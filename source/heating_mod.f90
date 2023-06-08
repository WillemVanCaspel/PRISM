module heating_Mod

    use Declarations_mod
    use Prognos_sigma_mod
    use splib_mod
    use forcing_mod

    contains

!=======================================================================
      SUBROUTINE PY_HEAT(DPGQSP,TIM, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX)
!=======================================================================
!     Assimilate Python-generated diabatic heating fields [K s-1]
!     represented by spherical harmonic coeffs on sigprim model lvls.
 
      IMPLICIT NONE
 
    !   INCLUDE 'mcons.inc'
    !   INCLUDE 'spcons.inc'
    !   INCLUDE 'mgrid.inc'
    !   INCLUDE 'spgrid.inc'
 
      REAL TIM
      COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      INTEGER WLUN
 
      INTEGER I, J, L, M, N, ND, ND0, K, X, Y
 
      REAL CTIMW(20000)
 
      REAL DAY, WT0, WT1, CRANG, DELT
      REAL TIME0
      COMMON /HTIM/ TIME0 ! offset from tim so forcing can be shifted
      LOGICAL FIRST /.TRUE./
      
      REAL WATERCOEFF, OZONECOEFF, MESOCOEFF
      INTEGER TROPINDEX, STRATINDEX
 
!     I/O parameters
      INTEGER BRECL, RECLEN
      INTEGER TDIM, MAXM, MAXN, MAXL, TOPLVL, MXM, MXN
      
      REAL READN0(200,200,200)
      REAL READN1(200,200,200)
      
      COMPLEX HEATIN0(0:M1MAX,0:N1MAX,L1MAX)
      COMPLEX HEATIN1(0:M1MAX,0:N1MAX,L1MAX)
      REAL INTERPF(K1MAX,K2MAX)
 
      DOUBLE PRECISION XD(200), YD(200) ! set lon, lat dimensions input grid here
      DOUBLE PRECISION ZD(200,200)
      INTEGER NIN ! number of interpolation points (= K1 x K2)
      DOUBLE PRECISION SLON(100000), SLAT(100000), ZOUT0(100000), ZOUT1(100000)
      COMPLEX SPECHEAT(0:M1MAX,0:N1MAX)
      INTEGER COUNTER
      
!     Input and output folders and files
      CHARACTER*50 INFOLD, INFILE, OUTFOLD, OUTFILE, TIDEFILE, TTIDEFILE
      CHARACTER*50 INITTEMP, INITZONAL, INITGRID, INITZGRID, GWFDUMP, CURT
      CHARACTER*50 SURFFILE
      COMMON /FILES/ INFOLD, INFILE, SURFFILE, OUTFOLD, OUTFILE, &
         INITTEMP, INITZONAL, INITGRID, GWFDUMP, TIDEFILE, TTIDEFILE,CURT
 
!-----------------------------------------------------------------------
 
      WLUN = 77
      
      IF (FIRST) THEN
        OPEN(WLUN, FILE=TRIM(INFOLD)//TRIM(TIDEFILE), STATUS='OLD',SHARED, &
	         FORM='UNFORMATTED', ACCESS='DIRECT', RECL=1000)
 
        READ (WLUN,REC=1) RECLEN, TDIM, MAXM, MAXN, MAXL, TOPLVL
        WRITE(6,*) ' HEAT DIMS', TDIM, MAXM, MAXN, MAXL, TOPLVL
        CLOSE(WLUN)
 
	    BRECL=MAXM*MAXN*MAXL*4 + 8
 
        OPEN(WLUN, FILE=TRIM(INFOLD)//TRIM(TIDEFILE), STATUS='OLD',SHARED, &
          FORM='UNFORMATTED', ACCESS='DIRECT', RECL=BRECL)
 
        READ (WLUN,REC=2) RECLEN, (CTIMW(I),I=1,TDIM) ! DOY
        READ (WLUN,REC=3) RECLEN, (XD(I),I=1,MAXM) ! LONGITUDE
        READ (WLUN,REC=4) RECLEN, (YD(I),I=1,MAXN) ! LATITUDE
        
        ! artificial time shift if desired
        DO I=1,TDIM
          CTIMW(I) = CTIMW(I) - 0/24.
        ENDDO
        
        CRANG=CTIMW(TDIM)-CTIMW(1)
        ND0=0
              
        TOPLVL = TOPLVL + 1 ! adjust for Python indexing
 
        MXN = N1
        MXM = M1
        IF ((MAXN - 1) .LE. N1) MXN = (MAXN - 1) ! -1 for Python index
        IF ((MAXM - 1) .LE. M1) MXM = (MAXM - 1)
                      
        NIN = K1*K2
                
        FIRST=.FALSE.
      ENDIF
 
!     interpolate to current time and cycle periodically over data time range
      DAY=(TIM)/86400.
      DELT=MOD(DAY-REAL(CTIMW(1)),CRANG)
      DO WHILE (DELT .LT. 0)
        DELT=DELT+CRANG
      ENDDO
      DAY=CTIMW(1)+DELT
 
      ND=1
      DO J=2,TDIM-1
        IF (DAY .GT. REAL(CTIMW(J))) ND=J
      ENDDO
            
      IF (ND .NE. ND0) THEN
        ! read in tstep, interpolate and convert to spectral space
        READ (WLUN,REC=4+ND) RECLEN, &
          (((READN0(I,J,L), I=1, MAXM), J=1, MAXN), L=1, MAXL)
             
        COUNTER = 1 ! sigprim lat and lon interpolation grid
        DO X=1,K1
          DO Y=1,K2
            SLON(COUNTER) = DBLE(DEGLON(X))
            SLAT(COUNTER) = DBLE(DEGLAT(Y))
            COUNTER = COUNTER + 1
          ENDDO
        ENDDO
                     
        DO L=1,MAXL
 
          DO X=1,MAXM ! interpolate 2d 
            DO Y=1,MAXN
              ZD(X,Y) = DBLE(READN0(X,Y,L))
            ENDDO
          ENDDO
          
          CALL pwl_interp_2d (MAXM,MAXN,XD,YD,ZD(1:MAXM,1:MAXN), &
                              NIN,SLON,SLAT,ZOUT0)
              
          COUNTER = 1
          DO X=1,K1
            DO Y=1,K2
              INTERPF(X,Y) = REAL(ZOUT0(COUNTER))
              COUNTER = COUNTER + 1
            ENDDO
          ENDDO
                     
          CALL ZEROSP(SPECHEAT,N1)
          CALL SPECTR(SPECHEAT,INTERPF)
      
          DO M=0,M1
            DO N=0,N1
              HEATIN0(M,N,L) = SPECHEAT(M,N)
            ENDDO
          ENDDO
        ENDDO
    
        READ (WLUN,REC=5+ND) RECLEN, &
          (((READN1(I,J,L), I=1, MAXM), J=1, MAXN), L=1, MAXL)
                
        COUNTER = 1 ! sigprim lat and lon interpolation grid
        DO X=1,K1
          DO Y=1,K2
            SLON(COUNTER) = DBLE(DEGLON(X))
            SLAT(COUNTER) = DBLE(DEGLAT(Y))
            COUNTER = COUNTER + 1
          ENDDO
        ENDDO
             
        DO L=1,MAXL
          
          DO X=1,MAXM ! interpolate 2d 
            DO Y=1,MAXN
              ZD(X,Y) = DBLE(READN1(X,Y,L))
            ENDDO
          ENDDO
          
          CALL pwl_interp_2d (MAXM,MAXN,XD,YD,ZD(1:MAXM,1:MAXN), &
                              NIN,SLON,SLAT,ZOUT1)
         
          COUNTER = 1
          DO X=1,K1
            DO Y=1,K2
              INTERPF(X,Y) = REAL(ZOUT1(COUNTER))
              COUNTER = COUNTER + 1
            ENDDO
          ENDDO
                     
          CALL ZEROSP(SPECHEAT,N1)
          CALL SPECTR(SPECHEAT,INTERPF)
      
          DO M=0,M1
            DO N=0,N1
              HEATIN1(M,N,L) = SPECHEAT(M,N)
            ENDDO
          ENDDO
        ENDDO
        
        WRITE(6,*) 'HTIM:', ND, DAY, CTIMW(ND)
        ND0=ND
      ENDIF
                      
!     Linearly interpolate in time
      DELT=CTIMW(ND+1)-CTIMW(ND)
      WT0=(CTIMW(ND+1)-DAY)/DELT
      WT1=(DAY-CTIMW(ND))/DELT
                              
      DO L=TOPLVL,L1 ! TOPLVL replaced wtih 113
        K = L - TOPLVL + 1 ! match sigprim and heating file indexing
        DO N=0,MXN
          DO M=2,2 ! or 1, MXM
            IF (L .LE. STRATINDEX) THEN
              DPGQSP(M,N,JPOT,L)=(DPGQSP(M,N,JPOT,L) &
                + MESOCOEFF*WT0*HEATIN0(M,N,K) + MESOCOEFF*WT1*HEATIN1(M,N,K))
            ELSE IF (L .LE. TROPINDEX) THEN
              DPGQSP(M,N,JPOT,L)=(DPGQSP(M,N,JPOT,L) &
                + OZONECOEFF*WT0*HEATIN0(M,N,K) + OZONECOEFF*WT1*HEATIN1(M,N,K))
            ELSE
              DPGQSP(M,N,JPOT,L)=(DPGQSP(M,N,JPOT,L) &
                + WATERCOEFF*WT0*HEATIN0(M,N,K) + WATERCOEFF*WT1*HEATIN1(M,N,K))
            ENDIF
	      END DO
        END DO
      END DO
 
    END subroutine PY_HEAT

end module heating_Mod