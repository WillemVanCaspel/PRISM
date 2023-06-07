module IO_mod
    
    use Declarations_mod
    use Config_mod          , only: INPATH, OUTPATH, GRIDP, GRIDZ 
    use Prognos_sigma_mod

    implicit none

    integer, private, parameter :: IO_IND = 18

    ! routines used by other modules
    public :: inialt_read, init_1d_damp, inittem
     
    contains 

     !=============================================
              subroutine inialt_read(PTHICK, PP, ZS, PHALFL, PBOT, PLID, NLEV)
        !=============================================
        !   read model pressure levels from file
        
              IMPLICIT NONE
       
        !       ! INCLUDE 'mcons.inc'
              REAL*4 PTHICK(L1MAX), PHALFL(0:L1MAX), PP(L1MAX+1), ZS(L1MAX) ! 4 byte now because of old input...
              REAL PBOT, PLID
              INTEGER NLEV, L 
              
              INTEGER RECLEN, BRECL, ZLVLS
        
              INTEGER IARG
              COMMON /OPENS/ IARG
        
        !   Input and output folders and files
              CHARACTER*50 INFOLD, INFILE, OUTFOLD, OUTFILE, TIDEFILE, TTIDEFILE
              CHARACTER*50 INITTEMP, INITZONAL, INITGRID, INITZGRID, GWFDUMP, CURT
              CHARACTER*50 SURFFILE
              COMMON /FILES/ INFOLD, INFILE, SURFFILE, OUTFOLD, OUTFILE, &
                 INITTEMP, INITZONAL, INITGRID, GWFDUMP, TIDEFILE, TTIDEFILE,CURT
                    
              OPEN(UNIT=IO_IND,STATUS='OLD',FILE=TRIM(INPATH)//TRIM(GRIDP),SHARED, &
                  FORM='UNFORMATTED',ACCESS='DIRECT',RECL=100)
        
              READ (IO_IND,REC=1) RECLEN, NLEV
              CLOSE(IO_IND)
                    
              BRECL=NLEV*4 + 8
                    
              OPEN(UNIT=IO_IND,STATUS='OLD',FILE=TRIM(INPATH)//TRIM(GRIDP),SHARED, &
                  FORM='UNFORMATTED',ACCESS='DIRECT',RECL=BRECL)
        
              READ (IO_IND,REC=2) RECLEN, (PP(L),L=1,NLEV)
              CLOSE(IO_IND)
              
              WRITE(6,*) 'Number of vertical sigma-levels: ', NLEV
        
              PHALFL(0) = PLID
              PHALFL(NLEV) = PBOT
        
              DO L=1, NLEV-1
                PHALFL(L) = .5*(PP(L)+PP(L+1))
                PTHICK(L) = PHALFL(L)-PHALFL(L-1)
              ENDDO
        
              PTHICK(NLEV) = PHALFL(NLEV)-PHALFL(NLEV-1)

              !   read in (hydrostatic) geometric altitude lvls 
              OPEN(IO_IND,FILE=TRIM(INPATH)//TRIM(GRIDZ),STATUS='OLD', &
                   SHARED,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=100)
        
              READ (IO_IND,REC=1) RECLEN, ZLVLS
              CLOSE(IO_IND)
                         
              BRECL=ZLVLS*4 + 8
                    
              OPEN(IO_IND,FILE=TRIM(INPATH)//TRIM(GRIDZ),STATUS='OLD', &
                   SHARED,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=BRECL)
        
              READ (IO_IND,REC=2) RECLEN,(ZS(L),L=1,ZLVLS)
              CLOSE(IO_IND)
                    
              RETURN
            end subroutine inialt_read

        subroutine init_1d_damp(UDAMP, TDAMP, AMPBOT, DAYINV, NFAC, NLEV, PP, PBOT)

            IMPLICIT NONE
       
              REAL UDAMP(L1MAX), TDAMP(L1MAX)
              real*4 PP(L1MAX+1) ! 4 byte now because of old input
              REAL PBOT, NFAC, Z, DAYINV, AMPBOT
              INTEGER NLEV, L 
              
              INTEGER RECLEN, BRECL, ZLVLS

             !   Rayleigh surface friction and Newtonian cooling profiles
              DO L=1, NLEV
                Z=PP(L)/PBOT ! sigma-coordinates
                
                IF (Z .GT. 0.7) THEN
                   UDAMP(L)= AMPBOT * DAYINV * (0.7 - Z) / (0.7 - 1.0)
                ENDIF
                
        !     standard newtonian cooling profile
                TDAMP(L) = NFAC * HCOOL(PP(L)) * DAYINV 
              ENDDO

        end subroutine init_1d_damp

        !===============================================
        SUBROUTINE INITTEM(TZMEAN,TRADEQ)
  !===============================================
  !*****************
  !.....TRADEQ=TEMP0
  !*****************
  !.....This version is for use in non-linear runs
  !.....and is meant to replace cira_midrad.for
  !c.....allows free specification of initial wind, global mean temp profile
  !c.....and radiation equilibreum temperature
  
  !.....Results are passed to INITVAR via COMMON /ZONE/
  !c.....Global mean temp profile TZMEAN and zonal mean radiative equilibrium
  !c.....state TRADEQ are to be sent to TINI
  !c.....Depends on prior call to VERINI
  
        IMPLICIT NONE
  
        ! INCLUDE 'mcons.inc'
        ! INCLUDE 'spcons.inc'
        ! INCLUDE 'mgrid.inc'
        ! INCLUDE 'spgrid.inc'
        ! INCLUDE 'sppoly.inc'
  
        REAL UH(K2MAX,L1MAX+1)
        REAL TZMEAN(L1MAX), TRADEQ(K2MAX,L1MAX)
        INTEGER I, J, K, L
        INTEGER WLUN, CLUN
        INTEGER IALT, ILAT
        REAL CTEMP(500), CWIND(500,500), CPRESW(500), CLATW(500)
        REAL CFLIP(500,500)
        REAL CPREST(500)
        REAL TRAD(500,500)
        INTEGER NALTW, NALTT, NLATW, NLATT
        REAL TEMP0(K2MAX,L1MAX), U0(K2MAX,L1MAX)
        COMPLEX VOR0(0:M1MAX,0:N1MAX,L1MAX)
        COMPLEX TMNLV(0:M1MAX,0:N1MAX,L1MAX)
        COMMON /ZONE/ TEMP0, U0, VOR0, TMNLV
  
        REAL PHIY(K2MAX,L1MAX+1)
        REAL P2SM(N1MAX,K2MAX), S2PM(K2MAX,N1MAX), YS2PM(K2MAX,N1MAX)
        REAL AD(N1MAX,N1MAX), ADI(N1MAX,N1MAX), ADIP(N1MAX,K2MAX)
        REAL ADERIV(K2MAX,K2MAX)
        REAL THETAY(K2MAX,L1MAX+1)
        REAL GEO(K2MAX)
        REAL GEOP(K1MAX, K2MAX)
        COMMON /GEOBOT/ GEO, GEOP
  
        REAL PH(L1MAX+1)
        INTEGER NL
        REAL Y2A(500,500), Y2(500)
        INTEGER IER, N
        REAL ZONW(K1MAX,K2MAX), MERW(K1MAX,K2MAX)
        COMPLEX VORT(0:M1MAX,0:N1MAX)
  
  !     Input and output folders and files
        CHARACTER*50 INFOLD, INFILE, OUTFOLD, OUTFILE, TIDEFILE, TTIDEFILE
        CHARACTER*50 INITTEMP, INITZONAL, INITGRID, INITZGRID, GWFDUMP, CURT
        CHARACTER*50 SURFFILE
        COMMON /FILES/ INFOLD, INFILE, SURFFILE, OUTFOLD, OUTFILE, &
           INITTEMP, INITZONAL, INITGRID, GWFDUMP, TIDEFILE, TTIDEFILE,CURT
  
  !===================================================
  
  !     Read in temperature profile file
        CLUN=123
        OPEN(UNIT=CLUN,STATUS='OLD',FILE=TRIM(INFOLD)//TRIM(INITTEMP),SHARED) 
  
  !     Read in zonal wind profile file
        WLUN=124
        OPEN(UNIT=WLUN,STATUS='OLD',FILE=TRIM(INFOLD)//TRIM(INITZONAL),SHARED)
  
  !     Input wind and temperatures on pressure grid (hPa)
        READ (CLUN,*) NALTT
        READ (CLUN,*) (CPREST(I),I=1,NALTT)
        READ (CLUN,*) (CTEMP(IALT),IALT=1,NALTT)
  
        READ (WLUN,*) NLATW, NALTW
        READ (WLUN,*) (CPRESW(I),I=1,NALTW), (CLATW(I),I=1,NLATW)
        READ (WLUN,*) ((CWIND(ILAT,IALT),ILAT=1,NLATW),IALT=1,NALTW)
        
        CLOSE(CLUN)
        CLOSE(WLUN)
        
  !      flip latitude to match python input file (index 0 = NP)
  !     DO I=1,NALTW
  !       DO J=1,NLATW
  !         K = NLATW - J + 1
  !         CFLIP(J,I) = CWIND(K,I)
  !       ENDDO
  !     ENDDO
  !     
  !     DO I=1,NALTW
  !       DO J=1,NLATW
  !         CWIND(J,I) = CFLIP(J,I)
  !       ENDDO
  !     ENDDO
        
  !.....Convert pressure to sigma coordinates (Psurf in Pa)
        DO J=1,NALTW
          CPRESW(J)=-ALOG(CPRESW(J)*100./PSURF)
        ENDDO
  
        DO J=1,NALTT
          CPREST(J)=-ALOG(CPREST(J)*100./PSURF)
        ENDDO
  
  !.....Interpolate global mean temperature profile onto model grid
        DO J=1,L1
          PH(J)=-ALOG(PLV(J))
        ENDDO
  
        CALL SPLINE(CPREST,CTEMP,NALTT,1.e30, 1.e30, Y2,IER) ! natural splines
  
        DO L=1,L1
          CALL SPLINT(CPREST,CTEMP,Y2,NALTT,PH(L),TZMEAN(L),IER)
        ENDDO
  
  !.....Interpolate wind onto model grid (half levels)
        DO J=1,L1+1
          PH(J)=-ALOG(PHLV(J-1))
        ENDDO
  
        CALL SPLIE2(CLATW,CPRESW,CWIND,NLATW,NALTW,500,500,Y2A,IER)
  
        DO L=1,L1+1
          DO K=1,K2
            CALL SPLIN2 (CLATW,CPRESW,CWIND,Y2A,NLATW,NALTW,500,500,DEGLAT(K), &
             PH(L),UH(K,L),IER)
          ENDDO
        ENDDO
  
  !.....y deriv of geopotential
        DO L=1,L1+1
          DO K=1,K2
            PHIY(K,L)=-(F0*MU(K)*A0 + UH(K,L) * TAN(PHI(K))) * UH(K,L)
          ENDDO
        ENDDO
  
  !.....Wind at full levels
      DO L=1,L1
        DO K=1,K2
          U0(K,L)=.5*(UH(K,L)+UH(K,L+1))
        END DO
      END DO
  
  !.....Physical to spectral & Spectral to physical matrices
  !C.....restrict # of legendre polynomials to 40 (smoothness constraint)
      NL=N1
      IF (N1 .GT. 40) NL=40
  
      DO I=1,NL+1
        DO J=1,K2
          P2SM(I,J)=P2(0,I-1,J)*G(J)
          S2PM(J,I)=P2(0,I-1,J)
          YS2PM(J,I)=-H2(0,I-1,J)*COSINV(J)
        ENDDO
      ENDDO
  
      CALL MATMUL(AD,N1MAX,P2SM,N1MAX,YS2PM(1,2),K2MAX,NL,K2,NL)
      CALL MATINV(ADI,N1MAX,AD,N1MAX,NL)
  
  !...anti-deriv wrt Y MATRIX
      CALL MATMUL(ADIP,N1MAX,ADI,N1MAX,P2SM,N1MAX,NL,NL,K2)
      CALL MATMUL(ADERIV,K2MAX,S2PM(1,2),K2MAX,ADIP,N1MAX,K2,NL,K2)
  
  !...get y-deriv of temperature
      DO L=1, L1
        DO K=1, K2
          THETAY(K,L) = -(PHIY(K,L+1)-PHIY(K,L)) / CP * PKLV(L) / &
       	 (PHLV(L)**KAPPA-PHLV(L-1)**KAPPA)
        ENDDO
      ENDDO
  
  !...y-integral of THETAY gives global temp field with zero mean
      CALL MATMUL(TEMP0,K2MAX,ADERIV,K2MAX,THETAY,K2MAX,K2,K2,L1)
  
  !...Add global mean temp
      DO K=1,K2
        DO L=1,L1
          TEMP0(K,L)=TEMP0(K,L)+TZMEAN(L)
        ENDDO
      ENDDO
  
  !...Geopotential at bottom boundary
      CALL MATMUL(GEO,K2MAX,ADERIV,K2MAX,PHIY(1,L1+1),K2MAX,K2,K2,1)
  
      write(6,*) 'Initializing tradeq = temp0'
      DO L=1,L1
        DO K=1,K2
          TRADEQ(K,L)=TEMP0(K,L)
        ENDDO
      ENDDO
  
  !.....Spectral components of zonal mean vorticity
  !.....used in frictional constraint (sponge layer)
      DO L=1,L1
        DO K=1,K2
          DO J=1,K1
            ZONW(J,K) = U0(K,L)
            MERW(J,K) = 0.
          ENDDO
        ENDDO
  
        CALL ZEROSP(VORT,N1)
        CALL CURLZ( VORT, ZONW, MERW)
        DO N=0,N1
          VOR0(0,N,L)=REAL(VORT(0,N))
        ENDDO
  
        DO N=41,N1
          VOR0(0,N,L)= (0.,0.)
        ENDDO
      END DO
  
      RETURN
    END subroutine INITTEM
  
        
        !============================================
              REAL FUNCTION HCOOL(PRESSURE)
        !============================================

              implicit none
              integer I
              real*4 PRESSURE

        !   from jung-hee, based on a hoskins paper; newtonian cooling due to O3 and CO2
        	real, dimension(81) :: COOL = (/                                            &
        	     0.05,                                                                  &
        	     0.1293,  0.1411,  0.1535,  0.1663,  0.1791,  0.1918,  0.2040,  0.2153, &
        	     0.2255,  0.2342,  0.2412,  0.2463,  0.2492,  0.2500,  0.2497,  0.2494, &
        	     0.2491,  0.2487,  0.2483,  0.2479,  0.2474,  0.2468,  0.2462,  0.2456, &
        	     0.2448,  0.2440,  0.2431,  0.2421,  0.2410,  0.2398,  0.2384,  0.2369, &
        	     0.2353,  0.2335,  0.2316,  0.2294,  0.2271,  0.2246,  0.2218,  0.2188, &
        	     0.2156,  0.2121,  0.2083,  0.2043,  0.2000,  0.1954,  0.1906,  0.1854, &
        	     0.1800,  0.1744,  0.1685,  0.1624,  0.1561,  0.1497,  0.1431,  0.1364, &
        	     0.1297,  0.1229,  0.1162,  0.1095,  0.1029,  0.0825,  0.0614,  0.0502, &
        	     0.0429,  0.0400,  0.0400,  0.0400,  0.0400,  0.0400,  0.0400,  0.0400, &
        	     0.0400,  0.0400,  0.0400,  0.0400,  0.0420,  0.0530,  0.0713,  0.1106/)
        
        	real, dimension(81) :: CPRES = (/                                  &
        	     0.0927,                                                       &
        	     0.2002,      0.2228,      0.2480,      0.2760,      0.3072,   &
        	     0.3420,      0.3807,      0.4237,      0.4716,      0.5250,   &
        	     0.5844,      0.6504,      0.7240,      0.8059,      0.8970,   &
        	     0.9985,      1.1114,      1.2371,      1.3770,      1.5327,   &
        	     1.7061,      1.8990,      2.1138,      2.3528,      2.6189,   &
        	     2.9151,      3.2448,      3.6117,      4.0202,      4.4749,   &
        	     4.9809,      5.5442,      6.1713,      6.8692,      7.6460,   &
        	     8.5108,      9.4733,     10.5447,     11.7372,     13.0646,   &
        	     14.5421,     16.1867,     18.0173,     20.0549,     22.3230,  &
        	     24.8476,     27.6577,     30.7856,     34.2672,     38.1426,  &
        	     42.4563,     47.2578,     52.6024,     58.5513,     65.1731,  &
        	     72.5437,     80.7479,     89.8800,    100.0450,    111.3590,  &
        	     123.9530,    137.9710,    153.5750,    170.9430,    190.2760, &
        	     211.7950,    235.7470,    262.4090,    292.0860,    325.1180, &
        	     361.8870,    402.8140,    448.3700,    499.0770,    555.5200, &
        	     618.3450,    688.2760,    766.1150,    852.7579,    949.1990/)
        
              IF (PRESSURE/1.E2 .LT. CPRES(1)) THEN
                HCOOL = .05
              ELSE IF (PRESSURE/1.E2 .GT. CPRES(81)) THEN
                HCOOL = .11
              ELSE
                DO I=1,80
                  IF ((PRESSURE/1.E2 .GT. CPRES(I)) .AND. &
                      (PRESSURE/1.E2 .LE. CPRES(I+1)))    &
                       HCOOL=COOL(I) + (PRESSURE/1.E2-CPRES(I))*(COOL(I+1)-COOL(I))/ &
                      (CPRES(I+1)-CPRES(I))
                ENDDO
              ENDIF
        
              RETURN
    END FUNCTION HCOOL
        
        ! !=============================================
        !       REAL FUNCTION CMAM_KZZ(Z, T)
        ! !=============================================
        ! !   Eddy + molecular Kzz profile (M**2 s**-1)
        ! !   Eddy Formula models profile in McLandress 2002
        ! !   Molec from Banks, as used in WACCM
        
        !       IMPLICIT NONE
              
        !       REAL Z, T, EDDY, SLOPE
        
        !       IF (Z .LT. 15) THEN
        !         EDDY=1.
        !       ELSE IF ((Z .GE. 15.) .AND. (Z .LT. 107.5)) THEN
        !         SLOPE=(ALOG(200.)-ALOG(1.))/92.5
        !         EDDY=EXP((Z-15.)*SLOPE)
        !       ELSE IF ((Z .GE. 107.5) .AND. (Z .LT. 115.)) THEN
        !         EDDY=200.
        !       ELSE IF (Z .GE. 115) THEN
        !         SLOPE=(ALOG(8.)-ALOG(200.))/15.
        !         EDDY=200.*EXP((Z-115)*SLOPE)
        !       ENDIF
        
        !       CMAM_KZZ=EDDY
        
        !       RETURN
        !       END
        
        ! !=========================================================
        ! !   SUBROUTINE SPEWPY(PGQSP, PSSP, OUTFILE)
        ! !=========================================================
        ! !   Write zonal and meridional wind, & vertical velocities,
        ! !   geopotential, potential temperature, for Python I/O.
        ! !
        ! !   IMPLICIT NONE
        ! !
        ! !   INCLUDE 'mcons.inc'
        ! !   INCLUDE 'spcons.inc'
        ! !   INCLUDE 'mgrid.inc'
        ! !   INCLUDE 'spgrid.inc'
        ! !
        ! !   COMPLEX PGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
        ! !   COMPLEX PSSP(0:M1MAX,0:N1MAX)
        ! !   INTEGER M, N, L, IREC
        ! !
        ! !   REAL GEO(K2MAX)
        ! !   REAL GEOB(K1MAX, K2MAX)
        ! !   COMMON /GEOBOT/ GEO, GEOB
        ! !   
        ! !   CHARACTER*50 OUTFILE
        ! !
        ! !   IREC=IREC+1
        ! !	  WRITE(6,*) 'Writing Python record: ', IREC, ' to  ', OUTFILE
        ! !
        ! !   WRITE(13) (((PGQSP(M,N,JVOR,L),M=0,M1),N=0,N1),L=1,L1)
        ! !   WRITE(13) (((PGQSP(M,N,JDIV,L),M=0,M1),N=0,N1),L=1,L1)
        ! !   WRITE(13) (((PGQSP(M,N,JPOT,L),M=0,M1),N=0,N1),L=1,L1)
        ! !   WRITE(13) ((PSSP(M,N),M=0,M1),N=0,N1)
        ! !   WRITE(13) ((GEOB(M,N),M=1,K1),N=1,K2)
        
        ! !   END
        
        ! !=========================================================
        !       SUBROUTINE SPEWPY(PGQSP, PSSP, OUTFILE)
        ! !=========================================================
        ! !   Write zonal (ZW) and meridional wind (MW), 
        ! !   surf geopotential (GEOB), temperature (TPR), surf pres (PRS)
        
        !       IMPLICIT NONE
        
        !       INCLUDE 'mcons.inc'
        !       INCLUDE 'spcons.inc'
        !       INCLUDE 'mgrid.inc'
        !       INCLUDE 'spgrid.inc'
        
        !       COMPLEX PGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
        !       COMPLEX PSSP(0:M1MAX,0:N1MAX)
        !       COMPLEX FSP0(0:M1MAX,0:N1MAX),FSP1(0:M1MAX,0:N1MAX)
        !       INTEGER M, N, L, IREC
              
        !       REAL ZW(K1MAX,K2MAX,L1MAX), MW(K1MAX,K2MAX,L1MAX)
        !       REAL TPR(K1MAX,K2MAX,L1MAX)
        !       REAL PRS(K1MAX,K2MAX)
        
        !       REAL GEO(K2MAX)
        !       REAL GEOB(K1MAX, K2MAX)
        !       COMMON /GEOBOT/ GEO, GEOB
              
        !       CHARACTER*50 OUTFILE
        
        !       IREC=IREC+1
        ! 	  WRITE(6,*) 'Writing Python record: ', IREC, ' to  ', OUTFILE
                       
        !       DO L=1,L1
        ! !...... compute streamfunction PSI and velocity potential CHI
        !         CALL IDELSQ( FSP0, PGQSP(0,0,JVOR,L), N1)
        !         CALL IDELSQ( FSP1, PGQSP(0,0,JDIV,L), N1)
        ! !.......transform to physical space
        !         CALL HVELOC( ZW(1,1,L), MW(1,1,L), FSP0, FSP1)
        !         CALL PHYSIC( TPR(1,1,L), PGQSP(0,0,JPOT,L))
        !       ENDDO
            
        !       CALL PHYSIC(PRS,PSSP)
              
        !       WRITE(13) (((ZW(M,N,L),M=1,K1),N=1,K2),L=1,L1)
        !       WRITE(13) (((MW(M,N,L),M=1,K1),N=1,K2),L=1,L1)
        !       WRITE(13) (((TPR(M,N,L),M=1,K1),N=1,K2),L=1,L1)
        !       WRITE(13) ((PRS(M,N),M=1,K1),N=1,K2)
        !       WRITE(13) ((GEOB(M,N),M=1,K1),N=1,K2)
              
        !       END
end module IO_mod