c===============================================
      SUBROUTINE INITTEM(TZMEAN,TRADEQ)
c===============================================
c*****************
c.....TRADEQ=TEMP0
C*****************
c.....This version is for use in non-linear runs
c.....and is meant to replace cira_midrad.for
c.....allows free specification of initial wind, global mean temp profile
c.....and radiation equilibreum temperature

c.....Results are passed to INITVAR via COMMON /ZONE/
c.....Global mean temp profile TZMEAN and zonal mean radiative equilibrium
c.....state TRADEQ are to be sent to TINI
c.....Depends on prior call to VERINI

      IMPLICIT NONE

      INCLUDE 'mcons.inc'
      INCLUDE 'spcons.inc'
      INCLUDE 'mgrid.inc'
      INCLUDE 'spgrid.inc'
      INCLUDE 'sppoly.inc'

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

C     Input and output folders and files
      CHARACTER*50 INFOLD, INFILE, OUTFOLD, OUTFILE, TIDEFILE, TTIDEFILE
      CHARACTER*50 INITTEMP, INITZONAL, INITGRID, INITZGRID, GWFDUMP, CURT
      CHARACTER*50 SURFFILE
      COMMON /FILES/ INFOLD, INFILE, SURFFILE, OUTFOLD, OUTFILE,
     1   INITTEMP, INITZONAL, INITGRID, GWFDUMP, TIDEFILE, TTIDEFILE,CURT

C===================================================

C     Read in temperature profile file
      CLUN=123
      OPEN(UNIT=CLUN,STATUS='OLD',FILE=TRIM(INFOLD)//TRIM(INITTEMP),SHARED) 

C     Read in zonal wind profile file
      WLUN=124
      OPEN(UNIT=WLUN,STATUS='OLD',FILE=TRIM(INFOLD)//TRIM(INITZONAL),SHARED)

C     Input wind and temperatures on pressure grid (hPa)
      READ (CLUN,*) NALTT
      READ (CLUN,*) (CPREST(I),I=1,NALTT)
      READ (CLUN,*) (CTEMP(IALT),IALT=1,NALTT)

      READ (WLUN,*) NLATW, NALTW
      READ (WLUN,*) (CPRESW(I),I=1,NALTW), (CLATW(I),I=1,NLATW)
      READ (WLUN,*) ((CWIND(ILAT,IALT),ILAT=1,NLATW),IALT=1,NALTW)
      
      CLOSE(CLUN)
      CLOSE(WLUN)
      
CC     flip latitude to match python input file (index 0 = NP)
C     DO I=1,NALTW
C       DO J=1,NLATW
C         K = NLATW - J + 1
C         CFLIP(J,I) = CWIND(K,I)
C       ENDDO
C     ENDDO
C     
C     DO I=1,NALTW
C       DO J=1,NLATW
C         CWIND(J,I) = CFLIP(J,I)
C       ENDDO
C     ENDDO
      
c.....Convert pressure to sigma coordinates (Psurf in Pa)
      DO J=1,NALTW
        CPRESW(J)=-ALOG(CPRESW(J)*100./PSURF)
      ENDDO

      DO J=1,NALTT
        CPREST(J)=-ALOG(CPREST(J)*100./PSURF)
      ENDDO

C.....Interpolate global mean temperature profile onto model grid
      DO J=1,L1
        PH(J)=-ALOG(PLV(J))
      ENDDO

      CALL SPLINE(CPREST,CTEMP,NALTT,1.e30, 1.e30, Y2,IER) ! natural splines

      DO L=1,L1
        CALL SPLINT(CPREST,CTEMP,Y2,NALTT,PH(L),TZMEAN(L),IER)
      ENDDO

C.....Interpolate wind onto model grid (half levels)
      DO J=1,L1+1
        PH(J)=-ALOG(PHLV(J-1))
      ENDDO

      CALL SPLIE2(CLATW,CPRESW,CWIND,NLATW,NALTW,500,500,Y2A,IER)

      DO L=1,L1+1
        DO K=1,K2
          CALL SPLIN2
     1    (CLATW,CPRESW,CWIND,Y2A,NLATW,NALTW,500,500,DEGLAT(K),
     2     PH(L),UH(K,L),IER)
        ENDDO
      ENDDO

C.....y deriv of geopotential
      DO L=1,L1+1
        DO K=1,K2
          PHIY(K,L)=-(F0*MU(K)*A0 + UH(K,L) * TAN(PHI(K))) * UH(K,L)
        ENDDO
      ENDDO

C.....Wind at full levels
	DO L=1,L1
	  DO K=1,K2
	    U0(K,L)=.5*(UH(K,L)+UH(K,L+1))
	  END DO
	END DO

C.....Physical to spectral & Spectral to physical matrices
C.....restrict # of legendre polynomials to 40 (smoothness constraint)
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

C...anti-deriv wrt Y MATRIX
	CALL MATMUL(ADIP,N1MAX,ADI,N1MAX,P2SM,N1MAX,NL,NL,K2)
	CALL MATMUL(ADERIV,K2MAX,S2PM(1,2),K2MAX,ADIP,N1MAX,K2,NL,K2)

C...get y-deriv of temperature
	DO L=1, L1
	  DO K=1, K2
	    THETAY(K,L) = -(PHIY(K,L+1)-PHIY(K,L)) / CP * PKLV(L) /
	1	 (PHLV(L)**KAPPA-PHLV(L-1)**KAPPA)
	  ENDDO
	ENDDO

C...y-integral of THETAY gives global temp field with zero mean
	CALL MATMUL(TEMP0,K2MAX,ADERIV,K2MAX,THETAY,K2MAX,K2,K2,L1)

C...Add global mean temp
	DO K=1,K2
	  DO L=1,L1
	    TEMP0(K,L)=TEMP0(K,L)+TZMEAN(L)
	  ENDDO
	ENDDO

C...Geopotential at bottom boundary
	CALL MATMUL(GEO,K2MAX,ADERIV,K2MAX,PHIY(1,L1+1),K2MAX,K2,K2,1)

	write(6,*) 'Initializing tradeq = temp0'
	DO L=1,L1
	  DO K=1,K2
	    TRADEQ(K,L)=TEMP0(K,L)
	  ENDDO
	ENDDO

C.....Spectral components of zonal mean vorticity
C.....used in frictional constraint (sponge layer)
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
	END

c===============================================
	SUBROUTINE INITVAR(PGQSP0, PSSP0)
c===============================================

      	IMPLICIT NONE

      	INCLUDE 'mcons.inc'
      	INCLUDE 'spcons.inc'
      	INCLUDE  'mgrid.inc'
      	INCLUDE  'spgrid.inc'

      REAL TEMPPH(K1MAX, K2MAX)
      	REAL ZONW(K1MAX,K2MAX), MERW(K1MAX,K2MAX)
      	COMPLEX PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      	COMPLEX PSSP0(0:M1MAX,0:N1MAX)
      	INTEGER L, K, J, IPGQ, M

      REAL TEMP0(K2MAX,L1MAX), U0(K2MAX,L1MAX)
      COMPLEX VOR0(0:M1MAX,0:N1MAX,L1MAX)
      COMPLEX TMNLV(0:M1MAX,0:N1MAX,L1MAX)
      COMMON /ZONE/ TEMP0, U0, VOR0, TMNLV

	REAL GEO(K2MAX)
	REAL GEOP(K1MAX, K2MAX)
	COMMON /GEOBOT/ GEO, GEOP
        REAL QB(K1MAX,K2MAX)

      	EXTERNAL YSP2LA

        integer jran
        real ran,rat(0:n1max),dwl
	REAL TEMPp(K2MAX,L1MAX), Up(K2MAX,L1MAX), vp(k2max,l1max)

C...Spectral components of temperature field
      	DO L=1,L1
       	  DO IPGQ=1,3+NTRACE
	    CALL ZEROSP( PGQSP0(0,0,IPGQ,L), N1)
       	  END DO

	  DO K=1,K2
	    DO J=1,K1
	      TEMPPH(J,K) = TEMP0(K,L)
	    ENDDO
	  ENDDO

          CALL SPECTR( PGQSP0(0,0,JPOT,L), TEMPPH )
      	END DO

C...Set q=log(p_surf)
        CALL ZEROSP (PSSP0,N1)

        DO J=1,K1
          QB(J,1) = ALOG(PSURF)
        ENDDO

        DO K=2,K2
          QB(1,K) = QB(1,K-1) +
	1      F0 * A0 * MU(K) * U0(K,l1) / RGAS / TEMP0(K,l1)
	1      * (PHI(K)-PHI(K-1))

          DO J=2,K1
            QB(J,K) = QB(1,K)
          ENDDO
        ENDDO

        PSSP0(0,0)=ALOG(PSURF)


C.....Spectral components of wind field
      DO L=1,L1
        DO K=1,K2
          DO J=1,K1
	        ZONW(J,K) = U0(K,L)
	        MERW(J,K) = 0.
	      ENDDO
         ENDDO

C........Compute divergence & vorticity
         CALL CURLZ( PGQSP0(0,0,JVOR,L), ZONW, MERW)
         CALL DIVERG(PGQSP0(0,0,JDIV,L), ZONW, MERW) ! PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      END DO

       	RETURN
      	END

