c===============================================
	SUBROUTINE DRY_ADJUST(DPGQSP,PGQSP1,U,TP,LENSQ)
c===============================================
	IMPLICIT NONE

	INCLUDE 'mcons.inc'
	INCLUDE 'spcons.inc'
	INCLUDE 'mgrid.inc'
	INCLUDE 'spgrid.inc'

	COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
	COMPLEX PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
	REAL U(K1MAX, K2MAX, JX:JY, L1MAX) !VELOCITY FIELD
	REAL TP(K1MAX, K2MAX, L1MAX) !TEMPERATURE FIELD

	REAL TSTRESS(K1MAX, K2MAX, 0:1)
	REAL USTRESS(K1MAX, K2MAX, 0:1)
	REAL VSTRESS(K1MAX, K2MAX, 0:1)
	REAL THETA  (K1MAX, K2MAX, 0:1)
	INTEGER L, M, N
	INTEGER LM, LP
	REAL THLF, THET, RHO, GRHO, NSQT, NSQ
	REAL DELU, DELV, SHSQ, RI, SHEAR, DIFF
	REAL RHOD, GR2D
	REAL TTEND(K1MAX, K2MAX)
	REAL UTEND(K1MAX, K2MAX)
	REAL VTEND(K1MAX, K2MAX)

	REAL LENSQ
c	PARAMETER (LENSQ=900.)	!MIXING LENGTH SQUARED (M^2)

C...HDELP(L), RMIN(L), G0SQ HAVE BEEN ADDED TO COMMON BLOCKS
C=========================================================

	LM=1
	LP=0

C...No stress at uppermost half-level
	DO N=1,K2
	  DO M=1,K1
	    TSTRESS(M,N,LP)=0.
	    USTRESS(M,N,LP)=0.
	    VSTRESS(M,N,LP)=0.
	    THETA(M,N,LP) = TP(M,N,1)/PKLV(1)
	  ENDDO
	ENDDO

	DO L=1,L1
	  LM = MOD(LM+1,2)
	  LP = MOD(LM+1,2)

C...No stress at surface (level L1+1/2)
	  IF (L .EQ. L1) THEN
	    DO N=1,K2
	      DO M=1,K1
		TSTRESS(M,N,LP)=0.
		USTRESS(M,N,LP)=0.
		VSTRESS(M,N,LP)=0.
	      ENDDO
	    ENDDO
	  ELSE
	    DO N=1,K2
	      DO M=1,K1
		THETA(M,N,LP) = TP(M,N,L+1) / PKLV(L+1)
		THLF =  .5 * (TP(M,N,L+1) + TP(M,N,L))
		THET =  .5*(THETA(M,N,LM) + THETA(M,N,LP))
		RHO  = PHLV(L) / (RGAS * THLF)
		GRHO = G0SQ * RHO / HDELP(L)
		NSQT = GRHO*(THETA(M,N,LP)-THETA(M,N,LM))
		NSQ  = -NSQT / THET
		DELU = (U(M,N,1,L+1) - U(M,N,1,L))
		DELV = (U(M,N,2,L+1) - U(M,N,2,L))

C...RICHARDSON NUMBER AND DIFFUSION COEFFICIENT
		SHSQ = G0SQ*(DELU**2+DELV**2)*(RHO/HDELP(L))**2
		SHSQ = MAX(SHSQ,1E-10)
		RI = NSQ / SHSQ
		SHEAR = SQRT(SHSQ)

		IF (RI .GT. RMIN(L)) THEN
		  DIFF = 0.
		ELSE IF (RI .GT. 0.) THEN
		  DIFF = LENSQ * SHEAR * SQRT(1.-RI/RMIN(L))
		ELSE
		  DIFF = LENSQ * SHEAR
		ENDIF

C...CALCULATE STRESS AT HALF-LEVELS
		RHOD = RHO * DIFF
		GR2D = GRHO * RHOD

		TSTRESS(M,N,LP) = RHOD * NSQT
		USTRESS(M,N,LP) = GR2D * DELU
		VSTRESS(M,N,LP) = GR2D * DELV
	      ENDDO
	    ENDDO
	  ENDIF

C...    Add contribution from viscous stress divergence
	  DO N=1,K2
	    DO M=1,K1
	      TTEND(M,N)= PKLV(L) *
	1	   (TSTRESS(M,N,LP)-TSTRESS(M,N,LM))/DP(L)

	      UTEND(M,N)=
	1	   (USTRESS(M,N,LP)-USTRESS(M,N,LM))/DP(L)

	      VTEND(M,N)=
	1	   (VSTRESS(M,N,LP)-VSTRESS(M,N,LM))/DP(L)
	    ENDDO
	  ENDDO

	  CALL SPECTR(DPGQSP(0,0,JPOT,L),TTEND)
	  CALL CURLZ (DPGQSP(0,0,JVOR,L),UTEND,VTEND)
	  CALL DIVERG(DPGQSP(0,0,JDIV,L),UTEND,VTEND)
	ENDDO

	END