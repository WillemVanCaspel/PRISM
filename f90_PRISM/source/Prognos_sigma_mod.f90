

module Prognos_sigma_mod

    use Declarations_mod
    use coolmd_mod
    use splib_mod 

    implicit none

    contains

    !========================================================================
      SUBROUTINE PLINI(OMEGA, R, CPRESS, GRAVIT)
    !========================================================================
!     PLINI initializes the common block /PLPARM/ which contains the
!     planetary parameters.
!     Input:  REAL OMEGA -- angular velocity of planetary rotation (in 1/s)
!             REAL R -- gas constant {in J/(kg*K)}
!             REAL CPRESS -- specific heat at const. pressure {in J/(kg*K)}
!             REAL GRAVIT -- gravitational acceleration (in m/s**2)
!     Output: Common block /PLPARM/ is initialized

      IMPLICIT NONE    
      REAL OMEGA, R, CPRESS, GRAVIT

      OMEGA0= OMEGA
      F0= 2.0 * OMEGA0
      FSP01= F0 / SQRT(3.0)
      RGAS= R
      CP= CPRESS
      KAPPA= RGAS / CP
      G0= GRAVIT
      G0SQ=G0*G0
      
      RETURN
      END
      
!========================================================================
      SUBROUTINE VERINI( NLEV, PTHICK, PHALFL, PP, DEL8DF, &
                         TDAMP, UDAMP, SFRIC)
!========================================================================
!
!     VERINI initializes the vertical truncation for a spectral pressure
!     coordinate model (using the Lorenz (1960) finite-differencing, with
!     vertically averaged divergence set to advection of surface pressure).
!
!     Note: Prior to calling VERINI, the horizontal truncation should have
!           been initialized by a call to routine SPINI in module sptran.F.
!           Also, the planetary parameters should have been initialized by
!           a call to routine PLINI in this module.
!
!     Note: This version of VERINI only sets up pressure grid and damping stuff.
!           Other things formerly set up in the old VERINI are now in TINI.
!           This makes it easier to use the vertical grid in initialization
!           subroutines.
!
!     NLEV is the number of layers in the vertical truncation, (NLEV >= 2).
!
!     PTHICK(NLEV) is the PRESSURE-thickness of each level.
!         Note: The upper boundary boundary is always assumed to be at zero
!         pressure.Let PSURF=sum(PTHICK).  Sigma coordinates are defined 
!         as sigma(i)-sigma(i-1) = pthick(i)/PSURF at the layer boundaries
!         (half-levels).)
!     PP(nlev) is the pressure at the full level.
!     NOTE: It is no longer assumed that this is midway between half-levels
!           DEL8DF is the del-8 (quad-harmonic) diffusion coefficient
!           (in units of A0**8/sec, where A0 is the planetary radius)
!     TDAMP(NLEV)is the sigma dependent radiative damping coefficient
!     (in units of 1/sec) that determines relaxation of temp back to
!     TRADEQ.
!     UDAMP(NLEV) is the sigma dependent Rayleigh friction coefficient
!     (in units of  1/sec) that damps velocity.
!     SFRIC is the surface rayleigh friction coefficient (in units of 1/sec)
!     that damps momentum in the lowermost level.
!
!     Lorenz, E.N. 1960: Energy and numerical weather prediction,
!     TELLUS,vol.12,pp.364-373.

      IMPLICIT NONE

      INTEGER NLEV
      REAL*4 PTHICK(NLEV), PP(NLEV), PHALFL(0:NLEV) 
      REAL DEL8DF
      REAL TDAMP(NLEV), UDAMP(NLEV), SFRIC, DELZ
      INTEGER L, L0, N
 
      IF (NLEV.GT.L1MAX) CALL SPERR( 'VERINI', 'Too many sigma levels')
      L1= NLEV
 
!     Initialize pressure related parameters
!     PLV and PHLV are now sigma-coord, NOT pressure-levels
      PSURF = PHALFL(L1) ! 1e5 Pa
     
      PHLV(0)= PHALFL(0) / PSURF
      DO 10 L=1,L1
        DP(L)= PTHICK(L)/PSURF
        HDP(L)= 0.5 * DP(L)
        PLV(L)= PP(L) / PSURF
        PHLV(L)= PHALFL(L) / PSURF
   10   CONTINUE

      DO L=1,L1-1
        HDELP(L) = HDP(L)+HDP(L+1)
      ENDDO

!     Standard vertical temp profile, and coordinate PK
      DO 20 L=1,L1
        PKLV(L)= PLV(L)**KAPPA
   20   CONTINUE

!     Compute PK-cap
      DO 30 L=1,L1-1
   30   PKCHLV(L)= 0.5 * (PKLV(L+1) - PKLV(L))
!
!     D to W conversion matrix (Z of Haltiner-Williams)
      DO 70 L=1,L1-1
        DO 70 L0=1,L1
          IF (L0.LE.L) THEN
            D2W(L,L0)= (PHLV(L)-1)*DP(L0)
          ELSE
            D2W(L,L0)= PHLV(L)*DP(L0)
          ENDIF
 70    CONTINUE

!     temp to GeoPotential conversion matrix (Hydrostatic relation)
       DO 90 L=1,L1
          DO L0=1,L1-1
             T2GPCP(L,L0)=-PHLV(L0)*CP*PKCHLV(L0)/PKLV(L0) + RGAS*DP(L0)
          ENDDO
          T2GPCP(L,L1)=RGAS*DP(L1)

          DO L0=1,L1-1
             T2GPCP(L,L0+1) = T2GPCP(L,L0+1)-PHLV(L0)*CP*PKCHLV(L0)/PKLV(L0+1)
          ENDDO
 90    CONTINUE

       DO 95 L=1, L1-1
          DO L0=L,L1-1
             T2GPCP(L,L0)   = T2GPCP(L,L0)   + CP*PKCHLV(L0)/PKLV(L0)
             T2GPCP(L,L0+1) = T2GPCP(L,L0+1) + CP*PKCHLV(L0)/PKLV(L0+1)
          ENDDO
 95    CONTINUE

!     Del**8 sub-grid-scale diffusion
      QHDAMP= DEL8DF
      DO 100 N=0,N1

!       Scale-dependent del**8 damping factor
!       (Note: QHDAMP is in units such that A0 = 1)
        TD8FAC(N) =  QHDAMP * FLOAT( N*(N+1)   )**4
        UD8FAC(N) =  QHDAMP * FLOAT( N*(N+1)-2 )**4
  100   CONTINUE
!
!     Mean thermal state, and Rayleigh/Newtonian damping coefficients
      DO 130 L=1,L1
        URLXLV(L)= UDAMP(L)
 130    TRLXLV(L)= TDAMP(L)

!     Surface drag (in units of (sigma/time) )
      USDRAG= DP(L1) * SFRIC

!     DRY ADJUSTMENT LIMIT ON RICHARDSON NUMBER
      DO L=1, L1
        DELZ = 7.*DP(L)/PLV(L)
        RMIN(L) = .25 * (1. + .1 * (10.*DELZ)**1.5)
      ENDDO
        
      RETURN
      END

!========================================================================
      SUBROUTINE TINI(TZSTD, TSTEEP, TRADEQ, VERVIS, DIFFC)
!========================================================================
!     TZSTD(NLEV) is the vertical profile of the mean temperature distribution.
!     Note: A more stably stratified version of TZSTD is used to define
!           the reference temp profile for the semi-implicit time-stepping.
!           TZSTD is also used to calculate standard height-levels
!           corresponding to the sigma levels)
!
!     TSTEEP is the steepness factor for the reference vertical temp
!     profile used for semi-implicit time-stepping. i.e. The reference
!     temp profile is steeper than TZSTD by a factor TSTEEP.
!     (One would typically choose TSTEEP > 1.0)
!
!     TRADEQ(K2MAX,NLEV) is the latitude-sigma field of the "radiative
!     equilibrium" temp distribution.
!
!     VERVIS is the vertical viscosity coefficient (in units of m**2/s) that
!     diffuses momentum.
!
!     DIFFC is thermal diffusion coefficient

      IMPLICIT NONE

      REAL TZSTD(L1MAX), TSTEEP
      REAL TRADEQ(K2MAX, L1MAX),VERVIS(l1max), DIFFC(l1max)
      REAL W2TT(L1MAX,L1MAX), FPH(K1MAX,K2MAX)
      COMPLEX FSP(0:M1MAX,0:N1MAX)
      INTEGER J, K, N, L, L0

      REAL WDAMP(L1MAX), MDAMP(L1MAX), TDIFF(L1MAX)
      REAL ALPHAT(L1MAX), ALPHAW(L1MAX)
      INTEGER SDTRUNC, DAMPTRUNC
      COMMON /WAVD/ WDAMP, MDAMP, TDIFF, ALPHAT, ALPHAW, SDTRUNC, DAMPTRUNC
	  
      REAL TEMP0(K2MAX,L1MAX), U0(K2MAX,L1MAX)
      COMPLEX VOR0(0:M1MAX,0:N1MAX,L1MAX)
      COMMON /ZONE/ TEMP0, U0, VOR0

      DO L=1,L1
        TSTDLV(L)= TZSTD(L)
      ENDDO

!     Compute standard height values, using the hydrostatic relation.
!     Also initialize reference temp profile(for semi-implicit scheme).
!     The reference temp profile is simply a steeper version of the
!     standard temp profile (steeper by a factor TRFFAC).
      TRFFAC= TSTEEP
      ZSTDLV(L1)= CP*(1.0-PKLV(L1))*TSTDLV(L1)/PKLV(L1) / G0

!     Reference temp profile coincides with standard profile at surface
      TREFLV(L1)= TSTDLV(L1)
      DO 40 L=L1-1,1,-1
        ZSTDLV(L)= ZSTDLV(L+1) + CP*PKCHLV(L) * &
                  (TSTDLV(L+1)/PKLV(L+1) +TSTDLV(L)/PKLV(L)) / G0
        TREFLV(L)= TRFFAC * TSTDLV(L) &
             + PKLV(L)/PKLV(L+1) * (TREFLV(L+1) - TRFFAC * TSTDLV(L+1))
   40   CONTINUE
!
!     W to temp tendency conversion matrix (Static stability)
!     NOTE: W2TT = - Y OF Haltiner-Williams
      DO 80 L=1,L1
        DO 80 L0=1,L1-1
          IF (L0.EQ.(L-1)) THEN
            W2TT(L,L0)=-0.5 * (TREFLV(L)-PKLV(L)/PKLV(L-1)*TREFLV(L-1)) / DP(L)
          ELSE IF (L0.EQ.L) THEN
            W2TT(L,L0)=-0.5 * (TREFLV(L+1)*PKLV(L)/PKLV(L+1)-TREFLV(L)) / DP(L)
          ELSE
            W2TT(L,L0)= 0.0
          ENDIF
   80     CONTINUE

!     D to temp-tendency conversion matrix
!     (Composite  W2TT # D2W - kappa * T#DP)
!     NNT2DT = -YZ OF HALTINER AND WILLIAMS
!     NOTE: THIS HAS OPPOSITE SIGN TO THAT USED IN DOCUMENTATION
      CALL MATMUL(NNT2DT, L1MAX, W2TT, L1MAX, D2W, L1MAX,L1, L1-1, L1)

      DO 85 L=1, L1
         DO 85, L0=1,L1
 85         D2TT(L,L0) = NNT2DT(L,L0) - KAPPA * TREFLV(L) * DP(L0)

!.....Compute spectral coefficients of zonally symmetric mean temp
      DO 130 L=1,L1
        DO 110 K=1,K2
          DO 110 J=1,K1
  110       FPH(J,K)= TRADEQ(K,L)
        CALL ZEROSP(FSP, N1)
        CALL SPECTR(FSP, FPH)
        DO 120 N=0,N1
  120     TMNLV(0,N,L)= REAL( FSP(0,N) )
  130   CONTINUE

!     Convert viscosity coefficient from z-coordinate to SIGMA-coordinate,
!     using the standard height values
      DO 140 L=1,L1-1
        VVISC(L)= VERVIS(L) * ( (HDP(L+1)  + HDP(L)     ) /  &
                             (ZSTDLV(L) - ZSTDLV(L+1))  )**2

        IF (VERVIS(L) .LT. 0) &
            VVISC(L)= VERVIS(L) * ( (HDP(L+1)  + HDP(L)     ) / &
                            (ZSTDLV(L) - ZSTDLV(L+1))  )**4

        TDIFF(L)= DIFFC(L) * ( (HDP(L+1)  + HDP(L)     ) /  & 
                            (ZSTDLV(L) - ZSTDLV(L+1))  )**2
!        TDIFF(L)= DIFFC(L) /PKLV(L) * ( (HDP(L+1)  + HDP(L)     ) /
!     1                       (ZSTDLV(L) - ZSTDLV(L+1))  )**2

        IF (DIFFC(L) .LT. 0) &
             TDIFF(L)= DIFFC(L) * ( (HDP(L+1)  + HDP(L)     ) / &
                             (ZSTDLV(L) - ZSTDLV(L+1))  )**4

  140   CONTINUE
      RETURN
      END
 
 
!========================================================================
      SUBROUTINE MATMUL(A, NA0, B, NB0, C, NC0, NROWB, NCOLB, NCOLC)
!========================================================================
! MATMUL multiplies two matrices together to produce a third matrix.
!  A = B . C
! Input:  INTEGER NA0, NB0, NC0 -- inner dimension of matrices A, B, C
!         INTEGER NROWB, NCOLB -- no. of rows and columns of matrix B
!         INTEGER NCOLC -- no. of columns of matrix C
!         REAL B(NB0,NCOLB), C(NC0,NCOLC) -- matrices to be multiplied
! Output: REAL A(NA0,NCOLC) -- matrix A = B . C (with NROWB rows)
!
      IMPLICIT NONE

      INTEGER NA0, NB0, NC0, NROWB, NCOLB, NCOLC
      REAL A(NA0,NCOLC), B(NB0,NCOLB), C(NC0,NCOLC)
      INTEGER I, J, K
      REAL TEMP
      DO 20 I=1,NROWB
        DO 20 K=1,NCOLC
          TEMP= 0.0
          DO 10 J=1,NCOLB
   10       TEMP= TEMP + B(I,J) * C(J,K)
          A(I,K)= TEMP
   20     CONTINUE
      RETURN
      END
 
!========================================================================
      SUBROUTINE MATINV (AINV, NAINV0, A, NA0, NA)
!========================================================================
! MATINV finds the inverse of a square matrix A.
! Input:  INTEGER NAINV0, NA0 -- inner dimension of matrices AINV, A
!         INTEGER NA -- no. of rows/columns in A
! (Note: NA must be <= L1MAX)
!         REAL A(NA0,NA) -- matrix A
! Output: REAL AINV(NAINV0,NA) -- inverse of A
!
      IMPLICIT NONE

      INTEGER MAXDIM
      PARAMETER (MAXDIM=300)

      INTEGER NAINV0, NA0, NA
      REAL AINV(NAINV0,NA), A(NA0,NA)
      REAL*8 A2(MAXDIM,MAXDIM), B(MAXDIM,MAXDIM), C(MAXDIM,MAXDIM)
      REAL*8 WKSPCE(MAXDIM), AA(MAXDIM,MAXDIM), BB(MAXDIM,MAXDIM)
      INTEGER I, J, IFAIL
!========================================================================
!
      IF (NA.GT.MAXDIM) CALL SPERR( 'MATINV', 'Matrix too large')

!     Initialize temporary identity matrix
      DO 10 I=1,NA
        DO 10 J=1,NA
!         Copy matrix to correct precision
          A2(I,J)= A(I,J)
!         Create Identity matrix for RHS
          IF (I.EQ.J) THEN
            B(I,J)= 1.0
          ELSE
            B(I,J)= 0.0
          ENDIF
   10     CONTINUE
!     Call Nag routine
      IFAIL= 0
      CALL F04AEF(A2, MAXDIM, B, MAXDIM, NA, NA, C, MAXDIM, &
                  WKSPCE, AA, MAXDIM, BB, MAXDIM, IFAIL)
      IF (IFAIL.NE.0) CALL SPERR( 'MATINV', 'Unable to invert matrix')
!     Copy inverse back to original precision
      DO 20 I=1,NA
        DO 20 J=1,NA
          AINV(I,J)= C(I,J)
   20     CONTINUE

      RETURN
      END
 
!========================================================================
      SUBROUTINE DDTINI( DT1, IMPLCT )
!========================================================================
! DDTINI sets up coefficients for time integration using routine
! DDTPGQ. DT1 is the time-step (in seconds). IMPLCT is the fractional
! implicitness (between 0.0 and 1.0) for time-stepping. The coefficients
! are stored in common module tmarch.h.
! **Note: Normally DDTPGQ uses a partially-implicit leap-frog scheme
!         for time-stepping (using PGQSP0, PGQSP1). But for start-up
!         purposes, you can set DT1 = 0.5 DT, and set PGQSP0 = PGQSP1
!         (i.e. half the actual time step) to get the effect of a
!         partially-implicit forward scheme.

      IMPLICIT NONE

      REAL DT1, IMPLCT
      REAL RTEMP, RTEMP2, TEMP
      REAL TMPMAT(L1MAX,L1MAX), TMPMA2(L1MAX,L1MAX)
      INTEGER L, L0, N
!========================================================================

!     Store time-step
      DT= DT1
!     Store implicitness factor
      IMPFAC= IMPLCT
!     Fractional implicit time-step
      DTFAC= IMPFAC * 2.0 * DT

!     Compute effective ( (-(Del**2)*D) -> (d**2)D/dt**2) matrix
!     TMPMAT is -CQ OF H-W
      CALL MATMUL( TMPMAT, L1MAX, T2GPCP, L1MAX, D2TT, L1MAX, &
                   L1, L1, L1)

!     Initialize useful wave-no dependent quantities
      DO 30 N=0,N1

!       Del**8 damping correction factor
        TD8COR(N) =  1.0 / (1.0 + DTFAC * TD8FAC(N))
        UD8COR(N) =  1.0 / (1.0 + DTFAC * UD8FAC(N))
!
!       Compute implicit correction matrix
        RTEMP = (A0QINV*N*(N+1)) * TD8COR(N) * UD8COR(N) * DTFAC**2
        RTEMP2 = RTEMP * RGAS

        DO 20 L=1,L1
           TEMP = RTEMP2 * TREFLV(L)

           DO 10 L0=1,L1
            TMPMA2(L,L0)= -RTEMP * TMPMAT(L,L0) + TEMP * DP(L0)
   10       CONTINUE

!         Add identity matrix
          TMPMA2(L,L)= TMPMA2(L,L) + 1.0
   20     CONTINUE

!       Invert matrix
        CALL MATINV( IMPCOR(1,1,N), L1MAX, TMPMA2, L1MAX, L1 )
   30   CONTINUE
 
      RETURN
      END
 











































!========================================================================
      SUBROUTINE DDTPGQ( DPGQSP, DPSSP, PGQSP0, PGQSP1, PSSP0, PSSP1 )
!========================================================================
! DDTPGQ calculates DPGQSP, the "adiabatic" time tendency of prognostic
! quantities in spectral space at some model time. (i.e. The only
! damping effect included in these calculations is the scale-selective
! Del**8 damping. Any other damping effects will need to be evaluated
! separately)
!
! Note: The coefficients required for time stepping must be initialized
!        by calling DDTINI before the first call to DDTPGQ.

! NOTE:  The -Laplacian of the surface geopotential must be added
!        in a user supplied forcing routine to the divergence tendency.

! Input:  COMPLEX PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX) -- spectral
!         representation of prognostic quantity at (current time) - DT1
!         COMPLEX PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX) -- spectral
!         representation of prognostic quantity at current model time
!         COMPLEX PSSP0(0:M1MAX,0:N1MAX) -- spectral
!         representation of ln(surface pressure) at previous model time
!         COMPLEX PSSSP1(0:M1MAX,0:N1MAX) -- spectral
!         representation of ln(surface pressure) at current model time
! Output: COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX) -- spectral
!                   "adiabatic" time-tendency of prognostic quantities
!         COMPLEX DPSSP(0:M1MAX,0:N1MAX) -- spectral
!                   "adiabatic" time-tendency of ln(surface pressure)
!
      IMPLICIT NONE

      COMPLEX PSSP0(0:M1MAX,0:N1MAX), PSSP1(0:M1MAX,0:N1MAX)
      COMPLEX DPSSP(0:M1MAX,0:N1MAX)
      COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      COMPLEX PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      COMPLEX PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      REAL U(K1MAX, K2MAX, JX:JY, L1MAX), TP(K1MAX, K2MAX, 0:1)
      REAL UCAP(K1MAX, K2MAX, JX:JY), TPCAP(K1MAX, K2MAX)
      REAL VUFLX(K1MAX, K2MAX, JX:JY, 0:1), VTPFLX(K1MAX, K2MAX, 0:1)
      REAL W(K1MAX, K2MAX, 0:1)
      REAL ABSVOR(K1MAX, K2MAX), DIV(K1MAX, K2MAX, L1MAX)
      REAL DELQ(K1MAX, K2MAX, JX:JY, L1MAX)
      REAL VDQ(K1MAX,K2MAX,L1MAX)
      REAL DBAR(K1MAX,K2MAX), GBAR(K1MAX,K2MAX)
      REAL FLXPH(K1MAX, K2MAX, JX:JY), FPH(K1MAX, K2MAX), TEMP, RT
      COMPLEX FSP1(0:M1MAX,0:N1MAX), FSP2(0:M1MAX,0:N1MAX)
      INTEGER LMHALF, LPZERO, LPHALF, LPONE
!     (Mnemonic for L-(1/2), L+0, L+(1/2), L+1)
      INTEGER I, J, K, L, L0, M, N, IPGQ
      REAL TRC(K1MAX, K2MAX, 0:1), VTRCFX(K1MAX, K2MAX, 0:1)
      REAL TRCCAP(K1MAX, K2MAX)
      REAL DAYINV, SURFD
      
      COMPLEX SURF1(0:100,0:100)
      COMMON /SURFACE/ SURF1
      
! ----------------------------------------------------------------------
!     Initialize all spectral time-tendencies to zero
      DO 10 L=1,L1
          CALL ZEROSP( DPSSP(0,0), N1)
        DO 10 IPGQ=1,3+NTRACE
          CALL ZEROSP( DPGQSP(0,0,IPGQ,L), N1)
   10     CONTINUE

      DO 15 K=1,K2
      DO 15 J=1,K1
         DBAR(J,K) = 0.
 15      GBAR(J,K) = 0.
! ----------------------------------------------------------------------
! **Note Values of prognostic quantities at current model time
!        (PGQSP1) are used in the following calculations. i.e. The
!        treatment is leap-frog in time.  These are the non-linear 
!        flux terms.
! ----------------------------------------------------------------------
!     Integrate continuity equation to find surface pressure tendency
!     (q=ln(p_surf))

      DO 50 L=1,L1

!       Compute streamfunction PSI and velocity potential CHI
           CALL IDELSQ( FSP1, PGQSP1(0,0,JVOR,L), N1)
           CALL IDELSQ( FSP2, PGQSP1(0,0,JDIV,L), N1)

!       Then compute horizontal velocity from PSI and CHI
           CALL HVELOC( U(1,1,JX,L), U(1,1,JY,L), FSP1, FSP2 )

!     Similarly, find Delq
           CALL GRAD (DELQ(1,1,JX,L), DELQ(1,1,JY,L), PSSP1(0,0))

!     Now find V.Delq and Div and accumulate in -GBAR and DBAR
           DO 20 K=1,K2
           DO 20 J=1,K1
              VDQ(J,K,L)=U(J,K,JX,L)*DELQ(J,K,JX,L)+ &
                         U(J,K,JY,L)*DELQ(J,K,JY,L)
 20           GBAR(J,K) = GBAR(J,K) - VDQ(J,K,L) * DP(L)


           CALL PHYSIC( DIV(1,1,L), PGQSP1(0,0,JDIV,L) )

           DO 25 K=1,K2
           DO 25 J=1,K1
 25           DBAR(J,K) = DBAR(J,K) + DIV(J,K,L) * DP(L)
 50   CONTINUE

      CALL SPECTR (DPSSP, GBAR)

! ----------------------------------------------------------------------
!     Start main altitude loop
      LMHALF= 0
      LPZERO= 0

      DO 300 L=0,L1
        LMHALF= MOD(LMHALF+1, 2)
        LPHALF= MOD(LMHALF+1, 2)
        LPZERO= MOD(LPZERO+1, 2)
        LPONE=  MOD(LPZERO+1, 2)

        IF (L.LT.L1) THEN
! ----------------------------------------------------------------------
!       Compute "perturbation" temp at level L+1
           CALL SPCOPY( FSP1, PGQSP1(0,0,JPOT,L+1), N1 )
           FSP1(0,0)= FSP1(0,0) - TREFLV(L+1)
           CALL PHYSIC( TP(1,1,LPONE), FSP1 )

!       Compute tracer concentration in physical space (at level L+1)
!           IF (NTRACE .GT. 0) THEN
!              CALL PHYSIC( TRC(1,1,LPONE), PGQSP1(0,0,JTR1,L+1) )
!           ENDIF
! ----------------------------------------------------------------------
        ENDIF
!
!       Compute omega and vertical fluxes at level L+1/2
        IF ( (L.EQ.0) .OR. (L.EQ.L1) ) THEN
! ----------------------------------------------------------------------
!         Omega is zero at the boundary half-levels 1/2 and L1+1/2
!         => All vertical fluxes are also zero
          DO 110 K=1,K2
            DO 110 J=1,K1
              W(J,K,LPHALF)= 0.0
              VTPFLX(J,K,LPHALF)= 0.0

	      IF (NTRACE .GT. 0) VTRCFX(J,K,LPHALF)= 0.0
  110         CONTINUE

          DO 120 I=JX,JY
            DO 120 K=1,K2
              DO 120 J=1,K1
  120           VUFLX(J,K,I,LPHALF)= 0.0
! ----------------------------------------------------------------------
        ELSE
! ----------------------------------------------------------------------
!         L+1/2 is not at the boundaries. Compute omega and vertical
!         fluxes at level L+1/2.
          DO 130 I=JX,JY
            DO 130 K=1,K2
              DO 130 J=1,K1
                UCAP(J,K,I)= 0.5 * (U(J,K,I,L+1) - U(J,K,I,L))
  130           CONTINUE

          DO 140 K=1,K2
            DO 140 J=1,K1
              TPCAP(J,K) = 0.5 * ( TP(J,K,LPONE)/PKLV(L+1) &
                                 - TP(J,K,LPZERO)/PKLV(L))

	      IF (NTRACE .GT. 0) TRCCAP(J,K)=0.5 * ( TRC(J,K,LPONE) - TRC(J,K,LPZERO))
  140         CONTINUE

!         Compute "pressure velocity" omega at L+1/2
          DO 150 K=1,K2
            DO 150 J=1,K1
              W(J,K,LPHALF)= W(J,K,LMHALF) - DP(L) * &
                    (DIV(J,K,L) + VDQ(J,K,L) - DBAR(J,K) + GBAR(J,K))
  150         CONTINUE

!         Compute vertical fluxes at level L+1/2
          DO 160 I=JX,JY
            DO 160 K=1,K2
              DO 160 J=1,K1
                VUFLX(J,K,I,LPHALF) = -W(J,K,LPHALF) * UCAP(J,K,I)
  160           CONTINUE

          DO 170 K=1,K2
            DO 170 J=1,K1
              VTPFLX(J,K,LPHALF) = -W(J,K,LPHALF) * TPCAP(J,K)

	      IF (NTRACE .GT. 0) VTRCFX(J,K,LPHALF) = -W(J,K,LPHALF) * TRCCAP(J,K)
  170         CONTINUE
! ----------------------------------------------------------------------
        ENDIF
!
        IF (L.GT.0) THEN
! ----------------------------------------------------------------------
!     Compute vorticity/divergence tendencies at level 
!     'Linear' contribution from geopotential added to D tendency later
! ----------------------------------------------------------------------
!
!         Absolute vorticity at level L
          CALL SPCOPY( FSP1, PGQSP1(0,0,JVOR,L), N1 )

!         Add planetary vorticity to relative vorticity
          FSP1(0,1)= FSP1(0,1) + CMPLX(FSP01, 0.)
          CALL PHYSIC( ABSVOR, FSP1 )

!         Horizontal "flux convergence" of momentum at level L
          DO 210 K=1,K2
            DO 210 J=1,K1
              TEMP =  RGAS * TP(J,K,LPZERO)
              FLXPH(J,K,JX)= + ABSVOR(J,K) * U(J,K,JY,L) &
                    - TEMP * DELQ(J,K,JX,L)
              FLXPH(J,K,JY)= - ABSVOR(J,K) * U(J,K,JX,L) &
                    - TEMP * DELQ(J,K,JY,L)
  210         CONTINUE

!         Add contribution from vertical flux divergence
          TEMP= 1.0 / DP(L)
          DO 220 I=JX,JY
            DO 220 K=1,K2
              DO 220 J=1,K1
                FLXPH(J,K,I)= FLXPH(J,K,I) + &
                     TEMP * (VUFLX(J,K,I,LPHALF)+VUFLX(J,K,I,LMHALF))
  220           CONTINUE

!         Vorticity/divergence tendencies from "flux conv."
          CALL CURLZ( DPGQSP(0,0,JVOR,L), FLXPH(1,1,JX), FLXPH(1,1,JY))
          CALL DIVERG(DPGQSP(0,0,JDIV,L), FLXPH(1,1,JX), FLXPH(1,1,JY))
!
!         Compute kinetic energy at level L
          DO 230 K=1,K2
            DO 230 J=1,K1
              FPH(J,K)= 0.5*(U(J,K,JX,L)**2 + U(J,K,JY,L)**2)
  230         CONTINUE

!         Convert to spectral space
          CALL ZEROSP( FSP1, N1 )
          CALL SPECTR( FSP1, FPH )

!         Compute the Laplacian, and add contribution to D tendency
          CALL DELSQ( FSP1, FSP1, N1 )
          DO 240 N=0,N1
            DO 240 M=0,MIN0(N,M1)
              DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) - FSP1(M,N)
  240         CONTINUE

! ----------------------------------------------------------------------
!         Compute temp tendency at level L
! ----------------------------------------------------------------------
!         Compute vertical "perturbation" temp flux convergence
          TEMP= PKLV(L) / DP(L)
          DO 250 K=1,K2
            DO 250 J=1,K1
              FPH(J,K)= TEMP * (VTPFLX(J,K,LPHALF)+VTPFLX(J,K,LMHALF))
  250         CONTINUE

!         Add horizontal "flux" of "perturbation" temp
          CALL GRAD (FLXPH(1,1,JX), FLXPH(1,1,JY), PGQSP1(0,0,JPOT,L))
           
          DO 260 I=JX,JY
            DO 260 K=1,K2
              DO 260 J=1,K1
                FPH(J,K)= FPH(J,K)- FLXPH(J,K,I) * U(J,K,I,L)
  260           CONTINUE

!     Surface pressure tendency terms. Recall NNT2DT = -YZ OF H-W
          DO 265 K=1,K2
             DO 265 J=1,K1
                TEMP = TREFLV(L)+ TP(J,K,LPZERO)
                TEMP = TEMP * (VDQ(J,K,L)+GBAR(J,K))
                TEMP = TEMP - TP(J,K,LPZERO) * DBAR(J,K)
                FPH(J,K) = FPH(J,K) + KAPPA*TEMP

                DO 265 L0=1,L1
 265               FPH(J,K) = FPH(J,K) + NNT2DT(L,L0) * VDQ(J,K,L0)

!         Convert to spectral space and add to temp tendency
          CALL SPECTR( DPGQSP(0,0,JPOT,L), FPH )
!
!     All adiabatic tendency terms, excluding the geopotential in the
!     D tendency equation and the linearized heating term in the temp
!     tendency equation, have now been computed for level L.
!
! ----------------------------------------------------------------------
!         Compute tracer time tendency at level L
! ----------------------------------------------------------------------
!
!	  IF (NTRACE .GT. 0) THEN

!         Compute vertical tracer flux convergence
!            TEMP= 1.0 / DP(L)
!            DO 270 K=1,K2
!              DO 270 J=1,K1
!                FPH(J,K)= TEMP * (VTRCFX(J,K,LPHALF)+VTRCFX(J,K,LMHALF))
!  270           CONTINUE
!
!         Add horizontal "flux" of tracer
!          CALL GRAD (FLXPH(1,1,JX), FLXPH(1,1,JY), PGQSP1(0,0,JTR1,L))
!
!            DO 280 I=JX,JY
!              DO 280 K=1,K2
!                DO 280 J=1,K1
!                  FPH(J,K)= FPH(J,K) - FLXPH(J,K,I) * U(J,K,I,L)
!  280             CONTINUE
!
!         Convert to spectral space and add to temp tendency
!            CALL SPECTR( DPGQSP(0,0,JTR1,L), FPH )
!	  ENDIF
! ----------------------------------------------------------------------
        ENDIF
!
  300   CONTINUE

!
!     Now we compute the remaining terms of the tendency equation,
!     including the damping terms. (These are the linear terms)
!     **Note: Values of prognostic quantities at model time TIM - DT1
!             (PGQSP0) are used in the following calculations. i.e. The
!             treatment is explicit in time.
!
      DO 500 L=1,L1
        RT = RGAS * TREFLV(L)

        DO 500 N=0,N1
! ----------------------------------------------------------------------
!         Add geopotential and q contribution to D tendency
          TEMP= A0QINV*N*(N+1)

          DO 410 L0=1,L1
            DO 410 M=0,MIN0(N,M1)
              DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) &
                    + TEMP * (T2GPCP(L,L0) * PGQSP0(M,N,JPOT,L0))
  410         CONTINUE

          DO 415 M=0,MIN0(N,M1)
  415         DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) + &
                  RT * TEMP * PSSP0(M,N)

!
!         Add contribution from reference static stability to temp
!         tendency (Recall D2TT = - Q of H-W)
          DO 420 L0=1,L1
            DO 420 M=0,MIN0(N,M1)
              DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L) &
                   + D2TT(L,L0) * PGQSP0(M,N,JDIV,L0)
  420         CONTINUE
!
!         Add Del**8 damping terms to VOR/DIV/POT tendencies
          DO 430 M=0,MIN0(N,M1)
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L) &
                  - UD8FAC(N) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) & 
                  - UD8FAC(N) * PGQSP0(M,N,JDIV,L)
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L) &
                  - TD8FAC(N) * PGQSP0(M,N,JPOT,L)
  430       CONTINUE

!	  IF (NTRACE .GT. 0) THEN
!
!         Add Del**8 damping term to tracer tendency
!            DO 440 M=0,MIN0(N,M1)
!              DPGQSP(M,N,JTR1,L)= DPGQSP(M,N,JTR1,L)
!     1            - TD8FAC(N) * PGQSP0(M,N,JTR1,L)
!  440         CONTINUE
!	  ENDIF

  500     CONTINUE

!     Add divergence contribution
!     and Del**8 damping to surface pressure tendency
      DO N=0,N1
         DO L=1,L1
            DO M=0,MIN0(N,M1)
               DPSSP(M,N)=DPSSP(M,N)-DP(L)*PGQSP0(M,N,JDIV,L)
            ENDDO
         ENDDO

         DO M=0,MIN0(N,M1)
            DPSSP(M,N)=DPSSP(M,N)-TD8FAC(N)*PSSP0(M,N)
         ENDDO
      ENDDO
      
!      DAYINV= 1.0 / 86400.
!      SURFD= 0.0 * DAYINV
!      
!     nudge surface temperature to assimilated state
!      DO N=0,N1
!        DO M=0,M1
!          DPSSP(M,N)= DPSSP(M,N) 
!     1      + SURFD * (SURF1(M,N) - PSSP0(M,N))
!        ENDDO
!      ENDDO
           
      RETURN
      END
 
































































!========================================================================
      SUBROUTINE GRAD(UPH,VPH,CHISP)
!========================================================================
!     GRAD  calculates the gradient dCHI/dx, dCHI/dy (UPH, VPH)  in physical
!     space given the spectral space representation of a function CHISP.
!     Input:  COMPLEX CHISP(0:M1MAX,0:N1MAX) -- spectral function
!     Output: REAL UPH(K1MAX,K2MAX) -- zonal component DCHI_DX
!             REAL VPH(K1MAX,K2MAX) -- meridional component DCHI_DY
!
!     This is the 'velocity potential' part of HVELOC, stripped down for
!     efficiency.

      IMPLICIT NONE

      REAL UPH(K1MAX,K2MAX), VPH(K1MAX,K2MAX)
      COMPLEX CHISP(0:M1MAX,0:N1MAX)
      COMPLEX ULA(0:M1MAX,K2MAX), VLA(0:M1MAX,K2MAX)
      COMPLEX FTMK1(0:M1MAX)
      REAL RTEMP1, RTEMP2
      INTEGER K, M

!     Calculate velocity components (ULA,VLA) in latitude space
      DO 10 K=1,K2
        RTEMP1= A0INV * COSINV(K)
        RTEMP2= A0INV * COSPHI(K)

         CALL XSP2LA(FTMK1,CHISP,P2(0,0,K),H2(0,0,K),N1,K)

         DO 11 M=0,M1
 11         ULA(M,K)= RTEMP1*FTMK1(M)

         CALL YSP2LA(FTMK1,CHISP,P2(0,0,K),H2(0,0,K),N1,K)

         DO 10 M=0,M1
 10         VLA(M,K)= RTEMP2*FTMK1(M)

!     Convert to physical space
      CALL L2P(UPH,ULA)
      CALL L2P(VPH,VLA)
      
      RETURN
      END
 
!========================================================================
      SUBROUTINE DDTFRC( DPGQSP, PGQSP0, TIM, RAYFAC, EDDIF, ZS)
!========================================================================
!     DDTFRC adds forcing/damping terms to the spectral time tendency DPGQSP
!     (A call to DDTFRC, if any, should immediately follow the call to
!     DDTPGQ that computes DPGQSP)
!
!     Input:  COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
!             "adiabatic" spectral time-tendency of prognostic quantities
!             COMPLEX PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
!             spectral representation of prognostic quantity at all
!             levels at (current model time) - DT1
!             COMPLEX PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
!             spectral representation of prognostic quantity at all
!             levels at current model time
!     Output: COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
!             spectral time-tendency, with effects of forcing/damping added
!
!     **Note: Values of prognostic quantities at current model time - DT1
!             (PGQSP0) are used in the following calculations. i.e., the
!             treatment is explicit in time. Numerical instabilities otherwise.

      IMPLICIT NONE

      COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      COMPLEX PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX)

      REAL TEMP
      COMPLEX STRESS(0:M1MAX, 0:N1MAX, 0:1)
      INTEGER L, M, N, IPGQ
      INTEGER LMHALF, LPHALF ! Mnemonic for L-(1/2), L+(1/2)
      INTEGER LMONE, LZERO

      REAL TEMP0(K2MAX,L1MAX), U0(K2MAX,L1MAX)
      COMPLEX VOR0(0:M1MAX,0:N1MAX,L1MAX)
      COMMON /ZONE/ TEMP0, U0, VOR0

      REAL WDAMP(L1MAX), MDAMP(L1MAX), TDIFF(L1MAX)
      REAL ALPHAT(L1MAX), ALPHAW(L1MAX)
      INTEGER SDTRUNC, DAMPTRUNC
      COMMON /WAVD/ WDAMP, MDAMP, TDIFF, ALPHAT, ALPHAW, SDTRUNC, DAMPTRUNC
      REAL DAYINV, TEMPAT, TEMPAW
      
      REAL DAMPFAC
      
      REAL TIME0
      COMMON /htim/ TIME0
      
      REAL EDDIF  
      REAL*4 ZS(L1MAX)
      COMPLEX DAMPU(L1MAX), DAMPT(L1MAX)
      
      REAL WCLAMP1, WCLAMP2, WCLAMP3
      INTEGER TROPLVL, STRATLVL
      COMMON /WAVECLAMP/ WCLAMP1, WCLAMP2, WCLAMP3, TROPLVL, STRATLVL   
      
      REAL ANDREWS, Z1, Z2, ETA0, A1, A2, Z
      REAL KZZTIM, B1, B2, B3, B4, B5, B6, B7, B8, B9, KZZPHASE
      REAL RHO, MOL, CONVFAC, GAUSS
      
      REAL TIM, RAYFAC ! dummies to make driver compatible with other versions of DDTFRC
      
!======================================================================================================================
!     Nudging to assimilated files, Newtonian cooling, and Rayleigh friction damping terms. 
!======================================================================================================================
      DAYINV=1.0/86400.
            
      DO L=1,L1      
        Z = ZS(L)
        
        CONVFAC = (2*3.14159265359/25000)**2
        
        RHO=1.34*EXP(-Z/7.)
        MOL=5.28e-13*EXP(Z/7)
        
        IF (Z .GT. 50) THEN
          MOL= MIN(MOL,10*DAYINV) * 0
        ELSE
          MOL=0
        ENDIF
                
        Z1= 13.1 * 7.0
        Z2= 14.3 * 7.0
        A1= 1.1 * 7.0   ! 4.0 * 7.0
        A2= 1.113 * 7.0 ! 2.0 * 7.0

        Z1= 85
        Z2= 100
        A1= 10   ! 4.0 * 7.0
        A2= 1   ! 2.0 * 7.0
        
        B1 = 4.06e-6
        B2 = -8.77e-7
        B3 = -2.28e-6
        B4 = 1.77e-6
        B5 = 2.15e-6
        B6 = -3.05e-7
        B7 = -2.66e-7
        B8 = 4.08e-7
        B9 = 1.59e-7 
                
        KZZTIM = (TIM/86400. - 80. - 31 + 30) ! start december 2012, 80d spin-up
        KZZPHASE = KZZTIM * 1 / 365 * 2 * 3.14159265359
                
        ETA0= (B1 + B2*SIN(KZZPHASE) + B3*COS(KZZPHASE) + B4*SIN(2*KZZPHASE)  &
               + B5*COS(2*KZZPHASE) + B6*SIN(3*KZZPHASE) + B7*COS(3*KZZPHASE) &
               + B8*SIN(4*KZZPHASE) + B9*COS(4*KZZPHASE))
             
        ETA0=ETA0*5000**2 
        
!        GAUSS = 100*EXP(-(KZZPHASE - (45. / 365. * 2 * 3.14159265359))**2 / (0.25**2))
!        ETA0 = ETA0 - GAUSS
        
        IF (ETA0 .LT. 0.) ETA0 = 0.
                
        IF (Z .LE. Z1) THEN
          ANDREWS = ETA0*EXP(-((Z - Z1)/A1)**2)
        ELSE
          IF (Z .LE. Z2) THEN
            ANDREWS = ETA0
          ELSE
            ANDREWS = ETA0*EXP(-((Z - Z2)/A2)**2) 
          ENDIF
        ENDIF
        
        IF (Z .GT. 20) THEN
          DAMPU(L)=COMPLEX(ANDREWS*CONVFAC*EDDIF + MOL,0)
          DAMPT(L)=COMPLEX(REAL(TRLXLV(L))+ANDREWS*CONVFAC*EDDIF+MOL,0)
        ELSE
          DAMPU(L)=URLXLV(L)
          DAMPT(L)=TRLXLV(L)
        ENDIF
                 
!        WRITE(6,'(9f8.2)') Z, DAMPU(L)/DAYINV, DAMPT(L)/DAYINV
      ENDDO
            
!     MLT nudging 
      DO L=1,STRATLVL-1
        DO N=0,N1
          DO M=0,DAMPTRUNC
            IF (M .GT. 0) THEN 
              TEMPAT = MAX(WCLAMP3*DAYINV,REAL(DAMPT(L)))
              TEMPAW = MAX(WCLAMP3*DAYINV,REAL(DAMPU(L)))
            ELSE
              TEMPAT = ALPHAT(L)
              TEMPAW = ALPHAW(L)
            ENDIF

            IF (M .GT. SDTRUNC) THEN
               DAMPFAC = 0.
            ELSE
               DAMPFAC = 1.
            ENDIF
                          
!           Nudge temperature to assimilated state with rate ALPHAT
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)   &
                  + TEMPAT * TMNLV(M,N,L) * DAMPFAC  &
                  - TEMPAT * PGQSP0(M,N,JPOT,L)
     
!           Nudge wind to assimilated state with rate ALPHAW
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)   &
                  + TEMPAW * VOR0(M,N,L) * DAMPFAC   &
                  - TEMPAW * PGQSP0(M,N,JVOR,L)
     
            IF (M .EQ. 0) THEN
!             damp mean state divergence field to zero with rate mdamp 
              DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) &
                  - MDAMP(L) * PGQSP0(M,N,JDIV,L)
            ENDIF
            
          ENDDO

          DO M=1,MIN0(N,M1)
!           Rayleigh friction
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L) &
                  - URLXLV(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) &
                  - URLXLV(L) * PGQSP0(M,N,JDIV,L)
     
!           Newtonian cooling
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L) &
                  - TRLXLV(L) * PGQSP0(M,N,JPOT,L)

!           Sponge layer damp vorticity & divergence
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L) &
                  - WDAMP(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) &
                  - WDAMP(L) * PGQSP0(M,N,JDIV,L)
          ENDDO
        ENDDO
      ENDDO

!     stratosphere nudging 
      DO L=STRATLVL,TROPLVL-1
        DO N=0,N1
          DO M=0,DAMPTRUNC
            IF (M .GT. 0) THEN 
              TEMPAT = MAX(WCLAMP2*DAYINV,REAL(DAMPT(L)))
              TEMPAW = MAX(WCLAMP2*DAYINV,REAL(DAMPU(L)))
            ELSE
              TEMPAT = ALPHAT(L)
              TEMPAW = ALPHAW(L)
            ENDIF
            
            IF (M .GT. SDTRUNC) THEN
               DAMPFAC = 0.
            ELSE
               DAMPFAC = 1.
            ENDIF
                          
!           Nudge temperature to assimilated state with rate ALPHAT
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)  &
                  + TEMPAT * TMNLV(M,N,L) * DAMPFAC &
                  - TEMPAT * PGQSP0(M,N,JPOT,L)
     
!           Nudge wind to assimilated state with rate ALPHAW
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)  &
                  + TEMPAW * VOR0(M,N,L) * DAMPFAC  &
                  - TEMPAW * PGQSP0(M,N,JVOR,L)
     
            IF (M .EQ. 0) THEN
!             damp mean state divergence field to zero with rate mdamp 
              DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)&
                - MDAMP(L) * PGQSP0(M,N,JDIV,L)
            ENDIF
            
          ENDDO

          DO M=1,MIN0(N,M1)
!           Rayleigh friction
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L) &
                  - URLXLV(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) &
                  - URLXLV(L) * PGQSP0(M,N,JDIV,L)
     
!           Newtonian cooling
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L) &
                  - TRLXLV(L) * PGQSP0(M,N,JPOT,L)

!           Sponge layer damp vorticity & divergence
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L) &
                  - WDAMP(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) &
                  - WDAMP(L) * PGQSP0(M,N,JDIV,L)
          ENDDO
        ENDDO
      ENDDO
      
!     troposphere nudging 
      DO L=TROPLVL,L1
        DO N=0,N1
          DO M=0,DAMPTRUNC
            IF (M .GT. 0) THEN 
              TEMPAT = MAX(WCLAMP1*DAYINV,REAL(DAMPT(L)))
              TEMPAW = MAX(WCLAMP1*DAYINV,REAL(DAMPU(L)))
            ELSE
              TEMPAT = ALPHAT(L)
              TEMPAW = ALPHAW(L)
            ENDIF
                    
            IF (M .GT. SDTRUNC) THEN
               DAMPFAC = 0.
            ELSE
               DAMPFAC = 1.
            ENDIF
            
!           Nudge temperature to assimilated state with rate ALPHAT
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)  &
                  + TEMPAT * TMNLV(M,N,L) * DAMPFAC &
                  - TEMPAT * PGQSP0(M,N,JPOT,L)
     
!           Nudge wind to assimilated state with rate ALPHAW
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)  &
                  + TEMPAW * VOR0(M,N,L) * DAMPFAC  &
                  - TEMPAW * PGQSP0(M,N,JVOR,L)
     
            IF (M .EQ. 0) THEN
!             damp mean state divergence field to zero with rate mdamp 
              DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) &
                  - MDAMP(L) * PGQSP0(M,N,JDIV,L)
            ENDIF
            
          ENDDO

          DO M=1,MIN0(N,M1)
!           Rayleigh friction
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L) &
                  - URLXLV(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) &
                  - URLXLV(L) * PGQSP0(M,N,JDIV,L)
     
!           Newtonian cooling
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L) &
                  - TRLXLV(L) * PGQSP0(M,N,JPOT,L)

!           Sponge layer damp vorticity & divergence
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L) &
                  - WDAMP(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) &
                  - WDAMP(L) * PGQSP0(M,N,JDIV,L)
          ENDDO
        ENDDO
      ENDDO

! ---------------------------------------------------------------
!     Compute vertical viscous damping (for vorticity/temperature)

      DO 800 IPGQ=JVOR,JPOT
        LMHALF= 1
        LPHALF= 0
        LZERO= 0
        LMONE= 1

        IF (IPGQ .EQ. JPOT) THEN
           IF (TDIFF(1) .LT. 0) GO TO 500
        ELSE
           IF (VVISC(1) .LT. 0) GO TO 500
        ENDIF

!       No stress at uppermost level
        CALL ZEROSP( STRESS(0,0,LPHALF), N1 )

        DO 350 L=1,L1
          LMHALF= MOD(LMHALF+1, 2)
          LPHALF= MOD(LMHALF+1, 2)

          IF (L.LT.L1) THEN
!            Compute stress at level L+1/2
!            SUBTRACT STRESS OF INITIAL STATE
             IF (IPGQ .EQ. JPOT) THEN
               TEMP = ABS(TDIFF(L))  / (HDP(L+1) + HDP(L)) ! ABS(TDIFF)
               DO 309 N=0,N1
                 DO 309 M=0,MIN0(N,M1)
!                    STRESS(M,N,LPHALF)= TEMP *
!     $                   (PGQSP0(M,N,IPGQ,L+1)/PKLV(L+1) -
!     $                   PGQSP0(M,N,IPGQ,L)/PKLV(L)) 
                   STRESS(M,N,LPHALF)= TEMP *     &
                         (PGQSP0(M,N,IPGQ,L+1) -  &
                          PGQSP0(M,N,IPGQ,L))
 309           CONTINUE

		  DO N=0,N1
		    STRESS(0,N,LPHALF) = STRESS(0,N,LPHALF) - TEMP * &
              (REAL(TMNLV(0,N,L+1)) - REAL(TMNLV(0,N,L)))
		  ENDDO
             ELSE
                TEMP= ABS(VVISC(L)) / (HDP(L+1) + HDP(L))
                DO 310 N=0,N1
                  DO 310 M=0,MIN0(N,M1)
                    STRESS(M,N,LPHALF)= TEMP * &
                         (PGQSP0(M,N,IPGQ,L+1) - PGQSP0(M,N,IPGQ,L))
 310              CONTINUE

		  IF (IPGQ .EQ. JVOR) THEN
		    DO N=0,N1
		      STRESS(0,N,LPHALF) = STRESS(0,N,LPHALF) - TEMP * &
                (REAL(VOR0(0,N,L+1)) - REAL(VOR0(0,N,L)))
		    ENDDO
		  ENDIF
             ENDIF

          ELSE
!           Compute stress at surface (level L1+1/2)
            DO 320 N=0,N1
              DO 320 M=0,MIN0(N,M1)
                 IF(IPGQ .EQ. JPOT) THEN
                    CALL ZEROSP( STRESS(0,0,LPHALF), N1 )
                 ELSE
                    STRESS(M,N,LPHALF) = -USDRAG * PGQSP0(M,N,IPGQ,L)
                 ENDIF
  320           CONTINUE
          ENDIF

!         Add contribution from viscous stress divergence
          DO 330 N=0,N1
            DO 330 M=0,MIN0(N,M1)
              IF(IPGQ .EQ. JPOT) THEN
!                DPGQSP(M,N,IPGQ,L)= DPGQSP(M,N,IPGQ,L) + PKLV(L) * 
!     1               (STRESS(M,N,LPHALF) - STRESS(M,N,LMHALF)) / DP(L)
                DPGQSP(M,N,IPGQ,L)= DPGQSP(M,N,IPGQ,L) +  &
                     (STRESS(M,N,LPHALF) - STRESS(M,N,LMHALF)) / DP(L)
              ELSE
                DPGQSP(M,N,IPGQ,L)= DPGQSP(M,N,IPGQ,L) +  &
                     (STRESS(M,N,LPHALF) - STRESS(M,N,LMHALF)) / DP(L)
              ENDIF
 330      CONTINUE

 350        CONTINUE

        GO TO 800

 500    CONTINUE

      WRITE(6,*) 'YOU FOOL!!!!'
      STOP

 800  CONTINUE

      RETURN
      END

!========================================================================
      SUBROUTINE DDTIMP( DPGQSP, DPSSP )
!========================================================================
!     DDTIMP makes implict corrections to the spectral time tendency DPGQSP.
!     (DDTIMP should be called after the calls to DDTPGQ and DDTFRC to
!     compute the "adiabatic" time-tendency and to include the effects of
!     forcing/damping)
!     Input:  COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
!             explicit spectral time-tendency of prognostic quantities
!     Output: COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
!             implicit spectral time-tendency of prognostic quantities
!
      IMPLICIT NONE

      COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      COMPLEX DPSSP(0:M1MAX,0:N1MAX)
      INTEGER L0, L, M, N
      REAL RTEMP, RT
      COMPLEX TEM(0:M1MAX,L1MAX)

      DO 200 N=0,N1
!       Making implicit corrections to D tendency
        DO 30 L=1,L1

!         Contribution from explicit D tendency, with Del**8 correction
          DO 10 M=0,MIN0(N,M1)
   10       TEM(M,L)= UD8COR(N) * DPGQSP(M,N,JDIV,L)

!         Add contribution from explicit temp tendency
          RTEMP= UD8COR(N) * (A0QINV*N*(N+1)) * TD8COR(N) * DTFAC
          RT= RTEMP * RGAS * TREFLV(L)

          DO 20 L0=1,L1
            DO 20 M=0,MIN0(N,M1)
   20         TEM(M,L)=TEM(M,L) + RTEMP*T2GPCP(L,L0)*DPGQSP(M,N,JPOT,L0)

!         Add contribution from explicit q tendency
          DO M=0,MIN0(N,M1)
             TEM(M,L) = TEM(M,L) + RT * DPSSP(M,N)
          ENDDO

   30     CONTINUE
!
!       Multiply TEM by implicit correction matrix to find implicit D
!       tendency
        DO 60 L=1,L1
          DO 40 M=0,MIN0(N,M1)
   40       DPGQSP(M,N,JDIV,L)= (0., 0.)

          DO 50 L0=1,L1
            DO 50 M=0,MIN0(N,M1)
   50         DPGQSP(M,N,JDIV,L)=  DPGQSP(M,N,JDIV,L) &
                                 + IMPCOR(L,L0,N)*TEM(M,L0)
   60     CONTINUE
!
        DO 100 L=1,L1
!
!         Make implicit correction to temp tendency due to
!         implicit D tendency (recall D2TT = -M(D->T))
          DO 70 L0=1,L1
            DO 70 M=0,MIN0(N,M1)
              DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L) &
                    + DTFAC * D2TT(L,L0) * DPGQSP(M,N,JDIV,L0)
   70         CONTINUE
!
!         Make implicit correction due to Del**8 damping of temp
          DO 80 M=0,MIN0(N,M1)
   80       DPGQSP(M,N,JPOT,L)= TD8COR(N) * DPGQSP(M,N,JPOT,L)

!	  IF (NTRACE .GT. 0) THEN
!
!         Make implicit correction due to Del**8 damping of tracer
!            DO 90 M=0,MIN0(N,M1)
!   90         DPGQSP(M,N,JTR1,L)= TD8COR(N) * DPGQSP(M,N,JTR1,L)
!	  ENDIF

  100     CONTINUE

!         Make implicit correction to q tendency due to 
!         implicit D tendency
          DO L0=1,L1
             DO M=0,MIN0(N,M1)
                DPSSP(M,N)= DPSSP(M,N)  &
                    - DTFAC * DP(L0) * DPGQSP(M,N,JDIV,L0)
             ENDDO
          ENDDO
!
!         Make implicit correction due to Del**8 damping of q
          DO 110 M=0,MIN0(N,M1)
  110       DPSSP(M,N)= TD8COR(N) * DPSSP(M,N)

  200   CONTINUE
  
!     nudge surface temperature to assimilation state
!      DO N=0,N1
!        DO M=0,M1
!          DPSSP(M,N)= DPSSP(M,N)
!     1      + SURFD * (SURF1(M,N) - PSSP1(M,N))
!     1      - SURFD * DPSSP(M,N) * 0
!        ENDDO
!      ENDDO

      RETURN
      END
 
!========================================================================
      SUBROUTINE ROBFIL( FSP0, FSP1, FSP2, ROBFAC )
!========================================================================
!     Given the spectral states FSP0, FSP1, FSP2 at three successive
!     time-steps, ROBFIL performs a Robert-filter operation on the
!     states, with filter weight ROBFAC.
!     Input:  COMPLEX FSP0(0:M1MAX,0:N1MAX) -- spectral
!             representation of prognostic quantity at model time TIM-DT
!             COMPLEX FSP1(0:M1MAX,0:N1MAX) -- spectral
!             representation of prognostic quantity at model time TIM
!             COMPLEX FSP2(0:M1MAX,0:N1MAX) -- spectral
!             representation of prognostic quantity at model time TIM+DT
!             REAL ROBFAC -- filter weighting factor
!     Output: (Filtered versions of FSP0, FSP1, FSP2)
!
      IMPLICIT NONE

      COMPLEX FSP0(0:M1MAX,0:N1MAX)
      COMPLEX FSP1(0:M1MAX,0:N1MAX), FSP2(0:M1MAX,0:N1MAX)
      REAL ROBFAC
      INTEGER M, N
      REAL TEMP

      TEMP= 1.0 - 2.0 * ROBFAC
      DO 10 N=0,N1
        DO 10 M=0,MIN0(N,M1)
          FSP1(M,N)=  TEMP  *  FSP1(M,N) &
                    + ROBFAC * (FSP0(M,N) + FSP2(M,N))
   10     CONTINUE
      RETURN
      END
      
!========================================================================
      SUBROUTINE MODROBFIL( FSP0, FSP1, FSP2, ROBFAC )
!========================================================================
!     Given the spectral states FSP0, FSP1, FSP2 at three successive
!     time-steps, ROBFIL performs a Robert-filter operation on the
!     states, with filter weight ROBFAC.
!     Input:  COMPLEX FSP0(0:M1MAX,0:N1MAX) -- spectral
!             representation of prognostic quantity at model time TIM-DT
!             COMPLEX FSP1(0:M1MAX,0:N1MAX) -- spectral
!             representation of prognostic quantity at model time TIM
!             COMPLEX FSP2(0:M1MAX,0:N1MAX) -- spectral
!             representation of prognostic quantity at model time TIM+DT
!             REAL ROBFAC -- filter weighting factor
!     Output: (Filtered versions of FSP0, FSP1, FSP2)
!
!     this version displaces the state at the future time slice as well 
!     as the current timeslice, based on a 2009 Paul D. Williams paper.
!
!     alpha = 0.53 would be a typical value for a semi-implicit scheme, 
!     but don't go below 0.5; will make it unconditionally unstable. 
!     alpha = 1.0 reduces this procedure to the regular Robert-filter

      IMPLICIT NONE

      COMPLEX FSP0(0:M1MAX,0:N1MAX)
      COMPLEX FSP1(0:M1MAX,0:N1MAX), FSP2(0:M1MAX,0:N1MAX)
      REAL ROBFAC, ALPHA
      INTEGER M, N
      
      COMPLEX TEMP(0:M1MAX,0:N1MAX)
      
      ALPHA= 0.53
      
      DO N=0,N1
        DO M=0,MIN0(N,M1)
!         displacement factor
          TEMP(M,N)= (ROBFAC / 2) * (FSP0(M,N) - 2 * FSP1(M,N) &
            + FSP2(M,N))
!         apply displacement to current and future time slice
          FSP1(M,N)= FSP1(M,N) + TEMP(M,N) * ALPHA
          FSP2(M,N)= FSP2(M,N) + TEMP(M,N) * (ALPHA - 1)
        ENDDO
      ENDDO
      
      RETURN
      END
    
!===============================================
	SUBROUTINE DRY_ADJUST(DPGQSP,PGQSP1,U,TP,LENSQ)
!===============================================
    ! NB: This is technically not a prognostic routine, but 
    ! kept here for simplicity, can be changed later
	IMPLICIT NONE

	! INCLUDE 'mcons.inc'
	! INCLUDE 'spcons.inc'
	! INCLUDE 'mgrid.inc'
	! INCLUDE 'spgrid.inc'

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
!	PARAMETER (LENSQ=900.)	!MIXING LENGTH SQUARED (M^2)

!...HDELP(L), RMIN(L), G0SQ HAVE BEEN ADDED TO COMMON BLOCKS
!=========================================================

	LM=1
	LP=0

!...No stress at uppermost half-level
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

!...No stress at surface (level L1+1/2)
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

!...RICHARDSON NUMBER AND DIFFUSION COEFFICIENT
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

!...CALCULATE STRESS AT HALF-LEVELS
		RHOD = RHO * DIFF
		GR2D = GRHO * RHOD

		TSTRESS(M,N,LP) = RHOD * NSQT
		USTRESS(M,N,LP) = GR2D * DELU
		VSTRESS(M,N,LP) = GR2D * DELV
	      ENDDO
	    ENDDO
	  ENDIF

!...    Add contribution from viscous stress divergence
	  DO N=1,K2
	    DO M=1,K1
	      TTEND(M,N)= PKLV(L) * &
	 	   (TSTRESS(M,N,LP)-TSTRESS(M,N,LM))/DP(L)

	      UTEND(M,N)= &
	 	   (USTRESS(M,N,LP)-USTRESS(M,N,LM))/DP(L)

	      VTEND(M,N)= &
	 	   (VSTRESS(M,N,LP)-VSTRESS(M,N,LM))/DP(L)
	    ENDDO
	  ENDDO

	  CALL SPECTR(DPGQSP(0,0,JPOT,L),TTEND)
	  CALL CURLZ (DPGQSP(0,0,JVOR,L),UTEND,VTEND)
	  CALL DIVERG(DPGQSP(0,0,JDIV,L),UTEND,VTEND)
	ENDDO

	END


end module Prognos_sigma_mod