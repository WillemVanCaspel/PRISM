C     Prognos_sigma - version 2.3, created January, 1999
C     A mechanistic spectral primitive equation model using sigma coordinates.
C     In this version, Temperature is a prognostic variable.
C     THIS PROGRAM IS BASED ON THE MATERIAL IN HALTINER & WILLIAMS, PG
c     257-262, AND SOME COMMENTS USE THEIR NOTATION
C
C     The prognostic quantities are vorticity, divergence,
C     q=log(surface pressure) and temperature at each level.
C
C     Author: R. Saravanan, modified by D. Ortland.
C
C     This software is distributed free for research purposes. It comes with
C     absolutely no warranties whatsoever. It may be freely copied, and may also
C     be modified, provided the modified software continues to be freely
C     available under the same terms. This software may not be used for any
C     commercial purpose.

C========================================================================
      SUBROUTINE PLINI(OMEGA, R, CPRESS, GRAVIT)
C========================================================================
C     PLINI initializes the common block /PLPARM/ which contains the
C     planetary parameters.
C     Input:  REAL OMEGA -- angular velocity of planetary rotation (in 1/s)
C             REAL R -- gas constant {in J/(kg*K)}
C             REAL CPRESS -- specific heat at const. pressure {in J/(kg*K)}
C             REAL GRAVIT -- gravitational acceleration (in m/s**2)
C     Output: Common block /PLPARM/ is initialized

      IMPLICIT NONE

      INCLUDE 'mcons.inc'
      INCLUDE 'spcons.inc'
      INCLUDE 'mgrid.inc'
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
      
C========================================================================
      SUBROUTINE VERINI( NLEV, PTHICK, PHALFL, PP, DEL8DF, 
     1                   TDAMP, UDAMP, SFRIC)
C========================================================================
C
C     VERINI initializes the vertical truncation for a spectral pressure
C     coordinate model (using the Lorenz (1960) finite-differencing, with
C     vertically averaged divergence set to advection of surface pressure).
C
C     Note: Prior to calling VERINI, the horizontal truncation should have
C           been initialized by a call to routine SPINI in module sptran.F.
C           Also, the planetary parameters should have been initialized by
C           a call to routine PLINI in this module.
C
C     Note: This version of VERINI only sets up pressure grid and damping stuff.
C           Other things formerly set up in the old VERINI are now in TINI.
C           This makes it easier to use the vertical grid in initialization
C           subroutines.
C
C     NLEV is the number of layers in the vertical truncation, (NLEV >= 2).
C
C     PTHICK(NLEV) is the PRESSURE-thickness of each level.
C         Note: The upper boundary boundary is always assumed to be at zero
C         pressure.Let PSURF=sum(PTHICK).  Sigma coordinates are defined 
C         as sigma(i)-sigma(i-1) = pthick(i)/PSURF at the layer boundaries
c         (half-levels).)
c     PP(nlev) is the pressure at the full level.
c     NOTE: It is no longer assumed that this is midway between half-levels
C           DEL8DF is the del-8 (quad-harmonic) diffusion coefficient
C           (in units of A0**8/sec, where A0 is the planetary radius)
C     TDAMP(NLEV)is the sigma dependent radiative damping coefficient
C     (in units of 1/sec) that determines relaxation of temp back to
C     TRADEQ.
C     UDAMP(NLEV) is the sigma dependent Rayleigh friction coefficient
C     (in units of  1/sec) that damps velocity.
C     SFRIC is the surface rayleigh friction coefficient (in units of 1/sec)
C     that damps momentum in the lowermost level.
C
C     Lorenz, E.N. 1960: Energy and numerical weather prediction,
C     TELLUS,vol.12,pp.364-373.

      IMPLICIT NONE

      INCLUDE 'mcons.inc'
      INCLUDE 'spcons.inc'
      INCLUDE 'mgrid.inc'
      INCLUDE 'spgrid.inc'
      INCLUDE 'tmarch.inc'

      INTEGER NLEV
      REAL PTHICK(NLEV), PP(NLEV), PHALFL(0:NLEV), DEL8DF
      REAL TDAMP(NLEV), UDAMP(NLEV), SFRIC, DELZ
      INTEGER L, L0, N
 
      IF (NLEV.GT.L1MAX)
     1  CALL SPERR( 'VERINI', 'Too many sigma levels')
      L1= NLEV
 
C     Initialize pressure related parameters
C     PLV and PHLV are now sigma-coord, NOT pressure-levels
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

C     Standard vertical temp profile, and coordinate PK
      DO 20 L=1,L1
        PKLV(L)= PLV(L)**KAPPA
   20   CONTINUE

C     Compute PK-cap
      DO 30 L=1,L1-1
   30   PKCHLV(L)= 0.5 * (PKLV(L+1) - PKLV(L))
C
C     D to W conversion matrix (Z of Haltiner-Williams)
      DO 70 L=1,L1-1
        DO 70 L0=1,L1
          IF (L0.LE.L) THEN
            D2W(L,L0)= (PHLV(L)-1)*DP(L0)
          ELSE
            D2W(L,L0)= PHLV(L)*DP(L0)
          ENDIF
 70    CONTINUE

C     temp to GeoPotential conversion matrix (Hydrostatic relation)
       DO 90 L=1,L1
          DO L0=1,L1-1
             T2GPCP(L,L0)=
     1            -PHLV(L0)*CP*PKCHLV(L0)/PKLV(L0) + RGAS*DP(L0)
          ENDDO
          T2GPCP(L,L1)=RGAS*DP(L1)

          DO L0=1,L1-1
             T2GPCP(L,L0+1) = T2GPCP(L,L0+1)
     1            -PHLV(L0)*CP*PKCHLV(L0)/PKLV(L0+1)
          ENDDO
 90    CONTINUE

       DO 95 L=1, L1-1
          DO L0=L,L1-1
             T2GPCP(L,L0)   = T2GPCP(L,L0)   + CP*PKCHLV(L0)/PKLV(L0)
             T2GPCP(L,L0+1) = T2GPCP(L,L0+1) + CP*PKCHLV(L0)/PKLV(L0+1)
          ENDDO
 95    CONTINUE

C     Del**8 sub-grid-scale diffusion
      QHDAMP= DEL8DF
      DO 100 N=0,N1

C       Scale-dependent del**8 damping factor
C       (Note: QHDAMP is in units such that A0 = 1)
        TD8FAC(N) =  QHDAMP * FLOAT( N*(N+1)   )**4
        UD8FAC(N) =  QHDAMP * FLOAT( N*(N+1)-2 )**4
  100   CONTINUE
C
C     Mean thermal state, and Rayleigh/Newtonian damping coefficients
      DO 130 L=1,L1
        URLXLV(L)= UDAMP(L)
 130    TRLXLV(L)= TDAMP(L)

C     Surface drag (in units of (sigma/time) )
      USDRAG= DP(L1) * SFRIC

C     DRY ADJUSTMENT LIMIT ON RICHARDSON NUMBER
      DO L=1, L1
        DELZ = 7.*DP(L)/PLV(L)
        RMIN(L) = .25 * (1. + .1 * (10.*DELZ)**1.5)
      ENDDO
        
      RETURN
      END

C========================================================================
      SUBROUTINE TINI(TZSTD, TSTEEP, TRADEQ, VERVIS, DIFFC)
C========================================================================
C     TZSTD(NLEV) is the vertical profile of the mean temperature distribution.
C     Note: A more stably stratified version of TZSTD is used to define
C           the reference temp profile for the semi-implicit time-stepping.
C           TZSTD is also used to calculate standard height-levels
C           corresponding to the sigma levels)
c
C     TSTEEP is the steepness factor for the reference vertical temp
C     profile used for semi-implicit time-stepping. i.e. The reference
C     temp profile is steeper than TZSTD by a factor TSTEEP.
C     (One would typically choose TSTEEP > 1.0)
c
C     TRADEQ(K2MAX,NLEV) is the latitude-sigma field of the "radiative
C     equilibrium" temp distribution.
c
C     VERVIS is the vertical viscosity coefficient (in units of m**2/s) that
C     diffuses momentum.
c
C     DIFFC is thermal diffusion coefficient

      IMPLICIT NONE

      INCLUDE 'mcons.inc'
      INCLUDE 'spcons.inc'
      INCLUDE 'mgrid.inc'
      INCLUDE 'spgrid.inc'
      INCLUDE 'tmarch.inc'

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

C     Compute standard height values, using the hydrostatic relation.
C     Also initialize reference temp profile(for semi-implicit scheme).
C     The reference temp profile is simply a steeper version of the
C     standard temp profile (steeper by a factor TRFFAC).
      TRFFAC= TSTEEP
      ZSTDLV(L1)= CP*(1.0-PKLV(L1))*TSTDLV(L1)/PKLV(L1) / G0

C     Reference temp profile coincides with standard profile at surface
      TREFLV(L1)= TSTDLV(L1)
      DO 40 L=L1-1,1,-1
        ZSTDLV(L)= ZSTDLV(L+1) + CP*PKCHLV(L) *
     1            (TSTDLV(L+1)/PKLV(L+1) +TSTDLV(L)/PKLV(L)) / G0
        TREFLV(L)= TRFFAC * TSTDLV(L)
     1       + PKLV(L)/PKLV(L+1) * (TREFLV(L+1) - TRFFAC * TSTDLV(L+1))
   40   CONTINUE
C
C     W to temp tendency conversion matrix (Static stability)
C     NOTE: W2TT = - Y OF Haltiner-Williams
      DO 80 L=1,L1
        DO 80 L0=1,L1-1
          IF (L0.EQ.(L-1)) THEN
            W2TT(L,L0)= 
     1        -0.5 * (TREFLV(L)-PKLV(L)/PKLV(L-1)*TREFLV(L-1)) / DP(L)
          ELSE IF (L0.EQ.L) THEN
            W2TT(L,L0)=
     1        -0.5 * (TREFLV(L+1)*PKLV(L)/PKLV(L+1)-TREFLV(L)) / DP(L)
          ELSE
            W2TT(L,L0)= 0.0
          ENDIF
   80     CONTINUE

C     D to temp-tendency conversion matrix
C     (Composite  W2TT # D2W - kappa * T#DP)
C     NNT2DT = -YZ OF HALTINER AND WILLIAMS
C     NOTE: THIS HAS OPPOSITE SIGN TO THAT USED IN DOCUMENTATION
      CALL MATMUL(NNT2DT, L1MAX, W2TT, L1MAX, D2W, L1MAX,
     1            L1, L1-1, L1)

      DO 85 L=1, L1
         DO 85, L0=1,L1
 85         D2TT(L,L0) = NNT2DT(L,L0) - KAPPA * TREFLV(L) * DP(L0)

c.....Compute spectral coefficients of zonally symmetric mean temp
      DO 130 L=1,L1
        DO 110 K=1,K2
          DO 110 J=1,K1
  110       FPH(J,K)= TRADEQ(K,L)
        CALL ZEROSP(FSP, N1)
        CALL SPECTR(FSP, FPH)
        DO 120 N=0,N1
  120     TMNLV(0,N,L)= REAL( FSP(0,N) )
  130   CONTINUE

C     Convert viscosity coefficient from z-coordinate to SIGMA-coordinate,
C     using the standard height values
      DO 140 L=1,L1-1
        VVISC(L)= VERVIS(L) * ( (HDP(L+1)  + HDP(L)     ) /
     1                       (ZSTDLV(L) - ZSTDLV(L+1))  )**2

        IF (VERVIS(L) .LT. 0)
     $       VVISC(L)= VERVIS(L) * ( (HDP(L+1)  + HDP(L)     ) /
     1                       (ZSTDLV(L) - ZSTDLV(L+1))  )**4

        TDIFF(L)= DIFFC(L) * ( (HDP(L+1)  + HDP(L)     ) /
     1                       (ZSTDLV(L) - ZSTDLV(L+1))  )**2
C        TDIFF(L)= DIFFC(L) /PKLV(L) * ( (HDP(L+1)  + HDP(L)     ) /
C     1                       (ZSTDLV(L) - ZSTDLV(L+1))  )**2

        IF (DIFFC(L) .LT. 0)
     $       TDIFF(L)= DIFFC(L) * ( (HDP(L+1)  + HDP(L)     ) /
     1                       (ZSTDLV(L) - ZSTDLV(L+1))  )**4

  140   CONTINUE
      RETURN
      END
 
 
C========================================================================
      SUBROUTINE MATMUL(A, NA0, B, NB0, C, NC0, NROWB, NCOLB, NCOLC)
C========================================================================
C MATMUL multiplies two matrices together to produce a third matrix.
C  A = B . C
C Input:  INTEGER NA0, NB0, NC0 -- inner dimension of matrices A, B, C
C         INTEGER NROWB, NCOLB -- no. of rows and columns of matrix B
C         INTEGER NCOLC -- no. of columns of matrix C
C         REAL B(NB0,NCOLB), C(NC0,NCOLC) -- matrices to be multiplied
C Output: REAL A(NA0,NCOLC) -- matrix A = B . C (with NROWB rows)
C
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
 
C========================================================================
      SUBROUTINE MATINV (AINV, NAINV0, A, NA0, NA)
C========================================================================
C MATINV finds the inverse of a square matrix A.
C Input:  INTEGER NAINV0, NA0 -- inner dimension of matrices AINV, A
C         INTEGER NA -- no. of rows/columns in A
C (Note: NA must be <= L1MAX)
C         REAL A(NA0,NA) -- matrix A
C Output: REAL AINV(NAINV0,NA) -- inverse of A
C
      IMPLICIT NONE

      INCLUDE 'mcons.inc'

      INTEGER MAXDIM
      PARAMETER (MAXDIM=300)

      INTEGER NAINV0, NA0, NA
      REAL AINV(NAINV0,NA), A(NA0,NA)
      REAL*8 A2(MAXDIM,MAXDIM), B(MAXDIM,MAXDIM), C(MAXDIM,MAXDIM)
      REAL*8 WKSPCE(MAXDIM), AA(MAXDIM,MAXDIM), BB(MAXDIM,MAXDIM)
      INTEGER I, J, IFAIL
C========================================================================
C
      IF (NA.GT.MAXDIM) CALL SPERR( 'MATINV', 'Matrix too large')

C     Initialize temporary identity matrix
      DO 10 I=1,NA
        DO 10 J=1,NA
C         Copy matrix to correct precision
          A2(I,J)= A(I,J)
C         Create Identity matrix for RHS
          IF (I.EQ.J) THEN
            B(I,J)= 1.0
          ELSE
            B(I,J)= 0.0
          ENDIF
   10     CONTINUE
C     Call Nag routine
      IFAIL= 0
      CALL F04AEF(A2, MAXDIM, B, MAXDIM, NA, NA, C, MAXDIM,
     1            WKSPCE, AA, MAXDIM, BB, MAXDIM, IFAIL)
      IF (IFAIL.NE.0)
     1  CALL SPERR( 'MATINV', 'Unable to invert matrix')
C     Copy inverse back to original precision
      DO 20 I=1,NA
        DO 20 J=1,NA
          AINV(I,J)= C(I,J)
   20     CONTINUE

      RETURN
      END
 
C========================================================================
      SUBROUTINE DDTINI( DT1, IMPLCT )
C========================================================================
C DDTINI sets up coefficients for time integration using routine
C DDTPGQ. DT1 is the time-step (in seconds). IMPLCT is the fractional
C implicitness (between 0.0 and 1.0) for time-stepping. The coefficients
C are stored in common module tmarch.h.
C **Note: Normally DDTPGQ uses a partially-implicit leap-frog scheme
C         for time-stepping (using PGQSP0, PGQSP1). But for start-up
C         purposes, you can set DT1 = 0.5 DT, and set PGQSP0 = PGQSP1
C         (i.e. half the actual time step) to get the effect of a
C         partially-implicit forward scheme.

      IMPLICIT NONE

      INCLUDE 'mcons.inc'
      INCLUDE 'spcons.inc'
      INCLUDE 'mgrid.inc'
      INCLUDE 'spgrid.inc'
      INCLUDE 'tmarch.inc'

      REAL DT1, IMPLCT
      REAL RTEMP, RTEMP2, TEMP
      REAL TMPMAT(L1MAX,L1MAX), TMPMA2(L1MAX,L1MAX)
      INTEGER L, L0, N
C========================================================================

C     Store time-step
      DT= DT1
C     Store implicitness factor
      IMPFAC= IMPLCT
C     Fractional implicit time-step
      DTFAC= IMPFAC * 2.0 * DT

C     Compute effective ( (-(Del**2)*D) -> (d**2)D/dt**2) matrix
C     TMPMAT is -CQ OF H-W
      CALL MATMUL( TMPMAT, L1MAX, T2GPCP, L1MAX, D2TT, L1MAX,
     1             L1, L1, L1)

C     Initialize useful wave-no dependent quantities
      DO 30 N=0,N1

C       Del**8 damping correction factor
        TD8COR(N) =  1.0 / (1.0 + DTFAC * TD8FAC(N))
        UD8COR(N) =  1.0 / (1.0 + DTFAC * UD8FAC(N))
C
C       Compute implicit correction matrix
        RTEMP = (A0QINV*N*(N+1)) * TD8COR(N) * UD8COR(N) * DTFAC**2
        RTEMP2 = RTEMP * RGAS

        DO 20 L=1,L1
           TEMP = RTEMP2 * TREFLV(L)

           DO 10 L0=1,L1
            TMPMA2(L,L0)= -RTEMP * TMPMAT(L,L0) + TEMP * DP(L0)
   10       CONTINUE

C         Add identity matrix
          TMPMA2(L,L)= TMPMA2(L,L) + 1.0
   20     CONTINUE

C       Invert matrix
        CALL MATINV( IMPCOR(1,1,N), L1MAX, TMPMA2, L1MAX, L1 )
   30   CONTINUE
 
      RETURN
      END
 
C========================================================================
      SUBROUTINE DDTPGQ( DPGQSP, DPSSP, PGQSP0, PGQSP1, PSSP0, PSSP1 )
C========================================================================
C DDTPGQ calculates DPGQSP, the "adiabatic" time tendency of prognostic
C quantities in spectral space at some model time. (i.e. The only
C damping effect included in these calculations is the scale-selective
C Del**8 damping. Any other damping effects will need to be evaluated
C separately)
C
C Note: The coefficients required for time stepping must be initialized
C        by calling DDTINI before the first call to DDTPGQ.

C NOTE:  The -Laplacian of the surface geopotential must be added
C        in a user supplied forcing routine to the divergence tendency.

C Input:  COMPLEX PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX) -- spectral
C         representation of prognostic quantity at (current time) - DT1
C         COMPLEX PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX) -- spectral
C         representation of prognostic quantity at current model time
C         COMPLEX PSSP0(0:M1MAX,0:N1MAX) -- spectral
C         representation of ln(surface pressure) at previous model time
C         COMPLEX PSSSP1(0:M1MAX,0:N1MAX) -- spectral
C         representation of ln(surface pressure) at current model time
C Output: COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX) -- spectral
C                   "adiabatic" time-tendency of prognostic quantities
C         COMPLEX DPSSP(0:M1MAX,0:N1MAX) -- spectral
C                   "adiabatic" time-tendency of ln(surface pressure)
C
      IMPLICIT NONE

      INCLUDE 'mcons.inc'
      INCLUDE 'spcons.inc'
      INCLUDE 'mgrid.inc'
      INCLUDE 'spgrid.inc'
      INCLUDE 'tmarch.inc'

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
C     (Mnemonic for L-(1/2), L+0, L+(1/2), L+1)
      INTEGER I, J, K, L, L0, M, N, IPGQ
      REAL TRC(K1MAX, K2MAX, 0:1), VTRCFX(K1MAX, K2MAX, 0:1)
      REAL TRCCAP(K1MAX, K2MAX)
      REAL DAYINV, SURFD
      
      COMPLEX SURF1(0:100,0:100)
      COMMON /SURFACE/ SURF1
      
C ----------------------------------------------------------------------
C     Initialize all spectral time-tendencies to zero
      DO 10 L=1,L1
          CALL ZEROSP( DPSSP(0,0), N1)
        DO 10 IPGQ=1,3+NTRACE
          CALL ZEROSP( DPGQSP(0,0,IPGQ,L), N1)
   10     CONTINUE

      DO 15 K=1,K2
      DO 15 J=1,K1
         DBAR(J,K) = 0.
 15      GBAR(J,K) = 0.
C ----------------------------------------------------------------------
C **Note Values of prognostic quantities at current model time
C        (PGQSP1) are used in the following calculations. i.e. The
C        treatment is leap-frog in time.  These are the non-linear 
C        flux terms.
C ----------------------------------------------------------------------
C     Integrate continuity equation to find surface pressure tendency
C     (q=ln(p_surf))

      DO 50 L=1,L1

C       Compute streamfunction PSI and velocity potential CHI
           CALL IDELSQ( FSP1, PGQSP1(0,0,JVOR,L), N1)
           CALL IDELSQ( FSP2, PGQSP1(0,0,JDIV,L), N1)

C       Then compute horizontal velocity from PSI and CHI
           CALL HVELOC( U(1,1,JX,L), U(1,1,JY,L), FSP1, FSP2 )

C     Similarly, find Delq
           CALL GRAD (DELQ(1,1,JX,L), DELQ(1,1,JY,L), PSSP1(0,0))

C     Now find V.Delq and Div and accumulate in -GBAR and DBAR
           DO 20 K=1,K2
           DO 20 J=1,K1
              VDQ(J,K,L)=U(J,K,JX,L)*DELQ(J,K,JX,L)+
     1                   U(J,K,JY,L)*DELQ(J,K,JY,L)
 20           GBAR(J,K) = GBAR(J,K) - VDQ(J,K,L) * DP(L)


           CALL PHYSIC( DIV(1,1,L), PGQSP1(0,0,JDIV,L) )

           DO 25 K=1,K2
           DO 25 J=1,K1
 25           DBAR(J,K) = DBAR(J,K) + DIV(J,K,L) * DP(L)
 50   CONTINUE

      CALL SPECTR (DPSSP, GBAR)

C ----------------------------------------------------------------------
C     Start main altitude loop
      LMHALF= 0
      LPZERO= 0

      DO 300 L=0,L1
        LMHALF= MOD(LMHALF+1, 2)
        LPHALF= MOD(LMHALF+1, 2)
        LPZERO= MOD(LPZERO+1, 2)
        LPONE=  MOD(LPZERO+1, 2)

        IF (L.LT.L1) THEN
C ----------------------------------------------------------------------
C       Compute "perturbation" temp at level L+1
           CALL SPCOPY( FSP1, PGQSP1(0,0,JPOT,L+1), N1 )
           FSP1(0,0)= FSP1(0,0) - TREFLV(L+1)
           CALL PHYSIC( TP(1,1,LPONE), FSP1 )

C       Compute tracer concentration in physical space (at level L+1)
C           IF (NTRACE .GT. 0) THEN
C              CALL PHYSIC( TRC(1,1,LPONE), PGQSP1(0,0,JTR1,L+1) )
C           ENDIF
C ----------------------------------------------------------------------
        ENDIF
C
C       Compute omega and vertical fluxes at level L+1/2
        IF ( (L.EQ.0) .OR. (L.EQ.L1) ) THEN
C ----------------------------------------------------------------------
C         Omega is zero at the boundary half-levels 1/2 and L1+1/2
C         => All vertical fluxes are also zero
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
C ----------------------------------------------------------------------
        ELSE
C ----------------------------------------------------------------------
C         L+1/2 is not at the boundaries. Compute omega and vertical
C         fluxes at level L+1/2.
          DO 130 I=JX,JY
            DO 130 K=1,K2
              DO 130 J=1,K1
                UCAP(J,K,I)= 0.5 * (U(J,K,I,L+1) - U(J,K,I,L))
  130           CONTINUE

          DO 140 K=1,K2
            DO 140 J=1,K1
              TPCAP(J,K) = 0.5 * ( TP(J,K,LPONE)/PKLV(L+1) 
     1                           - TP(J,K,LPZERO)/PKLV(L))

	      IF (NTRACE .GT. 0)
	1       TRCCAP(J,K)=   0.5 * ( TRC(J,K,LPONE) - TRC(J,K,LPZERO))
  140         CONTINUE

C         Compute "pressure velocity" omega at L+1/2
          DO 150 K=1,K2
            DO 150 J=1,K1
              W(J,K,LPHALF)= W(J,K,LMHALF) - DP(L) * 
     1              (DIV(J,K,L) + VDQ(J,K,L) - DBAR(J,K) + GBAR(J,K))
  150         CONTINUE

C         Compute vertical fluxes at level L+1/2
          DO 160 I=JX,JY
            DO 160 K=1,K2
              DO 160 J=1,K1
                VUFLX(J,K,I,LPHALF) = -W(J,K,LPHALF) * UCAP(J,K,I)
  160           CONTINUE

          DO 170 K=1,K2
            DO 170 J=1,K1
              VTPFLX(J,K,LPHALF) = -W(J,K,LPHALF) * TPCAP(J,K)

	      IF (NTRACE .GT. 0)
	1       VTRCFX(J,K,LPHALF) = -W(J,K,LPHALF) * TRCCAP(J,K)
  170         CONTINUE
C ----------------------------------------------------------------------
        ENDIF
C
        IF (L.GT.0) THEN
C ----------------------------------------------------------------------
C     Compute vorticity/divergence tendencies at level 
C     'Linear' contribution from geopotential added to D tendency later
C ----------------------------------------------------------------------
C
C         Absolute vorticity at level L
          CALL SPCOPY( FSP1, PGQSP1(0,0,JVOR,L), N1 )

C         Add planetary vorticity to relative vorticity
          FSP1(0,1)= FSP1(0,1) + CMPLX(FSP01, 0.)
          CALL PHYSIC( ABSVOR, FSP1 )

C         Horizontal "flux convergence" of momentum at level L
          DO 210 K=1,K2
            DO 210 J=1,K1
              TEMP =  RGAS * TP(J,K,LPZERO)
              FLXPH(J,K,JX)= + ABSVOR(J,K) * U(J,K,JY,L)
     1              - TEMP * DELQ(J,K,JX,L)
              FLXPH(J,K,JY)= - ABSVOR(J,K) * U(J,K,JX,L)
     1              - TEMP * DELQ(J,K,JY,L)
  210         CONTINUE

C         Add contribution from vertical flux divergence
          TEMP= 1.0 / DP(L)
          DO 220 I=JX,JY
            DO 220 K=1,K2
              DO 220 J=1,K1
                FLXPH(J,K,I)= FLXPH(J,K,I) +
     1               TEMP * (VUFLX(J,K,I,LPHALF)+VUFLX(J,K,I,LMHALF))
  220           CONTINUE

C         Vorticity/divergence tendencies from "flux conv."
          CALL CURLZ( DPGQSP(0,0,JVOR,L), FLXPH(1,1,JX), FLXPH(1,1,JY))
          CALL DIVERG(DPGQSP(0,0,JDIV,L), FLXPH(1,1,JX), FLXPH(1,1,JY))
C
C         Compute kinetic energy at level L
          DO 230 K=1,K2
            DO 230 J=1,K1
              FPH(J,K)= 0.5*(U(J,K,JX,L)**2 + U(J,K,JY,L)**2)
  230         CONTINUE

C         Convert to spectral space
          CALL ZEROSP( FSP1, N1 )
          CALL SPECTR( FSP1, FPH )

C         Compute the Laplacian, and add contribution to D tendency
          CALL DELSQ( FSP1, FSP1, N1 )
          DO 240 N=0,N1
            DO 240 M=0,MIN0(N,M1)
              DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) - FSP1(M,N)
  240         CONTINUE

C ----------------------------------------------------------------------
C         Compute temp tendency at level L
C ----------------------------------------------------------------------
C         Compute vertical "perturbation" temp flux convergence
          TEMP= PKLV(L) / DP(L)
          DO 250 K=1,K2
            DO 250 J=1,K1
              FPH(J,K)= TEMP * (VTPFLX(J,K,LPHALF)+VTPFLX(J,K,LMHALF))
  250         CONTINUE

C         Add horizontal "flux" of "perturbation" temp
          CALL GRAD (FLXPH(1,1,JX), FLXPH(1,1,JY), PGQSP1(0,0,JPOT,L))
           
          DO 260 I=JX,JY
            DO 260 K=1,K2
              DO 260 J=1,K1
                FPH(J,K)= FPH(J,K)- FLXPH(J,K,I) * U(J,K,I,L)
  260           CONTINUE

C     Surface pressure tendency terms. Recall NNT2DT = -YZ OF H-W
          DO 265 K=1,K2
             DO 265 J=1,K1
                TEMP = TREFLV(L)+ TP(J,K,LPZERO)
                TEMP = TEMP * (VDQ(J,K,L)+GBAR(J,K))
                TEMP = TEMP - TP(J,K,LPZERO) * DBAR(J,K)
                FPH(J,K) = FPH(J,K) + KAPPA*TEMP

                DO 265 L0=1,L1
 265               FPH(J,K) = FPH(J,K) + NNT2DT(L,L0) * VDQ(J,K,L0)

C         Convert to spectral space and add to temp tendency
          CALL SPECTR( DPGQSP(0,0,JPOT,L), FPH )
C
C     All adiabatic tendency terms, excluding the geopotential in the
C     D tendency equation and the linearized heating term in the temp
C     tendency equation, have now been computed for level L.
C
C ----------------------------------------------------------------------
C         Compute tracer time tendency at level L
C ----------------------------------------------------------------------
C
C	  IF (NTRACE .GT. 0) THEN

C         Compute vertical tracer flux convergence
C            TEMP= 1.0 / DP(L)
C            DO 270 K=1,K2
C              DO 270 J=1,K1
C                FPH(J,K)= TEMP * (VTRCFX(J,K,LPHALF)+VTRCFX(J,K,LMHALF))
C  270           CONTINUE
C
C         Add horizontal "flux" of tracer
C          CALL GRAD (FLXPH(1,1,JX), FLXPH(1,1,JY), PGQSP1(0,0,JTR1,L))
C
C            DO 280 I=JX,JY
C              DO 280 K=1,K2
C                DO 280 J=1,K1
C                  FPH(J,K)= FPH(J,K) - FLXPH(J,K,I) * U(J,K,I,L)
C  280             CONTINUE
C
C         Convert to spectral space and add to temp tendency
C            CALL SPECTR( DPGQSP(0,0,JTR1,L), FPH )
C	  ENDIF
C ----------------------------------------------------------------------
        ENDIF
C
  300   CONTINUE

C
C     Now we compute the remaining terms of the tendency equation,
C     including the damping terms. (These are the linear terms)
C     **Note: Values of prognostic quantities at model time TIM - DT1
C             (PGQSP0) are used in the following calculations. i.e. The
C             treatment is explicit in time.
C
      DO 500 L=1,L1
        RT = RGAS * TREFLV(L)

        DO 500 N=0,N1
C ----------------------------------------------------------------------
C         Add geopotential and q contribution to D tendency
          TEMP= A0QINV*N*(N+1)

          DO 410 L0=1,L1
            DO 410 M=0,MIN0(N,M1)
              DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1              + TEMP * (T2GPCP(L,L0) * PGQSP0(M,N,JPOT,L0))
  410         CONTINUE

          DO 415 M=0,MIN0(N,M1)
  415         DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L) +
     $            RT * TEMP * PSSP0(M,N)

C
C         Add contribution from reference static stability to temp
C         tendency (Recall D2TT = - Q of H-W)
          DO 420 L0=1,L1
            DO 420 M=0,MIN0(N,M1)
              DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)
     1             + D2TT(L,L0) * PGQSP0(M,N,JDIV,L0)
  420         CONTINUE
C
C         Add Del**8 damping terms to VOR/DIV/POT tendencies
          DO 430 M=0,MIN0(N,M1)
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)
     1            - UD8FAC(N) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1            - UD8FAC(N) * PGQSP0(M,N,JDIV,L)
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)
     1            - TD8FAC(N) * PGQSP0(M,N,JPOT,L)
  430       CONTINUE

C	  IF (NTRACE .GT. 0) THEN
C
C         Add Del**8 damping term to tracer tendency
C            DO 440 M=0,MIN0(N,M1)
C              DPGQSP(M,N,JTR1,L)= DPGQSP(M,N,JTR1,L)
C     1            - TD8FAC(N) * PGQSP0(M,N,JTR1,L)
C  440         CONTINUE
C	  ENDIF

  500     CONTINUE

C     Add divergence contribution
C     and Del**8 damping to surface pressure tendency
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
      
C      DAYINV= 1.0 / 86400.
C      SURFD= 0.0 * DAYINV
C      
C     nudge surface temperature to assimilated state
C      DO N=0,N1
C        DO M=0,M1
C          DPSSP(M,N)= DPSSP(M,N) 
C     1      + SURFD * (SURF1(M,N) - PSSP0(M,N))
C        ENDDO
C      ENDDO
           
      RETURN
      END
 
C========================================================================
      SUBROUTINE GRAD(UPH,VPH,CHISP)
C========================================================================
C     GRAD  calculates the gradient dCHI/dx, dCHI/dy (UPH, VPH)  in physical
C     space given the spectral space representation of a function CHISP.
C     Input:  COMPLEX CHISP(0:M1MAX,0:N1MAX) -- spectral function
C     Output: REAL UPH(K1MAX,K2MAX) -- zonal component DCHI_DX
C             REAL VPH(K1MAX,K2MAX) -- meridional component DCHI_DY
C
C     This is the 'velocity potential' part of HVELOC, stripped down for
C     efficiency.

      IMPLICIT NONE
      INCLUDE 'spcons.inc'
      INCLUDE 'spgrid.inc'
      INCLUDE 'sppoly.inc'

      REAL UPH(K1MAX,K2MAX), VPH(K1MAX,K2MAX)
      COMPLEX CHISP(0:M1MAX,0:N1MAX)
      COMPLEX ULA(0:M1MAX,K2MAX), VLA(0:M1MAX,K2MAX)
      COMPLEX FTMK1(0:M1MAX)
      REAL RTEMP1, RTEMP2
      INTEGER K, M

C     Calculate velocity components (ULA,VLA) in latitude space
      DO 10 K=1,K2
        RTEMP1= A0INV * COSINV(K)
        RTEMP2= A0INV * COSPHI(K)

         CALL XSP2LA(FTMK1,CHISP,P2(0,0,K),H2(0,0,K),N1,K)

         DO 11 M=0,M1
 11         ULA(M,K)= RTEMP1*FTMK1(M)

         CALL YSP2LA(FTMK1,CHISP,P2(0,0,K),H2(0,0,K),N1,K)

         DO 10 M=0,M1
 10         VLA(M,K)= RTEMP2*FTMK1(M)

C     Convert to physical space
      CALL L2P(UPH,ULA)
      CALL L2P(VPH,VLA)
      
      RETURN
      END
 
C========================================================================
      SUBROUTINE DDTFRC( DPGQSP, PGQSP0, TIM, RAYFAC, EDDIF, ZS)
C========================================================================
C     DDTFRC adds forcing/damping terms to the spectral time tendency DPGQSP
C     (A call to DDTFRC, if any, should immediately follow the call to
C     DDTPGQ that computes DPGQSP)
C
C     Input:  COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
C             "adiabatic" spectral time-tendency of prognostic quantities
C             COMPLEX PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
C             spectral representation of prognostic quantity at all
C             levels at (current model time) - DT1
C             COMPLEX PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
C             spectral representation of prognostic quantity at all
C             levels at current model time
C     Output: COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
C             spectral time-tendency, with effects of forcing/damping added
C
C     **Note: Values of prognostic quantities at current model time - DT1
C             (PGQSP0) are used in the following calculations. i.e., the
C             treatment is explicit in time. Numerical instabilities otherwise.

      IMPLICIT NONE

      INCLUDE 'mcons.inc'
      INCLUDE 'spcons.inc'
      INCLUDE 'mgrid.inc'
      INCLUDE 'spgrid.inc'
      INCLUDE 'tmarch.inc'

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
      
      REAL EDDIF, ZS(L1MAX)
      COMPLEX DAMPU(L1MAX), DAMPT(L1MAX)
      
      REAL WCLAMP1, WCLAMP2, WCLAMP3
      INTEGER TROPLVL, STRATLVL
      COMMON /WAVECLAMP/ WCLAMP1, WCLAMP2, WCLAMP3, TROPLVL, STRATLVL   
      
      REAL ANDREWS, Z1, Z2, ETA0, A1, A2, Z
      REAL KZZTIM, B1, B2, B3, B4, B5, B6, B7, B8, B9, KZZPHASE
      REAL RHO, MOL, CONVFAC, GAUSS
      
      REAL TIM, RAYFAC ! dummies to make driver compatible with other versions of DDTFRC
      
C======================================================================================================================
C     Nudging to assimilated files, Newtonian cooling, and Rayleigh friction damping terms. 
C======================================================================================================================
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
                
        ETA0= (B1 + B2*SIN(KZZPHASE) + B3*COS(KZZPHASE) + B4*SIN(2*KZZPHASE)
     1         + B5*COS(2*KZZPHASE) + B6*SIN(3*KZZPHASE) + B7*COS(3*KZZPHASE)
     1         + B8*SIN(4*KZZPHASE) + B9*COS(4*KZZPHASE))
             
        ETA0=ETA0*5000**2 
        
C        GAUSS = 100*EXP(-(KZZPHASE - (45. / 365. * 2 * 3.14159265359))**2 / (0.25**2))
C        ETA0 = ETA0 - GAUSS
        
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
            
C     MLT nudging 
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
                          
C           Nudge temperature to assimilated state with rate ALPHAT
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)
     1            + TEMPAT * TMNLV(M,N,L) * DAMPFAC
     1              - TEMPAT * PGQSP0(M,N,JPOT,L)
     
C           Nudge wind to assimilated state with rate ALPHAW
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)
     1            + TEMPAW * VOR0(M,N,L) * DAMPFAC
     1              - TEMPAW * PGQSP0(M,N,JVOR,L)
     
            IF (M .EQ. 0) THEN
C             damp mean state divergence field to zero with rate mdamp 
              DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1          - MDAMP(L) * PGQSP0(M,N,JDIV,L)
            ENDIF
            
          ENDDO

          DO M=1,MIN0(N,M1)
C           Rayleigh friction
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)
     1            - URLXLV(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1            - URLXLV(L) * PGQSP0(M,N,JDIV,L)
     
C           Newtonian cooling
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)
     1            - TRLXLV(L) * PGQSP0(M,N,JPOT,L)

C           Sponge layer damp vorticity & divergence
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)
     1            - WDAMP(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1            - WDAMP(L) * PGQSP0(M,N,JDIV,L)
          ENDDO
        ENDDO
      ENDDO

C     stratosphere nudging 
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
                          
C           Nudge temperature to assimilated state with rate ALPHAT
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)
     1            + TEMPAT * TMNLV(M,N,L) * DAMPFAC
     1              - TEMPAT * PGQSP0(M,N,JPOT,L)
     
C           Nudge wind to assimilated state with rate ALPHAW
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)
     1            + TEMPAW * VOR0(M,N,L) * DAMPFAC
     1              - TEMPAW * PGQSP0(M,N,JVOR,L)
     
            IF (M .EQ. 0) THEN
C             damp mean state divergence field to zero with rate mdamp 
              DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1          - MDAMP(L) * PGQSP0(M,N,JDIV,L)
            ENDIF
            
          ENDDO

          DO M=1,MIN0(N,M1)
C           Rayleigh friction
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)
     1            - URLXLV(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1            - URLXLV(L) * PGQSP0(M,N,JDIV,L)
     
C           Newtonian cooling
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)
     1            - TRLXLV(L) * PGQSP0(M,N,JPOT,L)

C           Sponge layer damp vorticity & divergence
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)
     1            - WDAMP(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1            - WDAMP(L) * PGQSP0(M,N,JDIV,L)
          ENDDO
        ENDDO
      ENDDO
      
C     troposphere nudging 
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
            
C           Nudge temperature to assimilated state with rate ALPHAT
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)
     1            + TEMPAT * TMNLV(M,N,L) * DAMPFAC
     1              - TEMPAT * PGQSP0(M,N,JPOT,L)
     
C           Nudge wind to assimilated state with rate ALPHAW
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)
     1            + TEMPAW * VOR0(M,N,L) * DAMPFAC
     1              - TEMPAW * PGQSP0(M,N,JVOR,L)
     
            IF (M .EQ. 0) THEN
C             damp mean state divergence field to zero with rate mdamp 
              DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1          - MDAMP(L) * PGQSP0(M,N,JDIV,L)
            ENDIF
            
          ENDDO

          DO M=1,MIN0(N,M1)
C           Rayleigh friction
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)
     1            - URLXLV(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1            - URLXLV(L) * PGQSP0(M,N,JDIV,L)
     
C           Newtonian cooling
            DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)
     1            - TRLXLV(L) * PGQSP0(M,N,JPOT,L)

C           Sponge layer damp vorticity & divergence
            DPGQSP(M,N,JVOR,L)= DPGQSP(M,N,JVOR,L)
     1            - WDAMP(L) * PGQSP0(M,N,JVOR,L)
            DPGQSP(M,N,JDIV,L)= DPGQSP(M,N,JDIV,L)
     1            - WDAMP(L) * PGQSP0(M,N,JDIV,L)
          ENDDO
        ENDDO
      ENDDO

C ---------------------------------------------------------------
C     Compute vertical viscous damping (for vorticity/temperature)

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

C       No stress at uppermost level
        CALL ZEROSP( STRESS(0,0,LPHALF), N1 )

        DO 350 L=1,L1
          LMHALF= MOD(LMHALF+1, 2)
          LPHALF= MOD(LMHALF+1, 2)

          IF (L.LT.L1) THEN
C            Compute stress at level L+1/2
C            SUBTRACT STRESS OF INITIAL STATE
             IF (IPGQ .EQ. JPOT) THEN
               TEMP = ABS(TDIFF(L))  / (HDP(L+1) + HDP(L)) ! ABS(TDIFF)
               DO 309 N=0,N1
                 DO 309 M=0,MIN0(N,M1)
c                    STRESS(M,N,LPHALF)= TEMP *
c     $                   (PGQSP0(M,N,IPGQ,L+1)/PKLV(L+1) -
c     $                   PGQSP0(M,N,IPGQ,L)/PKLV(L)) 
                   STRESS(M,N,LPHALF)= TEMP *
     $                   (PGQSP0(M,N,IPGQ,L+1) -
     $                    PGQSP0(M,N,IPGQ,L))
 309           CONTINUE

		  DO N=0,N1
		    STRESS(0,N,LPHALF) = STRESS(0,N,LPHALF) - TEMP *
     1        (REAL(TMNLV(0,N,L+1)) - REAL(TMNLV(0,N,L)))
		  ENDDO
             ELSE
                TEMP= ABS(VVISC(L)) / (HDP(L+1) + HDP(L))
                DO 310 N=0,N1
                  DO 310 M=0,MIN0(N,M1)
                    STRESS(M,N,LPHALF)= TEMP *
     1                   (PGQSP0(M,N,IPGQ,L+1) - PGQSP0(M,N,IPGQ,L))
 310              CONTINUE

		  IF (IPGQ .EQ. JVOR) THEN
		    DO N=0,N1
		      STRESS(0,N,LPHALF) = STRESS(0,N,LPHALF) - TEMP *
     1          (REAL(VOR0(0,N,L+1)) - REAL(VOR0(0,N,L)))
		    ENDDO
		  ENDIF
             ENDIF

          ELSE
C           Compute stress at surface (level L1+1/2)
            DO 320 N=0,N1
              DO 320 M=0,MIN0(N,M1)
                 IF(IPGQ .EQ. JPOT) THEN
                    CALL ZEROSP( STRESS(0,0,LPHALF), N1 )
                 ELSE
                    STRESS(M,N,LPHALF) = -USDRAG * PGQSP0(M,N,IPGQ,L)
                 ENDIF
  320           CONTINUE
          ENDIF

C         Add contribution from viscous stress divergence
          DO 330 N=0,N1
            DO 330 M=0,MIN0(N,M1)
              IF(IPGQ .EQ. JPOT) THEN
c                DPGQSP(M,N,IPGQ,L)= DPGQSP(M,N,IPGQ,L) + PKLV(L) * 
c     1               (STRESS(M,N,LPHALF) - STRESS(M,N,LMHALF)) / DP(L)
                DPGQSP(M,N,IPGQ,L)= DPGQSP(M,N,IPGQ,L) +  
     1               (STRESS(M,N,LPHALF) - STRESS(M,N,LMHALF)) / DP(L)
              ELSE
                DPGQSP(M,N,IPGQ,L)= DPGQSP(M,N,IPGQ,L) +
     1               (STRESS(M,N,LPHALF) - STRESS(M,N,LMHALF)) / DP(L)
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

C========================================================================
      SUBROUTINE DDTIMP( DPGQSP, DPSSP )
C========================================================================
C     DDTIMP makes implict corrections to the spectral time tendency DPGQSP.
C     (DDTIMP should be called after the calls to DDTPGQ and DDTFRC to
C     compute the "adiabatic" time-tendency and to include the effects of
C     forcing/damping)
C     Input:  COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
C             explicit spectral time-tendency of prognostic quantities
C     Output: COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX) --
C             implicit spectral time-tendency of prognostic quantities
C
      IMPLICIT NONE

      INCLUDE 'mcons.inc'
      INCLUDE 'spcons.inc'
      INCLUDE 'mgrid.inc'
      INCLUDE 'spgrid.inc'
      INCLUDE 'tmarch.inc'

      COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      COMPLEX DPSSP(0:M1MAX,0:N1MAX)
      INTEGER L0, L, M, N
      REAL RTEMP, RT
      COMPLEX TEM(0:M1MAX,L1MAX)

      DO 200 N=0,N1
C       Making implicit corrections to D tendency
        DO 30 L=1,L1

C         Contribution from explicit D tendency, with Del**8 correction
          DO 10 M=0,MIN0(N,M1)
   10       TEM(M,L)= UD8COR(N) * DPGQSP(M,N,JDIV,L)

C         Add contribution from explicit temp tendency
          RTEMP= UD8COR(N) * (A0QINV*N*(N+1)) * TD8COR(N) * DTFAC
          RT= RTEMP * RGAS * TREFLV(L)

          DO 20 L0=1,L1
            DO 20 M=0,MIN0(N,M1)
   20         TEM(M,L)=TEM(M,L) + RTEMP*T2GPCP(L,L0)*DPGQSP(M,N,JPOT,L0)

C         Add contribution from explicit q tendency
          DO M=0,MIN0(N,M1)
             TEM(M,L) = TEM(M,L) + RT * DPSSP(M,N)
          ENDDO

   30     CONTINUE
C
C       Multiply TEM by implicit correction matrix to find implicit D
C       tendency
        DO 60 L=1,L1
          DO 40 M=0,MIN0(N,M1)
   40       DPGQSP(M,N,JDIV,L)= (0., 0.)

          DO 50 L0=1,L1
            DO 50 M=0,MIN0(N,M1)
   50         DPGQSP(M,N,JDIV,L)=  DPGQSP(M,N,JDIV,L)
     1                           + IMPCOR(L,L0,N)*TEM(M,L0)
   60     CONTINUE
C
        DO 100 L=1,L1
C
C         Make implicit correction to temp tendency due to
C         implicit D tendency (recall D2TT = -M(D->T))
          DO 70 L0=1,L1
            DO 70 M=0,MIN0(N,M1)
              DPGQSP(M,N,JPOT,L)= DPGQSP(M,N,JPOT,L)
     1              + DTFAC * D2TT(L,L0) * DPGQSP(M,N,JDIV,L0)
   70         CONTINUE
C
C         Make implicit correction due to Del**8 damping of temp
          DO 80 M=0,MIN0(N,M1)
   80       DPGQSP(M,N,JPOT,L)= TD8COR(N) * DPGQSP(M,N,JPOT,L)

C	  IF (NTRACE .GT. 0) THEN
C
C         Make implicit correction due to Del**8 damping of tracer
C            DO 90 M=0,MIN0(N,M1)
C   90         DPGQSP(M,N,JTR1,L)= TD8COR(N) * DPGQSP(M,N,JTR1,L)
C	  ENDIF

  100     CONTINUE

C         Make implicit correction to q tendency due to 
C         implicit D tendency
          DO L0=1,L1
             DO M=0,MIN0(N,M1)
                DPSSP(M,N)= DPSSP(M,N)
     1              - DTFAC * DP(L0) * DPGQSP(M,N,JDIV,L0)
             ENDDO
          ENDDO
C
C         Make implicit correction due to Del**8 damping of q
          DO 110 M=0,MIN0(N,M1)
  110       DPSSP(M,N)= TD8COR(N) * DPSSP(M,N)

  200   CONTINUE
  
C     nudge surface temperature to assimilation state
C      DO N=0,N1
C        DO M=0,M1
C          DPSSP(M,N)= DPSSP(M,N)
C     1      + SURFD * (SURF1(M,N) - PSSP1(M,N))
C     1      - SURFD * DPSSP(M,N) * 0
c        ENDDO
c      ENDDO

      RETURN
      END
 
C========================================================================
      SUBROUTINE ROBFIL( FSP0, FSP1, FSP2, ROBFAC )
C========================================================================
C     Given the spectral states FSP0, FSP1, FSP2 at three successive
C     time-steps, ROBFIL performs a Robert-filter operation on the
C     states, with filter weight ROBFAC.
C     Input:  COMPLEX FSP0(0:M1MAX,0:N1MAX) -- spectral
C             representation of prognostic quantity at model time TIM-DT
C             COMPLEX FSP1(0:M1MAX,0:N1MAX) -- spectral
C             representation of prognostic quantity at model time TIM
C             COMPLEX FSP2(0:M1MAX,0:N1MAX) -- spectral
C             representation of prognostic quantity at model time TIM+DT
C             REAL ROBFAC -- filter weighting factor
C     Output: (Filtered versions of FSP0, FSP1, FSP2)
C
      IMPLICIT NONE

      INCLUDE 'spcons.inc'
      INCLUDE 'spgrid.inc'

      COMPLEX FSP0(0:M1MAX,0:N1MAX)
      COMPLEX FSP1(0:M1MAX,0:N1MAX), FSP2(0:M1MAX,0:N1MAX)
      REAL ROBFAC
      INTEGER M, N
      REAL TEMP

      TEMP= 1.0 - 2.0 * ROBFAC
      DO 10 N=0,N1
        DO 10 M=0,MIN0(N,M1)
          FSP1(M,N)=  TEMP  *  FSP1(M,N)
     1              + ROBFAC * (FSP0(M,N) + FSP2(M,N))
   10     CONTINUE
      RETURN
      END
      
C========================================================================
      SUBROUTINE MODROBFIL( FSP0, FSP1, FSP2, ROBFAC )
C========================================================================
C     Given the spectral states FSP0, FSP1, FSP2 at three successive
C     time-steps, ROBFIL performs a Robert-filter operation on the
C     states, with filter weight ROBFAC.
C     Input:  COMPLEX FSP0(0:M1MAX,0:N1MAX) -- spectral
C             representation of prognostic quantity at model time TIM-DT
C             COMPLEX FSP1(0:M1MAX,0:N1MAX) -- spectral
C             representation of prognostic quantity at model time TIM
C             COMPLEX FSP2(0:M1MAX,0:N1MAX) -- spectral
C             representation of prognostic quantity at model time TIM+DT
C             REAL ROBFAC -- filter weighting factor
C     Output: (Filtered versions of FSP0, FSP1, FSP2)
C
C     this version displaces the state at the future time slice as well 
C     as the current timeslice, based on a 2009 Paul D. Williams paper.
C
C     alpha = 0.53 would be a typical value for a semi-implicit scheme, 
C     but don't go below 0.5; will make it unconditionally unstable. 
C     alpha = 1.0 reduces this procedure to the regular Robert-filter

      IMPLICIT NONE

      INCLUDE 'spcons.inc'
      INCLUDE 'spgrid.inc'

      COMPLEX FSP0(0:M1MAX,0:N1MAX)
      COMPLEX FSP1(0:M1MAX,0:N1MAX), FSP2(0:M1MAX,0:N1MAX)
      REAL ROBFAC, ALPHA
      INTEGER M, N
      
      COMPLEX TEMP(0:M1MAX,0:N1MAX)
      
      ALPHA= 0.53
      
      DO N=0,N1
        DO M=0,MIN0(N,M1)
C         displacement factor
          TEMP(M,N)= (ROBFAC / 2) * (FSP0(M,N) - 2 * FSP1(M,N) 
     1      + FSP2(M,N))
C         apply displacement to current and future time slice
          FSP1(M,N)= FSP1(M,N) + TEMP(M,N) * ALPHA
          FSP2(M,N)= FSP2(M,N) + TEMP(M,N) * (ALPHA - 1)
        ENDDO
      ENDDO
      
      RETURN
      END