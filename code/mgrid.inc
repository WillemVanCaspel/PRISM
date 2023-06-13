C MGRID.INC
C Basic common module of multi-level pressure coordinate model
C  This file is included in the routines of PROGNOS.FOR and should also be
C included in a main driver.
C Must INCLUDE MCONS.INC and SPCONS.INC beforehand
C
C------------------------------------------------
C     Planetary parameters (set by call to PLINI)
C------------------------------------------------
C RADIUS -- planetary radius (in metres)
C OMEGA0 -- angular velocity of planetary rotation (in 1/sec)
C F0     -- coriolis parameter (= 2 * OMEGA0)
C FSP01  -- F0 / SQRT(3.0)
C          (F0 expressed as coefficient of spherical harmonic Y[0,1])
C RGAS   -- specific gas constant (in J/(kg*K))
C CP     -- specific heat at constant pressure (in J/(kg*K))
C KAPPA  -- RGAS / CP (an adiabatic exponent)
C G0     -- gravitational acceleration (in m/s**2)
C G0SQ     -- G0*G0, (used in dry adjustment)
C
      REAL RADIUS, OMEGA0, F0, FSP01, RGAS, CP, KAPPA, G0, G0SQ
      COMMON/PLPARM/  RADIUS, OMEGA0, F0, FSP01, RGAS, CP, KAPPA, G0, G0SQ
C
C------------------------------------------------
C     Vertical discretization details (set by call to VERINI)
C------------------------------------------------
C
C L1       -- number of pressure-levels ( >= 2, <=L1MAX)
C	      (set to NLEV by call to VERINI)
      INTEGER L1
      COMMON/VGRID1/ L1
C
C------------------------------------------------
C PSURF    -- surface pressure (in N/m**2)
C DP(L)    -- pressure-thickness of L-th level
C HDP(L)   -- (1/2) * DP(L)
C HDELP(L)   -- HDP(L+1) + HDP(L)
C PLV(L)   -- pressure in the middle of L-th level
C PHLV(L)  -- pressure values at half-levels (PLV(L)+HDP(L))
C             (PHLV(0)= 0.0, PHLV(L1)= PSURF)
C PKLV(L)  -- (PLV(L)/PSURF)**KAPPA at L-th pressure level
C PKCHLV(L)-- vertical differential of PKLV (PK-cap)
C             {i.e. (1/2)*(PKLV(L+1) - PKLV(L)) }
C TSTDLV(L)-- standard value of THETA at L-th level
C ZSTDLV(L)-- standard height above surface (in metres) of L-th level
C             (computed using the standard value of THETA)
C TREFLV(L)-- reference value of THETA at L-th level
C             (used to linearize THETA for semi-implicit stepping)
      REAL PSURF, DP, HDP, PLV, PHLV, PKLV, PKCHLV
      REAL TSTDLV, ZSTDLV, TREFLV, HDELP
      COMMON/VGRID2/ PSURF, DP(L1MAX), HDP(L1MAX), PLV(L1MAX),
     1               PHLV(0:L1MAX), PKLV(L1MAX), PKCHLV(L1MAX),
     2               TSTDLV(L1MAX), ZSTDLV(L1MAX), TREFLV(L1MAX+1),
     1  		HDELP(L1MAX)
C
C------------------------------------------------
C TRFFAC   -- steepness factor for reference THETA profile used in
C             semi-implicit scheme (as compared to the standard
C             THETA profile)
C QHDAMP   -- Del**8 diffusion coefficient (in units of A0**8/sec)
C UD8FAC(N)-- Del**8 damping factor for momentum
C TD8FAC(N)-- Del**8 damping factor for pot. temperature
C DCAP2D   -- DCAP to D conversion matrix
C D2DCAP   -- D to DCAP conversion matrix
C D2W      -- D to W conversion matrix
C D2TT     -- D to THETA tendency conversion matrix
C T2GPCP   -- THETA to GeoPotential-CAP conversion matrix
C NNT2DT   -- THETA to D-tendency (without Del-2 factor) conversion
C              matrix (= DCAP2D . T2GPCP)
C IMPCOR   -- implicit correction matrix (set in DDTINI)
C RMIN     -- minimum Richardson number for dry adjustment to kick in
      REAL TRFFAC, QHDAMP, UD8FAC, TD8FAC, RMIN
      REAL DCAP2D, D2DCAP, D2W, D2TT, T2GPCP, NNT2DT, IMPCOR
      COMMON/VGRID3/ TRFFAC, QHDAMP, UD8FAC(0:N1MAX), TD8FAC(0:N1MAX),
     1  DCAP2D(L1MAX,L1MAX-1), D2DCAP(L1MAX-1,L1MAX),
     2  D2W(L1MAX,L1MAX), D2TT(L1MAX,L1MAX),
     3  T2GPCP(L1MAX,L1MAX), NNT2DT(L1MAX,L1MAX),
     4  IMPCOR(L1MAX,L1MAX,0:N1MAX), RMIN(L1MAX)
C
      integer t_tracer
      real Z_LBC, NudgeAmp, Z_LVL
      integer nudge_coeff
      real NudgeFac(L1MAX)
