module Declarations_mod
 
    implicit none
    public

! MCONS.INC
!  This file is included in the routines of PROGNOS.FOR and should also be
! included in a main driver.
!
! Parameters (array dimensions):
! L1MAX    -- maximum no. of pressure-levels (>=2)
! NPGQ	   -- total # tracers (>=3) (VOR, DIV, POT, passive)
! NTRACE   -- number of "passive" tracers (<= NPGQ-3)
!             (Default value is 0). Set in main program.
!
! Array pointers:
! Prognostic quantities
! JVOR   -- denotes relative vorticity
! JDIV   -- denotes horizontal divergence
! JPOT   -- denotes potential temperature
! JTR1   -- denotes tracer 1
!
! Vector components
! JX     -- denotes zonal direction
! JY     -- denotes meridional direction

    integer, parameter :: L1MAX = 201

	integer, parameter :: NPGQ = 3

	INTEGER NTRACE
	!COMMON /NTR/ NTRACE

	integer, parameter :: JVOR=1, JDIV=2, JPOT=3, JTR1=4, JX=1, JY=2


! SPCONS.INC
!
!  This file is included in the routines of PROGNOS.FOR and SPLIB.FOR.
!  It should also be included in a main driver.
!
! Parameter definition module (array dimensions)
! N1MAX: This is the maximum order of triangular truncation (i.e.
!        order of Legendre polynomials). 
!
! M1MAX: This is the maximum value of zonal wavenumber resolved by the
!        truncation. 
!
! K1MAX: Number of equally spaced longitudes in the transform grid
!        [should be the lowest integral power of 2 >= (3*M1MAX+1) ].
!
! K2MAX: Number of gaussian latitudes in the transform grid [should
!        be the lowest odd integer >= (3*N1MAX+1)/2. 
! (K2MAX is preferably chosen to be odd, so that the equator can be
! one of the gaussian latitudes)
!
! MLOW:  Cut-off value of zonal wavenumber m that denotes "low order
!        zonal truncation". For truncations with very few zonal
!        wavenumbers (M1 <= MLOW), transforms can be performed more
!        efficiently by reversing the order of certain nested
!        do-loops at the heart of the Legendre transform. But the
!        choice of MLOW would be very much machine dependent. The
!        default value is 4.
!
! [Other typical choices of (N1MAX, M1MAX, K1MAX, K2MAX) would be
!     (10, 10,  32, 17) or (21, 21,  64,  33)
!  or (42, 42, 128, 65) or (85, 85, 256, 129) ]
!

	integer, parameter :: N1MAX=42, M1MAX=42, K1MAX=128, K2MAX=65, MLOW=4, N2MAX=N1MAX+1

! MGRID.INC
! Basic- common module of multi-level pressure coordinate model
!  This file is included in the routines of PROGNOS.FOR and should also be
! included in a main driver.
! Must INCLUDE MCONS.IN! and SPCONS.IN! beforehand
!
!----------------------------------------------
!     Planetary parameters (set by call to PLINI)
!----------------------------------------------
! RADIUS -- planetary radius (in metres)
! OMEGA0 -- angular velocity of planetary rotation (in 1/sec)
! F0     -- coriolis parameter (= 2 * OMEGA0)
! FSP01  -- F0 / SQRT(3.0)
!          (F0 expressed as coefficient of spherical harmoni! Y[0,1])
! RGAS   -- specifi! gas constant (in J/(kg*K))
! CP     -- specifi! heat at constant pressure (in J/(kg*K))
! KAPPA  -- RGAS / CP (an adiabati! exponent)
! G0     -- gravitational acceleration (in m/s**2)
! G0SQ     -- G0*G0, (used in dry adjustment)
!
    REAL RADIUS, OMEGA0, F0, FSP01, RGAS, CP, KAPPA, G0, G0SQ
    !COMMON/PLPARM/  RADIUS, OMEGA0, F0, FSP01, RGAS, CP, KAPPA, G0, G0SQ
!
!----------------------------------------------
!     Vertical discretization details (set by call to VERINI)
!----------------------------------------------
!
! L1       -- number of pressure-levels ( >= 2, <=L1MAX)
!	      (set to NLEV by call to VERINI)
    INTEGER L1
    !COMMON/VGRID1/ L1
!
!----------------------------------------------
! PSURF    -- surface pressure (in N/m**2)
! DP(L)    -- pressure-thickness of L-th level
! HDP(L)   -- (1/2) * DP(L)
! HDELP(L)   -- HDP(L+1) + HDP(L)
! PLV(L)   -- pressure in the middle of L-th level
! PHLV(L)  -- pressure values at half-levels (PLV(L)+HDP(L))
!             (PHLV(0)= 0.0, PHLV(L1)= PSURF)
! PKLV(L)  -- (PLV(L)/PSURF)**KAPPA at L-th pressure level
! PKCHLV(L)-- vertical differential of PKLV (PK-cap)
!             {i.e. (1/2)*(PKLV(L+1) - PKLV(L)) }
! TSTDLV(L)-- standard value of THETA at L-th level
! ZSTDLV(L)-- standard height above surface (in metres) of L-th level
!             (computed using the standard value of THETA)
! TREFLV(L)-- reference value of THETA at L-th level
!             (used to linearize THETA for semi-implicit stepping)
    REAL PSURF

    real, dimension(L1MAX)   :: DP, HDP, PLV, PKLV, PKCHLV, ZSTDLV, TSTDLV, HDELP
    real, dimension(L1MAX+1) :: TREFLV
    real, dimension(0:L1MAX) :: PHLV


! MPGRID
!----------------------------------------------
! TRFFA!   -- steepness factor for reference THETA profile used in
!             semi-implicit scheme (as compared to the standard
!             THETA profile)
! QHDAMP   -- Del**8 diffusion coefficient (in units of A0**8/sec)
! UD8FAC(N)-- Del**8 damping factor for momentum
! TD8FAC(N)-- Del**8 damping factor for pot. temperature
! DCAP2D   -- DCAP to D conversion matrix
! D2DCAP   -- D to DCAP conversion matrix
! D2W      -- D to W conversion matrix
! D2TT     -- D to THETA tendency conversion matrix
! T2GPCP   -- THETA to GeoPotential-CAP conversion matrix
! NNT2DT   -- THETA to D-tendency (without Del-2 factor) conversion
!              matrix (= DCAP2D . T2GPCP)
! IMPCOR   -- implicit correction matrix (set in DDTINI)
! RMIN     -- minimum Richardson number for dry adjustment to kick in
    real TRFFAC, QHDAMP

    real, dimension(0:N1MAX)       :: UD8FAC, TD8FAC
    real, dimension(L1MAX,L1MAX-1) :: DCAP2D, D2DCAP
    real, dimension(L1MAX,L1MAX)   :: D2W, D2TT, T2GPCP, NNT2DT
    real, dimension(L1MAX,L1MAX,0:N1MAX) :: IMPCOR
    real, dimension(L1MAX)         :: RMIN

! SPGRID.INC
!  This file is included in the routines of PROGNOS.FOR and SPLIB.FOR.
!  It should also be included in a main driver.
!  Must include SPCONS beforehand
!  Values set by call to SPLINI
!
!----------------------------------------------
! M1        = maximum value of zonal wavenumber M
! N1        = order of triangular truncation  (for scalars)
! K1        = no. of equally spaced longitudes   ( K1 >= 3*M1 + 1 )
! (For efficient F.F.T'S, K1 is always chosen to be a power of 2 )
! K2        = no. of gaussian latitudes          ( K2 >= (3*N1 + 1)/2 )
! (For "aesthetic" reasons, K2 is always chosen to be an odd number)
! N2        = N1 + 1  (order of truncation for vector components)
!
                   INTEGER M1, N1, K1, K2, N2
                   !COMMON/SPGRD1/ M1, N1, K1, K2, N2
             !
             !----------------------------------------------
             ! DEGLON(J) = longitude values (in degrees, 0 ... 360)
             ! DEGLAT(K) = gaussian latitude values (in degrees, -90 ... +90)
             ! LAMBDA(J) = longitude values (in radians, 0 ... 2*Pi)
             ! PHI(K)    = gaussian latitude values (in radians, -Pi/2 ... +Pi/2)
             ! MU(K)     = gaussian "latitude" values ( SIN(PHI(K)), -1 ... +1 )
             ! COSPHI(K) = cosine of latitude ( COS(PHI(K)), 0 ... 1 )
             ! COSINV(K) = reciprocal of cosine of latitude ( 1 / COS(PHI(K)) )
             ! G(K)      = gaussian quadrature weights (with factor of 1/2 included)
             ! D(M,N)    = coefficients of recurrence relation for H[M,N]
             ! CIM(M) = CMPLX(0.0, M)  (coefficient for evaluating LAMBDA-derivative)
             !
                  

                   real, dimension(K1MAX) :: DEGLON, LAMBDA
                   real, dimension(K2MAX) :: DEGLAT, PHI, MU, COSPHI, COSINV, G 
                   real, dimension(0:M1MAX, 0:N2MAX) :: D
                   complex, dimension(0:M1MAX) :: CIM
 

             !
             !----------------------------------------------
             ! A0        = planetary radius (in metres)
             ! A0Q       = A0**2
             ! A0INV     = 1 / A0
             ! A0QINV    = 1 / A0Q
             !
                   REAL :: A0, A0Q, A0INV, A0QINV
                   !COMMON/SPGRD3/ A0, A0Q, A0INV, A0QINV

! TMARCH.INC
!
! Must include MCONS & SPCONS beforehand. This is included in routines of
! PROGNOS.FOR
!
! Parameters used for time-marching
! Set in DDTINI:
                   real IMPFAC, DT, DTFAC
                   real, dimension(0:N1MAX) :: UD8COR, TD8COR
        
             !
             ! Set in VERINI:
                   REAL USDRAG
         
                   real, dimension(L1MAX) :: VVISC
                   complex, dimension(L1MAX) :: URLXLV, TRLXLV
                   complex, dimension(0:M1MAX,0:N1MAX,L1MAX) :: TMNLV
                                                 
!...heat code has 2 km grid. some RADPARAM values hardwired here

!...STUFF RELATED TO RADIATIVE HEATING CODE SETUP
                                 INTEGER KM
                                 PARAMETER (KM=76)
                             
                                 REAL*8 YKH2O(10)/1.0D-3,1.33D-2,4.22D-2,1.334D-1,4.217D-1,  & ! k-coefficients
                                                  1.334D0,5.623D0,31.62D0,177.8D0,1.0D3/          ! for water vapor
                                 REAL*8 DGKH2O(10)/3.0015D-1,5.019D-2,4.56D-2,3.827D-2,      & ! heating rate in 
                                                   2.968D-2,2.284D-2,2.324D-2,1.232D-2,5.16D-3,2.33D-3/ ! the troposphere
                                 COMMON /KDGH2O/ YKH2O, DGKH2O ! used by htuvmd.f solar heating
                             
                             !...f107 =(70,145,220): 60 < F107 < 270 ==> 10.7 cm solar flux index
                                 DOUBLE PRECISION F107 /70./
                                 DOUBLE PRECISION theta27 /0./! = phase of 27-day periodic variation: sin(theta27)
                                 DOUBLE PRECISION OME0 /.25/	!EFFECTIVE ALBEDO
                                 DOUBLE PRECISION DEC	!SOLAR DECLINATION
                             
                             !...mdlrd=1, no heat; 2 => diurnal avg heat rate; 3 => local heat rate;
                                 INTEGER MDLRD		
                                 INTEGER KLBD /6/		! effective cloud-top level
                                 DOUBLE PRECISION XMU/0.5D0/	! COS(solar zenith angle) = DUMMY VARIABLE
                                 INTEGER LTMU/1/		! 1 => input=DEC, PHI, TIME(hr); 0 => direct input = XMU
                             
                                 COMMON /RADPARAM/ F107, THETA27, OME0, DEC, MDLRD, KLBD ! used in heating codes

! SPPOLY.inc
!
! Common module for Legendre polynomial tables

      real, dimension(0:M1MAX, 0:N2MAX) :: PNORM, PRECC1, PRECC2, PK
      real, dimension(0:M1MAX, K2MAX)   :: PRECI
      real, dimension(0:M1MAX, 0:N1MAX) :: HRECC1, HRECC2
      real, dimension(0:M1MAX, 0:N1MAX) :: HK
      real, dimension(0:2*N2MAX+2)      :: FAC

      integer KCUR
      real, dimension(0:2*M1MAX+1,0:N2MAX) :: PK2
      real, dimension(0:2*M1MAX+1, 0:N1MAX):: HK2


      REAL P2(0:2*M1MAX+1, 0:N2MAX, K2MAX), H2(0:2*M1MAX+1, 0:N1MAX, K2MAX)
      !COMMON/SPPOL3/ P2, H2

! ftable.inc
!
! Declaration module for tables of trigonometric coefficients
    integer MAXDAT, MAXTAB
    parameter (MAXDAT=512)
    parameter (MAXTAB=512)
    real, dimension(4,2,MAXTAB) :: FFTAB
!
! SPFFTB.INC
!f
! Common module for f.f.t. table used by sptransforms
    INTEGER NFFTAB
    ! REAL FFTAB
    PARAMETER (NFFTAB=512)
    ! COMMON/SPFFTB/ FFTAB(4,2,NFFTAB)
    

end module Declarations_mod