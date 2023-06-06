!  
!  ...This version uses recurrence to compute normalized Legendre polynomials
!  ...without the use of factorials
!  
! Sptransforms -- A Spherical-harmonic transforms package
!                 version 3.3 created November 16, 1990
!
! Author: R. Saravanan, D.A.M.T.P., University of Cambridge,
!         Silver Street, Cambridge CB3 9EW, U.K.
!         E-mail: atm.amtp.cam.ac.uk::svn
! (Please send bug-reports/comments to above address)
!
! This software is distributed free for research purposes. It comes with
! absolutely no warranties whatsoever. It may be freely copied, and may also
! be modified, provided the modified software continues to be freely
! available under the same terms. This software may not be used for any
! commercial purpose.
!
! ----------- The Spherical Harmoni! Spectral Transforms --------------
!
! The notation and the method are based on "The Spectral Method" by
! Bennert Machenhauer (Chapter 3 of "Numerical Methods used in
! Atmospheric Models", vol.2, G.A.R.P. publications series no.17,
! W.M.O.(1979).
!
! The routines in this package use triangular truncation of order N1
! i.e. N = 0,...,N1  and M = -N,...,+N.
! the range of index M can be further restricted by specifying
! another parameter M1 that may be less than N1. Then we have
!                        M = -MIN(N,M1),...,+MIN(N,M1)
! e.g. The extreme case of M1=0 would correspond to axially
!      symmetric truncation.
!
! LAMBDA denotes longitude and PHI denotes latitude. Instead of PHI,
! we prefer to use MU = SIN(PHI). In the discussion that follows,
! we will use the term "latitude" to denote MU, rather than PHI.
! LAMBDA = 0 ... 2*PI, PHI = -PI/2 ... PI/2, MU = -1 ... +1
!
! We can expand a general function F(LAMBDA,MU) as
!   F(LAMBDA,MU) = Sigma(N=0 to N1) Sigma(M=-N to N) F[M,N] * Y[M,N]
! where F[M,N] are complex coefficients, Y[M,N] are the spherical
! harmonics,
!       Y[M,N](LAMBDA,MU) = P[M,N](MU) * EXP(I*M*LAMBDA),
! and P[M,N] are the associated legendre polynomials.
! The spherical harmonics satisfy the normalization condition
! (1/(4*PI)) * <<CONJG(Y[M,N])*Y[M2,N2]>> = Kronecker-delta[M,M2][N,N2]
! where << ... >> denotes double integration over MU, LAMBDA.
!
! Note: For components of vectors, The truncated series will be one
!       order higher in N as compared to scalars. i.e. Order N1+1
!
! For a real function F(LAMBDA,MU), we have the relation
!       F[-M,N] = (-1)**M * CONJG(F[M,N])
! So we need to store the coefficients F[M,N] only for M >= 0
!
!                    **** Transform grid ****
!
! To use these routines, we must choose the truncation order N1
! and the transform grid of K1 longitudes and K2 latitudes.
! i.e. We have the following ranges
! N = 0,...,N1        (upto N1+1 for vector components)
! M = 0,...,MIN(N,M1) (because we are interested in real functions)
! M2= 0,...,2*N+1     (This strange index is used only for storing the
! Legendre polynomials and their derivatives. Since they always occur
! in the form complex-number * real-polynomial, one can treat the
! complex multiplication as two real multiplications [for efficiency
! and vectorization]. So the polynomial values are stored in duplicate)
! LAMBDA[J] = 2*PI*(J-1)/K1, J = 1,...,K1  (longitudes)
! MU[K],                     K = 1,...,K2  (gaussian latitudes)
 
! Any real function F(LAMBDA,MU) defined on a sphere can therefore
! be represented in three different 2-D spaces:
! 1. Spectral space: Here the function is represented by its spectral
!    coefficients F[M,N]
!    (here F[-M,N] = (-1)**M * CONJG( F[M,N] ) )
! 2. Latitude space: Here the function is represented by its
!    latitude-dependent fourier coefficient F[M](MU[K])
!    (here F[-M](MU[K]) = CONJG( F[M](MU[K]) ) )
! 3. Physical space: Here the function is represented by its values
!    at the grid points F(LAMBDA[J],MU[K])
!
! The routines in this package basically enable you to move from
! one representation to another
!
! --------------- representation of a real function F ----------------
!
! For the transform operations, it would be convenient to declare the
! various representations of function F in the following manner:
! 1. Spectral space-
!       COMPLEX FSP(0:M1MAX,0:N1MAX)
!       (declare as FSP(0:M1MAX,0:N2MAX) for vector components)
! 2. Latitude space-
!       COMPLEX FLA(0:M1MAX,K2MAX)
!    For a single latitude, i.e., single value of K, declare as
!       COMPLEX FLAK(0:M1MAX)
! 3. Physical space-
!       REAL FPH(K1MAX,K2MAX)
!    For a single latitude, i.e., single value of K, declare as
!       REAL FPHK(K1MAX)
!
!
! Subroutines:
!     SPINI  -- Initializes the many tables required for the transforms.
!     LATINI -- Initializes table of Legendre polynomials for
!               for transforms on a single latitude circle.
!
! The following routines convert a function from spectral to physical
! space, and back, through latitude space. To save storage, they often
! operate on function values for one latitude at a time. Given the
! full spectral representation, one can concentrate on one latitude
! at a time for all operations except one. i.e. The conversion from
! latitude space to spectral space. This obviously requires function
! values at all latitudes. Subroutine LA2SP, which is responsible
! for this conversion, does not actually perform the full conversion,
! but accumulates contributions to the spectral coefficients from each
! latitude. To use LA2SP, one should first set the matrix of spectral
! coefficients to zero, using subroutine ZEROSP, and then successively
! call LA2SP for each latitude.
!
! The following routines operate on one latitude circle at a time:
!
! **NOTE** To use these routines at any given latitude, the table of
! Legendre polynomials at that latitude should have been initialized
! by a call to LATINI.
!
!     S2PK   -- converts from spectral space to physical space (S -> P)
!     P2SK   -- converts from physical space to spectral space (P -> S)
!               (Note: This routine actually accumulates)
!     HVELK  -- compute horizontal velocity from PSI and CHI   (S -> P)
!     DIVK   -- computes the divergence of a vector field      (P -> S)
!               (Note: This routine actually accumulates)
!     CURLK  -- computes the curl of a vector field            (P -> S)
!               (Note: This routine actually accumulates)
!
! The above routines are composites of the following more primitive
! routines, which the user may never need to use directly. They also
! operate on one latitude circle at a time.
!
!     SP2LA  -- converts from spectral space to latitude space (S -> L)
!     XSP2LA -- SP2LA of the LAMBDA-derivative                 (S -> L)
!     YSP2LA -- SP2LA of the MU-derivative                     (S -> L)
!     LA2SP  -- converts from latitude space to spectral space (L -> S)
!               (Note: This routine actually accumulates)
!     XLA2SP -- LA2SP of the LAMBDA-derivative                 (L -> S)
!     YLA2SP -- LA2SP of thE MU-derivative (by parts!)         (L -> S)
!     LA2PH  -- converts from latitude space to physical space (L -> P)
!     PH2LA  -- converts from physical space to latitude space (P -> L)
!
! The following routines are operators in spectral space. i.e. They do
! not transform to any other kind of space.
!
!     ZEROSP -- zeroes spectral space representation           (S -> S)
!     SPCOPY -- copies spectral space representations          (S -> S)
!     DELSQ  -- del-squared (Laplacian) operator               (S -> S)
!     IDELSQ -- inverse of DELSQ                               (S -> S)
!     SPDLM  -- caclulates d/dLAMBDA of a function             (S -> S)
!     HSPDMU -- calculates -(1-MU**2)*d/dMU of a function      (S -> S)
!
! The following set of routines are adapted from the book
! Numerical Recipes, W.H.Press et al, Cambridge(1986). They are used
! by the routines described above, but the user may never need to
! bother about them.
!
!    GAULEG -- This subroutine finds the abcissas and weights for the
!              n-point Gauss-Legendre quadrature.
!    FOUR0 -- Sets up a table of trigonometri! recurrence coefficients
!             for efficient fast fourier transforms on fixed length
!             arrays.
!    FOUR2 -- This subroutine is similar to FOUR1, but works only for
!             fixed array length. Routine FOUR0 sets up the table
!             needed for FOUR2.
!    REALF2 -- This subroutine is similar to REALFT, but works only for
!              fixed array lengths. Routine FOUR0 sets up the table
!              needed for REALF2.
!    FOURV -- Vectorizable version of FOUR2 (computes several
!             transforms simultaneously)
!    REALFV -- Vectorizable version of REALF2 (computes several
!             transforms simultaneously)
!
! The following routines consist of operators and mappings useful in
! solving the primitive equations for a "thin" planetary atmosphere.
! These are likely to be the most useful routines of all. These
! routines can be used only if the order of the spectral truncation is
! "low", because they involve all-latitudes-at-a-time transforms. These
! routines are based on the more primitive routines described above.
! Before using these routines, the spectral transform grid must be
! initialized by calling SPINI. (But calls to LATINI will not be
! necessary to use these routines.)
!
! The available routines are as follows:
!     PHCOPY -- copies physical representations                (P -> P)
!     MULCOS -- multiplies by a power of COSPHI                (P -> P)
!
!     PHYSIC -- convert from spectral to physical space        (S -> P)
!     JACOB1 -- computes the jacobian of one spectral function (S -> P)
!     JACOB2 -- computes the jacobian of two spectral functions(S -> P)
!     HVELOC -- calculate horizontal velocity from PSI and CHI (S -> P)
!
! The following routines all accumulate in spectral space. So take care
! to initialize the spectral coefficients using ZEROSP, if necessary.
!     SPECTR -- convert from physical to spectral space        (P -> S)
!     DIVERG -- computes the divergence of a vector field      (P -> S)
!     CURLZ  -- computes the curl of a vector field            (P -> S)
!
! The following routines offer more flexible transforms, but they may
! not need to be used very often.
!     S2P -- converts from spectral space to physical space
!     P2S -- converts from physical space to spectral space
!            (Note: This routine actually accumulates)
!     L2P -- converts from latitude space to physical space
!     P2L -- converts from physical space to latitude space
!
! The following are some subroutines useful in diagnostic
! calculations. Before using these routines, the spectral transform
! grid must be initialized by calling SPINI.
!
!     ZAVGE  -- calculate zonal average of a function in physical space
!     EDDY   -- determine eddy part of a function in physical space
!     GAVGE  -- calculate global average of a function in physical space
!     TOTAL  -- totals up an array of numbers
!     WTOTAL -- computes the weighted total of an array of numbers
!     WAVGE  -- computes the weighted average of an array of numbers
!
!
! Basic spectral transforms module

module splib_mod

    use Declarations_mod

    contains
!========================================================================
      SUBROUTINE SPINI(MMAX, NMAX, RADIUS)
        !========================================================================
        ! SPINI initializes the transform grid for triangular truncation of
        ! order NMAX. It chooses the minimum possible K1, K2 to avoid aliasing
        ! errors. For simplicity, K1 is chosen to be a power of 2. SPINI also
        ! initializes the tables required for efficient fast fourier
        ! transforms.
        ! Input:  INTEGER MMAX, NMAX -- maximum value of spherical harmonic
        !                               indices [M,N]
        !                 (MMAX should always be <= NMAX.
        !                  MMAX=NMAX for triangular truncation.
        !                  MMAX=0 for axially symmetric truncation.)
        !         REAL RADIUS -- planetary radius (in metres)
        !
            IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE  'spgrid.inc' 
            !   INCLUDE 'sppoly.inc'
            !   INCLUDE 'spfftb.inc'
          
              INTEGER MMAX, NMAX
              REAL RADIUS
              INTEGER I,J,K,M,N
              REAL DLAM, PI
              INTEGER M2
        
              PARAMETER (PI=3.14159265)
         
        !========================================================================
        ! Maximum zonal wavenumber
              M1= MIN0( MMAX, NMAX )
        ! No. of longitudinal intervals (always a power of 2)
              K1= 2**(INT( (ALOG(3*M1+1.)/ALOG(2.)) + 0.9999 ) )
        ! Maximum Legendre polynomial order (for scalars)
              N1= NMAX
        ! No. of gaussian latitudes
              K2= INT( (3*N1+1.)/2. + 0.9999 )
        ! Make it an odd number, so that the equator is always included
              IF (MOD(K2,2).EQ.0) K2=K2+1
        ! Maximum truncation order (for components of vectors)
              N2= N1+1
        ! Check if truncation fits into available storage
              IF (     (M1.GT.M1MAX) .OR. (N1.GT.N1MAX)   &
                  .OR. (K1.GT.K1MAX) .OR. (K2.GT.K2MAX) ) &
                  CALL SPERR('SPINI', 'Truncation too large for available storage')
         
         
        ! Initialize table of factorials
        !        FAC(0)= 1.D0
        !        DO 10 I=1,2*N2MAX+2
        !     10   FAC(I)= I * FAC(I-1)
        ! Set-up tables for F.F.T's (latitude space <-> physical space)
              IF (K1.GT.1) CALL FOUR0( FFTAB, NFFTAB, K1/2 )
        ! *** NOTE ***
        ! If the call to FOUR0 results in an error, double the value of
        ! NFFTAB, and try again.
        !
        ! Initialize longitude values
              DLAM= 2.*PI / K1
              DO 20 J=1,K1
                LAMBDA(J)= (J-1)*DLAM
                DEGLON(J)= LAMBDA(J) * 180./PI
           20   CONTINUE
        ! Initialize gaussian "latitudes" and gaussian quadrature weights
              CALL GAULEG( -1.0, 1.0, MU, G, K2)
        ! Initialize latitude angles etc.
              DO 30 K=1,K2
                PHI(K)= SIGN(ASIN( ABS(MU(K)) ), MU(K))
                DEGLAT(K)= PHI(K) * 180./PI
                COSPHI(K)= SQRT(1.0-MU(K)**2)
                COSINV(K)= 1.0/COSPHI(K)
        ! (Include factor of 1/2 in weight G(K) )
                G(K)= 0.5*G(K)
           30   CONTINUE
        ! Initialize coefficient for evaluating LAMBDA-derivative
              DO 40 M=0,M1
           40   CIM(M)= CMPLX( 0.0, FLOAT(M) )
        ! Initialize recurrence coeff. for MU-derivative of Legdr. polynomial
              D(0,0)= 0.0
              DO 50 N=1,N2
                DO 50 M=0,MIN0(N,M1)
                  D(M,N)= SQRT( FLOAT(N**2 - M**2) / (4 * N**2 - 1) )
           50     CONTINUE
        ! Initialize recurrence relations for computing Legendre polynomials
        ! (see Numerical Recipes, p.182)
              DO 60 K=1,K2
                PRECI(0,K)= 1.0
                DO 60 M=1,M1
        !              PRECI(M,K)= -(2*M-1) * COSPHI(K) * PRECI(M-1,K)
                    PRECI(M,K)= -SQRT(1+.5/M) * COSPHI(K) * PRECI(M-1,K)
           60     CONTINUE
        !
              DO 70 N=0,N2
                DO 70 M=0,MIN0(N-1,M1)
        !            PRECC1(M,N)=  FLOAT(2*N-1)/(N-M)
        !            PRECC2(M,N)= -FLOAT(N+M-1)/(N-M)
                  PRECC1(M,N)=  SQRT((4.*N*N-1.)/(N*N-M*M))
                  PRECC2(M,N)= -SQRT((2.*N+1.)*(N-1.+M)*(N-1.-M) &
                    /(2.*N-3.)/(N*N-M*M))
           70     CONTINUE
        ! Compute normalization factor for Legendre polynomial
        ! (as in Machenhauer 1979)
        !        DO 80 N=0,N2
        !             DO 80 M=0,MIN0(N,M1)
                   DO 80 M=0,M1
        !            PNORM(M,N)= ((-1)**M) * DSQRT(DBLE(2*N+1) * FAC(N-M)/FAC(N+M))
                     PNORM(M,1)= ((-1)**M)
           80     CONTINUE
        ! Compute coefficients for finding MU-derivative of Legendre polynomial
              DO 90 N=0,N1
                DO 90 M=0,MIN0(N,M1)
                    HRECC1(M,N)=      N * D(M,N+1)
                    HRECC2(M,N)= -(N+1) * D(M,N)
           90       CONTINUE
        ! Initialize planetary radius and related constants
              A0= RADIUS
              A0Q= A0**2
              A0INV= 1.0 / A0
              A0QINV= 1.0 / A0Q
        
        ! Initialize table of Associated Legendre polynomials and derivatives
              DO 120 K=1,K2
                 CALL POLINI( K )
                 DO 100 N=0,N2
                   DO 100 M2=0,2*MIN0(N,M1)+1
          100        P2(M2,N,K)= PK2(M2,N)
                 DO 110 N=0,N1
                   DO 110 M2=0,2*MIN0(N,M1)+1
          110        H2(M2,N,K)= HK2(M2,N)
          120  CONTINUE
        
        ! Clear polynomial buffer
              KCUR= 0
              RETURN
              END
         
        !========================================================================
              SUBROUTINE LATINI(K)
        !========================================================================
        ! LATINI initializes the tables required for spectral transforms
        ! at a single latitude circle MU=MU(K). It must be called before
        ! any calls to routines S2PK, P2SK, HVELK, DIVK or CURLK at that
        ! particular latitude are made. For efficiency, all such calls at
        ! the same latitude should be grouped together.
        ! Input:  INTEGER K -- Latitude number at which tables are to be
        !                      initialized.
        ! Output: (Tables for latitude MU=MU(K) are initialized.)
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              INTEGER K
              CALL POLINI(K)
        
              RETURN
              END
         
        !========================================================================
              SUBROUTINE POLINI(K)
        !========================================================================
        ! POLINI calculates the Legendre polynomials, and their MU-derivatives,
        ! at latitude MU(K).
        ! (Note: This is an auxiliary routine used by SPINI and LATINI)
        ! Input:  INTEGER K -- Latitude number at which polynomials and their
        !                      derivatives are to be evaluated (MU = MU(K))
        ! Output: (The Associated Legendre polynomials, with normalization as
        !          in Machenhauer(1979), are stored in array PK of common blocks
        !          /SPGRD2/, /SPGRD3/. The corresponding MU-derivatives are
        !          stored in array HK. PK2 and HK2 are also initialized.
        !          KCUR is set to K.)
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
            !   INCLUDE 'sppoly.inc'
        
              INTEGER K
              INTEGER M, N
              REAL*8 MUK
        ! Computing polynomials for latitude K
              KCUR= K
              MUK= MU(K)
        ! Initialize starting values for recurrence
              DO 10 M=0,M1
           10   PK(M,M)= PRECI(M,K)
        ! (N can vary between 0 and N2, but M can only vary between 0 and N1)
              DO 20 M=0,M1
        !     20   PK(M,M+1)= MUK*PRECC1(M,M+1)*PRECI(M,K)
           20   PK(M,M+1)= MUK*SQRT(2.*M+3.)*PRECI(M,K)
        ! Use recurrence relation to compute other values
              DO 30 N=2,N2
                DO 30 M=0,MIN0(N-2,M1)
                  PK(M,N)= MUK*PRECC1(M,N)*PK(M,N-1) &
                             + PRECC2(M,N)*PK(M,N-2)
           30     CONTINUE
        ! Normalize Associated Legendre polynomial (Machenhauer's convention)
              DO 70 N=0,N2
                DO 40 M=0,MIN0(N,M1)
        !     40     PK(M,N)= PNORM(M,N) * PK(M,N)
           40     PK(M,N)= PNORM(M,1) * PK(M,N)
                DO 50 M=0,MIN0(N,M1)
           50     PK2(2*M,  N)= PK(M,N)
                DO 60 M=0,MIN0(N,M1)
           60     PK2(2*M+1,N)= PK(M,N)
           70    CONTINUE
        ! Initialize table of MU-derivative of Associated Legendre polynomials
              HK(0,0)= 0.0
              DO 80 N=1,N1
                DO 80 M=0,MIN0(N,M1)
                  HK(M,N)=  HRECC1(M,N) * PK(M,N+1) &
                          + HRECC2(M,N) * PK(M,N-1)
           80     CONTINUE
        !
              DO 110 N=0,N1
                DO 90 M=0,MIN0(N,M1)
           90     HK2(2*M,  N)= HK(M,N)
                DO 100 M=0,MIN0(N,M1)
          100     HK2(2*M+1,N)= HK(M,N)
          110     CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE S2PK(FPHK,FSP,FSP2LA,K)
        !========================================================================
        ! S2PK converts the spectral space representation FSP to physical
        ! space representation FPH, for a single latitude MU=MU(K). FSP2LA is
        ! the intermediate subroutine to be used to convert from spectral to
        ! latitude space. i.e. FSP2LA = SP2LA or XSP2LA or YSP2LA
        ! Note 1: Truncation order is assumed to be N1 (i.e. as for a scalar)
        ! Note 2: Remember to declare FSP2LA as external to the calling routine
        !
        ! Input:  COMPLEX FSP(0:M1MAX,0:N1) -- spectral representation of F
        !         SUBROUTINE FSP2LA(FLAK,FSP,P0,H0,N0,K) -- spectral to latitude
        !                                        space conversion subroutine
        !         INTEGER K -- latitude number (MU = MU(K))
        ! Output: REAL FPHK(K1MAX) -- physical space representation oF F,
        !                             for MU=MU(K)
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
            !   INCLUDE 'sppoly.inc'
        
              INTEGER K
              REAL FPHK(K1MAX)
              COMPLEX FSP(0:M1MAX,0:N1)
              COMPLEX FTMKC(0:M1MAX)
        
        ! Dummy declaration (sometimes needed with the IMPLICIT NONE option)
        !      INTEGER FSP2LA
               EXTERNAL FSP2LA
        
        ! Convert to latitude space in temporary storage, using FSP2LA
        
              CALL FSP2LA(FTMKC,FSP,P2(0,0,K),H2(0,0,K),N1,K)
        
        ! Convert to physical space
              IF (K1.GT.1) THEN
                ! CALL LA2PH(FPHK,FTMKC)
              ELSE
        ! Axially symmetric case; handle separately
                FPHK(1)= REAL(FTMKC(0))
              ENDIF
              RETURN
              END
         
        !========================================================================
              SUBROUTINE P2SK(FSP,FPHK,FLA2SP,K)
        !========================================================================
        ! P2SK adds to spectral space representation FSP the contribution
        ! from physical space representation FPHK, for a single latitude
        ! MU=MU(K). FLA2SP is the intermediate subroutine to be used to convert
        ! from latitude to spectral space.
        ! i.e. FLA2SP = LA2SP or XLA2SP or YLA2SP
        ! Note 1: Truncation order is assumed to be N1 (i.e. as for a scalar)
        ! Note 2: Remember to declare FSP2LA as external to the calling routine
        !
        ! Input:  COMPLEX FSP(0:M1MAX,0:N1) -- old spectral coefficients
        !         REAL FPHK(K1MAX) -- physical space representation of F,
        !                             for MU=MU(K)
        !         SUBROUTINE FLA2SP(FSP,FLAK,P0,H0,N0,K) -- latitude to spectral
        !                                     space conversion subroutine
        !         INTEGER K -- latitude number (MU = MU(K))
        ! Output: COMPLEX FSP(0:M1MAX,0:N1) -- spectral representation of F
        !         (i.e. Old coefficients plus the contribution from latitude
        !                                                           MU=MU(K) )
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
            !   INCLUDE 'sppoly.inc'
        
              INTEGER K
              COMPLEX FSP(0:M1MAX,0:N1)
              REAL FPHK(K1MAX)
              COMPLEX FTMKC(0:M1MAX)
        
        ! Dummy declaration (sometimes needed with the IMPLICIT NONE option)
        !      INTEGER FLA2SP
              EXTERNAL FLA2SP
        
        ! Convert to latitude space in temporary storage
              IF (K1.GT.1) THEN
                CALL PH2LA(FTMKC,FPHK)
              ELSE
        ! Axially symmetric case; handle separately
                FTMKC(0)= CMPLX(FPHK(1),0.)
              ENDIF
        ! Convert to spectral space using FLA2SP (accumulate)
              CALL FLA2SP(FSP,FTMKC,P2(0,0,K),H2(0,0,K),N1,K)
        
              RETURN
              END

    
!========================================================================
              SUBROUTINE P2SK_XLA2SP(FSP,FPHK,K)
                !========================================================================
                ! P2SK adds to spectral space representation FSP the contribution
                ! from physical space representation FPHK, for a single latitude
                ! MU=MU(K). FLA2SP is the intermediate subroutine to be used to convert
                ! from latitude to spectral space.
                ! i.e. FLA2SP = LA2SP or XLA2SP or YLA2SP
                ! Note 1: Truncation order is assumed to be N1 (i.e. as for a scalar)
                ! Note 2: Remember to declare FSP2LA as external to the calling routine
                !
                ! Input:  COMPLEX FSP(0:M1MAX,0:N1) -- old spectral coefficients
                !         REAL FPHK(K1MAX) -- physical space representation of F,
                !                             for MU=MU(K)
                !         SUBROUTINE FLA2SP(FSP,FLAK,P0,H0,N0,K) -- latitude to spectral
                !                                     space conversion subroutine
                !         INTEGER K -- latitude number (MU = MU(K))
                ! Output: COMPLEX FSP(0:M1MAX,0:N1) -- spectral representation of F
                !         (i.e. Old coefficients plus the contribution from latitude
                !                                                           MU=MU(K) )
                !
                
                      IMPLICIT NONE
                    !   INCLUDE 'spcons.inc'
                    !   INCLUDE 'spgrid.inc'
                    !   INCLUDE 'sppoly.inc'
                
                      INTEGER K
                      COMPLEX FSP(0:M1MAX,0:N1)
                      REAL FPHK(K1MAX)
                      COMPLEX FTMKC(0:M1MAX)
                
                ! Dummy declaration (sometimes needed with the IMPLICIT NONE option)
                !      INTEGER FLA2SP
                      !EXTERNAL FLA2SP
                
                ! Convert to latitude space in temporary storage
                      IF (K1.GT.1) THEN
                        CALL PH2LA(FTMKC,FPHK)
                      ELSE
                ! Axially symmetric case; handle separately
                        FTMKC(0)= CMPLX(FPHK(1),0.)
                      ENDIF
                ! Convert to spectral space using FLA2SP (accumulate)
                      CALL XLA2SP(FSP,FTMKC,P2(0,0,K),H2(0,0,K),N1,K)
                
                      RETURN
                      END
         
!========================================================================
                      SUBROUTINE P2SK_YLA2SP(FSP,FPHK,K)
                        !========================================================================
                        ! P2SK adds to spectral space representation FSP the contribution
                        ! from physical space representation FPHK, for a single latitude
                        ! MU=MU(K). FLA2SP is the intermediate subroutine to be used to convert
                        ! from latitude to spectral space.
                        ! i.e. FLA2SP = LA2SP or XLA2SP or YLA2SP
                        ! Note 1: Truncation order is assumed to be N1 (i.e. as for a scalar)
                        ! Note 2: Remember to declare FSP2LA as external to the calling routine
                        !
                        ! Input:  COMPLEX FSP(0:M1MAX,0:N1) -- old spectral coefficients
                        !         REAL FPHK(K1MAX) -- physical space representation of F,
                        !                             for MU=MU(K)
                        !         SUBROUTINE FLA2SP(FSP,FLAK,P0,H0,N0,K) -- latitude to spectral
                        !                                     space conversion subroutine
                        !         INTEGER K -- latitude number (MU = MU(K))
                        ! Output: COMPLEX FSP(0:M1MAX,0:N1) -- spectral representation of F
                        !         (i.e. Old coefficients plus the contribution from latitude
                        !                                                           MU=MU(K) )
                        !
                        
                              IMPLICIT NONE
                            !   INCLUDE 'spcons.inc'
                            !   INCLUDE 'spgrid.inc'
                            !   INCLUDE 'sppoly.inc'
                        
                              INTEGER K
                              COMPLEX FSP(0:M1MAX,0:N1)
                              REAL FPHK(K1MAX)
                              COMPLEX FTMKC(0:M1MAX)
                        
                        ! Dummy declaration (sometimes needed with the IMPLICIT NONE option)
                        !      INTEGER FLA2SP
                              !EXTERNAL FLA2SP
                        
                        ! Convert to latitude space in temporary storage
                              IF (K1.GT.1) THEN
                                CALL PH2LA(FTMKC,FPHK)
                              ELSE
                        ! Axially symmetric case; handle separately
                                FTMKC(0)= CMPLX(FPHK(1),0.)
                              ENDIF
                        ! Convert to spectral space using FLA2SP (accumulate)
                              CALL YLA2SP(FSP,FTMKC,P2(0,0,K),H2(0,0,K),N1,K)
                        
                              RETURN
                              END
        !========================================================================
              SUBROUTINE HVELK(UPHK,VPHK,PSISP,CHISP,K)
        !========================================================================
        ! HVELK calculates the zonal and meridional components in physical
        ! space (UPHK, VPHK) of the horizontal velocity for a single latitude,
        ! given the spectral space representation of streamfunction PSISP and
        ! velocity-potential CHISP.
        ! Input:  INTEGER K -- latitude number ( MU = MU(K) )
        !         COMPLEX PSISP(0:M1MAX,0:N1) -- spectral streamfunction
        !         COMPLEX CHISP(0:M1MAX,0:N1) -- spectral velocity-potential
        ! Output: REAL UPHK(K1MAX), VPHK(K1MAX) -- zonal and meridional
        !         components of velocity in physical space at latitude MU(K)
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              INTEGER K
              COMPLEX UPHK(K1MAX), VPHK(K1MAX) ! changed from REAL
              COMPLEX PSISP(0:M1MAX,0:N1), CHISP(0:M1MAX,0:N1)
              COMPLEX ULAK(0:M1MAX), VLAK(0:M1MAX)
        ! Compute velocity components in latitude space (in temporary storage)
              CALL UVLAK(ULAK,VLAK,PSISP,CHISP,K)
        ! Convert to physical space
              !CALL LA2PH(UPHK,ULAK) LA2PH commented out by Wim because it's designed to pass reals to complex (broken)
              !CALL LA2PH(VPHK,VLAK)
              RETURN
              END
         
        !========================================================================
              SUBROUTINE UVLAK(ULAK,VLAK,PSISP,CHISP,K)
        !========================================================================
        ! UVLAK calculates the zonal and meridional components in latitude
        ! space (ULAK, VLAK) of the horizontal velocity for a single latitude,
        ! given the spectral space representation of streamfunction PSISP and
        ! velocity-potential CHISP.
        ! Note: This is an auxiliary routine used by routines HVELK and HVEL)
        ! Input:  INTEGER K -- latitude number ( MU = MU(K) )
        !         COMPLEX PSISP(0:M1MAX,0:N1) -- spectral streamfunction
        !         COMPLEX CHISP(0:M1MAX,0:N1) -- spectral velocity-potential
        ! Output: COMPLEX ULAK(0:M1MAX), VLAK(0:M1MAX) -- zonal and meridional
        !         components of velocity in latitude space at latitude MU(K)
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
            !   INCLUDE 'sppoly.inc'
        
              INTEGER K
              COMPLEX ULAK(0:M1MAX), VLAK(0:M1MAX)
              COMPLEX PSISP(0:M1MAX,0:N1), CHISP(0:M1MAX,0:N1)
              COMPLEX FTMK1(0:M1MAX), FTMK2(0:M1MAX)
              REAL RTEMP1, RTEMP2
              INTEGER M
              RTEMP1= A0INV * COSINV(K)
              RTEMP2= A0INV * COSPHI(K)
        ! Calculate scalar "velocity" U*COS(PHI)
              CALL XSP2LA(FTMK1,CHISP,P2(0,0,K),H2(0,0,K),N1,K)
              CALL YSP2LA(FTMK2,PSISP,P2(0,0,K),H2(0,0,K),N1,K)
        
              DO 10 M=0,M1
                ULAK(M)= RTEMP1*FTMK1(M) - RTEMP2*FTMK2(M)
           10   CONTINUE
        ! Calculate scalar "velocity" V*COS(PHI)
              CALL XSP2LA(FTMK1,PSISP,P2(0,0,K),H2(0,0,K),N1,K)
              CALL YSP2LA(FTMK2,CHISP,P2(0,0,K),H2(0,0,K),N1,K)
        
              DO 20 M=0,M1
                VLAK(M)= RTEMP1*FTMK1(M) + RTEMP2*FTMK2(M)
           20   CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE DIVK(FSP,FPHK1,FPHK2,K)
        !========================================================================
        ! DIVK adds to spectral space representation FSP the contribution
        ! from the divergence of a physical space vector field whose zonal and
        ! meridional components at a single latitude MU=MU(K) are given by
        ! FPHK1 and FPHK2. (When contributions from all the latitudes have
        ! been added up, the result will be the complete spectral space
        ! representation of the divergence of the vector field.)
        ! Note: Truncation order is assumed to be N1 (i.e. as for a scalar)
        ! Input:  COMPLEX FSP(0:M1MAX,0:N1MAX) -- old spectral coefficients
        !         REAL FPHK1(K1MAX), FPHK2(K1MAX) -- physical space
        !         representation of vector components (zonal, meridional)
        !         INTEGER K -- latitude number ( MU = MU(K) )
        ! Output: COMPLEX FSP(0:M1MAX,0:N1MAX) -- Spectral representation of F,
        !            with contribution to the divergence from latitude MU added.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              COMPLEX FSP(0:M1MAX,0:N1MAX)
              REAL FPHK1(K1MAX), FPHK2(K1MAX)
              INTEGER K
              REAL RTEMP, TEMPHK(K1MAX)
              INTEGER J
              !EXTERNAL XLA2SP, YLA2SP
        ! Accumulate zonal-derivative of FPHK1 / COSPHI
              RTEMP= A0INV * COSINV(K)
              DO 10 J=1,K1
           10   TEMPHK(J)= RTEMP * FPHK1(J)
              CALL P2SK_XLA2SP( FSP, TEMPHK, K)
        ! Accumulate (1/A0) * MU-derivative of COSPHI * FPHK2
              RTEMP= A0INV * COSPHI(K)
              DO 20 J=1,K1
           20   TEMPHK(J)= RTEMP * FPHK2(J)
              CALL P2SK_YLA2SP( FSP, TEMPHK, K)
              RETURN
              END
         
        !========================================================================
              SUBROUTINE CURLK(FSP,FPHK1,FPHK2,K)
        !========================================================================
        ! CURLK adds to spectral space representation FSP the contribution
        ! from curl-z of a physical space vector field whose zonal and
        ! meridional components at a single latitude MU=MU(K) are given by
        ! FPHK1 and FPHK2. (When contributions from all the latitudes have
        ! been added up, the result will be the complete spectral space
        ! representation of the z-component of the curl of the vector field.)
        ! Note: Truncation order is assumed to be N1 (i.e. as for a scalar)
        ! Input:  COMPLEX FSP(0:M1MAX,0:N1MAX) -- old spectral coefficients
        !         REAL FPHK1(K1MAX), FPHK2(K1MAX) -- physical space
        !         representation of vector components (zonal, meridional)
        !         INTEGER K -- latitude number ( MU = MU(K) )
        ! Output: COMPLEX FSP(0:M1MAX,0:N1MAX) -- Spectral representation of F,
        !            with contribution to curl-z from latitude MU added.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              COMPLEX FSP(0:M1MAX,0:N1MAX)
              REAL FPHK1(K1MAX), FPHK2(K1MAX)
              INTEGER K
              REAL RTEMP, TEMPHK(K1MAX)
              INTEGER J
              !EXTERNAL XLA2SP, YLA2SP
        ! Accumulate zonal-derivative of FPHK2 / COSPHI
              RTEMP= A0INV * COSINV(K)
              DO 10 J=1,K1
           10   TEMPHK(J)= RTEMP * FPHK2(J)
              CALL P2SK_XLA2SP( FSP, TEMPHK, K)
        ! Accumulate (1/A0) * MU-derivative of -COSPHI * FPHK1
              RTEMP= -A0INV * COSPHI(K)
              DO 20 J=1,K1
           20   TEMPHK(J)= RTEMP * FPHK1(J)
              CALL P2SK_YLA2SP( FSP, TEMPHK, K)
              RETURN
              END
         
        !========================================================================
              SUBROUTINE SP2LA(FLAK,FSP,P0,H0,N0,K)
        !========================================================================
        ! SP2LA converts the spectral space representation FSP to latitude
        ! space representation FLAK, for a single latitude MU.
        ! Input:  COMPLEX FSP(0:M1MAX,0:N0) -- Spectral representation of F.
        !         REAL P0(0:2*M1MAX+1,0:N0), H0(0:2*M1MAX+1,0:N0) --
        !         (Legendre polynomials and their MU-derivatives at some
        !          latitude MU. Usually {PK2, HK2} or {P2(0,0,K), H2(0,0,K)} )
        !         INTEGER N0 -- order of truncation (N1 or N2).
        !         INTEGER K -- latitude number ( MU = MU(K) )
        ! Output: COMPLEX FLAK(0:M1MAX) -- latitude representation of F,
        !                                  at latitude MU.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
        ! Since we are using the strange index M2, pretend that the complex
        ! arrays are real.
              INTEGER N0,K
              COMPLEX FLAK(0:2*M1MAX+1),     FSP(0:2*M1MAX+1,0:N0)
              REAL   P0(0:2*M1MAX+1,0:N0), H0(0:2*M1MAX+1,0:N0)
              INTEGER M2,N
        ! Zero latitude space coefficients
              DO 10 M2=0,2*M1+1
                FLAK(M2)= 0.0
           10   CONTINUE
        ! Multiply by associated legendre polynomial, and accumulate
              IF (M1.GT.MLOW) THEN
                DO 20 N=0,N0
                  DO 20 M2=0,2*MIN0(N,M1)+1
                    FLAK(M2)= FLAK(M2) + FSP(M2,N) * P0(M2,N)
           20       CONTINUE
              ELSE
                DO 30 N=0,N0
           30     FLAK(0)= FLAK(0) + FSP(0,N) * P0(0,N)
        !
                DO 40 M2=2,2*M1+1
                  DO 40 N=M2/2,N0
                    FLAK(M2)= FLAK(M2) + FSP(M2,N) * P0(M2,N)
           40       CONTINUE
              ENDIF
              RETURN
              END
         
        !========================================================================
              SUBROUTINE XSP2LA(FLAK,FSP,P0,H0,N0,K)
        !========================================================================
        ! XSP2LA returns FLAK, the latitude space representation of the
        ! LAMBDA-Derivative of a function F, given FSP, the spectral space
        ! represention of F. All this for a single latitude MU.
        ! i.e. XSP2LA is SP2LA combined with d/dLAMBDA.
        ! Input:  COMPLEX FSP(0:M1MAX,0:N0) -- Spectral representation of F.
        !         REAL P0(0:2*M1MAX+1,0:N0), H0(0:2*M1MAX+1,0:N0) --
        !         (Legendre polynomials and their MU-derivatives at some
        !          latitude MU. Usually {PK2, HK2} or {P2(0,0,K), H2(0,0,K)} )
        !         INTEGER N0 -- order of truncation (N1 or N2).
        !         INTEGER K -- latitude number ( MU = MU(K) )
        ! Output: COMPLEX FLAK(0:M1MAX) -- latitude representation of
        !                         LAMBDA-derivative of F, at latitude MU.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              INTEGER N0,K
              COMPLEX FLAK(0:M1MAX), FSP(0:M1MAX,0:N0)
              REAL P0(0:2*M1MAX+1,0:N0), H0(0:2*M1MAX+1,0:N0)
              INTEGER M
              IF (M1.GT.0) THEN
        ! Do a normal transform from spectral space to latitude space
                CALL SP2LA(FLAK,FSP,P0,H0,N0,K)
        ! Evaluate lambda-derivative
                DO 10 M=0,M1
                  FLAK(M)= CIM(M) * FLAK(M)
           10     CONTINUE
              ELSE
                FLAK(0)= (0., 0.)
              ENDIF
              RETURN
              END
         
        !========================================================================
              SUBROUTINE YSP2LA(FLAK,FSP,P0,H0,N0,K)
        !========================================================================
        ! YSP2LA returns FLAK, the latitude space representation of the
        ! MU-derivative of a function F, given FSP, the spectral space
        ! representation of F. All this for a single latitude MU.
        ! i.e. YSP2LA is SP2LA combined with d/dMU.
        ! Input:  COMPLEX FSP(0:M1MAX,0:N0) -- Spectral representation of F.
        !         REAL P0(0:2*M1MAX+1,0:N0), H0(0:2*M1MAX+1,0:N0) --
        !         (Legendre polynomials and their MU-derivatives at some
        !          latitude MU. Usually {PK2, HK2} or {P2(0,0,K), H2(0,0,K)} )
        !         INTEGER N0 -- order of truncation (N1 or N2).
        !         INTEGER K -- latitude number ( MU = MU(K) )
        ! Output: COMPLEX FLAK(0:M1MAX) -- latitude representation of
        !                             MU-derivative of F, at some latitude MU.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
        ! Since we are using the strange index M2, pretend that the complex
        ! arrays are real.
              INTEGER N0,K
              COMPLEX FLAK(0:2*M1MAX+1),     FSP(0:2*M1MAX+1,0:N0)
              REAL   P0(0:2*M1MAX+1,0:N0), H0(0:2*M1MAX+1,0:N0)
              INTEGER M2,N
              REAL RTEMP
        ! Zero latitude space coefficients
              DO 10 M2=0,2*M1+1
                FLAK(M2)= 0.0
           10   CONTINUE
        ! Multiply by -(1-MU**2) times the MU-derivative of associated
        ! legendre polynomial, and accumulate.
              IF (M1.GT.MLOW) THEN
                DO 20 N=0,N0
                  DO 20 M2=0,2*MIN0(N,M1)+1
                    FLAK(M2)= FLAK(M2) + FSP(M2,N) * H0(M2,N)
           20       CONTINUE
              ELSE
                DO 30 N=0,N0
           30     FLAK(0)= FLAK(0) + FSP(0,N) * H0(0,N)
        !
                DO 40 M2=2,2*M1+1
                  DO 40 N=M2/2,N0
                    FLAK(M2)= FLAK(M2) + FSP(M2,N) * H0(M2,N)
           40       CONTINUE
              ENDIF
        ! Divide by -(1-MU**2) to obtain MU-derivative
              RTEMP= - (COSINV(K)**2)
              DO 50 M2=0,2*M1+1
           50   FLAK(M2)= RTEMP * FLAK(M2)
              RETURN
              END
         
        !========================================================================
              SUBROUTINE LA2SP(FSP,FLAK,P0,H0,N0,K)
        !========================================================================
        ! LA2SP adds to spectral space representation FSP the contribution
        ! from latitude space representation FLAK, for a single latitude
        ! MU. (When contributions from all the latitudes have been added
        ! up, the result will be the complete spectral space representation.)
        ! Input:  COMPLEX FSP(0:M1MAX,0:N0) -- Spectral representation of F.
        !         COMPLEX FLAK(0:M1MAX) -- Latitude representation of F,
        !                                  at some latitude MU.
        !         REAL P0(0:2*M1MAX+1,0:N0), H0(0:2*M1MAX+1,0:N0) --
        !         (Legendre polynomials and their MU-derivatives at some
        !          latitude MU. Usually {PK2, HK2} or {P2(0,0,K), H2(0,0,K)} )
        !         INTEGER N0 -- order of truncation (N1 or N2).
        !         INTEGER K -- latitude number ( MU = MU(K) )
        ! Output: COMPLEX FSP(0:M1MAX,0:N0) -- Spectral representation of F,
        !                           with contribution from latitude MU added.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
        ! Since we are using the strange index M2, pretend that the complex
        ! arrays are real.
              INTEGER N0,K
              COMPLEX FSP(0:2*M1MAX+1,0:N0), FLAK(0:2*M1MAX+1)
              REAL  P0(0:2*M1MAX+1,0:N0),   H0(0:2*M1MAX+1,0:N0)
              INTEGER M2,N
              REAL GFLAK(0: 2*M1MAX+1)
        ! Evaluate temporary array
              DO 10 M2=0,2*M1+1
                GFLAK(M2)= G(K) * FLAK(M2)
           10   CONTINUE
        ! Evaluate gaussian quadrature contribution from latitude K.
        ! Remember, the factor of 1/2 is included in G(K)
              IF (M1.GT.MLOW) THEN
                DO 20 N=0,N0
                  DO 20 M2=0,2*MIN0(N,M1)+1
                    FSP(M2,N)= FSP(M2,N) + GFLAK(M2) * P0(M2,N)
           20       CONTINUE
              ELSE
                DO 30 N=0,N0
           30     FSP(0,N)= FSP(0,N) + GFLAK(0) * P0(0,N)
        !
                DO 40 M2=2,2*M1+1
                  DO 40 N=M2/2,N0
                    FSP(M2,N)= FSP(M2,N) + GFLAK(M2) * P0(M2,N)
           40       CONTINUE
              ENDIF
              RETURN
              END
         
        !========================================================================
              SUBROUTINE XLA2SP(FSP,FLAK,P0,H0,N0,K)
        !========================================================================
        ! XLA2SP does almost the same thing as LA2SP. The only difference
        ! is that instead of transforming the function F (given by FLAK), it
        ! transforms the LAMBDA-derivative of F to spectral space. i.e. It
        ! returns (accumulates) (0.,1.)*M times the result given by LA2SP.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              INTEGER N0,K
              COMPLEX FSP(0:M1MAX,0:N0), FLAK(0:M1MAX)
              REAL P0(0:2*M1MAX+1,0:N0), H0(0:2*M1MAX+1,0:N0)
              COMPLEX FTMKC(0:M1MAX)
              INTEGER M
              IF (M1.GT.0) THEN
        ! Evaluate temporary array (taking LAMBDA-derivative)
                DO 10 M=0,M1
           10     FTMKC(M)= CIM(M) * FLAK(M)
        ! Call LA2SP to do the rest of the job
                CALL LA2SP(FSP,FTMKC,P0,H0,N0,K)
              ENDIF
              RETURN
              END
         
        !========================================================================
              SUBROUTINE YLA2SP(FSP,FLAK,P0,H0,N0,K)
        !========================================================================
        ! YLA2SP does almost the same thing as LA2SP. the only difference
        ! is that instead of transforming the function F (given by FLAK), it
        ! transforms the MU-derivative of F to spectral space. It does this
        ! by integrating by parts, assuming the boundary terms vanish, and
        ! using H[M,N] instead of P[M,N] in the transform.
        ! i.e. YLA2SP returns (accumulates) dF/dMU for a single latitude.
        ! **CAUTION** YLA2SP will work only for functions which vanish at the
        !             poles. e.g. COS(PHI) * some-function
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
        ! Since we are using the strange index M2, pretend that the complex
        ! arrays are real.
              INTEGER N0,K
              COMPLEX FSP(0:2*M1MAX+1,0:N0), FLAK(0:2*M1MAX+1) ! changed from real
              REAL  P0(0:2*M1MAX+1,0:N0),   H0(0:2*M1MAX+1,0:N0)
              INTEGER M2,N
              REAL GFLAK(0: 2*M1MAX+1), RTEMP
        ! Evaluate temporary array (incorporates division by (1-MU**2) )
              RTEMP= G(K) * (COSINV(K)**2)
              DO 10 M2=0,2*M1+1
                GFLAK(M2)= RTEMP * FLAK(M2)
           10   CONTINUE
        ! Evaluate gaussian quadrature contribution from latitude K
        ! Remember, the factor of 1/2 is included in G(K)
              IF (M1.GT.MLOW) THEN
                DO 20 N=0,N0
                  DO 20 M2=0,2*MIN0(N,M1)+1
                    FSP(M2,N)= FSP(M2,N) + GFLAK(M2) * H0(M2,N)
           20       CONTINUE
              ELSE
                DO 30 N=0,N0
           30     FSP(0,N)= FSP(0,N) + GFLAK(0) * H0(0,N)
        !
                DO 40 M2=2,2*M1+1
                  DO 40 N=M2/2,N0
                    FSP(M2,N)= FSP(M2,N) + GFLAK(M2) * H0(M2,N)
           40       CONTINUE
              ENDIF
              RETURN
              END
         
        ! !========================================================================
        !       SUBROUTINE LA2PH(FPHK,FLAK)
        ! !========================================================================
        ! ! LA2PH converts a function F from latitude space representation FLAK
        ! ! to physical space representation FPHK, for any single latitude.
        ! !
        ! ! Note: LA2PH will not handle the case K1 = 1 (axially symmetric)
        ! ! Input:  COMPLEX FLAK(0:M1MAX) -- latitude space representation of F
        ! !                                  for some particular latitude.
        ! ! Output: REAL FPHK(K1MAX) -- physical space representation of F,
        ! !                             for some particular latitude.
        ! !
        
        !       IMPLICIT NONE
        !     !   INCLUDE 'spcons.inc'
        !     !   INCLUDE 'spgrid.inc'
        !     !   INCLUDE 'spfftb.inc'
        
        ! ! Although FPHK is real, pretend that it is complex
        !       COMPLEX FPHK( (K1MAX+1)/2 ),FLAK(0:M1MAX)
        !       INTEGER M
        !       IF (K1.LT.2) CALL SPERR( 'LA2PH', 'K1 less than 2')
        ! ! Copy the first M1+1 fourier coefficients to destination array
        ! ! (Take complex conjugate because the discrete fourier transform
        ! ! according to REALF2 is defined the "wrong" way)
        ! ! (The factor of 2 comes in because the inverse fourier transform
        ! ! according to REALF2 is divided by K1/2 and not by K1. Note that
        ! ! subroutine PH2LA divides by a factor K1)
        !       DO 10 M=0,M1
        !         FPHK(M+1)= CONJG( FLAK(M) ) * 2.
        !    10   CONTINUE
        ! ! Zero all extra coefficients, extending M1+1 values to K1/2 values.
        ! ! Since we are avoiding aliasing, it is safe to assume K1/2 >= M1+1
        !       DO 20 M=M1+1,K1/2 - 1
        !         FPHK(M+1)= (0.,0.)
        !    20   CONTINUE
        ! ! Since FLAK is the fourier transform of real function FPHK,
        ! ! invert it.
        !       CALL REALF2(FPHK,-1,FFTAB,NFFTAB,K1/2)
        !       RETURN
        !       END
         
        !========================================================================
              SUBROUTINE PH2LA(FLAK,FPHK)
        !========================================================================
        ! PH2LA converts a function F from physical space representation FPHK
        ! to latitude space representation FLAK, for any single latitude.
        !
        ! Note: PH2LA will not handle the case K1 = 1 (axially symmetric)
        ! Input:  REAL FPHK(K1MAX) -- physical space representation of F,
        !                             for some particular latitude.
        ! Output: COMPLEX FLAK(0:M1MAX) -- latitude space representation of F
        !                                  for some particular latitude.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
            !   INCLUDE 'spfftb.inc'
        
              REAL FPHK(K1MAX)
              COMPLEX FLAK(0:M1MAX)
              REAL FTMKR(K1MAX)
              COMPLEX FTMKC( (K1MAX+1)/2 )
              EQUIVALENCE (FTMKR,FTMKC)
              INTEGER J,M
              IF (K1.LT.2) CALL SPERR( 'PH2LA', 'K1 less than 2')
        ! Copy physical values into temporary storage
              DO 10 J=1,K1
           10   FTMKR(J)= FPHK(J)
        ! Fourier transform of real function FPHK
              CALL REALF2(FTMKR,+1,FFTAB,NFFTAB,K1/2)
        ! Retrieve the real first fourier coefficient
        ! (Divide by K1 because the discrete fourier transform according to
        ! REALF2 does not incorporate this factor)
              FLAK(0)= CMPLX( FTMKR(1), 0. ) / K1
        ! Retrieve the next M1 complex fourier coefficients
        ! (Take complex conjugate because the discrete fourier transform
        !  according to REALF2 is defined the "wrong" way)
              DO 20 M=1,M1
           20   FLAK(M)= CONJG( FTMKC(M+1) ) / K1
              RETURN
              END
         
        !========================================================================
              SUBROUTINE ZEROSP(FSP,N0)
        !========================================================================
        ! ZEROSP zeroes the spectral coefficient matrix FSP, for truncation
        ! order N0.
        ! Input:  INTEGER N0 -- order of truncation (N1 or N2).
        ! Output: COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              INTEGER N0
              COMPLEX FSP(0:M1MAX,0:N0)
              INTEGER M,N
              DO 10 N=0,N0
                DO 10 M=0,MIN0(N,M1)
                  FSP(M,N)= (0.,0.)
           10     CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE SPCOPY(FSP2,FSP1,N0)
        !========================================================================
        ! SPCOPY copies spectral coefficient matrix FSP1 to FSP2
        ! Input:  INTEGER N0 -- order of truncation (N1 or N2).
        !         COMPLEX FSP1(0:M1MAX,0:N0) -- spectral representation of F.
        ! Output: COMPLEX FSP2(0:M1MAX,0:N0) -- copy of FSP1.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              INTEGER N0
              COMPLEX FSP2(0:M1MAX,0:N0),FSP1(0:M1MAX,0:N0)
              INTEGER M,N
              DO 10 N=0,N0
                DO 10 M=0,MIN0(N,M1)
                  FSP2(M,N)= FSP1(M,N)
           10     CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE DELSQ(DSQFSP,FSP,N0)
        !========================================================================
        ! DELSQ calculates DSQFSP, del-squared (Laplacian) of a function F,
        ! given its spectral representation FSP.
        ! Input:  INTEGER N0 -- order of truncation (N1 or N2).
        !         COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F.
        ! Output: COMPLEX DSQFSP(0:M1MAX,0:N0) -- spectral representation
        !                                         of del-squared F.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              INTEGER N0
              COMPLEX DSQFSP(0:M1MAX,0:N0),FSP(0:M1MAX,0:N0)
              INTEGER M,N
              REAL TEMP
              DO 10 N=0,N0
                TEMP= -N * (N+1) * A0QINV
                DO 10 M=0,MIN0(N,M1)
                  DSQFSP(M,N)= FSP(M,N) * TEMP
           10     CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE IDELSQ(FSP,DSQFSP,N0)
        !========================================================================
        ! IDELSQ calculates FSP, the spectral representation of a function F,
        ! given the spectral representation of del-squared (Laplacian) of F.
        ! (Note: Since the inverse laplacian is not unique, FSP(0,0) is always
        !        set to zero)
        ! Input:  INTEGER N0 -- order of truncation (N1 or N2).
        !         COMPLEX DSQFSP(0:M1MAX,0:N0) -- spectral representation OF
        !                                         del-squared F.
        ! Output: COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              INTEGER N0
              COMPLEX FSP(0:M1MAX,0:N0),DSQFSP(0:M1MAX,0:N0)
              INTEGER M,N
              REAL TEMP
        ! Inverse laplacian is non-unique to upto a constant
              FSP(0,0)= (0.,0.)
              DO 10 N=1,N0
                TEMP= -A0Q / (N * (N+1) )
                DO 10 M=0,MIN0(N,M1)
                  FSP(M,N)= DSQFSP(M,N) * TEMP
           10     CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE SPDLM(DLMFSP,FSP,N0)
        !========================================================================
        ! SPDLM calculates DLMFSP, the LAMBDA-derivative of a function F,
        ! given its spectral representation FSP.
        ! Input:  INTEGER N0 -- order of truncation (N1 or N2).
        !         COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F.
        ! Output: COMPLEX DLMFSP(0:M1MAX,0:N0) -- spectral representation
        !                                         of d/dLAMBDA of F.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              INTEGER N0
              COMPLEX DLMFSP(0:M1MAX,0:N0),FSP(0:M1MAX,0:N0)
              INTEGER M,N
              DO 10 N=0,N0
                DO 10 M=0,MIN0(N,M1)
                  DLMFSP(M,N)= CIM(M) * FSP(M,N)
           10     CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE HSPDMU(DMUFSP,FSP,N0)
        !========================================================================
        ! HSPDMU calculates DLMFSP, -(1-MU**2) times the MU-derivative of
        ! a function F, given its spectral representation FSP.
        ! Input:  INTEGER N0 -- order of truncation (N1 or N2).
        !         COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F.
        ! **CAUTION** The input array FSP should be physically different
        !             from the output array DMUFSP
        ! Output: COMPLEX DMUFSP(0:M1MAX,0:N0) -- spectral representation
        !                 of -(1-MU**2) times d/dMU of F.
        ! NOTE: Since the -(1-MU**2)*d/dMU operation increases the order of
        !       truncated series by 1, HSPDMU should typically be called with
        !       N0=N2 for a scalar F (of truncation order N1), after setting
        !       the extra coefficients to zero. This way DMUFSP will contain
        !       all the terms of the derivative. Otherwise, terms may be lost.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              INTEGER N0
              COMPLEX DMUFSP(0:M1MAX,0:N0),FSP(0:M1MAX,0:N0)
              INTEGER M,N
              IF (M1.EQ.N0) THEN
                DO 20 M=0,M1-1
                  DMUFSP(M,M)=   - (M+2) * D(M,M+1)  * FSP(M,M+1)
                  DO 10 N=M+1,N0-1
                    DMUFSP(M,N)=   (N-1) * D(M,N)    * FSP(M,N-1) &
                                 - (N+2) * D(M,N+1)  * FSP(M,N+1)
           10       CONTINUE
                  DMUFSP(M,N0)=   (N0-1) * D(M,N0)   * FSP(M,N0-1)
           20     CONTINUE
                  DMUFSP(M1,N0)= (0.,0.)
              ELSE
                DO 40 M=0,M1
                  DMUFSP(M,M)=   - (M+2) * D(M,M+1)  * FSP(M,M+1)
                  DO 30 N=M+1,N0-1
                    DMUFSP(M,N)=   (N-1) * D(M,N)    * FSP(M,N-1) &
                                 - (N+2) * D(M,N+1)  * FSP(M,N+1)
           30       CONTINUE
                  DMUFSP(M,N0)=   (N0-1) * D(M,N0)   * FSP(M,N0-1)
           40     CONTINUE
              ENDIF
              RETURN
              END
         
        !========================================================================
              SUBROUTINE GAULEG(X1,X2,X,W,N)
        !========================================================================
        !
        ! Given the lower and upper limits of integration X1 and X2, and
        ! given N, this routine returns arrays X and W of length N, containing
        ! the abcissas and weights of the Gauss-Legendre n-point quadrature
        ! formula.
        ! (Numerical Recipes - W.H.Press et al, Cambridge(1986), p.125-126)
        !
        ! High precision is a good idea for this routine.
              IMPLICIT DOUBLE PRECISION (A-H,O-Z)
              IMPLICIT INTEGER (I-N)
              REAL X1,X2,X(N),W(N), P1,P2,P3
              PARAMETER (EPS=3.D-14)
        ! Increase if you don't have this precision.
        !
              M=(N+1)/2
        ! The roots are symmetric in the interval, so we only find half of them.
        !
              XM=0.5D0*(X2+X1)
              XL=0.5D0*(X2-X1)
        !
        ! Loop over the desired roots
              DO 12 I=1,M
                Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0)) 
        ! Starting with the above approximation to the I-th root, we enter the
        ! main loop of refinement by Newton's method.
        !
            1 CONTINUE ! 1 preceding here is for GO TO 1 below
                P1=1.D0
                P2=0.D0
        !
        ! Loop up the recurrence relation to get the Legendre polynomial
        ! evaluated at Z.
                DO 11 J=1,N
                  P3=P2
                  P2=P1
                  P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
           11     CONTINUE
        ! P1 is the now desired Legendre polynomial. We next compute PP, its
        ! derivative, by a standard relation involving also P2, the polynomial
        ! of one lower order.
        !
                PP=N*(Z*P1-P2)/(Z*Z-1.D0)
                Z1=Z
                Z=Z1-P1/PP
        ! Newton's method.
        !
                IF(ABS(Z-Z1).GT.EPS)GO TO 1
                X(I)=XM-XL*Z
                X(N+1-I)=XM+XL*Z
        ! Scale the root to the desired interval
        ! and put in its symmetric counterpart.
        !
                W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
                W(N+1-I)=W(I)
        ! Compute the weight and its symmetric counterpart.
        !
           12   CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE FOUR0(FFTAB, MAXTAB, NN)
        !========================================================================
        !
        ! Definitions:
        ! H(J) :  Discrete fourier transform of N points H(K)
        ! H(J) =  Sigma(K=0 TO N-1) H(K)*EXP(I*K*J*2*PI/N)
        ! H(K) :  Inverse discrete fourier transform of N points H(J)
        ! H(K) =  (1/N) * Sigma(J=0 TO N-1) H(J)*EXP(-I*K*J*2*PI/N)
        ! If H(K) is real, then the discrete fourier transform H(J)
        ! has the property --  H(N-J)=CONJG(H(J)), N=number of points.
        ! also, H(J=0) and H(J=N/2) are real and independent.
        !
        ! FOUR0 initializes a fast fourier transform table (FFTAB) used for
        ! transforming a complex array of length NN or, equivalently, a real
        ! array of length 2*NN. This table is used by routines FOUR2, FOURV,
        ! REALF2, and REALFV.
        ! Input:  INTEGER MAXTAB -- Declared outer dimension of array FFTAB
        !         INTEGER NN -- Length of complex array to be transformed
        !         (NN must be an integer power of 2; this is checked for)
        ! Output: REAL FFTAB(4,2,MAXTAB) -- Table of coefficients
        !         (Determining the size MAXTAB is not easy. A safe bet may
        !         be MAXTAB=256. FOUR0 aborts with an error message if NN
        !         is too large for a given MAXTAB. If that happens, increase
        !         the value of MAXTAB)
        ! (Adapted from Numerical Recipes - W.H.Press et al, Cambridge(1986))
        !
              !IMPLICIT REAL (A-H,O-Z)
              REAL A,B,C,D,E,F,G,H
              REAL O,P,Q,R,S,T,U,V,W,X,Y,Z
              !IMPLICIT INTEGER (I-N)
              INTEGER I,J,K,L,M,N
            !   INCLUDE 'ftable.inc'
              real, dimension(4,2,MAXTAB) :: FFTAB
              integer :: maxtab
              INTEGER NN
              INTEGER ISIGN, ISGN2, IW, MMAX, ISTEP
              DOUBLE PRECISION WR,WI,WPR,WPI,WTEMP,THETA
        ! Double precision for trigonometric recurrences.
        !
              N=2*NN
              J=1
        !
        ! Check whether NN is a power of two ????
            1 J= J*2
              IF (J.GT.NN) CALL SPERR('FOUR0', 'NN must be a power of 2')
              IF (J.LT.NN) GOTO 1
        ! Here begins the Danielson-Lanczos section of the initialization.
        ! Outer loop executed log (base 2) NN times. (Forward transform)
              ISIGN= 1
              ISGN2= 1
              MMAX=2
              IW= 1
          110 IF (N.GT.MMAX) THEN
                ISTEP=2*MMAX
                THETA=6.28318530717959D0/(ISIGN*MMAX)
        ! Initialize for trigonometric recurrence.
        !
                WPR=-2.D0*DSIN(0.5D0*THETA)**2
                WPI=DSIN(THETA)
                WR=1.D0
                WI=0.D0
                DO 115 M=1,MMAX,2
                  IF (IW.GT.MAXTAB) CALL SPERR('FOUR0', 'NN too large')
                  FFTAB(1,ISGN2,IW)= SNGL(WR)
                  FFTAB(2,ISGN2,IW)= SNGL(WI)
                  IW= IW+1
        ! Trigonometric recurrence.
                  WTEMP=WR
                  WR=WR*WPR-WI*WPI+WR
                  WI=WI*WPR+WTEMP*WPI+WI
          115     CONTINUE
                MMAX=ISTEP
                GOTO 110
        ! Not yet done.
              END IF
        ! All done.
        ! Outer loop executed log (base 2) NN times. (Inverse transform)
              ISIGN= -1
              ISGN2= 2
              MMAX=2
              IW= 1
          120 IF (N.GT.MMAX) THEN
                ISTEP=2*MMAX
                THETA=6.28318530717959D0/(ISIGN*MMAX)
        ! Initialize for trigonometric recurrence.
        !
                WPR=-2.D0*DSIN(0.5D0*THETA)**2
                WPI=DSIN(THETA)
                WR=1.D0
                WI=0.D0
                DO 125 M=1,MMAX,2
                  IF (IW.GT.MAXTAB) CALL SPERR('FOUR0', 'NN too large')
                  FFTAB(1,ISGN2,IW)= SNGL(WR)
                  FFTAB(2,ISGN2,IW)= SNGL(WI)
                  IW= IW+1
        ! Trigonometric recurrence.
                  WTEMP=WR
                  WR=WR*WPR-WI*WPI+WR
                  WI=WI*WPR+WTEMP*WPI+WI
          125     CONTINUE
                MMAX=ISTEP
                GOTO 120
        ! Not yet done.
              END IF
        ! All done.
        ! Initialize trigonometric coefficients for the real fourier transform
        ! (Forward transform)
              ISIGN= 1
              ISGN2= 1
              THETA=3.141592653589793D0/DBLE(ISIGN*NN)
              WPR=-2.0D0*DSIN(0.5D0*THETA)**2
              WPI=DSIN(THETA)
              WR=1.0D0+WPR
              WI=WPI
              IW= 1
              DO 210 I=2,NN/2+1
                IF (IW.GT.MAXTAB) CALL SPERR('FOUR0', 'NN too large')
                FFTAB(3,ISGN2,IW)= SNGL(WR)
                FFTAB(4,ISGN2,IW)= SNGL(WI)
                IW= IW+1
        !
        ! The recurrence.
                WTEMP=WR
                WR=WR*WPR-WI*WPI+WR
                WI=WI*WPR+WTEMP*WPI+WI
          210   CONTINUE
        !
        ! (Inverse transform)
              ISIGN= -1
              ISGN2= 2
              THETA=3.141592653589793D0/DBLE(ISIGN*NN)
              WPR=-2.0D0*DSIN(0.5D0*THETA)**2
              WPI=DSIN(THETA)
              WR=1.0D0+WPR
              WI=WPI
              IW= 1
              DO 220 I=2,NN/2+1
                IF (IW.GT.MAXTAB) CALL SPERR('FOUR0', 'NN too large')
                FFTAB(3,ISGN2,IW)= SNGL(WR)
                FFTAB(4,ISGN2,IW)= SNGL(WI)
                IW= IW+1
        !
        ! The recurrence.
                WTEMP=WR
                WR=WR*WPR-WI*WPI+WR
                WI=WI*WPR+WTEMP*WPI+WI
          220   CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE FOUR2(DATA, ISIGN, FFTAB, MAXTAB, NN)
        !========================================================================
        ! FOUR2 is a table-driven version of FOUR1 (for fixed length arrays).
        ! FOUR2 replaces DATA by its discrete fourier transform, if ISIGN is
        ! input as 1; or replaces DATA by NN times its inverse discrete fourier
        ! transform, if ISIGN is input as -1. DATA is a complex array of length
        ! NN or, equivalently, a real array of length 2*NN. FFTAB is a table
        ! of coefficients initialized by routine FOUR0, for complex transform
        ! length NN. (See routine FOUR0 for details)
        ! (Adapted from Numerical Recipes - W.H.Press et al, CAMBRIDGE(1986))
        !
              IMPLICIT REAL (A-H,O-Z)
              IMPLICIT INTEGER (I-N)
            !   INCLUDE 'ftable.inc'
        
              real, dimension(4,2,MAXTAB) :: FFTAB
              
              INTEGER NN
              !DIMENSION DATA(2*NN)
              real, dimension(2*NN) :: DATA
              REAL WRS, WIS
              INTEGER IW, ISGN2
        ! Initialize table look-up index
              IF (ISIGN.EQ.1) THEN
                ISGN2= 1
              ELSE
                ISGN2= 2
              ENDIF
              N=2*NN
              J=1
        !
        ! This is the bit reversal section of the routine.
              DO 11 I=1,N,2
                IF (J.GT.I) THEN
        ! Exchange two complex numbers.
                  TEMPR=DATA(J)
                  TEMPI=DATA(J+1)
                  DATA(J)=DATA(I)
                  DATA(J+1)=DATA(I+1)
                  DATA(I)=TEMPR
                  DATA(I+1)=TEMPI
                END IF
                M=N/2
            1   IF ((M.GE.2).AND.(J.GT.M)) THEN
                  J=J-M
                  M=M/2
                  GO TO 1
                END IF
                J=J+M
           11   CONTINUE
        ! Here begins the Danielson-Lanczos section of the routine.
        ! Outer loop executed log (base 2) NN times.
              MMAX=2
        ! Initialize trigonometric table index
        !
              IW= 1
            2 IF (N.GT.MMAX) THEN
                ISTEP=2*MMAX
                DO 13 M=1,MMAX,2
                  WRS= FFTAB(1,ISGN2,IW)
                  WIS= FFTAB(2,ISGN2,IW)
                  DO 12 I=M,N,ISTEP
        ! This is the Danielson-Lanczos formula.
        !
                    J=I+MMAX
                    TEMPR=WRS*DATA(J)-WIS*DATA(J+1)
                    TEMPI=WRS*DATA(J+1)+WIS*DATA(J)
                    DATA(J)=DATA(I)-TEMPR
                    DATA(J+1)=DATA(I+1)-TEMPI
                    DATA(I)=DATA(I)+TEMPR
                    DATA(I+1)=DATA(I+1)+TEMPI
           12       CONTINUE
        ! Increment trigonometric table index
                  IW= IW+1
           13     CONTINUE
                MMAX=ISTEP
                GOTO 2
        ! Not yet done.
              END IF
        ! All done.
              RETURN
              END
         
        !========================================================================
              SUBROUTINE REALF2(DATA, ISIGN, FFTAB, MAXTAB, N)
        !========================================================================
        ! REALF2 is a table-driven version of REALFT (for fixed length arrays).
        ! REALF2 calculates the fourier transform of a set of 2*N real-valued
        ! data points. It replaces this data (which is stored in array DATA)
        ! by the positive frequency half of its complex fourier transform. The
        ! real valued first and last components of the complex transform are
        ! returned as elements DATA(1) and DATA(2) respectively. N must  be a
        ! power of 2. This routine also calculates the inverse transform of a
        ! complex data array if it is the transform of real data. (Result in
        ! this case must be multiplied by 1/N). FFTAB is a table of
        ! coefficients initialized by routine FOUR0, for complex transform
        ! length N. (See routine FOUR0 for details)
        ! (Adapted from Numerical Recipes - W.H.Press et al, Cambridge(1986))
        !
              IMPLICIT REAL (A-H,O-Z)
              IMPLICIT INTEGER (I-N)
            !   INCLUDE 'ftable.inc'
        
              real, dimension(4,2,MAXTAB) :: FFTAB
              INTEGER N
              !DIMENSION DATA(2*N)
              real, dimension(2*N) :: DATA
              INTEGER IW, ISGN2
        ! Initialize the recurrence
              C1=0.5
              IF (ISIGN.EQ.1) THEN
                ISGN2= 1
                C2=-0.5
                CALL FOUR2(DATA,+1,FFTAB,MAXTAB,N)
        ! The forward transform is here
        !
              ELSE
        ! Otherwise set up for an inverse transform
        !
                ISGN2= 2
                C2=0.5
              END IF
        ! Initialize trigonometric table lookup index
              IW= 1
              N2P3=2*N+3
              DO 11 I=2,N/2+1
        ! Case I=1 done separately below
                I1=2*I-1
                I2=I1+1
                I3=N2P3-I2
                I4=I3+1
                WRS= FFTAB(3,ISGN2,IW)
                WIS= FFTAB(4,ISGN2,IW)
        ! The two separate transforms are separated out of Z.
        !
                H1R= C1*(DATA(I1)+DATA(I3))
                H1I= C1*(DATA(I2)-DATA(I4))
                H2R=-C2*(DATA(I2)+DATA(I4))
                H2I= C2*(DATA(I1)-DATA(I3))
        !
        ! Here they are recombined to form the true transform of the original
        ! real data.
                DATA(I1)= H1R+WRS*H2R-WIS*H2I
                DATA(I2)= H1I+WRS*H2I+WIS*H2R
                DATA(I3)= H1R-WRS*H2R+WIS*H2I
                DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        !
        ! Increment the trigonometric table look-up index
                IW= IW+1
           11   CONTINUE
              IF (ISIGN.EQ.1) THEN
                H1R=DATA(1)
                DATA(1)=H1R+DATA(2)
                DATA(2)=H1R-DATA(2)
        ! Squeeze the first and last data together to get them all within the
        ! original array.
        !
              ELSE
        ! This is the inverse transform for the case ISIGN=-1.
        !
                H1R=DATA(1)
                DATA(1)=C1*(H1R+DATA(2))
                DATA(2)=C1*(H1R-DATA(2))
                CALL FOUR2(DATA,-1,FFTAB,MAXTAB,N)
              END IF
              RETURN
              END
         
        !========================================================================
              SUBROUTINE FOURV(DATA,ISIGN,NDAT0,NDAT,FFTAB,MAXTAB,NN)
        !========================================================================
        ! FOURV is a "vectorizable" version of FOUR2 (for multiple transforms).
        ! FOURV does the same thing as FOUR2, but for NDAT arrays
        ! simultaneously. i.e. DATA is dimensioned as COMPLEX (NDAT0,NN),
        ! where NDAT0 is the declared inner dimesion, and NDAT is the operative
        ! inner dimension. FFTAB is a table of coefficients initialized by
        ! routine FOUR0, for complex transform length NN. (See routine FOUR0
        ! for details)
        ! (Adapted from Numerical Recipes - W.H.Press et al, Cambridge(1986))
        !
              IMPLICIT REAL (A-H,O-Z)
              IMPLICIT INTEGER (I-N)
            !   INCLUDE 'ftable.inc'
              real, dimension(4,2,MAXTAB) :: FFTAB
              INTEGER NN
              DIMENSION DATA(NDAT0,2*NN)
              REAL WRS, WIS, TEMPR(MAXDAT), TEMPI(MAXDAT)
              REAL TEMPR2(MAXDAT), TEMPI2(MAXDAT)
              INTEGER IW, ISGN2, IDAT
        ! Check if temporary storage is sufficient
              IF (NDAT.GT.MAXDAT) CALL SPERR('FOURN', &
               'NDAT too large; recompile with increased MAXDAT')
        ! Initialize table look-up index
              IF (ISIGN.EQ.1) THEN
                ISGN2= 1
              ELSE
                ISGN2= 2
              ENDIF
              N=2*NN
              J=1
        !
        ! This is the bit reversal section of the routine.
              DO 11 I=1,N,2
                IF (J.GT.I) THEN
        ! Exchange two complex numbers.
                  DO 4 IDAT=1,NDAT
                    TEMPR(IDAT)=DATA(IDAT,J)
                    TEMPI(IDAT)=DATA(IDAT,J+1)
                    TEMPR2(IDAT)=DATA(IDAT,I)
            4       TEMPI2(IDAT)=DATA(IDAT,I+1)
                  DO 5 IDAT=1,NDAT
                    DATA(IDAT,I)=TEMPR(IDAT)
            5       DATA(IDAT,I+1)=TEMPI(IDAT)
                  DO 6 IDAT=1,NDAT
                    DATA(IDAT,J)=TEMPR2(IDAT)
            6       DATA(IDAT,J+1)=TEMPI2(IDAT)
                END IF
                M=N/2
            1   IF ((M.GE.2).AND.(J.GT.M)) THEN
                  J=J-M
                  M=M/2
                  GO TO 1
                END IF
                J=J+M
           11   CONTINUE
        ! Here begins the Danielson-Lanczos section of the routine.
        ! Outer loop executed log (base 2) NN times.
              MMAX=2
        ! Initialize trigonometric table index
        !
              IW= 1
            2 IF (N.GT.MMAX) THEN
                ISTEP=2*MMAX
                DO 13 M=1,MMAX,2
                  WRS= FFTAB(1,ISGN2,IW)
                  WIS= FFTAB(2,ISGN2,IW)
                  DO 12 I=M,N,ISTEP
        ! This is the Danielson-Lanczos formula.
        !
                    J=I+MMAX
                    DO 7 IDAT=1,NDAT
                      TEMPR(IDAT)=WRS*DATA(IDAT,J)-WIS*DATA(IDAT,J+1)
            7         TEMPI(IDAT)=WRS*DATA(IDAT,J+1)+WIS*DATA(IDAT,J)
                    DO 8 IDAT=1,NDAT
                      DATA(IDAT,J)=DATA(IDAT,I)-TEMPR(IDAT)
            8         DATA(IDAT,J+1)=DATA(IDAT,I+1)-TEMPI(IDAT)
                    DO 9 IDAT=1,NDAT
                      DATA(IDAT,I)=DATA(IDAT,I)+TEMPR(IDAT)
            9         DATA(IDAT,I+1)=DATA(IDAT,I+1)+TEMPI(IDAT)
           12       CONTINUE
        ! Increment trigonometric table index
                  IW= IW+1
           13     CONTINUE
                MMAX=ISTEP
                GOTO 2
        ! Not yet done.
              END IF
        ! All done.
              RETURN
              END
         
        !========================================================================
              SUBROUTINE REALFV(DATA,ISIGN,NDAT0,NDAT,FFTAB,MAXTAB,N)
        !========================================================================
        ! REALFV is a "vectorizable" version of REALF2 (for multiple transforms)
        ! REALFV does the same thing as REALF2, but for NDAT arrays
        ! simultaneously. i.e. DATA is dimensioned as REAL (NDAT0,2*N),
        ! where NDAT0 is the declared inner dimesion, and NDAT is the operative
        ! inner dimension. FFTAB is a table of coefficients initialized by
        ! routine FOUR0, for complex transform length N. (See routine FOUR0
        ! for details)
        ! (Adapted from Numerical Recipes - W.H.Press et al, Cambridge(1986))
        !
              IMPLICIT REAL (A-H,O-Z)
              IMPLICIT INTEGER (I-N)
            !   INCLUDE 'ftable.inc'

              !integer MAXDAT, MAXTAB
              real, dimension(4,2,MAXTAB) :: FFTAB
        
              INTEGER N
              real, dimension(NDAT0,2*N) :: DATA
              REAL H1R(MAXDAT), H1I(MAXDAT), H2R(MAXDAT), H2I(MAXDAT)
              INTEGER IW, ISGN2
        ! Check if temporary storage is sufficient
              IF (NDAT.GT.MAXDAT) CALL SPERR('FOURN', &
               'NDAT too large; recompile with increased MAXDAT')
        ! Initialize the recurrence
              C1=0.5
              IF (ISIGN.EQ.1) THEN
                ISGN2= 1
                C2=-0.5
                CALL FOURV(DATA,+1,NDAT0,NDAT,FFTAB,MAXTAB,N)
        ! The forward transform is here
        !
              ELSE
        ! Otherwise set up for an inverse transform
        !
                ISGN2= 2
                C2=0.5
              END IF
        ! Initialize trigonometric table look-up index
              IW= 1
              N2P3=2*N+3
              DO 11 I=2,N/2+1
        ! Case I=1 done separately below
                I1=2*I-1
                I2=I1+1
                I3=N2P3-I2
                I4=I3+1
                WRS= FFTAB(3,ISGN2,IW)
                WIS= FFTAB(4,ISGN2,IW)
        ! The two separate transforms are separated out of Z.
        !
                DO 7 IDAT=1,NDAT
                  H1R(IDAT)= C1*(DATA(IDAT,I1)+DATA(IDAT,I3))
                  H1I(IDAT)= C1*(DATA(IDAT,I2)-DATA(IDAT,I4))
                  H2R(IDAT)=-C2*(DATA(IDAT,I2)+DATA(IDAT,I4))
            7     H2I(IDAT)= C2*(DATA(IDAT,I1)-DATA(IDAT,I3))
        !
        ! Here they are recombined to form the true transform of the original
        ! real data.
                DO 9 IDAT=1,NDAT
                  DATA(IDAT,I1)= H1R(IDAT)+WRS*H2R(IDAT)-WIS*H2I(IDAT)
                  DATA(IDAT,I2)= H1I(IDAT)+WRS*H2I(IDAT)+WIS*H2R(IDAT)
                  DATA(IDAT,I3)= H1R(IDAT)-WRS*H2R(IDAT)+WIS*H2I(IDAT)
            9     DATA(IDAT,I4)=-H1I(IDAT)+WRS*H2I(IDAT)+WIS*H2R(IDAT)
        !
        ! Increment the trigonometric table look-up index
                IW= IW+1
           11   CONTINUE
              IF (ISIGN.EQ.1) THEN
                DO 13 IDAT=1,NDAT
                  H1R(IDAT)=DATA(IDAT,1)
           13     H1I(IDAT)=DATA(IDAT,2)
                DO 14 IDAT=1,NDAT
           14     DATA(IDAT,1)=H1R(IDAT)+H1I(IDAT)
                DO 15 IDAT=1,NDAT
           15     DATA(IDAT,2)=H1R(IDAT)-H1I(IDAT)
        ! Squeeze the first and last data together to get them all within the
        ! original array.
        !
              ELSE
        ! This is the inverse transform for the case ISIGN=-1.
        !
                DO 17 IDAT=1,NDAT
                  H1R(IDAT)=DATA(IDAT,1)
           17     H1I(IDAT)=DATA(IDAT,2)
                DO 18 IDAT=1,NDAT
           18     DATA(IDAT,1)=C1*(H1R(IDAT)+H1I(IDAT))
                DO 19 IDAT=1,NDAT
           19     DATA(IDAT,2)=C1*(H1R(IDAT)-H1I(IDAT))
                CALL FOURV(DATA,-1,NDAT0,NDAT,FFTAB,MAXTAB,N)
              END IF
              RETURN
              END
         
        !========================================================================
              SUBROUTINE SPERR(PROC, MESG)
        !========================================================================
        ! SPERR aborts the program with a diagnostic error message.
        !   Input: CHARACTER*(*) PROC -- name of calling function/subroutine
        !         CHARACTER*(*) MESG -- diagnostic message
        !
              CHARACTER*(*) PROC, MESG
              WRITE(*,'(A,A,A,A)') '**Error in ', PROC, ': ', MESG
              STOP 'SPERR'
              END
         
        !
        !  All-latitudes-at-a-time spectral transforms begin here
        !
        !========================================================================
              SUBROUTINE PHCOPY(FPH2,FPH1)
        !========================================================================
        !  PHCOPY copies physical space representation FPH1 to FPH2
        !  Input:  REAL FPH1(K1MAX,K2MAX) -- physical representation oF F.
        !  Output: REAL FPH2(K1MAX,K2MAX) -- copy of FPH1.
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              REAL FPH2(K1MAX,K2MAX), FPH1(K1MAX,K2MAX)
              INTEGER J,K
              DO 10 K=1,K2
                DO 10 J=1,K1
                  FPH2(J,K)= FPH1(J,K)
           10     CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE MULCOS(FPH2,FPH1,C0,I0)
        !========================================================================
        !  MULCOS multiplies physical space representations by a coefficient
        !  C0 and by COS(PHI)**I0.
        !  Input:  REAL FPH1(K1MAX,K2MAX) -- physical representation of F
        !          REAL C0 -- coefficient to multiply F by
        !          INTEGER I0 -- exponent of COS(PHI)
        !  Output: REAL FPH2(K1MAX,K2MAX) -- C0*COS(PHI)**I0 times FPH1
        !  (Note: FPH2 may be the same array as FPH1)
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              REAL FPH2(K1MAX,K2MAX), FPH1(K1MAX,K2MAX), C0
              INTEGER I0
              REAL RTEMP(K2MAX)
              INTEGER J,K
              IF (I0.EQ.0) THEN
                DO 10 K=1,K2
           10     RTEMP(K)= C0
              ELSE IF (I0.EQ.1) THEN
                DO 20 K=1,K2
           20     RTEMP(K)= C0 * COSPHI(K)
              ELSE IF (I0.EQ.-1) THEN
                DO 30 K=1,K2
           30     RTEMP(K)= C0 * COSINV(K)
              ELSE
                DO 40 K=1,K2
           40     RTEMP(K)= C0 * (COSPHI(K)**I0)
              ENDIF
              DO 50 K=1,K2
                DO 50 J=1,K1
                  FPH2(J,K)= RTEMP(K) * FPH1(J,K)
           50     CONTINUE
              RETURN
              END
         
        !========================================================================
              SUBROUTINE PHYSIC(FPH,FSP)
        !========================================================================
        !  PHYSIC converts the spectral space representation FSP to physical
        !  space representation FPH.
        !  Note: Truncation order is assumed to be N1 (i.e. as for a scalar)
        !  Input:  COMPLEX FSP(0:M1MAX,0:N1MAX) -- spectral representation of F
        !  Output: REAL FPH(K1MAX,K2MAX) -- physical space representation of F
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              REAL FPH(K1MAX,K2MAX)
              COMPLEX FSP(0:M1MAX,0:N1MAX)
              EXTERNAL SP2LA
              CALL S2P_SP2LA( FPH, FSP, N1)
              RETURN
              END
         
        ! !========================================================================
        !       SUBROUTINE JACOB1( FPH, ASP, BLMPH, BMUPH )
        ! !========================================================================
        ! !  Given spectral space representation ASP of function A, and physical
        ! !  space representation BLMPH / BMUPH of the LAMBDA/MU derivatives of a
        ! !  function B, JACOB1 computes the jacobian of A & B, as defined by
        ! !  F = [1/(A0**2)]*[(dA/dLAMBDA)*(dB/dMU)-(dA/dMU)*(dB/dLAMBDA)]
        ! !  returning the jacobian FPH in physical space.
        ! !  Input:  COMPLEX ASP(0:M1MAX,0:N1MAX) -- spectral representation of A
        ! !          REAL BLMPH(K1MAX,K2MAX), BMUPH(K1MAX,K2MAX)-- physical space
        ! !               representation of dB/dLAMBDA, dB/dMU
        ! !  Output: REAL FPH(K1MAX,K2MAX) -- jacobian F in physical space
        ! !
        
        !       IMPLICIT NONE
        !     !   INCLUDE 'spcons.inc'
        !     !   INCLUDE 'spgrid.inc'
        
        !       REAL FPH(K1MAX,K2MAX)
        !       COMPLEX ASP(0:M1MAX,0:N1MAX)
        !       REAL BLMPH(K1MAX,K2MAX), BMUPH(K1MAX,K2MAX)
        !       REAL ALMPH(K1MAX,K2MAX), AMUPH(K1MAX,K2MAX)
        !       INTEGER J,K
        !       !EXTERNAL XSP2LA, YSP2LA
        ! !  Compute LAMBDA/MU derivatives of A
        !       CALL S2P( ALMPH, ASP, XSP2LA, N1)
        !       CALL S2P( AMUPH, ASP, YSP2LA, N1)
        ! !  Compute jacobian F in physical space
        !       DO 10 K=1,K2
        !         DO 10 J=1,K1
        !           FPH(J,K)= A0QINV * &
        !                     ( ALMPH(J,K)*BMUPH(J,K) - AMUPH(J,K)*BLMPH(J,K) )
        !    10     CONTINUE
        !       RETURN
        !       END
         
        ! !========================================================================
        !       SUBROUTINE JACOB2( FPH, ASP, BSP )
        ! !========================================================================
        ! !  Given spectral space representations ASP & BSP of functions A & B,
        ! !  JACOB2 computes the jacobian of A & B, as defined by
        ! !  F = [1/(A0**2)]*[(dA/dLAMBDA)*(dB/dMU)-(dA/dMU)*(dB/dLAMBDA)]
        ! !  returning the jacobian FPH in physical space.
        ! !  Input:  COMPLEX ASP(0:M1MAX,0:N1MAX), BSP(0:M1MAX,0:N1MAX) --
        ! !                     spectral representations of functions A, B
        ! !  Output: REAL FPH(K1MAX,K2MAX) -- jacobian F in physical space
        ! !
        
        !       IMPLICIT NONE
        !     !   INCLUDE 'spcons.inc'
        !     !   INCLUDE 'spgrid.inc'
        
        !       REAL FPH(K1MAX,K2MAX)
        !       COMPLEX ASP(0:M1MAX,0:N1MAX), BSP(0:M1MAX,0:N1MAX)
        !       REAL BLMPH(K1MAX,K2MAX), BMUPH(K1MAX,K2MAX)
        !       EXTERNAL XSP2LA, YSP2LA
        ! !  Compute LAMBDA/MU derivatives of B
        !       CALL S2P( BLMPH, BSP, XSP2LA, N1)
        !       CALL S2P( BMUPH, BSP, YSP2LA, N1)
        ! !  Call JACOB1 to do the rest of the job
        !       CALL JACOB1( FPH, ASP, BLMPH, BMUPH )
        !       RETURN
        !       END NB: commented out by Wim
         
        !========================================================================
              SUBROUTINE HVELOC(UPH,VPH,PSISP,CHISP)
        !========================================================================
        !  HVELOC calculates the horizontal velocity components U, V (UPH, VPH)
        !  in physical space, given the spectral space representation of
        !  streamfunction PSISP and velocity-potential CHISP.
        !  Input:  COMPLEX PSISP(0:M1MAX,0:N1MAX) -- spectral streamfunction
        !          COMPLEX CHISP(0:M1MAX,0:N1MAX) -- spectral velocity-potential
        !  Output: REAL UPH(K1MAX,K2MAX) -- zonal component of velocity
        !          REAL VPH(K1MAX,K2MAX) -- meridional component of velocity
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              REAL UPH(K1MAX,K2MAX), VPH(K1MAX,K2MAX)
              COMPLEX PSISP(0:M1MAX,0:N1MAX), CHISP(0:M1MAX,0:N1MAX)
              COMPLEX ULA(0:M1MAX,K2MAX), VLA(0:M1MAX,K2MAX)
              INTEGER K
        !  Calculate velocity components (ULA,VLA) in latitude space
              DO 10 K=1,K2
           10   CALL UVLAK(ULA(0,K),VLA(0,K),PSISP,CHISP,K)
        !  Convert to physical space
              CALL L2P(UPH,ULA)
              CALL L2P(VPH,VLA)
              RETURN
              END
         
        !========================================================================
              SUBROUTINE SPECTR(FSP,FPH)
        !========================================================================
        !  SPECTR converts the physical space representation FPH to spectral
        !  space, and accumulates in FSP.
        !  Note: Truncation order is assumed to be N1 (i.e. as for a scalar)
        !  Input:  COMPLEX FSP(0:M1MAX,0:N1MAX) -- old spectral coefficients
        !          REAL FPH(K1MAX,K2MAX) -- physical space representation of F
        !  Output: COMPLEX FSP(0:M1MAX,0:N1MAX) -- spectral representation of F
        !                                          plus old coefficients
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              COMPLEX FSP(0:M1MAX,0:N1MAX)
              REAL FPH(K1MAX,K2MAX)
              !EXTERNAL LA2SP
              CALL P2S_LA2SP( FSP, FPH, N1)
              RETURN
              END
         
        !========================================================================
              SUBROUTINE DIVERG(FSP,FPH1,FPH2)
        !========================================================================
        !  Given the component physical space representations FPH1 and FPH2
        !  of a vector field (say U,V), DIVERG computes the divergence in
        !  spectral space, and accumulates in FSP.
        !  Note: Truncation order is assumed to be N1 (i.e. as for a scalar)
        !  Input:  COMPLEX FSP(0:M1MAX,0:N1MAX) -- old spectral coefficients
        !          REAL FPH1(K1MAX,K2MAX), FPH2(K1MAX,K2MAX) -- physical
        !          representation of vector components (zonal, meridional)
        !  Output: COMPLEX FSP(0:M1MAX,0:N1MAX) -- spectral divergence of F
        !                                          plus old coefficients
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              COMPLEX FSP(0:M1MAX,0:N1MAX)
              REAL FPH1(K1MAX,K2MAX), FPH2(K1MAX,K2MAX)
              REAL TEMPH(K1MAX,K2MAX)
              EXTERNAL XLA2SP, YLA2SP
        !  Accumulate zonal-derivative of FPH1 / COSPHI
              CALL MULCOS( TEMPH, FPH1, A0INV, -1)
              CALL P2S_XLA2SP( FSP, TEMPH, N1 )
        !  Accumulate (1/A0) * MU-derivative of COSPHI * FPH2
              CALL MULCOS( TEMPH, FPH2, A0INV, +1)
              CALL P2S_YLA2SP( FSP, TEMPH, N1 )
              RETURN
              END
         
        !========================================================================
              SUBROUTINE CURLZ(FSP,FPH1,FPH2)
        !========================================================================
        !  Given the component physical space representations FPH1 and FPH2
        !  of a vector field (say U,V), CURLZ computes the curl (z-component) in
        !  spectral space, and accumulates in FSP.
        !  Note: Truncation order is assumed to be N1 (i.e. as for a scalar)
        !  Input:  COMPLEX FSP(0:M1MAX,0:N1MAX) -- old spectral coefficients
        !          REAL FPH1(K1MAX,K2MAX), FPH2(K1MAX,K2MAX) -- physical
        !          representation of vector components (zonal, meridional)
        !  Output: COMPLEX FSP(0:M1MAX,0:N1MAX) -- spectral curl of F
        !                                          plus old coefficients
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
        
              COMPLEX FSP(0:M1MAX,0:N1MAX)
              REAL FPH1(K1MAX,K2MAX), FPH2(K1MAX,K2MAX)
              REAL TEMPH(K1MAX,K2MAX)
              EXTERNAL XLA2SP, YLA2SP
        !  Accumulate zonal-derivative of FPH2 / COSPHI
              CALL MULCOS( TEMPH, FPH2, A0INV, -1)
              CALL P2S_XLA2SP( FSP, TEMPH, N1 )
        !  Accumulate (1/A0) * MU-derivative of -COSPHI * FPH1
              CALL MULCOS( TEMPH, FPH1, -A0INV, +1)
              CALL P2S_YLA2SP( FSP, TEMPH, N1 )
              RETURN
              END
         
        !========================================================================
              SUBROUTINE S2P(FPH,FSP,FSP2LA,N0)
        !========================================================================
        !  S2P converts the spectral space representation FSP to physical
        !  space representation FPH, for all latitudes. FSP2LA is the
        !  intermediate subroutine to be used to convert from spectral to
        !  latitude space. i.e. FSP2LA = SP2LA or XSP2LA or YSP2LA
        !  **NOTE** Remember to declare FSP2LA as external to the calling routine
        !
        !  Input:  COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F.
        !          SUBROUTINE FSP2LA(FLAK,FSP,P0,H0,N0,K) -- spectral to latitude
        !                                             space conversion subroutine
        !          INTEGER N0 -- order of truncation (N1 or N2).
        !  Output: REAL FPH(K1MAX,K2MAX) -- physical space representation of F
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
            !   INCLUDE 'sppoly.inc'
        
              INTEGER N0
              REAL FPH(K1MAX,K2MAX)
              COMPLEX FSP(0:M1MAX,0:N0)
              INTEGER K
              COMPLEX FLA(0:M1MAX,K2MAX)
        
        !  Dummy declaration (sometimes needed with the IMPLICIT NONE option)
        !       INTEGER FSP2LA
               EXTERNAL FSP2LA
        
        !  Convert to latitude space in temporary storage, using FSP2LA
              DO 10 K=1,K2
           10   CALL FSP2LA(FLA(0,K),FSP,P2(0,0,K),H2(0,0,K),N0,K)
        !  Convert to physical space
              CALL L2P(FPH,FLA)
              RETURN
              END

!========================================================================
              SUBROUTINE S2P_SP2LA(FPH,FSP,N0)
                !========================================================================
                !  S2P converts the spectral space representation FSP to physical
                !  space representation FPH, for all latitudes. FSP2LA is the
                !  intermediate subroutine to be used to convert from spectral to
                !  latitude space. i.e. FSP2LA = SP2LA or XSP2LA or YSP2LA
                !  **NOTE** Remember to declare FSP2LA as external to the calling routine
                !
                !  Input:  COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F.
                !          SUBROUTINE FSP2LA(FLAK,FSP,P0,H0,N0,K) -- spectral to latitude
                !                                             space conversion subroutine
                !          INTEGER N0 -- order of truncation (N1 or N2).
                !  Output: REAL FPH(K1MAX,K2MAX) -- physical space representation of F
                !
                
                      IMPLICIT NONE
                    !   INCLUDE 'spcons.inc'
                    !   INCLUDE 'spgrid.inc'
                    !   INCLUDE 'sppoly.inc'
                
                      INTEGER N0
                      REAL FPH(K1MAX,K2MAX)
                      COMPLEX FSP(0:M1MAX,0:N0)
                      INTEGER K
                      COMPLEX FLA(0:M1MAX,K2MAX)
                
                !  Dummy declaration (sometimes needed with the IMPLICIT NONE option)
                !       INTEGER FSP2LA
                !       EXTERNAL FSP2LA
                
                !  Convert to latitude space in temporary storage, using FSP2LA
                      DO 10 K=1,K2
                   10   CALL SP2LA(FLA(0,K),FSP,P2(0,0,K),H2(0,0,K),N0,K)
                !  Convert to physical space
                      CALL L2P(FPH,FLA)
                      RETURN
                      END
         
        !========================================================================
              SUBROUTINE P2S(FSP,FPH,FLA2SP,N0)
        !========================================================================
        !  P2S adds to spectral space representation FSP the contribution
        !  from physical space representation FPH, for all latitudes. FLA2SP is
        !  the intermediate subroutine to be used to convert from latitude to
        !  spectral space. i.e. FLA2SP = LA2SP or XLA2SP or YLA2SP
        !  **NOTE** Remember to declare FLA2SP as external to the calling routine
        !
        !  Input:  COMPLEX FSP(0:M1MAX,0:N0) -- old spectral coefficients
        !          REAL FPH(K1MAX,K2MAX) -- physical space representation of F
        !          SUBROUTINE FLA2SP(FSP,FLAK,P0,H0,N0,K) -- latitude to spectral
        !                                             space conversion subroutine
        !          INTEGER N0 -- order of truncation (N1 or N2)
        !  Output: COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F,
        !                                       plus the old coefficients
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
            !   INCLUDE 'sppoly.inc'
        
              INTEGER N0
              COMPLEX FSP(0:M1MAX,0:N0)
              REAL FPH(K1MAX,K2MAX)
              INTEGER K
              COMPLEX FLA(0:M1MAX,K2MAX)
        
        !  Dummy declaration (sometimes needed with the IMPLICIT NONE option)
        !       INTEGER FLA2SP
              EXTERNAL FLA2SP
        
        !  Convert to latitude space in temporary storage
              CALL P2L(FLA,FPH)
        !  Convert to spectral space using FLA2SP (accumulate)
              DO 10 K=1,K2
           10   CALL FLA2SP(FSP,FLA(0,K),P2(0,0,K),H2(0,0,K),N0,K)
              RETURN
              END

        !       CALL P2S( FSP, TEMPH, XLA2SP, N1 )
        ! !  Accumulate (1/A0) * MU-derivative of -COSPHI * FPH1
        !       CALL MULCOS( TEMPH, FPH1, -A0INV, +1)
        !       CALL P2S( FSP, TEMPH, YLA2SP, N1 )

    !========================================================================
    SUBROUTINE P2S_XLA2SP(FSP,FPH,N0)
    !========================================================================
    !  P2S adds to spectral space representation FSP the contribution
    !  from physical space representation FPH, for all latitudes. FLA2SP is
    !  the intermediate subroutine to be used to convert from latitude to
    !  spectral space. i.e. FLA2SP = LA2SP or XLA2SP or YLA2SP
    !  **NOTE** Remember to declare FLA2SP as external to the calling routine
    !
    !  Input:  COMPLEX FSP(0:M1MAX,0:N0) -- old spectral coefficients
    !          REAL FPH(K1MAX,K2MAX) -- physical space representation of F
    !          SUBROUTINE FLA2SP(FSP,FLAK,P0,H0,N0,K) -- latitude to spectral
    !                                             space conversion subroutine
    !          INTEGER N0 -- order of truncation (N1 or N2)
    !  Output: COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F,
    !                                       plus the old coefficients
    !
                
    IMPLICIT NONE
!   INCLUDE 'spcons.inc'
!   INCLUDE 'spgrid.inc'
!   INCLUDE 'sppoly.inc'
    
    INTEGER N0
    COMPLEX FSP(0:M1MAX,0:N0)
    REAL FPH(K1MAX,K2MAX)
    INTEGER K
    COMPLEX FLA(0:M1MAX,K2MAX)
    
    !  Dummy declaration (sometimes needed with the IMPLICIT NONE option)
    !       INTEGER FLA2SP
    !        EXTERNAL FLA2SP
    
!  Convert to latitude space in temporary storage
    CALL P2L(FLA,FPH)
!  Convert to spectral space using FLA2SP (accumulate)
    DO 10 K=1,K2
  10   CALL XLA2SP(FSP,FLA(0,K),P2(0,0,K),H2(0,0,K),N0,K)
    RETURN
    END

!========================================================================
    SUBROUTINE P2S_YLA2SP(FSP,FPH,N0)
        !========================================================================
        !  P2S adds to spectral space representation FSP the contribution
        !  from physical space representation FPH, for all latitudes. FLA2SP is
        !  the intermediate subroutine to be used to convert from latitude to
        !  spectral space. i.e. FLA2SP = LA2SP or XLA2SP or YLA2SP
        !  **NOTE** Remember to declare FLA2SP as external to the calling routine
        !
        !  Input:  COMPLEX FSP(0:M1MAX,0:N0) -- old spectral coefficients
        !          REAL FPH(K1MAX,K2MAX) -- physical space representation of F
        !          SUBROUTINE FLA2SP(FSP,FLAK,P0,H0,N0,K) -- latitude to spectral
        !                                             space conversion subroutine
        !          INTEGER N0 -- order of truncation (N1 or N2)
        !  Output: COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F,
        !                                       plus the old coefficients
        !
                    
        IMPLICIT NONE
    !   INCLUDE 'spcons.inc'
    !   INCLUDE 'spgrid.inc'
    !   INCLUDE 'sppoly.inc'
        
        INTEGER N0
        COMPLEX FSP(0:M1MAX,0:N0)
        REAL FPH(K1MAX,K2MAX)
        INTEGER K
        COMPLEX FLA(0:M1MAX,K2MAX)
        
        !  Dummy declaration (sometimes needed with the IMPLICIT NONE option)
        !       INTEGER FLA2SP
        !        EXTERNAL FLA2SP
        
    !  Convert to latitude space in temporary storage
        CALL P2L(FLA,FPH)
    !  Convert to spectral space using FLA2SP (accumulate)
        DO 10 K=1,K2
      10   CALL YLA2SP(FSP,FLA(0,K),P2(0,0,K),H2(0,0,K),N0,K)
        RETURN
        END

!========================================================================
        SUBROUTINE P2S_LA2SP(FSP,FPH,N0)
            !========================================================================
            !  P2S adds to spectral space representation FSP the contribution
            !  from physical space representation FPH, for all latitudes. FLA2SP is
            !  the intermediate subroutine to be used to convert from latitude to
            !  spectral space. i.e. FLA2SP = LA2SP or XLA2SP or YLA2SP
            !  **NOTE** Remember to declare FLA2SP as external to the calling routine
            !
            !  Input:  COMPLEX FSP(0:M1MAX,0:N0) -- old spectral coefficients
            !          REAL FPH(K1MAX,K2MAX) -- physical space representation of F
            !          SUBROUTINE FLA2SP(FSP,FLAK,P0,H0,N0,K) -- latitude to spectral
            !                                             space conversion subroutine
            !          INTEGER N0 -- order of truncation (N1 or N2)
            !  Output: COMPLEX FSP(0:M1MAX,0:N0) -- spectral representation of F,
            !                                       plus the old coefficients
            !
                        
            IMPLICIT NONE
        !   INCLUDE 'spcons.inc'
        !   INCLUDE 'spgrid.inc'
        !   INCLUDE 'sppoly.inc'
            
            INTEGER N0
            COMPLEX FSP(0:M1MAX,0:N0)
            REAL FPH(K1MAX,K2MAX)
            INTEGER K
            COMPLEX FLA(0:M1MAX,K2MAX)
            
            !  Dummy declaration (sometimes needed with the IMPLICIT NONE option)
            !       INTEGER FLA2SP
            !        EXTERNAL FLA2SP
            
        !  Convert to latitude space in temporary storage
            CALL P2L(FLA,FPH)
        !  Convert to spectral space using FLA2SP (accumulate)
            DO 10 K=1,K2
          10   CALL LA2SP(FSP,FLA(0,K),P2(0,0,K),H2(0,0,K),N0,K)
            RETURN
            END
             
         
        !========================================================================
              SUBROUTINE L2P(FPH,FLA)
        !========================================================================
        !  L2P converts a function F from latitude space representation FLA
        !  to physical space representation FPH, for all latitudes.
        !  Input:  COMPLEX FLA(0:M1MAX,K2MAX) -- latitude representation of F
        !  Output: REAL FPH(K1MAX,K2MAX) -- physical representation of F
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
            !   INCLUDE 'spfftb.inc'
        
              REAL FPH(K1MAX,K2MAX)
              COMPLEX FLA(0:M1MAX,K2MAX)
              COMPLEX FTMKC( (K1MAX+1)/2 )
              REAL FTMKR( K1MAX ), TEMPH( K2MAX, K1MAX )
        !  Note that TEMPH is in transpose form
              EQUIVALENCE (FTMKR, FTMKC)
              INTEGER J,K,M
              IF (K1.EQ.1) THEN
        !  Axially symmetric case; handle separately
                DO 10 K=1,K2
          10      FPH(1,K)= REAL( FLA(0,K) )
              ELSE
        !  Iterate for each latitude
                DO 50 K=1,K2
        !  Copy the first M1+1 fourier coefficients to temporary array
        !  (Take complex conjugate because the discrete fourier transform
        !  according to REALFV is defined the "wrong" way)
        !  (The factor of 2 comes in because the inverse fourier transform
        !  according to REALFV is divided by K1/2 and not by K1. Note that
        !  subroutine P2L divides by a factor K1)
        !
                  DO 20 M=0,M1
                    FTMKC(M+1)= CONJG( FLA(M,K) ) * 2.0
           20       CONTINUE
        !  Zero all extra coefficients, extending M1+1 values to K1/2 values.
        !  Since we are avoiding aliasing, it is safe to assume K1/2 >= M1+1
                  DO 30 M=M1+1,K1/2 - 1
                    FTMKC(M+1)= (0.,0.)
           30       CONTINUE
        !  Copy to temporary transposed physical representation
                  DO 40 J=1,K1
           40       TEMPH(K,J)= FTMKR(J)
           50       CONTINUE
        !  Since FLA(*,K) is the fourier transform of real function
        !  TEMPH(K,*), invert the transform.
                CALL REALFV(TEMPH,-1,K2MAX,K2,FFTAB,NFFTAB,K1/2)
        !  Transpose temporary matrix to obtain proper result
                DO 60 K=1,K2
                  DO 60 J=1,K1
           60     FPH(J,K)= TEMPH(K,J)
              ENDIF
              RETURN
              END
         
        !========================================================================
              SUBROUTINE P2L(FLA,FPH)
        !========================================================================
        !  P2L converts a function F from physical space representation FPH
        !  to latitude space representation FLA, for all latitudes.
        !  Input:  REAL FPH(K1MAX,K2MAX) -- physical representation of F
        !  Output: COMPLEX FLA(0:M1MAX,K2MAX) -- latitude representation of F
        !
        
              IMPLICIT NONE
            !   INCLUDE 'spcons.inc'
            !   INCLUDE 'spgrid.inc'
            !   INCLUDE 'spfftb.inc'
        
              REAL FPH(K1MAX,K2MAX)
              COMPLEX FLA(0:M1MAX,K2MAX)
              COMPLEX FTMKC( (K1MAX+1)/2 )
              REAL FTMKR(K1MAX), TEMPH( K2MAX, K1MAX)
        !  Note that TEMPH is in transpose form
              EQUIVALENCE (FTMKR, FTMKC)
              INTEGER J,K,M
              IF (K1.EQ.1) THEN
        !  Axially symmetric case; handle separately
                DO 10 K=1,K2
          10      FLA(0,K)= CMPLX( FPH(1,K), 0. )
              ELSE
        !  Copy physical values into temporary storage
                DO 20 K=1,K2
                  DO 20 J=1,K1
           20       TEMPH(K,J)= FPH(J,K)
        !  Fourier transforms of real functions TEMPH(K,*)
                CALL REALFV(TEMPH,+1,K2MAX,K2,FFTAB,NFFTAB,K1/2)
        !  Iterate for each latitude
                DO 50 K=1,K2
        !  Untranspose the result
                  DO 30 J=1,K1
           30       FTMKR(J)= TEMPH(K,J) / K1
        !  Retrieve the real first fourier coefficient
        !  (Divide by K1 because the discrete fourier transform according
        !  to REALFV does not incorporate this factor)
        !
                  FLA(0,K)= CMPLX( FTMKR(1), 0. )
        !
        !  Retrieve the next M1 complex fourier coefficients
        !  (Take complex conjugate because the discrete fourier transform
        !  according to REALFV is defined the "wrong" way)
                  DO 40 M=1,M1
           40       FLA(M,K)= CONJG( FTMKC(M+1) )
           50     CONTINUE
              ENDIF
              RETURN
              END
         
        !
        !  Diagnostic routines begin here
        !
        ! !========================================================================
        !       SUBROUTINE ZAVGE(FPROF,FPH) commented out by Wim
        ! !========================================================================
        ! !  ZAVGE computes the zonal average FPROF of a global field FPH.
        ! !  Input:  REAL FPH(K1MAX,K2MAX) -- values of a function on the globe
        ! !  Output: REAL FPROF(K2MAX) -- zonally averaged profile of the function
        ! !                            (FPROF(K) is the average at latitude MU(K) )
        ! !
        
        !       IMPLICIT NONE
        !     !   INCLUDE 'spcons.inc'
        !     !   INCLUDE 'spgrid.inc'
        
        !       REAL FPROF(K2MAX), FPH(K1MAX,K2MAX)
        !       INTEGER K
        !       REAL TOTAL
        !       DO 10 K=1,K2
        !         FPROF(K)= TOTAL(FPH(1,K),K1) / K1
        !    10   CONTINUE
        !       RETURN
        !       END
         
        ! !========================================================================
        !       SUBROUTINE EDDY(EPH,FPH) commented out by Wim
        ! !========================================================================
        ! !  EDDY returns the eddy component EPH of a global field FPH.
        ! !  i.e. It subtracts out the zonal mean part.
        ! !  Input:  REAL FPH(K1MAX,K2MAX) -- values of a function on the globe
        ! !  Output: REAL EPH(K1MAX,K2MAX) -- eddy part of the function FPH
        ! !
        
        !       IMPLICIT NONE
        !     !   INCLUDE 'spcons.inc'
        !     !   INCLUDE 'spgrid.inc'
        
        !       REAL EPH(K1MAX,K2MAX), FPH(K1MAX,K2MAX)
        !       INTEGER J,K
        !       REAL ZMEAN
        !       REAL TOTAL
        !       DO 10 K=1,K2
        !         ZMEAN= TOTAL(FPH(1,K),K1) / K1
        !         DO 10 J=1,K1
        !           EPH(J,K)= FPH(J,K) - ZMEAN
        !    10     CONTINUE
        !       RETURN
        !       END
         
        ! !========================================================================
        !       REAL FUNCTION GAVGE(FPH) commented out by Wim
        ! !========================================================================
        ! !  GAVGE returns the global average of a global field FPH.
        ! !  Input:  REAL FPH(K1MAX,K2MAX) -- values of a function on the globe
        ! !  Output: Global average of the function
        ! !
        
        !       IMPLICIT NONE
        !     !   INCLUDE 'spcons.inc'
        !     !   INCLUDE 'spgrid.inc'
        
        !       REAL FPH(K1MAX,K2MAX)
        !       REAL FPROF(K2MAX)
        !       REAL WTOTAL
        !       CALL ZAVGE(FPROF,FPH)
        !       GAVGE= WTOTAL(G,FPROF,K2)
        !       RETURN
        !       END
         
        !========================================================================
              REAL FUNCTION TOTAL(F,N)
        !========================================================================
        !  TOTAL computes the sum of elements of real array F, of length N
              INTEGER N
              REAL F(N)
              INTEGER J
              REAL SUM
              SUM= 0.0
              DO 10 J=1,N
           10   SUM= SUM + F(J)
              TOTAL=SUM
              RETURN
              END
         
        !========================================================================
              REAL FUNCTION WTOTAL(W,F,N)
        !========================================================================
        !  WTOTAL computes the weighted sum of elements of real array F, using
        !  the real array of weights W, both of length N.
              INTEGER N
              REAL W(N),F(N)
              INTEGER J
              REAL SUM
              SUM= 0.0
              DO 10 J=1,N
           10   SUM= SUM + W(J) * F(J)
              WTOTAL= SUM
              RETURN
              END
         
        ! !========================================================================
        !       REAL FUNCTION WAVGE(W,F,N) commented out by Wim
        ! !========================================================================
        ! !  WAVGE computes the weighted average of elements of real array F,
        ! !  using the real array of weights W, both of length N.
        !       INTEGER N
        !       REAL W(N),F(N)
        !       REAL WTOTAL, TOTAL
        !       WAVGE= WTOTAL(W,F,N) / TOTAL(W,N)
        !       RETURN
        !       END
         
end module splib_mod