C SPGRID.INC
C  This file is included in the routines of PROGNOS.FOR and SPLIB.FOR.
C  It should also be included in a main driver.
C  Must include SPCONS beforehand
C  Values set by call to SPLINI
C
C------------------------------------------------
C M1        = maximum value of zonal wavenumber M
C N1        = order of triangular truncation  (for scalars)
C K1        = no. of equally spaced longitudes   ( K1 >= 3*M1 + 1 )
C (For efficient F.F.T'S, K1 is always chosen to be a power of 2 )
C K2        = no. of gaussian latitudes          ( K2 >= (3*N1 + 1)/2 )
C (For "aesthetic" reasons, K2 is always chosen to be an odd number)
C N2        = N1 + 1  (order of truncation for vector components)
C
      INTEGER M1, N1, K1, K2, N2
      COMMON/SPGRD1/ M1, N1, K1, K2, N2
C
C------------------------------------------------
C DEGLON(J) = longitude values (in degrees, 0 ... 360)
C DEGLAT(K) = gaussian latitude values (in degrees, -90 ... +90)
C LAMBDA(J) = longitude values (in radians, 0 ... 2*Pi)
C PHI(K)    = gaussian latitude values (in radians, -Pi/2 ... +Pi/2)
C MU(K)     = gaussian "latitude" values ( SIN(PHI(K)), -1 ... +1 )
C COSPHI(K) = cosine of latitude ( COS(PHI(K)), 0 ... 1 )
C COSINV(K) = reciprocal of cosine of latitude ( 1 / COS(PHI(K)) )
C G(K)      = gaussian quadrature weights (with factor of 1/2 included)
C D(M,N)    = coefficients of recurrence relation for H[M,N]
C CIM(M) = CMPLX(0.0, M)  (coefficient for evaluating LAMBDA-derivative)
C
      REAL DEGLON, DEGLAT, LAMBDA, PHI, MU, COSPHI, COSINV, G, D
      COMPLEX CIM
      COMMON/SPGRD2/ DEGLON(K1MAX), DEGLAT(K2MAX),
     1 LAMBDA(K1MAX), PHI(K2MAX), MU(K2MAX), COSPHI(K2MAX),
     2 COSINV(K2MAX), G(K2MAX), D(0:M1MAX, 0:N2MAX), CIM(0:M1MAX)
C
C------------------------------------------------
C A0        = planetary radius (in metres)
C A0Q       = A0**2
C A0INV     = 1 / A0
C A0QINV    = 1 / A0Q
C
      REAL A0, A0Q, A0INV, A0QINV
      COMMON/SPGRD3/ A0, A0Q, A0INV, A0QINV
C
