C SPCONS.INC
c
C  This file is included in the routines of PROGNOS.FOR and SPLIB.FOR.
C  It should also be included in a main driver.
C
C Parameter definition module (array dimensions)
C N1MAX: This is the maximum order of triangular truncation (i.e.
C        order of Legendre polynomials). 
C
C M1MAX: This is the maximum value of zonal wavenumber resolved by the
C        truncation. 
C
C K1MAX: Number of equally spaced longitudes in the transform grid
C        [should be the lowest integral power of 2 >= (3*M1MAX+1) ].
C
C K2MAX: Number of gaussian latitudes in the transform grid [should
C        be the lowest odd integer >= (3*N1MAX+1)/2. 
C (K2MAX is preferably chosen to be odd, so that the equator can be
C one of the gaussian latitudes)
C
C MLOW:  Cut-off value of zonal wavenumber m that denotes "low order
C        zonal truncation". For truncations with very few zonal
C        wavenumbers (M1 <= MLOW), transforms can be performed more
C        efficiently by reversing the order of certain nested
C        do-loops at the heart of the Legendre transform. But the
C        choice of MLOW would be very much machine dependent. The
C        default value is 4.
C
C [Other typical choices of (N1MAX, M1MAX, K1MAX, K2MAX) would be
C     (10, 10,  32, 17) or (21, 21,  64,  33)
C  or (42, 42, 128, 65) or (85, 85, 256, 129) ]
C

	INTEGER N1MAX, M1MAX, K1MAX, K2MAX, MLOW, N2MAX
	PARAMETER (N1MAX=42, M1MAX=42, K1MAX=128, K2MAX=65, MLOW=4)
	PARAMETER (N2MAX=N1MAX+1)

