module IO_mod
    
    use Declarations_mod
    use Config_mod          , only: GRID 

    implicit none
     
    contains 

     !=============================================
              SUBROUTINE INIPRE_READ(PTHICK, PP, PHALFL, PBOT, PLID, NLEV)
        !=============================================
        !   read model pressure levels from file
        
              IMPLICIT NONE
       
        !       ! INCLUDE 'mcons.inc'
              REAL*4 PTHICK(L1MAX), PHALFL(0:L1MAX), PP(L1MAX+1) ! 4 byte now because of old input...
              REAL PBOT, PLID
              INTEGER NLEV, L 
              
              INTEGER RECLEN, BRECL
        
              INTEGER IARG
              COMMON /OPENS/ IARG
        
        !   Input and output folders and files
              CHARACTER*50 INFOLD, INFILE, OUTFOLD, OUTFILE, TIDEFILE, TTIDEFILE
              CHARACTER*50 INITTEMP, INITZONAL, INITGRID, INITZGRID, GWFDUMP, CURT
              CHARACTER*50 SURFFILE
              COMMON /FILES/ INFOLD, INFILE, SURFFILE, OUTFOLD, OUTFILE, &
                 INITTEMP, INITZONAL, INITGRID, GWFDUMP, TIDEFILE, TTIDEFILE,CURT

              write(*,*) 'durp durp darp ', GRID
                    
              OPEN(UNIT=18,STATUS='OLD',FILE=TRIM(INFOLD)//TRIM(INITGRID),SHARED, &
                  FORM='UNFORMATTED',ACCESS='DIRECT',RECL=100)
        
              READ (18,REC=1) RECLEN, NLEV
              CLOSE(18)
                    
              BRECL=NLEV*4 + 8
                    
              OPEN(UNIT=18,STATUS='OLD',FILE=TRIM(INFOLD)//TRIM(INITGRID),SHARED, &
                  FORM='UNFORMATTED',ACCESS='DIRECT',RECL=BRECL)
        
              READ (18,REC=2) RECLEN, (PP(L),L=1,NLEV)
              CLOSE(18)
              
              WRITE(6,*) 'Number of vertical sigma-levels: ', NLEV
        
              PHALFL(0) = PLID
              PHALFL(NLEV) = PBOT
        
              DO L=1, NLEV-1
                PHALFL(L) = .5*(PP(L)+PP(L+1))
                PTHICK(L) = PHALFL(L)-PHALFL(L-1)
              ENDDO
        
              PTHICK(NLEV) = PHALFL(NLEV)-PHALFL(NLEV-1)
                    
              RETURN
              END
        
        ! !============================================
        !       REAL FUNCTION HCOOL(PRESSURE)
        ! !============================================
        ! !   from jung-hee, based on a hoskins paper; newtonian cooling due to O3 and CO2
        ! 	REAL COOL(81) /
        ! 	6    0.05,
        ! 	1    0.1293,  0.1411,  0.1535,  0.1663,  0.1791,  0.1918,  0.2040,  0.2153,
        ! 	1    0.2255,  0.2342,  0.2412,  0.2463,  0.2492,  0.2500,  0.2497,  0.2494,
        ! 	2    0.2491,  0.2487,  0.2483,  0.2479,  0.2474,  0.2468,  0.2462,  0.2456,
        ! 	2    0.2448,  0.2440,  0.2431,  0.2421,  0.2410,  0.2398,  0.2384,  0.2369,
        ! 	3    0.2353,  0.2335,  0.2316,  0.2294,  0.2271,  0.2246,  0.2218,  0.2188,
        ! 	3    0.2156,  0.2121,  0.2083,  0.2043,  0.2000,  0.1954,  0.1906,  0.1854,
        ! 	4    0.1800,  0.1744,  0.1685,  0.1624,  0.1561,  0.1497,  0.1431,  0.1364,
        ! 	4    0.1297,  0.1229,  0.1162,  0.1095,  0.1029,  0.0825,  0.0614,  0.0502,
        ! 	5    0.0429,  0.0400,  0.0400,  0.0400,  0.0400,  0.0400,  0.0400,  0.0400,
        ! 	5    0.0400,  0.0400,  0.0400,  0.0400,  0.0420,  0.0530,  0.0713,  0.1106/
        
        ! 	REAL CPRES(81) /
        ! 	1    0.0927,
        ! 	1    0.2002,      0.2228,      0.2480,      0.2760,      0.3072,
        ! 	1    0.3420,      0.3807,      0.4237,      0.4716,      0.5250,
        ! 	2    0.5844,      0.6504,      0.7240,      0.8059,      0.8970,
        ! 	2    0.9985,      1.1114,      1.2371,      1.3770,      1.5327,
        ! 	3    1.7061,      1.8990,      2.1138,      2.3528,      2.6189,
        ! 	3    2.9151,      3.2448,      3.6117,      4.0202,      4.4749,
        ! 	4    4.9809,      5.5442,      6.1713,      6.8692,      7.6460,
        ! 	4    8.5108,      9.4733,     10.5447,     11.7372,     13.0646,
        ! 	5    14.5421,     16.1867,     18.0173,     20.0549,     22.3230,
        ! 	5    24.8476,     27.6577,     30.7856,     34.2672,     38.1426,
        ! 	6    42.4563,     47.2578,     52.6024,     58.5513,     65.1731,
        ! 	6    72.5437,     80.7479,     89.8800,    100.0450,    111.3590,
        ! 	7    123.9530,    137.9710,    153.5750,    170.9430,    190.2760,
        ! 	7    211.7950,    235.7470,    262.4090,    292.0860,    325.1180,
        ! 	8    361.8870,    402.8140,    448.3700,    499.0770,    555.5200,
        ! 	8    618.3450,    688.2760,    766.1150,    852.7579,    949.1990/
        
        !       IF (PRESSURE/1.E2 .LT. CPRES(1)) THEN
        !         HCOOL = .05
        !       ELSE IF (PRESSURE/1.E2 .GT. CPRES(81)) THEN
        !         HCOOL = .11
        !       ELSE
        !         DO I=1,80
        !           IF ((PRESSURE/1.E2 .GT. CPRES(I)) .AND.
        !      1      (PRESSURE/1.E2 .LE. CPRES(I+1)))
        !      1       HCOOL=COOL(I) + (PRESSURE/1.E2-CPRES(I))*(COOL(I+1)-COOL(I))/
        !      2      (CPRES(I+1)-CPRES(I))
        !         ENDDO
        !       ENDIF
        
        !       RETURN
        !       END
        
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