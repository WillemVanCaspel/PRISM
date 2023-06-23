c===============================================
      SUBROUTINE FORCING(DPGQSP,PGQSP1,TIM)
c===============================================
C Bottom boundary, thermal, gravity wave and
c zonally symmetric momentum forcing
c...no vorticity forcing
c******************************
c...***NOTE*** TT_OFF has been changed to midpoint of ramp up
c******************************

	IMPLICIT NONE

	INCLUDE 'mcons.inc'
	INCLUDE 'spcons.inc'
	INCLUDE 'mgrid.inc'
	INCLUDE 'spgrid.inc'
	INCLUDE 'spfftb.inc'
	INCLUDE 'tmarch.inc'

	REAL TIM
	REAL PI
	PARAMETER(PI=3.1415926) 
	COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
	COMPLEX PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
	COMPLEX GEOPIN(0:M1MAX,0:N1MAX)
	COMPLEX GEOPSP(0:M1MAX,0:N1MAX)
	COMPLEX GEOPSPB(0:M1MAX,0:N1MAX)
	COMPLEX CTEMP
	COMPLEX HOLD(0:M1MAX,0:N1MAX)
	REAL GEOP(K1MAX, K2MAX)
	REAL GEOG(K1MAX, K2MAX)
	REAL OMEG, Y, FACT, FACT0
      INTEGER LEV, ILON, ILAT, N, M, I, J, K, L, X, YT
	INTEGER NSPEC

	REAL GEO(K2MAX)
	COMMON /GEOBOT/ GEO, GEOP

	COMPLEX FORCESP(0:M1MAX,0:N1MAX)
	COMPLEX FORCESPT(0:M1MAX,0:N1MAX)
	REAL FU_PH(K1MAX, K2MAX), FV_PH(K1MAX, K2MAX)
	REAL FT_PH(K1MAX, K2MAX)
	REAL RLAT, Z, R2, FTEMP, RT2, FTTEMP
	REAL FM0, FT0, FTT0(6)
	REAL TFAC

	COMPLEX FSP(0:M1MAX,0:N1MAX)
	REAL FPH(K1MAX, K2MAX)

	REAL DTOR /0.0174533/
        integer fcall /0/

C.....Forcing parameters
	REAL FREQ, ZWN, GAMP, WID, CEN, W_ON, W_OFF
	REAL VAMP, VAC, VAW, VLC, VLW, VWN, VON, VOFF
	REAL VZW, VOMEG
	REAL TAMP, TAC, TAW, TLC, TLW, TWN, TON, TOFF
	REAL TZW, TOMEG
	REAL TT_ZW(6), TT_OMEG(6), TT_ON(6), TT_OFF(6)
	REAL TB_ON(6), TB_OFF(6)

	INTEGER NLAT, NALT, NSF, NSB, IS
	REAL HIN(200,200), PIN(200,200), LAT(200), PRES(200)
	REAL HEAT(K2MAX,L1MAX,6)
	REAL COP(K2MAX,L1MAX,6), SIP(K2MAX,L1MAX,6)
	REAL COB(K2MAX,L1MAX), SIB(K2MAX,L1MAX), PIB(200,200)!DUMMIES

	REAL CLATW(K2MAX)
	INTEGER NLATW, NTIMW, ND, ND0
	REAL    XNALTW, XNLATW, XNTIMW
	INTEGER INDX2(K2MAX), INDX3(L1MAX)
	REAL DAY, DELT, WT0, WT1, CRANG
	REAL TRAD(K2MAX,L1MAX)
    
      INTEGER ND0S, NDS
      REAL CRANGS, DAYS, DELTS, WT0S, WT1S
      LOGICAL FIRST /.TRUE./

C.....SPECTRUM PARAMS
	INTEGER MAXC,MAXA,MAXK
	PARAMETER(MAXC=300,MAXA=8,MAXK=3)
	INTEGER NAZ, NKH, NC
	REAL B0(MAXC,MAXA,MAXK,K2MAX) !SOURCE SPECTRA(C,NAZ,NKH,LAT)
	REAL C0(MAXC)		!PHASE SPEED SAMPLES
	REAL C(MAXC)		!PHASE SPEED SAMPLES
	REAL KH(MAXK,K2MAX)	!HORIZONTAL WAVENUMBER(NKH,NLAT)
	REAL EPS(MAXK,K2MAX)		!INTERMITTENCY
	INTEGER IZ0(MAXK,K2MAX)	!SOURCE GRID LEVEL
	REAL TLEV(K2MAX)	!TURBOPAUSE HEIGHT
        REAL CW                 !1/e-width of GW spectrum
	
C.....GW STUFF
	REAL KH1(MAXA), KH2(MAXA)
	REAL UCOMP(L1MAX)
	REAL GW(L1MAX)
	INTEGER IAZ
	REAL DAZ
	REAL RHO(L1MAX)
	REAL UGW(2,L1MAX)
	REAL BF(L1MAX), bff(k2max,l1max)
	REAL ZG(L1MAX)

	REAL BFRC(K2MAX,L1MAX,6)
	REAL DY
	INTEGER K1HOLD
	COMPLEX FSP1(0:M1MAX,0:N1MAX), FSP2(0:M1MAX,0:N1MAX)
	REAL U(K1MAX, K2MAX, JX:JY, L1MAX)
	REAL TP(K1MAX, K2MAX, L1MAX)
	REAL GWFRC(K1MAX,K2MAX,JX:JY,L1MAX)
	INTEGER IC0, MC		!B0(IC0:ICM+MC-1)>0
	INTEGER K1GW
	REAL C1, C2, DZ, DTDZ, DZ2
	INTEGER LP, LM, LL
	INTEGER ITEMP
	REAL ZF(L1MAX)		!alt grid for Bfrc

	integer FLAG !0 - center spectrum at ucomp, otherwise at 0

	real hr, hr0
	REAL UBAR(K2MAX,L1MAX), VBAR(K2MAX,L1MAX), TBAR(K2MAX,L1MAX)

	REAL ZONW(K1MAX,K2MAX)

	character*80 arg
	INTEGER IARG
	common /opens/ iarg
	real lensq

C.....LUNAR TIDE STUFF
      REAL FDAY, NDAY !DAY OF FULL MOON
      REAL MPHAS0, MPHSP, MPHAS, LOC0, MFAC
      DOUBLE PRECISION MFACM2, MFACN2, MFACS2, MFACL2
      REAL SIGLR, SIGSR, SIGZE, SIGMU
      DOUBLE PRECISION LOCM2, LOCN2, LOCS2, LOCL2
      REAL AMPM2, AMPN2, AMPS2, AMPL2
      COMPLEX LGEOPSP(0:M1MAX,0:N1MAX)
      DOUBLE PRECISION SPHASE, PPHASE, LTIM, SPIND, CONFAC, VPHASE, TPHASE
      DOUBLE PRECISION TZERO

c.....PYTHON SH ASSIMILATION DECLARATIONS
      INTEGER TDIM, MAXM, MAXN, MAXL, COUNTER
      INTEGER BRECL, RECLEN, MAXMAS, MAXNAS
      INTEGER MAXMS, MAXNS, TDIMS, BRECLS, NIN
      INTEGER MAXMSURF, MAXNSURF, TDIMSURF, BRECLSURF
      
      REAL CTIMW(20000)
      REAL CTIMWS(20000)
      
      DOUBLE PRECISION SLON(10000), SLAT(10000)
      DOUBLE PRECISION ZOUTU(100000), ZOUTV(100000), ZOUTT(100000)
      DOUBLE PRECISION ZDU(100,100), ZDV(100,100), ZDT(100,100)
      
      REAL INTERPU(K1MAX, K2MAX), INTERPV(K1MAX, K2MAX) 
      REAL INTERPT(K1MAX, K2MAX)

      COMPLEX TEMT0(0:M1MAX,0:N1MAX,L1MAX)
      COMPLEX TEMT1(0:M1MAX,0:N1MAX,L1MAX)
      COMPLEX VORT0(0:M1MAX,0:N1MAX,L1MAX)
      COMPLEX VORT1(0:M1MAX,0:N1MAX,L1MAX)
      
      DOUBLE PRECISION LOND(200),LATD(200)   ! assimilated winds and temps
      DOUBLE PRECISION LONDS(200),LATDS(200) ! ocean waves
      DOUBLE PRECISION LONDG(200),LATDG(200) ! surface geopotential
      
      COMPLEX SURF0(0:M1MAX,0:N1MAX)
      COMPLEX SURF1(0:M1MAX,0:N1MAX)
      COMPLEX SURF2(0:M1MAX,0:N1MAX)
      COMMON /SURFACE/ SURF1
      
      REAL READU0(200,200,200), READV0(200,200,200)
      REAL READT0(200,200,200)
      REAL READU1(200,200,200), READV1(200,200,200)
      REAL READT1(200,200,200)
      
      DOUBLE PRECISION SURFP0(200,200), SURFP1(200,200), SURFP2(200,200)

	  CHARACTER*80 FOLDER,FILENAME

      REAL TEMP0(K2MAX,L1MAX), U0(K2MAX,L1MAX)
      COMPLEX VOR0(0:M1MAX,0:N1MAX,L1MAX)
      COMPLEX VORTEMP(0:M1MAX,0:N1MAX)
      COMPLEX TEMPTEMP(0:M1MAX,0:N1MAX)
      COMMON /ZONE/ TEMP0, U0, VOR0
      
C     Input and output folders and files
      CHARACTER*50 INFOLD, INFILE, OUTFOLD, OUTFILE, TIDEFILE, TTIDEFILE
      CHARACTER*50 INITTEMP, INITZONAL, INITGRID, INITZGRID, GWFDUMP, CURT
      CHARACTER*50 OCEANWAVE
      COMMON /FILES/ INFOLD, INFILE, OCEANWAVE, OUTFOLD, OUTFILE,
     1   INITTEMP, INITZONAL, INITGRID, GWFDUMP, TIDEFILE, TTIDEFILE,CURT
c===============================================

	hr=tim/3600.

	IF (FCALL .EQ. 0) THEN
	  hr0=tim/3600.
	  FCALL=1

c...thermal forcing structures from file
c...multiple zonal wavenumbers and frequencies
	  WRITE(6,*) 'ENTER # HEAT SOURCE FILES:'
	  READ(5,*) NSF
	  write(6,*) NSF

	  open(unit=23,file=TRIM(OUTFOLD)//TRIM(GWFDUMP))
	  DO IS=1, NSF
	    WRITE(6,*) 'Enter zonal waven, frequency (1/day):'
	    READ (5,*) TT_ZW(IS), TT_OMEG(IS)
	    write(6,*) TT_ZW(IS), TT_OMEG(IS)
	    WRITE(6,*) 'Enter time on/off for forcing:'
	    READ (5,*) TT_ON(IS), TT_OFF(IS)
	    write(6,*) TT_ON(IS), TT_OFF(IS)

	    TT_OMEG(IS) = 2. * PI * (TT_OMEG(IS) / 86400.)
	    TT_ON(IS)=TT_ON(IS)*86400.
	    TT_OFF(IS)=TT_OFF(IS)*86400.

c	    CALL GETARG(IARG,ARG) !CONTINUE RUN FROM THIS FILE
c	    OPEN(UNIT=33,STATUS='OLD',FILE=ARG(1:(INDEX(ARG,' ')-1)),shared)
c	    iarg=iarg+1
	    write(6,'(A)') arg

	    READ(33,*) NLAT,NALT
	    READ(33,*) (LAT(I), I=1,NLAT)
	    READ(33,*) (PRES(I), I=1,NALT)
	    READ(33,*) ((HIN(I,J), I=1,NLAT),J=1,NALT)
	    READ(33,*) ((PIN(I,J), I=1,NLAT),J=1,NALT)
	    CLOSE (UNIT=33)

	    DO I=1,NALT
	      PRES(I)=PRES(I)/1.E3
	    ENDDO

	    CALL FINTERP(PRES, NALT, LAT, NLAT, HIN, PIN,
	1	 HEAT(1,1,IS), SIP(1,1,IS), COP(1,1,IS))

	    write(23,*) ((hin(i,j),i=1,nlat),j=1,nalt)
	    write(23,*) ((pin(i,j),i=1,nlat),j=1,nalt)
	    write(23,*) ((heat(i,j,is),i=1,k2),j=1,l1)
	    write(23,*) ((cop(i,j,is),i=1,k2),j=1,l1)
	  ENDDO
	  close(23)

c================================================================
c.....climatology
c================================================================
c	  CALL GETARG(IARG,ARG)
c.....READ IN FILE DIMENSIONS
      OPEN(55, FILE=TRIM(INFOLD)//TRIM(INFILE), STATUS='OLD',
	1      FORM='UNFORMATTED', ACCESS='DIRECT', RECL=1000)
c	  write(6,'(A)') arg
c	  IARG=IARG+1 ! from 5 to 6 at this point

      READ (55,REC=1) RECLEN, TDIM, MAXM, MAXN, MAXL
      
	  NTIMW=TDIM
	  CLOSE(55)

      BRECL=MAXM*MAXN*MAXL*4 + 8
      MAXMAS=MIN(MAXM - 1,M1)
      MAXNAS=MIN(MAXN - 1,N1)

c.....READ IN ALTITUDE AND TIME ARRAYS
      OPEN(55, FILE=TRIM(INFOLD)//TRIM(INFILE), STATUS='OLD',
	1      FORM='UNFORMATTED', ACCESS='DIRECT', RECL=BRECL)

      READ (55,REC=2) RECLEN, (CTIMW(I),I=1,TDIM)
      READ (55,REC=3) RECLEN, (LOND(I),I=1,MAXM)
      READ (55,REC=4) RECLEN, (LATD(I),I=1,MAXN)
                       
	  CRANG=CTIMW(NTIMW)-CTIMW(1)
      NIN = K1*K2
      
c...zonal mean body force from file (GW surrogate)
	  WRITE(6,*) 'ENTER # MOM SOURCE FILES:'
	  READ(5,*) NSB
	  write(6,*) NSB

	  DO IS=1, NSB
	    WRITE(6,*) 'Enter time on/off for forcing:'
	    READ(5,*) TB_ON(IS), TB_OFF(IS)
	    write(6,*) TB_ON(IS), TB_OFF(IS)

	    TB_ON(IS)=TB_ON(IS)*86400.
	    TB_OFF(IS)=TB_OFF(IS)*86400.

c	    CALL GETARG(IARG,ARG) !CONTINUE RUN FROM THIS FILE
c	    OPEN(UNIT=33,STATUS='OLD',FILE=ARG(1:(INDEX(ARG,' ')-1)),shared)
c	    iarg=iarg+1
c	    write(6,'(A)') arg
	    READ(33,*) NLAT,NALT
	    READ(33,*) (LAT(I), I=1,NLAT)
	    READ(33,*) (PRES(I), I=1,NALT)
	    READ(33,*) ((HIN(I,J), I=1,NLAT),J=1,NALT)
	    CLOSE (UNIT=33)

	    DO I=1,NALT
	      PRES(I)=PRES(I)/1.E3
	    ENDDO

        CALL FINTERP(PRES, NALT, LAT, NLAT, HIN, PIB,
     1    BFRC(1,1,IS), SIB(1,1), COB(1,1))
     
	  DO I=1,K2
	    DO J=1,L1
	      IF ((DEGLAT(I) .LT. -80.) .AND.
     1      (-7*ALOG(PLV(J)) .LT. 30.))
     1      BFRC(I,J,is)=0.
	    ENDDO
	  ENDDO

c	    write(23,*) ((bfrc(i,j,is),i=1,k2),j=1,l1)
c	    close(23)
	  ENDDO

c.....bottom boundary forcing
c.....GAMP, W_ON, W_OFF used to control realistic bottom topography
c	  OPEN(UNIT=44, FILE='GEOPSP.DAT',STATUS='OLD',SHARED,
c	1      FORM='UNFORMATTED')
c	  READ(44) NSPEC
c	  READ(44) ((geopin(M,N),M=0,NSPEC),N=0,NSPEC)
c	  CLOSE(44)

C------------------------------------------------------------------
C.....Bottom boundary forcing, use tanh to ramp up.

C.....Read in surface pressure file    
      OPEN(46, FILE=TRIM(INFOLD)//TRIM('surface_geopotential.dat'), STATUS='OLD',
	1      FORM='UNFORMATTED', ACCESS='DIRECT', RECL=1000)
    
      READ (46,REC=1) RECLEN, MAXMSURF, MAXNSURF
      CLOSE(46)
      
      BRECLSURF= MAXMSURF*MAXNSURF*8 + 8
      
      OPEN(46, FILE=TRIM(INFOLD)//TRIM('surface_geopotential.dat'), STATUS='OLD',
	1      FORM='UNFORMATTED', ACCESS='DIRECT', RECL=BRECLSURF)
    
      READ (46,REC=2)  RECLEN, (LONDG(M), M=1, MAXMSURF)
      READ (46,REC=3)  RECLEN, (LATDG(M), M=1, MAXNSURF)
    
C.....read single surface level field (1000 hpa)
      READ (46,REC=4)  RECLEN,
     1     ((SURFP2(M,N), M=1, MAXMSURF), N=1, MAXNSURF)
    
      NIN = K1*K2
     
      COUNTER = 1 ! sigprim lat and lon interpolation grid
      DO X=1,K1
        DO YT=1,K2
          SLON(COUNTER) = DBLE(DEGLON(X))
          SLAT(COUNTER) = DBLE(DEGLAT(YT))
          COUNTER = COUNTER + 1
        ENDDO
      ENDDO
 
      CALL pwl_interp_2d (MAXMSURF,MAXNSURF,LONDG,LATDG,SURFP2(1:MAXMSURF,1:MAXNSURF),
     1                    NIN,SLON,SLAT,ZOUTU)
                                
      COUNTER = 1
      DO X=1,K1
        DO YT=1,K2
          INTERPU(X,YT) = REAL(ZOUTU(COUNTER)) * 1
          COUNTER = COUNTER + 1
        ENDDO
      ENDDO
                  
      CALL ZEROSP(SURF2,N1)
      CALL SPECTR(SURF2,INTERPU)
      CLOSE(46)
      
C.....Read in ocean surface file    
      OPEN(45, FILE=TRIM(INFOLD)//TRIM(OCEANWAVE), STATUS='OLD',
	1      FORM='UNFORMATTED', ACCESS='DIRECT', RECL=1000)
    
      READ (45,REC=1) RECLEN, TDIMS, MAXMS, MAXNS
      CLOSE(45)
      
      BRECLS= MAXMS*MAXNS*8 + 8
      
      OPEN(45, FILE=TRIM(INFOLD)//TRIM(OCEANWAVE), STATUS='OLD',
	1      FORM='UNFORMATTED', ACCESS='DIRECT', RECL=BRECLS)
    
      READ (45,REC=2) RECLEN, (CTIMWS(I),I=1,TDIMS) ! DOY
      READ (45,REC=3) RECLEN, (LONDS(I),I=1,MAXMS) 
      READ (45,REC=4) RECLEN, (LATDS(I),I=1,MAXNS) 
      
C     Time shift ocean tide file 
      DO I=1,TDIMS
        CTIMWS(I) = CTIMWS(I) + 6.0 
      ENDDO
                      
      CRANGS=CTIMWS(TDIMS)-CTIMWS(1)
      ND0S=0
                     
	  WRITE(6,*)
	1      'ENTER FREQUENCY IN 1./DAYS (POSITIVE IS EASTWARD),'
	  WRITE(6,*)  'ZONAL WAVENUMBER S, GEOP AMPLITUDE,'
	  WRITE(6,*)
	1      'WIDTH AND CENTER OF FORCING REGION (DEGREES LAT):'
	  READ(5,*)  FREQ, ZWN, GAMP, WID, CEN ! line 12 in solar_tide.inp
	  write(6,*)  FREQ, ZWN, GAMP, WID, CEN

	  WRITE(6,*) 'ENTER TIME ON & OFF (DAYS) FOR WAVE FORCING:'
	  READ(5,*)  W_ON, W_OFF ! line 13 in solar_tide.inp
	  write(6,*)  W_ON, W_OFF

	  W_ON=W_ON*86400.
	  W_OFF=W_OFF*86400.

c.....gravity wave setup
C.....SET K1GW=-1,0 FOR NO GW FORCING
	  WRITE(6,*) 'enter longitude sampling, LENSQ'
	  READ(5,*) K1GW, LENSQ ! line 14 in solar_tide.inp
	  write(6,*) K1GW, LENSQ

	  IF (K1GW .GT. 0) THEN

	  WRITE(6,*) 'enter FLAG:'
	  READ(5,*) FLAG
	  write(6,*) FLAG

	  NC = 121 !# of wave packet
          NAZ = 4  !# of azimuth
	  NKH = 2  !# of GW zonal wavenumber
          CW  = 47 !1/e-width of GW spectrum, unit: m/s
	  DO I=1,NKH
             DO J=1,K2
               !ZG index is reversed in gwf.f
               IZ0(I,J) = 131 - 117    !GW source level at 261 hPa = ~10 km alt
               EPS(I,J) = 0.0005  !intermittency
             END DO
	  END DO 
          DO J=1,K2
	    TLEV(J) = 115. ! make a turbopause = 115 km
            !invariant wavenumber = 2pi/500000m and 2pi/1000000m
	    KH(1,J) = 2. * PI / 500000.
            KH(2,J) = 2. * PI / 1000000.
	  END DO
          DO I=1,NC
	    C0(I) = REAL(I)-1. ! define GW wave packet phase speed
          END DO
          DO I=1,NC
            DO J=1,NAZ
	      DO K=1,NKH
                DO L=1,K2
		  B0(I,J,K,L) = 0.004 * EXP(-(C0(I)/CW)**2) !#define GW source amp.
                END DO
	      END DO
            END DO
	  END DO
   
          !#DEFINE azimuth direction for gw propagation
          DO IAZ=1,NAZ
	    DAZ=360./NAZ*(IAZ-1)
            KH1(IAZ)=COSD(DAZ)
	    KH2(IAZ)=SIND(DAZ)
          END DO
C.....GW PARAMS FROM FILE
c	  CALL GETARG(IARG,ARG)	!CONTINUE RUN FROM THIS FILE
c	  OPEN(UNIT=33,STATUS='OLD',FILE=ARG(1:(INDEX(ARG,' ')-1)), SHARED)
c	  IARG=IARG+1
c	  write(6,'(A)') arg
C	    OPEN (UNIT=33,FILE=' ',STATUS='OLD', SHARED)
C	    READ(33,*) NAZ, NKH, NC
C	    READ(33,*) (TLEV(J),J=1,K2)
C	    READ(33,*) ((IZ0(I,J),I=1,NKH),J=1,K2)
C	    READ(33,*) ((EPS(I,J),I=1,NKH),J=1,K2)
C	    READ(33,*) ((KH(I,J),I=1,NKH),J=1,K2)
C	    READ(33,*) (C0(I),I=1,NC)
C	    READ(33,*) ((((B0(I,J,K,L),
C     1                 I=1,NC),J=1,NAZ),K=1,NKH),L=1,K2)
C	    CLOSE (UNIT=33)

C	    open(unit=71,file='gwdump.dat')
C
C	    DO IAZ=1,NAZ
C	      DAZ=360./NAZ*(IAZ-1)
C	      KH1(IAZ)=COSD(DAZ)
C	      KH2(IAZ)=SIND(DAZ)
C	    ENDDO
      ENDIF

c	  CMAX=150.

	  DO L=1,L1
	    RHO(L)=PLV(L1-L+1)*1.34
	    ZG(L)=-7.*ALOG(PLV(L1-L+1)*PSURF/1.E5)
	  ENDDO

	  C1=289./7.E6
	  C2=2./49.

C     Lunar tide stuff
      WRITE(6,*) 'Enter day of full moon' ! line 15 in solar_tide.inp
      READ(5,*) FDAY, NDAY, SPIND, MFAC, AMPM2, AMPN2, AMPS2, AMPL2
      WRITE(6,*) 'Gravitational tide parameters: ', FDAY, NDAY, SPIND,MFAC, AMPM2, AMPN2, AMPS2, AMPL2

	  MPHAS0=2*PI*FDAY/29.53059 ! lunar phase_0
            
      CONFAC = SQRT(24./5.)/2. ! see https://shtools.github.io/SHTOOLS/complex-spherical-harmonics.html
      
C     factor of sqrt(24/5.) is to convert from normalized P2 to 3*cos(lat)**2
      MFACM2 = -0.7933*(1 + 0.302 - 0.609)*CONFAC*AMPM2
      MFACN2 = -0.1518*(1 + 0.302 - 0.609)*CONFAC*AMPN2
      MFACS2 = -0.3700*(1 + 0.302 - 0.609)*CONFAC*AMPS2
      MFACL2 = -0.1005*(1 + 0.302 - 0.609)*CONFAC*AMPL2
      
      SIGLR = 2*PI/23.934469/3600. - 2*PI/27.321661/86400.
      SIGSR = 2*PI/23.934469/3600. - 2*PI/365.256363/86400.
      SIGZE = 2*PI/23.934469/3600.
      SIGMU = 2*PI/27.554551/86400.
      
      ENDIF ! fcall now set to 1
      
C-----FORCING ROUTINE BEGINS HERE------------------------------------------------------------------------------------------------------------

C     interpolate to current time and cycle periodically over data time range
      DAYS=TIM/86400.
      DELTS=MOD(DAYS-REAL(CTIMWS(1)),CRANGS)
      DO WHILE (DELTS .LT. 0)
        DELTS=DELTS+CRANGS
      ENDDO
      DAYS=CTIMWS(1)+DELTS
        
      NDS=1
      DO J=2,TDIMS-1
        IF (DAYS .GT. REAL(CTIMWS(J))) NDS=J
      ENDDO    
      
C     Linearly interpolate in time
      DELTS=CTIMWS(NDS+1)-CTIMWS(NDS)
      WT0S=(CTIMWS(NDS+1)-DAYS)/DELTS
      WT1S=(DAYS-CTIMWS(NDS))/DELTS
            
      IF (NDS .NE. ND0S) THEN
 
        NIN = K1*K2
        READ (45,REC=4+NDS) RECLEN,
     1    ((SURFP0(I,J), I=1, MAXMS), J=1, MAXNS)
     
        COUNTER = 1 ! sigprim lat and lon interpolation grid
        DO X=1,K1
          DO YT=1,K2
            SLON(COUNTER) = DBLE(DEGLON(X))
            SLAT(COUNTER) = DBLE(DEGLAT(YT))
            COUNTER = COUNTER + 1
          ENDDO
        ENDDO
 
        CALL pwl_interp_2d (MAXMS,MAXNS,LONDS,LATDS,SURFP0(1:MAXMS,1:MAXNS),
     1                      NIN,SLON,SLAT,ZOUTU)
                                
        COUNTER = 1
        DO X=1,K1
          DO YT=1,K2
            INTERPU(X,YT) = REAL(ZOUTU(COUNTER))
            COUNTER = COUNTER + 1
          ENDDO
        ENDDO
                
        CALL ZEROSP(SURF0,N1)
        CALL SPECTR(SURF0,INTERPU)
                           
        READ (45,REC=5+NDS) RECLEN,
     1    ((SURFP1(I,J), I=1, MAXMS), J=1, MAXNS)
                
        COUNTER = 1 ! sigprim lat and lon interpolation grid
        DO X=1,K1
          DO YT=1,K2
            SLON(COUNTER) = DBLE(DEGLON(X))
            SLAT(COUNTER) = DBLE(DEGLAT(YT))
            COUNTER = COUNTER + 1
          ENDDO
        ENDDO
                            
        CALL pwl_interp_2d (MAXMS,MAXNS,LONDS,LATDS,SURFP1(1:MAXMS,1:MAXNS),
     1                      NIN,SLON,SLAT,ZOUTU)
                        
        COUNTER = 1
        DO X=1,K1
          DO YT=1,K2
            INTERPU(X,YT) = REAL(ZOUTU(COUNTER))
            COUNTER = COUNTER + 1
          ENDDO
        ENDDO
                
        CALL ZEROSP(SURF1,N1)
        CALL SPECTR(SURF1,INTERPU)
        
        WRITE(6,*) 'TTIM:', NDS, DAYS, CTIMWS(NDS)
        ND0S=NDS
      ENDIF
            
c.....tranfer into array of correct size for physic, delsq
      DO N=0,N1
        DO M=0,M1
	      GEOPSP(M,N)=(WT0S*SURF0(M,N) + WT1S*SURF1(M,N))*MFAC + SURF2(M,N)
        ENDDO
      ENDDO

      DO ILAT=1,K2
        DO ILON=1,K1
          GEOP(ILON,ILAT)=GEO(ILAT)
        ENDDO
      ENDDO

      CALL PHYSIC(GEOG,GEOPSP)
      CALL SPECTR(GEOPSPB,GEOP)
      
      CALL DELSQ(GEOPSP,GEOPSP,N1)
      CALL DELSQ(GEOPSPB,GEOPSPB,N1)    
      
c.....boundary forcing
      TFAC=(TANH((TIM-W_ON)/(24*60*60*3)))
      FACT0=GAMP*TFAC
          
C.....PUT IN COMMON BLOCK FOR OUTPUT
      DO ILAT=1,K2
        DO ILON=1,K1
          GEOP(ILON,ILAT) = GEOG(ILON,ILAT) * FACT0 + GEO(ILAT) * 0
        ENDDO
      ENDDO
      
      DO N=0,N1
        DO M=0,MIN0(N,M1)
          CTEMP = FACT0 * GEOPSP(M,N) + GEOPSPB(M,N) * 0
          DO LEV=1,L1
            DPGQSP(M,N,JDIV,LEV)=DPGQSP(M,N,JDIV,LEV)-CTEMP
          END DO
        END DO
      END DO
      
c.....gravitational potential forcing, see doi:10.1029/2012JA017962

      TZERO = (NDAY + TIM/86400.)/36525. ! chapman & lindzen p. 78
            
      SPHASE = (270.43659 + 481267.89057*TZERO + 0.00198*TZERO**2)
      PPHASE = (334.32956 + 4069.03403*TZERO - 0.01032*TZERO**2 
     1         - 0.00001*TZERO**3)
      VPHASE = (-9.26009 + 445267.12165*TZERO + 0.00168*TZERO**2)
            
      SPHASE = SPHASE/360.*2.*PI
      PPHASE = PPHASE/360.*2.*PI
      VPHASE = VPHASE/360.*2.*PI
      TPHASE = TIM/86400.*2.*PI
                 
      MPHAS=TIM*MPHSP - MPHAS0
      
      LOC0=4*PI*TIM/86400. - 2.*MPHAS - 162.2*2*PI/360.
      
      LTIM = TIM - SPIND*86400. ! 80 day spin-up time used in input files
          
      LOCM2 = 2*SIGLR*LTIM - 2*MPHAS0
      LOCN2 = (2*SIGLR - SIGMU)*LTIM - 2*MPHAS0 - SPHASE + PPHASE
      LOCS2 = 2*(SIGSR*LTIM)
      LOCL2 = 2*(SIGZE*LTIM)
      
      LOCM2 = 2*(TPHASE - VPHASE)
      LOCN2 = 2*(TPHASE - VPHASE) - SPHASE + PPHASE
            
      CALL ZEROSP(LGEOPSP,N1)
      LGEOPSP(2,2)=(MFACM2*CMPLX(COS(LOCM2),SIN(LOCM2)) +
     1              MFACN2*CMPLX(COS(LOCN2),SIN(LOCN2)) +
     1              MFACS2*CMPLX(COS(LOCS2),SIN(LOCS2)) +
     1              MFACL2*CMPLX(COS(LOCL2),SIN(LOCL2)))
      CALL DELSQ(LGEOPSP,LGEOPSP,N1)

      DO LEV=1,L1
        DPGQSP(2,2,JDIV,LEV)=DPGQSP(2,2,JDIV,LEV) - LGEOPSP(2,2) ! minus for forcing
      END DO

c-----------------------------------------------------------------
c.....detailed thermal forcing
C      DO IS=1,NSF
C        FTT0(IS)=.25 * (1.+TANH((TIM-TT_ON(IS))/172800.))
C     1               * (1.-TANH((TIM-TT_OFF(IS))/172800.))
c	  FTT0(IS)=.5 * (1.+TANH((TIM-TT_ON(IS))/tt_off(is)))
C      ENDDO
C
c...this was originally in forcing_lunar.f. Above is taken from forcing_cycle.f
c	DO IS=1,NSF
cc	  FTT0(IS)=.25 * (1.+TANH((TIM-TT_ON(IS))/172800.))
cc	1              * (1.-TANH((TIM-TT_OFF(IS))/172800.))
c	  FTT0(IS)=.5 * (1.+TANH((TIM-TT_ON(IS))/tt_off(is)))
c	ENDDO

C      DO LEV=1,L1
C        DO K=1,K2
C          DO J=1,K1
C            FT_PH(J,K) = 0.
C            
C            DO IS=1,NSF
C              FTEMP=HEAT(K,LEV,IS) *
C     1        (COS(TT_ZW(IS)*LAMBDA(J)-TT_OMEG(IS)*TIM)*
C     1         COP(K,LEV,IS)
C     1        -SIN(TT_ZW(IS)*LAMBDA(J)-TT_OMEG(IS)*TIM)*
C     2         SIP(K,LEV,IS))
C
C               FT_PH(J,K) =FT_PH(J,K)+ FTEMP*FTT0(IS)
C            ENDDO
C            
C          ENDDO
C        ENDDO
C
C	  if (lev .eq. L1-4)  write(23,*) (ft_ph(j,30),j=1,k1)
C        CALL ZEROSP(FORCESPT,N1)
C        CALL SPECTR(FORCESPT,FT_PH)
C
C        DO N=0,N1
C          DO M=0,MIN0(N,M1)
C            DPGQSP(M,N,JPOT,LEV)=DPGQSP(M,N,JPOT,LEV)+FORCESPT(M,N)
C            WRITE(6,*) FORCESPT(M,N)
C          ENDDO
C        ENDDO
C      ENDDO

C-------------------------------------------------------------------
C.....zonal momentum forcing
C      DO IS=1,NSB
C        FTT0(IS)=.25 * (1.+TANH((TIM-TB_ON(IS))/172800.))
C     1               * (1.-TANH((TIM-TB_OFF(IS))/172800.))
C      ENDDO
C
C      DO LEV=1,L1
C        DO K=1,K2
C          DO J=1,K1
C            FU_PH(J, K)=0.0
C            FV_PH(J, K)=0.0
C
C            DO IS=1,NSB
C              FU_PH(J,K) =FU_PH(J,K)+ BFRC(K,LEV,IS)*FTT0(IS)
C            ENDDO
C          ENDDO
C        ENDDO
C
C        CALL ZEROSP(FORCESP,N1)
C        CALL CURLZ(FORCESP,FU_PH,FV_PH)
C
C        DO N=0,N1
C          DO M=0,MIN0(N,M1)
C            DPGQSP(M,N,JVOR,LEV)=DPGQSP(M,N,JVOR,LEV)+FORCESP(M,N)
C          ENDDO
C        ENDDO
C      ENDDO

C=====================================================================
C.....CLIMATOLOGY
C.....INTERPOLATE TO CURRENT TIME
C.....ASSUME PERIODIC OVER TIME INTERVAL CRANG

      DAY=TIM/86400.
      DELT=MOD(DAY-CTIMW(1),CRANG)
      DO WHILE (DELT .LT. 0)
        DELT=DELT+CRANG
      ENDDO
      DAY=CTIMW(1)+DELT

      ND=1
      DO J=2,NTIMW-1
        IF (DAY .GT. CTIMW(J)) ND=J
      ENDDO

C     Linear interpolation coefficients between T0 and T1
      DELT=CTIMW(ND+1)-CTIMW(ND)
      WT0=(CTIMW(ND+1)-DAY)/DELT
      WT1=(DAY-CTIMW(ND))/DELT

      if (ND0 .ne. ND) then
        WRITE(6,*) ND
        
        READ (55,REC=2+3*ND) RECLEN,
     1     (((READU0(M,N,L), M=1, MAXM), N=1, MAXN),L=1,MAXL)
        READ (55,REC=3+3*ND) RECLEN,
     1     (((READV0(M,N,L), M=1, MAXM), N=1, MAXN),L=1,MAXL)
        READ (55,REC=4+3*ND) RECLEN,
     1     (((READT0(M,N,L), M=1, MAXM), N=1, MAXN),L=1,MAXL)
                          
        DO L=1,MAXL ! interpolate 2d
   
          COUNTER = 1 ! sigprim lat and lon interpolation grid
          DO X=1,K1
            DO YT=1,K2
              SLON(COUNTER) = DBLE(DEGLON(X))
              SLAT(COUNTER) = DBLE(DEGLAT(YT))
              COUNTER = COUNTER + 1
            ENDDO
          ENDDO
          
          DO X=1,MAXM
            DO YT=1,MAXN
              ZDU(X,YT) = DBLE(READU0(X,YT,L))
              ZDV(X,YT) = DBLE(READV0(X,YT,L))
              ZDT(X,YT) = DBLE(READT0(X,YT,L))
            ENDDO
          ENDDO
                      
          CALL pwl_interp_2d (MAXM,MAXN,LOND,LATD,ZDU(1:MAXM,1:MAXN),NIN,SLON,
     1                        SLAT,ZOUTU)
          CALL pwl_interp_2d (MAXM,MAXN,LOND,LATD,ZDV(1:MAXM,1:MAXN),NIN,SLON,
     1                        SLAT,ZOUTV)
          CALL pwl_interp_2d (MAXM,MAXN,LOND,LATD,ZDT(1:MAXM,1:MAXN),NIN,SLON,
     1                        SLAT,ZOUTT)
          
          COUNTER = 1
          DO X=1,K1
            DO YT=1,K2
              INTERPU(X,YT) = REAL(ZOUTU(COUNTER))
              INTERPV(X,YT) = REAL(ZOUTV(COUNTER))
              INTERPT(X,YT) = REAL(ZOUTT(COUNTER))
              COUNTER = COUNTER + 1
            ENDDO
          ENDDO
          
C         IF (FIRST) THEN
C           WRITE(6,*) READU0(2,1:MAXN,1)
C           write(6,*) 'break'
C           WRITE(6,*) ZDU(2,1:MAXN)
C           write(6,*) 'break'
C           WRITE(6,*) INTERPU(2,1:K2)
C           WRITE(6,*) 'BREAK'
C           WRITE(6,*) K1, K2, NIN
C           FIRST = .FALSE.
C         ENDIF
                
          CALL ZEROSP(TEMPTEMP,N1)
          CALL SPECTR(TEMPTEMP,INTERPT)
                    
          CALL ZEROSP(VORTEMP,N1)
          CALL CURLZ(VORTEMP,INTERPU,INTERPV)
          
          DO M=0,M1
            DO N=0,N1
              VORT0(M,N,L) = VORTEMP(M,N)
              TEMT0(M,N,L) = TEMPTEMP(M,N)
            ENDDO
          ENDDO
          
        ENDDO
     
        READ (55,REC=5+3*ND) RECLEN,
     1     (((READU1(M,N,L), M=1, MAXM), N=1, MAXN),L=1,MAXL)
        READ (55,REC=6+3*ND) RECLEN,
     1     (((READV1(M,N,L), M=1, MAXM), N=1, MAXN),L=1,MAXL)
        READ (55,REC=7+3*ND) RECLEN,
     1     (((READT1(M,N,L), M=1, MAXM), N=1, MAXN),L=1,MAXL)
     
        COUNTER = 1 ! sigprim lat and lon interpolation grid
        DO X=1,K1
          DO YT=1,K2
            SLON(COUNTER) = DBLE(DEGLON(X))
            SLAT(COUNTER) = DBLE(DEGLAT(YT))
            COUNTER = COUNTER + 1
          ENDDO
        ENDDO
        
        DO L=1,MAXL ! interpolate 2d
          DO X=1,MAXM
            DO YT=1,MAXN
              ZDU(X,YT) = DBLE(READU1(X,YT,L))
              ZDV(X,YT) = DBLE(READV1(X,YT,L))
              ZDT(X,YT) = DBLE(READT1(X,YT,L))
            ENDDO
          ENDDO
          
          CALL pwl_interp_2d (MAXM,MAXN,LOND,LATD,ZDU(1:MAXM,1:MAXN),
     1                        NIN,SLON,SLAT,ZOUTU)
          CALL pwl_interp_2d (MAXM,MAXN,LOND,LATD,ZDV(1:MAXM,1:MAXN),
     1                        NIN,SLON,SLAT,ZOUTV)
          CALL pwl_interp_2d (MAXM,MAXN,LOND,LATD,ZDT(1:MAXM,1:MAXN),
     1                        NIN,SLON,SLAT,ZOUTT)
         
          COUNTER = 1
          DO X=1,K1
            DO YT=1,K2
              INTERPU(X,YT) = REAL(ZOUTU(COUNTER))
              INTERPV(X,YT) = REAL(ZOUTV(COUNTER))
              INTERPT(X,YT) = REAL(ZOUTT(COUNTER))
              COUNTER = COUNTER + 1
            ENDDO
          ENDDO
                     
          CALL ZEROSP(VORTEMP,N1)
          CALL CURLZ(VORTEMP,INTERPU,INTERPV)
          
          CALL ZEROSP(TEMPTEMP,N1)
          CALL SPECTR(TEMPTEMP,INTERPT)
                
          DO M=0,M1
            DO N=0,N1
              VORT1(M,N,L) = VORTEMP(M,N) 
              TEMT1(M,N,L) = TEMPTEMP(M,N) 
            ENDDO
          ENDDO
        ENDDO
        
        WRITE(6,*) 'BTIM:', ND, DAY, CTIMW(ND)
        ND0=ND
      ENDIF

      DO M=0,M1
        DO N=0,N1
          DO L=1,L1
            TMNLV(M,N,L)=(WT0*TEMT0(M,N,L)
     1       +    WT1*TEMT1(M,N,L))
            VOR0(M,N,L)=(WT0*VORT0(M,N,L)
     1       +    WT1*VORT1(M,N,L))
          ENDDO
        ENDDO
      ENDDO
      
c======================================================================
c.....k1gw flag:
c.....-1 => No adjustment, No GW
c.....0  => dry adjustment, No GW
c.....K1 => adjustment + GW, maintain longitude grid
c.....any positive => adjustment + GW on alternate lon grid
C.....finer grid may be required for GW interaction with PW
c.....at coarse resolution.  K1GW=1 gives GW on zonal mean only.

c.....toggle k1 temporarily to K1GW
      K1HOLD=K1

      IF ((K1 .NE. K1GW) .AND. (K1GW .GT. 0)) THEN
        K1=K1GW
        IF (K1.GT.1) CALL FOUR0( FFTAB, NFFTAB, K1/2 )
      ENDIF

      IF (K1GW .GT. -1) THEN

C.......expand horizontal velocity vector field and temperature
        DO L=1,L1
	      CALL IDELSQ( FSP1, PGQSP1(0,0,JVOR,L), N1)
	      CALL IDELSQ( FSP2, PGQSP1(0,0,JDIV,L), N1)
	      CALL HVELOC( U(1,1,JX,L), U(1,1,JY,L), FSP1, FSP2 )
	      CALL PHYSIC( TP(1,1,L), PGQSP1(0,0,JPOT,L) )
        ENDDO

C.......DRY CONVECTIVE ADJUSTMENT:
        CALL DRY_ADJUST(DPGQSP,PGQSP1,U,TP,LENSQ)

      ENDIF

c-------------------------------------------------------------------
C.....Gravity wave forcing
C.....SKIP IF GW FORCING TURNED OFF

      IF (K1GW .LT. 1) RETURN
      IF (HR .GE. HR0+24.) HR0=HR

C...compute GW forcing field
c	DO L=1,L1
c	  CALL ZAVGE(UBAR(1,L),U(1,1,JX,L))
c	  CALL ZAVGE(VBAR(1,L),U(1,1,JY,L))
c	  CALL ZAVGE(TBAR(1,L),TP(1,1,L))
c	ENDDO

	DO N=1,K2
	  DO M=1,K1
	    DO L=1,L1
	      LL=L1+1-L
	      UGW(1,L)=U(M,N,JX,LL)
	      UGW(2,L)=U(M,N,JY,LL)
cc	      UGW(1,L)=UBAR(N,LL)
cc	      UGW(2,L)=VBAR(N,LL)
	      LM=MAX(L-1, 1)
	      LP=MIN(L+1,L1)
	      DZ=ZG(LP)-ZG(LM)
	      LP=MAX(LL-1, 1)
	      LM=MIN(LL+1,L1)
	      DTDZ=(TP(M,N,LP)-TP(M,N,LM))/DZ
	      BF(L)= MAX(C1*(DTDZ + C2*TP(M,N,LL)),1.E-6)
cc	      DTDZ=(TBAR(N,LP)-TBAR(N,LM))/DZ
cc	      BF(L)= MAX(C1*(DTDZ + C2*TBAR(N,LL)),1.E-6)
	      BF(L)=SQRT(BF(L))
	      bff(n,l)=bf(l)

	      GWFRC(M,N,JX,L)=0.
	      GWFRC(M,N,JY,L)=0.
	    ENDDO

C...LOOP OVER AZIMUTH------------------
	    DO IAZ=1,NAZ
	      DO L=1,L1
		UCOMP(L)=UGW(1,L)*KH1(IAZ) + UGW(2,L)*KH2(IAZ)
	      ENDDO

	      DO J=1,NKH

C...FIND WHERE B0 NE 0
		IC0=NC
		MC=0
		DO L=NC,1,-1
		  IF (B0(L,IAZ,J,N) .GT. 0) THEN
		    IC0=L
		    MC=MC+1
		  ENDIF
		ENDDO

c...if flag eq 0 then parameterize spectrum with c = c0-ucomp(iz0)
		DO L=1,NC
		  C(L)=C0(L)
		  IF (FLAG .EQ. 0) C(L)=C0(L)+UCOMP(IZ0(J,N))
		ENDDO

		CALL GWF
	1	     (C(IC0), B0(IC0,IAZ,J,N), EPS(J,N), KH(J,N),
	1	  ZG, RHO, UCOMP, BF, L1, MC, IZ0(J,N), GW, TLEV(N))

		DO L=1,L1
cc		DO L=12,L1  !for gwup
cc		DO L=20,L1  !for gwup50 version of heating_diur
		  LL=L1+1-L
		  GWFRC(M,N,JX,LL)=GWFRC(M,N,JX,LL)+GW(L)*KH1(IAZ)
		  GWFRC(M,N,JY,LL)=GWFRC(M,N,JY,LL)+GW(L)*KH2(IAZ)
		ENDDO
	      ENDDO
	    ENDDO
	  ENDDO
	ENDDO

	DO L=1,L1
	  CALL CURLZ
	1      (DPGQSP(0,0,JVOR,L),GWFRC(1,1,JX,L),GWFRC(1,1,JY,L))
	  CALL DIVERG
	2      (DPGQSP(0,0,JDIV,L),GWFRC(1,1,JX,L),GWFRC(1,1,JY,L))

cc	  CALL ZEROSP(HOLD,N1)
cc	  CALL CURLZ
cc	1      (HOLD,GWFRC(1,1,JX,L),GWFRC(1,1,JY,L))
cc	  DO N=0,N1
cc	    DPGQSP(0,N,JVOR,L)=DPGQSP(0,N,JVOR,L)+HOLD(0,N)
cc	  END DO

cc	  CALL ZEROSP(HOLD,N1)
cc	  CALL DIVERG
cc	2      (HOLD,GWFRC(1,1,JX,L),GWFRC(1,1,JY,L))
cc	  DO N=0,N1
cc	    DPGQSP(0,N,JDIV,L)=DPGQSP(0,N,JDIV,L)+HOLD(0,N)
cc	  END DO
	ENDDO

	IF (K1 .NE. K1HOLD) THEN
	  K1=K1HOLD
	  CALL FOUR0( FFTAB, NFFTAB, K1/2 )
	ENDIF

C	IF (HR .EQ. HR0) THEN
C	  WRITE(71,*) ((GWFRC(1,K,JX,L),K=1,K2),L=1,L1)
C	  WRITE(6,*) 'GWF',HR
C	ENDIF

	RETURN
	END

c===============================================
      SUBROUTINE FINTERP (PRES, NALT, LAT, NLAT, HIN, PIN,
     1  HEAT, SIP, COP)
c===============================================

	IMPLICIT NONE

	INCLUDE 'mcons.inc'
	INCLUDE 'spcons.inc'
	INCLUDE 'mgrid.inc'
	INCLUDE 'spgrid.inc'

C...    INPUT
        INTEGER NLAT, NALT
        REAL HIN(200,200), PIN(200,200), LAT(200), PRES(200)

C...    OUTPUT
	REAL HEAT(K2MAX,L1MAX)
	REAL COP(K2MAX,L1MAX), SIP(K2MAX,L1MAX)

	INTEGER I,J,K,L,KP1,KP2,LP1,LP2
	INTEGER INDX1(L1MAX+1), INDX2(K2MAX)
	REAL  WT1(L1MAX+1), WT2(K2MAX+1)
	REAL Y1, Y2, Y3, Y4, W1, W2, W3, W4

c-----------------------------------------------------------------

	LP1=1
	LP2=L1
	DO I=1,L1
	  IF (PLV(I) .LT. PRES(NALT)) LP1=I+1
	  IF (PLV(L1-I+1) .GT. PRES(1)) LP2=L1-I
	ENDDO

	KP1=1
	KP2=K2
	DO I=1,K2
	  IF (DEGLAT(I) .LT. LAT(1)) KP1=I+1
	  IF (DEGLAT(K2-I+1) .GT. LAT(NLAT)) KP2=K2-I
	ENDDO

	DO I=LP1,LP2
	  INDX1(I)=1
	  DO J=1,NALT
	    IF (PRES(J) .GT. PLV(I)) INDX1(I)=J
	  ENDDO
	ENDDO

	DO I=KP1,KP2
	  INDX2(I)=1
	  DO J=1,NLAT
	    IF (LAT(J) .LT. DEGLAT(I)) INDX2(I)=J
	  ENDDO
	ENDDO

	DO K=LP1,LP2
	  WT1(K)=(ALOG(PLV(K))-ALOG(PRES(INDX1(K)))) /
	1      (ALOG(PRES(INDX1(K)+1))-ALOG(PRES(INDX1(K))))
	ENDDO

	DO K=KP1,KP2
	  WT2(K)=(DEGLAT(K)-LAT(INDX2(K))) /
	1      (LAT(INDX2(K)+1)-LAT(INDX2(K)))
	ENDDO

	DO L=LP1,LP2
	  DO K=KP1,KP2
	    W1=(1.-WT2(K)) * (1.-WT1(L))
	    W2=(WT2(K))    * (1.-WT1(L))
	    W3=(1.-WT2(K)) * (WT1(L))
	    W4=(WT2(K))    * (WT1(L))

	    Y1=HIN(INDX2(K),  INDX1(L))
	    Y2=HIN(INDX2(K)+1,INDX1(L))
	    Y3=HIN(INDX2(K),  INDX1(L)+1)
	    Y4=HIN(INDX2(K)+1,INDX1(L)+1)
	    HEAT(K,L)=(W1*Y1+W2*Y2+W3*Y3+W4*Y4)

	    Y1=COS(PIN(INDX2(K),  INDX1(L)))
	    Y2=COS(PIN(INDX2(K)+1,INDX1(L)))
	    Y3=COS(PIN(INDX2(K),  INDX1(L)+1))
	    Y4=COS(PIN(INDX2(K)+1,INDX1(L)+1))
	    COP(K,L)=(W1*Y1+W2*Y2+W3*Y3+W4*Y4)

	    Y1=SIN(PIN(INDX2(K),  INDX1(L)))
	    Y2=SIN(PIN(INDX2(K)+1,INDX1(L)))
	    Y3=SIN(PIN(INDX2(K),  INDX1(L)+1))
	    Y4=SIN(PIN(INDX2(K)+1,INDX1(L)+1))
	    SIP(K,L)=(W1*Y1+W2*Y2+W3*Y3+W4*Y4)
	  ENDDO
	ENDDO

	END
