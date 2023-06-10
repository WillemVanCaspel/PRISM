	SUBROUTINE zhu_HEAT (DPGQSP,PGQSP1, tim)

c...this version reads rad_gasses file from command line
	IMPLICIT NONE

	INCLUDE 'mcons.inc'
	INCLUDE 'spcons.inc'
	INCLUDE 'mgrid.inc'
	INCLUDE 'spgrid.inc'

	REAL TIM
	COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
	COMPLEX PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX)

	INTEGER KM
	PARAMETER (KM=76)

C...ATMOSPHERIC GASSES
	DOUBLE PRECISION CO2(KM,K2MAX), H2O(KM,K2MAX),
	1    O(KM,K2MAX), O2(KM,K2MAX), O3(KM,K2MAX)
	DOUBLE PRECISION ZM(KM), TEM(KM), BO3(KM,K2MAX),
	2    BO2(KM,K2MAX),QIR(KM),HSL(KM)
	REAL QNET(KM)
	DOUBLE PRECISION H2OLIM(KM),H2OTOP(KM)
	REAL CO2F(72,36), H2OF(72,36), OF(72,36),
	1    O2F(72,36), O3F(72,36)
	REAL CPRESW(172), CLATW(181)
	INTEGER NALTW, NLATW

	DOUBLE PRECISION F107	! =(70,145,220): 60 < F107 < 270
				!==> 10.7 cm solar flux index
	DOUBLE PRECISION theta27 ! = phase of 27-day periodic
				!variation: sin(theta27)
	DOUBLE PRECISION OME0 	!EFFECTIVE ALBEDO
	DOUBLE PRECISION DEC	!SOLAR DECLINATION
	INTEGER MDLRD	! 2 => diurnally averaged heating rate;
				!3 => local heating rate;
	INTEGER KLBD		! effective cloud-top level
	DOUBLE PRECISION XMU/0.5D0/	! COS(solar zenith angle) =
				!DUMMY VARIABLE
	INTEGER LTMU/1/		! 1 => input=DEC, PHI, TIME(hr);
				!0 => direct input = XMU

	COMMON /RADPARAM/ F107, THETA27, OME0, DEC, MDLRD, KLBD

	DOUBLE PRECISION DLAT, HOUR/0./	!USE HOUR WHEN LOCAL HEATING
				! RATE USED
	REAL Y2A(72,36), Y2(72), VAL
	INTEGER IER, WLUN

	REAL TEMP(KM), NUM(KM), TP(K1MAX,K2MAX,L1MAX),
	1    TBAR(K2MAX,L1MAX),TB(L1MAX)
	REAL ZG(L1MAX), ZR(KM), ZF(72)

	INTEGER INDX1(L1MAX), INDX2(KM), I1, I2
	REAL WT1(L1MAX), WT2(KM)
	INTEGER I, J, K, L, M, N
	LOGICAL FIRST /.TRUE./

	REAL HEAT(K1MAX,K2MAX,L1MAX),FORC(K1MAX, K2MAX)
	COMPLEX FORCESP(0:M1MAX,0:N1MAX,L1MAX)


	INTEGER KOUT
	DOUBLE PRECISION ZKM2K11,T2K11,XM2K11,P2K11,RHO2K11
	1    ,TOTN2K11,XCO2K11,XH2OK11,XO2K11,XOK11,XO3K11
	2    ,XCH4K11,XN2K11,XHK11,XH2K11,XARK11,XHEK11

	real hr, hr0,pi,HR2RAD
	COMPLEX IC,FACT,EXM(0:M1MAX)
	PARAMETER (PI=3.14159265)
	INTEGER IMAX

        character*80 arg
	INTEGER IARG
	common /opens/ iarg
C------------------------------------------------------------------
	if (mdlrd .eq. 1) return

	hr=tim/3600.

	IF (FIRST) THEN
c	  open(unit=21,file='d:/model/circ/output/zhudump.dat')
	  imax=k1

	  IC=CMPLX(0.,1.)
	  HR2RAD=PI/12.
	  hr0=tim/3600.

	  WLUN=224
	  CALL GETARG(IARG,ARG)
	  OPEN(UNIT=wlun,STATUS='OLD',FILE=ARG(1:(INDEX(ARG,' ')-1)),SHARED)
	  IARG=IARG+1
	  write(6,'(A)') arg

	  OPEN(UNIT=225,STATUS='OLD',
	1      FILE='input/numdens.dat',
	1      SHARED)
	  READ(225,*) (NUM(I),I=1,KM)
	  CLOSE(225)

	  DO I=1,KM
	    H2OTOP(I)=1E-5*NUM(I)
	  ENDDO

C...input gas mixing ratio
	  READ (WLUN,*) NLATW, NALTW
	  READ (WLUN,*) (CPRESW(I),I=1,NALTW), (CLATW(I),I=1,NLATW)
	  READ (WLUN,*) ((O3F(I,J), J=1,NLATW),I=1,NALTW)
	  READ (WLUN,*) ((O2F(I,J), J=1,NLATW),I=1,NALTW)
	  READ (WLUN,*) ((H2OF(I,J),J=1,NLATW),I=1,NALTW)
	  READ (WLUN,*) ((OF(I,J),  J=1,NLATW),I=1,NALTW)
	  READ (WLUN,*) ((CO2F(I,J),J=1,NLATW),I=1,NALTW)
	  CLOSE(WLUN)
	  CLOSE(WLUN)

C...interpolate mixing ratios onto heating alt grid and
c...model lat grid
	  DO L=1,NALTW
	    ZF(L)=-7.E3*ALOG(CPRESW(L)/1.E5)
	  ENDDO

	  DO I=1,KM
	    ZM(I)=(I-1)*2.D3
	    ZR(I)=(I-1)*2.E3
	  ENDDO

	  CALL  SPLIE2(ZF,CLATW,O3F,NALTW,NLATW,72,36,Y2A,IER)

	  DO K=1,K2
	    DO I=1,KM
	      CALL  SPLIN2(ZF,CLATW,O3F,Y2A,NALTW,NLATW,72,36,
	1	   ZR(I),DEGLAT(K),VAL,IER)
	      O3(I,K)=MAX(3E3,VAL*NUM(I))
	    ENDDO
	  ENDDO

	  CALL  SPLIE2(ZF,CLATW,O2F,NALTW,NLATW,72,36,Y2A,IER)

	  DO K=1,K2
	    DO I=1,KM
	      CALL  SPLIN2(ZF,CLATW,O2F,Y2A,NALTW,NLATW,72,36,
	1	   ZR(I),DEGLAT(K),VAL,IER)
	      O2(I,K)=VAL*NUM(I)
	    ENDDO
	  ENDDO

	  CALL  SPLIE2(ZF,CLATW,OF,NALTW,NLATW,72,36,Y2A,IER)

	  DO K=1,K2
	    DO I=1,KM
	      CALL  SPLIN2(ZF,CLATW,OF,Y2A,NALTW,NLATW,72,36,
	1	   ZR(I),DEGLAT(K),VAL,IER)
	      O(I,K)=VAL*NUM(I)
	    ENDDO
	  ENDDO

	  CALL  SPLIE2(ZF,CLATW,H2OF,NALTW,NLATW,72,36,Y2A,IER)

	  DO K=1,K2
	    DO I=1,KM
	      CALL  SPLIN2(ZF,CLATW,H2OF,Y2A,NALTW,NLATW,72,36,
	1	   ZR(I),DEGLAT(K),VAL,IER)
	      H2O(I,K)=VAL*NUM(I)
c	      H2O(I,K)=0.
	    ENDDO
	  ENDDO

	  CALL  SPLIE2(ZF,CLATW,CO2F,NALTW,NLATW,72,36,Y2A,IER)

	  DO K=1,K2
	    DO I=1,KM
	      CALL  SPLIN2(ZF,CLATW,CO2F,Y2A,NALTW,NLATW,72,36,
	1	   ZR(I),DEGLAT(K),VAL,IER)
	      CO2(I,K)=VAL*NUM(I)
	    ENDDO
	  ENDDO

C...MODEL ALT GRID (PUT IN INCREASING ORDER)
	  DO L=1,L1
	    ZG(L)=-7.E3*ALOG(PLV(L1+1-L)*PSURF/1.E5)

	    INDX1(L)=1
	    DO J=1,KM-1
	      IF (ZM(J) .LT. ZG(L)) INDX1(L)=J
	    ENDDO

	    I1=INDX1(L)
	    I2=I1+1
	    WT1(L)=(ZG(L)-ZM(I1)) / (ZM(I2)-ZM(I1))
	  ENDDO

	  DO L=1,KM
	    INDX2(L)=1
	    DO J=1,L1-1
	      IF (ZG(J) .LT. ZM(L)) INDX2(L)=J
	    ENDDO

	    I1=INDX2(L)
	    I2=I1+1
	    WT2(L)=(ZM(L)-ZG(I1)) / (ZG(I2)-ZG(I1))
	  ENDDO

C   Calculate (approximately) the column density of O2 and O3.
	  DO L=1,K2
	    BO2(KM,L)=O2(KM,L)*10.0D3 !  ZM(KM)>ZM(1)
	    BO3(KM,L)=O3(KM,L)*3.1D3
	    DO I=2,KM
	      K=KM-I+1
	      BO2(K,L)=BO2(K+1,L)+0.5D0*(O2(K,L)+O2(K+1,L))*
	1	   (ZM(K+1)-ZM(K))
	      BO3(K,L)=BO3(K+1,L)+0.5D0*(O3(K,L)+O3(K+1,L))*
	2	   (ZM(K+1)-ZM(K))
	    ENDDO
	  ENDDO

	  FIRST=.FALSE.
	ENDIF

c...update heating once a day
C...COMPUTED FOR LOCAL TIME=0 AT LONGITUDE=0
	if (hr .ge. hr0+24.) hr0=hr
	if (hr .ne. hr0) go to 100

C...interpolate zonal avg temperatures onto heating alt grid
	DO L=1,L1
	  K=L1+1-L
	  CALL PHYSIC( TP(1,1,K), PGQSP1(0,0,JPOT,L) )
	  CALL ZAVGE(TBAR(1,K),TP(1,1,K))
        END DO

	IF (IMAX .EQ. 1) THEN
	  DO K=1,K2
	    DO L=1,L1
	      TP(1,K,L)=TBAR(K,L)
	    ENDDO
	  ENDDO
	ENDIF

C...BEGIN LATITUDE LOOP
	DO K=1,K2
	  DLAT=DEGLAT(K)

c...takes too long here!!
C...BEGIN LONGITUDE/LOCAL TIME LOOP
	  DO I=1,IMAX
	    HOUR = LAMBDA(I)/HR2RAD

	  DO L=1,KM
	    IF (ZM(L) .GT. ZG(L1)) THEN
	      TEM(L)=TP(I,K,L1)
c	      TEM(L)=TBAR(K,L1)
	    ELSE
	      I1=INDX2(L)
	      I2=I1+1
	      TEM(L)=(1.-WT2(L))*TP(I,K,I1)+WT2(L)*TP(I,K,I2)
c	      TEM(L)=(1.-WT2(L))*TBAR(K,I1)+WT2(L)*TBAR(K,I2)
	    ENDIF
	  ENDDO

c...IR cooling rate QIR:
C...PREVENT H2O COOLING FROM GOING HAYWIRE
	  DO L=1,KM
	    H2OLIM(L)=MIN(H2O(L,K),H2OTOP(L))
	  ENDDO

c	  CALL QIRMID(QIR,TEM,CO2(1,K),O3(1,K),H2O(1,K),
	  CALL QIRMID(QIR,TEM,CO2(1,K),O3(1,K),H2OLIM,
	1      O(1,K),KM,KLBD)

C...BEGIN LONGITUDE/LOCAL TIME LOOP
c	  DO I=1,IMAX
c	    HOUR = LAMBDA(I)/HR2RAD

c... Solar heating rate HSL:
	    CALL HSLMID(HSL,TEM,O3(1,K),H2O(1,K),O2(1,K),
	1	 BO3(1,K),BO2(1,K),KM,DEC,DLAT,HOUR,XMU,
	1	 MDLRD,LTMU,OME0,F107,THETA27)

	    DO L=1,KM
	      QNET(L)=HSL(L)-QIR(L) !! net heating rate in K day^-1
c	      QNET(L)=HSL(L) !! net heating rate in K day^-1
c	      QNET(L)=-QIR(L) !! net heating rate in K day^-1
	    ENDDO

	    DO L=1,L1
	      IF (ZG(L) .GT. ZM(KM)) THEN
		HEAT(I,K,L)=QNET(KM)
	      ELSE
		I1=INDX1(L)
		I2=I1+1
		HEAT(I,K,L)=(1.-WT1(L))*QNET(I1)+WT1(L)*QNET(I2)
	      ENDIF

	      heat(i,k,l)=heat(i,k,l)*.5*(1-tanh((zg(L)-80e3)/5.))
	    ENDDO
	  ENDDO			!END LON LOOP
	ENDDO			!END LAT LOOP

c	write(6,*) 'zhu dump'
c	WRITE(21,*) (((HEAT(I,K,L),I=1,IMAX),K=1,K2),L=1,L1)
c	WRITE(21,*) ((o3(L,K),K=1,K2),L=1,KM)
c	write(6,*) hr

	DO L=1,L1
	  DO K=1,K2
            DO I=1,K1
	      J=I
c	      IF (MDLRD .EQ. 2) J=1
              FORC(I,K)=HEAT(J,K,L1+1-L)/86400.
            END DO
          END DO

	  CALL ZEROSP(FORCESP(0,0,L),N1)
	  CALL SPECTR(FORCESP(0,0,L),FORC)
        END DO

c	DO L=1,L1-7
c	  DO N=0,N1
c	    DO M=1,M1
c	      FORCESP(M,N,L)=0.
c	    ENDDO
c	  ENDDO
c	ENDDO

c	WRITE(21,*) (((FORCESP(I,K,L),I=0,M1),K=0,N1),L=1,L1)

 100	CONTINUE

C...SPECTRAL COMPONENENTS OF FORCING.  SHIFT LOCAL TIME
	FACT=IC*HR*HR2RAD
c	if (mdlrd .eq. 2) fact=0.

	exm(0)=1.
	DO M=1,MIN0(N,M1)
	  EXM(M)=CEXP(M*FACT)*(1.+TANH((TIM-432000.)/172800.))/2.
	ENDDO

	DO L=1,L1
          DO N=0,N1
            DO M=0,MIN0(N,M1)
              DPGQSP(M,N,JPOT,L)=DPGQSP(M,N,JPOT,L)+
	1	   FORCESP(M,N,L)*EXM(M)
            END DO
          END DO
        END DO

	END
