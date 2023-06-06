!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!######################### COOLMD.FOR ########################
!############# IMPLICIT REAL*8 (A-H,O-Z) ############
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    QIRMID,   RADCO2,   DETIJ,    INVERT,   AVT21,    AVT32,    
!    SPT,      RADO3,    DETIJR,   RADH2O,   INPUTC,   NLTEC,    
!
!  To calculate the IR radiative cooling rate for middle atmosphere modeling (8/95)


!  ==== IR radiative cooling by CO2, H2O, and O3 in the middle atmosphere ==== 

! NB: This module has loads of implicit types....

module coolmd_mod

    contains

      SUBROUTINE QIRMID(QIR,TEM,CO2,O3,H2O,OXY,IM,KLBD)
        !  To calculate the IR cooling rate in the middle atmosphere (K day^-1)
              IMPLICIT REAL*8 (A-H,O-Z)
              PARAMETER (KM=76)         !!! specified by users: KM=IM   ZM(KM)>ZM(1)
              COMMON /CTXZ23/ Z0(KM),P0(KM),T0(KM),GDP(KM)
              DIMENSION QIR(IM),TEM(IM),CO2(IM),O3(IM),H2O(IM),OXY(IM)
              DIMENSION WK1(KM),WK2(KM),WK3(KM),WK4(KM) &
             &         ,PRE(KM),QCO2(KM),QO3(KM),QH2O(KM)
              DO 10 K=1,KM
          10  WK1(K)=P0(K)/(1.380658D-23*TEM(K))        ! air density in m^-3
              DO 20 K=1,KM
              PRE(K)=P0(K)
              WK2(K)=1.5225D0*CO2(K)/WK1(K)             ! CO2 mass mixing ratio
              WK3(K)=1.661D0*O3(K)/WK1(K)               ! ozone mass mixing ratio
              WK4(K)=0.622D0*H2O(K)/WK1(K)              ! water vapor mass mixing ratio
          20  CONTINUE
              JCON=2       ! 1 => linear interpolation; 2 => quadrature interpolation
              JSEP=1       ! 0 => LTE; 1 => nonLTE in z>70 km; 2 => complete nonLTE
              CALL RADCO2(TEM,WK2,OXY,QCO2,IM,JCON,JSEP)
              CALL RADO3(KLBD,TEM,WK3,QO3,IM,JSEP)
              CALL RADH2O(TEM,WK4,QH2O,IM)
              DO 25 K=1,KM
              IF(Z0(K).GT.25.0D3) GO TO 25              ! CO2 heating at equatorial tropopause
              IF(QCO2(K).GT.0.0D0) QCO2(K)=0.0D0        ! is balanced by cirrus cooling
          25  CONTINUE
              DO 30 K=1,KM
          30  QIR(K)=-QCO2(K)-QO3(K)-QH2O(K)   ! total IR cooling rate in K/day (QIR > 0: cooling)
              RETURN
              END
              
              SUBROUTINE QIRMID_LTE(QIR,TEM,CO2,O3,H2O,OXY,IM,KLBD)
        !  To calculate the IR cooling rate in the middle atmosphere (K day^-1)
              IMPLICIT REAL*8 (A-H,O-Z)
              PARAMETER (KM=76)         !!! specified by users: KM=IM   ZM(KM)>ZM(1)
              COMMON /CTXZ23/ Z0(KM),P0(KM),T0(KM),GDP(KM)
              DIMENSION QIR(IM),TEM(IM),CO2(IM),O3(IM),H2O(IM),OXY(IM)
              DIMENSION WK1(KM),WK2(KM),WK3(KM),WK4(KM) &
             &         ,PRE(KM),QCO2(KM),QO3(KM),QH2O(KM)
              DO 10 K=1,KM
          10  WK1(K)=P0(K)/(1.380658D-23*TEM(K))        ! air density in m^-3
              DO 20 K=1,KM
              PRE(K)=P0(K)
              WK2(K)=1.5225D0*CO2(K)/WK1(K)             ! CO2 mass mixing ratio
              WK3(K)=1.661D0*O3(K)/WK1(K)               ! ozone mass mixing ratio
              WK4(K)=0.622D0*H2O(K)/WK1(K)              ! water vapor mass mixing ratio
          20  CONTINUE
              JCON=2       ! 1 => linear interpolation; 2 => quadrature interpolation
              JSEP=0       ! 0 => LTE; 1 => nonLTE in z>70 km; 2 => complete nonLTE
              CALL RADCO2(TEM,WK2,OXY,QCO2,IM,JCON,JSEP)
              CALL RADO3(KLBD,TEM,WK3,QO3,IM,JSEP)
              CALL RADH2O(TEM,WK4,QH2O,IM)
              DO 25 K=1,KM
              IF(Z0(K).GT.25.0D3) GO TO 25              ! CO2 heating at equatorial tropopause
              IF(QCO2(K).GT.0.0D0) QCO2(K)=0.0D0        ! is balanced by cirrus cooling
          25  CONTINUE
              DO 30 K=1,KM
          30  QIR(K)=-QCO2(K)-QO3(K)-QH2O(K)   ! total IR cooling rate in K/day (QIR > 0: cooling)
              RETURN
              END
              
        !  ========================= CO2 cooling : =================================   
        
              SUBROUTINE RADCO2(T,RR,BOXY,QCO2,IM,JCON,JSEP)
        !  To calculates the cooling rate (QCO2) by Curtis matrix CURT for
        !  given p, T, [O], and r[CO2] deviated from [CO2]0 around top levels. 
        !  CURTZ are the off-line calculated Curtis matrices. RR is mass mixing ratio.
        !  The units of [O] and p are  M^-3 and Pascal, respectively.
        !  JCON=1 --> linear interpolation; JCON=2 --> quadrature interpolation.
              IMPLICIT REAL*8 (A-H,O-Z)
              PARAMETER (KM=76)    !!! specified by users:  KM=IM,   ZM(KM)>ZM(1)
              COMMON /CTXZ23/ Z0(KM),P0(KM),T0(KM),GDP(KM)
              COMMON /CURTZ2/ R20(KM),H200(KM),CURT0(KM,KM) &
             & ,A10(KM,KM),A01(KM,KM),A20(KM,KM),A02(KM,KM),A11(KM,KM)
              DIMENSION CURT(KM,KM),THE1(KM),EE(KM,KM),BAIR(KM),DTJ(KM)
        !  EE(KM,KM) is used both for diagonal E matrix and non-LTE Curtix matrix
              DIMENSION T(IM),RR(IM),BOXY(IM),QCO2(IM)
              DATA AA21/1.43D0/,AA32/1.94D0/
              DO 5 J=1,IM
              QCO2(J)=0.0D0
              THE1(J)=47.612466D0/(DEXP(970.97D0/T(J))-1.0D0)  ! 47.612466=EXP(970.97/250)-1
           5  DTJ(J)=(T(J)-T0(J))
              DO 10 K=1,IM
              DO 10 J=1,IM
              EE(J,K)=0.0D0
          10  CURT(J,K)=CURT0(J,K)
              DO 20 K=1,IM        ! Get the Curtis matrix by interpolation
              DO 20 J=1,IM
              IF(DABS(CURT0(J,K)).LE.0.005D0) GO TO 20
              DELTA=DETIJ(T0,T,GDP,IM,J,K)
              CURT(J,K)=CURT0(J,K)+A10(J,K)*DTJ(J)+A01(J,K)*DELTA
              IF(JCON.EQ.2) CURT(J,K)=CURT(J,K)+0.5D0*(A20(J,K) &
             & *DTJ(J)**2+A02(J,K)*DELTA**2)+A11(J,K)*DELTA*DTJ(J)
          20  CONTINUE
              IF(JSEP.EQ.0) THEN
              CALL NLTEC(CURT,EE,IM,JSEP)
              GO TO 300
              ENDIF
              DO 28 J=1,IM
          28  BAIR(J)=P0(J)/(1.380658D-23*T(J))              ! air density in m^-3
              DO 30 J=1,IM
              PHIN=((SPT(T(J),1)+SPT(T(J),3)) &
             & *AVT21(J,BAIR,BOXY,T,IM)/AA21+SPT(T(J),2) &
             & *AVT32(J,BAIR,BOXY,T,IM)/AA32)/SPT(T(J),4)
          30  EE(J,J)=R20(J)/(RR(J)*2.0D0*H200(J)*PHIN)            ! EJJ=R0/(R*2*H00)
              CALL NLTEC(CURT,EE,IM,JSEP)
         300  CONTINUE
              DO 55 K=1,IM
              DO 55 J=1,IM
          55  QCO2(J)=QCO2(J)+EE(J,K)*THE1(K)
              DO 60 J=1,IM
              IF(Z0(J).LE.80.0D3) GO TO 60
              FACX=0.5D0*(1.0D0+TANH((Z0(J)-85.0D3)/7.0D3))
              QCO2(J)=QCO2(J)*((RR(J)/R20(J))*FACX+(1.0-FACX))
          60  CONTINUE
              IF(Z0(IM).LE.100.01D3) THEN   ! Aug. 95: correction around the top boundary
              DO 70 K=1,IM
              XXZ1=(Z0(IM)+1.0D3-Z0(K))/5.0D3
          70  QCO2(K)=QCO2(K)*(1.0D0-DEXP(-XXZ1**2))
              ENDIF
              RETURN
              END
          
              FUNCTION DETIJ(T0,T,GDP,IM,II,JJ)
        !   To calculate the DELTA(i,j) used for interpolations (FS81)
              IMPLICIT REAL*8 (A-H,O-Z)
              DIMENSION T0(IM),T(IM),GDP(IM)
              K1=MIN0(II,JJ)
              K2=MAX0(II,JJ)
              XX=0.0D0
              YY=0.0D0
              DO 10 K=K1,K2
              XX=XX+GDP(K)
          10  YY=YY+(T(K)-T0(K))*GDP(K)
              DETIJ=YY/XX
              RETURN
              END
          
              SUBROUTINE INVERT(D,N,M)
        !   To invert the matrix D,  M=N+1. Mth row can be any values. REAL*8
              IMPLICIT REAL*8 (A-H,O-Z)
              DIMENSION D(N,M)
              KK=0
              JJ=0
              DO 10 K=1,N
              DO 11 J=1,N
          11  D(J,M)=0.0D0
              D(K,M)=1.0D0
              JJ=KK+1
              LL=JJ
              KK=KK+1
          20  IF(DABS(D(JJ,KK)-1.0D-30)) 21,21,22
          21  JJ=JJ+1
              IF(JJ-N) 20,20,99
          99  WRITE(*,98)
          98  FORMAT('ERRORNEOUS INPUT')
              RETURN
          22  IF(LL-JJ) 23,24,23
          23  DO 25 MM=1,M
              DTEMP=D(JJ,MM)
              D(LL,MM)=D(JJ,MM)
          25  D(JJ,MM)=DTEMP
          24  DIV=D(K,K)
              DO 30 LJ=1,M
          30  D(K,LJ)=D(K,LJ)/DIV
              DO 12 I=1,N
              FAC=D(I,K)
              IF(I-K) 15,12,15
          15  DO 31 LJ=1,M
          31  D(I,LJ)=D(I,LJ)-FAC*D(K,LJ)
          12  CONTINUE
              DO 40 J=1,N
          40  D(J,K)=D(J,M)
          10  CONTINUE
              RETURN
              END
        
              FUNCTION AVT21(I,BAIR,BOXY,TEM,IM)
              IMPLICIT REAL*8 (A-H,O-Z)
              DIMENSION BAIR(IM),BOXY(IM),TEM(IM)
              IF(TEM(I).LT.200.0D0) VT1M=2.5D-21
              IF(TEM(I).GE.200.0D0) VT1M=2.5D-21*(1.0D0+0.03D0*(TEM(I)-200.D0))
              AVT21=BAIR(I)*VT1M+BOXY(I)*5.5D-18      ! revised k-rate for [O]
              RETURN
              END
              
              FUNCTION AVT32(I,BAIR,BOXY,TEM,IM)
              IMPLICIT REAL*8 (A-H,O-Z)
              DIMENSION BAIR(IM),BOXY(IM),TEM(IM)
              TRA=TEM(I)/273.3D0
              VT2M=1.24D-20*TRA*TRA
              AVT32=BAIR(I)*VT2M+BOXY(I)*5.5D-18      ! revised k-rate for [O]
              RETURN
              END
        
              FUNCTION SPT(T,N)
        !   To calculate the temperature dependence of band intensity S
        !   in Houghton's appedices. N=1,2,3 and 4 correspond to 626 
        !   fundamental band, 626 first hot band, isotopes fundamental 
        !   band and total 15 micron band, respectively. Only used in
        !   non-LTE calculations.
              IMPLICIT REAL*8 (A-H,O-Z)
              TX=T
              IF(TX.GT.320.0) TX=350.0D0
              IF(TX.LT.130.0) TX=130.0D0
              X=TX-250.0D0
              GO TO (10,20,30,40),N
          10  Y=5.0522D0-5.4693D-4*X-5.2873D-7*X*X
              GO TO 100
          20  Y=3.8067D0+7.0636D-3*X-2.1024D-5*X*X+1.964D-8*X*X*X
              GO TO 100
          30  Y=3.2338D0-5.5612D-4*X-5.3153D-7*X*X
              GO TO 100
          40  Y=5.0961D0+8.8837D-6*X+3.1662D-8*X*X
         100  SPT=10.0D0**Y
              RETURN
              END
        
        ! ======================== Ozone cooling : ==============================
        
              SUBROUTINE RADO3(KLBD,T,RR,QO3,IM,JSEP)
        !  To calculate the O3 cooling rate by Curtis matrix CURT.  !  ZM(M)>ZM(1)
        !  KLBD=lower boundary grid, KLBD=1 => surface as a lower boundary.
              IMPLICIT REAL*8 (A-H,O-Z)     !!! (KM,LBTRP)=(41,8),(60,12),(75,16)
              PARAMETER (KM=76,LBTRP=12)    !!! specified by users:  KM=IM,   ZM(KM)>ZM(1)
              COMMON /CTXZ23/ Z0(KM),P0(KM),T0(KM),GDP(KM)
              COMMON /CURTZ3/ R30(KM),H300(KM),ZCURT0(KM,KM) &
             & ,ZA10(KM,KM),ZA01(KM,KM),ZB01(KM,KM),CLB0(KM,LBTRP) &
             & ,A10LB(KM,LBTRP),A01LB(KM,LBTRP),B01LB(KM,LBTRP)
              DIMENSION CURT(KM,KM),THE1(KM),EE(KM,KM),CLB3(KM),DTJ(KM)
        !  EE(KM,KM) is used both for diagonal E matrix and non-LTE Curtix matrix
              DIMENSION T(IM),RR(IM),QO3(IM),WK1(KM)
              DO 5 J=1,IM
              THE1(J)=484.0D0/(DEXP(1546.0D0/T(J))-1.0D0)
              DTJ(J)=(T(J)-T0(J))
           5  CLB3(J)=CLB0(J,KLBD)+A10LB(J,KLBD)*DTJ(J) &
             & +A01LB(J,KLBD)*DETIJ(T0,T,GDP,KM,J,KLBD) &
             & +B01LB(J,KLBD)*DETIJR(R30,RR,GDP,KM,J,KLBD,0.2D0)
              DO 10 K=1,IM
              DO 10 J=1,IM
          10  CURT(J,K)=ZCURT0(J,K)
              DO 20 K=1,KM
              DO 20 J=1,KM        ! Get the Curtis matrix by interpolation
              IF(DABS(ZCURT0(J,K)).LE.0.005D0) GO TO 20
              DELTAT=DETIJ(T0,T,GDP,KM,J,K)
              DELTAM=DETIJR(R30,RR,GDP,KM,J,K,0.2D0)
              CURT(J,K)=ZCURT0(J,K)+ZA10(J,K)*DTJ(J) &
             &  +ZA01(J,K)*DETIJ(T0,T,GDP,KM,J,K)    &
             &  +ZB01(J,K)*DETIJR(R30,RR,GDP,KM,J,K,0.2D0)
          20  CONTINUE
              DO 25 J=1,KM
              WK1(J)=(RR(J)/R30(J))*H300(J)
          25  CLB3(J)=CLB3(J)*WK1(J)
              DO 30 K=1,KM
              DO 30 J=1,KM
          30  EE(J,K)=CURT(J,K)*WK1(J)
              DO 40 K=1,KM
              DO 40 J=1,KM
              CURT(J,K)=EE(J,K)
          40  EE(J,K)=0.0D0
              IF(JSEP.EQ.0) THEN
              CALL NLTEC(CURT,EE,IM,JSEP)
              GO TO 300
              ENDIF
              DO 50 J=1,IM
              WK1(J)=3.64D-21*P0(J)/(1.380658D-23*T(J))    ! PHIN=(4.3E-20/11.81)*air density
          50  EE(J,J)=R30(J)/(RR(J)*2.0D0*H300(J)*WK1(J))  ! EJJ=R0/(R*2*H00)
              CALL NLTEC(CURT,EE,IM,JSEP)
         300  CONTINUE
              DO 60 J=1,IM
          60  QO3(J)=CLB3(J)*THE1(KLBD)
              DO 70 K=KLBD,KM
              DO 70 J=1,IM
          70  QO3(J)=QO3(J)+EE(J,K)*THE1(K)
              DO 80 J=1,IM
              IF(J.LT.KLBD) QO3(J)=0.0D0
          80  CONTINUE
              RETURN
              END
          
              FUNCTION DETIJR(R0,RR,GDP,IM,II,JJ,DELR0)
        !   To calculate the DELTA(i,j) used for mixing ratio interpolations (ZHU, 1993)
              IMPLICIT REAL*8 (A-H,O-Z)
              DIMENSION R0(IM),RR(IM),GDP(IM)
              K1=MIN0(II,JJ)
              K2=MAX0(II,JJ)
              XX=0.0D0
              YY=0.0D0
              DO 10 K=K1,K2
              XX=XX+DELR0*R0(K)*GDP(K)
          10  YY=YY+(RR(K)-R0(K))*GDP(K)
              DETIJR=YY/XX
              RETURN
              END
              
        ! ======================== Water vapor cooling : ========================
        
              SUBROUTINE RADH2O(T,RR,QH2O,IM)
        !  To calculate the H2O cool-to-space cooling rate.  !  ZM(KM)>ZM(1)
              IMPLICIT REAL*8 (A-H,O-Z)
              PARAMETER (KM=76)     !!! specified by users:  KM=IM,   ZM(KM)>ZM(1)
              COMMON /CTXZ23/ Z0(KM),P0(KM),T0(KM),GDP(KM)
              COMMON /GAMWAT/ R10(KM),GAMR1(KM),GAMV1(KM) &
               ,YAR1(KM),YAV1(KM),YAR2(KM),YAV2(KM)
              DIMENSION GR2(KM),GV2(KM),THR(KM),THV(KM),PHIV(KM)
              DIMENSION T(IM),RR(IM),QH2O(IM),WK1(KM),WK2(KM)
              DO 10 K=1,KM
              WK1(K)=DETIJR(R10,RR,GDP,KM,K,KM,0.5D0)
              WK2(K)=RR(K)/3.0D-6         ! 3.0E-6 => reference mass mixing ratio
          10  PHIV(K)=2.1D-21*P0(K)/(1.380658D-23*T(K)) ! (4.3E-20/20.5)*air density in m^-3
              DO 20 K=1,KM 
              THR(K)=81.6D0/(DEXP(568.01D0/T(K))-1.0D0)          ! S*B in K/day
              THV(K)=2730.0D0/(DEXP(2300.8D0/T(K))-1.0D0)        ! S*B in K/day
              GR2(K)=GAMR1(K)+YAR1(K)*WK1(K)+YAR2(K)*WK1(K)**2
          20  GV2(K)=GAMV1(K)+YAV1(K)*WK1(K)+YAV2(K)*WK1(K)**2
              DO 60 K=1,KM
          60  QH2O(K)=-WK2(K)*(THR(K)*GR2(K) &
                      +THV(K)*GV2(K)/(1.0D0+0.5D0*GV2(K)/PHIV(K)))
              RETURN
              END
              
        
        !  ==================== Input of the Curtis matrices ====================
        
              SUBROUTINE INPUTC
        ! --- Input the Curtis matrices for cooling rate calculations
              IMPLICIT REAL*8 (A-H,O-Z)     !!! (KM,LBTRP)=(41,8),(60,12),(75,16)
              PARAMETER (KM=76,LBTRP=12)    !!! specified by users:  KM=IM,   ZM(KM)>ZM(1)
              COMMON /CTXZ23/ Z0(KM),P0(KM),T0(KM),GDP(KM)
              COMMON /GAMWAT/ R10(KM),GAMR1(KM),GAMV1(KM) &
             & ,YAR1(KM),YAV1(KM),YAR2(KM),YAV2(KM)
              COMMON /CURTZ2/ R20(KM),H200(KM),CURT0(KM,KM) &
             & ,A10(KM,KM),A01(KM,KM),A20(KM,KM),A02(KM,KM),A11(KM,KM)
              COMMON /CURTZ3/ R30(KM),H300(KM),ZCURT0(KM,KM) &
             & ,ZA10(KM,KM),ZA01(KM,KM),ZB01(KM,KM),CLB0(KM,LBTRP) &
             & ,A10LB(KM,LBTRP),A01LB(KM,LBTRP),B01LB(KM,LBTRP)
             
        !     Input and output folders and files
              CHARACTER*50 INFOLD, INFILE, OUTFOLD, OUTFILE, TIDEFILE, TTIDEFILE
              CHARACTER*50 INITTEMP, INITZONAL, INITGRID, INITZGRID, GWFDUMP, CURT
              CHARACTER*50 SURFFILE
              COMMON /FILES/ INFOLD, INFILE, SURFFILE, OUTFOLD, OUTFILE, &
                             INITTEMP, INITZONAL, INITGRID, GWFDUMP, TIDEFILE, TTIDEFILE,CURT
             
          52  FORMAT(8E12.5)
          53  FORMAT(2F12.6,5E12.3)
          55  FORMAT(7E14.5)
              OPEN(15,FILE=TRIM(INFOLD)//TRIM(CURT),STATUS='OLD')    ! All the Curtis matrices
              DO 10 K=1,KM
              READ(15,52) Z0(K),T0(K),P0(K),GDP(K), &
             &             R20(K),H200(K),R30(K),H300(K)
          10  READ(15,53) GAMR1(K),GAMV1(K),R10(K), &
             &             YAR1(K),YAV1(K),YAR2(K),YAV2(K)
              READ(15,55) ((CURT0(L1,L2),L1=1,KM),L2=1,KM)
              READ(15,55) ((A10(L1,L2),L1=1,KM),L2=1,KM)
              READ(15,55) ((A01(L1,L2),L1=1,KM),L2=1,KM)
              READ(15,55) ((A20(L1,L2),L1=1,KM),L2=1,KM)
              READ(15,55) ((A02(L1,L2),L1=1,KM),L2=1,KM)
              READ(15,55) ((A11(L1,L2),L1=1,KM),L2=1,KM)
              READ(15,55) ((ZCURT0(L1,L2),L1=1,KM),L2=1,KM)
              READ(15,55) ((ZA10(L1,L2),L1=1,KM),L2=1,KM)
              READ(15,55) ((ZA01(L1,L2),L1=1,KM),L2=1,KM)
              READ(15,55) ((ZB01(L1,L2),L1=1,KM),L2=1,KM)
              READ(15,55) ((CLB0(L1,L2),A10LB(L1,L2) &
             &   ,A01LB(L1,L2),B01LB(L1,L2),L1=1,KM),L2=1,LBTRP)
              CLOSE(UNIT=15)
              RETURN
              END SUBROUTINE INPUTC
          
              
        ! ==================== get non-LTE Curtis matrix : =======================
        
              SUBROUTINE NLTEC(CURT,EE,IM,JSEP)
        !  To calculate nonLTE curtis matrix from matrices C and E (JGR, 97, 12,790).
        !  The output nonLTE curtis matrix is stored in EE.        !  ZM(KM)>ZM(1)
        !  JSEP=0 => LTE; JSEP=1 => nonLTE in top K1 region; JSEP>1 => complete nonLTE
              IMPLICIT REAL*8 (A-H,O-Z)   !! {(KM,K1)=(41,21),(60,27),(75,30)}
              PARAMETER (KM=76,K1=36)    !!! specified by users:  KM=IM,   ZM(KM)>ZM(1)
              PARAMETER (KMP=KM+1,K1P=K1+1,K2=KM-K1)
              DIMENSION CURT(IM,IM),EE(IM,IM),AINV(KM,KMP)
              DIMENSION BINV(K1,K1P),WOK(K2,K1)
              IF(JSEP.EQ.0) THEN
              DO 26 K=1,IM
              DO 26 J=1,IM
          26  EE(J,K)=CURT(J,K)
              RETURN
              ENDIF
              IF(JSEP.EQ.1) GO TO 100
              DO 40 K=1,IM
              DO 40 J=1,IM
              XX=0.0D0
              DO 35 LL=1,IM
          35  XX=XX+CURT(J,LL)*EE(LL,K)
          40  AINV(J,K)=-XX
              DO 41 K=1,IM
          41  AINV(K,K)=AINV(K,K)+1.0D0
              CALL INVERT(AINV,KM,KMP)
              GO TO 200
         100  CONTINUE
              DO 140 K=1,K1
              DO 140 J=1,K1
              XX=0.0D0
              DO 135 LL=1,K1
         135  XX=XX+CURT(K2+J,K2+LL)*EE(K2+LL,K2+K)
         140  BINV(J,K)=-XX
              DO 141 K=1,K1
         141  BINV(K,K)=BINV(K,K)+1.0D0
              CALL INVERT(BINV,K1,K1P)
              DO 150 K=1,K2
              DO 150 J=1,IM
         150  AINV(J,K)=0.0D0
              DO 151 K=1,K2
         151  AINV(K,K)=1.0D0
              DO 160 K=1,K1
              DO 160 J=1,K1
         160  AINV(K2+J,K2+K)=BINV(J,K)
              DO 170 K=1,K1
              DO 170 J=1,K2
              XX=0.0D0
              DO 165 LL=1,K1
         165  XX=XX+CURT(J,K2+LL)*EE(K2+LL,K2+K)
              WOK(J,K)=XX
         170  CONTINUE
              DO 180 K=1,K1
              DO 180 J=1,K2
              XX=0.0D0
              DO 175 LL=1,K1
         175  XX=XX+WOK(J,LL)*BINV(LL,K)
              AINV(J,K2+K)=XX
         180  CONTINUE
         200  CONTINUE
              DO 50 K=1,IM
              DO 50 J=1,IM
              XX=0.0D0
              DO 45 LL=1,IM
          45  XX=XX+AINV(J,LL)*CURT(LL,K)
          50  EE(J,K)=XX
              RETURN
              END
            
end module coolmd_mod