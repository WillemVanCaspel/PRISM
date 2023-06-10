C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C########################## HEATMD.FOR ##########################
C############## IMPLICIT REAL*8 (A-H,O-Z) ##############
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    HSLMID,  HTH2O,   HEAT6,   HEAT3,   HTMU,    HT12,    EFFCY
C    xsrcq,
C
C  To calculate solar heating rate for middle atmosphere modeling (8/95)
c
c  The code was updated on 08/16/2000 to include the solar
c  variability dependence of the 11-year solar cycle (f107) and 
c  27-day solar rotational variability.


CC ==== solar heating by O3, O2, and tropospheric H2O (K day^-1) ==== 

      SUBROUTINE HSLMID(HSL,TEM,O3,H2O,O2,BO3,BO2,IM
     &   ,DEC,PHI,HOUR,XMU,MDEL,LTMU,OME0,f107,theta27)
CC  to calculate solar heating rate by O3, O2, and tropospheric H2O (K day^-1) 
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KM=76)     !!! specified by users:  KM=IM,   ZM(KM)>ZM(1)
      COMMON /CTXZ23/ Z0(KM),P0(KM),T0(KM),GDP(KM)
      DIMENSION HSL(IM),O3(IM),H2O(IM),O2(IM),BO3(IM),BO2(IM)
      DIMENSION HUVIR(KM),HTROP(KM)
      ITOP=0                       ! IM is the top boundary
      IF(MDEL.EQ.2) THEN           ! diurnally averaged heating rate
      DO 10 K=1,KM
  10  HTROP(K)=0.0D0 ! zonal averaged tropospheric H2O heating rate vanishes
      CALL HEAT6(HUVIR,DEC,PHI,O2,O3,BO2,BO3,P0,TEM
     &          ,IM,ITOP,OME0,f107,theta27)
      ENDIF
      IF(MDEL.EQ.3) THEN           ! local heating rate
      CALL HEAT3(HUVIR,HOUR,DEC,PHI,XMU,LTMU
     &           ,O2,O3,BO2,BO3,P0,TEM,IM,ITOP,OME0,f107,theta27)      
      RG0=0.2D0          ! surface albedo to solar near-IR flux
      CALL HTH2O(HTROP,HOUR,DEC,PHI,P0,TEM,H2O,IM,RG0)  ! tropospheric H2O heating rate
      ENDIF
      DO 20 K=1,KM
  20  HSL(K)=HUVIR(K)+HTROP(K)            ! heating rate in K/day
      RETURN
      END
      

C  ============== local H2O heating ========================

      SUBROUTINE HTH2O(HTW,TIME,DEC,PHI,PRE,TEM,H2O,IM,RG0)
CCCCC
C   Heating rate (MKS units) for H2O: partly follows parameterization
C   scheme by M. D. Chou and K. T. Lee (1995). H2O=number density (m^-3).
CCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (KM=76)     !!! specified by users:  KM=IM,   ZM(KM)>ZM(1)
      PARAMETER (MM=10)     !  PRE(KM) = top boundary pressure
      DIMENSION HTW(IM),PRE(IM),TEM(IM),H2O(IM),XXQ(KM)
      DIMENSION WRK1(KM),WRK2(KM),WRK5(MM)
      COMMON /KDGH2O/ YKH2O(10),DGKH2O(10)
      SOLAR0=1373.0D0     ! solar constant in W m^-2
      DO 10 I=1,IM
  10  HTW(I)=0.0D0
      RDEC=DEC/57.2958D0
      RPHI=PHI/57.2958D0
      ACB=DSIN(RDEC)*DSIN(RPHI)
      BCB=DCOS(RDEC)*DCOS(RPHI)
      HHX=(TIME-12.0D0)*0.261799D0               ! PI/12=0.261799
      BMU1=ACB+BCB*DCOS(HHX)
      IF(BMU1.LE.0.0D0) RETURN                   ! The sun does not rise
      PRER0=3.0D4                ! reference pressure in pascal
      TEMR0=240.0D0              ! reference temperature in K
      BC=1.380658D-23
      XRM=35.0D0/DSQRT(1224.0D0*BMU1**2+1.0D0)   ! magnification factor => 1/BMU1
      DO 20 I=1,IM
      FACT1=(PRE(I)/PRER0)**0.8                  ! [(p / pr)^m] 
      FACT2=DEXP(1.35D-3*(TEM(I)-TEMR0))         ! h(T, Tr)
      TOTN=PRE(I)/(BC*TEM(I))                    ! air number density (m^-3)
  20  XXQ(I)=FACT1*FACT2*0.622D0*H2O(I)/TOTN     ! weighted specific humidity
      WRK1(IM)=0.0D0
      DO 50 I=2,IM
      KI=IM+1-I
      DELP=PRE(KI)-PRE(KI+1)
  50  WRK1(KI)=WRK1(KI+1)+(XXQ(KI)+XXQ(KI+1))*DELP
      XFAC=XRM/19.62D0    ! 2.0*9.81 = 19.62,  slant path length in [kg m^-2]
      DO 55 I=1,IM
  55  WRK1(I)=WRK1(I)*XFAC
      WSPS=WRK1(1)         ! total slant path length in [kg m^-2]
      WRK2(1)=0.0D0
      DO 60 I=2,IM
      DELP=PRE(I-1)-PRE(I)
  60  WRK2(I)=WRK2(I-1)+(XXQ(I-1)+XXQ(I))*DELP
      YFAC=1.67D0/19.62D0
      DO 65 I=1,IM
  65  WRK2(I)=WSPS+WRK2(I)*YFAC  ! diffusive path length in [kg m^-2] ==> [w sup *](p)
      CONV12=0.1D0               ! conversion factor from [cm^2 g^-1] to [m^2 kg^-1]
      DO 80 M=1,MM
  80  WRK5(M)=YKH2O(M)*CONV12         ! weighted k coefficient in [m^2 kg^-1]
      DO 150 I=1,IM
      DO 150 M=1,MM
      HTW(I)=HTW(I)+(DEXP(-WRK5(M)*WRK1(I))+(1.67D0*RG0/XRM)*
     &               DEXP(-WRK5(M)*WRK2(I)))*WRK5(M)*DGKH2O(M)
 150  CONTINUE
      XCOE1=(SOLAR0/1004.0D0)*86400.0D0        ! heating rate in K day^-1
      DO 200 I=1,IM
 200  HTW(I)=XCOE1*XXQ(I)*HTW(I)               ! heating rate in K day^-1
      RETURN
      END
      

C  ================ diurnally averaged O2 and O3 heating ================

      SUBROUTINE HEAT6(QH,DEC,PHI,O2,O3,BO2,BO3,PRE,TEM
     &                ,IM,ITOP,OME0,f107,theta27)
CCCCC
C       To calculate the daily averaged heating rate (QH: K/day) for O3 
C       and O2, following the parameterization scheme by Strobel (1978).
C       DEC=declination (degree), PHI=latitude (degree).
C       O2,O3=number density (m^-3), BO2,BO3=column density (m^-2).
C       PRE=pressure (Pa), TEM=temperature. FCLEAR=Fraction of clear sky.
C       ITOP=1 => PRE(1)=TOP BOUNDARY;  ITOP=0 => PRE(IM)=TOP BOUNDARY
CCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION QH(IM),O2(IM),O3(IM),BO2(IM),BO3(IM),PRE(IM)
     &         ,TEM(IM)
      RDEC=DEC/57.2958D0
      RPHI=PHI/57.2958D0
      ACB=DSIN(RDEC)*DSIN(RPHI)
      BCB=DCOS(RDEC)*DCOS(RPHI)
      BCB2=BCB*BCB
      DO 10 I=1,IM
  10  QH(I)=0.0D0
      XA=0.0D0
      XB=ACB+BCB
      XC=ACB-BCB
      IF(XC.GT.0.0D0) XA=XC              !    The sun does not set
      IF(XB-XA.LE.0.0D0) RETURN          !    The sun does not rise
      XM=0.5D0*(XB+XA)
      XR=0.5D0*(XB-XA)
      DDX1=0.1252334D0*XR
      DDX2=0.3678315D0*XR
      DDX3=0.5873179D0*XR
      DDX4=0.7699026D0*XR
      DDX5=0.9041172D0*XR
      DDX6=0.9815606D0*XR
      WW1=0.2491470D0
      WW2=0.2334925D0
      WW3=0.2031674D0
      WW4=0.1600783D0
      WW5=0.1069393D0
      WW6=0.0471753D0
      BMU1A=XM-DDX1
      BMU1B=XM+DDX1
      BMU2A=XM-DDX2
      BMU2B=XM+DDX2
      BMU3A=XM-DDX3
      BMU3B=XM+DDX3
      BMU4A=XM-DDX4
      BMU4B=XM+DDX4
      BMU5A=XM-DDX5
      BMU5B=XM+DDX5
      BMU6A=XM-DDX6
      BMU6B=XM+DDX6
      X1A=1.0D0/DSQRT(BCB2-(BMU1A-ACB)**2)
      X1B=1.0D0/DSQRT(BCB2-(BMU1B-ACB)**2)
      X2A=1.0D0/DSQRT(BCB2-(BMU2A-ACB)**2)
      X2B=1.0D0/DSQRT(BCB2-(BMU2B-ACB)**2)
      X3A=1.0D0/DSQRT(BCB2-(BMU3A-ACB)**2)
      X3B=1.0D0/DSQRT(BCB2-(BMU3B-ACB)**2)
      X4A=1.0D0/DSQRT(BCB2-(BMU4A-ACB)**2)
      X4B=1.0D0/DSQRT(BCB2-(BMU4B-ACB)**2)
      X5A=1.0D0/DSQRT(BCB2-(BMU5A-ACB)**2)
      X5B=1.0D0/DSQRT(BCB2-(BMU5B-ACB)**2)
      X6A=1.0D0/DSQRT(BCB2-(BMU6A-ACB)**2)
      X6B=1.0D0/DSQRT(BCB2-(BMU6B-ACB)**2)
      PII=3.14159265D0
      IF(ITOP.EQ.1) TAUST=BO3(IM)
      IF(ITOP.EQ.0) TAUST=BO3(1)
      DO 30 I=1,IM
      P1=PRE(I)
      FAC1A=HTMU(BMU1A,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X1A
      FAC1B=HTMU(BMU1B,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X1B
      FAC2A=HTMU(BMU2A,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X2A
      FAC2B=HTMU(BMU2B,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X2B
      FAC3A=HTMU(BMU3A,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X3A
      FAC3B=HTMU(BMU3B,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X3B
      FAC4A=HTMU(BMU4A,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X4A
      FAC4B=HTMU(BMU4B,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X4B
      FAC5A=HTMU(BMU5A,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X5A
      FAC5B=HTMU(BMU5B,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X5B
      FAC6A=HTMU(BMU6A,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X6A
      FAC6B=HTMU(BMU6B,O2(I),O3(I),BO2(I),BO3(I),OME0
     &     ,TAUST,P1,f107,theta27)*X6B
      FACTT=(FAC1A+FAC1B)*WW1+(FAC2A+FAC2B)*WW2
     &     +(FAC3A+FAC3B)*WW3+(FAC4A+FAC4B)*WW4
     &     +(FAC5A+FAC5B)*WW5+(FAC6A+FAC6B)*WW6
      QH(I)=1.04D0*XR*FACTT/PII  ! 1.04=correction for numerical integration
      RHOI=PRE(I)/(287.0D0*TEM(I))
      QH(I)=86400.0D0*QH(I)/(1004.0D0*RHOI)
  30  CONTINUE
      RETURN
      END
      
C  ============== local O2 and O3 heating ========================

      SUBROUTINE HEAT3(QH,TIME,DEC,PHI,BMU,MODE,O2,O3,BO2,BO3
     &             ,PRE,TEM,IM,ITOP,OME0,f107,theta27)
CCCCC
C       To calculate the heating rate (QH: K/day) at TIME (hour) for O3 
C       and O2, following the parameterization scheme by Strobel (1978).
C       DEC=declination (degree), PHI=latitude (degree).
C       O2,O3=number density (m^-3), BO2,BO3=column density (m^-2).
C       PRE=pressure (Pa), TEM=temperature. FCLEAR=Fraction of clear sky.
C       ITOP=1 => PRE(1)=TOP BOUNDARY;  ITOP=0 => PRE(IM)=TOP BOUNDARY
C       MODE=1 => input=TIME,DEC,PHI; MODE=0 => input=BMU
CCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION QH(IM),O2(IM),O3(IM),BO2(IM),BO3(IM),PRE(IM)
     &     ,TEM(IM)
      DO 10 I=1,IM
  10  QH(I)=0.0D0
      IF(MODE.EQ.1) THEN
      RDEC=DEC/57.2958D0
      RPHI=PHI/57.2958D0
      ACB=DSIN(RDEC)*DSIN(RPHI)
      BCB=DCOS(RDEC)*DCOS(RPHI)
      HHX=(TIME-12.0D0)*0.261799D0             ! PI/12=0.261799
      BMU1=ACB+BCB*DCOS(HHX)
      ENDIF
      IF(MODE.EQ.0) BMU1=BMU
      IF(BMU1.LE.0.0D0) RETURN          !    The sun does not rise
      IF(ITOP.EQ.1) TAUST=BO3(IM)
      IF(ITOP.EQ.0) TAUST=BO3(1)
      DO 30 I=1,IM
      P1=PRE(I)
      QH(I)=HTMU(BMU1,O2(I),O3(I),BO2(I),BO3(I)
     &          ,OME0,TAUST,P1,f107,theta27)
      RHOI=PRE(I)/(287.0D0*TEM(I))
      QH(I)=86400.0D0*QH(I)/(1004.0D0*RHOI)
  30  CONTINUE
      RETURN
      END
      
      FUNCTION HTMU(BMU,O2,O3,BO2,BO3,OME0,BO3T,P1,f107,theta27)
CCCCC
C       Heating rate (MKS units) for O3 and O2, following the
C       parameterization scheme by Strobel (1978). BMU=COS(THETA)
C       O2,O3=number density (m^-3), BO2,BO3=column density (m^-2). 
C       Six bands: Chappius (7500-4500),Huggins(3550-2825),
C       Hadley(2725-2450),Herzberg continuum(2400-2060), 
C       Schumann-Runge continuum(1740-1260),Schumann-Runge(2025-1750).
CCCCC
      implicit real*8(a-h,o-z)
      XRM=35.0D0/DSQRT(1224.0D0*BMU**2+1.0D0)    ! magnification factor
      SNO3=BO3*XRM
      SNO2=BO2*XRM
C...OZONE ABSORPTION
c      HART=HT12(5.13D0,8.70D-22,SNO3,1.0D0,1)*EFFCY(P1,3)   !   Hartley band
c      HUGG2=HT12(-2.0D-2,2.470D-23,SNO3,0.01273D0,2)        !   Huggins band: long
c      HUGG3=HT12(-5.0D-2,3.573D-22,SNO3,0.01273D0,2)        !   Huggins band: short
c      HUBB=5.49882D0/SNO3+HUGG2+HUGG3
ccc inclusion of solar uv radiation variability:
c
      dsint5=0.5d0*dsin(theta27)
c
      chi1=1.d0+.02981d0*(f107-1.45d2)/1.5d2            !   Hartley band
      eps1=1.d0+(0.00098d0+1.400d-5*(f107-7.0d1))*dsint5
      xirrad1=chi1*eps1*5.13d0
      HART=HT12(xirrad1,8.70D-22,SNO3,1.0D0,1)*EFFCY(P1,3)
c
      chi1=1.d0+.00237d0*(f107-1.45d2)/1.5d2            !   Huggins band: long
      eps1=1.d0+(0.00000d0-1.723d-6*(f107-7.0d1))*dsint5
      xfact11=chi1*eps1
      xirrad1=xfact11*0.07d0
      chi2=1.d0+.01530d0*(f107-1.45d2)/1.5d2            !   Huggins band: short
      eps2=1.d0+(0.00062d0+1.720d-6*(f107-7.0d1))*dsint5
      xirrad2=chi2*eps2*0.05d0
      difi21=xirrad2-xirrad1
      HUGG2=HT12(difi21,2.470D-23,SNO3,0.01273D0,2)        !   Huggins band: long
      HUGG3=HT12(-xirrad2,3.573D-22,SNO3,0.01273D0,2)      !   Huggins band: short
      HUBB=5.49882D0*xfact11/SNO3+HUGG2+HUGG3
c
      chi1=1.d0+.04333d0*(f107-1.45d2)/1.5d2                !   Herzberg continuum
      eps1=1.d0+(0.00320d0+3.748d-5*(f107-7.0d1))*dsint5
      xirrad1=chi1*eps1*1.5d0
c
      chi5=1.d0+.0008d0*(f107-1.45d2)/1.5d2  ! Chappius band AND diffuse scattered radiation
      xirrad5=chi5*370.0d0
c
c      CHAP=HT12(370.0D0,2.85D-25,SNO3,1.0D0,1)         !   Chappius band
      CHAP=HT12(xirrad5,2.85D-25,SNO3,1.0D0,1)          !   Chappius band

      HTSUM1=O3*(CHAP+HART+HUBB)
      ATT=DEXP(-DMIN1(2.0D2,5.0D-28*SNO2+4.0D-22*SNO3))
c      HTSUM2=1.5D0*(5.0D-28*O2+4.0D-22*O3)*ATT              !   Herzberg continuum
      HTSUM2=xirrad1*(5.0D-28*O2+4.0D-22*O3)*ATT             !   Herzberg continuum

c
C...OXYGEN ABSORPTION
c      SRC1=HT12(0.65D-3,1.1D-21,SNO2,1.0D0,1)      !   Schumann-Runge continuum 1
c      SRC21=HT12(1.23D-3,3.0D-23,SNO2,1.0D0,2)     !   Schumann-Runge continuum 2
c      SRC22=HT12(-8.3D-4,2.0D-22,SNO2,1.0D0,2)     !   Schumann-Runge continuum 2
c      SRC23=HT12(-4.0D-4,1.5D-21,SNO2,1.0D0,2)     !   Schumann-Runge continuum 2
c      HTSRC=O2*(SRC1+SRC21+SRC22+SRC23)*EFFCY(P1,2)
      HTSRC=O2*xsrcq(SNO2,f107)*EFFCY(P1,2)
c
      eps1=1.d0+(0.01788d0+1.619d-4*(f107-7.0d1))*dsint5 !   Schumann-Runge continuum
      HTSRC=HTSRC*eps1   ! correction for 27-day period of uv variability only
c
      FACT1=DSQRT(1.0D0+1.734D-22*SNO2)        ! 4*sigma/(pi*y) = 4*2.07E-24/(pi*0.0152)
      FACT2=0.023876D0*(FACT1-1.0D0)           ! pi*y/2 = (pi*0.0152)/2 = 0.023876
      HTSRB=2.6496D-26*DEXP(-DMIN1(2.0D2,FACT2))/FACT1 ! F*sigma=0.0128*2.07E-24=2.6496E-26
      HTSRB=HTSRB*O2                  !   Schumann-Runge band
c
      chi1=1.d0+.08589d0*(f107-1.45d2)/1.5d2            !   Schumann-Runge band
      eps1=1.d0+(0.01124d0+6.589d-5*(f107-7.0d1))*dsint5
      HTSRB=HTSRB*chi1*eps1
c
C  Heating rate due to diffuse scattered solar radiation
      TAU=2.85D-25*BO3
      TAUS=2.85D-25*BO3T
      FAC1=TAUS/BMU+1.9D0*(TAUS-TAU)
      HTSCAT=2.109D-22*O3*OME0*BMU*DEXP(-DMIN1(2.0D2,FAC1))   ! 2.11E-22= 2*370*2.85E-25
      HTMU=HTSUM1+HTSUM2+HTSRC+HTSRB+HTSCAT*chi5
c  add heating rate due to Lyman-alpha line Nicolet (1985)
      q2eps1=0.49811d0*4.90246d-3*(1.0d0+0.2d-2*(f107-65.0d0))  ! 0.49811d0=heating efficient factor
c
      eps1=1.d0+(0.125d0+8.333d-4*(f107-7.0d1))*dsint5 !   Lyman-alpha
      q2eps1=q2eps1*eps1
c
      sno2cm=SNO2*1.0d-4    ! O2 column number density in cm^-2
      sigo2m=1.8728d-22/sno2cm**0.1145  ! O2 cross section in m^2: 1.8728=2.115*0.8855d0
      att=dexp(-DMIN1(2.0D2,2.115d-18*sno2cm**0.8855))
      HTMU=HTMU+q2eps1*sigo2m*o2*att  ! heating rate in J m^-3 s^-1
      RETURN
      END
      
      FUNCTION HT12(F1,SIG,SBN,SM,ICON)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(ICON.EQ.1) HT12=F1*SIG*DEXP(-DMIN1(2.0D2,SIG*SBN))
      IF(ICON.EQ.2) HT12=F1*DEXP(-DMIN1(2.0D2,SIG*SBN))/(SBN*SM)
      RETURN
      END
      
      FUNCTION EFFCY(P,MODE)
C  To calculate the heating efficiency for O3 Hartley band (MODE=3)
C  and O2 SR continuum (MODE=2). P = pressure in pascal. 
      IMPLICIT REAL*8 (A-H,O-Z)
      PMB=P*1.0D-2      ! pressure in mb
      IF(PMB.GT.1.0D0) PMB=1.0D0
      IF(PMB.LT.1.0D-4) PMB=1.0D-4
      XP=DLOG10(PMB)
      IF(MODE.EQ.3) THEN
        IF(PMB.LE.1.0E-2) THEN
	X=XP+3.0D0
	C0=0.66965D0
	C1=-0.009682D0
	C2=0.033093D0
	C3=0.017938D0
	ELSE
	X=XP+1.0D0
	C0=0.92621D0
	C1=0.13396D0
	C2=-0.076863D0
	C3=0.006897D0
	ENDIF
      ENDIF
      IF(MODE.EQ.2) THEN
        IF(PMB.LE.1.0D-2) THEN
	X=XP+3.0D0
	C0=0.75349D0
	C1=0.0036D0
	C2=0.059468D0
	C3=-0.022795D0
	ELSE
	X=XP+1.0D0
	C0=0.92621D0
	C1=0.13396D0
	C2=-0.076863D0
	C3=0.006897D0
	ENDIF
      ENDIF
      X2=X*X
      X3=X2*X
      EFFCY1=C0+C1*X+C2*X2+C3*X3
      IF(MODE.EQ.3) EFFCY=1.0D0-0.78831D0*(1.0D0-EFFCY1)
      IF(MODE.EQ.2) EFFCY=1.0D0-0.41D0*(1.0D0-EFFCY1)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  xsrcq  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function xsrcq(xn,f107)
c  net radiative heating rate (J s^-1) of Schumann-Runge continuum as a 
c  function of F107 and slant column number density in [molecules m^-2]
c  parameters determined from UARS/SOLSTICE solar flux
c  and Strobel's (1978) compilation of sigma(SRC).  02/08/1999
      implicit real*8 (a-h,o-z)
      f107sq=f107*f107
      uu0=6.794d-1+2.248d-4*f107-4.627d-7*f107sq
      uu1=2.281d-1-1.720d-4*f107+3.913d-7*f107sq
      uu2=8.253d-2-4.218d-5*f107+5.391d-8*f107sq
      uu3=3.569d-3-9.442d-6*f107+2.366d-8*f107sq
      xqinf=9.857d-25+2.782d-27*f107-4.503d-30*f107sq
      a0n=dmin1(1.094d-21*xn,200.0d0)
      a1n=dmin1(2.735d-22*xn,200.0d0)
      a2n=dmin1(6.8375d-23*xn,200.0d0)
      a3n=dmin1(1.709375d-23*xn,200.0d0)
      xsrcq=xqinf*(uu0*dexp(-a0n)+uu1*dexp(-a1n)
     &  +uu2*dexp(-a2n)+uu3*dexp(-a3n))
      return
      end
