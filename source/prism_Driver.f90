!      ===========================================================
      program prism
        !     ===========================================================
        !SIGMA-COORDINATE PRIMITIVE EQUATION MODEL DRIVER
        
              use Timing_mod
              use Declarations_mod
              use Prognos_sigma_mod
              use Config_mod,           only: Config_Constants
              use IO_mod,               only: INIPRE_READ
        
              IMPLICIT NONE
        
              !INCLUDE 'mcons_f90.inc'
              !INCLUDE 'spcons_f90.inc'
              !INCLUDE 'mgrid_f90.inc'
              !INCLUDE 'spgrid_f90.inc'
              !INCLUDE 'tmarch_f90.inc'
              !INCLUDE 'heating_f90.inc'
        
              REAL A, OMEGA, R, CPRESS, GRAVIT
              REAL*4 PTHICK(L1MAX)
              REAL TZMEAN(L1MAX), TRADEQ(K2MAX,L1MAX)
              REAL TSTEEP, DEL8DF, SFRIC, VERVIS(L1MAX), DIFFC(L1MAX)
              REAL ALPHA0, ALPHA1, TTOP, WTOP, SMAX, STEMP, SMIN
              REAL AMPBOT, EDDIF, HBOT, WBOT
              REAL TDAMP(L1MAX), UDAMP(L1MAX), DAMP8, RDAMP, NFAC
              REAL WAMP, WLEV, WWID
              REAL Teddy, Ueddy
              REAL RHO, MOL, DURP
        
              REAL*4 PHALFL(0:L1MAX)
              REAL*4 PP(L1MAX+1)
              REAL DAYINV, PBOT, PLID, IMPLCT, DT1, ROBFAC, TIM
              REAL ZBOT, ZLID
              REAL DZ0, DELZ, ZLO, WIDTHREF
              REAL ANDREWS, Z1, Z2, ETA0, A1, A2
              
              COMPLEX  PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
              COMPLEX  PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
              COMPLEX  PSSP0(0:M1MAX,0:N1MAX)
              COMPLEX  PSSP1(0:M1MAX,0:N1MAX)
              INTEGER MMAX, MMAX0, NMAX, NMAX0, NLEV, I,  L, M, N, IPGQ, IT
              REAL DUR, SKIP
              INTEGER NSTEP, NSKIP, ZLUN
              INTEGER TRMM, TIDE, ERATIDE
        
              REAL HCOOL, CMAM_KZZ
              EXTERNAL  HCOOL, CMAM_KZZ
        
              REAL Z
              INTEGER BRECL, ZLVLS
        
              INTEGER ihr1, ihr2, imin1, imin2, isec1, isec2, ihs1, ihs2
              REAL cpt1, cpt2
        
              INTEGER IARG
              COMMON /opens/ iarg
              
              REAL MEAN, SIG, LAM, TEMP1
              integer fsf,cff
              
              REAL TIME0
              COMMON /htim/ TIME0
              REAL WATERCOEFF, OZONECOEFF, MESOCOEFF
              INTEGER TROPINDEX, STRATINDEX
              REAL BFAC, ROBRAMP, ROBINIT
              
              REAL ZS(L1MAX)
        
        !   Damping and clamping coefficient profiles
              REAL WDAMP(L1MAX), MDAMP(L1MAX), TDIFF(L1MAX)
              REAL ALPHAT(L1MAX), ALPHAW(L1MAX), RAYFAC
              INTEGER SDTRUNC, DAMPTRUNC, RECLEN
              COMMON /WAVD/ WDAMP, MDAMP, TDIFF, ALPHAT, ALPHAW, SDTRUNC, DAMPTRUNC
        
        !   Input and output folders and files
              CHARACTER*50 INFOLD, INFILE, OUTFOLD, OUTFILE, TIDEFILE, TTIDEFILE
              CHARACTER*50 INITTEMP, INITZONAL, INITGRID, INITZGRID, GWFDUMP, CURT
              CHARACTER*50 SURFFILE
              COMMON /FILES/ INFOLD, INFILE, SURFFILE, OUTFOLD, OUTFILE, &
                             INITTEMP, INITZONAL, INITGRID, GWFDUMP, TIDEFILE, TTIDEFILE,CURT
             
              REAL WCLAMP1, WCLAMP2, WCLAMP3
              INTEGER TROPLVL, STRATLVL
              COMMON /WAVECLAMP/ WCLAMP1, WCLAMP2, WCLAMP3, TROPLVL, STRATLVL   
             
          !================================================================
          !read input from file
          !================================================================
              
              READ(5,*) INFOLD         ! path to input file
              READ(5,*) OUTFOLD        ! path to output file
              READ(5,*) INFILE         ! specified dynamics input file name
              READ(5,*) SURFFILE       ! surface geopotential forcing
              READ(5,*) OUTFILE        ! output file name
              READ(5,*) TIDEFILE       ! tide heating file
              READ(5,*) TTIDEFILE      ! tide heating file above sponge layer 
              READ(5,*) INITTEMP       ! initialized temperature field
              READ(5,*) INITZONAL      ! initialized zonal wind field
              READ(5,*) INITGRID       ! initialized model vertical pressure grid lvls
              READ(5,*) INITZGRID      ! initialized model vertical geom. alt grid lvls
              READ(5,*) GWFDUMP        ! dump file used in GWF routine
              READ(5,*) CURT           ! file containing curtis matrices
              
              WRITE(6,*) 'WIND/TEMP:    ', INFILE
              WRITE(6,*) 'SURF GEOP:    ', SURFFILE
              WRITE(6,*) 'TIDE FORC:    ', TIDEFILE
              WRITE(6,*) 'OUTPUT:       ', OUTFILE
              
              CALL CPU_TIME(cpt1)
              CALL GETTIM(ihr1,imin1,isec1,ihs1)
        
              CALL INPUTC ! call INPUTC once to read Curtis matrices at beginning
        
              write(6,*) 'use final save file? continue from file? (1/0=y/n)'
              read(5,*) fsf,cff
              write(6,*) fsf,cff
        
              NTRACE=0 ! no tracers. Tracer advection scheme does not work.
              iarg=1
        
        !     DZ0 and DELZ can be the same.  Set Dz0=0 for NCEP grid up to 30 km
              WRITE(6,*) 'ZBOT, ZTOP, DZ0, DELZ, Zlo'
              READ(5,*)  ZBOT, ZLID, DZ0, DELZ, Zlo
              WRITE(6,*)  ZBOT, ZLID, DZ0, DELZ, Zlo
        
              PBOT= 101325.00 ! era5 bottom model half level
              PLID= 1e-5 - 1e-6 ! ERAH: 5e-5 - 5e-6, ERAh : 1e-6 - 1e-7! ERAs: 1e-5 - 1e-6, ERAhhh: 2.5e-5 - 2.5e-6; ERAn: 1e-5 - 1e-6
              
              WRITE(6,*) 'TS (Hrs), dur (days), write interval (hrs), T0 (DOY):'
              READ(5,*) DT1, DUR, SKIP, TIM
              WRITE(6,*) DT1, DUR, SKIP, TIM
              NSTEP = IFIX(DUR * 24. / DT1) ! integer
              NSKIP = IFIX(SKIP / DT1)
              DT1= DT1 * 86400. / 24.
              TIM = TIM * 86400.
              TIME0 = TIM
        
              WRITE(6,*) 'HORIZONTAL TRUNCATION: MMAX, NMAX:'
              READ(5,*)  MMAX, NMAX, RAYFAC
              WRITE(6,*)  MMAX, NMAX
              NMAX = MIN(NMAX,N1MAX-1)
        
              WRITE(6,*) 'TROPOSPHERE STRATOSPHERE INDS HEIGHT INDEX: '
              READ(5,*) WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX
              write(6,*) WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX
        
              WRITE(6,*) 'DAMPING AT SMALLEST SCALE, ION drag 115 km (DAY-1)'
              WRITE(6,*)  'NEWT COOLING FACTOR,  & SURFACE FRICTION'
              READ(5,*)  DAMP8, RDAMP, NFAC, SFRIC
              WRITE(6,*)  DAMP8, RDAMP, NFAC, SFRIC
        
              WRITE(6,*) 'Lower and upper sponge params: '
              READ(5,*) ALPHA0, TTOP, ALPHA1, WTOP, WAMP, WLEV, WWID, SMAX, SMIN, BFAC, WIDTHREF
              WRITE(6,*) ALPHA0, TTOP, ALPHA1, WTOP, WAMP, WLEV, WWID, SMAX, SMIN, BFAC, WIDTHREF
              
              WRITE(6,*) 'WAVE FIELD CLAMP PARAMS:'
              READ(5,*) WCLAMP1, WCLAMP2, WCLAMP3, TROPLVL, STRATLVL
              WRITE(6,*) WCLAMP1, WCLAMP2, WCLAMP3, TROPLVL, STRATLVL
             
              WRITE(6,*) 'Bottom layer sponge params: '
              READ(5,*) AMPBOT, EDDIF, HBOT, WBOT
              WRITE(6,*) AMPBOT, EDDIF, HBOT, WBOT
        
              WRITE(6,*) 'Eddy diff factors for T/U (applied to CMAM profile)'
              READ(5,*) Teddy, Ueddy
              write(6,*) Teddy, Ueddy
        
              WRITE(6,*) 'HEAT FLAGS - RAD: 1(OFF) 2(ZSYM) 3(LOCAL), TRMM, TIDE'
              WRITE(6,*) '0/1=OFF/ON'
              READ(5,*)  MDLRD, TRMM, TIDE, ERATIDE
              WRITE(6,*) MDLRD, TRMM, TIDE, ERATIDE
        
        !     Initialize physical parameters
              DAYINV=1.0 / 86400.
              OMEGA=2.0 * 3.14159265359 * DAYINV ! angular velocity Earth
              R=287.058 ! gas constant
              CPRESS= 1013.2500 ! era5 model lvl global mean surface p
              GRAVIT=9.80665 ! gravitational acceleration m2 / s2
              A=6.3781E6 ! radius earth
              
              CALL PLINI( OMEGA, R, CPRESS, GRAVIT )
        
        !   Initialize horizontal truncation
              CALL SPINI( MMAX, NMAX, A )
              
              CALL Config_Constants()

        !   Initialize pressure-level values and compute pressure level thicknesses
              IF (DZ0 .EQ. 0.) THEN
                CALL INIPRE_READ(PTHICK, PP, PHALFL, PBOT, PLID, NLEV)
              ENDIF
              
        ! !   read in exact (hydrostatic) geometric altitude lvls 
        !       ZLUN = 27
        !       OPEN(ZLUN,FILE=TRIM(INFOLD)//TRIM(INITZGRID),STATUS='OLD',
        !      1    SHARED,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=100)
        
        !       READ (ZLUN,REC=1) RECLEN,ZLVLS
        !       CLOSE(ZLUN)
                         
        !       BRECL=ZLVLS*4 + 8
                    
        !       OPEN(ZLUN,FILE=TRIM(INFOLD)//TRIM(INITZGRID),STATUS='OLD',
        !      1    SHARED,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=BRECL)
        
        !       READ (ZLUN,REC=2) RECLEN,(ZS(L),L=1,ZLVLS)
        !       CLOSE(ZLUN)
              
        ! !   Rayleigh surface friction and Newtonian cooling profiles
        !       DO L=1, NLEV
        !         Z=PP(L)/PBOT ! sigma-coordinates
                
        !         IF (Z .GT. 0.7) THEN
        !            UDAMP(L)= AMPBOT * DAYINV * (0.7 - Z) / (0.7 - 1.0)
        !         ENDIF
                
        ! !     standard newtonian cooling profile
        !         TDAMP(L) = NFAC * HCOOL(PP(L)) * DAYINV 
        !       ENDDO
              
        ! !   Specified dynamics truncated at MMAX if input trunc exceeds MMAX
        !       READ(5,*) SDTRUNC, DAMPTRUNC
        !       IF (SDTRUNC .GT. MMAX) THEN
        !         SDTRUNC = MMAX
        !       ENDIF
        !       WRITE(6,*) 'Specified Dynamics truncated at M =', SDTRUNC, DAMPTRUNC
                          
        ! !   Lower sponge - clamp mean winds and temps. Allow for separate weights and layer top
        ! 	  ALPHA0=ALPHA0*DAYINV
        !       ALPHA1=ALPHA1*DAYINV
        
        !       DO L=1,NLEV
        !         Z=-7.*ALOG(PP(L)/1.E5)
        !         ALPHAT(L)=ALPHA0*(1.-TANH((Z-TTOP)/3.))/2. 
        !         ALPHAW(L)=ALPHA1*(1.-TANH((Z-WTOP)/3.))/2. 
        !       ENDDO
        
        ! !   wdamp damps waves, mdamp relaxes mean vorticity to VOR0 above WLEV
        !       WRITE(6,*) WAMP, WAMP, WAMP
        !       WRITE(6,*), WWID, WWID, WWID
        !       WRITE(6,*) WLEV
              
        ! 	  WAMP=WAMP*DAYINV
        !       WAMP=2.5e-4  ! 2.5e-4 for Hong and Lindzen solar maximum conditions
        !       SMAX=SMAX*DAYINV
        !       SMIN=SMIN*DAYINV
        !       AMPBOT=AMPBOT*DAYINV
                    
        ! !   upper sponge layer
        !       DO L=1,NLEV
        !         Z = ZS(L) 
        
        !         IF (Z .GT. 110) THEN
        !           STEMP = MAX(1e-4 * TANH((Z - 110) / 40), (1.+TANH((Z-WLEV)/WWID))*WAMP/2.)
        ! !        TDAMP(L)=TDAMP(L) + (1. + TANH((Z-WLEV)/WWID))*WAMP/2.  !+ 0.05*DAYINV
        !         ELSE
        !           STEMP = 0.
        !         ENDIF
                
                                
        !         IF (STEMP .GT. SMAX) STEMP = SMAX
        !         IF (STEMP .LT. SMIN) STEMP = 0.
                
        !         WDAMP(L)=STEMP
                        
        !         ALPHAW(L)=MAX(ALPHAW(L),WDAMP(L))
        !         ALPHAT(L)=MAX(ALPHAT(L),WDAMP(L))
                
        !         TDAMP(L)=TDAMP(L) + (1. + TANH((Z-WLEV)/WWID))*WAMP/2.  !+ 0.05*DAYINV
                
        !         MDAMP(L)=0
                             
        !         WDAMP(L)=WDAMP(L) !+ 0.05*DAYINV
        !         MDAMP(L)=WDAMP(L) 
                
        !         Z1= 12.0 * 7.0
        !         Z2= 13.3 * 7.0
        !         ETA0= 2.18e-5 * 0.0
        !         A1= 4.0 * 7.0
        !         A2= 2.0 * 7.0
                
        !         IF (Z .LE. Z1) THEN
        !           ANDREWS = ETA0*EXP(-((Z - Z1)/A1)**2)
        !         ELSE
        !           IF (Z .LE. Z2) THEN
        !             ANDREWS = ETA0
        !           ELSE
        !             ANDREWS = ETA0*EXP(-((Z - Z2)/A2)**2) 
        !           ENDIF
        !         ENDIF
                
        !         UDAMP(L)= UDAMP(L) + ANDREWS
        !       ENDDO
              
        ! !    WDAMP(1) = 25 * DAYINV
        
        !       DEL8DF= DAYINV / FLOAT(NMAX*(NMAX+1))**4 * DAMP8 ! Horiz diffusion coef
        
        ! !   Vertical grid stuff (output is in common blocks defined in MGRID.INC)
        !       SFRIC = SFRIC*DAYINV
        !       CALL VERINI
        !      1  (NLEV, PTHICK, PHALFL, PP, DEL8DF, TDAMP, UDAMP, SFRIC)
        
        ! !   Compute standard atmospheric mean and radiative equilibrium temperature profiles
        !       CALL INITTEM(TZMEAN, TRADEQ)
        
        ! !   Background eddy diffusion and molecular diffusion profiles
        !       WRITE(6,*) '      z','  vervis ',' diffc ','    udamp ',' tdamp ',
        !      1    '  wdamp ','   mdamp',  '  alphaw ',  ' alphat '
                  
        !       DO L=1,NLEV
        !         Z = ZS(L) 
              
        !         RHO=1.34*EXP(-Z/7.)
        !         MOL=3.5E-7*(TZMEAN(L)**.666)/RHO
        !         MOL= MIN(MOL,500.) ! molecular diffusion; crashes if too large
        
        !         VERVIS(L)= MIN(CMAM_KZZ(Z,TZMEAN(L))*Ueddy, 200.) + MOL !*Ueddy
        !         DIFFC(L)= MIN(CMAM_KZZ(Z,TZMEAN(L))*Teddy, 200.) + MOL !*Teddy
                
        !         IF (Z .GT. 0) THEN
        !           VERVIS(L)=VERVIS(L) 
        !           DIFFC(L)=DIFFC(L) 
        !         ELSE
        !           VERVIS(L)=0
        !           DIFFC(L)=0
        !         ENDIF
                
        !         Z1= 13.1 * 7.0
        !         Z2= 14.3 * 7.0
        !         ETA0= 250 * 0
        !         A1= 1.1 * 7.0   !4.0 * 7.0
        !         A2= 1.113 * 7.0 !2.0 * 7.0
                
        !         IF (Z .LE. Z1) THEN
        !           ANDREWS = ETA0*EXP(-((Z - Z1)/A1)**2)
        !         ELSE
        !           IF (Z .LE. Z2) THEN
        !             ANDREWS = ETA0
        !           ELSE
        !             ANDREWS = ETA0*EXP(-((Z - Z2)/A2)**2) 
        !           ENDIF
        !         ENDIF
                
        !         VERVIS(L)= VERVIS(L) + ANDREWS
        !         DIFFC(L)=DIFFC(L) + ANDREWS
                
        !         WRITE(6,'(9f8.2)') z, vervis(L), diffc(L),
        !      1      udamp(l)/dayinv, tdamp(l)/dayinv, wdamp(l)/dayinv,
        !      1      mdamp(l)/dayinv, alphaw(L)/dayinv, alphat(L)/dayinv
        !       ENDDO
              
        !       TSTEEP= 1.2 ! (never used anything else)
        !       CALL TINI(TZMEAN, TSTEEP, TRADEQ, VERVIS, DIFFC)
        
        ! !   Write meta-data to Python output file
        !       OPEN(13, FILE=TRIM(OUTFOLD)//TRIM(OUTFILE),FORM='UNFORMATTED')
        !       WRITE(13) K1, K2, M1, N1, L1, PBOT, PLID, DT1*NSKIP, NSTEP/NSKIP+1
        !       WRITE(13) (PLV(I),I=1,L1)
        !       WRITE(13) (DEGLAT(I), I=1,K2)
        !       WRITE(13) (ZSTDLV(I)/1000.+ZBOT,I=1,NLEV)
        
        !       WRITE(6,*) 'ENTER IMPLCT, ROBFAC' ! IMPLCT = .5 required for fast waves
        !       READ(5,*) IMPLCT, ROBFAC
        !       WRITE(6,*) IMPLCT, ROBFAC
              
        !       IF (CFF .EQ. 1) THEN
        ! !     Read in initial state if continue from file
        !         DO 110 L=1,NLEV
        !           DO 110 IPGQ=JVOR,JPOT
        ! 	        READ(11)((PGQSP0(M,N,IPGQ,L),M=0,MIN0(N,MMAX0)),N=0,NMAX0)
        !             READ(11)((PGQSP1(M,N,IPGQ,L),M=0,MIN0(N,MMAX0)),N=0,NMAX0)
        !   110   CONTINUE
        
        !         READ(11)((PSSP0(M,N),M=0,MIN0(N,MMAX0)),N=0,NMAX0)
        !         READ(11)((PSSP1(M,N),M=0,MIN0(N,MMAX0)),N=0,NMAX0)
        
        !         CLOSE(11)
        
        ! !     Save the initial state
        !         CALL SPEWPY (PGQSP1, PSSP1, OUTFILE)
        
        !       ELSE
        ! !     Initialize temperature distribution, wind field and surface condtions
        !         CALL INITVAR(PGQSP0,PSSP0)
        
        !         DO L=1,NLEV
        ! 	      DO IPGQ=JVOR,JPOT
        !             CALL SPCOPY( PGQSP1(0,0,IPGQ,L), PGQSP0(0,0,IPGQ,L), N1)
        ! 	      ENDDO
        !         ENDDO
        
        !         CALL SPCOPY( PSSP1,PSSP0, N1)
        
        ! !     Initialize backward scheme for first time-step
        !         CALL DDTINI( 0.5 * DT1, IMPLCT )
        
        ! !     First time step
        !         CALL TIME_STEP_INIT (PGQSP0, PGQSP1, PSSP0, PSSP1, TIM, .5*DT1, 0., 
        !      1         TRMM, TIDE, ERATIDE, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX, RAYFAC, EDDIF, ZS)
        
        ! !     Save the initial state
        !         CALL SPEWPY (PGQSP0, PSSP0, OUTFILE)
        !       ENDIF
        
        ! !   Initialize leap-frog time-stepping
        !       CALL DDTINI( DT1, IMPLCT )
        ! !   ROBFAC = .01
        
        ! !   Time-march
        !       DO IT=1,NSTEP
        !         TIM=TIM+DT1
                
        ! !     ramp down robfac to avoid crashing soon after initialization  
        !         ROBINIT = 1.0
        !         ROBRAMP = ROBINIT - 
        !      1            (ROBINIT - ROBFAC)*(1 - EXP(-(TIM-TIME0)/(5*86400.)))
             
        ! !     solar declination for zhu heating
        !         DEC=ASIND(COS((-ABS(TIM/86400./30.5)+5.6)/6.*3.14159)*SIND(23.))
              
        ! !     old robfil to damp physical mode after initialization to avoid crash
        !         IF (ROBRAMP .LT. ROBFAC + 0.01) THEN
        !           CALL TIME_STEP (PGQSP0, PGQSP1, PSSP0, PSSP1, TIM, DT1, ROBRAMP, 
        !      1      TRMM, TIDE, ERATIDE, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX, RAYFAC, EDDIF, ZS)
        !         ELSE
        !           CALL TIME_STEP_INIT (PGQSP0, PGQSP1, PSSP0, PSSP1, TIM, DT1, ROBRAMP, 
        !      1      TRMM, TIDE, ERATIDE, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX, RAYFAC, EDDIF, ZS)
        !         ENDIF
        
        !         IF (MOD(IT,NSKIP).EQ.0) THEN
        !           CALL SPEWPY (PGQSP1, PSSP1, OUTFILE)
        !         ENDIF
        !       ENDDO
        
        !       CLOSE( 12 )
        
        !       IF (FSF .eq. 1) THEN
        ! !     CALL GETARG(IARG,ARG)! SAVE FOR CONTINUATION
        ! !     OPEN(UNIT=11,FILE=ARG(1:(INDEX(ARG,' ')-1)), form='UNFORMATTED')
        ! !     iarg=iarg+1
        !         WRITE(11) MMAX,NMAX,NLEV,PBOT,PLID,TIM,DT1
        
        !         DO L=1,NLEV
        ! 	      DO IPGQ=JVOR,JPOT
        !             WRITE(11)((PGQSP0(M,N,IPGQ,L),M=0,MIN0(N,MMAX))
        !      1        ,N=0,NMAX)
        !             WRITE(11)((PGQSP1(M,N,IPGQ,L),M=0,MIN0(N,MMAX))
        !      1        ,N=0,NMAX)
        !           ENDDO
        !         ENDDO
        
        !         WRITE(11)((PSSP0(M,N),M=0,MIN0(N,MMAX))
        !      1   ,N=0,NMAX)
        !         WRITE(11)((PSSP1(M,N),M=0,MIN0(N,MMAX))
        !      1   ,N=0,NMAX)
        !         CLOSE(11)
        !       ENDIF
        
        !       CALL cpu_time(cpt2)
        !       CALL gettim(ihr2,imin2,isec2,ihs2)
        
        !       IF (ihr2-ihr1 .lt. 0) THEN
        ! 	    ihr2=ihr2+24
        !         isec2=(ihr2-ihr1)*3600 + (imin2-imin1)*60 + isec2-isec1
        !         ihr1=isec2/3600
        !         imin1=mod(isec2/60,60)
        ! 	    isec1=mod(isec2,60)
        
        !         WRITE(6,"(a13,i10,a9,i3,':',i2,':',i2,' hr:min:sec')")
        !      1     'elapsed time:',isec2,' seconds,',ihr1,imin1,isec1
        !         WRITE(6,*) 'CPU time: ',cpt2-cpt1
        !         WRITE(6,*) 'program complete'
        !       ENDIF
        
        !       END
        
        ! !   ===================================================
        !       SUBROUTINE TIME_STEP (PGQSP0, PGQSP1, PSSP0, PSSP1, TIM, DT1, ROBFAC, TRMM, TIDE, ERATIDE, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX, RAYFAC,
        !      1                      EDDIF, ZS)
        ! !   ===================================================
        !       IMPLICIT NONE
        
        !       INCLUDE 'mcons.inc'
        !       INCLUDE 'spcons.inc'
        !       INCLUDE 'mgrid.inc'
        !       INCLUDE 'spgrid.inc'
        !       INCLUDE 'tmarch.inc'
        
        !       COMPLEX  PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
        !       COMPLEX  PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
        !       COMPLEX  DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
        !       COMPLEX  PSSP0(0:M1MAX,0:N1MAX)
        !       COMPLEX  PSSP1(0:M1MAX,0:N1MAX)
        !       COMPLEX DPSSP(0:M1MAX,0:N1MAX)
        !       COMPLEX FSP1(0:M1MAX,0:N1MAX)
        !       REAL TIM, DT1, ROBFAC
        !       INTEGER TRMM, TIDE, ERATIDE
              
        !       REAL EDDIF, ZS(L1MAX)
        
        !       REAL WATERCOEFF, OZONECOEFF, MESOCOEFF, RAYFAC
        !       INTEGER TROPINDEX, STRATINDEX
        !       INTEGER L, M, N, IPGQ
              
        !       REAL TEMP0(K2MAX,L1MAX), U0(K2MAX,L1MAX)
        !       COMPLEX VOR0(0:M1MAX,0:N1MAX,L1MAX)
        !       COMMON /ZONE/ TEMP0, U0, VOR0
              
        !       COMPLEX SURF1(0:100,0:100)
        !       COMMON /SURFACE/ SURF1
              
        !       REAL WCLAMP1, WCLAMP2, WCLAMP3
        !       INTEGER TROPLVL, STRATLVL
        !       COMMON /WAVECLAMP/ WCLAMP1, WCLAMP2, WCLAMP3, TROPLVL, STRATLVL   
              
        ! !   Compute leap-frog "adiabatic" time tendency
        !       CALL DDTPGQ( DPGQSP, DPSSP, PGQSP0, PGQSP1, PSSP0, PSSP1 )
        
        ! !   Include effects of forcing/damping
        !       CALL DDTFRC( DPGQSP, PGQSP0, TIM, RAYFAC, EDDIF, ZS )
        
        ! !   Include bottom boundary, mechanical & thermal, GW forcing
        !       CALL FORCING(DPGQSP, PGQSP1, TIM)
        
        ! !   Radiative heating (no effect if mdlrd=1) ! hard coded in ZHU_HEAT
        !       CALL zhu_HEAT(DPGQSP, PGQSP1, TIM)
        
        ! !   trmm latent heating
        !       IF (TRMM .EQ. 1) CALL TRMM_HEAT(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, TROPINDEX)
        
        ! !   tide heat sources
        !       IF (tide .eq. 1) CALL TIDE_heat(DPGQSP, TIM)
              
        ! !   ERA (ECMWF) tide forcing
        !       IF (ERATIDE .EQ. 1) CALL PY_HEAT(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX)
        ! !    IF (ERATIDE .EQ. 1) CALL PY_HEATWAC(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX)
        
        ! !   Make implicit corrections to tendency
        !       CALL DDTIMP( DPGQSP,DPSSP )
        
        ! !   Leap-frog to next time-step
        !       DO 220 L=1,L1
        !         DO 220 IPGQ=JVOR,JPOT
        !           DO 210 N=0,N1
        !             DO 210 M=0,MIN0(N,M1),1
        !               FSP1(M,N)= PGQSP0(M,N,IPGQ,L)
        !      1        + 2.0 * DT1 * DPGQSP(M,N,IPGQ,L)
        !   210     CONTINUE
        
        ! !       Apply Robert filter
        ! 	      IF (ROBFAC .NE. 0.) THEN
        !             IF (L .LT. 21) THEN ! 16 ERAh
        !                 CALL ROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        !      1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L),
        !      1                       FSP1, 0.01 * 50)
        !                 CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        !      1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
             
        ! !           ELSE IF (L .LT. 10) THEN ! 31 ERAh
        ! !              CALL MODROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        ! !   1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L),
        ! !   1                       FSP1, ROBFAC * 2)
        ! !              CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        ! !   1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
                     
        !              ELSE
        ! !               CALL MODROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        ! !   1                           PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L),
        ! !   1                           FSP1, ROBFAC)
        !                  CALL ROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        !      1                           PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L),
        !      1                           FSP1, ROBFAC)
        !                  CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        !      1                        PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
        !              ENDIF
        !           ENDIF
                  
        ! !   uncomment if you want to hard-force 
        ! !        IF (L .GT. -1) THEN
        ! !          DO M=0,0
        ! !            DO N=0,N1       
        ! !                IF (IPGQ .EQ. JPOT) THEN
        ! !                  FSP1(M,N)= TMNLV(M,N,L)
        ! !                ELSE
        ! !                  FSP1(M,N)= VOR0(M,N,L)
        ! !                ENDIF
        ! !                PGQSP0(M,N,JDIV,L)= COMPLEX(0,0)
        ! !                PGQSP1(M,N,JDIV,L)= COMPLEX(0,0)
        ! !            ENDDO
        ! !          ENDDO
        ! !        ENDIF
        
        ! 	      CALL SPCOPY( PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), FSP1, N1)
        !   220 CONTINUE
        
        ! !   Surface pressure tendency
        !       DO 230 N=0,N1
        !         DO 230 M=0,MIN0(N,M1),1
        !           FSP1(M,N)= PSSP0(M,N)
        !      1      + 2.0 * DT1 * DPSSP(M,N)
             
        !   230 CONTINUE
        
        ! !   Apply Robert filter
        !       IF (ROBFAC .NE. 0.) THEN
        !         CALL MODROBFIL( PSSP0, PSSP1, FSP1, ROBFAC)
        !         CALL SPCOPY( PSSP0, PSSP1, N1)
        !       ENDIF
        
        !       CALL SPCOPY( PSSP1, FSP1, N1)
        
        !       RETURN
        !       END
        
        ! !   ===================================================
        !       SUBROUTINE TIME_STEP_INIT (PGQSP0, PGQSP1, PSSP0, PSSP1, TIM, DT1, ROBFAC, TRMM, TIDE, ERATIDE, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX, RAYFAC,
        !      1                           EDDIF, ZS)
        ! !   ===================================================
        !       IMPLICIT NONE
        
        !       INCLUDE 'mcons.inc'
        !       INCLUDE 'spcons.inc'
        !       INCLUDE  'mgrid.inc'
        !       INCLUDE  'spgrid.inc'
        
        !       COMPLEX  PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
        !       COMPLEX  PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
        !       COMPLEX  DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
        !       COMPLEX  PSSP0(0:M1MAX,0:N1MAX)
        !       COMPLEX  PSSP1(0:M1MAX,0:N1MAX)
        !       COMPLEX DPSSP(0:M1MAX,0:N1MAX)
        !       COMPLEX FSP1(0:M1MAX,0:N1MAX)
        !       REAL TIM, DT1, ROBFAC
        !       INTEGER TRMM, TIDE, ERATIDE
        
        !       REAL WATERCOEFF, OZONECOEFF, MESOCOEFF, RAYFAC
        !       INTEGER TROPINDEX, STRATINDEX
        !       INTEGER L, M, N, IPGQ
              
        !       REAL EDDIF, ZS(L1MAX)
              
        !       REAL WCLAMP1, WCLAMP2, WCLAMP3
        !       INTEGER TROPLVL, STRATLVL
        !       COMMON /WAVECLAMP/ WCLAMP1, WCLAMP2, WCLAMP3, TROPLVL, STRATLVL   
              
        ! !   Compute leap-frog "adiabatic" time tendency
        !       CALL DDTPGQ( DPGQSP, DPSSP, PGQSP0, PGQSP1, PSSP0, PSSP1 )
        
        ! !   Include effects of forcing/damping
        !       CALL DDTFRC( DPGQSP, PGQSP0, TIM, RAYFAC, EDDIF, ZS )
        
        ! !   Include bottom boundary, mechanical & thermal, GW forcing
        !       CALL FORCING(DPGQSP, PGQSP1, TIM)
        
        ! !   Radiative heating (no effect if mdlrd=1) ! hard coded in ZHU_HEAT
        !       CALL zhu_HEAT(DPGQSP, PGQSP1, TIM)
        
        ! !   trmm latent heating
        !       IF (TRMM .EQ. 1) CALL TRMM_HEAT(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, TROPINDEX)
        
        ! !   tide heat sources
        !       IF (tide .eq. 1) CALL TIDE_heat(DPGQSP, TIM)
              
        ! !   ERA (ECMWF) tide forcing. SH input file on model levels required.
        !       IF (ERATIDE .EQ. 1) CALL PY_HEAT(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX)
        ! !    IF (ERATIDE .EQ. 1) CALL PY_HEATWAC(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX)
        
        ! !   Make implicit corrections to tendency
        !       CALL DDTIMP( DPGQSP,DPSSP )
        
        ! !   Leap-frog to next time-step
        !       DO 220 L=1,L1
        !         DO 220 IPGQ=JVOR,JPOT
        !           DO 210 N=0,N1
        !             DO 210 M=0,MIN0(N,M1),1
        !               FSP1(M,N)= PGQSP0(M,N,IPGQ,L)
        !      1        + 2.0 * DT1 * DPGQSP(M,N,IPGQ,L)
        !   210     CONTINUE
        
        ! !       Apply Robert filter
        ! 	      IF (ROBFAC .NE. 0.) THEN
        !              IF (L .LT. 21) THEN ! 21 in ERAhh, 11 in ERAz
        !                 CALL ROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        !      1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L),
        !      1                       FSP1, 0.01 * 50)
        !                 CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        !      1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
             
        ! !           ELSE IF (L .LT. 31) THEN
        ! !              CALL ROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        ! !   1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L),
        ! !   1                       FSP1, 0.01 * 25)
        ! !              CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        ! !   1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
        !             ELSE
        !                 CALL ROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        !      1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L),
        !      1                       FSP1,  ROBFAC)
        !                 CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
        !      1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
        !             ENDIF
        !           ENDIF
        
        ! 	      CALL SPCOPY( PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), FSP1, N1)
        !   220 CONTINUE
        
        ! !   Surface pressure tendency
        !       DO 230 N=0,N1
        !         DO 230 M=0,MIN0(N,M1),1
        !           FSP1(M,N)= PSSP0(M,N)
        !      1      + 2.0 * DT1 * DPSSP(M,N)
             
        !   230 CONTINUE
        
        ! !   Apply Robert filter
        !       IF (ROBFAC .NE. 0.) THEN
        !         CALL ROBFIL( PSSP0, PSSP1, FSP1, ROBFAC)
        !         CALL SPCOPY( PSSP0, PSSP1, N1)
        !       ENDIF
        
        !       CALL SPCOPY( PSSP1, FSP1, N1)
        
        !       RETURN
            END program prism
              
