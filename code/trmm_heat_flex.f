      SUBROUTINE TRMM_HEAT ( DPGQSP,TIM,WATERCOEFF,OZONECOEFF,TROPINDEX )
C     this is the same as trmm_hspec_cycle, except subroutine name has changed
c     so it can be called from sigprim_3heat
C     interpolates to alt grid
c     Zonal mean heating not used - turn on if you want a Hadley cell

	IMPLICIT NONE

	INCLUDE 'mcons.inc'
	INCLUDE 'spcons.inc'
	INCLUDE 'mgrid.inc'
	INCLUDE 'spgrid.inc'

	REAL TIM
	COMPLEX DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)

	INTEGER IER, WLUN

	INTEGER I, J, K, L, M, N, ND, ND0
	LOGICAL FIRST /.TRUE./

	INTEGER NALTP, NTIMP
	PARAMETER(NALTP=42, NTIMP=20000)
	COMPLEX CTIMW(NTIMP), CALTW(NALTP)
	INTEGER NTRUN, NTIMW, NALTW
	COMPLEX XTRUN, XNTIMW, XNALTW

	COMPLEX HEATIN0(0:M1MAX,0:N1MAX,NALTP)
	COMPLEX HEATIN1(0:M1MAX,0:N1MAX,NALTP)
c	COMPLEX HEAT(0:M1MAX,0:N1MAX,L1MAX)

	REAL DAY, DT, WT0, WT1, day0, CRANG, DELT

        character*80 arg
	INTEGER IARG, BRECL
	common /opens/ iarg

	REAL RAMP, Z(L1MAX), WT(L1MAX)
	COMPLEX CO1, CO2
	INTEGER NALT, LMX, INDX(L1MAX)
	real time0
c	common /htim/ time0	!time a sequence of continuation runs started
	common /htim/ time0	!offset from tim so forcing can be shifted

        REAL MULTIPLIER(L1MAX)
        REAL WATERCOEFF, OZONECOEFF
        INTEGER TROPINDEX

C     Input and output folders and files
      CHARACTER*50 INFOLD, INFILE, OUTFOLD, OUTFILE, TIDEFILE, TTIDEFILE
      CHARACTER*50 INITTEMP, INITZONAL, INITGRID, INITZGRID, GWFDUMP, CURT
      CHARACTER*50 SURFFILE
      COMMON /FILES/ INFOLD, INFILE, SURFFILE, OUTFOLD, OUTFILE,
     1   INITTEMP, INITZONAL, INITGRID, GWFDUMP, TIDEFILE, TTIDEFILE,CURT

C------------------------------------------------------------------

	IF (FIRST) THEN
	  WLUN=224
	  CALL GETARG(IARG,ARG)
	  BRECL=32
	  OPEN(WLUN, FILE=TRIM(INFOLD)//TRIM(TIDEFILE), STATUS='OLD',SHARED,
	1      FORM='UNFORMATTED', ACCESS='DIRECT', RECL=BRECL)
	  IARG=IARG+1

c	  open(unit=23,file='d:/model/circ/output/dump.dat')

	  READ (WLUN,REC=1) XTRUN, XNALTW, XNTIMW
	  WRITE(6,*) ' HEAT DIMS', XTRUN, XNALTW, XNTIMW
	  CLOSE(WLUN)
	  NTRUN=XTRUN
	  NTIMW=XNTIMW
	  NALTW=XNALTW

	  IF (NALTW .GT. NALTP) then
	    WRITE(6,*) 'HEATING ARRAY SIZE TOO BIG'
	    STOP
	  ENDIF

	  BRECL=(NTRUN+1)*(NTRUN+1)*NALTW*8
	  OPEN(WLUN,  FILE=TRIM(INFOLD)//TRIM(TIDEFILE), STATUS='OLD',SHARED,
	1      FORM='UNFORMATTED', ACCESS='DIRECT', RECL=BRECL)


	  READ (WLUN,REC=2)
	1      (CALTW(I),I=1,NALTW), (CTIMW(I),I=1,NTIMW)

	  CRANG=CTIMW(NTIMW)-CTIMW(1)
	  ND0=0

	  DO L=1,L1
	    Z(L)=-7.*ALOG(PLV(L1-L+1)*PSURF/1.E5)
	  ENDDO

	  LMX=0
	  DO L=1,L1
	    IF (Z(L) .LT. REAL(CALTW(NALTW))) LMX=L
	  ENDDO

	  DO I=1,LMX
	    INDX(I)=1
	    DO J=1,NALTW-1
	      IF (REAL(CALTW(J)) .LT. Z(I)) INDX(I)=J
	    ENDDO
	  ENDDO

	  DO K=1,LMX
	    WT(K)=(Z(K)-CALTW(INDX(K))) / (CALTW(INDX(K)+1)-CALTW(INDX(K)))
	  ENDDO

          WRITE(6,*) 'Heating bottom, tropopause and top height (km): '
          WRITE(6,*) Z(1), Z(TROPINDEX), Z(LMX)

	  FIRST=.FALSE.
c	  time0=tim
	ENDIF

c	RAMP=.5*(1+TANH(((TIM-time0)/86400.-2)/2.))
	RAMP=1.

C...    INTERPOLATE TO CURRENT TIME
C...cycle periodically over the time range of the data
c	DAY=(TIM)/86400.
	DAY=(TIM+time0)/86400.
	DELT=MOD(DAY-REAL(CTIMW(1)),CRANG)
	DO WHILE (DELT .LT. 0)
	  DELT=DELT+CRANG
	ENDDO
	DAY=CTIMW(1)+DELT

	ND=1
	DO J=2,NTIMW-1
	  IF (DAY .GT. REAL(CTIMW(J))) ND=J
	ENDDO
	IF (ND .NE. ND0) THEN
	  READ (WLUN,REC=2+ND) (((HEATIN0(I,J,L), I=0, NTRUN), J=0, NTRUN), L=1, NALTW)
	  READ (WLUN,REC=3+ND) (((HEATIN1(I,J,L), I=0, NTRUN), J=0, NTRUN), L=1, NALTW)

	  WRITE(6,*) 'HTIM:',ND,DAY,CTIMW(ND)
	  ND0=ND
	ENDIF

	DELT=CTIMW(ND+1)-CTIMW(ND)
	WT0=(CTIMW(ND+1)-DAY)/DELT*ramp
	WT1=(DAY-CTIMW(ND))/DELT*ramp

C     heating multiplier as a function of model level

        DO K=1,L1MAX
          MULTIPLIER(K)=WATERCOEFF
          IF (K .GT. TROPINDEX) THEN
            MULTIPLIER(K)=OZONECOEFF
          ENDIF
        ENDDO
      WRITE(6,*) WT0, WT1
	DO L=1,LMX
	  K=L1+1-L
          DO N=0,N1

c...remove comment to include zonal mean heating
c            DO M=0,MIN0(N,M1)
            DO M=1,MIN0(N,M1)
	      CO1=WT0*HEATIN0(M,N,INDX(L))+WT1*HEATIN1(M,N,INDX(L))
	      CO2=WT0*HEATIN0(M,N,INDX(L)+1)+WT1*HEATIN1(M,N,INDX(L)+1)

	      DPGQSP(M,N,JPOT,K)=DPGQSP(M,N,JPOT,K) + MULTIPLIER(L)*((1.-WT(L))*CO1 + WT(L)*CO2)
c	      if (nd .ne. nd0) heat(m,n,k)=(1.-WT(L))*CO1 + WT(L)*CO2
c	      heat(m,n,k)=(1.-WT(L))*CO1 + WT(L)*CO2
            END DO
          END DO
        END DO

c	write(23,*) (((heat(m,n,k),m=0,m1),n=0,n1),k=1,l1)
c	if (nd .ne. nd0) then
c	  nd0=nd
c	  write(23,*) (((heat(m,n,k),m=0,m1),n=0,n1),k=1,l1)
c	endif

	END
