module timestep_mod

    use Declarations_mod
    use Prognos_sigma_mod
    use heating_mod
    use forcing_mod

    contains

    !     ===================================================
      SUBROUTINE TIME_STEP_INIT (PGQSP0, PGQSP1, PSSP0, PSSP1, TIM, DT1, ROBFAC, TRMM, &
                                TIDE, ERATIDE, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, &
                                STRATINDEX, RAYFAC, EDDIF, ZS)
!     ===================================================
      IMPLICIT NONE

    !   INCLUDE 'mcons.inc'
    !   INCLUDE 'spcons.inc'
    !   INCLUDE  'mgrid.inc'
    !   INCLUDE  'spgrid.inc'

      COMPLEX  PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      COMPLEX  PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      COMPLEX  DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
      COMPLEX  PSSP0(0:M1MAX,0:N1MAX)
      COMPLEX  PSSP1(0:M1MAX,0:N1MAX)
      COMPLEX DPSSP(0:M1MAX,0:N1MAX)
      COMPLEX FSP1(0:M1MAX,0:N1MAX)
      REAL TIM, DT1, ROBFAC
      INTEGER TRMM, TIDE, ERATIDE

      REAL WATERCOEFF, OZONECOEFF, MESOCOEFF, RAYFAC
      INTEGER TROPINDEX, STRATINDEX
      INTEGER L, M, N, IPGQ
      
      REAL EDDIF 
      REAL*4 ZS(L1MAX)
      
      REAL WCLAMP1, WCLAMP2, WCLAMP3
      INTEGER TROPLVL, STRATLVL
      COMMON /WAVECLAMP/ WCLAMP1, WCLAMP2, WCLAMP3, TROPLVL, STRATLVL   
!     Compute leap-frog "adiabatic" time tendency
      CALL DDTPGQ( DPGQSP, DPSSP, PGQSP0, PGQSP1, PSSP0, PSSP1 )
!     Include effects of forcing/damping
      CALL DDTFRC( DPGQSP, PGQSP0, TIM, RAYFAC, EDDIF, ZS )
!     Include bottom boundary, mechanical & thermal, GW forcing
      CALL FORCING(DPGQSP, PGQSP1, TIM)

      ! NB: These can be re-added again later
! !     Radiative heating (no effect if mdlrd=1) ! hard coded in ZHU_HEAT
!       CALL zhu_HEAT(DPGQSP, PGQSP1, TIM)

! !     trmm latent heating
!       IF (TRMM .EQ. 1) CALL TRMM_HEAT(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, TROPINDEX)

! !     tide heat sources
!       IF (tide .eq. 1) CALL TIDE_heat(DPGQSP, TIM)
      
!     ERA (ECMWF) tide forcing. SH input file on model levels required.
     IF (ERATIDE .EQ. 1) CALL PY_HEAT(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX)
!      IF (ERATIDE .EQ. 1) CALL PY_HEATWAC(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX)

!     Make implicit corrections to tendency
      CALL DDTIMP( DPGQSP,DPSSP )

!     Leap-frog to next time-step
      DO 220 L=1,L1
        DO 220 IPGQ=JVOR,JPOT
          DO 210 N=0,N1
            DO 210 M=0,MIN0(N,M1),1
              FSP1(M,N)= PGQSP0(M,N,IPGQ,L) &
              + 2.0 * DT1 * DPGQSP(M,N,IPGQ,L)
  210     CONTINUE

!         Apply Robert filter
	      IF (ROBFAC .NE. 0.) THEN
             IF (L .LT. 21) THEN ! 21 in ERAhh, 11 in ERAz
                CALL ROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), &
                             PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), &
                             FSP1, 0.01 * 50)
                CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), &
                             PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
     
!             ELSE IF (L .LT. 31) THEN
!                CALL ROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
!     1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L),
!     1                       FSP1, 0.01 * 25)
!                CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
!     1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
            ELSE
                CALL ROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), &
                             PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), &
                             FSP1,  ROBFAC)
                CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), &
                             PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
            ENDIF
          ENDIF

	      CALL SPCOPY( PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), FSP1, N1)
  220 CONTINUE

!     Surface pressure tendency
      DO 230 N=0,N1
        DO 230 M=0,MIN0(N,M1),1
          FSP1(M,N)= PSSP0(M,N) &
            + 2.0 * DT1 * DPSSP(M,N)
     
  230 CONTINUE

!     Apply Robert filter
      IF (ROBFAC .NE. 0.) THEN
        CALL ROBFIL( PSSP0, PSSP1, FSP1, ROBFAC)
        CALL SPCOPY( PSSP0, PSSP1, N1)
      ENDIF

      CALL SPCOPY( PSSP1, FSP1, N1)

      RETURN
    END subroutine TIME_STEP_INIT


!   ===================================================
    SUBROUTINE TIME_STEP (PGQSP0, PGQSP1, PSSP0, PSSP1, TIM, DT1, ROBFAC, TRMM, TIDE,        &
                          ERATIDE, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX, &
                          RAYFAC, EDDIF, ZS)
!   ===================================================
             IMPLICIT NONE
       
            !  INCLUDE 'mcons.inc'
            !  INCLUDE 'spcons.inc'
            !  INCLUDE 'mgrid.inc'
            !  INCLUDE 'spgrid.inc'
            !  INCLUDE 'tmarch.inc'
       
             COMPLEX  PGQSP0(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
             COMPLEX  PGQSP1(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
             COMPLEX  DPGQSP(0:M1MAX,0:N1MAX,NPGQ,L1MAX)
             COMPLEX  PSSP0(0:M1MAX,0:N1MAX)
             COMPLEX  PSSP1(0:M1MAX,0:N1MAX)
             COMPLEX DPSSP(0:M1MAX,0:N1MAX)
             COMPLEX FSP1(0:M1MAX,0:N1MAX)
             REAL TIM, DT1, ROBFAC
             INTEGER TRMM, TIDE, ERATIDE
             
             REAL EDDIF
             REAL*4 ZS(L1MAX)
       
             REAL WATERCOEFF, OZONECOEFF, MESOCOEFF, RAYFAC
             INTEGER TROPINDEX, STRATINDEX
             INTEGER L, M, N, IPGQ
             
             REAL TEMP0(K2MAX,L1MAX), U0(K2MAX,L1MAX)
             COMPLEX VOR0(0:M1MAX,0:N1MAX,L1MAX)
             COMMON /ZONE/ TEMP0, U0, VOR0
             
             COMPLEX SURF1(0:100,0:100)
             COMMON /SURFACE/ SURF1
             
             REAL WCLAMP1, WCLAMP2, WCLAMP3
             INTEGER TROPLVL, STRATLVL
             COMMON /WAVECLAMP/ WCLAMP1, WCLAMP2, WCLAMP3, TROPLVL, STRATLVL   
             
       !   Compute leap-frog "adiabatic" time tendency
             CALL DDTPGQ( DPGQSP, DPSSP, PGQSP0, PGQSP1, PSSP0, PSSP1 )
       
       !   Include effects of forcing/damping
         !    CALL DDTFRC( DPGQSP, PGQSP0, TIM, RAYFAC, EDDIF, ZS )
       
       !   Include bottom boundary, mechanical & thermal, GW forcing
            ! CALL FORCING(DPGQSP, PGQSP1, TIM)
       
    !    !   Radiative heating (no effect if mdlrd=1) ! hard coded in ZHU_HEAT
    !          CALL zhu_HEAT(DPGQSP, PGQSP1, TIM)
       
    !    !   trmm latent heating
    !          IF (TRMM .EQ. 1) CALL TRMM_HEAT(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, TROPINDEX)
       
    !    !   tide heat sources
    !          IF (tide .eq. 1) CALL TIDE_heat(DPGQSP, TIM)
             
       !   ERA (ECMWF) tide forcing
            ! IF (ERATIDE .EQ. 1) CALL PY_HEAT(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX)
       !    IF (ERATIDE .EQ. 1) CALL PY_HEATWAC(DPGQSP, TIM, WATERCOEFF, OZONECOEFF, MESOCOEFF, TROPINDEX, STRATINDEX)
       
       !   Make implicit corrections to tendency
             CALL DDTIMP( DPGQSP,DPSSP )
       
       !   Leap-frog to next time-step
             DO 220 L=1,L1
               DO 220 IPGQ=JVOR,JPOT
                 DO 210 N=0,N1
                   DO 210 M=0,MIN0(N,M1),1
                     FSP1(M,N)= PGQSP0(M,N,IPGQ,L) &
                     + 2.0 * DT1 * DPGQSP(M,N,IPGQ,L)
         210     CONTINUE
       
       !       Apply Robert filter
                 IF (ROBFAC .NE. 0.) THEN
                   IF (L .LT. 21) THEN ! 16 ERAh
                       CALL ROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), &
                                    PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), &
                                    FSP1, 0.01 * 50)
                       CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), &
                                    PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
            
       !           ELSE IF (L .LT. 10) THEN ! 31 ERAh
       !              CALL MODROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
       !   1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L),
       !   1                       FSP1, ROBFAC * 2)
       !              CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
       !   1                       PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
                    
                    ELSE
       !               CALL MODROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), 
       !   1                           PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L),
       !   1                           FSP1, ROBFAC)
                        CALL ROBFIL( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L),  & 
                                        PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), &
                                        FSP1, ROBFAC)
                        CALL SPCOPY( PGQSP0(0:M1MAX,0:N1MAX,IPGQ,L), &
                                     PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), N1)
                    ENDIF
                 ENDIF
                 
       !   uncomment if you want to hard-force 
       !        IF (L .GT. -1) THEN
       !          DO M=0,0
       !            DO N=0,N1       
       !                IF (IPGQ .EQ. JPOT) THEN
       !                  FSP1(M,N)= TMNLV(M,N,L)
       !                ELSE
       !                  FSP1(M,N)= VOR0(M,N,L)
       !                ENDIF
       !                PGQSP0(M,N,JDIV,L)= COMPLEX(0,0)
       !                PGQSP1(M,N,JDIV,L)= COMPLEX(0,0)
       !            ENDDO
       !          ENDDO
       !        ENDIF
       
                 CALL SPCOPY( PGQSP1(0:M1MAX,0:N1MAX,IPGQ,L), FSP1, N1)
         220 CONTINUE
       
       !   Surface pressure tendency
             DO 230 N=0,N1
               DO 230 M=0,MIN0(N,M1),1
                 FSP1(M,N)= PSSP0(M,N) &
                  + 2.0 * DT1 * DPSSP(M,N)
            
         230 CONTINUE
       
       !   Apply Robert filter
             IF (ROBFAC .NE. 0.) THEN
               CALL MODROBFIL( PSSP0, PSSP1, FSP1, ROBFAC)
               CALL SPCOPY( PSSP0, PSSP1, N1)
             ENDIF
       
             CALL SPCOPY( PSSP1, FSP1, N1)
       
             RETURN
             END

end module timestep_mod