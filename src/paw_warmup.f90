MODULE WARMUP_MODULE
!***********************************************************************
!***********************************************************************
!**  DATA TO BE SET BEFORE STARTING                                   **
!**    TON                                                            **
!**    NPULSE                                                         **
!***********************************************************************
!***********************************************************************
LOGICAL(4)            :: TON=.FALSE.
INTEGER(4)            :: NPULSE         ! #PULSES
INTEGER(4)            :: PULSE_LENGTH   ! #STEPS PER PULSE
REAL(8)               :: TARGET_TEMP    ! TARGET TEMPERATURE KELVIN
END MODULE WARMUP_MODULE
!
!     ..................................................................
      SUBROUTINE WARMUP$SETL4(ID,VAL)
!     ******************************************************************
!     ******************************************************************
      USE WARMUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'ON') THEN
        TON=VAL
      ELSE
        CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
        CALL ERROR$STOP('WARMUP$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WARMUP$SETI4(ID,VAL)
!     ******************************************************************
!     ******************************************************************
      USE WARMUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'NPULSE') THEN
        NPULSE=VAL
      ELSE IF(ID.EQ.'PULSE_LENGTH') THEN
        PULSE_LENGTH=VAL
      ELSE
        CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
        CALL ERROR$STOP('WARMUP$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WARMUP$SETR8(ID,VAL)
!     ******************************************************************
!     ******************************************************************
      USE WARMUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'TARGET_TEMP') THEN
        TARGET_TEMP=VAL
      ELSE
        CALL ERROR$MSG('KEYWORD NOT RECOGNIZED')
        CALL ERROR$STOP('WARMUP$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PAW_WARMUP_APPLY
!     ------------------------------------------------------------------------ 
!     THIS ROUTINE PREPS VARIOUS ARRAYS THAT ARE PASSED TO THE WARMUP ROUTINE
!     ------------------------------------------------------------------------ 
      USE WARMUP_MODULE
      IMPLICIT NONE
      LOGICAL(4),SAVE        :: FINISHED=.FALSE.
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      REAL(8)   ,ALLOCATABLE :: RM(:,:)
      REAL(8)   ,ALLOCATABLE :: RP(:,:)
      REAL(8)   ,ALLOCATABLE :: FORCE(:,:)
      REAL(8)                :: DELT
      REAL(8)   ,ALLOCATABLE :: RMASS(:)
!     ******************************************************************
      IF ((.NOT.TON).OR.FINISHED) RETURN
!
!     ==================================================================
!     == ALLOCATE ARRAYS                                              ==
!     ==================================================================
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      ALLOCATE(RM(3,NAT))
      ALLOCATE(RP(3,NAT))
      ALLOCATE(FORCE(3,NAT))
      ALLOCATE(RMASS(NAT))
!
!     ==================================================================
!     == COLLECT DATA                                                 ==
!     ==================================================================
      CALL ATOMS$GETR8('TIMESTEP',DELT)
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      CALL ATOMLIST$GETR8A('R(-)',0,3*NAT,RM)
      CALL ATOMLIST$GETR8A('R(+)',0,3*NAT,RP)
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCE)
      CALL ATOMLIST$GETR8A('MASS',0,NAT,RMASS)
!
!     ==================================================================
!     == RUN WARMUP                                                   ==
!     ==================================================================
      CALL PAW_WARMUP(NPULSE,PULSE_LENGTH,TARGET_TEMP,NAT,RM,R0,RP,FORCE &
     &               ,DELT,RMASS,FINISHED)
!
!     ==================================================================
!     == RETURN DATA                                                  ==
!     ==================================================================
      CALL ATOMLIST$SETR8A('R(+)',0,3*NAT,RP)
      CALL ATOMLIST$SETR8A('FORCE',0,3*NAT,FORCE)
!
!     ==================================================================
!     == CLEAN UP AND CLOSE DOWN                                      ==
!     ==================================================================
      DEALLOCATE(R0)
      DEALLOCATE(RM)
      DEALLOCATE(RP)
      DEALLOCATE(FORCE)
      DEALLOCATE(RMASS)
      RETURN
      END SUBROUTINE PAW_WARMUP_APPLY
!
!     .....................................................................
      SUBROUTINE PAW_WARMUP(NPULSE,NSTEP_IN,TFINAL,NAT,RM,R0,RP,FORCE  &
     &                     ,DELT,RMASS,TFINISHED)
!     *********************************************************************
!     ** SLOWLY WARM THE SYSTEM WHICH SINUSOIDAL PULSES.                 ** 
!     ** FIRST PULSE IS A RANDOM PULSE.  SUBSEQUENT PULSES ARE DESIGNED TO*
!     ** BE ORTHOGONAL TO CURRENT MOTION.                                **
!     ** NOTES:  WORKS BEST FROM AN OPTIMIZED STRUCTURE.                 **
!     **                                                                 **
!     ** P. M. MARGL                                                     **
!     *********************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(INOUT)   :: TFINISHED
      INTEGER(4),INTENT(IN)      :: NAT
      INTEGER(4),INTENT(IN)      ::  NSTEP_IN
      INTEGER(4),INTENT(INOUT)   :: NPULSE
      REAL(8)   ,INTENT(IN)      :: TFINAL, DELT
      REAL(8)   ,INTENT(IN)      :: RMASS(NAT)
      REAL(8)   ,INTENT(IN)      :: RM(3,NAT)
      REAL(8)   ,INTENT(IN)      :: R0(3,NAT)
      REAL(8)   ,INTENT(INOUT)   :: FORCE(3,NAT)
      REAL(8)   ,INTENT(INOUT)   :: RP(3,NAT)
      REAL(8)   ,PARAMETER       :: PI=4.D0*ATAN(1.D0)
      REAL(8)                    :: MASS(3,NAT)
      INTEGER(4),SAVE            :: NSTEP
      REAL(8)   ,ALLOCATABLE,SAVE:: DELTAV(:,:)
      REAL(8)   ,ALLOCATABLE,SAVE:: DELTAVO(:,:)
      REAL(8)   ,ALLOCATABLE,SAVE:: VPRO(:,:)
      REAL(8)   ,ALLOCATABLE,SAVE:: DWORK(:,:)
      REAL(8)   ,ALLOCATABLE,SAVE:: V(:,:)
      REAL(8)          :: XMIN, A, B, TEMPER, TEMPER1,RKAPPA, SVAR, RKB, SXVAR
      INTEGER(4)       :: ISTEP, I, J, NFILO, IPULSE, I1, I2, J1, J2, JC, IC
      INTEGER(4), SAVE :: ITOTAL = 0
      REAL(8)          :: TWOPI
      INTEGER(4)       :: NG
!     -----------------------------------------------------------------------
      TWOPI=2.D0*PI
      NG=3*NAT
!     ----------------------------------------------------
!     PUT RMASS ON ARRAY WHICH IS CONFORMABLE WITH V AND R
!     ----------------------------------------------------
      DO I=1,NAT
       DO J=1,3
        MASS(J,I)=RMASS(I)
       ENDDO
      ENDDO
!     -----------------------------------
!     TOTAL NUMBER OF TIMESTEPS IN WARMUP
!     -----------------------------------
      ITOTAL=ITOTAL+1
!     --------------
!     INITIALIZATION
!     --------------
      IF (ITOTAL .EQ. 1) THEN
       NSTEP = NSTEP_IN
       ALLOCATE(DELTAV(3,NAT))
       ALLOCATE(DELTAVO(3,NAT))
       ALLOCATE(VPRO(3*NAT,3*NAT))
       ALLOCATE(V(3,NAT))
       ALLOCATE(DWORK(3,NAT))
       IF (NPULSE .EQ. 0) THEN
        NPULSE=3*NAT-6
       ENDIF
!
       CALL FILEHANDLER$UNIT('PROT',NFILO)
       WRITE(NFILO,*)
       WRITE(NFILO,FMT='("WARMUP SETTINGS")')
       WRITE(NFILO,FMT='("VERSION 1.0")')
       WRITE(NFILO,FMT='("=============================")')
       WRITE(NFILO,FMT='("NUMBER OF PULSES:     ",I4)')NPULSE
       WRITE(NFILO,FMT='("NUMBER OF STEPS/PULSE:",I4)')NSTEP
       WRITE(NFILO,FMT='("TARGET TEMPERATURE: ",F6.1)')TFINAL
       WRITE(NFILO,FMT='("=============================")')
       WRITE(NFILO,*)
!
      ENDIF
!
      CALL CONSTANTS('KB',RKB)
!     -------------
!     CURRENT PULSE
!     -------------
      IPULSE=(ITOTAL-1)/NSTEP+1
!     ---------------------------------
!     CHECK FOR END OF PULSING SEQUENCE
!     ---------------------------------
      IF (IPULSE .GT. NPULSE) THEN
       CALL FILEHANDLER$UNIT('PROT',NFILO)
       WRITE(NFILO,*)'ION WARMUP FINISHED'
       TFINISHED=.TRUE.
       DEALLOCATE(DELTAV)
       DEALLOCATE(DELTAVO)
       DEALLOCATE(VPRO)
       DEALLOCATE(V)
       DEALLOCATE(DWORK)
       RETURN
      ENDIF
!     ------------------------------------------ 
!     NUMBER OF TIMESTEPS SPENT IN CURRENT PULSE
!     ------------------------------------------ 
      ISTEP=ITOTAL-(IPULSE-1)*NSTEP
!     ---------------------------
!     CALCULATE SHOULD-BE DELTA T 
!     MAKE A SINEWAVE
!     ---------------------------
      RKAPPA=TFINAL/(DBLE(NPULSE)*DBLE(NSTEP))*(DSIN(TWOPI*DBLE(ITOTAL)   &
&             /DBLE(NSTEP) - PI/2.D0)+1.D0)
!     -------------------------------
!     GET A NEW VECTOR FOR EACH PULSE 
!     -------------------------------
      IF (ISTEP .EQ. 1) THEN
       CALL FILEHANDLER$UNIT('PROT',NFILO)
       WRITE(NFILO,'("WARMING UP IONS - APPLYING PULSE NO.: ",I4)')IPULSE 
       CALL PAW_WARMUP_VECTOR(IPULSE,NAT,RM,RP,RMASS,DELTAVO)
!      -------------------------
! ---- COMPUTE NORMALIZED DELTAV
!      -------------------------
       SVAR=SUM(DELTAVO*DELTAVO)
       SVAR=1.D0/SQRT(SVAR)
       DELTAVO=DELTAVO*SVAR
       DELTAV=DELTAVO
      ENDIF
!     -------------------------------------------------------------
!     COMPUTE PROJECTOR UPON M|V><V|M, ABUSE RP FOR TEMPORARY SPACE
!     -------------------------------------------------------------
      VPRO=0.D0
      SVAR=0.D0
      RP=MASS*DELTAV
      SVAR=SUM(RP*RP)
!     --------------------------------------------------
!     NORMALIZE FOR CONSTRUCTION OF PROJECTOR
!     WHICH IS USED TO REMOVE ALL THE FORCES WHICH POINT
!     IN THE DIRECTION OF |M|DELTAV>
!     --------------------------------------------------
      SVAR=1.D0/SQRT(SVAR)
      RP=RP*SVAR
      IC=0
      DO I1=1,NAT
       DO J1=1,3
        IC=IC+1
        JC=0
        DO I2=1,NAT
         DO J2=1,3
          JC=JC+1
          VPRO(JC,IC)=RP(J1,I1)*RP(J2,I2)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!     ------------------
! ----PROJECT UPON FORCE
!     ------------------
      IC=0
      DO I1=1,NAT
       DO J1=1,3
        DWORK(J1,I1)=0.D0
        IC=IC+1
        JC=0
        DO I2=1,NAT
         DO J2=1,3
          JC=JC+1
          DWORK(J1,I1)=DWORK(J1,I1)+VPRO(JC,IC)*FORCE(J2,I2)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!     -------------------
!     SUBTRACT PROJECTION
!     -------------------
      FORCE=FORCE-DWORK
      SXVAR=SUM(RP*FORCE)
!     -------------------------------------------
!     CHECK ORTHOGONALITY BETWEEN FORCE AND M|DV>
!     -------------------------------------------
      IF (ABS(SXVAR).GT.1D-7) THEN
       CALL ERROR$MSG('FAILED TO ORTHOGONALIZE FORCES AND M|DV>')
       CALL ERROR$STOP('ATOMIC:MM_WARMUPAPPLY')
      ENDIF
!     -------------------------------------------------------------
!     PROPAGATE IONS DISREGARDING FORCES IN DIRECTION OF EXCITATION
!     -------------------------------------------------------------
      RP=2.D0*R0-RM+FORCE*DELT**2.D0/MASS
!     -------------------------------------------
!     COPY ORIGINAL DELTAV VECTOR INTO CURRENT DV
!     -------------------------------------------
      DELTAV=DELTAVO
!     ----------------------------------------------
!     COMPUTE V FROM PROJECTION-PROPAGATED POSITIONS
!     ----------------------------------------------
      V=(RP-R0)/DELT
      TEMPER=SUM(MASS*V**2.D0)/(DBLE(NG)*RKB)
!     ----------------------------------------------
!     CALCULATE XMIN BY SOLVING QUADRATIC SYSTEM
!     ----------------------------------------------
      A=2.D0/DBLE(NG)/RKB*SUM(MASS*DELTAVO*V)
      B=1.D0/DBLE(NG)/RKB*SUM(DELTAVO*MASS*DELTAVO)
      XMIN=-A/B/2.D0+SQRT((A/B)**2/4.D0+RKAPPA/B)
      PRINT*,'TEMPERATURE INCREMENT FACTOR ',XMIN
!
!     THIS IS THE OTHER SOLUTION, NOT USED
!     XMIN=-A/B/2.D0-SQRT((A/B)**2/4.D0+RKAPPA/B)
!     PRINT*,'XMIN FROM QUAD2 ',XMIN
!
!     ----------------------------
!     COMPUTE DELTAV NOW
!     ----------------------------
      DELTAV=DELTAV*XMIN
!
!     -------------------------------------
!     COMPUTE RP FROM R0 AND NEW VELOCITIES
!     -------------------------------------
      RP=R0+(V+DELTAV)*DELT
!
!     ----------------------------------------------------
!     STORE PREVIOUS TEMPERATURE AWAY FOR LATER COMPARISON 
!     ----------------------------------------------------
      TEMPER1=TEMPER
!     -----------------------------------------------------------
!     -------------------------------------
!     COMPUTE V FROM CURRENT TIME STEP ONLY
!     -------------------------------------
      V=(RP-R0)/DELT
      TEMPER=SUM(MASS*V**2.D0)/(DBLE(NG)*RKB)
!
      IF (ABS(TEMPER-TEMPER1-RKAPPA).GT.1.D-5) THEN
       CALL ERROR$MSG('WARMUP: DESIRED DELTA T NOT EFFECTED')
      CALL ERROR$STOP('IN WARMUP_APPLY')
      ENDIF
!
      RETURN
      END SUBROUTINE PAW_WARMUP 

!     .................................................................
      SUBROUTINE PAW_WARMUP_VECTOR(IPULSE,NAT,RM,RP,RMASS,DELTAV)
!     *****************************************************************
!     ** RETURNS A VELOCITY VECTOR DELTAV                            **
!     ** DELTAV IS CHOSEN RANDOMLY IN THE FIRST PULSE, ALL SUBSEQUENT**
!     ** DELTAV-S ARE GENERATED IN DIRECTIONS NORMAL TO THE CURRENT  **
!     ** VELOCITIES                                                  **
!     *****************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: IPULSE, NAT
      REAL(8)   ,INTENT(IN)   :: RM(3,NAT)
      REAL(8)   ,INTENT(IN)   :: RP(3,NAT)
      REAL(8)   ,INTENT(INOUT):: DELTAV(3,NAT)
      REAL(8)   ,INTENT(IN)   :: RMASS(NAT)
      REAL(8)   ,ALLOCATABLE  :: VPRO(:,:)
      REAL(8)   ,ALLOCATABLE  :: DWORK(:,:)
      REAL(8)   ,ALLOCATABLE  :: RCHECK(:,:)
      INTEGER(4)              :: I, J, I1, J1, I2, J2
      REAL(8)                 :: FACTOR, SVAR, RANF
!     *****************************************************************
!
!     ---------------------
!     GENERATE FIRST VECTOR
!     ---------------------
      IF (IPULSE .EQ. 1) THEN
       DO I=1,NAT
        DO J=1,3
         CALL RANDOM_NUMBER(RANF)
         DELTAV(J,I)=RANF-0.5D0
        ENDDO
       ENDDO
!      ---------------------------
!      GENERATE SUBSEQUENT VECTORS
!      ---------------------------
      ELSE
!      
       ALLOCATE(VPRO(3*NAT,3*NAT))
       ALLOCATE(DWORK(3,NAT))
       ALLOCATE(RCHECK(3,NAT))
!      --------------------------------------------------------
!      CALCULATE VELOCITIES WEIGHTED WITH SQRT(M)
!      BE AWARE THAT DELTAV CONTAINS NOW TEMPORARILY A VELOCITY
!      --------------------------------------------------------
       SVAR=0.D0
       DO I=1,NAT
        FACTOR=SQRT(RMASS(I))
        DO J=1,3
         DELTAV(J,I)=(RP(J,I)-RM(J,I))*FACTOR
         SVAR=SVAR+DELTAV(J,I)**2
        ENDDO
       ENDDO
       SVAR=SQRT(SVAR)
       DO I=1,NAT
        DO J=1,3
         DELTAV(J,I)=DELTAV(J,I)/SVAR
         RCHECK(J,I)=DELTAV(J,I)
        ENDDO
       ENDDO
!      -------------------------------------------
!      CONSTRUCT PROJECTOR UPON WEIGHTED VELOCITY
!      -------------------------------------------
       I=0
       J=0
       DO I1=1,NAT
        DO J1=1,3
         I=I+1
         J=0
         DO I2=1,NAT
          DO J2=1,3
           J=J+1
           VPRO(J,I)=DELTAV(J2,I2)*DELTAV(J1,I1)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
!      ------------------------------------------
!      INVENT SOME RANDOM VECTOR AND NORMALIZE IT
!      ------------------------------------------
       SVAR=0.D0
       DO I=1,NAT
        DO J=1,3
         CALL RANDOM_NUMBER(RANF)
         DELTAV(J,I)=RANF-0.5D0
         DWORK(J,I)=0.D0
         SVAR=SVAR+DELTAV(J,I)**2
        ENDDO
       ENDDO
       SVAR=SQRT(SVAR)
       DO I=1,NAT
        DO J=1,3
         DELTAV(J,I)=DELTAV(J,I)/SVAR
        ENDDO
       ENDDO
!      -----------------------------
!      PROJECT UPON PREVIOUS VECTORS
!      -----------------------------
       I=0
       DO I1=1,NAT
        DO J1=1,3
         J=0
         I=I+1
         DO I2=1,NAT
          DO J2=1,3
           J=J+1
           DWORK(J1,I1)=DWORK(J1,I1)+VPRO(J,I)*DELTAV(J2,I2)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       SVAR=0.D0
       DO I=1,NAT
        DO J=1,3
         DELTAV(J,I)=DELTAV(J,I)-DWORK(J,I)
         SVAR=SVAR+RCHECK(J,I)*DELTAV(J,I)
        ENDDO
       ENDDO
       IF (ABS(SVAR).GT.1.D-10) THEN
        CALL ERROR$MSG('FAILED TO FIND ORTHOGONAL VELOCITIES!')
        CALL ERROR$STOP('IN WARMUP')
       ENDIF
       DEALLOCATE(VPRO)
       DEALLOCATE(DWORK)
       DEALLOCATE(RCHECK)
      ENDIF
!     -------------------------------------------
!     DIVIDE BY SQRT(MASS) BEFORE HANDING IT BACK
!     -------------------------------------------
      DO I=1,NAT
       FACTOR=1.D0/SQRT(RMASS(I))
       DO J=1,3
        DELTAV(J,I)=DELTAV(J,I)*FACTOR
       ENDDO
      ENDDO

      RETURN

      END SUBROUTINE PAW_WARMUP_VECTOR
