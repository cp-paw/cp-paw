!*******************************************************************************
!*******************************************************************************
!*****                                                                      ****
!*****          THERMOSTAT OBJECT APPLIES THE NOSE THERMOSTAT               ****
!*****                                                                      ****
!*******************************************************************************
!*******************************************************************************
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE THERMOSTAT_MODULE
!*******************************************************************************
!**                                                                           **
!**  THERMOSTAT IMPLEMENTS THE NOSE THERMOSTAT                                **
!**                                                                           **
!**  FUNCTIONS :                                                              **
!**    NEW(ID)        ADD NEW THERMOSTAT                                      **
!**    SELECT(ID)     SELECT EXISTING THERMOSTAT                              **
!**    DELETE         DELETE SELECTED THERMOSTAT                              **
!**    COMPOSE                                                                **
!**    REPORT(NFIL)                                                           **
!**    PROPAGATE                                                              **
!**    SWITCH                                                                 **
!**    READ(NFIL,NFILO,TCHK)                                                  **
!**    WRITE(NFIL,NFILO,TCHK)                                                 **
!**                                                                           **
!**                                                                           **
!**    THERMOSTAT$NEW(ID)                                                     **
!**    THERMOSTAT$SETR8('TIMESTEP',DT)                                        **
!**    THERMOSTAT$SETR8('MASS',Q)                                             **
!**    THERMOSTAT$SETR8('FRICTION',ANNE)                                      **
!**    THERMOSTAT$SETR8('GFREE',GFREE)    (OPTIONAL,ONLY FOR REPORT)          **
!**                                                                           **
!**    THERMOSTAT$SELECT(ID)                                                  **
!**    THERMOSTAT$SETR8('EKIN(SYSTEM)')                                       **
!**    THERMOSTAT$PROPAGATE                                                   **
!**    THERMOSTAT$SWITCH                                                      **
!**                                                                           **
!**  REMARKS:                                                                 **
!**    1) PERIOD BECOMES EFFECTIVE ONLY THROUGH 'INITIALIZE'                  **
!**    2) GFREE AND TEMP ARE ACTUALIZED IN EACH ITERATION                     **
!**                                                                           **
!**                                                                           **
!**    3) IF THE 'COOLING' IS OBTAINED BEFORE CALLING THE FUNCTION            **
!**    THERMOSTAT$PROPAGATE THE VALUE IS OBTAINED FROM A POLYNOMIAL           **
!**    EXTRAPOLATION OF THE THERMOSTAT VARIABLE.                              **
!**    IF THE 'COOLING' IS OBTAIND AFTER CALLING THERMOSTAT$PROPAGATE         **
!**    IT IS OBTAINED FROM THE ACTUAL VALUES OF THE THERMOSTA VARIABLE        **
!**    4) THE THERMOSTAT CAN OPERATE IN TWO MODES. ONE IS A PURELY            **
!**    COOLING THERMOSTAT FOR THE WAVE FUNCTIONS. THIS OPTION IS SET          **
!**    BY SETL4('COOLONLY',.TRUE.)                                            **
!**    5) DURING READING  AND WRITING, THE THERMOSTAT ASSUMES THAT            **
!**    IT LIVES ON THE PROCESSOR GROUP MONOMER (SEE MPE OBJECT)               **
!**                                                                           **
!*******************************************************************************
TYPE THERMOSTAT_TYPE
! == INDIRECT SETTING: APPLIED ONLY THROUGH THERMOSTAT%INITIALIZE ======
  REAL(8)      :: GFREE       ! CONTROLLED SYSTEM IS 0.5*GFREE*TEMP
! == ACTUAL SETTING USED ===============================================
  CHARACTER(32):: ID          ! THERMOSTAT ID
  LOGICAL(4)   :: ON          ! ON/OFF SWITCH FOR THE THERMOSTAT
  LOGICAL(4)   :: TWAVE       ! USE SPECIAL WAVE FUNCTION THERMOSTAT
  LOGICAL(4)   :: STOP        ! INITIAL DX/DT IS SET TO ZERO 
  REAL(8)      :: DT          ! TIME STEP FOR X(T)
  REAL(8)      :: Q           ! MASS OF THE NOSE VARIABLE
  REAL(8)      :: FRICTION    ! FRICTION ACTING ON X(T)
! == ARGUMENTS =========================================================
  REAL(8)      :: EKIN_TARGET ! TARGET KINETIC ENERGY 
  REAL(8)      :: EKIN_ACTUAL ! ACTUAL KINETIC ENERGY
  REAL(8)      :: COOLING     ! FRICTION FOR THE CONTROLLED SYSTEM
  REAL(8)      :: EPOT        ! POTENTIAL THERMOSTAT ENERGY
  REAL(8)      :: EKIN        ! KINETIC THERMOSTAT ENERGY
  REAL(8)      :: EDISS       ! ENERGY DISSIPATED BY THE THERMOSTAT 
! == STATE =============================================================
  REAL(8)      :: XMMM        ! X(T-3*DELTA)
  REAL(8)      :: XMM         ! X(T-2*DELTA)
  REAL(8)      :: XM          ! X(T-DELTA)
  REAL(8)      :: X0          ! X(T)
  REAL(8)      :: XP          ! X(T+DELTA)
! == NEXT ==============================================================
  TYPE(THERMOSTAT_TYPE),POINTER :: NEXT
END TYPE THERMOSTAT_TYPE
LOGICAL(4)                    :: TINI=.FALSE.
TYPE(THERMOSTAT_TYPE),POINTER :: FIRST_THIS
TYPE(THERMOSTAT_TYPE),POINTER :: THIS
END MODULE THERMOSTAT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$NEW(ID)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
!     **************************************************************************
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        ALLOCATE(FIRST_THIS)
        THIS=>FIRST_THIS
      ELSE
        THIS=>FIRST_THIS
        DO 
          IF(THIS%ID.EQ.ID) THEN
            CALL ERROR$MSG('CANNOT CREATE THERMOSTAT WITH THE SAME NAME')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('THERMOSTAT$NEW')
          END IF
          IF(ASSOCIATED(THIS%NEXT)) THEN
            THIS=>THIS%NEXT
          ELSE
            ALLOCATE(THIS%NEXT)
            THIS=>THIS%NEXT
            EXIT
          END IF
        ENDDO
      END IF
      THIS%GFREE=0.D0
!   
      THIS%ID      =ID
      THIS%ON      =.FALSE.      
      THIS%TWAVE   =.FALSE.
      THIS%DT      =0.D0
      THIS%Q       =0.D0
      THIS%FRICTION=0.D0
      THIS%STOP    =.FALSE.      
!
      THIS%EKIN_TARGET=0.D0
      THIS%EKIN_ACTUAL=0.D0
      THIS%COOLING=0.D0
      THIS%EPOT=0.D0
      THIS%EKIN=0.D0
      THIS%EDISS=0.D0
!
      THIS%XMMM=0.D0
      THIS%XMM=0.D0
      THIS%XM=0.D0
      THIS%X0=0.D0
      THIS%XP=0.D0
!
      NULLIFY(THIS%NEXT)
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$DELETE
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      TYPE(THERMOSTAT_TYPE),POINTER :: PREVIOUS
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$DELETE')
      END IF
      IF(ASSOCIATED(THIS,FIRST_THIS)) THEN
        FIRST_THIS=THIS
        DEALLOCATE(THIS)
      ELSE
        PREVIOUS=>FIRST_THIS
        DO WHILE (.NOT.ASSOCIATED(PREVIOUS,THIS))
          PREVIOUS=>PREVIOUS%NEXT
        END DO
        PREVIOUS%NEXT=>THIS%NEXT
        DEALLOCATE(THIS)
      END IF 
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$SELECT(ID)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
!     **************************************************************************
      IF(.NOT.TINI) THEN
        CALL ERROR$MSG('CREATE THERMOSTATS BEFORE SELECTING')
        CALL ERROR$STOP('THERMOSTAT$SELECT')
      END IF
!
!     == UNSELECT ==============================================================
      IF(ID.EQ.'NONE') THEN
        NULLIFY(THIS)
      END IF
!
!     == FIND THERMOSTAT AND SELECT  ===========================================
      THIS=>FIRST_THIS
      DO WHILE (THIS%ID.NE.ID)
        IF(.NOT.ASSOCIATED(THIS%NEXT)) THEN
          CALL ERROR$MSG('THERMOSTAT DOES NOT EXIST')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('THERMOSTAT$SELECT')
        END IF
        THIS=>THIS%NEXT
      ENDDO
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$REPORT(NFIL)
!     **************************************************************************
!     **  WRITE REPORT                                                        **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: KELVIN
      REAL(8)               :: SVAR
      REAL(8)               :: OMEGA0
!     ******************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$REPORT')
      END IF
      IF(.NOT.THIS%ON) RETURN
!
      CALL CONSTANTS('KB',KELVIN)
      CALL REPORT$TITLE(NFIL,'THERMOSTAT: '//TRIM(THIS%ID))
      IF(THIS%TWAVE) THEN
        CALL REPORT$STRING(NFIL,'SPECIAL WAVE FUNCTION THERMOSTAT:COOLING ONLY')
      END IF
      IF(THIS%STOP) THEN
        CALL REPORT$CHVAL(NFIL,'INITIAL VELOCITIES SET TO','ZERO')
      END IF
      IF(THIS%GFREE.NE.0.D0) THEN
        SVAR=THIS%EKIN_TARGET/(0.5D0*THIS%GFREE*KELVIN)
        CALL REPORT$R8VAL(NFIL,'TEMPERATURE',SVAR,'K')
        CALL REPORT$R8VAL(NFIL,'#(DEGREES OF FREEDOM)',THIS%GFREE,' ')
      ELSE
        SVAR=THIS%EKIN_TARGET
        CALL REPORT$R8VAL(NFIL,'TARGET KINETIC ENERGY',SVAR,'H')
      END IF
      OMEGA0=SQRT(4.D0*THIS%EKIN_TARGET/THIS%Q)
      SVAR=2.D0*PI/OMEGA0
      CALL REPORT$R8VAL(NFIL,'OSCILATION PERIOD',SVAR,'A.U.')
      IF(THIS%FRICTION.NE.0.D0) THEN
        CALL REPORT$R8VAL(NFIL,'FRICTION ON NOSE VARIABLE',THIS%FRICTION,' ')
        SVAR=THIS%FRICTION/(OMEGA0*THIS%DT)
        CALL REPORT$R8VAL(NFIL,'SAME DIVIDED BY CRITICAL FRICTION',SVAR,' ')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$CONVERT(DT,TARGET,PERIOD,Q,FRICTION)
!     **************************************************************************
!     **  CREATED THE INPUT FOR THERMOSTAT FROM PHYSICAL INPUT                **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: DT       ! TIME STEP
      REAL(8)   ,INTENT(IN) :: TARGET   ! TARGET KINETIC ENERGY
      REAL(8)   ,INTENT(IN) :: PERIOD   ! PERIOD OF A NOSE OSCILLATION
      REAL(8)   ,INTENT(OUT):: Q        ! THERMOSTAT MASS
      REAL(8)   ,INTENT(OUT):: FRICTION ! THERMOSTAT FRICTION FOR
                                        ! CRITICAL DAMPING
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: OMEGA0
      REAL(8)               :: ALPHA
!     **************************************************************************
      OMEGA0=2.D0*PI/PERIOD
      Q=4.D0*TARGET/OMEGA0**2   ! OMEGA0=SQRT(4*<T>/Q)
      ALPHA=1.D0/OMEGA0         ! CRITICAL DAMPING
      FRICTION=ALPHA*DT/2.D0    ! FRICTION FOR CRITICAL DAMPING
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$SCALEGFREE(GFREE)
!     **************************************************************************
!     **  INTRODUCES THE NUMBER OF DEGREES OF FREEDOM                         **
!     **  SO THAT THE KINETIC ENERGY IS SCALED UP,                            **
!     **  WHILE TEMPERATURE, OSCILLATION PERIOD, AND DECAY TIME               **
!     **  REMAIN UNCHANGED                                                    **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: GFREE    ! #(DEGREES OF FREEDOM)
      REAL(8)               :: SCALE
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$SCALEGFREE')
      END IF
      IF(THIS%Q.LE.0.D0) THEN
        CALL ERROR$MSG('THERMOSTAT MASS IS ZERO')
        CALL ERROR$CHVAL('ID',THIS%ID)
        CALL ERROR$STOP('THERMOSTAT$SCALEGFREE')
      END IF
      IF(THIS%EKIN_TARGET.LE.0.D0) THEN
        CALL ERROR$MSG('TARGET KINETIC ENERGY IS ZERO')
        CALL ERROR$CHVAL('ID',THIS%ID)
        CALL ERROR$STOP('THERMOSTAT$SCALEGFREE')
      END IF
!     ==========================================================================
!     ==  NOW DO THE RESCALING                                                ==
!     ==========================================================================
      IF(THIS%GFREE.EQ.0.D0) THIS%GFREE=1.D0
      SCALE=GFREE/THIS%GFREE        
      THIS%GFREE      =THIS%GFREE      *SCALE
      THIS%EKIN_TARGET=THIS%EKIN_TARGET*SCALE
IF(THIS%ID.EQ.'ATOMS')PRINT*,'THERMOSTAT$SCALEGFREE',GFREE,THIS%Q,SCALE
      THIS%Q          =THIS%Q          *SCALE
IF(THIS%ID.EQ.'ATOMS')PRINT*,'THERMOSTAT$SCALEGFREE',THIS%Q
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$PROPAGATE()
!     **************************************************************************
!     **  CALCULATE XNOSP AND THE VELOCITY OF THE NOSE VARIABLE               **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: DT
      REAL(8)               :: X0,XM,XP
      REAL(8)               :: ANNEX  ! FRICTION FOR THE NOSE THERMOSTAT
      REAL(8)               :: QMASS  ! THERMOSTAT MASS
      REAL(8)               :: TARGET ! TARGET KINETIC ENERGY OF THE SYSTEM
      REAL(8)               :: EKIN
      REAL(8)               :: SVAR1,SVAR2,SVAR3
      REAL(8)               :: VNOS
      REAL(8)               :: ANNER  
      REAL(8)               :: PERIOD
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$PROPAGATE')
      END IF
!
      IF(.NOT.THIS%ON) RETURN
!
      IF(THIS%Q.LE.0.D0) THEN
        CALL ERROR$MSG('THERMOSTAT MASS IS ZERO')
        CALL ERROR$CHVAL('ID',THIS%ID)
        CALL ERROR$STOP('THERMOSTAT$PROPAGATE')
      END IF
      IF(THIS%EKIN_TARGET.LE.0.D0) THEN
        CALL ERROR$MSG('TARGET KINETIC ENERGY IS ZERO')
        CALL ERROR$CHVAL('ID',THIS%ID)
        CALL ERROR$STOP('THERMOSTAT$PROPAGATE')
      END IF
      IF(THIS%DT.LE.0.D0) THEN
        CALL ERROR$MSG('TIME STEP MUST BE POSITIVE AND FINITE')
        CALL ERROR$CHVAL('ID',THIS%ID)
        CALL ERROR$STOP('THERMOSTAT$PROPAGATE')
      END IF
!
!     ==========================================================================
!     == STOP NOSE COORDINATE (SET FRICTION TO ZERO)                          ==
!     ==========================================================================
      IF(THIS%STOP) THEN
        THIS%XM  =THIS%X0
        THIS%XMM =THIS%X0
        THIS%XMMM=THIS%X0
        THIS%STOP=.FALSE.
      END IF
!
!     ==========================================================================
!     ==  ADJUST DIRECT THERMOSTAT PARAMETERS                                 ==
!     ==========================================================================
      PERIOD=PI*SQRT(THIS%Q/THIS%EKIN_TARGET)
      IF(PERIOD.LT.3.D0*THIS%DT) THEN
        PRINT*,'WARNING! THERMOSTAT HAS TOO LARGE TIME STEP'
      END IF
!
!     ==========================================================================
!     ==  CALCULATE XNOSP AND THE VELOCITY OF THE NOSE VARIABLE               ==
!     ==========================================================================
      DT=THIS%DT
      X0=THIS%X0
      XM=THIS%XM
      EKIN=THIS%EKIN_ACTUAL
      ANNEX=THIS%FRICTION
      QMASS=THIS%Q  
      TARGET=THIS%EKIN_TARGET
!
!     ==  PROPAGATE THERMOSTAT VARIABLE ========================================
      SVAR1=2.D0/(1.D0+ANNEX)
      SVAR2=1.D0-SVAR1
      SVAR3=DT**2/QMASS/(1.D0+ANNEX)
      XP=SVAR1*X0+SVAR2*XM+SVAR3*2.D0*(EKIN-TARGET)
!
!     == UPDATE FRICTION =======================================================
      VNOS =(XP-XM)/(2.D0*DT)
      IF(THIS%TWAVE) THEN  
        IF(VNOS.LT.0) THEN !SWITCH OFF WHEN HEATING
          XP=X0
          THIS%XM  =X0
          THIS%XMM =X0
          THIS%XMMM=X0
          VNOS=0.D0
        END IF
      END IF
      ANNER=0.5D0*DT*VNOS
      ANNER=MAX(ANNER,-0.99999D0)
      ANNER=MIN(ANNER,1.D0)
!
!     ==========================================================================
!     ==  COLLECT RESULTS                                                     ==
!     ==========================================================================
      THIS%XP     = XP
      THIS%COOLING= ANNER
      THIS%EKIN   = 0.5D0*THIS%Q*VNOS**2 
      THIS%EPOT   = 2.D0*THIS%EKIN_TARGET*X0
      IF(THIS%TWAVE) THIS%EPOT=0.D0
      THIS%EDISS  = THIS%EDISS+ANNEX*0.5D0*QMASS*VNOS**2
THIS%EDISS=0.D0
!
!     ==========================================================================
!     ==  PRINTOUT FOR TEST                                                   ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,'("THERMOSTAT: XNOS :",F5.2," VNOS ",F5.2," ENOSEP ",F10.5)') &
     &          X0,VNOS,THIS%EKIN+THIS%EPOT
      END IF
PRINT*,'THERM COOLING',ANNER        
PRINT*,'THERM XP     ',XP,X0,XM
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$SWITCH
!     **************************************************************************
!     **  SWITCH NOSE COORDINATES                                             **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      REAL(8)            :: X0
      REAL(8)            :: XM
      REAL(8)            :: XMM,XMMM
      REAL(8)            :: XP
      REAL(8)            :: XDOT,XDOTLAST
      REAL(8)            :: DT
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$SWITCH')
      END IF
      IF(.NOT.THIS%ON) RETURN
!
!     ==========================================================================
!     ==  SWITCH VARIABLES                                                    ==
!     ==========================================================================
      THIS%XMMM=THIS%XMM
      THIS%XMM =THIS%XM
      THIS%XM  =THIS%X0
      THIS%X0  =THIS%XP
!
!     ==========================================================================
!     ==  PREDICT FRICTION FOR THE NEXT TIME STEP                             ==
!     ==  (REQUIRED IF THE PROPAGATION CANNOT BE REPEATED)                    ==
!     ==========================================================================
      DT =THIS%DT
      X0 =THIS%X0
      XM =THIS%XM
      XMM=THIS%XMM 
      XMMM=THIS%XMMM
!     XP=2.D0*X0-XM
!     XP=3.D0*X0-3.D0*XM+XMM
      XP=4.D0*X0-6.D0*XM+4.D0*XMM-XMMM
      THIS%XP     =XP
      XDOT=(XP-XM)/(2.D0*DT)
      IF(THIS%TWAVE.AND.XDOT.LT.0) THEN
        XDOTLAST=(X0-XMM)/(2.D0*DT)
!       ========================================================================
!       == USE THE AVERAGE FRICTION AMONG THE TIME STEPS                      ==
!       == ASSUMING A LINEAR INTERPOLATION OF XDOT,                           ==
!       == WHICH IS TRUNCATED AT ITS ZERO.                                    ==
!       == THE HISTOGRAMM SHAPED FRICTION OF HALF STEP IS SUBTRACTED.         ==
!       ========================================================================
        IF(XDOT.NE.-XDOTLAST) THEN
!         XDOT=0.5D0*XDOTLAST*XDOT/(XDOTLAST-XDOT)
!RINT*,'XDOT....',XDOT
!       ELSE
          XDOT=0.D0
        END IF
      END IF
      THIS%COOLING=XDOT*DT/2.D0
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$WRITE(NFIL,NFILO,TCHK)
!     **************************************************************************
!     **  WRITES DYNAMICAL VARIABLES ON THE RESTART FILE                      **
!     **************************************************************************
      USE RESTART_INTERFACE
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: NFIL  ! RESTART FILE UNIT
      INTEGER(4)            ,INTENT(IN) :: NFILO ! PROTOCOL FILE UNIT 
      LOGICAL(4)            ,INTENT(OUT):: TCHK  ! WRITTEN/NOT WRITTEN
      TYPE (SEPARATOR_TYPE),PARAMETER   :: MYSEPARATOR &
                =SEPARATOR_TYPE(1,'THERMOSTAT','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE)             :: SEPARATOR
      REAL(8)                           :: X0,XM,XMM
      INTEGER(4)                        :: NTASKS,THISTASK
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$WRITE')
      END IF
      TCHK=.FALSE. 
      IF(.NOT.THIS%ON) RETURN
!
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     == COLLECT DATA ==========================================================
      X0 =THIS%X0
      XM =THIS%XM
      XMM=THIS%XMM
!
!     == WRITE DATA TO FILE  ===================================================
      IF(THISTASK.EQ.1) THEN
        SEPARATOR=MYSEPARATOR
        SEPARATOR%NAME=THIS%ID
        CALL RESTART$WRITESEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
        WRITE(NFIL)X0,XM,XMM
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$READ(NFIL,NFILO,TCHK)
!     **************************************************************************
!     **  READS DYNAMICAL VARIABLES FROM THE RESTART FILE                     **
!     **************************************************************************
      USE RESTART_INTERFACE
      USE THERMOSTAT_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: NFIL 
      INTEGER(4)            ,INTENT(IN) :: NFILO
      LOGICAL(4)            ,INTENT(OUT):: TCHK
      TYPE (SEPARATOR_TYPE),PARAMETER   :: MYSEPARATOR &
                =SEPARATOR_TYPE(1,'THERMOSTAT','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE)             :: SEPARATOR
      REAL(8)                           :: X0,XM,XMM
      INTEGER(4)                        :: NTASKS,THISTASK
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$READ')
      END IF
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     == READ DATA FROM FILE ===================================================
      TCHK=THIS%ON
      SEPARATOR=MYSEPARATOR
      SEPARATOR%NAME=THIS%ID
      IF(THISTASK.EQ.1)CALL RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
!
!     == RETURN AND LEAVE THERMOSTAT VARIABLES UNCHANGED =======================
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(.NOT.TCHK)RETURN
!
!     == READ DATA FROM FILE ===================================================
      IF(SEPARATOR%VERSION.NE.MYSEPARATOR%VERSION) THEN
        CALL ERROR$STOP('THERMOSTAT$READ')
      END IF
      IF(THISTASK.EQ.1) THEN
        READ(NFIL)X0,XM,XMM
      END IF
!
!     == DISTRIBUTE ARRAY TO ALL PARALLEL TASKS ================================
      CALL MPE$BROADCAST('MONOMER',1,X0)
      CALL MPE$BROADCAST('MONOMER',1,XM)
      CALL MPE$BROADCAST('MONOMER',1,XMM)
!
!     == SET  DATA E ARRAY TO ALL PARALLEL TASKS ===============================
      THIS%X0=X0
      THIS%XM=XM
      THIS%XMM=XMM
      THIS%XMMM=XMM
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$SETR8(ID,VAL)
!     **************************************************************************
!     **  SET DATA                                                            **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$SETR8')
      END IF
      IF(ID.EQ.'TIMESTEP') THEN
        THIS%DT=VAL
      ELSE IF(ID.EQ.'FRICTION') THEN
        THIS%FRICTION=VAL
      ELSE IF(ID.EQ.'MASS') THEN
        THIS%Q=VAL
      ELSE IF(ID.EQ.'TARGET') THEN
        THIS%EKIN_TARGET=VAL
      ELSE IF(ID.EQ.'EKIN(SYSTEM)') THEN
        THIS%EKIN_ACTUAL=VAL
      ELSE IF(ID.EQ.'GFREE') THEN
        THIS%GFREE=VAL
      ELSE IF(ID.EQ.'XNOS0') THEN
        THIS%X0=VAL
      ELSE IF(ID.EQ.'XNOSM') THEN
        THIS%XM=VAL
      ELSE IF(ID.EQ.'XNOSMM') THEN
        THIS%XMM=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('THERMOSTAT$SETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$SETL4(ID,VAL)
!     **************************************************************************
!     **  SET DATA                                                            **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$SETL4')
      END IF
      IF(ID.EQ.'ON') THEN
         THIS%ON=VAL
      ELSE IF(ID.EQ.'STOP') THEN
         THIS%STOP=VAL
      ELSE IF(ID.EQ.'COOLONLY') THEN
         THIS%TWAVE=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('THERMOSTAT$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$GETR8(ID,VAL)
!     **************************************************************************
!     **  SET DATA                                                            **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$GETR8')
      END IF
      IF(ID.EQ.'COOLING') THEN
         VAL=THIS%COOLING
      ELSE IF(ID.EQ.'ENERGY') THEN
        VAL=THIS%EPOT+THIS%EKIN
      ELSE IF(ID.EQ.'EKIN') THEN
        VAL=THIS%EKIN
      ELSE IF(ID.EQ.'EPOT') THEN
        THIS%EPOT   = 2.D0*THIS%EKIN_TARGET*THIS%X0
        IF(THIS%TWAVE) THIS%EPOT=0.D0
        VAL=THIS%EPOT
      ELSE IF(ID.EQ.'EDISS') THEN
        VAL=THIS%EDISS
      ELSE IF(ID.EQ.'TARGET') THEN
        VAL=THIS%EKIN_TARGET
      ELSE IF(ID.EQ.'XNOS0') THEN
        VAL=THIS%X0
      ELSE IF(ID.EQ.'XNOSM') THEN
        VAL=THIS%XM
      ELSE IF(ID.EQ.'XNOSMM') THEN
        VAL=THIS%XMM
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('THERMOSTAT$GETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE THERMOSTAT$GETL4(ID,VAL)
!     **************************************************************************
!     **  SET DATA                                                            **
!     **************************************************************************
      USE THERMOSTAT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO THERMOSTAT SELECTED')
        CALL ERROR$STOP('THERMOSTAT$GETL4')
      END IF
      IF(ID.EQ.'ON') THEN
         VAL=THIS%ON
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('THERMOSTAT$GETL4')
      END IF
      RETURN
      END
