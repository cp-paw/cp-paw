!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAW
!     **************************************************************************
!     **  THIS IS THE MAIN PAW SUBROUTINE                                     **
!     **    1) READ AND INITIALIZE                                            **
!     **    2) ITERATE                                                        **
!     **    3) CLOSE DOWN                                                     **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      LOGICAL(4)   :: TNWSTR
      LOGICAL(4)   :: TSTOP
      LOGICAL(4)   :: TFIRST
      LOGICAL(4)   :: TLAST
      LOGICAL(4)   :: TPRINT
      LOGICAL(4)   :: TMERMN
      INTEGER(4)   :: NFILO     ! FILE UNIT OF PROTOCOL
      INTEGER(4)   :: NFI0      ! TIME STEP COUNTER AT START
      INTEGER(4)   :: NFI       ! TIME STEP COUNTER
      INTEGER(4)   :: NOMORE    ! 
      INTEGER(4)   :: IPRINT    ! 
      INTEGER(4)   :: NBEG      ! DECIDES IF RESTART FILE IS READ
      REAL(8)      :: DELT
      LOGICAL(4)   :: TCHK
!     **************************************************************************
                              CALL TRACE$PUSH('PAW')
                              CALL TIMING$CLOCKON('INITIALIZATION')
      CALL STOPIT$SETSTARTTIME
!
!     ==================================================================
!     ====  READ CONTROL INPUT DATA FILE "CNTL"                     ====
!     ==================================================================
      CALL READIN(NBEG,NOMORE,IPRINT,DELT,TMERMN,TNWSTR)
!
!     ==================================================================
!     ==  READ STRUCTURAL DATA FROM FILE "STRC"                       ==
!     ==================================================================
      CALL STRCIN

!     ==================================================================
!     ==  SET DIMER DIMENSION                                         ==
!     ==================================================================
      CALL DIMER$INITDIM()
!
!     ==================================================================
!     ==  INITIALIZE ATOMS OBJECT                                     ==
!     ==================================================================
      CALL ATOMS$INITIALIZE
!
!     ==================================================================
!     ==  GENERATE G-VECTORS                                          ==
!     ==================================================================
                              CALL TIMING$CLOCKON('WAVES$INITIALIZE')
      CALL WAVES$INITIALIZE
                              CALL TIMING$CLOCKOFF('WAVES$INITIALIZE')
!
!     ==================================================================
!     ==  READ RESTART FILE (WAVE FUNCTION COEFFICIENTS, ETC)         ==
!     ==================================================================
      IF(NBEG.GE.0) THEN
                              CALL TIMING$CLOCKON('RESTART I/O')
        CALL READRESTART
                              CALL TIMING$CLOCKOFF('RESTART I/O')
      END IF
!
!     ================================================================
!     ==  REPORT INPUT DATA                                         ==
!     ================================================================
      CALL IO$REPORT
!
!     ==================================================================
!     ==  RUN TIME STATISTICS FOR THE INITIALIZATION                  ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
                              CALL TIMING$CLOCKOFF('INITIALIZATION')
                              CALL TIMING$PRINT('MONOMER',NFILO)
!     == CLOCK WILL BE RESTARTED AGAIN AT THE END OF THE FIRST ITERATION TO
!     == ACCOUNT FOR THE FACT THAT MANY ROUTINES INIITIALIZE THEMSELFES IN THE 
!     == FIRST ITERATION. THE TIME TO THIS RESET IS NOT PRINTED 
                              CALL TIMING$START
!
!     ==================================================================
!     ==  RETURN IF NO ITERATIONS ARE REQUIRED                        ==
!     ==================================================================
      IF(NOMORE.EQ.0) THEN
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==================================================================
!     ==================================================================
!     == THE BASIC LOOP FOR MOLECULAR DYNAMICS STARTS HERE            ==
!     ==================================================================
!     ==================================================================
      CALL TIMESTEP$GETI4('ISTEP',NFI)
      NFI0=NFI
      TFIRST=.TRUE.
      TLAST=.FALSE.
      TSTOP=.FALSE.
1000  CONTINUE
!
!     ==================================================================
!     ==  ITERATION CONTROL (PROPER STOP ETC. )                       ==
!     ==================================================================
      NFI=NFI+1
      CALL STOPIT$GETL4('STOP',TSTOP)
      IF(TSTOP) TLAST=.TRUE.
      IF(NFI.GE.NFI0+NOMORE) TLAST=.TRUE.
      TPRINT=(MOD(NFI,IPRINT).EQ.0.OR.TFIRST.OR.TLAST)

      CALL HYPERFINE$SETL4('WAKE',TPRINT)
      CALL GRAPHICS$SETL4('WAKE',TPRINT.AND.TLAST)
      CALL OPTEELS$SETL4('ON',TLAST)
      CALL CORE$SETL4('ON',TLAST)
!
!     ==================================================================
!     ==   WRITE RESTART_OUT                                          ==
!     ==================================================================
      IF(TLAST.OR.(TPRINT.AND.(.NOT.TFIRST))) THEN
                              CALL TIMING$CLOCKON('RESTART I/O') 
        CALL WRITERESTART
        CALL WAVES$GETL4('WRITERHO',TCHK)
        IF(TCHK) CALL WAVES$FIXRHOWRITE()
                              CALL TIMING$CLOCKOFF('RESTART I/O')
!       CALL MM_PAW_WRITE_RESTART (NFI) ! CALGARY QM/MM IMPLEMENTATION
      END IF
!
!     ==================================================================
!     ==   PERFORM ONE TIME STEP                                      ==
!     ==================================================================
                              CALL TIMING$CLOCKON('TIMESTEP')
!     ==USE THE LINE WITH "NOT TSTOP" TO AVOID AN ADDITIONAL LAST TIME STEP
!     IF(.NOT.TSTOP) CALL TIMESTEP(DELT,TPRINT,NFI,TSTOP)
      CALL DYNOCC$GETL4('DYN',TMERMN)
      CALL TIMESTEP(DELT,TPRINT,NFI,TSTOP)
                              CALL TIMING$CLOCKOFF('TIMESTEP')
!
!     ==================================================================
!     ==   WRITE INFORMATION AND TRAJECTORIES                         ==
!     ==================================================================
!     __ADD TO TRAJECTORIES (TEMPORARY BUFFER)__________________________
      IF(.NOT.TLAST) THEN
        CALL WRITETRAJECTORY(NFI,DELT)
      END IF
!     __ WRITE TRAJECTORY FROM TEMPORARY BUFFER TO FILE_________________
      IF(TPRINT.OR.TLAST) THEN
        CALL TRAJECTORYIO$FLUSHALL
      END IF
!
!     ==================================================================
!     == RESET CLOCKS FOR TIMING                                      ==
!     == THIS IS DONE AFTER THE FIRST ITERATION TO AVOID COUNTING     ==
!     == THE TIME FOR SELF-INITIALIZATION OF OBJECTS                  ==
!     ==================================================================
      IF(TFIRST) THEN
        CALL TIMING$START
      ELSE
        CALL TIMING$COUNT
      END IF
!
!     ==================================================================
!     ==   STOP OR CONTINUE LOOP                                      ==
!     ==================================================================
      TFIRST=.FALSE.
      IF(.NOT.TLAST) GO TO 1000
!     ==================================================================
!     ==================================================================
!     ====   END OF TIME-STEP ITERATION                               ==
!     ==================================================================
!     ==================================================================
      CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      MODULE TIMESTEP_MODULE
      INTEGER(4) :: ISTEPNUMBER
      REAL(8)    :: DELTAT
      LOGICAL(4) :: TNEWTHERMOSTAT=.TRUE. !SWITCHES WAVES THERMOSTAT
      END MODULE TIMESTEP_MODULE
!
!     ..................................................................
      SUBROUTINE TIMESTEP$SETR8(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'DELTAT') THEN
       DELTAT=VAL_
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$GETR8(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'DELTAT') THEN
        VAL_=DELTAT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$GETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$SETI4(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'ISTEP') THEN
        ISTEPNUMBER=VAL_
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$GETI4(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'ISTEP') THEN
        VAL_=ISTEPNUMBER
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$SETL4(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'NEWTHERMOSTAT') THEN
        TNEWTHERMOSTAT=VAL_
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP$GETL4(ID_,VAL_)
      USE TIMESTEP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'NEWTHERMOSTAT') THEN
        VAL_=TNEWTHERMOSTAT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$GETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP(DELT,TPRINT,NFI,TSTOP)
!     ******************************************************************      
!     ******************************************************************      
      USE TIMESTEP_MODULE ,ONLY : DELTAT,ISTEPNUMBER,TNEWTHERMOSTAT
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: DELT   ! TIME STEP
      LOGICAL(4),INTENT(IN)   :: TPRINT ! FLAG FOR LONG PRINTOUT
      LOGICAL(4),INTENT(IN)   :: TSTOP  ! FLAG FOR LAST TIME STEP
      INTEGER(4),INTENT(INOUT):: NFI    ! TIME STEP COUNTER
      INTEGER(4)              :: NFILO
      LOGICAL(4)              :: TFOR   ! ON/OFF SWITCH FOR ATOMIC MOTION
      LOGICAL(4)              :: TGRA   ! ON/OFF SWITCH FOR GRADIENT CORRECTION
      REAL(8)                 :: ANNEE  ! FRICTION COEFFICIENT FOR WAVE FUNCTIONS
      REAL(8)                 :: ANNER  ! FRICTION COEFFICIENT FOR NUCLEI
      LOGICAL(4)              :: TCHK
      REAL(8)                 :: TEMPINST
      REAL(8)                 :: ENOSE  ! GENERIC THERMOSTAT ENERGY 
      REAL(8)                 :: EKIN   ! GENERIC KINETIC ENERGY
      REAL(8)                 :: SVAR
      LOGICAL(4)              :: TCHK1,TCHK2
      REAL(8)                 :: FAV,FMAX
!     ******************************************************************      
                              CALL TRACE$PUSH('TIMESTEP')
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!      
      DELTAT=DELT       !-> TIMESTEP_MODULE
      ISTEPNUMBER=NFI   !-> TIMESTEP_MODULE
      CALL ENERGYLIST$RESET
!
      CALL ATOMS$GETL4('MOVE',TFOR)
      CALL DFT$GETL4('GC',TGRA)
      CALL WAVES$SETL4('HAMILTON',TPRINT)
!
!     ==================================================================
!     ==  COMMUNICATE TIME STEP TO OBJECTS                            ==
!     ==================================================================
      CALL WAVES$SETR8('TIMESTEP',DELT)
      CALL ATOMS$SETR8('TIMESTEP',DELT)
      CALL CELL$SETR8('DT',DELT)
      CALL DYNOCC$SETR8('TIMESTEP',DELT)
      CALL QMMM$SETR8('TIMESTEP',DELT)
      CALL THERMOSTAT$SELECT('ATOMS') 
      CALL THERMOSTAT$SETR8('TIMESTEP',DELT)
      CALL THERMOSTAT$SELECT('WAVES') 
      CALL THERMOSTAT$SETR8('TIMESTEP',DELT)
      CALL DIALS$SETR8('TIMESTEP',DELT)
!
!     ==================================================================
!     ==================================================================
!     ==  CALCULATE TOTAL ENERGY AND FORCES                           ==
!     ==================================================================
!     ==================================================================
!     == LDA TOTAL ENERGY ==============================================
      CALL WAVES$ETOT
!
!     == EXTERNAL POTENTIAL ACTING ON ATOMS ============================
      CALL VEXT$APPLY
!
!     == OCCUPATIONS ===================================================
      CALL DYNOCC$GETR8('EPOT',SVAR)
      CALL ENERGYLIST$SET('OCCUPATIONAL ENTROPY TERM (-TS)',SVAR)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',SVAR)
!
!     == UNIT CELL ========= ===================================================
      CALL CELL$ETOT()
      CALL CELL$GETR8('EPOT',SVAR)
      CALL ENERGYLIST$SET('CELLOSTAT POTENTIAL',SVAR)     
      CALL ENERGYLIST$ADD('TOTAL ENERGY',SVAR)
!
!     ==================================================================
!     ==================================================================
!     ==  SET FRICTION FOR ANNEALING SCHEDULE                         ==
!     ==================================================================
!     ==================================================================
      CALL AUTOPI(TCHK)
      IF(.NOT.TSTOP.AND.TCHK) THEN
        CALL STOPIT$SETL4('STOP',TCHK)
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)'STOP SIGNAL FROM AUTOPILOT'
      END IF
! 
!     ==================================================================
!     ==  OBTAIN INSTANTANEOUS FRICTION VALUES FROM THERMOSTATS       ==
!     ==  AND COMMUNICATE TO WAVES AND ATOMS OBJECTS                  ==
!     ==================================================================
!
!     == SET THERMOSTAT FRICTION FOR ATOMS =============================
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$GETL4('ON',TCHK1)
      IF(TCHK1) THEN
        CALL THERMOSTAT$GETR8('COOLING',ANNER)
        CALL ATOMS$SETR8('FRICTION',ANNER)
      END IF
!
!     == SET THERMOSTAT FRICTION FOR WAVE FUNCTIONS ====================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$GETL4('ON',TCHK2)
      IF(TCHK2) THEN
        CALL THERMOSTAT$SELECT('WAVES')
        CALL THERMOSTAT$GETR8('COOLING',ANNEE)
        CALL WAVES$SETR8('FRICTION',ANNEE)
        CALL ATOMS$SETR8('ANNEE',ANNEE)
      ELSE
        CALL ATOMS$SETR8('ANNEE',0.D0)
      END IF
!!$!
!!$!     == DETERMINE IF CORRECTION FOR WAVE FUNCTION DRAG IS ON ==========
!!$      IF(TCHK1) THEN          ! ATOM THERMOSTAT OR BOTH THERMOSTATS
!!$        CALL WAVES$GETR8('FRICTION',ANNEE)
!!$        CALL ATOMS$SETR8('ANNEE',ANNEE)
!!$        IF(.NOT.TNEWTHERMOSTAT) THEN
!!$          CALL ATOMS$SETR8('ANNEE',0.D0)
!!$        END IF
!!$      ELSE ! FRICTION DYNAMICS (NO THERMOSTAT)
!!$        IF(TCHK2) THEN ! WAVE FUNCTION THERMOSTAT ALONE NOT ALLOWED
!!$          CALL ERROR$MSG('WAVE FUNCTION THERMOSTAT ALONE IS NOT ALLOWED')
!!$          CALL ERROR$STOP('TIMESTEP')
!!$        ELSE 
!!$          CALL ATOMS$GETR8('FRICTION',ANNER)
!!$          CALL ATOMS$SETR8('ANNEE',ANNER)
!!$        END IF
!!$      END IF
!
!     ==================================================================
!     ==================================================================
!     ==  PROPAGATE:                                                  ==
!     ==================================================================
!     ==================================================================
! 
!     ==================================================================
!     ==  PROPAGATE NUCLEI (CONSTRAINTS FOLLLOW LATER...)             ==
!     ==================================================================
!---DIMERMERGE FIX      
      CALL DIMER$GETL4('DIMER',TCHK)
      IF(TCHK) THEN
         CALL ATOMS$PROPAGATE_DIMER()
      ELSE
         CALL ATOMS$PROPAGATE()
      END IF
!---END DIMERMERGE FIX      
! 
!     ==================================================================
!     ==  PROPAGATE UNIT CELL                                         ==
!     ==================================================================
      CALL CELL$PROPAGATE()
! 
!     ==================================================================
!     ==  APPLY CONSTRAINTS TO ATOMIC COORDINATES                     ==
!     ==================================================================
      CALL ATOMS$CONSTRAINTS()
! 
!     ==================================================================
!     ==  PROPAGATE WAVE FUNCTIONS                                    ==
!     ==================================================================
      CALL WAVES$PROPAGATE()
      CALL WAVES$ORTHOGONALIZE()
!
!     ==================================================================
!     ==  PROPAGATE OCCUPATIONS                                       ==
!     ==================================================================
      CALL DYNOCC$PROPAGATE()
! 
!     ==================================================================
!     ==  PROPAGATE THERMOSTAT FOR THE NUCLEI                         ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
!       __EXTRACT KINETIC ENERGY (USING THE TRUE, EFFECTIVE MASSES)_____
        CALL ATOMS$EKIN(EKIN)
        CALL THERMOSTAT$SETR8('EKIN(SYSTEM)',EKIN)
!       __PROPAGATE THERMOSTAT AND RECORD ITS ENERGY____________________
        CALL THERMOSTAT$PROPAGATE()
      END IF
!
!     ==================================================================
!     ==  PROPAGATE THERMOSTAT FOR THE WAVE FUNCTIONS                 ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
        IF(TNEWTHERMOSTAT) THEN
          CALL ATOMS$EFFEKIN(EKIN)
          CALL THERMOSTAT$SETR8('TARGET',EKIN)
        END IF
        CALL WAVES$GETR8('EKIN(PSI)',EKIN)
        CALL THERMOSTAT$SETR8('EKIN(SYSTEM)',EKIN)
        CALL THERMOSTAT$PROPAGATE()
      END IF
!
!     ==================================================================
!     ==================================================================
!     ==  COLLECT KINETIC ENERGIES                                    ==
!     ==================================================================
!     ==================================================================
!
!     ==================================================================
!     __WAVE FUNCTIONS__________________________________________________
!     ==================================================================
      CALL WAVES$GETR8('EKIN(PSI)',EKIN)
      CALL ENERGYLIST$SET('WAVEFUNCTION KINETIC ENERGY',EKIN)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',EKIN)
!
!     ==================================================================
!     __ NUCLEI__________________________________________________________
!     ==================================================================
!     __CALCULATE KINETIC ENERGY OF THE ATOMS_ _________________________
      CALL ATOMS$EKIN(EKIN)
      CALL ENERGYLIST$SET('IONIC KINETIC ENERGY',EKIN)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',EKIN)
!     __SUBTRACT THE EFFECTIVE WAVE FUNCTION KINETIC ENERGY_____________
      CALL ATOMS$EFFEKIN(EKIN)
      CALL ENERGYLIST$SET('BO-WAVEFUNCTION KINETIC ENERGY',EKIN)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',-EKIN)
!     __STORE INSTANTANEOUS ATOMIC TEMPERATURE__________________________
      CALL ATOMS$TEMPERATURE(TEMPINST)
      CALL ENERGYLIST$SET('IONIC TEMPERATURE',TEMPINST)
!
!     ==================================================================
!     == UNIT CELL                                                    ==
!     ==================================================================
      CALL CELL$GETR8('EKIN',EKIN)
      CALL ENERGYLIST$SET('CELLOSTAT KINETIC',EKIN)     
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',EKIN)
!
!     ==================================================================
!     == OCCUPATIONS                                                  ==
!     ==================================================================
      CALL DYNOCC$GETR8('EKIN',SVAR)
      CALL ENERGYLIST$SET('OCCUPATION KINETIC ENERGY',SVAR)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',SVAR)
! 
!     ==================================================================
!     ==  THERMOSTAT ACTING ON THE NUCLEI                             ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
        CALL THERMOSTAT$GETR8('ENERGY',ENOSE)
        CALL ENERGYLIST$SET('ATOM THERMOSTAT',ENOSE)
        CALL ENERGYLIST$ADD('CONSTANT ENERGY',ENOSE)
      END IF
!
!     ==================================================================
!     ==  COLLECT ENERGY OF WAVE FUNCTION THERMOSTAT                  ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
        IF(TNEWTHERMOSTAT) THEN
          CALL THERMOSTAT$GETR8('EKIN',ENOSE)
        ELSE
          CALL THERMOSTAT$GETR8('ENERGY',ENOSE)
        END IF
        CALL ENERGYLIST$SET('ELECTRON THERMOSTAT',ENOSE)
        CALL ENERGYLIST$ADD('CONSTANT ENERGY',ENOSE)
      END IF
!     EFLUXE=-EKINC*(2.D0*ANNEE/DELT)/EMASS
!
!     ==================================================================
!     ==================================================================
!     == CHECK CONVERGENCE                                            ==
!     ==================================================================
!     ==================================================================
      CALL ATOMS$FORCECRITERION(FAV,FMAX)
!PRINT*,'FORCE FAV=',FAV,' FMAX=',FMAX
!
!     ==================================================================
!     ==================================================================
!     == WRITE INFORMATION                                            ==
!     ==================================================================
!     ==================================================================
!     __WRITE PROTOCOL__________________________________________________
      CALL PRINFO(TPRINT,TSTOP,NFI,DELT)
!
!     ==================================================================
!     ==================================================================
!     == UPDATE DYNAMICAL VARIABLES                                   ==
!     ==================================================================
!     ==================================================================
!     __ELECTRONIC WAVE FUNCTIONS_______________________________________
      CALL WAVES$SWITCH()
!     __ATOMIC POSITIONS________________________________________________
      CALL ATOMS$SWITCH()
!     __UNIT CELL_______________________________________________________________
!     __IMPORTANT! _____________________________________________________________
!     __CELL$SWITCH MUST BE CALLED AFTER WAVES$SWITCH AND ATOMS$SWITCH__________
      CALL CELL$SWITCH()
!     __OCCUPATIONS_____________________________________________________
      CALL DYNOCC$SWITCH() 
!     __WAVE FUNCTION THERMOSTAT________________________________________
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$SWITCH()
!     __ATOM THERMOSTAT_________________________________________________
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$SWITCH()
!
!     ==================================================================
!     == TURN DIALS                                                   ==
!     ==================================================================
      CALL DIALS$APPLY
                              CALL TRACE$POP()
      RETURN
      END SUBROUTINE TIMESTEP
!
!     ...................................................STOPIT.........
MODULE STOPIT_MODULE
!**                                                                   ** 
!**  SET SWITCH TO INITIATE SOFT-KILL                                 ** 
!**  SWITCH CAN BE SET BY                                             ** 
!**   A) CREATE A PREDEFINED EXITFILE                                 ** 
!**   B) BY THE CODE BY SETTING TSTOP EXPLICITELY                     ** 
!**                                                                   ** 
USE CLOCK_MODULE 
!USE CLOCK_MODULE , ONLY : DATE_TIME,(.LATER.)
LOGICAL(4),save :: TSTOP=.FALSE.   ! INITIATE SOFT-KILL
LOGICAL(4)      :: DISTRIBUTED=.FALSE. ! TSTOP=TRUE ON ALL NODES
LOGICAL(4)      :: TNOTIFY=.TRUE.  ! NOTIFICATION ABOUT STOP REQUIRED
LOGICAL(4)      :: EXITFILEREMOVED=.FALSE.
INTEGER(4)      :: RUNTIME=-1         ! RUNTIME IN SECONDS
TYPE(DATE_TIME) :: STARTTIME
END MODULE STOPIT_MODULE
!
!     ..................................................................
      SUBROUTINE STOPIT$SETL4(ID_,VAL_)
!     ****************************************************************** 
!     **  STOPIT$SET                                                  ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: ID_
      LOGICAL(4)  ,INTENT(IN) :: VAL_
!     ****************************************************************** 
      IF(TRIM(ID_).EQ.'STOP') THEN
        TSTOP=VAL_
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('STOPIT$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STOPIT$SETI4(ID_,VAL_)
!     ****************************************************************** 
!     **  STOPIT$SET                                                  ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: ID_
      INTEGER(4)  ,INTENT(IN) :: VAL_
!     ****************************************************************** 
      IF(TRIM(ID_).EQ.'RUNTIME') THEN 
        RUNTIME=VAL_    ! RUNTIME IN SECONDS
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('STOPIT$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE STOPIT$GETL4(ID_,VAL_)
!     ****************************************************************** 
!     **  STOPIT$GET                                                  ** 
!     **  USE STOPIT$GET('TSTOP',4,TSTOP) TO SEE IF STOP IS SET       ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(OUT):: VAL_
!     ****************************************************************** 
      IF(TRIM(ID_).EQ.'STOP') THEN
        CALL  STOPIT_UPDATE
        VAL_=TSTOP
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('STOPIT$GETL4')
      END IF
      RETURN
      END
!
!     ...................................................STOPIT.........
      SUBROUTINE STOPIT$SETSTARTTIME
!     ****************************************************************** 
!     **  STOPIT_UPDATE                                               ** 
!     **                                                              ** 
!     **  TESTS WHETHER TSTOP HAS BECOME TRUE IN THE MEANWHILE        ** 
!     **                                                              ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      IMPLICIT NONE
      CALL CLOCK$NOW(STARTTIME)
      RETURN
      END

!
!     ...................................................STOPIT.........
      SUBROUTINE STOPIT_UPDATE
!     ****************************************************************** 
!     **  STOPIT_UPDATE                                               ** 
!     **                                                              ** 
!     **  TESTS WHETHER TSTOP HAS BECOME TRUE IN THE MEANWHILE        ** 
!     **                                                              ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      LOGICAL(4)     :: TCHK
      INTEGER(4)     :: NVAL
      INTEGER(4)     :: NFILO
      CHARACTER(256) :: EXITFILE=' '
      CHARACTER(264) :: CMD
      INTEGER(4)     :: NTASKS,THISTASK
      TYPE(DATE_TIME) :: ENDTIME,NOW
      INTEGER(4)      :: ISVAR
!     ****************************************************************** 
!
!     ==================================================================
!     ==  REMOVE EXIT FILE IN THE FIRST REQUEST                       ==
!     ==================================================================
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      IF(THISTASK.EQ.1) THEN
        IF(.NOT.EXITFILEREMOVED) THEN
          CALL FILEHANDLER$DELETE('EXIT')
          EXITFILEREMOVED=.TRUE.
        END IF
      ELSE
        EXITFILEREMOVED=.TRUE.
      END IF

!     ==================================================================
!     ==  CHECK IF TIME EXCEEDS LIMIT                                 ==
!     ==================================================================
      IF(RUNTIME.GT.0) THEN
        ENDTIME=STARTTIME
        ISVAR=RUNTIME
        ENDTIME%SECOND=ENDTIME%SECOND+ISVAR
!
        ISVAR=INT(ENDTIME%SECOND/60)
        ENDTIME%SECOND=ENDTIME%SECOND-ISVAR*60
        ENDTIME%MINUTE=ENDTIME%MINUTE+ISVAR
!
        ISVAR=INT(ENDTIME%MINUTE/60)
        ENDTIME%MINUTE=ENDTIME%MINUTE-ISVAR*60
        ENDTIME%HOUR=ENDTIME%HOUR+ISVAR
!
        ISVAR=INT(ENDTIME%HOUR/24)
        ENDTIME%HOUR=ENDTIME%HOUR-ISVAR*24
        ENDTIME%DAY=ENDTIME%DAY+ISVAR
!
!       WARNING! THIS CHOICE STOPS AT THE LAST SECOND OF THE CURRENT MONTH!!!!!!!
        CALL CLOCK$NOW(NOW)
        IF(NOW.LATER.ENDTIME) TSTOP=.TRUE.
      END IF
!
!     ==================================================================
!     ==  CHECK IF EXITFILE EXISTS                                    ==
!     ==================================================================
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      IF(.NOT.TSTOP.AND.THISTASK.EQ.1) THEN
        CALL FILEHANDLER$FILENAME('EXIT',EXITFILE)
        INQUIRE(FILE=EXITFILE,EXIST=TSTOP)
      END IF
!
!     ==================================================================
!     ==  CHECK WHETHER TSTOP=T FOR ANY TASK                          ==
!     ==================================================================
      IF(.NOT.DISTRIBUTED) THEN
        IF(TSTOP) THEN 
          NVAL=1
        ELSE
          NVAL=0
        ENDIF
        CALL MPE$COMBINE('~','+',NVAL)
        TSTOP=(NVAL.NE.0) 
        DISTRIBUTED=TSTOP
      END IF
!
!     ==================================================================
!     == WRITE MESSAGE IF NOT DONE ALREADY =============================
!     ==================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(TSTOP.AND.TNOTIFY.AND.THISTASK.EQ.1) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)'STOP SIGNAL RECEIVED'
        FLUSH(NFILO)
        TNOTIFY=.FALSE.
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PRINFO(TPRINT,TLAST,NFI,DELT)
!     **************************************************************************
!     **  REPORTS ON THE PROCESS OF THE SIMULATION AND INVOKES                **
!     **  ANALYSIS ROUTINES                                                   **
!     **                                                                      **
!     **  PRINFO IS CALLED ONCE PER TIMESTEP                                  **
!     **************************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TPRINT  ! FLAG FOR EXTENSIVE PRINTOUT
      LOGICAL(4),INTENT(IN) :: TLAST   ! FLAG FOR ONE-TIME EXECUTION
      INTEGER(4),INTENT(IN) :: NFI     ! TIME STEP NUMBER
      REAL(8)   ,INTENT(IN) :: DELT    ! RIME STEP
      LOGICAL(4),SAVE       :: TFIRST=.TRUE.
      INTEGER(4)            :: NAT
      INTEGER(4)            :: NFILO
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: NTASKS,THISTASK
      REAL(8)               :: PICO
      REAL(8)               :: SECOND
      REAL(8)               :: TIME     ! ACTUAL TIME
      REAL(8)               :: TME1
      REAL(8)               :: ECONS
      REAL(8)               :: ETOT
      REAL(8)               :: EKINP
      REAL(8)               :: EKINC
      REAL(8)               :: ENOSEE
      REAL(8)               :: ENOSEP
      REAL(8)               :: EKINFC
      REAL(8)               :: HEAT
      REAL(8)               :: EEXT
      REAL(8)               :: OCCKIN
      REAL(8)               :: TEMPINST
      REAL(8)               :: CELVIN
      INTEGER(4)            :: ITEMP
      REAL(8)               :: ANNEE
      REAL(8)               :: ANNER
      REAL(8)               :: SVAR
      INTEGER(4)            :: ISVAR
      REAL(8)               :: EFFEKIN
      REAL(8)               :: ECELLPOT
      REAL(8)               :: ECELLKIN
      LOGICAL(4)            :: TCHK,TCHK1
      LOGICAL(4)            :: TQMMM=.FALSE.
      REAL(8)               :: QMMMKIN   ! EKIN OF QM-MM ENVIRONMENT
      REAL(8)               :: QMMMPOT   ! EPOT OF QM-MM ENVIRONMENT
      REAL(8)               :: QMMMTHERM ! ENERGY OF THE QM-MM THERMOSTAT
      REAL(8)               :: QMMMTEMP  ! TEMPERATURE OF THE QM-MM 
      LOGICAL(4)            :: CALGARY_QMMM
      REAL(8)               :: MM_KINETIC_ENERGY
      REAL(8)               :: MM_POT_ENERGY
      REAL(8)               :: MM_TEMP
      INTEGER(4)            :: IMM_TEMP
      REAL(8)               :: MM_NOSE_ENERGY
      REAL(8)               :: MM_FRIC
      LOGICAL(4)            :: TCOSMO
      REAL(8)               :: EKINCOSMO,EPOTCOSMO
      CHARACTER(8)          :: FFTYPE
!     **************************************************************************
                              CALL TRACE$PUSH('PRINFO')
      TIME=DBLE(NFI)*DELT
      CALL ATOMLIST$NATOM(NAT)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!   
!     ==========================================================================
!     ==   HYPERFINE PARAMETERS                                               ==
!     ==========================================================================
      IF(TPRINT) CALL HYPERFINE$PRINT
!   
!     ==========================================================================
!     ==  PLOT WAVE FUNCTIONS, DENSITIES                                      ==
!     ==========================================================================
!     == GRAPHICS OBJECT IS SET INTO THE WAKE STATE ONLY IF TLAST.AND.TPRINT
                             CALL TRACE$PASS('BEFORE GRAPHICS$PLOT')
      CALL GRAPHICS$PLOT()
                             CALL TRACE$PASS('AFTER GRAPHICS$PLOT')
!
!     ==========================================================================
      IF(TPRINT) THEN
!       == OPTEELS OBJECT IS SWITCHED ON ONLY IN THE LAST ITERATION ============
        CALL OPTEELS$PLOT()
      END IF
!   
!     ==========================================================================
!     ==   THE FOLLOWING IS ONLY EXECUTED ON THE FIRST NODE                   ==
!     ==========================================================================
!     IF(THISTASK.GT.1) THEN
!       CALL TRACE$POP
!       RETURN
!     END IF

!     CALL MM_RUN(CALGARY_QMMM)
      CALGARY_QMMM = .FALSE.
    
!     ==========================================================================
!     ==   WRITE HEADER FOR PROTOCOL FOR EACH TIME STEP                       ==
!     ==========================================================================
                             CALL TRACE$PASS('WRITE HEADER')
      IF(THISTASK.EQ.1.AND.TFIRST) THEN
        CALL COSMO$GETL4('ON',TCOSMO)
        IF(TQMMM) THEN
          WRITE(NFILO,FMT='()')
          WRITE(NFILO,FMT='(2X,A5,A9,1X,A4,1X,A10,2A13,2A6,3A10)') &
     &            'NFI','T[PSEC]','T[K]','EKIN(PSI)','E(RHO)','ECONS', &
     &            'ANNEE','ANNER','T(ENV)','EPOT(ENV)','ECNS(ENV)'
        ELSE IF (CALGARY_QMMM) THEN
          WRITE(NFILO,FMT='(/5X,"NFI",3X,"TIME",1X,"TEMP",3X,"EKINC"'   &
     &                  //',5X,"E(RHO)",6X,"ECONS",1X,"ANNEE",1X,"ANNER",' &
     &                  //'"   E_MM  TEMP FRIC")')
        ELSE IF (TCOSMO) THEN
          WRITE(NFILO,FMT='()')
          WRITE(NFILO,FMT='(2X,A5,A9,1X,A4,1X,A10,2A13,2A8,2A11)') &
     &            'NFI','T[PSEC]','T[K]','EKIN(PSI)','E(RHO)','ECONS', &
     &            'ANNEE','ANNER','COSMO-EPOT','COSMO-EKIN'
        ELSE
          WRITE(NFILO,FMT='()')
          WRITE(NFILO,FMT='(2X,A5,A9,1X,A4,1X,A11,2A14,2A8)') &
     &            'NFI','T[PSEC]','T[K]','EKIN(PSI)','E(RHO)','ECONS', &
     &            'ANNEE','ANNER'
        END IF
        TFIRST=.FALSE.
      END IF
      IF(TPRINT) TFIRST=.TRUE.
!   
!     ==========================================================================
!     ==   WRITE PROTOCOL FOR EACH TIME STEP                                  ==
!     ==========================================================================
                             CALL TRACE$PASS('WRITE PROTOCOL')
      IF(THISTASK.EQ.1) THEN

        CALL CONSTANTS('PICO',PICO)
        CALL CONSTANTS('SECOND',SECOND)
        CALL CONSTANTS('KB',CELVIN)
        TME1=TIME/(PICO*SECOND)
!
!       ========================================================================
!       == COMPARE TOTAL ENERGY FROM THE ENERGYLIST WITH THE SUM OF INDIVIDUAL==
!       == ENERGY CONTRIBUTIONS                                               ==
!       ========================================================================
        ETOT=0.D0
!       __PAW_WAVES1.F90________________________________________________________
        CALL ENERGYLIST$GET('PS  KINETIC',SVAR)
        ETOT=ETOT+SVAR
!       __PAW_AUGMENTATION.F90__________________________________________________
        CALL ENERGYLIST$GET('AE1-PS1 KINETIC',SVAR)
        ETOT=ETOT+SVAR
        CALL ENERGYLIST$GET('AE1 EXCHANGE-CORRELATION',SVAR)
        ETOT=ETOT+SVAR
        CALL ENERGYLIST$GET('PS1 EXCHANGE-CORRELATION',SVAR)
        ETOT=ETOT-SVAR
        CALL ENERGYLIST$GET('AE1 ELECTROSTATIC',SVAR)
        ETOT=ETOT+SVAR
        CALL ENERGYLIST$GET('PS1 ELECTROSTATIC',SVAR)
        ETOT=ETOT-SVAR
        CALL ENERGYLIST$GET('AE1 BACKGROUND',SVAR)
        ETOT=ETOT+SVAR
        CALL ENERGYLIST$GET('PS1 BACKGROUND',SVAR)
        ETOT=ETOT-SVAR
        CALL ENERGYLIST$GET('LDA+U EXCHANGE',SVAR)
        ETOT=ETOT+SVAR
        CALL ENERGYLIST$GET('CORE RELAXATION',SVAR)
        ETOT=ETOT+SVAR
        CALL ENERGYLIST$GET('EXTERNAL 1CENTER POTENTIAL',SVAR)
        ETOT=ETOT+SVAR
!       __PAW_POTENTIAL.F90_____________________________________________________
        CALL ENERGYLIST$GET('SOLVENT PAULI REPULSION',SVAR)
        ETOT=ETOT+SVAR
!       __PAW_POTENTIAL.F90_____________________________________________________
!       __NOTE: 'PAIRPOTENTIAL' IS PART OF 'PS  ELECTROSTATIC'__________________
        CALL ENERGYLIST$GET('PS  ELECTROSTATIC',SVAR)
        ETOT=ETOT+SVAR
!       __PAW_POTENTIAL.F90_____________________________________________________
        CALL ENERGYLIST$GET('PS  EXCHANGE-CORRELATION',SVAR)
        ETOT=ETOT+SVAR
!       __PAW_SIMPLELMTO.F90____________________________________________________
        CALL ENERGYLIST$GET('LOCAL CORRELATION',SVAR)
        ETOT=ETOT+SVAR
!       __SRC/PAW_DRIVER.F90____________________________________________________
        CALL ENERGYLIST$GET('OCCUPATIONAL ENTROPY TERM (-TS)',SVAR)
        ETOT=ETOT+SVAR
!       __SRC/PAW_DRIVER.F90____________________________________________________
        CALL ENERGYLIST$GET('CELLOSTAT POTENTIAL',SVAR)
        ETOT=ETOT+SVAR
!       __SRC/PAW_EXTPOT.F90____________________________________________________
        CALL ENERGYLIST$GET('EXTERNAL POTENTIAL',SVAR)
        ETOT=ETOT+SVAR
!       __SRC/PAW_ISOLATE.F90___________________________________________________
        CALL ENERGYLIST$GET('BACKGROUND ENERGY',SVAR)
        ETOT=ETOT+SVAR
!       __PAW_ISOLATE.F90_______________________________________________________
        CALL ENERGYLIST$GET('ISOLATE ENERGY',SVAR)
        ETOT=ETOT+SVAR
!       __PAW_ISOLATE.F90_______________________________________________________
        CALL ENERGYLIST$GET('QMMM POTENTIAL ENERGY',SVAR)
        ETOT=ETOT+SVAR
!       __PAW_ISOLATE.F90_______________________________________________________
        CALL ENERGYLIST$GET('COSMO POTENTIAL ENERGY',SVAR)
        ETOT=ETOT+SVAR
!       __PAW_ISOLATE.F90_______________________________________________________
        CALL ENERGYLIST$GET('VAN DER WAALS ENERGY',SVAR)
        ETOT=ETOT+SVAR
!
!!$SRC/PAW_DMFT.F90:      CALL ENERGYLIST$ADD('TOTAL ENERGY',ETOT)
!!$SRC/PAW_DMFT_BACKES.F90:      CALL ENERGYLIST$ADD('TOTAL ENERGY',ETOT)
!!$SRC/PAW_DMFT_PETER.F90:      CALL ENERGYLIST$ADD('TOTAL ENERGY',ETOT)
!!$SRC/PAW_LMTO.F90:      CALL ENERGYLIST$ADD('TOTAL ENERGY',ETOT)
!!$SRC/PAW_LMTO.F90:      CALL ENERGYLIST$ADD('TOTAL ENERGY',EXTOT)
!
        CALL ENERGYLIST$GET('TOTAL ENERGY',SVAR)
        IF(ABS(ETOT-SVAR).GT.1.D-7) THEN
          CALL ENERGYLIST$PRINT(6)
          CALL ERROR$MSG('INTERNAL CHECK FAILED')
          CALL ERROR$MSG('TOTAL ENERGY INCONSISTENT WITH INDIVIDUAL TERMS')
          CALL ERROR$R8VAL('TOTAL ENERGY FROM LIST',SVAR)
          CALL ERROR$R8VAL('SUM OF INDIVIDUAL TERMS',ETOT)
          CALL ERROR$STOP('PRINFO')
        END IF       
!
!       ========================================================================
!       == ADD KINETIC ENERGIES TO OBTAIN CONSERVED ENERGY ECONS              ==
!       ========================================================================
        CALL ENERGYLIST$GET('TOTAL ENERGY',ETOT)     
        ECONS=ETOT
!
!       == ATOMIC AND FICTITIOUS ELECTRONIC KINETIC ENERGY =====================
        CALL ENERGYLIST$GET('IONIC KINETIC ENERGY',EKINP)
        CALL ENERGYLIST$GET('WAVEFUNCTION KINETIC ENERGY',EKINC)     
        CALL ENERGYLIST$GET('BO-WAVEFUNCTION KINETIC ENERGY',EFFEKIN)
        ECONS=ECONS+EKINC-EFFEKIN+EKINP
!
!       == COSMO ===============================================================
        CALL ENERGYLIST$GET('COSMO KINETIC ENERGY',EKINCOSMO)
        ECONS=ECONS+EKINCOSMO
!
!       == CELLOSTAT ===========================================================
        CALL ENERGYLIST$GET('CELLOSTAT KINETIC',ECELLKIN)     
        ECONS=ECONS+ECELLKIN 
!
!       == ELECTRON AND ATOM THERMOSTATS =======================================
        CALL ENERGYLIST$GET('ATOM THERMOSTAT',ENOSEP)     
        CALL ENERGYLIST$GET('ELECTRON THERMOSTAT',ENOSEE)     
        ECONS=ECONS+ENOSEP+ENOSEE
!
        CALL ENERGYLIST$GET('CONSTRAINT KINETIC ENERGY',EKINFC)     
        ECONS=ECONS+EKINFC
!
!       == OCCUPATIONS =========================================================
        CALL ENERGYLIST$GET('OCCUPATION KINETIC ENERGY',OCCKIN)
        ECONS=ECONS+OCCKIN
!
!       == QM-MM ENVIRONMENT ===================================================
        CALL QMMM$GETL4('ON',TQMMM)
        IF(TQMMM) THEN
          CALL QMMM$GETR8('EKIN',QMMMKIN)
          CALL QMMM$GETR8('EPOT',QMMMPOT)
          ECONS=ECONS+QMMMKIN   !POTENTIAL ENERGY ALREADY INCLUDED IN ETOT
!
!         == THERMOSTAT FOR THE ENVIRONMENT ====================================
          CALL QMMM$GETR8('ETHERM',QMMMTHERM)
          ECONS=ECONS+QMMMTHERM
        END IF

!       == QMMM CALGARY IMPLEMENTATION  ========================================
        IF (CALGARY_QMMM) THEN
          CALL ENERGYLIST$GET('MM KINETIC ENERGY',MM_KINETIC_ENERGY)
          CALL ENERGYLIST$GET('MM POT ENERGY',MM_POT_ENERGY)
          CALL ENERGYLIST$GET('MM TEMPERATURE',MM_TEMP)
          CALL ENERGYLIST$GET('MM THERMOSTAT',MM_NOSE_ENERGY)
          CALL ENERGYLIST$GET('MM ATOM FRICTION',MM_FRIC)
          ECONS=ECONS + MM_KINETIC_ENERGY + MM_NOSE_ENERGY
        END IF
!
        CALL ENERGYLIST$GET('CONSTANT ENERGY',SVAR)  ! DO NO USE!
PRINT*,'CONSTANT ENERGY ',ECONS,SVAR
!
!       == SOME OTHER STUFF ====================================================
        CALL ENERGYLIST$GET('IONIC TEMPERATURE',TEMPINST)
        ITEMP=NINT(TEMPINST/CELVIN)
        CALL WAVES$GETR8('FRICTION',ANNEE)
        CALL ATOMS$GETR8('FRICTION',ANNER)
!
!       __REMARK: ALL DATA ENTRIES SHOULD BE SEPARATED BY A BLANK (1X),_________
!       __SO THAT THE PARSER PAW_SHOW.X DOES NOT GET CONFUSED IF VALUES_________
!       __ARE TOO LARGE. THIS HAS BEEN DONE SOFAR ONLY FOR THE GENERIC PRINTOUT_
!       __AND IT NEEDS TO BE DONE FOR TCOSMO, TQMMM, CALGARY_QMMM_______________
        CALL COSMO$GETL4('ON',TCOSMO)
        IF (TCOSMO) THEN
          CALL ENERGYLIST$GET('COSMO KINETIC ENERGY',EKINCOSMO)
          CALL ENERGYLIST$GET('COSMO POTENTIAL ENERGY',EPOTCOSMO)
          WRITE(NFILO,FMT='("!>",I6,F10.5,1X,I5,F10.5,2F13.5,2F8.5' &
     &                     //',2F11.5)') &
     &               NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER &
     &              ,EPOTCOSMO,EKINCOSMO
        ELSE IF(TQMMM) THEN
          CALL CONSTANTS('KB',CELVIN)
          CALL QMMM$GETI4('NAT:ENV',ISVAR)
          QMMMTEMP=2.D0*QMMMKIN/REAL(3*ISVAR,KIND=8)/CELVIN
          WRITE(NFILO,FMT='("!>",I5,F10.5,1X,I5,F10.5,2F13.5,2F6.3' &
     &                //',I10,2F10.5)') &
     &                NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER &
     &               ,NINT(QMMMTEMP),QMMMPOT,QMMMKIN+QMMMPOT
        ELSE IF(CALGARY_QMMM) THEN
          IMM_TEMP=NINT(MM_TEMP)
          WRITE(NFILO,FMT='("!>",I6,F8.4,1X,I4,F8.5,2F11.5,2F6.3,1X,F6.3,1X' &
     &       //',I4,1X,F5.2 )') &
     &                     NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER &
     &                    ,MM_POT_ENERGY, IMM_TEMP, MM_FRIC
        ELSE
          WRITE(NFILO,FMT='("!>",I6,1X,F9.5,1X,I5,1X,F10.6,1X,F13.6,1X,F13.6' &
     &                     //',1X,F7.4,1X,F7.4)') &
     &                NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER
        ENDIF
!
!       ========================================================================
!       ==   WRITE ENERGIES TO PROTOCOL                                       ==
!       ========================================================================
                             CALL TRACE$PASS('WRITE ENERGIES')
        IF(TPRINT) THEN
          CALL ENERGYLIST$PRINTHEADER(NFILO)     
!         ==  BASIC LDA ========================================================
          CALL ENERGYLIST$PRINTONE(NFILO,'TOTAL ENERGY')     
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  KINETIC')     
!         == AE ELECTROSTATIC ENERGY DOES NOT INCLUDE ISOLATE ENERGY  ==========
!         == (PREVIOUSLY IT WAS ADDED)                             =============
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  ELECTROSTATIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  EXCHANGE-CORRELATION')     
          CALL ENERGYLIST$PRINTONE(NFILO,'    LOCAL CORRELATIONS')     
          CALL ENERGYLIST$PRINTONE(NFILO,'BACKGROUND')     
          CALL ENERGYLIST$PRINTONE(NFILO,'ISOLATE ENERGY')     
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  KINETIC')     
!         == PS ELECTROSTATIC ENERGY DOES NOT INCLUDE ISOLATE ENERGY  ==========
!         == (PREVIOUSLY IT WAS ADDED)                             =============
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  ELECTROSTATIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  EXCHANGE-CORRELATION')     
          CALL ENERGYLIST$PRINTONE(NFILO,'CORE RELAXATION')     
!         == ATOM KINETIC ENERGY ===============================================
          CALL ENERGYLIST$PRINTONE(NFILO,'IONIC KINETIC ENERGY')
          CALL ENERGYLIST$PRINTONE(NFILO,'IONIC TEMPERATURE')
          CALL ENERGYLIST$PRINTONE(NFILO,'ATOM THERMOSTAT')     
!
          CALL ENERGYLIST$PRINTONE(NFILO,'WAVEFUNCTION KINETIC ENERGY')     
          CALL ENERGYLIST$PRINTONE(NFILO,'BO-WAVEFUNCTION KINETIC ENERGY')
          CALL ENERGYLIST$PRINTONE(NFILO,'ELECTRON THERMOSTAT')     
          CALL ENERGYLIST$PRINTONE(NFILO,'EXTERNAL 1CENTER POTENTIAL')     
          CALL ENERGYLIST$PRINTONE(NFILO,'LMTO INTERFACE')     
          CALL ENERGYLIST$PRINTONE(NFILO,'DMFT INTERFACE')     
!
!         == VAN DER WAALS ENERGY===============================================
          CALL VDW$GETL4('ON',TCHK)
          IF(TCHK) THEN
            CALL ENERGYLIST$PRINTONE(NFILO,'VAN DER WAALS ENERGY')
          END IF
!
!         == CELLOSTAT   =======================================================
          CALL CELL$GETL4('ON',TCHK)
          IF(TCHK) THEN
            CALL CELL$GETL4('MOVE',TCHK1)
            IF(TCHK1) THEN
              CALL ENERGYLIST$PRINTONE(NFILO,'CELLOSTAT KINETIC')
            END IF
            CALL ENERGYLIST$PRINTONE(NFILO,'CELLOSTAT POTENTIAL')
          END IF
!          
!         == OCCUPATIONS =======================================================
          CALL DYNOCC$GETL4('DYN',TCHK)
          IF(TCHK) THEN
            CALL ENERGYLIST$PRINTONE(NFILO,'OCCUPATIONAL ENTROPY TERM (-TS)')
            CALL ENERGYLIST$PRINTONE(NFILO,'OCCUPATION KINETIC ENERGY')
          END IF
!
!         == EXTERNAL POTENTIAL ================================================
          CALL ENERGYLIST$PRINTONE(NFILO,'EXTERNAL POTENTIAL')
!
!         == QM-MM =============================================================
          CALL QMMM$GETL4('ON',TQMMM)
          IF(TQMMM) THEN
            CALL ENERGYLIST$PRINTONE(NFILO,'QMMM KINETIC ENERGY')
            CALL ENERGYLIST$PRINTONE(NFILO,'QMMM POTENTIAL ENERGY')
          END IF
!
!         == COSMO =============================================================
          IF(TCOSMO) THEN
            CALL ENERGYLIST$PRINTONE(NFILO,'COSMO KINETIC ENERGY')
            CALL ENERGYLIST$PRINTONE(NFILO,'COSMO POTENTIAL ENERGY')
          END IF
          IF (CALGARY_QMMM) THEN
            WRITE(NFILO,*)'    ----- QM/MM MM CONTRIBUTIONS -----'
            CALL ENERGYLIST$PRINTONE(NFILO,'MM KINETIC ENERGY')
            CALL ENERGYLIST$PRINTONE(NFILO,'MM POT ENERGY')
            CALL ENERGYLIST$PRINTONE(NFILO,'MM TEMPERATURE')
            CALL ENERGYLIST$PRINTONE(NFILO,'MM THERMOSTAT')
            WRITE(NFILO,*)
!           CALL MM_PRINT_ENERGY_TERMS (NFILO,1)
!           CALL MM_PRINT_XYZ (NFILO,1)
            WRITE(NFILO,*)
          END IF

        END IF
      END IF
!   
!     ==========================================================================
!     ==   PROJECTED DENSITY OF STATES                                        ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE STATEANALYSIS')
      IF(TPRINT) THEN
        CALL WAVES$WRITEPDOS
        CALL ATOMLIST$REPORT(NFILO)
        CALL QMMM$REPORT(NFILO)
      ENDIF
                              CALL TRACE$PASS('AFTER STATEANALYSIS')
!   
!     ==========================================================================
!     ==   BANDDATA                                                           ==
!     ==========================================================================
      IF(TPRINT)THEN
        CALL BANDDATA$WRITEFILE
      ENDIF
!   
!     ==========================================================================
!     ==   WRITE FILE FOR COSMOTHERM                                          ==
!     ==========================================================================
      CALL COSMO$PRINTOUT()
!   
!     ==========================================================================
!     ==   CALCULATE OPTICAL MATRIXELEMENTS                                   ==
!     ==========================================================================
!     CALL OPTICS$EVALUATE
!   
!     ==========================================================================
!     ==   WRITE OCCUPATIONS AND ENERGIES TO PROTOCOL                         ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE OCCUPATIONS')
      IF(TPRINT) THEN
        WRITE(NFILO,FMT='()')
        CALL DYNOCC$GETL4('DYN',TCHK) 
        IF(TCHK) THEN
          CALL DYNOCC$REPORT(NFILO)
        ELSE
          CALL WAVES$REPORTEIG(NFILO)
        END IF
        CALL CORE$REPORT(NFILO)
      END IF
!   
!     ==========================================================================
!     ==   WRITE CONSTRAINT INFORMATION TO PROTOCOL                           ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE CONSTRAINTS')
      IF(TPRINT) THEN
        CALL CONSTRAINTS$REPORT(NFILO,'FULL')
      END IF
                              CALL TRACE$PASS('BEFORE FLUSH')
!   
!     ==========================================================================
!     ==   WRITE CONSTRAINT INFORMATION                                       ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('CONSTRAINTS',NFIL)
      IF(THISTASK.EQ.1) THEN
        WRITE(NFIL,FMT='(20("="),2X,"TIMESTEP: ",I10,2X,20("="))')NFI
        CALL CONSTRAINTS$REPORT(NFIL,'SHORT')
        FLUSH(NFIL)
      END IF
!   
!     ==========================================================================
!     ==   WRITE FILE STRC_OUT                                                ==
!     ==========================================================================
      CALL STRCOUT
      CALL QMMM$GETL4('ON',TQMMM)
      IF(TQMMM) THEN
        CALL FORCEFIELD$GETCH('FORCEFIELD',FFTYPE)
        IF(FFTYPE.EQ.'AMBER') THEN
           CALL FORCEFIELD$WRITE_MMSTRC
        ELSE IF(FFTYPE.EQ.'UFF') THEN
           CALL FORCEFIELD$WRITE_UFFSTRC
        END IF
      END IF
!   
!     ==========================================================================
!     ==   WRITE FILE STRC_OUT                                                ==
!     ==========================================================================
      CALL FILEHANDLER$CLOSE('PROT')
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITETRAJECTORY(NFI,DELT)
!     **************************************************************************
!     **  TRANSFERS TRAJECTORY DATA INTO THE TRAJECTORY OBJECT                **
!     **  DATA WILL BE ACCUMULATED AN BE WRITTEN IN BIGU CHUNKS               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFI
      REAL(8)   ,INTENT(IN)  :: DELT
      REAL(8)                :: CELVIN
      REAL(8)                :: GLIB
      REAL(8)                :: PICO
      REAL(8)                :: SECOND
      INTEGER(4)             :: NTASKS,THISTASK
      REAL(8)   ,ALLOCATABLE :: DWORK(:) 
      REAL(8)                :: TIME
      REAL(8)                :: SVAR
      INTEGER(4)             :: NAT,NATM
      REAL(8)                :: ECONS
      REAL(8)                :: ETOT
      REAL(8)                :: EKINP
      REAL(8)                :: EKINC
      REAL(8)                :: EFFEKIN
      REAL(8)                :: ENOSEP
      REAL(8)                :: ENOSEE
      REAL(8)                :: EKINFC
      REAL(8)                :: HEAT
      REAL(8)                :: OCCKIN
      LOGICAL(4)             :: TQMMM,TCALGARYQMMM,TCHK
      REAL(8)                :: QMMMKIN,QMMMPOT,QMMMTHERM
      REAL(8)                :: EEXT
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
                                 CALL TRACE$PUSH('WRITETRAJECTORY')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL CONSTANTS('PICO',PICO)
      CALL CONSTANTS('SECOND',SECOND)
      CALL CONSTANTS('KB',CELVIN)
      TIME=DBLE(NFI)*DELT
      CALL ATOMLIST$NATOM(NAT)
!
!     ==========================================================================
!     ==  COLLECT ENERGIES                                                    ==
!     ==========================================================================
!     ==========================================================================
!     == ADD UP CONSERVED ENERGY                                              ==
!     ==========================================================================
      ECONS=0.D0
!     
!     == BASIC LDA + ATOMIC AND FICTITIOUS ELECTRONIC KINETIC ENERGY ===========
      CALL ENERGYLIST$GET('TOTAL ENERGY',ETOT)     
      CALL ENERGYLIST$GET('IONIC KINETIC ENERGY',EKINP)
      CALL ENERGYLIST$GET('WAVEFUNCTION KINETIC ENERGY',EKINC)     
      CALL ENERGYLIST$GET('BO-WAVEFUNCTION KINETIC ENERGY',EFFEKIN)
      ECONS=ECONS+EKINC-EFFEKIN+EKINP+ETOT
!     
!     == ELECTRON AND ATOM THERMOSTATS =========================================
      CALL ENERGYLIST$GET('ATOM THERMOSTAT',ENOSEP)     
      CALL ENERGYLIST$GET('ELECTRON THERMOSTAT',ENOSEE)     
      ECONS=ECONS+ENOSEP+ENOSEE
!     
      CALL ENERGYLIST$GET('CONSTRAINT KINETIC ENERGY',EKINFC)     
      ECONS=ECONS+EKINFC
!     
!     == OCCUPATIONS ===========================================================
      CALL ENERGYLIST$GET('OCCUPATIONAL ENTROPY TERM (-TS)',HEAT) ! -T*S_MERMIN
      CALL ENERGYLIST$GET('OCCUPATION KINETIC ENERGY',OCCKIN)
      ECONS=ECONS+HEAT+OCCKIN
!     
!     == QM-MM ENVIRONMENT =====================================================
      CALL QMMM$GETL4('ON',TQMMM)
      IF(TQMMM) THEN
        CALL QMMM$GETR8('EKIN',QMMMKIN)
        CALL QMMM$GETR8('EPOT',QMMMPOT)
        ECONS=ECONS+QMMMKIN   !POTENTIAL ENERGY ALREADY INCLUDED IN ETOT
!     
!       == THERMOSTAT FOR THE ENVIRONMENT ======================================
        CALL QMMM$GETR8('ETHERM',QMMMTHERM)
        ECONS=ECONS+QMMMTHERM
      END IF
!     
!     == QMMM CALGARY IMPLEMENTATION  ==========================================
      TCALGARYQMMM = .FALSE.
      IF(TCALGARYQMMM) THEN
        CALL ENERGYLIST$GET('MM KINETIC ENERGY',QMMMKIN)
        CALL ENERGYLIST$GET('MM THERMOSTAT',QMMMTHERM)
        ECONS=ECONS+QMMMKIN+QMMMTHERM
      END IF
!     
!     == EXTERNAL POTENTIAL=====================================================
      CALL ENERGYLIST$GET('EXTERNAL POTENTIAL',EEXT)
      ECONS=ECONS+EEXT
!   
!     ==========================================================================
!     ==   WRITE ENERGY TRAJECTORY                                            ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE E-TRAJECTORY')
      CALL TRAJECTORYIO$SELECT('ENERGY-TRAJECTORY')
      CALL TRAJECTORYIO$GETL4('ON',TCHK)
      IF(TCHK) THEN
        ALLOCATE(DWORK(MAX(4*NAT,8)))
        CALL CONSTANTS('KB',CELVIN)
        GLIB=3*NAT-3
        IF(GLIB.NE.0.D0) THEN 
          SVAR=NINT(EKINP/(.5D0*GLIB*CELVIN))
        ELSE
          SVAR=0.D0
        END IF
        DWORK(1)=SVAR
        DWORK(2)=EKINC
        DWORK(3)=EKINP
        DWORK(4)=ETOT
        DWORK(5)=ECONS
        DWORK(6)=ENOSEE
        DWORK(7)=ENOSEP
        DWORK(8)=HEAT
        CALL TRAJECTORYIO$ADD(NFI,TIME,8,DWORK)
        DEALLOCATE(DWORK)
      END IF
      CALL TRAJECTORYIO$SELECT('NONE')
!   
!     ==========================================================================
!     ==   WRITE POSITION TRAJECTORY                                          ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE R-TRAJECTORY')
      CALL TRAJECTORYIO$SELECT('POSITION-TRAJECTORY')
      CALL TRAJECTORYIO$GETL4('ON',TCHK)
      IF(TCHK) THEN
        ALLOCATE(DWORK(9+8*NAT))
        CALL CELL$GETR8A('T(0)',9,DWORK(1:9))
        CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,DWORK(9+1:9+3*NAT))
        CALL ATOMLIST$GETR8A('Q',0,NAT,DWORK(9+3*NAT+1:9+3*NAT+NAT))
        CALL ATOMLIST$GETR8A('CHARGEANDMOMENTS',0,4*NAT &
     &                                  ,DWORK(9+3*NAT+NAT+1:9+3*NAT+NAT+4*NAT))
        CALL TRAJECTORYIO$ADD(NFI,TIME,9+8*NAT,DWORK)
        DEALLOCATE(DWORK)
      END IF
      CALL TRAJECTORYIO$SELECT('NONE')
!   
!     ==========================================================================
!     ==   WRITE FORCE TRAJECTORY                                             ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE F-TRAJECTORY')
      CALL TRAJECTORYIO$SELECT('FORCE-TRAJECTORY')
      CALL TRAJECTORYIO$GETL4('ON',TCHK)
      IF(TCHK) THEN
        ALLOCATE(DWORK(4*NAT))
        CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,DWORK)
        DWORK(3*NAT+1:4*NAT)=0.D0 ! SHALL CONTAIN IN FUTURE THE POTENTIALS
        CALL TRAJECTORYIO$ADD(NFI,TIME,4*NAT,DWORK)
        DEALLOCATE(DWORK)
      END IF
      CALL TRAJECTORYIO$SELECT('NONE')
!   
!     ==========================================================================
!     ==   WRITE POSITION TRAJECTORY FOR QMMM                                 ==
!     ==========================================================================
                              CALL TRACE$PASS('BEFORE QM-MM R-TRAJECTORY')
      CALL QMMM$GETL4('ON',TQMMM)
      IF(TQMMM) THEN
         CALL TRAJECTORYIO$SELECT('QMMM-POS-TRA')
         CALL TRAJECTORYIO$GETL4('ON',TCHK)
         IF(TCHK) THEN      
           CALL CLASSICAL$SELECT('QMMM')
           CALL CLASSICAL$GETI4('NAT',NATM)
           ALLOCATE(DWORK(9+4*NATM))
           CALL CELL$GETR8A('T(0)',9,DWORK(1:9))
           CALL CLASSICAL$GETR8A('R(0)',3*NATM,DWORK(10:9+3*NATM))
           CALL CLASSICAL$GETR8A('QEL',NATM,DWORK(10+3*NATM:))

           CALL TRAJECTORYIO$ADD(NFI,TIME,9+4*NATM,DWORK)
           DEALLOCATE(DWORK)
         END IF
         CALL TRAJECTORYIO$SELECT('NONE')
      END IF
!
                              CALL TRACE$POP
      RETURN
      END
!#END IF


