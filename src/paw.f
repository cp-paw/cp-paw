!#if defined(IBMLICENSE)
      PROGRAM MAIN
!     ******************************************************************
!     **     CP-PAW                                                   **
!     **     (C) COPYRIGHT   I B M   CORPORATION   1990-1997          **
!     **     LICENSED MATERIALS -  PROPERTY     OF     I B M          **
!     ******************************************************************
!     **                                                              **
!     **  THIS THE CAR-PARRINELLO FIRST PRINCIPLES MOLECULAR DYNAMICS **
!     **  PROGRAM BASED ON THE PROJECTOR-AUGMENTED PLANE WAVE METHOD  **
!     **  AND THE LOCAL DENSITY APPROXIMATION.                        **
!     **                                                              **
!     **  PLEASE READ THE PAW.README FILE CONTAINING DISCLAIMER       **
!     **  AND USER AGGREEMENT!                                        **
!     **                                                              **
!     **  AUTHOR: PETER E. BLOECHL                                    **
!     **                                                              **
!     ******************************************************************
      use clock_module
      IMPLICIT NONE
      CHARACTER(32) :: DATIME
      CHARACTER(256):: VERSIONTEXT
      CHARACTER(256):: VERSIONINFO
      INTEGER(4)    :: NFILO
      INTEGER(4)    :: NTASKS,THISTASK
      COMMON/VERSION/VERSIONTEXT
!     ******************************************************************
      VERSIONINFO = '@(#) PAW-VERSION %R% CREATED %U% %E%'
      VERSIONTEXT = 'PROGRAM VERSION %R% CREATED %U% %E%'
!
!     ==================================================================
!     == INITIALIZE MPE ROUTINE FOR PARALLEL PROCESSING               ==
!     ==================================================================
      CALL MPE$INIT
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     ==  ENTER CAR-PARRINELLO SIMULATION                             ==
!     ==================================================================
      CALL TRACE$PUSH('MAIN')
                              CALL TIMING$START
      CALL PAW
!
!     ==================================================================
!     ==  END OF CARPARRINELLO CALCULATION                            ==
!     ==================================================================
!     ==PRINTING IS ALLOWED ONLY ON THE FIRST TASK =====================
      IF(THISTASK.EQ.1) THEN      
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        CALL FILEHANDLER$REPORT(NFILO,'USED')
      END IF
!     == TIMING MUST BE CALLED BY ALL NODES =========================== 
      CALL TRACE$PASS('BEFORE TIMING')
      CALL TIMING$PRINT(NFILO,'ALL')
      CALL TRACE$PASS('AFTER TIMING')
      CALL USAGE$REPORT(NFILO)
      CALL TRACE$PASS('AFTER USAGE')
!
!     ==PRINTING IS ALLOWED ONLY ON THE FIRST TASK =====================
      IF(THISTASK.EQ.1) THEN
        CALL CLOCK$NOW(DATIME)        
        WRITE(NFILO,FMT='(72("="))')
        WRITE(NFILO,FMT='(72("="),T15,"  PROGRAM FINISHED ",A32,"  ")')DATIME
        WRITE(NFILO,FMT='(72("="))')
      END IF
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
      CALL TRACE$POP
      CALL ERROR$NORMALSTOP 
      STOP
      END
!
!     .....................................................PAW..........
      SUBROUTINE PAW
!     ******************************************************************
!     **                                                              **
!     **  THIS IS THE MAIN PAW SUBROUTINE                             **
!     **    1) READ AND INITIALIZE                                    **
!     **    2) ITERATE                                                **
!     **    3) CLOSE DOWN                                             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4)   :: TNWSTR
      LOGICAL(4)   :: TSTOPR
      LOGICAL(4)   :: TSTOP
      LOGICAL(4)   :: TRANP
      LOGICAL(4)   :: TFIRST
      LOGICAL(4)   :: TLAST
      LOGICAL(4)   :: TSTOPH
      LOGICAL(4)   :: TPRINT
      LOGICAL(4)   :: TOLATE,TMERMN
      CHARACTER(80):: TEXT
      INTEGER(4)   :: NFILO     ! FILE UNIT OF PROTOCOL
      INTEGER(4)   :: NFI0      ! TIME STEP COUNTER AT START
      INTEGER(4)   :: NFI       ! TIME STEP COUNTER
      INTEGER(4)   :: NOMORE    ! 
      INTEGER(4)   :: IPRINT    ! 
      INTEGER(4)   :: NBEG      ! DECIDES IF RESTART FILE IS READ
      REAL(8)      :: DELT
!     ******************************************************************
                              CALL TRACE$PUSH('PAW')
                              CALL TIMING$CLOCKON('INITIALIZATION')
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
!
!     ==================================================================
!     ==  INITIALIZE ATOMIC SETUPS                                    ==
!     ==================================================================
      CALL SETUP$READ
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
!       ==  OBTAIN OCCUPATIONS FOR MERMIN FUNCTIONAL ===================
        CALL DYNOCC$GETL4('DYN',TMERMN)
        IF(TMERMN) CALL DYNOCC$INIOCC
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
                              CALL TIMING$PRINT(NFILO,'ALL')
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
      CALL TRAJECTORYIO$INITIALIZE(IPRINT)
1000  CONTINUE
!
!     ==================================================================
!     ==  ITERATION CONTROL (PROPER STOP ETC. )                       ==
!     ==================================================================
      NFI=NFI+1
      IF(TSTOP) TLAST=.TRUE.
      CALL STOPIT$GETL4('STOP',TSTOP)
      IF(TSTOP) TLAST=.TRUE.
      IF(NFI.GE.NFI0+NOMORE) TLAST=.TRUE.
      TPRINT=(MOD(NFI,IPRINT).EQ.0.OR.TFIRST.OR.TLAST)

      CALL HYPERFINE$SETL4('WAKE',TPRINT)
      CALL GRAPHICS$SETL4('WAKE',TPRINT)
!
!     ==================================================================
!     ==   WRITE RESTART_OUT                                          ==
!     ==================================================================
      IF(TLAST.OR.(TPRINT.AND.(.NOT.TFIRST))) THEN
                              CALL TIMING$CLOCKON('RESTART I/O') 
        CALL WRITERESTART
                              CALL TIMING$CLOCKOFF('RESTART I/O')
!       CALL MM_PAW_WRITE_RESTART (NFI) ! CALGARY QM/MM IMPLEMENTATION
      END IF
!
!     ==================================================================
!     ==   PERFORM ONE TIME STEP                                      ==
!     ==================================================================
                              CALL TIMING$CLOCKON('TIMESTEP')
!     ==use the line with "not tstop" to avoid an additional last time step
!     if(.not.tstop) CALL TIMESTEP(DELT,TPRINT,NFI,TSTOP,TMERMN)
      CALL TIMESTEP(DELT,TPRINT,NFI,TSTOP,TMERMN)
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
        CALL TRAJECTORYIO$FLUSH
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
      SUBROUTINE TIMESTEP(DELT,TPRINT,NFI,TSTOP,TMERMN)
!     ******************************************************************      
!     ******************************************************************      
      USE TIMESTEP_MODULE ,ONLY : DELTAT,ISTEPNUMBER,TNEWTHERMOSTAT
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: DELT   ! TIME STEP
      LOGICAL(4),INTENT(IN)   :: TPRINT ! ON/OFF SWITCH FOR LONG PRINTOUT
      LOGICAL(4),INTENT(IN)   :: TSTOP  ! ON/OFF SWITCH FOR LAST TIME STEP
      LOGICAL(4),INTENT(IN)   :: TMERMN ! ON/OFF SWITCH FOR MERMIN FUNCTIONAL
      INTEGER(4),INTENT(INOUT):: NFI    ! TIME STEP COUNTER
      INTEGER(4)              :: NFILO
      LOGICAL(4)              :: TFOR   ! ON/OFF SWITCH FOR ATOMIC MOTION
      LOGICAL(4)              :: TGRA   ! ON/OFF SWITCH FOR GRADIENT CORRECTION
      LOGICAL(4)              :: TSPIN  ! ON/OFF SWITCH FOR SPIN POLARIZATION
      INTEGER(4)              :: NSPIN
      REAL(8)                 :: ANNEE  ! FRICTION COEFFICIENT FOR WAVE FUNCTIONS
      REAL(8)                 :: ANNER  ! FRICTION COEFFICIENT FOR NUCLEI
      REAL(8)   ,ALLOCATABLE  :: RHOE(:,:)  !(NNRX1,NSPIN)   
      LOGICAL(4)              :: TDIAG=.FALSE. !ON/OFF SWITCH FOR WAVE 
                                        ! FUNCTION DIAGONALIZATION
      INTEGER(4)              :: NR1,NR2,NR3,NNRX1
      LOGICAL(4)              :: TCHK
      REAL(8)                 :: TEMPINST
      REAL(8)                 :: ENOSE  ! GENERIC THERMOSTAT ENERGY 
      REAL(8)                 :: EKIN   ! GENERIC KINETIC ENERGY
      REAL(8)                 :: SVAR
      LOGICAL(4)              :: TCHK1,TCHK2
      real(8)                 :: stress(3,3)
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
      CALL CONTINUUM$SETR8('TIMESTEP',DELT)
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
!     == EXTERNAL POTENTIAL ACTING ON ATOMS ============================
      CALL VEXT$APPLY
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
      END IF
!
!     == DETERMINE IF CORRECTION FOR WAVE FUNCTION DRAG IS ON ==========
      IF(TCHK1) THEN          ! ATOM THERMOSTAT OR BOTH THERMOSTATS
        CALL WAVES$GETR8('FRICTION',ANNEE)
        CALL ATOMS$SETR8('ANNEE',ANNEE)
        IF(.NOT.TNEWTHERMOSTAT) THEN
          CALL ATOMS$SETR8('ANNEE',0.D0)
        END IF
      ELSE ! FRICTION DYNAMICS (NO THERMOSTAT)
        IF(TCHK2) THEN ! WAVE FUNCTION THERMOSTAT ALONE NOT ALLOWED
          CALL ERROR$MSG('WAVE FUNCTION THERMOSTAT ALONE IS NOT ALLOWED')
          CALL ERROR$STOP('TIMESTEP')
        ELSE 
          CALL ATOMS$GETR8('FRICTION',ANNER)
          CALL ATOMS$SETR8('ANNEE',ANNER)
        END IF
      END IF
!
!     ==================================================================
!     ==================================================================
!     ==  PROPAGATE:                                                  ==
!     ==================================================================
!     ==================================================================
! 
!     ==================================================================
!     ==  PROPAGATE NUCLEI                                            ==
!     ==================================================================
      CALL ATOMS$PROPAGATE()
      CALL ATOMS$CONSTRAINTS()
! 
!     ==================================================================
!     ==  PROPAGATE UNIT CELL                                         ==
!     ==================================================================
      CALL CELL$PROPAGATE()
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
!       CALL ENERGYLIST$RETURN('WAVEFUNCTION KINETIC ENERGY',EKIN)
        CALL WAVES$GETR8('EKIN(PSI)',EKIN)
        CALL THERMOSTAT$SETR8('EKIN(SYSTEM)',EKIN)
        CALL THERMOSTAT$PROPAGATE()
      END IF
!
!     ==================================================================
!     ==================================================================
!     ==  APPLY CONSTRAINTS                                           ==
!     ==================================================================
!     ==================================================================
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
      CALL CELL$GETR8('EPOT',SVAR)
      CALL ENERGYLIST$set('CELLOSTAT KINETIC',EKIN)     
      CALL ENERGYLIST$set('CELLOSTAT POTENTIAL',SVAR)     
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',SVAR+EKIN)
!
!     ==================================================================
!     == OCCUPATIONS                                                  ==
!     ==================================================================
      CALL DYNOCC$GETR8('HEAT',SVAR)
      CALL ENERGYLIST$SET('ELECTRONIC HEAT',SVAR)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',SVAR)
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
!     ==  collect energy of wave function thermostat                  ==
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
!     == WRITE INFORMATION                                            ==
!     ==================================================================
!     ==================================================================
!     __WRITE PROTOCOLL_________________________________________________
      CALL PRINFO(TPRINT,NFI,DELT)
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
!     __UNIT CELL_______________________________________________________
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
                              CALL TRACE$POP
      RETURN
      END
!
!     ...................................................STOPIT.........
MODULE STOPIT_MODULE
!**                                                                   ** 
!**  SET SWITCH TO INITIATE SOFT-KILL                                 ** 
!**  SWITCH CAN BE SET BY                                             ** 
!**   A) CREATE A PREDEFINED EXITfILE                                 ** 
!**   B) BY THE CODE BY SETTINg TSTOP EXPLICITELY                     ** 
!**                                                                   ** 
LOGICAL(4)     :: TSTOP=.FALSE.   ! INITIATE SOFT-KILL
LOGICAL(4)     :: DISTRIBUTED=.FALSE. ! TSTOP=TRUE ON ALL NODES
LOGICAL(4)     :: TNOTIFY=.TRUE.  ! NOTIFICATION ABOUT STOP REQUIRED
LOGICAL(4)     :: EXITFILEREMOVED
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
      SUBROUTINE STOPIT$SET(IDENT_,NBYTE_,VAL_)
!     ****************************************************************** 
!     **  STOPIT$SET                                                  ** 
!     ****************************************************************** 
      USE STOPIT_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: IDENT_
      INTEGER(4)  ,INTENT(IN) :: NBYTE_
      REAL(8)     ,INTENT(IN) :: VAL_
!     ****************************************************************** 
      CALL ERROR$MSG('ROUTINE MARKED FOR DELETION')
      CALL ERROR$STOP('STOPIT$SET')
!     CALL LINKEDLIST$SET($LIST,IDENT_,NBYTE_,VAL_)
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
!     ****************************************************************** 
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     ==  CHECK EXITFILE ONLY FROM THE FIRST TASK                     ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
!       ================================================================
!       ==  REMOVE EXIT FILE IN THE FIRST REQUEST                     ==
!       ================================================================
        IF(.NOT.EXITFILEREMOVED) THEN
          CALL FILEHANDLER$FILENAME('EXIT',EXITFILE)
          INQUIRE(FILE=EXITFILE,EXIST=TCHK)
          IF(TCHK) THEN
!           == CHAR(114)//CHAR(109)='RM'  (LOWERCASE)
            CMD=CHAR(114)//CHAR(109)//' '//EXITFILE
!           call error$msg('system call removed for absoft')
!           call error$stop('stopit$update')
            CALL SYSTEM(CMD)
          END IF
          EXITFILEREMOVED=.TRUE.
        END IF
!       ================================================================
!       ==  CHECK IF EXITFILE EXISTS                                  ==
!       ================================================================
        IF(.NOT.TSTOP) THEN
          CALL FILEHANDLER$FILENAME('EXIT',EXITFILE)
          INQUIRE(FILE=EXITFILE,EXIST=TSTOP)
          IF(TSTOP) THEN
            CALL FILEHANDLER$UNIT('PROT',NFILO)
            WRITE(NFILO,*)'EXITFILE EXISTS. PROGRAM WILL TERMINATE'
            CALL lib$FLUSHfile(NFILO)
            TNOTIFY=.FALSE.
          END IF
        END IF
      ELSE
        EXITFILEREMOVED=.TRUE.
      END IF

!     ==================================================================
!     ==  CHECK WHETHER TSTOP=T FOR ANY TASK                          ==
!     ==================================================================
      IF(.NOT.DISTRIBUTED) THEN
        IF(TSTOP) THEN 
          NVAL=1
        ELSE
          NVAL=0
        ENDIF
        CALL MPE$COMBINE('+',NVAL)
        TSTOP=(NVAL.NE.0) 
        DISTRIBUTED=TSTOP
      END IF
!
!     ==================================================================
!     == WRITE MESSAGE IF NOT DONE ALREADY =============================
!     ==================================================================
      IF(TSTOP.AND.TNOTIFY) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)'STOP SIGNAL RECEIVED'
        CALL lib$FLUSHfile(NFILO)
        TNOTIFY=.FALSE.
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PRINFO(TPRINT,NFI,DELT)
!     ******************************************************************
!     **  reports on the process of the simulation and invokes        ** 
!     **  analysis routines                                           ** 
!     **                                                              ** 
!     **  prinfo is called once per timestep                          ** 
!     ******************************************************************
      USE CONTINUUM_CONTROL_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TPRINT
      INTEGER(4),INTENT(IN) :: NFI
      REAL(8)   ,INTENT(IN) :: DELT
      LOGICAL(4),SAVE       :: TFIRST=.TRUE.
      INTEGER(4)            :: NAT
      INTEGER(4)            :: NFILO
      INTEGER(4)            :: NFIL
      REAL(8)               :: GLIB
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
      REAL(8)               :: Ecellpot
      REAL(8)               :: Ecellkin
      REAL(8)               :: ESOLV,EKINQ,QFRIC,QTOT
      LOGICAL(4)            :: TCHK,tchk1
      LOGICAL(4)            :: TCONTINUUM
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
      CHARACTER(126):: NAME
!     ******************************************************************
                              CALL TRACE$PUSH('PRINFO')
      TIME=DBLE(NFI)*DELT
      CALL ATOMLIST$NATOM(NAT)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      GLIB=3*NAT-3
      CALL MPE$QUERY(NTASKS,THISTASK)
!   
!     ==================================================================
!     ==   HYPERFINE PARAMETERS                                       ==
!     ==================================================================
      IF(TPRINT) CALL HYPERFINE$PRINT
!   
!     ==================================================================
!     ==   PLOT WAVE FUNCTIONS, DENSITIES                             ==
!     ==================================================================
      IF(TPRINT) THEN 
        CALL GRAPHICS$PLOT
      END IF
                             CALL TRACE$PASS('AFTER GRAPHICS$PLOT')
!   
!     ==================================================================
!     ==   THE FOLLOWING IS ONLY EXECUTED ON THE FIRST NODE           ==
!     ==================================================================
!     IF(THISTASK.GT.1) THEN
!       CALL TRACE$POP
!       RETURN
!     END IF

!     CALL MM_RUN(CALGARY_QMMM)
      CALGARY_QMMM = .FALSE.
    
!     ==================================================================
!     ==   WRITE HEADER FOR PROTOCOLL FOR EACH TIME STEP              ==
!     ==================================================================
                             CALL TRACE$PASS('WRITE HEADER')
      IF(THISTASK.EQ.1.AND.TFIRST) THEN
        IF(TQMMM) THEN
          WRITE(NFILO,FMT='()')
          WRITE(NFILO,FMT='(2X,A5,A9,1X,A4,1X,A9,2A11,2A6,3A10)') &
     &            'NFI','T[PSEC]','T[K]','EKIN(PSI)','E(RHO)','ECONS', &
     &            'ANNEE','ANNER','T(ENV)','EPOT(ENV)','ECNS(ENV)'
        ELSE IF (CALGARY_QMMM) THEN
          WRITE(NFILO,FMT='(/5X,"NFI",3X,"TIME",1X,"TEMP",3X,"EKINC"   &
     &                  ,5X,"E(RHO)",6X,"ECONS",1X,"ANNEE",1X,"ANNER", &
     &                  "   E_MM  TEMP FRIC")')
        ELSE
          WRITE(NFILO,FMT='()')
          WRITE(NFILO,FMT='(2X,A5,A9,1X,A4,1X,A9,2A11,2A8)') &
     &            'NFI','T[PSEC]','T[K]','EKIN(PSI)','E(RHO)','ECONS', &
     &            'ANNEE','ANNER'
        END IF
        TFIRST=.FALSE.
      END IF
      IF(TPRINT) TFIRST=.TRUE.
!   
!     ==================================================================
!     ==   WRITE PROTOCOLL FOR EACH TIME STEP                         ==
!     ==================================================================
                             CALL TRACE$PASS('WRITE PROTOCOLL')
      IF(THISTASK.EQ.1) THEN

        CALL CONSTANTS('PICO',PICO)
        CALL CONSTANTS('SECOND',SECOND)
        CALL CONSTANTS('KB',CELVIN)
        TME1=TIME/(PICO*SECOND)
!       ================================================================
!       == ADD UP CONSERVED ENERGY                                    ==
!       ================================================================
        ECONS=0.D0
!
!       == BASIC LDA + ATOMIC AND FICTITIOUS ELECTRONIC KINETIC ENERGY =
        CALL ENERGYLIST$RETURN('TOTAL ENERGY',ETOT)     
        CALL ENERGYLIST$RETURN('IONIC KINETIC ENERGY',EKINP)
        CALL ENERGYLIST$RETURN('WAVEFUNCTION KINETIC ENERGY',EKINC)     
        CALL ENERGYLIST$RETURN('BO-WAVEFUNCTION KINETIC ENERGY',EFFEKIN)
        ECONS=ECONS+EKINC-EFFEKIN+EKINP+ETOT
!
!       == ELECTRON AND ATOM THERMOSTATS ===============================
        CALL ENERGYLIST$RETURN('CELLOSTAT KINETIC',ECELLKIN)     
        CALL ENERGYLIST$RETURN('CELLOSTAT POTENTIAL',ECELLPOT)     
        ECONS=ECONS+ECELLKIN+ECELLPOT
print*,'ecellkin/pot ',ecellkin,ecellpot,ECELLKIN+ECELLPOT
!
!       == ELECTRON AND ATOM THERMOSTATS ===============================
        CALL ENERGYLIST$RETURN('ATOM THERMOSTAT',ENOSEP)     
        CALL ENERGYLIST$RETURN('ELECTRON THERMOSTAT',ENOSEE)     
        ECONS=ECONS+ENOSEP+ENOSEE
!
        CALL ENERGYLIST$RETURN('CONSTRAINT KINETIC ENERGY',EKINFC)     
        ECONS=ECONS+EKINFC
!
!       == OCCUPATIONS =================================================
        CALL ENERGYLIST$RETURN('ELECTRONIC HEAT',HEAT) ! -T*S_MERMIN
        CALL ENERGYLIST$RETURN('OCCUPATION KINETIC ENERGY',OCCKIN)
        ECONS=ECONS+HEAT+OCCKIN
!
!       == QM-MM ENVIRONMENT ===========================================
        CALL QMMM$GETL4('ON',TQMMM)
        IF(TQMMM) THEN
          CALL QMMM$GETR8('EKIN',QMMMKIN)
          CALL QMMM$GETR8('EPOT',QMMMPOT)
          ECONS=ECONS+QMMMKIN   !POTENTIAL ENERGY ALREADY INCLUDED IN ETOT
!
!         == THERMOSTAT FOR THE ENVIRONMENT ============================
          CALL QMMM$GETR8('ETHERM',QMMMTHERM)
          ECONS=ECONS+QMMMTHERM
        END IF

!       == QMMM CALGARY IMPLEMENTATION  ===============================
        IF (CALGARY_QMMM) THEN
          CALL ENERGYLIST$RETURN('MM KINETIC ENERGY',MM_KINETIC_ENERGY)
          CALL ENERGYLIST$RETURN('MM POT ENERGY',MM_POT_ENERGY)
          CALL ENERGYLIST$RETURN('MM TEMPERATURE',MM_TEMP)
          CALL ENERGYLIST$RETURN('MM THERMOSTAT',MM_NOSE_ENERGY)
          CALL ENERGYLIST$RETURN('MM ATOM FRICTION',MM_FRIC)
          ECONS=ECONS + MM_KINETIC_ENERGY + MM_NOSE_ENERGY
        END IF

!
!       == CONTINUUM (COSMO) ===========================================
        CALL CONTINUUM$RETURNPROTOCOL(TCONTINUUM,ESOLV,EKINQ,QFRIC,QTOT)
        IF(TCONTINUUM) THEN
          CALL ENERGYLIST$RETURN('SURFACE Q EKIN',EKINQ)    ! REPLACE LATER 
                                            !BY "CONTINUUM KINETIC ENERGY"
          ECONS=ECONS+EKINQ
        END IF
!
!       == EXTERNAL POTENTIAL============================================
        CALL ENERGYLIST$RETURN('EXTERNAL POTENTIAL',EEXT)
        ECONS=ECONS+EEXT
!
        CALL ENERGYLIST$RETURN('CONSTANT ENERGY',SVAR)
PRINT*,'CONSTANT ENERGY ',ECONS,SVAR
!
!       == SOME OTHER STUFF =============================================
        CALL ENERGYLIST$RETURN('IONIC TEMPERATURE',TEMPINST)
        ITEMP=NINT(TEMPINST/CELVIN)
        CALL WAVES$GETR8('FRICTION',ANNEE)
        CALL ATOMS$GETR8('FRICTION',ANNER)
!
        IF (TCONTINUUM) THEN
          WRITE(NFILO,FMT='("!>",I5,F9.5,1X,I5,F9.5,2F11.5,2F8.5 &
     &                      ,F11.5,1X,D7.2,1X,D7.2,2F6.3)') &
     &               NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER &
     &              ,ESOLV,EKINQ,QFRIC,QTOT
        ELSE IF(TQMMM) THEN
          CALL CONSTANTS('KB',CELVIN)
          CALL QMMM$GETI4('NAT:ENV',ISVAR)
          QMMMTEMP=2.D0*QMMMKIN/DBLE(ISVAR)/CELVIN
          WRITE(NFILO,FMT='("!>",I5,F9.5,1X,I5,F9.5,2F11.5,2F6.3 &
     &                ,I10,2F10.5)') &
     &                NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER &
     &               ,NINT(QMMMTEMP),QMMMPOT,QMMMKIN+QMMMPOT
        ELSE IF(CALGARY_QMMM) THEN
          IMM_TEMP=NINT(MM_TEMP)
          WRITE(NFILO,FMT='("!>",I5,F8.4,1X,I4,F8.5,2F11.5,2F6.3,1X,F6.3,1X, &
     &         I4,1X,F5.2 )') NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER,   &
     &                     MM_POT_ENERGY, IMM_TEMP, MM_FRIC
        ELSE
          WRITE(NFILO,FMT='("!>",I5,F9.5,1X,I5,F9.5,2F11.5,2F8.5)') &
     &                NFI,TME1,ITEMP,EKINC-EFFEKIN,ETOT,ECONS,ANNEE,ANNER
        ENDIF
!
!       ================================================================
!       ==   WRITE ENERGIES TO PROTOCOLL                              ==
!       ================================================================
                             CALL TRACE$PASS('WRITE ENERGIES')
        IF(TPRINT) THEN
          CALL ENERGYLIST$PRINTHEADER(NFILO)     
!         ==  BASIC LDA ================================================
          CALL ENERGYLIST$PRINTONE(NFILO,'TOTAL ENERGY')     
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  KINETIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  ELECTROSTATIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  EXCHANGE-CORRELATION')     
          CALL ENERGYLIST$PRINTONE(NFILO,'BACKGROUND')     
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  KINETIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  ELECTROSTATIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  EXCHANGE-CORRELATION')     
!         == ATOM KINETIC ENERGY =======================================
          CALL ENERGYLIST$PRINTONE(NFILO,'IONIC KINETIC ENERGY')
          CALL ENERGYLIST$PRINTONE(NFILO,'IONIC TEMPERATURE')
          CALL ENERGYLIST$PRINTONE(NFILO,'ATOM THERMOSTAT')     
!
          CALL ENERGYLIST$PRINTONE(NFILO,'WAVEFUNCTION KINETIC ENERGY')     
          CALL ENERGYLIST$PRINTONE(NFILO,'BO-WAVEFUNCTION KINETIC ENERGY')
          CALL ENERGYLIST$PRINTONE(NFILO,'ELECTRON THERMOSTAT')     
!
!         == CELLOSTAT   ================================================
          CALL ENERGYLIST$PRINTONE(NFILO,'CELLOSTAT KINETIC')
          CALL ENERGYLIST$PRINTONE(NFILO,'CELLOSTAT POTENTIAL')
          
!         == OCCUPATIONS ================================================
          CALL ENERGYLIST$PRINTONE(NFILO,'ELECTRONIC HEAT')
          CALL ENERGYLIST$PRINTONE(NFILO,'OCCUPATION KINETIC ENERGY')
          CALL ENERGYLIST$PRINTONE(NFILO,'EXTERNAL POTENTIAL')
          IF(TCONTINUUM) THEN
            WRITE(NFILO,*)'-----ELECTROSTATIC SOLVATION CONTRIBUTIONS-----'
            CALL ENERGYLIST$PRINTONE(NFILO,'INTER-TRIANGLE-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'INTRA-TRIANGLE-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'ION-TRIANGLE-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'CAVITY-HARDNESS-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'SURFACETENSION-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'TOTAL-SURFACE-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'SURFACE Q EKIN')
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
!     ==================================================================
!     ==   PROJECTED DENSITY OF STATES                                ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE STATEANALYSIS')
      IF(TPRINT) THEN
        CALL WAVES$WRITEPDOS
        CALL ATOMLIST$REPORT(NFILO)
        CALL QMMM$REPORT(NFILO)
      ENDIF
                              CALL TRACE$PASS('AFTER STATEANALYSIS')
!   
!     ==================================================================
!     ==   CALCULATE OPTICAL MATRIXELEMENTS                           ==
!     ==================================================================
!     CALL OPTICS$EVALUATE
!   
!     ==================================================================
!     ==   WRITE OCCUPATIONS AND ENERGIES TO PROTOCOLL                ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE OCCUPATIONS')
      IF(TPRINT) THEN
        WRITE(NFILO,FMT='()')
        CALL DYNOCC$GETL4('DYN',TCHK) 
        CALL WAVES$GETL4('SAFEORTHO',TCHK1)
        IF(TCHK.or. tchk1) THEN
          CALL DYNOCC$REPORT(NFILO)
        ELSE
          CALL WAVES$REPORTEIG(NFILO)
        END IF
      END IF
!   
!     ==================================================================
!     ==   WRITE CONSTRAINT INFORMATION TO PROTOCOLL                  ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE CONSTRAINTS')
      IF(TPRINT) THEN
        CALL CONSTRAINTS$REPORT(NFILO,'FULL')
      END IF
                              CALL TRACE$PASS('BEFORE FLUSH')
      CALL lib$FLUSHfile(NFILO)
!   
!     ==================================================================
!     ==   WRITE CONSTRAINT INFORMATION                               ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('CONSTRAINTS',NFIL)
      if(thistask.eq.1) then
        WRITE(NFIL,FMT='(20("="),2X,"TIMESTEP: ",I10,2X,20("="))')NFI
        CALL CONSTRAINTS$REPORT(NFIL,'SHORT')
        CALL lib$FLUSHfile(NFIL)
      end if
!   
!     ==================================================================
!     ==   WRITE FILE STRC_OUT                                        ==
!     ==================================================================
      CALL STRCOUT
                              CALL TRACE$POP
      RETURN
      END
!
!     ................................................................... 
      SUBROUTINE WRITETRAJECTORY(NFI,DELT)
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
      INTEGER(4)             :: NAT
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
      LOGICAL(4)             :: TQMMM,TCALGARYQMMM
      REAL(8)                :: QMMMKIN,QMMMPOT,QMMMTHERM
      LOGICAL(4)             :: TCONTINUUM
      REAL(8)                :: ESOLV,EKINQ,QFRIC,QTOT
      REAL(8)                :: EEXT
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
                                 CALL TRACE$PUSH('WRITETRAJECTORY')
      CALL MPE$QUERY(NTASKS,THISTASK)
      CALL CONSTANTS('PICO',PICO)
      CALL CONSTANTS('SECOND',SECOND)
      CALL CONSTANTS('KB',CELVIN)
      TIME=DBLE(NFI)*DELT
!
!     ==================================================================
!     ==  COLLECT ENERGIES                                            ==
!     ==================================================================
!     ================================================================
!     == ADD UP CONSERVED ENERGY                                    ==
!     ================================================================
      ECONS=0.D0
!     
!     == BASIC LDA + ATOMIC AND FICTITIOUS ELECTRONIC KINETIC ENERGY =
      CALL ENERGYLIST$RETURN('TOTAL ENERGY',ETOT)     
      CALL ENERGYLIST$RETURN('IONIC KINETIC ENERGY',EKINP)
      CALL ENERGYLIST$RETURN('WAVEFUNCTION KINETIC ENERGY',EKINC)     
      CALL ENERGYLIST$RETURN('BO-WAVEFUNCTION KINETIC ENERGY',EFFEKIN)
      ECONS=ECONS+EKINC-EFFEKIN+EKINP+ETOT
!     
!     == ELECTRON AND ATOM THERMOSTATS ===============================
      CALL ENERGYLIST$RETURN('ATOM THERMOSTAT',ENOSEP)     
      CALL ENERGYLIST$RETURN('ELECTRON THERMOSTAT',ENOSEE)     
      ECONS=ECONS+ENOSEP+ENOSEE
!     
      CALL ENERGYLIST$RETURN('CONSTRAINT KINETIC ENERGY',EKINFC)     
      ECONS=ECONS+EKINFC
!     
!     == OCCUPATIONS =================================================
      CALL ENERGYLIST$RETURN('ELECTRONIC HEAT',HEAT) ! -T*S_MERMIN
      CALL ENERGYLIST$RETURN('OCCUPATION KINETIC ENERGY',OCCKIN)
      ECONS=ECONS+HEAT+OCCKIN
!     
!     == QM-MM ENVIRONMENT ===========================================
      CALL QMMM$GETL4('ON',TQMMM)
      IF(TQMMM) THEN
        CALL QMMM$GETR8('EKIN',QMMMKIN)
        CALL QMMM$GETR8('EPOT',QMMMPOT)
        ECONS=ECONS+QMMMKIN   !POTENTIAL ENERGY ALREADY INCLUDED IN ETOT
!     
!       == THERMOSTAT FOR THE ENVIRONMENT ============================
        CALL QMMM$GETR8('ETHERM',QMMMTHERM)
        ECONS=ECONS+QMMMTHERM
      END IF
!     
!     == QMMM CALGARY IMPLEMENTATION  ===============================
      TCALGARYQMMM = .FALSE.
      IF(TCALGARYQMMM) THEN
        CALL ENERGYLIST$RETURN('MM KINETIC ENERGY',QMMMKIN)
        CALL ENERGYLIST$RETURN('MM THERMOSTAT',QMMMTHERM)
        ECONS=ECONS+QMMMKIN+QMMMTHERM
      END IF
!     
!     == CONTINUUM ===================================================
      CALL CONTINUUM$RETURNPROTOCOL(TCONTINUUM,ESOLV,EKINQ,QFRIC,QTOT)
      IF(TCONTINUUM) THEN
        CALL ENERGYLIST$RETURN('SURFACE Q EKIN',EKINQ)    ! REPLACE 
        ECONS=ECONS+EKINQ        ! LATER BY "CONTINUUM KINETIC ENERGY"
      END IF
!     
!     == EXTERNAL POTENTIAL============================================
      CALL ENERGYLIST$RETURN('EXTERNAL POTENTIAL',EEXT)
      ECONS=ECONS+EEXT
!   
!     ==================================================================
!     ==   WRITE ENERGY TRAJECTORY                                    ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE E-TRAJECTORY')
      CALL ATOMLIST$NATOM(NAT)
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
      CALL TRAJECTORYIO$ADD('ENERGY-TRAJECTORY',NFI,TIME,8,DWORK)
      DEALLOCATE(DWORK)
!   
!     ==================================================================
!     ==   WRITE POSITION TRAJECTORY                                  ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE R-TRAJECTORY')
      ALLOCATE(DWORK(9+4*NAT))
      CALL CELL$GETR8A('T(0)',9,DWORK(1:9))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,DWORK(10:9+3*NAT))
      CALL ATOMLIST$GETR8A('Q',0,NAT,DWORK(10+3*NAT:))
      CALL TRAJECTORYIO$ADD('POSITION-TRAJECTORY',NFI,TIME,9+4*NAT,DWORK)
      DEALLOCATE(DWORK)
!   
!     ==================================================================
!     ==   WRITE FORCE TRAJECTORY                                     ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE F-TRAJECTORY')
      ALLOCATE(DWORK(4*NAT))
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,DWORK)
      DWORK(3*NAT+1:4*NAT)=0.D0 ! SHALL CONTAIN IN FUTURE THE POTENTIALS
      CALL TRAJECTORYIO$ADD('FORCE-TRAJECTORY',NFI,TIME,4*NAT,DWORK)
      DEALLOCATE(DWORK)
!
                              CALL TRACE$POP
      RETURN
      END
!#end if


