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
      IMPLICIT NONE
      CHARACTER(26) :: DATIME
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
!     ==  COLLECT  SCCS INFO ABOUT PROGRAM VERSION                    ==
!     ==  SEE ALSO INPUT ROUTINE READIN                               ==
!     ==================================================================
      CALL TRACE$PUSH('MAIN')
!
!     ==================================================================
!     ==  ENTER CAR-PARRINELLO SIMULATION                             ==
!     ==================================================================
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
        CALL FDATE_(DATIME)
        WRITE(NFILO,FMT='(72("="))')
        WRITE(NFILO,FMT='(72("="),T15,"  PROGRAM FINISHED ",A24,"  ")')DATIME
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
!     **                                                              **
!     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL STOPIT$SIGNALCATCH
      LOGICAL(4)   :: TNWSTR
      LOGICAL(4)   :: TSTOPR
      LOGICAL(4)   :: TSTOP
      LOGICAL(4)   :: TRANP
      LOGICAL(4)   :: TFIRST
      LOGICAL(4)   :: TLAST
      LOGICAL(4)   :: TSTOPH
      LOGICAL(4)   :: TPRINT
      LOGICAL(4)   :: TOLATE,TMERMN,TSIC
      REAL(8)      :: CELVIN
      CHARACTER(80):: TEXT
!     ******************************************************************
                              CALL TRACE$PUSH('PAW')
                              CALL TIMING$CLOCKON('INITIALIZATION')
!
!     ==================================================================
!     ==================================================================
!     ====  READ CONTROL INPUT DATA FILE "CNTL"                     ====
!     ==================================================================
!     ==================================================================
      CALL READIN(NBEG,NOMORE,IPRINT,DELT,TMERMN,TNWSTR)
!
!     ==================================================================
!     ==  READ STRUCTURAL DATA FROM FILE "STRC"                       ==
!     ==================================================================
      CALL STRCIN
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NB)
      CALL DYNOCC$GETI4('NKPT',NKPT) 
      CALL DYNOCC$GETI4('NSPIN',NSPIN)
!
!     ==================================================================
!     ==  ASSIGN UNIT NUMBER TO PROTOCOLL                             ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==================================================================
!     ==  INITIALIZE ATOMS OBJECT                                     ==
!     ==================================================================
      CALL ATOMS$INITIALIZE
!
!     ==================================================================
!     ==  INITIALISE ROUTINES FOR ACCUMULATED DATA                    ==
!     ==================================================================
!     CALL PAIRCR(0,RBAS,NSP,NA,TAU0)
!     CALL ENFLUX(0,EFLUXE,EFLUXP)
!     CALL EDOS(0,NB,EIG)
!
!     ==================================================================
!     ==  INITIALIZE ATOMIC SETUPS                                    ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE READING SETUPS')
      CALL SETUP$READ
!
!     ==================================================================
!     ==  GENERATE G-VECTORS                                          ==
!     ==================================================================
                              CALL TIMING$CLOCKON('POTENTIAL$GVECTORS')
      CALL POTENTIAL$GVECTORS
                              CALL TIMING$CLOCKOFF('POTENTIAL$GVECTORS')
                              CALL TIMING$CLOCKON('WAVES$GVECTORS')
      CALL WAVES$GVECTORS
                              CALL TIMING$CLOCKOFF('WAVES$GVECTORS')
!
!     ==================================================================
!     ==  READ RESTART FILE (WAVE FUNCTION COEFFICIENTS, ETC)         ==
!     ==================================================================
      CALL WAVES$INITIALIZERANDOM
!
                              CALL TIMING$CLOCKON('READ PSI')
      IF(NBEG.GE.0) THEN
                              CALL TRACE$PASS('BEFORE READRESTRT')
        CALL READRESTART
                              CALL TRACE$PASS('AFTER READRESTRT')
!
!       ==  OBTAIN OCCUPATIONS FOR MERMIN FUNCTIONAL ===================
        CALL DYNOCC$GETL4('DYN',TMERMN)
        IF(TMERMN) CALL DYNOCC$INIOCC
      END IF
                              CALL TIMING$CLOCKOFF('READ PSI')
!
!     ================================================================
!     ==  REPORT INPUT DATA                                         ==
!     ================================================================
      CALL IO$REPORT
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
!     ==  CALCULATE PROJECTIONS                                       ==
!     ==================================================================
                              CALL TIMING$CLOCKON('PROJECTIONS')
      CALL WAVES$EIGR('R(-)')
      CALL WAVES$PROJECTIONS('PSI(-)',0)
      CALL WAVES$EIGR('R(0)')
      CALL WAVES$PROJECTIONS('PSI(0)',0)
                              CALL TIMING$CLOCKOFF('PROJECTIONS')
!
!     ==================================================================
!     == IMPOSE CONSTRAINTS                                           ==
!     ==================================================================
                              CALL TIMING$CLOCKON('GRAMM SCHMIDT')
      CALL WAVES$GRAMMSCHMIDT('PSI(-)')
      CALL WAVES$GRAMMSCHMIDT('PSI(0)')
                              CALL TIMING$CLOCKOFF('GRAMM SCHMIDT')
!
!     ==================================================================
!     ==  RUN TIME STATISTICS FOR THE INITIALIZATION                  ==
!     ==================================================================
                              CALL TIMING$CLOCKOFF('INITIALIZATION')
                              CALL TIMING$PRINT(NFILO,'ALL')
                              CALL TIMING$START
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
!     CALL LATER(MAXTIM,TOLATE)
!     IF(TOLATE) TLAST=.TRUE. 
      IF(TSTOP) TLAST=.TRUE.
      CALL SIGNAL(30,STOPIT$SIGNALCATCH)
      CALL STOPIT$GETL4('STOP',TSTOP)
      IF(TSTOP) TLAST=.TRUE.
      IF(NFI.GE.NFI0+NOMORE) TLAST=.TRUE.
      TPRINT=(MOD(NFI,IPRINT).EQ.0.OR.TFIRST.OR.TLAST)

      CALL HYPERFINE$SETL4('WAKE',TPRINT)
!
!     ==================================================================
!     ==   WRITE RESTART_OUT                                          ==
!     ==================================================================
                              CALL TIMING$CLOCKON('TIMESTEP')
      CALL TIMESTEP(DELT,TPRINT,NFI,TSTOP,TMERMN)
                              CALL TIMING$CLOCKOFF('TIMESTEP')
!
!     ==================================================================
!     ==   WRITE RESTART_OUT                                          ==
!     ==================================================================
      IF(TLAST.OR.(TPRINT.AND.(.NOT.TFIRST))) THEN
                              CALL TIMING$CLOCKON('WRITE WRITERESTART')
        CALL WRITERESTART
                              CALL TIMING$CLOCKOFF('WRITE WRITERESTART')
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
        CALL ERROR$MSG('IS NOT RECOGNIZED')
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
        CALL ERROR$MSG('IS NOT RECOGNIZED')
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
        CALL ERROR$MSG('IS NOT RECOGNIZED')
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
        CALL ERROR$MSG('IS NOT RECOGNIZED')
        CALL ERROR$STOP('TIMESTEP$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE TIMESTEP(DELT,TPRINT,NFI,TSTOP,TMERMN)
!     ******************************************************************      
!     ******************************************************************      
      USE TIMESTEP_MODULE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL(4),PARAMETER     :: TNEWTHERMOSTAT=.true.
      REAL(8)   ,INTENT(IN)    :: DELT   ! TIME STEP
      LOGICAL(4),INTENT(IN)    :: TPRINT ! ON/OFF SWITCH FOR LONG PRINTOUT
      LOGICAL(4),INTENT(IN)    :: TSTOP  ! ON/OFF SWITCH FOR LAST TIME STEP
      LOGICAL(4),INTENT(IN)    :: TMERMN ! ON/OFF SWITCH FOR MERMIN FUNCTIONAL
      INTEGER(4),INTENT(INOUT):: NFI    ! TIME STEP COUNTER
      LOGICAL(4)              :: TFOR   ! ON/OFF SWITCH FOR ATOMIC MOTION
      LOGICAL(4)              :: TGRA   ! ON/OFF SWITCH FOR GRADIENT CORRECTION
      LOGICAL(4)              :: TSPIN  ! ON/OFF SWITCH FOR SPIN POLARIZATION
      REAL(8)                 :: ANNEE  ! FRICTION COEFFICIENT FOR WAVE FUNCTIONS
      REAL(8)                 :: ANNER  ! FRICTION COEFFICIENT FOR NUCLEI
      REAL(8)  ,ALLOCATABLE   :: RHOE(:,:)  !(NNRX1,NSPIN)   
      LOGICAL(4)              :: TDIAG=.FALSE. !ON/OFF SWITCH FOR WAVE FUNCTION 
                                        !                   DIAGONALIZATION
      LOGICAL(4)              :: TCHK
      REAL(8)                 :: TEMPINST
!     ******************************************************************      
                              CALL TRACE$PUSH('TIMESTEP')
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!      
      DELTAT=DELT       !-> TIMETSEP_MODULE
      ISTEPNUMBER=NFI   !-> TIMESTEP_MODULE
      CALL POTENTIAL$RESET
      CALL ENERGYLIST$RESET
      CALL WAVES$SETR8('TIMESTEP',DELT)
      CALL ATOMS$SETR8('TIMESTEP',DELT)
      CALL DYNOCC$SETR8('TIMESTEP',DELT)
      CALL QMMM$SETR8('TIMESTEP',DELT)
      CALL CONTINUUM$SETR8('TIMESTEP',DELT)
      CALL THERMOSTAT$SELECT('ATOMS') 
           CALL THERMOSTAT$SETR8('TIMESTEP',DELT)
      CALL THERMOSTAT$SELECT('WAVES') 
           CALL THERMOSTAT$SETR8('TIMESTEP',DELT)
      CALL DIALS$SETR8('TIMESTEP',DELT)
!
      CALL ATOMS$GETL4('MOVE',TFOR)
      CALL DFT$GETL4('GC',TGRA)
      CALL WAVES$SETL4('HAMILTON',TPRINT)
!
!     ==================================================================
!     ==================================================================
!     ==  CALCULATE TOTAL ENERGY AND FORCES                           ==
!     ==================================================================
!     ==================================================================
!
!     ==================================================================
!     ==  CALCULATE EKIN(ELECTR.), <G|PSI> AND RHO(R)                ==
!     ==================================================================
      CALL WAVES$EKIN(EKIN)
!print*,'ekin+++++: ',ekin
        CALL ENERGYLIST$SET('PS  KINETIC',EKIN)
        CALL ENERGYLIST$ADD('AE  KINETIC',EKIN)
        CALL ENERGYLIST$ADD('TOTAL ENERGY',EKIN)
!
                              CALL TIMING$CLOCKON('WAVES$DENSITY')
      CALL DYNOCC$GETI4('NSPIN',NSPIN)
      CALL POTENTIAL$RETURNGRID(NR1,NR2,NR3,NNRX1)
      ALLOCATE(RHOE(NNRX1,NSPIN))
      CALL WAVES$DENSITY(NNRX1,NSPIN,RHOE)
                              CALL TIMING$CLOCKOFF('WAVES$DENSITY')
!
!     ==================================================================
!     ==  CALCULATE <O|PSI> AND <GRAD-O|PSI>                          ==
!     ==================================================================
                              CALL TIMING$CLOCKON('<PS-PSI|PROJECTOR>')
      CALL WAVES$PROJECTIONS('PSI(0)',0)
      IF(TFOR) THEN
        CALL WAVES$PROJECTIONS('PSI(0)',1)
      END IF
                              CALL TIMING$CLOCKOFF('<PS-PSI|PROJECTOR>')
!
!     ==================================================================
!     ==  CALCULATE CHARGEDENSITY AND EFFECTIVE POTENTIAL             ==
!     ==================================================================
      CALL MOMNTS
                              CALL TIMING$CLOCKON('VOFRHO')
      TSPIN=(NSPIN.EQ.2)
      CALL POTENTIAL$VOFRHO(NNRX1,NSPIN,RHOE,TSPIN,TGRA)
                              CALL TIMING$CLOCKOFF('VOFRHO')
!
                              CALL TIMING$CLOCKON('SPHERE')
      CALL SPHERE(TGRA)
      CALL WAVES$1CFORCE
                              CALL TIMING$CLOCKOFF('SPHERE')
!
!     ==================================================================
!     ==================================================================
!     ==  SET FRICTION FOR ANNEALING SCHEDULE                         ==
!     ==================================================================
!     ==================================================================
      CALL AUTOPI(TCHK)
!
      IF(.NOT.TSTOP.AND.TCHK) THEN
        CALL STOPIT$SETL4('STOP',TCHK)
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        WRITE(NFILO,*)'STOP SIGNAL FROM AUTOPILOT'
      END IF
!
!     ==================================================================
!     ==  DETERMINE FRICTION FOR THE ELECTRON DYNAMICS                ==
!     ==  DYNAMICALLY ACCORDING TO THE NOSE DYNAMICS                  ==
!     ==  THE WAVES FRICTION IS OVERWRITTEN                           ==
!     ==  THE ATOMS FRICTION IS CORRECTED                             ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
        CALL THERMOSTAT$GETR8('COOLING',ANNEE)
        CALL WAVES$SETR8('FRICTION',ANNEE)
        IF(TNEWTHERMOSTAT) THEN
          CALL ATOMS$SETR8('ANNEE',ANNEE)
        ELSE
          CALL ATOMS$SETR8('ANNEE',0.D0)
        ENDIF
      END IF
!
!     ==================================================================
!     ==================================================================
!     ==  PROPAGATE IONS :                                            ==
!     ==  ================                                            ==
!     ==  ANNER = 0.0 LEADS TO VERLET ALGHORITHM FOR NEWTON EQATIONS  ==
!     ==     R(+) = 2*R(0) - R(-) + DT**2/MASS * F(0)                 ==
!     ==  ANNER = 1.0 LEADS TO STEEPEST DESCENT                       ==
!     ==     R(+) = R(0) + DT**2/MASS * F(0)                          ==
!     ==================================================================
!     ==================================================================
      EKINP=0.D0
      CALL ATOMS$PROPAGATE
!
!     == CALCULATE KINETIC ENERGY OF THE ATOMS
      CALL ATOMS$EKIN(EKINP)
      CALL ENERGYLIST$SET('IONIC KINETIC ENERGY',EKINP)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',EKINP)
!
      CALL ATOMS$EFFEKIN(EFFEKIN)
      CALL ENERGYLIST$SET('BO-WAVEFUNCTION KINETIC ENERGY',EFFEKIN)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',-EFFEKIN)
!
      CALL ATOMS$TEMPERATURE(TEMPINST)
      CALL ENERGYLIST$SET('IONIC TEMPERATURE',TEMPINST)
!
!     ==================================================================
!     ==================================================================
!     ==  PROPAGATE ELECTRONIC WAVE FUNCTIONS :                       ==
!     ==  =====================================                       ==
!     ==  ANNEE = 0.0 LEADS TO VERLET ALGHORITHM FOR NEWTON EQATIONS  ==
!     ==     C(+) = 2*C(0) - C(-) + DT**2/EMASS * F(0)                ==
!     ==  ANNEE = 1.0 LEADS TO STEEPEST DESCENT                       ==
!     ==     C(+) = C(0) + DT**2/EMASS * F(0)                         ==
!     ==================================================================
!     ==================================================================
!
!     ==================================================================
!     ==  CALCULATE FICTITIOUS KINETIC ENERGY OF THE WAVE FUNCTIONS   ==
!     ==================================================================
      CALL WAVES$WAVEKINETIC(EKINCA)
      CALL ENERGYLIST$ADD('WAVEFUNCTION KINETIC ENERGY',0.5D0*EKINCA)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',0.5D0*EKINCA)
!
!     ==================================================================
!     ==  PROPAGATE WAVE FUNCTIONS WITHOUT CONSTRAINTS                ==
!     ===================================================================
                              CALL TIMING$CLOCKON('EFORCE')
      CALL WAVES$HPSI(NNRX1,NSPIN,RHOE)
      DEALLOCATE(RHOE)
                              CALL TIMING$CLOCKOFF('EFORCE')
      CALL WAVES$PROPAGATE
!
!     ==================================================================
!     ==  CALCULATE FORCE OF CONSTRAINT                               ==
!     ==================================================================
                              CALL TIMING$CLOCKON('OFORCE')
      CALL WAVES$OPSI
                              CALL TIMING$CLOCKOFF('OFORCE')
      IF(TFOR) THEN
                              CALL TIMING$CLOCKON('STRUCTURFACTORS')
        CALL WAVES$EIGR('R(+)')
                              CALL TIMING$CLOCKOFF('STRUCTURFACTORS')
      END IF
                              CALL TIMING$CLOCKON('<PS-PSI|PROJECTOR>')
      CALL WAVES$PROJECTIONS('PSI(+)',0)
      CALL WAVES$PROJECTIONS('O*PSI',0)
                              CALL TIMING$CLOCKOFF('<PS-PSI|PROJECTOR>')
!
!     ==================================================================
!     ==  ORTHOGONALIZE WAVEFUNCTIONS                                 ==
!     ==================================================================
                              CALL TIMING$CLOCKON('ORTHOGONALIZATION')
      CALL WAVES$ORTHOGONALIZE
                              CALL TIMING$CLOCKOFF('ORTHOGONALIZATION')
!
!     ==================================================================
!     ==  CALCULATE FICTITIOUS KINETIC ENERGY OF THE WAVE FUNCTIONS   ==
!     ==================================================================
!      CALL WAVES$WAVEKINETIC(EKINCB)
!      CALL ENERGYLIST$ADD('WAVEFUNCTION KINETIC ENERGY',EKINCB)
!      CALL ENERGYLIST$ADD('CONSTANT ENERGY',EKINCB)

      CALL WAVES$WAVEKINETIC(EKINCB)
      CALL ENERGYLIST$ADD('WAVEFUNCTION KINETIC ENERGY',0.5D0*EKINCB)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',0.5D0*EKINCB)
!
!     ==================================================================
!     ==  PROPAGATE THERMOSTAT FOR THE FICTITIOS ELECTRONIC VARIABLES ==
!     ==================================================================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$GETL4('ON',TCHK)
      IF(TCHK) THEN
        CALL ENERGYLIST$RETURN('WAVEFUNCTION KINETIC ENERGY',EKINC)
        IF(TNEWTHERMOSTAT) THEN
          CALL ATOMS$EFFEKIN(EFFEKIN)
          CALL THERMOSTAT$SETR8('TEMP',EFFEKIN)
        END IF
        CALL THERMOSTAT$SETR8('EKIN(SYSTEM)',EKINC)
        CALL THERMOSTAT$PROPAGATE
        IF(TNEWTHERMOSTAT) THEN
          CALL THERMOSTAT$GETR8('EKIN',ENOSEE)
print*,'thermostat: ',nfi,effekin,ekinc,enosee
        ELSE
          CALL THERMOSTAT$GETR8('ENERGY',ENOSEE)
        END IF
        CALL ENERGYLIST$SET('ELECTRON THERMOSTAT',ENOSEE)
        CALL ENERGYLIST$ADD('CONSTANT ENERGY',ENOSEE)
      END IF
!     EFLUXE=-EKINC*(2.D0*ANNEE/DELT)/EMASS
!
!     ==================================================================
!     ==  PROPAGATE OCCUPATIONS                                       ==
!     ==================================================================
      CALL DYNOCC$PROPAGATE
!
!     == ADD ENERGIES FROM DYNOCC OBJECT TO THE ENERGYLIST
      CALL DYNOCC$GETR8('HEAT',SVAR)
      CALL ENERGYLIST$SET('ELECTRONIC HEAT',SVAR)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',SVAR)
      CALL DYNOCC$GETR8('EKIN',SVAR)
      CALL ENERGYLIST$SET('OCCUPATION KINETIC ENERGY',SVAR)
      CALL ENERGYLIST$ADD('CONSTANT ENERGY',SVAR)
!
!     ==================================================================
!     ==================================================================
!     ==  PRINT ITERATION FOR THIS TIME STEP                          ==
!     ==================================================================
!     ==================================================================
      CALL PRINFO(TPRINT,NFI,DELT)
!
!     ==================================================================
!     ==  ACCUMULATE DATA                                             ==
!     ==================================================================
!     CALL PAIRCR(1,RBAS,NSP,NA,TAU0)
!     IF(NFI.GE.100) THEN
!       CALL ENFLUX(1,EFLUXE,EFLUXP)
!     END IF
!
!     ==================================================================
!     ==================================================================
!     == UPDATE ELECTRONIC AND IONIC VARIABLES                        ==
!     ==================================================================
!     ==================================================================
!     ==  ELECTRONIC WAVE FUNCTIONS  ===================================
!     __INTERCHANGE WAVE FUNCTIONS______________________________________
      CALL WAVES$SWITCH
!     ==  NOSE VARIABLE FOR THE ELECTRONS  ==============================
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$SWITCH
!     ==  ATOMIC POSITIONS =============================================
      CALL ATOMS$SWITCH
      CALL DYNOCC$SWITCH 
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
!     **                                                              ** 
!     **  SET SWITCH TO INITIATE SOFT-KILL                            ** 
!     **  SWITCH CAN BE SET BY                                        ** 
!     **   A) SEND A SIGNAL (KILL -30)                                ** 
!     **   B) CREATE A PREDEFINED EXITILE                            ** 
!     **   C) BY THE CODE BY SETTIN TSTOP EXPLICITELY                 ** 
!     **                                                              ** 
!     **  THE ROUTINE STOPIT$SIGNALCATCH NEEDS TO BE DECLARED         ** 
!     **  EXTERNAL AND THE SIGNAL NEEDS TO BE CAUGHT BY               ** 
!     **        CALL SIGNAL(30,STOPIT$SIGNALCATCH)                    ** 
!     **                                                              ** 
!     **                                                              ** 
!     INTEGER(4)     :: $LIST=0    ! TO BE REMOVED
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
      SUBROUTINE STOPIT$SIGNALCATCH
!     **                                                              ** 
!     **  THE ROUTINE STOPIT IS EXECUTED BY "CALL SIGNAL"             ** 
!     **  IT SETS THE PARAMETER TSTOP, WHICH CAUSES THE PROGRAM TO    **
!     **  TERMINATE AFTER THE NEXT ITERATION                          ** 
!     **                                                              ** 
      USE STOPIT_MODULE
      IMPLICIT NONE
      INTEGER(4) :: NFILO
!     ****************************************************************** 
      TSTOP=.TRUE.
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,*)'EXTERNAL STOP SIGNAL RECEIVED'
      CALL FLUSH_(NFILO)
      TNOTIFY=.FALSE.
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
            CALL FLUSH_(NFILO)
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
!       CALL MPE$OLDCOMBINE('ADD_I*4',4,NVAL)
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
        CALL FLUSH_(NFILO)
        TNOTIFY=.FALSE.
      END IF
      RETURN
      END
!
!     ...................................................PRINFO.........
      SUBROUTINE PRINFO(TPRINT,NFI,DELT)
!     ******************************************************************
!     **                                                              ** 
!     **                                                              ** 
!     ******************************************************************
      use continuum_control_module
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TPRINT
      INTEGER(4),INTENT(IN) :: NFI
      REAL(8)   ,INTENT(IN) :: DELT
      LOGICAL(4),SAVE       :: TFIRST=.TRUE.
      REAL(8)  ,ALLOCATABLE :: DWORK(:) 
      INTEGER(4)            :: NAT
      INTEGER(4)            :: NFILO
      INTEGER(4)            :: NFIL
      REAL(8)               :: GLIB
      INTEGER(4)            :: NTASKs,thisTASK
      REAL(8)               :: PICO
      REAL(8)               :: SECOND
      REAL(8)               :: Time     ! actual time
      REAL(8)               :: Tme1
      REAL(8)               :: ECONS
      REAL(8)               :: ETOT
      REAL(8)               :: EKINP
      REAL(8)               :: EKINC
      REAL(8)               :: ENOSEE
      REAL(8)               :: ENOSEP
      REAL(8)               :: EKINFC
      REAL(8)               :: HEAT
      REAL(8)               :: OCCKIN
      REAL(8)               :: TEMPINST
      REAL(8)               :: CELVIN
      INTEGER(4)            :: ITEMP
      REAL(8)               :: ANNEE
      REAL(8)               :: ANNER
      REAL(8)               :: SVAR
      REAL(8)               :: effekin
      real(8)               :: esolv,ekinq,qfric,qtot
      logical(4)            :: tcontinuum=.true.  ! switch for testing the continuum
      logical(4)            :: tqmmm=.false.
      REAL(8)               :: QMMMKIN   ! EKIN OF QM-MM ENVIRONMENT
      REAL(8)               :: QMMMPOT   ! EPOT OF QM-MM ENVIRONMENT
      REAL(8)               :: QMMMTHERM ! ENERGY OF THE QM-MM THERMOSTAT
!     ******************************************************************
      Time=DBLE(NFI)*DELT
!
                              CALL TRACE$PUSH('PRINFO')
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
      IF(THISTASK.GT.1) THEN
        CALL TRACE$POP
        RETURN
      END IF
!   
!     ==================================================================
!     ==   WRITE HEADER FOR PROTOCOLL FOR EACH TIME STEP              ==
!     ==================================================================
                             CALL TRACE$PASS('WRITE HEADER')
      IF(THISTASK.EQ.1.AND.TFIRST) THEN
        WRITE(NFILO,FMT='()')
          WRITE(NFILO,FMT='(2X,A5,A9,1X,A5,A8,2A11,2A8)') &
     &               'NFI','T[PSEC]','T[K]','EKIN(PSI)','E(RHO)','ECONS' &
     &              ,'ANNEE','ANNER'
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
        ECONS=0.D0
!
!       == BASIC LDA + ATOMIC AND FICTITIOUS ELECTRONIC KINETIC ENERGY =
        CALL ENERGYLIST$RETURN('TOTAL ENERGY',ETOT)     
        CALL ENERGYLIST$RETURN('IONIC KINETIC ENERGY',EKINP)
        CALL ENERGYLIST$RETURN('WAVEFUNCTION KINETIC ENERGY',EKINC)     
        CALL ENERGYLIST$return('BO-WAVEFUNCTION KINETIC ENERGY',EFFEKIN)
        ECONS=ECONS+EKINC-effekin+EKINP+ETOT
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
!
!       == CONTINUUM (COSMO) ===========================================
        IF(TCONTINUUM) THEN
          CALL ENERGYLIST$RETURN('SURFACE Q EKIN',EKINQ)    ! REPLACE LATER 
                                            !BY "CONTINUUM KINETIC ENERGY"
          CALL CONTINUUM$RETURNPROTOCOL(TCONTINUUM,ESOLV,EKINQ,QFRIC,QTOT)
          ECONS=ECONS+EKINQ
        END IF
!
!       == some other stuff =============================================
        CALL ENERGYLIST$RETURN('IONIC TEMPERATURE',TEMPINST)
        ITEMP=NINT(TEMPINST/CELVIN)
        CALL WAVES$GETR8('FRICTION',ANNEE)
        CALL ATOMS$GETR8('FRICTION',ANNER)
!
        IF (TCONTINUUM) THEN
          WRITE(NFILO,FMT='("!>",I5,F9.5,1X,I5,F9.5,2F11.5,2F8.5 &
     &                      ,F11.5,F8.5,F6.3,2F6.3)') &
     &               NFI,TME1,ITEMP,EKINC-effekin,ETOT,ECONS,ANNEE,ANNER &
     &              ,ESOLV,EKINQ,QFRIC,QTOT
        ELSE IF(TQMMM) THEN
          WRITE(NFILO,FMT='("!>",I5,F9.5,1X,I5,F9.5,2F11.5,2F6.3 &
     &                ,3f10.5)') &
     &                NFI,TME1,ITEMP,EKINC-effekin,ETOT,ECONS,ANNEE,ANNER &
     &               ,QMMMKIN,QMMMPOT,QMMMKIN+QMMMPOT
        ELSE
          WRITE(NFILO,FMT='("!>",I5,F9.5,1X,I5,F9.5,2F11.5,2F8.5)') &
     &                NFI,TME1,ITEMP,EKINC-effekin,ETOT,ECONS,ANNEE,ANNER
        ENDIF
!
!       ================================================================
!       ==   WRITE ENERGIES TO PROTOCOLL                              ==
!       ================================================================
                             CALL TRACE$PASS('WRITE ENERGIES')
        IF(TPRINT) THEN
          CALL ENERGYLIST$PRINTHEADER(NFILO)     
!         ==  basic lda ================================================
          CALL ENERGYLIST$PRINTONE(NFILO,'TOTAL ENERGY')     
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  KINETIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  ELECTROSTATIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'AE  EXCHANGE-CORRELATION')     
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  KINETIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  ELECTROSTATIC')     
          CALL ENERGYLIST$PRINTONE(NFILO,'PS  EXCHANGE-CORRELATION')     
!         == atom kinetic energy =======================================
          CALL ENERGYLIST$PRINTONE(NFILO,'IONIC KINETIC ENERGY')
          CALL ENERGYLIST$PRINTONE(NFILO,'IONIC TEMPERATURE')
          CALL ENERGYLIST$PRINTONE(NFILO,'ATOM THERMOSTAT')     
!
          CALL ENERGYLIST$PRINTONE(NFILO,'WAVEFUNCTION KINETIC ENERGY')     
          CALL ENERGYLIST$PRINTONE(NFILO,'BO-WAVEFUNCTION KINETIC ENERGY')
          CALL ENERGYLIST$PRINTONE(NFILO,'ELECTRON THERMOSTAT')     
!         == occupations ================================================
          CALL ENERGYLIST$PRINTONE(NFILO,'ELECTRONIC HEAT')
          CALL ENERGYLIST$PRINTONE(NFILO,'OCCUPATION KINETIC ENERGY')
          CALL ENERGYLIST$PRINTONE(NFILO,'OCCUPATION KINETIC ENERGY')
          if(tcontinuum) then
            WRITE(NFILO,*)'-----ELECTROSTATIC SOLVATION CONTRIBUTIONS-----'
            CALL ENERGYLIST$PRINTONE(NFILO,'INTER-TRIANGLE-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'INTRA-TRIANGLE-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'ION-TRIANGLE-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'CAVITY-HARDNESS-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'SURFACETENSION-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'TOTAL-SURFACE-POTENTIAL')
            CALL ENERGYLIST$PRINTONE(NFILO,'SURFACE Q EKIN')
          end if
        END IF
       END IF
!   
!     ==================================================================
!     ==   PROJECTED DENSITY OF STATES                                ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE STATEANALYSIS')
      IF(TPRINT) THEN
        CALL STATEANALYSIS
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
!       CALL OCCUPATION$REPORT(NFILO,'OCC')
        CALL OCCUPATION$REPORT(NFILO,'EIG')
        WRITE(NFILO,FMT='()')
        CALL DYNOCC$REPORT(NFILO)
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
      IF(NTASKS.EQ.1) THEN
        CALL FLUSH_(NFILO)
      ELSE
        CALL FLUSH_(NFILO)
      END IF
!   
!     ==================================================================
!     ==   WRITE CONSTRAINT INFORMATION                               ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('CONSTRAINTS',NFIL)
      WRITE(NFIL,FMT='(20("="),2X,"TIMESTEP: ",I10,2X,20("="))')NFI
      CALL CONSTRAINTS$REPORT(NFIL,'SHORT')
      CALL FLUSH_(NFIL)
!   
!     ==================================================================
!     ==   WRITE ENERGY TRAJECTORY                                    ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE E-TRAJECTORY')
      ALLOCATE(DWORK(MAX(4*NAT,8)))
      CALL CONSTANTS('KB',CELVIN)
      SVAR=NINT(EKINP/(.5D0*GLIB*CELVIN))
      DWORK(1)=SVAR
      DWORK(2)=EKINC
      DWORK(3)=EKINP
      DWORK(4)=ETOT
      DWORK(5)=ECONS
      DWORK(6)=ENOSEE
      DWORK(7)=ENOSEP
      DWORK(8)=HEAT
      CALL TRAJECTORYIO$ADD('ENERGY-TRAJECTORY',NFI,time,8,DWORK)
!   
!     ==================================================================
!     ==   WRITE POSITION TRAJECTORY                                  ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE R-TRAJECTORY')
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,DWORK(1:3*NAT))
      CALL ATOMLIST$GETR8A('Q',0,NAT,DWORK(1+3*NAT:))
      CALL TRAJECTORYIO$ADD('POSITION-TRAJECTORY',NFI,TIME,4*NAT,DWORK)
!   
!     ==================================================================
!     ==   WRITE FORCE TRAJECTORY                                     ==
!     ==================================================================
                              CALL TRACE$PASS('BEFORE F-TRAJECTORY')
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,DWORK)
      CALL TRAJECTORYIO$ADD('FORCE-TRAJECTORY',NFI,TIME,3*NAT,DWORK)
      DEALLOCATE(DWORK)
!   
!     ==================================================================
!     ==   WRITE FILE STRC_OUT                                        ==
!     ==================================================================
      CALL STRCOUT
!   
!     ==================================================================
!     ==   CLOSE ALL FILES                                            ==
!     ==================================================================
      IF(TPRINT.AND.NTASKS.EQ.1) CALL FILEHANDLER$CLOSEALL
                              CALL TRACE$POP
      RETURN
      END



