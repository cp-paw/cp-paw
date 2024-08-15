!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM MAIN
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)               :: LL_STP
      CHARACTER(256)              :: SETUPREPORT      
      CHARACTER(256)              :: STRING
      INTEGER(4)                  :: NSELECTION
      CHARACTER(256),ALLOCATABLE  :: OUTFILE(:)
      CHARACTER(16) ,ALLOCATABLE  :: SELECTION(:)
      INTEGER(4)                  :: NARGS
      INTEGER(4)                  :: IARG
      INTEGER(4)                  :: NFIL
      INTEGER(4)   ,PARAMETER     :: NCHOICEX=100
      CHARACTER(60)               :: CHOICE(2,NCHOICEX)
      INTEGER(4)                  :: I,NCHOICE
      INTEGER(4)                  :: ISELECTION
      LOGICAL(4)                  :: TCHK
!     **************************************************************************
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT
!
!     ==========================================================================
!     ==  DEFINE LIST OF SELECTIONS                                           ==
!     ==========================================================================
      I=1
      CHOICE(1,I)='SCATTERING' ; 
                  CHOICE(2,I)=-'PHASE SHIFTS (1ST SET: AE,2ND SET:PAW)'; I=I+1
      CHOICE(1,I)='NB'; CHOICE(2,I)=-'NR. OF WAVE FUNCTIONS'; I=I+1
      CHOICE(1,I)='NC'; CHOICE(2,I)=-'NR. OF CORE WAVE FUNCTIONS'; I=I+1
      CHOICE(1,I)='ATOM.L'; CHOICE(2,I)=-'MAIN ANGULAR MOMENTA OF WAVE FUNCTIONS'; I=I+1
      CHOICE(1,I)='ATOM.E'; CHOICE(2,I)=-'ENERGY EIGENVALUES IN H OF WAVE FUNCTIONS'; I=I+1
      CHOICE(1,I)='ATOM.F'; CHOICE(2,I)=-'OCCUPATIONS OF WAVE FUNCTIONS'; I=I+1
      CHOICE(1,I)='NPRO'; CHOICE(2,I)=-'NR. OF PROJECTOR FUNCTIONS'; I=I+1
      CHOICE(1,I)='LPRO'; CHOICE(2,I)=-'MAIN ANGULAR MOMENTA OF PROJECTOR FUNCTIONS'; I=I+1
      CHOICE(1,I)='ID'; CHOICE(2,I)=-'ID OF THIS SETUP CONSTRUCTION'; I=I+1
      CHOICE(1,I)='Z'; CHOICE(2,I)=-'ATOMIC NUMBER'; I=I+1
      CHOICE(1,I)='ZV'; CHOICE(2,I)=-'NR. OF VALENCE ELECTRONS'; I=I+1
      CHOICE(1,I)='AEPHI'; CHOICE(2,I)=-'ALL-ELECTRON PARTIAL WAVES'; I=I+1
      CHOICE(1,I)='PSPHI'; CHOICE(2,I)=-'AUXILIARY PARTIAL WAVES'; I=I+1
      CHOICE(1,I)='NLPHI'; CHOICE(2,I)=-'NODE-LESS PARTIAL WAVES'; I=I+1
      CHOICE(1,I)='QPHI'; CHOICE(2,I)=-'CORE-LESS PARTIAL WAVES'; I=I+1
      CHOICE(1,I)='PRO'; CHOICE(2,I)=-'PROJECTOR FUNCTIONS'; I=I+1
      CHOICE(1,I)='AEPHIDOT'; CHOICE(2,I)=-'ALL-ELECTRON SCATTERING PARTIAL WAVES'; I=I+1
      CHOICE(1,I)='PSPHIDOT'; CHOICE(2,I)=-'AUXILIARY SCATTERING PARTIAL WAVES'; I=I+1
      CHOICE(1,I)='NLPHIDOT'; CHOICE(2,I)=-'NODE-LESS SCATTERING PARTIAL WAVES'; I=I+1
      CHOICE(1,I)='PAWVALENCEPSI'; CHOICE(2,I)=-'VALENCE WAVE FUNCTIONS CALCULATED FROM PAW'; I=I+1
      CHOICE(1,I)='AEVALENCEPSI'; CHOICE(2,I)=-'ALL-ELECTRON VALENCE WAVE FUNCTIONS'; I=I+1
      CHOICE(1,I)='UPSI'; CHOICE(2,I)=-'NODE-LESS WAVE FUNCTIONS'; I=I+1
      CHOICE(1,I)='UPSISM'; CHOICE(2,I)=-'SMALL COMPONENT OF NODE-LESS WAVE FUNCTIONS'; I=I+1
      CHOICE(1,I)='AEPSI'; CHOICE(2,I)=-'ALL-ELECTRON WAVE FUNCTIONS'; I=I+1
      CHOICE(1,I)='AEPSISM'; CHOICE(2,I)=-'SMALL COMPONENT OF ALL-ELECTRON WAVE FUNCTIONS'; I=I+1
      CHOICE(1,I)='PARMS.PSI.RCL'; CHOICE(2,I)=-'METHOD FOR CONSTRUCTING AUXILIARY PARTIAL WAVES'; I=I+1
      CHOICE(1,I)='PARMS.PSI.TYPE'; CHOICE(2,I)=-'METHOD FOR CONSTRUCTING AUXILIARY PARTIAL WAVES'; I=I+1
      CHOICE(1,I)='PARMS.PSI.LAMBDA'; CHOICE(2,I)=-'DECAY PARAMETER FOR CONSTRUCTING AUXILIARY PARTIAL WAVES'; I=I+1
      CHOICE(1,I)='PARMS.CORE.RC'; CHOICE(2,I)=-'CUTOFF RADIUS FOR PSEUDIZED CORE'; I=I+1
      CHOICE(1,I)='PARMS.CORE.POW'; CHOICE(2,I)=-'LEADING POWER AT THE ORIGIN OF THE PSEUDIZED CORE'; I=I+1
      CHOICE(1,I)='PARMS.CORE.VAL0'; CHOICE(2,I)=-'VALUE AT THE ORIGIN OF THE PSEUDIZED CORE'; I=I+1
      CHOICE(1,I)='PARMS.POT.RC'; CHOICE(2,I)=-'CUTOFF RADIUS FOR PSEUDIZED POTENTIAL'; I=I+1
      CHOICE(1,I)='PARMS.POT.POW'; CHOICE(2,I)=-'LEADING POWER AT THE ORIGIN OF THE PSEUDIZED POTENTIAL'; I=I+1
      CHOICE(1,I)='PARMS.POT.VAL0'; CHOICE(2,I)=-'VALUE AT THE ORIGIN OF THE PSEUDIZED POTENTIAL'; I=I+1
      CHOICE(1,I)='PARMS.RCSM'; CHOICE(2,I)=-'NARROW COMPENSATION GAUSSIAN'; I=I+1
      CHOICE(1,I)='PARMS.RBOX'; CHOICE(2,I)=-'RBOX'; I=I+1
      CHOICE(1,I)='POT'; CHOICE(2,I)=-'POTENTIALS [AE,PS,V(PSRHO)]'; I=I+1
      CHOICE(1,I)='PROG'; CHOICE(2,I)=-'PROJECTOR FUNCTIONS IN G-SPACE'; I=I+1
      CHOICE(1,I)='VADDG'; CHOICE(2,I)=-'VADD IN G-SPACE'; I=I+1
      CHOICE(1,I)='CORE'; CHOICE(2,I)=-'AE AND PS CORE'; I=I+1
      NCHOICE=I-1
      CHOICE(1,I)='AECORE'; CHOICE(2,I)=-'ATOMIC CORE DENSITY'; I=I+1
      CHOICE(1,I)='PSCORE'; CHOICE(2,I)=-'PSEUDIZED CORE DENSITY'; I=I+1
!
!     ==========================================================================
!     ==  DETERMINE THE NUMBER OF ELEMENTS IN THE ARGUMENT LIST               ==
!     ==========================================================================
      NARGS=COMMAND_ARGUMENT_COUNT()
!
!     ==========================================================================
!     ==  SEARCH FOR HELP ARGUMENT                                            ==
!     ==========================================================================
      DO IARG=1,NARGS
        CALL GET_COMMAND_ARGUMENT(IARG,STRING)
        IF(STRING.EQ.'?'.OR.STRING.EQ.-'-H') THEN
          CALL ERRORMESSAGE(NCHOICE,CHOICE)
          STOP
        END IF
      ENDDO
!
!     ==========================================================================
!     ==========================================================================
!     ==  DETERMINE INPUT FILE READ IT INTO THE BUFFER                        ==
!     ==========================================================================
!     ==========================================================================
!     == LAST ARGUMENT IS THE INPUT  FILE FOR THE SETUP REPORT =================
      CALL GET_COMMAND_ARGUMENT(NARGS,SETUPREPORT)
!
!     ==========================================================================
!     ==========================================================================
!     ==  RESOLVE ARGUMENT LIST                                               ==
!     ==========================================================================
!     ==========================================================================
      NSELECTION=0
      DO IARG=1,NARGS-1
        CALL GET_COMMAND_ARGUMENT(IARG,STRING)
        IF(+STRING(1:2).EQ.'-S') NSELECTION=NSELECTION+1
      ENDDO
      IF(NSELECTION.EQ.0) THEN
        WRITE(*,*)'INPUT ERROR: OPTION -S IS MANDATORY'
        WRITE(*,*)
        CALL ERRORMESSAGE(NCHOICE,CHOICE)
        STOP
      END IF
      ALLOCATE(SELECTION(NSELECTION))
      ALLOCATE(OUTFILE(NSELECTION))
      SELECTION(:)=' '
      OUTFILE(:)=' '
      ISELECTION=0
      IARG=1
      DO WHILE(IARG.LT.NARGS)
        CALL GET_COMMAND_ARGUMENT(IARG,STRING)
!       == OPTION MUST START WITH A DASH =======================================
        IF(STRING(1:1).NE.'-') THEN
          WRITE(*,*)'INPUT ERROR: ARGUMENT MUST BEGIN WITH A DASH'
          WRITE(*,*)
          CALL ERRORMESSAGE(NCHOICE,CHOICE)
          STOP
        END IF       
!      
        IF(+STRING.EQ.+'-O') THEN
          IARG=IARG+1
          IF(ISELECTION.EQ.0) THEN
            CALL ERROR$MSG('SELECTION MUST BE SPECIFIED BEFORE THE OUTPUT FILE')
            CALL ERROR$STOP('MAIN')
          END IF
          CALL GET_COMMAND_ARGUMENT(IARG,OUTFILE(ISELECTION))
          IF(OUTFILE(ISELECTION)(1:1).EQ.'-'.OR.IARG.EQ.NARGS) THEN
            WRITE(*,*)'INPUT ERROR: OPTION -O MUST BE FOLLOWED BY OUTPUT FILE'
            WRITE(*,*)
            CALL ERRORMESSAGE(NCHOICE,CHOICE)
            STOP
          END IF
          IARG=IARG+1
!
!       == OPTION SPECIFY SELECTION ============================================
        ELSE IF(+STRING.EQ.+'-S') THEN
          IARG=IARG+1
          ISELECTION=ISELECTION+1
          CALL GET_COMMAND_ARGUMENT(IARG,SELECTION(ISELECTION))
          IF(SELECTION(ISELECTION)(1:1).EQ.'-'.OR.IARG.EQ.NARGS) THEN
            WRITE(*,*)'INPUT ERROR: OPTION -S MUST BE FOLLOWED BY SELECTION'
            WRITE(*,*)
            CALL ERRORMESSAGE(NCHOICE,CHOICE)
            STOP
          END IF
          SELECTION(ISELECTION)=+SELECTION(ISELECTION)
          TCHK=.FALSE.
          DO I=1,NCHOICE
            TCHK=TCHK.OR.(TRIM(SELECTION(ISELECTION)).EQ.+TRIM(CHOICE(1,I)))
            IF(TCHK) EXIT
          ENDDO
          IF(.NOT.TCHK) THEN
            WRITE(*,*)'INPUT ERROR: ILLEGAL VALUE FOR SELECTION'
            WRITE(*,*)'SELECTION=',TRIM(SELECTION(ISELECTION))
            WRITE(*,*)
            CALL ERRORMESSAGE(NCHOICE,CHOICE)
            STOP
          END IF
          IARG=IARG+1
        ELSE 
          WRITE(*,*)'INPUT ERROR: UNKNOWN ARGUMENT'
          WRITE(*,*)'SELECTION=',TRIM(STRING)
          WRITE(*,*)
          CALL ERRORMESSAGE(NCHOICE,CHOICE)
          STOP
        END IF
      ENDDO
!
!     ==========================================================================
!     == READ INPUT FILE                                                      ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('STP_REPORT',.FALSE.,SETUPREPORT)
      CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','FORM','FORMATTED')
!
!     ==========================================================================
!     == READ INPUT FILE                                                      ==
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_STP)
      CALL FILEHANDLER$UNIT('STP_REPORT',NFIL)
      CALL LINKEDLIST$READ(LL_STP,NFIL,'~')
!
!     ==========================================================================
!     == DEFINE OUTPUT FILE AND PREPARE OUTPUT                                ==
!     ==========================================================================
      DO ISELECTION=1,NSELECTION
        IF(OUTFILE(ISELECTION).NE.' ') THEN
          CALL FILEHANDLER$SETFILE('DAT',.FALSE.,OUTFILE(ISELECTION))
          CALL FILEHANDLER$SETSPECIFICATION('DAT','STATUS','UNKNOWN')
          CALL FILEHANDLER$SETSPECIFICATION('DAT','POSITION','REWIND')
          CALL FILEHANDLER$SETSPECIFICATION('DAT','ACTION','WRITE')
          CALL FILEHANDLER$SETSPECIFICATION('DAT','FORM','FORMATTED')
          CALL FILEHANDLER$UNIT('DAT',NFIL)
        ELSE 
          NFIL=6
        END IF
!       CALL LINKEDLIST$REPORT(LL_STP,NFIL)
!
!       ======================================================================
!       == SCAN FOR PARAMETERS                                              ==
!       ======================================================================
        CALL PARMSCONSTANTS(LL_STP,NFIL,SELECTION(ISELECTION),TCHK)
        IF(TCHK) CYCLE
!
!       ========================================================================
!       == SCAN FOR POTENTIALS                                                ==
!       ========================================================================
        CALL POTENTIALS(LL_STP,NFIL,SELECTION(ISELECTION),TCHK)
        IF(TCHK)CYCLE
!
!       ========================================================================
!       == SCAN FOR CORE DENSITIES                                            ==
!       ========================================================================
        CALL CORE(LL_STP,NFIL,SELECTION(ISELECTION),TCHK)
        IF(TCHK)CYCLE
!
!       ========================================================================
!       == SCAN FOR PROJECTOR FUNCTIONS AND VADD IN G-SPACE                   ==
!       ========================================================================
        CALL FOURIER(LL_STP,NFIL,SELECTION(ISELECTION),TCHK)
        IF(TCHK)CYCLE
!
!       ========================================================================
!       == TAKE CARE OF SCATTERING PROPERTIES                                 ==
!       ========================================================================
        IF(SELECTION(ISELECTION).EQ.'SCATTERING') THEN
          CALL SCATTERING(LL_STP,NFIL)
!
!       ========================================================================
!       == CONSTANTS                                                          ==
!       ========================================================================
        ELSE IF(SELECTION(ISELECTION).EQ.'NB') THEN
          CALL ATOMCONSTANTS(LL_STP,NFIL,'NB')
        ELSE IF(SELECTION(ISELECTION).EQ.'NC') THEN
          CALL ATOMCONSTANTS(LL_STP,NFIL,'NC')
        ELSE IF(SELECTION(ISELECTION).EQ.'ATOM.L') THEN
          CALL ATOMCONSTANTS(LL_STP,NFIL,'ATOM.L')
        ELSE IF(SELECTION(ISELECTION).EQ.'ATOM.E') THEN
          CALL ATOMCONSTANTS(LL_STP,NFIL,'ATOM.E')
        ELSE IF(SELECTION(ISELECTION).EQ.'ATOM.F') THEN
          CALL ATOMCONSTANTS(LL_STP,NFIL,'ATOM.F')
!   
        ELSE IF(SELECTION(ISELECTION).EQ.'NPRO') THEN
          CALL AUGMENTATIONCONSTANTS(LL_STP,NFIL,'NPRO')
        ELSE IF(SELECTION(ISELECTION).EQ.'LPRO') THEN
          CALL AUGMENTATIONCONSTANTS(LL_STP,NFIL,'LPRO')
!   
        ELSE IF(SELECTION(ISELECTION).EQ.'ID') THEN
          CALL GENERICCONSTANTS(LL_STP,NFIL,'ID')
        ELSE IF(SELECTION(ISELECTION).EQ.'Z') THEN
          CALL GENERICCONSTANTS(LL_STP,NFIL,'Z')
        ELSE IF(SELECTION(ISELECTION).EQ.'ZV') THEN
          CALL GENERICCONSTANTS(LL_STP,NFIL,'ZV')
!   
!       ========================================================================
!       == CONSTRUCT FILE FOR ATOMIC WAVE FUNCTIONS                           ==
!       ========================================================================
        ELSE IF(SELECTION(ISELECTION).EQ.'AEPSI') THEN
          CALL WAVEFUNCTIONS(LL_STP,NFIL,'AEPSI')
        ELSE IF(SELECTION(ISELECTION).EQ.'UPSI') THEN
          CALL WAVEFUNCTIONS(LL_STP,NFIL,'UPSI')!
        ELSE IF(SELECTION(ISELECTION).EQ.'UPSISM') THEN
          CALL WAVEFUNCTIONS(LL_STP,NFIL,'UPSISM')
!   
!       ========================================================================
!       == CONSTRUCT FILE FOR AUGMENTATION                                    ==
!       ========================================================================
        ELSE IF(SELECTION(ISELECTION).EQ.'AEPHI') THEN
          CALL AUGMENTATION(LL_STP,NFIL,'AEPHI')
        ELSE IF(SELECTION(ISELECTION).EQ.'PSPHI') THEN
          CALL AUGMENTATION(LL_STP,NFIL,'PSPHI')
        ELSE IF(SELECTION(ISELECTION).EQ.'NLPHI') THEN
          CALL AUGMENTATION(LL_STP,NFIL,'NLPHI')
        ELSE IF(SELECTION(ISELECTION).EQ.'QPHI') THEN
          CALL AUGMENTATION(LL_STP,NFIL,'QPHI')
        ELSE IF(SELECTION(ISELECTION).EQ.'PRO') THEN
          CALL AUGMENTATION(LL_STP,NFIL,'PRO')
        ELSE IF(SELECTION(ISELECTION).EQ.'AEPHIDOT') THEN
          CALL AUGMENTATION(LL_STP,NFIL,'AEPHIDOT')
        ELSE IF(SELECTION(ISELECTION).EQ.'PSPHIDOT') THEN
          CALL AUGMENTATION(LL_STP,NFIL,'PSPHIDOT')
        ELSE IF(SELECTION(ISELECTION).EQ.'NLPHIDOT') THEN
          CALL AUGMENTATION(LL_STP,NFIL,'NLPHIDOT')
!   
        ELSE IF(SELECTION(ISELECTION).EQ.'PAWVALENCEPSI') THEN
          CALL VALENCEWAVEFUNCTION(LL_STP,NFIL,'AUGPSI')
        ELSE IF(SELECTION(ISELECTION).EQ.'AEVALENCEPSI') THEN
          CALL VALENCEWAVEFUNCTION(LL_STP,NFIL,'AEPSI')
        ELSE IF(SELECTION(ISELECTION).EQ.'PSVALENCEPSI') THEN
          CALL VALENCEWAVEFUNCTION(LL_STP,NFIL,'PSPSI')
!   
        ELSE IF(SELECTION(ISELECTION).EQ.'AEPSISM') THEN
          CALL WAVEFUNCTIONS(LL_STP,NFIL,TRIM(SELECTION(ISELECTION)))
        ELSE IF(SELECTION(ISELECTION).EQ.'AEPSI') THEN
          CALL WAVEFUNCTIONS(LL_STP,NFIL,TRIM(SELECTION(ISELECTION)))
        ELSE IF(SELECTION(ISELECTION).EQ.'UPSI') THEN
          CALL WAVEFUNCTIONS(LL_STP,NFIL,TRIM(SELECTION(ISELECTION)))
        ELSE IF(SELECTION(ISELECTION).EQ.'UPSISM') THEN
          CALL WAVEFUNCTIONS(LL_STP,NFIL,TRIM(SELECTION(ISELECTION)))
        ELSE 
          CALL ERROR$MSG('SELECTION NOT RECOGNIZED')
          CALL ERROR$CHVAL('SELECTION',SELECTION(ISELECTION))
          CALL ERROR$STOP('MAIN')
        END IF
        IF(NFIL.NE.6) CALL FILEHANDLER$CLOSE('DAT')
      ENDDO
      CALL ERROR$NORMALSTOP()
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ERRORMESSAGE(NCHOICE,CHOICE)
!     **************************************************************************
!     **  WRITES AN ERROR MESSAGE TO STANDARD OUT                             **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NCHOICE
      CHARACTER(*),INTENT(IN) :: CHOICE(2,NCHOICE)
      CHARACTER(128)          :: STRING
      INTEGER(4)              :: I
!     **************************************************************************
      WRITE(*,FMT='(A)')'CALLING SEQUENCE:'
      WRITE(*,FMT='(T10,A)')-"PAW_STPA.X OPTIONS INFILE"
      WRITE(*,FMT='(A)')"OPTIONS CAN BE"
      WRITE(*,FMT='(T10,A)')-"-H"
      WRITE(*,FMT='(T10,A)')-"? "
      WRITE(*,FMT='(T10,A)')-"-S"//" SELECTION"
      WRITE(*,FMT='(T10,A)')-"-O"//" OUTFILE" 
      STRING='("THE OPTIONS '//-'-S'//' AND '//-'-O'
      STRING=TRIM(ADJUSTL(STRING))//' CAN BE SPECIFIED SEVERAL TIMES.")'
      WRITE(*,FMT=TRIM(STRING))
      WRITE(*,FMT=*)
       
      STRING='("INFILE IS THE NAME OF THE INPUT FILE, WHICH HAS THE FORM:")'
      WRITE(*,FMT=STRING)
      WRITE(*,FMT='(T10,A)')+"ROOT"//-"_STPFORZ"//+"NN"//-".MYXML"
      STRING='("WHERE NN IS THE ATOMIC NUMBER AND ROOT IS THE'
      STRING=TRIM(ADJUSTL(STRING))//' ROOT NAME OF THE PAW PROJECT")'
      WRITE(*,FMT=STRING)
      WRITE(*,FMT=*)

      STRING='("OUTFILE IS THE NAME OF THE OUTPUT FILE FOR THE PRECEEDING' 
      STRING=TRIM(ADJUSTL(STRING))//' SELECTION SPECIFIED BY '
      STRING=TRIM(ADJUSTL(STRING))//-' -S.")'
      WRITE(*,FMT=STRING)
      WRITE(*,FMT='("THE DEFAULT OUTFILE IS STOUT (STANDARD OUT, TERMINAL.)")') 
      WRITE(*,FMT=*)
      WRITE(*,FMT='("SELECTION CAN BE ONE OF:")')
!
      DO I=1,NCHOICE
        WRITE(*,FMT='(T2,A,T20,A)')-TRIM(CHOICE(1,I)),TRIM(CHOICE(2,I))
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PARMSCONSTANTS(LL_STP,NFIL,ID,TCHK)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*) ,INTENT(IN)    :: ID
      LOGICAL(4)   ,INTENT(OUT)   :: TCHK
      REAL(8)                     :: SVAR
      CHARACTER(64)               :: CHVAR
      INTEGER(4)                  :: ISVAR
      REAL(8)     ,ALLOCATABLE    :: ARR(:)
!     **************************************************************************
      TCHK=.TRUE.
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'PARAMETERS',1)
      IF(ID.EQ.'PARMS.RCSM') THEN
        CALL LINKEDLIST$GET(LL_STP,'RCSM',1,SVAR)
        WRITE(NFIL,*)SVAR
      ELSE IF(ID.EQ.'PARMS.RBOX') THEN
        CALL LINKEDLIST$GET(LL_STP,'RBOX',1,SVAR)
        WRITE(NFIL,*)SVAR
      ELSE IF(ID.EQ.'PARMS.PSI.TYPE') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'PSI',1)
        CALL LINKEDLIST$GET(LL_STP,'TYPE',1,CHVAR)
        WRITE(NFIL,*)CHVAR
      ELSE IF(ID.EQ.'PARMS.PSI.RCL') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'PSI',1)
        CALL LINKEDLIST$SIZE(LL_STP,'RCL',1,ISVAR)
        ALLOCATE(ARR(ISVAR))
        CALL LINKEDLIST$GET(LL_STP,'RCL',1,ARR)
        WRITE(NFIL,*)ARR
        DEALLOCATE(ARR)
      ELSE IF(ID.EQ.'PARMS.PSI.LAMBDA') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'PSI',1)
        CALL LINKEDLIST$SIZE(LL_STP,'LAMBDA',1,ISVAR)
        ALLOCATE(ARR(ISVAR))
        CALL LINKEDLIST$GET(LL_STP,'LAMBDA',1,ARR)
        WRITE(NFIL,*)ARR
        DEALLOCATE(ARR)
      ELSE IF(ID.EQ.'PARMS.CORE.RC') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'CORE',1)
        CALL LINKEDLIST$GET(LL_STP,'RC',1,SVAR)
        WRITE(NFIL,*)SVAR
      ELSE IF(ID.EQ.'PARMS.CORE.POW') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'CORE',1)
        CALL LINKEDLIST$GET(LL_STP,'POW',1,SVAR)
        WRITE(NFIL,*)SVAR
      ELSE IF(ID.EQ.'PARMS.CORE.VAL0') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'CORE',1)
        CALL LINKEDLIST$GET(LL_STP,'VAL0',1,SVAR)
        WRITE(NFIL,*)SVAR
      ELSE IF(ID.EQ.'PARMS.POT.RC') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'POT',1)
        CALL LINKEDLIST$GET(LL_STP,'RC',1,SVAR)
        WRITE(NFIL,*)SVAR
      ELSE IF(ID.EQ.'PARMS.POT.POW') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'POT',1)
        CALL LINKEDLIST$GET(LL_STP,'POW',1,SVAR)
        WRITE(NFIL,*)SVAR
      ELSE IF(ID.EQ.'PARMS.POT.VAL0') THEN
        CALL LINKEDLIST$SELECT(LL_STP,'POT',1)
        CALL LINKEDLIST$GET(LL_STP,'VAL0',1,SVAR)
        WRITE(NFIL,*)SVAR
      ELSE
        TCHK=.FALSE.
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GENERICCONSTANTS(LL_STP,NFIL,ID)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      REAL(8)                     :: SVAR
      CHARACTER(64)               :: CHVAR
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'GENERIC',1)
      IF(ID.EQ.'ID') THEN
        CALL LINKEDLIST$GET(LL_STP,'ID',1,CHVAR)
        WRITE(NFIL,*)CHVAR
      ELSE IF(ID.EQ.'Z') THEN
        CALL LINKEDLIST$GET(LL_STP,'Z',1,SVAR)
        WRITE(NFIL,*)SVAR
      ELSE IF(ID.EQ.'ZV') THEN
        CALL LINKEDLIST$GET(LL_STP,'ZV',1,SVAR)
        WRITE(NFIL,*)SVAR
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('GENERICCONSTANTS')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATIONCONSTANTS(LL_STP,NFIL,ID)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      INTEGER(4)                  :: ISVAR
      INTEGER(4)  ,ALLOCATABLE    :: IARR(:)
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION',1)
      IF(ID.EQ.'NPRO') THEN
        CALL LINKEDLIST$GET(LL_STP,'LNX',1,ISVAR)
        WRITE(NFIL,*)ISVAR
      ELSE IF(ID.EQ.'LPRO') THEN
        CALL LINKEDLIST$GET(LL_STP,'LNX',1,ISVAR)
        ALLOCATE(IARR(ISVAR))
        CALL LINKEDLIST$GET(LL_STP,'LOX',1,IARR)
        WRITE(NFIL,*)IARR
        DEALLOCATE(IARR)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AUGMENTATIONCONSTANTS')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMCONSTANTS(LL_STP,NFIL,ID)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      INTEGER(4)                  :: ISVAR
      INTEGER(4)  ,ALLOCATABLE    :: IARR(:)
      REAL(8)     ,ALLOCATABLE    :: ARR(:)
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'ATOM',1)
      IF(ID.EQ.'NB') THEN
        CALL LINKEDLIST$GET(LL_STP,'NB',1,ISVAR)
        WRITE(NFIL,*)ISVAR
      ELSE IF(ID.EQ.'NC') THEN
        CALL LINKEDLIST$GET(LL_STP,'NC',1,ISVAR)
        WRITE(NFIL,*)ISVAR
      ELSE IF(ID.EQ.'ATOM.L') THEN
        CALL LINKEDLIST$GET(LL_STP,'NB',1,ISVAR)
        ALLOCATE(IARR(ISVAR))
        CALL LINKEDLIST$GET(LL_STP,'L',1,IARR)
        WRITE(NFIL,*)IARR
        DEALLOCATE(IARR)
      ELSE IF(ID.EQ.'ATOM.E') THEN
        CALL LINKEDLIST$GET(LL_STP,'NB',1,ISVAR)
        ALLOCATE(ARR(ISVAR))
        CALL LINKEDLIST$GET(LL_STP,'E',1,ARR)
        WRITE(NFIL,*)ARR
        DEALLOCATE(ARR)
      ELSE IF(ID.EQ.'ATOM.F') THEN
        CALL LINKEDLIST$GET(LL_STP,'NB',1,ISVAR)
        ALLOCATE(ARR(ISVAR))
        CALL LINKEDLIST$GET(LL_STP,'F',1,ARR)
        WRITE(NFIL,*)ARR
        DEALLOCATE(ARR)
      ELSE IF(ID.EQ.'ATOM.SO') THEN
        CALL LINKEDLIST$GET(LL_STP,'NB',1,ISVAR)
        ALLOCATE(ARR(ISVAR))
        CALL LINKEDLIST$GET(LL_STP,'SO',1,ARR)
        WRITE(NFIL,*)ARR
        DEALLOCATE(ARR)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('ATOMCONSTANTS')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AUGMENTATION(LL_STP,NFIL,ID)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      INTEGER(4)                  :: NR
      REAL(8)      ,ALLOCATABLE   :: R(:)
      REAL(8)      ,ALLOCATABLE   :: PHI(:,:)
      INTEGER(4)                  :: LNX
      INTEGER(4)                  :: LN,IR
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'RGRID',1)
      CALL LINKEDLIST$GET(LL_STP,'NR',1,NR)
      ALLOCATE(R(NR))
      CALL LINKEDLIST$GET(LL_STP,'R',1,R)
!
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION',1)
      CALL LINKEDLIST$GET(LL_STP,'LNX',1,LNX)
      ALLOCATE(PHI(NR,LNX))
      DO LN=1,LNX
        IF(ID.EQ.'AEPHI') THEN
          CALL LINKEDLIST$GET(LL_STP,'AEPHI',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'PSPHI') THEN
          CALL LINKEDLIST$GET(LL_STP,'PSPHI',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'NLPHI') THEN
          CALL LINKEDLIST$GET(LL_STP,'NLPHI',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'QPHI') THEN
          CALL LINKEDLIST$GET(LL_STP,'QPHI',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'PRO') THEN
          CALL LINKEDLIST$GET(LL_STP,'PRO',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'AEPHIDOT') THEN
          CALL LINKEDLIST$GET(LL_STP,'AEPHIDOT',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'PSPHIDOT') THEN
          CALL LINKEDLIST$GET(LL_STP,'PSPHIDOT',LN,PHI(:,LN))
        ELSE IF(ID.EQ.'NLPHIDOT') THEN
          CALL LINKEDLIST$GET(LL_STP,'NLPHIDOT',LN,PHI(:,LN))
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('AUGMENTATION')
        END IF
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO IR=1,NR
        WRITE(NFIL,*)R(IR),PHI(IR,:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE VALENCEWAVEFUNCTION(LL_STP,NFIL,ID)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      INTEGER(4)                  :: NR
      REAL(8)      ,ALLOCATABLE   :: R(:)
      REAL(8)      ,ALLOCATABLE   :: PSI(:,:)
      INTEGER(4)                  :: NV
      INTEGER(4)                  :: IB,IR
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'RGRID',1)
      CALL LINKEDLIST$GET(LL_STP,'NR',1,NR)
      ALLOCATE(R(NR))
      CALL LINKEDLIST$GET(LL_STP,'R',1,R)
!
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION',1)
      CALL LINKEDLIST$GET(LL_STP,'NV',1,NV)
      ALLOCATE(PSI(NR,NV))
      DO IB=1,NV
        IF(ID.EQ.'AEPSI') THEN
          CALL LINKEDLIST$GET(LL_STP,'AEPSI',IB,PSI(:,IB))
        ELSE IF(ID.EQ.'PSPSI') THEN
          CALL LINKEDLIST$GET(LL_STP,'PSPSI',IB,PSI(:,IB))
        ELSE IF(ID.EQ.'AUGPSI') THEN
          CALL LINKEDLIST$GET(LL_STP,'AUGPSI',IB,PSI(:,IB))
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('VALENCEWAVEFUNCTIONS')
        END IF
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO IR=1,NR
        WRITE(NFIL,*)R(IR),PSI(IR,:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVEFUNCTIONS(LL_STP,NFIL,ID)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      INTEGER(4)                  :: NR
      REAL(8)      ,ALLOCATABLE   :: R(:)
      REAL(8)      ,ALLOCATABLE   :: PSI(:,:)
      INTEGER(4)                  :: NB
      INTEGER(4)                  :: IB,IR
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'RGRID',1)
      CALL LINKEDLIST$GET(LL_STP,'NR',1,NR)
      ALLOCATE(R(NR))
      CALL LINKEDLIST$GET(LL_STP,'R',1,R)
!
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'ATOM',1)
      CALL LINKEDLIST$GET(LL_STP,'NB',1,NB)
      ALLOCATE(PSI(NR,NB))
      DO IB=1,NB
        IF(ID.EQ.'AEPSI') THEN
          CALL LINKEDLIST$GET(LL_STP,'AEPSI',IB,PSI(:,IB))
        ELSE IF(ID.EQ.'AEPSISM') THEN
          CALL LINKEDLIST$GET(LL_STP,'AEPSISM',IB,PSI(:,IB))
        ELSE IF(ID.EQ.'UPSI') THEN
          CALL LINKEDLIST$GET(LL_STP,'UPSI',IB,PSI(:,IB))
!          PSI(:,IB)=PSI(:,IB)/MAXVAL(ABS(PSI(:,IB)))
        ELSE IF(ID.EQ.'UPSISM') THEN
          CALL LINKEDLIST$GET(LL_STP,'UPSI_SMALL',IB,PSI(:,IB))
        ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVEFUNCTIONS')
        END IF
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO IR=1,NR
        WRITE(NFIL,*)R(IR),PSI(IR,:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POTENTIALS(LL_STP,NFIL,ID,TCHK)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      LOGICAL(4)   ,INTENT(OUT)   :: TCHK
      INTEGER(4)                  :: NR
      REAL(8)      ,ALLOCATABLE   :: R(:)
      REAL(8)      ,ALLOCATABLE   :: POTS(:,:)
      INTEGER(4)                  :: IR
!     **************************************************************************
      TCHK=ID.EQ.'POT'
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'RGRID',1)
      CALL LINKEDLIST$GET(LL_STP,'NR',1,NR)
      ALLOCATE(R(NR))
      CALL LINKEDLIST$GET(LL_STP,'R',1,R)
!
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'ATOM',1)
      ALLOCATE(POTS(NR,3))
      CALL LINKEDLIST$GET(LL_STP,'AEPOT',1,POTS(:,1))
      CALL LINKEDLIST$GET(LL_STP,'PSPOT',1,POTS(:,2))
      CALL LINKEDLIST$GET(LL_STP,'VADD',1,POTS(:,3))
      POTS(:,3)=POTS(:,2)-POTS(:,3)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO IR=1,NR
        WRITE(NFIL,*)R(IR),POTS(IR,:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORE(LL_STP,NFIL,ID,TCHK)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      LOGICAL(4)   ,INTENT(OUT)   :: TCHK
      INTEGER(4)                  :: NR
      REAL(8)      ,ALLOCATABLE   :: R(:)
      REAL(8)      ,ALLOCATABLE   :: CORES(:,:)
      INTEGER(4)                  :: IR
!     **************************************************************************
      TCHK=ID.EQ.'CORE'
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'RGRID',1)
      CALL LINKEDLIST$GET(LL_STP,'NR',1,NR)
      ALLOCATE(R(NR))
      CALL LINKEDLIST$GET(LL_STP,'R',1,R)
!
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION',1)
      ALLOCATE(CORES(NR,2))
      CALL LINKEDLIST$GET(LL_STP,'AECORE',1,CORES(:,1))
      CALL LINKEDLIST$GET(LL_STP,'PSCORE',1,CORES(:,2))
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO IR=1,NR
        WRITE(NFIL,*)R(IR),CORES(IR,:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FOURIER(LL_STP,NFIL,ID,TCHK)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      CHARACTER(*),INTENT(IN)     :: ID
      LOGICAL(4)   ,INTENT(OUT)   :: TCHK
      INTEGER(4)                  :: NG
      INTEGER(4)                  :: NPRO
      INTEGER(4)                  :: ISVAR
      REAL(8)      ,ALLOCATABLE   :: G(:)
      REAL(8)      ,ALLOCATABLE   :: F(:,:)
      INTEGER(4)                  :: IG
!     **************************************************************************
      TCHK=(ID.EQ.'PROG').OR.(ID.EQ.'VADDG')
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'BESSELTRANSFORMED',1)
      CALL LINKEDLIST$GET(LL_STP,'NG',1,NG)
      ALLOCATE(G(NG))
      CALL LINKEDLIST$GET(LL_STP,'G',1,G)
      IF(ID.EQ.'PROG') THEN
        CALL LINKEDLIST$SIZE(LL_STP,'PRO',1,ISVAR)
        NPRO=NINT(REAL(ISVAR)/REAL(NG))
        IF(NPRO*NG.NE.ISVAR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
          CALL ERROR$I4VAL('NG',NG)
          CALL ERROR$I4VAL('ESTIMATED NPRO',NPRO)
          CALL ERROR$I4VAL('SIZE-NPRO*NG',ISVAR-NPRO*NG)
          CALL ERROR$STOP('FOURIER')
        ENDIF
        ALLOCATE(F(NG,NPRO))
        CALL LINKEDLIST$GET(LL_STP,'PRO',1,F)
!      
      ELSE IF(ID.EQ.'VADDG') THEN
        ALLOCATE(F(NG,1))
        CALL LINKEDLIST$GET(LL_STP,'VADD',1,F)
      END IF
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO IG=1,NG
        WRITE(NFIL,*)G(IG),F(IG,:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCATTERING(LL_STP,NFIL)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      INTEGER(4)   ,INTENT(IN)    :: NFIL
      REAL(8)                     :: EMIN,EMAX
      INTEGER(4)                  :: NE
      INTEGER(4)                  :: LX
      REAL(8)     ,ALLOCATABLE    :: AEPHASE(:,:)
      REAL(8)     ,ALLOCATABLE    :: PSPHASE(:,:)
      INTEGER(4)                  :: IE,L
      REAL(8)                     :: E
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'SCATTERING',1)
      CALL LINKEDLIST$GET(LL_STP,'EMIN',1,EMIN)
      CALL LINKEDLIST$GET(LL_STP,'EMAX',1,EMAX)
      CALL LINKEDLIST$GET(LL_STP,'NE',1,NE)
      CALL LINKEDLIST$GET(LL_STP,'LX',1,LX)
      ALLOCATE(AEPHASE(NE,LX+1))
      ALLOCATE(PSPHASE(NE,LX+1))
      DO L=0,LX
        CALL LINKEDLIST$GET(LL_STP,'AEPHASE',L+1,AEPHASE(:,L+1))
        CALL LINKEDLIST$GET(LL_STP,'PAWPHASE',L+1,PSPHASE(:,L+1))
      ENDDO
      DO IE=1,NE
        E=EMIN+(EMAX-EMIN)/REAL(NE-1,KIND=8)*REAL(IE-1,KIND=8)
        WRITE(NFIL,*)E,AEPHASE(IE,:),PSPHASE(IE,:)
      ENDDO
      DEALLOCATE(AEPHASE)
      DEALLOCATE(PSPHASE)
      RETURN
      END
