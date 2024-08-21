!*******************************************************************************
!**  TODO:                                                                    **
!**     - STANDARDFILES: DO CLEANUP                                           **
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
!     =====   INPUT AND OUTPUT ROUTINES                                     ====
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
!...............................................................................
MODULE IO_MODULE
USE LINKEDLIST_MODULE, ONLY: LL_TYPE
TYPE(LL_TYPE) :: LL_STRC
TYPE(LL_TYPE) :: LL_CNTL
CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE PUTHEADER(NFILO)
        USE CLOCK_MODULE
        IMPLICIT NONE
        INTEGER(4)   ,INTENT(IN) :: NFILO
        CHARACTER(32)            :: DATIME
        INTEGER(4)               :: STATUS
!       ****************************************************************
        WRITE(NFILO,FMT='()')
        WRITE(NFILO,FMT='(80("*"))')
        WRITE(NFILO,FMT='(80("*"),T15' &
     &              //',"                CP-PAW                     ")')
        WRITE(NFILO,FMT='(80("*"),T15' &
     &           //',"     FIRST PRINCIPLES MOLECULAR DYNAMICS      ")')
        WRITE(NFILO,FMT='(80("*"),T15' &
     &           //',"   WITH THE PROJECTOR AUGMENTED WAVE METHOD   ")')
        WRITE(NFILO,FMT='(80("*"))')
        WRITE(NFILO,FMT='(T10' &
     &       //',"P.E. BLOECHL, (C) CLAUSTHAL UNIVERSITY OF TECHNOLOGY (CUT)")')
        WRITE(NFILO,FMT='(T10' &
     &                      //',"DISTRIBUTED UNDER THE GNU PUBLIC LICENSE V3")')
        CALL CPPAW_WRITEVERSION(NFILO)

        CALL CLOCK$NOW(DATIME)
        WRITE(NFILO,FMT='("PROGRAM STARTED on ",A32)')DATIME
        END SUBROUTINE PUTHEADER
!       ................................................................
        SUBROUTINE WRITER8(NFIL,NAME,VALUE,UNIT)
        INTEGER(4)  ,INTENT(IN) :: NFIL
        CHARACTER(*),INTENT(IN) :: NAME
        REAL(8)     ,INTENT(IN) :: VALUE
        CHARACTER(*),INTENT(IN) :: UNIT
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,F10.5,A)')NAME,VALUE,UNIT
        RETURN
        END SUBROUTINE WRITER8
!       ................................................................
        SUBROUTINE WRITEI4(NFIL,NAME,VALUE,UNIT)
        INTEGER(4)  ,INTENT(IN) :: NFIL
        CHARACTER(*),INTENT(IN) :: NAME
        INTEGER(4)  ,INTENT(IN) :: VALUE
        CHARACTER(*),INTENT(IN) :: UNIT
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,I10,A)')NAME,VALUE,UNIT
        RETURN
        END SUBROUTINE WRITEI4
END MODULE IO_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE IO$REPORT
      USE IO_MODULE
      IMPLICIT NONE
      INTEGER(4) :: NFILO
      INTEGER(4)               :: NTASKS,THISTASK
!     ******************************************************************
                           CALL TRACE$PUSH('IO$REPORT')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==================================================================
!     == WAVES                                                        ==
!     ==================================================================
CALL TRACE$PASS('WAVES')
      IF(THISTASK.EQ.1)WRITE(NFILO,FMT='()')
      CALL WAVES$REPORT(NFILO)
!
!     ==================================================================
!     == POTENTIAL                                                    ==
!     ==================================================================
CALL TRACE$PASS('POTENTIAL')
      CALL POTENTIAL$REPORT(NFILO)
!
!     ==================================================================
!     == EXTERNAL ORBITAL POTENTIALS                                   ==
!     ==================================================================
CALL TRACE$PASS('EXT')
      CALL EXTERNAL1CPOT$REPORT(NFILO)
!
!     ==================================================================
!     == OCCUPATIONS                                                  ==
!     ==================================================================
CALL TRACE$PASS('DYNOCC')
      CALL DYNOCC$REPORT(NFILO)
!
!     ==================================================================
!     == ISOLATE                                                      ==
!     ==================================================================
CALL TRACE$PASS('ISOLATE')
      CALL ISOLATE$REPORT(NFILO)
!
!     ==================================================================
!     == K-POINTS                                                     ==
!     ==================================================================
CALL TRACE$PASS('OCCU')
      CALL OCCUPATION$REPORT(NFILO,'KPOINTS')
!
!     ==================================================================
!     == DFT FUNCTIONAL                                               ==
!     ==================================================================
CALL TRACE$PASS('DFT')
      CALL DFT$REPORT(NFILO)
!
!     ==================================================================
!     == ATOM SPECIES                                                 ==
!     ==================================================================
CALL TRACE$PASS('SPECIES')
      CALL SETUP$REPORT(NFILO)
!
!     ==================================================================
!     == LMTO OBJECT                                                  ==
!     ==================================================================
CALL TRACE$PASS('LMTO')
      CALL LMTO$REPORT(NFILO)
CALL TRACE$PASS('SIMPLELMTO')
      CALL SIMPLELMTO$REPORT(NFILO)
!
!     ==================================================================
!     == ATOMS                                                        ==
!     ==================================================================
CALL TRACE$PASS('ATOMS')
      CALL ATOMS$REPORT(NFILO)
CALL TRACE$PASS('ATOMLIST')
      CALL ATOMLIST$REPORT(NFILO)
!
!     ==================================================================
!     == CELL                                                         ==
!     ==================================================================
CALL TRACE$PASS('CELL')
      CALL CELL$REPORT(NFILO)
!
!     ==================================================================
!     == COSMO ENVIRONMENT                                           ==
!     ==================================================================
CALL TRACE$PASS('COSMO')
      CALL COSMO$REPORT(NFILO)
!
!     ==================================================================
!     == CLASSICAL ENVIRONMENT                                        ==
!     ==================================================================
CALL TRACE$PASS('QMMM')
      CALL QMMM$REPORT(NFILO)
!
!     ==================================================================
!     == QM-MM THERMOSTAT                                             ==
!     ==================================================================
CALL TRACE$PASS('QM-MM THERMOSTAT')
      CALL THERMOSTAT$SELECT('QM-MM')
      CALL THERMOSTAT$REPORT(NFILO)
!
!     ==================================================================
!     == ATOM THERMOSTAT                                              ==
!     ==================================================================
CALL TRACE$PASS('ATOM THERMOSTAT')
      CALL THERMOSTAT$SELECT('ATOMS')
      CALL THERMOSTAT$REPORT(NFILO)
!
!     ==================================================================
!     == WAVE FUNCTION THERMOSTAT                                     ==
!     ==================================================================
CALL TRACE$PASS('WAVES THERMOSTAT')
      CALL THERMOSTAT$SELECT('WAVES')
      CALL THERMOSTAT$REPORT(NFILO)
!
!     ==================================================================
!     == MINIMIZER                                                    ==
!     ==================================================================
CALL TRACE$PASS('AUTO')
      CALL AUTO$REPORT(NFILO)
!
!     ==================================================================
!     == GROUPS                                                       ==
!     ==================================================================
CALL TRACE$PASS('GROUPLIST')
      CALL GROUPLIST$REPORT(NFILO)
!
!     ==================================================================
!     == CONSTRAINTS                                                  ==
!     ==================================================================
CALL TRACE$PASS('CONSTRAINT')
      CALL CONSTRAINTS$REPORT(NFILO,'LONG')
!
!     ==================================================================
!     == DIALS                                                        ==
!     ==================================================================
CALL TRACE$PASS('DIALS')
      CALL DIALS$REPORT(NFILO)
!
!     ==================================================================
!     == PARALLELIZATION (MPE)                                        ==
!     ==================================================================
CALL TRACE$PASS('MPE')
      CALL MPE$REPORT(NFILO)
!
!     ==================================================================
!     == FILE HANDLER                                                 ==
!     ==================================================================
CALL TRACE$PASS('FILEHANDLER')
      CALL FILEHANDLER$REPORT(NFILO,'USED')
!
!     ==================================================================
!     == FUSH BUFFER OF PROTOCOL FILE                                 ==
!     ==================================================================
CALL TRACE$PASS('DONE')
      CALL FILEHANDLER$CLOSEALL
                            CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN(NBEG,NOMORE,IPRINT,DELT,TMERMIN,TNWSTR)
!     **************************************************************************
!     **  READ CONTROL INPUT FILE                                             **
!     **                                                                      **
!     **************************************************************************
      USE IO_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(OUT)  :: NBEG
      INTEGER(4)     ,INTENT(OUT)  :: NOMORE
      INTEGER(4)     ,INTENT(OUT)  :: IPRINT
      REAL(8)        ,INTENT(OUT)  :: DELT
      LOGICAL(4)     ,INTENT(OUT)  :: TNWSTR
      LOGICAL(4)     ,INTENT(OUT)  :: TMERMIN
      LOGICAL(4)                   :: TCHK
      INTEGER(4)                   :: NFIL
      INTEGER(4)                   :: NFILO
!      LOGICAL(4)                   :: TOLATE
!      INTEGER(4)                   :: MAXTIM(3)
      CHARACTER(512)               :: CNTLNAME
      CHARACTER(512)               :: ROOTNAME
      CHARACTER(32)                :: CH32SVAR
      LOGICAL(4)                   :: TPR=.FALSE.
      INTEGER(4)                   :: ISVAR
      INTEGER                      :: NTASKS,THISTASK
      INTEGER(4)                   :: LEN
      INTEGER(4)                   :: RUNTIME(3)
      INTEGER(4)   ,ALLOCATABLE    :: SPLITKEY(:)
      INTEGER                      :: ST
!     **************************************************************************
                          CALL TRACE$PUSH('READIN')
!
!     ==========================================================================
!     ==  SET CONTROLFILENAME AND STANDARD ROOT                               ==
!     ==  PAW_VERSION, WHICH RESOLVES THE OPTIONS OF THE EXECUTABLE           ==
!     ==  HAS BEEN CALLED ALREADY FROM THE MAIN PROGRAM. SEE PAW.F90          ==
!     ==========================================================================
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      IF(THISTASK.EQ.1) THEN
        ISVAR=COMMAND_ARGUMENT_COUNT()
        IF (ISVAR.LT.1) THEN
          CALL ERROR$MSG('THE NAME OF THE CONTROLFILE')
          CALL ERROR$MSG('MUST BE GIVEN AS ARGUMENT')
          CALL ERROR$STOP('READIN')
        END IF
        CALL GET_COMMAND_ARGUMENT(1,CNTLNAME,STATUS=ST)
        IF(ST.NE.0) THEN
          CALL ERROR$MSG('FAILURE COLLECTING COMMAND LINE ARGUMENT')
          CALL ERROR$I4VAL('STATUS',ST)
          CALL ERROR$STOP('READIN')
        END IF
      END IF
      CALL MPE$BROADCAST('~',1,CNTLNAME)
!
!     ==========================================================================
!     ==  DEFINE POLYMER IF NECESSARY                                         ==
!     ==========================================================================
      ALLOCATE(SPLITKEY(NTASKS))
      IF(INDEX(CNTLNAME,-'.POLYMER_CNTL',BACK=.TRUE.).NE.0) THEN
!       POLYMER_CNTL DEFINES SEVERAL MONOMERS AND THE ROOTNAMES FOR THE 
!       INDIVIDUAL MONOMERS. IT ALSO PROVIDES INFORMATION ON THE 
!       RELATIVE COMPUTATIONAL EFFORT, WHICH DIVIDE THE PROCESSORS AMONG 
!       INDIVIUDAL MONOMERS.

         !ALEXP-DIMER
         !=== THIS IS A QUICK-FIX FOR THE DIMER
         !=== LATER ONE COULD READ THE POLYMER DATA FROM A SPECIFIC CNTL FILE

         CALL DIMER$SETL4('DIMER',.TRUE.)
         !--- WE SPLIT IN 2 PARTS
         SPLITKEY(1:INT(NTASKS/2))=1
         SPLITKEY(INT(NTASKS/2)+1:NTASKS)=2

         !--- WE SET THE NAME OF THE CNTL-FILE BY HAND
         !    THIS IS ALMOST THE SAME CODE AS BELOW, BUT ENSURES THAT WE DO NOT
         !    HAVE TO CHANGE THE CODE BELOW FOR THE QUICK-FIX
         ISVAR=INDEX(CNTLNAME,-'.POLYMER_CNTL',BACK=.TRUE.)
         IF(ISVAR.GT.0) THEN      
            ROOTNAME=CNTLNAME(1:ISVAR-1)//-'_M'//.ITOS.SPLITKEY(THISTASK)
         ELSE
            CALL ERROR$MSG('ROOTNAME FOR POLYMER EMPTY')
            CALL ERROR$STOP('READIN')
         END IF
         CNTLNAME=TRIM(ADJUSTL(ROOTNAME))//-'.CNTL'
         !ALEXP-DIMER END

      ELSE
        SPLITKEY(:)=1
      END IF
      CALL MPE$NEW('~','MONOMER',NTASKS,SPLITKEY)
      DEALLOCATE(SPLITKEY)
!
!     ==========================================================================
!     ==  SET CONTROLFILENAME AND STANDARD ROOT                               ==
!     ==========================================================================
!     == IF ROOTNAME='-' USE THE ROOT OF THE CONTROLFILE =======================
      ISVAR=INDEX(CNTLNAME,-'.CNTL',BACK=.TRUE.)
      IF(ISVAR.GT.0) THEN      
        ROOTNAME=CNTLNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF

!     == CONNECT CONTROL FILE ==================================================
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL FILEHANDLER$SETFILE(+'CNTL',.FALSE.,CNTLNAME)
      CALL FILEHANDLER$SETSPECIFICATION(+'CNTL','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(+'CNTL','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(+'CNTL','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(+'CNTL','FORM','FORMATTED')
!    
!     ==========================================================================
!     ==  READ BUFFER CNTL                                                    ==
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL FILEHANDLER$UNIT('CNTL',NFIL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'MONOMER')
!    
!     ==========================================================================
!     ==  MARK ALL ELEMENTS AS READ FROM INPUT FILE                           ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$MARK(LL_CNTL,1)
!
!     ==========================================================================
!     ==  READ BLOCK !CONTROL!FILES!FILE                                      ==
!     ==========================================================================
      CALL READIN_FILES(LL_CNTL)
!
!     ==========================================================================
!     ==  WRITE HEADER                                                        ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL PUTHEADER(NFILO)
!
!     ==========================================================================
!     ==  READ BLOCK !DIMER  AND !CONTROL!GENERIC START= (FOR PLACEDIMER)     ==
!     ==========================================================================
      CALL DIMER$GETL4('DIMER',TCHK)
      IF(TCHK) THEN
         CALL READIN_DIMER(LL_CNTL)
      END IF
!    
!     ==========================================================================
!     ==  READ BLOCK !GENERIC                                                 ==
!     ==========================================================================
      CALL READIN_GENERIC(LL_CNTL,NBEG,NOMORE,IPRINT,DELT)
!    
!     ==================================================================
!     ==  READ BLOCK !DFT                                             ==
!     ==================================================================
      CALL READIN_DFT(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !FOURIER                                         ==
!     ==================================================================
      CALL READIN_FOURIER(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !PSIDYN                                         ==
!     ==================================================================
      CALL READIN_PSIDYN(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !RDYN                                            ==
!     ==================================================================
      CALL READIN_RDYN(LL_CNTL)
      CALL AUTO$SETR8('TOLERANCE',1.D-5)
!    
!     ==================================================================
!     ==  READ BLOCK !MERMIN                                          ==
!     ==================================================================
      CALL READIN_MERMIN(LL_CNTL,TMERMIN)
!    
!     ==================================================================
!     ==  READ BLOCK !CELL                                            ==
!     ==================================================================
      CALL READIN_CELL(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !CONTROL!SHADOW                                  ==
!     ==================================================================
      CALL READIN_SHADOW(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !CONTROL!QM-MM                                   ==
!     ==================================================================
      CALL READIN_QMMM(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !CONTROL!COSMO                                   ==
!     ==================================================================
      CALL READIN_COSMO(LL_CNTL)

!     ==================================================================
!     ==  READ BLOCK !CONTROL!QMMM (CALGARY IMPLEMENTATION)
!     ==================================================================
!     CALL READIN_QMMM_CALGARY(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !DATA                                            ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DATA')
!
!     == DEFAULT VALUES ================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NEWSTRUC',0,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'NEWSTRUC',0,.FALSE.)
!
!     == GET NON-DEFAULT VALUES=========================================
      CALL LINKEDLIST$GET(LL_CNTL,'NEWSTRUC',0,TNWSTR)
      IF(TNWSTR) THEN
        CALL ERROR$MSG('THE OPTION NEWSTRUC HAS BEEN REMOVED')
        CALL ERROR$STOP('READIN')
      END IF
!    
!     ==================================================================
!     ==  READ BLOCK !ANALYSE                                         ==
!     ==================================================================
!
!     ==================================================================
!     ==  READ BLOCK !ANALYSE!TRAJECTORIES                            ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
!
!     == DEFAULT VALUES ================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'PDOS',0,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'PDOS',0,.TRUE.)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'BALLSTICK',0,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'BALLSTICK',0,.TRUE.)
!
!     == GET NON-DEFAULT VALUES=========================================
!
      CALL LINKEDLIST$GET(LL_CNTL,'BALLSTICK',1,TCHK)
!     CALL EIGSSETBALLSTICK

      CALL LINKEDLIST$GET(LL_CNTL,'PDOS',1,TCHK)
!     CALL EIGSSETPDOS
!
!     ==================================================================
!     ==  READ BLOCK !ANALYSE!TRA                                     ==
!     ==================================================================
      CALL READIN_ANALYSE_TRAJECTORIES(LL_CNTL)
!
!     ==================================================================
!     ==  READ BLOCK !ANALYSE!WAVE                                    ==
!     ==================================================================
      CALL READIN_ANALYSE_WAVE(LL_CNTL)
!
!     ==================================================================
!     ==  READ BLOCK !ANALYSE!DENSITY                                 ==
!     ==================================================================
      CALL READIN_ANALYSE_DENSITY(LL_CNTL)
!
!     ==================================================================
!     ==  READ BLOCK !ANALYSE!POTENTIAL                               ==
!     ==================================================================
      CALL READIN_ANALYSE_POTENTIAL(LL_CNTL)
!
!     ==================================================================
!     ==  READ BLOCK !ANALYSE!1DPOT                                   ==
!     ==================================================================
      CALL READIN_ANALYSE_1DPOT(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !ANALYSE!POINTCHARGEPOT                          ==
!     ==================================================================
!     CALL READIN_ANALYSE_POTCHARGEPOT(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !ANALYSE!HYPERFINE                               ==
!     ==================================================================
      CALL READIN_ANALYSE_HYPERFINE(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !ANALYSE!HYPERFINE                               ==
!     ==================================================================
      CALL READIN_ANALYSE_CORELEVEL(LL_CNTL)
!    
!     ==================================================================
!     ==  READ BLOCK !ANALYSE!OPTIC                                   ==
!     ==================================================================
      CALL READIN_ANALYSE_OPTIC(LL_CNTL)
      CALL READIN_ANALYSE_OPTEELS(LL_CNTL)
      
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$REPORT_UNUSED(LL_CNTL,NFILO)
!
!     ==================================================================
!     ==  CHECK TIME TO STOP                                          ==
!     ==================================================================
!     IF(MAXTIM(1).EQ.1) WRITE(*,6900)'SUNDAY   ',MAXTIM(2),MAXTIM(3)
!     IF(MAXTIM(1).EQ.2) WRITE(*,6900)'MONDAY   ',MAXTIM(2),MAXTIM(3)
!     IF(MAXTIM(1).EQ.3) WRITE(*,6900)'TUESDAY  ',MAXTIM(2),MAXTIM(3)
!     IF(MAXTIM(1).EQ.4) WRITE(*,6900)'WEDNESDAY',MAXTIM(2),MAXTIM(3)
!     IF(MAXTIM(1).EQ.5) WRITE(*,6900)'THURSDAY ',MAXTIM(2),MAXTIM(3)
!     IF(MAXTIM(1).EQ.6) WRITE(*,6900)'FRIDAY   ',MAXTIM(2),MAXTIM(3)
!     IF(MAXTIM(1).EQ.7) WRITE(*,6900)'SATURDAY ',MAXTIM(2),MAXTIM(3)
!6900  FORMAT(1H ,'STOP ON ',A10,' AT ',I2,':',I2)
!     CALL LATER(MAXTIM,TOLATE)
!     IF(TOLATE) THEN
!       CALL ERROR$MSG('TIME OVER')
!       CALL ERROR$STOP('READIN')
!     END IF
!
!     ==================================================================
!     ==================================================================
!     ==  INTERNAL CHECKS AND PRINOUT                                 ==
!     ==================================================================
!     ==================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==================================================================
!     ==================================================================
!     ==  ELECTRONIC VARIABLES                                        ==
!     ==================================================================
!     ==================================================================
      WRITE(NFILO,FMT='()')
      WRITE(NFILO,FMT='("INFORMATION FROM CONTROL INPUT FILE"' &
     &             //'/"===================================")')
!
      IF(NBEG.LE.-1) THEN
        WRITE(NFILO,FMT='("START WITH RANDOM WAVE FUNCTIONS")')
      ELSE IF(NBEG.GE.0) THEN
        WRITE(NFILO,FMT='("WAVE FUNCTIONS ARE READ FROM FILE")')
      END IF
      IF(TNWSTR) THEN
        WRITE(NFILO,FMT='("ATOMIC STRUCTURE TAKEN FROM "' &
     &               //',"STRUCTURE INPUT FILE"' &
     &               //'/T10,"(STRUCTURE FROM RESTART FILE IGNORED)")')
      END IF
      WRITE(NFILO,FMT='(55("."),": "' &
     &             //',T1,"NUMBER OF ITERATIONS",T58,I10)')NOMORE
      WRITE(NFILO,FMT='(55("."),": "' &
     &             //',T1,"TIME STEP",T58,F10.5," A.U.")')DELT
!
!     ==================================================================
!     ==  DATA FILES                                                  ==
!     ==================================================================
      IF(NBEG.GE.0) THEN
        WRITE(NFILO,FMT='(55("."),": "' &
     &             //',T1,"START WITH WAVE FUNCTIONS FROM FILE"' &
     &             //',T58,"RESTART_IN")')
      END IF
      WRITE(NFILO,FMT='(55("."),": "' &
     &         //',T1,"DETAILED INFORMATION AFTER EACH"' &
     &         //',T58,I10," TIME STEPS")')IPRINT
!
                          CALL TRACE$POP
      RETURN
      END SUBROUTINE READIN
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_FILES(LL_CNTL_)
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      INTEGER(4)               :: NFILE
      INTEGER(4)               :: ITH
      CHARACTER(32)            :: ID
      CHARACTER(256)           :: NAME
      LOGICAL(4)               :: TCHK
      LOGICAL(4)               :: TMOD
      INTEGER(4)               :: THISTASK,NTASKS
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_FILES')  
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'FILES')
!
!     ==================================================================
!     == SET ROOTNAME AND PRCONNECT FILES                             ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ROOT',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'ROOT',1,NAME)
        CALL FILEHANDLER$SETROOT(TRIM(NAME))
      END IF
      CALL STANDARDFILES
!
!     == SCAN SUBLISTS =================================================
      CALL LINKEDLIST$NLISTS(LL_CNTL,'FILE',NFILE)
      DO ITH=1,NFILE
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',ITH)
!     
!       ==  READ ACTUAL VALUES  ========================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!CONTROL!FILES!FILE:ID NOT FOUND')
          CALL ERROR$MSG('NOTE: SYNTAX HAS CHANGED FROM IDENT TO ID')
          CALL ERROR$STOP('READIN_FILES')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,ID)
!
!       == COLLECT FILE NAME ===========================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,NAME)
        ELSE
          CALL ERROR$MSG('!CONTROL!FILES!FILE:NAME NOT FOUND')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('READIN_FILES')
        END IF
!
!       == DEFAULT FOR EXT=.FALSE.
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
        END IF
!     
!       == PERFORM ACTIONS =====================================================
        TMOD=.TRUE.
        IF(THISTASK.NE.1.AND.ID.EQ.'PROT') TMOD=.FALSE.
        IF(TMOD) THEN
          CALL FILEHANDLER$SETFILE(+ID,TCHK,NAME)
        END IF
!
!       == SPECIAL RULES =======================================================
        IF(ID.EQ.'AUGPARMS') THEN
          CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
          CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
        END IF

        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STANDARDFILES
!     **************************************************************************
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(32)        :: ID
      INTEGER(4)           :: NTASKS
      INTEGER(4)           :: THISTASK
!!$      CHARACTER(512)       :: PAWDIR=' '
      INTEGER(4)           :: ST
!     **************************************************************************
                                   CALL TRACE$PUSH('STANDARDFILES')
!!$      CALL GET_ENVIRONMENT_VARIABLE('PAWDIR',PAWDIR,STATUS=ST)
!!$      IF(ST.NE.0.AND.ST.NE.1) THEN
!!$        CALL ERROR$MSG('FAILURE COLLECTING ENVIRONMENT VARIABLE "PAWDIR"')
!!$        CALL ERROR$I4VAL('STATUS',ST)
!!$        IF(ST.EQ.-1) THEN
!!$          CALL ERROR$MSG('ENVIRONMENT VARIABLE DOES NOT FIT INTO STRING')
!!$        END IF
!!$        CALL ERROR$STOP('STANDARDFILES')
!!$      END IF
!  
!     ==========================================================================
!     == SET STANDARD FILENAMES                                               ==
!     ==========================================================================
!
!     ==  EXIT FILE ============================================================
      ID=+'EXIT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.EXIT')
!
!     ==  ERROR FILE ===========================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.ERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PROTOCOL FILE ========================================================
!     ==  EACH MONOMER HAS ITS OWN PROTOCOL FILE
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      ID=+'PROT'
      IF(THISTASK.GT.1)THEN
        CALL FILEHANDLER$SETFILE(ID,.FALSE.,-'/DEV/NULL')
      ELSE
        CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PROT')
      END IF
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  STRUCTURE FILE   =====================================================
      ID=+'STRC'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.STRC')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  STRUCTURE FILE   =====================================================
      ID=+'STRC_OUT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.STRC_OUT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  MOVIE FILE FOR QMMM MULTIPLE TIMESTEPS ===============================
      CALL FILEHANDLER$SETFILE('MMMOVIE',.TRUE.,-'.MMMOVIE.XYZ')
      CALL FILEHANDLER$SETSPECIFICATION('MMMOVIE','STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION('MMMOVIE','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('MMMOVIE','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('MMMOVIE','FORM','FORMATTED')
!
!     ==  MOVIE FILE FOR QMMM MULTIPLE TIMESTEPS ===============================
      CALL FILEHANDLER$SETFILE('SMOVIE',.TRUE.,-'.SMOVIE.XYZ')
      CALL FILEHANDLER$SETSPECIFICATION('SMOVIE','STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION('SMOVIE','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('SMOVIE','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('SMOVIE','FORM','FORMATTED')
!
!     ==  SCRAP FILE FOR TEST PURPOSES, THAT ALLOWS TO WRITE OUT TEST DATA======
      ID=+'INFO'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.INFO')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  ENTRY TO WHICH DIFFERENT FILES CAN BE TEMPORARILY ATTACHED ===========
      ID=+'HOOK'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.FORGOTTOASSIGNFILETOHOOKERROR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  RESTART_IN FILE ======================================================
      ID=+'RESTART_IN'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.RSTRT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
!
!     ==  RESTART_OUT FILE =====================================================
      ID=+'RESTART_OUT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.RSTRT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READWRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
!
!     ==  ENERGY TRAJECTORY ====================================================
      CALL TRAJECTORYIO$NEW(+'ENERGY-TRAJECTORY')
      CALL TRAJECTORYIO$SELECT('ENERGY-TRAJECTORY')
      CALL TRAJECTORYIO$SETCH('FILE',-'*_E.TRA')
      CALL TRAJECTORYIO$SELECT('NONE')
!
!     ==  BANDS  TRAJECTORY ====================================================
      CALL TRAJECTORYIO$NEW(+'BANDS-TRAJECTORY')
      CALL TRAJECTORYIO$SELECT('BANDS-TRAJECTORY')
      CALL TRAJECTORYIO$SETCH('FILE',-'*_B.TRA')
      CALL TRAJECTORYIO$SELECT('NONE')
!
!     ==  POSITION TRAJECTORY ==================================================
      CALL TRAJECTORYIO$NEW(+'POSITION-TRAJECTORY')
      CALL TRAJECTORYIO$SELECT('POSITION-TRAJECTORY')
      CALL TRAJECTORYIO$SETCH('FILE',-'*_R.TRA')
      CALL TRAJECTORYIO$SELECT('NONE')
!
!     ==  QM-MM POSITION TRAJECTORY ============================================
      CALL TRAJECTORYIO$NEW(+'QMMM-POS-TRA')
      CALL TRAJECTORYIO$SELECT('QMMM-POS-TRA')
      CALL TRAJECTORYIO$SETCH('FILE',-'*_RQMMM.TRA')
      CALL TRAJECTORYIO$SELECT('NONE')
!
!     ==  FORCE TRAJECTORY =====================================================
      CALL TRAJECTORYIO$NEW(+'FORCE-TRAJECTORY')
      CALL TRAJECTORYIO$SELECT('FORCE-TRAJECTORY')
      CALL TRAJECTORYIO$SETCH('FILE',-'*_F.TRA')
      CALL TRAJECTORYIO$SELECT('NONE')
!
!     ==  CONSTRAINTS TRAJECTORY ===============================================
      ID=+'CONSTRAINTS'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'_CONSTR.REPORT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PDOS  ================================================================
      ID=+'PDOS'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.PDOS')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
!
!     ==  BANDDATA==============================================================
      ID=+'BANDDATA'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.BANDDATA')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','UNFORMATTED')
!
!     ==  BALLSTICK ============================================================
      ID=+'BALLSTICK.DX'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'_BALLSTICK.DX')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  DENSITY  =============================================================
      ID=+'DENSITY.DX'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'_RHO.DX')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  DENSITY  =============================================================
      ID=+'POTENTIAL.DX'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'_POT.DX')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  DENSITY  =============================================================
      ID=+'PATH.DX'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'_PATH.DX')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
                                   CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_GENERIC(LL_CNTL_,NBEG,NOMORE,IPRINT,DELT)
!     **************************************************************************
!     ** 
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      REAL(8)      ,INTENT(OUT):: DELT
      INTEGER(4)   ,INTENT(OUT):: NOMORE
      INTEGER(4)   ,INTENT(OUT):: IPRINT
      INTEGER(4)   ,INTENT(OUT):: NBEG
      TYPE(LL_TYPE)            :: LL_CNTL
      LOGICAL(4)               :: TEQ
      LOGICAL(4)               :: TCHK
      CHARACTER(32)            :: CH32SVAR
      INTEGER(4)               :: LEN
      INTEGER(4)               :: ISVAR
      REAL(8)                  :: SVAR
      INTEGER(4)               :: RUNTIME(3)
!     **************************************************************************
                           CALL TRACE$PUSH('READIN_GENERIC')
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
!
!     ==========================================================================
!     == TAKE CARE OF TRACE FIRST ==============================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'TRACE',0,TCHK)
      IF(TCHK) THEN
         CALL LINKEDLIST$GET(LL_CNTL,'TRACE',1,TCHK)
         CALL TRACE$SETL4('ON',TCHK)
      ELSE
         CALL TRACE$SETL4('ON',.FALSE.)
      END IF
!
!     == PRECONDITION FILEHANDLER FOR LITTLE AND BIG ENDIAN =============
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ENDIAN',0,TCHK)
      IF(TCHK) THEN
         CALL LINKEDLIST$GET(LL_CNTL,'ENDIAN',1,CH32SVAR)
         IF(CH32SVAR.EQ.'INTEL') THEN
           CH32SVAR='LITTLE'
         ELSE IF(CH32SVAR.EQ.'IBM') THEN
           CH32SVAR='BIG'
         END IF
         IF(CH32SVAR.NE.'LITTLE'.AND.CH32SVAR.NE.'BIG') THEN
           CALL ERROR$MSG('ENDIAN MUST HAVE VALUES LITTLE,BIG,INTEL OR IBM')
           CALL ERROR$CHVAL('ENDIAN',CH32SVAR)
           CALL ERROR$STOP('READIN IN !CNTL!GENERIC:ENDIAN')
         END IF
         CALL FILEHANDLER$SETCH('ENDIAN',CH32SVAR)
      END IF
!
!     == DEFAULT VALUES
      CALL LINKEDLIST$EXISTD(LL_CNTL,'DT',0,TCHK)
      IF(.NOT.TCHK) CALL LINKEDLIST$SET(LL_CNTL,'DT',0,5.D0) 
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NSTEP',0,TCHK)
      IF(.NOT.TCHK) CALL LINKEDLIST$SET(LL_CNTL,'NSTEP',0,100) 
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NWRITE',0,TCHK)
      IF(.NOT.TCHK) CALL LINKEDLIST$SET(LL_CNTL,'NWRITE',0,100) 
      CALL LINKEDLIST$EXISTD(LL_CNTL,'START',0,TCHK)
      IF(.NOT.TCHK) CALL LINKEDLIST$SET(LL_CNTL,'START',0,.FALSE.) 
!
!     == CAPTURE USE OF OUTDATED SYNTAX ============================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STOREPSIR',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('SYNTAX HAS CHANGED')
        CALL ERROR$MSG('THIS OPTION HAS BEEN ELIMINATED')
        CALL ERROR$STOP('READIN')
      END IF
!
!     ==  READ ACTUAL VALUES  ==========================================
      CALL LINKEDLIST$GET(LL_CNTL,'DT',1,DELT)
      CALL LINKEDLIST$GET(LL_CNTL,'NSTEP',1,NOMORE)
      CALL LINKEDLIST$GET(LL_CNTL,'NWRITE',1,IPRINT)
      CALL LINKEDLIST$GET(LL_CNTL,'START',1,TEQ)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RSTRTTYPE',0,TCHK)
      IF(.NOT.TCHK) CALL LINKEDLIST$SET(LL_CNTL,'RSTRTTYPE',0,'DYNAMIC') 
      CALL LINKEDLIST$GET(LL_CNTL,'RSTRTTYPE',1,CH32SVAR)
      IF(CH32SVAR.NE.'STATIC'.AND.CH32SVAR.NE.'DYNAMIC') THEN
        CALL ERROR$MSG('ILLEGAL VALUE FOR RSTRTTYPE')
        CALL ERROR$CHVAL('RSTRTTYPE',CH32SVAR)
        CALL ERROR$MSG('ALLOWED VALUES ARE "DYNAMIC" AND "STATIC"')
        CALL ERROR$STOP('READIN')
      END IF
      CALL WAVES$SETCH('RSTRTTYPE',CH32SVAR)
!
!     ================================================================
      IF(TEQ) THEN
        NBEG=-1
      ELSE
        NBEG=0
      END IF
!
!     ==========================================================================
!     ==  READ RUNTIME                                                        ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RUNTIME',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SIZE(LL_CNTL,'RUNTIME',1,LEN)
        IF(LEN.LT.1.OR.LEN.GT.3) THEN
          CALL ERROR$MSG('ILLEGAL LENGTH FOR RUNTIME! ALLOWED LENGTH: ')
          CALL ERROR$MSG('1 - SECONDS')
          CALL ERROR$MSG('2 - MINUTES SECOUNDS')
          CALL ERROR$MSG('3 - HOURS MINUTES SECOUNDS')
          CALL ERROR$I4VAL('LENGTH:',LEN)
          CALL ERROR$STOP('READIN')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'RUNTIME',1,RUNTIME(1:LEN))
        IF(LEN.EQ.3)RUNTIME(1)=RUNTIME(3)+60*(RUNTIME(2)+60*RUNTIME(1))
        IF(LEN.EQ.2)RUNTIME(1)=RUNTIME(2)+60*RUNTIME(1)
        IF(LEN.EQ.1)RUNTIME(1)=RUNTIME(1)
        CALL STOPIT$SETI4('RUNTIME',RUNTIME(1))
      END IF
!
!     ==========================================================================
!     ==  READ AUTOCONV                                                       ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'AUTOCONV',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'AUTOCONV',1,ISVAR)
        IF(ISVAR.LT.1) THEN
          CALL ERROR$MSG('!CONTROL!GENERIC:AUTONV<1 NOT ALLOWED')
          CALL ERROR$I4VAL('AUTOCONV:',ISVAR)
          CALL ERROR$STOP('READIN')
        END IF
        CALL AUTO$SETI4('NTOL',ISVAR)
      END IF
!
!     ==========================================================================
!     ==  READ ETOL                                                           ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ETOL',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'ETOL',1,SVAR)
        CALL AUTO$SETR8('TOLERANCE',SVAR)
      END IF
!
!     ==========================================================================
!     ==  RESET ATOMIC STRUCTURE FROM STRUCTURE INPUT FILE                    ==
!     ==  I.E. IGNORE UNIT CELL AND ATOMIC OPOSITIONS FROM RESTART FILE       ==
!     == CAUTION! REDUNDANCY WITH !CNTL!RDYN:START                            ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NEWSTRC',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'NEWSTRC',1,TCHK)
        CALL CELL$SETL4('START',TCHK)
        CALL ATOMS$SETL4('START',TCHK)
      END IF
                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_DFT(LL_CNTL_)
!     **************************************************************************
!     **  READ INPUT DATA FROM THE !CONTROL!DFT BLOCK OF THE CONTROL FILE     **
!     **                                                                      **
!     **  CAUTION: VARIABLE CHTYPE OCCURS WITH SAME NAME IN LINKEDLIST_MODULE **
!     **           THAT IS THE REASON THE USE STATEMENT BEING VERY EXPLICIT   **
!     **************************************************************************
      USE LINKEDLIST_MODULE, ONLY : LL_TYPE &
     &                             ,LINKEDLIST$GET &
     &                             ,LINKEDLIST$EXISTD &
     &                             ,LINKEDLIST$EXISTL &
     &                             ,LINKEDLIST$SELECT &
     &                             ,LINKEDLIST$SET &
     &                             ,LINKEDLIST$SIZE 
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      INTEGER(4)               :: ILDA
      LOGICAL(4)               :: TCHK,TCHK1,TCHK2
      REAL(8)                  :: SVAR
      CHARACTER(32)            :: MODUS
      CHARACTER(32)            :: CHTYPE(3)  ! LIBXC IDENTIFIERS
      INTEGER(4)               :: LEN
      REAL(8)                  :: ANGSTROM  ! ANGSTROM
!     **************************************************************************
                          CALL TRACE$PUSH('READIN_DFT')
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
!
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DFT')
!
!     ==========================================================================
!     == FUNCTIONAL TYPE                                                      ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'LIBXC',1,TCHK)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'TYPE',1,TCHK1)
      IF(TCHK.AND.TCHK1) THEN
        CALL ERROR$MSG('TYPE AND LIBXC ARE MUTUALLY EXCLUSIVE')
        CALL ERROR$MSG('REMOVE ONE OF THEM FROM !CONTROL!DFT')
        CALL ERROR$STOP('READIN_DFT')
      END IF
!
!     == SET DEFAULT (PERDEW BURKE ERNZERHOF GGA) ==============================
      IF((.NOT.TCHK).AND.(.NOT.TCHK1)) THEN
        ILDA=10   ! PBE FUNCTIONAL, INTRINSIC IMPLEMENTATION 
        CALL DFT$SETI4('TYPE',ILDA)
      END IF
!
!     == READ FUNCTIONAL SPECIFICATION =========================================
      IF(TCHK) THEN
!       -- READ LIBXC INPUT ----------------------------------------------------
        CALL LINKEDLIST$SIZE(LL_CNTL,'LIBXC',1,LEN)
        IF(LEN.GT.SIZE(CHTYPE)) THEN
          CALL ERROR$MSG('NUMBER OF LIBXC IDENTIFIERS EXCEEDS MAXIMUM')
          CALL ERROR$STOP('READIN_DFT')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'LIBXC',1,CHTYPE(1:LEN))
        CALL DFT$SETCHA('TYPE',LEN,CHTYPE(1:LEN))
      ELSE IF(TCHK1) THEN
!       -- USE INTRINSIC FUNCTIONALS W/O LIBXC ---------------------------------
        CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,ILDA)
        CALL DFT$SETI4('TYPE',ILDA)
      END IF
!
!     ==========================================================================
!     == VAN DER WAALS INTERFACE                                              ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'VDW',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'VDW',1,TCHK)
      IF(TCHK) THEN
        CALL VDW$SETL4('ON',TCHK)
!       == SET ISOLATE OBJECT TO ON. DECOUPLE IS OFF PER DEFAULT 
        CALL ISOLATE$SETL4('ON',.TRUE.)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'VDW-3BODY',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'VDW-3BODY',1,TCHK)
        CALL VDW$SETL4('E63',TCHK)
        IF(ILDA.EQ.10) THEN
          CALL VDW$SETCH('FUNCTIONAL','PBE')
        ELSE
          CALL ERROR$MSG('GRIMMES VAN DER WAALS CORRECTION')
          CALL ERROR$MSG('INCLUDED ONLY FOR PBE0 FUNCTIONAL (ILDA=10)')
          CALL ERROR$STOP('READIN_DFT')
        END IF
      END IF
!
!     ==========================================================================
!     == NTBO INTERFACE                                                       ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTL(LL_CNTL,'NTBO',1,TCHK)
      IF(TCHK) THEN
!        CALL LMTO$SETL4('ON',.TRUE.)
CALL LMTO$SETL4('ON',.FALSE.)
        CALL SIMPLELMTO$SETL4('ON',.TRUE.)
        CALL LINKEDLIST$SELECT(LL_CNTL,'NTBO')
!
!       == MODUS  (CAN BE 'HYBRID', 'DMFT', OR 'ROBERT') ======================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'MODUS',1,TCHK1)
        MODUS='HYBRID' ! DEFAULT
        IF(TCHK1)CALL LINKEDLIST$GET(LL_CNTL,'MODUS',1,MODUS)
        IF(MODUS.NE.'HYBRID'.AND.MODUS.NE.'DMFT'.AND.MODUS.NE.'ROBERT') THEN
          CALL ERROR$MSG('INVALID VALUE FOR VARIABLE !CONTROL!DFT!NTBO:MODUS')
          CALL ERROR$MSG('ALLOWED VALUES ARE "HYBRID", "DMFT", "ROBERT"')
          CALL ERROR$CHVAL('MODUS',MODUS)
          CALL ERROR$STOP('READIN_DFT')
        END IF
        CALL LMTO$SETCH('MODUS',MODUS) 
        CALL SIMPLELMTO$SETCH('MODUS',MODUS) 
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'OFFSITE',1,TCHK1)
        IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'OFFSITE',1,TCHK2)
          IF(TCHK2.AND.MODUS.NE.'HYBRID') THEN
            CALL ERROR$MSG('INCONSISTENT SETTING')
            CALL ERROR$MSG('OFFSITE=T IS COMPATIBLE ONLY WITH MODUS "HYBRID"')
            CALL ERROR$L4VAL('OFFSITE',TCHK2)
            CALL ERROR$CHVAL('MODUS',MODUS)
            CALL ERROR$STOP('READIN_DFT')
          END IF
          CALL LMTO$SETL4('OFFSITE',TCHK2)
          CALL SIMPLELMTO$SETL4('OFFSITE',TCHK2)
        END IF
!
!!$        CALL LINKEDLIST$EXISTD(LL_CNTL,'DROP',1,TCHK1)
!!$        IF(TCHK1) THEN
!!$          CALL ERROR$MSG('OPTION DROP IS OBSOLETE')
!!$          CALL ERROR$STOP('READIN_DFT')
!!$!
!!$          CALL LINKEDLIST$GET(LL_CNTL,'DROP',1,TCHK2)
!!$          IF(TCHK2.AND.MODUS.NE.'OLDDMFT') THEN
!!$            CALL ERROR$MSG('INCONSISTENT SETTINGE')
!!$            CALL ERROR$MSG('DROP=T IS COMPATIBLE ONLY WITH MODUS "OLDDMFT"')
!!$            CALL ERROR$L4VAL('DROP',TCHK2)
!!$            CALL ERROR$CHVAL('MODUS',MODUS)
!!$            CALL ERROR$STOP('READIN_DFT')
!!$          END IF
!!$          CALL LMTO$SETL4('DROP',TCHK2)
!!$        END IF
!!$!
!!$        CALL LINKEDLIST$EXISTD(LL_CNTL,'PICK',1,TCHK1)
!!$        IF(TCHK1) THEN
!!$          CALL ERROR$MSG('OPTION PICK IS OBSOLETE')
!!$          CALL ERROR$STOP('READIN_DFT')
!!$          CALL LINKEDLIST$GET(LL_CNTL,'PICK',1,TCHK2)
!!$          IF(TCHK2.AND.MODUS.NE.'OLDDMFT') THEN
!!$            CALL ERROR$MSG('INCONSISTENT SETTINGE')
!!$            CALL ERROR$MSG('PICK=T IS COMPATIBLE ONLY WITH MODUS "OLDDMFT"')
!!$            CALL ERROR$L4VAL('DROP',TCHK2)
!!$            CALL ERROR$CHVAL('MODUS',MODUS)
!!$            CALL ERROR$STOP('READIN_DFT')
!!$          END IF
!!$          CALL LMTO$SETL4('PICK',TCHK2)
!!$          CALL LINKEDLIST$EXISTD(LL_CNTL,'DHOFK',1,TCHK2)
!!$          IF(TCHK2) THEN
!!$            CALL ERROR$MSG('OPTION DHOFK IS OBSOLETE')
!!$            CALL ERROR$STOP('READIN_DFT')
!!$            CALL LINKEDLIST$GET(LL_CNTL,'DHOFK',1,TCHK2)
!!$            CALL LMTO$SETL4('DHOFK',TCHK2)
!!$          END IF
!!$        END IF
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'HFWEIGHT',1,TCHK)
        IF(TCHK) THEN
          CALL ERROR$MSG('ENTRY HFWEIGHT IS OBSOLETE. CONSULT MANUAL')
          CALL ERROR$STOP('READIN_DFT')
        END IF
!!$
        CALL LINKEDLIST$EXISTD(LL_CNTL,'K2',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'K2',1,SVAR)
          CALL LMTO$SETR8('K2',SVAR)
          CALL SIMPLELMTO$SETR8('K2',SVAR)
        END IF
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'SCREENL[AA]',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'SCREENL[AA]',1,SVAR)
          IF(SVAR.GE.0.D0) THEN
            IF(SVAR.LT.2.D-2) THEN
!             __THIS AVOIDS UNDERFLOW IN SLATEC ROUTINE CALLED BY_______________
!             __RADIAL$YUKAWA IN PAW_SIMPLELMTO.F90_____________________________
              WRITE(*,*)'WARNING! SCREENL[AA] HAS BEEN CHANGED TO 0.02'
              SVAR=2.D-2
            END IF
          ELSE
            CALL ERROR$MSG('NEGATIVE VALUES FOR SCREENL[AA] ARE NOT ALLOWED.')
            CALL ERROR$MSG('TO SWITCH SCREENING OFF, DISABLE SCREENL[AA].')
            CALL ERROR$MSG('(!CONTROL!DFT!NTBO:SCREENL[AA])')
            CALL ERROR$R8VAL('SCREENL[AA]',SVAR)
            CALL ERROR$STOP('READIN_DFT')
          END IF
          SVAR=SVAR*ANGSTROM
          CALL SIMPLELMTO$SETR8('SCREENL',SVAR)
        END IF
!!$
        CALL LINKEDLIST$EXISTD(LL_CNTL,'SCALERCUT',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'SCALERCUT',1,SVAR)
          CALL LMTO$SETR8('SCALERCUT',SVAR)
          CALL WAVES$SETR8('SCALERCUT',SVAR)
        END IF

        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      END IF

                          CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_FOURIER(LL_CNTL_)
!     ******************************************************************
!     ** MODIFIES     : WAVES;POTENTIAL                               ** 
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK,TCHK1
      REAL(8)               :: EPWPSI    ! WAVE FUNCTION PLANE WAVE CUTOFF
      REAL(8)               :: EPWRHO    ! DENSITY PLANE WAVE CUTOFF
      REAL(8)               :: EPWBUCKET   ! MIN ENERGY FOR BUCKET POTENTIAL
      REAL(8)               :: BUCKETPAR ! PARAMETER FOR BUCKET POTENTIAL
      REAL(8)               :: CDUAL     ! EPWRHO=EPWPSI*CDUAL
      REAL(8)               :: RY        ! RYDBERG
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_FOURIER')
      LL_CNTL=LL_CNTL_
      CALL CONSTANTS('RY',RY)
!
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'FOURIER')
!
!     == EPWPSI ========================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EPWPSI',1,TCHK)
      IF(.NOT.TCHK) CALL LINKEDLIST$SET(LL_CNTL,'EPWPSI',0,30.D0)
      CALL LINKEDLIST$CONVERT(LL_CNTL,'EPWPSI',1,'R(8)')
      CALL LINKEDLIST$GET(LL_CNTL,'EPWPSI',1,EPWPSI)
      EPWPSI=EPWPSI*RY
      CALL WAVES$SETR8('EPWPSI',EPWPSI)
!
!     == EPWRHO ========================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EPWRHO',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$CONVERT(LL_CNTL,'EPWRHO',1,'R(8)')
        CALL LINKEDLIST$GET(LL_CNTL,'EPWRHO',1,EPWRHO)
        EPWRHO=EPWRHO*RY
      ELSE
        CALL LINKEDLIST$EXISTD(LL_CNTL,'CDUAL',1,TCHK)
        IF(.NOT.TCHK) THEN
          CDUAL=4.D0
          CALL LINKEDLIST$SET(LL_CNTL,'CDUAL',0,CDUAL)
        END IF
        CALL LINKEDLIST$CONVERT(LL_CNTL,'CDUAL',1,'R(8)')
        CALL LINKEDLIST$GET(LL_CNTL,'CDUAL',1,CDUAL)
        EPWRHO=EPWPSI*CDUAL
      END IF
      CALL WAVES$SETR8('EPWRHO',EPWRHO)
!
!     == EPWBUCKET =====================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EPWBUCKET',1,TCHK)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'BUCKETPAR',1,TCHK1)
      IF(TCHK.AND.TCHK1) THEN
        CALL LINKEDLIST$CONVERT(LL_CNTL,'EPWBUCKET',1,'R(8)')
        CALL LINKEDLIST$GET(LL_CNTL,'EPWBUCKET',1,EPWBUCKET)
        EPWBUCKET=EPWBUCKET*RY
        CALL LINKEDLIST$CONVERT(LL_CNTL,'BUCKETPAR',1,'R(8)')
        CALL LINKEDLIST$GET(LL_CNTL,'BUCKETPAR',1,BUCKETPAR)
        CALL WAVES$SETR8('EPWBUCKET',EPWBUCKET)
        CALL WAVES$SETR8('BUCKETPAR',BUCKETPAR)
      ELSE IF(TCHK.NEQV.TCHK1) THEN
        CALL ERROR$MSG('BUCKET POTENTIAL REQUIRES EPWBUCKET AND BUCKETPAR')
        CALL ERROR$L4VAL('EPWBUCKET PRESENT',TCHK)
        CALL ERROR$L4VAL('BUCKETPAR PRESENT',TCHK1)
        CALL ERROR$STOP('READIN_FOURIER')
      END IF
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_PSIDYN(LL_CNTL_)
!     **************************************************************************
!     ** READ INFORMATION FOR WAVE FUNCTION DYNAMICS FROM CONTROL FILE        **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      TYPE(LL_TYPE)         :: LL_TEST
      LOGICAL(4)            :: TCHK,TCHK1,TCHK2
      REAL(8)               :: FRICTION
      LOGICAL(4)            :: TSTOPE
      LOGICAL(4)            :: TRANE
      REAL(8)               :: AMPRE
      REAL(8)               :: EMASS
      REAL(8)               :: EMASSCG2
      REAL(8)               :: DT
      REAL(8)               :: PI
!     **************************************************************************
                           CALL TRACE$PUSH('READIN_PSIDYN')  
      LL_CNTL=LL_CNTL_
      PI=4.D0*ATAN(1.D0)
!
!     ==========================================================================
!     ==  READ BLOCK GENERIC                                                  ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
      CALL LINKEDLIST$GET(LL_CNTL,'DT',1,DT)
      CALL WAVES$SETR8('TIMESTEP',DT)
!
!     ==========================================================================
!     ==  READ BLOCK !CONTROL!PSIDYN                                          ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'PSIDYN')
!
!     ==========================================================================
!     == CAPTURE OUTDATED VERSIONS  ============================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'DIAG',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('SYNTAX HAS CHANGED')
        CALL ERROR$MSG('OPTION DIAG HAS BEEN REMOVED')
        CALL ERROR$STOP('READIN$PSIDYN')
      END IF
!
!     ==========================================================================
!     ==  BEGIN WITH ZERO VELOCITIES ===========================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STOP',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'STOP',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'STOP',1,TSTOPE)
      CALL WAVES$SETL4('STOP',TSTOPE)
!
!     ==========================================================================
!     ==  BEGIN WITH RANDOM VELOCITIES =========================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RANDOM',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'RANDOM',0,0.D0)
      CALL LINKEDLIST$CONVERT(LL_CNTL,'RANDOM',1,'R(8)')
      CALL LINKEDLIST$GET(LL_CNTL,'RANDOM',1,AMPRE)
      TRANE=(AMPRE.NE.0.D0)
      CALL WAVES$SETL4('RANDOMIZE',TRANE)
      CALL WAVES$SETR8('AMPRE',AMPRE)
!
!     ==========================================================================
!     ==  FRICTION  (MAY BE CHANGED BY !AUTO)                                 ==
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FRIC',0,0.D0)
      CALL LINKEDLIST$CONVERT(LL_CNTL,'FRIC',1,'R(8)')
      CALL LINKEDLIST$GET(LL_CNTL,'FRIC',1,FRICTION)
      CALL WAVES$SETR8('FRICTION',FRICTION)
!
!     ==========================================================================
!     ==  FICTITIOUS WAVE FUNCTION MASS                                       ==
      CALL LINKEDLIST$EXISTD(LL_CNTL,'MPSI',1,TCHK)
      IF(.NOT.TCHK) THEN
        EMASS=10.D0*DT**2
        CALL LINKEDLIST$SET(LL_CNTL,'MPSI',0,EMASS)
      ELSE
        CALL LINKEDLIST$CONVERT(LL_CNTL,'MPSI',1,'R(8)')
        CALL LINKEDLIST$GET(LL_CNTL,'MPSI',1,EMASS)
      END IF
      CALL WAVES$SETR8('EMASS',EMASS)
!
!     ==========================================================================
!     ==  ENHANCEMENT FACTOR FOR FICTITIOUS WAVE FUNCTION MASS                ==
      CALL LINKEDLIST$EXISTD(LL_CNTL,'MPSICG2',1,TCHK)
      IF(.NOT.TCHK) THEN
         EMASSCG2=0.5D0*(10.D0/(2.D0*PI))**2*DT**2/EMASS
         CALL LINKEDLIST$SET(LL_CNTL,'MPSICG2',0,EMASSCG2)
!        CALL LINKEDLIST$SET(LL_CNTL,'MPSICG2',0,0.D0)
      END IF
      CALL LINKEDLIST$CONVERT(LL_CNTL,'MPSICG2',1,'R(8)')
      CALL LINKEDLIST$GET(LL_CNTL,'MPSICG2',1,EMASSCG2)
      CALL WAVES$SETR8('EMASSCG2',EMASSCG2)
!
!     ==========================================================================
!     ==  STORE REAL SPACE WAVE FUNCTIONS TO AVOID ADDITIONAL FFT ==============
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STOREPSIR',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('THE PARAMETER STOREPSIR HAS BEEN ELIMINATED')
        CALL ERROR$STOP('READIN_WAVES')         
      END IF
!
!     ==========================================================================
!     ==  SAFEORTHO=T STRICTLY CONSERVES ENERGY,                              ==
!     ==  BUT DOES NOT PRODUCE EIGENSTATES =====================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'SAFEORTHO',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'SAFEORTHO',0,.TRUE.)
      CALL LINKEDLIST$GET(LL_CNTL,'SAFEORTHO',1,TCHK)
      CALL WAVES$SETL4('SAFEORTHO',TCHK)
!
!     ==========================================================================
!     ==  STRAIGHTEN MAKES ONCE A TRANSFORMATION TO EIGENSTATES               ==
!     ==                                                                      ==
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STRAIGHTEN',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'STRAIGHTEN',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'STRAIGHTEN',1,TCHK)
      CALL WAVES$SETL4('STRAIGHTEN',TCHK)
!
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'SWAPSTATES',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'SWAPSTATES',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'SWAPSTATES',1,TCHK)
      CALL WAVES$SETL4('SWAPSTATES',TCHK)
!
!     ==========================================================================
!     ==  READ       FIXRHO                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FIXRHO',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FIXRHO',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'FIXRHO',1,TCHK)
      CALL WAVES$SETL4('FIXRHO',TCHK)
      IF(TCHK) THEN
        LL_TEST=LL_CNTL
        CALL LINKEDLIST$SELECT(LL_TEST,'~')
        CALL LINKEDLIST$SELECT(LL_TEST,'CONTROL')
        CALL LINKEDLIST$EXISTL(LL_TEST,'RDYN',1,TCHK1)
        IF(TCHK1) THEN
          CALL ERROR$MSG('FROZEN POTENTIAL AND ATOM DYNAMICS ARE INCOMPATIBLE')
          CALL ERROR$MSG('!CONTROL!RDYN AND !CONTROL!PSIDYN:FIXRHO')
          CALL ERROR$MSG('ARE INCOMPATIBLE')
          CALL ERROR$STOP('READIN_PSIDYN')
        END IF
      END IF
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'WRITERHO',0,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'WRITERHO',1,TCHK)
        CALL WAVES$SETL4('WRITERHO',TCHK)
      END IF
!
!     ==========================================================================
!     ==  READ BLOCK !CONTROL!PSIDYN!AUTO                                     ==
!     ==========================================================================
      CALL READIN_AUTO(LL_CNTL,'WAVES',0.3D0,0.97D0,0.3D0,1.D0)
!
!     ==========================================================================
!     ==  READ BLOCK !CONTROL!PSIDYN!THERMOSTAT                               ==
!     ==========================================================================
      CALL READIN_THERMOSTAT(LL_CNTL,'WAVES',DT,1.D0,100.D0)
!
!     ==========================================================================
!     ==  CHECK FOR SIMULTANEOUS USE OF AUTO AND THERMOSTAT                   ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'THERMOSTAT',1,TCHK1)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'AUTO',1,TCHK2)
      IF(TCHK1.AND.TCHK2) THEN
        CALL ERROR$MSG('!AUTO AND !THERMOSTAT ARE MUTUALLY EXCLUSIVE')
        CALL ERROR$STOP('READIN_PSIDYN')
      END IF
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_RDYN(LL_CNTL_)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      LOGICAL(4)               :: TCHK,TCHK1,TCHK2
      REAL(8)                  :: FRICTION
      LOGICAL(4)               :: TFOR
      LOGICAL(4)               :: TSTOPR
      REAL(8)                  :: AMPRP
      LOGICAL(4)               :: TRANP
      REAL(8)                  :: CELVIN
      REAL(8)                  :: DT
!     **************************************************************************
                            CALL TRACE$PUSH('READIN_RDYN')
      LL_CNTL=LL_CNTL_
!
!     ==========================================================================
!     ==  READ BLOCK GENERIC                                                  ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
      CALL LINKEDLIST$GET(LL_CNTL,'DT',1,DT)
      CALL WAVES$SETR8('TIMESTEP',DT)
!
!     ==========================================================================
!     ==  GET INTO BLOCK !CONTROL!RDYN                                        ==
!     ==========================================================================
      CALL CONSTANTS('KB',CELVIN)
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'RDYN',1,TFOR)
      CALL LINKEDLIST$SELECT(LL_CNTL,'RDYN')
!
!     == SET ATOMS OBJECT ======================================================
      CALL ATOMS$SETL4('MOVE',TFOR)
!
!     == STOP ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STOP',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'STOP',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'STOP',1,TSTOPR)
      CALL ATOMS$SETL4('STOP',TSTOPR)
!
!     == RANDOMIZE =============================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RANDOM[K]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'RANDOM[K]',0,0.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'RANDOM[K]',1,AMPRP) 
      AMPRP=AMPRP*CELVIN
      TRANP=(AMPRP.NE.0.D0)
      CALL ATOMS$SETL4('RANDOMIZE',TRANP)
      CALL ATOMS$SETR8('AMPRE',AMPRP)
!
!     == FRICTION ==============================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FRIC',0,0.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'FRIC',1,FRICTION)
!     __ FRICTION MAY HAVE AN INFLUENCE ON THE ELECTRON FRICTION________________
!     __ EVEN IF ATOMS DO NOT MOVE______________________________________________
      IF(.NOT.TFOR) FRICTION=0.D0
      CALL ATOMS$SETR8('FRICTION',FRICTION)
!
!     == START WITH POSITIONS FROM STRUCTURE INPUT FILE 'STRC' =================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'START',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('OPTION !CONTROL!RDYN:START IS OBSOLETE')
        CALL ERROR$MSG('USE !CONTROL!GENERIC:NEWSTRC INSTEAD')
        CALL ERROR$STOP('READIN_RDYN')
!        CALL LINKEDLIST$GET(LL_CNTL,'START',1,TCHK)
      END IF
!      CALL ATOMS$SETL4('START',TCHK)
!
!     ==========================================================================
!     ==  ATTENTION                                                           ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'USEOPTFRIC',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'USEOPTFRIC',1,TCHK)
      END IF
      CALL ATOMS$SETL4('USEOPTFRIC',TCHK)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NONEGFRIC',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'NONEGFRIC',1,TCHK)
      END IF
      CALL ATOMS$SETL4('NONEGATIVEFRICTION',TCHK)
!
!     ==========================================================================
!     ==  READ BLOCK !CONTROL!RDYN!AUTO                                       ==
!     ==========================================================================
      CALL READIN_AUTO(LL_CNTL,'ATOMS',0.D0,1.D0,1.D-2,1.D0)
!
!     ==========================================================================
!     ==  READ BLOCK !CONTROL!RDYN!THERMOSTAT                                 ==
!     ==========================================================================
      CALL READIN_THERMOSTAT(LL_CNTL,'ATOMS',DT,293.D0*CELVIN,10.D0)
!
!     ==========================================================================
!     ==  READ BLOCK !CONTROL!RDYN!WARMUP                                     ==
!     ==========================================================================
      CALL READIN_RDYN_WARMUP(LL_CNTL)
!
!     ==========================================================================
!     ==  CHECK FOR SIMULTANEOUS USE OF AUTO AND THERMOSTAT                   ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'THERMOSTAT',1,TCHK1)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'AUTO',1,TCHK2)
      IF(TCHK1.AND.TCHK2) THEN
        CALL ERROR$MSG('AUTO AND THERMOSTAT MUST NOT')
        CALL ERROR$MSG('BE SPECIFIED SIMULTANEOUSLY')
        CALL ERROR$STOP('READIN_RDYN')
      END IF
                            CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_AUTO(LL_CNTL_,ID,ANNEL,FACL,ANNEU,FACU)
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      CHARACTER(*) ,INTENT(IN) :: ID       ! IDENTIFIER
      REAL(8)      ,INTENT(IN) :: ANNEL    ! DEFAULT: LOWER FRICTION
      REAL(8)      ,INTENT(IN) :: FACL     ! DEFAULT: SCALE FOR LOWER FRICTION
      REAL(8)      ,INTENT(IN) :: ANNEU    ! DEFAULT: UPPER FRICTION
      REAL(8)      ,INTENT(IN) :: FACU     ! DEFAULT: SCALE FOR UPPER FRICTION
      TYPE(LL_TYPE)            :: LL_CNTL
      LOGICAL(4)               :: TCHK
      REAL(8)                  :: FRICTION
      REAL(8)                  :: SCALE
      REAL(8)      ,PARAMETER  :: ANNEMIN=0.D0
!     ******************************************************************
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$EXISTL(LL_CNTL,'AUTO',1,TCHK)
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_CNTL,'AUTO')
      CALL AUTO$SELECT(ID)
      CALL AUTO$SETL4('ON',.TRUE.)
!     
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC(-)',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FRIC(-)',0,ANNEL)
      CALL LINKEDLIST$CONVERT(LL_CNTL,'FRIC(-)',1,'R(8)')
      CALL LINKEDLIST$GET(LL_CNTL,'FRIC(-)',1,FRICTION)
      CALL AUTO$SETR8('LOWERFRICTION',FRICTION)
!     
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FACT(-)',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FACT(-)',0,FACL)
      CALL LINKEDLIST$CONVERT(LL_CNTL,'FACT(-)',1,'R(8)')
      CALL LINKEDLIST$GET(LL_CNTL,'FACT(-)',1,SCALE)
      CALL AUTO$SETR8('LOWERFACTOR',SCALE)
!     
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC(+)',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FRIC(+)',0,ANNEU)
      CALL LINKEDLIST$CONVERT(LL_CNTL,'FRIC(+)',1,'R(8)')
      CALL LINKEDLIST$GET(LL_CNTL,'FRIC(+)',1,FRICTION)
      CALL AUTO$SETR8('UPPERFRICTION',FRICTION)
!     
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FACT(+)',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FACT(+)',0,FACU)
      CALL LINKEDLIST$CONVERT(LL_CNTL,'FACT(+)',1,'R(8)')
      CALL LINKEDLIST$GET(LL_CNTL,'FACT(+)',1,SCALE)
      CALL AUTO$SETR8('UPPERFACTOR',SCALE)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'MINFRIC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'MINFRIC',0,ANNEMIN)
      CALL LINKEDLIST$CONVERT(LL_CNTL,'MINFRIC',1,'R(8)')
      CALL LINKEDLIST$GET(LL_CNTL,'MINFRIC',1,FRICTION)
      CALL AUTO$SETR8('MINFRICTION',FRICTION)
!     
      CALL AUTO$SELECT('~')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_THERMOSTAT(LL_CNTL_,ID,DT,TARGET_,FREQ_)
!     **************************************************************************
!     ** CREATES A THERMOSTAT-INSTANCE NAMED BY ID AND READS SETTING          **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      CHARACTER(*) ,INTENT(IN) :: ID        ! ID FOR THE THERMOSTAT OBJECT
      REAL(8)      ,INTENT(IN) :: DT        ! TIME STEP
      REAL(8)      ,INTENT(IN) :: TARGET_   ! TARGET KINETIC ENERGY
      REAL(8)      ,INTENT(IN) :: FREQ_     ! TARGET FREQUENCY 
      TYPE(LL_TYPE)            :: LL_CNTL   ! LINKED LIST HOLDING INPUT DATA
      LOGICAL(4)               :: TCHK,TCHK1,TCHK2
      REAL(8)                  :: TERA
      REAL(8)                  :: SECOND
      REAL(8)                  :: CELVIN
      LOGICAL(4)               :: TON
      REAL(8)                  :: TARGET
      LOGICAL(4)               :: TSTOP
      REAL(8)                  :: FREQ
      REAL(8)                  :: MASS
      REAL(8)                  :: PERIOD
      REAL(8)                  :: FRICTION
!     **************************************************************************
                            CALL TRACE$PUSH('READIN_THERMOSTAT')
!
!     ==========================================================================
!     ==  INITIALIZE CONSTANTS                                                ==
!     ==========================================================================
      LL_CNTL=LL_CNTL_
      CALL CONSTANTS('TERA',TERA)
      CALL CONSTANTS('SECOND',SECOND)
      CALL CONSTANTS('KB',CELVIN)
!
!     ==========================================================================
!     == SELECT THERMOSTAT LIST OR EXIT IF IT DOES NOT EXIST                  ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTL(LL_CNTL,'THERMOSTAT',1,TON)
      CALL THERMOSTAT$NEW(ID)
      CALL THERMOSTAT$SETL4('ON',TON)
      CALL LINKEDLIST$SELECT(LL_CNTL,'THERMOSTAT')
      IF(.NOT.TON) THEN
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==========================================================================
!     == STOP INITIAL VELOCITIES                                              ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STOP',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'STOP',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'STOP',1,TSTOP)
      CALL THERMOSTAT$SETL4('STOP',TSTOP)
!
!     ==========================================================================
!     == TEMPERATURE                                                          ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'<EKIN>',1,TCHK1)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'T[K]',1,TCHK2)
      IF(TCHK1.AND.TCHK2) THEN
        CALL ERROR$MSG('<EKIN> AND T[K] MUST NOT BE SPECIFIED SIMULTANEOUSLY')
        CALL ERROR$STOP('READIN_THERMOSTAT')
      ELSE IF(TCHK1) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'<EKIN>',1,TARGET)
      ELSE IF(TCHK2) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'T[K]',1,TARGET)
        TARGET=0.5D0*TARGET*CELVIN
      ELSE
        TARGET=TARGET_
      END IF
!
!     ==================================================================
!     ==  "EIGENFREQUENCY" OF THE THERMOSTAT                          ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'PERIOD',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'PERIOD',1,PERIOD)
      ELSE 
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FREQ[THZ]',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FREQ[THZ]',0,FREQ_)
        CALL LINKEDLIST$GET(LL_CNTL,'FREQ[THZ]',1,FREQ)
        FREQ=FREQ*TERA/SECOND
        PERIOD=1.D0/FREQ       ! CONVERT TO OSCILLATION PERIOD
      END IF
!
!     ==================================================================
!     == INITIALIZE #(DEGREES OF FREEDOM) TO ONE                      ==
!     ==================================================================
      CALL THERMOSTAT$CONVERT(DT,TARGET,PERIOD,MASS,FRICTION)
      CALL THERMOSTAT$SETR8('TIMESTEP',DT)
      CALL THERMOSTAT$SETR8('TARGET',TARGET)
      CALL THERMOSTAT$SETR8('MASS',MASS)
!
!     ==================================================================
!     == FRICTION ON THE THERMOSTAT                                   ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FRIC',0,0.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'FRIC',1,FRICTION)
      CALL THERMOSTAT$SETR8('FRICTION',FRICTION) 
!
!     ==================================================================
!     == SWITCH WAVEFUNCTION THERMOSTAT BETWEEN OLD AND NEW VERSION   ==
!     ==================================================================
      IF(ID.EQ.'WAVES') THEN
        CALL TIMESTEP$SETL4('NEWTHERMOSTAT',.TRUE.)
!       == THIS  IS IMPORTANT FOR THE NEW WAVE FUNCTION THERMOSTAT =====
        CALL THERMOSTAT$SETL4('COOLONLY',.TRUE.)
      END IF
!
                            CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_RDYN_WARMUP(LL_CNTL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK
      REAL(8)               :: K_B    ! BOLTZMAN CONSTANT
      INTEGER(4)            :: ISVAR
      REAL(8)               :: SVAR
!     ******************************************************************
      LL_CNTL=LL_CNTL_
      CALL CONSTANTS('KB',K_B)
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'RDYN')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'WARMUP',1,TCHK)
      IF(.NOT.TCHK) RETURN
!
!     ----------------------------------------
!     INSERT WORKING VERSION OF WARMUP:: TKWOO
!     ----------------------------------------

      CALL LINKEDLIST$SELECT(LL_CNTL,'WARMUP')
      CALL WARMUP$SETL4('ON',TCHK)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NPULSES',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'NPULSES',0,30)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'NPULSES',1,ISVAR)
      CALL WARMUP$SETI4('NPULSE',ISVAR)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'PULSE_LENGTH',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'PULSE_LENGTH',0,20)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'PULSE_LENGTH',1,ISVAR)
      CALL WARMUP$SETI4('PULSE_LENGTH',ISVAR)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'TARGET_TEMP',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'TARGET_TEMP',0,300.D0)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'TARGET_TEMP',1,SVAR)
      CALL WARMUP$SETR8('TARGET_TEMP',SVAR)
!
      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_MERMIN(LL_CNTL_,TON)
!     ******************************************************************
!     **  READ BLOCK !CONTROL!MERMIN                                  **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK
      REAL(8)               :: KELVIN
      REAL(8)               :: EV
      LOGICAL(4),INTENT(OUT):: TON
      REAL(8)               :: MASS
      REAL(8)               :: FRICTION
      REAL(8)               :: TEMP
      REAL(8)               :: SVAR
      REAL(8)               :: FINAL,RATE
      INTEGER(4)            :: IDIAL,NDIAL
      CHARACTER(32)         :: DIALID
      CHARACTER(1)          :: CH1
      LOGICAL(4)            :: TADIABATIC
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_MERMIN')  
      LL_CNTL=LL_CNTL_
      CALL CONSTANTS('KB',KELVIN)
      CALL CONSTANTS('EV',EV)
!
!     ==================================================================
!     == SELECT MERMIN LIST OR EXIT IF IT DOES NOT EXIST              ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'MERMIN',1,TON)
      IF(.NOT.TON) THEN
        CALL DYNOCC$SETCH('STARTTYPE','N')
        CALL DYNOCC$SETL4('DYN',.FALSE.)
        CALL TRACE$POP
        RETURN
      END IF
      CALL WAVES$GETL4('SAFEORTHO',TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('!CONTROL!MERMIN REQUIRES !CONTROL!PSIDYN:SAFEORTHO=F')
        CALL ERROR$MSG('OTHERWISE STATES ARE NOT EIGENSTATES')
        CALL ERROR$STOP('READIN_MERMIN')
      END IF
!
!     == IF OCCUPATIONS ARE DYNAMICAL THE OCCUPATIONS ARE BY DEFAULT READ
!     == FROM THE RESTART FILE UNLESS THE PARAMETER STARTTYPE IS SPECIFIED 
!     ==  ON INPUT IN THIS !CNTL!MERMIN BLOCK
      CALL LINKEDLIST$SELECT(LL_CNTL,'MERMIN')
!
!     ==================================================================
!     == COLLECT DATA                                                 ==
!     ==================================================================
      CALL DYNOCC$SETL4('DYN',TON)
!
!     ==  CONSTANT CHARGE ENSEMBLE =====================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FIXQ',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FIXQ',0,.TRUE.)
      CALL LINKEDLIST$GET(LL_CNTL,'FIXQ',0,TCHK)
      CALL DYNOCC$SETL4('FIXQ',TCHK)
!
!     ==  CONSTANT SPIN ENSEMBLE =====================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FIXS',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FIXS',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'FIXS',0,TCHK)
      CALL DYNOCC$SETL4('FIXS',TCHK)
!
!     ==  TRUE ELECTRON TEMPERATURE =================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'T[K]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'T[K]',0,1000.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'T[K]',0,TEMP)
      TEMP=TEMP*KELVIN
      CALL DYNOCC$SETR8('TEMP',TEMP)
!
!     ==  FERMI LEVEL (ONLY USED IF FIXQ=.FALSE. =======================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EFERMI[EV]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EFERMI[EV]',0,0.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'EFERMI[EV]',0,SVAR)
      SVAR=SVAR*EV
      CALL DYNOCC$SETR8('EFERMI',SVAR)
!
!     ==  MAGNETIC FIELD (ONLY USED IF FIXS=.FALSE. =======================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'MAGFIELD[EV]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'MAGFIELD[EV]',0,0.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'MAGFIELD[EV]',0,SVAR)
      SVAR=SVAR*EV
      CALL DYNOCC$SETR8('MAGNETICFIELD',SVAR)
!
!     ==  ADIABATIC ==================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ADIABATIC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ADIABATIC',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'ADIABATIC',0,TADIABATIC)
      CALL DYNOCC$SETL4('ADIABATIC',TADIABATIC)
!
!     ==  RETARD   ==================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RETARD',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'RETARD',0,1.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'RETARD',0,SVAR)
      IF(SVAR.LE.0.D0) THEN
          CALL ERROR$MSG('THE VALUE OF RETARD MUST BE POSITIVE')
          CALL ERROR$R8VAL('RETARD',SVAR)
          CALL ERROR$STOP('READIN_MERMIN')
      END IF
      CALL DYNOCC$SETR8('RETARD',SVAR)
!
!     ==  BZITYPE  ==================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'TETRA+',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'TETRA+',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'TETRA+',0,TCHK)
      IF(TCHK) THEN
        IF(.NOT.TADIABATIC) THEN
          CALL ERROR$MSG('TETRAHEDRON METHOD WORKS ONLY IN QUASI-ADIABATIC MODE')
          CALL ERROR$MSG('SELECT ADIABATIC=T OR TETRA+=F')
          CALL ERROR$STOP('READIN_MERMIN')
        END IF
        IF(TEMP.GT.KELVIN) THEN
          CALL ERROR$MSG('TETRAHEDRON METHOD WORKS ONLY FOR ZERO KELVIN')
          CALL ERROR$MSG('SELECT T[K]=0. OR TETRA+=F')
          CALL ERROR$R8VAL('T[K]',TEMP/KELVIN)
          CALL ERROR$STOP('READIN_MERMIN')
        END IF
        CALL DYNOCC$SETCH('BZITYPE','TETRA+')
      ELSE
        CALL DYNOCC$SETCH('BZITYPE','SAMP')
      END IF
!
!     == RESTART WITH OCCUPATIONS FROM STRC FILE =======================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'MOVE',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'MOVE',0,.TRUE.)
      CALL LINKEDLIST$GET(LL_CNTL,'MOVE',0,TCHK)
      CALL DYNOCC$SETL4('DYN',TCHK)
!
!     == RESTART WITH OCCUPATIONS, ENERGIES OR FROM INITIALIZATION =====
!     == STARTTYPE OVERWRITES START
      IF(TADIABATIC) THEN
        CALL DYNOCC$SETCH('STARTTYPE','E')
      ELSE
        CALL DYNOCC$SETCH('STARTTYPE','X')
      END IF
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STARTTYPE',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'STARTTYPE',0,CH1)
        IF(TADIABATIC.AND.CH1.EQ.'X') THEN
          CALL ERROR$MSG('STARTTYPE="X" IS INCOMPATIBLE WITH ADIABATIC=T')
          CALL ERROR$STOP('READIN_MERMIN')
        END IF
        CALL DYNOCC$SETCH('STARTTYPE',CH1)
      ELSE
        CALL LINKEDLIST$EXISTD(LL_CNTL,'START',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'START',0,TCHK)
          IF(TCHK)CALL DYNOCC$SETCH('STARTTYPE','E')
        END IF
      END IF
!
!     ==  FICTITIOUS MASS FOR THE OCCUPATIONS=========================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'M',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'M',0,300.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'M',0,MASS)
      CALL DYNOCC$SETR8('MASS',MASS)
!
!     ==  FRICTION ON THE OCCUPATION DYNAMICS ========================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FRIC',0,0.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'FRIC',0,FRICTION)
      CALL DYNOCC$SETR8('FRICTION',FRICTION)
!
!     ==  INITIAL VELOCITIES OF OCCUPATIONS STOPPED =================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STOP',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'STOP',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'STOP',0,TCHK)
      CALL DYNOCC$SETL4('STOP',TCHK)
!
!     ==================================================================
!     ==  RESOLVE DIALS                                               ==
!     ==================================================================
      CALL LINKEDLIST$NLISTS(LL_CNTL,'DIAL',NDIAL)
      DO IDIAL=1,NDIAL
        CALL LINKEDLIST$SELECT(LL_CNTL,'DIAL',IDIAL)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!CNTL!MERMIN!DIAL:ID NOT FOUND')
          CALL ERROR$STOP('READIN_MERMIN')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'ID',0,DIALID)
        IF(TRIM(DIALID).EQ.'T[K]') THEN
          CALL DIALS$SELECT('TEMP(E)')
          CALL LINKEDLIST$GET(LL_CNTL,'RATE',1,RATE)
          CALL LINKEDLIST$GET(LL_CNTL,'FINAL',1,FINAL)
          RATE=RATE*KELVIN
          FINAL=FINAL*KELVIN
        ELSE IF(TRIM(DIALID).EQ.'SPIN[HBAR]') THEN
          CALL DIALS$SELECT('SPIN')
          CALL LINKEDLIST$GET(LL_CNTL,'RATE',1,RATE)
          CALL LINKEDLIST$GET(LL_CNTL,'FINAL',1,FINAL)
          RATE=2.D0*RATE
          FINAL=2.D0*FINAL
        ELSE IF(TRIM(DIALID).EQ.'CHARGE[E]') THEN
          CALL DIALS$SELECT('CHARGE')
          CALL LINKEDLIST$GET(LL_CNTL,'RATE',1,RATE)
          CALL LINKEDLIST$GET(LL_CNTL,'FINAL',1,FINAL)
          RATE=-RATE
          FINAL=-FINAL
        ELSE
          CALL ERROR$MSG('DIAL NAME NOT RECOGNIZED')
          CALL ERROR$CHVAL('DIALID',DIALID)
          CALL ERROR$STOP('READIN_MERMIN')
        END IF
        CALL DIALS$SETL4('ON',.TRUE.)
        CALL DIALS$SETR8('RATE',RATE)
        CALL DIALS$SETR8('FINAL',FINAL)
        CALL DIALS$SELECT('..')
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
      CALL READIN_AUTO(LL_CNTL,'OCC',9.6D-2,1.D0,0.5D0,1.D0)

                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_CELL(LL_CNTL_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_    ! CONTROL LINKED LIST
      TYPE(LL_TYPE)            :: LL_CNTL
      LOGICAL(4)               :: TCHK,TCHK1
      LOGICAL(4)               :: TMOVE       ! SWITCH TO DYNAMIC LATTICE
      REAL(8)                  :: P           ! EXTERNAL PRESSURE
      REAL(8)                  :: MASS        ! MASS FOR CELL DFYNAMICS
      REAL(8)                  :: FRIC        ! MASS FOR CELL DFYNAMICS
      CHARACTER(32)            :: CH32VAL
      REAL(8)                  :: DT          ! TIME STEP
      REAL(8)                  :: CONSTR(3,3)
      INTEGER(4)               :: I,NCONSTR
      REAL(8)                  :: GIGA,PASCAL
!     ******************************************************************
      LL_CNTL=LL_CNTL_
!
!     ==========================================================================
!     ==  READ BLOCK GENERIC                                                  ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
      CALL LINKEDLIST$GET(LL_CNTL,'DT',1,DT)
      CALL WAVES$SETR8('TIMESTEP',DT)
      CALL CELL$SETR8('DT',DT)
!
!     ==========================================================================
!     ==  NOW ENTER CELL BLOCK                                                ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'CELL',1,TCHK)
      CALL CELL$SETL4('ON',TCHK)
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_CNTL,'CELL')
!
!     == DYNAMICAL OR NOT ======================================================
      TMOVE=.TRUE.
      CALL LINKEDLIST$EXISTD(LL_CNTL,'MOVE',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'MOVE',1,TMOVE)
      END IF
      CALL CELL$SETL4('MOVE',TMOVE)
!
!     == CONSTRAINTS ON THE CELL DYNAMICS ======================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'CONSTRAINT',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('UNKNOWN DATA KEYWORD "CONSTRAINT" IN !CNTL!CELL')
        CALL ERROR$MSG('PROBABLY YOU MEAN "CONSTRAINTTYPE"')
        CALL ERROR$MSG('OR THE BLOCK "CONSTRAINT"')
        CALL ERROR$MSG('PLEASE CHECK THE MANUAL')
        CALL ERROR$STOP('READIN_CELL')
      END IF

      CALL LINKEDLIST$EXISTD(LL_CNTL,'CONSTRAINTTYPE',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'CONSTRAINTTYPE',1,CH32VAL)
        CALL CELL$SETCH('CONSTRAINTTYPE',CH32VAL)
      END IF

      CALL LINKEDLIST$NLISTS(LL_CNTL,'CONSTRAINT',NCONSTR)
      IF(TCHK.AND.(NCONSTR.GT.0)) THEN
        CALL ERROR$MSG('CONSTRAINTTYPE AND !CONSTRAINT ARE INCOMPATIBLE')
        CALL ERROR$STOP('READIN_CELL')
      END IF
      DO I=1,NCONSTR
        CALL LINKEDLIST$SELECT(LL_CNTL,'CONSTRAINT',I)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'C',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!CELL!CONSTRAINT:C IS MANDATORY')
          CALL ERROR$STOP('READIN_CELL')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'C',1,CONSTR)
        CALL CELL$SETR8A('CONSTRAINT',9,CONSTR)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO

!     == FRICTION ======================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'FRIC',1,0.D0)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'FRIC',1,FRIC)
      CALL CELL$SETR8('FRICTION',FRIC)
!
!     == MASS IF DYNAMICAL =============================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'M',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'M',1,MASS)
      ELSE
!       == NEGATIVE MASS TRIGGERS DEFAULT CHOICE IN CELL OBJECT
        MASS=-1.D0
      END IF
      CALL CELL$SETR8('MASS',MASS)
!     
!     == PRESSURE ======================================================
      P=0.D0
      CALL LINKEDLIST$EXISTD(LL_CNTL,'P',1,TCHK)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'P[GPA]',1,TCHK1)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'P',1,P)
        IF(TCHK1) THEN
          CALL ERROR$MSG('OPTIONS P AND P[GPA] IN !CONTROL!CELL')
          CALL ERROR$MSG('ARE MUTUALLY EXCLUSIVE')
          CALL ERROR$STOP('READIN_CELL')
        END IF
      ELSE
        IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'P[GPA]',1,P)
          CALL CONSTANTS$GET('GIGA',GIGA)
          CALL CONSTANTS$GET('PASCAL',PASCAL)
          P=P*GIGA*PASCAL       
        END IF
      END IF

      CALL CELL$SETR8('P',P)
!     
!     == STRESS ========================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STRESS',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('THE OPTION !CONTROL!CELL:STRESS HAS BEEN DISABLED')
        CALL ERROR$STOP('READIN_CELL')
!       CALL LINKEDLIST$GET(LL_CNTL,'STRESS',1,STRESS)
!       CALL CELL$SETR8A('STRESS',9,STRESS)
      END IF
!
      RETURN
      END SUBROUTINE READIN_CELL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_SHADOW(LL_CNTL_)
!     ******************************************************************
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK
      REAL(8)               :: SVAR,SVAR1,SVAR2
      INTEGER(4)            :: ISVAR
      INTEGER(4)            :: NFILO
      LOGICAL(4)            :: TON
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_SHADOW')  
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'SHADOW',1,TON)
      IF(.NOT.TON) THEN
        CALL TRACE$POP
        RETURN
      END IF
      CALL ERROR$MSG('!CONTROL!SHADOW IS DISABLED')
      CALL ERROR$STOP('READIN_SHADOW')
      CALL LINKEDLIST$SELECT(LL_CNTL,'SHADOW')

!     == SWITCH LONG-RANGE INTERACTIONS ON/OFF  ========================
!     CALL LINKEDLIST$EXISTD(LL_CNTL,'LONGRANGE',1,TCHK)
!     IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'LONGRANGE',0,.FALSE.)
!     CALL LINKEDLIST$GET(LL_CNTL,'LONGRANGE',0,TCHK)
!     CALL SHADOW$SET('TLONGRANGE',4,.TRUE.)
!
!     ==================================================================
!     ==  READ BLOCK !CONTROL!OPTIMIZE                                ==
!     ==================================================================
      CALL LINKEDLIST$EXISTL(LL_CNTL,'OPTIMIZE',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'OPTIMIZE')
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'TOL',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'TOL',0,1.D-6)
        CALL LINKEDLIST$GET(LL_CNTL,'TOL',0,SVAR)
!       CALL SHADOW$SET('TOL',8,SVAR)
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NSTEP',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'NSTEP',0,10000)
        CALL LINKEDLIST$GET(LL_CNTL,'NSTEP',0,ISVAR)
!       CALL SHADOW$SET('NSTEP',4,ISVAR)
!
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      END IF
!
!     ==================================================================
!     ==  READ BLOCK !CONTROL!SHADOW!PRECONDITION                     ==
!     ==================================================================
      CALL LINKEDLIST$EXISTL(LL_CNTL,'PRECONDITION',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'PRECONDITION')
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'OMEGA1[CM**-1]',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'OMEGA1[CM**-1]',0,300.D0)
        CALL LINKEDLIST$GET(LL_CNTL,'OMEGA1[CM**-1]',1,SVAR1)
        SVAR1=SVAR1*4.556D-6
!       CALL SHADOW$SET('OMEGA1',8,SVAR1)
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'OMEGA2[CM**-1]',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'OMEGA2[CM**-1]',0,3000.D0)
        CALL LINKEDLIST$GET(LL_CNTL,'OMEGA2[CM**-1]',1,SVAR2)
        SVAR2=SVAR2*4.556D-6
!       CALL SHADOW$SET('OMEGA2',8,SVAR2)
!
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      END IF
                           CALL TRACE$POP
!
!     ==================================================================
!     ==  REPORT  -> SHOULD BE MOVED INTO SHADOW OBJECT               ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,FMT='()')
      WRITE(NFILO,FMT='("SHADOW"/"====")')     
!     CALL SHADOW$GET('OPTIMIZE',4,TCHK)
      IF(TCHK) THEN
        WRITE(NFILO,FMT='(55("."),": ",T1,"OPTIMIZE STRUCTURE")')
!       CALL SHADOW$GET('NSTEP',4,ISVAR)
        WRITE(NFILO,FMT='(55("."),": "' &
     &           //',T1,"MAX. NUMBER OF ITERATIONS",T58,F10.5)')ISVAR
!       CALL SHADOW$GET('TOL',8,SVAR)
        WRITE(NFILO,FMT='(55("."),": ",T1,"TOLERANCE ON THE FORCE"' &
     &           //',T58,F10.5)')SVAR
      END IF
!     CALL SHADOW$GET('PRECONDITION',4,TCHK)
      IF(TCHK) THEN
        WRITE(NFILO,FMT='(55("."),": ",T1,"PRECONDITION DYNAMICS")')
!       CALL SHADOW$GET('OMEGA1',8,SVAR)
        WRITE(NFILO,FMT='(55("."),": ",T1,"OMEGA1",T58,F10.5)')SVAR
!       CALL SHADOW$GET('OMEGA2',8,SVAR)
        WRITE(NFILO,FMT='(55("."),": ",T1,"OMEGA2",T58,F10.5)')SVAR
      END IF
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_QMMM(LL_CNTL_)
!     ******************************************************************
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK
      LOGICAL(4)            :: TON
      LOGICAL(4)            :: TSTOP
      LOGICAL(4)            :: TMOVIE
      REAL(8)               :: AMPRE
      REAL(8)               :: FRICTION
      REAL(8)               :: KELVIN
      REAL(8)               :: DT
      INTEGER(4)            :: MULTIPLE
      INTEGER(4)            :: LOD, NSKIP
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_QMMM')  
      LL_CNTL=LL_CNTL_
      CALL CONSTANTS('KB',KELVIN)
!
!     ==================================================================
!     ==  READ BLOCK GENERIC                                          ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
      CALL LINKEDLIST$GET(LL_CNTL,'DT',1,DT)
!
!     ==================================================================
!     ==  READ BLOCK !CONTROL!QM-MM                                   ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'QM-MM',1,TON)
      CALL QMMM$SETL4('ON',TON)
!      CALL THERMOSTAT$NEW('QM-MM')   ! CREATED BUT NOT SWITCHED ON, IF NOT USED
      IF(.NOT.TON) THEN
!       == THERMOSTAT REMAINS OFF, BUT MUST EXIST FOR THE REPORT
        CALL THERMOSTAT$NEW('QM-MM')
        CALL TRACE$POP
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'QM-MM')
!
!     == DEFAULT VALUES ================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FREEZE',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FREEZE',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'FREEZE',1,TCHK)
      CALL QMMM$SETL4('MOVE',.NOT.TCHK)
!
!     == OVERSAMPLE WITH NMULTIPLE TIME STEPS ==========================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'MULTIPLE',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'MULTIPLE',0,1)
      CALL LINKEDLIST$GET(LL_CNTL,'MULTIPLE',1,MULTIPLE)
      CALL QMMM$SETI4('MULTIPLE',MULTIPLE)
!  
!     == CREATE MOVIE FILE FROM OVERSAMPLING ===========================
      IF(TCHK.AND.(MULTIPLE.GE.2)) THEN
         CALL LINKEDLIST$EXISTD(LL_CNTL,'MOVIE',1,TMOVIE)
         IF(.NOT.TMOVIE)CALL LINKEDLIST$SET(LL_CNTL,'MOVIE',0,.FALSE.)
         CALL LINKEDLIST$GET(LL_CNTL,'MOVIE',1,TMOVIE)
         CALL QMMM$SETL4('MOVIE',TMOVIE)
         CALL LINKEDLIST$EXISTD(LL_CNTL,'SKIP',1,TCHK)
         IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'SKIP',0,0)
         CALL LINKEDLIST$GET(LL_CNTL,'SKIP',1,NSKIP)
         CALL QMMM$SETI4('SKIP',NSKIP)
      END IF
!
!     == START WITH ZERO INITIAL VELOCITIES ============================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STOP',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'STOP',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'STOP',1,TSTOP)
      CALL QMMM$SETL4('STOP',TSTOP)
!
!     == RANDOMIZE INITIAL VELOCITIES ==================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RANDOM[K]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'RANDOM[K]',0,0.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'RANDOM[K]',1,AMPRE) 
      AMPRE=AMPRE*KELVIN
      IF(AMPRE.NE.0) CALL QMMM$SETR8('RANDOM',AMPRE)
!
!     == FRICTION OFR THE MM ATOMS ======================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FRIC',0,0.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'FRIC',0,FRICTION)
      CALL QMMM$SETR8('FRICTION',FRICTION)
!
!     == FRICTION OFR THE MM ATOMS ======================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ADIABATIC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ADIABATIC',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'ADIABATIC',0,TCHK)
      CALL QMMM$SETL4('ADIABATIC',TCHK)
!     == LEVEL OF DETAIL FOR THE OUTPUT =================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'LOD',1,TCHK)
      LOD=10        !DEFAULT IS 10. IF LOD.GE.10 THEN YOU GET THE FULL OUTPUT
      IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'LOD',1,LOD)
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$SETI4('LOD',LOD)
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$SETI4('LOD',LOD)
!      
!     ==================================================================
!     ==  !CONTROL!QM-MM!AUTO                                         ==
!     ==================================================================
      CALL READIN_AUTO(LL_CNTL,'QM-MM',0.0D0,1.D0,0.5D0,1.D0)
!      
!     ==================================================================
!     ==  !CONTROL!QM-MM!THERMOSTAT                                   ==
!     ==================================================================
      CALL READIN_THERMOSTAT(LL_CNTL,'QM-MM',DT,0.5D0*293.D0*KELVIN,10.D0)
                           CALL TRACE$POP 
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_ANALYSE_HYPERFINE(LL_CNTL_)
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: NHYP
      INTEGER(4)            :: IHYP
      CHARACTER(32)         :: ATOM
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_ANALYSE_HYPERFINE')  
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'HYPERFINE',NHYP)
      DO IHYP=1,NHYP
        CALL LINKEDLIST$SELECT(LL_CNTL,'HYPERFINE',IHYP)
        CALL LINKEDLIST$GET(LL_CNTL,'ATOM',0,ATOM)
        CALL HYPERFINE$SELECT('~')
        CALL HYPERFINE$SELECT(ATOM)
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EFG',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EFG',0,.FALSE.)
        CALL LINKEDLIST$GET(LL_CNTL,'EFG',1,TCHK)
        CALL HYPERFINE$SETL4('TEFG',TCHK)
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ISOMERSHIFT',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ISOMERSHIFT',0,.FALSE.)
        CALL LINKEDLIST$GET(LL_CNTL,'ISOMERSHIFT',1,TCHK)
        CALL HYPERFINE$SETL4('TIS',TCHK)
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FERMICONTACT',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FERMICONTACT',0,.FALSE.)
        CALL LINKEDLIST$GET(LL_CNTL,'FERMICONTACT',1,TCHK)
        CALL HYPERFINE$SETL4('TFC',TCHK)
!
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ANISOTROPIC',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ANISOTROPIC',0,.FALSE.)
        CALL LINKEDLIST$GET(LL_CNTL,'ANISOTROPIC',1,TCHK)
        CALL HYPERFINE$SETL4('TANIS',TCHK)
!
        CALL HYPERFINE$SELECT('~')
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_ANALYSE_CORELEVEL(LL_CNTL_)
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN)  :: LL_CNTL_
      TYPE(LL_TYPE)             :: LL_CNTL
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: N
      CHARACTER(32),ALLOCATABLE :: ATOMS(:)
      LOGICAL(4)                :: DEFAULT
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_ANALYSE_CORELEVEL')  
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'CORELEVELS',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL CORE$SETL4('DEFAULT',.FALSE.)
        CALL TRACE$POP
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'CORELEVELS')
      CALL LINKEDLIST$EXISTD(LL_CNTL,'DEFAULT',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'DEFAULT',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'DEFAULT',1,DEFAULT)
      CALL CORE$SETL4('DEFAULT',DEFAULT)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOMS',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SIZE(LL_CNTL,'ATOMS',1,N)
        ALLOCATE(ATOMS(N))
        CALL LINKEDLIST$GET(LL_CNTL,'ATOMS',1,ATOMS)
        CALL CORE$SETCHA('ATOMS',N,ATOMS)
        DEALLOCATE(ATOMS)
      END IF
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_COSMO(LL_CNTL_)
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      LOGICAL(4)               :: TCHK,TCHK1,TCHK2
      REAL(8)                  :: SVAR
      REAL(8)                  :: FRIC
      INTEGER(4)               :: MULTIPLE
      REAL(8)                  :: DT
!     ******************************************************************
      LL_CNTL=LL_CNTL_
!
!     ==================================================================
!     ==  READ BLOCK !CONTROL!GENERIC                                 ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
      CALL LINKEDLIST$GET(LL_CNTL,'DT',1,DT)
!     ==================================================================
!     ==  READ BLOCK !CONTROL!COSMO                                   ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'COSMO',1,TCHK)
      CALL COSMO$SETL4('ON',TCHK)
      IF(.NOT.TCHK) RETURN
                               CALL TRACE$PUSH('READIN_COSMO')
      CALL LINKEDLIST$SELECT(LL_CNTL,'COSMO')
!
!     == STOP =======================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'STOP',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'STOP',0,.FALSE.)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'STOP',1,TCHK)
      CALL COSMO$SETL4('STOP',TCHK)
!
!     == START =====================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'START',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'START',0,.FALSE.)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'START',1,TCHK)
      CALL COSMO$SETL4('START',TCHK)
!
!     == MASS ===================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'M',1,TCHK)
      IF(.NOT.TCHK) THEN
        PRINT*,'WARNING! NO MASS FOR SURFACE CHARGES SUPPLIED!'
        PRINT*,'DEFAULT VALUE OF 1000.D0 IS SUPPLIED BY PAW!'
        CALL LINKEDLIST$SET(LL_CNTL,'M',1,1000.D0)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'M',1,SVAR)
      CALL COSMO$SETR8('MASS',SVAR)
!
!     == MULTIPLE TIMESTEPS?  ========================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'MULTIPLE',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'MULTIPLE',1,MULTIPLE)
        CALL COSMO$SETI4('MULTIPLE',MULTIPLE)
      END IF
!
!     == ADIABATIC OR DYNAMIC?========================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ADIABATIC',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'ADIABATIC',0,.FALSE.)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'ADIABATIC',1,TCHK)
      CALL COSMO$SETL4('ADIABATIC',TCHK)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ETOL',1,TCHK1)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'QTOL',1,TCHK2)
      IF(.NOT.(TCHK1.OR.TCHK2)) THEN
        CALL COSMO$SETR8('ETOL',1.D-7)
      END IF
      IF(TCHK1) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'ETOL',1,SVAR)
        CALL COSMO$SETR8('ETOL',SVAR)
      END IF
!
      IF(TCHK2) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'QTOL',1,SVAR)
        CALL COSMO$SETR8('QTOL',SVAR)
      END IF
!
!     == EPSILON = DIELECTRIC CONSTANT ===========================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EPSILON',1,TCHK)
      IF(TCHK) THEN 
        CALL ERROR$MSG('EPSILON FOR COSMO SCREENING')
        CALL ERROR$MSG('MUST BE SPECIFIED IN THE STRUCTURE FILE')
        CALL ERROR$STOP('READIN_COSMO')
      END IF
!
!     ==  FRIC ===================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'FRIC',0,0.D0)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'FRIC',1,FRIC)
      CALL COSMO$SETR8('FRICTION',FRIC)
!
!     ==  VPAULI ===================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'VPAULI',1,TCHK)
      IF(TCHK) THEN
        CALL ERROR$MSG('VPAULI FOR COSMO SCREENING')
        CALL ERROR$MSG('MUST BE SPECIFIED IN THE STRUCTURE FILE')
        CALL ERROR$STOP('READIN_COSMO')
      END IF
!
!     ==  PRINT CHARGES ================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'CHARGES',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'CHARGES',0,.FALSE.)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'CHARGES',1,TCHK)
      CALL COSMO$SETL4('CHARGES',TCHK)
!
!     ==================================================================
!     ==  READ BLOCK !CONTROL!COSMO!OPTFRIC                           ==
!     ==================================================================
      CALL OPTFRIC$NEW('COSMO')
      CALL OPTFRIC$SELECT('COSMO')
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'OPTFRIC',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'OPTFRIC',0,.FALSE.)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'OPTFRIC',1,TCHK)
      CALL OPTFRIC$SETL4('ON',TCHK)
      CALL OPTFRIC$SETR8('STARTFRIC',FRIC)
!
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RETARD',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'RETARD',0,0.D0)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'RETARD',1,SVAR)
      CALL OPTFRIC$SETR8('RETARD',SVAR)
!
!     ==================================================================
!     ==  READ BLOCK !CONTROL!COSMO!AUTO                              ==
!     ==================================================================
!      CALL READIN_AUTO(LL_CNTL,'COSMO',0.3D0,0.96D0,1.D0,1.D0)
!
!     ==================================================================
!     ==  READ BLOCK !QNOSE                                           ==
!     ==================================================================
!      CALL CONSTANTS('KB',KELVIN)
!      CALL READIN_THERMOSTAT(LL_CNTL_,'COSMO',DT,0.5D0*293.D0*KELVIN,10.D0)
!
                                 CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_ANALYSE_TRAJECTORIES(LL_CNTL_)
!     ******************************************************************
!     **  ACTIVATES TRAJECTORY MODULE TO WRITE TRAJECTORIES           **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: NSKIP
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_ANALYSE_TRAJECTORIES')  
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TRA')

!
!     ==================================================================
!     ==  NUMBER OF SLICES TO BE SKIPPED FROM ALL TRAKJECTORIES       ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NSKIP',0,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'NSKIP',0,0)
      CALL LINKEDLIST$GET(LL_CNTL,'NSKIP',1,NSKIP)
!
!     ==================================================================
!     == POSITION TRAJECTORY                                          ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'R',0,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'R',0,.TRUE.)
      CALL LINKEDLIST$GET(LL_CNTL,'R',1,TCHK)
      CALL TRAJECTORYIO$SELECT('POSITION-TRAJECTORY')
      CALL TRAJECTORYIO$SETL4('ON',TCHK)
      CALL TRAJECTORYIO$SETI4('SKIP',NSKIP)
      CALL TRAJECTORYIO$SELECT('NONE')
!
!     ==================================================================
!     == QM-MM POSITION TRAJECTORY                                    ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'QMMM',0,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'QMMM',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'QMMM',1,TCHK)
      CALL TRAJECTORYIO$SELECT('QMMM-POS-TRA')
      CALL TRAJECTORYIO$SETL4('ON',TCHK)
      CALL TRAJECTORYIO$SETI4('SKIP',NSKIP)
      CALL TRAJECTORYIO$SELECT('NONE')
!
!     ==================================================================
!     ==  FORCE TRAJECTORY                                            ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FORCE',0,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FORCE',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'FORCE',1,TCHK)
      CALL TRAJECTORYIO$SELECT('FORCE-TRAJECTORY')
      CALL TRAJECTORYIO$SETL4('ON',TCHK)
      CALL TRAJECTORYIO$SETI4('SKIP',NSKIP)
      CALL TRAJECTORYIO$SELECT('NONE')
!
!     ==================================================================
!     ==  ENERGY TRAJECTORY                                           ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'E',0,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'E',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'E',1,TCHK)
      CALL TRAJECTORYIO$SELECT('ENERGY-TRAJECTORY')
      CALL TRAJECTORYIO$SETL4('ON',TCHK)
      CALL TRAJECTORYIO$SETI4('SKIP',NSKIP)
      CALL TRAJECTORYIO$SELECT('NONE')
!
!     ==================================================================
!     ==  ENERGY BANDS                                                ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'BANDS',0,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'BANDS',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'BANDS',1,TCHK)
      CALL TRAJECTORYIO$SELECT('BANDS-TRAJECTORY')
      CALL TRAJECTORYIO$SETL4('ON',TCHK)
      CALL TRAJECTORYIO$SETI4('SKIP',NSKIP)
      CALL TRAJECTORYIO$SELECT('NONE')
!
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_ANALYSE_WAVE(LL_CNTL_)
!     ******************************************************************
!     ==                                                              ==
!     ==  FILE=FILENAME                                               ==
!     ==  TITLE= IMAGE TITLE                                          ==
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: NWAVE
      INTEGER(4)            :: IWAVE
      INTEGER(4)            :: IB,IKPT,ISPIN
      REAL(8)               :: DR
      CHARACTER(256)        :: CH256SVAR1
      LOGICAL(4)            :: TIMAG
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_ANALYSE_WAVE')  
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'WAVE',NWAVE)
      CALL GRAPHICS$SETI4('NWAVE',NWAVE)
      DO IWAVE=1,NWAVE
        CALL LINKEDLIST$SELECT(LL_CNTL,'WAVE',IWAVE)
        CALL GRAPHICS$SETI4('IWAVE',IWAVE)
!
!       == SET TITLE OF THE IMAGE ======================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'TITLE',1,TCHK)
        IF(.NOT.TCHK) THEN
          CH256SVAR1=' '
          WRITE(CH256SVAR1,*)'WAVE',IWAVE
          CALL LINKEDLIST$SET(LL_CNTL,'TITLE',0,TRIM(CH256SVAR1))
        ELSE
          CALL LINKEDLIST$GET(LL_CNTL,'TITLE',1,CH256SVAR1)
        END IF
        CALL GRAPHICS$SETCH('TITLE',TRIM(CH256SVAR1))
!
!       == SET FILENAME FOR THE IMAGE===================================
        WRITE(CH256SVAR1,*)IWAVE
        CH256SVAR1=ADJUSTL(CH256SVAR1)
        CH256SVAR1='WAVE'//TRIM(CH256SVAR1(1:249))//'.WV'
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,CH256SVAR1)
        CALL GRAPHICS$SETCH('FILE',TRIM(CH256SVAR1))

!       == GRID SPACING
        DR=0.4D0
        CALL LINKEDLIST$EXISTD(LL_CNTL,'DR',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'DR',1,DR)
        CALL GRAPHICS$SETR8('DR',DR)
!
!       == BAND INDEX
        CALL LINKEDLIST$EXISTD(LL_CNTL,'B',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('BAND INDEX IS MANDATORY')
          CALL ERROR$CHVAL('FILE',CH256SVAR1)
          CALL ERROR$STOP('READIN_ANALAYSE_WAVE')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'B',1,IB)
        CALL GRAPHICS$SETI4('IB',IB)
!
!       == K-POINT INDEX
        IKPT=1
        CALL LINKEDLIST$EXISTD(LL_CNTL,'K',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'K',1,IKPT)
        CALL GRAPHICS$SETI4('IKPT',IKPT)
!
!       == SPIN INDEX
        ISPIN=1
        CALL LINKEDLIST$EXISTD(LL_CNTL,'S',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'S',1,ISPIN)
        CALL GRAPHICS$SETI4('ISPIN',ISPIN)
!
!       == SWITCH FOR IMAGINARY PART
        TIMAG=.FALSE.
        CALL LINKEDLIST$EXISTD(LL_CNTL,'IMAG',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'IMAG',1,TIMAG)
        CALL GRAPHICS$SETL4('IMAG',TIMAG)
!
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_ANALYSE_OPTEELS(LL_CNTL_)
!     **************************************************************************
!     **  INPUT DATA FOR OPTICAL ABSORPTION AND EELS SPECTRA                  **
!     **  CURRENTLY USED TO SET THE OPTEELS OBJECT                            **
!     **                                                                      **
!     **  FILE=FILENAME                                                       **
!     **                 MATTHE A. UIJTTEWAAL 2011 (ADJUSTED P.BLOECHL 2012)  **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_ !FROM IO_MODULE
      TYPE(LL_TYPE)            :: LL_CNTL
      LOGICAL(4)               :: TCHK
      LOGICAL(4)               :: TCH
      LOGICAL(4)               :: ORIG
      INTEGER(4)               :: NOPT
      INTEGER(4)               :: IOPT
      INTEGER(4)               :: NEELS
      INTEGER(4)               :: IEELS
      INTEGER(4)               :: IC
      INTEGER(4)               :: NR3
      INTEGER(4)               :: NI,NF !BI,BF
      CHARACTER(256)           :: CH256
      CHARACTER(32)            :: CH32
      CHARACTER(8)             :: CH8
      REAL(8)                  :: SVAR
      REAL(8)                  :: ZI,ZF
      REAL(8)                  :: EMAX
      REAL(8)                  :: EV     ! ELECTRON VOLT
      REAL(8)                  :: ANGSTROM  ! ANGSTROM
!     **************************************************************************
                           CALL TRACE$PUSH('READIN_ANALYSE_OPT')  
      CALL CONSTANTS$GET('EV',EV)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'EELS',NEELS)
      CALL OPTEELS$SETI4('NEELS',NEELS)
      DO IEELS=1,NEELS
        CALL LINKEDLIST$SELECT(LL_CNTL,'EELS',IEELS)
        CALL OPTEELS$SETI4('IEELS',IEELS)
!
!       == ATOM NAME FOR EELS ==================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOM',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!CONTROL!ANALYSE!EELS:ATOM IS A MANDATORY INPUT')
          CALL ERROR$STOP('READIN_ANALYSE_OPTEELS')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'ATOM',1,CH32)
        CALL OPTEELS$SETCH('ATOM',CH32)
!
!       == SHELL INDEX FOR EELS ================================================
        IC=1 !DEFAULT 1S
        CALL LINKEDLIST$EXISTD(LL_CNTL,'IC',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'IC',1,IC) 
        CALL OPTEELS$SETI4('IC',IC)
!
!       == SHELL INSTRUMENTAL BROADENING =======================================
        SVAR=1.D0
        CALL LINKEDLIST$EXISTD(LL_CNTL,'BROADENING[EV]',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'BROADENING[EV]',1,SVAR) 
        SVAR=SVAR*EV
        IF(SVAR.LE.0.D0) THEN
          CALL ERROR$MSG('ILLEGAL VALU OF !CONTROL!ANALYSE!EELS:BROADENING[EV]')
          CALL ERROR$MSG('BROADENING[EV] MUST BE LARGER THAN ZERO')
          CALL ERROR$R8VAL('BROADENING[EV]',SVAR)
          CALL ERROR$STOP('READIN_ANALYSE_OPTEELS')
        END IF
        CALL OPTEELS$SETR8('BROADENING',SVAR)
!
!       == SET FILENAME FOR THE DATA============================================
        WRITE(CH256,*)IC
        CH256=TRIM(ADJUSTL(CH32))//'_IC'//TRIM(ADJUSTL(CH256))//'.EELS'
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,CH256)
        CALL OPTEELS$SETCH('FILE',TRIM(CH256))
!
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
!
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
      CALL LINKEDLIST$NLISTS(LL_CNTL,'OPT',NOPT)
      CALL OPTEELS$SETI4('NOPT',NOPT)
      DO IOPT=1,NOPT
        CALL LINKEDLIST$SELECT(LL_CNTL,'OPT',IOPT)
        CALL OPTEELS$SETI4('IOPT',IOPT)
!
!       == SET FILENAME FOR THE DATA============================================
!       WRITE(CH256,*)IOPT ???
!       CH256=ADJUSTL(CH256)
        CH256='OPT'//ACHAR(IOPT)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,CH256)
        CALL OPTEELS$SETCH('FILE',TRIM(CH256))
!
!       == GET MAX.EXC.EN. =====================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EMAX[EV]',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('MAX. EXCITATION ENERGY (EMAX[EV]) IS MANDATORY')
          CALL ERROR$STOP('READIN_ANALAYSE_OPT')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'EMAX[EV]',1,EMAX)
        EMAX=EMAX*EV
        IF(EMAX.LT.0.) THEN
          CALL ERROR$MSG('EXCITATION ENERGY SHOULD BE POSITIVE')
          CALL ERROR$STOP('READIN_ANALAYSE_OPT')
        END IF
        CALL OPTEELS$SETR8('EMAX',EMAX)

!       == ORIGIN CENTERED STATES? =============================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'O',1,TCHK)
        TCH=.TRUE.
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'O',1,TCH)
        CALL OPTEELS$SETL4('ORIG',TCH)

!       == QUADRUPOLE CONTRIBUTION? ============================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'Q',1,TCHK)
        TCH=.FALSE.
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'Q',1,TCH)
        CALL OPTEELS$SETL4('QUAD',TCH)

!       == RESTRICT OPT TO SPECIFIED Z-RANGE ===================================
        ZI=-100. !DEFAULT NO RESTRICTION
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ZI[AA]',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'ZI[AA]',1,ZI) 
        ZI=ZI*ANGSTROM !INPUT IN ANGSTROM!!
        CALL OPTEELS$SETR8('ZI',ZI)
!
        ZF=100. !DEFAULT NO RESTRICTION
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ZF[A]',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'ZF[A]',1,ZF)
        ZF=ZF*ANGSTROM
        CALL OPTEELS$SETR8('ZF',ZF)
!
!       == RESTRICT OPT TO SPECIFIED GROUP =====================================
        CH256='ALL' !DEFAULT NO RESTRICTION
        CALL LINKEDLIST$EXISTD(LL_CNTL,'GROUP',1,TCHK) !RESTRICT OPT TO GROUP!
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'GROUP',1,CH256) 
        CALL OPTEELS$SETCH('GROUP',CH256)
!
!       == ATOM NAME FOR EELS ==================================================
        CH32=' ' !DEFAULT NOT SELECTED
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOM',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'ATOM',1,CH32)
        CALL OPTEELS$SETCH('ATOM',CH32)

!       == SHELL INDEX FOR EELS ================================================
        IC=1 !DEFAULT 1S
        CALL LINKEDLIST$EXISTD(LL_CNTL,'IC',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'IC',1,IC) !CHAR OR INT??
        CALL OPTEELS$SETI4('IC',IC)
!
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_ANALYSE_DENSITY(LL_CNTL_)
!     ******************************************************************
!     ==                                                              ==
!     ==  FILE=FILENAME                                               ==
!     ==  TITLE= IMAGE TITLE                                          ==

!     ******************************************************************
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK,TCHK1,TCHK2
      INTEGER(4)            :: NDEN
      INTEGER(4)            :: IDEN
      REAL(8)               :: DR
      CHARACTER(256)        :: CH256SVAR1
      CHARACTER(8)          :: TYPE
      LOGICAL(4)            :: TDIAG,TOCC,TCORE
      REAL(8)               :: EMIN,EMAX
      INTEGER(4)            :: IBMIN,IBMAX
      REAL(8)               :: EV
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_ANALYSE_DENSITY')  
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'DENSITY',NDEN)
      CALL GRAPHICS$SETI4('NDENSITY',NDEN)
      DO IDEN=1,NDEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'DENSITY',IDEN)
        CALL GRAPHICS$SETI4('IDENSITY',IDEN)
!
!       == SET TITLE OF THE IMAGE ======================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'TITLE',1,TCHK)
        IF(.NOT.TCHK) THEN
          CH256SVAR1=' '
          WRITE(CH256SVAR1,*)'DENSITY',IDEN
          CALL LINKEDLIST$SET(LL_CNTL,'TITLE',0,TRIM(CH256SVAR1))
        ELSE
          CALL LINKEDLIST$GET(LL_CNTL,'TITLE',1,CH256SVAR1)
        END IF
        CALL GRAPHICS$SETCH('TITLE',TRIM(CH256SVAR1))
!
!       == SET FILENAME FOR THE IMAGE===================================
        WRITE(CH256SVAR1,*)IDEN
        CH256SVAR1=ADJUSTL(CH256SVAR1)
        CH256SVAR1='DENSITY'//TRIM(CH256SVAR1(1:246))//'.WV'
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,CH256SVAR1)
        CALL GRAPHICS$SETCH('FILE',TRIM(CH256SVAR1))
!
!       == GRID SPACING
        DR=0.4D0
        CALL LINKEDLIST$EXISTD(LL_CNTL,'DR',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'DR',1,DR)
        CALL GRAPHICS$SETR8('DR',DR)
!
!       == TYPE CAN BE 'TOTAL','SPIN','UP','DOWN'
        TYPE='TOTAL'
        CALL LINKEDLIST$EXISTD(LL_CNTL,'TYPE',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,TYPE)
        TYPE=+TYPE
        IF(TYPE.NE.'TOTAL'.AND.TYPE.NE.'SPIN' &
       &                  .AND.TYPE.NE.'UP'.AND.TYPE.NE.'DOWN') THEN
          CALL ERROR$MSG('INVALID VALUE OF VARIABLE TYPE')
          CALL ERROR$MSG('TYPE MUST BE EITHER "TOTAL" OR "SPIN" OR')
          CALL ERROR$MSG('                            OR "UP" OR "DOWN"')
          CALL ERROR$CHVAL('TYPE',TYPE)
          CALL ERROR$STOP('READIN_ANALYSE_DENSITY')
        END IF
        CALL GRAPHICS$SETCH('TYPE',TYPE)
!
!       == TOCC
        TOCC=.TRUE.
        CALL LINKEDLIST$EXISTD(LL_CNTL,'OCC',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'OCC',1,TOCC)
        CALL GRAPHICS$SETL4('TOCC',TOCC)
!
!       == TDIAG (USE EIGENSTATES)
        TDIAG=.TRUE.
        CALL LINKEDLIST$EXISTD(LL_CNTL,'DIAG',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'DIAG',1,TDIAG)
        CALL GRAPHICS$SETL4('TDIAG',TDIAG)
!
!       == TCORE 
        TCORE=.FALSE.
        CALL LINKEDLIST$EXISTD(LL_CNTL,'CORE',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'CORE',1,TCORE)
        CALL GRAPHICS$SETL4('TCORE',TCORE)
!
!       == SELECTOR (NONE, EMINMAX ,BMINMAX)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EMIN[EV]',1,TCHK1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EMAX[EV]',1,TCHK)
        TCHK1=TCHK.OR.TCHK1
        CALL LINKEDLIST$EXISTD(LL_CNTL,'BMIN',1,TCHK2)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'BMAX',1,TCHK)
        TCHK2=TCHK.OR.TCHK2
        IF(TCHK1.AND.TCHK2) THEN
          CALL ERROR$MSG('SELECT EITHER BANDS OR ENERGY RANGE')
          CALL ERROR$STOP('READING_ANALYSE_DENSITY')
        END IF
        IF(TCHK1) THEN
          CALL GRAPHICS$SETCH('SELECT','EMINMAX')
          EMIN=-1.D+10
          EMAX=+1.D+10
          CALL LINKEDLIST$GET(LL_CNTL,'EMIN[EV]',1,EMIN)
          CALL LINKEDLIST$GET(LL_CNTL,'EMAX[EV]',1,EMAX)
          CALL CONSTANTS$GET('EV',EV)
          EMIN=EMIN*EV
          EMAX=EMAX*EV
          CALL GRAPHICS$SETR8('EMIN',EMIN)
          CALL GRAPHICS$SETR8('EMAX',EMAX)
        ELSE IF(TCHK2) THEN
          CALL GRAPHICS$SETCH('SELECT','BMINMAX')
          IBMIN=1
          IBMAX=10000000
          CALL LINKEDLIST$GET(LL_CNTL,'BMIN',1,IBMIN)
          CALL LINKEDLIST$GET(LL_CNTL,'BMAX',1,IBMAX)
          CALL GRAPHICS$SETI4('IBMIN',IBMIN)
          CALL GRAPHICS$SETI4('IBMAX',IBMAX)
        ELSE
          CALL GRAPHICS$SETCH('SELECT','NONE')
        END IF
!
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_ANALYSE_POTENTIAL(LL_CNTL_)
!     ******************************************************************
!     ==                                                              ==
!     ==  FILE=FILENAME                                               ==
!     ==  TITLE= IMAGE TITLE                                          ==
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK
      REAL(8)               :: DR
      CHARACTER(256)        :: CH256SVAR1
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_ANALYSE_POTENTIAL')  
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'POTENTIAL',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL TRACE$POP
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'POTENTIAL')
      CALL GRAPHICS$SETL4('TPOT',.TRUE.)
!
!     == SET TITLE OF THE IMAGE ======================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'TITLE',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'TITLE',0,'POTENTIAL')
      CALL LINKEDLIST$GET(LL_CNTL,'TITLE',1,CH256SVAR1)
      CALL GRAPHICS$SETCH('TITLE',TRIM(CH256SVAR1))
!
!     == SET FILENAME FOR TEH IMAGE===================================
      CH256SVAR1='POTENTIAL.WV'
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,CH256SVAR1)
      CALL GRAPHICS$SETCH('FILE',TRIM(CH256SVAR1))
!
!     == GRID SPACING
      DR=0.4D0
      CALL LINKEDLIST$EXISTD(LL_CNTL,'DR',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'DR',1,DR)
      CALL GRAPHICS$SETR8('DR',DR)
      CALL LINKEDLIST$SELECT(LL_CNTL,'..')
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_ANALYSE_1DPOT(LL_CNTL_)
!     **************************************************************************
!     **  ONE-DIMENSIONAL HARTREE POTENTIAL ALONG A SPECIFIED LATTICE DIRECTION*
!     **                                                                      **
!     **  FILE=FILENAME                                                       **
!     **  TITLE= IMAGE TITLE                                                  **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: IT(3)
      INTEGER(4)            :: I1DPOT
      INTEGER(4)            :: N1DPOT
      CHARACTER(256)        :: TITLE
      CHARACTER(256)        :: FILE
!     **************************************************************************
                           CALL TRACE$PUSH('READIN_ANALYSE_1DPOT')  
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'1DPOT',N1DPOT)
      CALL GRAPHICS$SETI4('N1DPOT',N1DPOT)
      DO I1DPOT=1,N1DPOT
        CALL LINKEDLIST$SELECT(LL_CNTL,'1DPOT',I1DPOT)
        CALL GRAPHICS$SETI4('I1DPOT',I1DPOT)
!
!       == SET TITLE OF THE IMAGE ==============================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'TITLE',1,TCHK)
        IF(.NOT.TCHK) THEN
          WRITE(TITLE,*)I1DPOT
          TITLE=-'ONEDPOT_'//TRIM(ADJUSTL(TITLE))
          CALL LINKEDLIST$SET(LL_CNTL,'TITLE',0,TRIM(TITLE))
        ELSE
          CALL LINKEDLIST$GET(LL_CNTL,'TITLE',1,TITLE)
        END IF
        CALL GRAPHICS$SETCH('TITLE',TRIM(TITLE))
!
!       == SET FILENAME FOR THE IMAGE===========================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FILE',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,FILE)
        ELSE
          WRITE(FILE,*)I1DPOT
          FILE=-'ONEDPOT_'//TRIM(ADJUSTL(TITLE))//-'.DAT'
          CALL LINKEDLIST$SET(LL_CNTL,'TITLE',0,TRIM(FILE))
        END IF
        CALL GRAPHICS$SETCH('FILE',TRIM(FILE))
!
!       == AXIS FOR 1D POTENTIAL AXIS=RBAS*IT ==================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'IT',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('MISSING MANDATORY ARGUMENT!CONTROL!ANALYSE!1DPLOT:IT')
          CALL ERROR$I4VAL('I1DPOT',I1DPOT)
          CALL ERROR$STOP('READIN_ANALYSE_1DPOT')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'IT',0,IT)
        CALL GRAPHICS$SETI4A('IT',3,IT)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!
!!$!     ..................................................................
!!$      SUBROUTINE READIN_ANALYSE_WAVE_OLD(LL_CNTL_)
!!$!     ******************************************************************
!!$!     ==                                                              ==
!!$!     ==  FILE=FILENAME                                               ==
!!$!     ==  TITLE= IMAGE TITLE                                          ==
!!$!     ==  TYPE='WAVE','DENSITY'                                       ==
!!$!     ==  WEIGHTING='NONE','SPIN','TOTAL'                             ==
!!$!     ==  FORMAT=ALL,RANGE,LIST                                       ==
!!$!     ==  STATE->IB,IKPT,ISPIN                                        ==
!!$!     ******************************************************************
!!$      USE LINKEDLIST_MODULE
!!$      IMPLICIT NONE
!!$      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
!!$      TYPE(LL_TYPE)         :: LL_CNTL
!!$      LOGICAL(4)            :: TCHK
!!$      TYPE(LL_TYPE)         :: LL_IMAGELIST
!!$      INTEGER(4)            :: NWAVE
!!$      INTEGER(4)            :: NSTATE
!!$      INTEGER(4)            :: IWAVE,ISTATE
!!$      INTEGER(4)            :: IB,IKPT,ISPIN
!!$      CHARACTER(256)        :: CH256SVAR1
!!$      CHARACTER(8)          :: CH8SVAR1
!!$!     ******************************************************************
!!$                           CALL TRACE$PUSH('READIN_ANALYSE_WAVE')  
!!$      LL_CNTL=LL_CNTL_
!!$      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
!!$      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
!!$      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
!!$      CALL LINKEDLIST$NLISTS(LL_CNTL,'WAVE',NWAVE)
!!$      DO IWAVE=1,NWAVE
!!$        CALL LINKEDLIST$SELECT(LL_CNTL,'WAVE',IWAVE)
!!$             CALL GRAPHICS$GETLIST(LL_IMAGELIST)    
!!$             CALL LINKEDLIST$SELECT(LL_IMAGELIST,'~')
!!$             CALL LINKEDLIST$SELECT(LL_IMAGELIST,'IMAGE',-1)
!!$!
!!$!       == SET TITLE OF THE IMAGE ======================================
!!$        CALL LINKEDLIST$EXISTD(LL_CNTL,'TITLE',1,TCHK)
!!$        IF(.NOT.TCHK) THEN
!!$          CH256SVAR1=' '
!!$          WRITE(CH256SVAR1,*)'WAVE',IWAVE
!!$          CALL LINKEDLIST$SET(LL_CNTL,'TITLE',0,TRIM(CH256SVAR1))
!!$        ELSE
!!$          CALL LINKEDLIST$GET(LL_CNTL,'TITLE',1,CH256SVAR1)
!!$        END IF
!!$        CALL LINKEDLIST$SET(LL_IMAGELIST,'TITLE',0,TRIM(CH256SVAR1))
!!$!
!!$!       == SET FILENAME FOR TEH IMAGE===================================
!!$        CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,CH256SVAR1)
!!$        CALL LINKEDLIST$SET(LL_IMAGELIST,'FILE',0,TRIM(CH256SVAR1))
!!$!
!!$!       __ GET TYPE_(DENSITY,WAVE)____________________________________
!!$        CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,CH8SVAR1)
!!$        CALL LINKEDLIST$SET(LL_IMAGELIST,'TYPE',0,TRIM(CH8SVAR1))
!!$!
!!$!       __ GET WEIGHTING_(NONE,SPIN,TOTAL)____________________________
!!$        IF(CH8SVAR1(1:4).EQ.'TYPE') THEN
!!$          CALL LINKEDLIST$SET(LL_CNTL,'WEIGHTING',0,'NONE')
!!$        END IF
!!$        CALL LINKEDLIST$EXISTD(LL_CNTL,'WEIGHTING',1,TCHK)
!!$        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'WEIGHTING',0,'TOTAL')
!!$        CALL LINKEDLIST$GET(LL_CNTL,'WEIGHTING',1,CH8SVAR1)
!!$        CALL LINKEDLIST$SET(LL_IMAGELIST,'WEIGHTING',0,TRIM(CH8SVAR1))
!!$!
!!$!         __ GET FORMAT_(ALL,RANGE,LIST)________________________________
!!$        CALL LINKEDLIST$EXISTD(LL_CNTL,'FORMAT',1,TCHK)
!!$        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'FORMAT',0,'ALL')
!!$        CALL LINKEDLIST$GET(LL_CNTL,'FORMAT',1,CH8SVAR1)
!!$        CALL LINKEDLIST$SET(LL_IMAGELIST,'FORMAT',0,TRIM(CH8SVAR1))
!!$!
!!$!       ================================================================
!!$!       ==  READ BLOCK !ANALYSE!WAVE!STATE                            ==
!!$!       ================================================================
!!$        CALL LINKEDLIST$NLISTS(LL_CNTL,'STATE',NSTATE)
!!$        DO ISTATE=1,NSTATE
!!$          CALL LINKEDLIST$SELECT(LL_CNTL,'STATE',ISTATE)
!!$          CALL LINKEDLIST$SELECT(LL_IMAGELIST,'STATE',-1)
!!$!         __ BAND INDEX________________________________________________
!!$          CALL LINKEDLIST$EXISTD(LL_CNTL,'B',1,TCHK)
!!$          IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'B',0,1)
!!$          CALL LINKEDLIST$GET(LL_CNTL,'B',1,IB)
!!$          CALL LINKEDLIST$SET(LL_IMAGELIST,'IB',0,IB)
!!$!         __ K-POINT INDEX_____________________________________________
!!$          CALL LINKEDLIST$EXISTD(LL_CNTL,'K',1,TCHK)
!!$          IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'K',0,1)
!!$          CALL LINKEDLIST$GET(LL_CNTL,'K',1,IKPT)
!!$          CALL LINKEDLIST$SET(LL_IMAGELIST,'IKPT',0,IKPT)
!!$!         __ SPIN INDEX _______________________________________________
!!$          CALL LINKEDLIST$EXISTD(LL_CNTL,'S',1,TCHK)
!!$          IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'S',0,1)
!!$          CALL LINKEDLIST$GET(LL_CNTL,'S',1,ISPIN)
!!$          CALL LINKEDLIST$SET(LL_IMAGELIST,'ISPIN',0,ISPIN)
!!$!         __ SELECT IMAGINARY OR REAL PART_____________________________
!!$          CALL LINKEDLIST$EXISTD(LL_CNTL,'IMAG',1,TCHK)
!!$          IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'IMAG',0,.FALSE.)
!!$          CALL LINKEDLIST$GET(LL_CNTL,'IMAG',1,TCHK)
!!$          CALL LINKEDLIST$SET(LL_IMAGELIST,'TIM',0,TCHK)
!!$!
!!$          CALL LINKEDLIST$SELECT(LL_IMAGELIST,'..')
!!$          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
!!$        ENDDO
!!$        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
!!$      ENDDO
!!$                           CALL TRACE$POP
!!$      RETURN
!!$      END
!!$!
!!$!     ..................................................................
!!$      SUBROUTINE READIN_ANALYSE_POTCHARGEPOT(LL_CNTL_)
!!$!     ******************************************************************
!!$!     ****************************************************************** 
!!$      USE LINKEDLIST_MODULE
!!$      IMPLICIT NONE
!!$      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
!!$      TYPE(LL_TYPE)         :: LL_CNTL
!!$      LOGICAL(4)            :: TCHK
!!$      LOGICAL(4)            :: TON
!!$      TYPE(LL_TYPE)         :: LL_IMAGELIST
!!$      CHARACTER(256)   :: CH256SVAR1
!!$      CHARACTER(8)     :: CH8SVAR1
!!$!     ******************************************************************
!!$                           CALL TRACE$PUSH('READIN_ANALYSE_POINTCHARGEPOT')!!!$      LL_CNTL=LL_CNTL_
!!$      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
!!$      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
!!$      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
!!$      CALL LINKEDLIST$EXISTL(LL_CNTL,'POINTCHARGEPOT',1,TON)
!!$      IF(.NOT.TON) THEN
!!$        CALL TRACE$POP
!!$        RETURN
!!$      END IF
!!$      CALL ERROR$MSG('OPTIN POINTCHARGEPOT DISABLED')
!!$      CALL ERROR$MSG('USE A PAW TOOL')
!!$      CALL ERROR$STOP('READIN_ANALYSE_POINTCHARGEPOT')
!!$      CALL LINKEDLIST$SELECT(LL_CNTL,'POINTCHARGEPOT')
!!$!
!!$      CALL GRAPHICS$GETLIST(LL_IMAGELIST)    
!!$      CALL LINKEDLIST$SELECT(LL_IMAGELIST,'~')
!!$      CALL LINKEDLIST$SELECT(LL_IMAGELIST,'IMAGE',-1)
!!$!     __ GET TITLE__________________________________________________
!!$      CH256SVAR1='POTENTIAL OF POINT CHARGE MODEL'
!!$      CALL LINKEDLIST$SET(LL_IMAGELIST,'TITLE',0,CH256SVAR1)
!!$!     __ GET FILE___________________________________________________
!!$      CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,CH256SVAR1)
!!$      CALL LINKEDLIST$SET(LL_IMAGELIST,'FILE',0,CH256SVAR1)
!!$!     __ GET FILE_(DENSITY,WAVE)____________________________________
!!$      CH8SVAR1='PQV'
!!$      CALL LINKEDLIST$SET(LL_IMAGELIST,'TYPE',0,CH8SVAR1)
!!$                           CALL TRACE$POP
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_ANALYSE_OPTIC(LL_CNTL_)
!     ******************************************************************
!     ** CALL THE INTERFACES TO WRITE THE DATA FOR THE OPTIC CODE     **
!     ****************************************************************** 
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      LOGICAL(4)               :: TON
!     ******************************************************************
                           CALL TRACE$PUSH('READIN_ANALYSE_OPTIC')  
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ANALYSE')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'OPTIC',1,TON)
      IF(.NOT.TON) THEN
        CALL TRACE$POP
        RETURN
      END IF
      CALL ERROR$MSG('OPTIC MODULE HAS BEEN REMOVED TEMPORARILY IN THIS VERSION')
      CALL ERROR$STOP('READIN_ANALYSE_OPTIC')

CALL ERROR$MSG('OPTIC MODULE IS NOT PARALLELIZED')
CALL ERROR$MSG('USE ONLY IN SCALAR MODE')
CALL ERROR$STOP('READIN_ANALYSE_OPTIC')
!! COMMENTED OUT
!     CALL OPTICS$SETL4('ON',TON)
      RETURN
      END
!

!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READIN_DIMER(LL_CNTL_)
!     ******************************************************************
!     **  READ BLOCK !CONTROL!MERMIN                                  **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)         :: LL_CNTL
      LOGICAL(4)            :: TCHK
      LOGICAL(4)            :: TON
      REAL(8)               :: SVARR
      INTEGER(4)            :: SVARI,ITH
      LOGICAL(4)            :: SVARL
      CHARACTER(32)         :: SVARCH
      REAL(8),ALLOCATABLE   :: SVARV(:)
      INTEGER(4)            :: NFIL
!     ******************************************************************

      CALL TRACE$PUSH('READIN_DIMER')  
      LL_CNTL=LL_CNTL_

!
!     ==================================================================
!     == SELECT DIMER LIST OR EXIT IF IT DOES NOT EXIST              ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'DIMER',1,TON)

      IF(.NOT.TON) THEN
         CALL DIMER$SETL4('DIMER',.FALSE.)
         CALL TRACE$POP
         RETURN
      ELSE !WE HAVE A DIMER BLOCK
         CALL DIMER$SETL4('DIMER',.TRUE.)

!     ========INITIALIZE DEFAULT VALUES =====================        
         CALL DIMER$INIT()

         !=== FIND OUT IF WE READ A RESTART FILE -> NO NEED TO PLACE DIMER
         CALL LINKEDLIST$SELECT(LL_CNTL,'~')
         CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
         CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
         
         
         CALL LINKEDLIST$EXISTD(LL_CNTL,'START',0,TCHK)
         IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'START',1,SVARL)
            CALL DIMER$SETL4('PLACEDIMER',SVARL)
         ELSE
            !THE DEAULT IS START=F -> RESTART-FILE -> DO NOT PLACE THE DIMER
            CALL DIMER$SETL4('PLACEDIMER',.FALSE.)
         END IF
         

         !=== NOW READ THE DIMER SETTINGS
         CALL LINKEDLIST$SELECT(LL_CNTL,'~')
         CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
         CALL LINKEDLIST$SELECT(LL_CNTL,'DIMER')
         

!     ========GENERAL SETTINGS ==============================        
         CALL LINKEDLIST$EXISTD(LL_CNTL,'DIMERDIST',0,TCHK)
         IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'DIMERDIST',1,SVARR)
            CALL DIMER$SETR8('D',SVARR)
         END IF
         
         CALL LINKEDLIST$EXISTD(LL_CNTL,'KDLENGTH',0,TCHK)
         IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'KDLENGTH',1,SVARL)
            CALL DIMER$SETL4('KDLENGTH',SVARL)
         END IF
         
         CALL LINKEDLIST$EXISTD(LL_CNTL,'STRETCHDIST',0,TCHK)
         IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'STRETCHDIST',1,SVARR)
            CALL DIMER$SETR8('STRETCHDIST',SVARR)
            CALL DIMER$SETL4('STRETCH',.TRUE.)
         END IF
         
!     ========CONSTRAINTS ======= ==============================        
         CALL LINKEDLIST$SELECT(LL_CNTL,'~')
         CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
         CALL LINKEDLIST$SELECT(LL_CNTL,'DIMER')
         CALL LINKEDLIST$EXISTL(LL_CNTL,'CONSTRAINTS',1,TON)
         IF(TON) THEN
            CALL LINKEDLIST$SELECT(LL_CNTL,'CONSTRAINTS')
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'CENTERID',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'CENTERID',1,SVARCH)
               CALL DIMER$SETCH('CENTER_ID',SVARCH)
            END IF
            CALL LINKEDLIST$EXISTD(LL_CNTL,'CENTERCOORD',0,TCHK)
            IF(TCHK) THEN
               IF(ALLOCATED(SVARV)) DEALLOCATE(SVARV)
               ALLOCATE(SVARV(3))
               CALL LINKEDLIST$GET(LL_CNTL,'CENTERCOORD',1,SVARV)
               CALL DIMER$SETV('CENTERCOORD',SVARV,3)
            END IF
            CALL LINKEDLIST$EXISTD(LL_CNTL,'CONSTRSTEP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'CONSTRSTEP',1,SVARR)
               CALL DIMER$SETR8('CONSTRSTEP',SVARR)
            END IF
            
            CALL LINKEDLIST$NLISTS(LL_CNTL,'COUPLE',SVARI)
            DO ITH=1,SVARI
               CALL LINKEDLIST$SELECT(LL_CNTL,'COUPLE',ITH)
               CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOM',1,TCHK)
               IF(TCHK) THEN
                  CALL LINKEDLIST$GET(LL_CNTL,'ATOM',1,SVARCH)
                  CALL DIMER$CONSTRLIST_ADD(SVARCH)
               END IF
               CALL LINKEDLIST$SELECT(LL_CNTL,'..')            
            END DO
         END IF
         
!     ========ITERATION CONTROL==============================        
         CALL LINKEDLIST$SELECT(LL_CNTL,'~')
         CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
         CALL LINKEDLIST$SELECT(LL_CNTL,'DIMER')
         CALL LINKEDLIST$EXISTL(LL_CNTL,'ITERATION',1,TON)
         IF(TON) THEN
            CALL LINKEDLIST$SELECT(LL_CNTL,'ITERATION')
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'LITERMAX',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'LITERMAX',1,SVARI)
               CALL DIMER$SETI4('CALCMULTIPLIERITERMAX',SVARI)
            END IF
            CALL LINKEDLIST$EXISTD(LL_CNTL,'DLAMBDA',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'DLAMBDA',1,SVARR)
               CALL DIMER$SETR8('DLAMBDA',SVARR)
            END IF

!DEPRECATED
!!$         CALL LINKEDLIST$EXISTD(LL_CNTL,'VITER',0,TCHK)
!!$         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'VITER',1,CALCVELOCITYITERMAX)
!!$
!!$         CALL LINKEDLIST$EXISTD(LL_CNTL,'DVELOCITY',0,TCHK)
!!$         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'DVELOCITY',1,DVELOCITY)
!!$
!!$         CALL LINKEDLIST$EXISTD(LL_CNTL,'DROT',0,TCHK)
!!$         IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'DROT',1,DROT)


         END IF

!     ========MOTION CONTROL ================================        
         CALL LINKEDLIST$SELECT(LL_CNTL,'~')
         CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
         CALL LINKEDLIST$SELECT(LL_CNTL,'DIMER')
         CALL LINKEDLIST$EXISTL(LL_CNTL,'MOTION',1,TON)
         IF(TON) THEN
            CALL LINKEDLIST$SELECT(LL_CNTL,'MOTION')
            
!ALEX: THIS IS FOR PCLIMB-TESTING!
            CALL LINKEDLIST$EXISTD(LL_CNTL,'CLIMBPERP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'CLIMBPERP',1,SVARL)
               CALL DIMER$SETL4('CLIMBPERP',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FORCEDSTEP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FORCEDSTEP',1,SVARR)
               CALL DIMER$SETR8('FORCEDSTEP',SVARR)
            END IF
            
!ALEX: THIS IS FOR PCLIMB-TESTING!



            CALL LINKEDLIST$EXISTD(LL_CNTL,'DFOLLOWDOWN',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'DFOLLOWDOWN',1,SVARL)
               CALL DIMER$SETL4('DIMERFOLLOWDOWN',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'INHIBITUP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'INHIBITUP',1,SVARL)
               CALL DIMER$SETL4('INHIBITUP',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'INHIBITPERP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'INHIBITPERP',1,SVARL)
               CALL DIMER$SETL4('INHIBITPERP',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'ONLYROT',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'ONLYROT',1,SVARL)
               CALL DIMER$SETL4('ONLYROT',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'ONLYPERP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'ONLYPERP',1,SVARL)
               CALL DIMER$SETL4('ONLYPERP',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'WDOWN',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'WDOWN',1,SVARL)
               CALL DIMER$SETL4('WDOWN',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'WDOWNFACT',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'WDOWNFACT',1,SVARR)
               CALL DIMER$SETR8('WDOWNFACT',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FMPARA',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FMPARA',1,SVARR)
               CALL DIMER$SETR8('FMPARA',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FMPERP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FMPERP',1,SVARR)
               CALL DIMER$SETR8('FMPERP',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FMROT',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FMROT',1,SVARR)
               CALL DIMER$SETR8('FMROT',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FRICPARA',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FRICPARA',1,SVARR)
               CALL DIMER$SETR8('FRICPARA',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FRICPERP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FRICPERP',1,SVARR)
               CALL DIMER$SETR8('FRICPERP',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FRICROT',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FRICROT',1,SVARR)
               CALL DIMER$SETR8('FRICROT',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'OPTFRICPARA',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'OPTFRICPARA',1,SVARL)
               CALL DIMER$SETL4('OPTFRICPARA',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'OPTFRICPERP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'OPTFRICPERP',1,SVARL)
               CALL DIMER$SETL4('OPTFRICPERP',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'OPTFRICROT',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'OPTFRICROT',1,SVARL)
               CALL DIMER$SETL4('OPTFRICROT',SVARL)
            END IF


            CALL LINKEDLIST$EXISTD(LL_CNTL,'TMAXPARA',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'TMAXPARA',1,SVARR)
               CALL DIMER$SETR8('TMAXPARA',SVARR)
            END IF

            CALL LINKEDLIST$EXISTD(LL_CNTL,'TMAXPERP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'TMAXPERP',1,SVARR)
               CALL DIMER$SETR8('TMAXPERP',SVARR)
            END IF

            CALL LINKEDLIST$EXISTD(LL_CNTL,'TMAXROT',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'TMAXROT',1,SVARR)
               CALL DIMER$SETR8('TMAXROT',SVARR)
            END IF



            !FORCEAUTOPILOT
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FAUTOPARA',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FAUTOPARA',1,SVARL)
               CALL DIMER$SETL4('FAUTOPARA',SVARL)
            END IF
            

            CALL LINKEDLIST$EXISTD(LL_CNTL,'FAUTOPERP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FAUTOPERP',1,SVARL)
               CALL DIMER$SETL4('FAUTOPERP',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FAUTOROT',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FAUTOROT',1,SVARL)
               CALL DIMER$SETL4('FAUTOROT',SVARL)
            END IF

            CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC(+)PARA',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FRIC(+)PARA',1,SVARR)
               CALL DIMER$SETR8('FRICAUTOPARA',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC(+)PERP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FRIC(+)PERP',1,SVARR)
               CALL DIMER$SETR8('FRICAUTOPERP',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'FRIC(+)ROT',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'FRIC(+)ROT',1,SVARR)
               CALL DIMER$SETR8('FRICAUTOROT',SVARR)
            END IF
            
            
         END IF

         
!     ========LENGTH CONTROL ================================        
         CALL LINKEDLIST$SELECT(LL_CNTL,'~')
         CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
         CALL LINKEDLIST$SELECT(LL_CNTL,'DIMER')
         CALL LINKEDLIST$EXISTL(LL_CNTL,'FLEXLENGTH',1,TON)
         IF(TON) THEN
            CALL LINKEDLIST$SELECT(LL_CNTL,'FLEXLENGTH')
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'DLFLEX',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'DLFLEX',1,SVARL)
               CALL DIMER$SETL4('DLFLEX',SVARL)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'LCS',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'LCS',1,SVARI)
               CALL DIMER$SETI4('LCS',SVARI)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'CENTERDIFFMIN',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'CENTERDIFFMIN',1,SVARR)
               CALL DIMER$SETR8('RCDIFFMIN',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'DSTEP',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'DSTEP',1,SVARR)
               CALL DIMER$SETR8('DSTEP',SVARR)
            END IF

            CALL LINKEDLIST$EXISTD(LL_CNTL,'NSTEPS',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'NSTEPS',1,SVARI)
               CALL DIMER$SETI4('NSTEPS',SVARI)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'DTOBE',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'DTOBE',1,SVARR)
               CALL DIMER$SETR8('DMIN',SVARR)
            END IF
         END IF
         



!     ========OUTPUT CONTROL ================================        
         CALL LINKEDLIST$SELECT(LL_CNTL,'~')
         CALL LINKEDLIST$SELECT(LL_CNTL,'CONTROL')
         CALL LINKEDLIST$SELECT(LL_CNTL,'DIMER')
         CALL LINKEDLIST$EXISTL(LL_CNTL,'OUTPUT',1,TON)
         IF(TON) THEN
            CALL LINKEDLIST$SELECT(LL_CNTL,'OUTPUT')
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'ENERGYTRA',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'ENERGYTRA',1,SVARR)
               CALL DIMER$SETR8('ENERGYTRA',SVARR)
            END IF
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'NATOMS',0,TCHK)
            IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'NATOMS',1,SVARI)
            SVARI=SVARI*3
            IF(ALLOCATED(SVARV)) DEALLOCATE(SVARV)
            ALLOCATE(SVARV(SVARI))
            
            CALL LINKEDLIST$EXISTD(LL_CNTL,'RTS',0,TCHK)
            IF(TCHK) THEN
               CALL LINKEDLIST$GET(LL_CNTL,'RTS',1,SVARV)
               CALL DIMER$SETV('RTS',SVARV,SVARI)
            END IF
         END IF
      END IF !FROM DIMER BLOCK
      
      CALL FILEHANDLER$UNIT('DPROT',NFIL)
      CALL DIMER$REPORT_SETTINGS(NFIL)
      CALL TRACE$POP
      RETURN
    END SUBROUTINE READIN_DIMER
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN
!     **************************************************************************
!     **                                                                      **
!     **  READ STRUCTURAL DATA                                                **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1995 ************
      USE IO_MODULE
      USE LINKEDLIST_MODULE
      USE CONSTANTS_MODULE
      IMPLICIT NONE
      INTEGER(4)            :: NFILO   ! PROTOCOL FILE UNIT
      INTEGER(4)            :: NFIL
      LOGICAL(4)            :: TCHK,TCHK1
      INTEGER(4)            :: NKPT
      REAL(8)               :: ANGSTROM
      REAL(8)               :: SVAR
!     **************************************************************************
                          CALL TRACE$PUSH('STRCIN')
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==========================================================================
!     ==  SELECT STRUCTURE FILE AND READ IT INTO THE BUFFER                   ==
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_STRC)
      CALL FILEHANDLER$UNIT('STRC',NFIL)
      CALL LINKEDLIST$READ(LL_STRC,NFIL,'MONOMER')
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
!    
!     ==========================================================================
!     ==  MARK ALL ELEMENTS AS READ FROM INPUT FILE                           ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$MARK(LL_STRC,1)
!
!     ==========================================================================
!     ==  ENTER BLOCK !STRUCTURE                                              ==
!     ==========================================================================
      
      CALL LINKEDLIST$EXISTL(LL_STRC,'STRUCTURE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('DATA FOR BLOCK STRUCTURE MISSING')
        CALL ERROR$STOP('STRCIN')
      END IF
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
!    
!     ==========================================================================
!     ==  READ BLOCK !STRUCTURE!GENERIC                                       ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')
!
!     ==  READ ACTUAL VALUES  ==================================================
!
!     ==  READ LUNIT ===========================================================
      CALL LINKEDLIST$EXISTD(LL_STRC,'LUNIT[AA]',1,TCHK1)
      IF(TCHK1) THEN 
        CALL LINKEDLIST$EXISTD(LL_STRC,'LUNIT',1,TCHK)
        IF(TCHK.AND.TCHK1) THEN
          CALL ERROR$MSG('LUNIT AND LUNIT[AA] ARE MUTUALLY EXCLUSIVE')
          CALL ERROR$STOP('STRCIN')
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'LUNIT[AA]',0,SVAR)
        CALL CONSTANTS$GET('ANGSTROM ',ANGSTROM)
        SVAR=SVAR*ANGSTROM
        CALL LINKEDLIST$SET(LL_STRC,'LUNIT',0,SVAR)
        CALL LINKEDLIST$RMDATA(LL_STRC,'LUNIT[AA]',1)
      END IF
      CALL LINKEDLIST$EXISTD(LL_STRC,'LUNIT',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'LUNIT',0,1.D0)
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
!    
!     ==========================================================================
!     ==  READ BLOCK !STRUCTURE!SPECIES                                       ==
!     ==========================================================================
      CALL STRCIN_SPECIES(LL_STRC)
!    
!     ==========================================================================
!     ==  READ BLOCK !STRUCTURE!LATTICE                                       ==
!     ==========================================================================
      CALL STRCIN_LATTICE(LL_STRC)
!    
!     ==========================================================================
!     ==  READ BLOCK !STRUCTURE!KPOINT                                        ==
!     ==========================================================================
      CALL STRCIN_KPOINT(LL_STRC,NKPT)
!
!     ==========================================================================
!     ==  ENTER BLOCK !ATOM                                                   ==
!     ==========================================================================
      CALL STRCIN_ATOM(LL_STRC)
!
!     ==================================================================
!     ==  ENTER BLOCK !GROUP                                          ==
!     ==================================================================
      CALL STRCIN_GROUP(LL_STRC)
!
!     ==================================================================
!     ==  ENTER BLOCK !STRC!OCCUP                                     ==
!     ==================================================================
      CALL STRCIN_OCCUP(LL_STRC,NKPT)
!
!     ==================================================================
!     ==  ENTER BLOCK !ISOLATE                                        ==
!     ==================================================================
      CALL STRCIN_ISOLATE(LL_STRC)
!
!     ==================================================================
!     ==  ENTER BLOCK !CONFINE                                        ==
!     ==================================================================
      CALL STRCIN_CONFINE(LL_STRC)
!
!     ==================================================================
!     ==  ENTER BLOCK !ORBPOT                                         ==
!     ==================================================================
      CALL STRCIN_ORBPOT(LL_STRC)
!
!     ==================================================================
!     ==  ENTER BLOCK !VEXT                                           ==
!     ==================================================================
      CALL STRCIN_VEXT(LL_STRC)
!
!     ==================================================================
!     ==  ENTER BLOCK !STRUCTURE!QM-MM                                ==
!     ==================================================================
      CALL STRCIN_SOLVENT(LL_STRC)
!
!     ==================================================================
!     ==  ENTER BLOCK !STRUCTURE!COSMO                                ==
!     ==================================================================
      CALL STRCIN_COSMO(LL_STRC)
!
!     ==================================================================
!     ==  ENTER BLOCK !STRUCTURE!CONSTRAINTS!RIGID                    ==
!     ================================================================== 
      CALL STRCIN_CONSTRAINTS_RIGID(LL_STRC)
      CALL STRCIN_CONSTRAINTS_FREEZE(LL_STRC)
      CALL STRCIN_CONSTRAINTS_LINEAR(LL_STRC)
      CALL STRCIN_CONSTRAINTS_BOND(LL_STRC)
      CALL STRCIN_CONSTRAINTS_MIDPLANE(LL_STRC)
      CALL STRCIN_CONSTRAINTS_TRANSLATION(LL_STRC)
      CALL STRCIN_CONSTRAINTS_ROTATION(LL_STRC)
      CALL STRCIN_CONSTRAINTS_ORIENTATION(LL_STRC)
      CALL STRCIN_CONSTRAINTS_COGSEP(LL_STRC)
      CALL STRCIN_CONSTRAINTS_ANGLE(LL_STRC)
      CALL STRCIN_CONSTRAINTS_TORSION(LL_STRC)
                          CALL TRACE$PASS('BLOCK !STRUCTURE FINISHED')
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$REPORT_UNUSED(LL_STRC,NFILO)

      CALL FILEHANDLER$CLOSE('STRC')
      CALL FILEHANDLER$CLOSE('PROT')
                          CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_LATTICE(LL_STRC_)
!     ****************************************************************** 
!     **  DEFINES LATTICE IN ATOMLIST OBJECT                          **
!     **                                                              **
!     **  REQUIRES PREDEFINED: NOTHING                                **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK
      REAL(8)                  :: UNIT
      REAL(8)                  :: RBAS(3,3)
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_LATTICE')  
      LL_STRC=LL_STRC_
!    
!     ==================================================================
!     ==  READ UNIT FROM BLOCK !STRUCTURE!GENERIC                     ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')
!
      CALL LINKEDLIST$GET(LL_STRC,'LUNIT',1,UNIT)
!    
!     ==================================================================
!     ==  READ BLOCK !STRUCTURE!LATTICE                               ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'LATTICE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('BLOCK !STRUCTURE!LATTICE IS MANDATORY')
        CALL ERROR$STOP('STRCIN_LATTICE')
      END IF
      CALL LINKEDLIST$SELECT(LL_STRC,'LATTICE')
!
!     == DEFAULT VALUES  ===============================================
      CALL LINKEDLIST$EXISTD(LL_STRC,'T',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('KEYWORD T= IS MANDATORY')
        CALL ERROR$STOP('STRCIN_LATTICE; !STRUCTURE!LATTICE')
      END IF
!
!     ==  READ ACTUAL VALUES  ==========================================
      CALL LINKEDLIST$GET(LL_STRC,'T',1,RBAS)
      RBAS(:,:)=RBAS(:,:)*UNIT
      CALL CELL$SETR8A('TREF',9,RBAS)
      CALL CELL$SETR8A('T0',9,RBAS)
      CALL CELL$SETR8A('TM',9,RBAS)
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_KPOINT(LL_STRC_,NKPT)
!     **************************************************************************
!     **  DEFINES THE K-POINT INFO OF OCCUPATIONS_MODULE                      **
!     **                                                                      **
!     **  REQUIRES PREDEFINED: NOTHING                                        ** 
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      INTEGER(4)   ,INTENT(OUT):: NKPT
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK,TCHK1
      INTEGER(4)               :: NKDIV(3)
      INTEGER(4)               :: ISHIFT(3)=0
      REAL(8)                  :: RBAS(3,3)
      REAL(8)                  :: RMAX
      REAL(8)    ,ALLOCATABLE  :: XK(:,:)
      REAL(8)    ,ALLOCATABLE  :: WGHT(:)
      LOGICAL(4)               :: TINV     ! EXPLOIT TIME INVERSION SYMMETRY
      INTEGER(4)               :: NSPIN,IKPT
!     **************************************************************************
                           CALL TRACE$PUSH('STRCIN_KPOINT')
      LL_STRC=LL_STRC_
!
!     ==========================================================================
!     == CHECK IF CALCULATION IS NON-COLLINEAR                                ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'OCCUPATIONS')
      CALL LINKEDLIST$EXISTD(LL_STRC,'NSPIN',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'NSPIN',0,1)
      CALL LINKEDLIST$GET(LL_STRC,'NSPIN',1,NSPIN)
      IF(NSPIN.NE.1.AND.NSPIN.NE.2.AND.NSPIN.NE.3) THEN
        CALL ERROR$MSG('NUMBER OF SPINS OUT OF RANGE')
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('STRCIN_KPOINT')
      END IF 
      IF(NSPIN.EQ.3) THEN 
        TINV=.FALSE.  ! TIME INVERSION SYMMETRY CANNOT BE EXPLOITED IN THE 
                      ! WITH THE PAULI EQUATION
      ELSE
        TINV=.TRUE. 
      END IF    
!
!     ==========================================================================
!     == FIRST SET FOR ISOLATED MOLECULE, THEN OVERWRITE FRO K-POINTS         ==
!     ==========================================================================
      CALL KPOINTS$SETL4('TINV',TINV)
      CALL KPOINTS$SETI4A('DIV',3,(/1,1,1/))
      CALL KPOINTS$SETI4A('SHIFT',3,(/0,0,0/))
!
!     ==========================================================================
!     == NOW ENTER K-POINT BLOCK                                              ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'KPOINTS',1,TCHK)
      CALL LINKEDLIST$SELECT(LL_STRC,'KPOINTS')
!     == SET DEFAULT ===========================================================
      IF(.NOT.TCHK) THEN
        NKDIV(:)=1
        CALL LINKEDLIST$SET(LL_STRC,'DIV',0,NKDIV) 
      END IF
!
!     ==  READ ACTUAL VALUES  ==================================================
      ISHIFT(:)=0
      CALL LINKEDLIST$EXISTD(LL_STRC,'SHIFT',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_STRC,'SHIFT',1,ISHIFT)
      CALL KPOINTS$SETI4A('SHIFT',3,ISHIFT)
!
      CALL LINKEDLIST$NDATA(LL_STRC,'K',NKPT)
      IF(NKPT.NE.0) THEN
        CALL ERROR$MSG('CHOICE OF INDIVIDUALLY SPECIFIED KPOINTS HAS BEEN DISABLED')
        CALL ERROR$MSG('OPTION STRUCTURE!KPOINT:K IS NO MORE ALLOWED')
        CALL ERROR$STOP('STRCIN_KPOINT')
      END IF
!
!     ==========================================================================
!     == DETERMINE DIVISIONS NKDIV OF RECIPROCAL UNIT CELL                    ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_STRC,'R',1,TCHK)
      CALL LINKEDLIST$EXISTD(LL_STRC,'DIV',1,TCHK1)
      IF(TCHK1.AND.TCHK) THEN
        CALL ERROR$MSG('!STRUCTURE!KPOINTS:R AND DIV ARE MUTUALLY EXCLUSIVE')
        CALL ERROR$STOP('STRCIN_KPOINT')
      ELSE IF(TCHK) THEN
!       == K-POINT GRID DEFINED BY R ===========================================
        CALL CELL$GETR8A('TREF',9,RBAS)
        CALL LINKEDLIST$GET(LL_STRC,'R',1,RMAX)
        CALL KPOINTS_NKDIV(RBAS,RMAX,NKDIV)
      ELSE IF(TCHK1) THEN
        CALL LINKEDLIST$GET(LL_STRC,'DIV',1,NKDIV)
        IF(NKDIV(1).LE.0.OR.NKDIV(2).LE.0.OR.NKDIV(3).LE.0) THEN
          CALL ERROR$MSG('!STRUCTURE!KPOINT:DIV MUST BE GREATER THAN ZERO')
          CALL ERROR$STOP('STRCIN_KPOINT')
        END IF
      ELSE
        NKDIV(:)=2
      END IF
      CALL KPOINTS$SETI4A('DIV',3,NKDIV)
!
!     ==========================================================================
!     == DETERMINE K-POINTS AND INTEGRATION WEIGHTS                           ==
!     ==========================================================================
      CALL KPOINTS$MAKEGRID()
      CALL KPOINTS$GETI4('NKPT',NKPT)

!      CALL KPOINTS_NKPT(TINV,NKDIV,ISHIFT,NKPT)
      ALLOCATE(XK(3,NKPT))
      ALLOCATE(WGHT(NKPT))
!      CALL KPOINTS_KPOINTS(TINV,NKDIV,ISHIFT,NKPT,XK,WGHT)
      CALL KPOINTS$GETR8A('XK',3*NKPT,XK)
      CALL KPOINTS$GETR8A('WKPT',NKPT,WGHT)
!
!     ==  PERFORM ACTIONS  =====================================================
      CALL DYNOCC$SETI4('NKPT',NKPT)
      CALL DYNOCC$SETI4A('NKDIV',3,NKDIV)
      CALL DYNOCC$SETI4A('ISHIFT',3,ISHIFT)
      CALL DYNOCC$SETR8A('XK',3*NKPT,XK) 
      CALL DYNOCC$SETR8A('WKPT',NKPT,WGHT)
      DEALLOCATE(XK)
      DEALLOCATE(WGHT)
                           CALL TRACE$POP
      RETURN
      END
!     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_SPECIES(LL_STRC_)
!     **************************************************************************
!     **                                                                      **
!     **  REQUIRES PREDEFINED: NOTHING                                        **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: NSP
      INTEGER(4)               :: ISP
!     **************************************************************************
                           CALL TRACE$PUSH('STRCIN_SPECIES')
      LL_STRC=LL_STRC_
!
!     ==========================================================================
!     == NOW ANALYZE THE SPECIES BLOCKS AND CONFIGURE SETUP INSTANCE
!     ==========================================================================
      CALL SETUP$READSTRCIN(LL_STRC)
!    
!     ==========================================================================
!     ==  READ BLOCK !STRUCTURE!NTBO                                          ==
!     ==========================================================================
      CALL STRCIN_LMTO(LL_STRC)
!    
!     ==========================================================================
!     ==  CONFIGURE SETUP OBJECT  (THE NAME "READ" IS MISLEADING)             ==
!     ==========================================================================
      CALL SETUP$READ()
!
!     ==========================================================================
!     == WARN BECAUSE OF OBSOLETE KEYWORDS LDAPLUSU AND HYBRID                ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$NLISTS(LL_STRC,'SPECIES',NSP)
      IF(NSP.EQ.0) THEN
        CALL ERROR$MSG('NO ATOM TYPES !STRUCTURE!SPECIES  SPECIFIED.')
        CALL ERROR$MSG('EVEN IF NOT USED, ONE ATOM TYPE MUST BE SPECIFIED')
        CALL ERROR$MSG('TO AVOID INTERNAL PROBLEMS')
!       == THE CODE RUNS OK UNTIL POTENTIAL_VOFRHO. PROBABLY THE CODE TRIES ====
!       == TO ACCESS AN ARRAY ELEMENT OF AN ARRAY OF ZERO SIZE =================
        CALL ERROR$STOP('STRCIN_SPECIES')
      END IF

      DO ISP=1,NSP
        CALL LINKEDLIST$SELECT(LL_STRC,'SPECIES',ISP)
        CALL LINKEDLIST$EXISTL(LL_STRC,'LDAPLUSU',1,TCHK)
        IF(TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!SPECIES:LDAPLUSU IS OBSOLETE')
          CALL ERROR$MSG('STRCIN_SPECIES')
        END IF
        CALL LINKEDLIST$EXISTL(LL_STRC,'HYBRID',1,TCHK)
        IF(TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!SPECIES:HYBRID IS OBSOLETE')
          CALL ERROR$MSG('STRCIN_SPECIES')
        END IF
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_LMTO(LL_STRC_)
!     **************************************************************************
!     **  DEFINES THE PARAMETERS FOR THE LMTO OBJECT                          **
!     **                                                                      **
!     **  !SPECIES MUST BE FINISHED BEFORE THIS IS CALLED                     **
!     **                                                                      **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)               :: ISP
      CHARACTER(32)            :: SPNAME     !SPECIES NAME
      LOGICAL(4)               :: TCHK,TCHK1
      INTEGER(4)               :: NSP        !#(SPECIES)
      REAL(8)                  :: SVAR
      INTEGER(4)               :: LENG
      REAL(8)   ,ALLOCATABLE   :: WORK(:)
      INTEGER(4),ALLOCATABLE   :: IWORK(:)
      REAL(8)                  :: AEZ       ! ATOMIC NUMBER
      REAL(8)                  :: RCOV      ! COVALENT RADIUS
!     **************************************************************************
                           CALL TRACE$PUSH('STRCIN_SPECIES_NTBO')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$NLISTS(LL_STRC,'SPECIES',NSP)
      DO ISP=1,NSP
        CALL LINKEDLIST$SELECT(LL_STRC,'SPECIES',ISP)
        CALL LMTO$SETI4('ISP ',ISP)
        CALL SIMPLELMTO$SETI4('ISP ',ISP)
!
!       ========================================================================
!       ==  SKIP IF NTBO BLOCK IS NOT PRESENT                                 ==
!       ========================================================================
        CALL LINKEDLIST$EXISTL(LL_STRC,'NTBO',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL LMTO$SETR8('LHFWEIGHT',0.D0)          ! SWITCH CONTRIBUTION OFF
          CALL SIMPLELMTO$SETR8('LHFWEIGHT',0.D0)    ! SWITCH CONTRIBUTION OFF
          CALL LINKEDLIST$SELECT(LL_STRC,'..')       ! LEAVE SPECIES BLOCK
          CYCLE
        END IF
!
!       ========================================================================
!       ==  SPECIES NAME / COLLECT INFORMATION FROM SETUP OBJECT              ==
!       ========================================================================
        CALL LINKEDLIST$GET(LL_STRC,'NAME',1,SPNAME)
        CALL SETUP$SELECT(SPNAME)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
        CALL SETUP$UNSELECT()
!
!       ========================================================================
!       ==  ENTER !SPECIES!NTBO BLOCK                                         ==
!       ========================================================================
        CALL LINKEDLIST$SELECT(LL_STRC,'NTBO')
!   
!       ========================================================================
!       == DEFINE LOCAL ORBITALS                                              ==
!       ========================================================================
!       ==  SELECTOR FOR LOCAL ORBITALS: #(ORBITALS PER L) =====================
        CALL LINKEDLIST$EXISTD(LL_STRC,'NOFL',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$SIZE(LL_STRC,'NOFL',1,LENG)
           ALLOCATE(IWORK(LENG))
           CALL LINKEDLIST$GET(LL_STRC,'NOFL',1,IWORK)
        ELSE
          LENG=1
          ALLOCATE(IWORK(LENG))   ! THE ARRAY WILL BE EXTENDED
          IWORK(1)=0           ! A SUM EQUAL ZERO SWITCHES TO TORB (SEE BELOW)
        END IF        
        CALL LMTO$SETI4A('NORBOFL',LENG,IWORK)
        CALL SIMPLELMTO$SETI4A('NORBOFL',LENG,IWORK)
        DEALLOCATE(IWORK)
!   
!       == RADIUS FOR MATCHING PARTIAL WAVES TO ENVELOPE FUNCTIONS =============
        CALL LINKEDLIST$EXISTD(LL_STRC,'RAUG/RCOV',1,TCHK)
        SVAR=1.D0
        IF(TCHK)CALL LINKEDLIST$GET(LL_STRC,'RAUG/RCOV',1,SVAR)
        SVAR=SVAR*RCOV
        CALL LMTO$SETR8('RAUG',SVAR)
        CALL SIMPLELMTO$SETR8('RAUG',SVAR)
!   
!       == RADIUS FOR MATCHING EXPONENTIAL TAILS ===============================
        SVAR=1.D0
        CALL LINKEDLIST$EXISTD(LL_STRC,'RTAIL/RCOV',1,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_STRC,'RTAIL/RCOV',1,SVAR)
        SVAR=SVAR*RCOV
        CALL LMTO$SETR8('RTAIL',SVAR)
!   
!       == DECAY CONSTANTS FOR THE EXPONENTIAL TAILS ===========================
        CALL LINKEDLIST$EXISTD(LL_STRC,'TAILLAMBDA',1,TCHK)
        IF(TCHK) THEN
          ALLOCATE(WORK(2))
          CALL LINKEDLIST$GET(LL_STRC,'TAILLAMBDA',1,WORK)
          IF(WORK(1).LE.0.D0.OR.WORK(2).LE.0.D0) THEN
            CALL ERROR$MSG('INVALID INPUT IN !STRUCTURE!SPECIES!NTBO')
            CALL ERROR$MSG('TAILLAMBDA MUST BE POSITIVE')
            CALL ERROR$STOP('STRCIN_LMTO')
          END IF
          CALL LMTO$SETR8('TAILLAMBDA1',WORK(1))
          CALL LMTO$SETR8('TAILLAMBDA2',WORK(2))
          DEALLOCATE(WORK)
        END IF
!
!       == TAIL-TRUNCATION RADIUS ==============================================
        SVAR=99.D0 ! DEFAULT IN UNITS OF RCOV
        CALL LINKEDLIST$EXISTD(LL_STRC,'RTAILCUT/RCOV',1,TCHK)
        IF(TCHK) CALL LINKEDLIST$GET(LL_STRC,'RTAILCUT/RCOV',1,SVAR)
        SVAR=SVAR*RCOV
        CALL LMTO$SETR8('RTAILCUT',SVAR)
!   
!       ========================================================================
!       == DEFINE INTERACTIONS                                                ==
!       ========================================================================
!       == LOCAL VALUE FOR THE ADMIXING FACTOR FOR THE LOCAL INTERACTION =======
        CALL LINKEDLIST$EXISTD(LL_STRC,'LHFWEIGHT',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'LHFWEIGHT',1,SVAR)
          CALL LMTO$SETR8('LHFWEIGHT',SVAR)
          CALL SIMPLELMTO$SETR8('LHFWEIGHT',SVAR)
        END IF
!
!       == CORE VALENCE EXCHANGE ===============================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'CV',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'CV',1,TCHK1)
          CALL LMTO$SETL4('COREVALENCE',TCHK1)
          CALL SIMPLELMTO$SETL4('COREVALENCE',TCHK1)
        END IF
!   
!       == FOCK SETUP ==========================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'FOCKSETUP',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'FOCKSETUP',1,TCHK1)
          CALL LMTO$SETL4('FOCKSETUP',TCHK1)
          CALL SIMPLELMTO$SETL4('FOCKSETUP',TCHK1)
        END IF
!   
!       == NEGLECT OF DIATOMIC DIFFERENTIAL OVERLAP ===========================
        CALL LINKEDLIST$EXISTD(LL_STRC,'NDDO',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'NDDO',1,TCHK1)
          CALL LMTO$SETL4('NDDO',TCHK1)
          CALL SIMPLELMTO$SETL4('NDDO',TCHK1)
        END IF
!   
!       ==  
        CALL LINKEDLIST$EXISTD(LL_STRC,'31',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'31',1,TCHK1)
          CALL LMTO$SETL4('31',TCHK1)
          CALL SIMPLELMTO$SETL4('31',TCHK1)
        END IF
!
!       == 
        CALL LINKEDLIST$EXISTD(LL_STRC,'BONDX',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'BONDX',1,TCHK1)
          CALL LMTO$SETL4('BONDX',TCHK1)
          CALL SIMPLELMTO$SETL4('BONDX',TCHK1)
        END IF
!   
!       ========================================================================
!       ==  OBSOLETE KEYWORDS                                                 ==
!       ========================================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'NTBO',1,TCHK)
        IF(TCHK) THEN
          CALL ERROR$MSG('NTBO IS OBSOLETE. USE NOFL TO SET NUMBER OF ORBITALS')
          CALL ERROR$STOP('STRCIN_LMTO')
        END IF

        CALL LINKEDLIST$EXISTD(LL_STRC,'TORB',1,TCHK)
        IF(TCHK) THEN
          CALL ERROR$MSG('TORB IS OBSOLETE. USE NOFL TO SET NUMBER OF ORBITALS')
          CALL ERROR$STOP('STRCIN_LMTO')
        END IF
!   
!       ========================================================================
!       ==  CLOSE LOOP                                                        ==
!       ========================================================================
        CALL LMTO$SETI4('ISP',0)             ! UNSELECT LMTO INSTANCE
        CALL SIMPLELMTO$SETI4('ISP',0)             ! UNSELECT LMTO INSTANCE
        CALL LINKEDLIST$SELECT(LL_STRC,'..') ! LEAVE NTBO BLOCK
        CALL LINKEDLIST$SELECT(LL_STRC,'..') ! LEAVE SPECIES BLOCK
      ENDDO
                                   CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_ATOM(LL_STRC_)
!     ******************************************************************
!     **  DEFINES THE ATOMLIST                                        **
!     **                                                              **
!     **  MODIFIES           : ATOMLIST                               **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: NAT
      INTEGER(4)               :: IAT,IAT1
      INTEGER(4)               :: ISP
      INTEGER(4)               :: NSP
      CHARACTER(32)            :: ATOM,SPECIES,STRING
      CHARACTER(32),ALLOCATABLE:: NAME(:)
      REAL(8)                  :: SVAR
      REAL(8)                  :: UNIT
      REAL(8)                  :: POS(3)
      REAL(8)                  :: PROTONMASS
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_ATOM')
      LL_STRC=LL_STRC_
      CALL CONSTANTS('U',PROTONMASS)
      CALL SETUP$GETI4('NSP',NSP)
!    
!     ==================================================================
!     ==  READ BLOCK !STRUCTURE!GENERIC                               ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')
      CALL LINKEDLIST$GET(LL_STRC,'LUNIT',1,UNIT)
!    
!     ==================================================================
!     ==  READ BLOCK !STRUCTURE!ATOM                                  ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NAT)
      ALLOCATE(NAME(NAT))
      CALL ATOMLIST$INITIALIZE(NAT)
      DO IAT=1,NAT
        CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT)
        CALL LINKEDLIST$GET(LL_STRC,'NAME',1,ATOM)
        CALL ATOMLIST$SETCH('NAME',IAT,ATOM)
        NAME(IAT)=ATOM
!  
        CALL LINKEDLIST$EXISTD(LL_STRC,'SP',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL LINKEDLIST$SET(LL_STRC,'SP',0,ATOM(1:2))
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'SP',1,SPECIES)
!
        TCHK=.FALSE.
        DO ISP=1,NSP
          CALL SETUP$ISELECT(ISP)
          CALL SETUP$GETCH('ID',STRING)
          IF(STRING.EQ.SPECIES) THEN
            TCHK=.TRUE.
            CALL ATOMLIST$SETI4('ISPECIES',IAT,ISP)
            CALL SETUP$GETR8('M',SVAR)
            CALL ATOMLIST$SETR8('MASS',IAT,SVAR)
            CALL SETUP$GETR8('PS<G2>',SVAR)
            CALL ATOMLIST$SETR8('PS<G2>',IAT,SVAR)
            CALL SETUP$GETR8('PS<G4>',SVAR)
            CALL ATOMLIST$SETR8('PS<G4>',IAT,SVAR)
          END IF
          CALL SETUP$UNSELECT()
        ENDDO 
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('NO SPECIES FOR ATOM')
          CALL ERROR$MSG('IN BLOCK !STRUCTURE!ATOM')
          CALL ERROR$CHVAL('ATOM NAME',ATOM)
          CALL ERROR$STOP('STRCIN_ATOM')
        END IF
!
!       __ATOMIC POSITION_______________________________________________
        CALL LINKEDLIST$GET(LL_STRC,'R',1,POS)
        POS(:)=POS(:)*UNIT
        CALL ATOMLIST$SETR8A('R(0)',IAT,3,POS) 
        CALL ATOMLIST$SETR8A('R(-)',IAT,3,POS)
        CALL ATOMLIST$SETR8A('R(+)',IAT,3,POS)
!
!       __ATOMIC MASS___________________________________________________
        CALL LINKEDLIST$EXISTD(LL_STRC,'M',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'M',1,SVAR)
          CALL ATOMLIST$SETR8('MASS',IAT,SVAR*PROTONMASS)
        END IF
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
!
!     ==  CHECKS AND TRANSFORMATIONS =================================
      CALL ATOMLIST$ORDER
      DO IAT=1,NAT
        DO IAT1=IAT+1,NAT
          IF(NAME(IAT).EQ.NAME(IAT1)) THEN
            CALL ERROR$MSG('ATOM NAMES MUST BE DIFFERENT')
            CALL ERROR$CHVAL('NAME',NAME(IAT))
            CALL ERROR$I4VAL('IAT1',IAT)
            CALL ERROR$I4VAL('IAT2',IAT1)
            CALL ERROR$STOP('STRCIN_ATOM')
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(NAME)
!
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_ATOM
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_ISOLATE(LL_STRC_)
!     ******************************************************************
!     **  DEFINES THE ISOLATE OBJECT                                  **
!     **                                                              **
!     **  REQUIRES PREDEFINED: NOTHING                                **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK
      LOGICAL(4)               :: TON
      INTEGER(4)               :: NFCT
      REAL(8)                  :: RC
      REAL(8)                  :: RCFAC
      REAL(8)                  :: GMAX2
      LOGICAL(4)               :: DECOUPLE
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_ISOLATE')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'ISOLATE',1,TON)
      IF(.NOT.TON) THEN
        CALL ISOLATE$SETL4('DECOUPLE',.FALSE.)
        CALL TRACE$POP
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_STRC,'ISOLATE',0)
      CALL ISOLATE$SETL4('ON',.TRUE.)
!
!     == SET DEFAULTS ==================================================
      CALL LINKEDLIST$EXISTD(LL_STRC,'NF',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'NF',0,3)
      CALL LINKEDLIST$GET(LL_STRC,'NF',1,NFCT)
      CALL ISOLATE$SETI4('NFCT',NFCT)

      CALL LINKEDLIST$EXISTD(LL_STRC,'RC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'RC',0,0.5D0)
      CALL LINKEDLIST$GET(LL_STRC,'RC',1,RC)
      CALL ISOLATE$SETR8('RC0',RC)

      CALL LINKEDLIST$EXISTD(LL_STRC,'RCFAC',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'RCFAC',0,1.5D0)
      CALL LINKEDLIST$GET(LL_STRC,'RCFAC',1,RCFAC)
      CALL ISOLATE$SETR8('RCFAC',RCFAC)

      CALL LINKEDLIST$EXISTD(LL_STRC,'GMAX2',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'GMAX2',0,3.D0)
      CALL LINKEDLIST$GET(LL_STRC,'GMAX2',1,GMAX2)
      CALL ISOLATE$SETR8('G2MAX',GMAX2)
!
      CALL LINKEDLIST$EXISTD(LL_STRC,'DECOUPLE',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'DECOUPLE',0,.TRUE.)
      CALL LINKEDLIST$GET(LL_STRC,'DECOUPLE',1,DECOUPLE)
      CALL ISOLATE$SETL4('DECOUPLE',DECOUPLE)
!
                           CALL TRACE$POP
      RETURN
      END SUBROUTINE STRCIN_ISOLATE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONFINE(LL_STRC_)
!     **************************************************************************
!     **  DEFINES THE CONFINING POTENTIAL MIMICING THE PAULI REPULSION        **
!     **  OF A SOLVENT                                                        **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK
      LOGICAL(4)               :: TON
      REAL(8)                  :: SVAR
!     **************************************************************************
                           CALL TRACE$PUSH('STRCIN_ISOLATE')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'CONFINE',1,TON)
      CALL POTENTIAL$SETL4('CONFINE',TON)
      IF(.NOT.TON) THEN
        CALL TRACE$POP
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_STRC,'CONFINE',0)
      
      CALL LINKEDLIST$EXISTD(LL_STRC,'POT',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'POT',0,20.D0)
      CALL LINKEDLIST$GET(LL_STRC,'POT',1,SVAR)
      CALL POTENTIAL$SETR8('VCONFINE',SVAR)
                           CALL TRACE$POP
      RETURN
      END SUBROUTINE STRCIN_CONFINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_OCCUP(LL_STRC_,NKPT)
!     ******************************************************************
!     **  DEFINES THE DYNOCC OBJECT                                   **
!     **                                                              **
!     **  REQUIRES PREDEFINED: FILEHANDLER,ATOMLIST                   **
!     **  MODIFIES           :: DFT,DYNOCC                            **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN):: LL_STRC_
      INTEGER(4)   ,INTENT(IN):: NKPT
      TYPE(LL_TYPE)           :: LL_STRC
      LOGICAL(4)              :: TCHK
      INTEGER(4)              :: NB
      INTEGER(4)              :: NSPIN
      INTEGER(4)              :: NTH
      INTEGER(4)              :: NAT
      INTEGER(4)              :: ITH,IB,IKPT,ISPIN,IAT
      INTEGER(4)              :: IK1,IK2,IS1,IS2
      REAL(8)                 :: F0      !OCCUPATION       
      REAL(8)                 :: QION    !IONIZATION STATE 
      REAL(8)                 :: TOTSPIN ! SPIN STATE
      REAL(8)                 :: SUMOFZ  ! SUM OF NUCLEAR CHARGES FROM ALL ATOMS
      REAL(8)                 :: SVAR
      REAL(8)                 :: FMAX
      INTEGER(4)              :: NFILO
      INTEGER(4)              :: NEMPTY,NOCC
      INTEGER(4)              :: ISP
      LOGICAL(4)              :: TNONCOLL
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_OCCUP')
      LL_STRC=LL_STRC_
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!    
!     ==================================================================
!     ==  READ BLOCK !STRUCTURE!OCCUPATIONS                           ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'OCCUPATIONS')
!
!     ==================================================================
!     == NSPIN  #SPINPONENTS IN THE DENSITY (1,2 OR 4)                ==
!     ==   (1,2,4) FOR SPIN-RESTRICTED, COLLINEAR, NON-COLLINEAR      ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_STRC,'NSPIN',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'NSPIN',0,1)
      CALL LINKEDLIST$GET(LL_STRC,'NSPIN',1,NSPIN)
      IF(NSPIN.NE.1.AND.NSPIN.NE.2.AND.NSPIN.NE.3) THEN
        CALL ERROR$MSG('NUMBER OF SPINS OUT OF RANGE')
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('STRCIN_OCCUP')
      END IF 
      TNONCOLL=(NSPIN.EQ.3)
      FMAX=1.D0
      IF(NSPIN.EQ.1) FMAX=2.D0
!
!     ==================================================================
!     == TOTAL SPIN                                                   ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_STRC,'SPIN[HBAR]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'SPIN[HBAR]',0,0.D0)
      CALL LINKEDLIST$CONVERT(LL_STRC,'SPIN[HBAR]',1,'R(8)')
      CALL LINKEDLIST$GET(LL_STRC,'SPIN[HBAR]',1,TOTSPIN) 
      TOTSPIN=2.D0*TOTSPIN   ! CONVERT INTO #(UP-E)-#(DOWN-E)
      IF(TNONCOLL) THEN
        IF(TCHK) THEN
          CALL ERROR$MSG('OPTION SPIN[HBAR] IS INCOMPATIBLE WITH NSPIN=3')
          CALL ERROR$STOP('STRCIN_OCCUP')
        END IF
        TOTSPIN=0.D0
      ELSE
        IF(ABS(TOTSPIN).GT.1.D-6.AND.NSPIN.EQ.1) THEN
          CALL ERROR$MSG('NONZERO SPIN[HBAR] IS INCOMPATIBLE WITH NSPIN=1')
          CALL ERROR$STOP('STRCIN_OCCUP')
        END IF
      END IF
!
!     ==================================================================
!     == COLLECT TOTAL CHARGE                                         ==
!     ==================================================================
      CALL LINKEDLIST$EXISTD(LL_STRC,'CHARGE[E]',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'CHARGE[E]',0,0.D0)
      CALL LINKEDLIST$CONVERT(LL_STRC,'CHARGE[E]',1,'R(8)')
      CALL LINKEDLIST$GET(LL_STRC,'CHARGE[E]',1,QION)
      QION =-QION
!    
!     ==================================================================
!     ==  COLLECT #(BANDS)                                            ==
!     ==================================================================
      NB=0
      CALL LINKEDLIST$EXISTD(LL_STRC,'NBAND',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_STRC,'NBAND',1,NB)
!
      NEMPTY=1
      CALL LINKEDLIST$EXISTD(LL_STRC,'EMPTY',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_STRC,'EMPTY',1,NEMPTY)
!     
!     ==================================================================
!     ==  DETERMINE TOTAL NUCLEAR CHARGE                              ==
!     ==================================================================
      CALL ATOMLIST$NATOM(NAT)
      SUMOFZ=0.D0
      DO IAT=1,NAT
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('ZV',SVAR)
        CALL SETUP$UNSELECT()
        SUMOFZ=SUMOFZ+SVAR
      ENDDO
!    
!     ==================================================================
!     ==  DETERMINE ACTUAL #(BANDS)                                   ==
!     ==================================================================
!     == ESTIMATE FROM #(ELECTRONS) AND #(EMPTY BANDS) REQUESTED ========
      IF(NSPIN.EQ.1.OR.NSPIN.EQ.3) THEN ! NON-SPIN-POLARIZED
        NOCC=1+INT((SUMOFZ+QION)/FMAX-1.D-6)
        NB=MAX(NB,NOCC+NEMPTY)
      ELSE IF(NSPIN.EQ.2) THEN
  !      NOCC=0.5D0*REAL(NOCC+1,KIND=8)
  !      NOCC=NOCC+INT(0.5D0*ABS(TOTSPIN)+1.D0-1.D-5)
         NOCC=1+INT(0.5D0*(SUMOFZ+QION)/FMAX-1.D-6)
         NB=MAX(NB,NOCC+NEMPTY)
         NOCC=1+INT(0.5D0*(SUMOFZ+QION+ABS(TOTSPIN))/FMAX-1.D-6)
         NB=MAX(NB,NOCC)
      END IF
!
!     ==  FIND NUMBER OF THE HIGHEST SPECIFIED BAND  ====================
      CALL LINKEDLIST$NLISTS(LL_STRC,'STATE',NTH)
      DO ITH=1,NTH
        CALL LINKEDLIST$SELECT(LL_STRC,'STATE',ITH)
        CALL LINKEDLIST$EXISTD(LL_STRC,'B',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'B',1,IB)
          NB=MAX(NB,IB)
        END IF
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
!
!     ==================================================================
!     ==  INITIALIZE DYNOCC OBJECT                                    ==
!     ==================================================================
      IF(TNONCOLL) NSPIN=1
      NB=2*INT((NB+1)/2)    !ONLY EVEN NUMBER OF STATES ALLOWED
      CALL DYNOCC$SETR8('SUMOFZ',SUMOFZ)
      CALL DYNOCC$SETR8('TOTCHA',QION)
      CALL DYNOCC$SETR8('SPIN',TOTSPIN)
      CALL DYNOCC$SETR8('FMAX',FMAX)
      CALL DYNOCC$CREATE(NB,NKPT,NSPIN) 
      CALL DYNOCC$INIOCC      ! GUESS OCCUPATIONS FOR EACH STATE
!
      IF(TNONCOLL) THEN
        CALL WAVES$SETI4('SPINORDIM',2)
        CALL DFT$SETL4('SPIN',.TRUE.)
      ELSE
        CALL WAVES$SETI4('SPINORDIM',1)
        CALL DFT$SETL4('SPIN',NSPIN.EQ.2)
      END IF        
!    
!     ==================================================================
!     ==  NOW ADAPT THE OCCUPATIONS                                   ==
!     ==================================================================
      CALL LINKEDLIST$NLISTS(LL_STRC,'STATE',NTH)
      DO ITH=1,NTH
        CALL LINKEDLIST$SELECT(LL_STRC,'STATE',ITH)
!       == COLLECT K-POINTS ============================================
        IKPT=0
        CALL LINKEDLIST$EXISTD(LL_STRC,'K',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'K',1,IKPT)
          IF(IKPT.GT.NKPT) THEN
            CALL ERROR$OVERFLOW('IKPT',IKPT,NKPT)
            CALL ERROR$STOP('STRCIN; !STRUCTURE!OCCUPATIONS!STATE:K')
          ENDIF
        END IF
!       == COLLECT SPIN ================================================
        ISPIN=0
        CALL LINKEDLIST$EXISTD(LL_STRC,'S',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'S',1,ISPIN)
          IF(ISPIN.NE.1.AND.ISPIN.NE.NSPIN) THEN
            CALL ERROR$MSG('!STRC!OCCUPATIONS!STATE:S OUT OF RANGE')
            CALL ERROR$MSG('ALLOWED VALUES ARE 1 AND 2')
            CALL ERROR$I4VAL('S',ISPIN)
            CALL ERROR$I4VAL('NSPIN',NSPIN)
            CALL ERROR$STOP('STRCIN; !STRUCTURE!OCCUPATIONS!STATE:S')
          ENDIF
        END IF
!       == COLLECT BAND INDEX ==========================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'B',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('BAND INDEX "B" HAS NOT BEEN SPECIFIED')
          CALL ERROR$I4VAL('STATE BLOCK NR.',ITH)
          CALL ERROR$STOP('STRCIN; !STRUCTURE!OCCUPATIONS!STATE:B')
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'B',1,IB)
!       == COLLECT OCCUPATION ==========================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'F',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('OCCUPATION HAS NOT BEEN SPECIFIED')
          CALL ERROR$I4VAL('STATE BLOCK NR.',ITH)
          CALL ERROR$STOP('STRCIN; !STRUCTURE!OCCUPATIONS!STATE:F')
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'F',1,F0)
        CALL DYNOCC$GETR8('FMAX',FMAX)
        IF(F0.LT.0.D0.OR.F0.GT.FMAX) THEN
          CALL ERROR$MSG('OCCUPATION OUT OF RANGE')
          CALL ERROR$I4VAL('STATE BLOCK NR.',ITH)
          CALL ERROR$R8VAL('SPECIFIED F',F0)
          CALL ERROR$R8VAL('MAX. F',FMAX)
          CALL ERROR$STOP('STRCIN; !STRUCTURE!OCCUPATIONS!STATE:F')
        END IF
!       == SELECT SPIN RANGE ===========================================
        IF(ISPIN.EQ.0) THEN
          IS1=1
          IS2=NSPIN
        ELSE
          IS1=ISPIN
          IS2=ISPIN
        END IF
!       == SELECT K-POINT RANGE ========================================
        IF(IKPT.EQ.0) THEN
          IK1=1
          IK2=NKPT
        ELSE
          IK1=IKPT
          IK2=IKPT
        END IF
!       == NOW SPECIFY OCCUPATIONS IN DYNOCC OBJECT ====================
        DO ISPIN=IS1,IS2
          DO IKPT=IK1,IK2
            CALL DYNOCC$MODOCC(IB,IKPT,ISPIN,F0)
          ENDDO
        ENDDO
!
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_OCCUP
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_ORBPOT(LL_STRC_)
!     **************************************************************************
!     **  DEFINE GROUPS IN GROUPLIST                                          **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE    ! +/-  MAKES STRING UPPERCASE/LOWERCASE
      USE LINKEDLIST_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK
      LOGICAL(4)               :: TON
      CHARACTER(32)            :: ATOM,NAME
      INTEGER(4)               :: NTH,ITH
      CHARACTER(32)            :: TYPE
      INTEGER(4)               :: ISPIN
      INTEGER(4)               :: NAT
      INTEGER(4)               :: I
      INTEGER(4)               :: IAT
      REAL(8)                  :: VALUE
      REAL(8)                  :: RC
      REAL(8)                  :: PWR
      REAL(8)                  :: AEZ
      REAL(8)                  :: RCOV
      INTEGER(4)               :: LM
!     **************************************************************************
                           CALL TRACE$PUSH('STRCIN_ORBPOT')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'ORBPOT',1,TON)
      IF(.NOT.TON) THEN
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==========================================================================
!     ==  COLLECT NSPIN FROM !STRUCTURE!OCCUPATIONS                           ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'ORBPOT')
      CALL LINKEDLIST$NLISTS(LL_STRC,'POT',NTH)
      DO ITH=1,NTH
        CALL LINKEDLIST$SELECT(LL_STRC,'POT',ITH)
!       == SPECIFY ATOM ========================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!ORBPOT!POT:ATOM NOT SPECIFIED')
          CALL ERROR$STOP('STRCIN_ORBPOT')
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'ATOM',1,ATOM)
!       == CHECK IF ATOM EXISTS ================================================
        CALL ATOMLIST$NATOM(NAT)
        TCHK=.FALSE.
        DO I=1,NAT
          IAT=I
          CALL ATOMLIST$GETCH('NAME',IAT,NAME)
          TCHK=TCHK.OR.(NAME.EQ.ATOM)
          IF(TCHK) EXIT
        ENDDO
        IF(.NOT.TCHK) THEN 
          CALL ERROR$MSG('!STRUCTURE!ORBPOT!POT:ATOM VALUE NOT RECOGNIZED')
          CALL ERROR$MSG('ATOM NOT IN ATOM LIST')
          CALL ERROR$CHVAL('ATOM',ATOM)
          CALL ERROR$MSG('SUGGESTION: CHECK UPPER-LOWER CASE')
          CALL ERROR$STOP('STRCIN_ORBPOT')
        END IF
        CALL ATOMLIST$GETR8('Z',IAT,AEZ)
        CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
!
!       == SPECIFY TYPE ========================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'TYPE',1,TCHK)
        IF(.NOT.TCHK) THEN
          TYPE='ALL'
        ELSE
          CALL LINKEDLIST$GET(LL_STRC,'TYPE',1,TYPE)
          TYPE=+TYPE   ! MAKE INPUT CASE INSENSITIVE BY UPPERCASING TYPE
        END IF
        TCHK=.FALSE.
        IF(TYPE.EQ.'P')   TCHK=.TRUE.
        IF(TYPE.EQ.'D')   TCHK=.TRUE.
        IF(TYPE.EQ.'F')   TCHK=.TRUE.
        IF(TYPE.EQ.'ALL') TCHK=.TRUE.
!       == PURE ANGULAR MOMENTA ================================================
!       == SEE SPHERICAL$LMBYNAME ==============================================
        IF(.NOT.TCHK) THEN
!         == SPHERICAL$LMBYNAME THROWS AN ERROR IF TYPE IS NOT RECOGNIZED ======
          CALL SPHERICAL$LMBYNAME(TYPE,LM)
        END IF
!       == SPECIFY SPIN DIRECTION ==============================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'S',1,TCHK)
        IF(.NOT.TCHK) THEN
          ISPIN=0          ! ISPIN=0 IS DEFAULT FOR ALL SPIN DIRECTIONS
        ELSE
          CALL LINKEDLIST$GET(LL_STRC,'S',1,ISPIN)
          IF(ISPIN.LT.1.OR.ISPIN.GT.3) THEN
            CALL ERROR$MSG('VALUE  OF !STRUCTURE!ORBPOT!POT:S NOT ALLOWED')
            CALL ERROR$I4VAL('ISPIN',ISPIN)
            CALL ERROR$STOP('STRCIN_ORBPOT')
          END IF
        END IF
!       == SPECIFY CUTOFF RADIUS ===============================================
        RC=RCOV
        CALL LINKEDLIST$EXISTD(LL_STRC,'RC',1,TCHK)
!!$        IF(.NOT.TCHK) THEN
!!$          CALL ERROR$MSG('!STRUCTURE!ORBPOT!POT:RC NOT SPECIFIED')
!!$          CALL ERROR$STOP('STRCIN_ORBPOT')
!!$        END IF
!!$        CALL LINKEDLIST$GET(LL_STRC,'RC',1,RC)
        IF(TCHK)CALL LINKEDLIST$GET(LL_STRC,'RC',1,RC)
!       == SPECIFY POWER =======================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'PWR',1,TCHK)
        IF(.NOT.TCHK) THEN
          PWR=2.D0
        ELSE
          CALL LINKEDLIST$GET(LL_STRC,'PWR',1,PWR)
        END IF
!       == SPECIFY VALUE =======================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'VALUE',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!ORBPOT!POT:VALUE NOT SPECIFIED')
          CALL ERROR$STOP('STRCIN_ORBPOT')
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'VALUE',1,VALUE)
!       ========================================================================
        CALL EXTERNAL1CPOT$SETPOT(ATOM,TYPE,ISPIN,VALUE,RC,PWR)
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_ORBPOT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_GROUP(LL_STRC_)
!     ******************************************************************
!     **  DEFINE GROUPS IN GROUPLIST                                  **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: NGRP
      INTEGER(4)               :: NAT
      INTEGER(4)               :: NSP
      INTEGER(4)               :: NATGR1
      INTEGER(4)               :: IAT,IGRP,IAT1,ISP
      CHARACTER(32)            :: GROUP
      CHARACTER(32)            :: NAME
      CHARACTER(32),ALLOCATABLE:: ATOMS(:)
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_GROUP')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$NLISTS(LL_STRC,'GROUP',NGRP)
      CALL ATOMLIST$NATOM(NAT)
      CALL SETUP$GETI4('NSP',NSP)
      CALL GROUPLIST$INITIALIZE(NAT,NGRP+1+NSP)
!
!     ==================================================================
!     ==  DEFINE DEFAULT GROUPS                                       ==
!     ==  1) GROUP ALL CONTAINING ALL ATOMS                           ==
!     ==  2) ONE GROUP FOR EACH SPECIES WITH THE SPECIES NAME         ==
!     ==     AS GROUP NAME                                            ==
!     ==================================================================
      DO IAT=1,NAT
        CALL GROUPLIST$ADD('ALL',IAT)
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETCH('ID',NAME)
        CALL SETUP$UNSELECT
        CALL GROUPLIST$ADD(NAME,IAT)
      ENDDO
!
!     ==  READ ACTUAL VALUES  ==========================================
      DO IGRP=1,NGRP
        CALL LINKEDLIST$SELECT(LL_STRC,'GROUP',IGRP)
        CALL LINKEDLIST$EXISTD(LL_STRC,'NAME',0,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRC!GROUP:NAME NOT FOUND')
          CALL ERROR$STOP('STRCIN_GROUP')
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'NAME',1,GROUP)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOMS',0,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRC!GROUP:ATOMS NOT FOUND')
          CALL ERROR$STOP('STRCIN_GROUP')
        END IF
        CALL LINKEDLIST$SIZE(LL_STRC,'ATOMS',1,NATGR1)
        ALLOCATE(ATOMS(NATGR1))
        CALL LINKEDLIST$GET(LL_STRC,'ATOMS',1,ATOMS(:))
        DO IAT=1,NATGR1        
          CALL ATOMLIST$INDEX(ATOMS(IAT),IAT1)
          CALL GROUPLIST$ADD(GROUP,IAT1)
        ENDDO
        DEALLOCATE(ATOMS)
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_GROUP
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_VEXT(LL_STRC_)
!     ******************************************************************
!     **  DEFINE GROUPS IN GROUPLIST                                  **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK,TCHK1,TCHK2
      INTEGER(4)               :: NAT
      LOGICAL(4)   ,ALLOCATABLE:: TIAT1(:)
      LOGICAL(4)   ,ALLOCATABLE:: TIAT2(:)
      CHARACTER(32)            :: NAME
      REAL(8)                  :: EV
      REAL(8)                  :: E0
      INTEGER(4)               :: I,N
      INTEGER(4)               :: IAT1,IAT2,IAT
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_VEXT')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'VEXT',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL TRACE$POP
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_STRC,'VEXT')
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(TIAT1(NAT))
      ALLOCATE(TIAT2(NAT))
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL LINKEDLIST$NLISTS(LL_STRC,'UNBIND',N)
      DO I=1,N
        CALL LINKEDLIST$SELECT(LL_STRC,'UNBIND',I)
!       == FIND FIRST GROUP OR ATOM
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM1',0,TCHK1)
        CALL LINKEDLIST$EXISTD(LL_STRC,'GROUP1',0,TCHK2)
        IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_STRC,'ATOM1',1,NAME)
          CALL ATOMLIST$INDEX(NAME,IAT)
          TIAT1(:)=.FALSE.
          TIAT1(IAT)=.TRUE.
        ELSE IF(TCHK2) THEN
          CALL LINKEDLIST$GET(LL_STRC,'GROUP1',1,NAME)
          CALL GROUPLIST$MEMBERS(NAME,NAT,TIAT1) 
        ELSE
          CALL ERROR$MSG('!STRC!GROUP:ATOM1 NOT FOUND')
          CALL ERROR$MSG('!STRC!GROUP:GROUP1 NOT FOUND')
          CALL ERROR$STOP('STRCIN_VEXT')
        END IF
!       == FIND SECOND GROUP OR ATOM
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM2',0,TCHK1)
        CALL LINKEDLIST$EXISTD(LL_STRC,'GROUP2',0,TCHK2)
        IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_STRC,'ATOM2',1,NAME)
          CALL ATOMLIST$INDEX(NAME,IAT)
          TIAT2(:)=.FALSE.
          TIAT2(IAT)=.TRUE.
        ELSE IF(TCHK2) THEN
          CALL LINKEDLIST$GET(LL_STRC,'GROUP2',1,NAME)
          CALL GROUPLIST$MEMBERS(NAME,NAT,TIAT2) 
        ELSE
          CALL ERROR$MSG('!STRC!GROUP:ATOM1 NOT FOUND')
          CALL ERROR$MSG('!STRC!GROUP:GROUP1 NOT FOUND')
          CALL ERROR$STOP('STRCIN_VEXT')
        END IF
!       ================================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'E[EV]',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL LINKEDLIST$SET(LL_STRC,'E[EV]',0,1.D0)
        END IF
        CALL CONSTANTS$GET('EV',EV)
        CALL LINKEDLIST$GET(LL_STRC,'E[EV]',1,E0)
        E0=E0*EV
!       ================================================================        
        DO IAT1=1,NAT
          DO IAT2=IAT1,NAT
            TCHK1=TIAT1(IAT1).AND.TIAT2(IAT2)
            TCHK2=TIAT2(IAT1).AND.TIAT1(IAT2)
            IF(TCHK1.OR.TCHK2) THEN
              CALL VEXT$SETUNBIND(IAT1,IAT2,E0)
            END IF
          ENDDO
        ENDDO
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      END DO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_VEXT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_LINEAR(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)               :: NUM,NUMJ
      INTEGER(4)               :: ITH,JTH
      INTEGER(4)               :: NAT
      INTEGER(4)               :: IAT
      CHARACTER(32)            :: CH32SVAR1
      REAL(8)     ,ALLOCATABLE :: VEC1(:,:) !(3,NAT)
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_LINEAR')
      LL_STRC=LL_STRC_
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(VEC1(3,NAT))

      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')

      CALL LINKEDLIST$NLISTS(LL_STRC,'LINEAR',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'LINEAR',ITH)
        VEC1(:,:)=0.D0
        CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NUMJ)
        IF(NUMJ.EQ.0) THEN
          CALL ERROR$MSG('BLOCK !ATOM IS MANDATORY')
          CALL ERROR$STOP('STRCIN; !STRUCTURE!CONSTRAINTS!LINEAR')
        END IF
!        
        DO JTH=1,NUMJ
          CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',JTH)
          CALL LINKEDLIST$GET(LL_STRC,'NAME',1,CH32SVAR1)
          CALL ATOMLIST$INDEX(CH32SVAR1,IAT)
          CALL LINKEDLIST$GET(LL_STRC,'R',1,VEC1(:,IAT))
          CALL LINKEDLIST$SELECT(LL_STRC,'..')
        ENDDO
!
        CALL CONSTRAINTS$OPEN('GENERALLINEAR','LINEAR CONSTRAINT')
        CALL CONSTRAINTS$DEFINER8A('VEC',3*NAT,VEC1)
        CALL READMOVABLECONSTRAINT(LL_STRC)
        CALL CONSTRAINTS$CLOSE
!
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
      DEALLOCATE(VEC1)
                             CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_LINEAR
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_BOND(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)               :: NUM
      INTEGER(4)               :: ITH
      CHARACTER(32)            :: CH32SVAR1
      CHARACTER(32)            :: CH32SVAR2
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_BOND')
!
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')

      CALL LINKEDLIST$NLISTS(LL_STRC,'BOND',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'BOND',ITH)
        CALL LINKEDLIST$GET(LL_STRC,'ATOM1',1,CH32SVAR1)
        CALL LINKEDLIST$GET(LL_STRC,'ATOM2',1,CH32SVAR2)
        CALL CONSTRAINTS$OPEN('BONDLENGTH' &
     &                  ,'BOND '//CH32SVAR1(1:10)//'-'//CH32SVAR2(1:10))
        CALL CONSTRAINTS$DEFINECH('ATOM1',CH32SVAR1(1:32))
        CALL CONSTRAINTS$DEFINECH('ATOM2',CH32SVAR2(1:32))
        CALL READMOVABLECONSTRAINT(LL_STRC)
        CALL CONSTRAINTS$CLOSE
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_BOND
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_RIGID(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)            :: NUM
      INTEGER(4)            :: ITH
      CHARACTER(32)         :: CH32SVAR1
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_RIGID')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')

      CALL LINKEDLIST$NLISTS(LL_STRC,'RIGID',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'RIGID',ITH)
        CALL LINKEDLIST$GET(LL_STRC,'GROUP',1,CH32SVAR1)
!
!       == PERFORM ACTIONS==============================================
        CALL CONSTRAINTS$OPEN('RIGID' &
     &                       ,'GROUP '//CH32SVAR1//' IS KEPT RIGID')
        CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
        CALL CONSTRAINTS$CLOSE
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO 
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_RIGID
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_FREEZE(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)            :: NUM
      INTEGER(4)            :: ITH
      CHARACTER(32)         :: CH32SVAR1
      LOGICAL(4)            :: TCHK,TCHK1
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_FREEZE')
      LL_STRC=LL_STRC_
!
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')
      CALL LINKEDLIST$NLISTS(LL_STRC,'FREEZE',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'FREEZE',ITH)
        CALL LINKEDLIST$EXISTD(LL_STRC,'GROUP',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'GROUP',1,CH32SVAR1)
          CALL CONSTRAINTS$OPEN('FIXGROUP','FIX GROUP:'//CH32SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$CLOSE
        END IF
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM',1,TCHK1)
        IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_STRC,'ATOM',1,CH32SVAR1)
          CALL CONSTRAINTS$OPEN('FIXATOM','FIX ATOM:'//CH32SVAR1)
          CALL CONSTRAINTS$DEFINECH('ATOM',CH32SVAR1)
          CALL CONSTRAINTS$CLOSE
        END IF
        IF(TCHK.EQV.TCHK1) THEN
          CALL ERROR$MSG('EXACTLY ONE OF THE KEYWORDS ATOM= OR GROUP=')
          CALL ERROR$MSG('NEED TO BE SPECIFIED')
          IF(TCHK1) THEN
            CALL ERROR$CHVAL('SPECIFIED','ATOM')
          END IF
          IF(TCHK) THEN
            CALL ERROR$CHVAL('SPECIFIED','GROUP')
          END IF
          CALL ERROR$STOP('STRCIN!STRUCTURE!CONSTRAINTS!FREEZE')
        END IF
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP

      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_FREEZE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_MIDPLANE(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)            :: NUM
      INTEGER(4)            :: ITH
      CHARACTER(32)    :: CH32SVAR1
      CHARACTER(32)    :: CH32SVAR2
      CHARACTER(32)    :: CH32SVAR3
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_MIDPLANE')
      LL_STRC=LL_STRC_
!
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')

      CALL LINKEDLIST$NLISTS(LL_STRC,'MIDPLANE',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'MIDPLANE',ITH)
        CALL LINKEDLIST$GET(LL_STRC,'ATOM1',1,CH32SVAR1)
        CALL LINKEDLIST$GET(LL_STRC,'ATOM2',1,CH32SVAR2)
        CALL LINKEDLIST$GET(LL_STRC,'ATOM3',1,CH32SVAR3)
        CALL CONSTRAINTS$OPEN('MIDPLANE' &
     &                    ,'MIDPLANE:'//CH32SVAR1(1:10)// &
     &                     ','//CH32SVAR2(1:10)//','//CH32SVAR3(1:10))
        CALL CONSTRAINTS$DEFINECH('ATOM1',CH32SVAR1)
        CALL CONSTRAINTS$DEFINECH('ATOM2',CH32SVAR2)
        CALL CONSTRAINTS$DEFINECH('ATOM3',CH32SVAR3)
        CALL READMOVABLECONSTRAINT(LL_STRC)
        CALL CONSTRAINTS$CLOSE
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_MIDPLANE

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_TRANSLATION(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)               :: NUM
      INTEGER(4)               :: ITH
      CHARACTER(32)            :: CH32SVAR1
      CHARACTER(256)           :: CH256SVAR1
      REAL(8)                  :: DIR(3) 
      LOGICAL(4)               :: TCHK
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_TRANSLATION')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')

      CALL LINKEDLIST$NLISTS(LL_STRC,'TRANSLATION',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'TRANSLATION',ITH)
        
        CALL LINKEDLIST$EXISTD(LL_STRC,'GROUP',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'GROUP',0,'ALL')
        CALL LINKEDLIST$GET(LL_STRC,'GROUP',1,CH32SVAR1)

        CALL LINKEDLIST$EXISTD(LL_STRC,'DIR',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'DIR',1,DIR)
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("CENTER OF GRAVITY OF GROUP ",A10' &
     &          //'," ALONG:",3F6.2," IS FIXED")') CH32SVAR1,DIR(:)
          CALL CONSTRAINTS$OPEN('TRANSLATIONALMOMENTUM',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('DIR',3,DIR)
          CALL CONSTRAINTS$CLOSE
        ELSE
          DIR(1)=1.D0   
          DIR(2)=0.D0   
          DIR(3)=0.D0   
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("CENTER OF GRAVITY OF GROUP ",A10' &
     &            //'," ALONG:",3F6.2," IS FIXED")') CH32SVAR1,DIR(:)
          CALL CONSTRAINTS$OPEN('TRANSLATIONALMOMENTUM',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('DIR',3,DIR)
          CALL CONSTRAINTS$CLOSE
          DIR(1)=0.D0   
          DIR(2)=1.D0   
          DIR(3)=0.D0   
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("CENTER OF GRAVITY OF GROUP ",A10' &
     &            //'," ALONG:",3F6.2," IS FIXED")') CH32SVAR1,DIR(:)
          CALL CONSTRAINTS$OPEN('TRANSLATIONALMOMENTUM',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('DIR',3,DIR)
          CALL CONSTRAINTS$CLOSE
          DIR(1)=0.D0   
          DIR(2)=0.D0   
          DIR(3)=1.D0   
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("CENTER OF GRAVITY OF GROUP ",A10' &
     &         //'," ALONG:",3F6.2," IS FIXED")') &
     &          CH32SVAR1,DIR(:)
          CALL CONSTRAINTS$OPEN('TRANSLATIONALMOMENTUM',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('DIR',3,DIR)
          CALL CONSTRAINTS$CLOSE
        END IF
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_TRANSLATION
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_ORIENTATION(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)               :: NUM
      INTEGER(4)               :: ITH
      CHARACTER(32)            :: CH32SVAR1
      CHARACTER(256)           :: CH256SVAR1
      REAL(8)                  :: AXIS(3) 
      LOGICAL                  :: TCHK
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_ORIENTATION')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')

      CALL LINKEDLIST$NLISTS(LL_STRC,'ORIENTATION',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'ORIENTATION',ITH)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'GROUP',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'GROUP',0,'ALL')
        CALL LINKEDLIST$GET(LL_STRC,'GROUP',1,CH32SVAR1)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'AXIS',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'AXIS',1,AXIS)
          CH256SVAR1=' '
            WRITE(CH256SVAR1,FMT='("ORIENTATION OF GROUP "' &
     &           //',A10," ABOUT:",3F6.2," IS FIXED")')CH32SVAR1,AXIS(:)
          CALL CONSTRAINTS$OPEN('ORIENTATION',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('PHI',3,AXIS)
          CALL CONSTRAINTS$CLOSE
        ELSE
          AXIS(1)=1.D0   
          AXIS(2)=0.D0   
          AXIS(3)=0.D0   
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("ORIENTATION OF GROUP "' &
     &         //',A10," ABOUT:",3F6.2," IS FIXED")')CH32SVAR1,AXIS(:)
          CALL CONSTRAINTS$OPEN('ORIENTATION',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('PHI',3,AXIS)
          CALL CONSTRAINTS$CLOSE
          AXIS(1)=0.D0   
          AXIS(2)=1.D0   
          AXIS(3)=0.D0   
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("ORIENTATION OF GROUP "' &
     &         //',A10," ABOUT:",3F6.2," IS FIXED")')CH32SVAR1,AXIS(:)
          CALL CONSTRAINTS$OPEN('ORIENTATION',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('PHI',3,AXIS)
          CALL CONSTRAINTS$CLOSE
          AXIS(1)=0.D0   
          AXIS(2)=0.D0   
          AXIS(3)=1.D0   
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("ORIENTATION OF GROUP "' &
     &         //',A10," ABOUT:",3F6.2," IS FIXED")')CH32SVAR1,AXIS(:)
          CALL CONSTRAINTS$OPEN('ORIENTATION',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('PHI',3,AXIS)
          CALL CONSTRAINTS$CLOSE
        END IF
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_ORIENTATION
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_ROTATION(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)               :: NUM
      INTEGER(4)               :: ITH
      CHARACTER(32)            :: CH32SVAR1
      CHARACTER(256)           :: CH256SVAR1
      REAL(8)                  :: AXIS(3) 
      LOGICAL                  :: TCHK
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_ROTATION')

      LL_STRC=LL_STRC_

      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')

      CALL LINKEDLIST$NLISTS(LL_STRC,'ROTATION',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'ROTATION',ITH)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'GROUP',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'GROUP',0,'ALL')
        CALL LINKEDLIST$GET(LL_STRC,'GROUP',1,CH32SVAR1)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'AXIS',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'AXIS',1,AXIS)
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("ROTATION OF GROUP ",A10' &
     &          //'," ABOUT:",3F6.2," IS FIXED")') &
     &          CH32SVAR1,AXIS(:)
          CALL CONSTRAINTS$OPEN('ANGULARMOMENTUM',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('PHI',3,AXIS)
          CALL CONSTRAINTS$CLOSE
        ELSE
          AXIS(1)=1.D0   
          AXIS(2)=0.D0   
          AXIS(3)=0.D0   
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("ROTATION OF GROUP ",A10' &
     &          //'," ABOUT:",3F6.2," IS FIXED")') &
     &          CH32SVAR1,AXIS(:)
          CALL CONSTRAINTS$OPEN('ANGULARMOMENTUM',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('PHI',3,AXIS)
          CALL CONSTRAINTS$CLOSE
          AXIS(1)=0.D0   
          AXIS(2)=1.D0   
          AXIS(3)=0.D0   
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("ROTATION OF GROUP ",A10' &
     &        //'," ABOUT:",3F6.2," IS FIXED")') &
     &        CH32SVAR1,AXIS(:)
          CALL CONSTRAINTS$OPEN('ANGULARMOMENTUM',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('PHI',3,AXIS)
          CALL CONSTRAINTS$CLOSE
          AXIS(1)=0.D0   
          AXIS(2)=0.D0   
          AXIS(3)=1.D0   
          CH256SVAR1=' '
          WRITE(CH256SVAR1,FMT='("ROTATION OF GROUP ",A10' &
     &          //'," ABOUT:",3F6.2," IS FIXED")') &
     &          CH32SVAR1,AXIS(:)
          CALL CONSTRAINTS$OPEN('ANGULARMOMENTUM',CH256SVAR1)
          CALL CONSTRAINTS$DEFINECH('GROUP',CH32SVAR1)
          CALL CONSTRAINTS$DEFINER8A('PHI',3,AXIS)
          CALL CONSTRAINTS$CLOSE
        END IF
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_ROTATION
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_COGSEP(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)               :: NUM
      INTEGER(4)               :: ITH
      CHARACTER(32)            :: NAME1,NAME2
      LOGICAL                  :: TCHK
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_ORIENTATION')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')

      CALL LINKEDLIST$NLISTS(LL_STRC,'COGSEP',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'COGSEP',ITH)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'GROUP1',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!CONSTRAINTS!COGSEP:GROUP1 NOT FOUND')
          CALL ERROR$STOP('STRCIN_CONSTRAINTS_COGSEP')
        ENDIF
        CALL LINKEDLIST$GET(LL_STRC,'GROUP1',1,NAME1)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'GROUP2',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!CONSTRAINTS!COGSEP:GROUP2 NOT FOUND')
          CALL ERROR$STOP('STRCIN_CONSTRAINTS_COGSEP')
        ENDIF
        CALL LINKEDLIST$GET(LL_STRC,'GROUP2',1,NAME2)
!
        CALL CONSTRAINTS$OPEN('COGSEP' &
     &                  ,'COGSEP '//NAME1(1:10)//'-'//NAME2(1:10))
        CALL CONSTRAINTS$DEFINECH('GROUP1',NAME1)
        CALL CONSTRAINTS$DEFINECH('GROUP2',NAME2)
        CALL READMOVABLECONSTRAINT(LL_STRC)
        CALL CONSTRAINTS$CLOSE
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_COGSEP
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_ANGLE(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)               :: NUM
      INTEGER(4)               :: ITH
      CHARACTER(32)            :: NAME1,NAME2,NAME3
      LOGICAL                  :: TCHK
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_ORIENTATION')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')

      CALL LINKEDLIST$NLISTS(LL_STRC,'ANGLE',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'ANGLE',ITH)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM1',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!CONSTRAINTS!ANGLE:ATOM1 NOT FOUND')
          CALL ERROR$STOP('STRCIN_CONSTRAINTS_ANGLE')
        ENDIF
        CALL LINKEDLIST$GET(LL_STRC,'ATOM1',1,NAME1)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM2',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!CONSTRAINTS!ANGLE:ATOM2 NOT FOUND')
          CALL ERROR$STOP('STRCIN_CONSTRAINTS_ANGLE')
        ENDIF
        CALL LINKEDLIST$GET(LL_STRC,'ATOM2',1,NAME2)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM3',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!CONSTRAINTS!ANGLE:ATOM3 NOT FOUND')
          CALL ERROR$STOP('STRCIN_CONSTRAINTS_ANGLE')
        ENDIF
        CALL LINKEDLIST$GET(LL_STRC,'ATOM3',1,NAME3)
!
        CALL CONSTRAINTS$OPEN('BONDANGLE' &
     &                  ,'ANGLE '//NAME1(1:10)//'-'//NAME2(1:10)//'-'//NAME3(1:10))
        CALL CONSTRAINTS$DEFINECH('ATOM1',NAME1)
        CALL CONSTRAINTS$DEFINECH('ATOM2',NAME2)
        CALL CONSTRAINTS$DEFINECH('ATOM3',NAME3)
        CALL READMOVABLECONSTRAINT(LL_STRC)
        CALL CONSTRAINTS$CLOSE
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_ANGLE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_CONSTRAINTS_TORSION(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      INTEGER(4)               :: NUM
      INTEGER(4)               :: ITH
      CHARACTER(32)            :: NAME1,NAME2,NAME3,NAME4
      LOGICAL                  :: TCHK
!     ******************************************************************
                           CALL TRACE$PUSH('STRCIN_CONSTRAINTS_ORIENTATION')
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'CONSTRAINTS')

      CALL LINKEDLIST$NLISTS(LL_STRC,'TORSION',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_STRC,'TORSION',ITH)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM1',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!CONSTRAINTS!TORSION:ATOM1 NOT FOUND')
          CALL ERROR$STOP('STRCIN_CONSTRAINTS_TORSION')
        ENDIF
        CALL LINKEDLIST$GET(LL_STRC,'ATOM1',1,NAME1)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM2',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!CONSTRAINTS!TORSION:ATOM2 NOT FOUND')
          CALL ERROR$STOP('STRCIN_CONSTRAINTS_TORSION')
        ENDIF
        CALL LINKEDLIST$GET(LL_STRC,'ATOM2',1,NAME2)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM3',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!CONSTRAINTS!TORSION:ATOM3 NOT FOUND')
          CALL ERROR$STOP('STRCIN_CONSTRAINTS_TORSION')
        ENDIF
        CALL LINKEDLIST$GET(LL_STRC,'ATOM3',1,NAME3)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'ATOM4',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!CONSTRAINTS!TORSION:ATOM4 NOT FOUND')
          CALL ERROR$STOP('STRCIN_CONSTRAINTS_TORSION')
        ENDIF
        CALL LINKEDLIST$GET(LL_STRC,'ATOM4',1,NAME4)
!
        CALL CONSTRAINTS$OPEN('TORSION' &
     &            ,'TORSION '//NAME1(1:10)//'-'//NAME2(1:10) &
     &                       //'-'//NAME3(1:10)//'-'//NAME4(1:10))
        CALL CONSTRAINTS$DEFINECH('ATOM1',NAME1)
        CALL CONSTRAINTS$DEFINECH('ATOM2',NAME2)
        CALL CONSTRAINTS$DEFINECH('ATOM3',NAME3)
        CALL CONSTRAINTS$DEFINECH('ATOM4',NAME4)
        CALL READMOVABLECONSTRAINT(LL_STRC)
        CALL CONSTRAINTS$CLOSE
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
                           CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCIN_CONSTRAINTS_TORSION
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READMOVABLECONSTRAINT(LL_STRC_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK
      REAL(8)                  :: SVAR
      REAL(8)                  :: SMASS
      REAL(8)                  :: ANNES1
      REAL(8)                  :: PI
      INTEGER(4)               :: ISVAR
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      LL_STRC=LL_STRC_
                       CALL TRACE$PUSH('READMOVABLECONSTRAINT')
!     
!     == DECIDE WHETHER PRINTOUT IS REQUESTED ========================
      CALL LINKEDLIST$EXISTD(LL_STRC,'SHOW',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'SHOW',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_STRC,'SHOW',1,TCHK)
      CALL CONSTRAINTS$DEFINEL4('SHOW',TCHK)
!        
!       == REQUEST WHETHER CONSTRAINT SHALL FLOAT ====================
      CALL LINKEDLIST$EXISTD(LL_STRC,'FLOAT',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'FLOAT',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_STRC,'FLOAT',1,TCHK)
      CALL CONSTRAINTS$DEFINEL4('FLOAT',TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$EXISTD(LL_STRC,'M',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'M',0,1.D+4)
        CALL LINKEDLIST$GET(LL_STRC,'M',1,SMASS)
        CALL CONSTRAINTS$DEFINER8('SMASS',SMASS)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'FRIC',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'FRIC',0,0.D0)
        CALL LINKEDLIST$GET(LL_STRC,'FRIC',1,ANNES1)
        CALL CONSTRAINTS$DEFINER8('ANNES',ANNES1)
      END IF
!        
!       == REQUEST WHETHER CONSTRAINT BE MOVED  ======================
      CALL LINKEDLIST$EXISTD(LL_STRC,'MOVE',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'MOVE',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_STRC,'MOVE',1,TCHK)
      CALL CONSTRAINTS$DEFINEL4('MOVE',TCHK)
!
      IF(TCHK) THEN 
        CALL LINKEDLIST$EXISTD(LL_STRC,'NSTEP',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_STRC,'NSTEP',0,100)
        CALL LINKEDLIST$GET(LL_STRC,'NSTEP',1,ISVAR)
        CALL CONSTRAINTS$DEFINEI4('NSTEP',ISVAR)
!
        CALL LINKEDLIST$EXISTD(LL_STRC,'VALUE[DEG]',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'VALUE[DEG]',1,SVAR)
          CALL CONSTRAINTS$DEFINER8('SFINAL',SVAR/180.D0*PI)
        ELSE
          CALL LINKEDLIST$EXISTD(LL_STRC,'VALUE',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'VALUE',1,SVAR)
            CALL CONSTRAINTS$DEFINER8('SFINAL',SVAR)
          END IF
        ENDIF
        CALL LINKEDLIST$EXISTD(LL_STRC,'VELOC',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'VELOC',1,SVAR)
          CALL CONSTRAINTS$DEFINER8('VFINAL',SVAR)
        END IF
      END IF
                       CALL TRACE$POP
      RETURN
    END SUBROUTINE READMOVABLECONSTRAINT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_SOLVENT(LL_STRC_)
!     ******************************************************************    
!     **                                                              **    
!     ******************************************************************    
      USE LINKEDLIST_MODULE
      USE FORCEFIELD_MODULE
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE ATOM_TYPE
        CHARACTER(32)     :: NAME
        REAL(8)           :: R(3)
        REAL(8)           :: Q
        REAL(8)           :: M
        CHARACTER(5)      :: FFTYPE
        CHARACTER(2)      :: ELEMENT
        INTEGER(4)        :: QMSATOM   ! M-ATOM FOR Q ATOM
                                       ! S ATOM FOR M-ATOM
                                       ! Q ATOM FOR S ATOM
      END TYPE ATOM_TYPE
      TYPE LINK_TYPE
        INTEGER(4)        :: MJOINT
        INTEGER(4)        :: QJOINT
        INTEGER(4)        :: SJOINT
        INTEGER(4)        :: MATOM
        INTEGER(4)        :: QATOM
        INTEGER(4)        :: SATOM
      END TYPE LINK_TYPE
      TYPE BOND_TYPE
        INTEGER(4)        :: ATOM1
        INTEGER(4)        :: ATOM2
        REAL(8)           :: BO
        INTEGER(4)        :: IT(3)
      END TYPE BOND_TYPE
      TYPE(LL_TYPE),INTENT(IN)    :: LL_STRC_
      TYPE(LL_TYPE)               :: LL_STRC
      REAL(8)                     :: UNIT
      REAL(8)                     :: PROTONMASS, ANGSTROM
      INTEGER(4)                  :: NATQ,NATM,NATS
      INTEGER(4)                  :: NBONDM,NBONDS,NCONECT
      INTEGER(4)                  :: NLINK
      INTEGER(4)                  :: NMAP
      TYPE(ATOM_TYPE),ALLOCATABLE :: QATOM(:)
      TYPE(ATOM_TYPE),ALLOCATABLE :: MATOM(:)
      TYPE(ATOM_TYPE),ALLOCATABLE :: SATOM(:)
      TYPE(LINK_TYPE),ALLOCATABLE :: LINK(:)
      TYPE(BOND_TYPE),ALLOCATABLE :: MBOND(:)
      TYPE(BOND_TYPE),ALLOCATABLE :: SBOND(:)
      INTEGER(4)                  :: IATQ,IATM,IATS,IATM1,IATM2,IATS1,IATS2,IAT
      INTEGER(4)                  :: IBONDM,IBONDS,ILINK,I,IRES,J,IVAR,ICONECT
      INTEGER(4)                  :: IZ
      LOGICAL(4)                  :: TCHK,TCHK2
      LOGICAL(4)     ,ALLOCATABLE :: TFREEZE(:)
      INTEGER(4)     ,ALLOCATABLE :: LINKARRAY(:,:)
      INTEGER(4)     ,ALLOCATABLE :: MAPARRAY(:,:)
      CHARACTER(32)               :: CH32SVAR1,DUMMY
      CHARACTER(40)               :: CH40SVAR
      CHARACTER(32)               :: FF
      CHARACTER(255)              :: PARMFILE, TOPFILE, MMFILE
      CHARACTER(4)                :: RESNAME
      REAL(8)                     :: MMS(3)            !ORIGIN OF THE MM SPHERE
      REAL(8)                     :: MSR(1)            !BOUNDARY RADIUS OF MM SPHERE
      REAL(8)                     :: T(9)
      INTEGER(4)                  :: ITVEC(3)
!     ******************************************************************    
                          CALL TRACE$PUSH('STRCIN_SOLVENT')
      LL_STRC=LL_STRC_
      CALL CONSTANTS('U',PROTONMASS)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
!    
!     ==================================================================
!     ==  READ UNIT FROM BLOCK !STRUCTURE!GENERIC                     ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')
      CALL LINKEDLIST$GET(LL_STRC,'LUNIT',1,UNIT)
!
!     ==================================================================
!     ==  READ BLOCK SOLVENT                                          ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'QM-MM',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL TRACE$POP()
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_STRC,'QM-MM')
!     CALL LINKEDLIST$REPORT(LL_STRC,6)
!
!     ==================================================================
!     ==  READ MMSPHERE                           HF   02/2007        ==
!     ==================================================================
      MMS=0.D0
      MSR=1000.D0
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'QM-MM')
      CALL LINKEDLIST$EXISTL(LL_STRC,'MMSPHERE',1,TCHK)
      IF(TCHK) THEN
         CALL LINKEDLIST$SELECT(LL_STRC,'MMSPHERE')
         CALL LINKEDLIST$EXISTD(LL_STRC,'O',1,TCHK)      
         IF(.NOT.TCHK) THEN
            PRINT*,"WARNING: NO MMSPHERE ORIGIN SPECIFIED. SWITCH TO 0,0,0!"
            MMS=0.D0
         ELSE
            CALL LINKEDLIST$GET(LL_STRC,'O',1,MMS)
            CALL CONSTANTS("ANGSTROM",ANGSTROM)
            MMS(:)=MMS(:)*ANGSTROM ! CONVERT TO ATOMIC UNITS
         END IF

         CALL LINKEDLIST$EXISTD(LL_STRC,'RADIUS[A]',1,TCHK)
         IF(.NOT.TCHK) THEN
            PRINT*,"WARNING: NO BOUNDARY RADIUS SPECIFIED. SWITCH TO 1000A0!"
            MSR=1000.D0
         ELSE
            CALL LINKEDLIST$GET(LL_STRC,'RADIUS[A]',1,MSR)
            CALL CONSTANTS("ANGSTROM",ANGSTROM)
            MSR=MSR*ANGSTROM ! CONVERT TO ATOMIC UNITS
            PRINT*,"************************************************"
            PRINT*,"ORIGIN [AU] :",MMS
            PRINT*,"ORIGIN [A] :",MMS/ANGSTROM
            PRINT*,"BOUNDARY RADIUS [AU] :",MSR
            PRINT*,"BOUNDARY RADIUS [A] :",MSR/ANGSTROM
            PRINT*,"************************************************"
         END IF
      END IF
! WARNING: CLASSICAL SYSTEM (QMMM OR SHADOW) MUST BE SELECTED BEFORE WRITING VARIABLES
!!!!      CALL CLASSICAL$SELECT('QMMM')
!!!!            CALL CLASSICAL$SETL4('MMSPHERE',.TRUE.)
!!!!      CALL CLASSICAL$SETR8A('MMS',3,MMS)
!!!!      CALL CLASSICAL$SETR8A('MSR',1,MSR)
!!!!      CALL CLASSICAL$SELECT('SHADOW')
!!!!            CALL CLASSICAL$SETL4('MMSPHERE',.TRUE.)
!!!!      CALL CLASSICAL$SETR8A('MMS',3,MMS)
!!!!      CALL CLASSICAL$SETR8A('MSR',1,MSR)
!STOP 'FORCED STOP'


!     ==================================================================
!     ==  READ FORCEFIELD AND INPUT FILE                              ==
!     ==================================================================
!     --- READ IN FORCEFIELD. IF NO FORCEFIELD IS GIVEN USE UFF ----
      CALL LINKEDLIST$SELECT(LL_STRC,'~')             !MHK
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')     !MHK
      CALL LINKEDLIST$SELECT(LL_STRC,'QM-MM')         !MHK
      CALL LINKEDLIST$EXISTL(LL_STRC,'FORCEFIELD',1,TCHK)
      IF(.NOT.TCHK) THEN
         PRINT*,"WARNING: NO FORCEFIELD SPECIFIED. SWITCH TO UFF."
         FF='UFF'
      ELSE
         CALL LINKEDLIST$SELECT(LL_STRC,'FORCEFIELD')
         CALL LINKEDLIST$EXISTD(LL_STRC,'FF',1,TCHK)
         IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('NO FORCEFIELD SPECIFIED')
            CALL ERROR$MSG('!STRUCTURE!QM-MM!FORCEFIELD:FF IS MANDATORY')
            CALL ERROR$STOP('STRCIN_SOLVENT')
         END IF
         CALL LINKEDLIST$GET(LL_STRC,'FF',1,FF)
!        ---- IS MM-FILE SPECIFIED? ---
         CALL LINKEDLIST$EXISTD(LL_STRC,'MMFILE',1,TCHK)
         IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'MMFILE',1,MMFILE)
            CALL FILEHANDLER$SETFILE('MMSTRC',.FALSE.,MMFILE)
         ELSE
            CALL FILEHANDLER$SETFILE('MMSTRC',.TRUE.,-'.PDB')
         END IF
         CALL FILEHANDLER$SETSPECIFICATION('MMSTRC','STATUS','OLD')
         CALL FILEHANDLER$SETSPECIFICATION('MMSTRC','POSITION','REWIND')
         CALL FILEHANDLER$SETSPECIFICATION('MMSTRC','ACTION','READ')
         CALL FILEHANDLER$SETSPECIFICATION('MMSTRC','FORM','FORMATTED')
!        ---- IS PARMAMETER-FILE SPECIFIED? ---
         CALL LINKEDLIST$EXISTD(LL_STRC,'PARMFILE',1,TCHK)
         IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'PARMFILE',1,PARMFILE)
            CALL FORCEFIELD$SETCH('PARMFILE',PARMFILE)
         END IF
!        ---- IS TOPOLOGY-FILE SPECIFIED? ---
         CALL LINKEDLIST$EXISTD(LL_STRC,'TOPFILE',1,TCHK)
         IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'TOPFILE',1,TOPFILE)
            CALL FORCEFIELD$SETCH('TOPFILE',TOPFILE)
         END IF
         CALL LINKEDLIST$SELECT(LL_STRC,'..')
      END IF
!   --- TRANSFER FORCE FIELD TYPE TO THE MODULES
      CALL FORCEFIELD$SETCH('FORCEFIELD',FF)  !IS THIS NECESSARY? CHECK THIS LATER!
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$SETCH('FF',FF)
      CALL CLASSICAL$SELECT('SHADOW')
      CALL CLASSICAL$SETCH('FF',FF)
!     ==================================================================
!     ==  READ TOP- AND PARM-FILE IF AMBER FORCEFIELD IS CHOSEN       ==
!     ==================================================================
      IF(FF.EQ.'AMBER') THEN
         CALL FORCEFIELD$READ_PARMFILE
         CALL FORCEFIELD$READ_TOPFILE
      END IF
!
!     ==================================================================
!     ==  READ CELL VECTORS                                           ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')             !MHK
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')     !MHK
      CALL LINKEDLIST$SELECT(LL_STRC,'QM-MM')         !MHK
      CALL LINKEDLIST$EXISTL(LL_STRC,'LATTICE',1,TCHK)
      IF(TCHK) THEN
         CALL LINKEDLIST$SELECT(LL_STRC,'LATTICE')
         CALL LINKEDLIST$EXISTD(LL_STRC,'T',1,TCHK)
         CALL LINKEDLIST$EXISTD(LL_STRC,'T[A]',1,TCHK2)
         IF(TCHK.AND.TCHK2) THEN
           CALL ERROR$MSG('SPECIFY EITHER T= OR T[A]=, BUT NOT BOTH')
           CALL ERROR$STOP('STRCIN_SOLVENT')
         END IF
         IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'T',1,T)
         ELSE IF(TCHK2) THEN
            CALL LINKEDLIST$GET(LL_STRC,'T[A]',1,T)
            T(:) = T(:) * ANGSTROM
         ELSE
            T=(/1.D+4,0.D0,0.D0,0.D0,1.D+4,0.D0,0.D0,0.D0,1.D+4/)
         END IF
         CALL CLASSICAL$SELECT('QMMM')
         CALL CLASSICAL$SETR8A('CELL(0)',9,T)
         CALL LINKEDLIST$SELECT(LL_STRC,'..')
      END IF
!
!     ==================================================================
!     ==  READ #(ATOMS) AND  #(BONDS)                                 ==
!     ==================================================================
      CALL ATOMLIST$NATOM(NATQ)
      NATS=NATQ
      IF(FF.EQ.'UFF') THEN
         CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NATM)
      ELSE IF(FF.EQ.'AMBER') THEN
         CALL FORCEFIELD$READ_MMSTRC
         CALL FORCEFIELD$GETI4('NMMATOM',NATM)
      END IF
      CALL LINKEDLIST$NLISTS(LL_STRC,'LINK',NLINK)
      ALLOCATE(MATOM(NATM))
      ALLOCATE(SATOM(NATS))
      ALLOCATE(QATOM(NATQ))
!
!     ==================================================================
!     ==  COLLECT QM ATOMS                                            ==
!     ==================================================================
      DO IATQ=1,NATQ
        CALL ATOMLIST$GETCH('NAME',IATQ,QATOM(IATQ)%NAME)
        CALL ATOMLIST$GETR8A('R(0)',IATQ,3,QATOM(IATQ)%R)
        CALL ATOMLIST$GETR8('MASS',IATQ,QATOM(IATQ)%M)
        QATOM(IATQ)%ELEMENT=QATOM(IATQ)%NAME(1:2)
        IF(QATOM(IATQ)%ELEMENT(2:2).EQ.'_') QATOM(IATQ)%ELEMENT(2:2)=' '
        QATOM(IATQ)%FFTYPE=' '
        QATOM(IATQ)%Q=0.D0
        QATOM(IATQ)%QMSATOM=0
      ENDDO
!
!     ==================================================================
!     ==================================================================
!     ==  READ ATOMS                                                  ==
!     ==================================================================
!     ==================================================================
      IF (FF.EQ.'UFF') THEN
        IATS=0
        ALLOCATE(TFREEZE(NATM))
        TFREEZE(:)=.FALSE.       !ALL ATOMS ARE FREE BY DEFAULT
        DO IATM=1,NATM
          CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IATM)
!    
!         ==  ATOM NAME  =================================================
          CALL LINKEDLIST$EXISTD(LL_STRC,'NAME',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('VARIABLE !QM-MM!ATOM:NAME NOT FOUND')
            CALL ERROR$STOP('STRCIN_SOLVENT') 
          END IF
          CALL LINKEDLIST$GET(LL_STRC,'NAME',1,MATOM(IATM)%NAME)
          MATOM(IATM)%ELEMENT=MATOM(IATM)%NAME(1:2)                          !*
          IF(MATOM(IATM)%ELEMENT(2:2).EQ.'_') MATOM(IATM)%ELEMENT(2:2)=' '   !*
!    
!         ==  QMATOM  =====================================================
          CALL LINKEDLIST$EXISTD(LL_STRC,'QMATOM',1,TCHK)
          IATQ=0
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'QMATOM',1,CH32SVAR1)
            DO I=1,NATQ
              IF(CH32SVAR1.EQ.QATOM(I)%NAME) THEN
                IATQ=I
                EXIT
              END IF
            ENDDO
            IF(IATQ.EQ.0) THEN
              CALL ERROR$MSG('QMATOM NOT FOUND')
              CALL ERROR$CHVAL('QMATOM',CH32SVAR1)
              CALL ERROR$CHVAL('MMATOM',MATOM(IATM)%NAME)
              CALL ERROR$STOP('STRCIN_SOLVENT')
            END IF
            IATS=IATS+1     ! INCREASE COUNTER FOR SHADOW ATOMS
            IF(IATS.GT.NATS) THEN
               CALL ERROR$MSG('#(SHADOW ATOMS) LARGER THAN EXPECTED')
               CALL ERROR$STOP('STRCIN_SOLVENT')
            END IF
            SATOM(IATS)%NAME=QATOM(IATQ)%NAME
!==================================
            SATOM(IATS)%ELEMENT= QATOM(IATQ)%ELEMENT
PRINT*,"SELEMENT= ",SATOM(IATS)%ELEMENT, IATS
!==================================
            QATOM(IATQ)%QMSATOM=IATM
            MATOM(IATM)%QMSATOM=IATS
            SATOM(IATS)%QMSATOM=IATQ
          ELSE
            MATOM(IATM)%QMSATOM=0
          END IF
!    
!         ==  FORCE FIELD ATOM TYPE========================================
          CALL LINKEDLIST$EXISTD(LL_STRC,'FFTYPE',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('VARIABLE !QM-MM!ATOM:FFTYPE NOT FOUND')
            CALL ERROR$STOP('STRCIN_SOLVENT') 
          END IF
          CALL LINKEDLIST$GET(LL_STRC,'FFTYPE',1,MATOM(IATM)%FFTYPE)
          IF(IATQ.NE.0) THEN
            SATOM(IATS)%FFTYPE=MATOM(IATM)%FFTYPE
          END IF
!    
!         ==  POSITION====================================================
          IF(IATQ.EQ.0) THEN
            CALL LINKEDLIST$EXISTD(LL_STRC,'R',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('KEYWORD QM-MM!ATOM:R NOT FOUND')
              CALL ERROR$STOP('STRCIN_SOLVENT') 
           END IF
           CALL LINKEDLIST$GET(LL_STRC,'R',1,MATOM(IATM)%R)
           MATOM(IATM)%R(:)=MATOM(IATM)%R(:)*UNIT
         ELSE
            CALL LINKEDLIST$EXISTD(LL_STRC,'R',1,TCHK)
            IF(TCHK) THEN
              CALL ERROR$MSG('R MUST NOT BE SPECIFIED')
              CALL ERROR$MSG('ATOMIC POSITION IS TAKEN FROM QM ATOM')
              CALL ERROR$STOP('STRCIN_SOLVENT')
            END IF
           MATOM(IATM)%R(:)=QATOM(IATQ)%R
           SATOM(IATS)%R(:)=QATOM(IATQ)%R
         END IF
!    
!        ==  MASS ========================================================
         IF(IATQ.EQ.0) THEN 
            CALL LINKEDLIST$EXISTD(LL_STRC,'M',1,TCHK)
            IF(TCHK) THEN
              CALL LINKEDLIST$GET(LL_STRC,'M',1,MATOM(IATM)%M)
              MATOM(IATM)%M=MATOM(IATM)%M*PROTONMASS
            ELSE
              CH32SVAR1=MATOM(IATM)%FFTYPE(1:2)
              IF(CH32SVAR1(2:2).EQ.'_') CH32SVAR1(2:2)=' '
              CALL PERIODICTABLE$GET(CH32SVAR1,'Z',IZ)
              CALL PERIODICTABLE$GET(IZ,'MASS',MATOM(IATM)%M)
              MATOM(IATM)%M=MATOM(IATM)%M
            END IF
          ELSE
            CALL LINKEDLIST$EXISTD(LL_STRC,'M',1,TCHK)
            IF(TCHK) THEN
              CALL ERROR$MSG('M MUST NOT BE SPECIFIED')
              CALL ERROR$MSG('MASS IS TAKEN FROM QM ATOM')
              CALL ERROR$STOP('STRCIN_SOLVENT')
            END IF
            MATOM(IATM)%M=QATOM(IATQ)%M
            SATOM(IATS)%M=QATOM(IATQ)%M
          END IF
!    
!         ==  CHARGE ======================================================
          CALL LINKEDLIST$EXISTD(LL_STRC,'Q',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'Q',1,MATOM(IATM)%Q)
            MATOM(IATM)%Q=-MATOM(IATM)%Q
          ELSE
            MATOM(IATM)%Q=0.D0
          END IF
          CALL LINKEDLIST$SELECT(LL_STRC,'..')
        ENDDO
!     ==================================================================
!     == READ AMBER ATOMS ==============================================
!     ==================================================================
      ELSE IF(FF.EQ.'AMBER') THEN
        IATS=0
        ALLOCATE(TFREEZE(NATM))
        TFREEZE(:)=.FALSE.       !ALL ATOMS ARE FREE BY DEFAULT
        DO IATM=1, NATM
!         ----- READ NAME OF MM-ATOM NAME= <ATOMNAME>_<ID>    
          MATOM(IATM)%NAME = TRIM(ADJUSTL(MMATOM(IATM)%NAME(1:LEN_TRIM(MMATOM(IATM)%NAME))))&
     &                     //'_'//TRIM(ADJUSTL(.ITOS.MMATOM(IATM)%ID))
          MATOM(IATM)%ELEMENT=TRIM(ADJUSTL(MMATOM(IATM)%ELEMENT))
!
!         ----- CHECK IF MM ATOM IS QM-ATOM
          IATQ=0
          IF(MMATOM(IATM)%FLAG.EQ.'Q') THEN
            DO I=1,NATQ
              IF(TRIM(ADJUSTL(MMATOM(IATM)%QMNAME)).EQ.TRIM(ADJUSTL(QATOM(I)%NAME))) THEN
                IATQ=I
                EXIT
              END IF
            END DO
            IF(IATQ.EQ.0) THEN
              CALL ERROR$MSG('QMATOM NOT FOUND')
              CALL ERROR$CHVAL('QMATOM',MMATOM(IATM)%QMNAME)
              CALL ERROR$CHVAL('MMATOM',MATOM(IATM)%NAME)
              CALL ERROR$STOP('STRCIN_SOLVENT')
            END IF
            IF(MMATOM(IATM)%FLAG.EQ.'Q') IATS=IATS+1     ! INCREASE COUNTER FOR SHADOW ATOMS
            IF(IATS.GT.NATS) THEN
               CALL ERROR$MSG('#(SHADOW ATOMS) LARGER THAN EXPECTED')
               CALL ERROR$STOP('STRCIN_SOLVENT')
            END IF
            SATOM(IATS)%NAME=MMATOM(IATM)%NAME !OR QMNAME ? 
            SATOM(IATS)%ELEMENT=TRIM(ADJUSTL(MMATOM(IATM)%ELEMENT))
            QATOM(IATQ)%QMSATOM=IATM
            MATOM(IATM)%QMSATOM=IATS
            SATOM(IATS)%QMSATOM=IATQ
          ELSE
            MATOM(IATM)%QMSATOM=0
          END IF
!         ----- READIN FORCEFIELD ATOM TYPE AND ATOM CHARGE
          CALL FORCEFIELD$GETRESNUMBER(MMATOM(IATM)%RESNAME,IRES)
          DO I=1,SIZE(TOP_RES(IRES)%ATOM)
            IF(TRIM(ADJUSTL(TOP_RES(IRES)%ATOM(I)%NAME)).EQ.TRIM(ADJUSTL(MMATOM(IATM)%NAME))) THEN
              MATOM(IATM)%FFTYPE =  TRIM(ADJUSTL(TOP_RES(IRES)%ATOM(I)%ATOM))
!             ------  HERE YOU CAN CHOOSE IF THE MM ATOMS GET ZERO CHARGE OR NOT ------------
              MATOM(IATM)%Q      = -TOP_RES(IRES)%ATOM(I)%CHARGE  
!             -------------------------------------------------------------------------------
              EXIT
            END IF
          END DO
          IF(IATQ.NE.0) THEN
            SATOM(IATS)%FFTYPE=MATOM(IATM)%FFTYPE
          END IF
!         ----- POSITIONS -----------------------------
          IF(IATQ.EQ.0) THEN
            MATOM(IATM)%R(:)= MMATOM(IATM)%R(:) * ANGSTROM
          ELSE
            MATOM(IATM)%R(:)=QATOM(IATQ)%R
            SATOM(IATS)%R(:)=QATOM(IATQ)%R
          END IF
!         ----- MASS ----------------------------------
          IF(IATQ.EQ.0) THEN
            CH32SVAR1 = ADJUSTL(MMATOM(IATM)%ELEMENT)
            CALL PERIODICTABLE$GET(CH32SVAR1,'Z',IZ)
            CALL PERIODICTABLE$GET(IZ,'MASS',MATOM(IATM)%M)
          ELSE
            MATOM(IATM)%M=QATOM(IATQ)%M
            SATOM(IATS)%M=QATOM(IATQ)%M
            SATOM(IATS)%Q=MATOM(IATM)%Q 
          END IF
!         ----- IS ATOM FIXED ? -----------------------
          IF(+MMATOM(IATM)%FLAG.EQ.'F') TFREEZE(IATM)=.TRUE. 
        END DO
      END IF

!
!     ==================================================================
!     ==================================================================
!     ==  READ LINKS                                                  ==
!     ==================================================================
!     ==================================================================
      IF(FF.EQ.'UFF') THEN
        CALL LINKEDLIST$NLISTS(LL_STRC,'LINK',NLINK)
        ALLOCATE(LINK(NLINK))
        DO ILINK=1,NLINK
          CALL LINKEDLIST$SELECT(LL_STRC,'LINK',ILINK)
!  
!         == JOINT ATOM ==================================================
!         == THE JOINT ATOM IS PRESENT IN THE QM AND THE MM SYSTEM =======
          CALL LINKEDLIST$EXISTD(LL_STRC,'MMJOINT',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('VARIABLE !QM-MM!LINK:MMJOINT NOT FOUND')
            CALL ERROR$STOP('STRCIN_SOLVENT') 
          END IF
          CALL LINKEDLIST$GET(LL_STRC,'MMJOINT',1,CH32SVAR1)
          LINK(ILINK)%MJOINT=0
          DO IATM=1,NATM
            IF(CH32SVAR1.EQ.MATOM(IATM)%NAME) THEN
              LINK(ILINK)%MJOINT=IATM
              EXIT
            END IF
          ENDDO
          IF(LINK(ILINK)%MJOINT.EQ.0) THEN
            CALL ERROR$MSG('!QM-MM!LINK:MMJOINT IS NOT A MM ATOM')
            CALL ERROR$CHVAL('MMJOINT',CH32SVAR1)
            CALL ERROR$STOP('STRCIN_SOLVENT') 
          END IF
          LINK(ILINK)%SJOINT=MATOM(LINK(ILINK)%MJOINT)%QMSATOM
          LINK(ILINK)%QJOINT=SATOM(LINK(ILINK)%SJOINT)%QMSATOM
!  
!         == MM ATOM =====================================================
          CALL LINKEDLIST$EXISTD(LL_STRC,'MMATOM',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('VARIABLE !QM-MM!LINK:MMATOM NOT FOUND')
            CALL ERROR$STOP('STRCIN_SOLVENT') 
          END IF
          CALL LINKEDLIST$GET(LL_STRC,'MMATOM',1,CH32SVAR1)
          LINK(ILINK)%MATOM=0
          DO IATM=1,NATM
            IF(CH32SVAR1.EQ.MATOM(IATM)%NAME) THEN
              LINK(ILINK)%MATOM=IATM
              EXIT
            END IF
          ENDDO
          IF(LINK(ILINK)%MATOM.EQ.0) THEN
            CALL ERROR$MSG('!QM-MM!LINK:MMATOM IS NOT A MM ATOM')
            CALL ERROR$CHVAL('MMATOM',CH32SVAR1)
            CALL ERROR$STOP('STRCIN_SOLVENT') 
          END IF
!            
!         == QM ATOM =====================================================
          CALL LINKEDLIST$EXISTD(LL_STRC,'QMATOM',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('VARIABLE !QM-MM!LINK:QMATOM NOT FOUND')
            CALL ERROR$STOP('STRCIN_SOLVENT') 
          END IF
          CALL LINKEDLIST$GET(LL_STRC,'QMATOM',1,CH32SVAR1)
          DO IATQ=1,NATQ
            IF(CH32SVAR1.EQ.QATOM(IATQ)%NAME) THEN
              LINK(ILINK)%QATOM=IATQ
              EXIT
            END IF
          ENDDO
          IF(LINK(ILINK)%QATOM.EQ.0) THEN
            CALL ERROR$MSG('!QM-MM!LINK:QMATOM IS NOT A QM ATOM')
            CALL ERROR$CHVAL('QMATOM',CH32SVAR1)
            CALL ERROR$STOP('STRCIN_SOLVENT') 
          END IF
!  
!         == SHADOW ATOM =================================================
          IATS=IATS+1
          IF(IATS.GT.NATS) THEN
             CALL ERROR$MSG('#(SHADOW ATOMS) LARGER THAN EXPECTED')
             CALL ERROR$STOP('STRCIN_SOLVENT')
          END IF
          LINK(ILINK)%SATOM=IATS
!  
!         == SHADOW ATOM FFTYPE ==========================================
          CALL LINKEDLIST$EXISTD(LL_STRC,'SHFFTYPE',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('VARIABLE !QM-MM!LINK!SHFFTYPE NOT FOUND')
            CALL ERROR$STOP('STRCIN_SOLVENT') 
          END IF
          CALL LINKEDLIST$GET(LL_STRC,'SHFFTYPE',1,SATOM(IATS)%FFTYPE)
!         
!         == SHADOW ATOM NAME,POSITION,MASS,CHARGE,QMSATOM================
          IATQ=LINK(ILINK)%QATOM
          SATOM(IATS)%NAME=QATOM(IATQ)%NAME
!=============================================0
          SATOM(IATS)%ELEMENT='H '
!=============================================0
          SATOM(IATS)%R=QATOM(IATQ)%R
          SATOM(IATS)%M=QATOM(IATQ)%M
          SATOM(IATS)%Q=QATOM(IATQ)%Q
!  
          IATS=LINK(ILINK)%SATOM
          IATQ=LINK(ILINK)%QATOM
          IATM=LINK(ILINK)%MATOM
          SATOM(IATS)%QMSATOM=IATQ
          QATOM(IATQ)%QMSATOM=IATM
          MATOM(IATM)%QMSATOM=IATS
!         
          CALL LINKEDLIST$SELECT(LL_STRC,'..')
        ENDDO
!     ===============================================================
!     =========== AMBER ==== READ LINKS =============================
!     ===============================================================
      ELSE IF(FF.EQ.'AMBER') THEN
                        CALL TRACE$PASS('READING LINK ATOMS IN AMBER')
!       == GET NUMBER OF LINK ATOMS ====================================
        NLINK=0
        DO I=1,SIZE(MMATOM)
          IF(MMATOM(I)%FLAG.EQ.'L') NLINK = NLINK + 1
        END DO
        ALLOCATE(LINK(NLINK))
!       == FILL LINK VARIABLE ==========================================
        ILINK=0
        DO I=1,SIZE(MMATOM)
          IF(MMATOM(I)%FLAG.EQ.'L') THEN
            ILINK= ILINK + 1
!           ---- JOINT ATOMS ------------
            DUMMY=TRIM(ADJUSTL(MMATOM(I)%QMNAME))
            DO J=1,NATQ
               IF(TRIM(ADJUSTL(QATOM(J)%NAME)).EQ.TRIM(ADJUSTL(DUMMY))) THEN
                  LINK(ILINK)%MJOINT=QATOM(J)%QMSATOM
                  EXIT
               END IF
            END DO
            IF(LINK(ILINK)%MJOINT.EQ.0) THEN
               CALL ERROR$MSG('!QM-MM!LINK:MMJOINT IS NOT A MM ATOM')
               CALL ERROR$CHVAL('MMJOINT',CH32SVAR1)
               CALL ERROR$STOP('STRCIN_SOLVENT') 
            END IF
            LINK(ILINK)%SJOINT=MATOM(LINK(ILINK)%MJOINT)%QMSATOM
            LINK(ILINK)%QJOINT=SATOM(LINK(ILINK)%SJOINT)%QMSATOM
!           ---- MM ATOM -----------   
            LINK(ILINK)%MATOM=I
!           ---- QM ATOM -----------
            DO IATQ=1,NATQ
               IF(TRIM(ADJUSTL(MMATOM(I)%LNAME)).EQ.TRIM(ADJUSTL(QATOM(IATQ)%NAME))) THEN
                  LINK(ILINK)%QATOM=IATQ
                  EXIT
               END IF
            ENDDO
            IF(LINK(ILINK)%QATOM.EQ.0) THEN
               CALL ERROR$MSG('!QM-MM!LINK:QMATOM IS NOT A QM ATOM')
               CALL ERROR$CHVAL('QMATOM',CH32SVAR1)
               CALL ERROR$STOP('STRCIN_SOLVENT') 
            END IF
            IATS=IATS+1
            IF(IATS.GT.NATS) THEN
               CALL ERROR$MSG('#(SHADOW ATOMS) LARGER THAN EXPECTED')
               CALL ERROR$I4VAL('IATS',IATS)
               CALL ERROR$I4VAL('NATS',NATS)
               CALL ERROR$STOP('STRCIN_SOLVENT')
            END IF
!           ---- SHADOW ATOM -------
            LINK(ILINK)%SATOM=IATS
!           ---- SHADOW ATOM FFTYPE-
            SATOM(IATS)%FFTYPE= 'L1'
!           ---- SHADOW ATOM NAME,POSITION,MASS,CHARGE,QMSATOM----------------
            IATQ=LINK(ILINK)%QATOM
            SATOM(IATS)%NAME='H'!MMATOM(I)%NAME
            SATOM(IATS)%ELEMENT='H '
            SATOM(IATS)%R=QATOM(IATQ)%R
            SATOM(IATS)%M=QATOM(IATQ)%M
            SATOM(IATS)%Q=QATOM(IATQ)%Q
!
            IATS=LINK(ILINK)%SATOM
            IATQ=LINK(ILINK)%QATOM
            IATM=LINK(ILINK)%MATOM
!
            SATOM(IATS)%QMSATOM=IATQ
            QATOM(IATQ)%QMSATOM=IATM
            MATOM(IATM)%QMSATOM=IATS
          END IF
        END DO
      END IF
      IF(IATS.NE.NATS) THEN
        CALL ERROR$MSG('CONSISTENCY CHECK IATS FAILED')
        CALL ERROR$I4VAL('IATS',IATS)
        CALL ERROR$I4VAL('NATS',NATS)
        CALL ERROR$STOP('STRCIN_SOLVENT')
      END IF
!
!     ==================================================================
!     ==================================================================
!     ==  READ BONDS                                                  ==
!     ==================================================================
!     ==================================================================
      NCONECT=0
      IF(FF.EQ.'UFF') THEN
        CALL LINKEDLIST$NLISTS(LL_STRC,'BOND',NBONDM)
        ALLOCATE(MBOND(NBONDM))
        DO IBONDM=1,NBONDM
          CALL LINKEDLIST$SELECT(LL_STRC,'BOND',IBONDM)
!         == FIRST ATOM ==================================================
          MBOND(IBONDM)%ATOM1=0
          CALL LINKEDLIST$GET(LL_STRC,'ATOM1',1,CH40SVAR)
          CALL STRCIN_RESOLVEEXTENDEDNAME(CH40SVAR,CH32SVAR1,ITVEC)
          MBOND(IBONDM)%IT(:)=-ITVEC(:)     
          DO IATM=1,NATM
            IF(CH32SVAR1.EQ.MATOM(IATM)%NAME) THEN
              MBOND(IBONDM)%ATOM1=IATM
              EXIT
            END IF
          ENDDO
          IF(MBOND(IBONDM)%ATOM1.EQ.0) THEN
            CALL ERROR$MSG('FIRST ATOM IN BOND NOT IN THE M-ATOM LIST')
            CALL ERROR$CHVAL('ATOM1',CH32SVAR1)
            CALL ERROR$STOP('STRCIN_SOLVENT; !STRUCTURE!QM-MM!BOND')
          END IF         
!   
!         == SECOND ATOM =================================================
          MBOND(IBONDM)%ATOM2=0
          CALL LINKEDLIST$GET(LL_STRC,'ATOM2',1,CH40SVAR)
          CALL STRCIN_RESOLVEEXTENDEDNAME(CH40SVAR,CH32SVAR1,ITVEC)
          MBOND(IBONDM)%IT(:)=MBOND(IBONDM)%IT(:)+ITVEC
          DO IATM=1,NATM
            IF(CH32SVAR1.EQ.MATOM(IATM)%NAME) THEN
              MBOND(IBONDM)%ATOM2=IATM
              EXIT
            END IF
          ENDDO
          IF(MBOND(IBONDM)%ATOM2.EQ.0) THEN
            CALL ERROR$MSG('SECOND ATOM IN BOND NOT IN THE M-ATOM LIST')
            CALL ERROR$CHVAL('ATOM2',CH32SVAR1)
            CALL ERROR$STOP('STRCIN_SOLVENT; !STRUCTURE!QM-MM!BOND')
          END IF         
!   
!         == BOND ORDER ==================================================
          CALL LINKEDLIST$EXISTD(LL_STRC,'BO',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_STRC,'BO',1,MBOND(IBONDM)%BO)
          ELSE
            MBOND(IBONDM)%BO=1
          END IF
          CALL LINKEDLIST$SELECT(LL_STRC,'..')
        ENDDO
!   == AMBER FORCEFIELD : BONDS ====================================================
      ELSE IF(FF.EQ.'AMBER') THEN
!       ----- GET NUMBER OF BONDS
        IAT=1
        NBONDM=0
!       ----- READ #BONDS FROM THE TOPOLOGY INFORMATIONS
        DO WHILE(IAT.LT.SIZE(MMATOM))
          RESNAME= MMATOM(IAT)%RESNAME
          CALL FORCEFIELD$GETRESNUMBER(RESNAME,IVAR)
          IBONDM= SIZE(TOP_RES(IVAR)%BOND)
          DO I=1,SIZE(TOP_RES(IVAR)%BOND)  !EXCLUDE PEPTIDE BONDS
            IF(TRIM(TOP_RES(IVAR)%BOND(I)%ATOM1).EQ.'+N  '.OR.TRIM(TOP_RES(IVAR)%BOND(I)%ATOM2)&
                 &.EQ.'+N  ') IBONDM = IBONDM -1
            IF(TRIM(TOP_RES(IVAR)%BOND(I)%ATOM1).EQ.'-C  '.OR.TRIM(TOP_RES(IVAR)%BOND(I)%ATOM2)&
                 &.EQ.'-C  ') IBONDM = IBONDM -1
          END DO
!         ----- ONLY IF THE ENTRY BELONGS TO AN AMINO ACID ADD THE PEPTIDE BOND
          IF(MMATOM(IAT)%KEYWORD.EQ.'ATOM  ') THEN
             NBONDM= NBONDM + IBONDM +1
          ELSE IF(MMATOM(IAT)%KEYWORD.EQ.'HETATM') THEN
             NBONDM = NBONDM + IBONDM
          END IF
!         -------------
          IVAR = SIZE(TOP_RES(IVAR)%ATOM)
          IAT  = IAT + IVAR
!         ----- CHECKS IF THE PEPTIDE CHAIN ID HAS CHANGED. ONLY FOR AMINO ACIDS
          IF(MMATOM(IAT-1)%KEYWORD.EQ.'ATOM  ') THEN
             IF(IAT.GE.SIZE(MMATOM)) THEN
               NBONDM = NBONDM - 1
             ELSE
               IF(MMATOM(IAT-1)%CHAINID.NE.MMATOM(IAT)%CHAINID) NBONDM = NBONDM - 1
             END IF
          END IF
        ENDDO
!             ----- READ #BONDS FROM THE CONECT INFORMATION
        DO I=1,SIZE(MMCONECT)
           IVAR=0
           IF(MMCONECT(I)%CATOM.EQ.0) CYCLE
           IF(MMCONECT(I)%ATOM4.NE.0) IVAR=IVAR+1
           IF(MMCONECT(I)%ATOM3.NE.0) IVAR=IVAR+1
           IF(MMCONECT(I)%ATOM2.NE.0) IVAR=IVAR+1
           IF(MMCONECT(I)%ATOM1.NE.0) IVAR=IVAR+1
           NCONECT=NCONECT+IVAR
        END DO
PRINT*,"FLAG: ",NBONDM,NCONECT
        ALLOCATE(MBOND(NBONDM+NCONECT))
        MBOND(:)%ATOM1=0
        MBOND(:)%ATOM2=0
        DO I=1,SIZE(MBOND)
           MBOND(I)%IT(:)=0
        END DO
        MBOND(:)%BO=0.D0!
!       ----- READ BONDS FROM MMATOM-ARRAY
        IAT=1
        IBONDM=1
        DO WHILE(IAT.LT.SIZE(MMATOM))
          RESNAME= MMATOM(IAT)%RESNAME  !WHICH RESIDUE
          CALL FORCEFIELD$GETRESNUMBER(RESNAME,IVAR)
          DO I=1,SIZE(TOP_RES(IVAR)%BOND)  !GO OVER ALL BONDS OF THE RESIDUE
             IF(TRIM(TOP_RES(IVAR)%BOND(I)%ATOM1).EQ.'+N  '.OR.TRIM(TOP_RES(IVAR)%BOND(I)%ATOM2).EQ.'+N  ') CYCLE
             IF(TRIM(TOP_RES(IVAR)%BOND(I)%ATOM1).EQ.'-C  '.OR.TRIM(TOP_RES(IVAR)%BOND(I)%ATOM2).EQ.'-C  ') CYCLE
             DO J=IAT, IAT+SIZE(TOP_RES(IVAR)%ATOM)-1 !CHECK ALL ATOMS OF THIS RESIDUE IN MMATOM
                IF(ADJUSTL(MMATOM(J)%NAME).EQ.TOP_RES(IVAR)%BOND(I)%ATOM1) THEN
                   MBOND(IBONDM)%ATOM1 = J
                END IF
                IF(ADJUSTL(MMATOM(J)%NAME).EQ.TOP_RES(IVAR)%BOND(I)%ATOM2) THEN
                   MBOND(IBONDM)%ATOM2 = J
                END IF
             END DO
             MBOND(IBONDM)%BO = 1.D0
             IBONDM = IBONDM + 1
          END DO
          IAT= IAT+SIZE(TOP_RES(IVAR)%ATOM)
          IF(IAT.GT.SIZE(MMATOM)) EXIT
!         ---- CONNECT DIFFERENT AMINOACIDS IF THEY ARE IN THE SAME CHAIN
          IF(MMATOM(IAT)%KEYWORD.EQ.'ATOM  ') THEN
            IF((MMATOM(IAT-1)%CHAINID.EQ.MMATOM(IAT)%CHAINID)) THEN

!             ---- GET LAST C -ATOM
              CALL FORCEFIELD$GETRESNUMBER(MMATOM(IAT-1)%RESNAME,IVAR)
              J=SIZE(TOP_RES(IVAR)%ATOM)
              DO I= IAT-J, IAT-1
                IF(TRIM(ADJUSTL(MMATOM(I)%NAME)).EQ.'C') THEN
                   MBOND(IBONDM)%ATOM1= I
                   PRINT*,"FLAG: ATOM FOUND! ",MMATOM(I)%NAME, I
                   EXIT
                END IF
              ENDDO
!             ---- GET NEXT (ACTUAL) N - ATOM  
              CALL FORCEFIELD$GETRESNUMBER(MMATOM(IAT)%RESNAME,IVAR)
              J=SIZE(TOP_RES(IVAR)%ATOM)
              DO I= IAT, IAT+J-1
                IF(TRIM(ADJUSTL(MMATOM(I)%NAME)).EQ.'N') THEN
                   MBOND(IBONDM)%ATOM2= I
                   PRINT*,"FLAG: ATOM FOUND! ",MMATOM(I)%NAME, I
                   EXIT
                END IF
              END DO
              MBOND(IBONDM)%BO = 1
              IBONDM = IBONDM + 1
            END IF
         END IF
        END DO

        ICONECT=NBONDM
        DO I=1,SIZE(MMCONECT)
           IF(MMCONECT(I)%CATOM.NE.0) THEN
              IF(MMCONECT(I)%ATOM4.NE.0) THEN
                 ICONECT=ICONECT+1
                 MBOND(ICONECT)%ATOM1=MMCONECT(I)%CATOM
                 MBOND(ICONECT)%ATOM2=MMCONECT(I)%ATOM4
                 MBOND(ICONECT)%BO=1
                 CYCLE
              END IF
              IF(MMCONECT(I)%ATOM3.NE.0) THEN
                 ICONECT=ICONECT+1
                 MBOND(ICONECT)%ATOM1=MMCONECT(I)%CATOM
                 MBOND(ICONECT)%ATOM2=MMCONECT(I)%ATOM3
                 MBOND(ICONECT)%BO=1
                 CYCLE
              END IF
              IF(MMCONECT(I)%ATOM2.NE.0) THEN
                 ICONECT=ICONECT+1
                 MBOND(ICONECT)%ATOM1=MMCONECT(I)%CATOM
                 MBOND(ICONECT)%ATOM2=MMCONECT(I)%ATOM2
                 MBOND(ICONECT)%BO=1
                 CYCLE
              END IF
              IF(MMCONECT(I)%ATOM1.NE.0) THEN
                 ICONECT=ICONECT+1
                 MBOND(ICONECT)%ATOM1=MMCONECT(I)%CATOM
                 MBOND(ICONECT)%ATOM2=MMCONECT(I)%ATOM1
                 MBOND(ICONECT)%BO=1
                 CYCLE
              END IF
           END IF
        END DO
      END IF
!
!     ==================================================================
!     ==  RESOLVE SHADOW BONDS                                        ==
!     ==================================================================
!     == COUNT #(SHADOW BONDS) =========================================
      NBONDS=0
      DO IBONDM=1,NBONDM+NCONECT
        IATM1=MBOND(IBONDM)%ATOM1
        IATM2=MBOND(IBONDM)%ATOM2
        IATS1=MATOM(IATM1)%QMSATOM
        IATS2=MATOM(IATM2)%QMSATOM
        IF(IATS1.NE.0.AND.IATS2.NE.0) NBONDS=NBONDS+1
      ENDDO
      ALLOCATE(SBOND(NBONDS))
!
!     == DETERMINE SHADOW BONDS ========================================
      IBONDS=0
      DO IBONDM=1,NBONDM+NCONECT
        IATM1=MBOND(IBONDM)%ATOM1
        IATM2=MBOND(IBONDM)%ATOM2
        IATS1=MATOM(IATM1)%QMSATOM
        IATS2=MATOM(IATM2)%QMSATOM
        IF(IATS1.NE.0.AND.IATS2.NE.0) THEN
          IBONDS=IBONDS+1   ! BOND LIES ENTIRELY IN THE REACTION CENTER
                            ! OR IS A LINK BOND
          SBOND(IBONDS)%ATOM1=IATS1
          SBOND(IBONDS)%ATOM2=IATS2
          SBOND(IBONDS)%BO=MBOND(IBONDM)%BO
          SBOND(IBONDS)%IT=MBOND(IBONDM)%IT
        END IF
      ENDDO
!
!     ==================================================================
!     ==================================================================
!     ==  TRANSFER DATA TO MODULES                                    ==
!     ==================================================================
!     ==================================================================
!
!     ==================================================================
!     ==  SET LINKS                                                   ==
!     ==================================================================
      ALLOCATE(LINKARRAY(6,NLINK))
      ALLOCATE(MAPARRAY(3,NATS))
      DO ILINK=1,NLINK
        LINKARRAY(1,ILINK)=LINK(ILINK)%QJOINT
        LINKARRAY(2,ILINK)=LINK(ILINK)%MJOINT
        LINKARRAY(3,ILINK)=LINK(ILINK)%SJOINT
        LINKARRAY(4,ILINK)=LINK(ILINK)%QATOM
        LINKARRAY(5,ILINK)=LINK(ILINK)%MATOM
        LINKARRAY(6,ILINK)=LINK(ILINK)%SATOM
      ENDDO
!     == MAPARRAY ONLY INCLUDES ATOMS NOT PARTICIPATING IN LINK BONDS
      NMAP=0
      DO IATS=1,NATS
        TCHK=.FALSE.
        DO ILINK=1,NLINK
          TCHK=TCHK.OR.(IATS.EQ.LINKARRAY(6,ILINK))
          TCHK=TCHK.OR.(IATS.EQ.LINKARRAY(3,ILINK))
          IF(TCHK) EXIT
        ENDDO
        IF(TCHK) CYCLE
        NMAP=NMAP+1
        MAPARRAY(1,NMAP)=SATOM(IATS)%QMSATOM
        MAPARRAY(2,NMAP)=QATOM(MAPARRAY(1,NMAP))%QMSATOM
        MAPARRAY(3,NMAP)=IATS
      ENDDO
      CALL QMMM$SETI4A('MAP',3*NMAP,MAPARRAY(:,1:NMAP))
      CALL QMMM$SETI4A('LINK',6*NLINK,LINKARRAY)
      DEALLOCATE(LINKARRAY)
      DEALLOCATE(MAPARRAY)
PRINT*,"FLAG: AFTER SETTING LINKS IN QMMM"
NBONDM=NBONDM+NCONECT
!
!     ==================================================================
!     ==  DEFINE CLASSICAL OBJECT                                     ==
!     ==================================================================
      CALL CLASSICAL$SELECT('QMMM')
      CALL CLASSICAL$SETL4A('TFREEZE',NATM,TFREEZE)
      CALL STRCIN_SOLVENT_SETM(NATM,MATOM,NBONDM,MBOND)
! PRINT*,"---------------------------"
! PRINT*,"   MATOM:"
! DO I=1,SIZE(MATOM)
!    WRITE(*,FMT='(A32, 3F10.4, 2F15.5, A10, A2, I6)')  MATOM(I)
! END DO
! PRINT*,"---------------------------"
! PRINT*,"   SATOM:",ALLOCATED(SATOM)
! DO I=1,SIZE(SATOM)
!    WRITE(*,FMT='(I6,A32, 3F10.4, 2F15.5, A10, A2, I6)') I,  SATOM(I)
! END DO
! PRINT*,"---------------------------"
! PRINT*,"   QATOM:"
! DO I=1,SIZE(QATOM)
!    WRITE(*,FMT='(A32, 3F10.4, 2F15.5, A10, A2, I6)')  QATOM(I)
! END DO
! PRINT*,"---------------------------"
! PRINT*,"   MBOND:",SIZE(MBOND)
! DO I=1,SIZE(MBOND)
!    WRITE(*,FMT='(2I8,F5.1)') MBOND(I)
! END DO
! PRINT*,"---------------------------"
! PRINT*,"   SBOND:"
! DO I=1,SIZE(SBOND)
!    WRITE(*,FMT='(2I8,F5.1)') SBOND(I)
! END DO
! PRINT*,"---------------------------"
! DO I=1,NLINK
!    WRITE(*,FMT='(6I8)') LINKARRAY(:,I)
! END DO
! !IF(ALLOCATED(LINKARRAY)) DEALLOCATE(LINKARRAY)
! STOP 'FORCED STOP IN READIN_SOLVENT'

      CALL CLASSICAL$SELECT('SHADOW')
      CALL STRCIN_SOLVENT_SETM(NATS,SATOM,NBONDS,SBOND)
!
!     ==================================================================
!     ==  DEALLOCATE AND RETURN                                       ==
!     ==================================================================
      DEALLOCATE(MATOM)
      DEALLOCATE(QATOM)
      DEALLOCATE(SATOM)
      DEALLOCATE(LINK)
      DEALLOCATE(MBOND)
      DEALLOCATE(SBOND)
                          CALL TRACE$POP
      RETURN
      CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
        SUBROUTINE STRCIN_SOLVENT_SETM(NAT,ATOM,NBOND,BOND)  
!       ****************************************************************    
!       **                                                            **    
!       ****************************************************************    
        IMPLICIT NONE
        INTEGER(4)     ,INTENT(IN) :: NAT
        INTEGER(4)     ,INTENT(IN) :: NBOND
        TYPE(ATOM_TYPE),INTENT(IN) :: ATOM(NAT)
        TYPE(BOND_TYPE),INTENT(IN) :: BOND(NBOND)
        INTEGER(4)                 :: IAT,IBOND
        CHARACTER(32)              :: ATOMNAME(NAT)
        CHARACTER(2)               :: ELEMENT(NAT)
        REAL(8)                    :: R(3,NAT)
        REAL(8)                    :: MASS(NAT)
        REAL(8)                    :: CHARGE(NAT)
        CHARACTER(5)               :: FFTYPE(NAT)
        REAL(8)                    :: BO(NBOND)
        INTEGER(4)                 :: I5VEC(5,NBOND)
!       ****************************************************************    
!
!       ================================================================
!       ==  SEND ATOM AND BOND INFORMATION TO CLASSICAL OBJECT        ==
!       ================================================================
        DO IAT=1,NAT
          R(:,IAT)   =ATOM(IAT)%R(:)
          MASS(IAT)  =ATOM(IAT)%M
          CHARGE(IAT)=ATOM(IAT)%Q
          FFTYPE(IAT)=ATOM(IAT)%FFTYPE
          ATOMNAME(IAT)=ATOM(IAT)%NAME
          ELEMENT(IAT)=ATOM(IAT)%ELEMENT
        ENDDO
        CALL CLASSICAL$SETI4('NAT',NAT)
        CALL CLASSICAL$SETR8A('R(0)',3*NAT,R)
        CALL CLASSICAL$SETR8A('R(-)',3*NAT,R)
        CALL CLASSICAL$SETR8A('MASS',NAT,MASS)
        CALL CLASSICAL$SETR8A('QEL',NAT,CHARGE)
        CALL CLASSICAL$SETCHA('TYPE',NAT,FFTYPE)
        CALL CLASSICAL$SETCHA('ATOMNAME',NAT,ATOMNAME)
        CALL CLASSICAL$SETCHA('ELEMENT',NAT,ELEMENT)
!
!       ================================================================
!       ==  SEND BOND INFORMATION TO CLASSICAL OBJECT                 ==
!       ================================================================
        DO IBOND=1,NBOND
          I5VEC(1,IBOND)=BOND(IBOND)%ATOM1
          I5VEC(2,IBOND)=BOND(IBOND)%ATOM2
          I5VEC(3:5,IBOND)=BOND(IBOND)%IT(:)
          BO(IBOND)=BOND(IBOND)%BO
        ENDDO
        CALL CLASSICAL$SETI4('NBOND',NBOND)
        CALL CLASSICAL$SETI4A('BOND',5*NBOND,I5VEC)
        CALL CLASSICAL$SETR8A('BONDORDER',NBOND,BO)
        CALL CLASSICAL$SETL4('LONGRANGE',.TRUE.)
        RETURN
        END SUBROUTINE STRCIN_SOLVENT_SETM
      END SUBROUTINE STRCIN_SOLVENT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_RESOLVEEXTENDEDNAME(XNAME,NAME,IT)
!     **************************************************************************
!     **  RESOLVES THE EXTENDED ATOM NAME NOTATION, WHICH INCLUDES            **
!     **  A LATTICE TRANSLATION                                               **
!     **                                                                      **
!     **  THE EXTENDED NOTATION INCLUDES AN INTEGER LATTICE TRANSLATIONS      **
!     **  IN THE ATOM NAME FOLLOWING A COLON                                  **
!     **                                                                      **
!     **   'O_23:+1-1+1'  ATOM 'O_23' SHIFTED BY RBAS(:,1)-RBAS(:,2)+RBAS(:,3)**
!     **                                                                      **
!     **   THE '+'SIGNS ARE NOT REQUIRED.                                     **
!     **   ONLY SINGLE-DIGIT TRANSLATIONS ARE PERMITTED                       **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: XNAME  ! EXTENDED ATOM NAME
      CHARACTER(*),INTENT(OUT):: NAME   ! NON-EXTENDED ATOM NAME
      INTEGER(4)  ,INTENT(OUT):: IT(3)  ! INTEGER LATTICE TRANSLATIONS
      INTEGER(4)              :: ICOLON ! POSITION OF THE COLON IN XNAME
      INTEGER(4)              :: IPOS,IND,SGN
      INTEGER(4)              :: ICH    ! ASCII NUMBER OF THE SELECTED LETTER
!     **************************************************************************
      ICOLON=INDEX(XNAME,':')
!     == RETURN IF NO TRANSLATION VECTOR GIVEN =================================
      IF(ICOLON.EQ.0) THEN
        NAME=XNAME
        IT(:)=0
        RETURN
      END IF
!
!     ==========================================================================
!     == RESOLVE EXTENDED ATOM NAME                                           ==
!     ==========================================================================
      NAME=XNAME(:ICOLON-1)
      IPOS=ICOLON+1
      IND=0
      SGN=+1
      DO WHILE(IND.LT.3) 
        ICH=IACHAR(XNAME(IPOS:IPOS))
!       ==  IACHAR('+')=43; IACHAR('-')=45; IACHAR('0')=48; IACHAR('1')=49;...
        IF(ICH.GE.48.AND.ICH.LE.57) THEN
          IND=IND+1
          IT(IND)=SGN*(ICH-48)
          SGN=+1
        ELSE IF(ICH.EQ.43) THEN
          SGN=+1
        ELSE IF(ICH.EQ.45) THEN
          SGN=-1
        ELSE
          CALL ERROR$MSG('ILLEGAL CHARACTER IN EXTENDED ATOM NOTATION')  
          CALL ERROR$CHVAL('EXT. NAME ',XNAME)
          CALL ERROR$CHVAL('ILLEGAL CHARACTER ',XNAME(IPOS:IPOS))
          CALL ERROR$STOP('STRCIN_RESOLVEEXTENDEDNAME')
        END IF
        IPOS=IPOS+1
      ENDDO
      IF(XNAME(IPOS:).NE.' ') THEN
        CALL ERROR$MSG('LETTERS FOUND BEYOND END OF EXTENDED ATOM NOTATION')  
        CALL ERROR$CHVAL('EXT. NAME ',XNAME)
        CALL ERROR$CHVAL('ADDITIONAL LETTERS ',XNAME(IPOS:))
        CALL ERROR$STOP('STRCIN_RESOLVEEXTENDEDNAME')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCIN_COSMO(LL_STRC_)
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)            :: LL_STRC
      LOGICAL(4)               :: TCHK
      LOGICAL(4)               :: TISO
      INTEGER(4)               :: NAT1
      INTEGER(4)               :: IAT1
      INTEGER(4)               :: IAT,NAT
      CHARACTER(32)            :: NAME
      CHARACTER(250)           :: TESSFILE
      REAL(8)      ,ALLOCATABLE:: RSOLV(:)
      REAL(8)                  :: SVAR
!     ******************************************************************
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'COSMO',1,TCHK)
      IF(.NOT.TCHK) RETURN
                                 CALL TRACE$PUSH('STRCIN_COSMO')
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(RSOLV(NAT))
      RSOLV(:)=0.D0
!
!     == CHECK WITH ISOLATE IF CALCULATION IS PERIODIC OR ISOLATED
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'ISOLATE',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_STRC,'ISOLATE')
        CALL LINKEDLIST$EXISTD(LL_STRC,'DECOUPLE',0,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'DECOUPLE',1,TISO)
        ELSE
          TISO=.TRUE.
        END IF
      ELSE
        TISO=.FALSE.
      END IF
      CALL COSMO$SETL4('PERIODIC',.NOT.TISO)
!
!     == NOW BACK TO COSMO BLOCK
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'COSMO')
!
!     == EPSILON = DIELECTRIC CONSTANT ===========================
      CALL LINKEDLIST$EXISTD(LL_STRC,'EPSILON',1,TCHK)
      IF(.NOT.TCHK) THEN ! USE METALLIC SCREENING AS DEFAULT
        CALL LINKEDLIST$SET(LL_STRC,'EPSILON',0,1.D12)
      END IF
      CALL LINKEDLIST$GET(LL_STRC,'EPSILON',1,SVAR)
      CALL COSMO$SETR8('EPSILON',SVAR)
!
!     ==  VPAULI ===================================================
      CALL LINKEDLIST$EXISTD(LL_STRC,'VPAULI',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_STRC,'VPAULI',0,0.D0)
      END IF
      CALL LINKEDLIST$GET(LL_STRC,'VPAULI',1,SVAR)
      CALL COSMO$SETR8('VPAULI',SVAR)
!
!     ==  TESSELATION FILE=====================================
      CALL LINKEDLIST$EXISTD(LL_STRC,'TESSFILE',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_STRC,'TESSFILE',1,TESSFILE)
        CALL FILEHANDLER$SETFILE(+'TESSFILE',.FALSE.,TRIM(TESSFILE))
        CALL FILEHANDLER$SETSPECIFICATION(+'TESSFILE','STATUS','OLD')
        CALL FILEHANDLER$SETSPECIFICATION(+'TESSFILE','POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(+'TESSFILE','ACTION','READ')
        CALL FILEHANDLER$SETSPECIFICATION(+'TESSFILE','FORM','FORMATTED')
        CALL COSMO$READTESSELATION('TESSFILE')
        CALL FILEHANDLER$CLOSE('TESSFILE')
      END IF
!
      CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NAT1)
      IF(NAT1.NE.NAT) THEN
        CALL ERROR$MSG('NUMBER OF ATOMS INCONSISTENT')
        CALL ERROR$STOP('STRCIN_COSMO')
      END IF
      DO IAT1=1,NAT1
        CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT1)
!
!       ==  NAME  ================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'NAME',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!STRUCTURE!COSMO!ATOM:NAME NOT SPECIFIED')
          CALL ERROR$STOP('STRCIN_COSMO')
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'NAME',1,NAME)
        CALL ATOMLIST$INDEX(NAME,IAT)
!
!       ==  RAD  ================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'RAD',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL LINKEDLIST$SET(LL_STRC,'RAD',1,3.D0)
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'RAD',1,RSOLV(IAT))
!
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
      CALL COSMO$SETR8A('RSOLV',NAT,RSOLV)
      DEALLOCATE(RSOLV)
                                          CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE STRCOUT
!     ******************************************************************
!     **                                                              **
!     **  UPDATES THE BUFFER STRC AND WRITES IT ON STRC_OUT           **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      USE IO_MODULE, ONLY : LL_STRC
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      INTEGER(4)               :: NFIL
      INTEGER(4)               :: NAT
      INTEGER(4)               :: NATMM
      REAL(8)                  :: R(3)
      REAL(8)                  :: Q
      REAL(8)   ,ALLOCATABLE   :: RMM(:,:)
      REAL(8)                  :: RUNIT
      REAL(8)                  :: RBAS(3,3)
      INTEGER(4)               :: IAT           ! RUNNING VARIABLES
      INTEGER(4)               :: IAT1          ! AUXILARY VARIABLES
      CHARACTER(32)            :: NAME
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: NLINK,ILINK
      INTEGER(4)               :: NTASKS,THISTASK
!     ******************************************************************
                          CALL TRACE$PUSH('STRCOUT')
!    
!     ==================================================================
!     ==  UPDATE BUFFER STRC                                          ==
!     ==================================================================
!     __PULL OUT LENGTH UNIT OF THE STRUCTURE FILE______________________
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')
      CALL LINKEDLIST$GET(LL_STRC,'LUNIT',1,RUNIT)
!
!     == LATTICE =======================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'LATTICE')
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL LINKEDLIST$SET(LL_STRC,'T',0,RBAS/RUNIT)

!     == ATOMS =========================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NAT)
      DO IAT=1,NAT
        CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT)
        CALL LINKEDLIST$GET(LL_STRC,'NAME',1,NAME)
        CALL ATOMLIST$INDEX(NAME,IAT1)
        CALL ATOMLIST$GETR8A('R(0)',IAT1,3,R)
        CALL LINKEDLIST$SET(LL_STRC,'R',0,R/RUNIT)
        CALL ATOMLIST$GETR8('Q',IAT1,Q)
        CALL LINKEDLIST$SET(LL_STRC,'Q',0,-Q)
        CALL LINKEDLIST$SET(LL_STRC,'INDEX',0,IAT1)
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO
!    
!     ==================================================================
!     ==  WRITE BUFFER STRC                                           ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$EXISTL(LL_STRC,'QM-MM',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_STRC,'QM-MM')
!       == REPLACE MM ATOMIC POSITIONS AND CHARGES
        CALL CLASSICAL$SELECT('QMMM')
        CALL CLASSICAL$GETI4('NAT',NATMM)
        ALLOCATE(RMM(3,NATMM))
        CALL CLASSICAL$GETR8A('R(0)',3*NATMM,RMM)
        CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NAT)
        IF(NAT.NE.NATMM) THEN
        END IF
        DO IAT=1,NAT
          CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',IAT)
          CALL LINKEDLIST$EXISTD(LL_STRC,'R',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$SET(LL_STRC,'R',1,RMM(:,IAT)/RUNIT)
          END IF
          CALL LINKEDLIST$SELECT(LL_STRC,'..')
        ENDDO  
        DEALLOCATE(RMM)
!
!       == REPLACE SHADOW ATOMIC POSITION IN THE LINKS
        CALL CLASSICAL$SELECT('SHADOW')
        CALL CLASSICAL$GETI4('NAT',NATMM)
        ALLOCATE(RMM(3,NATMM))
        CALL CLASSICAL$GETR8A('R(0)',3*NATMM,RMM)
        CALL LINKEDLIST$NLISTS(LL_STRC,'LINK',NLINK)
        DO ILINK=1,NLINK
          CALL LINKEDLIST$SELECT(LL_STRC,'LINK',ILINK)
          CALL LINKEDLIST$SELECT(LL_STRC,'SHADOW')
          CALL LINKEDLIST$SET(LL_STRC,'R',0,RMM(:,NATMM-NLINK+ILINK)/RUNIT)
          CALL LINKEDLIST$SELECT(LL_STRC,'..')
          CALL LINKEDLIST$SELECT(LL_STRC,'..')
        ENDDO
        DEALLOCATE(RMM)
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      END IF
!    
!     ==================================================================
!     ==  WRITE BUFFER STRC                                           ==
!     ==================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$UNIT('STRC_OUT',NFIL)
        CALL LINKEDLIST$SELECT(LL_STRC,'~')
        CALL LINKEDLIST$WRITE(LL_STRC,NFIL,'MONOMER')
        CALL FILEHANDLER$CLOSE('STRC_OUT')
      END IF
!
                                      CALL TRACE$POP
      RETURN
    END SUBROUTINE STRCOUT

