!........1.........2.........3.........4.........5.........6.........7.........8
MODULE TRAJECTORY_MODULE
TYPE TRA_TYPE
  INTEGER(4)             :: NSTEP
!CLEMENS
  REAL(8)   ,POINTER     :: CELL(:,:)!(9,NSTEP)
  REAL(8)   ,POINTER     :: R(:,:,:) !(3,NAT,NSTEP)  POSITION
  REAL(8)   ,POINTER     :: Q(:,:)   !(NAT,NSTEP)    CHARGE
  REAL(8)   ,POINTER     :: T(:)     !(NSTEP)        TIME
  REAL(8)   ,POINTER     :: CHARGE(:,:) !(NAT,NSTEP)  #(ELECTRONS IN SPHERE)
  REAL(8)   ,POINTER     :: SPIN(:,:,:) !(3,NAT,NSTEP)  SPIN MOMENT IN SPHERE
  INTEGER(4),POINTER     :: ISTEP(:) !(NSTEP)        TIME STEP COUNTER
END TYPE TRA_TYPE
TYPE MODE_TYPE
  CHARACTER(32)          :: ID
  REAL(8)      ,POINTER  :: VAL(:)
END TYPE MODE_TYPE
TYPE BOND_TYPE
  INTEGER(4) :: IAT1    ! CENTRAL ATOM
  INTEGER(4) :: IAT2    ! NEIGHBOR AT R(IAT2)=RBAS*IT
  INTEGER(4) :: IT(3)   ! LATTICE TRANSLATIONS
  REAL(8)    :: DIS     ! DISTANCE/ SUM OF COVALENT RADII
END TYPE BOND_TYPE
TYPE(TRA_TYPE)              :: TRA
CHARACTER(32)  ,ALLOCATABLE :: ATOM(:)
INTEGER(4)     ,ALLOCATABLE :: IZ(:)
REAL(8)        ,ALLOCATABLE :: MASS(:)
TYPE(MODE_TYPE),ALLOCATABLE :: MODE(:)
INTEGER(4)                  :: NMODES
INTEGER(4)                  :: NSTEP
INTEGER(4)                  :: NAT,QNAT
REAL(8)                     :: RBAS(3,3)
LOGICAL(4)                  :: TQMMM=.FALSE.
CHARACTER(16)               :: FFORMAT
CONTAINS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORY_NEWMODES
      IMPLICIT NONE
      INTEGER(4)    :: IMODE
!     **************************************************************************
      ALLOCATE(MODE(NMODES))
      NSTEP=TRA%NSTEP
      DO IMODE=1,NMODES
        MODE(IMODE)%ID=' '
        ALLOCATE(MODE(IMODE)%VAL(NSTEP))
      ENDDO
      RETURN
      END SUBROUTINE TRAJECTORY_NEWMODES
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORY_ATOMLOOKUP(NAME,IAT)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(OUT):: IAT
      INTEGER(4)              :: I
!     **************************************************************************
      IF(.NOT.TQMMM) THEN
         DO I=1,NAT
            IF(NAME.EQ.ATOM(I)) THEN
               IAT=I
               RETURN
            END IF
         ENDDO
      ELSE
         DO I=1,QNAT
            IF(NAME.EQ.ATOM(I)) THEN
               IAT=I
               RETURN
            END IF
         ENDDO
      END IF
      CALL ERROR$MSG('ATOM NAME NOT RECOGNIZED')
      CALL ERROR$CHVAL('NAME ',NAME)
      CALL ERROR$STOP('TRAJECTORY_ATOMLOOKUP')
      RETURN
      END SUBROUTINE TRAJECTORY_ATOMLOOKUP
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRAJECTORY_MODELOOKUP(NAME,IMODE)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(OUT):: IMODE
      INTEGER(4)              :: I
!     **************************************************************************
      DO I=1,NMODES
        IF(NAME.EQ.MODE(I)%ID) THEN
          IMODE=I
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('MODE NAME NOT RECOGNIZED')
      CALL ERROR$CHVAL('NAME ',NAME)
      CALL ERROR$STOP('TRAJECTORY_MODELOOKUP')
      RETURN
      END SUBROUTINE TRAJECTORY_MODELOOKUP
END MODULE TRAJECTORY_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM MAIN
!     **************************************************************************
!     **                                                                      **
!     **  MAIN ROUTINE OF TRAJECTORY TOOL                                     **
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      USE TRAJECTORY_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)  :: LL_CNTL
      TYPE(LL_TYPE)  :: LL_STRC
      CHARACTER(256) :: FILE
      CHARACTER(256) :: ROOT
      INTEGER(4)     :: I,NLISTS
      INTEGER(4)     :: NFIL
      INTEGER(4)     :: NFILO
      LOGICAL(4)     :: TCHK,TFILE
      CHARACTER(32)  :: ID
      REAL(8)        :: T1,T2,DT
      REAL(8)        :: SECOND
      REAL(8)        :: PICO
      REAL(8)        :: FEMTO
      CHARACTER(40)  :: STR
!     **************************************************************************
                          CALL TRACE$PUSH('MAIN')
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT

      CALL GET_COMMAND_ARGUMENT(1,FILE)
      I=INDEX(FILE,'.',BACK=.TRUE.)
      ROOT=FILE(1:I-1)
!
!     ==========================================================================
!     ==  DEFINE FILES                                                        ==
!     ==========================================================================
      CALL FILEHANDLER$SETROOT(ROOT)
      CALL FILEHANDLER$SETFILE('CNTL',.FALSE.,TRIM(FILE))
      CALL FILEHANDLER$SETSPECIFICATION('CNTL','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('CNTL','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('CNTL','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('CNTL','FORM','FORMATTED')
!
!     ==========================================================================
!     ==  READ CNTL FILE TO LINKEDLIST                                        ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('CNTL',NFIL)
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
!    
!     ==========================================================================
!     ==  MARK ALL ELEMENTS AS READ FROM INPUT FILE                           ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$MARK(LL_CNTL,1)
!
!     ==========================================================================
!     ==  RESET FILE NAME FOR STRUCTURE FILE IF REQUESTED                     ==
!     ==========================================================================
      CALL TRACE$PASS('RESET FILENAME FOR STRC FILE')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'FILES')
!
!     == ROOT NAME =============================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ROOT',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'ROOT',0,ROOT)
      CALL LINKEDLIST$GET(LL_CNTL,'ROOT',1,ROOT)
      CALL FILEHANDLER$SETROOT(ROOT)
!
!     ==========================================================================
!     ==  SET STANDARD FILE NAMES                                             ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('PROT',.TRUE.,-'.TPROT')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','FORM','FORMATTED')
      CALL FILEHANDLER$SETFILE('STRC',.TRUE.,-'.STRC_OUT')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('STRC','FORM','FORMATTED')
!
!     ==========================================================================
!SASCHA QM-MM
!     ==========================================================================
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')

      CALL LINKEDLIST$EXISTD(LL_CNTL,'QMMM',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'QMMM',0,TQMMM)
      CALL LINKEDLIST$GET(LL_CNTL,'QMMM',1,TCHK)
      IF(TCHK) THEN
         CALL FILEHANDLER$SETFILE('TRA',.TRUE.,-'_R.QMMMTRA')
         TQMMM=.TRUE.
         CALL LINKEDLIST$EXISTD(LL_CNTL,'FORMAT',1,TFILE)
         IF(.NOT.TFILE) CALL LINKEDLIST$SET(LL_CNTL,'FORMAT',0,'STRC')
         CALL LINKEDLIST$GET(LL_CNTL,'FORMAT',1,FFORMAT)
         IF(FFORMAT.EQ.'PDB') THEN
            CALL FILEHANDLER$SETFILE('MMSTRC',.TRUE.,-'.PDB')
            CALL FILEHANDLER$SETSPECIFICATION('MMSTRC','STATUS','OLD')
            CALL FILEHANDLER$SETSPECIFICATION('MMSTRC','POSITION','REWIND')
            CALL FILEHANDLER$SETSPECIFICATION('MMSTRC','ACTION','READ')
            CALL FILEHANDLER$SETSPECIFICATION('MMSTRC','FORM','FORMATTED')
         END IF
      ELSE
         CALL FILEHANDLER$SETFILE('TRA',.TRUE.,-'_R.TRA')
      END IF
      CALL FILEHANDLER$SETSPECIFICATION('TRA','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('TRA','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('TRA','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('TRA','FORM','UNFORMATTED')

      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'FILES')
!
!     ==========================================================================
!     ==  RENAME FILES IF DESIRED                                             ==
!     ==========================================================================
      CALL LINKEDLIST$NLISTS(LL_CNTL,'FILE',NLISTS)
      DO I=1,NLISTS
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',I)
!
!       == GET FILE ID =========================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('KEYWORD ID IS MANDATORY')
          CALL ERROR$STOP('MAIN')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,ID)
!
!       == GET FILE NAME =======================================================
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILE)
!
!       == GET EXTENSION FLAG ==================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
        END IF
        CALL FILEHANDLER$SETFILE(ID,TCHK,FILE)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
      CALL TRACE$PASS('FILENAME FOR STRC FILE SET')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
!
!     ==========================================================================
!     == HEADER                                                               ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,FMT='()')
      WRITE(NFILO,FMT='(72("*"))')
      WRITE(NFILO,FMT='(72("*"),T15 &
     &           ,"          TRAJECTORY ANALYSIS                ")')
      WRITE(NFILO,FMT='(72("*"),T15 &
     &           ,"   FOR THE PROJECTOR-AUGMENTED WAVE METHOD   ")')
      WRITE(NFILO,FMT='(72("*"))')
      WRITE(NFILO,FMT='(T10 &
     &         ,"P.E. BLOECHL, CLAUSTHAL UNIVERSITY OF TECHNOLOGY")')
      WRITE(NFILO,FMT='(T10 &
     &            ,"DISTRIBUTED UNDER THE GNU PUBLIC LICENSE V3")')
      WRITE(NFILO,*)
!
!     ==========================================================================
!     ==  READ STRC FILE TO LINKEDLIST                                        ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('STRC',NFIL)
      CALL LINKEDLIST$NEW(LL_STRC)
      CALL LINKEDLIST$READ(LL_STRC,NFIL,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
!    
!     ==========================================================================
!     ==  MARK ALL ELEMENTS AS READ FROM INPUT FILE                           ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$MARK(LL_STRC,1)
!
!     ==========================================================================
!     == GET ATOMIC NUMBER                                                    ==
!     ==========================================================================
      CALL READ_STRC(LL_STRC)
!
!     ==========================================================================
!     ==  READ TIME INTERVAL FOR ANALYSIS                                     ==
!     ==========================================================================
      CALL CONSTANTS$GET('SECOND',SECOND)
      CALL CONSTANTS$GET('FEMTO',FEMTO)
      CALL CONSTANTS$GET('PICO',PICO)
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'SEQUENCE')
!
!     == GET BEGINNING OF TIME SEQUENCE ========================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'T1[PS]',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'T1[PS]',0,-1.D+9)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'T1[PS]',1,T1)
      T1=T1*PICO*SECOND
!
!     == GET END OF TIME SEQUENCE ========================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'T2[PS]',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'T2[PS]',0,+1.D+9)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'T2[PS]',1,T2)
      T2=T2*PICO*SECOND
!
!     == GET STEP OF TIME SEQUENCE ========================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'DT[FS]',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SET(LL_CNTL,'DT[FS]',0,0.D0)
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'DT[FS]',1,DT)
      DT=DT*FEMTO*SECOND
      CALL LINKEDLIST$SELECT(LL_CNTL,'..')
!
!     ==========================================================================
!     ==  READ TRAJECTORY FILE                                                ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('TRA',NFIL)      
      CALL READTRA(NFIL,T1,T2,DT)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      CALL CLUSTERIT(LL_CNTL)
!
!     ==========================================================================
!     ==  PRODUCE OUTPUT                                                      ==
!     ==========================================================================
      CALL MKTRA(LL_CNTL,LL_STRC)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      CALL FILEHANDLER$REPORT(NFILO,'USED')
!      
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$REPORT_UNUSED(LL_CNTL,NFILO)
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$REPORT_UNUSED(LL_STRC,NFILO)
!
!     ==========================================================================
!     ==  CLOSING                                                             ==
!     ==========================================================================
      WRITE(NFILO,FMT='(72("="))')
      WRITE(NFILO,FMT='(72("="),T20," PAW_TRA TOOL FINISHED ")')
      WRITE(NFILO,FMT='(72("="))')
                          CALL TRACE$POP
      CALL ERROR$NORMALSTOP()
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MKTRA(LL_CNTL_,LL_STRC_)
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE), INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE), INTENT(IN) :: LL_STRC_
      TYPE(LL_TYPE)             :: LL_CNTL
      TYPE(LL_TYPE)             :: LL_STRC
      INTEGER(4)                :: NFIL
      LOGICAL(4)                :: TCHK
!     **************************************************************************
                          CALL TRACE$PUSH('MKTRA')
      LL_CNTL=LL_CNTL_
      LL_STRC=LL_STRC_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
!
!     ==========================================================================
!     ==  MAKE A MOVIE FILE                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTL(LL_CNTL,'MOVIE',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'MOVIE')
        CALL WRITETRA(LL_STRC,LL_CNTL)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      END IF    
!
!     ==========================================================================
!     ==  CALCULATE AVERAGE TEMPERATURE OF THE ATOMS                          ==
!     ==========================================================================
      CALL TEMPERATURE(LL_CNTL)
!
!     ==========================================================================
!     ==  CREATE SPAGHETTI PLOT                                               ==
!     ==========================================================================
       CALL SPAGHETTI(LL_CNTL)
!
!     ==========================================================================
!     ==  DEFINE MODES                                                        ==
!     ==========================================================================
      CALL MODES(LL_CNTL)
!
!     ==========================================================================
!     ==   PRINT MODES                                                        ==
!     ==========================================================================
      CALL WRITEMODES(LL_CNTL)
!
!     ==========================================================================
!     ==   PRINT NEIGHBORLIST                                                 ==
!     ==========================================================================
      CALL NEIGHBORS(LL_CNTL)
!
!     ==========================================================================
!     ==   WRITE ATOMIC POSITIONS FOR A GIVEN TIME STEP                       ==
!     ==========================================================================
      CALL SNAPSHOT(LL_CNTL)
!
!     ==========================================================================
!     ==   DISTANCE AND ANGLE CORRELATIONS                                    ==
!     ==========================================================================
      CALL CORRELATION(LL_CNTL)
!
!     ==========================================================================
!     ==   PRINT MASS WEIGHTED PATH LENGTH                                    ==
!     ==========================================================================
      CALL SOFT(LL_CNTL)
!
!     ==========================================================================
!     ==   UNIT CELL MATRIX AS FUNCTION OF TIME                               ==
!     ==========================================================================
      CALL CELL(LL_CNTL)
                                 CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CORRELATION(LL_CNTL_)
!     **************************************************************************
!     **                                                                      **
!     **  DISTANCES AND ANGLES AS FUNCTION OF TIME                            **
!     **  OR DENSITY AS FUNCTION OF DISTANCE OR ANGLE                         **
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      USE SORT_MODULE
      USE TRAJECTORY_MODULE, ONLY : TQMMM &
     &                             ,BOND_TYPE &
     &                             ,TRA_TYPE &
     &                             ,RBAS &
     &                             ,TRA &
     &                             ,ATOM &
     &                             ,NSTEP &
     &                             ,NAT &
     &                             ,QNAT   !#(ATOMS IN QM-MM) 
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      INTEGER(4)               :: NPLOT
      LOGICAL(4)               :: TCHK,TCHK1,TCHK2
      CHARACTER(1)             :: DA
      CHARACTER(16)            :: DEP
      INTEGER(4)               :: I,IAT,ISTEP,IATM
      INTEGER(4)               :: NFIL
      INTEGER(4)               :: NFILO
      CHARACTER(256)           :: FILE
      LOGICAL(4)   ,ALLOCATABLE:: TSELECTC(:)
      LOGICAL(4)   ,ALLOCATABLE:: TSELECTP(:)
      REAL(8)                  :: PICO,SECOND
      REAL(8)                  :: TIME
      REAL(8)                  :: T1,T2,T3
      REAL(8)                  :: DR(3)
      REAL(8)                  :: RMAX       ! DISTANCE CRITERION FOR BONDS
      REAL(8)                  :: SCALEX     ! 
      INTEGER(4)   ,PARAMETER  :: NBX=1000    ! X#(BONDS) -ARRAY DIMENSION
      INTEGER(4)               :: NB         ! ACTUAL NUMBER OF BONDS
      TYPE(BOND_TYPE)          :: BOND(NBX)
      INTEGER(4)   ,PARAMETER  :: INNX=1000    ! X#(BONDS) -ARRAY DIMENSION
      INTEGER(4)               :: INNMAX     ! ACTUAL NUMBER OF BONDS
      REAL(8)                  :: FIOFT(INNX,NSTEP)
      TYPE(BOND_TYPE)          :: BONDS(NBX) ! BOND LIST
      INTEGER(4)               :: NATM
      INTEGER(4)               :: IAT1,IAT2,IB,IX,INN
      REAL(8)      ,ALLOCATABLE:: POSM(:,:)
      REAL(8)      ,ALLOCATABLE:: RADM(:)
      LOGICAL(4)   ,ALLOCATABLE:: TSELECTCM(:)
      LOGICAL(4)   ,ALLOCATABLE:: TSELECTPM(:)
      REAL(8)                  :: SVAR,FAC,WEIGHT
      REAL(8)                  :: TM,T0,TP
      REAL(8)                  :: XMIN,XMAX,DX
      INTEGER(4)   ,PARAMETER  :: NX=1000
      REAL(8)      ,ALLOCATABLE:: P(:)
      INTEGER(4)               :: TO,FROM
      INTEGER(4)               :: NATOM
!     **************************************************************************
                               CALL TRACE$PUSH('CORRELATION')
      CALL CONSTANTS$GET('PICO',PICO)
      CALL CONSTANTS$GET('SECOND',SECOND)
      CALL FILEHANDLER$UNIT('PROT',NFILO)

      IF(TQMMM) THEN
         ALLOCATE(TSELECTP(QNAT))
         ALLOCATE(TSELECTC(QNAT))
         NATOM=QNAT
      ELSE
         ALLOCATE(TSELECTP(NAT))
         ALLOCATE(TSELECTC(NAT))
         NATOM=NAT
      END IF

      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'CORRELATION',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL TRACE$POP
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'CORRELATION')
      CALL REPORT$TITLE(NFILO,'CORRELATION')
!     
!     ==========================================================================
!     ==  SPECIFY FILE                                                        ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTL(LL_CNTL,'FILE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',0)
        FILE=-'.TRA.CORRELATION'
        CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,TRIM(FILE))
        TCHK=.TRUE.
        CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,TCHK)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ELSE
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,-'.TRA.CORRELATION')
          CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.TRUE.)
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILE)
!     
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
        CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      END IF
      CALL FILEHANDLER$SETFILE('CORRELATION',TCHK,TRIM(FILE))
      CALL FILEHANDLER$SETSPECIFICATION('CORRELATION','FORM','FORMATTED')
!     
!     ==========================================================================
!     ==  SELECT CENTRAL AND NEIGHBOR ATOMS                                   ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTL(LL_CNTL,'CENTER',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'CENTER',0)
        CALL SELECTION(LL_CNTL,NATOM,ATOM,TSELECTC)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ELSE
        TSELECTC(:)=.TRUE.
      END IF
      CALL LINKEDLIST$EXISTL(LL_CNTL,'PARTNER',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'PARTNER',0)
        CALL SELECTION(LL_CNTL,NATOM,ATOM,TSELECTP)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ELSE
        TSELECTP(:)=.TRUE.
      END IF
!     
!     ==========================================================================
!     == SELECT HISTOGRAM OR SPAGHETTI                                        ==
!     ========================================== ===============================
!     == OPTION SPAGHETTI IS NOT DOCUMENTED IN THE MANUAL.
!     == ITS IMPLEMENTATION IS NOT FUNCTIONAL
      CALL LINKEDLIST$EXISTD(LL_CNTL,'DEP',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'DEP',1,DEP)
        DEP=+DEP
        IF(DEP.NE.'SPAGHETTI'.AND.DEP.NE.'HISTOGRAM') THEN
          CALL ERROR$MSG('DEP CAN ONLY BE "SPAGHETTI" OR "HISTOGRAM"')
          CALL ERROR$CHVAL('DEP',DEP)
          CALL ERROR$STOP('CORRELATION')
        END IF
      ELSE
        DEP='HISTOGRAM'
      END IF
!     
!     ==========================================================================
!     == DEFINE POSITIONS TO BE CONSIDERED                                    ==
!     ==========================================================================
      NATM=0
      DO IAT=1,NATOM
        IF(.NOT.(TSELECTC(IAT).OR.TSELECTP(IAT))) CYCLE
        NATM=NATM+1
      ENDDO
      ALLOCATE(POSM(3,NATM))
      ALLOCATE(RADM(NATM))
      ALLOCATE(TSELECTCM(NATM))
      ALLOCATE(TSELECTPM(NATM))
      IATM=0
      DO IAT=1,NATOM
        IF(.NOT.(TSELECTC(IAT).OR.TSELECTP(IAT))) CYCLE
        IATM=IATM+1
        TSELECTCM(IATM)=TSELECTC(IAT)
        TSELECTPM(IATM)=TSELECTP(IAT)
      ENDDO
!     
!     ==========================================================================
!     == BRANCH OUT TO PAIR CORRELATION FUNCTION                              ==
!     ==========================================================================
      IF(DEP.EQ.'HISTOGRAM') THEN
        CALL FILEHANDLER$UNIT('CORRELATION',NFIL)
! FAST
        CALL PAIRCORRELATION1(NFIL,RBAS,NATOM,TSELECTC,TSELECTP,NSTEP,TRA%R)
! SECURE
!        CALL PAIRCORRELATION2(NFIL,RBAS,NATOM,TSELECTC,TSELECTP,NSTEP,TRA%R)
        CALL FILEHANDLER$CLOSE('CORRELATION')
        CALL REPORT$TITLE(NFILO,'CORRELATION')
        RETURN
      END IF
!     
!     ==========================================================================
!     == DEFINE DISTANCE CRITERION                                            ==
!     ==========================================================================
      RMAX=5.0
      RADM(:)=0.5D0*RMAX
      SCALEX=1.0D0
!     
!     ==========================================================================
!     == COLLECT TRAJECTORIES                                                 ==
!     ==========================================================================
      INNMAX=0
      DO ISTEP=1,NSTEP
        IATM=0
        DO IAT=1,NATOM
          IF(.NOT.(TSELECTC(IAT).OR.TSELECTP(IAT))) CYCLE
          IATM=IATM+1
          POSM(:,IATM)=TRA%R(:,IAT,ISTEP)
        ENDDO
        CALL BONDS_EVAL(RBAS,NATM,POSM,RADM,SCALEX,NBX,NB,BOND)
        FIOFT(:,ISTEP)=RMAX
        INN=0
        DO IB=1,NB
          IAT1=BOND(IB)%IAT1
          IAT2=BOND(IB)%IAT2
          TCHK1=TSELECTCM(IAT1).AND.TSELECTPM(IAT2)
          TCHK2=TSELECTCM(IAT2).AND.TSELECTPM(IAT1)
          IF(.NOT.(TCHK1.OR.TCHK2)) CYCLE
          T1=DBLE(BOND(IB)%IT(1))
          T2=DBLE(BOND(IB)%IT(2))
          T3=DBLE(BOND(IB)%IT(3))
          DR(:)=POSM(:,IAT2)-POSM(:,IAT1) &
     &         +RBAS(:,1)*T1+RBAS(:,2)*T2+RBAS(:,3)*T3
          INN=INN+1
          IF(INN.GT.INNX) THEN
            CALL ERROR$STOP('CORRELATION')
          END IF
          FIOFT(INN,ISTEP)=SQRT(DR(1)**2+DR(2)**2+DR(3)**2)
        ENDDO            
        INNMAX=MAX(INNMAX,INN)
      ENDDO   
!     
!     ==========================================================================
!     ==  ORDER DISTANCES                                                     ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('CORRELATION',NFIL)
      IF(DEP.EQ.'SPAGHETTI') THEN
        DO ISTEP=1,NSTEP
!         __ ORDER ARRAY________________________________________________________
          CALL SORT$SET(INNMAX,FIOFT(1,ISTEP))
          CALL SORT$RESTART
          CALL SORT$FLIP(FROM,TO)
          DO WHILE (FROM.NE.0.OR.TO.NE.0)
            IF(TO.EQ.0) THEN
              SVAR=FIOFT(FROM,ISTEP)
            ELSE IF(FROM.EQ.0) THEN
              FIOFT(TO,ISTEP)=SVAR
            ELSE
              FIOFT(TO,ISTEP)=FIOFT(FROM,ISTEP)
            END IF
            CALL SORT$FLIP(FROM,TO)
          ENDDO
          CALL SORT$UNSET
        ENDDO
!       == WRITE ===============================================================
        DO INN=1,INNMAX
          WRITE(NFIL,FMT='(2F10.5)')TRA%T(1)/(PICO*SECOND),RMAX
          DO ISTEP=1,NSTEP
            WRITE(NFIL,FMT='(2F10.5)')TRA%T(ISTEP)/(PICO*SECOND),FIOFT(INN,ISTEP)
          ENDDO
          WRITE(NFIL,FMT='(2F10.5)')TRA%T(NSTEP)/(PICO*SECOND),RMAX
        ENDDO
      ELSE IF(DEP.EQ.'HISTOGRAM') THEN
        XMAX=RMAX
        XMIN=0.D0
        DX=(XMAX-XMIN)/DBLE(NX-1)
        ALLOCATE(P(NX))
        TM=TRA%T(1)
        T0=TRA%T(1)
        SVAR=0.5D0/(TRA%T(NSTEP)-T0)/DX
        DO ISTEP=1,NSTEP-1
          TP=TRA%T(MIN(NSTEP,ISTEP+1))
          WEIGHT=SVAR*(TP-TM)
          DO INN=1,INNMAX
            IX=NINT((FIOFT(INN,ISTEP)-XMIN)/DX)+1
            IF(IX.LT.1) CYCLE
            IF(IX.GT.NX) CYCLE
            P(IX)=P(IX)+WEIGHT
          ENDDO
          TM=T0
          T0=TP
        ENDDO
        DO IX=1,NX
          WRITE(NFIL,*)XMIN+DBLE(IX-1)*DX,P(IX)
        ENDDO
        DEALLOCATE(P)
      END IF
      CALL FILEHANDLER$CLOSE('CORRELATION')
!     
!     ==========================================================================
!     ==  WRITE ATOMS TO PROTOCOL FILE                                        ==
!     ==========================================================================
!     
!     ==========================================================================
!     == CLOSE DOWN                                                           ==
!     ==========================================================================
      DEALLOCATE(POSM)
      DEALLOCATE(RADM)
      DEALLOCATE(TSELECTCM)
      DEALLOCATE(TSELECTPM)
                               CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAIRCORRELATION1(NFIL,RBAS,NAT,TSELECTC,TSELECTP,NSTEP,R)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE TRAJECTORY_MODULE, ONLY: BOND_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      INTEGER(4),INTENT(IN) :: NAT
      LOGICAL(4),INTENT(IN) :: TSELECTC(NAT)
      LOGICAL(4),INTENT(IN) :: TSELECTP(NAT)
      INTEGER(4),INTENT(IN) :: NSTEP
      REAL(8)   ,INTENT(IN) :: R(3,NAT,NSTEP)
      REAL(8)   ,PARAMETER  :: ANGSTROM=1.889726D0 ! ANGSTROM
      REAL(8)   ,PARAMETER  :: RMAX=5.D0*ANGSTROM ! OUTERMOST POINT OF HISTOGRAM
      REAL(8)   ,PARAMETER  :: RMAX2=RMAX**2
      REAL(8)   ,PARAMETER  :: SCALEX=2.D0  ! SELECT WITH RMAX*SCALEX
      INTEGER(4),PARAMETER  :: NBX=100000     ! #(BONDS CONSIDERED)
      INTEGER(4),PARAMETER  :: NP=1000  ! #(GRID POINTS ON HISTOGRAM)
      REAL(8)               :: HISTOGRAM(NP)
      INTEGER(4)            :: NB ! #(BONDS)
      REAL(8)   ,ALLOCATABLE:: RAD(:) !(NAT)
      TYPE(BOND_TYPE),ALLOCATABLE :: BOND(:)
      TYPE(BOND_TYPE)       :: BOND1
      REAL(8)               :: DR(3)
      REAL(8)               :: DIS2
      REAL(8)               :: SVAR ! SUPPORT VARIABLE
      INTEGER(4)            :: IATC,IATP,ISTEP,IP,I,IB
      INTEGER(4)            :: IT(3)
      INTEGER(4)            :: NATC  !#(CENTER ATOMS)
!     **************************************************************************
      ALLOCATE(RAD(NAT))
      RAD(:)=0.5D0*RMAX
      NATC=0
      DO IATC=1,NAT
        IF(TSELECTC(IATC))NATC=NATC+1
      ENDDO
      ALLOCATE(BOND(NBX))
      ISTEP=1
      NB=0
      CALL BONDS_EVAL(RBAS,NAT,R(:,:,ISTEP),RAD,SCALEX,NBX,NB,BOND)
!
!     __REMOVE UNDESIRED BONDS__________________________________________________
      IB=0
      DO I=1,NB
        IATC=BOND(I)%IAT1
        IF(.NOT.TSELECTC(IATC)) CYCLE
        IATP=BOND(I)%IAT2
        IF(.NOT.TSELECTP(IATP)) CYCLE
        IB=IB+1
        IF(IB.LT.I) BOND(IB)=BOND(I)
      ENDDO
      NB=IB
!
!     __ARRANGE BONDS SUITABLE FOR CASHING______________________________________
      DO IB=1,NB-1
        DO I=IB+1,NB
          IATC=BOND(IB)%IAT1
          IATP=BOND(IB)%IAT2
          IF(BOND(I)%IAT1.LT.IATC) THEN
            BOND1=BOND(IB)
            BOND(IB)=BOND(I)
            BOND(I)=BOND1
            IATC=BOND(IB)%IAT1  
            IATP=BOND(IB)%IAT2  
          ELSE IF(BOND(I)%IAT1.EQ.IATC) THEN
            IF(BOND(I)%IAT1.LT.IATP) THEN
              BOND1=BOND(IB)
              BOND(IB)=BOND(I)
              BOND(I)=BOND1
              IATP=BOND(IB)%IAT2
            END IF
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ACCUMULATE HISTOGRAM                                                ==
!     ==========================================================================
      SVAR=1.D0/RMAX*REAL(NP,KIND=8)
      HISTOGRAM(:)=0.D0
      DO ISTEP=1,NSTEP
        DO IB=1,NB
          IATC=BOND(IB)%IAT1
          IATP=BOND(IB)%IAT2
          DR(:)=R(:,IATP,ISTEP)-R(:,IATC,ISTEP)   
          IT=BOND(IB)%IT
          IF(SUM(IT**2).NE.0) DR(:)=DR(:)+MATMUL(RBAS,REAL(IT,KIND=8))
          DIS2=SUM((DR(:))**2)
          IF(DIS2.GE.RMAX2) CYCLE
          IP=1+INT(SQRT(DIS2)*SVAR)
          HISTOGRAM(IP)=HISTOGRAM(IP)+1.D0
        ENDDO
      ENDDO
      HISTOGRAM(:)=HISTOGRAM(:)/REAL(NSTEP,KIND=8)
!
!     __DIVIDE BY THE NUMBER OF CENTRAL ATOMS___________________________________
      HISTOGRAM(:)=HISTOGRAM(:)/REAL(NATC,KIND=8)
!
!     __CONVERT SUM INTO DENSITY________________________________________________
      HISTOGRAM(:)=HISTOGRAM(:)/(RMAX/REAL(NP,KIND=8))
!
!     ==========================================================================
!     ==  WRITE HISTOGRAM TO FILE                                             ==
!     ==========================================================================
      DO IP=1,NP
        SVAR=REAL(IP,KIND=8)/REAL(NP,KIND=8)*RMAX
        WRITE(NFIL,*)SVAR,HISTOGRAM(IP)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAIRCORRELATION2(NFIL,RBAS,NAT,TSELECTC,TSELECTP,NSTEP,R)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      INTEGER(4),INTENT(IN) :: NAT
      LOGICAL(4),INTENT(IN) :: TSELECTC(NAT)
      LOGICAL(4),INTENT(IN) :: TSELECTP(NAT)
      INTEGER(4),INTENT(IN) :: NSTEP
      REAL(8)   ,INTENT(IN) :: R(3,NAT,NSTEP)
      REAL(8)   ,PARAMETER  :: RMAX=5.D0
      REAL(8)   ,PARAMETER  :: RMAX2=RMAX**2
      INTEGER(4),PARAMETER  :: ITX=3
      INTEGER(4),PARAMETER  :: NP=1000  ! #(GRID POINTS ON HISTOGRAM)
      REAL(8)               :: HISTOGRAM(NP)
      INTEGER(4)            :: NT ! #(TRANSLATION VECTORS KEPT)
      REAL(8)   ,ALLOCATABLE:: TVEC(:,:) ! (3,NT)
      REAL(8)               :: R1(3)
      REAL(8)               :: DR(3)
      REAL(8)               :: DIS2
      REAL(8)               :: SVAR ! SUPPORT VARIABLE
      INTEGER(4)            :: IT1,IT2,IT3
      INTEGER(4)            :: IATC,IATP,IT,ISTEP,IP
!     **************************************************************************
!
!     ==========================================================================
!     == DETERMINE SET OF TRANSLATION VECTORS                                 ==
!     ==========================================================================
!     == SVAR ADDS MAIN DIAGONAL OF UNIT CELL TO MAXIMUM DISTANCE
      SVAR=SQRT(SUM((RBAS(:,1)+RBAS(:,2)+RBAS(:,3))**2))
      DO IT1=-ITX,ITX
        DO IT2=-ITX,ITX
          DO IT3=-ITX,ITX
            DR(:)=RBAS(:,1)*REAL(IT1,KIND=8) &
                 +RBAS(:,2)*REAL(IT2,KIND=8) &
                 +RBAS(:,3)*REAL(IT3,KIND=8) 
            IF(SQRT(SUM(DR**2)).LT.RMAX+SVAR)NT=NT+1
          ENDDO
        ENDDO
      ENDDO
      ALLOCATE(TVEC(3,NT))
      DO IT1=-ITX,ITX
        DO IT2=-ITX,ITX
          DO IT3=-ITX,ITX
            DR(:)=RBAS(:,1)*REAL(IT1,KIND=8) &
                 +RBAS(:,2)*REAL(IT2,KIND=8) &
                 +RBAS(:,3)*REAL(IT3,KIND=8) 
            IF(SQRT(SUM(DR**2)).LT.RMAX+SVAR) THEN
              IT=IT+1
              TVEC(:,IT)=DR(:)
            END IF
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ACCUMULATE HISTOGRAM                                                ==
!     ==========================================================================
      SVAR=1.D0/RMAX*REAL(NP,KIND=8)
      HISTOGRAM(:)=0.D0
      DO ISTEP=1,NSTEP
        DO IATC=1,NAT
          IF(.NOT.TSELECTC(IATC)) CYCLE
          R1(:)=R(:,IATC,ISTEP)
          DO IATP=1,NAT
            IF(.NOT.TSELECTP(IATP)) CYCLE
            DR(:)=R(:,IATP,ISTEP)-R1(:)
            DO IT=1,NT
              DIS2=SUM((DR(:)+TVEC(:,IT))**2)
              IF(DIS2.GE.RMAX2) CYCLE
              IP=1+INT(SQRT(DIS2)*SVAR)
              HISTOGRAM(IP)=HISTOGRAM(IP)+1.D0
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      HISTOGRAM=HISTOGRAM/REAL(NSTEP,KIND=8)
      HISTOGRAM=HISTOGRAM/REAL(NP,KIND=8)
!
!     ==========================================================================
!     ==  WRITE HISTOGRAM TO FILE                                             ==
!     ==========================================================================
      DO IP=1,NP
        SVAR=REAL(IP,KIND=8)/REAL(NP,KIND=8)*RMAX
        WRITE(NFIL,*)SVAR,HISTOGRAM(IP)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SNAPSHOT(LL_CNTL_)!
!     **************************************************************************
!     **  PICKS OUT ONE PARTICULAR TIME STEP FROM A TRAJECTORY FILE           **
!     **  AND PRODUCES PSEUDO STRUCTURE INPUT                                 **
!     **  AND A CSSR FILE                                                     **
!     **************************************************************************
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      USE TRAJECTORY_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      INTEGER(4)               :: NSNAPSHOTS
      INTEGER(4)               :: ISNAPSHOT
      LOGICAL(4)               :: TCHK,TCHK1,TCHK2
      INTEGER(4)               :: I,IAT,ISTEP,IATM
      INTEGER(4)               :: NFIL
      INTEGER(4)               :: NFILO
      CHARACTER(256)           :: FILE
      LOGICAL(4)   ,ALLOCATABLE:: TSELECT(:)
      REAL(8)                  :: PICO,SECOND
      REAL(8)                  :: TIME
      INTEGER(4)               :: NATM,NATOM
      CHARACTER(32),ALLOCATABLE:: ATOMM(:)
      REAL(8)      ,ALLOCATABLE:: POSM(:,:)
      REAL(8)      ,ALLOCATABLE:: QM(:)
!     **************************************************************************
                               CALL TRACE$PUSH('SNAPSHOT')
      IF(TQMMM) THEN
         ALLOCATE(TSELECT(QNAT))
         NATOM=QNAT
      ELSE
         ALLOCATE(TSELECT(NAT))
         NATOM=NAT
      END IF
      CALL CONSTANTS$GET('PICO',PICO)
      CALL CONSTANTS$GET('SECOND',SECOND)
      CALL FILEHANDLER$UNIT('PROT',NFILO)

      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'SNAPSHOT',NSNAPSHOTS)
      DO ISNAPSHOT=1,NSNAPSHOTS
        CALL LINKEDLIST$SELECT(LL_CNTL,'SNAPSHOT',ISNAPSHOT)
!     
!       ========================================================================
!       ==  SPECIFY FILE                                                      ==
!       ========================================================================
        CALL LINKEDLIST$EXISTL(LL_CNTL,'FILE',1,TCHK)
        IF(.NOT.TCHK) THEN
          IF(ISNAPSHOT.GT.1) THEN
            CALL ERROR$MSG('FILE MUST BE SPECIFIED FOR MULTIPLE SNAPSHOTS')
            CALL ERROR$STOP('SNAPSHOT')
          END IF
          FILE=-'.TRA.CSSR'
          TCHK=.TRUE.
          CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',0)
          CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,FILE)
          CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,TCHK)
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ELSE
          CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',1)
          CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SNAPSHOT!FILE:NAME IS NOT SPECIFIED')
            CALL ERROR$STOP('SNAPSHOT')
          END IF
          CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILE)
          TCHK=.FALSE.
          CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
          ELSE
            CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
          END IF
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        END IF
        CALL FILEHANDLER$SETFILE('SNAPSHOT',TCHK,TRIM(FILE))
        CALL FILEHANDLER$SETSPECIFICATION('SNAPSHOT','FORM','FORMATTED')
!     
!       ================================================================
!       ==  SPECIFY SNAPSHOT                                          ==
!       ================================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'T[PSEC]',1,TCHK1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'STEP',1,TCHK2)
        IF(TCHK1.AND.TCHK2) THEN
          CALL ERROR$MSG('KEYWORDS T[PSEC] AND STEP= ARE MUTUALLY EXCLUSIVE')
          CALL ERROR$STOP('SNAPSHOT')
        ENDIF
        IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'T[PSEC]',1,TIME)
          TIME=TIME*PICO*SECOND
          ISTEP=0
          DO I=1,NSTEP
            IF(TRA%T(I).LT.TIME) ISTEP=I
          ENDDO
          IF(ISTEP.EQ.0) THEN
            CALL ERROR$MSG('NO TRAJECTORY BEFORE SELECTED TIME')
            CALL ERROR$STOP('SNAPSHOT')
          ENDIF
        ELSE IF(TCHK2) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'STEP',1,ISTEP)
          IF(ISTEP.LE.0.OR.ISTEP.GT.NSTEP) THEN
            CALL ERROR$MSG('SELECTED TIMESTEP DOES NOT EXIST IN TRAJECTORY FILE')
            CALL ERROR$I4VAL('ISTEP',ISTEP)
            CALL ERROR$I4VAL('NSTEP',NSTEP)
            CALL ERROR$STOP('SNAPSHOT')
          ENDIF
        ELSE
          CALL ERROR$MSG('EITHER KEYWORDS T[PSEC] OR STEP MUST BE SELECTED')
          CALL ERROR$STOP('SNAPSHOT')
        END IF
!     
!       ==============================================================
!       ==  SELECT ATOMS                                            ==
!       ==============================================================
        CALL SELECTION(LL_CNTL,NATOM,ATOM,TSELECT)
        NATM=0
        DO IAT=1,NATOM
          IF(TSELECT(IAT)) NATM=NATM+1
        ENDDO
        ALLOCATE(ATOMM(NATM))
        ALLOCATE(POSM(3,NATM))
        ALLOCATE(QM(NATM))
        IATM=0
        DO IAT=1,NATOM
          IF(.NOT.TSELECT(IAT)) CYCLE
          IATM=IATM+1
          ATOMM(IATM)=ATOM(IAT)
          POSM(:,IATM)=TRA%R(:,IAT,ISTEP)
          QM(IATM)=TRA%Q(IAT,ISTEP)
        ENDDO
!     
!       ================================================================
!       ==  WRITE ATOMS TO PROTOCOL FILE                              ==
!       ================================================================
        CALL REPORT$TITLE(NFILO,'SNAPSHOT')
        CALL REPORT$R8VAL(NFILO,'TIME',TIME/(PICO*SECOND),'PSEC')
        WRITE(NFILO,FMT='("!LATTICE T=",3F15.8,"     ")')RBAS(:,1)
        WRITE(NFILO,FMT='("           ",3F15.8,"     ")')RBAS(:,2)
        WRITE(NFILO,FMT='("           ",3F15.8," !END")')RBAS(:,3)
        DO IATM=1,NATM
          WRITE(NFILO,FMT='("!ATOM NAME=''",A,"''",T30," R=",3F15.8," !END")') &
     &                          TRIM(ATOMM(IATM)),POSM(:,IATM)
        ENDDO
        WRITE(NFILO,FMT='("FILE ",A)')TRIM(FILE)
        WRITE(NFILO,*) 
!     
!       ================================================================
!       ==  WRITE ATOMS TO PROTOCOL FILE                              ==
!       ================================================================
        CALL FILEHANDLER$UNIT('SNAPSHOT',NFIL)
        CALL WRITECSSR(NFIL,'TEST',RBAS,NATM,ATOMM,POSM,QM,.TRUE.)
!     
!       ================================================================
!       == CLOSE DOWN                                                ==
!       ================================================================
        DEALLOCATE(ATOMM)
        DEALLOCATE(QM)
        DEALLOCATE(POSM)
        CALL FILEHANDLER$CLOSE('SNAPSHOT')
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                               CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SOFT(LL_CNTL_)
!     ******************************************************************
!     **                                                              **
!     **  PRODUCES PATH LENGTH IN MASS WEIGHTED COORDINATES           **
!     **      DS=SUM(M*DR^2)                                          **
!     **                                                              **
!     ******************************************************************
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      USE TRAJECTORY_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      INTEGER(4)               :: NPLOT
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: IPLOT,IAT,ISTEP,NATOM
      INTEGER(4)               :: NFIL
      INTEGER(4)               :: NFILO
      CHARACTER(256)           :: FILE
      REAL(8)                  :: DR(3)
      REAL(8)                  :: DIS
      LOGICAL(4)   ,ALLOCATABLE:: TSELECT(:)
      REAL(8)                  :: PICO,SECOND,ANGSTROM
      INTEGER(4)               :: N,I
      CHARACTER(16),ALLOCATABLE:: ATOMNAMES(:)
      REAL(8)                  :: S  ! PATH LENGTH 
      REAL(8)                  :: U  ! ATOMIC MASS UNIT: C12/12
!     ******************************************************************
                               CALL TRACE$PUSH('SOFT')
      IF(TQMMM) THEN
         NATOM=QNAT
         ALLOCATE(TSELECT(QNAT))
      ELSE
         NATOM=NAT
         ALLOCATE(TSELECT(NAT))
      END IF
      CALL CONSTANTS$GET('PICO',PICO)
      CALL CONSTANTS$GET('SECOND',SECOND)
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
      CALL CONSTANTS$GET('U',U)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'SOFT',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL TRACE$POP
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'SOFT')
!     
!     ================================================================
!     ==  SPECIFY FILE                                              ==
!     ================================================================
      CALL LINKEDLIST$EXISTL(LL_CNTL,'FILE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',0)
        FILE=-'.TRA.SOFT'
        CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,FILE)
        TCHK=.TRUE.
        CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,TCHK)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ELSE
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,-'.TRA.SOFT')
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILE)
!     
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
        CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      END IF
      CALL FILEHANDLER$SETFILE('SOFT',TCHK,TRIM(FILE))
      CALL FILEHANDLER$SETSPECIFICATION('SOFT','FORM','FORMATTED')
!     
!     ================================================================
!     ==  SELECT ATOMS                                              ==
!     ================================================================
      CALL SELECTION(LL_CNTL,NATOM,ATOM,TSELECT)
!     
!     ================================================================
!     ==  SOFT                                                      ==
!     ================================================================
      CALL FILEHANDLER$UNIT('SOFT',NFIL)
      REWIND(NFIL)
      S=0.D0
      WRITE(NFIL,FMT='(E20.5,10F10.5)')TRA%T(1)/(PICO*SECOND) &
     &                                ,S/(ANGSTROM*SQRT(U))
      DO ISTEP=2,NSTEP
        DIS=0.D0
        DO IAT=1,NATOM
          IF(.NOT.TSELECT(IAT)) CYCLE
          DR(:)=TRA%R(:,IAT,ISTEP)-TRA%R(:,IAT,ISTEP-1)
          DIS=DIS+SQRT(MASS(IAT)*DOT_PRODUCT(DR,DR))
        ENDDO
        S=S+DIS
        WRITE(NFIL,FMT='(E20.5,10F10.5)') &
     &       TRA%T(ISTEP)/(PICO*SECOND),S/(ANGSTROM*SQRT(U))
      ENDDO
      CALL FILEHANDLER$CLOSE('SOFT')
!     
      CALL LINKEDLIST$SELECT(LL_CNTL,'..')
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPAGHETTI(LL_CNTL_)
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      USE TRAJECTORY_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      INTEGER(4)               :: NPLOT
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: IPLOT,IAT,ISTEP,NATOM
      INTEGER(4)               :: NFIL
      INTEGER(4)               :: NFILO
      CHARACTER(256)           :: FILE
      REAL(8)                  :: DR(3)
      REAL(8)                  :: DIS
      LOGICAL(4)   ,ALLOCATABLE:: TSELECT(:)
      REAL(8)                  :: PICO,SECOND,ANGSTROM
      INTEGER(4)               :: N,I
      CHARACTER(16),ALLOCATABLE:: ATOMNAMES(:)
      CHARACTER(32)            :: FILEID
!     ******************************************************************
                               CALL TRACE$PUSH('SPAGHETTI')
      IF(TQMMM) THEN
         NATOM=QNAT
         ALLOCATE(TSELECT(QNAT))
      ELSE
         NATOM=NAT
         ALLOCATE(TSELECT(NAT))
      END IF
      CALL CONSTANTS$GET('PICO',PICO)
      CALL CONSTANTS$GET('SECOND',SECOND)
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'SPAGHETTI',NPLOT)
      DO IPLOT=1,NPLOT
        CALL LINKEDLIST$SELECT(LL_CNTL,'SPAGHETTI',IPLOT)
!
!       ================================================================
!       ==  SPECIFY FILE                                              ==
!       ================================================================
        CALL LINKEDLIST$EXISTL(LL_CNTL,'FILE',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('BLOCK !SPAGHETTI!FILE IS MANDATORY') 
          CALL ERROR$STOP('SPAGHETTI')
        END IF
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,-'.TRA.PASTA')
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILE)

        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
        CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
!
        IF(NPLOT.GT.1) THEN
          WRITE(FILEID,FMT=*)IPLOT
          FILEID='SPAGHETTI_'//ADJUSTL(FILEID)
        ELSE
          FILEID='SPAGHETTI'
        END IF
        CALL FILEHANDLER$SETFILE(FILEID,TCHK,TRIM(FILE))
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'FORM','FORMATTED')
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
!
!       ================================================================
!       ==  SELECT ATOMS                                              ==
!       ================================================================
        CALL SELECTION(LL_CNTL,NATOM,ATOM,TSELECT)
!
!       ================================================================
!       ==  SPAGHETTI                                                  ==
!       ================================================================
        CALL FILEHANDLER$UNIT(FILEID,NFIL)
        REWIND(NFIL)
        DO IAT=1,NATOM
          IF(.NOT.TSELECT(IAT)) CYCLE
          WRITE(NFIL,FMT='(E20.5,10F10.5)') &
     &         TRA%T(1)/(PICO*SECOND),0.D0
          DO ISTEP=2,NSTEP-1
            DR(:)=TRA%R(:,IAT,ISTEP)-TRA%R(:,IAT,1)
            DIS=SQRT(DOT_PRODUCT(DR,DR))
            WRITE(NFIL,FMT='(E20.5,10F10.5)') &
     &         TRA%T(ISTEP)/(PICO*SECOND),DIS/ANGSTROM
          ENDDO
          WRITE(NFIL,FMT='(E20.5,10F10.5)') &
     &         TRA%T(NSTEP)/(PICO*SECOND),0.D0
          
        ENDDO
        CALL FILEHANDLER$CLOSE(FILEID)
!
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TEMPERATURE(LL_CNTL_)
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      USE TRAJECTORY_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      INTEGER(4)               :: NPLOT
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: IPLOT,IAT,ISTEP
      INTEGER(4)               :: NFIL
      INTEGER(4)               :: NFILO
      CHARACTER(256)           :: FILE
      CHARACTER(32)            :: FILEID
      REAL(8)                  :: DR(3)
      REAL(8)                  :: DIS
      REAL(8)                  :: TAV,T2AV
      REAL(8)                  :: TSUM,T2SUM
      REAL(8)                  :: SVAR
      REAL(8)                  :: DELTAT
      LOGICAL(4)               :: TFILE   ! FILE SHALL BE WRITTEN OR NOT
      REAL(8)                  :: PICO,SECOND,ANGSTROM,KELVIN
      LOGICAL(4)   ,ALLOCATABLE:: TSELECT(:)
      INTEGER(4)               :: N,I
      CHARACTER(16),ALLOCATABLE:: ATOMNAMES(:)
      INTEGER(4)   ,PARAMETER  :: CHOICE=2
      REAL(8)                  :: SUM
      REAL(8)                  :: TRTRD,ALPHARTRD,TAVRTRD
      INTEGER(4)               :: NATOM
      CHARACTER(80)            :: STRING
!     ******************************************************************
                               CALL TRACE$PUSH('TEMPERATURE')
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'TEMPERATURE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL TRACE$POP
        RETURN
      END IF


      IF(TQMMM) THEN
        NATOM=QNAT
        ALLOCATE(TSELECT(QNAT))
      ELSE
        NATOM = NAT
        ALLOCATE(TSELECT(NAT))
      END IF
!
      CALL CONSTANTS$GET('PICO',PICO)
      CALL CONSTANTS$GET('SECOND',SECOND)
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)      
      CALL CONSTANTS$GET('KB',KELVIN)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ================================================================
!     ==  TEMPERATURE                                                  ==
!     ================================================================
      CALL REPORT$TITLE(NFILO,'TEMPERATURE OF INDIVIDUAL ATOMS')
      TSUM=0.D0
      T2SUM=0.D0
      DO IAT=1,NATOM
        TAV=0.D0
        T2AV=0.D0
        DO ISTEP=2,NSTEP-1
          DR(:)=TRA%R(:,IAT,ISTEP+1)-TRA%R(:,IAT,ISTEP-1)
          DELTAT=TRA%T(ISTEP+1)-TRA%T(ISTEP-1)
          DIS=DOT_PRODUCT(DR,DR)/DELTAT**2
          TAV =TAV +DIS
          T2AV=T2AV+DIS**2
        ENDDO
        TAV =TAV /DBLE(NSTEP-2)
        T2AV=T2AV/DBLE(NSTEP-2) 
        SVAR=0.5D0*MASS(IAT)/1.5D0
        TAV =SVAR*TAV
        T2AV=SVAR**2*T2AV
        TSUM=TSUM+TAV
        T2SUM=T2SUM+T2AV
        WRITE(NFILO,FMT='("ATOM ",A10," <T>:",F10.1,"K" &
     &                    ,"   SQR[<(T-<T>)^2>]",F15.0,"K")') &
     &    TRIM(ATOM(IAT)),TAV/KELVIN,SQRT(T2AV-TAV**2)/KELVIN
      ENDDO
      TSUM=TSUM/NATOM
      T2SUM=T2SUM/NATOM
      WRITE(NFILO,FMT='(80("-"))')
      WRITE(NFILO,FMT='("SUM",T17,"<T>:",F10.1,"K" &
     &                    ,"   SQR[<(T-<T>)^2>]",F15.0,"K")') &
     &    TSUM/KELVIN,SQRT(T2SUM-TSUM**2)/KELVIN
      STRING="NOTE THAT THE TEMPERATURE USES G=3*N INSTEAD OF" &
     &       //" G=3N-3 OR G=3N-6"
      WRITE(NFILO,FMT='(A)')TRIM(STRING)
      STRING="FOR THE NUMBER OF DEGREES OF FREEDOM."
      WRITE(NFILO,FMT='(A)')TRIM(STRING)
      STRING="THIS UNDERESTIMATES THE TRUE TEMPERATURE"
      WRITE(NFILO,FMT='(A)')TRIM(STRING)
!
!     ================================================================
!     ==  INDIVIDUAL PLOTS                                          ==
!     ================================================================
      CALL LINKEDLIST$NLISTS(LL_CNTL,'TEMPERATURE',NPLOT)
      DO IPLOT=1,NPLOT
        CALL LINKEDLIST$SELECT(LL_CNTL,'TEMPERATURE',IPLOT)
!
!       ================================================================
!       ==  SPECIFY RETARDATION                                       ==
!       ================================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'RETARD[PS]',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'RETARD[PS]',0,0.D0)
        CALL LINKEDLIST$GET(LL_CNTL,'RETARD[PS]',1,TRTRD)
        TRTRD=TRTRD*PICO*SECOND
!
!       ================================================================
!       ==  SPECIFY FILE                                              ==
!       ================================================================
! TODO: MANDATORY !FILE BUT NOT KEYS (NAME, EXT)
!       DEFAULT LEADS TO HIDDEN FILE ".TRA.TEMP"
        CALL LINKEDLIST$EXISTL(LL_CNTL,'FILE',1,TFILE)
        IF(.NOT.TFILE) THEN
          CALL ERROR$MSG('!TCNTL!TEMPERATURE!FILE NOT FOUND')
          CALL ERROR$STOP('TEMPERATURE')
        END IF
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,-'.TRA.TEMP')
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILE)
!       
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
        CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
!       
        IF(NPLOT.GT.1) THEN
          WRITE(FILEID,FMT=*)IPLOT
          FILEID='TEMPERATURE_'//ADJUSTL(FILEID)
        ELSE
          FILEID='TEMPERATURE'
        END IF
        CALL FILEHANDLER$SETFILE(FILEID,TCHK,TRIM(FILE))
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'FORM','FORMATTED')
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
!
!       ================================================================
!       ==  SELECT ATOMS                                              ==
!       ================================================================
        CALL SELECTION(LL_CNTL,NATOM,ATOM,TSELECT)
!
!       ================================================================
!       ==  TEMPERATURE                                               ==
!       ================================================================
        CALL FILEHANDLER$UNIT(FILEID,NFIL)
        TAVRTRD=0.D0
        REWIND(NFIL)
        IF(CHOICE.EQ.1) THEN
          DO IAT=1,NATOM
            IF(.NOT.TSELECT(IAT)) CYCLE
            WRITE(NFIL,*)TRA%T(1)/(PICO*SECOND),0.D0
            DO ISTEP=2,NSTEP-1
              DR(:)=TRA%R(:,IAT,ISTEP+1)-TRA%R(:,IAT,ISTEP-1)
              DELTAT=TRA%T(ISTEP+1)-TRA%T(ISTEP-1)
              DIS=DOT_PRODUCT(DR,DR)/DELTAT**2
              TAV=0.5D0*MASS(IAT)*DIS
              IF(TRTRD.LT.DELTAT.OR.ISTEP.EQ.2) THEN
                ALPHARTRD=1.D0
              ELSE
                ALPHARTRD=1.D0-EXP(-DELTAT/TRTRD)
              END IF
              TAVRTRD=ALPHARTRD*TAV+(1.D0-ALPHARTRD)*TAVRTRD
              WRITE(NFIL,*)TRA%T(ISTEP)/(PICO*SECOND),TAVRTRD/KELVIN
            ENDDO
            WRITE(NFIL,*)TRA%T(NSTEP)/(PICO*SECOND),0.D0
          ENDDO
        ELSE IF(CHOICE.EQ.2) THEN
          TAVRTRD=0.D0
          DO ISTEP=2,NSTEP-1
            TAV=0.D0
            SUM=0.D0
            DO IAT=1,NAT
              IF(.NOT.TSELECT(IAT)) CYCLE
              DR(:)=TRA%R(:,IAT,ISTEP+1)-TRA%R(:,IAT,ISTEP-1)
              DELTAT=TRA%T(ISTEP+1)-TRA%T(ISTEP-1)
              DIS=DOT_PRODUCT(DR,DR)/DELTAT**2
              TAV=TAV+0.5D0*MASS(IAT)*DIS
              SUM=SUM+1.D0
            ENDDO
            TAV=TAV/(1.5D0*SUM)
            IF(TRTRD.LT.DELTAT.OR.ISTEP.EQ.2) THEN
              ALPHARTRD=1.D0
            ELSE
              ALPHARTRD=1.D0-EXP(-DELTAT/TRTRD)
            END IF
            TAVRTRD=ALPHARTRD*TAV+(1.D0-ALPHARTRD)*TAVRTRD
            WRITE(NFIL,*)TRA%T(ISTEP)/(PICO*SECOND),TAVRTRD/KELVIN
          ENDDO
        ELSE 
          CALL ERROR$MSG('CHOICE NOT RECOGNIZED')
          CALL ERROR$STOP('TEMPERATURE')
        END IF
        CALL FILEHANDLER$CLOSE(FILEID)
!
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CELL(LL_CNTL_)
!     **************************************************************************
!     **  EXPORT UNIT CELL MATRIX AS FUNCTION OF TIME TO A FILE               **
!     **************************************************************************
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      USE TRAJECTORY_MODULE, ONLY: TQMMM &
     &                            ,NSTEP &
     &                            ,TRA
      IMPLICIT NONE
      TYPE(LL_TYPE), INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)             :: LL_CNTL
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: NFILO
      INTEGER(4)                :: NFIL
      CHARACTER(256)            :: FILE
      INTEGER(4)                :: I  
      REAL(8)                   :: ANGSTROM
      REAL(8)                   :: PICO
      REAL(8)                   :: SECOND
!     **************************************************************************
      CALL TRACE$PUSH('CELL')
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
      CALL CONSTANTS$GET('PICO',PICO)
      CALL CONSTANTS$GET('SECOND',SECOND)
      CALL FILEHANDLER$UNIT('PROT',NFILO)

      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'CELL',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL TRACE$POP
        RETURN
      END IF
      CALL LINKEDLIST$SELECT(LL_CNTL,'CELL')
!     
!     ==========================================================================
!     ==  SPECIFY FILE                                                        ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTL(LL_CNTL,'FILE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',0)
        FILE=-'.TRA.CELL'
        CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,TRIM(FILE))
        TCHK=.TRUE.
        CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,TCHK)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ELSE
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,-'.TRA.CELL')
          CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.TRUE.)
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILE)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(.NOT.TCHK) CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
        CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      END IF
      CALL FILEHANDLER$SETFILE('CELL',TCHK,TRIM(FILE))
      CALL FILEHANDLER$SETSPECIFICATION('CELL','FORM','FORMATTED')
!     
!     ==========================================================================
!     ==  CHECK FOR TIME FORMAT FLAG                                          ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'TIMEFORMAT',1,TCHK)
      IF(.NOT.TCHK) CALL LINKEDLIST$SET(LL_CNTL,'TIMEFORMAT',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'TIMEFORMAT',1,TCHK)
!     
!     ==========================================================================
!     ==  REPORT TO PROTOCOL FILE                                             ==
!     ==========================================================================
      CALL REPORT$TITLE(NFILO,'CELL')
      CALL REPORT$CHVAL(NFILO,'CELL FILE',FILE)
      CALL REPORT$L4VAL(NFILO,'TIMEFORMAT',TCHK)
!     
!     ==========================================================================
!     ==  WRITE DATA TO FILE                                                  ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('CELL',NFIL)
!       OUTPUT EITHER SIMULATION TIME IN PS OR ITERATION STEPS
      IF(TCHK) THEN
        WRITE(NFIL,FMT='("#",A9,9A10)')'TIME[PS]','T1X','T1Y','T1Z' &
     &                                           ,'T2X','T2Y','T2Z' &
     &                                           ,'T3X','T3Y','T3Z'
        DO I=1,NSTEP
          WRITE(NFIL,FMT='(F10.5,9F10.5)')TRA%T(I)/(PICO*SECOND) &
     &                                   ,TRA%CELL(:,I)/ANGSTROM
        ENDDO
      ELSE
        WRITE(NFIL,FMT='("#",A9,9A10)')'ITER','T1X','T1Y','T1Z' &
     &                                       ,'T2X','T2Y','T2Z' &
     &                                       ,'T3X','T3Y','T3Z'
        DO I=1,NSTEP
          WRITE(NFIL,FMT='(I10,9F10.5)')TRA%ISTEP(I),TRA%CELL(:,I)/ANGSTROM
        ENDDO
      END IF
!
      CALL FILEHANDLER$CLOSE('CELL')
      CALL TRACE$POP
      RETURN
      END SUBROUTINE CELL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEMODES(LL_CNTL_)
      USE LINKEDLIST_MODULE
      USE TRAJECTORY_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      INTEGER(4)               :: NPR,I
      CHARACTER(256)           :: NAME
      INTEGER(4)               :: IMODE,ISTEP
      INTEGER(4)               :: NFIL
      INTEGER(4)               :: NFILO
      REAL(8)      ,ALLOCATABLE:: ARRAY(:)
      REAL(8)                  :: PSEC
      REAL(8)                  :: SVAR
      LOGICAL(4)               :: TCHK
!     ******************************************************************
                                CALL TRACE$PUSH('WRITEMODES')
      LL_CNTL=LL_CNTL_
!
      CALL CONSTANTS$GET('SECOND',PSEC)
      CALL CONSTANTS$GET('PICO',SVAR)
      PSEC=PSEC*SVAR
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,FMT='("TRAJECTORY CONTAINS ",I10," STEPS")')NSTEP
!
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'OUTPUT',NPR)
      ALLOCATE(ARRAY(NSTEP))
      DO I=1,NPR      
        CALL LINKEDLIST$SELECT(LL_CNTL,'OUTPUT',I)
!
!       ================================================================
!       ==                                                            ==
!       ================================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'FILENAME',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!TCNTL!OUTPUT:FILENAME NOT SPECIFIED')
          CALL ERROR$STOP('WRITEMODES')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'FILENAME',1,NAME)
        CALL FILEHANDLER$SETFILE('MODEOUT',.FALSE.,NAME)
        CALL FILEHANDLER$SETSPECIFICATION('MODEOUT','STATUS','UNKNOWN')
        CALL FILEHANDLER$SETSPECIFICATION('MODEOUT','POSITION','ASIS')
        CALL FILEHANDLER$SETSPECIFICATION('MODEOUT','ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION('MODEOUT','FORM','FORMATTED')
        CALL FILEHANDLER$UNIT('MODEOUT',NFIL)
!
!       ================================================================
!       ==  GET MODE ID                                               ==
!       ================================================================
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,NAME)
        CALL TRAJECTORY_MODELOOKUP(NAME,IMODE)
!
!       ================================================================
!       ==  TRANSFORMATIONS                                           ==
!       ================================================================
        CALL LINKEDLIST$EXISTD(LL_CNTL,'TYPE',1,TCHK)
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'TYPE',1,NAME)
        IF(NAME.EQ.'VELOCITY') THEN
          CALL VELOCITY(NSTEP,TRA%T,MODE(IMODE)%VAL(:),ARRAY)
        ELSE
          ARRAY(:)=MODE(IMODE)%VAL(:)
        END IF
!
!       ================================================================
!       ==  WRITE TO FILE                                             ==
!       ================================================================
        DO ISTEP=1,NSTEP
          WRITE(NFIL,*)TRA%T(ISTEP)/PSEC,ARRAY(ISTEP)
        ENDDO
        CALL FILEHANDLER$CLOSE('MODEOUT')
                         
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
      DEALLOCATE(ARRAY)
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MODES(LL_CNTL_)
      USE LINKEDLIST_MODULE
      USE TRAJECTORY_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      INTEGER(4)               :: IMODE
      INTEGER(4)               :: I,N,J
      INTEGER(4)               :: IAT1,IAT2,IAT3,IAT4
      REAL(8)                  :: SCALE
      CHARACTER(32)            :: NAME
      REAL(8)      ,ALLOCATABLE:: ARRAY(:)
!     ******************************************************************
                                CALL TRACE$PUSH('MODES')
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'MODE',NMODES)
      CALL TRAJECTORY_NEWMODES
      ALLOCATE(ARRAY(NSTEP))
      DO IMODE=1,NMODES
        CALL LINKEDLIST$SELECT(LL_CNTL,'MODE',IMODE)
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,MODE(IMODE)%ID)
PRINT*,'PREPARING MODE:',MODE(IMODE)%ID
        MODE(IMODE)%VAL(:)=0.D0
        CALL LINKEDLIST$NLISTS(LL_CNTL,'BOND',N)
        DO I=1,N
          CALL LINKEDLIST$SELECT(LL_CNTL,'BOND',I)
          CALL LINKEDLIST$GET(LL_CNTL,'SCALE',1,SCALE)
          CALL LINKEDLIST$GET(LL_CNTL,'ATOM1',1,NAME)
          CALL TRAJECTORY_ATOMLOOKUP(NAME,IAT1)
PRINT*,'BOND: ATOM1=',NAME,IAT1
          CALL LINKEDLIST$GET(LL_CNTL,'ATOM2',1,NAME)
          CALL TRAJECTORY_ATOMLOOKUP(NAME,IAT2)
PRINT*,'BOND: ATOM1=',NAME,IAT2
          CALL BOND(IAT1,IAT2,ARRAY)

          MODE(IMODE)%VAL(:)=MODE(IMODE)%VAL(:)+ARRAY(:)*SCALE
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
!
!       ================================================================
!       ==  ANGLES                                                    ==
!       ================================================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'ANGLE',N)
        DO I=1,N
          CALL LINKEDLIST$SELECT(LL_CNTL,'ANGLE',I)
          CALL LINKEDLIST$GET(LL_CNTL,'SCALE',1,SCALE)
          CALL LINKEDLIST$GET(LL_CNTL,'ATOM1',1,NAME)
          CALL TRAJECTORY_ATOMLOOKUP(NAME,IAT1)
          CALL LINKEDLIST$GET(LL_CNTL,'ATOM2',1,NAME)
          CALL TRAJECTORY_ATOMLOOKUP(NAME,IAT2)
          CALL LINKEDLIST$GET(LL_CNTL,'ATOM3',1,NAME)
          CALL TRAJECTORY_ATOMLOOKUP(NAME,IAT3)
          CALL ANGLE(IAT1,IAT2,IAT3,ARRAY)
          MODE(IMODE)%VAL(:)=MODE(IMODE)%VAL(:)+ARRAY(:)*SCALE
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
!
!       ================================================================
!       ==  TORSIONS                                                  ==
!       ================================================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'TORSION',N)
        DO I=1,N
          CALL LINKEDLIST$SELECT(LL_CNTL,'TORSION',I)
          CALL LINKEDLIST$GET(LL_CNTL,'SCALE',1,SCALE)
          CALL LINKEDLIST$GET(LL_CNTL,'ATOM1',1,NAME)
          CALL TRAJECTORY_ATOMLOOKUP(NAME,IAT1)
          CALL LINKEDLIST$GET(LL_CNTL,'ATOM2',1,NAME)
          CALL TRAJECTORY_ATOMLOOKUP(NAME,IAT2)
          CALL LINKEDLIST$GET(LL_CNTL,'ATOM3',1,NAME)
          CALL TRAJECTORY_ATOMLOOKUP(NAME,IAT3)
          CALL LINKEDLIST$GET(LL_CNTL,'ATOM4',1,NAME)
          CALL TRAJECTORY_ATOMLOOKUP(NAME,IAT4)
          CALL TORSION(IAT1,IAT2,IAT3,IAT4,ARRAY)
          MODE(IMODE)%VAL(:)=MODE(IMODE)%VAL(:)+ARRAY(:)*SCALE
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
!
!       ================================================================
!       ==  PREVIOUSLY DEFINED MODE                                   ==
!       ================================================================
        CALL LINKEDLIST$NLISTS(LL_CNTL,'MODE',N)
        DO I=1,N
          CALL LINKEDLIST$SELECT(LL_CNTL,'MODE',I)
          CALL LINKEDLIST$GET(LL_CNTL,'SCALE',1,SCALE)
          CALL LINKEDLIST$GET(LL_CNTL,'ID',1,NAME)
          CALL TRAJECTORY_MODELOOKUP(NAME,J)
          MODE(IMODE)%VAL(:)=MODE(IMODE)%VAL(:)+MODE(J)%VAL(:)*SCALE
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
!
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
      DEALLOCATE(ARRAY)
                                CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READTRA(NFIL,T1,T2,DT)
!     **************************************************************************
!     **  READ TRAJECTORY FILE                                                **
!     **************************************************************************
      USE TRAJECTORY_MODULE, ONLY : TQMMM &
     &                             ,RBAS &
     &                             ,NAT &
     &                             ,QNAT &
     &                             ,ATOM &
     &                             ,NSTEP &
     &                             ,TRA 
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: T1
      REAL(8)   ,INTENT(IN)  :: T2
      REAL(8)   ,INTENT(IN)  :: DT
      REAL(8)                :: TLAST
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4)             :: LENG
      INTEGER(4)             :: ISTART,IREC
      INTEGER(4)             :: I,J
      INTEGER(4)             :: NFILO
      REAL(8)                :: TIME
      REAL(8)                :: PICO
      REAL(8)                :: SECOND
      REAL(8)                :: FIRSTONFILE,LASTONFILE
      LOGICAL(4)             :: TCHK,TSPINREAD
!     **************************************************************************
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL CONSTANTS$GET('SECOND',SECOND)
      CALL CONSTANTS$GET('PICO',PICO)
!
!     ==========================================================================
!     ==  CALCULATE #(ATOMS)                                                  ==
!     ==========================================================================
      REWIND(NFIL)
      READ(NFIL)I,TIME,LENG 
      FIRSTONFILE=TIME
!
      IF(.NOT.TQMMM) THEN
        IF(LENG.EQ.9+3*NAT+NAT+4*NAT) THEN
          TSPINREAD=.TRUE.              ! TRAJECTORY USES SPINS
        ELSE IF(LENG.EQ.9+3*NAT+NAT) THEN
          TSPINREAD=.FALSE.             ! NO SPIN INFORMATION IN TRAJECTORY 
        ELSE
          CALL ERROR$MSG('#(ATOMS) ON TRAJECTORY FILE INCONSISTENT')
          CALL ERROR$MSG('                          WITH STRC FILE')
          CALL ERROR$I4VAL('ARRAY SIZE',LENG)
          CALL ERROR$I4VAL('#(ATOMS) ON STRC FILE',NAT)
          CALL ERROR$STOP('READTRA')
        END IF           
      ELSE
        TSPINREAD=.FALSE.
        IF(LENG.NE.9+4*QNAT) THEN
          CALL ERROR$MSG('#(ATOMS) ON TRAJECTORY FILE INCONSISTENT')
          CALL ERROR$MSG('                          WITH STRC FILE')
          CALL ERROR$I4VAL('#(ATOMS) ON TRA FILE',LENG/4)
          CALL ERROR$I4VAL('#(ATOMS) ON STRC FILE',QNAT)
          CALL ERROR$STOP('READTRA')
        END IF
      END IF
!
!     ==========================================================================
!     ==  CALCULATE NUMBER OF RECORDS                                         ==
!     ==========================================================================
      REWIND(NFIL)
      NSTEP=0
      TLAST=-1.D+6
      IREC=0
      ISTART=0
      DO  
        IREC=IREC+1
        READ(NFIL,END=100)J,TIME
        LASTONFILE=TIME
        IF(TIME.LT.T1) THEN
          ISTART=IREC
          NSTEP=0
          TLAST=-1.D+6
          CYCLE
        END IF
        IF(TIME-TLAST.LT.DT) CYCLE
        IF(TIME.LE.T2) THEN
          TLAST=TIME
          NSTEP=NSTEP+1
        END IF
      ENDDO
 100  CONTINUE
      TRA%NSTEP=NSTEP
      IF(NSTEP.EQ.0) THEN
        CALL ERROR$MSG('TRAJECTORY IN THE SELECTED TIME WINDOW')
        CALL ERROR$R8VAL('TIME[PSEC] OF THE FIRST STEP ON FILE' &
     &                  ,FIRSTONFILE/(PICO*SECOND))
        CALL ERROR$R8VAL('TIME[PSEC] OF THE LAST TIME ON FILE' &
     &                  ,LASTONFILE/(PICO*SECOND))
        CALL ERROR$STOP('READTRA')
      END IF
!
!     ==========================================================================
!     ==  ALLOCATE TRAJECTORY AND READ                                        ==
!     ==========================================================================
!CLEMENS
      ALLOCATE(TRA%CELL(9,NSTEP))
      IF(.NOT.TQMMM) THEN
         ALLOCATE(TRA%R(3,NAT,NSTEP))
         ALLOCATE(TRA%Q(NAT,NSTEP))
      ELSE
         ALLOCATE(TRA%R(3,QNAT,NSTEP))
         ALLOCATE(TRA%Q(QNAT,NSTEP))
      END IF
      IF(TSPINREAD) THEN
         ALLOCATE(TRA%SPIN(3,NAT,NSTEP))
         ALLOCATE(TRA%CHARGE(NAT,NSTEP)) 
      END IF
      ALLOCATE(TRA%T(NSTEP))
      ALLOCATE(TRA%ISTEP(NSTEP))
!
!     == POSITION FILE AT BEGINNING OF SELECTED TIME-INTERVAL ==================
      REWIND(NFIL)
      DO IREC=1,ISTART
        READ(NFIL)
      ENDDO
!
!     == READ TRAJECTORY =======================================================
      TLAST=-1.D+9
      IREC=0
      I=1
      DO 
!CLEMENS
        IF(TSPINREAD) THEN
          READ(NFIL,END=200)TRA%ISTEP(I),TRA%T(I),LENG,TRA%CELL(:,I) &
    &                      ,TRA%R(:,:,I),TRA%Q(:,I) &
    &                      ,TRA%CHARGE(:,I),TRA%SPIN(:,:,I)
        ELSE
          READ(NFIL,END=200)TRA%ISTEP(I),TRA%T(I),LENG,TRA%CELL(:,I) &
    &                      ,TRA%R(:,:,I),TRA%Q(:,I)
        END IF
        IF(TRA%T(I)-TLAST.LT.DT) CYCLE
        TLAST=TRA%T(I)
        I=I+1
        IF(I.GT.NSTEP) EXIT
      ENDDO
 200  CONTINUE
!
!     ==========================================================================
!     ==  REPORT                                                              ==
!     ==========================================================================
      CALL REPORT$TITLE(NFILO,'TRAJECTORY INFORMATION')
      CALL REPORT$R8VAL(NFILO,'TRAJECTORY STARTS AT' &
     &                 ,FIRSTONFILE/(PICO*SECOND),'PSEC')
      CALL REPORT$R8VAL(NFILO,'TRAJECTORY ENDS AT' &
     &                       ,LASTONFILE/(PICO*SECOND),'PSEC')
      CALL REPORT$R8VAL(NFILO,'SEQUENCE READ STARTS AT' &
     &                        ,TRA%T(1)/(PICO*SECOND),'PSEC')
      CALL REPORT$R8VAL(NFILO,'SEQUENCE READ ENDS AT' &
     &                        ,TRA%T(NSTEP)/(PICO*SECOND),'PSEC') 
      CALL REPORT$I4VAL(NFILO,'SEQUENCE CONTAINS',NSTEP,'SNAPSHOTS') 
      WRITE(NFILO,*)
      WRITE(NFILO,FMT='("STRUCTURE ON THE FIRST FRAME:"/30("-"))')
      WRITE(NFILO,FMT='("T1 ",T10,3F10.5)')RBAS(:,1)
      WRITE(NFILO,FMT='("T2 ",T10,3F10.5)')RBAS(:,2)
      WRITE(NFILO,FMT='("T3 ",T10,3F10.5)')RBAS(:,3)
      WRITE(NFILO,FMT='(A10,T23,A5)')'ATOM','R(T1)'
      DO I=1,NAT
        WRITE(NFILO,FMT='(A,T10,3F10.5)')ATOM(I),TRA%R(:,I,1)
      ENDDO
      WRITE(NFILO,*)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITETRA(LL_STRC_,LL_CNTL_)
!     **************************************************************************
!     **  WRITE MOVIE FILE                                                    **
!     **************************************************************************
      USE TRAJECTORY_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE 
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN)  :: LL_STRC_
      TYPE(LL_TYPE),INTENT(IN)  :: LL_CNTL_
      INTEGER(4)                :: NFIL
      TYPE(LL_TYPE)             :: LL_STRC
      TYPE(LL_TYPE)             :: LL_CNTL
      LOGICAL(4)                :: TCHK
      LOGICAL(4)   ,PARAMETER   :: TTEST=.FALSE.
      LOGICAL(4)                :: TEXT
      INTEGER(4)                :: LENG
      REAL(8)                   :: SVAR
      REAL(8)                   :: BOXR0(3)
      REAL(8)                   :: BOXVEC(3,3)
      INTEGER(4)                :: ISTEP,IAT
      INTEGER(4)                :: NATM
      REAL(8)      ,ALLOCATABLE :: POSM(:,:)
      INTEGER(4)   ,ALLOCATABLE :: MAP(:)
      REAL(8)      ,ALLOCATABLE :: RAD(:)
      REAL(8)      ,ALLOCATABLE :: COLOR(:,:)
      INTEGER(4)                :: NBOND
      INTEGER(4)   ,ALLOCATABLE :: BOND(:,:)
      INTEGER(4)                :: IVEC(3)
      INTEGER(4)                :: IOBJ
      INTEGER(4)                :: ISKIP
      INTEGER(4)                :: IFRAME
      INTEGER(4)                :: NFILO
      CHARACTER(32)             :: FORMAT
      CHARACTER(2) ,ALLOCATABLE :: EL(:)
      CHARACTER(265)            :: FILE
      LOGICAL(4)   ,ALLOCATABLE :: TSELECT(:)
      INTEGER(4)                :: NAT0    ! #(SELECTED ATOMS)
      REAL(8)      ,ALLOCATABLE :: R0(:,:) ! POSITIONS OF SELECTED ATOMS
      INTEGER(4)   ,ALLOCATABLE :: IZ0(:)  ! ATOMIC NUMBER OF SELECTED ATOMS
      CHARACTER(32),ALLOCATABLE :: ATOM0(:)! NAMES OF SELECTED ATOMS
      REAL(8)      ,ALLOCATABLE :: Q0(:)   ! CHARGES OF SELECTED ATOMS
      INTEGER(4)                :: IZTHIS,IAT1,I
      CHARACTER(32),ALLOCATABLE :: ATOMM(:) ! 
      REAL(8)      ,ALLOCATABLE :: QM(:) 
      LOGICAL(4)                :: TBOX,TQMMM_OUT
      REAL(8)      ,ALLOCATABLE :: SCALEDRAD(:)
!     **************************************************************************
                           CALL TRACE$PUSH('WRITETRA')
      LL_STRC=LL_STRC_
      LL_CNTL=LL_CNTL_
      TQMMM_OUT=.FALSE.

      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')

      CALL LINKEDLIST$EXISTD(LL_CNTL,'QMMM',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'QMMM',0,TQMMM_OUT)
      CALL LINKEDLIST$GET(LL_CNTL,'QMMM',1,TQMMM_OUT)

      CALL LINKEDLIST$SELECT(LL_CNTL,'MOVIE')
      IF(TQMMM_OUT) THEN
         ALLOCATE(TSELECT(QNAT))
      ELSE
         ALLOCATE(TSELECT(NAT))
      END IF

!
!     ==================================================================
!     ==  SELECT FORMAT                                               ==
!     ==================================================================
      FORMAT='DX'
      CALL LINKEDLIST$EXISTD(LL_CNTL,'FORMAT',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'FORMAT',1,FORMAT)
      FORMAT=+FORMAT
      IF(FORMAT.NE.'DX'.AND.FORMAT.NE.'XYZ') THEN 
        CALL ERROR$MSG('MOVIE FILE FORMAT NOT RECOGNIZED')
        CALL ERROR$STOP('WRITETRA')
      END IF
!
!     ================================================================
!     ==  SPECIFY FILE                                              ==
!     ================================================================
      CALL LINKEDLIST$EXISTL(LL_CNTL,'FILE',1,TCHK)
      CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',0)
!     
!     == READ FLAG FOR EXTENSION =======================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
      IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
      CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TEXT)
!
!     == READ FILE NAME ================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NAME',1,TCHK)
      IF(.NOT.TCHK)THEN
        TEXT=.TRUE. 
        IF(FORMAT.EQ.'DX') THEN 
          CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,-'.TRA.MOVIE.DX')
        ELSE IF(FORMAT.EQ.'XYZ') THEN 
          CALL LINKEDLIST$SET(LL_CNTL,'NAME',0,-'.TRA.MOVIE.XYZ')
        ELSE
          CALL ERROR$MSG('MOVIE FILE FORMAT NOT RECOGNIZED')
          CALL ERROR$STOP('WRITETRA')
        END IF
      END IF
      CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILE)
!     
      CALL FILEHANDLER$SETFILE('MOVIE',TEXT,FILE)
      CALL FILEHANDLER$SETSPECIFICATION('MOVIE','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('MOVIE','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('MOVIE','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('MOVIE','FORM','FORMATTED')
!
      CALL LINKEDLIST$SELECT(LL_CNTL,'..')
!
!     ==================================================================
!     ==  GET VIEWBOX                                                 ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'VIEWBOX')
      TBOX=.FALSE.
      BOXR0(:)=0.D0
      BOXVEC(:,:)=RBAS(:,:)
      CALL LINKEDLIST$EXISTD(LL_CNTL,'T',1,TCHK)
      IF(TCHK) THEN
        TBOX=.TRUE.
        CALL LINKEDLIST$GET(LL_CNTL,'T',1,BOXVEC)
      END IF
      CALL LINKEDLIST$EXISTD(LL_CNTL,'C',1,TCHK)  ! CENTER OF THE BOX
      IF(TCHK) THEN
        TBOX=.TRUE.
        CALL LINKEDLIST$GET(LL_CNTL,'C',1,BOXR0)
        BOXR0(:)=BOXR0(:)-0.5D0*(BOXVEC(:,1)+BOXVEC(:,2)+BOXVEC(:,3))
        CALL LINKEDLIST$EXISTD(LL_CNTL,'O',1,TCHK)
        IF(TCHK) THEN
          CALL ERROR$MSG('OPTIONS O= AND C= MUST NOT BE SELECTED SIMULANEOUSLY')
          CALL ERROR$STOP('WRITETRA')
        END IF
      END IF
      CALL LINKEDLIST$EXISTD(LL_CNTL,'O',1,TCHK)
      IF(TCHK) THEN
        TBOX=.TRUE.
        CALL LINKEDLIST$GET(LL_CNTL,'O',1,BOXR0)
      END IF
!
!     ==================================================================
!     ==  DEFINE FRAMES TO BE SKIPPED                                 ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$SELECT(LL_CNTL,'MOVIE')
      ISKIP=0
      CALL LINKEDLIST$EXISTD(LL_CNTL,'SKIP',1,TCHK)
      IF(TCHK) CALL LINKEDLIST$GET(LL_CNTL,'SKIP',1,ISKIP)
!
!     ==================================================================
!     ==  SELECT ATOMS                                                ==
!     ==================================================================
      IF(.NOT.TQMMM_OUT)  CALL SELECTION(LL_CNTL,NAT,ATOM,TSELECT)
      IF(TQMMM_OUT)  CALL SELECTION(LL_CNTL,QNAT,ATOM,TSELECT)
      NAT0=0
      IF(.NOT.TQMMM_OUT) THEN
         DO IAT=1,NAT
            IF(TSELECT(IAT)) NAT0=NAT0+1
         ENDDO
      ELSE
         DO IAT=1,QNAT
            IF(TSELECT(IAT)) NAT0=NAT0+1
         ENDDO
      END IF

      ALLOCATE(R0(3,NAT0))
      ALLOCATE(IZ0(NAT0))
      ALLOCATE(ATOM0(NAT0))
      ALLOCATE(Q0(NAT0))
      IAT1=0
      IF(.NOT.TQMMM_OUT) THEN
         DO IAT=1,NAT
            IF(TSELECT(IAT)) THEN
               IAT1=IAT1+1
               IZ0(IAT1)=IZ(IAT)
               ATOM0(IAT1)=ATOM(IAT)
            END IF
         ENDDO
      ELSE
         DO IAT=1,QNAT
            IF(TSELECT(IAT)) THEN
               IAT1=IAT1+1
               IZ0(IAT1)=IZ(IAT)
               ATOM0(IAT1)=ATOM(IAT)
            END IF
         ENDDO
      END IF
!
!     ==================================================================
!     ==  REPORT INPUT                                                ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL REPORT$TITLE(NFILO,'MOVIE')
      WRITE(NFILO,FMT='("MOVIE FILE               : ",A)')TRIM(FILE)
      WRITE(NFILO,FMT='("CENTER OF VIEWBOX        :",3F10.1)') &
     &             BOXR0(:)+0.5D0*(BOXVEC(:,1)+BOXVEC(:,2)+BOXVEC(:,3))
      WRITE(NFILO,FMT='("CORNER OF VIEWBOX        :",3F10.1)')BOXR0
      WRITE(NFILO,FMT='("VECTOR SPANNING VIEWBOX  :",3F10.1)')BOXVEC(:,1)
      WRITE(NFILO,FMT='("VECTOR SPANNING VIEWBOX  :",3F10.1)')BOXVEC(:,2)
      WRITE(NFILO,FMT='("VECTOR SPANNING VIEWBOX  :",3F10.1)')BOXVEC(:,3)
      WRITE(NFILO,FMT='("ATOMS DISPLAYED:")')
      WRITE(NFILO,FMT='(10A8)')ATOM0
!     
!     ==================================================================
!     ==  BIG LOOP OVER TRAJECTORY                                    ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('MOVIE',NFIL)
      IFRAME=0
      DO ISTEP=1,NSTEP,ISKIP+1
        IFRAME=IFRAME+1
!     
!       ================================================================
!       ==  SELECT ATOMS AND MAP POSITION VECTOR INTO R0              ==
!       ================================================================
        IAT1=0
        IF(.NOT.TQMMM_OUT) THEN
           DO IAT=1,NAT
              IF(TSELECT(IAT)) THEN
                 IAT1=IAT1+1
                 R0(:,IAT1)=TRA%R(:,IAT,ISTEP)
                 Q0(IAT1)  =TRA%Q(IAT,ISTEP)
              END IF
           ENDDO
        ELSE
           DO IAT=1,QNAT
              IF(TSELECT(IAT)) THEN
                 IAT1=IAT1+1
                 R0(:,IAT1)=TRA%R(:,IAT,ISTEP)
                 Q0(IAT1)  =TRA%Q(IAT,ISTEP)
              END IF
           ENDDO
        END IF
           
!     
!       ================================================================
!       ==  MAP ATOMS INTO VIEWBOX                                    ==
!       ================================================================
!TBOX=.FALSE.
        IF(.NOT.TBOX) THEN  ! FIX DUE TO COMPILER BUG OF ABSOFT
          NATM=NAT0
        ELSE
          CALL MODEL$NATM(BOXR0,BOXVEC,NAT0,R0,RBAS,NATM)
        END IF
        ALLOCATE(POSM(3,NATM))
        ALLOCATE(MAP(NATM))
        IF(.NOT.TBOX) THEN  ! FIX DUE TO COMPILER BUG OF ABSOFT
          DO I=1,NAT0
            MAP(I)=I
            POSM(:,I)=R0(:,I)
          ENDDO
        ELSE
          CALL MODEL$ATOMS(BOXR0,BOXVEC,NAT0,R0,RBAS,NATM,POSM,MAP)
        END IF
        ALLOCATE(RAD(NATM))
        ALLOCATE(SCALEDRAD(NATM))
        ALLOCATE(COLOR(3,NATM))
        ALLOCATE(EL(NATM))
        ALLOCATE(ATOMM(NATM))
        ALLOCATE(QM(NATM))
        DO IAT=1,NATM
          IZTHIS=IZ0(MAP(IAT))
!PRINT*,"FLAG: ", IZTHIS, EL(IAT), IZ0(MAP(IAT)), MAP(IAT)
          CALL PERIODICTABLE$GET(IZTHIS,'R(COV)',RAD(IAT))
          CALL PERIODICTABLE$GET(IZTHIS,'SYMBOL',EL(IAT))
          CALL ATOMCOLOR(IZTHIS,IVEC)
          ATOMM(IAT)=ATOM0(MAP(IAT))
          COLOR(:,IAT)=REAL(IVEC,KIND=8)/200.D0
        ENDDO
!     
!       ================================================================
!       ==  CALCULATE BONDS                                           ==
!       ================================================================
        SCALEDRAD(:)=1.2D0*RAD(:)
        CALL MODEL$NBONDM(NATM,POSM,SCALEDRAD(:),NBOND)
        ALLOCATE(BOND(2,NBOND))
        CALL MODEL$BONDS(NATM,POSM,SCALEDRAD(:),NBOND,BOND)
!     
!       ================================================================
!       ==  WRITE CSSR FILE FOR TEST                                  ==
!       ================================================================
        IF(TTEST) THEN
          CALL FILEHANDLER$UNIT('PROT',NFIL)
          CALL WRITECSSR(NFIL,'TEST',RBAS,NATM,ATOMM,POSM,QM,.TRUE.)
        END IF
!     
!       ================================================================
!       ==  WRITE BALLSTICK MODEL TO DATAEXPLORER FILE                ==
!       ================================================================
        IF(FORMAT.EQ.'DX') THEN
          IOBJ=6*(IFRAME-1)
          SCALEDRAD(:)=0.5D0*RAD(:)
          CALL DXBALLSTICK(NFIL,IOBJ,NATM,COLOR,SCALEDRAD(:),POSM,NBOND,BOND &
      &         ,BOXR0,BOXVEC)
        ELSE IF(FORMAT.EQ.'XYZ') THEN
!PRINT*,"FLAG CALL WRITEXYZ(NATM): ",NATM
!          CALL WRITEXYZ(NFIL,IFRAME,NATM,EL,POSM)
          CALL WRITEEXTXYZ(NFIL,IFRAME,NATM,EL,POSM &
      &                   ,TRA%CELL(:,ISTEP),TRA%T(ISTEP))
        ELSE
          CALL ERROR$MSG('FORMAT NOT RECOGNIZED')
          CALL ERROR$STOP('WRITETRA')
        END IF
!
!       ================================================================
!       ==  DEALLOCATE                                                ==
!       ================================================================
        DEALLOCATE(EL)
        DEALLOCATE(POSM)
        DEALLOCATE(MAP)
        DEALLOCATE(RAD)
        DEALLOCATE(SCALEDRAD)
        DEALLOCATE(COLOR)
        DEALLOCATE(BOND)
        DEALLOCATE(ATOMM)
        DEALLOCATE(QM)
      ENDDO
!
!     ================================================================
!     ==  WRITE ENDING OF DATAEXPLORER FILE                         ==
!     ================================================================
      IF(FORMAT.EQ.'DX') THEN
        WRITE(NFIL,FMT='(A)')-"OBJECT ""SERIES"""
        WRITE(NFIL,FMT='(A)')-"CLASS SERIES"
        IFRAME=0
        DO ISTEP=1,NSTEP,ISKIP+1
          IFRAME=IFRAME+1
          IOBJ=6*(IFRAME-1)+6
          WRITE(NFIL,FMT='(A,1X,I10,1X,A,1X,I10,1X,A,1X,F20.3)') &
     &       -"MEMBER",IFRAME-1,-"VALUE",IOBJ,-"POSITION",TRA%T(ISTEP)
        ENDDO
        WRITE(NFIL,FMT='("#")')
        WRITE(NFIL,*)-'END'
      END IF
      CALL FILEHANDLER$CLOSE('MOVIE')
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL REPORT$CHVAL(NFILO,'FORMAT OF THE MOVIE FILE',FORMAT)
      CALL REPORT$I4VAL(NFILO,'#(FRAMES SKIPPED BETWEEN IMAGES)',ISKIP,' ')
      CALL REPORT$I4VAL(NFILO,'#(FRAMES)',IFRAME,' ')
      WRITE(NFILO,*)
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      DEALLOCATE(R0)
      DEALLOCATE(IZ0)
      DEALLOCATE(Q0)
      DEALLOCATE(ATOM0)
                             CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEXYZ(NFIL,FRAME,NAT,ID,R)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4)  ,INTENT(IN) :: FRAME
      INTEGER(4)  ,INTENT(IN) :: NAT
      CHARACTER(2),INTENT(IN) :: ID(NAT)
      REAL(8)     ,INTENT(IN) :: R(3,NAT)
      INTEGER(4)              :: IAT
      REAL(8)                 :: ANGSTROM
!     ******************************************************************
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
      IF(FRAME.EQ.1) REWIND NFIL
      WRITE(NFIL,*)NAT
      WRITE(NFIL,FMT='(A10,I10)')'NONAME',FRAME
      DO IAT=1,NAT
        WRITE(NFIL,FMT='(A2,2X,3(F10.5,1X))')ID(IAT),R(:,IAT)/ANGSTROM
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEEXTXYZ(NFIL,FRAME,NAT,ID,R,CELL,TIME)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE STRINGS_MODULE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4)  ,INTENT(IN) :: FRAME
      INTEGER(4)  ,INTENT(IN) :: NAT
      CHARACTER(2),INTENT(IN) :: ID(NAT)
      REAL(8)     ,INTENT(IN) :: R(3,NAT)
      REAL(8)     ,INTENT(IN) :: CELL(9)
      REAL(8)     ,INTENT(IN) :: TIME
      INTEGER(4)              :: IAT
      REAL(8)                 :: ANGSTROM, PICO, SECOND           
      CHARACTER(100)          :: STRING
      CHARACTER(200)          :: EXTXYZ
      CHARACTER(2)            :: SPECIES
!     ******************************************************************
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
      CALL CONSTANTS$GET('PICO', PICO)
      CALL CONSTANTS$GET('SECOND', SECOND)
      IF(FRAME.EQ.1) REWIND NFIL
      WRITE(NFIL,*)NAT
      WRITE(STRING,FMT='(9F10.5)')CELL/ANGSTROM
      EXTXYZ=+'L'//-'ATTICE="'//TRIM(ADJUSTL(STRING))//'" '//+'P'// &
     &       -'ROPERTIES=SPECIES:'//+'S'//-':1:POS:'//+'R:3 I'//-'TER='
      WRITE(STRING,FMT='(I10)')FRAME
      EXTXYZ=TRIM(EXTXYZ)//TRIM(ADJUSTL(STRING))//+' T'//-'IME='
      WRITE(STRING,FMT='(F10.5)')TIME/(PICO*SECOND)
      EXTXYZ=TRIM(EXTXYZ)//TRIM(ADJUSTL(STRING))
      WRITE(NFIL,*)EXTXYZ
      DO IAT=1,NAT
        ! TRANSFORM ELEMENT SYMBOL TO UPPERCASE+LOWERCASE
        SPECIES=+ID(IAT)
        SPECIES(2:2)=-SPECIES(2:2)
        WRITE(NFIL,FMT='(A2,2X,3(F10.5,1X))')SPECIES,R(:,IAT)/ANGSTROM
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MODEL$NATM(BOXR0,BOXVEC,NAT,POS,RBAS,NATM)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: BOXR0(3)
      REAL(8)   ,INTENT(IN) :: BOXVEC(3,3)
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: POS(3,NAT)
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      INTEGER(4),INTENT(OUT):: NATM
      INTEGER(4)            :: IAT,I1,I2,I3
      REAL(8)               :: BOXVECIN(3,3)
      REAL(8)               :: DET
      REAL(8)               :: R(3)
      REAL(8)               :: RBOXIN(3:3)
      REAL(8)               :: XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3
      INTEGER(4)            :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
      REAL(8)               :: T1,T2,T3
      LOGICAL(4)            :: TCHK
!     ******************************************************************
!     
!     ==================================================================
!     ==  INVERSE OF BOXVEC                                           ==
!     ==================================================================
      BOXVECIN(1,1)=BOXVEC(2,2)*BOXVEC(3,3)-BOXVEC(3,2)*BOXVEC(2,3)
      BOXVECIN(2,1)=BOXVEC(3,2)*BOXVEC(1,3)-BOXVEC(1,2)*BOXVEC(3,3)
      BOXVECIN(3,1)=BOXVEC(1,2)*BOXVEC(2,3)-BOXVEC(2,2)*BOXVEC(1,3)
      BOXVECIN(1,2)=BOXVEC(2,3)*BOXVEC(3,1)-BOXVEC(3,3)*BOXVEC(2,1)
      BOXVECIN(2,2)=BOXVEC(3,3)*BOXVEC(1,1)-BOXVEC(1,3)*BOXVEC(3,1)
      BOXVECIN(3,2)=BOXVEC(1,3)*BOXVEC(2,1)-BOXVEC(2,3)*BOXVEC(1,1)
      BOXVECIN(1,3)=BOXVEC(2,1)*BOXVEC(3,2)-BOXVEC(3,1)*BOXVEC(2,2)
      BOXVECIN(2,3)=BOXVEC(3,1)*BOXVEC(1,2)-BOXVEC(1,1)*BOXVEC(3,2)
      BOXVECIN(3,3)=BOXVEC(1,1)*BOXVEC(2,2)-BOXVEC(2,1)*BOXVEC(1,2)
      DET=DOT_PRODUCT(BOXVECIN(:,1),BOXVEC(:,1))
      BOXVECIN(:,:)=BOXVECIN(:,:)/DET
!     
!     ==================================================================
!     ==  CALCULATE #(ATOMS) IN THE MODEL                             ==
!     ==================================================================
      NATM=0
      DO IAT=1,NAT
        R(:)=POS(:,IAT)
        CALL BOXBOX(RBAS,BOXR0-R,BOXVEC,XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3)
        MIN1=NINT(XMIN1)-1
        MAX1=NINT(XMAX1)+1
        MIN2=NINT(XMIN2)-1
        MAX2=NINT(XMAX2)+1
        MIN3=NINT(XMIN3)-1
        MAX3=NINT(XMAX3)+1
        DO I1=MIN1,MAX1
          T1=REAL(I1,KIND=8)
          DO I2=MIN2,MAX2
            T2=REAL(I2,KIND=8)
            DO I3=MIN3,MAX3
              T3=REAL(I3,KIND=8)
              R=POS(:,IAT)+RBAS(:,1)*T1+RBAS(:,2)*T2+RBAS(:,3)*T3-BOXR0(:)
              R=MATMUL(BOXVECIN,R)
              TCHK=(R(1).GE.0.D0.AND.R(1).LE.1.D0).AND. &
                   (R(2).GE.0.D0.AND.R(2).LE.1.D0).AND. &
                   (R(3).GE.0.D0.AND.R(3).LE.1.D0)
              IF(TCHK) NATM=NATM+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MODEL$ATOMS(BOXR0,BOXVEC,NAT,POS,RBAS,NATM,POSM,MAP)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: BOXR0(3)
      REAL(8)   ,INTENT(IN) :: BOXVEC(3,3)
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: POS(3,NAT)
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      INTEGER(4),INTENT(IN) :: NATM
      REAL(8)   ,INTENT(OUT):: POSM(3,NATM)
      INTEGER(4),INTENT(OUT):: MAP(NATM)
      INTEGER(4)            :: IAT,IATM,I1,I2,I3
      REAL(8)               :: BOXVECIN(3,3)
      REAL(8)               :: DET
      REAL(8)               :: R(3),XR(3)
      REAL(8)               :: RBOXIN(3:3)
      REAL(8)               :: XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3
      INTEGER(4)            :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
      REAL(8)               :: T1,T2,T3
      LOGICAL(4)            :: TCHK
!     ******************************************************************
!     
!     ==================================================================
!     ==  INVERSE OF BOXVEC                                           ==
!     ==================================================================
      BOXVECIN(1,1)=BOXVEC(2,2)*BOXVEC(3,3)-BOXVEC(3,2)*BOXVEC(2,3)
      BOXVECIN(2,1)=BOXVEC(3,2)*BOXVEC(1,3)-BOXVEC(1,2)*BOXVEC(3,3)
      BOXVECIN(3,1)=BOXVEC(1,2)*BOXVEC(2,3)-BOXVEC(2,2)*BOXVEC(1,3)
      BOXVECIN(1,2)=BOXVEC(2,3)*BOXVEC(3,1)-BOXVEC(3,3)*BOXVEC(2,1)
      BOXVECIN(2,2)=BOXVEC(3,3)*BOXVEC(1,1)-BOXVEC(1,3)*BOXVEC(3,1)
      BOXVECIN(3,2)=BOXVEC(1,3)*BOXVEC(2,1)-BOXVEC(2,3)*BOXVEC(1,1)
      BOXVECIN(1,3)=BOXVEC(2,1)*BOXVEC(3,2)-BOXVEC(3,1)*BOXVEC(2,2)
      BOXVECIN(2,3)=BOXVEC(3,1)*BOXVEC(1,2)-BOXVEC(1,1)*BOXVEC(3,2)
      BOXVECIN(3,3)=BOXVEC(1,1)*BOXVEC(2,2)-BOXVEC(2,1)*BOXVEC(1,2)
      DET=DOT_PRODUCT(BOXVECIN(:,1),BOXVEC(:,1))
      BOXVECIN(:,:)=BOXVECIN(:,:)/DET
!     
!     ==================================================================
!     ==  CALCULATE #(ATOMS) IN THE MODEL                             ==
!     ==================================================================
      IATM=0
      DO IAT=1,NAT
        R(:)=POS(:,IAT)
        CALL BOXBOX(RBAS,BOXR0-R,BOXVEC,XMIN1,XMAX1,XMIN2,XMAX2,XMIN3,XMAX3)
        MIN1=NINT(XMIN1)-1
        MAX1=NINT(XMAX1)+1
        MIN2=NINT(XMIN2)-1
        MAX2=NINT(XMAX2)+1
        MIN3=NINT(XMIN3)-1
        MAX3=NINT(XMAX3)+1
        DO I1=MIN1,MAX1
          T1=REAL(I1,KIND=8)
          DO I2=MIN2,MAX2
            T2=REAL(I2,KIND=8)
            DO I3=MIN3,MAX3
              T3=REAL(I3,KIND=8)
              R=POS(:,IAT)+RBAS(:,1)*T1+RBAS(:,2)*T2+RBAS(:,3)*T3
              XR=MATMUL(BOXVECIN,R-BOXR0)
              TCHK=(XR(1).GE.0.D0.AND.XR(1).LE.1.D0).AND. &
                   (XR(2).GE.0.D0.AND.XR(2).LE.1.D0).AND. &
                   (XR(3).GE.0.D0.AND.XR(3).LE.1.D0)
              IF(TCHK) THEN
                IATM=IATM+1
                POSM(:,IATM)=R(:)
                MAP(IATM)=IAT
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MODEL$NBONDM(NAT,POS,RAD,NBOND)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: POS(3,NAT)
      REAL(8)   ,INTENT(IN) :: RAD(NAT)
      INTEGER(4),INTENT(OUT):: NBOND
      INTEGER(4)            :: IAT1,IAT2
      REAL(8)               :: DR(3)
      REAL(8)               :: DIS
!     ******************************************************************
      NBOND=0
      DO IAT1=1,NAT
        DO IAT2=IAT1+1,NAT
          DR(:)=POS(:,IAT1)-POS(:,IAT2)
          DIS=SQRT(DR(1)**2+DR(2)**2+DR(3)**2)
          IF(DIS.LT.RAD(IAT1)+RAD(IAT2)) NBOND=NBOND+1
        ENDDO
      ENDDO
      RETURN
      END
!     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MODEL$BONDS(NAT,POS,RAD,NBOND,BOND)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: POS(3,NAT)
      REAL(8)   ,INTENT(IN) :: RAD(NAT)
      INTEGER(4),INTENT(IN) :: NBOND
      INTEGER(4),INTENT(OUT):: BOND(2,NBOND)
      INTEGER(4)            :: IAT1,IAT2,IBOND
      REAL(8)               :: DR(3)
      REAL(8)               :: DIS
!     ******************************************************************
      IBOND=0
      DO IAT1=1,NAT
        DO IAT2=IAT1+1,NAT
          DR(:)=POS(:,IAT1)-POS(:,IAT2)
          DIS=SQRT(DOT_PRODUCT(DR,DR))
          IF(DIS.LT.RAD(IAT1)+RAD(IAT2)) THEN
            IBOND=IBOND+1
            BOND(1,IBOND)=IAT1
            BOND(2,IBOND)=IAT2
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE READ_STRC(LL_STRC_)
!      *************************************************************************
!      **  READS STRC FILE                                                    **
!      *************************************************************************
       USE LINKEDLIST_MODULE
       USE STRINGS_MODULE
       USE PERIODICTABLE_MODULE
       USE FORCEFIELD_MODULE, ONLY: MMATOM, PDB_FORM
       USE TRAJECTORY_MODULE, ONLY: NAT,IZ,MASS,ATOM,RBAS,QNAT,TQMMM,FFORMAT
       IMPLICIT NONE
       TYPE(LL_TYPE),INTENT(IN)  :: LL_STRC_
       TYPE(LL_TYPE)             :: LL_STRC
       CHARACTER(16),ALLOCATABLE :: SP(:)    !(NAT)
       CHARACTER(16)             :: SPNAME
       INTEGER(4)                :: I,IAT
       INTEGER(4)                :: NSP
       CHARACTER(2)              :: SYMBOL
       LOGICAL(4)                :: TCHK
       INTEGER(4)                :: IZ1
       REAL(8)                   :: MASS1
       REAL(8)                   :: MASSUNIT
       REAL(8)                   :: LUNIT
       INTEGER(4)                :: NFILO
INTEGER(4),ALLOCATABLE    :: IND(:)
INTEGER(4)                :: J,NATM
!      *************************************************************************
                               CALL TRACE$PUSH('READ_STRC')
       LL_STRC=LL_STRC_

       CALL CONSTANTS$GET('U',MASSUNIT)
!
!      =========================================================================
!      ==   GET NUMBER OF ATOMS AND ALLOCATE ARRAYS                           ==
!      =========================================================================
       IF(.NOT.TQMMM) THEN
          CALL LINKEDLIST$SELECT(LL_STRC,'~')
          CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
          CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',NAT)
          ALLOCATE(IZ(NAT))
          ALLOCATE(MASS(NAT))
          ALLOCATE(ATOM(NAT))
          ALLOCATE(SP(NAT))
       ELSE 
          IF(FFORMAT.EQ.'STRC') THEN
             CALL LINKEDLIST$SELECT(LL_STRC,'~')
             CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
             CALL LINKEDLIST$SELECT(LL_STRC,'QM-MM')
             CALL LINKEDLIST$NLISTS(LL_STRC,'ATOM',QNAT)
             CALL LINKEDLIST$SELECT(LL_STRC,'..')
          ELSE IF(FFORMAT.EQ.'PDB') THEN
             CALL FORCEFIELD$READ_MMSTRC
             CALL FORCEFIELD$GETI4('NMMATOM',QNAT)
! !-- REPORT ATOM DATA FROM PDB FILE
!          DO IAT=1,QNAT
!             WRITE(*,PDB_FORM) MMATOM(IAT)
!          ENDDO
! PRINT*,"FLAG: FORCED STOP - PRINTOUT OF MMATOMS. ",QNAT
! STOP         
         END IF
         ALLOCATE(IZ(QNAT))
         ALLOCATE(MASS(QNAT))
         ALLOCATE(ATOM(QNAT))
         ALLOCATE(SP(QNAT))
       END IF
       IZ(:)=0
       MASS(:)=0.D0
       ATOM(:)=' '
       SP(:)=' '
!
!      =========================================================================
!      ==   PICK UP ATOM SPECIES OF THE ATOMS                                 ==
!      =========================================================================
       CALL LINKEDLIST$SELECT(LL_STRC,'~')
       CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
       IF(.NOT.TQMMM) THEN
          CALL LINKEDLIST$SELECT(LL_STRC,'~')
          CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
          ALLOCATE(IND(NAT))
          DO I=1,NAT
             CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',I)
             CALL LINKEDLIST$GET(LL_STRC,'INDEX',1,IAT)
             IND(I)=IAT
             CALL LINKEDLIST$GET(LL_STRC,'NAME',1,ATOM(IAT))
             SP(IAT)=ATOM(IAT)(1:2)
             CALL LINKEDLIST$EXISTD(LL_STRC,'SP',1,TCHK)
             IF(TCHK)CALL LINKEDLIST$GET(LL_STRC,'SP',1,SP(IAT))
             CALL LINKEDLIST$SELECT(LL_STRC,'..')
          ENDDO
       ELSE 
          IF(FFORMAT.EQ.'STRC') THEN
             CALL LINKEDLIST$SELECT(LL_STRC,'~')
             CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
             ALLOCATE(IND(QNAT))
             CALL LINKEDLIST$SELECT(LL_STRC,'QM-MM')
             DO I=1,QNAT
                CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',I)
                CALL LINKEDLIST$EXISTD(LL_STRC,'INDEX',1,TCHK)
                IF(TCHK) THEN
                   CALL LINKEDLIST$GET(LL_STRC,'INDEX',1,IAT) !THERE IS NO INDEX IN QMMM PART
                ELSE 
                   IAT = I
                END IF
                IND(I)=IAT
                CALL LINKEDLIST$GET(LL_STRC,'NAME',1,ATOM(IAT))
                SP(IAT) = ATOM(IAT)(1:2)
                CALL LINKEDLIST$EXISTD(LL_STRC,'SP',1,TCHK)
                IF(TCHK)CALL LINKEDLIST$GET(LL_STRC,'SP',1,SP(IAT))
                CALL LINKEDLIST$SELECT(LL_STRC,'..')
             ENDDO
          END IF
          IF(FFORMAT.EQ.'PDB') THEN
             ALLOCATE(IND(QNAT))
             DO I=1,QNAT
                IAT=I
                IND(I)=IAT
                ATOM(IAT)=TRIM(ADJUSTL(MMATOM(IAT)%NAME))//'_'//.ITOS.IAT
               SP(IAT) = ATOM(IAT)(1:2)
            END DO
          END IF
       END IF



!      ==CHECK FOR IDENTICAL ATOM NAMES AND INDICES OUT OF RANGE ===============
       IF(.NOT.TQMMM) THEN
         DO I=1,NAT
            DO J=I+1,NAT
              IF(IND(I).EQ.IND(J)) THEN
                CALL ERROR$MSG('TWO ATOMS WITH THE SAME NAME ARE NOT ALLOWED')
                CALL ERROR$CHVAL('ATOM NAME',ATOM(IND(I)))
                CALL ERROR$I4VAL('POSITION OF FIRST ATOM IN STRC_OUT FILE ',I)
                CALL ERROR$I4VAL('POSITION OF SECOND ATOM IN STRC_OUT FILE ',J)
                CALL ERROR$STOP('READ_STRC')
              END IF
            ENDDO
         ENDDO
       END IF

       IF(ALLOCATED(IND)) DEALLOCATE(IND)
!
!      =========================================================================
!      ==   LOOK UP ATOMIC NUMBER OF SPECIES                                  ==
!      =========================================================================
       IF(.NOT.TQMMM.OR.FFORMAT.EQ.'STRC') THEN
          CALL LINKEDLIST$SELECT(LL_STRC,'~')
          CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
          CALL LINKEDLIST$NLISTS(LL_STRC,'SPECIES',NSP)
          DO I=1,NSP
             CALL LINKEDLIST$SELECT(LL_STRC,'SPECIES',I)
             CALL LINKEDLIST$GET(LL_STRC,'NAME',1,SPNAME)
             SYMBOL=SPNAME(1:2)
             PRINT*,"FLAG SPNAME: ",SPNAME
             IF(SYMBOL(2:2).EQ.'_')SYMBOL(2:2)=' '
             CALL PERIODICTABLE$GET(SYMBOL,'Z',IZ1)
             !
             CALL PERIODICTABLE$GET(SYMBOL,'MASS',MASS1)
             CALL LINKEDLIST$EXISTD(LL_STRC,'M',1,TCHK)
             IF(TCHK) THEN
                CALL LINKEDLIST$GET(LL_STRC,'M',1,MASS1)
                MASS1=MASS1*MASSUNIT
             END IF
             IF(.NOT.TQMMM) THEN
                DO IAT=1,NAT
                   IF(SPNAME.EQ.SP(IAT)) THEN
                      IZ(IAT)=IZ1
                      MASS(IAT)=MASS1
                   END IF
                ENDDO
             ELSE
                DO IAT=1,QNAT
                   IF(SPNAME.EQ.SP(IAT)) THEN
                      IZ(IAT)=IZ1
                      MASS(IAT)=MASS1
                   END IF
                ENDDO
             END IF
             CALL LINKEDLIST$SELECT(LL_STRC,'..')
          ENDDO
       ELSE  !FOR PDB FILES
          DO IAT=1, QNAT
             SYMBOL = TRIM(ADJUSTL(MMATOM(IAT)%ELEMENT))
             CALL PERIODICTABLE$GET(SYMBOL,'Z',IZ1)
             CALL PERIODICTABLE$GET(SYMBOL,'MASS',MASS1)
             IZ(IAT)=IZ1
             MASS(IAT)=MASS1
          END DO
       END IF
!
!      =========================================================================
!      ==  OVERWRITE DEFAULTS                                                 ==
!      =========================================================================
       IF(.NOT.TQMMM.OR.FFORMAT.EQ.'STRC') THEN
          CALL LINKEDLIST$SELECT(LL_STRC,'~')
          CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')

          IF(.NOT.TQMMM) THEN
             DO I=1,NAT
                CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',I)
                CALL LINKEDLIST$GET(LL_STRC,'INDEX',1,IAT)
                CALL LINKEDLIST$EXISTD(LL_STRC,'M',1,TCHK)
                IF(TCHK)THEN
                   CALL LINKEDLIST$GET(LL_STRC,'M',1,MASS(IAT))
                   MASS(IAT)=MASS(IAT)*MASSUNIT
                END IF
                CALL LINKEDLIST$SELECT(LL_STRC,'..')
             ENDDO
          ELSE
             DO I=1,QNAT
                CALL LINKEDLIST$SELECT(LL_STRC,'QM-MM')
                CALL LINKEDLIST$SELECT(LL_STRC,'ATOM',I)
!             CALL LINKEDLIST$GET(LL_STRC,'INDEX',1,IAT)  !NO INDEX IN QMMM PART
                IAT = I
                CALL LINKEDLIST$EXISTD(LL_STRC,'M',1,TCHK)
                IF(TCHK)THEN
                   CALL LINKEDLIST$GET(LL_STRC,'M',1,MASS(IAT))
                   MASS(IAT)=MASS(IAT)*MASSUNIT
                END IF
                CALL LINKEDLIST$SELECT(LL_STRC,'..')
                CALL LINKEDLIST$SELECT(LL_STRC,'..')
             ENDDO
          END IF
       END IF
!
!      =========================================================================
!      ==   CHECK                                                             ==
!      =========================================================================

       IF(.NOT.TQMMM) THEN
          DO I=1,NAT
             IF(IZ(I).EQ.0) THEN
                CALL ERROR$MSG('ATOMIC NUMBER NOT FOUND')
                CALL ERROR$I4VAL('IAT',I)
                CALL ERROR$CHVAL('NAME',ATOM(I))
                CALL ERROR$STOP('READ_STRC')
             END IF
          ENDDO
       ELSE
          DO I=1,NAT
             IF(IZ(I).EQ.0) THEN
                CALL ERROR$MSG('ATOMIC NUMBER NOT FOUND')
                CALL ERROR$I4VAL('IAT',I)
                CALL ERROR$CHVAL('NAME',ATOM(I))
                CALL ERROR$STOP('READ_STRC')
             END IF
          ENDDO
       END IF

!
!     ==========================================================================
!     ==  GET UNIT CELL                                                       ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$SELECT(LL_STRC,'GENERIC')
      CALL LINKEDLIST$EXISTD(LL_STRC,'LUNIT',1,TCHK)
      LUNIT=1.D0
      IF(TCHK)CALL LINKEDLIST$GET(LL_STRC,'LUNIT',1,LUNIT)
      CALL LINKEDLIST$SELECT(LL_STRC,'..')
      CALL LINKEDLIST$SELECT(LL_STRC,'LATTICE')
      CALL LINKEDLIST$GET(LL_STRC,'T',1,RBAS)
      RBAS=RBAS*LUNIT
!
!     ==========================================================================
!     ==  REPORT DATA                                                         ==
!     ==========================================================================
     CALL FILEHANDLER$UNIT('PROT',NFILO)
     CALL REPORT$TITLE(NFILO,'DATA FROM STRUCTURE FILE')
     WRITE(NFILO,FMT='(10A8)')ATOM
                                  CALL TRACE$POP
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DXBALLSTICK(NFIL,IOBJECT0,NAT,COLOR,RAD,POS &
     &                      ,NBOND,IBOND,BOXR0,BOXVEC)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: COLOR(3,NAT)
      REAL(8)   ,INTENT(IN) :: RAD(NAT)
      REAL(8)   ,INTENT(IN) :: POS(3,NAT)
      INTEGER(4),INTENT(IN) :: NBOND
      INTEGER(4),INTENT(IN) :: IBOND(2,NBOND)
      REAL(8)   ,INTENT(IN) :: BOXR0(3)
      REAL(8)   ,INTENT(IN) :: BOXVEC(3,3)
      INTEGER(4),INTENT(IN) :: IOBJECT0
      CHARACTER(256)        :: LINE
      INTEGER(4)            :: I1,I2,I3
      REAL(8)               :: T1,T2,T3
      REAL(8)               :: X,Y,Z
!     **************************************************************************
!     
!     ==========================================================================
!     ==   DATA ARRAY: SIZE OF THE SPHERES                                    ==
!     ==========================================================================
!     WRITE(NFIL,FMT='("#"/"#",T10,-"SPHERE SIZE"/"#")')
      WRITE(NFIL,FMT='(A,1X,I10)')-'OBJECT',IOBJECT0+1
      WRITE(NFIL,FMT='(A,1X,I10)')-"CLASS ARRAY TYPE FLOAT RANK 0 ITEMS ",NAT
      WRITE(NFIL,FMT='(A)')      -"DATA FOLLOWS"
      WRITE(NFIL,FMT='(10F10.5)')RAD(:)
      WRITE(NFIL,FMT='(A)')      -"ATTRIBUTE ""DEP"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='("#")')
!     
!     ==========================================================================
!     ==   COLORS                                                             ==
!     ==========================================================================
!     WRITE(NFIL,FMT='("#"/"#",T10,"COLORS"/"#")')
      WRITE(NFIL,FMT='(A,1X,I10)')-'OBJECT',IOBJECT0+2
      WRITE(NFIL,FMT='(A,1X,I10)') &
     &                        -"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS",NAT
      WRITE(NFIL,FMT='(A)')-"DATA FOLLOWS"
      WRITE(NFIL,FMT='(10F10.5)')COLOR(:,:)
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""DEP"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='("#")')
!     
!     ==========================================================================
!     ==   POSITIONS ARRAY: ATOMIC POSITIONS                                  ==
!     ==========================================================================
!     WRITE(NFIL,FMT='("#"/"#",T10,-"ATOMIC POSITIONS"/"#")')
      WRITE(NFIL,FMT='(A,1X,I10)')-"OBJECT ",IOBJECT0+3
      WRITE(NFIL,FMT='(A,1X,I10)') &
     &                       -"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS ",NAT
      WRITE(NFIL,FMT='(A)')      -"DATA FOLLOWS"
      WRITE(NFIL,FMT='(10F12.5)')POS(:,:)
      WRITE(NFIL,FMT='(A)')      -"ATTRIBUTE ""DEP"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='("#")')
!     
!     ==========================================================================
!     ==   CONNECTIONS ARRAY: BONDS                                           ==
!     ==========================================================================
      IF(NBOND.GT.0) THEN
!     WRITE(NFIL,FMT='("#"/"#",T10,"BONDS"/"#")')
      WRITE(NFIL,FMT='(A,1X,I10)')-"OBJECT ",IOBJECT0+4
      WRITE(NFIL,FMT='(A,1X,I10)') &
     &                       -"CLASS ARRAY TYPE INT RANK 1 SHAPE 2 ITEMS ",NBOND
      WRITE(NFIL,FMT='(A)')      -"DATA FOLLOWS"
      WRITE(NFIL,FMT='(10I5)')   IBOND(:,:)-1
      WRITE(NFIL,FMT='(A)')      -"ATTRIBUTE ""REF"" STRING ""POSITIONS"""
      WRITE(NFIL,FMT='(A)')      -"ATTRIBUTE ""ELEMENT TYPE"" STRING ""LINES"""
      WRITE(NFIL,FMT='("#")')
      END IF
!     
!     ==========================================================================
!     ==   BOX: LATTICE VECTORS                                               ==
!     ==========================================================================
!     WRITE(NFIL,FMT='("#"/"#",T10,-"LATTICE VECTORS"/"#")')
      WRITE(NFIL,FMT='(A,I10)')-"OBJECT ",IOBJECT0+5
      WRITE(NFIL,FMT='(A)')-"CLASS ARRAY TYPE FLOAT RANK 1 SHAPE 3 ITEMS 8"
      WRITE(NFIL,FMT='(A)')      -"DATA FOLLOWS"
      DO I3=0,1 
         T3=DBLE(I3)
        DO I2=0,1
          T2=DBLE(I2)
          DO I1=0,1
            T1=DBLE(I1)
            X=BOXR0(1)+BOXVEC(1,1)*T1+BOXVEC(1,2)*T2+BOXVEC(1,3)*T3
            Y=BOXR0(1)+BOXVEC(2,1)*T1+BOXVEC(2,2)*T2+BOXVEC(2,3)*T3
            Z=BOXR0(1)+BOXVEC(3,1)*T1+BOXVEC(3,2)*T2+BOXVEC(3,3)*T3
            WRITE(NFIL,FMT='(3F15.5)')X,Y,Z
          ENDDO
        ENDDO
      ENDDO
      WRITE(NFIL,FMT='("#")')
!     
!     ==========================================================================
!     ==   OBJECT MOLECULE:                                                   ==
!     ==========================================================================
      WRITE(NFIL,FMT='(A,1X,I10)')-"OBJECT",IOBJECT0+6
      WRITE(NFIL,FMT='(A)')      -"CLASS FIELD"
      WRITE(NFIL,FMT='(A,1X,I10)')-"COMPONENT ""DATA"" VALUE ",IOBJECT0+1
      WRITE(NFIL,FMT='(A,1X,I10)')-"COMPONENT ""POSITIONS"" VALUE ",IOBJECT0+3
      IF(NBOND.GT.0) THEN
      WRITE(NFIL,FMT='(A,1X,I10)')-"COMPONENT ""CONNECTIONS"" VALUE ",IOBJECT0+4
      END IF
      WRITE(NFIL,FMT='(A,1X,I10)')-"COMPONENT ""BOX"" VALUE ",IOBJECT0+5
      WRITE(NFIL,FMT='(A,1X,I10)')-"COMPONENT ""COLORS"" VALUE ",IOBJECT0+2
      WRITE(NFIL,FMT='(A)')-"ATTRIBUTE ""NAME"" STRING ""ATOMS"""
      WRITE(NFIL,FMT='("#")')
      RETURN
      END
!        
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMCOLOR(IZ,ICOLOR)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IZ
      INTEGER(4),INTENT(OUT) :: ICOLOR(3)
      INTEGER(4)             :: ICOLORSTANDARD(3,106)
      INTEGER(4)             :: I
!     ******************************************************************
      DATA (ICOLORSTANDARD(I,  1),I=1,3)/135,125,131/
      DATA (ICOLORSTANDARD(I,  2),I=1,3)/255,228,196/!BISQUE
      DATA (ICOLORSTANDARD(I,  3),I=1,3)/240,248,255/!ALICE BLUE
      DATA (ICOLORSTANDARD(I,  4),I=1,3)/200,200,  0/!YELLOW
      DATA (ICOLORSTANDARD(I,  5),I=1,3)/192,102,624/!CORAL
      DATA (ICOLORSTANDARD(I,  6),I=1,3)/ 50, 50, 50/
      DATA (ICOLORSTANDARD(I,  7),I=1,3)/ 10,200, 10/
      DATA (ICOLORSTANDARD(I,  8),I=1,3)/255,  0,  0/!RED
      DATA (ICOLORSTANDARD(I,  9),I=1,3)/ 50,205, 50/!LIME GREEN
      DATA (ICOLORSTANDARD(I, 10),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 11),I=1,3)/ 60,  1,  1/!REDISH BLACK
      DATA (ICOLORSTANDARD(I, 12),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 13),I=1,3)/180,  1,100/!PURPLE
      DATA (ICOLORSTANDARD(I, 14),I=1,3)/  1, 20,198/!DARK BLUE
      DATA (ICOLORSTANDARD(I, 15),I=1,3)/230,171, 17/
      DATA (ICOLORSTANDARD(I, 16),I=1,3)/240,240,  0/
      DATA (ICOLORSTANDARD(I, 17),I=1,3)/ 60,180,  0/!YELLOW GREEN
      DATA (ICOLORSTANDARD(I, 18),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 19),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 20),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 21),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 22),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 23),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 24),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 25),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 26),I=1,3)/176, 48, 96/
      DATA (ICOLORSTANDARD(I, 27),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 28),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 29),I=1,3)/278,134, 34/!FIREBRICK
      DATA (ICOLORSTANDARD(I, 30),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 31),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 32),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 33),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 34),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 35),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 36),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 37),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 38),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 39),I=1,3)/ 40, 40,140/
      DATA (ICOLORSTANDARD(I, 40),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 41),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 42),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 43),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 44),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 45),I=1,3)/230, 51, 41/
      DATA (ICOLORSTANDARD(I, 46),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 47),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 48),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 49),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 50),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 51),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 52),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 53),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 54),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 55),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 56),I=1,3)/100, 40,  0/
      DATA (ICOLORSTANDARD(I, 57),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 58),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 59),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 60),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 61),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 62),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 63),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 64),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 65),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 66),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 67),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 68),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 69),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 70),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 71),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 72),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 73),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 74),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 75),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 76),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 77),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 78),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 79),I=1,3)/255,215,  0/!GOLD
      DATA (ICOLORSTANDARD(I, 80),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 81),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 82),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 83),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 84),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 85),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 86),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 87),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 88),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 89),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 90),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 91),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 92),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 93),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 94),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 95),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 96),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 97),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 98),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I, 99),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,100),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,101),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,102),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,103),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,104),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,105),I=1,3)/  1,  1,  1/
      DATA (ICOLORSTANDARD(I,106),I=1,3)/  1,  1,  1/   
      DO I=1,3
        ICOLOR(I)=ICOLORSTANDARD(I,IZ)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BOND(IAT1,IAT2,ARRAY)
!     ******************************************************************
!     **   CALCULATE   BOND DISTANCE                                 ***
!     ******************************************************************
      USE TRAJECTORY_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT1
      INTEGER(4),INTENT(IN)  :: IAT2
      REAL(8)   ,INTENT(OUT) :: ARRAY(NSTEP)
      INTEGER(4)             :: ISTEP,I,J,K
      REAL(8)                :: DR(3)
      REAL(8)                :: DIS
      REAL(8)                :: RBASINV(3,3)
      INTEGER(4),PARAMETER   :: NAUX=98
      REAL(8)                :: RCOND,DET(2),AUX(NAUX)
      REAL(8)                :: DR1(3),DIS1,DT1(3)
      REAL(8)                :: DT0(3)
!     ******************************************************************
!RBASINV=RBAS
!CALL DGEICD(RBASINV,3,3,0,RCOND,DET,AUX,NAUX) !ESSL MATRIX INVERSION
      CALL LIB$INVERTR8(3,RBAS,RBASINV)
      DT0(:)=0.D0
      DO ISTEP=1,NSTEP
         DR=TRA%R(:,IAT1,ISTEP)-TRA%R(:,IAT2,ISTEP)+DT0(:)
         DIS=SQRT(DR(1)**2+DR(2)**2+DR(3)**2)
!        == LOOK FOR SHORTEST DISTANCE OF ALL PERIODIC IMAGES ==========
 1000    CONTINUE
         DO I=-1,1
           DO J=-1,1
             DO K=-1,1
               DT1=RBAS(:,1)*DBLE(I)+RBAS(:,2)*DBLE(J)+RBAS(:,3)*DBLE(K)
               DR1=DR+DT1
               DIS1=SQRT(DOT_PRODUCT(DR1,DR1))
               IF(DIS1.LT.DIS) THEN
                 DT0=DT0+DT1
                 DR=DR1
                 DIS=DIS1
                 GOTO 1000
               END IF
             ENDDO
           ENDDO
         ENDDO
!        == CALCULATE DISTANCE =========================================
         DIS=SQRT(DR(1)**2+DR(2)**2+DR(3)**2)
         ARRAY(ISTEP)=DIS
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ANGLE(IAT1,IAT2,IAT3,ARRAY)
!     ******************************************************************
!     **   CALCULATE   BOND ANGLE                                     **
!     ******************************************************************
      USE TRAJECTORY_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT1
      INTEGER(4),INTENT(IN)  :: IAT2
      INTEGER(4),INTENT(IN)  :: IAT3
      REAL(8)   ,INTENT(OUT) :: ARRAY(NSTEP)
      INTEGER(4)             :: ISTEP,I
      REAL(8)                :: DR1(3)
      REAL(8)                :: DR2(3)
      REAL(8)                :: D11,D22,D12
      REAL(8)                :: SVAR
      REAL(8)                :: RBASINV(3,3)
      INTEGER(4),PARAMETER   :: NAUX=98
      REAL(8)                :: RCOND,DET(2),AUX(NAUX)
!     ******************************************************************
!RBASINV=RBAS
!CALL DGEICD(RBASINV,3,3,0,RCOND,DET,AUX,NAUX) !ESSL MATRIX INVERSION
      CALL LIB$INVERTR8(3,RBAS,RBASINV)
      DO ISTEP=1,NSTEP
         DR1=TRA%R(:,IAT3,ISTEP)-TRA%R(:,IAT2,ISTEP)
         DR2=TRA%R(:,IAT1,ISTEP)-TRA%R(:,IAT2,ISTEP)
!        
!        == LOOK FOR SHORTEST DISTANCE OF ALL PERIODIC IMAGES ==========
         DR1=MATMUL(RBASINV,DR1)
         DR2=MATMUL(RBASINV,DR2)
         DO I=1,3
           DR1(I)=MODULO(DR1(I)+0.5D0,1.D0)-0.5D0
           DR2(I)=MODULO(DR2(I)+0.5D0,1.D0)-0.5D0
         ENDDO
         DR1=MATMUL(RBAS,DR1)
         DR2=MATMUL(RBAS,DR2)
!        == CALCULATE ANGLE ==============================================
         D11=DR1(1)**2+DR1(2)**2+DR1(3)**2
         D22=DR2(1)**2+DR2(2)**2+DR2(3)**2
         D12=DR1(1)*DR2(1)+DR1(2)*DR2(2)+DR1(3)*DR2(3)
         SVAR=D12/SQRT(D11*D22)
         ARRAY(ISTEP)=ACOS(SVAR)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TORSION(IAT1,IAT2,IAT3,IAT4,ARRAY)
!     ******************************************************************
!     **   CALCULATE BOND TORSION                                     **
!     **                                                              **
!     **   CALCULATE BOND TORSION                                     **
!     **   CALCULATE BOND TORSION                                     **
!     ******************************************************************
      USE TRAJECTORY_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT1   ! FIRST ATOM 
      INTEGER(4),INTENT(IN)  :: IAT2   ! 
      INTEGER(4),INTENT(IN)  :: IAT3
      INTEGER(4),INTENT(IN)  :: IAT4
      REAL(8)   ,INTENT(OUT) :: ARRAY(NSTEP)
      INTEGER(4)             :: ISTEP,I
      REAL(8)                :: DR12(3)
      REAL(8)                :: DR23(3)
      REAL(8)                :: DR34(3)
      REAL(8)                :: V1(3)
      REAL(8)                :: V2(3)
      REAL(8)                :: D11,D22,D12
      REAL(8)                :: SVAR
      REAL(8)                :: PREV
      REAL(8)                :: PI
      INTEGER(4)             :: N
      REAL(8)                :: RBASINV(3,3)
      INTEGER(4),PARAMETER   :: NAUX=98
      REAL(8)                :: RCOND,DET(2),AUX(NAUX)
!     ******************************************************************
!RBASINV=RBAS
!CALL DGEICD(RBASINV,3,3,0,RCOND,DET,AUX,NAUX) !ESSL MATRIX INVERSION
      CALL LIB$INVERTR8(3,RBAS,RBASINV)
      PI=4.D0*ATAN(1.D0)
      PREV=0.D0
      DO ISTEP=1,NSTEP
         DR12=TRA%R(:,IAT2,ISTEP)-TRA%R(:,IAT1,ISTEP)
         DR23=TRA%R(:,IAT3,ISTEP)-TRA%R(:,IAT2,ISTEP)
         DR34=TRA%R(:,IAT4,ISTEP)-TRA%R(:,IAT3,ISTEP)
!        
!        == LOOK FOR SHORTEST DISTANCE OF ALL PERIODIC IMAGES ==========
         DR12=MATMUL(RBASINV,DR12)
         DR23=MATMUL(RBASINV,DR23)
         DR34=MATMUL(RBASINV,DR34)
         DO I=1,3
           DR12(I)=MODULO(DR12(I)+0.5D0,1.D0)-0.5D0
           DR23(I)=MODULO(DR23(I)+0.5D0,1.D0)-0.5D0
           DR34(I)=MODULO(DR34(I)+0.5D0,1.D0)-0.5D0
         ENDDO
         DR12=MATMUL(RBAS,DR12)
         DR23=MATMUL(RBAS,DR23)
         DR34=MATMUL(RBAS,DR34)
!        == CALCULATE ANGLE ==============================================
         V1(1)=DR12(2)*DR23(3)-DR12(3)*DR23(2)
         V1(2)=DR12(3)*DR23(1)-DR12(1)*DR23(3)
         V1(3)=DR12(1)*DR23(2)-DR12(2)*DR23(1)
         V2(1)=DR23(2)*DR34(3)-DR23(3)*DR34(2)
         V2(2)=DR23(3)*DR34(1)-DR23(1)*DR34(3)
         V2(3)=DR23(1)*DR34(2)-DR23(2)*DR34(1)
         D11=V1(1)**2+V1(2)**2+V1(3)**2
         D22=V2(1)**2+V2(2)**2+V2(3)**2
         D12=V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
         SVAR=D12/SQRT(D11*D22)
         ARRAY(ISTEP)=ACOS(SVAR)
         SVAR=DR23(1)*(V1(2)*V2(3)-V1(3)*V2(2)) &
        &    +DR23(2)*(V1(3)*V2(1)-V1(1)*V2(3)) &
        &    +DR23(3)*(V1(1)*V2(2)-V1(2)*V2(1))
         IF(SVAR.LT.0.D0)ARRAY(ISTEP)=-ARRAY(ISTEP)
         N=NINT((ARRAY(ISTEP)-PREV)/(2.D0*PI))
         ARRAY(ISTEP)=ARRAY(ISTEP)+2.D0*PI*DBLE(N)
         PREV=ARRAY(ISTEP)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE VELOCITY(N,T,X,V)
!     ******************************************************************
!     **   CALCULATE TIME DERIVATIVE OF A COORDINATE TRAJECTORY       **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N
      REAL(8)   ,INTENT(IN)  :: T(N)
      REAL(8)   ,INTENT(IN)  :: X(N)
      REAL(8)   ,INTENT(OUT) :: V(N)
      INTEGER(4)             :: I
      REAL(8)                :: XP,XM,TP,TM
!     ******************************************************************
      DO I=2,N-1
        V(I)=0.D0
        TP=T(I+1)-T(I)
        IF(TP.EQ.0.D0) CYCLE
        TM=T(I-1)-T(I)
        IF(TM.EQ.0.D0) CYCLE
        XP=X(I+1)-X(I)
        XM=X(I-1)-X(I)
        V(I)=((TP/TM)*XM-(TM/TP)*XP)/(TP-TM)
      ENDDO
      IF(T(2)-T(1).NE.0) V(1)=(X(2)-X(1))/(T(2)-T(1))
      IF(T(N)-T(N-1).NE.0)V(N)=(X(N)-X(N-1))/(T(N)-T(N-1))
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NEIGHBORS(LL_CNTL_)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES NEIGHBORLIST AND ITS CHANGES DURING THE          **
!     **  TRAJECTORY.                                                 **
!     **                                                              **
!     ******************************************************************
      USE TRAJECTORY_MODULE
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      REAL(8)   ,PARAMETER  :: DISTX=1.2D0
      REAL(8)               :: DTMAX
      INTEGER(4)            :: NBX
      INTEGER(4)            :: NB(NSTEP)
      TYPE(BOND_TYPE),ALLOCATABLE:: BOND(:,:)
      REAL(8)               :: RCOV(NAT)
      REAL(8)               :: R(3,NAT)
      INTEGER(4),PARAMETER  :: NPX=8
      INTEGER(4),PARAMETER  :: NSKIP=1
      INTEGER(4)            :: IPOINT(NPX)
      INTEGER(4)            :: NPOINT
      INTEGER(4)            :: ISTEP,IAT,I,IB
      INTEGER(4)            :: ISTEP2,IAT2,IAT1
      LOGICAL(4)            :: TCHK
      REAL(8)               :: PICO,SECOND
      INTEGER(4)            :: IBRAKE,IMAKE
      INTEGER(4)            :: NFIL
      TYPE JUMP_TYPE
        LOGICAL(4)             :: ON
        TYPE(BOND_TYPE)        :: BOND
        LOGICAL(4)             :: FORMING
        INTEGER(4)             :: ISTEP
        REAL(8)                :: T
        TYPE(JUMP_TYPE),POINTER:: NEXT
      END TYPE JUMP_TYPE
      TYPE(JUMP_TYPE),POINTER  :: FIRSTJUMP
      TYPE(JUMP_TYPE),POINTER  :: JUMP
      TYPE(JUMP_TYPE),POINTER  :: JUMP1,NEXT
      CHARACTER(256)           :: FILENAME
!     **************************************************************************
      LL_CNTL=LL_CNTL_
      CALL CONSTANTS$GET('PICO',PICO)
      CALL CONSTANTS$GET('SECOND',SECOND)
!
!     ==========================================================================
!     ==  READ CONTROL FILE                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'NEIGHBORS',1,TCHK)
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_CNTL,'NEIGHBORS')
      CALL FILEHANDLER$SETFILE('NEIGHBORS',.TRUE.,-'.TRA.NEIGHBORS')
      CALL FILEHANDLER$SETSPECIFICATION('NEIGHBORS','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('NEIGHBORS','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('NEIGHBORS','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('NEIGHBORS','FORM','FORMATTED')
      CALL FILEHANDLER$UNIT('NEIGHBORS',NFIL)
!
!     ==========================================================================
!     == CALCULATE RADII                                                      ==
!     ==========================================================================
      DO IAT=1,NAT
        CALL PERIODICTABLE$GET(IZ(IAT),'R(COV)',RCOV(IAT))
      ENDDO
!
!     ==========================================================================
!     == COLLECT POTENTIAL NEIGHBORS                                          ==
!     ==========================================================================
PRINT*,'BONDS'
      NBX=NAT*NPX
      ALLOCATE(BOND(NBX,NSTEP))
      DO ISTEP=1,NSTEP,NSKIP
        R(:,:)=TRA%R(:,:,ISTEP)
        CALL BONDS_EVAL(RBAS,NAT,R,RCOV,DISTX,NBX,NB(ISTEP),BOND(:,ISTEP))
      ENDDO
!
!     ==========================================================================
!     == WRITE NEIGHBORHOOD                                                   ==
!     ==========================================================================
      WRITE(NFIL,FMT='("NEIGHBOR TABLE"/14("="))')
      ISTEP=1
      DO IAT=1,NAT
        CALL BONDS_NEIGHBORS(NFIL,IAT,NB(ISTEP),BOND(:,ISTEP))
      ENDDO
!
!     ==========================================================================
!     == COLLECT JUMPS                                                        ==
!     ==========================================================================
PRINT*,'COLLECT JUMPS'
      NULLIFY(FIRSTJUMP)
      NULLIFY(JUMP)
      DO ISTEP=1,NSTEP
!
!       ========================================================================
!       == LOOK FOR BOND FORMATION                                            ==
!       ========================================================================
        IF(ISTEP.GT.1) THEN
          DO IB=1,NB(ISTEP)
            CALL BONDS_FIND(BOND(IB,ISTEP),NB(ISTEP-1),BOND(1,ISTEP-1),I)
            IF(I.EQ.0) THEN
              IF(ASSOCIATED(JUMP)) THEN
                ALLOCATE(JUMP%NEXT)
                JUMP=>JUMP%NEXT
              ELSE
                ALLOCATE(FIRSTJUMP)
                JUMP=>FIRSTJUMP
              END IF
              NULLIFY(JUMP%NEXT)
              JUMP%ON=.TRUE.
              JUMP%BOND=BOND(IB,ISTEP)
              JUMP%FORMING=.TRUE.
              JUMP%ISTEP=ISTEP
              JUMP%T=TRA%T(ISTEP)
            END IF
          ENDDO
        END IF
!
!       ========================================================================
!       == LOOK FOR BOND BREAKING                                             ==
!       ========================================================================
        IF(ISTEP.LT.NSTEP) THEN
          DO IB=1,NB(ISTEP)
            CALL BONDS_FIND(BOND(IB,ISTEP),NB(ISTEP+1),BOND(1,ISTEP+1),I)
            IF(I.EQ.0) THEN
              IF(ASSOCIATED(JUMP)) THEN
                ALLOCATE(JUMP%NEXT)
                JUMP=>JUMP%NEXT
              ELSE
                ALLOCATE(FIRSTJUMP)
                JUMP=>FIRSTJUMP
              END IF
              NULLIFY(JUMP%NEXT)
              JUMP%ON=.TRUE.
              JUMP%BOND=BOND(IB,ISTEP)
              JUMP%FORMING=.FALSE.
              JUMP%ISTEP=ISTEP
              JUMP%T=TRA%T(ISTEP)
            END IF
          ENDDO        
        END IF
      ENDDO
!
!     ==========================================================================
!     == CLEAN JUMPS                                                          ==
!     ==========================================================================
!PRINT*,'COLLAPSING'
      IF(.NOT.ASSOCIATED(FIRSTJUMP)) GOTO 9999
      DTMAX=0.1D0*PICO*SECOND
      JUMP=>FIRSTJUMP
      DO 
        IF(.NOT.JUMP%ON) THEN
          IF(.NOT.ASSOCIATED(JUMP%NEXT)) EXIT 
          JUMP=>JUMP%NEXT
          CYCLE
        END IF
!
!       == FIND AND DELETE CORRELATED RETURN JUMPS =====================
        IAT1=JUMP%BOND%IAT1
        IAT2=JUMP%BOND%IAT2
        JUMP1=>JUMP
        LOOP2: DO WHILE (ASSOCIATED(JUMP1%NEXT))
          JUMP1=>JUMP1%NEXT
          IF(.NOT.JUMP1%ON)CYCLE
          IF(JUMP1%T-JUMP%T.GT.DTMAX) EXIT LOOP2
          IF(JUMP1%FORMING.EQV.JUMP%FORMING) CYCLE
          IF(JUMP1%BOND%IAT1.EQ.IAT1) THEN
            IF(JUMP1%BOND%IAT2.NE.IAT2) CYCLE
            IF(JUMP1%BOND%IT(1).NE.JUMP%BOND%IT(1)) CYCLE
            IF(JUMP1%BOND%IT(2).NE.JUMP%BOND%IT(2)) CYCLE
            IF(JUMP1%BOND%IT(3).NE.JUMP%BOND%IT(3)) CYCLE
          ELSE IF(JUMP1%BOND%IAT1.EQ.IAT2) THEN
            IF(JUMP1%BOND%IAT2.NE.IAT1) CYCLE
            IF(JUMP1%BOND%IT(1).NE.-JUMP%BOND%IT(1)) CYCLE
            IF(JUMP1%BOND%IT(2).NE.-JUMP%BOND%IT(2)) CYCLE
            IF(JUMP1%BOND%IT(3).NE.-JUMP%BOND%IT(3)) CYCLE
          ELSE
            CYCLE
          END IF
          JUMP%ON=.FALSE.
          JUMP1%ON=.FALSE.
          IF(JUMP%FORMING) THEN
            DO ISTEP2=JUMP%ISTEP,JUMP1%ISTEP
              CALL BONDS_DELETEBOND(JUMP%BOND,NB(ISTEP2),BOND(1,ISTEP2))
            ENDDO
          ELSE
            DO ISTEP2=JUMP%ISTEP+1,JUMP1%ISTEP-1
              CALL BONDS_ADDBOND(JUMP%BOND,RBAS,NAT,R,RCOV,NBX,NB(ISTEP2),BOND(1,ISTEP2))
            ENDDO
          END IF
          EXIT LOOP2
        ENDDO LOOP2
        IF(.NOT.ASSOCIATED(JUMP%NEXT)) EXIT
        JUMP=>JUMP%NEXT
      ENDDO 
!
!PRINT*,'MARKE X1'
!GOTO 1000
      JUMP=>FIRSTJUMP
      DO WHILE(ASSOCIATED(JUMP%NEXT))
        JUMP1=>JUMP%NEXT
        IF(JUMP1%ON) THEN
          JUMP=>JUMP1
        ELSE
          IF(ASSOCIATED(JUMP1%NEXT)) THEN
            JUMP%NEXT=>JUMP1%NEXT
          ELSE            
            NULLIFY(JUMP%NEXT)
          END IF          
          DEALLOCATE(JUMP1)
        END IF
      ENDDO
!
!PRINT*,'MARKE X2'
      IF(.NOT.FIRSTJUMP%ON) THEN
        IF(ASSOCIATED(FIRSTJUMP%NEXT)) THEN
          JUMP=>FIRSTJUMP
          FIRSTJUMP=>FIRSTJUMP%NEXT
          DEALLOCATE(JUMP)
        ELSE
          GOTO 9999
        END IF
      END IF
 1000 CONTINUE
!
!     ==========================================================================
!     == REPORT JUMPS                                                         ==
!     ==========================================================================
!PRINT*,'JUMPS AFTER CLEAN'
      JUMP=>FIRSTJUMP
      DO 
!PRINT*,'JUMP%ON ',JUMP%ON
        IF(.NOT.JUMP%ON) GOTO 2000
        IF(JUMP%FORMING) THEN
          WRITE(NFIL,FMT='("BOND-FORMATION AT TIME:",F15.3,"PSEC",I8,F15.0)') &
     &                  JUMP%T/(PICO*SECOND),JUMP%ISTEP,JUMP%T
        ELSE
          WRITE(NFIL,FMT='("BOND-BREAKING AT TIME :",F15.3,"PSEC",I8,F15.0)') &
     &                  JUMP%T/(PICO*SECOND),JUMP%ISTEP,JUMP%T
        END IF
        ISTEP=JUMP%ISTEP
        IAT1=JUMP%BOND%IAT1
        IAT2=JUMP%BOND%IAT2
        CALL BONDS_NEIGHBORS(NFIL,IAT1,NB(ISTEP),BOND(:,ISTEP))
        CALL BONDS_NEIGHBORS(NFIL,IAT2,NB(ISTEP),BOND(:,ISTEP))
!
 2000   CONTINUE
        IF(.NOT.ASSOCIATED(JUMP%NEXT)) EXIT
        JUMP=>JUMP%NEXT
      ENDDO

!PRINT*,'MARKE 4'
 9999 CONTINUE
!
!     ==========================================================================
!     == WRITE NEIGHBORHOOD                                                   ==
!     ==========================================================================
      WRITE(NFIL,FMT='("NEIGHBOR TABLE"/14("="))')
      ISTEP=NSTEP
      DO IAT=1,NAT
        CALL BONDS_NEIGHBORS(NFIL,IAT,NB(ISTEP),BOND(:,ISTEP))
      ENDDO
      DEALLOCATE(BOND)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BONDS_NEIGHBORS(NFIL,IAT0,NB,BOND)
!     ******************************************************************
!     **   WRITES OUT A NEIGHBOR LIST FOR A SPECIFIED ATOM            **
!     ******************************************************************
      USE TRAJECTORY_MODULE, ONLY: BOND_TYPE,ATOM
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: NFIL
      INTEGER(4)     ,INTENT(IN) :: IAT0
      INTEGER(4)     ,INTENT(IN) :: NB
      TYPE(BOND_TYPE),INTENT(IN) :: BOND(NB)
      INTEGER(4)                 :: IPOINT(NB)
      INTEGER(4)                 :: IB,I,IAT1,J,J1
      INTEGER(4)                 :: NPOINT
      CHARACTER(32)  ,ALLOCATABLE:: NEIGHBOR(:)
      CHARACTER(32)              :: STRING
      INTEGER(4)                 :: IT(3)
!     ******************************************************************
      IPOINT(:)=0
!
!     ==========================================================================
!     ==  COUNT NUMBER OF NEIGHBORS                                           ==
!     ==========================================================================
      NPOINT=0
      DO IB=1,NB
        IF(BOND(IB)%IAT1.EQ.IAT0.OR.BOND(IB)%IAT2.EQ.IAT0) THEN
          NPOINT=NPOINT+1
          IPOINT(NPOINT)=IB
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  COUNT NUMBER OF NEIGHBORS                                           ==
!     ==========================================================================
      ALLOCATE(NEIGHBOR(NPOINT))
      DO I=1,NPOINT
        IB=IPOINT(I)
        IAT1=BOND(IB)%IAT2
        IT(:)=BOND(IB)%IT(:)
        IF(IAT1.EQ.IAT0) THEN
          IAT1=BOND(IB)%IAT1
          IT(:)=-IT(:)
        END IF
        IF(IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) THEN
          WRITE(STRING,*)IT
          J1=0
          DO J=1,32
            IF(STRING(J:J).EQ.' ') CYCLE
            J1=J1+1
            STRING(J1:J1)=STRING(J:J)
            STRING(J:J)=' '
          ENDDO
          STRING=TRIM(ATOM(IAT1))//':'//STRING
        ELSE
          STRING=ATOM(IAT1)
        END IF
        NEIGHBOR(I)=STRING
      ENDDO  
      WRITE(NFIL,FMT='(A7,": ",8A15)') &
     &     TRIM(ATOM(IAT0)),(NEIGHBOR(I),I=1,NPOINT)
      DEALLOCATE(NEIGHBOR)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BONDS_FIND(BOND0,NB,BOND,IB0)
!     ******************************************************************
!     **   CALCULATE TIME DERIVATIVE OF A COORDINATE TRAJECTORY       **
!     ******************************************************************
      USE TRAJECTORY_MODULE, ONLY: BOND_TYPE
      IMPLICIT NONE
      TYPE(BOND_TYPE),INTENT(IN) :: BOND0
      INTEGER(4)     ,INTENT(IN) :: NB
      TYPE(BOND_TYPE),INTENT(IN) :: BOND(NB)
      INTEGER(4)     ,INTENT(OUT):: IB0
      INTEGER(4)                 :: IAT1,IAT2
      INTEGER(4)                 :: IB
!     ******************************************************************
      IB0=0
      IAT1=BOND0%IAT1
      IAT2=BOND0%IAT2
      DO IB=1,NB
        IF(BOND(IB)%IAT1.EQ.IAT1) THEN
          IF(BOND(IB)%IAT2.NE.IAT2) CYCLE
          IF(BOND(IB)%IT(1).NE.BOND0%IT(1)) CYCLE
          IF(BOND(IB)%IT(2).NE.BOND0%IT(2)) CYCLE
          IF(BOND(IB)%IT(3).NE.BOND0%IT(3)) CYCLE
        ELSE IF(BOND(IB)%IAT1.EQ.IAT2) THEN
          IF(BOND(IB)%IAT2.NE.IAT1) CYCLE
          IF(BOND(IB)%IT(1).NE.-BOND0%IT(1)) CYCLE
          IF(BOND(IB)%IT(2).NE.-BOND0%IT(2)) CYCLE
          IF(BOND(IB)%IT(3).NE.-BOND0%IT(3)) CYCLE
        ELSE
          CYCLE
        END IF
        IB0=IB
        EXIT
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BONDS_ADDBOND(BOND0,RBAS,NAT,R,RCOV,NBX,NB,BOND)
!     ******************************************************************
!     **   CALCULATE TIME DERIVATIVE OF A COORDINATE TRAJECTORY       **
!     ****************************************************************** 
      USE TRAJECTORY_MODULE, ONLY: BOND_TYPE
      IMPLICIT NONE
      TYPE(BOND_TYPE),INTENT(IN)   :: BOND0
      REAL(8)        ,INTENT(IN)   :: RBAS(3,3)
      INTEGER(4)     ,INTENT(IN)   :: NAT
      REAL(8)        ,INTENT(IN)   :: R(3,NAT)
      REAL(8)        ,INTENT(IN)   :: RCOV(NAT)
      INTEGER(4)     ,INTENT(IN)   :: NBX
      INTEGER(4)     ,INTENT(INOUT):: NB
      TYPE(BOND_TYPE),INTENT(INOUT):: BOND(NBX)
      INTEGER(4)                   :: IAT1,IAT2
      INTEGER(4)                   :: IB
      REAL(8)                      :: DR(3)
!     ******************************************************************
      NB=NB+1
      IF(NB.GT.NBX) THEN
        CALL ERROR$MSG('#(BONDS) EXCEEDS MAXIMUM')
        CALL ERROR$STOP('BONDS_FIND')
      END IF
      BOND(NB)=BOND0

!     == RECALCULATE DISTANCE ========
      IAT1=BOND(NB)%IAT1
      IAT2=BOND(NB)%IAT2
      DR(:)=R(:,IAT2)-R(:,IAT1) &
     &     +RBAS(:,1)*DBLE(BOND(NB)%IT(1)) &
     &     +RBAS(:,2)*DBLE(BOND(NB)%IT(2)) &
     &     +RBAS(:,3)*DBLE(BOND(NB)%IT(3))
      BOND(NB)%DIS=SQRT(DOT_PRODUCT(DR,DR))/(RCOV(IAT1)+RCOV(IAT2))
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BONDS_DELETEBOND(BOND0,NB,BOND)
!     ******************************************************************
!     **   CALCULATE TIME DERIVATIVE OF A COORDINATE TRAJECTORY       **
!     ****************************************************************** 
      USE TRAJECTORY_MODULE, ONLY: BOND_TYPE
      IMPLICIT NONE
      TYPE(BOND_TYPE),INTENT(IN)   :: BOND0
      INTEGER(4)     ,INTENT(INOUT):: NB
      TYPE(BOND_TYPE),INTENT(INOUT):: BOND(NB)
      INTEGER(4)                   :: IB,IB0
!     ******************************************************************
      CALL BONDS_FIND(BOND0,NB,BOND,IB0)
      IF(IB0.EQ.0) THEN
        CALL ERROR$MSG('SPECIFIED BOND DOES NOT EXIST')
        CALL ERROR$STOP('BONDS_DELETEBOND')
      END IF
      DO IB=IB0,NB-1
        BOND(IB)=BOND(IB+1)
      ENDDO
      BOND(NB)=BOND_TYPE(0,0,(/0,0,0/),0.D0)
      NB=NB-1
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BONDS_EVAL(RBAS,NAT,R,RCOV,SCALEX,NBX,NB,BOND)
!     **************************************************************************
!     **  LIST OF NEAREST NEIGHBOR BONDS                                      **
!     **************************************************************************
      USE TRAJECTORY_MODULE, ONLY: BOND_TYPE
      IMPLICIT NONE 
      INTEGER(4)     ,INTENT(IN) :: NAT              ! #(ATOMS)
      REAL(8)        ,INTENT(IN) :: R(3,NAT)         ! ATOMIC POSITION
      REAL(8)        ,INTENT(IN) :: RBAS(3,3)        ! ATOMIC POSITION
      REAL(8)        ,INTENT(IN) :: RCOV(NAT)        ! COVALENT RADIUS
      REAL(8)        ,INTENT(IN) :: SCALEX           ! CUTOFF CRITERION      
      INTEGER(4)     ,INTENT(IN) :: NBX              ! MAX#(BONDS PER ATOM)
      INTEGER(4)     ,INTENT(OUT):: NB               ! #(BONDS PER ATOM)
      TYPE(BOND_TYPE),INTENT(OUT):: BOND(NBX)
      TYPE(BOND_TYPE)            :: SBOND
      REAL(8)                    :: DIS,DIS0,DIS2,DIS2X
      REAL(8)                    :: DR(3)
      REAL(8)                    :: DR0(3)
      REAL(8)                    :: DR1(3)
      REAL(8)                    :: DR2(3)
      INTEGER(4)                 :: IAT1,IAT2,IT1,IT2,IT3,IB1,IB2
!     **************************************************************************
      DO IB1=1,NBX
        BOND(IB1)=BOND_TYPE(0,0,(/0,0,0/),0.D0)
      ENDDO
!
!     ==========================================================================
!     == COLLECT POTENTIAL NEIGHBORS                                          ==
!     ==========================================================================
      NB=0
      DO IAT1=1,NAT
        DO IAT2=IAT1,NAT
          DIS0=RCOV(IAT1)+RCOV(IAT2)
          DIS2X=( SCALEX*DIS0 )**2
          DR0(:)=R(:,IAT2)-R(:,IAT1) 
          DO IT1=-1,1
            DR1(:)=DR0(:)+RBAS(:,1)*DBLE(IT1) 
            DO IT2=-1,1
              DR2(:)=DR1(:)+RBAS(:,2)*DBLE(IT2) 
              DO IT3=-1,1
                IF(IAT1.EQ.IAT2) THEN
                  IF(IT1.LT.0) THEN
                    CYCLE
                  ELSE IF(IT1.EQ.0) THEN
                    IF(IT2.LT.0) THEN
                      CYCLE
                    ELSE IF(IT2.EQ.0) THEN
                      IF(IT3.LE.0) CYCLE
                    END IF
                  END IF  
                END IF
                DR(:)=DR2(:)+RBAS(:,3)*DBLE(IT3) 
                DIS2=DOT_PRODUCT(DR,DR)
                IF(DIS2.GT.DIS2X) CYCLE
                NB=NB+1
                IF(NB.GT.NBX) THEN
                  CALL ERROR$MSG('TO MANY NEIGHBORS')
                  CALL ERROR$I4VAL('NB',NB)
                  CALL ERROR$I4VAL('NBX',NBX)
                  CALL ERROR$STOP('BONDS_EVAL')
                END IF
                DIS=SQRT(DIS2)
                BOND(NB)=BOND_TYPE(IAT1,IAT2,(/IT1,IT2,IT3/),DIS/DIS0)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO 
!
!     ==========================================================================
!     == ORDER ACCORDING TO RELATIVE DISTANCE                                 ==
!     ==========================================================================
      DO IB1=1,NB
        DO IB2=IB1+1,NB
          IF(BOND(IB1)%DIS.LT.BOND(IB2)%DIS) THEN
            SBOND=BOND(IB1)
            BOND(IB1)=BOND(IB2)      
            BOND(IB2)=SBOND
          END IF
        ENDDO
      ENDDO  
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SELECTION(LL_CNTL_,NATOM,ATOM,TSELECT)
!     **************************************************************************
!     **                                                                      **
!     **  READS THE STANDARDIZED BLOCK !SELECT FROM THE                       **
!     **  INPUT LINKEDLIST STRUCTURE                                          **
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      INTEGER(4)   ,INTENT(IN) :: NATOM
      CHARACTER(*) ,INTENT(IN) :: ATOM(NATOM)
      LOGICAL(4)   ,INTENT(OUT):: TSELECT(NATOM)
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: N,I
      CHARACTER(16),ALLOCATABLE:: ATOMNAMES(:)
      LOGICAL(4)   ,PARAMETER  :: TPR=.FALSE.
      INTEGER(4)               :: IAT
!     **************************************************************************
      LL_CNTL=LL_CNTL_
      TSELECT(:)=.TRUE.
      CALL LINKEDLIST$EXISTL(LL_CNTL,'SELECT',1,TCHK)
      IF(.NOT.TCHK) RETURN
!
!     ==========================================================================
!     ==  READ OUT NAMES OF SELECTED ATOMS                                    ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'SELECT')
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOMS',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!SELECT:ATOMS IS MANDATORY')
        CALL ERROR$STOP('SELECTION')
      END IF
      CALL LINKEDLIST$SIZE(LL_CNTL,'ATOMS',1,N)  
      ALLOCATE(ATOMNAMES(N))
      CALL LINKEDLIST$GET(LL_CNTL,'ATOMS',1,ATOMNAMES(:))
      CALL LINKEDLIST$SELECT(LL_CNTL,'..')
!
!     ==========================================================================
!     ==  EVALUATE TSELECT VECTOR ENCODING SELECTED ATOMS                     ==
!     ==========================================================================
      TSELECT(:)=.FALSE.
      DO I=1,N
        TCHK=.FALSE.
        DO IAT=1,NATOM
          IF(ATOMNAMES(I).EQ.ATOM(IAT)) THEN
            IF(TCHK) THEN
              CALL ERROR$MSG('ATOM '//TRIM(ATOMNAMES(I))//' SELECTED TWICE')
              CALL ERROR$STOP('SELECTION')
            END IF
            TSELECT(IAT)=.TRUE.
            TCHK=.TRUE.
          END IF
        ENDDO
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('ATOM '//TRIM(ATOMNAMES(I))//' NOT IN THE LIST')
          CALL ERROR$STOP('SELECTION')
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  PRINT FOR TEST                                                      ==
!     ==========================================================================
      IF(TPR) THEN
        DO IAT=1,NATOM
          IF(TSELECT(IAT)) THEN
            I=I+1
            WRITE(*,FMT='(A8,L18,A8)')ATOM(IAT),TSELECT(IAT)
          ELSE
            WRITE(*,FMT='(A8,L18,A8)')ATOM(IAT),TSELECT(IAT)
          END IF
        ENDDO
      END IF
      DEALLOCATE(ATOMNAMES)
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITECSSR(NFIL,OBJECTNAME,RBAS,NAT,NAME,R,Q,TCRYSTAL)
!     **************************************************************************
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTERFACE OPERATOR (.DYAD.)
        FUNCTION DYADISCHES_PRODUCT(R1,R2) RESULT(R3)
          REAL(8), INTENT(IN) :: R1(3)
          REAL(8), INTENT(IN) :: R2(3)
          REAL(8)             :: R3(3)
        END FUNCTION DYADISCHES_PRODUCT
      END INTERFACE 
      INTEGER(4)  ,INTENT(IN) :: NFIL
      LOGICAL(4)  ,INTENT(IN) :: TCRYSTAL
      CHARACTER(*),INTENT(IN) :: OBJECTNAME
      INTEGER(4)  ,INTENT(IN) :: NAT
!     REAL(8)     ,INTENT(IN) :: RBAS(3,3)
      REAL(8)                 :: RBAS(3,3)
      CHARACTER(*),INTENT(IN) :: NAME(NAT)
      REAL(8)                 :: R(3,NAT)
!     REAL(8)     ,INTENT(IN) :: R(3,NAT)
      REAL(8)     ,INTENT(IN) :: Q(NAT)
      REAL(8)                 :: PI
      INTEGER(4)              :: NEIGH(8,NAT)
      REAL(8)                 :: RBASINV(3,3) ! RBAS**(-1)
      REAL(8)                 :: A,B,C            ! LENGTH OF LATTICE VECTORS
      REAL(8)                 :: ALPHA,BETA,GAMMA ! ANGLES BETWEEN LATTICE VECT
      REAL(8)                 :: DET
      INTEGER(4)              :: IAT,I,J
      REAL(8)                 :: VEC(3),SVAR,RBASNEU(3,3)
      REAL(8)                 :: ANGSTROM
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      NEIGH(:,:)=0
      CALL CONSTANTS$GET('ANGSTROM',ANGSTROM)
!
!     ==========================================================================
!     == CONVERT DATA TO UNITS OF LATTICE VECTORS                             ==
!     ==========================================================================
      IF(TCRYSTAL) THEN
        RBASINV(1,:)=RBAS(:,2).DYAD.RBAS(:,3)
        RBASINV(2,:)=RBAS(:,3).DYAD.RBAS(:,1)
        RBASINV(3,:)=RBAS(:,1).DYAD.RBAS(:,2)
        DET=DOT_PRODUCT(RBAS(:,1),RBASINV(1,:))
        RBASINV=RBASINV/DET
        A=SQRT(DOT_PRODUCT(RBAS(:,1),RBAS(:,1)))
        B=SQRT(DOT_PRODUCT(RBAS(:,2),RBAS(:,2)))
        C=SQRT(DOT_PRODUCT(RBAS(:,3),RBAS(:,3)))
        GAMMA =ACOS(DOT_PRODUCT(RBAS(:,1),RBAS(:,2))/(A*B))
        ALPHA =ACOS(DOT_PRODUCT(RBAS(:,2),RBAS(:,3))/(B*C))
        BETA  =ACOS(DOT_PRODUCT(RBAS(:,3),RBAS(:,1))/(C*A))
        ALPHA=180.D0/PI*ALPHA
        BETA =180.D0/PI*BETA
        GAMMA=180.D0/PI*GAMMA
        RBASNEU(:,:)=0.D0
        RBASNEU(3,3)=C
        RBASNEU(3,2)=B*COS(ALPHA*PI/180.D0)
        RBASNEU(2,2)=SQRT(B**2-RBASNEU(3,2)**2)
        RBASNEU(3,1)=A*COS(BETA*PI/180.D0)
        RBASNEU(2,1)=(A*B*COS(GAMMA*PI/180.D0)-RBASNEU(3,1)*RBASNEU(3,2)) &
       &             /RBASNEU(2,2)
        RBASNEU(1,1)=SQRT(A**2-RBASNEU(2,1)**2-RBASNEU(3,1)**2)
      END IF
!
!     ==========================================================================
!     == WRITE CSSR FILE                                                      ==
!     ==========================================================================
      IF(TCRYSTAL) THEN
        WRITE(NFIL,FMT='(T39,3F8.3 &
     &   /T22,3F8.3,T50,"SPGR = 1 P 1",T72,"OPT = 1" &
     &   /I4,''   1 CREATED BY PAW    '' &
     &   /"     0 ",A4,": ",A4)') &
     &   A/ANGSTROM,B/ANGSTROM,C/ANGSTROM,ALPHA,BETA,GAMMA &
     &                ,NAT,OBJECTNAME(1:4),OBJECTNAME(1:4)
      ELSE
        WRITE(NFIL,FMT='(//I4,''   1 CREATED BY PAW    '' &
     &   /"     0 ",A4,": ",A4)')NAT,OBJECTNAME(1:4),OBJECTNAME(1:4)
      END IF
      DO IAT=1,NAT
        IF(TCRYSTAL) THEN
          VEC=R(:,IAT)
          VEC=MATMUL(RBASINV,VEC)+100.D0
          VEC=MOD(VEC,1.D0)
          VEC=MATMUL(RBASNEU,VEC)
          WRITE(NFIL,FMT='(I4,1X,A5,3F10.5,1X,8I4,F8.3)') &
     &           IAT,NAME(IAT),VEC(:)/ANGSTROM,NEIGH(:,IAT),Q(IAT)

        ELSE
           WRITE(NFIL,FMT='(I4,1X,A5,3F10.5,1X,8I4,F8.3)') &
      &           IAT,NAME(IAT),R(:,IAT)/ANGSTROM,NEIGH(:,IAT),Q(IAT)
         END IF
      ENDDO

      RETURN
      END
!     ...1.........2.........3.........4.........5.........6.........7.........8
      FUNCTION DYADISCHES_PRODUCT(R1,R2) RESULT(R3)
        REAL(8), INTENT(IN) :: R1(3)
        REAL(8), INTENT(IN) :: R2(3)
        REAL(8)             :: R3(3)
        R3(1)=R1(2)*R2(3)-R1(3)*R2(2)
        R3(2)=R1(3)*R2(1)-R1(1)*R2(3)
        R3(3)=R1(1)*R2(2)-R1(2)*R2(1)
      END FUNCTION DYADISCHES_PRODUCT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CLUSTERIT(LL_CNTL_)
!     **************************************************************************
!     **                                                                      **
!     **  REDUCES THE TRAJECTORY TO ONLY A CLUSTER OF SPECIFIED ATOMS         **
!     **  ATOMS CAN BE REPEATED OUTSIDE THE UNIT CELL                         **
!     **  THESE ATOMS WILL GET NAMES IN THE EXTENDED ATOM NOTATION            **
!     **                                                                      **
!     **************************************************************************
      USE TRAJECTORY_MODULE
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)  ,INTENT(IN)  :: LL_CNTL_
      TYPE(LL_TYPE)               :: LL_CNTL
      TYPE(TRA_TYPE)              :: X_TRA
      CHARACTER(32)  ,ALLOCATABLE :: X_ATOM(:)
      INTEGER(4)     ,ALLOCATABLE :: X_IZ(:)
      REAL(8)        ,ALLOCATABLE :: X_MASS(:)
      TYPE(MODE_TYPE),ALLOCATABLE :: X_MODE(:)
      INTEGER(4)                  :: X_NAT
      REAL(8)                     :: X_RBAS(3,3)
      INTEGER(4)                  :: I2,IAT,X_IAT,IT,J,ISTEP,ISEP
      LOGICAL(4)                  :: TCHK
      REAL(8)                     :: TVEC(3)
      CHARACTER(32)               :: STRING
      INTEGER(4)                  :: NFILO
!     **************************************************************************
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'TCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'CLUSTER',1,TCHK)
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_CNTL,'CLUSTER')
      CALL LINKEDLIST$EXISTD(LL_CNTL,'ATOMS',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!TCNTL!CLUSTER:ATOMS NOT FOUND')
        CALL ERROR$STOP('CLUSTERIT')
      END IF
      CALL LINKEDLIST$SIZE(LL_CNTL,'ATOMS',1,X_NAT)  
!
!     ==========================================================================
!     ==  ALLOCATE ARRAYS                                                     ==
!     ==========================================================================
      ALLOCATE(X_ATOM(X_NAT))
      ALLOCATE(X_IZ(X_NAT))
      ALLOCATE(X_TRA%R(3,X_NAT,NSTEP))
      ALLOCATE(X_TRA%Q(X_NAT,NSTEP))
!
!     ==========================================================================
!     ==  COLLECT NEW DATA                                                    ==
!     ==========================================================================
      CALL LINKEDLIST$GET(LL_CNTL,'ATOMS',1,X_ATOM(:))
!
!     ==========================================================================
!     ==  REPORT                                                              ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL REPORT$TITLE(NFILO,'CLUSTER SELECTION')
      WRITE(NFILO,FMT='(5A15)')X_ATOM
!
!     ==========================================================================
!     ==  MAPPING                                                             ==
!     ==========================================================================
      DO X_IAT=1,X_NAT   
        ISEP=INDEX(X_ATOM(X_IAT),':')
        TVEC(:)=0.D0
        IF(ISEP.EQ.0) THEN    ! NO TRANSLATION
          STRING=X_ATOM(X_IAT)
        ELSE                  !TRANSLATION PRESENT
          STRING=X_ATOM(X_IAT)(ISEP+1:)
!         __RESOLVE TRANLATION VECTOR
          I2=1
          IF(STRING(1:1).EQ.'-'.OR.STRING(1:1).EQ.'+') I2=2
          READ(STRING(1:I2),*)IT
          TVEC(:)=TVEC(:)+RBAS(:,1)*DBLE(IT)
          STRING=STRING(I2+1:)
          I2=1
          IF(STRING(1:1).EQ.'-'.OR.STRING(1:1).EQ.'+') I2=2
          READ(STRING(1:I2),*)IT
          TVEC(:)=TVEC(:)+RBAS(:,2)*DBLE(IT)
          STRING=STRING(I2+1:)
          I2=1
          IF(STRING(1:1).EQ.'-'.OR.STRING(1:1).EQ.'+') I2=2
          READ(STRING(1:I2),*)IT
          TVEC(:)=TVEC(:)+RBAS(:,3)*DBLE(IT)
!         __ NOW PICK ATOM NAME_________________________________________________
          STRING=X_ATOM(X_IAT)(1:ISEP-1)
        END IF
!       ==  LOOK FOR ORIGINAL ATOM INDEX =======================================
        IAT=0
        DO J=1,NAT
          IF(STRING.EQ.ATOM(J)) THEN
            IAT=J
          END IF
        ENDDO
        IF(IAT.EQ.0) THEN
          CALL ERROR$MSG('ATOM NOT FOUND')
          CALL ERROR$CHVAL('ATOM',X_ATOM(X_IAT))
          CALL ERROR$STOP('CLUSTERIT')
        END IF
!       ==  COPY DATA INTO NEW ARRAY ===========================================
        X_IZ(X_IAT)=IZ(IAT)
        DO ISTEP=1,NSTEP
          X_TRA%R(:,X_IAT,ISTEP)=TRA%R(:,IAT,ISTEP)+TVEC(:)
          X_TRA%Q(X_IAT,ISTEP)=TRA%Q(IAT,ISTEP)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  MAP ARRAYS ONTO ORIGINAL ARRAYS                                     ==
!     ==========================================================================
      DEALLOCATE(ATOM)
      DEALLOCATE(IZ)
      DEALLOCATE(TRA%R)
      DEALLOCATE(TRA%Q)
!
      NAT=X_NAT
      ALLOCATE(ATOM(NAT))
      ALLOCATE(IZ(NAT))
      ALLOCATE(TRA%R(3,NAT,NSTEP))
      ALLOCATE(TRA%Q(NAT,NSTEP))
!
      ATOM(:)=X_ATOM(:)
      IZ(:)=X_IZ(:)
      TRA%R(:,:,:)=X_TRA%R(:,:,:)
      TRA%Q(:,:)=X_TRA%Q(:,:)
!
      DEALLOCATE(X_ATOM)
      DEALLOCATE(X_IZ)
      DEALLOCATE(X_TRA%R)
      DEALLOCATE(X_TRA%Q)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RESOLVEEXTENDEDNAME(XNAME,NAME,IT)
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
        IF(ICH.GE.48.AND.ICH.LE.57) THEN ! IF "0,1,...,9"
          IND=IND+1
          IT(IND)=SGN*(ICH-48)
          SGN=+1
        ELSE IF(ICH.EQ.43) THEN   ! IF "+"
          SGN=+1
        ELSE IF(ICH.EQ.45) THEN   ! IF "-"
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
