!************************************************************************
!***********************************************************************
!**                                                                   **
!**  NAME: PDOS                                                       **
!**                                                                   **
!**  PURPOSE: ANALYSIS TOOL FOR DENSITY OF STATES                     **
!**                                                                   **
!***********************************************************************
!***********************************************************************
MODULE SUBSTANCE_MODULE
TYPE SUBSTANCE_TYPE
  CHARACTER(32) :: NAME
  CHARACTER(256):: FILE
  REAL(8)       :: ETOT
END TYPE SUBSTANCE_TYPE
CONTAINS
!
!........................................................................
SUBROUTINE SUBSTANCE_ENERGY(NSUBSTANCE,SUBSTANCE,ID,ENERGY)
implicit none
INTEGER(4)          ,INTENT(IN) :: NSUBSTANCE
TYPE(SUBSTANCE_TYPE),INTENT(IN) :: SUBSTANCE(NSUBSTANCE)
CHARACTER(*)        ,INTENT(IN) :: ID
REAL(8)             ,INTENT(OUT):: ENERGY
integer(4)                      :: i
!**********************************************************************
DO I=1,NSUBSTANCE
 IF(ID.EQ.SUBSTANCE(I)%NAME) THEN
   ENERGY=SUBSTANCE(I)%ETOT
   RETURN
 END IF
ENDDO
CALL ERROR$MSG('SUBSTANCE NOT FOUND')
CALL ERROR$CHVAL('ID',ID)
CALL ERROR$STOP('SUBSTANCE_ENERGY')
STOP
END SUBROUTINE SUBSTANCE_ENERGY
END MODULE SUBSTANCE_MODULE
!
!     .................................................................
      PROGRAM GRAB
      USE LINKEDLIST_MODULE
      USE SUBSTANCE_MODULE
      IMPLICIT NONE
      LOGICAL(4)    ,PARAMETER  :: TPR=.FALSE.
      TYPE(LL_TYPE)             :: LL_CNTL
      LOGICAL(4)                :: TCHK,tfound
      INTEGER(4)                :: NFILO
      INTEGER(4)                :: NFIL
      INTEGER(4)                :: I,J
      INTEGER(4)                :: IPOS
      INTEGER(4)                :: NREACTION
      INTEGER(4)                :: NFROM
      INTEGER(4)                :: NTO
      INTEGER(4)                :: NSUBSTANCE
      TYPE(SUBSTANCE_TYPE),ALLOCATABLE :: SUBSTANCE(:)
      CHARACTER(256)            :: REACTIONSTRING
      CHARACTER(64)             :: REACTIONID
      CHARACTER(32)             :: ID
      REAL(8)                   :: REACTIONENERGY
      REAL(8)                   :: UNIT_ENERGY
      CHARACTER(16)             :: UNIT_ENERGY_NAME
      REAL(8)                   :: FAC
      REAL(8)                   :: ENERGY
      REAL(8)                   :: EREACT
      CHARACTER(16)             :: string
!     ******************************************************************
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT
!
!     ==================================================================
!     ==  RESOLVE ARGUMENTLIST AND INITIALIZE FILE HANDLER            ==
!     ==================================================================
      CALL INITIALIZEFILEHANDLER
!
!     ==================================================================
!     ==  READ CONTROL FILE                                           ==
!     ==================================================================
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL FILEHANDLER$UNIT('CNTL',NFIL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO) 
        CALL LINKEDLIST$REPORT(LL_CNTL,NFILO)
      END IF
!
!     ==================================================================
!     ==  WRITE HEADER                                                ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(80("*"),T15 &
     &             ,"           PROTOCOL GRAB TOOL                ")')
      WRITE(NFILO,FMT='(80("*"),T15 &
     &             ,"    FOR THE PROJECTOR-AUGMENTED WAVE METHOD  ")')
      WRITE(NFILO,FMT='(80("*"))')
      WRITE(NFILO,FMT='(T10 &
     &           ,"P.E. BLOECHL, CLAUSTHAL UNIVERSITY OF TECHNOLOGY ")')
      WRITE(NFILO,FMT='(T10 &
     &            ,"DISTRIBUTED UNDER THE GNU PUBLIC LICENSE V3")')
      WRITE(NFILO,*)
!
!     ==================================================================
!     ==  ANALYZE CONTROL FILE                                        ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'GENERIC',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
        CALL LINKEDLIST$EXISTD(LL_CNTL,'UNIT(E)',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL LINKEDLIST$SET(LL_CNTL,'UNIT(E)',1,'HARTREE')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'UNIT(E)',1,UNIT_ENERGY_NAME)
        IF(UNIT_ENERGY_NAME.NE.'HARTREE'.AND. &
     &     UNIT_ENERGY_NAME.NE.'JOULE'.AND. &
     &     UNIT_ENERGY_NAME.NE.'EV'.AND. &
     &     UNIT_ENERGY_NAME.NE.'KJ/MOL'.AND. &
     &     UNIT_ENERGY_NAME.NE.'KCAL/MOL'.AND. &
     &     UNIT_ENERGY_NAME.NE.'RY') THEN
           CALL ERROR$MSG('VALUE OF UNIT(E) NOT ALLOWED')
           CALL ERROR$CHVAL('UNIT(E)',UNIT_ENERGY_NAME)
           CALL ERROR$MSG('PERMITTED VALUES ARE:')
           CALL ERROR$MSG('"HARTREE","JOULE","EV","KJ/MOL","KCAL/MOL,RY"')
           CALL ERROR$STOP('GRAB')
        END IF
        CALL CONSTANTS(UNIT_ENERGY_NAME,UNIT_ENERGY)
      ELSE
        UNIT_ENERGY=1.D0
        UNIT_ENERGY_NAME='HARTREE'
      END IF
!
!     ==================================================================
!     ==  ANALYZE CONTROL FILE                                        ==
!     ==================================================================
      WRITE(NFILO,FMT='(/"SUBSTANCES"/"==========")')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'SUBSTANCE',NSUBSTANCE)
      ALLOCATE(SUBSTANCE(NSUBSTANCE))
      DO I=1,NSUBSTANCE
        CALL LINKEDLIST$SELECT(LL_CNTL,'SUBSTANCE',I)
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,SUBSTANCE(I)%NAME)
        tfound=.false.
!       == collect total energy from protocol file
        CALL LINKEDLIST$existd(LL_CNTL,'FILE',1,tchk)
        if(tchk) then
          CALL LINKEDLIST$GET(LL_CNTL,'FILE',1,SUBSTANCE(I)%FILE)
          SUBSTANCE(I)%ETOT=0.D0
          CALL GRABETOT(SUBSTANCE(I))
          tfound=.true.
        end if
!       == collect value directly
        CALL LINKEDLIST$EXISTD(LL_CNTL,'VALUE',1,TCHK)
        IF(TFOUND.AND.TCHK) THEN
          CALL ERROR$MSG('TWO VALUES MUST NOT BE GIVEN')
          CALL ERROR$MSG('CHECK SYNTAX OF !GCNTL!SUBSTANCE')
          CALL ERROR$STOP('GRAB')
        END IF
        IF(TCHK) THEN 
          TFOUND=.TRUE.
          CALL LINKEDLIST$GET(LL_CNTL,'VALUE',1,SUBSTANCE(I)%ETOT)
          CALL LINKEDLIST$EXISTD(LL_CNTL,'UNIT',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_CNTL,'UNIT',1,STRING)
            IF(STRING.NE.'HARTREE'.AND. &
     &        STRING.NE.'JOULE'.AND. &
     &        STRING.NE.'EV'.AND. &
     &        STRING.NE.'KJ/MOL'.AND. &
     &        STRING.NE.'KCAL/MOL'.AND. &
     &        STRING.NE.'RY') THEN
              CALL ERROR$MSG('VALUE OF UNIT(E) NOT ALLOWED')
              CALL ERROR$CHVAL('UNIT(E)',STRING)
              CALL ERROR$MSG('PERMITTED VALUES ARE:')
              CALL ERROR$MSG('"HARTREE","JOULE","EV","KJ/MOL","KCAL/MOL,RY"')
              CALL ERROR$STOP('GRAB')
            END IF
            CALL CONSTANTS(STRING,FAC)
            SUBSTANCE(I)%ETOT=SUBSTANCE(I)%ETOT*FAC
          END IF
        end if
        IF(.not.TFOUND) THEN
          CALL ERROR$MSG('No VALUES has been GIVEN')
          CALL ERROR$MSG('CHECK SYNTAX OF !GCNTL!SUBSTANCE')
          CALL ERROR$STOP('GRAB')
        END IF
        WRITE(NFILO,FMT='(60("."),T1,A,T60,": ",f20.10," ",A)') &
     &      TRIM(SUBSTANCE(I)%NAME),SUBSTANCE(I)%ETOT/UNIT_ENERGY,TRIM(UNIT_ENERGY_NAME)
        IF(TPR) THEN
          PRINT*,' NAME=  ',TRIM(SUBSTANCE(I)%NAME) &
     &        ,' FILE=  ',TRIM(SUBSTANCE(I)%FILE) &
     &        ,' ENERGY ',SUBSTANCE(I)%ETOT
        END IF
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
!
!     ==================================================================
!     ==  GRAB                                                        ==
!     ==================================================================
      WRITE(NFILO,FMT='(/"REACTIONS"/"=========")')
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'GCNTL')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'REACTION',NREACTION)
      DO I=1,NREACTION
        REACTIONSTRING=' '
        CALL LINKEDLIST$SELECT(LL_CNTL,'REACTION',I)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'ID',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL LINKEDLIST$SET(LL_CNTL,'ID',0,'NONE')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,reactionID)
        EREACT=0.D0
        IPOS=1
        CALL LINKEDLIST$NLISTS(LL_CNTL,'FROM',NFROM)
        DO J=1,NFROM
          CALL LINKEDLIST$SELECT(LL_CNTL,'FROM',J)
          CALL LINKEDLIST$GET(LL_CNTL,'ID',1,ID)
          CALL SUBSTANCE_ENERGY(NSUBSTANCE,SUBSTANCE,ID,ENERGY)
          CALL LINKEDLIST$GET(LL_CNTL,'FAC',1,FAC)
          EREACT=EREACT-FAC*ENERGY
          IF(ABS(FAC-1.D0).LT.1.D-2) THEN
            WRITE(REACTIONSTRING(IPOS:),FMT='(A)')TRIM(ID)//'+'
            IPOS=IPOS+LEN_TRIM(ID)+1
          ELSE
            WRITE(REACTIONSTRING(IPOS:),FMT='(F4.1,A)')FAC,'*'//TRIM(ID)//'+'
            IPOS=IPOS+4+1+LEN_TRIM(ID)+1
          END IF
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        IPOS=IPOS-1
        REACTIONSTRING(IPOS:)='->'
        IPOS=IPOS+2
        CALL LINKEDLIST$NLISTS(LL_CNTL,'TO',NTO)
        DO J=1,NTO
          CALL LINKEDLIST$SELECT(LL_CNTL,'TO',J)
          CALL LINKEDLIST$GET(LL_CNTL,'ID',1,ID)
          CALL SUBSTANCE_ENERGY(NSUBSTANCE,SUBSTANCE,ID,ENERGY)
          CALL LINKEDLIST$GET(LL_CNTL,'FAC',1,FAC)
          EREACT=EREACT+FAC*ENERGY!
          IF(ABS(FAC-1.D0).LT.1.D-2) THEN
            WRITE(REACTIONSTRING(IPOS:),FMT='(A)')TRIM(ID)//'+'
            IPOS=IPOS+LEN_TRIM(ID)+1
          ELSE
            WRITE(REACTIONSTRING(IPOS:),FMT='(F4.1,A)')FAC,'*'//TRIM(ID)//'+'
            IPOS=IPOS+4+1+LEN_TRIM(ID)+1
          END IF
          CALL LINKEDLIST$SELECT(LL_CNTL,'..')
        ENDDO
        IPOS=IPOS-1
        REACTIONSTRING(IPOS:IPOS)=' '
        IF(REACTIONID.NE.'NONE') THEN
          REACTIONSTRING=REACTIONID
        END IF
        IF(LEN_TRIM(REACTIONSTRING).LT.60) THEN
          IF(ABS(EREACT/UNIT_ENERGY).GT.999.999D0)THEN
            WRITE(NFILO,FMT='(60("."),T1,A,T60,"+ ",EN12.3," ",A)') &
     &       TRIM(REACTIONSTRING),EREACT/UNIT_ENERGY,TRIM(UNIT_ENERGY_NAME)
          ELSE
            WRITE(NFILO,FMT='(60("."),T1,A,T60,"+ ",F12.3,"     ",A)') &
     &       TRIM(REACTIONSTRING),EREACT/UNIT_ENERGY,TRIM(UNIT_ENERGY_NAME)
          END IF
        ELSE
          WRITE(NFILO,FMT='(A)')TRIM(REACTIONSTRING)
          WRITE(NFILO,FMT='(T60,"+ ",EN12.3," ",A)') &
     &       EREACT/UNIT_ENERGY,TRIM(UNIT_ENERGY_NAME)
        END IF
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
      WRITE(NFILO,*)
      CALL FILEHANDLER$REPORT(NFILO,'USED')
      WRITE(NFILO,FMT='(72("="))')
      WRITE(NFILO,FMT='(72("="),T20,"  PAW_GRAB TOOL FINISHED  ")')
      WRITE(NFILO,FMT='(72("="))')
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
      DEALLOCATE(SUBSTANCE)
      CALL FILEHANDLER$CLOSEALL
      CALL ERROR$NORMALSTOP()
      STOP
      END
!      
!     ...................................................................
      SUBROUTINE INITIALIZEFILEHANDLER
      use strings_module
      implicit none
      CHARACTER(256) :: ROOTNAME
      CHARACTER(256) :: PDOSINNAME
      INTEGER(4)     :: ISVAR
      INTEGER(4)     :: NARGS
!     *******************************************************************
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS.LT.1) THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$MSG('THE CONTROL FILE NAME MUST BE PROVIDED')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      END IF
      CALL GET_COMMAND_ARGUMENT(1,PDOSINNAME)
      ISVAR=INDEX(PDOSINNAME,-'.gCNTL',BACK=.TRUE.)
      IF(ISVAR.ne.0) THEN
        ROOTNAME=PDOSINNAME(1:ISVAR-1)
      ELSE
        ROOTNAME=' '
      END IF
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL STANDARDFILES
      CALL FILEHANDLER$SETFILE('CNTL',.FALSE.,PDOSINNAME)
      RETURN
      END
!
!     ..................................................................     
      SUBROUTINE STANDARDFILES
!     *****************************************************************
!     **                                                             **
!     *****************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(32)        :: CH32SVAR1
      CHARACTER(32)        :: ID
      INTEGER(4)           :: NFILO
!     *****************************************************************
                                   CALL TRACE$PUSH('STANDARDFILES')
!  
!     ==================================================================
!     == SET STANDARD FILENAMES                                       ==
!     ==================================================================
!
!     ==  ERROR FILE ===================================================
      ID=+'ERR'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.GERR')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','REPLACE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  PROTOCOL FILE ================================================
      ID=+'PROT'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.GPROT')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==  CONTROL FILE  == =============================================
      ID=+'CNTL'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,-'.GCNTL')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
                                   CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
MODULE READCNTLGRAB_MODULE
USE LINKEDLIST_MODULE
TYPE(LL_TYPE)   :: LL_CNTL
SAVE
END MODULE READCNTLGRAB_MODULE
!
!     ..................................................................
      SUBROUTINE READCNTL
      USE LINKEDLIST_MODULE
      USE READCNTLGRAB_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
      LOGICAL(4)           :: TCHK
      INTEGER(4)           :: NFIL
      CHARACTER(32)        :: ID 
      CHARACTER(256)       :: FILENAME
      INTEGER(4)           :: ITH
      INTEGER(4)           :: NUM
      INTEGER(4)           :: NFILO
!     ****************************************************************
                          CALL TRACE$PUSH('READCNTL')
!
!     ==================================================================
!     ==  READ CONTROL FILE                                           ==
!     ==================================================================
!
!     ==================================================================
!     ==  !PDOSIN!FILES!FILE                                          ==
!     ==================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'DCNTL')
      CALL LINKEDLIST$EXISTL(LL_CNTL,'FILES',1,TCHK)
      IF(.NOT.TCHK) RETURN
      CALL LINKEDLIST$SELECT(LL_CNTL,'FILES')
      CALL LINKEDLIST$NLISTS(LL_CNTL,'FILE',NUM)
      DO ITH=1,NUM
        CALL LINKEDLIST$SELECT(LL_CNTL,'FILE',ITH)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EXT',1,TCHK)
        IF(.NOT.TCHK)CALL LINKEDLIST$SET(LL_CNTL,'EXT',0,.FALSE.)
!       ==  READ ACTUAL VALUES  ======================================
        CALL LINKEDLIST$GET(LL_CNTL,'ID',1,ID)
        CALL LINKEDLIST$GET(LL_CNTL,'NAME',1,FILENAME)
        CALL LINKEDLIST$GET(LL_CNTL,'EXT',1,TCHK)
        CALL FILEHANDLER$SETFILE(ID,TCHK,FILENAME)
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
                          CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE GRABETOT(SUBSTANCE)
      USE STRINGS_MODULE
      USE SUBSTANCE_MODULE
      IMPLICIT NONE
      TYPE(SUBSTANCE_TYPE),INTENT(INOUT) :: SUBSTANCE
      INTEGER(4)                         :: NFIL
      CHARACTER(128)                     :: LINE
      CHARACTER(32)                      :: ID
      INTEGER(4)                         :: IPOS
!     ****************************************************************
                          CALL TRACE$PUSH('READCNTL')
!
!     ==================================================================
!     ==  SUBSTANCE FILE                                              ==
!     ==================================================================
      ID=+'SUBSTANCE'
      CALL FILEHANDLER$SETFILE(ID,.FALSE.,SUBSTANCE%FILE)
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('SUBSTANCE',NFIL)
      DO 
        READ(NFIL,FMT='(A)',END=9999)LINE
        IF(INDEX(LINE,+'TOTAL ENERGY',.FALSE.).NE.0) THEN
          IPOS=INDEX(LINE,':',.FALSE.)
          LINE=LINE(IPOS+1:)
          READ(LINE,*)SUBSTANCE%ETOT
        END IF
      ENDDO
 9999 CONTINUE
                          CALL TRACE$POP
      RETURN
      END
