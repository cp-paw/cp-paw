!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE  FILEHANDLER_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: FILEHANDLER                                                **
!**                                                                   **
!**  PURPOSE: TO OPEN, CLOSE AND ORGANIZE FILE-UNITS                  **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    FILEHANDLER$SETROOT(ROOTNAME)                                  **
!**    FILEHANDLER$UNIT(ID,STRING)                                    **
!**    FILEHANDLER$SETFILE(ID,EXT,FILENAME)                           **
!**    FILEHANDLER$SETSPECIFICATION(ID,SPEC,VALUE)                    **
!**    FILEHANDLER$CLOSE(ID)                                          **
!**    FILEHANDLER$CLOSEALL                                           **
!**    FILEHANDLER$PRINTFILEOFUNIT(UNIT)                              **
!**    FILEHANDLER$REPORT(NFIL,STRING)                                **
!**    FILEHANDLER$FILENAME(ID,NAME)                                  **
!**                                                                   **
!**  REMARKS                                                          **
!**   1) A FILE WITH NAME "OUT" IS PRECONNECTED AS FILE               **
!**      WITH ID "PROT"                                               **
!**   2) THE DEFAULT FOR POSITION IS 'APPEND' AND NOT 'ASIS'          **
!**   3) IF A FILE IS CLOSED WHILE NOT BEING POSITIONED               **
!**      EITHER AT THE BEGINNING OR THE END OF THE FILE               **
!**      IT WILL BE POSITIONED AT THE END WHEN OPENED AGAIN           **
!**   SIDE EFFECTS:                                                   **
!**   THE FILE-HANDLER ASSIGNS FILE UNITS ALMOST RANDOMLY.            **
!**   BUT ONLY WITH UNIT >1000                                        **
!**   IF FILES ARE CONNECTED NOT USING THE ONE NEEDS TO CHECK         **
!**   WHETHER THIS UNIT IS ALREADY CONNECTED                          **
!**                                                                   **
!***********************************************************************
IMPLICIT NONE                                                 
TYPE FILE_TYPE                                                
  CHARACTER(32) :: ID         ! INTERNAL IDENTIFIER
  CHARACTER(512):: NAME       ! FILENAME
  INTEGER(4)    :: UNIT       ! FORTRAN UNIT
  LOGICAL(4)    :: OPEN       ! A UNIT IS ASSIGNED TO THE FILE
  LOGICAL(4)    :: USED       ! A FILE HAS BEEN OPENED
  CHARACTER(8)  :: STATUS     ! UNKNOWN,OLD,NEW,SCRATCH,REPLACE
  CHARACTER(12) :: POSITION   ! REWIND,APPEND,ASIS
  CHARACTER(8)  :: PERMISSION ! R,W,RW,N, (N CHOOSES THE DEFAULT)
  LOGICAL(4)    :: FORMATTED  ! FORMATTED/UNFORMATTED SWITCH
END TYPE FILE_TYPE
TYPE (FILE_TYPE),ALLOCATABLE :: FILE(:)
INTEGER(4)                   :: NFILMAX=0
INTEGER(4)      ,PARAMETER   :: NFILMAXDEFAULT=20
INTEGER(4)      ,PARAMETER   :: NFILMAXINCREMENT=10
INTEGER(4)                   :: LASTFILE=0
CHARACTER(512)               :: BARENAME=' '
LOGICAL(4)                   :: TLITTLEENDIAN=.TRUE. !INTEL ARCHITECTURE USES LITTLE ENDIAN
                                                     !IBM USES BIG ENDIAN
!***********************************************************************
CONTAINS
!       ..................................................................     
        SUBROUTINE FILEHANDLER_CREATE
!       ****************************************************************
!       ****************************************************************
        USE STRINGS_MODULE
        IMPLICIT NONE
        INTEGER(4)    :: I
!       ****************************************************************
        IF(ALLOCATED(FILE)) RETURN
        NFILMAX=NFILMAXDEFAULT
        ALLOCATE(FILE(NFILMAX))
        DO I=1,NFILMAX
          CALL FILEHANDLER_DEFAULT(FILE(I))
        ENDDO
        FILE(1)%ID       ='PROT'
        FILE(1)%NAME     =-'OUT'
        FILE(1)%UNIT     =-1
        FILE(1)%OPEN     =.FALSE.
        FILE(1)%USED     =.FALSE.
        FILE(1)%STATUS   ='UNKNOWN'
        FILE(1)%POSITION ='APPEND'
        FILE(1)%PERMISSION ='W'
        FILE(1)%FORMATTED=.TRUE.
        LASTFILE=1
        END SUBROUTINE FILEHANDLER_CREATE
!       ..................................................................     
        SUBROUTINE FILEHANDLER_DEFAULT(FILE_)
!       ****************************************************************
!       ****************************************************************
        IMPLICIT NONE
        TYPE (FILE_TYPE),INTENT(OUT) :: FILE_
!       ****************************************************************
        FILE_%ID       ='NO ID'
        FILE_%NAME     =' '
        FILE_%UNIT     =-1
        FILE_%OPEN     =.FALSE.
        FILE_%USED     =.FALSE.
        FILE_%STATUS   ='UNKNOWN'
        FILE_%POSITION ='ASIS'
        FILE_%PERMISSION ='N'
        FILE_%FORMATTED=.FALSE.
        END SUBROUTINE FILEHANDLER_DEFAULT
!       ..................................................................     
        SUBROUTINE FILEHANDLER_ENLARGE
!       ****************************************************************
!       ****************************************************************
        TYPE(FILE_TYPE),ALLOCATABLE :: TEMPFILE(:)
        INTEGER(4)                  :: TEMPNFILMAX
        INTEGER(4)                  :: I
!       ****************************************************************
        IF(ALLOCATED(FILE)) THEN
          TEMPNFILMAX=LASTFILE
          ALLOCATE(TEMPFILE(TEMPNFILMAX))
          TEMPFILE(:)=FILE(:)
          DEALLOCATE(FILE)
        END IF
        NFILMAX=NFILMAX+NFILMAXINCREMENT
        ALLOCATE(FILE(NFILMAX))
        IF(ALLOCATED(TEMPFILE)) THEN
          FILE(1:LASTFILE)=TEMPFILE(1:LASTFILE)
          DEALLOCATE(TEMPFILE)
        END IF
        DO I=LASTFILE+1,NFILMAX
          CALL FILEHANDLER_DEFAULT(FILE(I))
        ENDDO
        END SUBROUTINE FILEHANDLER_ENLARGE
!
!       ................................................................
        SUBROUTINE FILEHANDLER_LOOKUP(ID_,IFIL_)
!       ****************************************************************
!       ****************************************************************
        CHARACTER(*),INTENT(IN)  :: ID_
        INTEGER(4)  ,INTENT(OUT) :: IFIL_
        INTEGER(4)               :: I
!       ****************************************************************
        IFIL_=0
        DO I=1,LASTFILE
          IF(TRIM(ID_).EQ.TRIM(FILE(I)%ID)) THEN
            IFIL_=I
            RETURN
          END IF
        ENDDO
        RETURN
        END SUBROUTINE FILEHANDLER_LOOKUP
END MODULE FILEHANDLER_MODULE
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FILEHANDLER$SETCH(ID,VAL)
!     ==================================================================
!     ==  SET TABLE MAXIMUM TABLE SIZE                                ==
!     ==================================================================
      USE FILEHANDLER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: ID
      CHARACTER(*),INTENT(IN):: VAL
!     ******************************************************************
      IF(ID.EQ.'ENDIAN') THEN
        IF(VAL.NE.'LITTLE'.AND.VAL.NE.'BIG') THEN
          CALL ERROR$MSG('VAL MUST BE  "LITTLE" OR "BIG"')
          CALL ERROR$CHVAL('VAL',TRIM(VAL))
          CALL ERROR$CHVAL('ID',TRIM(ID))
          CALL ERROR$STOP('FILEHANDLER$SETCH')
        END IF
        TLITTLEENDIAN=(VAL.EQ.'LITTLE')
      ELSE 
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',TRIM(ID))
        CALL ERROR$STOP('FILEHANDLER$SETCH')
      END IF
      RETURN
      END
!
!     .................................................................
      SUBROUTINE FILEHANDLER$SETROOT(STRING_)
!     ==================================================================
!     ==  SET TABLE MAXIMUM TABLE SIZE                                ==
!     ==================================================================
      USE FILEHANDLER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: STRING_
!     ******************************************************************
      BARENAME=STRING_
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FILEHANDLER$ATTACHED(ID_,TCHK)
!     **************************************************************************
!     **  CHECK IF A FILE IS ATTACHED FOR THE ENTRY WITH ID                   **
!     **************************************************************************
      USE FILEHANDLER_MODULE, ONLY : FILE,FILEHANDLER_CREATE,FILEHANDLER_LOOKUP
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(OUT):: TCHK
      INTEGER(4)              :: IFIL
!     **************************************************************************
      IF(.NOT.ALLOCATED(FILE))CALL FILEHANDLER_CREATE
      CALL FILEHANDLER_LOOKUP(ID_,IFIL)
      TCHK=(IFIL.NE.0) 
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      RECURSIVE SUBROUTINE FILEHANDLER$UNIT(ID_,UNIT_)
!     **************************************************************************
!     **  RETURN FILE UNIT FOR A GIVEN FILE (OPEN IF NECCESARY)               **
!     **************************************************************************
      USE FILEHANDLER_MODULE, ONLY : FILE,FILEHANDLER_CREATE,FILEHANDLER_LOOKUP
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(OUT):: UNIT_
      INTEGER(4)              :: IFIL
!     **************************************************************************
      IF(.NOT.ALLOCATED(FILE))CALL FILEHANDLER_CREATE
      CALL FILEHANDLER_LOOKUP(ID_,IFIL)
      IF(IFIL.EQ.0) THEN
        CALL ERROR$MSG('FILE IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',TRIM(ID_))
        CALL ERROR$STOP('FILEHANDLER$UNIT')
      END IF
      IF(.NOT.FILE(IFIL)%OPEN) THEN
        CALL FILEHANDLER_OPEN(FILE(IFIL))
      END IF
      UNIT_=FILE(IFIL)%UNIT
      RETURN
      END
!
!     .................................................................
      SUBROUTINE FILEHANDLER$SETFILE(ID_,EXT_,NAME_)
!     ==================================================================
!     ==  SET FILE BY EXTENSION                                       ==
!     ==================================================================
      USE FILEHANDLER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      CHARACTER(*),INTENT(IN) :: NAME_
      LOGICAL(4)  ,INTENT(IN) :: EXT_
      INTEGER(4)              :: IFIL
!     ******************************************************************
      IF(.NOT.ALLOCATED(FILE))CALL FILEHANDLER_CREATE
!
!     ==================================================================
!     ==  CHECK WHETHER FILE ALREADY EXISTS                           ==
!     ==================================================================
      CALL FILEHANDLER_LOOKUP(ID_,IFIL)
      IF(IFIL.EQ.0) THEN
        IF(LASTFILE.GE.NFILMAX) CALL FILEHANDLER_ENLARGE
        LASTFILE=LASTFILE+1
        FILE(LASTFILE)%ID=ID_
        IFIL=LASTFILE
      END IF
!
!     ==================================================================
!     ==  SET  FILENAME                                               ==
!     ==================================================================
      IF(EXT_) THEN
        FILE(IFIL)%NAME=TRIM(BARENAME)//NAME_
      ELSE
        FILE(IFIL)%NAME=NAME_
      END IF
!
!     ==================================================================
!     ==  CLOSE AND REOPEN FILE WITH SAME UNIT IF ALREADY OPENED      ==
!     ==================================================================
      IF(FILE(IFIL)%OPEN) THEN
        CALL FILEHANDLER_OPEN(FILE(IFIL))
      END IF
      RETURN      
      END
!
!     .................................................................
      SUBROUTINE FILEHANDLER$SETSPECIFICATION(ID_,SPEC_,VALUE_)
!     ==================================================================
!     ==  SET SPECIFICATIONS                                          ==
!     ==     FORM=FORMATTED / FORM=UNFORMATTED                        ==
!     ==     STATUS=UNKNOWN / STATUS=NEW / STATUS=OLD /STATUS=SCRATCH ==
!     ==================================================================
      USE FILEHANDLER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      CHARACTER(*),INTENT(IN) :: SPEC_
      CHARACTER(*),INTENT(IN) :: VALUE_
      INTEGER(4)              :: IFIL
!     ******************************************************************
      IF(.NOT.ALLOCATED(FILE))CALL FILEHANDLER_CREATE
      CALL FILEHANDLER_LOOKUP(ID_,IFIL)
      IF(IFIL.EQ.0) THEN
        CALL ERROR$MSG('IDENTIFIER NOT IN THE LIST')
        CALL ERROR$CHVAL('IDENTIFIER',ID_)
        CALL ERROR$STOP('FILEHANDLER$SETSPECIFICATION')
      END IF
!     ==================================================================
!     ==   STATUS=UNKNOWN,OLD,NEW,REPLACE                             ==
!     ==================================================================
      IF(SPEC_.EQ.'STATUS') THEN
        FILE(IFIL)%STATUS=VALUE_
        IF(VALUE_.EQ.'UNKNOWN') THEN
!         ==  CLEAR FILE IF IT EXISTS AND OVERWRITE =====================
        ELSE  IF(VALUE_.EQ.'OLD') THEN
!         ==  PROGRAM STOPS IF FILE DOES NOT EXIST ======================
        ELSE  IF(VALUE_.EQ.'NEW') THEN
!         ==  PROGRAM STOPS IF FILE DOES ALREADY EXIST ==================
        ELSE  IF(VALUE_.EQ.'SCRATCH') THEN
!         ==  LOOK UP MANUAL  ===========================================
        ELSE  IF(VALUE_.EQ.'REPLACE') THEN
!         ==  DELETE FILE IF ONE EXISTS =================================
!         ==  AND CREATE A NEW FILE =====================================
        ELSE
          CALL ERROR$MSG('VALUE FOR SPECIFICATION STATUS NOT RECOGNIZED')
          CALL ERROR$CHVAL('IDENTIFIER',VALUE_)
          CALL ERROR$CHVAL('FILE%ID',FILE(IFIL)%ID)
          CALL ERROR$CHVAL('FILE%NAME',FILE(IFIL)%NAME)
          CALL ERROR$STOP('FILEHANDLER$SETSPECIFICATION')
        END IF
!     ==================================================================
!     ==   POSITION=ASIS,REWIND,APPEND                                ==
!     ==================================================================
      ELSE IF(SPEC_.EQ.'POSITION') THEN
        FILE(IFIL)%POSITION=VALUE_
        IF(VALUE_.EQ.'ASIS') THEN
        ELSE  IF(VALUE_.EQ.'REWIND') THEN
        ELSE  IF(VALUE_.EQ.'APPEND') THEN
        ELSE  
          CALL ERROR$MSG('VALUE FOR SPECIFICATION POSITION NOT RECOGNIZED')
          CALL ERROR$CHVAL('IDENTIFIER',VALUE_)
          CALL ERROR$CHVAL('FILE%ID',FILE(IFIL)%ID)
          CALL ERROR$CHVAL('FILE%NAME',FILE(IFIL)%NAME)
          CALL ERROR$STOP('FILEHANDLER$SETSPECIFICATION')
        END IF
!   
!     ==================================================================
!     ==   ACTION='READ','WRITE','READWRITE'                          ==
!     ==================================================================
      ELSE IF(SPEC_.EQ.'ACTION') THEN
        IF(VALUE_.EQ.'READ') THEN
          FILE(IFIL)%PERMISSION='R'
        ELSE IF(VALUE_.EQ.'READWRITE') THEN
          FILE(IFIL)%PERMISSION='RW'
        ELSE IF(VALUE_.EQ.'WRITE') THEN
          FILE(IFIL)%PERMISSION='W'
        ELSE 
          CALL ERROR$MSG('VALUE FOR SPECIFICATION ACTION NOT RECOGNIZED')
          CALL ERROR$CHVAL('IDENTIFIER',VALUE_)
          CALL ERROR$CHVAL('FILE%ID',FILE(IFIL)%ID)
          CALL ERROR$CHVAL('FILE%NAME',FILE(IFIL)%NAME)
          CALL ERROR$STOP('FILEHANDLER$SETSPECIFICATION')
        END IF
!     ==================================================================
!     ==   FORM='FORMATTED'/'UNFORMATTED'                              ==
!     ==================================================================
      ELSE IF(SPEC_.EQ.'FORM') THEN
        IF(VALUE_.EQ.'FORMATTED') THEN
          FILE(IFIL)%FORMATTED=.TRUE.
        ELSE IF(VALUE_.EQ.'UNFORMATTED') THEN
          FILE(IFIL)%FORMATTED=.FALSE.
        ELSE
          CALL ERROR$MSG('VALUE FOR SPECIFICATION FORM NOT RECOGNIZED')
          CALL ERROR$CHVAL('IDENTIFIER',VALUE_)
          CALL ERROR$CHVAL('FILE%ID',FILE(IFIL)%ID)
          CALL ERROR$CHVAL('FILE%NAME',FILE(IFIL)%NAME)
          CALL ERROR$STOP('FILEHANDLER$SETSPECIFICATION')
        END IF
!     ==================================================================
!     ==   INCORRECT SPECIFICATION                                    ==
!     ==================================================================
      ELSE
        CALL ERROR$MSG('SPECIFICATION NOT RECOGNIZED')
        CALL ERROR$CHVAL('SPEC',SPEC_)
        CALL ERROR$CHVAL('FILE%ID',FILE(IFIL)%ID)
        CALL ERROR$CHVAL('FILE%NAME',FILE(IFIL)%NAME)
        CALL ERROR$STOP('FILEHANDLER$SETSPECIFICATION')
      END IF
      RETURN
      END
!
!     .................................................................
      SUBROUTINE FILEHANDLER$CLOSE(ID_)
!     ==================================================================
!     ==  CLOSE                                                       ==
!     ==================================================================
      USE FILEHANDLER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)                  :: IFIL
!     ******************************************************************
      IF(.NOT.ALLOCATED(FILE))CALL FILEHANDLER_CREATE
      CALL FILEHANDLER_LOOKUP(ID_,IFIL)
      IF(IFIL.EQ.0) THEN
        CALL ERROR$MSG('FILE IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('STRING=',ID_)
        CALL ERROR$STOP('FILEHANDLER$CLOSE')
      END IF
      CALL FILEHANDLER_CLOSE(FILE(IFIL))
      RETURN
      END
!
!     .................................................................
      SUBROUTINE FILEHANDLER$CLOSEALL
!     ==================================================================
!     ==  CLOSEALL                                                    ==
!     ==================================================================
      USE FILEHANDLER_MODULE
      IMPLICIT NONE
      INTEGER(4):: I
!     ******************************************************************
      IF(.NOT.ALLOCATED(FILE))CALL FILEHANDLER_CREATE
      DO I=1,LASTFILE
         CALL FILEHANDLER_CLOSE(FILE(I))
      ENDDO
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE FILEHANDLER$DELETE(ID_)
!     **************************************************************************
!     **  DELETE THE FILE                                                     **
!     **************************************************************************
      USE FILEHANDLER_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)              :: IFIL
      INTEGER(4)              :: I
      INTEGER(4)              :: NFIL
      INTEGER(4)              :: IOS
      CHARACTER(128)          :: IOMSG
      LOGICAL                 :: OD
      LOGICAL                 :: TCHK
!     **************************************************************************
      IF(.NOT.ALLOCATED(FILE))CALL FILEHANDLER_CREATE
!
!     ==========================================================================
!     == IDENTIFY ENTRY FOR THIS FILE-ID                                      ==
!     ==========================================================================
      CALL FILEHANDLER_LOOKUP(ID_,IFIL)
      IF(IFIL.EQ.0) THEN
        CALL ERROR$MSG('FILE IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$CHVAL('STRING=',ID_)
        CALL ERROR$STOP('FILEHANDLER$CLOSE')
      END IF
!
!     ==========================================================================
!     == CLOSE FILE IF IT IS OPEN                                             ==
!     ==========================================================================
      CALL FILEHANDLER_CLOSE(FILE(IFIL))
!
!     ==========================================================================
!     == FIND AVAILABLE FORTRAN UNIT                                          ==
!     ==========================================================================
      NFIL=-1
      DO I=1000,10000
        INQUIRE(UNIT=I,OPENED=OD,IOSTAT=IOS,IOMSG=IOMSG)
        IF(IOS.NE.0) THEN
          CALL ERROR$MSG('ERROR WHILE SCANNING FOR AVALIABLE FORTRAN FILE UNIT')
          CALL ERROR$MSG('ERROR INQUIRING ABOUT A FILE')
          CALL ERROR$I4VAL('UNIT',I)
          CALL ERROR$CHVAL('IO MESSAGE',TRIM(IOMSG))
          CALL ERROR$CHVAL('FILE ID',FILE(IFIL)%ID)
          CALL ERROR$CHVAL('FILENAME ',FILE(IFIL)%NAME)
          CALL ERROR$STOP('FILEHANDLER$DELETE')
        END IF
        IF(OD) CYCLE
        NFIL=I
        EXIT
      ENDDO
      IF(NFIL.EQ.-1) THEN
        CALL ERROR$MSG('NO FORTRAN FILE UNIT NUMBER AVAILABLE')
        CALL ERROR$CHVAL('FILE ID',FILE(IFIL)%ID)
        CALL ERROR$CHVAL('FILENAME ',FILE(IFIL)%NAME)
        CALL ERROR$STOP('FILEHANDLER$DELETE')
      END IF
!
!     ==========================================================================
!     == CLOSE FILE IF IT DOES NOT EXIST                                      ==
!     ==========================================================================
      INQUIRE(FILE=FILE(IFIL)%NAME,EXIST=TCHK)
      IF(.NOT.TCHK) RETURN
      OPEN(UNIT=NFIL,FILE=FILE(IFIL)%NAME,STATUS='OLD')
      CLOSE(UNIT=NFIL,STATUS='DELETE')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FILEHANDLER$PRINTFILEOFUNIT(UNIT_)
!     **************************************************************************
!     **  SAVE                                                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: UNIT_
      LOGICAL(4)            :: TCHK
      CHARACTER(512)        :: STRING
!     **************************************************************************
      INQUIRE(UNIT=UNIT_,OPENED=TCHK,NAME=STRING)
      IF(TCHK) THEN
        WRITE(*,FMT='("UNIT ",I5," IS CONNECTED TO FILE ",A)') &
     &         UNIT_,TRIM(STRING)
      ELSE
        WRITE(*,FMT='("UNIT ",I2," IS NOT CONNECTED TO ANY FILE ")') &
     &         UNIT_
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FILEHANDLER$REPORT(NFIL_,STRING_)
!     **************************************************************************
!     **   STRING CAN BE AN IDENTIFIER OR ,'ALL', OR 'USED'                   **
!     **************************************************************************
      USE FILEHANDLER_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL_
      CHARACTER(*),INTENT(IN) :: STRING_
      INTEGER(4)              :: IFIL
!     **************************************************************************
      IF(.NOT.ALLOCATED(FILE))CALL FILEHANDLER_CREATE
      WRITE(NFIL_,FMT='()')
      WRITE(NFIL_,FMT='("FILE REPORT"/"===========")')
      WRITE(NFIL_,FMT='("IDENTIFIER",T27,"FILENAME")')
      IF(STRING_.EQ.'USED') THEN
        DO IFIL=1,LASTFILE
          IF(FILE(IFIL)%USED) CALL MYPRINT
        ENDDO
      ELSE IF(STRING_.EQ.'ALL') THEN
        DO IFIL=1,LASTFILE
          CALL MYPRINT
        ENDDO
      ELSE
        CALL FILEHANDLER_LOOKUP(STRING_,IFIL)
        IF(IFIL.NE.0) CALL MYPRINT
      END IF
      RETURN
      CONTAINS
        SUBROUTINE MYPRINT
          WRITE(NFIL_,FMT='(A25," ",A)') &
     &                       FILE(IFIL)%ID,TRIM(FILE(IFIL)%NAME)
        RETURN
        END SUBROUTINE MYPRINT
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE FILEHANDLER$FILENAME(ID_,NAME_)
!     **************************************************************************
!     ** RETURNS THE NAME OF THE FILE CONNECTED TO A GIVEN file ID            **
!     **************************************************************************
      USE FILEHANDLER_MODULE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: ID_
      CHARACTER(*) ,INTENT(OUT):: NAME_
      INTEGER(4)               :: IFIL
!     **************************************************************************
      IF(.NOT.ALLOCATED(FILE))CALL FILEHANDLER_CREATE
      CALL FILEHANDLER_LOOKUP(ID_,IFIL)
      IF(IFIL.EQ.0) THEN
        PRINT*,'IDENTIFIER NOT IN THE LIST'
        PRINT*,'IDENTIFIER: ',ID_
        NAME_=' '
        RETURN
      END IF
      NAME_=FILE(IFIL)%NAME
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      RECURSIVE SUBROUTINE FILEHANDLER_OPEN(FILE_)
!     **************************************************************************
!     **  OPEN FILES UNLESS ALREADY OPEN                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      USE FILEHANDLER_MODULE
      IMPLICIT NONE
      TYPE (FILE_TYPE),INTENT(INOUT) :: FILE_
      INTEGER(4)                     :: IERR
      INTEGER(4)                     :: IOS
      LOGICAL(4)                     :: OD
      CHARACTER(9)                   :: ACTION
      CHARACTER(11)                  :: FORM
      CHARACTER(6)                   :: STDOUT
      CHARACTER(6)                   :: STDERR
      CHARACTER(5)                   :: STDIN
      CHARACTER(9)                   :: DEVNULL
      INTEGER(4)                     :: I
      LOGICAL(4)                     :: TOPENIBM ! IBM CHOICE BETWEEN 
                                                 ! LITTLE AND BIG ENDIAN
      CHARACTER(1)                   :: CONVERT  ! USED TO TEST FOR  
                                                 ! LITTLE AND BIG ENDIAN
      CHARACTER(128)                 :: IOMSG
!     **************************************************************************
      STDOUT =-'STDOUT'
      STDIN  =-'STDIN'
      STDERR =-'STDERR'
      DEVNULL=-'/DEV/NULL'
!
      IF(TRIM(FILE_%NAME).EQ.STDOUT) THEN
        FILE_%UNIT=6
        FILE_%USED=.TRUE.
        RETURN
      ELSE IF(TRIM(FILE_%NAME).EQ.STDIN) THEN
        FILE_%UNIT=5
        FILE_%USED=.TRUE.
        RETURN
      ELSE IF(TRIM(FILE_%NAME).EQ.STDERR) THEN
        FILE_%UNIT=0
        FILE_%USED=.TRUE.
        RETURN
      END IF
!
!     ==================================================================
!     ==  FIND UNUSED FORTRAN UNIT                                    ==
!     ==================================================================
      IERR=1
      IF(FILE_%OPEN) THEN
        IERR=2
        CLOSE(FILE_%UNIT,IOSTAT=IOS,IOMSG=IOMSG)
        IF(IOS.NE.0) THEN 
          CALL ERROR$MSG('ERROR CLOSING A FILE')
          CALL ERROR$CHVAL('IO MESSAGE',TRIM(IOMSG))
          CALL ERROR$I4VAL('INTERNAL ID TO IDENTIFY ERROR LOCATION',IERR)
          CALL ERROR$CHVAL('FILE ID',FILE_%ID)
          CALL ERROR$CHVAL('FILENAME ',FILE_%NAME)
          CALL ERROR$STOP('FILEHANDLER_OPEN')
        END IF
        FILE_%OPEN=.FALSE.
      ELSE
        IERR=3
        DO I=1000,10000
          IF(I.EQ.0) CYCLE
          IF(I.EQ.5) CYCLE
          IF(I.EQ.6) CYCLE
          INQUIRE(UNIT=I,OPENED=OD,IOSTAT=IOS,IOMSG=IOMSG)
          IF(IOS.NE.0) THEN
            CALL ERROR$MSG('ERROR INQUIRING ABOUT A FILE')
            CALL ERROR$CHVAL('IO MESSAGE',TRIM(IOMSG))
            CALL ERROR$I4VAL('INTERNAL ID TO IDENTIFY ERROR LOCATION',IERR)
            CALL ERROR$CHVAL('FILE ID',FILE_%ID)
            CALL ERROR$CHVAL('FILENAME ',FILE_%NAME)
            CALL ERROR$STOP('FILEHANDLER_OPEN')
          END IF
          IF(.NOT.OD) THEN
            FILE_%UNIT=I
            GOTO 200
          END IF
        ENDDO
        CALL ERROR$MSG('NO FORTRAN FILE UNIT NUMBER AVAILABLE')
        CALL ERROR$STOP('FILEHANDLER_OPEN')
 200    CONTINUE
      END IF
!  
!     ==================================================================
!     == CHECK WHETHER ANOTHER UNIT IS OPENED ON THE SAME FILE        ==
!     ==================================================================
      IF(TRIM(FILE_%NAME).EQ.DEVNULL) THEN
        OD=.FALSE.
      ELSE
        INQUIRE(FILE=FILE_%NAME,OPENED=OD) 
      END IF
      IF(OD) THEN
        CALL ERROR$MSG('ANOTHER FILE WITH THE SAME NAME')
        CALL ERROR$MSG('IS ALREADY OPENED')
        CALL ERROR$CHVAL('FILE ID',FILE_%ID)
        CALL ERROR$CHVAL('FILENAME',FILE_%NAME)
        CALL ERROR$STOP('FILEHANDLER_OPEN')
      END IF
!  
!     ==================================================================
!     == RESOLVE  PARAMETERS                                          ==
!     ==================================================================
      IF(FILE_%FORMATTED) THEN
        FORM='FORMATTED'
      ELSE
        FORM='UNFORMATTED'
      END IF
      IF(TRIM(FILE_%PERMISSION).EQ.'R') THEN
        ACTION='READ'
      ELSE IF(TRIM(FILE_%PERMISSION).EQ.'W') THEN
        ACTION='WRITE'
      ELSE IF(TRIM(FILE_%PERMISSION).EQ.'RW') THEN
        ACTION='READWRITE'
      ELSE IF(TRIM(FILE_%PERMISSION).EQ.'N') THEN
        ACTION=' '
      END IF
!  
!     =========================================================================
!     == OPEN FILE                                                           ==
!     =========================================================================
      IF(TRIM(FILE_%PERMISSION).EQ.'N') THEN
!       == PERMISSION 'N' CHOOSES DEFAULT PERMISSION OF (R,W,RW) ==============
!       == SPECIAL BECAUSE PARAMETER ACTION IS OMITTED  =======================
        IERR=4
        IF(FILE_%FORMATTED) THEN
!         == RECL IS EXPLICITELY SPECIFIED TO AVOID ARBITRARY BREAKINF OF =====
!         == LINES. RECL IS THE MAX NUMBER IF CHARACTERS PER LINE. ============
!         == SEE: HTTP://GCC.GNU.ORG/BUGZILLA/SHOW_BUG.CGI?ID=20257
          OPEN(UNIT=FILE_%UNIT,IOSTAT=IOS,IOMSG=IOMSG &
     &        ,FILE=FILE_%NAME &
     &        ,STATUS=FILE_%STATUS &
     &        ,FORM='FORMATTED' &
     &        ,POSITION=FILE_%POSITION &
     &        ,RECL=10000)    !MAX VALUE ACCEPTED BY GFORTRAN
        ELSE
!          == RECL NOT SPECIFIED FOR UNFORMATTED FILES. USES DEFAULT.         ==
!          == THE DEFAULT OR RECL IS 2^31 BYTES FOR DIRECT ACCESS FILES       ==
           OPEN(UNIT=FILE_%UNIT,IOSTAT=IOS,IOMSG=IOMSG &
     &        ,FILE=FILE_%NAME &
     &        ,STATUS=FILE_%STATUS &
     &        ,FORM='UNFORMATTED' &
     &        ,POSITION=FILE_%POSITION)
        END IF
        IF(IOS.NE.0) THEN
          CALL ERROR$MSG('ERROR OPENING A FILE')
          CALL ERROR$CHVAL('IO MESSAGE',TRIM(IOMSG))
          CALL ERROR$I4VAL('INTERNAL ID TO IDENTIFY ERROR LOCATION',IERR)
          CALL ERROR$CHVAL('FILE ID',FILE_%ID)
          CALL ERROR$CHVAL('FILENAME ',FILE_%NAME)
          CALL ERROR$STOP('FILEHANDLER_OPEN')
        END IF
      ELSE
!       == HERE FILES WITH SPECIFIED ACTION ARE OPENED ========================
        IERR=5
        IF(FILE_%FORMATTED) THEN
          OPEN(UNIT=FILE_%UNIT,IOSTAT=IOS,IOMSG=IOMSG &
     &        ,FILE=FILE_%NAME &
     &        ,STATUS=FILE_%STATUS &
     &        ,FORM='FORMATTED' &
     &        ,POSITION=FILE_%POSITION &
     &        ,ACTION=ACTION &
     &        ,RECL=1000)
          IF(IOS.NE.0) THEN
            CALL ERROR$MSG('ERROR OPENING A FILE')
            CALL ERROR$CHVAL('IO MESSAGE',TRIM(IOMSG))
            CALL ERROR$I4VAL('INTERNAL ID TO IDENTIFY ERROR LOCATION',IERR)
            CALL ERROR$CHVAL('FILE ID',FILE_%ID)
            CALL ERROR$CHVAL('FILENAME ',FILE_%NAME)
            CALL ERROR$STOP('FILEHANDLER_OPEN')
          END IF
        ELSE
#IF DEFINED(CPPVAR_ENDIANCHECK)
          CALL FILEHANDLER_READCONVERT(FILE_,FORM,CONVERT)
          IF(CONVERT.EQ.'U') TOPENIBM=(.NOT.TLITTLEENDIAN) !UNKNOWN
          IF(CONVERT.EQ.'B') TOPENIBM=.TRUE.
          IF(CONVERT.EQ.'L') TOPENIBM=.FALSE.
          IF((ACTION.EQ.'WRITE'.AND.FILE_%POSITION.EQ.'REWIND').AND.TLITTLEENDIAN) TOPENIBM=.FALSE.
          IF((ACTION.EQ.'WRITE'.AND.FILE_%POSITION.EQ.'REWIND').AND..NOT.TLITTLEENDIAN) TOPENIBM=.TRUE.
          IF(TOPENIBM) THEN 
            !FILE IS IBM COMPATIBLE AND HAS TO STAY SO
PRINT*,'FILEHANDLER: ATTENTION: FILE ',TRIM(FILE_%NAME),' IS OPENED IBM-COMPATIBLE'
            OPEN(UNIT=FILE_%UNIT,IOSTAT=IOS,IOMSG=IOMSG &
                  &    ,FILE=FILE_%NAME &
                  &    ,STATUS=FILE_%STATUS &
                  &    ,FORM=FORM &
                  &    ,POSITION=FILE_%POSITION &
                  &    ,ACTION=ACTION &
                  &    ,CONVERT='BIG_ENDIAN')            
          ELSE
            ! FILE IS INTEL COMPATIBLE
PRINT*,'FILEHANDLER: ATTENTION: FILE ',TRIM(FILE_%NAME),' IS OPENED INTEL-COMPATIBLE'
            OPEN(UNIT=FILE_%UNIT,IOSTAT=IOS,IOMSG=IOMSG &
                  &    ,FILE=FILE_%NAME &
                  &    ,STATUS=FILE_%STATUS &
                  &    ,FORM=FORM &
                  &    ,POSITION=FILE_%POSITION &
                  &    ,ACTION=ACTION &
                  &    ,CONVERT='LITTLE_ENDIAN')
          END IF
#ELSE
          OPEN(UNIT=FILE_%UNIT,IOSTAT=IOS,IOMSG=IOMSG &
     &        ,FILE=FILE_%NAME &
     &        ,STATUS=FILE_%STATUS &
     &        ,FORM='UNFORMATTED' &
     &        ,POSITION=FILE_%POSITION &
     &        ,ACTION=ACTION)
#ENDIF
          IF(IOS.NE.0) THEN
            CALL ERROR$MSG('ERROR OPENING A FILE')
            CALL ERROR$CHVAL('IO MESSAGE',TRIM(IOMSG))
            CALL ERROR$I4VAL('INTERNAL ID TO IDENTIFY ERROR LOCATION',IERR)
            CALL ERROR$CHVAL('FILE ID',FILE_%ID)
            CALL ERROR$CHVAL('FILENAME ',FILE_%NAME)
            CALL ERROR$STOP('FILEHANDLER_OPEN')
          END IF
        END IF
      END IF
      FILE_%OPEN=.TRUE.
      FILE_%USED=.TRUE.
!  
!     ==================================================================
!     == RESET PARAMETERS                                             ==
!     ==================================================================
!     PRINT*,'==========================================='
!     PRINT*,'FILE%ID       ',FILE_%ID
!     PRINT*,'FILE%STATUS   ',FILE_%STATUS
!     PRINT*,'FILE%POSITION ',FILE_%POSITION
!     PRINT*,'==========================================='
      IF(FILE_%STATUS.EQ.'NEW') FILE_%STATUS='OLD'
!     FILE_%POSITION='ASIS'
      FILE_%USED=.TRUE.
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE FILEHANDLER$IOSTATMESSAGE(IOS,X)
!     ******************************************************************
!     **                                                              **
!     **  RETURNS THE ERROR MESSAGES FOR A GIVEN IOSTAT VALUE "IOS"   **
!     **  (AIX XL FORTRAN COMPILER/6000)                              **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    THE MESSAGES ARE NOT A FORTRAN STANDARD AND MAY NOT BE    **
!     **    CORRECT ON A DIFFERENT PLATFORM                           **
!     **                                                              **
!     ******************************************************************
      INTEGER(4)  ,INTENT(IN) :: IOS
      CHARACTER(*),INTENT(OUT) :: X
!     ******************************************************************
      X='NO MESSAGE FOR THIS IOSTAT VALUE'
      IF(IOS.EQ.-1)X='END OF FILE'
      IF(IOS.EQ.-2)X='END OF INTERNAL FILE'
      IF(IOS.EQ. 6)X='FILE NOT FOUND AND STATUS=OLD SPECIFIED IN OPEN'
      IF(IOS.EQ. 7)X='INCORRECT FORMAT OF LIST-DIRECTED INPUT IN FILE'
      IF(IOS.EQ. 8) &
     &    X='INCORRECT FORMAT OF LIST-DIRECTED INPUT IN INTERNAL FILE'
      IF(IOS.EQ. 9)X='DATA ITEM TOO LONG FOR INTERNAL FILE'
      IF(IOS.EQ.10)X='READ ERROR ON DIRECT FILE'
      IF(IOS.EQ.11)X='WRITE ERROR ON DIRECT FILE'
      IF(IOS.EQ.12)X='READ ERROR ON SEQUENTIAL FILE'
      IF(IOS.EQ.13)X='WRITE ERROR ON SEQUENTIAL FILE'
      IF(IOS.EQ.14)X='ERROR OPENING A FILE'
      IF(IOS.EQ.15)X='PERMANENT I/O ERROR ENCOUNTERED ON A FILE'
      IF(IOS.EQ.16)X='RECORD NUMBER INVALID FOR DIRECT I/O'
      IF(IOS.EQ.17)X='I/O STATEMENT NOT ALLOWED FOR DIRECT FILE'
      IF(IOS.EQ.18)X='DIRECT I/O STATEMENT NON A FILE NOT OPEN'
      IF(IOS.EQ.19)X='UNFORMATTED I/O ON FORMATTED FILE'
      IF(IOS.EQ.20)X='FORMATTED I/O ON UNFORMATTED FILE'
      IF(IOS.EQ.21)X='SEQUENTIAL I/O ON DIRECT FILE'
      IF(IOS.EQ.22)X='DIRECT I/O ON SEQUENTIAL FILE'
      IF(IOS.EQ.23)X='FILE CONNECTED TO ANOTHER UNIT'
      IF(IOS.EQ.24)X='OPEN SPECIFIERS DO NOT MATCH FILE ATTRIBUTES'
      IF(IOS.EQ.25)X='RECL NOT GIVEN ON OPEN FOR DIRECT FILE'
      IF(IOS.EQ.26)X='NEGATIVE RECORD LENGTH IN OPEN'
      IF(IOS.EQ.27)X='OPEN ACCESS SPECIFIER NOT VALID'
      IF(IOS.EQ.28)X='OPEN FORMAT SPECIFIER NOT VALID'
      IF(IOS.EQ.31)X='OPEN FILE SPECIFIER NOT VALID'
      IF(IOS.EQ.35)X='RECURSIVE I/O OPERATION'
      IF(IOS.EQ.36)X='INVALID UNIT NUMBER'
      IF(IOS.EQ.38)X='REWIND ERROR'
      IF(IOS.EQ.39)X='ENDFILE ERROR'
      IF(IOS.EQ.40)X='BACKSPACE ERROR'
      IF(IOS.EQ.41)X='LOGICAL DATA NOT FOUND IN FILE'
      IF(IOS.EQ.42)X='LOGICAL DATA NOT FOUND IN INTERNAL FILE' 
      IF(IOS.EQ.43)X='CHARACTER DATA FOUND INSTEAD OF COMPLEX DATA'
      IF(IOS.EQ.44)X='CHARACTER DATA FOUND INSTEAD OF LOGICAL DATA'
      IF(IOS.EQ.45)X='CHARACTER DATA FOUND INSTEAD OF NUMERIC DATA'
      IF(IOS.EQ.46)X='COMPLEX DATA FOUND INSTEAD OF CHARACTER DATA'
      IF(IOS.EQ.47)X='COMPLEX DATA FOUND INSTEAD OF LOGICAL DATA'
      IF(IOS.EQ.48)X='COMPLEX DATA FOUND INSTEAD OF NUMERIC DATA'
      IF(IOS.EQ.49)X='NUMERIC DATA FOUND INSTEAD OF CHARACTER DATA'
      IF(IOS.EQ.50)X='NUMERIC DATA FOUND INSTEAD OF COMPLEX DATA'
      IF(IOS.EQ.56)X='INVALID DIGIT FOR BINARY, OCTAL OR HEXADECIMAL'
      IF(IOS.EQ.84)X='NAMELIST NAME NOT FOUND IN FILE'
      IF(IOS.EQ.85)X='NAMELIST NAME NOT FOUND IN INTERNAL FILE'
      IF(IOS.EQ.93)X='I/O STATEMENT NOT ALLOWED ON ERROR UNIT'
      IF(IOS.EQ.98)X='LENGTH OF DATA TOO LONG'
      IF(IOS.EQ.107)X='FILE EXISTS, STATUS=NEW SPECIFIED'
      IF(IOS.EQ.108) &
     &    X='THE I/O STATEMENT REFERS TO UNIT WHICH IS NOT CONNECTED'
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE IOSTATMESSAGE(IOS,NFIL)
!     ******************************************************************
!     **                                                              **
!     **  PRINTS THE ERROR MESSAGES FOR A GIVEN IOSTAT VALUE "IOS"    **
!     **  (AIX XL FORTRAN COMPILER/6000)                              **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    CALLS TO THIS ROUTINE SHOULD BE REPLACED BY A CALL TO     **
!     **    FILEHANDLER$IOSTATMESSAGE(IOS,X) AND AN EXPLICIT          **
!     **    PRINT STATEMENT OR CALL OF THE ERROR ROUTINE              **
!     **                                                              **
!     ******************************************************************
      INTEGER(4)    ,INTENT(IN):: IOS
      INTEGER(4)    ,INTENT(IN):: NFIL
      CHARACTER(128)           :: STRING
      CHARACTER(128)           :: X
!     ******************************************************************
      CALL FILEHANDLER$IOSTATMESSAGE(IOS,X)
      PRINT*,'ERROR WHILE READING OR WRITING ON FILE ',NFIL
      WRITE(*,FMT='('' AIX XL FORTRAN COMPILER MESSAGE FOR IOSTAT='',I3' &
     &     //'/A58)')IOS,TRIM(X)
      IF(NFIL.NE.0) THEN
        INQUIRE(NFIL,NAME=STRING)
        PRINT*,STRING
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE FILEHANDLER_CLOSE(FILE_)
!     ******************************************************************
!     ******************************************************************
      USE STRINGS_MODULE
      USE FILEHANDLER_MODULE, ONLY : FILE_TYPE
      IMPLICIT NONE
      TYPE (FILE_TYPE),INTENT(INOUT) :: FILE_
      INTEGER(4)                     :: NFIL
      CHARACTER(9)                   :: DEVNULL
!     ******************************************************************
      IF(.NOT.FILE_%OPEN) RETURN
      DEVNULL=-'/DEV/NULL'
      IF(FILE_%NAME.EQ.DEVNULL) RETURN  ! /DEV/NULL CANNOT BE CLOSED
      NFIL=FILE_%UNIT
      IF(NFIL.EQ.0.OR.NFIL.EQ.5.OR.NFIL.EQ.6) RETURN
!     INQUIRE(UNIT=NFIL,POSITION=FILE_%POSITION)
      IF(FILE_%POSITION.EQ.'ASIS')FILE_%POSITION='APPEND'
      CLOSE(NFIL)
      FILE_%OPEN=.FALSE.
      FILE_%UNIT=-1      
      RETURN
      END
!
#IF DEFINED(CPPVAR_ENDIANCHECK)
!     ..................................................................
      SUBROUTINE FILEHANDLER_READCONVERT(FILE_,FORM,CONVERT)
!     ******************************************************************
!     **  DETERMINES IF A FILE IS WRITTEN IN LITTLE OR BIG-ENDIAN.    **
!     **  TEST IF THE FIRST BYTE IS READABLE WITH A GIVEN SELECTION   **
!     **  IF READING IS NOT POSSIBLE (ERR=?), TRIES THE OTHER SELECTION**
!     ******************************************************************
      USE FILEHANDLER_MODULE, ONLY : FILE_TYPE
      IMPLICIT NONE
      TYPE (FILE_TYPE),INTENT(IN)  :: FILE_
      CHARACTER(11)   ,INTENT(IN)  :: FORM
      CHARACTER(1)    ,INTENT(OUT) :: CONVERT
      BYTE                         :: IVAR
      INTEGER(4)                   :: IOS
      CONVERT='L'
      IF(FORM.EQ.'FORMATTED') RETURN
      OPEN(UNIT=FILE_%UNIT,IOSTAT=IOS,ERR=9998 &
           &    ,FILE=FILE_%NAME &
           &    ,STATUS=FILE_%STATUS &
           &    ,FORM=FORM &
           &    ,POSITION='REWIND' &
           &    ,ACTION='READ' &
           &    ,CONVERT='LITTLE_ENDIAN')
      READ(FILE_%UNIT,ERR=200) IVAR ! A BYTE VARIABLE CAN BE READ 
                                    ! AT THE BEGINNING OF EVERY FILE
      CLOSE(FILE_%UNIT)
      RETURN ! READING HAPPEND WITHOUT ERROR, SO CONVERT=LITTLE=INTEL
200   CLOSE(FILE_%UNIT)
      OPEN(UNIT=FILE_%UNIT,IOSTAT=IOS,ERR=9998 &
           &    ,FILE=FILE_%NAME &
           &    ,STATUS=FILE_%STATUS &
           &    ,FORM=FORM &
           &    ,POSITION='REWIND' &
           &    ,ACTION='READ' &
           &    ,CONVERT='BIG_ENDIAN')
      READ(FILE_%UNIT,ERR=300) IVAR
      CONVERT='B'
      CLOSE(FILE_%UNIT)
      RETURN
300   PRINT*,'FILEHANDLER: ATTENTION: FILE COLD NOT BE READ IN ANY MODE: ',TRIM(FILE_%NAME)
      CONVERT='U'
      CLOSE(FILE_%UNIT)
9998  RETURN ! IN CASE OF OPENING ERROR GO BACK AND GIVE ERROR MESSAGE OF FILEHANDLER_OPEN
      END SUBROUTINE FILEHANDLER_READCONVERT
#ENDIF
