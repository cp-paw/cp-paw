!=======================================================================
!=======================================================================
!=====                                                             =====
!=====   INTERFACE USED BY ANY OBJECT TO REPORT ITS SETTINGS.      =====
!=====   ASSURES A CONSISTENT OUTPUT FORMAT THAT CAN EASILY        =====
!=====   BE ADAPTED                                                =====
!=====                                                             =====
!=======================================================================
!=======================================================================
!=======================================================================
!
!     ..................................................................
      SUBROUTINE REPORT$TITLE(NFIL,TITLE)
!     ******************************************************************
!     **  WRITES A TITLE                                              **
!     ******************************************************************
      INTEGER(4)   ,INTENT(IN) :: NFIL
      CHARACTER(*) ,INTENT(IN) :: TITLE
      CHARACTER(80)            :: UNDERLINE
      INTEGER(4)               :: I
!     ******************************************************************
      UNDERLINE=' '
      DO I=1,LEN(TITLE)
        UNDERLINE(I:I)='='
      ENDDO
      WRITE(NFIL,*)
      WRITE(NFIL,FMT='(A)')TITLE
      WRITE(NFIL,FMT='(A)')TRIM(UNDERLINE)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE REPORT$R8VAL(NFIL,NAME,VALUE,UNIT)
!     ******************************************************************
!     **  WRITES A REAL(8) VALUE                                      **
!     ******************************************************************
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: NAME
      REAL(8)     ,INTENT(IN) :: VALUE
      CHARACTER(*),INTENT(IN) :: UNIT
!     ******************************************************************
      IF(ABS(VALUE).LT.1.D+3.AND.ABS(VALUE).GT.1.D-3) THEN
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,F10.5," ",A)')NAME,VALUE,UNIT
      ELSE
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,ES10.2," ",A)')NAME,VALUE,UNIT
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE REPORT$I4VAL(NFIL,NAME,VALUE,UNIT)
!     ******************************************************************
!     **  WRITES AM INTEGER(4) VALUE                                   **
!     ******************************************************************
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(IN) :: VALUE
      CHARACTER(*),INTENT(IN) :: UNIT
!     ******************************************************************
      WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,I10," ",A)')NAME,VALUE,UNIT
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE REPORT$L4VAL(NFIL,NAME,VALUE)
!     ******************************************************************
!     **  WRITES A LOGICAL VALUE                                      **
!     ******************************************************************
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: NAME
      LOGICAL(4)  ,INTENT(IN) :: VALUE
!     ******************************************************************
      IF(VALUE) THEN
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,"TRUE")')NAME
      ELSE
        WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,"FALSE")')NAME
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE REPORT$CHVAL(NFIL,NAME,VALUE)
!     ******************************************************************
!     **  WRITES A CHARACTER VALUE                                    **
!     ******************************************************************
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: NAME
      CHARACTER(*),INTENT(IN) :: VALUE
!     ******************************************************************
      WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,A)')NAME,TRIM(VALUE)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE REPORT$STRING(NFIL,STRING)
!     ******************************************************************
!     **  WRITES A STRING                                             **
!     ******************************************************************
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: STRING
!     ******************************************************************
      WRITE(NFIL,FMT='(A)')TRIM(STRING)
      RETURN
      END
!
!.......................................................................
MODULE RESTART_INTERFACE
TYPE SEPARATOR_TYPE
  SEQUENCE
  INTEGER(4)    :: NREC      ! #(RECORDS WRITTEN EXCLUDING SEPARATOR)
  CHARACTER(32) :: ID        ! IDENTIFIER SPECIFIES THE TYPE OF DATA
  CHARACTER(32) :: NAME      ! SOME SPECIAL NAME 
  CHARACTER(64) :: VERSION   ! VERSION ID
  CHARACTER(256):: NOTES     ! OTHER REMARKS
END TYPE SEPARATOR_TYPE
TYPE (SEPARATOR_TYPE) :: HEADER &
     =SEPARATOR_TYPE(0,'HEADER','NONE','AUG1996','DEFAULT HEADER')
TYPE (SEPARATOR_TYPE) :: ENDOFFILE &
     =SEPARATOR_TYPE(0,'ENDOFFILE','NONE','AUG1996','DEFAULT EOF')
CONTAINS
!     ..................................................................
      SUBROUTINE RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
!     ******************************************************************
!     **  READS SEPARATOR FROM FILE NFIL                              **
!     **                                                              **
!     **  DATA ARE READ ONLY IF ID AND NAME COINCIDE                  **
!     **                                                              **
!     **                                                              **
!     **  TCHK ON INPUT : DATASET SHOULD BE READ OR SKIPPED           **
!     **  TCHK ON OUTPUT: DATASET SHOULD BE READ OR SKIPPED           **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)           ,INTENT(IN)   :: NFIL
      INTEGER(4)           ,INTENT(IN)   :: NFILO
      LOGICAL(4)           ,INTENT(INOUT):: TCHK
      TYPE (SEPARATOR_TYPE),INTENT(INOUT):: SEPARATOR
      TYPE (SEPARATOR_TYPE)              :: XSEPARATOR
      LOGICAL(4)                         :: TCHK1
      INTEGER(4)                         :: I
!     ******************************************************************
      READ(NFIL)XSEPARATOR
      TCHK1=(XSEPARATOR%ID.NE.SEPARATOR%ID) 
      TCHK1=TCHK1.OR.(XSEPARATOR%NAME.NE.SEPARATOR%NAME)
      IF(TCHK1) THEN
!       == NO MATCH : REPOSITION =====================================
        BACKSPACE(NFIL)
        TCHK=.FALSE.
      ELSE
        IF(TCHK) THEN
!       == MATCH: DECIDE TO READ =====================================
          TCHK=.TRUE.
          SEPARATOR=XSEPARATOR
          IF(TRIM(SEPARATOR%NAME).EQ.'NONE') THEN
            WRITE(NFILO,FMT='("READING: ",A)')TRIM(SEPARATOR%ID)
          ELSE
            WRITE(NFILO,FMT='("READING: ",A,":",A)') &
                TRIM(SEPARATOR%ID),TRIM(SEPARATOR%NAME)
          END IF
        ELSE
!       == MATCH: DECIDE NOT TO READ =================================
          TCHK=.FALSE.
          IF(TRIM(SEPARATOR%NAME).EQ.'NONE') THEN
            WRITE(NFILO,FMT='("SKIPPING: ",A)')TRIM(SEPARATOR%ID)
          ELSE
            WRITE(NFILO,FMT='("SKIPPING: ",A,":",A)') &
                TRIM(SEPARATOR%ID),TRIM(SEPARATOR%NAME)
          END IF
          DO I=1,XSEPARATOR%NREC
            READ(NFIL)
          ENDDO
        END IF
      END IF
      RETURN
      END SUBROUTINE RESTART$READSEPARATOR
!
!     ..................................................................
      SUBROUTINE RESTART$WRITESEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
!     ******************************************************************
!     **  READS SEPARATOR FROM FILE NFIL                              **
!     **  TCHK IS NOT USED                                            **
!     ******************************************************************
      IMPLICIT NONE
      TYPE (SEPARATOR_TYPE),INTENT(IN)   :: SEPARATOR
      INTEGER(4)           ,INTENT(IN)   :: NFIL
      INTEGER(4)           ,INTENT(IN)   :: NFILO
      LOGICAL(4)           ,INTENT(IN)   :: TCHK
!     ******************************************************************
      WRITE(NFIL)SEPARATOR
      IF(TRIM(SEPARATOR%NAME).EQ.'NONE') THEN
        WRITE(NFILO,FMT='("WRITING : ",A)')TRIM(SEPARATOR%ID) 
      ELSE
        WRITE(NFILO,FMT='("WRITING : ",A,":",A)') &
                 TRIM(SEPARATOR%ID),TRIM(SEPARATOR%NAME)
      END IF
      RETURN
      END SUBROUTINE RESTART$WRITESEPARATOR
!     ..................................................................
      SUBROUTINE RESTART$SKIP(NFIL,NFILO)
!     ******************************************************************
!     **  READS SEPARATOR FROM FILE NFIL                              **
!     **  TCHK ON INPUT : DATASET SHOULD BE READ OR SKIPPED           **
!     **  TCHK ON OUTPUT: DATASET SHOULD BE READ OR SKIPPED           **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)           ,INTENT(IN)   :: NFIL
      INTEGER(4)           ,INTENT(IN)   :: NFILO
      TYPE (SEPARATOR_TYPE)              :: SEPARATOR
      INTEGER(4)                         :: I,NREC
!     ******************************************************************
      READ(NFIL)SEPARATOR
      WRITE(NFILO,FMT='("SKIPPING : ",A)')TRIM(SEPARATOR%ID)
      NREC=SEPARATOR%NREC
      IF(NREC.LT.0.OR.NREC.GT.10000000) THEN
        CALL ERROR$MSG('FILE STRUCTURE CORRUPTED')
        CALL ERROR$CHVAL('NAME',SEPARATOR%NAME)
        CALL ERROR$I4VAL('NREC',NREC)
        CALL ERROR$STOP('RESTART$SKIP')
      END IF
      DO I=1,NREC
        READ(NFIL)
      ENDDO
      RETURN
      END SUBROUTINE RESTART$SKIP
!
!     ..................................................................
      SUBROUTINE RESTART$CHECK(NFIL)
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)           ,INTENT(IN)   :: NFIL
      TYPE (SEPARATOR_TYPE)              :: SEPARATOR
      INTEGER(4)                         :: NREC,I
      INTEGER(4)                         :: ITEM
!     ******************************************************************
      ITEM=0
      REWIND NFIL
      DO 
        ITEM=ITEM+1
        READ(NFIL,END=1000)SEPARATOR
        WRITE(*,FMT='("CHECKING RESTART FILE: ",A10,I5)')SEPARATOR%ID &
     &                                                  ,SEPARATOR%NREC
        IF(ITEM.EQ.1.AND.SEPARATOR%ID.NE.HEADER%ID) THEN
          CALL ERROR$MSG('FILE STRUCTURE IS CORRUPTED')
          CALL ERROR$MSG('FILE DOES NOT BEGIN WITH THE HEADER')
          CALL ERROR$CHVAL('ID',SEPARATOR%ID)
          CALL ERROR$STOP('RESTART$CHECK')
        END IF
        IF(SEPARATOR%ID.EQ.ENDOFFILE%ID) RETURN
        NREC=SEPARATOR%NREC
        IF(NREC.LT.0.OR.NREC.GT.10000000) THEN
          CALL ERROR$MSG('FILE STRUCTURE CORRUPTED')
          CALL ERROR$CHVAL('NAME',SEPARATOR%ID)
          CALL ERROR$I4VAL('NREC',NREC)
          CALL ERROR$STOP('RESTART$CHECK')
        END IF
        DO I=1,NREC
          READ(NFIL,END=2000)
        ENDDO
      ENDDO
1000  CONTINUE
      CALL ERROR$MSG('FILE STRUCTURE IS CORRUPTED')
      CALL ERROR$MSG('FILE ENDS WITHOUT AN "ENDOFFILE"')
      CALL ERROR$STOP('RESTART$CHECK')
2000  CONTINUE
      CALL ERROR$MSG('FILE STRUCTURE IS CORRUPTED')
      CALL ERROR$MSG('FILE ENDS BEFORE ALL RECORDS ARE FINISHED')
      CALL ERROR$STOP('RESTART$CHECK')
      RETURN
    END SUBROUTINE RESTART$CHECK
END MODULE RESTART_INTERFACE
