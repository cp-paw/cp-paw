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
!     ******************************************************************
      UNDERLINE=' '
      DO I=1,LEN(TITLE)
        UNDERLINE(I:I)='='
      ENDDO
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
      IF(DABS(VALUE).LT.1.D+3.AND.DABS(VALUE).GT.1.D-3) THEN
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
      WRITE(NFIL,FMT='(55("."),": ",T1,A,T58,A10," ",A)')NAME,TRIM(VALUE)
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
  INTEGER(4)    :: NREC      ! #(records written excluding separator)
  CHARACTER(32) :: ID        ! identifier specifies the type of data
  CHARACTER(32) :: NAME      ! some special name 
  CHARACTER(64) :: VERSION   ! version id
  CHARACTER(256):: NOTES     ! other remarks
END TYPE SEPARATOR_TYPE
TYPE (SEPARATOR_TYPE) :: HEADER &
     =SEPARATOR_TYPE(0,'HEADER','NONE','AUG1996','DEFAULT HEADER')
TYPE (SEPARATOR_TYPE) :: ENDOFFILE &
     =SEPARATOR_TYPE(0,'ENDOFFILE','NONE','AUG1996','DEFAULT EOF')
CONTAINS
!     ..................................................................
      SUBROUTINE READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
!     ******************************************************************
!     **  READS SEPARATOR FROM FILE NFIL                              **
!     **                                                              **
!     **  data are read only if id and name coincide                  **
!     **                                                              **
!     **                                                              **
!     **  TCHK ON INPUT : DATASET SHOULD BE READ OR SKIPPED           **
!     **  TCHK ON OUTPUT: DATASET SHOULD BE READ OR SKIPPED           **
!     ******************************************************************
      USE MPE_MODULE
      implicit none
      INTEGER(4)           ,INTENT(IN)   :: NFIL
      INTEGER(4)           ,INTENT(IN)   :: NFILO
      LOGICAL(4)           ,INTENT(INOUT):: TCHK
      TYPE (SEPARATOR_TYPE),INTENT(INOUT):: SEPARATOR
      TYPE (SEPARATOR_TYPE)              :: XSEPARATOR
      INTEGER(4)                         :: NTASKS,THISTASK
      LOGICAL(4)                         :: TCHK1
      integer(4)                         :: i
!     ******************************************************************
!
      CALL MPE$QUERY(NTASKS,THISTASK)
!
      IF(THISTASK.EQ.1) THEN
        READ(NFIL)XSEPARATOR
        TCHK1=(XSEPARATOR%ID.NE.SEPARATOR%ID) 
        TCHK1=TCHK1.OR.(XSEPARATOR%NAME.NE.SEPARATOR%NAME)
        IF(TCHK1) THEN
!         == no match : reposition =====================================
          BACKSPACE(NFIL)
          TCHK=.FALSE.
        ELSE
          IF(TCHK) THEN
!         == match: decide to read =====================================
            TCHK=.TRUE.
            SEPARATOR=XSEPARATOR
            IF(TRIM(SEPARATOR%NAME).EQ.'NONE') THEN
              WRITE(NFILO,FMT='("READING: ",A)')TRIM(SEPARATOR%ID)
            ELSE
              WRITE(NFILO,FMT='("READING: ",A,":",A)') &
                  TRIM(SEPARATOR%ID),TRIM(SEPARATOR%NAME)
            END IF
          ELSE
!         == match: decide not to read =================================
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
      END IF
!
      CALL MPE$BROADCAST(1,TCHK)
!
      RETURN
      END SUBROUTINE READSEPARATOR
!
!     ..................................................................
      SUBROUTINE WRITESEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
!     ******************************************************************
!     **  READS SEPARATOR FROM FILE NFIL                              **
!     **  TCHK ON INPUT : DATASET SHOULD BE READ OR SKIPPED           **
!     **  TCHK ON OUTPUT: DATASET SHOULD BE READ OR SKIPPED           **
!     ******************************************************************
      implicit none
      TYPE (SEPARATOR_TYPE),INTENT(IN)   :: SEPARATOR
      INTEGER(4)           ,INTENT(IN)   :: NFIL
      INTEGER(4)           ,INTENT(IN)   :: NFILO
      LOGICAL(4)           ,INTENT(INOUT):: TCHK
      INTEGER(4)                         :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
      IF(THISTASK.EQ.1) THEN
        WRITE(NFIL)SEPARATOR
        IF(TRIM(SEPARATOR%NAME).EQ.'NONE') THEN
          WRITE(NFILO,FMT='("WRITING : ",A)')TRIM(SEPARATOR%ID) 
        ELSE
          WRITE(NFILO,FMT='("WRITING : ",A,":",A)') &
                   TRIM(SEPARATOR%ID),TRIM(SEPARATOR%NAME)
        END IF
      END IF
      RETURN
      END SUBROUTINE WRITESEPARATOR
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
      INTEGER(4)                         :: NTASKS,THISTASK
      INTEGER(4)                         :: I,NREC
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
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
END MODULE RESTART_INTERFACE
