!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE STRINGS_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: STRINGS                                                            **
!**                                                                           **
!**  PURPOSE: DEFINE CASING OPERATORS FOR STRING MANIPULATIONS                **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**    '+'   UPPERCASE A STRING                                               **
!**    '-'   LOWERCASE A STRING                                               **
!**                                                                           **
!*******************************************************************************
PUBLIC
INTERFACE OPERATOR (-)
  MODULE PROCEDURE LOWER_CASE
END INTERFACE
INTERFACE OPERATOR (+) 
  MODULE PROCEDURE UPPER_CASE
END INTERFACE
INTERFACE OPERATOR (.ITOS.)
  MODULE PROCEDURE I_2_S
END INTERFACE

CONTAINS
!       .1.........2.........3.........4.........5.........6.........7.........8
        FUNCTION LOWER_CASE(OLD) RESULT(NEW)
          IMPLICIT NONE
          CHARACTER(*), INTENT(IN):: OLD
          CHARACTER(LEN(OLD))     :: NEW
          INTEGER(4)              :: I,ISVAR
!         ***************************************************************
          NEW=OLD
          DO I=1,LEN(TRIM(OLD))
            ISVAR=IACHAR(OLD(I:I))
            IF(ISVAR.GE.65.AND.ISVAR.LE.90) NEW(I:I)=ACHAR(ISVAR+32)
          ENDDO
        return
        END FUNCTION LOWER_CASE
!       .1.........2.........3.........4.........5.........6.........7.........8
        FUNCTION UPPER_CASE(OLD) RESULT(NEW)
          IMPLICIT NONE
          CHARACTER(*), INTENT(IN):: OLD
          CHARACTER(LEN(OLD))     :: NEW
          INTEGER(4)              :: I,ISVAR
!         **********************************************************************
          NEW=OLD 
          DO I=1,LEN(TRIM(OLD))
            ISVAR=IACHAR(OLD(I:I))
            IF(ISVAR.GE.97.AND.ISVAR.LE.122) NEW(I:I)=ACHAR(ISVAR-32)
          ENDDO
        END FUNCTION UPPER_CASE
!
!       .1.........2.........3.........4.........5.........6.........7.........8
        FUNCTION I_2_S(I) RESULT(STRING)
        IMPLICIT NONE
        INTEGER(4),INTENT(IN)        :: I
        CHARACTER(256)               :: STRING
        INTEGER(4)                   :: IBASE
        INTEGER(4)    ,PARAMETER     :: IBASEMAX=10
        INTEGER(4)                   :: IFL
        REAL(8)                      :: RBJ
        REAL(8)                      :: IACT
        INTEGER(4)                   :: IWRITE
        LOGICAL(4)                   :: TZERO
!       ************************************************************************
        STRING=' '
        IACT=ABS(REAL(I,KIND=8))
        IF(IACT.GT.10**(IBASEMAX-1)) THEN
           CALL ERROR$MSG('ERROR STOP IN INTEGER TO STRING')
           CALL ERROR$MSG('YOUR INTEGER NUMBER IS LARGER THAN IMPLEMENTED')
           CALL ERROR$STOP('I_2_S')
        END IF
        IF(I.LT.0) THEN
           STRING(1:1)='-'
           IWRITE=1
        ELSE
          IWRITE=0
        END IF
        TZERO=.TRUE.
        DO IBASE=1,IBASEMAX 
           RBJ=10**(IBASEMAX-IBASE)
           IFL=FLOOR(IACT/RBJ)
           IACT=IACT-REAL(IFL,KIND=8)*RBJ
           TZERO=TZERO.AND.(IFL.EQ.0)   ! DROP LEADING ZEROS
           IF(TZERO) CYCLE
           IWRITE=IWRITE+1
           STRING(IWRITE:IWRITE)=ACHAR(48+IFL)
        END DO
      END FUNCTION I_2_S
!
!       .1.........2.........3.........4.........5.........6.........7.........8
        FUNCTION I_2_S_ALEX(I) RESULT(STRING)
        IMPLICIT NONE
        INTEGER(4),INTENT(IN)        :: I
        CHARACTER(256)               :: STRING
        INTEGER(4)                   :: IBASE
        REAL(8)                      :: RBASE=10
        INTEGER(4)                   :: IBASEMAX=10
        INTEGER(4)                   :: IFL
        REAL(8)                      :: RBJ
        REAL(8)                      :: IACT
        INTEGER(4)                   :: IWRITE
        LOGICAL(4)                   :: TZERO
!       ************************************************************************
        CALL ERROR$MSG('ROUTINE IS MARKED FOR DELETION')
        CALL ERROR$MSG('I_2_S_ALEX IS REPLACED BY I_2_S')
        CALL ERROR$STOP('I_2_S_ALEX')
        STRING=' '
        IWRITE=1
        IACT=ABS(REAL(I,KIND=8))
        TZERO=.TRUE.
        IF(IACT.GT.RBASE**(IBASEMAX-1)) THEN
           PRINT*, 'ERROR STOP IN INTEGER TO STRING'
           PRINT*, 'YOUR INTEGER NUMBER IS LARGER THAN IMPLEMENTED'
        END IF
        IF(I.LT.0) THEN
           STRING(1:1)='-'
           IWRITE=2
        END IF
        DO IBASE=1,IBASEMAX 
           RBJ=RBASE**(IBASEMAX-IBASE)
           IFL=FLOOR(IACT/RBJ)
           IACT=IACT-REAL(IFL,KIND=8)*RBJ
           SELECT CASE (IFL)
             CASE (0)
               IF(TZERO) THEN !WE DO NOT KEEP LEADING ZEROS
               IF(RBJ.EQ.1.D0) THEN !WE KEEP THE 0 FOR 10**0
                 STRING(IWRITE:IWRITE)='0'
                 IWRITE=IWRITE+1
               END IF
             ELSE  !WE HAVE A VALUE .NEQ. 0
               STRING(IWRITE:IWRITE)='0'
               IWRITE=IWRITE+1
             END IF
           CASE (1)
             STRING(IWRITE:IWRITE)='1'
             IWRITE=IWRITE+1
             TZERO=.FALSE.
          CASE (2)
            STRING(IWRITE:IWRITE)='2'
            IWRITE=IWRITE+1
            TZERO=.FALSE.
          CASE (3)
            STRING(IWRITE:IWRITE)='3'
            IWRITE=IWRITE+1
            TZERO=.FALSE.
          CASE (4)
            STRING(IWRITE:IWRITE)='4'
            IWRITE=IWRITE+1
            TZERO=.FALSE.
          CASE (5)
            STRING(IWRITE:IWRITE)='5'
            IWRITE=IWRITE+1
            TZERO=.FALSE.
          CASE (6)
            STRING(IWRITE:IWRITE)='6'
            IWRITE=IWRITE+1
            TZERO=.FALSE.
          CASE (7)
            STRING(IWRITE:IWRITE)='7'
            IWRITE=IWRITE+1
            TZERO=.FALSE.
          CASE (8)
            STRING(IWRITE:IWRITE)='8'
            IWRITE=IWRITE+1
            TZERO=.FALSE.
          CASE (9)
            STRING(IWRITE:IWRITE)='9'
            IWRITE=IWRITE+1
            TZERO=.FALSE.
          END SELECT
        END DO
      END FUNCTION I_2_S_ALEX

END MODULE STRINGS_MODULE

