!
!.......................................................................
MODULE STRINGS_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: STRINGS                                                    **
!**                                                                   **
!**  PURPOSE: DEFINE CASING OPERATORS FOR STRING MANIPULATIONS        **
!**                                                                   **
!**  functions:                                                       **
!**    '+'   uppercase a string                                        **
!**    '-'   lowercase a string                                       **
!**                                                                   **
!***********************************************************************
PUBLIC
INTERFACE OPERATOR (-)
  MODULE PROCEDURE LOWER_CASE
END INTERFACE
INTERFACE OPERATOR (+) 
  MODULE PROCEDURE UPPER_CASE
END INTERFACE
CONTAINS
!       ..................................................................
        FUNCTION LOWER_CASE(OLD) RESULT(NEW)
          implicit none
          CHARACTER(*), INTENT(IN):: OLD
          CHARACTER(LEN(OLD))     :: NEW
          INTEGER(4)              :: I,ISVAR
!         ***************************************************************
          NEW=OLD 
          DO I=1,LEN(TRIM(OLD))
            ISVAR=IACHAR(OLD(I:I))
            IF(ISVAR.GE.65.AND.ISVAR.LE.90) NEW(I:I)=ACHAR(ISVAR+32)
          ENDDO
        END FUNCTION LOWER_CASE
!       ..................................................................
        FUNCTION UPPER_CASE(OLD) RESULT(NEW)
          implicit none
          CHARACTER(*), INTENT(IN):: OLD
          CHARACTER(LEN(OLD))     :: NEW
          INTEGER(4)              :: I,ISVAR
!         ***************************************************************
          NEW=OLD 
          DO I=1,LEN(TRIM(OLD))
            ISVAR=IACHAR(OLD(I:I))
            IF(ISVAR.GE.97.AND.ISVAR.LE.122) NEW(I:I)=ACHAR(ISVAR-32)
          ENDDO
        END FUNCTION UPPER_CASE
END MODULE STRINGS_MODULE

