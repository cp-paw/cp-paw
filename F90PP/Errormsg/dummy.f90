      MODULE STRINGS
      PUBLIC
      INTERFACE OPERATOR (+) 
        MODULE PROCEDURE UPPER_CASE
      END INTERFACE
      CONTAINS
!       ..................................................................
        FUNCTION UPPER_CASE(OLD) RESULT(new)
          CHARACTER(*), INTENT(IN):: OLD
          CHARACTER(LEN(OLD))     :: new
          new=OLD 
        END FUNCTION UPPER_CASE
      END MODULE STRINGS
!
!     ..................................................................
      SUBROUTINE UPPERCASE1(STRINGIN,STRING)
      USE STRINGS
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRINGIN
      CHARACTER(*),INTENT(OUT):: STRING
      STRING=+STRINGIN
      RETURN
      END
!     ..................................................................
      PROGRAM MAIN
      USE STRINGS
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: LINELENG=256
      INTEGER(4),PARAMETER :: NLINEX=100
      CHARACTER(LINELENG)  :: BUFFER(NLINEX)
      CHARACTER(LINELENG)  :: BUFFERdummy
      INTEGER(4)           :: I
!     ******************************************************************
      DO I=1,NLINEX 
        READ(*,FMT='(A)',END=1000)BUFFER(I)
print*,'test1a',trim(buffer(i))
!      call uppercase1(buffer(i),bufferdummy)
       bufferdummy=+buffer(i)
print*,'test1b'
      ENDDO
 1000 CONTINUE
print*,'done'
      STOP
      END
!Locate the critical line by looking for the string "+buffer" The +
!operator is defined in the module strings. It was intended as string
!operation to make text uppercase. Here it is stripped down to a simple
!copy operation. 
!
!The problem seems to be that the string length of the result is not
!recognized, because that is defined only in the +operation. I obtain
!an error message memory-access-error (in german "Speichezugriffsfehler")
!If I package the function into a subroutine, the error seems to disappear.
!Since that subroutine "uppercase1" still contains the same construction, 
!I suspect that it may still with different input. 
!
!The error is not systematic but occurs randomly. Therefore please use the 
!file dummy.in!
!
! Proceed as follows:
! f90 dummy.f90
! a.out <dummy.in
! If the message 'done' is not printed as last executable statement,
! the code terminated with error.


