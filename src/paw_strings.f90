!
!.......................................................................
MODULE STRINGS_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: STRINGS                                                    **
!**                                                                   **
!**  PURPOSE: DEFINE CASING OPERATORS FOR STRING MANIPULATIONS        **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    '+'   UPPERCASE A STRING                                        **
!**    '-'   LOWERCASE A STRING                                       **
!**                                                                   **
!***********************************************************************
PUBLIC
INTERFACE OPERATOR (-)
  MODULE PROCEDURE LOWER_CASE
END INTERFACE
INTERFACE OPERATOR (+) 
  MODULE PROCEDURE UPPER_CASE
END INTERFACE
INTERFACE OPERATOR (.itos.)
  MODULE PROCEDURE i_2_s
END INTERFACE

CONTAINS
!       ..................................................................
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
        END FUNCTION LOWER_CASE
!       ..................................................................
        FUNCTION UPPER_CASE(OLD) RESULT(NEW)
          IMPLICIT NONE
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
  function i_2_s(i) result(string)
    implicit none
    integer(4),intent(in)        :: i
    character(256)               :: string
    integer(4)                   :: ibase
    real(8)                      :: rbase=10
    integer(4)                   :: ibasemax=10
    integer(4)                   :: j,ifl
    real(8)                      :: rbj
    real(8)                      :: iact
    integer(4)                   :: iwrite
    logical(4)                   :: tzero
    
    string=''
    iwrite=1
    iact=abs(real(i,kind=8))
    tzero=.true.

    !checks
    if(i.gt.rbase**(ibasemax-1)) then
       print*, 'ERROR STOP IN INTEGER TO STRING'
       print*, 'YOUR INTEGER NUMBER IS LARGER THAN IMPLEMENTED'
    end if
    
    if(i.lt.0) then
       string(1:1)='-'
       iwrite=2
    end if
    
    do ibase=1,ibasemax 
       rbj=rbase**(ibasemax-ibase)
       ifl=floor(iact/rbj)
       iact=iact-real(ifl,kind=8)*rbj
       
       select case (ifl)
       case (0)
          if(tzero) then !we do not keep leading zeros
             if(rbj.eq.1.d0) then !we keep the 0 for 10**0
                string(iwrite:iwrite)='0'
                iwrite=iwrite+1
             end if
          else  !we have a value .neq. 0
             string(iwrite:iwrite)='0'
             iwrite=iwrite+1
          end if
       case (1)
          string(iwrite:iwrite)='1'
          iwrite=iwrite+1
          tzero=.false.
       case (2)
          string(iwrite:iwrite)='2'
          iwrite=iwrite+1
          tzero=.false.
       case (3)
          string(iwrite:iwrite)='3'
          iwrite=iwrite+1
          tzero=.false.
       case (4)
          string(iwrite:iwrite)='4'
          iwrite=iwrite+1
          tzero=.false.
       case (5)
          string(iwrite:iwrite)='5'
          iwrite=iwrite+1
          tzero=.false.
       case (6)
          string(iwrite:iwrite)='6'
          iwrite=iwrite+1
          tzero=.false.
       case (7)
          string(iwrite:iwrite)='7'
          iwrite=iwrite+1
          tzero=.false.
       case (8)
          string(iwrite:iwrite)='8'
          iwrite=iwrite+1
          tzero=.false.
       case (9)
          string(iwrite:iwrite)='9'
          iwrite=iwrite+1
          tzero=.false.
       end select
    end do
  end function i_2_s

END MODULE STRINGS_MODULE

