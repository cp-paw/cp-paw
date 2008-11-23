      PROGRAM MAIN
!     *************************************************************************
!     ** makes the names of all module files lowercase                       ** 
!     *************************************************************************
      IMPLICIT NONE
      INTEGER      ,PARAMETER :: LINELENG=512
      CHARACTER(LEN=LINELENG) :: LINE
      CHARACTER(4)            :: MODSTRING='.MOD'
      integer                 :: i1,i2,i
!     *************************************************************************
      call LOWERCASE(MODSTRING)
      DO
        LINE(:)=' '
        READ(*,FMT='(A)',END=1000)LINE
        I1=1
        DO
          I2=INDEX(LINE(I1:LINELENG),MODSTRING) 
          if(i2.eq.0) exit
          i=index(line(i1:i2),' ',back=.true.)
          if(i.ne.0) i1=i1+i
          i=index(line(i1:i2),'/',back=.true.)
          if(i.ne.0) i1=i1+i
          call lowercase(line(i1:i2-1)) 
          i1=i2+1
        ENDDO        
 
        WRITE(6,FMT='(A)')TRIM(LINE)
      ENDDO
1000  CONTINUE
      STOP
      END PROGRAM MAIN
!       ..................................................................
        subroutine  LOWERCASE(OLD)
          CHARACTER(*), INTENT(INout):: OLD
          CHARACTER(LEN(OLD))     :: NEW
          INTEGER(4)              :: I,ISVAR
          NEW=OLD 
          DO I=1,LEN(TRIM(NEW))
            ISVAR=IACHAR(NEW(I:I))
            IF(ISVAR.GE.65.AND.ISVAR.LE.90) NEW(I:I)=ACHAR(ISVAR+32)
          ENDDO
          old=new
          return
        END subroutine LOWERCASE
