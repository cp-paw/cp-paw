!
!......................................................................
MODULE CLOCK_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: CLOCK                                                      **
!**                                                                   **
!**  PURPOSE: ALLOWS TO TAKE, COMPARE AND PRINT TIMES.                **
!**    THE TIME IS TYPICALLY TAKEN IN THE FORM OF A DATE_TIME TYPE    **
!**    (WHICH IS SUPPLIED AS WELL). IT CAN BE COMPARED AND CONVERTED  **
!**    INTO A PRINTABLE FORMAT.                                       **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    TIME1.LATER.TIME2   COMPARES TO TIMES                          **
!**    CLOCK$NOW       RETURNS THE CURRENT TIME IN ONE OF THE THREE   **
!**                    FORMATS                                        **
!**    CLOCK$TIMESTAMP(DATE_TIME,STAMP) PRODUCES A PRINTABLE TIMESTAMP**
!**                    FROM THE  TIME SUPPLIED IN DATE_TIME_TYPE      **
!**                                                                   **
!**  REMARKS: IN FUTURE THERE SHOULD BE THE POSSIBILITY OF            **
!**    INTERCONVERSION INTO ANY OF THE THREE REPRESENTATIONS.         **
!**    AND HANDLING INTERNALLY FLOATING POINT NUMBERS                 **
!**    INSTEAD OF DATE_TIME                                           **
!**                                                                   **
!***********************************************************************
TYPE DATE_TIME
  INTEGER(4)     :: YEAR
  INTEGER(4)     :: MONTH
  INTEGER(4)     :: DAY
  INTEGER(4)     :: HOUR
  INTEGER(4)     :: MINUTE
  REAL(4)        :: SECOND
END TYPE DATE_TIME
REAL(8)          , SAVE       :: COUNTTURN=0.D0
REAL(8)          , SAVE       :: CURRENT
INTERFACE OPERATOR (.LATER.)
!!$   LOGICAL FUNCTION CLOCK_LATER1(TIME1,TIME2) 
!!$   TYPE (DATE_TIME),INTENT(IN) :: TIME1
!!$   TYPE (DATE_TIME),INTENT(IN) :: TIME2
!!$   END FUNCTION CLOCK_LATER1
  MODULE PROCEDURE CLOCK_LATER1
END INTERFACE
INTERFACE CLOCK$NOW
  MODULE PROCEDURE CLOCK_NOWSTAMP
  MODULE PROCEDURE CLOCK_COUNTER
  MODULE PROCEDURE CLOCK_WHATTIME
END INTERFACE
INTERFACE CLOCK$TIMESTAMP
  MODULE PROCEDURE CLOCK_TIMESTAMP
END INTERFACE
!***********************************************************************
CONTAINS
!       ...............................................................
        SUBROUTINE CLOCK_NOWSTAMP(STAMP)
!       ***************************************************************
!       ** TAKES THE CURRENT TIME AND PRODUCES A PRINTABLE TIMESTAMP **
!       ***************************************************************
        IMPLICIT NONE
        CHARACTER(32),INTENT(OUT) :: STAMP
        TYPE (DATE_TIME)          :: NOW
!       ***************************************************************
        CALL CLOCK_WHATTIME(NOW)
        CALL CLOCK_TIMESTAMP(NOW,STAMP)
        RETURN
        END SUBROUTINE CLOCK_NOWSTAMP
!
!       ...............................................................
        SUBROUTINE CLOCK_WHATTIME(NOW)
!       ***************************************************************
!       **  RETURNS CURRENT DATE AND TIME IN DATE_TIME FORMAT        **
!       ***************************************************************
        IMPLICIT NONE
        TYPE (DATE_TIME), INTENT(OUT) :: NOW
        INTEGER(4)                    :: VALUE(8)
!       ***************************************************************
        CALL DATE_AND_TIME(VALUES=VALUE)
        NOW%YEAR  =VALUE(1)
        NOW%MONTH =VALUE(2)
        NOW%DAY   =VALUE(3)
        NOW%HOUR  =VALUE(5)
        NOW%MINUTE=VALUE(6)
        NOW%SECOND=REAL(VALUE(7))+0.001*REAL(VALUE(8))
        RETURN
        END SUBROUTINE CLOCK_WHATTIME
!
!       ...............................................................
        SUBROUTINE CLOCK_TIMESTAMP(NOW,STAMP) 
!       ***************************************************************
!       ** CONVERTS DATE_TIME TO CHARACTER STRING                    **
!       ***************************************************************
        IMPLICIT NONE
        TYPE (DATE_TIME), INTENT(IN)  :: NOW
        CHARACTER(32),    INTENT(OUT) :: STAMP
        CHARACTER(3),     PARAMETER   :: MONTHS(12) &
     &                   =(/'JAN','FEB','MAR','APR','MAY','JUN' &
     &                     ,'JUL','AUG','SEP','OCT','NOV','DEC'/)
        CHARACTER(3),     PARAMETER    :: DAYS(7) &
     &             =(/'SUN','MON','TUE','WED','THU','FRI','SAT'/)
        INTEGER(4)                    :: Y,M,D
!       ***************************************************************
        IF(NOW%MONTH.LT.1.OR.NOW%MONTH.GT.12) THEN
           CALL ERROR$MSG('MONTH IS NOT IN A VALID RANGE')
           CALL ERROR$STOP('CLOCK: TIMESTAMP')
        END IF
!       ================================================================
!       == CALCULATE DAY OF THE WEEK                                  ==
!       ================================================================
        Y=NOW%YEAR
        M=NOW%MONTH
        D=NOW%DAY
        IF(M.LT.3) THEN
          M=M+12
          Y=Y-1
        END IF
        D=MOD((D+(13*M-27)/5+Y+Y/4-Y/100+Y/400),7)
!       ================================================================
!       ==  PRODUCE STAMP  "SUN 03 DEC 1995 13:15 (16.245S)"         ==
!       ================================================================
        WRITE(STAMP,FMT='(A3,1X,I2.2,1X,A3,1X,I4.4,1X' &
     &     //',I2.2,":",I2.2," (",F6.3,"S)")') &
     &     DAYS(D+1),NOW%DAY,MONTHS(NOW%MONTH),NOW%YEAR &
     &     ,NOW%HOUR,NOW%MINUTE,NOW%SECOND
        RETURN
        END SUBROUTINE CLOCK_TIMESTAMP
!
!       ................................................................
        LOGICAL FUNCTION CLOCK_LATER1(TIME1,TIME2) 
!       ***************************************************************
!       ** CHECKS WHETHER TIME1 IS LATER THAN TIME2                  **
!       ***************************************************************
        IMPLICIT NONE
        TYPE (DATE_TIME),INTENT(IN) :: TIME1
        TYPE (DATE_TIME),INTENT(IN) :: TIME2
        INTEGER(4)                   ::II
!       ***************************************************************
        II=TIME1%YEAR-TIME2%YEAR
        IF(II.NE.0) THEN
          CLOCK_LATER1=(II.GT.0); RETURN
        ELSE 
          II=TIME1%MONTH-TIME2%MONTH
          IF(II.NE.0) THEN
            CLOCK_LATER1=(II.GT.0); RETURN
          ELSE
            II=TIME1%DAY-TIME2%DAY
            IF(II.NE.0) THEN
              CLOCK_LATER1=(II.GT.0); RETURN
            ELSE
              II=TIME1%HOUR-TIME2%HOUR
              IF(II.NE.0) THEN
                CLOCK_LATER1=(II.GT.0); RETURN
              ELSE
                II=TIME1%MINUTE-TIME2%MINUTE
                IF(II.NE.0) THEN
                  CLOCK_LATER1=(II.GT.0); RETURN
                ELSE
                  CLOCK_LATER1=(TIME1%SECOND-TIME2%SECOND.GT.0.0)
                  RETURN
                END IF
              END IF
            END IF
          END IF
        END IF
        RETURN
        END FUNCTION CLOCK_LATER1
!
!       ...............................................................
        SUBROUTINE CLOCK_COUNTER(SECONDS)
!       ***************************************************************
        IMPLICIT NONE
        REAL(8) ,INTENT(OUT) :: SECONDS
        INTEGER(4)           :: COUNT
        INTEGER(4)           :: COUNTRATE
        INTEGER(4)           :: COUNTMAX
!       ***************************************************************
        CALL SYSTEM_CLOCK(COUNT,COUNTRATE,COUNTMAX)
        SECONDS=REAL(COUNT+COUNTTURN*COUNTMAX)/REAL(COUNTRATE)
        IF(SECONDS.LT.CURRENT) THEN
          COUNTTURN=COUNTTURN+1
        END IF
        CURRENT=SECONDS
        RETURN
        END SUBROUTINE CLOCK_COUNTER
END MODULE CLOCK_MODULE



