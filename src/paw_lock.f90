!
!     .....................................................LOCK.........
      MODULE LOCK_MODULE
!     ******************************************************************
!     **  LOCK                                                        **
!     **  IS USED TO PREVENT UNAUTHORIZED USE OF THE EXECUTABLE       **
!     **  1) IT ENFORCES AN EXPIRATION DATE                           **
!     **  2) BLOCKS COMPLETELY SELECTED OPTIONS                       **
!     **                                                              **
!     ** DATE_AND_TIME PROVIDES AN INTEGER ARRAY                      **
!     ** (1) YEAR,(2) MONTH,(3) DAY,                                  **
!     ** (4) HOUR,(5) MINUTES, (6) SECONDS                            **
!     ******************************************************************
      USE CLOCK_MODULE, ONLY : DATE_TIME
      IMPLICIT NONE
      TYPE(DATE_TIME),PARAMETER :: EXPIRATION=DATE_TIME(1997,6,31,0,0,0.0)
      LOGICAL(4)     ,PARAMETER :: TON=.FALSE.
      END MODULE LOCK_MODULE
!
!     .....................................................LOCK.........
      SUBROUTINE LOCK$DISABLE(MSG)
!     ******************************************************************
!     ** DISABLE STOPS THE PROGRAM IF THE LOCK IN SWITCHE ON          **
!     ** IT CAN BE USED TO PREVENT THE USE OF UNAUTHORIZED OPTIONS    **
!     ******************************************************************
      USE LOCK_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MSG
!     ******************************************************************
      IF(TON) THEN
        CALL ERROR$MSG('THIS OPTION HAS BEEN DISABLED')
        CALL ERROR$MSG(MSG)
        CALL ERROR$STOP('LOCK$DISABLE')
      END IF
      RETURN
      END
!
!     .....................................................LOCK.........
      SUBROUTINE LOCK$BREAKPOINT
!     ******************************************************************
!     **  CHECKS THE WHETHER THE EXPRIATION DATE HAS PASSED           **
!     **  AND, IF SO, STOPS THE EXECUTION                             **
!     ******************************************************************
      USE LOCK_MODULE
      USE CLOCK_MODULE
      IMPLICIT NONE
      TYPE(DATE_TIME) :: NOW
!     ******************************************************************
      IF(.NOT.TON) RETURN
      CALL CLOCK$NOW(NOW)      
      IF(NOW.LATER.EXPIRATION) THEN
        CALL ERROR$MSG('THIS VERSION OF THE CP-PAW IS EXPIRED')
        CALL ERROR$STOP('LOCK$BREAKPOINT')      
      END IF
      RETURN
      END
!
!     .....................................................LOCK.........
      SUBROUTINE LOCK$REPORT(NFIL)
!     ******************************************************************
!     ******************************************************************
      USE LOCK_MODULE
      USE CLOCK_MODULE
      IMPLICIT NONE
      CHARACTER(32)           :: DATE=' '
      INTEGER(4)   ,INTENT(IN):: NFIL
!     ******************************************************************
      IF(.NOT.TON) RETURN
      CALL CLOCK$TIMESTAMP(EXPIRATION,DATE)
      CALL REPORT$STRING(NFIL,'THIS CP-PAW VERSION IS ABOUT EXPIRE ON')
      CALL REPORT$STRING(NFIL,DATE)
      RETURN
      END



