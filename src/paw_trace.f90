!......................................................................
MODULE TRACE_MODULE
!**********************************************************************
!**                                                                  **
!**  NAME: TRACE                                                     **
!**                                                                  **
!**  PURPOSE: MAINTAINS A STACK OF SUBROUTINES FOR DEBUGGING PURPOSES**
!**    AND REPORTS THE ACTUAL POSITION OF THE EXECUTION.             **
!**    THE TRACE COMMANDS CAN BE SWITCHED OFF BY SETTING TOFF=.TRUE. **
!**                                                                  **
!**  FUNCTIONS:                                                      **
!**    TRACE$PUSH(SUBROUTINE_NAME)                                   **
!**    TRACE$PASS(MESSAGE)                                           **
!**    TRACE$POP                                                     **
!**    TRACE$WRITEHISTORY                                            **
!**                                                                  **
!**********************************************************************
IMPLICIT NONE
LOGICAL(4)              :: TOFF=.FALSE.
INTEGER(4)   ,PARAMETER :: MAXLEVEL=20
LOGICAL(4)   ,SAVE      :: TFIRST=.TRUE.
CHARACTER(32),SAVE      :: HISTORY(MAXLEVEL)
INTEGER(4)   ,SAVE      :: LEVEL=0
INTEGER(4)   ,SAVE      :: THISTASK=0
INTEGER(4)              :: NTASKS=1
CONTAINS
!**********************************************************************
!     .................................................................
      SUBROUTINE TRACE_FIRST
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      TFIRST=.TRUE.
      END SUBROUTINE TRACE_FIRST
END MODULE TRACE_MODULE
!
!     .................................................................
      SUBROUTINE TRACE$SETL4(ID,VAL)
!     ******************************************************************
!     ******************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'ON') THEN
        TOFF=.NOT.VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRAVCE$SETL4')
      END IF
      RETURN
      END
!     
!     .................................................................
      SUBROUTINE TRACE$PUSH(CURRENT_)
!     ******************************************************************
!     ******************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: CURRENT_
      REAL(8)                :: MAXMEM
      REAL(8)                :: MBYTE=2.D0**20
      CHARACTER(8)           :: DATE
      CHARACTER(10)          :: TIME
!     ******************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      CALL USAGE$GET('MAXMEM',MAXMEM)
      LEVEL=LEVEL+1
      IF(LEVEL.LT.MAXLEVEL)THEN
        HISTORY(LEVEL)=CURRENT_
      END IF
      WRITE(*,FMT='("TRACE-PUSH(",I3,"): LEVEL=",I3," INTO ",A)') &
     &           THISTASK,LEVEL,CURRENT_
      CALL DATE_AND_TIME(DATE,TIME)
      WRITE(*,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5,"TIME")') &
     &             THISTASK,MAXMEM/MBYTE,DATE,TIME
      RETURN 
      END
!
!     .................................................................
      SUBROUTINE TRACE$POP
!     ******************************************************************
!     ******************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      REAL(8)          :: MAXMEM
      REAL(8)          :: MBYTE=2.D0**20
      CHARACTER(8)           :: DATE
      CHARACTER(10)          :: TIME
!     ******************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      CALL USAGE$GET('MAXMEM',MAXMEM)
      CALL DATE_AND_TIME(DATE,TIME)
      IF(LEVEL.GE.1.AND.LEVEL.LE.MAXLEVEL)THEN
        WRITE(*,FMT='("TRACE-POP (",I3,"): LEVEL=",I3," FROM ",A)') &
     &              THISTASK,LEVEL,HISTORY(LEVEL)
        WRITE(*,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5,"TIME")') &
     &                THISTASK,MAXMEM/MBYTE,DATE,TIME
        HISTORY(LEVEL)=' '
      ELSE                  
        WRITE(*,FMT='("TRACE-POP(",I3,"): LEVEL=",I3)')THISTASK,LEVEL
        WRITE(*,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5,"TIME")') &
     &                  THISTASK,MAXMEM/MBYTE,DATE,TIME
      END IF
      LEVEL=LEVEL-1
      RETURN 
      END
!
!     .................................................................
      SUBROUTINE TRACE$PASS(MARK)
!     ******************************************************************
!     ******************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MARK
!     ******************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      WRITE(*,FMT='("TRACE-PASS(",I3,"): LEVEL=",I3," WHERE:",A)') &
     &          THISTASK,LEVEL,MARK
      RETURN 
      END
!
!     .................................................................
      SUBROUTINE TRACE$WRITEHISTORY(NFILERR)
!     ******************************************************************
!     ******************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      INTEGER               ::NFILERR
      INTEGER               ::I
!     ******************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      IF(LEVEL.EQ.0) RETURN
      WRITE(NFILERR,FMT='("TRACE-HISTORY:")')
      DO I=1,MIN(LEVEL,MAXLEVEL)
          WRITE(NFILERR,FMT='("LEVEL(",I3,")=",I3," : ",A)') &
     &          THISTASK,I,HISTORY(I)
      ENDDO
      IF(LEVEL.GT.MAXLEVEL) THEN
        WRITE(NFILERR,FMT='("CURRENT LEVEL(",I3,") IS ",I3)') &
     &        THISTASK,LEVEL
      END IF
      RETURN 
      END




