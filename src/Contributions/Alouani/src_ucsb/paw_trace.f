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
LOGICAL(4)              :: TOFF=.false.
INTEGER(4)   ,PARAMETER :: MAXLEVEL=20
LOGICAL(4)   ,SAVE      :: TFIRST=.TRUE.
CHARACTER(32),SAVE      :: HISTORY(MAXLEVEL)
INTEGER(4)   ,SAVE      :: LEVEL=0
INTEGER(4)   ,SAVE      :: ITASK=0
INTEGER(4)              :: NTASK=1
CONTAINS
!**********************************************************************
!     .................................................................
      SUBROUTINE TRACE_FIRST
      CALL MPE$QUERY(NTASK,ITASK)
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
      REAL(8)                :: MBYTE=2.d0**20
!     ******************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      CALL USAGE$GET('MAXMEM',MAXMEM)
      LEVEL=LEVEL+1
      IF(LEVEL.LT.MAXLEVEL)THEN
        HISTORY(LEVEL)=CURRENT_
      END IF
      WRITE(*,FMT='("TRACE-PUSH(",I3,"): LEVEL=",I3," INTO ",A)') &
     &           ITASK,LEVEL,CURRENT_
        WRITE(*,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5)')ITASK,MAXMEM/mbyte
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
      REAL(8)          :: MBYTE=2.d0**20
!     ******************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      CALL USAGE$GET('MAXMEM',MAXMEM)
      IF(level.ge.1.and.LEVEL.LE.MAXLEVEL)THEN
        WRITE(*,FMT='("TRACE-POP (",I3,"): LEVEL=",I3," FROM ",A)') &
     &              ITASK,LEVEL,HISTORY(LEVEL)
        WRITE(*,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5)')ITASK,MAXMEM/mbyte
        HISTORY(LEVEL)=' '
      ELSE                  
        WRITE(*,FMT='("TRACE-POP(",I3,"): LEVEL=",I3)')ITASK,LEVEL
        WRITE(*,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5)')ITASK,MAXMEM/mbyte
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
     &          ITASK,LEVEL,MARK
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
     &          ITASK,I,HISTORY(I)
      ENDDO
      IF(LEVEL.GT.MAXLEVEL) THEN
        WRITE(NFILERR,FMT='("CURRENT LEVEL(",I3,") IS ",I3)') &
     &        ITASK,LEVEL
      END IF
      RETURN 
      END




