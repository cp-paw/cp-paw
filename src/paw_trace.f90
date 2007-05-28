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
logical(4)   ,parameter :: tflush=.false.
CONTAINS
!**********************************************************************
!     .................................................................
      SUBROUTINE TRACE_FIRST
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(128)    :: TRACEFILE
      character(8)      :: id
!     *****************************************************************
      CALL MPE$QUERY('~',NTASKS,THISTASK)
!
!     ==================================================================
!     == DEFINE FILE FOR TRACE INFORMATION                          ==
!     ==================================================================
      WRITE(TRACEFILE,*)THISTASK
      TRACEFILE=-'.TRACE_'//TRIM(ADJUSTL(TRACEFILE))
      ID=+'TRACE'
      CALL FILEHANDLER$SETFILE(ID,.TRUE.,TRIM(TRACEFILE))
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
!      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')

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
      integer(4)              :: nfiltrace
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
      WRITE(*,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5,"TIME",A8," ",A10)') &
     &             THISTASK,MAXMEM/MBYTE,DATE,TIME
!
!     =======================================================================
!     == WRITE INFORMATION INTO TRACE FILE                                 ==
!     =======================================================================
      CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
      WRITE(NFILTRACE,FMT='("TRACE-PUSH(",I3,"): LEVEL=",I3," INTO ",A)') &
     &           THISTASK,LEVEL,CURRENT_
      CALL DATE_AND_TIME(DATE,TIME)
      WRITE(NFILTRACE,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5,"TIME",A8," ",A10)') &
     &             THISTASK,MAXMEM/MBYTE,DATE,TIME
!
!     =======================================================================
!     == flush files                                                       ==
!     =======================================================================
      IF(TFLUSH) THEN
        call filehandler$closeall()
      end if
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
      integer(4)              :: nfiltrace
!     ******************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      CALL USAGE$GET('MAXMEM',MAXMEM)
      CALL DATE_AND_TIME(DATE,TIME)
      IF(LEVEL.GE.1.AND.LEVEL.LE.MAXLEVEL)THEN
        WRITE(*,FMT='("TRACE-POP (",I3,"): LEVEL=",I3," FROM ",A)') &
     &              THISTASK,LEVEL,HISTORY(LEVEL)
        WRITE(*,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5,"TIME",a8," ",a10)') &
     &                THISTASK,MAXMEM/MBYTE,DATE,TIME
        HISTORY(LEVEL)=' '
      ELSE                  
        WRITE(*,FMT='("TRACE-POP(",I3,"): LEVEL=",I3)')THISTASK,LEVEL
        WRITE(*,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5,"TIME",a8," ",a10)') &
     &                  THISTASK,MAXMEM/MBYTE,DATE,TIME
      END IF
!
!     =======================================================================
!     == WRITE INFORMATION INTO TRACE FILE                                 ==
!     =======================================================================
      CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
      IF(LEVEL.GE.1.AND.LEVEL.LE.MAXLEVEL)THEN
        WRITE(nfiltrace,FMT='("TRACE-POP (",I3,"): LEVEL=",I3," FROM ",A)') &
     &              THISTASK,LEVEL,HISTORY(LEVEL)
        WRITE(nfiltrace,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5,"TIME",a8," ",a10)') &
     &                THISTASK,MAXMEM/MBYTE,DATE,TIME
        HISTORY(LEVEL)=' '
      ELSE                  
        WRITE(nfiltrace,FMT='("TRACE-POP(",I3,"): LEVEL=",I3)')THISTASK,LEVEL
        WRITE(nfiltrace,FMT='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5,"TIME",a8," ",a10)') &
     &                  THISTASK,MAXMEM/MBYTE,DATE,TIME
      END IF
!
!     =======================================================================
!     == flush files                                                       ==
!     =======================================================================
      if(tflush) then
        call filehandler$closeall()
      end if
!
!     =======================================================================
!     == reduce level indicator by one                                     ==
!     =======================================================================
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
      integer(4)              :: nfiltrace
!     ******************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      WRITE(*,FMT='("TRACE-PASS(",I3,"): LEVEL=",I3," WHERE:",A)') &
     &          THISTASK,LEVEL,MARK
!
!     =======================================================================
!     == WRITE INFORMATION INTO TRACE FILE                                 ==
!     =======================================================================
      CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
      WRITE(nfiltrace,fmT='("TRACE-PASS(",I3,"): LEVEL=",I3," WHERE:",A)') &
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




