!........1.........2.........3.........4.........5.........6.........7.........8
MODULE TRACE_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: TRACE                                                              **
!**                                                                           **
!**  PURPOSE: MAINTAINS A STACK OF SUBROUTINES FOR DEBUGGING PURPOSES         **
!**    AND REPORTS THE ACTUAL POSITION OF THE EXECUTION.                      **
!**    THE TRACE COMMANDS CAN BE SWITCHED OFF BY SETTING TOFF=.TRUE.          **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**    TRACE$PUSH(SUBROUTINE_NAME)                                            **
!**    TRACE$PASS(MESSAGE)                                                    **
!**    TRACE$POP                                                              **
!**    TRACE$WRITEHISTORY                                                     **
!**                                                                           **        
!*******************************************************************************
IMPLICIT NONE
TYPE TRACESWITCH_TYPE
  LOGICAL  :: WRITETRACEFILE=.FALSE.
  LOGICAL  :: SILENT=.FALSE.
  LOGICAL  :: FLUSH=.FALSE.
END TYPE TRACESWITCH_TYPE
TYPE(TRACESWITCH_TYPE) ,SAVE :: TRACESWITCH
LOGICAL(4)             ,SAVE :: TOFF=.FALSE.
LOGICAL(4)   ,PARAMETER      :: TFLUSH=.FALSE.
INTEGER(4)   ,PARAMETER      :: MAXLEVEL=50
LOGICAL(4)             ,SAVE :: TFIRST=.TRUE.
CHARACTER(32)          ,SAVE :: HISTORY(MAXLEVEL)
INTEGER(4)             ,SAVE :: LEVEL=0
INTEGER(4)             ,SAVE :: THISTASK=0
INTEGER(4)                   :: NTASKS=1
END MODULE TRACE_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE_FIRST
!      USE STRINGS_MODULE
      USE TRACE_MODULE
      IMPLICIT NONE
!     **************************************************************************
      CALL MPE$QUERY('~',NTASKS,THISTASK)
!
!     ==========================================================================
!     == DEFINE FILE FOR TRACE INFORMATION                                    ==
!     ==========================================================================
      TFIRST=.FALSE.
      END SUBROUTINE TRACE_FIRST
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE_ATTACHTRACEFILE()
!     **************************************************************************
!     ** ATTACH FILE FOR TRACE INFORMATION. THERE WILL BE ONE FILE PER TASK   **
!     ** FOR A PARALLEL JOB                                                   **
!     **************************************************************************
      USE STRINGS_MODULE
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(128)    :: TRACEFILE
      CHARACTER(128)    :: HOSTNAME
      CHARACTER(8)      :: ID
      INTEGER(4)        :: NFILTRACE
      INTEGER(4)        :: STATUS
!     **************************************************************************
      CALL MPE$QUERY('~',NTASKS,THISTASK)
      WRITE(TRACEFILE,*)THISTASK
      TRACEFILE=-'TRACE_'//TRIM(ADJUSTL(TRACEFILE))
      ID=+'TRACE'
      CALL FILEHANDLER$SETFILE(ID,.FALSE.,TRIM(TRACEFILE))
      CALL FILEHANDLER$SETSPECIFICATION(ID,'STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'POSITION','APPEND')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION(ID,'FORM','FORMATTED')
      CALL FILEHANDLER$UNIT(ID,NFILTRACE)
      REWIND NFILTRACE
      CALL GET_ENVIRONMENT_VARIABLE('HOSTNAME',HOSTNAME,STATUS=STATUS)
      WRITE(NFILTRACE,*)'TRACE FILE FROM HOST ',TRIM(HOSTNAME)
      TRACESWITCH%WRITETRACEFILE=.TRUE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE$SETL4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'ON') THEN
        TOFF=.NOT.VAL
!       __ ATTACH TRACEFILE IF IT DOES NOT EXIST YET
        IF(.NOT.TOFF.AND.(.NOT.TRACESWITCH%WRITETRACEFILE)) THEN
          CALL TRACE_ATTACHTRACEFILE()
        END IF
      ELSE IF(ID.EQ.'SILENT') THEN
        TRACESWITCH%SILENT=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRACE$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE$GETL4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ID.EQ.'TRACEFILE') THEN
!       == TRUE IF THE TRACEFILE IS ATTACHED AND CAN BE ACCESSED USING =========
!       == FILEHANDLER$UNIT('TRACE',NFIL) ======================================
        VAL=(.NOT.TOFF).AND.TRACESWITCH%WRITETRACEFILE
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('TRACE$GETL4')
      END IF
      RETURN
      END
!     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE$PUSH(CURRENT_)
!     **************************************************************************
!     **  REGISTER A SUBROUTINE INTO THE TRACE INFORMTAION                    **
!     **  TRACE$$PUSH IS CALLED IN THE BEGINNING OF A SUBROUTINE AND MUST BE  **
!     **  MATCHED BY A TRACE$POP BEFORE LEAVING                               **
!     **************************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: CURRENT_
      REAL(8)                :: MAXMEM
      REAL(8)                :: MBYTE=2.D0**20
      CHARACTER(8)           :: DATE
      CHARACTER(10)          :: TIME
      INTEGER(4)             :: NFILTRACE
      CHARACTER(128)         :: FMT_MEM
!     **************************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      FMT_MEM='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5'
      FMT_MEM=TRIM(ADJUSTL(FMT_MEM))//'" TIME",A8," ",A10)'
!
      CALL USAGE$GET('MAXMEM',MAXMEM)
      CALL DATE_AND_TIME(DATE,TIME)
      LEVEL=LEVEL+1
      IF(LEVEL.LE.MAXLEVEL)THEN
        HISTORY(LEVEL)=CURRENT_
      END IF
!
!     ==========================================================================
!     == WRITE INFORMATION INTO STANDARD OUT                                  ==
!     ==========================================================================
      IF(.NOT.TRACESWITCH%SILENT) THEN
        WRITE(*,FMT='("TRACE-PUSH(",I3,"): LEVEL=",I5," INTO ",A)') &
     &           THISTASK,LEVEL,CURRENT_
        WRITE(*,FMT=FMT_MEM)THISTASK,MAXMEM/MBYTE,DATE,TIME
      END IF
!
!     ==========================================================================
!     == WRITE INFORMATION INTO TRACE FILE                                    ==
!     ==========================================================================
      IF(TRACESWITCH%WRITETRACEFILE.AND..NOT.TRACESWITCH%SILENT) THEN
        CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
        WRITE(NFILTRACE,FMT='("TRACE-PUSH(",I3,"): LEVEL=",I5," INTO ",A)') &
     &           THISTASK,LEVEL,CURRENT_
        WRITE(NFILTRACE,FMT=FMT_MEM)THISTASK,MAXMEM/MBYTE,DATE,TIME
      END IF
!
!     ==========================================================================
!     == FLUSH FILES                                                          ==
!     ==========================================================================
      IF(TRACESWITCH%FLUSH) THEN
        CALL FILEHANDLER$CLOSEALL()
      END IF
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE$POP
!     **************************************************************************
!     ** CHECK A SUBROUTINE OUT OF THE TRACE STACK                            **
!     **************************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      REAL(8)          :: MAXMEM
      REAL(8)          :: MBYTE=2.D0**20
      CHARACTER(8)           :: DATE
      CHARACTER(10)          :: TIME
      INTEGER(4)              :: NFILTRACE
      CHARACTER(128)         :: FMT_MEM
!     **************************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      FMT_MEM='("TRACE-MEM(",I3,"): MAXMEM[MBYTE]=",F10.5'
      FMT_MEM=TRIM(ADJUSTL(FMT_MEM))//'" TIME",A8," ",A10)'
      CALL USAGE$GET('MAXMEM',MAXMEM)
      CALL DATE_AND_TIME(DATE,TIME)
!
!     ==========================================================================
!     == WRITE INFORMATION INTO STANDARD OUT                                  ==
!     ==========================================================================
      IF(.NOT.TRACESWITCH%SILENT) THEN
        IF(LEVEL.GE.1.AND.LEVEL.LE.MAXLEVEL)THEN
          WRITE(*,FMT='("TRACE-POP (",I3,"): LEVEL=",I5," FROM ",A)') &
     &                THISTASK,LEVEL,HISTORY(LEVEL)
          WRITE(*,FMT=FMT_MEM)THISTASK,MAXMEM/MBYTE,DATE,TIME
        ELSE                  
          WRITE(*,FMT='("TRACE-POP(",I3,"): LEVEL=",I5)')THISTASK,LEVEL
          WRITE(*,FMT=FMT_MEM)THISTASK,MAXMEM/MBYTE,DATE,TIME
        END IF
      END IF
!
!     ==========================================================================
!     == WRITE INFORMATION INTO TRACE FILE                                    ==
!     ==========================================================================
      IF(TRACESWITCH%WRITETRACEFILE.AND..NOT.TRACESWITCH%SILENT) THEN
        CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
        IF(LEVEL.GE.1.AND.LEVEL.LE.MAXLEVEL)THEN
          WRITE(NFILTRACE,FMT='("TRACE-POP (",I3,"): LEVEL=",I5," FROM ",A)') &
     &                THISTASK,LEVEL,HISTORY(LEVEL)
          WRITE(NFILTRACE,FMT=FMT_MEM)THISTASK,MAXMEM/MBYTE,DATE,TIME
        ELSE                  
          WRITE(NFILTRACE,FMT='("TRACE-POP(",I3,"): LEVEL=",I5)')THISTASK,LEVEL
          WRITE(NFILTRACE,FMT=FMT_MEM)THISTASK,MAXMEM/MBYTE,DATE,TIME
        END IF
      END IF
!
!     ==========================================================================
!     == FLUSH FILES                                                          ==
!     ==========================================================================
      IF(TRACESWITCH%FLUSH) THEN
        CALL FILEHANDLER$CLOSEALL()
      END IF
!
!     ==========================================================================
!     == REDUCE LEVEL INDICATOR BY ONE                                        ==
!     ==========================================================================
      IF(LEVEL.LE.MAXLEVEL.AND.LEVEL.GE.1)HISTORY(LEVEL)=' '
      LEVEL=LEVEL-1
      IF(LEVEL.LT.0) THEN
        CALL ERROR$MSG('LEVEL<0 ENCOUNTERED')
        CALL ERROR$MSG('LIKELY A TRACE$PUSH IS NOT PAIRED WITH A TRACE$POP')
        CALL ERROR$STOP('TRACE$POP')
      END IF
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE$PASS(MARK)
!     **************************************************************************
!     ** INCLUDE A COMMENT INTO THE TRACE STACK                               **
!     **************************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: MARK
      INTEGER(4)              :: NFILTRACE
!     **************************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
!
!     ==========================================================================
!     == WRITE INFORMATION INTO STANDARD OUT                                  ==
!     ==========================================================================
      IF(.NOT.TRACESWITCH%SILENT) THEN
        WRITE(*,FMT='("TRACE-PASS(",I3,"): LEVEL=",I3," WHERE:",A)') &
     &          THISTASK,LEVEL,MARK
      END IF
!
!     ==========================================================================
!     == WRITE INFORMATION INTO TRACE FILE                                    ==
!     ==========================================================================
      IF(TRACESWITCH%WRITETRACEFILE.AND..NOT.TRACESWITCH%SILENT) THEN
        CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
        WRITE(NFILTRACE,FMT='("TRACE-PASS(",I3,"): LEVEL=",I3," WHERE:",A)') &
     &          THISTASK,LEVEL,MARK
      END IF
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE$R8VAL(NAME,VAL)
!     **************************************************************************
!     ** REPORT A REAL NUMBER INSIDE THE TRACE INFORMATION                    **
!     **************************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      REAL(8)     ,INTENT(IN) :: VAL
      INTEGER(4)              :: NFILTRACE
!     **************************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      IF(TRACESWITCH%SILENT) RETURN
!
!     ==========================================================================
!     == WRITE INFORMATION INTO STANDARD OUT                                  ==
!     ==========================================================================
      WRITE(*,FMT='("TRACE-VALUE",A,": ",F30.20," ON TASK ",I5)') &
     &            TRIM(NAME),VAL,THISTASK 
!
!     ==========================================================================
!     == WRITE INFORMATION INTO TRACE FILE                                    ==
!     ==========================================================================
      IF(TRACESWITCH%WRITETRACEFILE) THEN
        CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
        WRITE(NFILTRACE,FMT='("TRACE-VALUE",A,": ",F30.20," ON TASK ",I5)') &
     &            TRIM(NAME),VAL,THISTASK 
      END IF
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE$C8VAL(NAME,VAL)
!     **************************************************************************
!     ** REPORT A COMPLEX NUMBER INSIDE THE TRACE INFORMATION                 **
!     **************************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      COMPLEX(8)  ,INTENT(IN) :: VAL
      INTEGER(4)              :: NFILTRACE
!     **************************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      IF(TRACESWITCH%SILENT) RETURN
      WRITE(*,FMT='("TRACE-VALUE",A,": (",F20.10,",",F20.10,") ON TASK ",I5)') &
     &            TRIM(NAME),VAL,THISTASK 
!
!     ==========================================================================
!     == WRITE INFORMATION INTO TRACE FILE                                    ==
!     ==========================================================================
      IF(TRACESWITCH%WRITETRACEFILE) THEN
        CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
        WRITE(NFILTRACE, &
     &        FMT='("TRACE-VALUE",A,": (",F20.10,",",F20.10,") ON TASK ",I5)') &
     &            TRIM(NAME),VAL,THISTASK 
      END IF
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE$I4VAL(NAME,VAL)
!     **************************************************************************
!     ** REPORT A INTEGER NUMBER INSIDE THE TRACE INFORMATION                 **
!     **************************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(IN) :: VAL
      INTEGER(4)              :: NFILTRACE
!     **************************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      IF(TRACESWITCH%SILENT) RETURN
      WRITE(*,FMT='("TRACE-VALUE",A,": ",I30," ON TASK ",I5)') &
     &            TRIM(NAME),VAL,THISTASK 
!
!     ==========================================================================
!     == WRITE INFORMATION INTO TRACE FILE                                    ==
!     ==========================================================================
      IF(TRACESWITCH%WRITETRACEFILE) THEN
        CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
        WRITE(NFILTRACE,FMT='("TRACE-VALUE",A,": ",I30," ON TASK ",I5)') &
     &            TRIM(NAME),VAL,THISTASK 
      END IF
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE$L4VAL(NAME,VAL)
!     **************************************************************************
!     ** REPORT A LOGICAL NUMBER INSIDE THE TRACE INFORMATION                 **
!     **************************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      LOGICAL(4)  ,INTENT(IN) :: VAL
      INTEGER(4)              :: NFILTRACE
!     **************************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      IF(TRACESWITCH%SILENT) RETURN
      WRITE(*,FMT='("TRACE-VALUE",A,": ",L1," ON TASK ",I5)') &
     &            TRIM(NAME),VAL,THISTASK 
!
!     ==========================================================================
!     == WRITE INFORMATION INTO TRACE FILE                                    ==
!     ==========================================================================
      IF(TRACESWITCH%WRITETRACEFILE) THEN
        CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
        WRITE(NFILTRACE,FMT='("TRACE-VALUE",A,": ",L1," ON TASK ",I5)') &
     &            TRIM(NAME),VAL,THISTASK 
      END IF
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TRACE$WRITEHISTORY(NFILERR)
!     **************************************************************************
!     ** WRITES THE STACK OF CURRENTLY REGISTERED SUBROUTINES TO THE FILE     **
!     ** WITH UNIT NFILERR                                                    **
!     **************************************************************************
      USE TRACE_MODULE
      IMPLICIT NONE
      INTEGER(4)              :: NFILERR
      INTEGER(4)              :: I
!     **************************************************************************
      IF(TOFF) RETURN
      IF(TFIRST) CALL TRACE_FIRST
      IF(LEVEL.EQ.0) RETURN
      WRITE(NFILERR,FMT='("TRACE-HISTORY:")')
      DO I=1,MIN(LEVEL,MAXLEVEL)
        WRITE(NFILERR,FMT='("LEVEL(",I3,")=",I3," : ",A)')THISTASK,I,HISTORY(I)
      ENDDO
      IF(LEVEL.GT.MAXLEVEL) THEN
        WRITE(NFILERR,FMT='("CURRENT LEVEL(",I3,") IS ",I3)')THISTASK,LEVEL
      END IF
      RETURN 
      END




