!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE TIMING_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: TIMING                                                             **
!**                                                                           **
!**  PURPOSE: USED TO TAKE THE TIME OF CERTAIN PROGRAM PARTS                  **
!**                                                                           **
!**  USAGE:                                                                   **
!**    1) DEFINE THE ZERO FOR THE TOTAL CLOCK TIME BY CALLING                 **
!**       TIMING$START                                                        **
!**    2) TAKE THE TIME OF A PROGRAM SEGMENT BY CALLING                       **
!**       TIMING$CLOCKON AND TIMING$CLOCKOFF BEFORE AND AFTER                 **
!**    3) PRINT OUT A REPORT USING TIMING$PRINT                               **
!**                                                                           **
!**  FUNCTIONS                                                                **
!**    TIMING$START                INITIALIZES TOTAL TIME                     **
!**    TIMING$CLOCKON(ID)          STARTS CLOCK FOR IDENTIFIER                **
!**    TIMING$CLOCKOFF(ID)         STOPS CLOCK FOR IDENTIFIER                 **
!**    TIMING$PRINT(ID)            PRINT REPORT FOR IDENTIFIER                **
!**    TIMING$PRINT('ALL')         PRINT FULL REPORT                          **
!**                                                                           **
!**  METHODS:                                                                 **
!**    TIMING_CLOCK                                                           **
!**    TIMING_CONVERT                                                         **
!**                                                                           **
!**  REMARKS:                                                                 **
!**    THE TOTAL TIME MUST BE STARTED BEFORE ANY OTHER                        **
!**    FUNCTION OF THE TIMING OBJECT CAN BE INVOKED                           **
!**                                                                           **
!**    THE TIME MAY BE STARTED SEVERAL TIMES; EACH START                      **
!**    RESETS ALL CLOCKS AND THE TOTAL TIME                                   **
!**                                                                           **
!**    ONLY THE FIRST 32 CHARACTERS OF IDENTIFIER ARE SIGNIFICANT             **
!**                                                                           **
!**    THE ROUTINE TIMING_CLOCK CALLS A SYSTEM ROUTINE                        **
!**                                                                           **
!**    WHILE THE CLOCK IS ACTIVE, THE START TIME IS SUBTRACTED                **
!**    THIS MAY BE MISLEADING WHILE TESTING                                   **
!**                                                                           **
!**    TIMING$PRINT PRODUCES A LISTING WITH                                   **
!**    1) THE WALLCLOCK TIME ON THE FIRST TASK OBTAINED FROM THE FORTAN       **
!**       INTRINSIC "SYSTEM_CLOCK"                                            **
!**    2) THE TOTAL CPU TIME SUMMED OVER ALL TASKS FROM THE FORTRAN INTRINSIC **
!**       "CPU_TIME"                                                          **
!**    FOR EACH CLOCK                                                         **
!**    3) THE TOTAL CPU TIME SUMMED OVER ALL TASKS                            **
!**    4) THE TOTAL CPU TIME PER ITERATION                                    **
!**    5) THE CPU TIME FOR THIS CLOCK RELATIVE TO THE TOTAL CPU TIME          **
!**                                                                           **
!*******************************************************************************
TYPE CLOCK_TYPE
  CHARACTER(32) :: NAME        
  LOGICAL(4)    :: RUNNING
  INTEGER(4)    :: COUNT
  REAL(8)       :: WALLCLOCK ! FROM FORTRAN INTRINSIC SYSTEM_CLOCK
  REAL(8)       :: CPUTIME   ! FROM FORTRAN INTRINSIC CPU_TIME
END TYPE CLOCK_TYPE
INTEGER(4),PARAMETER :: NENTRYX=100
TYPE(CLOCK_TYPE)     :: CLOCK(NENTRYX)
LOGICAL(4)           :: STARTED=.FALSE.
INTEGER(4)           :: NENTRY=0
REAL(8)              :: BEGINTIME=0.D0
END MODULE TIMING_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING_RESET(CLOCK_)
!     **************************************************************************
!     **  INITIALIZES A SPECIFIED CLOCK                                      **
!     **  NOTE THAT ALSO THE NAME IS EMPTY                                    **
!     **************************************************************************
      USE TIMING_MODULE ,ONLY : CLOCK_TYPE  
      IMPLICIT NONE
      TYPE(CLOCK_TYPE),INTENT(OUT) :: CLOCK_
!     **************************************************************************
      CLOCK_%NAME=' '
      CLOCK_%RUNNING=.FALSE.
      CLOCK_%COUNT=0
      CLOCK_%CPUTIME=0.D0
      CLOCK_%WALLCLOCK=0.D0
      RETURN
      END SUBROUTINE TIMING_RESET
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING$START
!     **************************************************************************
!     **  RESTARTS THE TIME OBJECT. ALL PRIOR INFORMATION IS LOST.            **
!     **************************************************************************
      USE TIMING_MODULE ,ONLY : STARTED &
     &                         ,NENTRY &
     &                         ,BEGINTIME 
      IMPLICIT NONE
      REAL(8) :: TIME
!     **************************************************************************
      STARTED=.TRUE.
      CALL TIMING_CLOCK(TIME)  ! FROM SYSTEM_CLOCK
      BEGINTIME=TIME           ! WALLCLOCK TIME WHEN ALL CLOCKS ARE STARTED
      NENTRY=0                 ! RESET THE NUMBER OF  CLOCKS
!
!     ==========================================================================
!     ==  THE TOTALCLOCK RUNS CONTINUOUSLY AND IS USED TO DETERMINE PERCENTAGES
!     ==========================================================================
      CALL TIMING$CLOCKON('TOTALCLOCK')
      RETURN
      END SUBROUTINE TIMING$START
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING$COUNT
!     **************************************************************************
!     **  ADVANCES THE COUNTER OF TOTAL CLOCK  (USED FOR ITERATIVE PROCEDURES)**
!     **************************************************************************
      IMPLICIT NONE
!     **************************************************************************
      CALL TIMING$CLOCKOFF('TOTALCLOCK')
      CALL TIMING$CLOCKON('TOTALCLOCK')
      RETURN
      END SUBROUTINE TIMING$COUNT
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING$CLOCKON(ID_)
!     **************************************************************************
!     ** ACTIVATE A SPECIFIED CLOCK                                           **
!     **************************************************************************
      USE TIMING_MODULE ,ONLY : STARTED &
     &                         ,NENTRYX &
     &                         ,NENTRY &
     &                         ,CLOCK_TYPE &
     &                         ,CLOCK
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      CHARACTER(32)           :: ID
      INTEGER(4)              :: IENTRY
      INTEGER(4)              :: I
      REAL(8)                 :: TIME
!     **************************************************************************
      IF(.NOT.STARTED) THEN
        CALL ERROR$MSG('CLOCK HAS NOT BEEN STARTED')
        CALL ERROR$MSG('CALL TIMING$START BEFORE ANY OTHER FUNCTON')
        CALL ERROR$STOP('TIMING$CLOCKON')
      END IF
!
!     ==========================================================================
!     ==  FIND CORRECT CLOCK                                                  ==
!     ==========================================================================
      ID=ID_
      IENTRY=0
      DO I=1,NENTRY
        IF(ID.EQ.CLOCK(I)%NAME) THEN
          IENTRY=I
          EXIT
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  CREATE NEW CLOCK IF NOT FOUND                                       ==
!     ==========================================================================
      IF(IENTRY.EQ.0) THEN
        IF(NENTRY.GE.NENTRYX) THEN
          CALL ERROR$MSG('MAXIMUM NUMBER OF CLOCKS EXCEEDED')
          CALL ERROR$MSG('INCREASE PARAMETER NENTRIX IN SUBROUTINE TIMING')
          CALL ERROR$STOP('TIMING$CLOCKON')
        END IF
        NENTRY=NENTRY+1
        IENTRY=NENTRY
        CALL TIMING_RESET(CLOCK(IENTRY))
        CLOCK(IENTRY)%NAME=ID
      END IF
!
!     ==========================================================================
!     ==  CHECK WHETHER CLOCK IS ALREADY ON                                   ==
!     ==========================================================================
      IF(CLOCK(IENTRY)%RUNNING) THEN
        CALL ERROR$MSG('ATTEMPT TO ACTIVATE WHILE IT IS STILL RUNNING')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$MSG('CALL TIMING$CLOCKOFF BEFORE TIMING$CLOCKON')
        CALL ERROR$STOP('TIMING$CLOCKON')
      END IF
!
!     ==========================================================================
!     ==  SWITCH CLOCK ON                                                     ==
!     ==========================================================================
      CLOCK(IENTRY)%RUNNING=.TRUE.          
      CLOCK(IENTRY)%COUNT=CLOCK(IENTRY)%COUNT+1
      CALL TIMING_CLOCK(TIME)
      CLOCK(IENTRY)%WALLCLOCK=CLOCK(IENTRY)%WALLCLOCK-TIME
      CALL CPU_TIME(TIME)
      CLOCK(IENTRY)%CPUTIME  =CLOCK(IENTRY)%CPUTIME  -TIME
      RETURN
      END SUBROUTINE TIMING$CLOCKON
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING$CLOCKOFF(ID_)
!     **************************************************************************
!     **  DEACTIVATE A SPECIFIED CLOCK                                        **
!     **************************************************************************
      USE TIMING_MODULE ,ONLY : STARTED &
     &                         ,NENTRY &
     &                         ,CLOCK_TYPE &
     &                         ,CLOCK
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      CHARACTER(32)           :: ID
      INTEGER(4)              :: IENTRY
      INTEGER(4)              :: I
      REAL(8)                 :: TIME
!     **************************************************************************
      IF(.NOT.STARTED) THEN
        CALL ERROR$MSG('CLOCK HAS NOT BEEN STARTED')
        CALL ERROR$MSG('CALL TIMING$START BEFORE ANY OTHER FUNCTON')
        CALL ERROR$STOP('TIMING$CLOCKOFF')
      END IF
!
!     ==========================================================================
!     ==  FIND CORRECT CLOCK                                                  ==
!     ==========================================================================
      ID=ID_
      IENTRY=0
      DO I=1,NENTRY
        IF(ID.EQ.CLOCK(I)%NAME) THEN
          IENTRY=I
          EXIT
        END IF
      ENDDO
      IF(IENTRY.EQ.0) THEN
        CALL ERROR$MSG('UNKNOWN CLOCK')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$STOP('TIMING$CLOCKOFF')
      END IF
!
!     == CLOCK IS NOT ON: ERROR EXIT ===========================================
      IF(.NOT.CLOCK(IENTRY)%RUNNING) THEN
        CALL ERROR$MSG('ATTEMPT TO STOP A CLOCK THAT IS NOT ACTIVE')
        CALL ERROR$CHVAL('ID',ID_)
        CALL ERROR$MSG('CALL TIMING$CLOCKON BEFORE TIMING$CLOCKOFF')
        CALL ERROR$STOP('TIMING$CLOCKOFF')
      END IF          
!
!     ==========================================================================
!     ==  SWITCH CLOCK OFF                                                    ==
!     ==========================================================================
      CLOCK(IENTRY)%RUNNING=.FALSE.
      CALL TIMING_CLOCK(TIME)  ! WALL CLOCK TIME IN SECONDS
      CLOCK(IENTRY)%WALLCLOCK=CLOCK(IENTRY)%WALLCLOCK+TIME
      CALL CPU_TIME(TIME)      ! CPU TIME IN SECONDS
      CLOCK(IENTRY)%CPUTIME=CLOCK(IENTRY)%CPUTIME+TIME
      RETURN
      END SUBROUTINE TIMING$CLOCKOFF
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING$PRINT(CID,NFIL)
!     **************************************************************************
!     **  REPORT TIMING INFORMATION                                           **
!     **************************************************************************
      USE TIMING_MODULE ,ONLY : NENTRYX &
     &                         ,STARTED &
     &                         ,CLOCK_TYPE &
     &                         ,CLOCK &
     &                         ,NENTRY
      USE CLOCK_MODULE  ,ONLY : CLOCK$NOW
      USE MPE_MODULE    ,ONLY : MPE$GATHER
      IMPLICIT NONE
      CHARACTER(*)    ,INTENT(IN)  :: CID  ! COMMUNICATOR ID (SEE MPE)
      INTEGER(4)      ,INTENT(IN)  :: NFIL
      INTEGER(4)                   :: I,J,K
      INTEGER(4)                   :: NTASKS,THISTASK  
      CHARACTER(15)                :: TIMESTRING(2)
      REAL(8)                      :: TIME
      REAL(8)                      :: PERCENT
      LOGICAL(4)                   :: ONCLOCK(NENTRYX)
      REAL(8)                      :: WALLCLOCK,WALLCLOCK1
      REAL(8)                      :: CPUTIME,CPUTIME1
      INTEGER(4)                   :: COUNT,COUNT1
      TYPE(CLOCK_TYPE),ALLOCATABLE :: CLOCKARRAY(:,:)
      CHARACTER(32)   ,ALLOCATABLE :: NAMEARRAY(:,:)
      CHARACTER(32)                :: NAMES(NENTRYX)
      REAL(8)         ,ALLOCATABLE :: TIMEARRAY(:,:)
      REAL(8)                      :: TIMES(NENTRYX)
      INTEGER(4)      ,ALLOCATABLE :: COUNTARRAY(:,:)
      INTEGER(4)                   :: COUNTS(NENTRYX)
      CHARACTER(30)                :: DOTSTRING,FMTSTRING
      CHARACTER(32)                :: NOWSTAMP
!     ************************  P.E. BLOECHL, CLAUSTHAL/GOSLAR 2005  ***********
      IF(.NOT.STARTED) THEN
        CALL ERROR$MSG('WARNING FROM TIMING: CLOCK HAS NOT BEEN STARTED')
        CALL ERROR$MSG('CALL TIMING$START BEFORE ANY OTHER FUNCTON')
        CALL ERROR$STOP('TIMIMG')
      END IF
      CALL MPE$QUERY(CID,NTASKS,THISTASK)
! 
!     ==========================================================================
!     == STOP ALL CLOCKS AND REMEMBER THOSE THAT HAVE BEEN ACTIVE             ==
!     ==========================================================================
      DO I=1,NENTRY
        ONCLOCK(I)=CLOCK(I)%RUNNING
        IF(ONCLOCK(I))CALL TIMING$CLOCKOFF(CLOCK(I)%NAME)
      ENDDO
! 
!     ==========================================================================
!     == COLLECT INFORMATION FROM OTHER PROCESSES                             ==
!     == DATA ARE FILLED INTO CLOCKARRAY                                      ==
!     ==========================================================================
      IF(THISTASK.EQ.1)ALLOCATE(CLOCKARRAY(NENTRYX,NTASKS))
      ALLOCATE(NAMEARRAY(NENTRYX,NTASKS))
      ALLOCATE(TIMEARRAY(NENTRYX,NTASKS))
      ALLOCATE(COUNTARRAY(NENTRYX,NTASKS))
      NAMES(:)=' '
      TIMES(:)=0.D0
      COUNTS(:)=0
!
      NAMES(1:NENTRY)=CLOCK(1:NENTRY)%NAME
      CALL MPE$GATHER(CID,1,NAMES,NAMEARRAY)
      IF(THISTASK.EQ.1)CLOCKARRAY(:,:)%NAME=NAMEARRAY(:,:)
!
      COUNTS(1:NENTRY)=CLOCK(1:NENTRY)%COUNT
      CALL MPE$GATHER(CID,1,COUNTS,COUNTARRAY)
      IF(THISTASK.EQ.1)CLOCKARRAY(:,:)%COUNT=COUNTARRAY(:,:)
!
      TIMES(1:NENTRY)=CLOCK(1:NENTRY)%WALLCLOCK
      CALL MPE$GATHER(CID,1,TIMES,TIMEARRAY)
      IF(THISTASK.EQ.1)CLOCKARRAY(:,:)%WALLCLOCK=TIMEARRAY(:,:)
!
      TIMES(1:NENTRY)=CLOCK(1:NENTRY)%CPUTIME
      CALL MPE$GATHER(CID,1,TIMES,TIMEARRAY)
      IF(THISTASK.EQ.1)CLOCKARRAY(:,:)%CPUTIME=TIMEARRAY(:,:)
!
      DEALLOCATE(NAMEARRAY)
      DEALLOCATE(TIMEARRAY)
      DEALLOCATE(COUNTARRAY)
!
!     ==========================================================================
!     ==  PRINT INFORMATION                                                   ==
!     ==========================================================================
      IF(THISTASK.EQ.1) THEN
        CALL CLOCK$NOW(NOWSTAMP)
        WRITE(NFIL,*)
        WRITE(NFIL,FMT='("RUN-TIME REPORT [",A,"]")')NOWSTAMP
        WRITE(NFIL,FMT='(50("="))')
!       == THE FIRST CLOCK IS 'TOTALCLOCK' =====================================
        DOTSTRING='(30("."),":",T1,A,T32,'
        FMTSTRING=TRIM(DOTSTRING)//'I10)'
        COUNT=CLOCKARRAY(1,1)%COUNT
        IF(COUNT.NE.1)WRITE(NFIL,FMT=FMTSTRING)'NUMBER OF ITERATIONS',COUNT
        IF(NTASKS.NE.1)WRITE(NFIL,FMT=FMTSTRING)'NUMBER OF TASKS',NTASKS
!
!       == ELAPSED WALL CLOCK TIME ON FIRST TASK ===============================
        FMTSTRING=TRIM(DOTSTRING)//'A15)'
        WALLCLOCK=CLOCKARRAY(1,1)%WALLCLOCK
        CALL TIMING_CONVERT(WALLCLOCK,TIMESTRING(1))
        WRITE(NFIL,FMT=FMTSTRING)'ELAPSED WALLCLOCK TIME',TIMESTRING(1)
! 
!       == CPUTIME SUMMED OVER ALL TASKS =======================================
        CPUTIME=SUM(CLOCKARRAY(1,:)%CPUTIME)
        CALL TIMING_CONVERT(CPUTIME,TIMESTRING(1))
        WRITE(NFIL,FMT=FMTSTRING)'TOTAL CPU TIME',TIMESTRING(1)
!
        WRITE(NFIL,*)
        WRITE(NFIL,FMT='(T1,A,T32,A15,T47,A15,T63,A7)') &
     &           'ID','TOTAL CPU','CPU/#ITER','CPU/TOT'
        WRITE(NFIL,FMT='(70("-"))')
        DO I=1,NENTRY
          COUNT1    =CLOCKARRAY(I,1)%COUNT
          WALLCLOCK1=CLOCKARRAY(I,1)%WALLCLOCK
          CPUTIME1  =CLOCKARRAY(I,1)%CPUTIME
          DO K=2,NTASKS
            DO J=1,NENTRYX  ! FIND CLOCK WITH THE CORRECT ID
              IF(CLOCKARRAY(I,1)%NAME.EQ.CLOCKARRAY(J,K)%NAME) THEN
                COUNT1    =COUNT1    +CLOCKARRAY(J,K)%COUNT
                WALLCLOCK1=WALLCLOCK1+CLOCKARRAY(J,K)%WALLCLOCK
                CPUTIME1  =CPUTIME1  +CLOCKARRAY(J,K)%CPUTIME
              END IF
            ENDDO
          ENDDO
          PERCENT=CPUTIME1/CPUTIME*100.D0 ! COMPARED TO 'TOTALCLOCK'
          TIME=CPUTIME1
          CALL TIMING_CONVERT(TIME,TIMESTRING(1)) ! SUMMED CPU TIME
          TIME=CPUTIME1/CLOCK(1)%COUNT
          CALL TIMING_CONVERT(TIME,TIMESTRING(2)) ! CPU TIME PER ITERATION
          WRITE(NFIL,FMT='(32("."),T1,A,T32,A15,T47,A15,T63,I3," %")') &
     &           TRIM(CLOCK(I)%NAME),(TIMESTRING(J),J=1,2),NINT(PERCENT)
        ENDDO
        DEALLOCATE(CLOCKARRAY)
      END IF
! 
!     ==========================================================================
!     == RESTART ALL CLOCKS THAT WERE RUNNING (RESET ORIGINAL STATE)          ==
!     ==========================================================================
      DO I=1,NENTRY
        IF(ONCLOCK(I)) THEN
          CALL TIMING$CLOCKON(CLOCK(I)%NAME)
          CLOCK(1)%COUNT=CLOCK(1)%COUNT-1  ! DO NOT COUNT INTERRUPTION
        END IF
      ENDDO
      RETURN
      END SUBROUTINE TIMING$PRINT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING_CLOCK(SECONDS)
!     **************************************************************************
!     ** RETURNS THE WALL CLOCK TIME IN SECONDS FROM SYSTEM_CLOCK             **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8) ,INTENT(OUT) :: SECONDS
      INTEGER(8)           :: COUNT
      INTEGER(8),SAVE      :: COUNTPREV=0
      INTEGER(8)           :: COUNTRATE
      INTEGER(8)           :: COUNTMAX
      INTEGER(8),SAVE      :: COUNTTURN=0
      REAL(8)              :: RCOUNT
!     **************************************************************************
      CALL SYSTEM_CLOCK(COUNT,COUNTRATE,COUNTMAX)
!     == CHECK REPEATCYCLE OF COUNT
      IF(COUNT.LT.COUNTPREV) THEN
        COUNTTURN=COUNTTURN+1
      END IF
      COUNTPREV=COUNT
! 
!     == OPERATIONS ARE DONE IN REAL MODE TO AVOID OVERFLOWS ===================
      RCOUNT=REAL(COUNT)+REAL(COUNTMAX,KIND=8)*REAL(COUNTTURN,KIND=8)
      SECONDS=RCOUNT/REAL(COUNTRATE,KIND=8)
      RETURN
      END SUBROUTINE TIMING_CLOCK
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING_CONVERT(TIME,TIMESTRING)
!     **************************************************************************
!     **  CONVERTS REAL(8) TIME IN SECONDS INTO A STRING                      **
!     **  OF TEH FORM               "000H00M00.0S"                            **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)      ,INTENT(IN) :: TIME
      CHARACTER(15),INTENT(OUT):: TIMESTRING 
      INTEGER(8)               :: HOURS,MINUTES,SECONDS,SECONDFRAC
      REAL(8)                  :: SVAR
!     **************************************************************************
      IF(TIME.GE.REAL(HUGE(HOURS)*3600.D0,KIND=8)) THEN
        CALL ERROR$MSG('TIME TOO LARGE FOR HANDLING WITH INTEGER(8)')
        CALL ERROR$R8VAL('TIME',TIME)
        CALL ERROR$R8VAL('HUGE(INTEGER*8)',REAL(HUGE(HOURS),KIND=8))
        CALL ERROR$STOP('TIMING_CONVERT')
      END IF
      SVAR   =TIME
      HOURS  =INT(SVAR/3600.D0,KIND=8)
      SVAR   =SVAR-REAL(3600.D0*HOURS,KIND=8)
      MINUTES=INT(SVAR/60.D0,KIND=8)
      SVAR   =SVAR-REAL(60.D0*MINUTES,KIND=8) 
      SECONDS=INT(SVAR,KIND=8)
      SVAR   =SVAR-REAL(SECONDS,KIND=8)
      SECONDFRAC=INT(SVAR*10.D0,KIND=8)
      WRITE(TIMESTRING,FMT='(I6,"H",I2,"M",I2,".",I1,"S")') &
     &               HOURS,MINUTES,SECONDS,SECONDFRAC
      RETURN
      END SUBROUTINE TIMING_CONVERT

#IF DEFINED(CPPVAR_ACCEL_PROFILE)
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE ACCELPROFILE_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: ACCELPROFILE                                                       **
!**                                                                           **
!**  PURPOSE: COLLECT LOW-OVERHEAD SIZE/TIME DATA FOR NUMERICAL HOTSPOTS      **
!**    THAT ARE CANDIDATES FOR NVPL, cuBLAS/cuSOLVER, cuFFT OR cuFFTMp.        **
!**                                                                           **
!*******************************************************************************
TYPE ACCELPROFILE_ENTRY_TYPE
  CHARACTER(32) :: NAME=' '
  INTEGER(8)    :: N1=0
  INTEGER(8)    :: N2=0
  INTEGER(8)    :: N3=0
  INTEGER(8)    :: N4=0
  INTEGER(8)    :: CALLS=0
  REAL(8)       :: SECONDS=0.D0
  REAL(8)       :: MAXSECONDS=0.D0
  REAL(8)       :: FLOPS=0.D0
  REAL(8)       :: BYTES=0.D0
END TYPE ACCELPROFILE_ENTRY_TYPE
INTEGER(4),PARAMETER              :: NENTRYX=512
TYPE(ACCELPROFILE_ENTRY_TYPE)     :: ENTRY(NENTRYX)
INTEGER(4)                        :: NENTRY=0
LOGICAL(4)                        :: INITIALIZED=.FALSE.
LOGICAL(4)                        :: ENABLED=.TRUE.
END MODULE ACCELPROFILE_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ACCELPROFILE_INIT
      USE ACCELPROFILE_MODULE
      IMPLICIT NONE
      CHARACTER(32) :: VAL
      INTEGER(4)    :: STATUS
!     **************************************************************************
      IF(INITIALIZED) RETURN
      ENABLED=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_ACCEL_PROFILE',VAL,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        VAL=ADJUSTL(VAL)
        IF(VAL(1:1).EQ.'0'.OR.VAL(1:1).EQ.'f'.OR.VAL(1:1).EQ.'F' &
     &     .OR.VAL(1:1).EQ.'n'.OR.VAL(1:1).EQ.'N') ENABLED=.FALSE.
      END IF
      INITIALIZED=.TRUE.
      RETURN
      END SUBROUTINE ACCELPROFILE_INIT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ACCELPROFILE$NOW(SECONDS)
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: SECONDS
      INTEGER(8)          :: COUNT
      INTEGER(8)          :: COUNTRATE
      INTEGER(8)          :: COUNTMAX
!     **************************************************************************
      CALL SYSTEM_CLOCK(COUNT,COUNTRATE,COUNTMAX)
      SECONDS=REAL(COUNT,KIND=8)/REAL(COUNTRATE,KIND=8)
      RETURN
      END SUBROUTINE ACCELPROFILE$NOW
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ACCELPROFILE$ADD(NAME,N1,N2,N3,N4,FLOPS,BYTES,SECONDS)
      USE ACCELPROFILE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(8)  ,INTENT(IN) :: N1
      INTEGER(8)  ,INTENT(IN) :: N2
      INTEGER(8)  ,INTENT(IN) :: N3
      INTEGER(8)  ,INTENT(IN) :: N4
      REAL(8)     ,INTENT(IN) :: FLOPS
      REAL(8)     ,INTENT(IN) :: BYTES
      REAL(8)     ,INTENT(IN) :: SECONDS
      INTEGER(4)              :: I
      INTEGER(4)              :: IENTRY
      CHARACTER(32)           :: NAME32
!     **************************************************************************
      IF(.NOT.INITIALIZED) CALL ACCELPROFILE_INIT
      IF(.NOT.ENABLED) RETURN
      NAME32=NAME
      IENTRY=0
      DO I=1,NENTRY
        IF(ENTRY(I)%NAME.EQ.NAME32 &
     &     .AND.ENTRY(I)%N1.EQ.N1 &
     &     .AND.ENTRY(I)%N2.EQ.N2 &
     &     .AND.ENTRY(I)%N3.EQ.N3 &
     &     .AND.ENTRY(I)%N4.EQ.N4) THEN
          IENTRY=I
          EXIT
        END IF
      ENDDO
      IF(IENTRY.EQ.0) THEN
        IF(NENTRY.GE.NENTRYX) RETURN
        NENTRY=NENTRY+1
        IENTRY=NENTRY
        ENTRY(IENTRY)%NAME=NAME32
        ENTRY(IENTRY)%N1=N1
        ENTRY(IENTRY)%N2=N2
        ENTRY(IENTRY)%N3=N3
        ENTRY(IENTRY)%N4=N4
      END IF
      ENTRY(IENTRY)%CALLS=ENTRY(IENTRY)%CALLS+1
      ENTRY(IENTRY)%SECONDS=ENTRY(IENTRY)%SECONDS+SECONDS
      ENTRY(IENTRY)%MAXSECONDS=MAX(ENTRY(IENTRY)%MAXSECONDS,SECONDS)
      ENTRY(IENTRY)%FLOPS=ENTRY(IENTRY)%FLOPS+FLOPS
      ENTRY(IENTRY)%BYTES=ENTRY(IENTRY)%BYTES+BYTES
      RETURN
      END SUBROUTINE ACCELPROFILE$ADD
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ACCELPROFILE$REPORT(CID,NFIL)
      USE ACCELPROFILE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: CID
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(256)          :: FILEBASE
      CHARACTER(256)          :: FILENAME
      CHARACTER(16)           :: RANKSTRING
      INTEGER(4)              :: STATUS
      INTEGER(4)              :: NFILCSV
      INTEGER(4)              :: NTASKS
      INTEGER(4)              :: THISTASK
      INTEGER(4)              :: I
      INTEGER(4)              :: J
      INTEGER(4)              :: BEST
      LOGICAL(4)              :: USED(NENTRYX)
      REAL(8)                 :: AVGSECONDS
      REAL(8)                 :: GFLOPS
      REAL(8)                 :: GBYTES
      REAL(8)                 :: RATE
      REAL(8)                 :: BANDWIDTH
!     **************************************************************************
      IF(.NOT.INITIALIZED) CALL ACCELPROFILE_INIT
      IF(.NOT.ENABLED) RETURN
      CALL MPE$QUERY(CID,NTASKS,THISTASK)
      FILEBASE='cppaw_accel_profile'
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_ACCEL_PROFILE_FILE',FILEBASE &
     &                             ,STATUS=STATUS)
      FILEBASE=ADJUSTL(FILEBASE)
      IF(NTASKS.GT.1) THEN
        WRITE(RANKSTRING,FMT='(".rank",I5.5,".csv")')THISTASK
        FILENAME=TRIM(FILEBASE)//TRIM(RANKSTRING)
      ELSE
        FILENAME=TRIM(FILEBASE)//'.csv'
      END IF
      OPEN(NEWUNIT=NFILCSV,FILE=TRIM(FILENAME),STATUS='REPLACE' &
     &    ,ACTION='WRITE',FORM='FORMATTED')
      WRITE(NFILCSV,FMT='(A)') &
     & 'op,n1,n2,n3,n4,calls,total_seconds,max_seconds,avg_seconds,'// &
     & 'gflop,gbyte,measured_gflop_per_s,measured_gbyte_per_s'
      DO I=1,NENTRY
        AVGSECONDS=ENTRY(I)%SECONDS/REAL(MAX(ENTRY(I)%CALLS,1_8),KIND=8)
        GFLOPS=ENTRY(I)%FLOPS*1.D-9
        GBYTES=ENTRY(I)%BYTES*1.D-9
        RATE=0.D0
        BANDWIDTH=0.D0
        IF(ENTRY(I)%SECONDS.GT.0.D0) THEN
          RATE=GFLOPS/ENTRY(I)%SECONDS
          BANDWIDTH=GBYTES/ENTRY(I)%SECONDS
        END IF
        WRITE(NFILCSV,FMT='(A,",",I0,",",I0,",",I0,",",I0,",",I0,'// &
     &     '",",ES16.8,",",ES16.8,",",ES16.8,",",ES16.8,",",ES16.8,'// &
     &     '",",ES16.8,",",ES16.8)') &
     &     TRIM(ENTRY(I)%NAME),ENTRY(I)%N1,ENTRY(I)%N2,ENTRY(I)%N3 &
     &     ,ENTRY(I)%N4,ENTRY(I)%CALLS,ENTRY(I)%SECONDS &
     &     ,ENTRY(I)%MAXSECONDS,AVGSECONDS,GFLOPS,GBYTES,RATE,BANDWIDTH
      ENDDO
      CLOSE(NFILCSV)
      IF(THISTASK.EQ.1) THEN
        WRITE(NFIL,*)
        WRITE(NFIL,FMT='(80("="))')
        WRITE(NFIL,FMT='("ACCELERATOR PORTING PROFILE")')
        WRITE(NFIL,FMT='(80("="))')
        IF(NTASKS.GT.1) THEN
          WRITE(NFIL,FMT='("CSV FILES: ",A,".rankNNNNN.csv")')TRIM(FILEBASE)
        ELSE
          WRITE(NFIL,FMT='("CSV FILE: ",A)')TRIM(FILENAME)
        END IF
        WRITE(NFIL,FMT='("Set CPPAW_ACCEL_PROFILE=0 to disable at run time.")')
        WRITE(NFIL,*)
        WRITE(NFIL,FMT='(T1,A,T28,A,T40,A,T52,A,T64,A)') &
     &       'HOTSPOT','CALLS','SECONDS','GFLOP','GBYTE'
        WRITE(NFIL,FMT='(80("-"))')
        USED(:)=.FALSE.
        DO I=1,MIN(20,NENTRY)
          BEST=0
          DO J=1,NENTRY
            IF(USED(J)) CYCLE
            IF(BEST.EQ.0) THEN
              BEST=J
            ELSE IF(ENTRY(J)%SECONDS.GT.ENTRY(BEST)%SECONDS) THEN
              BEST=J
            END IF
          ENDDO
          IF(BEST.EQ.0) EXIT
          USED(BEST)=.TRUE.
          WRITE(NFIL,FMT='(T1,A,T28,I10,T40,F10.3,T52,F10.3,T64,F10.3)') &
     &       TRIM(ENTRY(BEST)%NAME),ENTRY(BEST)%CALLS,ENTRY(BEST)%SECONDS &
     &       ,ENTRY(BEST)%FLOPS*1.D-9,ENTRY(BEST)%BYTES*1.D-9
        ENDDO
      END IF
      RETURN
      END SUBROUTINE ACCELPROFILE$REPORT
#ENDIF
