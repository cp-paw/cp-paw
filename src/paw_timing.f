!     
!.......................................................................
MODULE TIMING_MODULE
!***********************************************************************
!**                                                                   ** 
!**  NAME: TIMING                                                     ** 
!**                                                                   ** 
!**  PURPOSE: USED TO TAKE THE TIME OF CERTAIN PROGRAM PARTS          ** 
!**                                                                   ** 
!**  USAGE:                                                           ** 
!**    1) DEFINE THE ZERO FOR THE TOTAL CLOCK TIME BY CALLING         ** 
!**       TIMING$START                                                ** 
!**    2) TAKE THE TIME OF A PROGRAM SEGMENT BY CALLING               ** 
!**       TIMING$CLOCKON AND TIMING$CLOCKOFF BEFORE AND AFTER         ** 
!**    3) PRINT OUT A REPORT USING TIMING$PRINT                       ** 
!**                                                                   ** 
!**  FUNCTIONS                                                        **
!**    TIMING$START                INITIALIZES TOTAL TIME             **
!**    TIMING$CLOCKON(ID)          STARTS CLOCK FOR IDENTIFIER        **
!**    TIMING$CLOCKOFF(ID)         STOPS CLOCK FOR IDENTIFIER         **
!**    TIMING$PRINT(ID)            PRINT REPORT FOR IDENTIFIER        **
!**    TIMING$PRINT('ALL')         PRINT FULL REPORT                  **
!**                                                                   ** 
!**  METHODS:                                                         ** 
!**    TIMING_CLOCK                                                   ** 
!**    TIMING_CONVERT                                                 ** 
!**                                                                   ** 
!**  REMARKS:                                                         ** 
!**    THE TOTAL TIME MUST BE STARTED BEFORE ANY OTHER                ** 
!**    FUNCTION CAN BE INVOKED                                        ** 
!**                                                                   ** 
!**    THE TIME MAY BE STARTED SEVERAL TIMES; EACH START              ** 
!**    RESETS ALL CLOCKS AND THE TOTAL TIME                           ** 
!**                                                                   ** 
!**    ONLY THE FIRST 32 CHARACTERS OF IDENTIFIER ARE SIGNIFICANT     **   
!**                                                                   ** 
!**    THE ROUTINE TIMING_CLOCK CALLS A SYSTEM ROUTINE                ** 
!**                                                                   ** 
!***********************************************************************
TYPE CLOCK_TYPE
  CHARACTER(32) :: NAME        
  LOGICAL(4)    :: RUNNING
  REAL(8)       :: USED
  INTEGER(4)    :: COUNT
END TYPE CLOCK_TYPE
INTEGER(4),PARAMETER :: NENTRYX=100
TYPE(CLOCK_TYPE)     :: CLOCK(NENTRYX)
LOGICAL(4)           :: STARTED=.FALSE.
INTEGER(4)           :: NENTRY=0
REAL(8)              :: BEGINTIME=0.D0
END MODULE TIMING_MODULE
!
!     ................................................................
      SUBROUTINE TIMING_RESET(CLOCK_)
      USE TIMING_MODULE
      IMPLICIT NONE
      TYPE(CLOCK_TYPE),INTENT(OUT) :: CLOCK_
      CLOCK_%NAME=' '
      CLOCK_%RUNNING=.FALSE.
      CLOCK_%COUNT=0
      CLOCK_%USED=0.D0
      RETURN
      END 
!    
!     ..................................................................
      SUBROUTINE TIMING$START
!     ******************************************************************
!     ******************************************************************
      USE TIMING_MODULE
      IMPLICIT NONE
      REAL(8) :: TIME
!     ******************************************************************
      STARTED=.TRUE.
      CALL TIMING_CLOCK(TIME)
      BEGINTIME=TIME
      NENTRY=0
      RETURN
      END 
!    
!     ..................................................................
      SUBROUTINE TIMING$CLOCKON(ID_)
!     ******************************************************************
!     ******************************************************************
      USE TIMING_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      CHARACTER(32)           :: ID
      INTEGER(4)              :: IENTRY
      INTEGER(4)              :: I
      REAL(8)                 :: TIME
!     ******************************************************************
      IF(.NOT.STARTED) THEN
        CALL ERROR$MSG('WARNING FORM TIMING CLOCK HAS NOT BEEN STARTED')
        CALL ERROR$MSG('CALL TIMING$START BEFORE ANY OTHER FUNCTON')
        CALL ERROR$STOP('TIMING$CLOCKON')
      END IF
!
!     ==================================================================
!     ==  FIND CORRECT CLOCK                                          ==
!     ==================================================================
      ID=ID_
      IENTRY=0
      DO I=1,NENTRY
        IF(ID.EQ.CLOCK(I)%NAME) THEN
          IENTRY=I
          EXIT
        END IF
      ENDDO
!
!     ==================================================================
!     ==  CREATE NEW CLOCK IF NOT FOUND                               ==
!     ==================================================================
      IF(IENTRY.EQ.0) THEN
        IF(NENTRY.GE.NENTRYX) THEN
          CALL ERROR$MSG('WARNING FROM TIMING TOO MANY ENTRIES:')
          CALL ERROR$MSG('INCREASE PARAMETER NENTRIX IN SUBROUTINE TIMING')
          CALL ERROR$STOP('TIMING$CLOCKON')
          RETURN
        END IF
        NENTRY=NENTRY+1
        IENTRY=NENTRY
        CALL TIMING_RESET(CLOCK(IENTRY))
        CLOCK(IENTRY)%NAME=ID
      END IF
!
!     ==================================================================
!     ==  CHECK WHETHER CLOCK IS ALREADY ON                           ==
!     ==================================================================
      IF(CLOCK(IENTRY)%RUNNING) THEN
        PRINT*,'WARNING FROM TIMING'
        PRINT*,'IDENTIFIER WAS NOT CLOCKED OFF FOR CLOCKON'
        PRINT*,'IDENTIFIER:',ID_
        RETURN
      END IF
!
!     ==================================================================
!     ==  SWITCH CLOCK ON                                             ==
!     ==================================================================
      CLOCK(IENTRY)%RUNNING=.TRUE.          
      CLOCK(IENTRY)%COUNT=CLOCK(IENTRY)%COUNT+1
      CALL TIMING_CLOCK(TIME)
      CLOCK(IENTRY)%USED=CLOCK(IENTRY)%USED-TIME
      RETURN
      END
!    
!     ..................................................................
      SUBROUTINE TIMING$CLOCKOFF(ID_)
!     ******************************************************************
!     ==  CLOCK OFF                                                   ==
!     ******************************************************************
      USE TIMING_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      CHARACTER(32)           :: ID
      INTEGER(4)              :: IENTRY
      INTEGER(4)              :: I
      REAL(8)                 :: TIME
!     ******************************************************************
      CALL TIMING_CLOCK(TIME)
      IF(.NOT.STARTED) THEN
        CALL ERROR$MSG('CLOCK HAS NOT BEEN STARTED')
        CALL ERROR$MSG('CALL TIMING$START BEFORE ANY OTHER FUNCTON')
        CALL ERROR$STOP('TIMING$CLOCKOFF')
      END IF
!
!     ==================================================================
!     ==  FIND CORRECT CLOCK                                          ==
!     ==================================================================
      ID=ID_
      IENTRY=0
      DO I=1,NENTRY
        IF(ID.EQ.CLOCK(I)%NAME) THEN
          IENTRY=I
          EXIT
        END IF
      ENDDO
      IF(IENTRY.EQ.0) THEN
        PRINT*,'WARNING FROM TIMING'
        PRINT*,'IDENTIFIER NOT IN THE LIST FOR CLOCKOFF'
        PRINT*,'IDENTIFIER:',ID_
        RETURN
      END IF
!     == CLOCK TIME
      IF(.NOT.CLOCK(IENTRY)%RUNNING) THEN
        PRINT*,'WARNING FROM TIMING'
        PRINT*,'IDENTIFIER WAS NOT CLOCKED ON FOR CLOCKOFF'
        PRINT*,'IDENTIFIER:',ID_
        RETURN
      END IF          
      CLOCK(IENTRY)%RUNNING=.FALSE.
      CLOCK(IENTRY)%USED=CLOCK(IENTRY)%USED+TIME
      RETURN
      END
!    
!     ..................................................................
      SUBROUTINE TIMING$PRINT(NFIL,ID_)
!     ******************************************************************
!     **  report timing information                                   **
!     ******************************************************************
      USE TIMING_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(32)           :: ID
      INTEGER(4)              :: IENTRY
      INTEGER(4)              :: I,J
      REAL(8)                 :: TIME
      INTEGER(4)              :: NTASKS,THISTASK  
      REAL(8)                 :: TOTAL     ! TOTAL TIME SINCE LAST RESET ON THISTASK
      REAL(8)                 :: TSUM      ! TIME SUMMED OVER ALL TASKS
      REAL(8)                 :: TIDLE     ! AVERAGE IDLE TIME
      REAL(8)                 :: TLOCAL    ! TIME USED ON THISTASK
      REAL(8)                 :: TMAX      ! MAX(TIME ON ANY TASK)
      REAL(8)     ,ALLOCATABLE:: TIMEVAL(:)!(NTASKS) TIME ON EACH TASK
      REAL(8)                 :: SVAR
      INTEGER(4)              :: I1,I2
      INTEGER(4)              :: ISVAR
      CHARACTER(12)           :: TIMESTRING(2)
      REAL(8)                 :: PERCENT
      real(8)                 :: percentidle
!     ******************************************************************
      CALL TIMING_CLOCK(TIME)
      CALL MPE$QUERY(NTASKS,THISTASK)
      IF(.NOT.STARTED) THEN
        CALL ERROR$MSG('WARNING FORM TIMING CLOCK HAS NOT BEEN STARTED')
        CALL ERROR$MSG('CALL TIMING$START BEFORE ANY OTHER FUNCTON')
        CALL ERROR$STOP('TIMIMG')
      END IF
!
!     ==================================================================
!     ==  FIND CORRECT CLOCKS                                         ==
!     ==================================================================
      IF(ID_.EQ.'ALL') THEN
        I1=1
        I2=NENTRY
      ELSE
        ID=ID_
        IENTRY=0
        DO I=1,NENTRY
          IF(ID.EQ.CLOCK(I)%NAME) THEN
            IENTRY=I
            EXIT
          END IF
        ENDDO
        IF(IENTRY.EQ.0) THEN
          PRINT*,'WARNING FROM TIMING!'
          PRINT*,'ATTEMPT TO PRINT UNKNOWN CLOCK'
          RETURN
        END IF
        I1=IENTRY
        I2=IENTRY
      END IF
!
!     ==================================================================
!     ==  CHECK CONSISTENCY OF #(CLOCKS) ON ALL TASKS                 ==
!     ==================================================================
      ISVAR=NENTRY
      CALL MPE$BROADCAST(1,ISVAR)
      IF(ISVAR.NE.NENTRY) THEN
        CALL ERROR$MSG('INCONSISTENT TIMING TABLES ON PARALLEL TASKS')
        CALL ERROR$STOP('TIMING$REPORT')
      END IF
!
!     ==================================================================
!     ==  REPORT GLOBAL TIMING                                        ==
!     ==================================================================
      TLOCAL=TIME-BEGINTIME
      ALLOCATE(TIMEVAL(NTASKS))
      CALL MPE$GATHER(1,TLOCAL,TIMEVAL)
      IF(THISTASK.EQ.1) THEN
        TSUM=0.
        TMAX=0.D0
        DO I=1,NTASKS
          TSUM=TSUM+TIMEVAL(I)
          TMAX=MAX(TIMEVAL(I),TMAX)
        ENDDO
        TMAX=MAX(1.D-10,TMAX)    ! AVOID DIVIDE-BY-ZERO
        TIDLE=TMAX-TSUM/DBLE(NTASKS)
        CALL TIMING_CONVERT(TSUM,TIMESTRING(1))
        WRITE(NFIL,FMT='(/''RUN-TIME REPORT''/''==============='')')
        WRITE(NFIL,FMT='(''TOTAL CPU TIME USED  : '',A12)')TIMESTRING(1)
        WRITE(NFIL,FMT='(''PERCENTAGE IDLE      : '',F5.1)')100.D0*TIDLE/TMAX
        WRITE(NFIL,FMT='(''NUMBER OF PROCESSORS : '',I4)')NTASKS
      END IF
      TOTAL=Tmax*real(ntasks)  ! USED TO OBTAIN PERCENTAGE OIF INDIVIDUAL CLOCKS
!
!     ==================================================================
!     ==  REPORT ON INDIVIDUAL CLOCKS                                 ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        WRITE(NFIL,FMT='(A,T20," ",A6," ",A12," ",A12," ",A7," ",A7)') &
     &         "NAME","#CALLS","TOTAL","PER CALL","T/TOTAL","IDLE"
      END IF
      DO I=I1,I2
!
!       ================================================================
!       ==  COLLECT TIME USED FROM ALL TASKS (PROCESSORS)             ==
!       ================================================================
        IF(CLOCK(I)%RUNNING) THEN
          IF(THISTASK.EQ.1) THEN          
            WRITE(NFIL,FMT='(A," HAS NOT BEEN CLOCKED OFF")')CLOCK(I)%NAME
          END IF
          TLOCAL=CLOCK(I)%USED+TIME
        ELSE 
          TLOCAL=CLOCK(I)%USED
        END IF
        CALL MPE$GATHER(1,TLOCAL,TIMEVAL)
!
!       ================================================================
!       ==  ANALYSIS AND PRINTOUT ONLY ON TASK 1                      ==
!       ================================================================
        IF(THISTASK.EQ.1) THEN
          TMAX=0.D0
          TSUM = 0.D0
          DO J=1,NTASKS
            TSUM=TSUM+TIMEVAL(J)
            TMAX=MAX(TMAX,TIMEVAL(J))
          ENDDO
          TMAX=MAX(1.D-10,TMAX)    ! AVOID DIVIDE-BY-ZERO
          PERCENTIDLE=100*(TMAX-TSUM/REAL(NTASKS,KIND=8))/TMAX
!
!         ==============================================================
!         ==  PERCENTAGE OF IDLING TIME AND TOTAL CPU TIME TIME USED  ==
!         ==============================================================
          CALL TIMING_CONVERT(TSUM,TIMESTRING(1))
!
!         ==============================================================
!         ==  TIME PER CALL                                           ==
!         ==============================================================
          IF(CLOCK(I)%COUNT.NE.0) THEN
            SVAR=TSUM/DBLE(CLOCK(I)%COUNT)
          ELSE
            SVAR=0.D0
          END IF
          CALL TIMING_CONVERT(SVAR,TIMESTRING(2))
!
!         ==============================================================
!         ==  PERCENTAGE OF TOTAL TIME USED SINCE STARTING THE CLOCKS ==
!         ==============================================================
          IF(TSUM.NE.0.D0) THEN
            PERCENT=Tsum/TOTAL*100.D0
          ELSE
            PERCENT=0
          END IF
!
!         ==============================================================
!         ==  WRITE REPORT                                            ==
!         ==============================================================
          WRITE(NFIL, &
     &      FMT='(A20," ",I6," ",A12," ",A12," ",F5.1,"% ",F5.1,"%")') &
     &       CLOCK(I)%NAME,CLOCK(I)%COUNT,(TIMESTRING(J),J=1,2) &
     &      ,PERCENT,PERCENTIDLE
        END IF
!     
      ENDDO
      DEALLOCATE(TIMEVAL)
      RETURN
      END  
!    
!     ..................................................................
      SUBROUTINE TIMING_CLOCK(TIME)
!     **                                                              **
!     **  REPORTS CPU TIME IN SECONDS                                 **
!     **                                                              **
!     **  REMARK:                                                     **
!     **    CLOCKTIC IS THE TIME STEP OF THE CLOCK IN SECONDS         **
!     **    (MAY BE MACHINE DEPENDENT)                                **
!     **                                                              **
!     **    TIMES RETURNS IN UNITS OF CLOCKTICS                       **
!     **    1)  USER TIME OF THE PARENT PROCESS                       **
!     **    2)  SYSTEM TIME OF PARENT PROCESS                         **
!     **    3)  USER TIME OF CHILD PROCESSES                          **
!     **    4)  SYSTEM TIME OF CHILD PROCESSES                        **
!     **                                                              **
!     **  SUBROUTINES USED: TIMES (STANDARD C LIBRARY)                **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: TIME       ! USRTIME+SYSTIME IN SECONDS
      INTEGER             :: NTIME(4)
      REAL(8),PARAMETER   :: CLOCKTIC=0.01D0
!     ******************************************************************
      CALL lib$TIMES(NTIME)
      TIME=REAL(SUM(NTIME),KIND=8)*CLOCKTIC
      RETURN
      END 
!    
!     ..................................................................
      SUBROUTINE TIMING_CONVERT(TIME,TIMESTRING)
!     ******************************************************************
!     **  CONVERTS REAL(8) TIME IN SECONDS INTO A STRING             ***
!     **  OF TEH FORM               "000H00M00.0S"                   ***
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)      ,INTENT(IN) :: TIME
      CHARACTER(12),INTENT(OUT):: TIMESTRING 
      INTEGER(4)               :: HOURS,MINUTES,SECONDS,SECONDFRAC
      REAL(8)                  :: SVAR
!     ******************************************************************
      SVAR=TIME
      HOURS=INT(SVAR/3600.D0)
      SVAR=SVAR-DBLE(3600*HOURS)
      MINUTES=INT(SVAR/60.D0)
      SVAR=SVAR-DBLE(60*MINUTES) 
      SECONDS=INT(SVAR)
      SVAR=SVAR-DBLE(SECONDS)
      SECONDFRAC=INT(SVAR*10.D0)
      WRITE(TIMESTRING,FMT='(I3,''H'',I2,''M'',I2,''.'',I1,''S'')') &
     &               HOURS,MINUTES,SECONDS,SECONDFRAC
      RETURN
      END


