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
      SUBROUTINE TIMING$PRINT(cid,NFIL,ID_)
!     ******************************************************************
!     **  REPORT TIMING INFORMATION                                   **
!     ******************************************************************
      USE TIMING_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: cid  ! communicator id (see mpe)
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
      REAL(8)                 :: PERCENTIDLE
      character(64),allocatable:: idarr(:)
      real(8)      ,allocatable:: tarr(:)
      real(8)      ,allocatable:: tarrall(:,:)
      integer(4)   ,allocatable:: icountarr(:)
!     ************************  P.E. BLOECHL, TU-CLAUSTHAL 2005  *******
      CALL TIMING_CLOCK(TIME)
      CALL MPE$QUERY(cid,NTASKS,THISTASK)
CALL MPE$SYNC('MONOMER')
print*,thistask,'timing$print start',nfil,time,begintime,cid
      IF(.NOT.STARTED) THEN
        CALL ERROR$MSG('WARNING FroM TIMING: CLOCK HAS NOT BEEN STARTED')
        CALL ERROR$MSG('CALL TIMING$START BEFORE ANY OTHER FUNCTON')
        CALL ERROR$STOP('TIMIMG')
      END IF
!
!     ==================================================================
!     ==  collect timing information from all nodes                   ==
!     ==================================================================
!     == THE FIRST TASK DECIDES, WHICH CLOCKS ARE REPORTED
      ISVAR=NENTRY
      CALL MPE$BROADCAST(CID,1,ISVAR)
      ALLOCATE(IDARR(ISVAR))
      ALLOCATE(TARR(ISVAR))
      ALLOCATE(TARRALL(ISVAR,NTASKS))
      ALLOCATE(ICOUNTARR(ISVAR))
      IF(THISTASK.EQ.1) THEN
        DO I=1,NENTRY
          IDARR(I)=CLOCK(I)%NAME
        ENDDO
      END IF
      CALL MPE$BROADCAST(CID,1,IDARR)
      TARR(:)=0.D0
      ICOUNTARR(:)=0
      DO I=1,ISVAR
        DO J=1,NENTRY
          IF(CLOCK(J)%NAME.EQ.IDARR(I)) THEN
            TARR(I)=TARR(I)+CLOCK(J)%USED
            IF(CLOCK(J)%RUNNING)TARR(I)=TARR(I)+TIME
            ICOUNTARR(I)=ICOUNTARR(I)+CLOCK(I)%COUNT
          END IF
        ENDDO
      ENDDO
      CALL MPE$GATHER(CID,1,TARR,TARRALL)
      CALL MPE$COMBINE(CID,'+',ICOUNTARR)
      DEALLOCATE(IDARR)
      DEALLOCATE(TARR)
!
!     ==================================================================
!     ==  REPORT GLOBAL TIMING                                        ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        TLOCAL=TIME-BEGINTIME
        WRITE(NFIL,FMT='(/''RUN-TIME REPORT''/''==============='')')
        CALL TIMING_CONVERT(TLOCAL,TIMESTRING(1))
        WRITE(NFIL,FMT='(''CPU TIME ON FIRST TASK : '',A12)')TIMESTRING(1)
        CALL TIMING_CONVERT(TLOCAL*REAL(NTASKS,KIND=8),TIMESTRING(1))
        WRITE(NFIL,FMT='(''TOTAL CPU TIME         : '',A12)')TIMESTRING(1)
        WRITE(NFIL,FMT='(''NUMBER OF PROCESSORS   : '',I4)')NTASKS
        TOTAL=TLOCAL*REAL(NTASKS)  ! USED TO OBTAIN PERCENTAGE OIF INDIVIDUAL CLOCKS
!
!       ==================================================================
!       ==  REPORT ON INDIVIDUAL CLOCKS                                 ==
!       ==================================================================
        WRITE(NFIL,FMT='(A,T20," ",A6," ",A12," ",A12," ",A7," ",A7)') &
     &         "NAME","#CALLS","TOTAL","PER CALL","T/TOTAL","IDLE"
!
!       ================================================================
!       ==  ANALYSIS AND PRINTOUT ONLY ON TASK 1                      ==
!       ================================================================
        ALLOCATE(TIMEVAL(NTASKS))
        DO I=1,NENTRY
          TIMEVAL(:)=TARRALL(I,:)
!
          TMAX=0.D0
          TSUM = 0.D0
          DO J=1,NTASKS
            TSUM=TSUM+TIMEVAL(J)
            TMAX=MAX(TMAX,TIMEVAL(J))
          ENDDO
          IF(TSUM.EQ.0.D0) CYCLE  ! NO REPORT IF NO TIME HAS BEEN SPENT
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
            SVAR=TSUM/REAL(ICOUNTARR(I),KIND=8)
          ELSE
            SVAR=0.D0
          END IF
          CALL TIMING_CONVERT(SVAR,TIMESTRING(2))
!
!         ==============================================================
!         ==  PERCENTAGE OF TOTAL TIME USED SINCE STARTING THE CLOCKS ==
!         ==============================================================
          IF(TSUM.NE.0.D0) THEN
            PERCENT=TSUM/TOTAL*100.D0
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
        ENDDO
        DEALLOCATE(TIMEVAL)
!     
      END IF
      DEALLOCATE(TARRALL)
      DEALLOCATE(ICOUNTARR)
PRINT*,THISTASK,'TIMING$PRINT DONE'
      RETURN
      END  
!
!     ...............................................................
      SUBROUTINE TIMING_CLOCK(SECONDS)
!     ***************************************************************
      IMPLICIT NONE
      REAL(8) ,INTENT(OUT) :: SECONDS
      INTEGER(4)           :: COUNT
      INTEGER(4)           :: COUNTRATE
      INTEGER(4)           :: COUNTMAX
      INTEGER(4),SAVE      :: COUNTTURN=0
      REAL(8)              :: COUNTCYCLE
      REAL(8)   ,SAVE      :: CURRENT=0.D0
!     ***************************************************************
      CALL SYSTEM_CLOCK(COUNT,COUNTRATE,COUNTMAX)
      COUNTCYCLE=REAL(COUNTMAX)/REAL(COUNTRATE)
      SECONDS=REAL(COUNT)/REAL(COUNTRATE)+COUNTCYCLE*REAL(COUNTTURN)
      IF(SECONDS.LT.CURRENT) THEN
        COUNTTURN=COUNTTURN+1
      END IF
      CURRENT=SECONDS
      RETURN
      END SUBROUTINE TIMING_CLOCK
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


