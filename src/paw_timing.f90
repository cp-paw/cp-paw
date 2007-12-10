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
!**    FUNCTION of the timing object CAN BE INVOKED                   ** 
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
  REAL(8)       :: wallclock
  REAL(8)       :: systime
  REAL(8)       :: usrtime
  INTEGER(4)    :: COUNT
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
!     **  initializes a particular clock                                      **
!     **  note that also the name is empty                                    **
!     **************************************************************************
      USE TIMING_MODULE
      IMPLICIT NONE
      TYPE(CLOCK_TYPE),INTENT(OUT) :: CLOCK_
!     **************************************************************************
      CLOCK_%NAME=' '
      CLOCK_%RUNNING=.FALSE.
      CLOCK_%COUNT=0
      CLOCK_%USED=0.D0
      CLOCK_%usrtime=0.D0
      CLOCK_%systime=0.D0
      CLOCK_%wallclock=0.D0
      RETURN
      END 
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING$START
!     **************************************************************************
!     **  restarts the time object. All prior information is lost.            **
!     **************************************************************************
      USE TIMING_MODULE
      IMPLICIT NONE
      REAL(8) :: TIME
!     **************************************************************************
      STARTED=.TRUE.
      CALL TIMING_CLOCK(TIME)
      BEGINTIME=TIME
      NENTRY=0         !reset the number of  clocks
!
!     ==========================================================================
!     ==  the totalclock runs continuously and is used to determine percentages
!     ==========================================================================
      CALL TIMING$CLOCKON('TOTALCLOCK')
      RETURN
      END 
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING$count
!     **************************************************************************
!     **  advances the counter of total clock  (used for iterative procedures)**
!     **************************************************************************
      USE TIMING_MODULE
      IMPLICIT NONE
!     **************************************************************************
      CALL TIMING$CLOCKOff('TOTALCLOCK')
      CALL TIMING$CLOCKON('TOTALCLOCK')
      RETURN
      END 
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING$CLOCKON(ID_)
!     ******************************************************************
!     ******************************************************************
      USE TIMING_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      CHARACTER(32)           :: ID
      INTEGER(4)              :: IENTRY
      INTEGER(4)              :: I
      REAL(8)                 :: TIME,usrtime,systime
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
      CALL LIB$ETIME(USRTIME,SYSTIME)
      CLOCK(IENTRY)%WALLCLOCK=CLOCK(IENTRY)%WALLCLOCK-TIME
      CLOCK(IENTRY)%USRTIME=CLOCK(IENTRY)%USRTIME-USRTIME
      CLOCK(IENTRY)%SYSTIME=CLOCK(IENTRY)%SYSTIME-SYSTIME
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
      REAL(8)                 :: TIME,usrtime,systime
!     ******************************************************************
      CALL TIMING_CLOCK(TIME)
      CALL LIB$ETIME(USRTIME,SYSTIME)
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
      CLOCK(IENTRY)%WALLCLOCK=CLOCK(IENTRY)%WALLCLOCK+TIME
      CLOCK(IENTRY)%USRTIME=CLOCK(IENTRY)%USRTIME+USRTIME
      CLOCK(IENTRY)%SYSTIME=CLOCK(IENTRY)%SYSTIME+SYSTIME
      RETURN
      END
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TIMING$PRINT(cid,NFIL)
!     **************************************************************************
!     **  REPORT TIMING INFORMATION                                           **
!     **************************************************************************
      USE TIMING_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      CHARACTER(*)    ,INTENT(IN)  :: CID  ! COMMUNICATOR ID (SEE MPE)
      INTEGER(4)      ,INTENT(IN)  :: NFIL
      INTEGER(4)                   :: I,J,K
      INTEGER(4)                   :: NTASKS,THISTASK  
      CHARACTER(15)                :: TIMESTRING(2)
      REAL(8)                      :: TIME
      REAL(8)                      :: PERCENT
      REAL(8)                      :: PERCENTSYS
      REAL(8)         ,PARAMETER   :: R8SMALL=1.D-20
      LOGICAL(4)                   :: ONCLOCK(NENTRYX)
      REAL(8)                      :: WALLCLOCK,WALLCLOCK1
      REAL(8)                      :: USRTIME,USRTIME1
      REAL(8)                      :: SYSTIME,SYSTIME1
      INTEGER(4)                   :: COUNT1
      TYPE(CLOCK_TYPE),ALLOCATABLE :: CLOCKARRAY(:,:)
      CHARACTER(32)   ,ALLOCATABLE :: NAMEARRAY(:,:)
      CHARACTER(32)                :: NAMES(NENTRYX)
      REAL(8)         ,ALLOCATABLE :: TIMEARRAY(:,:)
      REAL(8)                      :: TIMES(NENTRYX)
      INTEGER(4)      ,ALLOCATABLE :: COUNTARRAY(:,:)
      INTEGER(4)                   :: COUNTS(NENTRYX)
!     ************************  P.E. BLOECHL, CLAUSTHAL/GOSLAR 2005  ***********
      IF(.NOT.STARTED) THEN
        CALL ERROR$MSG('WARNING FROM TIMING: CLOCK HAS NOT BEEN STARTED')
        CALL ERROR$MSG('CALL TIMING$START BEFORE ANY OTHER FUNCTON')
        CALL ERROR$STOP('TIMIMG')
      END IF
      CALL MPE$QUERY(CID,NTASKS,THISTASK)
! 
!     ==========================================================================
!     == STOP ALL CLOCKS                                                      ==
!     ==========================================================================
      DO I=1,NENTRY
        ONCLOCK(I)=CLOCK(I)%RUNNING
        IF(ONCLOCK(I))CALL TIMING$CLOCKOFF(CLOCK(I)%NAME)
      ENDDO
!     == CORRECT TOTAL ITERATION COUNT IF SET MANUALLY
      IF(CLOCK(1)%COUNT.GT.1)CLOCK(1)%COUNT=CLOCK(1)%COUNT-1
! 
!     ==========================================================================
!     == COLLECT INFORMATION FROM OTHER PROCESSES                             ==
!     == DATA ARE FILLED INTO CLOCKARRAY                                      ==
!     ==========================================================================
      IF(THISTASK.EQ.1)ALLOCATE(CLOCKARRAY(NENTRYX,NTASKS))
      ALLOCATE(NAMEARRAY(NENTRYX,NTASKS))
      ALLOCATE(TIMEARRAY(NENTRYX,NTASKS))
      ALLOCATE(COUNTARRAY(NENTRYX,NTASKS))
      NAMES(:)=''
      TIMES(:)=0.D0
      COUNTS(:)=0
!
      NAMES(1:NENTRY)=CLOCK(1:NENTRY)%NAME
      CALL MPE$GATHER(CID,1,NAMES,NAMEARRAY)  ! NEED A GATHER ROUTINE FOR STRINGS!!!
      IF(THISTASK.EQ.1)CLOCKARRAY(:,:)%NAME=NAMEARRAY(:,:)
!
      TIMES(1:NENTRY)=CLOCK(1:NENTRY)%USED
      CALL MPE$GATHER(CID,1,TIMES,TIMEARRAY)
      IF(THISTASK.EQ.1)CLOCKARRAY(:,:)%USED=TIMEARRAY(:,:)
!
      TIMES(1:NENTRY)=CLOCK(1:NENTRY)%WALLCLOCK
      CALL MPE$GATHER(CID,1,TIMES,TIMEARRAY)
      IF(THISTASK.EQ.1)CLOCKARRAY(:,:)%WALLCLOCK=TIMEARRAY(:,:)
!
      TIMES(1:NENTRY)=CLOCK(1:NENTRY)%USRTIME
      CALL MPE$GATHER(CID,1,TIMES,TIMEARRAY)
      IF(THISTASK.EQ.1)CLOCKARRAY(:,:)%USRTIME=TIMEARRAY(:,:)
!
      TIMES(1:NENTRY)=CLOCK(1:NENTRY)%SYSTIME
      CALL MPE$GATHER(CID,1,TIMES,TIMEARRAY)
      IF(THISTASK.EQ.1)CLOCKARRAY(:,:)%SYSTIME=TIMEARRAY(:,:)
!
      COUNTS(1:NENTRY)=CLOCK(1:NENTRY)%COUNT
      CALL MPE$GATHER(CID,1,COUNTS,COUNTARRAY)
      IF(THISTASK.EQ.1)CLOCKARRAY(:,:)%COUNT=COUNTARRAY(:,:)
!
      DEALLOCATE(NAMEARRAY)
      DEALLOCATE(TIMEARRAY)
      DEALLOCATE(COUNTARRAY)
!
!     ==================================================================
!     ==  PRINT INFORMATION                                           ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        WRITE(NFIL,FMT='(/''RUN-TIME REPORT''/''==============='')')
        WALLCLOCK=SUM(CLOCKARRAY(1,:)%WALLCLOCK)
        USRTIME=SUM(CLOCKARRAY(1,:)%USRTIME)
        SYSTIME=SUM(CLOCKARRAY(1,:)%SYSTIME)
        WRITE(NFIL,FMT='(''NUMBER OF ITERATIONS    '',I10)')CLOCKARRAY(1,1)%COUNT
        WRITE(NFIL,FMT='(''NUMBER OF PROCESSES     '',I10)')NTASKS
        CALL TIMING_CONVERT(CLOCKARRAY(1,1)%WALLCLOCK,TIMESTRING(1))
        WRITE(NFIL,FMT='(''ELAPSED WALLCLOCK TIME : '',A15)')TIMESTRING(1)
        CALL TIMING_CONVERT(WALLCLOCK,TIMESTRING(1))
        WRITE(NFIL,FMT='(''SUMMED WALLCLOCK TIME  : '',A15)')TIMESTRING(1)
        CALL TIMING_CONVERT(USRTIME+SYSTIME,TIMESTRING(1))
        WRITE(NFIL,FMT='(''TOTAL CPU TIME         : '',A15)')TIMESTRING(1)
        CALL TIMING_CONVERT(USRTIME,TIMESTRING(1))
        WRITE(NFIL,FMT='(''TOTAL USER CPU TIME    : '',A15)')TIMESTRING(1)
        CALL TIMING_CONVERT(SYSTIME,TIMESTRING(1))
        WRITE(NFIL,FMT='(''TOTAL SYSTEM CPU TIME  : '',A15)')TIMESTRING(1)
!
        WRITE(NFIL,FMT='(T1,A,T25,A15,T40,A15'// &
     &                   ',T56,A7,T64,A7,T72,A9)') &
     &           'ID','TOTAL CPU','CPU/#ITER','CPU/TOT' &
     &           ,'SYS/CPU','#CALLS/#ITER'
        WRITE(NFIL,FMT='(80("-"))')
        DO I=1,NENTRY
          WALLCLOCK1=CLOCKARRAY(I,1)%WALLCLOCK
          USRTIME1=CLOCKARRAY(I,1)%USRTIME
          SYSTIME1=CLOCKARRAY(I,1)%SYSTIME
          COUNT1=CLOCKARRAY(I,1)%COUNT
          DO K=2,NTASKS
            DO J=1,NENTRYX
              IF(CLOCKARRAY(I,1)%NAME.EQ.CLOCKARRAY(J,K)%NAME) THEN
                WALLCLOCK1=WALLCLOCK1+CLOCKARRAY(J,K)%WALLCLOCK
                USRTIME1=USRTIME1+CLOCKARRAY(J,K)%USRTIME
                SYSTIME1=SYSTIME1+CLOCKARRAY(J,K)%SYSTIME
                COUNT1=COUNT1+CLOCKARRAY(J,K)%COUNT
                EXIT
              END IF
            ENDDO
          ENDDO
          PERCENT=(USRTIME1+SYSTIME1)/(USRTIME+SYSTIME)*100.D0
          PERCENTSYS=SYSTIME1/(USRTIME1+SYSTIME1+R8SMALL)*100.D0
          TIME=USRTIME1+SYSTIME1
          CALL TIMING_CONVERT(TIME,TIMESTRING(1))
          TIME=(USRTIME1+SYSTIME1)/CLOCK(1)%COUNT
          CALL TIMING_CONVERT(TIME,TIMESTRING(2))
          WRITE(NFIL,FMT='(25("."),T1,A,T25,A15,T40,A15'// &
     &                   ',T58,I3," %",T66,I3," %",T72,I9)') &
     &           TRIM(CLOCK(I)%NAME),(TIMESTRING(J),J=1,2) &
     &          ,NINT(PERCENT),NINT(PERCENTSYS) &
     &          ,NINT(REAL(COUNT1)/REAL(CLOCK(1)%COUNT*NTASKS))
        ENDDO
        DEALLOCATE(CLOCKARRAY)
      END IF
! 
!     ==========================================================================
!     == RESTART ALL CLOCKS THAT WERE RUNNING (RESET ORIGINAL STATE)          ==
!     ==========================================================================
      DO I=1,NENTRY
        IF(ONCLOCK(I))CALL TIMING$CLOCKON(CLOCK(I)%NAME)
      ENDDO
      RETURN
      END  
!
!     ...............................................................
      SUBROUTINE TIMING_CLOCK(SECONDS)
!     ***************************************************************
      IMPLICIT NONE
      REAL(8) ,INTENT(OUT) :: SECONDS
      INTEGER              :: COUNT
      INTEGER,SAVE         :: COUNTPREV=0
      INTEGER              :: COUNTRATE
      INTEGER              :: COUNTMAX
      INTEGER,SAVE         :: COUNTTURN=0
      REAL(8)              :: RCOUNT
!     ***************************************************************
      CALL SYSTEM_CLOCK(COUNT,COUNTRATE,COUNTMAX)
!     == CHECK REPEATCYCLE OF COUNT
      IF(COUNT.LT.COUNTPREV) THEN
        COUNTTURN=COUNTTURN+1
      END IF
      COUNTPREV=COUNT
! 
!     == OPERATIONS ARE DONE IN REAL MODE TO AVOID OVERFLOWS ==========
      RCOUNT=REAL(COUNT)+REAL(COUNTMAX)*REAL(COUNTTURN)
      SECONDS=RCOUNT/REAL(COUNTRATE)
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
      CHARACTER(15),INTENT(OUT):: TIMESTRING 
      INTEGER(4)               :: HOURS,MINUTES,SECONDS,SECONDFRAC
      REAL(8)                  :: SVAR
!     ******************************************************************
      SVAR   =TIME
      HOURS  =INT(SVAR/3600.D0)
      SVAR   =SVAR-DBLE(3600*HOURS)
      MINUTES=INT(SVAR/60.D0)
      SVAR   =SVAR-DBLE(60*MINUTES) 
      SECONDS=INT(SVAR)
      SVAR   =SVAR-DBLE(SECONDS)
      SECONDFRAC=INT(SVAR*10.D0)
      WRITE(TIMESTRING,FMT='(I6,''H'',I2,''M'',I2,''.'',I1,''S'')') &
     &               HOURS,MINUTES,SECONDS,SECONDFRAC
      RETURN
      END


