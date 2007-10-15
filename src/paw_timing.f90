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
      CHARACTER(*),INTENT(IN) :: cid  ! communicator id (see mpe)
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4)              :: I,J,k,l
      REAL(8)                 :: TIME
      INTEGER(4)              :: NTASKS,THISTASK  
      REAL(8)                 :: TOTAL     ! TOTAL TIME SINCE LAST RESET ON THISTASK
      REAL(8)                 :: TSUM      ! TIME SUMMED OVER ALL TASKS
      REAL(8)                 :: TLOCAL    ! TIME USED ON THISTASK
      REAL(8)                 :: TMAX      ! MAX(TIME ON ANY TASK)
      REAL(8)     ,ALLOCATABLE:: TIMEVAL(:)!(NTASKS) TIME ON EACH TASK
      REAL(8)                 :: SVAR
      INTEGER(4)              :: ISVAR
      CHARACTER(15)           :: TIMESTRING(2)
      REAL(8)                 :: PERCENT
      REAL(8)                 :: PERCENTidle
      REAL(8)                 :: PERCENTsys
      CHARACTER(64),ALLOCATABLE:: IDARR(:)
      REAL(8)      ,ALLOCATABLE:: TARR(:)
      REAL(8)      ,ALLOCATABLE:: TARRALL(:,:)
      INTEGER(4)   ,ALLOCATABLE:: ICOUNTARR(:)
      LOGICAL(4)                   :: ONCLOCK(NENTRYX)
      real(8)                      :: wallclock,wallclock1
      real(8)                      :: usrtime,usrtime1
      real(8)                      :: systime,systime1
      integer(4)                   :: count,count1
      TYPE(CLOCK_TYPE),ALLOCATABLE :: CLOCKARRAY(:,:)
      character(32)   ,allocatable :: namearray(:,:)
      character(32)                :: names(nentryx)
      real(8)         ,allocatable :: timearray(:,:)
      real(8)                      :: times(nentryx)
      integer(4)      ,allocatable :: countarray(:,:)
      integer(4)                   :: counts(nentryx)
!     ************************  P.E. BLOECHL, CLAUSTHAL/Goslar 2005  ***********
      IF(.NOT.STARTED) THEN
        CALL ERROR$MSG('WARNING FroM TIMING: CLOCK HAS NOT BEEN STARTED')
        CALL ERROR$MSG('CALL TIMING$START BEFORE ANY OTHER FUNCTON')
        CALL ERROR$STOP('TIMIMG')
      END IF
      CALL MPE$QUERY(cid,NTASKS,THISTASK)
! 
!     ==========================================================================
!     == stop all clocks                                                      ==
!     ==========================================================================
      DO I=1,NENTRY
        ONCLOCK(I)=CLOCK(I)%RUNNING
        if(onclock(i))CALL TIMING$CLOCKOFF(CLOCK(I)%NAME)
      ENDDO
!     == CORRECT TOTAL ITERATION COUNT IF SET MANUALLY
      IF(CLOCK(1)%COUNT.GT.1)CLOCK(1)%COUNT=CLOCK(1)%COUNT-1
! 
!     ==========================================================================
!     == collect information from other processes                             ==
!     == data are filled into clockarray                                      ==
!     ==========================================================================
      if(thistask.eq.1)allocate(clockarray(nentryx,ntasks))
      allocate(namearray(nentryx,ntasks))
      allocate(timearray(nentryx,ntasks))
      allocate(countarray(nentryx,ntasks))
      names(:)=''
      TIMES(:)=0.D0
      counts(:)=0
!
      names(1:nentry)=clock(1:nentry)%name
      CALL MPE$GATHER(CID,1,names,namearray)  ! need a gather routine for strings!!!
      if(thistask.eq.1)CLOCKARRAY(:,:)%NAME=NAMEARRAY(:,:)
!
      TIMES(1:NENTRY)=CLOCK(1:NENTRY)%USED
      CALL MPE$GATHER(CID,1,TIMES,TIMEARRAY)
      if(thistask.eq.1)CLOCKARRAY(:,:)%USED=TIMEARRAY(:,:)
!
      TIMES(1:NENTRY)=CLOCK(1:NENTRY)%WALLCLOCK
      CALL MPE$GATHER(CID,1,TIMES,TIMEARRAY)
      if(thistask.eq.1)CLOCKARRAY(:,:)%WALLCLOCK=TIMEARRAY(:,:)
!
      TIMES(1:NENTRY)=CLOCK(1:NENTRY)%USRTIME
      CALL MPE$GATHER(CID,1,TIMES,TIMEARRAY)
      if(thistask.eq.1)CLOCKARRAY(:,:)%USRTIME=TIMEARRAY(:,:)
!
      TIMES(1:NENTRY)=CLOCK(1:NENTRY)%SYSTIME
      CALL MPE$GATHER(CID,1,TIMES,TIMEARRAY)
      if(thistask.eq.1)CLOCKARRAY(:,:)%SYSTIME=TIMEARRAY(:,:)
!
      COUNTS(1:NENTRY)=CLOCK(1:NENTRY)%COUNT
      CALL MPE$GATHER(CID,1,COUNTS,COUNTARRAY)
      if(thistask.eq.1)CLOCKARRAY(:,:)%COUNT=COUNTARRAY(:,:)
!
      DEALLOCATE(NAMEARRAY)
      DEALLOCATE(TIMEARRAY)
      DEALLOCATE(COUNTARRAY)
!
!     ==================================================================
!     ==  print information                                           ==
!     ==================================================================
      if(thistask.eq.1) then
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
     &           ,'SYS/CPU','#CALLS/#iter'
        WRITE(NFIL,FMT='(80("-"))')
        DO I=1,NENTRY
          WALLCLOCK1=CLOCKARRAY(I,1)%WALLCLOCK
          usrtime1=clockarray(i,1)%usrtime
          systime1=clockarray(i,1)%systime
          count1=clockarray(i,1)%count
          do k=2,ntasks
            do j=1,nentryx
              if(clockarray(i,1)%name.eq.clockarray(j,k)%name) then
                wallclock1=wallclock1+clockarray(j,k)%wallclock
                usrtime1=usrtime1+clockarray(j,k)%usrtime
                systime1=systime1+clockarray(j,k)%systime
                count1=count1+clockarray(j,k)%count
                exit
              end if
            enddo
          enddo
          percent=(usrtime1+systime1)/(usrtime+systime)*100.d0
          percentsys=systime1/(usrtime1+systime1+tiny(systime1))*100.d0
          time=usrtime1+systime1
          CALL TIMING_CONVERT(TIME,TIMESTRING(1))
          time=(usrtime1+systime1)/clock(1)%count
          CALL TIMING_CONVERT(TIME,TIMESTRING(2))
          WRITE(NFIL,FMT='(25("."),t1,A,t25,A15,t40,A15'// &
     &                   ',t58,i3," %",t66,i3," %",t72,i9)') &
     &           trim(CLOCK(I)%NAME),(TIMESTRING(J),J=1,2) &
     &          ,nint(PERCENT),nint(percentsys) &
     &          ,nint(real(count1)/real(clock(1)%count*NTASKS))
        enddo
        deallocate(clockarray)
      end if
!
!!$!     ==========================================================================
!!$!     ==  here the old routine                                                ==
!!$!     ==========================================================================
!!$!
!!$      CALL TIMING_CLOCK(TIME)
!!$!
!!$!     ==================================================================
!!$!     ==  collect timing information from all nodes                   ==
!!$!     ==================================================================
!!$!     == THE FIRST TASK DECIDES, WHICH CLOCKS ARE REPORTED
!!$      ISVAR=NENTRY
!!$      CALL MPE$BROADCAST(CID,1,ISVAR)
!!$      ALLOCATE(IDARR(ISVAR))
!!$      ALLOCATE(TARR(ISVAR))
!!$      ALLOCATE(TARRALL(ISVAR,NTASKS))
!!$      ALLOCATE(ICOUNTARR(ISVAR))
!!$      IF(THISTASK.EQ.1) THEN
!!$        DO I=1,NENTRY
!!$          IDARR(I)=CLOCK(I)%NAME
!!$        ENDDO
!!$      END IF
!!$      CALL MPE$BROADCAST(CID,1,IDARR)
!!$      TARR(:)=0.D0
!!$      ICOUNTARR(:)=0
!!$      DO I=1,ISVAR
!!$        DO J=1,NENTRY
!!$          IF(CLOCK(J)%NAME.EQ.IDARR(I)) THEN
!!$            TARR(I)=TARR(I)+CLOCK(J)%USED
!!$            IF(CLOCK(J)%RUNNING)TARR(I)=TARR(I)+TIME
!!$            ICOUNTARR(I)=ICOUNTARR(I)+CLOCK(I)%COUNT
!!$          END IF
!!$        ENDDO
!!$      ENDDO
!!$      CALL MPE$GATHER(CID,1,TARR,TARRALL)
!!$      CALL MPE$COMBINE(CID,'+',ICOUNTARR)
!!$      DEALLOCATE(IDARR)
!!$      DEALLOCATE(TARR)
!!$!
!!$!     ==================================================================
!!$!     ==  REPORT GLOBAL TIMING                                        ==
!!$!     ==================================================================
!!$      IF(THISTASK.EQ.1) THEN
!!$        TLOCAL=TIME-BEGINTIME
!!$        WRITE(NFIL,FMT='(/''RUN-TIME REPORT''/''==============='')')
!!$        CALL TIMING_CONVERT(TLOCAL,TIMESTRING(1))
!!$        WRITE(NFIL,FMT='(''CPU TIME ON FIRST TASK : '',A12)')TIMESTRING(1)
!!$        CALL TIMING_CONVERT(TLOCAL*REAL(NTASKS,KIND=8),TIMESTRING(1))
!!$        WRITE(NFIL,FMT='(''TOTAL CPU TIME         : '',A12)')TIMESTRING(1)
!!$        WRITE(NFIL,FMT='(''NUMBER OF PROCESSORS   : '',I4)')NTASKS
!!$        TOTAL=TLOCAL*REAL(NTASKS)  ! USED TO OBTAIN PERCENTAGE OIF INDIVIDUAL CLOCKS
!!$!
!!$!       ==================================================================
!!$!       ==  REPORT ON INDIVIDUAL CLOCKS                                 ==
!!$!       ==================================================================
!!$        WRITE(NFIL,FMT='(A,T20," ",A10," ",A15," ",A15," ",A7," ",A7)') &
!!$     &         "NAME","#CALLS","TOTAL","PER CALL","T/TOTAL","IDLE"
!!$!
!!$!       ================================================================
!!$!       ==  ANALYSIS AND PRINTOUT ONLY ON TASK 1                      ==
!!$!       ================================================================
!!$        ALLOCATE(TIMEVAL(NTASKS))
!!$        DO I=1,NENTRY
!!$          TIMEVAL(:)=TARRALL(I,:)
!!$!
!!$          TMAX=0.D0
!!$          TSUM = 0.D0
!!$          DO J=1,NTASKS
!!$            TSUM=TSUM+TIMEVAL(J)
!!$            TMAX=MAX(TMAX,TIMEVAL(J))
!!$          ENDDO
!!$          IF(TSUM.EQ.0.D0) CYCLE  ! NO REPORT IF NO TIME HAS BEEN SPENT
!!$          TMAX=MAX(1.D-10,TMAX)    ! AVOID DIVIDE-BY-ZERO
!!$          PERCENTIDLE=100*(TMAX-TSUM/REAL(NTASKS,KIND=8))/TMAX
!!$!
!!$!         ==============================================================
!!$!         ==  PERCENTAGE OF IDLING TIME AND TOTAL CPU TIME TIME USED  ==
!!$!         ==============================================================
!!$          CALL TIMING_CONVERT(TSUM,TIMESTRING(1))
!!$!
!!$!         ==============================================================
!!$!         ==  TIME PER CALL                                           ==
!!$!         ==============================================================
!!$          IF(CLOCK(I)%COUNT.NE.0) THEN
!!$            SVAR=TSUM/REAL(ICOUNTARR(I),KIND=8)
!!$          ELSE
!!$            SVAR=0.D0
!!$          END IF
!!$          CALL TIMING_CONVERT(SVAR,TIMESTRING(2))
!!$!
!!$!         ==============================================================
!!$!         ==  PERCENTAGE OF TOTAL TIME USED SINCE STARTING THE CLOCKS ==
!!$!         ==============================================================
!!$          IF(TSUM.NE.0.D0) THEN
!!$            PERCENT=TSUM/TOTAL*100.D0
!!$          ELSE
!!$            PERCENT=0
!!$          END IF
!!$!
!!$!         ==============================================================
!!$!         ==  WRITE REPORT                                            ==
!!$!         ==============================================================
!!$          WRITE(NFIL, &
!!$     &      FMT='(A20," ",I10," ",A15," ",A15," ",F5.1,"% ",F5.1,"%")') &
!!$     &       CLOCK(I)%NAME,CLOCK(I)%COUNT,(TIMESTRING(J),J=1,2) &
!!$     &      ,PERCENT,PERCENTIDLE
!!$        ENDDO
!!$        DEALLOCATE(TIMEVAL)
!!$!     
!!$      END IF
!!$      DEALLOCATE(TARRALL)
!!$      DEALLOCATE(ICOUNTARR)
!!$PRINT*,THISTASK,'TIMING$PRINT DONE'
! 
!     ==========================================================================
!     == restart all clocks that were running (reset original state)          ==
!     ==========================================================================
      DO I=1,NENTRY
        if(ONCLOCK(I))CALL TIMING$CLOCKON(CLOCK(I)%NAME)
      ENDDO
      RETURN
      END  
!
!     ...............................................................
      SUBROUTINE TIMING_CLOCK(SECONDS)
!     ***************************************************************
      IMPLICIT NONE
      REAL(8) ,intent(out) :: SECONDS
      INTEGER              :: COUNT
      INTEGER              :: COUNTRATE
      INTEGER              :: COUNTMAX
      INTEGER,SAVE         :: COUNTTURN=0
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
      CHARACTER(15),INTENT(OUT):: TIMESTRING 
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
      WRITE(TIMESTRING,FMT='(I6,''H'',I2,''M'',I2,''.'',I1,''S'')') &
     &               HOURS,MINUTES,SECONDS,SECONDFRAC
      RETURN
      END


