       PROGRAM MAIN
!      *************************************************************************
!      ** PROGRAM READS AN INPUT FILE AND PRODUCES A FORTRAN SUBROUTINE       **
!      ** WHICH PRINTS THE INPUT FILE TO A SPECIFIED FORTRAN CHANNEL          **
!      **                                                                     **
!      ** IT IS USED TO CONVERT THE PARMFILE FOR THE PAW INSTALLATION INTO    **
!      ** A FORTRAN ROUTINE, THAT WRITES THE PARMFILE TO A FILE               **
!      **                                                                     **
!      ********************************PETER BLOECHL, GOSLAR NOV. 9, 2010*******
       IMPLICIT NONE
       INTEGER,PARAMETER  :: MAXLEN=80
       INTEGER,PARAMETER  :: LINELEN=1024*4
       INTEGER,PARAMETER  :: OVERHEADLEN=25
       CHARACTER(LINELEN) :: LINE
       CHARACTER(MAXLEN)  :: LINEOUT
       INTEGER            :: LEN,I
       INTEGER,PARAMETER  :: NFILO=6    ! STANDARD IN
       INTEGER,PARAMETER  :: NFILIN=5   ! STANDARD OUT
!      *************************************************************************
!
!      =========================================================================
!      ==  PRODUCE SUBROUTINE VERSION$WRITEPARMFILE                           ==
!      =========================================================================
       WRITE(NFILO,FMT='(A)')'SUBROUTINE VERSION$WRITEPARMFILE(FILE)'
       WRITE(NFILO,FMT='(A)')'IMPLICIT NONE'
       WRITE(NFILO,FMT='(A)')'CHARACTER(*),INTENT(IN) :: FILE'
       WRITE(NFILO,FMT='(A)')'INTEGER     ,PARAMETER  :: NFILO=10'
       WRITE(NFILO,FMT='(A)')'OPEN(NFILO,FILE=FILE)'
       DO
         READ(NFILIN,FMT='(A)',ERR=100,END=100)LINE
         IF(LEN_TRIM(LINE).EQ.0) CYCLE
         LINE=ADJUSTL(LINE)
         LEN=LEN_TRIM(LINE)
         DO I=1,LEN ! REPLACE SINGLE QUOTES BY DOUBLE QUOTES
           IF(IACHAR(LINE(I:I)).EQ.39)LINE(I:I)='"'
         ENDDO
         DO WHILE (LEN.GT.0)
           IF(LEN.GT.MAXLEN-OVERHEADLEN) THEN
             LINEOUT="WRITE(NFILO,FMT='(A)')'" &
    &              //TRIM(LINE(1:MAXLEN-OVERHEADLEN))//"\'"
             LINE=LINE(MAXLEN-OVERHEADLEN+1:LEN)
           ELSE
             LINEOUT="WRITE(NFILO,FMT='(A)')'"//TRIM(LINE(1:LEN))//"'"
             LINE=' '
           END IF
           WRITE(NFILO,FMT='(A)')TRIM(LINEOUT)
           LEN=LEN_TRIM(LINE)
         ENDDO
       ENDDO  
100    CONTINUE
       WRITE(NFILO,FMT='(A)')'CLOSE(NFILO)'
       WRITE(NFILO,FMT='(A)')'RETURN'
       WRITE(NFILO,FMT='(A)')'END'
!
!      =========================================================================
!      ==  PRODUCE SUBROUTINE VERSION$GETPAWDIR                               ==
!      =========================================================================
       WRITE(NFILO,FMT='(A)')'SUBROUTINE VERSION$GETPAWDIR(PAWDIR)'
       WRITE(NFILO,FMT='(A)')'CHARACTER(*),INTENT(OUT):: PAWDIR'
       WRITE(NFILO,FMT='(A)')"PAWDIR='XXPAWDIRXX'" ! TO BE MODIFIED
       WRITE(NFILO,FMT='(A)')'RETURN'
       WRITE(NFILO,FMT='(A)')'END'
!
!      =========================================================================
!      ==  PRODUCE SUBROUTINE VERSION$GETPAWVERSION                           ==
!      =========================================================================
       WRITE(NFILO,FMT='(A)')'SUBROUTINE VERSION$GETPAWVERSIONR(PAWVERSION)'
       WRITE(NFILO,FMT='(A)')'CHARACTER(*),INTENT(OUT):: PAWVERSION'
       WRITE(NFILO,FMT='(A)')"PAWVERSION='XXPAWVERSIONXX'" ! TO BE MODIFIED
       WRITE(NFILO,FMT='(A)')'RETURN'
       WRITE(NFILO,FMT='(A)')'END'
!
       STOP
       END
