
!.......................................................................
      PROGRAM MAIN
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE)  :: LL_CNTL
      TYPE(LL_TYPE)  :: LL_STRC
      CHARACTER(256) :: FILE
      INTEGER(4)     :: NFIL1
      INTEGER(4)     :: NFIL2
      INTEGER(4)     :: I
      LOGICAL(4)     :: TCHK
      INTEGER(4)     :: NARGS
!     **************************************************************
                          CALL TRACE$PUSH('MAIN')
!     ==========================================================================
!     == MPE$INIT MUST BE CALLED ALSO FOR NON-PARALLEL CODES                  ==
!     ==========================================================================
      CALL MPE$INIT

      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS.LT.1)THEN
        CALL ERROR$MSG('ARGUMENT LIST OF EXECUTABLE IS EMPTY')
        CALL ERROR$STOP('INITIALIZEFILEANDLER')
      ENDIF
      CALL GET_COMMAND_ARGUMENT(1,FILE)
!
!     ==================================================================
!     ==  DEFINE FILES                                                ==
!     ==================================================================
      I=INDEX(FILE,'.',BACK=.TRUE.)
      CALL FILEHANDLER$SETROOT(FILE(1:I-1))
      CALL FILEHANDLER$SETFILE('TRAIN',.TRUE.,-'_R.TRA')
      CALL FILEHANDLER$SETSPECIFICATION('TRAIN','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('TRAIN','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('TRAIN','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('TRAIN','FORM','UNFORMATTED')
      CALL FILEHANDLER$SETFILE('TRAOUT',.TRUE.,-'_R.TRA-C')
      CALL FILEHANDLER$SETSPECIFICATION('TRAOUT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('TRAOUT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('TRAOUT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('TRAOUT','FORM','UNFORMATTED')
!
!     ==================================================================
!     ==  READ TRAJECTORY                                             ==
!     ==================================================================
      CALL FILEHANDLER$UNIT('TRAIN',NFIL1)      
      CALL FILEHANDLER$UNIT('TRAOUT',NFIL2)      
      CALL CORRECTTRA(NFIL1,NFIL2)
                          CALL TRACE$POP
      CALL ERROR$NORMALSTOP()
      STOP
      END
!
!     ..................................................................
      SUBROUTINE CORRECTTRA(NFIL1,NFIL2)
!     ******************************************************************
!     **  READ TRAJECTORY FILE *****************************************
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL1
      INTEGER(4),INTENT(IN)  :: NFIL2
      INTEGER(4)             :: LENG
      INTEGER(4)             :: I
      REAL(8)                :: TIME
      INTEGER(4)             :: NSTEP
      INTEGER(4)             :: ISTEP
      INTEGER(4)             :: NAT,NATX
      REAL(8)   ,ALLOCATABLE :: R(:,:)
      REAL(8)   ,ALLOCATABLE :: Q(:)
!     ******************************************************************
!
!     ==================================================================
!     ==  CALCULATE NUMBER OF RECORDS                                 ==
!     ==================================================================
      REWIND(NFIL1)
      NSTEP=0
      NATX=0
      DO 
        READ(NFIL1,END=100)I,TIME,LENG
        NATX=MAX(NATX,LENG/4)
        NSTEP=NSTEP+1
      ENDDO
 100  CONTINUE
      ALLOCATE(R(3,NATX))
      ALLOCATE(Q(NATX))
!
!     ==================================================================
!     ==  CALCULATE #(ATOMS)                                          ==
!     ==================================================================
      REWIND(NFIL1)
      REWIND(NFIL2)
      DO I=1,NSTEP
        READ(NFIL1)Istep,TIME,LENG 
        BACKSPACE(NFIL1)
        NAT=LENG/4
        READ(NFIL1)ISTEP,TIME,LENG,R(:,1:NAT),Q(1:NAT)
        IF(NATX.GT.NAT) THEN
          R(:,NAT+1:NATX)=0.D0
          Q(NAT+1:NATX)=0.D0
        END IF
        WRITE(NFIL2)ISTEP,TIME,4*NATX,R(:,:),Q(:)
      ENDDO
      DEALLOCATE(R)
      DEALLOCATE(Q)
      RETURN
      END
