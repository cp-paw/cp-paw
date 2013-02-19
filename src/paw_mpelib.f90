@PROCESS NOEXTCHK 
!*******************************************************************************
!**                                                                           **
!**  NAME: MPELIB                                                             **
!**                                                                           **
!**  PURPOSE: PERFORM COMMUNICATIONS AMONG PROCESSORS FOR                     **
!**   PARALLEL PROCESSING.                                                    **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**    MPE$INIT                                                               **
!**    MPE$QUERY(CID,NTASKS,THISTASK)                                         **
!**    MPE$SYNC(CID)                                                          **
!**    MPE$STOPALL(RETURNCODE)                                                **
!**    MPE$EXIT                                                               **
!**    MPE$ISPARALLEL(IPARALLEL) ??                                           **
!**    MPE$COMBINE(CID,OPERATION,VAL)                                         **
!**    MPE$BROADCAST(CID,FROMTASK,VAL)                                        **
!**    MPE$SEND(CID,TOTASK,TAG,VAL)                                           **
!**    MPE$RECEIVE(CID,FROMTASK,TAG,VAL)                                      **
!**    MPE$TRANSPOSE(CID,VAL)                                                 **
!**    MPE$GATHER(CID,TOTASK,INVAL,OUTVAL)                                    **
!**                                                                           **
!**  REMARKS:                                                                 **
!**    1) INTERFACE TO MPI (MESSAGE PASSING LIBRARY)                          **
!**    2) THE C-PREPROCESSOR IS USED TO CREATE BOTH THE SINGLE                **
!**       AND THE MULTIPROCESSOR VERSION OUT OF THIS FILE.                    **
!**       THUS THIS FILE MUST NOT BE MADE UPPERCASE                           **
!**    3) THIS MODULE REQUIRES THE FILE MPIF90.H IN THE PATH FOR              **
!**       INCLUDE FILES                                                       **
!**    4) CID STANDS FOR COMMUNICATOR ID AND DETERMINES THE GROUP OF          **
!**       PROCESSORS  ON WHICH THE COMMAND IS EXECUTED                        **
!**                                                                           **
!**                                          ERNST NUSTERER 1994/1996         **
!**                                MODIFIED: PETER E. BLOECHL, (1996)         **
!*******************************************************************************
!
!===============================================================================
!==INCLUDEFILE FOR MPI                                                        ==
!===============================================================================
MODULE MPI_MINE
#IFDEF CPPVARIABLE_PARALLEL
!  INCLUDE FILE IS UPPERCASE ON PURPOSE. 
   INCLUDE 'MPIF.H'    
#ENDIF
END MODULE MPI_MINE
!
MODULE MPE_MPIF_MODULE
#IFDEF CPPVARIABLE_PARALLEL
  USE MPI_MINE
!  USE MPI
#ENDIF
  TYPE THISTYPE
    CHARACTER(128)     :: ID
    INTEGER(4)         :: THISTASK
    INTEGER(4)         :: NTASKS
    INTEGER            :: COMM
    INTEGER   ,POINTER :: SENDTAG(:)
    INTEGER   ,POINTER :: RECEIVETAG(:)
    TYPE(THISTYPE),POINTER:: NEXT
  END TYPE THISTYPE
  TYPE MPECLOCKTYPE
    character(16) :: id
    INTEGER(4) :: COUNT
    REAL(8)    :: SIZE
    REAL(8)    :: SIZESQUARE
    REAL(8)    :: MINSIZE
    REAL(8)    :: MAXSIZE
    REAL(8)    :: TIME
    REAL(8)    :: TIMESQUARE
    REAL(8)    :: MINTIME
    REAL(8)    :: MAXTIME
    REAL(8)    :: STARTTIME
  END TYPE MPECLOCKTYPE
  LOGICAL(4)      :: TINI=.FALSE.
  INTEGER         :: COMM
  INTEGER         :: NTASKS
  INTEGER         :: THISTASK
  INTEGER         :: THISTASK_WORLD
  INTEGER,POINTER :: SENDTAG(:)
  INTEGER,POINTER :: RECEIVETAG(:)
  TYPE(THISTYPE),POINTER :: FIRSTTHIS
  TYPE(THISTYPE),POINTER :: THIS
  TYPE(MPECLOCKTYPE) :: COMBINE_CLOCK
  TYPE(MPECLOCKTYPE) :: BROADCAST_CLOCK
  TYPE(MPECLOCKTYPE) :: SENDRECEIVE_CLOCK
  TYPE(MPECLOCKTYPE) :: SEND_CLOCK
  TYPE(MPECLOCKTYPE) :: RECEIVE_CLOCK
  TYPE(MPECLOCKTYPE) :: TRANSPOSE_CLOCK
  TYPE(MPECLOCKTYPE) :: GATHER_CLOCK
  TYPE(MPECLOCKTYPE) :: SYNC_CLOCK
#IFDEF CPPVARIABLE_PARALLEL
!__DATATYPES CONNECTED WITH KIND PARAMETERS
!__FOR MPI1: USE INTEGER(4)=INTEGER*4, REAL(8)=INTEGER*8, COMPLEX(8)=COMPLEX*16
!__FOR MPI2: RECALCULATE VALUES IN MPE$INIT
  INTEGER,SAVE    :: MY_MPITYPE_INTEGER_KIND4=MPI_INTEGER4
  INTEGER,SAVE    :: MY_MPITYPE_REAL_KIND4=MPI_REAL4
  INTEGER,SAVE    :: MY_MPITYPE_REAL_KIND8=MPI_REAL8
  INTEGER,SAVE    :: MY_MPITYPE_COMPLEX_KIND8=MPI_DOUBLE_COMPLEX
CONTAINS
!     .............................................................................
      SUBROUTINE MPE_CHARTOASCII(LENCH,STRING,LENI,IASCII)
      IMPLICIT  NONE
      CHARACTER(*),INTENT(IN) :: STRING(LENCH)
      INTEGER     ,INTENT(IN) :: LENCH
      INTEGER     ,INTENT(IN) :: LENI
      INTEGER(4)  ,INTENT(OUT):: IASCII(LENI)
      INTEGER(4)              :: LENG,I,J,K
!     *******************************************************************************
      LENG=LEN(STRING)
      IF(LENG*LENCH.NE.LENI) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZES')
        CALL ERROR$STOP('MPE_CHARTOASCII')
      END IF
      K=0
      DO I=1,LENCH
        DO J=1,LENG
          K=K+1
          IASCII(K)=ICHAR(STRING(I)(J:J))
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MPE_CHARTOASCII
!
!     .............................................................................
      SUBROUTINE MPE_CHARFROMASCII(LENCH,STRING,LENI,IASCII)
      IMPLICIT  NONE
      CHARACTER(*),INTENT(OUT):: STRING(LENCH)
      INTEGER     ,INTENT(IN) :: LENCH
      INTEGER     ,INTENT(IN) :: LENI
      INTEGER(4)  ,INTENT(IN) :: IASCII(LENI)
      INTEGER(4)              :: LENG,I,J,K
!     *******************************************************************************
      LENG=LEN(STRING)
      IF(LENG*LENCH.NE.LENI) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZES')
        CALL ERROR$STOP('MPE_CHARFROMASCII')
      END IF
      K=0
      DO I=1,LENCH
        DO J=1,LENG
          K=K+1
          STRING(I)(J:J)=CHAR(IASCII(K))
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MPE_CHARFROMASCII
#ENDIF
END MODULE MPE_MPIF_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE MPE$SELECT(CID)
!      *************************************************************************
!      ** SELECT A SPECIFIED COMMUNICATOR (GROUP OF PROCESSES)                **
!      *************************************************************************
       USE MPE_MPIF_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: CID    ! COMMUNICATOR ID
!      *************************************************************************
       IF(.NOT.TINI) THEN
         CALL ERROR$MSG('MPI$INIT MUST BE CALLED BEFORE ANY OTHER MPE ROUTINE')
         CALL ERROR$STOP('MPE$SELECT')
       END IF
       IF(THIS%ID.EQ.CID) RETURN
       THIS=>FIRSTTHIS
       DO 
         IF(THIS%ID.EQ.CID) THEN
           NTASKS=THIS%NTASKS
           THISTASK=THIS%THISTASK
           COMM=THIS%COMM
           SENDTAG=>THIS%SENDTAG
           RECEIVETAG=>THIS%RECEIVETAG
           EXIT
         ELSE
           IF(ASSOCIATED(THIS%NEXT)) THEN
             THIS=>THIS%NEXT
           ELSE
             CALL ERROR$MSG('UNKNOWN ID')
             CALL ERROR$CHVAL('ID',TRIM(CID))
             CALL ERROR$STOP('MPE$SELECT')
           END IF
         END IF
       ENDDO
       RETURN
       END SUBROUTINE MPE$SELECT
!
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_COMBINE_MODULE
PUBLIC MPE$COMBINE
INTERFACE MPE$COMBINE    !(OP,VAL)
#   MODULE TEMPLATE MPE$COMBINE
END INTERFACE 
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TEMPLATE MPE$COMBINE
(<TYPEID><TYPE><MPI_TYPE>)
                =([C8],[COMPLEX(8)][MY_MPITYPE_COMPLEX_KIND8])
                 ([R8],[REAL(8)][MY_MPITYPE_REAL_KIND8])
                 ([I4],[INTEGER(4)][MY_MPITYPE_INTEGER_KIND4])
(<RANKID><SIZE><RESHAPE(VAL)><RANK>)
                 =([R0],[1],[VAL],[])
                 =([R1],[SIZE(VAL)],[RESHAPE(VAL,(/LENG/))],[(:)])
                 =([R2],[SIZE(VAL)],[RESHAPE(VAL,(/LENG/))],[(:,:)])
                 =([R3],[SIZE(VAL)],[RESHAPE(VAL,(/LENG/))],[(:,:,:)])
                 =([R4],[SIZE(VAL)],[RESHAPE(VAL,(/LENG/))],[(:,:,:,:)])
                 =([R5],[SIZE(VAL)],[RESHAPE(VAL,(/LENG/))],[(:,:,:,:,:)])
#BODY
!
!     ..................................................................
      SUBROUTINE MPE$COMBINE<TYPEID><RANKID>(CID,OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE$COMBINER8                                         **
!     **                                                              **
!     **  PERFORMS A VECTOR OPERATION ACCROSS ALL TASKS               **
!     **    KEY(THISTASK)=KEY(TASK1).OP. ... .OP.KEY(NTASKS)          **
!     **  WHERE KEY(ITASK) IS A VECTOR, AND .OP. IS AN OPERATION      **
!     **  SUCH AS OP='+' OR OP='*'                                    **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    IS IT ALLOWED TO REPLACE MPI_DOUBLE BY MPI_REAL8,ETC??    **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: CID
      CHARACTER(*),INTENT(IN)    :: OPERATION
      <TYPE>      ,INTENT(INOUT) :: VAL<RANK>
      <TYPE>      ,ALLOCATABLE   :: WORK(:)
      INTEGER                    :: LENG
      INTEGER                    :: IERR
      LOGICAL     ,PARAMETER     :: TTEST=.FALSE.
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE$SELECT(CID)
!
!     =================================================================
!     == ALLOCATE WORK ARRAY                                         ==
!     =================================================================
      LENG=<SIZE>
      ALLOCATE(WORK(LENG))
                                 ! WORK(:)=<RESHAPE(VAL)>   
      WORK(:)=TRANSFER(VAL,WORK) ! RESHAPE FAILED FOR THE PATHSCALE COMPILER
!     
!     =================================================================
!     == PERFORM OPERATION                                           ==
!     =================================================================
      IF(OPERATION.EQ.'+')THEN
        CALL MPI_ALLREDUCE(WORK,VAL,LENG,<MPI_TYPE>,MPI_SUM,COMM,IERR)
      ELSE IF(OPERATION.EQ.'*')THEN
        CALL MPI_ALLREDUCE(WORK,VAL,LENG,<MPI_TYPE>,MPI_PROD,COMM,IERR)
      ELSE IF(OPERATION.EQ.'MAX')THEN
        CALL MPI_ALLREDUCE(WORK,VAL,LENG,<MPI_TYPE>,MPI_MAX,COMM,IERR)
      ELSE IF(OPERATION.EQ.'MIN')THEN
        CALL MPI_ALLREDUCE(WORK,VAL,LENG,<MPI_TYPE>,MPI_MIN,COMM,IERR)
      ELSE
        CALL ERROR$MSG('OPERATION NOT SUPPORTED')
        CALL ERROR$STOP('MPE$COMBINE')
      END IF
      IF(TTEST) THEN
        WORK=<RESHAPE(VAL)>
        CALL MPE_TESTCOMBINE<TYPEID><RANKID>(CID,LENG,WORK)
      END IF
      DEALLOCATE(WORK)
#ENDIF
      RETURN
      END SUBROUTINE MPE$COMBINE<TYPEID><RANKID>
#IFDEF CPPVARIABLE_PARALLEL
!
!     .......................................................................................
      SUBROUTINE MPE_TESTCOMBINE<TYPEID><RANKID>(CID,LEN,VAL1)
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: CID
      INTEGER(4),INTENT(IN)    :: LEN
      <TYPE>    ,INTENT(IN)    :: VAL1(LEN)
      <TYPE>    ,ALLOCATABLE   :: VALARR(:,:)
      INTEGER                  :: ITASK
      INTEGER(4)               :: ISVARARR(1),ISVAR
      LOGICAL(4)               :: TCHK
      INTEGER                  :: LENSTD
      INTEGER                  :: IERR
      INTEGER                  :: STAT(MPI_STATUS_SIZE)
      INTEGER                  :: FROMTASK0,TOTASK0,TAG
!     =========================================================================================
      LENSTD=LEN
      ALLOCATE(VALARR(LEN,NTASKS))
      VALARR(:,1)=VAL1(:)
      DO ITASK=2,NTASKS
        TAG=2000+ITASK
        IF(ITASK.EQ.THISTASK) THEN
          LENSTD=LEN
          TOTASK0=0
          CALL MPI_SEND(VAL1,LENSTD,<MPI_TYPE>,TOTASK0,TAG,COMM,IERR)
        END IF
        IF(THISTASK.EQ.1) THEN
          FROMTASK0=ITASK-1
          CALL MPI_RECV(VALARR(:,ITASK),LENSTD,<MPI_TYPE>,FROMTASK0,TAG,COMM,STAT,IERR)        
        END IF
      ENDDO
!     =================================================================
!     == ANALYSE                                                     ==
!     =================================================================
      TCHK=.FALSE.
      IF(THISTASK.EQ.1) THEN
        DO ITASK=2,NTASKS
          TCHK=TCHK.OR.(MAXVAL(ABS(VALARR(:,ITASK)-VALARR(:,1))).GT.0.D0)
        ENDDO
        IF(TCHK) THEN
          PRINT*,'ERROR DETECTED BY MPE_TESTCOMBINE: <TYPEID><RANKID>;',LEN
          DO ITASK=2,NTASKS
            ISVARARR=MAXLOC(ABS(VALARR(:,ITASK)-VALARR(:,1)))
            ISVAR=ISVARARR(1)
            PRINT*,ITASK,ISVAR,VALARR(ISVAR,ITASK)-VALARR(ISVAR,1),VALARR(ISVAR,ITASK),VALARR(ISVAR,1)
          ENDDO
          CALL ERROR$MSG('ERROR IN MPE$COMBINE; RESULT DIFFERS ON DIFFERENT TASKS')
          CALL ERROR$STOP('MPE_TESTCOMBINE')
        END IF
      END IF
      RETURN
      END SUBROUTINE MPE_TESTCOMBINE<TYPEID><RANKID>
#ENDIF
#END TEMPLATE MPE$COMBINE
END MODULE MPE_COMBINE_MODULE
!
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_BROADCAST_MODULE
PUBLIC MPE$BROADCAST
INTERFACE MPE$BROADCAST  !(FROMTASK,VAL)
#   MODULE TEMPLATE MPE$BROADCAST
    MODULE PROCEDURE MPE$BROADCASTCHA
    MODULE PROCEDURE MPE$BROADCASTCH
END INTERFACE 
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TEMPLATE MPE$BROADCAST
(<TYPEID><TYPE><MPI_TYPE>)
                =([C8],[COMPLEX(8)][MY_MPITYPE_COMPLEX_KIND8])
                 ([R8],[REAL(8)][MY_MPITYPE_REAL_KIND8])
                 ([R4],[REAL(4)][MY_MPITYPE_REAL_KIND4])
                 ([I4],[INTEGER(4)][MY_MPITYPE_INTEGER_KIND4])
                 ([L4],[LOGICAL(4)][MPI_LOGICAL])
(<RANKID><SIZE><RANK>)
                 =([R0],[1],[])
                 =([R1],[SIZE(VAL)],[(:)])
                 =([R2],[SIZE(VAL)],[(:,:)])
                 =([R3],[SIZE(VAL)],[(:,:,:)])
                 =([R4],[SIZE(VAL)],[(:,:,:,:)])
                 =([R5],[SIZE(VAL)],[(:,:,:,:,:)])
                 =([R6],[SIZE(VAL)],[(:,:,:,:,:,:)])
#BODY
!
!     ..................................................................
      SUBROUTINE MPE$BROADCAST<TYPEID><RANKID>(CID,FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE$BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: CID
      <TYPE>      ,INTENT(INOUT) :: VAL<RANK>
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER                    :: FROMTASK0
      INTEGER                    :: IERR
      INTEGER                    :: LENG
      CHARACTER(128)             :: ERRORSTRING
      INTEGER                    :: ERRORSTRINGLEN
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(BROADCAST_CLOCK)
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      LENG=<SIZE>
!PRINT*,THISTASK_WORLD,TRIM(THIS%ID),THISTASK,'BEFORE  MPE$BROADCAST<TYPEID><RANKID> ' &
!&     ,VAL,<SIZE>,<MPI_TYPE>,FROMTASK0,COMM
      CALL MPI_BCAST(VAL,LENG,<MPI_TYPE>,FROMTASK0,COMM,IERR)
!PRINT*,THISTASK_WORLD,TRIM(THIS%ID),THISTASK,'AFTER  MPE$BROADCAST<TYPEID><RANKID> ' &
!&     ,VAL,<SIZE>,<MPI_TYPE>,FROMTASK0,COMM,IERR
      IF(IERR.NE.0) THEN
        CALL MPI_ERROR_STRING(IERR,ERRORSTRING,ERRORSTRINGLEN)
        CALL ERROR$MSG('MPI ERROR')
        CALL ERROR$CHVAL('ERRORSTRING',ERRORSTRING)
        CALL ERROR$STOP('MPE$BROADCAST<TYPEID><RANKID>')
      END IF
      CALL MPE_CLOCKOFF(BROADCAST_CLOCK,LENG,'<TYPEID>')
#ENDIF
      RETURN
      END SUBROUTINE MPE$BROADCAST<TYPEID><RANKID>
#END TEMPLATE MPE$BROADCAST
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE$BROADCASTCHA(CID,FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE$BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: CID
      CHARACTER(*),INTENT(INOUT) :: VAL(:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER                    :: FROMTASK0
      INTEGER                    :: LENG
      INTEGER                    :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(BROADCAST_CLOCK)
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      LENG=LEN(VAL)*SIZE(VAL)
      CALL MPI_BCAST(VAL,LENG,MPI_CHARACTER,FROMTASK0 &
     &                ,COMM,IERR)
      CALL MPE_CLOCKOFF(BROADCAST_CLOCK,LENG,'CH')
#ENDIF
      RETURN
      END SUBROUTINE MPE$BROADCASTCHA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE$BROADCASTCH(CID,FROMTASK,VAL) 
!     ******************************************************************
!     **                                                              **
!     ** MPE$BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: CID
      CHARACTER(*),INTENT(INOUT) :: VAL
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER                    :: FROMTASK0
      INTEGER                    :: LENG
      INTEGER                    :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(BROADCAST_CLOCK)
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      LENG=LEN(VAL)
      CALL MPI_BCAST(VAL,LENG,MPI_CHARACTER,FROMTASK0,COMM,IERR)
      CALL MPE_CLOCKOFF(BROADCAST_CLOCK,LENG,'CH')
#ENDIF
      RETURN
      END SUBROUTINE MPE$BROADCASTCH
END MODULE MPE_BROADCAST_MODULE
!
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_SENDRECEIVE_MODULE
PUBLIC MPE$SENDRECEIVE
INTERFACE MPE$SENDRECEIVE       !(TOTASK,TAG,VAL)
#   MODULE TEMPLATE MPE$SENDRECEIVE
    MODULE PROCEDURE MPE$SENDRECEIVECHA
    MODULE PROCEDURE MPE$SENDRECEIVECH
END INTERFACE 
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TEMPLATE MPE$SENDRECEIVE
(<TYPEID><TYPE><MPI_TYPE>)
                =([C8],[COMPLEX(8)][MY_MPITYPE_COMPLEX_KIND8])
                 ([R8],[REAL(8)][MY_MPITYPE_REAL_KIND8])
                 ([R4],[REAL(4)][MY_MPITYPE_REAL_KIND4])
                 ([I4],[INTEGER(4)][MY_MPITYPE_INTEGER_KIND4])
                 ([L4],[LOGICAL(4)][MPI_LOGICAL])
(<RANKID><SIZE><RANK>)
                 =([R0],[1],[])
                 =([R1],[SIZE(VAL)],[(:)])
                 =([R2],[SIZE(VAL)],[(:,:)])
                 =([R3],[SIZE(VAL)],[(:,:,:)])
                 =([R4],[SIZE(VAL)],[(:,:,:,:)])
                 =([R5],[SIZE(VAL)],[(:,:,:,:,:)])
                 =([R6],[SIZE(VAL)],[(:,:,:,:,:,:)])
#BODY
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE$SENDRECEIVE<TYPEID><RANKID>(CID,FROMTASK,TOTASK,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  NAME: MPE$SENDRECEIVE                                               **
!     **                                                                      **
!     **  PURPOSE: SEND VAL FROM FROMTASK TO TOTASK                           **
!     **     VAL(TOTASK)=VAL(FROMTASK)                                        **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **    TAG IS SOME INTEGER NUMBER WHICH ALLOWES TO DIFFERENTIATE         **
!     **    BETWEEN DIFFERENT MESSAGES.                                       **
!     **                                                                      **
!     **************************************************************************  
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: CID
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)  ,INTENT(IN)    :: TOTASK
      <TYPE>      ,INTENT(INOUT) :: VAL<RANK>
      CHARACTER(32)              :: MSGID
      INTEGER                    :: TOTASK0
      INTEGER                    :: FROMTASK0
      INTEGER                    :: TAGSTD
      INTEGER                    :: LENG
      INTEGER                    :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER                 :: STAT(MPI_STATUS_SIZE)
!     **************************************************************************  
!     == NO DATA NEEDS TO BE SEND IF RECEIVER = SENDER
      IF(FROMTASK.EQ.TOTASK) RETURN
!     == SKIP ON NODES THAT DO NOT PARTICIPATE IN DATA EXCHANGE
      IF(THISTASK.NE.FROMTASK.AND.THISTASK.NE.TOTASK) RETURN
      CALL MPE_CLOCKON(SENDRECEIVE_CLOCK)
      CALL MPE$SELECT(CID)
      LENG=<SIZE>
      IF(THISTASK.EQ.FROMTASK) THEN
        CALL MPE_COUNTTAGS('SEND',TOTASK,TAGSTD)
        TOTASK0=TOTASK-1
        CALL MPI_SEND(VAL,LENG,<MPI_TYPE>,TOTASK0,TAGSTD,COMM,IERR)
      ELSE IF(THISTASK.EQ.TOTASK) THEN
        CALL MPE_COUNTTAGS('RECEIVE',FROMTASK,TAGSTD)
        FROMTASK0=FROMTASK-1
        CALL MPI_RECV(VAL,LENG,<MPI_TYPE>,FROMTASK0,TAGSTD,COMM,STAT,IERR)        
      END IF
      CALL MPE_CLOCKOFF(SENDRECEIVE_CLOCK,LENG,'<TYPEID>')
#ENDIF
      RETURN
      END SUBROUTINE MPE$SENDRECEIVE<TYPEID><RANKID>
#END TEMPLATE MPE$SENDRECEIVE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE$SENDRECEIVECHA(CID,FROMTASK,TOTASK,VAL)
!     **************************************************************************  
!     **                                                                      **
!     **  NAME: MPE$SENDRECEIVE                                               **    
!     **                                                                      **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                          **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE$RECEIVE                      **
!     **     VAL(TOTASK)=VAL(THISTASK)                                        **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **    TAG IS SOME INTEGER NUMBER WHICH ALLOWES TO DIFFERENTIATE         **
!     **    BETWEEN DIFFERENT MESSAGES.                                       **
!     **                                                                      **
!     **************************************************************************  
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: CID
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)  ,INTENT(IN)    :: TOTASK
      CHARACTER(*),INTENT(IN)    :: VAL(:)
      INTEGER                    :: FROMTASK0
      INTEGER                    :: TOTASK0
      INTEGER                    :: TAGSTD
      INTEGER                    :: LENG
      INTEGER                    :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER                 :: STAT(MPI_STATUS_SIZE)
!     **************************************************************************  
!     == NO DATA NEEDS TO BE SEND IF RECEIVER = SENDER
      IF(FROMTASK.EQ.TOTASK) RETURN
!     == SKIP ON NODES THAT DO NOT PARTICIPATE IN DATA EXCHANGE
      IF(THISTASK.NE.FROMTASK.AND.THISTASK.NE.TOTASK) RETURN
!
      CALL MPE_CLOCKON(SENDRECEIVE_CLOCK)
      CALL MPE$SELECT(CID)
      CALL MPE$QUERY(CID,NTASKS,THISTASK)
      LENG=LEN(VAL)*SIZE(VAL)
      IF(THISTASK.EQ.FROMTASK) THEN
        CALL MPE_COUNTTAGS('SEND',TOTASK,TAGSTD)
        TOTASK0=TOTASK-1
        CALL MPI_SEND(VAL,LENG,MPI_CHARACTER,TOTASK0,TAGSTD &
     &               ,COMM,IERR)
      ELSE IF(THISTASK.EQ.TOTASK) THEN
        CALL MPE_COUNTTAGS('RECEIVE',FROMTASK,TAGSTD)
        FROMTASK0=FROMTASK-1
        CALL MPI_RECV(VAL,LENG,MPI_CHARACTER,FROMTASK0,TAGSTD &
     &               ,COMM,STAT,IERR)        
      END IF
      CALL MPE_CLOCKOFF(SENDRECEIVE_CLOCK,LENG,'CH')
#ENDIF
      RETURN
      END SUBROUTINE MPE$SENDRECEIVECHA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE$SENDRECEIVECH(CID,FROMTASK,TOTASK,VAL)
!     **************************************************************************  
!     **                                                                      **
!     **  NAME: MPE$SENDRECEIVE                                               **
!     **                                                                      **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                          **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE$RECEIVE                      **
!     **     VAL(TOTASK)=VAL(THISTASK)                                        **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **    TAG IS SOME INTEGER NUMBER WHICH ALLOWES TO DIFFERENTIATE         **
!     **    BETWEEN DIFFERENT MESSAGES.                                       **
!     **                                                                      **
!     **************************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: CID
      INTEGER(4)  ,INTENT(IN) :: FROMTASK
      INTEGER(4)  ,INTENT(IN) :: TOTASK
      CHARACTER(*),INTENT(IN) :: VAL
      INTEGER                 :: FROMTASK0
      INTEGER                 :: TOTASK0
      INTEGER                 :: TAGSTD
      INTEGER                 :: LENG
      INTEGER                 :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER                 :: STAT(MPI_STATUS_SIZE)
!     **************************************************************************  
!     == NO DATA NEEDS TO BE SEND IF RECEIVER = SENDER
      IF(FROMTASK.EQ.TOTASK) RETURN
!     == SKIP ON NODES THAT DO NOT PARTICIPATE IN DATA EXCHANGE
      IF(THISTASK.NE.FROMTASK.AND.THISTASK.NE.TOTASK) RETURN
!
      CALL MPE_CLOCKON(SENDRECEIVE_CLOCK)
      CALL MPE$SELECT(CID)
      CALL MPE$QUERY(CID,NTASKS,THISTASK)
      LENG=LEN(VAL)
      IF(THISTASK.EQ.FROMTASK) THEN
        CALL MPE_COUNTTAGS('SEND',TOTASK,TAGSTD)
        TOTASK0=TOTASK-1
        CALL MPI_SEND(VAL,LENG,MPI_CHARACTER,TOTASK0,TAGSTD,COMM,IERR)
      ELSE IF(THISTASK.EQ.TOTASK) THEN
        CALL MPE_COUNTTAGS('RECEIVE',FROMTASK,TAGSTD)
        FROMTASK0=FROMTASK-1
        CALL MPI_RECV(VAL,LENG,MPI_CHARACTER,FROMTASK0,TAGSTD,COMM,STAT,IERR)        
      END IF
      CALL MPE_CLOCKOFF(SENDRECEIVE_CLOCK,LENG,'CH')
#ENDIF
      RETURN
      END SUBROUTINE MPE$SENDRECEIVECH
END MODULE MPE_SENDRECEIVE_MODULE
!     **************************************************************************  
!
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_SEND_MODULE
PUBLIC MPE$SEND
INTERFACE MPE$SEND       !(TOTASK,TAG,VAL)
#   MODULE TEMPLATE MPE$SEND
    MODULE PROCEDURE MPE$SENDCHA
    MODULE PROCEDURE MPE$SENDCH
END INTERFACE 
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TEMPLATE MPE$SEND
(<TYPEID><TYPE><MPI_TYPE>)
                =([C8],[COMPLEX(8)][MY_MPITYPE_COMPLEX_KIND8])
                 ([R8],[REAL(8)][MY_MPITYPE_REAL_KIND8])
                 ([R4],[REAL(4)][MY_MPITYPE_REAL_KIND4])
                 ([I4],[INTEGER(4)][MY_MPITYPE_INTEGER_KIND4])
                 ([L4],[LOGICAL(4)][MPI_LOGICAL])
(<RANKID><SIZE><RANK>)
                 =([R0],[1],[])
                 =([R1],[SIZE(VAL)],[(:)])
                 =([R2],[SIZE(VAL)],[(:,:)])
                 =([R3],[SIZE(VAL)],[(:,:,:)])
                 =([R4],[SIZE(VAL)],[(:,:,:,:)])
                 =([R5],[SIZE(VAL)],[(:,:,:,:,:)])
                 =([R6],[SIZE(VAL)],[(:,:,:,:,:,:)])
#BODY
!     ..................................................................
      SUBROUTINE MPE$SEND<TYPEID><RANKID>(CID,TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE$SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE$RECEIVE              **
!     **     VAL(TOTASK)=VAL(THISTASK)                                **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    TAG IS SOME INTEGER NUMBER WHICH ALLOWES TO DIFFERENTIATE **
!     **    BETWEEN DIFFERENT MESSAGES.                               **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: CID
      INTEGER(4)  ,INTENT(IN)    :: TOTASK
      INTEGER(4)  ,INTENT(IN)    :: TAG
      <TYPE>      ,INTENT(IN)    :: VAL<RANK>
      INTEGER                    :: TOTASK0
      INTEGER                    :: TAGSTD
      INTEGER                    :: LENG
      INTEGER                    :: IERR
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(SEND_CLOCK)
      CALL MPE$SELECT(CID)
      TOTASK0=TOTASK-1
      TAGSTD=TAG
!PRINT*,'MPE$SEND<TYPEID><RANKID> ',THISTASK,TOTASK,TAG
      LENG=<SIZE>
      CALL MPI_SEND(VAL,LENG,<MPI_TYPE>,TOTASK0,TAGSTD,COMM,IERR)
      CALL MPE_CLOCKOFF(SEND_CLOCK,LENG,'<TYPEID>')
#ELSE
      CALL ERROR$MSG('MPE$SEND MUST NOT BE CALLED IN SCALAR MODE')
      CALL ERROR$STOP('MPE$SEND<TYPEID><RANKID>')
#ENDIF
      RETURN
      END SUBROUTINE MPE$SEND<TYPEID><RANKID>
#END TEMPLATE MPE$SEND
!
!     ..................................................................
      SUBROUTINE MPE$SENDCHA(CID,TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE$SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE$RECEIVE              **
!     **     VAL(TOTASK)=VAL(THISTASK)                                **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    TAG IS SOME INTEGER NUMBER WHICH ALLOWES TO DIFFERENTIATE **
!     **    BETWEEN DIFFERENT MESSAGES.                               **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: CID
      INTEGER(4)  ,INTENT(IN)    :: TOTASK
      INTEGER(4)  ,INTENT(IN)    :: TAG
      CHARACTER(*),INTENT(IN)    :: VAL(:)
      INTEGER                    :: TOTASK0
      INTEGER                    :: TAGSTD
      INTEGER                    :: LENG
      INTEGER                    :: IERR
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(SEND_CLOCK)
      CALL MPE$SELECT(CID)
      TOTASK0=TOTASK-1
      TAGSTD=TAG
      LENG=LEN(VAL)*SIZE(VAL)
      CALL MPI_SEND(VAL,LENG,MPI_CHARACTER,TOTASK0,TAGSTD,COMM,IERR)
      CALL MPE_CLOCKOFF(SEND_CLOCK,LENG,'CH')
#ELSE
      CALL ERROR$MSG('MPE$SEND MUST NOT BE CALLED IN SCALAR MODE')
      CALL ERROR$STOP('MPE$SENDCH')
#ENDIF
      RETURN
      END SUBROUTINE MPE$SENDCHA
!
!     ..................................................................
      SUBROUTINE MPE$SENDCH(CID,TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE$SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE$RECEIVE              **
!     **     VAL(TOTASK)=VAL(THISTASK)                                **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    TAG IS SOME INTEGER NUMBER WHICH ALLOWES TO DIFFERENTIATE **
!     **    BETWEEN DIFFERENT MESSAGES.                               **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: CID
      INTEGER(4)  ,INTENT(IN) :: TOTASK
      INTEGER(4)  ,INTENT(IN) :: TAG
      CHARACTER(*),INTENT(IN) :: VAL
      INTEGER                 :: TOTASK0
      INTEGER                 :: TAGSTD
      INTEGER                 :: LENG
      INTEGER                 :: IERR
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(SEND_CLOCK)
      CALL MPE$SELECT(CID)
      TOTASK0=TOTASK-1
      TAGSTD=TAG
      LENG=LEN(VAL)
      CALL MPI_SEND(VAL,LENG,MPI_CHARACTER,TOTASK0,TAGSTD,COMM,IERR)
      CALL MPE_CLOCKOFF(SEND_CLOCK,LENG,'CH')
#ELSE
      CALL ERROR$MSG('MPE$SEND MUST NOT BE CALLED IN SCALAR MODE')
      CALL ERROR$STOP('MPE$SENDCH')
#ENDIF
      RETURN
      END SUBROUTINE MPE$SENDCH
END MODULE MPE_SEND_MODULE
!***********************************************************************

!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_RECEIVE_MODULE
PUBLIC MPE$RECEIVE
INTERFACE MPE$RECEIVE       !(FOMTASK,TAG,VAL)
#   MODULE TEMPLATE MPE$RECEIVE
    MODULE PROCEDURE MPE$RECEIVECHA
    MODULE PROCEDURE MPE$RECEIVECH
END INTERFACE 
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TEMPLATE MPE$RECEIVE
(<TYPEID><TYPE><MPI_TYPE><ZERO>)
                =([C8],[COMPLEX(8)][MY_MPITYPE_COMPLEX_KIND8][(0.d0,0.d0)])
                 ([R8],[REAL(8)][MY_MPITYPE_REAL_KIND8][0.d0])
                 ([R4],[REAL(4)][MY_MPITYPE_REAL_KIND4][0.])
                 ([I4],[INTEGER(4)][MY_MPITYPE_INTEGER_KIND4][0])
                 ([L4],[LOGICAL(4)][MPI_LOGICAL][.false.])
(<RANKID><SIZE><RANK>)
                 =([R0],[1],[])
                 =([R1],[SIZE(VAL)],[(:)])
                 =([R2],[SIZE(VAL)],[(:,:)])
                 =([R3],[SIZE(VAL)],[(:,:,:)])
                 =([R4],[SIZE(VAL)],[(:,:,:,:)])
                 =([R5],[SIZE(VAL)],[(:,:,:,:,:)])
                 =([R6],[SIZE(VAL)],[(:,:,:,:,:,:)])
#BODY
!
!     ..................................................................
      SUBROUTINE MPE$RECEIVE<TYPEID><RANKID>(CID,FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE$RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE$SEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: CID
      INTEGER(4)  ,INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4)  ,INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      <TYPE>      ,INTENT(OUT):: VAL<RANK>
      INTEGER                 :: FROMTASK0
      INTEGER                 :: TAGSTD
      INTEGER                 :: LENG
      INTEGER                 :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER                 :: STAT(MPI_STATUS_SIZE)
!     ******************************************************************
      CALL MPE_CLOCKON(RECEIVE_CLOCK)
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      TAGSTD=TAG
!PRINT*,'MPE$RECEIVE<TYPEID><RANKID> ',THISTASK,FROMTASK,TAG
      LENG=<SIZE>
      CALL MPI_RECV(VAL,LENG,<MPI_TYPE>,FROMTASK0,TAGSTD,COMM,STAT,IERR)        
      CALL MPE_CLOCKOFF(RECEIVE_CLOCK,LENG,'<TYPEID>')
#ELSE
      CALL ERROR$MSG('MPE$RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
      CALL ERROR$STOP('MPE$RECEIVER8A')
      val=<ZERO>    ! only to make the compiler happy
#ENDIF
      RETURN
      END SUBROUTINE MPE$RECEIVE<TYPEID><RANKID>
#END TEMPLATE MPE$RECEIVE
!
!     ..................................................................
      SUBROUTINE MPE$RECEIVECHA(CID,FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE$RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE$BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: CID
      INTEGER(4)  ,INTENT(IN)  :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4)  ,INTENT(IN)  :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      CHARACTER(*),INTENT(OUT) :: VAL(:)   
      INTEGER                  :: FROMTASK0
      INTEGER                  :: TAGSTD
      INTEGER                  :: LENG
      INTEGER                  :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER(4)               :: STAT(MPI_STATUS_SIZE)   
!     ******************************************************************
      CALL MPE_CLOCKON(RECEIVE_CLOCK)
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      TAGSTD=TAG
      LENG=LEN(VAL)*SIZE(VAL)
      CALL MPI_RECV(VAL,LENG,MPI_CHARACTER,FROMTASK0,TAGSTD &
     &               ,COMM,STAT,IERR)        
      CALL MPE_CLOCKOFF(RECEIVE_CLOCK,LENG,'CH')
#ELSE
      VAL(:)=' '
      CALL ERROR$MSG('MPE$RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
      CALL ERROR$STOP('MPE$RECEIVECHA')
#ENDIF
      RETURN
      END SUBROUTINE MPE$RECEIVECHA
!
!     ..................................................................
      SUBROUTINE MPE$RECEIVECH(CID,FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE$RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE$BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: CID
      INTEGER(4)  ,INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4)  ,INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      CHARACTER(*),INTENT(OUT):: VAL
      INTEGER               :: FROMTASK0
      INTEGER               :: TAGSTD
      INTEGER               :: LENG
      INTEGER               :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER(4)              :: STAT(MPI_STATUS_SIZE)
!     ******************************************************************
      CALL MPE_CLOCKON(RECEIVE_CLOCK)
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      TAGSTD=TAG
      LENG=LEN(VAL)
      CALL MPI_RECV(VAL,LENG,MPI_CHARACTER,FROMTASK0,TAGSTD,COMM,STAT,IERR)        
      CALL MPE_CLOCKOFF(RECEIVE_CLOCK,LENG,'CH')
#ELSE
      VAL=' '
      CALL ERROR$MSG('MPE$RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
      CALL ERROR$STOP('MPE$RECEIVECHA')
#ENDIF
      RETURN
      END SUBROUTINE MPE$RECEIVECH
END MODULE MPE_RECEIVE_MODULE
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_TRANSPOSE_MODULE
PUBLIC MPE$TRANSPOSE
INTERFACE MPE$TRANSPOSE   !(VAL)
#  MODULE TEMPLATE MPE$TRANSPOSE
END INTERFACE 
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TEMPLATE MPE$TRANSPOSE
(<TYPEID><TYPE><MPI_TYPE>)
                =([C8],[COMPLEX(8)][MY_MPITYPE_COMPLEX_KIND8])
                 ([R8],[REAL(8)][MY_MPITYPE_REAL_KIND8])
                 ([R4],[REAL(4)][MY_MPITYPE_REAL_KIND4])
                 ([I4],[INTEGER(4)][MY_MPITYPE_INTEGER_KIND4])
                 ([L4],[LOGICAL(4)][MPI_LOGICAL])
(<RANKID><SIZE><RANK>)
                 =([R1],[SIZE(VAL)],[(:)])
                 =([R2],[SIZE(VAL)],[(:,:)])
                 =([R3],[SIZE(VAL)],[(:,:,:)])
                 =([R4],[SIZE(VAL)],[(:,:,:,:)])
                 =([R5],[SIZE(VAL)],[(:,:,:,:,:)])
                 =([R6],[SIZE(VAL)],[(:,:,:,:,:,:)])
#BODY
!
!     ..................................................................
      SUBROUTINE MPE$TRANSPOSE<TYPEID><RANKID>(CID,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE$TRANSPOSE                                         **
!     **                                                              **
!     **  PURPOSE: THE JTH BLOCK OF DATA SENT FROM TASK I IS          **
!     **    RECEIVED BY THE JTH PROCESS AND PLACED INTO THE ITH BLOCK.**
!     **                                                              **
!     **    IN OTHER WORDS:                                           **
!     **    THINK OF VAL AS VECTOR OF ARRAYS VAL(:,ITASK)             **
!     **    CONTAINING NTASKS PORTIONS OF SIZE (SIZE(VAL)/NTASKS).    **
!     **    THEN THE VECTORS FROM ALL TASKS FORM A MATRIX             **
!     **    (VAL(:,ITASK,JTASK)).                                     **
!     **    THIS OPERATION TRANSPOSES THIS MATRIX                     **
!     **        VAL(:,ITASK,THISTASK)= VAL(:,THISTASK,ITASK)          **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: CID
      <TYPE>      ,INTENT(INOUT) :: VAL<RANK>
      INTEGER                    :: IERR
      INTEGER                    :: LENG,LENG1
      <TYPE>      ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(TRANSPOSE_CLOCK)
      CALL MPE$SELECT(CID)
      LENG=<SIZE>
      LENG1=LENG/NTASKS
      IF(LENG1*NTASKS.NE.LENG) THEN
        CALL ERROR$MSG('ARRAY SIZE NOT DIVIDABLE BY #(PROCESSORS)')
        CALL ERROR$STOP('MPE$TRANSPOSE<TYPEID><RANKID>')
      END IF
      ALLOCATE(VAL1(LENG))
                                 ! VAL1=RESHAPE(VAL,(/LENG/))
      VAL1(:)=TRANSFER(VAL,VAL1) ! RESHAPE FAILED FOR THE PATHSCALE COMPILER
      CALL MPI_ALLTOALL(VAL1,LENG1,<MPI_TYPE>,VAL,LENG1,<MPI_TYPE> &
     &                  ,COMM,IERR)
      DEALLOCATE(VAL1)
      CALL MPE_CLOCKOFF(TRANSPOSE_CLOCK,LENG1,'<TYPEID>')
#ENDIF
      RETURN
      END SUBROUTINE MPE$TRANSPOSE<TYPEID><RANKID>
#END TEMPLATE MPE$TRANSPOSE
END MODULE MPE_TRANSPOSE_MODULE
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_GATHER_MODULE
PUBLIC MPE$GATHER
INTERFACE MPE$GATHER    !(TOTASK, INVAL,OUTVAL,VAL)
#  MODULE TEMPLATE MPE$GATHER
   MODULE PROCEDURE MPE$GATHERCHA
   MODULE PROCEDURE MPE$GATHERCH
END INTERFACE 
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TEMPLATE MPE$GATHER
(<TYPEID><TYPE><MPI_TYPE>)
                =([C8],[COMPLEX(8)][MY_MPITYPE_COMPLEX_KIND8])
                 ([R8],[REAL(8)][MY_MPITYPE_REAL_KIND8])
                 ([R4],[REAL(4)][MY_MPITYPE_REAL_KIND4])
                 ([I4],[INTEGER(4)][MY_MPITYPE_INTEGER_KIND4])
                 ([L4],[LOGICAL(4)][MPI_LOGICAL])
(<RANKID><SIZE><RANK><RANK+1><RANK,1>)
                 =([R0],[1]          ,[]             ,[(:)]            ,[(1)])            
                 =([R1],[SIZE(INVAL)],[(:)]          ,[(:,:)]          ,[(:,1)])          
                 =([R2],[SIZE(INVAL)],[(:,:)]        ,[(:,:,:)]        ,[(:,:,1)])        
                 =([R3],[SIZE(INVAL)],[(:,:,:)]      ,[(:,:,:,:)]      ,[(:,:,:,1)])      
                 =([R4],[SIZE(INVAL)],[(:,:,:,:)]    ,[(:,:,:,:,:)]    ,[(:,:,:,:,1)])    
                 =([R5],[SIZE(INVAL)],[(:,:,:,:,:)]  ,[(:,:,:,:,:,:)]  ,[(:,:,:,:,:,1)])  
                 =([R6],[SIZE(INVAL)],[(:,:,:,:,:,:)],[(:,:,:,:,:,:,:)],[(:,:,:,:,:,:,1)])
#BODY
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE$GATHER<TYPEID><RANKID>(CID,TOTASK,INVAL,OUTVAL)
!     **************************************************************************
!     **                                                                      **
!     ** MPE$GATHER                                                           **
!     **                                                                      **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF              **
!     ** THE WORLD                                                            **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                                    **
!     **                                                                      **
!     **************************************************************************
      USE MPE_MPIF_MODULE 
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: CID
      INTEGER(4)  ,INTENT(IN) :: TOTASK
      <TYPE>      ,INTENT(IN) :: INVAL<RANK>
      <TYPE>      ,INTENT(OUT):: OUTVAL<RANK+1>
      INTEGER                 :: TOTASK0
      INTEGER                 :: LENG
      INTEGER                 :: IERR
!     **************************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(GATHER_CLOCK)
      CALL MPE$SELECT(CID)
      TOTASK0=TOTASK-1
      LENG=<SIZE>
      CALL MPI_GATHER(INVAL,LENG,<MPI_TYPE>,OUTVAL,LENG,<MPI_TYPE>,TOTASK0 &
     &                                                              ,COMM,IERR)
      CALL MPE_CLOCKOFF(GATHER_CLOCK,LENG,'<TYPEID>')
#ELSE
      OUTVAL<RANK,1>=INVAL<RANK>
#ENDIF
      RETURN
      END SUBROUTINE MPE$GATHER<TYPEID><RANKID>
#END TEMPLATE MPE$GATHER
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE$GATHERCHA(CID,TOTASK,INVAL,OUTVAL)
!     **************************************************************************
!     **                                                                      **
!     ** MPE$GATHER                                                           **
!     **                                                                      **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF              **
!     ** THE WORLD                                                            **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                                    **
!     **                                                                      **
!     **************************************************************************
      USE MPE_MPIF_MODULE 
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN) :: CID
      INTEGER(4)      ,INTENT(IN) :: TOTASK
      CHARACTER(LEN=*),INTENT(IN) :: INVAL(:)
      CHARACTER(LEN=*),INTENT(OUT):: OUTVAL(:,:)
      INTEGER                     :: TOTASK0
      INTEGER                     :: LENG
      INTEGER                     :: IERR
!     **************************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(GATHER_CLOCK)
      CALL MPE$SELECT(CID)
      TOTASK0=TOTASK-1
      LENG=LEN(INVAL)*SIZE(INVAL)
      CALL MPI_GATHER(INVAL,LENG,MPI_CHARACTER,OUTVAL,LENG,MPI_CHARACTER &
     &                 ,TOTASK0,COMM,IERR)
      CALL MPE_CLOCKOFF(GATHER_CLOCK,LENG,'CH')
#ELSE
      OUTVAL(:,1)=INVAL(:)
#ENDIF
      RETURN
      END SUBROUTINE MPE$GATHERCHA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE$GATHERCH(CID,TOTASK,INVAL,OUTVAL)
!     **************************************************************************
!     **                                                                      **
!     ** MPE$GATHER                                                           **
!     **                                                                      **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF              **
!     ** THE WORLD                                                            **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                                    **
!     **                                                                      **
!     **************************************************************************
      USE MPE_MPIF_MODULE 
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN) :: CID
      INTEGER(4)      ,INTENT(IN) :: TOTASK
      CHARACTER(LEN=*),INTENT(IN) :: INVAL
      CHARACTER(LEN=*),INTENT(OUT):: OUTVAL(:)
      INTEGER                     :: TOTASK0
      INTEGER                     :: LENG
      INTEGER                     :: IERR
!     **************************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(GATHER_CLOCK)
      CALL MPE$SELECT(CID)
      TOTASK0=TOTASK-1
      LENG=LEN(INVAL)
      CALL MPI_GATHER(INVAL,LENG,MPI_CHARACTER,OUTVAL,LENG,MPI_CHARACTER &
     &                 ,TOTASK0,COMM,IERR)
      CALL MPE_CLOCKOFF(GATHER_CLOCK,LENG,'CH')
#ELSE
      OUTVAL(1)=INVAL
#ENDIF
      RETURN
      END SUBROUTINE MPE$GATHERCH
END MODULE MPE_GATHER_MODULE
!
!.......................................................................
MODULE MPE_MODULE
!=======================================================================
!==  THE FILE MPIF90.H IS THE INCLUDE FILE MPIF.H PROVIDED WITH THE 
!==  MESSAGE PASSING INTERFACE (MPI), WHICH IS CONVERTED 
!==  TO FORTRAN90 FREE FORMAT
!=======================================================================
USE MPE_GATHER_MODULE
USE MPE_TRANSPOSE_MODULE
USE MPE_SENDRECEIVE_MODULE
USE MPE_SEND_MODULE
USE MPE_RECEIVE_MODULE
USE MPE_COMBINE_MODULE
USE MPE_BROADCAST_MODULE
END MODULE MPE_MODULE
!    
!     ..................................................................
      SUBROUTINE MPE$INIT
!     ******************************************************************
!     **                                                              **
!     **  INITIALIZE THIS MODULE AND MPI                              **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **  THE COMMONBLOCK EXECUTABLEDESCRIPTION IS USED ONLY          **
!     **  TO BE ABLE ABLE TO DISTINGUISH A SCALAR EXECUTABLE          **
!     **  FROM A PARALLEL EXECUTABLE BECAUSE VALUE OF PARALLELVERSION **
!     **  WILL APPEAR LITERALLY IN THE EXECUTABLE                     **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER                  :: IERR
!      INTEGER     ,PARAMETER   :: MBYTE=2**20
!      INTEGER     ,PARAMETER   :: BUFFER_SIZE=6*MBYTE
!      CHARACTER(1),POINTER    :: BUFFER(:)
!      INTEGER(4),POINTER       :: BUFFER(:)
      REAL(4)                  :: XREAL4
      REAL(8)                  :: XREAL8
      COMPLEX(8)               :: XCOMPLEX8
      INTEGER(4)               :: XINTEGER4
      INTEGER                  :: SIZE
!     ******************************************************************
      IF(TINI) THEN
        CALL ERROR$MSG('MPE$INIT MUST ONLY BE CALLED ONCE')
        CALL ERROR$STOP('MPE$INIT')
      END IF     
      TINI=.TRUE.
!
!     ==================================================================
!     == GET NUMBER OF TASKS , ID OF THIS ONE                         ==
!     == IN THIS ROUTINE , TASKIDS RANGE FROM 0 TO NTASKNUM -1,       ==
!     == WHILE FOR A CALL FROM OUTSIDE, NTASKID S ARE TO BE           ==
!     == OUT OF [1,NTASKNUM]                                          ==
!     ==================================================================
#IFDEF CPPVARIABLE_PARALLEL
!     ================================================================
!     == INITIALIZE MPI                                             ==
!     ================================================================
      CALL MPI_INIT(IERR)
!     
!     ================================================================
!     ==  PROVIDE A BUFFER TO BUFFER MESSAGES                       ==
!     ==  SENT WITH MPI_NSEND AND MPI_IBSEND                        ==
!     ================================================================
!      SIZE=1+INT(REAL(8*BUFFER_SIZE)/REAL(BIT_SIZE(BUFFER)))
!      ALLOCATE(BUFFER(SIZE))
!      CALL MPI_BUFFER_ATTACH(BUFFER,BUFFER_SIZE,IERR)
!     
!     ================================================================
!     ==  DETERMINE #(TASKS) AND NUMBER OF THIS TASK                ==
!     ================================================================
      COMM=MPI_COMM_WORLD
      CALL MPI_COMM_SIZE(COMM,NTASKS,IERR)
      CALL MPI_COMM_RANK(COMM,THISTASK,IERR)
      THISTASK= THISTASK+1
      THISTASK_WORLD=THISTASK
      ALLOCATE(FIRSTTHIS)
      THIS=>FIRSTTHIS
      THIS%ID='~'
      THIS%THISTASK=THISTASK
      THIS%NTASKS=NTASKS
      THIS%COMM=COMM
      ALLOCATE(THIS%SENDTAG(NTASKS))
      ALLOCATE(THIS%RECEIVETAG(NTASKS))
      THIS%SENDTAG(:)=0
      THIS%RECEIVETAG(:)=0
      NULLIFY(THIS%NEXT)
      CALL MPE$SELECT('~')
      CALL MPI_PCONTROL(1,IERR)
#IFDEF CPPVAR_MPI2
      CALL MPI_SIZEOF(INTEGER4,SIZE,IERR)
      CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER,SIZE,MY_MPITYPE_INTEGER_KIND4,IERR)
      CALL MPI_SIZEOF(XREAL4,SIZE,IERR)
      CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,SIZE,MY_MPITYPE_REAL_KIND4,IERR)
      CALL MPI_SIZEOF(XREAL8,SIZE,IERR)
      CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,SIZE,MY_MPITYPE_REAL_KIND8,IERR)
      CALL MPI_SIZEOF(XCOMPLEX8,SIZE,IERR)
      CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_COMPLEX,SIZE,MY_MPITYPE_COMPLEX_KIND8,IERR)
#ENDIF
#ENDIF
!     
!     ================================================================
!     ==  RESET CLOCKS                                              ==
!     ================================================================
      BROADCAST_CLOCK%ID='BROADCAST'
      SENDRECEIVE_CLOCK%ID='SENDRECEIVE'
      SEND_CLOCK%ID='SEND'
      RECEIVE_CLOCK%ID='RECEIVE'
      GATHER_CLOCK%ID='GATHER'
      COMBINE_CLOCK%ID='COMBINE'
      TRANSPOSE_CLOCK%ID='TRANSPOSE'
      SYNC_CLOCK%ID='SYNC'
      CALL MPE_CLOCKRESET(BROADCAST_CLOCK)
      CALL MPE_CLOCKRESET(SENDRECEIVE_CLOCK)
      CALL MPE_CLOCKRESET(SEND_CLOCK)
      CALL MPE_CLOCKRESET(RECEIVE_CLOCK)
      CALL MPE_CLOCKRESET(GATHER_CLOCK)
      CALL MPE_CLOCKRESET(COMBINE_CLOCK)
      CALL MPE_CLOCKRESET(TRANSPOSE_CLOCK)
      CALL MPE_CLOCKRESET(SYNC_CLOCK)
      RETURN
      END
!    
!     ..................................................................
      SUBROUTINE MPE$NEW(CID,NEWCID,NTASKS_,ICOLOR_)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: CID
      CHARACTER(*),INTENT(IN)  :: NEWCID
      INTEGER     ,INTENT(IN)  :: NTASKS_
      INTEGER(4)  ,INTENT(IN)  :: ICOLOR_(NTASKS)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: I
      INTEGER(4)               :: NEWCOMM
      INTEGER(4)               :: ICOLOR(NTASKS_)
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE$SELECT(CID)
      IF(NTASKS_.NE.NTASKS) THEN
        CALL ERROR$MSG('INCONSISTENT SIZES')
        CALL ERROR$STOP('MPE$NEW')
      END IF
      ICOLOR(:)=ICOLOR_(:)
!
!     ===================================================================
!     == THE BROADCAST COMMAND IS TO ENSURE THAT ICOLOR IS IDENTICAL   ==
!     == ON ALL TASKS AND TO ENSURE THAT MPE$NEW IS EXECUTED ON ALL    ==
!     == NODES OF THE CORRESPONDING GROUP                              ==
!     ===================================================================
      CALL MPE$BROADCAST(CID,1,ICOLOR)
      ALLOCATE(THIS%NEXT)
      THIS=>THIS%NEXT
      THIS%ID=NEWCID
      NULLIFY(THIS%NEXT)
!PRINT*,'MPE$NEW',THIS%ID,THISTASK,ICOLOR
      CALL MPI_COMM_SPLIT(COMM,ICOLOR(THISTASK),THISTASK-1,NEWCOMM,IERR)
      COMM=NEWCOMM
      THIS%COMM=COMM
      CALL MPI_COMM_SIZE(COMM,NTASKS,IERR)
      THIS%NTASKS=NTASKS
      CALL MPI_COMM_RANK(COMM,THISTASK,IERR)
      THISTASK=THISTASK+1
!PRINT*,'MPE$NEW',THIS%ID,THISTASK
      THIS%THISTASK=THISTASK
      ALLOCATE(THIS%SENDTAG(NTASKS))
      ALLOCATE(THIS%RECEIVETAG(NTASKS))
      THIS%SENDTAG(:)=0
      THIS%RECEIVETAG(:)=0
      CALL MPE$SELECT(CID)
#ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MPE$EXIT 
!     ******************************************************************
!     ** MPE$EXIT                                                     **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER           :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPI_FINALIZE(IERR)
#ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MPE$STOPALL(NCODE_)
!     ******************************************************************
!     ** MPE$STOPALL                                                  **
!     ** STOP ALL TASKS IN PARTITION WITH EXITCODE NCODE_             **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCODE_
      INTEGER               :: IERR,IERROR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE$SELECT('~')
      WRITE(*,*) 'EXIT:',NCODE_
      CALL MPI_ABORT(COMM,IERROR,7777)
      STOP 'ERROR STOP'
#ELSE
      WRITE(*,*) 'EXIT:',NCODE_
      STOP 'ERROR STOP'
#ENDIF
      END
!
!     ..................................................................
      SUBROUTINE MPE$QUERY(CID,NTASKS_,THISTASK_)
!     ******************************************************************
!     ** MPE$QUERY                                                    **
!     == RETURN THE TOTAL NUMBER OF TASKS                             ==
!     == AND THE ID OF THIS TASK [1 ..NTASKS]                         ==
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE 
      CHARACTER(*),INTENT(IN)  :: CID   ! COMMUNICATOR ID
      INTEGER(4)  ,INTENT(OUT) :: NTASKS_   !#(TASKS IN THE GROUP)
      INTEGER(4)  ,INTENT(OUT) :: THISTASK_ ! NUMBER OF THE TASK IN THE GROUP
      INTEGER                  :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE$SELECT(CID)
      NTASKS_  =NTASKS
      THISTASK_=THISTASK
#ELSE
      NTASKS_  = 1
      THISTASK_= 1
#ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MPE$SYNC(CID)
!     ******************************************************************
!     **                                                              **
!     **  MPE$SYNC(CID)                                               **
!     **                                                              **
!     **  PURPOSE: SYNCHRONIZE ALL TASKS                              **
!     **  (WAITS UNTILL ALL TASKS ARRIVED AT THIS POINT OF THE CODE)  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: CID   ! COMMUNICATOR ID
      INTEGER                  :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE_CLOCKON(SYNC_CLOCK)
      CALL MPE$SELECT(CID)
      CALL MPI_BARRIER(COMM,IERR)
      CALL MPE_CLOCKOFF(SYNC_CLOCK,0,'I4')
#ENDIF
      RETURN
      END
!!$!
!!$!     ..................................................................
!!$      SUBROUTINE MPE$ISPARALLEL(IPARALLEL)
!!$!     ******************************************************************
!!$!     ** MPE$ISPARALLEL                                               **
!!$!     ******************************************************************
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(OUT) :: IPARALLEL      
!!$!     ******************************************************************
!!$#IFDEF CPPVARIABLE_PARALLEL
!!$      IPARALLEL=1
!!$#ELSE
!!$      IPARALLEL=0
!!$#ENDIF
!!$      RETURN
!!$      END
!!$!
!!$!     ..................................................................
!!$      SUBROUTINE MPE$INDEX(CID,INVAL,OUTVAL,SIZE)
!!$!     ******************************************************************
!!$!     **                                                              **
!!$!     ** MPE$INDEX                                                    **
!!$!     **                                                              **
!!$!     **                                                              **
!!$!     **                                                              **
!!$!     ** 2 TASKS AOUT(TASK1) := (1,8), AOUT(TASK2) := (4,9)           **
!!$!     **         ARES(TASK1)  = (1,4), AOUT(TASK2)  = (8,9)           **
!!$!     **                                                              **
!!$!     ******************************************************************
!!$      USE MPE_MPIF_MODULE
!!$      IMPLICIT NONE
!!$      CHARACTER(*),INTENT(IN)  :: CID   ! COMMUNICATOR ID
!!$      INTEGER(4)  ,INTENT(IN)  :: SIZE
!!$      CHARACTER(1),INTENT(IN)  :: INVAL(SIZE)
!!$      CHARACTER(1),INTENT(OUT) :: OUTVAL(SIZE)
!!$      INTEGER                  :: SIZESTD
!!$      INTEGER                  :: IERR
!!$!     ******************************************************************
!!$      CALL ERROR$MSG('ROUTINE IS MARKED FOR DELETION')
!!$      CALL ERROR$STOP('MPE$INDEX')
!!$#IFDEF CPPVARIABLE_PARALLEL
!!$      CALL MPE$SELECT(CID)
!!$      SIZESTD=SIZE
!!$      CALL MPI_ALLTOALL(INVAL,SIZESTD,MPI_BYTE,OUTVAL,SIZESTD,MPI_BYTE,COMM,IERR)
!!$#ELSE
!!$      OUTVAL(:)=INVAL(:)
!!$#ENDIF
!!$      RETURN
!!$      END
! 
!     .....................................................................
      SUBROUTINE MPE_COUNTTAGS(ID,OTHERTASK,TAG)
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4),INTENT(IN) :: OTHERTASK
      INTEGER   ,INTENT(OUT) :: TAG
      INTEGER(4)             :: TOTASK
      INTEGER(4)             :: FROMTASK
      INTEGER(4),PARAMETER   :: MAXTAG=10000
!     *********************************************************************
      IF(ID.EQ.'SEND') THEN
        TOTASK=OTHERTASK
        FROMTASK=THISTASK
        SENDTAG(TOTASK)=SENDTAG(TOTASK)+1
        IF(SENDTAG(TOTASK).GT.MAXTAG) SENDTAG(TOTASK)=1
        TAG=FROMTASK+NTASKS*(TOTASK-1+NTASKS*SENDTAG(TOTASK))
!PRINT*,'SENDTAG ',THISTASK,OTHERTASK,'::',SENDTAG,'--',TAG
      ELSE IF(ID.EQ.'RECEIVE') THEN
        FROMTASK=OTHERTASK
        TOTASK=THISTASK
        RECEIVETAG(FROMTASK)=RECEIVETAG(FROMTASK)+1
        IF(RECEIVETAG(FROMTASK).GT.MAXTAG) RECEIVETAG(FROMTASK)=1
        TAG=FROMTASK+NTASKS*(TOTASK-1+NTASKS*RECEIVETAG(FROMTASK))
!PRINT*,'RECEIVETAG ',THISTASK,OTHERTASK,'::',RECEIVETAG,'--',TAG
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED') 
        CALL ERROR$MSG('ALLOWED VALUES ARE SEND AND RECEIVE')
        CALL ERROR$STOP('MPE_COUNTTAGS')
      END IF
      RETURN
      END SUBROUTINE MPE_COUNTTAGS
! 
!     .....................................................................
      SUBROUTINE MPE$REPORT(NFIL)
      USE MPE_MODULE
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      TYPE(THISTYPE),POINTER :: CURRENT
      INTEGER(4)             :: NTASKS_WORLD
      INTEGER(4)             :: THISTASK_WORLD1
      INTEGER(4)             :: NTASKS_CURRENT
      INTEGER(4)             :: THISTASK_CURRENT
      INTEGER(4),ALLOCATABLE :: ICOLOR(:)
      INTEGER(4)             :: I,ISVAR
      CHARACTER(128)         :: CID
!     *********************************************************************
#IFDEF CPPVARIABLE_PARALLEL
                          CALL TRACE$PUSH('MPE$REPORT')
      CALL MPE$QUERY('~',NTASKS_WORLD,THISTASK_WORLD1)
      ALLOCATE(ICOLOR(NTASKS_WORLD))      
      DO I=1,NTASKS_WORLD
        ICOLOR(I)=I
      ENDDO
      CALL REPORT$TITLE(NFIL,'MPE-REPORT')
      WRITE(NFIL,FMT='(A10,20I4/T20,20I4)')'CID\TASKID',ICOLOR
      CURRENT=>FIRSTTHIS
      DO 
        CID=CURRENT%ID
        CALL MPE$QUERY(CID,NTASKS_CURRENT,THISTASK_CURRENT)
        ICOLOR(:)=0
        ICOLOR(THISTASK_WORLD)=THISTASK_WORLD
        CALL MPE$BROADCAST(CID,1,ICOLOR(THISTASK_WORLD))
        CALL MPE$COMBINE('~','+',ICOLOR)
        WRITE(NFIL,FMT='(A10,20I4/T20,20I4)')CID,ICOLOR
        IF(.NOT.ASSOCIATED(CURRENT%NEXT)) EXIT 
        CURRENT=>CURRENT%NEXT
      ENDDO
      DEALLOCATE(ICOLOR)
                                  CALL TRACE$POP()
#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE_CLOCKON(CLOCK)
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      TYPE(MPECLOCKTYPE),INTENT(INOUT):: CLOCK
      REAL(8)                         :: USRTIME
      REAL(8)                         :: SYSTIME
      REAL(8)                         :: CPUTIME
!     **************************************************************************
      CALL LIB$ETIME(USRTIME,SYSTIME)
      CPUTIME=USRTIME+SYSTIME
      CLOCK%STARTTIME=CPUTIME
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE_CLOCKOFF(CLOCK,dataSIZE,TYPEID)
      USE MPE_MPIF_MODULE, ONLY : MPECLOCKTYPE
      IMPLICIT NONE
      TYPE(MPECLOCKTYPE),INTENT(INOUT):: CLOCK
      integer           ,INTENT(IN)   :: dataSIZE                    
      CHARACTER(*)      ,INTENT(IN)   :: TYPEID
      REAL(8)                         :: SIZE1
      REAL(8)                         :: USRTIME
      REAL(8)                         :: SYSTIME
      REAL(8)                         :: CPUTIME
!     **************************************************************************
      IF(TYPEID.EQ.'R8') THEN
        SIZE1=real(DATASIZE*8,kind=8)
      ELSE IF(TYPEID.EQ.'I4') THEN
        SIZE1=real(DATASIZE*4,kind=8)
      ELSE IF(TYPEID.EQ.'C8') THEN
        SIZE1=real(DATASIZE*16,kind=8)
      ELSE IF(TYPEID.EQ.'L4') THEN
        SIZE1=real(DATASIZE*4,kind=8)
      ELSE IF(TYPEID.EQ.'R4') THEN
        SIZE1=real(DATASIZE*4,kind=8)
      ELSE IF(TYPEID.EQ.'C4') THEN
        SIZE1=real(DATASIZE*8,kind=8)
      ELSE IF(TYPEID.EQ.'CH') THEN
        SIZE1=real(DATASIZE,kind=8)
      ELSE
        CALL ERROR$MSG('TYPEID NOT RECOGNIZED')
        CALL ERROR$CHVAL('TYPEID',TYPEID)
        CALL ERROR$STOP('MPE_CLOCKOFF')
      END IF
      CLOCK%COUNT     =CLOCK%COUNT+1
      CLOCK%SIZE      =CLOCK%SIZE+SIZE1
      CLOCK%SIZESQUARE=CLOCK%SIZESQUARE+SIZE1**2
      CLOCK%MAXSIZE   =MAX(CLOCK%MAXSIZE,SIZE1)
      CLOCK%MINSIZE   =MIN(CLOCK%MINSIZE,SIZE1)
      CALL LIB$ETIME(USRTIME,SYSTIME)
      CPUTIME         =USRTIME+SYSTIME-CLOCK%STARTTIME
      CLOCK%TIME      =CLOCK%TIME+CPUTIME
      CLOCK%TIMESQUARE=CLOCK%TIMESQUARE+CPUTIME**2
      CLOCK%MAXTIME   =MAX(CLOCK%MAXTIME,CPUTIME)
      CLOCK%MINTIME   =MIN(CLOCK%MINTIME,CPUTIME)
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE$CLOCKREPORT(NFIL)
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4)        ,INTENT(IN)   :: NFIL
      INTEGER(4)               :: COUNT
      REAL(8)                  :: TTIME
      REAL(8)                  :: AVSIZE
      REAL(8)                  :: BANDWIDTH
      REAL(8)                  :: LATENCY   
!     **************************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      call mpe$query('~',ntasks,thistask)
      if(thistask.eq.1) then
        WRITE(NFIL,FMT='("TIMING REPORT OF MPE (MY PARALLEL ENVIRONMENT)")')
        WRITE(NFIL,FMT='("==============================================")')
        WRITE(NFIL,FMT='("ID",T17,"COUNT",T26,"T[S]",T39,"<SIZE/BYTE>",T52,"X(SIZE/BYTE)")')
      end if
      CALL MPE_CLOCKANALYZE(NFIL,BROADCAST_CLOCK)
      CALL MPE_CLOCKANALYZE(NFIL,SENDRECEIVE_CLOCK)
      CALL MPE_CLOCKANALYZE(NFIL,SEND_CLOCK)
      CALL MPE_CLOCKANALYZE(NFIL,RECEIVE_CLOCK)
      CALL MPE_CLOCKANALYZE(NFIL,GATHER_CLOCK)
      CALL MPE_CLOCKANALYZE(NFIL,COMBINE_CLOCK)
      CALL MPE_CLOCKANALYZE(NFIL,TRANSPOSE_CLOCK)
      CALL MPE_CLOCKANALYZE(NFIL,SYNC_CLOCK)
#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE_CLOCKANALYZE(nfil,CLOCK)
      USE MPE_MODULE
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4)                   :: nfil
      TYPE(MPECLOCKTYPE),INTENT(IN):: CLOCK
      INTEGER(4)                   :: COUNT
      REAL(8)                      :: TTIME
      REAL(8)                      :: AVSIZE
      REAL(8)                      :: BANDWIDTH
      REAL(8)                      :: LATENCY
      real(8)                      :: arr(3)
!      INTEGER     ,PARAMETER   :: MBYTE=2**20
!     **************************************************************************
      call mpe$query('~',ntasks,thistask)
      COUNT=CLOCK%COUNT
      CALL MPE$COMBINE('~','+',COUNT)
      ARR(1)=CLOCK%TIME
      ARR(2)=CLOCK%SIZE
      ARR(3)=CLOCK%MAXSIZE
      CALL MPE$COMBINE('~','+',ARR(1:2))
      CALL MPE$COMBINE('~','MAX',ARR(3))
      IF(COUNT.EQ.0)RETURN
      if(thistask.ne.1) return
      ARR(2)=ARR(2)/REAL(COUNT,KIND=8)
!!$      BANDWIDTH=SQRT((CLOCK%SIZESQUARE-CLOCK%SIZE**2) &
!!$     &              /(CLOCK%TIMESQUARE-CLOCK%TIME**2))
!!$      LATENCY  =(BANDWIDTH*CLOCK%SIZE-CLOCK%TIME)/REAL(COUNT,KIND=8)
      WRITE(NFIL,FMT='(a,T17,I9,6ES12.3)') &
     &            clock%id,clock%COUNT,arr(1:3)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MPE_CLOCKRESET(CLOCK)
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      TYPE(MPECLOCKTYPE),INTENT(OUT):: CLOCK
!     **************************************************************************
      CLOCK%COUNT      =0
      CLOCK%SIZE       =0.D0
      CLOCK%SIZESQUARE =0.D0
      CLOCK%MINSIZE    =0.D0
      CLOCK%MAXSIZE    =0.D0
      CLOCK%TIME       =0.D0
      CLOCK%TIMESQUARE =0.D0
      CLOCK%MINTIME    =0.D0
      CLOCK%MAXTIME    =0.D0
      CLOCK%STARTTIME  =0.D0
      RETURN
      END

