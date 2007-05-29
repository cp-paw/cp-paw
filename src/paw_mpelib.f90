@PROCESS NOEXTCHK 
!***********************************************************************
!**                                                                   **
!**  NAME: MPELIB                                                     **
!**                                                                   **
!**  PURPOSE: PERFORM COMMUNICATIONS AMONG PROCESSORS FOR             **
!**   PARALLEL PROCESSING.                                            **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    MPE$INIT                                                       **
!**    MPE$QUERY(CID,NTASKS,THISTASK)                                 **
!**    MPE$SYNC(CID)                                                  **
!**    MPE$STOPALL(RETURNCODE)                                        **
!**    MPE$EXIT                                                       **
!**    MPE$ISPARALLEL(IPARALLEL) ??                                   **
!**    MPE$COMBINE(CID,OPERATION,VAL)                                 **
!**    MPE$BROADCAST(CID,FROMTASK,VAL)                                **
!**    MPE$SEND(CID,TOTASK,TAG,VAL)                                   **
!**    MPE$RECEIVE(CID,FROMTASK,TAG,VAL)                              **
!**    MPE$TRANSPOSE(CID,VAL)                                         **
!**    MPE$GATHER(CID,TOTASK,INVAL,OUTVAL)                            **
!**                                                                   **
!**  REMARKS:                                                         **
!**    1) INTERFACE TO MPI (MESSAGE PASSING LIBRARY)                  **
!**    2) THE C-PREPROCESSOR IS USED TO CREATE BOTH THE SINGLE        **
!**       AND THE MULTIPROCESSOR VERSION OUT OF THIS FILE.            **
!**       THUS THIS FILE MUST NOT BE MADE UPPERCASE                   **
!**    3) THIS MODULE REQUIRES THE FILE MPIF90.H IN THE PATH FOR      **
!**       INCLUDE FILES                                               **
!**    4) CID STANDS FOR COMMUNICATOR ID AND DETERMINES THE GROUP OF  **
!**       PROCESSORS  ON WHICH THE COMMAND IS EXECUTED                **
!**                                                                   **
!**                                          ERNST NUSTERER 1994/1996 **
!**                                MODIFIED: PETER E. BLOECHL, (1996) **
!***********************************************************************
!
!=========================================================================
!==includefile for mpi                                                  ==
!=========================================================================
MODULE MPI
#IFDEF CPPVARIABLE_PARALLEL
!  include file is uppercase on purpose. 
   INCLUDE 'MPIF.H'    
#ENDIF
END MODULE MPI
!
MODULE MPE_MPIF_MODULE
#IFDEF CPPVARIABLE_PARALLEL
  use mpi
#ENDIF
  TYPE THISTYPE
    CHARACTER(128)     :: ID
    INTEGER(4)         :: THISTASK
    INTEGER(4)         :: NTASKS
    INTEGER            :: COMM
    INTEGER(4),POINTER :: SENDTAG(:)
    INTEGER(4),POINTER :: RECEIVETAG(:)
    TYPE(THISTYPE),POINTER:: NEXT
  END TYPE THISTYPE
  LOGICAL(4)      :: TINI=.FALSE.
  INTEGER         :: COMM
  INTEGER         :: NTASKS
  INTEGER         :: THISTASK
  INTEGER,POINTER :: SENDTAG(:)
  INTEGER,POINTER :: RECEIVETAG(:)
  TYPE(THISTYPE),POINTER :: FIRSTTHIS
  TYPE(THISTYPE),POINTER :: this
#IFDEF CPPVARIABLE_PARALLEL
!__DATATYPES CONNECTED WITH KIND PARAMETERS
!__FOR MPI1: USE INTEGER(4)=INTEGER*4, REAL(8)=INTEGER*8, COMPLEX(8)=COMPLEX*16
!__FOR MPI2: RECALCULATE VALUES IN MPE$INIT
  INTEGER,SAVE    :: MY_MPITYPE_INTEGER_KIND4=MPI_INTEGER4
  INTEGER,SAVE    :: MY_MPITYPE_REAL_KIND4=MPI_real4
  INTEGER,SAVE    :: MY_MPITYPE_REAL_KIND8=MPI_REAL8
  INTEGER,SAVE    :: MY_MPITYPE_COMPLEX_KIND8=MPI_double_complex
#ENDIF
END MODULE MPE_MPIF_MODULE
!      ..................................................................
       SUBROUTINE MPE$SELECT(CID)
       USE MPE_MPIF_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: CID
!      ******************************************************************
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
           exit
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
      character(128)             :: errorstring
      integer                    :: errorstringlen
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
!print*,'before  MPE$BROADCAST<TYPEID><RANKID> ',' <MPI_TYPE> ',<MPI_TYPE>,<SIZE>
      CALL MPI_BCAST(VAL,<SIZE>,<MPI_TYPE>,FROMTASK0,COMM,IERR)
      IF(IERR.NE.0) THEN
        CALL MPI_ERROR_STRING(IERR,ERRORSTRING,ERRORSTRINGLEN)
        CALL ERROR$MSG('MPI ERROR')
        CALL ERROR$CHVAL('ERRORSTRING',ERRORSTRING)
        CALL ERROR$STOP('MPE$BROADCAST<TYPEID><RANKID>')
      END IF
#ENDIF
      RETURN
      END SUBROUTINE MPE$BROADCAST<TYPEID><RANKID>
#END TEMPLATE MPE$BROADCAST
!
!     ..................................................................
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
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      LENG=LEN(VAL)*SIZE(VAL)
      CALL MPI_BCAST(VAL,LENG,MPI_CHARACTER,FROMTASK0 &
     &                ,COMM,IERR)
#ENDIF
      RETURN
      END SUBROUTINE MPE$BROADCASTCHA
!
!     ..................................................................
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
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      LENG=LEN(VAL)
      CALL MPI_BCAST(VAL,LENG,MPI_CHARACTER,FROMTASK0 &
     &                ,COMM,IERR)
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
!     ..................................................................
      SUBROUTINE MPE$SENDRECEIVE<TYPEID><RANKID>(CID,FROMTASK,TOTASK,VAL,msgid_)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE$SENDRECEIVE                                       **
!     **                                                              **
!     **  PURPOSE: SEND VAL FROM FROMTASK TO TOTASK                   **
!     **     VAL(TOTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    TAG IS SOME INTEGER NUMBER WHICH ALLOWES TO DIFFERENTIATE **
!     **    BETWEEN DIFFERENT MESSAGES.                               **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: CID
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)  ,INTENT(IN)    :: TOTASK
      <TYPE>      ,INTENT(INOUT) :: VAL<RANK>
      character(*),intent(in),optional :: msgid_
      character(32)              :: msgid
      INTEGER                    :: TOTASK0
      INTEGER                    :: FROMTASK0
      INTEGER                    :: TAGSTD
      INTEGER                    :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER                 :: STAT(MPI_STATUS_SIZE)
!     ******************************************************************
      CALL MPE$SELECT(CID)
      IF(PRESENT(MSGID_)) THEN
        MSGID=MSGID_
        CALL MPE$SENDRECEIVECH(CID,FROMTASK,TOTASK,MSGID)
        IF(trim(MSGID).NE.MSGID_) THEN
          CALL ERROR$I4VAL('THISTASK',THISTASK)
          CALL ERROR$I4VAL('FROMTASK',FROMTASK)
          CALL ERROR$I4VAL('TOTASK',TOTASK)
          CALL ERROR$CHVAL('MSGID_',MSGID_)
          CALL ERROR$CHVAL('MSGID',trim(MSGID))
          CALL ERROR$STOP('MPE$SENDRECEIVE<TYPEID><RANKID>')
        END IF
      END IF
      IF(FROMTASK.EQ.TOTASK) RETURN
      IF(THISTASK.EQ.FROMTASK) THEN
        CALL MPE_COUNTTAGS('SEND',TOTASK,TAGSTD)
        TOTASK0=TOTASK-1
!PRINT*,'SEND FROM ',FROMTASK,' TO ',TOTASK,' MSG: ',TAGSTD
        CALL MPI_SEND(VAL,<SIZE>,<MPI_TYPE>,TOTASK0,TAGSTD &
     &               ,COMM,IERR)
      ELSE IF(THISTASK.EQ.TOTASK) THEN
        CALL MPE_COUNTTAGS('RECEIVE',FROMTASK,TAGSTD)
        FROMTASK0=FROMTASK-1
!PRINT*,'RECEIVE FROM ',FROMTASK,' ON ',TOTASK,' MSG: ',TAGSTD
        CALL MPI_RECV(VAL,<SIZE>,<MPI_TYPE>,FROMTASK0,TAGSTD &
     &               ,COMM,STAT,IERR)        
      END IF
#ENDIF
      RETURN
      END SUBROUTINE MPE$SENDRECEIVE<TYPEID><RANKID>
#END TEMPLATE MPE$SENDRECEIVE
!
!     ..................................................................
      SUBROUTINE MPE$SENDRECEIVECHA(CID,FROMTASK,TOTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE$SENDRECEIVE                                              **
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
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)  ,INTENT(IN)    :: TOTASK
      CHARACTER(*),INTENT(IN)    :: VAL(:)
      INTEGER                    :: FROMTASK0
      INTEGER                    :: TOTASK0
      INTEGER                    :: TAGSTD
      INTEGER                    :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER                 :: STAT(MPI_STATUS_SIZE)
!     ******************************************************************
      CALL MPE$SELECT(CID)
      CALL MPE$QUERY(CID,NTASKS,THISTASK)
      IF(THISTASK.EQ.FROMTASK) THEN
        CALL MPE_COUNTTAGS('SEND',TOTASK,TAGSTD)
        TOTASK0=TOTASK-1
        CALL MPI_SEND(VAL,LEN(VAL)*SIZE(VAL),MPI_CHARACTER,TOTASK0,TAGSTD &
     &               ,COMM,IERR)
      ELSE IF(THISTASK.EQ.TOTASK) THEN
        CALL MPE_COUNTTAGS('RECEIVE',FROMTASK,TAGSTD)
        FROMTASK0=FROMTASK-1
        CALL MPI_RECV(VAL,LEN(VAL)*SIZE(VAL),MPI_CHARACTER,FROMTASK0,TAGSTD &
     &               ,COMM,STAT,IERR)        
      END IF
#ENDIF
      RETURN
      END SUBROUTINE MPE$SENDRECEIVECHA
!
!     ..................................................................
      SUBROUTINE MPE$SENDRECEIVECH(CID,FROMTASK,TOTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE$SENDRECEIVE                                              **
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
      INTEGER(4)  ,INTENT(IN) :: FROMTASK
      INTEGER(4)  ,INTENT(IN) :: TOTASK
      CHARACTER(*),INTENT(IN) :: VAL
      INTEGER                 :: FROMTASK0
      INTEGER                 :: TOTASK0
      INTEGER                 :: TAGSTD
      INTEGER                 :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER                 :: STAT(MPI_STATUS_SIZE)
!     ******************************************************************
      CALL MPE$SELECT(CID)
      CALL MPE$QUERY(CID,NTASKS,THISTASK)
      IF(THISTASK.EQ.FROMTASK) THEN
        CALL MPE_COUNTTAGS('SEND',TOTASK,TAGSTD)
        TOTASK0=TOTASK-1
        CALL MPI_SEND(VAL,LEN(VAL),MPI_CHARACTER,TOTASK0,TAGSTD &
     &               ,COMM,IERR)
      ELSE IF(THISTASK.EQ.TOTASK) THEN
        CALL MPE_COUNTTAGS('RECEIVE',FROMTASK,TAGSTD)
        FROMTASK0=FROMTASK-1
        CALL MPI_RECV(VAL,LEN(VAL),MPI_CHARACTER,FROMTASK0,TAGSTD &
     &               ,COMM,STAT,IERR)        
      END IF
#ENDIF
      RETURN
      END SUBROUTINE MPE$SENDRECEIVECH
END MODULE MPE_SENDRECEIVE_MODULE
!***********************************************************************
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
      INTEGER                    :: IERR
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
integer(4):: nfiltrace
CALL FILEHANDLER$UNIT('TRACE',NFILTRACE)
      CALL MPE$SELECT(CID)
      TOTASK0=TOTASK-1
      TAGSTD=TAG
write(nfiltrace,*)'before send',totask0,tagstd
       CALL MPI_SEND(VAL,<SIZE>,<MPI_TYPE>,TOTASK0,TAGSTD,COMM,IERR)
write(nfiltrace,*)'after send'
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
      INTEGER                    :: IERR
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
      CALL MPE$SELECT(CID)
      TOTASK0=TOTASK-1
      TAGSTD=TAG
      CALL MPI_SEND(VAL,LEN(VAL)*SIZE(VAL),MPI_CHARACTER,TOTASK0,TAGSTD &
     &               ,COMM,IERR)
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
      INTEGER                 :: IERR
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
      CALL MPE$SELECT(CID)
      TOTASK0=TOTASK-1
      TAGSTD=TAG
      CALL MPI_SEND(VAL,LEN(VAL),MPI_CHARACTER,TOTASK0,TAGSTD &
     &               ,COMM,IERR)
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
      SUBROUTINE MPE$RECEIVE<TYPEID><RANKID>(CID,FROMTASK,TAG,VAL)
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
      <TYPE>      ,INTENT(OUT):: VAL<RANK>
      INTEGER                 :: FROMTASK0
      INTEGER                 :: TAGSTD
      INTEGER                 :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER                 :: STAT(MPI_STATUS_SIZE)
!     ******************************************************************
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      TAGSTD=TAG
      CALL MPI_RECV(VAL,<SIZE>,<MPI_TYPE>,FROMTASK0,TAGSTD &
     &               ,COMM,STAT,IERR)        
#ELSE
      CALL ERROR$MSG('MPE$RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
      CALL ERROR$STOP('MPE$RECEIVER8A')
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
      INTEGER               :: FROMTASK0
      INTEGER               :: TAGSTD
      INTEGER               :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER(4)               :: STAT(MPI_STATUS_SIZE)   
!     ******************************************************************
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      TAGSTD=TAG
      CALL MPI_RECV(VAL,LEN(VAL)*SIZE(VAL),MPI_CHARACTER,FROMTASK0,TAGSTD &
     &               ,COMM,STAT,IERR)        
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
      INTEGER               :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER(4)              :: STAT(MPI_STATUS_SIZE)
!     ******************************************************************
      CALL MPE$SELECT(CID)
      FROMTASK0=FROMTASK-1
      TAGSTD=TAG
      CALL MPI_RECV(VAL,LEN(VAL),MPI_CHARACTER,FROMTASK0,TAGSTD &
     &               ,COMM,STAT,IERR)        
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
        CALL MPE$SYNC(CID)
        CALL MPI_ALLTOALL(VAL1,LENG1,<MPI_TYPE>,VAL,LENG1,<MPI_TYPE> &
     &                  ,COMM,IERR)
        DEALLOCATE(VAL1)
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
!     ..................................................................
      SUBROUTINE MPE$GATHER<TYPEID><RANKID>(CID,TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE$GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE 
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: CID
      INTEGER(4)  ,INTENT(IN) :: TOTASK
      <TYPE>      ,INTENT(IN) :: INVAL<RANK>
      <TYPE>      ,INTENT(OUT):: OUTVAL<RANK+1>
      INTEGER                 :: TOTASK0
      INTEGER                 :: LENG
      INTEGER                 :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE$SELECT(CID)
      TOTASK0=TOTASK-1
      LENG=<SIZE>
      CALL MPI_GATHER(INVAL,LENG,<MPI_TYPE>,OUTVAL,LENG,<MPI_TYPE> &
     &                 ,TOTASK0,COMM,IERR)
#ELSE
      OUTVAL<RANK,1>=INVAL<RANK>
#ENDIF
      RETURN
      END SUBROUTINE MPE$GATHER<TYPEID><RANKID>
#END TEMPLATE MPE$GATHER
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
      INTEGER     ,PARAMETER   :: MBYTE=2**20
      INTEGER     ,PARAMETER   :: BUFFER_SIZE=6*MBYTE
      CHARACTER(1),POINTER     :: BUFFER(:)
      real(4)                  :: xreal4
      real(8)                  :: xreal8
      complex(8)               :: xcomplex8
      integer(4)               :: xinteger4
      integer                  :: size
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
      ALLOCATE(BUFFER(BUFFER_SIZE))
      CALL MPI_BUFFER_ATTACH(BUFFER,BUFFER_SIZE,IERR)
!     
!     ================================================================
!     ==  DETERMINE #(TASKS) AND NUMBER OF THIS TASK                ==
!     ================================================================
      comm=mpi_comm_world
      CALL MPI_COMM_SIZE(COMM,NTASKS,IERR)
      CALL MPI_COMM_RANK(COMM,THISTASK,IERR)
      THISTASK= THISTASK+1
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
      call mpi_pcontrol(1,ierr)
#IFDEF CPPVAR_MPI2
      CALL MPI_SIZEOF(INTEGER4,SIZE,IERR)
      CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER,SIZE,My_MPITYPE_INTEGER_kind4,IERR)
      CALL MPI_SIZEOF(XREAL4,SIZE,IERR)
      CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,SIZE,MY_MPITYPE_REAL_kind4,IERR)
      CALL MPI_SIZEOF(XREAL8,SIZE,IERR)
      CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,SIZE,MY_MPITYPE_REAL_kind8,IERR)
      CALL MPI_SIZEOF(XCOMPLEX8,SIZE,IERR)
      CALL MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_COMPLEX,SIZE,MY_MPITYPE_COMPLEX_KIND8,IERR)
#ENDIF
#ENDIF
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
      INTEGER     ,INTENT(IN)  :: ICOLOR_(NTASKS)
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
!     == the broadcast command is to ensure that icolor is identical   ==
!     == on all tasks and to ensure that mpe$new is executed on all    ==
!     == nodes of the corresponding group                              ==
!     ===================================================================
      CALL MPE$BROADCAST(CID,1,ICOLOR)
      ALLOCATE(THIS%NEXT)
      THIS=>THIS%NEXT
      THIS%ID=NEWCID
      NULLIFY(THIS%NEXT)
      CALL MPI_COMM_SPLIT(COMM,ICOLOR(THISTASK),THISTASK-1,NEWCOMM,IERR)
      COMM=NEWCOMM
      THIS%COMM=COMM
      CALL MPI_COMM_SIZE(COMM,NTASKS,IERR)
      THIS%NTASKS=ntasks
      CALL MPI_COMM_RANK(COMM,THISTASK,IERR)
      thistask=thistask+1
      THIS%THISTASK=thistask
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
      CALL MPE$SELECT(CID)
      CALL MPI_BARRIER(COMM,IERR)
#ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MPE$ISPARALLEL(IPARALLEL)
!     ******************************************************************
!     ** MPE$ISPARALLEL                                               **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: IPARALLEL      
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      IPARALLEL=1
#ELSE
      IPARALLEL=0
#ENDIF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MPE$INDEX(CID,INVAL,OUTVAL,SIZE)
!     ******************************************************************
!     **                                                              **
!     ** MPE$INDEX                                                    **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     ** 2 TASKS AOUT(TASK1) := (1,8), AOUT(TASK2) := (4,9)           **
!     **         ARES(TASK1)  = (1,4), AOUT(TASK2)  = (8,9)           **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: CID   ! COMMUNICATOR ID
      INTEGER(4)  ,INTENT(IN)  :: SIZE
      CHARACTER(1),INTENT(IN)  :: INVAL(SIZE)
      CHARACTER(1),INTENT(OUT) :: OUTVAL(SIZE)
      INTEGER                  :: SIZESTD
      INTEGER                  :: IERR
!     ******************************************************************
      CALL ERROR$MSG('ROUTINE IS MARKED FOR DELETION')
      CALL ERROR$STOP('MPE$INDEX')
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE$SELECT(CID)
      SIZESTD=SIZE
      CALL MPI_ALLTOALL(INVAL,SIZESTD,MPI_BYTE,OUTVAL,SIZESTD &
     &                  ,MPI_BYTE,COMM,IERR)
#ELSE
      OUTVAL(:)=INVAL(:)
#ENDIF
      RETURN
      END
! 
!     .....................................................................
      SUBROUTINE MPE_COUNTTAGS(ID,OTHERTASK,TAG)
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4),INTENT(IN) :: OTHERTASK
      INTEGER(4),INTENT(OUT) :: TAG
      INTEGER(4)             :: TOTASK
      INTEGER(4)             :: FROMTASK
      INTEGER(4),PARAMETER   :: MAXTAG=100000
!     *********************************************************************
      IF(ID.EQ.'SEND') THEN
        TOTASK=OTHERTASK
        FROMTASK=THISTASK
        SENDTAG(TOTASK)=SENDTAG(TOTASK)+1
        IF(SENDTAG(TOTASK).GT.MAXTAG) SENDTAG(TOTASK)=1
        TAG=FROMTASK+NTASKS*(TOTASK-1+NTASKS*SENDTAG(TOTASK))
!print*,'sendtag ',thistask,othertask,'::',sendtag,'--',tag
      ELSE IF(ID.EQ.'RECEIVE') THEN
        FROMTASK=OTHERTASK
        TOTASK=THISTASK
        RECEIVETAG(FROMTASK)=RECEIVETAG(FROMTASK)+1
        IF(RECEIVETAG(FROMTASK).GT.MAXTAG) RECEIVETAG(FROMTASK)=1
        TAG=FROMTASK+NTASKS*(TOTASK-1+NTASKS*RECEIVETAG(FROMTASK))
!print*,'receivetag ',thistask,othertask,'::',receivetag,'--',tag
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
      TYPE(THISTYPE),pointer :: CURRENT
      INTEGER(4)             :: NTASKS_WORLD
      INTEGER(4)             :: THISTASK_WORLD
      INTEGER(4)             :: NTASKS_CURRENT
      INTEGER(4)             :: THISTASK_CURRENT
      INTEGER(4),ALLOCATABLE :: ICOLOR(:)
      INTEGER(4)             :: I
      CHARACTER(128)         :: CID
!     *********************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPE$QUERY('~',NTASKS_WORLD,THISTASK_WORLD)
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
#ENDIF
      RETURN
      END
