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
!**    MPE$QUERY(NTASKS,THISTASK)                                     **
!**    MPE$SYNC                                                       **
!**    MPE$STOPALL(RETURNCODE)                                        **
!**    MPE$EXIT                                                       **
!**    MPE$ISPARALLEL(IPARALLEL) ??                                   **
!**    MPE$COMBINE(OPERATION,VAL)                                     **
!**    MPE$BROADCAST(FROMTASK,VAL)                                    **
!**    MPE$SEND(TOTASK,TAG,VAL)                                       **
!**    MPE$RECEIVE(FROMTASK,TAG,VAL)                                  **
!**    MPE$TRANSPOSE(VAL)                                             **
!**    MPE$GATHER(TOTASK,INVAL,OUTVAL)                                **
!**                                                                   **
!**  REMARKS:                                                         **
!**    1) INTERFACE TO MPI (MESSAGE PASSING LIBRARY)                  **
!**    2) THE C-PREPROCESSOR IS USED TO CREATE BOTH THE SINGLE        **
!**       AND THE MULTIPROCESSOR VERSION OUT OF THIS FILE.            **
!**       THUS THIS FILE MUST NOT BE MADE UPPERCASE                   **
!**    3) THIS MODULE REQUIRES THE FILE MPIF90.H IN THE PATH FOR      **
!**       INCLUDE FILES                                               **
!**                                                                   **
!**           ERNST NUSTERER 1994/1996                                **
!** MODIFIED: PETER E. BLOECHL, IBM ZURICH RESEARCH LABORATORY (1996) **
!***********************************************************************
!
MODULE MPE_MPIF_MODULE
#IFDEF CPPVARIABLE_PARALLEL
  INCLUDE 'mpif.h'
#ENDIF
  INTEGER(4) :: NTASKS
  INTEGER(4) :: THISTASK
END MODULE MPE_MPIF_MODULE
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
                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
                 ([R8],[REAL(8)][MPI_REAL8])
                 ([I4],[INTEGER(4)][MPI_INTEGER4])
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
      SUBROUTINE MPE$COMBINE<TYPEID><RANKID>(OPERATION,VAL)
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
      CHARACTER(*),INTENT(IN)    :: OPERATION
      <TYPE>      ,INTENT(INOUT) :: VAL<RANK>
      <TYPE>      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
!
!     =================================================================
!     == ALLOCATE WORK ARRAY                                         ==
!     =================================================================
      LENG=<SIZE>
      ALLOCATE(WORK(LENG))
                                 ! WORK(:)=<RESHAPE(VAL)>   
      work(:)=transfer(val,work) ! reshape failed for the pathscale compiler
!     
!     =================================================================
!     == PERFORM OPERATION                                           ==
!     =================================================================
      IF(OPERATION.EQ.'+')THEN
        CALL MPI_ALLREDUCE(WORK,VAL,LENG &
     &           ,<MPI_TYPE>,MPI_SUM,MPI_COMM_WORLD,IERR)
      ELSE IF(OPERATION.EQ.'*')THEN
        CALL MPI_ALLREDUCE(WORK,VAL,LENG &
     &           ,<MPI_TYPE>,MPI_PROD,MPI_COMM_WORLD,IERR)
      ELSE IF(OPERATION.EQ.'MAX')THEN
        CALL MPI_ALLREDUCE(WORK,VAL,LENG &
     &           ,<MPI_TYPE>,MPI_MAX,MPI_COMM_WORLD,IERR)
      ELSE IF(OPERATION.EQ.'MIN')THEN
        CALL MPI_ALLREDUCE(WORK,VAL,LENG &
     &           ,<MPI_TYPE>,MPI_MIN,MPI_COMM_WORLD,IERR)
      ELSE
        CALL ERROR$MSG('OPERATION NOT SUPPORTED')
        CALL ERROR$STOP('MPE$COMBINE')
      END IF
      DEALLOCATE(WORK)
#ENDIF
      RETURN
      END SUBROUTINE MPE$COMBINE<TYPEID><RANKID>
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
                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
                 ([R8],[REAL(8)][MPI_REAL8])
                 ([R4],[REAL(4)][MPI_REAL4])
                 ([I4],[INTEGER(4)][MPI_INTEGER4])
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
      SUBROUTINE MPE$BROADCAST<TYPEID><RANKID>(FROMTASK,VAL)
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
      <TYPE>      ,INTENT(INOUT) :: VAL<RANK>
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPI_BCAST(VAL,<SIZE>,<MPI_TYPE>,FROMTASK-1 &
     &                ,MPI_COMM_WORLD,IERR)
#ENDIF
      RETURN
      END SUBROUTINE MPE$BROADCAST<TYPEID><RANKID>
#END TEMPLATE MPE$BROADCAST
!
!     ..................................................................
      SUBROUTINE MPE$BROADCASTCHA(FROMTASK,VAL)
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
      CHARACTER(*),INTENT(INOUT) :: VAL(:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
        LENG=LEN(VAL)*SIZE(VAL)
        CALL MPI_BCAST(VAL,LENG,MPI_CHARACTER,FROMTASK-1 &
     &                ,MPI_COMM_WORLD,IERR)
#ENDIF
      RETURN
      END SUBROUTINE MPE$BROADCASTCHA
!
!     ..................................................................
      SUBROUTINE MPE$BROADCASTCH(FROMTASK,VAL) 
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
      CHARACTER(*),INTENT(INOUT) :: VAL
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
        LENG=LEN(VAL)
        CALL MPI_BCAST(VAL,LENG,MPI_CHARACTER,FROMTASK-1 &
     &                ,MPI_COMM_WORLD,IERR)
#ENDIF
      RETURN
      END SUBROUTINE MPE$BROADCASTCH
END MODULE MPE_BROADCAST_MODULE
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
                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
                 ([R8],[REAL(8)][MPI_REAL8])
                 ([R4],[REAL(4)][MPI_REAL4])
                 ([I4],[INTEGER(4)][MPI_INTEGER4])
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
      SUBROUTINE MPE$SEND<TYPEID><RANKID>(TOTASK,TAG,VAL)
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
      INTEGER(4)  ,INTENT(IN)    :: TOTASK
      INTEGER(4)  ,INTENT(IN)    :: TAG
      <TYPE>      ,INTENT(IN)    :: VAL<RANK>
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
        CALL MPI_SEND(VAL,<SIZE>,<MPI_TYPE>,TOTASK-1,TAG &
     &               ,MPI_COMM_WORLD,IERR)
#ELSE
        CALL ERROR$MSG('MPE$SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR$STOP('MPE$SEND<TYPEID><RANKID>')
#ENDIF
      RETURN
      END SUBROUTINE MPE$SEND<TYPEID><RANKID>
#END TEMPLATE MPE$SEND
!
!     ..................................................................
      SUBROUTINE MPE$SENDCHA(TOTASK,TAG,VAL)
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
      INTEGER(4)  ,INTENT(IN)    :: TOTASK
      INTEGER(4)  ,INTENT(IN)    :: TAG
      CHARACTER(*),INTENT(IN)    :: VAL(:)
      INTEGER(4)                 :: IERR
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
        CALL MPI_SEND(VAL,LEN(VAL)*SIZE(VAL),MPI_CHARACTER,TOTASK-1,TAG &
     &               ,MPI_COMM_WORLD,IERR)
#ELSE
        CALL ERROR$MSG('MPE$SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR$STOP('MPE$SENDCH')
#ENDIF
      RETURN
      END SUBROUTINE MPE$SENDCHA
!
!     ..................................................................
      SUBROUTINE MPE$SENDCH(TOTASK,TAG,VAL)
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
      INTEGER(4)  ,INTENT(IN) :: TOTASK
      INTEGER(4)  ,INTENT(IN) :: TAG
      CHARACTER(*),INTENT(IN) :: VAL
      INTEGER(4)              :: IERR
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
      CALL MPI_SEND(VAL,LEN(VAL),MPI_CHARACTER,TOTASK-1,TAG &
     &               ,MPI_COMM_WORLD,IERR)
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
                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
                 ([R8],[REAL(8)][MPI_REAL8])
                 ([R4],[REAL(4)][MPI_REAL4])
                 ([I4],[INTEGER(4)][MPI_INTEGER4])
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
      SUBROUTINE MPE$RECEIVE<TYPEID><RANKID>(FROMTASK,TAG,VAL)
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
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      <TYPE>    ,INTENT(OUT):: VAL<RANK>
      INTEGER(4)            :: IERR
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER(4)            :: STAT(MPI_STATUS_SIZE)
!     ******************************************************************
      CALL MPI_RECV(VAL,<SIZE>,<MPI_TYPE>,FROMTASK-1,TAG &
     &               ,MPI_COMM_WORLD,STAT,IERR)        
#ELSE
      CALL ERROR$MSG('MPE$RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
      CALL ERROR$STOP('MPE$RECEIVER8A')
#ENDIF
      RETURN
      END SUBROUTINE MPE$RECEIVE<TYPEID><RANKID>
#END TEMPLATE MPE$RECEIVE
!
!     ..................................................................
      SUBROUTINE MPE$RECEIVECHA(FROMTASK,TAG,VAL)
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
      INTEGER(4)  ,INTENT(IN)  :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4)  ,INTENT(IN)  :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      CHARACTER(*),INTENT(OUT) :: VAL(:)   
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER(4)               :: STAT(MPI_STATUS_SIZE)   
      INTEGER(4)               :: IERR
!     ******************************************************************
      CALL MPI_RECV(VAL,LEN(VAL)*SIZE(VAL),MPI_CHARACTER,FROMTASK-1,TAG &
     &               ,MPI_COMM_WORLD,STAT,IERR)        
#ELSE
      VAL(:)=' '
      CALL ERROR$MSG('MPE$RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
      CALL ERROR$STOP('MPE$RECEIVECHA')
#ENDIF
      RETURN
      END SUBROUTINE MPE$RECEIVECHA
!
!     ..................................................................
      SUBROUTINE MPE$RECEIVECH(FROMTASK,TAG,VAL)
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
      INTEGER(4)  ,INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4)  ,INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      CHARACTER(*),INTENT(OUT):: VAL
#IFDEF  CPPVARIABLE_PARALLEL
      INTEGER(4)              :: STAT(MPI_STATUS_SIZE)
      INTEGER(4)              :: IERR
!     ******************************************************************
      CALL MPI_RECV(VAL,LEN(VAL),MPI_CHARACTER,FROMTASK-1,TAG &
     &               ,MPI_COMM_WORLD,STAT,IERR)        
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
                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
                 ([R8],[REAL(8)][MPI_REAL8])
                 ([R4],[REAL(4)][MPI_REAL4])
                 ([I4],[INTEGER(4)][MPI_INTEGER4])
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
      SUBROUTINE MPE$TRANSPOSE<TYPEID><RANKID>(VAL)
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
      <TYPE>    ,INTENT(INOUT) :: VAL<RANK>
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
!      INTEGER(4)               :: NTASKS,THISTASK
      <TYPE>    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
#IFDEF  CPPVARIABLE_PARALLEL
        LENG=<SIZE>
!       CALL MPE$QUERY(NTASKS,THISTASK)
        LENG1=LENG/NTASKS
        IF(LENG1*NTASKS.NE.LENG) THEN
          CALL ERROR$MSG('ARRAY SIZE NOT DIVIDABLE BY #(PROCESSORS)')
          CALL ERROR$STOP('MPE$TRANSPOSE<TYPEID><RANKID>')
        END IF
        ALLOCATE(VAL1(LENG))
                                   ! VAL1=RESHAPE(VAL,(/LENG/))
        val1(:)=transfer(val,val1) ! reshape failed for the pathscale compiler
	CALL MPE$SYNC
        CALL MPI_ALLTOALL(VAL1,LENG1,<MPI_TYPE>,VAL,LENG1,<MPI_TYPE> &
     &                  ,MPI_COMM_WORLD,IERR)
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
                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
                 ([R8],[REAL(8)][MPI_REAL8])
                 ([R4],[REAL(4)][MPI_REAL4])
                 ([I4],[INTEGER(4)][MPI_INTEGER4])
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
      SUBROUTINE MPE$GATHER<TYPEID><RANKID>(TOTASK,INVAL,OUTVAL)
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
      INTEGER(4),INTENT(IN) :: TOTASK
      <TYPE>    ,INTENT(IN) :: INVAL<RANK>
      <TYPE>    ,INTENT(OUT):: OUTVAL<RANK+1>
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
        LENG=<SIZE>
        CALL MPI_GATHER(INVAL,LENG,<MPI_TYPE>,OUTVAL,LENG,<MPI_TYPE> &
     &                 ,TOTASK-1,MPI_COMM_WORLD,IERR)
#ELSE
        OUTVAL<RANK,1>=INVAL<RANK>
#ENDIF
      RETURN
      END SUBROUTINE MPE$GATHER<TYPEID><RANKID>
#END TEMPLATE MPE$GATHER
END MODULE MPE_GATHER_MODULE
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
      INTEGER(4)               :: IERR
      INTEGER(4)  ,PARAMETER   :: MBYTE=2**20
      INTEGER(4)  ,PARAMETER   :: BUFFER_SIZE=6*MBYTE
      CHARACTER(1),POINTER     :: BUFFER(:)
!     ******************************************************************
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
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NTASKS,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,THISTASK,IERR)
      THISTASK= THISTASK+1
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
      INTEGER(4)        :: IERR
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
      INTEGER(4)            :: IERR,IERROR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      WRITE(*,*) 'EXIT:',NCODE_
      CALL MPI_ABORT(MPI_COMM_WORLD,IERROR,7777)
      STOP 'ERROR STOP'
#ELSE
      WRITE(*,*) 'EXIT:',NCODE_
      STOP 'ERROR STOP'
#ENDIF
      END
!
!     ..................................................................
      SUBROUTINE MPE$QUERY(NTASKS_,THISTASK_)
!     ******************************************************************
!     ** MPE$QUERY                                                    **
!     == RETURN THE TOTAL NUMBER OF TASKS                             ==
!     == AND THE ID OF THIS TASK [1 ..NTASKS]                         ==
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE 
      INTEGER(4),INTENT(OUT) :: NTASKS_
      INTEGER(4),INTENT(OUT) :: THISTASK_
      INTEGER(4)             :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
!     CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NTASKS,IERR)
!     CALL MPI_COMM_RANK(MPI_COMM_WORLD,THISTASK,IERR)
!     THISTASK= THISTASK+1
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
      SUBROUTINE MPE$SYNC
!     ******************************************************************
!     **                                                              **
!     **  MPE$SYNC                                                    **
!     **                                                              **
!     **  PURPOSE: SYNCHRONIZE ALL TASKS                              **
!     **  (WAITS UNTILL ALL TASKS ARRIVED AT THIS POINT OF THE CODE)  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4) :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
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
      SUBROUTINE MPE$INDEX(INVAL,OUTVAL,SIZE)
!     ******************************************************************
!     **                                                              **
!     ** MPE$INDEX                                                    **
!     **                                                              **
!     ** 2 TASKS AOUT(TASK1) := (1,8), AOUT(TASK2) := (4,9)           **
!     **         ARES(TASK1)  = (1,4), AOUT(TASK2)  = (8,9)           **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: SIZE
      CHARACTER(1),INTENT(IN)  :: INVAL(SIZE)
      CHARACTER(1),INTENT(OUT) :: OUTVAL(SIZE)
      INTEGER(4)               :: IERR
!     ******************************************************************
#IFDEF CPPVARIABLE_PARALLEL
        CALL MPI_ALLTOALL(INVAL,SIZE,MPI_BYTE,OUTVAL,SIZE &
     &                  ,MPI_BYTE,MPI_COMM_WORLD,IERR)
#ELSE
        OUTVAL(:)=INVAL(:)
#ENDIF
      RETURN
      END
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
USE MPE_SEND_MODULE
USE MPE_RECEIVE_MODULE
USE MPE_COMBINE_MODULE
USE MPE_BROADCAST_MODULE
END MODULE MPE_MODULE
