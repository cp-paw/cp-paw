!***********************************************************************
!**                                                                   **
!**  NAME: MPELIB                                                     **
!**                                                                   **
!**  PURPOSE: PERFORM COMMUNICATIONS AMONG PROCESSORS FOR             **
!**   PARALLEL PROCESSING.                                            **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    MPE__INIT                                                       **
!**    MPE__QUERY(NTASKS,THISTASK)                                     **
!**    MPE__SYNC                                                       **
!**    MPE__STOPALL(RETURNCODE)                                        **
!**    MPE__EXIT                                                       **
!**    MPE__ISPARALLEL(IPARALLEL) ??                                   **
!**    MPE__COMBINE(OPERATION,VAL)                                     **
!**    MPE__BROADCAST(FROMTASK,VAL)                                    **
!**    MPE__SEND(TOTASK,TAG,VAL)                                       **
!**    MPE__RECEIVE(FROMTASK,TAG,VAL)                                  **
!**    MPE__TRANSPOSE(VAL)                                             **
!**    MPE__GATHER(TOTASK,INVAL,OUTVAL)                                **
!**                                                                   **
!**  remarks:                                                         **
!**    1) INTERFACE to MPI (message passing library)                  **
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
END MODULE MPE_MPIF_MODULE
!
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_COMBINE_MODULE
PUBLIC MPE__COMBINE
INTERFACE MPE__COMBINE    !(OP,VAL)
  MODULE PROCEDURE MPE__COMBINEC8R0
  MODULE PROCEDURE MPE__COMBINER8R0
  MODULE PROCEDURE MPE__COMBINEI4R0
  MODULE PROCEDURE MPE__COMBINEC8R1
  MODULE PROCEDURE MPE__COMBINER8R1
  MODULE PROCEDURE MPE__COMBINEI4R1
  MODULE PROCEDURE MPE__COMBINEC8R2
  MODULE PROCEDURE MPE__COMBINER8R2
  MODULE PROCEDURE MPE__COMBINEI4R2
  MODULE PROCEDURE MPE__COMBINEC8R3
  MODULE PROCEDURE MPE__COMBINER8R3
  MODULE PROCEDURE MPE__COMBINEI4R3
  MODULE PROCEDURE MPE__COMBINEC8R4
  MODULE PROCEDURE MPE__COMBINER8R4
  MODULE PROCEDURE MPE__COMBINEI4R4
  MODULE PROCEDURE MPE__COMBINEC8R5
  MODULE PROCEDURE MPE__COMBINER8R5
  MODULE PROCEDURE MPE__COMBINEI4R5
END INTERFACE
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!#TEMPLATE MPE__COMBINE
!(<TYPEID><TYPE><MPI_TYPE>)
!                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
!                 ([R8],[REAL(8)][MPI_REAL8])
!                 ([I4],[INTEGER(4)][MPI_INTEGER4])
!(<RANKID><SIZE><RESHAPE(VAL)><RANK>)
!                 =([R0],[1],[VAL],[])
!                 =([R1],[SIZE(VAL)],[RESHAPE(VAL,(/LENG/))],[(:)])
!                 =([R2],[SIZE(VAL)],[RESHAPE(VAL,(/LENG/))],[(:,:)])
!                 =([R3],[SIZE(VAL)],[RESHAPE(VAL,(/LENG/))],[(:,:,:)])
!                 =([R4],[SIZE(VAL)],[RESHAPE(VAL,(/LENG/))],[(:,:,:,:)])
!                 =([R5],[SIZE(VAL)],[RESHAPE(VAL,(/LENG/))],[(:,:,:,:,:)])
!#BODY
!!
!!     ..................................................................
!      SUBROUTINE MPE__COMBINE<TYPEID><RANKID>(OPERATION,VAL)
!!     ******************************************************************
!!     **                                                              **
!!     **  NAME: MPE__COMBINER8                                         **
!!     **                                                              **
!!     **  PERFORMS A VECTOR OPERATION ACCROSS ALL TASKS               **
!!     **    KEY(THISTASK)=KEY(TASK1).OP. ... .OP.KEY(NTASKS)          **
!!     **  WHERE KEY(ITASK) IS A VECTOR, AND .OP. IS AN OPERATION      **
!!     **  SUCH AS OP='+' OR OP='*'                                    **
!!     **                                                              **
!!     **  REMARKS:                                                    **
!!     **    IS IT ALLOWED TO REPLACE MPI_DOUBLE BY MPI_REAL8,ETC??    **
!!     **                                                              **
!!     ******************************************************************
!      USE MPE_MPIF_MODULE
!      IMPLICIT NONE
!      CHARACTER(*),INTENT(IN)    :: OPERATION
!      <TYPE>      ,INTENT(INOUT) :: VAL<RANK>
!      <TYPE>      ,ALLOCATABLE   :: WORK(:)
!      INTEGER(4)                 :: LENG
!      INTEGER(4)                 :: IERR
!!     ******************************************************************
!#ifdef PARALLEL_MACHINE
!!
!!     =================================================================
!!     == ALLOCATE WORK ARRAY                                         ==
!!     =================================================================
!      LENG=<SIZE>
!      ALLOCATE(WORK(LENG))
!      WORK(:)=<RESHAPE(VAL)>
!!
!!     =================================================================
!!     == PERFORM OPERATION                                           ==
!!     =================================================================
!      IF(OPERATION.EQ.'+')THEN
!        CALL MPI_ALLREDUCE(WORK,VAL,LENG &
!     &           ,<MPI_TYPE>,MPI_SUM,MPI_COMM_WORLD,IERR)
!      ELSE IF(OPERATION.EQ.'*')THEN
!        CALL MPI_ALLREDUCE(WORK,VAL,LENG &
!     &           ,<MPI_TYPE>,MPI_PROD,MPI_COMM_WORLD,IERR)
!      ELSE IF(OPERATION.EQ.'MAX')THEN
!        CALL MPI_ALLREDUCE(WORK,VAL,LENG &
!     &           ,<MPI_TYPE>,MPI_MAX,MPI_COMM_WORLD,IERR)
!      ELSE IF(OPERATION.EQ.'MIN')THEN
!        CALL MPI_ALLREDUCE(WORK,VAL,LENG &
!     &           ,<MPI_TYPE>,MPI_MIN,MPI_COMM_WORLD,IERR)
!      ELSE
!        CALL ERROR__MSG('OPERATION NOT SUPPORTED')
!        CALL ERROR__STOP('MPE__COMBINE')
!      END IF
!      DEALLOCATE(WORK)
!#endif
!      RETURN
!      END SUBROUTINE MPE__COMBINE<TYPEID><RANKID>
!#INSTANCES
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEC8R0(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      COMPLEX(8)      ,INTENT(INOUT) :: VAL
      COMPLEX(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEC8R0
!
!     ..................................................................
      SUBROUTINE MPE__COMBINER8R0(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      REAL(8)      ,INTENT(INOUT) :: VAL
      REAL(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINER8R0
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEI4R0(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      INTEGER(4)      ,INTENT(INOUT) :: VAL
      INTEGER(4)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEI4R0
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEC8R1(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:)
      COMPLEX(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEC8R1
!
!     ..................................................................
      SUBROUTINE MPE__COMBINER8R1(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      REAL(8)      ,INTENT(INOUT) :: VAL(:)
      REAL(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINER8R1
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEI4R1(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEI4R1
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEC8R2(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:,:)
      COMPLEX(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEC8R2
!
!     ..................................................................
      SUBROUTINE MPE__COMBINER8R2(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      REAL(8)      ,INTENT(INOUT) :: VAL(:,:)
      REAL(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINER8R2
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEI4R2(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEI4R2
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEC8R3(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:,:,:)
      COMPLEX(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEC8R3
!
!     ..................................................................
      SUBROUTINE MPE__COMBINER8R3(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      REAL(8)      ,INTENT(INOUT) :: VAL(:,:,:)
      REAL(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINER8R3
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEI4R3(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEI4R3
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEC8R4(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:,:,:,:)
      COMPLEX(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEC8R4
!
!     ..................................................................
      SUBROUTINE MPE__COMBINER8R4(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      REAL(8)      ,INTENT(INOUT) :: VAL(:,:,:,:)
      REAL(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINER8R4
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEI4R4(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEI4R4
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEC8R5(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      COMPLEX(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEC8R5
!
!     ..................................................................
      SUBROUTINE MPE__COMBINER8R5(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      REAL(8)      ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      REAL(8)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINER8R5
!
!     ..................................................................
      SUBROUTINE MPE__COMBINEI4R5(OPERATION,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__COMBINER8                                         **
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
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)      ,ALLOCATABLE   :: WORK(:)
      INTEGER(4)                 :: LENG
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__COMBINEI4R5
!#END TEMPLATE MPE__COMBINE
END MODULE MPE_COMBINE_MODULE
!
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_BROADCAST_MODULE
PUBLIC MPE__BROADCAST
INTERFACE MPE__BROADCAST  !(FROMTASK,VAL)
  MODULE PROCEDURE MPE__BROADCASTC8R0
  MODULE PROCEDURE MPE__BROADCASTR8R0
  MODULE PROCEDURE MPE__BROADCASTR4R0
  MODULE PROCEDURE MPE__BROADCASTI4R0
  MODULE PROCEDURE MPE__BROADCASTL4R0
  MODULE PROCEDURE MPE__BROADCASTC8R1
  MODULE PROCEDURE MPE__BROADCASTR8R1
  MODULE PROCEDURE MPE__BROADCASTR4R1
  MODULE PROCEDURE MPE__BROADCASTI4R1
  MODULE PROCEDURE MPE__BROADCASTL4R1
  MODULE PROCEDURE MPE__BROADCASTC8R2
  MODULE PROCEDURE MPE__BROADCASTR8R2
  MODULE PROCEDURE MPE__BROADCASTR4R2
  MODULE PROCEDURE MPE__BROADCASTI4R2
  MODULE PROCEDURE MPE__BROADCASTL4R2
  MODULE PROCEDURE MPE__BROADCASTC8R3
  MODULE PROCEDURE MPE__BROADCASTR8R3
  MODULE PROCEDURE MPE__BROADCASTR4R3
  MODULE PROCEDURE MPE__BROADCASTI4R3
  MODULE PROCEDURE MPE__BROADCASTL4R3
  MODULE PROCEDURE MPE__BROADCASTC8R4
  MODULE PROCEDURE MPE__BROADCASTR8R4
  MODULE PROCEDURE MPE__BROADCASTR4R4
  MODULE PROCEDURE MPE__BROADCASTI4R4
  MODULE PROCEDURE MPE__BROADCASTL4R4
  MODULE PROCEDURE MPE__BROADCASTC8R5
  MODULE PROCEDURE MPE__BROADCASTR8R5
  MODULE PROCEDURE MPE__BROADCASTR4R5
  MODULE PROCEDURE MPE__BROADCASTI4R5
  MODULE PROCEDURE MPE__BROADCASTL4R5
  MODULE PROCEDURE MPE__BROADCASTC8R6
  MODULE PROCEDURE MPE__BROADCASTR8R6
  MODULE PROCEDURE MPE__BROADCASTR4R6
  MODULE PROCEDURE MPE__BROADCASTI4R6
  MODULE PROCEDURE MPE__BROADCASTL4R6
    MODULE PROCEDURE MPE__BROADCASTCHA
    MODULE PROCEDURE MPE__BROADCASTCH
END INTERFACE
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!#TEMPLATE MPE__BROADCAST
!(<TYPEID><TYPE><MPI_TYPE>)
!                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
!                 ([R8],[REAL(8)][MPI_REAL8])
!                 ([R4],[REAL(4)][MPI_REAL4])
!                 ([I4],[INTEGER(4)][MPI_INTEGER4])
!                 ([L4],[LOGICAL(4)][MPI_LOGICAL])
!(<RANKID><SIZE><RANK>)
!                 =([R0],[1],[])
!                 =([R1],[SIZE(VAL)],[(:)])
!                 =([R2],[SIZE(VAL)],[(:,:)])
!                 =([R3],[SIZE(VAL)],[(:,:,:)])
!                 =([R4],[SIZE(VAL)],[(:,:,:,:)])
!                 =([R5],[SIZE(VAL)],[(:,:,:,:,:)])
!                 =([R6],[SIZE(VAL)],[(:,:,:,:,:,:)])
!#BODY
!!
!!     ..................................................................
!      SUBROUTINE MPE__BROADCAST<TYPEID><RANKID>(FROMTASK,VAL)
!!     ******************************************************************
!!     **                                                              **
!!     ** MPE__BROADCASTR8                                              **
!!     **                                                              **
!!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!!     **                                                              **
!!     ******************************************************************
!      use mpe_mpif_module
!      IMPLICIT NONE
!      <type>      ,INTENT(INOUT) :: VAL<rank>
!      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
!      INTEGER(4)                 :: IERR
!!     ******************************************************************
!#ifdef PARALLEL_MACHINE
!      CALL MPI_BCAST(VAL,<SIZE>,<MPI_TYPE>,FROMTASK-1 &
!     &                ,MPI_COMM_WORLD,IERR)
!#endif
!      RETURN
!      END SUBROUTINE mpe__BROADCAST<TYPEID><RANKID>
!#INSTANCES
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTC8R0(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      COMPLEX(8)      ,INTENT(INOUT) :: VAL
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTC8R0
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR8R0(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(8)      ,INTENT(INOUT) :: VAL
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR8R0
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR4R0(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(4)      ,INTENT(INOUT) :: VAL
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR4R0
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTI4R0(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(INOUT) :: VAL
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTI4R0
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTL4R0(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      LOGICAL(4)      ,INTENT(INOUT) :: VAL
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTL4R0
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTC8R1(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTC8R1
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR8R1(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(8)      ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR8R1
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR4R1(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(4)      ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR4R1
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTI4R1(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTI4R1
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTL4R1(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      LOGICAL(4)      ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTL4R1
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTC8R2(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTC8R2
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR8R2(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(8)      ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR8R2
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR4R2(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(4)      ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR4R2
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTI4R2(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTI4R2
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTL4R2(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      LOGICAL(4)      ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTL4R2
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTC8R3(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTC8R3
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR8R3(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(8)      ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR8R3
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR4R3(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(4)      ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR4R3
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTI4R3(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTI4R3
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTL4R3(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      LOGICAL(4)      ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTL4R3
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTC8R4(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTC8R4
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR8R4(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(8)      ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR8R4
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR4R4(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(4)      ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR4R4
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTI4R4(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTI4R4
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTL4R4(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      LOGICAL(4)      ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTL4R4
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTC8R5(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTC8R5
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR8R5(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(8)      ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR8R5
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR4R5(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(4)      ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR4R5
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTI4R5(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTI4R5
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTL4R5(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      LOGICAL(4)      ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTL4R5
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTC8R6(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      COMPLEX(8)      ,INTENT(INOUT) :: VAL(:,:,:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTC8R6
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR8R6(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(8)      ,INTENT(INOUT) :: VAL(:,:,:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR8R6
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTR4R6(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      REAL(4)      ,INTENT(INOUT) :: VAL(:,:,:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTR4R6
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTI4R6(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(INOUT) :: VAL(:,:,:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTI4R6
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTL4R6(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      LOGICAL(4)      ,INTENT(INOUT) :: VAL(:,:,:,:,:,:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE mpe__BROADCASTL4R6
!#END TEMPLATE MPE__BROADCAST
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTCHA(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      CHARACTER(*),INTENT(INOUT) :: VAL(:)
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: leng
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      end subroutine mpe__broadcastcha
!
!     ..................................................................
      SUBROUTINE MPE__BROADCASTCH(FROMTASK,VAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__BROADCASTR8                                              **
!     **                                                              **
!     ** BROADCAST A MESSAGE FROM TASK NST TO THE REST OF THE WORLD   **
!     **   VAL(THISTASK)=VAL(FROMTASK)                                **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      CHARACTER(*),INTENT(INOUT) :: VAL
      INTEGER(4)  ,INTENT(IN)    :: FROMTASK
      INTEGER(4)                 :: leng
      INTEGER(4)                 :: IERR
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__BROADCASTCH
END MODULE MPE_broadcast_MODULE
!
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
module mpe_send_module
PUBLIC MPE__SEND
INTERFACE MPE__SEND       !(TOTASK,TAG,VAL)
  MODULE PROCEDURE MPE__SENDC8R0
  MODULE PROCEDURE MPE__SENDR8R0
  MODULE PROCEDURE MPE__SENDR4R0
  MODULE PROCEDURE MPE__SENDI4R0
  MODULE PROCEDURE MPE__SENDL4R0
  MODULE PROCEDURE MPE__SENDC8R1
  MODULE PROCEDURE MPE__SENDR8R1
  MODULE PROCEDURE MPE__SENDR4R1
  MODULE PROCEDURE MPE__SENDI4R1
  MODULE PROCEDURE MPE__SENDL4R1
  MODULE PROCEDURE MPE__SENDC8R2
  MODULE PROCEDURE MPE__SENDR8R2
  MODULE PROCEDURE MPE__SENDR4R2
  MODULE PROCEDURE MPE__SENDI4R2
  MODULE PROCEDURE MPE__SENDL4R2
  MODULE PROCEDURE MPE__SENDC8R3
  MODULE PROCEDURE MPE__SENDR8R3
  MODULE PROCEDURE MPE__SENDR4R3
  MODULE PROCEDURE MPE__SENDI4R3
  MODULE PROCEDURE MPE__SENDL4R3
  MODULE PROCEDURE MPE__SENDC8R4
  MODULE PROCEDURE MPE__SENDR8R4
  MODULE PROCEDURE MPE__SENDR4R4
  MODULE PROCEDURE MPE__SENDI4R4
  MODULE PROCEDURE MPE__SENDL4R4
  MODULE PROCEDURE MPE__SENDC8R5
  MODULE PROCEDURE MPE__SENDR8R5
  MODULE PROCEDURE MPE__SENDR4R5
  MODULE PROCEDURE MPE__SENDI4R5
  MODULE PROCEDURE MPE__SENDL4R5
  MODULE PROCEDURE MPE__SENDC8R6
  MODULE PROCEDURE MPE__SENDR8R6
  MODULE PROCEDURE MPE__SENDR4R6
  MODULE PROCEDURE MPE__SENDI4R6
  MODULE PROCEDURE MPE__SENDL4R6
    MODULE PROCEDURE MPE__SENDCHA
    MODULE PROCEDURE MPE__SENDCH
END INTERFACE
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!#TEMPLATE MPE__SEND
!(<TYPEID><TYPE><MPI_TYPE>)
!                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
!                 ([R8],[REAL(8)][MPI_REAL8])
!                 ([R4],[REAL(4)][MPI_REAL4])
!                 ([I4],[INTEGER(4)][MPI_INTEGER4])
!                 ([L4],[LOGICAL(4)][MPI_LOGICAL])
!(<RANKID><SIZE><RANK>)
!                 =([R0],[1],[])
!                 =([R1],[SIZE(VAL)],[(:)])
!                 =([R2],[SIZE(VAL)],[(:,:)])
!                 =([R3],[SIZE(VAL)],[(:,:,:)])
!                 =([R4],[SIZE(VAL)],[(:,:,:,:)])
!                 =([R5],[SIZE(VAL)],[(:,:,:,:,:)])
!                 =([R6],[SIZE(VAL)],[(:,:,:,:,:,:)])
!#BODY
!!     ..................................................................
!      SUBROUTINE MPE__SEND<TYPEID><RANKID>(TOTASK,TAG,VAL)
!!     ******************************************************************
!!     **                                                              **
!!     **  NAME: MPE__SEND                                              **
!!     **                                                              **
!!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
!!     **     VAL(TOTASK)=VAL(THISTASK)                                **
!!     **                                                              **
!!     **  REMARKS:                                                    **
!!     **    TAG IS SOME INTEGER NUMBER WHICH ALLOWES TO DIFFERENTIATE **
!!     **    BETWEEN DIFFERENT MESSAGES.                               **
!!     **                                                              **
!!     ******************************************************************
!      USE MPE_MPIF_MODULE
!      IMPLICIT NONE
!      INTEGER(4)  ,INTENT(IN)    :: TOTASK
!      INTEGER(4)  ,INTENT(IN)    :: TAG
!      <TYPE>      ,INTENT(IN)    :: VAL<RANK>
!      INTEGER(4)                 :: IERR,IDUM
!!     ******************************************************************
!#ifdef  PARALLEL_MACHINE
!        CALL MPI_SEND(VAL,<SIZE>,<MPI_TYPE>,TOTASK-1,TAG &
!     &               ,MPI_COMM_WORLD,IERR)
!#else
!        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
!        CALL ERROR__STOP('MPE__SEND<typeid><rankid>')
!#endif
!      RETURN
!      END SUBROUTINE mpe__SEND<typeid><rankid>
!#INSTANCES
!     ..................................................................
      SUBROUTINE MPE__SENDC8R0(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      COMPLEX(8)      ,INTENT(IN)    :: VAL
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDC8R0')
      RETURN
      END SUBROUTINE mpe__SENDC8R0
!     ..................................................................
      SUBROUTINE MPE__SENDR8R0(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(8)      ,INTENT(IN)    :: VAL
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR8R0')
      RETURN
      END SUBROUTINE mpe__SENDR8R0
!     ..................................................................
      SUBROUTINE MPE__SENDR4R0(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(4)      ,INTENT(IN)    :: VAL
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR4R0')
      RETURN
      END SUBROUTINE mpe__SENDR4R0
!     ..................................................................
      SUBROUTINE MPE__SENDI4R0(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      INTEGER(4)      ,INTENT(IN)    :: VAL
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDI4R0')
      RETURN
      END SUBROUTINE mpe__SENDI4R0
!     ..................................................................
      SUBROUTINE MPE__SENDL4R0(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      LOGICAL(4)      ,INTENT(IN)    :: VAL
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDL4R0')
      RETURN
      END SUBROUTINE mpe__SENDL4R0
!     ..................................................................
      SUBROUTINE MPE__SENDC8R1(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      COMPLEX(8)      ,INTENT(IN)    :: VAL(:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDC8R1')
      RETURN
      END SUBROUTINE mpe__SENDC8R1
!     ..................................................................
      SUBROUTINE MPE__SENDR8R1(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(8)      ,INTENT(IN)    :: VAL(:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR8R1')
      RETURN
      END SUBROUTINE mpe__SENDR8R1
!     ..................................................................
      SUBROUTINE MPE__SENDR4R1(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(4)      ,INTENT(IN)    :: VAL(:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR4R1')
      RETURN
      END SUBROUTINE mpe__SENDR4R1
!     ..................................................................
      SUBROUTINE MPE__SENDI4R1(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      INTEGER(4)      ,INTENT(IN)    :: VAL(:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDI4R1')
      RETURN
      END SUBROUTINE mpe__SENDI4R1
!     ..................................................................
      SUBROUTINE MPE__SENDL4R1(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      LOGICAL(4)      ,INTENT(IN)    :: VAL(:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDL4R1')
      RETURN
      END SUBROUTINE mpe__SENDL4R1
!     ..................................................................
      SUBROUTINE MPE__SENDC8R2(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      COMPLEX(8)      ,INTENT(IN)    :: VAL(:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDC8R2')
      RETURN
      END SUBROUTINE mpe__SENDC8R2
!     ..................................................................
      SUBROUTINE MPE__SENDR8R2(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(8)      ,INTENT(IN)    :: VAL(:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR8R2')
      RETURN
      END SUBROUTINE mpe__SENDR8R2
!     ..................................................................
      SUBROUTINE MPE__SENDR4R2(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(4)      ,INTENT(IN)    :: VAL(:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR4R2')
      RETURN
      END SUBROUTINE mpe__SENDR4R2
!     ..................................................................
      SUBROUTINE MPE__SENDI4R2(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      INTEGER(4)      ,INTENT(IN)    :: VAL(:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDI4R2')
      RETURN
      END SUBROUTINE mpe__SENDI4R2
!     ..................................................................
      SUBROUTINE MPE__SENDL4R2(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      LOGICAL(4)      ,INTENT(IN)    :: VAL(:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDL4R2')
      RETURN
      END SUBROUTINE mpe__SENDL4R2
!     ..................................................................
      SUBROUTINE MPE__SENDC8R3(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      COMPLEX(8)      ,INTENT(IN)    :: VAL(:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDC8R3')
      RETURN
      END SUBROUTINE mpe__SENDC8R3
!     ..................................................................
      SUBROUTINE MPE__SENDR8R3(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(8)      ,INTENT(IN)    :: VAL(:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR8R3')
      RETURN
      END SUBROUTINE mpe__SENDR8R3
!     ..................................................................
      SUBROUTINE MPE__SENDR4R3(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(4)      ,INTENT(IN)    :: VAL(:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR4R3')
      RETURN
      END SUBROUTINE mpe__SENDR4R3
!     ..................................................................
      SUBROUTINE MPE__SENDI4R3(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      INTEGER(4)      ,INTENT(IN)    :: VAL(:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDI4R3')
      RETURN
      END SUBROUTINE mpe__SENDI4R3
!     ..................................................................
      SUBROUTINE MPE__SENDL4R3(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      LOGICAL(4)      ,INTENT(IN)    :: VAL(:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDL4R3')
      RETURN
      END SUBROUTINE mpe__SENDL4R3
!     ..................................................................
      SUBROUTINE MPE__SENDC8R4(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      COMPLEX(8)      ,INTENT(IN)    :: VAL(:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDC8R4')
      RETURN
      END SUBROUTINE mpe__SENDC8R4
!     ..................................................................
      SUBROUTINE MPE__SENDR8R4(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(8)      ,INTENT(IN)    :: VAL(:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR8R4')
      RETURN
      END SUBROUTINE mpe__SENDR8R4
!     ..................................................................
      SUBROUTINE MPE__SENDR4R4(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(4)      ,INTENT(IN)    :: VAL(:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR4R4')
      RETURN
      END SUBROUTINE mpe__SENDR4R4
!     ..................................................................
      SUBROUTINE MPE__SENDI4R4(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      INTEGER(4)      ,INTENT(IN)    :: VAL(:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDI4R4')
      RETURN
      END SUBROUTINE mpe__SENDI4R4
!     ..................................................................
      SUBROUTINE MPE__SENDL4R4(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      LOGICAL(4)      ,INTENT(IN)    :: VAL(:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDL4R4')
      RETURN
      END SUBROUTINE mpe__SENDL4R4
!     ..................................................................
      SUBROUTINE MPE__SENDC8R5(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      COMPLEX(8)      ,INTENT(IN)    :: VAL(:,:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDC8R5')
      RETURN
      END SUBROUTINE mpe__SENDC8R5
!     ..................................................................
      SUBROUTINE MPE__SENDR8R5(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(8)      ,INTENT(IN)    :: VAL(:,:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR8R5')
      RETURN
      END SUBROUTINE mpe__SENDR8R5
!     ..................................................................
      SUBROUTINE MPE__SENDR4R5(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(4)      ,INTENT(IN)    :: VAL(:,:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR4R5')
      RETURN
      END SUBROUTINE mpe__SENDR4R5
!     ..................................................................
      SUBROUTINE MPE__SENDI4R5(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      INTEGER(4)      ,INTENT(IN)    :: VAL(:,:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDI4R5')
      RETURN
      END SUBROUTINE mpe__SENDI4R5
!     ..................................................................
      SUBROUTINE MPE__SENDL4R5(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      LOGICAL(4)      ,INTENT(IN)    :: VAL(:,:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDL4R5')
      RETURN
      END SUBROUTINE mpe__SENDL4R5
!     ..................................................................
      SUBROUTINE MPE__SENDC8R6(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      COMPLEX(8)      ,INTENT(IN)    :: VAL(:,:,:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDC8R6')
      RETURN
      END SUBROUTINE mpe__SENDC8R6
!     ..................................................................
      SUBROUTINE MPE__SENDR8R6(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(8)      ,INTENT(IN)    :: VAL(:,:,:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR8R6')
      RETURN
      END SUBROUTINE mpe__SENDR8R6
!     ..................................................................
      SUBROUTINE MPE__SENDR4R6(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      REAL(4)      ,INTENT(IN)    :: VAL(:,:,:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDR4R6')
      RETURN
      END SUBROUTINE mpe__SENDR4R6
!     ..................................................................
      SUBROUTINE MPE__SENDI4R6(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      INTEGER(4)      ,INTENT(IN)    :: VAL(:,:,:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDI4R6')
      RETURN
      END SUBROUTINE mpe__SENDI4R6
!     ..................................................................
      SUBROUTINE MPE__SENDL4R6(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
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
      LOGICAL(4)      ,INTENT(IN)    :: VAL(:,:,:,:,:,:)
      INTEGER(4)                 :: IERR,IDUM
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDL4R6')
      RETURN
      END SUBROUTINE mpe__SENDL4R6
!#END TEMPLATE MPE__SEND
!
!     ..................................................................
      SUBROUTINE MPE__SENDCHA(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
!     **     VAL(TOTASK)=VAL(THISTASK)                                **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    TAG IS SOME INTEGER NUMBER WHICH ALLOWES TO DIFFERENTIATE **
!     **    BETWEEN DIFFERENT MESSAGES.                               **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)    :: TOTASK
      INTEGER(4)  ,INTENT(IN)    :: TAG
      character(*),INTENT(IN)    :: VAL(:)
      INTEGER(4)                 :: IERR
!     ******************************************************************
        CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__SENDch')
      RETURN
      END SUBROUTINE mpe__SENDchA
!
!     ..................................................................
      SUBROUTINE MPE__SENDCH(TOTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__SEND                                              **
!     **                                                              **
!     **  PURPOSE: SEND BLOCKING. SEND VAL TO TOTASK                  **
!     **    AND WAIT UNTIL IT IS RECEIVED BY MPE__RECEIVE              **
!     **     VAL(TOTASK)=VAL(THISTASK)                                **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    TAG IS SOME INTEGER NUMBER WHICH ALLOWES TO DIFFERENTIATE **
!     **    BETWEEN DIFFERENT MESSAGES.                               **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: TOTASK
      INTEGER(4)  ,INTENT(IN) :: TAG
      CHARACTER(*),INTENT(IN) :: VAL
      INTEGER(4)              :: IERR
!     ******************************************************************
      CALL ERROR__MSG('MPE__SEND MUST NOT BE CALLED IN SCALAR MODE')
      CALL ERROR__STOP('MPE__SENDch')
      RETURN
      END SUBROUTINE mpe__SENDch
END MODULE MPE_send_moDULE
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_RECEIVE_MODULE
PUBLIC MPE__RECEIVE
INTERFACE MPE__RECEIVE       !(FOMTASK,TAG,VAL)
  MODULE PROCEDURE MPE__RECEIVEC8R0
  MODULE PROCEDURE MPE__RECEIVER8R0
  MODULE PROCEDURE MPE__RECEIVER4R0
  MODULE PROCEDURE MPE__RECEIVEI4R0
  MODULE PROCEDURE MPE__RECEIVEL4R0
  MODULE PROCEDURE MPE__RECEIVEC8R1
  MODULE PROCEDURE MPE__RECEIVER8R1
  MODULE PROCEDURE MPE__RECEIVER4R1
  MODULE PROCEDURE MPE__RECEIVEI4R1
  MODULE PROCEDURE MPE__RECEIVEL4R1
  MODULE PROCEDURE MPE__RECEIVEC8R2
  MODULE PROCEDURE MPE__RECEIVER8R2
  MODULE PROCEDURE MPE__RECEIVER4R2
  MODULE PROCEDURE MPE__RECEIVEI4R2
  MODULE PROCEDURE MPE__RECEIVEL4R2
  MODULE PROCEDURE MPE__RECEIVEC8R3
  MODULE PROCEDURE MPE__RECEIVER8R3
  MODULE PROCEDURE MPE__RECEIVER4R3
  MODULE PROCEDURE MPE__RECEIVEI4R3
  MODULE PROCEDURE MPE__RECEIVEL4R3
  MODULE PROCEDURE MPE__RECEIVEC8R4
  MODULE PROCEDURE MPE__RECEIVER8R4
  MODULE PROCEDURE MPE__RECEIVER4R4
  MODULE PROCEDURE MPE__RECEIVEI4R4
  MODULE PROCEDURE MPE__RECEIVEL4R4
  MODULE PROCEDURE MPE__RECEIVEC8R5
  MODULE PROCEDURE MPE__RECEIVER8R5
  MODULE PROCEDURE MPE__RECEIVER4R5
  MODULE PROCEDURE MPE__RECEIVEI4R5
  MODULE PROCEDURE MPE__RECEIVEL4R5
  MODULE PROCEDURE MPE__RECEIVEC8R6
  MODULE PROCEDURE MPE__RECEIVER8R6
  MODULE PROCEDURE MPE__RECEIVER4R6
  MODULE PROCEDURE MPE__RECEIVEI4R6
  MODULE PROCEDURE MPE__RECEIVEL4R6
    MODULE PROCEDURE MPE__RECEIVECHA
    MODULE PROCEDURE MPE__RECEIVECH
END INTERFACE
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!#TEMPLATE MPE__RECEIVE
!(<TYPEID><TYPE><MPI_TYPE>)
!                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
!                 ([R8],[REAL(8)][MPI_REAL8])
!                 ([R4],[REAL(4)][MPI_REAL4])
!                 ([I4],[INTEGER(4)][MPI_INTEGER4])
!                 ([L4],[LOGICAL(4)][MPI_LOGICAL])
!(<RANKID><SIZE><RANK>)
!                 =([R0],[1],[])
!                 =([R1],[SIZE(VAL)],[(:)])
!                 =([R2],[SIZE(VAL)],[(:,:)])
!                 =([R3],[SIZE(VAL)],[(:,:,:)])
!                 =([R4],[SIZE(VAL)],[(:,:,:,:)])
!                 =([R5],[SIZE(VAL)],[(:,:,:,:,:)])
!                 =([R6],[SIZE(VAL)],[(:,:,:,:,:,:)])
!#BODY
!!
!!     ..................................................................
!      SUBROUTINE MPE__RECEIVE<TYPEID><RANKID>(FROMTASK,TAG,VAL)
!!     ******************************************************************
!!     **                                                              **
!!     ** NAME: MPE__RECEIVE                                            **
!!     **                                                              **
!!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!!     **   MPE__BSEND                                                  **
!!     **                                                              **
!!     ******************************************************************
!      USE MPE_MPIF_MODULE
!      IMPLICIT NONE
!      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
!      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
!      <TYPE>    ,INTENT(OUT):: VAL<RANK>
!      INTEGER(4)            :: IERR
!#ifdef  PARALLEL_MACHINE
!      INTEGER(4)            :: STAT(MPI_STATUS_SIZE)
!!     ******************************************************************
!        CALL MPI_RECV(VAL,<SIZE>,<MPI_TYPE>,FROMTASK-1,TAG &
!     &               ,MPI_COMM_WORLD,STAT,IERR)
!#else
!        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
!        CALL ERROR__STOP('MPE__RECEIVER8A')
!#endif
!      RETURN
!      END subroutine MPE__RECEIVE<TYPEID><RANKID>
!#INSTANCES
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEC8R0(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      COMPLEX(8)    ,INTENT(OUT):: VAL
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEC8R0
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER8R0(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(8)    ,INTENT(OUT):: VAL
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER8R0
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER4R0(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(4)    ,INTENT(OUT):: VAL
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER4R0
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEI4R0(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      INTEGER(4)    ,INTENT(OUT):: VAL
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEI4R0
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEL4R0(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      LOGICAL(4)    ,INTENT(OUT):: VAL
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEL4R0
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEC8R1(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      COMPLEX(8)    ,INTENT(OUT):: VAL(:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEC8R1
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER8R1(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(8)    ,INTENT(OUT):: VAL(:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER8R1
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER4R1(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(4)    ,INTENT(OUT):: VAL(:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER4R1
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEI4R1(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      INTEGER(4)    ,INTENT(OUT):: VAL(:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEI4R1
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEL4R1(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      LOGICAL(4)    ,INTENT(OUT):: VAL(:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEL4R1
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEC8R2(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      COMPLEX(8)    ,INTENT(OUT):: VAL(:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEC8R2
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER8R2(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(8)    ,INTENT(OUT):: VAL(:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER8R2
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER4R2(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(4)    ,INTENT(OUT):: VAL(:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER4R2
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEI4R2(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      INTEGER(4)    ,INTENT(OUT):: VAL(:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEI4R2
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEL4R2(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      LOGICAL(4)    ,INTENT(OUT):: VAL(:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEL4R2
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEC8R3(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      COMPLEX(8)    ,INTENT(OUT):: VAL(:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEC8R3
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER8R3(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(8)    ,INTENT(OUT):: VAL(:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER8R3
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER4R3(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(4)    ,INTENT(OUT):: VAL(:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER4R3
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEI4R3(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      INTEGER(4)    ,INTENT(OUT):: VAL(:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEI4R3
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEL4R3(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      LOGICAL(4)    ,INTENT(OUT):: VAL(:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEL4R3
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEC8R4(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      COMPLEX(8)    ,INTENT(OUT):: VAL(:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEC8R4
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER8R4(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(8)    ,INTENT(OUT):: VAL(:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER8R4
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER4R4(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(4)    ,INTENT(OUT):: VAL(:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER4R4
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEI4R4(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      INTEGER(4)    ,INTENT(OUT):: VAL(:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEI4R4
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEL4R4(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      LOGICAL(4)    ,INTENT(OUT):: VAL(:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEL4R4
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEC8R5(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      COMPLEX(8)    ,INTENT(OUT):: VAL(:,:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEC8R5
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER8R5(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(8)    ,INTENT(OUT):: VAL(:,:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER8R5
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER4R5(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(4)    ,INTENT(OUT):: VAL(:,:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER4R5
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEI4R5(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      INTEGER(4)    ,INTENT(OUT):: VAL(:,:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEI4R5
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEL4R5(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      LOGICAL(4)    ,INTENT(OUT):: VAL(:,:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEL4R5
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEC8R6(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      COMPLEX(8)    ,INTENT(OUT):: VAL(:,:,:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEC8R6
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER8R6(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(8)    ,INTENT(OUT):: VAL(:,:,:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER8R6
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVER4R6(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      REAL(4)    ,INTENT(OUT):: VAL(:,:,:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVER4R6
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEI4R6(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      INTEGER(4)    ,INTENT(OUT):: VAL(:,:,:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEI4R6
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVEL4R6(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4),INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      LOGICAL(4)    ,INTENT(OUT):: VAL(:,:,:,:,:,:)
      INTEGER(4)            :: IERR
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVER8A')
      RETURN
      END subroutine MPE__RECEIVEL4R6
!#END TEMPLATE MPE__RECEIVE
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVECHA(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4)  ,INTENT(IN)  :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      CHARACTER(*),INTENT(OUT) :: VAL(:)
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVECHA')
      RETURN
      END subroutine MPE__RECEIVECHA
!
!     ..................................................................
      SUBROUTINE MPE__RECEIVECH(FROMTASK,TAG,VAL)
!     ******************************************************************
!     **                                                              **
!     ** NAME: MPE__RECEIVE                                            **
!     **                                                              **
!     ** PURPOSE: BLOCKED RECEIVE. RECEIVE A MESSAGE SENT BY          **
!     **   MPE__BSEND                                                  **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: FROMTASK !TASK NUMBER OF THE SENDING TASK
      INTEGER(4)  ,INTENT(IN) :: TAG      !IDENTIFIES A PARTICULAR MESSAGE
      CHARACTER(*),INTENT(OUT):: VAL
        CALL ERROR__MSG('MPE__RECEIVE MUST NOT BE CALLED IN SCALAR MODE')
        CALL ERROR__STOP('MPE__RECEIVECHA')
      RETURN
      END SUBROUTINE MPE__RECEIVECH
END MODULE MPE_RECEIVE_MODULE
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_TRANSPOSE_MODULE
PUBLIC MPE__TRANSPOSE
INTERFACE MPE__TRANSPOSE   !(VAL)
  MODULE PROCEDURE MPE__TRANSPOSEC8R1
  MODULE PROCEDURE MPE__TRANSPOSER8R1
  MODULE PROCEDURE MPE__TRANSPOSER4R1
  MODULE PROCEDURE MPE__TRANSPOSEI4R1
  MODULE PROCEDURE MPE__TRANSPOSEL4R1
  MODULE PROCEDURE MPE__TRANSPOSEC8R2
  MODULE PROCEDURE MPE__TRANSPOSER8R2
  MODULE PROCEDURE MPE__TRANSPOSER4R2
  MODULE PROCEDURE MPE__TRANSPOSEI4R2
  MODULE PROCEDURE MPE__TRANSPOSEL4R2
  MODULE PROCEDURE MPE__TRANSPOSEC8R3
  MODULE PROCEDURE MPE__TRANSPOSER8R3
  MODULE PROCEDURE MPE__TRANSPOSER4R3
  MODULE PROCEDURE MPE__TRANSPOSEI4R3
  MODULE PROCEDURE MPE__TRANSPOSEL4R3
  MODULE PROCEDURE MPE__TRANSPOSEC8R4
  MODULE PROCEDURE MPE__TRANSPOSER8R4
  MODULE PROCEDURE MPE__TRANSPOSER4R4
  MODULE PROCEDURE MPE__TRANSPOSEI4R4
  MODULE PROCEDURE MPE__TRANSPOSEL4R4
  MODULE PROCEDURE MPE__TRANSPOSEC8R5
  MODULE PROCEDURE MPE__TRANSPOSER8R5
  MODULE PROCEDURE MPE__TRANSPOSER4R5
  MODULE PROCEDURE MPE__TRANSPOSEI4R5
  MODULE PROCEDURE MPE__TRANSPOSEL4R5
  MODULE PROCEDURE MPE__TRANSPOSEC8R6
  MODULE PROCEDURE MPE__TRANSPOSER8R6
  MODULE PROCEDURE MPE__TRANSPOSER4R6
  MODULE PROCEDURE MPE__TRANSPOSEI4R6
  MODULE PROCEDURE MPE__TRANSPOSEL4R6
END INTERFACE
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!#TEMPLATE MPE__TRANSPOSE
!(<TYPEID><TYPE><MPI_TYPE>)
!                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
!                 ([R8],[REAL(8)][MPI_REAL8])
!                 ([R4],[REAL(4)][MPI_REAL4])
!                 ([I4],[INTEGER(4)][MPI_INTEGER4])
!                 ([L4],[LOGICAL(4)][MPI_LOGICAL])
!(<RANKID><SIZE><RANK>)
!                 =([R1],[SIZE(VAL)],[(:)])
!                 =([R2],[SIZE(VAL)],[(:,:)])
!                 =([R3],[SIZE(VAL)],[(:,:,:)])
!                 =([R4],[SIZE(VAL)],[(:,:,:,:)])
!                 =([R5],[SIZE(VAL)],[(:,:,:,:,:)])
!                 =([R6],[SIZE(VAL)],[(:,:,:,:,:,:)])
!#BODY
!!
!!     ..................................................................
!      SUBROUTINE MPE__TRANSPOSE<TYPEID><RANKID>(VAL)
!!     ******************************************************************
!!     **                                                              **
!!     **  NAME: MPE__TRANSPOSE                                         **
!!     **                                                              **
!!     **  PURPOSE: THE JTH BLOCK OF DATA SENT FROM TASK I IS          **
!!     **    RECEIVED BY THE JTH PROCESS AND PLACED INTO THE ITH BLOCK.**
!!     **                                                              **
!!     **    IN OTHER WORDS:                                           **
!!     **    THINK OF VAL AS VECTOR OF ARRAYS VAL(:,ITASK)             **
!!     **    CONTAINING NTASKS PORTIONS OF SIZE (SIZE(VAL)/NTASKS).    **
!!     **    THEN THE VECTORS FROM ALL TASKS FORM A MATRIX             **
!!     **    (VAL(:,ITASK,JTASK)).                                     **
!!     **    THIS OPERATION TRANSPOSES THIS MATRIX                     **
!!     **        VAL(:,ITASK,THISTASK)= VAL(:,THISTASK,ITASK)          **
!!     **                                                              **
!!     ******************************************************************
!      USE MPE_MPIF_MODULE
!      IMPLICIT NONE
!      <TYPE>    ,INTENT(INOUT) :: VAL<RANK>
!      INTEGER(4)               :: IERR
!      INTEGER(4)               :: LENG,LENG1
!      INTEGER(4)               :: NTASKS,THISTASK
!      <TYPE>    ,ALLOCATABLE   :: VAL1(:)
!!     ******************************************************************
!#ifdef  PARALLEL_MACHINE
!        LENG=<SIZE>
!        CALL MPE__QUERY(NTASKS,THISTASK)
!        LENG1=LENG/NTASKS
!        IF(LENG1*NTASKS.NE.LENG) THEN
!          CALL ERROR__MSG('ARRAY SIZE NOT DIVIDABLE BY #(PROCESSORS)')
!          CALL ERROR__STOP('MPE__TRANSPOSE<TYPEID><RANKID>')
!        END IF
!        ALLOCATE(VAL1(LENG))
!        VAL1=RESHAPE(VAL,(/LENG/))
!        CALL MPI_ALLTOALL(VAL1,LENG1,<MPI_TYPE>,VAL,LENG1,<MPI_TYPE> &
!     &                  ,MPI_COMM_WORLD,IERR)
!        DEALLOCATE(VAL1)
!#endif
!      RETURN
!      END SUBROUTINE MPE__TRANSPOSE<TYPEID><RANKID>
!#INSTANCES
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEC8R1(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      COMPLEX(8)    ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      COMPLEX(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEC8R1
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER8R1(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(8)    ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER8R1
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER4R1(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(4)    ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER4R1
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEI4R1(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      INTEGER(4)    ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEI4R1
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEL4R1(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      LOGICAL(4)    ,INTENT(INOUT) :: VAL(:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      LOGICAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEL4R1
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEC8R2(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      COMPLEX(8)    ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      COMPLEX(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEC8R2
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER8R2(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(8)    ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER8R2
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER4R2(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(4)    ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER4R2
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEI4R2(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      INTEGER(4)    ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEI4R2
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEL4R2(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      LOGICAL(4)    ,INTENT(INOUT) :: VAL(:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      LOGICAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEL4R2
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEC8R3(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      COMPLEX(8)    ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      COMPLEX(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEC8R3
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER8R3(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(8)    ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER8R3
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER4R3(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(4)    ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER4R3
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEI4R3(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      INTEGER(4)    ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEI4R3
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEL4R3(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      LOGICAL(4)    ,INTENT(INOUT) :: VAL(:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      LOGICAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEL4R3
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEC8R4(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      COMPLEX(8)    ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      COMPLEX(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEC8R4
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER8R4(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(8)    ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER8R4
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER4R4(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(4)    ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER4R4
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEI4R4(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      INTEGER(4)    ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEI4R4
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEL4R4(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      LOGICAL(4)    ,INTENT(INOUT) :: VAL(:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      LOGICAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEL4R4
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEC8R5(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      COMPLEX(8)    ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      COMPLEX(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEC8R5
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER8R5(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(8)    ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER8R5
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER4R5(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(4)    ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER4R5
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEI4R5(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      INTEGER(4)    ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEI4R5
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEL4R5(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      LOGICAL(4)    ,INTENT(INOUT) :: VAL(:,:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      LOGICAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEL4R5
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEC8R6(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      COMPLEX(8)    ,INTENT(INOUT) :: VAL(:,:,:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      COMPLEX(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEC8R6
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER8R6(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(8)    ,INTENT(INOUT) :: VAL(:,:,:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(8)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER8R6
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSER4R6(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      REAL(4)    ,INTENT(INOUT) :: VAL(:,:,:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      REAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSER4R6
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEI4R6(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      INTEGER(4)    ,INTENT(INOUT) :: VAL(:,:,:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEI4R6
!
!     ..................................................................
      SUBROUTINE MPE__TRANSPOSEL4R6(VAL)
!     ******************************************************************
!     **                                                              **
!     **  NAME: MPE__TRANSPOSE                                         **
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
      LOGICAL(4)    ,INTENT(INOUT) :: VAL(:,:,:,:,:,:)
      INTEGER(4)               :: IERR
      INTEGER(4)               :: LENG,LENG1
      INTEGER(4)               :: NTASKS,THISTASK
      LOGICAL(4)    ,ALLOCATABLE   :: VAL1(:)
!     ******************************************************************
      RETURN
      END SUBROUTINE MPE__TRANSPOSEL4R6
!#END TEMPLATE MPE__TRANSPOSE
END MODULE MPE_transpose_MODULE
!***********************************************************************
!***********************************************************************
@PROCESS NOEXTCHK
MODULE MPE_GATHER_MODULE
PUBLIC MPE__GATHER
INTERFACE MPE__GATHER    !(TOTASK, INVAL,OUTVAL,VAL)
  MODULE PROCEDURE MPE__GATHERC8R0
  MODULE PROCEDURE MPE__GATHERR8R0
  MODULE PROCEDURE MPE__GATHERR4R0
  MODULE PROCEDURE MPE__GATHERI4R0
  MODULE PROCEDURE MPE__GATHERL4R0
  MODULE PROCEDURE MPE__GATHERC8R1
  MODULE PROCEDURE MPE__GATHERR8R1
  MODULE PROCEDURE MPE__GATHERR4R1
  MODULE PROCEDURE MPE__GATHERI4R1
  MODULE PROCEDURE MPE__GATHERL4R1
  MODULE PROCEDURE MPE__GATHERC8R2
  MODULE PROCEDURE MPE__GATHERR8R2
  MODULE PROCEDURE MPE__GATHERR4R2
  MODULE PROCEDURE MPE__GATHERI4R2
  MODULE PROCEDURE MPE__GATHERL4R2
  MODULE PROCEDURE MPE__GATHERC8R3
  MODULE PROCEDURE MPE__GATHERR8R3
  MODULE PROCEDURE MPE__GATHERR4R3
  MODULE PROCEDURE MPE__GATHERI4R3
  MODULE PROCEDURE MPE__GATHERL4R3
  MODULE PROCEDURE MPE__GATHERC8R4
  MODULE PROCEDURE MPE__GATHERR8R4
  MODULE PROCEDURE MPE__GATHERR4R4
  MODULE PROCEDURE MPE__GATHERI4R4
  MODULE PROCEDURE MPE__GATHERL4R4
  MODULE PROCEDURE MPE__GATHERC8R5
  MODULE PROCEDURE MPE__GATHERR8R5
  MODULE PROCEDURE MPE__GATHERR4R5
  MODULE PROCEDURE MPE__GATHERI4R5
  MODULE PROCEDURE MPE__GATHERL4R5
  MODULE PROCEDURE MPE__GATHERC8R6
  MODULE PROCEDURE MPE__GATHERR8R6
  MODULE PROCEDURE MPE__GATHERR4R6
  MODULE PROCEDURE MPE__GATHERI4R6
  MODULE PROCEDURE MPE__GATHERL4R6
END INTERFACE
CONTAINS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!#TEMPLATE MPE__gather
!(<TYPEID><TYPE><MPI_TYPE>)
!                =([C8],[COMPLEX(8)][MPI_DOUBLE_COMPLEX])
!                 ([R8],[REAL(8)][MPI_REAL8])
!                 ([R4],[REAL(4)][MPI_REAL4])
!                 ([I4],[INTEGER(4)][MPI_INTEGER4])
!                 ([L4],[LOGICAL(4)][MPI_LOGICAL])
!(<RANKID><SIZE><RANK><RANK+1><RANK,1>)
!                 =([R0],[1]          ,[]             ,[(:)]            ,[(1)])
!                 =([R1],[SIZE(INVAL)],[(:)]          ,[(:,:)]          ,[(:,1)])
!                 =([R2],[SIZE(INVAL)],[(:,:)]        ,[(:,:,:)]        ,[(:,:,1)])
!                 =([R3],[SIZE(INVAL)],[(:,:,:)]      ,[(:,:,:,:)]      ,[(:,:,:,1)])
!                 =([R4],[SIZE(INVAL)],[(:,:,:,:)]    ,[(:,:,:,:,:)]    ,[(:,:,:,:,1)])
!                 =([R5],[SIZE(INVAL)],[(:,:,:,:,:)]  ,[(:,:,:,:,:,:)]  ,[(:,:,:,:,:,1)])
!                 =([R6],[SIZE(INVAL)],[(:,:,:,:,:,:)],[(:,:,:,:,:,:,:)],[(:,:,:,:,:,:,1)])
!#BODY
!!
!!     ..................................................................
!      SUBROUTINE MPE__GATHER<TYPEID><RANKID>(TOTASK,INVAL,OUTVAL)
!!     ******************************************************************
!!     **                                                              **
!!     ** MPE__GATHER                                                   **
!!     **                                                              **
!!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!!     ** THE WORLD                                                    **
!!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!!     **                                                              **
!!     ******************************************************************
!      USE MPE_MPIF_MODULE
!      IMPLICIT NONE
!      INTEGER(4),INTENT(IN) :: TOTASK
!      <TYPE>    ,INTENT(IN) :: INVAL<RANK>
!      <TYPE>    ,INTENT(OUT):: OUTVAL<RANK+1>
!      INTEGER(4)            :: LENG
!      INTEGER(4)            :: IERR
!!     ******************************************************************
!#ifdef PARALLEL_MACHINE
!        LENG=<SIZE>
!        CALL MPI_GATHER(INVAL,LENG,<MPI_TYPE>,OUTVAL,LENG,<MPI_TYPE> &
!     &                 ,TOTASK-1,MPI_COMM_WORLD,IERR)
!#else
!        OUTVAL<RANK,1>=INVAL<RANK>
!#endif
!      RETURN
!      END SUBROUTINE MPE__GATHER<TYPEID><RANKID>
!#INSTANCES
!
!     ..................................................................
      SUBROUTINE MPE__GATHERC8R0(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      COMPLEX(8)    ,INTENT(IN) :: INVAL
      COMPLEX(8)    ,INTENT(OUT):: OUTVAL(:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(1)=INVAL
      RETURN
      END SUBROUTINE MPE__GATHERC8R0
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR8R0(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(8)    ,INTENT(IN) :: INVAL
      REAL(8)    ,INTENT(OUT):: OUTVAL(:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(1)=INVAL
      RETURN
      END SUBROUTINE MPE__GATHERR8R0
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR4R0(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(4)    ,INTENT(IN) :: INVAL
      REAL(4)    ,INTENT(OUT):: OUTVAL(:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(1)=INVAL
      RETURN
      END SUBROUTINE MPE__GATHERR4R0
!
!     ..................................................................
      SUBROUTINE MPE__GATHERI4R0(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      INTEGER(4)    ,INTENT(IN) :: INVAL
      INTEGER(4)    ,INTENT(OUT):: OUTVAL(:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(1)=INVAL
      RETURN
      END SUBROUTINE MPE__GATHERI4R0
!
!     ..................................................................
      SUBROUTINE MPE__GATHERL4R0(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      LOGICAL(4)    ,INTENT(IN) :: INVAL
      LOGICAL(4)    ,INTENT(OUT):: OUTVAL(:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(1)=INVAL
      RETURN
      END SUBROUTINE MPE__GATHERL4R0
!
!     ..................................................................
      SUBROUTINE MPE__GATHERC8R1(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      COMPLEX(8)    ,INTENT(IN) :: INVAL(:)
      COMPLEX(8)    ,INTENT(OUT):: OUTVAL(:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,1)=INVAL(:)
      RETURN
      END SUBROUTINE MPE__GATHERC8R1
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR8R1(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(8)    ,INTENT(IN) :: INVAL(:)
      REAL(8)    ,INTENT(OUT):: OUTVAL(:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,1)=INVAL(:)
      RETURN
      END SUBROUTINE MPE__GATHERR8R1
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR4R1(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(4)    ,INTENT(IN) :: INVAL(:)
      REAL(4)    ,INTENT(OUT):: OUTVAL(:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,1)=INVAL(:)
      RETURN
      END SUBROUTINE MPE__GATHERR4R1
!
!     ..................................................................
      SUBROUTINE MPE__GATHERI4R1(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      INTEGER(4)    ,INTENT(IN) :: INVAL(:)
      INTEGER(4)    ,INTENT(OUT):: OUTVAL(:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,1)=INVAL(:)
      RETURN
      END SUBROUTINE MPE__GATHERI4R1
!
!     ..................................................................
      SUBROUTINE MPE__GATHERL4R1(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      LOGICAL(4)    ,INTENT(IN) :: INVAL(:)
      LOGICAL(4)    ,INTENT(OUT):: OUTVAL(:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,1)=INVAL(:)
      RETURN
      END SUBROUTINE MPE__GATHERL4R1
!
!     ..................................................................
      SUBROUTINE MPE__GATHERC8R2(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      COMPLEX(8)    ,INTENT(IN) :: INVAL(:,:)
      COMPLEX(8)    ,INTENT(OUT):: OUTVAL(:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,1)=INVAL(:,:)
      RETURN
      END SUBROUTINE MPE__GATHERC8R2
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR8R2(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(8)    ,INTENT(IN) :: INVAL(:,:)
      REAL(8)    ,INTENT(OUT):: OUTVAL(:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,1)=INVAL(:,:)
      RETURN
      END SUBROUTINE MPE__GATHERR8R2
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR4R2(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(4)    ,INTENT(IN) :: INVAL(:,:)
      REAL(4)    ,INTENT(OUT):: OUTVAL(:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,1)=INVAL(:,:)
      RETURN
      END SUBROUTINE MPE__GATHERR4R2
!
!     ..................................................................
      SUBROUTINE MPE__GATHERI4R2(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      INTEGER(4)    ,INTENT(IN) :: INVAL(:,:)
      INTEGER(4)    ,INTENT(OUT):: OUTVAL(:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,1)=INVAL(:,:)
      RETURN
      END SUBROUTINE MPE__GATHERI4R2
!
!     ..................................................................
      SUBROUTINE MPE__GATHERL4R2(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      LOGICAL(4)    ,INTENT(IN) :: INVAL(:,:)
      LOGICAL(4)    ,INTENT(OUT):: OUTVAL(:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,1)=INVAL(:,:)
      RETURN
      END SUBROUTINE MPE__GATHERL4R2
!
!     ..................................................................
      SUBROUTINE MPE__GATHERC8R3(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      COMPLEX(8)    ,INTENT(IN) :: INVAL(:,:,:)
      COMPLEX(8)    ,INTENT(OUT):: OUTVAL(:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,1)=INVAL(:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERC8R3
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR8R3(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(8)    ,INTENT(IN) :: INVAL(:,:,:)
      REAL(8)    ,INTENT(OUT):: OUTVAL(:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,1)=INVAL(:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERR8R3
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR4R3(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(4)    ,INTENT(IN) :: INVAL(:,:,:)
      REAL(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,1)=INVAL(:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERR4R3
!
!     ..................................................................
      SUBROUTINE MPE__GATHERI4R3(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      INTEGER(4)    ,INTENT(IN) :: INVAL(:,:,:)
      INTEGER(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,1)=INVAL(:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERI4R3
!
!     ..................................................................
      SUBROUTINE MPE__GATHERL4R3(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      LOGICAL(4)    ,INTENT(IN) :: INVAL(:,:,:)
      LOGICAL(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,1)=INVAL(:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERL4R3
!
!     ..................................................................
      SUBROUTINE MPE__GATHERC8R4(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      COMPLEX(8)    ,INTENT(IN) :: INVAL(:,:,:,:)
      COMPLEX(8)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,1)=INVAL(:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERC8R4
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR8R4(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(8)    ,INTENT(IN) :: INVAL(:,:,:,:)
      REAL(8)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,1)=INVAL(:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERR8R4
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR4R4(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(4)    ,INTENT(IN) :: INVAL(:,:,:,:)
      REAL(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,1)=INVAL(:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERR4R4
!
!     ..................................................................
      SUBROUTINE MPE__GATHERI4R4(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      INTEGER(4)    ,INTENT(IN) :: INVAL(:,:,:,:)
      INTEGER(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,1)=INVAL(:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERI4R4
!
!     ..................................................................
      SUBROUTINE MPE__GATHERL4R4(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      LOGICAL(4)    ,INTENT(IN) :: INVAL(:,:,:,:)
      LOGICAL(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,1)=INVAL(:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERL4R4
!
!     ..................................................................
      SUBROUTINE MPE__GATHERC8R5(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      COMPLEX(8)    ,INTENT(IN) :: INVAL(:,:,:,:,:)
      COMPLEX(8)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,:,1)=INVAL(:,:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERC8R5
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR8R5(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(8)    ,INTENT(IN) :: INVAL(:,:,:,:,:)
      REAL(8)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,:,1)=INVAL(:,:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERR8R5
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR4R5(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(4)    ,INTENT(IN) :: INVAL(:,:,:,:,:)
      REAL(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,:,1)=INVAL(:,:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERR4R5
!
!     ..................................................................
      SUBROUTINE MPE__GATHERI4R5(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      INTEGER(4)    ,INTENT(IN) :: INVAL(:,:,:,:,:)
      INTEGER(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,:,1)=INVAL(:,:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERI4R5
!
!     ..................................................................
      SUBROUTINE MPE__GATHERL4R5(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      LOGICAL(4)    ,INTENT(IN) :: INVAL(:,:,:,:,:)
      LOGICAL(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,:,1)=INVAL(:,:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERL4R5
!
!     ..................................................................
      SUBROUTINE MPE__GATHERC8R6(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      COMPLEX(8)    ,INTENT(IN) :: INVAL(:,:,:,:,:,:)
      COMPLEX(8)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,:,:,1)=INVAL(:,:,:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERC8R6
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR8R6(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(8)    ,INTENT(IN) :: INVAL(:,:,:,:,:,:)
      REAL(8)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,:,:,1)=INVAL(:,:,:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERR8R6
!
!     ..................................................................
      SUBROUTINE MPE__GATHERR4R6(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      REAL(4)    ,INTENT(IN) :: INVAL(:,:,:,:,:,:)
      REAL(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,:,:,1)=INVAL(:,:,:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERR4R6
!
!     ..................................................................
      SUBROUTINE MPE__GATHERI4R6(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      INTEGER(4)    ,INTENT(IN) :: INVAL(:,:,:,:,:,:)
      INTEGER(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,:,:,1)=INVAL(:,:,:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERI4R6
!
!     ..................................................................
      SUBROUTINE MPE__GATHERL4R6(TOTASK,INVAL,OUTVAL)
!     ******************************************************************
!     **                                                              **
!     ** MPE__GATHER                                                   **
!     **                                                              **
!     ** THE VALUE OUTVAL ON TOTASK RECEIVES THE INVAL VALUES OF      **
!     ** THE WORLD                                                    **
!     **       OUTVAL(TOTASK)=INVAL(ITASK)                            **
!     **                                                              **
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: TOTASK
      LOGICAL(4)    ,INTENT(IN) :: INVAL(:,:,:,:,:,:)
      LOGICAL(4)    ,INTENT(OUT):: OUTVAL(:,:,:,:,:,:,:)
      INTEGER(4)            :: LENG
      INTEGER(4)            :: IERR
!     ******************************************************************
        OUTVAL(:,:,:,:,:,:,1)=INVAL(:,:,:,:,:,:)
      RETURN
      END SUBROUTINE MPE__GATHERL4R6
!#END TEMPLATE MPE__GATHER
END MODULE MPE_gather_MODULE
!
!     ..................................................................
      SUBROUTINE MPE__INIT
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
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)               :: IERR
      integer(4)  ,PARAMETER   :: mbyte=2**20
      INTEGER(4)  ,PARAMETER   :: BUFFER_SIZE=3*mbyte
      CHARACTER(1),POINTER     :: BUFFER(:)
!     ******************************************************************
!
!     ==================================================================
!     == GET NUMBER OF TASKS , ID OF THIS ONE                         ==
!     == IN THIS ROUTINE , TASKIDS RANGE FROM 0 TO NTASKNUM -1,       ==
!     == WHILE FOR A CALL FROM OUTSIDE, NTASKID S ARE TO BE           ==
!     == OUT OF [1,NTASKNUM]                                          ==
!     ==================================================================
      RETURN
      end
!
!     ..................................................................
      SUBROUTINE MPE__EXIT
!     ******************************************************************
!     ** MPE__EXIT                                                     **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)        :: IERR
!     ******************************************************************
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MPE__STOPALL(NCODE_)
!     ******************************************************************
!     ** MPE__STOPALL                                                  **
!     ** STOP ALL TASKS IN PARTITION WITH EXITCODE NCODE_             **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCODE_
      INTEGER(4)            :: THISTASK
      INTEGER(4)            :: IERR,ierror
!     ******************************************************************
      WRITE(*,*) 'EXIT:',NCODE_
      STOP
      END
!
!     ..................................................................
      SUBROUTINE MPE__QUERY(NTASKS,THISTASK)
!     ******************************************************************
!     ** MPE__QUERY                                                    **
!     == RETURN THE TOTAL NUMBER OF TASKS                             ==
!     == AND THE ID OF THIS TASK [1 ..NTASKS]                         ==
!     ******************************************************************
      USE MPE_MPIF_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: Ntasks
      INTEGER(4),INTENT(OUT) :: thistask
      INTEGER(4)             :: IERR
!     ******************************************************************
      NtASKs   = 1
      thisTASK = 1
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MPE__SYNC
!     ******************************************************************
!     **                                                              **
!     **  MPE__SYNC                                                    **
!     **                                                              **
!     **  PURPOSE: SYNCHRONIZE ALL TASKS                              **
!     **  (WAITS UNTILL ALL TASKS ARRIVED AT THIS POINT OF THE CODE)  **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4) :: IERR
!     ******************************************************************
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MPE__ISPARALLEL(IPARALLEL)
!     ******************************************************************
!     ** MPE__ISPARALLEL                                               **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: IPARALLEL
!     ******************************************************************
      IPARALLEL=0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MPE__INDEX(INVAL,OUTVAL,SIZE)
!     ******************************************************************
!     **                                                              **
!     ** MPE__INDEX                                                    **
!     **                                                              **
!     ** 2 TASKS AOUT(TASK1) := (1,8), AOUT(TASK2) := (4,9)           **
!     **         ARES(TASK1)  = (1,4), AOUT(TASK2)  = (8,9)           **
!     **                                                              **
!     ******************************************************************
      use mpe_mpif_module
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: SIZE
      CHARACTER(1),INTENT(IN)  :: INVAL(SIZE)
      CHARACTER(1),INTENT(OUT) :: OUTVAL(SIZE)
      INTEGER(4)               :: IERR
!     ******************************************************************
        OUTVAL(:)=INVAL(:)
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
USE MPE_COMBINe_MODULE
USE MPE_BROADCAST_MODULE
END MODULE MPE_MODULE
