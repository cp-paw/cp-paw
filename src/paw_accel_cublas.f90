!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****  OPTIONAL NVIDIA cuBLAS/OpenACC ACCELERATION HELPERS                 *****
!****                                                                      *****
!****  This file is intentionally isolated from paw_library.f90 so the      *****
!****  generic BLAS/LAPACK wrappers stay mostly toolchain neutral.          *****
!****                                                                      *****
!*******************************************************************************
!*******************************************************************************
!
#IF DEFINED(CPPVAR_CUBLAS_ACC)
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE CPPAW_CUBLAS_ACC_MODULE
!     **************************************************************************
!     **  cuBLAS/OpenACC bridge for selected complex BLAS3 kernels.           **
!     **  The module is only compiled for NVIDIA HPC SDK accelerator targets. **
!     **************************************************************************
      USE CUBLAS_V2
      USE CUDAFOR
      USE OPENACC
      IMPLICIT NONE
      TYPE(CUBLASHANDLE) :: HANDLE
      LOGICAL(4)         :: HANDLE_READY=.FALSE.
      LOGICAL(4)         :: CONFIG_READY=.FALSE.
      LOGICAL(4)         :: ENABLED=.TRUE.
      LOGICAL(4)         :: SYNC_ENABLED=.TRUE.
      REAL(8)            :: MINFLOP=1.D7
      CONTAINS
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_INITCONFIG
      IMPLICIT NONE
      CHARACTER(128) :: VALUE
      INTEGER(4)     :: STATUS
      INTEGER(4)     :: IOS
!     **************************************************************************
      IF(CONFIG_READY) RETURN
      CONFIG_READY=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC',VALUE,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            ENABLED=.FALSE.
          CASE DEFAULT
            ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_MINFLOP',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        READ(VALUE,*,IOSTAT=IOS) MINFLOP
        IF(IOS.NE.0) MINFLOP=1.D7
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_SYNC',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            SYNC_ENABLED=.FALSE.
          CASE DEFAULT
            SYNC_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_INITCONFIG
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE(FLOPS)
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: FLOPS
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_SHOULD_USE=ENABLED.AND.(FLOPS.GE.MINFLOP)
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_ENSURE
      IMPLICIT NONE
      INTEGER(4)                :: ISTAT
      INTEGER(ACC_HANDLE_KIND)  :: STREAM
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      IF(.NOT.ENABLED) RETURN
      IF(.NOT.HANDLE_READY) THEN
        ISTAT=CUBLASCREATE(HANDLE)
        IF(ISTAT.NE.0) THEN
          CALL ERROR$MSG('CUBLASCREATE FAILED')
          CALL ERROR$I4VAL('STATUS',ISTAT)
          CALL ERROR$STOP('CPPAW_CUBLAS_ACC_ENSURE')
        END IF
        HANDLE_READY=.TRUE.
      END IF
      STREAM=ACC_GET_CUDA_STREAM(ACC_ASYNC_SYNC)
      ISTAT=CUBLASSETSTREAM(HANDLE,STREAM)
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUBLASSETSTREAM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUBLAS_ACC_ENSURE')
      END IF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ENSURE
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_FINISH(ISTAT)
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: ISTAT
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      IF(SYNC_ENABLED) THEN
        ISTAT=CUDADEVICESYNCHRONIZE()
      ELSE
        ISTAT=0
      END IF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_FINISH
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NN_COPY(N,M,L,A,B,C,USED)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      INTEGER(4),INTENT(IN)    :: M
      INTEGER(4),INTENT(IN)    :: L
      COMPLEX(8),INTENT(IN)    :: A(N,M)
      COMPLEX(8),INTENT(IN)    :: B(M,L)
      COMPLEX(8),INTENT(INOUT) :: C(N,L)
      LOGICAL(4),INTENT(OUT)   :: USED
      REAL(8)                  :: FLOPS
!     **************************************************************************
      FLOPS=8.D0*REAL(N,KIND=8)*REAL(M,KIND=8)*REAL(L,KIND=8)
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE(FLOPS)
      IF(.NOT.USED) RETURN
!$ACC DATA COPYIN(A(1:N,1:M),B(1:M,1:L)) COPY(C(1:N,1:L))
      CALL CPPAW_CUBLAS_ACC_ZGEMM_NN_PRESENT(N,M,L,A,B,C)
!$ACC END DATA
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NN_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NN_PRESENT(N,M,L,A,B,C)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      INTEGER(4),INTENT(IN)    :: M
      INTEGER(4),INTENT(IN)    :: L
      COMPLEX(8),INTENT(IN)    :: A(N,M)
      COMPLEX(8),INTENT(IN)    :: B(M,L)
      COMPLEX(8),INTENT(INOUT) :: C(N,L)
      COMPLEX(8)               :: ONE
      INTEGER(4)               :: ISTAT
!     **************************************************************************
      ONE=(1.D0,0.D0)
      CALL CPPAW_CUBLAS_ACC_ENSURE
!$ACC HOST_DATA USE_DEVICE(A,B,C)
      ISTAT=CUBLASZGEMM(HANDLE,CUBLAS_OP_N,CUBLAS_OP_N,N,L,M,ONE &
     &                 ,A,N,B,M,ONE,C,N)
!$ACC END HOST_DATA
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUBLASZGEMM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUBLAS_ACC_ZGEMM_NN_PRESENT')
      END IF
      CALL CPPAW_CUBLAS_ACC_FINISH(ISTAT)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NN_PRESENT
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_MATMUL_COPY(N,M,L,A,B,C,USED)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      INTEGER(4),INTENT(IN)    :: M
      INTEGER(4),INTENT(IN)    :: L
      COMPLEX(8),INTENT(IN)    :: A(N,M)
      COMPLEX(8),INTENT(IN)    :: B(M,L)
      COMPLEX(8),INTENT(OUT)   :: C(N,L)
      LOGICAL(4),INTENT(OUT)   :: USED
      REAL(8)                  :: FLOPS
!     **************************************************************************
      FLOPS=8.D0*REAL(N,KIND=8)*REAL(M,KIND=8)*REAL(L,KIND=8)
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE(FLOPS)
      IF(.NOT.USED) RETURN
!$ACC DATA COPYIN(A(1:N,1:M),B(1:M,1:L)) COPYOUT(C(1:N,1:L))
      CALL CPPAW_CUBLAS_ACC_ZGEMM_MATMUL_PRESENT(N,M,L,A,B,C)
!$ACC END DATA
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_MATMUL_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_MATMUL_PRESENT(N,M,L,A,B,C)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      INTEGER(4),INTENT(IN)    :: M
      INTEGER(4),INTENT(IN)    :: L
      COMPLEX(8),INTENT(IN)    :: A(N,M)
      COMPLEX(8),INTENT(IN)    :: B(M,L)
      COMPLEX(8),INTENT(OUT)   :: C(N,L)
      COMPLEX(8)               :: ONE
      COMPLEX(8)               :: ZERO
      INTEGER(4)               :: ISTAT
!     **************************************************************************
      ONE=(1.D0,0.D0)
      ZERO=(0.D0,0.D0)
      CALL CPPAW_CUBLAS_ACC_ENSURE
!$ACC HOST_DATA USE_DEVICE(A,B,C)
      ISTAT=CUBLASZGEMM(HANDLE,CUBLAS_OP_N,CUBLAS_OP_N,N,L,M,ONE &
     &                 ,A,N,B,M,ZERO,C,N)
!$ACC END HOST_DATA
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUBLASZGEMM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUBLAS_ACC_ZGEMM_MATMUL_PRESENT')
      END IF
      CALL CPPAW_CUBLAS_ACC_FINISH(ISTAT)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_MATMUL_PRESENT
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NC_COPY(LEN1,LEN2,N,PSI1 &
     &                                         ,PSI2,OPERATOR,USED)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: LEN1
      INTEGER(4),INTENT(IN)    :: LEN2
      INTEGER(4),INTENT(IN)    :: N
      COMPLEX(8),INTENT(IN)    :: PSI1(LEN1,N)
      COMPLEX(8),INTENT(IN)    :: PSI2(LEN2,N)
      COMPLEX(8),INTENT(OUT)   :: OPERATOR(LEN1,LEN2)
      LOGICAL(4),INTENT(OUT)   :: USED
      REAL(8)                  :: FLOPS
!     **************************************************************************
      FLOPS=8.D0*REAL(LEN1,KIND=8)*REAL(LEN2,KIND=8)*REAL(N,KIND=8)
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE(FLOPS)
      IF(.NOT.USED) RETURN
!$ACC DATA COPYIN(PSI1(1:LEN1,1:N),PSI2(1:LEN2,1:N)) &
!$ACC& COPYOUT(OPERATOR(1:LEN1,1:LEN2))
      CALL CPPAW_CUBLAS_ACC_ZGEMM_NC_PRESENT(LEN1,LEN2,N,PSI1 &
     &                                      ,PSI2,OPERATOR)
!$ACC END DATA
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NC_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NC_PRESENT(LEN1,LEN2,N,PSI1 &
     &                                            ,PSI2,OPERATOR)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: LEN1
      INTEGER(4),INTENT(IN)  :: LEN2
      INTEGER(4),INTENT(IN)  :: N
      COMPLEX(8),INTENT(IN)  :: PSI1(LEN1,N)
      COMPLEX(8),INTENT(IN)  :: PSI2(LEN2,N)
      COMPLEX(8),INTENT(OUT) :: OPERATOR(LEN1,LEN2)
      COMPLEX(8)             :: ONE
      COMPLEX(8)             :: ZERO
      INTEGER(4)             :: ISTAT
!     **************************************************************************
      ONE=(1.D0,0.D0)
      ZERO=(0.D0,0.D0)
      CALL CPPAW_CUBLAS_ACC_ENSURE
!$ACC HOST_DATA USE_DEVICE(PSI1,PSI2,OPERATOR)
      ISTAT=CUBLASZGEMM(HANDLE,CUBLAS_OP_N,CUBLAS_OP_C,LEN1,LEN2,N &
     &                 ,ONE,PSI1,LEN1,PSI2,LEN2,ZERO,OPERATOR,LEN1)
!$ACC END HOST_DATA
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUBLASZGEMM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUBLAS_ACC_ZGEMM_NC_PRESENT')
      END IF
      CALL CPPAW_CUBLAS_ACC_FINISH(ISTAT)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NC_PRESENT
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_COPY(TID,LEN,N1,PSI1 &
     &                                               ,N2,PSI2,OVERLAP,USED)
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)  :: TID
      INTEGER(4),INTENT(IN)  :: LEN
      INTEGER(4),INTENT(IN)  :: N1
      INTEGER(4),INTENT(IN)  :: N2
      COMPLEX(8),INTENT(IN)  :: PSI1(LEN,N1)
      COMPLEX(8),INTENT(IN)  :: PSI2(LEN,N2)
      COMPLEX(8),INTENT(OUT) :: OVERLAP(N1,N2)
      LOGICAL(4),INTENT(OUT) :: USED
      REAL(8)                :: FLOPS
!     **************************************************************************
      IF(TID) THEN
        FLOPS=4.D0*REAL(N1,KIND=8)*REAL(N1,KIND=8)*REAL(LEN,KIND=8)
      ELSE
        FLOPS=8.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)*REAL(LEN,KIND=8)
      END IF
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE(FLOPS)
      IF(.NOT.USED) RETURN
!$ACC DATA COPYIN(PSI1(1:LEN,1:N1),PSI2(1:LEN,1:N2)) &
!$ACC& COPYOUT(OVERLAP(1:N1,1:N2))
      CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT(TID,LEN,N1,PSI1 &
     &                                           ,N2,PSI2,OVERLAP)
!$ACC END DATA
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT(TID,LEN,N1,PSI1 &
     &                                                  ,N2,PSI2,OVERLAP)
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)  :: TID
      INTEGER(4),INTENT(IN)  :: LEN
      INTEGER(4),INTENT(IN)  :: N1
      INTEGER(4),INTENT(IN)  :: N2
      COMPLEX(8),INTENT(IN)  :: PSI1(LEN,N1)
      COMPLEX(8),INTENT(IN)  :: PSI2(LEN,N2)
      COMPLEX(8),INTENT(OUT) :: OVERLAP(N1,N2)
      COMPLEX(8)             :: ONE
      COMPLEX(8)             :: ZERO
      INTEGER(4)             :: ISTAT
      INTEGER(4)             :: I,J
!     **************************************************************************
      ONE=(1.D0,0.D0)
      ZERO=(0.D0,0.D0)
      CALL CPPAW_CUBLAS_ACC_ENSURE
      IF(TID) THEN
!$ACC HOST_DATA USE_DEVICE(PSI1,OVERLAP)
        ISTAT=CUBLASZHERK(HANDLE,CUBLAS_FILL_MODE_UPPER,CUBLAS_OP_C &
     &                   ,N1,LEN,1.D0,PSI1,LEN,0.D0,OVERLAP,N1)
!$ACC END HOST_DATA
        IF(ISTAT.NE.0) THEN
          CALL ERROR$MSG('CUBLASZHERK FAILED')
          CALL ERROR$I4VAL('STATUS',ISTAT)
          CALL ERROR$STOP('CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT')
        END IF
!$ACC PARALLEL LOOP PRESENT(OVERLAP)
        DO I=1,N1
          DO J=I+1,N1
            OVERLAP(J,I)=CONJG(OVERLAP(I,J))
          ENDDO
        ENDDO
!$ACC END PARALLEL LOOP
      ELSE
!$ACC HOST_DATA USE_DEVICE(PSI1,PSI2,OVERLAP)
        ISTAT=CUBLASZGEMM(HANDLE,CUBLAS_OP_C,CUBLAS_OP_N,N1,N2,LEN &
     &                   ,ONE,PSI1,LEN,PSI2,LEN,ZERO,OVERLAP,N1)
!$ACC END HOST_DATA
        IF(ISTAT.NE.0) THEN
          CALL ERROR$MSG('CUBLASZGEMM FAILED')
          CALL ERROR$I4VAL('STATUS',ISTAT)
          CALL ERROR$STOP('CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT')
        END IF
      END IF
      CALL CPPAW_CUBLAS_ACC_FINISH(ISTAT)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT
!
      END MODULE CPPAW_CUBLAS_ACC_MODULE
#ELSE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE CPPAW_CUBLAS_ACC_STUB_MODULE
!     **************************************************************************
!     **  Keeps CPU-only builds warning-free when cuBLAS support is disabled. **
!     **************************************************************************
      IMPLICIT NONE
      END MODULE CPPAW_CUBLAS_ACC_STUB_MODULE
#ENDIF
