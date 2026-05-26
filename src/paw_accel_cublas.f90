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
      LOGICAL(4)         :: RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: WAVE_OVERLAP_RESIDENT_ACTIVE=.FALSE.
      REAL(8)            :: MINFLOP=1.D7
      CONTAINS
!
!     ..........................................................................
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_BYTES(NAME,N1,N2,N3,N4,BYTES)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(IN) :: N1
      INTEGER(4)  ,INTENT(IN) :: N2
      INTEGER(4)  ,INTENT(IN) :: N3
      INTEGER(4)  ,INTENT(IN) :: N4
      REAL(8)     ,INTENT(IN) :: BYTES
!     **************************************************************************
      CALL ACCELPROFILE$ADD(NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                     ,INT(N3,KIND=8),INT(N4,KIND=8) &
     &                     ,0.D0,BYTES,0.D0)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_BYTES
#ENDIF
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
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_RESIDENCY',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_RESIDENCY' &
     &                               ,VALUE,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            RESIDENCY_ENABLED=.FALSE.
          CASE DEFAULT
            RESIDENCY_ENABLED=.TRUE.
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
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_RESIDENCY_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_RESIDENCY_ENABLED=ENABLED.AND.RESIDENCY_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_RESIDENCY_ENABLED
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_RESIDENT(ACTIVE)
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: ACTIVE
!     **************************************************************************
      WAVE_OVERLAP_RESIDENT_ACTIVE=ACTIVE
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_RESIDENT
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_WAVE_OVERLAP_RESIDENT_ACTIVE()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_WAVE_OVERLAP_RESIDENT_ACTIVE=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.WAVE_OVERLAP_RESIDENT_ACTIVE
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_WAVE_OVERLAP_RESIDENT_ACTIVE
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
      IF(RESIDENCY_ENABLED) THEN
!$ACC DATA PRESENT_OR_COPYIN(A(1:N,1:M),B(1:M,1:L)) COPY(C(1:N,1:L))
        CALL CPPAW_CUBLAS_ACC_ZGEMM_NN_PRESENT(N,M,L,A,B,C)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_ZGEMM_NN_RES' &
     &       ,N,M,L,0,16.D0*(REAL(N,KIND=8)*REAL(M,KIND=8) &
     &       +REAL(M,KIND=8)*REAL(L,KIND=8)))
#ENDIF
        RETURN
      END IF
!$ACC DATA COPYIN(A(1:N,1:M),B(1:M,1:L)) COPY(C(1:N,1:L))
      CALL CPPAW_CUBLAS_ACC_ZGEMM_NN_PRESENT(N,M,L,A,B,C)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_ZGEMM_NN' &
     &     ,N,M,L,0,16.D0*(REAL(N,KIND=8)*REAL(M,KIND=8) &
     &     +REAL(M,KIND=8)*REAL(L,KIND=8) &
     &     +2.D0*REAL(N,KIND=8)*REAL(L,KIND=8)))
#ENDIF
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_ZGEMM_MAT' &
     &     ,N,M,L,0,16.D0*(REAL(N,KIND=8)*REAL(M,KIND=8) &
     &     +REAL(M,KIND=8)*REAL(L,KIND=8) &
     &     +REAL(N,KIND=8)*REAL(L,KIND=8)))
#ENDIF
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
      SUBROUTINE CPPAW_CUBLAS_ACC_DGEMM_MATMUL_COPY(N,M,L,A,B,C,USED)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N
      INTEGER(4),INTENT(IN)  :: M
      INTEGER(4),INTENT(IN)  :: L
      REAL(8)   ,INTENT(IN)  :: A(N,M)
      REAL(8)   ,INTENT(IN)  :: B(M,L)
      REAL(8)   ,INTENT(OUT) :: C(N,L)
      LOGICAL(4),INTENT(OUT) :: USED
      REAL(8)                :: FLOPS
!     **************************************************************************
      FLOPS=2.D0*REAL(N,KIND=8)*REAL(M,KIND=8)*REAL(L,KIND=8)
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE(FLOPS)
      IF(.NOT.USED) RETURN
!$ACC DATA COPYIN(A(1:N,1:M),B(1:M,1:L)) COPYOUT(C(1:N,1:L))
      CALL CPPAW_CUBLAS_ACC_DGEMM_MATMUL_PRESENT(N,M,L,A,B,C)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_DGEMM_MAT' &
     &     ,N,M,L,0,8.D0*(REAL(N,KIND=8)*REAL(M,KIND=8) &
     &     +REAL(M,KIND=8)*REAL(L,KIND=8) &
     &     +REAL(N,KIND=8)*REAL(L,KIND=8)))
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_DGEMM_MATMUL_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_DGEMM_MATMUL_PRESENT(N,M,L,A,B,C)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N
      INTEGER(4),INTENT(IN)  :: M
      INTEGER(4),INTENT(IN)  :: L
      REAL(8)   ,INTENT(IN)  :: A(N,M)
      REAL(8)   ,INTENT(IN)  :: B(M,L)
      REAL(8)   ,INTENT(OUT) :: C(N,L)
      REAL(8)                :: ONE
      REAL(8)                :: ZERO
      INTEGER(4)             :: ISTAT
!     **************************************************************************
      ONE=1.D0
      ZERO=0.D0
      CALL CPPAW_CUBLAS_ACC_ENSURE
!$ACC HOST_DATA USE_DEVICE(A,B,C)
      ISTAT=CUBLASDGEMM(HANDLE,CUBLAS_OP_N,CUBLAS_OP_N,N,L,M,ONE &
     &                 ,A,N,B,M,ZERO,C,N)
!$ACC END HOST_DATA
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUBLASDGEMM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUBLAS_ACC_DGEMM_MATMUL_PRESENT')
      END IF
      CALL CPPAW_CUBLAS_ACC_FINISH(ISTAT)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_DGEMM_MATMUL_PRESENT
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_DGEMM_NT_COPY(LEN1,LEN2,N,PSI1 &
     &                                         ,PSI2,OPERATOR,USED)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: LEN1
      INTEGER(4),INTENT(IN)  :: LEN2
      INTEGER(4),INTENT(IN)  :: N
      REAL(8)   ,INTENT(IN)  :: PSI1(LEN1,N)
      REAL(8)   ,INTENT(IN)  :: PSI2(LEN2,N)
      REAL(8)   ,INTENT(OUT) :: OPERATOR(LEN1,LEN2)
      LOGICAL(4),INTENT(OUT) :: USED
      REAL(8)                :: FLOPS
!     **************************************************************************
      FLOPS=2.D0*REAL(LEN1,KIND=8)*REAL(LEN2,KIND=8)*REAL(N,KIND=8)
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE(FLOPS)
      IF(.NOT.USED) RETURN
!$ACC DATA COPYIN(PSI1(1:LEN1,1:N),PSI2(1:LEN2,1:N)) &
!$ACC& COPYOUT(OPERATOR(1:LEN1,1:LEN2))
      CALL CPPAW_CUBLAS_ACC_DGEMM_NT_PRESENT(LEN1,LEN2,N,PSI1 &
     &                                      ,PSI2,OPERATOR)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_DGEMM_NT' &
     &     ,LEN1,LEN2,N,0,8.D0*(REAL(LEN1,KIND=8)*REAL(N,KIND=8) &
     &     +REAL(LEN2,KIND=8)*REAL(N,KIND=8) &
     &     +REAL(LEN1,KIND=8)*REAL(LEN2,KIND=8)))
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_DGEMM_NT_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_DGEMM_NT_PRESENT(LEN1,LEN2,N,PSI1 &
     &                                            ,PSI2,OPERATOR)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: LEN1
      INTEGER(4),INTENT(IN)  :: LEN2
      INTEGER(4),INTENT(IN)  :: N
      REAL(8)   ,INTENT(IN)  :: PSI1(LEN1,N)
      REAL(8)   ,INTENT(IN)  :: PSI2(LEN2,N)
      REAL(8)   ,INTENT(OUT) :: OPERATOR(LEN1,LEN2)
      REAL(8)                :: ONE
      REAL(8)                :: ZERO
      INTEGER(4)             :: ISTAT
!     **************************************************************************
      ONE=1.D0
      ZERO=0.D0
      CALL CPPAW_CUBLAS_ACC_ENSURE
!$ACC HOST_DATA USE_DEVICE(PSI1,PSI2,OPERATOR)
      ISTAT=CUBLASDGEMM(HANDLE,CUBLAS_OP_N,CUBLAS_OP_T,LEN1,LEN2,N &
     &                 ,ONE,PSI1,LEN1,PSI2,LEN2,ZERO,OPERATOR,LEN1)
!$ACC END HOST_DATA
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUBLASDGEMM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUBLAS_ACC_DGEMM_NT_PRESENT')
      END IF
      CALL CPPAW_CUBLAS_ACC_FINISH(ISTAT)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_DGEMM_NT_PRESENT
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_ZGEMM_NC' &
     &     ,LEN1,LEN2,N,0,16.D0*(REAL(LEN1,KIND=8)*REAL(N,KIND=8) &
     &     +REAL(LEN2,KIND=8)*REAL(N,KIND=8) &
     &     +REAL(LEN1,KIND=8)*REAL(LEN2,KIND=8)))
#ENDIF
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
      SUBROUTINE CPPAW_CUBLAS_ACC_PROJECTION_PRESENT(NGL,NDIM,NB,LMNX &
     &                         ,LMNXX,IPRO,NPRO,PRO,PSI,GWEIGHT,WORK &
     &                         ,PROPSI)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NGL
      INTEGER(4),INTENT(IN)    :: NDIM
      INTEGER(4),INTENT(IN)    :: NB
      INTEGER(4),INTENT(IN)    :: LMNX
      INTEGER(4),INTENT(IN)    :: LMNXX
      INTEGER(4),INTENT(IN)    :: IPRO
      INTEGER(4),INTENT(IN)    :: NPRO
      COMPLEX(8),INTENT(IN)    :: PRO(NGL,LMNXX)
      COMPLEX(8),INTENT(IN)    :: PSI(NGL,NDIM,NB)
      REAL(8)   ,INTENT(IN)    :: GWEIGHT
      COMPLEX(8),INTENT(INOUT) :: WORK(LMNXX,NDIM*NB)
      COMPLEX(8),INTENT(OUT)   :: PROPSI(NDIM,NB,NPRO)
      COMPLEX(8)               :: ONE
      COMPLEX(8)               :: ZERO
      INTEGER(4)               :: ISTAT
      INTEGER(4)               :: IB
      INTEGER(4)               :: IDIM
      INTEGER(4)               :: LMN
      INTEGER(4)               :: ICOL
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                  :: ACCEL_T0
      REAL(8)                  :: ACCEL_T1
      REAL(8)                  :: ACCEL_FLOPS
      REAL(8)                  :: ACCEL_BYTES
#ENDIF
!     **************************************************************************
      ONE=(1.D0,0.D0)
      ZERO=(0.D0,0.D0)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC DATA PRESENT_OR_COPYIN(PRO(1:NGL,1:LMNX),PSI(1:NGL,1:NDIM,1:NB)) &
!$ACC& PRESENT(WORK(1:LMNXX,1:NDIM*NB),PROPSI(1:NDIM,1:NB,1:NPRO))
      CALL CPPAW_CUBLAS_ACC_ENSURE
!$ACC HOST_DATA USE_DEVICE(PRO,PSI,WORK)
      ISTAT=CUBLASZGEMM(HANDLE,CUBLAS_OP_C,CUBLAS_OP_N,LMNX,NDIM*NB &
     &                 ,NGL,ONE,PRO,NGL,PSI,NGL,ZERO,WORK,LMNXX)
!$ACC END HOST_DATA
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUBLASZGEMM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUBLAS_ACC_PROJECTION_PRESENT')
      END IF
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(WORK,PROPSI)
      DO IB=1,NB
        DO IDIM=1,NDIM
          DO LMN=1,LMNX
            ICOL=IDIM+(IB-1)*NDIM
            PROPSI(IDIM,IB,IPRO-1+LMN)=GWEIGHT*WORK(LMN,ICOL)
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL LOOP
      CALL CPPAW_CUBLAS_ACC_FINISH(ISTAT)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      ACCEL_FLOPS=8.D0*REAL(NGL,KIND=8)*REAL(LMNX,KIND=8) &
     &           *REAL(NDIM*NB,KIND=8)
      ACCEL_BYTES=16.D0*(REAL(NGL,KIND=8)*REAL(LMNX,KIND=8) &
     &                  +REAL(LMNX,KIND=8)*REAL(NDIM*NB,KIND=8))
      CALL ACCELPROFILE$ADD('CUBLAS_ZGEMM_PROJ_RES' &
     &     ,INT(NGL,KIND=8),INT(LMNX,KIND=8),INT(NDIM*NB,KIND=8) &
     &     ,0_8,ACCEL_FLOPS,ACCEL_BYTES,ACCEL_T1-ACCEL_T0)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_PROJ_RES' &
     &     ,NGL,LMNX,NDIM*NB,0,ACCEL_BYTES)
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_PROJECTION_PRESENT
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
      IF(RESIDENCY_ENABLED) THEN
!$ACC DATA PRESENT_OR_COPYIN(PSI1(1:LEN,1:N1),PSI2(1:LEN,1:N2)) &
!$ACC& COPYOUT(OVERLAP(1:N1,1:N2))
        CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT(TID,LEN,N1,PSI1 &
     &                                             ,N2,PSI2,OVERLAP)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_ZSPROD_RES' &
     &       ,LEN,N1,N2,0,16.D0*(REAL(LEN,KIND=8)*REAL(N1,KIND=8) &
     &       +REAL(N1,KIND=8)*REAL(N2,KIND=8)))
#ENDIF
        RETURN
      END IF
!$ACC DATA COPYIN(PSI1(1:LEN,1:N1),PSI2(1:LEN,1:N2)) &
!$ACC& COPYOUT(OVERLAP(1:N1,1:N2))
      CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT(TID,LEN,N1,PSI1 &
     &                                           ,N2,PSI2,OVERLAP)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_ZSPROD' &
     &     ,LEN,N1,N2,0,16.D0*(REAL(LEN,KIND=8)*REAL(N1,KIND=8) &
     &     +REAL(LEN,KIND=8)*REAL(N2,KIND=8) &
     &     +REAL(N1,KIND=8)*REAL(N2,KIND=8)))
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_RESIDENT_COPY(TID,LEN,N1 &
     &                                                       ,PSI1,N2,PSI2 &
     &                                                       ,OVERLAP,USED)
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
      REAL(8)                :: ACCEL_BYTES
#ENDIF
!     **************************************************************************
      IF(TID) THEN
        FLOPS=4.D0*REAL(N1,KIND=8)*REAL(N1,KIND=8)*REAL(LEN,KIND=8)
      ELSE
        FLOPS=8.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)*REAL(LEN,KIND=8)
      END IF
      USED=CPPAW_CUBLAS_ACC_WAVE_OVERLAP_RESIDENT_ACTIVE() &
     &     .AND.CPPAW_CUBLAS_ACC_SHOULD_USE(FLOPS)
      IF(.NOT.USED) RETURN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC DATA PRESENT(PSI1(1:LEN,1:N1),PSI2(1:LEN,1:N2)) &
!$ACC& COPYOUT(OVERLAP(1:N1,1:N2))
      CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT(TID,LEN,N1,PSI1 &
     &                                           ,N2,PSI2,OVERLAP)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      IF(TID) THEN
        ACCEL_BYTES=16.D0*(REAL(LEN,KIND=8)*REAL(N1,KIND=8) &
     &       +REAL(N1,KIND=8)*REAL(N2,KIND=8))
        CALL ACCELPROFILE$ADD('CUBLAS_ZHERK_OVL_RES' &
     &       ,INT(LEN,KIND=8),INT(N1,KIND=8),INT(N2,KIND=8),0_8 &
     &       ,FLOPS,ACCEL_BYTES,ACCEL_T1-ACCEL_T0)
      ELSE
        ACCEL_BYTES=16.D0*(REAL(LEN,KIND=8)*REAL(N1,KIND=8) &
     &       +REAL(LEN,KIND=8)*REAL(N2,KIND=8) &
     &       +REAL(N1,KIND=8)*REAL(N2,KIND=8))
        CALL ACCELPROFILE$ADD('CUBLAS_ZGEMM_OVL_RES' &
     &       ,INT(LEN,KIND=8),INT(N1,KIND=8),INT(N2,KIND=8),0_8 &
     &       ,FLOPS,ACCEL_BYTES,ACCEL_T1-ACCEL_T0)
      END IF
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_ZSPROD_OVL_RES' &
     &     ,LEN,N1,N2,0,16.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8))
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_RESIDENT_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_COPY(TID,LEN,N1,PSI1 &
     &                                                  ,N2,PSI2,OVERLAP &
     &                                                  ,USED)
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)  :: TID
      INTEGER(4),INTENT(IN)  :: LEN
      INTEGER(4),INTENT(IN)  :: N1
      INTEGER(4),INTENT(IN)  :: N2
      REAL(8)   ,INTENT(IN)  :: PSI1(LEN,N1)
      REAL(8)   ,INTENT(IN)  :: PSI2(LEN,N2)
      REAL(8)   ,INTENT(OUT) :: OVERLAP(N1,N2)
      LOGICAL(4),INTENT(OUT) :: USED
      REAL(8)                :: FLOPS
!     **************************************************************************
      IF(TID) THEN
        FLOPS=REAL(N1,KIND=8)*REAL(N1,KIND=8)*REAL(LEN,KIND=8)
      ELSE
        FLOPS=2.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)*REAL(LEN,KIND=8)
      END IF
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE(FLOPS)
      IF(.NOT.USED) RETURN
      IF(RESIDENCY_ENABLED) THEN
!$ACC DATA PRESENT_OR_COPYIN(PSI1(1:LEN,1:N1),PSI2(1:LEN,1:N2)) &
!$ACC& COPYOUT(OVERLAP(1:N1,1:N2))
        CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_PRESENT(TID,LEN,N1,PSI1 &
     &                                                ,N2,PSI2,OVERLAP)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_DSPROD_RES' &
     &       ,LEN,N1,N2,0,8.D0*(REAL(LEN,KIND=8)*REAL(N1,KIND=8) &
     &       +REAL(N1,KIND=8)*REAL(N2,KIND=8)))
#ENDIF
        RETURN
      END IF
!$ACC DATA COPYIN(PSI1(1:LEN,1:N1),PSI2(1:LEN,1:N2)) &
!$ACC& COPYOUT(OVERLAP(1:N1,1:N2))
      CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_PRESENT(TID,LEN,N1,PSI1 &
     &                                              ,N2,PSI2,OVERLAP)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_DSPROD' &
     &     ,LEN,N1,N2,0,8.D0*(REAL(LEN,KIND=8)*REAL(N1,KIND=8) &
     &     +REAL(LEN,KIND=8)*REAL(N2,KIND=8) &
     &     +REAL(N1,KIND=8)*REAL(N2,KIND=8)))
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_PRESENT(TID,LEN,N1,PSI1 &
     &                                                     ,N2,PSI2,OVERLAP)
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)  :: TID
      INTEGER(4),INTENT(IN)  :: LEN
      INTEGER(4),INTENT(IN)  :: N1
      INTEGER(4),INTENT(IN)  :: N2
      REAL(8)   ,INTENT(IN)  :: PSI1(LEN,N1)
      REAL(8)   ,INTENT(IN)  :: PSI2(LEN,N2)
      REAL(8)   ,INTENT(OUT) :: OVERLAP(N1,N2)
      REAL(8)                :: ONE
      REAL(8)                :: ZERO
      INTEGER(4)             :: ISTAT
      INTEGER(4)             :: I,J
!     **************************************************************************
      ONE=1.D0
      ZERO=0.D0
      CALL CPPAW_CUBLAS_ACC_ENSURE
      IF(TID) THEN
!$ACC HOST_DATA USE_DEVICE(PSI1,OVERLAP)
        ISTAT=CUBLASDSYRK(HANDLE,CUBLAS_FILL_MODE_UPPER,CUBLAS_OP_T &
     &                   ,N1,LEN,ONE,PSI1,LEN,ZERO,OVERLAP,N1)
!$ACC END HOST_DATA
        IF(ISTAT.NE.0) THEN
          CALL ERROR$MSG('CUBLASDSYRK FAILED')
          CALL ERROR$I4VAL('STATUS',ISTAT)
          CALL ERROR$STOP('CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_PRESENT')
        END IF
!$ACC PARALLEL LOOP PRESENT(OVERLAP)
        DO I=1,N1
          DO J=I+1,N1
            OVERLAP(J,I)=OVERLAP(I,J)
          ENDDO
        ENDDO
!$ACC END PARALLEL LOOP
      ELSE
!$ACC HOST_DATA USE_DEVICE(PSI1,PSI2,OVERLAP)
        ISTAT=CUBLASDGEMM(HANDLE,CUBLAS_OP_T,CUBLAS_OP_N,N1,N2,LEN &
     &                   ,ONE,PSI1,LEN,PSI2,LEN,ZERO,OVERLAP,N1)
!$ACC END HOST_DATA
        IF(ISTAT.NE.0) THEN
          CALL ERROR$MSG('CUBLASDGEMM FAILED')
          CALL ERROR$I4VAL('STATUS',ISTAT)
          CALL ERROR$STOP('CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_PRESENT')
        END IF
      END IF
      CALL CPPAW_CUBLAS_ACC_FINISH(ISTAT)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_PRESENT
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
