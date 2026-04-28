!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****  OPTIONAL NVIDIA cuSOLVER/OpenACC ACCELERATION HELPERS               *****
!****                                                                      *****
!****  This file is intentionally isolated from paw_library.f90 so the      *****
!****  generic LAPACK wrappers stay mostly toolchain neutral.               *****
!****                                                                      *****
!*******************************************************************************
!*******************************************************************************
!
#IF DEFINED(CPPVAR_CUSOLVER_ACC)
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE CPPAW_CUSOLVER_ACC_MODULE
!     **************************************************************************
!     **  cuSOLVER/OpenACC bridge for selected dense eigensolvers.           **
!     **  The module is only compiled for NVIDIA HPC SDK accelerator targets.**
!     **************************************************************************
      USE CUSOLVERDN
      USE CUBLAS_V2
      USE CUDAFOR
      USE OPENACC
      IMPLICIT NONE
      TYPE(CUSOLVERDNHANDLE) :: HANDLE
      LOGICAL(4)             :: HANDLE_READY=.FALSE.
      LOGICAL(4)             :: CONFIG_READY=.FALSE.
      LOGICAL(4)             :: ENABLED=.TRUE.
      INTEGER(4)             :: MIN_N=256
      CONTAINS
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUSOLVER_ACC_INITCONFIG
      IMPLICIT NONE
      CHARACTER(128) :: VALUE
      INTEGER(4)     :: STATUS
      INTEGER(4)     :: IOS
!     **************************************************************************
      IF(CONFIG_READY) RETURN
      CONFIG_READY=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUSOLVER_ACC',VALUE,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).EQ.0) RETURN
        SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
        CASE('0','no','NO','false','FALSE','off','OFF')
          ENABLED=.FALSE.
        CASE DEFAULT
          ENABLED=.TRUE.
        END SELECT
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUSOLVER_ACC_MIN_N',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        READ(VALUE,*,IOSTAT=IOS) MIN_N
        IF(IOS.NE.0) MIN_N=256
      END IF
      RETURN
      END SUBROUTINE CPPAW_CUSOLVER_ACC_INITCONFIG
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUSOLVER_ACC_SHOULD_USE(N)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
!     **************************************************************************
      CALL CPPAW_CUSOLVER_ACC_INITCONFIG
      CPPAW_CUSOLVER_ACC_SHOULD_USE=ENABLED.AND.(N.GE.MIN_N)
      RETURN
      END FUNCTION CPPAW_CUSOLVER_ACC_SHOULD_USE
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUSOLVER_ACC_ENSURE
      IMPLICIT NONE
      INTEGER(4)                :: ISTAT
      INTEGER(ACC_HANDLE_KIND)  :: STREAM
!     **************************************************************************
      CALL CPPAW_CUSOLVER_ACC_INITCONFIG
      IF(.NOT.ENABLED) RETURN
      IF(.NOT.HANDLE_READY) THEN
        ISTAT=CUSOLVERDNCREATE(HANDLE)
        IF(ISTAT.NE.0) THEN
          CALL ERROR$MSG('CUSOLVERDNCREATE FAILED')
          CALL ERROR$I4VAL('STATUS',ISTAT)
          CALL ERROR$STOP('CPPAW_CUSOLVER_ACC_ENSURE')
        END IF
        HANDLE_READY=.TRUE.
      END IF
      STREAM=ACC_GET_CUDA_STREAM(ACC_ASYNC_SYNC)
      ISTAT=CUSOLVERDNSETSTREAM(HANDLE,STREAM)
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUSOLVERDNSETSTREAM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUSOLVER_ACC_ENSURE')
      END IF
      RETURN
      END SUBROUTINE CPPAW_CUSOLVER_ACC_ENSURE
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUSOLVER_ACC_DSYEVD_COPY(N,H,E,U,USED,INFO)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N
      REAL(8)   ,INTENT(IN)  :: H(N,N)
      REAL(8)   ,INTENT(OUT) :: E(N)
      REAL(8)   ,INTENT(OUT) :: U(N,N)
      LOGICAL(4),INTENT(OUT) :: USED
      INTEGER(4),INTENT(OUT) :: INFO
      REAL(8)   ,ALLOCATABLE :: WORK(:)
      INTEGER(4)             :: LWORK
      INTEGER(4)             :: ISTAT
      INTEGER(4)             :: DEVINFO
!     **************************************************************************
      USED=.FALSE.
      INFO=0
      IF(N.LE.0) RETURN
      IF(.NOT.CPPAW_CUSOLVER_ACC_SHOULD_USE(N)) RETURN
      U=0.5D0*(H+TRANSPOSE(H))
      CALL CPPAW_CUSOLVER_ACC_ENSURE
!$ACC DATA COPYIN(U(1:N,1:N)) CREATE(E(1:N))
!$ACC HOST_DATA USE_DEVICE(U,E)
      ISTAT=CUSOLVERDNDSYEVD_BUFFERSIZE(HANDLE,CUSOLVER_EIG_MODE_VECTOR &
     &        ,CUBLAS_FILL_MODE_UPPER,N,U,N,E,LWORK)
!$ACC END HOST_DATA
!$ACC END DATA
      IF(ISTAT.NE.0.OR.LWORK.LE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      ALLOCATE(WORK(LWORK))
      DEVINFO=0
!$ACC DATA COPY(U(1:N,1:N)) COPYOUT(E(1:N)) CREATE(WORK(1:LWORK)) &
!$ACC& COPY(DEVINFO)
!$ACC HOST_DATA USE_DEVICE(U,E,WORK,DEVINFO)
      ISTAT=CUSOLVERDNDSYEVD(HANDLE,CUSOLVER_EIG_MODE_VECTOR &
     &        ,CUBLAS_FILL_MODE_UPPER,N,U,N,E,WORK,LWORK,DEVINFO)
!$ACC END HOST_DATA
!$ACC END DATA
      ISTAT=MAX(ISTAT,CUDADEVICESYNCHRONIZE())
      DEALLOCATE(WORK)
      IF(ISTAT.NE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      IF(DEVINFO.NE.0) THEN
        INFO=DEVINFO
        RETURN
      END IF
      USED=.TRUE.
      RETURN
      END SUBROUTINE CPPAW_CUSOLVER_ACC_DSYEVD_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUSOLVER_ACC_ZHEEVD_COPY(N,H,E,U,USED,INFO)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N
      COMPLEX(8),INTENT(IN)  :: H(N,N)
      REAL(8)   ,INTENT(OUT) :: E(N)
      COMPLEX(8),INTENT(OUT) :: U(N,N)
      LOGICAL(4),INTENT(OUT) :: USED
      INTEGER(4),INTENT(OUT) :: INFO
      COMPLEX(8),ALLOCATABLE :: WORK(:)
      INTEGER(4)             :: LWORK
      INTEGER(4)             :: ISTAT
      INTEGER(4)             :: DEVINFO
!     **************************************************************************
      USED=.FALSE.
      INFO=0
      IF(N.LE.0) RETURN
      IF(.NOT.CPPAW_CUSOLVER_ACC_SHOULD_USE(N)) RETURN
      U=0.5D0*(H+TRANSPOSE(CONJG(H)))
      CALL CPPAW_CUSOLVER_ACC_ENSURE
!$ACC DATA COPYIN(U(1:N,1:N)) CREATE(E(1:N))
!$ACC HOST_DATA USE_DEVICE(U,E)
      ISTAT=CUSOLVERDNZHEEVD_BUFFERSIZE(HANDLE,CUSOLVER_EIG_MODE_VECTOR &
     &        ,CUBLAS_FILL_MODE_UPPER,N,U,N,E,LWORK)
!$ACC END HOST_DATA
!$ACC END DATA
      IF(ISTAT.NE.0.OR.LWORK.LE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      ALLOCATE(WORK(LWORK))
      DEVINFO=0
!$ACC DATA COPY(U(1:N,1:N)) COPYOUT(E(1:N)) CREATE(WORK(1:LWORK)) &
!$ACC& COPY(DEVINFO)
!$ACC HOST_DATA USE_DEVICE(U,E,WORK,DEVINFO)
      ISTAT=CUSOLVERDNZHEEVD(HANDLE,CUSOLVER_EIG_MODE_VECTOR &
     &        ,CUBLAS_FILL_MODE_UPPER,N,U,N,E,WORK,LWORK,DEVINFO)
!$ACC END HOST_DATA
!$ACC END DATA
      ISTAT=MAX(ISTAT,CUDADEVICESYNCHRONIZE())
      DEALLOCATE(WORK)
      IF(ISTAT.NE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      IF(DEVINFO.NE.0) THEN
        INFO=DEVINFO
        RETURN
      END IF
      USED=.TRUE.
      RETURN
      END SUBROUTINE CPPAW_CUSOLVER_ACC_ZHEEVD_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUSOLVER_ACC_DSYGVD_COPY(N,H,S,E,U,USED,INFO)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N
      REAL(8)   ,INTENT(IN)  :: H(N,N)
      REAL(8)   ,INTENT(IN)  :: S(N,N)
      REAL(8)   ,INTENT(OUT) :: E(N)
      REAL(8)   ,INTENT(OUT) :: U(N,N)
      LOGICAL(4),INTENT(OUT) :: USED
      INTEGER(4),INTENT(OUT) :: INFO
      REAL(8)                :: B(N,N)
      REAL(8)   ,ALLOCATABLE :: WORK(:)
      INTEGER(4)             :: LWORK
      INTEGER(4)             :: ISTAT
      INTEGER(4)             :: DEVINFO
!     **************************************************************************
      USED=.FALSE.
      INFO=0
      IF(N.LE.0) RETURN
      IF(.NOT.CPPAW_CUSOLVER_ACC_SHOULD_USE(N)) RETURN
      U=H
      B=S
      CALL CPPAW_CUSOLVER_ACC_ENSURE
!$ACC DATA COPYIN(U(1:N,1:N),B(1:N,1:N)) CREATE(E(1:N))
!$ACC HOST_DATA USE_DEVICE(U,B,E)
      ISTAT=CUSOLVERDNDSYGVD_BUFFERSIZE(HANDLE,CUSOLVER_EIG_TYPE_1 &
     &        ,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_UPPER,N,U,N &
     &        ,B,N,E,LWORK)
!$ACC END HOST_DATA
!$ACC END DATA
      IF(ISTAT.NE.0.OR.LWORK.LE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      ALLOCATE(WORK(LWORK))
      DEVINFO=0
!$ACC DATA COPY(U(1:N,1:N),B(1:N,1:N)) COPYOUT(E(1:N)) &
!$ACC& CREATE(WORK(1:LWORK)) COPY(DEVINFO)
!$ACC HOST_DATA USE_DEVICE(U,B,E,WORK,DEVINFO)
      ISTAT=CUSOLVERDNDSYGVD(HANDLE,CUSOLVER_EIG_TYPE_1 &
     &        ,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_UPPER,N,U,N &
     &        ,B,N,E,WORK,LWORK,DEVINFO)
!$ACC END HOST_DATA
!$ACC END DATA
      ISTAT=MAX(ISTAT,CUDADEVICESYNCHRONIZE())
      DEALLOCATE(WORK)
      IF(ISTAT.NE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      IF(DEVINFO.NE.0) THEN
        INFO=DEVINFO
        RETURN
      END IF
      USED=.TRUE.
      RETURN
      END SUBROUTINE CPPAW_CUSOLVER_ACC_DSYGVD_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUSOLVER_ACC_ZHEGVD_COPY(N,H,S,E,VEC,SYM,USED,INFO)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N
      COMPLEX(8),INTENT(IN)  :: H(N,N)
      COMPLEX(8),INTENT(IN)  :: S(N,N)
      REAL(8)   ,INTENT(OUT) :: E(N)
      COMPLEX(8),INTENT(OUT) :: VEC(N,N)
      LOGICAL(4),INTENT(IN)  :: SYM
      LOGICAL(4),INTENT(OUT) :: USED
      INTEGER(4),INTENT(OUT) :: INFO
      COMPLEX(8)             :: B(N,N)
      COMPLEX(8),ALLOCATABLE :: WORK(:)
      INTEGER(4)             :: LWORK
      INTEGER(4)             :: ISTAT
      INTEGER(4)             :: DEVINFO
!     **************************************************************************
      USED=.FALSE.
      INFO=0
      IF(N.LE.0) RETURN
      IF(.NOT.CPPAW_CUSOLVER_ACC_SHOULD_USE(N)) RETURN
      IF(SYM) THEN
        VEC=0.5D0*(H+TRANSPOSE(CONJG(H)))
      ELSE
        VEC=H
      END IF
      B=S
      CALL CPPAW_CUSOLVER_ACC_ENSURE
!$ACC DATA COPYIN(VEC(1:N,1:N),B(1:N,1:N)) CREATE(E(1:N))
!$ACC HOST_DATA USE_DEVICE(VEC,B,E)
      ISTAT=CUSOLVERDNZHEGVD_BUFFERSIZE(HANDLE,CUSOLVER_EIG_TYPE_1 &
     &        ,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_UPPER,N,VEC,N &
     &        ,B,N,E,LWORK)
!$ACC END HOST_DATA
!$ACC END DATA
      IF(ISTAT.NE.0.OR.LWORK.LE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      ALLOCATE(WORK(LWORK))
      DEVINFO=0
!$ACC DATA COPY(VEC(1:N,1:N),B(1:N,1:N)) COPYOUT(E(1:N)) &
!$ACC& CREATE(WORK(1:LWORK)) COPY(DEVINFO)
!$ACC HOST_DATA USE_DEVICE(VEC,B,E,WORK,DEVINFO)
      ISTAT=CUSOLVERDNZHEGVD(HANDLE,CUSOLVER_EIG_TYPE_1 &
     &        ,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_UPPER,N,VEC,N &
     &        ,B,N,E,WORK,LWORK,DEVINFO)
!$ACC END HOST_DATA
!$ACC END DATA
      ISTAT=MAX(ISTAT,CUDADEVICESYNCHRONIZE())
      DEALLOCATE(WORK)
      IF(ISTAT.NE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      IF(DEVINFO.NE.0) THEN
        INFO=DEVINFO
        RETURN
      END IF
      USED=.TRUE.
      RETURN
      END SUBROUTINE CPPAW_CUSOLVER_ACC_ZHEGVD_COPY
!
      END MODULE CPPAW_CUSOLVER_ACC_MODULE
#ELSE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE CPPAW_CUSOLVER_ACC_STUB_MODULE
!     **************************************************************************
!     **  Keeps CPU-only builds warning-free when cuSOLVER support is disabled.**
!     **************************************************************************
      IMPLICIT NONE
      END MODULE CPPAW_CUSOLVER_ACC_STUB_MODULE
#ENDIF
