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
      LOGICAL(4)             :: CHECK_ENABLED=.FALSE.
      INTEGER(4)             :: MIN_N=256
      REAL(8)                :: CHECK_TOL=1.D-7
      CONTAINS
!
!     ..........................................................................
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      SUBROUTINE CPPAW_CUSOLVER_ACC_PROFILE_BYTES(NAME,N1,N2,N3,N4 &
     &                                           ,BYTES)
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
      END SUBROUTINE CPPAW_CUSOLVER_ACC_PROFILE_BYTES
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUSOLVER_ACC_PROFILE_TIME(NAME,N1,N2,N3,N4,T0)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(IN) :: N1
      INTEGER(4)  ,INTENT(IN) :: N2
      INTEGER(4)  ,INTENT(IN) :: N3
      INTEGER(4)  ,INTENT(IN) :: N4
      REAL(8)     ,INTENT(IN) :: T0
      REAL(8)                 :: T1
!     **************************************************************************
      CALL ACCELPROFILE$NOW(T1)
      CALL ACCELPROFILE$ADD(NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                     ,INT(N3,KIND=8),INT(N4,KIND=8) &
     &                     ,0.D0,0.D0,MAX(0.D0,T1-T0))
      RETURN
      END SUBROUTINE CPPAW_CUSOLVER_ACC_PROFILE_TIME
#ENDIF
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
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            ENABLED=.FALSE.
          CASE DEFAULT
            ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUSOLVER_ACC_CHECK',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            CHECK_ENABLED=.FALSE.
          CASE DEFAULT
            CHECK_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUSOLVER_ACC_CHECK_TOL',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUSOLVER_CHECK_TOL',VALUE &
     &                               ,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        READ(VALUE,*,IOSTAT=IOS) CHECK_TOL
        IF(IOS.NE.0) CHECK_TOL=1.D-7
        IF(CHECK_TOL.LE.0.D0) CHECK_TOL=1.D-7
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUSOLVER_ACC_MIN_N',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUSOLVER_MIN_N',VALUE &
     &                               ,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        READ(VALUE,*,IOSTAT=IOS) MIN_N
        IF(IOS.NE.0) MIN_N=256
        MIN_N=MAX(1,MIN_N)
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
      LOGICAL(4) FUNCTION CPPAW_CUSOLVER_ACC_DSYEVD_OK(N,H,E,U)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(IN) :: E(N)
      REAL(8)   ,INTENT(IN) :: U(N,N)
      REAL(8)   ,ALLOCATABLE:: A(:,:)
      REAL(8)   ,ALLOCATABLE:: R(:,:)
      REAL(8)               :: RERR
      REAL(8)               :: OERR
      REAL(8)               :: SCALE
      INTEGER(4)            :: I
!     **************************************************************************
      CALL CPPAW_CUSOLVER_ACC_INITCONFIG
      CPPAW_CUSOLVER_ACC_DSYEVD_OK=.TRUE.
      IF(.NOT.CHECK_ENABLED) RETURN
      ALLOCATE(A(N,N))
      ALLOCATE(R(N,N))
      A=0.5D0*(H+TRANSPOSE(H))
      R=MATMUL(A,U)
      DO I=1,N
        R(:,I)=R(:,I)-U(:,I)*E(I)
      ENDDO
      SCALE=MAX(1.D0,MAXVAL(ABS(A)))*MAX(1.D0,MAXVAL(ABS(U))) &
     &     *REAL(MAX(1,N),KIND=8)
      RERR=MAXVAL(ABS(R))/SCALE
      R=MATMUL(TRANSPOSE(U),U)
      DO I=1,N
        R(I,I)=R(I,I)-1.D0
      ENDDO
      OERR=MAXVAL(ABS(R))
      CPPAW_CUSOLVER_ACC_DSYEVD_OK=(RERR.LE.CHECK_TOL) &
     &                       .AND.(OERR.LE.10.D0*CHECK_TOL)
      DEALLOCATE(R)
      DEALLOCATE(A)
      RETURN
      END FUNCTION CPPAW_CUSOLVER_ACC_DSYEVD_OK
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUSOLVER_ACC_ZHEEVD_OK(N,H,E,U)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(IN) :: E(N)
      COMPLEX(8),INTENT(IN) :: U(N,N)
      COMPLEX(8),ALLOCATABLE:: A(:,:)
      COMPLEX(8),ALLOCATABLE:: R(:,:)
      REAL(8)               :: RERR
      REAL(8)               :: OERR
      REAL(8)               :: SCALE
      INTEGER(4)            :: I
!     **************************************************************************
      CALL CPPAW_CUSOLVER_ACC_INITCONFIG
      CPPAW_CUSOLVER_ACC_ZHEEVD_OK=.TRUE.
      IF(.NOT.CHECK_ENABLED) RETURN
      ALLOCATE(A(N,N))
      ALLOCATE(R(N,N))
      A=0.5D0*(H+TRANSPOSE(CONJG(H)))
      R=MATMUL(A,U)
      DO I=1,N
        R(:,I)=R(:,I)-U(:,I)*E(I)
      ENDDO
      SCALE=MAX(1.D0,MAXVAL(ABS(A)))*MAX(1.D0,MAXVAL(ABS(U))) &
     &     *REAL(MAX(1,N),KIND=8)
      RERR=MAXVAL(ABS(R))/SCALE
      R=MATMUL(CONJG(TRANSPOSE(U)),U)
      DO I=1,N
        R(I,I)=R(I,I)-(1.D0,0.D0)
      ENDDO
      OERR=MAXVAL(ABS(R))
      CPPAW_CUSOLVER_ACC_ZHEEVD_OK=(RERR.LE.CHECK_TOL) &
     &                       .AND.(OERR.LE.10.D0*CHECK_TOL)
      DEALLOCATE(R)
      DEALLOCATE(A)
      RETURN
      END FUNCTION CPPAW_CUSOLVER_ACC_ZHEEVD_OK
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUSOLVER_ACC_DSYGVD_OK(N,H,S,E,U)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(IN) :: S(N,N)
      REAL(8)   ,INTENT(IN) :: E(N)
      REAL(8)   ,INTENT(IN) :: U(N,N)
      REAL(8)   ,ALLOCATABLE:: R(:,:)
      REAL(8)   ,ALLOCATABLE:: SU(:,:)
      REAL(8)               :: RERR
      REAL(8)               :: OERR
      REAL(8)               :: SCALE
      INTEGER(4)            :: I
!     **************************************************************************
      CALL CPPAW_CUSOLVER_ACC_INITCONFIG
      CPPAW_CUSOLVER_ACC_DSYGVD_OK=.TRUE.
      IF(.NOT.CHECK_ENABLED) RETURN
      ALLOCATE(R(N,N))
      ALLOCATE(SU(N,N))
      SU=MATMUL(S,U)
      R=MATMUL(H,U)
      DO I=1,N
        R(:,I)=R(:,I)-SU(:,I)*E(I)
      ENDDO
      SCALE=MAX(1.D0,MAXVAL(ABS(H)),MAXVAL(ABS(S))) &
     &     *MAX(1.D0,MAXVAL(ABS(U))) &
     &     *MAX(1.D0,MAXVAL(ABS(E)))*REAL(MAX(1,N),KIND=8)
      RERR=MAXVAL(ABS(R))/SCALE
      R=MATMUL(TRANSPOSE(U),SU)
      DO I=1,N
        R(I,I)=R(I,I)-1.D0
      ENDDO
      OERR=MAXVAL(ABS(R))
      CPPAW_CUSOLVER_ACC_DSYGVD_OK=(RERR.LE.CHECK_TOL) &
     &                       .AND.(OERR.LE.10.D0*CHECK_TOL)
      DEALLOCATE(SU)
      DEALLOCATE(R)
      RETURN
      END FUNCTION CPPAW_CUSOLVER_ACC_DSYGVD_OK
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUSOLVER_ACC_ZHEGVD_OK(N,H,S,E,VEC,SYM)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      COMPLEX(8),INTENT(IN) :: S(N,N)
      REAL(8)   ,INTENT(IN) :: E(N)
      COMPLEX(8),INTENT(IN) :: VEC(N,N)
      LOGICAL(4),INTENT(IN) :: SYM
      COMPLEX(8),ALLOCATABLE:: A(:,:)
      COMPLEX(8),ALLOCATABLE:: R(:,:)
      COMPLEX(8),ALLOCATABLE:: SV(:,:)
      REAL(8)               :: RERR
      REAL(8)               :: OERR
      REAL(8)               :: SCALE
      INTEGER(4)            :: I
!     **************************************************************************
      CALL CPPAW_CUSOLVER_ACC_INITCONFIG
      CPPAW_CUSOLVER_ACC_ZHEGVD_OK=.TRUE.
      IF(.NOT.CHECK_ENABLED) RETURN
      ALLOCATE(A(N,N))
      ALLOCATE(R(N,N))
      ALLOCATE(SV(N,N))
      IF(SYM) THEN
        A=0.5D0*(H+TRANSPOSE(CONJG(H)))
      ELSE
        A=H
      END IF
      SV=MATMUL(S,VEC)
      R=MATMUL(A,VEC)
      DO I=1,N
        R(:,I)=R(:,I)-SV(:,I)*E(I)
      ENDDO
      SCALE=MAX(1.D0,MAXVAL(ABS(A)),MAXVAL(ABS(S))) &
     &     *MAX(1.D0,MAXVAL(ABS(VEC))) &
     &     *MAX(1.D0,MAXVAL(ABS(E)))*REAL(MAX(1,N),KIND=8)
      RERR=MAXVAL(ABS(R))/SCALE
      R=MATMUL(CONJG(TRANSPOSE(VEC)),SV)
      DO I=1,N
        R(I,I)=R(I,I)-(1.D0,0.D0)
      ENDDO
      OERR=MAXVAL(ABS(R))
      CPPAW_CUSOLVER_ACC_ZHEGVD_OK=(RERR.LE.CHECK_TOL) &
     &                       .AND.(OERR.LE.10.D0*CHECK_TOL)
      DEALLOCATE(SV)
      DEALLOCATE(R)
      DEALLOCATE(A)
      RETURN
      END FUNCTION CPPAW_CUSOLVER_ACC_ZHEGVD_OK
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
#ENDIF
!     **************************************************************************
      USED=.FALSE.
      INFO=0
      IF(N.LE.0) RETURN
      IF(.NOT.CPPAW_CUSOLVER_ACC_SHOULD_USE(N)) RETURN
      U=0.5D0*(H+TRANSPOSE(H))
      CALL CPPAW_CUSOLVER_ACC_ENSURE
      LWORK=0
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC DATA COPYIN(U(1:N,1:N)) CREATE(E(1:N))
!$ACC HOST_DATA USE_DEVICE(U,E)
      ISTAT=CUSOLVERDNDSYEVD_BUFFERSIZE(HANDLE,CUSOLVER_EIG_MODE_VECTOR &
     &        ,CUBLAS_FILL_MODE_UPPER,N,U,N,E,LWORK)
!$ACC END HOST_DATA
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUSOLVER_ACC_PROFILE_TIME('ACC_SETUP_CUSOLVER_DSYEVD' &
     &                                   ,N,MAX(LWORK,0),0,0,ACCEL_T0)
#ENDIF
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUSOLVER_ACC_PROFILE_BYTES('ACC_COPY_CUSOLVER_DSYEVD' &
     &     ,N,LWORK,0,0,8.D0*(3.D0*REAL(N,KIND=8)*REAL(N,KIND=8) &
     &     +REAL(N,KIND=8)))
#ENDIF
      DEALLOCATE(WORK)
      IF(ISTAT.NE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      IF(DEVINFO.NE.0) THEN
        INFO=DEVINFO
        RETURN
      END IF
      IF(.NOT.CPPAW_CUSOLVER_ACC_DSYEVD_OK(N,H,E,U)) RETURN
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
#ENDIF
!     **************************************************************************
      USED=.FALSE.
      INFO=0
      IF(N.LE.0) RETURN
      IF(.NOT.CPPAW_CUSOLVER_ACC_SHOULD_USE(N)) RETURN
      U=0.5D0*(H+TRANSPOSE(CONJG(H)))
      CALL CPPAW_CUSOLVER_ACC_ENSURE
      LWORK=0
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC DATA COPYIN(U(1:N,1:N)) CREATE(E(1:N))
!$ACC HOST_DATA USE_DEVICE(U,E)
      ISTAT=CUSOLVERDNZHEEVD_BUFFERSIZE(HANDLE,CUSOLVER_EIG_MODE_VECTOR &
     &        ,CUBLAS_FILL_MODE_UPPER,N,U,N,E,LWORK)
!$ACC END HOST_DATA
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUSOLVER_ACC_PROFILE_TIME('ACC_SETUP_CUSOLVER_ZHEEVD' &
     &                                   ,N,MAX(LWORK,0),0,0,ACCEL_T0)
#ENDIF
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUSOLVER_ACC_PROFILE_BYTES('ACC_COPY_CUSOLVER_ZHEEVD' &
     &     ,N,LWORK,0,0,16.D0*3.D0*REAL(N,KIND=8)*REAL(N,KIND=8) &
     &     +8.D0*REAL(N,KIND=8))
#ENDIF
      DEALLOCATE(WORK)
      IF(ISTAT.NE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      IF(DEVINFO.NE.0) THEN
        INFO=DEVINFO
        RETURN
      END IF
      IF(.NOT.CPPAW_CUSOLVER_ACC_ZHEEVD_OK(N,H,E,U)) RETURN
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
#ENDIF
!     **************************************************************************
      USED=.FALSE.
      INFO=0
      IF(N.LE.0) RETURN
      IF(.NOT.CPPAW_CUSOLVER_ACC_SHOULD_USE(N)) RETURN
      U=H
      B=S
      CALL CPPAW_CUSOLVER_ACC_ENSURE
      LWORK=0
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC DATA COPYIN(U(1:N,1:N),B(1:N,1:N)) CREATE(E(1:N))
!$ACC HOST_DATA USE_DEVICE(U,B,E)
      ISTAT=CUSOLVERDNDSYGVD_BUFFERSIZE(HANDLE,CUSOLVER_EIG_TYPE_1 &
     &        ,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_UPPER,N,U,N &
     &        ,B,N,E,LWORK)
!$ACC END HOST_DATA
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUSOLVER_ACC_PROFILE_TIME('ACC_SETUP_CUSOLVER_DSYGVD' &
     &                                   ,N,MAX(LWORK,0),0,0,ACCEL_T0)
#ENDIF
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUSOLVER_ACC_PROFILE_BYTES('ACC_COPY_CUSOLVER_DSYGVD' &
     &     ,N,LWORK,0,0,8.D0*(6.D0*REAL(N,KIND=8)*REAL(N,KIND=8) &
     &     +REAL(N,KIND=8)))
#ENDIF
      DEALLOCATE(WORK)
      IF(ISTAT.NE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      IF(DEVINFO.NE.0) THEN
        INFO=DEVINFO
        RETURN
      END IF
      IF(.NOT.CPPAW_CUSOLVER_ACC_DSYGVD_OK(N,H,S,E,U)) RETURN
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
#ENDIF
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
      LWORK=0
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC DATA COPYIN(VEC(1:N,1:N),B(1:N,1:N)) CREATE(E(1:N))
!$ACC HOST_DATA USE_DEVICE(VEC,B,E)
      ISTAT=CUSOLVERDNZHEGVD_BUFFERSIZE(HANDLE,CUSOLVER_EIG_TYPE_1 &
     &        ,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_UPPER,N,VEC,N &
     &        ,B,N,E,LWORK)
!$ACC END HOST_DATA
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUSOLVER_ACC_PROFILE_TIME('ACC_SETUP_CUSOLVER_ZHEGVD' &
     &                                   ,N,MAX(LWORK,0),0,0,ACCEL_T0)
#ENDIF
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUSOLVER_ACC_PROFILE_BYTES('ACC_COPY_CUSOLVER_ZHEGVD' &
     &     ,N,LWORK,0,0,16.D0*6.D0*REAL(N,KIND=8)*REAL(N,KIND=8) &
     &     +8.D0*REAL(N,KIND=8))
#ENDIF
      DEALLOCATE(WORK)
      IF(ISTAT.NE.0) THEN
        INFO=ISTAT
        RETURN
      END IF
      IF(DEVINFO.NE.0) THEN
        INFO=DEVINFO
        RETURN
      END IF
      IF(.NOT.CPPAW_CUSOLVER_ACC_ZHEGVD_OK(N,H,S,E,VEC,SYM)) RETURN
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
