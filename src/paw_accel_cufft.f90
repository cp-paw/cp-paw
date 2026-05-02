!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****  OPTIONAL NVIDIA cuFFT/OpenACC ACCELERATION HELPERS                  *****
!****                                                                      *****
!****  This file is intentionally isolated from paw_library.f90 so the      *****
!****  generic FFT wrappers stay mostly toolchain neutral.                  *****
!****                                                                      *****
!*******************************************************************************
!*******************************************************************************
!
#IF DEFINED(CPPVAR_CUFFT_ACC)
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE CPPAW_CUFFT_ACC_MODULE
!     **************************************************************************
!     **  cuFFT/OpenACC bridge for complex batched 1-D FFT kernels.          **
!     **  The module is only compiled for NVIDIA HPC SDK accelerator targets.**
!     **************************************************************************
      USE CUFFT
      USE OPENACC
      IMPLICIT NONE
      TYPE CPPAW_CUFFT_ACC_PLAN_TYPE
        INTEGER(4) :: LEN=0
        INTEGER(4) :: NFFT=0
        INTEGER(4) :: PLAN=0
      END TYPE CPPAW_CUFFT_ACC_PLAN_TYPE
      TYPE CPPAW_CUFFT_ACC_PLAN3D_TYPE
        INTEGER(4) :: N1=0
        INTEGER(4) :: N2=0
        INTEGER(4) :: N3=0
        INTEGER(4) :: PLAN=0
      END TYPE CPPAW_CUFFT_ACC_PLAN3D_TYPE
      INTEGER(4),PARAMETER            :: NPLANX=512
      TYPE(CPPAW_CUFFT_ACC_PLAN_TYPE) :: PLANS(NPLANX)
      TYPE(CPPAW_CUFFT_ACC_PLAN3D_TYPE) :: PLANS3D(NPLANX)
      INTEGER(4)                      :: NPLAN=0
      INTEGER(4)                      :: NPLAN3D=0
      LOGICAL(4)                      :: CONFIG_READY=.FALSE.
      LOGICAL(4)                      :: ENABLED=.FALSE.
      LOGICAL(4)                      :: ENABLED_3D=.FALSE.
      REAL(8)                         :: MINELEMENTS=0.D0
      REAL(8)                         :: MIN3DELEMENTS=0.D0
      CONTAINS
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUFFT_ACC_INITCONFIG
      IMPLICIT NONE
      CHARACTER(128) :: VALUE
      INTEGER(4)     :: STATUS
      INTEGER(4)     :: IOS
!     **************************************************************************
      IF(CONFIG_READY) RETURN
      CONFIG_READY=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUFFT_ACC',VALUE,STATUS=STATUS)
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
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUFFT_ACC_MIN_ELEMENTS',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        READ(VALUE,*,IOSTAT=IOS) MINELEMENTS
        IF(IOS.NE.0) MINELEMENTS=0.D0
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUFFT_ACC_3D',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            ENABLED_3D=.FALSE.
          CASE DEFAULT
            ENABLED_3D=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUFFT_ACC_3D_MIN_ELEMENTS' &
     &                             ,VALUE,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        READ(VALUE,*,IOSTAT=IOS) MIN3DELEMENTS
        IF(IOS.NE.0) MIN3DELEMENTS=0.D0
      END IF
      RETURN
      END SUBROUTINE CPPAW_CUFFT_ACC_INITCONFIG
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUFFT_ACC_SHOULD_USE(LEN,NFFT)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN
      INTEGER(4),INTENT(IN) :: NFFT
      REAL(8)               :: ELEMENTS
!     **************************************************************************
      CALL CPPAW_CUFFT_ACC_INITCONFIG
      ELEMENTS=REAL(LEN,KIND=8)*REAL(NFFT,KIND=8)
      CPPAW_CUFFT_ACC_SHOULD_USE=ENABLED.AND.(ELEMENTS.GE.MINELEMENTS)
      RETURN
      END FUNCTION CPPAW_CUFFT_ACC_SHOULD_USE
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUFFT_ACC_SHOULD_USE_3D(N1,N2,N3)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      INTEGER(4),INTENT(IN) :: N3
      REAL(8)               :: ELEMENTS
!     **************************************************************************
      CALL CPPAW_CUFFT_ACC_INITCONFIG
      ELEMENTS=REAL(N1,KIND=8)*REAL(N2,KIND=8)*REAL(N3,KIND=8)
      CPPAW_CUFFT_ACC_SHOULD_USE_3D=ENABLED.AND.ENABLED_3D &
     &                             .AND.(ELEMENTS.GE.MIN3DELEMENTS)
      RETURN
      END FUNCTION CPPAW_CUFFT_ACC_SHOULD_USE_3D
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUFFT_ACC_GETPLAN(LEN,NFFT,PLAN,FOUND)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: LEN
      INTEGER(4),INTENT(IN)  :: NFFT
      INTEGER(4),INTENT(OUT) :: PLAN
      LOGICAL(4),INTENT(OUT) :: FOUND
      INTEGER(4)             :: I
      INTEGER(4)             :: ISTAT
!     **************************************************************************
      FOUND=.TRUE.
      DO I=1,NPLAN
        IF(PLANS(I)%LEN.EQ.LEN.AND.PLANS(I)%NFFT.EQ.NFFT) THEN
          PLAN=PLANS(I)%PLAN
          RETURN
        END IF
      ENDDO
      IF(NPLAN.GE.NPLANX) THEN
        FOUND=.FALSE.
        PLAN=0
        RETURN
      END IF
      NPLAN=NPLAN+1
      ISTAT=CUFFTPLAN1D(PLANS(NPLAN)%PLAN,LEN,CUFFT_Z2Z,NFFT)
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUFFTPLAN1D FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$I4VAL('LEN',LEN)
        CALL ERROR$I4VAL('NFFT',NFFT)
        CALL ERROR$STOP('CPPAW_CUFFT_ACC_GETPLAN')
      END IF
      PLANS(NPLAN)%LEN=LEN
      PLANS(NPLAN)%NFFT=NFFT
      PLAN=PLANS(NPLAN)%PLAN
      RETURN
      END SUBROUTINE CPPAW_CUFFT_ACC_GETPLAN
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUFFT_ACC_GETPLAN3D(N1,N2,N3,PLAN,FOUND)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N1
      INTEGER(4),INTENT(IN)  :: N2
      INTEGER(4),INTENT(IN)  :: N3
      INTEGER(4),INTENT(OUT) :: PLAN
      LOGICAL(4),INTENT(OUT) :: FOUND
      INTEGER(4)             :: I
      INTEGER(4)             :: ISTAT
!     **************************************************************************
      FOUND=.TRUE.
      DO I=1,NPLAN3D
        IF(PLANS3D(I)%N1.EQ.N1.AND.PLANS3D(I)%N2.EQ.N2 &
     &     .AND.PLANS3D(I)%N3.EQ.N3) THEN
          PLAN=PLANS3D(I)%PLAN
          RETURN
        END IF
      ENDDO
      IF(NPLAN3D.GE.NPLANX) THEN
        FOUND=.FALSE.
        PLAN=0
        RETURN
      END IF
      NPLAN3D=NPLAN3D+1
      ISTAT=CUFFTPLAN3D(PLANS3D(NPLAN3D)%PLAN,N3,N2,N1,CUFFT_Z2Z)
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUFFTPLAN3D FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$I4VAL('N1',N1)
        CALL ERROR$I4VAL('N2',N2)
        CALL ERROR$I4VAL('N3',N3)
        CALL ERROR$STOP('CPPAW_CUFFT_ACC_GETPLAN3D')
      END IF
      PLANS3D(NPLAN3D)%N1=N1
      PLANS3D(NPLAN3D)%N2=N2
      PLANS3D(NPLAN3D)%N3=N3
      PLAN=PLANS3D(NPLAN3D)%PLAN
      RETURN
      END SUBROUTINE CPPAW_CUFFT_ACC_GETPLAN3D
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUFFT_ACC_FFTC8_COPY(DIR,LEN,NFFT,X,Y,USED)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      LOGICAL(4)  ,INTENT(OUT):: USED
      INTEGER(4)              :: PLAN
      INTEGER(4)              :: DIRECTION
      INTEGER(4)              :: ISTAT
      INTEGER(4)              :: I
      INTEGER(4)              :: IFFT
      REAL(8)                 :: SCALE
      LOGICAL(4)              :: FOUND
!     **************************************************************************
      USED=.FALSE.
      IF(LEN.LE.0.OR.NFFT.LE.0) RETURN
      IF(.NOT.CPPAW_CUFFT_ACC_SHOULD_USE(LEN,NFFT)) RETURN
      IF(DIR.EQ.'RTOG') THEN
        DIRECTION=CUFFT_FORWARD
      ELSE IF(DIR.EQ.'GTOR') THEN
        DIRECTION=CUFFT_INVERSE
      ELSE
        RETURN
      END IF
      CALL CPPAW_CUFFT_ACC_GETPLAN(LEN,NFFT,PLAN,FOUND)
      IF(.NOT.FOUND) RETURN
      ISTAT=CUFFTSETSTREAM(PLAN,ACC_GET_CUDA_STREAM(ACC_ASYNC_SYNC))
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUFFTSETSTREAM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUFFT_ACC_FFTC8_COPY')
      END IF
      Y=X
!$ACC DATA COPY(Y(1:LEN,1:NFFT))
!$ACC HOST_DATA USE_DEVICE(Y)
      ISTAT=CUFFTEXECZ2Z(PLAN,Y,Y,DIRECTION)
!$ACC END HOST_DATA
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUFFTEXECZ2Z FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUFFT_ACC_FFTC8_COPY')
      END IF
      IF(DIR.EQ.'RTOG') THEN
        SCALE=1.D0/REAL(LEN,KIND=8)
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(Y)
        DO IFFT=1,NFFT
          DO I=1,LEN
            Y(I,IFFT)=Y(I,IFFT)*SCALE
          ENDDO
        ENDDO
!$ACC END PARALLEL LOOP
      END IF
!$ACC END DATA
      USED=.TRUE.
      RETURN
      END SUBROUTINE CPPAW_CUFFT_ACC_FFTC8_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUFFT_ACC_3DFFTC8_COPY(DIR,N1,N2,N3,X,Y,USED)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: N1
      INTEGER(4)  ,INTENT(IN) :: N2
      INTEGER(4)  ,INTENT(IN) :: N3
      COMPLEX(8)  ,INTENT(IN) :: X(N1,N2,N3)
      COMPLEX(8)  ,INTENT(OUT):: Y(N1,N2,N3)
      LOGICAL(4)  ,INTENT(OUT):: USED
      INTEGER(4)              :: PLAN
      INTEGER(4)              :: DIRECTION
      INTEGER(4)              :: ISTAT
      INTEGER(4)              :: I1
      INTEGER(4)              :: I2
      INTEGER(4)              :: I3
      REAL(8)                 :: SCALE
      LOGICAL(4)              :: FOUND
!     **************************************************************************
      USED=.FALSE.
      IF(N1.LE.0.OR.N2.LE.0.OR.N3.LE.0) RETURN
      IF(.NOT.CPPAW_CUFFT_ACC_SHOULD_USE_3D(N1,N2,N3)) RETURN
      IF(DIR.EQ.'RTOG') THEN
        DIRECTION=CUFFT_FORWARD
      ELSE IF(DIR.EQ.'GTOR') THEN
        DIRECTION=CUFFT_INVERSE
      ELSE
        RETURN
      END IF
      CALL CPPAW_CUFFT_ACC_GETPLAN3D(N1,N2,N3,PLAN,FOUND)
      IF(.NOT.FOUND) RETURN
      ISTAT=CUFFTSETSTREAM(PLAN,ACC_GET_CUDA_STREAM(ACC_ASYNC_SYNC))
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUFFTSETSTREAM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUFFT_ACC_3DFFTC8_COPY')
      END IF
      Y=X
!$ACC DATA COPY(Y(1:N1,1:N2,1:N3))
!$ACC HOST_DATA USE_DEVICE(Y)
      ISTAT=CUFFTEXECZ2Z(PLAN,Y,Y,DIRECTION)
!$ACC END HOST_DATA
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUFFTEXECZ2Z FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUFFT_ACC_3DFFTC8_COPY')
      END IF
      IF(DIR.EQ.'RTOG') THEN
        SCALE=1.D0/(REAL(N1,KIND=8)*REAL(N2,KIND=8)*REAL(N3,KIND=8))
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(Y)
        DO I3=1,N3
          DO I2=1,N2
            DO I1=1,N1
              Y(I1,I2,I3)=Y(I1,I2,I3)*SCALE
            ENDDO
          ENDDO
        ENDDO
!$ACC END PARALLEL LOOP
      END IF
!$ACC END DATA
      USED=.TRUE.
      RETURN
      END SUBROUTINE CPPAW_CUFFT_ACC_3DFFTC8_COPY
!
      END MODULE CPPAW_CUFFT_ACC_MODULE
#ELSE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      MODULE CPPAW_CUFFT_ACC_STUB_MODULE
!     **************************************************************************
!     **  Keeps CPU-only builds warning-free when cuFFT support is disabled.  **
!     **************************************************************************
      IMPLICIT NONE
      END MODULE CPPAW_CUFFT_ACC_STUB_MODULE
#ENDIF
