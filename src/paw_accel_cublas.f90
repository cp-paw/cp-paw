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
#IF DEFINED(CPPVAR_GPU_RESIDENCY_PROFILE)
      LOGICAL(4)         :: RESIDENCY_ENABLED=.TRUE.
      LOGICAL(4)         :: PRO_EXPANSION_ENABLED=.TRUE.
      LOGICAL(4)         :: ADDPRO_CACHE_ENABLED=.TRUE.
      LOGICAL(4)         :: PROJ_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: ADDPRO_CACHE_HPSI_ENABLED=.TRUE.
      LOGICAL(4)         :: ADDPRO_CACHE_OPSI_ENABLED=.TRUE.
      LOGICAL(4)         :: FORCE_PSI_RESIDENCY_ENABLED=.TRUE.
      LOGICAL(4)         :: HPSI_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: PSIM_PROPAGATE_ENABLED=.FALSE.
      LOGICAL(4)         :: PSIM_PHASE_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: OPSI_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: ORTHO_CONST_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: ORTHO_X_RESIDENCY_ENABLED=.TRUE.
      LOGICAL(4)         :: ONECENTER_OVERLAP_ENABLED=.TRUE.
      LOGICAL(4)         :: DENMAT_ENERGY_ENABLED=.FALSE.
#ELSE
      LOGICAL(4)         :: RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: PRO_EXPANSION_ENABLED=.FALSE.
      LOGICAL(4)         :: ADDPRO_CACHE_ENABLED=.FALSE.
      LOGICAL(4)         :: PROJ_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: ADDPRO_CACHE_HPSI_ENABLED=.FALSE.
      LOGICAL(4)         :: ADDPRO_CACHE_OPSI_ENABLED=.FALSE.
      LOGICAL(4)         :: FORCE_PSI_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: HPSI_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: PSIM_PROPAGATE_ENABLED=.FALSE.
      LOGICAL(4)         :: PSIM_PHASE_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: OPSI_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: ORTHO_CONST_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: ORTHO_X_RESIDENCY_ENABLED=.FALSE.
      LOGICAL(4)         :: ONECENTER_OVERLAP_ENABLED=.FALSE.
      LOGICAL(4)         :: DENMAT_ENERGY_ENABLED=.FALSE.
#ENDIF
      LOGICAL(4)         :: INVERSION_BATCH_ENABLED=.TRUE.
      LOGICAL(4)         :: WAVE_OVERLAP_RESIDENT_ACTIVE=.FALSE.
      LOGICAL(4)         :: WAVE_OVERLAP_FORCE_ACTIVE=.FALSE.
      REAL(8)            :: MINFLOP=1.D7
      REAL(8)            :: MINFLOP_PROJECTION=1.D7
      REAL(8)            :: MINFLOP_OVERLAP=1.D7
      REAL(8)            :: MINFLOP_ADDPRODUCT=1.D7
      REAL(8)            :: MINFLOP_MATMUL=1.D7
      REAL(8)            :: MINFLOP_DENMAT=1.D8
      REAL(8)            :: MINFLOP_OFFDEN=1.D8
      CHARACTER(16)      :: OVERLAP_PROFILE_ID=''
      CONTAINS
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_SET_OVERLAP_PROFILE(PROFILE_ID)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: PROFILE_ID
!     **************************************************************************
      OVERLAP_PROFILE_ID=ADJUSTL(PROFILE_ID)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_SET_OVERLAP_PROFILE
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_CLEAR_OVERLAP_PROFILE()
      IMPLICIT NONE
!     **************************************************************************
      OVERLAP_PROFILE_ID=''
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_CLEAR_OVERLAP_PROFILE
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
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_ZSPROD_NAMES(SUFFIX,PRESENT_NAME &
     &                                        ,COPY_NAME)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: SUFFIX
      CHARACTER(32),INTENT(OUT):: PRESENT_NAME
      CHARACTER(32),INTENT(OUT):: COPY_NAME
!     **************************************************************************
      IF(LEN_TRIM(OVERLAP_PROFILE_ID).GT.0) THEN
        PRESENT_NAME='ACC_PRESENT_ZSP_'//TRIM(OVERLAP_PROFILE_ID)//'_' &
     &       //TRIM(SUFFIX)
        COPY_NAME='ACC_COPY_ZSP_'//TRIM(OVERLAP_PROFILE_ID)//'_' &
     &       //TRIM(SUFFIX)
        IF(TRIM(SUFFIX).EQ.'P1') COPY_NAME=TRIM(COPY_NAME)//'_IN'
        IF(TRIM(SUFFIX).EQ.'P2') COPY_NAME=TRIM(COPY_NAME)//'_IN'
      ELSE IF(TRIM(SUFFIX).EQ.'P1') THEN
        PRESENT_NAME='ACC_PRESENT_CUBLAS_ZSPROD_PSI1'
        COPY_NAME='ACC_COPY_CUBLAS_ZSPROD_PSI1_IN'
      ELSE IF(TRIM(SUFFIX).EQ.'P2') THEN
        PRESENT_NAME='ACC_PRESENT_CUBLAS_ZSPROD_PSI2'
        COPY_NAME='ACC_COPY_CUBLAS_ZSPROD_PSI2_IN'
      ELSE
        PRESENT_NAME='ACC_PRESENT_CUBLAS_ZSPROD_OUT'
        COPY_NAME='ACC_COPY_CUBLAS_ZSPROD_OUT'
      END IF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ZSPROD_NAMES
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_ZSPROD_COPY_NAME(SUFFIX,COPY_NAME)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: SUFFIX
      CHARACTER(32),INTENT(OUT) :: COPY_NAME
!     **************************************************************************
      IF(LEN_TRIM(OVERLAP_PROFILE_ID).GT.0) THEN
        COPY_NAME='ACC_COPY_ZSP_'//TRIM(OVERLAP_PROFILE_ID)//'_' &
     &       //TRIM(SUFFIX)
      ELSE IF(TRIM(SUFFIX).EQ.'OVL') THEN
        COPY_NAME='ACC_COPY_CUBLAS_ZSPROD_OVL_RES'
      ELSE
        COPY_NAME='ACC_COPY_CUBLAS_ZSPROD'
      END IF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ZSPROD_COPY_NAME
#ENDIF
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D(PRESENT_NAME &
     &                                                 ,COPY_NAME &
     &                                                 ,N1,N2,N3,ARRAY)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: PRESENT_NAME
      CHARACTER(*),INTENT(IN) :: COPY_NAME
      INTEGER(4)  ,INTENT(IN) :: N1
      INTEGER(4)  ,INTENT(IN) :: N2
      INTEGER(4)  ,INTENT(IN) :: N3
      COMPLEX(8)              :: ARRAY(N1,N2,N3)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      LOGICAL(4)              :: ISPRESENT
      REAL(8)                 :: BYTES
!     **************************************************************************
      ISPRESENT=ACC_IS_PRESENT(ARRAY)
      BYTES=16.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)*REAL(N3,KIND=8)
      IF(ISPRESENT) THEN
        CALL ACCELPROFILE$ADD(PRESENT_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,INT(N3,KIND=8),0_8,0.D0,0.D0,0.D0)
      ELSE
        CALL ACCELPROFILE$ADD(COPY_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,INT(N3,KIND=8),0_8,0.D0,BYTES,0.D0)
      END IF
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D_IO(PRESENT_NAME &
     &                                                    ,COPY_NAME &
     &                                                    ,N1,N2,N3,ARRAY)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: PRESENT_NAME
      CHARACTER(*),INTENT(IN) :: COPY_NAME
      INTEGER(4)  ,INTENT(IN) :: N1
      INTEGER(4)  ,INTENT(IN) :: N2
      INTEGER(4)  ,INTENT(IN) :: N3
      COMPLEX(8)              :: ARRAY(N1,N2,N3)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      LOGICAL(4)              :: ISPRESENT
      REAL(8)                 :: BYTES
!     **************************************************************************
      ISPRESENT=ACC_IS_PRESENT(ARRAY)
      BYTES=32.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)*REAL(N3,KIND=8)
      IF(ISPRESENT) THEN
        CALL ACCELPROFILE$ADD(PRESENT_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,INT(N3,KIND=8),0_8,0.D0,0.D0,0.D0)
      ELSE
        CALL ACCELPROFILE$ADD(COPY_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,INT(N3,KIND=8),0_8,0.D0,BYTES,0.D0)
      END IF
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D_IO
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D_INOUT(PRESENT_NAME &
     &                                                       ,COPY_IN_NAME &
     &                                                       ,COPY_OUT_NAME &
     &                                                       ,N1,N2,N3,ARRAY)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: PRESENT_NAME
      CHARACTER(*),INTENT(IN) :: COPY_IN_NAME
      CHARACTER(*),INTENT(IN) :: COPY_OUT_NAME
      INTEGER(4)  ,INTENT(IN) :: N1
      INTEGER(4)  ,INTENT(IN) :: N2
      INTEGER(4)  ,INTENT(IN) :: N3
      COMPLEX(8)              :: ARRAY(N1,N2,N3)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      LOGICAL(4)              :: ISPRESENT
      REAL(8)                 :: BYTES
!     **************************************************************************
      ISPRESENT=ACC_IS_PRESENT(ARRAY)
      BYTES=16.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)*REAL(N3,KIND=8)
      IF(ISPRESENT) THEN
        CALL ACCELPROFILE$ADD(PRESENT_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,INT(N3,KIND=8),0_8,0.D0,0.D0,0.D0)
      ELSE
        CALL ACCELPROFILE$ADD(COPY_IN_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,INT(N3,KIND=8),0_8,0.D0,BYTES,0.D0)
        CALL ACCELPROFILE$ADD(COPY_OUT_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,INT(N3,KIND=8),0_8,0.D0,BYTES,0.D0)
      END IF
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D_INOUT
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D(PRESENT_NAME &
     &                                                 ,COPY_NAME &
     &                                                 ,N1,N2,ARRAY)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: PRESENT_NAME
      CHARACTER(*),INTENT(IN) :: COPY_NAME
      INTEGER(4)  ,INTENT(IN) :: N1
      INTEGER(4)  ,INTENT(IN) :: N2
      COMPLEX(8)              :: ARRAY(N1,N2)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      LOGICAL(4)              :: ISPRESENT
      REAL(8)                 :: BYTES
!     **************************************************************************
      ISPRESENT=ACC_IS_PRESENT(ARRAY)
      BYTES=16.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)
      IF(ISPRESENT) THEN
        CALL ACCELPROFILE$ADD(PRESENT_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,0_8,0_8,0.D0,0.D0,0.D0)
      ELSE
        CALL ACCELPROFILE$ADD(COPY_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,0_8,0_8,0.D0,BYTES,0.D0)
      END IF
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D_IO(PRESENT_NAME &
     &                                                    ,COPY_NAME &
     &                                                    ,N1,N2,ARRAY)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: PRESENT_NAME
      CHARACTER(*),INTENT(IN) :: COPY_NAME
      INTEGER(4)  ,INTENT(IN) :: N1
      INTEGER(4)  ,INTENT(IN) :: N2
      COMPLEX(8)              :: ARRAY(N1,N2)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      LOGICAL(4)              :: ISPRESENT
      REAL(8)                 :: BYTES
!     **************************************************************************
      ISPRESENT=ACC_IS_PRESENT(ARRAY)
      BYTES=32.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)
      IF(ISPRESENT) THEN
        CALL ACCELPROFILE$ADD(PRESENT_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,0_8,0_8,0.D0,0.D0,0.D0)
      ELSE
        CALL ACCELPROFILE$ADD(COPY_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,0_8,0_8,0.D0,BYTES,0.D0)
      END IF
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D_IO
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_1D(PRESENT_NAME &
     &                                                 ,COPY_NAME,N1,ARRAY)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: PRESENT_NAME
      CHARACTER(*),INTENT(IN) :: COPY_NAME
      INTEGER(4)  ,INTENT(IN) :: N1
      COMPLEX(8)              :: ARRAY(N1)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      LOGICAL(4)              :: ISPRESENT
      REAL(8)                 :: BYTES
!     **************************************************************************
      ISPRESENT=ACC_IS_PRESENT(ARRAY)
      BYTES=16.D0*REAL(N1,KIND=8)
      IF(ISPRESENT) THEN
        CALL ACCELPROFILE$ADD(PRESENT_NAME,INT(N1,KIND=8),0_8,0_8,0_8 &
     &                       ,0.D0,0.D0,0.D0)
      ELSE
        CALL ACCELPROFILE$ADD(COPY_NAME,INT(N1,KIND=8),0_8,0_8,0_8 &
     &                       ,0.D0,BYTES,0.D0)
      END IF
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_1D
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_R8_2D(PRESENT_NAME &
     &                                                ,COPY_NAME,N1,N2,ARRAY)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: PRESENT_NAME
      CHARACTER(*),INTENT(IN) :: COPY_NAME
      INTEGER(4)  ,INTENT(IN) :: N1
      INTEGER(4)  ,INTENT(IN) :: N2
      REAL(8)                 :: ARRAY(N1,N2)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      LOGICAL(4)              :: ISPRESENT
      REAL(8)                 :: BYTES
!     **************************************************************************
      ISPRESENT=ACC_IS_PRESENT(ARRAY)
      BYTES=8.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)
      IF(ISPRESENT) THEN
        CALL ACCELPROFILE$ADD(PRESENT_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,0_8,0_8,0.D0,0.D0,0.D0)
      ELSE
        CALL ACCELPROFILE$ADD(COPY_NAME,INT(N1,KIND=8),INT(N2,KIND=8) &
     &                       ,0_8,0_8,0.D0,BYTES,0.D0)
      END IF
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_PROFILE_PRESENT_R8_2D
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_READ_REAL_ENV(NAME,VALUE)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: NAME
      REAL(8)     ,INTENT(INOUT) :: VALUE
      CHARACTER(128) :: TEXT
      REAL(8)        :: TMP
      INTEGER(4)     :: STATUS
      INTEGER(4)     :: IOS
!     **************************************************************************
      CALL GET_ENVIRONMENT_VARIABLE(NAME,TEXT,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        READ(TEXT,*,IOSTAT=IOS) TMP
        IF(IOS.EQ.0) VALUE=TMP
      END IF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_READ_REAL_ENV
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV(NAME,VALUE,FOUND)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: NAME
      LOGICAL(4)  ,INTENT(INOUT) :: VALUE
      LOGICAL(4)  ,INTENT(OUT)   :: FOUND
      CHARACTER(128)             :: TEXT
      INTEGER(4)                 :: STATUS
!     **************************************************************************
      FOUND=.FALSE.
      CALL GET_ENVIRONMENT_VARIABLE(NAME,TEXT,STATUS=STATUS)
      IF(STATUS.NE.0) RETURN
      TEXT=ADJUSTL(TEXT)
      IF(LEN_TRIM(TEXT).EQ.0) RETURN
      FOUND=.TRUE.
      SELECT CASE(TEXT(1:MIN(LEN(TEXT),LEN_TRIM(TEXT))))
      CASE('0','no','NO','false','FALSE','off','OFF')
        VALUE=.FALSE.
      CASE DEFAULT
        VALUE=.TRUE.
      END SELECT
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_INITCONFIG
      IMPLICIT NONE
      CHARACTER(128) :: VALUE
      INTEGER(4)     :: STATUS
      INTEGER(4)     :: IOS
      LOGICAL(4)     :: FOUND
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
      MINFLOP_PROJECTION=MINFLOP
      MINFLOP_OVERLAP=MINFLOP
      MINFLOP_ADDPRODUCT=MINFLOP
      MINFLOP_MATMUL=MINFLOP
      CALL CPPAW_CUBLAS_ACC_READ_REAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_PROJECTION_MINFLOP',MINFLOP_PROJECTION)
      CALL CPPAW_CUBLAS_ACC_READ_REAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_OVERLAP_MINFLOP',MINFLOP_OVERLAP)
      CALL CPPAW_CUBLAS_ACC_READ_REAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_ADDPRODUCT_MINFLOP',MINFLOP_ADDPRODUCT)
      CALL CPPAW_CUBLAS_ACC_READ_REAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_MATMUL_MINFLOP',MINFLOP_MATMUL)
      CALL CPPAW_CUBLAS_ACC_READ_REAL_ENV &
     &    ('CPPAW_GPU_DENMAT_MINFLOP',MINFLOP_DENMAT)
      CALL CPPAW_CUBLAS_ACC_READ_REAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_DENMAT_MINFLOP',MINFLOP_DENMAT)
      CALL CPPAW_CUBLAS_ACC_READ_REAL_ENV &
     &    ('CPPAW_GPU_OFFDEN_MINFLOP',MINFLOP_OFFDEN)
      CALL CPPAW_CUBLAS_ACC_READ_REAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_OFFDEN_MINFLOP',MINFLOP_OFFDEN)
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
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_PRO_EXPANSION',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_PRO_EXPANSION' &
     &                               ,VALUE,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            PRO_EXPANSION_ENABLED=.FALSE.
          CASE DEFAULT
            PRO_EXPANSION_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_ADDPRO_CACHE',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_ADDPRO_CACHE' &
     &                               ,VALUE,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            ADDPRO_CACHE_ENABLED=.FALSE.
          CASE DEFAULT
            ADDPRO_CACHE_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_GPU_PROJ_RESIDENCY' &
     &    ,PROJ_RESIDENCY_ENABLED,FOUND)
      IF(.NOT.FOUND) CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_PROJ_RESIDENCY' &
     &    ,PROJ_RESIDENCY_ENABLED,FOUND)
      CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_GPU_ADDPRO_CACHE_HPSI' &
     &    ,ADDPRO_CACHE_HPSI_ENABLED,FOUND)
      IF(.NOT.FOUND) CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_ADDPRO_CACHE_HPSI' &
     &    ,ADDPRO_CACHE_HPSI_ENABLED,FOUND)
      CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_GPU_ADDPRO_CACHE_OPSI' &
     &    ,ADDPRO_CACHE_OPSI_ENABLED,FOUND)
      IF(.NOT.FOUND) CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_ADDPRO_CACHE_OPSI' &
     &    ,ADDPRO_CACHE_OPSI_ENABLED,FOUND)
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_FORCE_PSI_RESIDENCY',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_FORCE_PSI_RESIDENCY' &
     &                               ,VALUE,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            FORCE_PSI_RESIDENCY_ENABLED=.FALSE.
          CASE DEFAULT
            FORCE_PSI_RESIDENCY_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_HPSI_RESIDENCY',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_HPSI_RESIDENCY' &
     &                               ,VALUE,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            HPSI_RESIDENCY_ENABLED=.FALSE.
          CASE DEFAULT
            HPSI_RESIDENCY_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_OPSI_RESIDENCY',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_OPSI_RESIDENCY' &
     &                               ,VALUE,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            OPSI_RESIDENCY_ENABLED=.FALSE.
          CASE DEFAULT
            OPSI_RESIDENCY_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_GPU_PSIM_PROPAGATE' &
     &    ,PSIM_PROPAGATE_ENABLED,FOUND)
      IF(.NOT.FOUND) CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_GPU_PSIM_RESIDENCY' &
     &    ,PSIM_PROPAGATE_ENABLED,FOUND)
      IF(.NOT.FOUND) CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_PSIM_RESIDENCY' &
     &    ,PSIM_PROPAGATE_ENABLED,FOUND)
      CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_GPU_PSIM_PHASE_RESIDENCY' &
     &    ,PSIM_PHASE_RESIDENCY_ENABLED,FOUND)
      IF(.NOT.FOUND) CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_GPU_PSIM_KEEP_RESIDENT' &
     &    ,PSIM_PHASE_RESIDENCY_ENABLED,FOUND)
      IF(.NOT.FOUND) CALL CPPAW_CUBLAS_ACC_READ_LOGICAL_ENV &
     &    ('CPPAW_CUBLAS_ACC_PSIM_PHASE_RESIDENCY' &
     &    ,PSIM_PHASE_RESIDENCY_ENABLED,FOUND)
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_ORTHO_CONST_RESIDENCY',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_ORTHO_CONST_RESIDENCY' &
     &                               ,VALUE,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            ORTHO_CONST_RESIDENCY_ENABLED=.FALSE.
          CASE DEFAULT
            ORTHO_CONST_RESIDENCY_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_ORTHO_X_RESIDENCY',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_ORTHO_X_RESIDENCY' &
     &                               ,VALUE,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            ORTHO_X_RESIDENCY_ENABLED=.FALSE.
          CASE DEFAULT
            ORTHO_X_RESIDENCY_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_1COVERLAP',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_1COVERLAP' &
     &                               ,VALUE,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            ONECENTER_OVERLAP_ENABLED=.FALSE.
          CASE DEFAULT
            ONECENTER_OVERLAP_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_DENMAT_ENERGY',VALUE &
     &                             ,STATUS=STATUS)
      IF(STATUS.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_DENMAT_ENERGY' &
     &                               ,VALUE,STATUS=STATUS)
      END IF
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            DENMAT_ENERGY_ENABLED=.FALSE.
          CASE DEFAULT
            DENMAT_ENERGY_ENABLED=.TRUE.
          END SELECT
        END IF
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_INVERSION_BATCH' &
     &                             ,VALUE,STATUS=STATUS)
      IF(STATUS.EQ.0) THEN
        VALUE=ADJUSTL(VALUE)
        IF(LEN_TRIM(VALUE).GT.0) THEN
          SELECT CASE(VALUE(1:MIN(LEN(VALUE),LEN_TRIM(VALUE))))
          CASE('0','no','NO','false','FALSE','off','OFF')
            INVERSION_BATCH_ENABLED=.FALSE.
          CASE DEFAULT
            INVERSION_BATCH_ENABLED=.TRUE.
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
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE_PROJECTION(FLOPS)
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: FLOPS
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_SHOULD_USE_PROJECTION=ENABLED &
     &                                      .AND.(FLOPS.GE.MINFLOP_PROJECTION)
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE_PROJECTION
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE_OVERLAP(FLOPS)
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: FLOPS
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_SHOULD_USE_OVERLAP=ENABLED &
     &                                   .AND.(FLOPS.GE.MINFLOP_OVERLAP)
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE_OVERLAP
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE_ADDPRODUCT(FLOPS)
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: FLOPS
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_SHOULD_USE_ADDPRODUCT=ENABLED &
     &                                      .AND.(FLOPS.GE.MINFLOP_ADDPRODUCT)
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE_ADDPRODUCT
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE_MATMUL(FLOPS)
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: FLOPS
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_SHOULD_USE_MATMUL=ENABLED &
     &                                  .AND.(FLOPS.GE.MINFLOP_MATMUL)
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE_MATMUL
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE_OFFDEN(FLOPS)
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: FLOPS
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_SHOULD_USE_OFFDEN=ENABLED &
     &                                  .AND.(FLOPS.GE.MINFLOP_OFFDEN)
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_SHOULD_USE_OFFDEN
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
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_PRO_EXPANSION_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_PRO_EXPANSION_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.PRO_EXPANSION_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_PRO_EXPANSION_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_ADDPRO_CACHE_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_ADDPRO_CACHE_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.PRO_EXPANSION_ENABLED &
     &     .AND.ADDPRO_CACHE_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_ADDPRO_CACHE_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_ADDPRO_CACHE_CONTEXT_ENABLED &
     &                                    (PROFILE_ID)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: PROFILE_ID
      LOGICAL(4)              :: TCONTEXT
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      TCONTEXT=.TRUE.
      IF(INDEX(TRIM(PROFILE_ID),'HPSI').GT.0) THEN
        TCONTEXT=ADDPRO_CACHE_HPSI_ENABLED
      ELSE IF(INDEX(TRIM(PROFILE_ID),'OPSI').GT.0) THEN
        TCONTEXT=ADDPRO_CACHE_OPSI_ENABLED
      END IF
      CPPAW_CUBLAS_ACC_ADDPRO_CACHE_CONTEXT_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.PRO_EXPANSION_ENABLED &
     &     .AND.ADDPRO_CACHE_ENABLED.AND.TCONTEXT
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_ADDPRO_CACHE_CONTEXT_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_PROJ_RESIDENCY_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_PROJ_RESIDENCY_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.PROJ_RESIDENCY_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_PROJ_RESIDENCY_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_FORCE_PSI_RESIDENCY_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_FORCE_PSI_RESIDENCY_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.FORCE_PSI_RESIDENCY_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_FORCE_PSI_RESIDENCY_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_HPSI_RESIDENCY_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_HPSI_RESIDENCY_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.HPSI_RESIDENCY_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_HPSI_RESIDENCY_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_PSIM_PROPAGATE_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_PSIM_PROPAGATE_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.PSIM_PROPAGATE_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_PSIM_PROPAGATE_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_PSIM_PHASE_RESIDENCY_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_PSIM_PHASE_RESIDENCY_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.PSIM_PROPAGATE_ENABLED &
     &     .AND.PSIM_PHASE_RESIDENCY_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_PSIM_PHASE_RESIDENCY_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_OPSI_RESIDENCY_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_OPSI_RESIDENCY_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.OPSI_RESIDENCY_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_OPSI_RESIDENCY_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_ORTHO_CONST_RESIDENCY_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_ORTHO_CONST_RESIDENCY_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.ORTHO_CONST_RESIDENCY_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_ORTHO_CONST_RESIDENCY_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_ORTHO_X_RESIDENCY_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_ORTHO_X_RESIDENCY_ENABLED=ENABLED &
     &     .AND.RESIDENCY_ENABLED.AND.ORTHO_X_RESIDENCY_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_ORTHO_X_RESIDENCY_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_1COVERLAP_ENABLED(FLOPS)
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: FLOPS
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_1COVERLAP_ENABLED=ENABLED &
     &     .AND.ONECENTER_OVERLAP_ENABLED.AND.(FLOPS.GE.MINFLOP_OVERLAP)
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_1COVERLAP_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_DENMAT_ENERGY_ENABLED &
     &                                  (FLOPS)
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: FLOPS
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_DENMAT_ENERGY_ENABLED=ENABLED &
     &     .AND.DENMAT_ENERGY_ENABLED.AND.(FLOPS.GE.MINFLOP_DENMAT)
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_DENMAT_ENERGY_ENABLED
!
!     ..........................................................................
      LOGICAL(4) FUNCTION CPPAW_CUBLAS_ACC_INVERSION_BATCH_ENABLED()
      IMPLICIT NONE
!     **************************************************************************
      CALL CPPAW_CUBLAS_ACC_INITCONFIG
      CPPAW_CUBLAS_ACC_INVERSION_BATCH_ENABLED=ENABLED &
     &                                      .AND.INVERSION_BATCH_ENABLED
      RETURN
      END FUNCTION CPPAW_CUBLAS_ACC_INVERSION_BATCH_ENABLED
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
      SUBROUTINE CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_FORCE(ACTIVE)
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: ACTIVE
!     **************************************************************************
      WAVE_OVERLAP_FORCE_ACTIVE=ACTIVE
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_FORCE
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
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE_ADDPRODUCT(FLOPS)
      IF(.NOT.USED) RETURN
      IF(RESIDENCY_ENABLED) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D &
     &      ('ACC_PRESENT_CUBLAS_ZGEMM_NN_A' &
     &      ,'ACC_COPY_CUBLAS_ZGEMM_NN_A_IN',N,M,A)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D &
     &      ('ACC_PRESENT_CUBLAS_ZGEMM_NN_B' &
     &      ,'ACC_COPY_CUBLAS_ZGEMM_NN_B_IN',M,L,B)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D_IO &
     &      ('ACC_PRESENT_CUBLAS_ZGEMM_NN_C' &
     &      ,'ACC_COPY_CUBLAS_ZGEMM_NN_C_IO',N,L,C)
#ENDIF
!$ACC DATA PRESENT_OR_COPYIN(A(1:N,1:M),B(1:M,1:L)) &
!$ACC& PRESENT_OR_COPY(C(1:N,1:L))
        CALL CPPAW_CUBLAS_ACC_ZGEMM_NN_PRESENT(N,M,L,A,B,C)
!$ACC END DATA
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
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE_MATMUL(FLOPS)
      IF(.NOT.USED) RETURN
      IF(RESIDENCY_ENABLED) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D &
     &      ('ACC_PRESENT_ZGEMM_MAT_A','ACC_COPY_ZGEMM_MAT_A_IN' &
     &      ,N,M,A)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D &
     &      ('ACC_PRESENT_ZGEMM_MAT_B','ACC_COPY_ZGEMM_MAT_B_IN' &
     &      ,M,L,B)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D &
     &      ('ACC_PRESENT_ZGEMM_MAT_C','ACC_COPY_ZGEMM_MAT_C_OUT' &
     &      ,N,L,C)
#ENDIF
!$ACC DATA PRESENT_OR_COPYIN(A(1:N,1:M),B(1:M,1:L)) &
!$ACC& PRESENT_OR_COPYOUT(C(1:N,1:L))
        CALL CPPAW_CUBLAS_ACC_ZGEMM_MATMUL_PRESENT(N,M,L,A,B,C)
!$ACC END DATA
        RETURN
      END IF
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
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE_MATMUL(FLOPS)
      IF(.NOT.USED) RETURN
      IF(RESIDENCY_ENABLED) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_R8_2D &
     &      ('ACC_PRESENT_DGEMM_MAT_A','ACC_COPY_DGEMM_MAT_A_IN' &
     &      ,N,M,A)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_R8_2D &
     &      ('ACC_PRESENT_DGEMM_MAT_B','ACC_COPY_DGEMM_MAT_B_IN' &
     &      ,M,L,B)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_R8_2D &
     &      ('ACC_PRESENT_DGEMM_MAT_C','ACC_COPY_DGEMM_MAT_C_OUT' &
     &      ,N,L,C)
#ENDIF
!$ACC DATA PRESENT_OR_COPYIN(A(1:N,1:M),B(1:M,1:L)) &
!$ACC& PRESENT_OR_COPYOUT(C(1:N,1:L))
        CALL CPPAW_CUBLAS_ACC_DGEMM_MATMUL_PRESENT(N,M,L,A,B,C)
!$ACC END DATA
        RETURN
      END IF
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
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE_OVERLAP(FLOPS)
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
      SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NT_COPY(LEN1,LEN2,N,PSI1 &
     &                                         ,PSI2,OPERATOR,USED)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: LEN1
      INTEGER(4),INTENT(IN)  :: LEN2
      INTEGER(4),INTENT(IN)  :: N
      COMPLEX(8),INTENT(IN)  :: PSI1(LEN1,N)
      COMPLEX(8),INTENT(IN)  :: PSI2(LEN2,N)
      COMPLEX(8),INTENT(OUT) :: OPERATOR(LEN1,LEN2)
      LOGICAL(4),INTENT(OUT) :: USED
      REAL(8)                :: FLOPS
!     **************************************************************************
      FLOPS=8.D0*REAL(LEN1,KIND=8)*REAL(LEN2,KIND=8)*REAL(N,KIND=8)
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE_OFFDEN(FLOPS)
      IF(.NOT.USED) RETURN
!$ACC DATA COPYIN(PSI1(1:LEN1,1:N),PSI2(1:LEN2,1:N)) &
!$ACC& COPYOUT(OPERATOR(1:LEN1,1:LEN2))
      CALL CPPAW_CUBLAS_ACC_ZGEMM_NT_PRESENT(LEN1,LEN2,N,PSI1 &
     &                                      ,PSI2,OPERATOR)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_ZGEMM_NT' &
     &     ,LEN1,LEN2,N,0,16.D0*(REAL(LEN1,KIND=8)*REAL(N,KIND=8) &
     &     +REAL(LEN2,KIND=8)*REAL(N,KIND=8) &
     &     +REAL(LEN1,KIND=8)*REAL(LEN2,KIND=8)))
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NT_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NT_PRESENT(LEN1,LEN2,N,PSI1 &
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
      ISTAT=CUBLASZGEMM(HANDLE,CUBLAS_OP_N,CUBLAS_OP_T,LEN1,LEN2,N &
     &                 ,ONE,PSI1,LEN1,PSI2,LEN2,ZERO,OPERATOR,LEN1)
!$ACC END HOST_DATA
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUBLASZGEMM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUBLAS_ACC_ZGEMM_NT_PRESENT')
      END IF
      CALL CPPAW_CUBLAS_ACC_FINISH(ISTAT)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_ZGEMM_NT_PRESENT
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
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE_OVERLAP(FLOPS)
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
      CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D &
     &    ('ACC_PRESENT_PROJ_PRO','ACC_COPY_PROJ_PRO_IN' &
     &    ,NGL,LMNX,PRO)
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
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CHARACTER(32)          :: ACC_PRESENT_PSI1
      CHARACTER(32)          :: ACC_COPY_PSI1
      CHARACTER(32)          :: ACC_PRESENT_PSI2
      CHARACTER(32)          :: ACC_COPY_PSI2
      CHARACTER(32)          :: ACC_PRESENT_OUT
      CHARACTER(32)          :: ACC_COPY_OUT
      CHARACTER(32)          :: ACC_COPY_AGG
#ENDIF
!     **************************************************************************
      IF(TID) THEN
        FLOPS=4.D0*REAL(N1,KIND=8)*REAL(N1,KIND=8)*REAL(LEN,KIND=8)
      ELSE
        FLOPS=8.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)*REAL(LEN,KIND=8)
      END IF
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE_OVERLAP(FLOPS)
      IF(.NOT.USED) RETURN
      IF(RESIDENCY_ENABLED) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL CPPAW_CUBLAS_ACC_ZSPROD_NAMES('P1',ACC_PRESENT_PSI1 &
     &      ,ACC_COPY_PSI1)
        CALL CPPAW_CUBLAS_ACC_ZSPROD_NAMES('P2',ACC_PRESENT_PSI2 &
     &      ,ACC_COPY_PSI2)
        CALL CPPAW_CUBLAS_ACC_ZSPROD_NAMES('OUT',ACC_PRESENT_OUT &
     &      ,ACC_COPY_OUT)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D &
     &      (ACC_PRESENT_PSI1,ACC_COPY_PSI1,LEN,N1,PSI1)
        IF(.NOT.TID) THEN
          CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D &
     &        (ACC_PRESENT_PSI2,ACC_COPY_PSI2,LEN,N2,PSI2)
        END IF
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_2D &
     &      (ACC_PRESENT_OUT,ACC_COPY_OUT,N1,N2,OVERLAP)
#ENDIF
        IF(TID) THEN
!$ACC DATA PRESENT_OR_COPYIN(PSI1(1:LEN,1:N1)) &
!$ACC& PRESENT_OR_COPYOUT(OVERLAP(1:N1,1:N2))
          CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT(TID,LEN,N1,PSI1 &
     &                                               ,N2,PSI2,OVERLAP)
!$ACC END DATA
        ELSE
!$ACC DATA PRESENT_OR_COPYIN(PSI1(1:LEN,1:N1),PSI2(1:LEN,1:N2)) &
!$ACC& PRESENT_OR_COPYOUT(OVERLAP(1:N1,1:N2))
          CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT(TID,LEN,N1,PSI1 &
     &                                               ,N2,PSI2,OVERLAP)
!$ACC END DATA
        END IF
        RETURN
      END IF
      IF(TID) THEN
!$ACC DATA COPYIN(PSI1(1:LEN,1:N1)) COPYOUT(OVERLAP(1:N1,1:N2))
        CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_PRESENT(TID,LEN,N1,PSI1 &
     &                                             ,N2,PSI2,OVERLAP)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL CPPAW_CUBLAS_ACC_ZSPROD_COPY_NAME('GEN',ACC_COPY_AGG)
        CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES(ACC_COPY_AGG &
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
      CALL CPPAW_CUBLAS_ACC_ZSPROD_COPY_NAME('GEN',ACC_COPY_AGG)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES(ACC_COPY_AGG &
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
      CHARACTER(32)          :: ACC_COPY_OVL_RES
#ENDIF
!     **************************************************************************
      IF(TID) THEN
        FLOPS=4.D0*REAL(N1,KIND=8)*REAL(N1,KIND=8)*REAL(LEN,KIND=8)
      ELSE
        FLOPS=8.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8)*REAL(LEN,KIND=8)
      END IF
      USED=CPPAW_CUBLAS_ACC_WAVE_OVERLAP_RESIDENT_ACTIVE() &
     &     .AND.(WAVE_OVERLAP_FORCE_ACTIVE &
     &           .OR.CPPAW_CUBLAS_ACC_SHOULD_USE_OVERLAP(FLOPS))
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
      CALL CPPAW_CUBLAS_ACC_ZSPROD_COPY_NAME('OVL',ACC_COPY_OVL_RES)
      CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES(ACC_COPY_OVL_RES &
     &     ,LEN,N1,N2,0,16.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8))
#ENDIF
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_SCALARPRODUCT_RESIDENT_COPY
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_GAMMA_CORRECTION_RESIDENT(NGL,NDIM &
     &                                                       ,N1,PSI1 &
     &                                                       ,N2,PSI2 &
     &                                                       ,NGAMMA &
     &                                                       ,OVERLAP)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NGL
      INTEGER(4),INTENT(IN)  :: NDIM
      INTEGER(4),INTENT(IN)  :: N1
      INTEGER(4),INTENT(IN)  :: N2
      INTEGER(4),INTENT(IN)  :: NGAMMA
      COMPLEX(8),INTENT(IN)  :: PSI1(NGL,NDIM,N1)
      COMPLEX(8),INTENT(IN)  :: PSI2(NGL,NDIM,N2)
      COMPLEX(8),INTENT(INOUT):: OVERLAP(N1,N2)
      COMPLEX(8),ALLOCATABLE :: GAMMA1(:,:)
      COMPLEX(8),ALLOCATABLE :: GAMMA2(:,:)
      INTEGER(4)             :: I1
      INTEGER(4)             :: I2
      INTEGER(4)             :: IDIM
!     **************************************************************************
      IF(NGAMMA.EQ.0) RETURN
      ALLOCATE(GAMMA1(NDIM,N1))
      ALLOCATE(GAMMA2(NDIM,N2))
!$ACC DATA PRESENT(PSI1(1:NGL,1:NDIM,1:N1),PSI2(1:NGL,1:NDIM,1:N2)) &
!$ACC& COPYOUT(GAMMA1(1:NDIM,1:N1),GAMMA2(1:NDIM,1:N2))
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(PSI1,GAMMA1)
      DO I1=1,N1
        DO IDIM=1,NDIM
          GAMMA1(IDIM,I1)=PSI1(NGAMMA,IDIM,I1)
        ENDDO
      ENDDO
!$ACC END PARALLEL LOOP
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(PSI2,GAMMA2)
      DO I2=1,N2
        DO IDIM=1,NDIM
          GAMMA2(IDIM,I2)=PSI2(NGAMMA,IDIM,I2)
        ENDDO
      ENDDO
!$ACC END PARALLEL LOOP
!$ACC END DATA
      DO I1=1,N1
        DO I2=1,N2
          DO IDIM=1,NDIM
            OVERLAP(I1,I2)=OVERLAP(I1,I2) &
     &                    -CONJG(GAMMA1(IDIM,I1))*GAMMA2(IDIM,I2)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(GAMMA1)
      DEALLOCATE(GAMMA2)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_GAMMA_CORRECTION_RESIDENT
!
!     ..........................................................................
      SUBROUTINE CPPAW_CUBLAS_ACC_INVERSION_RESIDENT_COPY(NGL,NDIM,N1 &
     &                                                    ,PSI1,N2,PSI2 &
     &                                                    ,MINUSG,OVERLAP &
     &                                                    ,USED)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NGL
      INTEGER(4),INTENT(IN)  :: NDIM
      INTEGER(4),INTENT(IN)  :: N1
      INTEGER(4),INTENT(IN)  :: N2
      COMPLEX(8),INTENT(IN)  :: PSI1(NGL,NDIM,N1)
      COMPLEX(8),INTENT(IN)  :: PSI2(NGL,NDIM,N2)
      INTEGER(4),INTENT(IN)  :: MINUSG(NGL)
      COMPLEX(8),INTENT(OUT) :: OVERLAP(N1,N2)
      LOGICAL(4),INTENT(OUT) :: USED
      COMPLEX(8),ALLOCATABLE :: PSI2M(:,:)
      COMPLEX(8)             :: ONE
      COMPLEX(8)             :: ZERO
      REAL(8)                :: FLOPS
      INTEGER(4)             :: ISTAT
      INTEGER(4)             :: I1
      INTEGER(4)             :: I2
      INTEGER(4)             :: IDIM
      INTEGER(4)             :: IG
      INTEGER(4)             :: I
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
      REAL(8)                :: ACCEL_BYTES
#ENDIF
!     **************************************************************************
      FLOPS=8.D0*REAL(NGL*NDIM,KIND=8)*REAL(N1,KIND=8) &
     &          *REAL(N2,KIND=8)
      USED=CPPAW_CUBLAS_ACC_WAVE_OVERLAP_RESIDENT_ACTIVE() &
     &     .AND.CPPAW_CUBLAS_ACC_INVERSION_BATCH_ENABLED() &
     &     .AND.(WAVE_OVERLAP_FORCE_ACTIVE &
     &           .OR.CPPAW_CUBLAS_ACC_SHOULD_USE_OVERLAP(FLOPS))
      IF(.NOT.USED) RETURN
      ONE=(1.D0,0.D0)
      ZERO=(0.D0,0.D0)
      ALLOCATE(PSI2M(NGL*NDIM,N2))
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC DATA PRESENT(PSI1(1:NGL,1:NDIM,1:N1),PSI2(1:NGL,1:NDIM,1:N2)) &
!$ACC& PRESENT_OR_COPYIN(MINUSG(1:NGL)) CREATE(PSI2M(1:NGL*NDIM,1:N2)) &
!$ACC& COPYOUT(OVERLAP(1:N1,1:N2))
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(I) PRESENT(PSI2,PSI2M,MINUSG)
      DO I2=1,N2
        DO IDIM=1,NDIM
          DO IG=1,NGL
            I=(IDIM-1)*NGL+IG
            PSI2M(I,I2)=CONJG(PSI2(MINUSG(IG),IDIM,I2))
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL LOOP
      CALL CPPAW_CUBLAS_ACC_ENSURE
!$ACC HOST_DATA USE_DEVICE(PSI1,PSI2M,OVERLAP)
      ISTAT=CUBLASZGEMM(HANDLE,CUBLAS_OP_C,CUBLAS_OP_N,N1,N2,NGL*NDIM &
     &                  ,ONE,PSI1,NGL*NDIM,PSI2M,NGL*NDIM,ZERO &
     &                  ,OVERLAP,N1)
!$ACC END HOST_DATA
      IF(ISTAT.NE.0) THEN
        CALL ERROR$MSG('CUBLASZGEMM FAILED')
        CALL ERROR$I4VAL('STATUS',ISTAT)
        CALL ERROR$STOP('CPPAW_CUBLAS_ACC_INVERSION_RESIDENT_COPY')
      END IF
      CALL CPPAW_CUBLAS_ACC_FINISH(ISTAT)
!$ACC END DATA
      DO I1=1,N1
        DO I2=1,N2
          OVERLAP(I1,I2)=CONJG(OVERLAP(I1,I2))
        ENDDO
      ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      ACCEL_BYTES=16.D0*(REAL(NGL*NDIM,KIND=8)*REAL(N2,KIND=8) &
     &     +REAL(N1,KIND=8)*REAL(N2,KIND=8)) &
     &     +4.D0*REAL(NGL,KIND=8)
      CALL ACCELPROFILE$ADD('CUBLAS_ZGEMM_OVL_RES_INV' &
     &     ,INT(NGL*NDIM,KIND=8),INT(N1,KIND=8),INT(N2,KIND=8),0_8 &
     &     ,FLOPS,ACCEL_BYTES,ACCEL_T1-ACCEL_T0)
#ENDIF
      DEALLOCATE(PSI2M)
      RETURN
      END SUBROUTINE CPPAW_CUBLAS_ACC_INVERSION_RESIDENT_COPY
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
      USED=CPPAW_CUBLAS_ACC_SHOULD_USE_OVERLAP(FLOPS)
      IF(.NOT.USED) RETURN
      IF(RESIDENCY_ENABLED) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_R8_2D &
     &      ('ACC_PRESENT_CUBLAS_DSPROD_PSI1' &
     &      ,'ACC_COPY_CUBLAS_DSPROD_PSI1_IN',LEN,N1,PSI1)
        IF(.NOT.TID) THEN
          CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_R8_2D &
     &        ('ACC_PRESENT_CUBLAS_DSPROD_PSI2' &
     &        ,'ACC_COPY_CUBLAS_DSPROD_PSI2_IN',LEN,N2,PSI2)
        END IF
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_R8_2D &
     &      ('ACC_PRESENT_CUBLAS_DSPROD_OUT' &
     &      ,'ACC_COPY_CUBLAS_DSPROD_OUT',N1,N2,OVERLAP)
#ENDIF
        IF(TID) THEN
!$ACC DATA PRESENT_OR_COPYIN(PSI1(1:LEN,1:N1)) &
!$ACC& PRESENT_OR_COPYOUT(OVERLAP(1:N1,1:N2))
          CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_PRESENT(TID,LEN,N1,PSI1 &
     &                                                  ,N2,PSI2,OVERLAP)
!$ACC END DATA
        ELSE
!$ACC DATA PRESENT_OR_COPYIN(PSI1(1:LEN,1:N1),PSI2(1:LEN,1:N2)) &
!$ACC& PRESENT_OR_COPYOUT(OVERLAP(1:N1,1:N2))
          CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_PRESENT(TID,LEN,N1,PSI1 &
     &                                                  ,N2,PSI2,OVERLAP)
!$ACC END DATA
        END IF
        RETURN
      END IF
      IF(TID) THEN
!$ACC DATA COPYIN(PSI1(1:LEN,1:N1)) COPYOUT(OVERLAP(1:N1,1:N2))
        CALL CPPAW_CUBLAS_ACC_SCALARPRODUCT_R8_PRESENT(TID,LEN,N1,PSI1 &
     &                                                ,N2,PSI2,OVERLAP)
!$ACC END DATA
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL CPPAW_CUBLAS_ACC_PROFILE_BYTES('ACC_COPY_CUBLAS_DSPROD' &
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
