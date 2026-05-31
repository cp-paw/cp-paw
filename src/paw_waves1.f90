!........1.........2.........3.........4.........5.........6.........7.........8
!*******************************************************************************
!*******************************************************************************
!**                                                                           **
!**  WAVES OBJECT                                                             **
!**  THIS OBJECT IS CONTINUED IN PAW_WAVES2.F90                               **
!**                                                                           **
!**  PURPOSE: OPERATIONS ON THE PLANE WAVE PART OF THE WAVE FUNCTIONS         **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**    WAVES$SETR8A(ID,LEN,VAL)                                               **
!**    WAVES$GETR8A(ID,LEN,VAL)                                               **
!**    WAVES$INITIALIZE                                                       **
!**    WAVES$ETOT                                                             **
!**    WAVES$PROPAGATE                                                        **
!**    WAVES$ORTHOGONALIZE                                                    **
!**    WAVES$SWITCH                                                           **
!**    WAVES$REPORT(NFIL)                                                     **
!**                                                                           **
!**  DEPENDECIES:                                                             **
!**    ERROR                                                                  **
!**    MPELIB                                                                 **
!**    PLANEWAVE                                                              **
!**    SETUP                                                                  **
!**    ...                                                                    **
!**                                                                           **
!**  DEFINITION OF SPIN DIMENSIONS:                                           **
!**                           | NSPIN | NDIM | NDIMD |                        **
!**        -------------------------------------------                        **
!**        NON-SPIN-POLARIZED |   1   |   1  |   1   |                        **
!**        SPIN POLARIZED     |   2   |   1  |   2   |                        **
!**        NON-COLLINEAR      |   1   |   2  |   4   |                        **
!**                                                                           **
!**     NSPIN IS THE NUMBER OF SLATER DETERMINANTS                            **
!**     NDIM IS THE NUMBER OF SPINOR DIMENSIONS OF THE WAVE FUNCTION          **
!**     NDIMD IS THE NUMBER OF DENSITY COMPONENTS. NDIMD=NSPIN*NDIM**2        **
!**     TWO REPRESENTATIONS ARE USED:                                         **
!**           NSPIN=2,NDIM=1:  UP/DOWN AND TOTAL/SPIN                         **
!**           NSPIN=1,NDIM=2:  UPUP/UPDOWN/DOWNUP/DOWNDOWN AND                **
!**                            TOTAL/MX/MY/MZ                                 **
!**                                                                           **
!*******************************************************************************
!*******************************************************************************
MODULE RSPACEOP_MODULE
TYPE RSPACEMAT_TYPE
  INTEGER(4)      :: IAT1        ! 1ST ATOM (LINKED TO THE LEFT INDEX OF MAT)
  INTEGER(4)      :: IAT2        ! 2ND ATOM (LINKED TO THE RIGHT INDEX OF MAT)
  INTEGER(4)      :: IT(3)       ! LATTICE TRANSLATIONS TO BE ADDED TO ATOM 2
  INTEGER(4)      :: N1          ! LEFT DIMENSION OF MAT
  INTEGER(4)      :: N2          ! RIGHT DIMENSION OF MAT
  INTEGER(4)      :: N3          ! THIRD INDEX OF MAT
  REAL(8),ALLOCATABLE :: MAT(:,:,:)  ! (N1,N2,N3)
END TYPE RSPACEMAT_TYPE
CONTAINS
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE RSPACEOP$COPY(IN,OUT)
!      *************************************************************************
!      ** COPIES ONE ITEM OF RSPACEOP
!      *************************************************************************
       IMPLICIT NONE
       TYPE(RSPACEMAT_TYPE),INTENT(IN)    :: IN
       TYPE(RSPACEMAT_TYPE),INTENT(INOUT) :: OUT
!      *************************************************************************
       OUT%IAT1=IN%IAT1
       OUT%IAT2=IN%IAT2
       OUT%IT  =IN%IT
       IF(OUT%N1.NE.IN%N1.OR.OUT%N2.NE.IN%N2.OR.OUT%N3.NE.IN%N3) THEN
         IF(ALLOCATED(OUT%MAT))DEALLOCATE(OUT%MAT)
       END IF
       IF(.NOT.ALLOCATED(OUT%MAT)) THEN
         OUT%N1  =IN%N1
         OUT%N2  =IN%N2
         OUT%N3  =IN%N3
         ALLOCATE(OUT%MAT(IN%N1,IN%N2,IN%N3))
       END IF
       OUT%MAT=IN%MAT
       RETURN
       END SUBROUTINE RSPACEOP$COPY
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE RSPACEOP$DELETE(IN)
!      *************************************************************************
!      ** COPIES ONE ITEM OF RSPACEOP
!      *************************************************************************
       IMPLICIT NONE
       TYPE(RSPACEMAT_TYPE),INTENT(INOUT) :: IN
!      *************************************************************************
       IF(ALLOCATED(IN%MAT))DEALLOCATE(IN%MAT)
       IN%N1=-1
       IN%N2=-1
       IN%N3=-1
       IN%IAT1=-1
       IN%IAT2=-1
       IN%IT  =(/0,0,0/)
       RETURN
       END SUBROUTINE RSPACEOP$DELETE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RSPACEOP$WRITEMAT(NFIL,TITLE,NN,MAT)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN) :: NFIL
      CHARACTER(*)        ,INTENT(IN) :: TITLE
      INTEGER(4)          ,INTENT(IN) :: NN
      TYPE(RSPACEMAT_TYPE),INTENT(IN) :: MAT(NN)
      INTEGER(4)                      :: IN,I1,I3
      CHARACTER(128)                  :: FOMT
!     **************************************************************************
      WRITE(NFIL,FMT='(82("="))')
      WRITE(NFIL,FMT='(82("="),T10,"  ",A,"   ")')TRIM(TITLE)
      WRITE(NFIL,FMT='(82("="))')
      DO IN=1,NN
        IF(MAT(IN)%N1*MAT(IN)%N2.EQ.0) CYCLE
        FOMT='(82("="),T10,"IN=",I8," IAT1=",I5," IAT2=",I5," IT=",3I3," ")'
        WRITE(NFIL,FMT=FOMT)IN,MAT(IN)%IAT1,MAT(IN)%IAT2,MAT(IN)%IT
        DO I3=1,MAT(IN)%N3
          WRITE(NFIL,FMT='(82("-"),T10," COMPONENT ",I1,"  ")')I3
          DO I1=1,MAT(IN)%N1
            WRITE(NFIL,FMT='(20F10.5)')MAT(IN)%MAT(I1,:,I3)
          ENDDO
        ENDDO
        WRITE(NFIL,FMT='(82("="))')
      ENDDO
      RETURN
      END SUBROUTINE RSPACEOP$WRITEMAT
END MODULE RSPACEOP_MODULE
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE WAVES_MODULE                                                
!*******************************************************************************
!**                                                                           **
!**                                                                           **
!************************************************P.E. BLOECHL, (1995)***********
USE LINKEDLIST_MODULE
USE RSPACEOP_MODULE, ONLY: RSPACEMAT_TYPE  ! (IN PAW_WAVES1.F90)
! ------------------------------------------------------------------------------
! THE TYPE MAP_TYPE DESCRIBES THE ARRANGEMENT OF PROJECTOR FUNCTIONS 
! AND THEIR LINKS TO ATOMS AN ANGULAR MOMENTA
TYPE MAP_TYPE
  INTEGER(4)         :: NAT         ! #(ATOMS)
  INTEGER(4)         :: NSP         ! #(ATOM TYPES)
  INTEGER(4)         :: NBAREPRO    ! #(LN*NAT)
  INTEGER(4)         :: NRL         ! #(LOCAL REAL SPACE GRID POINTS)
  INTEGER(4)         :: NPRO        ! #(LMN*NAT)
  INTEGER(4)         :: LMX
  INTEGER(4)         :: LNXX        ! X(LNX)
  INTEGER(4),POINTER :: LNX(:)      !(NSP)
  INTEGER(4),POINTER :: LMNX(:)     !(NSP)
  INTEGER(4),POINTER :: LOX(:,:)    !(LNXX,NSP) MAIN ANGULAR MOMENTUM OF AN LN
  INTEGER(4),POINTER :: ISP(:)      !(NAT) POINTER TO ATOM TYPES
END TYPE MAP_TYPE
! ------------------------------------------------------------------------------
! THE TYPE GSET_TYPE KEEPS PROJECTOR FUNCTIONS SPHERICAL HARMONICS AND
! STRUCTURE FACTORS FOR A GIVEN K-POINT
! GSET_TYPE REFERS TO MAP_TYPE FOR DIMENSIONING OF ITS ARRAYS
TYPE GSET_TYPE 
  CHARACTER(16)      :: ID         ! GSET IDENTIFIER
  LOGICAL(4)         :: TINV       ! SWITCH: TIME INVERSION SYMMETRY
  INTEGER(4)         :: NGL        ! #(PLANE WAVES(LOCAL))
  REAL(8)   ,POINTER :: PRO(:,:)   ! (NGL,NBAREPRO) PROJECTOR
  REAL(8)   ,POINTER :: DPRO(:,:)  ! (NGL,NBAREPRO) PROJECTOR |G|DP/DG-L*P
  REAL(8)   ,POINTER :: YLM(:,:)   ! (NGL,LMXX) REAL SPHERICAL HARMONICS
  REAL(8)   ,POINTER :: SYLM(:,:,:)! (NGL,LMXX,6) STRAINED SPHERICAL HARMONICS
  REAL(8)   ,POINTER :: MPSI(:)    ! (NGL) WAVE FUNCTION MASS
!PB070802  REAL(8)   ,POINTER :: DMPSI(:)   ! (NGL) DMPSI/DT
  REAL(8)   ,POINTER :: BUCKET(:)  ! BUCKET POTENTIAL
  REAL(8)   ,POINTER :: DBUCKET(:) ! 1/G* DBUCKET/DG
END TYPE GSET_TYPE
TYPE WVSET_TYPE  !==============================================================
  TYPE(GSET_TYPE),POINTER :: GSET
  INTEGER(4)         :: NB
  INTEGER(4)         :: NBH
  COMPLEX(8),POINTER :: PSI0(:,:,:)=>NULL()     !(NGL,NDIM,NBH)  PSPSI(0)
  COMPLEX(8),POINTER :: PSIM(:,:,:)=>NULL()     !(NGL,NDIM,NBH)  PSPSI(-,+)(G)
  COMPLEX(8),POINTER :: PROJ(:,:,:)=>NULL()     !(NDIM,NBH,NPRO) <PSPSI|P>
  COMPLEX(8),POINTER :: HPROJ(:,:,:)=>NULL()    !(NDIM,NBH,NPRO) 
  COMPLEX(8),POINTER :: TBC_NEW(:,:,:)=>NULL()  !(NDIM,NBH,NORB) |PSI>=|CHI>*TBC
  COMPLEX(8),POINTER :: HTBC_NEW(:,:,:)=>NULL()     !(NDIM,NBH,NPRO) DE/DTBC
  COMPLEX(8),POINTER :: HPSI(:,:,:)=>NULL()     !(NGWLX,NB,IDIM)
                                        ! +(WAVES$HPSI)-(WAVES$PROPAGATE)
  COMPLEX(8),POINTER :: OPSI(:,:,:)     !(NGWLX,NB,IDIM)
                                        ! +(WAVES$HPSI)-(WAVES$PROPAGATE)
  COMPLEX(8),POINTER :: RLAM0(:,:)      !(NB,NB)
  COMPLEX(8),POINTER :: RLAMM(:,:)      !(NB,NB)
  COMPLEX(8),POINTER :: RLAM2M(:,:)     !(NB,NB)
  COMPLEX(8),POINTER :: RLAM3M(:,:)     !(NB,NB)
  COMPLEX(8),POINTER :: EIGVEC(:,:)     !(NB,NB) ENERGY EIGENVALUES AFTER DIAG
  REAL(8)   ,POINTER :: EIGVAL(:)       !(NB)    EIGENVALUES AFTER DIAG
  REAL(8)   ,POINTER :: EXPECTVAL(:)    !(NB) !<PSI_N|H|PSI_N>
  LOGICAL(4)         :: HPSI_ACC_RESIDENT=.FALSE.
  LOGICAL(4)         :: PSIM_ACC_RESIDENT=.FALSE.
  LOGICAL(4)         :: PROJ_ACC_RESIDENT=.FALSE.
END TYPE WVSET_TYPE
TYPE EXTERNALPOINTER_TYPE !=====================================================
  INTEGER(4) :: IB
  INTEGER(4) :: IKPT
  INTEGER(4) :: ISPIN
  INTEGER(4) :: IAT
  LOGICAL(4) :: TIM
END TYPE EXTERNALPOINTER_TYPE
!===============================================================================
!== DATA THAT DESCRIBE FUNCTIONALITY OF WAVES OBJECT                          ==
!==  THESE DATA HAVE TO BE SET EXPLICITELY                                    ==
!===============================================================================
INTEGER(4)  :: NKPT          ! #(K-POINTS) SEE ALSO NKPTL
INTEGER(4)  :: NSPIN=1       ! #(SPINS)
INTEGER(4)  :: NDIM=1        ! #(WAVE FUNCTION COMPONENTS)
INTEGER(4)  :: NDIMD=1       ! #(WAVE FUNCTION COMPONENTS)
REAL(8)     :: EPWPSI=0.D0   ! PLANE WAVE CUTOFF GMAX**2/2
LOGICAL(4)  :: TBUCKET=.FALSE.
REAL(8)     :: EPWPSI0=0.D0  ! PARAMETER FOR PLANE WAVE CUTOFF POTENTIAL
REAL(8)     :: D2EPWPSI=0.D0 ! PARAMETER FOR PLANE WAVE CUTOFF POTENTIAL
REAL(8)     :: EPWRHO=0.D0   ! PLANE WAVE CUTOFF GMAX**2/2
REAL(8)     :: EMASS=0.D0    ! MPSI(G)=EMASS*(1+EMASSCG2*G2)
REAL(8)     :: EMASSCG2=0.D0 ! MPSI(G)=EMASS*(1+EMASSCG2*G2)
REAL(8)     :: ANNEE=0.D0    ! FRICTION 
LOGICAL(4)  :: TSTOP=.FALSE. ! INITIAL VELOCITY SET TO ZERO
CHARACTER(32):: RSTRTTYPE='DYNAMIC'
LOGICAL(4)  :: TRANDOM=.FALSE.    ! INITIAL VELOCITIES RANDOMIZED
REAL(8)     :: AMPRANDOM=0.D0
LOGICAL(4)  :: TSAFEORTHO=.TRUE.  ! CHOICE OR ORTHOGONALIZATION
LOGICAL(4)  :: TSWAPSTATES=.FALSE. ! CHOICE OR SWITCHING STATES IN A BAND CROSSING
LOGICAL(4)  :: TSTRAIGHTEN=.TRUE. ! TRANSFORM TO EIGENSTATES
REAL(8)     :: DELT=0.D0          ! TIME STEP IN A.U.
LOGICAL(4)  :: THAMILTON=.FALSE.  ! HAMILTON MATRIX AVAILABLE
LOGICAL(4)  :: TRAWSTATES=.FALSE. ! PROVIDES NON-DIAGONALIZED WAVE FUNCTIONS THROUGH $GETR
LOGICAL(4)  :: TFORCEX=.TRUE.
LOGICAL(4)  :: TSTRESSX=.FALSE.
REAL(8)     :: SCALERCUT=2.D0  ! NEIGHBORLIST FOR OFFSITE DENSITY MATRIX CONTAINS
                             ! PAIRS WITH DIS <(RCOV1+RCOV2)*SCALERCUT
!===============================================================================
!== PERMANENT DATA, WHICH ARE ORGANIZED BY THE ROUTINES ITSELF                ==
!===============================================================================
INTEGER(4)      ,SAVE     :: NKPTL         ! #(LOCAL K-POINTS)
INTEGER(4)      ,POINTER  :: KMAP(:)
REAL(8)                   :: WAVEEKIN1=0.D0   ! CALCULATED IN WAVES$ETOT
REAL(8)                   :: WAVEEKIN2=0.D0   ! CALULATED IN WAVES$ORTHOGONALIZE
TYPE(WVSET_TYPE),POINTER  :: THISARRAY(:,:)   ! (NKPTL,NSPIN)
TYPE(WVSET_TYPE),POINTER  :: THIS            ! CURRENT SET OF WAVES
TYPE(GSET_TYPE) ,POINTER  :: GSET            ! CURRENT SET OF GSET
TYPE(MAP_TYPE)            :: MAP
TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: OSDENMAT(:)  !OFFSITE DENSITY MATRIX      
TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: OSHAMIL(:)      
LOGICAL(4)                :: TPR=.FALSE.
LOGICAL(4)                :: TFIRST=.TRUE.
TYPE(EXTERNALPOINTER_TYPE):: EXTPNTR
LOGICAL(4)                :: TFIXRHO=.FALSE.
LOGICAL(4)                :: TWRITERHO=.FALSE.
CHARACTER(8)              :: OPTIMIZERTYPE  ! SWITCH FOR CONJUGATE GRADIENT OR DYNAMICS
#IF DEFINED(CPPVAR_CUBLAS_ACC)
LOGICAL(4)                :: ACC_PRO_CACHE_READY=.FALSE.
INTEGER(4)                :: ACC_PRO_CACHE_NGL=0
INTEGER(4)                :: ACC_PRO_CACHE_NAT=0
INTEGER(4)                :: ACC_PRO_CACHE_NPRO=0
INTEGER(4)                :: ACC_PRO_CACHE_LMNXX=0
CHARACTER(16)             :: ACC_PRO_CACHE_GSET_ID=' '
REAL(8)     ,ALLOCATABLE :: ACC_PRO_CACHE_R(:,:)
COMPLEX(8)  ,ALLOCATABLE :: ACC_PRO_CACHE_PRO(:,:,:)
COMPLEX(8)  ,ALLOCATABLE :: ACC_PRO_CACHE_EIGR(:,:)
#ENDIF
CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_SELECTWV(IKPT,ISPIN)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT
      INTEGER(4),INTENT(IN) :: ISPIN
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THISARRAY)) THEN
        CALL ERROR$MSG('THISARRAY DOES NOT EXIST') 
        CALL ERROR$STOP('WAVES_SELECTWV')
      END IF 
      IF(IKPT.GT.NKPTL) THEN
        CALL ERROR$MSG('IKPT OUT OF RANGE')
        CALL ERROR$I4VAL('IKPT',IKPT)
        CALL ERROR$I4VAL('NKPTL',NKPTL)
        CALL ERROR$STOP('WAVES_SELECTWV')
      END IF
      IF(ISPIN.GT.NSPIN) THEN
        CALL ERROR$MSG('ISPIN OUT OF RANGE')
        CALL ERROR$I4VAL('ISPIN',ISPIN)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('WAVES_SELECTWV')
      END IF
      THIS=>THISARRAY(IKPT,ISPIN)
      GSET=>THIS%GSET
      RETURN
      END SUBROUTINE WAVES_SELECTWV
!
#IF DEFINED(CPPVAR_CUBLAS_ACC)
!     ..........................................................................
      SUBROUTINE WAVES_ACC_THIS_PROJ_DELETE()
      USE OPENACC
      IMPLICIT NONE
!     **************************************************************************
      IF(.NOT.ASSOCIATED(THIS)) RETURN
      IF(.NOT.ASSOCIATED(THIS%PROJ)) RETURN
      IF(ACC_IS_PRESENT(THIS%PROJ)) THEN
!$ACC EXIT DATA DELETE(THIS%PROJ(1:NDIM,1:THIS%NBH,1:MAP%NPRO))
      END IF
      THIS%PROJ_ACC_RESIDENT=.FALSE.
      RETURN
      END SUBROUTINE WAVES_ACC_THIS_PROJ_DELETE
!
!     ..........................................................................
      SUBROUTINE WAVES_ACC_THIS_PROJ_COPYIN()
      USE OPENACC
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &      CPPAW_CUBLAS_ACC_PROJ_RESIDENCY_ENABLED
      IMPLICIT NONE
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                 :: ACCEL_T0
      REAL(8)                 :: ACCEL_T1
      REAL(8)                 :: ACCEL_BYTES
#ENDIF
!     **************************************************************************
      IF(.NOT.CPPAW_CUBLAS_ACC_PROJ_RESIDENCY_ENABLED()) RETURN
      IF(.NOT.ASSOCIATED(THIS)) RETURN
      IF(.NOT.ASSOCIATED(THIS%PROJ)) RETURN
      IF(THIS%NBH.LE.0.OR.MAP%NPRO.LE.0) RETURN
      IF(ACC_IS_PRESENT(THIS%PROJ)) THEN
        THIS%PROJ_ACC_RESIDENT=.TRUE.
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$ADD('ACC_PRESENT_THIS_PROJ' &
     &      ,INT(NDIM,KIND=8),INT(THIS%NBH,KIND=8) &
     &      ,INT(MAP%NPRO,KIND=8),0_8,0.D0,0.D0,0.D0)
#ENDIF
        RETURN
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      ACCEL_BYTES=16.D0*REAL(NDIM,KIND=8)*REAL(THIS%NBH,KIND=8) &
     &          *REAL(MAP%NPRO,KIND=8)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC ENTER DATA COPYIN(THIS%PROJ(1:NDIM,1:THIS%NBH,1:MAP%NPRO))
!$ACC WAIT
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('ACC_COPY_THIS_PROJ_IN' &
     &    ,INT(NDIM,KIND=8),INT(THIS%NBH,KIND=8) &
     &    ,INT(MAP%NPRO,KIND=8),0_8,0.D0,ACCEL_BYTES &
     &    ,ACCEL_T1-ACCEL_T0)
#ENDIF
      THIS%PROJ_ACC_RESIDENT=.TRUE.
      RETURN
      END SUBROUTINE WAVES_ACC_THIS_PROJ_COPYIN
!
!     ..........................................................................
      SUBROUTINE WAVES_ACC_PRO_CACHE_INVALIDATE()
      IMPLICIT NONE
!     **************************************************************************
      IF(ACC_PRO_CACHE_READY) THEN
!$ACC EXIT DATA DELETE(ACC_PRO_CACHE_PRO,ACC_PRO_CACHE_EIGR)
      END IF
      IF(ALLOCATED(ACC_PRO_CACHE_PRO)) DEALLOCATE(ACC_PRO_CACHE_PRO)
      IF(ALLOCATED(ACC_PRO_CACHE_EIGR)) DEALLOCATE(ACC_PRO_CACHE_EIGR)
      IF(ALLOCATED(ACC_PRO_CACHE_R)) DEALLOCATE(ACC_PRO_CACHE_R)
      ACC_PRO_CACHE_READY=.FALSE.
      ACC_PRO_CACHE_NGL=0
      ACC_PRO_CACHE_NAT=0
      ACC_PRO_CACHE_NPRO=0
      ACC_PRO_CACHE_LMNXX=0
      ACC_PRO_CACHE_GSET_ID=' '
      RETURN
      END SUBROUTINE WAVES_ACC_PRO_CACHE_INVALIDATE
!
!     ..........................................................................
      LOGICAL(4) FUNCTION WAVES_ACC_PRO_CACHE_IS_VALID(GSETID,NGL,NAT &
     &                                                ,NPRO,LMNXX,R)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: GSETID
      INTEGER(4)  ,INTENT(IN) :: NGL
      INTEGER(4)  ,INTENT(IN) :: NAT
      INTEGER(4)  ,INTENT(IN) :: NPRO
      INTEGER(4)  ,INTENT(IN) :: LMNXX
      REAL(8)     ,INTENT(IN) :: R(3,NAT)
      REAL(8)     ,PARAMETER  :: RTOL=1.D-12
!     **************************************************************************
      WAVES_ACC_PRO_CACHE_IS_VALID=.FALSE.
      IF(.NOT.ACC_PRO_CACHE_READY) RETURN
      IF(ACC_PRO_CACHE_NGL.NE.NGL) RETURN
      IF(ACC_PRO_CACHE_NAT.NE.NAT) RETURN
      IF(ACC_PRO_CACHE_NPRO.NE.NPRO) RETURN
      IF(ACC_PRO_CACHE_LMNXX.NE.LMNXX) RETURN
      IF(TRIM(ACC_PRO_CACHE_GSET_ID).NE.TRIM(GSETID)) RETURN
      IF(.NOT.ALLOCATED(ACC_PRO_CACHE_R)) RETURN
      IF(MAXVAL(ABS(ACC_PRO_CACHE_R(1:3,1:NAT)-R(1:3,1:NAT))).GT.RTOL) RETURN
      WAVES_ACC_PRO_CACHE_IS_VALID=.TRUE.
      RETURN
      END FUNCTION WAVES_ACC_PRO_CACHE_IS_VALID
!
!     ..........................................................................
      SUBROUTINE WAVES_ACC_PRO_CACHE_BUILD(MAP,GSET,NAT,R,NGL,LMNXX)
      IMPLICIT NONE
      TYPE(MAP_TYPE) ,INTENT(IN) :: MAP
      TYPE(GSET_TYPE),INTENT(IN) :: GSET
      INTEGER(4)     ,INTENT(IN) :: NAT
      REAL(8)        ,INTENT(IN) :: R(3,NAT)
      INTEGER(4)     ,INTENT(IN) :: NGL
      INTEGER(4)     ,INTENT(IN) :: LMNXX
      INTEGER(4)                 :: IAT,ISP,LNX,LMNX,IBPRO,LNXX
      INTEGER(4)     ,ALLOCATABLE:: LOX(:)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                    :: ACCEL_T0
      REAL(8)                    :: ACCEL_T1
      REAL(8)                    :: ACCEL_BYTES
#ENDIF
!     **************************************************************************
      IF(WAVES_ACC_PRO_CACHE_IS_VALID(GSET%ID,NGL,NAT,MAP%NPRO,LMNXX,R)) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$ADD('ACC_PRESENT_PRO_CACHE_REUSE' &
     &      ,INT(NGL,KIND=8),INT(LMNXX,KIND=8),INT(NAT,KIND=8),0_8 &
     &      ,0.D0,0.D0,0.D0)
#ENDIF
        RETURN
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      CALL WAVES_ACC_PRO_CACHE_INVALIDATE()
      ACC_PRO_CACHE_NGL=NGL
      ACC_PRO_CACHE_NAT=NAT
      ACC_PRO_CACHE_NPRO=MAP%NPRO
      ACC_PRO_CACHE_LMNXX=LMNXX
      ACC_PRO_CACHE_GSET_ID=GSET%ID
      ALLOCATE(ACC_PRO_CACHE_R(3,NAT))
      ALLOCATE(ACC_PRO_CACHE_PRO(NGL,LMNXX,NAT))
      ALLOCATE(ACC_PRO_CACHE_EIGR(NGL,NAT))
      ACC_PRO_CACHE_R(1:3,1:NAT)=R(1:3,1:NAT)
!$ACC ENTER DATA CREATE(ACC_PRO_CACHE_PRO(1:NGL,1:LMNXX,1:NAT) &
!$ACC&                  ,ACC_PRO_CACHE_EIGR(1:NGL,1:NAT))
      LNXX=MAXVAL(MAP%LNX)
      ALLOCATE(LOX(LNXX))
!$ACC DATA PRESENT_OR_COPYIN(GSET%PRO(1:NGL,1:MAP%NBAREPRO) &
!$ACC&       ,GSET%YLM(1:NGL,1:MAP%LMX)) &
!$ACC& PRESENT(ACC_PRO_CACHE_PRO(1:NGL,1:LMNXX,1:NAT) &
!$ACC&        ,ACC_PRO_CACHE_EIGR(1:NGL,1:NAT))
      DO IAT=1,NAT
        ISP=MAP%ISP(IAT)
        LNX=MAP%LNX(ISP)
        LMNX=MAP%LMNX(ISP)
        LOX=MAP%LOX(:,ISP)
        IBPRO=1+SUM(MAP%LNX(1:ISP-1))
        CALL PLANEWAVE$STRUCTUREFACTOR(R(1,IAT),NGL,ACC_PRO_CACHE_EIGR(1,IAT))
!$ACC UPDATE DEVICE(ACC_PRO_CACHE_EIGR(1:NGL,IAT))
        CALL WAVES_EXPANDPRO_ACC(LNX,LOX,LMNX,NGL &
     &      ,GSET%PRO(:,IBPRO:IBPRO+LNX-1),MAP%LMX,GSET%YLM &
     &      ,ACC_PRO_CACHE_EIGR(:,IAT),ACC_PRO_CACHE_PRO(:,:,IAT))
      ENDDO
!$ACC END DATA
      DEALLOCATE(LOX)
      ACC_PRO_CACHE_READY=.TRUE.
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      ACCEL_BYTES=8.D0*REAL(NGL,KIND=8)*REAL(MAP%NBAREPRO,KIND=8) &
     &          +8.D0*REAL(NGL,KIND=8)*REAL(MAP%LMX,KIND=8) &
     &          +16.D0*REAL(NGL,KIND=8)*REAL(NAT,KIND=8) &
     &          +16.D0*REAL(NGL,KIND=8)*REAL(LMNXX,KIND=8) &
     &          *REAL(NAT,KIND=8)
      CALL ACCELPROFILE$ADD('ACC_BUILD_PRO_CACHE' &
     &      ,INT(NGL,KIND=8),INT(LMNXX,KIND=8),INT(NAT,KIND=8),0_8 &
     &      ,0.D0,ACCEL_BYTES,ACCEL_T1-ACCEL_T0)
#ENDIF
      RETURN
      END SUBROUTINE WAVES_ACC_PRO_CACHE_BUILD
#ENDIF
END MODULE WAVES_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$STATESELECTED(TCHK)
!     **************************************************************************
!     **  WAVES$STATESELECTED                                                 **
!     **  TESTS IF THE EXTERNAL POINTER SELECTS A STATE THAT IS               **
!     **  AVAILABLE ON THIS TASK                                              **
!     **************************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(OUT):: TCHK
!     **************************************************************************
      TCHK=.TRUE.
      TCHK=TCHK.AND.EXTPNTR%IKPT.NE.0
      TCHK=TCHK.AND.EXTPNTR%ISPIN.NE.0
      TCHK=TCHK.AND.EXTPNTR%IB.NE.0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$SELECTSTATEPOINTER(IB,IKPT,ISPIN,TCHK)
!     **************************************************************************
!     **  SELECTS A SPECIFIC STATE OR WAVES$GET.. ROUTINES.                   **
!     **  IF THE PRESENT TASK DOES NOT HAVE INFORMATION ON THIS STATE         **
!     **  TCHK WILL BE FALSE. OTHERWISE IT WILL BE TRUE                       **
!     **************************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IB
      INTEGER(4),INTENT(IN) :: IKPT
      INTEGER(4),INTENT(IN) :: ISPIN
      LOGICAL(4),INTENT(OUT):: TCHK
!     **************************************************************************
      CALL WAVES$SETI4('IKPT',IKPT)
      CALL WAVES$SETI4('ISPIN',ISPIN)
      CALL WAVES$SETI4('IB',IB)
      TCHK=.TRUE.
      TCHK=TCHK.AND.EXTPNTR%IKPT.NE.0
      TCHK=TCHK.AND.EXTPNTR%ISPIN.NE.0
      TCHK=TCHK.AND.EXTPNTR%IB.NE.0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$SETR8(ID,VAL)
!     **************************************************************************
!     **  WAVES$SETR8A                                                        **
!     **************************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'EPWRHO') THEN
        EPWRHO=VAL
      ELSE IF(ID.EQ.'EPWPSI') THEN
        EPWPSI=VAL
      ELSE IF(ID.EQ.'EMASS') THEN
        EMASS=VAL
      ELSE IF(ID.EQ.'EMASSCG2') THEN
        EMASSCG2=VAL
      ELSE IF(ID.EQ.'EPWBUCKET') THEN
        EPWPSI0=VAL
      ELSE IF(ID.EQ.'BUCKETPAR') THEN
        D2EPWPSI=VAL
        TBUCKET=(D2EPWPSI.NE.0.D0)
      ELSE IF(ID.EQ.'FRICTION') THEN
        ANNEE=VAL
      ELSE IF(ID.EQ.'TIMESTEP') THEN
        DELT=VAL
      ELSE IF(ID.EQ.'AMPRE') THEN
        AMPRANDOM=VAL
      ELSE IF(ID.EQ.'SCALERCUT') THEN
        SCALERCUT=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$SETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$GETR8(ID,VAL)
!     **************************************************************************
!     **  WAVES$GETR8                                                         **
!     **************************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ID.EQ.'EPWRHO') THEN
        VAL=EPWRHO
      ELSE IF(ID.EQ.'EPWPSI') THEN
        VAL=EPWPSI
      ELSE IF(ID.EQ.'EMASS') THEN
        VAL=EMASS
      ELSE IF(ID.EQ.'EMASSCG2') THEN
        VAL=EMASSCG2
      ELSE IF(ID.EQ.'FRICTION') THEN
        VAL=ANNEE
      ELSE IF(ID.EQ.'EKIN(PSI)') THEN
        IF(WAVEEKIN2.EQ.0.D0) THEN
          VAL=WAVEEKIN1
        ELSE 
          VAL=0.5D0*(WAVEEKIN1+WAVEEKIN2)
        END IF
      ELSE IF(ID.EQ.'SCALERCUT') THEN
        VAL=SCALERCUT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$GETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$GETR8A(ID,LEN,VAL)
!     **************************************************************************
!     **  WAVES$GETR8A                                                        **
!     **************************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
      COMPLEX(8)  ,ALLOCATABLE:: CWORK1(:)
      COMPLEX(8)  ,ALLOCATABLE:: CWORK2(:)
      COMPLEX(8)  ,ALLOCATABLE:: CWORK3(:)
      INTEGER(4)              :: IKPT,ISPIN,IB,IDIM,IBH,IB1,IB2,IG,IR
      INTEGER(4)              :: IPRO,IPRO1,IPRO2,IAT,ISP,L,LN,LMNX,LMN
      INTEGER(4)              :: NB,NBH,NRL,NAT,NGL
      COMPLEX(8)              :: CSVAR1,CSVAR2
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)              :: NTASKS,THISTASK
      LOGICAL(4)              :: TINV
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$GETR8A')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(ID.EQ.'XXXXXXXXXXXXXX') THEN
!     
!     ==========================================================================
!     ==  REAL SPACE PS WAVE FUNCTION                               ==
!     ==  NOTE: IB,IKPT,ISPIN MUST BE SET ON LINKEDLIST             ==
!     ==========================================================================
!     ELSE IF(ID.EQ.'<PSI|H|PSI>') THEN
!
!     ==========================================================================
!     ==  ENERGY EIGENVALUES                                        ==
!     ==  NOTE: IKPT,ISPIN MUST BE SET ON LINKEDLIST                ==
!     ==========================================================================
      ELSE IF(ID.EQ.'EIGVAL') THEN
        IKPT=EXTPNTR%IKPT
        ISPIN=EXTPNTR%ISPIN
        IF(IKPT.EQ.0) THEN
          CALL ERROR$MSG('STATE NOT AVAILABLE ON THIS TASK')
          CALL ERROR$MSG('THIS ERROR OCCURS ONLY FOR PARALLEL JOBS')
          CALL ERROR$MSG("USE WAVES$GETL4('AVAILABLESTATE',TCHK) TO EXPLORE")
          CALL ERROR$I4VAL('THISTASK(MOMOMER)',THISTASK)
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('ISPIN',ISPIN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        CALL WAVES_SELECTWV(IKPT,ISPIN)
        CALL PLANEWAVE$SELECT(GSET%ID)
        IF(.NOT.ASSOCIATED(THIS%EIGVAL)) THEN
          CALL ERROR$MSG('EIGENVALUES NOT AVAILABLE')
          CALL ERROR$L4VAL('THAMILTON',THAMILTON)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        IF(LEN.NE.THIS%NB) THEN
          CALL ERROR$MSG('INCONSISTENT LENGTH')
          CALL ERROR$I4VAL('NB',THIS%NB)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        VAL(:)=THIS%EIGVAL(:)
!
!     ==========================================================================
!     ==  ENERGY EXPECTATIONS VALUES OF THE ONE-PARTCLE STATES      ==
!     ==  NOTE: IKPT,ISPIN MUST BE SET ON LINKEDLIST                ==
!     ==========================================================================
      ELSE IF(ID.EQ.'<PSI|H|PSI>') THEN
        IKPT=EXTPNTR%IKPT
        IF(IKPT.EQ.0) THEN
          CALL ERROR$MSG('STATE NOT AVAILABLE ON THIS TASK')
          CALL ERROR$MSG('THIS ERROR OCCURS ONLY FOR PARALLEL JOBS')
          CALL ERROR$MSG("USE WAVES$GETL4('AVAILABLESTATE',TCHK) TO EXPLORE")
          CALL ERROR$I4VAL('THISTASK(MOMOMER)',THISTASK)
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('ISPIN',ISPIN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        ISPIN=EXTPNTR%ISPIN
        CALL WAVES_SELECTWV(IKPT,ISPIN)
        CALL PLANEWAVE$SELECT(GSET%ID)
        IF(.NOT.ASSOCIATED(THIS%EXPECTVAL)) THEN
          CALL ERROR$MSG('ENERGY EXPECTATION VALUES NOT AVAILABLE')
          CALL ERROR$L4VAL('THAMILTON',THAMILTON)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        IF(LEN.NE.THIS%NB) THEN
          CALL ERROR$MSG('INCONSISTENT LENGTH')
          CALL ERROR$I4VAL('NB',THIS%NB)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        VAL(:)=THIS%EXPECTVAL(:)
!
!     ==========================================================================
!     ==  REAL SPACE PS WAVE FUNCTION                               ==
!     ==  NOTE: IB,IKPT,ISPIN MUST BE SET ON LINKEDLIST             ==
!     ==========================================================================
      ELSE IF(ID.EQ.'PSPSI') THEN
        IKPT=EXTPNTR%IKPT
        IF(IKPT.EQ.0) THEN
          CALL ERROR$MSG('STATE NOT AVAILABLE ON THIS TASK')
          CALL ERROR$MSG('THIS ERROR OCCURS ONLY FOR PARALLEL JOBS')
          CALL ERROR$MSG("USE WAVES$GETL4('AVAILABLESTATE',TCHK) TO EXPLORE")
          CALL ERROR$I4VAL('THISTASK(MOMOMER)',THISTASK)
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('ISPIN',ISPIN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        ISPIN=EXTPNTR%ISPIN
        IB=EXTPNTR%IB
        IF(NSPIN.EQ.1.AND.NDIM.EQ.2) THEN
          IDIM=ISPIN
          ISPIN=1
        ELSE
          IDIM=1
        END IF
        CALL WAVES_SELECTWV(IKPT,ISPIN)
        CALL PLANEWAVE$SELECT(GSET%ID)
        IF(.NOT.ASSOCIATED(THIS%EIGVEC)) THEN
          CALL ERROR$MSG('EIGENSTATES NOT AVAILABLE')
          CALL ERROR$L4VAL('THAMILTON',THAMILTON)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        NGL=GSET%NGL
        NB=THIS%NB
        NBH=THIS%NBH
        TINV=GSET%TINV
        ALLOCATE(CWORK1(NGL))
        CWORK1(:)=(0.D0,0.D0)
        IF(TINV) THEN
          ALLOCATE(CWORK2(NGL))
          CWORK2(:)=(0.D0,0.D0)
          DO IBH=1,NBH
            IB1=2*IBH-1
            IB2=2*IBH
            CSVAR1=0.5D0*(THIS%EIGVEC(IB1,IB)-CI*THIS%EIGVEC(IB2,IB))
            CSVAR2=0.5D0*(THIS%EIGVEC(IB1,IB)+CI*THIS%EIGVEC(IB2,IB))
            CSVAR2=CONJG(CSVAR2)
            DO IG=1,NGL
              CWORK1(IG)=CWORK1(IG)+THIS%PSI0(IG,IDIM,IBH)*CSVAR1
              CWORK2(IG)=CWORK2(IG)+THIS%PSI0(IG,IDIM,IBH)*CSVAR2
            ENDDO
          ENDDO
          ALLOCATE(CWORK3(NGL))
          CALL PLANEWAVE$INVERTG(NGL,CWORK2,CWORK3)
          DO IG=1,NGL
            CWORK1(IG)=CWORK1(IG)+CWORK3(IG)
          ENDDO
          DEALLOCATE(CWORK2)
          DEALLOCATE(CWORK3)
        ELSE
          DO IB1=1,NB
            CSVAR1=THIS%EIGVEC(IB1,IB)
            DO IG=1,NGL
              CWORK1(IG)=CWORK1(IG)+THIS%PSI0(IG,IDIM,IB1)*CSVAR1
            ENDDO
          ENDDO
        END IF
        NRL=MAP%NRL
        IF(NRL.NE.LEN) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        ALLOCATE(CWORK2(NRL))
        CALL PLANEWAVE$FFT('GTOR',1,NGL,CWORK1,NRL,CWORK2)
        DEALLOCATE(CWORK1)
        IF(EXTPNTR%TIM) THEN
          DO IR=1,NRL
            VAL(IR)=AIMAG(CWORK2(IR))
          ENDDO
        ELSE
          DO IR=1,NRL
            VAL(IR)=REAL(CWORK2(IR),KIND=8)
          ENDDO
        END IF
        DEALLOCATE(CWORK2)
!
!     ==========================================================================
!     ==  GET PROJECTIONS                                                     ==
!     ==  FOR A SPECIFIED KPOINT(IKPT), SPIN (ISPIN), BAND(IB),               ==
!     ==  REAL OR IMAGINARY PART (TIM), ATOM(IAT) AS SPECIFIED IN EXTPTR.     ==
!     ==  THE PROJECTION IS OBTAINED FOR EIGENSTATES                          ==
!     ==========================================================================
      ELSE IF(ID.EQ.'<PSPSI|PRO>') THEN
!       == THE RESULT IS <PRO|PSPSI>. THE CURRENT ID IS MISLEADING! ============
        IKPT=EXTPNTR%IKPT   ! IKPT REFERS TO LOCAL KPOINTS
        IF(IKPT.EQ.0) THEN
          CALL ERROR$MSG('STATE NOT AVAILABLE ON THIS TASK')
          CALL ERROR$MSG('THIS ERROR OCCURS ONLY FOR PARALLEL JOBS')
          CALL ERROR$MSG("USE WAVES$GETL4('AVAILABLESTATE',TCHK) TO EXPLORE")
          CALL ERROR$I4VAL('THISTASK(MOMOMER)',THISTASK)
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('ISPIN',ISPIN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        ISPIN=EXTPNTR%ISPIN
        IB=EXTPNTR%IB
        IF(NSPIN.EQ.1.AND.NDIM.EQ.2) THEN
          IDIM=ISPIN
          ISPIN=1
        ELSE
          IDIM=1
        END IF
        CALL WAVES_SELECTWV(IKPT,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%EIGVEC)) THEN
          CALL ERROR$MSG('EIGENSTATES NOT AVAILABLE')
          CALL ERROR$L4VAL('THAMILTON',THAMILTON)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        TINV=GSET%TINV
        NB=THIS%NB
        NBH=THIS%NBH
        NAT=MAP%NAT
!
        IPRO=0
        DO IAT=1,NAT
          IF(IAT.EQ.EXTPNTR%IAT)IPRO1=IPRO+1
          ISP=MAP%ISP(IAT)
          DO LN=1,MAP%LNX(ISP)
            L=MAP%LOX(LN,ISP)
            IPRO=IPRO+2*L+1
          ENDDO
          IF(IAT.EQ.EXTPNTR%IAT) THEN
            IPRO2=IPRO
            EXIT
          END IF
        ENDDO
        LMNX=IPRO2-IPRO1+1
!
!       == TRANSFORM TO EIGENSTATES ============================================
        ALLOCATE(CWORK1(LMNX))
        CWORK1(:)=(0.D0,0.D0)
        IF(TINV) THEN
          DO IBH=1,NBH
            IB1=2*IBH-1
            IB2=2*IBH
            CSVAR1=0.5D0*(THIS%EIGVEC(IB1,IB)-CI*THIS%EIGVEC(IB2,IB))
            CSVAR2=0.5D0*(THIS%EIGVEC(IB1,IB)+CI*THIS%EIGVEC(IB2,IB))
            DO IPRO=IPRO1,IPRO2
              LMN=IPRO-IPRO1+1
              CWORK1(LMN)=CWORK1(LMN)+THIS%PROJ(IDIM,IBH,IPRO)*CSVAR1 &
                               +CONJG(THIS%PROJ(IDIM,IBH,IPRO))*CSVAR2
            ENDDO
          ENDDO
        ELSE
          DO IB1=1,NB
            CSVAR1=THIS%EIGVEC(IB1,IB)
            DO IPRO=IPRO1,IPRO2
              LMN=IPRO-IPRO1+1
              CWORK1(LMN)=CWORK1(LMN)+THIS%PROJ(IDIM,IB1,IPRO)*CSVAR1
            ENDDO
          ENDDO             
        END IF
        IF(LMNX.NE.LEN) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$MSG('LMNX AND LEN MUST BE EQUAL')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('LMNX',LMNX)
          CALL ERROR$I4VAL('EXTPNTR%IAT',EXTPNTR%IAT)
          CALL ERROR$I4VAL('EXTPNTR%IKPT',EXTPNTR%IKPT)
          CALL ERROR$I4VAL('EXTPNTR%ISPIN',EXTPNTR%ISPIN)
          CALL ERROR$I4VAL('EXTPNTR%IB',EXTPNTR%IB)
          CALL ERROR$L4VAL('EXTPNTR%TIM',EXTPNTR%TIM)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        VAL(:)=0.D0
        IF(EXTPNTR%TIM) THEN
          DO LMN=1,LMNX
            VAL(LMN)=AIMAG(CWORK1(LMN))
          ENDDO
        ELSE
          DO LMN=1,LMNX
            VAL(LMN)=REAL(CWORK1(LMN),KIND=8)
          ENDDO
        END IF
        DEALLOCATE(CWORK1)
!!$PRINT*,'EXTPNTR',EXTPNTR%IAT,EXTPNTR%IB,' VAL ',VAL
!!$PRINT*,'EIGVEC ',THIS%EIGVEC(:,IB)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$GETR8A')
      END IF
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$GETC8A(ID,LEN,VAL)
!     **************************************************************************
!     **  WAVES$GETC8A                                                        **
!     **************************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      COMPLEX(8)  ,INTENT(OUT):: VAL(LEN)
      COMPLEX(8)  ,ALLOCATABLE:: CWORK1(:)
      COMPLEX(8)  ,ALLOCATABLE:: CWORK2(:)
      COMPLEX(8)  ,ALLOCATABLE:: CWORK3(:)
      INTEGER(4)              :: IKPT,ISPIN,IB,IDIM,IBH,IB1,IB2,IG,IR
      INTEGER(4)              :: IPRO,IPRO1,IPRO2,IAT,ISP,L,LN,LMNX,LMN
      INTEGER(4)              :: NB,NBH,NRL,NAT,NGL
      COMPLEX(8)              :: CSVAR1,CSVAR2
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL(4)              :: TINV
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$GETR8A')
!
!     ==========================================================================
!     ==  REAL SPACE PS WAVE FUNCTION                                         ==
!     ==  NOTE: IB,IKPT,ISPIN MUST BE SET ON LINKEDLIST                       ==
!     ==========================================================================
      IF(ID.EQ.'PSPSI') THEN
        IKPT=EXTPNTR%IKPT
        IF(IKPT.EQ.0) THEN
          CALL ERROR$MSG('STATE NOT AVAILABLE ON THIS TASK')
          CALL ERROR$MSG('THIS ERROR OCCURS ONLY FOR PARALLEL JOBS')
          CALL ERROR$MSG("USE WAVES$GETL4('AVAILABLESTATE',TCHK) TO EXPLORE")
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        ISPIN=EXTPNTR%ISPIN
        IB=EXTPNTR%IB
        IF(NSPIN.EQ.1.AND.NDIM.EQ.2) THEN
          IDIM=ISPIN
          ISPIN=1
        ELSE
          IDIM=1
        END IF
        CALL WAVES_SELECTWV(IKPT,ISPIN)
        CALL PLANEWAVE$SELECT(GSET%ID)
        IF(.NOT.ASSOCIATED(THIS%EIGVEC)) THEN
          CALL ERROR$MSG('EIGENSTATES NOT AVAILABLE')
          CALL ERROR$L4VAL('THAMILTON',THAMILTON)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        NGL=GSET%NGL
        NB=THIS%NB
        NBH=THIS%NBH
        TINV=GSET%TINV
        ALLOCATE(CWORK1(NGL))
        CWORK1(:)=(0.D0,0.D0)
!!$WRITE(*,*)"IB=",IB
!!$WRITE(*,FMT='("EIGVEC=",100("(",2F10.5,") "))')THIS%EIGVEC(:,IB)
        IF(TINV) THEN
          ALLOCATE(CWORK2(NGL))
          CWORK2(:)=(0.D0,0.D0)
          DO IBH=1,NBH
            IB1=2*IBH-1
            IB2=2*IBH
            CSVAR1=0.5D0*(THIS%EIGVEC(IB1,IB)-CI*THIS%EIGVEC(IB2,IB))
            CSVAR2=0.5D0*(THIS%EIGVEC(IB1,IB)+CI*THIS%EIGVEC(IB2,IB))
            CSVAR2=CONJG(CSVAR2)
            DO IG=1,NGL
              CWORK1(IG)=CWORK1(IG)+THIS%PSI0(IG,IDIM,IBH)*CSVAR1
              CWORK2(IG)=CWORK2(IG)+THIS%PSI0(IG,IDIM,IBH)*CSVAR2
            ENDDO
          ENDDO
          ALLOCATE(CWORK3(NGL))
          CALL PLANEWAVE$INVERTG(NGL,CWORK2,CWORK3)
          DO IG=1,NGL
            CWORK1(IG)=CWORK1(IG)+CWORK3(IG)
          ENDDO
          DEALLOCATE(CWORK2)
          DEALLOCATE(CWORK3)
        ELSE
          DO IB1=1,NB
            CSVAR1=THIS%EIGVEC(IB1,IB)
            DO IG=1,NGL
              CWORK1(IG)=CWORK1(IG)+THIS%PSI0(IG,IDIM,IB1)*CSVAR1
            ENDDO
          ENDDO
        END IF
        NRL=MAP%NRL
        IF(NRL.NE.LEN) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        CALL PLANEWAVE$FFT('GTOR',1,NGL,CWORK1,NRL,VAL)
        DEALLOCATE(CWORK1)
!
!     ==========================================================================
!     ==  GET PROJECTIONS                                                     ==
!     ==  FOR A SPECIFIED KPOINT(IKPT), SPIN (ISPIN), BAND(IB),               ==
!     ==  REAL OR IMAGINARY PART (TIM), ATOM(IAT) AS SPECIFIED IN EXTPTR.     ==
!     ==  THE PROJECTION IS OBTAINED FOR EIGENSTATES                          ==
!     ==========================================================================
      ELSE IF(ID.EQ.'<PSPSI|PRO>') THEN
        IKPT=EXTPNTR%IKPT   ! IKPT REFERS TO LOCAL KPOINTS
        IF(IKPT.EQ.0) THEN
          CALL ERROR$MSG('STATE NOT AVAILABLE ON THIS TASK')
          CALL ERROR$MSG('THIS ERROR OCCURS ONLY FOR PARALLEL JOBS')
          CALL ERROR$MSG("USE WAVES$GETL4('AVAILABLESTATE',TCHK) TO EXPLORE")
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        ISPIN=EXTPNTR%ISPIN
        IB=EXTPNTR%IB
        IF(NSPIN.EQ.1.AND.NDIM.EQ.2) THEN
          IDIM=ISPIN
          ISPIN=1
        ELSE
          IDIM=1
        END IF
        CALL WAVES_SELECTWV(IKPT,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%EIGVEC)) THEN
          CALL ERROR$MSG('EIGENSTATES NOT AVAILABLE')
          CALL ERROR$L4VAL('THAMILTON',THAMILTON)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        TINV=GSET%TINV
        NB=THIS%NB
        NBH=THIS%NBH
        NAT=MAP%NAT
!
        IPRO=0
        DO IAT=1,NAT
          IF(IAT.EQ.EXTPNTR%IAT)IPRO1=IPRO+1
          ISP=MAP%ISP(IAT)
          DO LN=1,MAP%LNX(ISP)
            L=MAP%LOX(LN,ISP)
            IPRO=IPRO+2*L+1
          ENDDO
          IF(IAT.EQ.EXTPNTR%IAT) THEN
            IPRO2=IPRO
            EXIT
          END IF
        ENDDO
        LMNX=IPRO2-IPRO1+1
!
!       == TRANSFORM TO EIGENSTATES ============================================
        ALLOCATE(CWORK1(LMNX))
        CWORK1(:)=(0.D0,0.D0)
        IF(TINV) THEN
          DO IBH=1,NBH
            IB1=2*IBH-1
            IB2=2*IBH
            CSVAR1=0.5D0*(THIS%EIGVEC(IB1,IB)-CI*THIS%EIGVEC(IB2,IB))
            CSVAR2=0.5D0*(THIS%EIGVEC(IB1,IB)+CI*THIS%EIGVEC(IB2,IB))
            DO IPRO=IPRO1,IPRO2
              LMN=IPRO-IPRO1+1
              CWORK1(LMN)=CWORK1(LMN)+THIS%PROJ(IDIM,IBH,IPRO)*CSVAR1 &
                               +CONJG(THIS%PROJ(IDIM,IBH,IPRO))*CSVAR2
            ENDDO
          ENDDO
        ELSE
          DO IB1=1,NB
            CSVAR1=THIS%EIGVEC(IB1,IB)
            DO IPRO=IPRO1,IPRO2
              LMN=IPRO-IPRO1+1
              CWORK1(LMN)=CWORK1(LMN)+THIS%PROJ(IDIM,IB1,IPRO)*CSVAR1
            ENDDO
          ENDDO             
        END IF
        IF(LMNX.GT.LEN) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        VAL(:)=0.D0
        VAL(:LMNX)=CWORK1(:LMNX)
        DEALLOCATE(CWORK1)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$GETC8A')
      END IF
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$SETI4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
      INTEGER(4)              :: IKPT,IKPTL,I
      INTEGER(4)              :: THISTASK,NTASKS
!     **************************************************************************
      IF(ID.EQ.'SPINORDIM') THEN
        NDIM=VAL
!
      ELSE IF(ID.EQ.'IKPT') THEN   ! GLOBAL IKPT
        IKPT=VAL
!       == CHECK IF K-POINT IS PRESENT =========================================
        CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
        IF(KMAP(IKPT).NE.THISTASK) THEN
          IKPTL=0
        ELSE
!         == CONVERT GLOBAL IKPT INTO LOCAL IKPT ===============================
          IKPTL=0
          DO I=1,IKPT
            IF(KMAP(I).EQ.THISTASK) IKPTL=IKPTL+1
          ENDDO
        END IF
        CALL MPE$BROADCAST('K',1,IKPTL)
        EXTPNTR%IKPT=IKPTL
!
      ELSE IF(ID.EQ.'ISPIN') THEN
        EXTPNTR%ISPIN=VAL
      ELSE IF(ID.EQ.'IB') THEN
        EXTPNTR%IB=VAL
      ELSE IF(ID.EQ.'IAT') THEN
        EXTPNTR%IAT=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$GETI4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ID.EQ.'SPINORDIM') THEN
        VAL=NDIM
      ELSE IF(ID.EQ.'NSPIN') THEN
        VAL=NSPIN
      ELSE IF(ID.EQ.'NDIM') THEN
        VAL=NDIM
      ELSE IF(ID.EQ.'NDIMD') THEN
        VAL=NDIMD
!
      ELSE IF(ID.EQ.'NB') THEN
        IF(EXTPNTR%ISPIN.EQ.0) EXTPNTR%ISPIN=1
        IF(EXTPNTR%IKPT.EQ.0)  EXTPNTR%IKPT=1
        CALL WAVES_SELECTWV(EXTPNTR%IKPT,EXTPNTR%ISPIN)
        VAL=THIS%NB
!
      ELSE IF(ID.EQ.'NKPT') THEN
        VAL=NKPT   !GLOBAL NKPT
!
      ELSE IF(ID.EQ.'NR1') THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR1',VAL)
!
      ELSE IF(ID.EQ.'NR1L') THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR1L',VAL)
!
      ELSE IF(ID.EQ.'NR2') THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR2',VAL)
!
      ELSE IF(ID.EQ.'NR3') THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR3',VAL)
!
      ELSE IF(ID.EQ.'NR1START') THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR1START',VAL)
!
      ELSE IF(ID.EQ.'NND') THEN
        IF(ALLOCATED(OSDENMAT)) THEN
          VAL=SIZE(OSDENMAT)
        ELSE
          VAL=-1
        END IF
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$SETL4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'STOP') THEN
        TSTOP=VAL
      ELSE IF(ID.EQ.'SAFEORTHO') THEN
        TSAFEORTHO=VAL
      ELSE IF(ID.EQ.'SWAPSTATES') THEN
        TSWAPSTATES=VAL
      ELSE IF(ID.EQ.'RANDOMIZE') THEN
        TRANDOM=VAL
!
!     == UNITARY TRANSFORMATION ONTO EIGENSTATES (PERFORMED ONCE) ==============
      ELSE IF(ID.EQ.'STRAIGHTEN') THEN
        TSTRAIGHTEN=VAL
      ELSE IF(ID.EQ.'HAMILTON') THEN
        THAMILTON=VAL
      ELSE IF(ID.EQ.'FORCE') THEN
        TFORCEX=VAL
      ELSE IF(ID.EQ.'STRESS') THEN
        TSTRESSX=VAL
      ELSE IF(ID.EQ.'RAWSTATES') THEN
        TRAWSTATES=VAL
      ELSE IF(ID.EQ.'TIM') THEN
        EXTPNTR%TIM=VAL
      ELSE IF(ID.EQ.'FIXRHO') THEN
        TFIXRHO=VAL
      ELSE IF(ID.EQ.'WRITERHO') THEN
        TWRITERHO=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$GETL4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      LOGICAL(4)  ,INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'RAWSTATES') THEN
        VAL=TRAWSTATES
      ELSE IF(ID.EQ.'TIM') THEN
        VAL=EXTPNTR%TIM
      ELSE IF(ID.EQ.'SAFEORTHO') THEN
        VAL=TSAFEORTHO
      ELSE IF(ID.EQ.'SWAPSTATES') THEN
        VAL=TSWAPSTATES
      ELSE IF(ID.EQ.'WRITERHO') THEN
        VAL=TWRITERHO
      ELSE IF(ID.EQ.'AVAILABLESTATE') THEN
        VAL=(EXTPNTR%IKPT.NE.0)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$GETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$SETCH(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'RSTRTTYPE') THEN
        IF(VAL.NE.'DYNAMIC'.AND.VAL.NE.'STATIC') THEN
          CALL ERROR$MSG('VALUE NOT RECOGNIZED')
          CALL ERROR$MSG('ALLOWED VALUES ARE "STATIC" AND "DYNAMIC"')
          CALL ERROR$CHVAL('VAL',VAL)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$SETCH')
        END IF
        RSTRTTYPE=VAL
      ELSE IF(ID.EQ.'OPTIMIZER') THEN
        IF(VAL.NE.'CG'.AND.VAL.NE.'MD') THEN
          CALL ERROR$MSG('VALUE NOT RECOGNIZED')
          CALL ERROR$MSG('ALLOWED VALUES ARE "CG" AND "MD"')
          CALL ERROR$CHVAL('VAL',VAL)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$SETCH')
        END IF
        OPTIMIZERTYPE=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$SETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$GETRSPACEMATA(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE RSPACEOP_MODULE, ONLY : RSPACEMAT_TYPE &
     &                           ,RSPACEOP$COPY 
      USE WAVES_MODULE, ONLY: OSDENMAT
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)          :: ID
      INTEGER(4)  ,INTENT(IN)          :: LEN
      TYPE(RSPACEMAT_TYPE),INTENT(OUT) :: VAL(LEN)
      INTEGER(4)                       :: I
!     **************************************************************************
      IF(ID.EQ.'DENMAT') THEN
        IF(.NOT.ALLOCATED(OSDENMAT)) THEN
          CALL ERROR$MSG('OSDENMAT IS NOT AVAILABLE')
          CALL ERROR$MSG('CODE ERROR')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETRSPACEMATA')
        END IF
        IF(LEN.NE.SIZE(OSDENMAT)) THEN
          CALL ERROR$MSG('LEN INCONSISTENT WITH SIZE OF OSDENMAT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('SIZE(OSDENMAT)',SIZE(OSDENMAT))
          CALL ERROR$STOP('WAVES$GETRSPACEMATA')
        END IF
        DO I=1,LEN
          CALL RSPACEOP$COPY(OSDENMAT(I),VAL(I))
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$GETRSPACEMATA')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$SETRSPACEMATA(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE RSPACEOP_MODULE, ONLY: RSPACEMAT_TYPE &
     &                          ,RSPACEOP$COPY &
     &                          ,RSPACEOP$WRITEMAT
      USE WAVES_MODULE   , ONLY: OSHAMIL &
     &                          ,OSDENMAT 
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)          :: ID
      INTEGER(4)  ,INTENT(IN)          :: LEN
      TYPE(RSPACEMAT_TYPE),INTENT(IN)  :: VAL(LEN)
      INTEGER(4)                       :: I
!     **************************************************************************
      IF(ID.EQ.'HAMIL') THEN
        IF(.NOT.ALLOCATED(OSDENMAT)) THEN
          CALL ERROR$MSG('OSDENMAT IS NOT AVAILABLE')
          CALL ERROR$MSG('MUST BE PRESENT TO CHECK SIZE CONSISTENCY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$SETRSPACEMATA')
        END IF
        IF(LEN.NE.SIZE(OSDENMAT)) THEN
          CALL ERROR$MSG('LEN INCONSISTENT WITH SIZE OF OSDENMAT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('SIZE(OSDENMAT)',SIZE(OSDENMAT))
          CALL ERROR$STOP('WAVES$SETRSPACEMATA')
        END IF
        IF(.NOT.ALLOCATED(OSHAMIL)) THEN
          ALLOCATE(OSHAMIL(LEN))
        END IF
        DO I=1,LEN
          CALL RSPACEOP$COPY(VAL(I),OSHAMIL(I))
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$SETRSPACEMATA')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$INITIALIZE
!     **************************************************************************
!     **  GENERATE G-VECTORS AND OTHER INITIALIZATION                         **
!     **************************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      REAL(8)                :: RBAS(3,3) ! LATTICE VECTORS
      REAL(8)                :: GBAS(3,3) ! RECIPROCAL LATTICE VECTORS
      REAL(8)                :: CELLVOL   ! UNIT CELL  VOLUME
      REAL(8)   ,ALLOCATABLE :: XK(:,:)   ! K-POINTS IN RELATIVE COORDINATES
      REAL(8)   ,ALLOCATABLE :: GVEC(:,:)     ! G VECTOR
      INTEGER(4)             :: LN,ISP,IKPT,ISPIN,IAT,IKPTG,IKPTL
      INTEGER(4)             :: NB
      INTEGER(4)             :: NBH
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NRL,NR1L,NR1START,NR2,NR3
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: NFILO
      INTEGER(4)             :: LNX
      INTEGER(4)             :: NTASKS,THISTASK
      INTEGER(4),ALLOCATABLE :: ICOLOR(:)
      LOGICAL(4),ALLOCATABLE :: TINVARR(:)
      LOGICAL(4)             :: TKGROUP
      INTEGER(4)             :: I,J,K
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$INITIALIZE')
!     
!     ==========================================================================
!     ==  DETERMINE MAPPRO AND MAPBAREPRO                                     ==
!     ==========================================================================
      CALL ATOMLIST$NATOM(MAP%NAT)
      CALL SETUP$GETI4('NSP',MAP%NSP)   !FORMER CALL SETUP$NSPECIES(MAP%NSP)
      CALL SETUP$GETI4('LNXX',MAP%LNXX) !FORMER CALL SETUP$LNXX(MAP%LNXX)
      ALLOCATE(MAP%ISP(MAP%NAT))
      ALLOCATE(MAP%LNX(MAP%NSP))
      ALLOCATE(MAP%LMNX(MAP%NSP))
      ALLOCATE(MAP%LOX(MAP%LNXX,MAP%NSP))
      CALL ATOMLIST$GETI4A('ISPECIES',0,MAP%NAT,MAP%ISP)
      MAP%LMX=0
      MAP%LOX(:,:)=0
      MAP%NBAREPRO=0
      DO ISP=1,MAP%NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX) !FORMER CALL SETUP$LNX(ISP,LNX)
        MAP%LNX(ISP)=LNX
        CALL SETUP$GETI4A('LOX',LNX,MAP%LOX(1:LNX,ISP))
                       !FORMER CALL SETUP$LOFLN(ISP,LNX,MAP%LOX(1:LNX,ISP))
        MAP%LMNX(ISP)=0
        DO LN=1,LNX
          MAP%LMNX(ISP)=MAP%LMNX(ISP)+2*MAP%LOX(LN,ISP)+1
          MAP%LMX=MAX(MAP%LMX,(MAP%LOX(LN,ISP)+1)**2)
        ENDDO
        MAP%NBAREPRO=MAP%NBAREPRO+LNX
        CALL SETUP$UNSELECT()
      ENDDO
      MAP%NPRO=0
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
        MAP%NPRO=MAP%NPRO+MAP%LMNX(ISP)
      ENDDO
!     
!     ==================================================================
!     ==  K-POINT PARALLELIZATION                                     ==
!     ==================================================================
      CALL DYNOCC$GETI4('NKPT',NKPT) 
      ALLOCATE(TINVARR(NKPT))
      ALLOCATE(XK(3,NKPT))
      CALL DYNOCC$GETR8A('XK',3*NKPT,XK)
      DO IKPT=1,NKPT
        CALL PLANEWAVE$CHECKINVERSION(XK(1,IKPT),TINVARR(IKPT))
        IF(NDIM.NE.1)TINVARR(IKPT)=.FALSE.
      ENDDO
      DEALLOCATE(XK)
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      ALLOCATE(KMAP(NKPT))
      ALLOCATE(ICOLOR(NTASKS))
      CALL WAVES_KDISTRIBUTE(NTASKS,NKPT,TINVARR,ICOLOR,KMAP)
      DEALLOCATE(TINVARR)
      CALL MPE$NEW('MONOMER','K',NTASKS,ICOLOR)
      DEALLOCATE(ICOLOR)
      NKPTL=0
      DO IKPT=1,NKPT
        IF(KMAP(IKPT).EQ.THISTASK) NKPTL=NKPTL+1
      ENDDO
      CALL MPE$BROADCAST('K',1,NKPTL)
!     
!     ==================================================================
!     ==  ALLOCATE THISARRAY                                          ==
!     ==================================================================
      CALL DYNOCC$GETI4('NSPIN',NSPIN)
      ALLOCATE(THISARRAY(NKPTL,NSPIN))
      DO IKPT=1,NKPTL
        ALLOCATE(THISARRAY(IKPT,1)%GSET)
        DO ISPIN=2,NSPIN
          THISARRAY(IKPT,ISPIN)%GSET=>THISARRAY(IKPT,1)%GSET
        ENDDO
        DO ISPIN=1,NSPIN
          NULLIFY(THISARRAY(IKPT,ISPIN)%PSI0)
          NULLIFY(THISARRAY(IKPT,ISPIN)%PSIM)
          NULLIFY(THISARRAY(IKPT,ISPIN)%PROJ)
          THISARRAY(IKPT,ISPIN)%PROJ_ACC_RESIDENT=.FALSE.
          NULLIFY(THISARRAY(IKPT,ISPIN)%TBC_NEW)
          NULLIFY(THISARRAY(IKPT,ISPIN)%HTBC_NEW)
          NULLIFY(THISARRAY(IKPT,ISPIN)%HPSI)
          NULLIFY(THISARRAY(IKPT,ISPIN)%OPSI)
          NULLIFY(THISARRAY(IKPT,ISPIN)%RLAM0)
          NULLIFY(THISARRAY(IKPT,ISPIN)%RLAMM)
          NULLIFY(THISARRAY(IKPT,ISPIN)%RLAM2M)
          NULLIFY(THISARRAY(IKPT,ISPIN)%RLAM3M)
          NULLIFY(THISARRAY(IKPT,ISPIN)%EIGVEC)
          NULLIFY(THISARRAY(IKPT,ISPIN)%EIGVAL)
          NULLIFY(THISARRAY(IKPT,ISPIN)%EXPECTVAL)
        ENDDO
      ENDDO
!     
!     ==================================================================
!     ==  DEFINE REAL SPACE GRID                                      ==
!     ==================================================================
      CALL CELL$GETR8A('TREF',9,RBAS)
      CALL PLANEWAVE$DIVIDERGRIDONTASKS(EPWRHO,RBAS,'MONOMER',NR1START,NR1L,NR2,NR3)
      CALL POTENTIAL$INITIALIZE(EPWRHO,NR1START,NR1L,NR2,NR3)
!
      CALL PLANEWAVE$DIVIDERGRIDONTASKS(EPWRHO,RBAS,'K',NR1START,NR1L,NR2,NR3)
      MAP%NRL=NR1L*NR2*NR3
!     
!     ==================================================================
!     ==  INITIALIZE PLANE WAVE OBJECT                                ==
!     ==================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      IF(NDIM.EQ.1) THEN
        NDIMD=NSPIN
      ELSE
        NDIMD=NDIM**2
      END IF
      IKPTL=0
      DO IKPTG=1,NKPT
        TKGROUP=THISTASK.EQ.KMAP(IKPTG)
        CALL MPE$BROADCAST('K',1,TKGROUP)
        IF(.NOT.TKGROUP) CYCLE
        IKPTL=IKPTL+1
        CALL WAVES_SELECTWV(IKPTL,1)
        WRITE(GSET%ID,FMT=*)IKPTG
        GSET%ID='WAVE '//ADJUSTL(GSET%ID)
        IF(NDIM.EQ.1) THEN
          CALL PLANEWAVE$CHECKINVERSION(XK(1,IKPTL),TINV)
        ELSE
          TINV=.FALSE.
        END IF
!!$CALL FILEHANDLER$UNIT('PROT',NFILO)
!!$WRITE(NFILO,*)'TINV FORCED TO BE FALSE IN WAVES$INITIALIZE!!!'
!!$PRINT*,'TINV FORCED TO BE FALSE IN WAVES$INITIALIZE!!!'
!!$TINV=.FALSE.
!       ==  NOT THAT XK REFERS TO THE LOCAL K-POINT INDEX =================
!PRINT*,THISTASK,'BEFORE ',XK(:,IKPTL),TINV
        CALL PLANEWAVE$INITIALIZE(GSET%ID,'K',RBAS,XK(1,IKPTL),TINV,EPWPSI &
     &                           ,NR1START,NR1L,NR2,NR3)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETL4('TINV',GSET%TINV)
!       == GSPACE AND RSPACE DIMENSIONS ================================
        CALL PLANEWAVE$GETI4('NGL',NGL)
        CALL PLANEWAVE$GETI4('NRL',NRL)
        GSET%NGL=NGL
        NULLIFY(GSET%PRO)
        NULLIFY(GSET%DPRO)
        NULLIFY(GSET%YLM)
        NULLIFY(GSET%SYLM)
        NULLIFY(GSET%MPSI)
!PB070802        NULLIFY(GSET%DMPSI)
        NULLIFY(GSET%BUCKET)
        NULLIFY(GSET%DBUCKET)
      ENDDO
      DEALLOCATE(XK)
      CALL WAVES_UPDATEGSET()
!
!     ==================================================================
!     ==  EVALUATE NUMBER OF G-VECTORS ETC. ON LOCAL PROCESSOR        ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NB)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
!         == #(BANDS) ==================================================
          IF(GSET%TINV) THEN
            NBH=(NB+1)/2
          ELSE
            NBH=NB
          END IF
          THIS%NB=NB
          THIS%NBH=NBH
!         == ALLOCATE PROJECTIONS ======================================
          ALLOCATE(THIS%PROJ(NDIM,NBH,MAP%NPRO))
          THIS%PROJ_ACC_RESIDENT=.FALSE.
!         == ALLOCATE WAVE FUNCTIONS ===================================
          NGL=GSET%NGL
          ALLOCATE(THIS%PSI0(NGL,NDIM,NBH))
          ALLOCATE(THIS%PSIM(NGL,NDIM,NBH))
          THIS%PSI0(:,:,:)=(0.D0,0.D0)
          THIS%PSIM(:,:,:)=(0.D0,0.D0)
          ALLOCATE(GVEC(3,NGL))
          CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
          CALL WAVES_INITIALIZERANDOM(NGL,NDIM,NBH,GVEC,THIS%PSI0)
          DEALLOCATE(GVEC)
          DO I=1,NBH
            DO J=1,NDIM
              DO K=1,NGL
                THIS%PSIM(K,J,I)=THIS%PSI0(K,J,I)
              ENDDO
            ENDDO
          ENDDO
!          THIS%PSIM=THIS%PSI0 ! THIS STATEMENT GAVE PROBLEMS WITH IFC10
!         == ALLOCATE LAGRANGE MULTIPLYERS =============================
          NULLIFY(THIS%HPSI)
          NULLIFY(THIS%OPSI)
          ALLOCATE(THIS%RLAM0(NB,NB))
          THIS%RLAM0(:,:)=(0.D0,0.D0)
          NULLIFY(THIS%RLAMM)
          NULLIFY(THIS%RLAM2M)
          NULLIFY(THIS%RLAM3M)
          NULLIFY(THIS%EIGVAL)
          NULLIFY(THIS%EIGVEC)
          NULLIFY(THIS%EXPECTVAL)
        ENDDO
      ENDDO
!     
!     ==================================================================
!     ==  REPORT DATA ABOUT FFTS                                      ==
!     ==================================================================
!!$      CALL FILEHANDLER$UNIT('PROT',NFILO)
!!$      CALL PLANEWAVE$REPORT(NFILO)
!
!      DO IKPT=1,NKPTL
!        CALL WAVES_SELECTWV(IKPT,1)
!        CALL PLANEWAVE$SELECT(GSET%ID)
!        CALL PLANEWAVE$REPORT(NFILO)
!      ENDDO
!     
!     ================================================================
!     ==  SEND DATA TO OPTIC CODE                                   ==
!     ================================================================
!      CALL OPTICS$GETL4('ON',TCHK)
!      IF(TCHK) THEN
!        ALLOCATE(XK(3,NKPT))
!        CALL DYNOCC$GETR8A('XK',3*NKPT,XK)
!        TCHK=.FALSE.
!        DO IKPT=1,NKPT
!          TCHK=XK(1,IKPT).EQ.0.AND.XK(2,IKPT).EQ.0.AND.XK(3,IKPT).EQ.0
!          IF(TCHK) THEN
!            CALL WAVES_SELECTWV(IKPT,1)
!            CALL PLANEWAVE$SELECT(GSET%ID)
!            CALL OPTICS$LATTICE
!            EXIT
!          END IF
!        ENDDO
!        IF(.NOT.TCHK) THEN
!          CALL ERROR$MSG('OPTICS CODE REQUIRES GAMMA POINT IN THE K-POINT SET')
!          CALL ERROR$STOP('WAVES$INITIALIZE')
!        END IF
!        DEALLOCATE(XK)
!      END IF
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_DYNOCCGETR8A(ID,LEN,VAL)
!     **                                                              **
!     **  GET A REAL(8) ARRAY FROM DYNOCC. THIS ROUTINE IS AN         **
!     **  INTERFACE BETWEEN WAVES OBJECT AND DYNOCC OBJECT.           **
!     **  IT IS REQUIRED FOR THE K-POINT PARALLELIZATION, BECAUSE     **
!     **  EACH NODE ONLY MAINTAINS A FRACTION OF ALL K-POINTS         **
!     **  WHILE THE DYNOCC OBJECT WORKS ON ALL K-POINTS               **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   ASSUMES THAT THE INPUT ARRAY IS IDENTICAL ON ALL           **
!     **      TASKS OF THE LOCAL K-GROUP                              **
!     **   KMAP(NKPT) POINTS TO THE FIRST TASK IN THE K-GROUP, WHERE  **
!     **     THE K-POINT RESIDES (THE TASK IS RELATIVE TO THE MONOMER **
!     **     GROUP)                                                   **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE MPE_MODULE
      USE WAVES_MODULE  ,ONLY : NKPT,NKPTL,KMAP
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
      REAL(8)     ,ALLOCATABLE:: VALG(:)
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: IKPTL,IKPT
      INTEGER(4)              :: NSPIN,NBX
      INTEGER(4)              :: ISVAR1L,ISVAR2L,ISVAR1G,ISVAR2G
      INTEGER(4)              :: ISPIN
!     ******************************************************************
      IF(ID.EQ.'OCC'.OR.ID.EQ.'EPSILON') THEN
        CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
        CALL DYNOCC$GETI4('NB',NBX)
        CALL DYNOCC$GETI4('NSPIN',NSPIN)
        ALLOCATE(VALG(NBX*NKPT*NSPIN))
        CALL DYNOCC$GETR8A(ID,NBX*NKPT*NSPIN,VALG)
        IKPTL=0
        DO IKPT=1,NKPT
           IF(KMAP(IKPT).EQ.THISTASK) THEN
            IKPTL=IKPTL+1
            DO ISPIN=1,NSPIN
!             == DIMENSIONS: OCC(NBX,NKPT,NSPIN)
              ISVAR1L=1+NBX*(IKPTL-1+NKPTL*(ISPIN-1))
              ISVAR2L=ISVAR1L+NBX-1
              ISVAR1G=1+NBX*(IKPT-1+NKPT*(ISPIN-1))
              ISVAR2G=ISVAR1G+NBX-1
              VAL(ISVAR1L:ISVAR2L)=VALG(ISVAR1G:ISVAR2G)
            ENDDO
          END IF
        ENDDO
        DEALLOCATE(VALG)
        CALL MPE$BROADCAST('K',1,VAL)
!
      ELSE IF(ID.EQ.'XK') THEN
        CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
        ALLOCATE(VALG(3*NKPT))
        CALL DYNOCC$GETR8A('XK',3*NKPT,VALG)
        IKPTL=0
        DO IKPT=1,NKPT
          IF(KMAP(IKPT).EQ.THISTASK) THEN
            IKPTL=IKPTL+1
            ISVAR1L=1+3*(IKPTL-1)
            ISVAR2L=ISVAR1L+3-1
            ISVAR1G=1+3*(IKPT-1)
            ISVAR2G=ISVAR1G+3-1
            VAL(ISVAR1L:ISVAR2L)=VALG(ISVAR1G:ISVAR2G)
          END IF
        ENDDO                        
        DEALLOCATE(VALG)
        CALL MPE$BROADCAST('K',1,VAL)
!
      ELSE IF(ID.EQ.'WKPT') THEN
!       == THIS OPTION IS USED IN DMFT MODULE  =================================
        CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
        ALLOCATE(VALG(NKPT))
        CALL DYNOCC$GETR8A('WKPT',NKPT,VALG)
        IKPTL=0
        DO IKPT=1,NKPT
          IF(KMAP(IKPT).EQ.THISTASK) THEN
            IKPTL=IKPTL+1
            VAL(IKPTL)=VALG(IKPT)
          END IF
        ENDDO                         
        DEALLOCATE(VALG)
        CALL MPE$BROADCAST('K',1,VAL)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES_DYNOCCGETR8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_DYNOCCSETR8A(ID,LEN,VAL)
!     **************************************************************************
!     **  SET A REAL(8) ARRAY TO DYNOCC. THIS ROUTINE IS AN                   **
!     **  INTERFACE BETWEEN WAVES OBJECT AND DYNOCC OBJECT.                   **
!     **  IT IS REQUIRED FOR THE K-POINT PARALLELIZATION, BECAUSE             **
!     **  EACH NODE ONLY MAINTAINS A FRACTION OF ALL K-POINTS                 **
!     **  WHILE THE DYNOCC OBJECT WORKS ON ALL K-POINTS                       **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **   ASSUMES THAT THE INPUT ARRAY IS IDENTICAL ON ALL                   **
!     **      TASKS OF THE LOCAL K-GROUP                                      **
!     **   KMAP(NKPT) POINTS TO THE FIRST TASK IN THE K-GROUP WHERE           **
!     **     THE K-POINT RESIDES (THE TASK IS RELATIVE TO THE MONOMER         **
!     **     GROUP)                                                           **
!     **                                                                      **
!     **                                                                      **
!     *********************P.E. BLOECHL, TU-CLAUSTHAL (2005)********************
      USE MPE_MODULE
      USE WAVES_MODULE  ,ONLY : NKPT,NKPTL,KMAP
!     -- KMAP(NKPT) POINTS TO THE FIRST TASK IN THE K-GROUP WHERE THE 
!     --    K-POINT RESIDES (THE TASK IS RELATIVE TO THE MONOMER GROUP)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: VAL(LEN)
      REAL(8)     ,ALLOCATABLE:: VALG(:)
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: IKPTL,IKPT
      INTEGER(4)              :: NSPIN,NBX
      INTEGER(4)              :: ISVAR1L,ISVAR2L,ISVAR1G,ISVAR2G
      INTEGER(4)              :: ISPIN
!     **************************************************************************
      IF(ID.EQ.'EIG'.OR.ID.EQ.'M<PSIDOT|PSIDOT>'.OR.ID.EQ.'EPSILON') THEN
        CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
        CALL DYNOCC$GETI4('NB',NBX)
        CALL DYNOCC$GETI4('NSPIN',NSPIN)
        ALLOCATE(VALG(NBX*NKPT*NSPIN))
        VALG(:)=0.D0
        IKPTL=0
        DO IKPT=1,NKPT
          IF(KMAP(IKPT).EQ.THISTASK) THEN
            IKPTL=IKPTL+1
            DO ISPIN=1,NSPIN
!             == DIMENSIONS: EIG(NBX,NKPT,NSPIN)
              ISVAR1L=1+NBX*(IKPTL-1+NKPTL*(ISPIN-1))
              ISVAR2L=ISVAR1L+NBX-1
              ISVAR1G=1+NBX*(IKPT-1+NKPT*(ISPIN-1))
              ISVAR2G=ISVAR1G+NBX-1
              VALG(ISVAR1G:ISVAR2G)=VAL(ISVAR1L:ISVAR2L)
            ENDDO
          END IF
        ENDDO
        CALL MPE$COMBINE('MONOMER','+',VALG)
        CALL DYNOCC$SETR8A(ID,NBX*NKPT*NSPIN,VALG)
        DEALLOCATE(VALG)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES_DYNOCCSETR8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$ETOT()
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATE PS KINETIC ENERGY IN G-SPACE                               **
!     **  EVALUATE NUMBER OF ELECTRONS IN G-SPACE                             **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (2000)***********
      USE MPE_MODULE
      USE WAVES_MODULE
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &       CPPAW_CUBLAS_ACC_HPSI_RESIDENCY_ENABLED &
     &      ,CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_RESIDENT &
     &      ,CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_FORCE &
     &      ,CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D
#ENDIF
      IMPLICIT NONE
      INTEGER(4)             :: LMRXX
      REAL(8)   ,ALLOCATABLE :: QLM(:,:)        !(LMRXX) MULTIPOLE MOMENTS
      REAL(8)   ,ALLOCATABLE :: VQLM(:,:)       !(LMRXX) MULTIPOLE POTENTIALS
      REAL(8)   ,ALLOCATABLE :: RHO(:,:)        ! CHARGE DENSITY
      REAL(8)   ,ALLOCATABLE :: RHOKIN(:,:)     ! KINETIC ENERGY DENSITY
      COMPLEX(8),ALLOCATABLE :: DENMAT(:,:,:,:) ! 1CENTER DENSITY MATRIX
      COMPLEX(8),ALLOCATABLE :: EDENMAT(:,:,:,:)! ENERGY-WEIGHTED DENSITY MATRIX
      COMPLEX(8),ALLOCATABLE :: DH(:,:,:,:)     ! 1CENTER HAMILTONIAN
      COMPLEX(8),ALLOCATABLE :: DH1(:,:,:,:)    ! 1CENTER HAMILTONIAN
      REAL(8)   ,ALLOCATABLE :: DO(:,:,:,:)     ! 1CENTER OVERLAP
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      REAL(8)   ,ALLOCATABLE :: R(:,:)
      REAL(8)   ,ALLOCATABLE :: FORCE(:,:)
      REAL(8)   ,ALLOCATABLE :: FORCET(:,:)
      REAL(8)                :: STRESS1(3,3),STRESS(3,3)
      REAL(8)                :: STRESSKIN(3,3)
      REAL(8)                :: RBAS(3,3) ! REAL SPACE LATTICE VECTORS
      REAL(8)                :: RHOB      ! BACKGROUND DENSITY
      REAL(8)                :: POTB      ! AVERAGE ELECTROSTATIC POTENTIAL
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NRL
      REAL(8)                :: EKIN
      INTEGER(4)             :: IKPT,ISPIN,IAT,ISP,IB,IR
      INTEGER(4)             :: L,M,LN
      INTEGER(4)             :: LMN1,LMN2
      INTEGER(4)             :: LMNXX
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NBH,NBX,NB
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: ISVAR
      REAL(8)                :: SVAR1,SVAR2
      COMPLEX(8)             :: CSVAR1,CSVAR2
      COMPLEX(8),ALLOCATABLE :: HAMILTON(:,:)
      REAL(8)   ,ALLOCATABLE :: EIG(:,:,:)
      LOGICAL(4)             :: TSTRESS
      LOGICAL(4)             :: TFORCE
      LOGICAL(4)             :: TCHK
      COMPLEX(8),ALLOCATABLE :: QMAT(:,:)   
      INTEGER(4)             :: NFILO
      LOGICAL(4)             :: TCONV ! MIXER SAYS THAT WAVE FUNCTIONS ARE CONVERGED !KAESTNERCG
      REAL(8)                :: CONVPSI ! CONVERGENCE CRITERION FOR WAVE FUNCTIONS !KAESTNERCG
      LOGICAL(4)             :: TRHOKIN=.FALSE. !KINETIC ENERGY DENSITY REQUIRED
      INTEGER(4) ::NTASKS_W,THISTASK_W
      REAL(8) :: RBASM(3,3)
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      LOGICAL(4)             :: THPSIRESIDENT
#ENDIF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_ETOT_T0
      REAL(8)                :: ACCEL_ETOT_T1
#ENDIF
!     **************************************************************************
      CALL MPE$QUERY('~',NTASKS_W,THISTASK_W)
                              CALL TRACE$PUSH('WAVES$ETOT')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == CHECK CONSISTENCY WITH OCCUPATIONS OBJECT                            ==
!     ==========================================================================
      CALL DYNOCC$GETL4('DYN',TCHK)
      IF(TCHK) THEN
        IF(TSAFEORTHO) THEN
          CALL ERROR$MSG('FOR DYNAMICAL OCCUPATIONS, SAFEORTHO MUST BE FALSE')
          CALL ERROR$STOP('WAVES$ETOT')
        END IF
      END IF
!
!     ==========================================================================
!     == SWITCHES FOR FORCE AND STRESS CALCULATION                            ==
!     ==========================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      CALL ATOMS$GETL4('MOVE',TFORCE)
      TSTRESS=TSTRESS.OR.TSTRESSX
      TFORCE=TFORCE.OR.TFORCEX
      IF(TSTRESS)CALL POTENTIAL$SETL4('STRESS',TSTRESS)
      IF(TFORCE)CALL POTENTIAL$SETL4('FORCE',TFORCE)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SETUP_FLAGS' &
     &    ,0_8,0_8,0_8,0_8,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == COLLECT VARIABLES                                                    ==
!     ==========================================================================
      NAT=MAP%NAT
!     == UNIT CELL =============================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      STRESS=0.D0
!     == ATOMIC POSITIONS ======================================================
      ALLOCATE(R(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
      ALLOCATE(FORCE(3,NAT))
      FORCE(:,:)=0.D0
!     == NUMBER OF BANDS =======================================================
      CALL DYNOCC$GETI4('NB',NBX)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SETUP_COLLECT' &
     &    ,INT(NAT,KIND=8),INT(NBX,KIND=8),0_8,0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == INITIALIZE GSET: YLM, PRO, EIGR                                      ==
!     ==========================================================================
      IF(TFIRST.OR.TSTRESS) THEN
        CALL WAVES_UPDATEGSET()
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SETUP_GSET' &
     &    ,INT(NAT,KIND=8),INT(NBX,KIND=8),0_8,0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == DEALLOCATE EIGENVALUES AND EIGENVECTORS                              ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(ASSOCIATED(THIS%EIGVAL))DEALLOCATE(THIS%EIGVAL)
          IF(ASSOCIATED(THIS%EIGVEC))DEALLOCATE(THIS%EIGVEC)
        ENDDO
      ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SETUP_DEALLOC' &
     &    ,INT(NKPTL,KIND=8),INT(NSPIN,KIND=8),INT(NBX,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
      IF(OPTIMIZERTYPE.EQ.'CG') THEN   !KAESTNERCG
        CALL WAVES$KAESTNERCG1(TFIRST)
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SETUP_CG' &
     &    ,INT(NBX,KIND=8),INT(NKPTL,KIND=8),INT(NSPIN,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == RANDOMIZE INITIAL WAVE FUNCTIONS                                     ==
!     ==========================================================================
      IF(TFIRST.AND.TRANDOM) THEN
        CALL WAVES$RANDOMIZE()
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SETUP_RANDOM' &
     &    ,INT(NBX,KIND=8),INT(NKPTL,KIND=8),INT(NSPIN,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == GRAMM-SCHMIDT ORTHOGONALIZATION OF INITIAL WAVE FUNCTIONS            ==
!     == ATTENTION THIS REORTHOGONALIZATION MAY GIVE THE                      ==
!     == AN UNCONTROLLED KICK AFTER RESTARTING                                ==
!     ==========================================================================
      IF(TFIRST) THEN
        CALL WAVES$GRAMMSCHMIDT()
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SETUP_GRAM' &
     &    ,INT(NBX,KIND=8),INT(NKPTL,KIND=8),INT(NSPIN,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == CALCULATE PROJECTIONS                                                ==
!     ==========================================================================
      IF(TFIRST) THEN
        CALL WAVES$PROJECTIONS('PSI0')
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SETUP_PROJECTIONS' &
     &    ,INT(NBX,KIND=8),INT(NKPTL,KIND=8),INT(NSPIN,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == CALCULATE LMTO STRUCTURE CONSTANTS                                   ==
!     ==========================================================================
      IF(TFIRST.OR.TFORCE.OR.TSTRESS) THEN
                               CALL TIMING$CLOCKON('STRUCTURECONSTANTS')
!IF(TFIRST) THEN
!!$        CALL LMTO$MAKESTRUCTURECONSTANTS()
!END IF
                               CALL TIMING$CLOCKOFF('STRUCTURECONSTANTS')
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SETUP_LMTO' &
     &    ,INT(NAT,KIND=8),0_8,0_8,0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
                               CALL TIMING$CLOCKON('WAVES$ETOT')
!!$                               CALL TIMING$CLOCKON('WAVES$TONTBO')
!!$      CALL WAVES$TONTBO()
!!$                               CALL TIMING$CLOCKOFF('WAVES$TONTBO')
      IF(TFIRST) TFIRST=.FALSE.
!
!     ==========================================================================
!     == KINETIC ENERGY                                                       ==
!     ==========================================================================
      CALL WAVES$EKIN(EKIN,STRESSKIN)   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_EKIN' &
     &    ,0_8,0_8,0_8,0_8,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
      STRESS(:,:)=STRESS(:,:)+STRESSKIN(:,:)
      CALL ENERGYLIST$SET('PS  KINETIC',EKIN)
      CALL ENERGYLIST$ADD('AE  KINETIC',EKIN)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EKIN)
IF(1.EQ.0) THEN
      CALL CELL$GETR8A('TM',9,RBASM)
WRITE(*,FMT='("STRESSTEST:",66("="))')
WRITE(*,FMT='("STRESSTEST: EKIN=",F20.10)')EKIN
WRITE(*,FMT='("STRESSTEST: STRESS*RBAS=",F20.10)')SUM(RBAS*STRESSKIN)
WRITE(*,FMT='("STRESSTEST: RBAS1=",3F20.10)')RBAS(1,:)
WRITE(*,FMT='("STRESSTEST: RBAS2=",3F20.10)')RBAS(2,:)
WRITE(*,FMT='("STRESSTEST: RBAS3=",3F20.10)')RBAS(3,:)
WRITE(*,FMT='("STRESSTEST: STRESS1=",3F20.10)')STRESSKIN(1,:)
WRITE(*,FMT='("STRESSTEST: STRESS2=",3F20.10)')STRESSKIN(2,:)
WRITE(*,FMT='("STRESSTEST: STRESS3=",3F20.10)')STRESSKIN(3,:)
END IF
!
!     ==========================================================================
!     == ONE-CENTER DENSITY MATRICES                                          ==
!     == SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                              ==
!     == SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)                       ==
!     == NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z)         ==
!     ==========================================================================
      NAT=MAP%NAT
      LMNXX=0
      DO ISP=1,MAP%NSP
        LMNXX=MAX(LMNXX,MAP%LMNX(ISP))
      ENDDO
      ALLOCATE(DENMAT(LMNXX,LMNXX,NDIMD,NAT))
      ALLOCATE(EDENMAT(LMNXX,LMNXX,NDIMD,NAT))
      CALL WAVES$DENMAT(LMNXX,NDIMD,NAT,DENMAT,EDENMAT) !<<<<<<<<<<<<<<<<<<<<<<<
      CALL WAVES$OFFSITEDENMAT()  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_DENMAT' &
     &    ,INT(LMNXX,KIND=8),INT(NDIMD,KIND=8),INT(NAT,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == PSEUDO DENSITY STILL WITHOUT PSEUDOCORE                              ==
!     == SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                              ==
!     == SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)                       ==
!     == NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z)         ==
!     ==========================================================================
      NRL=MAP%NRL
      ALLOCATE(RHO(NRL,NDIMD))
      ALLOCATE(RHOKIN(NRL,NDIMD))
      CALL WAVES$RHO(NRL,NDIMD,RHO,TRHOKIN,RHOKIN)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_RHO' &
     &    ,INT(NRL,KIND=8),INT(NDIMD,KIND=8),0_8,0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     ==========================================================================
IF(1.EQ.0) THEN
  SVAR1=0.D0
  DO IB=1,NRL
    SVAR1=SVAR1+RHO(IB,1)
  ENDDO
  CALL PLANEWAVE$GETR8('RWEIGHT',SVAR2)
  SVAR1=SVAR1*SVAR2
  CALL MPE$COMBINE('MONOMER','+',SVAR1)   !LADUNGSDICHTE EINSAMMELN
  PRINT*,'TOTAL CHARGE IN PSEUDO WAVE FUNCTIONS W/O PSCORE ',SVAR1
END IF
!
!     ==========================================================================
!     == MIX RHO AND CHECK IF CONVERGED                                       ==
!     ==========================================================================
      IF(OPTIMIZERTYPE.EQ.'CG') THEN                                 !KAESTNERCG
         CALL TIMING$CLOCKON('MIXER')                                !KAESTNERCG
         TSTOP=.TRUE. ! DO NOT PROPAGATE THE WAVE FUNCTION           !KAESTNERCG
              ! MUST BE SET IN EACH ITERATION, SEE WAVES$PROPAGATE   !KAESTNERCG
 PRINT*,'----------------- NEW CG SCF ITERATION------------------'   !KAESTNERCG
CALL ERROR$MSG('MIXER SWITCHED OFF BECAUSE IT IS STILL INCOMPATIBLE')
CALL ERROR$MSG('WITH A COMPLEX DENSITY MATRIX')
CALL ERROR$STOP('WAVES$ETOT')
!         CALL WAVES_MIXRHO_PUL(NRL*NDIM,LMNXX*LMNXX*NDIMD*NAT,&     !KAESTNERCG
!              RHO,REAL(DENMAT,KIND=8),TCONV,CONVPSI)                !KAESTNERCG
         CALL TIMING$CLOCKOFF('MIXER')                               !KAESTNERCG
      END IF                                                         !KAESTNERCG
!
!     ==========================================================================
!     == MULTIPOLE MOMENTS OF ONE-CENTER PART                                 ==
!     == CONTAINS CORE AND PSEUDOCORE                                         ==
!     ==========================================================================
      CALL SETUP$GETI4('LMRXX',LMRXX) !FORMER CALL SETUP$LMRXX(LMRXX)
      ALLOCATE(QLM(LMRXX,NAT))
      NAT=MAP%NAT
      CALL WAVES$MOMENTS(NAT,LMRXX,NDIMD,LMNXX,DENMAT,QLM)
!
!     ==========================================================================
!     == PREPARE SNAPSHOT OF SPIN TRAJECTORY                                  ==
!     ==========================================================================
!
!     ==========================================================================
!     == OVERWRITE DENSITY TO FOR FIXED POTENTIAL CALC.                       ==
!     ==========================================================================
      IF(TFIXRHO) THEN
        CALL WAVES_FIXRHOGET(NRL,NDIMD,LMRXX,NAT,QLM,RHO,DENMAT)
      END IF
      IF(TWRITERHO) THEN
        CALL WAVES_FIXRHOSET(NRL,NDIMD,LMRXX,NAT,QLM,RHO,DENMAT)
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_MOMENTS' &
     &    ,INT(LMRXX,KIND=8),INT(NAT,KIND=8),INT(LMNXX,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == ANALYSE SPIN DENSITY                                                 ==
!     ==========================================================================
!     __I SUSPECT THAT THIS ROUTINE IS RESPONSIBLE FOR A SEGMENTATION FAULT.____
!     __THEREFORE IT IS TEMPORARILY DISABLED.___________________________________
!      CALL WAVES$SPINS(NRL,NDIMD,RHO,NAT,LMNXX,DENMAT)
!
!     ==========================================================================
!     == POTENTIAL (POTENTIAL IS STORED BACK INTO THE DENSITY ARRAY!)         ==
!     ==========================================================================
      ALLOCATE(VQLM(LMRXX,NAT))
      CALL WAVES_VOFRHO(NRL,NDIMD,RHO,RHOB,NAT,LMRXX,QLM,VQLM)
      DEALLOCATE(QLM)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_VOFRHO' &
     &    ,INT(NRL,KIND=8),INT(NDIMD,KIND=8),INT(NAT,KIND=8) &
     &    ,INT(LMRXX,KIND=8),0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF

!      ALLOCATE(FORCET(3,NAT))
!      ALLOCATE(VQLM(LMRXX,NAT))
!      CALL POTENTIAL$VOFRHO(NRL,NDIMD,RHO,LMRXX,NAT,QLM,VQLM &
!     &                     ,R,FORCET,RBAS,STRESST,RHOB)
!      FORCE=FORCE+FORCET
!      STRESS=STRESS+STRESST
!      DEALLOCATE(QLM)
!      DEALLOCATE(FORCET)
!
!     ==========================================================================
!     == AUGMENTATION                                                         ==
!     ==========================================================================
      ALLOCATE(DH(LMNXX,LMNXX,NDIMD,NAT))
      ALLOCATE(DO(LMNXX,LMNXX,NDIMD,NAT))
      CALL WAVES$SPHERE(LMNXX,NDIMD,NAT,LMRXX,RHOB,DENMAT,EDENMAT &
     &                 ,VQLM,DH,DO,POTB)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SPHERE' &
     &    ,INT(LMNXX,KIND=8),INT(NDIMD,KIND=8),INT(NAT,KIND=8) &
     &    ,INT(LMRXX,KIND=8),0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == INTERFACE TO NTBO BASIS                                              ==
!     == INPUT:   OSDENMAT (OFFSITE-DENSITY MATRIX)                           ==
!     == RETURNS: OSHAMIL (OFFSITE-HAMILTONIAN)                               ==
!     ==========================================================================
!!$      ALLOCATE(DH1(LMNXX,LMNXX,NDIMD,NAT))
!!$      DH1=(0.D0,0.D0)
!!$      CALL LMTO$ETOT(LMNXX,NDIMD,NAT,DENMAT,DH1)
!!$      DH=DH+DH1
!!$      DEALLOCATE(DH1)
!
      CALL SIMPLELMTO$ETOT()
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_SIMPLELMTO' &
     &    ,0_8,0_8,0_8,0_8,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     ==  SUBTRACT AVERAGE ELECTROSTATIC POTENTIAL                            ==
!     ==  TO ACCOUNT FOR BACKGROUND DENSITY                                   ==
!     ==  POTB=-1/OMEGA*INT D^3R(\TILDE{V}+V^1-\TILDE{V}^1)                   ==
!     ==  NOTE THAT \INT D^3R \TILDE{V}=0                                     ==
!     ==========================================================================
      CALL WAVES$ADDCONSTANTPOT(NRL,LMNXX,NDIMD,NAT,POTB,RHO,DH,DO)
      CALL GRAPHICS$SETR8('POTSHIFT',POTB)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_ADDCONSTANT' &
     &    ,INT(NRL,KIND=8),INT(LMNXX,KIND=8),INT(NDIMD,KIND=8) &
     &    ,INT(NAT,KIND=8),0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == COMMUNICATE DATA WITH OPTICS MODULE                                  ==
!     ==========================================================================
!     CALL OPTICS$GETL4('ON',TCHK)
!     IF(TCHK) THEN
!        CALL MPE$COMBINE('NONE','+',DO) ! DO CURRENTLY USED ONLY FOR OPTICS
!        IF(NDIM.EQ.2) THEN
!          CALL ERROR$MSG('OPTICS CODE DOES NOT WORK WITH NONCOLLINEAR MODE')
!          CALL ERROR$STOP('WAVES$ETOT') 
!        END IF
!        CALL OPTICS$WRITEOUT_PARTIALS(NSPIN,NAT,LMNXX,DH,DO)
!      END IF
!
!     ==========================================================================
!     == PRINTOUT FOR TESTING                                                 ==
!     ==========================================================================
!!$IF(NSPIN.EQ.2) THEN
!!$  DO IAT=1,NAT
!!$    ISP=MAP%ISP(IAT)
!!$    LMNX=MAP%LMNX(ISP)
!!$    DO ISPIN=1,NDIMD
!!$      WRITE(*,FMT='("D ======== ATOM ",I2," SPIN ",I2," =========")')IAT,ISPIN
!!$      LMN1=0.D0
!!$      DO LN=1,MAP%LNX(ISP)
!!$        L=MAP%LOX(LN,ISP)
!!$        IF(L.NE.2) THEN
!!$          LMN1=LMN1+2*L+1
!!$          CYCLE
!!$        END IF
!!$        WRITE(*,FMT='("=",I1,10E10.3)') &
!!$             LN,(REAL(DENMAT(LMN1+M,LMN1+M,ISPIN,IAT),KIND=8),M=1,2*L+1)
!!$        LMN1=LMN1+2*L+1
!!$      ENDDO
!!$    ENDDO
!!$    DO ISPIN=1,NDIMD
!!$      WRITE(*,FMT='("H======== ATOM",I2," SPIN ",I2," =========")')IAT,ISPIN
!!$      LMN1=0.D0
!!$      DO LN=1,MAP%LNX(ISP)
!!$        L=MAP%LOX(LN,ISP)
!!$        IF(L.NE.2) THEN
!!$          LMN1=LMN1+2*L+1
!!$          CYCLE
!!$        END IF
!!$        WRITE(*,FMT='("=",I2,10E10.3)') &
!!$             LN,(REAL(DH(LMN1+M,LMN1+M,ISPIN,IAT)),M=1,2*L+1)
!!$        LMN1=LMN1+2*L+1
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO 
!!$ELSE
!!$  DO IAT=1,NAT
!!$    ISP=MAP%ISP(IAT)
!!$    LMNX=MAP%LMNX(ISP)
!!$    DO ISPIN=1,NDIMD
!!$!      WRITE(*,FMT='("D======== ATOM ",I2," SPIN ",I2," =========")')IAT,ISPIN
!!$      LMN1=0.D0
!!$      DO LN=1,MAP%LNX(ISP)
!!$        L=MAP%LOX(LN,ISP)
!!$        IF(L.NE.2) THEN
!!$          LMN1=LMN1+2*L+1
!!$          CYCLE
!!$        END IF
!!$!        WRITE(*,FMT='("=",I2,10E10.3)') &
!!$!             LN,(REAL(DENMAT(LMN1+M,LMN1+M,ISPIN,IAT),KIND=8),M=1,2*L+1)
!!$        LMN1=LMN1+2*L+1
!!$      ENDDO
!!$    ENDDO
!!$    DO ISPIN=1,NDIMD
!!$!      WRITE(*,FMT='("H ======== ATOM ",I2," SPIN ",I2," =========")')IAT,ISPIN
!!$      LMN1=0.D0
!!$      DO LN=1,MAP%LNX(ISP)
!!$        L=MAP%LOX(LN,ISP)
!!$        IF(L.NE.2) THEN
!!$!          LMN1=LMN1+2*L+1
!!$!          CYCLE
!!$        END IF
!!$!        WRITE(*,FMT='("=",I2,10E10.3)') &
!!$!             LN,(REAL(DH(LMN1+M,LMN1+M,ISPIN,IAT)),M=1,2*L+1)
!!$        LMN1=LMN1+2*L+1
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO 
!!$END IF
      DEALLOCATE(VQLM)
      DEALLOCATE(DENMAT)
      DEALLOCATE(EDENMAT)
!
!     ==========================================================================
!     == DECOMPOSE INTO SPIN-UP AND SPIN DOWN FOR NSPIN=2                     ==
!     ==========================================================================
      IF(NSPIN.EQ.2) THEN
        DO IR=1,NRL
          SVAR1=RHO(IR,1)
          SVAR2=RHO(IR,2)
          RHO(IR,1)=SVAR1+SVAR2
          RHO(IR,2)=SVAR1-SVAR2
        ENDDO
        DO IAT=1,NAT
          ISP=MAP%ISP(IAT)
          LMNX=MAP%LMNX(ISP)
          DO LMN1=1,LMNX
            DO LMN2=1,LMNX
              CSVAR1=DH(LMN1,LMN2,1,IAT)
              CSVAR2=DH(LMN1,LMN2,2,IAT)
              DH(LMN1,LMN2,1,IAT)=CSVAR1+CSVAR2
              DH(LMN1,LMN2,2,IAT)=CSVAR1-CSVAR2
              CSVAR1=DO(LMN1,LMN2,1,IAT)
              CSVAR2=DO(LMN1,LMN2,2,IAT)
              DO(LMN1,LMN2,1,IAT)=REAL(CSVAR1+CSVAR2,KIND=8)
              DO(LMN1,LMN2,2,IAT)=REAL(CSVAR1-CSVAR2,KIND=8)
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  SEND DATA TO BANDDATA MODULE                                        ==
!     ==  THIS HAS TO BE DONE HERE, BECAUSE VOFR, DH AND DO ARE ONLY          ==
!     ==  AVAILABLE HERE                                                      ==
!     ==========================================================================
      CALL BANDDATA$SETI4('LMNXX',LMNXX)
      CALL BANDDATA$SETI4('NAT',NAT)
      CALL BANDDATA$SETI4('NRL',NRL)
      CALL BANDDATA$SETI4('NDIMD',NDIMD)
      CALL BANDDATA$SETR8A('VOFRL',NRL*NDIMD,RHO(:,:))
      CALL BANDDATA$SETC8A('DH',LMNXX*LMNXX*NDIMD*NAT,DH(:,:,:,:))
      CALL BANDDATA$SETR8A('DO',LMNXX*LMNXX*NDIMD*NAT,DO(:,:,:,:))
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_BANDDATA' &
     &    ,INT(NRL,KIND=8),INT(LMNXX,KIND=8),INT(NDIMD,KIND=8) &
     &    ,INT(NAT,KIND=8),0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     ==  RECEIVE POTENTIALS FROM NTBO INTERFACE                              ==
!     ==========================================================================
                               CALL TIMING$CLOCKON('WAVES$FROMNTBO')
!      CALL WAVES$FROMNTBO()
      CALL WAVES$OFFSITEHAMIL() ! CALCULATE THIS$HPROJ
                               CALL TIMING$CLOCKOFF('WAVES$FROMNTBO')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_OFFSITEHAMIL' &
     &    ,INT(LMNXX,KIND=8),INT(NDIMD,KIND=8),INT(NAT,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     == FORCES AND STRESSES                                                  ==
!     ==========================================================================
      FORCE=0.D0
      IF(TFORCE.OR.TSTRESS) THEN
        STRESS1(:,:)=0.D0
        CALL WAVES$FORCE(NAT,LMNXX,NDIMD,DH,FORCE,STRESS1)
        STRESS=STRESS+STRESS1
!WRITE(*,FMT='("PRO STRESS ",3F15.7)')STRESS1(1,:)
!WRITE(*,FMT='("PRO STRESS ",3F15.7)')STRESS1(2,:)
!WRITE(*,FMT='("PRO STRESS ",3F15.7)')STRESS1(3,:)
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_FORCE' &
     &    ,INT(NAT,KIND=8),INT(LMNXX,KIND=8),INT(NDIMD,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     ==  CALL GONJUGATE GRADIENT                                             ==
!     ==========================================================================
!!$      IF(OPTIMIZERTYPE.EQ.'CG') THEN                            !KAESTNERCG
!!$         CALL TIMING$CLOCKON('CG')                              !KAESTNERCG
!!$         CALL CG$STATE_BY_STATE(MAP%NRL,NDIMD,RHO(:,:),CONVPSI,NAT,LMNXX,DH) !KAESTNERCG
!!$         CALL TIMING$CLOCKOFF('CG')                             !KAESTNERCG
!!$         TCONV=.FALSE. ! TCONV HAS NOT BEEN SET!!!
!!$         IF(TCONV) CALL STOPIT$SETL4('STOP',.TRUE.)             !KAESTNERCG
!!$      END IF                                                    !KAESTNERCG
!
!     ==========================================================================
!     ==  EVALUATE H*PSI                                                      ==
!     ==========================================================================
!PRINT*,'RHO',(SUM(ABS(RHO)).GT.0.D0.OR.SUM(ABS(RHO)).LE.0.D0)
      CALL WAVES$HPSI(NRL,NDIMD,NAT,LMNXX,RHO,DH)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_HPSI' &
     &    ,INT(NRL,KIND=8),INT(NDIMD,KIND=8),INT(NAT,KIND=8) &
     &    ,INT(LMNXX,KIND=8),0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
      DEALLOCATE(RHO)
      DEALLOCATE(DH)
      DEALLOCATE(DO)
!
!     ==========================================================================
!     ==  EVALUATE ENERGY EXPECTATION VALUES                                  ==
!     ==========================================================================
CALL TIMING$CLOCKON('W:EXPECT')
      ALLOCATE(EIG(NBX,NKPTL,NSPIN))
      ALLOCATE(HAMILTON(2,2))
!PRINT*,'==================BAND MISMATCHES========================='
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          NB=THIS%NB
          NBH=THIS%NBH
          NGL=GSET%NGL
          IF(.NOT.ASSOCIATED(THIS%EXPECTVAL))ALLOCATE(THIS%EXPECTVAL(NB))
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          THPSIRESIDENT=CPPAW_CUBLAS_ACC_HPSI_RESIDENCY_ENABLED() &
     &                   .AND.THIS%HPSI_ACC_RESIDENT
          IF(THPSIRESIDENT) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &          ('ACC_PRESENT_EXPECT_PSI0','ACC_COPY_EXPECT_PSI0_IN' &
     &          ,NGL,NDIM,NBH,THIS%PSI0)
            CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &          ('ACC_PRESENT_EXPECT_HPSI','ACC_COPY_EXPECT_HPSI_IN' &
     &          ,NGL,NDIM,NBH,THIS%HPSI)
#ENDIF
            CALL CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_RESIDENT(.TRUE.)
            CALL CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_FORCE(.TRUE.)
          END IF
!$ACC DATA PRESENT_OR_COPYIN(THIS%PSI0(1:NGL,1:NDIM,1:NBH) &
!$ACC&                       ,THIS%HPSI(1:NGL,1:NDIM,1:NBH)) IF(THPSIRESIDENT)
#ENDIF
          DO IB=1,NBH
            IF(GSET%TINV) THEN
              CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,1,2 &
     &              ,THIS%PSI0(:,:,IB),THIS%HPSI(:,:,IB),HAMILTON &
     &              ,'EXPECT_INV')
              EIG(2*IB-1,IKPT,ISPIN)=REAL(HAMILTON(1,1),KIND=8)
              EIG(2*IB  ,IKPT,ISPIN)=REAL(HAMILTON(2,2),KIND=8)
            ELSE 
              CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,1,1 &
     &            ,THIS%PSI0(:,:,IB),THIS%HPSI(:,:,IB),HAMILTON &
     &            ,'EXPECT')
              EIG(IB,IKPT,ISPIN)=REAL(HAMILTON(1,1),KIND=8)
            END IF
          ENDDO
#IF DEFINED(CPPVAR_CUBLAS_ACC)
!$ACC END DATA
          IF(THPSIRESIDENT) THEN
            CALL CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_FORCE(.FALSE.)
            CALL CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_RESIDENT(.FALSE.)
            IF(.NOT.THAMILTON) THEN
!$ACC EXIT DATA COPYOUT(THIS%HPSI(1:NGL,1:NDIM,1:NBH))
              THIS%HPSI_ACC_RESIDENT=.FALSE.
            END IF
          END IF
#ENDIF
          THIS%EXPECTVAL(:)=EIG(1:NB,IKPT,ISPIN)
CALL TIMESTEP$GETI4('ISTEP',ISVAR)
!!$DO IB=1,NBH-1
!!$  IF(EIG(IB+1,IKPT,ISPIN).LT.EIG(IB,IKPT,ISPIN)) THEN
!!$    WRITE(*,FMT='("BAND MISMATCH:",I10,3I5,2F10.5)')ISVAR,IB,IKPT,ISPIN &
!!$&                                              ,EIG(IB:IB+1,IKPT,ISPIN)
!!$  END IF
!!$ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(HAMILTON)
!     == OCCUPATIONS ===========================================================
      CALL WAVES_DEDFFROMNTBO(NB,NKPTL,NSPIN,EIG) ! NO FUNCTION 
      CALL WAVES_DYNOCCSETR8A('EPSILON',NB*NKPTL*NSPIN,EIG)
      DEALLOCATE(EIG)
!
!     ==========================================================================
!     ==  EVALUATE <S^2>                                                      ==
!     ==========================================================================
      IF ((NDIM.EQ.2).AND.THAMILTON) THEN   ! THAMILTON SHOULD BE REPLACED
        CALL TIMING$CLOCKON('S2')
        ALLOCATE(OCC(NBX,NKPTL,NSPIN))
        CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
        CALL WAVES_TOTALSPINRESET()
        DO IKPT=1,NKPTL
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          NB=THIS%NB
          NBH=THIS%NBH
          NGL=GSET%NGL         
          ALLOCATE(QMAT(2*NB*NSPIN,2*NB*NSPIN))
          CALL WAVES_SPINOROVERLAP(NBH,NB,IKPT,QMAT)
          CALL WAVES_TOTALSPIN(NBX,NB,NKPT,IKPT,NSPIN,OCC,QMAT)
          DEALLOCATE(QMAT)
        END DO
        DEALLOCATE(OCC)
        CALL TIMING$CLOCKOFF('S2')
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        CALL WAVES$REPORTSPIN(NFILO)
      END IF
!
!     ==========================================================================
!     ==  EVALUATE <PSI|H|PSI>                                                ==
!     ==========================================================================
      IF(THAMILTON) THEN
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            CALL PLANEWAVE$SELECT(GSET%ID)      
            NB=THIS%NB
            NBH=THIS%NBH
            NGL=GSET%NGL
            ALLOCATE(HAMILTON(NB,NB))
#IF DEFINED(CPPVAR_CUBLAS_ACC)
            THPSIRESIDENT=CPPAW_CUBLAS_ACC_HPSI_RESIDENCY_ENABLED() &
       &                   .AND.THIS%HPSI_ACC_RESIDENT
            IF(THPSIRESIDENT) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
              CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
       &          ('ACC_PRESENT_HAMILTON_PSI0','ACC_COPY_HAMILTON_PSI0_IN' &
       &          ,NGL,NDIM,NBH,THIS%PSI0)
              CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
       &          ('ACC_PRESENT_HAMILTON_HPSI','ACC_COPY_HAMILTON_HPSI_IN' &
       &          ,NGL,NDIM,NBH,THIS%HPSI)
#ENDIF
              CALL CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_RESIDENT(.TRUE.)
            END IF
!$ACC DATA PRESENT_OR_COPYIN(THIS%PSI0(1:NGL,1:NDIM,1:NBH) &
!$ACC&                       ,THIS%HPSI(1:NGL,1:NDIM,1:NBH)) IF(THPSIRESIDENT)
#ENDIF
            CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,NBH,NB,THIS%PSI0,THIS%HPSI &
     &                        ,HAMILTON,'HAMILTON')
#IF DEFINED(CPPVAR_CUBLAS_ACC)
!$ACC END DATA
            IF(THPSIRESIDENT) THEN
              CALL CPPAW_CUBLAS_ACC_SET_WAVE_OVERLAP_RESIDENT(.FALSE.)
!$ACC EXIT DATA COPYOUT(THIS%HPSI(1:NGL,1:NDIM,1:NBH))
              THIS%HPSI_ACC_RESIDENT=.FALSE.
            END IF
#ENDIF
            IF(.NOT.ASSOCIATED(THIS%EIGVAL))ALLOCATE(THIS%EIGVAL(NB))
            IF(.NOT.ASSOCIATED(THIS%EIGVEC))ALLOCATE(THIS%EIGVEC(NB,NB))
            CALL LIB$DIAGC8(NB,HAMILTON,THIS%EIGVAL,THIS%EIGVEC)
            DEALLOCATE(HAMILTON)
!
!           ====================================================================
!           == TRANSFORMS STATES ONTO HAMILTON EIGENSTATES. THIS ACCELERATES  ==
!           == CONERGENCE FOR SAFEORTHO=F. IT DOES NOT AFFECT THE OCCUPATIONS.==
!           == THEREFORE IT IS RECOMMENDED TO RESTART THOSE AS WELL.          ==
!           ==                                                                ==
!           == A SECOND CALL GIVES VELOCITY TO THE WAVE FUNCTIONS, WHICH IS   ==
!           == NOT UNDERSTOOD. THEREFORE IT IS EXECUTED ONLY ONCE AND         ==
!           == SWITCHED OFF BELOW.                                            ==
!           ====================================================================
            IF(TSTRAIGHTEN) THEN
PRINT*,'STRAIGHTENING STATES...'
              CALL WAVES_STRAIGHTEN(NDIM,NBH,NB,NGL,MAP%NPRO,THIS%EIGVEC &
    &                              ,THIS%PSI0,THIS%PSIM,THIS%HPSI,THIS%PROJ &
    &                              ,THIS%RLAM0)
!             == DISCARD DATA THAT IS INCONSISTENT WITH NEW STATES =============
              IF(ASSOCIATED(THIS%RLAMM))DEALLOCATE(THIS%RLAMM)
              IF(ASSOCIATED(THIS%RLAM2M))DEALLOCATE(THIS%RLAM2M)
              IF(ASSOCIATED(THIS%RLAM3M))DEALLOCATE(THIS%RLAM3M)
PRINT*,'......STATES STRAIGHTENED'
            END IF !TSTRAIGHTEN
          ENDDO
        ENDDO
!
!       == ENSURE THAT STRAIGTHENING IS DONE ONLY ONCE... ======================
        IF(TSTRAIGHTEN) TSTRAIGHTEN=.FALSE. ! REMOVED OCT.23, 2024
      ELSE
        IF(OPTIMIZERTYPE.NE.'CG') THEN    !KAESTNERCG
          IF(ASSOCIATED(THIS%EIGVAL))DEALLOCATE(THIS%EIGVAL)
          IF(ASSOCIATED(THIS%EIGVEC))DEALLOCATE(THIS%EIGVEC)
        END IF                         !KAESTNERCG
      END IF
!
CALL TIMING$CLOCKOFF('W:EXPECT')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_EXPECT' &
     &    ,INT(NBX,KIND=8),INT(NKPTL,KIND=8),INT(NSPIN,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     ==  FIRST HALF OF WAVE FUNCTIONS KINETIC ENERGY                         ==
!     ==========================================================================
      CALL WAVES_WAVEKINETIC(WAVEEKIN1)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_WAVEKINETIC' &
     &    ,0_8,0_8,0_8,0_8,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T0)
#ENDIF
!
!     ==========================================================================
!     ==  CLOSE DOWN                                                          ==
!     ==========================================================================
!     == FORCES ================================================================
      ALLOCATE(FORCET(3,NAT))
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCET)
      FORCE=FORCET+FORCE
      CALL ATOMLIST$SETR8A('FORCE',0,3*NAT,FORCE)
      DEALLOCATE(FORCET)
!     == STRESS ================================================================
      CALL CELL$GETR8A('STRESS_I',9,STRESS1)
      STRESS=STRESS1-STRESS  ! IN THIS ROUTINE STRESS=+DE/DEPSILON!
      CALL CELL$SETR8A('STRESS_I',9,STRESS)
!     == DEALLOCATE ARRAYS =====================================================
      DEALLOCATE(FORCE)
      DEALLOCATE(R)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_ETOT_T1)
      CALL ACCELPROFILE$ADD('PAW_ETOT_CLOSE' &
     &    ,0_8,0_8,0_8,0_8,0.D0,0.D0,ACCEL_ETOT_T1-ACCEL_ETOT_T0)
#ENDIF
                              CALL TIMING$CLOCKOFF('WAVES$ETOT')
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_STRAIGHTEN(NDIM,NBH,NB,NGL,NPRO,UMAT &
    &                             ,PSI0,PSIM,HPSI,PROJ,RLAM0)
!     **************************************************************************
!     ** PERFORMS A UNITARY TRANSFORM OF THE ACTIAL STATE OF WAVE FUNCTIONS
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NDIM         ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN)    :: NBH          ! #(SUPER WAVE FUNCTIONS)
      INTEGER(4),INTENT(IN)    :: NB           ! #(BANDS)
      INTEGER(4),INTENT(IN)    :: NGL          ! #(LOCAL G-VECTORS)
      INTEGER(4),INTENT(IN)    :: NPRO         ! #(PROJECTIONS)
      COMPLEX(8),INTENT(INOUT) :: UMAT(NB,NB)  ! TRANSFORMATION
      COMPLEX(8),INTENT(INOUT) :: PSI0(NGL,NDIM,NBH)
      COMPLEX(8),INTENT(INOUT) :: PSIM(NGL,NDIM,NBH)
      COMPLEX(8),INTENT(INOUT) :: HPSI(NGL,NDIM,NBH)
      COMPLEX(8),INTENT(INOUT) :: PROJ(NDIM,NBH,NPRO)
      COMPLEX(8),INTENT(INOUT) :: RLAM0(NB,NB)
      COMPLEX(8),ALLOCATABLE   :: AUX(:,:,:)
      INTEGER(4)               :: I
!     **************************************************************************
!
!     ==========================================================================
!     ==  TRANSFORM WAVE FUNCTIONS                                            ==
!     ==========================================================================
      ALLOCATE(AUX(NGL,NDIM,NBH))
      AUX=PSI0
      PSI0=(0.D0,0.D0)
      CALL WAVES_ADDPSI(NGL,NDIM,NBH,NB,PSI0,AUX,UMAT)
      AUX=PSIM
      PSIM=(0.D0,0.D0)
      CALL WAVES_ADDPSI(NGL,NDIM,NBH,NB,PSIM,AUX,UMAT)
      AUX=HPSI
      HPSI=(0.D0,0.D0)
      CALL WAVES_ADDPSI(NGL,NDIM,NBH,NB,HPSI,AUX,UMAT)
      DEALLOCATE(AUX)
!
!     ==========================================================================
!     ==  TRANSFORM PROJECTIONS                                               ==
!     ==========================================================================
      ALLOCATE(AUX(NDIM,NBH,NPRO))
      AUX=PROJ
      PROJ=(0.D0,0.D0)
      CALL WAVES_ADDOPROJ(NPRO,NDIM,NBH,NB,PROJ,AUX,UMAT)
      DEALLOCATE(AUX)
!
!     ==========================================================================
!     ==  TRANSFORM LAGRANGE MULTIPLIERS                                      ==
!     ==========================================================================
      RLAM0=MATMUL(TRANSPOSE(CONJG(UMAT)),MATMUL(RLAM0,UMAT))
!
!     ==========================================================================
!     == FINALLY, TRANSFORM EIGENSTATES TO A UNIT MATRIX.                     ==
!     ==========================================================================
      UMAT(:,:)=(0.D0,0.D0)
      DO I=1,NB
        UMAT(I,I)=(1.D0,0.D0)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_VOFRHO(NRL,NDIMD,RHO,RHOB,NAT,LMRXX,QLM,VQLM)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
!      USE WAVES_MODULE, ONLY : WAVES_SELECTWV,GSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NDIMD  !#(DENSITY COMPONENTS)
      INTEGER(4),INTENT(IN)   :: NRL    !#(LOCAL R-SPACE GRID POINTS)
      REAL(8)   ,INTENT(INOUT):: RHO(NRL,NDIMD)  
      REAL(8)   ,INTENT(OUT)  :: RHOB   ! BACKGROUND DENSITY
      INTEGER(4),INTENT(IN)   :: NAT    !#(ATOMS)
      INTEGER(4),INTENT(IN)   :: LMRXX   !#(ANGULAR MOMENTA FOR 1-C DENSITY)
      REAL(8)   ,INTENT(IN)   :: QLM(LMRXX,NAT)  ! MULTIPOLE MOMENTS
      REAL(8)   ,INTENT(OUT)  :: VQLM(LMRXX,NAT) ! "MULTIPOLE" POTENTIALS
      REAL(8)   ,ALLOCATABLE  :: RHO_V(:,:)     ! CHARGE DENSITY
      REAL(8)                 :: RBAS(3,3)      ! LATTICE VECTORS
      REAL(8)                 :: R(3,NAT)       ! ATOMIC POSITIONS
      REAL(8)                 :: STRESS(3,3)    ! STRESS TENSOR
      REAL(8)                 :: STRESST(3,3)    ! STRESS TENSOR
      REAL(8)                 :: FORCE(3,NAT)
      REAL(8)                 :: FORCET(3,NAT)
      INTEGER(4)              :: NR1L,NR1L_V,NR2,NR3
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES_VOFRHO')

!     == ASSUMES THAT NR2 AND NR3 ARE IDENTICAL FOR DENSITY AND FOR
!     == ALL WAVES IN THE GROUP/ NR1L ON THE OTHER HAND MAY DIFFER
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$GETI4('NR1L',NR1L_V)
      CALL PLANEWAVE$GETI4('NR2',NR2)
      CALL PLANEWAVE$GETI4('NR3',NR3)
      NR1L=NRL/(NR2*NR3)  ! NRL=NR1L*NR2*NR3  REFERS TO THE WAVE FUNCTION!
!
!     ==========================================================================
!     == PERFORM SUM OVER K-POINT GROUPS AND REDISTRIBUTE                     ==
!     ==========================================================================
      ALLOCATE(RHO_V(NR1L_V*NR2*NR3,NDIMD))
      CALL WAVES_MAPPSITOPOT('PSITOPOT',NR1L,NR1L_V,NR2,NR3,NDIMD,RHO,RHO_V)
!
!     ==========================================================================
!     == CALCULATE POTENTIAL FROM THE CHARGE DENSITY                          ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
      FORCE(:,:)=0.D0
      STRESS(:,:)=0.D0
      VQLM(:,:)=0.D0
      CALL POTENTIAL$VOFRHO(NR1L_V*NR2*NR3,NDIMD,RHO_V,LMRXX,NAT,QLM,VQLM &
     &                     ,R,FORCE,RBAS,STRESS,RHOB)
!
!     ==========================================================================
!     == MAP POTENTIAL BACK ONTO THE WAVE FUNCTION GRID                       ==
!     ==========================================================================
      CALL WAVES_MAPPSITOPOT('POTTOPSI',NR1L,NR1L_V,NR2,NR3,NDIMD,RHO,RHO_V)
      DEALLOCATE(RHO_V)
!
!     ==========================================================================
!     ==  STORE FORCES AND STRESSES BACK TO THE OWNING OBJECTS                ==
!     ==========================================================================
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCET)
      FORCE=FORCET+FORCE
      CALL ATOMLIST$SETR8A('FORCE',0,3*NAT,FORCE)
!
      CALL CELL$GETR8A('STRESS_I',9,STRESST)
      STRESS=STRESST-STRESS  ! IN THIS ROUTINE STRESS=+DE/DEPSILON!
      CALL CELL$SETR8A('STRESS_I',9,STRESS)
!
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES_VOFRHO
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_MAPPSITOPOT(ID,NR1L_P,NR1L_V,NR2,NR3,NDIMD,RHO_P,RHO_V)
!     **************************************************************************
!     **                                                                      **
!     **  ADDS UP THE DENSITY OVER K-POINTS AND MAPS THE DENSITY              **
!     **  FROM THE SHEETS OF THE WAVE FUNCTIONS INTO THOSE OF THE             **
!     **  POTENTIAL AND VICE VERSA.                                           **
!     **                                                                      **
!     **  THE DENSITY IS DIVIDED INTO SHEETS, WHERE EACH SHEET IS             **
!     **  KEPT BY ONE NODE. IN THE K-POINT PARALLELIZATION, THE               **
!     **  SHEETS DIFFER FOR THE WAVE FUNCTIONS AND THE DENSITY,               **
!     **  BECAUSE A K-GROUP MUST KEEP ALL GRID POINTS OF THE DENSITY          **
!     **  WHILE FOR THE POTENTIAL, THEY ARE DISTRIBUTED OVER ALL              **
!     **  NODES OF THE MONOMER GROUP.                                         **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **    THE DIVISION IN SHEETS IS THE SAME FOR ALL KPOINTS WITHIN         **
!     **    ONE K-GROUP                                                       **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE WAVES_MODULE, ONLY : WAVES_SELECTWV,GSET,NKPTL
      USE MPE_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: NDIMD
      INTEGER(4)  ,INTENT(IN) :: NR1L_P
      INTEGER(4)  ,INTENT(IN) :: NR1L_V
      INTEGER(4)  ,INTENT(IN) :: NR2
      INTEGER(4)  ,INTENT(IN) :: NR3
      REAL(8)     ,INTENT(INOUT) :: RHO_P(NR1L_P,NR2,NR3,NDIMD)
      REAL(8)     ,INTENT(INOUT) :: RHO_V(NR1L_V,NR2,NR3,NDIMD)
      REAL(8)     ,ALLOCATABLE   :: ARR(:,:,:,:)
      INTEGER(4)  ,ALLOCATABLE   :: NR1FIRST_P(:)
      INTEGER(4)  ,ALLOCATABLE   :: NR1LAST_P(:)
      INTEGER(4)  ,ALLOCATABLE   :: NR1FIRST_V(:)
      INTEGER(4)  ,ALLOCATABLE   :: NR1LAST_V(:)
      INTEGER(4)  ,ALLOCATABLE   :: NR1LARR(:)
      INTEGER(4)            :: NTASKS,THISTASK
      INTEGER(4)            :: ISVAR
      INTEGER(4)            :: NR1START_P
      INTEGER(4)            :: ITASK,ITASK_P,ITASK_V
      INTEGER(4)            :: IWORK2(2)
      INTEGER(4)            :: IRSTART,IREND
      INTEGER(4)            :: IR1_P,IR2_P,IR1_V,IR2_V
INTEGER(4) :: IR1,IR2,IR3
REAL(8)    :: SVAR
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES_MAPPSITOPOT')
      IF(ID.NE.'PSITOPOT'.AND.ID.NE.'POTTOPSI') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES_MAPPSITOPOT')
      ENDIF
!
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)      
!
!     ==========================================================================
!     == DETERMINE WHICH PLANES OF RHO_V ARE HELD BY EACH PROCESSOR           ==
!     ==========================================================================
      CALL PLANEWAVE$SELECT('DENSITY')
!     __THE GRID POINTS OF DENSITY AND POTENTIAL ARE DIVIDED UP INTO SHEETS.
!     __EACH SHEAT HAS A RANGE OF NR1L INDICES OF THE FIRST INDEX.
!     __THE LOCAL VALUES OF NR1L ARE STORED ON THE ARRAY NR1LARR
!     __WHICH SATISFIES SUM(NR1LARR(1:NTASKS))=NR1.
!     __THE FIRST INDEX LOCAL INDEX IR ON TASK 'ITASK' IS 
!     __     IR_FIRST=1+SUM(NR1LARR(1:ITASK-1) AND THE LAST VALUE IS
!     __     IR_LAST=IR_FIRST-1+NR1LARR(ITASK)
!PRINT*,'MAPPSITOPOT TASK ',THISTASK,' ID(1)=',ID
      ALLOCATE(NR1LARR(NTASKS))
      NR1LARR(:)=0 !JUST INITIALIZATION
      CALL PLANEWAVE$GETI4A('NR1LARR',NTASKS,NR1LARR)
!PRINT*,'MAPPSITOPOT TASK ',THISTASK,' NR1LARR=',NR1LARR
!
      ALLOCATE(NR1FIRST_V(NTASKS))
      ALLOCATE(NR1LAST_V(NTASKS))
      ISVAR=0
      DO ITASK=1,NTASKS
        NR1FIRST_V(ITASK)=ISVAR+1
        ISVAR=ISVAR+NR1LARR(ITASK)
        NR1LAST_V(ITASK)=ISVAR
      ENDDO
!
      DEALLOCATE(NR1LARR)
!
!     ==========================================================================
!     == COMMUNICATE WHICH PLANES OF RHO_P ARE HELD BY EACH PROCESSOR         ==
!     ==========================================================================
!     __ IT IS POSSIBLE THAT THERE IS NO K-POINT ON A TASK
      IF(NKPTL.GT.0) THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR1START',NR1START_P)
      ELSE
        NR1START_P=1
      END IF
!PRINT*,'MAPPSITOPOT TASK ',THISTASK,' NR1START_P=',NR1START_P,NKPTL
      ALLOCATE(NR1FIRST_P(NTASKS))
      ALLOCATE(NR1LAST_P(NTASKS))
      DO ITASK=1,NTASKS
        IWORK2(1)=NR1START_P
        IWORK2(2)=NR1START_P-1+NR1L_P
        CALL MPE$BROADCAST('MONOMER',ITASK,IWORK2(:))
        NR1FIRST_P(ITASK)=IWORK2(1)
        NR1LAST_P(ITASK)=IWORK2(2)
      ENDDO
!PRINT*,'MAPPSITOPOT TASK ',THISTASK,' ID(2)=',ID
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      IF(ID.EQ.'PSITOPOT') THEN
        ALLOCATE(ARR(NR1L_V,NR2,NR3,NDIMD))
        RHO_V(:,:,:,:)=0.D0
      END IF
!PRINT*,'MAPPSITOPOT TASK ',THISTASK,'BEFORE DISTRIBUTE ',ID
!
!     ==========================================================================
!     ==  DISTRIBUTE DENSITY/POTENTIAL                                        ==
!     ==========================================================================
      DO ITASK_P=1,NTASKS
        DO ITASK_V=1,NTASKS
          IF(.NOT.(ITASK_P.EQ.THISTASK.OR.ITASK_V.EQ.THISTASK)) CYCLE
!
          IRSTART=MAX(NR1FIRST_V(ITASK_V),NR1FIRST_P(ITASK_P)) 
          IREND  =MIN(NR1LAST_V(ITASK_V),NR1LAST_P(ITASK_P)) 
          IF(IRSTART.GT.IREND) CYCLE  ! NO MESSAGE, GO TO NEXT....
!
!         ======================================================================
!         ==                                                                  ==
!         ======================================================================
          IR1_P=IRSTART-NR1FIRST_P(ITASK_P)+1
          IR2_P=IREND-NR1FIRST_P(ITASK_P)+1
          IR1_V=IRSTART-NR1FIRST_V(ITASK_V)+1
          IR2_V=IREND-NR1FIRST_V(ITASK_V)+1
          IF(ITASK_P.EQ.ITASK_V) THEN
!           == BOTH ON THE SAME TASK; JUST COPY
            IF(ID.EQ.'PSITOPOT') THEN
              RHO_V(IR1_V:IR2_V,:,:,:)=RHO_V(IR1_V:IR2_V,:,:,:) &
     &                                +RHO_P(IR1_P:IR2_P,:,:,:)
            ELSE
              RHO_P(IR1_P:IR2_P,:,:,:)=RHO_V(IR1_V:IR2_V,:,:,:)
            END IF
          ELSE
            IF(ID.EQ.'PSITOPOT') THEN
              IF(THISTASK.EQ.ITASK_P) THEN
!               __ SEND ____
!PRINT*,'SEND   : ',THISTASK,':',ITASK_P,'->',ITASK_V,' ; ',IR1_P,IR2_P
                CALL MPE$SENDRECEIVE('MONOMER',ITASK_P,ITASK_V &
     &                              ,RHO_P(IR1_P:IR2_P,:,:,:))
              ELSE
!               __ RECEIVE ____
!PRINT*,'RECEIVE: ',THISTASK,':',ITASK_P,'->',ITASK_V,' ; ',IR1_V,IR2_V
                CALL MPE$SENDRECEIVE('MONOMER',ITASK_P,ITASK_V &
     &                              ,ARR(1:IR2_P-IR1_P+1,:,:,:))
                RHO_V(IR1_V:IR2_V,:,:,:)=RHO_V(IR1_V:IR2_V,:,:,:) &
     &                                +ARR(1:IR2_P-IR1_P+1,:,:,:)
              END IF
            ELSE
              IF(THISTASK.EQ.ITASK_V) THEN
!               ____SEND____
!PRINT*,'SEND   : ',THISTASK,':',ITASK_V,'->',ITASK_P
                CALL MPE$SENDRECEIVE('MONOMER',ITASK_V,ITASK_P &
      &                             ,RHO_V(IR1_V:IR2_V,:,:,:))
              ELSE
!               __ RECEIVE______
!PRINT*,'RECEIVE: ', THISTASK,':',ITASK_V,'->',ITASK_P
                CALL MPE$SENDRECEIVE('MONOMER',ITASK_V,ITASK_P,RHO_P(IR1_P:IR2_P,:,:,:))
              END IF
            END IF
          END IF
        ENDDO
      ENDDO      
!PRINT*,'MAPPSITOPOT TASK ',THISTASK,'AFTER DISTRIBUTE'


!PRINT*,'MAPPSITOPOT: ID=',TRIM(ID),' TASK=',THISTASK,': DONE'
      IF(ALLOCATED(ARR))DEALLOCATE(ARR)

IF(ID.EQ.'PSITOPOT_X') THEN
! PRINT*,THISTASK,'NR1FIRST_V',NR1FIRST_V  
! PRINT*,THISTASK,'NR1LAST_V ',NR1LAST_V  
! PRINT*,THISTASK,'NR1FIRST_P',NR1FIRST_P  
! PRINT*,THISTASK,'NR1LAST_P ',NR1LAST_P  
! PRINT*,'RHO_V(X,1,1) ',RHO_V(1:10,1,1,1)
 SVAR=0.D0
 DO IR3=1,NR3
   DO IR2=1,NR2
     DO IR1=1,NR1L_V
        SVAR=SVAR+ABS(RHO_V(IR1,IR2,IR2,1))
     ENDDO
   ENDDO
 ENDDO
 CALL MPE$COMBINE('MONOMER','+',SVAR)
 IR1=NR1L_V*NR2*NR3
 CALL MPE$COMBINE('MONOMER','+',IR1)
! PRINT*,THISTASK,'SUM ',SVAR/REAL(IR1,KIND=8)
END IF
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES_MAPPSITOPOT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$KAESTNERCG1(TFIRST_)
!     **************************************************************************
!     **                                                                      **
!     **  THIS IS FOR KAESTNERS CONJUGATE GRADIENT IMPLEMENTATION             **
!     **  HANDLE WITH SPECIAL CARE TO AVOID PROBLEMS WITH ALLOCATION          **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TFIRST_
      INTEGER(4)            :: IKPT,ISPIN
      INTEGER(4)            :: NB
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$KAESTNERCG1')
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NB=THIS%NB
          IF(.NOT.ASSOCIATED(THIS%EIGVAL))ALLOCATE(THIS%EIGVAL(NB))
          IF(.NOT.ASSOCIATED(THIS%EIGVEC))ALLOCATE(THIS%EIGVEC(NB,NB))
!         KAESTNER PROBABLY NEEDS THE VALUES OF EIGVAL AND EIGVEC FROM THE 
!         PREVIOUS ITERATIONS. HOWVEER THE ARRAYS HAVE BEEN DEALLOCATED
!         BEFORE THEY ARE ALLOCATED AGAIN HERE, SO THAT THIS INFORMATION
!         IS LOST. FIND THE DEALLOCATION IN WAVES$ETOT.
          CALL ERROR$MSG('CORRECT CODING ERROR BEFORE USING CONJUGATE GRADIENT OPTION')
          CALL ERROR$MSG('LOOK INTO THIS ROUTINE FOR MORE INFORMATION')
          CALL ERROR$STOP('WAVES$KAESTNERCG1')
          IF(TFIRST_) THEN
            ALLOCATE(THIS%EIGVAL(NB))          !KAESTNERCG
            THIS%EIGVAL(:)=1.D100              !KAESTNERCG
          END IF
        ENDDO
      ENDDO
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$KAESTNERCG1
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$RANDOMIZE()
!     **************************************************************************
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE WAVES_MODULE !, ONLY : NSPIN,NKPTL,NDIM,GSET,THIS,AMPRANDOM
      IMPLICIT NONE
      INTEGER(4)             :: IKPT,ISPIN
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NBH
      REAL(8)   ,ALLOCATABLE :: G2(:)
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$RANDOMIZE')
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          ALLOCATE(G2(NGL))
          CALL PLANEWAVE$GETR8A('G2',NGL,G2)
          CALL WAVES_RANDOMIZE(NGL,NDIM,NBH,AMPRANDOM,G2,THIS%PSIM)
          DEALLOCATE(G2)
        ENDDO
      ENDDO
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$RANDOMIZE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$GRAMMSCHMIDT()
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **  REMARK: I IS IMPORTANT TO DO THE ORTHOGONALIZATION WITH THE         **
!     **    CONSISTENT SET OF ATOMIC POSITIONS.                               **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE WAVES_MODULE, ONLY : NKPTL,NSPIN,NDIM,MAP,GSET,THIS &
     &                        ,WAVES_SELECTWV
      IMPLICIT NONE
      INTEGER(4)             :: IKPT,ISPIN
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NB
      INTEGER(4)             :: NBH
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      REAL(8)   ,ALLOCATABLE :: RM(:,:)
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$GRAMMSCHMIDT')
                              CALL TIMING$CLOCKON('W:GRAMSCHMIDT')
      NAT=MAP%NAT
      ALLOCATE(R0(3,NAT))
      ALLOCATE(RM(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      CALL ATOMLIST$GETR8A('R(-)',0,3*NAT,RM)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NB=THIS%NB
          NBH=THIS%NBH
          CALL WAVES_GRAMSCHMIDT(MAP,GSET,NAT,R0,NGL,NDIM,NBH,NB,THIS%PSI0 &
     &                          ,'PSI0')
          CALL WAVES_GRAMSCHMIDT(MAP,GSET,NAT,RM,NGL,NDIM,NBH,NB,THIS%PSIM &
     &                          ,'PSIM')
        ENDDO
      ENDDO
                              CALL TIMING$CLOCKOFF('W:GRAMSCHMIDT')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$GRAMMSCHMIDT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$PROJECTIONS(ID)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: ID
      INTEGER(4)             :: IKPT,ISPIN
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NBH
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R(:,:)
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$PROJECTIONS')
                              CALL TIMING$CLOCKON('W:PROJ')
      NAT=MAP%NAT
      ALLOCATE(R(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          CALL WAVES_ACC_THIS_PROJ_DELETE()
#ENDIF
          IF(ID.EQ.'PSI0') THEN
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,MAP%NPRO &
     &                            ,THIS%PSI0,THIS%PROJ,'SETUP0')
            CALL MPE$COMBINE('K','+',THIS%PROJ)
#IF DEFINED(CPPVAR_CUBLAS_ACC)
            CALL WAVES_ACC_THIS_PROJ_COPYIN()
#ENDIF
          ELSE
            CALL ERROR$MSG('ID NOT RECOGNIZED')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('WAVES$PROJECTIONS')
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(R)
                              CALL TIMING$CLOCKOFF('W:PROJ')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$PROJECTIONS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$EKIN(EKIN,STRESS)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES PSEUDO-KINETIC ENERGY                                     **
!     **  AND THE CORRESPONDING STRESS                                        **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: EKIN
      REAL(8)   ,INTENT(OUT) :: STRESS(3,3)
      INTEGER(4)             :: IKPT,ISPIN
      INTEGER(4)             :: NGL,NBH,NB,NBX
      REAL(8)                :: GWEIGHT
      REAL(8)                :: EKIN1
      REAL(8)                :: STRESS1(3,3)
      LOGICAL(4)             :: TSTRESS
!REDUCEKPOINTS UP TO HERE
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:) !(NBX,NKPTL,NSPIN)
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$EKIN')
                              CALL TIMING$CLOCKON('W:EKIN')
!
!     ==========================================================================
!     ==  GET OCCUPATIONS FROM DYNOCC OBJECT                                  ==
!     ==========================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      EKIN=0.D0
      STRESS(:,:)=0.D0
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETR8('GWEIGHT',GWEIGHT)
          NGL=GSET%NGL
          NBH=THIS%NBH
          NB=THIS%NB
          CALL WAVES_EKIN(NGL,NDIM,NBH,NB,OCC(1,IKPT,ISPIN),GWEIGHT &
      &                         ,THIS%PSI0,EKIN1,TSTRESS,STRESS1 &
      &                         ,TBUCKET,GSET%BUCKET,GSET%DBUCKET)
          EKIN  =EKIN  +EKIN1
          STRESS=STRESS+STRESS1
        ENDDO
      ENDDO
      DEALLOCATE(OCC)
      CALL MPE$COMBINE('MONOMER','+',EKIN)
      CALL MPE$COMBINE('MONOMER','+',STRESS)
                              CALL TIMING$CLOCKOFF('W:EKIN')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$EKIN
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_EKIN(NGL,NDIM,NBH,NB,F,GWEIGHT,PSI,EKIN &
     &                         ,TSTRESS,STRESS,TBUCKET,BUCKET,DBUCKET)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATE PS KINETIC ENERGY IN G-SPACE                               **
!     **  EVALUATE NUMBER OF ELECTRONS IN G-SPACE                             **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **    REQUIRES PLANEWAVE OBJECT TO BE SET PROPERLY                      **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1999)***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB         ! #(STATES)
      INTEGER(4),INTENT(IN) :: NBH        ! #(WAVE FUNCTIONS)
      INTEGER(4),INTENT(IN) :: NDIM       ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN) :: NGL        ! #(PLANE WAVES)
      REAL(8)   ,INTENT(IN) :: F(NB)      ! OCCUPATION
      REAL(8)   ,INTENT(IN) :: GWEIGHT
      COMPLEX(8),INTENT(IN) :: PSI(NGL,NDIM,NBH) ! PS-WAVE FUNCTION
      REAL(8)   ,INTENT(OUT):: EKIN       ! KINETIC ENERGY
      LOGICAL(4),INTENT(IN) :: TSTRESS    ! SWITCH STRESS ON/OFF
      REAL(8)   ,INTENT(OUT):: STRESS(3,3)! STRESS
      LOGICAL(4),INTENT(IN) :: TBUCKET    ! BUCKET POTENTIAL PRESENT
      REAL(8)   ,INTENT(IN) :: BUCKET(NGL) ! BUCKET POTENTIAL
      REAL(8)   ,INTENT(IN) :: DBUCKET(NGL) ! 1/G *DBUCKET/DG
      REAL(8)   ,ALLOCATABLE:: G2(:)      ! G**2
      REAL(8)   ,ALLOCATABLE:: GVEC(:,:)  ! G   
      COMPLEX(8),ALLOCATABLE:: PSI1(:,:)
      REAL(8)   ,ALLOCATABLE:: DMAT(:)
      LOGICAL(4)            :: TINV
      INTEGER(4)            :: IB,IG,IDIM
      REAL(8)               :: FP,FM
      REAL(8)   ,PARAMETER  :: DSMALL=1.D-12
      REAL(8)               :: EBUCKET
      REAL(8)               :: SVAR
!     **************************************************************************
!
!     ==========================================================================
!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT        ==
!     ==========================================================================
      CALL PLANEWAVE$GETL4('TINV',TINV)
      IF(TINV) THEN
        IF(NBH.NE.(NB+1)/2) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_EKIN')
        END IF
      ELSE 
        IF(NBH.NE.NB) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_EKIN')
        END IF
      END IF
!
!     ==========================================================================
!     ==  CALCULATE DMAT(G)=F(IB)*PSI^*(G,IB)*PSI(G,IB)                       ==
!     ==========================================================================
      ALLOCATE(DMAT(NGL))
      ALLOCATE(PSI1(NGL,NDIM))
      DMAT(:)=0.D0
      DO IB=1,NBH
!       == DETERMINE OCCUPATIONS ===============================================
        IF(TINV) THEN
          FP=0.5D0*(F(2*IB-1)+F(2*IB))
          FM=0.5D0*(F(2*IB-1)-F(2*IB))
        ELSE
          FP=F(IB)
          FM=0.D0
        END IF
!
!       ========================================================================
!       == GENERAL WAVE FUNCTION / FIRST PART OF SUPER WAVE FUNCTIONS         ==
!       == <PSI_+|G^2|PSI_+> FOR SUPER WAVE FUNCTIONS                         ==
!       ========================================================================
        IF(FP.EQ.0.D0.AND.FM.EQ.0.D0) CYCLE
        DO IDIM=1,NDIM
          DO IG=1,NGL
            DMAT(IG)=DMAT(IG) &
    &                 +FP*REAL(CONJG(PSI(IG,IDIM,IB))*PSI(IG,IDIM,IB),KIND=8)
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  <PSI_+|G^2|PSI_->                                                 ==
!       ========================================================================
        IF(.NOT.TINV.OR.ABS(FM).LT.DSMALL) CYCLE
        DO IDIM=1,NDIM
          CALL PLANEWAVE$INVERTG(NGL,PSI(1,IDIM,IB),PSI1(1,IDIM))
        ENDDO
        DO IDIM=1,NDIM
          DO IG=1,NGL
            DMAT(IG)=DMAT(IG) &
     &                +FM*REAL(CONJG(PSI(IG,IDIM,IB))*PSI1(IG,IDIM),KIND=8)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(PSI1)
!
!     ==========================================================================
!     ==  CALCULATE KINETIC ENERGY AND STRESS                                 ==
!     ==========================================================================
      IF(TSTRESS) THEN
        ALLOCATE(GVEC(3,NGL))
        CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
!       == KINETIC ENERGY AND KINETIC STRESS ===================================
        STRESS(:,:)=0.D0
        DO IG=1,NGL
          STRESS(1,1)=STRESS(1,1)+GVEC(1,IG)*GVEC(1,IG)*DMAT(IG)
          STRESS(1,2)=STRESS(1,2)+GVEC(1,IG)*GVEC(2,IG)*DMAT(IG)
          STRESS(1,3)=STRESS(1,3)+GVEC(1,IG)*GVEC(3,IG)*DMAT(IG)
          STRESS(2,2)=STRESS(2,2)+GVEC(2,IG)*GVEC(2,IG)*DMAT(IG)
          STRESS(2,3)=STRESS(2,3)+GVEC(2,IG)*GVEC(3,IG)*DMAT(IG)
          STRESS(3,3)=STRESS(3,3)+GVEC(3,IG)*GVEC(3,IG)*DMAT(IG)
        ENDDO
        EKIN=0.5D0*(STRESS(1,1)+STRESS(2,2)+STRESS(3,3))
        STRESS=-STRESS
!       == ENERGY AND STRESS DUE TO THE BUCKET POTENTIAL =======================
        IF(TBUCKET) THEN
          EBUCKET=0.D0
          DO IG=1,NGL
            IF(BUCKET(IG).EQ.0.D0) CYCLE
            EBUCKET=EBUCKET+BUCKET(IG)*DMAT(IG)
            SVAR=DMAT(IG)*DBUCKET(IG)
            STRESS(1,1)=STRESS(1,1)-GVEC(1,IG)*GVEC(1,IG)*SVAR
            STRESS(1,2)=STRESS(1,2)-GVEC(1,IG)*GVEC(2,IG)*SVAR
            STRESS(1,3)=STRESS(1,3)-GVEC(1,IG)*GVEC(3,IG)*SVAR
            STRESS(2,2)=STRESS(2,2)-GVEC(2,IG)*GVEC(2,IG)*SVAR
            STRESS(2,3)=STRESS(2,3)-GVEC(2,IG)*GVEC(3,IG)*SVAR
            STRESS(3,3)=STRESS(3,3)-GVEC(3,IG)*GVEC(3,IG)*SVAR
          ENDDO
          EKIN=EKIN+EBUCKET
        END IF
!       == PARALLELIZE AND CLOSE DOWN ==========================================
        STRESS(2,1)=STRESS(1,2)
        STRESS(3,1)=STRESS(1,3)
        STRESS(3,2)=STRESS(2,3)
        EKIN=EKIN*GWEIGHT
        STRESS=STRESS*GWEIGHT
        DEALLOCATE(GVEC)
        DEALLOCATE(DMAT)
      ELSE
        ALLOCATE(G2(NGL))
        CALL PLANEWAVE$GETR8A('G2',NGL,G2)
!       == KINETIC ENERGY ======================================================
        EKIN=0.D0
        DO IG=1,NGL
          EKIN=EKIN+G2(IG)*DMAT(IG)
        ENDDO
        EKIN=0.5D0*EKIN
!       == BUCKET POTENTIAL ====================================================
        IF(TBUCKET) THEN
          EBUCKET=0.D0
          DO IG=1,NGL
            EBUCKET=EBUCKET+BUCKET(IG)*DMAT(IG)
          ENDDO
          EKIN=EKIN+EBUCKET
        END IF
!       ==  WEIGHT RESULT AND CLOSE DOWN =======================================
        EKIN=EKIN*GWEIGHT
        DEALLOCATE(DMAT)
        DEALLOCATE(G2)
        STRESS(:,:)=0.D0
      END IF
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$DENMAT(LMNXX,NDIMD_,NAT,DENMAT,EDENMAT)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES ONE-CENTER DENSITY MATRIX FROM THE ACTUAL                 **
!     **  PROJECTIONS <PRO|PSI> IN THE WAVES OBJECT                           **
!     **                                                                      **
!     **  ONE-CENTER DENSITY MATRICES                                         **
!     **  SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                             **
!     **  SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)                      **
!     **  NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z)        **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE MPE_MODULE
      USE WAVES_MODULE
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &        CPPAW_CUBLAS_ACC_DENMAT_ENERGY_ENABLED
#ENDIF
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: LMNXX
      INTEGER(4),INTENT(IN)  :: NDIMD_
      INTEGER(4),INTENT(IN)  :: NAT    
      COMPLEX(8),INTENT(OUT) :: DENMAT(LMNXX,LMNXX,NDIMD,NAT)
      COMPLEX(8),INTENT(OUT) :: EDENMAT(LMNXX,LMNXX,NDIMD,NAT)
      INTEGER(4)             :: NBH   !#(SUPER WAVE FUNCTIONS)
      INTEGER(4)             :: NB    !#(WAVE FUNCTIONS)
      INTEGER(4)             :: IPRO  ! PROJECTOR INDEX
      INTEGER(4)             :: ISP   ! SPECIES INDEX
      INTEGER(4)             :: LMNX  ! 
      INTEGER(4)             :: NBX   ! 
      INTEGER(4)             :: IAT,ISPIN,IKPT,LMN1,LMN2,IDIM
      COMPLEX(8),ALLOCATABLE :: PROJ(:,:,:) !(NDIM,NBH,LMNX) <PRO|PSPSI>
      COMPLEX(8),ALLOCATABLE :: DENMAT1(:,:,:)
      COMPLEX(8),ALLOCATABLE :: EDENMAT1(:,:,:)
      COMPLEX(8),ALLOCATABLE :: LAGR(:,:)
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      COMPLEX(8)             :: CSVAR1,CSVAR2
      INTEGER(4)             :: NTASKS,THISTASK
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IB
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      INTEGER(4)             :: LMNXMAX
      LOGICAL(4)             :: TUSEACCLAGR
      REAL(8)                :: ACCEL_DENMAT_FLOPS
#ENDIF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
#ENDIF
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$DENMAT')
                              CALL TIMING$CLOCKON('W:DENMAT')
      IF(NDIMD_.NE.NDIMD.OR.NAT.NE.MAP%NAT) THEN 
        CALL ERROR$MSG('ARRAY SIZE INCONSISTENT')
        CALL ERROR$STOP('WAVES$DENMAT')
      END IF
!
!     ==========================================================================
!     ==  GET OCCUPATIONS FROM DYNOCC OBJECT                                  ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_DENMAT_OCC_SETUP' &
     &    ,INT(NBX,KIND=8),INT(NKPTL,KIND=8),INT(NSPIN,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      CALL MPE$QUERY('K',NTASKS,THISTASK)
      DENMAT(:,:,:,:)=(0.D0,0.D0)
      EDENMAT(:,:,:,:)=(0.D0,0.D0)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_DENMAT_INIT' &
     &    ,INT(LMNXX,KIND=8),INT(NDIMD,KIND=8),INT(NAT,KIND=8),0_8 &
     &    ,0.D0,16.D0*REAL(LMNXX,KIND=8)*REAL(LMNXX,KIND=8) &
     &    *REAL(NDIMD,KIND=8)*REAL(NAT,KIND=8)*2.D0 &
     &    ,ACCEL_T1-ACCEL_T0)
#ENDIF
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          NBH=THIS%NBH
          NB=THIS%NB
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          ALLOCATE(LAGR(NB,NB))
          DO IB=1,NB
            LAGR(:,IB)=THIS%RLAM0(:,IB)*OCC(IB,IKPT,ISPIN)
          ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          CALL ACCELPROFILE$ADD('PAW_DENMAT_LAGR_SETUP' &
     &      ,INT(NB,KIND=8),INT(NB,KIND=8),1_8,0_8 &
     &      ,0.D0,16.D0*REAL(NB,KIND=8)*REAL(NB,KIND=8) &
     &      ,ACCEL_T1-ACCEL_T0)
#ENDIF
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          LMNXMAX=0
          DO IAT=1,NAT
            ISP=MAP%ISP(IAT)
            LMNXMAX=MAX(LMNXMAX,MAP%LMNX(ISP))
          ENDDO
          ACCEL_DENMAT_FLOPS=32.D0*REAL(NBH,KIND=8)*REAL(NBH,KIND=8) &
     &       *REAL(LMNXMAX,KIND=8)*REAL(NDIM,KIND=8) &
     &       +16.D0*REAL(NBH,KIND=8)*REAL(LMNXMAX,KIND=8) &
     &       *REAL(LMNXMAX,KIND=8)*REAL(NDIM,KIND=8)*REAL(NDIM,KIND=8)
          TUSEACCLAGR=TINV &
     &       .AND.CPPAW_CUBLAS_ACC_DENMAT_ENERGY_ENABLED &
     &            (ACCEL_DENMAT_FLOPS)
          IF(TUSEACCLAGR) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL ACCELPROFILE$ADD('ACC_COPY_DENMAT_LAGR_IN' &
     &        ,INT(NB,KIND=8),INT(NB,KIND=8),1_8,0_8 &
     &        ,0.D0,16.D0*REAL(NB,KIND=8)*REAL(NB,KIND=8),0.D0)
#ENDIF
!$ACC ENTER DATA COPYIN(LAGR(1:NB,1:NB))
          END IF
#ENDIF
          IPRO=1
          DO IAT=1,NAT
            ISP=MAP%ISP(IAT)
            LMNX=MAP%LMNX(ISP)
            IF(MOD(IAT-1,NTASKS)+1.NE.THISTASK) THEN  
!             -- PARALLELIZATION. MAY NOT BE MOST EFFICIENT, 
!             -- BUT ALLOWS COMBINE OVER K-POINTS AND "G-GROUP"
              IPRO=IPRO+LMNX
              CYCLE
            END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
            ALLOCATE(PROJ(NDIM,NBH,LMNX))
            PROJ(:,:,:)=THIS%PROJ(:,:,IPRO:IPRO-1+LMNX)
            ALLOCATE(DENMAT1(LMNX,LMNX,NDIMD))
            ALLOCATE(EDENMAT1(LMNX,LMNX,NDIMD))
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL ACCELPROFILE$NOW(ACCEL_T1)
            CALL ACCELPROFILE$ADD('PAW_DENMAT_SITE_SETUP' &
     &        ,INT(LMNX,KIND=8),INT(NBH,KIND=8),INT(NDIM,KIND=8) &
     &        ,INT(NDIMD,KIND=8),0.D0 &
     &        ,16.D0*REAL(LMNX,KIND=8)*REAL(NBH,KIND=8) &
     &        *REAL(NDIM,KIND=8) &
     &        ,ACCEL_T1-ACCEL_T0)
            CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!!$PRINT*,"========  WAVES$DENMAT  ======IKPT=",IKPT,' IAT=',IAT,' IPRO=',IPRO
!!$DO IB=1,NB
!!$  IF(OCC(IB,IKPT,ISPIN).LT.1.D-5) CYCLE
!!$  WRITE(*,FMT='(I3,40("(",2F10.5,")"))')IB,THIS%PROJ(:,IB,IPRO:IPRO-1+LMNX)
!!$ENDDO
            CALL WAVES_DENMAT(NDIM,NBH,NB,LMNX,OCC(1,IKPT,ISPIN),LAGR &
     &                       ,PROJ,DENMAT1,EDENMAT1)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL ACCELPROFILE$NOW(ACCEL_T1)
            CALL ACCELPROFILE$ADD('PAW_DENMAT_SITE_KERNEL' &
     &        ,INT(LMNX,KIND=8),INT(NBH,KIND=8),INT(NB,KIND=8) &
     &        ,INT(NDIM,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
            CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
            IF(NDIM.EQ.1) THEN
              DENMAT(1:LMNX,1:LMNX,ISPIN,IAT) &
     &                   =DENMAT(1:LMNX,1:LMNX,ISPIN,IAT)+DENMAT1(:,:,1)
              EDENMAT(1:LMNX,1:LMNX,ISPIN,IAT) &
     &                   =EDENMAT(1:LMNX,1:LMNX,ISPIN,IAT)+EDENMAT1(:,:,1)
            ELSE
              DENMAT(1:LMNX,1:LMNX,:,IAT) &
     &                       =DENMAT(1:LMNX,1:LMNX,:,IAT)+DENMAT1(:,:,:)
              EDENMAT(1:LMNX,1:LMNX,:,IAT) &
     &                       =EDENMAT(1:LMNX,1:LMNX,:,IAT)+EDENMAT1(:,:,:)
            ENDIF
            DEALLOCATE(DENMAT1)
            DEALLOCATE(EDENMAT1)
            DEALLOCATE(PROJ)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL ACCELPROFILE$NOW(ACCEL_T1)
            CALL ACCELPROFILE$ADD('PAW_DENMAT_SITE_ACCUM' &
     &        ,INT(LMNX,KIND=8),INT(NBH,KIND=8),INT(NDIMD,KIND=8) &
     &        ,0_8,0.D0 &
     &        ,16.D0*REAL(LMNX,KIND=8)*REAL(LMNX,KIND=8) &
     &        *REAL(NDIMD,KIND=8)*4.D0 &
     &        ,ACCEL_T1-ACCEL_T0)
#ENDIF
            IPRO=IPRO+LMNX
          ENDDO
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          IF(TUSEACCLAGR) THEN
!$ACC EXIT DATA DELETE(LAGR(1:NB,1:NB))
          END IF
#ENDIF
          DEALLOCATE(LAGR)
        ENDDO
      ENDDO
!     == THE PROJECTIONS ARE IDENTICAL AND COMPLETE FOR EACH K-GROUP
!     == EACH K-GROUP HOLDS ONLY THE WAVE FUNCTIONS BELONGING TO IT.
!     == EACH PROCESSOR OF EACH K-GROUP ONLY ADDS UP A FRACTION OF THE PROJECTIONS
!     == THEREFORE THERE IS NO DOUBLE COUNTING BY SUMMING OVER THE MONOMER
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      CALL MPE$COMBINE('MONOMER','+',DENMAT)
      CALL MPE$COMBINE('MONOMER','+',EDENMAT)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_DENMAT_COMBINE' &
     &    ,INT(LMNXX,KIND=8),INT(NDIMD,KIND=8),INT(NAT,KIND=8) &
     &    ,INT(NTASKS,KIND=8),0.D0 &
     &    ,16.D0*REAL(LMNXX,KIND=8)*REAL(LMNXX,KIND=8) &
     &    *REAL(NDIMD,KIND=8)*REAL(NAT,KIND=8)*2.D0 &
     &    ,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  CONVERT SPIN-UP AND SPIN-DOWN DENSITY MATRIX INTO                   ==
!     ==  TOTAL AND SPIN DENSITY MATRICES                                     ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      IF(NSPIN.EQ.2) THEN
        DO IAT=1,NAT
          ISP=MAP%ISP(IAT)
          LMNX=MAP%LMNX(ISP)
          DO LMN1=1,LMNX
             DO LMN2=1,LMNX
              CSVAR1=DENMAT(LMN1,LMN2,1,IAT)
              CSVAR2=DENMAT(LMN1,LMN2,2,IAT)
              DENMAT(LMN1,LMN2,1,IAT)=CSVAR1+CSVAR2   ! TOTAL 
              DENMAT(LMN1,LMN2,2,IAT)=CSVAR1-CSVAR2   ! SPIN
!
              CSVAR1=EDENMAT(LMN1,LMN2,1,IAT)
              CSVAR2=EDENMAT(LMN1,LMN2,2,IAT)
              EDENMAT(LMN1,LMN2,1,IAT)=CSVAR1+CSVAR2   ! TOTAL 
              EDENMAT(LMN1,LMN2,2,IAT)=CSVAR1-CSVAR2   ! SPIN
            ENDDO
          ENDDO
        ENDDO
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_DENMAT_SPIN_CONVERT' &
     &    ,INT(LMNXX,KIND=8),INT(NDIMD,KIND=8),INT(NAT,KIND=8) &
     &    ,INT(NSPIN,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  PRINT FOR TEST                                                      ==
!     ==========================================================================
      IF(TPRINT) THEN
        WRITE(*,FMT='("TEST PRINT FROM WAVES$DENMAT")')
        DO IAT=1,NAT
          ISP=MAP%ISP(IAT)
          LMNX=MAP%LMNX(ISP)
          DO IDIM=1,NDIMD
            WRITE(*,FMT='("DENMAT FOR ATOM ",I3," AND SPIN ",I3)')IAT,IDIM
            DO LMN1=1,LMNX
              WRITE(*,FMT='("R",I3,100F10.5)')LMN1,REAL(DENMAT(LMN1,:LMNX,IDIM,IAT))
              WRITE(*,FMT='("I",I3,100F10.5)')LMN1,AIMAG(DENMAT(LMN1,:LMNX,IDIM,IAT))
            ENDDO
          ENDDO
        ENDDO
      END IF
!!$!
!!$!     =======================================================================
!!$!     ==  NOW CALCULATE OFF-SITE DENSITY MATRIX                            ==
!!$!     =======================================================================
!!$      CALL LOCALIZE$SWITCHCHIDENMAT(.TRUE.)
!!$      ALLOCATE(XK(3,NKPTL))
!!$      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
!!$      DO IKPT=1,NKPTL
!!$        DO ISPIN=1,NSPIN
!!$          CALL WAVES_SELECTWV(IKPT,ISPIN)
!!$          CALL PLANEWAVE$SELECT(GSET%ID)
!!$          NBH=THIS%NBH
!!$          NB=THIS%NB
!!$          TINV=GSET%TINV
!!$          CALL LOCALIZE$CHIDENMAT(XK(:,IKPT),TINV,NB,OCC(1,IKPT,ISPIN),NDIM,NBH,MAP%NPRO,THIS%PROJ)
!!$        ENDDO
!!$      ENDDO
!!$      DEALLOCATE(XK)
!!$      CALL LOCALIZE$REPORTCHIDENMAT(6)
!!$      CALL LOCALIZE$NLEXCHANGE()
!!$STOP 'FORCED STOP IN WAVES$DENMAT'
!
      DEALLOCATE(OCC)
                              CALL TIMING$CLOCKOFF('W:DENMAT')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$DENMAT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$TONTBO()
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES ONE-CENTER DENSITY MATRIX FROM THE ACTUAL                 **
!     **  PROJECTIONS <PRO|PSI> IN THE WAVES OBJECT                           **
!     **                                                                      **
!     **  ONE-CENTER DENSITY MATRICES                                         **
!     **  SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                             **
!     **  SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)                      **
!     **  NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z)        **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      INTEGER(4)             :: NBH    !#(SUPER WAVE FUNCTIONS)
      INTEGER(4)             :: NPRO   !#(PROJECTOR FUNCTIONS)
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      INTEGER(4)             :: ISPIN,IKPT,IB
      LOGICAL(4)             :: TON
      INTEGER(4)             :: NTASKS,THISTASK,COUNT
      INTEGER(4)             :: NORB
!     **************************************************************************
      CALL LMTO$GETL4('ON',TON)
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('WAVES$TONTBO')
                              CALL TIMING$CLOCKON('W:TONTBO')
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      NPRO=MAP%NPRO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      CALL LMTO$GETI4('NLOCORB',NORB)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          NBH=THIS%NBH
          IF(.NOT.ASSOCIATED(THIS%TBC_NEW))ALLOCATE(THIS%TBC_NEW(NDIM,NBH,NORB))
          CALL LMTO$PROJTONTBO_NEW('FWRD',XK(:,IKPT),NDIM,NBH,NPRO,THIS%PROJ &
     &                                                       ,NORB,THIS%TBC_NEW)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  PRINT FOR TESTING                                                   ==
!     ==========================================================================
      IF(TPRINT) THEN
        WRITE(*,FMT='(82("="),T20," TIGHT-BINDING ORBITAL COEFFICIENTS ")')
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            NBH=THIS%NBH
            WRITE(*,FMT='(82("="),T20,"  IKPT ",I5,"ISPIN=",I2,"  ")')IKPT,ISPIN
            DO IB=1,NBH
              WRITE(*,FMT='("PROJ",I5,100F10.5)')2*IB-1,REAL(THIS%PROJ(1,IB,:))
              WRITE(*,FMT='("PROJ",I5,100F10.5)')2*IB,AIMAG(THIS%PROJ(1,IB,:))
            ENDDO
            DO IB=1,NBH
              WRITE(*,FMT='("TBC ",I5,100F10.5)')2*IB-1,REAL(THIS%TBC_NEW(1,IB,:))
              WRITE(*,FMT='("TBC ",I5,100F10.5)')2*IB,AIMAG(THIS%TBC_NEW(1,IB,:))
            ENDDO
          ENDDO
        ENDDO
        CALL ERROR$MSG('REGULAR TERMINATION AFTER PRINTOUT')
        CALL ERROR$MSG('UNDO BY SETTING TPRINT=.FALSE.')
        CALL ERROR$STOP('WAVES$TONTBO')
      END IF
                              CALL TIMING$CLOCKOFF('W:TONTBO')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$TONTBO
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$FROMNTBO()
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES ONE-CENTER DENSITY MATRIX FROM THE ACTUAL                 **
!     **  PROJECTIONS <PRO|PSI> IN THE WAVES OBJECT                           **
!     **                                                                      **
!     **  ONE-CENTER DENSITY MATRICES                                         **
!     **  SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                             **
!     **  SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)                      **
!     **  NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z)        **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TOLD=.FALSE.
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      INTEGER(4)             :: NBH   !#(SUPER WAVE FUNCTIONS)
      INTEGER(4)             :: NPRO   !#(PROJECTOR FUNCTIONS)
      INTEGER(4)             :: NORB  
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      COMPLEX(8),ALLOCATABLE :: HPROJ(:,:,:)
      INTEGER(4)             :: ISPIN,IKPT,IB
      LOGICAL(4)             :: TON
!     **************************************************************************
      CALL LMTO$GETL4('ON',TON)
      IF(.NOT.TON) RETURN
      CALL LMTO$GETL4('THTBC',TON) !HTBC HAS NOT BEEN CALCULATED
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('WAVES$FROMNTBO')
                              CALL TIMING$CLOCKON('W:FROMNTBO')
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      NPRO=MAP%NPRO
!
!     ==========================================================================
!     ==  TRANSFORM THIS$HTBC FROM NTBOS TO PAW PARTIAL WAVES                 ==
!     ==========================================================================
!!!! PARALLELIZE LOOP WITH RESPECT TO STATES.
      CALL LMTO$GETI4('NLOCORB',NORB)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          NBH=THIS%NBH
          ALLOCATE(HPROJ(NDIM,NBH,NPRO))
          IF(.NOT.ASSOCIATED(THIS%HTBC_NEW)) THEN
!         == THIS%HTBC IS EVALUATED IN LMTO_NTBODENMATDER_NEW
            CALL ERROR$MSG('INTERNAL ERROR:')
            CALL ERROR$MSG('LMTO$PROJTONTBO USES NON-ASSOCIATED POINTER')
            CALL ERROR$MSG('CHECK CALL IN PAW_WAVES1.F90')
            CALL ERROR$STOP('WAVES$FROMNTBO')
          END IF
          CALL LMTO$PROJTONTBO_NEW('BACK',XK(:,IKPT),NDIM,NBH,NPRO,HPROJ &
     &                                                      ,NORB,THIS%HTBC_NEW)
          DEALLOCATE(THIS%HTBC_NEW)
          ALLOCATE(THIS%HTBC_NEW(NDIM,NBH,NPRO))
          THIS%HTBC_NEW=HPROJ
          DEALLOCATE(HPROJ)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  PRINT FOR TESTING                                                   ==
!     ==========================================================================
      IF(TPRINT) THEN
        WRITE(*,FMT='(82("="),T20," TIGHT-BINDING ORBITAL POTENTIALS ")')
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            NBH=THIS%NBH
            WRITE(*,FMT='(80("="),T20," IKPT=",I5," ISPIN=",I2," ")')IKPT,ISPIN
            DO IB=1,NBH
              WRITE(*,FMT='(80("-"),T20,"  IB=",I5,"  ")')IB
              WRITE(*,FMT='("RE(HTBC) ",100F10.5)')REAL(THIS%HTBC_NEW(1,IB,:))
              WRITE(*,FMT='("IM(HTBC) ",100F10.5)')AIMAG(THIS%HTBC_NEW(1,IB,:))
            ENDDO
          ENDDO
        ENDDO
        WRITE(*,*)'FORCED STOP IN WAVES$FROMNTBO'
        STOP 'FORCED STOP IN WAVES$FROMNTBO'
      END IF
                              CALL TIMING$CLOCKOFF('W:FROMNTBO')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$FROMNTBO
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_DEDFFROMNTBO(NB,NKPTL,NSPIN,EIG)
!     **************************************************************************
!     ** ADDS TERMS TO DEDF, WHICH ARE NOT TAKEN CARE OF WITH <PSI|HPSI>      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB
      INTEGER(4),INTENT(IN)    :: NKPTL
      INTEGER(4),INTENT(IN)    :: NSPIN
      REAL(8)   ,INTENT(INOUT) :: EIG(NB,NKPTL,NSPIN)
      REAL(8)                  :: DEIG(NB,NKPTL,NSPIN)
      LOGICAL(4)               :: TON
!     **************************************************************************
      CALL LMTO$GETL4('ON',TON)
      IF(.NOT.TON) RETURN
      CALL LMTO$GETL4('THTBC',TON) !HTBC HAS NOT BEEN CALCULATED
      IF(.NOT.TON) RETURN
!
      CALL LMTO$DEDF(NB,NKPTL,NSPIN,DEIG)
      EIG=EIG+DEIG
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_DENMAT(NDIM,NBH,NB,LMNX,OCC,LAGR,PROPSI &
     &                       ,DENMAT,EDENMAT)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATE ONE-CENTER DENSITY MATRIX                                  **
!     **             <P_I|PSI_N>F_N<PSI_N|P_J>                                **
!     **  FROM PROJECTIONS <P|PSI> AND THE OCCUPATIONS                        **
!     **                                                                      **
!     **  FOR NDIM=2,NDIMD=4:  SPINOR WAVE FUNCTIONS                          **
!     **  THE DENSITY MATRIX ON RETURN CONTAINS THE TOTAL DENSITY             **
!     **  AND THE X,Y,Z SPIN DENSITY                                          **
!     **                                                                      **
!     **                                                                      **
!     **  SUPERWAVE FUNCTIONS ARE DEFINED AS: PSI=PSI1+I*PSI2                 **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1999)***********
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &        CPPAW_CUBLAS_ACC_DENMAT_ENERGY_ENABLED
#ENDIF
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NDIM   ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN) :: NBH    ! #(WAVE FUNCTIONS
      INTEGER(4),INTENT(IN) :: NB     ! #(STATES)
      INTEGER(4),INTENT(IN) :: LMNX   ! #(PROJECTORS ON THIS SITE)
      REAL(8)   ,INTENT(IN) :: OCC(NB)! OCCUPATIONS
      COMPLEX(8),INTENT(IN) :: LAGR(NB,NB)     ! LAGRANGE/F TIMES OCCUPATION
      COMPLEX(8),INTENT(IN) :: PROPSI(NDIM,NBH,LMNX) !<PRO|PSI>
      COMPLEX(8),INTENT(OUT):: DENMAT(LMNX,LMNX,NDIM**2)
      COMPLEX(8),INTENT(OUT):: EDENMAT(LMNX,LMNX,NDIM**2)
      COMPLEX(8)            :: DENMAT1(LMNX,LMNX,NDIM,NDIM)
      COMPLEX(8)            :: EDENMAT1(LMNX,LMNX,NDIM,NDIM)
      COMPLEX(8)            :: FUNC(LMNX,NDIM)
      LOGICAL(4)            :: TINV
      INTEGER(4)            :: LMN1,LMN2,IDIM1,IDIM2,IB,IB1,IB2
      REAL(8)               :: SVAR1,SVAR2
      COMPLEX(8)            :: CSVAR,CFACR,CFACI,CFAC1,CFAC2
      INTEGER(4)            :: IDIM,NDIMD
      INTEGER(4)            :: NFILO
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      LOGICAL(4)             :: TUSEACCENERGY
      REAL(8)                :: ACCEL_ENERGY_FLOPS
#ENDIF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
#ENDIF
!     **************************************************************************
      NDIMD=NDIM**2
      DENMAT(:,:,:)=(0.D0,0.D0)
      EDENMAT(:,:,:)=(0.D0,0.D0)
!
!     ==========================================================================
!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT        ==
!     ==========================================================================
      CALL PLANEWAVE$GETL4('TINV',TINV)
      IF(TINV) THEN
        IF(NBH.NE.(NB+1)/2) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_DENMAT')
        END IF
      ELSE 
        IF(NBH.NE.NB) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_DENMAT')
        END IF
      END IF
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      ACCEL_ENERGY_FLOPS=32.D0*REAL(NBH,KIND=8)*REAL(NBH,KIND=8) &
     &                  *REAL(LMNX,KIND=8)*REAL(NDIM,KIND=8) &
     &                  +16.D0*REAL(NBH,KIND=8) &
     &                  *REAL(LMNX,KIND=8)*REAL(LMNX,KIND=8) &
     &                  *REAL(NDIM,KIND=8)*REAL(NDIM,KIND=8)
      TUSEACCENERGY=TINV &
     &             .AND.CPPAW_CUBLAS_ACC_DENMAT_ENERGY_ENABLED &
     &                  (ACCEL_ENERGY_FLOPS)
#ENDIF
!
!     ==========================================================================
!     ==  SUM UP THE DENSITY MATRIX                                           ==
!     ==========================================================================
!     == SUPER WAVE FUNCTIONS FOR TINV=TRUE
!     ==   <P|PSI1+I*PSI2>(F1+F2)/2<PSI1-I*PSI2|P>
!     == + <P|PSI1+I*PSI2>(F1-F2)/2<PSI1+I*PSI2|P>
!     == = <P|PSI1>F1<PSI1|P>+<P|PSI2>F2<PSI2|P>
!     == + I[<P|PSI2>F1<PSI1|P>-<P|PSI1>F2<PSI2|P>]
!     == IMAGINARY PART IS DROPPED
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      DENMAT1(:,:,:,:)=(0.D0,0.D0)
      DO IB=1,NBH
!       == FIND OCCUPATION OF THE STATE ========================================
        IF(TINV) THEN
          SVAR1=0.5D0*(OCC(2*IB-1)+OCC(2*IB))
          SVAR2=0.5D0*(OCC(2*IB-1)-OCC(2*IB))
        ELSE
          SVAR1=OCC(IB)
        END IF
        DO LMN1=1,LMNX
          DO IDIM1=1,NDIM
            IF(TINV) THEN
!             == FUNC=(F1+F2)/2<PSI1-I*PSI2|P>+(F1-F2)/2<PSI1+I*PSI2|P>
              FUNC(LMN1,IDIM1)=SVAR1*CONJG(PROPSI(IDIM1,IB,LMN1)) &
     &                        +SVAR2*PROPSI(IDIM1,IB,LMN1)
            ELSE 
              FUNC(LMN1,IDIM1)=SVAR1*CONJG(PROPSI(IDIM1,IB,LMN1))
            END IF
          ENDDO
        ENDDO
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            DO IDIM1=1,NDIM
              DO IDIM2=1,NDIM
                DENMAT1(LMN1,LMN2,IDIM1,IDIM2)=DENMAT1(LMN1,LMN2,IDIM1,IDIM2) &
     &                      +PROPSI(IDIM1,IB,LMN1)*FUNC(LMN2,IDIM2)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     == DENSITY MATRIX FOR REAL WAVE FUNCTIONS IS REAL. IMAGINARY PART
!     == CONTAINS CRAP DUE TO SUPER WAVE FUNCTIONS
      IF(TINV) THEN
        DENMAT1(:,:,:,:)=REAL(DENMAT1(:,:,:,:),KIND=8)
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_DENMAT_DENSITY_LOOP' &
     &    ,INT(LMNX,KIND=8),INT(NBH,KIND=8),INT(NDIM,KIND=8) &
     &    ,INT(NB,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     == NOW SUM UP EDENMAT  <P|PSITILDE>*LAMBDA*<PSITILDE|P>                 ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      EDENMAT1(:,:,:,:)=(0.D0,0.D0)
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      IF(TUSEACCENERGY) THEN
        CALL WAVES_DENMAT_ENERGY_TINV_ACC(NDIM,NBH,NB,LMNX,LAGR &
     &                                   ,PROPSI,EDENMAT1)
      ELSE
#ENDIF
      IF(TINV) THEN
        DO IB1=1,NBH
          FUNC(:,:)=(0.D0,0.D0)
          DO IB2=1,NBH
            CFACR=LAGR(2*IB1-1,2*IB2-1)+CI*LAGR(2*IB1,2*IB2-1)
            CFACI=LAGR(2*IB1-1,2*IB2)+CI*LAGR(2*IB1,2*IB2)
            CFAC1=0.5D0*(CFACR-CI*CFACI)
            CFAC2=0.5D0*(CFACR+CI*CFACI)
            DO IDIM=1,NDIM
              FUNC(:,IDIM)=FUNC(:,IDIM)+CFAC1*PROPSI(IDIM,IB2,:) &
     &                                 +CFAC2*CONJG(PROPSI(IDIM,IB2,:))
            ENDDO
          ENDDO
          DO IDIM2=1,NDIM
            DO IDIM1=1,NDIM
              DO LMN2=1,LMNX
                DO LMN1=1,LMNX
                  EDENMAT1(LMN1,LMN2,IDIM1,IDIM2)=EDENMAT1(LMN1,LMN2,IDIM1,IDIM2) &
     &                              +CONJG(PROPSI(IDIM1,IB1,LMN1))*FUNC(LMN2,IDIM2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!       == DENSITY MATRIX FOR REAL WAVE FUNCTIONS IS REAL. IMAGINARY PART
!       == CONTAINS CRAP DUE TO SUPER WAVE FUNCTIONS
        EDENMAT1(:,:,:,:)=REAL(EDENMAT1(:,:,:,:),KIND=8)
      ELSE
        DO IB1=1,NBH
          DO IB2=1,NBH
            DO IDIM=1,NDIM
              FUNC(:,IDIM)=FUNC(:,IDIM)+LAGR(IB1,IB2)*CONJG(PROPSI(IDIM,IB2,:))
            ENDDO
          ENDDO
          DO IDIM2=1,NDIM
            DO IDIM1=1,NDIM
              DO LMN2=1,LMNX
                DO LMN1=1,LMNX
                  EDENMAT1(LMN1,LMN2,IDIM1,IDIM2)=EDENMAT1(LMN1,LMN2,IDIM1,IDIM2) &
     &                                 +PROPSI(IDIM1,IB1,LMN1)*FUNC(LMN2,IDIM2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      END IF
#ENDIF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_DENMAT_ENERGY_LOOP' &
     &    ,INT(LMNX,KIND=8),INT(NBH,KIND=8),INT(NDIM,KIND=8) &
     &    ,INT(NB,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     == MAP DENSITY MATRIX ONTO TOTAL AND SPIN DENSITY                       ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      IF(NDIM.EQ.1) THEN  !== TOTAL DENSITY ====================================
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            DENMAT(LMN1,LMN2,1)=DENMAT1(LMN1,LMN2,1,1)
            EDENMAT(LMN1,LMN2,1)=EDENMAT1(LMN1,LMN2,1,1)
          ENDDO
        ENDDO
      ELSE IF(NDIM.EQ.2) THEN  !== TOTAL DENSITY, X,Y,Z SPIN DENSITY ===========
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            CSVAR=DENMAT1(LMN1,LMN2,1,1)+DENMAT1(LMN1,LMN2,2,2)
            DENMAT(LMN1,LMN2,1) =CSVAR
            CSVAR=DENMAT1(LMN1,LMN2,1,2)+DENMAT1(LMN1,LMN2,2,1)
            DENMAT(LMN1,LMN2,2) =CSVAR
!           == DENMAT=TR[SIGMA*DENMAT1]=SUM_{I,J} SIGMA_{I,J}*DENMAT{J,I}
!           == HENCE THE MINUS SIGN FOR THE ANTISYMMETRIC SIGMA
            CSVAR=CI*(DENMAT1(LMN1,LMN2,1,2)-DENMAT1(LMN1,LMN2,2,1))
            DENMAT(LMN1,LMN2,3) =CSVAR
            CSVAR=DENMAT1(LMN1,LMN2,1,1)-DENMAT1(LMN1,LMN2,2,2)
            DENMAT(LMN1,LMN2,4) =CSVAR
!
            CSVAR=EDENMAT1(LMN1,LMN2,1,1)+EDENMAT1(LMN1,LMN2,2,2)
            EDENMAT(LMN1,LMN2,1) =CSVAR
            CSVAR=EDENMAT1(LMN1,LMN2,1,2)+EDENMAT1(LMN1,LMN2,2,1)
            EDENMAT(LMN1,LMN2,2) =CSVAR
            CSVAR=CI*(EDENMAT1(LMN1,LMN2,1,2)-EDENMAT1(LMN1,LMN2,2,1))
            EDENMAT(LMN1,LMN2,3) =CSVAR
            CSVAR=EDENMAT1(LMN1,LMN2,1,1)-EDENMAT1(LMN1,LMN2,2,2)
            EDENMAT(LMN1,LMN2,4) =CSVAR
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == SYMMETRIZE DENSITY MATRIX (DENSITY MATRIX MUST BE HERMITEAN)         ==
!     ==========================================================================
      DO IDIM=1,NDIMD
        DO LMN1=1,LMNX
          DO LMN2=LMN1+1,LMNX
            CSVAR=0.5D0*(DENMAT(LMN1,LMN2,IDIM)+CONJG(DENMAT(LMN2,LMN1,IDIM)))
            DENMAT(LMN1,LMN2,IDIM)=CSVAR
            DENMAT(LMN2,LMN1,IDIM)=CONJG(CSVAR)
!
            CSVAR=0.5D0*(EDENMAT(LMN1,LMN2,IDIM)+CONJG(EDENMAT(LMN2,LMN1,IDIM)))
            EDENMAT(LMN1,LMN2,IDIM)=CSVAR
            EDENMAT(LMN2,LMN1,IDIM)=CONJG(CSVAR)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == SYMMETRIZE DENSITY MATRIX WITH RESPECT TO TIME INVERSION             ==
!     == FOR COLLINEAR CALCULATIONS PSI(K)=CONJG(PSI(-K)). THEREFORE THE      ==                     
!     == DENSITY MATRIX IS REAL AFTER SUMMING OVER K-POINTS. THIS IS NOT TRUE ==
!     == FOR TRANSPORT CALCULATIONS. THIS SYMMETRY FAILS FRO SPIN-ORBIT       ==
!     == AND NON-COLLINEAR CALCULATIONS WHERE AN EXPLICIT MAGNETIC FIELD      ==
!     == IS PRESENT.
!     ==========================================================================
      IF(NDIM.EQ.1) THEN
        DENMAT(:,:,1)=REAL(DENMAT(:,:,1),KIND=8)
        EDENMAT(:,:,1)=REAL(EDENMAT(:,:,1),KIND=8)
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_DENMAT_MAP_SYM' &
     &    ,INT(LMNX,KIND=8),INT(NBH,KIND=8),INT(NDIMD,KIND=8) &
     &    ,INT(NB,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     == PRINTOUT FOR TEST                                                    ==
!     ==========================================================================
      IF(TPR) THEN
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        DO IDIM=1,NDIMD
          WRITE(NFILO,FMT='("WAVES_DENMAT: DENMAT FOR IDIM= ",I2)') IDIM
          DO LMN1=1,LMNX
            WRITE(NFILO,FMT='(9F10.6)')REAL(DENMAT(LMN1,:LMNX,IDIM),KIND=8)
          ENDDO
        ENDDO
        DO IB=1,NBH
          WRITE(NFILO,FMT='("PROPSI",I5,50F10.6)')IB,PROPSI(1,IB,:LMNX)
        ENDDO
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      SUBROUTINE WAVES_DENMAT_ENERGY_TINV_ACC(NDIM,NBH,NB,LMNX,LAGR &
     &                                       ,PROPSI,EDENMAT1)
!     **************************************************************************
!     **  OpenACC diagnostic for the time-inversion one-center energy matrix. **
!     **  It is intentionally opt-in because it copies the Lambda block per   **
!     **  site; the profiling rows show whether a broader resident rewrite is **
!     **  worth the extra data-structure work.                                **
!     **************************************************************************
      USE OPENACC
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NDIM
      INTEGER(4),INTENT(IN) :: NBH
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LMNX
      COMPLEX(8),INTENT(IN) :: LAGR(NB,NB)
      COMPLEX(8),INTENT(IN) :: PROPSI(NDIM,NBH,LMNX)
      COMPLEX(8),INTENT(OUT):: EDENMAT1(LMNX,LMNX,NDIM,NDIM)
      COMPLEX(8),ALLOCATABLE :: FUNCACC(:,:,:)
      INTEGER(4)            :: LMN1,LMN2,IDIM1,IDIM2,IB1,IB2
      REAL(8)               :: RSUM,ISUM
      COMPLEX(8)            :: CFACR,CFACI,CFAC1,CFAC2,CTMP
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)               :: ACCEL_T0
      REAL(8)               :: ACCEL_T1
      REAL(8)               :: ACCEL_KERNEL_T0
      REAL(8)               :: ACCEL_KERNEL_T1
      REAL(8)               :: ACCEL_BYTES
      REAL(8)               :: ACCEL_FLOPS
      REAL(8)               :: ACCEL_COPY_TIME
      LOGICAL(4)            :: TLAGRPRESENT
      LOGICAL(4)            :: TPROPSIPRESENT
#ENDIF
!     **************************************************************************
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
      ACCEL_FLOPS=32.D0*REAL(NBH,KIND=8)*REAL(NBH,KIND=8) &
     &            *REAL(LMNX,KIND=8)*REAL(NDIM,KIND=8) &
     &            +16.D0*REAL(NBH,KIND=8) &
     &            *REAL(LMNX,KIND=8)*REAL(LMNX,KIND=8) &
     &            *REAL(NDIM,KIND=8)*REAL(NDIM,KIND=8)
      TLAGRPRESENT=ACC_IS_PRESENT(LAGR)
      TPROPSIPRESENT=ACC_IS_PRESENT(PROPSI)
      ACCEL_BYTES=16.D0*REAL(LMNX,KIND=8)*REAL(LMNX,KIND=8) &
     &            *REAL(NDIM,KIND=8)*REAL(NDIM,KIND=8)
      IF(.NOT.TLAGRPRESENT) THEN
        ACCEL_BYTES=ACCEL_BYTES &
     &      +16.D0*REAL(NB,KIND=8)*REAL(NB,KIND=8)
      ELSE
        CALL ACCELPROFILE$ADD('ACC_PRESENT_DENMAT_LAGR' &
     &      ,INT(NB,KIND=8),INT(NB,KIND=8),1_8,0_8,0.D0,0.D0,0.D0)
      END IF
      IF(.NOT.TPROPSIPRESENT) THEN
        ACCEL_BYTES=ACCEL_BYTES &
     &      +16.D0*REAL(NDIM,KIND=8)*REAL(NBH,KIND=8) &
     &      *REAL(LMNX,KIND=8)
      END IF
#ENDIF
      ALLOCATE(FUNCACC(LMNX,NDIM,NBH))
!$ACC DATA PRESENT_OR_COPYIN(LAGR(1:NB,1:NB) &
!$ACC& ,PROPSI(1:NDIM,1:NBH,1:LMNX)) &
!$ACC& CREATE(FUNCACC(1:LMNX,1:NDIM,1:NBH)) &
!$ACC& COPYOUT(EDENMAT1(1:LMNX,1:LMNX,1:NDIM,1:NDIM))
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_KERNEL_T0)
#ENDIF
!$ACC PARALLEL LOOP COLLAPSE(3) GANG &
!$ACC& PRIVATE(IB2,RSUM,ISUM,CFACR,CFACI,CFAC1,CFAC2,CTMP) &
!$ACC& PRESENT(LAGR,PROPSI,FUNCACC)
      DO IB1=1,NBH
        DO IDIM2=1,NDIM
          DO LMN2=1,LMNX
            RSUM=0.D0
            ISUM=0.D0
!$ACC LOOP VECTOR REDUCTION(+:RSUM,ISUM)
            DO IB2=1,NBH
              CFACR=LAGR(2*IB1-1,2*IB2-1) &
     &             +CI*LAGR(2*IB1,2*IB2-1)
              CFACI=LAGR(2*IB1-1,2*IB2) &
     &             +CI*LAGR(2*IB1,2*IB2)
              CFAC1=0.5D0*(CFACR-CI*CFACI)
              CFAC2=0.5D0*(CFACR+CI*CFACI)
              CTMP=CFAC1*PROPSI(IDIM2,IB2,LMN2) &
     &            +CFAC2*CONJG(PROPSI(IDIM2,IB2,LMN2))
              RSUM=RSUM+REAL(CTMP,KIND=8)
              ISUM=ISUM+AIMAG(CTMP)
            ENDDO
            FUNCACC(LMN2,IDIM2,IB1)=CMPLX(RSUM,ISUM,KIND=8)
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL LOOP
!$ACC PARALLEL LOOP COLLAPSE(4) GANG &
!$ACC& PRIVATE(IB1,RSUM,ISUM,CTMP) &
!$ACC& PRESENT(PROPSI,FUNCACC,EDENMAT1)
      DO IDIM2=1,NDIM
        DO IDIM1=1,NDIM
          DO LMN2=1,LMNX
            DO LMN1=1,LMNX
              RSUM=0.D0
              ISUM=0.D0
!$ACC LOOP VECTOR REDUCTION(+:RSUM,ISUM)
              DO IB1=1,NBH
                CTMP=CONJG(PROPSI(IDIM1,IB1,LMN1)) &
     &              *FUNCACC(LMN2,IDIM2,IB1)
                RSUM=RSUM+REAL(CTMP,KIND=8)
                ISUM=ISUM+AIMAG(CTMP)
              ENDDO
              EDENMAT1(LMN1,LMN2,IDIM1,IDIM2) &
     &          =CMPLX(RSUM,ISUM,KIND=8)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL LOOP
!$ACC PARALLEL LOOP COLLAPSE(4) PRESENT(EDENMAT1)
      DO IDIM2=1,NDIM
        DO IDIM1=1,NDIM
          DO LMN2=1,LMNX
            DO LMN1=1,LMNX
              EDENMAT1(LMN1,LMN2,IDIM1,IDIM2) &
     &          =REAL(EDENMAT1(LMN1,LMN2,IDIM1,IDIM2),KIND=8)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL LOOP
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
!$ACC WAIT
      CALL ACCELPROFILE$NOW(ACCEL_KERNEL_T1)
      CALL ACCELPROFILE$ADD('ACC_KERNEL_DENMAT_ENERGY_TINV' &
     &    ,INT(LMNX,KIND=8),INT(NBH,KIND=8),INT(NDIM,KIND=8) &
     &    ,INT(NB,KIND=8),ACCEL_FLOPS,0.D0 &
     &    ,ACCEL_KERNEL_T1-ACCEL_KERNEL_T0)
#ENDIF
!$ACC END DATA
      DEALLOCATE(FUNCACC)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      ACCEL_COPY_TIME=(ACCEL_T1-ACCEL_T0) &
     &               -(ACCEL_KERNEL_T1-ACCEL_KERNEL_T0)
      IF(ACCEL_COPY_TIME.LT.0.D0) ACCEL_COPY_TIME=0.D0
      CALL ACCELPROFILE$ADD('ACC_COPY_DENMAT_ENERGY_TINV' &
     &    ,INT(LMNX,KIND=8),INT(NBH,KIND=8),INT(NDIM,KIND=8) &
     &    ,INT(NB,KIND=8),0.D0,ACCEL_BYTES,ACCEL_COPY_TIME)
#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
#ENDIF
      SUBROUTINE WAVES$OFFSITEDENMAT()
!     **************************************************************************
!     ** EVALUATES THE DENSITY MATRIX FOR A NEIGHBORLIST                      **
!     **                                                                      **
!     ** NEIGHBORLIST INCLUDES ALL ATOM PAIRS WITH DISTANCE SMALLER THAN      **
!     ** THE SUM OF COVALENT RADII TIME SCALERCUT                             **
!     ** SCALERCUT=5. IS GOOD FOR THE HUBBARD MODEL WITH LATTICE CONSTANT=3\AA**
!     **************************************************************************
      USE MPE_MODULE
      USE PERIODICTABLE_MODULE
      USE WAVES_MODULE, ONLY: MAP &
     &                       ,NDIMD &
     &                       ,OSDENMAT &
     &                       ,OSHAMIL &
     &                       ,SCALERCUT !RADIUS SCALE FACTOR FOR NEIGHBORLIST
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPRNN=.FALSE.  !PRINT NEIGHBORLIST ENTRIES
      INTEGER(4),PARAMETER   :: NNXPERATOM=400
      INTEGER(4)             :: NNX
      INTEGER(4)             :: NND
      INTEGER(4),ALLOCATABLE :: NNLIST(:,:) !(5,NNX)
      INTEGER(4)             :: NAT
      REAL(8)                :: RBAS(3,3) !LATTICE VECTORS
      REAL(8)   ,ALLOCATABLE :: R0(:,:)   !(3,NAT) ATOMIC POSITIONS
      REAL(8)   ,ALLOCATABLE :: RC(:)     !(NAT) RADIUS OF NEIGHBORLIST
      REAL(8)                :: SVAR
      INTEGER(4)             :: IAT1,IAT2,N1,N2
      INTEGER(4)             :: IAT,ISP,NN
      INTEGER(4)             :: NTASKS,THISTASK
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
#ENDIF
!     **************************************************************************
                                          CALL TRACE$PUSH('WAVES_OFFSITEDENMAT')
      NAT=MAP%NAT
!
!     ==========================================================================
!     == CLEAN DENSITY MATRIX STRUCTURE                                       ==
!     ==========================================================================
      IF(ALLOCATED(OSDENMAT)) THEN
        NND=SIZE(OSDENMAT)
        DO NN=1,NND
          IF(ALLOCATED(OSDENMAT(NN)%MAT)) DEALLOCATE(OSDENMAT(NN)%MAT)
        ENDDO
        DEALLOCATE(OSDENMAT)
      END IF
      IF(ALLOCATED(OSHAMIL)) THEN
        NND=SIZE(OSHAMIL)
        DO NN=1,NND
          IF(ALLOCATED(OSHAMIL(NN)%MAT)) DEALLOCATE(OSHAMIL(NN)%MAT)
        ENDDO
        DEALLOCATE(OSHAMIL)
      END IF
!
!     ==========================================================================
!     == OBTAIN ATOMIC STRUCTURE                                              ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     == SET UP NEIGHBORLIST (ENCODED IN NNLIST)                              ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      ALLOCATE(RC(NAT))
      DO IAT=1,NAT
        ISP=MAP%ISP(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',SVAR)
        CALL SETUP$UNSELECT()
        CALL PERIODICTABLE$GET(SVAR,'R(COV)',RC(IAT))
        RC(IAT)=RC(IAT)*SCALERCUT
      ENDDO
      NNX=NNXPERATOM*NAT
      ALLOCATE(NNLIST(5,NNX))
!     == SYNCHRONIZE NEIGHBORLIST OVER ALL MPI PROCESSES. NEIGBORLIST MAY ======
!     == DIFFER ON DIFFERENT TASKS DUE TO AMBIGUOUS DISTANCE CRITERION. ========
      IF(THISTASK.EQ.1)CALL LMTO$NEIGHBORLIST(RBAS,NAT,R0,RC,NNX,NND,NNLIST)
      CALL MPE$BROADCAST('MONOMER',1,NND)
      CALL MPE$BROADCAST('MONOMER',1,NNLIST(:,:NND))
      DEALLOCATE(RC)
      DEALLOCATE(R0)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OFFDEN_NEIGHBOR' &
     &    ,INT(NAT,KIND=8),INT(NND,KIND=8),INT(NTASKS,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  REPORT NEIGHBORLIST                                                 ==
!     ==========================================================================
      IF(TPRNN) THEN
        WRITE(*,FMT= &
     &       '(82("="),T5," NEIGHBORLIST ENTRIES IN WAVES$OFFSITEDENMAT ")')
        DO NN=1,NND
          WRITE(*,FMT='(I10," ATOM1=",I4," ATOM2=",I4," IT=",3I3)') &
     &                NN,NNLIST(:,NN)
        ENDDO
      END IF
!
!     ==========================================================================
!     == ALLOCATE DENSITY MATRIX
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      ALLOCATE(OSDENMAT(NND))
      DO NN=1,NND
        IAT1=NNLIST(1,NN)
        IAT2=NNLIST(2,NN)
        N1=MAP%LMNX(MAP%ISP(IAT1))
        N2=MAP%LMNX(MAP%ISP(IAT2))
        OSDENMAT(NN)%IAT1=IAT1
        OSDENMAT(NN)%IAT2=IAT2
        OSDENMAT(NN)%IT=NNLIST(3:5,NN)
        OSDENMAT(NN)%N1=N1
        OSDENMAT(NN)%N2=N2
        OSDENMAT(NN)%N3=NDIMD  !(TOTAL),(TOTAL,Z),(TOTAL,X,Y,Z)
        ALLOCATE(OSDENMAT(NN)%MAT(N1,N2,NDIMD))
        OSDENMAT(NN)%MAT(:,:,:)=0.D0
      ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OFFDEN_ALLOC' &
     &    ,INT(NAT,KIND=8),INT(NND,KIND=8),INT(NDIMD,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     == EVALUATE DENSITY MATRIX
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      CALL WAVES_SUMMUPOFFSITEDENMAT()
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OFFDEN_SUM' &
     &    ,INT(NAT,KIND=8),INT(NND,KIND=8),INT(NDIMD,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
                                                                CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_SUMMUPOFFSITEDENMAT()
!     **************************************************************************
!     **  CONSTRUCT DENSITY MATRIX IN A NTBO BASIS                            **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY: NKPTL &
     &                       ,NSPIN &
     &                       ,NDIM  &
     &                       ,NDIMD &
     &                       ,THIS  &
     &                       ,MAP   &
     &                       ,WAVES_SELECTWV &
     &                       ,GSET  &
     &                       ,OSDENMAT  ! OFF-SITE DENSITY MATRIX
!     &                       ,RSPACEOP ! (FROM RSPACEOP_MODULE)
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)             :: NAT
      INTEGER(4)             :: N1,N2
      INTEGER(4)             :: NND
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: NB,NBH,NBX
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      INTEGER(4),ALLOCATABLE :: IPRO1(:)
      INTEGER(4),ALLOCATABLE :: NPROAT(:)
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      INTEGER(4)             :: IAT,NN,II,ISP,IPRO,IKPT,ISPIN,I,J,IBH,IB
      REAL(8)                :: SVAR
      REAL(8)                :: F1,F2
      LOGICAL(4)             :: TINV
      LOGICAL(4)             :: TOFFDENBLAS
      LOGICAL(4)             :: TOFFDENCUBLAS
      LOGICAL(4)             :: TOFFDENBATCH
      LOGICAL(4)             :: TOFFDENDEVICEPACK
      LOGICAL(4)             :: TOFFDENDEVICEACCUM
      INTEGER(4)             :: OFFDENBATCHSIZE
      INTEGER(4)             :: IAT1,IAT2,IT(3),I0,J0,IDIM,JDIM
      COMPLEX(8)             :: EIKR,C1(NDIM),C2(NDIM),CSVAR22(NDIM,NDIM)
      INTEGER(4)             :: NTASKS,THISTASK,ICOUNT
      CHARACTER(16)          :: ENVVAL
      INTEGER(4)             :: ENVSTAT
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
#ENDIF
!     **************************************************************************
                                    CALL TRACE$PUSH('WAVES_SUMMUPOFFSITEDENMAT')
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
      TOFFDENBLAS=.FALSE.
      TOFFDENCUBLAS=.FALSE.
      TOFFDENBATCH=.FALSE.
      TOFFDENDEVICEPACK=.FALSE.
      TOFFDENDEVICEACCUM=.FALSE.
      OFFDENBATCHSIZE=64
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_OFFDEN_LOCAL',ENVVAL &
     &                             ,STATUS=ENVSTAT)
      IF(ENVSTAT.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_OFFDEN_BLAS',ENVVAL &
     &                               ,STATUS=ENVSTAT)
      END IF
      IF(ENVSTAT.EQ.0) THEN
        SELECT CASE(TRIM(ADJUSTL(ENVVAL)))
        CASE('1','T','t','TRUE','true','True','YES','yes','ON','on')
          TOFFDENBLAS=.TRUE.
        CASE DEFAULT
          TOFFDENBLAS=.FALSE.
        END SELECT
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_OFFDEN_CUBLAS',ENVVAL &
     &                             ,STATUS=ENVSTAT)
      IF(ENVSTAT.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_OFFDEN',ENVVAL &
     &                               ,STATUS=ENVSTAT)
      END IF
      IF(ENVSTAT.EQ.0) THEN
        SELECT CASE(TRIM(ADJUSTL(ENVVAL)))
        CASE('1','T','t','TRUE','true','True','YES','yes','ON','on')
          TOFFDENBLAS=.TRUE.
          TOFFDENCUBLAS=.TRUE.
        CASE DEFAULT
          TOFFDENCUBLAS=.FALSE.
        END SELECT
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_OFFDEN_CUBLAS_BATCH' &
     &                             ,ENVVAL,STATUS=ENVSTAT)
      IF(ENVSTAT.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_OFFDEN_BATCH' &
     &                               ,ENVVAL,STATUS=ENVSTAT)
      END IF
      IF(ENVSTAT.EQ.0) THEN
        SELECT CASE(TRIM(ADJUSTL(ENVVAL)))
        CASE('1','T','t','TRUE','true','True','YES','yes','ON','on')
          TOFFDENBLAS=.TRUE.
          TOFFDENCUBLAS=.TRUE.
          TOFFDENBATCH=.TRUE.
        CASE DEFAULT
          TOFFDENBATCH=.FALSE.
        END SELECT
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_OFFDEN_DEVICE_PACK' &
     &                             ,ENVVAL,STATUS=ENVSTAT)
      IF(ENVSTAT.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_OFFDEN_DEVICE_PACK' &
     &                               ,ENVVAL,STATUS=ENVSTAT)
      END IF
      IF(ENVSTAT.EQ.0) THEN
        SELECT CASE(TRIM(ADJUSTL(ENVVAL)))
        CASE('1','T','t','TRUE','true','True','YES','yes','ON','on')
          TOFFDENBLAS=.TRUE.
          TOFFDENCUBLAS=.TRUE.
          TOFFDENBATCH=.TRUE.
          TOFFDENDEVICEPACK=.TRUE.
        CASE DEFAULT
          TOFFDENDEVICEPACK=.FALSE.
        END SELECT
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_OFFDEN_DEVICE_ACCUM' &
     &                             ,ENVVAL,STATUS=ENVSTAT)
      IF(ENVSTAT.NE.0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('CPPAW_CUBLAS_ACC_OFFDEN_DEVICE_ACCUM' &
     &                               ,ENVVAL,STATUS=ENVSTAT)
      END IF
      IF(ENVSTAT.EQ.0) THEN
        SELECT CASE(TRIM(ADJUSTL(ENVVAL)))
        CASE('1','T','t','TRUE','true','True','YES','yes','ON','on')
          TOFFDENBLAS=.TRUE.
          TOFFDENCUBLAS=.TRUE.
          TOFFDENBATCH=.TRUE.
          TOFFDENDEVICEPACK=.TRUE.
          TOFFDENDEVICEACCUM=.TRUE.
        CASE DEFAULT
          TOFFDENDEVICEACCUM=.FALSE.
        END SELECT
      END IF
      CALL GET_ENVIRONMENT_VARIABLE('CPPAW_GPU_OFFDEN_BATCH_SIZE' &
     &                             ,ENVVAL,STATUS=ENVSTAT)
      IF(ENVSTAT.EQ.0) READ(ENVVAL,*,ERR=101,END=101) OFFDENBATCHSIZE
  101 CONTINUE
      IF(OFFDENBATCHSIZE.LT.1) OFFDENBATCHSIZE=64
!
!     ==========================================================================
!     ==  CONSTRUCT INDEX ARRAYS                                              ==
!     ==========================================================================
      NAT=MAP%NAT
      ALLOCATE(IPRO1(NAT))
      ALLOCATE(NPROAT(NAT))
      IPRO=1
      DO IAT=1,NAT
        ISP=MAP%ISP(IAT)
        IPRO1(IAT)=IPRO
        NPROAT(IAT)=MAP%LMNX(ISP)
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OFFDEN_SUM_SETUP' &
     &    ,INT(NKPTL,KIND=8),INT(NSPIN,KIND=8),INT(NAT,KIND=8),0_8 &
     &    ,0.D0,8.D0*REAL(NBX,KIND=8)*REAL(NKPTL,KIND=8) &
     &    *REAL(NSPIN,KIND=8),ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  ADD UP DENSITY MATRIX                                               ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      NND=SIZE(OSDENMAT)
      DO NN=1,NND
        OSDENMAT(NN)%MAT=0.D0
      ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OFFDEN_SUM_ZERO' &
     &    ,INT(NND,KIND=8),INT(NDIMD,KIND=8),0_8,0_8 &
     &    ,0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
      NPRO=MAP%NPRO
      CALL MPE$QUERY('K',NTASKS,THISTASK)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          NBH=THIS%NBH
          NB=THIS%NB
          IF(TOFFDENBATCH.AND.TINV.AND.NDIM.EQ.1) THEN
            CALL WAVES_OFFDEN_TINV_NDIM1_CUBLAS_BATCH(NND,NPRO &
     &           ,NPROAT,IPRO1,NTASKS,THISTASK,NBH,NB,NDIMD &
     &           ,NSPIN,ISPIN,OFFDENBATCHSIZE,OCC(1,IKPT,ISPIN) &
     &           ,XK(1,IKPT),THIS%PROJ(1,1,1),OSDENMAT &
     &           ,TOFFDENDEVICEPACK,TOFFDENDEVICEACCUM)
            CYCLE
          END IF
          ICOUNT=0
          DO NN=1,NND
            ICOUNT=ICOUNT+1
            IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) CYCLE
            IAT1=OSDENMAT(NN)%IAT1
            IAT2=OSDENMAT(NN)%IAT2
            IT  =OSDENMAT(NN)%IT
!           SVAR=-2.D0*PI*SUM(XK(:,IKPT)*REAL(IT,KIND=8))
            SVAR=2.D0*PI*SUM(XK(:,IKPT)*REAL(IT,KIND=8))
            EIKR=EXP(CI*SVAR)  !<P_{R+T}|PSI>=<P_R|PSI>*EIKR
            I0=IPRO1(IAT1)-1
            J0=IPRO1(IAT2)-1
            IF(TOFFDENBLAS.AND.TINV.AND.NDIM.EQ.1) THEN
              CALL WAVES_OFFDEN_TINV_NDIM1_BLAS(NPROAT(IAT1) &
     &             ,NPROAT(IAT2),NBH,NB,NDIMD,NSPIN,ISPIN &
     &             ,OCC(1,IKPT,ISPIN),EIKR,THIS%PROJ(1,1,I0+1) &
     &             ,THIS%PROJ(1,1,J0+1),OSDENMAT(NN)%MAT &
     &             ,TOFFDENCUBLAS)
              CYCLE
            END IF
            DO I=1,NPROAT(IAT1)
              DO J=1,NPROAT(IAT2)
!
                IF(TINV) THEN
                  CSVAR22(:,:)=(0.D0,0.D0)
                  DO IBH=1,NBH
                    F1=OCC(2*IBH-1,IKPT,ISPIN)
                    F2=OCC(2*IBH,IKPT,ISPIN)
                    C1(:)=THIS%PROJ(:,IBH,I0+I)
                    C2(:)=THIS%PROJ(:,IBH,J0+J)*EIKR   ! EXP(-I*K*T)
                    DO JDIM=1,NDIM
                      CSVAR22(:,JDIM)=CSVAR22(:,JDIM) &
     &                               +0.5D0*((F1+F2)*C1(:)*CONJG(C2(JDIM)) &
     &                                      +(F1-F2)*C1(:)*C2(JDIM))
                    ENDDO
                  ENDDO
                  CSVAR22=REAL(CSVAR22) ! IMAG(CSVAR) CONTAINS CRAP 
                                        !  DUE TO SUPER WAVE FUNCTIONS
                ELSE
                  CSVAR22(:,:)=(0.D0,0.D0)
                  DO IB=1,NB
                    F1=OCC(IB,IKPT,ISPIN)
                    C1(:)=THIS%PROJ(:,IB,I0+I)
                    C2(:)=THIS%PROJ(:,IB,J0+J)*EIKR ! EXP(-I*K*T)
                    DO JDIM=1,NDIM
                      CSVAR22(:,JDIM)=CSVAR22(:,JDIM)+F1*C1(:)*CONJG(C2(JDIM))
                    ENDDO
                  ENDDO
                END IF
!
!           == DISTRIBUTE ONTO DENSITY MATRIX ENTRIES ==========================
!           == D(IDIMD)=SUM_{IDIM,JDIM} D(IDIM,JDIM)*PAULI_{IDIMD}(JDIM,IDIM) ==
!           == IDIMD IN {TOTAL,X,Y,Z}; IDIM IN {UP,DOWN} =======================
!           == TRANSFORMATION MUST BE CONSISTENT WITH WAVES_DENMAT =============
                IF(NSPIN.EQ.1) THEN
                  IF(NDIM.EQ.1) THEN !NON-SPIN-POLARIZED
                    OSDENMAT(NN)%MAT(I,J,1)=OSDENMAT(NN)%MAT(I,J,1) &
     &                                       +REAL(CSVAR22(1,1))
                  ELSE ! NONCOLLINEAR
                    OSDENMAT(NN)%MAT(I,J,1)=OSDENMAT(NN)%MAT(I,J,1) &
     &                                  +REAL(CSVAR22(1,1)+CSVAR22(2,2))
                    OSDENMAT(NN)%MAT(I,J,2)=OSDENMAT(NN)%MAT(I,J,2) &
     &                                  +REAL(CSVAR22(1,2)+CSVAR22(2,1))
                    OSDENMAT(NN)%MAT(I,J,3)=OSDENMAT(NN)%MAT(I,J,3) &
     &                                  -AIMAG(CSVAR22(1,2)-CSVAR22(2,1))
                    OSDENMAT(NN)%MAT(I,J,4)=OSDENMAT(NN)%MAT(I,J,4) &
     &                                  +REAL(CSVAR22(1,1)-CSVAR22(2,2))
                  END IF
                ELSE IF(NSPIN.EQ.2) THEN !SPIN POLARIZED
                  IF(ISPIN.EQ.1) THEN
                    OSDENMAT(NN)%MAT(I,J,1)=OSDENMAT(NN)%MAT(I,J,1) &
     &                                       +REAL(CSVAR22(1,1))
                    OSDENMAT(NN)%MAT(I,J,2)=OSDENMAT(NN)%MAT(I,J,2) &
     &                                       +REAL(CSVAR22(1,1))
                  ELSE
                    OSDENMAT(NN)%MAT(I,J,1)=OSDENMAT(NN)%MAT(I,J,1) &
     &                                       +REAL(CSVAR22(1,1))
                    OSDENMAT(NN)%MAT(I,J,2)=OSDENMAT(NN)%MAT(I,J,2) &
     &                                       -REAL(CSVAR22(1,1))
                  END IF
                END IF
              ENDDO   !J  (PROJECTORS)
            ENDDO   !I  (PROJECTORS)
          ENDDO   !NN
        ENDDO   !ISPIN
      ENDDO  !IKPT
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OFFDEN_SUM_LOCAL' &
     &    ,INT(NND,KIND=8),INT(NPRO,KIND=8),INT(NB,KIND=8) &
     &    ,INT(NTASKS,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  SUM OVER MONOMER INCLUDES ALSO THE KPOINT SUM                       ==
!     ==========================================================================
!     == THE FOLLOWING STATEMENT CAUSED MPI PROBLEMS (MPI_MESSAGE_TRUNCATE) ====
!     == BECAUSE THE SHAPE OF OSDENMAT(NN)%MAT DIFFERED ON DIFFERENT TASKS =====
!     == THIS WAS DUE TO DIFFERING NEIGHBORLISTS ON DIFFERENT TASKS ============
!     == THE PROBLEM IS FIXED BY BROADCASTING THE NEIGHBORLIST =================
!     == IN WAVES$OFFSITEDENMAT(). =============================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      DO NN=1,NND
!       -- ONCE, THE CODE FAILED WITHIN MPI IN THE FOLLOWING CALL. THE FAILURE -
!       -- WAS NOT DETERMINISTIC -----------------------------------------------
        CALL MPE$COMBINE('MONOMER','+',OSDENMAT(NN)%MAT)
      ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OFFDEN_SUM_COMBINE' &
     &    ,INT(NND,KIND=8),INT(NDIMD,KIND=8),INT(NTASKS,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,FMT='(82("="),T10," NON-LOCAL DENSITY MATRIX  ")')
        WRITE(*,FMT='(82("="),T10," FROM WAVES_SUMMUPOFFSITEDENMAT ")')
        DO NN=1,NND
          IAT1=OSDENMAT(NN)%IAT1
          IAT2=OSDENMAT(NN)%IAT2
          IT=OSDENMAT(NN)%IT
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
     &                                                         IAT1,IAT2,IT
          N1=OSDENMAT(NN)%N1
          N2=OSDENMAT(NN)%N2
          DO I=1,1 !OSDENMAT(NN)%N3
            DO J=1,N1 
              WRITE(*,FMT='(I3,300F10.3)')I,OSDENMAT(NN)%MAT(J,:,I)
            ENDDO
            WRITE(*,FMT='(82("-"))')
          ENDDO
        ENDDO
        CALL ERROR$MSG('FORCED STOP AFTER PRINTING DENSITY MATRIX')
        CALL ERROR$STOP('WAVES_SUMMUPOFFSITEDENMAT')
      END IF
                                                                CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_OFFDEN_TINV_NDIM1_BLAS(N1,N2,NBH,NB,NDIMD &
     &                                       ,NSPIN,ISPIN,OCC,EIKR &
     &                                       ,PROJ1,PROJ2,MAT &
     &                                       ,TCUBLAS)
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &        CPPAW_CUBLAS_ACC_ZGEMM_NT_COPY
#ENDIF
!     **************************************************************************
!     **  BLAS diagnostic for the inversion-symmetric, scalar off-site        **
!     **  density-matrix contraction. It is opt-in and intentionally limited  **
!     **  to NDIM=1, where the original tensor-to-spin mapping is scalar.     **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      INTEGER(4),INTENT(IN) :: NBH
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: NSPIN
      INTEGER(4),INTENT(IN) :: ISPIN
      REAL(8)   ,INTENT(IN) :: OCC(NB)
      COMPLEX(8),INTENT(IN) :: EIKR
      COMPLEX(8),INTENT(IN) :: PROJ1(NBH,N1)
      COMPLEX(8),INTENT(IN) :: PROJ2(NBH,N2)
      REAL(8)   ,INTENT(INOUT) :: MAT(N1,N2,NDIMD)
      LOGICAL(4),INTENT(IN) :: TCUBLAS
      COMPLEX(8),ALLOCATABLE :: A(:,:)
      COMPLEX(8),ALLOCATABLE :: B(:,:)
      COMPLEX(8),ALLOCATABLE :: WORK(:,:)
      COMPLEX(8)            :: ONE
      COMPLEX(8)            :: ZERO
      REAL(8)               :: FPLUS
      REAL(8)               :: FMINUS
      INTEGER(4)            :: I,J,IBH
      LOGICAL(4)            :: TCUBLAS_USED
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)               :: ACCEL_T0
      REAL(8)               :: ACCEL_T1
      REAL(8)               :: ACCEL_FLOPS
      CHARACTER(32)         :: ACCEL_GEMM_NAME
#ENDIF
!     **************************************************************************
      ONE=(1.D0,0.D0)
      ZERO=(0.D0,0.D0)
      TCUBLAS_USED=.FALSE.
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      ALLOCATE(A(N1,NBH))
      ALLOCATE(B(N2,NBH))
      ALLOCATE(WORK(N1,N2))
      DO IBH=1,NBH
        FPLUS=0.5D0*(OCC(2*IBH-1)+OCC(2*IBH))
        FMINUS=0.5D0*(OCC(2*IBH-1)-OCC(2*IBH))
        DO I=1,N1
          A(I,IBH)=PROJ1(IBH,I)
        ENDDO
        DO J=1,N2
          B(J,IBH)=FPLUS*CONJG(PROJ2(IBH,J))*CONJG(EIKR) &
     &            +FMINUS*PROJ2(IBH,J)*EIKR
        ENDDO
      ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OFFDEN_BLAS_PACK' &
     &    ,INT(N1,KIND=8),INT(N2,KIND=8),INT(NBH,KIND=8),0_8 &
     &    ,0.D0,16.D0*REAL(NBH,KIND=8)*REAL(N1+N2,KIND=8) &
     &    ,ACCEL_T1-ACCEL_T0)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      IF(TCUBLAS) THEN
        CALL CPPAW_CUBLAS_ACC_ZGEMM_NT_COPY(N1,N2,NBH,A,B,WORK &
     &                                     ,TCUBLAS_USED)
      END IF
#ENDIF
      IF(.NOT.TCUBLAS_USED) THEN
        CALL ZGEMM('N','T',N1,N2,NBH,ONE,A,N1,B,N2,ZERO,WORK,N1)
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      ACCEL_FLOPS=8.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8) &
     &           *REAL(NBH,KIND=8)
      IF(TCUBLAS_USED) THEN
        ACCEL_GEMM_NAME='CUBLAS_ZGEMM_OFFDEN_TINV_NDIM1'
      ELSE
        ACCEL_GEMM_NAME='ZGEMM_OFFDEN_TINV_NDIM1'
      END IF
      CALL ACCELPROFILE$ADD(ACCEL_GEMM_NAME &
     &    ,INT(N1,KIND=8),INT(N2,KIND=8),INT(NBH,KIND=8),0_8 &
     &    ,ACCEL_FLOPS,16.D0*(REAL(NBH,KIND=8)*REAL(N1+N2,KIND=8) &
     &    +REAL(N1,KIND=8)*REAL(N2,KIND=8)),ACCEL_T1-ACCEL_T0)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      IF(NSPIN.EQ.1) THEN
        DO J=1,N2
          DO I=1,N1
            MAT(I,J,1)=MAT(I,J,1)+REAL(WORK(I,J),KIND=8)
          ENDDO
        ENDDO
      ELSE IF(NSPIN.EQ.2) THEN
        IF(ISPIN.EQ.1) THEN
          DO J=1,N2
            DO I=1,N1
              MAT(I,J,1)=MAT(I,J,1)+REAL(WORK(I,J),KIND=8)
              MAT(I,J,2)=MAT(I,J,2)+REAL(WORK(I,J),KIND=8)
            ENDDO
          ENDDO
        ELSE
          DO J=1,N2
            DO I=1,N1
              MAT(I,J,1)=MAT(I,J,1)+REAL(WORK(I,J),KIND=8)
              MAT(I,J,2)=MAT(I,J,2)-REAL(WORK(I,J),KIND=8)
            ENDDO
          ENDDO
        END IF
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OFFDEN_BLAS_ACCUM' &
     &    ,INT(N1,KIND=8),INT(N2,KIND=8),INT(NDIMD,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
      DEALLOCATE(WORK)
      DEALLOCATE(B)
      DEALLOCATE(A)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_OFFDEN_TINV_NDIM1_CUBLAS_BATCH(NND,NPRO &
     &          ,NPROAT,IPRO1,NTASKS,THISTASK,NBH,NB,NDIMD,NSPIN &
     &          ,ISPIN,BATCHSIZE,OCC,XK,PROJ,OSDENMAT,TDEVICEPACK &
     &          ,TDEVICEACCUM)
      USE RSPACEOP_MODULE, ONLY: RSPACEMAT_TYPE
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      USE OPENACC
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &        CPPAW_CUBLAS_ACC_SHOULD_USE_OFFDEN &
     &       ,CPPAW_CUBLAS_ACC_ZGEMM_NT_COPY &
     &       ,CPPAW_CUBLAS_ACC_ZGEMM_NT_PRESENT
#ENDIF
!     **************************************************************************
!     **  Chunked cuBLAS diagnostic for scalar, inversion-symmetric off-site  **
!     **  density-matrix contractions. Neighbors sharing the same first atom  **
!     **  and second-projector size are stacked into one wide ZGEMM.          **
!     **  The optional device-pack mode keeps PROJ/OCC in one OpenACC region, **
!     **  packs A/B on the GPU, and can accumulate WORK into a flat real      **
!     **  off-site matrix buffer before one host-visible result copy.         **
!     **************************************************************************
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      INTEGER(4),INTENT(IN)  :: NND
      INTEGER(4),INTENT(IN)  :: NPRO
      INTEGER(4),INTENT(IN)  :: NPROAT(*)
      INTEGER(4),INTENT(IN)  :: IPRO1(*)
      INTEGER(4),INTENT(IN)  :: NTASKS
      INTEGER(4),INTENT(IN)  :: THISTASK
      INTEGER(4),INTENT(IN)  :: NBH
      INTEGER(4),INTENT(IN)  :: NB
      INTEGER(4),INTENT(IN)  :: NDIMD
      INTEGER(4),INTENT(IN)  :: NSPIN
      INTEGER(4),INTENT(IN)  :: ISPIN
      INTEGER(4),INTENT(IN)  :: BATCHSIZE
      REAL(8)   ,INTENT(IN)  :: OCC(NB)
      REAL(8)   ,INTENT(IN)  :: XK(3)
      COMPLEX(8),INTENT(IN)  :: PROJ(NBH,NPRO)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT) :: OSDENMAT(NND)
      LOGICAL(4),INTENT(IN)  :: TDEVICEPACK
      LOGICAL(4),INTENT(IN)  :: TDEVICEACCUM
      LOGICAL(4),ALLOCATABLE :: DONE(:)
      INTEGER(4),ALLOCATABLE :: IDX(:)
      INTEGER(4),ALLOCATABLE :: J0ARR(:)
      INTEGER(4),ALLOCATABLE :: MOFFARR(:)
      INTEGER(4),ALLOCATABLE :: MATOFF(:)
      COMPLEX(8),ALLOCATABLE :: A(:,:)
      COMPLEX(8),ALLOCATABLE :: B(:,:)
      COMPLEX(8),ALLOCATABLE :: WORK(:,:)
      COMPLEX(8),ALLOCATABLE :: EIKRARR(:)
      REAL(8)   ,ALLOCATABLE :: MATFLAT(:)
      COMPLEX(8)            :: ONE
      COMPLEX(8)            :: ZERO
      COMPLEX(8)            :: EIKR
      REAL(8)               :: SVAR
      REAL(8)               :: FPLUS
      REAL(8)               :: FMINUS
      REAL(8)               :: OFFDENFLOPS
      INTEGER(4)            :: NN0
      INTEGER(4)            :: NN
      INTEGER(4)            :: K
      INTEGER(4)            :: I
      INTEGER(4)            :: J
      INTEGER(4)            :: JDIM
      INTEGER(4)            :: IBH
      INTEGER(4)            :: IAT1
      INTEGER(4)            :: IAT2
      INTEGER(4)            :: I0
      INTEGER(4)            :: J0
      INTEGER(4)            :: N1
      INTEGER(4)            :: N2
      INTEGER(4)            :: COUNT
      INTEGER(4)            :: MAXBATCH
      INTEGER(4)            :: IOFF
      INTEGER(4)            :: MOFF
      INTEGER(4)            :: NTOT
      LOGICAL(4)            :: TDEVICEACTIVE
      LOGICAL(4)            :: TDEVICEACCUMACTIVE
      LOGICAL(4)            :: TDEVICEFLAT
      LOGICAL(4)            :: TCUBLAS_USED
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)               :: ACCEL_T0
      REAL(8)               :: ACCEL_T1
      REAL(8)               :: ACCEL_FLOPS
      CHARACTER(32)         :: ACCEL_GEMM_NAME
#ENDIF
!     **************************************************************************
      ONE=(1.D0,0.D0)
      ZERO=(0.D0,0.D0)
      MAXBATCH=MAX(1,BATCHSIZE)
      ALLOCATE(DONE(NND))
      ALLOCATE(IDX(MAXBATCH))
      DONE=.FALSE.
      TDEVICEFLAT=.FALSE.
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      TDEVICEFLAT=TDEVICEPACK.AND.TDEVICEACCUM
#ENDIF
      NTOT=0
      IF(TDEVICEFLAT) THEN
        ALLOCATE(MATOFF(NND))
        MATOFF=0
        DO NN=1,NND
          IF(MOD(NN-1,NTASKS).NE.THISTASK-1) CYCLE
          IAT1=OSDENMAT(NN)%IAT1
          IAT2=OSDENMAT(NN)%IAT2
          N1=NPROAT(IAT1)
          N2=NPROAT(IAT2)
          MATOFF(NN)=NTOT+1
          NTOT=NTOT+N1*N2*NDIMD
        ENDDO
        IF(NTOT.GT.0) THEN
          ALLOCATE(MATFLAT(NTOT))
          MATFLAT=0.D0
!$ACC ENTER DATA COPYIN(MATFLAT(1:NTOT))
        ELSE
          TDEVICEFLAT=.FALSE.
          DEALLOCATE(MATOFF)
        END IF
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      IF(TDEVICEPACK) THEN
#IF DEFINED(CPPVAR_CUBLAS_ACC)
        IF(ACC_IS_PRESENT(PROJ)) THEN
          CALL ACCELPROFILE$ADD('ACC_PRESENT_OFFDEN_DPACK_PROJ' &
     &        ,INT(NBH,KIND=8),INT(NPRO,KIND=8),INT(NB,KIND=8),0_8 &
     &        ,0.D0,0.D0,0.D0)
        ELSE
#ENDIF
          CALL ACCELPROFILE$ADD('ACC_COPY_OFFDEN_DPACK_PROJ_IN' &
     &        ,INT(NBH,KIND=8),INT(NPRO,KIND=8),INT(NB,KIND=8),0_8 &
     &        ,0.D0,16.D0*REAL(NBH,KIND=8)*REAL(NPRO,KIND=8) &
     &        +8.D0*REAL(NB,KIND=8),0.D0)
#IF DEFINED(CPPVAR_CUBLAS_ACC)
        END IF
#ENDIF
      END IF
#ENDIF
!$ACC DATA PRESENT_OR_COPYIN(PROJ(1:NBH,1:NPRO),OCC(1:NB)) IF(TDEVICEPACK)
      DO NN0=1,NND
        IF(DONE(NN0)) CYCLE
        IF(MOD(NN0-1,NTASKS).NE.THISTASK-1) THEN
          DONE(NN0)=.TRUE.
          CYCLE
        END IF
        IAT1=OSDENMAT(NN0)%IAT1
        IAT2=OSDENMAT(NN0)%IAT2
        N1=NPROAT(IAT1)
        N2=NPROAT(IAT2)
        COUNT=0
        DO NN=NN0,NND
          IF(DONE(NN)) CYCLE
          IF(MOD(NN-1,NTASKS).NE.THISTASK-1) THEN
            DONE(NN)=.TRUE.
            CYCLE
          END IF
          IF(OSDENMAT(NN)%IAT1.NE.IAT1) CYCLE
          IF(NPROAT(OSDENMAT(NN)%IAT2).NE.N2) CYCLE
          COUNT=COUNT+1
          IDX(COUNT)=NN
          DONE(NN)=.TRUE.
          IF(COUNT.EQ.MAXBATCH) EXIT
        ENDDO
        IF(COUNT.EQ.0) CYCLE
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
        ALLOCATE(A(N1,NBH))
        ALLOCATE(B(N2*COUNT,NBH))
        ALLOCATE(WORK(N1,N2*COUNT))
        OFFDENFLOPS=8.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8) &
     &             *REAL(NBH,KIND=8)*REAL(COUNT,KIND=8)
        TDEVICEACTIVE=.FALSE.
        TDEVICEACCUMACTIVE=.FALSE.
#IF DEFINED(CPPVAR_CUBLAS_ACC)
        TDEVICEACTIVE=TDEVICEPACK &
     &      .AND.CPPAW_CUBLAS_ACC_SHOULD_USE_OFFDEN(OFFDENFLOPS)
#ENDIF
        I0=IPRO1(IAT1)-1
        IF(TDEVICEACTIVE) THEN
          TDEVICEACCUMACTIVE=TDEVICEFLAT
          ALLOCATE(J0ARR(COUNT))
          ALLOCATE(EIKRARR(COUNT))
          IF(TDEVICEACCUMACTIVE) ALLOCATE(MOFFARR(COUNT))
          DO K=1,COUNT
            NN=IDX(K)
            IAT2=OSDENMAT(NN)%IAT2
            J0ARR(K)=IPRO1(IAT2)-1
            SVAR=2.D0*PI*SUM(XK(:)*REAL(OSDENMAT(NN)%IT,KIND=8))
            EIKRARR(K)=EXP(CI*SVAR)
            IF(TDEVICEACCUMACTIVE) MOFFARR(K)=MATOFF(NN)
          ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          IF(TDEVICEACCUMACTIVE) THEN
!$ACC ENTER DATA CREATE(A(1:N1,1:NBH),B(1:N2*COUNT,1:NBH) &
!$ACC&                  ,WORK(1:N1,1:N2*COUNT)) &
!$ACC& COPYIN(J0ARR(1:COUNT),EIKRARR(1:COUNT),MOFFARR(1:COUNT))
          ELSE
!$ACC ENTER DATA CREATE(A(1:N1,1:NBH),B(1:N2*COUNT,1:NBH) &
!$ACC&                  ,WORK(1:N1,1:N2*COUNT)) &
!$ACC& COPYIN(J0ARR(1:COUNT),EIKRARR(1:COUNT))
          END IF
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(PROJ,A)
          DO IBH=1,NBH
            DO I=1,N1
              A(I,IBH)=PROJ(IBH,I0+I)
            ENDDO
          ENDDO
!$ACC END PARALLEL LOOP
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(IOFF,FPLUS,FMINUS) &
!$ACC& PRESENT(PROJ,OCC,J0ARR,EIKRARR,B)
          DO K=1,COUNT
            DO IBH=1,NBH
              DO J=1,N2
                IOFF=(K-1)*N2
                FPLUS=0.5D0*(OCC(2*IBH-1)+OCC(2*IBH))
                FMINUS=0.5D0*(OCC(2*IBH-1)-OCC(2*IBH))
                B(IOFF+J,IBH)=FPLUS &
     &              *CONJG(PROJ(IBH,J0ARR(K)+J))*CONJG(EIKRARR(K)) &
     &              +FMINUS*PROJ(IBH,J0ARR(K)+J)*EIKRARR(K)
              ENDDO
            ENDDO
          ENDDO
!$ACC END PARALLEL LOOP
!$ACC WAIT
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          CALL ACCELPROFILE$ADD('PAW_OFFDEN_DEVICE_PACK' &
     &        ,INT(N1,KIND=8),INT(N2,KIND=8),INT(NBH,KIND=8) &
     &        ,INT(COUNT,KIND=8),0.D0 &
     &        ,16.D0*REAL(NBH,KIND=8) &
     &        *(REAL(N1,KIND=8)+REAL(N2,KIND=8)*REAL(COUNT,KIND=8)) &
     &        ,ACCEL_T1-ACCEL_T0)
          CALL ACCELPROFILE$ADD('ACC_COPY_OFFDEN_DPACK_META_IN' &
     &        ,INT(COUNT,KIND=8),0_8,0_8,0_8,0.D0 &
     &        ,REAL(4*COUNT,KIND=8)+16.D0*REAL(COUNT,KIND=8),0.D0)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          TCUBLAS_USED=.TRUE.
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          CALL CPPAW_CUBLAS_ACC_ZGEMM_NT_PRESENT(N1,N2*COUNT,NBH &
     &        ,A,B,WORK)
#ENDIF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          CALL ACCELPROFILE$ADD('CUBLAS_ZGEMM_OFFDEN_TINV_DPACK' &
     &        ,INT(N1,KIND=8),INT(N2,KIND=8),INT(NBH,KIND=8) &
     &        ,INT(COUNT,KIND=8),OFFDENFLOPS &
     &        ,16.D0*(REAL(NBH,KIND=8)*REAL(N1+N2,KIND=8) &
     &        +REAL(N1,KIND=8)*REAL(N2,KIND=8))*REAL(COUNT,KIND=8) &
     &        ,ACCEL_T1-ACCEL_T0)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          IF(TDEVICEACCUMACTIVE) THEN
!$ACC WAIT
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(IOFF,MOFF) &
!$ACC& PRESENT(WORK,MATFLAT,MOFFARR)
            DO K=1,COUNT
              DO J=1,N2
                DO I=1,N1
                  IOFF=(K-1)*N2
                  MOFF=MOFFARR(K)+(J-1)*N1+I-1
                  IF(NSPIN.EQ.1) THEN
                    MATFLAT(MOFF)=MATFLAT(MOFF) &
     &                  +REAL(WORK(I,IOFF+J),KIND=8)
                  ELSE IF(NSPIN.EQ.2) THEN
                    IF(ISPIN.EQ.1) THEN
                      MATFLAT(MOFF)=MATFLAT(MOFF) &
     &                    +REAL(WORK(I,IOFF+J),KIND=8)
                      MATFLAT(MOFF+N1*N2)=MATFLAT(MOFF+N1*N2) &
     &                    +REAL(WORK(I,IOFF+J),KIND=8)
                    ELSE
                      MATFLAT(MOFF)=MATFLAT(MOFF) &
     &                    +REAL(WORK(I,IOFF+J),KIND=8)
                      MATFLAT(MOFF+N1*N2)=MATFLAT(MOFF+N1*N2) &
     &                    -REAL(WORK(I,IOFF+J),KIND=8)
                    END IF
                  END IF
                ENDDO
              ENDDO
            ENDDO
!$ACC END PARALLEL LOOP
!$ACC WAIT
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL ACCELPROFILE$NOW(ACCEL_T1)
            CALL ACCELPROFILE$ADD('PAW_OFFDEN_DEVICE_ACCUM' &
     &          ,INT(N1,KIND=8),INT(N2,KIND=8),INT(NDIMD,KIND=8) &
     &          ,INT(COUNT,KIND=8),0.D0 &
     &          ,8.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8) &
     &          *REAL(NDIMD,KIND=8)*REAL(COUNT,KIND=8),ACCEL_T1-ACCEL_T0)
            CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC EXIT DATA DELETE(A(1:N1,1:NBH),B(1:N2*COUNT,1:NBH) &
!$ACC&                 ,WORK(1:N1,1:N2*COUNT),J0ARR(1:COUNT) &
!$ACC&                 ,EIKRARR(1:COUNT),MOFFARR(1:COUNT))
          ELSE
!$ACC UPDATE SELF(WORK(1:N1,1:N2*COUNT))
!$ACC WAIT
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL ACCELPROFILE$NOW(ACCEL_T1)
            CALL ACCELPROFILE$ADD('ACC_COPY_OFFDEN_DPACK_WORK_OUT' &
     &          ,INT(N1,KIND=8),INT(N2,KIND=8),INT(COUNT,KIND=8),0_8 &
     &          ,0.D0,16.D0*REAL(N1,KIND=8)*REAL(N2,KIND=8) &
     &          *REAL(COUNT,KIND=8),ACCEL_T1-ACCEL_T0)
            CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC EXIT DATA DELETE(A(1:N1,1:NBH),B(1:N2*COUNT,1:NBH) &
!$ACC&                 ,WORK(1:N1,1:N2*COUNT),J0ARR(1:COUNT) &
!$ACC&                 ,EIKRARR(1:COUNT))
          END IF
          IF(TDEVICEACCUMACTIVE) DEALLOCATE(MOFFARR)
          DEALLOCATE(EIKRARR)
          DEALLOCATE(J0ARR)
        ELSE
          DO IBH=1,NBH
            DO I=1,N1
              A(I,IBH)=PROJ(IBH,I0+I)
            ENDDO
          ENDDO
          DO K=1,COUNT
            NN=IDX(K)
            IAT2=OSDENMAT(NN)%IAT2
            J0=IPRO1(IAT2)-1
            IOFF=(K-1)*N2
            SVAR=2.D0*PI*SUM(XK(:)*REAL(OSDENMAT(NN)%IT,KIND=8))
            EIKR=EXP(CI*SVAR)
            DO IBH=1,NBH
              FPLUS=0.5D0*(OCC(2*IBH-1)+OCC(2*IBH))
              FMINUS=0.5D0*(OCC(2*IBH-1)-OCC(2*IBH))
              DO J=1,N2
                B(IOFF+J,IBH)=FPLUS*CONJG(PROJ(IBH,J0+J)) &
     &              *CONJG(EIKR)+FMINUS*PROJ(IBH,J0+J)*EIKR
              ENDDO
            ENDDO
          ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          CALL ACCELPROFILE$ADD('PAW_OFFDEN_BATCH_PACK' &
     &        ,INT(N1,KIND=8),INT(N2,KIND=8),INT(NBH,KIND=8) &
     &        ,INT(COUNT,KIND=8),0.D0 &
     &        ,16.D0*REAL(NBH,KIND=8) &
     &        *(REAL(N1,KIND=8)+REAL(N2,KIND=8)*REAL(COUNT,KIND=8)) &
     &        ,ACCEL_T1-ACCEL_T0)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          TCUBLAS_USED=.FALSE.
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          CALL CPPAW_CUBLAS_ACC_ZGEMM_NT_COPY(N1,N2*COUNT,NBH &
     &       ,A,B,WORK,TCUBLAS_USED)
#ENDIF
          IF(.NOT.TCUBLAS_USED) THEN
            CALL ZGEMM('N','T',N1,N2*COUNT,NBH,ONE,A,N1 &
     &                ,B,N2*COUNT,ZERO,WORK,N1)
          END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          ACCEL_FLOPS=OFFDENFLOPS
          IF(TCUBLAS_USED) THEN
            ACCEL_GEMM_NAME='CUBLAS_ZGEMM_OFFDEN_TINV_STACK'
          ELSE
            ACCEL_GEMM_NAME='ZGEMM_OFFDEN_TINV_STACK'
          END IF
          CALL ACCELPROFILE$ADD(ACCEL_GEMM_NAME &
     &        ,INT(N1,KIND=8),INT(N2,KIND=8),INT(NBH,KIND=8) &
     &        ,INT(COUNT,KIND=8),ACCEL_FLOPS &
     &        ,16.D0*(REAL(NBH,KIND=8)*REAL(N1+N2,KIND=8) &
     &        +REAL(N1,KIND=8)*REAL(N2,KIND=8))*REAL(COUNT,KIND=8) &
     &        ,ACCEL_T1-ACCEL_T0)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
        END IF
        IF(.NOT.TDEVICEACCUMACTIVE) THEN
          DO K=1,COUNT
            NN=IDX(K)
            IOFF=(K-1)*N2
            IF(NSPIN.EQ.1) THEN
              DO J=1,N2
                DO I=1,N1
                  OSDENMAT(NN)%MAT(I,J,1)=OSDENMAT(NN)%MAT(I,J,1) &
     &                +REAL(WORK(I,IOFF+J),KIND=8)
                ENDDO
              ENDDO
            ELSE IF(NSPIN.EQ.2) THEN
              IF(ISPIN.EQ.1) THEN
                DO J=1,N2
                  DO I=1,N1
                    OSDENMAT(NN)%MAT(I,J,1)=OSDENMAT(NN)%MAT(I,J,1) &
     &                  +REAL(WORK(I,IOFF+J),KIND=8)
                    OSDENMAT(NN)%MAT(I,J,2)=OSDENMAT(NN)%MAT(I,J,2) &
     &                  +REAL(WORK(I,IOFF+J),KIND=8)
                  ENDDO
                ENDDO
              ELSE
                DO J=1,N2
                  DO I=1,N1
                    OSDENMAT(NN)%MAT(I,J,1)=OSDENMAT(NN)%MAT(I,J,1) &
     &                  +REAL(WORK(I,IOFF+J),KIND=8)
                    OSDENMAT(NN)%MAT(I,J,2)=OSDENMAT(NN)%MAT(I,J,2) &
     &                  -REAL(WORK(I,IOFF+J),KIND=8)
                  ENDDO
                ENDDO
              END IF
            END IF
          ENDDO
        END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T1)
        CALL ACCELPROFILE$ADD('PAW_OFFDEN_BATCH_ACCUM' &
     &      ,INT(N1,KIND=8),INT(N2,KIND=8),INT(NDIMD,KIND=8) &
     &      ,INT(COUNT,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
        DEALLOCATE(WORK)
        DEALLOCATE(B)
        DEALLOCATE(A)
      ENDDO
      IF(TDEVICEFLAT) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC UPDATE SELF(MATFLAT(1:NTOT))
!$ACC WAIT
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T1)
        CALL ACCELPROFILE$ADD('ACC_COPY_OFFDEN_DPACK_FLAT_OUT' &
     &      ,INT(NTOT,KIND=8),INT(NDIMD,KIND=8),0_8,0_8,0.D0 &
     &      ,8.D0*REAL(NTOT,KIND=8),ACCEL_T1-ACCEL_T0)
        CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
!$ACC EXIT DATA DELETE(MATFLAT(1:NTOT))
        DO NN=1,NND
          IF(MATOFF(NN).EQ.0) CYCLE
          IAT1=OSDENMAT(NN)%IAT1
          IAT2=OSDENMAT(NN)%IAT2
          N1=NPROAT(IAT1)
          N2=NPROAT(IAT2)
          DO JDIM=1,NDIMD
            DO J=1,N2
              DO I=1,N1
                MOFF=MATOFF(NN)+(JDIM-1)*N1*N2+(J-1)*N1+I-1
                OSDENMAT(NN)%MAT(I,J,JDIM) &
     &            =OSDENMAT(NN)%MAT(I,J,JDIM)+MATFLAT(MOFF)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T1)
        CALL ACCELPROFILE$ADD('PAW_OFFDEN_FLAT_ACCUM_SCATTER' &
     &      ,INT(NTOT,KIND=8),INT(NND,KIND=8),INT(NDIMD,KIND=8),0_8 &
     &      ,0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
        DEALLOCATE(MATFLAT)
        DEALLOCATE(MATOFF)
      END IF
!$ACC END DATA
      DEALLOCATE(IDX)
      DEALLOCATE(DONE)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$OFFSITEHAMIL()
!     **************************************************************************
!     **  CONSTRUCT HPROJ = H<P|PSI> = DE/D<PSI|P> FROM THE OFF-SITE          **
!     **  HAMILTONIAN AND THE PROJECTIONS                                     **
!     **                                                                      **
!     **  THE OFF-SITE HAMILTONIAN IS IN THE SPIN REPRESENTATION              **
!     **    NSPIN=1,NDIM=1: TOTAL                        (NDIMD=1)            **
!     **    NSPIN=2,NDIM=1: TOTAL,SPIN_Z                 (NDIMD=2)            **
!     **    NSPIN=1,NDIM=2: TOTAL,SPIN_X,SPIN_Y,SPIN_Z   (NDIMD=3)            **
!     **                                                                      **
!     **  OSHAMIL=DE/DRHODAGGER                                               **
!     **************************************************************************
      USE WAVES_MODULE, ONLY: WVSET_TYPE &
     &                       ,NKPTL &
     &                       ,NSPIN &
     &                       ,NDIM &
     &                       ,NDIMD &
     &                       ,THIS &
     &                       ,MAP &
     &                       ,WAVES_SELECTWV &
     &                       ,GSET &
     &                       ,OSHAMIL
      USE RSPACEOP_MODULE, ONLY: RSPACEMAT_TYPE &
     &                          ,RSPACEOP$COPY &
     &                          ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)             :: NAT
      INTEGER(4)             :: N1,N2
      INTEGER(4)             :: NND
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: NB,NBH
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      INTEGER(4),ALLOCATABLE :: IPRO1(:)
      INTEGER(4),ALLOCATABLE :: NPROAT(:)
      INTEGER(4)             :: IAT,NN,ISP,IPRO,IKPT,ISPIN,I,J,IBH,IB
      REAL(8)                :: SVAR
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IAT1,IAT2,IT(3),I0,J0,IDIM,JDIM
      COMPLEX(8)             :: EIKR,C1,C2,CSVAR22(NDIM,NDIM)
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
!     **************************************************************************
      IF(.NOT.ALLOCATED(OSHAMIL)) RETURN
                                 CALL TRACE$PUSH('WAVES$OFFSITEHAMIL')
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
!
!     ==========================================================================
!     ==  CONSTRUCT INDEX ARRAYS                                              ==
!     ==========================================================================
      NAT=MAP%NAT
      ALLOCATE(IPRO1(NAT))
      ALLOCATE(NPROAT(NAT))
      IPRO=1
      DO IAT=1,NAT
        ISP=MAP%ISP(IAT)
        IPRO1(IAT)=IPRO
        NPROAT(IAT)=MAP%LMNX(ISP)
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
      NPRO=SUM(NPROAT(:))
!
!     ==========================================================================
!     ==  ADD UP DENSITY MATRIX                                               ==
!     ==========================================================================
      NND=SIZE(OSHAMIL)
      NPRO=MAP%NPRO
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          NBH=THIS%NBH
          NB=THIS%NB
          IF(.NOT.ASSOCIATED(THIS%HPROJ)) ALLOCATE(THIS%HPROJ(NDIM,NBH,NPRO))
          THIS%HPROJ(:,:,:)=(0.D0,0.D0)
          DO NN=1,NND
            IAT1=OSHAMIL(NN)%IAT1
            IAT2=OSHAMIL(NN)%IAT2
            IT=OSHAMIL(NN)%IT
            SVAR=2.D0*PI*SUM(XK(:,IKPT)*REAL(IT))
            EIKR=EXP(-CI*SVAR)  !<P_{R+T}|PSI>=<P_R|PSI>*EIKR
            I0=IPRO1(IAT1)-1
            J0=IPRO1(IAT2)-1
            DO I=1,NPROAT(IAT1)
              DO J=1,NPROAT(IAT2)
!
!               == CONVERT FROM TOTAL/SPIN INTO UP-DOWN REPRESENTATION =====
                IF(NSPIN.EQ.1) THEN
                  IF(NDIM.EQ.1) THEN !NON-SPIN-POLARIZED
                    CSVAR22(1,1)=CMPLX(OSHAMIL(NN)%MAT(I,J,1),0.D0,KIND=8)
                  ELSE ! NONCOLLINEAR
                    CSVAR22(1,1)=CMPLX(OSHAMIL(NN)%MAT(I,J,1) &
     &                                +OSHAMIL(NN)%MAT(I,J,4),0.D0,KIND=8)
                    CSVAR22(1,2)=CMPLX(OSHAMIL(NN)%MAT(I,J,2) &
     &                               ,-OSHAMIL(NN)%MAT(I,J,3),KIND=8)
                    CSVAR22(2,1)=CMPLX(OSHAMIL(NN)%MAT(I,J,2) &
     &                               ,+OSHAMIL(NN)%MAT(I,J,3),KIND=8)
                    CSVAR22(2,2)=CMPLX(OSHAMIL(NN)%MAT(I,J,1) &
     &                                -OSHAMIL(NN)%MAT(I,J,4),0.D0,KIND=8)
                  END IF
                ELSE IF(NSPIN.EQ.2) THEN !SPIN POLARIZED
                  IF(ISPIN.EQ.1) THEN
                    CSVAR22(1,1)=CMPLX(OSHAMIL(NN)%MAT(I,J,1) &
     &                                +OSHAMIL(NN)%MAT(I,J,2),0.D0,KIND=8)
                  ELSE
                    CSVAR22(1,1)=CMPLX(OSHAMIL(NN)%MAT(I,J,1) &
     &                                -OSHAMIL(NN)%MAT(I,J,2),0.D0,KIND=8)
                  END IF
                END IF
                CSVAR22(:,:)=CSVAR22(:,:)*EIKR
!
                DO IB=1,NBH
                  DO IDIM=1,NDIM
                    DO JDIM=1,NDIM
!===============================================================================
!       E=E_0 + SUM_{A,B} H_{A,B}D_{A,B} +... = E_0 + TR[ H * D-DAGGER ] +... ==
!                                                                             ==
!         => DE=<DPSI_N| [SUM_{A,B} |P_B>H_{A,B}<P_A|PSI_N> F_N + ...         ==
!===============================================================================
                      THIS%HPROJ(JDIM,IB,J0+J)=THIS%HPROJ(JDIM,IB,J0+J) &
      &                           +CSVAR22(JDIM,IDIM)*THIS%PROJ(IDIM,IB,I0+I)
!
! THIS IS THE OLD VERSION (THIS IS CORRECT: SEE METHODS SECTION 'SECOND QUANT..'
!                      THIS%HPROJ(JDIM,IB,J0+J)=THIS%HPROJ(JDIM,IB,J0+J) &
!      &                           +CSVAR22(JDIM,IDIM)*THIS%PROJ(IDIM,IB,I0+I)
! THIS SHOULD BE CORRECT.(NO!)
!!$                      THIS%HPROJ(IDIM,IB,I0+I)=THIS%HPROJ(IDIM,IB,I0+I) &
!!$      &                           +CSVAR22(IDIM,JDIM)*THIS%PROJ(JDIM,IB,J0+J)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == PRINT FOR TESTING                                                    ==
!     ==========================================================================
      IF(TPR) THEN
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            NBH=THIS%NBH
            WRITE(*,FMT='(80("="),T10," HPROJ FOR XK=",3F10.5,"  ")') &
     &                                                                XK(:,IKPT)
            DO IBH=1,NBH
              DO IDIM=1,NDIM
                WRITE(*,FMT='(80("-"),T10," IBH=",I5," IDIM=",I2,"  ")')IBH,IDIM
                WRITE(*,FMT='("RE:",10F20.5)')REAL(THIS%HPROJ(IDIM,IBH,:))
                WRITE(*,FMT='("IM:",10F20.5)')AIMAG(THIS%HPROJ(IDIM,IBH,:))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CALL ERROR$MSG('FORCED STOP AFTER WRITING DIAGNOSTIC RESULTS')
        CALL ERROR$STOP('WAVES$OFFSITEHAMIL')
      END IF
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$RHO(NRL,NDIMD_,RHO,TRHOKIN,RHOKIN)
!     **************************************************************************
!     **  EVALUATES PSEUDO-DENSITY FROM THE ACTUAL PSEUDO WAVE                **
!     **  FUNCTIONS                                                           **
!     **                                                                      **
!     **  RHO1 IS MADE ALLOCATABLE, BECAUSE IT IS TOO LARGE FOR THE STACK     **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NRL             !#(LOCAL REAL SPACE GRID POINTS)
      INTEGER(4),INTENT(IN)  :: NDIMD_          !#(SPIN COMPONENTS IN DENSITY)
      REAL(8)   ,INTENT(OUT) :: RHO(NRL,NDIMD_) !ELECTRON DENSITY
      LOGICAL(4),INTENT(IN)  :: TRHOKIN         !KINETIC ENERGY DENSITY REQUIRED
      REAL(8)   ,INTENT(OUT) :: RHOKIN(NRL,NDIMD_) !KINETIC ENERGY DENSITY
      REAL(8)   ,ALLOCATABLE :: RHO1(:,:)    ! (NRL,NDIM**2)
      REAL(8)   ,ALLOCATABLE :: RHOKIN1(:,:) ! (NRL,NDIM**2)
      INTEGER(4)             :: IKPT,ISPIN,IR
      INTEGER(4)             :: NBX        
      REAL(8)  ,ALLOCATABLE  :: OCC(:,:,:) 
      REAL(8)                :: SVAR1,SVAR2
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
      REAL(8)                :: ACCEL_TOTAL_T0
      REAL(8)                :: ACCEL_TOTAL_T1
      REAL(8)                :: ACCEL_DENSITY_T
      REAL(8)                :: ACCEL_COMBINE_T
      REAL(8)                :: ACCEL_SPIN_T
#ENDIF
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$RHO')
                              CALL TIMING$CLOCKON('W:RHO')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      ACCEL_DENSITY_T=0.D0
      ACCEL_COMBINE_T=0.D0
      ACCEL_SPIN_T=0.D0
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T0)
#ENDIF
      IF(NDIMD_.NE.NDIMD.OR.NRL.NE.MAP%NRL) THEN
        CALL ERROR$MSG('ARRAY SIZE INCONSISTEN')
        CALL ERROR$STOP('WAVES$RHO')
      END IF
!
!     ==========================================================================
!     ==  GET OCCUPATIONS FROM DYNOCC OBJECT                                  ==
!     ==========================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
!
!     ==========================================================================
!     ==  CALCULATE DENSITY                                                   ==
!     ==========================================================================
      ALLOCATE(RHO1(NRL,NDIM**2))
      ALLOCATE(RHOKIN1(NRL,NDIM**2)) 
      RHO(:,:)=0.D0
      RHOKIN(:,:)=0.D0
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          CALL WAVES_DENSITY(GSET%NGL,MAP%NRL,NDIM,THIS%NB,THIS%NBH &
     &             ,OCC(:THIS%NB,IKPT,ISPIN),THIS%PSI0,RHO1,TRHOKIN,RHOKIN1)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          ACCEL_DENSITY_T=ACCEL_DENSITY_T+ACCEL_T1-ACCEL_T0
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          IF(NDIM.EQ.1) THEN
            RHO(:,ISPIN)=RHO(:,ISPIN)+RHO1(:,1)
            IF(TRHOKIN) RHOKIN(:,ISPIN)=RHOKIN(:,ISPIN)+RHOKIN1(:,1)
          ELSE
            RHO=RHO+RHO1
            IF(TRHOKIN) RHOKIN=RHOKIN+RHOKIN1
          END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          ACCEL_COMBINE_T=ACCEL_COMBINE_T+ACCEL_T1-ACCEL_T0
#ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(RHO1)
      DEALLOCATE(RHOKIN1)  
!
!     ==========================================================================
!     == CONVERT INTO TOTAL AND SPIN DENSITY                                  ==
!     ==========================================================================
      IF(NSPIN.EQ.2) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
        DO IR=1,NRL
          SVAR1=RHO(IR,1)
          SVAR2=RHO(IR,2)
          RHO(IR,1)=SVAR1+SVAR2   ! TOTAL DENSITY
          RHO(IR,2)=SVAR1-SVAR2   ! SPIN DENSITY
        ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T1)
        ACCEL_SPIN_T=ACCEL_T1-ACCEL_T0
#ENDIF
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$ADD('PAW_RHO_DENSITY' &
     &    ,INT(NRL,KIND=8),INT(NDIMD_,KIND=8),INT(NKPTL,KIND=8) &
     &    ,INT(NSPIN,KIND=8),0.D0,0.D0,ACCEL_DENSITY_T)
      CALL ACCELPROFILE$ADD('PAW_RHO_COMBINE' &
     &    ,INT(NRL,KIND=8),INT(NDIMD_,KIND=8),INT(NKPTL,KIND=8) &
     &    ,INT(NSPIN,KIND=8),0.D0,0.D0,ACCEL_COMBINE_T)
      IF(NSPIN.EQ.2) THEN
        CALL ACCELPROFILE$ADD('PAW_RHO_SPIN' &
     &      ,INT(NRL,KIND=8),INT(NDIMD_,KIND=8),INT(NKPTL,KIND=8) &
     &      ,INT(NSPIN,KIND=8),0.D0,0.D0,ACCEL_SPIN_T)
      END IF
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T1)
      CALL ACCELPROFILE$ADD('PAW_RHO_TOTAL' &
     &    ,INT(NRL,KIND=8),INT(NDIMD_,KIND=8),INT(NKPTL,KIND=8) &
     &    ,INT(NSPIN,KIND=8),0.D0,0.D0,ACCEL_TOTAL_T1-ACCEL_TOTAL_T0)
#ENDIF
                              CALL TIMING$CLOCKOFF('W:RHO')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$RHO
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$MOMENTS(NAT,LMRXX,NDIMD_,LMNXX,DENMAT,QLM)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES MULTIPOLE MOMENTS OF THE ONE-CENTER DENSITY               **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: LMRXX
      INTEGER(4),INTENT(IN)  :: NAT
      INTEGER(4),INTENT(IN)  :: LMNXX
      INTEGER(4),INTENT(IN)  :: NDIMD_
      COMPLEX(8),INTENT(IN)  :: DENMAT(LMNXX,LMNXX,NDIMD,NAT)
      REAL(8)   ,INTENT(OUT) :: QLM(LMRXX,NAT)
      INTEGER(4)             :: ISP
      INTEGER(4)             :: IAT
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: LMRX
      INTEGER(4)             :: NTASKS,THISTASK
      COMPLEX(8),ALLOCATABLE :: DENMAT1(:,:)
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$MOMENTS')
                              CALL TIMING$CLOCKON('W:MOMENTS')
      IF(NDIMD_.NE.NDIMD.OR.NAT.NE.MAP%NAT) THEN
        CALL ERROR$MSG('ARRAY SIZE INCONSISTEN')
        CALL ERROR$STOP('WAVES$RHO')
      END IF
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
!     == DENMAT IS IDENTICAL AND COMPLETE ON EACH PROCESS OF A MONOMER
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      QLM(:,:)=0.D0
      DO IAT=THISTASK,NAT,NTASKS   ! DISTRIBUTE WORK ACCROSS TASKS
        ISP=MAP%ISP(IAT)
        CALL SETUP$ISELECT(ISP)
        LMNX=MAP%LMNX(ISP)
        CALL SETUP$GETI4('LMRX',LMRX) !FORMER CALL SETUP$LMRX(ISP,LMRX)
        CALL SETUP$ISELECT(0)
        ALLOCATE(DENMAT1(LMNX,LMNX))
        DENMAT1(:,:)=DENMAT(1:LMNX,1:LMNX,1,IAT)
!       ==  SETUP MUST BE UNSELECTED BEFORE CALLING AUGMENTATION$MOMENTS
        CALL AUGMENTATION$MOMENTS(ISP,LMNX,DENMAT1,LMRX,QLM(1,IAT))
        DEALLOCATE(DENMAT1)
      ENDDO
      CALL MPE$COMBINE('MONOMER','+',QLM)
!PI=4.D0*ATAN(1.D0)
!SVAR1=0.D0
!DO IAT=1,NAT
!  SVAR1=SVAR1+QLM(1,IAT)*SQRT(4.D0*PI)
!ENDDO
!PRINT*,'TOTAL CHARGE IN AUGMENTATION',SVAR1
                              CALL TIMING$CLOCKOFF('W:MOMENTS')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$MOMENTS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$SPINS(NRL,NDIMD,RHO,NAT,LMNXX,DENMAT)
!      USE WAVES_MODULE, ONLY : MAP 
      USE MPE_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NRL
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: LMNXX
      REAL(8)   ,INTENT(IN) :: RHO(NRL,NDIMD)      
      COMPLEX(8),INTENT(IN) :: DENMAT(LMNXX,LMNXX,NDIMD,NAT)
      INTEGER(4)            :: NTASKS,THISTASK
      REAL(8)               :: PSSPIN(3)
      REAL(8)               :: AUGSPIN(3)
      REAL(8)               :: TOTSPIN(3)
      INTEGER(4)            :: GID     ! GRID ID FOR RADIAL GRID
      INTEGER(4)            :: NR      ! #(RADIAL GRID POINTS)
      REAL(8)   ,ALLOCATABLE:: R(:)    ! RADIAL GRID
      REAL(8)               :: RBAS(3,3),GBAS(3,3),CELLVOL
      INTEGER(4)            :: IAT
      REAL(8)               :: AEZ     ! ATOMIC NUMBER
      INTEGER(4)            :: ISP
      INTEGER(4)            :: LNX
      INTEGER(4),ALLOCATABLE:: LOX(:)
      REAL(8)   ,ALLOCATABLE:: DOVER(:,:)
      INTEGER(4)            :: L1,L2,LN1,LN2,LMN1,LMN2,M,IDIMD
      REAL(8)               :: CM(NDIMD,NAT)
      REAL(8)               :: CM1(NDIMD)   ! CHARGE AND MOMENT
      REAL(8)   ,ALLOCATABLE:: WORK(:), WORK1(:)
      REAL(8)   ,ALLOCATABLE:: AEPHI(:,:)
      REAL(8)               :: RAD  !APPROXIMATE ASA RADIUS OF THE ATOM
      INTEGER(4)            :: NFILO
!     **************************************************************************
                            CALL TRACE$PUSH('WAVES$SPINS')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     ==========================================================================
!     ==  DETERMINE SPINS ON THE INDIVIDUAL ATOMS                             ==
!     ==========================================================================
      CM(:,:)=0.D0
      DO IAT=1,NAT
        IF(MOD(IAT-1,NTASKS)+1.NE.THISTASK)  CYCLE
!PRINT*,'MARKE 1',THISTASK,IAT
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX) !FORMER CALL SETUP$LNX(ISP,LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX) !FORMER CALL SETUP$LOFLN(ISP,LNX,LOX)
!CONSISTENCY CHECK IN SEARCH OF A BUG
IF(SUM(2*LOX+1).GT.LMNXX) THEN
  CALL ERROR$MSG('LMNX EXCEEDS LMNXX')
  CALL ERROR$I4VAL('IAT',IAT)
  CALL ERROR$I4VAL('LNX',LNX)
  CALL ERROR$I4VAL('LMNXX',LMNXX)
  CALL ERROR$STOP('WAVES$SPINS')
END IF
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(AEZ,'R(ASA)',RAD)
!PRINT*,'MARKE 2',THISTASK,IAT
!
!       ========================================================================
!       ==  DETERMINE OVERLAP OF PARTIALWAVES WITHIN ASA SPHERE               ==
!       ========================================================================
        CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        ALLOCATE(AEPHI(NR,LNX))
        CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
                             !FORMER CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
!PRINT*,'MARKE 3',THISTASK,IAT
        ALLOCATE(DOVER(LNX,LNX))
        ALLOCATE(WORK(NR))
        ALLOCATE(WORK1(NR))
        DOVER(:,:)=0.D0
!PRINT*,'MARKE 4',THISTASK,IAT
        DO LN1=1,LNX
          DO LN2=LN1,LNX
            IF(LOX(LN1).NE.LOX(LN2)) CYCLE
            WORK(:)=AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
            CALL RADIAL$INTEGRATE(GID,NR,WORK,WORK1)
            CALL RADIAL$VALUE(GID,NR,WORK1,RAD,DOVER(LN1,LN2))
            DOVER(LN2,LN1)=DOVER(LN1,LN2)
          ENDDO
        ENDDO
!PRINT*,'MARKE 5',THISTASK,IAT
        DEALLOCATE(WORK1)
        DEALLOCATE(WORK)
        DEALLOCATE(AEPHI)
        DEALLOCATE(R)
!PRINT*,'MARKE 6',THISTASK,IAT
!
!       ========================================================================
!       ==  DETERMINE CHARGE AND SPIN PROJECTIONS                             ==
!       ========================================================================
        CM1(:)=0.D0
        LMN1=0
        DO LN1=1,LNX
          L1=LOX(LN1)
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            IF(L1.EQ.L2) THEN
              DO M=1,2*L1+1
                CM1(:)=CM1(:)+REAL(DENMAT(LMN1+M,LMN2+M,:,IAT))*DOVER(LN1,LN2)
              ENDDO
            END IF
            LMN2=LMN2+2*L2+1
          ENDDO
          LMN1=LMN1+2*L1+1
        ENDDO
!PRINT*,'MARKE 7',THISTASK,IAT
        CM(:,IAT)=CM(:,IAT)+CM1(:)
        DEALLOCATE(DOVER)
        DEALLOCATE(LOX)
!!$DO IDIMD=1,NDIMD
!!$  PRINT*,'DENMAT FOR ATOM ',AEZ,' AND IDIMD=',IDIMD
!!$  DO LMN1=1,SUM(2*LOX(:)+1)
!!$    WRITE(*,FMT='(15F10.3)')REAL(DENMAT(LMN1,:,IDIMD,IAT))
!!$  ENDDO
!!$ENDDO
        CALL SETUP$ISELECT(0)
!PRINT*,'MARKE 8',THISTASK,IAT
      ENDDO
      CALL MPE$COMBINE('MONOMER','+',CM)
      CALL ATOMLIST$SETR8A('CHARGEANDMOMENTS',0,4*NAT,CM)
!!$PRINT*, 'CHARGE AND SPINS '
!!$DO IAT=1,NAT
!!$  CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
!!$  CALL SETUP$ISELECT(ISP)
!!$  CALL SETUP$GETR8('AEZ',AEZ)
!!$  WRITE(*,*)IAT,AEZ,CM(:,IAT)
!!$ENDDO
                                   CALL TRACE$POP()
RETURN

!
!     ==========================================================================
!     ==  DETERMINE TOTAL MAGNETIZATION AS INTEGRATED MOMENT DENSITY          ==
!     ==========================================================================
      PSSPIN(:)=0.D0
      IF(NDIMD.EQ.2) THEN
        PSSPIN(3)=SUM(RHO(:,2))
      ELSE IF(NDIMD.EQ.4) THEN
        PSSPIN(1)=SUM(RHO(:,2))
        PSSPIN(2)=SUM(RHO(:,3))
        PSSPIN(3)=SUM(RHO(:,4))
      END IF
      CALL MPE$COMBINE('MONOMER','+',PSSPIN)
      NR=NRL
      CALL MPE$COMBINE('MONOMER','+',NR)
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      PSSPIN(:)=PSSPIN(:)*CELLVOL/REAL(NR,KIND=8)
!
!     ==========================================================================
!     ==  DETERMINE AUGMENTATION PART OF THE SPIN DENSITY                     ==
!     ==========================================================================
      AUGSPIN(:)=0.D0
      DO IAT=1,NAT
        IF(MOD(IAT-1,NTASKS)+1.NE.THISTASK)  CYCLE
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX) !FORMER CALL SETUP$LNX(ISP,LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX) !FORMER CALL SETUP$LOFLN(ISP,LNX,LOX)
        ALLOCATE(DOVER(LNX,LNX))
        CALL SETUP$GETR8A('DO',LNX*LNX,DOVER)
                                     !FORMER CALL SETUP$1COVERLAP(ISP,LNX,DOVER)
        LMN1=0
        DO LN1=1,LNX
          L1=LOX(LN1)
          LMN2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            IF(L1.EQ.L2) THEN
              DO M=1,2*L1+1
                IF(NDIMD.EQ.2) THEN
                  AUGSPIN(3)=AUGSPIN(3)+REAL(DENMAT(LMN1+M,LMN2+M,2,IAT))*DOVER(LN1,LN2)
                ELSE IF(NDIMD.EQ.4) THEN
                  AUGSPIN(1)=AUGSPIN(1)+REAL(DENMAT(LMN1+M,LMN2+M,2,IAT))*DOVER(LN1,LN2)
                  AUGSPIN(2)=AUGSPIN(2)+REAL(DENMAT(LMN1+M,LMN2+M,3,IAT))*DOVER(LN1,LN2)
                  AUGSPIN(3)=AUGSPIN(3)+REAL(DENMAT(LMN1+M,LMN2+M,4,IAT))*DOVER(LN1,LN2)
                END IF
              ENDDO
            END IF
            LMN2=LMN2+2*L2+1
          ENDDO
          LMN1=LMN1+2*L1+1
        ENDDO
        DEALLOCATE(LOX)
        DEALLOCATE(DOVER)
      ENDDO
      CALL SETUP$ISELECT(0)
      CALL MPE$COMBINE('MONOMER','+',AUGSPIN)
!
!     ==========================================================================
!     ==  REPORT TOTAL SPIN                                                   ==
!     ==========================================================================
      TOTSPIN=PSSPIN+AUGSPIN
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,FMT='("ABS.VALUE OF INTEGRATED SPIN MOMENT DENSITY",F10.5)') &
     &                 SQRT(SUM(TOTSPIN(:)**2))
      WRITE(NFILO,FMT='("DIRECTION OF INTEGRATED SPIN MOMENT DENSITY",3F10.5)')&
     &                 TOTSPIN/SQRT(SUM(TOTSPIN(:)**2))
                                   CALL TRACE$POP()
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$SPHERE(LMNXX,NDIMD_,NAT,LMRXX,RHOB,DENMAT,EDENMAT &
     &                       ,VQLM,DH,DO,POTB)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES ONE-CENTER HAMILTONIAN AN OVERLAPMATRIX                   **
!     **                                                                      **
!     **  ONE-CENTER DENSITY MATRICES                                         **
!     **  SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                             **
!     **  SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)                      **
!     **  NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z)        **
!     **                                                                      **
!     **  ON RETURN, POTB IS MINUS THE AVERAGED ELECTROSTATIC                 **
!     **  POTENTIAL FROM THE AUGMENTATION,                                    **
!     **  I.E. -(AEPOTHARTREE-(PSPOTHARTREE-VADD))                            **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: LMNXX
      INTEGER(4),INTENT(IN)  :: NDIMD_
      INTEGER(4),INTENT(IN)  :: NAT    
      INTEGER(4),INTENT(IN)  :: LMRXX
      REAL(8)   ,INTENT(IN)  :: RHOB
      COMPLEX(8),INTENT(IN)  :: DENMAT(LMNXX,LMNXX,NDIMD_,NAT)
      COMPLEX(8),INTENT(IN)  :: EDENMAT(LMNXX,LMNXX,NDIMD,NAT)
      REAL(8)   ,INTENT(INOUT):: VQLM(LMRXX,NAT)
      COMPLEX(8),INTENT(OUT) :: DH(LMNXX,LMNXX,NDIMD_,NAT)
      REAL(8)   ,INTENT(OUT) :: DO(LMNXX,LMNXX,NDIMD_,NAT)
      REAL(8)   ,INTENT(OUT) :: POTB 
      REAL(8)                :: POTB1
      INTEGER(4)             :: ISP   ! SPECIES INDEX
      INTEGER(4)             :: IAT   ! 
      INTEGER(4)             :: LMNX  ! 
      INTEGER(4)             :: LMRX  ! 
      INTEGER(4)             :: NTASKS,THISTASK
      COMPLEX(8),ALLOCATABLE :: DENMAT1(:,:,:)
      COMPLEX(8),ALLOCATABLE :: EDENMAT1(:,:,:)
      COMPLEX(8),ALLOCATABLE :: DH1(:,:,:)  !(LMNX,LMNX,NDIMD)
      REAL(8)   ,ALLOCATABLE :: DOV1(:,:,:) !(LMNX,LMNX,NDIMD)
      REAL(8)                :: RBAS(3,3)
      REAL(8)                :: GBAS(3,3)
      REAL(8)                :: SVAR
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
#ENDIF
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$SPHERE')
                              CALL TIMING$CLOCKON('W:SPHERE')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      IF(NDIMD_.NE.NDIMD.OR.NAT.NE.MAP%NAT) THEN
        CALL ERROR$MSG('ARRAY SIZE INCONSISTENT')
        CALL ERROR$STOP('WAVES$DENMAT')
      END IF
!
!     ==========================================================================
!     ==  EVALUATE ONE-CENTER HAMILTONIAN                                     ==
!     ==========================================================================
      DH(:,:,:,:)=(0.D0,0.D0)
      DO(:,:,:,:)=0.D0
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      POTB=0
      DO IAT=THISTASK,NAT,NTASKS   ! DISTRIBUTE WORK ACCROSS TASKS
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LMRX',LMRX) !FORMER CALL SETUP$LMRX(ISP,LMRX)
        CALL SETUP$UNSELECT()
        ALLOCATE(DENMAT1(LMNX,LMNX,NDIMD))
        ALLOCATE(EDENMAT1(LMNX,LMNX,NDIMD))
        DENMAT1(:,:,:)=DENMAT(1:LMNX,1:LMNX,:,IAT)
        EDENMAT1(:,:,:)=EDENMAT(1:LMNX,1:LMNX,:,IAT)
        ALLOCATE(DH1(LMNX,LMNX,NDIMD))
        ALLOCATE(DOV1(LMNX,LMNX,NDIMD))
        CALL AUGMENTATION$SPHERE(ISP,IAT,LMNX,NDIMD,DENMAT1,EDENMAT1 &
     &               ,LMRX,VQLM(1,IAT),RHOB,POTB1,DH1,DOV1)
        POTB=POTB+POTB1
        DH(1:LMNX,1:LMNX,:,IAT)=DH1(:,:,:)
        DO(1:LMNX,1:LMNX,:,IAT)=DOV1(:,:,:)
        DEALLOCATE(DH1)
        DEALLOCATE(DOV1)
        DEALLOCATE(DENMAT1)
        DEALLOCATE(EDENMAT1)
      END DO
      CALL MPE$COMBINE('MONOMER','+',DH)
      CALL MPE$COMBINE('MONOMER','+',DO)
      CALL MPE$COMBINE('MONOMER','+',POTB)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL GBASS(RBAS,GBAS,SVAR)
      POTB=POTB/SVAR
      CALL AUGMENTATION$SYNC
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_SPHERE_TOTAL' &
     &    ,INT(NAT,KIND=8),INT(LMNXX,KIND=8),INT(NDIMD_,KIND=8) &
     &    ,INT(LMRXX,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
                              CALL TIMING$CLOCKOFF('W:SPHERE')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$SPHERE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$ADDCONSTANTPOT(NRL,LMNXX,NDIMD_,NAT,POTB,RHO,DH,DO)
!     **************************************************************************
!     **                                                                      **
!     **  ADDS A CONSTANT POTENTIAL TO ENSURE THAT THE AVERAGED               **
!     **  ELECTROSTATIC POTENTIAL IS ZERO.                                    **
!     **                                                                      **
!     **  POTB IS MINUS THE INTEGRATED ELECTROSTATIC POTENTIAL FROM           **
!     **  THE AUGMENTATION AVERAGED OVER THE UNIT CELL                        **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NRL
      REAL(8)   ,INTENT(IN)   :: POTB
      INTEGER(4),INTENT(IN)   :: LMNXX
      INTEGER(4),INTENT(IN)   :: NDIMD_
      INTEGER(4),INTENT(IN)   :: NAT
      REAL(8)   ,INTENT(IN)   :: DO(LMNXX,LMNXX,NDIMD_,NAT)
      COMPLEX(8),INTENT(INOUT):: DH(LMNXX,LMNXX,NDIMD_,NAT)
      REAL(8)   ,INTENT(INOUT):: RHO(NRL,NDIMD_)
      INTEGER(4)              :: ISP
      INTEGER(4)              :: IAT
      INTEGER(4)              :: LMN1,LMN2
      INTEGER(4)              :: LN1,LN2
      INTEGER(4)              :: L1,L2,M
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$ADDCONSTANTPOT')
                              CALL TIMING$CLOCKON('W:ADDCONSTANTPOT')
      RHO(:,1)=RHO(:,1)+POTB
      DO IAT=1,NAT
        ISP=MAP%ISP(IAT)
        LMN1=0
        DO LN1=1,MAP%LNX(ISP)
          L1=MAP%LOX(LN1,ISP)
          LMN2=0
          DO LN2=1,MAP%LNX(ISP)
            L2=MAP%LOX(LN2,ISP)
            IF(L2.EQ.L1) THEN
              DO M=1,2*L1+1
                DH(LMN1+M,LMN2+M,1,IAT)=DH(LMN1+M,LMN2+M,1,IAT) &
     &                               +CMPLX(POTB*DO(LMN1+M,LMN2+M,1,IAT),KIND=8)
              ENDDO
            END IF
            LMN2=LMN2+2*L2+1
          ENDDO
          LMN1=LMN1+2*L1+1
        ENDDO
      ENDDO
                              CALL TIMING$CLOCKOFF('W:ADDCONSTANTPOT')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$ADDCONSTANTPOT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$FORCE(NAT,LMNXX,NDIMD_,DH,FORCE,STRESS)
!     **************************************************************************
!     **  EVALUATES FORCES FROM THE AUGMENTATION PART                         **
!     **                                                                      **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*****************************
      USE MPE_MODULE
      USE WAVES_MODULE
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &       CPPAW_CUBLAS_ACC_FORCE_PSI_RESIDENCY_ENABLED &
     &      ,CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D
#ENDIF
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NAT
      INTEGER(4),INTENT(IN)  :: LMNXX
      INTEGER(4),INTENT(IN)  :: NDIMD_
      COMPLEX(8),INTENT(IN)  :: DH(LMNXX,LMNXX,NDIMD,NAT)
      REAL(8)   ,INTENT(OUT) :: FORCE(3,NAT)
      REAL(8)   ,INTENT(OUT) :: STRESS(3,3)
      INTEGER(4)             :: IKPT,ISPIN
      LOGICAL(4)             :: TSTRESS     ! CALCLUATE STRESSES OR NOT
      INTEGER(4)             :: NBX
      REAL(8)                :: GWEIGHT    
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)  ! (NBX,NKPTL,NSPIN)    
      LOGICAL(4)             :: TINV        ! SUPER WAVE FUNCTIONS OR NOT
      INTEGER(4)             :: NB          ! #(WAVE FUNCTIONS)
      INTEGER(4)             :: NBH         ! #(SUPER WAVE FUNCTIONS)
      INTEGER(4)             :: NGL         ! LOCAL #(G-VECTORS) 
      INTEGER(4)             :: IJ,I,J,IPRO,ISP,IAT,IBPRO,LMNX,LNX,IG
      REAL(8)   ,ALLOCATABLE :: GVEC(:,:)   ! (3,NGL)    
      REAL(8)   ,ALLOCATABLE :: GIJ(:,:)    ! (6,NGL)    
      REAL(8)   ,ALLOCATABLE :: DO1(:,:)    ! (LNX,LNX)    
      COMPLEX(8),ALLOCATABLE :: DH1(:,:,:)  ! (LMNX,LMNX,NDIM**2)    
      COMPLEX(8),ALLOCATABLE :: DEDPROJ(:,:,:) ! (NDIM,NBH,LMNX)
      COMPLEX(8),ALLOCATABLE :: DEDPROJ1(:,:,:) ! (NDIM,NBH,LMNX)
      COMPLEX(8),ALLOCATABLE :: DEDPRO(:,:) ! (NGL,LMNX)
      COMPLEX(8),ALLOCATABLE :: EIGR(:)     ! (NGL)
      REAL(8)                :: FORCE1(3)
      REAL(8)                :: STRESS1(3,3)
      REAL(8)                :: SVAR
      REAL(8)                :: R(3,NAT)
      REAL(8)   ,PARAMETER   :: RSMALL=1.D-20
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      LOGICAL(4)             :: TRESIDENTFORCEPSI
#ENDIF
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$FORCE')
                              CALL TIMING$CLOCKON('W:FORCE')
      IF(NDIMD_.NE.NDIMD.OR.NAT.NE.MAP%NAT) THEN
        CALL ERROR$MSG('ARRAY SIZE INCONSISTENT')
        CALL ERROR$STOP('WAVES$FORCE')
      END IF
      IF(NDIMD.NE.NSPIN*NDIM**2) THEN
        CALL ERROR$MSG('VALUE OF NDIMD INCONSISTENT WITH ITS DEFINITION')
        CALL ERROR$STOP('WAVES$FORCE')
      END IF
!
!     ==========================================================================
!     ==  GET OCCUPATIONS FROM DYNOCC OBJECT                                  ==
!     ==========================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      FORCE(:,:)=0.D0
      STRESS(:,:)=0.D0
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          CALL PLANEWAVE$GETR8('GWEIGHT',GWEIGHT)
          NB=THIS%NB
          NBH=THIS%NBH
          NGL=GSET%NGL
!
!         ======================================================================
!         ==                                                                  ==
!         ======================================================================
          ALLOCATE(GVEC(3,NGL))
          CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
          ALLOCATE(GIJ(6,NGL))
          IF(TSTRESS) THEN
            DO IG=1,NGL
              SVAR=GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2
              SVAR=1.D0/(SVAR+RSMALL)
              IJ=0
              DO I=1,3
                DO J=I,3
                  IJ=IJ+1
                  GIJ(IJ,IG)=GVEC(I,IG)*GVEC(J,IG)*SVAR
                ENDDO
              ENDDO
            ENDDO
          ELSE
            GIJ(:,:)=0.D0
          END IF
!
!         ======================================================================
!         ==                                                                  ==
!         ======================================================================
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          TRESIDENTFORCEPSI=CPPAW_CUBLAS_ACC_FORCE_PSI_RESIDENCY_ENABLED()
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          IF(TRESIDENTFORCEPSI) THEN
            CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &          ('ACC_PRESENT_FORCE_PSI0','ACC_COPY_FORCE_PSI0_IN' &
     &          ,NGL,NDIM,NBH,THIS%PSI0)
          END IF
#ENDIF
!$ACC DATA PRESENT_OR_COPYIN(THIS%PSI0(1:NGL,1:NDIM,1:NBH)) &
!$ACC& IF(TRESIDENTFORCEPSI)
#ENDIF
          IPRO=1
          DO IAT=1,NAT
            ISP=MAP%ISP(IAT)
            CALL SETUP$ISELECT(ISP)
            IBPRO=1+SUM(MAP%LNX(1:ISP-1))
            LMNX=MAP%LMNX(ISP)
            LNX=MAP%LNX(ISP)
!           == DEDPROJ= (DH-LAMBDA*DO)<P|PSPSI> ================================
            ALLOCATE(DO1(LNX,LNX))
            CALL SETUP$GETR8A('DO',LNX*LNX,DO1)
                                       !FORMER CALL SETUP$1COVERLAP(ISP,LNX,DO1)
            ALLOCATE(DEDPROJ(NDIM,NBH,LMNX))
            ALLOCATE(DH1(LMNX,LMNX,NDIM**2))
            IF(NDIM.EQ.1) THEN
              DH1(:,:,1)=DH(1:LMNX,1:LMNX,ISPIN,IAT)
            ELSE 
              DH1(:,:,:)=DH(1:LMNX,1:LMNX,:,IAT)
            END IF
!           ==  DEDPROJ=DE/D<PSPSI|PRO>=DH*<PRO|PSPSI>
            CALL WAVES_DEDPROJ(NDIM,NBH,NB,LNX,MAP%LOX(1:LNX,ISP),LMNX &
     &                        ,OCC(:,IKPT,ISPIN) &
     &                        ,THIS%PROJ(:,:,IPRO:IPRO+LMNX-1),DH1,DO1 &
     &                        ,THIS%RLAM0,DEDPROJ)
            DEALLOCATE(DH1)
            DEALLOCATE(DO1)
!
!           == ADD CONTRIBUTION FROM NTBOS  ====================================
            IF(ASSOCIATED(THIS%HPROJ)) THEN
              CALL WAVES_FORCE_ADDHTBC(NDIM,NBH,NB,LMNX,OCC(:,IKPT,ISPIN) &
     &                                ,THIS%HPROJ(:,:,IPRO:IPRO+LMNX-1),DEDPROJ)
            END IF
!
            IF(ASSOCIATED(THIS%HTBC_NEW)) THEN
              CALL WAVES_FORCE_ADDHTBC(NDIM,NBH,NB,LMNX,OCC(:,IKPT,ISPIN) &
     &                             ,THIS%HTBC_NEW(:,:,IPRO:IPRO+LMNX-1),DEDPROJ)
            END IF

!           == |DE/DPRO>=|PSPSI>DEDPROJ ========================================
            ALLOCATE(DEDPRO(NGL,LMNX))
            CALL WAVES_DEDPRO(GSET%TINV,NGL,NDIM,NBH,THIS%PSI0,LMNX &
     &                       ,DEDPROJ,DEDPRO)
            DEALLOCATE(DEDPROJ)
!           == DE= <DPRO|DEDPRO> ===============================================
            ALLOCATE(EIGR(NGL))
            CALL PLANEWAVE$STRUCTUREFACTOR(R(1,IAT),NGL,EIGR)
!           == F=2*RE[ <PSI|DPRO/DR>*DH*<PRO|PSI> ]
            CALL WAVES_PROFORCE(LNX,LMNX,MAP%LOX(1:LNX,ISP),NGL,GWEIGHT,GVEC &
     &           ,GIJ &
     &           ,GSET%PRO(:,IBPRO:IBPRO+LNX-1),GSET%DPRO(:,IBPRO:IBPRO+LNX-1) &
     &                 ,MAP%LMX,GSET%YLM,GSET%SYLM &
     &                 ,EIGR,DEDPRO,FORCE1,TSTRESS,STRESS1)
            DEALLOCATE(EIGR)
            DEALLOCATE(DEDPRO)
            FORCE(:,IAT)=FORCE(:,IAT)+FORCE1(:)
            STRESS(:,:) =STRESS(:,:) +STRESS1(:,:)
            IPRO=IPRO+LMNX
            CALL SETUP$UNSELECT()
          ENDDO ! END OF LOOP OVER IAT
#IF DEFINED(CPPVAR_CUBLAS_ACC)
!$ACC END DATA
#ENDIF
          DEALLOCATE(GIJ)
          DEALLOCATE(GVEC)
        ENDDO  ! END OF LOOP OVER ISPIN
      ENDDO    ! END OF LOOP OVER IKPT
      DEALLOCATE(OCC)
!     == EACH PROCESSOR CALCULATES THE CONTRIBUTION TO THE FORCE FROM THE 
!     == LOCAL WAVE FUNCTION COMPONENTS
      CALL MPE$COMBINE('MONOMER','+',FORCE)
      IF(TSTRESS)CALL MPE$COMBINE('MONOMER','+',STRESS)
!
                              CALL TIMING$CLOCKOFF('W:FORCE')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$FORCE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_FORCE_ADDHTBC(NDIM,NBH,NB,LMNX,OCC,HTBC,DEDPROJ)
!     **************************************************************************
!     ** ADDS THE CONTRIBUTION FROM THE NTB-ORBITALS TO DE/DPROJ              **
!     ** THE OCCUPATIONS ARE INCLUDED AT THIS POINT AND NOT BEFORE            **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NDIM
      INTEGER(4),INTENT(IN)    :: NBH
      INTEGER(4),INTENT(IN)    :: NB
      INTEGER(4),INTENT(IN)    :: LMNX
      REAL(8)   ,INTENT(IN)    :: OCC(NB)
      COMPLEX(8),INTENT(IN)    :: HTBC(NDIM,NBH,LMNX)
      COMPLEX(8),INTENT(INOUT) :: DEDPROJ(NDIM,NBH,LMNX)
      INTEGER(4)               :: IB,IB1,IB2,LMN,IDIM
      LOGICAL(4)               :: TINV
      REAL(8)                  :: F1,F2
      COMPLEX(8)               :: CSVAR
!     **************************************************************************
      CALL PLANEWAVE$GETL4('TINV',TINV)
!
!     ==========================================================================
!     ==  HPROJ = DH<P|PSI>*F                                                 ==
!     ==========================================================================
      IF(TINV) THEN
        DO IB=1,NBH
          IB1=2*IB-1       
          IB2=IB1+1
          F1=0.5D0*(OCC(IB1)+OCC(IB2))
          F2=0.5D0*(OCC(IB1)-OCC(IB2))
          DO LMN=1,LMNX
            DO IDIM=1,NDIM
              CSVAR=HTBC(IDIM,IB,LMN)
              DEDPROJ(IDIM,IB,LMN)=DEDPROJ(IDIM,IB,LMN)+F1*CSVAR+F2*CONJG(CSVAR)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO IB=1,NB
          F1=OCC(IB)
          DO LMN=1,LMNX
            DO IDIM=1,NDIM
              CSVAR=HTBC(IDIM,IB,LMN)
              DEDPROJ(IDIM,IB,LMN)=DEDPROJ(IDIM,IB,LMN)+F1*CSVAR
            ENDDO
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$HPSI(NRL,NDIMD_,NAT,LMNXX,RHO,DH)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES FORCES FROM THE AUGMENTATION PART                         **
!     **                                                                      **
!     ********************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : NDIMD,NDIM,NKPTL,NSPIN,MAP,GSET,THIS &
     &                        ,WAVES_SELECTWV
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &       CPPAW_CUBLAS_ACC_HPSI_RESIDENCY_ENABLED &
     &      ,CPPAW_CUBLAS_ACC_SHOULD_USE_ADDPRODUCT &
     &      ,CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D
#ENDIF
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NRL
      INTEGER(4),INTENT(IN) :: NDIMD_
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RHO(NRL,NDIMD_) ! PS-POTENTIAL
      COMPLEX(8),INTENT(IN) :: DH(LMNXX,LMNXX,NDIMD_,NAT)!ONE-CENTER HAMILTONIAN
      INTEGER(4)            :: IKPT,ISPIN
      INTEGER(4)            :: NGL
      INTEGER(4)            :: NBH
      REAL(8)               :: R(3,NAT)
!     --------------------------------------------------------------------------
      INTEGER(4)             :: IPRO
      INTEGER(4)             :: ISP
      INTEGER(4)             :: IAT
      INTEGER(4)             :: LMNX
      COMPLEX(8),ALLOCATABLE :: DH1(:,:,:)      ! 1CENTER HAMILTONIAN
      COMPLEX(8),ALLOCATABLE :: HPROJ(:,:,:)    ! DH*PROJ
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      LOGICAL(4)             :: THPSIRESIDENT
      REAL(8)                :: ADDPROFLOPS
#ENDIF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
      REAL(8)                :: ACCEL_TOTAL_T0
      REAL(8)                :: ACCEL_TOTAL_T1
#ENDIF
!     **************************************************************************
                              CALL TRACE$PUSH('WAVES$HPSI')
                              CALL TIMING$CLOCKON('W:HPSI')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T0)
#ENDIF
      IF(NDIMD_.NE.NDIMD.OR.NRL.NE.MAP%NRL.OR.NAT.NE.MAP%NAT) THEN
        CALL ERROR$MSG('ARRAY SIZE INCONSISTENT')
        CALL ERROR$I4VAL('NDIMD_',NDIMD_)
        CALL ERROR$I4VAL('NDIMD',NDIMD)
        CALL ERROR$I4VAL('NRL',NRL)
        CALL ERROR$I4VAL('MAP%NRL',MAP%NRL)
        CALL ERROR$I4VAL('NAT',NAT)
        CALL ERROR$I4VAL('MAP%NAT',MAP%NAT)
        CALL ERROR$STOP('WAVES$HPSI')
      END IF
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)

      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)      
!         ======================================================================
!         ==  MULTIPLY WAVE FUNCTION WITH LOCAL POTENTIAL AND                 ==
!         ==  ADD KINETIC ENERGY                                              ==
!         ======================================================================
          NGL=GSET%NGL
          NBH=THIS%NBH
          IF(.NOT.ASSOCIATED(THIS%HPSI))ALLOCATE(THIS%HPSI(NGL,NDIM,NBH))
          THIS%HPSI_ACC_RESIDENT=.FALSE.
!         == NOTE: THE ARRAY RHO CONTAINS THE POTENTIAL ========================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
	CALL TIMING$CLOCKON('W:HPSI.VPSI')
!===============================================================================
	IF(1.EQ.1)THEN !OLD VERSION CHANGE WAS REQUIRED FOR KAESTNERS CONJUGATE GRADIENT
!===============================================================================
	         CALL WAVES_VPSI(GSET,NGL,NDIM,NBH,NRL,THIS%PSI0,RHO(1,ISPIN) &
	     &                   ,THIS%HPSI)
	CALL TIMING$CLOCKOFF('W:HPSI.VPSI')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          CALL ACCELPROFILE$ADD('PAW_HPSI_VPSI' &
     &        ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &        ,INT(NRL,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!         ======================================================================
!         ==  EVALUATE  DH<P|PSI>                                             ==
!         ======================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
	CALL TIMING$CLOCKON('W:HPSI.HPROJ')
          ALLOCATE(HPROJ(NDIM,NBH,MAP%NPRO))
          HPROJ(:,:,:)=(0.D0,0.D0)
          IPRO=1
          DO IAT=1,NAT
            ISP=MAP%ISP(IAT)
            LMNX=MAP%LMNX(ISP)
            ALLOCATE(DH1(LMNX,LMNX,NDIM**2))
            IF(NDIM.EQ.1) THEN
              DH1(:,:,1)=DH(1:LMNX,1:LMNX,ISPIN,IAT)
            ELSE 
              DH1(:,:,:)=DH(1:LMNX,1:LMNX,:,IAT)
            END IF
            CALL WAVES_HPROJ(NDIM,NBH,LMNX &
     &         ,DH1,THIS%PROJ(:,:,IPRO:IPRO+LMNX-1),HPROJ(:,:,IPRO:IPRO+LMNX-1))
            DEALLOCATE(DH1)
            IPRO=IPRO+LMNX
          ENDDO
	CALL TIMING$CLOCKOFF('W:HPSI.HPROJ')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          CALL ACCELPROFILE$ADD('PAW_HPSI_HPROJ' &
     &        ,INT(NBH,KIND=8),INT(MAP%NPRO,KIND=8),INT(NAT,KIND=8) &
     &        ,INT(NDIM,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!         ======================================================================
!         ==  ADD POTENTIAL FROM LMTO INTERFACE                               ==
!         ======================================================================
!         == HTBC CONTAINS ALREADY HPROJ IN TERMS OF PARTIAL WAVE 
!         == PROJECTOR FUNCTIONS.
!!$          IF(ASSOCIATED(THIS%HTBC_NEW)) THEN
!!$            HPROJ(:,:,:)=HPROJ(:,:,:)+THIS%HTBC_NEW(:,:,:)
!!$            DEALLOCATE(THIS%HTBC_NEW)
!!$            NULLIFY(THIS%HTBC_NEW)
!!$          END IF
          IF(ASSOCIATED(THIS%HPROJ)) THEN
            HPROJ(:,:,:)=HPROJ(:,:,:)+THIS%HPROJ(:,:,:)
            DEALLOCATE(THIS%HPROJ)
            NULLIFY(THIS%HPROJ)
          END IF
!
!         ======================================================================
!         ==  ADD  |P>DH<P|PSI>                                               ==
!         ======================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          CALL TIMING$CLOCKON('W:HPSI.ADDPROJ')
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          THPSIRESIDENT=CPPAW_CUBLAS_ACC_HPSI_RESIDENCY_ENABLED()
          IF(THPSIRESIDENT) THEN
            DO IAT=1,NAT
              ISP=MAP%ISP(IAT)
              LMNX=MAP%LMNX(ISP)
              ADDPROFLOPS=8.D0*REAL(NGL,KIND=8)*REAL(LMNX,KIND=8) &
     &                   *REAL(NDIM*NBH,KIND=8)
              THPSIRESIDENT=THPSIRESIDENT &
     &             .AND.CPPAW_CUBLAS_ACC_SHOULD_USE_ADDPRODUCT(ADDPROFLOPS)
            ENDDO
          END IF
          IF(THPSIRESIDENT) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &          ('ACC_PRESENT_HPSI_ADDPRO','ACC_COPY_HPSI_ADDPRO_IN' &
     &          ,NGL,NDIM,NBH,THIS%HPSI)
#ENDIF
!$ACC ENTER DATA COPYIN(THIS%HPSI(1:NGL,1:NDIM,1:NBH))
            THIS%HPSI_ACC_RESIDENT=.TRUE.
          END IF
#ENDIF
          CALL WAVES_ADDPRO(MAP,GSET,NAT,R,NGL,NDIM,NBH,MAP%NPRO &
     &                     ,THIS%HPSI,HPROJ,'HPSI')
	CALL TIMING$CLOCKOFF('W:HPSI.ADDPROJ')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          CALL ACCELPROFILE$ADD('PAW_HPSI_ADDPRO' &
     &        ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &        ,INT(MAP%NPRO,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
          DEALLOCATE(HPROJ)
!===============================================================================
	ELSE
!===============================================================================
         CALL WAVES_HPSI(MAP,GSET,ISPIN,NGL,NDIM,NDIMD,NBH,MAP%NPRO,LMNXX,NAT,NRL&
     &                  ,THIS%PSI0,RHO(1,ISPIN),R,THIS%PROJ,DH,THIS%HPSI)
!!LMTO INTERFACE MISSING!!!
!===============================================================================
	END IF
!===============================================================================
	        ENDDO
	      ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T1)
      CALL ACCELPROFILE$ADD('PAW_HPSI_TOTAL' &
     &    ,INT(NRL,KIND=8),INT(NDIM,KIND=8),INT(NAT,KIND=8) &
     &    ,INT(MAP%NPRO,KIND=8),0.D0,0.D0,ACCEL_TOTAL_T1-ACCEL_TOTAL_T0)
#ENDIF
                              CALL TIMING$CLOCKOFF('W:HPSI')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$HPSI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_HPROJ(NDIM,NB,LMNX,DH,PROJ,HPROJ)
!     **************************************************************************
!     **  HPROJ = DH*PROC= DH*<P|PSI>                                         **
!     **  DH IS THE ON-SITE AUGMENTATION HAMILTONIAN                          **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1991)***********
      USE WAVES_MODULE, ONLY : MAP_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NDIM
      INTEGER(4),INTENT(IN)  :: NB
      INTEGER(4),INTENT(IN)  :: LMNX
      COMPLEX(8),INTENT(IN)  :: DH(LMNX,LMNX,NDIM**2)
      COMPLEX(8),INTENT(IN)  :: PROJ(NDIM,NB,LMNX)
      COMPLEX(8),INTENT(OUT) :: HPROJ(NDIM,NB,LMNX)
      INTEGER(4)             :: IB
      INTEGER(4)             :: LMN1,LMN2
      COMPLEX(8)             :: DHUPUP,DHDNDN
      COMPLEX(8)             :: DHUPDN,DHDNUP
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
#ENDIF
!     **************************************************************************
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      HPROJ(:,:,:)=(0.D0,0.D0)
!
!     ==========================================================================
!     ==  EVALUATE  DH<P|PSI>                                                 ==
!     ==========================================================================
      IF(NDIM.EQ.1) THEN
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            DHUPUP=DH(LMN1,LMN2,1)
            DO IB=1,NB
              HPROJ(1,IB,LMN1)=HPROJ(1,IB,LMN1)+DHUPUP*PROJ(1,IB,LMN2)
            ENDDO
          ENDDO
        ENDDO
      ELSE 
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            DHUPUP=DH(LMN1,LMN2,1)+DH(LMN1,LMN2,4)
            DHDNDN=DH(LMN1,LMN2,1)-DH(LMN1,LMN2,4)
!           == HERE, DH=DH/DD AND NOT CONJG(DE/DD) !! 
!JSC+PEB    DHUPDN=CMPLX(DH(LMN1,LMN2,2),DH(LMN1,LMN2,3),KIND=8)
!            DHUPDN=CMPLX(DH(LMN1,LMN2,2),-DH(LMN1,LMN2,3),KIND=8)
            DHUPDN=DH(LMN1,LMN2,2)-CI*DH(LMN1,LMN2,3)
            DHDNUP=CONJG(DHUPDN)
            DO IB=1,NB
              HPROJ(1,IB,LMN1)=HPROJ(1,IB,LMN1)+DHUPUP*PROJ(1,IB,LMN2) &
     &                                         +DHUPDN*PROJ(2,IB,LMN2)
              HPROJ(2,IB,LMN1)=HPROJ(2,IB,LMN1)+DHDNUP*PROJ(1,IB,LMN2) &
     &                                         +DHDNDN*PROJ(2,IB,LMN2)
            ENDDO
          ENDDO
        ENDDO
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_HPROJ_TOTAL' &
     &    ,INT(NB,KIND=8),INT(LMNX,KIND=8),INT(NDIM,KIND=8),0_8 &
     &    ,0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_HPSI(MAP,GSET,ISPIN,NGL,NDIM,NDIMD,NBH,NPRO,LMNXX,NAT &
     &                     ,NRL,PSI,POT,R,PROJ,DH,HPSI)
!     **************************************************************************
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1991)***********
      USE WAVES_MODULE, ONLY: MAP_TYPE,GSET_TYPE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: ISPIN
      TYPE(MAP_TYPE) ,INTENT(IN) :: MAP  !NAT,NSP,NBAREPRO,NRL,NPRO 
                                         !,LMX,LNXX,LNX,LMNX,LOX,ISP
      TYPE(GSET_TYPE),INTENT(IN) :: GSET
      INTEGER(4)     ,INTENT(IN) :: NGL
      INTEGER(4)     ,INTENT(IN) :: NDIM
      INTEGER(4)     ,INTENT(IN) :: NDIMD
      INTEGER(4)     ,INTENT(IN) :: NBH
      INTEGER(4)     ,INTENT(IN) :: NRL
      INTEGER(4)     ,INTENT(IN) :: NPRO
      INTEGER(4)     ,INTENT(IN) :: LMNXX
      INTEGER(4)     ,INTENT(IN) :: NAT
      COMPLEX(8)     ,INTENT(IN) :: PSI(NGL,NDIM,NBH)
      REAL(8)        ,INTENT(IN) :: POT(NRL,NDIM**2) !RHO(1,ISPIN)
      REAL(8)        ,INTENT(IN) :: R(3,NAT)
      COMPLEX(8)     ,INTENT(IN) :: PROJ(NDIM,NBH,NPRO)!(NDIM,NBH,NPRO)<PSPSI|P>
      COMPLEX(8)     ,INTENT(IN) :: DH(LMNXX,LMNXX,NDIMD,NAT)
      COMPLEX(8)     ,INTENT(OUT):: HPSI(NGL,NDIM,NBH)
      COMPLEX(8)     ,ALLOCATABLE:: HPROJ(:,:,:)  
      COMPLEX(8)     ,ALLOCATABLE:: DH1(:,:,:)
      INTEGER(4)                 :: IPRO,IAT,ISP,LMNX
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                    :: ACCEL_T0
      REAL(8)                    :: ACCEL_T1
      REAL(8)                    :: ACCEL_TOTAL_T0
      REAL(8)                    :: ACCEL_TOTAL_T1
#ENDIF
!     **************************************************************************
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T0)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
	CALL TIMING$CLOCKON('W:HPSI.VPSI')
      CALL WAVES_VPSI(GSET,NGL,NDIM,NBH,NRL,PSI,POT,HPSI)
	CALL TIMING$CLOCKOFF('W:HPSI.VPSI')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_HPSI_VPSI' &
     &    ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &    ,INT(NRL,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  EVALUATE  DH<P|PSI>                                                 ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
	CALL TIMING$CLOCKON('W:HPSI.HPROJ')
      ALLOCATE(HPROJ(NDIM,NBH,NPRO))
      HPROJ(:,:,:)=(0.D0,0.D0)
      IPRO=1
      DO IAT=1,NAT
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        ALLOCATE(DH1(LMNX,LMNX,NDIM**2))
        IF(NDIM.EQ.1) THEN
          DH1(:,:,1)=DH(1:LMNX,1:LMNX,ISPIN,IAT)
        ELSE 
          DH1(:,:,:)=DH(1:LMNX,1:LMNX,:,IAT)
        END IF
        CALL WAVES_HPROJ(NDIM,NBH,LMNX,DH1 &
       &       ,PROJ(:,:,IPRO:IPRO+LMNX-1),HPROJ(:,:,IPRO:IPRO+LMNX-1))
        DEALLOCATE(DH1)
        IPRO=IPRO+LMNX
      ENDDO
	CALL TIMING$CLOCKOFF('W:HPSI.HPROJ')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_HPSI_HPROJ' &
     &    ,INT(NBH,KIND=8),INT(NPRO,KIND=8),INT(NAT,KIND=8) &
     &    ,INT(NDIM,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  ADD  |P>DH<P|PSI>                                                   ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
	CALL TIMING$CLOCKON('W:HPSI.ADDPRO')
      CALL WAVES_ADDPRO(MAP,GSET,NAT,R,NGL,NDIM,NBH,NPRO,HPSI,HPROJ &
     &                 ,'HPSI')
	CALL TIMING$CLOCKOFF('W:HPSI.ADDPRO')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_HPSI_ADDPRO' &
     &    ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &    ,INT(NPRO,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
      DEALLOCATE(HPROJ)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T1)
      CALL ACCELPROFILE$ADD('PAW_HPSI_TOTAL' &
     &    ,INT(NRL,KIND=8),INT(NDIM,KIND=8),INT(NAT,KIND=8) &
     &    ,INT(NPRO,KIND=8),0.D0,0.D0,ACCEL_TOTAL_T1-ACCEL_TOTAL_T0)
#ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_OPSI(NB,NBH,NPRO,NAT,NGL,R0,PROJ,OPSI)
!     **************************************************************************
!     **                                                                      **
!     **  |PSI>+|P>DO<P|PSI>                                                  **
!     **                                                                      **
!     *********************************JOHANNES KAESTNER 2004*******************
      USE WAVES_MODULE, ONLY : MAP &
     &                        ,GSET &
     &                        ,NDIM
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)   :: NB
      INTEGER(4) ,INTENT(IN)   :: NBH
      INTEGER(4) ,INTENT(IN)   :: NPRO
      INTEGER(4) ,INTENT(IN)   :: NAT
      INTEGER(4) ,INTENT(IN)   :: NGL
      REAL(8)    ,INTENT(IN)   :: R0(3,NAT)
      COMPLEX(8) ,INTENT(IN)   :: PROJ(NDIM,NBH,NPRO)  !<P|PSI>
      COMPLEX(8) ,INTENT(INOUT):: OPSI(NGL,NDIM,NBH) !PSI ON ETRY, OPSI ON EXIT
      INTEGER(4)               :: IPRO,IAT,ISP,LNX,LMNX
      COMPLEX(8) ,ALLOCATABLE  :: OPROJ(:,:,:)
      REAL(8)    ,ALLOCATABLE  :: DO(:,:)      ! 1-CENTER OVERLAP
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                  :: ACCEL_T0
      REAL(8)                  :: ACCEL_T1
      REAL(8)                  :: ACCEL_TOTAL_T0
      REAL(8)                  :: ACCEL_TOTAL_T1
#ENDIF
!     **************************************************************************
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  EVALUATE  DO<P|PSI>  ASSUMING <PRO|PSI> IS STILL VALID              ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      ALLOCATE(OPROJ(NDIM,NBH,NPRO))
      OPROJ(:,:,:)=(0.D0,0.D0)
      IPRO=1
      DO IAT=1,NAT
         ISP=MAP%ISP(IAT)
         LNX=MAP%LNX(ISP)
         LMNX=MAP%LMNX(ISP)
         ALLOCATE(DO(LNX,LNX))
         CALL SETUP$ISELECT(ISP)
         CALL SETUP$GETR8A('DO',LNX*LNX,DO)
                                    !FORMER CALL SETUP$1COVERLAP(ISP,LNX,DO)
         CALL WAVES_OPROJ(LNX,MAP%LOX(1:LNX,ISP),DO,NDIM,LMNX,NBH &
     &                  ,PROJ(:,:,IPRO:IPRO+LMNX-1),OPROJ(:,:,IPRO:IPRO+LMNX-1))
         DEALLOCATE(DO)
         IPRO=IPRO+LMNX
         CALL SETUP$UNSELECT()
      ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OPSI_OPROJ' &
     &    ,INT(NBH,KIND=8),INT(NPRO,KIND=8),INT(NAT,KIND=8) &
     &    ,INT(NDIM,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  ADD  |PSI>+|P>DO<P|PSI>                                             ==
!     ==========================================================================
      !ALLOCATE(THIS%OPSI(NGL,NDIM,NBH))
      !THIS%OPSI(:,:,:)=THIS%PSI0(:,:,:)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      CALL WAVES_ADDPRO(MAP,GSET,NAT,R0,NGL,NDIM,NBH,NPRO,OPSI,OPROJ &
     &                 ,'OPSI')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      CALL ACCELPROFILE$ADD('PAW_OPSI_ADDPRO' &
     &    ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &    ,INT(NPRO,KIND=8),0.D0,0.D0,ACCEL_T1-ACCEL_T0)
#ENDIF
      DEALLOCATE(OPROJ)
      ! THIS THIS%OPSI IS MY O|PSI> FOR CG
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T1)
      CALL ACCELPROFILE$ADD('PAW_OPSI_TOTAL' &
     &    ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &    ,INT(NPRO,KIND=8),0.D0,0.D0,ACCEL_TOTAL_T1-ACCEL_TOTAL_T0)
#ENDIF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_VPSI(GSET,NGL,NDIM,NBH,NRL,PSI,V,HPSI)
!     **************************************************************************
!     **                                                                      **
!     **  EVALUATES H*PSI WITHOUT THE AUGMENTATION PART                       **
!     **  |HPSI>=(-0.5*NABLA**2+PS-V)|PSPSI(0)>                               **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1991)***********
      USE WAVES_MODULE, ONLY : GSET_TYPE,TBUCKET
      IMPLICIT NONE
      TYPE(GSET_TYPE),INTENT(IN) :: GSET
      INTEGER(4),INTENT(IN)  :: NGL         ! #(G-VECTORS)
      INTEGER(4),INTENT(IN)  :: NDIM        ! #(BANDS)
      INTEGER(4),INTENT(IN)  :: NBH         ! #(BANDS)
      INTEGER(4),INTENT(IN)  :: NRL         ! #(R-SPACE POINTS) 
      COMPLEX(8),INTENT(IN)  :: PSI(NGL,NDIM,NBH)  ! PSPSI(0)
      REAL(8)   ,INTENT(IN)  :: V(NRL,NDIM**2) ! PSPOT
      COMPLEX(8),INTENT(OUT) :: HPSI(NGL,NDIM,NBH)     ! PSH|PSPSI>
      REAL(8)   ,ALLOCATABLE :: G2(:)          ! G**2
      INTEGER(4)             :: IB,IR,IG,IDIM
      REAL(8)   ,ALLOCATABLE :: VUPUP(:),VDNDN(:)
      COMPLEX(8),ALLOCATABLE :: VUPDN(:)
      COMPLEX(8),ALLOCATABLE :: PSIOFR(:,:)
      COMPLEX(8)             :: PSIUP,PSIDN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
      REAL(8)                :: ACCEL_TOTAL_T0
      REAL(8)                :: ACCEL_TOTAL_T1
      REAL(8)                :: ACCEL_FFT_GTOR_T
      REAL(8)                :: ACCEL_FFT_RTOG_T
      REAL(8)                :: ACCEL_POT_T
      REAL(8)                :: ACCEL_KIN_T
      REAL(8)                :: ACCEL_BUCKET_T
#ENDIF
!     **************************************************************************
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      ACCEL_FFT_GTOR_T=0.D0
      ACCEL_FFT_RTOG_T=0.D0
      ACCEL_POT_T=0.D0
      ACCEL_KIN_T=0.D0
      ACCEL_BUCKET_T=0.D0
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  MULTIPLY WAVE FUNCTIONS WITH THE POTENTIAL                          ==
!     ==                                                                      ==
!     ==  THE FFTS ARE DONE FOR EACH WAVE FUNCTION INDIVIUDALLY TO AVOID      ==
!     ==  A MEMORY SPIKE. THERE MAY BE A LOSS OF SPEED BECAUSE DOING THE      ==
!     ==  FFTS AND THE MULTIPLICATION IN ONE SHOT WOULD BE MORE EFFICIENT.    ==
!     ==  IN THE PREVIOUS IMPLEMENTATION THIS HAD NOT BEEN EXPLOITED ANYWAY.  ==
!     ==========================================================================
      ALLOCATE(PSIOFR(NRL,NDIM))
      IF(NDIM.EQ.1) THEN
        DO IB=1,NBH
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          CALL PLANEWAVE$FFT('GTOR',NDIM,NGL,PSI(:,:,IB),NRL,PSIOFR)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          ACCEL_FFT_GTOR_T=ACCEL_FFT_GTOR_T+ACCEL_T1-ACCEL_T0
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          DO IR=1,NRL
            PSIOFR(IR,1)=V(IR,1)*PSIOFR(IR,1)
          ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          ACCEL_POT_T=ACCEL_POT_T+ACCEL_T1-ACCEL_T0
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          CALL PLANEWAVE$FFT('RTOG',NDIM,NGL,HPSI(:,:,IB),NRL,PSIOFR)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          ACCEL_FFT_RTOG_T=ACCEL_FFT_RTOG_T+ACCEL_T1-ACCEL_T0
#ENDIF
        ENDDO
      ELSE
        ALLOCATE(VUPUP(NRL))
        ALLOCATE(VDNDN(NRL))
        ALLOCATE(VUPDN(NRL))
        DO IR=1,NRL
          VUPUP(IR)=V(IR,1)+V(IR,4)
          VDNDN(IR)=V(IR,1)-V(IR,4)
          VUPDN(IR)=CMPLX(V(IR,2),-V(IR,3),KIND=8)
        ENDDO
        DO IB=1,NBH
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          CALL PLANEWAVE$FFT('GTOR',NDIM,NGL,PSI(:,:,IB),NRL,PSIOFR)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          ACCEL_FFT_GTOR_T=ACCEL_FFT_GTOR_T+ACCEL_T1-ACCEL_T0
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          DO IR=1,NRL
            PSIUP=PSIOFR(IR,1)
            PSIDN=PSIOFR(IR,2)
            PSIOFR(IR,1)=VUPUP(IR)*PSIUP+      VUPDN(IR) *PSIDN
            PSIOFR(IR,2)=VDNDN(IR)*PSIDN+CONJG(VUPDN(IR))*PSIUP
          ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          ACCEL_POT_T=ACCEL_POT_T+ACCEL_T1-ACCEL_T0
          CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
          CALL PLANEWAVE$FFT('RTOG',NDIM,NGL,HPSI(:,:,IB),NRL,PSIOFR)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL ACCELPROFILE$NOW(ACCEL_T1)
          ACCEL_FFT_RTOG_T=ACCEL_FFT_RTOG_T+ACCEL_T1-ACCEL_T0
#ENDIF
        ENDDO
        DEALLOCATE(VUPUP)
        DEALLOCATE(VUPDN)
        DEALLOCATE(VDNDN)
      END IF
      DEALLOCATE(PSIOFR)
!
!     ==========================================================================
!     ==  ADD KINETIC ENERGY CONTRIBUTION                                     ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      ALLOCATE(G2(NGL))
      CALL PLANEWAVE$GETR8A('G2',NGL,G2)
      DO IB=1,NBH
        DO IDIM=1,NDIM
          DO IG=1,NGL
            HPSI(IG,IDIM,IB)=HPSI(IG,IDIM,IB)+0.5D0*G2(IG)*PSI(IG,IDIM,IB)
          ENDDO
        ENDDO
      ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      ACCEL_KIN_T=ACCEL_T1-ACCEL_T0
#ENDIF
!
!     ==========================================================================
!     ==  ADD BUCKET POTENTIAL                                                ==
!     ==========================================================================
      IF(TBUCKET) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
        DO IB=1,NBH
          DO IDIM=1,NDIM
            DO IG=1,NGL
              HPSI(IG,IDIM,IB)=HPSI(IG,IDIM,IB) &
     &                        +GSET%BUCKET(IG)*PSI(IG,IDIM,IB)
            ENDDO
          ENDDO
        ENDDO
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T1)
        ACCEL_BUCKET_T=ACCEL_T1-ACCEL_T0
#ENDIF
      END IF
      DEALLOCATE(G2)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$ADD('PAW_VPSI_FFT_GTOR' &
     &    ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &    ,INT(NRL,KIND=8),0.D0,0.D0,ACCEL_FFT_GTOR_T)
      CALL ACCELPROFILE$ADD('PAW_VPSI_POT' &
     &    ,INT(NRL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &    ,0_8,0.D0,0.D0,ACCEL_POT_T)
      CALL ACCELPROFILE$ADD('PAW_VPSI_FFT_RTOG' &
     &    ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &    ,INT(NRL,KIND=8),0.D0,0.D0,ACCEL_FFT_RTOG_T)
      CALL ACCELPROFILE$ADD('PAW_VPSI_KIN' &
     &    ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &    ,0_8,0.D0,0.D0,ACCEL_KIN_T)
      IF(TBUCKET) THEN
        CALL ACCELPROFILE$ADD('PAW_VPSI_BUCKET' &
     &      ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &      ,0_8,0.D0,0.D0,ACCEL_BUCKET_T)
      END IF
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T1)
      CALL ACCELPROFILE$ADD('PAW_VPSI_TOTAL' &
     &    ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &    ,INT(NRL,KIND=8),0.D0,0.D0,ACCEL_TOTAL_T1-ACCEL_TOTAL_T0)
#ENDIF
      RETURN
      END
!!$!
!!$!     ..................................................................
!!$      SUBROUTINE WAVES_EKIN_OLD(NGL,NDIM,NBH,NB,F,GWEIGHT,PSI,EKIN &
!!$     &                         ,TSTRESS,STRESS,TBUCKET,BUCKET,DBUCKET)
!!$!     ******************************************************************
!!$!     **                                                              **
!!$!     **  EVALUATE PS KINETIC ENERGY IN G-SPACE                       **
!!$!     **  EVALUATE NUMBER OF ELECTRONS IN G-SPACE                     **
!!$!     **                                                              **
!!$!     **  REMARKS:                                                    **
!!$!     **    REQUIRES PLANEWAVE OBJECT TO BE SET PROPERLY              **
!!$!     **                                                              **
!!$!     *******************************************P.E. BLOECHL, (1999)***
!!$      USE MPE_MODULE
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: NB         ! #(STATES)
!!$      INTEGER(4),INTENT(IN) :: NBH        ! #(WAVE FUNCTIONS)
!!$      INTEGER(4),INTENT(IN) :: NDIM       ! #(SPINOR COMPONENTS)
!!$      INTEGER(4),INTENT(IN) :: NGL        ! #(PLANE WAVES)
!!$      REAL(8)   ,INTENT(IN) :: F(NB)      ! OCCUPATION
!!$      REAL(8)   ,INTENT(IN) :: GWEIGHT
!!$      COMPLEX(8),INTENT(IN) :: PSI(NGL,NDIM,NBH) ! PS-WAVE FUNCTION
!!$      REAL(8)   ,INTENT(OUT):: EKIN       ! KINETIC ENERGY
!!$      LOGICAL(4),INTENT(IN) :: TSTRESS    ! SWITCH STRESS ON/OFF
!!$      REAL(8)   ,INTENT(OUT):: STRESS(3,3)! STRESS
!!$      LOGICAL(4),INTENT(IN) :: TBUCKET    ! BUCKET POTENTIAL PRESENT
!!$      REAL(8)   ,INTENT(IN) :: BUCKET(NGL) ! BUCKET POTENTIAL
!!$      REAL(8)   ,INTENT(IN) :: DBUCKET(NGL) ! 1/G *DBUCKET/DG
!!$      REAL(8)   ,ALLOCATABLE:: G2(:)      ! G**2
!!$      REAL(8)   ,ALLOCATABLE:: GVEC(:,:)  ! G   
!!$      COMPLEX(8),ALLOCATABLE:: PSI1(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: DMAT(:)
!!$      LOGICAL(4)            :: TINV
!!$      INTEGER(4)            :: IB,IG,IDIM1,IDIM2
!!$      REAL(8)               :: FP,FM
!!$      REAL(8)   ,PARAMETER  :: DSMALL=1.D-12
!!$      REAL(8)               :: EBUCKET
!!$      REAL(8)               :: SVAR
!!$!     ******************************************************************
!!$      CALL ERROR$MSG('ROUTINE MARKED FOR DELETION. CONTAINS ERRORS!')
!!$      CALL ERROR$STOP('WAVES_EKIN_OLD')
!!$!
!!$!     ==================================================================
!!$!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT==
!!$!     ==================================================================
!!$      CALL PLANEWAVE$GETL4('TINV',TINV)
!!$      IF(TINV) THEN
!!$        IF(NBH.NE.(NB+1)/2) THEN
!!$          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
!!$          CALL ERROR$STOP('WAVES_EKIN')
!!$        END IF
!!$      ELSE 
!!$        IF(NBH.NE.NB) THEN
!!$          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
!!$          CALL ERROR$STOP('WAVES_EKIN')
!!$        END IF
!!$      END IF
!!$!
!!$!     ==================================================================
!!$!     ==  CALCULATE DMAT(G)=F(IB)*PSI^*(G,IB)*PSI(G,IB)               ==
!!$!     ==================================================================
!!$      ALLOCATE(DMAT(NGL))
!!$      ALLOCATE(PSI1(NGL,NDIM))
!!$      DMAT(:)=0.D0
!!$      DO IB=1,NBH
!!$!       == DETERMINE OCCUPATIONS =======================================
!!$        IF(TINV) THEN
!!$          FP=0.5D0*(F(2*IB-1)+F(2*IB))
!!$          FM=0.5D0*(F(2*IB-1)-F(2*IB))
!!$        ELSE
!!$          FP=F(IB)
!!$          FM=0.D0
!!$        END IF
!!$!
!!$!       ================================================================
!!$!       == GENERAL WAVE FUNCTION / FIRST PART OF SUPER WAVE FUNCTIONS ==
!!$!       == <PSI_+|G^2|PSI_+> FOR SUPER WAVE FUNCTIONS                 ==
!!$!       ================================================================
!!$        IF(FP.EQ.0.D0.AND.FM.EQ.0.D0) CYCLE
!!$        DO IDIM1=1,NDIM
!!$          DO IDIM2=IDIM1,NDIM
!!$            SVAR=FP
!!$            IF(IDIM1.NE.IDIM2)SVAR=2.D0*SVAR
!!$            DO IG=1,NGL
!!$              DMAT(IG)=DMAT(IG) &
!!$    &              +SVAR*REAL(CONJG(PSI(IG,IDIM1,IB))*PSI(IG,IDIM2,IB),KIND=8)
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$!
!!$!       ================================================================
!!$!       ==  <PSI_+|G^2|PSI_->                                         ==
!!$!       ================================================================
!!$        IF(.NOT.TINV.OR.ABS(FM).LT.DSMALL) CYCLE
!!$        DO IDIM1=1,NDIM
!!$          CALL PLANEWAVE$INVERTG(NGL,PSI(1,IDIM1,IB),PSI1(1,IDIM1))
!!$        ENDDO
!!$        DO IDIM1=1,NDIM
!!$          DO IDIM2=IDIM1,NDIM
!!$            SVAR=FM
!!$            IF(IDIM1.NE.IDIM2)SVAR=2.D0*SVAR
!!$            DO IG=1,NGL
!!$              DMAT(IG)=DMAT(IG) &
!!$     &           +SVAR*REAL(CONJG(PSI(IG,IDIM1,IB))*PSI1(IG,IDIM2),KIND=8)
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$      DEALLOCATE(PSI1)
!!$!
!!$!     ==================================================================
!!$!     ==  CALCULATE KINETIC ENERGY AND STRESS                         ==
!!$!     ==================================================================
!!$      IF(TSTRESS) THEN
!!$        ALLOCATE(GVEC(3,NGL))
!!$        CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
!!$!       == KINETIC ENERGY AND KINETIC STRESS ===========================
!!$        STRESS(:,:)=0.D0
!!$        DO IG=1,NGL
!!$          STRESS(1,1)=STRESS(1,1)+GVEC(1,IG)*GVEC(1,IG)*DMAT(IG)
!!$          STRESS(1,2)=STRESS(1,2)+GVEC(1,IG)*GVEC(2,IG)*DMAT(IG)
!!$          STRESS(1,3)=STRESS(1,3)+GVEC(1,IG)*GVEC(3,IG)*DMAT(IG)
!!$          STRESS(2,2)=STRESS(2,2)+GVEC(2,IG)*GVEC(2,IG)*DMAT(IG)
!!$          STRESS(2,3)=STRESS(2,3)+GVEC(2,IG)*GVEC(3,IG)*DMAT(IG)
!!$          STRESS(3,3)=STRESS(3,3)+GVEC(3,IG)*GVEC(3,IG)*DMAT(IG)
!!$        ENDDO
!!$        EKIN=0.5D0*(STRESS(1,1)+STRESS(2,2)+STRESS(3,3))
!!$        STRESS=-STRESS
!!$!       == ENERGY AND STRESS DUE TO THE BUCKET POTENTIAL ===============
!!$        IF(TBUCKET) THEN
!!$          EBUCKET=0.D0
!!$          DO IG=1,NGL
!!$            IF(BUCKET(IG).EQ.0.D0) CYCLE
!!$            EBUCKET=EBUCKET+BUCKET(IG)*DMAT(IG)
!!$            SVAR=DMAT(IG)*DBUCKET(IG)
!!$            STRESS(1,1)=STRESS(1,1)-GVEC(1,IG)*GVEC(1,IG)*SVAR
!!$            STRESS(1,2)=STRESS(1,2)-GVEC(1,IG)*GVEC(2,IG)*SVAR
!!$            STRESS(1,3)=STRESS(1,3)-GVEC(1,IG)*GVEC(3,IG)*SVAR
!!$            STRESS(2,2)=STRESS(2,2)-GVEC(2,IG)*GVEC(2,IG)*SVAR
!!$            STRESS(2,3)=STRESS(2,3)-GVEC(2,IG)*GVEC(3,IG)*SVAR
!!$            STRESS(3,3)=STRESS(3,3)-GVEC(3,IG)*GVEC(3,IG)*SVAR
!!$          ENDDO
!!$          EKIN=EKIN+EBUCKET
!!$        END IF
!!$!       == PARALLELIZE AND CLOSE DOWN ==================================
!!$        STRESS(2,1)=STRESS(1,2)
!!$        STRESS(3,1)=STRESS(1,3)
!!$        STRESS(3,2)=STRESS(2,3)
!!$        EKIN=EKIN*GWEIGHT
!!$        STRESS=STRESS*GWEIGHT
!!$        CALL MPE$COMBINE('NONE','+',EKIN)
!!$        CALL MPE$COMBINE('NONE','+',STRESS)
!!$        DEALLOCATE(GVEC)
!!$        DEALLOCATE(DMAT)
!!$      ELSE
!!$        ALLOCATE(G2(NGL))
!!$        CALL PLANEWAVE$GETR8A('G2',NGL,G2)
!!$!       == KINETIC ENERGY ==============================================
!!$        EKIN=0.D0
!!$        DO IG=1,NGL
!!$          EKIN=EKIN+G2(IG)*DMAT(IG)
!!$        ENDDO
!!$        EKIN=0.5D0*EKIN
!!$!       == BUCKET POTENTIAL ============================================
!!$        IF(TBUCKET) THEN
!!$          EBUCKET=0.D0
!!$          DO IG=1,NGL
!!$            EBUCKET=EBUCKET+BUCKET(IG)*DMAT(IG)
!!$          ENDDO
!!$          EKIN=EKIN+EBUCKET
!!$        END IF
!!$!       ==  PARALLELIZE AND CLOSE DOWN =================================
!!$        EKIN=EKIN*GWEIGHT
!!$        CALL MPE$COMBINE('NONE','+',EKIN)
!!$        DEALLOCATE(DMAT)
!!$        DEALLOCATE(G2)
!!$      END IF
!!$!
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_DENSITY(NGL,NRL,NDIM,NB,NBH,F,PSIOFG,RHO,TRHOKIN,RHOKIN)
!     **************************************************************************
!     **  CALCULATE ELECTRON DENSITY RHO IN REAL SPACE                        **
!     **  IF(TRHOKIN) CALCULATE ALSO KINETIC-ENERGY DENSITY RHOKIN IN R. SPACE**
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1991-2021)******
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NGL         ! MAX # PLANE WAVES
      INTEGER(4),INTENT(IN)  :: NRL         ! # R-SPACE POINTS
      INTEGER(4),INTENT(IN)  :: NDIM        ! DIMENSION OF THE WAVE FUNCTION
      INTEGER(4),INTENT(IN)  :: NB          ! # BANDS
      INTEGER(4),INTENT(IN)  :: NBH
      REAL(8)   ,INTENT(IN)  :: F(NB)       ! OCCUPATIONS
      COMPLEX(8),INTENT(IN)  :: PSIOFG(NGL,NDIM,NBH)
      REAL(8)   ,INTENT(OUT) :: RHO(NRL,NDIM**2) ! DENSITY IN R-SPACE
      LOGICAL(4),INTENT(IN)  :: TRHOKIN
      REAL(8)   ,INTENT(OUT) :: RHOKIN(NRL,NDIM**2) ! KINETIC-ENERGY DENSITY
      COMPLEX(8),ALLOCATABLE :: PSIOFR(:,:,:)
      COMPLEX(8),ALLOCATABLE :: PSIKINOFR(:,:,:)
      COMPLEX(8),ALLOCATABLE :: EI2KR(:)      !(NRL) SQUARED BLOCH PHASE FACTOR
      COMPLEX(8),ALLOCATABLE :: PSI1(:)          !(NRL)
      REAL(8)   ,ALLOCATABLE :: G2(:)            !(NGL)
      COMPLEX(8),ALLOCATABLE :: PSIKINOFG(:,:,:) !(NGL,NDIM,NBH)
      COMPLEX(8)             :: CSVAR
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IBH,IR,IDIM
      REAL(8)                :: F1,F2
      REAL(8)                :: RE,IM
      REAL(8)                :: SVAR1,SVAR2
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                :: ACCEL_T0
      REAL(8)                :: ACCEL_T1
      REAL(8)                :: ACCEL_TOTAL_T0
      REAL(8)                :: ACCEL_TOTAL_T1
      REAL(8)                :: ACCEL_FFT_T
      REAL(8)                :: ACCEL_KIN_PREP_T
      REAL(8)                :: ACCEL_KIN_FFT_T
      REAL(8)                :: ACCEL_ACCUM_T
#ENDIF
!     **************************************************************************
                               CALL TRACE$PUSH('WAVES_DENSITY')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      ACCEL_FFT_T=0.D0
      ACCEL_KIN_PREP_T=0.D0
      ACCEL_KIN_FFT_T=0.D0
      ACCEL_ACCUM_T=0.D0
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T0)
#ENDIF
!RELEASED 8.OCT.99 
!
!     ==========================================================================
!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT        ==
!     ==========================================================================
      CALL PLANEWAVE$GETL4('TINV',TINV)
      IF(TINV) THEN
        IF(NBH.NE.(NB+1)/2) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_DENSITY')
        END IF
      ELSE 
        IF(NBH.NE.NB) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_DENSITY')
        END IF
      END IF
!
!     ==========================================================================
!     ==  FFT                                                                 ==
!     ==========================================================================
      ALLOCATE(PSIOFR(NRL,NDIM,NBH))
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      CALL PLANEWAVE$FFT('GTOR',NBH*NDIM,NGL,PSIOFG,NRL,PSIOFR)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      ACCEL_FFT_T=ACCEL_T1-ACCEL_T0
#ENDIF
!
!     == CONSTRUCT <R| 1/2 * G^2 |PSITILDE> ====================================
      IF(TRHOKIN) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
        ALLOCATE(G2(NGL))
        CALL PLANEWAVE$GETR8A('G2',NGL,G2)
!       == THIS WILL GIVE A SPIKE IN THE MEMORY REQUIREMENT, ===================
!       == BECAUSE AN ADDITIONAL SET OF WAVE FUNCTIONS IS ALLOCATED ============
!       == BOTH IN A PLANE WAVE BASIS AND IN REAL SPACE ========================
        ALLOCATE(PSIKINOFG(NGL,NDIM,NBH))
        DO IBH=1,NBH
          DO IDIM=1,NDIM
            PSIKINOFG(:,IDIM,IBH)=0.5D0*G2(:)*PSIOFG(:,IDIM,IBH)
          ENDDO
        ENDDO
        ALLOCATE(PSIKINOFR(NRL,NDIM,NBH))
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T1)
        ACCEL_KIN_PREP_T=ACCEL_T1-ACCEL_T0
        CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
        CALL PLANEWAVE$FFT('GTOR',NBH*NDIM,NGL,PSIKINOFG,NRL,PSIKINOFR)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        CALL ACCELPROFILE$NOW(ACCEL_T1)
        ACCEL_KIN_FFT_T=ACCEL_T1-ACCEL_T0
#ENDIF
        DEALLOCATE(PSIKINOFG)
        DEALLOCATE(G2)
      END IF
!
!     ==========================================================================
!     ==  CALCULATE CHARGE DENSITY                                            ==
!     ==========================================================================
!
!     ==========================================================================
!     ==  CALCULATE DENSITY FOR SUPER WAVE FUNCTIONS                          ==
!     ==========================================================================
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
      IF(TINV) THEN
        IF(NDIM.EQ.2) THEN
          CALL ERROR$MSG('SPINOR WAVE FUNCTIONS ARE NOT COMPATIBLE')
          CALL ERROR$MSG('WITH SUPER WAVE FUNCTIONS')
          CALL ERROR$STOP('WAVES_DENSITY')
        END IF
        ALLOCATE(EI2KR(NRL))
        CALL PLANEWAVE$GETC8A('EIKR',NRL,EI2KR)
        DO IR=1,NRL
          EI2KR(IR)=EI2KR(IR)**2
        ENDDO
        ALLOCATE(PSI1(NRL))
        RHO(:,:)=0.D0
        IF(TRHOKIN) RHOKIN(:,:)=0.D0
        DO IBH=1,NBH
          F1=F(2*IBH-1)
          F2=0.D0
          IF(2*IBH.LE.NB) F2=F(2*IBH)
          SVAR1=0.5D0*(F1+F2)
          SVAR2=0.5D0*(F1-F2)
!         == NOTE THAT OOCUPATIONS CAN BE NEGATIVE FROM K-INTEGRATION ==========
          IF(SVAR1.EQ.0.D0) THEN
            IF(SVAR2.EQ.0.D0) CYCLE
          END IF
          DO IR=1,NRL
            PSI1(IR)=SVAR1*CONJG(PSIOFR(IR,1,IBH))
          ENDDO
          IF(SVAR2.NE.0.D0) THEN
            DO IR=1,NRL
              PSI1(IR)=PSI1(IR)+SVAR2*PSIOFR(IR,1,IBH)*EI2KR(IR)
            ENDDO
          END IF
          DO IR=1,NRL
            RHO(IR,1)=RHO(IR,1)+REAL(PSIOFR(IR,1,IBH)*PSI1(IR),KIND=8)
          ENDDO
          IF(TRHOKIN) THEN
            DO IR=1,NRL
              RHOKIN(IR,1)=RHOKIN(IR,1) &
       &                  +REAL(PSIKINOFR(IR,1,IBH)*PSI1(IR),KIND=8)
            ENDDO
          END IF
        ENDDO
        DEALLOCATE(PSI1)
        DEALLOCATE(EI2KR)
!
!     ==========================================================================
!     ==  CALCULATE DENSITY FOR GENERAL WAVE FUNCTIONS                        ==
!     ==========================================================================
      ELSE
        IF(NDIM.EQ.1) THEN
!       == SCALAR WAVE FUNCTIONS ===============================================
          RHO(:,:)=0.D0
          DO IBH=1,NB
            F1=F(IBH)
            IF(F1.EQ.0.D0) CYCLE
            DO IR=1,NRL
              RE= REAL(PSIOFR(IR,1,IBH),KIND=8)
              IM=AIMAG(PSIOFR(IR,1,IBH))
              RHO(IR,1)=RHO(IR,1)+F1*(RE**2+IM**2)
            ENDDO
          ENDDO
          IF(TRHOKIN) THEN
            RHOKIN(:,:)=0.D0
            DO IBH=1,NB
              F1=F(IBH)
              IF(F1.EQ.0.D0) CYCLE
              DO IR=1,NRL
                RHO(IR,1)=RHO(IR,1) &
                         +F1*REAL(PSIOFR(IR,1,IBH)*CONJG(PSIKINOFR(IR,1,IBH)),8)
              ENDDO
            ENDDO
          END IF
        ELSE
!       == SPINOR WAVE FUNCTIONS ===============================================
          RHO(:,:)=0.D0
          DO IBH=1,NBH
            F1=F(IBH)
            IF(F1.EQ.0.D0) CYCLE
!           == REAL(RHO11), REAL(RHO22), RE(RHO12), IM(RHO12) ==================
            DO IR=1,NRL
              RHO(IR,1)=RHO(IR,1) &
    &               +F1*REAL(PSIOFR(IR,1,IBH)*CONJG(PSIOFR(IR,1,IBH)),KIND=8)
              RHO(IR,4)=RHO(IR,4) &
    &               +F1*REAL(PSIOFR(IR,2,IBH)*CONJG(PSIOFR(IR,2,IBH)),KIND=8)
              CSVAR=PSIOFR(IR,1,IBH)*CONJG(PSIOFR(IR,2,IBH))
              RHO(IR,2)=RHO(IR,2)+F1*REAL(CSVAR,KIND=8)
              RHO(IR,3)=RHO(IR,3)+F1*AIMAG(CSVAR)
            ENDDO
          ENDDO
!         == TRANSFORM TO NT,NX,NY,NZ ==========================================
          DO IR=1,NRL
            SVAR1=RHO(IR,1)
            SVAR2=RHO(IR,4)
            RHO(IR,1)=SVAR1+SVAR2
            RHO(IR,4)=SVAR1-SVAR2
            RHO(IR,2)= 2.D0*RHO(IR,2)
            RHO(IR,3)=-2.D0*RHO(IR,3)
          ENDDO
          IF(TRHOKIN) THEN
            RHOKIN(:,:)=0.D0
            DO IBH=1,NBH
              F1=F(IBH)
              IF(F1.EQ.0.D0) CYCLE
!             == REAL(RHO11), REAL(RHO22), RE(RHO12), IM(RHO12) ================
              DO IR=1,NRL
                RHOKIN(IR,1)=RHOKIN(IR,1) &
    &               +F1*REAL(PSIOFR(IR,1,IBH)*CONJG(PSIKINOFR(IR,1,IBH)),KIND=8)
                RHOKIN(IR,4)=RHOKIN(IR,4) &
    &               +F1*REAL(PSIOFR(IR,2,IBH)*CONJG(PSIKINOFR(IR,2,IBH)),KIND=8)
                CSVAR=PSIOFR(IR,1,IBH)*CONJG(PSIKINOFR(IR,2,IBH))
                RHOKIN(IR,2)=RHOKIN(IR,2)+F1*REAL(CSVAR,KIND=8)
                RHOKIN(IR,3)=RHOKIN(IR,3)+F1*AIMAG(CSVAR)
              ENDDO
            ENDDO
!           == TRANSFORM TO NT,NX,NY,NZ ========================================
            DO IR=1,NRL
              SVAR1=RHOKIN(IR,1)
              SVAR2=RHOKIN(IR,4)
              RHOKIN(IR,1)=SVAR1+SVAR2
              RHOKIN(IR,4)=SVAR1-SVAR2
              RHOKIN(IR,2)= 2.D0*RHOKIN(IR,2)
              RHOKIN(IR,3)=-2.D0*RHOKIN(IR,3)
            ENDDO
          END IF
        ENDIF
      END IF
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_T1)
      ACCEL_ACCUM_T=ACCEL_T1-ACCEL_T0
#ENDIF
      DEALLOCATE(PSIOFR)
      IF(TRHOKIN) DEALLOCATE(PSIKINOFR)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$ADD('PAW_DENSITY_FFT' &
     &    ,INT(NGL,KIND=8),INT(NRL,KIND=8),INT(NDIM,KIND=8) &
     &    ,INT(NB,KIND=8),0.D0,0.D0,ACCEL_FFT_T)
      IF(TRHOKIN) THEN
        CALL ACCELPROFILE$ADD('PAW_DENSITY_KIN_PREP' &
     &      ,INT(NGL,KIND=8),INT(NRL,KIND=8),INT(NDIM,KIND=8) &
     &      ,INT(NB,KIND=8),0.D0,0.D0,ACCEL_KIN_PREP_T)
        CALL ACCELPROFILE$ADD('PAW_DENSITY_KIN_FFT' &
     &      ,INT(NGL,KIND=8),INT(NRL,KIND=8),INT(NDIM,KIND=8) &
     &      ,INT(NB,KIND=8),0.D0,0.D0,ACCEL_KIN_FFT_T)
      END IF
      CALL ACCELPROFILE$ADD('PAW_DENSITY_ACCUM' &
     &    ,INT(NRL,KIND=8),INT(NDIM,KIND=8),INT(NB,KIND=8) &
     &    ,0_8,0.D0,0.D0,ACCEL_ACCUM_T)
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T1)
      CALL ACCELPROFILE$ADD('PAW_DENSITY_TOTAL' &
     &    ,INT(NGL,KIND=8),INT(NRL,KIND=8),INT(NDIM,KIND=8) &
     &    ,INT(NB,KIND=8),0.D0,0.D0,ACCEL_TOTAL_T1-ACCEL_TOTAL_T0)
#ENDIF
                               CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_DEDPROJ(NDIM,NBH,NB,LNX,LOX,LMNX,OCC &
     &                       ,PROJ,DH,DO,EPS,DEDPROJ)
!     **************************************************************************
!     **  CALCULATES THE DERIVATIVE OF THE ONE-CENTER                         **
!     **  ENERGIES WITH RESPECT TO THE PROJECTOR FUNCTIONS                    **
!     **    DEDPROJ = DE/(DPROJ)                                              **
!     **            = DH<P|PSITILDE>F-DO<P|PSITILDE>*LAMBDA                   **
!     **  LAMBDA IS CALCULATED FROM THIS%RLAM0 BY MULTIPLICATION WITH         **
!     **  0.5*(FI+FJ)                                                         **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1999)***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NDIM      ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN)   :: NB        ! #(BANDS)
      INTEGER(4),INTENT(IN)   :: NBH       ! #(WAVE FUNCTIONS)
      INTEGER(4),INTENT(IN)   :: LMNX      ! #(PROJECTORS ON THIS SITE)
      REAL(8)   ,INTENT(IN)   :: OCC(NB)   ! OCCUPATIONS
      COMPLEX(8),INTENT(IN)   :: PROJ(NDIM,NBH,LMNX) ! <P|PSI>
      COMPLEX(8),INTENT(IN)   :: DH(LMNX,LMNX,NDIM**2) ! DE/DD
      INTEGER(4),INTENT(IN)   :: LNX
      INTEGER(4),INTENT(IN)   :: LOX(LNX)
      REAL(8)   ,INTENT(IN)   :: DO(LNX,LNX) ! DO/DD
      COMPLEX(8),INTENT(IN)   :: EPS(NB,NB)
      COMPLEX(8),INTENT(OUT)  :: DEDPROJ(NDIM,NBH,LMNX)    ! DE/D<P|
      LOGICAL(4)              :: TINV
      COMPLEX(8),ALLOCATABLE  :: OPROJ(:,:,:)
      COMPLEX(8),ALLOCATABLE  :: LAMBDA(:,:)
      REAL(8)                 :: F1,F2
      INTEGER(4)              :: IB,IB1,IB2,LMN,IDIM
      COMPLEX(8)              :: CSVAR
!     **************************************************************************
                               CALL TRACE$PUSH('WAVES_DEDPROJ')
!
!     ==========================================================================
!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT        ==
!     ==========================================================================
      CALL PLANEWAVE$GETL4('TINV',TINV)
      IF(TINV) THEN
        IF(NBH.NE.(NB+1)/2) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_DEDPROJ')
        END IF
      ELSE 
        IF(NBH.NE.NB) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_DEDPROJ')
        END IF
      END IF
!       
!     ==========================================================================
!     ==  HPROJ = DH<P|PSI>                                                   ==
!     ==========================================================================
      CALL WAVES_HPROJ(NDIM,NBH,LMNX,DH,PROJ,DEDPROJ)
!
!     ==========================================================================
!     ==  HPROJ = DH<P|PSI>*F                                                 ==
!     ==========================================================================
      IF(TINV) THEN
        DO IB=1,NBH
          IB1=2*IB-1       
          IB2=IB1+1
          F1=0.5D0*(OCC(IB1)+OCC(IB2))
          F2=0.5D0*(OCC(IB1)-OCC(IB2))
          DO LMN=1,LMNX
            DO IDIM=1,NDIM
              CSVAR=DEDPROJ(IDIM,IB,LMN)
              DEDPROJ(IDIM,IB,LMN)=F1*CSVAR+F2*CONJG(CSVAR)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO IB=1,NB
          F1=OCC(IB)
          DO LMN=1,LMNX
            DO IDIM=1,NDIM
              DEDPROJ(IDIM,IB,LMN)=DEDPROJ(IDIM,IB,LMN)*F1
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  ADD -DO<P|PSI>LAMBDA*(FI+FJ)/2                                      ==
!     ==========================================================================
      ALLOCATE(LAMBDA(NB,NB))
      DO IB1=1,NB
        F1=OCC(IB1)
        DO IB2=1,NB
          F2=OCC(IB2)
          LAMBDA(IB1,IB2)=-EPS(IB1,IB2)*0.5D0*(F1+F2)
        ENDDO
      ENDDO
      ALLOCATE(OPROJ(NDIM,NBH,LMNX))
      CALL WAVES_OPROJ(LNX,LOX,DO,NDIM,LMNX,NBH,PROJ,OPROJ)
      CALL WAVES_ADDOPROJ(LMNX,NDIM,NBH,NB,DEDPROJ,OPROJ,LAMBDA)
      DEALLOCATE(LAMBDA)
      DEALLOCATE(OPROJ)
                               CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_DEDPRO(TINV,NGL,NDIM,NBH,PSI,LMNX,DEDPROJ,DEDPRO)
!     **************************************************************************
!     **  CALCULATES THE FUNCTION |DEDPRO_I>=DE/D<PRO_I|                      **
!     **  USED TO CALCULATE FORCES AND STRESSES AS                            **
!     **    DE = SUM_I <DPRO_I|DEDPRO_I> + <DEDPRO_I|DPRO_I>                  **
!     **  WHERE                                                               **
!     **    |DEDPRO_I>=|DEDPRO_II> + SUM_N |PSI_N> F_N <PSI_N|PRO_I>          **
!     **                   +(DH-DO EPS)                                       **
!     **                                                                      **
!     **   REMARKS:                                                           **
!     **    - CALL FOR EACH ATOM INDIVIDUALLY                                 **
!     **    - INITIALIZE DEDPRO=(0.D0,0.D0) BEFORE CALLING                    **
!     **    - FOR LSD INITIALIZE ONCE AND CALL FOR SPIN UP AND DOWN           **
!     **    - DEDPRO SHOULD BE A REAL FUNCTION IN REAL SPACE                  **
!     **       BUT THAT THE SYMMETRY IS NOT ENFORCES YET                      **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1999)***********
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)   :: TINV
      INTEGER(4),INTENT(IN)   :: NGL       ! (#(G-VECTORS)
      INTEGER(4),INTENT(IN)   :: NDIM      ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN)   :: NBH       ! #(WAVE FUNCTIONS)
      COMPLEX(8),INTENT(IN)   :: PSI(NGL,NDIM,NBH)   ! |PSI>
      INTEGER(4),INTENT(IN)   :: LMNX      ! #(PROJECTORS ON THIS SITE)
      COMPLEX(8),INTENT(IN)   :: DEDPROJ(NDIM,NBH,LMNX) ! <P|PSI>
      COMPLEX(8),INTENT(OUT)  :: DEDPRO(NGL,LMNX)    ! DE/D<P|
      INTEGER(4)              :: LMN,IG
      COMPLEX(8),ALLOCATABLE  :: PSIM(:)
      LOGICAL(4),PARAMETER    :: TESSL=.TRUE.
      COMPLEX(8)              :: DEDPROJ1(NDIM,NBH,LMNX) ! CONJG(<P|PSI>)
!     **************************************************************************
!
!     ==========================================================================
!     ==  SUPERPOSE WAVE FUNCTION TO OBTAIN DE/D<P|                           ==
!     ==========================================================================
      DEDPROJ1=CONJG(DEDPROJ)   !TRANSPOSE ASSUMING A HERMITEAN MATRIX
!     ==  DEPRO=DEPRO+PSI*CONJG(DEPROJ)
      CALL LIB$MATMULC8(NGL,NDIM*NBH,LMNX,PSI,DEDPROJ1,DEDPRO)
      IF(TINV) THEN
        ALLOCATE(PSIM(NGL))
        DO LMN=1,LMNX
          CALL PLANEWAVE$INVERTG(NGL,DEDPRO(1,LMN),PSIM)
          DO IG=1,NGL
            DEDPRO(IG,LMN)=0.5D0*(DEDPRO(IG,LMN)+PSIM(IG))
          ENDDO
        ENDDO
        DEALLOCATE(PSIM)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NB,NPRO,PSI &
     &                             ,PROPSI,PROFILE_ID)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATE PROJECTIONS                                               **
!     **                                                                      **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1998)***********
      USE WAVES_MODULE, ONLY : MAP_TYPE,GSET_TYPE
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      USE OPENACC
      USE WAVES_MODULE, ONLY: &
     &         ACC_PRO_CACHE_PRO &
     &        ,WAVES_ACC_PRO_CACHE_BUILD
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &         CPPAW_CUBLAS_ACC_SHOULD_USE_PROJECTION &
     &        ,CPPAW_CUBLAS_ACC_RESIDENCY_ENABLED &
     &        ,CPPAW_CUBLAS_ACC_PRO_EXPANSION_ENABLED &
     &        ,CPPAW_CUBLAS_ACC_PROJECTION_PRESENT &
     &        ,CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D
#ENDIF
      IMPLICIT NONE
      TYPE(MAP_TYPE) ,INTENT(IN) :: MAP
      TYPE(GSET_TYPE),INTENT(IN) :: GSET
      INTEGER(4)     ,INTENT(IN) :: NAT
      REAL(8)        ,INTENT(IN) :: R(3,NAT)
      INTEGER(4)     ,INTENT(IN) :: NGL        ! #(PLANE WAVES)
      INTEGER(4)     ,INTENT(IN) :: NDIM       ! #(SPINOR DIMENSIONS)
      INTEGER(4)     ,INTENT(IN) :: NB         ! #(STATES)
      INTEGER(4)     ,INTENT(IN) :: NPRO       ! #(PROJECTIONS)
      COMPLEX(8)     ,INTENT(IN) :: PSI(NGL,NDIM,NB) !WAVEFUNCTIONS
      COMPLEX(8)     ,INTENT(OUT):: PROPSI(NDIM,NB,NPRO) !
      CHARACTER(*)   ,INTENT(IN) :: PROFILE_ID
      REAL(8)        ,ALLOCATABLE:: GVEC(:,:)  !(3,NGL)
      COMPLEX(8)     ,ALLOCATABLE:: PRO(:,:)   !(NGL,LMNXX)
      COMPLEX(8)     ,ALLOCATABLE:: PROPSI1(:) !(LMNX,NDIM,NB)
      COMPLEX(8)     ,ALLOCATABLE:: PROPSIACC(:,:) !(LMNXX,NDIM*NB)
      COMPLEX(8)     ,ALLOCATABLE:: EIGR(:)    !(NGL)
      INTEGER(4)     ,ALLOCATABLE:: LOX(:)     !(LNXX)
      INTEGER(4)                 :: IPRO,IBPRO,IAT,LMN,IB,IDIM,II
      INTEGER(4)                 :: LMNXX,LNXX
      INTEGER(4)                 :: LNX,LMNX,ISP

      INTEGER(4)     ,SAVE       :: BLOCKSIZE=128
      INTEGER(4)                 :: BS,IGB
      COMPLEX(8)     ,ALLOCATABLE:: EIGRALL(:,:)    !(NGL,NAT)
      COMPLEX(8)     ,ALLOCATABLE:: PSITMP(:,:) !(LMNX,NDIM,NB)
      COMPLEX(8)     ,ALLOCATABLE:: PROTMP(:,:) !(LMNX,NDIM,NB)
      LOGICAL(4)                 :: TUSECUBLASPROJ
      LOGICAL(4)                 :: TUSERESIDENTPROJ
      LOGICAL(4)                 :: TUSEGPUEXPANDPRO
      LOGICAL(4)                 :: TSUPER
      REAL(8)                    :: PROJFLOPS
      REAL(8)                    :: GWEIGHT
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CHARACTER(32)               :: ACC_PRESENT_PSI
      CHARACTER(32)               :: ACC_COPY_PSI
#ENDIF
!     **************************************************************************
                                CALL TIMING$CLOCKON('WAVES_PROJECTIONS')
!
!     ==========================================================================
!     ==  CONSISTENCY CHECKS                                                  ==
!     ==========================================================================
      IF(NPRO.NE.MAP%NPRO) THEN
        CALL ERROR$MSG('INCONSISTENT NPRO')
        CALL ERROR$STOP('WAVES_PROJECTIONS')
      END IF
      IF(GSET%NGL.NE.NGL) THEN
        CALL ERROR$MSG('INCONSISTENT NGL')
        CALL ERROR$STOP('WAVES_PROJECTIONS')
      END IF
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      IF(ACC_IS_PRESENT(PROPSI)) THEN
!$ACC EXIT DATA DELETE(PROPSI(1:NDIM,1:NB,1:NPRO))
      END IF
#ENDIF
!
!     ==========================================================================
!     ==  COLLECT G-VECTORS                                                   ==
!     ==========================================================================
      ALLOCATE(GVEC(3,NGL))
      CALL PLANEWAVE$GETR8A('GVEC',3*GSET%NGL,GVEC)
!
!     ==========================================================================
!     ==  <PRO|PSI>                                                           ==
!     ==========================================================================
      LMNXX=MAXVAL(MAP%LMNX)
      LNXX=MAXVAL(MAP%LNX)
      ALLOCATE(PRO(NGL,LMNXX))
      ALLOCATE(EIGR(NGL))
      ALLOCATE(LOX(LNXX))
      TUSECUBLASPROJ=.FALSE.
      TUSERESIDENTPROJ=.FALSE.
      TUSEGPUEXPANDPRO=.FALSE.
      TSUPER=.TRUE.
      GWEIGHT=1.D0
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      TUSERESIDENTPROJ=CPPAW_CUBLAS_ACC_RESIDENCY_ENABLED()
      IF(TUSERESIDENTPROJ) THEN
        CALL PLANEWAVE$GETL4('SUPER',TSUPER)
        CALL PLANEWAVE$GETR8('GWEIGHT',GWEIGHT)
        TUSERESIDENTPROJ=.NOT.TSUPER
      END IF
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        PROJFLOPS=8.D0*REAL(NGL,KIND=8)*REAL(LMNX,KIND=8) &
     &           *REAL(NDIM*NB,KIND=8)
        TUSECUBLASPROJ=TUSECUBLASPROJ &
     &           .OR.CPPAW_CUBLAS_ACC_SHOULD_USE_PROJECTION(PROJFLOPS)
        TUSERESIDENTPROJ=TUSERESIDENTPROJ &
     &           .AND.CPPAW_CUBLAS_ACC_SHOULD_USE_PROJECTION(PROJFLOPS)
      ENDDO
#ENDIF
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      IF(TUSERESIDENTPROJ) THEN
!       == DEVICE-RESIDENT PROJECTION PATH =====================================
        ALLOCATE(PROPSIACC(LMNXX,NDIM*NB))
        IPRO=1
        TUSEGPUEXPANDPRO=CPPAW_CUBLAS_ACC_PRO_EXPANSION_ENABLED()
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        ACC_PRESENT_PSI='ACC_PRESENT_PROJ_'//TRIM(PROFILE_ID)//'_PSI'
        ACC_COPY_PSI='ACC_COPY_PROJ_'//TRIM(PROFILE_ID)//'_PSI_IN'
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &      (ACC_PRESENT_PSI,ACC_COPY_PSI,NGL,NDIM,NB,PSI)
        CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &      ('ACC_PRESENT_PROJ_PROPSI','ACC_COPY_PROJ_PROPSI_OUT' &
     &      ,NDIM,NB,NPRO,PROPSI)
#ENDIF
        IF(TUSEGPUEXPANDPRO) THEN
          CALL WAVES_ACC_PRO_CACHE_BUILD(MAP,GSET,NAT,R,NGL,LMNXX)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &        ('ACC_PRESENT_PROJ_PRO_CACHE','ACC_COPY_PROJ_PRO_CACHE_IN' &
     &        ,NGL,LMNXX,NAT,ACC_PRO_CACHE_PRO)
#ENDIF
!$ACC DATA PRESENT_OR_COPYIN(PSI(1:NGL,1:NDIM,1:NB)) &
!$ACC& PRESENT(ACC_PRO_CACHE_PRO(1:NGL,1:LMNXX,1:NAT)) &
!$ACC& COPYOUT(PROPSI(1:NDIM,1:NB,1:NPRO)) &
!$ACC& CREATE(PROPSIACC(1:LMNXX,1:NDIM*NB))
          DO IAT=1,MAP%NAT
            ISP=MAP%ISP(IAT)
            LNX=MAP%LNX(ISP)
            LMNX=MAP%LMNX(ISP)
            CALL CPPAW_CUBLAS_ACC_PROJECTION_PRESENT(NGL,NDIM,NB,LMNX &
     &          ,LMNXX,IPRO,NPRO,ACC_PRO_CACHE_PRO(:,:,IAT),PSI &
     &          ,GWEIGHT,PROPSIACC,PROPSI)
            IPRO=IPRO+LMNX
          ENDDO
!$ACC END DATA
        ELSE
!$ACC DATA PRESENT_OR_COPYIN(PSI(1:NGL,1:NDIM,1:NB)) &
!$ACC& COPYOUT(PROPSI(1:NDIM,1:NB,1:NPRO)) &
!$ACC& CREATE(PROPSIACC(1:LMNXX,1:NDIM*NB))
          DO IAT=1,MAP%NAT
            ISP=MAP%ISP(IAT)
            LNX=MAP%LNX(ISP)
            LMNX=MAP%LMNX(ISP)
            LOX=MAP%LOX(:,ISP)
            IBPRO=1+SUM(MAP%LNX(1:ISP-1))
            CALL PLANEWAVE$STRUCTUREFACTOR(R(:,IAT),NGL,EIGR)
            CALL WAVES_EXPANDPRO(LNX,LOX,LMNX,NGL,GVEC &
     &                          ,GSET%PRO(:,IBPRO:IBPRO+LNX-1),MAP%LMX &
     &                          ,GSET%YLM,EIGR,PRO)
            CALL CPPAW_CUBLAS_ACC_PROJECTION_PRESENT(NGL,NDIM,NB,LMNX &
     &              ,LMNXX,IPRO,NPRO,PRO,PSI,GWEIGHT,PROPSIACC,PROPSI)
            IPRO=IPRO+LMNX
          ENDDO
!$ACC END DATA
        END IF
      ELSE IF(.TRUE.) THEN
#ELSE
      IF(.TRUE.) THEN
#ENDIF
!       == UNBLOCKED ORIGINAL CODE SEGMENT =====================================
        ALLOCATE(PROPSI1(LMNXX*NDIM*NB))
        IPRO=1
#IF DEFINED(CPPVAR_CUBLAS_ACC)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
        IF(TUSECUBLASPROJ) THEN
          ACC_PRESENT_PSI='ACC_PRESENT_PROJ_'//TRIM(PROFILE_ID)//'_PSI'
          ACC_COPY_PSI='ACC_COPY_PROJ_'//TRIM(PROFILE_ID)//'_PSI_IN'
          CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &        (ACC_PRESENT_PSI,ACC_COPY_PSI,NGL,NDIM,NB,PSI)
        END IF
#ENDIF
#ENDIF
!$ACC DATA COPYIN(PSI(1:NGL,1:NDIM,1:NB)) IF(TUSECUBLASPROJ)
        DO IAT=1,MAP%NAT
          ISP=MAP%ISP(IAT)
          LNX=MAP%LNX(ISP)
          LMNX=MAP%LMNX(ISP)
          LOX=MAP%LOX(:,ISP)
          IBPRO=1+SUM(MAP%LNX(1:ISP-1))
          CALL PLANEWAVE$STRUCTUREFACTOR(R(:,IAT),NGL,EIGR)
          CALL WAVES_EXPANDPRO(LNX,LOX,LMNX,NGL,GVEC &
     &                        ,GSET%PRO(:,IBPRO:IBPRO+LNX-1),MAP%LMX &
     &                        ,GSET%YLM,EIGR,PRO)
          CALL PLANEWAVE$SCALARPRODUCT(' ',NGL,1,LMNX,PRO,NDIM*NB,PSI,PROPSI1)
!!$PRINT*,'IAT ',IAT,'R ',R
!!$PRINT*,'IAT ',IAT,'NGL ',NGL,' LNX=',LNX,' LOX=',LOX &
!!$       ,' LMNX=',LMNX,' LMX=',MAP%LMX
!!$PRINT*,'IAT ',IAT,'EIGR-1 ',MAXVAL(ABS(EIGR-(1.D0,0.D0)))
!!$PRINT*,'IAT ',IAT,'PROPSI1 ',PROPSI1 
          II=0
          DO IB=1,NB
            DO IDIM=1,NDIM
              DO LMN=1,LMNX
                II=II+1    ! II=(LMN,IDIM,IB)
                PROPSI(IDIM,IB,IPRO-1+LMN)=PROPSI1(II)
              ENDDO
            ENDDO
          ENDDO
          IPRO=IPRO+LMNX
        ENDDO
!$ACC END DATA
      ELSE
!       == BLOCKED CODE SEGMENT ================================================
        ALLOCATE(EIGRALL(NGL,MAP%NAT))
        DO IAT=1,MAP%NAT
          CALL PLANEWAVE$STRUCTUREFACTOR(R(1,IAT),NGL,EIGRALL(:,IAT))
        ENDDO
        ALLOCATE(PROPSI1(LMNXX*NB))
        PROPSI(:,:,:)=0
        DO IDIM=1,NDIM
          IF(.NOT.ALLOCATED(PSITMP))ALLOCATE(PSITMP(BLOCKSIZE,NB))
          IF(.NOT.ALLOCATED(PROTMP))ALLOCATE(PROTMP(BLOCKSIZE,LMNXX))
          DO IGB=1,NGL,BLOCKSIZE
            BS=MIN(NGL-IGB+1,BLOCKSIZE)
            IF(BS.NE.BLOCKSIZE)THEN
              DEALLOCATE(PSITMP)
              ALLOCATE(PSITMP(BS,NB))
              DEALLOCATE(PROTMP)
              ALLOCATE(PROTMP(BS,LMNXX))
            END IF
            PSITMP(1:BS,1:NB)=PSI(IGB:IGB+BS-1,IDIM,1:NB)

            IPRO=1
            DO IAT=1,MAP%NAT
              ISP=MAP%ISP(IAT)
              LNX=MAP%LNX(ISP)
              LMNX=MAP%LMNX(ISP)
              LOX=MAP%LOX(:,ISP)
              IBPRO=1+SUM(MAP%LNX(1:ISP-1))
              CALL WAVES_EXPANDPRO(LNX,LOX,LMNX,BS,GVEC(1:3,IGB:IGB+BS-1) &
     &                      ,GSET%PRO(IGB:IGB+BS-1,IBPRO:IBPRO+LNX-1) &
     &                      ,MAP%LMX,GSET%YLM(IGB:IGB+BS-1,:) &
     &                      ,EIGRALL(IGB:IGB+BS-1,IAT),PROTMP)
              CALL PLANEWAVE$SCALARPRODUCT(' ',BS,1,LMNX,PROTMP,NB,PSITMP &
     &                                                            ,PROPSI1)
              II=0
              DO IB=1,NB
                DO LMN=1,LMNX
                  II=II+1    ! II=(LMN,IDIM,IB)
                  PROPSI(IDIM,IB,IPRO-1+LMN)=PROPSI(IDIM,IB,IPRO-1+LMN) &
     &                                      +PROPSI1(II)
                ENDDO
              ENDDO
              IPRO=IPRO+LMNX
            ENDDO
          ENDDO
          DEALLOCATE(PSITMP)
          DEALLOCATE(PROTMP)
        ENDDO
        DEALLOCATE(EIGRALL)
!       == END OF BLOCKED CODE SEGMENT =========================================
      END IF
      DEALLOCATE(LOX)
      DEALLOCATE(PRO)
      IF(ALLOCATED(PROPSI1)) DEALLOCATE(PROPSI1)
      IF(ALLOCATED(PROPSIACC)) DEALLOCATE(PROPSIACC)
      DEALLOCATE(EIGR)
      DEALLOCATE(GVEC)

!!$DO IB=1,NB
!!$  PRINT*,'PSI '
!!$  DO LMN=1,NPRO
!!$     PRINT*,'PROPSI ',IB,LMN,PROPSI(1,IB,LMN)
!!$  ENDDO
!!$ENDDO
!!$STOP 'FORCED IN PAW_PROJECTIONS'
                               CALL TIMING$CLOCKOFF('WAVES_PROJECTIONS')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_ADDPRO(MAP,GSET,NAT,R,NGL,NDIM,NB,NPRO,PSI &
     &                        ,PROPSI,PROFILE_ID)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATE PROJECTIONS                                               **
!     **                                                                      **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1998)***********
      USE WAVES_MODULE, ONLY : MAP_TYPE,GSET_TYPE
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      USE WAVES_MODULE, ONLY: &
     &         ACC_PRO_CACHE_PRO &
     &        ,WAVES_ACC_PRO_CACHE_BUILD
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &         CPPAW_CUBLAS_ACC_SHOULD_USE_ADDPRODUCT &
     &        ,CPPAW_CUBLAS_ACC_ADDPRO_CACHE_CONTEXT_ENABLED &
     &        ,CPPAW_CUBLAS_ACC_ZGEMM_NN_PRESENT &
     &        ,CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &        ,CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D_INOUT
#ENDIF
      IMPLICIT NONE
      TYPE(MAP_TYPE) ,INTENT(IN) :: MAP
      TYPE(GSET_TYPE),INTENT(IN) :: GSET
      INTEGER(4)     ,INTENT(IN) :: NAT        ! #(ATOMS)
      REAL(8)        ,INTENT(IN) :: R(3,NAT)   ! ATOMIC POSITIONS
      INTEGER(4)     ,INTENT(IN) :: NGL        ! #(PLANE WAVES)
      INTEGER(4)     ,INTENT(IN) :: NDIM       ! #(SPINOR DIMENSIONS)
      INTEGER(4)     ,INTENT(IN) :: NB         ! #(STATES)
      INTEGER(4)     ,INTENT(IN) :: NPRO       ! #(PROJECTIONS)
      COMPLEX(8)     ,INTENT(INOUT) :: PSI(NGL,NDIM,NB)  !WAVEFUNCTIONS
      COMPLEX(8)     ,INTENT(IN) :: PROPSI(NDIM,NB,NPRO) !
      CHARACTER(*)   ,INTENT(IN) :: PROFILE_ID
      REAL(8)        ,ALLOCATABLE:: GVEC(:,:)  !(3,NGL)
      COMPLEX(8)     ,ALLOCATABLE:: PRO(:,:)   !(NGL,NDIM,LBLOCK)
      COMPLEX(8)     ,ALLOCATABLE:: PROPSI1(:) !(LMNX*NDIM*NB)
      COMPLEX(8)     ,ALLOCATABLE:: EIGR(:)    !(NGL)
      INTEGER(4)     ,ALLOCATABLE:: LOX(:)     !(LNXX)
      LOGICAL(4)    ,PARAMETER   :: TESSL=.TRUE.
      INTEGER(4)                 :: IPRO,IBPRO,IAT,LMN,IB,IDIM,II
      INTEGER(4)                 :: LMNXX,LNXX
      INTEGER(4)                 :: LNX,LMNX,ISP
!
      INTEGER(4)     ,SAVE       :: BLOCKSIZE=128
      INTEGER(4)                 :: BS,IGB
      COMPLEX(8)     ,ALLOCATABLE:: EIGRALL(:,:)    !(NGL,NAT)
      COMPLEX(8)     ,ALLOCATABLE:: PSITMP(:,:) !(LMNX,NDIM,NB)
      COMPLEX(8)     ,ALLOCATABLE:: PROTMP(:,:) !(LMNX,NDIM,NB)
      LOGICAL(4)                 :: TUSECUBLASADDPRO
      LOGICAL(4)                 :: TUSEGPUADDPROCACHE
      REAL(8)                    :: ADDPROFLOPS
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      REAL(8)                    :: ACCEL_T0
      REAL(8)                    :: ACCEL_T1
      REAL(8)                    :: ACCEL_TOTAL_T0
      REAL(8)                    :: ACCEL_TOTAL_T1
      CHARACTER(32)              :: ACC_PRESENT_PSI
      CHARACTER(32)              :: ACC_COPY_PSI_IN
      CHARACTER(32)              :: ACC_COPY_PSI_OUT
      CHARACTER(32)              :: ACC_PRESENT_PROPSI
      CHARACTER(32)              :: ACC_COPY_PROPSI
      CHARACTER(32)              :: ACC_PRESENT_CACHE
      CHARACTER(32)              :: ACC_COPY_CACHE
#ENDIF
!     **************************************************************************
                                 CALL TRACE$PUSH('WAVES_ADDPRO')
                                 CALL TIMING$CLOCKON('WAVES_ADDPRO')
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T0)
#ENDIF
!
!     ==========================================================================
!     ==  CONSISTENCY CHECKS                                                  ==
!     ==========================================================================
      IF(NPRO.NE.MAP%NPRO) THEN
        CALL ERROR$MSG('INCONSISTENT NPRO')
        CALL ERROR$STOP('WAVES_ADDPRO')
      END IF
      IF(GSET%NGL.NE.NGL) THEN
        CALL ERROR$MSG('INCONSISTENT NGL')
        CALL ERROR$STOP('WAVES_ADDPRO')
      END IF
!
!     ==========================================================================
!     ==  COLLECT G-VECTORS                                                   ==
!     ==========================================================================
      ALLOCATE(GVEC(3,NGL))
      CALL PLANEWAVE$GETR8A('GVEC',3*GSET%NGL,GVEC)
!
!     ==========================================================================
!     ==  COLLECT G-VECTORS                                                   ==
!     ==========================================================================
      LMNXX=MAXVAL(MAP%LMNX)
      LNXX=MAXVAL(MAP%LNX)
      ALLOCATE(PRO(NGL,LMNXX))
      ALLOCATE(PROPSI1(LMNXX*NDIM*NB))
      ALLOCATE(EIGR(NGL))
      ALLOCATE(LOX(LNXX))
      TUSEGPUADDPROCACHE=.FALSE.
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      TUSECUBLASADDPRO=.TRUE.
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        ADDPROFLOPS=8.D0*REAL(NGL,KIND=8)*REAL(LMNX,KIND=8) &
     &             *REAL(NDIM*NB,KIND=8)
        TUSECUBLASADDPRO=TUSECUBLASADDPRO &
     &             .AND.CPPAW_CUBLAS_ACC_SHOULD_USE_ADDPRODUCT(ADDPROFLOPS)
      ENDDO
      IF(TUSECUBLASADDPRO) THEN
#ELSE
      IF(.FALSE.) THEN
#ENDIF
!       == UNBLOCKED ORIGINAL CODE SEGMENT =====================================
        IPRO=1
#IF DEFINED(CPPVAR_CUBLAS_ACC)
        TUSEGPUADDPROCACHE= &
     &       CPPAW_CUBLAS_ACC_ADDPRO_CACHE_CONTEXT_ENABLED(PROFILE_ID)
#ENDIF
        IF(TUSEGPUADDPROCACHE) THEN
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          CALL WAVES_ACC_PRO_CACHE_BUILD(MAP,GSET,NAT,R,NGL,LMNXX)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
          ACC_PRESENT_PSI='ACC_PRESENT_ADDPRO_'//TRIM(PROFILE_ID)//'_PSI'
          ACC_COPY_PSI_IN='ACC_COPY_ADDPRO_'//TRIM(PROFILE_ID)//'_PSI_IN'
          ACC_COPY_PSI_OUT='ACC_COPY_ADDPRO_'//TRIM(PROFILE_ID)//'_PSI_OUT'
          ACC_PRESENT_PROPSI='ACC_PRESENT_ADDPRO_'//TRIM(PROFILE_ID) &
     &        //'_PROPSI'
          ACC_COPY_PROPSI='ACC_COPY_ADDPRO_'//TRIM(PROFILE_ID) &
     &        //'_PROPSI_IN'
          ACC_PRESENT_CACHE='ACC_PRESENT_ADDPRO_'//TRIM(PROFILE_ID) &
     &        //'_CACHE'
          ACC_COPY_CACHE='ACC_COPY_ADDPRO_'//TRIM(PROFILE_ID) &
     &        //'_CACHE_IN'
          CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D_INOUT &
     &        (ACC_PRESENT_PSI,ACC_COPY_PSI_IN,ACC_COPY_PSI_OUT &
     &        ,NGL,NDIM,NB,PSI)
          CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &        (ACC_PRESENT_PROPSI,ACC_COPY_PROPSI &
     &        ,NDIM,NB,NPRO,PROPSI)
          CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &        (ACC_PRESENT_CACHE,ACC_COPY_CACHE,NGL,LMNXX,NAT &
     &        ,ACC_PRO_CACHE_PRO)
#ENDIF
!$ACC DATA PRESENT_OR_COPY(PSI(1:NGL,1:NDIM,1:NB)) &
!$ACC& PRESENT_OR_COPYIN(PROPSI(1:NDIM,1:NB,1:NPRO)) &
!$ACC& PRESENT(ACC_PRO_CACHE_PRO(1:NGL,1:LMNXX,1:NAT)) &
!$ACC& CREATE(PROPSI1(1:LMNXX*NDIM*NB))
          DO IAT=1,MAP%NAT
            ISP=MAP%ISP(IAT)
            LMNX=MAP%LMNX(ISP)
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(II) PRESENT(PROPSI,PROPSI1)
            DO IB=1,NB
              DO IDIM=1,NDIM
                DO LMN=1,LMNX
                  II=LMN+(IDIM-1)*LMNX+(IB-1)*LMNX*NDIM
                  PROPSI1(II)=PROPSI(IDIM,IB,IPRO-1+LMN)
                ENDDO
              ENDDO
            ENDDO
!$ACC END PARALLEL LOOP
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL ACCELPROFILE$NOW(ACCEL_T0)
#ENDIF
            CALL CPPAW_CUBLAS_ACC_ZGEMM_NN_PRESENT(NGL,LMNX,NDIM*NB &
     &          ,ACC_PRO_CACHE_PRO(:,:,IAT),PROPSI1,PSI)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL ACCELPROFILE$NOW(ACCEL_T1)
            CALL ACCELPROFILE$ADD('CUBLAS_ZGEMM_ADDPRO_CACHE' &
     &          ,INT(NGL,KIND=8),INT(LMNX,KIND=8),INT(NDIM*NB,KIND=8) &
     &          ,0_8,8.D0*REAL(NGL,KIND=8)*REAL(LMNX,KIND=8) &
     &          *REAL(NDIM*NB,KIND=8) &
     &          ,16.D0*(REAL(LMNX,KIND=8)*REAL(NDIM*NB,KIND=8)) &
     &          ,ACCEL_T1-ACCEL_T0)
#ENDIF
            IPRO=IPRO+LMNX
          ENDDO
!$ACC END DATA
#ENDIF
        ELSE
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          ACC_PRESENT_PSI='ACC_PRESENT_ADDPRO_'//TRIM(PROFILE_ID)//'_PSI'
          ACC_COPY_PSI_IN='ACC_COPY_ADDPRO_'//TRIM(PROFILE_ID)//'_PSI_IN'
          ACC_COPY_PSI_OUT='ACC_COPY_ADDPRO_'//TRIM(PROFILE_ID)//'_PSI_OUT'
          CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D_INOUT &
     &        (ACC_PRESENT_PSI,ACC_COPY_PSI_IN,ACC_COPY_PSI_OUT &
     &        ,NGL,NDIM,NB,PSI)
#ENDIF
#ENDIF
!$ACC DATA COPY(PSI(1:NGL,1:NDIM,1:NB))
        DO IAT=1,MAP%NAT
          ISP=MAP%ISP(IAT)
          LNX=MAP%LNX(ISP)
          LMNX=MAP%LMNX(ISP)
          LOX=MAP%LOX(:,ISP)
          IBPRO=1+SUM(MAP%LNX(1:ISP-1))
          CALL PLANEWAVE$STRUCTUREFACTOR(R(1,IAT),NGL,EIGR)
          CALL WAVES_EXPANDPRO(LNX,LOX,LMNX,NGL,GVEC &
       &              ,GSET%PRO(:,IBPRO:IBPRO+LNX-1),MAP%LMX,GSET%YLM,EIGR,PRO)
!    
!         ======================================================================
!         ==  |PSI>=|PSI>+|PRO>PROPSI                                         ==
!         ======================================================================
          II=0
          DO IB=1,NB
            DO IDIM=1,NDIM
              DO LMN=1,LMNX
                II=II+1    ! II=(LMN,IDIM,IB)
                PROPSI1(II)=PROPSI(IDIM,IB,IPRO-1+LMN)
              ENDDO
            ENDDO
          ENDDO
!         == PSI=PSI+PRO*PROPSI1
          CALL LIB$ADDPRODUCTC8(.FALSE.,NGL,LMNX,NDIM*NB,PRO,PROPSI1,PSI)
          IPRO=IPRO+LMNX
        ENDDO
!$ACC END DATA
        END IF
!
      ELSE
!       == BLOCKED CODE SEGMENT ================================================
        ALLOCATE(EIGRALL(NGL,MAP%NAT))
        DO IAT=1,MAP%NAT
          CALL PLANEWAVE$STRUCTUREFACTOR(R(1,IAT),NGL,EIGRALL(:,IAT))
        ENDDO
       
        DO IDIM=1,NDIM
          IF(.NOT.ALLOCATED(PSITMP))ALLOCATE(PSITMP(BLOCKSIZE,NB))
          IF(.NOT.ALLOCATED(PROTMP))ALLOCATE(PROTMP(BLOCKSIZE,LMNXX))
          DO IGB=1,NGL,BLOCKSIZE
            BS=MIN(NGL-IGB+1,BLOCKSIZE)
            IF(BS.NE.BLOCKSIZE)THEN
              DEALLOCATE(PSITMP)
              ALLOCATE(PSITMP(BS,NB))
              DEALLOCATE(PROTMP)
              ALLOCATE(PROTMP(BS,LMNXX))
            ENDIF
            PSITMP(1:BS,1:NB)=PSI(IGB:IGB+BS-1,IDIM,1:NB)
            
            IPRO=1
            DO IAT=1,MAP%NAT
              ISP=MAP%ISP(IAT)
              LNX=MAP%LNX(ISP)
              LMNX=MAP%LMNX(ISP)
              LOX=MAP%LOX(:,ISP)
              IBPRO=1+SUM(MAP%LNX(1:ISP-1))
              CALL WAVES_EXPANDPRO(LNX,LOX,LMNX,BS,GVEC(1:3,IGB:IGB+BS-1) &
           &                      ,GSET%PRO(IGB:IGB+BS-1,IBPRO:IBPRO+LNX-1) &
           &                      ,MAP%LMX,GSET%YLM(IGB:IGB+BS-1,:) &
           &                      ,EIGRALL(IGB:IGB+BS-1,IAT),PROTMP)
!
!             ==================================================================
!             ==  |PSI>=|PSI>+|PRO>PROPSI                                     ==
!             ==================================================================
              II=0
              DO IB=1,NB
                !DO IDIM=1,NDIM
                  DO LMN=1,LMNX
                    II=II+1    ! II=(LMN,IDIM,IB)
                    PROPSI1(II)=PROPSI(IDIM,IB,IPRO-1+LMN)
                  ENDDO
                !ENDDO
              ENDDO
!             == PSI=PSI+PRO*PROPSI1 ===========================================
              CALL LIB$ADDPRODUCTC8(.FALSE.,BS,LMNX,NB,PROTMP,PROPSI1,PSITMP)
              IPRO=IPRO+LMNX
            ENDDO
            PSI(IGB:IGB+BS-1,IDIM,1:NB)=PSITMP(1:BS,1:NB)
          ENDDO
          DEALLOCATE(PROTMP)
          DEALLOCATE(PSITMP)
        ENDDO
        DEALLOCATE(EIGRALL)
!       == END OF BLOCKED CODE SEGMENT =========================================
      END IF  
      DEALLOCATE(LOX)
      DEALLOCATE(PRO)
      DEALLOCATE(GVEC)
      DEALLOCATE(EIGR)
      DEALLOCATE(PROPSI1)
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
      CALL ACCELPROFILE$NOW(ACCEL_TOTAL_T1)
      CALL ACCELPROFILE$ADD('PAW_ADDPRO_TOTAL' &
     &    ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NB,KIND=8) &
     &    ,INT(NPRO,KIND=8),0.D0,0.D0,ACCEL_TOTAL_T1-ACCEL_TOTAL_T0)
#ENDIF
                                 CALL TIMING$CLOCKOFF('WAVES_ADDPRO')
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_EXPANDPRO(LNX,LOX,LMNX,NGL,GVEC &
     &                          ,BAREPRO,LMX,YLM,EIGR,PRO)
!     **************************************************************************
!     **                                                                      **
!     **  TAKES THE BARE PROJECTOR FUNCTIONS AND CALCULATES FULL              **
!     **  PROJECTOR FUNCTIONS WITH STRUCTURE FACTOR AND I**L                  **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1998)***********
      IMPLICIT NONE
      INTEGER(4)    ,INTENT(IN) :: LNX       ! #(PROJECTORS PER ATOM)
      INTEGER(4)    ,INTENT(IN) :: LMNX      ! #(PROJECTORS PER ATOM)
      INTEGER(4)    ,INTENT(IN) :: LOX(LNX)  ! ANGULAR MOMENTA
      INTEGER(4)    ,INTENT(IN) :: NGL           ! #(G-VECTORS(LOCAL))
      REAL(8)       ,INTENT(IN) :: GVEC(3,NGL)   ! G-VECTORS
      REAL(8)       ,INTENT(IN) :: BAREPRO(NGL,LNX)  ! PROJECTOR FUNCTIONS
      INTEGER(4)    ,INTENT(IN) :: LMX           ! X#(ANGULAR MOMENTA IN YLM))
      REAL(8)       ,INTENT(IN) :: YLM(NGL,LMX)  ! CI**(-L)*SPHERICAL HARMONICS
      COMPLEX(8)    ,INTENT(IN) :: EIGR(NGL)     ! STRUCTURE FACTORS
      COMPLEX(8)    ,INTENT(OUT):: PRO(NGL,LMNX) !PROJECTOR FUNCTION
      COMPLEX(8)    ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)                :: LMN,LN,L,LM,IM,IG
      REAL(8)                   :: SVAR
      COMPLEX(8)                :: CSVAR
!     **************************************************************************
      LMN=0
      DO LN=1,LNX
        L=LOX(LN)
        LM=L**2
        DO IM=1,2*L+1
          LMN=LMN+1
          LM=LM+1
          CSVAR=(-CI)**L
          DO IG=1,NGL
            SVAR=YLM(IG,LM)*BAREPRO(IG,LN)
            PRO(IG,LMN)=SVAR*CSVAR*EIGR(IG)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_EXPANDPRO_ACC(LNX,LOX,LMNX,NGL,BAREPRO,LMX,YLM &
     &                              ,EIGR,PRO)
!     **************************************************************************
!     **  OpenACC projector expansion for the resident cuBLAS projection path.  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)    ,INTENT(IN) :: LNX
      INTEGER(4)    ,INTENT(IN) :: LMNX
      INTEGER(4)    ,INTENT(IN) :: LOX(LNX)
      INTEGER(4)    ,INTENT(IN) :: NGL
      REAL(8)       ,INTENT(IN) :: BAREPRO(NGL,LNX)
      INTEGER(4)    ,INTENT(IN) :: LMX
      REAL(8)       ,INTENT(IN) :: YLM(NGL,LMX)
      COMPLEX(8)    ,INTENT(IN) :: EIGR(NGL)
      COMPLEX(8)    ,INTENT(OUT):: PRO(NGL,LMNX)
      INTEGER(4)                :: LN,L,LM,LMN,LMN0,IM,IG,ILN
      REAL(8)                   :: SVAR
      COMPLEX(8)                :: CSVAR
!     **************************************************************************
!$ACC DATA PRESENT_OR_COPYIN(BAREPRO(1:NGL,1:LNX),YLM(1:NGL,1:LMX) &
!$ACC&       ,EIGR(1:NGL)) PRESENT(PRO(1:NGL,1:LMNX)) COPYIN(LOX(1:LNX))
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(LN,L,LM,LMN0,IM,ILN,SVAR,CSVAR)
      DO LMN=1,LMNX
        DO IG=1,NGL
          LMN0=0
          LN=1
          DO ILN=1,LNX
            IF(LMN.LE.LMN0+2*LOX(ILN)+1) THEN
              LN=ILN
              EXIT
            END IF
            LMN0=LMN0+2*LOX(ILN)+1
          ENDDO
          L=LOX(LN)
          IM=LMN-LMN0
          LM=L**2+IM
          SELECT CASE(MOD(L,4))
          CASE(0)
            CSVAR=(1.D0,0.D0)
          CASE(1)
            CSVAR=(0.D0,-1.D0)
          CASE(2)
            CSVAR=(-1.D0,0.D0)
          CASE DEFAULT
            CSVAR=(0.D0,1.D0)
          END SELECT
          SVAR=YLM(IG,LM)*BAREPRO(IG,LN)
          PRO(IG,LMN)=SVAR*CSVAR*EIGR(IG)
        ENDDO
      ENDDO
!$ACC END PARALLEL LOOP
!$ACC END DATA
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_PROFORCE(LNX,LMNX,LOX,NGL,GWEIGHT,GVEC,GIJ &
     &          ,BAREPRO,DBAREPRO,LMX,YLM,SYLM,EIGR,DEDPRO,FORCE,TSTRESS,STRESS)
!     **************************************************************************
!     **                                                                      **
!     **  TAKES THE BARE PROJECTOR FUNCTIONS AND CALCULATES FULL              **
!     **  PROJECTOR FUNCTIONS WITH STRUCTURE FACTOR AND I**L                  **
!     **                                                                      **
!     **  FOR NDIM.NE.3 FIRST DERIVATIVES ARE CALCULATED                      **
!     **  FOR NDIM.NE.6 SECOND DERIVATIVES ARE CALCULATED                     **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1998)***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LNX           ! #(BARE PROJECTORSS)
      INTEGER(4),INTENT(IN) :: LMNX          ! #(PROJECTORSS)
      INTEGER(4),INTENT(IN) :: LOX(LNX)      ! ANGULAR MOMENTA
      INTEGER(4),INTENT(IN) :: NGL           ! #(G-VECTORS(LOCAL))
      REAL(8)   ,INTENT(IN) :: GWEIGHT
      REAL(8)   ,INTENT(IN) :: GVEC(3,NGL)   ! G-VECTORS
      REAL(8)   ,INTENT(IN) :: GIJ(6,NGL)    ! GI*GJ/G**2
      REAL(8)   ,INTENT(IN) :: BAREPRO(NGL,LNX) ! BARE PROJECTOR FUNCTIONS
      REAL(8)   ,INTENT(IN) :: DBAREPRO(NGL,LNX) ! |G|DBAREPRO/DG
      INTEGER(4),INTENT(IN) :: LMX           ! X#(ANGULAR MOMENTA IN YLM))
      REAL(8)   ,INTENT(IN) :: YLM(NGL,LMX)  ! SPHERICAL HARMONICS
      REAL(8)   ,INTENT(IN) :: SYLM(NGL,LMX,6) ! STRAINED SPHERICAL HARMONICS
      COMPLEX(8),INTENT(IN) :: EIGR(NGL)     ! STRUCTURE FACTORS
      COMPLEX(8),INTENT(IN) :: DEDPRO(NGL,LMNX) !PROJECTOR FUNCTION
      REAL(8)   ,INTENT(OUT):: FORCE(3)
      LOGICAL(4),INTENT(IN) :: TSTRESS       ! SWITCH FOR STRESS CALC.
      REAL(8)   ,INTENT(OUT):: STRESS(3,3)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: CSVAR,CSVAR1,CSVAR2,CIL
      COMPLEX(8)            :: CWORK1(NGL)
      COMPLEX(8)            :: CWORK2(NGL)
      REAL(8)               :: SVAR,SVAR1,SVAR2
      REAL(8)               :: F(3)
      REAL(8)               :: S(6),SDIAG
      INTEGER(4)            :: LN,LMN,LM,L,M
      INTEGER(4)            :: IG,I,J,IJ
!     **************************************************************************
      F(:)=0.D0
      STRESS=0.D0
      S(:)=0.D0
      SDIAG=0.D0
!   
!     ==========================================================================
!     ==  START THE LOOP                                                      ==
!     ==========================================================================
      LMN=0
      DO LN=1,LNX
        L=LOX(LN)
        CIL=(-CI)**L
        IF(TSTRESS) THEN
          DO IG=1,NGL
            CSVAR=CIL*EIGR(IG)
            CWORK1(IG)= BAREPRO(IG,LN)*CSVAR
            CWORK2(IG)=DBAREPRO(IG,LN)*CSVAR
          ENDDO
        ELSE
          DO IG=1,NGL
            CWORK1(IG)=CI*BAREPRO(IG,LN)*CIL*EIGR(IG)
          ENDDO
        END IF
!
        LM=L**2
        DO M=1,2*L+1
          LM=LM+1
          LMN=LMN+1
!           
!         ======================================================================
!         ==  COMPOSE PROJECTOR FUNCTIONS                                     ==
!         ======================================================================
!         ==  MULTIPLY PROJECTOR WITH STRUCTURE FACTOR ===========
          IF(TSTRESS) THEN
            DO IG=1,NGL
!             == THE USE OF CONJG IN THE FOLLOWING LINE IS UNCLEAR
!             == PREVIOUSLY CONJG WAS USED WHENEVER TINV=.TRUE.
              CSVAR =CONJG(DEDPRO(IG,LMN))
              CSVAR1=CSVAR*CWORK1(IG)
              CSVAR2=CSVAR*CWORK2(IG)
              SVAR  =REAL(CI*CSVAR1,KIND=8)*YLM(IG,LM)
              DO I=1,3
                F(I)=F(I)+SVAR*GVEC(I,IG)
              ENDDO
              SVAR1=REAL(CSVAR1,KIND=8)
              SVAR2=REAL(CSVAR2,KIND=8)*YLM(IG,LM)
              DO IJ=1,6
                S(IJ)=S(IJ)-SVAR1*SYLM(IG,LM,IJ)-SVAR2*GIJ(IJ,IG)
              ENDDO
              SDIAG=SDIAG-SVAR1*YLM(IG,LM)
            ENDDO
          ELSE
            DO IG=1,NGL
              SVAR=REAL(CWORK1(IG)*CONJG(DEDPRO(IG,LMN)),KIND=8)*YLM(IG,LM)
              DO I=1,3
                F(I)=F(I)+SVAR*GVEC(I,IG)
              ENDDO
            ENDDO
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == <P|DE/D<P| + DE/D|P> |P>                                             ==
!     ==========================================================================
      DO I=1,3
        FORCE(I)=GWEIGHT*2.D0*F(I)
      ENDDO
      IF(TSTRESS) THEN
        IJ=0
        DO I=1,3
          DO J=I,3
            IJ=IJ+1
            IF(I.EQ.J) THEN
              STRESS(I,J)=GWEIGHT*(2.D0*S(IJ)+SDIAG)
            ELSE
              STRESS(I,J)=GWEIGHT*2.D0*S(IJ)
              STRESS(J,I)=STRESS(I,J)
            END IF
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_UPDATEGSET()
!     **************************************************************************
!     **                                                                      **
!     **  THIS ROUTINE MUST BE CALLED INITIALLY AND AFTER EACH                **
!     **  CHANGE OF THE CELL-SHAPE.                                           **
!     **                                                                      **
!     **    GSET%PRO PROJECTOR FUNCTIONS  GRADIENT                            **
!     **    GSET%DPRO                                                         **
!     **    GSET%YLM                                                          **
!     **    GSET%SYLM                                                         **
!     **    GSET%BUCKET                                                      **
!     **    GSET%DBUCKET                                                      **
!     **                                                                      **
!     *******************************************P.E. BLOECHL, (1999)***********
      USE WAVES_MODULE, ONLY : GSET_TYPE,GSET,MAP,TBUCKET,EPWPSI0,D2EPWPSI &
     &                         ,NKPTL,WAVES_SELECTWV 
      IMPLICIT NONE
      REAL(8)        ,ALLOCATABLE   :: G2(:)
      REAL(8)        ,ALLOCATABLE   :: GVEC(:,:)
      REAL(8)        ,ALLOCATABLE   :: YLM(:)
      INTEGER(4)                    :: IG,ISP,IND,LN,IKPT
      INTEGER(4)                    :: NBAREPRO,NGL
      REAL(8)                       :: RBAS(3,3),GBAS(3,3),CELLVOL
      INTEGER(4)                    :: LMX
      LOGICAL(4)                    :: TSTRESS
      REAL(8)                       :: MING2,MING,DG,ABSG
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES_UPDATEGSET')
!
!     ==========================================================================
!     ==  UPDATE RADIAL PROJECTOR FUNCTIONS                                   ==
!     ==========================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      CALL CELL$GETL4('ON',TSTRESS)
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL POTENTIAL$SETR8A('RBAS',9,RBAS)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      LMX=MAP%LMX
      DO IKPT=1,NKPTL
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL=GSET%NGL
        NBAREPRO=MAP%NBAREPRO
!
!       ========================================================================
!       == UPDATE PLANEWAVE OBJECT                                            ==
!       ========================================================================
        CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
!
!       ========================================================================
!       ==  EVALUATE PROJECTOR FUNCTIONS IN G-SPACE                           ==
!       ========================================================================
        IF(.NOT.ASSOCIATED(GSET%PRO)) ALLOCATE(GSET%PRO(NGL,NBAREPRO))
        IF(.NOT.ASSOCIATED(GSET%DPRO))ALLOCATE(GSET%DPRO(NGL,NBAREPRO))
        ALLOCATE(G2(NGL))
        CALL PLANEWAVE$GETR8A('G2',NGL,G2)
        IND=0
        DO ISP=1,MAP%NSP
          CALL SETUP$ISELECT(ISP)
          DO LN=1,MAP%LNX(ISP)
            IND=IND+1
            CALL SETUP$GETFOFG('PRO',.FALSE.,LN,NGL,G2,CELLVOL,GSET%PRO(:,IND))
            CALL SETUP$GETFOFG('PRO',.TRUE.,LN,NGL,G2,CELLVOL,GSET%DPRO(:,IND))
          ENDDO
          CALL SETUP$UNSELECT()
        ENDDO
!
!       ========================================================================
!       ==  UPDATE SPHERICAL HARMONICS                                        ==
!       ========================================================================
        IF(.NOT.ASSOCIATED(GSET%YLM)) ALLOCATE(GSET%YLM(NGL,LMX))
        IF(.NOT.ASSOCIATED(GSET%SYLM))ALLOCATE(GSET%SYLM(NGL,LMX,6))
        ALLOCATE(GVEC(3,NGL))
        CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
        ALLOCATE(YLM(LMX))
        DO IG=1,NGL
          CALL GETYLM(LMX,GVEC(1,IG),YLM)
          GSET%YLM(IG,:)=YLM(:) 
        ENDDO        
        DEALLOCATE(YLM)
!       == NOW THE STRAINED SPHERICAL HARMONICS ================================
        CALL WAVES_STRAINEDYLM(NGL,LMX,GVEC,GSET%YLM,GSET%SYLM)
        DEALLOCATE(GVEC)
!
!       ========================================================================
!       ==  EVALUATE BUCKET POTENTIAL                                         ==
!       ========================================================================
        IF(.NOT.ASSOCIATED(GSET%BUCKET)) ALLOCATE(GSET%BUCKET(NGL))
        IF(.NOT.ASSOCIATED(GSET%DBUCKET)) ALLOCATE(GSET%DBUCKET(NGL))
        IF(TBUCKET) THEN
          MING2=2.D0*EPWPSI0
          MING=SQRT(MING2)
          DO IG=1,NGL
            IF(G2(IG).LT.MING2) THEN
              GSET%BUCKET(IG)=0.D0
              GSET%DBUCKET(IG)=0.D0
            ELSE
              ABSG=SQRT(G2(IG))
              DG=ABSG-MING
              GSET%BUCKET(IG) =     D2EPWPSI*DG**2
              GSET%DBUCKET(IG)=2.D0*D2EPWPSI*DG/ABSG
            END IF 
          ENDDO
        ELSE
          GSET%BUCKET(:)=0.D0
          GSET%DBUCKET(:)=0.D0
        END IF
        DEALLOCATE(G2)
      ENDDO
!
!     ==========================================================================
!     ==  EVALUATE WAVE FUNCTION MASS                                         ==
!     ==========================================================================
      CALL WAVES_PSIMASS
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_STRAINEDYLM(NGL,LMX,GVEC,YLM,SYLM)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES STRAINED SPHERICAL HARMONICS ON THE G-SPACE GRID  **
!     **                                                              **
!     *******************************************P.E. BLOECHL, (1999)***
      USE WAVES_MODULE, ONLY : GSET_TYPE,MAP_TYPE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN)    :: NGL  ! #(G-VECTORS)
      INTEGER(4)     ,INTENT(IN)    :: LMX ! X#(SPHERICAL HARMONICS
      REAL(8)        ,INTENT(IN)    :: GVEC(3,NGL)
      REAL(8)        ,INTENT(IN)    :: YLM(NGL,LMX)
      REAL(8)        ,INTENT(OUT)   :: SYLM(NGL,LMX,6)
      INTEGER(4)                    :: IG,IND,I,J,LM1,LM2,IJ,L
      REAL(8)                       :: CC(3,3)
      REAL(8)                       :: SVAR
      REAL(8)                       :: GIJ(NGL,6)
      REAL(8)        ,PARAMETER     :: RSMALL=1.D-20
!     ******************************************************************
!
!     ==========================================================================
!     ==  PREPARE GI*GJ/G2                                                    ==
!     ==========================================================================
      DO IG=1,NGL
        SVAR=1.D0/(GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2+RSMALL)
        IJ=0
        DO I=1,3
          DO J=I,3
            IJ=IJ+1
            GIJ(IG,IJ)=SVAR*GVEC(I,IG)*GVEC(J,IG)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  NOW STRAIN THE SPHERICAL HARMONICS                                  ==
!     ==========================================================================
      SYLM=0.D0
      DO LM1=1,LMX
!
!       ========================================================================
!       ==                                                                    ==
!       ========================================================================
        DO IND=1,10
          CALL SPHERICAL$CCMAT0(LM1,IND,LM2,CC)
          IF(LM2.EQ.0) EXIT
          IJ=0
          DO I=1,3
            DO J=I,3
              IJ=IJ+1
              DO IG=1,NGL
                SYLM(IG,LM1,IJ)=SYLM(IG,LM1,IJ)+CC(I,J)*YLM(IG,LM2)
              ENDDO
            ENDDO
          ENDDO
        ENDDO        
!
!       ========================================================================
!       ==                                                                    ==
!       ========================================================================
        DO IND=1,10
          CALL SPHERICAL$CCMATM(LM1,IND,LM2,CC)
          IF(LM2.EQ.0) EXIT
          IJ=0
          DO I=1,3
            DO J=I,3
              IJ=IJ+1
              DO IG=1,NGL
                SYLM(IG,LM1,IJ)=SYLM(IG,LM1,IJ)+CC(I,J)*YLM(IG,LM2)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  NOW THE TERM -L*GI*GJ/G**2                                        ==
!       ========================================================================
        L=INT(SQRT(REAL(LM1-1,KIND=8))+1.D-5)
        SVAR=-REAL(L,KIND=8)
        DO IJ=1,6
          DO IG=1,NGL
            SYLM(IG,LM1,IJ)=SYLM(IG,LM1,IJ)+SVAR*YLM(IG,LM1)*GIJ(IG,IJ)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$STOP()
!     **************************************************************************
!     **  SETS THE VELOCITY OF WAVE FUNCTIONS TO ZERO                         **
!     **************************************************************************
      USE WAVES_MODULE, ONLY : NKPTL &
     &                        ,NSPIN &
     &                        ,WAVES_SELECTWV &
     &                        ,THIS &
     &                        ,WAVEEKIN1 
      IMPLICIT NONE
      INTEGER(4)           :: IKPT
      INTEGER(4)           :: ISPIN
!     **************************************************************************
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
!!$PRINT*,'WAVES$STOP ',IKPT,ISPIN,ASSOCIATED(THIS%PSIM),ASSOCIATED(THIS%PSI0)
!!$PRINT*,'WAVES$STOP PSI0:',SHAPE(THIS%PSI0)
!!$PRINT*,'WAVES$STOP PSIM:',SHAPE(THIS%PSIM)
          THIS%PSIM=THIS%PSI0
!!$PRINT*,'WAVES$STOP PASS1'
          THIS%PSIM(:,:,:)=THIS%PSI0(:,:,:)
!!$PRINT*,'WAVES$STOP PASS'

          IF(ASSOCIATED(THIS%RLAMM))DEALLOCATE(THIS%RLAMM)
          IF(ASSOCIATED(THIS%RLAM2M))DEALLOCATE(THIS%RLAM2M)
          IF(ASSOCIATED(THIS%RLAM3M))DEALLOCATE(THIS%RLAM3M)
        ENDDO
      ENDDO
      WAVEEKIN1=0.D0
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES$PROPAGATE()
!     **************************************************************************
!     **  PROPAGATES WAVE FUNCTIONS                                           **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **    PARAMETERS DELT,EMASS,EMASSCG2,ANNEE AND TSTOP MUST BE SET        **
!     **                                                                      **
!     **************************************************************************
      USE WAVES_MODULE
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      USE CPPAW_CUBLAS_ACC_MODULE, ONLY: &
     &       CPPAW_CUBLAS_ACC_PSIM_PROPAGATE_ENABLED &
     &      ,CPPAW_CUBLAS_ACC_PSIM_PHASE_RESIDENCY_ENABLED &
     &      ,CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &      ,CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D_INOUT
#ENDIF
      IMPLICIT NONE
      INTEGER(4)            :: IKPT,ISPIN,IG,IDIM,IB
      INTEGER(4)            :: NGL
      INTEGER(4)            :: NBH
      REAL(8)   ,ALLOCATABLE:: ARR1(:)
      REAL(8)   ,ALLOCATABLE:: ARR2(:)
      REAL(8)   ,ALLOCATABLE:: ARR3(:)
      REAL(8)               :: SVAR1,SVAR2,SVAR3
      LOGICAL(4)            :: TSTRESS
      REAL(8)               :: FRICMAT(3,3)
      REAL(8)  ,ALLOCATABLE :: GVEC(:,:)
      REAL(8)  ,ALLOCATABLE :: ANNEEVEC(:)
#IF DEFINED(CPPVAR_CUBLAS_ACC)
      LOGICAL(4)            :: TRESIDENTPSIM
      LOGICAL(4)            :: TPHASERESIDENTPSIM
#ENDIF
!     **************************************************************************
      IF(OPTIMIZERTYPE.EQ.'CG') RETURN                           !KAESTNERCG
                              CALL TRACE$PUSH('WAVES$PROPAGATE')
!
!     ==========================================================================
!     ==  STOP WAVE FUNCTIONS                                                 ==
!     ==========================================================================
      IF(TSTOP) THEN
        CALL WAVES$STOP()
        TSTOP=.FALSE.
      END IF
!
!     ==========================================================================
!     ==  PROPAGATE WAVE FUNCTIONS (PUT PSI(+) INTO PSIM)                     ==
!     ==========================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      DO IKPT=1,NKPTL
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL=GSET%NGL
!       == EVALUATE FACTORS ====================================================
        ALLOCATE(ARR1(NGL))
        ALLOCATE(ARR2(NGL))
        ALLOCATE(ARR3(NGL))
        IF(TSTRESS) THEN
          CALL CELL$GETR8A('PSIANNEMAT',9,FRICMAT) !\ALPHADOT*DELTA/2
! ==============================================================================
! THERE IS A FRICTION LIKE TERM WHEN A UNIT CELL CHANGES, WHICH CAN BE
! REMOVED WITHOUT CHANGING THE ENERGY CONSERVATION, IF AN ANALOGOUS
! TERM IN THE STRESSES IS REMOVED AS WELL. THIS FRICTION LIKE TERM IS
! SWITCHED OFF BY THE PARAMETER "TNORED" IN THE CELL OBJECT.
! ==============================================================================
          ALLOCATE(GVEC(3,NGL))
          CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
          ALLOCATE(ANNEEVEC(NGL))
!!$          CALL WAVES_DMPSI(.FALSE.,ANNEE,FRICMAT,NGL,GVEC &
!!$     &                      ,EMASS,EMASSCG2,GSET%DBUCKET,ANNEEVEC)
          CALL WAVES_DMPSI(TSTRESS,ANNEE,FRICMAT,NGL,GVEC &
     &                      ,EMASS,EMASSCG2,GSET%DBUCKET,ANNEEVEC)
          DEALLOCATE(GVEC)
          ARR1(:)=2.D0/(1.D0+ANNEEVEC(:))
          ARR2(:)=1.D0-ARR1(:)
          ARR3(:)=-DELT**2/(1.D0+ANNEEVEC(:))/GSET%MPSI(:)
          DEALLOCATE(ANNEEVEC)
!PB070802        IF(ASSOCIATED(GSET%DMPSI)) THEN
!PB070802          DO IG=1,NGL
!PB070802!PB         SVAR1=1.D0+ANNEE+GSET%DMPSI(IG)
!PB070802            SVAR1=1.D0+ANNEE
!PB070802            ARR1(IG)=2.D0/SVAR1
!PB070802            ARR2(IG)=1.D0-ARR1(IG)
!PB070802            ARR3(IG)=-DELT**2/GSET%MPSI(IG)/SVAR1
!PB070802          ENDDO
        ELSE
          SVAR1=2.D0/(1.D0+ANNEE)
          SVAR2=1.D0-SVAR1
          SVAR3=-DELT**2/(1.D0+ANNEE)
          DO IG=1,NGL
            ARR1(IG)=SVAR1
            ARR2(IG)=SVAR2
            ARR3(IG)=SVAR3/GSET%MPSI(IG)
          ENDDO
        END IF
!       ==  NOW PROPAGATE ======================================================
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(.NOT.ASSOCIATED(THIS%HPSI)) THEN
            CALL ERROR$MSG('HPSI NOT AVAILABLE')
            CALL ERROR$MSG('EITHER WAVES$ETOT HAS NOT BEEN CALLED')
            CALL ERROR$MSG('OR WAVES$PROPAGATE HAS BEEN CALLED PREVIOUSLY')
            CALL ERROR$I4VAL('IKPT',IKPT)
            CALL ERROR$I4VAL('ISPIN',ISPIN)
            CALL ERROR$STOP('WAVES$PROPAGATE')
          END IF
          NBH=THIS%NBH
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          TRESIDENTPSIM=CPPAW_CUBLAS_ACC_PSIM_PROPAGATE_ENABLED() &
     &                  .AND.(.NOT.TSTRESS)
          TPHASERESIDENTPSIM= &
     &       CPPAW_CUBLAS_ACC_PSIM_PHASE_RESIDENCY_ENABLED() &
     &       .AND.(.NOT.TSTRESS)
          IF(TRESIDENTPSIM) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
            CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &          ('ACC_PRESENT_PROP_PSI0','ACC_COPY_PROP_PSI0_IN' &
     &          ,NGL,NDIM,NBH,THIS%PSI0)
            IF(TPHASERESIDENTPSIM) THEN
              CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &            ('ACC_PRESENT_PROP_PSIM','ACC_COPY_PROP_PSIM_IN' &
     &            ,NGL,NDIM,NBH,THIS%PSIM)
            ELSE
              CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D_INOUT &
     &            ('ACC_PRESENT_PROP_PSIM','ACC_COPY_PROP_PSIM_IN' &
     &            ,'ACC_COPY_PROP_PSIM_OUT',NGL,NDIM,NBH,THIS%PSIM)
            END IF
            CALL CPPAW_CUBLAS_ACC_PROFILE_PRESENT_C8_3D &
     &          ('ACC_PRESENT_PROP_HPSI','ACC_COPY_PROP_HPSI_IN' &
     &          ,NGL,NDIM,NBH,THIS%HPSI)
            CALL ACCELPROFILE$ADD('ACC_COPY_PROP_COEFF_IN' &
     &          ,INT(NGL,KIND=8),3_8,0_8,0_8 &
     &          ,0.D0,8.D0*3.D0*REAL(NGL,KIND=8),0.D0)
#ENDIF
            IF(TPHASERESIDENTPSIM.AND.(.NOT.THIS%PSIM_ACC_RESIDENT)) THEN
!$ACC ENTER DATA COPYIN(THIS%PSIM(1:NGL,1:NDIM,1:NBH))
              THIS%PSIM_ACC_RESIDENT=.TRUE.
            END IF
            CALL WAVES_PROPAGATE_PSIM_ACC(NGL,NDIM,NBH,ARR1,ARR2,ARR3 &
     &                                    ,THIS%PSI0,THIS%PSIM,THIS%HPSI)
            IF(THIS%HPSI_ACC_RESIDENT) THEN
!$ACC EXIT DATA DELETE(THIS%HPSI(1:NGL,1:NDIM,1:NBH))
            END IF
            THIS%HPSI_ACC_RESIDENT=.FALSE.
          ELSE
            IF(THIS%HPSI_ACC_RESIDENT) THEN
#IF DEFINED(CPPVAR_ACCEL_PROFILE)
              CALL ACCELPROFILE$ADD('ACC_COPY_PROP_HPSI_HOST_OUT' &
     &            ,INT(NGL,KIND=8),INT(NDIM,KIND=8),INT(NBH,KIND=8) &
     &            ,0_8,0.D0,16.D0*REAL(NGL,KIND=8) &
     &            *REAL(NDIM,KIND=8)*REAL(NBH,KIND=8),0.D0)
#ENDIF
!$ACC EXIT DATA COPYOUT(THIS%HPSI(1:NGL,1:NDIM,1:NBH))
              THIS%HPSI_ACC_RESIDENT=.FALSE.
            END IF
#ENDIF
            DO IB=1,NBH
              DO IDIM=1,NDIM
                DO IG=1,NGL
                  THIS%PSIM(IG,IDIM,IB)=ARR1(IG)*THIS%PSI0(IG,IDIM,IB) &
       &                               +ARR2(IG)*THIS%PSIM(IG,IDIM,IB) &
       &                               +ARR3(IG)*THIS%HPSI(IG,IDIM,IB)
                ENDDO
              ENDDO
            ENDDO
#IF DEFINED(CPPVAR_CUBLAS_ACC)
          END IF
#ENDIF
          DEALLOCATE(THIS%HPSI)
        ENDDO
        DEALLOCATE(ARR1)
        DEALLOCATE(ARR2)
        DEALLOCATE(ARR3)
      ENDDO
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_PROPAGATE_PSIM_ACC(NGL,NDIM,NBH,ARR1,ARR2,ARR3 &
     &                                    ,PSI0,PSIM,HPSI)
!     **************************************************************************
!     **  GPU PROTOTYPE FOR THE PROPAGATION UPDATE OF PSIM.                  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NGL
      INTEGER(4),INTENT(IN)    :: NDIM
      INTEGER(4),INTENT(IN)    :: NBH
      REAL(8)   ,INTENT(IN)    :: ARR1(NGL)
      REAL(8)   ,INTENT(IN)    :: ARR2(NGL)
      REAL(8)   ,INTENT(IN)    :: ARR3(NGL)
      COMPLEX(8),INTENT(IN)    :: PSI0(NGL,NDIM,NBH)
      COMPLEX(8),INTENT(INOUT) :: PSIM(NGL,NDIM,NBH)
      COMPLEX(8),INTENT(IN)    :: HPSI(NGL,NDIM,NBH)
      INTEGER(4)               :: IG
      INTEGER(4)               :: IDIM
      INTEGER(4)               :: IB
!     **************************************************************************
!$ACC DATA PRESENT_OR_COPY(PSIM(1:NGL,1:NDIM,1:NBH)) &
!$ACC& PRESENT_OR_COPYIN(PSI0(1:NGL,1:NDIM,1:NBH) &
!$ACC&       ,ARR1(1:NGL),ARR2(1:NGL),ARR3(1:NGL)) &
!$ACC& PRESENT_OR_COPYIN(HPSI(1:NGL,1:NDIM,1:NBH))
!$ACC PARALLEL LOOP COLLAPSE(3) &
!$ACC& PRESENT(PSI0(1:NGL,1:NDIM,1:NBH),PSIM(1:NGL,1:NDIM,1:NBH) &
!$ACC&        ,HPSI(1:NGL,1:NDIM,1:NBH) &
!$ACC&        ,ARR1(1:NGL),ARR2(1:NGL),ARR3(1:NGL))
      DO IB=1,NBH
        DO IDIM=1,NDIM
          DO IG=1,NGL
            PSIM(IG,IDIM,IB)=ARR1(IG)*PSI0(IG,IDIM,IB) &
     &                       +ARR2(IG)*PSIM(IG,IDIM,IB) &
     &                       +ARR3(IG)*HPSI(IG,IDIM,IB)
          ENDDO
        ENDDO
      ENDDO
!$ACC END PARALLEL LOOP
!$ACC END DATA
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_PSIMASS
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE WAVE FUNCTION MASS AND ITS CHANGES                   **
!     **                                                                      **
!     **************************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)   ,ALLOCATABLE :: GVEC(:,:)
!PB070802      REAL(8)   ,ALLOCATABLE :: DG2(:)
!PB070802      REAL(8)                :: FRICMAT(3,3)
      REAL(8)                :: G2
      REAL(8)                :: FAC1
      INTEGER(4)             :: IKPT,IG
      INTEGER(4)             :: NGL
      LOGICAL(4)             :: TSTRESS
!     **************************************************************************
      CALL CELL$GETL4('ON',TSTRESS)
!PB070802      IF(TSTRESS) THEN
!PB070802        CALL CELL$GETR8A('FRICMAT',9,FRICMAT) !\ALPHADOT*DELTA/2
!PB070802      ELSE
!PB070802        FRICMAT=0.D0
!PB070802      END IF
      FAC1=EMASS*2.D0*EMASSCG2
      DO IKPT=1,NKPTL
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL=GSET%NGL
        IF(TSTRESS) THEN
          ALLOCATE(GVEC(3,NGL))
          CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
          IF(.NOT.ASSOCIATED(GSET%MPSI)) ALLOCATE(GSET%MPSI(NGL))
!PB070802    IF(.NOT.ASSOCIATED(GSET%DMPSI))ALLOCATE(GSET%DMPSI(NGL))
!PB070802    ALLOCATE(DG2(NGL))
!PB070802    DO IG=1,NGL
!PB070802      DG2(IG)=FRICMAT(1,1)     *GVEC(1,IG)*GVEC(1,IG) &
!PB070802   &        +FRICMAT(1,2)*2.D0*GVEC(1,IG)*GVEC(2,IG) &
!PB070802   &        +FRICMAT(1,3)*2.D0*GVEC(1,IG)*GVEC(3,IG) &
!PB070802   &        +FRICMAT(2,2)     *GVEC(2,IG)*GVEC(2,IG) &
!PB070802   &        +FRICMAT(2,3)*2.D0*GVEC(2,IG)*GVEC(3,IG) &
!PB070802   &        +FRICMAT(3,3)     *GVEC(3,IG)*GVEC(3,IG)
!PB070802 ENDDO
!         == MASS ==============================================================
          DO IG=1,NGL
            G2=GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2
            GSET%MPSI(IG) =EMASS*(1.D0+EMASSCG2*G2)
!PB070802   GSET%DMPSI(IG)=EMASS*EMASSCG2*2.D0*DG2(IG)
          END DO
!         == ADD MASS FOR BUCKET POTENTIAL =====================================
          IF(TBUCKET) THEN
            DO IG=1,NGL
              GSET%MPSI(IG) =GSET%MPSI(IG) +FAC1*GSET%BUCKET(IG) 
!PB070802     GSET%DMPSI(IG)=GSET%DMPSI(IG)+FAC1*GSET%DBUCKET(IG)*DG2(IG)
            ENDDO
          END IF
!PB070802 DEALLOCATE(DG2)
          DEALLOCATE(GVEC)
        ELSE 
!PB070802 IF(ASSOCIATED(GSET%DMPSI))DEALLOCATE(GSET%DMPSI)   
          IF(.NOT.ASSOCIATED(GSET%MPSI)) THEN
            ALLOCATE(GVEC(3,NGL))
            CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
            ALLOCATE(GSET%MPSI(NGL))
            DO IG=1,NGL
              G2=GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2
              GSET%MPSI(IG)=EMASS*(1.D0+EMASSCG2*G2)
            ENDDO
            DEALLOCATE(GVEC)
            IF(TBUCKET) THEN
              DO IG=1,NGL
                GSET%MPSI(IG)=GSET%MPSI(IG)+FAC1*GSET%BUCKET(IG)
              ENDDO
            END IF
          END IF
        END IF
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_DMPSI(TSTRESS,ANNEE,FRICMAT,NG,GVEC &
     &                      ,EMASS,EMASSCG2,DBUCKET,ANNEEVEC)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)  :: TSTRESS             
      REAL(8)   ,INTENT(IN)  :: ANNEE
      REAL(8)   ,INTENT(IN)  :: FRICMAT(3,3)
      INTEGER(4),INTENT(IN)  :: NG
      REAL(8)   ,INTENT(IN)  :: GVEC(3,NG)
      REAL(8)   ,INTENT(IN)  :: EMASS
      REAL(8)   ,INTENT(IN)  :: EMASSCG2
      REAL(8)   ,INTENT(IN)  :: DBUCKET(NG)
      REAL(8)   ,INTENT(OUT) :: ANNEEVEC(NG)
      REAL(8)                :: FAC1,FAC2,SVAR
      INTEGER(4)             :: IG
!     **************************************************************************
      IF(.NOT.TSTRESS) THEN
        ANNEEVEC(:)=ANNEE
        RETURN
      END IF
      FAC1=ANNEE-0.5D0*(FRICMAT(1,1)+FRICMAT(2,2)+FRICMAT(3,3))
      FAC2=-0.5D0*EMASS*2.D0*EMASSCG2
      DO IG=1,NG
        SVAR=FRICMAT(1,1)     *GVEC(1,IG)*GVEC(1,IG) &
     &      +FRICMAT(1,2)*2.D0*GVEC(1,IG)*GVEC(2,IG) &
     &      +FRICMAT(1,3)*2.D0*GVEC(1,IG)*GVEC(3,IG) &
     &      +FRICMAT(2,2)     *GVEC(2,IG)*GVEC(2,IG) &
     &      +FRICMAT(2,3)*2.D0*GVEC(2,IG)*GVEC(3,IG) &
     &      +FRICMAT(3,3)     *GVEC(3,IG)*GVEC(3,IG)
        ANNEEVEC(IG)=FAC1+FAC2*SVAR*(1.D0+DBUCKET(IG))
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_WAVEKINETIC(EKIN)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES FICTITIOUS KINETIC ENERGY OF THE WAVE FUNCTIONS          **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **    PARAMETERS DELT,EMASS,EMASSCG2 MUST BE SET                        **
!     **                                                                      **
!     **************************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT):: EKIN
      INTEGER(4)            :: IKPT,ISPIN,IG,IDIM,IB
      INTEGER(4)            :: NGL
      INTEGER(4)            :: NBX
      COMPLEX(8),ALLOCATABLE:: TPSIP(:)
      COMPLEX(8),ALLOCATABLE:: TPSIM(:)
      REAL(8)               :: EKIN1,SUM
      REAL(8)   ,ALLOCATABLE:: OCC(:,:,:)
      REAL(8)   ,ALLOCATABLE:: EIG(:,:,:)
      REAL(8)               :: RBAS(3,3),GBAS(3,3),CELLVOL
      REAL(8)               :: V1,V2
      REAL(8)               :: F1,F2
      LOGICAL(4)            :: TINV
      COMPLEX(8)            :: CSVAR
      REAL(8)   ,PARAMETER  :: DSMALL=1.D-12
      INTEGER(4)            :: IB1,IB2
!     **************************************************************************
!
!     ==========================================================================
!     == COLLECT UNIT CELL VOLUME                                             ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL GBASS(RBAS,GBAS,CELLVOL)
!
!     ==========================================================================
!     == COLLECT OCCUPATIONS                                                  ==
!     ==========================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      ALLOCATE(EIG(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
      EIG(:,:,:)=0.D0
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      EKIN=0.D0
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          TINV=GSET%TINV
          ALLOCATE(TPSIP(NGL))
          ALLOCATE(TPSIM(NGL))
          EKIN1=0.D0
          DO IB=1,THIS%NBH
            IF(TINV) THEN 
              IB1=2*IB-1
              IB2=2*IB
              F1=OCC(IB1,IKPT,ISPIN)
              F2=OCC(IB2,IKPT,ISPIN)
            ELSE
              F1=OCC(IB,IKPT,ISPIN)
              F2=0.D0
            END IF
            IF(F1.EQ.0.D0.AND.F2.EQ.0.D0) CYCLE
!
!           ====================================================================
!           ==  WAVE FUNCTION KINETIC ENERGY FOR GENERAL WAVE FUNCTIONS       ==
!           ==  AND FIRST PART FOR SUPER WAVE FUNCTIONS                       ==
!           ====================================================================
            SUM=0.D0
            DO IDIM=1,NDIM
              DO IG=1,NGL
                CSVAR=THIS%PSI0(IG,IDIM,IB)-THIS%PSIM(IG,IDIM,IB)
                V1=REAL(CSVAR,KIND=8)
                V2=AIMAG(CSVAR)
                SUM=SUM+GSET%MPSI(IG)*(V1**2+V2**2)
              ENDDO
            ENDDO
            SUM=SUM/DELT**2
!
            IF(.NOT.TINV) THEN
              EKIN1=EKIN1+SUM*F1
              EIG(IB,IKPT,ISPIN)=SUM
              CYCLE
            ELSE 
              EIG(IB1,IKPT,ISPIN)=0.5D0*SUM              
              EIG(IB2,IKPT,ISPIN)=0.5D0*SUM              
              EKIN1=EKIN1+SUM*0.5D0*(F1+F2)
            END IF
!
!           ====================================================================
!           ==  NOW CONTINUE WITH SUPERWAVE FUNCTIONS ONLY                    ==
!           ====================================================================
!           == PSI- ONLY REQUIRED IF THE OCCUPATIONS DIFFER ==========
            IF(F1-F2.EQ.0.D0) CYCLE
            SUM=0.D0
            DO IDIM=1,NDIM
              TPSIP(:)=THIS%PSI0(:,IDIM,IB)-THIS%PSIM(:,IDIM,IB)
              CALL PLANEWAVE$INVERTG(NGL,TPSIP,TPSIM)
              DO IG=1,NGL
                SUM=SUM+GSET%MPSI(IG) &
     &                 *(REAL(TPSIP(IG),KIND=8)*REAL(TPSIM(IG),KIND=8) &
     &                   +AIMAG(TPSIP(IG))*AIMAG(TPSIM(IG)))
              ENDDO            
            ENDDO
            SUM=SUM/DELT**2
            EIG(IB1,IKPT,ISPIN)=EIG(IB1,IKPT,ISPIN)+0.5D0*SUM              
            EIG(IB2,IKPT,ISPIN)=EIG(IB2,IKPT,ISPIN)-0.5D0*SUM              
            EKIN1=EKIN1+SUM*0.5D0*(F1-F2)
          ENDDO
          DEALLOCATE(TPSIP)
          DEALLOCATE(TPSIM)
          EKIN=EKIN+EKIN1*CELLVOL
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == COMMUNICATE                                                          ==
!     ==========================================================================
      CALL MPE$COMBINE('MONOMER','+',EKIN)
      CALL MPE$COMBINE('K','+',EIG)
!
!     ==========================================================================
!     == SET CONTRIBUTION TO EIGENVALUES                                      ==
!     ==========================================================================
      CALL WAVES_DYNOCCSETR8A('M<PSIDOT|PSIDOT>',NBX*NKPTL*NSPIN,EIG)
      DEALLOCATE(EIG)
      DEALLOCATE(OCC)
      RETURN
      END
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WAVES_KDISTRIBUTE(NTASKS,NKPT,TINV,ICOLOR,KMAP)
!     **************************************************************************
!     **   ICOLOR GIVES EACH TASK AN INTEGER. TASKS WITH EQUAL                **
!     **   VALUE OF ICOLOR WILL BE IN THE SAME GROUP                          **
!     **                                                                      **
!     **   KMAP POINTS TO THE TASK, WHICH IS THE FIRST TASK IN THE            **
!     **   CORRESPONDING K-GROUP                                              **
!     **                                                                      **
!     ******************************* PETER BLOECHL, TU CLAUSTJAL 2005 *********
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NTASKS      ! #(TASKS)
       INTEGER(4),INTENT(IN) :: NKPT        ! #(K-POINTS)
       LOGICAL(4),INTENT(IN) :: TINV(NKPT)  ! TRUE FOR SPECIAL K-POINTS
       INTEGER(4),INTENT(OUT):: ICOLOR(NTASKS) !USED TO SPLIT PROCESSORS
       INTEGER(4),INTENT(OUT):: KMAP(NKPT)   ! FIRST TASK IN K-GROUP
       REAL(8)   ,PARAMETER  :: WA=1.D0,WB=0.5D0,WA0=0.2,WB0=0.2
       REAL(8)               :: WK(NKPT)
       REAL(8)               :: WK0(NKPT)
       INTEGER(4)            :: SK(NKPT)
       INTEGER(4)            :: SKSAVE(NKPT)
       REAL(8)   ,ALLOCATABLE:: XLOAD(:)
       REAL(8)               :: XLOAD1
       REAL(8)   ,PARAMETER  :: R8LARGE=1.D+20
       INTEGER(4)            :: I,J,IKPT,IPOS,IGROUP
       INTEGER(4)            :: ISVAR
       INTEGER(4)            :: IND(1)
       INTEGER(4)            :: IFIRST(NTASKS)
       INTEGER(4),ALLOCATABLE:: NA(:),NB(:)
       INTEGER(4),ALLOCATABLE:: MG(:),MGSAVE(:)
       INTEGER(4)            :: NGROUPS
!      *************************************************************************
       DO I=1,NKPT
         WK(I)=WA
         WK0(I)=WA0
         IF(TINV(I))WK(I)=WB
         IF(TINV(I))WK0(I)=WB0
       ENDDO
!
!      =========================================================================
!      ==  OPTMIZATION OF SK. SK US A POINTER IKPT->IGROUP                    ==
!      =========================================================================
       XLOAD1=R8LARGE
!      == TEST DIFFERENT NUMBER OF SUBGROUPS ===================================
       DO NGROUPS=1,NTASKS       
         IF(NGROUPS.GT.NKPT) EXIT  ! THERE MUST NOT BE AN EMPTY GROUP
!        =======================================================================
!        == DISTRIBUTE ALL K-POINTS INTO "NGROUPS" GROUPS                     ==
!        == NB(IGROUP) IS THE NUMBER OF KPOINTS WITH REAL WAVE                ==
!        == FUNCTIONS IN THE GROUP "IGROUP"; NA(IGROUP IS THE                 ==
!        == NUMBER OF COMPLEX WAVE FUNCTIONS                                  ==
!        =======================================================================
         ALLOCATE(NA(NGROUPS)) !#(K-POINTS WITH COMPLEX PSI)
         ALLOCATE(NB(NGROUPS)) !#(K-POINTS WITH REAL PSI)
         NA(:)=0
         NB(:)=0
         DO IKPT=1,NKPT
           IPOS=MODULO(IKPT-1,NGROUPS)+1
           IF(TINV(IKPT)) THEN
             NB(IPOS)=NB(IPOS)+1
           ELSE
             NA(IPOS)=NA(IPOS)+1
           END IF
         ENDDO
!
!        === RESHUFFLE NA AND NB TO DISTRIBUTE LOAD EQUALLY AMONG GROUPS =======
         CALL WAVES_PREBALANCE(NGROUPS,NA,NB,WA,WB)
!
!        == DETERMINE POINTER SK(IKPT)=IGROUP
         DO IKPT=1,NKPT
           DO IGROUP=1,NGROUPS
             IF(TINV(IKPT).AND.NB(IGROUP).NE.0) THEN
               NB(IGROUP)=NB(IGROUP)-1   ! FOR CONSISTENCY CHECK
               SK(IKPT)=IGROUP
               EXIT
             ELSE IF(.NOT.TINV(IKPT).AND.NA(IGROUP).NE.0) THEN
               NA(IGROUP)=NA(IGROUP)-1   ! FOR CONSISTENCY CHECK
               SK(IKPT)=IGROUP
               EXIT
             END IF
           ENDDO
         ENDDO
!
!        == CONSISTENCY CHECK ==================================================
         DO IGROUP=1,NGROUPS
           IF(NB(IGROUP).NE.0.OR.NA(IGROUP).NE.0) THEN
             CALL ERROR$MSG('INCONSISTENCY DETECTED WITH NA,NB')
             CALL ERROR$I4VAL('NGROUPS',NGROUPS)
             CALL ERROR$I4VAL('NKPT',NKPT)
             CALL ERROR$I4VAL('NA: #(TINV)',NA(IGROUP))
             CALL ERROR$I4VAL('NA: #(.NOT.TINV)',NB(IGROUP))
             CALL ERROR$STOP('WAVES_KDISTRIBUTE')
           END IF
         ENDDO
         DEALLOCATE(NA)
         DEALLOCATE(NB)
!
!        =======================================================================
!        == CHECK IF THIS #(GROUPS) IMPROVED THE LOADBALANCE                  ==
!        == AND DISCARD IF NOT.                                               ==
!        =======================================================================
         ALLOCATE(MG(NGROUPS))
         ALLOCATE(XLOAD(NGROUPS))
         CALL WAVES_LOADPERPROC(NGROUPS,NTASKS,NKPT,WK,WK0,SK,MG,XLOAD)
         IND=MAXLOC(XLOAD)
         IF(XLOAD(IND(1)).LE.XLOAD1) THEN
           IF(ALLOCATED(MGSAVE))DEALLOCATE(MGSAVE)
           ALLOCATE(MGSAVE(NGROUPS))
           MGSAVE=MG
           SKSAVE=SK
           XLOAD1=XLOAD(IND(1))
         END IF
         DEALLOCATE(MG)
         DEALLOCATE(XLOAD)
       ENDDO
       NGROUPS=SIZE(MGSAVE)
       ALLOCATE(MG(NGROUPS))
       MG(:)=MGSAVE(:)  ! #(PROCS PER GROUP)
       SK=SKSAVE        ! POINTER: IKPT-> IGROUP
       DEALLOCATE(MGSAVE)
!
!      =========================================================================
!      ==  WRAP UP  AND PREPARE OUTPUT                                        ==
!      =========================================================================
       ISVAR=0
       DO I=1,NGROUPS
         IFIRST(I)=ISVAR+1   ! FIRST TASK FOR THIS GROUP
         ISVAR=ISVAR+MG(I)
         ICOLOR(IFIRST(I):ISVAR)=I   
       ENDDO
       DO IKPT=1,NKPT 
         KMAP(IKPT)=IFIRST(SK(IKPT))
       ENDDO
!
!      =========================================================================
!      == CHECKS                                                              ==
!      =========================================================================
       IF(SUM(MG).NE.NTASKS) THEN
         CALL ERROR$MSG('INCONSISTENT DIVISION')
         CALL ERROR$STOP('WAVES_KDISTRIBUTE')
       END IF
       ALLOCATE(NA(NTASKS))
       NA(:)=0
       DO IKPT=1,NKPT
         IGROUP=SK(IKPT)
         I=IFIRST(IGROUP)
         J=IFIRST(IGROUP)-1+MG(IGROUP)
         NA(I:J)=NA(I:J)+1
       ENDDO
       DO I=1,NTASKS
         IF(NA(I).EQ.0) THEN
           CALL ERROR$MSG('TASK WITH NKPTL=0 IS NOT ALLOWED')
           CALL ERROR$MSG('REDUCE NUMBER OF PROCESSORS')
           CALL ERROR$STOP('WAVES_KDISTRIBUTE')
         END IF
       ENDDO              

       DO I=1,NTASKS
         IF(ICOLOR(I).EQ.0) THEN
           CALL ERROR$MSG('ICOLOR HAS AN ELEMENT ZERO')
           CALL ERROR$I4VAL('I',I)
           CALL ERROR$I4VAL('ICOLOR(I)',ICOLOR(I))
           CALL ERROR$STOP('WAVES_KDISTRIBUTE')
         END IF
       ENDDO
       DO I=1,NKPT
         IF(KMAP(I).EQ.0) THEN
           CALL ERROR$MSG('KMAP HAS AN ELEMENT ZERO')
           CALL ERROR$I4VAL('I',I)
           CALL ERROR$I4VAL('KMAP(I)',KMAP(I))
           CALL ERROR$STOP('WAVES_KDISTRIBUTE')
         END IF
       ENDDO

       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE WAVES_LOADPERPROC(NGROUPS,NTASKS,NKPT,WK,WK0,SK,MG,XLOAD)
!      *************************************************************************
!      **  ASSIGN TO EACH GROUP A NUMBER OF PROCESSORS                        **
!      **                                                                     **
!      ********************* PETER BLOECHL, TU CLAUSTJAL 2005 ******************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NGROUPS
       INTEGER(4),INTENT(IN) :: NTASKS
       INTEGER(4),INTENT(IN) :: NKPT
       REAL(8)   ,INTENT(IN) :: WK(NKPT)     ! WEIGHT PER K-POINT
       REAL(8)   ,INTENT(IN) :: WK0(NKPT)     ! WEIGHT PER K-POINT
       INTEGER(4),INTENT(IN) :: SK(NKPT)     ! GROUP INDEX
       INTEGER(4),INTENT(OUT):: MG(NGROUPS)   ! #PROC FOR K-GROUP
       REAL(8)   ,INTENT(OUT):: XLOAD(NGROUPS) ! LOAD PER PROCFOR EACH GROUP
       REAL(8)               :: WG(NGROUPS) ! TOTAL LOAD OF THE K-GROUP
       REAL(8)               :: WG0(NGROUPS) ! TOTAL LOAD OF THE K-GROUP
       INTEGER(4)            :: NODESLEFT
       INTEGER(4)            :: IPOS1,IPOS2
       REAL(8)               :: XLOADTEST(NGROUPS) ! LOAD PER PROCESSOR
       REAL(8)               :: XLOAD1
       INTEGER(4)            :: IND(1)
       INTEGER(4)            :: I
       REAL(8)   ,PARAMETER  :: R8LARGE=1.D+20
!      *************************************************************************
!
!      =========================================================================
!      == START UP                                                            ==
!      =========================================================================
!      == DETERMINE WORK FOR EACH GROUP
       WG(:)=0.D0
       WG0(:)=0.D0
       DO I=1,NKPT
         IF(SK(I).LT.1.OR.SK(I).GT.NGROUPS) THEN
           CALL ERROR$MSG('SK OUT OF RANGE')
           CALL ERROR$STOP('WAVES_LOADPERPROC')
         END IF
         WG(SK(I))=WG(SK(I))+WK(I)
         WG0(SK(I))=WG0(SK(I))+WK0(I)
       ENDDO
       IF(NGROUPS.GT.NTASKS) THEN
         CALL ERROR$MSG('NGROUPS MUST BE SMALLER OR EQUAL TO NTASKS')
         CALL ERROR$STOP('WAVES_LOADPERPROC')
         XLOAD(:)=0
         XLOAD(1:NTASKS)=1.D0
         MG(:)=0
         MG(1:NTASKS)=1
         RETURN
       END IF
       MG(:)=1  ! GIVE EACH GROUP ONE PROCESSOR
! 
!      =========================================================================
!      == ENFORCE SUM RULE                                                    ==
!      =========================================================================
       NODESLEFT=NTASKS-SUM(MG)  ! #(UNOCCUPIED NODES)
       IF(NODESLEFT.LT.0) THEN
         CALL ERROR$MSG('NODESLEFT<0')
         CALL ERROR$STOP('WAVES_LOADPERPROC')
       END IF
       DO I=1,NODESLEFT
         XLOAD(:)=WG0(:)+WG(:)/MG(:)  !DETERMINE LOAD FOR EACH GROUP
         IND=MAXLOC(XLOAD)
         IPOS1=IND(1)
         MG(IPOS1)=MG(IPOS1)+1  ! GIVE MAXIMUM LOADED GROUP AN ADDITIONAL NODE
       ENDDO  
!
!      =========================================================================
!      == RESHUFFLE NODES TO REDUCE MAXIMUM LOAD                              ==
!      =========================================================================
1000   CONTINUE
       XLOAD(:)=WG0(:)+WG(:)/MG(:)  ! CURRENT WORK PER PROC IN EACH GROUP
       IND=MAXLOC(XLOAD)
       IPOS1=IND(1)
       XLOAD1=XLOAD(IPOS1)          ! MAX WORK PER PROC
       XLOADTEST(IPOS1)=WG0(IPOS1)+WG(IPOS1)/REAL(MG(IPOS1)+1)
       DO I=1,NGROUPS
         IF(I.EQ.IPOS1) CYCLE
         IF(MG(I).GT.1) THEN
           XLOADTEST(I)=WG0(I)+WG(I)/REAL(MG(I)-1)
         ELSE
           XLOADTEST(I)=R8LARGE
         END IF
       ENDDO
       IND=MINLOC(XLOADTEST)
       IPOS2=IND(1)
       IF(IPOS2.NE.IPOS1.AND.XLOADTEST(IPOS2).LT.XLOAD1) THEN
         MG(IPOS1)=MG(IPOS1)+1
         MG(IPOS2)=MG(IPOS2)-1
         GOTO 1000
       END IF
!
!      =========================================================================
!      == DETERMINE LOAD FOR EACH PROCESSOR                                   ==
!      =========================================================================
       XLOAD(:)=WG0(:)+WG(:)/MG(:)
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE WAVES_PREBALANCE(NGROUPS,NA,NB,WA,WB)
!      *************************************************************************
!      **  TRIES DO DISTRIBUTE THE COMPUTATIONAL EFFORT EQUALLY               **
!      **  ONTO NGROUPS GROUPS. IT ASSUMES THAT THERE ARE TWO TYPES           **
!      **  OF KPOINTS WITH WEIGHTS WA AND WB. NA AND NB SPECIFIES             **
!      **  THE NUMBER OF K-POINTS FROM EACH TYPE IN THE CORRESPONDING         **
!      **  GROUPS                                                             **
!      **                                                                     **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)    :: NGROUPS
       INTEGER(4),INTENT(INOUT) :: NA(NGROUPS)
       INTEGER(4),INTENT(INOUT) :: NB(NGROUPS)
       REAL(8)   ,INTENT(IN)    :: WA
       REAL(8)   ,INTENT(IN)    :: WB
       REAL(8)                  :: XLOAD(NGROUPS)
       REAL(8)                  :: XLOADTESTA(NGROUPS)
       REAL(8)                  :: XLOADTESTB(NGROUPS)
       INTEGER(4)               :: IND(1)
       REAL(8)                  :: XLOAD1,XLOAD2A,XLOAD2B
       INTEGER(4)               :: IPOS1,IPOS2A,IPOS2B
       INTEGER(4)               :: J,ICOUNT
       REAL(8)   ,PARAMETER     :: R8LARGE=1.D20
!      ***************************************************************
       ICOUNT=0
       DO ICOUNT=1,100
!        == XLOAD QUANTIFIES THE ESTIMATED COMPUTATIONAL EFFORT FOR EACH GROUP
         XLOAD(:)=REAL(NA(:),KIND=8)*WA+REAL(NB(:),KIND=8)*WB
!        == FIND THE GROUP HAVING THE LARGEST LOAD
         IND=MAXLOC(XLOAD)
         IPOS1=IND(1)
         XLOAD1=XLOAD(IPOS1)
!        == IT IS NOT ALLOWED TO REDUCE THE NUMBER OF KPOINTS IN A GROUP TO ZERO
         IF(NA(IPOS1)+NB(IPOS1).LE.1) EXIT
!        ==  TRY TO RESHUFFLE FROM THE GROUP WITH THE LARGEST LOAD TO THE OTHERS
!        ==  XLOADTESTA ESTIMATES THE LOAD WHEN A K-POINT OF TYPE A IS ADDED
!        ==  XLOADTESTB ESTIMATES THE LOAD WHEN A K-POINT OF TYPE B IS ADDED
!        ==  THE FACTOR IN THE FOLLOWING LINE IS BECAUSE THE PGI COMPILER DOES
!        ==  NOT HANDLE MINLOC WITH THE LARGEST POSSIBLE VALUE.
         XLOADTESTA(IPOS1)=R8LARGE
         XLOADTESTB(IPOS1)=R8LARGE
         DO J=1,NGROUPS
           IF(J.EQ.IPOS1) CYCLE
           XLOADTESTA(J)=XLOAD(J)+WA
           XLOADTESTB(J)=XLOAD(J)+WB
         ENDDO
         IF(NA(IPOS1).EQ.0) XLOADTESTA=R8LARGE
         IF(NB(IPOS1).EQ.0) XLOADTESTB=R8LARGE
         IND=MINLOC(XLOADTESTA)
         IPOS2A=IND(1)       
         XLOAD2A=XLOADTESTA(IPOS2A)
         IND=MINLOC(XLOADTESTB)
         IPOS2B=IND(1)
         XLOAD2B=XLOADTESTB(IPOS2B)
         IF(MIN(XLOAD2A,XLOAD2B).LT.XLOAD1) THEN 
           IF(XLOAD2A.LE.XLOAD1.AND.NA(IPOS1).GT.0) THEN
             NA(IPOS1)=NA(IPOS1)-1        
             NA(IPOS2A)=NA(IPOS2A)+1        
             CYCLE
           ELSE IF(XLOAD2A.GE.XLOAD1.AND.NB(IPOS1).GT.0) THEN
             NB(IPOS1)=NB(IPOS1)-1        
             NB(IPOS2A)=NB(IPOS2A)+1        
             CYCLE
           END IF
         END IF
         EXIT !EXIT LOOK IF THERE IS NO MORE CHANGE
       ENDDO
       RETURN
       END
