!**********************************************************************
!***********************************************************************
!**                                                                   **
!**  WAVES OBJECT                                                     **
!**  This object is continued in paw_waves2.f90                       **
!**                                                                   **
!**  PURPOSE: OPERATIONS ON THE PLANE WAVE PART OF THE WAVE FUNCTIONS **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    WAVES$SETR8A(ID,LEN,VAL)                                       **
!**    WAVES$GETR8A(ID,LEN,VAL)                                       **
!**    WAVES$INITIALIZE                                               **
!**    WAVES$ETOT                                                     **
!**    WAVES$PROPAGATE                                                **
!**    WAVES$ORTHOGONALIZE                                            **
!**    WAVES$SWITCH                                                   **
!**    WAVES$REPORT(NFIL)                                             **
!**                                                                   **
!**  DEPENDECIES:                                                     **
!**    ERROR                                                          **
!**    MPELIB                                                         **
!**    PLANEWAVE                                                      **
!**    SETUP                                                          **
!**    ...                                                            **
!**                                                                   **
!**  definition of spin dimensions:                                   **
!**                           | nspin | ndim | ndimd |                **
!**        -------------------------------------------                **
!**        non-spin-polarized |   1   |   1  |   1   |                **
!**        spin polarized     |   2   |   1  |   2   |                **
!**        non-collinear      |   1   |   2  |   4   |                **
!**                                                                   **
!**     nspin is the number of slater determinants                    **
!**     ndim is the number of spinor dimensions of the wave function  **
!**     ndimd is the number of density components. Two representations**
!**           are used:                                               **
!**           nspin=2,ndim=1:  up/down and total/spin                 **
!**           nspin=1,ndim=2:  upup/updown/downup/downdown and        **
!**                            total/mx/my/mz                         **
!**                                                                   **
!***********************************************************************
!***********************************************************************
!
!......................................................WAVES............
MODULE WAVES_MODULE                                                
!***********************************************************************
!**                                                                   **
!**                                                                   **
!******************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1995)**
USE LINKEDLIST_MODULE
! ----------------------------------------------------------------------
! THE TYPE MAP_TYPE DESCRIBES THE ARRANGEMENT OF PROJECTOR FUNCTIONS 
! AND THEIR LINKS TO ATOMS AN ANGULAR MOMENTA
TYPE MAP_TYPE
  INTEGER(4)         :: NAT         ! #(ATOms)
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
! ----------------------------------------------------------------------
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
  REAL(8)   ,POINTER :: DMPSI(:)   ! (NGL) DMPSI/DT
  REAL(8)   ,POINTER :: BUCKET(:)  ! BUCKET POTENTIAL
  REAL(8)   ,POINTER :: DBUCKET(:) ! 1/G* DBUCKET/DG
END TYPE GSET_TYPE
TYPE WVSET_TYPE  !======================================================
  TYPE(GSET_TYPE),POINTER :: GSET
  INTEGER(4)         :: NB
  INTEGER(4)         :: NBH
  COMPLEX(8),POINTER :: PSI0(:,:,:)     !(NGL,NDIM,NBH)  PSPSI(0)
  COMPLEX(8),POINTER :: PSIM(:,:,:)     !(NGL,NDIM,NBH)  PSPSI(-,+)(G)
  COMPLEX(8),POINTER :: PROJ(:,:,:)     !(NDIM,NBH,NPRO) <PSPSI|P>
  COMPLEX(8),POINTER :: HPSI(:,:,:)     !(NGWLX,NB,IDIM)
                                        ! +(WAVES$HPSI)-(WAVES$PROPAGATE)
  COMPLEX(8),POINTER :: OPSI(:,:,:)     !(NGWLX,NB,IDIM)
                                        ! +(WAVES$HPSI)-(WAVES$PROPAGATE)
  COMPLEX(8),POINTER :: RLAM0(:,:)      !(NB,NB)
  COMPLEX(8),POINTER :: RLAMM(:,:)      !(NB,NB)
  COMPLEX(8),POINTER :: RLAM2M(:,:)     !(NB,NB)
  COMPLEX(8),POINTER :: RLAM3M(:,:)     !(NB,NB)
  COMPLEX(8),POINTER :: EIGVEC(:,:)     !(NB,NB) ENERGY EIGENVALUES AFTER DIAG
  REAL(8)   ,POINTER :: EIGVAL(:)       !(NB)    EIGENValues AFTER DIAG
  REAL(8)   ,POINTER :: EXPECTVAL(:)    !(NB) !<PSI_N|H|PSI_N>
END TYPE WVSET_TYPE
TYPE EXTERNALPOINTER_TYPE !=============================================
  INTEGER(4) :: IB
  INTEGER(4) :: IKPT
  INTEGER(4) :: ISPIN
  INTEGER(4) :: IAT
  LOGICAL(4) :: TIM
END TYPE EXTERNALPOINTER_TYPE
!========================================================================
!== DATA THAT DESCRIBE FUNCTIONALITY OF WAVES OBJECT                   ==
!==  THESE DATA HAVE TO BE SET EXPLICITELY                             ==
!========================================================================
INTEGER(4)  :: NKPT          ! #(K-POINTS)
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
REAL(8)     :: DELT=0.D0          ! TIME STEP IN A.U.
LOGICAL(4)  :: THAMILTON=.FALSE.  ! HAMILTON MATRIX AVAILABLE
LOGICAL(4)  :: TRAWSTATES=.FALSE. ! PROVIDES NON-DIAGONALIZED WAVE FUNCTIONS THROUGH $GETR
LOGICAL(4)  :: TFORCEX=.TRUE.
LOGICAL(4)  :: TSTRESSX=.FALSE.
!========================================================================
!== PERMANENT DATA, WHICH ARE ORGANIZED BY THE ROUTINES ITSELF         == 
!========================================================================
REAL(8)                   :: WAVEEKIN1=0.D0
REAL(8)                   :: WAVEEKIN2=0.D0
TYPE(WVSET_TYPE),POINTER  :: THISARRAY(:,:)   ! (NKPT,NSPIN)
TYPE(WVSET_TYPE),POINTER  :: THIS            ! CURRENT SET OF WAVES
TYPE(GSET_TYPE) ,POINTER  :: GSET            ! CURRENT SET OF GSET
TYPE(MAP_TYPE)            :: MAP
LOGICAL(4)                :: TPR=.FALSE.
LOGICAL(4)                :: TFIRST=.TRUE.
TYPE(EXTERNALPOINTER_TYPE):: EXTPNTR
LOGICAL(4)                :: TFIXRHO=.FALSE.
LOGICAL(4)                :: TWRITERHO=.FALSE.
CONTAINS
!***********************************************************************
      SUBROUTINE WAVES_SELECTWV(IKPT,ISPIN)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT
      INTEGER(4),INTENT(IN) :: ISPIN
!     ******************************************************************
      IF(.NOT.ASSOCIATED(THISARRAY)) THEN
        CALL ERROR$MSG('THISARRAY DOES NOT EXIST') 
        CALL ERROR$STOP('WAVES_SELECTWV')
      END IF 
      THIS=>THISARRAY(IKPT,ISPIN)
      GSET=>THIS%GSET
      RETURN
      END SUBROUTINE WAVES_SELECTWV
END MODULE WAVES_MODULE
!
!     ..................................................................
      SUBROUTINE WAVES$SETR8(ID,VAL)
!     ******************************************************************
!     **  WAVES$GVECTORS                                              **
!     **  GENERATE G-VECTORS AND OTHER INITIALIZATION                 **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     ******************************************************************
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
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$GETR8(ID,VAL)
!     ******************************************************************
!     **  WAVES$GVECTORS                                              **
!     **  GENERATE G-VECTORS AND OTHER INITIALIZATION                 **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
!     ******************************************************************
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
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$GETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$GETR8A(ID,LEN,VAL)
!     ******************************************************************
!     **  WAVES$GVECTORS                                              **
!     **  GENERATE G-VECTORS AND OTHER INITIALIZATION                 **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
      REAL(8)     ,ALLOCATABLE:: EIGVAL(:)
      COMPLEX(8)  ,ALLOCATABLE:: CWORK1(:)
      COMPLEX(8)  ,ALLOCATABLE:: CWORK2(:)
      COMPLEX(8)  ,ALLOCATABLE:: CWORK3(:)
      INTEGER(4)              :: IKPT,ISPIN,IB,IDIM,IBH,IB1,IB2,IG,IR
      INTEGER(4)              :: IPRO,IPRO1,IPRO2,IAT,ISP,L,LN,LMNX,LMN
      INTEGER(4)              :: NB,NBH,NRL,NPRO,NAT,NGL
      COMPLEX(8)              :: CSVAR1,CSVAR2
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL(4)              :: TINV
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES$GETR8A')
      IF(ID.EQ.'XXXXXXXXXXXXXX') THEN
!     
!     ================================================================
!     ==  REAL SPACE PS WAVE FUNCTION                               ==
!     ==  NOTE: IB,IKPT,ISPIN MUST BE SET ON LINKEDLIST             ==
!     ================================================================
!     ELSE IF(ID.EQ.'<PSI|H|PSI>') THEN
!
!     ================================================================
!     ==  ENERGY EIGENVALUES                                        ==
!     ==  NOTE: IKPT,ISPIN MUST BE SET ON LINKEDLIST                ==
!     ================================================================
      ELSE IF(ID.EQ.'EIGVAL') THEN
        IKPT=EXTPNTR%IKPT
        ISPIN=EXTPNTR%ISPIN
        CALL WAVES_SELECTWV(IKPT,ISPIN)
        CALL PLANEWAVE$SELECT(GSET%ID)
        IF(.NOT.ASSOCIATED(THIS%EIGVAL)) THEN
          CALL ERROR$MSG('EIGENVALUES NOT AVAILABLE')
          CALL ERROR$L4VAL('THAMILTON',THAMILTON)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        IF(LEN.NE.THIS%NB) THEN
          CALL ERROR$MSG('INCONSISTENT LENGTH')
          CALL ERROR$I4VAL('NB',THIS%NB)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
        VAL(:)=THIS%Eigval(:)
!
!     ================================================================
!     ==  energy expectations values of the one-partcle states      ==
!     ==  NOTE: IKPT,ISPIN MUST BE SET ON LINKEDLIST                ==
!     ================================================================
      ELSE IF(ID.EQ.'<PSI|H|PSI>') THEN
        IKPT=EXTPNTR%IKPT
        ISPIN=EXTPNTR%ISPIN
        CALL WAVES_SELECTWV(IKPT,ISPIN)
        CALL PLANEWAVE$SELECT(GSET%ID)
        IF(.NOT.ASSOCIATED(THIS%EXPECTVAL)) THEN
          CALL ERROR$MSG('ENERGY EXPECTATION VALUES NOT AVAILABLE')
          CALL ERROR$L4VAL('THAMILTON',THAMILTON)
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
!     ================================================================
!     ==  REAL SPACE PS WAVE FUNCTION                               ==
!     ==  NOTE: IB,IKPT,ISPIN MUST BE SET ON LINKEDLIST             ==
!     ================================================================
      ELSE IF(ID.EQ.'PSPSI') THEN
        IKPT=EXTPNTR%IKPT
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
        deALLOCATE(CWORK1)
        IF(EXTPNTR%TIM) THEN
          DO IR=1,NRL
            VAL(IR)=AIMAG(CWORK2(IR))
          ENDDO
        ELSE
          DO IR=1,NRL
            VAL(IR)=REAL(CWORK2(IR))
          ENDDO
        END IF
        DEALLOCATE(CWORK2)
!
!     ================================================================
!     ==  GET PROJECTIONS                                           ==
!     ================================================================
      ELSE IF(ID.EQ.'<PSPSI|PRO>') THEN
        IKPT=EXTPNTR%IKPT
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
        NPRO=MAP%NPRO
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
! DO IPRO=1,NPRO
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
        IF(EXTPNTR%TIM) THEN
          DO LMN=1,LMNX
            VAL(LMN)=AIMAG(CWORK1(LMN))
          ENDDO
        ELSE
          DO LMN=1,LMNX
            VAL(LMN)=REAL(CWORK1(LMN))
          ENDDO
        END IF
        DEALLOCATE(CWORK1)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$GETR8A')
      END IF
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SETI4(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'SPINORDIM') THEN
        NDIM=VAL
      ELSE IF(ID.EQ.'IKPT') THEN
        EXTPNTR%IKPT=VAL
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
!     ..................................................................
      SUBROUTINE WAVES$GETI4(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'SPINORDIM') THEN
        VAL=NDIM
      ELSE IF(ID.EQ.'NSPIN') THEN
        VAL=NSPIN
      ELSE IF(ID.EQ.'NB') THEN
        IF(EXTPNTR%ISPIN.EQ.0) EXTPNTR%ISPIN=1
        IF(EXTPNTR%IKPT.EQ.0)  EXTPNTR%IKPT=1
        CALL WAVES_SELECTWV(EXTPNTR%IKPT,EXTPNTR%ISPIN)
        VAL=THIS%NB
      ELSE IF(ID.EQ.'NKPT') THEN
        VAL=NKPT
      ELSE IF(ID.EQ.'NR1') THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR1',VAL)
      ELSE IF(ID.EQ.'NR1L') THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR1L',VAL)
      ELSE IF(ID.EQ.'NR2') THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR2',VAL)
      ELSE IF(ID.EQ.'NR3') THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR3',VAL)
      ELSE IF(ID.EQ.'NR1START') THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR1START',VAL)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SETL4(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'STOP') THEN
        TSTOP=VAL
      ELSE IF(ID.EQ.'SAFEORTHO') THEN
        TSAFEORTHO=VAL
      ELSE IF(ID.EQ.'SWAPSTATES') THEN
        TSWAPSTATES=VAL
      ELSE IF(ID.EQ.'RANDOMIZE') THEN
        TRANDOM=VAL
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
!     ..................................................................
      SUBROUTINE WAVES$GETL4(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      LOGICAL(4)  ,INTENT(OUT) :: VAL
!     ******************************************************************
      IF(ID.EQ.'RAWSTATES') THEN
        VAL=TRAWSTATES
      ELSE IF(ID.EQ.'TIM') THEN
        VAL=EXTPNTR%TIM
      ELSE IF(ID.EQ.'SAFEORTHO') THEN
        VAL=TSAFEORTHO
      ELSE IF(ID.EQ.'SWAPSTATES') THEN
        VAL=TSWAPSTATES
      ELSE IF(ID.EQ.'WRITERHO') THEN
        val=TWRITERHO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$getL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SETCH(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'RSTRTTYPE') THEN
        RSTRTTYPE=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES$SETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$INITIALIZE
!     ******************************************************************
!     **  GENERATE G-VECTORS AND OTHER INITIALIZATION                 **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL(4)             :: TCHK
      REAL(8)                :: RBAS(3,3) ! LATTICE VECTORS
      REAL(8)                :: GBAS(3,3) ! RECIPROCAL LATTICE VECTORS
      REAL(8)                :: CELLVOL   ! UNIT CELL  VOLUME
      REAL(8)   ,ALLOCATABLE :: XK(:,:)   ! K-POINTS IN RELATIVE COORDINATES
      REAL(8)   ,ALLOCATABLE :: G2(:)     ! G**2
      INTEGER(4)             :: LN,ISP,IKPT,ISPIN,IAT
      INTEGER(4)             :: NB
      INTEGER(4)             :: NBH
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NRL,NR1L,NR1START,NR2,NR3
      INTEGER(4)             :: LNXX
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: NFILO
      integer(4)             :: lnx
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$GVECTORS')
!     
!     ================================================================
!     ==  DETERMINE MAPPRO AND MAPBAREPRO                           ==
!     ================================================================
      CALL ATOMLIST$NATOM(MAP%NAT)
      CALL SETUP$NSPECIES(MAP%NSP)
      CALL SETUP$LNXX(MAP%LNXX)
      ALLOCATE(MAP%ISP(MAP%NAT))
      ALLOCATE(MAP%LNX(MAP%NSP))
      ALLOCATE(MAP%LMNX(MAP%NSP))
      ALLOCATE(MAP%LOX(MAP%LNXX,MAP%NSP))
      CALL ATOMLIST$GETI4A('ISPECIES',0,MAP%NAT,MAP%ISP)
      MAP%LMX=0
      MAP%LOX(:,:)=0
      MAP%NBAREPRO=0
      DO ISP=1,MAP%NSP
        CALL SETUP$LNX(ISP,LNX)
        MAP%LNX(ISP)=lnx
        CALL SETUP$LOFLN(ISP,LNX,MAP%LOX(1:LNX,ISP))
        MAP%LMNX(ISP)=0
        DO LN=1,lnx
          MAP%LMNX(ISP)=MAP%LMNX(ISP)+2*MAP%LOX(LN,ISP)+1
          MAP%LMX=MAX(MAP%LMX,(MAP%LOX(LN,ISP)+1)**2)
        ENDDO
        MAP%NBAREPRO=MAP%NBAREPRO+lnx
      ENDDO
      MAP%NPRO=0
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
        MAP%NPRO=MAP%NPRO+MAP%LMNX(ISP)
      ENDDO
!     
!     ==================================================================
!     ==  ALLOCATE THISARRAY                                          ==
!     ==================================================================
      CALL DYNOCC$GETI4('NKPT',NKPT) 
      CALL DYNOCC$GETI4('NSPIN',NSPIN)
      ALLOCATE(THISARRAY(NKPT,NSPIN))
      DO IKPT=1,NKPT
        ALLOCATE(THISARRAY(IKPT,1)%GSET)
        DO ISPIN=2,NSPIN
          THISARRAY(IKPT,ISPIN)%GSET=>THISARRAY(IKPT,1)%GSET
        ENDDO
      ENDDO
!     
!     ==================================================================
!     ==  DEFINE REAL SPACE GRID                                      ==
!     ==================================================================
      CALL CELL$GETR8A('TREF',9,RBAS)
      CALL PLANEWAVE$DIVIDERGRIDONTASKS(EPWRHO,RBAS,NR1START,NR1L,NR2,NR3)
      CALL POTENTIAL$INITIALIZE(EPWRHO,NR1START,NR1L,NR2,NR3)
      MAP%NRL=NR1L*NR2*NR3
!     
!     ==================================================================
!     ==  INITIALIZE PLANE WAVE OBJECT                                ==
!     ==================================================================
      ALLOCATE(XK(3,NKPT))
      CALL DYNOCC$GETR8A('XK',3*NKPT,XK)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      IF(NDIM.EQ.1) THEN
        NDIMD=NSPIN
      ELSE
        NDIMD=NDIM**2
      END IF
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        WRITE(GSET%ID,FMT=*)IKPT
        GSET%ID='WAVE '//ADJUSTL(GSET%ID)
        IF(NDIM.EQ.1) THEN
          CALL PLANEWAVE$CHECKINVERSION(XK(1,IKPT),TINV)
        ELSE
          TINV=.FALSE.
        END IF
!print*,'tinv forced to be false in waves$initialize!!!'
!tinv=.false.
        CALL PLANEWAVE$INITIALIZE(GSET%ID,RBAS,XK(1,IKPT),TINV,EPWPSI &
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
        NULLIFY(GSET%DMPSI)
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
      DO IKPT=1,NKPT
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
!         == ALLOCATE PROJECTIONSONS ===================================
          ALLOCATE(THIS%PROJ(NDIM,NBH,MAP%NPRO))
!         == ALLOCATE WAVE FUNCTIONS ===================================
          NGL=GSET%NGL
          ALLOCATE(THIS%PSI0(NGL,NDIM,NBH))
          ALLOCATE(THIS%PSIM(NGL,NDIM,NBH))
          ALLOCATE(G2(NGL))
          CALL PLANEWAVE$GETR8A('G2',NGL,G2)
          CALL WAVES_INITIALIZERANDOM(NGL,NDIM,NBH,G2,THIS%PSI0)
          DEALLOCATE(G2)
          THIS%PSIM=THIS%PSI0
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
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$REPORT(NFILO)
      ENDDO
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
!     ..................................................................
      SUBROUTINE WAVES$ETOT()
!     ******************************************************************
!     **                                                              **
!     **  EVALUATE PS KINETIC ENERGY IN G-SPACE                       **
!     **  EVALUATE NUMBER OF ELECTRONS IN G-SPACE                     **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (2000)***
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      COMPLEX(8)             :: CSUM,CSUM1,CSUM2
      INTEGER(4)             :: LMRXX
      REAL(8)   ,ALLOCATABLE :: QLM(:,:)  !(LMRXX) MULTIPOLE MOMENTS
      REAL(8)   ,ALLOCATABLE :: VQLM(:,:) !(LMRXX) MULTIPOLE POTENTIALS
      REAL(8)   ,ALLOCATABLE :: RHO(:,:)  ! CHARGE DENSITY
      REAL(8)   ,ALLOCATABLE :: RHO1(:,:) ! AUXILIARY DENSITY ARRAY
      COMPLEX(8),ALLOCATABLE :: PROJ(:,:,:)     ! <PRO|PSI0>
      REAL(8)   ,ALLOCATABLE :: DENMAT(:,:,:,:) ! 1CENTER DENSITY MATRIX
      REAL(8)   ,ALLOCATABLE :: DENMATI(:,:,:,:) ! 1CENTER DENSITY MATRIX
!                     (IMAGINARY PART FOR SPIN ORBIT AND CURRENT DENSITIES)
      REAL(8)   ,ALLOCATABLE :: DH(:,:,:,:)     ! 1CENTER HAMILTONIAN
      REAL(8)   ,ALLOCATABLE :: DO(:,:,:,:)     ! 1CENTER OVERLAP
      REAL(8)   ,ALLOCATABLE :: DENMAT1(:,:,:)  ! 1CENTER DENSITY MATRIX
      REAL(8)   ,ALLOCATABLE :: DENMATI1(:,:,:) ! 1CENTER DENSITY MATRIX (IMAG)
      REAL(8)   ,ALLOCATABLE :: DH1(:,:,:)      ! 1CENTER DENSITY MATRIX
      REAL(8)   ,ALLOCATABLE :: DOV1(:,:,:)      ! 1CENTER DENSITY MATRIX
      COMPLEX(8),ALLOCATABLE :: DEDPRO(:,:)     ! DE/DPRO
      COMPLEX(8),ALLOCATABLE :: DEDPROJ(:,:,:)  ! DE/DPROJ
      COMPLEX(8),ALLOCATABLE :: HPROJ(:,:,:)    ! DH*PROJ
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      REAL(8)   ,ALLOCATABLE :: R(:,:)
      REAL(8)   ,ALLOCATABLE :: FORCE(:,:)
      REAL(8)   ,ALLOCATABLE :: FORCET(:,:)
      REAL(8)                :: FORCE1(3)
      REAL(8)                :: STRESS1(3,3),STRESS(3,3)
      REAL(8)                :: STRESST(3,3)
      REAL(8)                :: RBAS(3,3) ! real space lattice vectors
      REAL(8)                :: gbas(3,3) ! reciprocal space lattice vectors
      REAL(8)                :: RHOB      ! background density
      REAL(8)                :: potb,potb1 ! average electrostatic potential
      REAL(8)   ,PARAMETER   :: TINY=1.D-300
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NRL
      INTEGER(4)             :: LMRX
      REAL(8)                :: EKIN,EKIN1
      INTEGER(4)             :: IKPT,ISPIN,IAT,ISP,IBH,IB,IG,IR
      INTEGER(4)             :: L,M,LN,I,J,IJ
      INTEGER(4)             :: IPRO,IPRO1,IPRO2,IBPRO
      INTEGER(4)             :: LMN1,LMN2
      INTEGER(4)             :: LMNXX,LNX
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NBH,NBX,NB
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: THISTASK,NTASKS
      REAL(8)                :: svar,SVAR1,SVAR2
      REAL(8)   ,ALLOCATABLE :: DO1(:,:)
      COMPLEX(8),ALLOCATABLE :: HAMILTON(:,:)
      COMPLEX(8),ALLOCATABLE :: EIGR(:)
      REAL(8)   ,ALLOCATABLE :: EIG(:,:,:)
      REAL(8)                :: GWEIGHT
      REAL(8)   ,ALLOCATABLE :: G2(:)      !G**2
      REAL(8)   ,ALLOCATABLE :: GVEC(:,:)  !GI
      REAL(8)   ,ALLOCATABLE :: GIJ(:,:)   !GI*GJ/G2 (only filled if tstress!)
      LOGICAL(4)             :: TINV
      LOGICAL(4)             :: TSTRESS
      LOGICAL(4)             :: TFORCE
      LOGICAL(4)             :: TCHK
      LOGICAL(4)             :: TSO=.FALSE. ! EVALUATE DENMATI
      REAL(8)                :: PI
      COMPLEX(8),ALLOCATABLE :: QMAT(:,:)   
      INTEGER(4)             :: NFILO
      INTEGER(4)             :: L1,L2,LN1,LN2
      integer(4)             :: ngamma
!     ******************************************************************      
                             CALL TRACE$PUSH('WAVES$ETOT')
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     == CHECK CONSISTENCY WITH OCCUPATIONS OBJECT                    ==
!     ==================================================================
      CALL DYNOCC$GETL4('DYN',TCHK)
      IF(TCHK) THEN
        IF(TSAFEORTHO) THEN
          CALL ERROR$MSG('FOR DYNAMICAL OCCUPATIONS, SAFEORTHO MUST BE FALSE')
          CALL ERROR$STOP('WAVES$ETOT')
        END IF
      END IF
!
!     ==================================================================
!     == SWITCHES FOR FORCE AND STRESS CALCULATION                    ==
!     ==================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      CALL ATOMS$GETL4('MOVE',TFORCE)
      TSTRESS=TSTRESS.OR.TSTRESSX
      TFORCE=TFORCE.OR.TFORCEX
      IF(TSTRESS)CALL POTENTIAL$SETL4('STRESS',TSTRESS)
      IF(TFORCE)CALL POTENTIAL$SETL4('FORCE',TFORCE)
!
!     ==================================================================
!     == COLLECT VARIABLES                                            ==
!     ==================================================================
      NAT=MAP%NAT
!     == UNIT CELL =====================================================
      CALL CELL$GETR8A('T0',9,RBAS)
!     == ATOMIC POSITIONS ==============================================
      ALLOCATE(R(3,NAT))
      ALLOCATE(FORCE(3,NAT))
      FORCE(:,:)=0.D0
      STRESS=0.D0
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
!     == OCCUPATIONS ===================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPT,NSPIN))
      ALLOCATE(EIG(NBX,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NBX*NKPT*NSPIN,OCC)
!
!     ==================================================================
!     == INITIALIZE GSET: YLM, PRO, EIGR                              ==
!     ==================================================================
      IF(TFIRST.OR.TSTRESS) THEN
        CALL WAVES_UPDATEGSET()
      END IF
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(ASSOCIATED(THIS%EIGVAL))DEALLOCATE(THIS%EIGVAL)
          IF(ASSOCIATED(THIS%EIGVEC))DEALLOCATE(THIS%EIGVEC)
          NGL=GSET%NGL
          NB=THIS%NB
          NBH=THIS%NBH
          IF(TFIRST) THEN
            IF(TRANDOM) THEN
              ALLOCATE(G2(NGL))
              CALL PLANEWAVE$GETR8A('G2',NGL,G2)
              CALL WAVES_RANDOMIZE(NGL,NDIM,NBH,AMPRANDOM,G2,THIS%PSIM)
              DEALLOCATE(G2)
            END IF
!call waves_comparepsi('before gramschmidt',ngl,ndim,nbh,this%psi0,THIS%psim)
            CALL WAVES_GRAMSCHMIDT(MAP,GSET,NAT,R,NGL,NDIM,NBH,NB,THIS%PSI0)
            CALL WAVES_GRAMSCHMIDT(MAP,GSET,NAT,R,NGL,NDIM,NBH,NB,THIS%PSIM)
!call waves_comparepsi('after gramschmidt',ngl,ndim,nbh,this%psi0,THIS%psim)
!
!           ============================================================
!           == PROJECTIONS (REPLACEMENT BY PREDICTION IS FRAGILE!)    ==
!           ============================================================
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,MAP%NPRO &
     &                            ,THIS%PSI0,THIS%PROJ)
          END IF
        ENDDO
      ENDDO
      IF(TFIRST) TFIRST=.FALSE.
                              CALL TIMING$CLOCKON('WAVES$ETOT')
!
!     ==================================================================
!     == KINETIC ENERGY                                               ==
!     ==================================================================
CALL TIMING$CLOCKON('W:EKIN')
      EKIN=0.D0
      DO IKPT=1,NKPT
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
      CALL ENERGYLIST$SET('PS  KINETIC',EKIN)
      CALL ENERGYLIST$ADD('AE  KINETIC',EKIN)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EKIN)
CALL TIMING$CLOCKOFF('W:EKIN')
!D WRITE(*,FMT='("KIN STRESS ",3F15.7)')STRESS(1,:)
!D WRITE(*,FMT='("KIN STRESS ",3F15.7)')STRESS(2,:)
!D WRITE(*,FMT='("KIN STRESS ",3F15.7)')STRESS(3,:)
!
!     ==================================================================
!     == ONE-CENTER DENSITY MATRICES                                  ==
!     == SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                      ==
!     == SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)               ==
!     == NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z) ==
!     ==================================================================
CALL TRACE$PASS('BEFORE ONE-CENTER DENSITY MATRIX')     
CALL TIMING$CLOCKON('W:1CD')
      NAT=MAP%NAT
      LMNXX=0
      DO ISP=1,MAP%NSP
        LMNXX=MAX(LMNXX,MAP%LMNX(ISP))
      ENDDO
      ALLOCATE(DENMAT(LMNXX,LMNXX,NDIMD,NAT))
      DENMAT(:,:,:,:)=0.D0
      IF(TSO) THEN
        ALLOCATE(DENMATI(LMNXX,LMNXX,NDIMD,NAT))
        DENMATI(:,:,:,:)=0.D0
      END IF
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NBH=THIS%NBH
          NB=THIS%NB
          IPRO=1
          DO IAT=1,NAT
            ISP=MAP%ISP(IAT)
            LMNX=MAP%LMNX(ISP)
            ALLOCATE(PROJ(NDIM,NBH,LMNX))
            PROJ(:,:,:)=THIS%PROJ(:,:,IPRO:IPRO-1+LMNX)
            ALLOCATE(DENMAT1(LMNX,LMNX,NDIMD))
            ALLOCATE(DENMATI1(LMNX,LMNX,NDIMD))
            CALL WAVES_DENMAT(NDIM,NBH,NB,LMNX,OCC(1,IKPT,ISPIN),PROJ &
     &                       ,DENMAT1,DENMATI1)
            IF(NDIM.EQ.1) THEN
              DENMAT(1:LMNX,1:LMNX,ISPIN,IAT) &
     &                   =DENMAT(1:LMNX,1:LMNX,ISPIN,IAT)+DENMAT1(:,:,1)
            ELSE
              DENMAT(1:LMNX,1:LMNX,:,IAT) &
     &                       =DENMAT(1:LMNX,1:LMNX,:,IAT)+DENMAT1(:,:,:)
              IF(TSO) THEN
                DENMATI(1:LMNX,1:LMNX,:,IAT) &
     &                       =DENMATI(1:LMNX,1:LMNX,:,IAT)+DENMATI1(:,:,:)
              END IF
            ENDIF
            DEALLOCATE(DENMATI1)
            DEALLOCATE(DENMAT1)
            DEALLOCATE(PROJ)
            IPRO=IPRO+LMNX
          ENDDO
        ENDDO
      ENDDO
      IF(NSPIN.EQ.2) THEN
        DO IAT=1,NAT
          ISP=MAP%ISP(IAT)
          LMNX=MAP%LMNX(ISP)
          DO LMN1=1,LMNX
             DO LMN2=1,LMNX
              SVAR1=DENMAT(LMN1,LMN2,1,IAT)
              SVAR2=DENMAT(LMN1,LMN2,2,IAT)
              DENMAT(LMN1,LMN2,1,IAT)=SVAR1+SVAR2
              DENMAT(LMN1,LMN2,2,IAT)=SVAR1-SVAR2
            ENDDO
          ENDDO
        ENDDO
      END IF
CALL TIMING$CLOCKOFF('W:1CD')
!
!     ==================================================================
!     == PSEUDO DENSITY STILL WITHOUT PSEUDOCORE                      ==
!     == SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                      ==
!     == SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)               ==
!     == NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z) ==
!     ==================================================================
CALL TRACE$PASS('BEFORE PSEUDO DENSITY')     
CALL TIMING$CLOCKON('W:PSRHO')
      NRL=MAP%NRL
      ALLOCATE(RHO(NRL,NDIMD))
      ALLOCATE(RHO1(NRL,NDIMD))
      RHO(:,:)=0.D0
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL WAVES_DENSITY(GSET%NGL,MAP%NRL,NDIM,THIS%NB,THIS%NBH &
     &             ,OCC(1,IKPT,ISPIN),THIS%PSI0,RHO1)
          IF(NDIM.EQ.1) THEN
            RHO(:,ISPIN)=RHO(:,ISPIN)+RHO1(:,1)
          ELSE
            RHO=RHO+RHO1
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(RHO1)
      IF(NSPIN.EQ.2) THEN
        DO IR=1,NRL
          SVAR1=RHO(IR,1)
          SVAR2=RHO(IR,2)
          RHO(IR,1)=SVAR1+SVAR2
          RHO(IR,2)=SVAR1-SVAR2
        ENDDO
      END IF
CALL TIMING$CLOCKOFF('W:PSRHO')
!
!     ==================================================================
!     ==================================================================
IF(1.EQ.0) THEN
  SVAR1=0.D0
  DO IB=1,NRL
    SVAR1=SVAR1+RHO(IB,1)
  ENDDO
  CALL PLANEWAVE$GETR8('RWEIGHT',SVAR2)
  SVAR1=SVAR1*SVAR2
  CALL MPE$COMBINE('+',SVAR1)
  PRINT*,'TOTAL CHARGE IN PSEUDO WAVE FUNCTIONS W/O PSCORE ',SVAR1
END IF
!
!     ==================================================================
!     == MULTIPOLE MOMENTS OF ONE-CENTER PART                         ==
!     == CONTAINS CORE AND PSEUDOCORE                                 ==
!     ==================================================================
CALL TRACE$PASS('BEFORE MULTIPOLE MOMENTS')     
CALL TIMING$CLOCKON('W:MOMENTS')
      CALL SETUP$LMRXX(LMRXX)
      ALLOCATE(QLM(LMRXX,NAT))
      QLM(:,:)=0.D0
      DO IAT=THISTASK,NAT,NTASKS   ! DISTRIBUTE WORK ACCROSS TASKS
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        CALL SETUP$LMRX(ISP,LMRX)
        ALLOCATE(DENMAT1(LMNX,LMNX,1))
        DENMAT1(:,:,1)=DENMAT(1:LMNX,1:LMNX,1,IAT)
        CALL AUGMENTATION$MOMENTS(ISP,LMNX,DENMAT1,LMRX,QLM(1,IAT))
        DEALLOCATE(DENMAT1)
      ENDDO
      CALL MPE$COMBINE('+',QLM)
PI=4.D0*DATAN(1.D0)
SVAR1=0.D0
DO IAT=1,NAT
  SVAR1=SVAR1+QLM(1,IAT)*SQRT(4.D0*PI)
ENDDO
PRINT*,'TOTAL CHARGE IN AUGMENTATION',SVAR1
CALL TIMING$CLOCKOFF('W:MOMENTS')
!
!     ==================================================================
!     == overwrite density to for fixed potential calc.               ==
!     ==================================================================
      IF(TFIXRHO) THEN
        CALL WAVES_FIXRHOGET(NRL,NDIMD,LMRXX,NAT,QLM,RHO,DENMAT)
      END if
      IF(TWRITERHO) then
        CALL WAVES_FIXRHOSET(NRL,NDIMD,LMRXX,NAT,QLM,RHO,DENMAT)
      END IF
!
!     ==================================================================
!     == POTENTIAL                                                    ==
!     ==================================================================
      ALLOCATE(FORCET(3,NAT))
      ALLOCATE(VQLM(LMRXX,NAT))
      CALL POTENTIAL$VOFRHO(NRL,NDIMD,RHO,LMRXX,NAT,QLM,VQLM &
     &                     ,R,FORCET,RBAS,STRESST,RHOB)
      FORCE=FORCE+FORCET
      STRESS=STRESS+STRESST
      DEALLOCATE(QLM)
      DEALLOCATE(FORCET)
!
!     ==================================================================
!     == AUGMENTATION                                                 ==
!     ================================================================== 
CALL TRACE$PASS('BEFORE AUGMENTATION')     
CALL TIMING$CLOCKON('W:SPHERE')
      ALLOCATE(DH(LMNXX,LMNXX,NDIMD,NAT))
      DH(:,:,:,:)=0.D0
      ALLOCATE(DO(LMNXX,LMNXX,NDIMD,NAT))
      DO(:,:,:,:)=0.D0
      potb=0
      DO IAT=THISTASK,NAT,NTASKS   ! DISTRIBUTE WORK ACCROSS TASKS
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        CALL SETUP$LMRX(ISP,LMRX)
        ALLOCATE(DENMAT1(LMNX,LMNX,NDIMD))
        DENMAT1(:,:,:)=DENMAT(1:LMNX,1:LMNX,:,IAT)
        ALLOCATE(DENMATI1(LMNX,LMNX,NDIMD))
        IF(TSO) THEN
          DENMATI1(:,:,:)=DENMATI(1:LMNX,1:LMNX,:,IAT)
        ELSE
          DENMATI1(:,:,:)=0.D0
        END IF
        ALLOCATE(DH1(LMNX,LMNX,NDIMD))
        ALLOCATE(DOV1(LMNX,LMNX,NDIMD))
        CALL AUGMENTATION$SPHERE(ISP,IAT,LMNX,NDIMD,DENMAT1,DENMATI1 &
     &               ,LMRX,VQLM(1,IAT),RHOB,potb1,DH1,DOV1)
        potb=potb+potb1
        DH(1:LMNX,1:LMNX,:,IAT)=DH1(:,:,:)
        DO(1:LMNX,1:LMNX,:,IAT)=DOV1(:,:,:)
        DEALLOCATE(DH1)
        DEALLOCATE(DOV1)
        DEALLOCATE(DENMATI1)
        DEALLOCATE(DENMAT1)
      ENDDO
      DEALLOCATE(VQLM)
      CALL MPE$COMBINE('+',DH)
      CALL MPE$COMBINE('+',DO)
      CALL MPE$COMBINE('+',potb)
!     ==  SPREAD INFO FROM SPHERE OVER ALL TASKS AND UPDATE ENERGYLIST
      CALL AUGMENTATION$SYNC
!
!     ==  subtract average electrostatic potential =====================
!     ==  to account for background density ============================
!     ==  potb=-1/omega*int d^3r(\tilde{v}+v^1-\tilde{v}^1) ============
!     ==  note that \int d^3r \tilde{v}=0
      call gbass(rbas,gbas,svar)
      potb=potb/svar
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
     &                                 +POTB*DO(LMN1+M,LMN2+M,1,IAT)
              ENDDO
            END IF
            LMN2=LMN2+2*L2+1
          ENDDO
          LMN1=LMN1+2*L1+1
        ENDDO
      ENDDO
CALL TIMING$CLOCKOFF('W:SPHERE')
!
!     ==================================================================
!     == COMMUNICATE DATA WITH OPTICS MODULE                          ==
!     ==================================================================
!     CALL OPTICS$GETL4('ON',TCHK)
!     IF(TCHK) THEN
!        CALL MPE$COMBINE('+',DO) ! DO CURRENTLY USED ONLY FOR OPTICS
!        IF(NDIM.EQ.2) THEN
!          CALL ERROR$MSG('OPTICS CODE DOES NOT WORK WITH NONCOLLINEAR MODE')
!          CALL ERROR$STOP('WAVES$ETOT') 
!        END IF
!        CALL OPTICS$WRITEOUT_PARTIALS(NSPIN,NAT,LMNXX,DH,DO)
!      END IF
      DEALLOCATE(DO)  ! DO NOT USED YET....; EXCEPT IN OPTICS
!
!     ==================================================================
!     == PRINTOUT FOR TESTING 
!     ==================================================================
IF(NSPIN.EQ.2) THEN
  DO IAT=1,NAT
    ISP=MAP%ISP(IAT)
    LMNX=MAP%LMNX(ISP)
    DO ISPIN=1,NDIMD
      WRITE(*,FMT='("D ======== ATOM ",I2," SPIN ",I2," =========")')IAT,ISPIN
      LMN1=0.D0
      DO LN=1,MAP%LNX(ISP)
        L=MAP%LOX(LN,ISP)
        IF(L.NE.2) THEN
          LMN1=LMN1+2*L+1
          CYCLE
        END IF
        WRITE(*,FMT='("=",I1,10E10.3)') &
             LN,(DENMAT(LMN1+M,LMN1+M,ISPIN,IAT),M=1,2*L+1)
        LMN1=LMN1+2*L+1
      ENDDO
    ENDDO
    DO ISPIN=1,NDIMD
      WRITE(*,FMT='("H======== ATOM",I2," SPIN ",I2," =========")')IAT,ISPIN
      LMN1=0.D0
      DO LN=1,MAP%LNX(ISP)
        L=MAP%LOX(LN,ISP)
        IF(L.NE.2) THEN
          LMN1=LMN1+2*L+1
          CYCLE
        END IF
        WRITE(*,FMT='("=",I2,10E10.3)') &
             LN,(DH(LMN1+M,LMN1+M,ISPIN,IAT),M=1,2*L+1)
        LMN1=LMN1+2*L+1
      ENDDO
    ENDDO
  ENDDO 
ELSE
  DO IAT=1,NAT
    ISP=MAP%ISP(IAT)
    LMNX=MAP%LMNX(ISP)
    DO ISPIN=1,NDIMD
!      WRITE(*,FMT='("D======== ATOM ",I2," SPIN ",I2," =========")')IAT,ISPIN
      LMN1=0.D0
      DO LN=1,MAP%LNX(ISP)
        L=MAP%LOX(LN,ISP)
        IF(L.NE.2) THEN
          LMN1=LMN1+2*L+1
          CYCLE
        END IF
!        WRITE(*,FMT='("=",I2,10E10.3)') &
!             LN,(DENMAT(LMN1+M,LMN1+M,ISPIN,IAT),M=1,2*L+1)
        LMN1=LMN1+2*L+1
      ENDDO
    ENDDO
    DO ISPIN=1,NDIMD
!      WRITE(*,FMT='("H ======== ATOM ",I2," SPIN ",I2," =========")')IAT,ISPIN
      LMN1=0.D0
      DO LN=1,MAP%LNX(ISP)
        L=MAP%LOX(LN,ISP)
        IF(L.NE.2) THEN
!          LMN1=LMN1+2*L+1
!          CYCLE
        END IF
!        WRITE(*,FMT='("=",I2,10E10.3)') &
!             LN,(DH(LMN1+M,LMN1+M,ISPIN,IAT),M=1,2*L+1)
        LMN1=LMN1+2*L+1
      ENDDO
    ENDDO
  ENDDO 
END IF
      DEALLOCATE(DENMAT)
      IF(TSO)DEALLOCATE(DENMATI)
!
!     ==================================================================
!     == DECOMPOSE INTO SPIN-UP AND SPIN DOWN FOR NSPIN=2             ==
!     ==================================================================
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
              SVAR1=DH(LMN1,LMN2,1,IAT)
              SVAR2=DH(LMN1,LMN2,2,IAT)
              DH(LMN1,LMN2,1,IAT)=SVAR1+SVAR2
              DH(LMN1,LMN2,2,IAT)=SVAR1-SVAR2
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==================================================================
!     == FORCES AND STRESSES                                          ==
!     ==================================================================
CALL TRACE$PASS('BEFORE FORCES AND STRESSES')     
!CHECK IF ATOMS ARE MOVED
WRITE(*,FMT='("BEFORE 1C STRESS ",3F15.7)')STRESS(1,:)
WRITE(*,FMT='("BEFORE 1C STRESS ",3F15.7)')STRESS(2,:)
WRITE(*,FMT='("BEFORE 1C STRESS ",3F15.7)')STRESS(3,:)
CALL TIMING$CLOCKON('W:FORCE')
      IF(TFORCE.OR.TSTRESS) THEN
        DO IKPT=1,NKPT
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            CALL PLANEWAVE$SELECT(GSET%ID)
            CALL PLANEWAVE$GETL4('TINV',TINV)
            CALL PLANEWAVE$GETR8('GWEIGHT',GWEIGHT)
            NB=THIS%NB
            NBH=THIS%NBH
            NGL=GSET%NGL
            ALLOCATE(GVEC(3,NGL))
            CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
            ALLOCATE(GIJ(6,NGL))
            IF(TSTRESS) THEN
              DO IG=1,NGL
                SVAR=GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2
                SVAR=1.D0/(SVAR+TINY)
                IJ=0
                DO I=1,3
                  DO J=I,3
                    IJ=IJ+1
                    GIJ(IJ,IG)=GVEC(I,IG)*GVEC(J,IG)*SVAR
                  ENDDO
                ENDDO
              ENDDO
            else
              gij(:,:)=0
            END IF
            IPRO=1
            DO IAT=1,NAT
              ISP=MAP%ISP(IAT)
              IBPRO=1+SUM(MAP%LNX(1:ISP-1))
              LMNX=MAP%LMNX(ISP)
              LNX=MAP%LNX(ISP)
!             == DEDPROJ= (DH-LAMBDA*DO)<P|PSPSI> ======================
              ALLOCATE(DO1(LNX,LNX))
              CALL SETUP$1COVERLAP(ISP,LNX,DO1)
              ALLOCATE(DEDPROJ(NDIM,NBH,LMNX))
              ALLOCATE(DH1(LMNX,LMNX,NDIM**2))
              if(ndim.eq.1) then
                DH1(:,:,1)=DH(1:LMNX,1:LMNX,ispin,IAT)
              else 
                DH1(:,:,:)=DH(1:LMNX,1:LMNX,:,IAT)
              end if
              CALL WAVES_DEDPROJ(NDIM,NBH,NB,LNX,MAP%LOX(1:lnx,ISP),LMNX &
     &                       ,OCC(:,IKPT,ISPIN) &
     &                       ,THIS%PROJ(:,:,ipro:ipro+lmnx-1),DH1,DO1 &
     &                       ,THIS%RLAM0,DEDPROJ)
              DEALLOCATE(DH1)
              DEALLOCATE(DO1)
!             == |DEDPRO>=|PSPSI>DEDPROJ ===============================
              ALLOCATE(DEDPRO(NGL,LMNX))
              CALL WAVES_DEDPRO(GSET%TINV,NGL,NDIM,NBH,THIS%PSI0,LMNX,DEDPROJ,DEDPRO)
              DEALLOCATE(DEDPROJ)
!             == DE= <DPRO|DEDPRO> =====================================
              ALLOCATE(EIGR(NGL))
              CALL PLANEWAVE$STRUCTUREFACTOR(R(1,IAT),NGL,EIGR)
              CALL WAVES_PROFORCE(LNX,LMNX,MAP%LOX(1:lnx,ISP),NGL,GWEIGHT,GVEC,GIJ &
     &                   ,GSET%PRO(:,IBPRO:ibpro+lnx-1),GSET%DPRO(:,IBPRO:ibpro+lnx-1) &
     &                   ,MAP%LMX,GSET%YLM,GSET%SYLM &
     &                   ,EIGR,DEDPRO,FORCE1,TSTRESS,STRESS1)
              CALL MPE$COMBINE('+',FORCE1)
              IF(TSTRESS)CALL MPE$COMBINE('+',STRESS1)
              DEALLOCATE(EIGR)
              DEALLOCATE(DEDPRO)
              FORCE(:,IAT)=FORCE(:,IAT)+FORCE1(:)
              STRESS(:,:) =STRESS(:,:) +STRESS1(:,:)
              IPRO=IPRO+LMNX
            ENDDO
            DEALLOCATE(GIJ)
            DEALLOCATE(GVEC)
          ENDDO
        ENDDO
      END IF
CALL TIMING$CLOCKOFF('W:FORCE')

!
!     ==================================================================
!     ==  EVALUATE H*PSI                                              ==
!     ================================================================== 
CALL TRACE$PASS('BEFORE HPSI')     
CALL TIMING$CLOCKON('W:HPSI')
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)      
!         ==============================================================
!         ==  MULTIPLY WAVE FUNCTION WITH LOCAL POTENTIAL AND         ==
!         ==  ADD KINETIC ENERGY                                      ==
!         ==============================================================
          NGL=GSET%NGL
          NRL=MAP%NRL
          NBH=THIS%NBH
          IF(.NOT.ASSOCIATED(THIS%HPSI))ALLOCATE(THIS%HPSI(NGL,NDIM,NBH))
!         == NOTE: THE ARRAY RHO CONTAINS THE POTENTIAL ================
CALL TIMING$CLOCKON('W:HPSI.1')
          CALL WAVES_VPSI(GSET,NGL,NDIM,NBH,NRL,THIS%PSI0,RHO(1,ISPIN) &
     &                   ,THIS%HPSI)
CALL TIMING$CLOCKOFF('W:HPSI.1')
!
!         ==============================================================
!         ==  EVALUATE  DH<P|PSI>                                     ==
!         ==============================================================
CALL TIMING$CLOCKON('W:HPSI.2')
          ALLOCATE(HPROJ(NDIM,NBH,MAP%NPRO))
          HPROJ(:,:,:)=(0.D0,0.D0)
          IPRO=1
          DO IAT=1,NAT
            ISP=MAP%ISP(IAT)
            LMNX=MAP%LMNX(ISP)
            ALLOCATE(DH1(LMNX,LMNX,NDIM**2))
            if(ndim.eq.1) then
              DH1(:,:,1)=DH(1:LMNX,1:LMNX,ispin,IAT)
            else 
              DH1(:,:,:)=DH(1:LMNX,1:LMNX,:,IAT)
            end if
            CALL WAVES_HPROJ(NDIM,NBH,LMNX &
     &           ,dh1,THIS%PROJ(:,:,ipro:ipro+lmnx-1),HPROJ(:,:,ipro:ipro+lmnx-1))
            deallocate(dh1)
            IPRO=IPRO+LMNX
          ENDDO
CALL TIMING$CLOCKOFF('W:HPSI.2')
!
!         ==============================================================
!         ==  ADD  |P>DH<P|PSI>                                       ==
!         ==============================================================
CALL TIMING$CLOCKON('W:HPSI.3')
          CALL WAVES_ADDPRO(MAP,GSET,NAT,R,NGL,NDIM,NBH,MAP%NPRO &
     &                     ,THIS%HPSI,HPROJ)
CALL TIMING$CLOCKOFF('W:HPSI.3')
          DEALLOCATE(HPROJ)
        ENDDO
      ENDDO
      DEALLOCATE(RHO)
      DEALLOCATE(DH)
CALL TIMING$CLOCKOFF('W:HPSI')
!
!     ==================================================================
!     ==  EVALUATE ENERGY EXPECTATION VALUES                          ==
!     ================================================================== 
CALL TIMING$CLOCKON('W:EXPECT')
      ALLOCATE(HAMILTON(2,2))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          NB=THIS%NB
          NBH=THIS%NBH
          NGL=GSET%NGL
          IF(.NOT.ASSOCIATED(THIS%EXPECTVAL))ALLOCATE(THIS%EXPECTVAL(NB))
          DO IB=1,NBH
            IF(GSET%TINV) THEN
              CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,1,2 &
     &              ,THIS%PSI0(:,:,IB),THIS%HPSI(:,:,IB),HAMILTON)
              EIG(2*IB-1,IKPT,ISPIN)=REAL(HAMILTON(1,1))
              EIG(2*IB  ,IKPT,ISPIN)=REAL(HAMILTON(2,2))
            ELSE 
              CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,1,1 &
     &            ,THIS%PSI0(:,:,IB),THIS%HPSI(:,:,IB),HAMILTON)
              EIG(IB,IKPT,ISPIN)=REAL(HAMILTON(1,1))
            END IF
          ENDDO
          THIS%EXPECTVAL(:)=EIG(1:NB,IKPT,ISPIN)
!PRINT*,'EIG ',EIG(:,IKPT,ISPIN)
        ENDDO
      ENDDO
      DEALLOCATE(HAMILTON)
!     == OCCUPATIONS ===================================================
      CALL DYNOCC$SETR8A('EPSILON',NB*NKPT*NSPIN,EIG)
!
!     ==================================================================
!     ==  EVALUATE <S^2>                                              ==
!     ==================================================================
      IF ((NDIM.EQ.2).AND.THAMILTON) THEN   ! THAMILTON SHOULD BE REPLACED
         CALL TIMING$CLOCKON('S2')
         DO IKPT=1,NKPT
            CALL WAVES_SELECTWV(IKPT,1)
            CALL PLANEWAVE$SELECT(GSET%ID)      
            NB=THIS%NB
            NBH=THIS%NBH
            NGL=GSET%NGL         
            ALLOCATE(QMAT(2*NB*NSPIN,2*NB*NSPIN))
            CALL WAVES_SPINOROVERLAP(NBH,NB,IKPT,QMAT)
            CALL WAVES_TOTALSPIN(NB,NKPT,IKPT,NSPIN,OCC,QMAT)
            DEALLOCATE(QMAT)
         END DO
         CALL TIMING$CLOCKOFF('S2')
         CALL FILEHANDLER$UNIT('PROT',NFILO)
         CALL WAVES$REPORTSPIN(NFILO)
      END IF
!
!     ==================================================================
!     ==  EVALUATE <PSI|H|PSI>                                        ==
!     ================================================================== 
      IF(THAMILTON) THEN
        DO IKPT=1,NKPT
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            CALL PLANEWAVE$SELECT(GSET%ID)      
            NB=THIS%NB
            NBH=THIS%NBH
            NGL=GSET%NGL
            ALLOCATE(HAMILTON(NB,NB))
            CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,NBH,NB,THIS%PSI0,THIS%HPSI,HAMILTON)
            IF(.NOT.ASSOCIATED(THIS%EIGVAL))ALLOCATE(THIS%EIGVAL(NB))
            IF(.NOT.ASSOCIATED(THIS%EIGVEC))ALLOCATE(THIS%EIGVEC(NB,NB))
            CALL LIB$DIAGC8(NB,HAMILTON,THIS%EIGVAL,THIS%EIGVEC)
            DEALLOCATE(HAMILTON)
          ENDDO
        ENDDO
      ELSE
        IF(ASSOCIATED(THIS%EIGVAL))DEALLOCATE(THIS%EIGVAL)
        IF(ASSOCIATED(THIS%EIGVEC))DEALLOCATE(THIS%EIGVEC)
      END IF
CALL TIMING$CLOCKOFF('W:EXPECT')
!
!     ==================================================================
!     ==  FIRST HALF OF WAVE FUNCTIONS KINETIC ENERGY                 ==
!     ================================================================== 
      CALL WAVES_WAVEKINETIC(WAVEEKIN1)
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ================================================================== 
!     == FORCES ========================================================
      ALLOCATE(FORCET(3,NAT))
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCET)
      FORCE=FORCET+FORCE
      CALL ATOMLIST$SETR8A('FORCE',0,3*NAT,FORCE)
      DEALLOCATE(FORCET)
!     == STRESS ========================================================
!D WRITE(*,FMT='("RBAS ",3F15.7)')RBAS(1,:)
!D WRITE(*,FMT='("RBAS ",3F15.7)')RBAS(2,:)
!D WRITE(*,FMT='("RBAS ",3F15.7)')RBAS(3,:)
!D WRITE(*,FMT='("AFTER WAVES STRESS ",3F15.7)')STRESS(1,:)
!D WRITE(*,FMT='("AFTER WAVES STRESS ",3F15.7)')STRESS(2,:)
!D WRITE(*,FMT='("AFTER WAVES STRESS ",3F15.7)')STRESS(3,:)
      CALL CELL$GETR8A('STRESS_I',9,STRESS1)
      STRESS=STRESS1-STRESS  ! IN THIS ROUTINE STRESS=+DE/DEPSILON!
      CALL CELL$SETR8A('STRESS_I',9,STRESS)
!     == DEALLOCATE ARRAYS =============================================
      DEALLOCATE(FORCE)
      DEALLOCATE(OCC)
      DEALLOCATE(R)
      DEALLOCATE(EIG)
                              CALL TIMING$CLOCKOFF('WAVES$ETOT')
                                   CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_HPROJ(NDIM,NB,LMNX,dh,PROJ,HPROJ)
!     ******************************************************************
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      USE WAVES_MODULE, ONLY : MAP_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NDIM
      INTEGER(4),INTENT(IN)  :: NB
      INTEGER(4),INTENT(IN)  :: LMNX
      REAL(8)   ,INTENT(IN)  :: DH(LMNX,LMNX,NDIM**2)
      COMPLEX(8),INTENT(IN)  :: PROJ(NDIM,NB,LMNX)
      COMPLEX(8),INTENT(OUT) :: HPROJ(NDIM,NB,LMNX)
      INTEGER(4)             :: IB
      INTEGER(4)             :: LMN1,LMN2
      REAL(8)                :: DHUPUP,DHDNDN
      COMPLEX(8)             :: DHUPDN,DHDNUP
!     ******************************************************************
      HPROJ(:,:,:)=(0.D0,0.D0)
!
!     ==============================================================
!     ==  EVALUATE  DH<P|PSI>                                     ==
!     ==============================================================
      IF(NDIM.EQ.1) THEN
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            DHUPUP=DH(LMN1,LMN2,1)
            DO IB=1,NB
              HPROJ(1,IB,LMN1)=HPROJ(1,IB,LMN1) &
     &                        +DHUPUP*PROJ(1,IB,LMN2)
            ENDDO
          ENDDO
        ENDDO
      ELSE 
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            DHUPUP=DH(LMN1,LMN2,1)+DH(LMN1,LMN2,4)
            DHDNDN=DH(LMN1,LMN2,1)-DH(LMN1,LMN2,4)
!           == HERE, DH=DH/DD AND NOT CONJG(DE/DD) !! 
!JSC+PEB    DHUPDN=CMPLX(DH(LMN1,LMN2,2),DH(LMN1,LMN2,3),8)
            DHUPDN=CMPLX(DH(LMN1,LMN2,2),-DH(LMN1,LMN2,3),8)
            DHDNUP=CONJG(DHUPDN)
            DO IB=1,NB
              HPROJ(1,IB,LMN1)=HPROJ(1,IB,LMN1) &
     &                        +DHUPUP*PROJ(1,IB,LMN2) &
     &                        +DHUPDN*PROJ(2,IB,LMN2)
              HPROJ(2,IB,LMN1)=HPROJ(2,IB,LMN1) &
     &                        +DHDNUP*PROJ(1,IB,LMN2) &
     &                        +DHDNDN*PROJ(2,IB,LMN2)
            ENDDO
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_VPSI(GSET,NGL,NDIM,NBH,NRL,PSI,V,HPSI)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES H*PSI WITHOUT THE AUGMENTATION PART               **
!     **  |HPSI>=(-0.5*NABLA**2+PS-V)|PSPSI(0)>                       **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
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
      INTEGER(4)             :: IC,IDIM1,IDIM2,IB,IR,IG,IDIM
      REAL(8)   ,ALLOCATABLE :: VUPUP(:),VDNDN(:)
      COMPLEX(8),ALLOCATABLE :: VUPDN(:)
      COMPLEX(8),ALLOCATABLE :: PSIOFR(:,:,:)
      COMPLEX(8)             :: PSIUP,PSIDN
      REAL(8)                :: SVAR
!     ******************************************************************
!
!     ==================================================================
!     ==  FOURIER TRANSFORM WAVE FUNCTIONS TO REAL SPACE              ==
!     ==================================================================
!     -- REMARK: CONSIDER BLOCKING THIS INTO SMALLER BUNCHES OF BANDS
      ALLOCATE(PSIOFR(NRL,NDIM,NBH))
      CALL PLANEWAVE$FFT('GTOR',NBH*NDIM,NGL,PSI,NRL,PSIOFR)
!
!     ==================================================================
!     ==  MULTIPLY WAVE FUNCTIONS WITH THE POTENTIAL                  ==
!     ==================================================================
      IF(NDIM.EQ.1) THEN
        DO IB=1,NBH
          DO IR=1,NRL
            PSIOFR(IR,1,IB)=V(IR,1)*PSIOFR(IR,1,IB)
          ENDDO
        ENDDO
      ELSE
        ALLOCATE(VUPUP(NRL))
        ALLOCATE(VDNDN(NRL))
        ALLOCATE(VUPDN(NRL))
        DO IR=1,NRL
          VUPUP(IR)=V(IR,1)+V(IR,4)
          VDNDN(IR)=V(IR,1)-V(IR,4)
          VUPDN(IR)=CMPLX(V(IR,2),-V(IR,3),8)
        ENDDO
        DO IB=1,NBH
          DO IR=1,NRL
            PSIUP=PSIOFR(IR,1,IB)
            PSIDN=PSIOFR(IR,2,IB)
            PSIOFR(IR,1,IB)=VUPUP(IR)*PSIUP+      VUPDN(IR) *PSIDN
            PSIOFR(IR,2,IB)=VDNDN(IR)*PSIDN+CONJG(VUPDN(IR))*PSIUP
          ENDDO
        ENDDO
        DEALLOCATE(VUPUP)
        DEALLOCATE(VUPDN)
        DEALLOCATE(VDNDN)
      END IF
!
!     ==================================================================
!     ==  FOURIER TRANSFORM WAVE FUNCTIONS TO RECIPROCAL SPACE SPACE  ==
!     ==================================================================
      CALL PLANEWAVE$FFT('RTOG',NBH*NDIM,NGL,HPSI,NRL,PSIOFR)
      DEALLOCATE(PSIOFR)
!
!     ==================================================================
!     ==  ADD KINETIC ENERGY CONTRIBUTION                             ==
!     ==================================================================
      ALLOCATE(G2(NGL))
      CALL PLANEWAVE$GETR8A('G2',NGL,G2)
      DO IB=1,NBH
        DO IDIM=1,NDIM
          DO IG=1,NGL
            HPSI(IG,IDIM,IB)=HPSI(IG,IDIM,IB)+0.5D0*G2(IG)*PSI(IG,IDIM,IB)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  ADD BUCKET POTENTIAL                                        ==
!     ==================================================================
      IF(TBUCKET) THEN
        DO IB=1,NBH
          DO IDIM=1,NDIM
            DO IG=1,NGL
              HPSI(IG,IDIM,IB)=HPSI(IG,IDIM,IB) &
     &                        +GSET%BUCKET(IG)*PSI(IG,IDIM,IB)
            ENDDO
          ENDDO
        ENDDO
      END IF
      DEALLOCATE(G2)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_DENMAT(NDIM,NBH,NB,LMNX,OCC,PROPSI,DENMAT,DENMATI)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATE ONE-CENTER DENSITY MATRIX                          **
!     **             <P_I|PSI_N>F_N<PSI_N|P_J>                        **
!     **  FROM PROJECTIONS <P|PSI> AND THE OCCUPATIONS                **
!     **                                                              **
!     **  FOR NDIM=2,NDIMD=4:  SPINOR WAVE FUNCTIONS                  **
!     **  THE DENSITY MATRIX ON RETURN CONTAINS THE TOTAL DENSITY     **
!     **  AND THE X,Y,Z SPIN DENSITY                                  **
!     **                                                              **
!     **                                                              **
!     **  SUPERWAVE FUNCTIONS ARE DEFINED AS: PSI=PSI1+I*PSI2         **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
!RELEASED: 8.OCT.99
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NDIM   ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN) :: NBH    ! #(WAVE FUNCTIONS
      INTEGER(4),INTENT(IN) :: NB     ! #(STATES)
      INTEGER(4),INTENT(IN) :: LMNX   ! #(PROJECTORS ON THIS SITE)
      REAL(8)   ,INTENT(IN) :: OCC(NB)! OCCUPATIONS
      COMPLEX(8),INTENT(IN) :: PROPSI(NDIM,NBH,LMNX) !<PRO|PSI>
      REAL(8)   ,INTENT(OUT):: DENMAT(LMNX,LMNX,NDIM**2)
      REAL(8)   ,INTENT(OUT):: DENMATI(LMNX,LMNX,NDIM**2)
      COMPLEX(8)            :: DENMAT1(LMNX,LMNX,NDIM,NDIM)
      COMPLEX(8)            :: FUNC(LMNX,NDIM)
      LOGICAL(4)            :: TINV
      INTEGER(4)            :: LMN1,LMN2,IDIM1,IDIM2,IB
      REAL(8)               :: SUM,SVAR1,SVAR2,SVAR
      COMPLEX(8)            :: CSVAR
      INTEGER(4)            :: NDIMD,IDIM
!     ******************************************************************
      NDIMD=NDIM**2
!
!     ==================================================================
!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT==
!     ==================================================================
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
!
!     ==================================================================
!     ==  SUM UP THE DENSITY MATRIX                                   ==
!     ==================================================================
      DENMAT1(:,:,:,:)=0.D0
      DO IB=1,NBH
!       == FIND OCCUPATION OF THE STATE ===============================
        IF(TINV) THEN
          SVAR1=0.5D0*(OCC(2*IB-1)+OCC(2*IB))
          SVAR2=0.5D0*(OCC(2*IB-1)-OCC(2*IB))
        ELSE
          SVAR1=OCC(IB)
        END IF
        DO LMN1=1,LMNX
          DO IDIM1=1,NDIM
            IF(TINV) THEN
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
!
!     ==================================================================
!     == MAP DENSITY MATRIX ONTO TOTAL AND SPIN DENSITY               ==
!     ==================================================================
      IF(NDIM.EQ.1) THEN  !== TOTAL DENSITY ===========================
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            DENMAT(LMN1,LMN2,1)=REAL(DENMAT1(LMN1,LMN2,1,1))
          ENDDO
        ENDDO
      ELSE IF(NDIM.EQ.2) THEN  !== TOTAL DENSITY, X,Y,Z SPIN DENSITY ==
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            CSVAR=DENMAT1(LMN1,LMN2,1,1)+DENMAT1(LMN1,LMN2,2,2)
            DENMAT(LMN1,LMN2,1) = REAL(CSVAR)
            DENMATI(LMN1,LMN2,1)=AIMAG(CSVAR)
            CSVAR=DENMAT1(LMN1,LMN2,1,2)+DENMAT1(LMN1,LMN2,2,1)
            DENMAT(LMN1,LMN2,2) = REAL(CSVAR)
            DENMATI(LMN1,LMN2,2)=AIMAG(CSVAR)
            CSVAR=DENMAT1(LMN1,LMN2,1,2)-DENMAT1(LMN1,LMN2,2,1)
            DENMAT(LMN1,LMN2,3) =-AIMAG(CSVAR)
            DENMATI(LMN1,LMN2,3)=+REAL(CSVAR)
            CSVAR=DENMAT1(LMN1,LMN2,1,1)-DENMAT1(LMN1,LMN2,2,2)
            DENMAT(LMN1,LMN2,4) = REAL(CSVAR)
            DENMATI(LMN1,LMN2,4)=AIMAG(CSVAR)
          ENDDO
        ENDDO
      END IF
!
!     ==================================================================
!     == SYMMETRIZE DENSITY MATRIX                                    ==
!     ==================================================================
      DO IDIM=1,NDIMD
        DO LMN1=1,LMNX
          DO LMN2=LMN1+1,LMNX
            SVAR=0.5D0*(DENMAT(LMN1,LMN2,IDIM)+DENMAT(LMN2,LMN1,IDIM))
            DENMAT(LMN1,LMN2,IDIM)=SVAR
            DENMAT(LMN2,LMN1,IDIM)=SVAR
            SVAR=0.5D0*(DENMATI(LMN1,LMN2,IDIM)-DENMATI(LMN2,LMN1,IDIM))
            DENMATI(LMN1,LMN2,IDIM)= SVAR
            DENMATI(LMN2,LMN1,IDIM)=-SVAR
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_EKIN(NGL,NDIM,NBH,NB,F,GWEIGHT,PSI,EKIN &
     &                         ,TSTRESS,STRESS,TBUCKET,BUCKET,DBUCKET)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATE PS KINETIC ENERGY IN G-SPACE                       **
!     **  EVALUATE NUMBER OF ELECTRONS IN G-SPACE                     **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    REQUIRES PLANEWAVE OBJECT TO BE SET PROPERLY              **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
      USE MPE_MODULE
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
      INTEGER(4)            :: IB,IG,IDIM,I,J
      INTEGER(4)            :: NDIMD
      REAL(8)               :: FP,FM
      REAL(8)   ,PARAMETER  :: DSMALL=1.D-12
      REAL(8)               :: EBUCKET
      REAL(8)               :: SVAR
!     ******************************************************************
!
!     ==================================================================
!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT==
!     ==================================================================
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
!     ==================================================================
!     ==  CALCULATE DMAT(G)=F(IB)*PSI^*(G,IB)*PSI(G,IB)               ==
!     ==================================================================
      NDIMD=(NDIM*(NDIM+1))/2
      ALLOCATE(DMAT(NGL))
      ALLOCATE(PSI1(NGL,NDIM))
      DMAT(:)=0.D0
      DO IB=1,NBH
!       == DETERMINE OCCUPATIONS =======================================
        IF(TINV) THEN
          FP=0.5D0*(F(2*IB-1)+F(2*IB))
          FM=0.5D0*(F(2*IB-1)-F(2*IB))
        ELSE
          FP=F(IB)
          FM=0.D0
        END IF
!
!       ================================================================
!       == GENERAL WAVE FUNCTION / FIRST PART OF SUPER WAVE FUNCTIONS ==
!       == <PSI_+|G^2|PSI_+> FOR SUPER WAVE FUNCTIONS                 ==
!       ================================================================
        IF(FP.EQ.0.D0.AND.FM.EQ.0.D0) CYCLE
        DO IDIM=1,NDIM
          DO IG=1,NGL
            DMAT(IG)=DMAT(IG) &
    &                 +FP*REAL(CONJG(PSI(IG,IDIM,IB))*PSI(IG,IDIM,IB))
          ENDDO
        ENDDO
!
!       ================================================================
!       ==  <PSI_+|G^2|PSI_->                                         ==
!       ================================================================
        IF(.NOT.TINV.OR.ABS(FM).LT.DSMALL) CYCLE
        DO IDIM=1,NDIM
          CALL PLANEWAVE$INVERTG(NGL,PSI(1,IDIM,IB),PSI1(1,IDIM))
        ENDDO
        DO IDIM=1,NDIM
          DO IG=1,NGL
            DMAT(IG)=DMAT(IG) &
     &                +FM*REAL(CONJG(PSI(IG,IDIM,IB))*PSI1(IG,IDIM))
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(PSI1)
!
!     ==================================================================
!     ==  CALCULATE KINETIC ENERGY AND STRESS                         ==
!     ==================================================================
      IF(TSTRESS) THEN
        ALLOCATE(GVEC(3,NGL))
        CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
!       == KINETIC ENERGY AND KINETIC STRESS ===========================
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
!       == ENERGY AND STRESS DUE TO THE BUCKET POTENTIAL ===============
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
!       == PARALLELIZE AND CLOSE DOWN ==================================
        STRESS(2,1)=STRESS(1,2)
        STRESS(3,1)=STRESS(1,3)
        STRESS(3,2)=STRESS(2,3)
        EKIN=EKIN*GWEIGHT
        STRESS=STRESS*GWEIGHT
        CALL MPE$COMBINE('+',EKIN)
        CALL MPE$COMBINE('+',STRESS)
        DEALLOCATE(GVEC)
        DEALLOCATE(DMAT)
      ELSE
        ALLOCATE(G2(NGL))
        CALL PLANEWAVE$GETR8A('G2',NGL,G2)
!       == KINETIC ENERGY ==============================================
        EKIN=0.D0
        DO IG=1,NGL
          EKIN=EKIN+G2(IG)*DMAT(IG)
        ENDDO
        EKIN=0.5D0*EKIN
!       == BUCKET POTENTIAL ============================================
        IF(TBUCKET) THEN
          EBUCKET=0.D0
          DO IG=1,NGL
            EBUCKET=EBUCKET+BUCKET(IG)*DMAT(IG)
          ENDDO
          EKIN=EKIN+EBUCKET
        END IF
!       ==  PARALLELIZE AND CLOSE DOWN =================================
        EKIN=EKIN*GWEIGHT
        CALL MPE$COMBINE('+',EKIN)
        DEALLOCATE(DMAT)
        DEALLOCATE(G2)
        STRESS(:,:)=0.D0
      END IF
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_EKIN_OLD(NGL,NDIM,NBH,NB,F,GWEIGHT,PSI,EKIN &
     &                         ,TSTRESS,STRESS,TBUCKET,BUCKET,DBUCKET)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATE PS KINETIC ENERGY IN G-SPACE                       **
!     **  EVALUATE NUMBER OF ELECTRONS IN G-SPACE                     **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    REQUIRES PLANEWAVE OBJECT TO BE SET PROPERLY              **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
      USE MPE_MODULE
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
      INTEGER(4)            :: IB,IG,IDIM,I,J,IDIM1,IDIM2
      INTEGER(4)            :: NDIMD
      REAL(8)               :: FP,FM
      REAL(8)   ,PARAMETER  :: DSMALL=1.D-12
      REAL(8)               :: EBUCKET
      REAL(8)               :: SVAR
!     ******************************************************************
      CALL ERROR$MSG('ROUTINE MARKED FOR DELETION. CONTAINS ERRORS!')
      CALL ERROR$STOP('WAVES_EKIN_OLD')
!
!     ==================================================================
!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT==
!     ==================================================================
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
!     ==================================================================
!     ==  CALCULATE DMAT(G)=F(IB)*PSI^*(G,IB)*PSI(G,IB)               ==
!     ==================================================================
      NDIMD=(NDIM*(NDIM+1))/2
      ALLOCATE(DMAT(NGL))
      ALLOCATE(PSI1(NGL,NDIM))
      DMAT(:)=0.D0
      DO IB=1,NBH
!       == DETERMINE OCCUPATIONS =======================================
        IF(TINV) THEN
          FP=0.5D0*(F(2*IB-1)+F(2*IB))
          FM=0.5D0*(F(2*IB-1)-F(2*IB))
        ELSE
          FP=F(IB)
          FM=0.D0
        END IF
!
!       ================================================================
!       == GENERAL WAVE FUNCTION / FIRST PART OF SUPER WAVE FUNCTIONS ==
!       == <PSI_+|G^2|PSI_+> FOR SUPER WAVE FUNCTIONS                 ==
!       ================================================================
        IF(FP.EQ.0.D0.AND.FM.EQ.0.D0) CYCLE
        DO IDIM1=1,NDIM
          DO IDIM2=IDIM1,NDIM
            SVAR=FP
            IF(IDIM1.NE.IDIM2)SVAR=2.D0*SVAR
            DO IG=1,NGL
              DMAT(IG)=DMAT(IG) &
    &                 +SVAR*REAL(CONJG(PSI(IG,IDIM1,IB))*PSI(IG,IDIM2,IB))
            ENDDO
          ENDDO
        ENDDO
!
!       ================================================================
!       ==  <PSI_+|G^2|PSI_->                                         ==
!       ================================================================
        IF(.NOT.TINV.OR.ABS(FM).LT.DSMALL) CYCLE
        DO IDIM1=1,NDIM
          CALL PLANEWAVE$INVERTG(NGL,PSI(1,IDIM1,IB),PSI1(1,IDIM1))
        ENDDO
        DO IDIM1=1,NDIM
          DO IDIM2=IDIM1,NDIM
            SVAR=FM
            IF(IDIM1.NE.IDIM2)SVAR=2.D0*SVAR
            DO IG=1,NGL
              DMAT(IG)=DMAT(IG) &
     &                +SVAR*REAL(CONJG(PSI(IG,IDIM1,IB))*PSI1(IG,IDIM2))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(PSI1)
!
!     ==================================================================
!     ==  CALCULATE KINETIC ENERGY AND STRESS                         ==
!     ==================================================================
      IF(TSTRESS) THEN
        ALLOCATE(GVEC(3,NGL))
        CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
!       == KINETIC ENERGY AND KINETIC STRESS ===========================
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
!       == ENERGY AND STRESS DUE TO THE BUCKET POTENTIAL ===============
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
!       == PARALLELIZE AND CLOSE DOWN ==================================
        STRESS(2,1)=STRESS(1,2)
        STRESS(3,1)=STRESS(1,3)
        STRESS(3,2)=STRESS(2,3)
        EKIN=EKIN*GWEIGHT
        STRESS=STRESS*GWEIGHT
        CALL MPE$COMBINE('+',EKIN)
        CALL MPE$COMBINE('+',STRESS)
        DEALLOCATE(GVEC)
        DEALLOCATE(DMAT)
      ELSE
        ALLOCATE(G2(NGL))
        CALL PLANEWAVE$GETR8A('G2',NGL,G2)
!       == KINETIC ENERGY ==============================================
        EKIN=0.D0
        DO IG=1,NGL
          EKIN=EKIN+G2(IG)*DMAT(IG)
        ENDDO
        EKIN=0.5D0*EKIN
!       == BUCKET POTENTIAL ============================================
        IF(TBUCKET) THEN
          EBUCKET=0.D0
          DO IG=1,NGL
            EBUCKET=EBUCKET+BUCKET(IG)*DMAT(IG)
          ENDDO
          EKIN=EKIN+EBUCKET
        END IF
!       ==  PARALLELIZE AND CLOSE DOWN =================================
        EKIN=EKIN*GWEIGHT
        CALL MPE$COMBINE('+',EKIN)
        DEALLOCATE(DMAT)
        DEALLOCATE(G2)
      END IF
!
      RETURN
      END
!
!     .....................................................RHOOFR.......
      SUBROUTINE WAVES_DENSITY(NGL,NRL,NDIM,NB,NBH,F,PSIOFG,RHO)
!     ******************************************************************
!     **                                                              **
!     **  THE  ELECTRON DENSITY RHOE IN REAL SPACE                    **
!     **                                                              **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NGL         ! MAX # PLANE WAVES
      INTEGER(4),INTENT(IN)  :: NRL         ! # R-SPACE POINTS
      INTEGER(4),INTENT(IN)  :: NDIM        ! DIMENSION OF THE WAVE FUNCTION
      INTEGER(4),INTENT(IN)  :: NB          ! # BANDS
      INTEGER(4),INTENT(IN)  :: NBH
      REAL(8)   ,INTENT(IN)  :: F(NB)       ! OCCUPATIONS
      COMPLEX(8),INTENT(IN)  :: PSIOFG(NGL,NDIM,NBH)
      REAL(8)   ,INTENT(OUT) :: RHO(NRL,NDIM**2) ! DENSITY IN R-SPACE
      COMPLEX(8),ALLOCATABLE :: PSIOFR(:,:,:)
      COMPLEX(8)             :: EI2KR(NRL)   !(NRL) SQUARED BLOCH PHASE FACTOR 
      COMPLEX(8)             :: PSI1(NRL)
      COMPLEX(8)             :: CSVAR
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IB,IBH,IDIM1,IDIM2,IR,IC
      REAL(8)                :: F1,F2
      REAL(8)                :: RE,IM
      REAL(8)                :: SVAR1,SVAR2
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
!     ******************************************************************
!RELEASED 8.OCT.99 
!
!     ==================================================================
!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT==
!     ==================================================================
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
!
!     ==================================================================
!     ==  FFT                                                         ==
!     ==================================================================
      ALLOCATE(PSIOFR(NRL,NDIM,NBH))
      CALL PLANEWAVE$FFT('GTOR',NBH*NDIM,NGL,PSIOFG,NRL,PSIOFR)
!
!     ==================================================================
!     ==  CALCULATE CHARGE DENSITY                                    ==
!     ==================================================================
!
!     ==================================================================
!     ==  CALCULATE DENSITY FOR SUPER WAVE FUNCTIONS                  ==
!     ==================================================================
      IF(TINV) THEN
        IF(NDIM.EQ.2) THEN
          CALL ERROR$MSG('SPINOR WAVE FUNCTIONS ARE NOT COMPATIBLE')
          CALL ERROR$MSG('WITH SUPER WAVE FUNCTIONS')
          CALL ERROR$STOP('WAVES_DENSITY')
        END IF
        CALL PLANEWAVE$GETC8A('EIKR',NRL,EI2KR)
        DO IR=1,NRL
          EI2KR(IR)=EI2KR(IR)**2
        ENDDO
        RHO(:,:)=0.D0
        DO IBH=1,NBH
          F1=F(2*IBH-1)
          F2=0.D0
          IF(2*IBH.LE.NB) F2=F(2*IBH)
          SVAR1=0.5D0*(F1+F2)
          SVAR2=0.5D0*(F1-F2)
!         == NOTE THAT OOCUPATIONS CAN BE NEGATIVE FROM K-INTEGRATION
          IF(SVAR1.EQ.0.D0) THEN
            IF(SVAR2.EQ.0) CYCLE
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
            RHO(IR,1)=RHO(IR,1)+REAL(PSIOFR(IR,1,IBH)*PSI1(IR))
          ENDDO
        ENDDO
!
!     ==================================================================
!     ==  CALCULATE DENSITY FOR GENERAL WAVE FUNCTIONS                ==
!     ==================================================================
      ELSE
        IF(NDIM.EQ.1) THEN
!       == SCALAR WAVE FUNCTIONS =======================================
          RHO(:,:)=0.D0
          DO IBH=1,NB
            F1=F(IBH)
            IF(F1.EQ.0.D0) CYCLE
            DO IR=1,NRL
              RE= REAL(PSIOFR(IR,1,IBH))
              IM=AIMAG(PSIOFR(IR,1,IBH))
              RHO(IR,1)=RHO(IR,1)+F1*(RE**2+IM**2)
            ENDDO
          ENDDO
        ELSE
!       == SPINOR WAVE FUNCTIONS =======================================
          RHO(:,:)=0.D0
          DO IBH=1,NBH
            F1=F(IBH)
            IF(F1.EQ.0.D0) CYCLE
!           == REAL(RHO11), REAL(RHO22), RE(RHO12), IM(RHO12) ==========
            DO IR=1,NRL
              RHO(IR,1)=RHO(IR,1) &
    &                 +F1*REAL(PSIOFR(IR,1,IBH)*CONJG(PSIOFR(IR,1,IBH)))
              RHO(IR,4)=RHO(IR,4) &
    &                 +F1*REAL(PSIOFR(IR,2,IBH)*CONJG(PSIOFR(IR,2,IBH)))
              CSVAR=PSIOFR(IR,1,IBH)*CONJG(PSIOFR(IR,2,IBH))
              RHO(IR,2)=RHO(IR,2)+F1*REAL(CSVAR)
              RHO(IR,3)=RHO(IR,3)+F1*AIMAG(CSVAR)
            ENDDO
          ENDDO
!         == TRANSFORM TO NT,NX,NY,NZ ==================================
          DO IR=1,NRL
            SVAR1=RHO(IR,1)
            SVAR2=RHO(IR,4)
            RHO(IR,1)=SVAR1+SVAR2
            RHO(IR,4)=SVAR1-SVAR2
            RHO(IR,2)= 2.D0*RHO(IR,2)
            RHO(IR,3)=-2.D0*RHO(IR,3)
          ENDDO
        ENDIF
      END IF
      DEALLOCATE(PSIOFR)
      RETURN
      END
!
!     .....................................................NLSMX........
      SUBROUTINE WAVES_DEDPROJ(NDIM,NBH,NB,LNX,LOX,LMNX,OCC &
     &                       ,PROJ,DH,DO,EPS,DEDPROJ)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE DERIVATIVE OF THE ONE-CENTER                 **
!     **  ENERGIES WITH RESPECT TO THE PROJECTOR FUNCTIONS            **
!     **    DEDPROJ = DE/(DPROJ)                                      **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NDIM      ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN)   :: NB        ! #(BANDS)
      INTEGER(4),INTENT(IN)   :: NBH       ! #(WAVE FUNCTIONS)
      INTEGER(4),INTENT(IN)   :: LMNX      ! #(PROJECTORS ON THIS SITE)
      REAL(8)   ,INTENT(IN)   :: OCC(NB)   ! OCCUPATIONS
      COMPLEX(8),INTENT(IN)   :: PROJ(NDIM,NBH,LMNX) ! <P|PSI>
      REAL(8)   ,INTENT(IN)   :: DH(LMNX,LMNX,NDIM**2) ! DE/DD
      INTEGER(4),INTENT(IN)   :: LNX
      INTEGER(4),INTENT(IN)   :: LOX(LNX)
      REAL(8)   ,INTENT(IN)   :: DO(LNX,LNX) ! DO/DD
      COMPLEX(8),INTENT(IN)   :: EPS(NB,NB)
      COMPLEX(8),INTENT(OUT)  :: DEDPROJ(NDIM,NBH,LMNX)    ! DE/D<P|
      LOGICAL(4)              :: TINV
      COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
      COMPLEX(8),ALLOCATABLE  :: OPROJ(:,:,:)
      COMPLEX(8),ALLOCATABLE  :: LAMBDA(:,:)
      REAL(8)                 :: F1,F2
      INTEGER(4)              :: IB,IB1,IB2,LMN,IDIM
      COMPLEX(8)              :: CSVAR
!     ******************************************************************
!
!     ==================================================================
!     ==  CHECK IF SUPERWAVEFUNCTIONS ARE USED AND IF #(BANDS) CORRECT==
!     ==================================================================
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
!     ==================================================================
!     ==  HPROJ = DH<P|PSI>                                           ==
!     ==================================================================
      CALL WAVES_HPROJ(NDIM,NBH,LMNX,dh,PROJ,DEDPROJ)
!
!     ==================================================================
!     ==  HPROJ = DH<P|PSI>*F                                         ==
!     ==================================================================
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
!     ==================================================================
!     ==  ADD -DO<P|PSI>LAMBDA*(FI+FJ)/2                              ==
!     ==================================================================
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
      RETURN
      END
!
!     .....................................................NLSMX........
      SUBROUTINE WAVES_DEDPRO(TINV,NGL,NDIM,NBH,PSI,LMNX,DEDPROJ,DEDPRO)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE FUNCTION |DEDPRO_I>=DE/D<PRO_I|              **
!     **  USED TO CALCULATE FORCES AND STRESSES AS                    **
!     **    DE = SUM_I <DPRO_I|DEDPRO_I> + <DEDPRO_I|DPRO_I>          **
!     **  WHERE                                                       **
!     **    |DEDPRO_I>=|DEDPRO_II> + SUM_N |PSI_N> F_N <PSI_N|PRO_I>  **
!     **                   +(DH-DO EPS)                               **
!     **                                                              **
!     **   REMARKS:                                                   **
!     **    - CALL FOR EACH ATOM INDIVIDUALLY                         **
!     **    - INITIALIZE DEDPRO=(0.D0,0.D0) BEFORE CALLING            **
!     **    - FOR LSD INITIALIZE ONCE AND CALL FOR SPIN UP AND DOWN   **
!     **    - DEDPRO SHOULD BE A REAL FUNCTION IN REAL SPACE          **
!     **       BUT THAT THE SYMMETRY IS NOT ENFORCES YET              **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)   :: TINV
      INTEGER(4),INTENT(IN)   :: NGL       ! (#(G-VECTORS)
      INTEGER(4),INTENT(IN)   :: NDIM      ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN)   :: NBH       ! #(WAVE FUNCTIONS)
      COMPLEX(8),INTENT(IN)   :: PSI(NGL,NDIM,NBH)   ! |PSI>
      INTEGER(4),INTENT(IN)   :: LMNX      ! #(PROJECTORS ON THIS SITE)
      COMPLEX(8),INTENT(IN)   :: DEDPROJ(NDIM,NBH,LMNX) ! <P|PSI>
      COMPLEX(8),INTENT(OUT)  :: DEDPRO(NGL,LMNX)    ! DE/D<P|
      INTEGER(4)              :: IB,IDIM,LMN,IG
      COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
      COMPLEX(8),ALLOCATABLE  :: PSIM(:)
      LOGICAL(4),PARAMETER    :: TESSL=.TRUE.
      COMPLEX(8)              :: DEDPROJ1(NDIM,NBH,LMNX) ! CONJG(<P|PSI>)
      COMPLEX(8)              :: CSVAR
!     ******************************************************************
!
!     ==================================================================
!     ==  SUPERPOSE WAVE FUNCTION TO OBTAIN DE/D<P|                   ==
!     ==================================================================
      DEDPROJ1=CONJG(DEDPROJ)
      CALL LIB$MATMULC8(NGL,NDIM*NBH,LMNX,PSI,DEDPROJ1,DEDPRO)
!     IF(TESSL) THEN
!        DEDPROJ1=CONJG(DEDPROJ)
!        CALL ZGEMUL(PSI,NGL,'N',DEDPROJ1,NDIM*NBH,'N',DEDPRO,NGL &
!    &              ,NGL,NDIM*NBH,LMNX)
!     ELSE
!       DEDPRO(:,:)=(0.D0,0.D0)
!       DO LMN=1,LMNX
!         DO IB=1,NBH
!           DO IDIM=1,NDIM
!             CSVAR=CONJG(DEDPROJ(IDIM,IB,LMN))
!             DO IG=1,NGL
!               DEDPRO(IG,LMN)=DEDPRO(IG,LMN)+PSI(IG,IDIM,IB)*CSVAR
!             ENDDO
!           ENDDO
!         ENDDO
!       ENDDO!

!     END IF
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
!     ..................................................................
      SUBROUTINE WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NB,NPRO,PSI,PROPSI)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE PROJECTIONS                                       **
!     **                                                              **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1998)***
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : MAP_TYPE,GSET_TYPE
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
      REAL(8)        ,ALLOCATABLE:: GVEC(:,:)  !(3,NGL)
      COMPLEX(8)     ,ALLOCATABLE:: PRO(:,:)   !(NGL,LMNXX)
      COMPLEX(8)     ,ALLOCATABLE:: PROPSI1(:) !(LMNX,NDIM,NB)
      COMPLEX(8)     ,ALLOCATABLE:: EIGR(:)    !(NGL)
      INTEGER(4)     ,ALLOCATABLE:: LOX(:)     !(LNXX)
      INTEGER(4)                 :: IPRO,IBPRO,IAT,LMN,IB,IDIM,II
      INTEGER(4)                 :: LMNXX,LNXX
      INTEGER(4)                 :: LNX,LMNX,ISP
!     ******************************************************************
                                CALL TIMING$CLOCKON('WAVES_PROJECTIONS')
!
!     ==================================================================
!     ==  CONSISTENCY CHECKS                                          ==
!     ==================================================================
      IF(NPRO.NE.MAP%NPRO) THEN
        CALL ERROR$MSG('INCONSISTENT NPRO')
        CALL ERROR$STOP('WAVES_PROJECTIONS')
      END IF
      IF(GSET%NGL.NE.NGL) THEN
        CALL ERROR$MSG('INCONSISTENT NGL')
        CALL ERROR$STOP('WAVES_PROJECTIONS')
      END IF
!
!     ==================================================================
!     ==  COLLECT G-VECTORS                                           ==
!     ==================================================================
      ALLOCATE(GVEC(3,NGL))
      CALL PLANEWAVE$GETR8A('GVEC',3*GSET%NGL,GVEC)
!
!     ==================================================================
!     ==  <PRO|PSI>                                                   ==
!     ==================================================================
      LMNXX=MAXVAL(MAP%LMNX)
      LNXX=MAXVAL(MAP%LNX)
      ALLOCATE(PRO(NGL,LMNXX))
      ALLOCATE(EIGR(NGL))
      ALLOCATE(LOX(LNXX))
      ALLOCATE(PROPSI1(LMNXX*NDIM*NB))
      IPRO=1
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
        LNX=MAP%LNX(ISP)
        LMNX=MAP%LMNX(ISP)
        LOX=MAP%LOX(:,ISP)
        IBPRO=1+SUM(MAP%LNX(1:ISP-1))
        CALL PLANEWAVE$STRUCTUREFACTOR(R(1,IAT),NGL,EIGR)
        CALL WAVES_EXPANDPRO(LNX,LOX,LMNX,NGL,GVEC &
     &                      ,GSET%PRO(:,IBPRO:ibpro+lnx-1),MAP%LMX,GSET%YLM,EIGR,PRO)
        CALL PLANEWAVE$SCALARPRODUCT(' ',NGL,1,LMNX,PRO,NDIM*NB,PSI &
     &                              ,PROPSI1)
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
       DEALLOCATE(LOX)
      DEALLOCATE(PRO)
      DEALLOCATE(PROPSI1)
      DEALLOCATE(EIGR)
      DEALLOCATE(GVEC)
      CALL MPE$COMBINE('+',PROPSI)
                               CALL TIMING$CLOCKOFF('WAVES_PROJECTIONS')
      RETURN
      END
!
!     .....................................................NLSMX........
      SUBROUTINE WAVES_ADDPRO(MAP,GSET,NAT,R,NGL,NDIM,NB,NPRO,PSI,PROPSI)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE PROJECTIONS                                       **
!     **                                                              **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1998)***
      USE WAVES_MODULE, ONLY : MAP_TYPE,GSET_TYPE
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
      REAL(8)        ,ALLOCATABLE:: GVEC(:,:)  !(3,NGL)
      COMPLEX(8)     ,ALLOCATABLE:: PRO(:,:)   !(NGL,NDIM,LBLOCK)
      COMPLEX(8)     ,ALLOCATABLE:: PROPSI1(:) !(LMNX*NDIM*NB)
      COMPLEX(8)     ,ALLOCATABLE:: EIGR(:)    !(NGL)
      INTEGER(4)     ,ALLOCATABLE:: LOX(:)     !(LNXX)
      COMPLEX(8)                 :: CSVAR
      REAL(8)                    :: SVAR
      LOGICAL(4)    ,PARAMETER   :: TESSL=.TRUE.
      INTEGER(4)                 :: IPRO,IBPRO,IAT,LMN,IB,IDIM,IG,II
      INTEGER(4)                 :: LMNXX,LNXX
      INTEGER(4)                 :: LNX,LMNX,ISP
!     ******************************************************************
                                 CALL TIMING$CLOCKON('WAVES_ADDPRO')
!
!     ==================================================================
!     ==  CONSISTENCY CHECKS                                          ==
!     ==================================================================
      IF(NPRO.NE.MAP%NPRO) THEN
        CALL ERROR$MSG('INCONSISTENT NPRO')
        CALL ERROR$STOP('WAVES_ADDPRO')
      END IF
      IF(GSET%NGL.NE.NGL) THEN
        CALL ERROR$MSG('INCONSISTENT NGL')
        CALL ERROR$STOP('WAVES_ADDPRO')
      END IF
!
!     ==================================================================
!     ==  COLLECT G-VECTORS                                           ==
!     ==================================================================
      ALLOCATE(GVEC(3,NGL))
      CALL PLANEWAVE$GETR8A('GVEC',3*GSET%NGL,GVEC)
!
!     ==================================================================
!     ==  COLLECT G-VECTORS                                           ==
!     ==================================================================
      LMNXX=MAXVAL(MAP%LMNX)
      LNXX=MAXVAL(MAP%LNX)
      ALLOCATE(PRO(NGL,LMNXX))
      ALLOCATE(PROPSI1(LMNXX*NDIM*NB))
      ALLOCATE(EIGR(NGL))
      ALLOCATE(LOX(LNXX))
      IPRO=1
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
        LNX=MAP%LNX(ISP)
        LMNX=MAP%LMNX(ISP)
        LOX=MAP%LOX(:,ISP)
        IBPRO=1+SUM(MAP%LNX(1:ISP-1))
        CALL PLANEWAVE$STRUCTUREFACTOR(R(1,IAT),NGL,EIGR)
        CALL WAVES_EXPANDPRO(LNX,LOX,LMNX,NGL,GVEC &
     &                      ,GSET%PRO(:,IBPRO:ibpro+lnx-1),MAP%LMX,GSET%YLM,EIGR,PRO)
!
!       ===============================================================
!       ==  |PSI>=|PSI>+|PRO>PROPSI                                  ==
!       ===============================================================
        II=0
        DO IB=1,NB
          DO IDIM=1,NDIM
            DO LMN=1,LMNX
              II=II+1    ! II=(LMN,IDIM,IB)
              PROPSI1(II)=PROPSI(IDIM,IB,IPRO-1+LMN)
            ENDDO
          ENDDO
        ENDDO
        CALL LIB$ADDPRODUCTC8(.FALSE.,NGL,LMNX,NDIM*NB,PRO,PROPSI1,PSI)
!       IF(TESSL) THEN
!         CALL ZGEMM('N','T',NGL,NDIM*NB,LMNX,(1.D0,0.D0),PRO,NGL &
!   &               ,PROPSI1,NDIM*NB,(1.D0,0.D0),PSI,NGL)
!         CALL ZGEMM('N','N',NGL,NDIM*NB,LMNX,(1.D0,0.D0),PRO,NGL &
!   &               ,PROPSI1,LMNX,(1.D0,0.D0),PSI,NGL)
!       ELSE
!         II=0
!         DO IB=1,NB
!           DO IDIM=1,NDIM
!             DO LMN=1,LMNX
!               II=II+1  ! II=(LMN,IDIM,IB)
!               CSVAR=PROPSI1(II)
!               DO IG=1,NGL
!                 PSI(IG,IDIM,IB)=PSI(IG,IDIM,IB)+PRO(IG,LMN)*CSVAR     
!              ENDDO
!             ENDDO
!           ENDDO
!         ENDDO
!       END IF
        IPRO=IPRO+LMNX
      ENDDO
      DEALLOCATE(LOX)
      DEALLOCATE(PRO)
      DEALLOCATE(GVEC)
      DEALLOCATE(EIGR)
      DEALLOCATE(PROPSI1)
                                 CALL TIMING$CLOCKOFF('WAVES_ADDPRO')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_EXPANDPRO(LNX,LOX,LMNX,NGL,GVEC &
     &                          ,BAREPRO,LMX,YLM,EIGR,PRO)
!     ******************************************************************
!     **                                                              **
!     **  TAKES THE BARE PROJECTOR FUNCTIONS AND CALCULATES FULL      **
!     **  PROJECTOR FUNCTIONS WITH STRUCTURE FACTOR AND I**L          **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1998)***
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
!     ******************************************************************
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
!     ..................................................................
      SUBROUTINE WAVES_EXPANDPRO_OLD(IPRO1,NPRO &
     &                          ,NSP,LNX,LNXX,LOX,NAT,ISPECIES,NGL,GVEC &
     &                          ,NBAREPRO,BAREPRO,LMX,YLM,EIGR,PRO)
!     ******************************************************************
!     **                                                              **
!     **  TAKES THE BARE PROJECTOR FUNCTIONS AND CALCULATES FULL      **
!     **  PROJECTOR FUNCTIONS WITH STRUCTURE FACTOR AND I**L          **
!     **                                                              **
!     **  FOR NDIM.NE.3 FIRST DERIVATIVES ARE CALCULATED              **
!     **  FOR NDIM.NE.6 SECOND DERIVATIVES ARE CALCULATED             **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1998)***
      IMPLICIT NONE
      INTEGER(4)    ,INTENT(IN) :: IPRO1    ! FIRST PROJECTOR
      INTEGER(4)    ,INTENT(IN) :: NPRO     ! X#(NUMBER OF PROJECTORS)
      INTEGER(4)    ,INTENT(IN) :: NSP      ! #(SPECIES)
      INTEGER(4)    ,INTENT(IN) :: LNX(NSP) ! #(PROJECTORS PER ATOM)
      INTEGER(4)    ,INTENT(IN) :: LNXX     ! #(PROJECTORS PER ATOM)
      INTEGER(4)    ,INTENT(IN) :: LOX(LNXX,NSP) ! ANGULAR MOMENTA
      INTEGER(4)    ,INTENT(IN) :: NAT           ! #(ATOMS)
      INTEGER(4)    ,INTENT(IN) :: ISPECIES(NAT) ! SPECIES OF THE ATOM
      INTEGER(4)    ,INTENT(IN) :: NGL           ! #(G-VECTORS(LOCAL))
      REAL(8)       ,INTENT(IN) :: GVEC(3,NGL)   ! G-VECTORS
      INTEGER(4)    ,INTENT(IN) :: NBAREPRO      ! X#(NUMBER OF PROJECTORS)
      REAL(8)       ,INTENT(IN) :: BAREPRO(NGL,NBAREPRO)  ! PROJECTOR FUNCTIONS
      INTEGER(4)    ,INTENT(IN) :: LMX           ! X#(ANGULAR MOMENTA IN YLM))
      REAL(8)       ,INTENT(IN) :: YLM(NGL,LMX)  ! CI**(-L)*SPHERICAL HARMONICS
      COMPLEX(8)    ,INTENT(IN) :: EIGR(NGL,NAT)     ! STRUCTURE FACTORS
      COMPLEX(8)    ,INTENT(OUT):: PRO(NGL,NPRO) !PROJECTOR FUNCTION
      COMPLEX(8)    ,PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)                :: CSVAR
      INTEGER(4)                :: NDER
      INTEGER(4)                :: ICOUNT,IPRO
      INTEGER(4)                :: IAT,ISP,LMN,LN,IM,L,IG,LM,IBAREPRO
      LOGICAL(4)                :: TGRAD,TSTRESS
      REAL(8)                   :: SVAR
!     ******************************************************************
!
!     ==================================================================
!     ==  START THE LOOP                                              ==
!     ==================================================================
      ICOUNT=0
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        IBAREPRO=SUM(LNX(1:ISP-1))
        DO LN=1,LNX(ISP)
          IBAREPRO=IBAREPRO+1
          L=LOX(LN,ISP)
          LM=L**2
          DO IM=1,2*L+1
            LM=LM+1
            IPRO=IPRO+1
            IF(IPRO.LT.IPRO1) CYCLE
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.NPRO) RETURN
!           
!           ==========================================================
!           ==  COMPOSE PROJECTOR FUNCTIONS                         ==
!           ==========================================================
!           ==  MULTIPLY PROJECTOR WITH STRUCTURE FACTOR ===========
            CSVAR=(-CI)**L
            DO IG=1,NGL
              SVAR=YLM(IG,LM)*BAREPRO(IG,IBAREPRO)
              PRO(IG,ICOUNT)=SVAR*CSVAR*EIGR(IG,IAT)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_PROFORCE(LNX,LMNX,LOX,NGL,GWEIGHT,GVEC,GIJ &
     &           ,BAREPRO,DBAREPRO,LMX,YLM,SYLM,EIGR,DEDPRO,FORCE,TSTRESS,STRESS)
!     ******************************************************************
!     **                                                              **
!     **  TAKES THE BARE PROJECTOR FUNCTIONS AND CALCULATES FULL      **
!     **  PROJECTOR FUNCTIONS WITH STRUCTURE FACTOR AND I**L          **
!     **                                                              **
!     **  FOR NDIM.NE.3 FIRST DERIVATIVES ARE CALCULATED              **
!     **  FOR NDIM.NE.6 SECOND DERIVATIVES ARE CALCULATED             **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1998)***
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
      REAL(8)   ,PARAMETER  :: TINY=1.D-300
      COMPLEX(8)            :: CSVAR,CSVAR1,CSVAR2,CIL
      COMPLEX(8)            :: CWORK1(NGL)
      COMPLEX(8)            :: CWORK2(NGL)
      REAL(8)               :: SVAR,SVAR1,SVAR2
      REAL(8)               :: F(3)
      REAL(8)               :: S(6),SDIAG
      INTEGER(4)            :: LN,LMN,LM,L,M
      INTEGER(4)            :: IG,I,J,IJ
!     ******************************************************************
      F(:)=0.D0
      STRESS=0.D0
      S(:)=0.D0
      SDIAG=0.D0
!   
!     ==================================================================
!     ==  START THE LOOP                                              ==
!     ==================================================================
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
!         ==========================================================
!         ==  COMPOSE PROJECTOR FUNCTIONS                         ==
!         ==========================================================
!         ==  MULTIPLY PROJECTOR WITH STRUCTURE FACTOR ===========
          IF(TSTRESS) THEN
            DO IG=1,NGL
!             == THE USE OF CONJG IN THE FOLLOWING LINE IS UNCLEAR
!             == PREVIOUSLY CONJG WAS USED WHENEVER TINV=.TRUE.
              CSVAR =CONJG(DEDPRO(IG,LMN))
              CSVAR1=CSVAR*CWORK1(IG)
              CSVAR2=CSVAR*CWORK2(IG)
              SVAR  =REAL(CI*CSVAR1)*YLM(IG,LM)
              DO I=1,3
                F(I)=F(I)+SVAR*GVEC(I,IG)
              ENDDO
              SVAR1=REAL(CSVAR1)
              SVAR2=REAL(CSVAR2)*YLM(IG,LM)
              DO IJ=1,6
                S(IJ)=S(IJ)-SVAR1*SYLM(IG,LM,IJ)-SVAR2*GIJ(IJ,IG)
              ENDDO
              SDIAG=SDIAG-SVAR1*YLM(IG,LM)
            ENDDO
          ELSE
            DO IG=1,NGL
              SVAR=REAL(CWORK1(IG)*CONJG(DEDPRO(IG,LMN)))*YLM(IG,LM)
              DO I=1,3
                F(I)=F(I)+SVAR*GVEC(I,IG)
              ENDDO
            ENDDO
          END IF
        ENDDO
      ENDDO
!
!     ==================================================================
!     == <P|DE/D<P| + DE/D|P> |P>                                     ==
!     ==================================================================
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
!     ..................................................................
      SUBROUTINE WAVES_UPDATEGSET()
!     ******************************************************************
!     **                                                              **
!     **  THIS ROUTINE MUST BE CALLED INITIALLY AND AFTER EACH        **
!     **  CHANGE OF THE CELL-SHAPE.                                   **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
      USE WAVES_MODULE
      IMPLICIT NONE
      REAL(8)        ,ALLOCATABLE   :: G2(:)
      REAL(8)        ,ALLOCATABLE   :: GVEC(:,:)
      REAL(8)        ,ALLOCATABLE   :: YLM(:)
      REAL(8)        ,ALLOCATABLE   :: R(:,:)
      INTEGER(4)                    :: IG,ISP,IND,LN,IKPT,ISPIN
      INTEGER(4)                    :: NBAREPRO,NGL
      REAL(8)                       :: RBAS(3,3),GBAS(3,3),CELLVOL
      INTEGER(4)                    :: LMX
      LOGICAL(4)                    :: TSTRESS
      REAL(8)       ,ALLOCATABLE    :: WORK(:,:)
      REAL(8)                       :: MING2,MING,DG,ABSG
!     ******************************************************************
                            CALL TRACE$PUSH('WAVES_UPDATEGSET')
!
!     ==================================================================
!     ==  UPDATE RADIAL PROJECTOR FUNCTIONS                           ==
!     ==================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      CALL CELL$GETL4('ON',TSTRESS)
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL POTENTIAL$SETR8A('RBAS',9,RBAS)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      LMX=MAP%LMX
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL=GSET%NGL
        NBAREPRO=MAP%NBAREPRO
!
!       ================================================================
!       == UPDATE PLANEWAVE OBJECT                                    ==
!       ================================================================
        CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
!
!       ================================================================
!       ==  EVALUATE PROJECTOR FUNCTIONS IN G-SPACE                   ==
!       ================================================================
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
        ENDDO
!
!       ================================================================
!       ==  UPDATE SPHERICAL HARMONICS                                ==
!       ================================================================
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
!       == NOW THE STRAINED SPHERICAL HARMONICS ==========================
        CALL WAVES_STRAINEDYLM(NGL,LMX,GVEC,GSET%YLM,GSET%SYLM)
        DEALLOCATE(GVEC)
!
!       ================================================================
!       ==  EVALUATE BUCKET POTENTIAL                                 ==
!       ================================================================
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
!     ================================================================
!     ==  EVALUATE WAVE FUNCTION MASS                               ==
!     ================================================================
      CALL WAVES_PSIMASS
                               CALL TRACE$POP
      RETURN
      END
!
!     .....................................................NLSMX........
      SUBROUTINE WAVES_STRAINEDYLM(NGL,LMX,GVEC,YLM,SYLM)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES STRAINED SPHERICAL HARMONICS ON THE G-SPACE GRID  **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
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
      REAL(8)        ,PARAMETER     :: TINY=1.D-300
!     ******************************************************************
!
!     ==================================================================
!     ==  PREPARE GI*GJ/G2                                            ==
!     ==================================================================
      DO IG=1,NGL
        SVAR=1.D0/(GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2+TINY)
        IJ=0
        DO I=1,3
          DO J=I,3
            IJ=IJ+1
            GIJ(IG,IJ)=SVAR*GVEC(I,IG)*GVEC(J,IG)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  NOW STRAIN THE SPHERICAL HARMONICS                          ==
!     ==================================================================
      SYLM=0.D0
      DO LM1=1,LMX
!
!       ================================================================
!       ==                                                            ==
!       ================================================================
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
!       ================================================================
!       ==                                                            ==
!       ================================================================
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
!       ================================================================
!       ==  NOW THE TERM -L*GI*GJ/G**2                                ==
!       ================================================================
        L=INT(SQRT(REAL(LM1-1)+0.1))
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
!     ..................................................................
      SUBROUTINE WAVES$PROPAGATE()
!     ******************************************************************
!     **                                                              **
!     **  PROPAGATES WAVE FUNCTIONS                                   **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    PARAMETERS DELT,EMASS,EMASSCG2,ANNEE AND TSTOP MUST BE SET**
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)            :: IKPT,ISPIN,IG,IDIM,IB
      INTEGER(4)            :: NGL
      INTEGER(4)            :: NBH
      REAL(8)   ,ALLOCATABLE:: ARR1(:)
      REAL(8)   ,ALLOCATABLE:: ARR2(:)
      REAL(8)   ,ALLOCATABLE:: ARR3(:)
      REAL(8)   ,ALLOCATABLE:: GVEC(:,:)
      REAL(8)               :: SVAR1,SVAR2,SVAR3
!     ******************************************************************
                            CALL TRACE$PUSH('WAVES$PROPAGATE')
!
!     ==================================================================
!     ==  STOP WAVE FUNCTIONS                                         ==
!     ==================================================================
      IF(TSTOP) THEN
        DO IKPT=1,NKPT
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            THIS%PSIM(:,:,:)=THIS%PSI0(:,:,:)
          ENDDO
        ENDDO
        waveekin1=0.d0
        TSTOP=.FALSE.
      END IF
!
!     ==================================================================
!     ==  PROPAGATE WAVE FUNCTIONS (PUT PSI(+) INTP PSIM)             ==
!     ==================================================================
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL=GSET%NGL
!       == EVALUATE FACTORS ============================================
        ALLOCATE(ARR1(NGL))
        ALLOCATE(ARR2(NGL))
        ALLOCATE(ARR3(NGL))
        IF(ASSOCIATED(GSET%DMPSI)) THEN
          DO IG=1,NGL
!PB         SVAR1=1.D0+ANNEE+GSET%DMPSI(IG)
            SVAR1=1.D0+ANNEE
            ARR1(IG)=2.D0/SVAR1
            ARR2(IG)=1.D0-ARR1(IG)
            ARR3(IG)=-DELT**2/GSET%MPSI(IG)/SVAR1
          ENDDO
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
!       ==  NOW PROPAGATE ==============================================
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
          DO IB=1,NBH
            DO IDIM=1,NDIM
              DO IG=1,NGL
                THIS%PSIM(IG,IDIM,IB)=ARR1(IG)*THIS%PSI0(IG,IDIM,IB) &
       &                             +ARR2(IG)*THIS%PSIM(IG,IDIM,IB) &
       &                             +ARR3(IG)*THIS%HPSI(IG,IDIM,IB) 
              ENDDO
            ENDDO
          ENDDO
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
!     ..................................................................
      SUBROUTINE WAVES_PSIMASS
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE WAVE FUNCTION MASS AND ITS CHANGES           **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)   ,ALLOCATABLE :: GVEC(:,:)
      REAL(8)   ,ALLOCATABLE :: DG2(:)
      REAL(8)                :: FRICMAT(3,3)
      REAL(8)                :: G2
      REAL(8)                :: FAC1
      INTEGER(4)             :: IKPT,IG
      INTEGER(4)             :: NGL
      LOGICAL(4)             :: TSTRESS
!     ******************************************************************
      CALL CELL$GETL4('ON',TSTRESS)
      IF(TSTRESS) THEN
        CALL CELL$GETR8A('FRICMAT',9,FRICMAT) !\ALPHADOT*DELTA/2
      ELSE
        FRICMAT=0.D0
      END IF
      FAC1=EMASS*2.D0*EMASSCG2
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL=GSET%NGL
        IF(TSTRESS) THEN
          ALLOCATE(GVEC(3,NGL))
          CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
          IF(.NOT.ASSOCIATED(GSET%MPSI)) ALLOCATE(GSET%MPSI(NGL))
          IF(.NOT.ASSOCIATED(GSET%DMPSI))ALLOCATE(GSET%DMPSI(NGL))
          ALLOCATE(DG2(NGL))
          DO IG=1,NGL
            DG2(IG)=FRICMAT(1,1)     *GVEC(1,IG)*GVEC(1,IG) &
    &              +FRICMAT(1,2)*2.D0*GVEC(1,IG)*GVEC(2,IG) &
    &              +FRICMAT(1,3)*2.D0*GVEC(1,IG)*GVEC(3,IG) &
    &              +FRICMAT(2,2)     *GVEC(2,IG)*GVEC(2,IG) &
    &              +FRICMAT(2,3)*2.D0*GVEC(2,IG)*GVEC(3,IG) &
    &              +FRICMAT(3,3)     *GVEC(3,IG)*GVEC(3,IG)
          ENDDO
!         == MASS ======================================================
          DO IG=1,NGL
            G2=GVEC(1,IG)**2+GVEC(2,IG)**2+GVEC(3,IG)**2
            GSET%MPSI(IG) =EMASS*(1.D0+EMASSCG2*G2)
            GSET%DMPSI(IG)=EMASS*EMASSCG2*2.D0*DG2(IG)
          END DO
!         == ADD MASS FOR BUCKET POTENTIAL =============================
          IF(TBUCKET) THEN
            DO IG=1,NGL
              GSET%MPSI(IG) =GSET%MPSI(IG) +FAC1*GSET%BUCKET(IG)
              GSET%DMPSI(IG)=GSET%DMPSI(IG)+FAC1*GSET%DBUCKET(IG)*DG2(IG)
            ENDDO
          END IF
          DEALLOCATE(DG2)
          DEALLOCATE(GVEC)
        ELSE 
          IF(ASSOCIATED(GSET%DMPSI))DEALLOCATE(GSET%DMPSI)
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
!     ..................................................................
      SUBROUTINE WAVES_WAVEKINETIC(EKIN)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES FICTITIOUS KINETIC ENERGY OF THE WAVE FUNCTIONS  **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    PARAMETERS DELT,EMASS,EMASSCG2 MUST BE SET                **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT):: EKIN
      INTEGER(4)            :: IKPT,ISPIN,IG,IDIM,IB
      INTEGER(4)            :: NGL
      INTEGER(4)            :: NBX
!REAL(8)   ,ALLOCATABLE:: MARR(:)
      COMPLEX(8),ALLOCATABLE:: TPSIP(:)
      COMPLEX(8),ALLOCATABLE:: TPSIM(:)
      REAL(8)               :: EKIN1,SUM
      REAL(8)   ,ALLOCATABLE:: OCC(:,:,:)
      REAL(8)   ,ALLOCATABLE:: EIG(:,:,:)
      REAL(8)               :: RBAS(3,3),GBAS(3,3),CELLVOL
      REAL(8)               :: V1,V2
      REAL(8)               :: F1,F2
      REAL(8)               :: DOCC
      LOGICAL(4)            :: TINV
      COMPLEX(8)            :: CSVAR
      REAL(8)   ,PARAMETER  :: DSMALL=1.D-12
      INTEGER(4)            :: IB1,IB2
!     ******************************************************************
!
!     ==================================================================
!     == COLLECT UNIT CELL VOLUME                                     ==
!     ==================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL GBASS(RBAS,GBAS,CELLVOL)
!
!     ==================================================================
!     == COLLECT OCCUPATIONS                                          ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPT,NSPIN))
      ALLOCATE(EIG(NBX,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NBX*NKPT*NSPIN,OCC)
      EIG(:,:,:)=0.D0
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      EKIN=0.D0
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          TINV=GSET%TINV
!
!ALLOCATE(MARR(NGL))
!CALL PLANEWAVE$GETR8A('G2',NGL,MARR)
!DO IG=1,NGL
!MARR(IG)=EMASS*(1.D0+EMASSCG2*MARR(IG))/DELT**2
!ENDDO
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
!           ============================================================
!           ==  WAVE FUNCTION KINETIC ENERGY FOR GENERAL WAVE FUNCTIONS=
!           ==  AND FIRST PART FOR SUPER WAVE FUNCTIONS                =
!           ============================================================
            SUM=0.D0
            DO IDIM=1,NDIM
              DO IG=1,NGL
                CSVAR=THIS%PSI0(IG,IDIM,IB)-THIS%PSIM(IG,IDIM,IB)
                V1=REAL(CSVAR)
                V2=AIMAG(CSVAR)
!               SUM=SUM+MARR(IG)*(V1**2+V2**2)
                SUM=SUM+GSET%MPSI(IG)*(V1**2+V2**2)
              ENDDO
            ENDDO
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
!           ============================================================
!           ==  NOW CONTINUE WITH SUPERWAVE FUNCTIONS ONLY            ==
!           ============================================================
!           == PSI- ONLY REQUIRED IF THE OCCUPATIONS DIFFER ==========
            IF(F1-F2.EQ.0.D0) CYCLE
            SUM=0.D0
            DO IDIM=1,NDIM
              TPSIP(:)=THIS%PSI0(:,IDIM,IB)-THIS%PSIM(:,IDIM,IB)
              CALL PLANEWAVE$INVERTG(NGL,TPSIP,TPSIM)
              DO IG=1,NGL
!               SUM=SUM+MARR(IG)*(REAL(TPSIP(IG))*REAL(TPSIM(IG)) &
!    &                         +AIMAG(TPSIP(IG))*AIMAG(TPSIM(IG)))
                SUM=SUM+GSET%MPSI(IG)*(REAL(TPSIP(IG))*REAL(TPSIM(IG)) &
     &                               +AIMAG(TPSIP(IG))*AIMAG(TPSIM(IG)))
              ENDDO            
            ENDDO
            EIG(IB1,IKPT,ISPIN)=EIG(IB1,IKPT,ISPIN)+0.5D0*SUM              
            EIG(IB2,IKPT,ISPIN)=EIG(IB2,IKPT,ISPIN)-0.5D0*SUM              
            EKIN1=EKIN1+SUM*0.5D0*(F1-F2)
          ENDDO
          DEALLOCATE(TPSIP)
          DEALLOCATE(TPSIM)
!DEALLOCATE(MARR)
          EKIN=EKIN+EKIN1*CELLVOL/DELT**2
        ENDDO
      ENDDO
!
!     ==================================================================
!     == COMMUNICATE                                                  ==
!     ==================================================================
      CALL MPE$COMBINE('+',EKIN)
      CALL MPE$COMBINE('+',EIG)
!
!     ==================================================================
!     == SET CONTRIBUTION TO EIGENVALUES                              ==
!     ==================================================================
      CALL DYNOCC$SETR8A('M<PSIDOT|PSIDOT>',NBX*NKPT*NSPIN,EIG)
      DEALLOCATE(EIG)
      DEALLOCATE(OCC)
      RETURN
      END
