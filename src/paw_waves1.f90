!***********************************************************************
!***********************************************************************
!**                                                                   **
!**  WAVES OBJECT                                                     **
!**  THIS OBJECT IS CONTINUED IN PAW_WAVES2.F90                       **
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
!**  DEFINITION OF SPIN DIMENSIONS:                                   **
!**                           | NSPIN | NDIM | NDIMD |                **
!**        -------------------------------------------                **
!**        NON-SPIN-POLARIZED |   1   |   1  |   1   |                **
!**        SPIN POLARIZED     |   2   |   1  |   2   |                **
!**        NON-COLLINEAR      |   1   |   2  |   4   |                **
!**                                                                   **
!**     NSPIN IS THE NUMBER OF SLATER DETERMINANTS                    **
!**     NDIM IS THE NUMBER OF SPINOR DIMENSIONS OF THE WAVE FUNCTION  **
!**     NDIMD IS THE NUMBER OF DENSITY COMPONENTS. NDIMD=NSPIN*NDIM**2**
!**     TWO REPRESENTATIONS ARE USED:                                 **
!**           NSPIN=2,NDIM=1:  UP/DOWN AND TOTAL/SPIN                 **
!**           NSPIN=1,NDIM=2:  UPUP/UPDOWN/DOWNUP/DOWNDOWN AND        **
!**                            TOTAL/MX/MY/MZ                         **
!**                                                                   **
!***********************************************************************
!***********************************************************************
!
!......................................................WAVES............
MODULE WAVES_MODULE                                                
!***********************************************************************
!**                                                                   **
!**                                                                   **
!************************************************P.E. BLOECHL, (1995)***
USE LINKEDLIST_MODULE
! ----------------------------------------------------------------------
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
  REAL(8)   ,POINTER :: EIGVAL(:)       !(NB)    EIGENVALUES AFTER DIAG
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
REAL(8)     :: DELT=0.D0          ! TIME STEP IN A.U.
LOGICAL(4)  :: THAMILTON=.FALSE.  ! HAMILTON MATRIX AVAILABLE
LOGICAL(4)  :: TRAWSTATES=.FALSE. ! PROVIDES NON-DIAGONALIZED WAVE FUNCTIONS THROUGH $GETR
LOGICAL(4)  :: TFORCEX=.TRUE.
LOGICAL(4)  :: TSTRESSX=.FALSE.
!========================================================================
!== PERMANENT DATA, WHICH ARE ORGANIZED BY THE ROUTINES ITSELF         == 
!========================================================================
INTEGER(4)      ,SAVE     :: NKPTL         ! #(LOCAL K-POINTS)
INTEGER(4)      ,POINTER  :: KMAP(:)
REAL(8)                   :: WAVEEKIN1=0.D0
REAL(8)                   :: WAVEEKIN2=0.D0
TYPE(WVSET_TYPE),POINTER  :: THISARRAY(:,:)   ! (NKPTL,NSPIN)
TYPE(WVSET_TYPE),POINTER  :: THIS            ! CURRENT SET OF WAVES
TYPE(GSET_TYPE) ,POINTER  :: GSET            ! CURRENT SET OF GSET
TYPE(MAP_TYPE)            :: MAP
LOGICAL(4)                :: TPR=.FALSE.
LOGICAL(4)                :: TFIRST=.TRUE.
TYPE(EXTERNALPOINTER_TYPE):: EXTPNTR
LOGICAL(4)                :: TFIXRHO=.FALSE.
LOGICAL(4)                :: TWRITERHO=.FALSE.
CHARACTER(8)              :: OPTIMIZERTYPE  ! SWITCH FOR CONJUGATE GRADIENT OR DYNAMICS

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
END MODULE WAVES_MODULE
!
!     ..................................................................
      SUBROUTINE WAVES$STATESELECTED(IB,IKPT,ISPIN,TCHK)
!     ******************************************************************
!     **  WAVES$STATESELECTED                                         **
!     **  TESTS IF THE EXTERNAL POINTER SELECTS A STATE THAT IS       **
!     **  AVAILABLE ON THIS TASK                                      **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IB
      INTEGER(4),INTENT(IN) :: IKPT
      INTEGER(4),INTENT(IN) :: ISPIN
      LOGICAL(4),INTENT(OUT):: TCHK
!     ******************************************************************
      TCHK=.TRUE.
      TCHK=TCHK.AND.EXTPNTR%IKPT.NE.0
      TCHK=TCHK.AND.EXTPNTR%ISPIN.NE.0
      TCHK=TCHK.AND.EXTPNTR%IB.NE.0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SETR8(ID,VAL)
!     ******************************************************************
!     **  WAVES$SETR8A                                                **
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
!     **  WAVES$GETR8A                                                **
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
!     **  WAVES$GETR8A                                                **
!     ******************************************************************
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
        IF(IKPT.EQ.0) THEN
          CALL ERROR$MSG('STATE NOT AVAILABLE ON THIS TASK')
          CALL ERROR$MSG('THIS ERROR OCCURS ONLY FOR PARALLEL JOBS')
          CALL ERROR$MSG("USE WAVES$GETL4('AVAILABLESTATE',TCHK) TO EXPLORE")
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
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
        VAL(:)=THIS%EIGVAL(:)
!
!     ================================================================
!     ==  ENERGY EXPECTATIONS VALUES OF THE ONE-PARTCLE STATES      ==
!     ==  NOTE: IKPT,ISPIN MUST BE SET ON LINKEDLIST                ==
!     ================================================================
      ELSE IF(ID.EQ.'<PSI|H|PSI>') THEN
        IKPT=EXTPNTR%IKPT
        IF(IKPT.EQ.0) THEN
          CALL ERROR$MSG('STATE NOT AVAILABLE ON THIS TASK')
          CALL ERROR$MSG('THIS ERROR OCCURS ONLY FOR PARALLEL JOBS')
          CALL ERROR$MSG("USE WAVES$GETL4('AVAILABLESTATE',TCHK) TO EXPLORE")
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF
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
!     ================================================================
!     ==  GET PROJECTIONS                                           ==
!     ================================================================
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
            VAL(LMN)=REAL(CWORK1(LMN),KIND=8)
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
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
      INTEGER(4)              :: IKPT,IKPTL,I
      INTEGER(4)              :: THISTASK,NTASKS
!     ******************************************************************
      IF(ID.EQ.'SPINORDIM') THEN
        NDIM=VAL
!
      ELSE IF(ID.EQ.'IKPT') THEN
        IKPT=VAL
!       == CHECK IF K-POINT IS PRESENT ===============================
        CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
        IF(KMAP(IKPT).NE.THISTASK) THEN
          IKPTL=0
        ELSE
!         == CONVERT GLOBAL IKPT INTO LOCAL IKPT =====================
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
        IF(VAL.NE.'DYNAMIC'.AND.VAL.NE.'STATIC') THEN
          CALL ERROR$MSG('VALUE NOT RECOGNIZED. ALLOWED VALUES ARE "STATIC" AND "DYNAMIC"')
          CALL ERROR$CHVAL('VAL',VAL)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('WAVES$SETCH')
        END IF
        RSTRTTYPE=VAL
      ELSE IF(ID.EQ.'OPTIMIZER') THEN
        IF(VAL.NE.'CG'.AND.VAL.NE.'MD') THEN
          CALL ERROR$MSG('VALUE NOT RECOGNIZED. ALLOWED VALUES ARE "CG" AND "MD"')
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
!     ..................................................................
      SUBROUTINE WAVES$INITIALIZE
!     ******************************************************************
!     **  GENERATE G-VECTORS AND OTHER INITIALIZATION                 **
!     ******************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      REAL(8)                :: RBAS(3,3) ! LATTICE VECTORS
      REAL(8)                :: GBAS(3,3) ! RECIPROCAL LATTICE VECTORS
      REAL(8)                :: CELLVOL   ! UNIT CELL  VOLUME
      REAL(8)   ,ALLOCATABLE :: XK(:,:)   ! K-POINTS IN RELATIVE COORDINATES
      REAL(8)   ,ALLOCATABLE :: G2(:)     ! G**2
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
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES$INITIALIZE')
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
        MAP%LNX(ISP)=LNX
        CALL SETUP$LOFLN(ISP,LNX,MAP%LOX(1:LNX,ISP))
        MAP%LMNX(ISP)=0
        DO LN=1,LNX
          MAP%LMNX(ISP)=MAP%LMNX(ISP)+2*MAP%LOX(LN,ISP)+1
          MAP%LMX=MAX(MAP%LMX,(MAP%LOX(LN,ISP)+1)**2)
        ENDDO
        MAP%NBAREPRO=MAP%NBAREPRO+LNX
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
CALL FILEHANDLER$UNIT('PROT',NFILO)
!WRITE(NFILO,*)'TINV FORCED TO BE FALSE IN WAVES$INITIALIZE!!!'
!PRINT*,'TINV FORCED TO BE FALSE IN WAVES$INITIALIZE!!!'
!TINV=.FALSE.
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
      CALL PLANEWAVE$REPORT(NFILO)
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
!     ..................................................................
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
      IF(ID.EQ.'OCC') THEN
        CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
        CALL DYNOCC$GETI4('NB',NBX)
        CALL DYNOCC$GETI4('NSPIN',NSPIN)
        ALLOCATE(VALG(NBX*NKPT*NSPIN))
        CALL DYNOCC$GETR8A('OCC',NBX*NKPT*NSPIN,VALG)
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
      ELSE IF(ID.EQ.'WKPT') THEN
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
!     ..................................................................
      SUBROUTINE WAVES_DYNOCCSETR8A(ID,LEN,VAL)
!     **                                                              **
!     **  SET A REAL(8) ARRAY TO DYNOCC. THIS ROUTINE IS AN         **
!     **  INTERFACE BETWEEN WAVES OBJECT AND DYNOCC OBJECT.           **
!     **  IT IS REQUIRED FOR THE K-POINT PARALLELIZATION, BECAUSE     **
!     **  EACH NODE ONLY MAINTAINS A FRACTION OF ALL K-POINTS         **
!     **  WHILE THE DYNOCC OBJECT WORKS ON ALL K-POINTS               **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   ASSUMES THAT THE INPUT ARRAY IS IDENTICAL ON ALL           **
!     **      TASKS OF THE LOCAL K-GROUP                              **
!     **   KMAP(NKPT) POINTS TO THE FIRST TASK IN THE K-GROUP WHERE   **
!     **     THE K-POINT RESIDES (THE TASK IS RELATIVE TO THE MONOMER **
!     **     GROUP)                                                   **
!     **                                                              **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
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
!     ******************************************************************
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
!     ..................................................................
      SUBROUTINE WAVES$ETOT()
!     ******************************************************************
!     **                                                              **
!     **  EVALUATE PS KINETIC ENERGY IN G-SPACE                       **
!     **  EVALUATE NUMBER OF ELECTRONS IN G-SPACE                     **
!     **                                                              **
!     *******************************************P.E. BLOECHL, (2000)***
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: LMRXX
      REAL(8)   ,ALLOCATABLE :: QLM(:,:)  !(LMRXX) MULTIPOLE MOMENTS
      REAL(8)   ,ALLOCATABLE :: VQLM(:,:) !(LMRXX) MULTIPOLE POTENTIALS
      REAL(8)   ,ALLOCATABLE :: RHO(:,:)  ! CHARGE DENSITY
      COMPLEX(8),ALLOCATABLE :: DENMAT(:,:,:,:) ! 1CENTER DENSITY MATRIX
      COMPLEX(8),ALLOCATABLE :: EDENMAT(:,:,:,:)! ENERGY-WEIGHTED DENSITY MATRIX
      REAL(8)   ,ALLOCATABLE :: DH(:,:,:,:)     ! 1CENTER HAMILTONIAN
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
      REAL(8)   ,PARAMETER   :: TINY=1.D-300
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
      REAL(8)                :: SVAR1,SVAR2
      COMPLEX(8),ALLOCATABLE :: HAMILTON(:,:)
      REAL(8)   ,ALLOCATABLE :: EIG(:,:,:)
      LOGICAL(4)             :: TSTRESS
      LOGICAL(4)             :: TFORCE
      LOGICAL(4)             :: TCHK
      COMPLEX(8),ALLOCATABLE :: QMAT(:,:)   
      INTEGER(4)             :: NFILO
      LOGICAL(4)             :: TCONV ! MIXER SAYS THAT WAVE FUNCTIONS ARE CONVERGED !KAESTNERCG
      REAL(8)                :: CONVPSI ! CONVERGENCE CRITERION FOR WAVE FUNCTIONS !KAESTNERCG
!     ******************************************************************      
INTEGER(4) ::NTASKS_W,THISTASK_W
CALL MPE$QUERY('~',NTASKS_W,THISTASK_W)
                              CALL TRACE$PUSH('WAVES$ETOT')
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
!     == NUMBER OF BANDS ===============================================
      CALL DYNOCC$GETI4('NB',NBX)
!
!     ==================================================================
!     == INITIALIZE GSET: YLM, PRO, EIGR                              ==
!     ==================================================================
      IF(TFIRST.OR.TSTRESS) THEN
        CALL WAVES_UPDATEGSET()
      END IF
!
!     ==================================================================
!     == DEALLOCATE EIGENVALUES AND EIGENVECTORS                      ==
!     ==================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(ASSOCIATED(THIS%EIGVAL))DEALLOCATE(THIS%EIGVAL)
          IF(ASSOCIATED(THIS%EIGVEC))DEALLOCATE(THIS%EIGVEC)
        ENDDO
      ENDDO
      IF(OPTIMIZERTYPE.EQ.'CG') THEN   !KAESTNERCG
        CALL WAVES$KAESTNERCG1(TFIRST)
      END IF
!
!     ==================================================================
!     == RANDOMIZE INITIAL WAVE FUNCTIONS                             ==
!     ==================================================================
      IF(TFIRST.AND.TRANDOM) THEN
        CALL WAVES$RANDOMIZE()
      END IF
!
!     ==================================================================
!     == GRAMM-SCHMIDT ORTHOGONALIZATION OF INITIAL WAVE FUNCTIONS    ==
!     ==================================================================
      IF(TFIRST) THEN
        CALL WAVES$GRAMMSCHMIDT()
      END IF
!
!     ==================================================================
!     == CALCULATE PROJECTIONS                                        ==
!     ==================================================================
      IF(TFIRST) THEN
        CALL WAVES$PROJECTIONS('PSI0')
      END IF


      IF(TFIRST) TFIRST=.FALSE.
                              CALL TIMING$CLOCKON('WAVES$ETOT')
!
!     ==================================================================
!     == KINETIC ENERGY                                               ==
!     ==================================================================
      CALL WAVES$EKIN(EKIN,STRESSKIN)   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      STRESS(:,:)=STRESS(:,:)+STRESSKIN(:,:)
      CALL ENERGYLIST$SET('PS  KINETIC',EKIN)
      CALL ENERGYLIST$ADD('AE  KINETIC',EKIN)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EKIN)
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
      NAT=MAP%NAT
      LMNXX=0
      DO ISP=1,MAP%NSP
        LMNXX=MAX(LMNXX,MAP%LMNX(ISP))
      ENDDO
      ALLOCATE(DENMAT(LMNXX,LMNXX,NDIMD,NAT))
      ALLOCATE(EDENMAT(LMNXX,LMNXX,NDIMD,NAT))
      CALL WAVES$DENMAT(LMNXX,NDIMD,NAT,DENMAT,EDENMAT) !<<<<<<<<<<<<<<<
!
!     ==================================================================
!     == PSEUDO DENSITY STILL WITHOUT PSEUDOCORE                      ==
!     == SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                      ==
!     == SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)               ==
!     == NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z) ==
!     ==================================================================
      NRL=MAP%NRL
      ALLOCATE(RHO(NRL,NDIMD))
      CALL WAVES$RHO(NRL,NDIMD,RHO)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
  CALL MPE$COMBINE('MONOMER','+',SVAR1)   !LADUNGSDICHTE EINSAMMELN
  PRINT*,'TOTAL CHARGE IN PSEUDO WAVE FUNCTIONS W/O PSCORE ',SVAR1
END IF
!
!     ==================================================================
!     == MIX RHO AND CHECK IF CONVERGED                               ==
!     ==================================================================
      IF(OPTIMIZERTYPE.EQ.'CG') THEN                                             !KAESTNERCG
         CALL TIMING$CLOCKON('MIXER')                                        !KAESTNERCG
         TSTOP=.TRUE. ! DO NOT PROPAGATE THE WAVE FUNCTION                   !KAESTNERCG
                      ! MUST BE SET IN EACH ITERATION, SEE WAVES$PROPAGATE   !KAESTNERCG
 PRINT*,'----------------- NEW CG SCF ITERATION------------------'           !KAESTNERCG
CALL ERROR$MSG('MIXER SWITCHED OFF BECAUSE IT IS STILL INCOMPATIBLE')
CALL ERROR$MSG('WITH A COMPLEX DENSITY MATRIX')
CALL ERROR$STOP('WAVES$ETOT')
!         CALL WAVES_MIXRHO_PUL(NRL*NDIM,LMNXX*LMNXX*NDIMD*NAT,&              !KAESTNERCG
!              RHO,REAL(DENMAT,KIND=8),TCONV,CONVPSI)                         !KAESTNERCG
         CALL TIMING$CLOCKOFF('MIXER')                                       !KAESTNERCG
      END IF                                                                 !KAESTNERCG
!
!     ==================================================================
!     == MULTIPOLE MOMENTS OF ONE-CENTER PART                         ==
!     == CONTAINS CORE AND PSEUDOCORE                                 ==
!     ==================================================================
      CALL SETUP$LMRXX(LMRXX)
      ALLOCATE(QLM(LMRXX,NAT))
      NAT=MAP%NAT
      CALL WAVES$MOMENTS(NAT,LMRXX,NDIMD,LMNXX,DENMAT,QLM)
!
!     ==================================================================
!     == OVERWRITE DENSITY TO FOR FIXED POTENTIAL CALC.               ==
!     ==================================================================
      IF(TFIXRHO) THEN
        CALL WAVES_FIXRHOGET(NRL,NDIMD,LMRXX,NAT,QLM,RHO,DENMAT)
      END IF
      IF(TWRITERHO) THEN
        CALL WAVES_FIXRHOSET(NRL,NDIMD,LMRXX,NAT,QLM,RHO,DENMAT)
      END IF
!
!     ==========================================================================
!     == ANALYSE SPIN DENSITY                                                 ==
!     ==========================================================================
      IF(THAMILTON) THEN
        CALL WAVES$SPINS(NRL,NDIMD,RHO,NAT,LMNXX,DENMAT)
      END IF
!
!     ==================================================================
!     == POTENTIAL (POTENTIAL IS STORED BACK INTO THE DENSITY ARRAY!) ==
!     ==================================================================
      ALLOCATE(VQLM(LMRXX,NAT))
      CALL WAVES_VOFRHO(NRL,NDIMD,RHO,RHOB,NAT,LMRXX,QLM,VQLM)
      DEALLOCATE(QLM)

!      ALLOCATE(FORCET(3,NAT))
!      ALLOCATE(VQLM(LMRXX,NAT))
!      CALL POTENTIAL$VOFRHO(NRL,NDIMD,RHO,LMRXX,NAT,QLM,VQLM &
!     &                     ,R,FORCET,RBAS,STRESST,RHOB)
!      FORCE=FORCE+FORCET
!      STRESS=STRESS+STRESST
!      DEALLOCATE(QLM)
!      DEALLOCATE(FORCET)
!
!     ==================================================================
!     == AUGMENTATION                                                 ==
!     ================================================================== 
      ALLOCATE(DH(LMNXX,LMNXX,NDIMD,NAT))
      ALLOCATE(DO(LMNXX,LMNXX,NDIMD,NAT))
      CALL WAVES$SPHERE(LMNXX,NDIMD,NAT,LMRXX,RHOB,DENMAT,EDENMAT &
     &                 ,VQLM,DH,DO,POTB)
!
!     ===================================================================
!     ==  SUBTRACT AVERAGE ELECTROSTATIC POTENTIAL =====================
!     ==  TO ACCOUNT FOR BACKGROUND DENSITY ============================
!     ==  POTB=-1/OMEGA*INT D^3R(\TILDE{V}+V^1-\TILDE{V}^1) ============
!     ==  NOTE THAT \INT D^3R \TILDE{V}=0
!     ===================================================================
print*,'before waves$addconstantpot'
      CALL WAVES$ADDCONSTANTPOT(NRL,LMNXX,NDIMD,NAT,POTB,RHO,DH,DO)
      CALL GRAPHICS$SETR8('POTSHIFT',POTB)
!
!     ==================================================================
!     == COMMUNICATE DATA WITH OPTICS MODULE                          ==
!     ==================================================================
!     CALL OPTICS$GETL4('ON',TCHK)
!     IF(TCHK) THEN
!        CALL MPE$COMBINE('NONE','+',DO) ! DO CURRENTLY USED ONLY FOR OPTICS
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
             LN,(REAL(DENMAT(LMN1+M,LMN1+M,ISPIN,IAT),KIND=8),M=1,2*L+1)
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
!             LN,(REAL(DENMAT(LMN1+M,LMN1+M,ISPIN,IAT),KIND=8),M=1,2*L+1)
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
      DEALLOCATE(VQLM)
      DEALLOCATE(DENMAT)
      DEALLOCATE(EDENMAT)
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
      IF(TFORCE.OR.TSTRESS) THEN
        CALL WAVES$FORCE(NAT,LMNXX,NDIMD,DH,FORCE,STRESS)
      END IF
!
!     ==================================================================
!     ==  CALL GONJUGATE GRADIENT                                     ==
!     ================================================================== 
      IF(OPTIMIZERTYPE.EQ.'CG') THEN                                             !KAESTNERCG
         CALL TIMING$CLOCKON('CG')                                           !KAESTNERCG
         CALL CG$STATE_BY_STATE(MAP%NRL,NDIMD,RHO(:,:),CONVPSI,NAT,LMNXX,DH) !KAESTNERCG
         CALL TIMING$CLOCKOFF('CG')                                          !KAESTNERCG
         TCONV=.FALSE. ! TCONV HAS NOT BEEN SET!!!
         IF(TCONV) CALL STOPIT$SETL4('STOP',.TRUE.)                          !KAESTNERCG
      END IF                                                                 !KAESTNERCG
!
!     ==================================================================
!     ==  EVALUATE H*PSI                                              ==
!     ================================================================== 
!PRINT*,'RHO',(SUM(ABS(RHO)).GT.0.D0.OR.SUM(ABS(RHO)).LE.0.D0)
      CALL WAVES$HPSI(NRL,NDIMD,NAT,LMNXX,RHO,DH)
      DEALLOCATE(RHO)
      DEALLOCATE(DH)
!
!     ==================================================================
!     ==  EVALUATE ENERGY EXPECTATION VALUES                          ==
!     ================================================================== 
CALL TIMING$CLOCKON('W:EXPECT')
      ALLOCATE(EIG(NBX,NKPTL,NSPIN))
      ALLOCATE(HAMILTON(2,2))
      DO IKPT=1,NKPTL
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
              EIG(2*IB-1,IKPT,ISPIN)=REAL(HAMILTON(1,1),KIND=8)
              EIG(2*IB  ,IKPT,ISPIN)=REAL(HAMILTON(2,2),KIND=8)
            ELSE 
              CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,1,1 &
     &            ,THIS%PSI0(:,:,IB),THIS%HPSI(:,:,IB),HAMILTON)
              EIG(IB,IKPT,ISPIN)=REAL(HAMILTON(1,1),KIND=8)
            END IF
          ENDDO
          THIS%EXPECTVAL(:)=EIG(1:NB,IKPT,ISPIN)
        ENDDO
      ENDDO
      DEALLOCATE(HAMILTON)
!     == OCCUPATIONS ===================================================
      CALL WAVES_DYNOCCSETR8A('EPSILON',NB*NKPTL*NSPIN,EIG)
      DEALLOCATE(EIG)
!
!     ==================================================================
!     ==  EVALUATE <S^2>                                              ==
!     ==================================================================
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
          CALL WAVES_TOTALSPIN(NB,NKPT,IKPT,NSPIN,OCC,QMAT)
          DEALLOCATE(QMAT)
        END DO
        DEALLOCATE(OCC)
        CALL TIMING$CLOCKOFF('S2')
        CALL FILEHANDLER$UNIT('PROT',NFILO)
        CALL WAVES$REPORTSPIN(NFILO)
      END IF
!
!     ==================================================================
!     ==  EVALUATE <PSI|H|PSI>                                        ==
!     ================================================================== 
      IF(THAMILTON) THEN
        DO IKPT=1,NKPTL
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
        IF(OPTIMIZERTYPE.NE.'CG') THEN    !KAESTNERCG
          IF(ASSOCIATED(THIS%EIGVAL))DEALLOCATE(THIS%EIGVAL)
          IF(ASSOCIATED(THIS%EIGVEC))DEALLOCATE(THIS%EIGVEC)
        END IF                         !KAESTNERCG
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
      DEALLOCATE(R)
                              CALL TIMING$CLOCKOFF('WAVES$ETOT')
                              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_VOFRHO(NRL,NDIMD,RHO,RHOB,NAT,LMRXX,QLM,VQLM)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
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
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES_VOFRHO')

!     == ASSUMES THAT NR2 AND NR3 ARE IDENTICAL FOR DENSITY AND FOR
!     == ALL WAVES IN THE GROUP/ NR1L ON THE OTHER HAND MAY DIFFER
      CALL PLANEWAVE$SELECT('DENSITY')
      CALL PLANEWAVE$GETI4('NR1L',NR1L_V)
      CALL PLANEWAVE$GETI4('NR2',NR2)
      CALL PLANEWAVE$GETI4('NR3',NR3)
      NR1L=NRL/(NR2*NR3)  ! NRL=NR1L*NR2*NR3  REFERS TO THE WAVE FUNCTION!
!
!     ==================================================================
!     == PERFORM SUM OVER K-POINT GROUPS AND REDISTRIBUTE             ==
!     ==================================================================
      ALLOCATE(RHO_V(NR1L_V*NR2*NR3,NDIMD))
      CALL WAVES_MAPPSITOPOT('PSITOPOT',NR1L,NR1L_V,NR2,NR3,NDIMD,RHO,RHO_V)
!
!     ==================================================================
!     == CALCULATE POTENTIAL FROM THE CHARGE DENSITY                  ==
!     ==================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
      FORCE(:,:)=0.D0
      STRESS(:,:)=0.D0
      VQLM(:,:)=0.D0
      CALL POTENTIAL$VOFRHO(NR1L_V*NR2*NR3,NDIMD,RHO_V,LMRXX,NAT,QLM,VQLM &
     &                     ,R,FORCE,RBAS,STRESS,RHOB)
!
!     ==================================================================
!     == MAP POTENTIAL BACK ONTO THE WAVE FUNCTION GRID               ==
!     ==================================================================
      CALL WAVES_MAPPSITOPOT('POTTOPSI',NR1L,NR1L_V,NR2,NR3,NDIMD,RHO,RHO_V)
      DEALLOCATE(RHO_V)
!
!     ==================================================================
!     ==  STORE FORCES AND STRESSES BACK TO THE OWNING OBJECTS        ==
!     ==================================================================
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
!     ..................................................................
      SUBROUTINE WAVES_MAPPSITOPOT(ID,NR1L_P,NR1L_V,NR2,NR3,NDIMD,RHO_P,RHO_V)
!     ******************************************************************
!     **                                                              **
!     **  ADDS UP THE DENSITY OVER K-POINTS AND MAPS THE DENSITY      **
!     **  FROM THE SHEETS OF THE WAVE FUNCTIONS INTO THOSE OF THE     **
!     **  POTENTIAL AND VICE VERSA.                                   **
!     **                                                              **
!     **  THE DENSITY IS DIVIDED INTO SHEETS, WHERE EACH SHEET IS     **
!     **  KEPT BY ONE NODE. IN THE K-POINT PARALLELIZATION, THE       **
!     **  SHEETS DIFFER FOR THE WAVE FUNCTIONS AND THE DENSITY,       **
!     **  BECAUSE A K-GROUP MUST KEEP ALL GRID POINTS OF THE DENSITY  **
!     **  WHILE FOR THE POTENTIAL, THEY ARE DISTRIBUTED OVER ALL      **
!     **  NODES OF THE MONOMER GROUP.                                 **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **    THE DIVISION IN SHEETS IS THE SAME FOR ALL KPOINTS WITHIN **
!     **    ONE K-GROUP                                               **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
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
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES_MAPPSITOPOT')
      IF(ID.NE.'PSITOPOT'.AND.ID.NE.'POTTOPSI') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('WAVES_MAPPSITOPOT')
      ENDIF
!
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)      
!
!     ====================================================================
!     == DETERMINE WHICH PLANES OF RHO_V ARE HELD BY EACH PROCESSOR   ==
!     ====================================================================
      CALL PLANEWAVE$SELECT('DENSITY')
      ALLOCATE(NR1LARR(NTASKS))
      CALL PLANEWAVE$GETI4A('NR1LARR',NTASKS,NR1LARR)
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
!     ====================================================================
!     == COMMUNICATE WHICH PLANES OF RHO_P ARE HELD BY EACH PROCESSOR   ==
!     ====================================================================
!     __ IT IS POSSIBLE THAT THERE IS NO K-POINT ON A TASK
      IF(NKPTL.GT.0) THEN
        CALL WAVES_SELECTWV(1,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETI4('NR1START',NR1START_P)
      ELSE
        NR1START_P=1
      END IF
      ALLOCATE(NR1FIRST_P(NTASKS))
      ALLOCATE(NR1LAST_P(NTASKS))
      DO ITASK=1,NTASKS
        IWORK2(1)=NR1START_P
        IWORK2(2)=NR1START_P-1+NR1L_P
        CALL MPE$BROADCAST('MONOMER',ITASK,IWORK2(:))
        NR1FIRST_P(ITASK)=IWORK2(1)
        NR1LAST_P(ITASK)=IWORK2(2)
      ENDDO
!
!     ====================================================================
!     ==                                                                ==
!     ====================================================================
      IF(ID.EQ.'PSITOPOT') THEN
        ALLOCATE(ARR(NR1L_V,NR2,NR3,NDIMD))
        RHO_V(:,:,:,:)=0.D0
      END IF
!
!     ====================================================================
!     ==  DISTRIBUETE THE DENSITY/POTENTIAL                             ==
!     ====================================================================
      DO ITASK_P=1,NTASKS
        DO ITASK_V=1,NTASKS
          IF(.NOT.(ITASK_P.EQ.THISTASK.OR.ITASK_V.EQ.THISTASK)) CYCLE
!
          IRSTART=MAX(NR1FIRST_V(ITASK_V),NR1FIRST_P(ITASK_P)) 
          IREND  =MIN(NR1LAST_V(ITASK_V),NR1LAST_P(ITASK_P)) 
          IF(IRSTART.GT.IREND) CYCLE  ! NO MESSAGE, GO TO NEXT....
!
!         ==============================================================
!         ==                                                          ==
!         ==============================================================
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
                CALL MPE$SENDRECEIVE('MONOMER',ITASK_V,ITASK_P,RHO_V(IR1_V:IR2_V,:,:,:))
              ELSE
!               __ RECEIVE______
!PRINT*,'RECEIVE: ', THISTASK,':',ITASK_V,'->',ITASK_P
                CALL MPE$SENDRECEIVE('MONOMER',ITASK_V,ITASK_P,RHO_P(IR1_P:IR2_P,:,:,:))
              END IF
            END IF
          END IF
        ENDDO
      ENDDO      
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
!     ..................................................................
      SUBROUTINE WAVES$KAESTNERCG1(TFIRST_)
!     ******************************************************************
!     **                                                              **
!     **  THIS IS FOR KAESTNER'S CONJUGATE GRADIENT IMPLEMENTATION    **
!     **  HANDLE WITH SPECIAL CARE TO AVOID PROBLEMS WITH ALLOCATION  **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TFIRST_
      INTEGER(4)            :: IKPT,ISPIN
      INTEGER(4)            :: NB
!     ******************************************************************
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
!     ..................................................................
      SUBROUTINE WAVES$RANDOMIZE()
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE WAVES_MODULE !, ONLY : NSPIN,NKPTL,NDIM,GSET,THIS,AMPRANDOM
      IMPLICIT NONE
      INTEGER(4)             :: IKPT,ISPIN
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NBH
      REAL(8)   ,ALLOCATABLE :: G2(:)
!     ******************************************************************
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
!     ..................................................................
      SUBROUTINE WAVES$GRAMMSCHMIDT()
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: IKPT,ISPIN
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NB
      INTEGER(4)             :: NBH
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R(:,:)
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES$GRAMMSCHMIDT')
                              CALL TIMING$CLOCKON('W:GRAMSCHMIDT')
      NAT=MAP%NAT
      ALLOCATE(R(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NB=THIS%NB
          NBH=THIS%NBH
          CALL WAVES_GRAMSCHMIDT(MAP,GSET,NAT,R,NGL,NDIM,NBH,NB,THIS%PSI0)
          CALL WAVES_GRAMSCHMIDT(MAP,GSET,NAT,R,NGL,NDIM,NBH,NB,THIS%PSIM)
        ENDDO
      ENDDO
                              CALL TIMING$CLOCKOFF('W:GRAMSCHMIDT')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$GRAMMSCHMIDT
!
!     ..................................................................
      SUBROUTINE WAVES$PROJECTIONS(ID)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN):: ID
      INTEGER(4)             :: IKPT,ISPIN
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NBH
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R(:,:)
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES$PROJECTIONS')
                              CALL TIMING$CLOCKON('W:PROJ')
      NAT=MAP%NAT
      ALLOCATE(R(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          IF(ID.EQ.'PSI0') THEN
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,MAP%NPRO &
     &                            ,THIS%PSI0,THIS%PROJ)
            CALL MPE$COMBINE('K','+',THIS%PROJ)
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
!     ..................................................................
      SUBROUTINE WAVES$EKIN(EKIN,STRESS)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES PSEUDO-KINETIC ENERGY                             **
!     **  AND THE CORRESPONDING STRESS                                **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
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
!     ******************************************************************      
                              CALL TRACE$PUSH('WAVES$EKIN')
                              CALL TIMING$CLOCKON('W:EKIN')
!
!     ==================================================================
!     ==  GET OCCUPATIONS FROM DYNOCC OBJECT                          ==
!     ==================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
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
!     *******************************************P.E. BLOECHL, (1999)***
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
    &                 +FP*REAL(CONJG(PSI(IG,IDIM,IB))*PSI(IG,IDIM,IB),KIND=8)
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
     &                +FM*REAL(CONJG(PSI(IG,IDIM,IB))*PSI1(IG,IDIM),KIND=8)
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
!       ==  WEIGHT RESULT AND CLOSE DOWN ===============================
        EKIN=EKIN*GWEIGHT
        DEALLOCATE(DMAT)
        DEALLOCATE(G2)
        STRESS(:,:)=0.D0
      END IF
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$DENMAT(LMNXX,NDIMD_,NAT,DENMAT,EDENMAT)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES ONE-CENTER DENSITY MATRIX FROM THE ACTUAL         **
!     **  PROJECTIONS <PRO|PSI> IN THE WAVES OBJECT                   **
!     **                                                              **
!     **  ONE-CENTER DENSITY MATRICES                                 **
!     **  SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                     **
!     **  SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)              **
!     **  NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z)**
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE MPE_MODULE
      USE WAVES_MODULE
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
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      COMPLEX(8)             :: CSVAR1,CSVAR2
      INTEGER(4)             :: NTASKS,THISTASK
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES$DENMAT')
                              CALL TIMING$CLOCKON('W:DENMAT')
      IF(NDIMD_.NE.NDIMD.OR.NAT.NE.MAP%NAT) THEN 
        CALL ERROR$MSG('ARRAY SIZE INCONSISTENT')
        CALL ERROR$STOP('WAVES$DENMAT')
      END IF
!
!     ==================================================================
!     ==  GET OCCUPATIONS FROM DYNOCC OBJECT                          ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      CALL MPE$QUERY('K',NTASKS,THISTASK)
      DENMAT(:,:,:,:)=(0.D0,0.D0)
      EDENMAT(:,:,:,:)=(0.D0,0.D0)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NBH=THIS%NBH
          NB=THIS%NB
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
            ALLOCATE(PROJ(NDIM,NBH,LMNX))
            PROJ(:,:,:)=THIS%PROJ(:,:,IPRO:IPRO-1+LMNX)
            ALLOCATE(DENMAT1(LMNX,LMNX,NDIMD))
            ALLOCATE(EDENMAT1(LMNX,LMNX,NDIMD))
            CALL WAVES_DENMAT(NDIM,NBH,NB,LMNX,OCC(1,IKPT,ISPIN),THIS%RLAM0 &
     &                       ,PROJ,DENMAT1,EDENMAT1)
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
            IPRO=IPRO+LMNX
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(OCC)
!     == THE PROJECTIONS ARE IDENTICAL AND COMPLETE FOR EACH K-GROUP
!     == EACH K-GROUP HOLDS ONLY THE WAVE FUNCTIONS BELONGING TO IT.
!     == EACH PROCESSOR OF EACH K-GROUP ONLY ADDS UP A FRACTION OF THE PROJECTIONS
!     == THEREFORE THERE IS NO DOUBLE COUNTING BY SUMMING OVER THE MONOMER
      CALL MPE$COMBINE('MONOMER','+',DENMAT)
      CALL MPE$COMBINE('MONOMER','+',EDENMAT)
!
!     ==================================================================
!     ==  CONVERT SPIN-UP AND SPIN-DOWN DENSITY MATRIX INTO           ==
!     ==  TOTAL AND SPIN DENSITY MATRICES                             ==
!     ==================================================================
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
!
!     ==================================================================
!     ==  PRINT FOR TEST                                              ==
!     ==================================================================
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
!
                              CALL TIMING$CLOCKOFF('W:DENMAT')
                              CALL TRACE$POP
     RETURN
     END SUBROUTINE WAVES$DENMAT
!
!     ..................................................................
      SUBROUTINE WAVES_DENMAT(NDIM,NBH,NB,LMNX,OCC,LAMBDA,PROPSI &
     &                       ,DENMAT,EDENMAT)
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
!     *******************************************P.E. BLOECHL, (1999)***
!RELEASED: 8.OCT.99
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NDIM   ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN) :: NBH    ! #(WAVE FUNCTIONS
      INTEGER(4),INTENT(IN) :: NB     ! #(STATES)
      INTEGER(4),INTENT(IN) :: LMNX   ! #(PROJECTORS ON THIS SITE)
      REAL(8)   ,INTENT(IN) :: OCC(NB)! OCCUPATIONS
      COMPLEX(8),INTENT(IN) :: LAMBDA(NB,NB)   !LAGRANGE/F
      COMPLEX(8),INTENT(IN) :: PROPSI(NDIM,NBH,LMNX) !<PRO|PSI>
      COMPLEX(8),INTENT(OUT):: DENMAT(LMNX,LMNX,NDIM**2)
      COMPLEX(8),INTENT(OUT):: EDENMAT(LMNX,LMNX,NDIM**2)
      COMPLEX(8)            :: DENMAT1(LMNX,LMNX,NDIM,NDIM)
      COMPLEX(8)            :: EDENMAT1(LMNX,LMNX,NDIM,NDIM)
      COMPLEX(8)            :: FUNC(LMNX,NDIM)
      COMPLEX(8)            :: LAGR(NB,NB)
      LOGICAL(4)            :: TINV
      INTEGER(4)            :: LMN1,LMN2,IDIM1,IDIM2,IB,IB1,IB2
      REAL(8)               :: SVAR1,SVAR2
      COMPLEX(8)            :: CSVAR,CFACR,CFACI,CFAC1,CFAC2
      INTEGER(4)            :: IDIM,NDIMD
      INTEGER(4)            :: NFILO
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
!     ******************************************************************
      NDIMD=NDIM**2
      DENMAT(:,:,:)=(0.D0,0.D0)
      EDENMAT(:,:,:)=(0.D0,0.D0)
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
!     == SUPER WAVE FUNCTIONS FOR TINV=TRUE
!     ==   <P|PSI1+I*PSI2>(F1+F2)/2<PSI1-I*PSI2|P>
!     == + <P|PSI1+I*PSI2>(F1-F2)/2<PSI1+I*PSI2|P>
!     == = <P|PSI1>F1<PSI1|P>+<P|PSI2>F2<PSI2|P>
!     == + I[<P|PSI2>F1<PSI1|P>-<P|PSI1>F2<PSI2|P>]
!     == IMAGINARY PART IS DROPPED
      DENMAT1(:,:,:,:)=(0.D0,0.D0)
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
!
!     ==================================================================
!     == NOW SUM UP EDENMAT  <P|PSITILDE>*LAMBDA*<PSITILDE|P>         ==
!     ==================================================================
      DO IB2=1,NB
        LAGR(:,IB2)=LAMBDA(:,IB2)*OCC(IB2)
      ENDDO
      EDENMAT1(:,:,:,:)=(0.D0,0.D0)
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
     &                         +CFAC2*CONJG(PROPSI(IDIM,IB2,:))
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
     &                    +PROPSI(IDIM1,IB1,LMN1)*FUNC(LMN2,IDIM2)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==================================================================
!     == MAP DENSITY MATRIX ONTO TOTAL AND SPIN DENSITY               ==
!     ==================================================================
      IF(NDIM.EQ.1) THEN  !== TOTAL DENSITY ===========================
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            DENMAT(LMN1,LMN2,1)=DENMAT1(LMN1,LMN2,1,1)
            EDENMAT(LMN1,LMN2,1)=EDENMAT1(LMN1,LMN2,1,1)
          ENDDO
        ENDDO
      ELSE IF(NDIM.EQ.2) THEN  !== TOTAL DENSITY, X,Y,Z SPIN DENSITY ==
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
!     ==================================================================
!     == SYMMETRIZE DENSITY MATRIX                                    ==
!     ==================================================================
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
!     ==================================================================
!     == PRINOUT FOR TEST                                             ==
!     ==================================================================
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
!     ..................................................................
      SUBROUTINE WAVES$RHO(NRL,NDIMD_,RHO)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES PSEUDO-DENSITY FROM THE ACTUAL PSEUDO WAVE        **
!     **  FUNCTIONS                                                   **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NRL
      INTEGER(4),INTENT(IN)  :: NDIMD_
      REAL(8)   ,INTENT(OUT) :: RHO(NRL,NDIMD_)
      REAL(8)                :: RHO1(NRL,NDIMD_)
      INTEGER(4)             :: IKPT,ISPIN,IR
      INTEGER(4)             :: NBX        
      REAL(8)  ,ALLOCATABLE  :: OCC(:,:,:) 
      REAL(8)                :: SVAR1,SVAR2
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES$RHO')
                              CALL TIMING$CLOCKON('W:RHO')
      IF(NDIMD_.NE.NDIMD.OR.NRL.NE.MAP%NRL) THEN
        CALL ERROR$MSG('ARRAY SIZE INCONSISTEN')
        CALL ERROR$STOP('WAVES$RHO')
      END IF
!
!     ==================================================================
!     ==  GET OCCUPATIONS FROM DYNOCC OBJECT                          ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
!
!     ==================================================================
!     ==  CALCULATE DENSITY                                           ==
!     ==================================================================
      RHO(:,:)=0.D0
      DO IKPT=1,NKPTL
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
!
!     ===================================================================
!     == CONVERT INTO TOTAL AND SPIN DENSITY                           ==
!     ===================================================================
      IF(NSPIN.EQ.2) THEN
        DO IR=1,NRL
          SVAR1=RHO(IR,1)
          SVAR2=RHO(IR,2)
          RHO(IR,1)=SVAR1+SVAR2   ! TOTAL DENSITY
          RHO(IR,2)=SVAR1-SVAR2   ! SPIN DENSITY
        ENDDO
      END IF
                              CALL TIMING$CLOCKOFF('W:RHO')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$RHO
!
!     ..................................................................
      SUBROUTINE WAVES$MOMENTS(NAT,LMRXX,NDIMD_,LMNXX,DENMAT,QLM)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES MULTIPOLE MOMENTS OF THE ONE-CENTER DENSITY       **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
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
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES$MOMENTS')
                              CALL TIMING$CLOCKON('W:MOMENTS')
      IF(NDIMD_.NE.NDIMD.OR.NAT.NE.MAP%NAT) THEN
        CALL ERROR$MSG('ARRAY SIZE INCONSISTEN')
        CALL ERROR$STOP('WAVES$RHO')
      END IF
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
!     == DENMAT IS IDENTICAL AND COMPLETE ON EACH PROCESS OF A MONOMER
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      QLM(:,:)=0.D0
      DO IAT=THISTASK,NAT,NTASKS   ! DISTRIBUTE WORK ACCROSS TASKS
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        CALL SETUP$LMRX(ISP,LMRX)
        ALLOCATE(DENMAT1(LMNX,LMNX))
        DENMAT1(:,:)=DENMAT(1:LMNX,1:LMNX,1,IAT)
        CALL AUGMENTATION$MOMENTS(ISP,LMNX,DENMAT1,LMRX,QLM(1,IAT))
        DEALLOCATE(DENMAT1)
      ENDDO
      CALL MPE$COMBINE('MONOMER','+',QLM)
!PI=4.D0*DATAN(1.D0)
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
      USE MPE_MODULE
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
      INTEGER(4)            :: NR
      REAL(8)               :: RBAS(3,3),GBAS(3,3),CELLVOL
      INTEGER(4)            :: IAT
      INTEGER(4)            :: ISP
      INTEGER(4)            :: LNX
      INTEGER(4),ALLOCATABLE:: LOX(:)
      REAL(8)   ,ALLOCATABLE:: DOVER(:,:)
      INTEGER(4)            :: L1,L2,LN1,LN2,LMN1,LMN2,M
      INTEGER(4)            :: NFILO
!     ***************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     ===========================================================================
!     ==  DETERMINE TOTAL MAGNETIZATION AS INTEGRATED MOMENT DENSITY           ==
!     ===========================================================================
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
!     ===========================================================================
!     ==  DETERMINE AUGMENTATION PART OF THE SPIN DENSITY                      ==
!     ===========================================================================
      AUGSPIN(:)=0.D0
      DO IAT=1,NAT
        IF(MOD(IAT-1,NTASKS)+1.NE.THISTASK)  CYCLE
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$LNX(ISP,LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$LOFLN(ISP,LNX,LOX)
        ALLOCATE(DOVER(LNX,LNX))
        CALL SETUP$1COVERLAP(ISP,LNX,DOVER)
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
                ELSE IF(NDIMD.EQ.3) THEN
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
      CALL MPE$COMBINE('MONOMER','+',AUGSPIN)
!
!     ===========================================================================
!     ==  REPORT TOTAL SPIN                                                    ==
!     ===========================================================================
      TOTSPIN=PSSPIN+AUGSPIN
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,FMT='("ABS.VALUE OF INTEGRATED SPIN MOMENT DENSITY",F10.5)') &
     &                 SQRT(SUM(TOTSPIN(:)**2))
      WRITE(NFILO,FMT='("DIRECTION OF INTEGRATED SPIN MOMENT DENSITY",3F10.5)') &
     &                 TOTSPIN/SQRT(SUM(TOTSPIN(:)**2))
      RETURN
      END

!
!     ..................................................................
      SUBROUTINE WAVES$SPHERE(LMNXX,NDIMD_,NAT,LMRXX,RHOB,DENMAT,EDENMAT &
     &                       ,VQLM,DH,DO,POTB)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES ONE-CENTER HAMILTONIAN AN OVERLAPMATRIX           **
!     **                                                              **
!     **  ONE-CENTER DENSITY MATRICES                                 **
!     **  SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                     **
!     **  SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)              **
!     **  NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z)**
!     **                                                              **
!     **  ON RETURN, POTB IS MINUS THE AVERAGED ELECTROSTATIC         **
!     **  POTENTIAL FROM THE AUGMENTATION,                            **
!     **  I.E. -(AEPOTHARTREE-(PSPOTHARTREE-VADD))                    **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
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
      REAL(8)   ,INTENT(OUT) :: DH(LMNXX,LMNXX,NDIMD_,NAT)
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
      REAL(8)   ,ALLOCATABLE :: DH1(:,:,:)  !(LMNX,LMNX,NDIMD)
      REAL(8)   ,ALLOCATABLE :: DOV1(:,:,:) !(LMNX,LMNX,NDIMD)
      REAL(8)                :: RBAS(3,3)
      REAL(8)                :: GBAS(3,3)
      REAL(8)                :: SVAR
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES$SPHERE')
                              CALL TIMING$CLOCKON('W:SPHERE')
      IF(NDIMD_.NE.NDIMD.OR.NAT.NE.MAP%NAT) THEN
        CALL ERROR$MSG('ARRAY SIZE INCONSISTENT')
        CALL ERROR$STOP('WAVES$DENMAT')
      END IF
!
!     ===================================================================
!     ==  EVALUATE ONE-CENTER HAMILTONIAN                              ==
!     ===================================================================
      DH(:,:,:,:)=0.D0
      DO(:,:,:,:)=0.D0
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      POTB=0
      DO IAT=THISTASK,NAT,NTASKS   ! DISTRIBUTE WORK ACCROSS TASKS
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        CALL SETUP$LMRX(ISP,LMRX)
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
!     ===================================================================
!     ==                                                               ==
!     ===================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL GBASS(RBAS,GBAS,SVAR)
      POTB=POTB/SVAR
      CALL AUGMENTATION$SYNC
                              CALL TIMING$CLOCKOFF('W:SPHERE')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$SPHERE
!
!     ..................................................................
      SUBROUTINE WAVES$ADDCONSTANTPOT(NRL,LMNXX,NDIMD_,NAT,POTB,RHO,DH,DO)
!     ******************************************************************
!     **                                                              **
!     **  ADDS A CONSTANT POTENTIAL TO ENSURE THAT THE AVERAGED       **
!     **  ELECTROSTATIC POTENTIAL IS ZERO.                            **
!     **                                                              **
!     **  POTB IS MINUS THE INTEGRATED ELECTROSTATIC POTENTIAL FROM   **
!     **  THE AUGMENTATION AVERAGED OVER THE UNIT CELL                **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NRL
      REAL(8)   ,INTENT(IN)   :: POTB
      INTEGER(4),INTENT(IN)   :: LMNXX
      INTEGER(4),INTENT(IN)   :: NDIMD_
      INTEGER(4),INTENT(IN)   :: NAT
      REAL(8)   ,INTENT(IN)   :: DO(LMNXX,LMNXX,NDIMD_,NAT)
      REAL(8)   ,INTENT(INOUT):: DH(LMNXX,LMNXX,NDIMD_,NAT)
      REAL(8)   ,INTENT(INOUT):: RHO(NRL,NDIMD_)
      INTEGER(4)              :: ISP
      INTEGER(4)              :: IAT
      INTEGER(4)              :: LMN1,LMN2
      INTEGER(4)              :: LN1,LN2
      INTEGER(4)              :: L1,L2,M
!     ******************************************************************
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
     &                                 +POTB*DO(LMN1+M,LMN2+M,1,IAT)
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
!     ...................................................................
      SUBROUTINE WAVES$FORCE(NAT,LMNXX,NDIMD_,DH,FORCE,STRESS)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES FORCES FROM THE AUGMENTATION PART                 **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NAT
      INTEGER(4),INTENT(IN)  :: LMNXX
      INTEGER(4),INTENT(IN)  :: NDIMD_
      REAL(8)   ,INTENT(IN)  :: DH(LMNXX,LMNXX,NDIMD,NAT)
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
      REAL(8)   ,ALLOCATABLE :: DH1(:,:,:)  ! (LMNX,LMNX,NDIM**2)    
      COMPLEX(8),ALLOCATABLE :: DEDPROJ(:,:,:) ! (NDIM,NBH,LMNX)
      COMPLEX(8),ALLOCATABLE :: DEDPRO(:,:) ! (NGL,LMNX)
      COMPLEX(8),ALLOCATABLE :: EIGR(:)     ! (NGL)
      REAL(8)                :: FORCE1(3)
      REAL(8)                :: STRESS1(3,3)
      REAL(8)                :: SVAR
      REAL(8)                :: R(3,NAT)
      REAL(8)   ,PARAMETER   :: TINY=1.D-300
!     ******************************************************************
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
!     ==================================================================
!     ==  GET OCCUPATIONS FROM DYNOCC OBJECT                          ==
!     ==================================================================
      CALL CELL$GETL4('MOVE',TSTRESS)
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
!
!     ===================================================================
!     ==                                                               ==
!     ===================================================================
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
!         ===============================================================
!         ==                                                           ==
!         ===============================================================
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
          ELSE
            GIJ(:,:)=0
          END IF
!
!         ===============================================================
!         ==                                                           ==
!         ===============================================================
          IPRO=1
          DO IAT=1,NAT
            ISP=MAP%ISP(IAT)
            IBPRO=1+SUM(MAP%LNX(1:ISP-1))
            LMNX=MAP%LMNX(ISP)
            LNX=MAP%LNX(ISP)
!           == DEDPROJ= (DH-LAMBDA*DO)<P|PSPSI> =========================
            ALLOCATE(DO1(LNX,LNX))
            CALL SETUP$1COVERLAP(ISP,LNX,DO1)
            ALLOCATE(DEDPROJ(NDIM,NBH,LMNX))
            ALLOCATE(DH1(LMNX,LMNX,NDIM**2))
            IF(NDIM.EQ.1) THEN
              DH1(:,:,1)=DH(1:LMNX,1:LMNX,ISPIN,IAT)
            ELSE 
              DH1(:,:,:)=DH(1:LMNX,1:LMNX,:,IAT)
            END IF
!           ==  DEDPROJ=DE/D<PSPSI|PRO>=DH*<PRO|PSPSI>
            CALL WAVES_DEDPROJ(NDIM,NBH,NB,LNX,MAP%LOX(1:LNX,ISP),LMNX &
     &                     ,OCC(:,IKPT,ISPIN) &
     &                     ,THIS%PROJ(:,:,IPRO:IPRO+LMNX-1),DH1,DO1 &
     &                     ,THIS%RLAM0,DEDPROJ)
            DEALLOCATE(DH1)
            DEALLOCATE(DO1)
!           == |DE/DPRO>=|PSPSI>DEDPROJ ==================================
            ALLOCATE(DEDPRO(NGL,LMNX))
            CALL WAVES_DEDPRO(GSET%TINV,NGL,NDIM,NBH,THIS%PSI0,LMNX,DEDPROJ,DEDPRO)
            DEALLOCATE(DEDPROJ)
!           == DE= <DPRO|DEDPRO> ========================================
            ALLOCATE(EIGR(NGL))
            CALL PLANEWAVE$STRUCTUREFACTOR(R(1,IAT),NGL,EIGR)
!           == F=2*RE[ <PSI|DPRO/DR>*DH*<PRO|PSI> ]
            CALL WAVES_PROFORCE(LNX,LMNX,MAP%LOX(1:LNX,ISP),NGL,GWEIGHT,GVEC,GIJ &
     &                 ,GSET%PRO(:,IBPRO:IBPRO+LNX-1),GSET%DPRO(:,IBPRO:IBPRO+LNX-1) &
     &                 ,MAP%LMX,GSET%YLM,GSET%SYLM &
     &                 ,EIGR,DEDPRO,FORCE1,TSTRESS,STRESS1)
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
!     ...................................................................
      SUBROUTINE WAVES$HPSI(NRL,NDIMD_,NAT,LMNXX,RHO,DH)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES FORCES FROM THE AUGMENTATION PART                 **
!     **                                                              **
!     ************P.E. BLOECHL, TU-CLAUSTHAL (2005)*********************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NRL
      INTEGER(4),INTENT(IN)  :: NDIMD_
      INTEGER(4),INTENT(IN)  :: LMNXX
      INTEGER(4),INTENT(IN)  :: NAT
      REAL(8)   ,INTENT(IN)  :: RHO(NRL,NDIMD_) ! PS-POTENTIAL
      REAL(8)   ,INTENT(IN)  :: DH(LMNXX,LMNXX,NDIMD_,NAT) ! ONE-CENTER HAMILTONIAN
      INTEGER(4)             :: IKPT,ISPIN
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NBH
      REAL(8)                :: R(3,NAT)
!     ------------------------------------------------------------------
      INTEGER(4)             :: IPRO
      INTEGER(4)             :: ISP
      INTEGER(4)             :: IAT
      INTEGER(4)             :: LMNX
      REAL(8)   ,ALLOCATABLE :: DH1(:,:,:)      ! 1CENTER HAMILTONIAN
      COMPLEX(8),ALLOCATABLE :: HPROJ(:,:,:)    ! DH*PROJ
!     ******************************************************************
                              CALL TRACE$PUSH('WAVES$HPSI')
                              CALL TIMING$CLOCKON('W:HPSI')
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
!         ==============================================================
!         ==  MULTIPLY WAVE FUNCTION WITH LOCAL POTENTIAL AND         ==
!         ==  ADD KINETIC ENERGY                                      ==
!         ==============================================================
          NGL=GSET%NGL
          NBH=THIS%NBH
          IF(.NOT.ASSOCIATED(THIS%HPSI))ALLOCATE(THIS%HPSI(NGL,NDIM,NBH))
!         == NOTE: THE ARRAY RHO CONTAINS THE POTENTIAL ================
CALL TIMING$CLOCKON('W:HPSI.VPSI')
!======================================================================== 
IF(1.EQ.1)THEN !OLD VERSION CHANGE WAS REQUIRED FOR KAESTNERS CONJUGATE GRADIENT
!======================================================================== 
         CALL WAVES_VPSI(GSET,NGL,NDIM,NBH,NRL,THIS%PSI0,RHO(1,ISPIN) &
     &                   ,THIS%HPSI)
CALL TIMING$CLOCKOFF('W:HPSI.VPSI')
!
!         ==============================================================
!         ==  EVALUATE  DH<P|PSI>                                     ==
!         ==============================================================
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
     &           ,DH1,THIS%PROJ(:,:,IPRO:IPRO+LMNX-1),HPROJ(:,:,IPRO:IPRO+LMNX-1))
            DEALLOCATE(DH1)
            IPRO=IPRO+LMNX
          ENDDO
CALL TIMING$CLOCKOFF('W:HPSI.HPROJ')
!
!         ==============================================================
!         ==  ADD  |P>DH<P|PSI>                                       ==
!         ==============================================================
CALL TIMING$CLOCKON('W:HPSI.ADDPROJ')
          CALL WAVES_ADDPRO(MAP,GSET,NAT,R,NGL,NDIM,NBH,MAP%NPRO &
     &                     ,THIS%HPSI,HPROJ)
CALL TIMING$CLOCKOFF('W:HPSI.ADDPROJ')
          DEALLOCATE(HPROJ)
!======================================================================== 
ELSE
!======================================================================== 
         CALL WAVES_HPSI(MAP,GSET,ISPIN,NGL,NDIM,NDIMD,NBH,MAP%NPRO,LMNXX,NAT,NRL&
     &                  ,THIS%PSI0,RHO(1,ISPIN),R,THIS%PROJ,DH,THIS%HPSI)
!======================================================================== 
END IF
!======================================================================== 
        ENDDO
      ENDDO
                              CALL TIMING$CLOCKOFF('W:HPSI')
                              CALL TRACE$POP
      RETURN
      END SUBROUTINE WAVES$HPSI
!
!     .................................................................
      SUBROUTINE WAVES_HPROJ(NDIM,NB,LMNX,DH,PROJ,HPROJ)
!     *****************************************************************
!     **                                                             **
!     *******************************************P.E. BLOECHL, (1991***
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
!     *****************************************************************
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
!     ....................................................................
      SUBROUTINE WAVES_HPSI(MAP,GSET,ISPIN,NGL,NDIM,NDIMD,NBH,NPRO,LMNXX,NAT,NRL &
     &                     ,PSI,POT,R,PROJ,DH,HPSI)
      USE WAVES_MODULE, ONLY: MAP_TYPE,GSET_TYPE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: ISPIN
      TYPE(MAP_TYPE) ,INTENT(IN) :: MAP  !NAT,NSP,NBAREPRO,NRL,NPRO,LMX,LNXX,LNX,LMNX,LOX,ISP
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
      COMPLEX(8)     ,INTENT(IN) :: PROJ(NDIM,NBH,NPRO)     !(NDIM,NBH,NPRO) <PSPSI|P>
      REAL(8)        ,INTENT(IN) :: DH(LMNXX,LMNXX,NDIMD,NAT)
      COMPLEX(8)     ,INTENT(OUT):: HPSI(NGL,NDIM,NBH)
      COMPLEX(8)     ,ALLOCATABLE:: HPROJ(:,:,:)  
      REAL(8)        ,ALLOCATABLE:: DH1(:,:,:)
      INTEGER(4)                 :: IPRO,IAT,ISP,LMNX
!     ****************************************************************
CALL TIMING$CLOCKON('W:HPSI.VPSI')
      CALL WAVES_VPSI(GSET,NGL,NDIM,NBH,NRL,PSI,POT,HPSI)
CALL TIMING$CLOCKOFF('W:HPSI.VPSI')
!
!     ================================================================
!     ==  EVALUATE  DH<P|PSI>                                       ==
!     ================================================================
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
!
!     ==============================================================
!     ==  ADD  |P>DH<P|PSI>                                       ==
!     ==============================================================
CALL TIMING$CLOCKON('W:HPSI.ADDPRO')
      CALL WAVES_ADDPRO(MAP,GSET,NAT,R,NGL,NDIM,NBH,NPRO,HPSI,HPROJ)
CALL TIMING$CLOCKOFF('W:HPSI.ADDPRO')
      DEALLOCATE(HPROJ)
      RETURN
      END
!
!     .................................................................
      SUBROUTINE WAVES_OPSI(NB,NBH,NPRO,NAT,NGL,R0,PROJ,OPSI)
!     *****************************************************************
!     **                                                             **
!     **  |PSI>+|P>DO<P|PSI>                                         **
!     **                                                             **
!     *********************************JOHANNES KAESTNER 2004**********
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)   :: NB
      INTEGER(4) ,INTENT(IN)   :: NBH
      INTEGER(4) ,INTENT(IN)   :: NPRO
      INTEGER(4) ,INTENT(IN)   :: NAT
      INTEGER(4) ,INTENT(IN)   :: NGL
      REAL(8)    ,INTENT(IN)   :: R0(3,NAT)
      COMPLEX(8) ,INTENT(IN)   :: PROJ(NDIM,NBH,NPRO)
      COMPLEX(8) ,INTENT(INOUT):: OPSI(NGL,NDIM,NBH) !PSI ON ETRY, OPSI ON EXIT
      INTEGER(4)               :: IPRO,IAT,ISP,LNX,LMNX
      COMPLEX(8) ,ALLOCATABLE  :: OPROJ(:,:,:)
      REAL(8)    ,ALLOCATABLE  :: DO(:,:)      ! 1-CENTER OVERLAP
!     *****************************************************************
!
!     ==============================================================
!     ==  EVALUATE  DO<P|PSI>  ASSUMING <PRO|PSI> IS STILL VALID  ==
!     ==============================================================
      ALLOCATE(OPROJ(NDIM,NBH,NPRO))
      OPROJ(:,:,:)=(0.D0,0.D0)
      IPRO=1
      DO IAT=1,NAT
         ISP=MAP%ISP(IAT)
         LNX=MAP%LNX(ISP)
         LMNX=MAP%LMNX(ISP)
         ALLOCATE(DO(LNX,LNX))
         CALL SETUP$1COVERLAP(ISP,LNX,DO)
         CALL WAVES_OPROJ(LNX,MAP%LOX(1:LNX,ISP),DO,NDIM,LMNX,NBH &
              &          ,PROJ(:,:,IPRO:IPRO+LMNX-1),OPROJ(:,:,IPRO:IPRO+LMNX-1))
         DEALLOCATE(DO)
         IPRO=IPRO+LMNX
      ENDDO
!
!     ==============================================================
!     ==  ADD  |PSI>+|P>DO<P|PSI>                                 ==
!     ==============================================================
      !ALLOCATE(THIS%OPSI(NGL,NDIM,NBH))
      !THIS%OPSI(:,:,:)=THIS%PSI0(:,:,:)
      CALL WAVES_ADDPRO(MAP,GSET,NAT,R0,NGL,NDIM,NBH,NPRO,OPSI,OPROJ)
      DEALLOCATE(OPROJ)
      ! THIS THIS%OPSI IS MY O|PSI> FOR CG
      RETURN
      END 
!
!     .................................................................
      SUBROUTINE WAVES_VPSI(GSET,NGL,NDIM,NBH,NRL,PSI,V,HPSI)
!     *****************************************************************
!     **                                                             **
!     **  EVALUATES H*PSI WITHOUT THE AUGMENTATION PART              **
!     **  |HPSI>=(-0.5*NABLA**2+PS-V)|PSPSI(0)>                      **
!     **                                                             **
!     *******************************************P.E. BLOECHL, (1991)**
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
      COMPLEX(8),ALLOCATABLE :: PSIOFR(:,:,:)
      COMPLEX(8)             :: PSIUP,PSIDN
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
!     *******************************************P.E. BLOECHL, (1999)***
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
      INTEGER(4)            :: IB,IG,IDIM1,IDIM2
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
    &              +SVAR*REAL(CONJG(PSI(IG,IDIM1,IB))*PSI(IG,IDIM2,IB),KIND=8)
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
     &           +SVAR*REAL(CONJG(PSI(IG,IDIM1,IB))*PSI1(IG,IDIM2),KIND=8)
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
        CALL MPE$COMBINE('NONE','+',EKIN)
        CALL MPE$COMBINE('NONE','+',STRESS)
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
        CALL MPE$COMBINE('NONE','+',EKIN)
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
!     *******************************************P.E. BLOECHL, (1991)***
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
      INTEGER(4)             :: IBH,IR
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
          CALL ERROR$STOP('WAVES_DENSITY')
        END IF
      ELSE 
        IF(NBH.NE.NB) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF BANDS')
          CALL ERROR$STOP('WAVES_DENSITY')
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
            RHO(IR,1)=RHO(IR,1)+REAL(PSIOFR(IR,1,IBH)*PSI1(IR),KIND=8)
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
              RE= REAL(PSIOFR(IR,1,IBH),KIND=8)
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
    &               +F1*REAL(PSIOFR(IR,1,IBH)*CONJG(PSIOFR(IR,1,IBH)),KIND=8)
              RHO(IR,4)=RHO(IR,4) &
    &               +F1*REAL(PSIOFR(IR,2,IBH)*CONJG(PSIOFR(IR,2,IBH)),KIND=8)
              CSVAR=PSIOFR(IR,1,IBH)*CONJG(PSIOFR(IR,2,IBH))
              RHO(IR,2)=RHO(IR,2)+F1*REAL(CSVAR,KIND=8)
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
!     **            = DH<P|PSITILDE>F-DO<P|PSITILDE>*LAMBDA           **
!     **  LAMBDA IS CALCULATED FROM THIS%RLAM0 BY MULTIPLICATION WITH **
!     **  0.5*(FI+FJ)                                                 **
!     **                                                              **
!     *******************************************P.E. BLOECHL, (1999)***
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
                               CALL TRACE$PUSH('WAVES_DEDPROJ')
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
      CALL WAVES_HPROJ(NDIM,NBH,LMNX,DH,PROJ,DEDPROJ)
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
                               CALL TRACE$Pop
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
!     *******************************************P.E. BLOECHL, (1999)***
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
      COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
      COMPLEX(8),ALLOCATABLE  :: PSIM(:)
      LOGICAL(4),PARAMETER    :: TESSL=.TRUE.
      COMPLEX(8)              :: DEDPROJ1(NDIM,NBH,LMNX) ! CONJG(<P|PSI>)
!     ******************************************************************
!
!     ==================================================================
!     ==  SUPERPOSE WAVE FUNCTION TO OBTAIN DE/D<P|                   ==
!     ==================================================================
      DEDPROJ1=CONJG(DEDPROJ)
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
!     ..................................................................
      SUBROUTINE WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NB,NPRO,PSI,PROPSI)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE PROJECTIONS                                       **
!     **                                                              **
!     **                                                              **
!     *******************************************P.E. BLOECHL, (1998)***
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
     &                      ,GSET%PRO(:,IBPRO:IBPRO+LNX-1),MAP%LMX,GSET%YLM,EIGR,PRO)
        CALL PLANEWAVE$SCALARPRODUCT(' ',NGL,1,LMNX,PRO,NDIM*NB,PSI &
     &                              ,PROPSI1)
        II=0
        DO IB=1,NB
          DO IDIM=1,NDIM
            DO LMN=1,LMNX
              II=II+1    ! II=(LMN,IDIM,IB)
              PROPSI(IDIM,IB,IPRO-1+LMN)=PROPSI1(II)
!PRINT*,'IAT,LMN,PROPSI',IAT,LMN,IB,PROPSI1(II)
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
!DO IB=1,NB
!  DO LMN=1,LMNX
!     PRINT*,'PROPSI ',IB,LMN,PROPSI(1,IB,LMN)
!  ENDDO
!ENDDO
!STOP
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
!     *******************************************P.E. BLOECHL, (1998)***
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
      LOGICAL(4)    ,PARAMETER   :: TESSL=.TRUE.
      INTEGER(4)                 :: IPRO,IBPRO,IAT,LMN,IB,IDIM,II
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
     &                      ,GSET%PRO(:,IBPRO:IBPRO+LNX-1),MAP%LMX,GSET%YLM,EIGR,PRO)
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
!       == PSI=PSI+PRO*PROPSI1
        CALL LIB$ADDPRODUCTC8(.FALSE.,NGL,LMNX,NDIM*NB,PRO,PROPSI1,PSI)
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
!     *******************************************P.E. BLOECHL, (1998)***
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
!     *******************************************P.E. BLOECHL, (1998)***
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
!     *******************************************P.E. BLOECHL, (1999)***
      USE WAVES_MODULE
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
!     ==================================================================
!     ==  UPDATE RADIAL PROJECTOR FUNCTIONS                           ==
!     ==================================================================
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
      REAL(8)               :: SVAR1,SVAR2,SVAR3
!     ******************************************************************
      IF(OPTIMIZERTYPE.EQ.'CG') RETURN                           !KAESTNERCG
                              CALL TRACE$PUSH('WAVES$PROPAGATE')
!
!     ==================================================================
!     ==  STOP WAVE FUNCTIONS                                         ==
!     ==================================================================
      IF(TSTOP) THEN
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            THIS%PSIM(:,:,:)=THIS%PSI0(:,:,:)
          ENDDO
        ENDDO
        WAVEEKIN1=0.D0
        TSTOP=.FALSE.
      END IF
!
!     ==================================================================
!     ==  PROPAGATE WAVE FUNCTIONS (PUT PSI(+) INTO PSIM)             ==
!     ==================================================================
      DO IKPT=1,NKPTL
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
      DO IKPT=1,NKPTL
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
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      ALLOCATE(EIG(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
      EIG(:,:,:)=0.D0
!
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
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
!           ============================================================
!           ==  WAVE FUNCTION KINETIC ENERGY FOR GENERAL WAVE FUNCTIONS=
!           ==  AND FIRST PART FOR SUPER WAVE FUNCTIONS                =
!           ============================================================
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
!     ==================================================================
!     == COMMUNICATE                                                  ==
!     ==================================================================
      CALL MPE$COMBINE('MONOMER','+',EKIN)
      CALL MPE$COMBINE('K','+',EIG)
!
!     ==================================================================
!     == SET CONTRIBUTION TO EIGENVALUES                              ==
!     ==================================================================
      CALL WAVES_DYNOCCSETR8A('M<PSIDOT|PSIDOT>',NBX*NKPTL*NSPIN,EIG)
      DEALLOCATE(EIG)
      DEALLOCATE(OCC)
      RETURN
      END
! 
!      ...............................................................
       SUBROUTINE WAVES_KDISTRIBUTE(NTASKS,NKPT,TINV,ICOLOR,KMAP)
!      **                                                           **
!      **   ICOLOR GIVES EACH TASK AN INTEGER. TASKS WITH EQUAL     **
!      **   VALUE OF ICOLOR WILL BE IN THE SAME GROUP               **
!      **                                                           **
!      **   KMAP POINTS TO THE TASK, WHICH IS THE FIRST TASK IN THE **
!      **   CORRESPONDING K-GROUP                                   **
!      **                                                           **
!      ********************* PETER BLOECHL, TU CLAUSTJAL 2005 *********
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
       INTEGER(4)            :: I,J,IKPT,IPOS,IGROUP
       INTEGER(4)            :: ISVAR
       INTEGER(4)            :: IND(1)
       INTEGER(4)            :: IFIRST(NTASKS)
       INTEGER(4),ALLOCATABLE:: NA(:),NB(:)
       INTEGER(4),ALLOCATABLE:: MG(:),MGSAVE(:)
       INTEGER(4)            :: NGROUPS
!      ***************************************************************
       DO I=1,NKPT
         WK(I)=WA
         WK0(I)=WA0
         IF(TINV(I))WK(I)=WB
         IF(TINV(I))WK0(I)=WB0
       ENDDO
!
!      ==============================================================
!      ==  OPTMIZATION OF SK. SK US A POINTER IKPT->IGROUP         ==
!      ==============================================================
       XLOAD1=HUGE(XLOAD1)
!      == TEST DIFFERENT NUMBER OF SUBGROUPS ========================
       DO NGROUPS=1,NTASKS       
         IF(NGROUPS.GT.NKPT) EXIT  ! THERE MUST NOT BE AN EMPTY GROUP
!        ============================================================
!        == DISTRIBUTE ALL K-POINTS INTO "NGROUPS" GROUPS          ==
!        == NB(IGROUP) IS THE NUMBER OF KPOINTS WITH REAL WAVE     ==
!        == FUNCTIONS IN THE GROUP "IGROUP"; NA(IGROUP IS THE      == 
!        == NUMBER OF COMPLEX WAVE FUNCTIONS                       ==
!        ============================================================
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
!        === RESHUFFLE NA AND NB TO DISTRIBUTE LOAD EQUALLY AMONG GROUPS
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
!        == CONSISTENCY CHECK =======================================
         DO IGROUP=1,NGROUPS
           IF(NB(IGROUP).NE.0.OR.NA(IGROUP).NE.0) THEN
             CALL ERROR$MSG('INCONSISTENCY DETECTED WITH NA,NB')
             CALL ERROR$STOP('WAVES_KDISTRIBUTE')
           END IF
         ENDDO
         DEALLOCATE(NA)
         DEALLOCATE(NB)
!
!        =============================================================
!        == CHECK IF THIS #(GROUPS) IMPROVED THE LOADBALANCE        ==
!        == AND DISCARD IF NOT.                                     ==
!        =============================================================
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
!      ==============================================================
!      ==  WRAP UP  AND PREPARE OUTPUT                             ==
!      ==============================================================
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
!      ==============================================================
!      == CHECKS ====================================================
!      ==============================================================
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

       RETURN
       END
!
!      ...............................................................
       SUBROUTINE WAVES_LOADPERPROC(NGROUPS,NTASKS,NKPT,WK,WK0,SK,MG,XLOAD)
!      **  ASSIGN TO EACH GROUP A NUMBER OF PROCESSORS              **
!      **                                                           **
!      **                                                           **
!      ********************* PETER BLOECHL, TU CLAUSTJAL 2005 *********
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
!      ***************************************************************
!
!      ===============================================================
!      == START UP                                                  ==
!      ===============================================================
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
!      ===============================================================
!      == ENFORCE SUM RULE                                          ==
!      ===============================================================
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
!      ===============================================================
!      == RESHUFFLE NODES TO REDUCE MAXIMUM LOAD                    ==
!      ===============================================================
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
           XLOADTEST(I)=HUGE(XLOADTEST)
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
!      ===============================================================
!      == DETERMINE LOAD FOR EACH PROCESSOR                         ==
!      ===============================================================
       XLOAD(:)=WG0(:)+WG(:)/MG(:)
       RETURN
       END
!
!      ...............................................................
       SUBROUTINE WAVES_PREBALANCE(NGROUPS,NA,NB,WA,WB)
!      **                                                             **
!      **  TRIES DO DISTRIBUTE THE COMPUTATIONAL EFFORT EQUALLY       **
!      **  ONTO NGROUPS GROUPS. IT ASSUMES THAT THERE ARE TWO TYPES   **
!      **  OF KPOINTS WITH WEIGHTS WA AND WB. NA AND NB SPECIFIES     **
!      **  THE NUMBER OF K-POINTS FROM EACH TYPE IN THE CORRESPONDING **
!      **  GROUPS                                                     **
!      **                                                             **
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
!      ***************************************************************
       ICOUNT=0
       DO ICOUNT=1,100
!        == XLOAD QUANTIFIES THE ESTIMATED COMPUTATIONAL EFFORT FOR EACH GROUP
         XLOAD(:)=REAL(NA(:),KIND=8)*WA+REAL(NB(:),KIND=8)*NB
!        == FIND THE GROUP HAVING THE LARGEST LOAD
         IND=MAXLOC(XLOAD)
         IPOS1=IND(1)
         XLOAD1=XLOAD(IPOS1)
!        == IT IS NOT ALLOWED TO REDUCE THE NUMBER OF KPOINTS IN A GROUP TO ZERO
         IF(NA(IPOS1)+NB(IPOS1).EQ.1) EXIT
!        ==  TRY TO RESHUFFLE FROM THE GROUP WITH THE LARGEST LOAD TO THE OTHERS
!        ==  XLOADTESTA ESTIMATES THE LOAD WHEN A K-POINT OF TYPE A IS ADDED
!        ==  XLOADTESTB ESTIMATES THE LOAD WHEN A K-POINT OF TYPE B IS ADDED
         XLOADTESTA(IPOS1)=HUGE(XLOADTESTA)
         XLOADTESTB(IPOS1)=HUGE(XLOADTESTA)
         DO J=1,NGROUPS
           IF(J.EQ.IPOS1) CYCLE
           XLOADTESTA(J)=XLOAD(J)+WA
           XLOADTESTB(J)=XLOAD(J)+WB
         ENDDO
         IF(NA(IPOS1).EQ.0) XLOADTESTA=HUGE(XLOADTESTA)
         IF(NB(IPOS1).EQ.0) XLOADTESTB=HUGE(XLOADTESTB)
         IND=MINLOC(XLOADTESTA)
         IPOS2A=IND(1)       
         XLOAD2A=XLOADTESTA(IPOS2A)
         IND=MINLOC(XLOADTESTB)
         IPOS2B=IND(1)
         XLOAD2B=XLOADTESTB(IPOS2B)
         IF(MIN(XLOAD2A,XLOAD2B).LT.XLOAD1) THEN 
          IF(XLOAD2A.LE.XLOAD1) THEN
             NA(IPOS1)=NA(IPOS1)-1        
             NA(IPOS2A)=NA(IPOS2A)+1        
           ELSE
             NB(IPOS1)=NB(IPOS1)-1        
             NB(IPOS2A)=NB(IPOS2A)+1        
           END IF
           CYCLE
         ELSE
           EXIT
         END IF
         IF(ICOUNT.EQ.100) THEN
           CALL ERROR$MSG('LOOP NOT CONVERGED')
           CALL ERROR$STOP('WAVES_PREBALANCE')
         END IF
       ENDDO
       RETURN
       END
