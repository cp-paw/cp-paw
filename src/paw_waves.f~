!**********************************************************************
!***********************************************************************
!**                                                                   **
!**  WAVES OBJECT                                                     **
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
!**                                                                   **
!**                                                                   **
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
  INTEGER(4)         :: NAT         ! #(ATO
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
!     ==  energy expectations values of the one-partcle states      ==
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
        VAL(:)=THIS%EXPECTVAL(:)
!
!     ================================================================
!     ==  ENERGY EIGENVALUES                                        ==
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
        VAL(:)=THIS%EIGVAL(:)
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
        CALL SETUP$LNX(ISP,MAP%LNX(ISP))
        CALL SETUP$LOFLN(ISP,MAP%LNX(ISP),MAP%LOX(1,ISP))
        MAP%LMNX(ISP)=0
        DO LN=1,MAP%LNX(ISP)
          MAP%LMNX(ISP)=MAP%LMNX(ISP)+2*MAP%LOX(LN,ISP)+1
          MAP%LMX=MAX(MAP%LMX,(MAP%LOX(LN,ISP)+1)**2)
        ENDDO
        MAP%NBAREPRO=MAP%NBAREPRO+MAP%LNX(ISP)
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
      CALL OPTICS$GETL4('ON',TCHK)
      IF(TCHK) THEN
        ALLOCATE(XK(3,NKPT))
        CALL DYNOCC$GETR8A('XK',3*NKPT,XK)
        TCHK=.FALSE.
        DO IKPT=1,NKPT
          TCHK=XK(1,IKPT).EQ.0.AND.XK(2,IKPT).EQ.0.AND.XK(3,IKPT).EQ.0
          IF(TCHK) THEN
            CALL WAVES_SELECTWV(IKPT,1)
            CALL PLANEWAVE$SELECT(GSET%ID)
            CALL OPTICS$LATTICE
            EXIT
          END IF
        ENDDO
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('OPTICS CODE REQUIRES GAMMA POINT IN THE K-POINT SET')
          CALL ERROR$STOP('WAVES$INITIALIZE')
        END IF
        DEALLOCATE(XK)
      END IF
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
      REAL(8)   ,ALLOCATABLE :: GIJ(:,:)   !GI*GJ/G2
      LOGICAL(4)             :: TINV
      LOGICAL(4)             :: TSTRESS
      LOGICAL(4)             :: TFORCE
      LOGICAL(4)             :: TCHK
      LOGICAL(4)             :: TSO=.FALSE. ! EVALUATE DENMATI
      REAL(8)                :: PI
      COMPLEX(8),ALLOCATABLE :: QMAT(:,:)   
      INTEGER(4)             :: NFILO
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
call waves_comparepsi('before gramschmidt',ngl,ndim,nbh,this%psi0,THIS%psim)
            CALL WAVES_GRAMSCHMIDT(MAP,GSET,NAT,R,NGL,NDIM,NBH,NB,THIS%PSI0)
            CALL WAVES_GRAMSCHMIDT(MAP,GSET,NAT,R,NGL,NDIM,NBH,NB,THIS%PSIM)
call waves_comparepsi('after gramschmidt',ngl,ndim,nbh,this%psi0,THIS%psim)
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
!print*,ikpt,gset%tinv
!print*,'occ',occ(:,ikpt,ispin)
          CALL WAVES_EKIN(NGL,NDIM,NBH,NB,OCC(1,IKPT,ISPIN),GWEIGHT &
      &                         ,THIS%PSI0,EKIN1,TSTRESS,STRESS1 &
      &                         ,TBUCKET,GSET%BUCKET,GSET%DBUCKET)
!        CALL WAVES_EKIN_OLD(GSET%NGL,NDIM,THIS%NBH,THIS%NB,OCC(1,IKPT,ISPIN) &
!     &                   ,THIS%PSI0,EKIN1,TSTRESS,STRESS1)
          EKIN  =EKIN  +EKIN1
          STRESS=STRESS+STRESS1
        ENDDO
      ENDDO
!stop
PRINT*,'EKIN ',EKIN
      CALL ENERGYLIST$SET('PS  KINETIC',EKIN)
      CALL ENERGYLIST$ADD('AE  KINETIC',EKIN)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EKIN)
WRITE(*,FMT='("KIN STRESS ",3F15.7)')STRESS(1,:)
WRITE(*,FMT='("KIN STRESS ",3F15.7)')STRESS(2,:)
WRITE(*,FMT='("KIN STRESS ",3F15.7)')STRESS(3,:)
!
!     ==================================================================
!     == ONE-CENTER DENSITY MATRICES                                  ==
!     == SPIN RESTRICTED NSPIN=1;NDIM=1: (TOTAL)                      ==
!     == SPIN POLARIZED  NSPIN=2;NDIM=1: (TOTAL,SPIN_Z)               ==
!     == NONCOLLINEAR    NSPIN=1;NDIM=2: (TOTAL,SPIN_X,SPIN_Y,SPIN_Z) ==
!     ==================================================================
CALL TRACE$PASS('BEFORE ONE-CENTER DENSITY MATRIX')     
CALL TIMING$CLOCKOFF('W:EKIN')
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
      CALL MPE$COMBINE('+',potb)
!     ==  SPREAD INFO FROM SPHERE OVER ALL TASKS AND UPDATE ENERGYLIST
      CALL AUGMENTATION$SYNC
!
!     ==  subtract average electrostatic potential =====================
!     ==  to account for background density
      call gbass(rbas,gbas,svar)
      potb=potb/svar
      rho(:,:)=rho(:,:)+potb
CALL TIMING$CLOCKOFF('W:SPHERE')
!
!     ==================================================================
!     == COMMUNICATE DATA WITH OPTICS MODULE                          ==
!     ==================================================================
      CALL OPTICS$GETL4('ON',TCHK)
      IF(TCHK) THEN
        CALL MPE$COMBINE('+',DO) ! DO CURRENTLY USED ONLY FOR OPTICS
        IF(NDIM.EQ.2) THEN
          CALL ERROR$MSG('OPTICS CODE DOES NOT WORK WITH NONCOLLINEAR MODE')
          CALL ERROR$STOP('WAVES$ETOT') 
        END IF
        CALL OPTICS$WRITEOUT_PARTIALS(NSPIN,NAT,LMNXX,DH,DO)
      END IF
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
            IF(TSTRESS) THEN
              ALLOCATE(GIJ(6,NGL))
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
              ALLOCATE(DH1(LMNX,LMNX,NDIMD))
              DH1(:,:,:)=DH(1:LMNX,1:LMNX,:,IAT)
              CALL WAVES_DEDPROJ(NDIM,NDIMD,NBH,NB,LNX,MAP%LOX(1,ISP),LMNX &
     &                       ,OCC(1,IKPT,ISPIN) &
     &                       ,THIS%PROJ(1,1,IPRO),DH1,DO1 &
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
              CALL WAVES_PROFORCE(LNX,LMNX,MAP%LOX(1,ISP),NGL,GWEIGHT,GVEC,GIJ &
     &                   ,GSET%PRO(1,IBPRO),GSET%DPRO(1,IBPRO) &
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
            IF(TSTRESS) DEALLOCATE(GIJ)
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
            CALL WAVES_HPROJ(LMNXX,NDIMD,DH(1,1,ISPIN,IAT) &
      &              ,NDIM,NBH,LMNX,THIS%PROJ(1,1,IPRO),HPROJ(1,1,IPRO))
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
     &                    ,THIS%PSI0(1,1,IB),THIS%HPSI(1,1,IB),HAMILTON)
              EIG(2*IB-1,IKPT,ISPIN)=REAL(HAMILTON(1,1))
              EIG(2*IB  ,IKPT,ISPIN)=REAL(HAMILTON(2,2))
            ELSE 
              CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,1,1 &
     &                    ,THIS%PSI0(1,1,IB),THIS%HPSI(1,1,IB),HAMILTON)
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
WRITE(*,FMT='("RBAS ",3F15.7)')RBAS(1,:)
WRITE(*,FMT='("RBAS ",3F15.7)')RBAS(2,:)
WRITE(*,FMT='("RBAS ",3F15.7)')RBAS(3,:)
WRITE(*,FMT='("AFTER WAVES STRESS ",3F15.7)')STRESS(1,:)
WRITE(*,FMT='("AFTER WAVES STRESS ",3F15.7)')STRESS(2,:)
WRITE(*,FMT='("AFTER WAVES STRESS ",3F15.7)')STRESS(3,:)
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
      SUBROUTINE WAVES_HPROJ(LMNXX,NDIMD,DH,NDIM,NB,LMNX,PROJ,HPROJ)
!     ******************************************************************
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      USE WAVES_MODULE, ONLY : MAP_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: LMNXX
      INTEGER(4),INTENT(IN)  :: NDIMD
      REAL(8)   ,INTENT(IN)  :: DH(LMNXX,LMNXX,NDIMD)
      INTEGER(4),INTENT(IN)  :: NDIM
      INTEGER(4),INTENT(IN)  :: NB
      INTEGER(4),INTENT(IN)  :: LMNX
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
      SUBROUTINE WAVES_DEDPROJ(NDIM,NDIMD,NBH,NB,LNX,LOX,LMNX,OCC &
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
      INTEGER(4),INTENT(IN)   :: NDIMD     ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN)   :: NB        ! #(BANDS)
      INTEGER(4),INTENT(IN)   :: NBH       ! #(WAVE FUNCTIONS)
      INTEGER(4),INTENT(IN)   :: LMNX      ! #(PROJECTORS ON THIS SITE)
      REAL(8)   ,INTENT(IN)   :: OCC(NB)   ! OCCUPATIONS
      COMPLEX(8),INTENT(IN)   :: PROJ(NDIM,NBH,LMNX) ! <P|PSI>
      REAL(8)   ,INTENT(IN)   :: DH(LMNX,LMNX,NDIMD) ! DE/DD
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
      CALL WAVES_HPROJ(LMNX,NDIMD,DH,NDIM,NBH,LMNX,PROJ,DEDPROJ)
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
     &                      ,GSET%PRO(1,IBPRO),MAP%LMX,GSET%YLM,EIGR,PRO)
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
     &                      ,GSET%PRO(1,IBPRO),MAP%LMX,GSET%YLM,EIGR,PRO)
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
            CALL SETUP$GETFOFG('PRO',.FALSE.,LN,NGL,G2,CELLVOL,GSET%PRO(1,IND))
            CALL SETUP$GETFOFG('PRO',.TRUE.,LN,NGL,G2,CELLVOL,GSET%DPRO(1,IND))
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
!
!     ..................................................................
      SUBROUTINE WAVES$ORTHOGONALIZE()
!     ******************************************************************
!     ** ENFORCES THE ORTHONORMALITY CONDITION OF THE WAVE FUNCTIONS  **
!     ******************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      COMPLEX(8),ALLOCATABLE :: OPROJ(:,:,:)
      REAL(8)   ,ALLOCATABLE :: MARR(:)
      REAL(8)   ,ALLOCATABLE :: R0(:,:)      !CURRENT ATOMIC POSITIONS
      REAL(8)   ,ALLOCATABLE :: RP(:,:)      ! NEXT ATOMIC POSITIONS
      REAL(8)   ,ALLOCATABLE :: RM(:,:)      ! LAST ATOMIC POSITIONS
      COMPLEX(8)             :: CSUM
      COMPLEX(8),ALLOCATABLE :: MAT(:,:)
      COMPLEX(8),ALLOCATABLE :: OMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: OOMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: AUXMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: LAMBDA(:,:)
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DO(:,:)      ! 1-CENTER OVERLAP
      REAL(8)   ,ALLOCATABLE :: G2(:)        ! SQUARE OF REC. LATTICE VECTORS
      REAL(8)   ,ALLOCATABLE :: GVEC(:,:)    ! RECIPROCAL LATTICE VECTORS
      REAL(8)                :: OCCI,OCCJ
      REAL(8)                :: FAC
      INTEGER(4)             :: NBX
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: IKPT,ISPIN,IPRO,IAT,ISP,M
      INTEGER(4)             :: LN1,L1,M1,LMN1,IPRO1
      INTEGER(4)             :: LN2,L2,M2,LMN2,IPRO2
      INTEGER(4)             :: IBH,IB,IDIM,IG,I,J
      INTEGER(4)             :: NGL,NBH,NB,LMNX,LNX,LMX,LN
      INTEGER(4)             :: NAT,IND
      REAL(8)   ,PARAMETER   :: DSMALL=1.D-12
      REAL(8)                :: SVAR,DOVER1
      REAL(8)                :: RBAS(3,3)      ! UNIT CELL
      REAL(8)                :: GBAS(3,3)      ! RECIPROCAL UNIT CELL
      REAL(8)                :: CELLVOL        ! UNIT CELL VOLUME
      REAL(8)                :: MAPTOCELL(3,3) ! R(+)=MAPTOCELL*X(+)
      REAL(8)                :: CELLSCALE      ! PSI(+)=CELLSCALE*PSI(+)
      REAL(8)   ,ALLOCATABLE :: YLM(:)
      LOGICAL(4)             :: TSTRESS
      LOGICAL(4)             :: TINV
      LOGICAL(4),PARAMETER   :: TTEST=.FALSE.
      COMPLEX(8)             :: CSVAR
      REAL(8)   ,ALLOCATABLE :: NORM(:)
      REAL(8)   ,ALLOCATABLE :: RMAT(:,:),ROMAT(:,:),ROOMAT(:,:),RLAMBDA(:,:)
      INTEGER(4),ALLOCATABLE :: SMAP(:)
      INTEGER(4)             :: I1,J1,K
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$ORTHOGONALIZE')
                             CALL TIMING$CLOCKON('WAVES$ORTHOGONALIZE')
      NPRO=MAP%NPRO
      NAT=MAP%NAT
      CALL CELL$GETL4('MOVE',TSTRESS)
!
!     ==================================================================
!     == COLLECT OCCUPATIONS                                          ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NBX*NKPT*NSPIN,OCC)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==================================================================
!     ==  CALCULATE FORCE OF CONSTRAINT                               ==
!     ==================================================================
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          NB=THIS%NB
!
!         ==============================================================
!         ==  EVALUATE  DO<P|PSI>  ASSUMING <PRO|PSI> IS STILL VALID  ==
!         ==============================================================
          ALLOCATE(OPROJ(NDIM,NBH,NPRO))
          OPROJ(:,:,:)=(0.D0,0.D0)
          IPRO=1
          DO IAT=1,NAT
            ISP=MAP%ISP(IAT)
            LNX=MAP%LNX(ISP)
            LMNX=MAP%LMNX(ISP)
            ALLOCATE(DO(LNX,LNX))
            CALL SETUP$1COVERLAP(ISP,LNX,DO)
            CALL WAVES_OPROJ(LNX,MAP%LOX(1,ISP),DO,NDIM,LMNX,NBH &
      &                      ,THIS%PROJ(1,1,IPRO),OPROJ(1,1,IPRO))
            DEALLOCATE(DO)
            IPRO=IPRO+LMNX
          ENDDO
!
!         ==============================================================
!         ==  ADD  |PSI>+|P>DO<P|PSI>                                 ==
!         ==============================================================
          ALLOCATE(THIS%OPSI(NGL,NDIM,NBH))
          THIS%OPSI(:,:,:)=THIS%PSI0(:,:,:)
          CALL WAVES_ADDPRO(MAP,GSET,NAT,R0,NGL,NDIM,NBH,NPRO,THIS%OPSI,OPROJ)
          DEALLOCATE(OPROJ)
!
!         ==============================================================
!         ==  DIVIDE BY WAVE FUNCTION MASS                            ==
!         ==============================================================
          ALLOCATE(MARR(NGL))
          CALL PLANEWAVE$GETR8A('G2',NGL,MARR)
          IF(ASSOCIATED(GSET%DMPSI)) THEN
            DO IG=1,NGL
!             SVAR=1.D0+ANNEE+GSET%DMPSI(IG)
              SVAR=1.D0+ANNEE 
              MARR(IG)=DELT**2/(SVAR*GSET%MPSI(IG))
            ENDDO
          ELSE
            SVAR=DELT**2/(1.D0+ANNEE)
            DO IG=1,NGL
!MARR(IG)=EMASS*(1.D0+EMASSCG2*MARR(IG))  !OLD VERSION
!MARR(IG)=DELT**2/(MARR(IG)*(1.D0+ANNEE))
              MARR(IG)=SVAR/GSET%MPSI(IG)
            ENDDO
          END IF
          DO IB=1,NBH
            DO IDIM=1,NDIM
              DO IG=1,NGL
                THIS%OPSI(IG,IDIM,IB)=MARR(IG)*THIS%OPSI(IG,IDIM,IB)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(MARR)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  UPDATE STRUCTUREFACTOR, ETC                                 ==
!     ==================================================================
      ALLOCATE(RP(3,NAT))
      ALLOCATE(RM(3,NAT))
      CALL ATOMLIST$GETR8A('R(+)',0,3*NAT,RP)
      IF(TSTRESS) THEN
!       == PREDICT NEW POSITIONS =======================================
        CALL CELL$GETR8A('TP',9,RBAS)
        CALL CELL$GETR8A('MAPTOCELL',9,MAPTOCELL)
        CALL ATOMLIST$GETR8A('R(-)',0,3*NAT,RM)
        DO IAT=1,NAT
          RP(:,IAT)=RP(:,IAT)-RM(:,IAT)+MATMUL(MAPTOCELL,RM(:,IAT))
        ENDDO
!       ==  SCALING FACTOR FOR WAVE FUNCTIONS ==========================
        CALL GBASS(RBAS,GBAS,CELLVOL)
        CELLSCALE=MAPTOCELL(1,1)*(MAPTOCELL(2,2)*MAPTOCELL(3,3)  &
     &                           -MAPTOCELL(3,2)*MAPTOCELL(2,3)) &
     &           +MAPTOCELL(2,1)*(MAPTOCELL(3,2)*MAPTOCELL(1,3)  &
     &                           -MAPTOCELL(1,2)*MAPTOCELL(3,3)) &
     &           +MAPTOCELL(3,1)*(MAPTOCELL(1,2)*MAPTOCELL(2,3)  &
     &                           -MAPTOCELL(2,2)*MAPTOCELL(1,3))
        CELLSCALE=1.D0/SQRT(CELLSCALE)
!       == NOW K-POINT DEPENDENT DATA ==================================
        DO IKPT=1,NKPT
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
!
!         ==  UPDATE PROJECTOR FUNCTIONS ===============================
          ALLOCATE(G2(NGL))
          CALL PLANEWAVE$GETR8A('G2',NGL,G2)
          IND=0
          DO ISP=1,MAP%NSP
            CALL SETUP$ISELECT(ISP)
            DO LN=1,MAP%LNX(ISP)
              IND=IND+1
              CALL SETUP$GETFOFG('PRO',.FALSE.,LN,NGL,G2,CELLVOL,GSET%PRO(1,IND))
              CALL SETUP$GETFOFG('PRO',.TRUE.,LN,NGL,G2,CELLVOL,GSET%DPRO(1,IND))
            ENDDO
          ENDDO
          DEALLOCATE(G2)
!
!         ==  UPDATE SPHERICAL HARMONICS  ==============================
          ALLOCATE(GVEC(3,NGL))
          CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVEC)
          LMX=MAP%LMX
          ALLOCATE(YLM(LMX))
          DO IG=1,NGL
            CALL GETYLM(LMX,GVEC(1,IG),YLM)
            GSET%YLM(IG,:)=YLM(:) 
          ENDDO        
          DEALLOCATE(YLM)
!
!         == NOW THE STRAINED SPHERICAL HARMONICS ======================
          CALL WAVES_STRAINEDYLM(NGL,LMX,GVEC,GSET%YLM,GSET%SYLM)
          DEALLOCATE(GVEC)
!
!         == RESCALE WAVE FUNCTIONS ====================================
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            THIS%PSIM(:,:,:)=THIS%PSIM(:,:,:)*CELLSCALE
            THIS%OPSI(:,:,:)=THIS%OPSI(:,:,:)*CELLSCALE 
          ENDDO
        ENDDO
      ELSE
        CELLSCALE=1.D0
      END IF
!
!     ==================================================================
!     ==  NOW ORTHOGONALIZE                                           ==
!     ==================================================================
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          NGL=GSET%NGL
          NBH=THIS%NBH
          NB=THIS%NB
!
!         ==============================================================
!         ==  CALCULATE PROJECTIONS FOR THE NEW POSITIONS             ==
!         ==============================================================
          CALL WAVES_PROJECTIONS(MAP,GSET,NAT,RP,NGL,NDIM,NBH,NPRO &
     &                                             ,THIS%PSIM,THIS%PROJ)
          ALLOCATE(OPROJ(NDIM,NBH,NPRO))
          CALL WAVES_PROJECTIONS(MAP,GSET,NAT,RP,NGL,NDIM,NBH,NPRO &
     &                                                 ,THIS%OPSI,OPROJ)
!
!         ==============================================================
!         ==  1C-OVERLAP OF <PSI0|PSI0>, <OPSI|PSI0> AND <OPSI|OPSI>  ==
!         ==============================================================
          ALLOCATE(MAT(NB,NB))
          CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIM,NBH,NB,NPRO &
         &                    ,THIS%PROJ,THIS%PROJ,MAT)
          ALLOCATE(OMAT(NB,NB))
          CALL WAVES_1COVERLAP(.FALSE.,MAP,NDIM,NBH,NB,NPRO &
       &                      ,OPROJ,THIS%PROJ,OMAT)
          ALLOCATE(OOMAT(NB,NB))
          CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIM,NBH,NB,NPRO &
       &                      ,OPROJ,OPROJ,OOMAT)
!
!         ==============================================================
!         ==  NOW ADD OVERLAP OF PSEUDO WAVE FUNCTIONS                ==
!         ==============================================================
          ALLOCATE(AUXMAT(NB,NB))
          CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB &
      &                     ,THIS%PSIM,THIS%PSIM,AUXMAT)
          DO I=1,NB
            DO J=1,NB
              MAT(I,J)=MAT(I,J)+AUXMAT(I,J)
            ENDDO
          ENDDO
          CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,NBH,NB &
      &                     ,THIS%OPSI,THIS%PSIM,AUXMAT)
          DO I=1,NB
            DO J=1,NB
              OMAT(I,J)=OMAT(I,J)+AUXMAT(I,J)
            ENDDO
          ENDDO
          CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB &
       &                    ,THIS%OPSI,THIS%OPSI,AUXMAT)
          DO I=1,NB
            DO J=1,NB
              OOMAT(I,J)=OOMAT(I,J)+AUXMAT(I,J)
            ENDDO
          ENDDO
          DEALLOCATE(AUXMAT)
!
!         ==================================================================
!         ==  CALCULATE LAGRANGE PARAMETERS                               ==
!         ==================================================================
          ALLOCATE(LAMBDA(NB,NB))
          LAMBDA(:,:)=THIS%RLAM0(:,:)
          DO I=1,NB
            DO J=I+1,NB
              CSVAR=0.5D0*(LAMBDA(I,J)+CONJG(LAMBDA(J,I)))
              LAMBDA(I,J)=CSVAR
              LAMBDA(J,I)=CONJG(CSVAR)
            ENDDO
          ENDDO
!===============================================================          
          IF(.NOT.TSAFEORTHO) THEN
            ALLOCATE(SMAP(NB))
            DO I=1,NB
              SMAP(I)=I
            ENDDO
            IF(TSWAPSTATES) THEN
              SVAR=0.D0
              DO I=1,NB
                DO J=I+1,NB
                  SVAR=MAX(SVAR,ABS(LAMBDA(I,J)))
                  SVAR=MAX(SVAR,ABS(LAMBDA(J,I)))
                  I1=SMAP(I)
                  J1=SMAP(J)
                  IF(REAL(LAMBDA(J1,J1)).LT.REAL(LAMBDA(I1,I1))) THEN
                    K=SMAP(I)
                    SMAP(I)=SMAP(J)
                    SMAP(J)=K
                  END IF
                ENDDO
              ENDDO  
!PRINT*,'SMAP ',SMAP
!PRINT*,'LAMBDA ',(REAL(LAMBDA(I,I))*27.211D0,I=1,NB)
!PRINT*,'MAX ',SVAR*27.211D0
              DO I=1,NB-1
                I1=SMAP(I)
                J1=SMAP(I+1)
                IF(REAL(LAMBDA(I1,I1)).GT.REAL(LAMBDA(J1,J1))) THEN
                  CALL ERROR$MSG('STATE ORDERING FAILED')
                  CALL ERROR$STOP('WAVES$ORTHOGONALIZE')
                END IF
              ENDDO
!!$PRINT*,'---------- LAMBDA SPIN: ',ISPIN,' ------------------'
!!$WRITE(*,FMT='("LAMBDA SMAP:",20I4)')SMAP
!!$WRITE(*,FMT='("LAMBDA IMAP:",20I4)') (I,I=1,NB)
!!$PRINT*,'MAX NON-DIAGNAL LAMBDA',SVAR*27.211D0
!!$WRITE(*,FMT='("LAMBDA BEFORE:",10F8.3)') (REAL(LAMBDA(I,I))*27.211D0,I=1,NB)
!!$WRITE(*,FMT='("LAMBDA OCC:   ",10F8.3)')OCC(:,IKPT,ISPIN)
!!$ELSE
!!$   PRINT*,'LAMBDA SMAP SWITCHED OFF!!'
            END IF
          END IF
!===============================================================          
          DO I=1,NB
            DO J=1,NB
              OCCI=OCC(I,IKPT,ISPIN)
              OCCJ=OCC(J,IKPT,ISPIN)
              IF(OCCI+OCCJ.LT.DSMALL)OCCJ=DSMALL
              LAMBDA(I,J)=LAMBDA(I,J)*0.5D0*OCCI/(OCCI+OCCJ)
            ENDDO
          ENDDO
!
          IF(TSAFEORTHO) THEN
            CALL PLANEWAVE$GETL4('TINV',TINV)
            IF(TINV) THEN
              ALLOCATE(RMAT(NB,NB))
              RMAT=REAL(MAT)
              ALLOCATE(ROMAT(NB,NB))
              ROMAT=REAL(OMAT)
              ALLOCATE(ROOMAT(NB,NB))
              ROOMAT=REAL(OOMAT)
              ALLOCATE(RLAMBDA(NB,NB))
              RLAMBDA=REAL(LAMBDA)
              CALL WAVES_ORTHO_X(NB,OCC(1,IKPT,ISPIN) &
       &                        ,ROOMAT,RMAT,ROMAT,RLAMBDA)
              LAMBDA=CMPLX(RLAMBDA,0.D0,8)
              DEALLOCATE(RLAMBDA)
              DEALLOCATE(RMAT)
              DEALLOCATE(ROMAT)
              DEALLOCATE(ROOMAT)
            ELSE
              CALL WAVES_ORTHO_X_C(NB,OCC(1,IKPT,ISPIN),OOMAT,MAT,OMAT,LAMBDA)
            END IF
          ELSE
            CALL WAVES_ORTHO_Y_C(NB,MAT,OMAT,OOMAT,LAMBDA,SMAP)
          END IF
          DEALLOCATE(MAT)
          DEALLOCATE(OMAT)
          DEALLOCATE(OOMAT)
          IF(.NOT.TSAFEORTHO)DEALLOCATE(SMAP)
!
!         ==================================================================
!         ==  CALCULATE |PSI(+)>=|PSI>+|CHI>LAMBDA                        ==
!         ==================================================================
          CALL WAVES_ADDOPSI(NGL,NDIM,NBH,NB,THIS%PSIM,THIS%OPSI,LAMBDA)
          DEALLOCATE(THIS%OPSI)
PRINT*,'WARNING FROM WAVES$ORTHOGONALIZE:'
PRINT*,'MAKE SURE THAT PDOS AND GRAPHICS PICK UP A CONSISTENT SET OF '
PRINT*,'WAVE FUNCTIONS  AND PROJECTOR FUNCTIONS'
          CALL WAVES_ADDOPROJ(NPRO,NDIM,NBH,NB,THIS%PROJ,OPROJ,LAMBDA)
          DEALLOCATE(OPROJ)
!
!         ==================================================================
!         ==  RESCALE GAMMA                                               ==
!         ==================================================================
!         == MASS IS NOT INCLUDED HERE (MASS TENSOR IS TAKEN CARE OF IN 
!         == PSIBAR AND (1/M)O|PSI> 
          DO I=1,NB
            DO J=I,NB
              LAMBDA(I,J)=0.5D0*(LAMBDA(I,J)+CONJG(LAMBDA(J,I)))
              LAMBDA(J,I)=CONJG(LAMBDA(I,J))
            ENDDO
          ENDDO
!PRINT*,'LAMBDA[EV] ',(REAL(LAMBDA(I,I))*27.211D0,I=1,NB)
          IF(.NOT.ASSOCIATED(THIS%RLAM0))ALLOCATE(THIS%RLAM0(NB,NB))
          THIS%RLAM0(:,:)=LAMBDA(:,:)
          DEALLOCATE(LAMBDA)
!
!         ==================================================================
!         ==  TEST ORTHONORMALITY                                         ==
!         ==================================================================
          IF(TTEST) THEN
            ALLOCATE(OPROJ(NDIM,NBH,NPRO))
            CALL WAVES_PROJECTIONS(MAP,GSET,NAT,RP,NGL,NDIM,NBH,NPRO,THIS%PSIM,OPROJ)
            ALLOCATE(AUXMAT(NB,NB))
            CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIM,NBH,NB,NPRO,OPROJ,OPROJ,AUXMAT)
            DEALLOCATE(OPROJ)
            CSUM=(0.D0,0.D0)
            DO I=1,NB
              CSUM=CSUM+AUXMAT(I,I)*OCC(I,IKPT,ISPIN)
            ENDDO
!PRINT*,'1C-CHARGE AFTER ORTHOGONALIZATION ',CSUM
            ALLOCATE(MAT(NB,NB))
            CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,THIS%PSIM,THIS%PSIM,MAT)
            CSUM=(0.D0,0.D0)
            DO I=1,NB
              CSUM=CSUM+MAT(I,I)*OCC(I,IKPT,ISPIN)
            ENDDO
!PRINT*,'PS-CHARGE AFTER ORTHOGONALIZATION ',CSUM
            DO J=1,NB
              DO I=1,NB
                MAT(I,J)=MAT(I,J)+AUXMAT(I,J)
              ENDDO
            ENDDO
            ALLOCATE(NORM(NB))
            CALL LIB$DIAGC8(NB,MAT,NORM,AUXMAT)
!CALL CDIAG(NB,NB,MAT,NORM,AUXMAT)
            DEALLOCATE(MAT)
            DEALLOCATE(AUXMAT)
!PRINT*,'NORM ',NORM
            DO I=1,NB
              IF(ABS(NORM(I)-1.D0).GT.1.D-4) THEN
                CALL ERROR$MSG('ORTHOGONALIZATION FAILED')
                CALL ERROR$I4VAL('STATE',I)
                CALL ERROR$R8VAL('NORM',NORM(I))
                CALL ERROR$STOP('WAVES$ORTHOGONALIZE')
              END IF
            ENDDO
            DEALLOCATE(NORM)
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(OCC)
!
!     ==================================================================
!     ==  NOW TRANSFORM MACK INTO ORIGINAL CELL                       ==
!     ==================================================================
      IF(TSTRESS) THEN
        CALL CELL$GETR8A('T0',9,RBAS)
        DO IKPT=1,NKPT
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            THIS%PSIM(:,:,:)=THIS%PSIM(:,:,:)/CELLSCALE
          ENDDO   
        ENDDO
      END IF
      DEALLOCATE(RP)
      DEALLOCATE(R0)
      DEALLOCATE(RM)
!
!     ==================================================================
!     ==  CALCULATE SECOND PART OF WAVE FUNCTION KINETIC ENERGY       ==
!     ==================================================================
      CALL WAVES_WAVEKINETIC(WAVEEKIN2)
                             CALL TIMING$CLOCKOFF('WAVES$ORTHOGONALIZE')
                                    CALL TRACE$POP
      RETURN
      END
!
!      .................................................................
       SUBROUTINE WAVES_OPROJ(LNX,LOX,DO,NDIM,LMNX,NB,PROJ,OPROJ)
!      *****************************************************************
!      **                                                             **
!      **  DO<PRO|PSI>                                                **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: LNX
       INTEGER(4),INTENT(IN) :: LOX(LNX)
       REAL(8)   ,INTENT(IN) :: DO(LNX,LNX)
       INTEGER(4),INTENT(IN) :: LMNX
       INTEGER(4),INTENT(IN) :: NDIM
       INTEGER(4),INTENT(IN) :: NB
       COMPLEX(8),INTENT(IN) :: PROJ(NDIM,NB,LMNX)
       COMPLEX(8),INTENT(OUT):: OPROJ(NDIM,NB,LMNX)
       INTEGER(4)            :: LMN1,LMN2
       INTEGER(4)            :: LMN10,LMN20
       INTEGER(4)            :: LN1,LN2
       INTEGER(4)            :: L1,L2
       INTEGER(4)            :: IB,IDIM,M
       REAL(8)               :: DOVER1
!      *****************************************************************
       OPROJ(:,:,:)=(0.D0,0.D0)
       LMN10=0
       DO LN1=1,LNX
         L1=LOX(LN1)
         LMN20=0
         DO LN2=1,LNX
           L2=LOX(LN2)
           IF(L1.EQ.L2) THEN
             DOVER1=DO(LN1,LN2)
             DO M=1,2*L1+1
               LMN1=LMN10+M
               LMN2=LMN20+M
               DO IDIM=1,NDIM
                 DO IB=1,NB
                   OPROJ(IDIM,IB,LMN1)=OPROJ(IDIM,IB,LMN1) &
     &                          +DOVER1*PROJ(IDIM,IB,LMN2)
                 ENDDO
               ENDDO
             END DO
           END IF
           LMN20=LMN20+2*L2+1
         ENDDO
         LMN10=LMN10+2*L1+1
       ENDDO
       RETURN
       END
!
!      .................................................................
       SUBROUTINE WAVES_ADDOPSI(NGL,NDIM,NBH,NB,PSIBAR,OPSI,LAMBDA)
!      *****************************************************************
!      **                                                             **
!      **  |PSI(+)>=|PSIBAR>+O|PSI(0)>*LAMBDA                         **
!      **                                                             **
!      **                                                             **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)   :: NGL
       INTEGER(4),INTENT(IN)   :: NDIM
       INTEGER(4),INTENT(IN)   :: NBH
       INTEGER(4),INTENT(IN)   :: NB
       COMPLEX(8),INTENT(INOUT):: PSIBAR(NGL*NDIM,NBH)
       COMPLEX(8),INTENT(INOUT):: OPSI(NGL*NDIM,NBH)
       COMPLEX(8),INTENT(IN)   :: LAMBDA(NB,NB)
       LOGICAL(4),PARAMETER    :: TESSL=.TRUE.
       LOGICAL(4)              :: TINV
       INTEGER(4)              :: IBH1,IBH2,I,IDIM,IBH
       INTEGER(4)              :: IB1A,IB1B,IB2A,IB2B
       COMPLEX(8),ALLOCATABLE  :: LAMBDA1(:,:)
       COMPLEX(8),ALLOCATABLE  :: LAMBDA2(:,:)
       COMPLEX(8),ALLOCATABLE  :: TPSI(:)
       COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
       COMPLEX(8)              :: CSVAR,CSVAR1,CSVAR2
       INTEGER(4)              :: NGLNDIM
!      *****************************************************************
                               CALL TIMING$CLOCKON('WAVES_ADDOPSI')
       TINV=NBH.NE.NB
       NGLNDIM=NGL*NDIM
       IF(.NOT.TINV) THEN
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA,PSIBAR)
!        IF(TESSL) THEN
!          CALL ZGEMM('N','N',NGLNDIM,NBH,NBH,(1.D0,0.D0),OPSI,NGLNDIM &
!                     ,LAMBDA,NBH,(1.D0,0.D0),PSIBAR,NGLNDIM)
!        ELSE
!          DO IBH1=1,NBH
!            DO IBH2=1,NBH
!              CSVAR=LAMBDA(IBH2,IBH1)
!              DO I=1,NGLNDIM
!                PSIBAR(I,IBH1)=PSIBAR(I,IBH1)+OPSI(I,IBH2)*CSVAR
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDIF
       ELSE
         ALLOCATE(LAMBDA1(NBH,NBH))
         ALLOCATE(LAMBDA2(NBH,NBH))
         DO IBH1=1,NBH
           IB1A=2*IBH1-1
           IB1B=2*IBH1
           DO IBH2=1,NBH
             IB2A=2*IBH2-1
             IB2B=2*IBH2
             CSVAR1=   LAMBDA(IB1A,IB2A)+CI*LAMBDA(IB1A,IB2B)
             CSVAR2=CI*LAMBDA(IB1B,IB2A)-   LAMBDA(IB1B,IB2B)
             LAMBDA1(IBH1,IBH2)=0.5D0*(CSVAR1-CSVAR2)
             LAMBDA2(IBH1,IBH2)=0.5D0*(CSVAR1+CSVAR2)
           ENDDO
         ENDDO
!        == ADD O|PSI_+>LAMBDA1 =======================================
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA1,PSIBAR)
!        IF(TESSL) THEN
!          CALL ZGEMM('N','N',NGLNDIM,NBH,NBH,(1.D0,0.D0),OPSI,NGLNDIM &
!                     ,LAMBDA1,NBH,(1.D0,0.D0),PSIBAR,NGLNDIM)
!        ELSE 
!          DO IBH1=1,NBH
!            DO IBH2=1,NBH
!              CSVAR=LAMBDA1(IBH2,IBH1)
!              DO I=1,NGLNDIM
!                PSIBAR(I,IBH1)=PSIBAR(I,IBH1)+OPSI(I,IBH2)*CSVAR
!              ENDDO
!            ENDDO
!          ENDDO
!        END IF
         DEALLOCATE(LAMBDA1)
!        == INVERT OPSI ================================================
         ALLOCATE(TPSI(NGLNDIM))
         DO IBH1=1,NBH
           I=1
           DO IDIM=1,NDIM
             CALL PLANEWAVE$INVERTG(NGL,OPSI(I,IBH1),TPSI(I))
             I=I+NGL
           ENDDO
           OPSI(:,IBH1)=TPSI(:)
         ENDDO
         DEALLOCATE(TPSI)
!
!        == ADD O|PSI_+>LAMBDA1 =======================================
         CALL LIB$ADDPRODUCTC8(.FALSE.,NGLNDIM,NBH,NBH,OPSI,LAMBDA2,PSIBAR)
!        IF(TESSL) THEN
!          CALL ZGEMM('N','N',NGLNDIM,NBH,NBH,(1.D0,0.D0),OPSI,NGLNDIM &
!                     ,LAMBDA2,NBH,(1.D0,0.D0),PSIBAR,NGLNDIM)
!        ELSE 
!          DO IBH1=1,NBH
!            DO IBH2=1,NBH
!              CSVAR=LAMBDA2(IBH2,IBH1)
!              DO I=1,NGLNDIM
!                PSIBAR(I,IBH1)=PSIBAR(I,IBH1)+OPSI(I,IBH2)*CSVAR
!              ENDDO
!            ENDDO
!          ENDDO
!        END IF
         DEALLOCATE(LAMBDA2)
       END IF
                               CALL TIMING$CLOCKOFF('WAVES_ADDOPSI')
       RETURN
       END
!
!      .................................................................
       SUBROUTINE WAVES_ADDOPROJ(NPRO,NDIM,NBH,NB,PROJ,OPROJ,LAMBDA)
!      *****************************************************************
!      **                                                             **
!      **  |PSI(+)>=|PSIBAR>+O|PSI(0)>*LAMBDA                         **
!      **                                                             **
!      **                                                             **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)   :: NDIM
       INTEGER(4),INTENT(IN)   :: NBH
       INTEGER(4),INTENT(IN)   :: NB
       INTEGER(4),INTENT(IN)   :: NPRO
       COMPLEX(8),INTENT(INOUT):: PROJ(NDIM,NBH,NPRO)
       COMPLEX(8),INTENT(IN)   :: OPROJ(NDIM,NBH,NPRO)
       COMPLEX(8),INTENT(IN)   :: LAMBDA(NB,NB)
       LOGICAL(4),PARAMETER    :: TESSL=.TRUE.
       LOGICAL(4)              :: TINV
       INTEGER(4)              :: IBH1,IBH2,I,IDIM,IBH
       INTEGER(4)              :: IB1A,IB1B,IB2A,IB2B
       COMPLEX(8),ALLOCATABLE  :: LAMBDA1(:,:)
       COMPLEX(8),ALLOCATABLE  :: LAMBDA2(:,:)
       COMPLEX(8),ALLOCATABLE  :: TPSI(:)
       COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
       COMPLEX(8)              :: CSVAR,CSVAR1,CSVAR2
       INTEGER(4)              :: IPRO
!      *****************************************************************
       TINV=NBH.NE.NB
       IF(.NOT.TINV) THEN
         DO IBH1=1,NBH
           DO IBH2=1,NBH
             CSVAR=LAMBDA(IBH2,IBH1)
             DO IPRO=1,NPRO
               DO IDIM=1,NDIM
                 PROJ(IDIM,IBH1,IPRO)=PROJ(IDIM,IBH1,IPRO) &
     &                              +OPROJ(IDIM,IBH2,IPRO)*CSVAR
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ELSE
         ALLOCATE(LAMBDA1(NBH,NBH))
         ALLOCATE(LAMBDA2(NBH,NBH))
         DO IBH1=1,NBH
           IB1A=2*IBH1-1
           IB1B=2*IBH1
           DO IBH2=1,NBH
             IB2A=2*IBH2-1
             IB2B=2*IBH2
             CSVAR1=   LAMBDA(IB1A,IB2A)+CI*LAMBDA(IB1A,IB2B)
             CSVAR2=CI*LAMBDA(IB1B,IB2A)-   LAMBDA(IB1B,IB2B)
             LAMBDA1(IBH1,IBH2)=0.5D0*(CSVAR1-CSVAR2)
             LAMBDA2(IBH1,IBH2)=0.5D0*(CSVAR1+CSVAR2)
           ENDDO
         ENDDO
!        == ADD O|PSI_+>LAMBDA1 =======================================
         DO IBH1=1,NBH
           DO IBH2=1,NBH
             CSVAR1=LAMBDA1(IBH2,IBH1)
             CSVAR2=LAMBDA2(IBH2,IBH1)
             DO IPRO=1,NPRO
               DO IDIM=1,NDIM
                 PROJ(IDIM,IBH1,IPRO)=PROJ(IDIM,IBH1,IPRO) &
     &                             +OPROJ(IDIM,IBH2,IPRO) *CSVAR1 &
     &                       +CONJG(OPROJ(IDIM,IBH2,IPRO))*CSVAR2
               ENDDO
             ENDDO
           ENDDO
         ENDDO
         DEALLOCATE(LAMBDA1)
         DEALLOCATE(LAMBDA2)
       END IF
       RETURN
       END
!
!      ..............................................................
       SUBROUTINE WAVES_OVERLAP(TID,NGL,NDIM,NBH,NB,PSI1,PSI2,MAT)
!      **                                                          **
!      **  CALCULATES <PSI1|PSI2>                                  **
!      **                                                          **
       USE MPE_MODULE
       IMPLICIT NONE
       LOGICAL(4),INTENT(IN) :: TID !INDICATES THAT PSI1=PSI2
       INTEGER(4),INTENT(IN) :: NGL
       INTEGER(4),INTENT(IN) :: NDIM
       INTEGER(4),INTENT(IN) :: NBH
       INTEGER(4),INTENT(IN) :: NB
       COMPLEX(8),INTENT(IN) :: PSI1(NGL,NDIM,NBH)
       COMPLEX(8),INTENT(IN) :: PSI2(NGL,NDIM,NBH)
       COMPLEX(8),INTENT(OUT):: MAT(NB,NB)
       INTEGER(4)            :: IBH1,IBH2,IB1,IB2,IDIM,IG
       INTEGER(4)            :: IB1A,IB1B,IB2A,IB2B,I,J
       COMPLEX(8)            :: CSVARPP,CSVARPM,CSVARP,CSVARM
       COMPLEX(8)            :: CSVAR
       COMPLEX(8),ALLOCATABLE:: TMAT(:,:)
       REAL(8)               :: RE,IM
       REAL(8)               :: MAT2(2,2)
       LOGICAL(4)            :: TINV
!      **************************************************************
                             CALL TIMING$CLOCKON('WAVES_OVERLAP')
       TINV=(NB.NE.NBH)
       IF(.NOT.TINV) THEN
         IF(TID) THEN
           CALL PLANEWAVE$SCALARPRODUCT('=',NGL,NDIM,NB,PSI1,NB,PSI2,MAT)
         ELSE
           CALL PLANEWAVE$SCALARPRODUCT(' ',NGL,NDIM,NB,PSI1,NB,PSI2,MAT)
         END IF
         CALL MPE$COMBINE('+',MAT)
                             CALL TIMING$CLOCKOFF('WAVES_OVERLAP')
         RETURN
       ENDIF
!
!      =================================================================
!      == NOW DEAL WITH SUPER WAVE FUNCTIONS                          ==
!      =================================================================
!      == <PSI_+|PSI_+> ================================================
       ALLOCATE(TMAT(NBH,NBH))
       IF(TID) THEN
         CALL PLANEWAVE$SCALARPRODUCT('=',NGL,NDIM,NBH,PSI1,NBH,PSI2,TMAT)
       ELSE
         CALL PLANEWAVE$SCALARPRODUCT(' ',NGL,NDIM,NBH,PSI1,NBH,PSI2,TMAT)
       END IF
       DO IBH1=1,NBH
         IB1A=2*IBH1-1
         IB1B=MIN(2*IBH1,NB)
         DO IBH2=1,NBH
           IB2A=2*IBH2-1
           IB2B=MIN(2*IBH2,NB)
           RE=0.5D0* REAL(TMAT(IBH1,IBH2))
           IM=0.5D0*AIMAG(TMAT(IBH1,IBH2))
!          == WATCH THE ORDER OF THE FOLLOWING 4 LINES!!! =============
!          == BECAUSE IB1B=IB1A FOR 2*IBH1>NB =========================
!          == AND     IB2B=IB2A FOR 2*IBH2>NB =========================
           MAT2(1,1)=+RE
           MAT2(1,2)=+IM
           MAT2(2,1)=-IM
           MAT2(2,2)=+RE
           DO I=IB1A,IB1B
             DO J=IB2A,IB2B
               MAT(I,J)=CMPLX(MAT2(I-IB1A+1,J-IB2A+1),0.D0,8)
             ENDDO
           ENDDO
         ENDDO
       ENDDO
!      == NOW <PSI-|PSI+> ==============================================
       CALL PLANEWAVE$SCALARPRODUCT('-',NGL,NDIM,NBH,PSI1,NBH,PSI2,TMAT)
       DO IBH1=1,NBH
         IB1A=2*IBH1-1
         IB1B=MIN(2*IBH1,NB)
         DO IBH2=1,NBH
!          == <PSI-(I)|PSI+(J)> ========================================
           IB2A=2*IBH2-1
           IB2B=MIN(2*IBH2,NB)
           RE=0.5D0* REAL(TMAT(IBH1,IBH2))
           IM=0.5D0*AIMAG(TMAT(IBH1,IBH2))
           MAT2(1,1)=+RE
           MAT2(1,2)=+IM
           MAT2(2,1)=+IM
           MAT2(2,2)=-RE
           DO I=IB1A,IB1B
             DO J=IB2A,IB2B
               MAT(I,J)=MAT(I,J)+CMPLX(MAT2(I-IB1A+1,J-IB2A+1),0.D0,8)
             ENDDO
           ENDDO
         ENDDO
       ENDDO
       DEALLOCATE(TMAT)
       CALL MPE$COMBINE('+',MAT)
                             CALL TIMING$CLOCKOFF('WAVES_OVERLAP')
       RETURN
       END
!
!     ..................................................................
      SUBROUTINE WAVES_1COVERLAP(TID,MAP,NDIM,NBH,NB,NPRO,PROJ1,PROJ2,MAT)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE 1C CONTRIBUTION TO THE OVERLAP MATRIX        **
!     **  ( RESULT IS ADDED TO "MAT"! )                               **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY: MAP_TYPE
      IMPLICIT NONE
      LOGICAL(4)  ,INTENT(IN)   :: TID
      TYPE(MAP_TYPE),INTENT(IN) :: MAP
      INTEGER(4)  ,INTENT(IN)   :: NDIM
      INTEGER(4)  ,INTENT(IN)   :: NBH
      INTEGER(4)  ,INTENT(IN)   :: NB
      INTEGER(4)  ,INTENT(IN)   :: NPRO
      COMPLEX(8)  ,INTENT(IN)   :: PROJ1(NDIM,NBH,NPRO) ! <P|PSI1>
      COMPLEX(8)  ,INTENT(IN)   :: PROJ2(NDIM,NBH,NPRO) ! <P|PSI2>
      COMPLEX(8)  ,INTENT(OUT)  :: MAT(NB,NB)
      REAL(8)     ,ALLOCATABLE  :: DOVER(:,:)
      INTEGER(4)                :: NTASKS,THISTASK
      INTEGER(4)                :: ICOUNT,IAT,ISP
      INTEGER(4)                :: LMN1,LN1,L1
      INTEGER(4)                :: LMN2,LN2,L2
      INTEGER(4)                :: IB1,IB2,IPRO,IPRO1,IPRO2,M
      INTEGER(4)                :: IBH1,IBH2,IDIM,LNX
      INTEGER(4)                :: IB1A,IB1B,IB2A,IB2B,I,J
      COMPLEX(8)                :: CSVAR,CSVAR1,CSVAR2
      REAL(8)                   :: RE1,IM1,RE2,IM2
      REAL(8)                   :: RMAT(2,2)
      LOGICAL(4)                :: TINV
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
      MAT(:,:)=(0.D0,0.D0)
      TINV=NB.NE.NBH
      ICOUNT=0
      IPRO=0
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
!       __ SELECTION FOR PARALLEL PROCESSING____________________________
        ICOUNT=ICOUNT+1
        IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) THEN
          IPRO=IPRO+MAP%LMNX(ISP)
          CYCLE
        END IF
!       __ NOW CONTINUE_________________________________________________
        LNX=MAP%LNX(ISP)
        ALLOCATE(DOVER(LNX,LNX))
        CALL SETUP$1COVERLAP(ISP,LNX,DOVER)
!
        DO IBH1=1,NBH
          DO IBH2=1,NBH
            CSVAR1=(0.D0,0.D0)
            CSVAR2=(0.D0,0.D0)
            IPRO1=IPRO
            DO LN1=1,LNX
              L1=MAP%LOX(LN1,ISP)
              IPRO2=IPRO
              DO LN2=1,LNX
                L2=MAP%LOX(LN2,ISP)
                IF(L1.EQ.L2) THEN
                  CSVAR=(0.D0,0.D0)
                  DO IDIM=1,NDIM
                    DO M=1,2*L1+1
                      CSVAR=CSVAR + CONJG(PROJ1(IDIM,IBH1,IPRO1+M)) &
     &                            *       PROJ2(IDIM,IBH2,IPRO2+M)
                    ENDDO
                  ENDDO
                  CSVAR1=CSVAR1+CSVAR*DOVER(LN1,LN2)
                  IF(TINV) THEN
                    CSVAR=(0.D0,0.D0)
                    DO IDIM=1,NDIM
                      DO M=1,2*L1+1
                        CSVAR=CSVAR + PROJ1(IDIM,IBH1,IPRO1+M) &
     &                              * PROJ2(IDIM,IBH2,IPRO2+M)
                      ENDDO
                    ENDDO
                    CSVAR2=CSVAR2+CSVAR*DOVER(LN1,LN2)
                  END IF
                END IF
                IPRO2=IPRO2+2*L2+1
              ENDDO
              IPRO1=IPRO1+2*L1+1
            ENDDO

            IF(.NOT.TINV) THEN
              MAT(IBH1,IBH2)=MAT(IBH1,IBH2)+CSVAR1
            ELSE
              RE1=  REAL(CSVAR1)
              IM1= AIMAG(CSVAR1)
              RE2=  REAL(CSVAR2)
              IM2=-AIMAG(CSVAR2)  ! COMPLEX CONJUGATE OF CSVAR2 USED
              RMAT(1,1)=0.5D0*( RE1+RE2)
              RMAT(1,2)=0.5D0*( IM1-IM2)
              RMAT(2,1)=0.5D0*(-IM1-IM2)
              RMAT(2,2)=0.5D0*( RE1-RE2)
              IB1A=2*IBH1-1
              IB1B=MIN(2*IBH1,NB)
              IB2A=2*IBH2-1
              IB2B=MIN(2*IBH2,NB)
              DO I=IB1A,IB1B
                DO J=IB2A,IB2B
                  MAT(I,J)=MAT(I,J)+CMPLX(RMAT(I-IB1A+1,J-IB2A+1),0.D0,8)
                ENDDO
              ENDDO
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(DOVER)
        IPRO=IPRO+MAP%LMNX(ISP)
      ENDDO
      CALL MPE$COMBINE('+',MAT)
      RETURN
      END
!
!      .................................................................
       SUBROUTINE WAVES_ORTHO_Y_C(NB,PHIPHI,CHIPHI,CHICHI,X,MAP)
!      **                                                             **
!      **  CALCULATE LAGRANGE MULTIPLIERS FOR ORTHOGONALIZATION       **
!      **    |PHI(I)>=|PHI(I)>+SUM_J |CHI(J)>X(J,I)                   **
!      **  WITH                                                       **
!      **    X(I>J)=0                                                 **
!      **                                                             **
!      **  ATTENTION!! CHIPHI0=<CHI|O|PHI>                            **
!      **        AND   PHIPHI0=<PHI|O|PHI>                            **
!      **        CONVERSION IS DONE IN THE INITIALIZATION             **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NB
       COMPLEX(8),INTENT(IN) :: PHIPHI(NB,NB) !<PHI_I|O|PHI_J>
       COMPLEX(8),INTENT(IN) :: CHIPHI(NB,NB) !<PHI_I|O|CHI>
       COMPLEX(8),INTENT(IN) :: CHICHI(NB,NB) !<CHI_I|O|CHI_J>
       COMPLEX(8),INTENT(OUT):: X(NB,NB)      ! X(I>J)=0
       INTEGER(4),INTENT(IN) :: MAP(NB)
       COMPLEX(8)            :: A(NB,NB)
       COMPLEX(8)            :: B(NB,NB)
       COMPLEX(8)            :: C(NB,NB)
       COMPLEX(8)            :: ALPHA(NB,NB)   ! (I>J)=0
       COMPLEX(8)            :: WORK(NB,NB)    
       COMPLEX(8)            :: Z(NB)
       INTEGER(4)            :: I,J,K,L,N,M
       INTEGER(4)            :: N0
       COMPLEX(8)            :: CSVAR
       REAL(8)               :: SVAR,SVAR1,SVAR2
       REAL(8)               :: MAXDEV
       LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
       REAL(8)   ,PARAMETER  :: TOL=1.D-10
       INTEGER(4)            :: NU,NU0,IU,JU
!      *****************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_Y_C')
!
!      =================================================================
!      ==  INITIALIZE                                                 ==
!      =================================================================
       DO I=1,NB
         DO J=1,NB
           A(I,J)=PHIPHI(I,J)
           B(I,J)=CHIPHI(I,J)
           C(I,J)=CHICHI(I,J)
           X(I,J)=(0.D0,0.D0)
           ALPHA(I,J)=(0.D0,0.D0)
         ENDDO
         A(I,I)=A(I,I)-(1.D0,0.D0)
         ALPHA(I,I)=(1.D0,0.D0)
       ENDDO
!
!      =================================================================
!      ==  ORTHOGONALIZATION LOOP                                     ==
!      =================================================================
!                            CALL TRACE$PASS('BEFORE ORTHOGONALIZATION LOOP')
       DO NU=1,NB
         N=MAP(NU)
!
!        ===============================================================
!        == NORMALIZE PHI(N)                                          ==
!        == PHI(N)=PHI(N)+CHI(N)*Z(N)                                 ==
!        ===============================================================
!                            CALL TRACE$PASS('NORMALIZE')
         SVAR=CONJG(B(N,N))*B(N,N)-A(N,N)*C(N,N)-AIMAG(B(N,N))**2
         SVAR1=-REAL(B(N,N))
         IF(SVAR.GE.0.D0) THEN
           SVAR2=SQRT(SVAR)
           IF(SVAR1*SVAR2.GE.0.D0) THEN
             Z(N)=SVAR1-SVAR2
           ELSE
             Z(N)=SVAR1+SVAR2
           END IF
         ELSE
           PRINT*,'ORTHOGONALIZATION FAILED! TRYING BEST APPROXIMATION...'
           Z(N)=SVAR1
         END IF
         Z(N)=Z(N)/REAL(C(N,N))
!
!        ===============================================================
!        == NOW UPDATE MATRICES                                       ==
!        ===============================================================
         NU0=NU   !SET N0=N FOR FAST CALCULATION AND N0=1 FOR TESTS
!N0=1
         DO I=1,NB
!          == A(N,M)+B(N,N)*DELTA(M)=0  ======================
           X(I,N)=X(I,N)+ALPHA(I,N)*Z(N)
         ENDDO           
         DO J=1,NB
           A(N,J)=A(N,J)+CONJG(Z(N))*B(N,J)
         ENDDO
         DO I=1,NB
           A(I,N)=A(I,N)+CONJG(B(N,I))*Z(N) 
         ENDDO
         A(N,N)=A(N,N)+CONJG(Z(N))*C(N,N)*Z(N)
         DO I=1,NB
           B(I,N)=B(I,N)+C(I,N)*Z(N)
         ENDDO
         IF(TTEST) THEN
           CALL TESTA(N,N,CSVAR)
           IF(ABS(CSVAR).GT.TOL) THEN
             WRITE(*,FMT='("NORMALIZATION OF PHI(",I4,")")')N
             PRINT*,N,N,CSVAR
           END IF
         END IF
!
!        ===============================================================
!        == ORTHOGONALIZE HIGHER PHI'S TO THIS PHI                    ==
!        == PHI(J)=PHI(J)+CHI(N)*Z(J)       J>N                       ==
!        ===============================================================
         DO IU=1,NU
           I=MAP(IU)
           Z(I)=(0.D0,0.D0)
         ENDDO
         DO IU=NU+1,NB
           I=MAP(IU)
           Z(I)=-CONJG(A(I,N)/B(N,N))
         ENDDO
!               CALL TRACE$PASS('ORTHOGONALIZE HIGHER PHIS TO THIS PHI')
!
!        ===============================================================
!        == NOW UPDATE MATRICES                                       ==
!        ===============================================================
         NU0=NU+1   !SET N0=N FOR FAST CALCULATION AND N0=1 FOR TESTS
!N0=1
         DO I=1,NB
!          == A(N,M)+B(N,N)*DELTA(M)=0  ======================
           DO JU=NU0,NB
             J=MAP(JU)
             X(I,J)=X(I,J)+ALPHA(I,N)*Z(J)
           ENDDO
         ENDDO           
         DO IU=NU0,NB
           I=MAP(IU)
           DO J=1,NB
             A(I,J)=A(I,J)+CONJG(Z(I))*B(N,J)
           ENDDO
         ENDDO
         DO I=1,NB
           DO JU=NU0,NB
             J=MAP(JU)
             A(I,J)=A(I,J)+CONJG(B(N,I))*Z(J) 
           ENDDO
         ENDDO
         DO IU=NU0,NB
           I=MAP(IU)
           DO JU=NU0,NB
             J=MAP(JU)
             A(I,J)=A(I,J)+CONJG(Z(I))*C(N,N)*Z(J)
           ENDDO
         ENDDO
         DO I=1,NB
           DO JU=NU0,NB
             J=MAP(JU)
             B(I,J)=B(I,J)+C(I,N)*Z(J)
           ENDDO
         ENDDO
         IF(TTEST) THEN
           DO IU=1,NU
             I=MAP(IU)
             DO JU=NU,NB
               J=MAP(JU)
               CALL TESTA(I,J,CSVAR)
               IF(ABS(CSVAR).GT.TOL) THEN
                 WRITE(*,FMT='("HIGHER PHIS ORTHOGONALIZED TO PHI(",I4,")")')N
                 PRINT*,I,J,CSVAR
               END IF
             ENDDO
           ENDDO
         END IF
!
!        ===============================================================
!        == ORTHOGONALIZE HIGHER CHI'S TO THIS PHI                    ==
!        == CHI(M)=CHI(M)+CHI(N)*DELTA(M)   M>N                       ==
!        ===============================================================
!               CALL TRACE$PASS('ORTHOGONALIZE HIGHER CHIS TO THIS PHI')
         DO IU=1,NU
           I=MAP(IU)
           Z(I)=(0.D0,0.D0)
         ENDDO
         DO IU=NU+1,NB
           I=MAP(IU)
!          == |CHI(J)>=|CHI(J)>+|CHI(N)>*Z(J) ==========================
!          == B(M,N)+B(N,N)*DELTA(M)=0
           Z(I)=-CONJG(B(I,N)/B(N,N))
         ENDDO
         NU0=NU+1   !SET N0=N+1 FOR FAST CALCULATION AND N0=1 FOR TESTS
!N0=1
         DO I=1,NB
           DO JU=NU0,NB
             J=MAP(JU)
             ALPHA(I,J)=ALPHA(I,J)+ALPHA(I,N)*Z(J)
           ENDDO
         ENDDO
         DO IU=NU0,NB
           I=MAP(IU)
           DO J=1,NB
             B(I,J)=B(I,J)+CONJG(Z(I))*B(N,J)
           ENDDO 
         ENDDO
         WORK(:,:)=0.D0
         DO IU=NU0,NB
           I=MAP(IU)
           DO J=1,NB
             WORK(I,J)=WORK(I,J)+CONJG(Z(I))*C(N,J) 
           ENDDO
         ENDDO
         DO I=1,NB
           DO JU=NU0,NB
             J=MAP(JU)
             WORK(I,J)=WORK(I,J)+C(I,N)*Z(J)
           ENDDO
         ENDDO
         DO IU=NU0,NB
           I=MAP(IU)
           DO JU=NU0,NB
             J=MAP(JU)
             WORK(I,J)=WORK(I,J)+CONJG(Z(I))*C(N,N)*Z(J)
           ENDDO
         ENDDO
         C(:,:)=C(:,:)+WORK(:,:)
         IF(TTEST) THEN
           DO IU=NU+1,NB
             I=MAP(IU)
             DO JU=1,NU
               J=MAP(JU)
               CALL TESTB(I,J,CSVAR)
               IF(ABS(CSVAR).GT.TOL) THEN
!                WRITE(*,FMT='("HIGHER CHIS ORTHOGONALIZED TO PHI(",I4,")")')N
!                PRINT*,I,J,CSVAR
               END IF
             ENDDO
            ENDDO
         END IF
       ENDDO
!
!      =================================================================
!      == TEST ORTHOGONALITY                                          ==
!      =================================================================
       IF(TTEST) THEN
         MAXDEV=0.D0
         DO I=1,NB
           DO J=1,NB
             CSVAR=PHIPHI(I,J)
             DO K=1,NB
               CSVAR=CSVAR+CONJG(X(K,I))*CHIPHI(K,J) &
      &                   +CONJG(CHIPHI(K,I))*X(K,J)
               DO L=1,NB
                 CSVAR=CSVAR+CONJG(X(K,I))*CHICHI(K,L)*X(L,J)
               ENDDO
             ENDDO
             IF(I.EQ.J) CSVAR=CSVAR-(1.D0,0.D0)
             MAXDEV=MAX(MAXDEV,ABS(CSVAR))
             IF(ABS(CSVAR).GT.TOL) THEN
               CALL ERROR$MSG('ORTHOGONALIZATION FAILED')
               CALL ERROR$I4VAL('I',I)
               CALL ERROR$I4VAL('J',J)
               CALL ERROR$C8VAL('<PHI(+)|O|PHI(+)>-1',CSVAR)
               CALL ERROR$STOP('WAVES_ORTHO_Y')
             END IF
           ENDDO
         ENDDO
         PRINT*,'MAX. DEVIATION IN ORTHO_Y_C',MAXDEV
       END IF
                             CALL TRACE$POP
       RETURN
       CONTAINS
!      .................................................................
       SUBROUTINE TESTA(I,J,CSVAR)
       INTEGER(4),INTENT(IN) :: I
       INTEGER(4),INTENT(IN) :: J
       COMPLEX(8),INTENT(OUT):: CSVAR
!      *****************************************************************
       CSVAR=PHIPHI(I,J)
       DO K=1,NB
         CSVAR=CSVAR+CONJG(X(K,I))*CHIPHI(K,J)+CONJG(CHIPHI(K,I))*X(K,J)
         DO L=1,NB
           CSVAR=CSVAR+CONJG(X(K,I))*CHICHI(K,L)*X(L,J)
         ENDDO
       ENDDO
       IF(I.EQ.J) CSVAR=CSVAR-(1.D0,0.D0)
       RETURN
       END SUBROUTINE TESTA
!      .................................................................
       SUBROUTINE TESTB(I,J,CSVAR)
       INTEGER(4),INTENT(IN) :: I
       INTEGER(4),INTENT(IN) :: J
       COMPLEX(8),INTENT(OUT):: CSVAR
!      *****************************************************************
       CSVAR=CHIPHI(I,J)
       DO K=1,NB
         CSVAR=CSVAR+CONJG(ALPHA(K,I))*CHIPHI(K,J)+CHIPHI(I,K)*X(K,J)
         DO L=1,NB
           CSVAR=CSVAR+CONJG(ALPHA(K,I))*CHICHI(K,L)*X(L,J)
         ENDDO
       ENDDO
       RETURN
       END SUBROUTINE TESTB
      END
!
!      .................................................................
       SUBROUTINE WAVES_ORTHO_Y(NB,PHIPHI0,CHIPHI0,CHICHI0,RLAMBDA)
!      **                                                             **
!      **  CALCULATE LAGRANGE MULTIPLIERS FOR ORTHOGONALIZATION       **
!      **    |PHI_N(+)>=|PHIBAR_N>+SUM_M |CHI_M(0)>LAMBDA(M,N)        **
!      **  WITH                                                       **
!      **    |CHI>=O|PHI>                                             **
!      **    |PHIBAR>=2|PHI(0)>-PHI(0)+H|PSI(0)>DT^2/M_PSI            **
!      **    LAMDA(I>J)=0                                             **
!      **                                                             **
!      **  ATTENTION!! CHIPHI0=<CHI|O|PHI>                            **
!      **        AND   PHIPHI0=<PHIBAR|O|PHIBAR>                      **
!      **        CONVERSION IS DONE IN THE INITIALIZATION             **
!      **                                                             **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NB
       REAL(8)   ,INTENT(IN) :: PHIPHI0(NB,NB) !<PHIBAR_I|O|PHIBAR_J>
       REAL(8)   ,INTENT(IN) :: CHIPHI0(NB,NB) !<PHIBAR_I|O|CHI>
       REAL(8)   ,INTENT(IN) :: CHICHI0(NB,NB) !<CHI_I|O|CHI_J>
       REAL(8)   ,INTENT(OUT):: RLAMBDA(NB,NB) ! (I>J)=0
       REAL(8)               :: PHIPHI(NB,NB)
       REAL(8)               :: CHIPHI(NB,NB)
       REAL(8)               :: CHICHI(NB,NB)
       REAL(8)               :: ALPHA(NB,NB)   ! (I>J)=0
       REAL(8)               :: WORK(NB,NB)   ! (I>J)=0
       REAL(8)               :: DELTA(NB)
       INTEGER(4)            :: I,J,K,L,N,M,M1,M2
       INTEGER(4)            :: NU,MU,M1U,IU
       REAL(8)               :: SVAR
       LOGICAL   ,PARAMETER  :: TPR=.FALSE.
       LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
       REAL(8)   ,PARAMETER  :: TOL=1.D-10
       REAL(8)               :: TEST(NB,NB)
       REAL(8)               :: TEST1(NB,NB)
       REAL(8)               :: TEST2(NB,NB)
       INTEGER(4)            :: MAP(NB)
!      *****************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_Y')
!
!      =================================================================
!      ==  INITIALIZE                                                 ==
!      =================================================================
       DO I=1,NB
         MAP(I)=I
         DO J=1,NB
           PHIPHI(I,J)=PHIPHI0(I,J)
           CHIPHI(I,J)=CHIPHI0(I,J)
           CHICHI(I,J)=CHICHI0(I,J)
           RLAMBDA(I,J)=0.D0
           ALPHA(I,J)=0.D0
         ENDDO
         PHIPHI(I,I)=PHIPHI(I,I)-1.D0
         ALPHA(I,I)=1.D0
       ENDDO
                             CALL TRACE$PASS('BEFORE TESTME')
       IF(TPR)CALL TESTME
!
!      =================================================================
!      ==  ORTHOGONALIZATION LOOP                                     ==
!      =================================================================
                             CALL TRACE$PASS('BEFORE ORTHOGONALIZATION LOOP')
       DO NU=1,NB
         N=MAP(NU)
!
!        ===============================================================
!        == NORMALIZE PHI(N)                                          ==
!        == PHI(N)=PHI(N)+CHI(N)*SVAR                                 ==
!        ===============================================================
                             CALL TRACE$PASS('NORMALIZE')
         SVAR=1.D0-PHIPHI(N,N)*CHICHI(N,N)/CHIPHI(N,N)**2
         IF(SVAR.GT.0.D0) THEN
           SVAR=CHIPHI(N,N)/CHICHI(N,N)*(DSQRT(SVAR)-1.D0)             
         ELSE
           PRINT*,'ORTHOGONALIZATION FAILED! TRYING BEST APPROXIMATION...'
           SVAR=-CHIPHI(N,N)/CHICHI(N,N)
         END IF       
         DO I=1,NB
           RLAMBDA(I,N)=RLAMBDA(I,N)+ALPHA(I,N)*SVAR
         ENDDO
         PHIPHI(N,N)=PHIPHI(N,N)+CHICHI(N,N)*SVAR**2
         DO I=1,NB
           PHIPHI(I,N)=PHIPHI(I,N)+CHIPHI(N,I)*SVAR
           PHIPHI(N,I)=PHIPHI(N,I)+CHIPHI(N,I)*SVAR
         ENDDO
         DO I=1,NB
           CHIPHI(I,N)=CHIPHI(I,N)+SVAR*CHICHI(N,I)
         ENDDO
         IF(TPR) THEN
           PRINT*,'AFTER NORMALIZATION'
           CALL TESTME
         END IF
!
!        ===============================================================
!        == ORTHOGONALIZE HIGHER PHI'S TO THIS PHI                    ==
!        == PHI(M)=PHI(M)+CHI(N)*DELTA(M)   M>N                       ==
!        ===============================================================
                CALL TRACE$PASS('ORTHOGONALIZE HIGHER PHIS TO THIS PHI')
         DELTA=0.D0
         DO MU=NU+1,NB
           M=MAP(MU)
!          == PHIPHI(N,M)+CHIPHI(N,N)*DELTA(M)=0  ======================
           DELTA(M)=-PHIPHI(N,M)/CHIPHI(N,N)
           DO IU=1,NU
             I=MAP(IU)
             RLAMBDA(I,M)=RLAMBDA(I,M)+ALPHA(I,N)*DELTA(M)
           ENDDO
         ENDDO           
         DO M1=1,NB
           DO M2=1,NB
             PHIPHI(M1,M2)=PHIPHI(M1,M2)+DELTA(M1)*CHIPHI(N,M2) &
     &                    +CHIPHI(N,M1)*DELTA(M2) &
     &                    +DELTA(M1)*CHICHI(N,N)*DELTA(M2)
           ENDDO
         ENDDO
         DO M1U=NU+1,NB
           M1=MAP(M1U)
           DO M2=1,NB
             CHIPHI(M2,M1)=CHIPHI(M2,M1)+DELTA(M1)*CHICHI(N,M2)
           ENDDO
         ENDDO
         IF(TPR) THEN
           WRITE(*,FMT='("HIGHER PHIS ORTHOGONALIZED TO PHI(",I4,")")')N
           CALL TESTME
         END IF
!
!        ===============================================================
!        == ORTHOGONALIZE HIGHER CHI'S TO THIS PHI                    ==
!        == CHI(M)=CHI(M)+CHI(N)*DELTA(M)   M>N                       ==
!        ===============================================================
                CALL TRACE$PASS('ORTHOGONALIZE HIGHER CHIS TO THIS PHI')
         DELTA=0.D0
         DO MU=NU+1,NB
           M=MAP(MU)
!          == |CHI(M)>=|CHI(M)>+|CHI(N)>*DELTA(M) ============================
!          == CHIPHI(M,N)+CHIPHI(N,N)*DELTA(M)=0
           DELTA(M)=-CHIPHI(M,N)/CHIPHI(N,N)
         ENDDO
         DO MU=NU+1,NB
           M=MAP(MU)
           DO I=1,NB
             ALPHA(I,M)=ALPHA(I,M)+ALPHA(I,N)*DELTA(M)
           ENDDO
           DO I=1,NB
             CHIPHI(M,I)=CHIPHI(M,I)+DELTA(M)*CHIPHI(N,I)
           ENDDO 
         ENDDO
         WORK(:,:)=0.D0
         DO M1=1,NB
           DO M2=1,NB
             WORK(M1,M2)=WORK(M1,M2)+DELTA(M1)*CHICHI(N,M2) &
        &                                    +CHICHI(M1,N)*DELTA(M2) &
        &                          +DELTA(M1)*CHICHI(N,N) *DELTA(M2)
           ENDDO
         ENDDO
         CHICHI(:,:)=CHICHI(:,:)+WORK(:,:)
         IF(TPR) THEN
           WRITE(*,FMT='("HIGHER CHIS ORTHOGONALIZED TO PHI(",I4,")")')
           CALL TESTME
         END IF
       ENDDO
!
!      =================================================================
!      == TEST ORTHOGONALITY                                          ==
!      =================================================================
       IF(TTEST) THEN
         CALL TESTME
         DO I=1,NB
           DO J=1,NB
             SVAR=PHIPHI0(I,J)
             DO K=1,NB
               SVAR=SVAR+CHIPHI0(K,I)*RLAMBDA(K,J)+RLAMBDA(K,I)*CHIPHI0(K,J)
               DO L=1,NB
                 SVAR=SVAR+RLAMBDA(K,I)*CHICHI0(K,L)*RLAMBDA(L,J)
               ENDDO
             ENDDO
             IF(I.EQ.J) SVAR=SVAR-1.D0
             IF(ABS(SVAR).GT.TOL) THEN
               CALL ERROR$MSG('ORTHOGONALIZATION FAILED')
               CALL ERROR$I4VAL('I',I)
               CALL ERROR$I4VAL('J',J)
               CALL ERROR$R8VAL('<PHI(+)|O|PHI(+)>-1',SVAR)
               CALL ERROR$STOP('WAVES_ORTHO_Y')
             END IF
           ENDDO
         ENDDO
       END IF
                             CALL TRACE$POP
       RETURN
       CONTAINS
!      ............................................................
       SUBROUTINE TESTME
       REAL(8)              :: SVAR,SVAR1,SVAR2
       WRITE(*,FMT='("STARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTARTSTART")')
       DO I=1,NB
         DO J=1,NB
           SVAR=PHIPHI0(I,J)
           SVAR1=0.D0
           SVAR2=0.D0
           DO K=1,NB
             SVAR=SVAR+RLAMBDA(K,I)*CHIPHI0(K,J)+CHIPHI0(K,I)*RLAMBDA(K,J)
             SVAR1=SVAR1+CHIPHI0(K,I)*ALPHA(K,J)
             DO L=1,NB
               SVAR=SVAR+RLAMBDA(K,I)*CHICHI0(K,L)*RLAMBDA(L,J)
               SVAR1=SVAR1+RLAMBDA(K,I)*CHICHI0(K,L)*ALPHA(L,J)
               SVAR2=SVAR2+ALPHA(K,I)*CHICHI0(K,L)*ALPHA(L,J)
             ENDDO
           ENDDO
           TEST(I,J)=SVAR
           TEST1(I,J)=SVAR1
           TEST2(I,J)=SVAR2
         ENDDO
         TEST(I,I)=TEST(I,I)-1.D0
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("TEST  =",8E10.2)')TEST(I,:)
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("PHIPHI=",8E10.2)')PHIPHI(I,:)
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("CHIPHI=",8E10.2)')CHIPHI(I,:)
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("CHICHI=",8E10.2)')CHICHI(I,:)
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("LAMBDA=",8F10.3)')RLAMBDA(I,:)
       ENDDO
       DO I=1,NB
         WRITE(*,FMT='("ALPHA =",8F10.3)')ALPHA(I,:)
       ENDDO
       DO I=1,NB
         DO J=1,NB
           IF(DABS(TEST(I,J)-PHIPHI(I,J)).GT.1.D-7) THEN
!             PRINT*,'ERROR IN PHIPHI FOR ',I,J,PHIPHI(I,J),TEST(I,J)
           END IF
           IF(DABS(TEST1(I,J)-CHIPHI(J,I)).GT.1.D-7) THEN
!             PRINT*,'ERROR IN CHIPHI FOR ',I,J,CHIPHI(J,I),TEST1(I,J)
           END IF
           IF(DABS(TEST2(I,J)-CHICHI(I,J)).GT.1.D-7) THEN
!             PRINT*,'ERROR IN CHICHI FOR ',I,J,CHICHI(I,J),TEST2(I,J)
           END IF
         ENDDO
       ENDDO
       WRITE(*,FMT='("ENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDEND")')
       RETURN
       END SUBROUTINE TESTME
      END
!
!     .....................................................ORTHO........
      SUBROUTINE WAVES_ORTHO_X_C(NB,OCC,CHICHI,PSIPSI,CHIPSI,LAMBDA)
!     ******************************************************************
!     **                                                              **
!     **  IMPOSES THE ORTHOGONALITY CONSTRAINT ONTO THE ELECTRONS     **
!     **                                                              **
!     **  NEW VERSION WITH DIAGONALIZATION FOR CHIPSI                    **
!     **                                                              **
!     **                                                              **
!     **  THE METHOD IS DESCRIBED IN :                                **
!     **    R.CAR AND M.PARRINELLO, IN "SIMPLE MOLECULAR SYSTEMS      **
!     **    AT VERY HIGH DENSITY", PAGE 455                           **
!     **    ED. A.POLIAN, PLOUBEYRE AND N.BOCCARA                     **
!     **    (PLENUM PUBLISHING CORPORATION,1989)                      **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
!     IMPLICIT NONE
      REAL(8)   ,PARAMETER     :: EPS    = 1.D-8
      REAL(8)   ,PARAMETER     :: DSMALL = 1.D-12
      INTEGER(4),PARAMETER     :: ITERX    = 200
      INTEGER(4),INTENT(IN)    :: NB
      REAL(8)   ,INTENT(IN)    :: OCC(NB)
      COMPLEX(8),INTENT(INOUT) :: LAMBDA(NB,NB)
      COMPLEX(8),INTENT(IN)    :: PSIPSI(NB,NB)
      COMPLEX(8),INTENT(IN)    :: CHIPSI(NB,NB)   !
      COMPLEX(8),INTENT(IN)    :: CHICHI(NB,NB)
      COMPLEX(8),ALLOCATABLE   :: GAMN(:,:)  
      REAL(8)                  :: EIG(NB)
      INTEGER(4)               :: IND,ITER,I,J,K ! RUNNING VARIABLES
      INTEGER(4)               :: IMAX,I0,J0   ! AUXILARY VARIABLES
      REAL(8)                  :: DIGAM,SVAR,FI,FJ,EIGI ! AUXILARY VARIABLES
      COMPLEX(8)               :: HAUX(NB,NB)    
      COMPLEX(8)               :: U(NB,NB)       
      LOGICAL(4)               :: TCONVERGED
      LOGICAL(4)               :: TESSL=.TRUE.
      COMPLEX(8)               :: CSVAR
      REAL(8)                  :: OCCI,OCCJ
      INTEGER(4),EXTERNAL      :: IDAMAX
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_X_C')
      ALLOCATE(GAMN(NB,NB))
!
!     ==================================================================
!     ==  CALCULATE  PSIPSI(I,J)= <PSIBAR(I)|PSIBAR(J)>-1             ==
!     ==        AND  CHIPSI(I,J)   = <PSI0(I)|PSIBAR(J)>              ==
!     ==================================================================
!
!     ==================================================================
!     ==  DIAGONALIZE 0.5*(CHIPSI(I,J)+CHIPSI(J,I))                   ==
!     ==================================================================
      DO I=1,NB
        DO J=1,NB
           GAMN(I,J)=0.5D0*(CHIPSI(I,J)+CONJG(CHIPSI(J,I)))
        ENDDO
      ENDDO
      CALL LIB$DIAGC8(NB,GAMN,EIG,U)
!CALL CDIAG(NB,NB,GAMN,EIG,U)
!WRITE(*,FMT='("EIG",20E10.3)')EIG
!DO I=1,NB
!  WRITE(*,FMT='("U",I2,20E10.3)')I,U(I,:)
!ENDDO

!
!     ==================================================================
!     ==================================================================
!     ==  ITERATIVE CALCULATION OF GAMMA                              ==
!     ==================================================================
!     ==================================================================
      DO ITER=1,ITERX
!PRINTA*,'==================',ITER,'==========================='
!       ================================================================
!       ==  CALCULATE <PHI(+)|PHI(+)>-1 WITH PRESENT LAMBDA           ==
!       ==  GAMN(I,J)=PSIPSI(I,J)+LAMBDA(K,I)*CHIPSI(K,J)             ==
!       ==                       +CHIPSI(K,I)*LAMBDA(K,J)             == 
!       ==           +LAMBDA(K,I)*CHICHI(K,L)*LAMBDA(L,J)-1(I,J)      ==
!       ================================================================
        IF(TESSL) THEN
!         __GAMN(I,J) = CHICHI(I,K)*LAMBDA(K,J)___________________________
          CALL LIB$MATMULC8(NB,NB,NB,CHICHI,LAMBDA,HAUX)
!CALL ZGEMUL(CHICHI,NB,'N',LAMBDA,NB,'N',HAUX,NB,NB,NB,NB)
!         __HAUX(I,J) = HAUX(I,J)+2*CHIPSI(I,J)___________________________
          HAUX=HAUX+2.D0*CHIPSI
!CALL ZAXPY(NB*NB,(2.D0,0.D0),CHIPSI,1,HAUX,1)
!         __GAMN(I,J) = LAMBDA(K,I)*HAUX(K,J)_____________________________
          CALL LIB$SCALARPRODUCTC8(.FALSE.,NB,NB,LAMBDA,NB,HAUX,GAMN)
!CALL ZGEMUL(LAMBDA,NB,'C',HAUX,NB,'N',GAMN,NB,NB,NB,NB)
!         __HAUX(I,J) = HAUX(I,J)+2*CHIPSI(I,J)___________________________
          GAMN(:,:)=GAMN(:,:)+PSIPSI(:,:)
!         __GAMN(I,J) = GAMN(I,J)-1_______________________________________
          DO I=1,NB
            DO J=I,NB
              CSVAR=0.5D0*(GAMN(I,J)+CONJG(GAMN(J,I)))
              GAMN(I,J)=CSVAR
              GAMN(J,I)=CONJG(CSVAR)
            ENDDO
            GAMN(I,I)=GAMN(I,I)-1.D0
          ENDDO
        ELSE
          DO I=1,NB
            DO J=1,NB
              GAMN(I,J)=PSIPSI(I,J)
              IF(I.EQ.J) GAMN(I,J)=GAMN(I,J)-(1.D0,0.D0)
              HAUX(I,J)=CHIPSI(I,J)
              DO K=1,NB
                GAMN(I,J)=GAMN(I,J)+CONJG(CHIPSI(K,I))*LAMBDA(K,J)
                HAUX(I,J)=HAUX(I,J)+CHICHI(I,K)*LAMBDA(K,J)
              ENDDO
            ENDDO
          ENDDO
          DO I=1,NB
            DO J=1,NB
              DO K=1,NB
                GAMN(I,J)=GAMN(I,J)+CONJG(LAMBDA(K,I))*HAUX(K,J)
              ENDDO
            ENDDO
          ENDDO
        END IF
!DO I=1,NB
!WRITE(*,FMT='("A",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       == FIND LARGEST ELEMENT OF THE OVERLAP MATRIX                 ==
!       ================================================================
        DIGAM=0.D0
        TCONVERGED=.TRUE.
        DO I=1,NB
          DO J=I,NB
            IF(ABS(GAMN(I,J)).GT.EPS) THEN
              TCONVERGED=.FALSE.
              DIGAM=MAX(DIGAM,ABS(GAMN(I,J)))
            END IF
          ENDDO
        ENDDO
        IF(TCONVERGED) GOTO 9000
PRINT*,'ITER ',ITER,DIGAM
!
!       ==================================================================
!       ==  OBTAIN CHANGE OF THE LAMBDA MATRIX                          ==
!       ==================================================================
!       == TRANSFORM OVERLAP MATRIX GAMN
        IF(TESSL) THEN
!         ----  HAUX(I,L)=U(K,I)*H0(K,L)
          CALL LIB$SCALARPRODUCTC8(.FALSE.,NB,NB,U,NB,GAMN,HAUX)
!CALL ZGEMUL(U,NB,'C',GAMN,NB,'N',HAUX,NB,NB,NB,NB)
!         ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
          CALL LIB$MATMULC8(NB,NB,NB,HAUX,U,GAMN)
!CALL ZGEMUL(HAUX,NB,'N',U,NB,'N',GAMN,NB,NB,NB,NB)
!
!         ==  MULTIPLY WITH 1/(EIG(I)+EIG(J))
          DO I=1,NB
            EIGI=EIG(I)
            DO J=1,NB
              GAMN(I,J)=GAMN(I,J)/(EIGI+EIG(J))
            ENDDO
          ENDDO
!
!         == TRANSFORM OVERLAP MATRIX GAMN BACK
!         ----  HAUX(I,L)=U(K,I)*H0(K,L)
          CALL LIB$MATMULC8(NB,NB,NB,U,GAMN,HAUX)
!CALL ZGEMUL(U,NB,'N',GAMN,NB,'N',HAUX,NB,NB,NB,NB)
!         ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
          CALL LIB$DYADSUMC8(NB,NB,NB,HAUX,U,GAMN)
!CALL ZGEMUL(HAUX,NB,'N',U,NB,'C',GAMN,NB,NB,NB,NB)
        ELSE
          DO I=1,NB
            DO J=1,NB
              HAUX(I,J)=(0.D0,0.D0)
              DO K=1,NB
                HAUX(I,J)=HAUX(I,J)+CONJG(U(K,I))*GAMN(K,J)
              ENDDO
            ENDDO
          ENDDO
          DO I=1,NB
            DO J=1,NB
              GAMN(I,J)=(0.D0,0.D0)
              DO K=1,NB
                GAMN(I,J)=GAMN(I,J)+HAUX(I,K)*U(K,J)
              ENDDO
              GAMN(I,J)=GAMN(I,J)/(EIG(I)+EIG(J))
            ENDDO
          ENDDO
          DO I=1,NB
            DO J=1,NB
              HAUX(I,J)=(0.D0,0.D0)
              DO K=1,NB
                HAUX(I,J)=HAUX(I,J)+U(I,K)*GAMN(K,J)
              ENDDO
            ENDDO
          ENDDO
          DO I=1,NB
            DO J=1,NB
              GAMN(I,J)=(0.D0,0.D0)
              DO K=1,NB
                GAMN(I,J)=GAMN(I,J)+HAUX(I,K)*CONJG(U(J,K))
              ENDDO
            ENDDO
          ENDDO
        END IF
!
!       ================================================================
!       ==  PROPAGATE GAMMA                                           ==
!       ================================================================
        DO I=1,NB
          OCCI=OCC(I)+DSMALL
          DO J=1,NB
            OCCJ=OCC(J)+DSMALL
            SVAR=2.D0*OCCI/(OCCI+OCCJ)
            LAMBDA(I,J)=LAMBDA(I,J)-SVAR*GAMN(I,J)
          ENDDO
        ENDDO
!DO I=1,NB
!WRITE(*,FMT='("DL",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       == SYMMETRIZE LAMBDA                                          ==
!       ================================================================
        DO I=1,NB
          DO J=1,NB
            IF(OCC(I).LE.OCC(J)) THEN
              LAMBDA(I,J)=CONJG(LAMBDA(J,I))*(OCC(I)+DSMALL)/(OCC(J)+DSMALL)
            END IF
          ENDDO
        ENDDO
!DO I=1,NB
!WRITE(*,FMT='("L",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       ==  ALTERNATIVE                                               ==
!       ================================================================
!       DO I=1,NB
!         FI=F(I)
!         DO J=I+1,NB
!           FJ=F(J)
!           IF(FI+FJ.GT.1.D-6) THEN
!             FI=1.D0
!             FJ=0.D0
!             LAMBDA(J,I)=0.D0
!           ELSE
!             IF(FI.LE.FJ) THEN
!               LAMBDA(I,J)=LAMBDA(J,I)*FI/FJ
!             ELSE IF(FI.GT.FJ) THEN
!               LAMBDA(J,I)=LAMBDA(I,J)*FJ/FI
!             END IF
!           ENDIF
!           SVAR=1.D0/(FI*EIG(I)+FJ*EIG(J))
!           LAMBDA(I,J)=LAMBDA(I,J)-SVAR*GAMN(I,J)*FI
!           LAMBDA(J,I)=LAMBDA(J,I)-SVAR*GAMN(J,I)*FJ
!         ENDDO
!       ENDDO
!
!       DO I=1,NB
!         LAMBDA(I,I)=LAMBDA(I,I)-GAMN(I,I)/(2.D0*EIG(I))
!       ENDDO
!
      ENDDO
      CALL ERROR$MSG('LOOP FOR ORTHOGONALIZATION IS NOT CONVERGED')
      CALL ERROR$STOP('WAVES_ORTHO_X_C')
!
9000  CONTINUE
      DEALLOCATE(GAMN)
                             CALL TRACE$POP
      RETURN
      END
!
!     .....................................................ORTHO........
      SUBROUTINE WAVES_ORTHO_X(NB,OCC,CHICHI,PSIPSI,CHIPSI,LAMBDA)
!     ******************************************************************
!     **                                                              **
!     **  IMPOSES THE ORTHOGONALITY CONSTRAINT ONTO THE ELECTRONS     **
!     **                                                              **
!     **  NEW VERSION WITH DIAGONALIZATION FOR CHIPSI                    **
!     **                                                              **
!     **                                                              **
!     **  THE METHOD IS DESCRIBED IN :                                **
!     **    R.CAR AND M.PARRINELLO, IN "SIMPLE MOLECULAR SYSTEMS      **
!     **    AT VERY HIGH DENSITY", PAGE 455                           **
!     **    ED. A.POLIAN, PLOUBEYRE AND N.BOCCARA                     **
!     **    (PLENUM PUBLISHING CORPORATION,1989)                      **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
!     IMPLICIT NONE
      REAL(8)   ,PARAMETER     :: EPS    = 1.D-8
      REAL(8)   ,PARAMETER     :: DSMALL = 1.D-12
      INTEGER(4),PARAMETER     :: MAX    = 200
      INTEGER(4),INTENT(IN)    :: NB
      REAL(8)   ,INTENT(IN)    :: OCC(NB)
      REAL(8)   ,INTENT(INOUT) :: LAMBDA(NB,NB)
      REAL(8)   ,INTENT(IN)    :: PSIPSI(NB,NB)
      REAL(8)   ,INTENT(IN)    :: CHIPSI(NB,NB)   !
      REAL(8)   ,INTENT(IN)    :: CHICHI(NB,NB)
      REAL(8)   ,ALLOCATABLE   :: GAMN(:,:)  
      REAL(8)                  :: EIG(NB)
      INTEGER(4)               :: IND,ITER,I,J ! RUNNING VARIABLES
      INTEGER(4)               :: IMAX,I0,J0   ! AUXILARY VARIABLES
      REAL(8)                  :: DIGAM,SVAR,FI,FJ,EIGI ! AUXILARY VARIABLES
      REAL(8)                  :: HAUX(NB,NB)    
      REAL(8)                  :: U(NB,NB)       
      REAL(8)                  :: OCCI,OCCJ
      INTEGER(4),EXTERNAL      :: IDAMAX
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_X')
      ALLOCATE(GAMN(NB,NB))
!
!     ==================================================================
!     ==  CALCULATE  PSIPSI(I,J)= <PSIBAR(I)|PSIBAR(J)>-1             ==
!     ==        AND  CHIPSI(I,J)   = <PSI0(I)|PSIBAR(J)>              ==
!     ==================================================================
!
!     ==================================================================
!     ==  DIAGONALIZE 0.5*(CHIPSI(I,J)+CHIPSI(J,I))                   ==
!     ==================================================================
      CALL LIB$DIAGR8(NB,CHIPSI,EIG,U)
!CALL DIAG(NB,NB,CHIPSI,EIG,U)
!WRITE(*,FMT='("EIG",20E10.3)')EIG
!DO I=1,NB
!  WRITE(*,FMT='("U",I2,20E10.3)')I,U(I,:)
!ENDDO
!
!     ==================================================================
!     ==================================================================
!     ==  ITERATIVE CALCULATION OF GAMMA                              ==
!     ==================================================================
!     ==================================================================
      DO ITER=1,MAX
!PRINT*,'==================',ITER,'==========================='
!       ================================================================
!       ==  CALCULATE <PHI(+)|PHI(+)>-1 WITH PRESENT LAMBDA           ==
!       ==  GAMN(I,J)=PSIPSI(I,J)+LAMBDA(K,I)*CHIPSI(K,J)             ==
!       ==                       +CHIPSI(K,I)*LAMBDA(K,J)             == 
!       ==           +LAMBDA(K,I)*CHICHI(K,L)*LAMBDA(L,J)-1(I,J)      ==
!       ================================================================
!       __GAMN(I,J) = CHICHI(I,K)*LAMBDA(K,J)___________________________
        CALL LIB$MATMULR8(NB,NB,NB,CHICHI,LAMBDA,HAUX)
!CALL DGEMUL(CHICHI,NB,'N',LAMBDA,NB,'N',HAUX,NB,NB,NB,NB)
!       __HAUX(I,J) = HAUX(I,J)+2*CHIPSI(I,J)___________________________
        HAUX=HAUX+2.D0*CHIPSI
!CALL DAXPY(NB*NB,2.D0,CHIPSI,1,HAUX,1)
!       __GAMN(I,J) = LAMBDA(K,I)*HAUX(K,J)_____________________________
        CALL LIB$SCALARPRODUCTR8(.FALSE.,NB,NB,LAMBDA,NB,HAUX,GAMN)
!CALL DGEMUL(LAMBDA,NB,'T',HAUX,NB,'N',GAMN,NB,NB,NB,NB)
!       __GAMN(I,J) = GAMN(I,J)-1_______________________________________
        DO I=1,NB
          DO J=I,NB
            SVAR=0.5D0*(GAMN(I,J)+GAMN(J,I))+PSIPSI(I,J)
            GAMN(J,I)=SVAR
            GAMN(I,J)=SVAR
          ENDDO
          GAMN(I,I)=GAMN(I,I)-1.D0
        ENDDO
!DO I=1,NB
!WRITE(*,FMT='("A",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       == FIND LARGEST ELEMENT OF THE OVERLAP MATRIX                 ==
!       ================================================================
        IMAX=IDAMAX(NB*NB,GAMN(1,1),1)
        J0=(IMAX-1)/NB+1
        I0=IMAX-(J0-1)*NB
        DIGAM=DABS(GAMN(I0,J0))
!       PRINT*,'ITER ',ITER,I0,J0,DIGAM,NCON
        IF(DIGAM.LT.EPS) GOTO 9000
PRINT*,'ITER ',ITER,DIGAM
!
!       ==================================================================
!       ==  OBTAIN CHANGE OF THE LAMBDA MATRIX                          ==
!       ==================================================================
!       == TRANSFORM OVERLAP MATRIX GAMN
!       ----  HAUX(I,L)=U(K,I)*H0(K,L)
        CALL LIB$SCALARPRODUCTR8(.FALSE.,NB,NB,U,NB,GAMN,HAUX)
!CALL DGEMUL(U,NB,'T',GAMN,NB,'N',HAUX,NB,NB,NB,NB)
!       ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
        CALL LIB$MATMULR8(NB,NB,NB,HAUX,U,GAMN)
!CALL DGEMUL(HAUX,NB,'N',U,NB,'N',GAMN,NB,NB,NB,NB)
!
!       ==  MULTIPLY WITH 1/(EIG(I)+EIG(J))
        DO I=1,NB
          EIGI=EIG(I)
          DO J=1,NB
            GAMN(I,J)=GAMN(I,J)/(EIGI+EIG(J))
          ENDDO
        ENDDO
!
!       == TRANSFORM OVERLAP MATRIX GAMN BACK
!       ----  HAUX(I,L)=U(K,I)*H0(K,L)
        CALL LIB$MATMULR8(NB,NB,NB,U,GAMN,HAUX)
!CALL DGEMUL(U,NB,'N',GAMN,NB,'N',HAUX,NB,NB,NB,NB)
!       ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
        CALL LIB$DYADSUMR8(NB,NB,NB,HAUX,U,GAMN)
!CALL DGEMUL(HAUX,NB,'N',U,NB,'T',GAMN,NB,NB,NB,NB)
!DO I=1,NB
!WRITE(*,FMT='("DL",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       ==  PROPAGATE GAMMA                                           ==
!       ================================================================
        DO I=1,NB
          OCCI=OCC(I)+DSMALL
          DO J=1,NB
            OCCJ=OCC(J)+DSMALL
            SVAR=2.D0*OCCI/(OCCI+OCCJ)
            LAMBDA(I,J)=LAMBDA(I,J)-SVAR*GAMN(I,J)
          ENDDO
        ENDDO
!DO I=1,NB
!WRITE(*,FMT='("L",I2,20E10.3)')I,GAMN(I,:)
!ENDDO
!
!       ================================================================
!       == SYMMETRIZE LAMBDA                                          ==
!       ================================================================
        DO I=1,NB
          DO J=1,NB
            IF(OCC(I).LT.OCC(J)) THEN
              LAMBDA(I,J)=LAMBDA(J,I)*(OCC(I)+DSMALL)/(OCC(J)+DSMALL)
            END IF
          ENDDO
        ENDDO
!
!       ================================================================
!       ==  ALTERNATIVE                                               ==
!       ================================================================
!       DO I=1,NB
!         FI=F(I)
!         DO J=I+1,NB
!           FJ=F(J)
!           IF(FI+FJ.GT.1.D-6) THEN
!             FI=1.D0
!             FJ=0.D0
!             LAMBDA(J,I)=0.D0
!           ELSE
!             IF(FI.LE.FJ) THEN
!               LAMBDA(I,J)=LAMBDA(J,I)*FI/FJ
!             ELSE IF(FI.GT.FJ) THEN
!               LAMBDA(J,I)=LAMBDA(I,J)*FJ/FI
!             END IF
!           ENDIF
!           SVAR=1.D0/(FI*EIG(I)+FJ*EIG(J))
!           LAMBDA(I,J)=LAMBDA(I,J)-SVAR*GAMN(I,J)*FI
!           LAMBDA(J,I)=LAMBDA(J,I)-SVAR*GAMN(J,I)*FJ
!         ENDDO
!       ENDDO
!
!       DO I=1,NB
!         LAMBDA(I,I)=LAMBDA(I,I)-GAMN(I,I)/(2.D0*EIG(I))
!       ENDDO
!
      ENDDO
!      PRINT*,'EIG ',EIG
!      PRINT*,'OCC ',OCC
!      DO I=1,NB
!        DO J=1,NB
!          WRITE(*,FMT='(2I3,7F10.5)')I,J,LAMBDA(I,J),U(I,J),PSIPSI(I,J) &
!    &                               ,CHIPSI(I,J),CHICHI(I,J)
!        ENDDO
!      ENDDO
      CALL ERROR$MSG('LOOP FOR ORTHOGONALIZATION IS NOT CONVERGED')
      CALL ERROR$STOP('WAVES_ORTHO_X')
      
!
9000  CONTINUE
      DEALLOCATE(GAMN)
                             CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_GRAMSCHMIDT(MAP,GSET,NAT,R,NGL,NDIM,NBH,NB,PSI)
!     ******************************************************************
!     **                                                              **
!     **  GRAM-SCHMIDT ORTHOGONALIZATION OF A SET OF WAVE FUNCTIONS   **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY: MAP_TYPE,GSET_TYPE
      IMPLICIT NONE
      TYPE(MAP_TYPE)  ,INTENT(IN) :: MAP
      TYPE(GSET_TYPE) ,INTENT(IN) :: GSET
      INTEGER(4)      ,INTENT(IN) :: NAT
      REAL(8)         ,INTENT(IN) :: R(3,NAT)
      INTEGER(4)      ,INTENT(IN) :: NGL
      INTEGER(4)      ,INTENT(IN) :: NDIM
      INTEGER(4)      ,INTENT(IN) :: NBH
      INTEGER(4)      ,INTENT(IN) :: NB
      COMPLEX(8)      ,INTENT(INOUT):: PSI(NGL,NDIM,NBH)
      INTEGER(4)                  :: NPRO
      INTEGER(4)                  :: IDIM,IG,I,J
      INTEGER(4)                  :: IBH1,IBH2,IB1A,IB1B,IB2A,IB2B
      COMPLEX(8)                  :: CSVAR
      COMPLEX(8)      ,ALLOCATABLE:: X(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: PROJ(:,:,:)
      COMPLEX(8)      ,ALLOCATABLE:: OVERLAP(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: AUXMAT(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: X1(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: X2(:,:)
      COMPLEX(8)      ,ALLOCATABLE:: PSIINV(:,:,:)
      COMPLEX(8)      ,ALLOCATABLE:: CWORK(:,:,:)
      COMPLEX(8)                  :: CSVAR1,CSVAR2
      COMPLEX(8)      ,PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL(4)                  :: TINV
      REAL(8)                     :: NORM(NB),SVAR
      COMPLEX(8)                  :: XTWOBYTWO(2,2)
      LOGICAL(4)      ,PARAMETER  :: TTEST=.false.
      INTEGER(4)      ,ALLOCATABLE:: SMAP(:)
!     ******************************************************************
!
!     ==================================================================
!     ==  CALCULATE PROJECTIONS FOR THE NEW POSITIONS                 ==
!     ==================================================================
      CALL PLANEWAVE$SELECT(GSET%ID)
      TINV=GSET%TINV
      NPRO=MAP%NPRO
      ALLOCATE(PROJ(NDIM,NBH,NPRO))
      CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,NPRO,PSI,PROJ)
!     
!     ==================================================================
!     ==  OVERLAP OF <PSI0|PSI0>,                                     ==
!     ==================================================================
      ALLOCATE(OVERLAP(NB,NB))
      ALLOCATE(AUXMAT(NB,NB))
      CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,PSI,PSI,OVERLAP)
      CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIM,NBH,NB,NPRO,PROJ,PROJ,AUXMAT)
      DO J=1,NB
        DO I=1,NB
          OVERLAP(I,J)=OVERLAP(I,J)+AUXMAT(I,J)
        ENDDO
      ENDDO
      DEALLOCATE(AUXMAT)
      DEALLOCATE(PROJ)
!     
!     =================================================================
!     ==  OBTAIN ORTHOGONALIZATION TRANSFORM                         ==
!     =================================================================
      ALLOCATE(SMAP(NB))
      DO I=1,NB
        SMAP(I)=I
      ENDDO
      ALLOCATE(X(NB,NB))
      CALL WAVES_ORTHO_Y_C(NB,OVERLAP,OVERLAP,OVERLAP,X,SMAP)
      DEALLOCATE(OVERLAP)
      DEALLOCATE(SMAP)
!     
!     =================================================================
!     ==  TRANSFORM WAVE FUNCTIONS  |PSI>=|PSI>X                     ==
!     =================================================================
      IF(TINV) THEN
        ALLOCATE(X1(NBH,NBH))
        ALLOCATE(X2(NBH,NBH))
!
!       == FIRST FOLD DOWN OVERLAP MATRIX TO SUPER WAVE FUNCTIONS ======
        DO IBH1=1,NBH
          IB1A=2*IBH1-1
          IB1B=MIN(2*IBH1,NB)
          DO IBH2=1,NBH
            IB2A=2*IBH2-1
            IB2B=MIN(2*IBH2,NB)
            XTWOBYTWO(:,:)=(0.D0,0.D0)
            DO I=IB1A,IB1B
              DO J=IB2A,IB2B
                XTWOBYTWO(I-IB1A+1,J-IB2A+1)=X(I,J)
              ENDDO
            ENDDO
            CSVAR1=    XTWOBYTWO(1,1)+CI*XTWOBYTWO(1,2)
            CSVAR2=-CI*XTWOBYTWO(2,1)   +XTWOBYTWO(2,2)
            X1(IBH1,IBH2)=0.5D0*(CSVAR1+CSVAR2)
            X2(IBH1,IBH2)=0.5D0*(CSVAR1-CSVAR2)
          ENDDO
        ENDDO
        DEALLOCATE(X)
        ALLOCATE(PSIINV(NGL,NDIM,NBH))
        PSIINV=PSI
        CALL PLANEWAVE$ADDPRODUCT(' ',NGL,NDIM,NBH,PSI,NBH,PSIINV,X1)
        CALL PLANEWAVE$ADDPRODUCT('-',NGL,NDIM,NBH,PSI,NBH,PSIINV,X2)
        DEALLOCATE(PSIINV)
!        DO I=NBH,1,-1
!          DO J=1,I       !WORKS ONLY FOR TRIANGULAR X
!            DO IDIM=1,NDIM
!              DO IG=1,NGL
!                PSI(IG,IDIM,I)=PSI(IG,IDIM,I) &
!     &                        +PSI(IG,IDIM,J)   *X1(J,I) &
!     &                        +PSIINV(IG,IDIM,J)*X2(J,I)
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDDO         
        DEALLOCATE(X1)
        DEALLOCATE(X2)
      ELSE
        ALLOCATE(PSIINV(NGL,NDIM,NB))
        PSIINV=PSI
        CALL PLANEWAVE$ADDPRODUCT(' ',NGL,NDIM,NB,PSI,NB,PSIINV,X)
        DEALLOCATE(PSIINV)
!        DO I=NB,1,-1
!          DO J=1,I
!            DO IDIM=1,NDIM
!              DO IG=1,NGL
!                PSI(IG,IDIM,I)=PSI(IG,IDIM,I)+PSI(IG,IDIM,J)*X(J,I)
!              ENDDO
!            ENDDO
!          ENDDO
!       ENDDO         
        DEALLOCATE(X)
      ENDIF
! 
!     =================================================================
!     ==                                                             ==
!     =================================================================
      IF(TTEST) THEN
        ALLOCATE(PROJ(NDIM,NBH,NPRO))
        CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,NPRO,PSI,PROJ)
        ALLOCATE(AUXMAT(NB,NB))
        CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIM,NBH,NB,NPRO,PROJ,PROJ,AUXMAT)
        DEALLOCATE(PROJ)
        ALLOCATE(OVERLAP(NB,NB))
        CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH,NB,PSI,PSI,OVERLAP)
        DO J=1,NB
          DO I=1,NB
            OVERLAP(I,J)=OVERLAP(I,J)+AUXMAT(I,J)
          ENDDO
        ENDDO
        DEALLOCATE(AUXMAT)
        ALLOCATE(X(NB,NB))
        CALL LIB$DIAGC8(NB,OVERLAP,NORM,X)
!CALL CDIAG(NB,NB,OVERLAP,NORM,X)
        DEALLOCATE(OVERLAP)
        DO I=1,NB
          IF(ABS(NORM(I)-1.D0).GT.1.D-4) THEN
            CALL ERROR$MSG('GRAM-SCHMIDT ORTHOGONALIZATION FAILED')
            CALL ERROR$I4VAL('STATE',I)
            CALL ERROR$R8VAL('NORM',NORM(I))
            CALL ERROR$STOP('WAVES_GRAMSCHMIDT')
          END IF
        ENDDO
        DEALLOCATE(X)
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SWITCH()
!     ******************************************************************
!     **  STEP FORWARD                                                **
!     **                                                              **
!     **  REMARKS: THE PROPAGATION MUST BE CONSISTENT WITH THAT IN    **
!     **   ORTHOGONALIZE                                              **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: IKPT
      INTEGER(4)             :: ISPIN
      COMPLEX(8),POINTER     :: TPSI(:,:,:)
      LOGICAL(4)             :: TSTRESS
      REAL(8)                :: CELLSCALE
      REAL(8)                :: MAPTOCELL(3,3)
      INTEGER(4)             :: NLAMBDA
      INTEGER(4),PARAMETER   :: NLAMBDAX=2
      INTEGER(4)             :: NB
      COMPLEX(8),ALLOCATABLE :: RLAMP(:,:)
      COMPLEX(8),POINTER     :: RLAMPTR(:,:)
!     ******************************************************************
      CALL CELL$GETL4('MOVE',TSTRESS)
      WAVEEKIN1=0.D0
      WAVEEKIN2=0.D0
      IF(TSTRESS) THEN
        CALL CELL$GETR8A('MAPTOCELL',9,MAPTOCELL)
        CELLSCALE=MAPTOCELL(1,1)*(MAPTOCELL(2,2)*MAPTOCELL(3,3)  &
     &                           -MAPTOCELL(3,2)*MAPTOCELL(2,3)) &
     &           +MAPTOCELL(2,1)*(MAPTOCELL(3,2)*MAPTOCELL(1,3)  &
     &                           -MAPTOCELL(1,2)*MAPTOCELL(3,3)) &
     &           +MAPTOCELL(3,1)*(MAPTOCELL(1,2)*MAPTOCELL(2,3)  &
     &                           -MAPTOCELL(2,2)*MAPTOCELL(1,3))
        CELLSCALE=1.D0/SQRT(CELLSCALE)
      ENDIF
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          TPSI=>THIS%PSI0
          THIS%PSI0=>THIS%PSIM
          THIS%PSIM=>TPSI
          IF(TSTRESS) THEN
            THIS%PSI0=THIS%PSI0*CELLSCALE
            THIS%PSIM=THIS%PSIM*CELLSCALE
          END IF
!
!         ==============================================================
!         == EXTRAPOLATE LAGRANGE MULTIPLIERS                         ==
!         ==============================================================
          NLAMBDA=4
          IF(.NOT.ASSOCIATED(THIS%RLAM3M)) NLAMBDA=3
          IF(.NOT.ASSOCIATED(THIS%RLAM2M)) NLAMBDA=2
          IF(.NOT.ASSOCIATED(THIS%RLAMM)) NLAMBDA=1
          IF(.NOT.ASSOCIATED(THIS%RLAM0)) NLAMBDA=0
          NLAMBDA=MIN(NLAMBDA,NLAMBDAX)
          NB=THIS%NB
          ALLOCATE(RLAMP(NB,NB))
          IF(NLAMBDA.EQ.0) THEN
            RLAMP=0.D0
          ELSE IF(NLAMBDA.EQ.1) THEN
            RLAMP=THIS%RLAM0
            THIS%RLAMM=>THIS%RLAM0
            NULLIFY(THIS%RLAM0)
            IF(NLAMBDA.EQ.NLAMBDAX) THEN
              THIS%RLAM0=>THIS%RLAMM
              NULLIFY(THIS%RLAMM)
            END IF 
          ELSE IF(NLAMBDA.EQ.2) THEN
            RLAMP=2.D0*THIS%RLAM0-THIS%RLAMM
            THIS%RLAM2M=>THIS%RLAMM
            THIS%RLAMM=>THIS%RLAM0
            NULLIFY(THIS%RLAM0)
            IF(NLAMBDA.EQ.NLAMBDAX) THEN  ! REUSE MEMORY 
              THIS%RLAM0=>THIS%RLAM2M     
              NULLIFY(THIS%RLAM2M)        
            END IF 
          ELSE IF(NLAMBDA.EQ.3) THEN
            RLAMP=3.D0*THIS%RLAM0-3.D0*THIS%RLAMM+THIS%RLAM2M
            THIS%RLAM3M=>THIS%RLAM2M
            THIS%RLAM2M=>THIS%RLAMM
            THIS%RLAMM=>THIS%RLAM0
            NULLIFY(THIS%RLAM0)
            IF(NLAMBDA.EQ.NLAMBDAX) THEN
              THIS%RLAM0=>THIS%RLAM3M
              NULLIFY(THIS%RLAM3M)
            END IF 
          ELSE IF(NLAMBDA.EQ.4) THEN
            RLAMP=4.D0*THIS%RLAM0-6.D0*THIS%RLAMM+4.D0*THIS%RLAM2M-THIS%RLAM3M
            RLAMPTR=>THIS%RLAM3M
            THIS%RLAM3M=>THIS%RLAM2M
            THIS%RLAM2M=>THIS%RLAMM
            THIS%RLAMM=>THIS%RLAM0
            THIS%RLAM0=>RLAMPTR
            NULLIFY(RLAMPTR)
          END IF
          IF(.NOT.ASSOCIATED(THIS%RLAM0))ALLOCATE(THIS%RLAM0(NB,NB))
          THIS%RLAM0=RLAMP             
          DEALLOCATE(RLAMP)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DELETE HAMILTONIAN AND EIGENVALUES                          ==
!     ==================================================================
      IF(THAMILTON) THEN
        DO IKPT=1,NKPT
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            CALL PLANEWAVE$SELECT(GSET%ID)
            DEALLOCATE(THIS%EIGVAL)         
            DEALLOCATE(THIS%EIGVEC)         
            DEALLOCATE(THIS%EXPECTVAL)
          ENDDO
        ENDDO
        THAMILTON=.FALSE.
      END IF
      RETURN
      END
!
!     ......................................................RANWAV......
      SUBROUTINE WAVES_RANDOMIZE(NG,NDIM,NB,AMPLITUDE,G2,PSI)
!     ******************************************************************
!     **                                                              **
!     **  CREATE RANDOM WAVE FUNCTIONS                                **
!     **                                                              **
!     **  THE MAXIMUM WEIGHT OF THE WAVE FUNCTIONS AT EPW[RY]=GC2     **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB              ! #(BANDS)
      INTEGER(4),INTENT(IN)    :: NG              ! #(PLANE WAVES),MAX
      INTEGER(4),INTENT(IN)    :: NDIM            ! #(PLANE WAVES),MAX
      REAL(8)   ,INTENT(IN)    :: AMPLITUDE       ! SCALE FACTOR
      REAL(8)   ,INTENT(IN)    :: G2(NG)          ! G**2
      COMPLEX(8),INTENT(INOUT) :: PSI(NG,NDIM,NB) ! PS-WAVE FUNCTION
      INTEGER(4)               :: IB,IG,IDIM
      REAL(8)                  :: PI,GC,FAC
      REAL(8)   ,PARAMETER     :: GC2=10.D0
      REAL(8)                  :: SCALE(NG)
      REAL(8)                  :: REC,RIM
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      FAC=2.D0*SQRT(PI*GC2)
      FAC=FAC**3/REAL(NDIM,8)*(2.D0/3.D0)
      FAC=AMPLITUDE/FAC
      DO IG=1,NG
        SCALE(IG)=FAC*EXP(-0.5D0*G2(IG)/GC2)
      ENDDO
      CALL LIB$RANDOMSEED
      DO IB=1,NB
        DO IDIM=1,NDIM
          DO IG=1,NG
            CALL LIB$RANDOM(REC)
            CALL LIB$RANDOM(RIM)
            REC=2.D0*REC-1.D0
            RIM=2.D0*RIM-1.D0
            PSI(IG,IDIM,IB)=PSI(IG,IDIM,IB)+CMPLX(REC,RIM)*SCALE(IG)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ......................................................RANWAV......
      SUBROUTINE WAVES_INITIALIZERANDOM(NG,NDIM,NB,G2,PSI)
!     ******************************************************************
!     **                                                              **
!     **  CREATE RANDOM WAVE FUNCTIONS                                **
!     **                                                              **
!     **  THE MAXIMUM WEIGHT OF THE WAVE FUNCTIONS AT EPW[RY]=GC2     **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB              ! #(BANDS)
      INTEGER(4),INTENT(IN)    :: NG              ! #(PLANE WAVES),MAX
      INTEGER(4),INTENT(IN)    :: NDIM            ! #(PLANE WAVES),MAX
      REAL(8)   ,INTENT(IN)    :: G2(NG)          ! G**2
      COMPLEX(8),INTENT(OUT)   :: PSI(NG,NDIM,NB) ! PS-WAVE FUNCTION
      INTEGER(4)               :: IB,IG,IDIM
      REAL(8)                  :: PI,GC,FAC
      REAL(8)   ,PARAMETER     :: GC2=10.D0
      REAL(8)                  :: SCALE(NG)
      REAL(8)                  :: REC,RIM
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      FAC=2.D0*SQRT(PI*GC2)
      FAC=FAC**3/REAL(NDIM,8)*(2.D0/3.D0)
      FAC=1.D0/FAC
      DO IG=1,NG
        SCALE(IG)=FAC*EXP(-0.5D0*G2(IG)/GC2)
      ENDDO
      CALL LIB$RANDOMSEED
      DO IB=1,NB
        DO IDIM=1,NDIM
          DO IG=1,NG
            CALL LIB$RANDOM(REC)
            CALL LIB$RANDOM(RIM)
            REC=2.D0*REC-1.D0
            RIM=2.D0*RIM-1.D0
            PSI(IG,IDIM,IB)=CMPLX(REC,RIM)*SCALE(IG)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!======================================================================
!======================================================================
!======================================================================
!======================================================================
!======================================================================
!======================================================================
!
!     ..................................................................
      SUBROUTINE WAVES$REPORT(NFIL)
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)   ,PARAMETER  :: MBYTE=2.D0**20
      INTEGER(4)            :: IKPT
      REAL(8)               :: RY
      REAL(8)               :: SVAR
      REAL(8)               :: MEMORY
      INTEGER(4)            :: NG
!     ******************************************************************
      CALL CONSTANTS('RY',RY)
      CALL REPORT$TITLE(NFIL,'WAVE FUNCTIONS')
      CALL WAVES_SELECTWV(1,1)
      CALL REPORT$I4VAL(NFIL,'NUMBER OF BANDS',THIS%NB,' ')

      CALL REPORT$I4VAL(NFIL,'NUMBER OF K-POINTS',NKPT,' ')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF SPINS',NSPIN,' ')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF SPINOR COMPONENTS',NDIM,' ')
      CALL REPORT$R8VAL(NFIL,'PLANE WAVE CUTOFF',EPWPSI/RY,'RY')
      CALL REPORT$R8VAL(NFIL,'WAVE FUNCTION MASS',EMASS,'A.U.')
      CALL REPORT$R8VAL(NFIL,'G**2 ENHANCEMENT OF FUNCTION MASS',EMASSCG2,' ')
      CALL REPORT$R8VAL(NFIL,'BUCKET POTENTIAL STARTS AT 0.5G^2=',EPWPSI0/RY,'RY')
      CALL REPORT$R8VAL(NFIL,'BUCKET POTENTIAL PREFACTOR',D2EPWPSI,'H')
      IF(ANNEE.NE.0) THEN
        CALL REPORT$R8VAL(NFIL,'FRICTION',ANNEE,' ')
      END IF
      IF(TSTOP) THEN
        CALL REPORT$CHVAL(NFIL,'INITIAL VELOCITY IS SET TO','ZERO')
      END IF
!     IF(TRANDOMIZE) THEN
!       CALL REPORT$R8VAL(NFIL &
!    &      ,'INITIAL VELOCITIES ARE RANDOMIZED WITH ENERGY',AMPRE,'H')
!     END IF
      IF(.NOT.TSAFEORTHO) THEN
        WRITE(NFIL,FMT='("EIGENSTATES ARE CALCULATED."' &
     &                 //'," (NO STRICT ENERGY CONSERVATION)")')
      END IF
!     
!     ================================================================
!     ==  REPORT INFORMATION ABOUT G-VECTORS                        ==
!     ================================================================
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        NG=GSET%NGL
        CALL MPE$COMBINE('+',NG)
        CALL REPORT$I4VAL(NFIL &
     &       ,'NUMBER OF (GLOBAL) PLANE WAVES FOR WAVE FUNCTION',NG,' ')
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$WRITEPDOS
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE WAVES_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: NFIL
      REAL(8)   ,ALLOCATABLE :: VAL(:,:)
      REAL(8)   ,ALLOCATABLE :: DER(:,:)
      REAL(8)   ,ALLOCATABLE :: OV(:,:,:)
      REAL(8)   ,ALLOCATABLE :: R(:,:)
      REAL(8)   ,ALLOCATABLE :: RAD(:)
      REAL(8)   ,ALLOCATABLE :: WORK(:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)                :: AEZ
      INTEGER(4),ALLOCATABLE :: IZ(:)
      INTEGER(4)             :: NAT,NSP
      REAL(8)                :: R1,DEX,XEXP,RI
      INTEGER(4)             :: NR
      COMPLEX(8)             :: CSVAR
      INTEGER(4)             :: NB,NBH
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      COMPLEX(8),ALLOCATABLE :: VEC(:,:)      
      COMPLEX(8),ALLOCATABLE :: VECTOR1(:,:,:)
      INTEGER(4)             :: LNXX,LNX,NPRO
      INTEGER(4)             :: NGL
      INTEGER(4),ALLOCATABLE :: LOX(:)
      LOGICAL(4)             :: TINV
      REAL(8)                :: RBAS(3,3)
      INTEGER(4)             :: ISP,IR
      INTEGER(4)             :: NTASKS,THISTASK
      INTEGER(4)             :: IB1,IB2,IBH,LN1,LN2,IDIM,IKPT,ISPIN,IPRO
      COMPLEX(8),ALLOCATABLE :: PROJ(:,:,:)
      CHARACTER(16),ALLOCATABLE :: ATOMID(:)
      REAL(8)                :: SVAR      
      CHARACTER(32)          :: FLAG='011004'
      INTEGER(4)             :: NBX
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$WRITEPDOS')
      IF(.NOT.THAMILTON) THEN
        CALL ERROR$MSG('EIGENVALUES NOT PRESENT')
        CALL ERROR$STOP('WAVES$WRITEPDOS')
      END IF
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     == GENERAL QUANTITIES                                           ==
!     ==================================================================
      NSP=MAP%NSP
      NAT=MAP%NAT
      NPRO=MAP%NPRO
      LNXX=MAP%LNXX
      CALL PDOS$SETI4('NAT',NAT)
      CALL PDOS$SETI4('NSP',NSP)
      CALL PDOS$SETI4('NKPT',NKPT)
      CALL PDOS$SETI4('NSPIN',NSPIN)
      CALL PDOS$SETI4('NDIM',NDIM)
      CALL PDOS$SETI4('NPRO',NPRO)
      CALL PDOS$SETI4('LNXX',LNXX)
      CALL PDOS$SETI4A('LNX',NSP,MAP%LNX)
      CALL PDOS$SETI4A('LOX',LNXX*NSP,MAP%LOX)
      CALL PDOS$SETI4A('ISPECIES',NAT,MAP%ISP)
      IF(THISTASK.EQ.1) THEN
        CALL FILEHANDLER$UNIT('PDOS',NFIL)
        REWIND NFIL
        WRITE(NFIL)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,FLAG
        WRITE(NFIL)MAP%LNX(:),MAP%LOX(:,:),MAP%ISP(:)
      END IF
!
!     ==================================================================
!     == ATOMIC STRUCTURE                                             ==
!     ==================================================================
      ALLOCATE(R(3,NAT))
      ALLOCATE(ATOMID(NAT))
      CALL CELL$GETR8A('T(0)',9,RBAS)
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R)
      CALL ATOMLIST$GETCHA('NAME',0,NAT,ATOMID)
      CALL PDOS$SETR8A('RBAS',9,RBAS)
      CALL PDOS$SETR8A('R',3*NAT,R)
!     CALL PDOS$SETCHA('ATOMID',NAT,ATOMID)
      IF(THISTASK.EQ.1) THEN
        WRITE(NFIL)RBAS,R,ATOMID
      END IF
      DEALLOCATE(ATOMID)
!
!     ==================================================================
!     == ELEMENT SPECIFIC QUANTITIES                                  ==
!     ==================================================================
      ALLOCATE(VAL(LNXX,NSP))
      ALLOCATE(DER(LNXX,NSP))
      ALLOCATE(OV(LNXX,LNXX,NSP))
      ALLOCATE(LOX(LNXX))
      ALLOCATE(IZ(NSP))
      ALLOCATE(RAD(NSP))
      CALL SETUP$GETI4('NR',NR)
      ALLOCATE(AEPHI(NR,LNXX))
      ALLOCATE(WORK(NR))
      CALL SETUP$GETR8('R1',R1)
      CALL SETUP$GETR8('DEX',DEX)
      XEXP=EXP(DEX)
      DO ISP=1,NSP
        OV(:,:,:)=0.D0
        CALL SETUP$ISELECT(ISP)
        LNX=MAP%LNX(ISP)
        LOX=MAP%LOX(:,ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        IZ(ISP)=NINT(AEZ)
        CALL PERIODICTABLE$GET(IZ(ISP),'R(ASA)',RAD(ISP))
        CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,AEPHI)
        DO LN1=1,LNX
          CALL RADIAL$VALUE(R1,DEX,NR,AEPHI(1,LN1),RAD(ISP),VAL(LN1,ISP))
          CALL RADIAL$DERIVATIVE(R1,DEX,NR,AEPHI(1,LN1),RAD(ISP),DER(LN1,ISP))
          DO LN2=LN1,LNX
            IF(LOX(LN1).NE.LOX(LN2)) CYCLE
            RI=R1/XEXP
            DO IR=1,NR
              RI=RI*XEXP
              WORK(IR)=RI**2*AEPHI(IR,LN1)*AEPHI(IR,LN2)
            ENDDO
            CALL RADIAL$INTEGRAL1(R1,DEX,NR,WORK,RAD(ISP),OV(LN1,LN2,ISP))
            OV(LN2,LN1,ISP)=OV(LN1,LN2,ISP)
          ENDDO
        ENDDO
        IF(THISTASK.EQ.1) THEN
          WRITE(NFIL)IZ(ISP),RAD(ISP),VAL(1:LNX,ISP),DER(1:LNX,ISP) &
     &              ,OV(1:LNX,1:LNX,ISP)
        END IF
      ENDDO
      CALL PDOS$SETI4A('IZ',NSP,IZ)
      CALL PDOS$SETR8A('RAD',NSP,RAD)
      CALL PDOS$SETR8A('PHI',LNXX*NSP,VAL)
      CALL PDOS$SETR8A('DPHIDR',LNXX*NSP,DER)
      CALL PDOS$SETR8A('OVERLAP',LNXX*LNXX*NSP,OV)
      DEALLOCATE(VAL)
      DEALLOCATE(IZ)
      DEALLOCATE(RAD)
      DEALLOCATE(DER)
      DEALLOCATE(OV)
      DEALLOCATE(WORK)
      DEALLOCATE(AEPHI)
      DEALLOCATE(LOX)
!
!     ==================================================================
!     ==  NOW WRITE PROJECTIONS                                       ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NBX*NKPT*NSPIN,OCC)
      ALLOCATE(XK(3,NKPT))
      CALL DYNOCC$GETR8A('XK',3*NKPT,XK)
      CALL PDOS$SETR8A('XK',3*NKPT,XK)
      ALLOCATE(VEC(NDIM,NPRO))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(.NOT.ASSOCIATED(THIS%EIGVEC)) THEN
            CALL ERROR$MSG('EIGENVALUES NOT PRESENT')
            CALL ERROR$STOP('WAVES$WRITEPDOS')
          END IF
          NB=THIS%NB
          NBH=THIS%NBH
          NGL=GSET%NGL
          ALLOCATE(VECTOR1(NDIM,NPRO,NB))
          TINV=GSET%TINV
          IF(THISTASK.EQ.1) THEN
            WRITE(NFIL)XK(:,IKPT),NB
          END IF
!
!         =============================================================
!         ==  CALCULATE PROJECTIONS                                  ==
!         =============================================================
          ALLOCATE(PROJ(NDIM,NBH,NPRO))
          CALL WAVES_PROJECTIONS(MAP,GSET,NAT,R,NGL,NDIM,NBH,NPRO &
     &                            ,THIS%PSI0,PROJ)
!
!         ==============================================================
!         == UNRAVEL SUPER WAVE FUNCTIONS                             ==
!         ==============================================================
          IF(TINV) THEN
            DO IBH=1,NBH
              IB1=1+2*(IBH-1)
              IB2=2+2*(IBH-1)
              DO IPRO=1,NPRO
                DO IDIM=1,NDIM
                  VECTOR1(IDIM,IPRO,IB1)=REAL(PROJ(IDIM,IBH,IPRO))
                  VECTOR1(IDIM,IPRO,IB2)=AIMAG(PROJ(IDIM,IBH,IPRO))
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO IB1=1,NB
              DO IPRO=1,NPRO
                DO IDIM=1,NDIM
                  VECTOR1(IDIM,IPRO,IB1)=PROJ(IDIM,IB1,IPRO)
                ENDDO
              ENDDO
            ENDDO
          END IF
          DEALLOCATE(PROJ)
!
!         ==============================================================
!         ==  TRANSFORM TO EIGENSTATES IF SAFEORTHO=.TRUE.,           ==
!         ==  AND WRITE                                               ==
!         ==============================================================
          IF(TSAFEORTHO) THEN
            DO IB1=1,NB
              VEC=0.D0
              DO IB2=1,NB
                CSVAR=THIS%EIGVEC(IB2,IB1)
                DO IPRO=1,NPRO
                  DO IDIM=1,NDIM
                    VEC(IDIM,IPRO)=VEC(IDIM,IPRO)+VECTOR1(IDIM,IPRO,IB2)*CSVAR
                  ENDDO
                ENDDO
              ENDDO
              IF(THISTASK.EQ.1) THEN
                IF(FLAG.EQ.'011004') THEN
                  WRITE(NFIL)THIS%EIGVAL(IB1),OCC(IB1,IKPT,ISPIN),VEC
                ELSE
                  WRITE(NFIL)THIS%EIGVAL(IB1),VEC
                END IF
!PRINT*,'E ',THIS%EIGVAL(IB1)
!WRITE(*,FMT='(10F10.5)')VEC
              END IF
            ENDDO
          ELSE
            DO IB1=1,NB
              IF(THISTASK.EQ.1) THEN
                IF(FLAG.EQ.'011004') THEN
                  WRITE(NFIL)THIS%EXPECTVAL(IB1),OCC(IB1,IKPT,ISPIN) &
    &                       ,VECTOR1(:,:,IB1)
                ELSE
                  WRITE(NFIL)THIS%EXPECTVAL(IB1),VECTOR1(:,:,IB1)
                END IF
              END IF
            ENDDO
          END IF
          DEALLOCATE(VECTOR1)
        ENDDO
      ENDDO
      DEALLOCATE(XK)
      DEALLOCATE(VEC)
      DEALLOCATE(R)
      DEALLOCATE(OCC)
      IF(THISTASK.EQ.1) THEN
        CALL LIB$FLUSHFILE(NFIL)
        CALL FILEHANDLER$CLOSE('PDOS')
      END IF
                             CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$REPORTEIG(NFIL)
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: IKPT
      INTEGER(4)            :: ISPIN
      INTEGER(4)            :: IB
      INTEGER(4)            :: ITEN
      REAL(8)               :: EV
      CHARACTER(64)         :: STRING
      INTEGER(4)            :: NB
!     ******************************************************************
      CALL CONSTANTS('EV',EV)
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(.NOT.ASSOCIATED(THIS%EIGVAL)) THEN
            CALL ERROR$MSG('EIGENVALUES NOT PRESENT')
            CALL ERROR$STOP('WAVES$REPORTEIG')
          END IF
          NB=THIS%NB
          IF(NSPIN.EQ.1) THEN
            WRITE(STRING,FMT='("EIGENVALUES [EV] FOR K-POINT ",I4)')IKPT
          ELSE
            WRITE(STRING,FMT='("EIGENVALUES [EV] FOR K-POINT ",I4' &
     &                        //'," AND SPIN ",I1)')IKPT,ISPIN
          END IF
          CALL REPORT$TITLE(NFIL,STRING)
          ITEN=0
          DO WHILE (NB.GT.ITEN)
            WRITE(NFIL,FMT='(I3,":",10F8.3)') &
     &           ITEN,(THIS%EIGVAL(IB)/EV,IB=ITEN+1,MIN(ITEN+10,NB))
            ITEN=ITEN+10
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$WRITE(NFIL,NFILO,TCHK)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: NFIL
      INTEGER(4)            ,INTENT(IN) :: NFILO
      LOGICAL(4)            ,INTENT(OUT):: TCHK
      TYPE (SEPARATOR_TYPE),PARAMETER  :: SEP_WAVES &
           =SEPARATOR_TYPE(0,'WAVES','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE)            :: SEPARATOR
      INTEGER(4)                       :: THISTASK,NTASKS
      INTEGER(4)                       :: IKPT,ISPIN
      INTEGER(4)                       :: NB,NBH
!     ******************************************************************
              CALL TRACE$PUSH('WAVES$WRITE')
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     == WRITE WAVE FUNCTION FOR THE THIS TIME STEP ==================
      TCHK=.TRUE.
      SEPARATOR=SEP_WAVES
      SEPARATOR%NREC=-1
      CALL WAVES_WRITEPSI(NFIL,SEPARATOR%NREC)
      CALL WRITESEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      CALL WAVES_WRITEPSI(NFIL,SEPARATOR%NREC)
              CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$READ(NFIL,NFILO,TCHK)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **  --TCHK RETURNS IF ANYTHING HAS BEEN READ. AS RETURN CODE    **
!     **    IT SHOULD BETTER BE THAT ALL RELEVANT QUANTITIES HAVE     **
!     **    BEEN READ                                                 **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)           ,INTENT(IN) :: NFIL
      INTEGER(4)           ,INTENT(IN) :: NFILO
      LOGICAL(4)           ,INTENT(OUT):: TCHK  ! SOMETHING HAS BEEN READ
      TYPE (SEPARATOR_TYPE),PARAMETER  :: SEP_WAVES &
           =SEPARATOR_TYPE(0,'WAVES','NONE','AUG1996','NONE')
      TYPE (SEPARATOR_TYPE)            :: SEPARATOR
      LOGICAL(4)                       :: TREAD
      REAL(8)              ,ALLOCATABLE:: EIG(:,:,:)
      INTEGER(4)                       :: IKPT,ISPIN,IB
      INTEGER(4)                       :: NKPT1,NSPIN1,NB1,NBX
      INTEGER(4)                       :: NB
      COMPLEX(8)           ,ALLOCATABLE:: TMP(:,:)
      REAL(8)              ,ALLOCATABLE:: TMPR8(:,:)
      INTEGER(4)                       :: THISTASK,NTASKS
!     ******************************************************************
                                  CALL TRACE$PUSH('WAVES$READ')
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     ==  READ WAVE FUNCTIONS FOR THIS TIME STEP                      ==
!     ==================================================================
      TCHK=.TRUE.
      SEPARATOR=SEP_WAVES
      CALL READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      IF(.NOT.TCHK) RETURN
      CALL WAVES_READPSI(NFIL)
!
!     ==================================================================
!     == SET OCCUPATIONS FOR DYNOCC OBJECT                            ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(EIG(NBX,NKPT,NSPIN))
      EIG(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          DO IB=1,THIS%NB
            EIG(IB,IKPT,ISPIN)=REAL(THIS%RLAM0(IB,IB))
          ENDDO
        ENDDO
      ENDDO
      CALL DYNOCC$SETR8A('EPSILON',NBX*NKPT*NSPIN,EIG)
      DEALLOCATE(EIG)
!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
                                  CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_WRITEPSI(NFIL,NREC)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL 
      INTEGER(4)  ,INTENT(INOUT):: NREC
      COMPLEX(8)  ,ALLOCATABLE:: PSIG(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSI1(:,:)
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: IOS
      INTEGER(4)              :: IKPT,ISPIN,IB,IDIM,IWAVE
      INTEGER(4)              :: LEN
      INTEGER(4)              :: NG1,NGG,NGL,NBH,NB
      LOGICAL(4)              :: TSUPER
      CHARACTER(64)           :: IOSTATMSG
      REAL(8)                 :: XK(3)
      REAL(8)     ,ALLOCATABLE:: GVECL(:,:)
      REAL(8)     ,ALLOCATABLE:: GVECG(:,:)
      REAL(8)                 :: GBAS(3,3)
      INTEGER(4)              :: NREC1,ISVAR
      CHARACTER(8)            :: KEY
      REAL(8)                 :: RBAS(3,3)
      INTEGER(4)  ,ALLOCATABLE:: IGVECG(:,:)
      INTEGER(4)  ,ALLOCATABLE:: IGVECL(:,:)
      INTEGER(4)              :: NWAVE=2
!     ******************************************************************
              CALL TRACE$PUSH('WAVES_WRITEPSI')
      CALL MPE$QUERY(NTASKS,THISTASK)
      IF(RSTRTTYPE.EQ.'STATIC') THEN
        NWAVE=1
      ELSE IF(RSTRTTYPE.EQ.'DYNAMIC') THEN
        NWAVE=2
      ELSE
        CALL ERROR$STOP('WAVES_WRITEPSI')
      END IF
!
!     ==================================================================
!     ==  COUNT #(RECORDS)                                            ==
!     ==================================================================
      IF(NREC.EQ.-1) THEN
        NREC=1
        DO IKPT=1,NKPT
          NREC=NREC+2
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          TSUPER=GSET%TINV
          NREC=NREC+NWAVE*THIS%NBH*NSPIN
          ISVAR=-1
          CALL WAVES_WRITELAMBDA(NFIL,IKPT,ISVAR)
          NREC=NREC+ISVAR
        ENDDO  
        CALL TRACE$POP
        RETURN
      END IF
      NREC1=0
!
!     ==================================================================
!     ==  GET DIMENSIONS ETC.                                         ==
!     ==================================================================
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     ==  WRITE SIZES                                                 ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
!       WRITE(NFIL)NKPT,NSPIN !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        CALL CELL$GETR8A('T(0)',9,RBAS)
        WRITE(NFIL)NKPT,NSPIN,RBAS,NWAVE !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        NREC1=NREC1+1
      END IF
!
!     ==================================================================
!     ==  LOOP OVER K-POINTS AND SPINS                                ==
!     ==================================================================
              CALL TRACE$PASS('MARKE 5')
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NGL=GSET%NGL
        NBH=THIS%NBH
        NB=THIS%NB
        CALL PLANEWAVE$GETI4('NGG',NGG)
!       
!       ================================================================
!       == WRITE SIZE AND TYPE OF WAVE FUNCTION                       ==
!       ================================================================
        CALL PLANEWAVE$GETL4('TINV',TSUPER)
        IF(THISTASK.EQ.1) THEN
          KEY='PSI'
          WRITE(NFIL)KEY,NGG,NDIM,NB,TSUPER  !<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC1=NREC1+1
        END IF
!       
!       ================================================================
!       == WRITE K-POINTS AND G-VECTORS                               ==
!       ================================================================
        CALL PLANEWAVE$GETR8A('GBAS',9,GBAS)
        CALL PLANEWAVE$GETR8A('XK',3,XK)
        ALLOCATE(IGVECL(3,NGL))
        CALL PLANEWAVE$GETI4A('IGVEC',3*NGL,IGVECL)
        ALLOCATE(IGVECG(3,NGG))
        CALL PLANEWAVE$COLLECTI4(3,NGL,IGVECL,NGG,IGVECG)
        IF(THISTASK.EQ.1) THEN
          WRITE(NFIL)XK,IGVECG !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC1=NREC1+1
        END IF
        DEALLOCATE(IGVECG)
        DEALLOCATE(IGVECL)
        ALLOCATE(PSIG(NGG,NDIM))
        ALLOCATE(PSI1(NGL,NDIM))
        DO IWAVE=1,NWAVE
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            CALL PLANEWAVE$SELECT(GSET%ID)
            DO IB=1,NBH
              DO IDIM=1,NDIM
                IF(IWAVE.EQ.1) THEN
                  PSI1(:,IDIM)=THIS%PSI0(:,IDIM,IB)
                ELSE IF(IWAVE.EQ.2) THEN
                  PSI1(:,IDIM)=THIS%PSIM(:,IDIM,IB)
                END IF
              ENDDO
              DO IDIM=1,NDIM
                CALL PLANEWAVE$COLLECTC8(1,NGL,PSI1(1,IDIM),NGG,PSIG(1,IDIM))
              ENDDO
              IF(THISTASK.EQ.1) THEN
                WRITE(NFIL,ERR=9999,IOSTAT=IOS)PSIG !<<<<<<<<<<<<<<<<<<<
                NREC1=NREC1+1
              END IF
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(PSI1)
        DEALLOCATE(PSIG)
        ISVAR=0
        CALL WAVES_WRITELAMBDA(NFIL,IKPT,ISVAR)
        NREC1=NREC1+ISVAR
      ENDDO 
 !
!     ==================================================================
!     ==  DEALLOCATE                                                  ==
!     ==================================================================
      IF(THISTASK.EQ.1.AND.NREC1.NE.NREC) THEN
        CALL ERROR$MSG('#(RECORDS WRITTEN DIFFERENT FROM TARGET')
        CALL ERROR$I4VAL('TARGET',NREC)
        CALL ERROR$I4VAL('ACTUAL',NREC1)
        CALL ERROR$STOP('WAVES_WRITEPSI')
      END IF
              CALL TRACE$POP
      RETURN
 9999 CONTINUE
      CALL FILEHANDLER$IOSTATMESSAGE(IOS,IOSTATMSG)
      CALL ERROR$MSG('ERROR WRITING WAVE FUNCTION TO RESTART FILE')
      CALL ERROR$I4VAL('IOS',IOS)
      CALL ERROR$CHVAL('IOSTATMSG',IOSTATMSG)
      CALL ERROR$I4VAL('IB',IB)
      CALL ERROR$I4VAL('IKPT',IKPT)
      CALL ERROR$I4VAL('ISPIN',ISPIN)
      CALL ERROR$I4VAL('NGG',NGG)
      CALL ERROR$STOP('WAVES_WRITEPSI')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_READPSI(NFIL)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     ******************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL 
      INTEGER(4)              :: NTASKS,THISTASK
      INTEGER(4)              :: NKPT_,NSPIN_
      INTEGER(4)              :: NGG_,NDIM_,NB_,NBH_
      LOGICAL(4)              :: TSUPER_,TSUPER
      INTEGER(4)              :: IKPT_,ISPIN_
      INTEGER(4)              :: IKPT1,IKPT,IB,IDIM,IG,IB1,IB2,ISPIN,I
      INTEGER(4)              :: IKPT0
      INTEGER(4)              :: IWAVE
      LOGICAL(4)              :: TREAD(NKPT)
      REAL(8)                 :: KREAD(3,NKPT)
      REAL(8)                 :: K(3)        ! ACTUAL K-POINT 
      REAL(8)                 :: K_(3)       ! K-POINT ON FILE
      REAL(8)                 :: GBAS_(3,3)  ! REC. LATT. VECT. ON FILE
      REAL(8)     ,ALLOCATABLE:: GVECG_(:,:) ! G-VECTORS ON FILE
      REAL(8)     ,ALLOCATABLE:: GVECG(:,:)  ! G-VECTORS (GLOBAL)
      REAL(8)     ,ALLOCATABLE:: GVECL(:,:)  ! G-VECTORS (LOCAL)
      REAL(8)                 :: RBAS(3,3)
      REAL(8)                 :: SVAR,SVAR1
      INTEGER(4)              :: NGG,NGL,NBH,NB
      INTEGER(4)  ,ALLOCATABLE:: MAPG(:)
      COMPLEX(8)  ,ALLOCATABLE:: PSI(:,:,:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSIL(:,:,:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSIIN(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: PSIG(:,:)
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)              :: ISVAR1
      INTEGER(4)              :: IOS
      CHARACTER(8)            :: KEY
      LOGICAL(4)              :: TCYCLE
      LOGICAL(4)              :: GBASFIX
      INTEGER(4)              :: IFORMAT
      INTEGER(4) ,ALLOCATABLE :: IGVECG_(:,:)
      REAL(8)                 :: XG1,XG2,XG3
      INTEGER(4)              :: NWAVE
      INTEGER(4) ,ALLOCATABLE :: MINUSG(:)
      LOGICAL(4)              :: TCHK
      COMPLEX(8)              :: F1,F2
      COMPLEX(8),ALLOCATABLE  :: PSIINSUPER(:,:)
      COMPLEX(8)              :: csvar,cmat(1,1)
!     ******************************************************************
                               CALL TRACE$PUSH('WAVES_READPSI')
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     ==  READ SIZES                                                  ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
!       IFORMAT=1
!        READ(NFIL,ERR=100)NKPT_,NSPIN_,RBAS
!        IFORMAT=2
! 100    CONTINUE
        NWAVE=2
        IFORMAT=0
        READ(NFIL,ERR=100)NKPT_,NSPIN_
        IFORMAT=1
        BACKSPACE(NFIL)
        READ(NFIL,ERR=100)NKPT_,NSPIN_,RBAS
        IFORMAT=2
        BACKSPACE(NFIL)
        READ(NFIL,ERR=100)NKPT_,NSPIN_,RBAS,NWAVE
 100    CONTINUE
        IF(IFORMAT.EQ.0) THEN
          CALL ERROR$MSG('FORMAT NOT RECOGNIZED')
          CALL ERROR$STOP('WAVES_READPSI')
        END IF
      END IF
      CALL MPE$BROADCAST(1,NSPIN_)
      CALL MPE$BROADCAST(1,NKPT_)
      CALL MPE$BROADCAST(1,IFORMAT)
      CALL MPE$BROADCAST(1,NWAVE)
      IF(IFORMAT.EQ.2) THEN
        CALL MPE$BROADCAST(1,RBAS)
        DO IKPT=1,NKPT
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
        ENDDO
        CALL GBASS(RBAS,GBAS_,SVAR)
      END IF
!
!     ==================================================================
!     ==  LOOP OVER K-POINTS AND SPINS                                ==
!     ==================================================================
              CALL TRACE$PASS('MARKE 5')
      GBASFIX=.FALSE.
      TREAD(:)=.FALSE.
      DO IKPT_=1,NKPT_
!
!       ================================================================
!       ==  READ COORDINATES OF THE WAVE FUNCTIONS                    ==
!       ================================================================
        IF(THISTASK.EQ.1) THEN
          READ(NFIL)KEY,NGG_,NDIM_,NB_,TSUPER_   !<<<<<<<<<<<<<<<<<<<<<<
        END IF
        CALL MPE$BROADCAST(1,KEY)
        CALL MPE$BROADCAST(1,NGG_)
        CALL MPE$BROADCAST(1,NDIM_)
        CALL MPE$BROADCAST(1,NB_)
        CALL MPE$BROADCAST(1,TSUPER_)
        IF(KEY.NE.'PSI') THEN
          CALL ERROR$MSG('ID IS NOT "PSI"')
          CALL ERROR$MSG('FILE IS CORRUPTED')
          CALL ERROR$I4VAL('IKPT_',IKPT_)
          CALL ERROR$STOP('WAVES_READPSI')
        END IF
        ALLOCATE(GVECG_(3,NGG_))
        IF(THISTASK.EQ.1) THEN
          IF(IFORMAT.EQ.1) THEN
            READ(NFIL)K_,GBAS_,GVECG_ !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ELSE
            ALLOCATE(IGVECG_(3,NGG_))
            READ(NFIL)K_,IGVECG_ !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            DO IG=1,NGG_
              XG1=REAL(IGVECG_(1,IG),8)+K_(1)
              XG2=REAL(IGVECG_(2,IG),8)+K_(2)
              XG3=REAL(IGVECG_(3,IG),8)+K_(3)
              DO I=1,3
                GVECG_(I,IG)=GBAS_(I,1)*XG1+GBAS_(I,2)*XG2+GBAS_(I,3)*XG3
              ENDDO
            ENDDO
            IF(TSUPER_) THEN
              ALLOCATE(MINUSG(NGG_))
              CALL PLANEWAVE$MINUSG(K_,NGG_,IGVECG_,MINUSG)
            END IF
            DEALLOCATE(IGVECG_)
          END IF
        END IF
        CALL MPE$BROADCAST(1,K_)
        CALL MPE$BROADCAST(1,GVECG_)
!
        IF(IFORMAT.EQ.1.AND.(.NOT.GBASFIX)) THEN
          GBASFIX=.TRUE.
          CALL MPE$BROADCAST(1,GBAS_)
          CALL GBASS(GBAS_,RBAS,SVAR)
          DO IKPT1=1,NKPT
            CALL WAVES_SELECTWV(IKPT1,1)
            CALL PLANEWAVE$SELECT(GSET%ID)
            CALL PLANEWAVE$SETR8A('RBAS',9,RBAS)
          ENDDO
        END IF
!       
!       ==================================================================
!       ==  FIND NEAREST K-POINT                                        ==
!       ==================================================================
        SVAR=1.D+10
        IKPT=1
        DO IKPT1=1,NKPT
          CALL WAVES_SELECTWV(IKPT1,1)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETR8A('XK',3,K)
          SVAR1=(K(1)-K_(1))**2+(K(2)-K_(2))**2+(K(3)-K_(3))**2
          IF(SVAR1.LT.SVAR) THEN
            IKPT=IKPT1
            SVAR=SVAR1
          END IF
        ENDDO
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
!       
!       ==================================================================
!       ==  FIND NEAREST K-POINT                                        ==
!       ==================================================================
        IF(TREAD(IKPT)) THEN
          CALL PLANEWAVE$GETR8A('XK',3,K)
          SVAR1=(K(1)-K_(1))**2-(K(1)-KREAD(1,IKPT))**2 &
     &         +(K(2)-K_(2))**2-(K(2)-KREAD(2,IKPT))**2 &
     &         +(K(3)-K_(3))**2-(K(3)-KREAD(3,IKPT))**2 
          IF(SVAR1.GT.0.D0) THEN  ! PREVIOUS CHOICE WAS BETTER
            NBH_=NB_
            IF(TSUPER_)NBH_=INT(0.5D0*REAL(NB+1,8))
            DO IWAVE=1,NWAVE
              DO ISPIN_=1,NSPIN_
                DO IB=1,NBH_
                  IF(THISTASK.EQ.1)READ(NFIL)
                ENDDO
              ENDDO
            ENDDO
            CALL WAVES_READLAMBDA(NFIL,IKPT,.TRUE.)
            DEALLOCATE(GVECG_)
            IF(ALLOCATED(MINUSG))DEALLOCATE(MINUSG)
            CYCLE
          END IF
        END IF
        TREAD(IKPT)=.TRUE.
        KREAD(:,IKPT)=K_(:)
        CALL TRACE$PASS('MARKE 6')
!       == SET DATA FOR FURTHER USE ==================================
        NGL=GSET%NGL
        CALL PLANEWAVE$GETI4('NGG',NGG)
!       
!       ==============================================================
!       ==  DEFINE MAPPING OF THE ARRAYS                            ==
!       ==============================================================
        ALLOCATE(GVECL(3,NGL))
        ALLOCATE(GVECG(3,NGG))
        CALL PLANEWAVE$GETR8A('GVEC',3*NGL,GVECL)
        CALL PLANEWAVE$COLLECTR8(3,NGL,GVECL,NGG,GVECG)
        DEALLOCATE(GVECL)
        ALLOCATE(MAPG(NGG))
        IF(THISTASK.EQ.1) THEN
          CALL WAVES_MAPG(NGG_,GBAS_,GVECG_,NGG,GVECG,MAPG)
        END IF
        CALL MPE$BROADCAST(1,MAPG)
!!$DO IG=1,NGG
!!$  IF(MAPG(IG).NE.IG) THEN
!!$     PRINT*,'MAPG ',IKPT_,IG,MAPG(IG),NGG,NGG_
!!$  END IF
!!$ENDDO
        DEALLOCATE(GVECG_)
        DEALLOCATE(GVECG)
              CALL TRACE$PASS('MARKE 7')
!       
!       ==============================================================
!       ==  COLLECT DATA                                            ==
!       ==============================================================
IF(TSUPER_.AND.NDIM_.EQ.2) THEN
  CALL ERROR$MSG('NONCOLLINEAR WAVE FUNCTIONS AND SUPER WAVE FUNCTIONS')
  CALL ERROR$MSG('ARE NOT PROPERLY IMPLEMENTED BELOW')
  CALL ERROR$STOP('WAVES_READPSI')
END IF
        NBH=THIS%NBH
        NB=THIS%NB
        ALLOCATE(PSIIN(NGG_,NDIM_))
        ALLOCATE(PSIG(NGG,NDIM_))
        ALLOCATE(PSIL(NGL,NDIM_,NB_,NSPIN_))
        ALLOCATE(PSI(NGL,NDIM,NB,NSPIN))
        IF(TSUPER_) ALLOCATE(PSIINSUPER(NGG_,NDIM_))
        DO IWAVE=1,NWAVE
          DO ISPIN=1,NSPIN_
            TCHK=.FALSE.
            DO IB=1,NB_
              IF(THISTASK.EQ.1) THEN 
                IF(.NOT.TSUPER_) THEN
                  READ(NFIL,ERR=9999,IOSTAT=IOS)PSIIN !<<<<<<<<<<<<<<<<<<<
                ELSE  ! THIS FOR READING SUPER WAVE FUNCTIONS
                  IF(.NOT.TCHK) THEN   
                    READ(NFIL,ERR=9999,IOSTAT=IOS)PSIINSUPER !<<<<<<<<<<<<
                    DO IG=1,NGG_
                      IF(MINUSG(IG).LT.IG) CYCLE
                      DO IDIM=1,NDIM_
                        F1=PSIINSUPER(IG,IDIM)
                        F2=CONJG(PSIINSUPER(MINUSG(IG),IDIM))
                        PSIIN(IG,IDIM)=0.5D0*(F1+F2)
                        PSIIN(MINUSG(IG),IDIM)=CONJG(PSIIN(IG,IDIM))
                      ENDDO
                    ENDDO
                    TCHK=.TRUE.
                  ELSE
                    DO IG=1,NGG_
                      IF(MINUSG(IG).LT.IG) CYCLE
                      DO IDIM=1,NDIM_
                        F1=PSIINSUPER(IG,IDIM)
                        F2=CONJG(PSIINSUPER(MINUSG(IG),IDIM))
                        PSIIN(IG,IDIM)=0.5D0*CI*(F1-F2)
                        PSIIN(MINUSG(IG),IDIM)=CONJG(PSIIN(IG,IDIM))
                      END DO
                    ENDDO
                    TCHK=.FALSE.
                  END IF
                END IF              
                DO IDIM=1,NDIM_
                  DO IG=1,NGG
                    IF(MAPG(IG).NE.0) THEN
                      PSIG(IG,IDIM)=PSIIN(MAPG(IG),IDIM)
                    ELSE
                      PSIG(IG,IDIM)=(0.D0,0.D0)
                    END IF
                  ENDDO
                ENDDO
              END IF
              DO IDIM=1,NDIM_
                CALL PLANEWAVE$DISTRIBUTEC8(1,NGG,PSIG(1,IDIM),NGL,PSIL(1,IDIM,IB,ISPIN))
              ENDDO
            ENDDO
          ENDDO
!       
!         ==============================================================
!         ==  TRANSFORM INTO THE CORRECT FORMAT                       ==
!         ==============================================================
!         == SAME REPRESENTATION =========================================
          CALL TRACE$PASS('MARKE 9')
          IF(NDIM.EQ.NDIM_.AND.NSPIN.EQ.NSPIN_) THEN
            DO IB=1,NB
              IF(IB.GT.NB_) EXIT
              PSI(:,:,IB,:)=PSIL(:,:,IB,:)
            ENDDO
!         == FROM NON-SPIN POLARIZED CALCULATION =========================
          ELSE IF(NSPIN_.EQ.1.AND.NDIM_.EQ.1) THEN
            IF(NSPIN.EQ.2) THEN
              DO ISPIN=1,NSPIN
                DO IB=1,NB
                  IF(IB.GT.NB_) EXIT
                  PSI(:,1,IB,ISPIN)=PSIL(:,1,IB,1)
                ENDDO
              ENDDO
            ELSE IF(NDIM.EQ.2) THEN
              DO IB=1,NB
                IB1=INT(0.5*REAL(IB+1)) ! IB1 =1,1,2,2,3,3,...
                IF(IB1.GT.NB_) EXIT
                IDIM=2+IB-2*IB1         ! IDIM=1,2,1,2,1,2,...
                PSI(:,IDIM,IB,1)=PSIL(:,1,IB1,1)
              ENDDO
            ELSE
              CALL ERROR$MSG('TRANSFORMATION NOT IMPLEMENTED')
              CALL ERROR$STOP('WAVES$READPSI')
            END IF
!         == FROM SPIN POLARIZED CALCULATION ===========================
          ELSE IF(NSPIN_.EQ.2.AND.NDIM_.EQ.1) THEN
            IF(NSPIN.EQ.1.AND.NDIM.EQ.1) THEN
              DO IB=1,NB
                PSI(:,1,IB,1)=PSIL(:,1,IB,1)
              ENDDO
            ELSE IF(NSPIN.EQ.1.AND.NDIM.EQ.2) THEN
              DO IB=1,(NB+1)/2
                IB1=2*IB-1
                IB2=2*IB
                PSI(:,1,IB1,1)=PSIL(:,1,IB,1)
                PSI(:,2,IB1,1)=(0.D0,0.D0)
                IF(IB2.GT.NB) EXIT
                PSI(:,1,IB2,1)=(0.D0,0.D0)
                PSI(:,2,IB2,1)=PSIL(:,1,IB,2)
              ENDDO
            ELSE
              CALL ERROR$MSG('TRANSFORMATION NOT IMPLEMENTED')
              CALL ERROR$STOP('WAVES$READPSI')
            END IF
!         == FROM NONCOLLINEAR CALCULATION  ==============================
          ELSE IF(NSPIN_.EQ.2.AND.NDIM_.EQ.1) THEN
            IF(NSPIN.EQ.1.AND.NDIM.EQ.1) THEN
              DO IB=1,NB 
                PSI(:,1,IB,1)=PSIL(:,1,2*IB-1,1)
              ENDDO
            ELSE IF(NSPIN.EQ.1.AND.NDIM.EQ.1) THEN
              DO IB=1,NB 
                IB1=2*IB-1
                IB2=2*IB
                PSI(:,1,IB,1)=PSIL(:,1,IB1,1)
                PSI(:,1,IB,2)=PSIL(:,1,IB2,1)
              ENDDO
            ELSE
              CALL ERROR$MSG('TRANSFORMATION NOT IMPLEMENTED')
              CALL ERROR$STOP('WAVES$READPSI')
            END IF
          END IF           
          CALL TRACE$PASS('MARKE 10')
!       
!         ==============================================================
!         ==  MAP ONTO SUPER WAVE FUNCTIONS                           ==
!         ==============================================================
          CALL PLANEWAVE$GETL4('TINV',TSUPER)
          IF(TSUPER) THEN
            DO ISPIN=1,NSPIN
!             == remove the factor exp(i*phi) from the wave function
              do ib=1,nb
                call PLANEWAVE$SCALARPRODUCT('-',NGL,1,1 &
        &                   ,psi(:,1,ib,ispin),1,psi(:,1,ib,ispin),cmat)
                csvar=cmat(1,1)
                csvar=conjg(sqrt(csvar/sqrt(csvar*conjg(csvar))))
                psi(:,1,ib,ispin)=psi(:,1,ib,ispin)*csvar
              enddo
              DO IB=1,NBH
                IB1=2*IB-1
                IB2=2*IB
                PSI(:,1,IB,ISPIN)=PSI(:,1,IB1,ISPIN)+CI*PSI(:,1,IB2,ISPIN)
              ENDDO
            ENDDO
          END IF
!       
!         ==============================================================
!         ==  MAP BACK                                                ==
!         ==============================================================
          CALL TRACE$PASS('MARKE 11')
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            IF(IWAVE.EQ.1) THEN
              DO IB=1,NBH
                DO IDIM=1,NDIM
                  DO IG=1,NGL
                    THIS%PSI0(IG,IDIM,IB)=PSI(IG,IDIM,IB,ISPIN)
                  ENDDO
                ENDDO
              ENDDO
              IF(NWAVE.EQ.1) THEN
                THIS%PSIM=THIS%PSI0
              END IF
            ELSE IF(IWAVE.EQ.2) THEN
              DO IB=1,NBH
                DO IDIM=1,NDIM
                  DO IG=1,NGL
                    THIS%PSIM(IG,IDIM,IB)=PSI(IG,IDIM,IB,ISPIN)
                  ENDDO
                ENDDO
              ENDDO
            END IF
          ENDDO
        ENDDO 
        DEALLOCATE(PSIL)
        DEALLOCATE(PSIG)
        DEALLOCATE(MAPG)
        DEALLOCATE(PSI)
        DEALLOCATE(PSIIN)
        IF(ALLOCATED(PSIINSUPER))DEALLOCATE(PSIINSUPER)
        IF(ALLOCATED(MINUSG)) DEALLOCATE(MINUSG)
        CALL WAVES_READLAMBDA(NFIL,IKPT,.FALSE.)
      ENDDO
!
!     ==================================================================
!     ==  COMPLETE K-POINTS                                           ==
!     ==================================================================
      CALL TRACE$PASS('MARKE 14')
      DO IKPT=1,NKPT
        IF(TREAD(IKPT)) CYCLE
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        CALL PLANEWAVE$GETR8A('K',3,K)
        SVAR=1.D+10
        DO IKPT_=1,NKPT_
          IF(.NOT.TREAD(IKPT_)) CYCLE
          SVAR1=(K(1)-KREAD(1,IKPT_))**2+(K(2)-KREAD(2,IKPT_))**2 &
     &                                  +(K(3)-KREAD(3,IKPT_))**2
          IF(SVAR1.LT.SVAR) THEN
            IKPT0=IKPT_
            SVAR=SVAR1
          END IF
        ENDDO       
        CALL WAVES_COPYPSI(IKPT,IKPT0)
        CALL WAVES_COPYLAMBDA(IKPT0,IKPT)
      ENDDO 

!
!     ==================================================================
!     ==  CLOSE DOWN                                                  ==
!     ==================================================================
                               CALL TRACE$POP
      RETURN
 9999 CONTINUE
      CALL ERROR$MSG('ERROR WHILE READING FROM FILE')
      CALL ERROR$I4VAL('IOS',IOS)
      CALL ERROR$STOP('WAVES_READPSI')
      END
!
!     ..................................................................
      SUBROUTINE WAVES_WRITELAMBDA(NFIL,IKPT,NREC)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NFIL
      INTEGER(4)   ,INTENT(IN) :: IKPT
      INTEGER(4)   ,INTENT(INOUT) :: NREC
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)               :: ISPIN,I
      INTEGER(4)               :: NB
      CHARACTER(8)             :: KEY
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA(:,:,:)
      INTEGER(4)               :: NLAMBDA
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
      CALL WAVES_SELECTWV(IKPT,1)
!
!     ==================================================================
!     == COUNT DEPTH OF HISTORY FOR LAGRANGE PARAMETERS               ==
!     ==================================================================
      NLAMBDA=4
      IF(.NOT.ASSOCIATED(THIS%RLAM3M)) NLAMBDA=3
      IF(.NOT.ASSOCIATED(THIS%RLAM2M)) NLAMBDA=2
      IF(.NOT.ASSOCIATED(THIS%RLAMM)) NLAMBDA=1
      IF(.NOT.ASSOCIATED(THIS%RLAM0)) NLAMBDA=0
      IF(RSTRTTYPE.EQ.'STATIC') NLAMBDA=0
!
!     ==================================================================
!     == RETURN #(RECORDS) OR IF ON TASK OTHER THAN THE FIRST         ==
!     ==================================================================
      IF(NREC.EQ.-1) THEN
        NREC=1+NSPIN*NLAMBDA
        RETURN
      END IF
      IF(THISTASK.NE.1) RETURN
      NREC=0
!
!     ==================================================================
!     == NOW WRITE TO FILE                                            ==
!     ==================================================================
      KEY='LAMBDA'
      NB=THIS%NB
      WRITE(NFIL)KEY,NB,NDIM,NSPIN,NLAMBDA  !<<<<<<<<<<<<<<<<<<<<<<<<<<<
      NREC=NREC+1
      IF(NLAMBDA.GE.1) THEN
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          WRITE(NFIL)THIS%RLAM0 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC=NREC+1
        ENDDO
      END IF
      IF(NLAMBDA.GE.2) THEN
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          WRITE(NFIL)THIS%RLAMM !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC=NREC+1
        ENDDO
      END IF
      IF(NLAMBDA.GE.3) THEN
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          WRITE(NFIL)THIS%RLAM2M !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC=NREC+1
        ENDDO
      END IF
      IF(NLAMBDA.GE.4) THEN
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          WRITE(NFIL)THIS%RLAM3M !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          NREC=NREC+1
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_READLAMBDA(NFIL,IKPT,TSKIP)
!     ******************************************************************
!     **                                                              **
!     **  WRITE WAVE FUNCTIONS TO RESTART FILE                        **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NFIL
      INTEGER(4)   ,INTENT(IN) :: IKPT
      LOGICAL(4)   ,INTENT(IN) :: TSKIP
      INTEGER(4)               :: NTASKS,THISTASK
      INTEGER(4)               :: ISPIN,I,IB1,IB2
      INTEGER(4)               :: NB_,NDIM_,NSPIN_
      INTEGER(4)               :: NB,NBA
      INTEGER(4)               :: NLAMBDA
      CHARACTER(8)             :: KEY
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA(:,:,:)
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA2(:,:,:)
!     ******************************************************************
                           CALL TRACE$PUSH('WAVES_READLAMBDA')
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     == READ AND BROADCAST HISTORY-DEPTH OF LAGRANGE PARAMETERS      ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        READ(NFIL)KEY,NB_,NDIM_,NSPIN_,NLAMBDA  !<<<<<<<<<<<<<<<<<<<<<<<<<<
        IF(KEY.NE.'LAMBDA') THEN
          CALL ERROR$MSG('ID IS NOT "LAMBDA"')
          CALL ERROR$MSG('FILE IS CORRUPTED')
          CALL ERROR$STOP('WAVES_READLAMBDA')
        END IF
      END IF
      IF(TSKIP) THEN
        DO I=1,NLAMBDA
          DO ISPIN=1,NSPIN_
           IF(THISTASK.EQ.1)READ(NFIL)
          ENDDO
        ENDDO
        CALL TRACE$POP
        RETURN
      END IF
      CALL MPE$BROADCAST(1,NLAMBDA)
!
!     ==================================================================
!     ==  READ AND FOLD DOWN                                          ==
!     ==================================================================
      IF(THISTASK.EQ.1) THEN
        CALL WAVES_SELECTWV(IKPT,1)
        NB=THIS%NB
        ALLOCATE(LAMBDA(NB_,NB_,NSPIN_))
        ALLOCATE(LAMBDA2(NB,NB,NSPIN))
        NBA=MIN(NB,NB_)
        DO I=1,NLAMBDA
          DO ISPIN=1,NSPIN_
            READ(NFIL)LAMBDA(:,:,ISPIN)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ENDDO
!
!         ==============================================================
!         ==  TRANSFORM BETWEEN DATA MODELS                           ==
!         ==============================================================
          LAMBDA2=(0.D0,0.D0)
          IF(NSPIN.EQ.NSPIN_.AND.NDIM.EQ.NDIM_) THEN
            LAMBDA2(1:NBA,1:NBA,:)=LAMBDA(1:NBA,1:NBA,:)
          ELSE
!           == FROM NON-SPIN POLARIZED CALCULATION =====================
            IF(NDIM_.EQ.1.AND.NSPIN_.EQ.1) THEN
              IF(NDIM.EQ.1.AND.NSPIN.EQ.2) THEN
                LAMBDA2(1:NBA,1:NBA,1)=LAMBDA(1:NBA,1:NBA,1)
                LAMBDA2(1:NBA,1:NBA,2)=LAMBDA(1:NBA,1:NBA,1)
              ELSE IF(NDIM.EQ.2.AND.NSPIN.EQ.1) THEN
                DO IB1=1,NB_
                  DO IB2=1,NB_
                    IF(2*IB1-1.GT.NB.OR.2*IB2-1.GT.NB)CYCLE
                    LAMBDA2(2*IB1-1,2*IB2-1,1)=LAMBDA(IB1,IB2,1)
                    IF(2*IB1.GT.NB.OR.2*IB2.GT.NB)CYCLE
                    LAMBDA2(2*IB1  ,2*IB2  ,1)=LAMBDA(IB1,IB2,1)
                  ENDDO
                ENDDO
              ELSE
                CALL ERROR$MSG('INVALID OPTION TRANSFORMING FROM')
                CALL ERROR$MSG('NON-SPINPOLARIZED DATA')
                CALL ERROR$I4VAL('NDIM',NDIM)
                CALL ERROR$I4VAL('NSPIN',NSPIN)
                CALL ERROR$I4VAL('NDIM_',NDIM_)
                CALL ERROR$I4VAL('NSPIN_',NSPIN_)
                CALL ERROR$STOP('WAVES_READLAMBDA')
              END IF
!           == FROM SPIN POLARIZED CALCULATION =========================
            ELSE IF(NDIM_.EQ.1.AND.NSPIN_.EQ.2) THEN
              IF(NDIM.EQ.1.AND.NSPIN.EQ.1) THEN
                LAMBDA2(1:NBA,1:NBA,1)=LAMBDA(1:NBA,1:NBA,1)
              ELSE IF(NDIM.EQ.2.AND.NSPIN.EQ.1) THEN
                DO IB1=1,NB_
                  DO IB2=1,NB_
                    IF(2*IB1-1.GT.NB.OR.2*IB2-1.GT.NB)CYCLE
                    LAMBDA2(2*IB1-1,2*IB2-1,1)=LAMBDA(IB1,IB2,1)
                    IF(2*IB1.GT.NB.OR.2*IB2.GT.NB)CYCLE
                    LAMBDA2(2*IB1  ,2*IB2  ,1)=LAMBDA(IB1,IB2,2)
                  ENDDO
                ENDDO
              ELSE
                CALL ERROR$MSG('INVALID OPTION TRANSFORMING FROM')
                CALL ERROR$MSG('SPINPOLARIZED DATA')
                CALL ERROR$I4VAL('NDIM',NDIM)
                CALL ERROR$I4VAL('NSPIN',NSPIN)
                CALL ERROR$I4VAL('NDIM_',NDIM_)
                CALL ERROR$I4VAL('NSPIN_',NSPIN_)
                CALL ERROR$STOP('WAVES_READLAMBDA')
              END IF
!           == FROM NONCOLLINEAR WAVE FUNCTION =========================
            ELSE IF(NDIM_.EQ.2.AND.NSPIN_.EQ.1) THEN
              IF(NDIM.EQ.1.AND.NSPIN.EQ.1) THEN
                DO IB1=1,NB
                  DO IB2=1,NB
                    IF(2*IB1-1.GT.NB_.OR.2*IB2-1.GT.NB_)CYCLE
                    LAMBDA2(IB1,IB2,1)=LAMBDA(2*IB1-1,2*IB2-1,1)
                  ENDDO
                ENDDO
              ELSE IF(NDIM.EQ.2.AND.NSPIN.EQ.1) THEN
                DO IB1=1,NB
                  DO IB2=1,NB
                    IF(2*IB1-1.GT.NB_.OR.2*IB2-1.GT.NB_)CYCLE
                    LAMBDA2(IB1,IB2,1)=LAMBDA(2*IB1-1,2*IB2-1,1)
                    IF(2*IB1.GT.NB_.OR.2*IB2.GT.NB_)CYCLE
                    LAMBDA2(IB1,IB2,2)=LAMBDA(2*IB1,2*IB2,1)
                  ENDDO
                ENDDO
              ELSE
                CALL ERROR$MSG('INVALID OPTION TRANSFORMING FROM')
                CALL ERROR$MSG('NON-COLLINEAR DATA')
                CALL ERROR$I4VAL('NDIM',NDIM)
                CALL ERROR$I4VAL('NSPIN',NSPIN)
                CALL ERROR$I4VAL('NDIM_',NDIM_)
                CALL ERROR$I4VAL('NSPIN_',NSPIN_)
                CALL ERROR$STOP('WAVES_READLAMBDA')
              END IF
            ELSE
              CALL ERROR$MSG('INVALID OPTION')
              CALL ERROR$I4VAL('NDIM',NDIM)
              CALL ERROR$I4VAL('NSPIN',NSPIN)
              CALL ERROR$I4VAL('NDIM_',NDIM_)
              CALL ERROR$I4VAL('NSPIN_',NSPIN_)
              CALL ERROR$STOP('WAVES_READLAMBDA')
            END IF
          END IF
!
!         ==============================================================
!         ==  MAP ONTO THIS                                           ==
!         ==============================================================
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            NB=THIS%NB
            IF(I.EQ.1) THEN
              IF(.NOT.ASSOCIATED(THIS%RLAM0))ALLOCATE(THIS%RLAM0(NB,NB))
              THIS%RLAM0=LAMBDA2(:,:,ISPIN)
            ELSE IF(I.EQ.2) THEN
              IF(.NOT.ASSOCIATED(THIS%RLAMM))ALLOCATE(THIS%RLAMM(NB,NB))
              THIS%RLAMM=LAMBDA2(:,:,ISPIN)
            ELSE IF(I.EQ.3) THEN 
              IF(.NOT.ASSOCIATED(THIS%RLAM2M))ALLOCATE(THIS%RLAM2M(NB,NB))
              THIS%RLAM2M=LAMBDA2(:,:,ISPIN)
            ELSE IF(I.EQ.4) THEN
              IF(.NOT.ASSOCIATED(THIS%RLAM3M))ALLOCATE(THIS%RLAM3M(NB,NB))
              THIS%RLAM3M=LAMBDA2(:,:,ISPIN)
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(LAMBDA)
        DEALLOCATE(LAMBDA2)
      END IF
!
!     ==================================================================
!     ==  BROADCAST LAMBDA                                            ==
!     ==================================================================
      DO ISPIN=1,NSPIN
        CALL WAVES_SELECTWV(IKPT,ISPIN)
        NB=THIS%NB
        IF(NLAMBDA.GE.1)THEN
          IF(.NOT.ASSOCIATED(THIS%RLAM0))ALLOCATE(THIS%RLAM0(NB,NB))
          CALL MPE$BROADCAST(1,THIS%RLAM0)
        END IF
        IF(NLAMBDA.GE.2)THEN
          IF(.NOT.ASSOCIATED(THIS%RLAMM))ALLOCATE(THIS%RLAMM(NB,NB))
          CALL MPE$BROADCAST(1,THIS%RLAMM)
        END IF
        IF(NLAMBDA.GE.3)THEN
          IF(.NOT.ASSOCIATED(THIS%RLAM2M))ALLOCATE(THIS%RLAM2M(NB,NB))
          CALL MPE$BROADCAST(1,THIS%RLAM2M)
        END IF
        IF(NLAMBDA.GE.4)THEN
          IF(.NOT.ASSOCIATED(THIS%RLAM3M))ALLOCATE(THIS%RLAM3M(NB,NB))
          CALL MPE$BROADCAST(1,THIS%RLAM3M)
        END IF
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_MAPG(NGA,GBASA,GVECA,NGB,GVECB,MAP)
!     ******************************************************************
!     **  ARRAY TO MAP WAVE FUNCTIONS READ FROM FILE INTO             **
!     **  AN EXISTING GRID                                            **
!     **  NGA,GBASA,GVECA DESCRIBE THE ARRAY READ FROM FILE           **
!     **  NGB,GVECB DESCRIBE THE EXISTING G-ARRAY                     **
!     **                                                              **
!     **  DO IG=1,NG                                                  **
!     **    IF(MAP(IG).EQ.0) THEN                                     **
!     **      PSI(IG)=(0.D0,0.D0)                                     **
!     **    ELSE                                                      **
!     **      PSI(IG)=PSIIN(MAP(IG))                                  **
!     **    END IF                                                    **
!     **  ENDDO                                                       **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NGA
      REAL(8)   ,INTENT(IN) :: GBASA(3,3)
      REAL(8)   ,INTENT(IN) :: GVECA(3,NGA)
      INTEGER(4),INTENT(IN) :: NGB
      REAL(8)   ,INTENT(IN) :: GVECB(3,NGB)
      INTEGER(4),INTENT(OUT):: MAP(NGB)
      INTEGER(4)            :: IG,I
      REAL(8)   ,ALLOCATABLE:: MAP3D(:,:,:)
      REAL(8)               :: GBASIN(3,3)
      INTEGER(4)            :: IVEC1(3),IVECA(3,NGA)
      REAL(8)               :: KA(3)
      REAL(8)               :: G(3)
      INTEGER(4)            :: MING(3),MAXG(3)
!     ******************************************************************
      CALL LIB$INVERTR8(3,GBASA,GBASIN)
!
!     ==================================================================
!     == DETERMINE K-POINT                                            ==
!     ==================================================================
      G(:)=MATMUL(GBASIN,GVECA(:,1))
      IVEC1=NINT(G)
      DO I=1,3
        IVEC1(I)=NINT(G(I))
      ENDDO
      G=REAL(NINT(G))
      G(:)=MATMUL(GBASA,REAL(IVEC1,KIND=8))
      KA(:)=GVECA(:,1)-G(:)
!
!     ==================================================================
!     == DIMENSION AND ALLOCATE 3-D GRID                              ==
!     ==================================================================
      MING=0
      MAXG=0
      DO IG=1,NGA
        IVECA(:,IG)=NINT(MATMUL(GBASIN,GVECA(:,IG)-KA(:)))
        DO I=1,3
          MING(I)=MIN(MING(I),IVECA(I,IG))
          MAXG(I)=MAX(MAXG(I),IVECA(I,IG))
        ENDDO
      ENDDO
      ALLOCATE(MAP3D(MING(1):MAXG(1),MING(2):MAXG(2),MING(3):MAXG(3)))
!
!     ==================================================================
!     == MAP FIRST ARRAY ONTO 3-D GRID                                ==
!     ==================================================================
      MAP3D(:,:,:)=0
      DO IG=1,NGA
        IF(MAP3D(IVECA(1,IG),IVECA(2,IG),IVECA(3,IG)).EQ.0) THEN
          MAP3D(IVECA(1,IG),IVECA(2,IG),IVECA(3,IG))=IG
        ELSE
          CALL ERROR$MSG('TWO G-VECTORS ARE MAPPED ONTO THE SAME POINT')
          CALL ERROR$STOP('WAVES_MAPG')
        END IF
      ENDDO
!
!     ==================================================================
!     ==  FOR EACH VECTOR (B) PICK OUT THE CLOSEST GRID (A) POINT     ==
!     ==================================================================
      LOOP1:DO IG=1,NGB
        MAP(IG)=0
        G(:)=MATMUL(GBASIN,GVECB(:,IG)-KA(:))
        DO I=1,3
          IVEC1(I)=NINT(G(I))
          IF(IVEC1(I).LT.MING(I).OR.IVEC1(I).GT.MAXG(I)) CYCLE LOOP1
        ENDDO
        MAP(IG)=MAP3D(IVEC1(1),IVEC1(2),IVEC1(3))
      ENDDO LOOP1
      DEALLOCATE(MAP3D)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_COPYPSI(IKPT1,IKPT2)
!     ******************************************************************
!     **                                                              **
!     **  MAP WAVE FUNCTIONS FROM ONE K-POINT IKPT2 ONTO IKPT1        **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IKPT1
      INTEGER(4),INTENT(IN) :: IKPT2
      INTEGER(4)            :: NGL1,NGL2,NGG1,NGG2
      REAL(8)   ,ALLOCATABLE:: GVECL(:,:)
      REAL(8)   ,ALLOCATABLE:: GVECG1(:,:)
      REAL(8)   ,ALLOCATABLE:: GVECG2(:,:)
      INTEGER(4),ALLOCATABLE:: MAPG(:)
      LOGICAL(4)            :: TSUPER1,TSUPER2
      INTEGER(4)            :: NB1,NB2,NBH1,NBH2,NBN,NBHN
      COMPLEX(8),ALLOCATABLE:: PSI1(:,:,:)
      COMPLEX(8),ALLOCATABLE:: PSI2(:,:,:)
      COMPLEX(8),ALLOCATABLE:: PSIG1(:)
      COMPLEX(8),ALLOCATABLE:: PSIG2(:)
      COMPLEX(8),ALLOCATABLE:: PSITMP2(:,:)
      INTEGER(4)            :: IWAVE,ISPIN,IB,IB1,IB2,IDIM,IG
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)               :: GBAS(3,3)
      INTEGER(4)            :: NTASKS,THISTASK
!     ******************************************************************
      CALL MPE$QUERY(NTASKS,THISTASK)
!
!     ==================================================================
!     ==  MAPPING FROM IKPT2 TO IKPT1                                 ==
!     ==================================================================
!     == GET GLOBAL GVECTORS FROM IKPT1 ================================
      CALL WAVES_SELECTWV(IKPT1,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      NGL1=GSET%NGL
      ALLOCATE(GVECL(3,NGL1))
      CALL PLANEWAVE$GETR8A('GVEC',3*NGL1,GVECL)
      CALL PLANEWAVE$GETI4('NGG',NGG1)
      IF(THISTASK.NE.1) NGG1=1 
      ALLOCATE(GVECG1(3,NGG1))
      CALL PLANEWAVE$COLLECTR8(3,NGL1,GVECL,NGG1,GVECG1)
      DEALLOCATE(GVECL)
!
!     == GET GLOBAL GVECTORS FROM IKPT2 ================================
      CALL WAVES_SELECTWV(IKPT2,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      NGL2=GSET%NGL
      ALLOCATE(GVECL(3,NGL2))
      CALL PLANEWAVE$GETR8A('GVEC',3*NGL2,GVECL)
      CALL PLANEWAVE$GETI4('NGG',NGG2)
      IF(THISTASK.NE.1) NGG2=1 
      ALLOCATE(GVECG2(3,NGG2))
      CALL PLANEWAVE$COLLECTR8(3,NGL2,GVECL,NGG2,GVECG2)
      DEALLOCATE(GVECL)
!
!     == COPY GBAS      ================================================
      CALL WAVES_SELECTWV(IKPT2,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      CALL PLANEWAVE$GETR8A('GBAS',9,GBAS)
      CALL WAVES_SELECTWV(IKPT2,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      CALL PLANEWAVE$SETR8A('GBAS',9,GBAS)
!
!     == DEFINE MAPPING ================================================
      IF(THISTASK.EQ.1) THEN
        ALLOCATE(MAPG(NGG1))
        CALL WAVES_MAPG(NGG2,GBAS,GVECG2,NGG1,GVECG1,MAPG)
      END IF
      DEALLOCATE(GVECG1)
      DEALLOCATE(GVECG2)
!
!     ==================================================================
!     ==  COPY                                                        ==
!     ==================================================================
      CALL WAVES_SELECTWV(IKPT1,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      CALL PLANEWAVE$GETL4('TINV',TSUPER1)
      NB1=THIS%NB
      NBH1=THIS%NBH
      CALL WAVES_SELECTWV(IKPT2,1)
      CALL PLANEWAVE$SELECT(GSET%ID)
      CALL PLANEWAVE$GETL4('TINV',TSUPER2)
      NB2=THIS%NB
      NBH2=THIS%NBH

      NBN=MIN(NB1,NB2)
      IF(TSUPER1) THEN
        NBHN=MIN(NBN/2,NBH1)
      ELSE
        NBHN=NBN
      END IF
      ALLOCATE(PSI1(NGL1,NDIM,NB1))
      ALLOCATE(PSI2(NGL2,NDIM,NB2))
      ALLOCATE(PSIG1(NGG1))
      ALLOCATE(PSIG2(NGG2))
      ALLOCATE(PSITMP2(NGL2,2))
!
!     ==================================================================
!     == MAP DATA ONTO TEMP ARRAY                                     ==
!     ==================================================================
      DO IWAVE=1,2
        DO ISPIN=1,NSPIN
!
!         ==============================================================
!         == MAP DATA ONTO TEMP ARRAY                                 ==
!         ==============================================================
          CALL WAVES_SELECTWV(IKPT2,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(IWAVE.EQ.1) THEN
            DO IB=1,NBH2
              PSI2(:,:,IB)=THIS%PSI0(:,:,IB)
            ENDDO
          ELSE IF(IWAVE.EQ.2) THEN
            DO IB=1,NBH2
              PSI2(:,:,IB)=THIS%PSIM(:,:,IB)
            ENDDO
          END IF
!
!         ==============================================================
!         ==  EXPAND TO REGULAR WAVE FUNCTIONS                        ==
!         ==============================================================
          IF(TSUPER2) THEN
            DO IB=NBH2,1,-1
              IB1=2*IB-1
              IB2=2*IB
              DO IDIM=1,NDIM
                PSITMP2(:,1)=PSI2(:,IDIM,IB)
                CALL PLANEWAVE$INVERTG(NGL2,PSI2(1,IDIM,IB),PSITMP2(1,2))
                PSI2(:,IDIM,IB1)= 0.5D0   *(PSITMP2(:,1)+PSITMP2(:,2))
                PSI2(:,IDIM,IB2)=-0.5D0*CI*(PSITMP2(:,1)-PSITMP2(:,2))
              END DO
            ENDDO
          END IF
!
!         ==============================================================
!         ==  MAP ONTO IKPT1                                          ==
!         ==============================================================
          DO IB=1,NBN
            DO IDIM=1,NDIM
              CALL WAVES_SELECTWV(IKPT2,ISPIN)
              CALL PLANEWAVE$SELECT(GSET%ID)
              CALL PLANEWAVE$COLLECTC8(1,NGL2,PSI2(1,IDIM,IB),NGG2,PSIG2)
              IF(THISTASK.EQ.1) THEN
                DO IG=1,NGG1
                  IF(MAPG(IG).NE.0) THEN
                    PSIG1(IG)=PSIG2(MAPG(IG))
                  ELSE
                    PSIG1(IG)=(0.D0,0.D0)
                  END IF
                ENDDO
              END IF
              CALL WAVES_SELECTWV(IKPT1,ISPIN)
              CALL PLANEWAVE$SELECT(GSET%ID)
              CALL PLANEWAVE$DISTRIBUTEC8(1,NGG1,PSIG1,NGL1,PSI1(1,IDIM,IB))
            ENDDO
          ENDDO
!
!         ==============================================================
!         == MAP ONTO SUPER WAVE FUNCTIONS                            ==
!         ==============================================================
          IF(TSUPER1) THEN
            DO IB=1,NBHN
              IB1=2*IB-1
              IB2=2*IB
              PSI1(:,1,IB)=PSI1(:,1,IB1)+CI*PSI1(:,1,IB2)
            ENDDO
          END IF
!
!         ==============================================================
!         == SAVE INTO THIS                                           ==
!         ==============================================================
          CALL WAVES_SELECTWV(IKPT1,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(IWAVE.EQ.1) THEN
            DO IB=1,NBHN
              THIS%PSI0(:,:,IB)=PSI1(:,:,IB)
            ENDDO
          ELSE IF(IWAVE.EQ.2) THEN
            DO IB=1,NBHN
              THIS%PSIM(:,:,IB)=PSI1(:,:,IB)
            ENDDO
          END IF
        ENDDO
      ENDDO
      IF(THISTASK.EQ.1)DEALLOCATE(MAPG)
      DEALLOCATE(PSI1)
      DEALLOCATE(PSI2)
      DEALLOCATE(PSIG1)
      DEALLOCATE(PSIG2)
      DEALLOCATE(PSITMP2)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_COPYLAMBDA(IKPT1,IKPT2)
!     ******************************************************************
!     **                                                              **
!     **  COPY LAGRANGE MULTIPLIERS FROM K-POINT IKPT1 TO IKPT2       **
!     **                                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: IKPT1
      INTEGER(4)   ,INTENT(IN) :: IKPT2
      INTEGER(4)               :: NB1,NB2,NBA
      INTEGER(4)               :: ISPIN
      COMPLEX(8)   ,ALLOCATABLE:: LAMBDA(:,:)
!     ******************************************************************
      CALL WAVES_SELECTWV(IKPT1,1)
      NB1=THIS%NB
      CALL WAVES_SELECTWV(IKPT2,1)
      NB2=THIS%NB
      NBA=MIN(NB1,NB2)
      ALLOCATE(LAMBDA(NBA,NBA))
      DO ISPIN=1,NSPIN
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM0)) RETURN
        LAMBDA=THIS%RLAM0(1:NBA,1:NBA)
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM0)) ALLOCATE(THIS%RLAM0(NB2,NB2))
        THIS%RLAM0=(0.D0,0.D0)        
        THIS%RLAM0(1:NBA,1:NBA)=LAMBDA
      ENDDO
      DO ISPIN=1,NSPIN
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAMM)) RETURN
        LAMBDA=THIS%RLAMM(1:NBA,1:NBA)
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAMM))ALLOCATE(THIS%RLAMM(NB2,NB2))
        THIS%RLAMM=(0.D0,0.D0)        
        THIS%RLAMM(1:NBA,1:NBA)=LAMBDA
      ENDDO
      DO ISPIN=1,NSPIN
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM2M)) RETURN
        LAMBDA=THIS%RLAM2M(1:NBA,1:NBA)
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM2M)) ALLOCATE(THIS%RLAM2M(NB2,NB2))
        THIS%RLAM2M=(0.D0,0.D0)        
        THIS%RLAM2M(1:NBA,1:NBA)=LAMBDA
      ENDDO
      DO ISPIN=1,NSPIN
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM3M)) RETURN
        LAMBDA=THIS%RLAM3M(1:NBA,1:NBA)
        CALL WAVES_SELECTWV(IKPT2,ISPIN)
        IF(.NOT.ASSOCIATED(THIS%RLAM3M)) ALLOCATE(THIS%RLAM3M(NB2,NB2))
        THIS%RLAM2M=(0.D0,0.D0)        
        THIS%RLAM3M(1:NBA,1:NBA)=LAMBDA
      ENDDO
      DEALLOCATE(LAMBDA)
      RETURN
      END SUBROUTINE WAVES_COPYLAMBDA!
!      .................................................................
       SUBROUTINE WAVES_SPINOROVERLAP(NBH,NB,IKPT,QMAT)
!      *****************************************************************
!      **  CALCULATES Q_K,I,L,J=1/2*<PSI_I,K|PSI_J,L>                 **
!      **  I AND J ARE BAND-INDICES                                   **
!      **  K AND L ARE 1 OR 2 AND ARE THE SPINOR PART OF THE WAVEFUNCTION
!      *****************************************************************
       USE WAVES_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NBH,NB,IKPT
       COMPLEX(8),INTENT(OUT):: QMAT(2*NB*NSPIN,2*NB*NSPIN)
       COMPLEX(8),ALLOCATABLE :: AUXMAT(:,:)
       COMPLEX(8),ALLOCATABLE :: AUXMAT2(:,:)
       COMPLEX(8),ALLOCATABLE :: OPROJ(:,:,:)
       COMPLEX(8),ALLOCATABLE :: PSI0(:,:,:)
       INTEGER(4)            :: NGL,NBD,I,J,NDIMHALF,ISPIN
       IF (NDIM.EQ.1.AND.NSPIN.EQ.1) THEN
          CALL ERROR$MSG('S^2 ONLY POSSIBLE FOR NOT SPIN RESTRICTED CALCULATION')
          CALL ERROR$I4VAL('NDIM',NDIM)
          CALL ERROR$I4VAL('NSPIN',NSPIN)
          CALL ERROR$STOP('WAVES_SPINOROVERLAP')
       END IF
       IF(NDIM.EQ.2) THEN   !NON-COLLINEAR
          CALL WAVES_SELECTWV(IKPT,1)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          NGL=GSET%NGL   
          NBD=2*NB
          NDIMHALF=1
          ALLOCATE(AUXMAT(NBD,NBD))
          CALL WAVES_1COVERLAP(.TRUE.,MAP,NDIMHALF,NBD,NBD,MAP%NPRO,THIS%PROJ,THIS%PROJ,AUXMAT)
          CALL WAVES_OVERLAP(.TRUE.,NGL,NDIMHALF,NBD,NBD,THIS%PSI0,THIS%PSI0,QMAT)
          DO J=1,NBD
             DO I=1,NBD
                QMAT(I,J)=(QMAT(I,J)+AUXMAT(I,J))*0.5D0
             ENDDO
          ENDDO
          DEALLOCATE(AUXMAT)
       ELSE  ! COLLINEAR SPIN POLARIZED
          ALLOCATE(AUXMAT(NB*2,NB*2))
          DO ISPIN=1,NSPIN
             CALL WAVES_SELECTWV(IKPT,ISPIN)
             CALL PLANEWAVE$SELECT(GSET%ID)      
             NGL=GSET%NGL              
             NBD=2*NB
             IF(ISPIN.EQ.1) ALLOCATE(OPROJ(1,NBH*2,MAP%NPRO))
             IF(ISPIN.EQ.1) THEN
                OPROJ(:,1:NBH,:)=THIS%PROJ
             ELSE
                OPROJ(:,NBH+1:2*NBH,:)=THIS%PROJ
             END IF
          END DO
          CALL WAVES_1COVERLAP(.TRUE.,MAP,1,NBH*2,NB*2,MAP%NPRO,OPROJ,OPROJ,AUXMAT)
          DEALLOCATE(OPROJ)
          ALLOCATE(AUXMAT2(NB*2,NB*2))
          DO ISPIN=1,NSPIN
             CALL WAVES_SELECTWV(IKPT,ISPIN)
             CALL PLANEWAVE$SELECT(GSET%ID)      
             NGL=GSET%NGL
             IF(ISPIN.EQ.1) ALLOCATE(PSI0(NGL,1,2*NBH))
             IF(ISPIN.EQ.1) THEN
                PSI0(:,:,1:NBH)=THIS%PSI0
             ELSE
                PSI0(:,:,NBH+1:2*NBH)=THIS%PSI0
             END IF
          END DO
          CALL WAVES_OVERLAP(.TRUE.,NGL,NDIM,NBH*2,NB*2,PSI0,PSI0,AUXMAT2)
          DEALLOCATE(PSI0)

          QMAT=(0.D0,0.D0)
          DO I=1,NB  ! UP UP
             DO J=1,NB
                QMAT(2*I-1,2*J-1)=(AUXMAT(I,J)+AUXMAT2(I,J))*0.5D0
             ENDDO
          ENDDO
          DO I=NB+1,2*NB  ! DOWN DOWN
             DO J=NB+1,2*NB
                QMAT(2*I,2*J)=(AUXMAT(I,J)+AUXMAT2(I,J))*0.5D0
             ENDDO
          ENDDO
          DO I=1,NB  ! UP DOWN
             DO J=NB+1,2*NB
                QMAT(2*I-1,2*J)=(AUXMAT(I,J)+AUXMAT2(I,J))*0.5D0
             ENDDO
          ENDDO         
          DO I=NB+1,2*NB  ! DOWN UP
             DO J=1,NB
                QMAT(2*I,2*J-1)=(AUXMAT(I,J)+AUXMAT2(I,J))*0.5D0
             ENDDO
          ENDDO
          DEALLOCATE(AUXMAT)
          DEALLOCATE(AUXMAT2)
       END IF
       RETURN
     END SUBROUTINE WAVES_SPINOROVERLAP
!
!.....................................................................
module wavesfixrho_module
! this module is used to keep the density fixed
  real(8), allocatable    :: qlm(:,:)
  real(8), allocatable    :: rho(:,:)
  real(8), allocatable    :: denmat(:,:,:,:)
end module wavesfixrho_module
!
!      .................................................................
       subroutine waves_fixrhoget(nrl,ndimd,lmrxx,nat,qlm_,rho_,denmat_)
       use wavesfixrho_module 
       implicit none
       integer(4) ,intent(in)  :: nrl,ndimd,lmrxx,nat
       real(8)    ,intent(out) :: qlm_(lmrxx,nat),rho_(nrl,ndimd)
       real(8)    ,intent(out) :: denmat_(lmrxx,lmrxx,ndimd,nat)
!      ******************************************************************
       if (.not.allocated(qlm)) then
          allocate(qlm(lmrxx,nat))
          allocate(rho(nrl,ndimd))
          allocate(denmat(lmrxx,lmrxx,ndimd,nat))
          call waves$fixrhoread()
       end if
       qlm_(:,:)=qlm(:,:)
       rho_(:,:)=rho(:,:)
       denmat_(:,:,:,:)=denmat(:,:,:,:)
       end
!
!      .................................................................
       subroutine waves_fixrhoset(nrl,ndimd,lmrxx,nat,qlm_,rho_,denmat_)
       use wavesfixrho_module 
       implicit none
       integer(4) ,intent(in)  :: nrl,ndimd,lmrxx,nat
       real(8)    ,intent(out) :: qlm_(lmrxx,nat),rho_(nrl,ndimd)
       real(8)    ,intent(out) :: denmat_(lmrxx,lmrxx,ndimd,nat)
!      ******************************************************************
       if (.not.allocated(qlm)) then
          allocate(qlm(lmrxx,nat))
          allocate(rho(nrl,ndimd))
          allocate(denmat(lmrxx,lmrxx,ndimd,nat))
       end if
       qlm(:,:)=qlm_(:,:)
       rho(:,:)=rho_(:,:)
       denmat(:,:,:,:)=denmat_(:,:,:,:)
       end
!
!      .................................................................
       subroutine waves$fixrhoread()
       use wavesfixrho_module
       implicit none
       integer(4)                 :: nfil
!      ******************************************************************
       CALL FILEHANDLER$SETFILE('FIXRHO',.FALSE.,'fixrho.bin')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','FORM','UNFORMATTED')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','POSITION','REWIND')
       CALL FILEHANDLER$UNIT('FIXRHO',NFIL)
       read(nfil)qlm(:,:)
       read(nfil)rho(:,:)
       read(nfil)denmat(:,:,:,:)
       call filehandler$close('FIXRHO')
       return
       end
!
!      .................................................................
       subroutine waves$fixrhowrite()
       use wavesfixrho_module
       implicit none
       integer(4)                 :: nfil
!      ******************************************************************
       CALL FILEHANDLER$SETFILE('FIXRHO',.FALSE.,'fixrho.bin')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','FORM','UNFORMATTED')
       CALL FILEHANDLER$SETSPECIFICATION('FIXRHO','POSITION','REWIND')
       CALL FILEHANDLER$UNIT('FIXRHO',NFIL)
       write(nfil)qlm(:,:)
       write(nfil)rho(:,:)
       write(nfil)denmat(:,:,:,:)
       call filehandler$close('FIXRHO')
       return
       end
       
       

!
!.....................................................................
MODULE TOTALSPIN_MODULE
REAL(8)               :: TOTSPIN(4) ! DIFFERS FROM JOHANNES' VERSION
END MODULE TOTALSPIN_MODULE
!
!      .................................................................
       SUBROUTINE WAVES_TOTALSPIN(NB,NKPT,IKPT,NSPIN,OCC,QMAT)
!      *****************************************************************
!      **  CALCULATES <S^2> IN UNITS OF HBAR^2                        **
!      **                                                             **
!      **                                                             **
!      *****************************************************************
       USE MPE_MODULE
       USE TOTALSPIN_MODULE, ONLY: TOTSPIN ! DIFFERS FROM JOHANNES' VERSION
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NB,NKPT,IKPT,NSPIN
       REAL(8)   ,INTENT(IN) :: OCC(NB,NKPT,NSPIN) 
       COMPLEX(8),INTENT(IN) :: QMAT(2,NB*NSPIN,2,NB*NSPIN) ! Q_I,J,K,L=1/2*<PSI_I,K|PSI_J,K>  
       COMPLEX(8)            :: SUM,PART,SPINX,SPINY,SPINZ
       REAL(8)               :: IOCC(NB*NSPIN),SVAR
       INTEGER(4)            :: I,J,NTASKS,THISTASK,NBD
       
       CALL MPE$QUERY(NTASKS,THISTASK)

       DO I=1,NB
          IOCC(I)=OCC(I,IKPT,1)
       END DO
       IF (NSPIN.EQ.2) THEN
          DO I=NB+1,2*NB
             IOCC(I)=OCC(I-NB,IKPT,2)
          END DO
       END IF
       NBD=NSPIN*NB !ONE BAND FOR ONE ELECTRON

       SUM=0.D0
       DO I=1,NBD
             PART=QMAT(1,I,2,I)+QMAT(2,I,1,I)
             SUM=SUM+PART*IOCC(I)
       END DO
       SPINX=SUM

       SUM=0.D0
       DO I=1,NBD
             PART=QMAT(1,I,2,I)-QMAT(2,I,1,I)
             SUM=SUM+PART*IOCC(I)
       END DO
       SPINY=SUM/(0,1)

       SUM=0.D0
       DO I=1,NBD
             PART=QMAT(1,I,1,I)-QMAT(2,I,2,I)
             SUM=SUM+PART*IOCC(I)
       END DO
       SPINZ=SUM


       SUM=0.D0
       DO I=1,NBD
          SUM=SUM+IOCC(I)*0.75D0
       END DO
       
       SUM=SUM+SPINX**2+SPINY**2+SPINZ**2
       PRINT*,'TOTALSPIN: PART 3/4+X+Y+Z: ',SUM
       
       SVAR=SUM
       DO I=1,NBD
          DO J=1,NBD
             PART=-(QMAT(1,I,1,J)-QMAT(2,I,2,J))*(QMAT(1,J,1,I)-QMAT(2,J,2,I))
             SUM=SUM+PART*SQRT(IOCC(I)*IOCC(J))
          END DO
       END DO
       PRINT*,'TOTALSPIN: PART UP UP    : ',SUM-SVAR

       SVAR=SUM
       DO I=1,NBD
          DO J=1,NBD
             PART=-4.D0*QMAT(1,I,2,J)*QMAT(2,J,1,I)
             SUM=SUM+PART*SQRT(IOCC(I)*IOCC(J))
          END DO
       END DO       
       PRINT*,'TOTALSPIN: PART DOWN UP  : ',SUM-SVAR

       IF (THISTASK.EQ.1) THEN
          PRINT*,'SPIN: 2*<S_X>: ',2.D0*REAL(SPINX),' HBAR'
          PRINT*,'SPIN: 2*<S_Y>: ',2.D0*REAL(SPINY),' HBAR'
          PRINT*,'SPIN: 2*<S_Z>: ',2.D0*REAL(SPINZ),' HBAR'
          PRINT*,'TOTALSPIN: <S^2>: ',REAL(SUM),' HBAR^2'
          PRINT*,'SPIN QUATUM NUMBER S=',-0.5+SQRT(0.25+REAL(SUM))
          IF (IKPT.EQ.1) TOTSPIN=0.D0 ! SUM UP TOTAL SPIN OVER K-POINTS
          TOTSPIN(1)=TOTSPIN(1)+DBLE(SUM)
          TOTSPIN(2)=TOTSPIN(2)+DBLE(SPINX)
          TOTSPIN(3)=TOTSPIN(3)+DBLE(SPINY)
          TOTSPIN(4)=TOTSPIN(4)+DBLE(SPINZ)
       END IF
       RETURN
     END SUBROUTINE WAVES_TOTALSPIN
!
!     ..................................................................
      SUBROUTINE WAVES$REPORTSPIN(NFIL)
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE TOTALSPIN_MODULE, ONLY: TOTSPIN ! DIFFERS FROM JOHANNES' VERSION
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)               :: SVAR
!     ******************************************************************
      IF (NSPIN.EQ.1.AND.NDIM.EQ.1) THEN
         RETURN
      ELSE
         CALL REPORT$TITLE(NFIL,'TOTAL SPIN ANALYSIS')
         IF (NDIM.EQ.2) CALL REPORT$R8VAL(NFIL,'<S_X>',TOTSPIN(2),'HBAR')
         IF (NDIM.EQ.2) CALL REPORT$R8VAL(NFIL,'<S_Y>',TOTSPIN(3),'HBAR')
         CALL REPORT$R8VAL(NFIL,'<S_Z>',TOTSPIN(4),'HBAR')
         SVAR=DSQRT(TOTSPIN(2)**2+TOTSPIN(3)**2+TOTSPIN(4)**2)
         IF (NDIM.EQ.2) CALL REPORT$R8VAL(NFIL,'|S|',SVAR,'HBAR')
!        CALL REPORT$R8VAL(NFIL,'TOTAL SPIN <S^2>',TOTSPIN(1),'HBAR^2')
!        SVAR=-0.5+SQRT(0.25+TOTSPIN(1))
!        CALL REPORT$R8VAL(NFIL,'SPIN QUANTUM NUMBER S',SVAR,'')
         RETURN
      END IF
    END SUBROUTINE WAVES$REPORTSPIN
!
!     ....................................................................
      subroutine waves_comparepsi(id,ngl,ndim,nbh,psi1,psi2)
      USE WAVES_MODULE, only : gset,delt
      implicit none
      character(*)    ,intent(in) :: id
      INTEGER(4)      ,INTENT(IN) :: NGL
      INTEGER(4)      ,INTENT(IN) :: NDIM
      INTEGER(4)      ,INTENT(IN) :: NBH
      COMPLEX(8)      ,INTENT(IN) :: PSI1(NGL,NDIM,NBH)
      COMPLEX(8)      ,INTENT(IN) :: PSI2(NGL,NDIM,NBH)
      complex(8)                  :: csvar
      integer(4)                  :: ig,idim,ib
      real(8)                     :: sum,sumtot
      real(8)                     :: rbas(3,3),gbas(3,3),cellvol
!     *********************************************************************
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL GBASS(RBAS,GBAS,CELLVOL)
      sumtot=0.d0
      do ib=1,nbh
        sum=0.d0
        do ig=1,ngl
          do idim=1,ndim
            csvar=psi1(ig,idim,ib)-psi2(ig,idim,ib)
            sum=sum+GSET%MPSI(IG)*(real(csvar)**2+aimag(csvar)**2)
          enddo
        enddo
        sumtot=sumtot+sum
        write(*,*)' kineticenergytest ',ib,sum*CELLVOL/DELT**2,id
      enddo
        write(*,*)' kineticenergytest total',sumtot*CELLVOL/DELT**2,id
      return
      end
