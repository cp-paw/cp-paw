!TODO:
!  SET EPWRHO IN WAVES AND PASS THROUGH FROM WAVES INTO POTENTIAL
!  IF REQUIRED
!
!......................................................WAVES............
MODULE WAVES_MODULE                                                
!***********************************************************************
!**                                                                   **
!**  NAME: WAVES                                                      **
!**                                                                   **
!**  PURPOSE: OPERATIONS ON THE PLANE WAVE PART OF THE WAVE FUNCTIONS **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    WAVES$SETR8A(ID,LEN,VAL)                                       **
!**    WAVES$GETR8A(ID,LEN,VAL)                                       **
!**    WAVES$GVECTORS                                                 **
!**    WAVES$ETOT                                                     **
!**    WAVES$PROPAGATE                                                **
!**    WAVES$CONSTRAINTS                                              **
!**    WAVES$REPORT(NFIL)                                             **
!
!**    WAVES$INITIALIZERANDOM                                         **
!**    WAVES$RANDOMIZEVELOCITY                                        **
!**    WAVES$STOP                                                     **
!**    WAVES$EIGR(STRING)                                             **
!**    WAVES$EKIN(EKIN)                                               **
!**    WAVES$DENSITY(NR,NSPIN,RHO)                                    **
!**    WAVES$HPSI(NR,NSPIN,POT)                                       **
!**    ...                                                            **
!**                                                                   **
!**  DEPENDECIES:                                                     **
!**    ERROR                                                          **
!**    MPELIB                                                         **
!**    PLANEWAVE                                                      **
!**    SETUP                                                          **
!**    ...                                                            **
!**                                                                   **
!**   REMARKS                                                         **
!**   1) HPSI IS REUSED IN WAVES$PROPAGATE TO HOLD PSI(-).            **
!**     NEEDED TO GET THE PROPER KINETIC ENERGY OF THE WAVE FUNCTIONS **
!**                                                                   **
!**     NGW(IKPT)    NUMBER OF PLANE WAVES                            **
!**     NGWL(IKPT)   NUMBER OF PLANE WAVES ON LOCAL PROCESSOR         **
!**     NGWX         MAX.(ALL K-POINTS) NUMBER OF PLANE WAVES         **
!**     NGWLX        MAX.(ALL K-POINTS) NUMBER OF PLANE WAVES         **
!**                  ON LOCAL PROCESSOR                               **
!**     NR1,NR2,NR3  GRID DIMENSIONS FOR THE DENSITY                  **
!**     NR1L         NUMBER OF PLANES ON THIS PROCESSOR               **
!**     NR1START     FIRST PLANE ON LOCAL PROCESSOR                   **
!**     NRSTORE      FIRST DIMENSION FOR THE REAL SPACE WAVE          **
!**                  WAVE FUNCTIONS, WHICH ARE OVERLAYED WITH         **
!**                  ONE SET OF REC. SPACE WAVE FUNCTIONS             **
!**     NNR          NUMBER OF R-GRID POINTS                          **
!**     NNRL         NUMBER OF R-GRID POINTS ON LOCAL PROCESSOR       **
!**                                                                   **
!**                                                                   **
!
!  A TYPE GVECTORS_TYPE DESCRIBES A K-POINT AND ITS G-GRID
!  A TYPE PSIBLOCK_TYPE DESCRIBES THE WAVE FUNCTIONS 
!     AND POINTS TO A GVECTORS_TYPE
!
! REMARKS PARALLELIZE PSIPRO ACCORING TO ATOMS
! COMMUNICATE AFTER CALCULATING PROJECTIONS
! AND BEFORE AQDDING PROJECTORS
!
! SET EPWRHO FOR WAVES INSTEAD OF POTENTIAL
!
!
!******************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1995)**
USE LINKEDLIST_MODULE
TYPE MAP_TYPE
  INTEGER(4)         :: NAT         ! #(ATOMS)
  INTEGER(4)         :: NSP         ! #(ATOM TYPES)
  INTEGER(4)         :: NBAREPRO    ! #(LN*NAT)
  INTEGER(4)         :: NPRO        ! #(LMN*NAT)
  INTEGER(4),POINTER :: LNX(:)      !(NSP)
  INTEGER(4),POINTER :: LMNX(:)     !(NSP)
  INTEGER(4),POINTER :: LOX(:,:)    !(LNX,NSP) MAIN ANGULAR MOMENTUM OF AN LN
  INTEGER(4),POINTER :: ISP(:)      !(NAT) POINTER TO ATOM TYPES
  INTEGER(4),POINTER :: IBAREPRO(:) !(NSP) POINTER TO FIRST LN OF AN ATOM TYPE
  INTEGER(4),POINTER :: IPRO(:)     !(NAT) POINTER TO FIRST LMN OF AN ATOM
END TYPE MAP_TYPE
TYPE GSET_TYPE 
  CHARACTER(16)      :: ID         ! GSET IDENTIFIER
  REAL(8)            :: WKPT       ! K-POINT WEIGHT
  REAL(8)            :: XK(3)      !(3,NKPT) K-POINT RELATIVE COORDINATES
  LOGICAL(4)         :: TINV       ! SWITCH: TIME INVERSION SYMMETRY
  INTEGER(4)         :: NGL        !#(PLANE WAVES(LOCAL))
  REAL(8)   ,POINTER :: PRO(:,:)   !(NGL,NBAREPRO) PROJECTOR
  REAL(8)   ,POINTER :: YLM(:,:)   !(NGL,LMXX) REAL SPHERICAL HARMONICS
  REAL(8)   ,POINTER :: EIGR(:,:)  !(NGL,NAT) STRUCTURE FACTORS
END TYPE GSET_TYPE
TYPE WVSET_TYPE  !======================================================
  TYPE(GSET_TYPE),POINTER :: GSET
  INTEGER(4)         :: NBH
  COMPLEX(8),POINTER :: PSI0(:,:,:)     !(NGL,NBH,IDIM)  PSPSI(0)
  COMPLEX(8),POINTER :: PSIM(:,:,:)     !(NGL,NBH,IDIM)  PSPSI(-,+)(G)
  COMPLEX(8),POINTER :: PSIR(:,:,:)     !(NRL,NBH,IDIM)  PSPSI(-,+)(R)
  COMPLEX(8),POINTER :: PROJ(:,:,:)     !(NPRO,NBH,IDIM) <PSPSI|P>
  COMPLEX(8),POINTER :: DPROJ(:,:,:,:)  !(NPRO,NBH,3,IDIM) GRAD<PSPSI|P>
  COMPLEX(8),POINTER :: HPROJ(:,:,:)    !(NPRO,NBH,IDIM)
  COMPLEX(8),POINTER :: OPROJ(:,:,:)    !(NPRO,NBH,IDIM)
  COMPLEX(4),POINTER :: STORE(:,:,:)    !(NRSTORE,NB,IDIM)
                                        ! +(WAVES$DENSITY) -(WAVES$HPSI)
  COMPLEX(8),POINTER :: HPSI(:,:,:)     !(NGWLX,NB,IDIM)
                                        ! +(WAVES$HPSI)-(WAVES$PROPAGATE)
  COMPLEX(8),POINTER :: OPSI(:,:,:)     !(NGWLX,NB,IDIM)
                                        ! +(WAVES$OPSI) -(WAVES$ORTHOGONALIZE)
  REAL(8)   ,POINTER :: RLAMP(:,:)      !(NB,NB)
  REAL(8)   ,POINTER :: RLAM0(:,:)      !(NB,NB)
  REAL(8)   ,POINTER :: RLAMM(:,:)      !(NB,NB)
  REAL(8)   ,POINTER :: RLAM2M(:,:)     !(NB,NB)
  REAL(8)   ,POINTER :: RLAM3M(:,:)     !(NB,NB)
  REAL(8)   ,POINTER :: HAME(:,:)        !(NB)
  REAL(8)   ,POINTER :: HAMU(:,:)       !(NB,NB)
END TYPE WVSET_TYPE
TYPE EXTERNALPOINTER_TYPE !=============================================
  INTEGER(4) :: IB
  INTEGER(4) :: IKPT
  INTEGER(4) :: ISPIN
  INTEGER(4) :: IAT
END TYPE EXTERNALPOINTER_TYPE
!========================================================================
!== DATA                                                               ==
!========================================================================
REAL(8) ,PARAMETER :: EPSILONGAMMA=1.D-7
!========================================================================
!== DATA THAT DESCRIBE FUNCTIONALITY OF WAVES OBJECT                   ==
!==  THESE DATA HAVE TO BE SET EXPLICITELY                             ==
!========================================================================
INTEGER(4)  :: NKPT          ! #(K-POINTS)
INTEGER(4)  :: NSPIN         ! #(SPINS)
INTEGER(4)  :: NDIM          ! #(WAVE FUNCTION COMPONENTS)
INTEGER(4)  :: NCOMP         ! #(DENSITY COMPONENTS) (NSPIN OR (NDIM+1)*NDIM/2)
INTEGER(4)  :: NB            ! #(BANDS)
REAL(8)     :: EPWPSI=0.D0   ! PLANE WAVE CUTOFF GMAX**2/2
REAL(8)     :: EPWRHO=0.D0   ! PLANE WAVE CUTOFF GMAX**2/2
REAL(8)     :: EMASS=0.D0    ! MPSI(G)=EMASS*(1+EMASSCG2*G2)
REAL(8)     :: EMASSCG2=0.D0 ! MPSI(G)=EMASS*(1+EMASSCG2*G2)
REAL(8)     :: ANNEE=0.D0    ! FRICTION 
LOGICAL(4)  :: TSTOP=.FALSE. ! INITIAL VELOCITY SET TO ZERO
LOGICAL(4)  :: TSAFEORTHO=.TRUE.  ! CHOICE OR ORTHOGONALIZATION
LOGICAL(4)  :: TRANDOMIZE=.FALSE. ! RANDOMIZE INITIAL WAVE FUNCTIONS
REAL(8)     :: AMPRE=0.D0    ! KINETIC ENERGY FOR INITIAL RANDOMIZATION
LOGICAL     :: TSTORE=.TRUE. ! STORE REAL SPACE WAVE FUNCTIONS
REAL(8)     :: DELT=0.D0     ! TIME STEP IN A.U.
LOGICAL(4)  :: THAMILTON=.FALSE.    ! HAMILTON MATRIX AVAILABLE
LOGICAL(4)  :: TRAWSTATES=.FALSE.   ! PROVIDES NON-DIAGONALIZED WAVE FUNCTIONS THROUGH $GETR
!========================================================================
!== PERMANENT DATA, WHICH ARE ORGANIZED BY THE ROUTINES ITSELF         == 
!========================================================================
INTEGER(4)  :: NRSTORE
INTEGER(4)  :: NRL           !#(R-SPACE GRID POINTS(LOCAL))
INTEGER(4)  :: NR1,NR2,NR3   ! GLOBAL R-SPACE GRID
INTEGER(4)  :: NR1L,NR1START ! LOCAL SLAB OF R-SPACE GRID
REAL(8)     :: WKPT          ! K-POINT WEIGHT
TYPE(WVSET_TYPE),POINTER :: THISARRAY(:,:)   ! (NKPT,NSPIN)
TYPE(WVSET_TYPE),POINTER :: THIS            ! CURRENT SET OF WAVES
TYPE(GSET_TYPE) ,POINTER :: GSET            ! CURRENT SET OF GSET
TYPE(MAP_TYPE)           :: MAP
LOGICAL(4)               :: TPR=.FALSE.
CONTAINS
!***********************************************************************
      SUBROUTINE WAVES_SELECTWV(IKPT,ISPIN)
      IF(.NOT.ASSOCIATED(THISARRAY)) THEN
        CALL ERROR$MSG('THISARRAY DOES NOT EXIST') 
        CALL ERROR$STOP('WAVES_SELECTWV')
      END IF 
      THIS=>THISARRAY(NKPT,NSPIN)
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
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('WAVES$EPWRHO')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$GVECTORS
!     ******************************************************************
!     **  WAVES$GVECTORS                                              **
!     **  GENERATE G-VECTORS AND OTHER INITIALIZATION                 **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL(4)             :: TCHK
      REAL(8)                :: RBAS(3,3) ! LATTICE VECTORS
      REAL(8)                :: GBAS(3,3) ! RECIPROCAL LATTICE VECTORS
      REAL(8)                :: CELLVOL   ! UNIT CELL  VOLUME
      REAL(8)   ,ALLOCATABLE :: XK(:,:)   ! K-POINTS IN RELATIVE COORDINATES
      REAL(8)   ,ALLOCATABLE :: GVEC(:,:) ! G-VECTORS IN CARTESIAN COORDINATES
      REAL(8)   ,ALLOCATABLE :: G2(:)     ! G**2
      INTEGER(4)             :: LN,ISP,IKPT,ISPIN,IAT
      INTEGER(4)             :: NBH
      INTEGER(4)             :: NGL
      INTEGER(4)             :: LNXX
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$GVECTORS')
!     
!     ==================================================================
!     ==  DEFINE REAL SPACE GRID                                      ==
!     ==================================================================
      CALL CELL$GETR8A('T',9,RBAS)
      CALL PLANEWAVE$DIVIDERGRIDONTASKS(EPWRHO,RBAS,NR1START,NR1L,NR2,NR3)
      CALL POTENTIAL$SETGRID(EPWRHO,RBAS,NR1START,NR1L,NR2,NR3)
!     
!     ==================================================================
!     ==  GET KPOINTS AND DIMENSIONS                                  ==
!     ==================================================================
      CALL DYNOCC$GETI4('NB',NB)
      CALL DYNOCC$GETI4('NKPT',NKPT) 
      CALL DYNOCC$GETI4('NSPIN',NSPIN)
      ALLOCATE(XK(3,NKPT))
      CALL DYNOCC$GETR8A('XK',3*NKPT,XK)
!     
!     ================================================================
!     ==  DETERMINE MAPPRO AND MAPBAREPRO                           ==
!     ================================================================
      CALL ATOMLIST$NATOM(MAP%NAT)
      CALL SETUP$NSPECIES(MAP%NSP)
      CALL SETUP$LNXX(LNXX)
      ALLOCATE(MAP%ISP(MAP%NAT))
      ALLOCATE(MAP%LNX(MAP%NSP))
      ALLOCATE(MAP%LMNX(MAP%NSP))
      ALLOCATE(MAP%LOX(LNXX,MAP%NSP))
      ALLOCATE(MAP%IBAREPRO(MAP%NSP))
      ALLOCATE(MAP%IPRO(MAP%NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,MAP%NAT,MAP%ISP)
      MAP%LOX(:,:)=0
      MAP%NBAREPRO=0
      DO ISP=1,MAP%NSP
        CALL SETUP$LNX(ISP,MAP%LNX(ISP))
        CALL SETUP$LOFLN(ISP,MAP%LNX(ISP),MAP%LOX(1,ISP))
        MAP%IBAREPRO(ISP)=MAP%NBAREPRO+1
        MAP%LMNX(ISP)=0
        DO LN=1,MAP%LNX(ISP)
          MAP%LMNX(ISP)=MAP%LMNX(ISP)+2*MAP%LOX(LN,ISP)
          MAP%NBAREPRO=MAP%NBAREPRO+1
        ENDDO
      ENDDO
      MAP%NPRO=0
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
        MAP%IPRO(IAT)=MAP%NPRO+1
        DO LN=1,MAP%LNX(ISP)
          MAP%NPRO=MAP%NPRO+2*MAP%LOX(LN,ISP)+1
        ENDDO
      ENDDO
!     
!     ================================================================
!     ==  ALLOCATE THISARRAY                                        ==
!     ================================================================
      ALLOCATE(THISARRAY(NKPT,NSPIN))
      DO IKPT=1,NKPT
        ALLOCATE(THISARRAY(IKPT,1)%GSET)
        DO ISPIN=2,NSPIN
          THISARRAY(IKPT,ISPIN)%GSET=>THISARRAY(IKPT,1)%GSET
        ENDDO
      ENDDO
!     
!     ================================================================
!     ==  INITIALIZE PLANE WAVE OBJECT                              ==
!     ================================================================
      CALL GBASS(RBAS,GBAS,CELLVOL)
      DO IKPT=1,NKPT
        CALL WAVES$SELECTWV(IKPT,1)
        WRITE(GSET%ID,FMT=*)IKPT
        GSET%ID='WAVE '//ADJUSTL(GSET%ID)
        CALL PLANEWAVE$INITIALIZE(GSET%ID,RBAS,XK,EPWPSI,NR1START,NR1L,NR2,NR3)
      ENDDO
!     
!     ==================================================================
!     ==  EVALUATE NUMBER OF G-VECTORS ETC. ON LOCAL PROCESSOR        ==
!     ==================================================================
      DO IKPT=1,NKPT
        CALL WAVES$SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
!       == K-POINT WEIGHT  (GENERAL K-POINTS ARE PAIRS OF G AND -G)
        GSET%WKPT=1.D0/DBLE(NKPT)
        IF(.NOT.GSET%TINV)WKPT=2.D0*WKPT
!       == GSPACE AND RSPACE DIMENSIONS ================================
        CALL PLANEWAVE$GETI4A('NGL',NGL)
        CALL PLANEWAVE$GETI4A('NRL',NRL)
        GSET%NGL=NGL
!       == #(DIMENSIONS) ===============================================
        NDIM=1    ! CAN BE GENERALIZED LATER TO SPINOR REPRESENTATION
!       == NUMBER OF DENSITY COMPONENTS
        IF(NSPIN.EQ.2) THEN
          NCOMP=NSPIN
        ELSE
          NCOMP=((NDIM+1)*NDIM)/2
        END IF
        DO ISPIN=1,NSPIN
          CALL WAVES$SELECTWV(IKPT,ISPIN)
          
!         == #(BANDS) ====================================================
          CALL PLANEWAVE$GETL4('TINV',GSET%TINV)
          IF(GSET%TINV) THEN
            NBH=(NB+1)/2
          ELSE
            NBH=NB
          END IF
          THIS%NBH=NBH
!         == ALLOCATE WAVE FUNCTIONS =====================================
          ALLOCATE(THIS%PSI0(NGL,NBH,NDIM))
          ALLOCATE(THIS%PSIM(NGL,NBH,NDIM))
!         == ALLOCATE PROJECTIONSONS =====================================
          ALLOCATE(THIS%PROJ(MAP%NPRO,NBH,NDIM))
          ALLOCATE(THIS%HPROJ(MAP%NPRO,NBH,NDIM))
          ALLOCATE(THIS%OPROJ(MAP%NPRO,NBH,NDIM))    
          ALLOCATE(THIS%DPROJ(MAP%NPRO,NBH,NDIM,3))
!         == ALLOCATE LAGRANGE MULTIPLIERS ===============================
          ALLOCATE(THIS%RLAMP(NB,NB))
          ALLOCATE(THIS%RLAM0(NB,NB)) 
          ALLOCATE(THIS%RLAMM(NB,NB))
          THIS%RLAMP(:,:)=0.D0
          THIS%RLAMM(:,:)=0.D0
          THIS%RLAM0(:,:)=0.D0
        ENDDO
      ENDDO
!     
!     ================================================================
!     ==  GET G-SPACE PROJECTOR FUNCTIONS                           ==
!     ================================================================
                           CALL TRACE$PASS('GPROJECTORS')
      DO IKPT=1,NKPT
        CALL WAVES$SELECTWV(IKPT,1)
        NGL=GSET%NGL
        ALLOCATE(GSET%PRO(NGL,MAP%NBAREPRO))
        ALLOCATE(GVEC(3,NGL))
        ALLOCATE(G2(NGL))
        CALL PLANEWAVE$GETR8A('GVEC',3*GSET%NGL,GVEC)
        CALL PLANEWAVE$GETR8A('G2',GSET%NGL,G2)
        DO ISP=1,MAP%NSP
          CALL SETUP$GPROJECTORS(ISP,MAP%LMNX(ISP),CELLVOL,NGL,NGL &
     &                          ,G2,GVEC,GSET%PRO)
        ENDDO
        DEALLOCATE(GVEC)
        DEALLOCATE(G2)
      ENDDO
      DEALLOCATE(XK)
!     
!     ================================================================
!     ==  SEND DATA TO OPTIC CODE                                   ==
!     ================================================================
!     CALL OPTIC3$LATTICE(RBAS,GBAS,CELLVOL,NGW,NGWX,NKPT,IND1T,IND2T,IND3T)
                           CALL TRACE$POP
      RETURN
      END
!
!     .............................................WAVES_PSEKIN..........
      SUBROUTINE WAVES$ETOT
!     ******************************************************************
!     **                                                              **
!     **  EVALUATE PS KINETIC ENERGY IN G-SPACE                       **
!     **  EVALUATE NUMBER OF ELECTRONS IN G-SPACE                     **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      USE WAVES_MODULE
      IMPLICIT NONE
      REAL(8)                :: F(NB,NKPT,NSPIN)
      COMPLEX(8)             :: CSUM
      REAL(8)   ,ALLOCATABLE :: QLM(:,:)  !(LMRXX,NAT) MULTIPOLE MOMENTS
      REAL(8)   ,ALLOCATABLE :: VQLM(:,:) !(LMRXX,NAT) MULTIPOLE POTENTIALS
      REAL(8)   ,ALLOCATABLE :: RHO(:,:)  ! CHARGE DENSITY
      REAL(8)   ,ALLOCATABLE :: RHO1(:,:) ! AUXILIARY DENSITY ARRAY
      COMPLEX(8),ALLOCATABLE :: PROJ(:,:,:)    ! <PRO|PSI0>
      COMPLEX(8),ALLOCATABLE :: dPROJ(:,:,:,:)   ! <PRO|PSI0>
      REAL(8)   ,ALLOCATABLE :: DENMAT(:,:,:,:)! ONE-CENTER DENSITY MATRIX
      REAL(8)   ,ALLOCATABLE :: DENMAT1(:,:,:) ! AUX. ONE-CENTER DENSITY MATRIX
      REAL(8)   ,ALLOCATABLE :: DATH(:,:,:) ! AUX.ONE-CENTER DENSITY MATRIX
      REAL(8)   ,ALLOCATABLE :: DO(:,:,:) ! AUX.ONE-CENTER DENSITY MATRIX
      real(8)                :: rhob
      INTEGER(4)             :: NGL
      INTEGER(4)             :: lmrx
      REAL(8)                :: EKIN,EKIN1
      INTEGER(4)             :: IKPT,ISPIN,IAT,ISP,ibh,ib,idim,ipro
      INTEGER(4)             :: LMN1,LMN2
      INTEGER(4)             :: LMNXX
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NBH
      INTEGER(4)             :: LMNX
      integer(4)             :: thistask,ntasks
!     ******************************************************************      
      call mpe$query(thistask,ntasks)
      ALLOCATE(RHO(NRL,NCOMP))
!
!     ==================================================================
!     == KINETIC ENERGY,PSEUDO DENSITY AND PROJECTIONS  ================
!     ==================================================================
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,F)
      IF(NDIM.EQ.1) THEN
        ALLOCATE(RHO1(NRL,NSPIN))
      ELSE
        ALLOCATE(RHO1(NRL,NDIM*(NDIM+1)/2))
      END IF
      EKIN=0.D0
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES$SELECTWV(IKPT,ISPIN)
          NGL=GSET%NGL        
          NBH=THIS%NBH
!         == KINETIC ENERGY ============================================
          CALL WAVES_EKIN(NGL,NDIM,NBH,NB,F(1,IKPT,ISPIN),THIS%PSI0,EKIN1)
          EKIN=EKIN+GSET%WKPT*EKIN1
!         == DENSITY ===================================================
          CALL WAVES_DENSITY(NGL,NRL,NDIM,NB,NBH,F(1,IKPT,ISPIN) &
     &              ,THIS%PSI0,RHO1,THIS%PSIR)
          IF(NDIM.EQ.1) THEN
            RHO(:,ISPIN)=RHO(:,ISPIN)+GSET%WKPT*RHO1(:,1)
          ELSE
            RHO=RHO+GSET%WKPT*RHO1
          END IF
!         == PROJECTIONS ===============================================
          CALL WAVES_PROJECTIONS(NGL,NDIM,NBH,MAP%NPRO,MAP%NBAREPRO,1,LMX,MAP &
     &          ,GSET%PRO,GSET%YLM,THIS%PSI0,THIS%PROJ)
        ENDDO
      ENDDO
      CALL ENERGYLIST$SET('PS  KINETIC',EKIN)
      CALL ENERGYLIST$ADD('AE  KINETIC',EKIN)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EKIN)
!
!     ==================================================================
!     == ONE-CENTER DENSITY MATRICES                                  ==
!     ==================================================================
      NAT=MAP%NAT
      LMNXX=0
      DO ISP=1,map%NSP
        LMNXX=MAX(LMNXX,MAP%LMNX(ISP))
      ENDDO
      ALLOCATE(DENMAT(LMNXX,LMNXX,NCOMP,NAT))
      DO IAT=1,NAT
        IPRO=MAP%IPRO(IAT)
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        ALLOCATE(PROJ(LMNX,NBH,NDIM))
        PROJ(:,:,:)=THIS%PROJ(IPRO:IPRO-1+LMNX,:,:)
        ALLOCATE(DENMAT1(LMNX,LMNX,NCOMP))
        CALL WAVES_DENMAT(NDIM,NBH,NB,LMNX,F,PROJ,DENMAT1)
        DO LMN1=1,LMNX
          DO LMN2=1,LMNX
            DENMAT(LMN1,LMN2,:,IAT)=DENMAT1(LMN1,LMN2,:)
          ENDDO
        ENDDO
        DEALLOCATE(DENMAT1)
        DEALLOCATE(PROJ)
      ENDDO
!
!     ==================================================================
!     == MULTIPOLE MOMENTS OF ONE-CENTER PART                         ==
!     ==================================================================
! define lmrx
      DO IAT=1,NAT
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        ALLOCATE(DENMAT1(LMNX,LMNX,NCOMP))
        DENMAT1(:,:,:)=DENMAT(1:LMNX,1:LMNX,:,IAT)
! QLM HAS NO ATOM INDEX!
        CALL AUGMENTATION_MOMNTS(ISP,NSPIN,LMNX,DENMAT1,LMRX,QLM)
        DEALLOCATE(DENMAT1)
      ENDDO
!
!     ==================================================================
!     == POTENTIAL                                                    ==
!     ==================================================================
!TAKE CARE OF TGRA!!!
      IF(NSPIN.EQ.2) THEN
        CALL POTENTIAL$VOFRHO(NRL,NSPIN,RHO,NDIM)
      ELSE
        CALL POTENTIAL$VOFRHO(NRL,NCOMP,RHO,NDIM)
      END IF
!
!     ==================================================================
!     == AUGMENTATION                                                 ==
!     ================================================================== 
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES$SELECTWV(IKPT,ISPIN)
          NBH=THIS%NBH
          ALLOCATE(THIS%HPROJ(map%NPRO,NBH,NDIM))
        ENDDO
      ENDDO
!define rhob
      DO IAT=1,NAT,NTASKS
        ISP=MAP%ISP(IAT)
        LMNX=MAP%LMNX(ISP)
        CALL AUGMENTATION$SPHERE(ISP,IAT,LMNX,Ncomp,DENMAT,LMRX,VQLM,RHOB &
     &                              ,DATH,DO)
        CALL WAVES_HPRO(THIS)
        DO IKPT=1,NKPT
          DO ISPIN=1,NSPIN
            CALL WAVES$SELECTWV(IKPT,ISPIN)
            NBH=THIS%NBH
            ALLOCATE(THIS%HPROJ(map%NPRO,NBH,NDIM))
            DO IBH=1,NBH
              DO LMN1=1,LMNX
                CSUM=(0.D0,0D0)
                DO LMN2=1,LMNX
                  CSUM=CSUM+DATH(LMN1,LMN2,IDIM)*THIS%PROJ(LMN2,IB,IDIM)
                ENDDO
                THIS%HPROJ(IPRO,IBH,LMN1)=CSUM
              ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO
!
!     CALL MPE$COMBINE('+',ONL)
!     CALL MPE$COMBINE('+',HNL)
!     ==================================================================
!     == HPRO                                                         ==
!     ==================================================================
!
!     ==================================================================
!     == FORCES                                                       ==
!     ==================================================================
       DO IKPT=1,NKPT
         DO ISPIN=1,NSPIN
           CALL WAVES$SELECTWV(IKPT,ISPIN)
           NGL=GSET%NGL        
           NBH=THIS%NBH
           ALLOCATE(DPROJ(NDIM,NBH,3,map%NPRO))
           CALL WAVES_PROJECTIONS(NGL,NDIM,NBH,map%NPRO,MAP%NBAREPRO,3,LMX,MAP &
      &          ,GSET%PRO,GSET%YLM,THIS%PSI0,DPROJ)
 !CALCULATE ONE-CENTER FORCES
           DEALLOCATE(DPROJ)          
         ENDDO
       ENDDO
!
!     ==================================================================
!     == HPSI                                                         ==
!     ==================================================================
      
      RETURN
      END
!
!     .............................................WAVES_PSEKIN..........
      SUBROUTINE WAVES_DENMAT(NDIM,NBH,NB,LMNX,F,PSIPRO,DENMAT)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATE 1C DENSITY MATRIX                                  **
!     **                                                              **
!     ** SUPERWAVE FUNCTIONS: PSI=PSI1+I*PSI2                         **
!     ** (PSI1-I*PSI2)*(PSI1+I*PSI2) = PSI1**1+PSI2**2                **
!     ** (PSI1+I*PSI2)*(PSI1+I*PSI2) = PSI1**1-PSI2**2+2I*PSI1PSI2    **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NDIM   ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN) :: NBH    ! #(WAVE FUNCTIONS
      INTEGER(4),INTENT(IN) :: NB     ! #(STATES)
      INTEGER(4),INTENT(IN) :: LMNX   ! #(PROJECTORS ON THIS SITE)
      REAL(8)   ,INTENT(IN) :: F(NB)  ! OCCUPATIONS
      COMPLEX(8),INTENT(IN) :: PSIPRO(NDIM,NBH,LMNX)
      REAL(8)   ,INTENT(OUT):: DENMAT(LMNX,LMNX,NDIM*(NDIM+1)/2)
      LOGICAL(4)            :: TINV
      INTEGER(4)            :: LMN1,LMN2,IDIM1,IDIM2,IDIM3
      INTEGER(4)            :: IB
      REAL(8)               :: SUM,SVAR1,SVAR2
!     ******************************************************************
      TINV=NBH.NE.NB
      DO LMN1=1,LMNX
        DO LMN2=1,LMNX
          IDIM3=0
          DO IDIM1=1,NDIM
            DO IDIM2=IDIM1,1,-1
              IDIM3=IDIM3+1
              SUM=0.D0
              DO IB=1,NBH
                SVAR1=REAL(CONJG(PSIPRO(IDIM1,IB,LMN1))*PSIPRO(IDIM2,IB,LMN2))
                IF(.NOT.TINV) THEN
                  SUM=SUM+SVAR1*F(IB)
                ELSE
!                 == special for super wave functions =================
                  SVAR2=REAL(PSIPRO(IDIM1,IB,LMN1)*PSIPRO(IDIM2,IB,LMN2))
                  SUM=SUM+0.5D0*(F(2*IB-1)*(SVAR1+SVAR2) &
     &                          +F(2*IB)  *(SVAR1-SVAR2))
                ENDIF
              ENDDO
              DENMAT(LMN1,LMN2,IDIM3)=SUM
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     .............................................WAVES_PSEKIN..........
      SUBROUTINE WAVES_EKIN(NGL,NDIM,NBH,NB,F,PSI,EKIN)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATE PS KINETIC ENERGY IN G-SPACE                       **
!     **  EVALUATE NUMBER OF ELECTRONS IN G-SPACE                     **
!     **                                                              **
!     **                                                              **
!     **  remarks:                                                    **
!     **    requires planewave object to be set properly              **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1999)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB         ! #(states)
      INTEGER(4),INTENT(IN) :: NBH        ! #(wave functions)
      INTEGER(4),INTENT(IN) :: Ndim       ! #(spinor components)
      INTEGER(4),INTENT(IN) :: NGL        ! #(PLANE WAVES)
      REAL(8)   ,INTENT(IN) :: F(NB)      ! OCCUPATION
      COMPLEX(8),INTENT(IN) :: PSI(NGL,NDIM,NBH) ! PS-WAVE FUNCTION
      REAL(8)   ,INTENT(OUT):: EKIN       ! KINETIC ENERGY
      REAL(8)               :: G2(NGL)    ! G**2
      COMPLEX(8)            :: PSI1(NGL)
      COMPLEX(8)            :: PSI2(NGL)
      LOGICAL(4)            :: TINV
      COMPLEX(8)            :: CSVAR1,CSVAR2
      INTEGER(4)            :: IB,IG,iDIM
!     ******************************************************************
      CALL PLANEWAVE$GETRL4('TINV',TINV)
      CALL PLANEWAVE$GETR8A('G2',NGL,G2)
      EKIN=0.D0
      DO IB=1,NBH
        DO IDIM=1,NDIM
          DO IG=1,NGL
            PSI1(IG)=G2(IG)*PSI(NGL,IDIM,IB)
          ENDDO
          CALL PLANEWAVE$SCALARPRODUCT(NGL,1,PSI1,1,PSI(1,IDIM,IB) &
     &                                ,CSVAR1,.FALSE.)
          IF(.NOT.TINV) THEN
            EKIN=EKIN+REAL(CSVAR1,8)*F(IB)
          ELSE
            CALL PLANEWAVE$INVERTG(NGL,PSI1,PSI2)
            CALL PLANEWAVE$SCALARPRODUCT(NGL,1,PSI2,1,PSI(1,idim,IB),CSVAR2,.FALSE.)
            EKIN=0.5D0*(F(2*IB-1)*REAL(CSVAR1+CSVAR2) &
     &                 +F(2*IB)  *REAL(CSVAR1-CSVAR2))   
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END
!
!     .....................................................RHOOFR.......
      SUBROUTINE WAVES_DENSITY(NGL,NRL,NDIM,NB,NBH,F,PSIOFG,RHO,PSIOFR)
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
      COMPLEX(8),INTENT(OUT) :: PSIOFR(NRL,NDIM,NBH)
      REAL(8)   ,INTENT(OUT) :: RHO(NRL,NDIM*(NDIM+1)/2) ! DENSITY IN R-SPACE
      COMPLEX(8),ALLOCATABLE :: EIKR(:)     !(nrl) bloch phase factor e**(ikr)
      REAL(8)   ,ALLOCATABLE :: PSI1(:,:)   
      REAL(8)   ,ALLOCATABLE :: PSI2(:,:)
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IB,IBH,IDIM1,IDIM2,ir,ic
      REAL(8)                :: FAC1,FAC2
      complex(8)             :: csvar
!     ******************************************************************
!
!     ==================================================================
!     ==  FFT                                                         ==
!     ==================================================================
      CALL PLANEWAVE$FFT('GTOR',NBH,NGL,PSIOFG,NRL,PSIOFR)
!
!     ================================================================
!     ==  CALCULATE CHARGE DENSITY                                  ==
!     ================================================================
      CALL PLANEWAVE$GETL4('TINV',TINV)
      IF(.NOT.TINV) THEN
        DO IB=1,NB
          FAC1=F(IB)
          ic=0
          DO IDIM1=1,NDIM
            DO IDIM2=IDIM1,NDIM
              ic=ic+1
              DO IR=1,NRL
                RHO(IR,ic)=RHO(IR,ic) &
     &                    +FAC1*PSIOFR(IR,IDIM1,IB)*PSIOFR(IR,IDIM2,IB)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        ALLOCATE(EIKR(NRL))
        ALLOCATE(PSI1(NRL,NDIM))
        ALLOCATE(PSI2(NRL,NDIM))
        CALL PLANEWAVE$GETC8('EIKR',NGL,EIKR)
!     ==  THE WAVE FUNCTION IS THE REAL AND IMAGINARY PART OF PSI, NOT U
        DO IBH=1,NBH
          FAC1=F(2*IB-1)
          FAC2=F(2*IB)
          DO IDIM1=1,NDIM
            DO IR=1,NRL
              CSVAR=PSIOFR(IR,IDIM1,IB)*EIKR(IR)
              PSI1(IR,IDIM1)=REAL(CSVAR)
              PSI2(IR,IDIM1)=AIMAG(CSVAR)
            ENDDO
          ENDDO
          Ic=0
          DO IDIM1=1,NDIM
            DO IDIM2=NDIM,IDIM1,-1
              Ic=Ic+1
              DO IR=1,NRL
                RHO(IR,Ic)=RHO(IR,Ic) &
     &                    +FAC1*PSI1(IR,IDIM1)*PSI1(IR,IDIM2) &
     &                    +FAC2*PSI2(IR,IDIM1)*PSI2(IR,IDIM2)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(EIKR)
        DEALLOCATE(PSI1)
        DEALLOCATE(PSI2)
      ENDIF          
      RETURN
      END
!
!     .....................................................NLSMX........
      SUBROUTINE WAVES_ADDPRO(NGL,NDIM,NB,NPRO,NBAREPRO,LMX,MAP &
     &          ,BAREPRO,YLM,PSI,PSIPRO)
!     ******************************************************************
!     **                                                              **
!     **  adds projectors to wave functions                           **
!     **      |Psi>=|psi>+ sum_i |p_i> c_i                            **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1998)***
      use waves_module, only : map_type
      IMPLICIT NONE
      INTEGER(4)    ,INTENT(IN) :: NGL      ! #(PLANE WAVES)
      INTEGER(4)    ,INTENT(IN) :: NDIM     ! #(SPINOR DIMENSIONS)
      INTEGER(4)    ,INTENT(IN) :: NB       ! #(STATES)
      INTEGER(4)    ,INTENT(IN) :: NPRO     ! #(PROJECTORS)
      INTEGER(4)    ,INTENT(IN) :: NBAREPRO ! #(PROJECTORS)
      INTEGER(4)    ,INTENT(IN) :: LMX      ! #(ANGULAR MOMENTA ON YLM)
      TYPE(MAP_TYPE),INTENT(IN) :: MAP
      REAL(8)       ,INTENT(IN) :: BAREPRO(NGL,NBAREPRO)
      REAL(8)       ,INTENT(IN) :: YLM(NGL,LMX)   ! ANGULAR MOMENTA
      COMPLEX(8)    ,INTENT(INOUT) :: PSI(NGL,NDIM,NB) !WAVEFUNCTIONS
      COMPLEX(8)    ,INTENT(IN) :: PSIPRO(NDIM,NB,NPRO) !
      COMPLEX(8)    ,ALLOCATABLE:: PRO(:,:)     !(NGL,NPRO)
      COMPLEX(8)    ,ALLOCATABLE:: EIGR(:,:)      !(NGL,NAT)   STRUCTURE FACTORS
      REAL(8)       ,ALLOCATABLE:: GVEC(:,:)      !(3,NGL)
      INTEGER(4)    ,PARAMETER  :: NPROTARGET=64
      INTEGER(4)                :: IPRO1            !FIRST PROJECTOR
      INTEGER(4)                :: NPRO1            !LAST PROJECTOR CALCULATED
      INTEGER(4)                :: NPRO2
      INTEGER(4)                :: nat
      INTEGER(4)                :: idim,ib,ig,ipro2,ipro
      complex(8)                :: csvar
!     ******************************************************************
      NPRO1=int(real(NPROTARGET)/real(NDIM))+1
      ALLOCATE(PRO(NGL,NPRO1))
      NAT=MAP%NAT
      IPRO1=0
      DO IPRO=1,NPRO,NPRO1
        CALL WAVES_EXPANDPRO(1,NBAREPRO,NAT,LMX,NGL,IPRO1+1,NPRO1 &
     &                          ,MAP,BAREPRO,YLM,EIGR,GVEC,PRO)
!
!       ===============================================================
!       ==  <PSI|PRO>                                                ==
!       ===============================================================
        NPRO2=MIN(NPRO1,NPRO-IPRO1+1)
        DO IPRO2=1,NPRO2
          DO IB=1,NB
            DO IDIM=1,NDIM
              CSVAR=PSIPRO(IDIM,IB,IPRO1+IPRO2-1)
              DO IG=1,NGL
                PSI(IG,IDIM,IB)=PSI(IG,IDIM,IB)+PRO(IG,IPRO2)*CSVAR     
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(PRO)
      RETURN
      END
!
!     .....................................................NLSMX........
      SUBROUTINE WAVES_PROJECTIONS(NGL,NDIM,NB,NPRO,NBAREPRO,NmDER,LMX,MAP &
     &          ,BAREPRO,YLM,PSI,PSIPRO)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE PROJECTIONS                                       **
!     **                                                              **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1998)***
      use waves_module, only : map_type
      IMPLICIT NONE
      INTEGER(4)    ,INTENT(IN) :: NGL      ! #(PLANE WAVES)
      INTEGER(4)    ,INTENT(IN) :: NDIM     ! #(SPINOR DIMENSIONS)
      INTEGER(4)    ,INTENT(IN) :: NB       ! #(STATES)
      INTEGER(4)    ,INTENT(IN) :: NPRO     ! #(PROJECTORS)
      INTEGER(4)    ,INTENT(IN) :: NBAREPRO ! #(PROJECTORS)
      INTEGER(4)    ,INTENT(IN) :: LMX      ! #(ANGULAR MOMENTA ON YLM)
      INTEGER(4)    ,INTENT(IN) :: NMDER     !#(DERIVATIVES: 1,3,6)
      TYPE(MAP_TYPE),INTENT(IN) :: MAP
      REAL(8)       ,INTENT(IN) :: BAREPRO(NGL,NBAREPRO)
      REAL(8)       ,INTENT(IN) :: YLM(NGL,LMX)   ! ANGULAR MOMENTA
      COMPLEX(8)    ,INTENT(IN) :: PSI(NGL,NDIM,NB) !WAVEFUNCTIONS
      COMPLEX(8)    ,INTENT(OUT):: PSIPRO(NDIM,NB,NMDER,NPRO) !
      COMPLEX(8)    ,ALLOCATABLE:: PRO(:,:,:)     !(NGL,NDIM,NPRO)
      COMPLEX(8)    ,ALLOCATABLE:: EIGR(:,:)      !(NGL,NAT)   STRUCTURE FACTORS
      REAL(8)       ,ALLOCATABLE:: GVEC(:,:)      !(3,NGL)
      INTEGER(4)    ,PARAMETER  :: NPROTARGET=64
      INTEGER(4)                :: IPRO1            !FIRST PROJECTOR
      INTEGER(4)                :: NPRO1            !LAST PROJECTOR CALCULATED
      INTEGER(4)                :: NPRO2
      INTEGER(4)                :: ipro
!     ******************************************************************
      NPRO1=NPROTARGET/(NDIM*NMDER)+1
      ALLOCATE(PRO(NGL,NDIM,NPRO1))
      IPRO1=0
      DO IPRO=1,NPRO,NPRO1
        CALL WAVES_EXPANDPRO(NMDER,NBAREPRO,map%NAT,LMX,NGL,IPRO1+1,NPRO1 &
     &                          ,MAP,BAREPRO,YLM,EIGR,GVEC,PRO)
!
!       ===============================================================
!       ==  <PSI|PRO>                                                ==
!       ===============================================================
        NPRO2=MIN(NPRO1,NPRO-IPRO1+1)
        CALL PLANEWAVE$SCALARPRODUCT(NGL,NDIM*NB,PSI,NMDER*NPRO2,PRO &
     &         ,PSIPRO(1,1,1,IPRO1),.FALSE.)
      ENDDO
      DEALLOCATE(PRO)
      RETURN
      END
!
!     .....................................................NLSMX........
      SUBROUTINE WAVES_EXPANDPRO(NMDER,npro,NBAREPRO,NAT,LMX,NGL,IPRO1,NPRO1 &
     &                          ,MAP,BAREPRO,YLM,EIGR,GVEC,PRO)
!     ******************************************************************
!     **                                                              **
!     **  TAKES THE BARE PROJECTOR FUNCTIONS AND CALCULATES FULL      **
!     **  PROJECTOR FUNCTIONS WITH STRUCTURE FACTOR AND I**L          **
!     **                                                              **
!     **  FOR NDIM.NE.3 FIRST DERIVATIVES ARE CALCULATED              **
!     **  FOR NDIM.NE.6 SECOND DERIVATIVES ARE CALCULATED             **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1998)***
      use waves_module, only : map_type
      IMPLICIT NONE
      INTEGER(4)    ,INTENT(IN) :: NMDER    ! #(DERIVATIVES 1,3,6)
      INTEGER(4)    ,INTENT(IN) :: NAT      ! #(ATOMS)
      INTEGER(4)    ,INTENT(IN) :: NGL      !#(G-VECTORS(LOCAL))
      INTEGER(4)    ,INTENT(IN) :: LMX      !X#(ANGULAR MOMENTA))
      INTEGER(4)    ,INTENT(IN) :: IPRO1    ! FIRST PROJECTOR
      INTEGER(4)    ,INTENT(IN) :: NPRO1    !X#(NUMBER OF PROJECTORS)
      INTEGER(4)    ,INTENT(IN) :: Nbarepro !X#(NUMBER OF PROJECTORS)
      INTEGER(4)    ,INTENT(IN) :: Npro     !X#(NUMBER OF PROJECTORS)
      TYPE(MAP_TYPE),INTENT(IN) :: MAP
      REAL(8)       ,INTENT(IN) :: BAREPRO(NGL,NBAREPRO)
      REAL(8)       ,INTENT(IN) :: YLM(NGL,LMX) !CI**(-L)*SPHERICAL HARMONICS
      COMPLEX(8)    ,INTENT(IN) :: EIGR(NGL,NAT)      ! STRUCTURE FACTORS
      REAL(8)       ,INTENT(IN) :: GVEC(3,NGL)
      COMPLEX(8)    ,INTENT(OUT):: PRO(NGL,NMDER,NPRO)
      COMPLEX(8)    ,PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)                :: CSVAR
      INTEGER(4)                :: NDER
      INTEGER(4)                :: ICOUNT,IPRO
      INTEGER(4)                :: IAT,ISP,LMN,LN,IM,L,IG,lm,ibarepro
!     ******************************************************************
      IF(NMDER.EQ.1) THEN
        NDER=0
      ELSE IF(NMDER.EQ.3) THEN
        NDER=1
      ELSE IF(NMDER.EQ.6) THEN
        NDER=2
      ELSE
        CALL ERROR$MSG('NMDER MUST BE 1,3, OR 6')
        CALL ERROR$STOP('WAVES_PROJECTIONS')
      END IF
!
!     ==================================================================
!     ==  START THE LOOP                                              ==
!     ==================================================================
      ICOUNT=0
      IPRO=0
      DO IAT=1,MAP%NAT
        ISP=MAP%ISP(IAT)
        IBAREPRO=MAP%IBAREPRO(IAT)-1
        LMN=0
        DO LN=1,MAP%LNX(ISP)
          IBAREPRO=IBAREPRO+1
          L=MAP%LOX(LN,ISP)
          LM=L**2
          DO IM=1,2*L+1
            LM=LM+1
            IPRO=IPRO+1
            IF(IPRO.LT.IPRO1) CYCLE
            ICOUNT=ICOUNT+1
            IF(ICOUNT.GT.NPRO1) EXIT
!           
!           ==========================================================
!           ==  COMPOSE PROJECTOR FUNCTIONS                         ==
!           ==========================================================
!           ==  MULTIPLY PROJECTOR WITH STRUCTURE FACTOR ===========
            CSVAR=(-CI)**L
            DO IG=1,NGL
              PRO(IG,1,IPRO)=YLM(IG,LM)*BAREPRO(IG,IBAREPRO) &
           &                           *(CSVAR*EIGR(IG,IAT))
            ENDDO
!           ==  MULTIPLY WITH GRADIENT =============================
            IF(NDER.EQ.1) THEN
              DO IG=1,NGl
                PRO(IG,3,IPRO)=-GVEC(3,IG)*CI*PRO(IG,1,IPRO)
                PRO(IG,2,IPRO)=-GVEC(2,IG)*CI*PRO(IG,1,IPRO)
                PRO(IG,1,IPRO)=-GVEC(1,IG)*CI*PRO(IG,1,IPRO)
              ENDDO
            ELSE IF(NDER.EQ.2) THEN
              DO IG=1,NGl
                PRO(IG,6,IPRO)=-GVEC(3,IG)*GVEC(3,IG)*PRO(IG,1,IPRO)
                PRO(IG,5,IPRO)=-GVEC(2,IG)*GVEC(3,IG)*PRO(IG,1,IPRO)
                PRO(IG,4,IPRO)=-GVEC(2,IG)*GVEC(2,IG)*PRO(IG,1,IPRO)
                PRO(IG,3,IPRO)=-GVEC(1,IG)*GVEC(3,IG)*PRO(IG,1,IPRO)
                PRO(IG,2,IPRO)=-GVEC(1,IG)*GVEC(2,IG)*PRO(IG,1,IPRO)
                PRO(IG,1,IPRO)=-GVEC(1,IG)*GVEC(1,IG)*PRO(IG,1,IPRO)
              ENDDO
            END IF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
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
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)   ,PARAMETER  :: MBYTE=2.D0**20
      INTEGER(4)            :: IKPT
      REAL(8)               :: RY
      REAL(8)               :: SVAR
      REAL(8)               :: MEMORY
!     ******************************************************************
      CALL CONSTANTS('RY',RY)
      CALL REPORT$TITLE(NFIL,'WAVE FUNCTIONS')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF BANDS',NB,' ')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF K-POINTS',NKPT,' ')
      CALL REPORT$I4VAL(NFIL,'NUMBER OF SPINS',NSPIN,' ')
      CALL REPORT$R8VAL(NFIL,'PLANE WAVE CUTOFF',EPWPSI/RY,'RY')
      CALL REPORT$R8VAL(NFIL,'WAVE FUNCTION MASS',EMASS,'A.U.')
      CALL REPORT$R8VAL(NFIL,'G**2 ENHANCEMENT OF FUNCTION MASS',EMASSCG2,' ')
      IF(ANNEE.NE.0) THEN
        CALL REPORT$R8VAL(NFIL,'FRICTION',ANNEE,' ')
      END IF
      IF(TSTOP) THEN
        CALL REPORT$CHVAL(NFIL,'INITIAL VELOCITY IS SET TO','ZERO')
      END IF
      IF(TRANDOMIZE) THEN
        CALL REPORT$R8VAL(NFIL &
     &      ,'INITIAL VELOCITIES ARE RANDOMIZED WITH ENERGY',AMPRE,'H')
      END IF
      IF(.NOT.TSAFEORTHO) THEN
        WRITE(NFIL,FMT='("EIGENSTATES ARE CALCULATED." &
     &                  ," (NO STRICT ENERGY CONSERVATION)")')
      END IF
!     
!     ================================================================
!     ==  REPORT INFORMATION ABOUT G-VECTORS                        ==
!     ================================================================
      DO IKPT=1,NKPT
        CALL WAVES_SELECTWV(IKPT,1)
        CALL REPORT$I4VAL(NFIL &
     &       ,'NUMBER OF PLANE WAVES FOR WAVE FUNCTION',GSET%NGL,' ')
      ENDDO
      CALL REPORT$I4VAL(NFIL,'#R-POINTS IN 1. DIRECTION',NR1,' ')
      CALL REPORT$I4VAL(NFIL,'#R-POINTS IN 2. DIRECTION',NR2,' ')
      CALL REPORT$I4VAL(NFIL,'#R-POINTS IN 3. DIRECTION',NR3,' ')
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$GETR8(ID_,VAL_)
!     ******************************************************************
!     **  WAVES_GETR8                                                 **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(OUT):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'EPWPSI') THEN
        VAL_=EPWPSI
      ELSE IF(ID_.EQ.'EMASS') THEN
        VAL_=EMASS
      ELSE IF(ID_.EQ.'EMASSCG2') THEN
        VAL_=EMASSCG2
      ELSE IF(ID_.EQ.'FRICTION') THEN
        VAL_=ANNEE
      ELSE
        CALL ERROR$MSG('IDENT_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('IDENT_',ID_)
        CALL ERROR$STOP('WAVES$GETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$GETR8A(ID_,LEN_,VAL_)
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LEN_
      REAL(8)     ,INTENT(OUT):: VAL_(LEN_)
      REAL(8)     ,ALLOCATABLE:: DWORK(:)
      REAL(8)     ,ALLOCATABLE:: DWORK2(:,:)
      COMPLEX(8)  ,ALLOCATABLE:: CWORK(:)
      INTEGER(4)              :: IB,IKPT,ISPIN  ! STATE SPECIFIERS
      INTEGER(4)              :: IAT
      INTEGER(4)              :: IDUMMY
      INTEGER(4)              :: NNRL
      INTEGER(4)              :: LMN            ! DO LOOP COUNTER
      REAL(8)                 :: SVAR1,SVAR2
      INTEGER(4)              :: I
!     ******************************************************************
                                         CALL TRACE$PUSH('WAVES$GETR8A')
      CALL WAVES_NEWLIST
      CALL LINKEDLIST$GET(LL_WAVE,'IB',1,IB)
      CALL LINKEDLIST$GET(LL_WAVE,'IKPT',1,IKPT)
      CALL LINKEDLIST$GET(LL_WAVE,'ISPIN',1,ISPIN)
      CALL LINKEDLIST$GET(LL_WAVE,'IAT',1,IAT)
!     == CHECK DIMENSIONS
      IF(IB.LE.0.OR.IB.GT.NB) THEN
        CALL ERROR$MSG('BAND DIMENSION INCONSISTENT')
        CALL ERROR$I4VAL('IB ',IB)
        CALL ERROR$I4VAL('NB ',NB)
        CALL ERROR$STOP('WAVES$GETR8A)')
      END IF
      IF(IKPT.LE.0.OR.IKPT.GT.NKPT) THEN
        CALL ERROR$MSG('K-POINT DIMENSION  INCONSISTENT')
        CALL ERROR$I4VAL('IKPT ',IKPT)
        CALL ERROR$I4VAL('NKPT ',NKPT)
        CALL ERROR$STOP('WAVES$GETR8A')
      END IF
      IF(ISPIN.LE.0.OR.ISPIN.GT.NSPIN) THEN
        CALL ERROR$MSG('SPIN-DIMENSION  INCONSISTENT')
        CALL ERROR$I4VAL('ISPIN ',ISPIN)
        CALL ERROR$I4VAL('NSPIN ',NSPIN)
        CALL ERROR$STOP('WAVES$GETR8A')
      END IF
!
!     =================================================================
!     ==  HAMILTON MATRIX                                            ==
!     =================================================================
      IF(ID_.EQ.'<PSI|H|PSI>') THEN
        IF(LEN_.NE.NB*NB) THEN
          CALL ERROR$MSG('LENGTH INCONSISTENT')
          CALL ERROR$CHVAL('ID_',ID_)
          CALL ERROR$I4VAL('LEN_',LEN_)
          CALL ERROR$I4VAL('NB*NB',NB*NB)
          CALL ERROR$STOP('WAVES$GETR8A')
        END IF        
           
        IF(ALLOCATED(HAMU)) THEN
          ALLOCATE(DWORK2(NB,NB))
          DO I=1,NB
            DWORK2(:,I)=HAMU(:,I,IKPT,ISPIN)*HAME(I,IKPT,ISPIN)
          ENDDO        
          DWORK2=MATMUL(DWORK2,TRANSPOSE(HAMU(:,:,IKPT,ISPIN)))
          VAL_=RESHAPE(DWORK2,(/NB*NB/))
          DEALLOCATE(DWORK2)
        ELSE
          CALL ERROR$MSG('HAMILTON MATRIX NOT AVAILABLE')
          CALL ERROR$STOP('WAVES$GETR8A(<PSI|H|PSI>)')
        END IF
!     
!     ================================================================
!     ==  REAL SPACE PS WAVE FUNCTION                               ==
!     ==  NOTE: IB,IKPT,ISPIN MUST BE SET ON LINKEDLIST             ==
!     ================================================================
      ELSE IF(ID_.EQ.'PSPSI') THEN
!     
!       == EVALUATE DENSITY AND REAL SPACE WAVE FUNCTIONS ============== 
        CALL PLANEWAVE$SELECT('WAVE',IKPT)
        CALL PLANEWAVE$RSIZE(IDUMMY,NR1L,IDUMMY)
        NNRL=NR1L*NR2*NR3
        IF(LEN_.NE.NNRL) THEN
          CALL ERROR$MSG('R-DIMENSIONS OF DENSITY INCONSISTENT')
          CALL ERROR$I4VAL('LEN_',LEN_)
          CALL ERROR$I4VAL('NNRL_',NNRL)
          CALL ERROR$STOP('WAVES$GET(PSPSI)')
        END IF
!
!       ================================================================
!       == RETURN RAW WAVE FUNCTIONS                                  ==
!       ================================================================
        IF(TRAWSTATES) THEN
!PRINT*,'RAWSTATES PSWAVE'
          ALLOCATE(DWORK2(NNRL,2))
          ALLOCATE(CWORK(NGWLX))
          VAL_(:)=0.D0
          CALL PLANEWAVE$SUPFFT('GTOR',' ',1,NNRL,NR1L,NR2,NR3,NGWLX &
      &     ,DWORK2(1,1),DWORK2(1,2),C0(:,IB,IKPT,ISPIN),CWORK)
          VAL_(:)=VAL_(:)+DWORK2(:,1)
          DEALLOCATE(DWORK2)
          DEALLOCATE(CWORK)
                          CALL TRACE$POP
          RETURN
        ELSE
!
!       ================================================================
!       == DIAGONALIZE HAMILTONIAN                                    ==
!       ================================================================
!PRINT*,'EIGSTATES PSWAVE'
          IF(.NOT.ALLOCATED(HAMU)) THEN
            CALL ERROR$MSG('TRANSFORMATION TO EIGENSTATES NOT AVAILABLE')
            CALL ERROR$STOP('WAVES$GETR8A')
          END IF
          ALLOCATE(DWORK2(NNRL,2))
          VAL_(:)=0.D0
          DO I=1,NB,2
            IF(I.NE.NB) THEN
              CALL PLANEWAVE$SUPFFT('GTOR',' ',2,NNRL,NR1L,NR2,NR3,NGWLX &
      &       ,DWORK2(1,1),DWORK2(1,2),C0(:,I,IKPT,ISPIN),C0(:,I+1,IKPT,ISPIN))
              VAL_(:)=VAL_(:)+DWORK2(:,1)*HAMU(I,IB,IKPT,ISPIN) &
      &                      +DWORK2(:,2)*HAMU(I+1,IB,IKPT,ISPIN)
            ELSE
              CALL PLANEWAVE$SUPFFT('GTOR',' ',1,NNRL,NR1L,NR2,NR3,NGWLX &
      &       ,DWORK2(1,1),DWORK2(1,2),C0(:,I,IKPT,ISPIN),C0(:,I,IKPT,ISPIN))
              VAL_(:)=VAL_(:)+DWORK2(:,1)*HAMU(I,IB,IKPT,ISPIN)
            ENDIF
          ENDDO
          DEALLOCATE(DWORK2)
!         CALL PLANEWAVE$SUPFFT('GTOR',' ',1,NNRL,NR1L,NR2,NR3,NGWLX &
!    &                       ,VAL_(:),SVAR1,C0(:,IB,IKPT,ISPIN),SVAR2)
       END IF
!     
!     ================================================================
!     ==  GET PROJECTIONS                                           ==
!     ================================================================
      ELSE IF(ID_.EQ.'<PSPSI|PRO>') THEN
        CALL SETUP$LMNXX(LMNXX)
        IF(LEN_.NE.LMNXX) THEN
          CALL ERROR$MSG('LMN-DIMENSIONS INCONSISTENT')
          CALL ERROR$I4VAL('LEN_',LEN_)
          CALL ERROR$I4VAL('LMNXX',LMNXX)
          CALL ERROR$STOP('WAVES$GET(<PSPSI|PRO>)')
        END IF
!
        IF(TRAWSTATES) THEN
!PRINT*,'RAWSTATES <PSI|PRO>'
          VAL_(:)=FNL(IAT,IB,:,IKPT,ISPIN)
        ELSE
!PRINT*,'EIGSTATES <PSI|PRO>'
          ALLOCATE(DWORK(LMNXX))
          DWORK=0.D0
          DO LMN=1,LMNXX
            DO I=1,NB
              DWORK(LMN)=DWORK(LMN) &
     &                  +FNL(IAT,I,LMN,IKPT,ISPIN)*HAMU(I,IB,IKPT,ISPIN)
            ENDDO
          ENDDO        
          VAL_(:)=DWORK(:)
          DEALLOCATE(DWORK)
        END IF
!
!     =================================================================
!     ==  UNKNOWN SPECIFIER                                          ==
!     =================================================================
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('WAVES$GETR8A')
      END IF
           CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$GETL4(ID_,VAL_)
!     ******************************************************************
!     **  WAVES_GETL4                                                 **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID_
      LOGICAL(4)  ,INTENT(OUT) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'SAFEORTHO') THEN
        VAL_=TSAFEORTHO
      ELSE IF(ID_.EQ.'STOP') THEN
        VAL_=TSTOP
      ELSE IF(ID_.EQ.'RANDOMIZE') THEN
        VAL_=TRANDOMIZE
      ELSE IF(ID_.EQ.'STOREPSIR') THEN
        VAL_=TSTORE
      ELSE IF(ID_.EQ.'RAWSTATES') THEN
        VAL_=TRAWSTATES
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('WAVES$GETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$GETI4(ID_,VAL_)
!     ******************************************************************
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(OUT):: VAL_
!     ******************************************************************
      CALL WAVES_NEWLIST
      IF(ID_.EQ.'NSPIN') THEN
        VAL_=NSPIN
      ELSE IF(ID_.EQ.'NKPT') THEN
        VAL_=NKPT
      ELSE IF(ID_.EQ.'NB') THEN
        VAL_=NB
      ELSE IF(ID_.EQ.'NR1') THEN
        VAL_=NR1
      ELSE IF(ID_.EQ.'NR1L') THEN
        VAL_=NR1L
      ELSE IF(ID_.EQ.'NR1START') THEN
        CALL LINKEDLIST$GET(LL_WAVE,'NR1START',1,VAL_)
      ELSE IF(ID_.EQ.'NR2') THEN
        VAL_=NR2
      ELSE IF(ID_.EQ.'NR3') THEN
        VAL_=NR3
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('WAVES$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SETL4(ID_,VAL_)
!     ******************************************************************
!     **  WAVES_SET                                                   **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      LOGICAL(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'SAFEORTHO') THEN
        TSAFEORTHO=VAL_
      ELSE IF(ID_.EQ.'STOP') THEN
        TSTOP=VAL_
      ELSE IF(ID_.EQ.'RANDOMIZE') THEN
        TRANDOMIZE=VAL_
      ELSE IF(ID_.EQ.'STOREPSIR') THEN
        TSTORE=VAL_
      ELSE IF(ID_.EQ.'HAMILTON') THEN
        THAMILTON=VAL_
      ELSE IF(ID_.EQ.'RAWSTATES') THEN
        TRAWSTATES=VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('WAVES$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SETI4(ID_,VAL_)
!     ******************************************************************
!     **  WAVES_SET                                                   **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: VAL_
!     ******************************************************************
      IF(ID_.EQ.'IB') THEN
        CALL LINKEDLIST$SET(LL_WAVE,ID_,0,VAL_)
      ELSE IF(ID_.EQ.'IKPT') THEN
        CALL LINKEDLIST$SET(LL_WAVE,ID_,0,VAL_)
      ELSE IF(ID_.EQ.'ISPIN') THEN
        CALL LINKEDLIST$SET(LL_WAVE,ID_,0,VAL_)
      ELSE IF(ID_.EQ.'IAT') THEN
        CALL LINKEDLIST$SET(LL_WAVE,ID_,0,VAL_)
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('WAVES$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SETR8(ID_,VAL_)
!     ******************************************************************
!     **  WAVES_SETR8                                                 **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      REAL(8)     ,INTENT(IN):: VAL_
!     ******************************************************************
      IF(ID_.EQ.'EPWPSI') THEN
        EPWPSI=VAL_
      ELSE IF(ID_.EQ.'EMASS') THEN
        EMASS=VAL_
      ELSE IF(ID_.EQ.'EMASSCG2') THEN
        EMASSCG2=VAL_
      ELSE IF(ID_.EQ.'FRICTION') THEN
        ANNEE=VAL_
      ELSE IF(ID_.EQ.'AMPRE') THEN
        AMPRE=VAL_
      ELSE IF(ID_.EQ.'TIMESTEP') THEN
        DELT=VAL_
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('WAVES$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SETR8A(ID_,LEN_,VAL_)
!     ******************************************************************
!     **  WAVES_SETR8A                                                **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LEN_
      REAL(8)     ,INTENT(IN) :: VAL_(LEN_)
!     ******************************************************************
      IF(ID_.EQ.'RLAM(0)') THEN
        IF(LEN_.NE.NB*NB*NKPT*NSPIN) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$I4VAL('LEN_',LEN_)
          CALL ERROR$I4VAL('NB*NB*NKPT*NSPIN',NB*NB*NKPT*NSPIN)
          CALL ERROR$STOP('WAVES$SETR8A')
        END IF
        RLAM0=RESHAPE(VAL_,(/NB,NB,NKPT,NSPIN/))
      ELSE IF(ID_.EQ.'RLAM(-)') THEN
        IF(LEN_.NE.NB*NB*NKPT*NSPIN) THEN
          CALL ERROR$MSG('SIZE INCONSISTENT')
          CALL ERROR$I4VAL('LEN_',LEN_)
          CALL ERROR$I4VAL('NB*NB*NKPT*NSPIN',NB*NB*NKPT*NSPIN)
          CALL ERROR$STOP('WAVES$SETR8A')
        END IF
        RLAMM=RESHAPE(VAL_,(/NB,NB,NKPT,NSPIN/))
      ELSE
        CALL ERROR$MSG('ID_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID_',ID_)
        CALL ERROR$STOP('WAVES$SETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$INITIALIZERANDOM
!     ******************************************************************
!     **  WAVES$INITIALIZERANDOM                                      **
!     **  INITIALIZE RANDOM                                           **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL              :: TGAMMA
      INTEGER(4)           :: ISPIN,IKPT  ! DO LOOP COUNTERS
!     ******************************************************************
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          CALL WAVES$SELECTWV(IKPT,ISPIN)
          NGL=GSET%NGL
          NDIM=THIS%NDIM
          NBH=THIS%NBH
          ALLOCATE(G2(NGL))
          CALL PLANEWAVE$GETR8('G2',NGL,G2)
          CALL WAVES_INITIALIZERANDOM(NGL,NBH*NDIM,G2,THIS%PSI0)
          THIS%PSIM=THIS%PSI0
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$RANDOMIZEVELOCITY
!     ******************************************************************
!     **  WAVES$RANDOMIZEVELOCITY                                     **
!     **  RANDOMIZE VELOCITY OF THE WAVE FUNCTIONS                    **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL             :: TGAMMA
      REAL(8)             :: OCC(NB,NKPT,NSPIN)  
      INTEGER(4)          :: ISPIN,IKPT,IB ! DO LOOP COUNTERS
      INTEGER(4)          :: IDUMMY2
      INTEGER(4)          :: NGWL1,NGW1
      REAL(8)             :: FAC
      REAL(8)             :: AMPRE1
      REAL(8)             :: TOTEL      ! #(ELECTRONS)
      REAL(8)             :: TOTSTATE   ! #(STATES)
!     ******************************************************************
      IF(.NOT.TRANDOMIZE) RETURN
!     
!     ==================================================================
!     == SCALE UP TO ACCOUNT FOR PARTIAL OCCUPATION OF ALL STATES     ==
!     ==================================================================
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
      TOTEL=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IB=1,NB
            TOTEL=TOTEL+OCC(IB,IKPT,ISPIN)
          ENDDO
        ENDDO
      ENDDO
      TOTSTATE=DBLE(NSPIN*NKPT*NB)
!     
!     ==================================================================
!     == CALCULATE TOTAL NUMBER OF ELECTRONS                          ==
!     ==================================================================
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
!     
!         +++ RESCALE AMPRE ++++++++++++++++++++++++++++++++++++++++++
          CALL PLANEWAVE$SELECT('WAVE',IKPT)
          CALL PLANEWAVE$GSIZE(NGWL1,NGW1,IDUMMY2)
!         +++ END ++++++++++++++++++++++++++++++++++++++++++++++++++++
!         __SCALE UP TO ACCOUNT FOR PARTIAL OCCUPATION
          FAC=TOTSTATE/TOTEL
!         __SCALE UP TO ACCOUNT FOR CONSTRAINTS OF ORTHOGONALITY
          FAC=FAC*DBLE(NGW1-NB+1)/DBLE(NGW1)
!         __ FIND ENERGY FOR THIS FRACTION OF STATES
          FAC=FAC/DBLE(NSPIN*NKPT)*DBLE(NGWL1)/DBLE(NGW1)
!         __AMPRE IS THE AVERAGE KINETIC ENERGY FOR THIS FRACTION_____
          AMPRE1=AMPRE*FAC
          CALL WAVES_RANDOMIZEVELOCITY1(NGWLX,NGWL(IKPT),NB &
     &              ,AMPRE1,EMASS,EMASSCG2,CELLVOL,DELT &
     &              ,HSG(1,IKPT),CM(1,1,IKPT,ISPIN))
          TGAMMA=(DSQRT(HSG(1,IKPT)).LT.EPSILONGAMMA)
          CALL WAVES_CORRECTGAMMA(TGAMMA,NGWLX,NB,CM(1,1,IKPT,ISPIN))
        ENDDO
      ENDDO
      TRANDOMIZE=.FALSE.
      RETURN
      CONTAINS
!     .................................................................
      SUBROUTINE WAVES_RANDOMIZEVELOCITY1(NGWX,NGW,NB &
     &                ,AMPRE,EMASS,EMASSCG2,CELLVOL,DELT,G2,PSI)
!     ******************************************************************
!     **                                                              **
!     **  RANDOMIZE WAVE FUNCTIONS OF THE PREVIOUS TIME STEP          **
!     **  IN ORDER TO RANDOMIZE THE VELOCITIES                        **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1996)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NGW,NGWX ! # PLANE WAVES, MAX
      INTEGER(4),INTENT(IN)   :: NB       ! # BANDS
      REAL(8)   ,INTENT(IN)   :: AMPRE    ! TARGET FICT. KINETIC ENERGY
      REAL(8)   ,INTENT(IN)   :: EMASS    ! WAVE FUNCTION MASS
      REAL(8)   ,INTENT(IN)   :: EMASSCG2 ! WAVE FUNCTION MASS (GDEPFAC)
      REAL(8)   ,INTENT(IN)   :: CELLVOL  ! UNIT CELL VOLUME   
      REAL(8)   ,INTENT(IN)   :: DELT     ! TIME STEP
      REAL(8)   ,INTENT(IN)   :: G2(NGWX) ! G**2
      COMPLEX(8),INTENT(INOUT):: PSI(NGWX,NB) ! PS WAVE
      INTEGER(4)              :: IB,IG
      REAL(8)                 :: SVAR,FAC,RAN1,RAN2
      REAL(8)                 :: SUM,SUM1,SUM2
      INTEGER(4)              :: NFILO
!     ******************************************************************
      CALL FILEHANDLER$UNIT('PROT',NFILO)    
!     == THE FIRST FACTOR 2 COMES FROM REAL AND IMAGINARY PART
!     == THE SECOND FACTOR 2 COMES FROM ONLY COUNTING G BUT NOT -G
!     == THE THIRD FACTOR 2 THE FACT THAT <RAN^2>=0.5
      FAC=DELT*DSQRT(AMPRE/(CELLVOL*2.D0*2.D0*EMASS*DBLE(NB*NGW)))
      SUM=0.D0
      DO IB=1,NB
        DO IG=1,NGW
!         __ EVALUATE GAUSSIAN DISTRIBUTED RANDOM NUMBERS_______________
          CALL GAUSS_RANDOM_NUMBER(RAN1)
          CALL GAUSS_RANDOM_NUMBER(RAN2)
          SVAR=FAC/DSQRT(1.D0+EMASSCG2*G2(IG))
          SUM=SUM+SVAR**2*(RAN1**2+RAN2**2)
          PSI(IG,IB)=PSI(IG,IB)+SVAR*CMPLX(RAN1,RAN2,KIND=8)
        ENDDO
      ENDDO
      WRITE(NFILO,*)'SUM FROM RANDOMIZEVELOCITY ',SUM*2.D0/DELT**2*EMASS*CELLVOL*2.D0
      RETURN
      END SUBROUTINE WAVES_RANDOMIZEVELOCITY1
      END
!
!     ..................................................................
      SUBROUTINE WAVES$STOP
!     ******************************************************************
!     **  SET WAVE FUNCTION VELOCITIES TO ZERO                        **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)   :: IKPT,ISPIN   ! DO LOOP COUNTERS
!     ******************************************************************
      DO ISPIN=1,NSPIN
         DO IKPT=1,NKPT
           CALL WAVES_SELECTWV(IKPT,ISPIN)
           THIS%PSIM=THIS%PSI0
         ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$EIGR(STRING)
!     ******************************************************************
!     **  GENERATE STRUCTURE FACTORS                                  **
!     **  EXP(-I*(G+K)*R), WHERE R IS THE ATOMIC COORDINATE           **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4), PARAMETER   :: NMAX=500
      CHARACTER(*),INTENT(IN) :: STRING
      REAL(8)                 :: RBAS(3,3)
      REAL(8)                 :: GBAS(3,3)
      COMPLEX(8)              :: DCWORK(2,NMAX)  
      REAL(8)                 :: RPOS(3)
      INTEGER(4)              :: IKPT,ISPIN,IAT ! DO LOOP COUNTERS
!     ******************************************************************
      IF(STRING.NE.'R(0)'.AND.STRING.NE.'R(-)'.AND.STRING.NE.'R(+)') THEN
        CALL ERROR$MSG('STRING MUST BE EITHER R(0),R(-), OR R(+)')
        CALL ERROR$CHVAL('STRING',STRING)
        CALL ERROR$STOP('WAVES$EIGR')
      END IF
      DO IKPT=1,NKPT
        CALL PLANEWAVE$SETECTWV(IKPT,ISPIN)
        IF(.NOT.ALLOCATED(GSET%EIGR)) THEN
          ALLOCATE(GSET%EIGR(NGL,NAT))
        END IF
        NGL=GSET%NGL
        DO IAT=1,NAT
          CALL ATOMLIST$GETR8A(STRING(1:4),IAT,3,RPOS)
          CALL PLANEWAVE$STRUCTUREFACTOR(RPOS,NGL,EIGR(1,IAT))
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$EKIN(EKIN)
!     ******************************************************************
!     **  CALCULATE KINETIC ENERGY AND CHECK NUMBER OF ELECTRONS      **
!     ******************************************************************
      USE MPE_MODULE
      USE WAVES_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: EKIN          ! ELECTRONIC KINETIC ENERGY
      REAL(8)                :: OCC(NB,NKPT,NSPIN) ! OCCUPATIONS 
      REAL(8)                :: DHELP(2)
      REAL(8)                :: RSUM
      REAL(8)                :: EKIN1,RSUM1
      INTEGER(4)             :: IKPT,ISPIN ! DO LOOP COUNTERS
      INTEGER(4)             :: NGL
      REAL(8)   ,ALLOCATABLE :: G2
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$EKIN')
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
      EKIN=0.D0
      RSUM=0.D0
      DO IKPT=1,NKPT
        CALL WAVES$SELECTWV(IKPT,1)
        NGL=GSET%NGL
        CALL PLANEWAVE$SELECT(GSET%ID)
        ALLOCATE(G2(NGL))
        CALL PLANEWAVE$GETR8A('G2',NGL,G2)
        WKPT=GSET%WKPT
        DO ISPIN=1,NSPIN
          CALL WAVES$SELECTWV(IKPT,ISPIN)
          NBH=THIS%NBH
          CALL WAVES_EKIN(NGL,NBH,OCC(1,IKPT,ISPIN),WKPT,G2 &
     &                   ,CELLVOL,EKIN1,RSUM1,THIS%PSI0)
          EKIN=EKIN+EKIN1
          RSUM=RSUM+RSUM1
        ENDDO
      ENDDO
!     
!     ++  PARALLEL BEGIN +++++++++++++++++++++++++++++++++++++++++++++
      DHELP(1)=EKIN
      DHELP(2)=RSUM
      CALL MPE$COMBINE('+',DHELP)        
      EKIN = DHELP(1)
      RSUM = DHELP(2)
!     ++  PARALLEL END +++++++++++++++++++++++++++++++++++++++++++++++
                           CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$DENSITY(NNRX_,NSPIN_,RHO)
!     ******************************************************************
!     **  CALCULATE ELECTRON DENSITY                                  **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NNRX_
      INTEGER(4),INTENT(IN) :: NSPIN_
      REAL(8)   ,INTENT(OUT):: RHO(NNRX_,NSPIN_)
      REAL(8)               :: OCC(NB,NKPT,NSPIN)  
      INTEGER(4)            :: IKPT,ISPIN,IR    ! DO LOOP VARIABLES
      INTEGER(4)            :: IDUMMY
      INTEGER(4)            :: NNRL
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$DENSITY')
!
      IF(NNRX_.LT.NR1L*NR2*NR3) THEN
        CALL ERROR$MSG('R-DIMENSIONS OF DENSITY INCONSISTENT')
        CALL ERROR$I4VAL('NNRX_',NNRX_)
        CALL ERROR$I4VAL('NR1L*NR2*NR3',NR1L*NR2*NR3)
        CALL ERROR$STOP('WAVES$DENSITY')
      END IF
      IF(NSPIN_.LT.NSPIN) THEN
        CALL ERROR$MSG('SPIN-DIMENSION OF DENSITY INCONSISTENT')
        CALL ERROR$STOP('WAVES$DENSITY')
      END IF
!     
!     == ALLOCATE STORE ==============================================        
      IF(ALLOCATED(STORE)) THEN
        DEALLOCATE(STORE)
      END IF
      ALLOCATE(STORE(NRSTORE,NB,NKPT,NSPIN)) 
!     
!     == EVALUATE DENSITY AND REAL SPACE WAVE FUNCTIONS ============== 
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
      DO ISPIN=1,NSPIN
        DO IR=1,NNRX_
          RHO(IR,ISPIN)=0.D0
        ENDDO
        DO IKPT=1,NKPT
          CALL PLANEWAVE$SELECT('WAVE',IKPT)
          CALL PLANEWAVE$RSIZE(IDUMMY,NR1L,IDUMMY)
          NNRL=NR1L*NR2*NR3
          CALL WAVES_DENSITY(NGWLX,NRSTORE,NR1L,NR2,NR3,NB &
     &         ,OCC(1,IKPT,ISPIN),WKPT &
     &         ,C0(1,1,IKPT,ISPIN),RHO(1,ISPIN),STORE(1,1,IKPT,ISPIN))
        ENDDO
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$1CDENSITYMATRIX(IAT_,LMNXX_,NSPINX_,DENMAT)
!     ******************************************************************
!     **  CALCULATE ELECTRON DENSITY AND CHECK NUMBER OF ELECTRONS    **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT_      ! ATOM INDEX
      INTEGER(4),INTENT(IN) :: LMNXX_    ! 
      INTEGER(4),INTENT(IN) :: NSPINX_
      REAL(8)   ,INTENT(OUT):: DENMAT(LMNXX_,LMNXX_,NSPINX_)
      REAL(8)               :: OCC(NB,NKPT,NSPIN)  
      INTEGER(4),ALLOCATABLE:: LMNX(:)             !(NSP) 
      INTEGER(4)            :: ISP
      INTEGER(4)            :: IKPT,ISPIN !DO LOOP COUNTERS
      CHARACTER*1           :: OP
!     ******************************************************************
      CALL ATOMLIST$GETI4('ISPECIES',IAT_,ISP)
      CALL SETUP$NSPECIES(NSP)
      ALLOCATE(LMNX(NSP)) 
      CALL SETUP$LMNX(ISP,LMNX(ISP))
!     
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
      DO ISPIN=1,NSPIN
        OP=' '
        DO IKPT=1,NKPT
          CALL WAVES_1CDENSITYMATRIX(NAT,LMNXX_,LMNX(ISP),NB,IAT_ &
     &        ,WKPT,OCC(1,IKPT,ISPIN),FNL(1,1,1,IKPT,ISPIN) &
     &        ,DENMAT(1,1,ISPIN),OP)
          OP='+'
        ENDDO
      ENDDO
      DEALLOCATE(LMNX) 
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SET1CHAMILTON(IAT_,LMNXX_,NSPINX_,DENMAT)
!     ******************************************************************
!     **  GET (<AEPHI|AE-H1|AEPHI>-<PSPHI|PS-H1|PSPHI>)<PS-P|PS-PSI>  **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT_
      INTEGER(4),INTENT(IN) :: LMNXX_
      INTEGER(4),INTENT(IN) :: NSPINX_
      REAL(8)   ,INTENT(IN) :: DENMAT(LMNXX_,LMNXX_,NSPINX_)
      INTEGER(4)            :: ISP
      INTEGER(4)            :: LMNX1
      INTEGER(4)            :: ISPIN,IKPT
!     ******************************************************************
      CALL ATOMLIST$GETI4('ISPECIES',IAT_,ISP)
      CALL SETUP$LMNX(ISP,LMNX1)
      IF(LMNX1.GT.LMNXX_) THEN
        CALL ERROR$STOP('WAVES$1CHAMILTONIAN')
      END IF
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          CALL WAVES_1CHAMILTONIAN(NAT,LMNXX_,LMNX1,NB,IAT_ &
     &              ,FNL(1,1,1,IKPT,ISPIN),HNL(1,1,1,IKPT,ISPIN) &
     &              ,DENMAT(1,1,ISPIN))
        ENDDO
      ENDDO
!     
!     ================================================================
!     ==  SEND DATA TO OPTIC CODE                                   ==
!     ================================================================
      CALL OPTIC3$DH(NSPINX_,LMNXX_,NAT,IAT_,DENMAT)
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SET1COVERLAP(IAT_,LMNXX_,NSPINX_,DENMAT)
!     ******************************************************************
!     **  GET (<AEPHI|AEPHI>-<PSPHI|PSPHI>)<PS-P|PS-PSI>              **
!     ******************************************************************
      USE WAVES_MODULE
      INTEGER(4),INTENT(IN) :: IAT_
      INTEGER(4),INTENT(IN) :: LMNXX_
      INTEGER(4),INTENT(IN) :: NSPINX_
      REAL(8)   ,INTENT(IN) :: DENMAT(LMNXX_,LMNXX_,NSPINX_)
      INTEGER(4)            :: ISPIN,IKPT ! DO LOOP COUNTERS
      INTEGER(4)            :: ISP
      INTEGER(4)            :: LMNX1
!     ******************************************************************
      CALL ATOMLIST$GETI4('ISPECIES',IAT_,ISP)
      CALL SETUP$LMNX(ISP,LMNX1)
      IF(LMNX1.GT.LMNXX_) THEN
        CALL ERROR$STOP('WAVES$1COVERLAP')
      END IF
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          CALL WAVES_1CHAMILTONIAN(NAT,LMNXX_,LMNX1,NB,IAT_ &
     &              ,FNL(1,1,1,IKPT,ISPIN),ONL(1,1,1,IKPT,ISPIN) &
     &              ,DENMAT(1,1,ISPIN))
        ENDDO
      ENDDO
      CALL OPTIC3$DO(NSPINX_,LMNXX_,NAT,IAT_,DENMAT)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$RESETONLHNL
!     ******************************************************************
!     **  RESET FNL, HNL                                              **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)    :: ISPIN,IKPT,IB,IAT,LMN  ! DO LOOP COUNTERS
!     ******************************************************************
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IB = 1,NB
            DO IAT = 1,NAT
              DO LMN= 1,LMNXX
                ONL(IAT,IB,LMN,IKPT,ISPIN)=0.D0
                HNL(IAT,IB,LMN,IKPT,ISPIN)=0.D0
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$COMMUNICATEONLHNL
!     ******************************************************************
!     **  COMMUNICATE                                                 **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$COMMUNICATEONLHNL')
      CALL MPE$COMBINE('+',ONL)
      CALL MPE$COMBINE('+',HNL)
                           CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$1CFORCE
!     ******************************************************************
!     **  CALCULATE 1-CENTER CONTRIBUTION TO THE FORCE                **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)                :: OCC(NB,NKPT,NSPIN)        
      REAL(8)   ,ALLOCATABLE :: FORCE(:,:)   !(3,NAT)     
      REAL(8)   ,ALLOCATABLE :: FORCET(:,:)  !(3,NAT)     
      INTEGER(4),ALLOCATABLE :: ISPECIES(:)  !(NAT)       
      INTEGER(4),ALLOCATABLE :: LMNX(:)      !(NSP)       
      INTEGER(4)             :: ISP,I,IKPT,ISPIN,IAT
      INTEGER(4)             :: NTASKS,THISTASK
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$1CFORCE')
      CALL SETUP$NSPECIES(NSP)
      ALLOCATE(LMNX(NSP))      
      DO ISP=1,NSP
        CALL SETUP$LMNX(ISP,LMNX(ISP))
      ENDDO
      ALLOCATE(ISPECIES(NAT))  
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
!     == GET OCCUPATIONS =============================================
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
!     
!     ==================================================================
!     ==  CALCULATE FORCES DUE TO 1CENTER TERMS                       ==
!     ==================================================================
      ALLOCATE(FORCET(3,NAT)) 
      DO IAT=1,NAT
        DO I = 1,3
          FORCET(I,IAT)=0.D0
        ENDDO
      ENDDO
      CALL MPE$QUERY(NTASKS,THISTASK)
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IAT=THISTASK,NAT,NTASKS
            ISP=ISPECIES(IAT)              
            CALL WAVES_1CFORCE(IAT,NAT,NB,LMNXX &
     &                ,NB,LMNX(ISP),WKPT,OCC(:,IKPT,ISPIN) &
     &                ,RLAM0(:,:,IKPT,ISPIN) &
     &                ,DFNL(:,:,:,:,IKPT,ISPIN),HNL(:,:,:,IKPT,ISPIN) &
     &                ,ONL(:,:,:,IKPT,ISPIN),FORCET(:,IAT))
          ENDDO
        ENDDO
      ENDDO
      CALL MPE$COMBINE('+',FORCET)
!     
!     ==================================================================
!     ==  ADD THESE FORCES TO ATOMLIST                                ==
!     ==================================================================
      ALLOCATE(FORCE(3,NAT))   
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCE)
      DO IAT=1,NAT
        DO I=1,3
          FORCE(I,IAT)=FORCE(I,IAT)+FORCET(I,IAT)
        ENDDO
      ENDDO
      CALL ATOMLIST$SETR8A('FORCE',0,3*NAT,FORCE)
!     
!     ==================================================================
!     ==  PRINT FOR TEST                                              ==
!     ==================================================================
      IF(TPR) THEN
        PRINT*,'FORCE AFTER WAVES_1CFORCE' 
        DO IAT=1,NAT
          WRITE(*,FMT='(3F20.10)')(FORCE(I,IAT),I=1,3)
        ENDDO
      END IF
!     
!     ==================================================================
!     ==  CLOSE ROUTINE                                               ==
!     ==================================================================
      DEALLOCATE(FORCE)    
      DEALLOCATE(FORCET)   
      DEALLOCATE(ISPECIES) 
      DEALLOCATE(LMNX)     
                           CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$WAVEKINETIC(EKINC)
!     ******************************************************************
!     **  CALCULATE THE KINETIC ENERGY OF THE WAVE FUNCTIONS          **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT):: EKINC
      LOGICAL(4)            :: TGAMMA
      REAL(8)               :: EKINC1
      INTEGER(4)            :: IB,ISPIN,IKPT
      REAL(8)               :: OCC(NB,NKPT,NSPIN)  
      REAL(8)               :: EIG(NB,NKPT,NSPIN)  
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$WAVEKINETIC')
!     
!     ================================================================
!     ==  KINETIC ENERGY OF EACH BAND                               ==
!     ================================================================
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          TGAMMA=(DSQRT(HSG(1,IKPT)).LT.EPSILONGAMMA)
!         == AFTER PROPAGATE, HPSI CONTAINS PSI(-)
          IF(ALLOCATED(HPSI)) THEN
            CALL WAVES_1WAVEKINETIC(NGWLX,NGWL(IKPT),NB &
     &         ,TGAMMA,WKPT,CELLVOL,HPSI(1,1,IKPT,ISPIN),CM(1,1,IKPT,ISPIN) &
     &         ,2.D0*DELT,EMASS,EMASSCG2,HSG(1,IKPT),EIG(:,IKPT,ISPIN))
          ELSE
            CALL WAVES_1WAVEKINETIC(NGWLX,NGWL(IKPT),NB &
     &         ,TGAMMA,WKPT,CELLVOL,C0(1,1,IKPT,ISPIN),CM(1,1,IKPT,ISPIN) &
     &         ,DELT,EMASS,EMASSCG2,HSG(1,IKPT),EIG(:,IKPT,ISPIN))
          END IF
        ENDDO
      ENDDO
!     
!     ++ COMMUNICATE EKINC +++++++++++++++++++++++++++++++++++++++++++
      CALL MPE$COMBINE('+',EIG)
!     ++ END +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     
!     ================================================================
!     ==  PROVIDE KINETIC ENERGY TO OCCUPATIONS MODULE              ==
!     ================================================================
!EIG=0.D0
      CALL DYNOCC$SETR8A('M<PSIDOT|PSIDOT>',NB*NKPT*NSPIN,EIG)
!     
!     ================================================================
!     ==  SUM UP KINETIC ENERGY                                     ==
!     ================================================================
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
      EKINC=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IB=1,NB
            EKINC=EKINC+EIG(IB,IKPT,ISPIN)*OCC(IB,IKPT,ISPIN)*WKPT
          ENDDO
        ENDDO
      ENDDO
                           CALL TRACE$POP
      RETURN
      END
!
!     ....................................WAVES_1WAVEKINETIC............
      SUBROUTINE WAVES_1WAVEKINETIC(NGWX,NGW,NB &
     &                             ,TGAMMA,WKPT,CELLVOL,C0,CM &
     &                             ,DELT,EMASS,EMASSCG2,G2,EIG)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE KINETIC ENERGY OF THE FICTITIOUS ELECTRONIC       **
!     **  VARIABLES BETWEEN C(0) AND C(-)                             **
!     **  REMARK : THIS FORMULA IS NOT CONSISTENT WITH THE VERLET     **
!     **  ALGORITHM                                                   **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NGWX
      INTEGER(4),INTENT(IN)    :: NGW              ! NUMBER OF G-VECTORS
      INTEGER(4),INTENT(IN)    :: NB               ! NUMBER OF STATES
      LOGICAL(4),INTENT(IN)    :: TGAMMA           !
      REAL(8)   ,INTENT(IN)    :: WKPT             ! K-POINT  WEIGHT
      REAL(8)   ,INTENT(IN)    :: CELLVOL            ! UNIT-CELL VOLUME
      REAL(8)   ,INTENT(IN)    :: DELT             ! TIME STEP
      REAL(8)   ,INTENT(IN)    :: EMASS            ! WAVE FUNCTION MASS
      REAL(8)   ,INTENT(IN)    :: EMASSCG2         
      REAL(8)   ,INTENT(IN)    :: G2(NGWX)         ! G**2
      COMPLEX(8),INTENT(IN)    :: C0(NGWX,NB)      ! WAVE FUNCTION (0)
      COMPLEX(8),INTENT(IN)    :: CM(NGWX,NB)      ! WAVE FUNCTION (-)
      REAL(8)   ,INTENT(OUT)   :: EIG(NB)          ! KINETIC ENERGY
      REAL(8)                  :: SVAR,DEKINC
      COMPLEX(8)               :: SPEED            
      INTEGER(4)               :: IB,IG            ! RUNNING INDICES
!     ******************************************************************
!
!     ==================================================================
!     ==  PLANE WAVE CONTRIBUTION                                     ==
!     ==================================================================
      DO IB=1,NB
        DEKINC=0.D0
        IF(TGAMMA) THEN
          SPEED=C0(1,IB)-CM(1,IB)
          DEKINC=-0.5D0*DBLE(SPEED*CONJG(SPEED))*(1.D0+EMASSCG2*G2(1))
        END IF
        DO IG=1,NGW
          SPEED=C0(IG,IB)-CM(IG,IB)
          DEKINC=DEKINC+DBLE(SPEED*CONJG(SPEED))*(1.D0+EMASSCG2*G2(IG))
        ENDDO
        EIG(IB)=DEKINC*CELLVOL*(2.D0*EMASS/DELT**2)
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$HPSI(NNRX_,NSPIN_,POT)
!     ******************************************************************
!     **  EVALUATE HPSI                                               **
!     ******************************************************************
      USE WAVES_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NNRX_
      INTEGER(4),INTENT(IN) :: NSPIN_
      REAL(8)   ,INTENT(IN) :: POT(NNRX_,NSPIN_)
      LOGICAL               :: TGAMMA
      REAL(8)   ,ALLOCATABLE:: DWORK(:,:,:,:) !(NB,NB,NKPT,NSPIN)
      INTEGER(4),ALLOCATABLE:: LOX(:,:)       !(LNXX,NSP)           
      INTEGER(4),ALLOCATABLE:: ISPECIES(:)    !(NAT)          
      INTEGER(4),ALLOCATABLE:: LNX(:)         !(NSP)          
      REAL(8)               :: EIG(NB,NKPT,NSPIN)
      REAL(8)   ,ALLOCATABLE:: HAMILTON(:,:,:,:)
      INTEGER(4)            :: ISP,ISPIN,IKPT,IB
      INTEGER(4)            :: LNXX
      REAL(8)               :: SVAR
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$HPSI')
      IF(NNRX_.LT.NR1L*NR2*NR3) THEN
        CALL ERROR$MSG('R-DIMENSIONS OF DENSITY INCONSISTENT')
        CALL ERROR$STOP('WAVES$HPSI')
      END IF
      IF(NSPIN_.LT.NSPIN) THEN
        CALL ERROR$MSG('SPIN-DIMENSION OF DENSITY INCONSISTENT')
        CALL ERROR$STOP('WAVES$HPSI')
      END IF
!     
!     ==================================================================
!     ==  HERE THE ARRAY STORE AND HPSI ARE OVERLAYED===================
!     ==================================================================
      IF(.NOT.ALLOCATED(STORE)) THEN
        CALL ERROR$MSG('ARRAY STORE DOES NOT EXIST')
        CALL ERROR$STOP('WAVES$HPSI')
      END IF
      ALLOCATE(HPSI(NGWLX,NB,NKPT,NSPIN))
!     
!     ==================================================================
!     ==  NOW COLLECT INFORMATION ABOUT SETUPS                        ==
!     ==================================================================
      CALL SETUP$NSPECIES(NSP)
      CALL SETUP$LMNXX(LMNXX)
      CALL SETUP$LNXX(LNXX)
      ALLOCATE(LNX(NSP))        
      ALLOCATE(LOX(LNXX,NSP))   
      DO ISP=1,NSP
        CALL SETUP$LNX(ISP,LNX(ISP))
        CALL SETUP$LOFLN(ISP,LNX(ISP),LOX(1,ISP))
      ENDDO
      ALLOCATE(ISPECIES(NAT))   
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
!     
!     ==================================================================
!     ==  SEND DATA TO OPTIC CODE                                     ==
!     ==================================================================
      CALL OPTIC3$WAVES(NGWLX,NB,NKPT,NSPIN,C0)   
!     
!     ==================================================================
!     ==  CALCULATE H-TILDE|PSI-TILDE>                                ==
!     ==================================================================
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          CALL PLANEWAVE$SELECT('WAVE',IKPT)
          CALL WAVES_HPSI(NGWLX,NRSTORE &
     &              ,NB,NGWL(IKPT),NR1L,NR2,NR3,HSG(1,IKPT) &
     &              ,C0(1,1,IKPT,ISPIN),POT(1,ISPIN) &
     &             ,HPSI(1,1,IKPT,ISPIN),STORE(1,1,IKPT,ISPIN))
           CALL WAVES_ADDPRO(NGWLX,LMNXX,LNXX &
     &              ,NSP,NAT,ISPECIES,NB,NGWL(IKPT),LNX,LOX &
     &              ,HPSI(1,1,IKPT,ISPIN),WNL(1,1,1,IKPT) &
     &              ,HNL(1,1,1,IKPT,ISPIN),EIGR(1,1,IKPT))
          TGAMMA=(DSQRT(HSG(1,IKPT)).LT.EPSILONGAMMA)
          CALL WAVES_CORRECTGAMMA(TGAMMA,NGWLX,NB,HPSI(1,1,IKPT,ISPIN))
        ENDDO
      ENDDO
      DEALLOCATE(ISPECIES) 
      DEALLOCATE(LOX)      
      DEALLOCATE(LNX)      
      DEALLOCATE(STORE)
!     
!     ================================================================
!     ==  ENERGY EXPECTATION VALUES OF THE HAMILTON OPERATOR        ==
!     ================================================================
      IF(THAMILTON) THEN
        IF(.NOT.ALLOCATED(HAMU))ALLOCATE(HAMU(NB,NB,NKPT,NSPIN))
        IF(.NOT.ALLOCATED(HAME))ALLOCATE(HAME(NB,NKPT,NSPIN))
        ALLOCATE(HAMILTON(NB,NB,NKPT,NSPIN))
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            TGAMMA=(DSQRT(HSG(1,IKPT)).LT.EPSILONGAMMA)
            CALL WAVES_SIGSET(TGAMMA,0,C0(1,1,IKPT,ISPIN) &
     &             ,HPSI(1,1,IKPT,ISPIN),NB,NB,NGWLX,NGWL(IKPT),CELLVOL &
     &             ,HAMILTON(1,1,IKPT,ISPIN))
          ENDDO
        ENDDO
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
!           == EIG ARE THE EXPECTATION VALUES OF THE HAMILTONIAN,     ==
!           ==  NOT ITS EIGNEVALUES
            DO IB=1,NB
              EIG(IB,IKPT,ISPIN)=HAMILTON(IB,IB,IKPT,ISPIN)
            ENDDO
          ENDDO
        ENDDO
        CALL MPE$COMBINE('+',HAMILTON)
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            CALL DIAG(NB,NB,HAMILTON(1,1,IKPT,ISPIN) &
     &               ,HAME(1,IKPT,ISPIN),HAMU(1,1,IKPT,ISPIN))
          ENDDO
        ENDDO
        DEALLOCATE(HAMILTON)
      ELSE
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB=1,NB
              TGAMMA=(DSQRT(HSG(1,IKPT)).LT.EPSILONGAMMA)
              CALL WAVES_SIGSET(TGAMMA,0,C0(1,IB,IKPT,ISPIN) &
     &             ,HPSI(1,IB,IKPT,ISPIN),1,1,NGWLX,NGWL(IKPT),CELLVOL &
     &             ,EIG(IB,IKPT,ISPIN))
            ENDDO
          ENDDO
        ENDDO
      END IF
      CALL MPE$COMBINE('+',EIG)
!PRINT*,'THAMILTON',THAMILTON
!CALL CONSTANTS$GET('EV',SVAR)
!WRITE(*,FMT='("EIG",10F8.3)')EIG/SVAR
      CALL DYNOCC$SETR8A('EPSILON',NB*NKPT*NSPIN,EIG)
!     
!     ================================================================
!     ==  REORDER STATES                                            ==
!     ================================================================
      IF(.NOT.TSAFEORTHO) THEN
        CALL WAVES_ORDERSTATES
      END IF
                           CALL TRACE$POP
      RETURN
      CONTAINS
!     ..................................................................
      SUBROUTINE WAVES_ORDERSTATES
!     ******************************************************************
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,PARAMETER :: TOL=0.D0
      INTEGER(4)           :: ISPIN
      INTEGER(4)           :: IKPT
      INTEGER(4)           :: LMN
      INTEGER(4)           :: FROM,TO
      REAL(8)              :: TMP1(NB),TMP2(NB),TMP3(NB),TMP4(NB),TMP(NB)
      INTEGER(4)           :: I
!     ******************************************************************
PRINT*,'ORDERING STATES'
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
!
!         ==============================================================
!         == DEFINE ORDERING                                          ==
!         ==============================================================
          TMP(:)=EIG(:,IKPT,ISPIN)
!         WRITE(*,FMT='("EIG",9F8.5)')EIG(:,IKPT,ISPIN)*27.211D0
          DO I=1,NB-1
            IF(DABS(TMP(I+1)-TMP(I)).LT.TOL) TMP(I+1)=TMP(I)
          ENDDO
!         WRITE(*,FMT='("TMP",9F8.5)')TMP(:)*27.211D0
          CALL SORT$SET(NB,TMP(:))
!
!         ==============================================================
!         == ORDER WAVE FUNCTIONS AND RELATED                         ==
!         ==============================================================
          CALL SORT$ORDERC16(NGWLX,NB,C0(:,:,IKPT,ISPIN))
          CALL SORT$ORDERC16(NGWLX,NB,CM(:,:,IKPT,ISPIN))
          CALL SORT$ORDERC16(NGWLX,NB,HPSI(:,:,IKPT,ISPIN))
          IF(ALLOCATED(OPSI)) THEN
            CALL SORT$ORDERC16(NGWLX,NB,OPSI(:,:,IKPT,ISPIN))
          END IF
!
!         ==============================================================
!         == ORDER PROJECTIONS AND RELATED                            ==
!         ==============================================================
          DO LMN=1,LMNXX
            CALL SORT$ORDERR8(NAT,NB,FNL(:,:,LMN,IKPT,ISPIN))
            CALL SORT$ORDERR8(NAT,NB,ONL(:,:,LMN,IKPT,ISPIN))
            CALL SORT$ORDERR8(NAT,NB,DFNL(:,:,LMN,1,IKPT,ISPIN))
            CALL SORT$ORDERR8(NAT,NB,DFNL(:,:,LMN,2,IKPT,ISPIN))
            CALL SORT$ORDERR8(NAT,NB,DFNL(:,:,LMN,3,IKPT,ISPIN))
          ENDDO
!
!         ==============================================================
!         == ORDER LAGRANGE MULTIPLIERS                               ==
!         ==============================================================
          CALL SORT$ORDERR8(NB,NB,RLAM0(:,:,IKPT,ISPIN))
          CALL SORT$ORDERR8(NB,NB,RLAMM(:,:,IKPT,ISPIN))
          IF(ALLOCATED(RLAM2M)) THEN
            CALL SORT$ORDERR8(NB,NB,RLAM2M(:,:,IKPT,ISPIN))
            IF(ALLOCATED(RLAM3M)) THEN
              CALL SORT$ORDERR8(NB,NB,RLAM3M(:,:,IKPT,ISPIN))
            END IF
          END IF
          CALL SORT$RESTART
          CALL SORT$FLIP(FROM,TO)
          DO WHILE(FROM.NE.0.OR.TO.NE.0) 
            PRINT*,'FLIP FROM STATE ',FROM,' TO STATE ',TO
            IF(FROM.EQ.0) THEN
              RLAM0(TO,:,IKPT,ISPIN)=TMP1(:)
              RLAMM(TO,:,IKPT,ISPIN)=TMP2(:)
              IF(ALLOCATED(RLAM2M)) THEN
                RLAM2M(TO,:,IKPT,ISPIN)=TMP3(:)
                IF(ALLOCATED(RLAM3M)) THEN
                  RLAM3M(TO,:,IKPT,ISPIN)=TMP4(:)
                END IF
              END IF
            ELSE IF(TO.EQ.0) THEN
              TMP1(:)=RLAM0(FROM,:,IKPT,ISPIN)
              TMP2(:)=RLAMM(FROM,:,IKPT,ISPIN)
              IF(ALLOCATED(RLAM2M)) THEN
                TMP3(:)=RLAM2M(FROM,:,IKPT,ISPIN)
                IF(ALLOCATED(RLAM3M)) THEN
                  TMP4(:)=RLAM3M(FROM,:,IKPT,ISPIN)
                END IF
              END IF
            ELSE
              RLAM0(TO,:,IKPT,ISPIN)=RLAM0(FROM,:,IKPT,ISPIN)
              RLAMM(TO,:,IKPT,ISPIN)=RLAMM(FROM,:,IKPT,ISPIN)
              IF(ALLOCATED(RLAM2M)) THEN
                RLAM2M(TO,:,IKPT,ISPIN)=RLAM2M(FROM,:,IKPT,ISPIN)
                IF(ALLOCATED(RLAM3M)) THEN
                  RLAM3M(TO,:,IKPT,ISPIN)=RLAM3M(FROM,:,IKPT,ISPIN)
                END IF
              END IF
            END IF
            CALL SORT$FLIP(FROM,TO)
          ENDDO
!         == REARRANGE BANDS IN DYNOCC OBJECT ==========================
          CALL DYNOCC$ORDER(IKPT,ISPIN)
!         == UNSET SORT OBJECT =========================================
          CALL SORT$UNSET
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE WAVES_ORDERSTATES
      END
!
!     ..................................................................
      SUBROUTINE WAVES$OPSI
!     ******************************************************************
!     **  EVALUATE OPSI                                               **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL(4)              :: TGAMMA
      INTEGER(4),ALLOCATABLE  :: LOX(:,:)    !(LNXX,NSP)   
      INTEGER(4),ALLOCATABLE  :: ISPECIES(:) !(NAT)   
      INTEGER(4),ALLOCATABLE  :: LNX(:)      !(NSP)   
      INTEGER(4)              :: LNXX
      INTEGER(4)              :: ISP,ISPIN,IKPT
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$OPSI')
!     
      CALL SETUP$NSPECIES(NSP)
      CALL SETUP$LNXX(LNXX)
      CALL SETUP$LMNXX(LMNXX)
      IF(.NOT.ALLOCATED(OPSI)) THEN
        ALLOCATE(OPSI(NGWLX,NB,NKPT,NSPIN))
      ELSE
        CALL ERROR$MSG('ARRAY OPSI ALREADY OCCUPIED')
        CALL ERROR$STOP('WAVES$OPSI')
      END IF
      ALLOCATE(LNX(NSP))
      ALLOCATE(LOX(LNXX,NSP))
      DO ISP=1,NSP
        CALL SETUP$LNX(ISP,LNX(ISP))
        CALL SETUP$LOFLN(ISP,LNX(ISP),LOX(1,ISP))
      ENDDO
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          CALL WAVES_SETEQUAL(NGWLX,NGWL(IKPT) &
     &              ,NB,C0(1,1,IKPT,ISPIN),OPSI(1,1,IKPT,ISPIN))
          CALL WAVES_ADDPRO(NGWLX,LMNXX,LNXX &
     &              ,NSP,NAT,ISPECIES,NB,NGWL(IKPT),LNX,LOX &
     &              ,OPSI(1,1,IKPT,ISPIN),WNL(1,1,1,IKPT) &
     &              ,ONL(1,1,1,IKPT,ISPIN),EIGR(1,1,IKPT))
          TGAMMA=(DSQRT(HSG(1,IKPT)).LT.EPSILONGAMMA)
          CALL WAVES_CORRECTGAMMA(TGAMMA,NGWLX,NB &
     &                           ,OPSI(1,1,IKPT,ISPIN))
          CALL WAVES_BYPSIMASS(NGWLX,NGWL(IKPT),NB,EMASS,EMASSCG2 &
     &                        ,HSG(1,IKPT),OPSI(1,1,IKPT,ISPIN))
        ENDDO
      ENDDO
      DEALLOCATE(ISPECIES)
      DEALLOCATE(LNX)
      DEALLOCATE(LOX)
                           CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$PROPAGATE
!     ******************************************************************
!     **  PROPAGATE WAVE FUNCTIONS WITHOUT ORTHOGONALIZATION FORCES   **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL    :: TGAMMA
      INTEGER(4) :: IKPT,ISPIN
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$PROPAGATE')
!
!     ==================================================================
!     == STOP WAVE FUNCTIONS AND RANDOMIZE                            ==
!     ==================================================================
      IF(TSTOP) THEN
        CALL WAVES$STOP
        RLAMM(:,:,:,:)=RLAM0(:,:,:,:)
        IF(ALLOCATED(RLAM2M)) THEN
          RLAM2M(:,:,:,:)=RLAMM(:,:,:,:)
          IF(ALLOCATED(RLAM3M)) THEN
            RLAM3M(:,:,:,:)=RLAM2M(:,:,:,:)
          END IF            
        END IF
        TSTOP=.FALSE.
      END IF
      CALL WAVES$RANDOMIZEVELOCITY
!
!     ==================================================================
!     == PROPAGATE WAVE FUNCTIONS                                     ==
!     ==================================================================
!     == NOTE FORCES HPSI MUST BE CALCULATED BEFORE!!  =================
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
!         == ON OUTPUT: HPSI CONTAINS PSI(-)
          CALL WAVES_PROPAGATE(NGWLX,NGWL(IKPT),NB &
     &              ,EMASS,EMASSCG2,HSG(1,IKPT),DELT,ANNEE &
     &              ,C0(1,1,IKPT,ISPIN),CM(1,1,IKPT,ISPIN) &
     &              ,HPSI(1,1,IKPT,ISPIN))
!         ==  C(G=0) IS MADE REAL ======================================
          TGAMMA=(DSQRT(HSG(1,IKPT)).LT.EPSILONGAMMA)
          CALL WAVES_CORRECTGAMMA(TGAMMA,NGWLX,NB,CM(1,1,IKPT,ISPIN))
        ENDDO
      ENDDO
      DEALLOCATE(HPSI)
!
      CALL TRACE$POP 
      RETURN
      END
!
!     ......................................................RANWAV......
      SUBROUTINE WAVES_PROPAGATE(NGWX,NGW,NB &
     &                ,EMASS,EMASSCG2,G2,DELT,ANNEE,PSI0,PSI1,HPSI)
!     ******************************************************************
!     **                                                              **
!     **  PROPAGATES WAVE FUNCTION ACCORDING TO FORCES                **
!     **  (NO CONSTRAINTS)                                            **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NGW,NGWX       ! #(G-VECTORS); (MAX)
      INTEGER(4),INTENT(IN)    :: NB             ! #(BANDS)
      REAL(8)   ,INTENT(IN)    :: G2(NGWX)       ! G**2
      REAL(8)   ,INTENT(IN)    :: ANNEE          ! FRICTION
      REAL(8)   ,INTENT(IN)    :: EMASS          ! WAVE FUNCTION MASS
      REAL(8)   ,INTENT(IN)    :: EMASSCG2      
      REAL(8)   ,INTENT(IN)    :: DELT           ! TIME STEP
      COMPLEX(8),INTENT(IN)    :: PSI0(NGWX,NB)  ! |PSI(0)>
      COMPLEX(8),INTENT(INOUT) :: HPSI(NGWX,NB)  ! H|PSI(0)>
      COMPLEX(8),INTENT(INOUT) :: PSI1(NGWX,NB)  ! IN:|PSI(-)>; OUT: |PSI-BAR>
      REAL(8)                  :: SVAR1,SVAR2,SVAR3,FAC ! AUXILIARY VARIABLES
      INTEGER(4)               :: IB,IG                 ! RUNNING INDICES
      COMPLEX(8)               :: CSVAR
!     ******************************************************************
      SVAR1=2.D0/(1.D0+ANNEE)
      SVAR2=1.D0-SVAR1
      FAC=-DELT**2 / EMASS / (1.D0+ANNEE)
      DO IB=1,NB
        DO IG=1,NGW
          SVAR3=FAC/(1.D0+EMASSCG2*G2(IG))
          CSVAR=PSI1(IG,IB)
          PSI1(IG,IB) = SVAR1*PSI0(IG,IB) + SVAR2*PSI1(IG,IB) &
     &                                    + SVAR3*HPSI(IG,IB)
          HPSI(IG,IB)=CSVAR
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$PROJECTIONS(STRING,IDIM)
!     ******************************************************************
!     **  GENERATE STRUCTURE FACTORS                                  **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(*),INTENT(IN) :: STRING
      INTEGER(4)  ,INTENT(IN) :: IDIM
      REAL(8)                 :: OCC(NB,NKPT,NSPIN) 
      INTEGER(4)  ,ALLOCATABLE:: LNX(:)      !(NSP)      
      INTEGER(4)  ,ALLOCATABLE:: LOX(:,:)    !(LNXX,NSP) 
      INTEGER(4)  ,ALLOCATABLE:: ISPECIES(:) !(NAT)      
      INTEGER(4)  ,ALLOCATABLE:: LMNX(:)     !(NSP)      
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$PROJECTIONS')
!       == SELECT C0,CM
!       == GRADIENT OR NOT
        CALL SETUP$NSPECIES(NSP)
        CALL SETUP$LNXX(LNXX)
        CALL SETUP$LMNXX(LMNXX)
        ALLOCATE(LNX(NSP))
        ALLOCATE(LMNX(NSP))
        ALLOCATE(LOX(LNXX,NSP))
        DO ISP=1,NSP
          CALL SETUP$LNX(ISP,LNX(ISP))
          CALL SETUP$LMNX(ISP,LMNX(ISP))
          CALL SETUP$LOFLN(ISP,LNX(ISP),LOX(1,ISP))
        ENDDO
        CALL ATOMLIST$NATOM(NAT)
        ALLOCATE(ISPECIES(NAT))
        CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
        IF(STRING.EQ.'PSI(0)'.AND.IDIM.EQ.0) THEN
          IGRAD=0
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              CALL WAVES_PROJECTIONS(LNXX,LMNXX,NGWLX &
     &                 ,NSP,NAT,ISPECIES,NB,LNX,LOX,CELLVOL &
     &                 ,NGWL(IKPT),HSG(1,IKPT),GSK(1,1,IKPT) &
     &                 ,WNL(1,1,1,IKPT),EIGR(1,1,IKPT) &
     &                 ,C0(1,1,IKPT,ISPIN),IGRAD,FNL(1,1,1,IKPT,ISPIN))
            ENDDO
          ENDDO
!         ++ SUM PROJECTIONS +++++++++++++++++++++++++++++++++++++++++++
!         MLEN = 8*NAT*NB*LMNXX*NKPT*NSPIN
!         CALL MPE$COMBINE('ADD_R*8',MLEN,FNL)
          CALL WAVES_PROJECTIONCOMBINE(1,NSPIN,NKPT,NB,NAT,ISPECIES &
     &                                  ,LMNXX,NSP,LMNX,FNL)
!         ++ END  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!        == TRANSFER PROJECTIONS TO PLOTTING ROUTINE ===================
          CALL STATEANALYSIS$SETPROJECTIONS(NAT,NB,LMNXX &
     &                  ,NKPT,NSPIN,FNL)
        ELSE IF(STRING.EQ.'PSI(0)'.AND.IDIM.EQ.1) THEN
          IGRAD=1
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              CALL WAVES_PROJECTIONS(LNXX,LMNXX,NGWLX &
     &                 ,NSP,NAT,ISPECIES,NB,LNX,LOX,CELLVOL &
     &                 ,NGWL(IKPT),HSG(1,IKPT),GSK(1,1,IKPT) &
     &                 ,WNL(1,1,1,IKPT),EIGR(1,1,IKPT) &
     &                ,C0(1,1,IKPT,ISPIN),IGRAD,DFNL(1,1,1,1,IKPT,ISPIN))
            ENDDO
          ENDDO
!         ++ SUM PROJECTIONS +++++++++++++++++++++++++++++++++++++++++++
!         MLEN = 8*3*NAT*NB*LMNXX*NKPT*NSPIN
!         CALL MPE$COMBINE('ADD_R*8',MLEN,DFNL)
          CALL WAVES_PROJECTIONCOMBINE(3,NSPIN,NKPT,NB,NAT,ISPECIES &
     &                                  ,LMNXX,NSP,LMNX,DFNL)
!         ++ END  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ELSE IF(STRING.EQ.'PSI(-)'.AND.IDIM.EQ.0) THEN
          IGRAD=0
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              CALL WAVES_PROJECTIONS(LNXX,LMNXX,NGWLX &
     &                 ,NSP,NAT,ISPECIES,NB,LNX,LOX,CELLVOL &
     &                 ,NGWL(IKPT),HSG(1,IKPT),GSK(1,1,IKPT) &
     &                 ,WNL(1,1,1,IKPT),EIGR(1,1,IKPT) &
     &                 ,CM(1,1,IKPT,ISPIN),IGRAD,ONL(1,1,1,IKPT,ISPIN))
            ENDDO
          ENDDO
!         ++ SUM PROJECTIONS +++++++++++++++++++++++++++++++++++++++++++
!         MLEN = 8*NAT*NB*LMNXX*NKPT*NSPIN
!         CALL MPE$COMBINE('ADD_R*8',MLEN,ONL)
          CALL WAVES_PROJECTIONCOMBINE(1,NSPIN,NKPT,NB,NAT,ISPECIES &
     &                                  ,LMNXX,NSP,LMNX,ONL)
!         ++ END  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ELSE IF(STRING.EQ.'PSI(+)'.AND.IDIM.EQ.0) THEN
          IGRAD=0
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              CALL WAVES_PROJECTIONS(LNXX,LMNXX,NGWLX &
     &                 ,NSP,NAT,ISPECIES,NB,LNX,LOX,CELLVOL &
     &                 ,NGWL(IKPT),HSG(1,IKPT),GSK(1,1,IKPT) &
     &                 ,WNL(1,1,1,IKPT),EIGR(1,1,IKPT) &
     &                 ,CM(1,1,IKPT,ISPIN),IGRAD,FNL(1,1,1,IKPT,ISPIN))
            ENDDO
          ENDDO
!         ++ SUM PROJECTIONS +++++++++++++++++++++++++++++++++++++++++++
!         MLEN = 8*NAT*NB*LMNXX*NKPT*NSPIN
!         CALL MPE$COMBINE('ADD_R*8',MLEN,FNL)
          CALL WAVES_PROJECTIONCOMBINE(1,NSPIN,NKPT,NB,NAT,ISPECIES &
     &                                  ,LMNXX,NSP,LMNX,FNL)
!         ++ END  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ELSE IF(STRING.EQ.'O*PSI'.AND.IDIM.EQ.0) THEN
          IF(.NOT.ALLOCATED(OPSI)) THEN
            CALL ERROR$MSG('OPSI IS NOT ALLOCATED')
            CALL ERROR$STOP('WAVES$PROJECTIONS')
          END IF
          IGRAD=0
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              CALL WAVES_PROJECTIONS(LNXX,LMNXX,NGWLX &
     &                 ,NSP,NAT,ISPECIES,NB,LNX,LOX,CELLVOL &
     &                 ,NGWL(IKPT),HSG(1,IKPT),GSK(1,1,IKPT) &
     &                 ,WNL(1,1,1,IKPT),EIGR(1,1,IKPT) &
     &                 ,OPSI(1,1,IKPT,ISPIN),IGRAD,ONL(1,1,1,IKPT,ISPIN))
            ENDDO
          ENDDO
!         ++ SUM PROJECTIONS +++++++++++++++++++++++++++++++++++++++++++
!         MLEN = 8*NAT*NB*LMNXX*NKPT*NSPIN
!         CALL MPE$COMBINE('ADD_R*8',MLEN,ONL)
          CALL WAVES_PROJECTIONCOMBINE(1,NSPIN,NKPT,NB,NAT,ISPECIES &
     &                                  ,LMNXX,NSP,LMNX,ONL)
!         ++ END  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ELSE
          CALL ERROR$MSG('SELECTION NOT ALLOWED')
          CALL ERROR$STOP('WAVES$PROJECTIONS')
        END IF
        DEALLOCATE(ISPECIES)
        DEALLOCATE(LNX)
        DEALLOCATE(LOX)
        DEALLOCATE(LMNX)
                             CALL TRACE$POP
        RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$ORTHOGONALIZE
!     ******************************************************************
!     **  WAVES$ORHOGONALIZE                                          **
!     **  SWITCH THE TWO WAVE FUNCTIONS                               **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      LOGICAL(4)             :: TGAMMA
      REAL(8)                :: OCC(NB,NKPT,NSPIN)      
      REAL(8)  ,ALLOCATABLE  :: DOVER(:,:,:)  !(LNXX,LNXX,NSP)
      INTEGER(4),ALLOCATABLE :: LOX(:,:)      !(LNXX,NSP)
      INTEGER(4),ALLOCATABLE :: ISPECIES(:)   !(NAT)
      INTEGER(4),ALLOCATABLE :: LMNX(:)       !(NSP)
      INTEGER(4),ALLOCATABLE :: LNX(:)        !(NSP)
      INTEGER(4)             :: LNXX
      INTEGER(4)             :: ISPIN,IKPT,ISP
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$ORTHOGONALIZE')
      IF(.NOT.ALLOCATED(OPSI)) THEN
        CALL ERROR$MSG('OPSI NOT ALLOCATED')
        CALL ERROR$STOP('WAVES$ORTHOGONALIZE')
      END IF
!     
      CALL SETUP$NSPECIES(NSP)
      CALL SETUP$LNXX(LNXX)
      ALLOCATE(DOVER(LNXX,LNXX,NSP))
      ALLOCATE(LNX(NSP))
      ALLOCATE(LMNX(NSP))
      ALLOCATE(LOX(LNXX,NSP))
      DO ISP=1,NSP
        CALL SETUP$LNX(ISP,LNX(ISP))
        CALL SETUP$LMNX(ISP,LMNX(ISP))
        CALL SETUP$LOFLN(ISP,LNX(ISP),LOX(1,ISP))
        CALL SETUP$1COVERLAP(ISP,LNXX,DOVER(1,1,ISP))
      ENDDO
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
!     == GET OCCUPATIONS =============================================
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
!     
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          TGAMMA=(DSQRT(HSG(1,IKPT)).LT.EPSILONGAMMA)
          CALL WAVES_ORTHO(NB,LMNXX,LNXX,NGWLX &
     &              ,NSP,NAT,ISPECIES,LNX,LMNX,LOX,CELLVOL &
     &              ,TGAMMA,NB,OCC(1,IKPT,ISPIN),NGWL(IKPT) &
     &              ,OPSI(1,1,IKPT,ISPIN),CM(1,1,IKPT,ISPIN) &
     &              ,ONL(1,1,1,IKPT,ISPIN),FNL(1,1,1,IKPT,ISPIN) &
     &              ,RLAM0(1,1,IKPT,ISPIN),DOVER,DELT,EMASS,ANNEE,TSAFEORTHO)
        ENDDO
      ENDDO
      DEALLOCATE(ISPECIES)
      DEALLOCATE(LOX)
      DEALLOCATE(LMNX)
      DEALLOCATE(LNX)
      DEALLOCATE(DOVER)
      DEALLOCATE(OPSI)
!
!     ==================================================================
!     == PROPAGATE LAGRANGE PARAMETERS FOR ORTHOGONALIZATION          ==
!     == RLAM IS NOT A DYNAMICAL QUANTITY BUT IS RECALCULATED         ==
!     == DURING THE ORTHOGONALIZATION                                 ==
!     ==================================================================
      IF(ALLOCATED(RLAM2M)) THEN
        IF(ALLOCATED(RLAM3M)) THEN
          RLAMP(:,:,:,:)=4.D0*RLAM0-6.D0*RLAMM+4.D0*RLAM2M-RLAM3M
        ELSE
          RLAMP(:,:,:,:)=3.D0*RLAM0-3.D0*RLAMM+RLAM2M
        END IF
      ELSE
        RLAMP(:,:,:,:)=2.D0*RLAM0-RLAMM
      END IF
!
                           CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$GRAMMSCHMIDT(STRING)
!     ******************************************************************
!     **  GRAMM SCHMIDT ORTHOGONALIZATION                             **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRING
      LOGICAL                     :: TGAMMA
      REAL(8)         ,ALLOCATABLE:: DOVER(:,:,:) !(LNXX,LNXX,NSP)
      INTEGER(4)       ,ALLOCATABLE:: LOX(:,:)     !(LNXX,NSP)
      INTEGER(4)       ,ALLOCATABLE:: ISPECIES(:)  !(NAT)
      INTEGER(4)       ,ALLOCATABLE:: LMNX(:)      !(NSP)
      INTEGER(4)       ,ALLOCATABLE:: LNX(:)       !(NSP)
      INTEGER(4)                   :: LNXX
      INTEGER(4)                   :: ISP,ISPIN,IKPT
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$GRAMMSCHMIDT')
        CALL SETUP$NSPECIES(NSP)
        CALL SETUP$LNXX(LNXX)
        ALLOCATE(LNX(NSP))             
        ALLOCATE(LMNX(NSP))            
        ALLOCATE(LOX(LNXX,NSP))        
        ALLOCATE(DOVER(LNXX,LNXX,NSP)) 
        DO ISP=1,NSP
          CALL SETUP$LNX(ISP,LNX(ISP))
          CALL SETUP$LMNX(ISP,LMNX(ISP)) 
          CALL SETUP$LOFLN(ISP,LNX(ISP),LOX(1,ISP))
          CALL SETUP$1COVERLAP(ISP,LNXX,DOVER(1,1,ISP))
        ENDDO
        ALLOCATE(ISPECIES(NAT)) 
        CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)        
        IF(STRING.EQ.'PSI(0)') THEN
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              TGAMMA=(DSQRT(HSG(1,IKPT)).LT.EPSILONGAMMA) 
              CALL WAVES_GRAMMSCHMIDT(TGAMMA,NB,NGWL(IKPT) &
     &                  ,C0(1,1,IKPT,ISPIN),FNL(1,1,1,IKPT,ISPIN),DOVER &
     &                  ,NGWLX,LMNXX,LNXX &
     &                  ,NSP,NAT,ISPECIES,LNX,LMNX,LOX,CELLVOL)
           ENDDO
         ENDDO
       ELSE IF(STRING.EQ.'PSI(-)') THEN
          DO ISPIN=1,NSPIN
            DO IKPT=1,NKPT
              TGAMMA=(DSQRT(HSG(1,IKPT)).LT.EPSILONGAMMA) 
              CALL WAVES_GRAMMSCHMIDT(TGAMMA,NB,NGWL(IKPT) &
     &                  ,CM(1,1,IKPT,ISPIN),ONL(1,1,1,IKPT,ISPIN),DOVER &
     &                  ,NGWLX,LMNXX,LNXX &
     &                  ,NSP,NAT,ISPECIES,LNX,LMNX,LOX,CELLVOL)
           ENDDO
         ENDDO
       ELSE
         CALL ERROR$MSG('OPTION FOR STRING NOT ALLOWED')
         CALL ERROR$STOP('WAVES$GRAMMSCHMIDT')
       END IF
       DEALLOCATE(ISPECIES) 
       DEALLOCATE(LOX)      
       DEALLOCATE(LMNX)     
       DEALLOCATE(LNX)      
       DEALLOCATE(DOVER)    
                             CALL TRACE$POP
       RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$SWITCH
!     ******************************************************************
!     **  SWITCH THE TWO WAVE FUNCTIONS                               **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: ISPIN
      INTEGER(4)             :: IKPT
      REAL(4)   ,ALLOCATABLE :: DWORK(:,:,:,:)
!     ******************************************************************
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          CALL WAVES_SWITCH(NGWLX,NGWL(IKPT),NB &
     &              ,C0(1,1,IKPT,ISPIN),CM(1,1,IKPT,ISPIN))
        ENDDO
      ENDDO
!
!     == SWITCH LAGRANGE PARAMETERS ====================================
      IF(.NOT.ALLOCATED(RLAM2M)) THEN
        ALLOCATE(RLAM2M(NB,NB,NKPT,NSPIN))
      ELSE
        IF(.NOT.ALLOCATED(RLAM3M)) THEN
!         ALLOCATE(RLAM3M(NB,NB,NKPT,NSPIN))
        END IF
!       RLAM3M(:,:,:,:)=RLAM2M(:,:,:,:)
      END IF
      RLAM2M(:,:,:,:)=RLAMM(:,:,:,:)
      RLAMM(:,:,:,:)=RLAM0(:,:,:,:)
      RLAM0(:,:,:,:)=RLAMP(:,:,:,:)
      RLAMP(:,:,:,:)=0.D0
!SECOND ORDER LAGRANGE PREDICTION
!      ALLOCATE(DWORK(NB,NB,NKPT,NSPIN))
!      DWORK(:,:,:,:)=RLAM0(:,:,:,:)
!      RLAM0(:,:,:,:)=RLAMM(:,:,:,:)
!      RLAM0(:,:,:,:)=DWORK(:,:,:,:)
!      DEALLOCATE(DWORK)
      IF(ALLOCATED(HPSI)) DEALLOCATE(HPSI)
!
!     == DEALLOCATE HAMILTON MATRIX ====================================
      IF(ALLOCATED(HAMU)) DEALLOCATE(HAMU)
      IF(ALLOCATED(HAME)) DEALLOCATE(HAME)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES$DIAG
!     ******************************************************************
!     **  ROTATE TRAJECTORY                                           **
!     ******************************************************************
      USE WAVES_MODULE
      IMPLICIT NONE
      INTEGER(4),ALLOCATABLE :: ISPECIES(:)    !(NAT)
      INTEGER(4),ALLOCATABLE :: LMNX(:)        !(NSP)
      REAL(8)  ,ALLOCATABLE :: H0(:,:,:,:)    !(NB,NB,NKPT,NSPIN)   
      REAL(8)  ,ALLOCATABLE :: HP(:,:,:,:)    !(NB,NB,NKPT,NSPIN)   
      INTEGER(4)             :: ISP,ISPIN,IKPT
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES$DIAG')
      ALLOCATE(H0(NB,NB,NKPT,NSPIN))    
      ALLOCATE(HP(NB,NB,NKPT,NSPIN))    
!     
      CALL ERROR$MSG('ROUTINE NEEDS RECONSTRUCTION')
      CALL ERROR$STOP('WAVES$DIAG')
      CALL OCCUPATION$GET('H(0)',NB*NB*NKPT*NSPIN,H0)
      CALL OCCUPATION$GET('H(+)',NB*NB*NKPT*NSPIN,HP)
      CALL SETUP$NSPECIES(NSP)
      ALLOCATE(LMNX(NSP)) 
      DO ISP=1,NSP
        CALL SETUP$LMNX(ISP,LMNX(ISP))
      ENDDO
      ALLOCATE(ISPECIES(NAT)) 
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
           CALL WAVES_DIAG(NB,NGWLX,LMNXX &
     &               ,NSP,NAT,ISPECIES,LMNX,NB,NGWL(IKPT) &
     &               ,H0(1,1,IKPT,ISPIN),HP(1,1,IKPT,ISPIN) &
     &               ,FNL(1,1,1,IKPT,ISPIN) &
     &               ,CM(1,1,IKPT,ISPIN),C0(1,1,IKPT,ISPIN) &
     &               ,RLAMM(:,:,IKPT,ISPIN),RLAM0(:,:,IKPT,ISPIN))
        ENDDO
      ENDDO
      CALL OCCUPATION$SET('H(0)',NB*NB*NKPT*NSPIN,H0)
      CALL OCCUPATION$SET('H(+)',NB*NB*NKPT*NSPIN,HP)
!
      DEALLOCATE(H0)
      DEALLOCATE(HP)
      DEALLOCATE(ISPECIES)  
      DEALLOCATE(ISPECIES)  
                           CALL TRACE$POP
      RETURN
      END 
!
!     ==================================================================
!     ==  GENERATE G-VECTORS                                          ==
!     ==  INITIALIZE RANDOM                                           ==
!     ==  RANDOMIZE VELOCITY                                          ==
!     ==  FFT                                                         ==
!     ==  ORTHOGONALIZE                                               ==
!     ==  GRAMM SCHMITT                                               ==
!     ==  STRUCTURE FACTORS                                           ==
!     ==  KINETIC ENERGY                                              ==
!     ==  DENSITY                                                     ==
!     ==  ADDLOCAL                                                    ==
!     ==  PROJECTRIONS                                                ==
!     ==  ADDSPARABLE                                                 ==
!     ==  STOP                                                        ==
!     ==  READ                                                        ==
!     ==  WRITE                                                       ==
!     ==  PLOTRHO                                                     ==
!     ==================================================================
!
!     ..................................................................
      SUBROUTINE WAVES_BYPSIMASS(NGWX,NGW,NB,EMASS,EMASSCG2,G2,PSI)
!     ******************************************************************
!     **                                                              **
!     **  DIVIDE BY THE G-DEPENDENT MASS FOR THE WAVE FUNCTIONS       **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4), INTENT(IN)     :: NGW,NGWX
      INTEGER(4), INTENT(IN)     :: NB
      REAL(8),    INTENT(IN)     :: EMASS,EMASSCG2
      REAL(8),    INTENT(IN)     :: G2(NGWX)
      COMPLEX(8),INTENT(INOUT)  :: PSI(NGWX,NB)
      INTEGER(4)                 :: IB,IG       ! RUNNING VARIABLES
      REAL(8)                   :: SVAR1       ! AUXILARY VARIABLES
      SVAR1=1.D0 / EMASS 
      DO IB=1,NB
        DO IG=1,NGW
          PSI(IG,IB)=PSI(IG,IB)*SVAR1/(1.D0+EMASSCG2*G2(IG))
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_CORRECTGAMMA(TGAMMA,NGW,NB,PSI)
!     ******************************************************************
!     **                                                              **
!     **  SETS THE IMAGINARY PART OF THE WAVE FUNCTION FOR THE FIRST  **
!     **  G=VECTOR TO ZERO IF TGAMMA=.TRUE.                           **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      LOGICAL,   INTENT(IN)    :: TGAMMA
      INTEGER(4), INTENT(IN)    :: NGW         ! #(G-VECTORS)
      INTEGER(4), INTENT(IN)    :: NB          ! #(BANDS)
      COMPLEX(8),INTENT(INOUT) :: PSI(NGW,NB) ! PSEUDO WAVE FUNCTION
      INTEGER(4)                :: IB
      REAL(8)                  :: SVAR
!     ******************************************************************
      IF(TGAMMA) THEN
        DO IB=1,NB
          SVAR=DBLE(PSI(1,IB))
          PSI(1,IB)=CMPLX(SVAR,0.D0,KIND=8)
        ENDDO
      END IF
      RETURN
      END
!
!     ......................................................RANWAV......
      SUBROUTINE WAVES_INITIALIZERANDOM(NGW,NB,G2,PSI)
!     ******************************************************************
!     **                                                              **
!     **  CREATE RANDOM WAVE FUNCTIONS                                **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NB            ! #(BANDS)
      INTEGER(4),INTENT(IN)    :: NGW           ! #(PLANE WAVES),MAX
      REAL(8)   ,INTENT(IN)    :: G2(NGW)       ! G**2
      COMPLEX(8),INTENT(INOUT) :: PSI(NGW,NB)   ! PS-WAVE FUNCTION
      INTEGER(4)               :: IB,IG,IG1     !
      REAL(8)                  :: REC,RIMC
!     ******************************************************************
      DO IB=1,NB
        IG1=1
        IF(DABS(G2(1)).LT.1.D-6) THEN
          PSI(1,IB)=(1.D0,0.D0)
          IG1=2
        END IF
        DO IG=IG1,NGW
          CALL RANDOM_NUMBER(REC)
          CALL RANDOM_NUMBER(RIMC)
          PSI(IG,IB)=CMPLX(REC,RIMC,KIND=8)/G2(IG)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ......................................................RANWAV......
      SUBROUTINE WAVES_SETEQUAL(NGW,NB,COLD,CNEW)
!     ******************************************************************
!     **                                                              **
!     **  COPIES THE PS WAVE COLD INTO CNEW                           **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NB            ! # BANDS
      INTEGER(4),INTENT(IN)   :: NGW           ! # PLANE WAVES, MAX
      COMPLEX(8),INTENT(OUT)  :: CNEW(NGW,NB) ! WAVE FUNCTIONS
      COMPLEX(8),INTENT(IN)   :: COLD(NGW,NB) ! WAVE FUNCTIONS
      INTEGER(4)              :: IB,IG
!     ******************************************************************
      DO IB=1,NB
        DO IG=1,NGW
          CNEW(IG,IB)=COLD(IG,IB)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ......................................................RANWAV......
      SUBROUTINE WAVES_SETZERO(NGWX,NGW,NB,CNEW)
!     ******************************************************************
!     **                                                              **
!     **  SETS PSWAVE TO ZERO                                         **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NB             ! #(BANDS)
      INTEGER(4),INTENT(IN)  :: NGW,NGWX       ! #(G-VECTORS)
      COMPLEX(8),INTENT(OUT) :: CNEW(NGWX,NB)  ! WAVE FUNCTIONS
      INTEGER(4)             :: IB,IG          ! DO LOOP COUNTERS
!     ******************************************************************
      DO IB=1,NB
        DO IG=1,NGW
          CNEW(IG,IB)=(0.D0,0.D0)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ................................................WAVES_EIKR........
      SUBROUTINE WAVES_EIKR(NR1,NR2,NR3,NR1START,NR1END,IXK,EIKR)
!     ******************************************************************
!     **                                                              **
!     ** CALCULATES PHASE  FACTOR EXP( I*K*R ) ON THE REAL SPACE GRID **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NR1,NR2,NR3
      INTEGER(4),INTENT(IN) :: NR1START,NR1END
      INTEGER(4),INTENT(IN) :: IXK(3)
      COMPLEX(8),INTENT(OUT):: EIKR((NR1END-NR1START+1)*NR2*NR3)
      COMPLEX(8)            :: DEIK1,DEIK2,DEIK3
      COMPLEX(8)            :: EIKR1,EIKR2,EIKR3
      REAL(8)               :: PI
      REAL(8)               :: PHAS1,PHAS2,PHAS3
      INTEGER(4)            :: IRR,IR1,IR2,IR3
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      PHAS1=2.D0*PI * DBLE(IXK(1))/2.D0 / DBLE(NR1)
      PHAS2=2.D0*PI * DBLE(IXK(2))/2.D0 / DBLE(NR2)
      PHAS3=2.D0*PI * DBLE(IXK(3))/2.D0 / DBLE(NR3)
      DEIK1=CMPLX(DCOS(PHAS1),DSIN(PHAS1),KIND=8)
      DEIK2=CMPLX(DCOS(PHAS2),DSIN(PHAS2),KIND=8)
      DEIK3=CMPLX(DCOS(PHAS3),DSIN(PHAS3),KIND=8)
      EIKR3=(1.D0,0.D0)
      IRR=0
      DO IR3=1,NR3
        EIKR2=EIKR3
        DO IR2=1,NR2
          EIKR1=EIKR2*DEIK1**(NR1START-1)
          DO IR1=NR1START,NR1END
            IRR=IRR+1
            EIKR(IRR)=EIKR1
            EIKR1=EIKR1*DEIK1
          ENDDO
          EIKR2=EIKR2*DEIK2
        ENDDO
        EIKR3=EIKR3*DEIK3
      ENDDO
      RETURN
      END
!
!     .....................................................RHOOFR.......
      SUBROUTINE WAVES_DENSITY(NGWX,NRSTORE,NR1,NR2,NR3,NB,F,WKPT &
     &                ,PSIOFG,RHO,PSIOFR)
!     ******************************************************************
!     **                                                              **
!     **  THE  ELECTRON DENSITY RHOE IN REAL SPACE                    **
!     **                                                              **
!     **  INPUT:                                                      **
!     **    (NR1,NR2,NR3) NUMBER OF POINTS OF THE FOURIER GRID        **
!     **                  FOR THE THREE DIRECTIONS                    **
!     **    NGW           NUMBER OF G-VEKTORS FOR TEH WAVE FUNCTION   **
!     **                  DEPENDENT ON THE K-POINT                    **
!     **    N             NUMBER OF BANDS                             **
!     **    F             OCCUPATIONS                                 **
!     **    WKPT          GEOMETRICAL WEIGHT OF ONE K-POINT           **
!     **    C              WAVE FUNCTIONS IN G-SPACE                  **
!     **                                                              **
!     **  OUTPUT:                                                     **
!     **    RHOE          CHARGE DENSITY IN REAL SPACE                **
!     **    PSIOFR        WAVE FUNCTIONS IN REAL SPACE (REAL(4)!!)     **
!     **                                                              **
!     **  SAVE: NONE                                                  **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NRSTORE     ! MAX # PLANE WAVES
      INTEGER(4),INTENT(IN)  :: NGWX        ! MAX # PLANE WAVES
      INTEGER(4),INTENT(IN)  :: NB          ! # BANDS
      INTEGER(4),INTENT(IN)  :: NR1,NR2,NR3 ! # R-SPACE POINTS
      REAL(8)   ,INTENT(IN)  :: F(NB)       ! OCCUPATIONS
      REAL(8)   ,INTENT(IN)  :: WKPT        ! GEOMETRICAL K-POINT WEIGHT
      COMPLEX(8),INTENT(IN)  :: PSIOFG(NGWX,NB)
      REAL(4)   ,INTENT(OUT) :: PSIOFR(NRSTORE,NB)
      REAL(8)   ,INTENT(OUT) :: RHO(NR1*NR2*NR3) ! DENSITY IN R-SPACE
      LOGICAL                :: TSTORE
      INTEGER(4)             :: NNR          ! =NR1*NR2*NR3
      INTEGER(4)             :: IBPAIR,I,IR,IB
      INTEGER(4)             :: NFFT   
      INTEGER(4)             :: IB1,IB2
      REAL(8)                :: FAC,SVAR
      REAL(8)   ,ALLOCATABLE :: DWORK(:,:) !(NR1*NR2*NR3,2) !(IP)
!     ******************************************************************
      NNR=NR1*NR2*NR3
      TSTORE=(NRSTORE.GT.1)
      ALLOCATE(DWORK(NR1*NR2*NR3,2)) 
!
!     ==================================================================
!     ==  LOOP OVER PAIRS OF WAVE FUNCTIONS                           ==
!     ==================================================================
      DO IBPAIR=1,NB,2
        NFFT=MIN(2,NB-IBPAIR+1)
        IB1=IBPAIR
        IB2=IB1+NFFT-1
!
!       ================================================================
!       ==  FOURIERTRANSFORM TO REAL SPACE                            ==
!       ================================================================
        CALL PLANEWAVE$SUPFFT('GTOR',' ',NFFT,NNR,NR1,NR2,NR3,NGWX &
     &              ,DWORK(1,1),DWORK(1,2),PSIOFG(1,IB1),PSIOFG(1,IB2))
!
!       ================================================================
!       ==  CALCULATE CHARGE DENSITY                                  ==
!       ================================================================
        I=0
        DO IB=IB1,IB2
          I=I+1
          FAC=F(IB)*WKPT
          IF(TSTORE) THEN
!           __STORE REAL SPACE WAVE FUNCTIONS___________________________ 
            DO IR=1,NNR
              SVAR=DWORK(IR,I)
              RHO(IR)=RHO(IR)+FAC*SVAR**2
              PSIOFR(IR,IB)=REAL(SVAR,KIND=4)
            ENDDO
          ELSE
!           __OR DONT___________________________________________________
            DO IR=1,NNR
              SVAR=DWORK(IR,I)
              RHO(IR)=RHO(IR)+FAC*SVAR**2
            ENDDO
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(DWORK)
      RETURN
      END
!
!     .....................................................MOMNTS.......
      SUBROUTINE WAVES_1CDENSITYMATRIX(NAT,LMNXX,LMNX,NB,IAT &
     &          ,WKPT,F,FNL,DENMAT,OP)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE ONE-CENTER DENSITY MATRIX FOR A GIVEN        **
!     **  SET OF OCCUPATIONS                                          **
!     **      N^1=<R|PHI_I>DENMAT(I,J)<PHI(J)|R>                      **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),  INTENT(IN)   :: IAT   ! ATOM INDEX
      INTEGER(4),  INTENT(IN)   :: NAT   ! #(ATOMS)
      INTEGER(4),  INTENT(IN)   :: LMNX,LMNXX 
      INTEGER(4),  INTENT(IN)   :: NB    ! #(BANDS)
      CHARACTER*1,INTENT(IN)   :: OP    ! SWITCH: RESET DENMAT OR ADD RESULT
      REAL(8),     INTENT(IN)   :: WKPT  ! GEOMETRICAL K-POINT WEIGHT
      REAL(8),     INTENT(IN)   :: F(NB) ! OCCUPATIONS
      REAL(8),     INTENT(IN)   :: FNL(NAT,NB,LMNX)    ! PROJECTIONS
      REAL(8),     INTENT(INOUT):: DENMAT(LMNXX,LMNXX) ! 1C DENSITY MATRIX
      REAL(8)                  :: SVAR
      INTEGER(4)                :: LMN1,LMN2,IB
!     ******************************************************************
      IF(OP.NE.' '.AND.OP.NE.'+') THEN
        CALL ERROR$STOP('WAVES_1CDENSITYMATRIX')
      END IF
      IF(LMNXX.LT.LMNX) THEN
        CALL ERROR$STOP('WAVES_1CDENSITYMATRIX')
      END IF
!
!     ==================================================================
!     ==   RESET 1C DENSITY MATRIX                                    ==
!     ==================================================================
      IF(OP.EQ.' ') THEN
        DO LMN1=1,LMNXX
          DO LMN2=1,LMNXX
            DENMAT(LMN1,LMN2)=0.D0
          ENDDO
        ENDDO
      END IF
!
!     ==================================================================
!     ==   CALCULATE ONE CENTER CHARGE DENSITY                        ==
!     ==================================================================
      DO LMN1=1,LMNX
        DO LMN2=1,LMNX
          SVAR=0.D0
          DO IB=1,NB
            SVAR=SVAR+F(IB)*WKPT*FNL(IAT,IB,LMN1)*FNL(IAT,IB,LMN2)
          ENDDO
          DENMAT(LMN1,LMN2)=DENMAT(LMN1,LMN2)+SVAR
        ENDDO
      ENDDO
      RETURN
      END
!
!     .....................................................MOMNTS.......
      SUBROUTINE WAVES_1CHAMILTONIAN(NAT,LMNXX,LMNX,NB,IAT &
     &          ,FNL,HNL,DATH)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE ONE-CENTER DENSITY MATRIX FOR A GIVEN        **
!     **  SET OF OCCUPATIONS                                          **
!     **      N^1=<R|PHI_I>DENMAT(I,J)<PHI(J)|R>                      **
!     **                                                              **
!     **  INPUT:                                                      **
!     **    NB       NUMBER OF BANDS                                  **
!     **    LOX      MAIN ANGULAR MOMENTUM OF A RADIAL PARTIAL WAVE   **
!     **    F        OCCUPATIONS                                      **
!     **    FNL      PROJECTIONS                                      **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),  INTENT(IN)   :: IAT   ! ATOM INDEX
      INTEGER(4),  INTENT(IN)   :: NAT   ! #(ATOMS)
      INTEGER(4),  INTENT(IN)   :: LMNX,LMNXX 
      INTEGER(4),  INTENT(IN)   :: NB    ! #(BANDS)
      REAL(8),     INTENT(IN)   :: DATH(LMNXX,LMNXX) ! 1C DIFF-HAMILTONIAN
      REAL(8),     INTENT(IN)   :: FNL(NAT,NB,LMNX)    ! <PS-P|PS-PSI>
      REAL(8),     INTENT(OUT)  :: HNL(NAT,NB,LMNX)    ! 
      INTEGER(4)                :: LMN,LMN1,LMN2,IB
!     ******************************************************************
!
      DO LMN=1,LMNX
        DO IB=1,NB
          HNL(IAT,IB,LMN)=0.D0
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==   CALCULATE PREFACTOR OF PROJECTORS IN PS-H|PS-PSI>          ==
!     ==================================================================
      DO LMN1=1,LMNX
        DO LMN2=1,LMNX
          DO IB=1,NB
            HNL(IAT,IB,LMN1)=HNL(IAT,IB,LMN1) &
     &                      +DATH(LMN1,LMN2)*FNL(IAT,IB,LMN2)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
@PROCESS NOEXTCHK
!     .....................................................NLSMX........
      SUBROUTINE WAVES_PROJECTIONS(LNXX,LMNXX,NGWX &
     &                 ,NSP,NAT,ISPECIES,NB,LNX,LOX,CELLVOL &
     &                 ,NGW,HSG,GSK &
     &                 ,WNL,EIGR,PSIG,IGRAD,FNL)
!     ******************************************************************
!     **                                                              **
!     **  THE ARRAY FNL WHICH IS USED IN NLRH                         **
!     **           FNL = <PSI|PRO>                                    **
!     **  OR THE GRADIENT ( FOR KK.NE.0 )                             **
!     **           FNL=<PSI|GRAD(PRO)>                                **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      INTEGER(4),INTENT(IN) :: IGRAD    ! SWITCH: GRADIENT CALCULATION OR NOT
      INTEGER(4),INTENT(IN) :: NAT      ! #(ATOMS)
      INTEGER(4),INTENT(IN) :: NSP      ! #(ATOM TYPES)
      INTEGER(4),INTENT(IN) :: LNXX
      INTEGER(4),INTENT(IN) :: NGW,NGWX ! #(PLANE WAVES), MAX
      INTEGER(4),INTENT(IN) :: NB       ! #(BANDS)
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: LNX(NSP)
      INTEGER(4),INTENT(IN) :: LOX(LNXX,NSP)
      INTEGER(4),INTENT(IN) :: ISPECIES(NAT)
      COMPLEX(8),INTENT(IN) :: PSIG(NGWX,NB)       ! PS WAVE
      COMPLEX(8),INTENT(IN) :: EIGR(NGWX,NAT)      ! STRUCTURE FACTORS
      REAL(8)   ,INTENT(IN) :: WNL(NGWX,LMNXX,NSP) ! FORM FACTORS
      REAL(8)   ,INTENT(OUT):: FNL(NAT,NB,LMNXX,1+2*IGRAD)
      REAL(8)   ,INTENT(IN) :: GSK(3,NGWX)
      REAL(8)   ,INTENT(IN) :: CELLVOL               ! R-SPACE UNIT CELL VOLUME
      REAL(8)   ,INTENT(IN) :: HSG(NGWX)           ! G**2
      COMPLEX(8)            :: CSVAR
      COMPLEX(8),ALLOCATABLE:: AUXC(:,:)
      REAL(8)   ,ALLOCATABLE:: RNL(:)
      INTEGER(4)            :: NPROS
      INTEGER(4)            :: KKX
      INTEGER(4)            :: IB,IG
      INTEGER(4)            :: KK,LMN,IAT,I,NPRO,LN,ISP,L,IM
      INTEGER(4)            :: KK1,IAT1,ISP1,LMN1,LN1,L1,IM1,IPRO1
      INTEGER(4)            :: IPRO,IDPRO,ISWEEP,IPROLAST
      INTEGER(4)            :: IFIRST,ILAST,ICOUNT
      INTEGER(4)            :: II
!     ******************************************************************
      NPROS=NB
      KKX=1+2*IGRAD
      IF(IGRAD.NE.0.AND.IGRAD.NE.1) THEN
        CALL ERROR$MSG('IGRAD IN NLSMX MUST BE 0 OR 1')
        CALL ERROR$STOP('WAVES_PROJECTIONS')
      END IF
!     ==================================================================
!     ==  INITIALIZE ARRAYS                                           ==
!     ==================================================================
      DO I=1,KKX
        DO LMN=1,LMNXX
          DO IB=1,NB
            DO IAT=1,NAT
              FNL(IAT,IB,LMN,I)=0.D0
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      NPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          NPRO=NPRO+2*L+1
        ENDDO
      ENDDO
      ALLOCATE(AUXC(NGWX,NB))
      ALLOCATE(RNL(NPROS*NB))
!     PRINT*,'NPRO,NPROS,KKX,NB ',NPRO,NPROS,KKX,NB
!
!     ==================================================================
!     ==  START THE LOOP                                              ==
!     ==================================================================
      IPRO=0
      IDPRO=0
      ISWEEP=0
      IPROLAST=MIN(NPROS,NPRO*KKX-IDPRO)
      DO KK=1,KKX
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          LMN=0
          DO LN=1,LNX(ISP)
            L=LOX(LN,ISP)
            DO IM=1,2*L+1
              LMN=LMN+1
              IPRO=IPRO+1
!
!             ==========================================================
!             ==  COMPOSE PROJECTOR FUNCTIONS                         ==
!             ==========================================================
!             ==  MULTIPLY PROJECTOR WITH STRUCTURE FACTOR ===========
!             ==  THE FACTORS 2 AND .5 ARE TO INCLUDE G<0  ===========
              CSVAR=2.D0*(-CI)**L
              DO IG=1,NGW
                AUXC(IG,IPRO)=WNL(IG,LMN,ISP)*(CSVAR*EIGR(IG,IAT))
              ENDDO
              IF(DABS(HSG(1)).LT.1.D-6) THEN
                AUXC(1,IPRO)=.5D0*AUXC(1,IPRO)
              END IF
!             ==  MULTIPLY WITH GRADIENT =============================
              IF(IGRAD.EQ.1) THEN
                DO IG=1,NGW
                  AUXC(IG,IPRO)=-GSK(KK,IG)*CI*AUXC(IG,IPRO)
                ENDDO
              END IF
!             PRINT*,'IPRO ',IPRO,LMN,IAT,KK,IDPRO+IPRO
!
!             ==========================================================
!             ==  MULTIPLY WITH WAVE FUNCTIONS                        ==
!             ==========================================================
              IF(IPRO.GE.IPROLAST) THEN
!               == EVALUATE PROJECTIONS ================================
                              CALL TIMING$CLOCKON('WAVES_PROJ.DGEMUL')
!               CALL DGEMUL(TRANSFER(AUXC,RNL,2*NGWX*IPRO),2*NGWX,'T' &
!    &                     ,TRANSFER(PSIG,RNL,2*NGWX*NB),2*NGWX,'N' &
!    &                     ,RNL,IPRO,IPRO,2*NGW,NB)
! JUST ANOTHER VERSION OF THE DGEMUL CALL FOR TEST
!               DO I=1,IPRO
!                 DO IB=1,NB
!                   II=I+IPRO*(IB-1)
!                   RNL(II)=DOT_PRODUCT(AUXC(:,I),PSIG(:,IB))
!                 ENDDO
!               ENDDO
! THE PREVIOUS VERSION WITH THE PROBLEM OF INCORRECT INTERFACE
                CALL DGEMUL(AUXC,2*NGWX,'T' &
     &                     ,PSIG,2*NGWX,'N' &
     &                     ,RNL,IPRO,IPRO,2*NGW,NB)
                              CALL TIMING$CLOCKOFF('WAVES_PROJ.DGEMUL')
!               == MAP PROJECTIONS BACK ================================           
                IFIRST=IDPRO+1
                ILAST=IDPRO+IPRO
                ICOUNT=0
                DO KK1=1,KKX
                  DO IAT1=1,NAT
                    ISP1=ISPECIES(IAT1)
                    LMN1=0
                    DO LN1=1,LNX(ISP1)
                      L1=LOX(LN1,ISP1)
                      DO IM1=1,2*L1+1
                        LMN1=LMN1+1
                        ICOUNT=ICOUNT+1
                        IF(ICOUNT.GE.IFIRST.AND.ICOUNT.LE.ILAST) THEN
                          IPRO1=ICOUNT-IFIRST+1
!                         PRINT*,'MAP ',IPRO1,LMN1,IAT1,KK1,ICOUNT
                          DO IB=1,NB
                            II=IPRO1+IPRO*(IB-1)
                            FNL(IAT1,IB,LMN1,KK1)=RNL(II)*CELLVOL
                          ENDDO
                        END IF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
!               ==  RESET POINTERS =====================================
                IDPRO=IDPRO+IPRO
                ISWEEP=ISWEEP+1
                IPRO=0
                IPROLAST=MIN(NPROS,NPRO*KKX-IDPRO)
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(RNL)
      DEALLOCATE(AUXC)
!
!     ==================================================================
!     ==  OUTPUT OF OCCUPATIONS                                       ==
!     ==================================================================
      IF(TPR) THEN
        IF(IGRAD.EQ.1) THEN
          CALL WAVES_PRFNL(NAT,NB,LMNXX,FNL(1,1,1,1),'<DPRO/DX|PSI>')
          CALL WAVES_PRFNL(NAT,NB,LMNXX,FNL(1,1,1,2),'<DPRO/DY|PSI>')
          CALL WAVES_PRFNL(NAT,NB,LMNXX,FNL(1,1,1,3),'<DPRO/DZ|PSI>')
        ELSE IF(IGRAD.EQ.0) THEN
          CALL WAVES_PRFNL(NAT,NB,LMNXX,FNL(1,1,1,1),'<PRO|PSI>')
        ENDIF
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_PROJECTIONCOMBINE(NDIM,NSPIN,NKPT,NB,NAT,ISPECIES &
     &                                  ,LMNXX,NSP,LMNX,FNL)
!     ******************************************************************
!     **                                                              **
!     **  ADDS THE PROJECTIONS OVER ALL TASKS, BY COMPRESSING AND     **
!     **  DECOMPRESSING THE ARRAY                                     **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NKPT          ! #(K-POINTS)
      INTEGER(4),INTENT(IN)   :: NSPIN         ! #(SPINS)
      INTEGER(4),INTENT(IN)   :: NB            ! #(BANDS) 
      INTEGER(4),INTENT(IN)   :: NAT           ! #(ATOMS)
      INTEGER(4),INTENT(IN)   :: NSP           ! #(SPECIES)
      INTEGER(4),INTENT(IN)   :: NDIM          ! 
      INTEGER(4),INTENT(IN)   :: LMNXX         ! MAX #(PROJECTIONS;L,M,N)
      INTEGER(4),INTENT(IN)   :: ISPECIES(NAT) ! ATOM TYPE INDEX ARRAY
      INTEGER(4),INTENT(IN)   :: LMNX(NSP)     ! #(PROJECTIONS)
      REAL(8)   ,INTENT(INOUT):: FNL(NAT,NB,LMNXX,NDIM,NKPT,NSPIN)
      INTEGER(4)              :: ICOUNT,NWORK
      INTEGER(4)              :: ISPIN,IKPT,IAT,ISP,LMN,IB,IDIM
      REAL(8)   ,ALLOCATABLE  :: WORK(:)
!     ******************************************************************
!
!     ==================================================================
!     == EVALUATE SIZE OF WORK ARRAY                                  ==
!     ==================================================================
      ICOUNT=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        ICOUNT=ICOUNT+LMNX(ISP)
      ENDDO
      NWORK=ICOUNT*NDIM*NSPIN*NKPT*NB
      ALLOCATE(WORK(NWORK))
!
!     ==================================================================
!     == MAP ONTO WORK ARRAY                                          ==
!     ==================================================================
      ICOUNT=0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IDIM=1,NDIM
            DO IAT=1,NAT
              ISP=ISPECIES(IAT)
              DO LMN=1,LMNX(ISP)
                DO IB=1,NB
                  ICOUNT=ICOUNT+1
                  WORK(ICOUNT)=FNL(IAT,IB,LMN,IDIM,IKPT,ISPIN)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     == ADD WORK ARRAY                                               ==
!     ==================================================================
      CALL MPE$COMBINE('+',WORK)
!
!     ==================================================================
!     == MAP BACK FROM WORK ARRAY                                     ==
!     ==================================================================
      ICOUNT=0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IDIM=1,NDIM
            DO IAT=1,NAT
              ISP=ISPECIES(IAT)
              DO LMN=1,LMNX(ISP)
                DO IB=1,NB
                  ICOUNT=ICOUNT+1
                  FNL(IAT,IB,LMN,IDIM,IKPT,ISPIN)=WORK(ICOUNT)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(WORK)
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE WAVES_PRFNL(NAT,NB,LMNXX,FNL,TEXT)
!     ******************************************************************
!     **                                                              **
!     ** PRINT PROJECTIONS     <P|PSI>                                **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      CHARACTER*(*),INTENT(IN) :: TEXT
      REAL(8)      ,INTENT(IN) :: FNL(NAT,NB,LMNXX)  !<PS-P|PS-PSI>
      INTEGER(4)   ,INTENT(IN) :: NAT                ! #(ATOMS)
      INTEGER(4)   ,INTENT(IN) :: NB                 ! #(BANDS)
      INTEGER(4)   ,INTENT(IN) :: LMNXX  
      INTEGER(4)               :: IAT,IB,LMN
!     ******************************************************************
      WRITE(*,FMT='("PRFNL: ",A)')TEXT
      DO IAT=1,NAT
        WRITE(*,FMT='("PROJECTIONS FOR ATOM ",I5)')IAT
        DO IB=1,NB
          WRITE(*,FMT='(I5,10F10.5)') &
     &          IB,(FNL(IAT,IB,LMN),LMN=1,LMNXX)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ......................................................EKINFE......
      SUBROUTINE WAVES_WAVEKINETIC(NGWX,NGW,NB &
     &           ,TGAMMA,F,WKPT,CELLVOL,C0,CM &
     &           ,DELT,EMASS,EMASSCG2,G2,EKINC)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATE KINETIC ENERGY OF THE FICTITIOUS ELECTRONIC       **
!     **  VARIABLES BETWEEN C(0) AND C(-)                             **
!     **  REMARK : THIS FORMULA IS NOT CONSISTENT WITH THE VERLET     **
!     **  ALGORITHM                                                   **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NGWX
      INTEGER(4),INTENT(IN) :: NGW              ! NUMBER OF G-VECTORS
      INTEGER(4),INTENT(IN) :: NB               ! NUMBER OF STATES
      LOGICAL   ,INTENT(IN) :: TGAMMA           !
      REAL(8)   ,INTENT(IN) :: F(NB)            ! OCCUPATIONS
      REAL(8)   ,INTENT(IN) :: WKPT             ! K-POINT  WEIGHT
      REAL(8)   ,INTENT(IN) :: CELLVOL            ! UNIT-CELL VOLUME
      REAL(8)   ,INTENT(IN) :: DELT             ! TIME STEP
      REAL(8)   ,INTENT(IN) :: EMASS            ! WAVE FUNCTION MASS
      REAL(8)   ,INTENT(IN) :: EMASSCG2         
      REAL(8)   ,INTENT(IN) :: G2(NGWX)         ! G**2
      COMPLEX(8),INTENT(IN) :: C0(NGWX,NB)      ! WAVE FUNCTION (0)
      COMPLEX(8),INTENT(IN) :: CM(NGWX,NB)      ! WAVE FUNCTION (-)
      REAL(8)   ,INTENT(OUT):: EKINC            ! KINETIC ENERGY
      REAL(8)               :: SVAR,DEKINC,SVAR1
      COMPLEX(8)            :: SPEED            
      INTEGER(4)            :: IB,IG            ! RUNNING INDICES
!     ******************************************************************
!
!     ==================================================================
!     ==  PLANE WAVE CONTRIBUTION                                     ==
!     ==================================================================
      SVAR=2.D0*EMASS/DELT**2
      EKINC=0.D0
      DO IB=1,NB
        DEKINC=0.D0
        IF(TGAMMA) THEN
          SPEED=C0(1,IB)-CM(1,IB)
          DEKINC=-0.5D0*DBLE(SPEED*CONJG(SPEED))
        END IF
        DO IG=1,NGW
          SVAR1=(1.D0+EMASSCG2*G2(IG))
          SPEED=C0(IG,IB)-CM(IG,IB)
          DEKINC=DEKINC+SVAR1*DBLE(SPEED*CONJG(SPEED))
        ENDDO
        EKINC=EKINC+DEKINC*F(IB)*WKPT*CELLVOL*SVAR
      ENDDO
      RETURN
      END
!
!     ................................................WAVES_SWITCH......
      SUBROUTINE WAVES_SWITCH(NGWX,NGW,NB,PSI1,PSI2)
!     ******************************************************************
!     **                                                              **
!     **  INTERCHANGES PSI1 AND PSI2                                  **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: NB       ! #(BANDS)
      INTEGER(4) ,INTENT(IN)    :: NGW,NGWX ! #(PLANE WAVES),MAX
      COMPLEX(8),INTENT(INOUT)  :: PSI1(NGWX,NB)
      COMPLEX(8),INTENT(INOUT)  :: PSI2(NGWX,NB)
      COMPLEX(8)                :: SPEED
      INTEGER(4)                :: IB,IG
!     ******************************************************************
      DO IB=1,NB
        DO IG=1,NGW
          SPEED=PSI1(IG,IB)
          PSI1(IG,IB)=PSI2(IG,IB)
          PSI2(IG,IB)=SPEED
        ENDDO
      ENDDO
      RETURN
      END
!
!     .....................................................DFORCE.......
      SUBROUTINE WAVES_HPSI(NGWX,NRSTORE,NB,NGW,NR1,NR2,NR3,HSG &
     &                  ,PSI,V,HPSI,STORE)
!     ******************************************************************
!     **                                                              **
!     **  EVALUATES H*PSI WITHOUT THE AUGMENTATION PART               **
!     **  |HPSI>=(-0.5*NABLA**2+PS-V)|PSPSI(0)>                       **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NGW,NGWX    ! #(PLANE WAVES),MAX
      INTEGER(4),INTENT(IN)  :: NB          ! #(BANDS)
      INTEGER(4),INTENT(IN)  :: NRSTORE     ! #(R-SPACE GRID POINTS) OR 1
      INTEGER(4),INTENT(IN)  :: NR1,NR2,NR3 ! #(R-SPACE POINTS) / EACH LATTICE VECTOR
      COMPLEX(8),INTENT(IN)  :: PSI(NGWX,NB)      ! PSPSI(0)
      REAL(8)   ,INTENT(IN)  :: V(NR1*NR2*NR3)    ! PSPOT
      REAL(8)   ,INTENT(IN)  :: HSG(NGW)          ! G**2
      REAL(4)   ,INTENT(IN)  :: STORE(NRSTORE,NB) ! REAL SPACE PSPSI(0)
      COMPLEX(8),INTENT(OUT) :: HPSI(NGWX,NB)     ! PSH|PSPSI>
      LOGICAL                :: TSTORE
      INTEGER(4)             :: NNR
      INTEGER(4)             :: IB,IB1,IB2,IBPAIR,NFFT
      INTEGER(4)             :: I,IR,IG
      REAL(8)                :: SVAR
      REAL(8)   ,ALLOCATABLE :: DWORK(:,:)     !(NR1*NR2*NR3,2)
!     ******************************************************************
      NNR=NR1*NR2*NR3
      TSTORE=(NRSTORE.GT.1)
      ALLOCATE(DWORK(NR1*NR2*NR3,2))
!
!     ==================================================================
!     ==  LOOP OVER PAIRS OF WAVE FUNCTIONS                           ==
!     ==================================================================
      DO IBPAIR=1,NB,2
        IB1=IBPAIR
        IB2=MIN(IB1+1,NB)
        NFFT=IB2-IB1+1
!
!       ================================================================
!       ==  MULTIPLY WAVE FUNCTIONS WITH THE POTENTIAL                ==
!       ================================================================
        IF(TSTORE) THEN
!         __USE STORED REAL SPACE WAVE FUNCTIONS IF THEY EXIST__________
          IF(NRSTORE*4.LT.NGWX*16) THEN
            CALL ERROR$MSG('OVERLAYING OF STORE AND HPSI FAILED')
            CALL ERROR$MSG('4*NRSTORE<16*NGWX')
            CALL ERROR$STOP('WAVES_HPSI')
          END IF
          I=0
          DO IB=IB1,IB2
            I=I+1
            DO IR=1,NNR
              SVAR=DBLE(STORE(IR,IB))
              DWORK(IR,I)=V(IR)*SVAR
            ENDDO
          ENDDO
        ELSE 
!         __AND RECALCULATE THEM IF THEY DONT___________________________
          CALL PLANEWAVE$SUPFFT('GTOR',' ',NFFT,NNR,NR1,NR2,NR3,NGW &
     &              ,DWORK(1,1),DWORK(1,2),PSI(1,IB1),PSI(1,IB2))
          DO I=1,NFFT
            DO IR=1,NNR
              DWORK(IR,I)=V(IR)*DWORK(IR,I)
            ENDDO
          ENDDO          
        END IF
!
!       ================================================================
!       ==  FOURIER TRANSFORM TO RECIPROCAL SPACE                     ==
!       ================================================================
        CALL PLANEWAVE$SUPFFT('RTOG',' ',NFFT,NNR,NR1,NR2,NR3,NGW &
     &              ,DWORK(1,1),DWORK(1,2),HPSI(1,IB1),HPSI(1,IB2))
      ENDDO
      DEALLOCATE(DWORK)
!
!     ==================================================================
!     ==  ADD KINETIC ENERGY CONTRIBUTION                             ==
!     ==================================================================
      DO IB=1,NB
        DO IG=1,NGW
          HPSI(IG,IB)=HPSI(IG,IB)+0.5D0*HSG(IG)*PSI(IG,IB)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..............................................WAVES_ADDPRO........
      SUBROUTINE WAVES_ADDPRO(NGWX,LMNXX,LNXX &
     &                       ,NSP,NAT,ISPECIES,NB,NGW,LNX,LOX &
     &                       ,PSIG,PROG,XNL,EIGR)
!     ******************************************************************
!     **                                                              **
!     **  ADDS PROJECTORS IN G-SPACE TO A WAVE FUNCTION               **
!     **                                                              **
!     **    PSIG(IG,IB) = PSIG(IG,IB)+PROG(IG,RLMN)*(-I)**L*XNL(RLMN) **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
      IMPLICIT NONE
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)
      INTEGER(4) ,INTENT(IN)    :: NGW,NGWX ! NUMBER OF PLANE WAVES,MAX
      INTEGER(4) ,INTENT(IN)    :: NB       ! NUMBER OF BANDS
      INTEGER(4) ,INTENT(IN)    :: LMNXX
      INTEGER(4) ,INTENT(IN)    :: NSP      ! NUMBER OF ATOM TYPES
      INTEGER(4) ,INTENT(IN)    :: NAT      ! NUMBER OF ATOMS
      INTEGER(4) ,INTENT(IN)    :: LNXX     ! MAX #(DIFFERENT PROJECTOR FUNCTIONS)
      INTEGER(4) ,INTENT(IN)    :: LNX(NSP) ! #(DIFFERENT PROJECTOR FUNCTIONS)
      INTEGER(4) ,INTENT(IN)    :: LOX(LNXX,NSP) ! #(DIFFERENT PROJECTOR FUNCTIONS)
      INTEGER(4) ,INTENT(IN)    :: ISPECIES(NAT)        ! ATOM TYPE INDEX ARRAY
      REAL(8)   ,INTENT(IN)    :: PROG(NGWX,LMNXX,NSP) ! <PS-P|G>
      COMPLEX(8),INTENT(IN)    :: EIGR(NGWX,NAT)       ! STRUCTURE FACTOR
      REAL(8)   ,INTENT(IN)    :: XNL(NAT,NB,LMNXX)
      COMPLEX(8),INTENT(INOUT) :: PSIG(NGWX,NB)        ! PS-PSI
      REAL(8)   ,ALLOCATABLE   :: AUXR(:),AUXI(:)      ! WORK ARRAYS
      COMPLEX(8)               :: CILIN
      COMPLEX(8)               :: CSVAR
      INTEGER(4)                :: IB,IAT,IG,LMN,LN,IM  ! RUNNING VARIABLES
      INTEGER(4)                :: L     ! MAIN ANGULAR MOMENTUM
      INTEGER(4)                :: ISP   ! ATOM TYPE INDEX
      INTEGER(4)                :: ISVAR ! AUXILIARY VARIABLE
      REAL(8)                  :: FAC,SVAR
!     ******************************************************************
      ALLOCATE(AUXR(NGWX))
      ALLOCATE(AUXI(NGWX))
      DO IB=1,NB
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          DO IG=1,NGW
            AUXR(IG)=0.D0
            AUXI(IG)=0.D0
          ENDDO
          LMN=0
!         ==  WNL=SUM{G|PRO(G,LMN)*(-I)**L*XNL(LMN)}
          DO LN=1,LNX(ISP)
            L=LOX(LN,ISP)
            CILIN=(-CI)**L
            ISVAR=L/2
            IF(L.EQ.2*ISVAR) THEN
              FAC=DBLE(CILIN)
              DO IM=1,2*L+1
                LMN=LMN+1
                SVAR=FAC*XNL(IAT,IB,LMN)
!               ==  AUXR(I)=AUXR(I)+PROG(I)*SVAR
                IF(DABS(SVAR).GT.1.D-6) THEN
                  CALL DAXPY(NGW,SVAR,PROG(1,LMN,ISP),1,AUXR,1)
                END IF
              ENDDO
            ELSE
              FAC=AIMAG(CILIN)
              DO IM=1,2*L+1
                LMN=LMN+1
                SVAR=FAC*XNL(IAT,IB,LMN)
!               ==  AUXI(I)=AUXI(I)+PROG(I)*SVAR
                IF(DABS(SVAR).GT.1.D-6) THEN
                  CALL DAXPY(NGW,SVAR,PROG(1,LMN,ISP),1,AUXI,1)
                END IF
              ENDDO
            END IF
          ENDDO
!         ==  REMARK: CAUX IS EQUIVALENT TO AUX
          DO IG=1,NGW
            CSVAR=CMPLX(AUXR(IG),AUXI(IG),KIND=8)
            PSIG(IG,IB)=PSIG(IG,IB)+EIGR(IG,IAT)*CSVAR
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(AUXR)
      DEALLOCATE(AUXI)
      RETURN
      END
!
!     .....................................................ORTHO........
      SUBROUTINE WAVES_ORTHO(NX,LMNXX,LNXX,NGWX &
     &                ,NSP,NAT,ISPECIES,LNX,LMNX,LOX,CELLVOL &
     &                ,TGAMMA,NB,F,NGW,CHI,PSI,ONL,FNL,LAMBDA &
     &                ,DOVER,DELT,EMASS,ANNEE,TSAFEORTHO)
!     ******************************************************************
!     **                                                              **
!     **  IMPOSES THE ORTHOGONALITY CONSTRAINT ONTO THE ELECTRONS     **
!     **                                                              **
!     **  NEW VERSION WITH DIAGONALIZATION FOR RHO                    **
!     **                                                              **
!     **  THE METHOD IS DESCRIBED IN :                                **
!     **    R.CAR AND M.PARRINELLO, IN "SIMPLE MOLECULAR SYSTEMS      **
!     **    AT VERY HIGH DENSITY", PAGE 455                           **
!     **    ED. A.POLIAN, PLOUBEYRE AND N.BOCCARA                     **
!     **    (PLENUM PUBLISHING CORPORATION,1989)                      **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
      USE MPE_MODULE
      IMPLICIT NONE
      REAL(8)    ,PARAMETER     :: EPS=1.D-6
      REAL(8)    ,PARAMETER     :: DSMALL=1.D-12
      LOGICAL(4) ,PARAMETER     :: TTEST=.FALSE.
      INTEGER(4) ,INTENT(IN)    :: NB,NX           ! NUMBER OF STATES ,MAX
      INTEGER(4) ,INTENT(IN)    :: NAT             ! NUMBER OF ATOMS
      INTEGER(4) ,INTENT(IN)    :: NSP             ! NUMBER OF ATOM-TYPESS
      INTEGER(4) ,INTENT(IN)    :: NGW,NGWX        ! NUMBER OF PLANE WAVES, MAX
      INTEGER(4) ,INTENT(IN)    :: LMNX(NSP),LMNXX !
      INTEGER(4) ,INTENT(IN)    :: LNX(NSP),LNXX
      INTEGER(4) ,INTENT(IN)    :: LOX(LNXX,NSP)
      INTEGER(4) ,INTENT(IN)    :: ISPECIES(NAT)
      LOGICAL(4) ,INTENT(IN)    :: TGAMMA 
      REAL(8)    ,INTENT(IN)    :: CELLVOL           ! CELL VOLUME
      REAL(8)    ,INTENT(IN)    :: EMASS           ! FICTITIOUS ELECTRON MASS
      REAL(8)    ,INTENT(IN)    :: DELT            ! TIME STEP
      REAL(8)    ,INTENT(IN)    :: ANNEE           ! ELECTON FRICTION
      REAL(8)    ,INTENT(IN)    :: F(NB)           ! OCCUPATIONS
      REAL(8)    ,INTENT(IN)    :: DOVER(LNXX,LNXX,NSP)
      REAL(8)    ,INTENT(INOUT) :: LAMBDA(NX,NX)
      REAL(8)    ,INTENT(IN)    :: ONL(NAT,NB,LMNXX)
      REAL(8)    ,INTENT(INOUT) :: FNL(NAT,NB,LMNXX)
      COMPLEX(8) ,INTENT(IN)    :: CHI(NGWX,NB)     ! [PS-O]|PSI(0)>/M_PSI
      COMPLEX(8) ,INTENT(INOUT) :: PSI(NGWX,NB)     ! |PSI-BAR> / |PSI(+)>
      LOGICAL    ,INTENT(IN)    :: TSAFEORTHO
      REAL(8)                   :: PSIPSI(NX,NX)    
      REAL(8)                   :: CHIPSI(NX,NX)    
      REAL(8)                   :: CHICHI(NX,NX)    
      INTEGER(4)                :: IB1,IB2,I,J,IAT,ICOUNT,LMN
      INTEGER(4)                :: NTASKS,THISTASK
      REAL(8)                   :: SVAR,FAC
      REAL(8)                   :: FI,FJ
      INTEGER(4)                :: N2X,J0,I0,ISP,IMAX
      REAL(8)                   :: DIGAM
      INTERFACE 
        FUNCTION IDAMAX(LEN,ARRAY,STRIDE) RESULT(IPOS)
!       ==  RETURNS POSITION OF THE ELEMENT WITH THE LARGEST          ==
!       ==  ABSOLUTE VALUE (ESSL ROUTINE)                             ==
        INTEGER(4) :: LEN
        INTEGER(4) :: STRIDE
        REAL(8)    :: ARRAY(LEN)
        INTEGER(4) :: IPOS
        END FUNCTION IDAMAX
      END INTERFACE
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO')
!
!     ==================================================================
!     ==  CALCULATE  PSIPSI(I,J)=     <PSIBAR(I)|PSIBAR(J)>-1            ==
!     ==        AND  CHIPSI(I,J)=     <PSI0(I)|PSIBAR(J)>                ==
!     ==================================================================
      CALL WAVES_SIGSET(TGAMMA,1,CHI,CHI,NX,NB,NGWX,NGW,CELLVOL,CHICHI)
      CALL WAVES_SIGSET(TGAMMA,1,PSI,PSI,NX,NB,NGWX,NGW,CELLVOL,PSIPSI)
      CALL WAVES_SIGSET(TGAMMA,0,CHI,PSI,NX,NB,NGWX,NGW,CELLVOL,CHIPSI)
!
!     == ADD AUGMENTATION CONTRIBUTION TO THE OVERLAP ==================
      CALL WAVES_SIGSI(LNXX,LMNXX,NX,NAT,NSP,ISPECIES &
     &                ,LNX,LOX,NB,ONL,ONL,DOVER,CHICHI)
      CALL WAVES_SIGSI(LNXX,LMNXX,NX,NAT,NSP,ISPECIES &
     &                ,LNX,LOX,NB,FNL,FNL,DOVER,PSIPSI)
      CALL WAVES_SIGSI(LNXX,LMNXX,NX,NAT,NSP,ISPECIES &
     &                ,LNX,LOX,NB,ONL,FNL,DOVER,CHIPSI)
      CALL MPE$COMBINE('+',CHIPSI)
      CALL MPE$COMBINE('+',CHICHI)
      CALL MPE$COMBINE('+',PSIPSI)
!
!     ==================================================================
!     ==  CALCULATE LAGRANGE PARAMETERS                               ==
!     ==================================================================
      SVAR=2.D0*DELT**2/(1.D0+ANNEE)
      DO I=1,NB
        DO J=1,NB
          FI=F(I)
          FJ=F(J)
          IF(FI+FJ.LT.DSMALL)FJ=DSMALL
          LAMBDA(I,J)=LAMBDA(I,J)*SVAR*FI/(FI+FJ)
        ENDDO
      ENDDO
!
      IF (TSAFEORTHO) THEN
        CALL WAVES_ORTHO_X(NX,NB,F,CHICHI,PSIPSI,CHIPSI,LAMBDA)
      ELSE
        CALL WAVES_ORTHO_Y(NX,NB,PSIPSI,CHIPSI,CHICHI,LAMBDA)
      END IF
!
!     ==================================================================
!     ==  CALCULATE |PSI(+)>=|PSI>+|CHI>LAMBDA                        ==
!     ==================================================================
      DO J=1,NB
        DO I=1,NB
          CALL DZAXPY(2*NGW,LAMBDA(J,I),CHI(1,J),1,PSI(1,I),1,PSI(1,I),1)
        ENDDO
      ENDDO
      CALL MPE$QUERY(NTASKS,THISTASK)
      ICOUNT=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO IB1=1,NB
          ICOUNT=ICOUNT+1
!         __ SELECTION FOR PARALLEL PROCESSING__________________________
          IF(MOD(ICOUNT-1,NTASKS).EQ.THISTASK-1) THEN
            DO LMN=1,LMNX(ISP)
              SVAR=0.D0
              DO IB2=1,NB
                SVAR=SVAR+ONL(IAT,IB2,LMN)*LAMBDA(IB2,IB1)
              ENDDO
              FNL(IAT,IB1,LMN)=FNL(IAT,IB1,LMN)+SVAR
            ENDDO
!         __ELSE SELECTION FOR PARALLEL PROCESSING______________________
          ELSE
            DO LMN=1,LMNX(ISP)
              FNL(IAT,IB1,LMN)=0.D0
            ENDDO             
!         __SELECTION OF PARALLELPROCESSING FINISHED____________________
          END IF
        ENDDO
      ENDDO
      CALL MPE$COMBINE('+',FNL)
!
!     ==================================================================
!     ==  RESCALE GAMMA                                               ==
!     ==================================================================
!     == MASS IS NOT INCLUDED HERE (MASS TENSOR IS TAKEN CARE OF IN 
!     == PSIBAR AND (1/M)O|PSI> 
      FAC=(1.D0+ANNEE)/(2.D0*DELT**2)
      DO I=1,NB
        DO J=I,NB
          LAMBDA(I,J)=(LAMBDA(I,J)+LAMBDA(J,I))*FAC
          LAMBDA(J,I)=LAMBDA(I,J)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  TEST ORTHONORMALITY                                         ==
!     ==================================================================
      IF(TTEST) THEN
        CALL WAVES_SIGSET(TGAMMA,1,PSI,PSI,NX,NB,NGWX,NGW,CELLVOL,PSIPSI)
        CALL WAVES_SIGSI(LNXX,LMNXX,NX,NAT,NSP,ISPECIES &
     &                  ,LNX,LOX,NB,FNL,FNL,DOVER,PSIPSI)
        CALL MPE$COMBINE('+',PSIPSI)
        DO I=1,NB
          PSIPSI(I,I)=PSIPSI(I,I)-1.D0
        ENDDO
        N2X=NX*NX
        IMAX=IDAMAX(N2X,PSIPSI,1)
        J0=(IMAX-1)/NX+1
        I0=IMAX-(J0-1)*NX
        DIGAM=DABS(PSIPSI(I0,J0))
        IF(DIGAM.GT.EPS) THEN
          CALL ERROR$MSG('ORTHOGONALIZATION TEST FAILED')
          CALL ERROR$I4VAL('I',I0)
          CALL ERROR$I4VAL('J',J0)
          CALL ERROR$R8VAL('PSIPSI(I,J)',PSIPSI(I0,J0))
          CALL ERROR$STOP('WAVES_ORTHO')
        END IF
      END IF
                             CALL TRACE$POP
      RETURN
      END
!
!      .................................................................
       SUBROUTINE WAVES_ORTHO_Y(NX,NB,PHIPHI0,CHIPHI0,CHICHI0,RLAMBDA)
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
       INTEGER(4),INTENT(IN) :: NB,NX
       REAL(8)   ,INTENT(IN) :: PHIPHI0(NX,NB) !<PHIBAR_I|O|PHIBAR_J>
       REAL(8)   ,INTENT(IN) :: CHIPHI0(NX,NB) !<PHIBAR_I|O|CHI>
       REAL(8)   ,INTENT(IN) :: CHICHI0(NX,NB) !<CHI_I|O|CHI_J>
       REAL(8)   ,INTENT(OUT):: RLAMBDA(NX,NB) ! (I>J)=0
       REAL(8)               :: PHIPHI(NB,NB)
       REAL(8)               :: CHIPHI(NB,NB)
       REAL(8)               :: CHICHI(NB,NB)
       REAL(8)               :: ALPHA(NB,NB)   ! (I>J)=0
       REAL(8)               :: WORK(NB,NB)   ! (I>J)=0
       REAL(8)               :: DELTA(NB)
       INTEGER(4)            :: I,J,K,L,N,M,M1,M2
       REAL(8)               :: SVAR
       LOGICAL   ,PARAMETER  :: TPR=.FALSE.
       LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
       REAL(8)   ,PARAMETER  :: TOL=1.D-6
!      *****************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_Y')
!
!      =================================================================
!      ==  INITIALIZE                                                 ==
!      =================================================================
       DO I=1,NB
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
       DO N=1,NB
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
         DO M=N+1,NB
!          == PHIPHI(N,M)+CHIPHI(N,N)*DELTA(M)=0  ======================
           DELTA(M)=-PHIPHI(N,M)/CHIPHI(N,N)
           DO I=1,N
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
         DO M1=N+1,NB
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
         DO M=N+1,NB
!          == |CHI(M)>=|CHI(M)>+|CHI(N)>*DELTA(M) ============================
!          == CHIPHI(M,N)+CHIPHI(N,N)*DELTA(M)=0
           DELTA(M)=-CHIPHI(M,N)/CHIPHI(N,N)
         ENDDO
         DO M=N+1,NB
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
       REAL(8)              :: TEST(NB,NB)
       REAL(8)              :: TEST1(NB,NB)
       REAL(8)              :: TEST2(NB,NB)
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
      SUBROUTINE WAVES_ORTHO_X(NX,NB,F,CHICHI,PSIPSI,CHIPSI,LAMBDA)
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
      INTEGER(4),PARAMETER     :: MAX    = 100
      INTEGER(4),INTENT(IN)    :: NB,NX
      REAL(8)   ,INTENT(IN)    :: F(NX)
      REAL(8)   ,INTENT(INOUT) :: LAMBDA(NX,NX)
      REAL(8)   ,INTENT(IN)    :: PSIPSI(NX,NX)
      REAL(8)   ,INTENT(IN)    :: CHIPSI(NX,NX)   !
      REAL(8)   ,INTENT(IN)    :: CHICHI(NX,NX)
      REAL(8)   ,ALLOCATABLE   :: GAMN(:,:)  
      REAL(8)                  :: EIG(NX)
      INTEGER(4)               :: IND,ITER,I,J ! RUNNING VARIABLES
      INTEGER(4)               :: IMAX,I0,J0   ! AUXILARY VARIABLES
      REAL(8)                  :: DIGAM,SVAR,FI,FJ,EIGI ! AUXILARY VARIABLES
      REAL(8)                  :: HAUX(NX,NX)    
      REAL(8)                  :: U(NX,NX)       
!     ******************************************************************
                             CALL TRACE$PUSH('WAVES_ORTHO_X')
      ALLOCATE(GAMN(NX,NX))
!
!     ==================================================================
!     ==  CALCULATE  PSIPSI(I,J)= <PSIBAR(I)|PSIBAR(J)>-1             ==
!     ==        AND  CHIPSI(I,J)   = <PSI0(I)|PSIBAR(J)>                 ==
!     ==================================================================
!
!     ==================================================================
!     ==  DIAGONALIZE 0.5*(CHIPSI(I,J)+CHIPSI(J,I))                         ==
!     ==================================================================
      CALL DIAG(NX,NB,CHIPSI,EIG,U)
!
!     ==================================================================
!     ==================================================================
!     ==  ITERATIVE CALCULATION OF GAMMA                              ==
!     ==================================================================
!     ==================================================================
      DO ITER=1,MAX
!       ================================================================
!       ==  CALCULATE <PHI(+)|PHI(+)>-1 WITH PRESENT LAMBDA           ==
!       ==  GAMN(I,J)=PSIPSI(I,J)+LAMBDA(K,I)*CHIPSI(K,J)             ==
!       ==                       +CHIPSI(K,I)*LAMBDA(K,J)             == 
!       ==           +LAMBDA(K,I)*CHICHI(K,L)*LAMBDA(L,J)-1(I,J)      ==
!       ================================================================
!       __GAMN(I,J) = CHICHI(I,K)*LAMBDA(K,J)___________________________
        CALL DGEMUL(CHICHI,NX,'N',LAMBDA,NX,'N',HAUX,NX,NB,NB,NB)
!       __HAUX(I,J) = HAUX(I,J)+2*CHIPSI(I,J)___________________________
        CALL DAXPY(NX*NB,2.D0,CHIPSI,1,HAUX,1)
!       __GAMN(I,J) = LAMBDA(K,I)*HAUX(K,J)_____________________________
        CALL DGEMUL(LAMBDA,NX,'T',HAUX,NX,'N',GAMN,NX,NB,NB,NB)
!       __GAMN(I,J) = GAMN(I,J)-1_______________________________________
        DO I=1,NB
          DO J=I,NB
            SVAR=0.5D0*(GAMN(I,J)+GAMN(J,I))+PSIPSI(I,J)
            GAMN(J,I)=SVAR
            GAMN(I,J)=SVAR
          ENDDO
          GAMN(I,I)=GAMN(I,I)-1.D0
        ENDDO
!
!       ================================================================
!       == FIND LARGEST ELEMENT OF THE OVERLAP MATRIX                 ==
!       ================================================================
        IMAX=IDAMAX(NX*NX,GAMN(1,1),1)
        J0=(IMAX-1)/NX+1
        I0=IMAX-(J0-1)*NX
        DIGAM=DABS(GAMN(I0,J0))
!       PRINT*,'ITER ',ITER,I0,J0,DIGAM,NCON
        IF(DIGAM.LT.EPS) GOTO 9000
!
!       ==================================================================
!       ==  OBTAIN CHANGE OF THE LAMBDA MATRIX                          ==
!       ==================================================================
!       == TRANSFORM OVERLAP MATRIX GAMN
!       ----  HAUX(I,L)=U(K,I)*H0(K,L)
        CALL DGEMUL(U,NX,'T',GAMN,NX,'N',HAUX,NX,NB,NB,NB)
!       ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
        CALL DGEMUL(HAUX,NX,'N',U,NX,'N',GAMN,NX,NB,NB,NB)
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
        CALL DGEMUL(U,NX,'N',GAMN,NX,'N',HAUX,NX,NB,NB,NB)
!       ----  GAMN(I,J)=HAUX(I,L)*U(L,I)
        CALL DGEMUL(HAUX,NX,'N',U,NX,'T',GAMN,NX,NB,NB,NB)
!
!       ================================================================
!       ==  PROPAGATE GAMMA                                           ==
!       ================================================================
        DO I=1,NB
          FI=F(I)+DSMALL
          DO J=1,NB
            FJ=F(J)+DSMALL
            SVAR=2.D0*FI/(FI+FJ)
            LAMBDA(I,J)=LAMBDA(I,J)-SVAR*GAMN(I,J)
          ENDDO
        ENDDO
!
!       ================================================================
!       == SYMMETRIZE LAMBDA                                          ==
!       ================================================================
        DO I=1,NB
          DO J=1,NB
            IF(F(I).LT.F(J)) THEN
              LAMBDA(I,J)=LAMBDA(J,I)*(F(I)+DSMALL)/(F(J)+DSMALL)
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
      CALL ERROR$MSG('LOOP FOR ORTHOGONALIZATION IS NOT CONVERGED')
      CALL ERROR$STOP('WAVES_ORTHO_X')
!
9000  CONTINUE
      DEALLOCATE(GAMN)
                             CALL TRACE$POP
      RETURN
      END
!
!     .....................................................SIGSET.......
@PROCESS NOEXTCHK
      SUBROUTINE WAVES_SIGSET(TGAMMA,IDEN,VECT1,VECT2,NX,NB,NGWX,NGW &
     &                       ,CELLVOL,SIG)
!     ******************************************************************
!     **                                                              **
!     **  COMPUTES THE MATRIX                                         **
!     **     SIG(I,J) =  <PHI1(I)|PHI2(J)>                            **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
      IMPLICIT NONE
      LOGICAL   ,INTENT(IN) :: TGAMMA
      INTEGER(4),INTENT(IN) :: IDEN     ! SWITCH (1/0) VECT1=VECT2
      INTEGER(4),INTENT(IN) :: NB,NX    ! #(BANDS),MAX
      INTEGER(4),INTENT(IN) :: NGW,NGWX ! #(G-VECTORS),MAX
      COMPLEX(8),INTENT(IN) :: VECT1(NGWX,NB)
      COMPLEX(8),INTENT(IN) :: VECT2(NGWX,NB)
      REAL(8)   ,INTENT(OUT):: SIG(NX,NX)
      REAL(8)   ,INTENT(IN) :: CELLVOL
      REAL(8)               :: TWOOM,TR1,TR2
      INTEGER(4)            :: I,J,ID,IB1,IB2
!     ******************************************************************
!
!     =================================================================
!     ==  SIG(I,J)=VECT1(I)*VECT2(J)                                 ==
!     =================================================================
      IF(IDEN.EQ.1) THEN
        DO J=1,NB
          ID=NB-J+1
          CALL DGEMUL(VECT1(1,J),2*NGWX,'T',VECT2(1,J),2*NGWX,'N' &
     &                 ,SIG(J,J),NX,ID,2*NGW,1)
!         CALL DGEMUL(TRANSFER(VECT1(:,J:NB),SIG,2*NGWX*ID),2*NGWX,'T' &
!    &               ,TRANSFER(VECT2(:,J),SIG,2*NGWX),2*NGWX,'N' &
!    &               ,SIG(J,J),NX,ID,2*NGW,1)
!         CALL DGEMUL(VECT1(:,J:NB),2*NGWX,'T',VECT2(:,J),2*NGWX,'N' &
!    &                 ,SIG(J:NB,J),ID,ID,2*NGW,1)
          DO I=J,NB
            SIG(J,I)=SIG(I,J)
          ENDDO
        ENDDO
      ELSE
        CALL DGEMUL(VECT1,2*NGWX,'T' &
     &             ,VECT2,2*NGWX,'N' &
     &             ,SIG,NX,NB,2*NGW,NB)
!       CALL DGEMUL(TRANSFER(VECT1,SIG,2*NGWX*NB),2*NGWX,'T' &
!    &             ,TRANSFER(VECT2,SIG,2*NGWX*NB),2*NGWX,'N' &
!    &             ,SIG,NX,NB,2*NGW,NB)
      END IF
!
!     =================================================================
!     == RESCALE                                                     ==
!     =================================================================
      TWOOM=2.D0*CELLVOL
      DO IB2=1,NB
        DO IB1=1,NB
          SIG(IB1,IB2)=SIG(IB1,IB2)*TWOOM
        ENDDO
      ENDDO
!
!     =================================================================
!     == CORRECT THE GAMMA POINT                                     ==
!     =================================================================
      IF(TGAMMA) THEN
!       ==   REMARK: IMAGINARY PART OF VECT(G=0)=0
        DO J=1,NB
          DO I=1,NB
            TR1=DBLE(VECT1(1,I))
            TR2=DBLE(VECT2(1,J))
            SIG(I,J)=SIG(I,J)-TR1*TR2*CELLVOL
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_SIGSI(LNXX,LMNXX,NX,NAT,NSP,ISPECIES &
     &                      ,LNX,LOX,NB,FNL1,FNL2,DOVER,DOV)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE 1C CONTRIBUTION TO THE OVERLAP MATRIX        **
!     **  ( RESULT IS ADDED TO "DOV" )                                **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: LNXX
      INTEGER(4),INTENT(IN)   :: LMNXX
      INTEGER(4),INTENT(IN)   :: NAT            ! #(ATOMS)
      INTEGER(4),INTENT(IN)   :: NSP            ! #(ATOM TYPES)
      INTEGER(4),INTENT(IN)   :: NB,NX          ! #(BANDS),MAX
      INTEGER(4),INTENT(IN)   :: ISPECIES(NAT)  ! ATOM TYPE INDEX ARRAY
      INTEGER(4),INTENT(IN)   :: LNX(NSP)
      INTEGER(4),INTENT(IN)   :: LOX(LNXX,NSP)
      REAL(8)   ,INTENT(IN)   :: FNL1(NAT,NB,LMNXX)
      REAL(8)   ,INTENT(IN)   :: FNL2(NAT,NB,LMNXX)
      REAL(8)   ,INTENT(IN)   :: DOVER(LNXX,LNXX,NSP)
      REAL(8)   ,INTENT(INOUT):: DOV(NX,NX)
      INTEGER(4)              :: NUMTASK,ITASK
      INTEGER(4)              :: ICOUNT,IAT,ISP
      INTEGER(4)              :: LMN1,LN1,L1
      INTEGER(4)              :: LMN2,LN2,L2
      REAL(8)                 :: DOVER1,SVAR
      INTEGER(4)              :: IB1,IB2,IM
!     ******************************************************************
!
!     PRINT*,'====  DOV  BEFORE ====='
!     DO I=1,NB
!       WRITE(*,6000)(DOV(I,J),J=1,NB)
!     ENDDO
!
      CALL MPE$QUERY(NUMTASK,ITASK)
!     PRINT*,'NUMTASK,ITASK ',NUMTASK,ITASK
      ICOUNT=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LMN1=0
        DO LN1=1,LNX(ISP)
          L1=LOX(LN1,ISP)
          ICOUNT=ICOUNT+1
!         __ SELECTION FOR PARALLEL PROCESSING__________________________
          IF(MOD(ICOUNT-1,NUMTASK).EQ.ITASK-1) THEN
            LMN2=0
            DO LN2=1,LNX(ISP)
              L2=LOX(LN2,ISP)
              IF(L1.EQ.L2) THEN
                DOVER1=DOVER(LN1,LN2,ISP)
                DO IB1=1,NB
                  DO IB2=1,NB
                    SVAR=0.D0
                    DO IM=1,2*L1+1
                      SVAR=SVAR &
     &                    +FNL1(IAT,IB1,LMN1+IM)*FNL2(IAT,IB2,LMN2+IM)
                    ENDDO
                    DOV(IB1,IB2)=DOV(IB1,IB2)+DOVER1*SVAR
                  ENDDO
                ENDDO
              END IF
              LMN2=LMN2+2*L2+1
            ENDDO
!         __SELECTION OF PARALLELPROCESSING FINISHED____________________
          END IF
          LMN1=LMN1+2*L1+1
        ENDDO
      ENDDO
!
!     PRINT*,'====  DOV AFTER ====='
!     DO I=1,NB
!       WRITE(*,6000)(DOV(I,J),J=1,NB)
!6000   FORMAT(10F7.3)
!     ENDDO
!
      RETURN
      END
!
!     .....................................................GRAHAM.......
      SUBROUTINE WAVES_GRAMMSCHMIDT(TGAMMA,NB,NGW,PSI,FNL,DOVER &
     &                ,NGWX,LMNXX,LNXX &
     &                ,NSP,NAT,ISPECIES,LNX,LMNX,LOX,CELLVOL)
!     ******************************************************************
!     **                                                              **
!     ** GRAM-SCHMIDT ORTHOGONALIZATION OF THE WAVEFUNCTIONS PSI      **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
      USE MPE_MODULE
      IMPLICIT NONE
      INTERFACE
        FUNCTION ZDOTC(LENGTH,F1,STRIDE1,F2,STRIDE2)
!         == ESSL ROUTINE: SCALAR PRODUCT OF TWO COMPLEX(8) VECTORS   ==
          INTEGER(4),INTENT(IN) :: LENGTH
          COMPLEX(8),INTENT(IN) :: F1
          INTEGER(4),INTENT(IN) :: STRIDE1
          COMPLEX(8),INTENT(IN) :: F2
          INTEGER(4),INTENT(IN) :: STRIDE2
          COMPLEX(8)            :: ZDOTC
        END FUNCTION ZDOTC
        FUNCTION DNRM2(LENGTH,F,STRIDE)
!         == ESSL ROUTINE: NORM OF A REAL(8)VECTOR =====================
          INTEGER(4),INTENT(IN) :: LENGTH
          COMPLEX(8),INTENT(IN) :: F
          INTEGER(4),INTENT(IN) :: STRIDE
          REAL(8)               :: DNRM2         
        END FUNCTION DNRM2
      END INTERFACE
      LOGICAL   ,INTENT(IN)   :: TGAMMA
      INTEGER(4),INTENT(IN)   :: NB
      INTEGER(4),INTENT(IN)   :: NGW
      COMPLEX(8),INTENT(INOUT):: PSI(NGWX,NB)
      REAL(8)   ,INTENT(INOUT):: FNL(NAT,NB,LMNXX)
      REAL(8)   ,INTENT(IN)   :: DOVER(LNXX,LNXX,NSP)
      INTEGER(4),INTENT(IN)   :: NGWX
      INTEGER(4),INTENT(IN)   :: LMNXX
      INTEGER(4),INTENT(IN)   :: LNXX
      INTEGER(4),INTENT(IN)   :: NSP
      INTEGER(4),INTENT(IN)   :: NAT
      INTEGER(4),INTENT(IN)   :: ISPECIES(NAT)
      INTEGER(4),INTENT(IN)   :: LNX(NSP)
      INTEGER(4),INTENT(IN)   :: LMNX(NSP)
      INTEGER(4),INTENT(IN)   :: LOX(LNXX,NSP)
      REAL(8)   ,INTENT(IN)   :: CELLVOL
      COMPLEX(8)              :: SIGSIG,CSIG
      REAL(8)                 :: SIG
      REAL(8)                 :: DEVORTHO
      REAL(8)                 :: DEVNORM
      INTEGER(4)              :: IB1,IB2,LMN1,LMN2,LN1,LN2,L1,L2,I
      INTEGER(4)              :: IAT,ISP,LMN,IM
      REAL(8)                 :: DOVER1,SVAR,ANORM
      INTEGER(4)              :: NFILO
!     ******************************************************************
      DEVORTHO=0.D0
      DEVNORM=0.D0
      SIG=0.D0
!
!     ==================================================================
!     ==  ORTHOGONALIZATION                                           ==
!     ==================================================================
      DO IB1=1,NB
        IF(IB1.GT.1) THEN
          DO IB2=1,IB1-1
            CSIG=ZDOTC(NGW,PSI(1,IB2),1,PSI(1,IB1),1)
            SIG=2.D0*DBLE(CSIG)
            IF(TGAMMA) THEN
              SIG=SIG-DBLE(PSI(1,IB2)*PSI(1,IB1))
            END IF
!           ++  PARALLEL BEGIN +++++++++++++++++++++++++++++++++++++++++
            CALL MPE$COMBINE('+',SIG)
!           ++  PARALLEL END +++++++++++++++++++++++++++++++++++++++++++
            SIG=SIG*CELLVOL
!
            DO IAT=1,NAT
              ISP=ISPECIES(IAT)
              LMN1=0
              DO LN1=1,LNX(ISP)
                L1=LOX(LN1,ISP)
                LMN2=0
                DO LN2=1,LNX(ISP)
                  L2=LOX(LN2,ISP)
                  IF(L1.EQ.L2) THEN
                    DOVER1=DOVER(LN1,LN2,ISP)
                    SVAR=0.D0
                    DO IM=1,2*L1+1
                      SVAR=SVAR &
     &                    +FNL(IAT,IB1,LMN1+IM)*FNL(IAT,IB2,LMN2+IM)
                    ENDDO
                    SIG=SIG+DOVER1*SVAR
                  END IF
                  LMN2=LMN2+2*L2+1
                ENDDO
                LMN1=LMN1+2*L1+1
              ENDDO
            ENDDO
!
            DEVORTHO=MAX(DEVORTHO,SIG)
            SIGSIG=CMPLX(SIG,0.D0,KIND=8)
            CALL ZAXPY(NGW,-SIGSIG,PSI(1,IB2),1,PSI(1,IB1),1)
!
            DO IAT=1,NAT
              ISP=ISPECIES(IAT)
              DO LMN=1,LMNX(ISP)
                FNL(IAT,IB1,LMN)=FNL(IAT,IB1,LMN)-SIG*FNL(IAT,IB2,LMN)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
!       ================================================================
!       ==  NORMALIZATION                                             ==
!       ================================================================
        ANORM=DNRM2(2*NGW,PSI(1,IB1),1)
        ANORM=2.D0*ANORM*ANORM
        IF(TGAMMA) THEN
          ANORM=ANORM-DBLE(PSI(1,IB1)*CONJG(PSI(1,IB1)))
        END IF
        ANORM=ANORM*CELLVOL  
!       ++  PARALLEL BEGIN +++++++++++++++++++++++++++++++++++++++++++++
        CALL MPE$COMBINE('+',ANORM)
!       ++  PARALLEL END +++++++++++++++++++++++++++++++++++++++++++++++
!
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          LMN1=0
          DO LN1=1,LNX(ISP)
            L1=LOX(LN1,ISP)
            LMN2=0
            DO LN2=1,LNX(ISP)
              L2=LOX(LN2,ISP)
              IF(L1.EQ.L2) THEN
                DOVER1=DOVER(LN1,LN2,ISP)
                SVAR=0.D0
                DO IM=1,2*L1+1
                  SVAR=SVAR &
     &                +FNL(IAT,IB1,LMN1+IM)*FNL(IAT,IB1,LMN2+IM)
                ENDDO
                ANORM=ANORM+DOVER1*SVAR
              END IF
              LMN2=LMN2+2*L2+1
            ENDDO
            LMN1=LMN1+2*L1+1
          ENDDO
        ENDDO
!
        IF(ANORM.LE.0.D0) THEN
          CALL ERROR$MSG('NORM OF WAVE FUNCTION IS LESS OR EQUAL ZERO.')
          CALL ERROR$I4VAL('IB1',IB1)
          CALL ERROR$R8VAL('ANORM',ANORM)
          CALL ERROR$STOP('WAVES_GRAMMSCHMIDT')
        END IF
        ANORM=DSQRT(ANORM)
        DEVNORM=MAX(DEVNORM,SIG)
        DO I=1,NGW
          PSI(I,IB1)=PSI(I,IB1)/ANORM
        ENDDO
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          DO LMN=1,LMNX(ISP)
            FNL(IAT,IB1,LMN)=FNL(IAT,IB1,LMN)/ANORM
          ENDDO
        ENDDO
      ENDDO
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,FMT='("DEVIATIONS IN GRAMM-SCHMIDT:",2E10.2)') &
     &      DEVORTHO,DEVNORM
      RETURN
      END
!
!     ...................................................WAVES_DIAG.....
      SUBROUTINE WAVES_DIAG(NX,NGWX,LMNXX,NSP,NAT,ISPECIES,LMNX &
     &                ,NB,NGW,H0,HP,FNL0,CP,C0,RLAMP,RLAM0)
!     ******************************************************************
!     **                                                              **
!     **  APPLIES A UNITARY TRANSFORMATION TO THE ELCTRONIC           **
!     **  TRAJECTORIES SUCH THAT THE HAMILTONIAN FOR THE              **
!     **  NEXT TIME STEP IS DIAGONAL                                  **
!     **                                                              **
!     **  OUTPUT:                                                     **
!     **    EIG        EIGENVALUES OF THE HAMILTONIAN FOR T+DELT      **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL(4)  :: TFIRST
      INTEGER(4)  :: ISPECIES(NAT)
      INTEGER(4)  :: LMNX(NSP)
      REAL(8)     :: FNL0(NAT,NB,LMNXX)
      COMPLEX(8)  :: CP(NGWX,NB)
      COMPLEX(8)  :: C0(NGWX,NB)
      REAL(8)     :: H0(NX,NX)
      REAL(8)     :: HP(NX,NX)
      REAL(8)     :: RLAMP(NX,NX)
      REAL(8)     :: RLAM0(NX,NX)
      REAL(8)   ,ALLOCATABLE :: U(:,:)      !(NX,NX)
      REAL(8)   ,ALLOCATABLE :: FNLI(:)     !(NX)
      COMPLEX(8),ALLOCATABLE :: C3(:,:)     !(NGWX,NX)
      REAL(8)   ,ALLOCATABLE :: HAUX(:,:)   !(NX,NX)
      REAL(8)   ,ALLOCATABLE :: WORK2(:)
!     ******************************************************************
!
!     ==================================================================
!     ====  DIAGONALIZE M<PSIDOT|PSIDOT>+<PSI|H|PSI> FOR T=+DELT      ==
!     ==================================================================
      ALLOCATE(U(NX,NX))
!     IF(TFIRST) THEN
        ALLOCATE(WORK2(NX))
        CALL DIAG(NX,NB,HP,WORK2,U)
        DEALLOCATE(WORK2)
!
!     DO I=1,NB
!       DO J=1,NB
!         SVAR=0.D0
!         IF(I.EQ.J) SVAR=-1.D0
!         DO K=1,NB
!           SVAR=SVAR+U(I,K)*U(J,K)
!         ENDDO
!         IF(DABS(SVAR).GT.1.D-7) THEN
!           PRINT*,'ERROR ',I,J,SVAR
!         END IF
!       ENDDO
!     ENDDO
!
!     DO IB1=1,NB
!       DO IB2=1,NB
!         U(IB1,IB2)=0.D0
!       ENDDO
!       U(IB1,IB1)=1.D0
!     ENDDO
!     PHI=1.0D0 
!     I1=1
!     I2=5
!     U(I1,I1)=DCOS(PHI)
!     U(I2,I2)=DCOS(PHI)
!     U(I1,I2)=DSIN(PHI)
!     U(I2,I1)=-DSIN(PHI)
!
!     CALL MEMORY$ALLOCATE(8,$SAVE,NX*NX)
!     CALL BYTECOPY(8*NX*NX,U,SAVE)
!     DATA TFIRST/.TRUE./
!     TFIRST=.FALSE.
!     ELSE
!       CALL BYTECOPY(8*NX*NX,SAVE,U)
!     END IF
!     DO IB=1,NB
!       WRITE(*,FMT='(8F10.5)')(U(IB,IB2),IB2=1,NB)
!     ENDDO
!
!     ==================================================================
!     ====  APH(I,J)=U(K,I)*APH(K,L)*U(L,J)                           ==
!     ====  H0(I,J)=U(K,I)*H0(K,L)*U(L,J)                             ==
!     ==================================================================
!
      ALLOCATE(HAUX(NX,NX))
!
!     ==  TRANSFORM H0 ===============================================
!     ---- HAUX(I,L)=U(K,I)*H0(K,L)
      CALL DGEMUL(U,NX,'T',H0,NX,'N',HAUX,NX,NB,NB,NB)
!     ---- H0(I,J)=HAUX(I,L)*U(L,J)
      CALL DGEMUL(HAUX,NX,'N',U,NX,'N',H0,NX,NB,NB,NB)
!
!     ==  TRANSFORM HP   ===============================================
!     ---- HAUX(I,L)=U(K,I)*HP(K,L)
      CALL DGEMUL(U,NX,'T',HP,NX,'N',HAUX,NX,NB,NB,NB)
!     ---- HP(I,J)=HAUX(I,L)*U(L,J)
      CALL DGEMUL(HAUX,NX,'N',U,NX,'N',HP,NX,NB,NB,NB)
!
!     ==  TRANSFORM RLAMP  =============================================
!     ---- HAUX(I,L)=U(K,I)*RLAMP(K,L)
      CALL DGEMUL(U,NX,'T',RLAMP,NX,'N',HAUX,NX,NB,NB,NB)
!     ---- RLAMP(I,J)=HAUX(I,L)*U(L,J)
      CALL DGEMUL(HAUX,NX,'N',U,NX,'N',RLAMP,NX,NB,NB,NB)
!
!     ==  TRANSFORM RLAM0  =============================================
!     ---- HAUX(I,L)=U(K,I)*RLAM0(K,L)
      CALL DGEMUL(U,NX,'T',RLAM0,NX,'N',HAUX,NX,NB,NB,NB)
!     ---- RLAM0(I,J)=HAUX(I,L)*U(L,J)
      CALL DGEMUL(HAUX,NX,'N',U,NX,'N',RLAM0,NX,NB,NB,NB)
!
      DEALLOCATE(HAUX)
!
      DO I=1,NB
        DO J=I+1,NB
          H0(I,J)=0.5D0*(H0(I,J)+H0(J,I))
          H0(J,I)=H0(I,J)
          HP(I,J)=0.5D0*(HP(I,J)+HP(J,I))
          HP(J,I)=HP(I,J)
          RLAM0(I,J)=0.5D0*(RLAM0(I,J)+RLAM0(J,I))
          RLAM0(J,I)=RLAM0(I,J)
          RLAMP(I,J)=0.5D0*(RLAMP(I,J)+RLAMP(J,I))
          RLAMP(J,I)=RLAMP(I,J)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ====  ROTATE WAVE FUNCTIONS |PSI(I)>=|PSI(J)>*U(J,I)            ==
!     ==================================================================
      ALLOCATE(FNLI(NX))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LMN=1,LMNX(ISP)
          DO IB=1,NB
            FNLI(IB)=FNL0(IAT,IB,LMN)
            FNL0(IAT,IB,LMN)=0.D0
          ENDDO
          DO IB1=1,NB
            DO IB2=1,NB
              FNL0(IAT,IB1,LMN)=FNL0(IAT,IB1,LMN)+FNLI(IB2)*U(IB2,IB1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(FNLI)
!
      ALLOCATE(C3(NGWX,NX))
      CALL WAVES_SETEQUAL(NGWX,NGW,NB,C0,C3)
      CALL WAVES_SETZERO(NGWX,NGW,NB,C0)
!
      DO I=1,NB
        DO J=1,NB
!         ---- C0(I)=C0(I)+C3(J)*U(J,I) -------
          CALL DZAXPY(2*NGW,U(J,I),C3(1,J),1,C0(1,I),1,C0(1,I),1)
        ENDDO
      ENDDO
!
      CALL WAVES_SETEQUAL(NGWX,NGW,NB,CP,C3)
      CALL WAVES_SETZERO(NGWX,NGW,NB,CP)
!
      DO I=1,NB
        DO J=1,NB
!         ---- CP(I)=CP(I)+C3(J)*U(J,I) -------
          CALL DZAXPY(2*NGW,U(J,I),C3(1,J),1,CP(1,I),1,CP(1,I),1)
        ENDDO
      ENDDO
      DEALLOCATE(C3)
      DEALLOCATE(U)
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE WAVES_PRHAMILTONIAN(STRING,NFIL,NX,NB,HAMILTON)
!     ******************************************************************
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRING ! DESCRIPTION
      INTEGER(4)  ,INTENT(IN) :: NFIL   ! FORTRAN FILE NUMBER
      INTEGER(4)  ,INTENT(IN) :: NB,NX  ! #(BANDS),MAX
      REAL(8)     ,INTENT(IN) :: HAMILTON(NX,NX) ! HAMILTON MATRIX
      INTEGER(4)              :: IB1,IB2,I1,I2
!     ******************************************************************
      WRITE(NFIL,FMT='(A)')STRING
!     WRITE(NFIL,*)'NX,NB,HAM',NX,NB,HAMILTON(3,3)
      DO IB1=1,NB
        I1=1
        DO WHILE (I1.LE.NB)
          I2=MIN(I1+9,NB)
          WRITE(NFIL,FMT='(I5,10E19.2)') &
     &         IB1,(HAMILTON(IB1,IB2),IB2=I1,I2)
          I1=I2+1
        ENDDO
      ENDDO
      RETURN
      END
!
!     .....................................................INVOV........
      SUBROUTINE WAVES_1CFORCE(IAT,NAT,NX,LMNXX &
     &                     ,NB,LMNX,WKPT,F,H,DFNL,HNL,ONL,FION)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE FORCE DUE TO THE POSITION DEPENDENCE         **
!     **  OF THE OVERLAP OPERATOR                                     **
!     **                                                              **
!     ************P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1992)***
      IMPLICIT NONE
      INTEGER(4),PARAMETER  :: IPR=0
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: NX
      INTEGER(4),INTENT(IN) :: LMNXX
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(IN) :: WKPT
      REAL(8)   ,INTENT(IN) :: F(NB)
      REAL(8)   ,INTENT(IN) :: H(NX,NX)   ! ORTHO-LAGRANGE MULTIPLIER
      REAL(8)   ,INTENT(IN) :: DFNL(NAT,NB,LMNXX,3)
      REAL(8)   ,INTENT(IN) :: HNL(NAT,NB,LMNXX)  
      REAL(8)   ,INTENT(IN) :: ONL(NAT,NB,LMNXX)
      REAL(8)   ,INTENT(INOUT):: FION(3)
      INTEGER(4)            :: IDIR,IB,LMN,IB1,IB2,I
      REAL(8)               :: FAC,SVAR
      REAL(8)               :: F1,F2
!     ******************************************************************
!
!     ==================================================================
!     ==  PRINTOUT FOR TEST                                           ==
!     ==================================================================
      IF(IPR.EQ.1) THEN
        PRINT*,'OVFOR START'
        WRITE(*,FMT='("FORCE ON ATOM ",I3,":",3F15.5)') &
     &               IAT,(FION(I),I=1,3)
      ENDIF
!
!     ==================================================================
!     ==  2 * TR[ <PSI|P> DH <GRAD*P|PSI> ]                           ==
!     ==================================================================
      DO IDIR=1,3
        DO IB=1,NB
          DO LMN=1,LMNX
            FAC=2.D0*F(IB)*WKPT
            FION(IDIR)=FION(IDIR) &
     &                -FAC*DFNL(IAT,IB,LMN,IDIR)*HNL(IAT,IB,LMN)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  2*TR[ <PSI|P>DO<GRAD*P|PSI> * LAMBDA ]                      ==
!     ==================================================================
      DO LMN=1,LMNX
        DO IB1=1,NB
          F1=F(IB1)
          SVAR=0.D0
          DO IB2=1,NB
            F2=F(IB2)
            SVAR=SVAR + H(IB1,IB2)* ONL(IAT,IB2,LMN) &
     &                * 2.D0*F1*F2/(F1+F2+1.D-12)
          ENDDO
          SVAR=SVAR*2.D0*WKPT
          FION(1)=FION(1)+DFNL(IAT,IB1,LMN,1)*SVAR
          FION(2)=FION(2)+DFNL(IAT,IB1,LMN,2)*SVAR
          FION(3)=FION(3)+DFNL(IAT,IB1,LMN,3)*SVAR
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  PRINTOUT FOR TEST                                           ==
!     ==================================================================
      IF(IPR.EQ.1) THEN
        PRINT*,'AFTER OVFOR'
        WRITE(*,FMT='("FORCE ON ATOM ",I3,":",3F15.5)') &
     &               IAT,(FION(I),I=1,3)
      ENDIF 
      RETURN
      END
