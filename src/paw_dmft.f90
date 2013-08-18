!........1.........2.........3.........4.........5.........6.........7.........8
MODULE DMFT_MODULE
!*******************************************************************************
!** IPROOFCHI(NCHI) IS A POINTER THAT SELECTS A PARTICULAR LOCAL ORBITAL      **
!**                                                                           **
!**  DESCRIPTION OF VARIABLES                                                 **
!**    DIAGSLOC:   |CHI_NON-ORTHONORMAL>*DIAGSLOC = |CHI_ORTHONORMAL>         **
!**    GLOCLAUR1,1,2 EXPANSION OF SPIN-AVERAGED GLOC IN 1/(CI*HBAR*OMEGA)     **
!**                                                                           **
!*******************************************************************************
TYPE DENMAT_TYPE
INTEGER(4)          :: LNX
INTEGER(4),POINTER  :: LOX(:) !(LNX)
INTEGER(4),POINTER  :: LMN(:) !(NLOC)
INTEGER(4)          :: LMNX
COMPLEX(8),POINTER  :: MAT(:,:,:) => NULL() ! DENSITY MATRIX
COMPLEX(8),POINTER  :: H(:,:,:)   => NULL() ! HAMILTONIAN FROM DOUBLE COUNTING
END TYPE DENMAT_TYPE
TYPE ATOMSET_TYPE
  INTEGER(4)           :: NLOC
  INTEGER(4)           :: ICHI1
  INTEGER(4)           :: ICHI2
  REAL(8)   ,POINTER   :: U(:,:,:,:)        => NULL() !(NLOC,NLOC,NLOC,NLOC) 
  COMPLEX(8),POINTER   :: GLOC(:,:,:,:)     => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: GLOCLAUR(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,3)
  COMPLEX(8),POINTER   :: SLOC(:,:,:,:)     => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: SLOCLAUR(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,3)
  COMPLEX(8),POINTER   :: DIAGSLOC(:,:,:)   => NULL() !(NLOC,NLOC,NDIMD)
  TYPE(DENMAT_TYPE)    :: DENMAT
END TYPE ATOMSET_TYPE
TYPE KSET_TYPE
  REAL(8)            :: WKPTL
  LOGICAL(4)         :: TADDMINUSK    !ADD  THE TERM FOR -K
  REAL(8)   ,POINTER :: ERHO(:,:)     !(NB,NSPIN)
  COMPLEX(8),POINTER :: PIPSI(:,:,:,:)!(NDIM,NCHI,NB,NSPIN)
  COMPLEX(8),POINTER :: RHO(:,:,:)    !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: H0(:,:,:)     !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: HRHO(:,:,:)   !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: SINV(:,:,:)   !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: SMAT(:,:,:)   !(NCHI,NCHI,NDIMD)
END TYPE KSET_TYPE
LOGICAL(4),PARAMETER   :: TON=.TRUE.
LOGICAL(4),SAVE        :: TINI=.FALSE.
REAL(8)   ,PARAMETER   :: AMIX=1.D-1
INTEGER(4)             :: NOMEGA
INTEGER(4)             :: NCHI          ! #(CORRELATED ORBITALS)
INTEGER(4)             :: NB            ! #(BAND STATES PER K-POINT)
INTEGER(4)             :: NKPTL         ! #(KPOINTS ON THIS TASK)
INTEGER(4)             :: NSPIN         ! #(SPIN COMPONENTS)
INTEGER(4)             :: NDIM          ! #(SPINOR COMPONENTS)
INTEGER(4)             :: NDIMD         ! CAN BE 1,2,4
INTEGER(4)             :: NAT           ! #(ATOMS)
REAL(8)   ,ALLOCATABLE :: OMEGA(:)      ! MATSUBARA FREQUENCIES
INTEGER(4),ALLOCATABLE :: IPROOFCHI(:)  !(NCHI) MAP ICHI TO IPRO
REAL(8)                :: KBT           ! TEMPERATURE (K_B*T)
REAL(8)                :: MU            ! CHEMICAL POTENTIAL
REAL(8)                :: DELTAT        ! TIMESTEP
!== KSET =======================================================================
TYPE(KSET_TYPE)   ,ALLOCATABLE :: KSET(:)  !(NKPTL)
TYPE(ATOMSET_TYPE),ALLOCATABLE :: ATOMSET(:)  !(NAT)
END MODULE DMFT_MODULE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_INI()
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TINI,NDIM,NDIMD,NSPIN,NKPTL,NB,NAT,NCHI &
     &                       ,NOMEGA,KBT,DELTAT,MU,OMEGA,IPROOFCHI &
     &                       ,KSET,ATOMSET
      USE WAVES_MODULE, ONLY : KMAP,NDIM_W=>NDIM,NKPTL_W=>NKPTL,NSPIN_W=>NSPIN
      IMPLICIT NONE
      REAL(8)                :: PI
      INTEGER(4)             :: NTASKS_K,THISTASK_K
      INTEGER(4)             :: NTASKS_M,THISTASK_M
      INTEGER(4)             :: NKPT
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: NLOC
      INTEGER(4)             :: L
      INTEGER(4)             :: NU,ISP,IAT,LN,IM,IKPTL,IKPT,ICHI,IPRO,I
      INTEGER(4)             :: LMN
      INTEGER(4)             :: I1,I2 ! BOUNDS ON CHI-ARRAY FOR EACH ATOM
      REAL(8)   ,ALLOCATABLE :: WKPT(:) !(NKPT) K-POINT WEIGHTS
      INTEGER(4),ALLOCATABLE :: LOX(:)  !(LNX)
      LOGICAL(4),ALLOCATABLE :: TORB(:) !(LNX)
      REAL(8)                :: EV
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
      PI=4.D0*ATAN(1.D0)
      CALL CONSTANTS$GET('EV',EV)
      CALL MPE$QUERY('K',NTASKS_K,THISTASK_K)
      CALL MPE$QUERY('MONOMER',NTASKS_M,THISTASK_M)
!
!     ==========================================================================
!     == OTHER VARIABLES                                                      ==
!     ==========================================================================
!THIS IS HARD WIRED
      NOMEGA=40
      KBT=0.333D0*EV
      DELTAT=10.D0
      MU=0.D0
!
!     ==========================================================================
!     == COLLECT PERMANENT DATA                                               ==
!     ==========================================================================
      NDIM=NDIM_W   ! IMPORT NDIM FROM WAVES_MODULE INTO DMFT MODULE
      NSPIN=NSPIN_W ! IMPORT NSPIN FROM WAVES_MODULE INTO DMFT MODULE
      IF(NDIM.EQ.1.AND.NSPIN.EQ.1) THEN
        NDIMD=1
      ELSE IF(NDIM.EQ.1.AND.NSPIN.EQ.2) THEN
        NDIMD=2
      ELSE IF(NDIM.EQ.2.AND.NSPIN.EQ.1) THEN
        NDIMD=4
      ELSE
        CALL ERROR$MSG('INCONSISTENT SPIN VARABLED INHERITED FROM PAW_WAVES.')
        CALL ERROR$MSG('(NDIM,NSPIN) MAY BE (1,1), (1,2) OR (2,1)')
        CALL ERROR$I4VAL('NDIM',NDIM)
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('DMFT_INI')
      END IF
!
      NKPTL=NKPTL_W ! IMPORT NKPTL FROM WAVES_MODULE INTO DMFT MODULE
      CALL DYNOCC$GETI4('NB',NB)
!
!     ==========================================================================
!     ==  CALCULATE MATSUBARA FREQUENCIES                                     ==
!     ==========================================================================
      ALLOCATE(OMEGA(NOMEGA))
      DO NU=1,NOMEGA
        OMEGA(NU)=REAL(2*NU-1,KIND=8)*PI*KBT
      ENDDO
!
!     ==========================================================================
!     == SELECT CORRELATED ORBITALS (->NCHI, ATOMSET)                         ==
!     ==========================================================================
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(ATOMSET(NAT))
      ATOMSET(:)%ICHI1=1
      ATOMSET(:)%ICHI2=0
!
!     == GET NCHI AND BOUNDS ON ATOMSET ========================================
      NCHI=0
      DO IAT=1,NAT
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(TORB(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETL4A('TORB',LNX,TORB)
        CALL SETUP$UNSELECT()
        I1=NCHI+1   ! FIRST ICHI FOR THIS ATOM
        DO LN=1,LNX
          L=LOX(LN)
          IF(TORB(LN)) THEN
            NCHI=NCHI+2*L+1
          END IF
        ENDDO
        I2=NCHI
PRINT*,'TORB ',IAT,TORB
        DEALLOCATE(LOX)
        DEALLOCATE(TORB)
!
!       == SAVE BOUNDS TO ATOMSET ==============================================
        NLOC=I2-I1+1
        ATOMSET(IAT)%NLOC=NLOC
        ATOMSET(IAT)%ICHI1=I1
        ATOMSET(IAT)%ICHI2=I2
!
!       == ALLOCATE ATOMSET SUBARRAYS ==========================================
        ALLOCATE(ATOMSET(IAT)%U(NLOC,NLOC,NLOC,NLOC))
        ALLOCATE(ATOMSET(IAT)%GLOC(NLOC,NLOC,NDIMD,NOMEGA))
        ALLOCATE(ATOMSET(IAT)%GLOCLAUR(NLOC,NLOC,NDIMD,3))
        ALLOCATE(ATOMSET(IAT)%SLOC(NLOC,NLOC,NDIMD,NOMEGA))
        ALLOCATE(ATOMSET(IAT)%SLOCLAUR(NLOC,NLOC,NDIMD,3))
        ALLOCATE(ATOMSET(IAT)%DIAGSLOC(NLOC,NLOC,NDIMD))
        ATOMSET(IAT)%U=0.D0
        ATOMSET(IAT)%GLOC=(0.D0,0.D0)
        ATOMSET(IAT)%GLOCLAUR=(0.D0,0.D0)
        ATOMSET(IAT)%SLOC=(0.D0,0.D0)
        ATOMSET(IAT)%SLOCLAUR=(0.D0,0.D0)
        ATOMSET(IAT)%DIAGSLOC=(0.D0,0.D0)
      ENDDO
PRINT*,'NCHI ',NCHI
!
!     == ACCUMULATE IPROOFCHI ==================================================
      ALLOCATE(IPROOFCHI(NCHI))
      ICHI=0
      IPRO=0
      DO IAT=1,NAT
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(TORB(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETL4A('TORB',LNX,TORB)
        DO LN=1,LNX
          L=LOX(LN)
          IF(TORB(LN)) THEN
            DO IM=1,2*L+1
              ICHI=ICHI+1
              IPRO=IPRO+1
              IPROOFCHI(ICHI)=IPRO
            ENDDO
          ELSE
            IPRO=IPRO+2*L+1
          END IF
        ENDDO
        DEALLOCATE(LOX)
        DEALLOCATE(TORB)
        CALL SETUP$UNSELECT()
      ENDDO
PRINT*,'IPROOFCHI ',IPROOFCHI
!
!     ==========================================================================
!     == COLLECT K-POINT WEIGHTS                                              ==
!     ==========================================================================
      ALLOCATE(KSET(NKPTL))
      DO IKPTL=1,NKPTL
        ALLOCATE(KSET(IKPTL)%ERHO(NB,NSPIN))
        ALLOCATE(KSET(IKPTL)%PIPSI(NDIM,NCHI,NB,NSPIN))
        ALLOCATE(KSET(IKPTL)%RHO(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%H0(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%HRHO(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%SINV(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%SMAT(NCHI,NCHI,NDIMD))
        KSET(IKPTL)%WKPTL=0.D0
        KSET(IKPTL)%ERHO=0.D0
        KSET(IKPTL)%PIPSI=(0.D0,0.D0)
        KSET(IKPTL)%RHO=(0.D0,0.D0)
        KSET(IKPTL)%HRHO=(0.D0,0.D0)
        KSET(IKPTL)%H0=(0.D0,0.D0)
        KSET(IKPTL)%SINV=(0.D0,0.D0)
        KSET(IKPTL)%SMAT=(0.D0,0.D0)
      ENDDO
! 
!     == COLLECT K-POINT WEIGHTS ===============================================
      CALL DYNOCC$GETI4('NKPT',NKPT)
      ALLOCATE(WKPT(NKPT))
      CALL DYNOCC$GETR8A('WKPT',NKPT,WKPT)
      IKPTL=0
      DO IKPT=1,NKPT
        IF(KMAP(IKPT).EQ.THISTASK_M) THEN
          IKPTL=IKPTL+1
          IF(IKPTL.GT.NKPTL) THEN
            CALL ERROR$MSG('KMAP INCONSISTENT WITH NKPTL FROM PAW_WAVES')
            CALL ERROR$I4VAL('NKPTL',NKPTL)
            CALL ERROR$STOP('DMFT_INI')
          END IF
          KSET(IKPTL)%WKPTL=WKPT(IKPT)
        END IF
      ENDDO
      IF(IKPTL.NE.NKPTL) THEN
        CALL ERROR$MSG('KMAP INCONSISTENT WITH NKPTL FROM PAW_WAVES')
        CALL ERROR$I4VAL('NKPTL',NKPTL)
        CALL ERROR$I4VAL('IKPTL',IKPTL)
        CALL ERROR$STOP('DMFT_INI')
      END IF
!
!     ==========================================================================
!     ==  ALLOCATE PERMANENT ARRAYS
!     ==========================================================================
!
!     == ATOMSET ===============================================================
      DO IAT=1,NAT
      ENDDO
!
!     == ATOMSET%DENMAT%... FOR DOUBLE COUNTING CORRECTION =====================
      DO IAT=1,NAT
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        ALLOCATE(TORB(LNX))
        CALL SETUP$GETL4A('TORB',LNX,TORB)
        CALL SETUP$UNSELECT()
!
        ATOMSET(IAT)%DENMAT%LNX=LNX
        ALLOCATE(ATOMSET(IAT)%DENMAT%LOX(LNX))
        ATOMSET(IAT)%DENMAT%LOX=LOX
        LMNX=SUM(2*LOX+1)
        ATOMSET(IAT)%DENMAT%LMNX=LMNX
        NLOC=ATOMSET(IAT)%NLOC
        ALLOCATE(ATOMSET(IAT)%DENMAT%LMN(NLOC))
        LMN=0
        I=0
        DO LN=1,LNX
          L=LOX(LN)
          IF(TORB(LN)) THEN
            DO IM=1,2*L+1
              I=I+1
              LMN=LMN+1
              ATOMSET(IAT)%DENMAT%LMN(I)=LMN
            ENDDO
          ELSE
            LMN=LMN+2*L+1
          END IF
        ENDDO
        ALLOCATE(ATOMSET(IAT)%DENMAT%MAT(LMNX,LMNX,NDIMD))
        ALLOCATE(ATOMSET(IAT)%DENMAT%H(LMNX,LMNX,NDIMD))
        ATOMSET(IAT)%DENMAT%MAT=(0.D0,0.D0)
        ATOMSET(IAT)%DENMAT%H  =(0.D0,0.D0)
        DEALLOCATE(LOX)
        DEALLOCATE(TORB)
      ENDDO        
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$GREEN()
!     **************************************************************************
!     ** CALCULATES THE LOCAL INTERACTING GREENS FUNCTION                     **
!     **                                                                      **
!     ** PIPSI <PI_A|PSI_N> IS THE PRE-FACTOR OF LOCAL ORBITAL |CHI_A> IN     **
!     **       THE LOCAL ORBITAL EXPANSION OF |PSI_N>                         **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON
      USE MPE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: ITER
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT$GREEN')
      CALL DMFT_INI()
!
!     ==========================================================================
!     ==  COLLECT DFT HAMILTONIAN  (HAMILTONIAN INITIALIZED TO ERHO)          ==
!     ==========================================================================
      CALL DMFT$COLLECTHAMILTONIAN_WITHKSET()  
      CALL DMFT$COLLECTFULLDENMAT()  
!
!     ==========================================================================
!     ==  OBTAIN BARE U-TENSOR FROM LMTO OBJECT.                              ==
!     ==     SHAPE OF ORBITALS ARE DEFINED BY NPRO AND ATOMIC STRUCTRE        ==
!     ==     ORBITALS ARE SELECTED BY TORB.                                   ==
!     ==========================================================================
PRINT*,'SETTING UP U-TENSOR...'
      CALL DMFT_UCHI_WITHATOMSET() 
PRINT*,'... U-TENSOR DONE'
!
!     ==========================================================================
!     ==  TRANSFORMATION ONTO A ORTHONORMAL CORRELATED BASIS SET              ==
!     ==    |CHIORTHO>   =|CHINONORTHO>*DIAGSLOC                              ==
!     ==========================================================================
      CALL DMFT_DIAGSLOC_WITHKSET() 
!
!     ==========================================================================
!     ==  CONSTRUCT NON-INTERACTING HAMILTONIAN THAT PRODUCES THE CORRECT     ==
!     ==  ONE-PARTICLE DENSITY MATRIX                                         ==
!     ==========================================================================
PRINT*,'SETTING UP NON-INTERACTING HAMILTONIAN HRHO...'
      CALL DMFT_CONSTRAINTS_WITHKSET('HRHO') !NEW
PRINT*,'... HRHO DONE'
!
!     ==========================================================================
!     == ITERATION TO ENFORCE CONSTRAINTS                                     ==
!     ==========================================================================
      DO ITER=1,10
WRITE(*,FMT='(82("="),T20," ITERATION ",I5)')ITER
        CALL DMFT_GLOC_WITHATOMSET() 
!
!       ========================================================================
!       ==  CALL THE SOLVER                                                   ==
!       ========================================================================
        CALL DMFT$SOLVER_WITHKSET() 
!
!       ========================================================================
!       ==  CONSTRAINTS                                                       ==
!       ========================================================================
        CALL DMFT_CONSTRAINTS_WITHKSET('H0')
!PRINT*,'FORCED STOP AFTER CONSTRAINTS_WITHKSET("H0")'
!STOP
        PRINT*,'ITERATION COMPLETED ',ITER
      ENDDO ! END OF LOOP OVER ITERATIONS TO ENFORCE CONSTRAINT
!
      CALL DMFT$ADDTOHPSI_WITHKSET()

!STOP 'END OF LOOP. STOPPING.'
                                       CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_DIAGSLOC_WITHKSET()
!     **************************************************************************
!     **  DETERMINE OVERLAP MATRIX SMAT(NCHI,NCHI,ISPIN) OF LOCAL ORBITALS,   **
!     **  THE K-DEPENDENT INVERSE SINV(NCHI,NCHI,IKPT,ISPIN)                  **
!     **  AND THE MATRIX THAT CONVERTS THE LOCAL ORBITALS                     **
!     **  INTO AN ORTHONORMAL BASISSET OF CORRELATED ORBITALS                 **
!     **                                                                      **
!     **   |PI_NEU>  = |PI_ALT>  DIAGSLOC                                     **
!     **   |PI_ALT>  = |PI_NEU>  INVERT(DIAGSLOC)                             **
!     **   |CHI_ALT> = |CHI_NEU> TRANSPOSE(DIAGSLOC)                          **
!     **                                                                      **
!     **  WHERE THE PI ARE PROJECTOR FUNCTIONS AND CHI ARE ORBITALS           **
!     **                                                                      **
!     **   SINV(K) =SUM_N <PI|PSI(N,K)><PSI(N,K)|PI>                          **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NSPIN,NDIM,NDIMD,KSET,ATOMSET,NAT
      IMPLICIT NONE
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
      INTEGER(4)            :: IKPT,ISPIN,IDIM1,IDIM2,IDIMD,I
      COMPLEX(8),ALLOCATABLE:: MAT1(:,:) !(2*NCHI,2*NCHI)
      COMPLEX(8),ALLOCATABLE:: MAT2(:,:) !(2*NCHI,2*NCHI)
      REAL(8)   ,ALLOCATABLE:: FLOC(:)  !(NLOC)
      INTEGER(4)            :: I1,I2,NLOC,IAT
!     **************************************************************************
!
!     ==========================================================================
!     == CALCULATE INVERSE OVERLAP MATRIX                                     ==
!     ==   SINV(A,S1,B,S2,K)=SUM_N <PI(A,S1)|PSI(N,K)><PSI(N,K)|PI(B,S2)>     ==
!     ==                                                                      ==
!     ==   NON-SPIN POLARIZED:                                                ==
!     ==     SINV(S1,S2)=1/2*[ SINV(1)*(1,0/0,1) ]                            ==
!     ==   COLLINEAR SPIN-POLARIZED                                           ==
!     ==     SINV(S1,S2)=1/2*[ SINV(1)*(1,0/0,1)+SINV(2)*(1,0/0,-1) ]         ==
!     ==   NON-COLLINEAR SPIN-POLARIZED                                       ==
!     ==     SINV(S1,S2)=1/2*[ SINV(1)*(1,0/0,1)+SINV(2)*(0,1/1,0)            ==
!     ==                      +SINV(3)*(0,-I/I,0)+SINV(2)*(1,0/0,-1) ]        ==
!     ==   WHERE                                                              ==
!     ==     SINV(1)=SINV(1,1)+  SINV(2,2)                                    ==
!     ==     SINV(2)=SINV(1,2)+I*SINV(2,1)                                    ==
!     ==     SINV(3)=SINV(1,2)-I*SINV(2,1)                                    ==
!     ==     SINV(4)=SINV(1,1)-  SINV(2,2)  !FOR NSPIN=2 THIS IS SIMV(2)      ==
!     ==                                                                      ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        KSET(IKPT)%SINV(:,:,:)=(0.D0,0.D0)
        DO ISPIN=1,NSPIN
          DO IDIM2=1,NDIM
            DO IDIM1=1,NDIM
!             == NON-COLLINEAR               IDIMD=(1)
!             == COLLINEAR SPIN-POLARIED     IDIMD=(1,2)
!             == NON-COLLINEAR SPIN-POLARIED IDIMD=(1,2,3,4)
              IDIMD=IDIM1+NDIM*(IDIM2-1)+ISPIN-1
              KSET(IKPT)%SINV(:,:,IDIMD) &
     &                     =MATMUL(KSET(IKPT)%PIPSI(IDIM1,:,:,ISPIN) &
     &                     ,CONJG(TRANSPOSE(KSET(IKPT)%PIPSI(IDIM2,:,:,ISPIN))))
              
            ENDDO
          ENDDO
        ENDDO
        CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%SINV)
      ENDDO
!
!     ==========================================================================
!     == OVERLAP MATRIX                                                       ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        CALL SPINOR$INVERT(NDIMD,NCHI,KSET(IKPT)%SINV,KSET(IKPT)%SMAT)
      ENDDO
!
!     ==========================================================================
!     == OBTAIN MATRIX THAT MAKES THE OVERLAP EQUAL TO THE UNIT MATRIX        ==
!     ==========================================================================
      DO IAT=1,NAT
        NLOC=ATOMSET(IAT)%NLOC
        IF(NLOC.LE.0) CYCLE
        I1=ATOMSET(IAT)%ICHI1
        I2=ATOMSET(IAT)%ICHI2
        ATOMSET(IAT)%DIAGSLOC(:,:,:)=(0.D0,0.D0)
        DO IKPT=1,NKPTL
          ATOMSET(IAT)%DIAGSLOC=ATOMSET(IAT)%DIAGSLOC &
    &                           +KSET(IKPT)%WKPTL*KSET(IKPT)%SINV(I1:I2,I1:I2,:)
        ENDDO
!
!       == MAKE REAL IF TIME-INVERSION SYMMETRY IS VALID =======================
        IF(NDIMD.LT.4) ATOMSET(IAT)%DIAGSLOC=REAL(ATOMSET(IAT)%DIAGSLOC)
!
!       == CONVERT FROM (TXYZ) TO (UP-DOWN REPRESENTATION) =====================
        CALL SPINOR$CONVERT('BACK',NLOC,NDIMD,ATOMSET(IAT)%DIAGSLOC)
!
!       == DIAGONALIZE AND FORM TRANSFORMATION MATRIX ==========================
        IF(NDIMD.EQ.1.OR.NDIMD.EQ.2) THEN
          ALLOCATE(MAT1(NLOC,NLOC))
          ALLOCATE(FLOC(NLOC))
          DO IDIMD=1,NDIMD
            MAT1=ATOMSET(IAT)%DIAGSLOC(:,:,IDIMD)
!           == COMPLEX DIAGONALIZATION USED. (CAN PROBABLY BE REAL)
            CALL LIB$DIAGC8(NLOC,MAT1,FLOC,ATOMSET(IAT)%DIAGSLOC(:,:,IDIMD))
            DO I=1,NCHI
              ATOMSET(IAT)%DIAGSLOC(:,I,IDIMD) &
      &                   =ATOMSET(IAT)%DIAGSLOC(:,I,IDIMD)/SQRT(FLOC(I))
            ENDDO
          ENDDO
          DEALLOCATE(MAT1)
          DEALLOCATE(FLOC)
        ELSE IF(NDIMD.EQ.4) THEN  ! NON-COLLINEAR CASE
          ALLOCATE(MAT1(2*NLOC,2*NLOC))
          ALLOCATE(MAT2(2*NLOC,2*NLOC))
          ALLOCATE(FLOC(2*NLOC))  ! (COMPLEX ARRAY)
          MAT1(1:NLOC ,1:NLOC )=ATOMSET(IAT)%DIAGSLOC(:,:,1)      
          MAT1(1:NLOC ,NLOC+1:)=ATOMSET(IAT)%DIAGSLOC(:,:,2)      
          MAT1(NLOC+1:,1:NLOC )=ATOMSET(IAT)%DIAGSLOC(:,:,3)      
          MAT1(NLOC+1:,NLOC+1:)=ATOMSET(IAT)%DIAGSLOC(:,:,4)      
          CALL LIB$DIAGC8(2*NLOC,MAT1,FLOC,MAT2)
          DO I=1,2*NLOC
            MAT2(:,I)=MAT2(:,I)/SQRT(FLOC(I))
          ENDDO
          ATOMSET(IAT)%DIAGSLOC(:,:,1)=MAT2(1:NLOC ,1:NLOC )
          ATOMSET(IAT)%DIAGSLOC(:,:,2)=MAT2(1:NLOC ,NLOC+1:)      
          ATOMSET(IAT)%DIAGSLOC(:,:,3)=MAT2(NLOC+1:,1:NLOC )      
          ATOMSET(IAT)%DIAGSLOC(:,:,4)=MAT2(NLOC+1:,NLOC+1:)
          DEALLOCATE(MAT1)
          DEALLOCATE(MAT2)
          DEALLOCATE(FLOC)
        ELSE
          CALL ERROR$MSG('ILLEGAL VALUE OF NDIMD')
          CALL ERROR$STOP('DMFT_DIAGSLOC')
        END IF
!
!       == CONVERT FROM (UP-DOWN REPRESENTATION TO (TXYZ) ======================
        CALL SPINOR$CONVERT('FWRD',NLOC,NDIMD,ATOMSET(IAT)%DIAGSLOC)
      ENDDO  !END OF LOOP OVER ATOMS
!
!     ==========================================================================
!     == OPTIONAL TESTS
!     ==========================================================================
      IF(TTEST) THEN
        DO IAT=1,NAT
          NLOC=ATOMSET(IAT)%NLOC
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          ALLOCATE(MAT1(NLOC,NLOC))
          ALLOCATE(MAT2(NLOC,NLOC))
          DO IDIMD=1,NDIMD
            MAT1=(0.D0,0.D0)
            DO IKPT=1,NKPTL
              MAT1=MAT1+KSET(IKPT)%WKPTL*KSET(IKPT)%SINV(I1:I2,I1:I2,IDIMD)
            END DO
            IF(NDIMD.LT.4)MAT1=REAL(MAT1)
!
            MAT2=ATOMSET(IAT)%DIAGSLOC(:,:,IDIMD)
            MAT1=MATMUL(CONJG(TRANSPOSE(MAT2)),MATMUL(MAT1,MAT2))
!
            WRITE(*,FMT='(82("="),T10,"  TEST DIAGSLOC FOR IDIMD=",I2,"  ")') &
     &                                            IDIMD
            DO I=1,NLOC
              WRITE(*,FMT='("TEST",10("(",2F10.5,")"))')MAT1(I,:)
            ENDDO
          ENDDO ! END OF LOOP OVER IDIMD
        ENDDO   !END OF LOOP OVER ATOMS
        CALL ERROR$MSG('SCHEDULED STOP AFTER TEST')
        CALL ERROR$STOP('DMFT_DIAGSLOC')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$COLLECTHAMILTONIAN_WITHKSET()
!     **************************************************************************
!     ** COLLECTS THE HAMILTONIAN AND STORES IT ON THE MODULE                 **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,NDIMD,IPROOFCHI,KBT &
     &                      ,NAT,KSET,ATOMSET
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      COMPLEX(8),ALLOCATABLE :: RHO(:,:,:) !(NCHI,NCHI,NDIMD)
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IPRO,IB,J,IAT,IDIM1,IDIM2
      INTEGER(4)             :: IDIMD
      INTEGER(4)             :: I1,I2
      REAL(8)                :: F(NB,NKPTL,NSPIN)
      REAL(8)                :: EMIN,EMAX
      COMPLEX(8)             :: MAT(NCHI,NCHI)
      COMPLEX(8)             :: CSVAR
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                      CALL TRACE$PUSH('DMFT$COLLECTHAMILTONIAN')
!
!     ==========================================================================
!     ==  DETERMINE ERHO 
!     ==========================================================================
      CALL DMFT$COLLECTOCCUPATIONS(NKPTL,NSPIN,NB,F)
      EMIN=+HUGE(EMIN)
      EMAX=-HUGE(EMAX)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          DO IB=1,NB
            CALL DMFT_EOFF(KBT,F(IB,IKPT,ISPIN),KSET(IKPT)%ERHO(IB,ISPIN))
            EMIN=MIN(EMIN,KSET(IKPT)%ERHO(IB,ISPIN))
            EMAX=MAX(EMAX,KSET(IKPT)%ERHO(IB,ISPIN))
          ENDDO
        ENDDO
      ENDDO 
PRINT*,'BOUNDS OF SPECTRUM [EV] ',EMIN*27.211D0,EMAX*27.211D0
!
!     ==========================================================================
!     ==  EXTRACT <PSI|PSI>
!     ==========================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          IF(THIS%NB.NE.NB) THEN
            CALL ERROR$MSG('INCONSISTENT NUMBER OF STATES IN WAVES AND DYNOCC')
            CALL ERROR$I4VAL('NB IN DYNOCC',NB)
            CALL ERROR$I4VAL('NB IN WAVES ',THIS%NB)
            CALL ERROR$STOP('DMFT$COLLECTHAMILTONIAN')
          END IF
          NBH=THIS%NBH
!
!         ======================================================================
!         ==  DETERMINE ORBITAL PROJECTIONS                                   ==
!         ======================================================================
!         == SPECIFY IF SECOND K-POINT IS OBTAINED FROM TIME-INVERSION.
!         == NON-COLLINEAR: NO, BECAUSE THERE IS NO TIME INVERSION SYMMETRY
!         == TINV: NO, K-POINT FALLS ONTO ITSELF UNDER TIME-INVERSION SYMMETRY.
          KSET(IKPT)%TADDMINUSK=(NDIMD.NE.4).AND.(.NOT.GSET%TINV)
!
          DO ICHI=1,NCHI
            IPRO=IPROOFCHI(ICHI)
            DO IBH=1,NBH
              IF(NBH.NE.NB) THEN
                KSET(IKPT)%PIPSI(:,ICHI,2*IBH-1,ISPIN) &
     &                                              = REAL(THIS%TBC(:,IBH,IPRO))
                KSET(IKPT)%PIPSI(:,ICHI,2*IBH  ,ISPIN)& 
     &                                              =AIMAG(THIS%TBC(:,IBH,IPRO))
              ELSE
                KSET(IKPT)%PIPSI(:,ICHI,IBH,ISPIN)=THIS%TBC(:,IBH,IPRO)
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==   RHO(A,S1,B,S2,K)=SUM_N <PI(A,S1)|PSI(N,K)>F_N <PSI(N,K)|PI(B,S2)>  ==
!     ==                                                                      ==
!     ==   NON-SPIN POLARIZED:                                                ==
!     ==     RHO(S1,S2)=1/2*[ RHO(1)*(1,0/0,1) ]                              ==
!     ==   COLLINEAR SPIN-POLARIZED                                           ==
!     ==     RHO(S1,S2)=1/2*[ RHO(1)*(1,0/0,1)+RHO(2)*(1,0/0,-1) ]            ==
!     ==   NON-COLLINEAR SPIN-POLARIZED                                       ==
!     ==     RHO(S1,S2)=1/2*[ RHO(1)*(1,0/0,1)+RHO(2)*(0,1/1,0)               ==
!     ==                      +RHO(3)*(0,-I/I,0)+RHO(2)*(1,0/0,-1) ]          ==
!     ==   WHERE                                                              ==
!     ==     RHO(1)=   RHO(1,1)+RHO(2,2)                                      ==
!     ==     RHO(2)=   RHO(2,1)+RHO(1,2)                                      ==
!     ==     RHO(3)=-I(RHO(2,1)-RHO(1,2))                                     ==
!     ==     RHO(4)=   RHO(1,1)-RHO(2,2)  !FOR NSPIN=2 THIS IS SIMV(2)        ==
!     ==                                                                      ==
!     ==  CAUTION! MATRICES ARE ARRANGED IN A FORTRAN-STYLE NOTATION WITH THE ==
!     ==     FIRST INDEX RUNNING FIRST. THE NOTATION (A,B/C,D) HOWEVER IS A   ==
!     ==     LATEX STYLE NOTATION WITH THE SECOND INDEX RUNNING FIRST!        ==
!     ==                                                                      ==
!     ==========================================================================
!     ==  NDIM=1,NSPIN=1  NDIMD=1: (T)                                        ==
!     ==  NDIM=1,NSPIN=2  NDIMD=2: (T,Z)                                      ==
!     ==  NDIM=2,NSPIN=1  NDIMD=4: (T,X,Y,Z)                                  ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        KSET(IKPT)%RHO=(0.D0,0.D0)
        DO ISPIN=1,NSPIN
          DO IDIM1=1,NDIM
            DO IDIM2=1,NDIM
              MAT(:,:)=(0.D0,0.D0)
              DO IB=1,NB
                DO J=1,NCHI
                  CSVAR=F(IB,IKPT,ISPIN) &
       &               *CONJG(KSET(IKPT)%PIPSI(IDIM2,J,IB,ISPIN))
                  MAT(:,J)=MAT(:,J) &
       &                  +KSET(IKPT)%PIPSI(IDIM1,:,IB,ISPIN)*CSVAR
                ENDDO
              ENDDO
!
!             ==================================================================
!             == DISTRIBUTE DENSITY MATRIX ONTO VARIOUS SPIN COMPONENTS       ==
!             == USING THE UP-DOWN NOTATION (1=UU,2=DU,3=UD,4=DD)             ==
!             == NSPIN=1,2: IDIM1=IDIM2=NDIM=1 =>IDIMD=ISPIN                  ==
!             == NSPIN=1,NDIM=2,ISPIN=1: (1,2,3,4)                            ==
!             ==================================================================
              IDIMD=IDIM1+NDIM*(IDIM2-1)+ISPIN-1
              KSET(IKPT)%RHO(:,:,IDIMD)=MAT
            ENDDO
          ENDDO
        ENDDO
!
!       ========================================================================
!       == CONVERT FROM UP-DOWN REPRESENTATION INTO (TXYZ)                    ==
!       ========================================================================
        CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%RHO)
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE LOCAL DENSITY MATRIX FOR TESTING 
!     ==========================================================================
      IF(TTEST) THEN
!       == CALCULATE DENSITY MATRIX FOR TESTING ================================
        ALLOCATE(RHO(NCHI,NCHI,NDIMD))
        RHO(:,:,:)=(0.D0,0.D0)
        DO IKPT=1,NKPTL
          RHO=RHO+KSET(IKPT)%WKPTL*KSET(IKPT)%RHO
        ENDDO
!       == MAKE HERMITEAN ======================================================
        DO IDIMD=1,NDIMD
          RHO(:,:,IDIMD)=0.5D0*(RHO(:,:,IDIMD)+CONJG(TRANSPOSE(RHO(:,:,IDIMD))))
        ENDDO
!       == PRINT ===============================================================
        PRINT*,' IN DMFT$COLLECTHAMILTONIAN_WITHKSET'
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          CALL SPINOR$PRINTMATRIX(6,'DENSITY MATRIX(UPDOWN)','UPDOWN' &
      &                          ,I1,I2,NDIMD,NCHI,RHO)
        ENDDO
        DEALLOCATE(RHO)
      END IF
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$COLLECTFULLDENMAT()
!     **************************************************************************
!     ==  EXTRACT ONSITE DENSITY MATRICES WITH ALL PROJECTOR FUNCTIONS        ==
!     ==  FOR DOUBLE COUNTING                                                 ==
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NKPTL,NSPIN,NAT,NDIM,NDIMD,NB,ATOMSET
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: I1,I2
      COMPLEX(8)             :: CSVAR
      INTEGER(4)             :: NB_,NSPIN_
      INTEGER(4)             :: IAT,IKPT,ISPIN,IPRO,IBH,IDIM1,IDIM2,IDIMD,I
      INTEGER(4)             :: LMN
      REAL(8)                :: OCC(NB,NKPTL,NSPIN)
!     **************************************************************************
!
!     ==========================================================================
!     == COLLECT OCCUPATIONS                                                  ==
!     ==========================================================================
      CALL DYNOCC$GETI4('NB',NB_)
      IF(NB_.NE.NB) THEN
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$I4VAL('NB',NB)
        CALL ERROR$I4VAL('NB_',NB_)
        CALL ERROR$STOP('DMFT$COLLECTOCCUPATIONS')
      END IF
!
      CALL DYNOCC$GETI4('NSPIN',NSPIN_)
      IF(NSPIN_.NE.NSPIN) THEN
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$I4VAL('NSPIN_',NSPIN_)
        CALL ERROR$STOP('DMFT$COLLECTOCCUPATIONS')
      END IF
      CALL WAVES_DYNOCCGETR8A('OCC',NB*NKPTL*NSPIN,OCC)
!
!     ==========================================================================
!     == SET DENSITY MATRIX TO ZERO
!     ==========================================================================
      DO IAT=1,NAT
        ATOMSET(IAT)%DENMAT%MAT=(0.D0,0.D0)
      ENDDO        
!
!     ==========================================================================
!     == ADD UP DENSITY MATRIX                                                ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          IF(THIS%NB.NE.NB) THEN
            CALL ERROR$MSG('INCONSISTENT NUMBER OF STATES IN WAVES AND DYNOCC')
            CALL ERROR$I4VAL('NB IN DYNOCC',NB)
            CALL ERROR$I4VAL('NB IN WAVES ',THIS%NB)
            CALL ERROR$STOP('DMFT$COLLECTHAMILTONIAN')
          END IF
          NBH=THIS%NBH
!
!         ======================================================================
!         ==  DETERMINE ORBITAL PROJECTIONS                                   ==
!         ======================================================================
          IPRO=0
          DO IAT=1,NAT
            LMNX=ATOMSET(IAT)%DENMAT%LMNX
            IF(ATOMSET(IAT)%NLOC.LE.0) THEN
              IPRO=IPRO+LMNX
              CYCLE
            END IF   
            I1=IPRO+1         
            I2=IPRO+LMNX
            DO IBH=1,NBH
              IF(NBH.NE.NB) THEN
!               == THIS IS FOR GENERAL SPECIAL K-POINTS ========================
!               == TBC IS FROM A SUPER WAVE FUNCTION AND CONTAINS TWO WAVE =====
!               == FUNCTIONS.
                DO IDIM2=1,NDIM
                  DO IDIM1=1,NDIM
                    IDIMD=IDIM1+NDIM*(IDIM2-1)+(ISPIN-1)
                    DO LMN=1,LMNX
                      I=IPRO+LMN
                      CSVAR=REAL(THIS%TBC(IDIM1,IBH,I))*OCC(2*IBH-1,IKPT,ISPIN)
                      ATOMSET(IAT)%DENMAT%MAT(:,LMN,IDIMD) &
     &                               =ATOMSET(IAT)%DENMAT%MAT(:,LMN,IDIMD) &
     &                               +REAL(THIS%TBC(IDIM2,IBH,I1:I2))*CSVAR
!
                      CSVAR=AIMAG(THIS%TBC(IDIM1,IBH,I))*OCC(2*IBH,IKPT,ISPIN)
                      ATOMSET(IAT)%DENMAT%MAT(:,LMN,IDIMD) &
     &                             =ATOMSET(IAT)%DENMAT%MAT(:,LMN,IDIMD) &
     &                             +AIMAG(THIS%TBC(IDIM2,IBH,I1:I2))*CSVAR
                    ENDDO
                  ENDDO
                ENDDO
              ELSE
!               == THIS IS FOR GENERAL K-POINTS ================================
                DO IDIM2=1,NDIM
                  DO IDIM1=1,NDIM
                    IDIMD=IDIM1+NDIM*(IDIM2-1)+(ISPIN-1)
                    DO LMN=1,LMNX
                      I=IPRO+LMN
                      CSVAR=CONJG(THIS%TBC(IDIM1,IBH,I))*OCC(IBH,IKPT,ISPIN)
                      ATOMSET(IAT)%DENMAT%MAT(:,LMN,IDIMD) &
     &                             =ATOMSET(IAT)%DENMAT%MAT(:,LMN,IDIMD) &
     &                             +THIS%TBC(IDIM2,IBH,I1:I2)*CSVAR
                    ENDDO
                  ENDDO
                ENDDO
              END IF
            ENDDO  ! END OF LOOP OVER BANDS
            IPRO=IPRO+LMNX
          ENDDO ! END OF LOOP OVER ATOMS
        ENDDO ! END OF LOOP OVER SPIN
      ENDDO ! END OF LOOP OVER K-POINTS
!
!     ==========================================================================
!     ==  INCLUDE CONTRIBUTION FROM -K FOR NON-SPINPOLARIZED AND COLLINEAR    ==
!     ==  FOR TIME INVERSION W/O SPIN PSI(-K)=CONJG(PSI(+K))                  ==
!     ==========================================================================
      IF(NDIMD.LT.4) THEN
        DO IAT=1,NAT
          DO IDIMD=1,NDIMD
            ATOMSET(IAT)%DENMAT%MAT(:,:,IDIMD) &
     &                           =REAL(ATOMSET(IAT)%DENMAT%MAT(:,:,IDIMD))
          ENDDO        
        ENDDO   
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM TO (T,X,Y,Z) REPRESENTATION AND MAKE HERMITEAN            ==
!     ==========================================================================
      DO IAT=1,NAT
        IF(ATOMSET(IAT)%NLOC.LE.0) CYCLE
        LMNX=ATOMSET(IAT)%DENMAT%LMNX
        CALL SPINOR$CONVERT('FWRD',LMNX,NDIMD,ATOMSET(IAT)%DENMAT%MAT)
        DO IDIMD=1,NDIMD
          ATOMSET(IAT)%DENMAT%MAT(:,:,IDIMD) &
    &              =0.5D0*(ATOMSET(IAT)%DENMAT%MAT(:,:,IDIMD) &
    &                     +CONJG(TRANSPOSE(ATOMSET(IAT)%DENMAT%MAT(:,:,IDIMD))))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  PRINT
!     ==========================================================================
      IF(TPRINT) THEN
        DO IAT=1,NAT
          IF(ATOMSET(IAT)%NLOC.LE.0) CYCLE
          WRITE(*,FMT='(82("="),T10," DENSITY MATRIX FOR ATOM ",I3," ")')IAT
          DO IDIMD=1,NDIMD
            DO LMN=1,LMNX
              WRITE(*,FMT='("IDIMD=",I1,":",100("(",2F10.5,")"))')IDIMD &
       &              ,ATOMSET(IAT)%DENMAT%MAT(LMN,:,IDIMD)
            ENDDO
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$COLLECTOCCUPATIONS(NKPTL,NSPIN,NB,F)
!     **************************************************************************
!     ** COLLECTS THE OCCUPATIONS FOR THIS NODE AND REMOVES THE K-POINT       **
!     ** WEIGHT                                                               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NKPTL
      INTEGER(4),INTENT(IN)  :: NSPIN
      INTEGER(4),INTENT(IN)  :: NB
      REAL(8)   ,INTENT(OUT) :: F(NB,NKPTL,NSPIN)
      INTEGER(4)             :: NB_,NSPIN_
      REAL(8)                :: WKPT(NKPTL)
      INTEGER(4)             :: IKPT,ISPIN
!     **************************************************************************
      CALL DYNOCC$GETI4('NB',NB_)
      IF(NB_.NE.NB) THEN
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$I4VAL('NB',NB)
        CALL ERROR$I4VAL('NB_',NB_)
        CALL ERROR$STOP('DMFT$COLLECTOCCUPATIONS')
      END IF
!
      CALL DYNOCC$GETI4('NSPIN',NSPIN_)
      IF(NSPIN_.NE.NSPIN) THEN
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$I4VAL('NSPIN_',NSPIN_)
        CALL ERROR$STOP('DMFT$COLLECTOCCUPATIONS')
      END IF
!
      CALL WAVES_DYNOCCGETR8A('OCC',NB*NKPTL*NSPIN,F)
      CALL WAVES_DYNOCCGETR8A('WKPT',NKPTL,WKPT)
!
!     == MULTIPLY WITH SPIN MULTIPLICITY =======================================
      IF(NSPIN.EQ.1) WKPT=2.D0*WKPT
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          F(:,IKPT,ISPIN)=F(:,IKPT,ISPIN)/WKPT(IKPT)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE DMFT_EOFF(KBT,F,E)
!     **************************************************************************
!     ** CALCULATES ENERGY FROM F=1/(1+EXP(E/KBT))                            **
!     ** SMALL AND LARGE OCCUPATIONS ARE SIMPLIFIED BY LINEAR EXTRAPOLATION   **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: KBT
      REAL(8),INTENT(IN)  :: F
      REAL(8),INTENT(OUT) :: E
      REAL(8),PARAMETER   :: FMIN=1.D-6
      REAL(8)             :: F0,SLOPE,VAL
!     **************************************************************************
      IF(F.GT.FMIN.AND.F.LT.1.D0-FMIN) THEN
        E=LOG((1.D0-F)/F)
      ELSE 
        F0=FMIN
        VAL=LOG( (1.D0-F0)/F0 )
        SLOPE=-1.D0/( F0*(1.D0-F0) )
        IF(F.GT.0.5D0) THEN
          VAL=-VAL
          F0=1.D0-F0
        END IF
        E=VAL+SLOPE*(F-F0)
      END IF
      E=E*KBT
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$SOLVER_WITHKSET()
!     **************************************************************************
!     ** MIMICKS A DMFT SOLVER                                                **
!     **  E=1/BETA * SUM_NU E(IOMEGA_NU0+) TR[ DV * GLOC^\DAGGER]             **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NAT,NDIMD,NOMEGA,OMEGA,KBT,ATOMSET
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      REAL(8)   ,PARAMETER   :: DMFTHFWEIGHT=0.25D0
      REAL(8)                :: EV
      REAL(8)                :: PHILW
      REAL(8)                :: ETOT
      REAL(8)                :: EDC
      INTEGER(4)             :: NLOC  !#(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4)             :: LMNX  !#(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4)             :: IAT
      INTEGER(4)             :: IDIMD,I,J,LMN1,LMN2,NU
      COMPLEX(8),ALLOCATABLE :: H(:,:,:)
!     **************************************************************************
PRINT*,'ENTERING DMFT$SOLVER...'
      CALL CONSTANTS('EV',EV)
!
!     ==========================================================================
!     == LOOP OVER ATOMS                                                      ==
!     ==========================================================================
      ETOT=0.D0
      DO IAT=1,NAT
        NLOC=ATOMSET(IAT)%NLOC
        IF(NLOC.LE.0) CYCLE
!
!       ========================================================================
!       == COLLECT LOCAL ARRAYS GLOC1 AND ULOC1                               ==
!       ========================================================================
PRINT*,'CONSTRUCTING SELF ENERGY FOR ATOM ',IAT,' ...'
        CALL DMFT$HFSOLVER_WITHSET(NLOC,NDIMD,NOMEGA,KBT,OMEGA &
     &                            ,ATOMSET(IAT)%GLOC,ATOMSET(IAT)%GLOCLAUR &
     &                            ,ATOMSET(IAT)%U,PHILW &
     &                            ,ATOMSET(IAT)%SLOC,ATOMSET(IAT)%SLOCLAUR)
        ETOT=ETOT+PHILW
!!$WRITE(*,FMT='(82("="),T10," SLOC FROM HFSOLVER_WITHKSET ",I4,"  ")'),IAT
!!$DO NU=1,NOMEGA
!!$  DO IDIMD=1,NDIMD
!!$    DO I=1,NLOC
!!$      WRITE(*,FMT='(4I5,100F10.5)')IAT,NU,IDIMD,I &
!!$                              ,ATOMSET(IAT)%SLOC(I,:,IDIMD,NU)
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
!!$WRITE(*,FMT='(82("="),T10," SLOCLAUR FROM HFSOLVER_WITHKSET ",I4,"  ")'),IAT
!!$DO NU=1,3
!!$  DO IDIMD=1,NDIMD
!!$    DO I=1,NLOC
!!$      WRITE(*,FMT='(4I5,100F10.5)')IAT,NU,IDIMD,I &
!!$                              ,ATOMSET(IAT)%SLOCLAUR(I,:,IDIMD,NU)
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
!
!       ========================================================================
!       == SUBTRACT DOUBLE COUNTING TERM                                      ==
!       ========================================================================
        LMNX=ATOMSET(IAT)%DENMAT%LMNX
        CALL DMFT_DC(IAT,LMNX,NDIMD,ATOMSET(IAT)%DENMAT%MAT &
       &                       ,EDC,ATOMSET(IAT)%DENMAT%H)
        ETOT=ETOT+EDC
!
!       ========================================================================
!       == SHIFT STATIC SELF ENERGY TO THE DC TERM                            ==
!       ========================================================================
        NLOC=ATOMSET(IAT)%NLOC
        ALLOCATE(H(NLOC,NLOC,NDIMD))
!       == SHIFT ONLY REAL PART BECAUSE ATOMSET%DENMAT%H IS REAL ===============
        H=ATOMSET(IAT)%SLOCLAUR(:,:,:,1)
!
!       == REMOVE FROM SLOC ====================================================
        DO NU=1,NOMEGA
          ATOMSET(IAT)%SLOC(:,:,:,NU)=ATOMSET(IAT)%SLOC(:,:,:,NU)-H
        ENDDO
        ATOMSET(IAT)%SLOCLAUR(:,:,:,1)=ATOMSET(IAT)%SLOCLAUR(:,:,:,1)-H
!
!       == ADD TO DOUBLE COUNTING HAMILTONIAN===================================
        DO I=1,NLOC
          LMN1=ATOMSET(IAT)%DENMAT%LMN(I)
          DO J=1,NLOC
            LMN2=ATOMSET(IAT)%DENMAT%LMN(J)
            ATOMSET(IAT)%DENMAT%H(LMN1,LMN2,:) &
   &                     =ATOMSET(IAT)%DENMAT%H(LMN1,LMN2,:)+H(I,J,:)
          ENDDO
        ENDDO
        DEALLOCATE(H)
!
!       ========================================================================
!       == SCALE                                                              ==
!       ========================================================================
        ATOMSET(IAT)%SLOC    =DMFTHFWEIGHT*ATOMSET(IAT)%SLOC
        ATOMSET(IAT)%SLOCLAUR=DMFTHFWEIGHT*ATOMSET(IAT)%SLOCLAUR
        ATOMSET(IAT)%DENMAT%H=DMFTHFWEIGHT*ATOMSET(IAT)%DENMAT%H
!
!       ========================================================================
!       == PRINT                                                              ==
!       ========================================================================
        IF(TPRINT) THEN
          WRITE(*,FMT='(82("="),T10,A," IAT=",I4,"  ")') &
     &                     " HAMILTONIAN FROM HFSOLVER_WITHKSET ",IAT
          DO NU=1,NOMEGA
            DO IDIMD=1,NDIMD
              DO I=1,ATOMSET(IAT)%NLOC
                WRITE(*,FMT='(4I5,100("(",2F10.5,")"))')IAT,NU,IDIMD,I &
     &                              ,ATOMSET(IAT)%SLOC(I,:,IDIMD,NU)
              ENDDO
            ENDDO
          ENDDO
        END IF
!
      ENDDO ! END OF LOOP OVER ATOMS (IAT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_DC(IAT,LMNX,NDIMD,RHO,EDC,HAM)
!     **************************************************************************
!     ** CALCULATES THE DFT EXCHANGE ENERGY FOR THE CORRELATED ORBITALS.      **
!     ** EDC IS TO BE ADDED! TO THE TOTAL ENERGY (MINUS SIGN IS INCLUDED)     **
!     **                                                                      **
!     ** REMARK: SOME ROUTINES ONLY WORK IN THE NON-COLLINEAR MODE            **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      INTEGER(4),INTENT(IN)  :: LMNX          !  #(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4),INTENT(IN)  :: NDIMD
      REAL(8)   ,INTENT(OUT) :: EDC           ! -E_(DFT-EXCHANGE)
      COMPLEX(8),INTENT(IN)  :: RHO(LMNX,LMNX,NDIMD) ! DENSITY MATRIX
      COMPLEX(8),INTENT(OUT) :: HAM(LMNX,LMNX,NDIMD) ! HAMILTONIAN CONTRIB.
      REAL(8)                :: D(LMNX,LMNX,4)
      REAL(8)                :: H(LMNX,LMNX,4)
      REAL(8)                :: HALL(LMNX,LMNX,4)
      REAL(8)   ,ALLOCATABLE :: DT(:,:,:)    !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: DTALL(:,:,:) !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: HT(:,:,:)    !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: HTALL(:,:,:) !(LMNXT,LMNXT,NDIMD)
      REAL(8)                :: EX
      INTEGER(4)             :: IDFTTYPE
      INTEGER(4)             :: LMNXT   ! #(VALENCE+SCATTERING WAVES)
      INTEGER(4)             :: LNXT    ! #(VALENCE+SCATTERING WAVES W/O M)
      INTEGER(4),ALLOCATABLE :: LOXT(:)  
      INTEGER(4)             :: LRX
      INTEGER(4)             :: LMRX    ! #(SPHERICAL HARMONICS IN THE DENSITY)
      INTEGER(4)             :: GID     ! GRID ID
      INTEGER(4)             :: NR      ! #(GRID POINTS)
      REAL(8)   ,ALLOCATABLE :: AECORE(:) !(NR) CORE DENSITY
      INTEGER(4)             :: ISP
!     **************************************************************************
      HAM=0.D0
      CALL DFT$GETI4('TYPE',IDFTTYPE)
      PRINT*,'IDFTTYPE ',IDFTTYPE
      IF(IDFTTYPE.EQ.5002) RETURN
!
!     ==========================================================================
!     == COLLECT DATA                                                         ==
!     ==========================================================================
      ISP=ISPECIES(IAT)
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(AECORE(NR))
      CALL SETUP$GETR8A('AECORE',NR,AECORE)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
      LMNXT=POTPAR(ISP)%TAILED%LMNX  ! SIZE OF U-TENSOR ON POTPAR
      LNXT=POTPAR(ISP)%TAILED%LNX  
      ALLOCATE(LOXT(LNXT))
      LOXT=POTPAR(ISP)%TAILED%LOX  
!
!     ==========================================================================
!     ==  MAP INTO NON-COLLINEAR ARRAYS, KEEP ONLY REAL PART                  ==
!     ==========================================================================
      D(:,:,:)=0.D0
      D(:,:,1)=REAL(RHO(:,:,1))
      IF(NDIMD.EQ.2) THEN 
        D(:,:,4)=REAL(RHO(:,:,2))
      ELSE IF(NDIMD.EQ.4) THEN
        D(:,:,2:4)=REAL(RHO(:,:,2:4))
      END IF
!
!     ==========================================================================
!     ==  MAP ONTO REAL ARRAY. (T,X,Y,Z) REPRESENTATION                       ==
!     ==========================================================================
      ALLOCATE(DTALL(LMNXT,LMNXT,4))
      ALLOCATE(DT(LMNXT,LMNXT,4))
      ALLOCATE(HTALL(LMNXT,LMNXT,4))
      ALLOCATE(HT(LMNXT,LMNXT,4))
      CALL LMTO_BLOWUPDENMATNL(IAT,IAT,4,LMNX,LMNX,D,LMNXT,LMNXT,DT)
      POTPAR(ISP)%TALLORB=.TRUE.
      CALL LMTO_BLOWUPDENMATNL(IAT,IAT,4,LMNX,LMNX,D,LMNXT,LMNXT,DTALL)
      CALL LMTO_SIMPLEDC(GID,NR,LMNXT,LNXT,LOXT,POTPAR(ISP)%TAILED%AEF &
     &                  ,LRX,AECORE,DT,DTALL,EX,HT,HTALL)
      CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,4,LMNXT,LMNXT,HTALL,LMNX,LMNX,HALL)
      POTPAR(ISP)%TALLORB=.FALSE.  ! DO NOT FORGET THIS!!!!!
      CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,4,LMNXT,LMNXT,HT,LMNX,LMNX,H)
      EDC=-EX
      H=-(H+HALL)
PRINT*,'DOUBLE COUNTING CORRECTION ENERGY FOR ATOM=',IAT,-EX
!
!     ==========================================================================
!     ==  MAP BACK FROM  NON-COLLINEAR ARRAYS                                 ==
!     ==========================================================================
      HAM=0.D0
      HAM(:,:,1)=CMPLX(H(:,:,1))
      IF(NDIMD.EQ.2) THEN
         HAM(:,:,2)=CMPLX(H(:,:,4))
      ELSE IF(NDIMD.EQ.4) THEN
         HAM(:,:,2:4)=CMPLX(H(:,:,2:4))
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_UCHI_WITHATOMSET()
!     **************************************************************************
!     ** CALCULATE THE U-TENSOR OF THE CORRELATED ORBITALS IN THE SELECTED SET
!     **************************************************************************
      USE LMTO_MODULE, ONLY: ISPECIES,LNX,LOX,POTPAR
      USE DMFT_MODULE, ONLY: NAT,ATOMSET
      IMPLICIT NONE
      REAL(8)  ,ALLOCATABLE :: UB(:,:,:,:)
      INTEGER(4)            :: ISP ! ATOM TYPE
      INTEGER(4)            :: L,LN,LMN,ICHI
      INTEGER(4)            :: LMNX
      INTEGER(4)            :: NCHI1        ! #(SELECTED ORBITALS ON THIS ATOM)
      INTEGER(4)            :: IAT
!     **************************************************************************
                                          CALL TRACE$PUSH('DMFT_UCHI')
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
!
!       == CALCULATE LOCAL DIMENSIONS LMNX AND NCHI1 ===========================
        LMN=0
        ICHI=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(POTPAR(ISP)%TORB(LN)) ICHI=ICHI+2*L+1
          LMN=LMN+2*L+1
        ENDDO
        LMNX=LMN
        NCHI1=ICHI
        IF(NCHI1.EQ.0) CYCLE
        IF(NCHI1.NE.ATOMSET(IAT)%NLOC) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
          CALL ERROR$STOP('DMFT_UCHI')
        END IF
!
!       ========================================================================
!       ==  CALCULATE U-TENSOR IN THE BASIS OF LOCAL ORBITALS IGNORING TORB   ==
!       ========================================================================
        ALLOCATE(UB(LMNX,LMNX,LMNX,LMNX))
        CALL DMFT_ULOCAL(IAT,LMNX,UB)
!
!       ========================================================================
!       ==  THROW OUT ELEMENTS THAT ARE NOT USED                              ==
!       ========================================================================
        LMN=0
        ICHI=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(POTPAR(ISP)%TORB(LN)) THEN
            UB(:,:,:,ICHI+1:ICHI+2*L+1)=UB(:,:,:,LMN+1:LMN+2*L+1)
            ICHI=ICHI+2*L+1
          END IF
          LMN=LMN+2*L+1
        ENDDO
!  
!       ========================================================================
        LMN=0
        ICHI=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(POTPAR(ISP)%TORB(LN)) THEN
            UB(:,:,ICHI+1:ICHI+2*L+1,:NCHI1)=UB(:,:,LMN+1:LMN+2*L+1,:NCHI1)
            ICHI=ICHI+2*L+1
          END IF
          LMN=LMN+2*L+1
        ENDDO
!  
!       ========================================================================
        LMN=0
        ICHI=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(POTPAR(ISP)%TORB(LN)) THEN
            UB(:,ICHI+1:ICHI+2*L+1,:NCHI1,:NCHI1) &
      &                                     =UB(:,LMN+1:LMN+2*L+1,:NCHI1,:NCHI1)
            ICHI=ICHI+2*L+1
          END IF
          LMN=LMN+2*L+1
        ENDDO
!  
!       ========================================================================
        LMN=0
        ICHI=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(POTPAR(ISP)%TORB(LN)) THEN
            UB(ICHI+1:ICHI+2*L+1,:NCHI1,:NCHI1,:NCHI1) &
       &                               =UB(LMN+1:LMN+2*L+1,:NCHI1,:NCHI1,:NCHI1)
            ICHI=ICHI+2*L+1
          END IF
          LMN=LMN+2*L+1
        ENDDO  
!
!       ========================================================================
!       ==  PUT ONTO ATOMSET                                                  ==
!       ========================================================================
        ATOMSET(IAT)%U=UB(:NCHI1,:NCHI1,:NCHI1,:NCHI1)
        DEALLOCATE(UB)
      ENDDO  ! END LOOP OVER ATOMS
                                          CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_ULOCAL(IAT,LMNX,U)
!     **************************************************************************
!     **  CALCULATES THE U-TENSOR OF THE NATURAL TIGHT-BINDING ORBITALS       **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR,SBAR,SBARLI1,LNX,LOX
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      INTEGER(4),INTENT(IN)  :: LMNX          !  #(LOCAL ORBITALS ON THIS SITE)
      REAL(8)   ,INTENT(OUT) :: U(LMNX,LMNX,LMNX,LMNX)
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      REAL(8)   ,ALLOCATABLE :: UT(:,:,:,:)
      INTEGER(4)             :: ISP     ! ATOM TYPE
      INTEGER(4)             :: LMNXT   ! #(VALENCE+SCATTERING WAVES)
      INTEGER(4)             :: LMX     ! #(SCATTERING WAVES)
      INTEGER(4)             :: NNS     ! #(PAIRS FOR STRUCTURE CONSTANTS)
      REAL(8)                :: SVAR
      INTEGER(4)             :: LMN,L,LM,LN,IM,NN,NN0,J
!     **************************************************************************
      ISP=ISPECIES(IAT)
      LMNXT=POTPAR(ISP)%TAILED%LMNX  ! SIZE OF U-TENSOR ON POTPAR
!
!     ==========================================================================
!     == COLLECT U-TENSOR IN EXPANDED BASIS                                   ==
!     ==========================================================================
      ALLOCATE(UT(LMNXT,LMNXT,LMNXT,LMNXT))
      UT=POTPAR(ISP)%TAILED%U
!
!     ==========================================================================
!     == TEST CONSISTENCY OF ARRAY DIMENSIONS                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          LMN=LMN+2*L+1 
        ENDDO
        IF(LMN.NE.LMNX) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
          CALL ERROR$I4VAL('CALCULATED LMNX',LMN)
          CALL ERROR$I4VAL('INPUT LMNX',LMNX)
          CALL ERROR$STOP('DMFT_ULOCAL')
        END IF
        LM=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          LM=MAX(LM,SBARLI1(L+1,ISP)-1+(2*L+1))
        ENDDO
        IF(LMNX+LM.NE.LMNXT) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
          CALL ERROR$I4VAL('CALCULATED LMNXT',LMNX+LM)
          CALL ERROR$I4VAL('LMNX FROM POTPAR%TAILED',LMNXT)
          CALL ERROR$STOP('DMFT_ULOCAL')
        END IF
      END IF
!
!     ==========================================================================
!     == COLLECT STRUCTURE CONSTANTS                                          ==
!     ==========================================================================
      LMX=LMNXT-LMNX
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
        IF(SBAR(NN)%IAT2.NE.IAT) CYCLE
        IF(SUM(SBAR(NN)%IT(:)**2).NE.0) CYCLE
        IF(SBAR(NN)%N1.NE.LMX.OR.SBAR(NN)%N2.NE.LMX) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS FOR SBAR')
          CALL ERROR$I4VAL('LMX=LMNXT-LMNX',LMX)
          CALL ERROR$I4VAL('SBAR%N1',SBAR(NN)%N1)
          CALL ERROR$I4VAL('SBAR%N2',SBAR(NN)%N2)
          CALL ERROR$STOP('DMFT_ULOCAL')
        END IF
        NN0=NN
        EXIT
      ENDDO
      NN=NN0
!
!     ==========================================================================
!     == DOWNFOLD UTENSOR
!     ==   |KBAR_I>=|K_I>-\SUM_J |JBAR_J> SBAR(I,J)                           ==
!     ==========================================================================
!      U(:,:,::)=POTPAR(ISP)%TAILED%U(:LMNX,:LMNX,:LMNX,:LMNX)
      LMN=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        LM=SBARLI1(L+1,ISP)-1
        DO IM=1,2*L+1 
          LMN=LMN+1
          LM=LM+1
          DO J=1,LMX
            SVAR=SBAR(NN)%MAT(LM,J)   
            UT(:,:,:,LMN)=UT(:,:,:,LMN)-SVAR*UT(:,:,:,LMNX+J)
          ENDDO
        ENDDO
      ENDDO
!
      LMN=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        LM=SBARLI1(L+1,ISP)-1
        DO IM=1,2*L+1 
          LMN=LMN+1
          LM=LM+1
          DO J=1,LMX
            SVAR=SBAR(NN)%MAT(LM,J)   
            UT(:,:,LMN,:LMNX)=UT(:,:,LMN,:LMNX)-SVAR*UT(:,:,LMNX+J,:LMNX)
          ENDDO
        ENDDO
      ENDDO
!
      LMN=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        LM=SBARLI1(L+1,ISP)-1
        DO IM=1,2*L+1 
          LMN=LMN+1
          LM=LM+1
          DO J=1,LMX
            SVAR=SBAR(NN)%MAT(LM,J)   
            UT(:,LMN,:LMNX,:LMNX)=UT(:,LMN,:LMNX,:LMNX) &
      &                         -SVAR*UT(:,LMNX+J,:LMNX,:LMNX)
          ENDDO
        ENDDO
      ENDDO
!
      LMN=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        LM=SBARLI1(L+1,ISP)-1
        DO IM=1,2*L+1 
          LMN=LMN+1
          LM=LM+1
          DO J=1,LMX
            SVAR=SBAR(NN)%MAT(LM,J)   
            UT(LMN,:LMNX,:LMNX,:LMNX)=UT(LMN,:LMNX,:LMNX,:LMNX) &
                                     -SVAR*UT(LMNX+J,:LMNX,:LMNX,:LMNX)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == MAP ON SMALLER ARRAY                                                 ==
!     ==========================================================================
      U=UT(:LMNX,:LMNX,:LMNX,:LMNX)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$HFSOLVER_WITHSET(NLOC,NDIMD,NOMEGA,KBT,OMEGA &
     &                        ,GLOC,GLOCLAUR,UCHI,ETOT,SLOC,SLOCLAUR)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NLOC
      INTEGER(4),INTENT(IN)  :: NDIMD
      INTEGER(4),INTENT(IN)  :: NOMEGA
      REAL(8)   ,INTENT(IN)  :: KBT
      REAL(8)   ,INTENT(IN)  :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN)  :: GLOC(NLOC,NLOC,NDIMD,NOMEGA)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR(NLOC,NLOC,NDIMD,3)
      REAL(8)   ,INTENT(OUT) :: UCHI(NLOC,NLOC,NLOC,NLOC)
      REAL(8)   ,INTENT(OUT) :: ETOT
      COMPLEX(8),INTENT(OUT) :: SLOC(NLOC,NLOC,NDIMD,NOMEGA)
      COMPLEX(8),INTENT(OUT) :: SLOCLAUR(NLOC,NLOC,NDIMD,3)
      LOGICAL(4),PARAMETER   :: TPRINT=.TRUE.
      COMPLEX(8)             :: RHO(NLOC,NLOC,NDIMD)
      COMPLEX(8)             :: H(NLOC,NLOC,NDIMD)
      INTEGER(4)             :: NU,I,J,K,L
      REAL(8)                :: SVAR
!     **************************************************************************
                                     CALL TRACE$PUSH('DMFT$HFSOLVER_WITHSET')
!
!     ==========================================================================
!     == CALCULATE ONE-PARTICLE DENSITY MATRIX                                ==
!     ==========================================================================
      CALL DMFT_RHO_WITHSET(NLOC,NDIMD,NOMEGA,KBT,OMEGA,GLOC,GLOCLAUR,RHO)
!
      IF(TPRINT) THEN
        CALL SPINOR$PRINTMATRIX(6 &
     &    ,'DENSITY MATRIX (UPDOWN) IN DMFT$HFSOLVER_WITHSET','UPDOWN',1,NLOC &
     &                       ,NDIMD,NLOC,RHO)
      END IF
!
!     ==========================================================================
!     == CALCULATE HARTREE FOCK POTENTIAL                                     ==
!     ==========================================================================
      ETOT=0.D0
      H=(0.D0,0.D0)
      DO I=1,NLOC
        DO J=1,NLOC
          DO K=1,NLOC
            DO L=1,NLOC
              SVAR=-0.25D0*UCHI(I,J,K,L) 
              IF(SVAR.EQ.0.D0) CYCLE
              ETOT=ETOT+SVAR*REAL(SUM(RHO(K,J,:)*RHO(L,I,:)))
              H(K,J,:)=H(K,J,:)+SVAR*RHO(I,L,:)
              H(I,L,:)=H(I,L,:)+SVAR*RHO(K,J,:)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO NU=1,NOMEGA
        SLOC(:,:,:,NU)=H(:,:,:)
      ENDDO
      SLOCLAUR(:,:,:,1)=H
      SLOCLAUR(:,:,:,2)=(0.D0,0.D0)
      SLOCLAUR(:,:,:,3)=(0.D0,0.D0)
PRINT*,'EX IN DMFT$HFSOLVER_WITHSET:',ETOT
!
      IF(TPRINT) THEN
        CALL SPINOR$PRINTMATRIX(6 &
     &                      ,'H_XC (UPDOWN) IN DMFT$HFSOLVER_WITHSET','UPDOWN' &
     &                      ,1,NLOC,NDIMD,NLOC,H)
      END IF
!
                                     CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_CONSTRAINTS_WITHKSET(TYPE)
!     **************************************************************************
!     **  ADJUSTS H0(:,:,IKPT,IPSIN) SUCH THAT                                **
!     **                                                                      **
!     **  G(K,NU)=<PI(K)|PSI(K)> &                                            **
!     **         /[I*OMEGA_NU+MU-<PSI(K)|PI>(H0(K)-SIGMA(NU))<PI(K)|PSI(K)>] &**
!     **         *<PSI(K)|PI(K)>                                              **
!     **                                                                      **
!     **  PRODUCES THE SAME K-DEPENDENT DENSITY MATRIX AS GRHO                **
!     **                                                                      **
!     **  GRHO(K,NU)=<PI(K)|PSI(K)>/[I*OMEGA_NU+MU-HRHO(K)]<PI(K)|PSI(K)>     **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NSPIN,NDIMD,NAT,NOMEGA &
     &                      ,OMEGA,KBT,MU,AMIX,KSET,ATOMSET
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: TYPE  ! CAN BE 'HRHO','H0'
      REAL(8)   ,PARAMETER     :: TOL=1.D-6
      INTEGER(4),PARAMETER     :: NITER=1000
      LOGICAL(4),PARAMETER     :: TMIX=.FALSE.
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT3(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: G(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR1(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR3(NCHI,NCHI,NDIMD)
!      COMPLEX(8)               :: S(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: SLAUR1(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: SLAUR2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: SLAUR3(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: DEVRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),ALLOCATABLE   :: X4(:,:,:)
      COMPLEX(8)               :: DH0(NCHI,NCHI,NDIMD)
      LOGICAL(4)               :: CONVG
      REAL(8)                  :: MAXDEV
      COMPLEX(8)               :: CSVAR
      INTEGER(4)               :: IKPT,ITER,IDIMD,NU,IAT
      INTEGER(4)               :: I1,I2
      LOGICAL(4)               :: TH0
      INTEGER(4)               :: LX4 ! FIRST DIMENSION OF X4
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_CONSTRAINTS_WITHKSET')
      TH0=(TYPE.EQ.'H0') 
      IF(.NOT.(TH0.OR.TYPE.EQ.'HRHO')) THEN
        CALL ERROR$MSG('ILLEGAL VALUE OF TYPE. (MAY BE "H0" OR "HRHO")')
        CALL ERROR$CHVAL('TYPE',TYPE)
        CALL ERROR$STOP('DMFT_CONSTRAINTS_WITHKSET')     
      END IF
      IF(.NOT.TH0)  THEN
        DO IKPT=1,NKPTL
          KSET(IKPT)%H0=(0.D0,0.D0)
        ENDDO
      END IF
!
!     == ALLOCATE X4 ===========================================================
      IF(NDIMD.LE.2) THEN
        LX4=NCHI**2
      ELSE IF(NDIMD.EQ.4) THEN
        LX4=(2*NCHI)**2
      ENDIF
      ALLOCATE(X4(LX4,LX4,NSPIN))
!       
!     ==========================================================================
!     ==  ADJUST H0 TO OBEY DENSITY MATRIX CONSTRAINT                         ==
!     ==========================================================================
      DO ITER=1,NITER
!       
!       ========================================================================
!       ==  DEVIATION FROM TARGET DENSITY MATRIX                              ==
!       ========================================================================
        MAXDEV=0.D0
        DO IKPT=1,NKPTL
!         == LAURENT EXPANSION FOR THE GREENS FUNCTION =========================
          SLAUR1=(0.D0,0.D0)
          SLAUR2=(0.D0,0.D0)
          SLAUR3=(0.D0,0.D0)
          IF(TH0) THEN
            DO IAT=1,NAT
              I1=ATOMSET(IAT)%ICHI1
              I2=ATOMSET(IAT)%ICHI2
              SLAUR1(I1:I2,I1:I2,:)=ATOMSET(IAT)%SLOCLAUR(:,:,:,1)
              SLAUR2(I1:I2,I1:I2,:)=ATOMSET(IAT)%SLOCLAUR(:,:,:,2)
              SLAUR3(I1:I2,I1:I2,:)=ATOMSET(IAT)%SLOCLAUR(:,:,:,3)
            ENDDO
          END IF
          GLAUR1=KSET(IKPT)%SINV
!
          MAT=-MU*KSET(IKPT)%SMAT
          MAT=MAT+KSET(IKPT)%H0
          IF(TH0) MAT=MAT+SLAUR1(:,:,:)
!         == GLAUR2=SINV*MAT*SINV ==============================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT,KSET(IKPT)%SINV,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,MAT2,GLAUR2)
!
          MAT=-MU*KSET(IKPT)%SMAT
          MAT=MAT+KSET(IKPT)%H0
          IF(TH0) MAT=MAT+SLAUR1(:,:,:)
!         == MAT3=MAT*(SINV*MAT) ===============================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT,KSET(IKPT)%SINV,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT,MAT2,MAT3)
          IF(TH0) MAT3=MAT3+SLAUR2(:,:,:)
!         == GLAUR3=SINV*(MAT3*SINV) ===========================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT3,KSET(IKPT)%SINV,MAT)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,MAT,GLAUR3)
!
          X4=(0.D0,0.D0)        
          DEVRHO=(0.D0,0.D0)
          DO NU=1,NOMEGA
        
!           == CONSTRUCT LATTICE GREENS FUNCTION =============================
            MAT=(CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT
            MAT=MAT-KSET(IKPT)%H0
            IF(TH0) THEN
              DO IAT=1,NAT
                I1=ATOMSET(IAT)%ICHI1
                I2=ATOMSET(IAT)%ICHI2
                MAT(I1:I2,I1:I2,:)=MAT(I1:I2,I1:I2,:) &
     &                            -ATOMSET(IAT)%SLOC(:,:,:,NU)
              ENDDO
            END IF
            CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
!           == SUBTRACT LAURENT EXPANSION TO IMPROVE OMEGA-CONVERGENCE =========
            CSVAR=1.D0/(CI*OMEGA(NU))
            DEVRHO=DEVRHO+KBT*(G-CSVAR*(GLAUR1+CSVAR*(GLAUR2+CSVAR*GLAUR3)))
!
!           == ACCUMULATE SECOND ORDER TERM ====================================
            CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,G) ! CONVERT TO UPDOWN
            CALL DMFT_ADDX4(NDIMD,NCHI,G,LX4,NSPIN,X4)
          ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
          X4=X4*KBT
!
!         == INCLUDE NEGATIVE FREQUENCIES ======================================
!         == ALREADY INCLUDED FOR X4
          DO IDIMD=1,NDIMD
            DEVRHO(:,:,IDIMD)=DEVRHO(:,:,IDIMD) &
     &                       +CONJG(TRANSPOSE(DEVRHO(:,:,IDIMD)))
          ENDDO
!         == ADD TAILS (GLAUR3 DOES NOT CONTRIBUTE) ============================
          DEVRHO=DEVRHO+0.5D0*GLAUR1-0.25D0/KBT*GLAUR2
!
!         == SUBTRACT TARGET DENSITY =========================================
          DEVRHO=DEVRHO-KSET(IKPT)%RHO
!
!!$          MAXDEV=MAX(MAXDEV,MAXVAL(ABS(DEVRHO)))

CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,DEVRHO) ! CONVERT TO UPDOWN
MAXDEV=MAX(MAXDEV,MAXVAL(ABS(DEVRHO)))
DO IDIMD=1,NDIMD
    PRINT*,'ITER=',ITER,' IKPT=',IKPT,' IDIMD=',IDIMD &
  &       ,' MAXVAL(X4)  ',MAXVAL(ABS(X4(:,:,IDIMD))) &
  &       ,' MAX(DEVRHO) ',MAXVAL(ABS(DEVRHO(:,:,IDIMD)))
ENDDO
CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,DEVRHO) ! CONVERT TO TXYZ
!
          IF(TMIX) THEN
            DH0=AMIX*DEVRHO
          ELSE
            CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,DEVRHO) ! CONVERT TO UPDOWN
!           == DEVRHO+X4*DH0=0 -> DH0
            CALL DMFT_SOLVEX4(NDIMD,NCHI,DEVRHO,LX4,NSPIN,X4,DH0)
            CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,DH0)    ! CONVERT TO TXYZ
          END IF
!
!         =================================================================== 
!         == ADD DH
!         ====================================================================
          KSET(IKPT)%H0=KSET(IKPT)%H0+DH0
        ENDDO ! END OF LOOP OVER K-POINTS
!
!       ========================================================================
!       == CHECK CONVERGENCE
!       ========================================================================
        PRINT*,'MAXDEV',ITER,MAXDEV
        CONVG=MAXDEV.LT.TOL
        IF(CONVG) EXIT
      ENDDO   !END OF LOOP OVER ITERATIONS
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('LOOP NOT CONVERGED')
        CALL ERROR$STOP('DMFT_CONSTRAINTS')
      END IF
!
!     ==========================================================================
!     == MAP ONTO HRHO                                                        ==
!     ==========================================================================
      IF(.NOT.TH0) THEN
        DO IKPT=1,NKPTL
          KSET(IKPT)%HRHO=KSET(IKPT)%H0
        ENDDO
      END IF
PRINT*,'TH0 ',TH0

                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVEX4(NDIMD,NCHI,DRHO,LX4,NSPIN,X4,DH)
!     **************************************************************************
!     ** DRHO+X4*DH=0
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NDIMD
      INTEGER(4),INTENT(IN)    :: NCHI
      INTEGER(4),INTENT(IN)    :: LX4
      INTEGER(4),INTENT(IN)    :: NSPIN
      COMPLEX(8),INTENT(IN)    :: DRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN)    :: X4(LX4,LX4,NSPIN)
      COMPLEX(8),INTENT(OUT)   :: DH(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: BVEC(LX4)
      COMPLEX(8)               :: XVEC(LX4)
      INTEGER(4)               :: IA,IB
      INTEGER(4)               :: IDA,IDB
      INTEGER(4)               :: ICHIA,ICHIB
      INTEGER(4)               :: IAB
      INTEGER(4)               :: IDIMAB
      INTEGER(4)               :: ISPIN 
!     **************************************************************************
!
!     ==========================================================================
!     == COLLINEAR CASE: SPIN INDICES ARE DECOUPLED                           ==
!     ==========================================================================
      IF(NDIMD.LE.2) THEN
        IF(LX4.NE.NCHI**2.OR.NSPIN.NE.NDIMD) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$I4VAL('NSPIN',NSPIN)
          CALL ERROR$I4VAL('NCHI',NCHI)
          CALL ERROR$I4VAL('LX4',LX4)
          CALL ERROR$STOP('DMFT_ADDX4')
        END IF
        DO ISPIN=1,NSPIN
          DO IA=1,NCHI
            DO IB=1,NCHI
              IAB=IA+NCHI*(IB-1)
              BVEC(IAB)=DRHO(IA,IB,ISPIN)
            ENDDO
          ENDDO              
!         == A*X=B -> X (THEREFORE MINUS SIGN NEEDED)===========================
          CALL LIB$MATRIXSOLVEC8(LX4,LX4,1,X4(:,:,ISPIN),XVEC,BVEC)
          DO IA=1,NCHI
            DO IB=1,NCHI
              IAB=IA+NCHI*(IB-1)
              DH(IA,IB,ISPIN)=-XVEC(IAB)
            ENDDO
          ENDDO              
        ENDDO
!
!     ==========================================================================
!     == NON-COLLINEAR CASE:                                                  ==
!     ==========================================================================
      ELSE IF(NDIMD.EQ.4) THEN
        IF(LX4.NE.(2*NCHI)**2.OR.NSPIN.NE.1) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$I4VAL('NSPIN',NSPIN)
          CALL ERROR$I4VAL('NCHI',NCHI)
          CALL ERROR$I4VAL('LX4',LX4)
          CALL ERROR$STOP('DMFT_ADDX4')
        END IF
        DO IDA=1,2
          DO ICHIA=1,NCHI
            IA=ICHIA+NCHI*(IDA-1)
            DO IDB=1,2
              DO ICHIB=1,NCHI
                IB=ICHIB+NCHI*(IDB-1)
                IAB=IA+2*NCHI*(IB-1)
                IDIMAB=IDA+2*(IDB-1)
                BVEC(IAB)=DRHO(ICHIA,ICHIB,IDIMAB)
              ENDDO
            ENDDO              
          ENDDO              
        ENDDO              
!         == A*X=B -> X (THEREFORE MINUS SIGN NEEDED)===========================
        CALL LIB$MATRIXSOLVEC8(LX4,LX4,1,X4,XVEC,BVEC)
        DO IDA=1,2
          DO ICHIA=1,NCHI
            IA=ICHIA+NCHI*(IDA-1)
            DO IDB=1,2
              DO ICHIB=1,NCHI
                IB=ICHIB+NCHI*(IDB-1)
                IAB=IA+2*NCHI*(IB-1)
                IDIMAB=IDA+2*(IDB-1)
                DH(ICHIA,ICHIB,IDIMAB)=-XVEC(IAB)
              ENDDO
            ENDDO              
          ENDDO              
        ENDDO              
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_ADDX4(NDIMD,NCHI,G,LX4,NSPIN,X4)
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NDIMD
      INTEGER(4),INTENT(IN)    :: NCHI
      INTEGER(4),INTENT(IN)    :: LX4
      INTEGER(4),INTENT(IN)    :: NSPIN
      COMPLEX(8),INTENT(IN)    :: G(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(INOUT) :: X4(LX4,LX4,NSPIN)
      INTEGER(4)               :: IA,IB,IC,ID
      INTEGER(4)               :: IDA,IDB,IDC,IDD
      INTEGER(4)               :: ICHIA,ICHIB,ICHIC,ICHID
      INTEGER(4)               :: IAB,ICD
      INTEGER(4)               :: IDIMAC,IDIMDB,IDIMCA,IDIMBD
      INTEGER(4)               :: ISPIN 
!     **************************************************************************
!
!     ==========================================================================
!     == COLLINEAR CASE: SPIN INDICES ARE DECOUPLED                           ==
!     ==========================================================================
      IF(NDIMD.LE.2) THEN
        IF(LX4.NE.NCHI**2.OR.NSPIN.NE.NDIMD) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$I4VAL('NSPIN',NSPIN)
          CALL ERROR$I4VAL('NCHI',NCHI)
          CALL ERROR$I4VAL('LX4',LX4)
          CALL ERROR$STOP('DMFT_ADDX4')
        END IF
        DO ISPIN=1,NDIMD
          IAB=0
          DO IA=1,NCHI
            DO IB=1,NCHI
              IAB=IAB+1
              ICD=0
              DO IC=1,NCHI
                DO ID=1,NCHI
                  ICD=ICD+1
!                 == THE SECOND TERM IS FOR -I*OMEGA_\NU =======================
!                 == G(-I*OMEGA)=G^\DAGGER(I*OMEGA) ============================
                  X4(IAB,ICD,ISPIN)=X4(IAB,ICD,ISPIN) &
     &                             +G(IA,IC,ISPIN)*G(ID,IB,ISPIN) &
     &                             +CONJG(G(IC,IA,ISPIN)*G(IB,ID,ISPIN))
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
!     ==========================================================================
!     == NON-COLLINEAR CASE:                                                  ==
!     ==========================================================================
      ELSE
        IF(LX4.NE.(2*NCHI)**2.OR.NSPIN.NE.1) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$I4VAL('NSPIN',NSPIN)
          CALL ERROR$I4VAL('NCHI',NCHI)
          CALL ERROR$I4VAL('LX4',LX4)
          CALL ERROR$STOP('DMFT_ADDX4')
        END IF
        IAB=0
        DO IDA=1,2
        DO ICHIA=1,NCHI
          DO IDB=1,2
          DO ICHIB=1,NCHI
            IAB=IAB+1
            ICD=0
            DO IDC=1,2
            DO ICHIC=1,NCHI
              DO IDD=1,2
              DO ICHID=1,NCHI
                ICD=ICD+1
!
                IDIMAC=IDA+2*(IDC-1)
                IDIMCA=IDC+2*(IDA-1)
                IDIMDB=IDD+2*(IDB-1)
                IDIMBD=IDB+2*(IDD-1)
                X4(IAB,ICD,1)=X4(IAB,ICD,1) &
     &                       +G(ICHIA,ICHIC,IDIMAC)*G(ICHID,ICHIB,IDIMDB) &
     &                 +CONJG(G(ICHIC,ICHIA,IDIMCA)*G(ICHIB,ICHID,IDIMBD))
              ENDDO
              ENDDO
            ENDDO
            ENDDO
          ENDDO
          ENDDO
        ENDDO
        ENDDO
      END IF      
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_RHO_WITHSET(NCHI,NDIMD,NOMEGA,KBT,OMEGA,G,GLAUR,RHO)
!     **************************************************************************
!     **  CALCULATES THE DENSITY MATRIX FROM A GREENS FUNCTION                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: KBT
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN) :: G(NCHI,NCHI,NDIMD,NOMEGA)
      COMPLEX(8),INTENT(IN) :: GLAUR(NCHI,NCHI,NDIMD,3)
      COMPLEX(8),INTENT(OUT):: RHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)  ! SQRT(-1)
      INTEGER(4)            :: IDIMD,NU
      COMPLEX(8)            :: CSVAR
!     **************************************************************************
      RHO=(0.D0,0.D0)
      DO NU=1,NOMEGA
        CSVAR=1.D0/(CI*OMEGA(NU))
        RHO=RHO+G(:,:,:,NU)-GLAUR(:,:,:,1)*CSVAR &
     &                     -GLAUR(:,:,:,2)*CSVAR**2 &
     &                     -GLAUR(:,:,:,3)*CSVAR**3
      ENDDO
      DO IDIMD=1,NDIMD
        RHO(:,:,IDIMD)=RHO(:,:,IDIMD)+CONJG(TRANSPOSE(RHO(:,:,IDIMD)))
      ENDDO
      RHO=RHO*KBT
      RHO=RHO+0.5D0*GLAUR(:,:,:,1)-0.25D0/KBT*GLAUR(:,:,:,2)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_GLOC_WITHATOMSET()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NDIMD,NAT,NOMEGA,OMEGA,MU,KSET,ATOMSET
      IMPLICIT NONE
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER     :: TPRINT=.FALSE.
      LOGICAL(4),PARAMETER     :: TTEST=.FALSE.
      COMPLEX(8)               :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT3(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR1(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR3(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: SLAUR1(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: SLAUR2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: SLAUR3(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: G(NCHI,NCHI,NDIMD)
      REAL(8)                  :: WKPTL
      INTEGER(4)               :: IKPT,NU,IAT,I1,I2,I,IDIMD
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_GLOC')
      DO IAT=1,NAT
        ATOMSET(IAT)%GLOC=(0.D0,0.D0)
        ATOMSET(IAT)%GLOCLAUR=(0.D0,0.D0)
      ENDDO
!
      DO IKPT=1,NKPTL
        WKPTL=KSET(IKPT)%WKPTL
!       == LAURENT EXPANSION FOR THE GREENS FUNCTION =========================
        SLAUR1=(0.D0,0.D0)
        SLAUR2=(0.D0,0.D0)
        SLAUR3=(0.D0,0.D0)
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          SLAUR1(I1:I2,I1:I2,:)=ATOMSET(IAT)%SLOCLAUR(:,:,:,1)
          SLAUR2(I1:I2,I1:I2,:)=ATOMSET(IAT)%SLOCLAUR(:,:,:,2)
          SLAUR3(I1:I2,I1:I2,:)=ATOMSET(IAT)%SLOCLAUR(:,:,:,3)
        ENDDO
        GLAUR1=KSET(IKPT)%SINV
!
        MAT=-MU*KSET(IKPT)%SMAT
        MAT=MAT+KSET(IKPT)%H0
        MAT=MAT+SLAUR1
!       == GLAUR2=SINV*MAT*SINV ==============================================
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT,KSET(IKPT)%SINV,MAT2)
        CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,MAT2,GLAUR2)
!
        MAT=-MU*KSET(IKPT)%SMAT
        MAT=MAT+KSET(IKPT)%H0
        MAT=MAT+SLAUR1
!       == MAT3=MAT*(SINV*MAT) ===============================================
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT,KSET(IKPT)%SINV,MAT2)
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT,MAT2,MAT3)
        MAT3=MAT3+SLAUR2
!       == GLAUR3=SINV*(MAT3*SINV) ===========================================
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT3,KSET(IKPT)%SINV,MAT)
        CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,MAT,GLAUR3)
        IF(KSET(IKPT)%TADDMINUSK) THEN
          DO IDIMD=1,NDIMD
              GLAUR1(:,:,IDIMD)=0.5D0*(GLAUR1(:,:,IDIMD) &
      &                               +CONJG(TRANSPOSE(GLAUR1(:,:,IDIMD))))
              GLAUR2(:,:,IDIMD)=0.5D0*(GLAUR2(:,:,IDIMD) &
      &                               +CONJG(TRANSPOSE(GLAUR2(:,:,IDIMD))))
              GLAUR3(:,:,IDIMD)=0.5D0*(GLAUR3(:,:,IDIMD) &
      &                               +CONJG(TRANSPOSE(GLAUR3(:,:,IDIMD))))
          ENDDO
        END IF
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          ATOMSET(IAT)%GLOCLAUR(:,:,:,1)=ATOMSET(IAT)%GLOCLAUR(:,:,:,1) &
     &                                  +WKPTL*GLAUR1(I1:I2,I1:I2,:)
          ATOMSET(IAT)%GLOCLAUR(:,:,:,2)=ATOMSET(IAT)%GLOCLAUR(:,:,:,2) &
     &                                  +WKPTL*GLAUR2(I1:I2,I1:I2,:)
          ATOMSET(IAT)%GLOCLAUR(:,:,:,3)=ATOMSET(IAT)%GLOCLAUR(:,:,:,3) &
     &                                  +WKPTL*GLAUR3(I1:I2,I1:I2,:)
        ENDDO
!
        DO NU=1,NOMEGA
!         == CONSTRUCT LATTICE GREENS FUNCTION =================================
          MAT=(CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT
          MAT=MAT-KSET(IKPT)%H0
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            MAT(I1:I2,I1:I2,:)=MAT(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
          ENDDO
          CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
!         == ACCOUNT FOR K-POINT (-K) ==========================================
          IF(KSET(IKPT)%TADDMINUSK) THEN
            DO IDIMD=1,NDIMD
              MAT(:,:,IDIMD)=2.D0*CI*OMEGA(NU) &
      &                          *TRANSPOSE(CONJG(KSET(IKPT)%SMAT(:,:,IDIMD))) &
      &                          -TRANSPOSE(CONJG(MAT(:,:,IDIMD)))
            ENDDO
            CALL SPINOR$INVERT(NDIMD,NCHI,MAT,MAT2)
            G=0.5D0*(G+MAT2)
          END IF
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            ATOMSET(IAT)%GLOC(:,:,:,NU)=ATOMSET(IAT)%GLOC(:,:,:,NU) &
     &                                 +WKPTL*G(I1:I2,I1:I2,:)
          ENDDO
!
        ENDDO ! END OF LOOP OVER MATSUBARA FREQUENCIES
      ENDDO   ! END OF LOOP OVER KPOINTS
!
!     ==========================================================================
!     == TAKE K-POINTS INTO ACCOUNT THAT HAVE BEEN LEFT OUT DUE TO            ==
!     == TIME INVERSION SYMMETRY                                              ==
!     ==========================================================================
!     == ALL LAURENT EXPANSION COEFFICIENTS ARE HERMITIAN ======================
      DO IAT=1,NAT
        I1=ATOMSET(IAT)%ICHI1
        I2=ATOMSET(IAT)%ICHI2
        DO I=1,3
          DO IDIMD=1,NDIMD
            ATOMSET(IAT)%GLOCLAUR(:,:,IDIMD,I)=0.5D0 &
     &                   *(ATOMSET(IAT)%GLOCLAUR(:,:,IDIMD,I) &
     &                    +CONJG(TRANSPOSE(ATOMSET(IAT)%GLOCLAUR(:,:,IDIMD,I))))
          ENDDO
        ENDDO
      ENDDO

                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$ADDTOHPSI_WITHKSET()
!     **************************************************************************
!     ** 
!     ** 
!     ** 
!     ** 
!     ** 
!     **************************************************************************
      USE DMFT_MODULE, ONLY  : TON,NKPTL,NSPIN,NDIM,NDIMD,NB,NCHI,NAT &
     &                        ,IPROOFCHI,ATOMSET,KSET
      USE WAVES_MODULE, ONLY : GSET,WAVES_SELECTWV,THIS,MAP
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      COMPLEX(8),ALLOCATABLE :: DHPIPSI(:,:,:) !(NDIM,NCHI,NB)
      COMPLEX(8),ALLOCATABLE :: MAT(:,:,:)     !(NCHI,NCHI,NDIMD)
      LOGICAL(4)             :: TCHK,TRESET
      INTEGER(4)             :: NBH
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IPRO,IDIM1,IDIM2,IDIMD,IAT
      INTEGER(4)             :: I1,I2
!     **************************************************************************
      IF(.NOT.TON) RETURN
PRINT*,'ENTERING DMFT$ADDTOHPSI_WITHKSET'
                              CALL TRACE$PUSH('DMFT$ADDTOHPSI_WITHKSET')
!
!     ==========================================================================
!     == CHECK IF HTBC ALREADY CONTAINS INFORMATION                           ==
!     ==========================================================================
      CALL LMTO$GETL4('ON',TCHK)
      IF(TCHK) CALL LMTO$GETL4('THTBC',TCHK)
      TRESET=.NOT.TCHK
!
!     ==========================================================================
!     ==  COLLECT CONSTANTS
!     ==========================================================================
      NPRO=MAP%NPRO
      ALLOCATE(DHPIPSI(NDIM,NCHI,NB))
      ALLOCATE(MAT(NCHI,NCHI,NDIMD))
!
!     ==========================================================================
!     ==  CONVERT ON-SITE TERMS FROM (TXYZ) TO UPDOWN                         ==
!     ==========================================================================
      DO IAT=1,NAT
        LMNX=ATOMSET(IAT)%DENMAT%LMNX
        CALL SPINOR$CONVERT('BACK',LMNX,NDIMD,ATOMSET(IAT)%DENMAT%H)
        ATOMSET(IAT)%DENMAT%H=2.D0*ATOMSET(IAT)%DENMAT%H
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          IF(THIS%NB.NE.NB) THEN
            CALL ERROR$MSG('INCONSISTENT NUMBER OF STATES IN WAVES AND DYNOCC')
            CALL ERROR$I4VAL('NB IN DYNOCC',NB)
            CALL ERROR$I4VAL('NB IN WAVES ',THIS%NB)
            CALL ERROR$STOP('DMFT$ADDTOHPSI')
          END IF
          NBH=THIS%NBH
!
!         ======================================================================
!         == ALLOCATE THIS%HTBC IF NECESSARY                                  ==
!         ======================================================================
          IF(TRESET) THEN
            IF(.NOT.ASSOCIATED(THIS%HTBC))ALLOCATE(THIS%HTBC(NDIM,NBH,NPRO))
            THIS%HTBC=(0.D0,0.D0)
          END IF
!
!         ======================================================================
!         == DF/DRHO=H0-HRHO; DHDPIPSI= (H0-HRHO)<PI|PSI>                     ==
!         ======================================================================
          MAT=KSET(IKPT)%H0(:,:,:)-KSET(IKPT)%HRHO(:,:,:)
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,MAT)
          MAT=2.D0*MAT   
          IDIMD=0
          DO IDIM2=1,NDIM
            DO IDIM1=1,NDIM
              IDIMD=IDIM1+NDIM*(IDIM2-1)+ISPIN-1
              DHPIPSI(IDIM1,:,:)=MATMUL(MAT(:,:,IDIMD) &
       &                                     ,KSET(IKPT)%PIPSI(IDIM2,:,:,ISPIN))
            ENDDO
          ENDDO
!
!         ======================================================================
!         ==  CONVERT TO SUPER WAVE FUNCTIONS                                 ==
!         ==  PIPSI(ICHI,2*IBH-1,IKPT,ISPIN)= REAL(THIS%TBC(1,IBH,IPRO))      ==
!         ==  PIPSI(ICHI,2*IBH  ,IKPT,ISPIN)=AIMAG(THIS%TBC(1,IBH,IPRO))      ==
!         ======================================================================
          IF(NBH.NE.NB) THEN
            DO IBH=1,NBH
              DHPIPSI(:,:,IBH)=REAL(DHPIPSI(:,:,2*IBH-1)) &
        &                 -CI*AIMAG(DHPIPSI(:,:,2*IBH)) 
            ENDDO
          END IF
!
!         ======================================================================
!         ==  DETERMINE CONTRIBUTION TO PROJECTOR PART OF H|PSI>              ==
!         ==  EXPAND TO ALL PROJECTOR FUNCTIONS                               ==
!         ======================================================================
          DO ICHI=1,NCHI
            IPRO=IPROOFCHI(ICHI)
!           == PIPSI(ICHI,IBH,IKPT,ISPIN)  =THIS%TBC(1,IBH,IPRO) ===============
            THIS%HTBC(:,:,IPRO)=THIS%HTBC(:,:,IPRO)+DHPIPSI(:,ICHI,:NBH)
          ENDDO
!
!         ======================================================================
!         ==  NOW THE SITE-LOCAL TERM FROM THE DOUBLE COUNTING                ==
!         ======================================================================
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            DO IDIM2=1,NDIMD
              DO IDIM1=1,NDIMD
                IDIMD=IDIM1+NDIM*(IDIM2-1)+ISPIN-1
!               == HTBC=HAM*TBC ================================================
!               == HTBC(IB,I)=HTBC(IB,I)+HAM(I,J)*TBC(IB,J) ====================
                THIS%HTBC(IDIM1,:,I1:I2)=THIS%HTBC(IDIM1,:,I1:I2) &
         &                +MATMUL(THIS%TBC(IDIM1,:,I1:I2) &
         &                       ,TRANSPOSE(ATOMSET(IAT)%DENMAT%H(:,:,IDIMD)))
              ENDDO
            ENDDO
          ENDDO
        ENDDO ! END OF LOOP OVER NSPIN
      ENDDO ! END OF LOOP OVER K-POINTS
PRINT*,'.... DMFT$ADDTOHPSI DONE'
                              CALL TRACE$POP()
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE DMFT_PRECONDITION(M,C,NOMEGA,OMEGA,MOMEGA)
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: M
      REAL(8)   ,INTENT(IN) :: C
      INTEGER(4),INTENT(IN) :: NOMEGA
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(OUT):: MOMEGA(NOMEGA)
      INTEGER(4)            :: NU
      REAL(8)               :: PI,PIHALF
      REAL(8)               :: FACTOR1,FACTOR2
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      PIHALF=0.5D0*PI
      FACTOR1=M/PIHALF**2
      FACTOR2=PIHALF*SQRT(C/M)
      DO NU=1,NOMEGA
        MOMEGA(NU)=FACTOR1*ATAN(FACTOR2/OMEGA(NU))**2
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TESTBESSEL()
      INTEGER(4)           :: L=2
      INTEGER(4),PARAMETER :: NKAPPA=5
      INTEGER(4)           :: I
      REAL(8)              :: RAD
      REAL(8)              :: K2(NKAPPA)=(/-0.2D0,-0.1D0,0.D0,1.D0,2.D0/)
      REAL(8)              :: JVAL(NKAPPA),JDER(NKAPPA)
      REAL(8)              :: KVAL(NKAPPA),KDER(NKAPPA)
      INTEGER(4)           :: NFIL1=12,NFIL2=14
      OPEN(NFIL1,FILE='JVAL.DAT')
      OPEN(NFIL2,FILE='KVAL.DAT')
      DO I=1,1000
        RAD=1.D-2*REAL(I)
        DO IKAPPA=1,NKAPPA
          CALL LMTO$SOLIDBESSELRAD(L,RAD,K2(IKAPPA),JVAL(IKAPPA),JDER(IKAPPA))
          CALL LMTO$SOLIDHANKELRAD(L,RAD,K2(IKAPPA),KVAL(IKAPPA),KDER(IKAPPA))
        ENDDO
        WRITE(NFIL1,*)RAD,JVAL,JDER
        WRITE(NFIL2,*)RAD,KVAL,KDER
      ENDDO
      CLOSE(NFIL1)
      CLOSE(NFIL2)
PRINT*,'.... BESSELTEST DONE'
STOP 'FORCED AFTER BESSELTEST'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINOR$TEST()
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NCHI=3
      INTEGER(4),PARAMETER :: NDIMD=4
      COMPLEX(8)           :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8)           :: B(NCHI,NCHI,NDIMD)
      COMPLEX(8)           :: C(NCHI,NCHI,NDIMD)
      COMPLEX(8)           :: D(NCHI,NCHI,NDIMD)
      COMPLEX(8)           :: MAT1(2*NCHI,2*NCHI)
      COMPLEX(8)           :: MAT2(2*NCHI,2*NCHI)
      INTEGER(4),PARAMETER :: NFIL=6
      REAL(8)              :: RAN1,RAN2
      INTEGER(4)           :: IDIMD,I,J
      REAL(8)  ,PARAMETER  :: TOL=1.D-7
      REAL(8)              :: DEV
!     **************************************************************************
                                             CALL TRACE$PUSH('SPINOR$TEST')
PRINT*,' STARTING SPINOR TEST'
      DO IDIMD=1,NDIMD
        DO I=1,NCHI
          DO J=1,NCHI
            CALL RANDOM_NUMBER(RAN1)
            CALL RANDOM_NUMBER(RAN2)
            A(I,J,IDIMD)=CMPLX(RAN1,RAN2)
          ENDDO
        ENDDO
      ENDDO
      DO IDIMD=1,NDIMD
        DO I=1,NCHI
          DO J=1,NCHI
            CALL RANDOM_NUMBER(RAN1)
            CALL RANDOM_NUMBER(RAN2)
            B(I,J,IDIMD)=CMPLX(RAN1,RAN2)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  TEST FORWARD AND BACK TRANSFORMATION                                ==
!     ==========================================================================
      C=A
      CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,A)
      CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,A)
      DEV=MAXVAL(ABS(C-A))
      IF(DEV.GT.TOL) THEN 
        CALL SPINOR$PRINTMATRIX(NFIL,'CONVERSION TEST (ZERO)','DIRECT',1,NCHI &
     &                       ,NDIMD,NCHI,C-A)
        CALL ERROR$MSG('CONVERSION TEST FAILED')
        CALL ERROR$STOP('SPINOR$TEST')
      END IF
!
!     ==========================================================================
!     ==  TEST MULTIPLICATION                                                 ==
!     ==========================================================================
      C=A
      D=B
      CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,C)
      CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,D)
      MAT1=(0.D0,0.D0)
      MAT2=(0.D0,0.D0)
      IF(NDIMD.EQ.1) THEN
        MAT1(:NCHI  ,:NCHI)  =C(:,:,1)
        MAT1(NCHI+1:,NCHI+1:)=C(:,:,1)
        MAT2(:NCHI  ,:NCHI)  =D(:,:,1)
        MAT2(NCHI+1:,NCHI+1:)=D(:,:,1)
      ELSE IF(NDIMD.EQ.2) THEN
        MAT1(:NCHI  ,:NCHI)  =C(:,:,1)
        MAT1(NCHI+1:,NCHI+1:)=C(:,:,2)
        MAT2(:NCHI  ,:NCHI)  =D(:,:,1)
        MAT2(NCHI+1:,NCHI+1:)=D(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        MAT1(:NCHI  ,:NCHI)  =C(:,:,1)
        MAT1(NCHI+1:,:NCHI)  =C(:,:,2)
        MAT1(:NCHI  ,NCHI+1:)=C(:,:,3)
        MAT1(NCHI+1:,NCHI+1:)=C(:,:,4)
        MAT2(:NCHI  ,:NCHI)  =D(:,:,1)
        MAT2(NCHI+1:,:NCHI)  =D(:,:,2)
        MAT2(:NCHI  ,NCHI+1:)=D(:,:,3)
        MAT2(NCHI+1:,NCHI+1:)=D(:,:,4)
      END IF
      MAT1=MATMUL(MAT1,MAT2)
      IF(NDIMD.EQ.1) THEN
        C(:,:,1)=MAT1(:NCHI  ,:NCHI) 
      ELSE IF(NDIMD.EQ.2) THEN
        C(:,:,1)=MAT1(:NCHI  ,:NCHI)
        C(:,:,2)=MAT1(NCHI+1:,NCHI+1:)
      ELSE IF(NDIMD.EQ.4) THEN
        C(:,:,1)=MAT1(:NCHI  ,:NCHI)
        C(:,:,2)=MAT1(NCHI+1:,:NCHI)
        C(:,:,3)=MAT1(:NCHI  ,NCHI+1:)
        C(:,:,4)=MAT1(NCHI+1:,NCHI+1:)
      END IF
      CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,C)
!
      CALL SPINOR$MATMUL(NDIMD,NCHI,A,B,D)
      DEV=MAXVAL(ABS(D-C))
      IF(DEV.GT.TOL) THEN 
        CALL SPINOR$PRINTMATRIX(NFIL,'MATMUL TEST (ZERO)','DIRECT',1,NCHI &
     &                       ,NDIMD,NCHI,C-D)
        CALL ERROR$MSG('MULTIPLICATION  TEST FAILED')
        CALL ERROR$STOP('SPINOR$TEST')
      END IF
!
!     ==========================================================================
!     ==  TEST INVERSION AND MULTIPLICATION                                   ==
!     ==========================================================================
PRINT*,' BEFORE SPINOR$INVERT'
      CALL SPINOR$INVERT(NDIMD,NCHI,A,B)
!
PRINT*,' BEFORE SPINOR$MATMUL'
      CALL SPINOR$MATMUL(NDIMD,NCHI,A,B,C)
!
PRINT*,' BEFORE SPINOR$PRINTL'
!      CALL SPINOR$PRINTMATRIX(NFIL,'IDENTITY TEST','TXYZ',1,NCHI,NDIMD,NCHI,C)
      CALL SPINOR$PRINTMATRIX(NFIL,'IDENTITY TEST','UPDOWN',1,NCHI,NDIMD,NCHI,C)
                                             CALL TRACE$POP()
      STOP
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINOR$MATMUL(NDIMD,NCHI,A,B,C)
!     **************************************************************************
!     **  C=A*B                                                               **
!     **                                                                      **
!     **  DERIVATION IN METHODS: SECTION "SECOND QUANTIZATION WITH            **
!     **  NON-ORTHONORMAL BVAIS SETS, SUBSECTION "INVERSION OF A MATRIX IN    **
!     **  SPINOR REPRESENTATION                                               **
!     **************************************PETER BLOECHL GOSLAR 2013***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      COMPLEX(8),INTENT(IN) :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: B(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(OUT) :: C(NCHI,NCHI,NDIMD)
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      INTEGER(4)             :: IDIMD,I,IP,IPP
!     **************************************************************************
      C(:,:,:)=0.D0
      C(:,:,1)=MATMUL(A(:,:,1),B(:,:,1))
      DO IDIMD=2,NDIMD
        C(:,:,1)=C(:,:,1)+MATMUL(A(:,:,IDIMD),B(:,:,IDIMD))
        C(:,:,IDIMD)=C(:,:,IDIMD)+MATMUL(A(:,:,1),B(:,:,IDIMD)) &
     &                           +MATMUL(A(:,:,IDIMD),B(:,:,1)) 
      ENDDO
!
!     == THIS TERM IS ONLY PRESENT IN THE NON-COLLINEAR CASE ===================
      IF(NDIMD.EQ.4) THEN
        DO I=1,3
!         == (X,Y,Z)->(2,3,4)
          IP=2+MOD(I,3)
          IPP=2+MOD(I+1,3)
          C(:,:,1+I)=C(:,:,1+I)+CI*(MATMUL(A(:,:,IP),B(:,:,IPP)) &
     &                             -MATMUL(A(:,:,IPP),B(:,:,IP)))
        ENDDO
      END IF
!
!     == FINALLY APPLY FACTOR 1/2 ==============================================
      C(:,:,:)=0.5D0*C(:,:,:)
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINOR$INVERT(NDIMD,NCHI,A,AINV)
!     **************************************************************************
!     ** INVERTS A MATRIX GIVE IN SPINOR REPRESENTATION                       **
!     **    A(S1,S2)=0.5 SUM_{J=0}^3 A(J) SIGMAJ(S1,S2)                       **
!     ** WHERE SIGMAJ IS THE 2X2 UNIT MATRIX FOR J=0 AND                      **
!     **                 THE PAULI MATRIX FOR J>0                             **
!     **                                                                      **
!     **  RATIONALE FOR THE APPARENTLY COMPLICATED TREATMENT:                 **
!     **  BECAUSE THIS INVERSION IS DONE FREQUENTY, I TRY TO KEEP ARRAYS IN   **
!     **  IN THE STACK. FURTHERMORE, I TRY TO MAINTAIN AN ANALOGOUS TREATMENT **
!     **  OF NON-COLLINEAR, COLLINEAR SPIN-POLARIZED AND NON-COLLINEAR CASES. **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      COMPLEX(8),INTENT(IN) :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(OUT):: AINV(NCHI,NCHI,NDIMD)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
      COMPLEX(8)            :: MAT(NCHI,NCHI,NDIMD)
      REAL(8)               :: SVAR
      INTEGER(4)            :: IDIMD,I
!     **************************************************************************
      IF(NDIMD.EQ.1) THEN
         CALL LIB$INVERTC8(NCHI,A(:,:,1),AINV(:,:,1))
         AINV(:,:,1)=4.D0*AINV(:,:,1)
      ELSE IF(NDIMD.EQ.2) THEN
         AINV(:,:,1)=0.5D0*(A(:,:,1)+A(:,:,2))
         AINV(:,:,2)=0.5D0*(A(:,:,1)-A(:,:,2))
         CALL LIB$INVERTC8(NCHI,AINV(:,:,1),MAT(:,:,1))
         CALL LIB$INVERTC8(NCHI,AINV(:,:,2),MAT(:,:,2))
         AINV(:,:,1)=MAT(:,:,1)+MAT(:,:,2)
         AINV(:,:,2)=MAT(:,:,1)-MAT(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
!        =======================================================================
!        == NON-COLLINEAR CASE                                                ==
!        == 4 INVERSIONS + 6 MULTIPLICATIONS + ORDER N**2                     ==
!        =======================================================================
!        == A(T,X,Y,Z)-> AINV(11,21,12,22) =====================================
         MAT=A
         CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,MAT)
!
!        == INVERSION IN UP-DOWN REPRESENTATION ================================
         CALL LIB$INVERTC8(NCHI,MAT(:,:,1),AINV(:,:,1))   ! 1/A11
         CALL LIB$INVERTC8(NCHI,MAT(:,:,4),AINV(:,:,4))   ! 1/A22
         AINV(:,:,3)=-MATMUL(AINV(:,:,1),MAT(:,:,3))      ! -(1/A11)*A12
         AINV(:,:,2)=-MATMUL(AINV(:,:,4),MAT(:,:,2))      ! -(1/A22)*A21
                                                          !A22-A21*(1/A11)*A12
         MAT(:,:,4)=MAT(:,:,4)+MATMUL(MAT(:,:,2),AINV(:,:,3))
                                                          ! A11-A12*(1/A22)*A21
         MAT(:,:,1)=MAT(:,:,1)+MATMUL(MAT(:,:,3),AINV(:,:,2))
                                                  ! B11=1/[A22-A21*(1/A11)*A12]
         CALL LIB$INVERTC8(NCHI,MAT(:,:,4),AINV(:,:,4))
                                                  !B22=1/[A11-A12*(1/A22)*A21]
         CALL LIB$INVERTC8(NCHI,MAT(:,:,1),AINV(:,:,1)) 
         AINV(:,:,3)=MATMUL(AINV(:,:,3),AINV(:,:,4))    ! B12
         AINV(:,:,2)=MATMUL(AINV(:,:,2),AINV(:,:,1)) ! B21
!
!        == MAT(11,21,12,22) -> AINV(T,X,Y,Z) ==================================
         CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,AINV)
         IF(TTEST) THEN
           CALL SPINOR$MATMUL(NDIMD,NCHI,A,AINV,MAT)
           CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,MAT)
           DO IDIMD=1,NDIMD,NDIMD-1 ! (1), (1,2), (1,4)
             DO I=1,NCHI
               MAT(I,I,IDIMD)=MAT(I,I,IDIMD)-CMPLX(1.D0,0.D0)
             ENDDO
           ENDDO
           CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,MAT)
           SVAR=MAXVAL(ABS(MAT))
           IF(SVAR.GT.1.D-6) THEN
             CALL SPINOR_PRINTMATRIX(6,'INVERSION TEST (0)',1,NCHI &
        &                           ,NDIMD,NCHI,MAT)
             CALL ERROR$MSG('INVERSION TEST FAILED') 
             CALL ERROR$R8VAL('MAX DEV',SVAR)
             CALL ERROR$STOP('DMFT_INVERTSPINORMATRIX')
           END IF
         END IF
      ELSE
        CALL ERROR$MSG('ILLEGAL VALUE FOR ARGUMENT NDIMD') 
        CALL ERROR$I4VAL('NDIMD',NDIMD)
        CALL ERROR$STOP('DMFT_INVERTSPINORMATRIX')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINOR$CONVERT(ID,NCHI,NDIMD,A)
!     **************************************************************************
!     **  TRANSFORMATION OF A MATRIX FROM THE UP/DOWN-SPIN REPRESENTATION     **
!     **  (UU,DU,UD,DD) TO THE (TOTAL,X,Y,Z) REPRESENTATION AND BACK.         **
!     **                                                                      **
!     **  CAUTION! THE UP-DOWN REPRESENTATION USES THE INDEX ORDER OF FORTRAN **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4)  ,INTENT(IN)    :: NCHI
      INTEGER(4)  ,INTENT(IN)    :: NDIMD
      COMPLEX(8)  ,INTENT(INOUT) :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8)  ,PARAMETER     :: CI=(0.D0,1.D0)
      COMPLEX(8)                 :: B(NCHI,NCHI,NDIMD)
      LOGICAL(4)                 :: TOUPDN
!     **************************************************************************
      TOUPDN=ID.EQ.'BACK'
      IF(.NOT.(TOUPDN.OR.ID.EQ.'FWRD')) THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED. MUST BE EITHER "FWRD" OR "BACK"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('DMFT_TXYZTOUPDOWN')
      END IF
!
      IF(NDIMD.EQ.1) THEN
        IF(TOUPDN) THEN
          A(:,:,1)=0.5D0*A(:,:,1)
        ELSE
          A(:,:,1)=2.D0*A(:,:,1)
        END IF
      ELSE IF(NDIMD.EQ.2) THEN
        IF(TOUPDN) THEN
          B(:,:,1)=0.5D0*(A(:,:,1)+A(:,:,2))
          B(:,:,2)=0.5D0*(A(:,:,1)-A(:,:,2))
        ELSE
          B(:,:,1)=A(:,:,1)+A(:,:,2)
          B(:,:,2)=A(:,:,1)-A(:,:,2)
        END IF
        A=B
      ELSE IF(NDIMD.EQ.4) THEN
        IF(TOUPDN) THEN ! A(T,X,Y,Z)-> B(UU,UD,DU,DD)
          B(:,:,1)=0.5D0*(A(:,:,1)+A(:,:,4))
          B(:,:,2)=0.5D0*(A(:,:,2)+CI*A(:,:,3))
          B(:,:,3)=0.5D0*(A(:,:,2)-CI*A(:,:,3))
          B(:,:,4)=0.5D0*(A(:,:,1)-A(:,:,4))
        ELSE ! A(UU,UD,DU,DD)-> B(T,X,Y,Z)
          B(:,:,1)=A(:,:,1)+A(:,:,4)
          B(:,:,2)=A(:,:,2)+A(:,:,3)
          B(:,:,3)=-CI*(A(:,:,2)-A(:,:,3))
          B(:,:,4)=A(:,:,1)-A(:,:,4)
        END IF
        A=B
      ELSE
        CALL ERROR$MSG('ILLEGAL VALUE FOR ARGUMENT NDIMD') 
        CALL ERROR$I4VAL('NDIMD',NDIMD)
        CALL ERROR$STOP('SPINOR_CONVERT')
      END IF   
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINOR$PRINTMATRIX(NFIL,NAME,TYPE,I1,I2,NDIMD,NCHI,A)
!     **************************************************************************
!     ** PRINT A SPINOR MATRIX FOR THE SECTION (I1:I2)X(I1:I2) IN THE         **
!     ** REPRESENTATION SPECIFIED BY THE VARIABLE TYPE.                       **
!     ** TYPE MAY BE 'UPDOWN', 'TXYZ' OR 'DIRECT'
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: TYPE  ! MAY BE
      INTEGER(4)  ,INTENT(IN) :: I1,I2
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(IN) :: NDIMD
      INTEGER(4)  ,INTENT(IN) :: NCHI
      COMPLEX(8)  ,INTENT(IN) :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)              :: B(NCHI,NCHI,NDIMD)
!     **************************************************************************
      IF(I2.LT.I1) RETURN
!
!     ==========================================================================
!     == ARGUMENT CHECKS                                                      ==
!     ==========================================================================
      IF(I1.LT.1.OR.I1.GT.NCHI) THEN
        CALL ERROR$MSG('LOWER BOUND I1 OUT OF RANGE')
        CALL ERROR$CHVAL('I1',I1)
        CALL ERROR$CHVAL('I2',I2)
        CALL ERROR$CHVAL('NCHI',NCHI)
        CALL ERROR$CHVAL('NAME',NAME)
        CALL ERROR$STOP('SPINOR$PRINTMATRIX')
      END IF
      IF(I2.LT.1.OR.I2.GT.NCHI) THEN
        CALL ERROR$MSG('UPPER BOUND I2 OUT OF RANGE')
        CALL ERROR$CHVAL('I1',I1)
        CALL ERROR$CHVAL('I2',I2)
        CALL ERROR$CHVAL('NCHI',NCHI)
        CALL ERROR$CHVAL('NAME',NAME)
        CALL ERROR$STOP('SPINOR$PRINTMATRIX')
      END IF
!
!     ==========================================================================
!     == CONVERT MATRIX IF REQUIRED                                           ==
!     ==========================================================================
      B=A
      IF(TYPE.EQ.'UPDOWN') THEN
!       == TRANSFORM FROM (T,X,Y,Z) TO (11,12,21,22) ===========================
        CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,B)
      ELSE IF(TYPE.EQ.'TXYZ'.OR.TYPE.EQ.'DIRECT') THEN
      ELSE
        CALL ERROR$MSG('TYPE ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('TYPE',TYPE)
        CALL ERROR$CHVAL('NAME',NAME)
        CALL ERROR$STOP('SPINOR$PRINTMATRIX')
      END IF
!
!     ==========================================================================
!     == PRINT MATRIX                                                         ==
!     ==========================================================================
      CALL SPINOR_PRINTMATRIX(NFIL,NAME,I1,I2,NDIMD,NCHI,B)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINOR_PRINTMATRIX(NFIL,NAME,I1,I2,NDIMD,NCHI,A)
!     **************************************************************************
!     ** PRINT THE SUBBLOCK I1:I2 X I1:I2 FROM THE MATRIX IN SPINOR           **
!     ** REPRESENTATION                                                       **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: NAME
      INTEGER(4)  ,INTENT(IN) :: NDIMD
      INTEGER(4)  ,INTENT(IN) :: NCHI
      INTEGER(4)  ,INTENT(IN) :: I1,I2
      COMPLEX(8)  ,INTENT(IN) :: A(NCHI,NCHI,NDIMD)
      INTEGER(4)              :: I,IDIMD
!     **************************************************************************
      WRITE(NFIL,FMT='(82("="),T10," ",A," ")')NAME
      DO IDIMD=1,NDIMD
        DO I=I1,I2
          WRITE(NFIL,FMT='("IDIMD=",I1,":",100("(",2F10.5,")"))') &
    &                    IDIMD,A(I,I1:I2,IDIMD)
        ENDDO
      ENDDO
      RETURN
      END
