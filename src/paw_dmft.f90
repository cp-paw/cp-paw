!........1.........2.........3.........4.........5.........6.........7.........8
MODULE DMFT_MODULE
!*******************************************************************************
!** IPROOFCHI(NCHI) IS A POINTER THAT SELECTS A PARTICULAR LOCAL ORBITAL      **
!**                                                                           **
!**  DESCRIPTION OF VARIABLES                                                 **
!**    GLOCLAUR1,1,2 EXPANSION OF SPIN-AVERAGED GLOC IN 1/(CI*HBAR*OMEGA)     **
!**                                                                           **
!*******************************************************************************
TYPE DENMAT_TYPE   ! CONSIDERS ALL ORBITALS ON THIS SITE
INTEGER(4)          :: LNX              !#(RADIAL PARTIAL WAVES)
INTEGER(4),POINTER  :: LOX(:) => NULL() !(LNX)
INTEGER(4),POINTER  :: LMN(:) => NULL() !(NLOC) MAPS LOCAL ORBITALS 
                                        !ONTO PROJECTOR ARRAY
INTEGER(4)          :: LMNX
COMPLEX(8),POINTER  :: RHO(:,:,:) => NULL() ! DENSITY MATRIX
COMPLEX(8),POINTER  :: H(:,:,:)   => NULL() ! HAMILTONIAN FROM DOUBLE COUNTING
END TYPE DENMAT_TYPE
TYPE NATORB_TYPE                             ! NATURAL ORBITALS
COMPLEX(8),POINTER  :: CHIPHI(:,:) => NULL() !(2*NLOC,2*NLOC) <CHI|NATORB>
COMPLEX(8),POINTER  :: PIPHI(:,:)  => NULL() !(2*NLOC,2*NLOC) <PI|NATORB>
END TYPE NATORB_TYPE
TYPE ATOMSET_TYPE
  INTEGER(4)           :: NLOC
  INTEGER(4)           :: ICHI1
  INTEGER(4)           :: ICHI2
  REAL(8)              :: LHFWEIGHT
  REAL(8)   ,POINTER   :: U(:,:,:,:)        => NULL() !(NLOC,NLOC,NLOC,NLOC) 
  REAL(8)   ,POINTER   :: DEDU(:,:,:,:)     => NULL() !(NLOC,NLOC,NLOC,NLOC) 
  COMPLEX(8),POINTER   :: GLOC(:,:,:,:)     => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: GLOCLAUR(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,3)
  COMPLEX(8),POINTER   :: SLOC(:,:,:,:)     => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: SLOCLAUR(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,3)
  COMPLEX(8),POINTER   :: RHO(:,:,:)        => NULL() !LOCAL DENSITY MATRIX
  COMPLEX(8),POINTER   :: SMAT(:,:,:)       => NULL() !OVERLAP OF LOCAL ORBITALS
  TYPE(DENMAT_TYPE)    :: DENMAT
  TYPE(NATORB_TYPE)    :: NATORB
END TYPE ATOMSET_TYPE
TYPE KSET_TYPE
  REAL(8)            :: WKPT
  LOGICAL(4)         :: TADDMINUSK    !ADD  THE TERM FOR -K
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
REAL(8)   ,PARAMETER   :: MU=0.D0       ! CHEMICAL POTENTIAL
!== KSET =======================================================================
TYPE(KSET_TYPE)   ,ALLOCATABLE :: KSET(:)  !(NKPTL)
TYPE(ATOMSET_TYPE),ALLOCATABLE :: ATOMSET(:)  !(NAT)
END MODULE DMFT_MODULE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_INI()
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TINI,NDIM,NDIMD,NSPIN,NKPTL,NB,NAT,NCHI &
     &                       ,NOMEGA,KBT,MU,OMEGA,IPROOFCHI &
     &                       ,KSET,ATOMSET
      USE WAVES_MODULE, ONLY : KMAP,NDIM_W=>NDIM,NKPTL_W=>NKPTL,NSPIN_W=>NSPIN
      USE LMTO_MODULE, ONLY: HYBRIDSETTING,HFWEIGHT
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
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
      PI=4.D0*ATAN(1.D0)
      CALL MPE$QUERY('K',NTASKS_K,THISTASK_K)
      CALL MPE$QUERY('MONOMER',NTASKS_M,THISTASK_M)
!
!     ==========================================================================
!     == HARDWIRED VARIABLES                                                  ==
!     ==========================================================================
      NOMEGA=50
!
!     ==========================================================================
!     == INHERIT KBT FROM OCCUPATIONS OBJECT                                  ==
!     ==   (SPECIFIED IN !CONTROL!MERMIN:T[K])                                ==
!     ==========================================================================
      CALL DYNOCC$GETR8('TEMP',KBT)
      IF(KBT.LT.1.D-5) THEN
        CALL ERROR$MSG('TEMPERATURE MUST BE FINITE')
        CALL ERROR$MSG('SET !CONTROL!MERMIN:T[K]')
        CALL ERROR$R8VAL('KBT',KBT)
        CALL ERROR$STOP('DMFT_INI')
      END IF
PRINT*,'KBT=',KBT,' KBT[EV]=',KBT*27.211D0
!!$KBT=0.333D0*27.211D0
!!$print*,'kbt=',kbt,' kbt[ev]=',kbt*27.211d0
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
!
!       == SAVE BOUNDS TO ATOMSET ==============================================
        NLOC=I2-I1+1
        ATOMSET(IAT)%NLOC=NLOC
        ATOMSET(IAT)%ICHI1=I1
        ATOMSET(IAT)%ICHI2=I2
!
!       == NOW DETERMINE SUBSTRUCTURE DENMAT ==================================
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
!
!       ========================================================================
!       == INHERIT SCREENING FACTOR  FROM LMTO_MODULE                         ==
!       ========================================================================
        IF(HYBRIDSETTING(ISP)%LHFWEIGHT.GE.0.D0) THEN
          ATOMSET(IAT)%LHFWEIGHT=HYBRIDSETTING(ISP)%LHFWEIGHT
        ELSE
          ATOMSET(IAT)%LHFWEIGHT=HFWEIGHT
        END IF
PRINT*,'IAT=',IAT,' LOCAL HFWEIGHT=',ATOMSET(IAT)%LHFWEIGHT
!
!       == ALLOCATE ATOMSET SUBARRAYS ==========================================
        ALLOCATE(ATOMSET(IAT)%U(NLOC,NLOC,NLOC,NLOC))
        ALLOCATE(ATOMSET(IAT)%DEDU(NLOC,NLOC,NLOC,NLOC))
        ALLOCATE(ATOMSET(IAT)%GLOC(NLOC,NLOC,NDIMD,NOMEGA))
        ALLOCATE(ATOMSET(IAT)%GLOCLAUR(NLOC,NLOC,NDIMD,3))
        ALLOCATE(ATOMSET(IAT)%SLOC(NLOC,NLOC,NDIMD,NOMEGA))
        ALLOCATE(ATOMSET(IAT)%SLOCLAUR(NLOC,NLOC,NDIMD,3))
        ALLOCATE(ATOMSET(IAT)%SMAT(LMNX,LMNX,NDIMD))
        ALLOCATE(ATOMSET(IAT)%DENMAT%RHO(LMNX,LMNX,NDIMD))
        ALLOCATE(ATOMSET(IAT)%DENMAT%H(LMNX,LMNX,NDIMD))
        ATOMSET(IAT)%U=0.D0
        ATOMSET(IAT)%GLOC=(0.D0,0.D0)
        ATOMSET(IAT)%GLOCLAUR=(0.D0,0.D0)
        ATOMSET(IAT)%SLOC=(0.D0,0.D0)
        ATOMSET(IAT)%SLOCLAUR=(0.D0,0.D0)
        ATOMSET(IAT)%SMAT=(0.D0,0.D0)
        ATOMSET(IAT)%DENMAT%RHO =(0.D0,0.D0)
        ATOMSET(IAT)%DENMAT%H   =(0.D0,0.D0)
!
        DEALLOCATE(LOX)
        DEALLOCATE(TORB)
      ENDDO
PRINT*,'NCHI ',NCHI
!
!     ==========================================================================
!     == ACCUMULATE IPROOFCHI                                                 ==
!     ==========================================================================
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
        ALLOCATE(KSET(IKPTL)%PIPSI(NDIM,NCHI,NB,NSPIN))
        ALLOCATE(KSET(IKPTL)%RHO(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%H0(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%HRHO(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%SINV(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%SMAT(NCHI,NCHI,NDIMD))
        KSET(IKPTL)%WKPT=0.D0
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
          KSET(IKPTL)%WKPT=WKPT(IKPT)
        END IF
      ENDDO
      IF(IKPTL.NE.NKPTL) THEN
        CALL ERROR$MSG('KMAP INCONSISTENT WITH NKPTL FROM PAW_WAVES')
        CALL ERROR$I4VAL('NKPTL',NKPTL)
        CALL ERROR$I4VAL('IKPTL',IKPTL)
        CALL ERROR$STOP('DMFT_INI')
      END IF
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
      REAL(8)                :: ETOT
      REAL(8)                :: SVAR
      INTEGER(4)             :: ITER
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT$GREEN')
WRITE(*,FMT='(82("="),T20," ENTERING DMFT$GREEN ")')
      CALL DMFT_INI()
!
!     ==========================================================================
!     ==  COLLECT DFT HAMILTONIAN                                             ==
!     ==========================================================================
      CALL DMFT_COLLECTHAMILTONIAN()  
      CALL DMFT_COLLECTFULLDENMAT()  
!
!     ==========================================================================
!     ==  OBTAIN BARE U-TENSOR FROM LMTO OBJECT.                              ==
!     ==     SHAPE OF ORBITALS ARE DEFINED BY NPRO AND ATOMIC STRUCTRE        ==
!     ==     ORBITALS ARE SELECTED BY TORB.                                   ==
!     ==========================================================================
      CALL DMFT_UTENSOR() 
!
!     ==========================================================================
!     ==  OVERLAP MATRIX IN K-SPACE                                           ==
!     ==========================================================================
      CALL DMFT_SMAT()
!
!     ==========================================================================
!     ==  CONSTRUCT LOCAL NATURAL ORBITALS                                    ==
!     ==========================================================================
      CALL DMFT_NATORB()
!
!     ==========================================================================
!     ==  CONSTRUCT NON-INTERACTING HAMILTONIAN THAT PRODUCES THE CORRECT     ==
!     ==  ONE-PARTICLE DENSITY MATRIX                                         ==
!     ==========================================================================
      CALL DMFT_HRHO()
      CALL DMFT_CONSTRAINTS('HRHO')
!
!     ==========================================================================
!     == ITERATION TO ENFORCE CONSTRAINTS                                     ==
!     ==========================================================================
      DO ITER=1,3
WRITE(*,FMT='(82("="),T20," ITERATION ",I5)')ITER
        CALL DMFT_GLOC_WITHATOMSET() 
!
!       ========================================================================
!       ==  CALL THE SOLVER                                                   ==
!       ========================================================================
        CALL DMFT_SOLVER(ETOT) 
PRINT*,'ETOT AFTER DMFT_SOLVER ',ETOT
!
!       ========================================================================
!       ==  CONSTRAINTS                                                       ==
!       ========================================================================
        CALL DMFT_CONSTRAINTS('H0')
!        CALL DMFT_CONSTRAINTS_WITHKSET('H0')
        PRINT*,'ITERATION COMPLETED ',ITER
      ENDDO ! END OF LOOP OVER ITERATIONS TO ENFORCE CONSTRAINT
!
!     ==========================================================================
!     ==  CALCULATE MISSING TOTAL ENERGY CONTRIBUTION                         ==
!     ==========================================================================
      CALL DMFT_DETOT(SVAR)
      ETOT=ETOT+SVAR      
      CALL ENERGYLIST$SET('DMFT INTERFACE',ETOT)
      CALL ENERGYLIST$ADD('LOCAL CORRELATION',ETOT)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',ETOT)
!
!     ==========================================================================
!     ==  ADD HAMILTONIAN THIS%HTBC                                           ==
!     ==========================================================================
      CALL DMFT_ADDTOHPSI()
WRITE(*,FMT='(82("="),T20," LEAVING DMFT$GREEN ")')
                                       CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SMAT()
!     **************************************************************************
!     **  DETERMINE OVERLAP MATRIX SMAT(NCHI,NCHI,ISPIN) OF LOCAL ORBITALS,   **
!     **  THE K-DEPENDENT INVERSE SINV(NCHI,NCHI,IKPT,ISPIN)                  **
!     **                                                                      **
!     **   SINV(K) =SUM_N <PI|PSI(N,K)><PSI(N,K)|PI>                          **
!     **   WHERE THE PI ARE PROJECTOR FUNCTIONS AND CHI ARE ORBITALS          **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NSPIN,NDIM,NDIMD,KSET
      IMPLICIT NONE
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
      INTEGER(4)            :: IKPT,ISPIN,IDIM1,IDIM2,IDIMD
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
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_NATORB()
!     **************************************************************************
!     ** NATURAL ORBITALS IN A TWO COMPONENT PAULI-SPINOR REPRESENTATION      **
!     **  USES KSET%RHO AND KSET%SINV 
!     **************************************************************************
      USE DMFT_MODULE, ONLY : NAT,NKPTL,NDIMD,ATOMSET,KSET
      IMPLICIT NONE
      LOGICAL(4)             :: NATORB=.TRUE.
      INTEGER(4)             :: NLOC
      COMPLEX(8),ALLOCATABLE :: RHO(:,:,:)
      COMPLEX(8),ALLOCATABLE :: SINV(:,:,:)
      COMPLEX(8),ALLOCATABLE :: RHOB(:,:)
      COMPLEX(8),ALLOCATABLE :: SINVB(:,:)
      REAL(8)   ,ALLOCATABLE :: X(:)
      REAL(8)                :: SVAR
      INTEGER(4)             :: I1,I2
      INTEGER(4)             :: IAT,ISPIN,IKPT,I
!     **************************************************************************
      DO IAT=1,NAT
        NLOC=ATOMSET(IAT)%NLOC
        IF(NLOC.LE.0) CYCLE
!
!       ========================================================================
!       == ALLOCATE ARRAYS                                                    ==
!       ========================================================================
        IF(.NOT.ASSOCIATED(ATOMSET(IAT)%NATORB%CHIPHI)) THEN
          ALLOCATE(ATOMSET(IAT)%NATORB%CHIPHI(2*NLOC,2*NLOC))
        END IF
        IF(.NOT.ASSOCIATED(ATOMSET(IAT)%NATORB%PIPHI)) THEN
          ALLOCATE(ATOMSET(IAT)%NATORB%PIPHI(2*NLOC,2*NLOC))
        END IF
!
!       ========================================================================
!       == OBTAIN LOCAL DENSITY MATRIX AND INVERSE OVERLAP FROM KSET          ==
!       ========================================================================
        ALLOCATE(RHO(NLOC,NLOC,NDIMD))
        ALLOCATE(SINV(NLOC,NLOC,NDIMD))
        RHO(:,:,:)=(0.D0,0.D0)
        SINV(:,:,:)=(0.D0,0.D0)
        DO IKPT=1,NKPTL
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          SVAR=KSET(IKPT)%WKPT
          RHO =RHO +SVAR*KSET(IKPT)%RHO(I1:I2,I1:I2,:)
          SINV=SINV+SVAR*KSET(IKPT)%SINV(I1:I2,I1:I2,:)
        ENDDO
!
!       ========================================================================
!       == CONVERT TO UP-DOWN REPRESENTATION                                  ==
!       ========================================================================
        CALL SPINOR$CONVERT('BACK',NLOC,NDIMD,RHO) 
        CALL SPINOR$CONVERT('BACK',NLOC,NDIMD,SINV) 
!
!       ========================================================================
!       == DETERMINE OCCUPATIONS AND EIGENVALUES                              ==
!       ========================================================================
        ATOMSET(IAT)%NATORB%CHIPHI=(0.D0,0.D0)
        ATOMSET(IAT)%NATORB%PIPHI=(0.D0,0.D0)
        IF(NDIMD.LE.2) THEN
          ALLOCATE(X(NLOC))
          DO ISPIN=1,NDIMD
            I1=1+NLOC*(ISPIN-1)
            I2=NLOC*ISPIN
            IF(NATORB) THEN
!             == A*U=B*U*X  WITH U^DAGGER*S*U=1 ================================
!             == X IS THE OCCUPATION ===========================================
              CALL LIB$GENERALEIGENVALUEC8(NLOC,RHO(:,:,ISPIN),SINV(:,:,ISPIN) &
      &                              ,X,ATOMSET(IAT)%NATORB%CHIPHI(I1:I2,I1:I2))
              ATOMSET(IAT)%NATORB%PIPHI(I1:I2,I1:I2) &
      &         =MATMUL(SINV(:,:,ISPIN),ATOMSET(IAT)%NATORB%CHIPHI(I1:I2,I1:I2))
            ELSE
!CALL SPINOR_PRINTMATRIX(6,'XXX',1,5,2,NLOC,SINV)
!             == EIGENSTATES OF SINV ===========================================
              CALL LIB$DIAGC8(NLOC,SINV(:,:,ISPIN),X &
      &                                ,ATOMSET(IAT)%NATORB%CHIPHI(I1:I2,I1:I2))
              ATOMSET(IAT)%NATORB%PIPHI(I1:I2,I1:I2) &
      &                                 =ATOMSET(IAT)%NATORB%CHIPHI(I1:I2,I1:I2)
              DO I=I1,I2
                SVAR=SQRT(X(I-I1+1))
                ATOMSET(IAT)%NATORB%PIPHI(I1:I2,I) &
      &                                 =ATOMSET(IAT)%NATORB%PIPHI(I1:I2,I)*SVAR
                ATOMSET(IAT)%NATORB%CHIPHI(I1:I2,I) &
      &                                =ATOMSET(IAT)%NATORB%CHIPHI(I1:I2,I)/SVAR
              ENDDO
            END IF
          ENDDO
          DEALLOCATE(X)
!         == COMPLETE SECOND SPIN DIRECTION FOR NON-SPIN-POLARIZED CASE ========
          IF(NDIMD.EQ.1) THEN
            I1=NLOC+1
            I2=2*NLOC
            ATOMSET(IAT)%NATORB%CHIPHI(I1:I2,I1:I2) &
      &                               =ATOMSET(IAT)%NATORB%CHIPHI(1:NLOC,1:NLOC)
            ATOMSET(IAT)%NATORB%PIPHI(I1:I2,I1:I2) &
      &                               =ATOMSET(IAT)%NATORB%PIPHI(1:NLOC,1:NLOC)
          END IF
        ELSE 
          ALLOCATE(RHOB(2*NLOC,2*NLOC))
          ALLOCATE(SINVB(2*NLOC,2*NLOC))
          ALLOCATE(X(2*NLOC))
          RHOB(  :NLOC,  :NLOC)=RHO(:,:,1)
          RHOB(NLOC+1:,  :NLOC)=RHO(:,:,2)
          RHOB(  :NLOC,NLOC+1:)=RHO(:,:,3)
          RHOB(NLOC+1:,NLOC+1:)=RHO(:,:,4)
          SINVB(  :NLOC,  :NLOC)=SINV(:,:,1)
          SINVB(NLOC+1:,  :NLOC)=SINV(:,:,2)
          SINVB(  :NLOC,NLOC+1:)=SINV(:,:,3)
          SINVB(NLOC+1:,NLOC+1:)=SINV(:,:,4)
          IF(NATORB) THEN
!           == A*U=B*U*F  WITH U^DAGGER*S*U=1 ==================================
            CALL LIB$GENERALEIGENVALUEC8(2*NLOC,RHOB,SINVB &
      &                                           ,X,ATOMSET(IAT)%NATORB%CHIPHI)
            ATOMSET(IAT)%NATORB%PIPHI=MATMUL(SINVB,ATOMSET(IAT)%NATORB%CHIPHI)
          ELSE
!           == EIGENSTATES OF SINV ===========================================
            CALL LIB$DIAGC8(2*NLOC,SINVB,X,ATOMSET(IAT)%NATORB%CHIPHI)
            ATOMSET(IAT)%NATORB%PIPHI=ATOMSET(IAT)%NATORB%CHIPHI
            DO I=1,2*NLOC
              ATOMSET(IAT)%NATORB%PIPHI(I,I1:I2) &
      &                               =ATOMSET(IAT)%NATORB%PIPHI(I,:)*SQRT(X(I))
              ATOMSET(IAT)%NATORB%CHIPHI(I1:I2,I) &
      &                              =ATOMSET(IAT)%NATORB%CHIPHI(:,I)/SQRT(X(I))
            ENDDO
          END IF
          DEALLOCATE(X)
          DEALLOCATE(RHOB)
          DEALLOCATE(SINVB)
        END IF
        DEALLOCATE(RHO)
        DEALLOCATE(SINV)
      ENDDO ! END OF LOOP OVER ATOMS
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_COLLECTHAMILTONIAN()
!     **************************************************************************
!     ** COLLECTS THE HAMILTONIAN AND STORES IT ON THE MODULE                 **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,NDIMD,IPROOFCHI &
     &                      ,NAT,KSET,ATOMSET
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER   :: TTEST=.FALSE.
      COMPLEX(8),ALLOCATABLE :: RHO(:,:,:) !(NCHI,NCHI,NDIMD)
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IPRO,IB,J,IAT,IDIM1,IDIM2
      INTEGER(4)             :: IDIMD
      INTEGER(4)             :: I1,I2
      REAL(8)                :: F(NB,NKPTL,NSPIN)
      COMPLEX(8)             :: MAT(NCHI,NCHI)
      COMPLEX(8)             :: CSVAR
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                      CALL TRACE$PUSH('DMFT_COLLECTHAMILTONIAN')
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
            CALL ERROR$STOP('DMFT_COLLECTHAMILTONIAN')
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
!
!     ==  COLLECT OCCUPATIONS W/O K-POINT WEIGHT ===============================
      CALL DMFT_COLLECTOCCUPATIONS(NKPTL,NSPIN,NB,F)
!
!     == SUM UP DENSITY MATRIX =================================================
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
          RHO=RHO+KSET(IKPT)%WKPT*KSET(IKPT)%RHO
        ENDDO
!       == MAKE HERMITEAN ======================================================
        DO IDIMD=1,NDIMD
          RHO(:,:,IDIMD)=0.5D0*(RHO(:,:,IDIMD)+CONJG(TRANSPOSE(RHO(:,:,IDIMD))))
        ENDDO
!       == PRINT ===============================================================
        PRINT*,' IN DMFT_COLLECTHAMILTONIAN'
        DO IAT=1,NAT
          IF(ATOMSET(IAT)%NLOC.LE.0) CYCLE
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          CALL SPINOR$PRINTMATRIX(6,'DENSITY MATRIX(UPDOWN)','UPDOWN' &
      &                          ,I1,I2,NDIMD,NCHI,RHO)
          CALL SPINOR$PRINTMATRIX(6,'DENSITY MATRIX(TXYZ)','TXYZ' &
      &                          ,I1,I2,NDIMD,NCHI,RHO)
        ENDDO
        DEALLOCATE(RHO)
      END IF
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_COLLECTFULLDENMAT()
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
        CALL ERROR$STOP('DMFT_COLLECTFULLDENMAT')
      END IF
!
      CALL DYNOCC$GETI4('NSPIN',NSPIN_)
      IF(NSPIN_.NE.NSPIN) THEN
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$I4VAL('NSPIN_',NSPIN_)
        CALL ERROR$STOP('DMFT_COLLECTFULLDENMAT')
      END IF
      CALL WAVES_DYNOCCGETR8A('OCC',NB*NKPTL*NSPIN,OCC)
!
!     ==========================================================================
!     == SET DENSITY MATRIX TO ZERO
!     ==========================================================================
      DO IAT=1,NAT
        ATOMSET(IAT)%DENMAT%RHO=(0.D0,0.D0)
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
            CALL ERROR$STOP('DMFT_COLLECTFULLDENMAT')
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
                      ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
     &                               =ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
     &                               +REAL(THIS%TBC(IDIM2,IBH,I1:I2))*CSVAR
!
                      CSVAR=AIMAG(THIS%TBC(IDIM1,IBH,I))*OCC(2*IBH,IKPT,ISPIN)
                      ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
     &                             =ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
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
                      ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
     &                             =ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
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
            ATOMSET(IAT)%DENMAT%RHO(:,:,IDIMD) &
     &                           =REAL(ATOMSET(IAT)%DENMAT%RHO(:,:,IDIMD))
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
        CALL SPINOR$CONVERT('FWRD',LMNX,NDIMD,ATOMSET(IAT)%DENMAT%RHO)
        DO IDIMD=1,NDIMD
          ATOMSET(IAT)%DENMAT%RHO(:,:,IDIMD) &
    &              =0.5D0*(ATOMSET(IAT)%DENMAT%RHO(:,:,IDIMD) &
    &                     +CONJG(TRANSPOSE(ATOMSET(IAT)%DENMAT%RHO(:,:,IDIMD))))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  PRINT
!     ==========================================================================
      IF(TPRINT) THEN
        PRINT*,'DENSITY MATRIX REPOST FROM DMFT_COLLECTFULLDENMAT'
        DO IAT=1,NAT
          IF(ATOMSET(IAT)%NLOC.LE.0) CYCLE
          LMNX=ATOMSET(IAT)%DENMAT%LMNX
          WRITE(*,FMT='(82("="),T10," DENSITY MATRIX FOR ATOM ",I3," ")')IAT
          DO IDIMD=1,NDIMD
            DO LMN=1,LMNX
              WRITE(*,FMT='("IDIMD=",I1,":",100("(",2F10.5,")"))')IDIMD &
       &              ,ATOMSET(IAT)%DENMAT%RHO(LMN,:,IDIMD)
            ENDDO
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_COLLECTOCCUPATIONS(NKPTL,NSPIN,NB,F)
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
        CALL ERROR$STOP('DMFT_COLLECTOCCUPATIONS')
      END IF
!
      CALL DYNOCC$GETI4('NSPIN',NSPIN_)
      IF(NSPIN_.NE.NSPIN) THEN
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$I4VAL('NSPIN_',NSPIN_)
        CALL ERROR$STOP('DMFT_COLLECTOCCUPATIONS')
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_HRHO()
!     **************************************************************************
!     **  CALCULATES A STATIC HAMILTONIAN HRHO SUCH THAT THE RESULTING        **
!     **  DENSITY MATRIX IS RHO=[1+EXP(BETA*HRHO)]^(-1)                       **
!     **************************************************************************
      USE DMFT_MODULE ,ONLY: NKPTL,NCHI,NDIMD,NOMEGA,OMEGA,KBT,KSET
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      REAL(8)   ,PARAMETER   :: TOL=1.D-5  ! TOLERANCE FOR THE DENSITY MATRIX
      COMPLEX(8),ALLOCATABLE :: MAT(:,:,:)
      COMPLEX(8),ALLOCATABLE :: A(:,:)
      COMPLEX(8),ALLOCATABLE :: B(:,:)
      REAL(8)   ,ALLOCATABLE :: F(:)   !OCCUPATIONS
      COMPLEX(8),ALLOCATABLE :: U(:,:)
      INTEGER(4)             :: N
      REAL(8)                :: SVAR
      INTEGER(4)             :: IKPT,ISPIN,I
!     **************************************************************************
      IF(NDIMD.LE.2) THEN  ! NON-SPIN AND COLLINEAR
        N=NCHI
      ELSE IF(NDIMD.EQ.4) THEN ! NON-COLLINEAR
        N=2*NCHI
        ALLOCATE(A(N,N))
        ALLOCATE(MAT(NCHI,NCHI,NDIMD))
      ELSE
        CALL ERROR$STOP('DMFT_HRHO')
      END IF
      ALLOCATE(F(N))
      ALLOCATE(U(N,N))
      ALLOCATE(B(N,N))
!
!     ==========================================================================
!     == CALCULATE HRHO                                                       ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        IF(NDIMD.LE.2) THEN
!         == CONVERT TO UP-DOWN REPRESENTATION =================================
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,KSET(IKPT)%RHO) 
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,KSET(IKPT)%SINV) 
          DO ISPIN=1,NDIMD
!           == A*U=B*U*F  WITH U^DAGGER*S*U=1 ==================================
            CALL LIB$GENERALEIGENVALUEC8(NCHI,KSET(IKPT)%RHO(:,:,ISPIN) &
     &                                       ,KSET(IKPT)%SINV(:,:,ISPIN),F,U)
            B=TRANSPOSE(CONJG(U))
            DO I=1,NCHI
              CALL DMFT_EOFF(KBT,F(I),SVAR)
              B(I,:)=SVAR*B(I,:)
            ENDDO
            KSET(IKPT)%HRHO(:,:,ISPIN)=MATMUL(U,B)
          ENDDO
          CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%RHO)
          CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%SINV)
          CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%HRHO)
        ELSE
          ALLOCATE(MAT(NCHI,NCHI,NDIMD))
          MAT=KSET(IKPT)%RHO
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,MAT) 
          A(  :NCHI,  :NCHI)=MAT(:,:,1)
          A(NCHI+1:,  :NCHI)=MAT(:,:,2)
          A(  :NCHI,NCHI+1:)=MAT(:,:,3)
          A(NCHI+1:,NCHI+1:)=MAT(:,:,4)
          MAT=KSET(IKPT)%SINV
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,MAT) 
          B(  :NCHI,  :NCHI)=MAT(:,:,1)
          B(NCHI+1:,  :NCHI)=MAT(:,:,2)
          B(  :NCHI,NCHI+1:)=MAT(:,:,3)
          B(NCHI+1:,NCHI+1:)=MAT(:,:,4)
!         == A*U=B*U*F  WITH U^DAGGER*S*U=1 ====================================
          CALL LIB$GENERALEIGENVALUEC8(N,A,B,F,U)
          B=TRANSPOSE(CONJG(U))
          DO I=1,N
            CALL DMFT_EOFF(KBT,F(I),SVAR)
            B(I,:)=SVAR*B(I,:)
          ENDDO
          A=MATMUL(U,B)
          MAT(:,:,1)=A(  :NCHI,  :NCHI)
          MAT(:,:,2)=A(NCHI+1:,  :NCHI)
          MAT(:,:,3)=A(  :NCHI,NCHI+1:)
          MAT(:,:,4)=A(NCHI+1:,NCHI+1:)
          CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,MAT) 
          KSET(IKPT)%HRHO=MAT
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  TEST                                                                ==
!     ==========================================================================
      IF(TTEST) THEN 
        IF(.NOT.ALLOCATED(MAT))ALLOCATE(MAT(NCHI,NCHI,NDIMD))
        DO IKPT=1,NKPTL
!         == CALCULATE THE DENSITY MATRIX FROM THE GREENS FUNCTION WITH HRHO ===
          CALL DMFT_HRHO_TEST(NCHI,NDIMD,NOMEGA,OMEGA,KBT &
     &                       ,KSET(IKPT)%SMAT,KSET(IKPT)%HRHO,MAT) !MAT=RHO
!         == CALCULATE DEVIATION FROM TARGET DENSITY MATRIX ====================
          MAT=MAT-KSET(IKPT)%RHO
          IF(MAXVAL(ABS(MAT)).GT.TOL) THEN
            CALL SPINOR_PRINTMATRIX(6,'RHO[HRHO]-RHO',1,NCHI &
     &                           ,NDIMD,NCHI,MAT)
            CALL ERROR$MSG('TEST OF HRHO FAILED')
            CALL ERROR$R8VAL('MAX DEV OF RHO',MAXVAL(ABS(MAT)))
            CALL ERROR$STOP('DMFT_HRHO')
          END IF
        ENDDO
      END IF
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
      SUBROUTINE DMFT_HRHO_TEST(NCHI,NDIMD,NOMEGA,OMEGA,KBT,S,H,RHO)
!     **************************************************************************
!     ** CALCULATES THE DENSITY MATRIX FROM THE GREENS FUNCTION FOR           **
!     ** A SPECIFIED HAMILTONIAN H AND OVERLAP MATRIX S                       **
!     **                                                                      **
!     **  DOES NOT REFER TO ANY MODULES                                       **
!     **  DOES OT KNOW K-POINTS                                               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NCHI
      INTEGER(4),INTENT(IN)    :: NDIMD
      INTEGER(4),INTENT(IN)    :: NOMEGA
      REAL(8)   ,INTENT(IN)    :: KBT
      REAL(8)   ,INTENT(IN)    :: OMEGA(NOMEGA)        ! MATSUBARA FREQUENCIES
      COMPLEX(8),INTENT(IN)    :: S(NCHI,NCHI,NDIMD)   ! OVERLAP MATRIX
      COMPLEX(8),INTENT(IN)    :: H(NCHI,NCHI,NDIMD)   ! HAMILTON MATRIX
      COMPLEX(8),INTENT(OUT)   :: RHO(NCHI,NCHI,NDIMD) !DENSITY MATRIX
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: SINV(NCHI,NCHI,NDIMD)
      REAL(8)                  :: FN(3)
      INTEGER(4)               :: NU,IDIMD
!     **************************************************************************
!       
!     ==========================================================================
!     ==  MATSUBARA SUM                                                       ==
!     ==========================================================================
      RHO=(0.D0,0.D0)
      DO NU=1,NOMEGA
!       == CONSTRUCT LATTICE GREENS FUNCTION ===================================
        MAT1=CI*OMEGA(NU)*S-H
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,MAT2)  !MAT2=G
        RHO=RHO+KBT*MAT2
      ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!     == INCLUDE NEGATIVE FREQUENCIES ==========================================
      DO IDIMD=1,NDIMD
        RHO(:,:,IDIMD)=RHO(:,:,IDIMD)+CONJG(TRANSPOSE(RHO(:,:,IDIMD)))
      ENDDO
!     
!     ==========================================================================
!     ==  REGULARIZATION                                                      ==
!     ==========================================================================
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,3,FN)
      CALL SPINOR$INVERT(NDIMD,NCHI,S,SINV)
!     == GLAUR1=SINV ===========================================================
      RHO=RHO+FN(1)*SINV
!     == GLAUR2=SINV*H*SINV ====================================================
      CALL SPINOR$MATMUL(NDIMD,NCHI,H,SINV,MAT1)
      CALL SPINOR$MATMUL(NDIMD,NCHI,SINV,MAT1,MAT2)
      RHO=RHO+FN(2)*MAT2
!     == GLAUR3=SINV*H*SINV*H*SINV =============================================
      CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,H,MAT1)
      CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,SINV,MAT2)
      RHO=RHO+FN(3)*MAT2
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVER(ETOT)
!     **************************************************************************
!     ** MIMICKS A DMFT SOLVER                                                **
!     **  E=1/BETA * SUM_NU E(IOMEGA_NU0+) TR[ DV * GLOC^\DAGGER]             **
!     **************************************************************************
      USE DMFT_MODULE ,ONLY: NAT,NDIMD,NOMEGA,KBT,ATOMSET
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: ETOT
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      REAL(8)                :: EV
      REAL(8)                :: PHILW
      REAL(8)                :: LHFWEIGHT
      REAL(8)                :: EDC
      INTEGER(4)             :: NLOC  !#(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4)             :: LMNX  !#(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4)             :: IAT
      INTEGER(4)             :: IDIMD,I,J,LMN1,LMN2,NU
      COMPLEX(8),ALLOCATABLE :: RHO(:,:,:)
      COMPLEX(8),ALLOCATABLE :: HAM(:,:,:)
!     **************************************************************************
                                      CALL TRACE$PUSH('DMFT_SOLVER')
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
        ATOMSET(IAT)%DENMAT%H=(0.D0,0.D0)
        ATOMSET(IAT)%SLOC=(0.D0,0.D0)
        ATOMSET(IAT)%SLOCLAUR=(0.D0,0.D0)
!
!       ========================================================================
!       == DETERMINE HARTREE FOCK CONTRIBUTION                                ==
!       ========================================================================
        ALLOCATE(RHO(NLOC,NLOC,NDIMD))
        ALLOCATE(HAM(NLOC,NLOC,NDIMD))
        DO I=1,NLOC
          DO J=1,NLOC
            LMN1=ATOMSET(IAT)%DENMAT%LMN(I)
            LMN2=ATOMSET(IAT)%DENMAT%LMN(J)
            RHO(I,J,:)=ATOMSET(IAT)%DENMAT%RHO(LMN1,LMN2,:)
          ENDDO
        ENDDO
        CALL DMFT_FOCK(NLOC,NDIMD,RHO,ATOMSET(IAT)%U,PHILW,HAM)
        ETOT=ETOT+PHILW  !SCREENING DONE ALREADY IN U-TENSOR
PRINT*,'AFTER FOCK ',IAT,ETOT,PHILW
        DO I=1,NLOC
          DO J=1,NLOC
            LMN1=ATOMSET(IAT)%DENMAT%LMN(I)
            LMN2=ATOMSET(IAT)%DENMAT%LMN(J)
            ATOMSET(IAT)%DENMAT%H(LMN1,LMN2,:) &
     &                            =ATOMSET(IAT)%DENMAT%H(LMN1,LMN2,:)+HAM(I,J,:)
          ENDDO
        ENDDO
!
!       == PRINT IF DESIRED ====================================================
        IF(TPRINT) THEN
          WRITE(*,FMT='(82("="),T10," ",A," IAT=",I4,"  ")') &
       &                        'FOCK HAMILTONIAN',IAT
          DO IDIMD=1,NDIMD
            DO I=1,NLOC
              WRITE(*,FMT='(4I5,100F10.5)')IAT,NU,IDIMD,I,HAM(I,:,IDIMD)
            ENDDO
          ENDDO
        END IF
!
        DEALLOCATE(RHO)
        DEALLOCATE(HAM)
!
!       ========================================================================
!       == SUBTRACT DOUBLE COUNTING TERM                                      ==
!       ========================================================================
!       == COLLECT LOCAL HF WEIGHT =============================================
        LHFWEIGHT=ATOMSET(IAT)%LHFWEIGHT
        LMNX=ATOMSET(IAT)%DENMAT%LMNX
        ALLOCATE(HAM(LMNX,LMNX,NDIMD))
        CALL DMFT_DC(IAT,LMNX,NDIMD,ATOMSET(IAT)%DENMAT%RHO,EDC,HAM)
        ETOT=ETOT+LHFWEIGHT*EDC
PRINT*,'AFTER DC ',IAT,ETOT,LHFWEIGHT*EDC
        ATOMSET(IAT)%DENMAT%H=ATOMSET(IAT)%DENMAT%H+LHFWEIGHT*HAM
!
!       == PRINT IF DESIRED ====================================================
        IF(TPRINT) THEN
          WRITE(*,FMT='(82("="),T10," ",A," IAT=",I4,"  ")') &
       &                        'DBLE-CNTNG HAMILTONIAN',IAT
          DO IDIMD=1,NDIMD
            DO I=1,LMNX
              WRITE(*,FMT='(4I5,100F10.5)')IAT,NU,IDIMD,I,HAM(I,:,IDIMD) 
            ENDDO
          ENDDO
        END IF
!
        DEALLOCATE(HAM)
!
!       ========================================================================
!       == ADD DYNAMIC CONTRIBUTIONS                                          ==
!       == HFWEIGHT IS ABSORBED IN THE LOCAL U-TENSOR
!       ========================================================================
        CALL DMFT_DYNAMICSOLVER(NLOC,NDIMD,NOMEGA,KBT &
     &                         ,ATOMSET(IAT)%GLOC,ATOMSET(IAT)%GLOCLAUR &
     &                         ,ATOMSET(IAT)%U,PHILW &
     &                         ,ATOMSET(IAT)%SLOC,ATOMSET(IAT)%SLOCLAUR &
     &                         ,ATOMSET(IAT)%DEDU &
     &                         ,ATOMSET(IAT)%NATORB%PIPHI &
     &                         ,ATOMSET(IAT)%NATORB%CHIPHI)
        ETOT=ETOT+PHILW !SCREENING ALREADY CONTAINED IN U-TENSOR
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
      SUBROUTINE DMFT_DETOT(ETOT)
!     **************************************************************************
!     ** CALCULATES THE ENERGY CONTRIBUTION FROM THE LOGARISHM TERM ETC       **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NAT,NKPTL,NCHI,NDIMD,NOMEGA,OMEGA,KBT,MU &
     &                      ,ATOMSET,KSET
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: ETOT
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      COMPLEX(8)             :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)             :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)             :: MAT3(NCHI,NCHI,NDIMD)
      REAL(8)                :: VEC(NCHI)
      COMPLEX(8),ALLOCATABLE :: BMAT1(:,:)
      COMPLEX(8),ALLOCATABLE :: BMAT2(:,:)
      REAL(8)   ,ALLOCATABLE :: BVEC(:)
      COMPLEX(8)             :: GRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8)             :: GFULL(NCHI,NCHI,NDIMD)
      REAL(8)                :: XSUM
      REAL(8)                :: ETOT1
      INTEGER(4)             :: NU,I,IAT,I1,I2,IDIMD,IKPT
!     **************************************************************************
      IF(NDIMD.EQ.4) THEN
        ALLOCATE(BMAT1(2*NCHI,2*NCHI))
        ALLOCATE(BMAT2(2*NCHI,2*NCHI))
        ALLOCATE(BVEC(2*NCHI))
      END IF
      ETOT=0.D0
      DO IKPT=1,NKPTL
        ETOT1=0.D0
        DO NU=1,NOMEGA
!         == CALCULATE GRHO ====================================================
          MAT1=(CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT-KSET(IKPT)%HRHO
          CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,GRHO)
!     
!         == HPRIME+SIGMA-HRHO =================================================
          MAT2=KSET(IKPT)%H0-KSET(IKPT)%HRHO
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            MAT2(I1:I2,I1:I2,:)=MAT2(I1:I2,I1:I2,:)+ATOMSET(IAT)%SLOC(:,:,:,NU)
          ENDDO
!     
!         == CALCULATE G =======================================================
          MAT1=MAT1-MAT2
          CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,GFULL)
!     
!         ======================================================================
!         == CONSTRAINT TERM
!         ======================================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,GRHO-GFULL &
    &                       ,KSET(IKPT)%H0-KSET(IKPT)%HRHO,MAT1) 
          XSUM=0.D0
          DO I=1,NCHI
            XSUM=XSUM+2.D0*REAL(MAT1(I,I,1))  ! FACTOR TWO COMES FROM (-OMEGA)
          ENDDO
          ETOT1=ETOT1+XSUM
!     
!         ======================================================================
!         == (HPRIME+SIGMA-HRHO)*G                                            ==
!         ======================================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,GFULL,MAT1) 
          XSUM=0.D0
          DO I=1,NCHI
            XSUM=XSUM+2.D0*REAL(MAT1(I,I,1))  ! FACTOR TWO COMES FROM (-OMEGA)
          ENDDO
          ETOT1=ETOT1+XSUM
!     
!         ======================================================================
!         == LOGARITHM TERM                                                   ==
!         ==                                                                  ==
!         == I AM USING THE TRICK OF DAHLEN, PHYS.REV.A 73,P12511 (2006) EQ. B7=
!         == (THE MATH BEHIND THE TRICK NEEDS TO BE CHECKED).                 ==
!         ==                                                                  ==
!         == NOT USED IS THE POTENTIALLY MORE EFFICIENT ROUTE USING           ==
!         == TR[LN(A)]=LN(DET[A]). ACCORDING TO ROBERT SCHADE THE DETERMINANT ==
!         == IS CALCULATED EFFICIENTLY VIA THE LU DECOMPOSITION IN LAPACK     ==
!         ======================================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,GRHO,MAT2,MAT1) !MAT1=GRHO*MAT1
          DO I=1,NCHI
            MAT1(I,I,1)=MAT1(I,I,1)+(2.D0,0.D0) !(UNIT IN SPINOR REPRESENTATION)
          ENDDO
          DO IDIMD=1,NDIMD  ! MAT2(OMEGA)=MAT1(-OMEGA)
            MAT2(:,:,IDIMD)=TRANSPOSE(CONJG(MAT1(:,:,IDIMD)))
          ENDDO
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,MAT2,MAT3) !MAT1=GRHO*MAT1
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,MAT3) !MAT1=GRHO*MAT1
          IF(NDIMD.EQ.4) THEN
            BMAT1(:NCHI,:NCHI)    =MAT3(:,:,1)
            BMAT1(NCHI+1:,:NCHI)  =MAT3(:,:,2)
            BMAT1(:NCHI,NCHI+1:)  =MAT3(:,:,3)
            BMAT1(NCHI+1:,NCHI+1:)=MAT3(:,:,4)
            CALL LIB$DIAGC8(2*NCHI,BMAT1,BVEC,BMAT2)
            XSUM=0.D0
            DO I=1,2*NCHI
              XSUM=XSUM+LOG(BVEC(I))
            ENDDO
            ETOT1=ETOT1+XSUM
          ELSE
            XSUM=0.D0
            DO IDIMD=1,NDIMD
              CALL LIB$DIAGC8(NCHI,MAT3(:,:,IDIMD),VEC,MAT1(:,:,IDIMD))
              DO I=1,NCHI
                XSUM=XSUM+LOG(VEC(I))  !SLOW BUT AVOIDS OVERFLOWS              
              ENDDO
            ENDDO
            ETOT1=ETOT1+XSUM
          END IF
        ENDDO ! END OF LOOP OVER MATSUBARA FREQUENCIES
        ETOT=ETOT+KSET(IKPT)%WKPT*ETOT1
      ENDDO ! END OF LOOK OVER K-POINTS
      ETOT=-KBT*ETOT
PRINT*,'ETOT FROM DMFT_DETOT: ',ETOT
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_UTENSOR()
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
                                          CALL TRACE$PUSH('DMFT_UTENSOR')
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
!
!       ========================================================================
!       == CALCULATE LOCAL DIMENSIONS LMNX AND NCHI1 ===========================
!       ========================================================================
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
          CALL ERROR$STOP('DMFT_UTENSOR')
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
!
!       ========================================================================
!       == SCREEN U-TENSOR BY LOCAL HF-WEIGHT                                 ==
!       ========================================================================
        ATOMSET(IAT)%U=ATOMSET(IAT)%LHFWEIGHT*ATOMSET(IAT)%U
!
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
      SUBROUTINE DMFT_DYNAMICSOLVER(NLOC,NDIMD,NOMEGA,KBT &
     &                        ,GLOC,GLOCLAUR,UCHI,ETOT,SLOC,SLOCLAUR,DEDUCHI &
     &                        ,PIPHI,CHIPHI)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NLOC
      INTEGER(4),INTENT(IN)  :: NDIMD
      INTEGER(4),INTENT(IN)  :: NOMEGA
      REAL(8)   ,INTENT(IN)  :: KBT
      COMPLEX(8),INTENT(IN)  :: PIPHI(2*NLOC,2*NLOC)  !<PI|NATORB>
      COMPLEX(8),INTENT(IN)  :: CHIPHI(2*NLOC,2*NLOC) !<CHI|NATORB>
      COMPLEX(8),INTENT(IN)  :: GLOC(NLOC,NLOC,NDIMD,NOMEGA)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR(NLOC,NLOC,NDIMD,3)
      REAL(8)   ,INTENT(IN)  :: UCHI(NLOC,NLOC,NLOC,NLOC)
      REAL(8)   ,INTENT(OUT) :: ETOT
      COMPLEX(8),INTENT(OUT) :: SLOC(NLOC,NLOC,NDIMD,NOMEGA)
      COMPLEX(8),INTENT(OUT) :: SLOCLAUR(NLOC,NLOC,NDIMD,3)
      REAL(8)   ,INTENT(OUT) :: DEDUCHI(NLOC,NLOC,NLOC,NLOC)
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      COMPLEX(8)             :: U(2*NLOC,2*NLOC,2*NLOC,2*NLOC)
      COMPLEX(8)             :: DEDU(2*NLOC,2*NLOC,2*NLOC,2*NLOC)
      COMPLEX(8)             :: V(2*NLOC,2*NLOC,2*NLOC,2*NLOC)
      COMPLEX(8)             :: G(2*NLOC,2*NLOC,NOMEGA)
      COMPLEX(8)             :: GLAUR(2*NLOC,2*NLOC,3)
      COMPLEX(8)             :: S(2*NLOC,2*NLOC,NOMEGA)
      COMPLEX(8)             :: SLAUR(2*NLOC,2*NLOC,3)
      INTEGER(4)             :: NU,I,J
      INTEGER(4)             :: NFIL
!     **************************************************************************
                                     CALL TRACE$PUSH('DMFT$HFSOLVER')
!
!     ==========================================================================
!     == CONVERT U-TENSOR                                                     ==
!     ==========================================================================
      U(  :NLOC,  :NLOC,  :NLOC,  :NLOC)=UCHI                              
      U(NLOC+1:,NLOC+1:,NLOC+1:,NLOC+1:)=UCHI                              
      DO I=1,2*NLOC
        V=(0.D0,0.D0)
        DO J=1,2*NLOC
          V(I,:,:,:)=V(I,:,:,:)+CONJG(PIPHI(J,I))*U(J,:,:,:)
        ENDDO
      ENDDO
      DO I=1,2*NLOC
        U=(0.D0,0.D0)
        DO J=1,2*NLOC
          U(:,I,:,:)=U(:,I,:,:)+CONJG(PIPHI(J,I))*V(:,J,:,:)
        ENDDO
      ENDDO
      DO I=1,2*NLOC
        V=(0.D0,0.D0)
        DO J=1,2*NLOC
          V(:,:,I,:)=V(:,:,I,:)+U(:,:,J,:)*PIPHI(J,I)
        ENDDO
      ENDDO
      DO I=1,2*NLOC
        U=(0.D0,0.D0)
        DO J=1,2*NLOC
          U(:,:,:,I)=U(:,:,:,I)+V(:,:,:,J)*PIPHI(J,I)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == CONVERT GREENS FUNCTION                                              ==
!     ==========================================================================
      DO NU=1,NOMEGA
        CALL SPINOR$BLOWUP(NDIMD,NLOC,GLOC(:,:,:,NU),G(:,:,NU))
        G(:,:,NU)=MATMUL(TRANSPOSE(CONJG(CHIPHI)),MATMUL(G(:,:,NU),CHIPHI))
      ENDDO
      DO I=1,3
        CALL SPINOR$BLOWUP(NDIMD,NLOC,GLOCLAUR(:,:,:,I),GLAUR(:,:,I))
        GLAUR(:,:,I)=MATMUL(TRANSPOSE(CONJG(CHIPHI)) &
     &                      ,MATMUL(GLAUR(:,:,I),CHIPHI))
      ENDDO
!
!     ==========================================================================
!     == WRITE DATA
!     ==========================================================================
      CALL DMFT_SOLVERIO(2*NLOC,NOMEGA,KBT,G,GLAUR,U,ETOT,S,SLAUR,DEDU)
!
!     ==========================================================================
!     == CONVERT GREEN'S FUNCTION                                             ==
!     ==========================================================================
      DO NU=1,NOMEGA
        S(:,:,NU)=MATMUL(CHIPHI,MATMUL(S(:,:,NU),CONJG(TRANSPOSE(CHIPHI))))
        CALL SPINOR$SHRINKDOWN(NDIMD,NLOC,S(:,:,NU),SLOC(:,:,:,NU))
      ENDDO
      DO I=1,3
        SLAUR(:,:,I)=MATMUL(CHIPHI,MATMUL(SLAUR(:,:,I) &
     &                                   ,CONJG(TRANSPOSE(CHIPHI))))
        CALL SPINOR$SHRINKDOWN(NDIMD,NLOC,SLAUR(:,:,I),SLOCLAUR(:,:,:,I))
      ENDDO
!
!     ==========================================================================
!     == CONVERT DERIVATIVE OF U-TENSOR                                       ==
!     ==========================================================================
      U=DEDU
      DO I=1,2*NLOC
        V=(0.D0,0.D0)
        DO J=1,2*NLOC
          V(I,:,:,:)=V(I,:,:,:)+CONJG(PIPHI(I,J))*U(J,:,:,:)
        ENDDO
      ENDDO
      DO I=1,2*NLOC
        U=(0.D0,0.D0)
        DO J=1,2*NLOC
          U(:,I,:,:)=U(:,I,:,:)+CONJG(PIPHI(I,J))*V(:,J,:,:)
        ENDDO
      ENDDO
      DO I=1,2*NLOC
        V=(0.D0,0.D0)
        DO J=1,2*NLOC
          V(:,:,I,:)=V(:,:,I,:)+U(:,:,J,:)*PIPHI(I,J)
        ENDDO
      ENDDO
      DO I=1,2*NLOC
        U=(0.D0,0.D0)
        DO J=1,2*NLOC
          U(:,:,:,I)=U(:,:,:,I)+V(:,:,:,J)*PIPHI(I,J)
        ENDDO
      ENDDO
      DEDUCHI=REAL(U(  :NLOC,  :NLOC,  :NLOC,  :NLOC) &
     &            +U(NLOC+1:,NLOC+1:,NLOC+1:,NLOC+1:))
                                     CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVERIO(NORB,NOMEGA,KBT,G,GLAUR,U,ETOT,S,SLAUR,DEDU)
!     **************************************************************************
!     ** INTERFACE ROUTINE FOR THE COMMUNICATION WITH THE DMFT SOLVER         **
!     ** THE SPIN-UP ORBITALS ARE STORED FIRST AND THEN THE SPIN-DOWN         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB     !#(SPIN ORBITALS)
      INTEGER(4),INTENT(IN) :: NOMEGA   !#(POSITIVE MATSUBARA FREQUENCIES
      REAL(8)   ,INTENT(IN) :: KBT      
      COMPLEX(8),INTENT(IN) :: G(NORB,NORB,NOMEGA)    !GREENS FUNCTION
      COMPLEX(8),INTENT(IN) :: GLAUR(NORB,NORB,3)     !LAURENT EXPANSION OF G
      COMPLEX(8),INTENT(IN) :: U(NORB,NORB,NORB,NORB) !U-TENSOR
      REAL(8)   ,INTENT(OUT):: ETOT                   !LUTTINGER WARD FUNCTIONAL
      COMPLEX(8),INTENT(OUT):: S(NORB,NORB,NOMEGA)    !SELF ENERGY
      COMPLEX(8),INTENT(OUT):: SLAUR(NORB,NORB,3)     !LAURENT EXPANSION OF S
      COMPLEX(8),INTENT(OUT):: DEDU(NORB,NORB,NORB,NORB) 
      INTEGER(4)            :: NORB_,NOMEGA_
      INTEGER(4)            :: NFIL
!     **************************************************************************
CALL TESTG(NORB,NOMEGA,KBT,G,GLAUR)
      ETOT=0.D0
      S=(0.D0,0.D0)
      SLAUR=(0.D0,0.D0)
      DEDU=(0.D0,0.D0)
      RETURN
!
!     ==========================================================================
!     == WRITE DATA
!     ==========================================================================
!WARNING!! THE TWO FILES ARE NOT YET DEFINED
!CAREFUL: IN FORTRAN, THE FIRST INDEX RUNS FIRST
      CALL FILEHANDLER$UNIT('DMFTOSOLVER',NFIL)
      REWIND NFIL
      WRITE(NFIL,*)NORB,NOMEGA,KBT
      WRITE(NFIL,FMT='(2F20.10)')G
      WRITE(NFIL,FMT='(2F20.10)')GLAUR
      WRITE(NFIL,FMT='(2F20.10)')U
      CALL FILEHANDLER$CLOSE(NFIL)
!
!     ==========================================================================
!     == READ DATA DATA
!     ==========================================================================
      CALL FILEHANDLER$UNIT('SOLVERTODMFT',NFIL)
      REWIND NFIL
      READ(NFIL,*)NORB_,NOMEGA_,ETOT
      READ(NFIL,FMT='(2F20.10)')S
      READ(NFIL,FMT='(2F20.10)')SLAUR
      READ(NFIL,FMT='(2F20.10)')DEDU
      CALL FILEHANDLER$CLOSE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TESTG(NORB,NOMEGA,KBT,G,GLAUR)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB     !#(SPIN ORBITALS)
      INTEGER(4),INTENT(IN) :: NOMEGA   !#(POSITIVE MATSUBARA FREQUENCIES
      REAL(8)   ,INTENT(IN) :: KBT      
      COMPLEX(8),INTENT(IN) :: G(NORB,NORB,NOMEGA)    !GREENS FUNCTION
      COMPLEX(8),INTENT(IN) :: GLAUR(NORB,NORB,3)     !LAURENT EXPANSION OF G
      REAL(8)               :: PI
      INTEGER(4)            :: NU,I
      REAL(8)               :: OMEGA(NOMEGA)
      COMPLEX(8)            :: RHO(NORB,NORB)
      REAL(8)               :: FN(3)
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      RHO=(0.D0,0.D0)
      DO NU=1,NOMEGA
        OMEGA(NU)=REAL(2*NU-1,KIND=8)*PI*KBT
        RHO=RHO+KBT*G(:,:,NU)
      ENDDO
      RHO=RHO+TRANSPOSE(CONJG(RHO))
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,3,FN)
      RHO=RHO+GLAUR(:,:,1)*FN(1) &
     &       +GLAUR(:,:,2)*FN(2) &
     &       +GLAUR(:,:,3)*FN(3) 
!
      WRITE(*,FMT='(100("="),T30,"  ",A,"  ")')'DENSITY MATRIX IN TESTG'
      DO I=1,NORB
        WRITE(*,FMT='(100("(",2F10.5,")"))')RHO(I,:)
      ENDDO
      return
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_FOCK(NLOC,NDIMD,RHO,UCHI,ETOT,HAM)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NLOC
      INTEGER(4),INTENT(IN)  :: NDIMD
      COMPLEX(8),INTENT(OUT) :: RHO(NLOC,NLOC,NDIMD)
      REAL(8)   ,INTENT(OUT) :: UCHI(NLOC,NLOC,NLOC,NLOC)
      REAL(8)   ,INTENT(OUT) :: ETOT
      COMPLEX(8),INTENT(OUT) :: HAM(NLOC,NLOC,NDIMD)
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      INTEGER(4)             :: I,J,K,L
      REAL(8)                :: SVAR
!     **************************************************************************
                                     CALL TRACE$PUSH('DMFT_FOCK')
!
!     ==========================================================================
!     == CALCULATE HARTREE FOCK POTENTIAL                                     ==
!     ==========================================================================
      ETOT=0.D0
      HAM=(0.D0,0.D0)
      DO I=1,NLOC
        DO J=1,NLOC
          DO K=1,NLOC
            DO L=1,NLOC
              SVAR=-0.25D0*UCHI(I,J,K,L) 
              IF(SVAR.EQ.0.D0) CYCLE
              ETOT=ETOT+SVAR*REAL(SUM(RHO(K,J,:)*RHO(L,I,:)))
              HAM(K,J,:)=HAM(K,J,:)+SVAR*RHO(I,L,:)
              HAM(I,L,:)=HAM(I,L,:)+SVAR*RHO(K,J,:)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == PRINT IF DESIRED                                                     ==
!     ==========================================================================
      IF(TPRINT) THEN
        CALL SPINOR$PRINTMATRIX(6 &
     &                      ,'DENMAT (TXYZ) IN DMFT$HFSOLVER_WITHSET','TXYZ' &
     &                      ,1,NLOC,NDIMD,NLOC,RHO)
        CALL SPINOR$PRINTMATRIX(6 &
     &                      ,'H_XC (TXYZ) IN DMFT$HFSOLVER_WITHSET','TXYZ' &
     &                      ,1,NLOC,NDIMD,NLOC,HAM)
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
      REAL(8)   ,PARAMETER     :: TOL=1.D-5
      INTEGER(4),PARAMETER     :: NITER=500
      CHARACTER(16)            :: MIXTYPE
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT3(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: G(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR(NCHI,NCHI,NDIMD,3)
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
COMPLEX(8),ALLOCATABLE   :: X4SAVE(:,:,:)
COMPLEX(8),ALLOCATABLE   :: X4VEC(:)
REAL(8) :: AMP=1.D-3
COMPLEX(8)               :: DHAM1(NCHI,NCHI,NDIMD)
COMPLEX(8)               :: DHAM2(NCHI,NCHI,NDIMD)
COMPLEX(8)               :: DRHO1(NCHI,NCHI,NDIMD)
COMPLEX(8)               :: DRHO2(NCHI,NCHI,NDIMD)
INTEGER(4) :: ISPIN
REAL(8)    :: EIG(NCHI)
COMPLEX(8) :: U(NCHI,NCHI)
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

DO ISPIN=1,NSPIN
  DO IKPT=1,NKPTL
    CALL LIB$DIAGC8(NCHI,KSET(IKPT)%RHO(:,:,ISPIN),EIG(:NCHI),U)
    PRINT*,'EIGENVALUES OF RHO ',ISPIN,EIG(:NCHI)
  ENDDO
ENDDO


!!$DHAM1=(0.D0,0.D0)
!!$DHAM1(3,5,1)=AMP
!!$DHAM1(2,5,2)=AMP
!!$CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,DHAM1) ! CONVERT TO 0XYZ
!!$CALL DMFT_TESTRHO(TYPE,NCHI,NDIMD,DHAM1,DRHO1)
!!$CALL DMFT_TESTRHO(TYPE,NCHI,NDIMD,-DHAM1,DRHO2)
!!$DHAM2=(DRHO2-DRHO1)/(2.D0*AMP)
!!$CALL SPINOR_PRINTMATRIX(6,'HAM NUMERICAL ',1,NCHI,NDIMD,NCHI,DHAM2)
!       
!     ========================================================================
!     ==  DEVIATION FROM TARGET DENSITY MATRIX                              ==
!     ========================================================================
      DO IKPT=1,NKPTL
!       == LAURENT EXPANSION FOR THE GREENS FUNCTION =========================
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
!       ========================================================================
!       ==  ADJUST H0 TO OBEY DENSITY MATRIX CONSTRAINT                       ==
!       ========================================================================
        DO ITER=1,NITER
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
!!$CALL SPINOR_PRINTMATRIX(6,'G DIRECT ',1,NCHI,NDIMD,NCHI,G)
!!$CALL SPINOR_PRINTMATRIX(6,'G TEST ',1,NCHI,NDIMD,NCHI &
!!$&                   ,KSET(IKPT)%SINV/(CI*OMEGA(NU)))
!!$PRINT*,'KBT ',KBT
!!$STOP 'FORCED'
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
!RHO0=0.5D0*KET(IKPT)%SINV

!!$CALL SPINOR_PRINTMATRIX(6,'DEVRHO DIRECT ',1,NCHI,NDIMD,NCHI,DEVRHO)
!!$CALL SPINOR_PRINTMATRIX(6,'DEVRHO TEST ',1,NCHI,NDIMD,NCHI &
!!$&                   ,KSET(IKPT)%SINV/2.D0)
!!$PRINT*,'KBT ',KBT
!!$STOP 'FORCED'
!
!         == SUBTRACT TARGET DENSITY ===========================================
!!$CALL SPINOR_PRINTMATRIX(6,'RHO FROM G ',1,NCHI,NDIMD,NCHI,DEVRHO)
!!$CALL SPINOR_PRINTMATRIX(6,'KSET%RHO   ',1,NCHI,NDIMD,NCHI,KSET(IKPT)%RHO)
          DEVRHO=DEVRHO-KSET(IKPT)%RHO
!
!         == REGULARIZE X4 =====================================================
          GLAUR(:,:,:,1)=GLAUR1
          GLAUR(:,:,:,2)=GLAUR2
          GLAUR(:,:,:,3)=GLAUR3
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,GLAUR(:,:,:,1)) ! CNVRTTO UPDOWN
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,GLAUR(:,:,:,2)) ! CNVRTTO UPDOWN
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,GLAUR(:,:,:,3)) ! CNVRTTO UPDOWN
          CALL DMFT_ADDX4LAUR(NDIMD,NCHI,KBT,NOMEGA,OMEGA,GLAUR,LX4,NSPIN,X4)
IF(.NOT.ALLOCATED(X4SAVE)) THEN
ALLOCATE(X4SAVE(LX4,LX4,NSPIN))
  X4SAVE=X4
!!$  DO ISPIN=1,NSPIN
!!$    DO I1=1,NCHI**2
!!$      DO I2=1,NCHI**2
!!$        IF(ABS(X4(I1,I2,ISPIN)).GT.1.D-3) THEN
!!$          PRINT*,'X4 ',I1,I2,ISPIN,X4(I1,I2,ISPIN)
!!$        END IF
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO
ELSE
  X4=X4SAVE
END IF

!
DO ISPIN=1,NSPIN
    CALL LIB$DIAGC8(NCHI,KSET(IKPT)%RHO(:,:,ISPIN)+DEVRHO(:,:,ISPIN),EIG(:NCHI),U)
    PRINT*,'EIGENVALUES OF RHO ',ISPIN,EIG(:NCHI)
IF(MINVAL(EIG).LT.0.D0.OR.MAXVAL(EIG).GT.1.D0) THEN
   CALL ERROR$MSG('OCCUPATIONS OUT OF RANGE')
   CALL ERROR$R8VAL('MAX(F)',MAXVAL(EIG))
   CALL ERROR$R8VAL('MIN(F)',MINVAL(EIG))
   CALL ERROR$STOP('DMFT_CONSTRAINTS_WITHKSET')
END IF
ENDDO
!
CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,DEVRHO) ! CONVERT TO UPDOWN
MAXDEV=MAXVAL(ABS(DEVRHO))
DO IDIMD=1,NDIMD
    PRINT*,'ITER=',ITER,' IKPT=',IKPT,' IDIMD=',IDIMD &
  &       ,' MAXVAL(ABS(X4))  ',MAXVAL(ABS(X4(:,:,IDIMD))) &
  &       ,' MAX(DEVRHO) ',MAXVAL(ABS(DEVRHO(:,:,IDIMD)))
ENDDO
CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,DEVRHO) ! CONVERT TO TXYZ
!
!          MIXTYPE='MIX'
!          MIXTYPE='LINEAR'
          MIXTYPE='LINEAR0'
          IF(MIXTYPE.EQ.'MIX') THEN
            DH0=AMIX*DEVRHO
          ELSE IF(MIXTYPE.EQ.'LINEAR') THEN
!
!!$ALLOCATE(X4VEC(LX4))
!!$CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,DHAM1) ! CONVERT TO UPDOWN
!!$DO ISPIN=1,NSPIN
!!$  DO I1=1,NCHI
!!$    DO I2=1,NCHI
!!$      X4VEC(I1+NCHI*(I2-1))=DHAM1(I1,I2,ISPIN)
!!$    ENDDO
!!$  ENDDO
!!$  X4VEC=MATMUL(X4(:,:,ISPIN),X4VEC)
!!$  DO I1=1,NCHI
!!$    DO I2=1,NCHI
!!$      DRHO2(I1,I2,ISPIN)=X4VEC(I1+NCHI*(I2-1))
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
!!$DRHO2=DRHO2/AMP
!!$CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,DRHO2) ! CONVERT TO 0XYZ
!!$CALL SPINOR_PRINTMATRIX(6,'DRHO ANALYTICAL ',1,NCHI,NDIMD,NCHI,DRHO2)

            CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,DEVRHO) ! CONVERT TO UPDOWN
!           == DEVRHO+X4*DH0=0 -> DH0
            CALL DMFT_SOLVEX4(NDIMD,NCHI,DEVRHO,LX4,NSPIN,X4,DH0)
            CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,DH0)    ! CONVERT TO TXYZ
!            DH0=0.5D0*DH0   !PART OF THE CONVERSION FOR HAMILTONIANS
          ELSE IF(MIXTYPE.EQ.'LINEAR0') THEN
!!$PRINT*,'X4SAVE ',X4SAVE(1,1,1)
!!$CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,KSET(IKPT)%SINV) ! CONVERT TO UPDOWN
!!$PRINT*,'X4NEW ',-1.D0/(4.D0*KBT)*KSET(IKPT)%SINV(1,1,1)**2
!!$STOP
            CALL SPINOR$MATMUL(NDIMD,NCHI,DEVRHO,KSET(IKPT)%SMAT,MAT)
            CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SMAT,MAT,DH0)
            DH0=4.D0*KBT*DH0
   DH0=0.5D0*DH0
!!$CALL SPINOR_PRINTMATRIX(6,'DH0 A',1,NCHI,NDIMD,NCHI,DH0)
!!$CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,DEVRHO) ! CONVERT TO UPDOWN
!!$CALL DMFT_SOLVEX4(NDIMD,NCHI,DEVRHO,LX4,NSPIN,X4,DH0)
!!$CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,DH0)    ! CONVERT TO TXYZ
!!$CALL SPINOR_PRINTMATRIX(6,'DH0 B',1,NCHI,NDIMD,NCHI,DH0)
          ELSE
            CALL ERROR$MSG('MIXTYPE NOT RECOGNIZED')
            CALL ERROR$STOP('DMFT_CONSTRAINTS_WITHKSET')
          END IF
PRINT*,'MAXVAL OF DH ',MAXVAL(ABS(DH0)),MAXLOC(ABS(DH0))
!
!         ======================================================================
!         == ADD DH                                                           ==
!         ======================================================================
          KSET(IKPT)%H0=KSET(IKPT)%H0+DH0
!
!         ======================================================================
!         == CHECK CONVERGENCE                                                ==
!         ======================================================================
          MAXDEV=MAXVAL(ABS(DEVRHO))
          PRINT*,'MAXDEV',ITER,IKPT,MAXDEV
          CONVG=MAXDEV.LT.TOL
          IF(CONVG) EXIT
        ENDDO ! END OF LOOP OVER ITERATIONS
        IF(.NOT.CONVG) THEN
          CALL ERROR$MSG('LOOP NOT CONVERGED')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('ITER',ITER)
          CALL ERROR$R8VAL('MAX. DEVIATION',MAXDEV)
          CALL ERROR$R8VAL('TOLERANCE',TOL)
          CALL ERROR$STOP('DMFT_CONSTRAINTS_WITHKSET')
        END IF
      ENDDO   !END OF LOOP OVER K-POINTS
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
      SUBROUTINE DMFT_CONSTRAINTS(TYPE)
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
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NDIMD,NAT,NOMEGA &
     &                      ,OMEGA,KBT,MU,AMIX,KSET,ATOMSET
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: TYPE  ! CAN BE 'HRHO','H0'
      REAL(8)   ,PARAMETER     :: TOL=1.D-5
      INTEGER(4),PARAMETER     :: NITER=100000
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
      COMPLEX(8)               :: RHO(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: DH0(NCHI,NCHI,NDIMD)
      LOGICAL(4)               :: CONVG
      REAL(8)                  :: MAXDEV
      INTEGER(4)               :: IKPT,ITER,IDIMD,NU,IAT
      INTEGER(4)               :: I1,I2
      LOGICAL(4)               :: TH0
      REAL(8)                  :: FN(2)
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_CONSTRAINTS')
      TH0=(TYPE.EQ.'H0') 
      IF(.NOT.(TH0.OR.TYPE.EQ.'HRHO')) THEN
        CALL ERROR$MSG('ILLEGAL VALUE OF TYPE. (MAY BE "H0" OR "HRHO")')
        CALL ERROR$CHVAL('TYPE',TYPE)
        CALL ERROR$STOP('DMFT_CONSTRAINTS')     
      END IF
      IF(.NOT.TH0)  THEN
        DO IKPT=1,NKPTL
          KSET(IKPT)%H0=KSET(IKPT)%HRHO
        ENDDO
      END IF
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,2,FN)
!       
!     ========================================================================
!     ==  DEVIATION FROM TARGET DENSITY MATRIX                              ==
!     ========================================================================
      DO IKPT=1,NKPTL
!       == LAURENT EXPANSION FOR THE GREENS FUNCTION =========================
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
!       ========================================================================
!       ==  ADJUST H0 TO OBEY DENSITY MATRIX CONSTRAINT                       ==
!       ========================================================================
        DO ITER=1,NITER
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
!         == LOOP OVER OMEGA ===================================================
          RHO=(0.D0,0.D0)
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
            RHO=RHO+KBT*G
          ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!         == INCLUDE NEGATIVE FREQUENCIES ======================================
          DO IDIMD=1,NDIMD
            RHO(:,:,IDIMD)=RHO(:,:,IDIMD)+CONJG(TRANSPOSE(RHO(:,:,IDIMD)))
          ENDDO
!         == ADD TAILS (GLAUR3 DOES NOT CONTRIBUTE) ============================
          RHO=RHO+FN(1)*GLAUR1+FN(2)*GLAUR2
!RHO0=0.5D0*KET(IKPT)%SINV
          DEVRHO=RHO-KSET(IKPT)%RHO
!
!         ======================================================================
!         == ADJUST HAMILTONIAN (LAGRANGE MULTIPLIER)                         ==
!         ======================================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,DEVRHO,KSET(IKPT)%SMAT,MAT)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SMAT,MAT,DH0)
          DH0=4.D0*KBT*DH0
          KSET(IKPT)%H0=KSET(IKPT)%H0+DH0
!
!         ======================================================================
!         == CHECK CONVERGENCE                                                ==
!         ======================================================================
IF(MOD(ITER,100).EQ.0)PRINT*,'MAXVAL OF DH /DRHO ',ITER,MAXVAL(ABS(DH0)),MAXVAL(ABS(DEVRHO))
          MAXDEV=MAXVAL(ABS(DEVRHO))
          CONVG=MAXDEV.LT.TOL
          IF(CONVG) EXIT
        ENDDO ! END OF LOOP OVER ITERATIONS
!!$CALL SPINOR_PRINTMATRIX(6,'FINAL KSET%H0 ',1,NCHI,NDIMD,NCHI,KSET(IKPT)%H0)
!!$CALL SPINOR_PRINTMATRIX(6,'FINAL DEVRHO ',1,NCHI,NDIMD,NCHI,DEVRHO)
        IF(.NOT.CONVG) THEN
          CALL ERROR$MSG('LOOP NOT CONVERGED')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('ITER',ITER)
          CALL ERROR$R8VAL('MAX. DEVIATION',MAXDEV)
          CALL ERROR$R8VAL('TOLERANCE',TOL)
          CALL ERROR$STOP('DMFT_CONSTRAINTS')
        END IF
      ENDDO   !END OF LOOP OVER K-POINTS
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
      SUBROUTINE DMFT_TESTRHO(TYPE,NCHI_,NDIMD_,DHAM,DRHO)
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
     &                      ,OMEGA,KBT,MU,KSET,ATOMSET
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: TYPE  ! CAN BE 'HRHO','H0'
      INTEGER(4),INTENT(IN)    :: NCHI_
      INTEGER(4),INTENT(IN)    :: NDIMD_
      COMPLEX(8),INTENT(IN)    :: DHAM(NCHI_,NCHI_,NDIMD_)
      COMPLEX(8),INTENT(OUT)   :: DRHO(NCHI_,NCHI_,NDIMD_)
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT3(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: G(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR1(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR3(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: SLAUR1(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: SLAUR2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: SLAUR3(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: CSVAR
      INTEGER(4)               :: IKPT,IDIMD,NU,IAT
      INTEGER(4)               :: I1,I2
      LOGICAL(4)               :: TH0
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_TESTRHO')
      TH0=(TYPE.EQ.'H0') 
      IF(.NOT.(TH0.OR.TYPE.EQ.'HRHO')) THEN
        CALL ERROR$MSG('ILLEGAL VALUE OF TYPE. (MAY BE "H0" OR "HRHO")')
        CALL ERROR$CHVAL('TYPE',TYPE)
        CALL ERROR$STOP('DMFT_TESTRHO')     
      END IF
      IF(.NOT.TH0)  THEN
        DO IKPT=1,NKPTL
          KSET(IKPT)%H0=(0.D0,0.D0)
        ENDDO
      END IF
!       
!     ========================================================================
!     ==  DEVIATION FROM TARGET DENSITY MATRIX                              ==
!     ========================================================================
      DO IKPT=1,NKPTL
!       == LAURENT EXPANSION FOR THE GREENS FUNCTION =========================
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
        MAT=-MU*KSET(IKPT)%SMAT
        MAT=MAT+KSET(IKPT)%H0+DHAM
        IF(TH0) MAT=MAT+SLAUR1(:,:,:)
!       == GLAUR2=SINV*MAT*SINV ==============================================
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT,KSET(IKPT)%SINV,MAT2)
        CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,MAT2,GLAUR2)
!
        MAT=-MU*KSET(IKPT)%SMAT
        MAT=MAT+KSET(IKPT)%H0+DHAM
        IF(TH0) MAT=MAT+SLAUR1(:,:,:)
!       == MAT3=MAT*(SINV*MAT) ===============================================
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT,KSET(IKPT)%SINV,MAT2)
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT,MAT2,MAT3)
        IF(TH0) MAT3=MAT3+SLAUR2(:,:,:)
!       == GLAUR3=SINV*(MAT3*SINV) ===========================================
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT3,KSET(IKPT)%SINV,MAT)
        CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,MAT,GLAUR3)
!
        DRHO=(0.D0,0.D0)
        DO NU=1,NOMEGA
!         == CONSTRUCT LATTICE GREENS FUNCTION =============================
          MAT=(CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT
          MAT=MAT-KSET(IKPT)%H0-DHAM
          IF(TH0) THEN
            DO IAT=1,NAT
              I1=ATOMSET(IAT)%ICHI1
              I2=ATOMSET(IAT)%ICHI2
              MAT(I1:I2,I1:I2,:)=MAT(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
            ENDDO
          END IF
          CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
!         == SUBTRACT LAURENT EXPANSION TO IMPROVE OMEGA-CONVERGENCE =========
          CSVAR=1.D0/(CI*OMEGA(NU))
          DRHO=DRHO+KBT*(G-CSVAR*(GLAUR1+CSVAR*(GLAUR2+CSVAR*GLAUR3)))
        ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!       == INCLUDE NEGATIVE FREQUENCIES ======================================
!       == ALREADY INCLUDED FOR X4
        DO IDIMD=1,NDIMD
          DRHO(:,:,IDIMD)=DRHO(:,:,IDIMD) &
     &                   +CONJG(TRANSPOSE(DRHO(:,:,IDIMD)))
        ENDDO
!       == ADD TAILS (GLAUR3 DOES NOT CONTRIBUTE) ============================
        DRHO=DRHO+0.5D0*GLAUR1-0.25D0/KBT*GLAUR2
!
!       == SUBTRACT TARGET DENSITY ===========================================
        DRHO=DRHO-KSET(IKPT)%RHO
      ENDDO   !END OF LOOP OVER K-POINTS
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
      COMPLEX(8),INTENT(INOUT)    :: X4(LX4,LX4,NSPIN)
      COMPLEX(8),INTENT(OUT)   :: DH(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: BVEC(LX4)
      COMPLEX(8)               :: XVEC(LX4)
      INTEGER(4)               :: IA,IB
      INTEGER(4)               :: IDA,IDB
      INTEGER(4)               :: ICHIA,ICHIB
      INTEGER(4)               :: IAB
      INTEGER(4)               :: IDIMAB
      INTEGER(4)               :: ISPIN 
REAL(8) :: EIG(LX4)
COMPLEX(8) :: U(LX4,LX4)
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
! =============================================
!!$DO IA=1,LX4
!!$  X4(IA,IA,ISPIN)=X4(IA,IA,ISPIN)-1.D-2
!!$ENDDO

!         ====
PRINT*,'LX4 ',LX4,ISPIN
CALL SPINOR_PRINTMATRIX(6,'DRHO (UPDOWN)',1,NCHI,1,NCHI,DRHO(:,:,ISPIN))
CALL LIB$DIAGC8(NCHI,DRHO(:,:,ISPIN),EIG(:NCHI),U)
PRINT*,'EIGENVALUES OF DRHO ',ISPIN,EIG(:NCHI)
CALL LIB$DIAGC8(LX4,X4(:,:,ISPIN),EIG,U)
PRINT*,'EIGENVALUES OF X4 ',ISPIN,EIG
IF(MAXVAL(EIG).GT.0.D0) PRINT*,'WARNING!!!',MAXVAL(EIG)
!
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
          DO IA=1,NCHI
            DO IB=1,NCHI
              IAB=IA+NCHI*(IB-1)
              DO IC=1,NCHI
                DO ID=1,NCHI
                  ICD=IC+NCHI*(ID-1)
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
            IAB=ICHIA+NCHI*(IDA-1+2*(ICHIB+NCHI*(IDB-1)))
            DO IDC=1,2
            DO ICHIC=1,NCHI
              DO IDD=1,2
              DO ICHID=1,NCHI
                ICD=ICHIC+NCHI*(IDC-1+2*(ICHID+NCHI*(IDD-1)))
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
      SUBROUTINE DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,N,FN)
!     **************************************************************************
!     ** CALCULATES THE DIFFERENCE OF THE INFINITE AND THE TRUNCATED          **
!     ** MATSUBARA SUM OF THE LAURENT TERMS 1/(I*OMEGA_NU)**J                 **
!     **                                                                      **
!     ** THIS METHOD SUFFERS FROM THE FACT THAT THE LARGEST CONTRIBUTION      **
!     ** COMES FROM THE POINTS NEXT TO THE ORIGIN. THUS ONE CALCULATES THE    **
!     ** TIN DIFFERENCE OF TWO VERY LARGE NUMBERS. THEREFORE, IT ONLY MAKES   **
!     ** SENCE TO GO TO N=6 AT MOST. ALL OTHERS ARE PROBABLY WORSE THAN       **
!     ** LEAVING THE CORRECTION OUT.                                          **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)    :: KBT
      INTEGER(4),INTENT(IN)    :: NOMEGA
      REAL(8)   ,INTENT(IN)    :: OMEGA(NOMEGA)
      INTEGER(4),INTENT(IN)    :: N
      REAL(8)   ,INTENT(OUT)   :: FN(N)
!     ==  PARTIAL_X^(J-1):1/(1+EXP(X))
      REAL(8),PARAMETER :: ARR(20) &
     &   =(/0.5D0,-0.25D0,0.D0,0.125D0,0.D0,-0.25D0,0.D0,1.0625D0,0.D0,-7.75D0 &
     &     ,0.D0,86.375D0,0.D0,1365.25D0,0.D0,29049.03125D0,0.D0,800572.75D0 &
     &     ,0.D0,27741322.625D0/)
      REAL(8)                  :: FACTOR(20)
      REAL(8)                  :: SVAR
      INTEGER(4)               :: J,NU
!     **************************************************************************
      IF(N.GT.6) THEN
        CALL ERROR$MSG('N EXCEEDS INTERNAL MAXIMUM')
        CALL ERROR$STOP('DMFT_REGMATSUBARA')
      END IF
!
!     ==========================================================================
!     == ATTACH 1/[KBT**(J-1) * (J-1)!]
!     ==========================================================================
!     == FACTOR=1/[(J-1)!]*PARTIAL^(J-1) [1/(1+EXP(X))]
      SVAR=1.D0
      DO J=1,N
        FACTOR(J)=ARR(J)*SVAR            ! SVAR=1/[KBT**(J-1) * (J-1)!]
        SVAR=SVAR/(KBT*REAL(J,KIND=8))   ! SVAR=1/[KBT**J * J!]
      ENDDO
!
!     ==========================================================================
!     == SUBTRACT TRUNCATED MATSUBARA SUM                                     ==
!     ==========================================================================
      FN(:)=FACTOR(:N)
      DO NU=1,NOMEGA              ! FN(J)=SUM_NU: 1/(I*OMEGA_NU)**JJ
        SVAR=-1.D0/OMEGA(NU)**2   ! [1.D0/(CI*OMEGA(NU))]**2
        DO J=1,N/2                ! REAL PART FROM -OMEGA
          FN(2*J)=FN(2*J)-2.D0*KBT*SVAR**J ! ODD POWERS VANISH BECAUSE OF -OMEGA
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_ADDX4LAUR(NDIMD,NCHI,KBT,NOMEGA,OMEGA,GLAUR,LX4,NSPIN,X4)
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NDIMD
      INTEGER(4),INTENT(IN)    :: NCHI
      INTEGER(4),INTENT(IN)    :: NOMEGA
      REAL(8)   ,INTENT(IN)    :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(IN)    :: KBT
      INTEGER(4),INTENT(IN)    :: LX4
      INTEGER(4),INTENT(IN)    :: NSPIN
      COMPLEX(8),INTENT(IN)    :: GLAUR(NCHI,NCHI,NDIMD,3)
      COMPLEX(8),INTENT(INOUT) :: X4(LX4,LX4,NSPIN)
      REAL(8)                  :: FN(6),F2,F4,F6
      INTEGER(4)               :: IA,IB,IC,ID
      INTEGER(4)               :: IDA,IDB,IDC,IDD
      INTEGER(4)               :: ICHIA,ICHIB,ICHIC,ICHID
      INTEGER(4)               :: IAB,ICD
      INTEGER(4)               :: IDIMAC,IDIMDB,IDIMCA,IDIMBD
      INTEGER(4)               :: ISPIN 
!     **************************************************************************
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,6,FN)
      F2=FN(2)
      F4=FN(4)
      F6=FN(6)
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
          DO IA=1,NCHI
            DO IB=1,NCHI
              IAB=IA+NCHI*(IB-1)
              DO IC=1,NCHI
                DO ID=1,NCHI
                  ICD=IC+NCHI*(ID-1)
!                 == THE SECOND TERM IS FOR -I*OMEGA_\NU =======================
!                 == G(-I*OMEGA)=G^\DAGGER(I*OMEGA) ============================
                  X4(IAB,ICD,ISPIN)=X4(IAB,ICD,ISPIN) &
     &                       +F2* GLAUR(IA,IC,ISPIN,1)*GLAUR(ID,IB,ISPIN,1) &
     &                       +F4*(GLAUR(IA,IC,ISPIN,1)*GLAUR(ID,IB,ISPIN,3) &
     &                           +GLAUR(IA,IC,ISPIN,2)*GLAUR(ID,IB,ISPIN,2) &
     &                           +GLAUR(IA,IC,ISPIN,3)*GLAUR(ID,IB,ISPIN,1)) &
     &                       +F6* GLAUR(IA,IC,ISPIN,3)*GLAUR(ID,IB,ISPIN,3) 
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
            IAB=ICHIA+NCHI*(IDA-1+2*(ICHIB+NCHI*(IDB-1)))
            DO IDC=1,2
            DO ICHIC=1,NCHI
              DO IDD=1,2
              DO ICHID=1,NCHI
                ICD=ICHIC+NCHI*(IDC-1+2*(ICHID+NCHI*(IDD-1)))
!
                IDIMAC=IDA+2*(IDC-1)
                IDIMCA=IDC+2*(IDA-1)
                IDIMDB=IDD+2*(IDB-1)
                IDIMBD=IDB+2*(IDD-1)
                X4(IAB,ICD,1)=X4(IAB,ICD,1) &
     &           +F2* GLAUR(ICHIA,ICHIC,IDIMAC,1)*GLAUR(ICHID,ICHIB,IDIMDB,1) &
     &           +F4*(GLAUR(ICHIA,ICHIC,IDIMAC,1)*GLAUR(ICHID,ICHIB,IDIMDB,3) &
     &               +GLAUR(ICHIA,ICHIC,IDIMAC,2)*GLAUR(ICHID,ICHIB,IDIMDB,2) &
     &               +GLAUR(ICHIA,ICHIC,IDIMAC,3)*GLAUR(ICHID,ICHIB,IDIMDB,1)) &
     &           +F6* GLAUR(ICHIA,ICHIC,IDIMAC,3)*GLAUR(ICHID,ICHIB,IDIMDB,3) 
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
        WKPTL=KSET(IKPT)%WKPT
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
      SUBROUTINE DMFT_ADDTOHPSI()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY  : TON,NKPTL,NSPIN,NDIM,NDIMD,NB,NCHI,NAT &
     &                        ,IPROOFCHI,ATOMSET,KSET
      USE WAVES_MODULE, ONLY : GSET,WAVES_SELECTWV,THIS,MAP
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      COMPLEX(8),ALLOCATABLE :: DHPIPSI(:,:,:) !(NDIM,NCHI,NB)
      COMPLEX(8),ALLOCATABLE :: MAT(:,:,:)     !(NCHI,NCHI,NDIMD)
      LOGICAL(4)             :: TCHK,TRESET
      INTEGER(4)             :: NBH
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: NFIL
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IPRO,IDIM1,IDIM2,IDIMD,IAT
      INTEGER(4)             :: I1,I2
!     **************************************************************************
      IF(.NOT.TON) RETURN
PRINT*,'ENTERING DMFT_ADDTOHPSI'
                              CALL TRACE$PUSH('DMFT_ADDTOHPSI')
!
!     ==========================================================================
!     == CHECK IF HTBC ALREADY CONTAINS INFORMATION                           ==
!     ==========================================================================
      CALL LMTO$GETL4('ON',TCHK)
PRINT*,'LMTO ON? ',TCHK
      IF(TCHK) CALL LMTO$GETL4('THTBC',TCHK)
      TRESET=.NOT.TCHK
      CALL LMTO$SETL4('THTBC',.TRUE.)
PRINT*,'RESET? ',TRESET
!
!     ==========================================================================
!     ==  PRINT LOCAL HAMILTONIAN (STATIC AND DOUBLE COUNTING)                ==
!     ==========================================================================
      IF(TPRINT) THEN
        NFIL=6 
        WRITE(NFIL,FMT='(82("-"),T10," ",A," ")')'LOCAL HAMILTONIAN (FOCK+DC)' 
        WRITE(NFIL,FMT='(82("-"),T10," ",A," ")')'FROM DMFT_ADDTOHPSI'
        WRITE(NFIL,*)'LOCAL HAMILTONIAN (FOCK+DC)' 
        DO IAT=1,NAT
          IF(ATOMSET(IAT)%NLOC.LE.0) CYCLE ! NO DOUBLE COUNTING          
          LMNX=ATOMSET(IAT)%DENMAT%LMNX
          PRINT*,'ATOM ',IAT,' LMNX=',LMNX
          CALL SPINOR_PRINTMATRIX(NFIL,'DC HAMILTONIAN(TXYZ)',1,LMNX &
    &                            ,NDIMD,LMNX,ATOMSET(IAT)%DENMAT%H)
        ENDDO
      END IF
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
        IF(ATOMSET(IAT)%NLOC.LE.0) CYCLE ! NO DOUBLE COUNTING          
        LMNX=ATOMSET(IAT)%DENMAT%LMNX
!CALL SPINOR_PRINTMATRIX(6,'ATOMSET%H',1,LMNX,NDIMD,LMNX,ATOMSET(IAT)%DENMAT%H)
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
            CALL ERROR$STOP('DMFT_ADDTOHPSI')
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
!!$DO IAT=1,NAT
!!$  IF(ATOMSET(IAT)%NLOC.LE.0) CYCLE ! NO DOUBLE COUNTING          
!!$  I1=ATOMSET(IAT)%ICHI1
!!$  I2=ATOMSET(IAT)%ICHI2
!!$  CALL SPINOR_PRINTMATRIX(6,'CORR. HAMILTONIAN',1,I2-I1+1,NDIMD,I2-I1+1,MAT(I1:I2,I1:I2,:))
!!$ENDDO
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
          IPRO=0
          DO IAT=1,NAT
            LMNX=ATOMSET(IAT)%DENMAT%LMNX
            I1=IPRO+1
            I2=IPRO+LMNX
            IPRO=IPRO+LMNX
            IF(ATOMSET(IAT)%NLOC.LE.0) CYCLE ! NO DOUBLE COUNTING          
!
            DO IDIM2=1,NDIM
              DO IDIM1=1,NDIM
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
      IMPLICIT NONE
      INTEGER(4)           :: L=2
      INTEGER(4),PARAMETER :: NKAPPA=5
      INTEGER(4)           :: I
      REAL(8)              :: RAD
      REAL(8)              :: K2(NKAPPA)=(/-0.2D0,-0.1D0,0.D0,1.D0,2.D0/)
      REAL(8)              :: JVAL(NKAPPA),JDER(NKAPPA)
      REAL(8)              :: KVAL(NKAPPA),KDER(NKAPPA)
      INTEGER(4)           :: NFIL1=12,NFIL2=14,IKAPPA
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
      SUBROUTINE SPINOR$BLOWUP(NDIMD,NCHI,A,B)
!     **************************************************************************
!     ** CONVERT A MATRIX FROM THE 0XYZ REPRESENTATION INTO THE               **
!     ** (UU,UD/DU,DD) REPRESENTATION                                         **
!     **************************************PETER BLOECHL GOSLAR 2013***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      COMPLEX(8),INTENT(IN) :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(OUT):: B(2*NCHI,2*NCHI)
      COMPLEX(8)            :: A1(NCHI,NCHI,NDIMD)
!     **************************************************************************
      A1=A
      CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,A1)
      B=(0.D0,0.D0)
      IF(NDIMD.EQ.1) THEN
        B(  :NCHI,  :NCHI)=A1(:,:,1)
        B(NCHI+1:,NCHI+1:)=A1(:,:,1)
      ELSE IF(NDIMD.EQ.2) THEN
        B(  :NCHI,  :NCHI)=A1(:,:,1)
        B(NCHI+1:,NCHI+1:)=A1(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        B(  :NCHI,  :NCHI)=A1(:,:,1)
        B(NCHI+1:,  :NCHI)=A1(:,:,2)
        B(  :NCHI,NCHI+1:)=A1(:,:,3)
        B(NCHI+1:,NCHI+1:)=A1(:,:,4)
      ELSE
        CALL ERROR$MSG('ILLEGAL VALUE FOR NDIMD. (MUST BE 1,2 OR 4)')
        CALL ERROR$I4VAL('NDIMD',NDIMD)
        CALL ERROR$STOP('SPINOR$BLOWUP')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINOR$SHRINKDOWN(NDIMD,NCHI,B,A)
!     **************************************************************************
!     ** CONVERT MATRIX B FROM THE (UU,UD/DU,DD) REPRESENTATION               **
!     ** INTO MATRIX A IN 0XYZ REPRESENTATION                                 **
!     **************************************PETER BLOECHL GOSLAR 2013***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      COMPLEX(8),INTENT(IN) :: B(2*NCHI,2*NCHI)
      COMPLEX(8),INTENT(OUT):: A(NCHI,NCHI,NDIMD)
!     **************************************************************************
      IF(NDIMD.EQ.1) THEN
        A(:,:,1)=0.5D0*(B(  :NCHI,  :NCHI)+B(NCHI+1:,NCHI+1:))
      ELSE IF(NDIMD.EQ.2) THEN
        A(:,:,1)=B(  :NCHI,  :NCHI)
        A(:,:,2)=B(NCHI+1:,NCHI+1:)
      ELSE IF(NDIMD.EQ.4) THEN
        A(:,:,1)=B(  :NCHI,  :NCHI)
        A(:,:,2)=B(NCHI+1:,  :NCHI)
        A(:,:,3)=B(  :NCHI,NCHI+1:)
        A(:,:,4)=B(NCHI+1:,NCHI+1:)
      ELSE
        CALL ERROR$MSG('ILLEGAL VALUE FOR NDIMD. (MUST BE 1,2 OR 4)')
        CALL ERROR$I4VAL('NDIMD',NDIMD)
        CALL ERROR$STOP('SPINOR$BLOWUP')
      END IF
      CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,A)
      RETURN
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
