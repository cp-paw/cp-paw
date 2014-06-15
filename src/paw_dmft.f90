!........1.........2.........3.........4.........5.........6.........7.........8
MODULE DMFT_MODULE
!*******************************************************************************
!**                                                                           **
!**  DESCRIPTION OF VARIABLES                                                 **
!**    GLOCLAUR1,1,2 EXPANSION OF SPIN-AVERAGED GLOC IN 1/(CI*HBAR*OMEGA)     **
!**                                                                           **
!*******************************************************************************
TYPE DENMAT_TYPE   ! CONSIDERS ALL ORBITALS ON THIS SITE
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
  COMPLEX(8),POINTER   :: GLOCLAUR(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,NLAU+1)
  COMPLEX(8),POINTER   :: DPHIDG(:,:,:,:)   => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: DPHIDGLAUR(:,:,:,:)=> NULL() !(NLOC,NLOC,NDIMD,NLAU)
  COMPLEX(8),POINTER   :: SLOC(:,:,:,:)     => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: SLOCLAUR(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,NLAU)
  COMPLEX(8),POINTER   :: RHO(:,:,:)        => NULL() !LOCAL DENSITY MATRIX
  COMPLEX(8),POINTER   :: SMAT(:,:,:)       => NULL() !OVERLAP OF LOCAL ORBITALS
  TYPE(DENMAT_TYPE)    :: DENMAT
  TYPE(NATORB_TYPE)    :: NATORB
END TYPE ATOMSET_TYPE
TYPE KSET_TYPE
  REAL(8)            :: WKPT
  LOGICAL(4)         :: TADDMINUSK    !ADD  THE TERM FOR -K
  REAL(8)   ,POINTER :: F(:,:)        ! (NB,NSPIN)         
  COMPLEX(8),POINTER :: PIPSI(:,:,:,:)!(NDIM,NCHI,NB,NSPIN)
  COMPLEX(8),POINTER :: RHO(:,:,:)    !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: GAMMA(:,:,:)  !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: HRHO(:,:,:)   !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: SINV(:,:,:)   !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: SMAT(:,:,:)   !(NCHI,NCHI,NDIMD)
END TYPE KSET_TYPE
LOGICAL(4),PARAMETER   :: TON=.TRUE.
LOGICAL(4),SAVE        :: TINI=.FALSE.
REAL(8)   ,PARAMETER   :: AMIX=1.D-1
INTEGER(4)             :: NOMEGA
INTEGER(4)             :: NLAU          ! #(LAURENT EXPANSION TERMS SELF ENERGY
INTEGER(4)             :: NCHI          ! #(CORRELATED ORBITALS)
INTEGER(4)             :: NB            ! #(BAND STATES PER K-POINT)
INTEGER(4)             :: NKPTL         ! #(KPOINTS ON THIS TASK)
INTEGER(4)             :: NSPIN         ! #(SPIN COMPONENTS)
INTEGER(4)             :: NDIM          ! #(SPINOR COMPONENTS)
INTEGER(4)             :: NDIMD         ! CAN BE 1,2,4
INTEGER(4)             :: NAT           ! #(ATOMS)
REAL(8)   ,ALLOCATABLE :: OMEGA(:)      ! MATSUBARA FREQUENCIES
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
     &                       ,NOMEGA,NLAU,KBT,MU,OMEGA &
     &                       ,KSET,ATOMSET
      USE WAVES_MODULE, ONLY : KMAP,NDIM_W=>NDIM,NKPTL_W=>NKPTL,NSPIN_W=>NSPIN
      USE LMTO_MODULE, ONLY: HYBRIDSETTING,HFWEIGHT,POTPAR1
      IMPLICIT NONE
      REAL(8)                :: PI
      INTEGER(4)             :: NTASKS_K,THISTASK_K
      INTEGER(4)             :: NTASKS_M,THISTASK_M
      INTEGER(4)             :: NKPT
      INTEGER(4)             :: NLOC
      INTEGER(4)             :: NU,ISP,IAT,IKPTL,IKPT
      REAL(8)   ,ALLOCATABLE :: WKPT(:) !(NKPT) K-POINT WEIGHTS
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
      NOMEGA=100 ! #(POSITIVE MATSUBARA FREQUENCIES)
      NLAU=1     ! NUMBER OF LAURENT EXPANSION TERMS OF THE SELF ENERGY
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
!!$PRINT*,'KBT=',KBT,' KBT[EV]=',KBT*27.211D0
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
!
!       == SAVE BOUNDS TO ATOMSET ==============================================
        NLOC=SUM(2*POTPAR1(ISP)%LOFH+1)
        ATOMSET(IAT)%NLOC=NLOC
        ATOMSET(IAT)%ICHI1=NCHI+1      ! FIRST ICHI FOR THIS ATOM
        ATOMSET(IAT)%ICHI2=NCHI+NLOC   ! LAST ICHI FOR THIS ATOM
        NCHI=NCHI+NLOC
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
        ALLOCATE(ATOMSET(IAT)%GLOCLAUR(NLOC,NLOC,NDIMD,NLAU+1))
        ALLOCATE(ATOMSET(IAT)%DPHIDG(NLOC,NLOC,NDIMD,NOMEGA))
        ALLOCATE(ATOMSET(IAT)%DPHIDGLAUR(NLOC,NLOC,NDIMD,NLAU))
        ALLOCATE(ATOMSET(IAT)%SLOC(NLOC,NLOC,NDIMD,NOMEGA))
        ALLOCATE(ATOMSET(IAT)%SLOCLAUR(NLOC,NLOC,NDIMD,NLAU))
        ALLOCATE(ATOMSET(IAT)%SMAT(NLOC,NLOC,NDIMD))
        ALLOCATE(ATOMSET(IAT)%DENMAT%RHO(NLOC,NLOC,NDIMD))
        ALLOCATE(ATOMSET(IAT)%DENMAT%H(NLOC,NLOC,NDIMD))
        ATOMSET(IAT)%U=0.D0
        ATOMSET(IAT)%GLOC=(0.D0,0.D0)
        ATOMSET(IAT)%GLOCLAUR=(0.D0,0.D0)
        ATOMSET(IAT)%SLOC=(0.D0,0.D0)
        ATOMSET(IAT)%SLOCLAUR=(0.D0,0.D0)
        ATOMSET(IAT)%SMAT=(0.D0,0.D0)
        ATOMSET(IAT)%DENMAT%RHO =(0.D0,0.D0)
        ATOMSET(IAT)%DENMAT%H   =(0.D0,0.D0)
      ENDDO
PRINT*,'NCHI ',NCHI
!
!     ==========================================================================
!     == COLLECT K-POINT WEIGHTS                                              ==
!     ==========================================================================
      ALLOCATE(KSET(NKPTL))
      DO IKPTL=1,NKPTL
        ALLOCATE(KSET(IKPTL)%F(NB,NSPIN))
        ALLOCATE(KSET(IKPTL)%PIPSI(NDIM,NCHI,NB,NSPIN))
        ALLOCATE(KSET(IKPTL)%RHO(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%GAMMA(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%HRHO(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%SINV(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%SMAT(NCHI,NCHI,NDIMD))
        KSET(IKPTL)%WKPT=0.D0
        KSET(IKPTL)%F=0.D0
        KSET(IKPTL)%PIPSI=(0.D0,0.D0)
        KSET(IKPTL)%RHO=(0.D0,0.D0)
        KSET(IKPTL)%HRHO=(0.D0,0.D0)
        KSET(IKPTL)%GAMMA=(0.D0,0.D0)
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
      CALL DMFT_COLLECTNATORB()
      CALL DMFT_RHOOFK()
      CALL DMFT_RHOLOCAL()

!     CALL DMFT_COLLECTHAMILTONIAN()   ! ABSORBED IN COLLECTNATORB AND RHOOFK
!     CALL DMFT_COLLECTFULLDENMAT()   !ABSORBED IN RHOLOCAL
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
!      CALL DMFT_SMAT()   ! ABSORBED IN DMFT_RHOOFK
!
!     ==========================================================================
!     ==  CONSTRUCT NON-INTERACTING HAMILTONIAN THAT PRODUCES THE CORRECT     ==
!     ==  ONE-PARTICLE DENSITY MATRIX                                         ==
!     ==========================================================================
!CALL DMFT_EXPLOREMODULE()
      CALL DMFT_HRHO()
      CALL DMFT_CONSTRAINTS()
!
!     ==========================================================================
!     ==  CONSTRUCT LOCAL NATURAL ORBITALS                                    ==
!     ==========================================================================
      CALL DMFT_NATORB()
!
!     ==========================================================================
!     == ITERATION TO ENFORCE CONSTRAINTS                                     ==
!     ==========================================================================
      DO ITER=1,3
WRITE(*,FMT='(82("="),T20," ITERATION ",I5)')ITER
        CALL DMFT_GLOC() 
!
!       ========================================================================
!       ==  CALL THE SOLVER                                                   ==
!       ========================================================================
        CALL DMFT_SOLVER(ETOT) 
PRINT*,'ETOT AFTER DMFT_SOLVER ',ETOT
!
!       ========================================================================
!       ==  MIXING                                                            ==
!       ========================================================================
        CALL DMFT_MIX()
!
!       ========================================================================
!       ==  CONSTRAINTS                                                       ==
!       ========================================================================
        CALL DMFT_CONSTRAINTS()
        PRINT*,'ITERATION COMPLETED ',ITER
      ENDDO ! END OF LOOP OVER ITERATIONS TO ENFORCE CONSTRAINT
!
!     ==========================================================================
!     ==  ADD HARTREE FOCK AND DFT DOUBLE COUNTING                            ==
!     ==========================================================================
      CALL DMFT_STATICSOLVER(SVAR)
      ETOT=ETOT+SVAR
! 
!     ==========================================================================
!     ==  CALCULATE MISSING TOTAL ENERGY CONTRIBUTION                         ==
!     ==========================================================================
PRINT*,'WARNING! DMFT_DETOT SWITCHED OFF!!!!'
!      CALL DMFT_DETOT(SVAR)
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
      SUBROUTINE DMFT_COLLECTNATORB()
!     **************************************************************************
!     ** COLLECTS THE PROJECTIONS OF THE EXTENDED NATURAL ORBITALS AND THEIR  **
!     ** OCCUPATION. THIS IS ALL THE INFORMATION REQUIRED FROM THE WAVES OBJECT*
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,NDIMD &
     &                      ,KSET
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IB,J
      REAL(8)                :: F(NB,NKPTL,NSPIN)
      REAL(8)                :: WKPT(NKPTL)
      INTEGER(4)             :: NB_,NSPIN_
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                      CALL TRACE$PUSH('DMFT_COLLECTNATORB')
!
!     ==========================================================================
!     == COLLECT ARRAY SIZES                                                  ==
!     == ASSUMES THAT NKPTL,NB,NSPIN ARE KNOWN                                ==
!     ==========================================================================
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
!     ==========================================================================
!     == COLLECT OCCUPATIONS K-POINT WEIGHTS                                  ==
!     ==========================================================================
      CALL WAVES_DYNOCCGETR8A('WKPT',NKPTL,WKPT)
      DO IKPT=1,NKPTL
        KSET(IKPT)%WKPT=WKPT(IKPT)
      ENDDO
!
!     ==========================================================================
!     == COLLECT OCCUPATIONS WITHOUT K-POINT WEIGHT AND SPIN MULTIPLICITY     ==
!     ==========================================================================
!     == OCCUPATIONS CONTAINS SPIN MULTIPLICITY AND K-POINT WEIGHT =============
      CALL WAVES_DYNOCCGETR8A('OCC',NB*NKPTL*NSPIN,F)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          KSET(IKPT)%F(:,ISPIN)=F(:,IKPT,ISPIN)/KSET(IKPT)%WKPT
!         == OCCUPATIONS ARE IN [0,1] ==========================================
          IF(NSPIN.EQ.1) KSET(IKPT)%F(:,ISPIN)=0.5D0*KSET(IKPT)%F(:,ISPIN)
        ENDDO
      ENDDO
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
            DO IBH=1,NBH
              IF(NBH.NE.NB) THEN
                KSET(IKPT)%PIPSI(:,ICHI,2*IBH-1,ISPIN) &
     &                                          = REAL(THIS%TBC_NEW(:,IBH,ICHI))
                KSET(IKPT)%PIPSI(:,ICHI,2*IBH  ,ISPIN)& 
     &                                          =AIMAG(THIS%TBC_NEW(:,IBH,ICHI))
              ELSE
                KSET(IKPT)%PIPSI(:,ICHI,IBH,ISPIN)=THIS%TBC_NEW(:,IBH,ICHI)
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_RHOOFK()
!     **************************************************************************
!     **  CALCULATES THE K-DEPENDENT DENSITY MATRIX AND THE OVERLAP MATRIX    **
!     **  FROM THE TB-PROJECTIONS OF THE NATURAL ORBITALS AND THEIR OCCUPATIONS* 
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,NDIMD &
     &                      ,NAT,KSET,ATOMSET
      USE MPE_MODULE
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER   :: TTEST=.FALSE.
      COMPLEX(8),ALLOCATABLE :: RHO(:,:,:) !(NCHI,NCHI,NDIMD)
      INTEGER(4)             :: IKPT,ISPIN,ICHI,IB,J,IAT,IDIM1,IDIM2
      INTEGER(4)             :: IDIMD
      INTEGER(4)             :: I1,I2
      COMPLEX(8)             :: MAT(NCHI,NCHI)
      COMPLEX(8)             :: CSVAR
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                      CALL TRACE$PUSH('DMFT_RHOOFK')
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
!     == SUM UP DENSITY MATRIX =================================================
      DO IKPT=1,NKPTL
        KSET(IKPT)%RHO=(0.D0,0.D0)
        KSET(IKPT)%SINV(:,:,:)=(0.D0,0.D0)
        DO ISPIN=1,NSPIN
          DO IDIM1=1,NDIM
            DO IDIM2=1,NDIM
              MAT(:,:)=(0.D0,0.D0)
              DO IB=1,NB
                DO J=1,NCHI
                  CSVAR=KSET(IKPT)%F(IB,ISPIN) &
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
              KSET(IKPT)%SINV(:,:,IDIMD) &
     &                     =MATMUL(KSET(IKPT)%PIPSI(IDIM1,:,:,ISPIN) &
     &                     ,CONJG(TRANSPOSE(KSET(IKPT)%PIPSI(IDIM2,:,:,ISPIN))))
            ENDDO
          ENDDO
        ENDDO
!
!       ========================================================================
!       == CONVERT FROM UP-DOWN REPRESENTATION INTO (TXYZ)                    ==
!       ========================================================================
        CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%RHO)
        CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%SINV)
!
!       ========================================================================
!       == OVERLAP MATRIX                                                     ==
!       ========================================================================
        CALL SPINOR$INVERT(NDIMD,NCHI,KSET(IKPT)%SINV,KSET(IKPT)%SMAT)
      ENDDO

                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_RHOLOCAL()
!     **************************************************************************
!     ==  EXTRACT ONSITE DENSITY MATRICES WITH ALL PROJECTOR FUNCTIONS        ==
!     ==  FOR DOUBLE COUNTING                                                 ==
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NKPTL,NSPIN,NAT,NDIM,NDIMD,NB,ATOMSET,KSET
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: NLOC
      INTEGER(4)             :: I1,I2
      COMPLEX(8)             :: CSVAR
      INTEGER(4)             :: NB_,NSPIN_
      INTEGER(4)             :: IAT,IKPT,ISPIN,IPRO,IDIM1,IDIM2,IDIMD,I
      INTEGER(4)             :: LMN
!     **************************************************************************
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
          IPRO=0
          DO IAT=1,NAT
            NLOC=ATOMSET(IAT)%NLOC
            IF(NLOC.LE.0) THEN
              IPRO=IPRO+NLOC
              CYCLE
            END IF   
            I1=IPRO+1         
            I2=IPRO+NLOC
            ATOMSET(IAT)%DENMAT%RHO(:,:,:)=ATOMSET(IAT)%DENMAT%RHO(:,:,:) &
     &                           +KSET(IKPT)%WKPT*KSET(IKPT)%RHO(I1:I2,I1:I2,:)
            IPRO=IPRO+NLOC
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
        NLOC=ATOMSET(IAT)%NLOC
        IF(NLOC.LE.0) CYCLE
        CALL SPINOR$CONVERT('FWRD',NLOC,NDIMD,ATOMSET(IAT)%DENMAT%RHO)
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
          NLOC=ATOMSET(IAT)%NLOC
          IF(NLOC.LE.0) CYCLE
          WRITE(*,FMT='(82("="),T10," DENSITY MATRIX FOR ATOM ",I3," ")')IAT
          DO IDIMD=1,NDIMD
            DO LMN=1,NLOC
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
      SUBROUTINE DMFT_COLLECTHAMILTONIAN()
!     **************************************************************************
!     ** COLLECTS THE HAMILTONIAN AND STORES IT ON THE MODULE                 **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,NDIMD &
     &                      ,NAT,KSET,ATOMSET
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER   :: TTEST=.FALSE.
      COMPLEX(8),ALLOCATABLE :: RHO(:,:,:) !(NCHI,NCHI,NDIMD)
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IB,J,IAT,IDIM1,IDIM2
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
            DO IBH=1,NBH
              IF(NBH.NE.NB) THEN
                KSET(IKPT)%PIPSI(:,ICHI,2*IBH-1,ISPIN) &
     &                                          = REAL(THIS%TBC_NEW(:,IBH,ICHI))
                KSET(IKPT)%PIPSI(:,ICHI,2*IBH  ,ISPIN)& 
     &                                          =AIMAG(THIS%TBC_NEW(:,IBH,ICHI))
              ELSE
                KSET(IKPT)%PIPSI(:,ICHI,IBH,ISPIN)=THIS%TBC_NEW(:,IBH,ICHI)
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
      INTEGER(4)             :: NLOC
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
            NLOC=ATOMSET(IAT)%NLOC
            IF(NLOC.LE.0) THEN
              IPRO=IPRO+NLOC
              CYCLE
            END IF   
            I1=IPRO+1         
            I2=IPRO+NLOC
            DO IBH=1,NBH
              IF(NBH.NE.NB) THEN
!               == THIS IS FOR GENERAL SPECIAL K-POINTS ========================
!               == TBC IS FROM A SUPER WAVE FUNCTION AND CONTAINS TWO WAVE =====
!               == FUNCTIONS.
                DO IDIM2=1,NDIM
                  DO IDIM1=1,NDIM
                    IDIMD=IDIM1+NDIM*(IDIM2-1)+(ISPIN-1)
                    DO LMN=1,NLOC
                      I=IPRO+LMN
                      CSVAR=REAL(THIS%TBC_NEW(IDIM1,IBH,I))*OCC(2*IBH-1,IKPT,ISPIN)
                      ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
     &                               =ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
     &                               +REAL(THIS%TBC_NEW(IDIM2,IBH,I1:I2))*CSVAR
!
                      CSVAR=AIMAG(THIS%TBC_NEW(IDIM1,IBH,I))*OCC(2*IBH,IKPT,ISPIN)
                      ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
     &                             =ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
     &                             +AIMAG(THIS%TBC_NEW(IDIM2,IBH,I1:I2))*CSVAR
                    ENDDO
                  ENDDO
                ENDDO
              ELSE
!               == THIS IS FOR GENERAL K-POINTS ================================
                DO IDIM2=1,NDIM
                  DO IDIM1=1,NDIM
                    IDIMD=IDIM1+NDIM*(IDIM2-1)+(ISPIN-1)
                    DO LMN=1,NLOC
                      I=IPRO+LMN
                      CSVAR=CONJG(THIS%TBC_NEW(IDIM1,IBH,I))*OCC(IBH,IKPT,ISPIN)
                      ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
     &                             =ATOMSET(IAT)%DENMAT%RHO(:,LMN,IDIMD) &
     &                             +THIS%TBC_NEW(IDIM2,IBH,I1:I2)*CSVAR
                    ENDDO
                  ENDDO
                ENDDO
              END IF
            ENDDO  ! END OF LOOP OVER BANDS
            IPRO=IPRO+NLOC
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
        NLOC=ATOMSET(IAT)%NLOC
        IF(NLOC.LE.0) CYCLE
        CALL SPINOR$CONVERT('FWRD',NLOC,NDIMD,ATOMSET(IAT)%DENMAT%RHO)
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
          NLOC=ATOMSET(IAT)%NLOC
          IF(NLOC.LE.0) CYCLE
          WRITE(*,FMT='(82("="),T10," DENSITY MATRIX FOR ATOM ",I3," ")')IAT
          DO IDIMD=1,NDIMD
            DO LMN=1,NLOC
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
      SUBROUTINE DMFT_UTENSOR()
!     **************************************************************************
!     ** CALCULATE THE U-TENSOR OF THE CORRELATED ORBITALS IN THE SELECTED SET
!     **************************************************************************
      USE LMTO_MODULE, ONLY: ISPECIES
      USE DMFT_MODULE, ONLY: NAT,ATOMSET
      IMPLICIT NONE
      INTEGER(4)            :: ISP ! ATOM TYPE
      INTEGER(4)            :: IAT
      INTEGER(4)            :: NH
!     **************************************************************************
                                          CALL TRACE$PUSH('DMFT_UTENSOR')
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        NH=ATOMSET(IAT)%NLOC
!
!       ========================================================================
!       ==  CALCULATE U-TENSOR IN THE BASIS OF LOCAL ORBITALS IGNORING TORB   ==
!       ========================================================================
        CALL DMFT_ULOCAL(IAT,NH,ATOMSET(IAT)%U)
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
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR1,SBAR_NEW
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      INTEGER(4),INTENT(IN)  :: LMNX          !  #(LOCAL ORBITALS ON THIS SITE)
      REAL(8)   ,INTENT(OUT) :: U(LMNX,LMNX,LMNX,LMNX)
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      REAL(8)   ,ALLOCATABLE :: UT(:,:,:,:)
      INTEGER(4)             :: ISP     ! ATOM TYPE
      INTEGER(4)             :: LMNXT   ! #(VALENCE+SCATTERING WAVES)
      INTEGER(4)             :: NT      ! #(TAIL FUNCTIONS(LMN))
      INTEGER(4)             :: NH      ! #(HEAD FUNCTIONS (LMN))
      INTEGER(4)             :: NNS     ! #(PAIRS FOR STRUCTURE CONSTANTS)
      INTEGER(4)             :: IH,IT,NN,NN0
!     **************************************************************************
      ISP=ISPECIES(IAT)
      LMNXT=POTPAR1(ISP)%TAILED%LMNX  ! SIZE OF U-TENSOR ON POTPAR
      NH=SUM(2*POTPAR1(ISP)%LOFH+1)     
      NT=SUM(2*POTPAR1(ISP)%LOFT+1)     
      IF(NH+NT.NE.LMNXT) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
        CALL ERROR$MSG('LMNXT=NH+NT MUST BE OBEYED')
        CALL ERROR$I4VAL('LMNXT',LMNXT)
        CALL ERROR$I4VAL('NH',NH)
        CALL ERROR$I4VAL('NT',NT)
        CALL ERROR$I4VAL('ISP',ISP)
        CALL ERROR$STOP('DMFT_ULOCAL')
      END IF
      IF(NH.NE.LMNX) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
        CALL ERROR$MSG('LMNX=NH MUST BE OBEYED')
        CALL ERROR$I4VAL('LMNX',LMNX)
        CALL ERROR$I4VAL('NH',NH)
        CALL ERROR$I4VAL('ISP',ISP)
        CALL ERROR$STOP('DMFT_ULOCAL')
      END IF
!
!     ==========================================================================
!     == COLLECT U-TENSOR IN EXPANDED BASIS                                   ==
!     ==========================================================================
      ALLOCATE(UT(LMNXT,LMNXT,LMNXT,LMNXT))
      UT=POTPAR1(ISP)%TAILED%U
!
!     ==========================================================================
!     == COLLECT STRUCTURE CONSTANTS                                          ==
!     ==========================================================================
      NNS=SIZE(SBAR_NEW)
      DO NN=1,NNS
        IF(SBAR_NEW(NN)%IAT1.NE.IAT) CYCLE
        IF(SBAR_NEW(NN)%IAT2.NE.IAT) CYCLE
        IF(SUM(SBAR_NEW(NN)%IT(:)**2).NE.0) CYCLE
        NN0=NN
        EXIT
      ENDDO
      NN=NN0
!
!     __ CHECK DIMENSIONS_______________________________________________________
      IF(SBAR_NEW(NN)%N1.NE.NH.OR.SBAR_NEW(NN)%N2.NE.NT) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS FOR SBAR')
        CALL ERROR$MSG('NH=SBAR%N1 AND NT=SBAR%N2 MUST BE OBEYD')
        CALL ERROR$I4VAL('NH',NH) ! INDEX OF HEAD FUNCTIONS
        CALL ERROR$I4VAL('NT',NT) ! INDEX OF HEAD FUNCTIONS
        CALL ERROR$I4VAL('SBAR%N1',SBAR_NEW(NN)%N1) 
        CALL ERROR$I4VAL('SBAR%N2',SBAR_NEW(NN)%N2) 
        CALL ERROR$STOP('DMFT_ULOCAL')
      END IF
!
!     ==========================================================================
!     == DOWNFOLD UTENSOR
!     ==   |KBAR_I>=|K_I>-\SUM_J |JBAR_J> SBAR(I,J)                           ==
!     ==========================================================================
!     U(:,:,::)=POTPAR(ISP)%TAILED%U(:LMNX,:LMNX,:LMNX,:LMNX)

      DO IH=1,NH
        DO IT=1,NT
          UT(:,:,:,IH)=UT(:,:,:,IH)-SBAR_NEW(NN)%MAT(IH,IT)*UT(:,:,:,NH+IT)
        ENDDO
      ENDDO
      DO IH=1,NH
        DO IT=1,NT
          UT(:,:,IH,:NH)=UT(:,:,IH,:NH) &
     &                   -SBAR_NEW(NN)%MAT(IH,IT)*UT(:,:,NH+IT,:NH)
        ENDDO
      ENDDO
      DO IH=1,NH
        DO IT=1,NT
          UT(:,IH,:NH,:NH)=UT(:,IH,:NH,:NH) &
     &                    -SBAR_NEW(NN)%MAT(IH,IT)*UT(:,NH+IT,:NH,:NH)
        ENDDO
      ENDDO
      DO IH=1,NH
        DO IT=1,NT
          UT(IH,:NH,:NH,:NH)=UT(IH,:NH,:NH,:NH) &
     &                     -SBAR_NEW(NN)%MAT(IH,IT)*UT(NH+IT,:NH,:NH,:NH)
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
      CHARACTER(16),PARAMETER:: ORBTYPE='QUAMBO' ! CAN BE 'NAT','QUAMBO','ORTHO'
      LOGICAL(4)             :: NATORB
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
      IF(ORBTYPE.EQ.'NATORB') THEN
        NATORB=.TRUE.
      ELSE IF(ORBTYPE.EQ.'QUAMBO') THEN
        NATORB=.TRUE.
      ELSE IF(ORBTYPE.EQ.'ORTHO') THEN
        NATORB=.FALSE.
      ELSE
        CALL ERROR$MSG('ORBTYPE NOT RECOGNIZED')
        CALL ERROR$MSG('ALLOWED VALUES ARE "NATORB", "QUAMBO" AND "ORTHO"')
        CALL ERROR$CHVAL('ORBTYPE',ORBTYPE)
        CALL ERROR$STOP('DMFT_NATORB')
      END IF

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
!       == OVERWRITE RHO FOR QUAMBO CONSTRUCTION                              ==
!       ========================================================================
        IF(ORBTYPE.EQ.'QUAMBO') THEN
          RHO=(0.D0,0.D0)
          DO I=1,NLOC
            RHO(I,I,1)=CMPLX(I,KIND=8)
          ENDDO
          IF(NDIMD.GT.1) THEN
            DO I=1,NLOC
              RHO(I,I,NDIMD)=-CMPLX(I,KIND=8)
            ENDDO
          END IF
        END IF
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
      SUBROUTINE DMFT_HRHO()
!     **************************************************************************
!     **  CALCULATES A STATIC HAMILTONIAN HRHO SUCH THAT THE RESULTING        **
!     **  DENSITY MATRIX IS RHO=[1+EXP(BETA*HRHO)]^(-1)                       **
!     **************************************************************************
      USE DMFT_MODULE ,ONLY: NKPTL,NCHI,NDIMD,NOMEGA,NLAU,OMEGA,KBT,MU,KSET
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
              SVAR=MU+SVAR
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
            SVAR=MU+SVAR
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
          CALL DMFT_HRHO_TEST(NCHI,NDIMD,NOMEGA,NLAU,OMEGA,KBT &
     &                       ,KSET(IKPT)%SMAT,KSET(IKPT)%HRHO,MAT) !MAT=RHO
!         == CALCULATE DEVIATION FROM TARGET DENSITY MATRIX ====================
          MAT=MAT-KSET(IKPT)%RHO
          IF(MAXVAL(ABS(MAT)).GT.TOL) THEN
            CALL SPINOR_PRINTMATRIX(6,'RHO[HRHO]-RHO',1,NCHI &
     &                           ,NDIMD,NCHI,MAT)
            CALL ERROR$MSG('TEST OF HRHO FAILED')
            CALL ERROR$MSG('INCREASE THE NUMBER OF MATSUBARA FREQUENCIES')
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
      SUBROUTINE DMFT_HRHO_TEST(NCHI,NDIMD,NOMEGA,NLAU,OMEGA,KBT,S,H,RHO)
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
      INTEGER(4),INTENT(IN)    :: NLAU
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
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,NLAU+1,FN)
      CALL SPINOR$INVERT(NDIMD,NCHI,S,SINV)
!     == GLAUR1=SINV ===========================================================
      RHO=RHO+FN(1)*SINV
      IF(NLAU.GE.1) THEN
!       == GLAUR2=SINV*H*SINV ==================================================
        CALL SPINOR$MATMUL(NDIMD,NCHI,H,SINV,MAT1)
        CALL SPINOR$MATMUL(NDIMD,NCHI,SINV,MAT1,MAT2)
        RHO=RHO+FN(2)*MAT2
      END IF
      IF(NLAU.GE.2) THEN
!       == GLAUR3=SINV*H*SINV*H*SINV ===========================================
!       == GLAUR3 CONTRIBUTES NOTHING
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,H,MAT1)
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,SINV,MAT2)
        RHO=RHO+FN(3)*MAT2
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_GLOC()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NDIMD,NAT,NOMEGA,NLAU,OMEGA,MU &
     &                      ,KSET,ATOMSET
      IMPLICIT NONE
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER     :: TPRINT=.FALSE.
      LOGICAL(4),PARAMETER     :: TTEST=.FALSE.
      COMPLEX(8)               :: MATX(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR(NCHI,NCHI,NDIMD,NLAU+1)
      COMPLEX(8)               :: SLAUR(NCHI,NCHI,NDIMD,NLAU)
      COMPLEX(8)               :: G(NCHI,NCHI,NDIMD)
      REAL(8)                  :: WKPTL
      INTEGER(4)               :: IKPT,NU,IAT,I1,I2,IDIMD,ILAU
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_GLOC')
      IF(NLAU.GT.1) THEN
        CALL ERROR$MSG('NLAU EXCEEDS MAXIMUM VALUE OF 1')
        CALL ERROR$I4VAL('NLAU',NLAU)
        CALL ERROR$STOP('DMFT_GLOC')
      END IF
!
      DO IAT=1,NAT
        ATOMSET(IAT)%GLOC=(0.D0,0.D0)
        ATOMSET(IAT)%GLOCLAUR=(0.D0,0.D0)
      ENDDO
!
!     ==========================================================================
!     == CONSTRUCT LAURENT EXPANSION SELF-ENERGY, WHICH IS K-INDEPENDENT      ==
!     ==========================================================================
      SLAUR=(0.D0,0.D0)
      DO IAT=1,NAT
        I1=ATOMSET(IAT)%ICHI1
        I2=ATOMSET(IAT)%ICHI2
        DO ILAU=1,NLAU
          SLAUR(I1:I2,I1:I2,:,ILAU)=ATOMSET(IAT)%SLOCLAUR(:,:,:,ILAU)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == PERFORM BRILLOUIN-ZONE SAMPLING                                      ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        WKPTL=KSET(IKPT)%WKPT
!       == LAURENT EXPANSION FOR THE GREENS FUNCTION ===========================
        GLAUR(:,:,:,1)=KSET(IKPT)%SINV
!
        MATX=-MU*KSET(IKPT)%SMAT+KSET(IKPT)%HRHO-KSET(IKPT)%GAMMA+SLAUR(:,:,:,1)
!       == GLAUR(2)=SINV*MAT*SINV ==============================================
        CALL SPINOR$MATMUL(NDIMD,NCHI,MATX,KSET(IKPT)%SINV,MAT)
        CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,MAT,GLAUR(:,:,:,2))
        IF(KSET(IKPT)%TADDMINUSK) THEN
          DO IDIMD=1,NDIMD
            DO ILAU=1,NLAU+1
              GLAUR(:,:,IDIMD,ILAU)=0.5D0*(GLAUR(:,:,IDIMD,ILAU) &
      &                   +CONJG(TRANSPOSE(GLAUR(:,:,IDIMD,ILAU))))
            ENDDO
          ENDDO
        END IF
!       == PROJECT LAURENT TERMS ONTO LOCAL HILBERT SPACE ======================
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          DO ILAU=1,NLAU
            ATOMSET(IAT)%GLOCLAUR(:,:,:,ILAU)=ATOMSET(IAT)%GLOCLAUR(:,:,:,ILAU)&
     &                                     +WKPTL*GLAUR(I1:I2,I1:I2,:,ILAU)
          ENDDO
        ENDDO
!
        DO NU=1,NOMEGA
!         == CONSTRUCT LATTICE GREENS FUNCTION =================================
          MAT=(CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT
          MAT=MAT-KSET(IKPT)%HRHO+KSET(IKPT)%GAMMA
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            MAT(I1:I2,I1:I2,:)=MAT(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
          ENDDO
          CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
!         == ACCOUNT FOR K-POINT (-K) ==========================================
          IF(KSET(IKPT)%TADDMINUSK) THEN
!           == TRANSPOSE SMAT,HRHO,GAMMA, BUT NEITHER CI NOR SLOC ==============
            MAT=(-CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT
            MAT=MAT-KSET(IKPT)%HRHO+KSET(IKPT)%GAMMA
            DO IDIMD=1,NDIMD
              MAT(:,:,IDIMD)=TRANSPOSE(CONJG(MAT(:,:,IDIMD)))
            ENDDO
            DO IAT=1,NAT
              I1=ATOMSET(IAT)%ICHI1
              I2=ATOMSET(IAT)%ICHI2
              MAT(I1:I2,I1:I2,:)=MAT(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
            ENDDO
            CALL SPINOR$INVERT(NDIMD,NCHI,MAT,MAT2)
            G=0.5D0*(G+MAT2)
          END IF
!         == PROJECT GREENS FUNCTION ONTO INFIVIDUAL SITES =====================
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
        DO ILAU=1,NLAU+1
          DO IDIMD=1,NDIMD
            ATOMSET(IAT)%GLOCLAUR(:,:,IDIMD,ILAU)=0.5D0 &
     &                   *(ATOMSET(IAT)%GLOCLAUR(:,:,IDIMD,ILAU) &
     &                 +CONJG(TRANSPOSE(ATOMSET(IAT)%GLOCLAUR(:,:,IDIMD,ILAU))))
          ENDDO
        ENDDO
      ENDDO
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVER(ETOT)
!     **************************************************************************
!     **  CALL DMFT INTERFACE FOR DYNAMIC CORRELATIONS (LOOP OVER ATOMS)      **
!     **  PRODUCES ETOT                                                       **
!     **           ATOMSET(IAT)%DPHIDG                                        **
!     **           ATOMSET(IAT)%DPHIDGLAUR                                    **
!     **           ATOMSET(IAT)%DEDU                                          **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE ,ONLY: NAT,NDIMD,NOMEGA,NLAU,KBT,ATOMSET
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: ETOT
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      REAL(8)                :: PHILW
      INTEGER(4)             :: NLOC  !#(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4)             :: IAT
!     **************************************************************************
                                      CALL TRACE$PUSH('DMFT_SOLVER')
      ETOT=0.D0
      DO IAT=1,NAT
        NLOC=ATOMSET(IAT)%NLOC
        IF(NLOC.LE.0) CYCLE
!
!       ========================================================================
!       == ADD DYNAMIC CONTRIBUTIONS  (WITHOUT HF CONTRIBUTION)               ==
!       == HFWEIGHT IS ABSORBED IN THE LOCAL U-TENSOR                         ==
!       ========================================================================
!IF(IAT.EQ.2)CALL DMFT_SOLVERTEST(IAT)
        CALL DMFT_DYNAMICSOLVER(NLOC,NDIMD,NOMEGA,NLAU,KBT &
     &                         ,ATOMSET(IAT)%GLOC,ATOMSET(IAT)%GLOCLAUR &
     &                         ,ATOMSET(IAT)%U,PHILW &
     &                         ,ATOMSET(IAT)%DPHIDG,ATOMSET(IAT)%DPHIDGLAUR &
     &                         ,ATOMSET(IAT)%DEDU &
     &                         ,ATOMSET(IAT)%NATORB%PIPHI &
     &                         ,ATOMSET(IAT)%NATORB%CHIPHI)
        ETOT=ETOT+PHILW !SCREENING ALREADY CONTAINED IN U-TENSOR
      ENDDO ! END OF LOOP OVER ATOMS (IAT)
      WRITE(*,FMT='(50("."),T1'&
     &        //',"LW FUNCTIONAL FOR DYNAMIC CORRELATION OF ATOM ",I3,T50,":"' &
     &        //',F20.10)')IAT,ETOT
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVERTEST(IAT)
      USE DMFT_MODULE ,ONLY: NDIMD,NOMEGA,NLAU,KBT,ATOMSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),PARAMETER  :: NSCALE=3
      REAL(8)               :: SCALE(NSCALE)
      REAL(8)               :: Y0 
      REAL(8)               :: PHILW
      INTEGER(4)            :: NLOC       
      REAL(8)               :: SIGN
      COMPLEX(8),ALLOCATABLE:: GLOCSAVE(:,:,:,:)
      COMPLEX(8),ALLOCATABLE:: DGLOC(:,:,:,:)
      COMPLEX(8),ALLOCATABLE:: SIGMA(:,:,:,:)
      INTEGER(4)            :: I,IPM
!     **************************************************************************
      NLOC=ATOMSET(IAT)%NLOC
      ALLOCATE(GLOCSAVE(NLOC,NLOC,NDIMD,NOMEGA))
      ALLOCATE(SIGMA(NLOC,NLOC,NDIMD,NOMEGA))
      CALL DMFT_DYNAMICSOLVER(NLOC,NDIMD,NOMEGA,NLAU,KBT &
     &                         ,ATOMSET(IAT)%GLOC,ATOMSET(IAT)%GLOCLAUR &
     &                         ,ATOMSET(IAT)%U,PHILW &
     &                         ,ATOMSET(IAT)%DPHIDG,ATOMSET(IAT)%DPHIDGLAUR &
     &                         ,ATOMSET(IAT)%DEDU &
     &                         ,ATOMSET(IAT)%NATORB%PIPHI &
     &                         ,ATOMSET(IAT)%NATORB%CHIPHI)
      SIGMA=ATOMSET(IAT)%DPHIDG
      GLOCSAVE=ATOMSET(IAT)%GLOC
!
      ALLOCATE(DGLOC(NLOC,NLOC,NDIMD,NOMEGA))
      DGLOC=(0.D0,0.D0)
      DGLOC(1,1,1,1)=(0.D0,1.D0)

      Y0=2.D0*KBT*REAL(SUM(SIGMA*DGLOC))
WRITE(*,FMT='("Y(",I3,")=",2F20.10)')0,Y0
!
      SCALE=(/1.D0,2.D0,4.D0/)
      DO I=1,NSCALE
        Y0=0.D0
        DO IPM=-1,1,2
          SIGN=REAL(IPM,KIND=8)
          ATOMSET(IAT)%GLOC=GLOCSAVE+DGLOC*SCALE(I)*SIGN
          CALL DMFT_DYNAMICSOLVER(NLOC,NDIMD,NOMEGA,NLAU,KBT &
     &                         ,ATOMSET(IAT)%GLOC,ATOMSET(IAT)%GLOCLAUR &
     &                         ,ATOMSET(IAT)%U,PHILW &
     &                         ,ATOMSET(IAT)%DPHIDG,ATOMSET(IAT)%DPHIDGLAUR &
     &                         ,ATOMSET(IAT)%DEDU &
     &                         ,ATOMSET(IAT)%NATORB%PIPHI &
     &                         ,ATOMSET(IAT)%NATORB%CHIPHI)
          Y0=Y0+PHILW*SIGN/(2.D0*SCALE(I))
        ENDDO
WRITE(*,FMT='("Y(",I3,")=",2F20.10)')I,Y0
      ENDDO
STOP 'FORCED'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DMFT_DYNAMICSOLVER(NLOC,NDIMD,NOMEGA,NLAU,KBT &
     &                        ,GLOC,GLOCLAUR,UCHI,ETOT,SLOC,SLOCLAUR,DEDUCHI &
     &                        ,PIPHI,CHIPHI)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NLOC
      INTEGER(4),INTENT(IN)  :: NDIMD
      INTEGER(4),INTENT(IN)  :: NOMEGA
      INTEGER(4),INTENT(IN)  :: NLAU   !#(LAURENT EXP. TERMS FOR SELF ENERGY)
      REAL(8)   ,INTENT(IN)  :: KBT
      COMPLEX(8),INTENT(IN)  :: PIPHI(2*NLOC,2*NLOC)  !<PI|NATORB>
      COMPLEX(8),INTENT(IN)  :: CHIPHI(2*NLOC,2*NLOC) !<CHI|NATORB>
      COMPLEX(8),INTENT(IN)  :: GLOC(NLOC,NLOC,NDIMD,NOMEGA)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR(NLOC,NLOC,NDIMD,NLAU+1)
      REAL(8)   ,INTENT(IN)  :: UCHI(NLOC,NLOC,NLOC,NLOC)
      REAL(8)   ,INTENT(OUT) :: ETOT
      COMPLEX(8),INTENT(OUT) :: SLOC(NLOC,NLOC,NDIMD,NOMEGA)
      COMPLEX(8),INTENT(OUT) :: SLOCLAUR(NLOC,NLOC,NDIMD,NLAU)
      REAL(8)   ,INTENT(OUT) :: DEDUCHI(NLOC,NLOC,NLOC,NLOC)
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      COMPLEX(8)             :: U(2*NLOC,2*NLOC,2*NLOC,2*NLOC)
      COMPLEX(8)             :: DEDU(2*NLOC,2*NLOC,2*NLOC,2*NLOC)
      COMPLEX(8)             :: V(2*NLOC,2*NLOC,2*NLOC,2*NLOC)
      COMPLEX(8)             :: G(2*NLOC,2*NLOC,NOMEGA)
      COMPLEX(8)             :: GLAUR(2*NLOC,2*NLOC,NLAU+1)
      COMPLEX(8)             :: S(2*NLOC,2*NLOC,NOMEGA)
      COMPLEX(8)             :: SLAUR(2*NLOC,2*NLOC,NLAU)
      INTEGER(4)             :: NU,I,J
!     **************************************************************************
                                     CALL TRACE$PUSH('DMFT_DYNAMICSOLVER')
!
!     ==========================================================================
!     == CONVERT U-TENSOR                                                     ==
!     ==========================================================================
      U(  :NLOC,  :NLOC,  :NLOC,  :NLOC)=CMPLX(UCHI,KIND=8)
      U(NLOC+1:,NLOC+1:,NLOC+1:,NLOC+1:)=CMPLX(UCHI,KIND=8)
      U(  :NLOC,NLOC+1:,  :NLOC,NLOC+1:)=CMPLX(UCHI,KIND=8)
      U(NLOC+1:,  :NLOC,NLOC+1:,  :NLOC)=CMPLX(UCHI,KIND=8)
      V=(0.D0,0.D0)
      DO I=1,2*NLOC
        DO J=1,2*NLOC
          V(I,:,:,:)=V(I,:,:,:)+CONJG(PIPHI(J,I))*U(J,:,:,:)
        ENDDO
      ENDDO
!
      U=(0.D0,0.D0)
      DO I=1,2*NLOC
        DO J=1,2*NLOC
          U(:,I,:,:)=U(:,I,:,:)+CONJG(PIPHI(J,I))*V(:,J,:,:)
        ENDDO
      ENDDO
!
      V=(0.D0,0.D0)
      DO I=1,2*NLOC
        DO J=1,2*NLOC
          V(:,:,I,:)=V(:,:,I,:)+U(:,:,J,:)*PIPHI(J,I)
        ENDDO
      ENDDO
!
      U=(0.D0,0.D0)
      DO I=1,2*NLOC
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
      DO I=1,NLAU+1
        CALL SPINOR$BLOWUP(NDIMD,NLOC,GLOCLAUR(:,:,:,I),GLAUR(:,:,I))
        GLAUR(:,:,I)=MATMUL(TRANSPOSE(CONJG(CHIPHI)) &
     &                      ,MATMUL(GLAUR(:,:,I),CHIPHI))
      ENDDO
!
!     ==========================================================================
!     == WRITE DATA
!     ==========================================================================
      CALL DMFT_SOLVERIO(2*NLOC,NOMEGA,NLAU,KBT,G,GLAUR,U,ETOT,S,SLAUR,DEDU)
!
!     ==========================================================================
!     == CONVERT GREENS FUNCTION                                              ==
!     ==========================================================================
!     == THE FACTOR 1/2 IS NEEDED SO THAT 
!     ==        DPHI=BETA*SUM_NU TR[SUM_IDIM SIGMA(IDIMD)*DG(IDIMD)]
      DO NU=1,NOMEGA
        S(:,:,NU)=MATMUL(CHIPHI,MATMUL(S(:,:,NU),CONJG(TRANSPOSE(CHIPHI))))
        CALL SPINOR$SHRINKDOWN(NDIMD,NLOC,S(:,:,NU),SLOC(:,:,:,NU))
      ENDDO
      SLOC=0.5D0*SLOC
      DO I=1,NLAU
        SLAUR(:,:,I)=MATMUL(CHIPHI,MATMUL(SLAUR(:,:,I) &
     &                                   ,CONJG(TRANSPOSE(CHIPHI))))
        CALL SPINOR$SHRINKDOWN(NDIMD,NLOC,SLAUR(:,:,I),SLOCLAUR(:,:,:,I))
      ENDDO
      SLAUR=0.5D0*SLAUR
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
     &            +U(NLOC+1:,NLOC+1:,NLOC+1:,NLOC+1:) &
     &            +U(NLOC+1:,  :NLOC,NLOC+1:,  :NLOC) &
     &            +U(  :NLOC,NLOC+1:,  :NLOC,NLOC+1:) )
                                     CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVERIO(NORB,NOMEGA,NLAU,KBT,G,GLAUR,U,ETOT,S,SLAUR,DEDU)
!     **************************************************************************
!     ** INTERFACE ROUTINE FOR THE COMMUNICATION WITH THE DMFT SOLVER         **
!     ** THE SPIN-UP ORBITALS ARE STORED FIRST AND THEN THE SPIN-DOWN         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB     !#(SPIN ORBITALS)
      INTEGER(4),INTENT(IN) :: NOMEGA   !#(POSITIVE MATSUBARA FREQUENCIES
      INTEGER(4),INTENT(IN) :: NLAU     !#(LAURENT EXP. TERMS IN SELF ENERGY)
      REAL(8)   ,INTENT(IN) :: KBT      
      COMPLEX(8),INTENT(IN) :: G(NORB,NORB,NOMEGA)    !GREENS FUNCTION
      COMPLEX(8),INTENT(IN) :: GLAUR(NORB,NORB,NLAU+1)!LAURENT EXPANSION OF G
      COMPLEX(8),INTENT(IN) :: U(NORB,NORB,NORB,NORB) !U-TENSOR
      REAL(8)   ,INTENT(OUT):: ETOT                   !LUTTINGER WARD FUNCTIONAL
      COMPLEX(8),INTENT(OUT):: S(NORB,NORB,NOMEGA)    !SELF ENERGY
      COMPLEX(8),INTENT(OUT):: SLAUR(NORB,NORB,NLAU)  !LAURENT EXPANSION OF S
      COMPLEX(8),INTENT(OUT):: DEDU(NORB,NORB,NORB,NORB) 
      INTEGER(4)            :: NORB_,NOMEGA_
      INTEGER(4)            :: NFIL
      CHARACTER(32),PARAMETER :: TYPE='LINEAR' ! 'NONE','LINEAR'
!     **************************************************************************
                                     CALL TRACE$PUSH('DMFT_SOLVERIO')
!CALL TESTG(NORB,NOMEGA,NLAU,KBT,G,GLAUR)
      ETOT=0.D0
      S=(0.D0,0.D0)
      SLAUR=(0.D0,0.D0)
      DEDU=(0.D0,0.D0)
      IF(TYPE.EQ.'NONE') THEN
        PRINT*,'USING STATIC LUTTINGER WARD FUNCTIONAL'
        PRINT*,'(NO DYNAMIC CORRECTION. FOR TESTING PURPOSES ONLY)'
        CALL TRACE$POP()
        RETURN
      ELSE IF(TYPE.EQ.'LINEAR') THEN
        PRINT*,'USING LINEAR LUTTINGER-WARD FUNCTIONAL'
        PRINT*,'(USE FOR TESTING PURPOSES ONLY)'
        CALL DMFT_SOLVER_LINEAR(NORB,NOMEGA,KBT,G,ETOT,S)  
        CALL TRACE$POP()
        RETURN
      END IF
!
!     ==========================================================================
!     == WRITE DATA
!     ==========================================================================
!WARNING!! THE TWO FILES ARE NOT YET DEFINED
!CAREFUL: IN FORTRAN, THE FIRST INDEX RUNS FIRST
      CALL FILEHANDLER$UNIT('DMFTOSOLVER',NFIL)
      REWIND NFIL
      WRITE(NFIL,*)NORB,NOMEGA,NLAU,KBT
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
                                     CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVER_LINEAR(NORB,NOMEGA,KBT,G,ETOT,S)  
!     **************************************************************************
!     ** MIMICKS A SOLVER FOR THE LUTTINGER WARD FUNCTIONAL AND SELF ENERGY   **
!     ** ON THE BASIS OF A LINEAR LUTTINGER WARD FUNCTIONAL, WHICH IS         **
!     ** DEFINED BY A FROZEN SELF ENERGY.                                     **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB     !#(SPIN ORBITALS)
      INTEGER(4),INTENT(IN) :: NOMEGA   !#(POSITIVE MATSUBARA FREQUENCIES
      REAL(8)   ,INTENT(IN) :: KBT      
      COMPLEX(8),INTENT(IN) :: G(NORB,NORB,NOMEGA)    !GREENS FUNCTION
      REAL(8)   ,INTENT(OUT):: ETOT                   !LUTTINGER WARD FUNCTIONAL
      COMPLEX(8),INTENT(OUT):: S(NORB,NORB,NOMEGA)    !SELF ENERGY
      REAL(8)               :: PI
      REAL(8)               :: OMEGA
      INTEGER(4)            :: NU,I
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)               :: SVAR
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      S=(0.D0,0.D0)
      ETOT=0.D0
      DO NU=1,NOMEGA
        OMEGA=REAL(2*NU-1,KIND=8)*PI*KBT
        SVAR=0.D0
        DO I=1,NORB
          S(I,I,NU)=0.007D0*EXP(-0.1D0/ABS(OMEGA))/(CI*OMEGA)
!         == FACTOR 2 FROM -OMEGA_NU ===========================================
          SVAR=SVAR+2.D0*REAL(S(I,I,NU)*G(I,I,NU),KIND=8)
        ENDDO
        ETOT=ETOT+KBT*SVAR
      ENDDO      
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_CONSTRAINTS()
!     **************************************************************************
!     **  ADJUSTS GAMMA(:,:,IKPT,IPSIN) SUCH THAT                             **
!     **                                                                      **
!     **  G(K,NU)=<PI(K)|PSI(K)> &                                            **
!     **         /[I*OMEGA_NU+MU-<PSI(K)|PI>(H0(K)-SIGMA(NU))<PI(K)|PSI(K)>] &**
!     **         *<PSI(K)|PI(K)>                                              **
!     **                                                                      **
!     **  PRODUCES THE SAME K-DEPENDENT DENSITY MATRIX AS GRHO                **
!     **                                                                      **
!     **  GRHO(K,NU)=<PI(K)|PSI(K)>/[I*OMEGA_NU+MU-HRHO(K)]<PI(K)|PSI(K)>     **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NDIMD,NAT,NOMEGA,NLAU &
     &                      ,OMEGA,KBT,MU,AMIX,KSET,ATOMSET
      IMPLICIT NONE
      REAL(8)   ,PARAMETER     :: TOL=1.D-5
      INTEGER(4),PARAMETER     :: NITER=100000
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MATX(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: DMATX(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: G(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR(NCHI,NCHI,NDIMD,NLAU+1)
      COMPLEX(8)               :: SLAUR(NCHI,NCHI,NDIMD,NLAU)
      COMPLEX(8)               :: DRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: DH0(NCHI,NCHI,NDIMD)
      LOGICAL(4)               :: CONVG
      REAL(8)                  :: MAXDEV
      INTEGER(4)               :: IKPT,ITER,IDIMD,NU,IAT
      INTEGER(4)               :: I1,I2,ILAU
      REAL(8)                  :: FN(2)
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_CONSTRAINTS')
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,2,FN)
!
      IF(NLAU.GT.2) THEN
        CALL ERROR$MSG('NLAU EXCEEDS MAXIMUM PERMITTED VALUE OF 2')
        CALL ERROR$I4VAL('NLAU',NLAU)
        CALL ERROR$STOP('DMFT_CONSTRAINTS')
      END IF
      IF(NLAU.LT.1) THEN
        CALL ERROR$MSG('NLAU MUST BE AT LEAST EQUAL TO ONE')
        CALL ERROR$I4VAL('NLAU',NLAU)
        CALL ERROR$STOP('DMFT_CONSTRAINTS')
      END IF
!
!     ==========================================================================
!     == CONSTRUCT LAURENT EXPANSION OF SELF-ENERGY, WHICH IS K-INDEPENDENT   ==
!     ==========================================================================
      SLAUR=(0.D0,0.D0)
      DO IAT=1,NAT
        I1=ATOMSET(IAT)%ICHI1
        I2=ATOMSET(IAT)%ICHI2
        DO ILAU=1,NLAU
          SLAUR(I1:I2,I1:I2,:,ILAU)=ATOMSET(IAT)%SLOCLAUR(:,:,:,ILAU)
        ENDDO
      ENDDO
!       
!     ==========================================================================
!     ==  DEVIATION FROM TARGET DENSITY MATRIX                                ==
!     ==========================================================================
      DO IKPT=1,NKPTL
!       
!       ========================================================================
!       ==  ADJUST H0 TO OBEY DENSITY MATRIX CONSTRAINT                       ==
!       ========================================================================
        DO ITER=1,NITER
!         == LAURENT EXPANSION FOR THE GREENS FUNCTION DIFFERENCE ==============
          GLAUR=(0.D0,0.D0)
!         == MATX=-MU*S+HRHO ; DMATX=SLAUR(0)-GAMMA ============================
          MATX=-MU*KSET(IKPT)%SMAT+KSET(IKPT)%HRHO
          DMATX=SLAUR(:,:,:,1)-KSET(IKPT)%GAMMA 
!         == GLAUR2=SINV*MATX*SINV =============================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,DMATX,KSET(IKPT)%SINV,MAT)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,MAT,GLAUR(:,:,:,2))
!         == THE NEXT TERM GLAUR3 DOES NOT CONTRIBUTE TO THE DENSITY MATRIX. ===
!         == GLAUR4 IS NOT IMPLEMENTED BECAUSE NLAU<3 ==========================
!
!         == LOOP OVER OMEGA ===================================================
          DRHO=(0.D0,0.D0)
          DO NU=1,NOMEGA
!           == CONSTRUCT LATTICE GREENS FUNCTION =============================
            MAT=(CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT-KSET(IKPT)%HRHO 
            CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
            DRHO=DRHO-KBT*G  ! SUBTRACT NONINTERACTING G
            MAT=MAT+KSET(IKPT)%GAMMA
            DO IAT=1,NAT
              I1=ATOMSET(IAT)%ICHI1
              I2=ATOMSET(IAT)%ICHI2
              MAT(I1:I2,I1:I2,:)=MAT(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
            ENDDO
            CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
            DRHO=DRHO+KBT*G
          ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!         == INCLUDE NEGATIVE FREQUENCIES ======================================
          DO IDIMD=1,NDIMD
            DRHO(:,:,IDIMD)=DRHO(:,:,IDIMD)+CONJG(TRANSPOSE(DRHO(:,:,IDIMD)))
          ENDDO
!         == ADD TAILS (GLAUR3 DOES NOT CONTRIBUTE) ============================
          DO ILAU=1,NLAU+1
            DRHO=DRHO+FN(ILAU)*GLAUR(:,:,:,ILAU)
          ENDDO
!
!         ======================================================================
!         == ADJUST HAMILTONIAN (LAGRANGE MULTIPLIER)                         ==
!         ======================================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,DRHO,KSET(IKPT)%SMAT,MAT)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SMAT,MAT,DH0)
          KSET(IKPT)%GAMMA=KSET(IKPT)%GAMMA-4.D0*KBT*DH0
!
!         ======================================================================
!         == CHECK CONVERGENCE                                                ==
!         ======================================================================
          IF(MOD(ITER,100).EQ.0)PRINT*,'MAXVAL OF DH /DRHO ',ITER &
      &                        ,4.D0*KBT*MAXVAL(ABS(DH0)),MAXVAL(ABS(DRHO))
          MAXDEV=MAXVAL(ABS(DRHO))*REAL(NKPTL,KIND=8)
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
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_MIX()
!     **************************************************************************
!     **  MIXES DPHIDG INTO THE SELF ENERGY                                   **
!     **  DPHIDG IS THE OUTPUT SELF ENERGY OF THE SOLVER                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NAT,NOMEGA,OMEGA,ATOMSET
      IMPLICIT NONE
      REAL(8)  ,PARAMETER :: MIX1=1.D-2
      REAL(8)  ,PARAMETER :: MIX2=1.D0
      INTEGER  ,PARAMETER :: NUH=50
      REAL(8)             :: SVAR
      INTEGER(4)          :: NLOC
      INTEGER(4)          :: IAT,NU
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_MIX')
      IF(NUH.GT.NOMEGA) THEN
        CALL ERROR$MSG('NUH MUST NOT BE LARGER THAN NOMEGA')
        CALL ERROR$STOP('DMFT_MIX')
      END IF
      DO IAT=1,NAT
        NLOC=ATOMSET(IAT)%NLOC
        IF(NLOC.LE.0) CYCLE
!
!       == DIFFERNCE OUT-IN, WHICH WILL BE SCALED AND ADDED TO CURRENT VALUES ==
        ATOMSET(IAT)%DPHIDG    =ATOMSET(IAT)%DPHIDG    -ATOMSET(IAT)%SLOC
        ATOMSET(IAT)%DPHIDGLAUR=ATOMSET(IAT)%DPHIDGLAUR-ATOMSET(IAT)%SLOCLAUR
!
!       == MIX FREQUENCY DEPENDENT SELF ENERGY =================================
        DO NU=1,NOMEGA
!          SVAR=1.D0/(1.D0-MIX1/(MIX1-1.D0)+MIX2*OMEGA(NU)**2)
          SVAR=1.D0/(1.D0+(1.D0-MIX1)/(MIX1+(OMEGA(NU)/OMEGA(NUH))**2))
          ATOMSET(IAT)%SLOC(:,:,:,NU)=ATOMSET(IAT)%SLOC(:,:,:,NU) &
     &                          +SVAR*ATOMSET(IAT)%DPHIDG(:,:,:,NU)
        ENDDO
!
!       == MIX IN WITH FREQUENCY INDEPENDENT MASS ==============================
        SVAR=1.D0
        ATOMSET(IAT)%SLOCLAUR=ATOMSET(IAT)%SLOCLAUR+SVAR*ATOMSET(IAT)%DPHIDGLAUR
!
        ATOMSET(IAT)%DPHIDG=(0.D0,0.D0)
        ATOMSET(IAT)%DPHIDGLAUR=(0.D0,0.D0)
      ENDDO
                               CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_STATICSOLVER(ETOT)
!     **************************************************************************
!     **  ENERGY AND SELF ENERGY FROM HARTREE FOCK AND DFT DOUBLE COUNTING    **
!     **  PRODUCES ETOT                                                       **
!     **           ATOMSET(IAT)%DENMAT%H                                      **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE ,ONLY: NAT,NDIMD,ATOMSET
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: ETOT
      LOGICAL(4),PARAMETER   :: STATICOFF=.FALSE.
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      REAL(8)                :: PHILW
      REAL(8)                :: LHFWEIGHT
      REAL(8)                :: EDC
      INTEGER(4)             :: NLOC  !#(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4)             :: LMNX  !#(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4)             :: IAT
      INTEGER(4)             :: IDIMD,I,NU
      COMPLEX(8),ALLOCATABLE :: HAM(:,:,:)
!     **************************************************************************
                                      CALL TRACE$PUSH('DMFT_STATICSOLVER')
      IF(STATICOFF) THEN
        ETOT=0.D0
        DO IAT=1,NAT
          ATOMSET(IAT)%DENMAT%H=(0.D0,0.D0)
        ENDDO
                                      CALL TRACE$POP()
        RETURN
      END IF
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
!       == DETERMINE HARTREE-FOCK CONTRIBUTION                                ==
!       ========================================================================
        CALL DMFT_FOCK(NLOC,NDIMD,ATOMSET(IAT)%DENMAT%RHO,ATOMSET(IAT)%U &
     &                                             ,PHILW,ATOMSET(IAT)%DENMAT%H)
        ETOT=ETOT+PHILW  !SCREENING DONE ALREADY IN U-TENSOR
        WRITE(*,FMT='(50("."),T1,"HARTREE-FOCK ENERGY FOR ATOM ",I3,T50,":"' &
       &             //',F20.10)')IAT,PHILW
!
!       == PRINT IF DESIRED ====================================================
        IF(TPRINT) THEN
          WRITE(*,FMT='(82("="),T10," ",A," IAT=",I4,"  ")') &
       &                        'FOCK HAMILTONIAN',IAT
          DO IDIMD=1,NDIMD
            DO I=1,NLOC
              WRITE(*,FMT='(4I5,100F10.5)')IAT,NU,IDIMD,I &
      &                                   ,ATOMSET(IAT)%DENMAT%H(I,:,IDIMD)
            ENDDO
          ENDDO
        END IF
!
!       ========================================================================
!       == SUBTRACT DFT DOUBLE COUNTING TERM                                  ==
!       ========================================================================
!       == COLLECT LOCAL HF WEIGHT =============================================
        LHFWEIGHT=ATOMSET(IAT)%LHFWEIGHT
        ALLOCATE(HAM(NLOC,NLOC,NDIMD))
        CALL DMFT_DC(IAT,NLOC,NDIMD,ATOMSET(IAT)%DENMAT%RHO,EDC,HAM)
        ETOT=ETOT+LHFWEIGHT*EDC
        ATOMSET(IAT)%DENMAT%H=ATOMSET(IAT)%DENMAT%H+LHFWEIGHT*HAM
        WRITE(*,FMT='(50("."),T1,"DOUBLE COUNTING FOR ATOM ",I3,T50,":"' &
       &             //',F20.10)')IAT,LHFWEIGHT*EDC
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
      ENDDO ! END OF LOOP OVER ATOMS (IAT)
                                      CALL TRACE$POP()
      RETURN
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
      COMPLEX(8),INTENT(IN)  :: RHO(NLOC,NLOC,NDIMD)
      REAL(8)   ,INTENT(IN)  :: UCHI(NLOC,NLOC,NLOC,NLOC)
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
      SUBROUTINE DMFT_DC(IAT,LMNX,NDIMD,RHO,EDC,HAM)
!     **************************************************************************
!     ** CALCULATES THE DFT EXCHANGE ENERGY FOR THE CORRELATED ORBITALS.      **
!     ** EDC IS TO BE ADDED! TO THE TOTAL ENERGY (MINUS SIGN IS INCLUDED)     **
!     **                                                                      **
!     ** REMARK: SOME ROUTINES ONLY WORK IN THE NON-COLLINEAR MODE            **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR=>POTPAR1,SBAR=>SBAR_NEW
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      INTEGER(4),INTENT(IN)  :: LMNX          !  #(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4),INTENT(IN)  :: NDIMD
      COMPLEX(8),INTENT(IN)  :: RHO(LMNX,LMNX,NDIMD) ! DENSITY MATRIX
      REAL(8)   ,INTENT(OUT) :: EDC           ! -E_(DFT-EXCHANGE)
      COMPLEX(8),INTENT(OUT) :: HAM(LMNX,LMNX,NDIMD) ! HAMILTONIAN CONTRIB.
      REAL(8)   ,ALLOCATABLE :: D(:,:,:)     !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: DT(:,:,:)    !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: DTALL(:,:,:) !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: HT(:,:,:)    !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: HTALL(:,:,:) !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: H(:,:,:)  !(LMNX,LMNX,NDIMD)
      REAL(8)   ,ALLOCATABLE :: HALL(:,:,:)  !(LMNX,LMNX,NDIMD)
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
      INTEGER(4)             :: NNS
      INTEGER(4)             :: INS
      INTEGER(4)             :: NN
!     **************************************************************************
      EDC=0.D0
      HAM=0.D0
      CALL DFT$GETI4('TYPE',IDFTTYPE)
!      PRINT*,'IDFTTYPE ',IDFTTYPE
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
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
        IF(SBAR(NN)%IAT2.NE.IAT) CYCLE
        IF(SUM(ABS(SBAR(NN)%IT)).NE.0) CYCLE
        INS=NN
        EXIT
      ENDDO
!
!     ==========================================================================
!     ==  MAP ONTO REAL ARRAY. (T,X,Y,Z) REPRESENTATION                       ==
!     ==========================================================================
!  CAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTION
!
!  THE DOUBLE COUNTING TERM IS NOT CALCULATED CORRECTLY. DTALL SHALL DESCRIBE 
!  THE FULL DENSITY INCLUDING PARTIAL WAVES THAT HAVE NO CORRESPONDING LOCAL 
!  ORBITALS. CHECK WITH PBLOECHL C611D05
!
!  CAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTIONCAUTION

      ALLOCATE(D(LMNX,LMNX,NDIMD))
      ALLOCATE(H(LMNX,LMNX,NDIMD))
      ALLOCATE(HALL(LMNX,LMNX,NDIMD))
      ALLOCATE(DT(LMNXT,LMNXT,NDIMD))
      ALLOCATE(DTALL(LMNXT,LMNXT,NDIMD))
      ALLOCATE(HT(LMNXT,LMNXT,NDIMD))
      ALLOCATE(HTALL(LMNXT,LMNXT,NDIMD))
      D=REAL(RHO)
      CALL LMTO_EXPANDLOCAL('FWRD',NDIMD,LMNX,LMNXT,SBAR(INS)%MAT,D,DT)
      DTALL=DT
      CALL LMTO_SIMPLEDC_NEW(GID,NR,NDIMD,LMNXT,LNXT,LOXT &
     &                  ,POTPAR(ISP)%TAILED%AEF &
     &                  ,LRX,AECORE,DT,DTALL,EX,HT,HTALL)
      CALL LMTO_EXPANDLOCAL('BACK',NDIMD,LMNX,LMNXT,SBAR(INS)%MAT,HALL,HTALL)
      CALL LMTO_EXPANDLOCAL('BACK',NDIMD,LMNX,LMNXT,SBAR(INS)%MAT,H,HT)
      EDC=-EX
      HAM=-CMPLX(H+HALL,KIND=8)
      DEALLOCATE(D)
      DEALLOCATE(DT)
      DEALLOCATE(DTALL)
      DEALLOCATE(HTALL)
      DEALLOCATE(HT)
      DEALLOCATE(H)
      DEALLOCATE(HALL)
!PRINT*,'DOUBLE COUNTING CORRECTION ENERGY FOR ATOM=',IAT,-EX
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
!         == MAT1=1/GRHO =======================================================
          MAT1=(CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT-KSET(IKPT)%HRHO
!         == CALCULATE GRHO ====================================================
          CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,GRHO)
!     
!         == MAT2=SIGMA-GAMMA (=HPRIME+SIGMA-HRHO) =============================
          MAT2=-KSET(IKPT)%GAMMA
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
          CALL SPINOR$MATMUL(NDIMD,NCHI,GRHO-GFULL,-KSET(IKPT)%GAMMA,MAT1)
          XSUM=0.D0
          DO I=1,NCHI
            XSUM=XSUM+2.D0*REAL(MAT1(I,I,1))  ! FACTOR TWO COMES FROM (-OMEGA)
          ENDDO
          ETOT1=ETOT1+XSUM
!     
!         ======================================================================
!         == (SIGMA-GAMMA)*G (=(HPRIME+SIGMA-HRHO)*G)                         ==
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
      WRITE(*,FMT='(50("."),T1,"ENERGY FROM MATSUBARA SUM FOR ATOM ",I3' &
     &                  //',T50,":",F20.10)')IAT,ETOT
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TESTG(NORB,NOMEGA,NLAU,KBT,G,GLAUR)
!     **************************************************************************
!     ** PRINTS THE LOCAL DENSITY MATRIX FROM THE LOCAL GREENS FUNCTION       **
!     ** FOR TESTING PURPOSES ONLY                                            **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB     !#(SPIN ORBITALS)
      INTEGER(4),INTENT(IN) :: NOMEGA   !#(POSITIVE MATSUBARA FREQUENCIES
      INTEGER(4),INTENT(IN) :: NLAU   !#(LAUENT EXPANSION TERMS SELF ENERGY)
      REAL(8)   ,INTENT(IN) :: KBT      
      COMPLEX(8),INTENT(IN) :: G(NORB,NORB,NOMEGA)     !GREENS FUNCTION
      COMPLEX(8),INTENT(IN) :: GLAUR(NORB,NORB,NLAU+1) !LAURENT EXPANSION OF G
      REAL(8)               :: PI
      REAL(8)               :: SVAR
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
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,NLAU+1,FN)
      DO I=1,NLAU+1
        RHO=RHO+GLAUR(:,:,I)*FN(I)
      ENDDO
!
      WRITE(*,FMT='(100("="),T30,"  ",A,"  ")')'DENSITY MATRIX IN TESTG'
      DO I=1,NORB
        WRITE(*,FMT='(100("(",2F10.5,")"))')RHO(I,:)
      ENDDO
      SVAR=0.D0
      DO I=1,NORB
        SVAR=SVAR+REAL(RHO(I,I))
      ENDDO
      WRITE(*,FMT='("NUMBER OF ELECTRONS IN G: ",F10.5)')SVAR
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_ADDTOHPSI()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY  : TON,NKPTL,NSPIN,NDIM,NDIMD,NB,NCHI,NAT &
     &                        ,ATOMSET,KSET
      USE WAVES_MODULE, ONLY : GSET,WAVES_SELECTWV,THIS,MAP
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      COMPLEX(8),ALLOCATABLE :: DHPIPSI(:,:,:) !(NDIM,NCHI,NB)
      COMPLEX(8),ALLOCATABLE :: MAT(:,:,:)     !(NCHI,NCHI,NDIMD)
      LOGICAL(4)             :: TCHK,TRESET
      INTEGER(4)             :: NBH
      INTEGER(4)             :: NLOC
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
          NLOC=ATOMSET(IAT)%NLOC
          IF(NLOC.LE.0) CYCLE ! NO DOUBLE COUNTING          
          PRINT*,'ATOM ',IAT,' NLOC=',NLOC
          CALL SPINOR_PRINTMATRIX(NFIL,'DC HAMILTONIAN(TXYZ)',1,NLOC &
    &                            ,NDIMD,NLOC,ATOMSET(IAT)%DENMAT%H)
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
        NLOC=ATOMSET(IAT)%NLOC
        IF(NLOC.LE.0) CYCLE ! NO DOUBLE COUNTING          
!CALL SPINOR_PRINTMATRIX(6,'ATOMSET%H',1,NLOC,NDIMD,NLOC,ATOMSET(IAT)%DENMAT%H)
        CALL SPINOR$CONVERT('BACK',NLOC,NDIMD,ATOMSET(IAT)%DENMAT%H)
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
            IF(.NOT.ASSOCIATED(THIS%HTBC_NEW)) &
      &               ALLOCATE(THIS%HTBC_NEW(NDIM,NBH,NPRO))
            THIS%HTBC_NEW=(0.D0,0.D0)
          END IF
!
!         ======================================================================
!         == DF/DRHO=H0-HRHO; DHDPIPSI= (H0-HRHO)<PI|PSI>                     ==
!         ======================================================================
          MAT=KSET(IKPT)%GAMMA(:,:,:)
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
!           == PIPSI(ICHI,IBH,IKPT,ISPIN)=THIS%TBC(1,IBH,ICHI) =================
            THIS%HTBC_NEW(:,:,ICHI)=THIS%HTBC_NEW(:,:,ICHI)+DHPIPSI(:,ICHI,:NBH)
          ENDDO
!
!         ======================================================================
!         ==  NOW THE SITE-LOCAL TERM FROM THE DOUBLE COUNTING                ==
!         ======================================================================
          IPRO=0
          DO IAT=1,NAT
            NLOC=ATOMSET(IAT)%NLOC
            I1=IPRO+1
            I2=IPRO+NLOC
            IPRO=IPRO+NLOC
            IF(NLOC.LE.0) CYCLE ! NO DOUBLE COUNTING          
!
            DO IDIM2=1,NDIM
              DO IDIM1=1,NDIM
                IDIMD=IDIM1+NDIM*(IDIM2-1)+ISPIN-1
!               == HTBC=HAM*TBC ================================================
!               == HTBC(IB,I)=HTBC(IB,I)+HAM(I,J)*TBC(IB,J) ====================
                THIS%HTBC_NEW(IDIM1,:,I1:I2)=THIS%HTBC_NEW(IDIM1,:,I1:I2) &
         &                +MATMUL(THIS%TBC_NEW(IDIM1,:,I1:I2) &
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
            A(I,J,IDIMD)=CMPLX(RAN1,RAN2,KIND=8)
          ENDDO
        ENDDO
      ENDDO
      DO IDIMD=1,NDIMD
        DO I=1,NCHI
          DO J=1,NCHI
            CALL RANDOM_NUMBER(RAN1)
            CALL RANDOM_NUMBER(RAN2)
            B(I,J,IDIMD)=CMPLX(RAN1,RAN2,KIND=8)
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
!     **  DERIVATION IN METHODS: SECTION: WORKING WITH SPIN ORBITALS,         **
!     **  SUBSECTION: INVERSION OF A MATRIX IN SPINOR REPRESENTATION          **
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
               MAT(I,I,IDIMD)=MAT(I,I,IDIMD)-(1.D0,0.D0)
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_EXPLOREMODULE()
!     **************************************************************************
!     ** TESTR ROUTINE TO EXPLORE THE CONTENT OF THE DMFT MODULE              **
!     ** THIS ROUTINE SHALL NOT CHANGE THE STATE OF THE OBJECT.               **
!     **************************************************************************
      USE DMFT_MODULE, ONLY :NAT,NKPTL,NDIMD,NCHI,KSET,ATOMSET
      INTEGER(4)  :: NLOC
      INTEGER(4)  :: IKPT,ISPIN,IAT
      INTEGER(4)  :: I,J,K,L
!     **************************************************************************
      DO IKPT=1,NKPTL
        DO ISPIN=1,NDIMD
          DO I=1,NCHI
            DO J=1,NCHI
              WRITE(*,FMT='("IK=",I4," IS=",I2," IJ=",2I4," RHO=",2F20.10)') &
     &                            IKPT,ISPIN,I,J,KSET(IKPT)%RHO(I,J,ISPIN)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO IAT=1,NAT
        NLOC=ATOMSET(IAT)%NLOC
        DO I=1,NLOC
          DO J=1,NLOC
            DO K=1,NLOC
              DO L=1,NLOC
                WRITE(*,FMT='("IAT=",I4," IJKL=",4I4," U=",F20.10)') &
     &                   IAT,I,J,K,L,ATOMSET(IAT)%U(I,J,K,L)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO     
      RETURN
      END
