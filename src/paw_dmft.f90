! FUDGES:
! REGMATSUBARA: THE 2ND LAURENT TERM HAS BEEN CHANGED TO AVOID  
!               A NON-MONOTONEOUS FERMI FUNCTION
!
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
!
TYPE NATORB_TYPE                             ! NATURAL ORBITALS
  COMPLEX(8),POINTER  :: CHIPHI(:,:) => NULL() !(2*NLOC,2*NLOC) <CHI|NATORB>
  COMPLEX(8),POINTER  :: PIPHI(:,:)  => NULL() !(2*NLOC,2*NLOC) <PI|NATORB>
END TYPE NATORB_TYPE
!
TYPE ATOMSET_TYPE
  !**  nloc       =#(local orbotals on this atom)
  !**  u          = U-tensor
  !**  gni        = non-interacting local Greens function
  !**  dphidg     = self energy from luttinger-ward functional
  !**  rho        = local density matrix 
  INTEGER(4)           :: NLOC
  INTEGER(4)           :: ICHI1
  INTEGER(4)           :: ICHI2
  REAL(8)              :: LHFWEIGHT
  REAL(8)   ,POINTER   :: U(:,:,:,:)        => NULL() !(NLOC,NLOC,NLOC,NLOC) 
  REAL(8)   ,POINTER   :: DEDU(:,:,:,:)     => NULL() !(NLOC,NLOC,NLOC,NLOC) 
  COMPLEX(8),POINTER   :: GLOC(:,:,:,:)     => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: GLOCLAUR(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,NLAU+1)
  COMPLEX(8),POINTER   :: GNI(:,:,:,:)      => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: GNILAUR(:,:,:,:)  => NULL() !(NLOC,NLOC,NDIMD,NLAU+1)
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
  !** pipsi   = <pi|psi> projection of band state onto local orbital
  !** rho     = density matrix
  !** rho     = density matrix made n-representable
  !** gamma   = lagrange multipliers for density matrix constraint
  REAL(8)            :: WKPT
  LOGICAL(4)         :: TADDMINUSK    !ADD THE TERM FOR -K
  REAL(8)   ,POINTER :: F(:,:)        !(NB,NSPIN)  occupation       
  REAL(8)   ,POINTER :: E(:,:)        !(NB,NSPIN)  energy       
  COMPLEX(8),POINTER :: PIPSI(:,:,:,:)!(NDIM,NCHI,NB,NSPIN) 
  COMPLEX(8),POINTER :: RHO(:,:,:)    !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: RHOADAPTED(:,:,:) !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: GAMMA(:,:,:)  !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: HRHO(:,:,:)   !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: SINV(:,:,:)   !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: SMAT(:,:,:)   !(NCHI,NCHI,NDIMD)
  COMPLEX(8),ALLOCATABLE :: CHIPHI(:,:,:) !(NCHI,NCHI,NDIMD)
  COMPLEX(8),ALLOCATABLE :: PIPHI(:,:,:)  !(NCHI,NCHI,NDIMD)
END TYPE KSET_TYPE
LOGICAL(4),PARAMETER   :: TON=.TRUE.
LOGICAL(4),PARAMETER   :: TONATORB=.FALSE.
LOGICAL(4),SAVE        :: TINI=.FALSE.
logical(4),PARAMETER   :: tsimpledmft=.true. ! no density matrix constraint
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
      NOMEGA=300 ! #(POSITIVE MATSUBARA FREQUENCIES)
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
        ALLOCATE(ATOMSET(IAT)%GNI(NLOC,NLOC,NDIMD,NOMEGA))
        ALLOCATE(ATOMSET(IAT)%GNILAUR(NLOC,NLOC,NDIMD,NLAU+1))
        ALLOCATE(ATOMSET(IAT)%DPHIDG(NLOC,NLOC,NDIMD,NOMEGA))
        ALLOCATE(ATOMSET(IAT)%DPHIDGLAUR(NLOC,NLOC,NDIMD,NLAU))
        ALLOCATE(ATOMSET(IAT)%SLOC(NLOC,NLOC,NDIMD,NOMEGA))
        ALLOCATE(ATOMSET(IAT)%SLOCLAUR(NLOC,NLOC,NDIMD,NLAU))
        ALLOCATE(ATOMSET(IAT)%SMAT(NLOC,NLOC,NDIMD))
        ALLOCATE(ATOMSET(IAT)%DENMAT%RHO(NLOC,NLOC,NDIMD))
        ALLOCATE(ATOMSET(IAT)%DENMAT%H(NLOC,NLOC,NDIMD))
        ATOMSET(IAT)%U=0.D0
        ATOMSET(IAT)%DEDU=0.D0
        ATOMSET(IAT)%GLOC=(0.D0,0.D0)
        ATOMSET(IAT)%GLOCLAUR=(0.D0,0.D0)
        ATOMSET(IAT)%GNI=(0.D0,0.D0)
        ATOMSET(IAT)%GNILAUR=(0.D0,0.D0)
        ATOMSET(IAT)%DPHIDG=(0.D0,0.D0)
        ATOMSET(IAT)%DPHIDGLAUR=(0.D0,0.D0)
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
        ALLOCATE(KSET(IKPTL)%E(NB,NSPIN))
        ALLOCATE(KSET(IKPTL)%PIPSI(NDIM,NCHI,NB,NSPIN))
        ALLOCATE(KSET(IKPTL)%RHO(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%RHOADAPTED(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%GAMMA(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%HRHO(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%SINV(NCHI,NCHI,NDIMD))
        ALLOCATE(KSET(IKPTL)%SMAT(NCHI,NCHI,NDIMD))
        KSET(IKPTL)%WKPT=0.D0
        KSET(IKPTL)%F=0.D0
        KSET(IKPTL)%E=0.D0
        KSET(IKPTL)%PIPSI=(0.D0,0.D0)
        KSET(IKPTL)%RHO=(0.D0,0.D0)
        KSET(IKPTL)%RHOADAPTED=(0.D0,0.D0)
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
      LOGICAL(4)               :: TADJUST=.FALSE.
      REAL(8)                  :: FACTOR(20)
      REAL(8)                  :: SVAR
      REAL(8)                  :: PI
      INTEGER(4)               :: J,NU
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
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
!
!     ==========================================================================
!     == ADJUST REGULARIZATION TO AVOID NON-MONOTONEOUS FERMI FUNCTION        ==
!     ==========================================================================
      IF(TADJUST) THEN
        IF(N.GE.2) THEN
          FN(2)=-1.D0/(2.D0*PI**2*REAL(NOMEGA,KIND=8))
          FN(2)=FN(2)-1.D-2
        END IF
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
      USE DMFT_MODULE, ONLY: TON &
     &                      ,tsimpledmft
      USE MPE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER   :: NITERX=500
      REAL(8)   ,PARAMETER   :: ETOL=1.D-6
      REAL(8)   ,PARAMETER   :: SIGMATOL=1.D-5
      REAL(8)                :: XDEV ! X(DEVIATION OF THE SELF ENERGY)
      REAL(8)                :: ETOT
      REAL(8)                :: ELAST
      REAL(8)                :: SVAR
      INTEGER(4)             :: ITER
      LOGICAL(4) :: TREPEAT=.false. ! USED FOR GRADIENT TEST
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
      CALL DMFT_DOSHIST() ! DOS HISTROGRAM FROM 1P-ENERGIES 
      CALL DMFT_GNI()     ! NON INTERACTING GREENS FUNCTION (GBAR)
!
!     == CAUTION: STOREKSET SETS A LOT OF DATA TO ZERO ON OUTPUT
1000  IF(TREPEAT) CALL DMFT_STOREKSET() ! STORE AND RESET COMPLETE INPUT 
!
!     ==========================================================================
!     ==  DENSITY MATRIX AND OVERLAP MATRIX IN K-SPACE                        ==
!     ==========================================================================
      CALL DMFT_RHOOFK()
!
!     ==========================================================================
!     ==  OVERLAP MATRIX IN K-SPACE
!     ==========================================================================
!This apparently overwrites sinv constructed in dmft_rhoofk
!      CALL DMFT_OVERLAPOFK()
!
      IF(TREPEAT)CALL DMFT_CHECKGRADIENT1()
!
!     ==========================================================================
!     ==  OBTAIN BARE U-TENSOR FROM LMTO OBJECT.                              ==
!     ==     SHAPE OF ORBITALS ARE DEFINED BY NPRO AND ATOMIC STRUCTRE        ==
!     ==     ORBITALS ARE SELECTED BY TORB.                                   ==
!     ==========================================================================
      CALL DMFT_UTENSOR() 
!
!     ==========================================================================
!     ==  CONSTRUCT NON-INTERACTING HAMILTONIAN THAT PRODUCES THE CORRECT     ==
!     ==  ONE-PARTICLE DENSITY MATRIX                                         ==
!     ==========================================================================
!CALL DMFT_EXPLOREMODULE()
      CALL DMFT_HRHO()
!      CALL DMFT_TONATORB('FWRD')

!
!     ==========================================================================
!     ==  CONSTRUCT LOCAL NATURAL ORBITALS                                    ==
!     ==========================================================================
      CALL DMFT_NATORB()
!      CALL DMFT_NININTSPECTALDENSITY() ! USED FOR TESTING
!      CALL DMFT_GLOCNI()
!      CALL DMFT_LOCALSPECTRALFUNCTION()
!
!     ==========================================================================
!     == ITERATION TO ENFORCE CONSTRAINTS                                     ==
!     ==========================================================================
      ELAST=0.D0
      DO ITER=1,NITERX
        WRITE(*,FMT='(82("="),T20," ITERATION ",I5)')ITER
        CALL DMFT_GLOC() 
!
print*,'before printing spectral function iter =',iter
call DMFT_fullLOCALSPECTRALFUNCTION()
print*,'after printing spectral function'

!CALL DMFT_LOCALSPECTRALFUNCTION()
!
!       ========================================================================
!       ==  CALL THE SOLVER                                                   ==
!       ========================================================================
        CALL DMFT_SOLVER(ETOT) 
        WRITE(*,FMT='(60("."),T1'&
     &        //',"--DMFT--ENERGY AFTER DMFT_SOLVER"' &
     &        //',T60,":",F20.10)')ETOT
!
!       ========================================================================
!       ==  MIXING                                                            ==
!       ========================================================================
        CALL DMFT_MIX(XDEV)
!
!       ========================================================================
!       ==  CONSTRAINTS                                                       ==
!       ========================================================================
        if(.not.tsimpledmft) then
          CALL DMFT_CONSTRAINTS()
        end if
print*,'constraints done ',tsimpledmft
        PRINT*,'ITERATION COMPLETED ',ITER,XDEV,ETOT,ETOT-ELAST
        IF(ABS(ETOT-ELAST).LT.ETOL.AND.XDEV.LT.SIGMATOL) EXIT
        ELAST=ETOT
      ENDDO ! END OF LOOP OVER ITERATIONS TO ENFORCE CONSTRAINT
!
!     ==========================================================================
!     ==  ADD HARTREE FOCK AND DFT DOUBLE COUNTING                            ==
!     ==========================================================================
!     == SWITCH STATICSOLVER OFF BECAUSE THIS DOES NOT ENTER GAMMA      
      CALL DMFT_RHOLOCAL() ! LOCAL DENSITY MATRIX USED IN STATICSOLVER
      CALL DMFT_STATICSOLVER(SVAR)
      IF(.NOT.TREPEAT) ETOT=ETOT+SVAR
      WRITE(*,FMT='(60("."),T1'&
     &        //',"--DMFT--ENERGY AFTER DMFT_STATICSOLVER"' &
     &        //',T60,":",F20.10)')ETOT
! 
!     ==========================================================================
!     ==  CALCULATE MISSING TOTAL ENERGY CONTRIBUTION                         ==
!     ==========================================================================
      CALL DMFT_DETOT(SVAR)
      ETOT=ETOT+SVAR      
      WRITE(*,FMT='(60("."),T1'&
     &        //',"--DMFT--ENERGY AFTER DMFT_DETOT"' &
     &        //',T60,":",F20.10)')ETOT
      CALL ENERGYLIST$SET('DMFT INTERFACE',ETOT)
      CALL ENERGYLIST$ADD('LOCAL CORRELATION',ETOT)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',ETOT)
!
!     ==========================================================================
!     ==  ADD HAMILTONIAN THIS%HTBC                                           ==
!     ==========================================================================
      IF(TREPEAT) THEN
        CALL DMFT_CHECKGRADIENT2(ETOT)
        GOTO 1000
      END IF
IF(TREPEAT) STOP 'FORCED BECAUSE TREPEAT'
!
!      CALL DMFT_TONATORB('BACK')
      CALL DMFT_ADDTOHPSI()
WRITE(*,FMT='(82("="),T20," LEAVING DMFT$GREEN ")')

call error$stop('forced in dmft$green')
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
     &                      ,KSET &
&  ,KBT
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IB,J
      REAL(8)                :: F(NB,NKPTL,NSPIN)
      REAL(8)                :: F1(NB,NKPTL,NSPIN)
      REAL(8)                :: E(NB,NKPTL,NSPIN)
      REAL(8)                :: WKPT(NKPTL)
      INTEGER(4)             :: NB_,NSPIN_
      INTEGER(4)             :: ISPINDEG
      REAL(8)                :: SVAR
      REAL(8)                :: MU1
      REAL(8)                :: EMERMN
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
        CALL ERROR$STOP('DMFT_COLLECTNATORB')
      END IF
!
      CALL DYNOCC$GETI4('NSPIN',NSPIN_)
      IF(NSPIN_.NE.NSPIN) THEN
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$I4VAL('NSPIN_',NSPIN_)
        CALL ERROR$STOP('DMFT_COLLECTNATORB')
      END IF
!
!     ==========================================================================
!     == COLLECT OCCUPATIONS K-POINT WEIGHTS                                  ==
!     ==========================================================================
!     == SPECIAL POINTS HAVE HALF THE WEIGHT OF GENERAL K-POINTS
!     == SUM(WKPT)=1
      CALL WAVES_DYNOCCGETR8A('WKPT',NKPTL,WKPT)
      DO IKPT=1,NKPTL
        KSET(IKPT)%WKPT=WKPT(IKPT)
      ENDDO
!
!     ==========================================================================
!     == COLLECT OCCUPATIONS WITHOUT K-POINT WEIGHT AND SPIN MULTIPLICITY     ==
!     ==========================================================================
!     == OCC FROM WAVES OBJECT CONTAIN SPIN MULTIPLICITY AND K-POINT WEIGHT ====
!     == THE OCCPATIONS STORED OBEY:  0<KSET(IKPT)%F<1 =========================
      CALL WAVES_DYNOCCGETR8A('OCC',NB*NKPTL*NSPIN,F)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          KSET(IKPT)%F(:,ISPIN)=F(:,IKPT,ISPIN)/KSET(IKPT)%WKPT
!         == OCCUPATIONS ARE IN [0,1] ==========================================
          IF(NSPIN.EQ.1) KSET(IKPT)%F(:,ISPIN)=0.5D0*KSET(IKPT)%F(:,ISPIN)
!!$WRITE(*,FMT='(78("="),T10," IKPT=",I5," ISPIN=",I2,"  ")')IKPT,ISPIN
!!$WRITE(*,FMT='("OCC(1) =",10F10.5)')KSET(IKPT)%F(:,ISPIN)
!!$WRITE(*,FMT='("KINDOFE=",10F10.5)') &
!!$      27.211D0*KBT*LOG((1.D0-KSET(IKPT)%F(:,ISPIN))/KSET(IKPT)%F(:,ISPIN))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == COLLECT ENERGY EIGENVALUES
!     ==========================================================================
      CALL WAVES_DYNOCCGETR8A('EPSILON',NB*NKPTL*NSPIN,E)
!     == GET NUMBER OF ELECTRONS TO EXTRACT THE FERMI LEVEL ====================
      SVAR=0.D0
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          SVAR=SVAR+SUM(KSET(IKPT)%F(:,ISPIN))*KSET(IKPT)%WKPT
        ENDDO
      ENDDO
      ISPINDEG=3-NSPIN
      SVAR=SVAR*REAL(ISPINDEG,KIND=8)
PRINT*,'NUMBER OF ELECTRONS=',SVAR
!     == FOR PARALLEL CALCULATION THERE MAY BE AN INDIVIDUAL FERMI LEVEL FOR  ==
!     == EACH  PROCESSOR. THE ENERGY LEVELS ON THIS TASK IS USED CONSISTENTLY ==
!     == WITH THE NUMBER OF ELECTRONS ON THIS TASK.                           ==
PRINT*,'ISPINDEG ',ISPINDEG
PRINT*,'NB=',NB,' NEL=',SVAR,' KBT[EV] ',KBT*27.211D0,' NSPIN ',NSPIN
      CALL DMFT_MERMIN(NB,NKPTL,NSPIN &
     &                 ,SVAR,ISPINDEG,KBT,KSET(:)%WKPT,E,F1,MU1,EMERMN)
PRINT*,'MU FROM MERMIN[EV]',MU1*27.211D0
!     == SHIFT ENERGIES RELATIVE TO THE FERMI-LEVEL ============================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          KSET(IKPT)%E(:,ISPIN)=E(:,IKPT,ISPIN)-MU1
!!$DO IB=1,NB
!!$  IF(ABS(KSET(IKPT)%E(IB,ISPIN)).LT.2.D0/27.211D0) THEN
!!$    WRITE(*,FMT='(3I3,2F10.5)')IKPT,ISPIN,IB,KSET(IKPT)%E(IB,ISPIN)*27.211 &
!!$ &                                          ,E(IB,IKPT,ISPIN)*27.211
!!$  END IF
!!$ENDDO
        ENDDO
      ENDDO
!!$PRINT*,'START HERE ======================'
!!$DO IKPT=1,NKPTL
!!$  DO ISPIN=1,NSPIN
!!$WRITE(*,FMT='(4F10.5)')-40.D0,0.D0,0.D0,0.D0
!!$WRITE(*,FMT='(4F10.5)')-40.D0,1.D0,1.D0,1.D0
!!$    DO IB=1,NB
!!$SVAR=KSET(IKPT)%E(IB,ISPIN)
!!$SVAR=1.D0/(1.D0+EXP(SVAR/KBT))
!!$WRITE(*,FMT='(4F10.5)') &
!!$&    KSET(IKPT)%E(IB,ISPIN)*27.211D0,SVAR,F1(IB,IKPT,ISPIN),KSET(IKPT)%F(IB,ISPIN) 
!!$    ENDDO
!!$WRITE(*,FMT='(4F10.5)')10.D0,0.D0,0.D0,0.D0
!!$  ENDDO
!!$ENDDO
!!$STOP
!
!     ==========================================================================
!     ==  EXTRACT <PI|PSI>                                                    ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          IF(THIS%NB.NE.NB) THEN
            CALL ERROR$MSG('INCONSISTENT NUMBER OF STATES IN WAVES AND DYNOCC')
            CALL ERROR$I4VAL('NB IN DYNOCC',NB)
            CALL ERROR$I4VAL('NB IN WAVES ',THIS%NB)
            CALL ERROR$STOP('DMFT_COLLECTNATORB')
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
      SUBROUTINE DMFT_GNI() 
!     **************************************************************************
!     ** CALCULATE NON-INTERACTING GREENS FUNCTION FOR TESTING PURPOSES       **
!     **                                                                      **
!     ** THIS IS THE MATRIX GRHO WHICH WILL HAVE A DEGENERATE SPECTRAL FNCTN  **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI &
     &                      ,NB &
     &                      ,NKPTL &
     &                      ,NSPIN &
     &                      ,NDIM &
     &                      ,NDIMD &
     &                      ,NAT &
     &                      ,NOMEGA &
     &                      ,OMEGA &
     &                      ,KSET &
     &                      ,ATOMSET &
     &                      ,MU 
      IMPLICIT NONE
      COMPLEX(8),PARAMETER      :: CI=(0.D0,1.D0)
      COMPLEX(8)                :: CSVAR
      INTEGER(4)                :: ICHI1,ICHI2
      INTEGER(4)                :: NLOC
      INTEGER(4)                :: NU
      REAL(8)                   :: WKPTL ! K-POINT WEIGHT
      REAL(8)                   :: E     ! BAND ENERGY
      COMPLEX(8),ALLOCATABLE    :: VEC(:,:)
      COMPLEX(8),ALLOCATABLE    :: MAT(:,:)
      INTEGER(4)                :: IAT,IKPT,ISPIN,IB,I
      INTEGER(4)                :: IDIM1,IDIM2,IDIMD
      INTEGER(4)                :: NFIL
      CHARACTER(64)             :: FILE='FORBACKES.DAT'
!     **************************************************************************
      DO IAT=1,NAT
        ICHI1=ATOMSET(IAT)%ICHI1
        ICHI2=ATOMSET(IAT)%ICHI2
        NLOC=ATOMSET(IAT)%NLOC
        ALLOCATE(VEC(NDIM,NLOC))
        ALLOCATE(MAT(NLOC,NLOC))
        ATOMSET(IAT)%GNI=(0.D0,0.D0)    !(NLOC,NLOC,NDIMD,NOMEGA)
!+KSET(IKPT)%GAMMA
        DO IKPT=1,NKPTL
          WKPTL=KSET(IKPT)%WKPT
          DO ISPIN=1,NSPIN
             DO IB=1,NB
              E=KSET(IKPT)%E(IB,ISPIN)+MU  !KSET%E IS RELATIVE TO EF
              VEC(:,:)=KSET(IKPT)%PIPSI(:,ICHI1:ICHI2,IB,ISPIN)
              DO IDIM1=1,NDIM
                DO IDIM2=1,NDIM
                  IDIMD=IDIM1+NDIM*(IDIM2-1)+ISPIN-1
                  DO I=1,NLOC
                    MAT(I,:)=VEC(IDIM1,I)*CONJG(VEC(IDIM2,:))
                  ENDDO
                  IF(KSET(IKPT)%TADDMINUSK) THEN
                    MAT=MAT+CONJG(TRANSPOSE(MAT(:,:)))
                  END IF
                  MAT=MAT*WKPTL
                  DO NU=1,NOMEGA
                    CSVAR=1.D0/(CI*OMEGA(NU)+MU-E)
                    ATOMSET(IAT)%GNI(:,:,IDIMD,NU) &
    &                             =ATOMSET(IAT)%GNI(:,:,IDIMD,NU)+MAT(:,:)*CSVAR
                  ENDDO 
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(VEC)
        DEALLOCATE(MAT)
      ENDDO
!
!     ==========================================================================
!     == PRINT DATA FILE FOR XMGRACE                                          ==
!     ==========================================================================
!!$      PRINT*,'STARTING TO WRITE IN DMFT_GNI',MU
!!$      OPEN(UNIT=NFIL,FILE=TRIM(FILE))
!!$      DO NU=1,NOMEGA
!!$        WRITE(NFIL,*)OMEGA(NU),(REAL(ATOMSET(IAT)%GNI(I,I,1,NU)) &
!!$                      ,AIMAG(ATOMSET(IAT)%GNI(I,I,1,NU)),I=1,NLOC)
!!$      ENDDO
!!$      CLOSE(NFIL)
!!$STOP 'FORCED AFTER DMFT_GNI'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_DOSHIST() 
!     **************************************************************************
!     ** DOS-HISTOGRAM                                                        **
!     ** DOS IS CALCULATED DIRECTLY FROM THE ENERGY VALUES OBTAINED FROM      **
!     ** DYNOCC SHIFTED RELATIVE TO THE FERMI LEVEL.                          **
!     ** -) BOTH SPINS ARE ADDED                                              **
!     ** -) ENERGY AXIS IN EV AND RELATIVE TO THE FERMI LEVEL                 **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NB,NKPTL,NSPIN,NDIM,NDIMD,NAT,NOMEGA,OMEGA &
     &                      ,KSET,ATOMSET,MU
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TON=.FALSE.
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      INTEGER(4),PARAMETER   :: NE=500
      REAL(8)   ,PARAMETER   :: EMIN=-2.D0/27.211D0
      REAL(8)   ,PARAMETER   :: EMAX=+2.D0/27.211D0
      COMPLEX(8)             :: CSVAR
      INTEGER(4)             :: ICHI1,ICHI2
      INTEGER(4)             :: NLOC
      INTEGER(4)             :: NU
      REAL(8)                :: WKPTL ! K-POINT WEIGHT
      REAL(8)                :: E     ! BAND ENERGY
      COMPLEX(8),ALLOCATABLE :: VEC(:,:)
      COMPLEX(8),ALLOCATABLE :: MAT(:,:)
      INTEGER(4)             :: IAT,IKPT,ISPIN,IB,I,IE
      INTEGER(4)             :: IDIM1,IDIM2,IDIMD
      INTEGER(4)             :: NFIL=10015
      CHARACTER(64)          :: FILE='FORDOS.DAT'
      REAL(8)   ,ALLOCATABLE :: DOS(:,:)
      REAL(8)                :: EV
      INTEGER(4)             :: NEG
      REAL(8)   ,ALLOCATABLE :: SMEAR(:)
      REAL(8)   ,ALLOCATABLE :: DOSCOPY(:,:)
      REAL(8)                :: X
!     **************************************************************************
      IF(.NOT.TON) RETURN
      CALL CONSTANTS('EV',EV)
      ALLOCATE(DOS(NE,NCHI))
      DOS=0.D0
      DO IAT=1,NAT
        ICHI1=ATOMSET(IAT)%ICHI1
        ICHI2=ATOMSET(IAT)%ICHI2
        NLOC=ATOMSET(IAT)%NLOC
PRINT*,'IAT,ICHI1,ICHI2,NLOC,NDIM,NDIMD,NSPIN ' &
&      ,IAT,ICHI1,ICHI2,NLOC,NDIM,NDIMD,NSPIN
IF(NLOC.LE.0) PRINT*,'SKIPPING ATOM=',IAT
IF(NLOC.LE.0) CYCLE
        ALLOCATE(VEC(NDIM,NLOC))
        ALLOCATE(MAT(NLOC,NLOC))
        DO IKPT=1,NKPTL
          WKPTL=KSET(IKPT)%WKPT
          DO ISPIN=1,NSPIN
            DO IB=1,NB
              E=KSET(IKPT)%E(IB,ISPIN)
              IE=1+INT((E-EMIN)/(EMAX-EMIN)*REAL(NE-1,KIND=8))
              IF(IE.LT.1.OR.IE.GT.NE) CYCLE
              VEC(:,:)=KSET(IKPT)%PIPSI(:,ICHI1:ICHI2,IB,ISPIN)
              DO IDIM1=1,NDIM
                DO IDIM2=1,NDIM
                  IDIMD=IDIM1+NDIM*(IDIM2-1)+ISPIN-1
                  DO I=1,NLOC
                    MAT(I,:)=VEC(IDIM1,I)*CONJG(VEC(IDIM2,:))
                  ENDDO
                  IF(KSET(IKPT)%TADDMINUSK) THEN
                    MAT=MAT+CONJG(TRANSPOSE(MAT(:,:)))
                  END IF
                  MAT=MAT*WKPTL
                  DO I=1,NLOC
                    DOS(IE,ICHI1-1+I)=DOS(IE,ICHI1-1+I)+REAL(MAT(I,I),KIND=8)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(VEC)
        DEALLOCATE(MAT)
      ENDDO
!
!     ==========================================================================
!     == SMEAR OUT                                                            ==
!     ==========================================================================
      NEG=50
      ALLOCATE(SMEAR(-NEG:NEG))
      DO IE=-NEG,NEG
        X=REAL(IE,KIND=8)/REAL(NEG,KIND=8)  ! [-1,1]
        X=X*4.D0
        SMEAR(IE)=EXP(-X**2)
      ENDDO
      SMEAR=SMEAR/SUM(SMEAR)
      ALLOCATE(DOSCOPY(1-NEG:NE+NEG,NCHI))
      DOSCOPY=0.D0
      DO IE=-NEG,NEG
        DOSCOPY(1+IE:NE+IE,:)=DOSCOPY(1+IE:NE+IE,:)+SMEAR(IE)*DOS(:,:)
      ENDDO
!
!     ==========================================================================
!     == PRINT DATA FILE FOR XMGRACE                                          ==
!     ==========================================================================
      PRINT*,'STARTING TO WRITE IN DMFT_DOSHIST'
      OPEN(UNIT=NFIL,FILE=TRIM(FILE))
      IAT=2
      NLOC=ATOMSET(IAT)%NLOC
      DO IE=-NEG,NE+NEG
        E=EMIN+(EMAX-EMIN)/REAL(NE-1,KIND=8)*REAL(IE-1,KIND=8)
        WRITE(NFIL,*)E/EV,DOSCOPY(IE,:)
      ENDDO
      CLOSE(NFIL)
STOP 'FORCED STOP AFTER DMFT_DOSHIST'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_MERMIN(NBANDS,NKPT,NSPIN &
     &                 ,TOTCHA,ISPINDEG,TEMP,WKPT,EIG,F,CHMPOT,EMERMN)
!      SUBROUTINE DYNOCC_MERMIN(NBANDS,NKPT,NSPIN &
!     &                 ,TOTCHA,ISPINDEG,TEMP,WKPT,EIG,F,CHMPOT,EMERMN)
!     ******************************************************************
!     **  CALCULATES THE OCCUPATIONS OF THE ELECTRONIC LEVELS         **
!     **  ACCORDING TO THE FERMI DISTRIBUTION;                        **
!     **  AND CALCULATES THE ENERGY -T*S RELATED TO THE ENTROPY OF    **
!     **  THE ELECTRONS.                                              **
!     **                                                              **
!     **  THE ROUTINE IS A COPY OF DYNOCC_MERMIN                      **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      LOGICAL(4) ,PARAMETER  :: TPR=.FALSE.
      INTEGER(4) ,PARAMETER  :: ITERX=1000   ! MAX #(ITERATIONS)
      REAL(8)    ,PARAMETER  :: TOL=1.D-10    ! TOLERANCE IN #(ELECTRONS)
      INTEGER(4) ,INTENT(IN) :: NBANDS       ! #(BANDS)
      INTEGER(4) ,INTENT(IN) :: NKPT         ! #(K-POINTS)    
      INTEGER(4) ,INTENT(IN) :: NSPIN        ! #(SPINS)    
      REAL(8)    ,INTENT(IN) :: TOTCHA       ! TOTAL CHARGE
      INTEGER(4) ,INTENT(IN) :: ISPINDEG     ! SPIN DEGENERACY (1 OR 2)
      REAL(8)    ,INTENT(IN) :: TEMP         ! K_B*TEMPERATURE IN HARTREE
      REAL(8)    ,INTENT(IN) :: WKPT(NKPT)         ! K-POINT WEIGHTS
      REAL(8)    ,INTENT(IN) :: EIG(NBANDS,NKPT,NSPIN) ! EIGENVALUES
      REAL(8)    ,INTENT(OUT):: F(NBANDS,NKPT,NSPIN)   ! OCCUPATIONS
      REAL(8)    ,INTENT(OUT):: CHMPOT               ! CHEMICAL POTENTIAL
      REAL(8)    ,INTENT(OUT):: EMERMN               ! HEAT OF THE ELECTRONS
      INTEGER(4)             :: ISTART
      INTEGER(4)             :: IB,I,ISPIN,IKPT
      REAL(8)                :: X0,DX,Y0,XM,YM
      REAL(8)                :: SVAR
      REAL(8)                :: SUMV
      INTEGER(4)             :: ITER       ! ITERATION COUNT
      REAL(8)                :: DQ         ! DEVIATION IN TOTAL CHARGE
      REAL(8)                :: EV         ! ELECTRON VOLT
      REAL(8)                :: DE
      REAL(8)                :: F1
      INTEGER(4)             :: IBI
      REAL(8)                :: FMAX         ! #(ELECTRONS PER STATE) 1 OR 2
!     ******************************************************************
                           CALL TRACE$PUSH('DYNOCC_MERMIN')
!
      IF(ISPINDEG.NE.1.AND.ISPINDEG.NE.2.D0) THEN
        CALL ERROR$MSG('SPIN-DEGENERACY CAN BE EITHER ONE OR TWO')
        CALL ERROR$STOP('DYNOCC_MERMIN')
      END IF
      FMAX=REAL(ISPINDEG,KIND=8)
!
!     ==================================================================
!     ==  ESTIMATE CHEMICAL POTENTIAL BY AVERAGING THE ONE-PARTICLE   ==
!     ==  ENERGIES OF THE HIGHEST OCCUPIED BAND                       ==
!     ==================================================================
      IB=INT(TOTCHA/FMAX)  
      IF(IB.GE.NBANDS*NSPIN) THEN
        CALL ERROR$MSG('TO FEW BANDS FOR THE NUMBER OF ELECTRONS')
        CALL ERROR$STOP('DYNOCC_MERMIN')
      END IF
      IB=MAX(IB,1)
      I=0
      CHMPOT=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          I=I+1
          CHMPOT=CHMPOT+EIG(IB,IKPT,ISPIN)
        ENDDO
      ENDDO
      CHMPOT=CHMPOT/DBLE(I) 
                           CALL TRACE$PASS('A')
!
!     ==================================================================
!     ==  FIND CHEMICAL POTENTIAL BY BISECTION                        ==
!     ==================================================================
      X0=CHMPOT
      DX=1.D-2
      ISTART=1
      CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
      CHMPOT=X0
      DO ITER=1,ITERX
        SUMV=0.D0
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB=1,NBANDS
              SVAR=(EIG(IB,IKPT,ISPIN)-CHMPOT)/TEMP
              IF(SVAR.GT.+50.D0)SVAR=+50.D0
              IF(SVAR.LT.-50.D0)SVAR=-50.D0
              F(IB,IKPT,ISPIN)=1.D0/(1.D0+EXP(SVAR))*FMAX
              SUMV=SUMV+F(IB,IKPT,ISPIN)*WKPT(IKPT)
            ENDDO
          ENDDO
        ENDDO
        DQ=SUMV-TOTCHA
        IF(ABS(DQ).LT.TOL) GOTO 110
        X0=CHMPOT
        Y0=DQ
        CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
        CHMPOT=X0
      ENDDO
      CALL ERROR$MSG('OCCUPATIONS NOT CONVERGED')
      CALL ERROR$MSG('PROBABLY THE TEMPERATURE IS ZERO')
      CALL ERROR$STOP('DYNOCC_MERMIN')
 110  CONTINUE
                           CALL TRACE$PASS('B')
!
!     ==================================================================
!     ==  CALCULATE HEAT OF THE ELECTRONS                             ==
!     ==================================================================
      EMERMN=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IB=1,NBANDS
            F1=F(IB,IKPT,ISPIN)/FMAX
            IF(F1.NE.0.D0.AND.1.D0-F1.NE.0.D0) THEN
              DE=+TEMP*(F1*LOG(F1)+(1.D0-F1)*LOG(1.D0-F1))
              EMERMN=EMERMN+DE*WKPT(IKPT)*FMAX
            END IF
          ENDDO
        ENDDO
      ENDDO
                           CALL TRACE$PASS('C')
!
!     ==================================================================
!     ==  PRINT FOR CHECK                                             ==
!     ==================================================================
      IF(TPR) THEN
        CALL CONSTANTS('EV',EV)
        WRITE(*,FMT='("#ELECTRONS( IN)=",F10.5' &
     &               //'," CHEMICAL POTENTIAL=",F10.3' &
     &               //'/"# ELECTRONS(OUT)=",F10.5)')TOTCHA,CHMPOT/EV,TOTCHA+DQ
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            WRITE(*,FMT='(5("(",F8.3,";",F5.2,")"))') &
     &         (EIG(IB,IKPT,ISPIN)/EV,F(IB,IKPT,ISPIN),IB=1,NBANDS)
          ENDDO
        ENDDO
      END IF
                         CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
MODULE DMFT_CHECKGRADIENT_MODULE
INTEGER(4) :: IKPT=1
INTEGER(4) :: ICHI1=1
INTEGER(4) :: ICHI2=1
INTEGER(4) :: IDIMD=1
REAL(8)    :: DEL=1.D-2
INTEGER(4) :: ISTEP=0
REAL(8)    :: ETOT0
REAL(8)    :: ETOT(4)
REAL(8)    :: DETOT
END MODULE DMFT_CHECKGRADIENT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_CHECKGRADIENT1()
      USE DMFT_MODULE, ONLY: KSET_TYPE,KSET,NKPTL,NB,NSPIN,NDIM,NCHI,NDIMD
      USE DMFT_CHECKGRADIENT_MODULE
      IMPLICIT NONE
      INTEGER(4) ::I
!     **************************************************************************
                                      CALL TRACE$PUSH('DMFT_CHECKGRADIENT1')
WRITE(*,FMT='(80("+"))')
WRITE(*,FMT='(80("+"),T10," CHECKGRADIENT 1 ")')
WRITE(*,FMT='(80("+"))')
      IF(ISTEP.EQ.0) THEN
        CONTINUE
      ELSE IF(ISTEP.EQ.1) THEN
        KSET(IKPT)%RHO(ICHI1,ICHI2,IDIMD)=KSET(IKPT)%RHO(ICHI1,ICHI2,IDIMD)-DEL
      ELSE IF(ISTEP.EQ.2) THEN
        KSET(IKPT)%RHO(ICHI1,ICHI2,IDIMD)=KSET(IKPT)%RHO(ICHI1,ICHI2,IDIMD)+DEL
      ELSE IF(ISTEP.EQ.3) THEN
        KSET(IKPT)%RHO(ICHI1,ICHI2,IDIMD)=KSET(IKPT)%RHO(ICHI1,ICHI2,IDIMD) &
     &                                -2.D0*DEL
      ELSE IF(ISTEP.EQ.4) THEN
        KSET(IKPT)%RHO(ICHI1,ICHI2,IDIMD)=KSET(IKPT)%RHO(ICHI1,ICHI2,IDIMD) &
     &                                +2.D0*DEL
      END IF
      DO I=1,NKPTL
        CALL SPINOR_PRINTMATRIX(6,'RHO IN CHECKGRADIENT 1 ' &
    &                      ,1,NCHI,NDIMD,NCHI,KSET(I)%RHO)
      ENDDO
PRINT*,'WKPT ',KSET(:)%WKPT
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_CHECKGRADIENT2(ETOTIN)
      USE DMFT_MODULE, ONLY: KSET_TYPE,KSET,NKPTL,NB,NSPIN,NDIM,NCHI,NDIMD
      USE DMFT_CHECKGRADIENT_MODULE
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: ETOTIN
      INTEGER(4) ::I
!     **************************************************************************
                                      CALL TRACE$PUSH('DMFT_CHECKGRADIENT1')
      DO I=1,NKPTL
        CALL SPINOR_PRINTMATRIX(6,'RHO IN CHECKGRADIENT 2 ' &
    &                      ,1,NCHI,NDIMD,NCHI,KSET(I)%RHO)
        CALL SPINOR_PRINTMATRIX(6,'GAMMA IN CHECKGRADIENT 2 ' &
    &                      ,1,NCHI,NDIMD,NCHI,KSET(I)%GAMMA)
      ENDDO
      IF(ISTEP.EQ.0) THEN
!       == TR[A*B] IN SPINOR REPRESENTATION IS 1/2 SUM_IDIMD TR[A(IDIM)*B(IDIM)]
!       == A IS GAMMA AND B IS DG ==============================================
        DETOT=0.5D0*REAL(KSET(IKPT)%GAMMA(ICHI1,ICHI2,IDIMD),KIND=8)
        ETOT0=ETOTIN
      ELSE IF(ISTEP.EQ.1) THEN
        ETOT(1)=ETOTIN
      ELSE IF(ISTEP.EQ.2) THEN
        ETOT(2)=ETOTIN
      ELSE IF(ISTEP.EQ.3) THEN
        ETOT(3)=ETOTIN
      ELSE IF(ISTEP.EQ.4) THEN
        ETOT(4)=ETOTIN
        PRINT*,'ANALYTIC ',DETOT &
   &    ,' NUMERIC ',(ETOT(2)-ETOT(1))/(2.D0*DEL)/KSET(IKPT)%WKPT &
   &    ,' NUMERIC ',(ETOT(4)-ETOT(3))/(4.D0*DEL)/KSET(IKPT)%WKPT
        PRINT*,'DISPLACEMENTS, ENERGIES ETOT0=',ETOT0
        PRINT*,-2.D0*DEL,(ETOT(3)-ETOT0),-2.D0*DEL*DETOT
        PRINT*,-1.D0*DEL,(ETOT(1)-ETOT0),-1.D0*DEL*DETOT
        PRINT*,+1.D0*DEL,(ETOT(2)-ETOT0),+1.D0*DEL*DETOT
        PRINT*,+2.D0*DEL,(ETOT(4)-ETOT0),+2.D0*DEL*DETOT

WRITE(*,FMT='(80("+"))')
WRITE(*,FMT='(80("+"),T10," CHECKGRADIENT 2 DONE ")')
WRITE(*,FMT='(80("+"))')
        STOP 'IN CHECKGRADIENT2'
      END IF
      ISTEP=ISTEP+1
WRITE(*,FMT='(80("+"))')
WRITE(*,FMT='(80("+"),T10," CHECKGRADIENT 2 DONE ")')
WRITE(*,FMT='(80("+"))')
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_STOREKSET()
      USE DMFT_MODULE, ONLY: KSET_TYPE,KSET,NKPTL,NB,NSPIN,NDIM,NCHI,NDIMD &
     &                      ,NAT,ATOMSET
      IMPLICIT NONE
      TYPE(KSET_TYPE),SAVE,ALLOCATABLE :: COPY(:) 
      LOGICAL(4)     ,SAVE :: TINI=.FALSE.
      INTEGER(4)           :: IKPT,IAT
!     **************************************************************************
                                      CALL TRACE$PUSH('DMFT_STOREKSET')
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        ALLOCATE(COPY(NKPTL))
        DO IKPT=1,NKPTL
          ALLOCATE(COPY(IKPT)%F(NB,NSPIN))
          ALLOCATE(COPY(IKPT)%PIPSI(NDIM,NCHI,NB,NSPIN))
          ALLOCATE(COPY(IKPT)%RHO(NCHI,NCHI,NDIMD))
          ALLOCATE(COPY(IKPT)%GAMMA(NCHI,NCHI,NDIMD))
          ALLOCATE(COPY(IKPT)%HRHO(NCHI,NCHI,NDIMD))
          ALLOCATE(COPY(IKPT)%SINV(NCHI,NCHI,NDIMD))
          ALLOCATE(COPY(IKPT)%SMAT(NCHI,NCHI,NDIMD))
          COPY(IKPT)%WKPT=KSET(IKPT)%WKPT
          COPY(IKPT)%TADDMINUSK=KSET(IKPT)%TADDMINUSK
          COPY(IKPT)%F=KSET(IKPT)%F
          COPY(IKPT)%PIPSI=KSET(IKPT)%PIPSI
          COPY(IKPT)%RHO=KSET(IKPT)%RHO
          COPY(IKPT)%GAMMA=KSET(IKPT)%GAMMA
          COPY(IKPT)%HRHO=KSET(IKPT)%HRHO
          COPY(IKPT)%SINV=KSET(IKPT)%SINV
          COPY(IKPT)%SMAT=KSET(IKPT)%SMAT
        ENDDO
      ELSE
        DO IKPT=1,NKPTL
          KSET(IKPT)%WKPT=COPY(IKPT)%WKPT
          KSET(IKPT)%TADDMINUSK=COPY(IKPT)%TADDMINUSK
          KSET(IKPT)%F=COPY(IKPT)%F
          KSET(IKPT)%PIPSI=COPY(IKPT)%PIPSI
          KSET(IKPT)%RHO=COPY(IKPT)%RHO
          KSET(IKPT)%GAMMA=COPY(IKPT)%GAMMA
          KSET(IKPT)%HRHO=COPY(IKPT)%HRHO
          KSET(IKPT)%SINV=COPY(IKPT)%SINV
          KSET(IKPT)%SMAT=COPY(IKPT)%SMAT
        ENDDO        
      END IF
PRINT*,'CAUTION!!!!: DMFT_STOREKSET SETS GAMMA, HRHO ETC. TO ZERO'
DO IKPT=1,NKPTL
   KSET(IKPT)%GAMMA=(0.D0,0.D0)
   KSET(IKPT)%HRHO=(0.D0,0.D0)
   KSET(IKPT)%SINV=(0.D0,0.D0)
   KSET(IKPT)%SMAT=(0.D0,0.D0)
ENDDO
!!$DO IAT=1,NAT
!!$  ATOMSET(IAT)%SLOC=(0.D0,0.D0)
!!$  ATOMSET(IAT)%SLOCLAUR=(0.D0,0.D0)
!!$ENDDO

                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_RHOOFK()
!     **************************************************************************
!     **  CALCULATES THE K-DEPENDENT DENSITY MATRIX AND THE OVERLAP MATRIX    **
!     **  FROM TB-PROJECTIONS OF THE NATURAL ORBITALS AND THEIR OCCUPATIONS   **
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
!!$IF(IKPT.EQ.1) CALL SPINOR_PRINTMATRIX(6,'RHO(UPDOWN)  IN RHOOFK',1,NCHI,NDIMD,NCHI &
!!$  &    ,KSET(IKPT)%RHO)
!!$CALL SPINOR_PRINTMATRIX(6,'SINV(UPDOWN) IN RHOOFK',1,NCHI,NDIMD,NCHI &
!!$   &    ,KSET(IKPT)%SINV)
!
!       ========================================================================
!       == CONVERT FROM UP-DOWN REPRESENTATION INTO (TXYZ)                    ==
!       ========================================================================
        CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%RHO)
        CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%SINV)
!!$IF(IKPT.EQ.1) THEN
!!$PRINT*,'NSPIN ',NSPIN,' NDIM ',NDIM,' NDIMD ',NDIMD
!!$CALL SPINOR_PRINTMATRIX(6,'RHO(T,S)  IN RHOOFK',1,NCHI,NDIMD,NCHI &
!!$   &    ,KSET(IKPT)%RHO)
!!$CALL SPINOR_PRINTMATRIX(6,'SINV(T,S) IN RHOOFK',1,NCHI,NDIMD,NCHI &
!!$   &    ,KSET(IKPT)%SINV)
!!$END IF
!
!       ========================================================================
!       == OVERLAP MATRIX                                                     ==
!       ========================================================================
        CALL SPINOR$INVERT(NDIMD,NCHI,KSET(IKPT)%SINV,KSET(IKPT)%SMAT)
!!$IF(IKPT.EQ.1) THEN
!!$CALL SPINOR_PRINTMATRIX(6,'SMAT(T,S) IN RHOOFK',1,NCHI,NDIMD,NCHI &
!!$   &    ,KSET(IKPT)%SMAT)
!!$CALL SPINOR_PRINTMATRIX(6,'MULTIPLICATIONTEST',1,NCHI,NDIMD,NCHI &
!!$   &    ,MATMUL(KSET(IKPT)%SINV(:,:,1),KSET(IKPT)%SMAT(:,:,1)))
!!$END IF
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
      LOGICAL(4),PARAMETER   :: TPRINT=.TRUE.
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
     &                            +KSET(IKPT)%WKPT*KSET(IKPT)%RHO(I1:I2,I1:I2,:)
          IPRO=IPRO+NLOC
        ENDDO ! END OF LOOP OVER ATOMS
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
!     ==  PRINT
!     ==========================================================================
      IF(TPRINT) THEN
        PRINT*,'DENSITY MATRIX REPORT FROM DMFT_RHOLOCAL'
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
      SUBROUTINE DMFT_OVERLAPOFK()
!     **************************************************************************
!     **  CALCULATES THE K-DEPENDENT OVERLAP MATRIX                           **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON &
     &                      ,NCHI &
     &                      ,NKPTL &
     &                      ,NSPIN &
     &                      ,NDIM &
     &                      ,NDIMD &
     &                      ,NAT &
     &                      ,KSET &
     &                      ,ATOMSET
      USE LMTO_MODULE, ONLY : OVERLAP
      USE MPE_MODULE
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER   :: TTEST=.FALSE.
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      REAL(8)                :: PI
      INTEGER(4)             :: NND
      REAL(8)                :: SVAR
      COMPLEX(8)             :: EIKR
      INTEGER(4)             :: IAT1,IAT2,IT(3)
      INTEGER(4)             :: I1,I2,J1,J2
      INTEGER(4)             :: IKPT,NN,IDIMD,I
INTEGER(4) :: NFIL=1500,NN2,LM1
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                      CALL TRACE$PUSH('DMFT_OVERLAPOFK')
      PI=4.D0*ATAN(1.D0)
      IF(NDIM.EQ.1) THEN
        NDIMD=NSPIN
      ELSE IF(NDIM.EQ.2) THEN
        NDIMD=4
      END IF
!
!     ==========================================================================
!     == CALCULATES OVERLAP MATRIX ON THE NEIGHBORLIST. 
!     == CHECKS IF OFFSITEX IS TRUE. STOPS, IF NOT
!     ==========================================================================
!THE ROUTINE LMTO_OVERLAPEVAL DOES NOT WORK BECAUSE THE ARRAY OFFSITEX 
!IS NOT SET UP DURING THE INITIALIZATION, WHEN GOING THROUGH THE DMFT OPTION.
CALL ERROR$MSG('FORCED STOP AFTER CALLING LMTO_OVERLAPEVAL IN THE DMFT OBJECT')
CALL ERROR$STOP('DMFT_OVERLAPOFK')
      CALL LMTO_OVERLAPEVAL()
CALL DMFT_REPORTOVERLAP(6)
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
!
!     ==========================================================================
!     ==  ADD UP DENSITY MATRIX                                               ==
!     ==========================================================================
      NND=SIZE(OVERLAP)
      DO IKPT=1,NKPTL
!!$OPEN(NFIL,FILE='X1.DAT')
!!$REWIND(NFIL)
!!$DO IDIMD=1,NDIMD
!!$  WRITE(NFIL,FMT='(82("="),T10," OLD OVERLAP FOR IK=",I2,"  ")')IKPT
!!$  DO I=1,NCHI
!!$    WRITE(NFIL,FMT='(200F10.5)')KSET(IKPT)%SMAT(I,:,IDIMD)
!!$  ENDDO
!!$ENDDO
!!$CLOSE(NFIL)
        IF(.NOT.ASSOCIATED(KSET(IKPT)%SMAT)) &
     &                        ALLOCATE(KSET(IKPT)%SMAT(NCHI,NCHI,NDIMD))
        KSET(IKPT)%SMAT(:,:,:)=(0.D0,0.D0)
        DO NN=1,NND
          IAT1=OVERLAP(NN)%IAT1
          IAT2=OVERLAP(NN)%IAT2
          IT=OVERLAP(NN)%IT
          SVAR=2.D0*PI*SUM(XK(:,IKPT)*REAL(IT))
          EIKR=EXP(CI*SVAR)  !<P_{R+T}|PSI>=<P_R|PSI>*EIKR
          I1=ATOMSET(IAT1)%ICHI1
          I2=ATOMSET(IAT1)%ICHI2
          J1=ATOMSET(IAT2)%ICHI1
          J2=ATOMSET(IAT2)%ICHI2
          IF(I2-I1+1.NE.OVERLAP(NN)%N1) THEN
            CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
            CALL ERROR$STOP('DMFT_OVERLAPOFK')
          END IF
          IF(J2-J1+1.NE.OVERLAP(NN)%N2) THEN
            CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
            CALL ERROR$STOP('DMFT_OVERLAPOFK')
          END IF
          KSET(IKPT)%SMAT(I1:I2,J1:J2,1)=KSET(IKPT)%SMAT(I1:I2,J1:J2,1) &
      &                               +EIKR*OVERLAP(NN)%MAT
        ENDDO  ! END OF LOOP OVER NEIGHBORLIST
!
!       ==  MAP INTO SPINOR FORM ===============================================
!       == CONVERT OVERLAP INTO TXYZ REPRESENTATION
        CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%SMAT)
!
!       == INVERT OVERLAP MATRIX ===============================================
        IF(.NOT.ASSOCIATED(KSET(IKPT)%SINV)) &
     &                        ALLOCATE(KSET(IKPT)%SINV(NCHI,NCHI,NDIMD))
        CALL SPINOR$INVERT(NDIMD,NCHI,KSET(IKPT)%SMAT,KSET(IKPT)%SINV)
!!$OPEN(NFIL,FILE='X2.DAT')
!!$REWIND(NFIL)
!!$DO IDIMD=1,NDIMD
!!$  WRITE(NFIL,FMT='(82("="),T10," NEW OVERLAP FOR IK=",I2,"  ")')IKPT
!!$  DO I=1,NCHI
!!$    WRITE(NFIL,FMT='(200F10.5)')KSET(IKPT)%SMAT(I,:,IDIMD)
!!$  ENDDO
!!$ENDDO
!!$CLOSE(NFIL)

      ENDDO ! END OF LOOP OVER KPOINTS
!STOP 'FORCED IN DMFT_OVERLAPOFK'
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_REPORTOVERLAP(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : OVERLAP
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NNS
      INTEGER(4)            :: NN,LM1
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTSBAR')
      NNS=SIZE(OVERLAP)
      WRITE(NFIL,FMT='(82("="))')
      WRITE(NFIL,FMT='(82("="),T10,"   SCREENED STRUCTURE CONSTANTS   ")')
      WRITE(NFIL,FMT='(82("="))')
      DO NN=1,NNS
        WRITE(NFIL,FMT='(82("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 OVERLAP(NN)%IAT1,OVERLAP(NN)%IAT2,OVERLAP(NN)%IT
        DO LM1=1,OVERLAP(NN)%N1
          WRITE(NFIL,FMT='(20F10.5)')OVERLAP(NN)%MAT(LM1,:)
        ENDDO
        WRITE(NFIL,FMT='(82("="))')
      ENDDO
                             CALL TRACE$POP()
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
print*,'---alpha--- ',atomset(iat)%lhfweight
print*,'---u------- ',atomset(iat)%u
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
      USE LMTO_MODULE, ONLY : ISPECIES &
     &                       ,POTPAR1 &
     &                       ,TCTE &
     &                       ,CTE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      INTEGER(4),INTENT(IN)  :: LMNX          !  #(LOCAL ORBITALS ON THIS SITE)
      REAL(8)   ,INTENT(OUT) :: U(LMNX,LMNX,LMNX,LMNX)
      REAL(8)   ,ALLOCATABLE :: U1(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: U2(:,:,:,:)
      INTEGER(4)             :: ISP     ! ATOM TYPE
      INTEGER(4)             :: LMNXT   ! #(VALENCE+SCATTERING WAVES)
      INTEGER(4)             :: IH,IT
!     **************************************************************************
      IF(.NOT.TCTE) THEN
        CALL DMFT_ULOCAL_NONCTE(IAT,LMNX,U)
        RETURN
      END IF
      ISP=ISPECIES(IAT)
      LMNXT=CTE(IAT)%LMNXT ! SIZE OF U-TENSOR ON POTPAR
      IF(LMNX.NE.CTE(IAT)%LMNXH) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
        CALL ERROR$MSG('LMNX MUST BE EQUAL TO THE NUMBER OF HEAD FUNCTIONS')
        CALL ERROR$I4VAL('LMNX',LMNX)
        CALL ERROR$I4VAL('LMNXH',CTE(IAT)%LMNXH)
        CALL ERROR$I4VAL('ISP',ISP)
        CALL ERROR$STOP('DMFT_ULOCAL')
      END IF
!
!     ==========================================================================
!     == COLLECT U-TENSOR IN EXPANDED BASIS                                   ==
!     ==========================================================================
      ALLOCATE(U1(LMNXT,LMNXT,LMNXT,LMNXT))
      U1=POTPAR1(ISP)%TAILED%U
!
!     ==========================================================================
!     == COLLECT U-TENSOR IN EXPANDED BASIS                                   ==
!     ==========================================================================
      ALLOCATE(U2(LMNXT,LMNXT,LMNXT,LMNX))
      U2(:,:,:,:)=0.D0
      DO IH=1,LMNX
        DO IT=1,LMNXT
          U2(:,:,:,IH)=U2(:,:,:,IH)+U1(:,:,:,IT)*CTE(IAT)%MAT(IT,IH)
        ENDDO
      ENDDO
      U1(:,:,:,:)=0.D0
      DO IH=1,LMNX
        DO IT=1,LMNXT
          U1(:,:,IH,:LMNX)=U1(:,:,IH,:LMNX) &
     &                    +U2(:,:,IT,:LMNX)*CTE(IAT)%MAT(IT,IH)
        ENDDO
      ENDDO
      U2(:,:,:,:)=0.D0
      DO IH=1,LMNX
        DO IT=1,LMNXT
          U2(:,IH,:LMNX,:LMNX)=U2(:,IH,:LMNX,:LMNX) &
     &                        +U1(:,IT,:LMNX,:LMNX)*CTE(IAT)%MAT(IT,IH)
        ENDDO
      ENDDO
      U(:,:,:,:)=0.D0
      DO IH=1,LMNX
        DO IT=1,LMNXT
          U(IH,:,:,:)=U(IH,:,:,:) &
     &               +U2(IT,:LMNX,:LMNX,:LMNX)*CTE(IAT)%MAT(IT,IH)
        ENDDO
      ENDDO
      DEALLOCATE(U1)
      DEALLOCATE(U2)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_ULOCAL_NONCTE(IAT,LMNX,U)
!     **************************************************************************
!     **  CALCULATES THE U-TENSOR OF THE NATURAL TIGHT-BINDING ORBITALS       **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES &
     &                       ,POTPAR1 &
     &                       ,SBAR_NEW
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
      SUBROUTINE DMFT_NATORB()
!     **************************************************************************
!     ** NATURAL ORBITALS IN A TWO COMPONENT PAULI-SPINOR REPRESENTATION      **
!     **  USES KSET%RHO AND KSET%SINV 
!     **************************************************************************
      USE DMFT_MODULE, ONLY : NAT,NKPTL,NDIMD,ATOMSET,KSET
      IMPLICIT NONE
!     __________________________________ORBTYPE CAN BE 'NATORB','QUAMBO','ORTHO'
      CHARACTER(16),PARAMETER:: ORBTYPE='QUAMBO' 
!      CHARACTER(16),PARAMETER:: ORBTYPE='NATORB'
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
!!$CALL SPINOR_PRINTMATRIX(6,'RHO(T,S)  IN NATORB',1,NLOC,NDIMD,NLOC &
!!$   &    ,RHO)
!!$CALL SPINOR_PRINTMATRIX(6,'SINV(T,S) IN NATORB',1,NLOC,NDIMD,NLOC &
!!$   &    ,SINV)
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
!       == THE CHOICE OF THE MATRIX RHO SHOULD BE ADAPTED TO THE NORM OF 
!       == THE ORBITALS
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
!
WRITE(*,FMT='(100("="),T30,"  ",A,I3,"  ")')' ATOMSET%NATORB%PIPHI ATOM=',IAT
DO I=1,2*NLOC
  WRITE(*,FMT='(100("(",2F10.5,")"))')ATOMSET(IAT)%NATORB%PIPHI(I,:)
ENDDO
WRITE(*,FMT='(100("="),T30,"  ",A,I3,"  ")')' ATOMSET%NATORB%CHIPHI ATM=',IAT
DO I=1,2*NLOC
  WRITE(*,FMT='(100("(",2F10.5,")"))')ATOMSET(IAT)%NATORB%CHIPHI(I,:)
ENDDO

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
      USE DMFT_MODULE ,ONLY: NKPTL,NCHI,NDIMD,NOMEGA,NLAU,OMEGA,KBT,MU,KSET &
     &                     ,TONATORB
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      REAL(8)   ,PARAMETER   :: TOL=1.D-5  ! TOLERANCE FOR THE DENSITY MATRIX
!      REAL(8)   ,PARAMETER   :: EMAXBYKT=4.D0
      REAL(8)   ,PARAMETER   :: EMAXBYKT=6.D0
      COMPLEX(8),ALLOCATABLE :: MAT(:,:,:)
      COMPLEX(8),ALLOCATABLE :: A(:,:)
      COMPLEX(8),ALLOCATABLE :: B(:,:)
      REAL(8)   ,ALLOCATABLE :: F(:)   !OCCUPATIONS
      REAL(8)   ,ALLOCATABLE :: E(:)   !ENERGIES
      COMPLEX(8),ALLOCATABLE :: U(:,:)
      INTEGER(4)             :: N
      REAL(8)                :: SVAR
      INTEGER(4)             :: IKPT,ISPIN,I
REAL(8)    :: ev
INTEGER(4) :: IDIMD,ICHI
!     **************************************************************************
      CALL CONSTANTS('EV',EV)
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
      ALLOCATE(E(N))
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
!
!           == KEEP TRANSFORM TO NATURAL ORBITALS ==============================
            IF(TONATORB) THEN
              IF(.NOT.ALLOCATED(KSET(IKPT)%CHIPHI)) &
      &                          ALLOCATE(KSET(IKPT)%CHIPHI(NCHI,NCHI,NDIMD))
              KSET(IKPT)%CHIPHI(:,:,ISPIN)=U
            END IF
!
!           == CALCULATE HRHO ==================================================
            B=TRANSPOSE(CONJG(U))
            DO I=1,NCHI
              CALL DMFT_EOFF(KBT,F(I),SVAR)
              SVAR=MAX(-EMAXBYKT*KBT,MIN(EMAXBYKT*KBT,SVAR))
              SVAR=MU+SVAR
              E(I)=SVAR
              B(I,:)=SVAR*B(I,:)
WRITE(*,FMT='("IN HRHO: K=",I5," S=",I1," E[EV]=",F15.5," OCC=",2F10.5)') &
&   IKPT,ISPIN,E(I)/EV,F(I),1.D0/(1.D0+EXP((E(I)-MU)/KBT))
            ENDDO ! ICHI
            KSET(IKPT)%HRHO(:,:,ISPIN)=MATMUL(U,B)
!
!           == CALCULATE N-REPRESENTABLE DENSITY MATRIX ========================
            U=MATMUL(KSET(IKPT)%SINV(:,:,ISPIN),U)
            B=TRANSPOSE(CONJG(U))
            DO I=1,NCHI
              SVAR=1.D0/(1.D0+EXP((E(I)-MU)/KBT))
              B(I,:)=SVAR*B(I,:)
            ENDDO ! ICHI
            KSET(IKPT)%RHOADAPTED(:,:,ISPIN)=MATMUL(U,B)
          ENDDO ! ISPIN
          CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%RHO)
          CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%RHOADAPTED)
          CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%SINV)
          CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%HRHO)
          IF(TONATORB) THEN
            CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%CHIPHI)
          END IF
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
!
!         == KEEP TRANSFORM TO NATURAL ORBITALS ================================
          IF(TONATORB) THEN
            IF(.NOT.ALLOCATED(KSET(IKPT)%CHIPHI)) &
      &                             ALLOCATE(KSET(IKPT)%CHIPHI(NCHI,NCHI,NDIMD))
            KSET(IKPT)%CHIPHI(:,:,ISPIN)=U
            CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%CHIPHI)
          END IF
!         === CONSTRUCT HRHO
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
!CALL DMFT_HRHOTEST2()
!
!     ==========================================================================
!     ==  TEST                                                                ==
!     ==========================================================================
!!$      IF(TTEST) THEN 
!!$        IF(.NOT.ALLOCATED(MAT))ALLOCATE(MAT(NCHI,NCHI,NDIMD))
!!$        DO IKPT=1,NKPTL
!!$!         == CALCULATE THE DENSITY MATRIX FROM THE GREENS FUNCTION WITH HRHO ===
!!$          CALL DMFT_HRHO_TEST(NCHI,NDIMD,NOMEGA,NLAU,OMEGA,KBT &
!!$     &                       ,KSET(IKPT)%SMAT,KSET(IKPT)%HRHO,MAT) !MAT=RHO
!!$!         == CALCULATE DEVIATION FROM TARGET DENSITY MATRIX ====================
!!$          IF(MAXVAL(ABS(MAT-KSET(IKPT)%RHO)).GT.TOL) THEN
!!$IF(IKPT.EQ.1) THEN
!!$  DO IDIMD=1,NDIMD
!!$    WRITE(*,FMT='(82("="),T10," RHOADAPTED FOR IDIMD=",I2,"  ")')IDIMD
!!$    DO ICHI=1,6
!!$      WRITE(*,FMT='(200F10.5)')KSET(IKPT)%RHOADAPTED(ICHI,:6,IDIMD)
!!$    ENDDO
!!$    WRITE(*,FMT='(82("="),T10," RHO1 FOR IDIMD=",I2,"  ")')IDIMD
!!$    DO ICHI=1,6
!!$      WRITE(*,FMT='(200F10.5)')KSET(IKPT)%RHO(ICHI,:6,IDIMD)
!!$    ENDDO
!!$    WRITE(*,FMT='(82("="),T10," RHO2 FOR IDIMD=",I2,"  ")')IDIMD
!!$    DO ICHI=1,6
!!$      WRITE(*,FMT='(200F10.5)')MAT(ICHI,:6,IDIMD)
!!$    ENDDO
!!$  ENDDO
!!$END IF
!!$            CALL SPINOR_PRINTMATRIX(6,'RHO[HRHO]-RHO',1,NCHI &
!!$     &                           ,NDIMD,NCHI,MAT-KSET(IKPT)%RHO)
!!$            CALL ERROR$MSG('TEST OF HRHO FAILED')
!!$            CALL ERROR$MSG('INCREASE THE NUMBER OF MATSUBARA FREQUENCIES')
!!$            CALL ERROR$R8VAL('MAX DEV OF RHO',MAXVAL(ABS(MAT-KSET(IKPT)%RHO)))
!!$            CALL ERROR$STOP('DMFT_HRHO')
!!$          END IF
!!$        ENDDO
!!$      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_TONATORB(ID)
!     **************************************************************************
!     **  CALCULATES A STATIC HAMILTONIAN HRHO SUCH THAT THE RESULTING        **
!     **  DENSITY MATRIX IS RHO=[1+EXP(BETA*HRHO)]^(-1)                       **
!     **  THIS VERSION WORKS WITH THE EIGENVALUES FROM THE DYNOCC OBJECT      **
!     **  AND IS USED FOR TESTING ONLY
!     **************************************************************************
      USE DMFT_MODULE ,ONLY: TONATORB,NKPTL,NCHI,NDIMD,KSET
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)              :: IKPT
      COMPLEX(8)              :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)              :: MAT2(NCHI,NCHI,NDIMD)
!     **************************************************************************
      IF(.NOT.TONATORB) RETURN
!
!     ==========================================================================
!     ==  COMPLETE TRANSFORMATION                                             ==
!     ==========================================================================
      DO IKPT=1,NKPTL
!       == PREPARE PIPHI=SINV*U ================================================
        IF(.NOT.ALLOCATED(KSET(IKPT)%CHIPHI)) THEN
          CALL ERROR$MSG('TRANSFORMATION NOT AVAILABLE')
          CALL ERROR$STOP('DMFT_TONATORB')
        END IF
        IF(.NOT.ALLOCATED(KSET(IKPT)%PIPHI)) &
     &                               ALLOCATE(KSET(IKPT)%PIPHI(NCHI,NCHI,NDIMD))
        CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,KSET(IKPT)%CHIPHI &
     &                                                 ,KSET(IKPT)%PIPHI)
      ENDDO
!      
!     ==========================================================================
!     ==  PERFORM TRANSFORMATION                                              ==
!     ==========================================================================
      IF(ID.EQ.'FWRD') THEN      
        DO IKPT=1,NKPTL
!         == TRANSFORM HRHO ====================================================
          CALL SPINOR$CONJUGATE(NDIMD,NCHI,KSET(IKPT)%PIPHI,MAT1)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,KSET(IKPT)%HRHO,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,KSET(IKPT)%PIPHI,KSET(IKPT)%HRHO)
!         == TRANSFORM GAMMA ===================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,KSET(IKPT)%GAMMA,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,KSET(IKPT)%PIPHI,KSET(IKPT)%GAMMA)
!         == TRANSFORM SMAT ====================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,KSET(IKPT)%SMAT,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,KSET(IKPT)%PIPHI,KSET(IKPT)%SMAT)
!         == TRANSFORM RHO =====================================================
          CALL SPINOR$CONJUGATE(NDIMD,NCHI,KSET(IKPT)%CHIPHI,MAT1)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,KSET(IKPT)%RHO,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,KSET(IKPT)%CHIPHI,KSET(IKPT)%RHO)
!         == TRANSFORM RHOADAPTED ==============================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,KSET(IKPT)%RHOADAPTED,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,KSET(IKPT)%CHIPHI &
     &                                                 ,KSET(IKPT)%RHOADAPTED)
!         == TRANSFORM SINV ====================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,KSET(IKPT)%SINV,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,KSET(IKPT)%CHIPHI,KSET(IKPT)%SINV)
!         == PRINT TRANSFORMED MATRICES
!!$          CALL SPINOR$PRINTMATRIX(6,'RHO','TXYZ' &
!!$     &                           ,1,NCHI,NDIMD,NCHI,KSET(IKPT)%RHO)
!!$          CALL SPINOR$PRINTMATRIX(6,'RHOADAPTED','TXYZ' &
!!$     &                           ,1,NCHI,NDIMD,NCHI,KSET(IKPT)%RHOADAPTED)
!!$          CALL SPINOR$PRINTMATRIX(6,'HRHO','TXYZ' &
!!$     &                           ,1,NCHI,NDIMD,NCHI,KSET(IKPT)%HRHO)
!!$          CALL SPINOR$PRINTMATRIX(6,'SMAT','TXYZ' &
!!$     &                           ,1,NCHI,NDIMD,NCHI,KSET(IKPT)%SMAT)
!!$          CALL SPINOR$PRINTMATRIX(6,'SINV','TXYZ' &
!!$     &                           ,1,NCHI,NDIMD,NCHI,KSET(IKPT)%SINV)
        ENDDO
      ELSE IF(ID.EQ.'BACK') THEN
        DO IKPT=1,NKPTL
!         == TRANSFORM HRHO ====================================================
          CALL SPINOR$CONJUGATE(NDIMD,NCHI,KSET(IKPT)%CHIPHI,MAT1)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%HRHO,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%CHIPHI,MAT2,KSET(IKPT)%HRHO)
!         == TRANSFORM GAMMA ===================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%GAMMA,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%CHIPHI,MAT2,KSET(IKPT)%GAMMA)
!         == TRANSFORM SMAT ====================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SMAT,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%CHIPHI,MAT2,KSET(IKPT)%SMAT)
!         == TRANSFORM RHO =====================================================
          CALL SPINOR$CONJUGATE(NDIMD,NCHI,KSET(IKPT)%PIPHI,MAT1)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%RHO,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%PIPHI,MAT2,KSET(IKPT)%RHO)
!         == TRANSFORM RHOADAPTED ==============================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%RHOADAPTED,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%PIPHI,MAT2 &
     &                                                   ,KSET(IKPT)%RHOADAPTED)
!         == TRANSFORM SINV ====================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%PIPHI,MAT2,KSET(IKPT)%SINV)
!         == PRINT TRANSFORMED MATRICES
!!$          CALL SPINOR$PRINTMATRIX(6,'RHO','TXYZ' &
!!$     &                           ,1,NCHI,NDIMD,NCHI,KSET(IKPT)%RHO)
!!$          CALL SPINOR$PRINTMATRIX(6,'RHOADAPTED','TXYZ' &
!!$     &                           ,1,NCHI,NDIMD,NCHI,KSET(IKPT)%RHOADAPTED)
!!$          CALL SPINOR$PRINTMATRIX(6,'HRHO','TXYZ' &
!!$     &                           ,1,NCHI,NDIMD,NCHI,KSET(IKPT)%HRHO)
!!$          CALL SPINOR$PRINTMATRIX(6,'SMAT','TXYZ' &
!!$     &                           ,1,NCHI,NDIMD,NCHI,KSET(IKPT)%SMAT)
!!$          CALL SPINOR$PRINTMATRIX(6,'SINV','TXYZ' &
!!$     &                           ,1,NCHI,NDIMD,NCHI,KSET(IKPT)%SINV)
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('DMFT_TONATORB')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_HRHOTEST2()
!     **************************************************************************
!     **  CALCULATES A STATIC HAMILTONIAN HRHO SUCH THAT THE RESULTING        **
!     **  DENSITY MATRIX IS RHO=[1+EXP(BETA*HRHO)]^(-1)                       **
!     **  THIS VERSION WORKS WITH THE EIGENVALUES FROM THE DYNOCC OBJECT      **
!     **  AND IS USED FOR TESTING ONLY
!     **************************************************************************
      USE DMFT_MODULE ,ONLY: NKPTL,NSPIN,NB,NCHI,NDIM,NDIMD,KSET
      IMPLICIT NONE
      LOGICAL(4),PARAMETER :: TON=.TRUE. 
      INTEGER(4)           :: IKPT,IB,ICHI1,ICHI2,IDIM1,IDIM2,IDIMD,ISPIN
!     **************************************************************************
      IF(.NOT.TON) RETURN
!
!     ==========================================================================
!     ==  PRINT RESULT                                                        ==
!     ==========================================================================
      DO IKPT=1,1
        DO IDIMD=1,NDIMD
          WRITE(*,FMT='(82("="),T10,"  ORIGINAL HRHO IDIMD=",I2," K=",I3,"  ")')
          DO ICHI1=1,NCHI
            WRITE(*,FMT='(200F10.5)')KSET(IKPT)%HRHO(ICHI1,:,IDIMD)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  COMPOSE NEW HRHO
!     ==========================================================================
      DO IKPT=1,1
        KSET(IKPT)%HRHO=(0.D0,0.D0)
        DO ISPIN=1,NSPIN
          DO IB=1,NB
            DO ICHI1=1,NCHI
              DO ICHI2=1,NCHI
                DO IDIM1=1,NDIM
                  DO IDIM2=1,NDIM
                    IDIMD=IDIM1+NDIM*(IDIM2-1)+ISPIN-1
                    KSET(IKPT)%HRHO(ICHI1,ICHI2,IDIMD) &
     &                      =KSET(IKPT)%HRHO(ICHI1,ICHI2,IDIMD) &
     &                          +CONJG(KSET(IKPT)%PIPSI(IDIM1,ICHI1,IB,ISPIN)) &
     &                               * KSET(IKPT)%E(IB,ISPIN) &
     &                               * KSET(IKPT)%PIPSI(IDIM2,ICHI2,IB,ISPIN)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  PRINT RESULT                                                        ==
!     ==========================================================================
      DO IKPT=1,1
        DO IDIMD=1,NDIMD
          WRITE(*,FMT='(82("="),T10,"  TEST HRHO IDIMD=",I2," K=",I3,"  ")')
          DO ICHI1=1,NCHI
            WRITE(*,FMT='(200F10.5)')KSET(IKPT)%HRHO(ICHI1,:,IDIMD)
          ENDDO
        ENDDO
      ENDDO
      CALL ERROR$MSG('FORCED STOP')
      CALL ERROR$STOP('DMFT_HRHOTEST2')
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
        IF(F.LT.0.5D0) THEN
          F0=FMIN
        ELSE
          F0=1.D0-FMIN
        END IF
        VAL=LOG( (1.D0-F0)/F0 )
        SLOPE=-1.D0/( F0*(1.D0-F0) )
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
      REAL(8)                  :: FN(NLAU+1)
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
      COMPLEX(8)               :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GLAUR(NCHI,NCHI,NDIMD,NLAU+1)
      COMPLEX(8)               :: SLAUR(NCHI,NCHI,NDIMD,NLAU)
      COMPLEX(8)               :: G(NCHI,NCHI,NDIMD)
      REAL(8)                  :: WKPTL
      INTEGER(4)               :: IKPT,NU,IAT,I1,I2,IDIMD,ILAU
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_GLOC')
      IF(NLAU.GT.2) THEN
        CALL ERROR$MSG('NLAU EXCEEDS MAXIMUM VALUE OF 2')
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
!
!       ========================================================================
!       == LAURENT EXPANSION FOR THE GREENS FUNCTION ===========================
!       ========================================================================
        GLAUR(:,:,:,1)=KSET(IKPT)%SINV
        IF(NLAU+1.GE.2) THEN
          MAT2=-MU*KSET(IKPT)%SMAT+KSET(IKPT)%HRHO &
     &                            -KSET(IKPT)%GAMMA+SLAUR(:,:,:,1)
!         == GLAUR(2)=SINV*MAT2*SINV ===========================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,KSET(IKPT)%SINV,MAT)
          CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,MAT,GLAUR(:,:,:,2))
          IF(NLAU+1.GE.3) THEN
!           == GLAUR(3)=SINV*MAT2*SINV*MAT2*SINV+SINV*SLAUR(2)*SINV
            CALL SPINOR$MATMUL(NDIMD,NCHI,GLAUR(:,:,:,2),MAT,GLAUR(:,:,:,3))
            CALL SPINOR$MATMUL(NDIMD,NCHI,KSET(IKPT)%SINV,SLAUR(:,:,:,2),MAT2)
            CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,KSET(IKPT)%SINV,MAT)
            GLAUR(:,:,:,3)=GLAUR(:,:,:,3)+MAT
          END IF
        END IF

        IF(KSET(IKPT)%TADDMINUSK) THEN
          DO ILAU=1,NLAU+1
            CALL SPINOR$CONJUGATE(NDIMD,NCHI,GLAUR(:,:,:,ILAU),MAT)
            GLAUR(:,:,:,ILAU)=0.5D0*(GLAUR(:,:,:,ILAU)+MAT)
          ENDDO
        END IF
!
!       == PROJECT LAURENT TERMS ONTO LOCAL HILBERT SPACE ======================
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          DO ILAU=1,NLAU+1
            ATOMSET(IAT)%GLOCLAUR(:,:,:,ILAU)=ATOMSET(IAT)%GLOCLAUR(:,:,:,ILAU)&
     &                                       +WKPTL*GLAUR(I1:I2,I1:I2,:,ILAU)
          ENDDO
        ENDDO
!
!       ========================================================================
!       == PERFORM MATSUBARA SUM                                              ==
!       ========================================================================
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
            CALL SPINOR$CONJUGATE(NDIMD,NCHI,MAT,MAT2)
            MAT=MAT2
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
          CALL SPINOR$CONJUGATE(NDIMD,I2-I1+1 &
     &                    ,ATOMSET(IAT)%GLOCLAUR(:,:,:,ILAU),MAT(I1:I2,I1:I2,:))
          ATOMSET(IAT)%GLOCLAUR(:,:,:,ILAU)=0.5D0 &
     &                   *(ATOMSET(IAT)%GLOCLAUR(:,:,:,ILAU)+MAT(I1:I2,I1:I2,:))
        ENDDO
      ENDDO
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_fullLOCALSPECTRALFUNCTION()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NDIMD,NAT,MU,KSET,ATOMSET,nomega,omega
      use strings_module
      IMPLICIT NONE
      TYPE SPECTR_TYPE
        complex(8),ALLOCATABLE :: MAT(:,:,:,:) ! NORB,NORB,NE 
      END TYPE SPECTR_TYPE
      LOGICAL(4),PARAMETER :: TON=.true.
      COMPLEX(8),PARAMETER :: CI=(0.D0,1.D0)  ! SQRT(-1)
      REAL(8)   ,PARAMETER :: EMIN=-15.D0/27.211D0
      REAL(8)   ,PARAMETER :: EMAX=+15.D0/27.211D0
      INTEGER(4),PARAMETER :: NE=1000
      INTEGER(4)           :: NORB
      REAL(8)   ,allocatable  :: E(:) !(NE)
      TYPE(SPECTR_TYPE)     :: SPECTR(NAT)
      complex(8),ALLOCATABLE:: spectral(:,:,:) !(NORB,NORB,NOMEGA)
      complex(8),ALLOCATABLE:: GREEN(:,:,:) !(NORB,NORB,NOMEGA)
      integer(4),save       :: icallcount
      integer(4)            :: nfil
      character(128)        :: file
      character(128)        :: string
      INTEGER(4)            :: IE,IAT,IDIMD,I,nu
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_FULLLOCALSPECTRALFUNCTION')
      IF(.NOT.TON) RETURN
      DO IAT=1,NAT
        NORB=ATOMSET(IAT)%ICHI2-ATOMSET(IAT)%ICHI1+1
        ALLOCATE(SPECTR(IAT)%MAT(NORB,NORB,NDIMD,NE))
        SPECTR(IAT)%MAT=(0.D0,0.d0)
      ENDDO
!     == ENERGY GRID ===========================================================
      allocate(e(ne))
      DO IE=1,NE
        E(IE)=EMIN+(EMAX-EMIN)/REAL(NE-1,KIND=8)*REAL(IE-1,KIND=8)
      ENDDO
!
!     ==========================================================================
!     == PERFORM BRILLOUIN-ZONE SAMPLING                                      ==
!     ==========================================================================
      DO IAT=1,NAT
        NORB=ATOMSET(IAT)%ICHI2-ATOMSET(IAT)%ICHI1+1
        ALLOCATE(GREEN(NORB,NORB,NOMEGA))
        ALLOCATE(SPECTRAL(NORB,NORB,NE))
        DO IDIMD=1,NDIMD
          DO NU=1,NOMEGA
            GREEN(:,:,NU)=ATOMSET(IAT)%GLOC(:,:,IDIMD,NU)
          ENDDO
          CALL DMFT_ANALYTICCONTINUATION(NORB,NOMEGA,OMEGA,GREEN &
     &                                  ,NE,EMIN,EMAX,SPECTRAL)
          DO IE=1,NE
            SPECTR(IAT)%MAT(:,:,IDIMD,IE)=SPECTRAL(:,:,IE)
          ENDDO
        ENDDO
        DEALLOCATE(GREEN)
        DEALLOCATE(SPECTRAL)
      ENDDO
!
!     ==========================================================================
!     ==  WRITE SPECTRAL FUNCTION                                             ==
!     ==========================================================================
      ICALLCOUNT=ICALLCOUNT+1
      DO IAT=1,NAT
        WRITE(STRING,*)IAT
        FILE='_fullSPECTRALFUNCTION_AT'//TRIM(ADJUSTL(STRING))
        WRITE(STRING,*)ICALLCOUNT
        FILE=TRIM(ADJUSTL(FILE))//'_I'//TRIM(ADJUSTL(STRING))
        FILE=TRIM(ADJUSTL(FILE))//'.DAT'
        CALL FILEHANDLER$SETFILE('SPECTRALFILE',.TRUE.,-FILE)
        IF(ICALLCOUNT.EQ.1) THEN
          CALL FILEHANDLER$SETSPECIFICATION('SPECTRALFILE','STATUS','UNKNOWN')
          CALL FILEHANDLER$SETSPECIFICATION('SPECTRALFILE','ACTION','WRITE')
          CALL FILEHANDLER$SETSPECIFICATION('SPECTRALFILE','FORM','FORMATTED')
        END IF
        CALL FILEHANDLER$UNIT('SPECTRALFILE',NFIL)

        DO IE=1,NE
          NORB=ATOMSET(IAT)%ICHI2-ATOMSET(IAT)%ICHI1+1
          WRITE(NFIL,FMT='(200F15.5)')E(IE)*27.211D0 &
!     &       ,(REAL(SPECTR(IAT)%MAT(I,I,1,IE)),I=1,NORB) &
     &       ,(AIMAG(SPECTR(IAT)%MAT(I,I,1,IE)),I=1,NORB)
        ENDDO
        CALL FILEHANDLER$CLOSE('SPECTRALFILE')
      ENDDO
STOP 'FORCED AFTER WRITING SPECTRAL FUNCTION'
!
!     ==========================================================================
!     ==  CLEAN                                                               ==
!     ==========================================================================
      DO IAT=1,NAT
        DEALLOCATE(SPECTR(IAT)%MAT)
      ENDDO
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_LOCALSPECTRALFUNCTION()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NDIMD,NAT,MU,KSET,ATOMSET
      use strings_module
      IMPLICIT NONE
      LOGICAL(4),PARAMETER :: TON=.true.
      COMPLEX(8),PARAMETER :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER :: TPRINT=.FALSE.
      LOGICAL(4),PARAMETER :: TTEST=.FALSE.
      COMPLEX(8)           :: MATX(NCHI,NCHI,NDIMD)
      COMPLEX(8)           :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)           :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)           :: G(NCHI,NCHI,NDIMD)
      REAL(8)              :: WKPTL
      INTEGER(4)           :: IKPT,IE,IAT,I1,I2,IDIMD,I
      INTEGER(4)           :: NORB
      real(8)   ,parameter :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER :: EMIN=-15.D0/27.211D0
      REAL(8)   ,PARAMETER :: EMAX=+15.D0/27.211D0
      INTEGER(4),PARAMETER :: NE=1000
      REAL(8)              :: E(NE)
      REAL(8)   ,PARAMETER :: DELTA=0.1D0/27.211D0
      TYPE SPECTR_TYPE
        REAL(8),ALLOCATABLE :: MAT(:,:,:,:) ! NORB,NORB,NE 
      END TYPE SPECTR_TYPE
      TYPE(SPECTR_TYPE)     :: SPECTR(NAT)
      integer(4),save       :: icallcount
      integer(4)            :: nfil
      character(128)        :: file
      character(128)        :: string
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_LOCALSPECTRALFUNCTION')
      IF(.NOT.TON) RETURN
      DO IAT=1,NAT
        NORB=ATOMSET(IAT)%ICHI2-ATOMSET(IAT)%ICHI1+1
PRINT*,'NORB',NORB,IAT
        ALLOCATE(SPECTR(IAT)%MAT(NORB,NORB,NDIMD,NE))
        SPECTR(IAT)%MAT=0.D0
      ENDDO
!     == ENERGY GRID ===========================================================
      DO IE=1,NE
        E(IE)=EMIN+(EMAX-EMIN)/REAL(NE-1,KIND=8)*REAL(IE-1,KIND=8)
      ENDDO
PRINT*,'MU ',MU
!
!     ==========================================================================
!     == PERFORM BRILLOUIN-ZONE SAMPLING                                      ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        WKPTL=KSET(IKPT)%WKPT
        DO IE=1,NE
!         == CONSTRUCT LATTICE GREENS FUNCTION =================================
          MAT=(E(IE)+MU-CI*DELTA)*KSET(IKPT)%SMAT
          MAT=MAT-KSET(IKPT)%HRHO+KSET(IKPT)%GAMMA
          CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
!
!         == PROJECT GREENS FUNCTION ONTO INFIVIDUAL SITES =====================
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            SPECTR(IAT)%MAT(:,:,:,IE)=SPECTR(IAT)%MAT(:,:,:,IE) &
     &                               +WKPTL*AIMAG(G(I1:I2,I1:I2,:))/PI
          ENDDO
        ENDDO ! END OF LOOP OVER MATSUBARA FREQUENCIES
      ENDDO   ! END OF LOOP OVER KPOINTS
!
!     ==========================================================================
!     ==  WRITE SPECTRAL FUNCTION                                             ==
!     ==========================================================================
      ICALLCOUNT=ICALLCOUNT+1
      DO IAT=1,NAT
        WRITE(STRING,*)IAT
        FILE='_SPECTRALFUNCTION_AT'//TRIM(ADJUSTL(STRING))
        WRITE(STRING,*)ICALLCOUNT
        FILE=TRIM(ADJUSTL(FILE))//'_I'//TRIM(ADJUSTL(STRING))
        FILE=TRIM(ADJUSTL(FILE))//'.DAT'
        CALL FILEHANDLER$SETFILE('SPECTRALFILE',.TRUE.,-FILE)
        IF(ICALLCOUNT.EQ.1) THEN
          CALL FILEHANDLER$SETSPECIFICATION('SPECTRALFILE','STATUS','UNKNOWN')
          CALL FILEHANDLER$SETSPECIFICATION('SPECTRALFILE','ACTION','WRITE')
          CALL FILEHANDLER$SETSPECIFICATION('SPECTRALFILE','FORM','FORMATTED')
        END IF
        CALL FILEHANDLER$UNIT('SPECTRALFILE',NFIL)
        DO IE=1,NE
          NORB=ATOMSET(IAT)%ICHI2-ATOMSET(IAT)%ICHI1+1
          WRITE(NFIL,FMT='(200F15.5)')E(IE)*27.211D0 &
     &       ,(SPECTR(IAT)%MAT(I,I,1,IE),I=1,NORB)
        ENDDO
        CALL FILEHANDLER$CLOSE('SPECTRALFILE')
      ENDDO
STOP 'FORCED AFTER WRITING SPECTRAL FUNCTION'
!
!     ==========================================================================
!     ==  CLEAN                                                               ==
!     ==========================================================================
      DO IAT=1,NAT
        DEALLOCATE(SPECTR(IAT)%MAT)
      ENDDO
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_GLOCNI()
!     **************************************************************************
!     ** THIS ROUTINE IS DERIEVED FROM DMFT_GLOC                              **
!     ** IN CONTRAST TO THE LATTER IT DOES NOT INCLUDE THE SELF ENERGY.       **
!     ** THE RESULT IS STORED IN ATOMSET%GNI AND ATOMSET%GNILAUR
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
      COMPLEX(8)               :: G(NCHI,NCHI,NDIMD)
      REAL(8)                  :: WKPTL
      INTEGER(4)               :: IKPT,NU,IAT,I1,I2,IDIMD,ILAU
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_GLOCNI')
      IF(NLAU.GT.1) THEN
        CALL ERROR$MSG('NLAU EXCEEDS MAXIMUM VALUE OF 1')
        CALL ERROR$I4VAL('NLAU',NLAU)
        CALL ERROR$STOP('DMFT_GLOCNI')
      END IF
!
      DO IAT=1,NAT
        ATOMSET(IAT)%GNI=(0.D0,0.D0)
        ATOMSET(IAT)%GNILAUR=(0.D0,0.D0)
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
        MATX=-MU*KSET(IKPT)%SMAT+KSET(IKPT)%HRHO-KSET(IKPT)%GAMMA
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
            ATOMSET(IAT)%GNILAUR(:,:,:,ILAU)=ATOMSET(IAT)%GNILAUR(:,:,:,ILAU) &
     &                                     +WKPTL*GLAUR(I1:I2,I1:I2,:,ILAU)
          ENDDO
        ENDDO
!
        DO NU=1,NOMEGA
!         == CONSTRUCT LATTICE GREENS FUNCTION =================================
          MAT=(CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT
          MAT=MAT-KSET(IKPT)%HRHO+KSET(IKPT)%GAMMA
          CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
!         == ACCOUNT FOR K-POINT (-K) ==========================================
          IF(KSET(IKPT)%TADDMINUSK) THEN
!           == TRANSPOSE SMAT,HRHO,GAMMA, BUT NEITHER CI NOR SLOC ==============
            MAT=(-CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT
            MAT=MAT-KSET(IKPT)%HRHO+KSET(IKPT)%GAMMA
            DO IDIMD=1,NDIMD
              MAT(:,:,IDIMD)=TRANSPOSE(CONJG(MAT(:,:,IDIMD)))
            ENDDO
            CALL SPINOR$INVERT(NDIMD,NCHI,MAT,MAT2)
            G=0.5D0*(G+MAT2)
          END IF
!         == PROJECT GREENS FUNCTION ONTO INFIVIDUAL SITES =====================
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            ATOMSET(IAT)%GNI(:,:,:,NU)=ATOMSET(IAT)%GNI(:,:,:,NU) &
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
            ATOMSET(IAT)%GNILAUR(:,:,IDIMD,ILAU)=0.5D0 &
     &                                 *(ATOMSET(IAT)%GNILAUR(:,:,IDIMD,ILAU) &
     &                  +CONJG(TRANSPOSE(ATOMSET(IAT)%GNILAUR(:,:,IDIMD,ILAU))))
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
!CALL DMFT_SOLVERTEST(IAT)
!!$STOP
!IF(IAT.EQ.2)CALL DMFT_SOLVERTEST(IAT)
!CALL TESTRHOOFG(NLOC,NDIMD,NOMEGA,NLAU,KBT,ATOMSET(IAT)%GLOC,ATOMSET(IAT)%GLOCLAUR)
!       __SELF ENERGY IS RETURNED AS ATOMSET%DPHIDG_____________________________
        CALL DMFT_DYNAMICSOLVER(NLOC,NDIMD,NOMEGA,NLAU,KBT &
     &                         ,ATOMSET(IAT)%GLOC,ATOMSET(IAT)%GLOCLAUR &
     &                         ,ATOMSET(IAT)%U,PHILW &
     &                         ,ATOMSET(IAT)%DPHIDG,ATOMSET(IAT)%DPHIDGLAUR &
     &                         ,ATOMSET(IAT)%DEDU &
     &                         ,ATOMSET(IAT)%NATORB%PIPHI &
     &                         ,ATOMSET(IAT)%NATORB%CHIPHI)
        IF(SUM(ABS(ATOMSET(IAT)%DPHIDG(:,:,:,1))).GT.1.D-6) THEN
          CALL SPINOR_PRINTMATRIX(6,'DPHIDG(NU=1) IN SOLVER ' &
    &                      ,1,NLOC,NDIMD,NLOC,ATOMSET(IAT)%DPHIDG(:,:,:,1))
        END IF
        ETOT=ETOT+PHILW !SCREENING ALREADY CONTAINED IN U-TENSOR
        WRITE(*,FMT='(60("."),T1'&
     &        //',"--DMFT--LW FUNCTIONAL FOR DYNAMIC CORRELATION OF ATOM "' &
     &        //',I3,T60,":",F20.10)')IAT,PHILW
      ENDDO ! END OF LOOP OVER ATOMS (IAT)
      WRITE(*,FMT='(60("."),T1'&
     &        //',"--DMFT--LW FUNCTIONAL FOR DYNAMIC CORRELATION (TOTAL)"' &
     &        //',T60,":",F20.10)')ETOT
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVERTEST(IAT)
!     **************************************************************************
!     ** TEST NUMERICAL DERIVATIVE OF DYNAMICSOLVER                           **
!     **                                                                      **
!     **  SIGMA_AB(NU) = BETA*DPHI/DG_BA(NU)                                  **
!     **  DPHI=KBT*TRACE[SIGMA(NU)*DG(NU)+SIGMA(-NU)*DG(-NU)]                 **
!     **************************************************************************
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
      INTEGER(4)            :: I,IPM,NU
      COMPLEX(8)            :: CSVAR
!     **************************************************************************
!
!     ==========================================================================
!     == REFERENCE CALCULATION AND CALCULATION OF SIGMA                       ==
!     ==========================================================================
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
!     ==========================================================================
!     == DEFINE DISPLACEMENT DGLOC                                            ==
!     ==========================================================================
      ALLOCATE(DGLOC(NLOC,NLOC,NDIMD,NOMEGA))
      DGLOC=(0.D0,0.D0)
      DGLOC(1,1,1,1)=(1.D0,2.D0)
!
!     ==========================================================================
!     == ANALYTIC DERIVATIVE FOR A GIVEN DISPLACEMENT DGLOC*SCALE             ==
!     ==========================================================================
      Y0=0.D0
      DO NU=1,NOMEGA
        CALL SPINOR$TRACEAB(NDIMD,NLOC,SIGMA(:,:,:,NU),DGLOC(:,:,:,NU),CSVAR)
        Y0=Y0+2.D0*KBT*REAL(CSVAR)  ! FACTOR 2 INCLUDES NEGATIVE FREQUENCIES
      ENDDO

WRITE(*,FMT='("SOLVERTEST Y(",I3,")=",2F20.10)')0,Y0
!
!     ==========================================================================
!     == FINITE DISPLACEMENTS                                                 ==
!     ==========================================================================
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
        ENDDO  ! IPM: PLUS AND MINUS
WRITE(*,FMT='("SOLVERTEST Y(",I3,")=",2F20.10)')I,Y0
      ENDDO ! I: DISPLACEMENTS SIZE
STOP 'FORCED STOP IN DMFT_SOLVERTEST'
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
      INTEGER(4)             :: NU,I,J,IDIMD,ILAU
!     **************************************************************************
                                     CALL TRACE$PUSH('DMFT_DYNAMICSOLVER')
!
!     ==========================================================================
!     == CONVERT U-TENSOR                                                     ==
!     ==========================================================================
      U=(0.D0,0.D0)
      U(  :NLOC,  :NLOC,  :NLOC,  :NLOC)=CMPLX(UCHI,KIND=8)
      U(NLOC+1:,NLOC+1:,NLOC+1:,NLOC+1:)=CMPLX(UCHI,KIND=8)
      U(  :NLOC,NLOC+1:,  :NLOC,NLOC+1:)=CMPLX(UCHI,KIND=8)
      U(NLOC+1:,  :NLOC,NLOC+1:,  :NLOC)=CMPLX(UCHI,KIND=8)
!
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
      DO NU=1,NOMEGA
        S(:,:,NU)=MATMUL(CHIPHI,MATMUL(S(:,:,NU),CONJG(TRANSPOSE(CHIPHI))))
        CALL SPINOR$SHRINKDOWN(NDIMD,NLOC,S(:,:,NU),SLOC(:,:,:,NU))
      ENDDO
      DO I=1,NLAU
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
      USE STRINGS_MODULE
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
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: NU,I,J,K,L
      COMPLEX(8)            :: UCOPY(NORB,NORB,NORB,NORB) !U-TENSOR
!      CHARACTER(32),PARAMETER :: TYPE='LINEAR' ! 'NONE','LINEAR','STATIC'
!      CHARACTER(32),PARAMETER :: TYPE='NONE' ! 'NONE','LINEAR',
!      CHARACTER(32),PARAMETER :: TYPE='STATIC' ! 'NONE','LINEAR',
      CHARACTER(32),PARAMETER :: TYPE='2NDORDER' ! 'NONE','LINEAR',
      REAL(8)                 :: PI
      LOGICAL(4),PARAMETER    :: TMYSECONDORDER=.true.
!     **************************************************************************
                                     CALL TRACE$PUSH('DMFT_SOLVERIO')
      PI=4.D0*ATAN(1.D0)
      ETOT=0.D0
      S=(0.D0,0.D0)
      SLAUR=(0.D0,0.D0)
      DEDU=(0.D0,0.D0)
!
!     ==========================================================================
!     == SKIP IF U-TENSOR VANISHES                                            ==
!     ==========================================================================
      IF(MAXVAL(ABS(U)).LT.1.D-8) THEN
        PRINT*,'NO DYNAMIC LUTTINGER WARD FUNCTIONAL CALCULATED:'
        PRINT*,'U-TENSOR VANISHES'
        CALL TRACE$POP()
        RETURN 
      END IF
!
!     ==========================================================================
!     == EVALUATE DYNMICAL LUTTINGER WARD FUNCTIONAL AND SELF ENERGY          ==
!     ==========================================================================
      IF(TYPE.EQ.'NONE') THEN
        PRINT*,'USING STATIC LUTTINGER WARD FUNCTIONAL'
        PRINT*,'(NO DYNAMIC CORRECTION. FOR TESTING PURPOSES ONLY)'
!       == SOMETHING MUST BE THERE TO AVOID DIVIDE BY ZEROS ===================
        DO I=1,NORB
          S(I,I,1)=(1.D-8,1.D-8)
          ETOT=ETOT+2.D0*REAL(G(I,I,1)*S(I,I,1))
        ENDDO
        CALL TRACE$POP()
        RETURN
      ELSE IF(TYPE.EQ.'LINEAR') THEN
        PRINT*,'USING LINEAR LUTTINGER-WARD FUNCTIONAL'
        PRINT*,'(USE FOR TESTING PURPOSES ONLY)'
        CALL DMFT_SOLVER_LINEAR(NORB,NOMEGA,KBT,G,ETOT,S)  
        CALL TRACE$POP()
        RETURN
      ELSE IF(TYPE.EQ.'STATIC') THEN
        PRINT*,'USING HF AS LUTTINGER-WARD FUNCTIONAL'
        PRINT*,'(USE FOR TESTING PURPOSES ONLY)'
        CALL DMFT_SOLVER_STATIC(NORB,NOMEGA,KBT,G,U,ETOT,S)  
        CALL TRACE$POP()
        RETURN
      ELSE IF(TYPE.EQ.'2NDORDER') THEN
!!$UCOPY=U
!!$DO I=1,6
!!$  IF(I.EQ.1) J=1   ! S
!!$  IF(I.EQ.2) J=2   ! EG
!!$  IF(I.EQ.3) J=4   ! EG
!!$  IF(I.EQ.4) J=7   ! S  
!!$  IF(I.EQ.5) J=8
!!$  IF(I.EQ.6) J=10
!!$  UCOPY(J,:,:,:)=0.D0
!!$  UCOPY(:,J,:,:)=0.D0
!!$  UCOPY(:,:,J,:)=0.D0
!!$  UCOPY(:,:,:,J)=0.D0
!!$ENDDO
!!$        CALL MSIGMA$2NDORDER(NORB,NOMEGA,NLAU,KBT,G,UCOPY,ETOT,S,SLAUR)

!       =====================================================================
!       == STEFFEN ONLY USES THE DIAGRAM U_IJIJ BUT NOT U_IJJI  =============
!       == NEVERTHELESS THE RESULTS ARE SIMILAR AS WITH INCLUDING EXCHANGE ==
!       =====================================================================
        PRINT*,'ENTERING CALCULATION OF LW FUNCTIONAL'
        CALL MSIGMA$2NDORDER(NORB,NOMEGA,NLAU,KBT,G,U,ETOT,S,SLAUR)
        PRINT*,'CALCULATION OF LW FUNCTIONAL DONE'

        IF(TMYSECONDORDER) THEN
          CALL DMFT_SOLVERIO_WRITEGF('S',NORB,NOMEGA,KBT,S)
          PRINT*,'FROM STEFFEN: ETOT',ETOT
          PRINT*,'FROM STEFFEN: SIGMA',S(:,:,1)
!         == MY RESULTS HAVE A STRANGE BEHAVIOR IN THE SELF ENERGY FOR LARGE 
!         ==  MATSUBARA FREQUENCIES. PROBABLY REGULARIZATION NEEDED
!         PRINT*,'ENTERING CALCULATION OF LW FUNCTIONAL'
          CALL DMFT_SOLVER_2NDORDER(NORB,NOMEGA,NLAU,KBT,G,GLAUR,U,ETOT,S,SLAUR)  
          PRINT*,'CALCULATION OF LW FUNCTIONAL DONE'
          CALL DMFT_SOLVERIO_WRITEGF('G',NORB,NOMEGA,KBT,G)
          CALL DMFT_SOLVERIO_WRITEGF('S',NORB,NOMEGA,KBT,S)
          PRINT*,'MINE: ETOT',ETOT
          PRINT*,'MINE: SIGMA',S(:,:,1)
          STOP 'MINE'
        END IF
        CALL TRACE$POP()
        RETURN
      ELSE
        CALL ERROR$MSG('TYPE NOT RECOGNIZED')
        CALL ERROR$CHVAL('TYPE',TYPE)
        CALL ERROR$STOP('DMFT_SOLVERIO')
      END IF
                                     CALL TRACE$POP()
      RETURN
      END
!   
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVERIO_WRITEGF(ID,NORB,NOMEGA,KBT,GF)
!     **************************************************************************
!     ** WRITE GREENS FUNCTION OR SELF ENERGY TO FILE
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: NORB
      INTEGER(4)  ,INTENT(IN) :: NOMEGA
      REAL(8)     ,INTENT(IN) :: KBT
      COMPLEX(8)  ,INTENT(IN) :: GF(NORB,NORB,NOMEGA)
      LOGICAL(4)  ,SAVE       :: TINI=.FALSE.
      CHARACTER(32)           :: FILEID
      INTEGER(4)              :: NFIL
      INTEGER(4)              :: M,N
      REAL(8)                 :: OMEGA
      REAL(8)                 :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
!
!     ==========================================================================
!     == DEFINE FILES FOR FILEHANDLER                                         ==
!     ==========================================================================
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        FILEID=+'DMFT_GREENF_M'
        CALL FILEHANDLER$SETFILE(FILEID,.TRUE.,-'.DMFTGFIN_MINE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'STATUS','REPLACE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'FORM','FORMATTED')
        FILEID=+'DMFT_SIGMA_M'
        CALL FILEHANDLER$SETFILE(FILEID,.TRUE.,-'.DMFTSIGMA_MINE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'STATUS','REPLACE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'FORM','FORMATTED')
      END IF
!
!     ==========================================================================
!     == OPEN FILE AND DETERMINE CONVERSION FACTOR                            ==
!     ==========================================================================
      IF(ID.EQ.'G') THEN
        CALL FILEHANDLER$UNIT('DMFT_GREENF_M',NFIL)
      ELSE IF(ID.EQ.'S') THEN
        CALL FILEHANDLER$UNIT('DMFT_SIGMA_M',NFIL)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED (MUST BE "G" OR "S")')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('MSIGMA_WRITEGF')
      END IF
!
!     ==========================================================================
!     == WRITE FILE                                                           ==
!     ==========================================================================
      DO N=1,NOMEGA
        OMEGA=REAL(2*N-1)*PI*KBT
        WRITE(NFIL,'(ES24.16E3,5X)',ADVANCE='NO') OMEGA
        DO M=1,NORB
          WRITE(NFIL,'(ES24.16E3,3X,ES24.16E3,5X)',ADVANCE='NO') &
     &                      REAL( GF(M,M,N),KIND=8), AIMAG( GF(M,M,N))
        ENDDO !END LOOP OVER ORBITALS
        WRITE(NFIL,*)
      ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!
!     ==========================================================================
!     == CLOSE FILES                                                          ==
!     ==========================================================================
      CALL FILEHANDLER$CLOSE('DMFT_GREENF')
      CALL FILEHANDLER$CLOSE('DMFT_SIGMA')
      RETURN
      END SUBROUTINE DMFT_SOLVERIO_WRITEGF
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_EXPANDG(NORB,NOMEGA,NLAU,KBT,G,GLAU,NOMEGA2,GBIG)  
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB     !#(SPIN ORBITALS)
      INTEGER(4),INTENT(IN) :: NOMEGA   !#(POSITIVE MATSUBARA FREQUENCIES
      INTEGER(4),INTENT(IN) :: NOMEGA2   !#(POSITIVE MATSUBARA FREQUENCIES
      INTEGER(4),INTENT(IN) :: NLAU   !#(POSITIVE MATSUBARA FREQUENCIES
      REAL(8)   ,INTENT(IN) :: KBT      
      COMPLEX(8),INTENT(IN) :: G(NORB,NORB,NOMEGA)    !GREENS FUNCTION
      COMPLEX(8),INTENT(IN) :: GLAU(NORB,NORB,NLAU+1)    !GREENS FUNCTION
      COMPLEX(8),INTENT(OUT):: GBIG(NORB,NORB,NOMEGA)    !SELF ENERGY
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4) :: NU,ILAU
      REAL(8)    :: OMEGA,PI
      COMPLEX(8) :: CSVAR
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      GBIG(:,:,:NOMEGA)=G
      DO NU=NOMEGA+1,NOMEGA2
        GBIG(:,:,NU)=(0.D0,0.D0)
        OMEGA=REAL(2*NU-1)*PI*KBT
        CSVAR=1.D0/(CI*OMEGA)
        DO ILAU=1,NLAU+1
          GBIG(:,:,NU)=GBIG(:,:,NU)+GLAU(:,:,ILAU)*CSVAR**ILAU
        ENDDO
      ENDDO
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
      INTEGER(4)            :: NU,I,J
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)   ,PARAMETER  :: A=2.D0*7.D-3,B=1.D-3,LAMBDA=1.D-2
!      REAL(8)   ,PARAMETER  :: A=0.D0,B=1.D-1,LAMBDA=1.D-1
!      REAL(8)   ,PARAMETER  :: A=-7.D-3,B=1.D0,LAMBDA=1.D-1
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      S=(0.D0,0.D0)
      DO NU=1,NOMEGA
        OMEGA=REAL(2*NU-1,KIND=8)*PI*KBT
!       == DEFINE SELF ENERGY ==================================================
        DO I=1,NORB
          S(I,I,NU)=EXP(-LAMBDA/ABS(OMEGA))*(A/(CI*OMEGA)+B*(1.D0,0.D0))
        ENDDO
      ENDDO
!
!     == CALCULATE ENERGY=======================================================
      ETOT=0.D0
      DO I=1,NORB
        DO J=1,NORB
!         ==  FACTOR 2 FROM -OMEGA_NU 
          ETOT=ETOT+2.D0*KBT*REAL(SUM(S(I,J,:)*G(J,I,:)))
        ENDDO
      ENDDO  
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVER_STATIC(NORB,NOMEGA,KBT,G,U,ETOT,S)  
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
      COMPLEX(8),INTENT(IN) :: U(NORB,NORB,NORB,NORB) !U-TENSOR
      REAL(8)   ,INTENT(OUT):: ETOT                   !LUTTINGER WARD FUNCTIONAL
      COMPLEX(8),INTENT(OUT):: S(NORB,NORB,NOMEGA)    !SELF ENERGY
      INTEGER(4),PARAMETER  :: NLAU=0
      REAL(8)               :: PI
      REAL(8)               :: OMEGA(NOMEGA)
      INTEGER(4)            :: NU,I,J,K,L
      COMPLEX(8)            :: RHO(NORB,NORB)
      COMPLEX(8)            :: HAM(NORB,NORB)
      COMPLEX(8)            :: CSVAR
      REAL(8)               :: FN(NLAU+1)      
      REAL(8)               :: SVAR
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      S=(0.D0,0.D0)
      ETOT=0.D0
!
!     ==========================================================================
!     == CALCULATE DENSITY MATRIX 
!     ==========================================================================
      RHO=(0.D0,0.D0)
      DO NU=1,NOMEGA
        OMEGA(NU)=REAL(2*NU-1,KIND=8)*PI*KBT
        RHO=RHO+KBT*G(:,:,NU)
      ENDDO
      RHO=RHO+TRANSPOSE(CONJG(RHO))
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,NLAU+1,FN)
      DO I=1,NLAU+1
        DO J=1,NORB
          RHO(J,J)=RHO(J,J)+FN(I)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == PRINT DENSITY MATRIX FOR TESTING
!     ==========================================================================
      WRITE(*,FMT='(100("="),T30,"  ",A,"  ")')'DENSITY MATRIX IN SOLVER_STATIC'
      DO I=1,NORB
        WRITE(*,FMT='(100("(",2F10.5,")"))')RHO(I,:)
      ENDDO
      SVAR=0.D0
      DO I=1,NORB
        SVAR=SVAR+REAL(RHO(I,I))
      ENDDO
      WRITE(*,FMT='("NUMBER OF ELECTRONS IN G: ",F10.5)')SVAR
!
!     ==========================================================================
!     == CALCULATE HARTREE FOCK POTENTIAL                                     ==
!     ==========================================================================
      ETOT=0.D0
      HAM=(0.D0,0.D0)
      DO I=1,NORB
        DO J=1,NORB
          DO K=1,NORB
            DO L=1,NORB
              CSVAR=-0.5D0*U(I,J,K,L) 
              IF(ABS(CSVAR).EQ.0.D0) CYCLE
              ETOT=ETOT+REAL(CSVAR*RHO(K,J)*RHO(L,I))
              HAM(K,J)=HAM(K,J)+CSVAR*RHO(I,L)
              HAM(I,L)=HAM(I,L)+CSVAR*RHO(K,J)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
PRINT*,'ETOT',ETOT
PRINT*,'HAM',HAM
!
!     ==========================================================================
!     == CALCULATE SELF ENERGY
!     ==========================================================================
      DO NU=1,NOMEGA
        S(:,:,NU)=HAM(:,:)
      ENDDO      
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SOLVER_2NDORDER(NORB,NOMEGA,NLAU,KBT,G,GLAU,UIN &
     &                               ,ETOT,SIGMA,SIGMALAU)  
!     **************************************************************************
!     ** MIMICKS A SOLVER FOR THE LUTTINGER WARD FUNCTIONAL AND SELF ENERGY   **
!     ** ON THE BASIS OF A LINEAR LUTTINGER WARD FUNCTIONAL, WHICH IS         **
!     ** DEFINED BY A FROZEN SELF ENERGY.                                     **
!     **                                                                      **
!     ** SIGMA_{A,B}(NU)=BETA*DPHI_LW/DG_{B,A}(NU)                            **
!     **                                                                      **
!     **  Y(I1O1 I2O2)=SUM_O3I3 U(O1O3 I1I3)*KBT*SUM_NU GNU(O2I3)*GNU(O3I2)   **
!     **  PHI=-0.25*KBT*SUM_{I1O1 I2O2} Y(I1O1 I2O2)*Y(I2O2,I1O1              **
!     **                                                                      **
!     ** NO REGULARIZATION!                                                   **
!     **   RESTRICT BOSONIC MATSUBARA FREQUENCIES TO AVOID INCREASING TAIL!   **
!     **                                                                      **
!     **   TAILS NOT COMPLETED!!!!!!                                          **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************PETER BLOECHL GOSLAR 2015***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB     !#(SPIN ORBITALS)
      INTEGER(4),INTENT(IN) :: NOMEGA   !#(POSITIVE MATSUBARA FREQUENCIES
      INTEGER(4),INTENT(IN) :: NLAU     !#(LAURENT COEFFICIENTS-1)
      REAL(8)   ,INTENT(IN) :: KBT      
      COMPLEX(8),INTENT(IN) :: G(NORB,NORB,NOMEGA)    !GREENS FUNCTION
      COMPLEX(8),INTENT(IN) :: GLAU(NORB,NORB,NLAU+1)    !GREENS FUNCTION
      COMPLEX(8),INTENT(IN) :: UIN(NORB,NORB,NORB,NORB) !U-TENSOR
      REAL(8)   ,INTENT(OUT):: ETOT                    !LUTTINGER WARD FUNCTIONAL
      COMPLEX(8),INTENT(OUT):: SIGMA(NORB,NORB,NOMEGA) !SELF ENERGY
      COMPLEX(8),INTENT(OUT):: SIGMALAU(NORB,NORB,NLAU) !SELF ENERGY
      LOGICAL(4),PARAMETER  :: TINCLUDEEXCHANGE=.TRUE. !INCLUDE EXCHANGE DIAGRAMS
      LOGICAL(4),PARAMETER  :: TDIAGG=.TRUE.   ! USE ONLY DIAGONAL OF GREENSF
      INTEGER(4),PARAMETER  :: NUBX=100   ! CAN BE UP TO 2*NOMEGA
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)               :: PI
      INTEGER(4)            :: NUB,NU1,NU2
      INTEGER(4)            :: ILAU
      INTEGER(4)            :: I1,O1,I2,O2,I3,O3
      COMPLEX(8)            :: POLE(NUBX,NLAU)
      COMPLEX(8)            :: U(NORB,NORB,NORB,NORB)
      COMPLEX(8)            :: Y(NORB,NORB,NORB,NORB)
      COMPLEX(8)            :: YU(NORB,NORB,NORB,NORB)
      COMPLEX(8)            :: CSVAR
      REAL(8)               :: SVAR,EPREV
      INTEGER(4)            :: ISVAR1,ISVAR2
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      SIGMA=(0.D0,0.D0)
      SIGMALAU=(0.D0,0.D0)
      ETOT=0.D0
!
!     ==========================================================================
!     == ANTISYMMETRIZE U-TENSOR TO INCLUDE EXCHANGE DIAGRAMS                 ==
!     ==========================================================================
      IF(TINCLUDEEXCHANGE) THEN
        DO I2=1,NORB
          DO I1=1,NORB
            DO O2=1,NORB
              DO O1=1,NORB
                U(O1,O2,I1,I2)=0.25D0*(UIN(O1,O2,I1,I2)-UIN(O1,O2,I2,I1) &
     &                                -UIN(O2,O1,I1,I2)+UIN(O2,O1,I2,I1))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        U=UIN
      END IF
!
!     ==========================================================================
!     == PREPARE TAILS OF GREENS FUNCTION                                     ==
!     ==========================================================================
      DO NU1=NOMEGA+1,NOMEGA+NUBX
        CSVAR=1.D0/(CI*REAL(2*NU1-1)*PI*KBT)  !1/(I*OMEGA)
        DO ILAU=1,NLAU
          POLE(NU1-NOMEGA,ILAU)=CSVAR**ILAU
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == BOSONIC MATSUBARA SUM
!     ==========================================================================
      DO NUB=0,NUBX
        Y(:,:,:,:)=(0.D0,0.D0)
        DO I2=1,NORB
          DO O2=1,NORB
            DO I3=1,NORB
              DO O3=1,NORB
!               == IMPLEMENT DIAGONAL APPROXIMATION ============================ 
                IF(TDIAGG) THEN
                  IF(I3.NE.O2) CYCLE
                  IF(I2.NE.O3) CYCLE
                END IF
!
!               == SKIP IF U-TENSOR ELEMENTS VANISH ============================
                IF(MAXVAL(ABS(U(:,O3,:,I3))).LT.1.D-8) CYCLE
!               ================================================================
!               == PERFORM FERMIONIC MATSUBARA SUM                            ==
!               == DO NU2=-NOMEGA+1,NOMEGA-NUB
!               ==   NU1=NU2+NUB
!               ==   CSVAR=CSVAR+G(O2,I3,NU1)*G(O3,I2,NU2)
!               == ENDDO
!               ================================================================
                CSVAR=(0.D0,0.D0)
!
!               == NU2,NU1=NU2+NUB NEGATIVE OR ZERO ============================
                DO NU2=-NOMEGA-NUB+1,-NOMEGA+1
                  NU1=NU2+NUB
                  DO ILAU=1,NLAU+1
                    CSVAR=CSVAR+CONJG(G(I3,O2,-NU1+1) &
        &                            *GLAU(I2,O3,ILAU)*POLE(-NU2+1-NOMEGA,ILAU))
                  ENDDO
                ENDDO
                DO NU2=-NOMEGA+1,-NUB
                  NU1=NU2+NUB
                  CSVAR=CSVAR+CONJG(G(I3,O2,-NU1+1)*G(I2,O3,-NU2+1))
                ENDDO
!
!               == NU2 NEGATIVE, NU1=NU2+NUB POSITIVE ==========================
!               == NU2=-NUB+1,...,0
                ISVAR1=MAX(-NUB+1,-NOMEGA+1)
                ISVAR2=MIN(0,NOMEGA-NUB)
                DO NU2=ISVAR1,ISVAR2
                  NU1=NU2+NUB
                  CSVAR=CSVAR+G(O2,I3,NU1)*CONJG(G(I2,O3,-NU2+1))
                ENDDO
!
!               == NU2,NU1+MUB POSITIVE ========================================
                DO NU2=1,NOMEGA-NUB
                  NU1=NU2+NUB
                  CSVAR=CSVAR+G(O2,I3,NU1)*G(O3,I2,NU2)
                ENDDO

                DO NU2=NOMEGA-NUB+1,NOMEGA
                  NU1=NU2+NUB
                  DO ILAU=1,NLAU+1
                    CSVAR=CSVAR+GLAU(O2,I3,ILAU)*POLE(NU1-NOMEGA,ILAU) &
        &                      *G(O3,I2,NU2)
                  ENDDO
                ENDDO

!               ================================================================
                CSVAR=CSVAR*KBT
                DO I1=1,NORB
                  DO O1=1,NORB
                    Y(I1,O1,I2,O2)=Y(I1,O1,I2,O2)+U(O1,O3,I1,I3)*CSVAR
                  ENDDO 
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
!       ========================================================================
!       == ADD TO LUTTINGER WARD FUNCTIONAL                                   ==
!       ========================================================================
!       == NOTE THAT ETOT IS THE LW FUNCTIONAL AND NOT THE INTERACTION ENERGY!
        SVAR=-0.25*KBT
        IF(NUB.NE.0) SVAR=SVAR*2.D0  !ACCOUNT FOR NEGATIVE BOSONIC FREQUENCIES
!EPREV=ETOT
        DO I1=1,NORB
          DO I2=1,NORB
            DO O1=1,NORB
              DO O2=1,NORB
                ETOT=ETOT+SVAR*REAL(Y(I1,O1,I2,O2)*Y(I2,O2,I1,O1))
              ENDDO
            ENDDO
          ENDDO
        ENDDO 
!PRINT*,'ETOT ',ETOT,ETOT-EPREV,NUB
!
!       ========================================================================
!       == CALCULATE SELF ENERGY                                              ==
!       ========================================================================
!       == PERFORM MULTIPLICATION WITH U-TENSOR ================================
        YU=(0.D0,0D0)
        DO I1=1,NORB
          DO O1=1,NORB
            DO I2=1,NORB
              DO O2=1,NORB
                DO I3=1,NORB
                  DO O3=1,NORB
                    YU(I2,O2,I3,O3)=YU(I2,O2,I3,O3)+Y(I2,O2,I1,O1)*U(O1,O3,I1,I3)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO 
!       == FACTOR -0.25*KBT PREFACTOR IN FRONT OF PHI=-/(2BETA)\SUM_NUB TR[Y^2]
!       == FACTOR 4 BECAUSE OF SECOND ORDER (4 GREENS FUNCTION) 
!       == FINAL 1/BETA FROM MATSUBARA SUM DROPS WITH DEFINITION
        YU=-KBT*YU
!
!       == CALCULATE SELF ENERGY ===============================================
        DO I2=1,NORB
          DO O2=1,NORB
            DO I3=1,NORB
              DO O3=1,NORB
!               == IMPLEMENT DIAGONAL APPROXIMATION ============================ 
                IF(TDIAGG) THEN
                  IF(I3.NE.O2) CYCLE
                  IF(I2.NE.O3) CYCLE
                END IF
!
!               == SKIP IF PREFACTOR VANISHES ===== ============================
                CSVAR=YU(I2,O2,I3,O3)
                IF(ABS(CSVAR).LT.1.D-8) CYCLE
!               ================================================================
!               == PERFORM FERMIONIC MATSUBARA SUM                            ==
!               == DO NU2=-NOMEGA+1,NOMEGA-NUB
!               ==   NU1=NU2+NUB
!               ==   CSVAR=CSVAR+G(O2,I3,NU1)*G(O3,I2,NU2)
!               == ENDDO
!               ================================================================
!               == NUB, NU1=NU2+NUB POSITIVE
                DO NU2=1,NOMEGA-NUB
                  NU1=NU2+NUB
!                 == SIGMA(A,B)=BETA*DPHI/DG(B,A)
                  SIGMA(I2,O3,NU2)=SIGMA(I2,O3,NU2)+CSVAR*G(O2,I3,NU1)
                ENDDO
!               == PATCH UP WITH TAILS =========================================
                DO NU2=NOMEGA-NUB+1,NOMEGA
                  NU1=NU2+NUB
!                 == SIGMA(A,B)=BETA*DPHI/DG(B,A)
                  DO ILAU=1,NLAU+1
                    SIGMA(I2,O3,NU2)=SIGMA(I2,O3,NU2) &
           &                        +CSVAR*GLAU(O2,I3,ILAU)*POLE(NU1-NOMEGA,ILAU)
                  ENDDO
                ENDDO
!               
!               == NUB NEGATIVE ================================================
                IF(NUB.NE.0) THEN  ! DO NOT COUNT NUB=0 TWICE!
!                 == NUB NEGATIVE, NU1=NU2+NUB NEGATIVE ========================
                  ISVAR1=MAX(1,-NOMEGA+1+NUB)
                  ISVAR2=MIN(NUB,NOMEGA)
!!$                  DO NU2=1,ISVAR1-1
!!$                    NU1=NU2-NUB
!!$                    SIGMA(I2,O3,NU2)=SIGMA(I2,O3,NU2) &
!!$    &                                   +CONJG(CSVAR*GTAIL(I3,O2,-NU1+1-NOMEGA))
!!$                  ENDDO
                  DO NU2=ISVAR1,ISVAR2
                    NU1=NU2-NUB
                    SIGMA(I2,O3,NU2)=SIGMA(I2,O3,NU2) &
    &                                              +CONJG(CSVAR*G(I3,O2,-NU1+1))
                  ENDDO
!                 == NUB NEGATIVE, NU1=NU2+NUB POSITIVE
                  DO NU2=NUB+1,NOMEGA
                    NU1=NU2-NUB
                    SIGMA(I2,O3,NU2)=SIGMA(I2,O3,NU2)+CONJG(CSVAR)*G(O2,I3,NU1)
                  ENDDO
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == FACTOR 2 FROM INCLUSION OF EXCHANGE DIAGRAMS                         ==
!     ==========================================================================
      IF(TINCLUDEEXCHANGE) THEN
        ETOT=2.D0*ETOT
        SIGMA=2.D0*SIGMA
      END IF
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
     &                      ,OMEGA,KBT,MU,KSET,ATOMSET
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER     :: NITER1=1000     ! X#(CONJUGATE GRADIENT STEPS)
      INTEGER(4),PARAMETER     :: NITER2=1000    ! X#(STEPS) FOR LINE SARCH
!     == TOLERANCE FOR LINE SEARCH MUST BE MUCH SMALLER THAN THAT FOR THE 
!     == OUTER LOOP. IF THE TOLERANCE IS TOO SMALL IT MAY GET STUCK
!     == BECAUSE OF NOISE BEING DOMINANT ON THUS SCALE.
      REAL(8)   ,PARAMETER     :: TOL1=1.D-4     ! TOLERANCE FOR CG ITERATION
      REAL(8)   ,PARAMETER     :: TOL2=1.D-10     ! TOLERANCE FOR LINE SEARCH
      REAL(8)   ,PARAMETER     :: AMIX=1.D-3     ! INITIAL MIXING
      LOGICAL(4),PARAMETER     :: TPRDRHO=.FALSE.! PRINT DEV OF DENSITY MATRIX
      LOGICAL(4),PARAMETER     :: TPRLINE=.TRUE.! PRINT LINE-SEARCH PROGRESS
      INTEGER(4),PARAMETER     :: NMEMX=20
      INTEGER(4)               :: NMEM
      COMPLEX(8)               :: SK(NCHI,NCHI,NDIMD,NMEMX)
      COMPLEX(8)               :: YK(NCHI,NCHI,NDIMD,NMEMX)
      COMPLEX(8)               :: DGAMMA(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GAMMA0(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GAMMAM(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: GAMMAP(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: FORCE0(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: FORCEM(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: DRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: CSVAR
      LOGICAL(4)               :: CONVG
      REAL(8)                  :: S0,SM,SP
      REAL(8)                  :: F0,FM,FP
      REAL(8)                  :: Y0,YM,YP,YM1
      REAL(8)                  :: D2YDS2
      REAL(8)                  :: DLEN
      INTEGER(4)               :: IKPT
      INTEGER(4)               :: ITER1,ITER2
      REAL(8)                  :: FN(NLAU+1)
      REAL(8)                  :: SVAR
      REAL(8)                  :: PI
      INTEGER(4)               :: I,IDIMD
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_CONSTRAINTS')
      PI=4.D0*ATAN(1.D0)
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,NLAU+1,FN)
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
!     ==  DEVIATION FROM TARGET DENSITY MATRIX                                ==
!     ==========================================================================
      DO IKPT=1,NKPTL
!
!       == THIS SYMMETRIZATION IS NOT REQUIRED. IT IS INTRODUCED TO AVOID
!       == POTENTIAL PROBLEMS DURING DEBUGGING
        CALL SPINOR$CONJUGATE(NDIMD,NCHI,KSET(IKPT)%GAMMA,MAT1)
        KSET(IKPT)%GAMMA=0.5D0*(KSET(IKPT)%GAMMA+MAT1)
!       
!       ========================================================================
!       ==  ADJUST H0 TO OBEY DENSITY MATRIX CONSTRAINT                       ==
!       ==  USE CONJUGATE GRADIENT METHOD TO MINIMIZE THE ABSOLUTE SQUARE     ==
!       ==  OF THE DEVIATION OF THE DENSITY MATRIX FROM THE TARGET            ==
!       ========================================================================
!       == INITIATE CONJUGATE GRADIENT ITERATION ===============================
!       == DGAMMA=DMAXDEV/DGAMMA ===========================================
        GAMMAM=KSET(IKPT)%GAMMA
        CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
     &              ,GAMMAM,YM,FORCEM)
        FORCEM=0.5D0*FORCEM  ! CONVERT TO WIRTINGER DERIVATIVE
        FORCEM=-FORCEM       ! DOWNHILL DIRECTION
        GAMMA0=GAMMAM+AMIX*FORCEM   ! FIRST SEARCH DIRECTION IS DOWNHILL
!
!IF(YM.GT.100.D0) THEN
!!$ CALL PLOTDOS(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
!!$&                    ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
!!$&                    ,GAMMAM)
!!$  STOP 'FORCED'
!END IF

        CONVG=.FALSE.
        NMEM=0
        DO ITER1=1,NITER1
          CONVG=.FALSE.
!
!         ======================================================================
!         == LINE SEARCH                                                      ==
!         ======================================================================
          DGAMMA=GAMMA0-GAMMAM
          CALL SPINOR$TRACEAB(NDIMD,NCHI,DGAMMA,DGAMMA,CSVAR)
          DLEN=SQRT(2.D0*REAL(CSVAR))
          IF(DLEN.EQ.0.D0) THEN
            CONVG=.TRUE.
            EXIT
!!$            CALL ERROR$MSG('NO SEARCH DIRECTION PROVIDED')
!!$            CALL ERROR$I4VAL('ITER1 ',ITER1)
!!$            CALL ERROR$STOP('DMFT_CONSTRAINTS')
          END IF
          DGAMMA=DGAMMA/DLEN
          SM=0.D0
          YM1=YM  ! PREVIOUS VALUE IN LINE SEARCH
          CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCEM,DGAMMA,CSVAR)
          FM=2.D0*REAL(CSVAR,KIND=8)
          IF(TPRLINE)PRINT*,'LINESEARCH',SM,YM,FM,0.D0
 
          CONVG=(ABS(FM).LT.1.D-8) 
          IF(CONVG) THEN
            PRINT*,'FORCE VANISHES: EXIT CONSTRAINT LOOP ITER1=',ITER1
            EXIT
          END IF
!
          YM1=YM
          S0=SM+AMIX*FM
          DO ITER2=1,NITER2 ! LINE SEARCH
            GAMMA0=GAMMAM+DGAMMA*S0
            CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                  ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
     &                  ,GAMMA0,Y0,FORCE0)
            FORCE0=0.5D0*FORCE0  ! CONVERT TO WIRTINGER DERIVATIVE
            FORCE0=-FORCE0       ! DOWNHILL DIRECTION
            CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCE0,DGAMMA,CSVAR)
            F0=2.D0*REAL(CSVAR,KIND=8)     ! DY/DS
            D2YDS2=-(F0-FM)/(S0-SM)   ! D^2Y/DS^2
!
!           == KEEP EACH SUCCESSFUL VALUE OF GAMMA =============================
            IF(Y0.LT.YM) KSET(IKPT)%GAMMA=GAMMA0
!
!           ====================================================================
!           == CHECK CONVERGENCE ===============================================
!           ====================================================================
            IF(TPRLINE)PRINT*,'LINESEARCH',S0,Y0,F0,D2YDS2
!IF(ITER2.GT.100)WRITE(*,'("LINESEARCH",I10,5E20.5)')ITER2,S0,Y0,F0,D2YDS2
            CONVG=.FALSE.
!           ====== TRANSLATE FORCE INTO DEVIATION IN DENSITY ===================
!           == DELTA-RHO=BETA*E
!           == (DELTA_RHO)**2=BETA**2*E**2
!           == F=D(DELTA-RHO)**2/DE=2*BETA**2*E
!           == E=F/(2*BETA**2)
!           == DELTA-RHO=F/(2*BETA)=KBT*F/2<TOL
            CONVG=(0.5D0*KBT*ABS(F0).LT.MIN(0.5D0*TOL1,TOL2)) 
!           CONVG=CONVG.OR.(ABS(S0-SM).LT.1.D-8)
!            CONVG=CONVG.OR.(SQRT(Y0).LT.TOL2) 
            IF(CONVG) EXIT
!
!           ====================================================================
!           == DETERMINE SIZE OF NEXT STEP =====================================
!           ====================================================================
            IF(Y0.GT.YM1) THEN 
!             == VALUE OF YM IN INCREASED. =====================================
!             == UNDO PREVIOUS STEP AND RETRY WITH SMALLER STEP SIZE ===========
              IF(FM*F0.LE.0.D0) THEN
                S0=(FM*S0-F0*SM)/(FM-F0)
                IF(TPRLINE)PRINT*,'INTERRUPT OK',Y0-YM1
              ELSE
                S0=0.5D0*(S0+SM) !FORLATER: -(S0-SM)+(S0+SM)/2 = -S0/2-SM*3/2
                IF(TPRLINE)PRINT*,'INTERRUPT BISECT',Y0-YM1
              END IF
              CYCLE
            END IF
            IF(D2YDS2.GT.0.D0) THEN  ! NEWTON-RAPHSON FOR POSITIVE CURVARTURE
              SP=(FM*S0-F0*SM)/(FM-F0)
              IF(ABS(SP-S0).GT.KBT) THEN
                IF(TPRLINE)PRINT*,'STEP SIZE LIMITED',KBT
                SP=S0+KBT*(SP-S0)/ABS(SP-S0)
              END IF
            ELSE ! SWITCH TO STEPPING MODE, IF CURVATURE IS NEGATIVE ========
              SP=S0+KBT*SIGN(1.D0,F0)
            END IF
!           == SWITCH ==========================================================
            SM=S0
            YM1=Y0
            FM=F0
            S0=SP
          ENDDO ! END OF LINE SEARCH
          IF(.NOT.CONVG) THEN
            CALL ERROR$MSG('LINE SEARCH NOT CONVERGED')
            CALL ERROR$R8VAL('CONV CRITERIUM',Y0)
            CALL ERROR$R8VAL('TOLERANCE',TOL2)
            CALL ERROR$STOP('DMFT_CONSTRAINTS')
          END IF
!
!         ======================================================================
!         == CHECK CONVERGENCE                                                ==
!         ======================================================================
PRINT*,'IKPT=',IKPT,' ITER=',ITER1,' |DEV^2|=',SQRT(Y0),' NITER2 ',ITER2
          CONVG=SQRT(Y0).LT.TOL1
          IF(CONVG) EXIT
!
!         ======================================================================
!         == PREPARE NEXT SEARCH DIRECTION DGAMMA                             ==
!         ======================================================================
          CALL DMFT_LBFGS(NCHI,NDIMD,NMEMX,NMEM,AMIX,YK,SK &
      &                  ,GAMMA0,GAMMAM,FORCE0,FORCEM,GAMMAP)
!
!!$CALL DMFT_SPINORANGLE(NCHI,NDIMD,GAMMAP-GAMMA0,GAMMA0-GAMMAM,SVAR)
!!$PRINT*,'ANGLE BETWEEN OLD AND NEW SEARCH DIRECTION   ',SVAR/PI*180.D0
!!$CALL DMFT_SPINORANGLE(NCHI,NDIMD,GAMMA0-GAMMAM,FORCE0,SVAR)
!!$PRINT*,'ANGLE BETWEEN FORCE AND OLD SEARCH DIRECTION ',SVAR/PI*180.D0
!!$CALL DMFT_SPINORANGLE(NCHI,NDIMD,GAMMAP-GAMMA0,FORCE0,SVAR)
!!$PRINT*,'ANGLE BETWEEN FORCE AND NEW SEARCH DIRECTION ',SVAR/PI*180.D0
!!$CALL DMFT_SPINORANGLE(NCHI,NDIMD,FORCEM,FORCE0,SVAR)
!!$PRINT*,'ANGLE BETWEEN OLD AND CURRENT FORCE          ',SVAR/PI*180.D0
!!$DO I=1,NMEM
!!$  CALL DMFT_SPINORANGLE(NCHI,NDIMD,GAMMAP-GAMMA0,SK(:,:,:,I),SVAR)
!!$  PRINT*,'ANGLE BETWEEN SEARCH DIRECTIONS  ',I,SVAR/PI*180.D0
!!$ENDDO
!
!         __ NEXT SEARCH DIRECTION______________________________________________
          GAMMAM=GAMMA0
          FORCEM=FORCE0
          YM=Y0
          GAMMA0=GAMMAP
        ENDDO ! END OF CONJUGATE GRADIENT LOOP
        IF(.NOT.CONVG) THEN
          CALL ERROR$MSG('LOOP NOT CONVERGED')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('ITER1',ITER1)
          CALL ERROR$R8VAL('MAX. DEVIATION',Y0)
          CALL ERROR$R8VAL('TOLERANCE',TOL1)
          CALL ERROR$STOP('DMFT_CONSTRAINTS')
        END IF
!
!       ========================================================================
!       ==  PRINT THE DEVIATION OF RHO FROM TARGET
!       ========================================================================
        IF(TPRDRHO) THEN
          PRINT*,'VIOLATION FOR DENSITY MATRIX CONSTRAINT FOR K-POINT ',IKPT
          CALL DMFT_PRINTDRHO(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
     &              ,KSET(IKPT)%GAMMA)
        END IF
      ENDDO   !END OF LOOP OVER K-POINTS
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_CONSTRAINTS_A()
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
     &                      ,OMEGA,KBT,MU,KSET,ATOMSET
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER     :: NITER1=1000     ! X#(CONJUGATE GRADIENT STEPS)
      INTEGER(4),PARAMETER     :: NITER2=10000   ! X#(STEPS) FOR LINE SARCH
      REAL(8)   ,PARAMETER     :: TOL1=1.D-3     ! TOLERANCE FOR CG ITERATION
      REAL(8)   ,PARAMETER     :: TOL2=1.D-8     ! TOLERANCE FOR LINE SEARCH
      REAL(8)   ,PARAMETER     :: SCALEINI=1.D-3 ! INITIAL MIXING
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0) ! SQRT(-1)
      CHARACTER(2),PARAMETER   :: CGID='PR'      ! CAN BE 'FR', 'PR' 
      LOGICAL(4),PARAMETER     :: TPRLINE=.FALSE. ! PRINT LINE-SEARCH PROGRESS
      LOGICAL(4),PARAMETER     :: TPRDRHO=.TRUE. ! PRINT DEV OF DENSITY MATRIX
      LOGICAL(4),PARAMETER     :: TPRFILE=.FALSE. ! PRINT DMFTCONSTRAINTS FILE
      COMPLEX(8)               :: DGAMMA(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: FORCE(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: FORCEM(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: CSVAR
      LOGICAL(4)               :: CONVG
      REAL(8)                  :: MAXDEV,MAXDEV0
      REAL(8)                  :: S0,SM,SP
      REAL(8)                  :: F0,FM,FP
      REAL(8)                  :: Y0,YM,YP
      REAL(8)                  :: D2YDS2
      REAL(8)                  :: FMFM
      INTEGER(4)               :: IKPT,IDIMD
      INTEGER(4)               :: ITER1,ITER2,ITER,I,J
      REAL(8)                  :: FN(NLAU+1)
      REAL(8)                  :: BETA
      REAL(8)                  :: SVAR,SVAR1,SVAR2
      REAL(8)                  :: COSALPHA,COSBETA,COSGAMMA
      CHARACTER(128)           :: STRING
      REAL(8)                  :: DLEN,FLEN
      REAL(8)                  :: PI
      INTEGER(4),SAVE          :: COUNTER=0
      CHARACTER(32),PARAMETER  :: FILEID='DMFTCLOOP'
      CHARACTER(128)           :: FILE
      INTEGER(4)               :: NFIL
      INTEGER(4),PARAMETER     :: NDIFF=8
      INTEGER(4),PARAMETER     :: NX=100
      INTEGER(4)               :: N
      REAL(8)                  :: DY(NDIFF+1)
      REAL(8)                  :: S
      COMPLEX(8)               :: Z(NDIFF)
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_CONSTRAINTS')
      PI=4.D0*ATAN(1.D0)
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,NLAU+1,FN)
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

      IF(TPRFILE) THEN
        COUNTER=COUNTER+1
        WRITE(FILE,*)COUNTER
        FILE='DMFTCONSTRAINTS_'//TRIM(ADJUSTL(FILE))//'.DAT'
        CALL FILEHANDLER$SETFILE(FILEID,.FALSE.,-FILE)
        IF(COUNTER.EQ.1) THEN
          CALL FILEHANDLER$SETSPECIFICATION(FILEID,'STATUS','REPLACE')
          CALL FILEHANDLER$SETSPECIFICATION(FILEID,'POSITION','REWIND')
          CALL FILEHANDLER$SETSPECIFICATION(FILEID,'ACTION','WRITE')
          CALL FILEHANDLER$SETSPECIFICATION(FILEID,'FORM','FORMATTED')
        END IF
        CALL FILEHANDLER$UNIT(FILEID,NFIL)
        REWIND(NFIL)
      END IF
!
!     ==========================================================================
!     ==  DEVIATION FROM TARGET DENSITY MATRIX                                ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        IF(TPRFILE) THEN
          WRITE(NFIL,*)'# ,ITER1,Y0,FLEN/DLEN,D2YDS2/DLEN**2' &
       &               //',S0,DLEN,COSALPHA,COSGAMMA,IKPT'
        END IF
!
!       == THIS SYMMETRIZATION IS NOT REQUIRED. IT IS INTRODUCED TO AVOID
!       == POTENTIAL PROBLEMS DURING DEBUGGING
        CALL SPINOR$CONJUGATE(NDIMD,NCHI,KSET(IKPT)%GAMMA,MAT1)
        KSET(IKPT)%GAMMA=0.5D0*(KSET(IKPT)%GAMMA+MAT1)
!       
!       ========================================================================
!       ==  ADJUST H0 TO OBEY DENSITY MATRIX CONSTRAINT                       ==
!       ==  USE CONJUGATE GRADIENT METHOD TO MINIMIZE THE ABSOLUTE SQUARE     ==
!       ==  OF THE DEVIATION OF THE DENSITY MATRIX FROM THE TARGET            ==
!       ========================================================================
!       == INITIATE CONJUGATE GRADIENT ITERATION ===============================
!       == DGAMMA=DMAXDEV/DGAMMA ===========================================
        CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
     &              ,KSET(IKPT)%GAMMA,MAXDEV,FORCE)
        MAXDEV0=MAXDEV
        FORCE=0.5D0*FORCE  ! CONVERT TO WIRTINGER DERIVATIVE
        FORCE=-FORCE       ! DOWNHILL DIRECTION
        DGAMMA=FORCE       ! FIRST SEARCH DIRECTION IS DOWNHILL
!
        COSALPHA=0.D0      ! COSINUS OF ANGLE BETWEEN SEARCH DIRECTIONS
        BETA=0.D0          ! ADMIXTURE OF THE PREVIOUS SEARCH DIRECTION
        DO ITER1=1,NITER1
          CONVG=.FALSE.
          FORCEM=FORCE     ! STORE PREVIOUS FORCE
!         __ENSURE THAT DGAMMA REMAINS HERMITEAN________________________________
          CALL SPINOR$CONJUGATE(NDIMD,NCHI,DGAMMA,MAT1)
          DGAMMA=0.5D0*(DGAMMA+MAT1)
!
!         ======================================================================
!         == LINE SEARCH                                                      ==
!         ======================================================================
          CALL SPINOR$TRACEAB(NDIMD,NCHI,DGAMMA,DGAMMA,CSVAR)
          DLEN=SQRT(REAL(CSVAR,KIND=8))
          IF(ABS(DLEN).LT.1.D-9)DLEN=1.D-9  !AVOID DIVIDE BY ZERO
!
          SM=0.D0
          YM=MAXDEV
          CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCEM,DGAMMA,CSVAR)
          FM=2.D0*REAL(CSVAR,KIND=8)
          IF(TPRLINE)PRINT*,'LINESEARCH',SM*DLEN,YM,FM/DLEN,0.D0
          CONVG=(ABS(FM/DLEN).LT.TOL2) 
          IF(CONVG) THEN
            PRINT*,'FORCE VANISHES: EXIT CONSTRAINT LOOP'
            EXIT
          END IF
!
          S0=SM+SCALEINI*SIGN(1.D0,FM)/DLEN
          DO ITER2=1,NITER2 ! LINE SEARCH

            KSET(IKPT)%GAMMA=KSET(IKPT)%GAMMA+DGAMMA*(S0-SM)
            CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                  ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
     &                  ,KSET(IKPT)%GAMMA,MAXDEV,FORCE)
            FORCE=0.5D0*FORCE  ! CONVERT TO WIRTINGER DERIVATIVE
            FORCE=-FORCE       ! DOWNHILL DIRECTION

            Y0=MAXDEV
            CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCE,DGAMMA,CSVAR)
            F0=2.D0*REAL(CSVAR,KIND=8)     ! DY/DS
            D2YDS2=-(F0-FM)/(S0-SM)   ! D^2Y/DS^2
            CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCE,FORCE,CSVAR)
            FLEN=SQRT(REAL(CSVAR,KIND=8))
!
!           ====================================================================
!           == TEST DERIVATIVE                                                ==
!           ====================================================================
!!$            CALL DMFT_TESTDER(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
!!$     &                  ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
!!$     &                  ,KSET(IKPT)%GAMMA)
!
!           == CHECK CONVERGENCE ===============================================
            IF(TPRLINE)PRINT*,'LINESEARCH',S0*DLEN,Y0,F0/DLEN,D2YDS2/DLEN**2
            CONVG=(ABS(F0/DLEN).LT.TOL2) 
            CONVG=CONVG.OR.(ABS(S0-SM)*DLEN.LT.1.D-8)
            CONVG=CONVG.OR.(SQRT(Y0).LT.1.D-8) 
!!$            IF(.FALSE.) THEN
!!$!            IF(CONVG.AND.SQRT(Y0).GT.1.D-8) THEN
!!$PRINT*,'CONVERGENCE WHILE Y0 TOO LARGE: Y0=',Y0
!!$EXIT
!!$              CALL DMFT_SUMDRHO2DERIVS(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
!!$     &              ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
!!$     &              ,KSET(IKPT)%GAMMA,DGAMMA,NDIFF,DY)
!!$              PRINT*,'DY ',DY
!!$              PRINT*,'F ',FLEN
!!$              DO N=1,NX
!!$                S=-0.2+0.4D0*REAL(N-1,KIND=8)/REAL(NX-1,KIND=8)
!!$                SVAR=(0.D0,0.D0)
!!$                DO J=1,NDIFF+1
!!$                  SVAR=SVAR+DY(J)*S**(J-1)
!!$                ENDDO
!!$               CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
!!$     &                    ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
!!$     &                    ,KSET(IKPT)%GAMMA+DGAMMA*S,MAXDEV,FORCE)
!!$               PRINT*,S,SVAR,MAXDEV
!!$              ENDDO
!!$PRINT*,'KSET(IKPT)%HRHO   ',KSET(IKPT)%HRHO 
!!$PRINT*,'KSET(IKPT)%HRHO/SMAT/EV   ',KSET(IKPT)%HRHO/KSET(IKPT)%SMAT*27.211D0
!!$PRINT*,'KSET(IKPT)%GAMMA  ',KSET(IKPT)%GAMMA 
!!$PRINT*,'DGAMMA            ',DGAMMA 
!!$              CALL PLOTDOS(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
!!$     &                    ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
!!$     &                    ,KSET(IKPT)%GAMMA)
!!$STOP
!!$
!!$
!!$              DO I=1,NDIFF
!!$                DY(I)=DY(I+1)*REAL(I,KIND=8)
!!$              ENDDO
!!$              PRINT*,'DY ',DY(:NDIFF)
!!$              CALL POLYNOM$ZEROS(3,DY,Z)
!!$              DO I=1,3
!!$                CSVAR=(0.D0,0.D0)
!!$                DO J=1,4
!!$                  CSVAR=CSVAR+DY(J)*Z(I)**(J-1)
!!$                ENDDO
!!$                PRINT*,'ZERO ',Z(I),CSVAR
!!$              ENDDO
!!$              STOP 'FORCED'
!!$              CONVG=.FALSE.
!!$              CYCLE
!!$            END IF
            IF(CONVG) EXIT
!
!           ====================================================================
!           == DETERMINE SIZE OF NEXT STEP =====================================
!           ====================================================================
            IF(Y0.GT.YM) THEN 
!             == VALUE OF YM IN INCREASED. =====================================
!             == UNDO PREVIOUS STEP AND RETRY WITH SMALLER STEP SIZE ===========
              KSET(IKPT)%GAMMA=KSET(IKPT)%GAMMA-DGAMMA*(S0-SM)
              S0=0.5D0*(S0+SM) !FORLATER: -(S0-SM)+(S0+SM)/2 = -S0/2-SM*3/2
              IF(TPRLINE)PRINT*,'INTERRUPT'
              CYCLE
            END IF
            IF(D2YDS2.GT.0.D0) THEN  ! NEWTON-RAPHSON FOR POSITIVE CURVARTURE
              SP=(FM*S0-F0*SM)/(FM-F0)
            ELSE ! SWITCH TO STEEPEST DESCENT, IF CURVATURE IS NEGATIVE ========
              SP=S0+SCALEINI*SIGN(1.D0,F0)/DLEN
            END IF
!           == SWITCH ==========================================================
            SM=S0
            YM=Y0
            FM=F0
            S0=SP
          ENDDO ! END OF LINE SEARCH
          IF(.NOT.CONVG) THEN
            CALL ERROR$MSG('LINE SEARCH NOT CONVERGED')
            CALL ERROR$R8VAL('CONV CITERIUM',Y0)
            CALL ERROR$R8VAL('TOLERANCE',TOL2)
            CALL ERROR$STOP('DMFT_CONSTRAINTS')
          END IF
!
!         ======================================================================
!         == TEST DERIVATIVE                                                  ==
!         ======================================================================
!!$          CALL DMFT_TESTDER(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
!!$     &                  ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
!!$     &                  ,KSET(IKPT)%GAMMA)
!
!         ======================================================================
!         == REPORT PROGRESS                                                  ==
!         ======================================================================
          IF(MOD(ITER1,1).EQ.0) THEN
!           __ ANGLE BETWEEN CURRENT AND PREVIOUS FORCE_________________________
            CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCEM,FORCE,CSVAR)
            COSGAMMA=REAL(CSVAR,KIND=8)
            CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCE,FORCE,CSVAR)
            COSGAMMA=COSGAMMA/SQRT(REAL(CSVAR,KIND=8))
            CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCEM,FORCEM,CSVAR)
            COSGAMMA=COSGAMMA/SQRT(REAL(CSVAR,KIND=8))
!
!           __ REPORT___________________________________________________________
            STRING='("DMFT_CONSTRAINTS ITER=",I6," |DRHO|=",F10.5'
            STRING=TRIM(STRING)//'," COS(PHI)=",2F10.5," LINENITER=",I5'
            STRING=TRIM(STRING)//',E10.2'
            STRING=TRIM(STRING)//')'
            WRITE(*,FMT=TRIM(STRING))ITER1,SQRT(Y0),COSALPHA,COSGAMMA,ITER2 &
       &                            ,D2YDS2/DLEN**2
            IF(TPRFILE) THEN
              WRITE(NFIL,*)ITER1,Y0,FLEN/DLEN,D2YDS2/DLEN**2 &
       &                    ,S0,DLEN,COSALPHA,COSGAMMA,IKPT
            END IF
          END IF
!
!         ======================================================================
!         == CHECK CONVERGENCE                                                ==
!         ======================================================================
          CONVG=SQRT(MAXDEV).LT.TOL1
          IF(CONVG) EXIT
!
!         ======================================================================
!         == PREPARE NEXT SEARCH DIRECTION DGAMMA                             ==
!         ======================================================================
          IF(CGID.EQ.'PR') THEN ! POLAK-RIBIERE. BETTER FOR NONLINEARITIES_____
            CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCEM,FORCEM,CSVAR)
            FMFM=REAL(CSVAR,KIND=8)
            IF(FMFM.GT.1.D-8) THEN
              CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCE-FORCEM,FORCE,CSVAR)
              BETA=REAL(CSVAR,KIND=8)/FMFM
!             ==HTTP://EN.WIKIPEDIA.ORG/WIKI/NONLINEAR_CONJUGATE_GRADIENT_METHOD
              IF(BETA.LT.0.D0) BETA=0.D0
            ELSE
              BETA=0.D0
            END IF
          ELSE IF(CGID.EQ.'FR') THEN ! FLETCHER_REEVES__________________________
            CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCEM,FORCEM,CSVAR)
            FMFM=REAL(CSVAR,KIND=8)
            IF(FMFM.GT.1.D-8) THEN
              CALL SPINOR$TRACEAB(NDIMD,NCHI,FORCE,FORCE,CSVAR)
              BETA=REAL(CSVAR,KIND=8)/FMFM
            ELSE
              BETA=0.D0
            END IF
          ELSE
            CALL ERROR$MSG('ID FOR CONJUGATE GRADIENT METHOD NOT RECOGNIZED')
            CALL ERROR$CHVAL('CGID',CGID)
            CALL ERROR$STOP('DMFT_CONSTRAINTS')
          END IF
!
!         ======================================================================
!         == ANGLE BETWEEN THE THE NEW AND OLD SEARCH DIRECTIONS              ==
!         == COSALPHA= <FORCE+DGAMMA*SVAR|DGAMMA>=<DGAMMA|DGAMMA>*BETA        ==
!         ======================================================================
          COSALPHA=1.D0/SQRT(1.D0+(FLEN/(DLEN*BETA))**2)
!
!         __ NEXT SEARCH DIRECTION______________________________________________
          DGAMMA=FORCE+DGAMMA*BETA   ! NEW SEARCH DIRECTION
        ENDDO ! END OF CONJUGATE GRADIENT LOOP
!!$        IF(.NOT.CONVG) THEN
!!$          CALL ERROR$MSG('LOOP NOT CONVERGED')
!!$          CALL ERROR$I4VAL('IKPT',IKPT)
!!$          CALL ERROR$I4VAL('ITER1',ITER1)
!!$          CALL ERROR$R8VAL('MAX. DEVIATION',MAXDEV)
!!$          CALL ERROR$R8VAL('TOLERANCE',TOL1)
!!$          CALL ERROR$STOP('DMFT_CONSTRAINTS')
!!$        END IF
!
!       ========================================================================
!       ==  PRINT THE DEVIATION OF RHO FROM TARGET
!       ========================================================================
        IF(TPRDRHO) THEN
          PRINT*,'VIOLATION FOR DENSITY MATRIX CONSTRAINT FOR K-POINT ',IKPT
PRINT*,'GAMMA',KSET(IKPT)%GAMMA
          CALL DMFT_PRINTDRHO(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
     &              ,KSET(IKPT)%GAMMA)
        END IF
      ENDDO   !END OF LOOP OVER K-POINTS
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SPINORANGLE(NCHI,NDIMD,A,B,PHI)
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NCHI
      INTEGER(4),INTENT(IN)    :: NDIMD
      COMPLEX(8),INTENT(IN)    :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN)    :: B(NCHI,NCHI,NDIMD)
      REAL(8)                  :: PHI
      COMPLEX(8)               :: CSVAR
!     **************************************************************************
      CALL SPINOR$TRACEAB(NDIMD,NCHI,A,B,CSVAR)
      PHI=REAL(CSVAR)
      CALL SPINOR$TRACEAB(NDIMD,NCHI,A,A,CSVAR)
      PHI=PHI/SQRT(REAL(CSVAR))
      CALL SPINOR$TRACEAB(NDIMD,NCHI,B,B,CSVAR)
      PHI=PHI/SQRT(REAL(CSVAR))
      PHI=ACOS(PHI)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_LBFGS(NCHI,NDIMD,NMEMX,NMEM,AMIX,YK,SK,X0,XM,F0,FM,XP)
!     **************************************************************************
!     ** LIMITED-MEMORY BROYDEN-FLETSCHER-GOLDFARB-SHANNO  (L-BFGS) ALGORITHM **
!     ** FOR OPTIMIZING NEARLY QUADRATIC FUNCTIONS
!     ** 
!     ** SEE P 779 OF UPDATING QUASI-NEWTON MATRICES IWTH LIMITED STORAGE,    **
!     ** JORGE NOCEDAL, MATHEMATICS OF COMPUTATION, 35. P773 (1980)           **
!     **                                                                      **
!     ** INITIALIZE NMEM=0 BEFORE THE FIRST CALL OR TO RESTART HISTORY.       **
!     ** DURING OPTIMIZATION LEAVE NMEM,YK,SK UNTOUCHED. THEY ARE UPDATED     **
!     ** INSIDE THIS ROUTINE.                                                 **
!     ** FOLLOWING THE ROUTINE,                                               **
!     **  (1) SHIFT POSITIONS AND FORCES, I.E. XM=X0,FM=X0,X0=XP              **
!     **  (2) PERFORM A LINE SEARCH ALONG X(S)=X0+(X0-XM)*S                   **
!     **      UNTIL F(S)*(X0-XM)=0                                            **
!     **                                                                      **
!     **************************************PETER BLOECHL GOSLAR 2015***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NCHI
      INTEGER(4),INTENT(IN)    :: NDIMD
      INTEGER(4),INTENT(IN)    :: NMEMX ! MAX #(STEPS STORED IN HISTORY)
      INTEGER(4),INTENT(INOUT) :: NMEM  ! ACTUAL #(STEPS STORED IN HISTORY)
      COMPLEX(8),INTENT(INOUT) :: YK(NCHI,NCHI,NDIMD,NMEMX)
      COMPLEX(8),INTENT(INOUT) :: SK(NCHI,NCHI,NDIMD,NMEMX)
      REAL(8)   ,INTENT(IN)    :: AMIX  ! MIXING FACTOR
      COMPLEX(8),INTENT(IN)    :: X0(NCHI,NCHI,NDIMD) ! ACTUAL COORDINATES
      COMPLEX(8),INTENT(IN)    :: XM(NCHI,NCHI,NDIMD) ! PREVIOUS COORDINATES
      COMPLEX(8),INTENT(IN)    :: F0(NCHI,NCHI,NDIMD) ! ACTUAL FORCE
      COMPLEX(8),INTENT(IN)    :: FM(NCHI,NCHI,NDIMD) ! PREVIOUS FORCE
      COMPLEX(8),INTENT(OUT)   :: XP(NCHI,NCHI,NDIMD) ! NEXT COORDINATES
      COMPLEX(8)               :: ALPHA(NMEMX)
      COMPLEX(8)               :: BETA(NMEMX)
      COMPLEX(8)               :: RHO(NMEMX)
      COMPLEX(8)               :: Q(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: CSVAR
      INTEGER(4)               :: I
!     **************************************************************************
!
!     ==========================================================================
!     ==  UPDATE MEMORY                                                       ==
!     ==========================================================================
      Q=X0-XM
      CALL SPINOR$TRACEAB(NDIMD,NCHI,Q,Q,CSVAR)
      IF(REAL(CSVAR).GT.1.D-10) THEN
        NMEM=MIN(NMEMX,NMEM+1)
        DO I=NMEM-1,1,-1
          SK(:,:,:,I+1)=SK(:,:,:,I)
          YK(:,:,:,I+1)=YK(:,:,:,I)
        ENDDO
        SK(:,:,:,1)=X0-XM
        YK(:,:,:,1)=F0-FM
      END IF
!
!     ==========================================================================
!     ==  PREDICT NEXT POSITION                                               ==
!     ==========================================================================
      Q=F0
      DO I=1,NMEM
        CALL SPINOR$TRACEAB(NDIMD,NCHI,SK(:,:,:,I),YK(:,:,:,I),CSVAR)
        RHO(I)=1.D0/CMPLX(CSVAR,KIND=8)
        CALL SPINOR$TRACEAB(NDIMD,NCHI,SK(:,:,:,I),Q,CSVAR)
        ALPHA(I)=RHO(I)*CSVAR
        Q=Q-YK(:,:,:,I)*ALPHA(I)
      ENDDO
      Q=-AMIX*Q
      DO I=NMEM,1,-1
        CALL SPINOR$TRACEAB(NDIMD,NCHI,YK(:,:,:,I),Q,CSVAR)
        BETA(I)=RHO(I)*CSVAR
        Q=Q+SK(:,:,:,I)*(ALPHA(I)-BETA(I))
      ENDDO
      XP=X0-Q
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_TESTDER(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                  ,SMAT,SINV,HRHO,GAMMA0)
!     **************************************************************************
!     ** TEST DERIVATIVES OF THE SQUARE DEVIATION OF THE DENSITY MATRIX       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NLAU
      REAL(8)   ,INTENT(IN) :: KBT
      REAL(8)   ,INTENT(IN) :: MU
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(IN) :: FN(NLAU+1)
      COMPLEX(8),INTENT(IN) :: SMAT(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: SINV(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: HRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: GAMMA0(NCHI,NCHI,NDIMD)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: GAMMA(NCHI,NCHI,NDIMD)
      REAL(8)               :: MAXDEV0,MAXDEVP,MAXDEVM
      COMPLEX(8)            :: FORCE(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: DIR(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: DIRDAGGER(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: CSVAR
      REAL(8)               :: DER,DERP,DERM
      REAL(8)               :: SVARP,SVARM
      INTEGER(4)            :: IDIS,IDIMD
      REAL(8)               :: DF
      COMPLEX(8)            :: FORCE0(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: FORCEM(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: FORCEP(NCHI,NCHI,NDIMD)
!     **************************************************************************
!
!     ==========================================================================
!     == TEST DERIVATIVES WITHOUT WIRTINGER CALCULUS                          ==
!     ==========================================================================
      CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                  ,SMAT,SINV,HRHO,GAMMA0,MAXDEV0,FORCE)
      FORCE=-FORCE
!
!     == DISPLACEMENT ==========================================================
      DIR=FORCE
      DIR=DIR/SQRT(SUM(ABS(DIR)**2))
!
      WRITE(*,FMT='(82("="),T20," TEST DERIVATIVES ")')
      WRITE(*,FMT='(5A20)')' DISPLACEMENT ',' FORCE ',' FORCE ' &
     &                                     ,' VALUE ',' VALUE '
      CALL SPINOR$TRACEAB(NDIMD,NCHI,DIR,FORCE,CSVAR)
      DER=REAL(CSVAR,KIND=8)
      WRITE(*,FMT='(5F20.10)')0.D0,DER,DER,MAXDEV0,MAXDEV0
      DO IDIS=1,5
        SVARP=+1.D-5*REAL(IDIS,8)
        SVARM=-SVARP
        GAMMA=GAMMA0+SVARP*DIR 
        CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA,MAXDEVP,FORCE)
        FORCE=-FORCE
        CALL SPINOR$TRACEAB(NDIMD,NCHI,DIR,FORCE,CSVAR)
        DER=REAL(CSVAR,KIND=8)
        GAMMA=GAMMA0+SVARM*DIR 
        CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA,MAXDEVM,FORCE)
        FORCE=-FORCE
        CALL SPINOR$TRACEAB(NDIMD,NCHI,DIR,FORCE,CSVAR)
        DER=0.5D0*(DER+REAL(CSVAR,KIND=8))
        WRITE(*,FMT='(5F20.10)')SVARP,-(MAXDEVP-MAXDEVM)/(SVARP-SVARM),DER &
     &                          ,MAXDEVP,MAXDEVM
      ENDDO
!
!     ==========================================================================
!     ==  TEST WIRTINGER DERIVATIVES                                          ==
!     ==========================================================================
      CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA0,MAXDEV0,FORCE)
      FORCE=0.5D0*FORCE  ! CONVERT TO WIRTINGER DERIVATIVES
      FORCE=-FORCE       ! CONVERT GRADIENT INTO FORCE
!
!     == DISPLACEMENT ==========================================================
      DIR=FORCE
      DIR=DIR/SQRT(SUM(ABS(DIR)**2))
      CALL SPINOR$CONJUGATE(NDIMD,NCHI,DIR,DIRDAGGER)
!
!     ==
      WRITE(*,FMT='(82("="),T20," TEST WIRTINGER DERIVATIVES ")')
      WRITE(*,FMT='(5A20)')' DISPLACEMENT ',' FORCE ',' FORCE ' &
     &                                     ,' VALUE ',' VALUE '
!
!     == THE FACTOR 2 COMES FROM THE WIRTINGER DERIVATIVES =====================
!     == DE=SUM(CONJG(F)*DIR+CONJG(DIR)*F)*SCALE
      CALL SPINOR$MATMUL(NDIMD,NCHI,DIR,FORCE,MAT)
      CALL SPINOR$TRACE(NDIMD,NCHI,MAT,CSVAR)
      DER=REAL(CSVAR+CONJG(CSVAR),KIND=8)
      WRITE(*,FMT='(5F20.10)')0.D0,DER,DER,MAXDEV0,MAXDEV0

      DO IDIS=1,5
        SVARP=+1.D-5*REAL(IDIS,8)
        SVARM=-SVARP
!
        GAMMA=GAMMA0+SVARP*DIR
        CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA,MAXDEVP,FORCE)
        FORCE=0.5D0*FORCE ! CONVERT TO WIRTINGER DERIVATIVES
        FORCE=-FORCE      ! CONVERT GRADIENT INTO FORCE
        CALL SPINOR$MATMUL(NDIMD,NCHI,DIR,FORCE,MAT)
        CALL SPINOR$TRACE(NDIMD,NCHI,MAT,CSVAR)
        DERP=REAL(CSVAR+CONJG(CSVAR),KIND=8)
!
        GAMMA=GAMMA0+SVARM*DIR
        CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA,MAXDEVM,FORCE)
        FORCE=0.5D0*FORCE ! CONVERT TO WIRTINGER DERIVATIVES
        FORCE=-FORCE      ! CONVERT GRADIENT INTO FORCE
        CALL SPINOR$MATMUL(NDIMD,NCHI,DIR,FORCE,MAT)
        CALL SPINOR$TRACE(NDIMD,NCHI,MAT,CSVAR)
        DERM=REAL(CSVAR+CONJG(CSVAR),KIND=8)
        WRITE(*,FMT='(5F20.10)')SVARP,-(MAXDEVP-MAXDEVM)/(SVARP-SVARM) &
     &                      ,0.5D0*(DERP+DERM),MAXDEVP,MAXDEVM
      ENDDO
!STOP 'FORCED STOP IN DMFT_TESTDER'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PLOTDOS(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU0,OMEGA,FN &
     &                        ,SMAT,SINV,HRHO,GAMMA)
!     **************************************************************************
!     ** F=TR[DRHO*S*DRHO*S]
!     ** DF=TR[DFDGAMMA*DGAMMA]
!     **************************************************************************
      USE STRINGS_MODULE
      USE DMFT_MODULE, ONLY : NAT &
     &                       ,ATOMSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NLAU
      REAL(8)   ,INTENT(IN) :: KBT
      REAL(8)   ,INTENT(IN) :: MU0
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(IN) :: FN(NLAU+1)
      COMPLEX(8),INTENT(IN) :: SMAT(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: SINV(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: HRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: GAMMA(NCHI,NCHI,NDIMD)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL(4),PARAMETER  :: TPRDRHO=.FALSE.
      COMPLEX(8)            :: GINV(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: G(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: CSVAR
      INTEGER(4)            :: NU,IAT,I1,I2,IDIMD,ILAU,I
      REAL(8)               :: SVAR
      INTEGER(4),PARAMETER  :: NE=1000
      INTEGER(4)            :: IE
      REAL(8)               :: EMAX,EMIN
      REAL(8)               :: DOS,DOSBAR
      REAL(8)               :: DOS1,DOSBAR1
      REAL(8)               :: MU,X
      REAL(8)               :: EV,PI
      INTEGER(4),PARAMETER  :: NFIL=88
!     **************************************************************************
      CALL CONSTANTS$GET('EV',EV)
      PI=4.D0*ATAN(1.D0)
      EMAX=+1.D+3*KBT
      EMIN=-EMAX
      OPEN(NFIL,FILE=-'DOSTEST.DAT')
      DO IE=1,NE
        X=EMIN+(EMAX-EMIN)/REAL(NE-1,KIND=8)*REAL(IE-1,KIND=8)      
        MU=MU0
        DOS=0.D0
        DOSBAR=0.D0
!
!       ========================================================================
!       == CALCULATE DRHO=RHO-RHOBAR                                          ==
!       ========================================================================
        DO NU=1,NOMEGA
!         == CONSTRUCT LATTICE GREENS FUNCTION GBAR ============================
          GINV=(CI*OMEGA(NU)+MU)*SMAT-HRHO 
          CALL SPINOR$INVERT(NDIMD,NCHI,GINV,G)
!         == ADD UP DERIVATIVE =================================================
          CALL SPINOR$TRACEAB(NDIMD,NCHI,G,SMAT,CSVAR)
          DOSBAR=DOSBAR+2.D0*KBT*REAL(CSVAR,KIND=8) !2=NEGATIVE FREQ
!         == CONSTRUCT LATTICE GREENS FUNCTION G=== ============================
          GINV=GINV+GAMMA+X*SMAT
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            GINV(I1:I2,I1:I2,:)=GINV(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
          ENDDO
          CALL SPINOR$INVERT(NDIMD,NCHI,GINV,G)
!         == ADD UP DERIVATIVE =================================================
          CALL SPINOR$TRACEAB(NDIMD,NCHI,G,SMAT,CSVAR)
          DOS=DOS+2.D0*KBT*REAL(CSVAR,KIND=8) !2=NEGATIVE FREQ
        ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!       == LAURENT EXPANSION ===================================================
        CALL SPINOR$TRACEAB(NDIMD,NCHI,SINV,SMAT,CSVAR)
        DOSBAR=DOSBAR+FN(1)*REAL(CSVAR,KIND=8)
        DOS   =DOS   +FN(1)*REAL(CSVAR,KIND=8)
        IF(NLAU.GE.1) THEN
          MAT1=HRHO-MU*SMAT
          CALL SPINOR$TRACEAB(NDIMD,NCHI,MAT1,SINV,CSVAR)
          DOSBAR=DOSBAR+FN(2)*REAL(CSVAR,KIND=8)
          DOS   =DOS   +FN(2)*REAL(CSVAR,KIND=8)
          MAT1=-GAMMA-X*SMAT
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
              MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:)+ATOMSET(IAT)%SLOCLAUR(:,:,:,1)
          ENDDO
          CALL SPINOR$TRACEAB(NDIMD,NCHI,MAT1,SINV,CSVAR)
          DOS=DOS+FN(2)*REAL(CSVAR,KIND=8)
        END IF
        WRITE(NFIL,*)X,DOS,DOSBAR
      ENDDO ! END LOOP OVER ENERGIES
      CLOSE(NFIL)
PRINT*,'XX',-1/(PI*OMEGA(NOMEGA)),' BETA=',1/KBT
PRINT*,' NOMEGA=',NOMEGA, 'OMEGAX=',OMEGA(NOMEGA)
PRINT*,'NLAU',NLAU
PRINT*,'FN()',FN(:)        
PRINT*,'H/S ',HRHO,SMAT,SINV
CALL SPINOR$MATMUL(NDIMD,NCHI,HRHO,SINV,MAT1)
PRINT*,'H*SINV ',MAT1
CALL SPINOR$INVERT(NDIMD,NCHI,SMAT,MAT1)
CALL SPINOR$MATMUL(NDIMD,NCHI,HRHO,MAT1,MAT2)
PRINT*,'H/SMAT ',MAT2
        CALL SPINOR$TRACEAB(NDIMD,NCHI,SINV,SMAT,CSVAR)
PRINT*,'TR(S*SINV) ',CSVAR
PRINT*,'NCHI ',NCHI,NDIMD
STOP 'IN PLOTDOS'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PLOTDOS1(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU0,OMEGA,FN &
     &                        ,SMAT,SINV,HRHO,GAMMA)
!     **************************************************************************
!     ** F=TR[DRHO*S*DRHO*S]
!     ** DF=TR[DFDGAMMA*DGAMMA]
!     **************************************************************************
      USE STRINGS_MODULE
      USE DMFT_MODULE, ONLY : NAT &
     &                       ,ATOMSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NLAU
      REAL(8)   ,INTENT(IN) :: KBT
      REAL(8)   ,INTENT(IN) :: MU0
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(IN) :: FN(NLAU+1)
      COMPLEX(8),INTENT(IN) :: SMAT(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: SINV(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: HRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: GAMMA(NCHI,NCHI,NDIMD)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL(4),PARAMETER  :: TPRDRHO=.FALSE.
      COMPLEX(8)            :: GINV(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: G(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: CSVAR
      INTEGER(4)            :: NU,IAT,I1,I2,IDIMD,ILAU,I
      REAL(8)               :: SVAR
      INTEGER(4),PARAMETER  :: NE=1000
      INTEGER(4)            :: IE
      REAL(8)               :: EMAX,EMIN
      REAL(8)               :: DOS,DOSBAR
      REAL(8)               :: DOS1,DOSBAR1
      REAL(8)               :: MU,X
      REAL(8)               :: EV
      INTEGER(4),PARAMETER  :: NFIL=88
!     **************************************************************************
      CALL CONSTANTS$GET('EV',EV)
      EMAX=+1.D+2*KBT
      EMIN=-EMAX
      OPEN(NFIL,FILE=-'DOSTEST.DAT')
      DO IE=1,NE
        X=EMIN+(EMAX-EMIN)/REAL(NE-1,KIND=8)*REAL(IE-1,KIND=8)      
        MU=MU0+X
        DOS=0.D0
        DOSBAR=0.D0
!
!       ========================================================================
!       == CALCULATE DRHO=RHO-RHOBAR                                          ==
!       ========================================================================
        DO NU=1,NOMEGA
!         == CONSTRUCT LATTICE GREENS FUNCTION GBAR ============================
          GINV=(CI*OMEGA(NU)+MU)*SMAT-HRHO 
          CALL SPINOR$INVERT(NDIMD,NCHI,GINV,G)
!         == ADD UP DERIVATIVE =================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,G,SMAT,MAT1)
          CALL SPINOR$TRACEAB(NDIMD,NCHI,MAT1,MAT1,CSVAR)
          DOSBAR=DOSBAR-2.D0*REAL(CSVAR,KIND=8) !2=NEGATIVE FREQ;-1=FROM DERIVATIVE
!         == CONSTRUCT LATTICE GREENS FUNCTION G=== ============================
          GINV=GINV+GAMMA
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            GINV(I1:I2,I1:I2,:)=GINV(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
          ENDDO
          CALL SPINOR$INVERT(NDIMD,NCHI,GINV,G)
!         == ADD UP DERIVATIVE =================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,G,SMAT,MAT1)
          CALL SPINOR$TRACEAB(NDIMD,NCHI,MAT1,MAT1,CSVAR)
          DOS=DOS-2.D0*REAL(CSVAR,KIND=8) !2=NEGATIVE FREQ; -1=FROM DERIVATIVE
        ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!       == LAURENT EXPANSION ===================================================
        CALL SPINOR$TRACEAB(NDIMD,NCHI,SINV,SMAT,CSVAR)
        DOSBAR1=DOSBAR-FN(2)*REAL(CSVAR,KIND=8)
        DOS1   =DOS   -FN(2)*REAL(CSVAR,KIND=8)
!        WRITE(NFIL,*)MU/EV,DOS*EV,DOSBAR*EV,(DOS-DOSBAR)*EV
        WRITE(NFIL,*)X,DOS,DOSBAR,(DOS-DOSBAR),DOS1,DOSBAR1
      ENDDO ! END LOOP OVER ENERGIES
      CLOSE(NFIL)
PRINT*,'FN()',FN(:)        
PRINT*,'OMEGA(1/NOMEGA) ',OMEGA(1),OMEGA(NOMEGA)
PRINT*,'H/S ',HRHO,SMAT,SINV
CALL SPINOR$MATMUL(NDIMD,NCHI,HRHO,SINV,MAT1)
PRINT*,'H*SINV ',MAT1
CALL SPINOR$INVERT(NDIMD,NCHI,SMAT,MAT1)
CALL SPINOR$MATMUL(NDIMD,NCHI,HRHO,MAT1,MAT2)
PRINT*,'H/SMAT ',MAT2
        CALL SPINOR$TRACEAB(NDIMD,NCHI,SINV,SMAT,CSVAR)
PRINT*,'TR(S*SINV) ',CSVAR
PRINT*,'NCHI ',NCHI,NDIMD
STOP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_DRHO(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN,SMAT,SINV &
     &                    ,HRHO,GAMMA,DRHO)
!     **************************************************************************
!     ** CALCULATE THE DEVIATION OF THE ACTUAL DENSITTY MATRIX FROM THE       **
!     ** TARGET DENSITY MATRIX FOR A SPECIFIC K-POINT                         **
!     **************************************************************************
      USE DMFT_MODULE, ONLY : NAT &
     &                       ,ATOMSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NLAU
      REAL(8)   ,INTENT(IN) :: KBT
      REAL(8)   ,INTENT(IN) :: MU
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(IN) :: FN(NLAU+1)
      COMPLEX(8),INTENT(IN) :: SMAT(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: SINV(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: HRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: GAMMA(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(OUT):: DRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT2(NCHI,NCHI,NDIMD)
      INTEGER(4)            :: NU,IAT,I1,I2
!     **************************************************************************
!
      DRHO=(0.D0,0.D0)
!     ==========================================================================
!     == MATSUBARA SUM OF  DRHO=RHO-RHOBAR                                    ==
!     ==========================================================================
      DO NU=1,NOMEGA
        MAT1=(CI*OMEGA(NU)+MU)*SMAT-HRHO   !1/G0(IOMEGA)
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,MAT2)
        DRHO=DRHO-KBT*MAT2  ! SUBTRACT NONINTERACTING G
!
        MAT1=MAT1+GAMMA
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
        ENDDO
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,MAT2) !1/GFULL(IOMEGA)
        DRHO=DRHO+KBT*MAT2
      ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!
!     == INCLUDE NEGATIVE FREQUENCIES ==========================================
      CALL SPINOR$CONJUGATE(NDIMD,NCHI,DRHO,MAT1)
      DRHO=DRHO+MAT1
!
!     ==========================================================================
!     == REGULARIZE                                                           ==
!     ==========================================================================
      IF(NLAU.GE.1) THEN
        MAT1=(0.D0,0.D0)
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:)+ATOMSET(IAT)%SLOCLAUR(:,:,:,1)
        ENDDO
        MAT1=MAT1-GAMMA
        CALL SPINOR$MATMUL(NDIMD,NCHI,SINV,MAT1,MAT2)
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,SINV,MAT1)
        DRHO=DRHO+FN(2)*MAT1
        IF(NLAU.GE.2) THEN            
!         == SUBTRACT CONTRIIB. FROM BARE GREEN FUNCTION=======================
          MAT1=HRHO-MU*SMAT
          CALL SPINOR$MATMUL(NDIMD,NCHI,SINV,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,MAT2,MAT1)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,SINV,MAT2)
          DRHO=DRHO-FN(3)*MAT2
!         == ADD CONTRIB FROM FULL GREENS FUNCTION =============================
          MAT1=HRHO-MU*SMAT-GAMMA
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:) &
       &                       +ATOMSET(IAT)%SLOCLAUR(:,:,:,1)
          ENDDO
          CALL SPINOR$MATMUL(NDIMD,NCHI,SINV,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,MAT2,MAT1)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,SINV,MAT2)
          DRHO=DRHO+FN(3)*MAT2
          MAT1=(0.D0,0.D0)
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:) &
     &                         +ATOMSET(IAT)%SLOCLAUR(:,:,:,2)
          ENDDO
          CALL SPINOR$MATMUL(NDIMD,NCHI,SINV,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,SINV,MAT1)
          DRHO=DRHO+FN(3)*MAT1
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                        ,SMAT,SINV,HRHO,GAMMA,F,DFDGAMMA)
!     **************************************************************************
!     ** F=TR[DRHO*S*DRHO*S]                                                  **
!     ** DF=TR[DFDGAMMA*DGAMMA]                                               **
!     **************************************************************************
      USE DMFT_MODULE, ONLY : NAT &
     &                       ,ATOMSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NLAU
      REAL(8)   ,INTENT(IN) :: KBT
      REAL(8)   ,INTENT(IN) :: MU
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(IN) :: FN(NLAU+1)
      COMPLEX(8),INTENT(IN) :: SMAT(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: SINV(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: HRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: GAMMA(NCHI,NCHI,NDIMD)
      REAL(8)   ,INTENT(OUT):: F
      COMPLEX(8),INTENT(OUT):: DFDGAMMA(NCHI,NCHI,NDIMD)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL(4),PARAMETER  :: TPRDRHO=.FALSE.
      COMPLEX(8)            :: DRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: SDRHOS(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: G(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: CSVAR
      INTEGER(4)            :: NU,IAT,I1,I2,IDIMD,I
      REAL(8)               :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
!
!     ==========================================================================
!     == CALCULATE F=DEVIATION OF THE DENSITY MATRIX                          ==
!     ==========================================================================
      CALL DMFT_DRHO(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN,SMAT,SINV &
     &                    ,HRHO,GAMMA,DRHO)
!
      IF(TPRDRHO) THEN
PRINT*,'HRHO  ',HRHO
PRINT*,'GAMMA ',GAMMA
PRINT*,'DRHO  ',DRHO
PRINT*,'SMAT  ',SMAT
        CALL SPINOR$MATMUL(NDIMD,NCHI,DRHO,SMAT,MAT1)
PRINT*,'MAT1 ',MAT1
        CALL SPINOR$TRACEAB(NDIMD,NCHI,MAT1,MAT1,CSVAR)
PRINT*,'CSVAR ',CSVAR,REAL(CSVAR),SQRT(REAL(CSVAR))
        DO IDIMD=1,NDIMD
          WRITE(*,FMT='(82("="),T10," DMFT_SUMDRHO2 IDIMD=",I5," ")')IDIMD
          WRITE(*,FMT='(82("="),T10," REAL(DRHO) IDIMD=",I5," ")')IDIMD
          DO I=1,NCHI
            WRITE(*,FMT='(100F10.5)')REAL(DRHO(I,:,IDIMD),8)
            WRITE(*,FMT='(100F10.5)')REAL(MAT1(I,:,IDIMD),8)
          ENDDO
          WRITE(*,FMT='(82("="),T10," IMAG(DRHO) ",I5," ")')IDIMD
          DO I=1,NCHI
            WRITE(*,FMT='(100F10.5)')AIMAG(DRHO(I,:,IDIMD))
            WRITE(*,FMT='(100F10.5)')AIMAG(MAT1(I,:,IDIMD))
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == CALCULATE F=TR[DRHO*S*DRHO*S]                                        ==
!     == INCORPORATE OVERLAP MATRIX  SO THAT DRHO <- S*DRHO*S                 ==
!     ==========================================================================
      CALL SPINOR$MATMUL(NDIMD,NCHI,DRHO,SMAT,MAT1)
      CALL SPINOR$MATMUL(NDIMD,NCHI,SMAT,MAT1,SDRHOS)
!
      CALL SPINOR$TRACEAB(NDIMD,NCHI,MAT1,MAT1,CSVAR)
      F=REAL(CSVAR,KIND=8)
!
!     ==========================================================================
!     == CALCULATE DERIVATIVE DFDGAMMA OF F=TR[DRHO*S*DRHO*S]                 ==
!     == F(DRHO+DIR)-F(DRHO)=TR[DFDGAMMA*DIR]=SUM[CONJG(DFDGAMMA)*DIR]        ==
!     ==========================================================================
      DFDGAMMA=(0.D0,0.D0)
      DO NU=1,NOMEGA
!       == CONSTRUCT LATTICE GREENS FUNCTION ===================================
        MAT1=(CI*OMEGA(NU)+MU)*SMAT-HRHO+GAMMA
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
        ENDDO
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,G)
!       == ADD UP DERIVATIVE ===================================================
        CALL SPINOR$MATMUL(NDIMD,NCHI,SDRHOS,G,MAT1)
        CALL SPINOR$MATMUL(NDIMD,NCHI,G,MAT1,MAT2)
        DFDGAMMA=DFDGAMMA-2.D0*KBT*MAT2
      ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!     == ADD NEGATIVE FREQUENCIES ==============================================
      CALL SPINOR$CONJUGATE(NDIMD,NCHI,DFDGAMMA,MAT1)
      DFDGAMMA=DFDGAMMA+MAT1
!
!     == ADD LAURENT TERMS =====================================================   
      IF(NLAU.GE.1) THEN
        DFDGAMMA=DFDGAMMA-2.D0*FN(2)*DRHO
        IF(NLAU.GE.2) THEN            
          MAT1=HRHO-MU*SMAT-GAMMA
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:)+ATOMSET(IAT)%SLOCLAUR(:,:,:,1)
          ENDDO
          CALL SPINOR$MATMUL(NDIMD,NCHI,DRHO,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,SINV,MAT1)
          CALL SPINOR$CONJUGATE(NDIMD,NCHI,MAT1,MAT2)
          MAT1=MAT1+MAT2
          DFDGAMMA=DFDGAMMA-FN(3)*MAT1
        END IF
      END IF

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SUMDRHO2DERIVS(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                        ,SMAT,SINV,HRHO,GAMMA,DGAMMA,NDIFF,DY)
!     **************************************************************************
!     ** TAYOR EXPANSION COEFFICIENTS OF F=TR[DRHO*S*DRHO*S]
!     ** FOR THE SEARCH DIRECTION DGAMMA 
!     **************************************************************************
      USE DMFT_MODULE, ONLY : NAT &
     &                       ,ATOMSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NLAU
      REAL(8)   ,INTENT(IN) :: KBT
      REAL(8)   ,INTENT(IN) :: MU
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(IN) :: FN(NLAU+1)
      COMPLEX(8),INTENT(IN) :: SMAT(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: SINV(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: HRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: GAMMA(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: DGAMMA(NCHI,NCHI,NDIMD)
      INTEGER(4),INTENT(IN) :: NDIFF
      REAL(8)   ,INTENT(OUT):: DY(NDIFF+1) ! DY(I+1)=D^IY/DY^I; I=0,...,NDIFF
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: DRHO(NCHI,NCHI,NDIMD,NDIFF+1)
      COMPLEX(8)            :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT3(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: G(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: CSVAR
      INTEGER(4)            :: NU,IAT,I1,I2,I,J,N
      REAL(8)               :: SVAR
!     **************************************************************************
!
!     ==========================================================================
!     == CALCULATE DRHO=RHO-RHOBAR                                            ==
!     ==========================================================================
      DRHO=(0.D0,0.D0)
      DO NU=1,NOMEGA
!       == CONSTRUCT LATTICE GREENS FUNCTION ===================================
        MAT1=(CI*OMEGA(NU)+MU)*SMAT-HRHO 
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,G)
        DRHO(:,:,:,1)=DRHO(:,:,:,1)-KBT*G  ! SUBTRACT NONINTERACTING G
        MAT1=MAT1+GAMMA
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
        ENDDO
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,G)
        DRHO(:,:,:,1)=DRHO(:,:,:,1)+KBT*G
        CALL SPINOR$MATMUL(NDIMD,NCHI,G,DGAMMA,MAT1)
        MAT2=MAT1
        DO I=1,NDIFF
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,G,MAT3)
          DRHO(:,:,:,I+1)=DRHO(:,:,:,I+1)+KBT*MAT3
          IF(I.EQ.NDIFF) EXIT
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,MAT1,MAT3)
          MAT2=MAT3
        ENDDO
      ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!
!     == INCLUDE NEGATIVE FREQUENCIES AND MULTIPLY WITH OVERLAP MATRIX =========
      DO I=1,NDIFF+1
        CALL SPINOR$CONJUGATE(NDIMD,NCHI,DRHO(:,:,:,I),MAT1)
        MAT1=DRHO(:,:,:,I)+MAT1
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,SMAT,DRHO(:,:,:,I))
      ENDDO
!
!     ==========================================================================
!     == CALCULATE F=TR[DRHO*S*DRHO*S]                                        ==
!     ==========================================================================
      DY(:)=0.D0
      DO N=0,NDIFF
        DO J=0,N 
          I=N-J
          CALL SPINOR$TRACEAB(NDIMD,NCHI,DRHO(:,:,:,I+1),DRHO(:,:,:,J+1),CSVAR)
          DY(N+1)=DY(N+1)+REAL(CSVAR,KIND=8)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_PRINTDRHO(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                        ,SMAT,SINV,HRHO,GAMMA)
!     **************************************************************************
!     ** CALCULATES THE DEVIATION OF THE DENSITY MATRIX FROM THE              **
!     ** TARGET DENSITY MATRIX                                                **
!     **************************************************************************
      USE DMFT_MODULE, ONLY : NAT &
     &                       ,ATOMSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NLAU
      REAL(8)   ,INTENT(IN) :: KBT
      REAL(8)   ,INTENT(IN) :: MU
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(IN) :: FN(NLAU+1)
      COMPLEX(8),INTENT(IN) :: SMAT(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: SINV(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: HRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN) :: GAMMA(NCHI,NCHI,NDIMD)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: DRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: G(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: CSVAR
      INTEGER(4)            :: NU,IAT,I1,I2,IDIMD,I
!     **************************************************************************
!
!     ==========================================================================
!     == CALCULATE DRHO=RHO-RHOBAR                                            ==
!     ==========================================================================
      DRHO=(0.D0,0.D0)
      DO NU=1,NOMEGA
!       == CONSTRUCT LATTICE GREENS FUNCTION ===================================
        MAT1=(CI*OMEGA(NU)+MU)*SMAT-HRHO 
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,G)
        DRHO=DRHO-KBT*G  ! SUBTRACT NONINTERACTING G
        MAT1=MAT1+GAMMA
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
        ENDDO
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT1,G)
        DRHO=DRHO+KBT*G
      ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!
!     == INCLUDE NEGATIVE FREQUENCIES ==========================================
      CALL SPINOR$CONJUGATE(NDIMD,NCHI,DRHO,MAT1)
      DRHO=DRHO+MAT1
!
!     == ADD LAURENT TERMS ====================================================   
      IF(NLAU.GE.1) THEN
        MAT1=(0.D0,0.D0)
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:)+ATOMSET(IAT)%SLOCLAUR(:,:,:,1)
        ENDDO
        MAT1=MAT1-GAMMA
        CALL SPINOR$MATMUL(NDIMD,NCHI,SINV,MAT1,MAT2)
        CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,SINV,MAT1)
        DRHO=DRHO+FN(2)*MAT1
        IF(NLAU.GE.2) THEN            
          MAT1=HRHO-MU*SMAT
          CALL SPINOR$MATMUL(NDIMD,NCHI,SINV,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,MAT2,MAT1)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT1,SINV,MAT2)
          DRHO=DRHO-FN(3)*MAT2
!         == ADD CONTRIB FROM FULL GREENS FUNCTION =============================
          MAT1=HRHO-MU*SMAT-GAMMA
          DO IAT=1,NAT
            I1=ATOMSET(IAT)%ICHI1
            I2=ATOMSET(IAT)%ICHI2
            MAT1(I1:I2,I1:I2,:)=MAT1(I1:I2,I1:I2,:)+ATOMSET(IAT)%SLOCLAUR(:,:,:,2)
          ENDDO
          CALL SPINOR$MATMUL(NDIMD,NCHI,SINV,MAT1,MAT2)
          CALL SPINOR$MATMUL(NDIMD,NCHI,MAT2,SINV,MAT1)
          DRHO=DRHO+FN(3)*MAT1
        END IF
      END IF

      DO IDIMD=1,NDIMD
        WRITE(*,FMT='(82("="),T10," REAL(DRHO) IDIMD=",I5," ")')IDIMD
        DO I=1,NCHI
          WRITE(*,FMT='(100F10.5)')REAL(DRHO(I,:,IDIMD),8)
        ENDDO
        WRITE(*,FMT='(82("="),T10," IMAG(DRHO) ",I5," ")')IDIMD
        DO I=1,NCHI
          WRITE(*,FMT='(100F10.5)')AIMAG(DRHO(I,:,IDIMD))
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_MIX(XDEV)
!     **************************************************************************
!     **  MIXES DPHIDG INTO THE SELF ENERGY                                   **
!     **  DPHIDG IS THE OUTPUT SELF ENERGY OF THE SOLVER                      **
!     **                                                                      **
!     **  ON INPUT, DPHIDG IS THE SELF ENERGY RETURNED FROM THE SOLVER. ON    **
!     **  OUTPUT DPHIDG IS SET TO ZERO, BECAUSE IT HAS BEEN ABSORBED IN SLOC  **
!     **                                                                      **
!     **  ON INPUT, SLOC IS THE SELF ENERGY USED IN THE PREVIOUS ITERATION.   **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NAT,NDIMD,NLAU,NOMEGA,OMEGA,ATOMSET
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT):: XDEV   ! MAX DEVIATION OF THE SELF ENERGY
!!$      REAL(8)  ,PARAMETER :: MIX1=0.1D0
!!$      REAL(8)  ,PARAMETER :: MIX2=0.1D0
      REAL(8)   ,PARAMETER :: MIX1=1.D-1
      REAL(8)   ,PARAMETER :: MIX2=0.0D0
      INTEGER(4),PARAMETER :: NUH=50
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
      REAL(8)              :: SVAR
      INTEGER(4)           :: NLOC
      INTEGER(4)           :: IAT,NU,ILAU,IDIMD,I
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_MIX')
      IF(NUH.GT.NOMEGA) THEN
        CALL ERROR$MSG('NUH MUST NOT BE LARGER THAN NOMEGA')
        CALL ERROR$STOP('DMFT_MIX')
      END IF
      XDEV=0.D0
      DO IAT=1,NAT
        NLOC=ATOMSET(IAT)%NLOC
        IF(NLOC.LE.0) CYCLE
!
!       == DIFFERNCE OUT-IN, WHICH WILL BE SCALED AND ADDED TO CURRENT VALUES ==
        ATOMSET(IAT)%DPHIDG    =ATOMSET(IAT)%DPHIDG    -ATOMSET(IAT)%SLOC
        ATOMSET(IAT)%DPHIDGLAUR=ATOMSET(IAT)%DPHIDGLAUR-ATOMSET(IAT)%SLOCLAUR
        XDEV=MAX(XDEV,MAXVAL(ABS(ATOMSET(IAT)%DPHIDG)))
        XDEV=MAX(XDEV,MAXVAL(ABS(ATOMSET(IAT)%DPHIDGLAUR)))
!
!       == MIX FREQUENCY DEPENDENT SELF ENERGY =================================
        DO NU=1,NOMEGA
!          SVAR=1.D0/(1.D0-MIX1/(MIX1-1.D0)+MIX2*OMEGA(NU)**2)
          SVAR=1.D0/(1.D0+(1.D0-MIX1)/(MIX1+(OMEGA(NU)/OMEGA(NUH))**2))
SVAR=MIX1
          ATOMSET(IAT)%SLOC(:,:,:,NU)=ATOMSET(IAT)%SLOC(:,:,:,NU) &
     &                               +SVAR*ATOMSET(IAT)%DPHIDG(:,:,:,NU)
        ENDDO
!
!       == MIX IN WITH FREQUENCY INDEPENDENT MASS ==============================
        SVAR=1.D0
 SVAR=MIX1
        ATOMSET(IAT)%SLOCLAUR=ATOMSET(IAT)%SLOCLAUR+SVAR*ATOMSET(IAT)%DPHIDGLAUR
!
        ATOMSET(IAT)%DPHIDG=(0.D0,0.D0)
        ATOMSET(IAT)%DPHIDGLAUR=(0.D0,0.D0)
      ENDDO
!
!     ==========================================================================
!     == PRINT FOR TESTING                                                    ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,FMT='(82("="),T10," INFO FROM DMFT_MIX  ")')
        DO IAT=1,NAT
          IF(SUM(ABS(ATOMSET(IAT)%SLOC))+SUM(ABS(ATOMSET(IAT)%SLOCLAUR)).LT.1.D-8) CYCLE
          WRITE(*,FMT='(82("="),T10," ATOM",I4," ")')IAT
          NLOC=ATOMSET(IAT)%NLOC
          DO IDIMD=1,NDIMD
            WRITE(*,FMT='(82("="),T10," REAL(SIGMA)  ")')
            DO I=1,NLOC
              WRITE(*,FMT='(100F10.5)')REAL(ATOMSET(IAT)%SLOC(I,:,IDIMD,1),8)
            ENDDO
            WRITE(*,FMT='(82("="),T10," IMAG(SIGMA)  ")')
            DO I=1,NLOC
              WRITE(*,FMT='(100F10.5)')AIMAG(ATOMSET(IAT)%SLOC(I,:,IDIMD,1))
            ENDDO
            DO ILAU=1,NLAU
              WRITE(*,FMT='(82("="),T10," REAL(SLAUR(",I5,"))  ")')ILAU
              DO I=1,NLOC
                WRITE(*,FMT='(100F10.5)')REAL(ATOMSET(IAT)%SLOCLAUR(I,:,IDIMD,ILAU),8)
              ENDDO
              WRITE(*,FMT='(82("="),T10," IMAG(SLAUR(",I5,"))  ")')ILAU
              DO I=1,NLOC
                WRITE(*,FMT='(100F10.5)')AIMAG(ATOMSET(IAT)%SLOCLAUR(I,:,IDIMD,ILAU))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
!     
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
      LOGICAL(4),PARAMETER   :: TPRINT=.TRUE.
      REAL(8)                :: PHILW
      REAL(8)                :: LHFWEIGHT
      REAL(8)                :: EDC
      INTEGER(4)             :: NLOC  !#(LOCAL ORBITALS ON THIS SITE)
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
        WRITE(*,FMT='(60("."),T1'&
     &        //',"--DMFT--HF ENERGY FOR ATOM  "' &
     &        //',I3,T60,":",F20.10)')IAT,PHILW
!
!       == PRINT IF DESIRED ====================================================
        IF(TPRINT) THEN
          WRITE(*,FMT='(82("="),T10," ",A," IAT=",I4,"  ")') &
       &                        'LOCAL DENSITY MATRIX',IAT
          DO IDIMD=1,NDIMD
            DO I=1,NLOC
              WRITE(*,FMT='(4I5,100F10.5)')IAT,NU,IDIMD,I &
      &                                   ,ATOMSET(IAT)%DENMAT%RHO(I,:,IDIMD)
            ENDDO
          ENDDO
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
!       == COLLECT LOCAL HF WEIGHT (NOT REQUIRED FOR THE OTHER TERMS, BECAUSE ==
!       == THE U-TENSOR HAS ALREADY BEEN SCALED) ===============================
        LHFWEIGHT=ATOMSET(IAT)%LHFWEIGHT
        ALLOCATE(HAM(NLOC,NLOC,NDIMD))
        CALL DMFT_DC(IAT,NLOC,NDIMD,ATOMSET(IAT)%DENMAT%RHO,LHFWEIGHT,EDC,HAM)
!       __SCALING WITH LHFWEIGHT IS DONE INSIDE DMFT_DC_________________________
!FUDGE FACTOR FOR HUBBARD MODEL
PRINT*,'CAUTION! FUDGE FACTOR 1.1 FOR DOUBLE COUNTING ENERGY INCLUDED'
! FUDGE FACTOR DEPENDS IN LHFWEIGHT
EDC=EDC*1.1D0
        ETOT=ETOT+EDC
!       == FACTOR TWO IS FROM THE SPINOR REPRESENTATION ========================
!       == TR[H*DRHO] IN UP-DOWN TRANSLATES INTO                              ==
!       == 1/2 SUM_IDIMD TR[H(IDIMD)*DRHO(IDIMD)] IN THE SPINOR REPRESENTATION==
        HAM=2.D0*HAM
        ATOMSET(IAT)%DENMAT%H=ATOMSET(IAT)%DENMAT%H+HAM
        WRITE(*,FMT='(60("."),T1'&
     &        //',"--DMFT--DOUBLE COUNTING FOR ATOM  "' &
     &        //',I3,T60,":",F20.10)')IAT,LHFWEIGHT*EDC
!
!       == PRINT IF DESIRED ====================================================
        IF(TPRINT) THEN
          WRITE(*,FMT='(82("="),T10," ",A," IAT=",I4,"  ")') &
       &                        'DBLE-CNTNG HAMILTONIAN',IAT
          DO IDIMD=1,NDIMD
            DO I=1,NLOC
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
!     == DE = TR[H*DRHO] = 1/2 * \SUM_IDIM TR[ H(IDIM)*DRHO(IDIM) ]
      HAM=2.D0*HAM
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
      SUBROUTINE DMFT_DC(IAT,LMNX,NDIMD,RHO,LHFWEIGHT,EDC,HAM)
!     **************************************************************************
!     ** CALCULATES THE DFT EXCHANGE ENERGY FOR THE CORRELATED ORBITALS.      **
!     ** EDC IS TO BE ADDED! TO THE TOTAL ENERGY (MINUS SIGN IS INCLUDED)     **
!     **                                                                      **
!     ** REMARK: SOME ROUTINES ONLY WORK IN THE NON-COLLINEAR MODE            **
!     **                                                                      **
!     ** REMARK: THE ENERGY SCALING WITH LHFWEIGHT IS TAKEN CARE OF INSIDE    **
!     **         THIS ROUTINE. IT IS ALSO PASSED ON TO LMTOAUGMENTATION.      **
!     **                                                                      **
!     ** REMARK: IN CONTRAST TO THE HYBRID FUNCTIONALS, ALSO THE CORRELATION  **
!     **         CONTRIBUTION IS ACCOUNTED FOR.                               **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR=>POTPAR1,SBAR=>SBAR_NEW
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      INTEGER(4),INTENT(IN)  :: LMNX          !  #(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4),INTENT(IN)  :: NDIMD
      REAL(8)   ,INTENT(IN)  :: LHFWEIGHT  !SCREENING FACTOR FOR U-TENSOR
      COMPLEX(8),INTENT(IN)  :: RHO(LMNX,LMNX,NDIMD) ! DENSITY MATRIX
      REAL(8)   ,INTENT(OUT) :: EDC           ! -E_(DFT-EXCHANGE)
      COMPLEX(8),INTENT(OUT) :: HAM(LMNX,LMNX,NDIMD) ! HAMILTONIAN CONTRIB.
      REAL(8)   ,ALLOCATABLE :: D(:,:,:)     !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: DT(:,:,:)    !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: HT(:,:,:)    !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: H(:,:,:)  !(LMNX,LMNX,NDIMD)
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
      ALLOCATE(D(LMNX,LMNX,NDIMD))
      ALLOCATE(H(LMNX,LMNX,NDIMD))
      ALLOCATE(DT(LMNXT,LMNXT,NDIMD))
      ALLOCATE(HT(LMNXT,LMNXT,NDIMD))
      D=REAL(RHO)
      CALL LMTO_EXPANDLOCALWITHCTE('FWRD',IAT,NDIMD,LMNX,LMNXT,D,DT)
      CALL LMTOAUGMENTATION$SETI4('IAT',IAT)
      CALL LMTO_SIMPLEDC_NEW_NEW(GID,NR,NDIMD,LMNXT,LNXT,LOXT &
     &                    ,POTPAR(ISP)%TAILED%AEF &
     &                    ,LRX,AECORE,DT,LHFWEIGHT,EX,HT)
      CALL LMTOAUGMENTATION$SETI4('IAT',0) !UNSET ATOM INDEX
      CALL LMTO_EXPANDLOCALWITHCTE('BACK',IAT,NDIMD,LMNX,LMNXT,H,HT)
      EDC=-EX*LHFWEIGHT
      HAM=-CMPLX(H,KIND=8)*LHFWEIGHT
      DEALLOCATE(D)
      DEALLOCATE(DT)
      DEALLOCATE(HT)
      DEALLOCATE(H)
!PRINT*,'DOUBLE COUNTING CORRECTION ENERGY FOR ATOM=',IAT,-EX
!
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
      COMPLEX(8)             :: CVEC(NCHI)
      COMPLEX(8),ALLOCATABLE :: BMAT1(:,:)
      COMPLEX(8)             :: GRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8)             :: GFULL(NCHI,NCHI,NDIMD)
      COMPLEX(8)             :: CSVAR
      REAL(8)                :: XSUM
      REAL(8)                :: ETOT1
      REAL(8)                :: TERM1,TERM2,TERM3
      REAL(8)                :: TERM1S,TERM2S,TERM3S
      INTEGER(4)             :: NU,I,IAT,I1,I2,IDIMD,IKPT
!     **************************************************************************
      IF(NDIMD.EQ.4) THEN
        ALLOCATE(BMAT1(2*NCHI,2*NCHI))
      END IF
      TERM1S=0.D0
      TERM2S=0.D0
      TERM3S=0.D0
      ETOT=0.D0
      DO IKPT=1,NKPTL
        TERM1=0.D0
        TERM2=0.D0
        TERM3=0.D0
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
          CALL SPINOR$TRACEAB(NDIMD,NCHI,GFULL-GRHO,KSET(IKPT)%GAMMA,CSVAR)
          ETOT1=ETOT1+2.D0*REAL(CSVAR)  !FACTOR 2 FROM -OMEGA
          TERM1=TERM1+2.D0*REAL(CSVAR)  !FACTOR 2 FROM -OMEGA
!     
!         ======================================================================
!         == (SIGMA-GAMMA)*G (=(HPRIME+SIGMA-HRHO)*G)                         ==
!         ======================================================================
!         == MAT2=SIGMA-GAMMA
          CALL SPINOR$TRACEAB(NDIMD,NCHI,MAT2,GFULL,CSVAR)
          ETOT1=ETOT1+2.D0*REAL(CSVAR)  !FACTOR 2 FROM -OMEGA
          TERM2=TERM2+2.D0*REAL(CSVAR)  !FACTOR 2 FROM -OMEGA
!     
!         ======================================================================
!         == LOGARITHM TERM                                                   ==
!         ==                                                                  ==
!         == NOT USED IS THE POTENTIALLY MORE EFFICIENT ROUTE USING           ==
!         == TR[LN(A)]=LN(DET[A]). ACCORDING TO ROBERT SCHADE THE DETERMINANT ==
!         == IS CALCULATED EFFICIENTLY VIA THE LU DECOMPOSITION IN LAPACK     ==
!         ======================================================================
          CALL SPINOR$MATMUL(NDIMD,NCHI,GRHO,MAT2,MAT1) !MAT1=GRHO*MAT1
          CALL SPINOR$UNITY(NDIMD,NCHI,MAT2)            !MAT2=1
          MAT1=MAT2-MAT1    ! 1-GRHO(SIGMA+GAMMA)
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,MAT1)
          IF(NDIMD.EQ.4) THEN
            BMAT1(:NCHI,:NCHI)    =MAT1(:,:,1)
            BMAT1(NCHI+1:,:NCHI)  =MAT1(:,:,2)
            BMAT1(:NCHI,NCHI+1:)  =MAT1(:,:,3)
            BMAT1(NCHI+1:,NCHI+1:)=MAT1(:,:,4)
            CALL LIB$EIGVALNONHERMITEANC8(2*NCHI,BMAT1,CVEC)
            XSUM=0.D0
            DO I=1,2*NCHI
              XSUM=XSUM+2.D0*LOG(ABS(CVEC(I))) ! FACTOR 2 FROM -NU
            ENDDO
            ETOT1=ETOT1+XSUM 
            TERM3=TERM3+XSUM 
          ELSE
            DO IDIMD=1,NDIMD
              CALL LIB$EIGVALNONHERMITEANC8(NCHI,MAT1(:,:,IDIMD),CVEC)
              XSUM=0.D0
              DO I=1,NCHI
                XSUM=XSUM+2.D0*LOG(ABS(CVEC(I))) ! FACTOR TWO FROM -NU
              ENDDO
              IF(NDIMD.EQ.1) XSUM=XSUM*2.D0 !SPIN DEGENERACY
              ETOT1=ETOT1+XSUM
              TERM3=TERM3+XSUM
            ENDDO
          END IF
        ENDDO ! END OF LOOP OVER MATSUBARA FREQUENCIES
        ETOT=ETOT+KSET(IKPT)%WKPT*ETOT1
        TERM1S=TERM1S+KSET(IKPT)%WKPT*TERM1
        TERM2S=TERM2S+KSET(IKPT)%WKPT*TERM2
        TERM3S=TERM3S+KSET(IKPT)%WKPT*TERM3
      ENDDO ! END OF LOOK OVER K-POINTS
      ETOT=-KBT*ETOT
      TERM1S=-KBT*TERM1S
      TERM2S=-KBT*TERM2S
      TERM3S=-KBT*TERM3S
      WRITE(*,FMT='(60("."),T1'&
     &        //',"--DMFT--ENERGY FROM MATSUBARA SUM"' &
     &        //',T60,":",F20.10)')ETOT
      PRINT*,'TERM1S=',TERM1S,' TERM2S=',TERM2S,' TERM3S=',TERM3S
      PRINT*,'+TERM1S+TERM2S+TERM3S=',TERM1S+TERM2S+TERM3S
      PRINT*,'-TERM1S-TERM2S+TERM3S=',-TERM1S-TERM2S+TERM3S
      RETURN
      END




!
!     ...1.........2.........3.........4.........5.........6.........7.........8
MODULE POLES_MODULE
TYPE POLE_TYPE 
  REAL(8)                :: F
  COMPLEX(8),ALLOCATABLE :: U(:,:)
END TYPE POLE_TYPE
END MODULE POLES_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_GOFPOLE(NDIMD,NCHI,NPOLE,POLES,NOMEGA,G,NLAU,GLAU)
!     **************************************************************************
!     ** CONSTRUCT MATSUBARA GREENS FUNCTION FROM A POLE REPRESENTATION       **
!     **************************************************************************
      USE POLES_MODULE , ONLY: POLE_TYPE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: NDIMD
      INTEGER(4)     ,INTENT(IN) :: NCHI
      INTEGER(4)     ,INTENT(IN) :: NPOLE
      INTEGER(4)     ,INTENT(IN) :: NLAU
      TYPE(POLE_TYPE),INTENT(IN) :: POLES(NPOLE)
      INTEGER(4)     ,INTENT(IN) :: NOMEGA
      COMPLEX(8)     ,INTENT(OUT):: G(NCHI,NCHI,NDIMD,NOMEGA)
      COMPLEX(8)     ,INTENT(OUT):: GLAU(NCHI,NCHI,NDIMD,NLAU+1)
      COMPLEX(8)     ,PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)                    :: PI
      INTEGER(4)                 :: NDIM
      INTEGER(4)                 :: NSPIN
      INTEGER(4)                 :: IP,ISPIN,IDIMD,IDIM1,IDIM2
      INTEGER(4)                 :: I,J,NU,ILAU
      REAL(8)                    :: OMEGA(NOMEGA)
      REAL(8)                    :: FPOLE,EPOLE
      COMPLEX(8)                 :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)                 :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)                 :: MAT3(NCHI,NCHI,NDIMD)
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      NDIM=2
      NSPIN=1
      IF(NDIMD.LE.2) NDIM=1
      IF(NDIMD.EQ.2) NSPIN=2
!
!     ==========================================================================
!     == CALCULATE HBAR*OMEGA*BETA
!     ==========================================================================
      DO NU=1,NOMEGA
        OMEGA(NU)=REAL(2*NU-1,KIND=8)*PI
      ENDDO
!
!     ==========================================================================
!     == ADD UP POLES TO GREENS FUNCTION
!     ==========================================================================
      G(:,:,:,:)=(0.D0,0.D0)
      DO IP=1,NPOLE
        FPOLE=POLES(IP)%F
        EPOLE=LOG((1.D0-FPOLE)/FPOLE)
        IF(FPOLE.EQ.0.D0.OR.FPOLE.EQ.1.D0) CYCLE
        IDIMD=0
        DO ISPIN=1,NSPIN
          DO IDIM1=1,NDIM
            DO IDIM2=1,NDIM
              IDIMD=IDIMD+1
              DO I=1,NCHI
                DO J=1,NCHI
                  MAT1(I,J,IDIMD)=POLES(IP)%U(I,IDIM1)*CONJG(POLES(IP)%U(J,IDIM2))
                ENDDO !JCHI
              ENDDO  !ICHI
            ENDDO  !IDIM2
          ENDDO  ! IDIM1
        ENDDO ! ISPIN
        DO NU=1,NOMEGA
          G(:,:,:,NU)=G(:,:,:,NU)+MAT1(:,:,:)/(CI*OMEGA(NU)-EPOLE)
        ENDDO
        CALL SPINOR$UNITY(NCHI,NDIMD,MAT2)
        GLAU(:,:,:,1)=G(:,:,:,1)+MAT2(:,:,:)
        DO ILAU=2,NLAU+1
          CALL SPINOR$MATMUL(NCHI,NDIMD,MAT2,MAT1,MAT3)
          MAT2=MAT3
          GLAU(:,:,:,ILAU)=G(:,:,:,ILAU)+MAT2(:,:,:)/REAL(ILAU,KIND=8)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_POLEFORCE(NDIMD,NCHI,NPOLE,POLES &
     &                         ,NOMEGA,SIGMA,NLAU,SIGMALAU,POLEFORCES)
!     **************************************************************************
!     ** CONSTRUCT FORCES ON THE POLE REPRESENTATION                          **
!     **************************************************************************
      USE POLES_MODULE , ONLY: POLE_TYPE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: NDIMD
      INTEGER(4)     ,INTENT(IN) :: NCHI
      INTEGER(4)     ,INTENT(IN) :: NPOLE
      INTEGER(4)     ,INTENT(IN) :: NLAU
      TYPE(POLE_TYPE),INTENT(IN) :: POLES(NPOLE)
      INTEGER(4)     ,INTENT(IN) :: NOMEGA
      COMPLEX(8)     ,INTENT(IN) :: SIGMA(NCHI,NCHI,NDIMD,NOMEGA)
      COMPLEX(8)     ,INTENT(IN) :: SIGMALAU(NCHI,NCHI,NDIMD,NLAU)
      TYPE(POLE_TYPE),INTENT(OUT):: POLEFORCES(NPOLE)
      COMPLEX(8)     ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)                 :: NDIM
      INTEGER(4)                 :: NSPIN
      REAL(8)                    :: PI
      INTEGER(4)                 :: IP,ISPIN,IDIMD,IDIM1,IDIM2,IDIM
      INTEGER(4)                 :: ICHI,JCHI,NU
      REAL(8)                    :: OMEGA(NOMEGA)
      REAL(8)                    :: FPOLE,EPOLE
      COMPLEX(8)                 :: MAT1(NCHI,NCHI,NDIMD)
      COMPLEX(8)                 :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)                 :: MAT3(NCHI,NCHI,NDIMD)
      COMPLEX(8)                 :: VEC(NCHI)
      COMPLEX(8)                 :: CSVAR
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      NDIM=2
      NSPIN=1
      IF(NDIMD.LE.2) NDIM=1
      IF(NDIMD.EQ.2) NSPIN=2
!
!     ==========================================================================
!     == CALCULATE HBAR*OMEGA*BETA
!     ==========================================================================
      DO NU=1,NOMEGA
        OMEGA(NU)=REAL(2*NU-1,KIND=8)*PI
      ENDDO
!
!     ==========================================================================
!     == INITIALIZE POLEFORCES TO ZERO
!     ==========================================================================
      DO IP=1,NPOLE
        POLEFORCES(IP)%F=POLES(IP)%F
        IF(.NOT.ALLOCATED(POLEFORCES(IP)%U)) THEN
           ALLOCATE(POLEFORCES(IP)%U(NCHI,NDIM))
        END IF
        POLEFORCES(IP)%U=(0.D0,0.D0)
      ENDDO
!
!     ==========================================================================
!     == ADD UP POLES TO GREENS FUNCTION
!     ==========================================================================
      DO IP=1,NPOLE
        FPOLE=POLES(IP)%F
        EPOLE=LOG((1.D0-FPOLE)/FPOLE)
        IF(FPOLE.EQ.0.D0.OR.FPOLE.EQ.1.D0) CYCLE
        MAT1=(0.D0,0.D0)
        MAT2=(0.D0,0.D0)
        DO NU=1,NOMEGA
          MAT1=MAT1+SIGMA(:,:,:,NU)/(CI*OMEGA(NU)-EPOLE)
          MAT2=MAT2+SIGMA(:,:,:,NU)/(CI*OMEGA(NU)-EPOLE)**2
        ENDDO
        CALL SPINOR$CONJUGATE(NDIMD,NCHI,MAT1,MAT3)
        MAT1=MAT1+MAT3
        CALL SPINOR$CONJUGATE(NDIMD,NCHI,MAT2,MAT3)
        MAT2=MAT2+MAT3
!
!       == DERIVATIVE WRT. U ===================================================
        CALL SPINOR$MATVECMUL(NDIMD,NCHI,NDIM,MAT1,POLES(IP)%U,POLEFORCES(IP)%U)
!
!       == DERIVATIVE WRT. F ===================================================
        CSVAR=0.D0
        DO IDIM=1,NDIM
          CALL SPINOR$MATVECMUL(NDIMD,NCHI,1,MAT2,POLES(IP)%U(:,IDIM),VEC)
          CSVAR=CSVAR+DOT_PRODUCT(POLES(IP)%U(:,IDIM),VEC)
        ENDDO
!       __REAL(Z) INHERITS KIND PARAMETER FROM COMPLEX ARGUMENT Z
        POLEFORCES(IP)%F=REAL(CSVAR)/(FPOLE*(1.D0-FPOLE))
      ENDDO
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_POLELOC(NDIMD,NCHI,NPOLE,POLES,ICHI1,ICHI2,LOCPOLES)
!     **************************************************************************
!     ** CONSTRUCT FORCES ON THE POLE REPRESENTATION                          **
!     **************************************************************************
      USE POLES_MODULE , ONLY: POLE_TYPE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: NDIMD
      INTEGER(4)     ,INTENT(IN) :: NCHI
      INTEGER(4)     ,INTENT(IN) :: NPOLE
      INTEGER(4)     ,INTENT(IN) :: ICHI1
      INTEGER(4)     ,INTENT(IN) :: ICHI2
      TYPE(POLE_TYPE),INTENT(IN) :: POLES(NPOLE)
      TYPE(POLE_TYPE),INTENT(OUT):: LOCPOLES(NPOLE)
      INTEGER(4)                 :: NDIM
      INTEGER(4)                 :: NSPIN
      INTEGER(4)                 :: IP
!     **************************************************************************
      NDIM=2
      NSPIN=1
      IF(NDIMD.LE.2) NDIM=1
      IF(NDIMD.EQ.2) NSPIN=2
      DO IP=1,NPOLE
        LOCPOLES(IP)%F=POLES(IP)%F
        IF(.NOT.ALLOCATED(LOCPOLES(IP)%U)) THEN
          ALLOCATE(LOCPOLES(IP)%U(ICHI2-ICHI1+1,NDIM))
        END IF
        LOCPOLES(IP)%U(:,:)=POLES(IP)%U(ICHI1:ICHI2,:)
      ENDDO
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
      REAL(8)               :: FN(NLAU+1)
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
STOP  'FORCED STOP IN TESTG'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TESTRHOOFG(NORB,NDIMD,NOMEGA,NLAU,KBT,G,GLAUR)
!     **************************************************************************
!     ** PRINTS THE LOCAL DENSITY MATRIX FROM THE LOCAL GREENS FUNCTION       **
!     ** FOR TESTING PURPOSES ONLY                                            **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB     !#(SPIN ORBITALS)
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: NOMEGA   !#(POSITIVE MATSUBARA FREQUENCIES
      INTEGER(4),INTENT(IN) :: NLAU   !#(LAUENT EXPANSION TERMS SELF ENERGY)
      REAL(8)   ,INTENT(IN) :: KBT      
      COMPLEX(8),INTENT(IN) :: G(NORB,NORB,NDIMD,NOMEGA)    !GREENS FUNCTION
      COMPLEX(8),INTENT(IN) :: GLAUR(NORB,NORB,NDIMD,NLAU+1)!LAURENT EXPANSION
      REAL(8)               :: PI
      REAL(8)               :: SVAR
      INTEGER(4)            :: NU,I
      REAL(8)               :: OMEGA(NOMEGA)
      COMPLEX(8)            :: RHO(NORB,NORB,NDIMD)
      REAL(8)               :: FN(NLAU+1)
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      RHO=(0.D0,0.D0)
      DO NU=1,NOMEGA
        OMEGA(NU)=REAL(2*NU-1,KIND=8)*PI*KBT
        RHO=RHO+KBT*G(:,:,:,NU)
      ENDDO
!     == ADDING -NU. ASSUMING (T,S) REPRESENTATION OF RHO
      DO I=1,NDIMD
        RHO(:,:,I)=RHO(:,:,I)+TRANSPOSE(CONJG(RHO(:,:,I)))
      ENDDO
      CALL DMFT_REGMATSUBARA(KBT,NOMEGA,OMEGA,NLAU+1,FN)
      DO I=1,NLAU+1
        RHO=RHO+GLAUR(:,:,:,I)*FN(I)
      ENDDO
!
      CALL SPINOR_PRINTMATRIX(6,'DENSITY MATRIX (TS) IN TESTRHOOFG' &
    &                        ,1,NORB,NDIMD,NORB,RHO)
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
!WHAT IS THE FACTOR 2?
!        ATOMSET(IAT)%DENMAT%H=2.D0*ATOMSET(IAT)%DENMAT%H
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
!         __WHAT IS THE FUNNY FACTOR 2???_______________________________________
!??          MAT=2.D0*MAT   
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
      SUBROUTINE SPINOR$CONJUGATE(NDIMD,NCHI,A,ADAGGER)
!     **************************************************************************
!     **  ADAGGER IS THE HERMITEAN CONJUGATE OF A                             **
!     **                                                                      **
!     **************************************PETER BLOECHL GOSLAR 2015***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NDIMD
      INTEGER(4),INTENT(IN)  :: NCHI
      COMPLEX(8),INTENT(IN)  :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(OUT) :: ADAGGER(NCHI,NCHI,NDIMD)
      INTEGER(4)             :: IDIMD
!     **************************************************************************
      DO IDIMD=1,NDIMD
        ADAGGER(:,:,IDIMD)=TRANSPOSE(CONJG(A(:,:,IDIMD)))
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINOR$TRACE(NDIMD,NCHI,A,VAL)
!     **************************************************************************
!     **  VAL=TRACE[A]
!     **                                                                      **
!     **************************************PETER BLOECHL GOSLAR 2015***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NDIMD
      COMPLEX(8),INTENT(IN)  :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(OUT) :: VAL
      INTEGER(4)             :: I
!     **************************************************************************
      VAL=0.D0
      DO I=1,NCHI
        VAL=VAL+A(I,I,1)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINOR$TRACEAB(NDIMD,NCHI,A,B,VAL)
!     **************************************************************************
!     **  VAL=TRACE[A*B]
!     **                                                                      **
!     **************************************PETER BLOECHL GOSLAR 2015***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NDIMD
      COMPLEX(8),INTENT(IN)  :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN)  :: B(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(OUT) :: VAL
      INTEGER(4)             :: IDIMD,I,J
!     **************************************************************************
      VAL=(0.D0,0.D0)
      DO IDIMD=1,NDIMD
        DO I=1,NCHI
          DO J=1,NCHI
            VAL=VAL+A(I,J,IDIMD)*B(J,I,IDIMD)
          ENDDO
        ENDDO
      ENDDO
      VAL=0.5D0*VAL
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
      SUBROUTINE SPINOR$MATVECMUL(NDIMD,NCHI,NDIM,A,B,C)
!     **************************************************************************
!     **  MULTIPLY A MATRIX IN SPINOR REPRESENTATION WITH A VECTOR            **
!     **  
!     **  FOR NDIMD=1, NDIM=1 AND 2 IS ALLOWED. 
!     **  FOR NDIMD>1, NDIM MUST BE 2.
!     **                                                                      **
!     **  FOR NDIM=1 THE VECTORS CAN BE OF EITHER SPIN DIRECTION              **
!     **                                                                      **
!     **************************************PETER BLOECHL GOSLAR 2015***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NDIMD
      INTEGER(4),INTENT(IN)  :: NDIM
      COMPLEX(8),INTENT(IN)  :: A(NCHI,NCHI,NDIMD)
      COMPLEX(8),INTENT(IN)  :: B(NCHI,NDIM)
      COMPLEX(8),INTENT(OUT) :: C(NCHI,NDIM)
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      INTEGER(4)             :: IDIM,IDIM2
      REAL(8)                :: SVAR
      COMPLEX(8)             :: CSVAR
!     **************************************************************************
      DO IDIM=1,NDIM
        C(:,IDIM)=0.5D0*MATMUL(A(:,:,1),B(:,IDIM))
      ENDDO
      IF(NDIMD.LE.1) RETURN
!
!     == SPIN POLARIZED =======================================================
      DO IDIM=1,NDIM
        SVAR=0.5D0*REAL(3-2*IDIM,KIND=8)
        C(:,IDIM)=C(:,IDIM)+SVAR*MATMUL(A(:,:,NDIMD),B(:,IDIM))
      ENDDO
      IF(NDIMD.LE.2) RETURN
!
!     == NON-COLLINEAR ========================================================
      DO IDIM=1,NDIM
        IDIM2=3-IDIM  ! OPPOSITE SPIN DIRECTION
        CSVAR=-CI*REAL(3-2*IDIM,KIND=8)
        C(:,IDIM)=C(:,IDIM)+0.5D0*MATMUL(A(:,:,2)+CSVAR*A(:,:,3),B(:,IDIM2))
      ENDDO
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
      SUBROUTINE SPINOR$UNITY(NCHI,NDIMD,A)
!     **************************************************************************
!     **  UNITY                                                               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: NCHI
      INTEGER(4)  ,INTENT(IN)  :: NDIMD
      COMPLEX(8)  ,INTENT(OUT) :: A(NCHI,NCHI,NDIMD)
      INTEGER(4)               :: I
!     **************************************************************************
      A(:,:,:)=(0.D0,0.D0)
      DO I=1,NCHI
        A(I,I,1)=(2.D0,0.D0)
      ENDDO
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
        CALL ERROR$I4VAL('I1',I1)
        CALL ERROR$I4VAL('I2',I2)
        CALL ERROR$I4VAL('NCHI',NCHI)
        CALL ERROR$CHVAL('NAME',NAME)
        CALL ERROR$STOP('SPINOR$PRINTMATRIX')
      END IF
      IF(I2.LT.1.OR.I2.GT.NCHI) THEN
        CALL ERROR$MSG('UPPER BOUND I2 OUT OF RANGE')
        CALL ERROR$I4VAL('I1',I1)
        CALL ERROR$I4VAL('I2',I2)
        CALL ERROR$I4VAL('NCHI',NCHI)
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
      DO IDIMD=1,NDIMD
        WRITE(NFIL,FMT='(82("="),T10," ",A,"(REAL) IDIMD=",I2," ")') &
     &               TRIM(NAME),IDIMD
        DO I=I1,I2
          WRITE(NFIL,FMT='(100F10.5)') &
    &                    REAL(A(I,I1:I2,IDIMD),KIND=8)
        ENDDO
        WRITE(NFIL,FMT='(82("="),T10," ",A,"(IMAG) IDIMD=",I2," ")') &
     &               TRIM(NAME),IDIMD
        DO I=I1,I2
          WRITE(NFIL,FMT='(100F10.5)') &
    &                    AIMAG(A(I,I1:I2,IDIMD))
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
      IMPLICIT NONE
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_NININTSPECTRALDENSITY()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NDIMD,NAT,MU &
     &                      ,KSET,ATOMSET
      IMPLICIT NONE
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER     :: TPRINT=.FALSE.
      LOGICAL(4),PARAMETER     :: TTEST=.FALSE.
      INTEGER(4),PARAMETER     :: NE=2000
      REAL(8)   ,PARAMETER     :: EMIN=-2.5D0/27.211D0,EMAX=2.5D0/27.211D0
      REAL(8)   ,PARAMETER     :: DELTA=1.D-3
      REAL(8)                  :: EGRID(NE)
      REAL(8)                  :: A(NCHI,NCHI,NDIMD) ! SPECTRAL FUNCTION
                                ! LOCAL SPECTRAL FUNCTION
      REAL(8)   ,ALLOCATABLE   :: ALOC(:,:,:)  !(NLOC,NLOCM,NDIMD) 
      COMPLEX(8)               :: MATX(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: G(NCHI,NCHI,NDIMD)
      REAL(8)                  :: WKPTL
      INTEGER(4)               :: IKPT,NU,IAT,I1,I2,IDIMD
      REAL(8)                  :: PI
      INTEGER(4)               :: NFIL
      INTEGER(4)               :: IE,I
      INTEGER(4)               :: NLOC
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_NINTSPECTRALDENSITY')
      PI=4.D0*ATAN(1.D0)
      DO IE=1,NE
        EGRID(IE)=EMIN+(EMAX-EMIN)/REAL(NE-1,KIND=8)*REAL(IE-1,KIND=8)
      ENDDO
      NFIL=986
      OPEN(UNIT=NFIL,FILE='SPECTR.DAT',FORM='FORMATTED')
!
!     ==========================================================================
!     == STEP THROUGH ENERGY GRID                                             ==
!     ==========================================================================
      DO IE=1,NE
!
!       ========================================================================
!       == PERFORM BRILLOUIN-ZONE SAMPLING                                      
!       ========================================================================
        A=0.D0
        DO IKPT=1,NKPTL
          WKPTL=KSET(IKPT)%WKPT
!         == CONSTRUCT LATTICE GREENS FUNCTION =================================
          MAT=(EGRID(IE)+CI*DELTA+MU)*KSET(IKPT)%SMAT
          MAT=MAT-KSET(IKPT)%HRHO
          CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
!         == ACCOUNT FOR K-POINT (-K) ==========================================
          IF(KSET(IKPT)%TADDMINUSK) THEN
!           == TRANSPOSE SMAT,HRHO,GAMMA, BUT NEITHER CI NOR SLOC ==============
            MAT=((EGRID(IE)-CI*DELTA)+MU)*KSET(IKPT)%SMAT
            MAT=MAT-KSET(IKPT)%HRHO
            DO IDIMD=1,NDIMD
              MAT(:,:,IDIMD)=TRANSPOSE(CONJG(MAT(:,:,IDIMD)))
            ENDDO
            CALL SPINOR$INVERT(NDIMD,NCHI,MAT,MAT2)
            G=0.5D0*(G+MAT2)
          END IF
!         == MAP ONTO SPECTRAL FUNCTION ========================================
          A=A-WKPTL*AIMAG(G)/PI
        ENDDO   ! END OF LOOP OVER KPOINTS

        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          NLOC=I2-I1+1
          IF(NLOC.EQ.0) CYCLE
          ALLOCATE(ALOC(NLOC,NLOC,NDIMD))
          ALOC(:,:,:)=A(I1:I2,I1:I2,:)
          WRITE(NFIL,*)EGRID(IE),IAT,(ALOC(I,I,1),I=1,NLOC)
          DEALLOCATE(ALOC)
        ENDDO
!
      ENDDO ! END OF LOOP OVER ENERGIES
      CLOSE(NFIL)
STOP 'FORCED IN SPECTRALFUNCTION'
                              CALL TRACE$POP()
      RETURN
      END
!
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*****2ND ORDER SELF ENERGY FROM STEFFEN BACKES ********************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MSIGMA$2NDORDER(NORB,NOMEGA,NLAU,KBT,GF &
     &                          ,UTENSOR,LWENERGY,SIGMA,SLAU_OUT)
!     **************************************************************************
!     ** THIS FUNCTION SHOULD BE CALLED FROM OUTSIDE                          **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NORB
      INTEGER(4),INTENT(IN)  :: NOMEGA
      INTEGER(4),INTENT(IN)  :: NLAU
      REAL(8)   ,INTENT(IN)  :: KBT
      COMPLEX(8),INTENT(IN)  :: GF(NORB,NORB,NOMEGA)
      COMPLEX(8),INTENT(IN)  :: UTENSOR(NORB,NORB,NORB,NORB)
      REAL(8)   ,INTENT(OUT) :: LWENERGY
      COMPLEX(8),INTENT(OUT) :: SIGMA(NORB,NORB,NOMEGA)
      COMPLEX(8),INTENT(OUT) :: SLAU_OUT(NORB,NORB,NLAU)
      LOGICAL(4),PARAMETER   :: TPR=.TRUE.
      REAL(8)                :: OMEGA(NOMEGA)
      REAL(8)   ,ALLOCATABLE :: GLAU(:,:,:)
      REAL(8)   ,ALLOCATABLE :: SLAU(:,:,:)
      COMPLEX(8),ALLOCATABLE :: GFFULL(:,:,:)
      INTEGER(4)             :: SUM_MAX
      INTEGER(4)             :: NMININDEX
      INTEGER(4)             :: NMAXINDEX
      INTEGER(4)             :: NO_LAUR
      INTEGER(4)             :: N,M,NU,I,J
      REAL(8)                :: BETA
      REAL(8)                :: PI
!     **************************************************************************
                               CALL TRACE$PUSH('MSIGMA$2NDORDER')
      WRITE(*,FMT='(82("="),T10," IN MSIGMA$2NDORDER ")')
      PI=4.D0*ATAN(1.D0)
!
!     ==========================================================================
!     == CALCULATE MATSUBARA FREQUENCIES                                      ==
!     ==========================================================================
      BETA=1.D0/KBT
      DO NU=1,NOMEGA
!       __WN = (2*N+1)*PI/BETA  STEFFEN STARTS FROM ZERO, PETER STARTS FROM 1
        OMEGA(NU)=REAL(2*NU-1,KIND=8)*PI*KBT
      ENDDO
!
!     ==========================================================================
!     == PASS DIMENSIONS TO MSIGMA_MODULE AND ALLOCATE ARRAYS                 ==
!     ==========================================================================
      NO_LAUR=3   ! NO_LAUR TERMS IN THE LAURENT EXPANSION
      SUM_MAX=NOMEGA+300
SUM_MAX=NOMEGA+0
      CALL MSIGMA_DETERMINEBOUNDS(NOMEGA,SUM_MAX,NMININDEX,NMAXINDEX)
!
!     ==========================================================================
!     == REPORT DATA                                                          ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,*)'CALLED SIGMA_2NDORDER WITH THE FOLLOWING PARAMETERS:'
        WRITE(*,*)'NORB=.............',NORB
        WRITE(*,*)'NOMEGA=...........',NOMEGA
        WRITE(*,*)'NO_LAUR=..........',NO_LAUR
        WRITE(*,*)'BETA=.............',BETA
        WRITE(*,*)'SUM_MAX(NO_W+300)=',SUM_MAX
        WRITE(*,*)'NMAXINDEX........=',NMAXINDEX
        WRITE(*,*)'NMININDEX........=',NMININDEX
        WRITE(*,FMT='(T20,"  U-TENSOR  ")')
        DO I=1,NORB
          WRITE(6,FMT='(200F11.5)') (REAL(UTENSOR(I,J,I,J)),J=1,NORB)
        END DO
!
!       __ WRITE GREENS FUNCTION TO FILE________________________________________
        CALL MSIGMA_WRITEGF('G',NORB,NOMEGA,OMEGA,GF)
      END IF
!
!     ==========================================================================
!     == FIT LAURENT COEFFICIENTS TO GF.                                      ==
!     == USED TO EXPAND GREENS FUNCTION TO A LARGER MATSUBARA GRID            ==
!     ==========================================================================
      ALLOCATE(GLAU(NORB,NORB,NO_LAUR+1))
      CALL MSIGMA_SETGLAU(NORB,NOMEGA,OMEGA,GF,NO_LAUR+1,GLAU)
!
!     __REPORT FILLING__________________________________________________________
      CALL MSIGMA_SETFILLING(NORB,NOMEGA,NO_LAUR,BETA,OMEGA,GF,GLAU)
!
!     ==========================================================================
!     == EXPAND GREENS FUNCTION TO A LARGER ARRAY OF MATSUBARA FREQUENCIES    ==
!     ==========================================================================
      ALLOCATE(GFFULL(NORB,NORB,NMAXINDEX-NMININDEX+1))
      CALL MSIGMA_FILLGFFULL(NORB,NOMEGA,NLAU,KBT,GF,GLAU &
     &                      ,NMAXINDEX,NMININDEX,GFFULL)
!
!     ==========================================================================
!     == CALCULATE SELF ENERGY                                                ==
!     ==========================================================================
      SIGMA = (0.D0,0.D0)
      DO N=1,NOMEGA
!        WRITE(*,'(F6.2,A6)') (N-1)*100.D0/NOMEGA,'% DONE'
        DO M=1,NORB
!         __SIGMA(M,M,N+1) = GETSIGMA(M,N,UTENSOR)______________________________
          CALL MSIGMA_GETSIGMA(NORB,SUM_MAX,NMININDEX,NMAXINDEX &
    &                          ,M,N,BETA,UTENSOR,GFFULL,SIGMA(M,M,N))
        END DO
      END DO
!
!     ==========================================================================
!     == CALCULATE ENERGY  E=G*SIGMA?                                         ==
!     ==========================================================================
!     __LWENERGY = GETLW(GF,SIGMA)___________________________________________
      CALL MSIGMA_GETLW(NORB,NOMEGA,BETA,GF,SIGMA,LWENERGY)
!
!     ==========================================================================
!     == OBTAIN LAURENT COEFFICIENTS BY FITTING TO SIGMA 
!     ==========================================================================
      ALLOCATE(SLAU(NORB,NORB,NLAU))
      CALL MSIGMA_SETGLAU(NORB,NOMEGA,OMEGA,SIGMA,NLAU,SLAU)
      IF(NLAU.GE.1) SLAU_OUT(:,:,1)=CMPLX(SLAU(:,:,1),0.D0,KIND=8)
      IF(NLAU.GE.2) SLAU_OUT(:,:,2)=CMPLX(0.D0,-SLAU(:,:,2),KIND=8)
      IF(NLAU.GE.3) SLAU_OUT(:,:,3)=CMPLX(-SLAU(:,:,3),0.D0,KIND=8)
      IF(NLAU.GE.4) SLAU_OUT(:,:,4)=CMPLX(0.D0,SLAU(:,:,4),KIND=8)
!
!     ==========================================================================
!     == WRITE SELF ENERGIES TO FILE                                          ==
!     ==========================================================================
      IF(TPR) THEN
        CALL MSIGMA_WRITEGF('S',NORB,NOMEGA,OMEGA,SIGMA) 
      END IF
      WRITE(*,FMT='(82("="),T10," MSIGMA$2NDORDER DONE ")')
                               CALL TRACE$POP()
      RETURN
      END SUBROUTINE MSIGMA$2NDORDER
!   
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MSIGMA_WRITEGF(ID,NORB,NOMEGA,OMEGA,GF)
!     **************************************************************************
!     ** WRITE GREENS FUNCTION OR SELF ENERGY TO FILE
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: NORB
      INTEGER(4)  ,INTENT(IN) :: NOMEGA
      REAL(8)     ,INTENT(IN) :: OMEGA(NOMEGA)
      COMPLEX(8)  ,INTENT(IN) :: GF(NORB,NORB,NOMEGA)
      LOGICAL(4)  ,SAVE       :: TINI=.FALSE.
      CHARACTER(32)           :: FILEID
      INTEGER(4)              :: NFIL
      INTEGER(4)              :: NFILEV
      INTEGER(4)              :: M,N
      REAL(8)                 :: EV
      REAL(8)                 :: FACTOR
!     **************************************************************************
!
!     ==========================================================================
!     == DEFINE FILES FOR FILEHANDLER                                         ==
!     ==========================================================================
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        FILEID=+'DMFT_GREENF'
        CALL FILEHANDLER$SETFILE(FILEID,.TRUE.,-'.DMFTGFIN')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'STATUS','REPLACE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'FORM','FORMATTED')
        FILEID=+'DMFT_GREENF_EV'
        CALL FILEHANDLER$SETFILE(FILEID,.TRUE.,-'.DMFTGFIN_EV')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'STATUS','REPLACE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'FORM','FORMATTED')
        FILEID=+'DMFT_SIGMA'
        CALL FILEHANDLER$SETFILE(FILEID,.TRUE.,-'.DMFTSIGMA')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'STATUS','REPLACE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'FORM','FORMATTED')
        FILEID=+'DMFT_SIGMA_EV'
        CALL FILEHANDLER$SETFILE(FILEID,.TRUE.,-'.DMFTSIGMA_EV')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'STATUS','REPLACE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION(FILEID,'FORM','FORMATTED')
      END IF
!
!     ==========================================================================
!     == OPEN FILE AND DETERMINE CONVERSION FACTOR                            ==
!     ==========================================================================
      CALL CONSTANTS$GET('EV',EV)
      IF(ID.EQ.'G') THEN
        FACTOR=EV
        CALL FILEHANDLER$UNIT('DMFT_GREENF',NFIL)
        CALL FILEHANDLER$UNIT('DMFT_GREENF_EV',NFILEV)
      ELSE IF(ID.EQ.'S') THEN
        FACTOR=1.D0/EV
        CALL FILEHANDLER$UNIT('DMFT_SIGMA',NFIL)
        CALL FILEHANDLER$UNIT('DMFT_SIGMA_EV',NFILEV)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED (MUST BE "G" OR "S")')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('MSIGMA_WRITEGF')
      END IF
!
!     ==========================================================================
!     == WRITE FILE                                                           ==
!     ==========================================================================
      DO N=1,NOMEGA
        WRITE(NFIL,'(ES24.16E3,5X)',ADVANCE='NO') OMEGA(N)
        WRITE(NFILEV,'(ES24.16E3,5X)',ADVANCE='NO')OMEGA(N)/EV
        DO M=1,NORB
          WRITE(NFIL,'(ES24.16E3,3X,ES24.16E3,5X)',ADVANCE='NO') &
     &                      REAL( GF(M,M,N),KIND=8), AIMAG( GF(M,M,N))
          WRITE(NFILEV,'(ES24.16E3,3X,ES24.16E3,5X)',ADVANCE='NO') &
     &                      REAL( GF(M,M,N),KIND=8)*FACTOR &
     &                     ,AIMAG( GF(M,M,N))*FACTOR
        ENDDO !END LOOP OVER ORBITALS
        WRITE(NFIL,*)
        WRITE(NFILEV,*)
      ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!
!     ==========================================================================
!     == CLOSE FILES                                                          ==
!     ==========================================================================
      CALL FILEHANDLER$CLOSE('DMFT_GREENF')
      CALL FILEHANDLER$CLOSE('DMFT_GREENF_EV')
      CALL FILEHANDLER$CLOSE('DMFT_SIGMA')
      CALL FILEHANDLER$CLOSE('DMFT_SIGMA_EV')
      RETURN
      END SUBROUTINE MSIGMA_WRITEGF
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MSIGMA_SETGLAU(NORB,NOMEGA,OMEGA,GF_IN,NO_LAUR,LAUARRAY)
!     **************************************************************************
!     **  FIT LAURENT COEFF UP TO 3RD ORDER
!     **  
!     **  DETERMINE THE COEFFICIENTS OF 
!     **   REG(OMEGA)=A+B/(CI*OMEGA)**2 
!     **   IMG(OMEGA)=C/(CI*OMEGA) +D/(CI*OMEGA)**3 
!     **  FROM A PAIR OF OMEGA-VALUES. THESE COEFFICIENTS ARE THAN AVERAGED   **
!     **  A SET OF PAIRS
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NORB
      INTEGER(4),INTENT(IN)   :: NOMEGA
      INTEGER(4),INTENT(IN)   :: NO_LAUR
      REAL(8)   ,INTENT(IN)   :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN)   :: GF_IN(NORB,NORB,NOMEGA)
      REAL(8)   ,INTENT(OUT)  :: LAUARRAY(NORB,NORB,NO_LAUR)
      INTEGER(4),PARAMETER    :: N_AVG=25
      INTEGER(4),PARAMETER    :: SHIFT=25
      REAL(8)                 :: C2_TMP(NORB,NORB)
      REAL(8)                 :: COEFFS(NORB,NORB,4)
      REAL(8)                 :: FN(NORB,NORB)
      REAL(8)                 :: FM(NORB,NORB)
      REAL(8)                 :: FNI(NORB,NORB)
      REAL(8)                 :: FMI(NORB,NORB)
      REAL(8)                 :: WN,WM
      INTEGER(4)              :: NUN,NUM
      INTEGER(4)              :: M, N, NL
!     **************************************************************************
      IF(NO_LAUR.GT.4) THEN
        CALL ERROR$MSG('NO_LAUR LIMITED TO VALUES BELOW 4')
        CALL ERROR$I4VAL('NO_LAUR',NO_LAUR)
        CALL ERROR$STOP('MSIGMA_SETGLAU')
      END IF
      COEFFS = 0.D0
      DO N=1,N_AVG
        NUN = NOMEGA-N+1
        NUM = NOMEGA-N-SHIFT+1
        WN  = OMEGA(NUN)
        WM  = OMEGA(NUM)
        FN  = REAL( GF_IN(:,:,NUN),KIND=8)
        FM  = REAL( GF_IN(:,:,NUM),KIND=8)
        FNI = AIMAG(GF_IN(:,:,NUN))
        FMI = AIMAG(GF_IN(:,:,NUM))
        C2_TMP = (FN-FM) * (WN*WM)**2  / ( WN**2 - WM**2 )
        COEFFS(:,:,1)=COEFFS(:,:,1) + FN + C2_TMP/WN**2
        COEFFS(:,:,2)=COEFFS(:,:,2) + (FMI*WM**3-FNI*WN**3) / (WN**2-WM**2)
        COEFFS(:,:,3)=COEFFS(:,:,3) + C2_TMP
        COEFFS(:,:,4)=COEFFS(:,:,4) + (FMI*WM-FNI*WN)*(WN*WM)**2 /(WN**2-WM**2)
      END DO
      COEFFS = COEFFS/N_AVG
!
!     ==========================================================================
!     == MAP INTO RESULT ARRAY                                                ==
!     ==========================================================================
      WRITE(6,'(A15,I3,A15)') 'CALCULATE ONLY ',NO_LAUR,' LAURENT COEFF.'
      LAUARRAY = 0.D0
      DO NL=1,NO_LAUR
        LAUARRAY(:,:,NL) = COEFFS(:,:,NL)
      END DO
!
!     ==========================================================================
!     ==  NOW PRINT
!     ==========================================================================
      DO M=1,NORB
        WRITE(*,'(A15,I2,A13,4F8.4)')'LAURENT : ORB=',M &
     &                               ,', C(0,1,2,3)=',REAL(LAUARRAY(M,M,:))
      END DO   
      RETURN
      END SUBROUTINE MSIGMA_SETGLAU
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MSIGMA_DETERMINEBOUNDS(NO_W,SUM_MAX,NMININDEX,NMAXINDEX)
!     **************************************************************************
!     ** DETERMINE THE ARRAY BOUNDS FOR DIAGRAM SUMMATION 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NO_W
      INTEGER(4),INTENT(IN) :: SUM_MAX
      INTEGER(4),INTENT(OUT):: NMININDEX
      INTEGER(4),INTENT(OUT):: NMAXINDEX
      INTEGER(4)            :: N,N1,N2
      INTEGER(4)            :: LOWER_POLE
      INTEGER(4)            :: UPPER_POLE
!     **************************************************************************
      NMININDEX = 0
      NMAXINDEX = 0
      DO N=0,NO_W-1
         DO N1=-SUM_MAX,SUM_MAX-1
            LOWER_POLE = MIN(0,N1-N)
            UPPER_POLE = MAX(0,N1-N)
            DO N2=LOWER_POLE-SUM_MAX,UPPER_POLE+SUM_MAX-1
               IF (N1<NMININDEX) NMININDEX = N1
               IF (N2<NMININDEX) NMININDEX = N2
               IF (N+N2-N1<NMININDEX) NMININDEX = N+N2-N1
               
               IF (N1>NMAXINDEX) NMAXINDEX = N1
               IF (N2>NMAXINDEX) NMAXINDEX = N2
               IF (N+N2-N1>NMAXINDEX) NMAXINDEX = N+N2-N1
            END DO
         END DO
      END DO
      WRITE(6,*) 'MIN. INDEX=',NMININDEX
      WRITE(6,*) 'MAX. INDEX=',NMAXINDEX
      RETURN
      END SUBROUTINE MSIGMA_DETERMINEBOUNDS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MSIGMA_SETFILLING(NORB,NOMEGA,NLAU,BETA,OMEGA,GF,GLAU)
!     **************************************************************************
!     ** REPORT DIAGONAL ORBITAL OCCUPATIONS
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NLAU
      REAL(8)   ,INTENT(IN) :: BETA
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN) :: GF(NORB,NORB,NOMEGA)
      REAL(8)   ,INTENT(IN) :: GLAU(NORB,NORB,NLAU+1)
      INTEGER(4)            :: ORB,NU
      REAL(8)               :: FILL
!     **************************************************************************
      DO ORB=1,NORB
        FILL = 0.D0
        DO NU=1,NOMEGA
          FILL=FILL+REAL(GF(ORB,ORB,NU),KIND=8)
        END DO
        IF(NLAU+1.GE.3) THEN
          DO NU=1,NOMEGA
            FILL=FILL+GLAU(ORB,ORB,3)/OMEGA(NU)**2
          END DO
        END IF
        FILL = FILL*2.0/BETA 
        IF(NLAU+1.GE.2) FILL=FILL + GLAU(ORB,ORB,2)*0.5D0 
        IF(NLAU+1.GE.3) FILL=FILL - GLAU(ORB,ORB,3)*BETA/4.0
        WRITE(*,FMT='("ORB=",I3," FILLING=",F9.6)')ORB,FILL
      END DO
      RETURN      
      END SUBROUTINE MSIGMA_SETFILLING
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MSIGMA_FILLGFFULL(NORB,NOMEGA,NLAU,KBT,GF,GLAU &
     &                            ,NMAXINDEX,NMININDEX,GFFULL)
!     **************************************************************************
!     ** EXPAND GREENS FUNCTION TO A LARGER MATSUBARA GRID                    **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NLAU
      REAL(8)   ,INTENT(IN) :: KBT
      COMPLEX(8),INTENT(IN) :: GF(NORB,NORB,NOMEGA)
      REAL(8)   ,INTENT(IN) :: GLAU(NORB,NORB,NLAU+1)
      INTEGER(4),INTENT(IN) :: NMININDEX
      INTEGER(4),INTENT(IN) :: NMAXINDEX
      COMPLEX(8),INTENT(OUT):: GFFULL(NORB,NORB,NMAXINDEX-NMININDEX+1)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)            :: N,ILAU,NU
      REAL(8)               :: OMEGANU
      COMPLEX(8)            :: CSVAR
      COMPLEX(8)            :: GF1(NORB,NORB)
      REAL(8)               :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      GFFULL=(0.D0, 0.D0)
      DO N=NMININDEX,NMAXINDEX
        NU=N+1
!       == THE CONTENT HAS BEEN THE FUNCTION GETGF ===========================
!       ==  IF WE HAVE DATA, RETURN IT =======================================
        IF (N<NOMEGA .AND. -NOMEGA<=N) THEN
          IF (N>=0) THEN
            GF1 = GF(:,:,NU)
          ELSE 
            GF1 = CONJG(TRANSPOSE(GF(:,:,ABS(N)) ))
          END IF 
        ELSE
!         __RETURN EXTRAPOLATION
          OMEGANU=REAL(2*NU-1,KIND=8)*PI*KBT
          GF1=(0.D0,0.D0)
          CSVAR=(1.D0,0.D0)
          DO ILAU=2,NLAU+1
            CSVAR=CSVAR/(CI*OMEGANU)   
            GF1 = GF1+CSVAR*GLAU(:,:,ILAU)
          ENDDO
        END IF
        GFFULL(:,:,N-NMININDEX+1)=GF1
      ENDDO
      RETURN
      END SUBROUTINE MSIGMA_FILLGFFULL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MSIGMA_GETSIGMA(NORB,SUM_MAX,NMININDEX,NMAXINDEX &
     &                          ,ORB,N_IN,BETA,U,GFFULL,S)
!     **************************************************************************
!     ** CALCULATE SIGMA IN 2ND ORDER                                         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NORB ! #(ORBITALS)
      INTEGER(4),INTENT(IN) :: SUM_MAX 
      INTEGER(4),INTENT(IN) :: NMAXINDEX
      INTEGER(4),INTENT(IN) :: NMININDEX
      INTEGER(4),INTENT(IN) :: ORB  ! ORBITAL INDEX
      INTEGER(4),INTENT(IN) :: N_IN ! MATSUBARA FREQUENCY INDEX (STARTS WITH 1)
      REAL(8)   ,INTENT(IN) :: BETA
      COMPLEX(8),INTENT(IN) :: U(NORB,NORB,NORB,NORB+1)
      COMPLEX(8),INTENT(IN) :: GFFULL(NORB,NORB,NMAXINDEX-NMININDEX+1)
      COMPLEX(8),INTENT(OUT):: S    ! SELF ENERGY
      INTEGER(4)            :: M, N1, LOWER_POLE, UPPER_POLE,N
      INTEGER(4)            :: IND1A,IND1B,IND2A,IND2B
      COMPLEX(8)            :: S_PRE_N1
!     **************************************************************************
      S = (0.D0,0.D0)
      N=N_IN-1
      !2ND ORDER DIAGRAM (D+E) EQ. 63
      DO M=1,NORB
        IF(M.NE.ORB) THEN 
          S_PRE_N1 = (0.D0,0.D0)
          DO N1=-SUM_MAX,SUM_MAX-1
            !DIAGRAM (D+E)
            !G(N2) HAS SYMMETRY AXIS AT N2 = 0
            !G(N+N2-N1) HAS SYMMETRY AXIS AT N2 = N1-N
            !CENTER IS AT (N1-N-1)/2
            !SUM FROM LEFT POLE-SUM_MAX UP TO RIGHT POLE+SUM_MAX
            LOWER_POLE = MIN(0, N1-N);
            UPPER_POLE = MAX(0, N1-N);
            IND1A=LOWER_POLE-SUM_MAX-NMININDEX+1
            IND1B=UPPER_POLE+SUM_MAX-NMININDEX+0
            IND2A=IND1A+N-N1
            IND2B=IND1B+N-N1
            S_PRE_N1=S_PRE_N1 + GFFULL(ORB,ORB,N1-NMININDEX+1) &
     &              * DOT_PRODUCT(CONJG(GFFULL(M,M,IND1A:IND1B))  &
     &                                 ,GFFULL(M,M,IND2A:IND2B) )
          END DO
          S = S - REAL(U(ORB,M,ORB,M),KIND=8)**2 * S_PRE_N1/BETA**2
        END IF
      END DO
      RETURN
      END SUBROUTINE MSIGMA_GETSIGMA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MSIGMA_GETLW(NORB,NOMEGA,BETA,GF,SIGMA,LW)
!     **************************************************************************
!     ** CALCULATE LUTTINGER-WARD ENERGY BY G*SIGMA/4 WHICH IS VALID FOR      **
!     **  SECOND ORDER                                                        **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NORB
      INTEGER(4),INTENT(IN)  :: NOMEGA
      REAL(8)   ,INTENT(IN)  :: BETA
      COMPLEX(8),INTENT(IN)  :: GF(NORB,NORB,NOMEGA)
      COMPLEX(8),INTENT(IN)  :: SIGMA(NORB,NORB,NOMEGA)
      REAL(8)   ,INTENT(OUT) :: LW
      INTEGER(4)             :: N,M
!     **************************************************************************
      LW = 0.D0
      DO N=1,NOMEGA
        DO M=1,NORB
!         == 2 TIMES REAL ACCOUNTS FOR NEGATIVE MATSUBARA FREQUENCIES ==========
          LW = LW +2.D0*REAL(GF(M,M,N)*SIGMA(M,M,N),KIND=8)
         END DO
      END DO
      LW=LW/(4.D0*BETA)
      WRITE(6,'(A11,ES14.6E3)') 'LW-ENERGY: ',LW
      RETURN
      END SUBROUTINE MSIGMA_GETLW
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_ANALYTICCONTINUATION(NORB,NOMEGA,OMEGA,GREEN &
     &                                    ,NE,EMIN,EMAX,SPECTRAL)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NORB
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN) :: GREEN(NORB,NORB,NOMEGA)
      INTEGER(4),INTENT(IN) :: NE
      REAL(8)   ,INTENT(IN) :: EMIN
      REAL(8)   ,INTENT(IN) :: EMAX
      COMPLEX(8),INTENT(OUT):: SPECTRAL(NORB,NORB,NE)
      INTEGER(4),PARAMETER  :: NPOLE=10
      INTEGER(4),PARAMETER  :: NP=2*NPOLE
      INTEGER(4),PARAMETER  :: NIMAGX=1
      real(8)   ,parameter  :: imagx=5.d-3
      COMPLEX(8),PARAMETER  :: CONE=(1.D0,0.D0)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)               :: e,gamma
      INTEGER(4)            :: NO
      INTEGER(4)            :: NUP(NP)   ! INDEX FOR MATSUBARA FREQS
      COMPLEX(8)            :: ZP(NP)
      COMPLEX(8)            :: FP(NP)
      COMPLEX(8),ALLOCATABLE:: ZO(:)
      COMPLEX(8),ALLOCATABLE:: FO(:)
      COMPLEX(8)            :: csvar
      INTEGER(4)            :: I,J,IE,IO1,IO2,IND
integer(4):: nfil1,nfil2
!     **************************************************************************
!     ==========================================================================
!     == SELECT POINTS TO BE FITTED                                       
!     ==========================================================================
      DO I=1,NP
        NUP(I)=I
      ENDDO
      NUP(NP)=NOMEGA
      DO I=1,3
        NUP(NP-I)=NINT(REAL(NUP(NP-I+1),KIND=8)/2.D0)
      ENDDO      

      DO I=1,NP
        IF(NUP(I).GT.0) THEN
          ZP(I)=CI*OMEGA(NUP(I))
        ELSE
          ZP(I)=-CI*OMEGA(-NUP(I))
        END IF
      ENDDO
!
!     ==========================================================================
!     == ENERGY GRID
!     ==========================================================================
      NO=NIMAGX*NE
      ALLOCATE(ZO(NO))
      ALLOCATE(FO(NO))
      IND=0
      DO I=1,NE
        e=EMIN+(EMAX-EMIN)*REAL(I-1,KIND=8)/REAL(NE-1,KIND=8)
        DO J=1,NIMAGX
          gamma=IMAGX*REAL(J,KIND=8)/REAL(NIMAGX,KIND=8)
          IND=IND+1
          ZO(IND)=cmplx(e,gamma)
        ENDDO
      ENDDO
!!$ind=0
!!$do i=1,ne
!!$  DO J=1,NIMAGX
!!$    e=omega(nomega)/real(ne)*real(i)
!!$    ind=ind+1
!!$    ZO(IND)=cmplx(0.d0,e)
!!$  enddo
!!$enddo
!
!     ==========================================================================
!     == ENERGY GRID
!     ==========================================================================
      DO IO1=1,NORB
        DO IO2=IO1,NORB
!
          DO I=1,NP
            IF(NUP(I).GT.0) THEN
              FP(I)=GREEN(IO1,IO2,NUP(I))
            ELSE
              FP(I)=CONJG(GREEN(IO2,IO1,-NUP(I)))
            END IF
          ENDDO
!
!!$print*,'before dmft_padels'
          CALL DMFT_PADELS(NPOLE,NP,ZP,FP,NO,ZO,FO)
!!$print*,'after dmft_padels'
!!$open(newunit=nfil1,file='out1.dat')
!!$do i=1,np
!!$  write(nfil1,*)aimag(zp(i)),real(fp(i)),aimag(fp(i))
!!$enddo
!!$close(nfil1)
!!$open(newunit=nfil1,file='out2.dat')
!!$do i=1,no
!!$  write(nfil1,*)aimag(zo(i)),real(fo(i)),aimag(fo(i))
!!$enddo
!!$close(nfil1)
!!$stop 'forced'
!
          IND=0
          DO I=1,NE
            DO J=1,NIMAGX
              IND=IND+1
              csvar=FO(IND)
            ENDDO
            SPECTRAL(IO1,IO2,I)=csvar
            SPECTRAL(IO2,IO1,I)=CONJG(csvar)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_PADELS(NPOLE,NP,ZP,FP,NO,ZO,FO)
!     **************************************************************************
!     **  LEAST-SQUARE PADE APPROXIMATION                                     **
!     **  METHOD TAKEN FROM SCHOETT ARXIV1511.03496 EQS. 8,9,10,11            **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NPOLE
      INTEGER(4),INTENT(IN) :: NP       ! NUMBER OF DATA POINTS TO BE FITTED
      COMPLEX(8),INTENT(IN) :: ZP(NP)   ! COMPLEX ARGUMENTS OF DATA POINTS
      COMPLEX(8),INTENT(IN) :: FP(NP)   ! FUNCTION VALUES OF DATA POINTS
      INTEGER(4),INTENT(IN) :: NO       ! NUMBER OF INTERPOLATED VALUES 
      COMPLEX(8),INTENT(IN) :: ZO(NO)   ! COMPLEX ARGUMENTS OF THE INTERPOLATION
      COMPLEX(8),INTENT(OUT):: FO(NO)   ! INTERPOLATED FUNCTION VALUES
      COMPLEX(8),PARAMETER  :: CONE=(1.D0,0.D0)
      COMPLEX(8),PARAMETER  :: Ci=(0.D0,1.D0)
      COMPLEX(8),PARAMETER  :: Cnull=(0.D0,0.D0)
      COMPLEX(8)            :: Amat(NP,2*NPOLE)
      COMPLEX(8)            :: B(NP)
      COMPLEX(8)            :: X(2*NPOLE)
      COMPLEX(8)            :: cdenom(nPOLE+1)
      COMPLEX(8)            :: apole(NPOLE)
      COMPLEX(8)            :: zpole(NPOLE)
      COMPLEX(8),allocatable:: P(:),Q(:)
      real(8)               :: re,im
      INTEGER(4)            :: J,k
      logical(4),parameter  :: tprint=.true.
!     **************************************************************************
!
!     ==========================================================================
!     == set up system of equations for least-square fit
!     ==========================================================================
!     == MATRIX K (EQ 11)
      Amat(:,1)=CONE
      DO J=2,NPOLE
        Amat(:,J)=Amat(:,J-1)*ZP(:)
      ENDDO
      Amat(:,NPOLE+1)=-FP
      DO J=NPOLE+2,2*NPOLE
        Amat(:,J)=Amat(:,J-1)*ZP(:)
      ENDDO
!     ==  VECTOR V (EQ 11)
      B=-Amat(:,2*NPOLE)*ZP(:)
!
!     ==========================================================================
!     == solve an over-determined equation system A*x=b via Sing Val. decomp. ==
!     ==========================================================================
      CALL DMFT_MYSVDEQS(NP,2*NPOLE,Amat,B,X)
!
!     ==========================================================================
!     == determine poles                                                      ==
!     ==========================================================================
      cdenom(:npole)=x(npole+1:)
      cdenom(npole+1)=cone
      call polynom$zerosc8(npole+1,cdenom,zpole)

!!$apole=(0.D0,0.D0)
!!$do j=1,npole+1
!!$  apole(:)=apole(:)+cdenom(j)*zpole(:)**(j-1)
!!$enddo
!!$print*,'cdenom ',cdenom
!!$print*,'zpole ',zpole
!!$print*,'nonzero ',maxval(abs(apole))
!!$stop

      apole=(0.D0,0.D0)
      DO J=1,NPOLE
        apole(:)=apole(:)+X(J)*Zpole(:)**(J-1)
      ENDDO
      do j=1,npole
        do k=1,npole
          if(k.eq.j) cycle
          apole(j)=apole(j)/(zpole(j)-zpole(k))
        enddo
      enddo
!
!     ==========================================================================
!     == shift poles                                                          ==
!     ==========================================================================
      do j=1,npole
        re=real(zpole(j),kind=8)
        im=aimag(zpole(j))
        im=min(1.d-8,im)
        zpole(j)=cmplx(re,im)
      enddo
      IF(TPRINT) THEN
        DO J=1,NPOLE
          WRITE(*,FMT='("POLE=",2e15.5," AMPLITUDE=",2e15.5)')ZPOLE(J),APOLE(J)
        ENDDO
      END if
!apole=real(apole)
!
!     ==========================================================================
!     == evaluate rational polynomial                                         ==
!     ==========================================================================
      fo(:)=cnull
      do j=1,npole
        fo(:)=fo(:)+apole(j)/(zo(:)-zpole(j))
      enddo

!!$      allocate(p(no))
!!$      allocate(q(no))
!!$      P=(0.D0,0.D0)
!!$      Q=(0.D0,0.D0)
!!$      DO J=1,NPOLE
!!$        P(:)=P(:)+X(J)      *ZO(:)**(J-1)
!!$        Q(:)=Q(:)+X(NPOLE+J)*ZO(:)**(J-1)
!!$      ENDDO
!!$      q(:)=q(:)+zo(:)**npole
!!$      FO(:)=P(:)/Q(:)
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE DMFT_PADELSpoles(NPOLE,NP,ZP,FP,NO,ZO,FO)
!!$!     **************************************************************************
!!$!     **  LEAST-SQUARE PADE APPROXIMATION                                     **
!!$!     **  METHOD TAKEN FROM SCHOETT ARXIV1511.03496 EQS. 8,9,10,11            **
!!$!     **************************************************************************
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: NPOLE
!!$      INTEGER(4),INTENT(IN) :: NP       ! NUMBER OF DATA POINTS TO BE FITTED
!!$      COMPLEX(8),INTENT(IN) :: ZP(NP)   ! COMPLEX ARGUMENTS OF DATA POINTS
!!$      COMPLEX(8),INTENT(IN) :: FP(NP)   ! FUNCTION VALUES OF DATA POINTS
!!$      COMPLEX(8),INTENT(IN) :: Zpole(Npole) ! poles of the Greens function
!!$      COMPLEX(8),INTENT(OUT):: FO(NO)   ! INTERPOLATED FUNCTION VALUES
!!$      COMPLEX(8),PARAMETER  :: CONE=(1.D0,0.D0)
!!$      COMPLEX(8)            :: A(NP,2*NPOLE)
!!$      COMPLEX(8)            :: B(NP)
!!$      COMPLEX(8)            :: X(2*NPOLE)
!!$      COMPLEX(8),allocatable:: P(:),Q(:)
!!$      INTEGER(4)            :: J
!!$!     **************************************************************************
!!$      allocate(p(no))
!!$      allocate(q(no))
!!$!
!!$!     ==========================================================================
!!$!     == set up system of equations for least-square fit
!!$!     ==========================================================================
!!$!     == MATRIX K (EQ 11)
!!$      A(:,1)=CONE
!!$      DO J=2,NPOLE
!!$        A(:,J)=A(:,J-1)*ZP(:)
!!$      ENDDO
!!$      A(:,NPOLE+1)=-FP
!!$      DO J=NPOLE+2,2*NPOLE
!!$        A(:,J)=A(:,J-1)*ZP(:)
!!$      ENDDO
!!$!     ==  VECTOR V (EQ 11)
!!$      B=-A(:,2*NPOLE)*ZP(:)
!!$!
!!$!     ==========================================================================
!!$!     == solve an over-determined equation system A*x=b via Sing Val. decomp. ==
!!$!     ==========================================================================
!!$      CALL DMFT_MYSVDEQS(NP,2*NPOLE,A,B,X)
!!$!
!!$!     ==========================================================================
!!$!     == evaluate rational polynomial                                         ==
!!$!     ==========================================================================
!!$      call polynomial$zerosc8(npole+1,x(:npole+1),znumerator
!!$      P=(0.D0,0.D0)
!!$      Q=(0.D0,0.D0)
!!$      DO J=1,NPOLE
!!$        P(:)=P(:)+X(J)      *ZO(:)**(J-1)
!!$        Q(:)=Q(:)+X(NPOLE+J)*ZO(:)**(J-1)
!!$      ENDDO
!!$      q(:)=q(:)+zo(:)**npole
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_MYSVDEQS(M,N,A,B,X)
!     **************************************************************************
!     ** SOLVES A LINEAR EQUATION SYSTEM A*X=B FOR X  IN THE LEAST-SQUARE     **
!     ** SENSE. THAT IS THE NUMBER M OF EQUATION CAN BE MUCH LARGER THEN      **
!     ** THE NUMBER OF UNKNOWNS X.
!     ** M IS LARGER OR EQUAL TO N
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: A(M,N)
      COMPLEX(8),INTENT(IN) :: B(M)
      COMPLEX(8),INTENT(OUT):: X(N)
      REAL(8)   ,allocatable:: S(:)    !(m)    !singular values of A
      COMPLEX(8),allocatable:: U(:,:)  !(M,m)  !(LDU,N) left vectors
      COMPLEX(8),allocatable:: VT(:,:) !(N,N) !(LDVT,N) right vectors
!     **************************************************************************
      IF(M.LT.N) THEN
        CALL ERROR$MSG('ERROR: M MUST NOT BE SMALLER THAN N')
        CALL ERROR$STOP('DMFT_MYSVDEQS')
      END IF
      IF(N.LE.1) THEN
        CALL ERROR$MSG('N MUST BE LARGER THAN 1')
        CALL ERROR$STOP('DMFT_MYSVDEQS')
      END IF
      allocate(s(m))
      allocate(u(m,m))
      allocate(vt(n,n))

!     ==========================================================================
!     == SINGULAR VALUE DECOMPOSITION A=U*S*V^DAGGER
!     ==========================================================================
!     CAUTION: ROUTINE PRODUCES THE FULL M*M MATRIX U INCLUDING NULL-VECTORS
      CALL LIB$SVDC8(M,N,A,U,S,VT)
!
!     ==========================================================================
!     == A*X=B => U*S*VT*X=B => X=V*1/S*UT*B
!     ==========================================================================
      X=MATMUL(CONJG(TRANSPOSE(U(:,:N))),B)
      X=X/S(:n)
      X=MATMUL(CONJG(TRANSPOSE(VT)),X)
!
      deallocate(s)
      deallocate(u)
      deallocate(vt)
      RETURN
      END
   
   
!!
!! EXAMPLE
!!
!!$PROGRAM MAIN
!!$   USE PARAMS
!!$   USE MSIGMA, ONLY : SIGMA_2NDORDER
!!$   
!!$  INTEGER,PARAMETER                          :: NO_ORB=10, NO_W=100, NLAU=1
!!$  REAL(8), PARAMETER                        :: KBT=0.0009500036733
!!$  COMPLEX(8), DIMENSION(NO_ORB,NO_ORB,NO_W) :: GF
!!$  REAL(8), DIMENSION(NO_ORB,NO_ORB,NLAU)    :: GLAU
!!$  COMPLEX(8), DIMENSION(NO_ORB,NO_ORB,NO_ORB,NO_ORB) :: U
!!$
!!$  REAL(8)                                   :: LWENERGY
!!$  COMPLEX(8), DIMENSION(NO_ORB,NO_ORB,NO_W) :: SIGMA
!!$  REAL(8), DIMENSION(NO_ORB,NO_ORB,NLAU)    :: SLAU
!!$  
!!$  INTEGER  :: M,M1,M2,N
!!$  REAL(8) :: READTMP(NO_ORB), READTMP2(2*NO_ORB+1)
!!$  
!!$  U = (0.0,0.0)
!!$  OPEN(UNIT=1001,FILE="UMATRIX.DAT")
!!$  DO M1=1,NO_ORB
!!$     READ(1001,*) READTMP 
!!$     DO M2=1,NO_ORB
!!$        U(M1,M2,M1,M2) = READTMP(M2)
!!$     END DO
!!$  END DO
!!$  CLOSE(1001)
!!$  OPEN(UNIT=1002,FILE="G_LOC.DAT")
!!$  DO N=1,NO_W
!!$     READ(1002,*) READTMP2 
!!$     DO M=0,NO_ORB-1
!!$        GF(M+1,M+1,N) = CMPLX(READTMP2(2*M+2),READTMP2(2*M+3),KIND=8)
!!$     END DO
!!$  END DO
!!$  CLOSE(1002)
!!$  
!!$  CALL SIGMA_2NDORDER(NO_ORB, NO_W, NLAU, KBT, GF, &
!!$                    & GLAU, U, LWENERGY, SIGMA, SLAU)
!!$   
!!$END PROGRAM MAIN
