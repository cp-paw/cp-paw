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
  COMPLEX(8),POINTER   :: GNI(:,:,:,:)     => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: GNILAUR(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,NLAU+1)
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
  REAL(8)   ,POINTER :: E(:,:)        ! (NB,NSPIN)         
  COMPLEX(8),POINTER :: PIPSI(:,:,:,:)!(NDIM,NCHI,NB,NSPIN)
  COMPLEX(8),POINTER :: RHO(:,:,:)    !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: RHOADAPTED(:,:,:)    !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: GAMMA(:,:,:)  !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: HRHO(:,:,:)   !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: SINV(:,:,:)   !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: SMAT(:,:,:)   !(NCHI,NCHI,NDIMD)
END TYPE KSET_TYPE
LOGICAL(4),PARAMETER   :: TON=.TRUE.
LOGICAL(4),SAVE        :: TINI=.FALSE.
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
      INTEGER(4),PARAMETER   :: NITERX=2000
      REAL(8)   ,PARAMETER   :: ETOL=1.D-8
      REAL(8)   ,PARAMETER   :: SIGMATOL=1.D-8
      REAL(8)                :: XDEV ! X(DEVIATION OF THE SELF ENERGY)
      REAL(8)                :: ETOT
      REAL(8)                :: ELAST
      REAL(8)                :: SVAR
      INTEGER(4)             :: ITER
      LOGICAL(4) :: TREPEAT=.FALSE. ! USED FOR GRADIENT TEST
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
      CALL DMFT_OVERLAPOFK()
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
!      CALL DMFT_CONSTRAINTS()
!
!     ==========================================================================
!     ==  CONSTRUCT LOCAL NATURAL ORBITALS                                    ==
!     ==========================================================================
      CALL DMFT_NATORB()
!      CALL DMFT_NININTSPECTALDENSITY() ! USED FOR TESTING
!      CALL DMFT_GLOCNI()
      CALL DMFT_LOCALSPECTRALFUNCTION()
!
!     ==========================================================================
!     == ITERATION TO ENFORCE CONSTRAINTS                                     ==
!     ==========================================================================
      ELAST=0.D0
      DO ITER=1,NITERX
        WRITE(*,FMT='(82("="),T20," ITERATION ",I5)')ITER
        CALL DMFT_GLOC() 
!
!       ========================================================================
!       ==  CALL THE SOLVER                                                   ==
!       ========================================================================
        CALL DMFT_SOLVER(ETOT) 
!
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
        CALL DMFT_CONSTRAINTS()
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
!
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
      USE DMFT_MODULE, ONLY: NCHI,NB,NKPTL,NSPIN,NDIM,NDIMD,NAT,NOMEGA,OMEGA &
     &                      ,KSET,ATOMSET,MU
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
!     **                                                                      **
!     ** ENERGY AXIS IN EV AND RELATIVE TO THE FERMI LEVEL                    **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NB,NKPTL,NSPIN,NDIM,NDIMD,NAT,NOMEGA,OMEGA &
     &                      ,KSET,ATOMSET,MU
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TON=.FALSE.
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
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
      INTEGER(4),PARAMETER   :: NE=500
      REAL(8)   ,PARAMETER   :: EMIN=-2.D0/27.211D0
      REAL(8)   ,PARAMETER   :: EMAX=+2.D0/27.211D0
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
PRINT*,'IAT,ICHI1,ICHI2,NLOC,NDIM,NDIMD ',IAT,ICHI1,ICHI2,NLOC,NDIM,NDIMD
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
                    DOS(IE,ICHI1-1+I)=DOS(IE,ICHI1-1+I)+MAT(I,I)
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
        DETOT=REAL(KSET(IKPT)%GAMMA(ICHI1,ICHI2,IDIMD),KIND=8)
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
PRINT*,'ENERGIES ',ETOT
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
      CALL LMTO_OVERLAPEVAL()
!CALL DMFT_REPORTOVERLAP(6)
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
CALL SPINOR_PRINTMATRIX(6,'RHO(T,S)  IN NATORB',1,NLOC,NDIMD,NLOC &
   &    ,RHO)
CALL SPINOR_PRINTMATRIX(6,'SINV(T,S) IN NATORB',1,NLOC,NDIMD,NLOC &
   &    ,SINV)
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
      USE DMFT_MODULE ,ONLY: NKPTL,NCHI,NDIMD,NOMEGA,NLAU,OMEGA,KBT,MU,KSET
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      REAL(8)   ,PARAMETER   :: TOL=1.D-5  ! TOLERANCE FOR THE DENSITY MATRIX
      COMPLEX(8),ALLOCATABLE :: MAT(:,:,:)
      COMPLEX(8),ALLOCATABLE :: A(:,:)
      COMPLEX(8),ALLOCATABLE :: B(:,:)
      REAL(8)   ,ALLOCATABLE :: F(:)   !OCCUPATIONS
      REAL(8)   ,ALLOCATABLE :: E(:)   !ENERGIES
      COMPLEX(8),ALLOCATABLE :: U(:,:)
      INTEGER(4)             :: N
      REAL(8)                :: SVAR
      INTEGER(4)             :: IKPT,ISPIN,I,IE
INTEGER(4) :: IDIMD,ICHI1,ICHI2,ICHI
      INTEGER(4),PARAMETER   :: NE=1000
      REAL(8)   ,PARAMETER   :: EMIN=-10.D0/27.211D0,EMAX=10.D0/27.211D0
      REAL(8)   ,PARAMETER   :: DE=(EMAX-EMIN)/REAL(NE-1,KIND=8)
      REAL(8)                :: HIST(NE)
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
      ALLOCATE(E(N))
      ALLOCATE(U(N,N))
      ALLOCATE(B(N,N))
!
!     ==========================================================================
!     == CALCULATE HRHO                                                       ==
!     ==========================================================================
      HIST=0.D0
      DO IKPT=1,NKPTL
        IF(NDIMD.LE.2) THEN
!         == CONVERT TO UP-DOWN REPRESENTATION =================================
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,KSET(IKPT)%RHO) 
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,KSET(IKPT)%SINV) 
          DO ISPIN=1,NDIMD
!           == A*U=B*U*F  WITH U^DAGGER*S*U=1 ==================================
CALL SPINOR$PRINTMATRIX(6,'RHO IN HRHO','TXYZ',1,NCHI,NDIMD,NCHI,KSET(IKPT)%RHO)
CALL SPINOR$PRINTMATRIX(6,'SINV IN HRHO','TXYZ',1,NCHI,NDIMD,NCHI,KSET(IKPT)%SINV)
            CALL LIB$GENERALEIGENVALUEC8(NCHI,KSET(IKPT)%RHO(:,:,ISPIN) &
     &                                       ,KSET(IKPT)%SINV(:,:,ISPIN),F,U)
!!$WRITE(*,FMT='(78("="),T10," IKPT=",I5," ISPIN=",I2,"  ")')IKPT,ISPIN
!!$WRITE(*,FMT='("OCC=",10F10.5)')F
!!$WRITE(*,FMT='("SUM(OCC)=",10F10.5)')SUM(F)

!!$DO ICHI=1,NCHI
!!$WRITE(*,FMT='(I3,200("(",2F10.5,")"))')ICHI,U(2:6,ICHI)
!!$ENDDO
!
!           == CALCULATE HRHO ==================================================
            B=TRANSPOSE(CONJG(U))
            DO I=1,NCHI
              CALL DMFT_EOFF(KBT,F(I),SVAR)
              SVAR=MAX(-100.D0*KBT,MIN(100.D0*KBT,SVAR))
              SVAR=MU+SVAR
              E(I)=SVAR
              B(I,:)=SVAR*B(I,:)
!!$WRITE(*,FMT='("IN HRHO: E[EV]=",F15.5," OCC=",2F10.5)')E(I)*27.211D0,F(I) &
!!$& ,1.D0/(1.D0+EXP((E(I)-MU)/KBT))
            ENDDO
            KSET(IKPT)%HRHO(:,:,ISPIN)=MATMUL(U,B)
!
!           == CALCULATE N-REPRESENTABLE DENSITY MATRIX ========================
            U=MATMUL(KSET(IKPT)%SINV(:,:,ISPIN),U)
            B=TRANSPOSE(CONJG(U))
            DO I=1,NCHI
!!$IE=NINT((E(I)-EMIN)/DE)
!!$IE=MAX(1,MIN(NE,IE))
!!$HIST(IE)=HIST(IE)+ABS(U(6,I))**2 !ABS(U(7,I))**2
              SVAR=1.D0/(1.D0+EXP((E(I)-MU)/KBT))
              B(I,:)=SVAR*B(I,:)
            ENDDO
            KSET(IKPT)%RHOADAPTED(:,:,ISPIN)=MATMUL(U,B)
!!$WRITE(*,FMT='(78("="),T10," IKPT=",I5," ISPIN=",I2,"  ")')IKPT,ISPIN
!!$WRITE(*,FMT='("E[EV]=",10F10.5)')E*27.211D0
!!$WRITE(*,FMT='("OCCN =",10F10.5)')E*27.211D0
          ENDDO
          CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%RHO)
          CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,KSET(IKPT)%RHOADAPTED)
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


!!$OPEN(20,FILE='HIST.DAT')
!!$DO IE=1,NE
!!$WRITE(20,*)(EMIN+DE*REAL(IE-1,KIND=8))*27.211,HIST(IE)
!!$ENDDO
!!$CLOSE(20)
!!$STOP 'FORCED IN DMFT_HRHO'
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
      SUBROUTINE DMFT_LOCALSPECTRALFUNCTION()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NDIMD,NAT,MU,KSET,ATOMSET
      IMPLICIT NONE
      LOGICAL(4),PARAMETER     :: TON=.FALSE.
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER     :: TPRINT=.FALSE.
      LOGICAL(4),PARAMETER     :: TTEST=.FALSE.
      COMPLEX(8)               :: MATX(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: G(NCHI,NCHI,NDIMD)
      REAL(8)                  :: WKPTL
      INTEGER(4)               :: IKPT,IE,IAT,I1,I2,IDIMD,I
      INTEGER(4)               :: NORB
      REAL(8)   ,PARAMETER :: EMIN=-15.D0/27.211D0,EMAX=15.D0/27.211D0
      INTEGER(4),PARAMETER :: NE=1000
      REAL(8)              :: E(NE)
      REAL(8)   ,PARAMETER :: DELTA=0.1D0/27.211D0
      TYPE SPECTR_TYPE
        REAL(8),ALLOCATABLE :: MAT(:,:,:,:) ! NORB,NORB,NE 
      END TYPE SPECTR_TYPE
      TYPE(SPECTR_TYPE)     :: SPECTR(NAT)
      REAL(8)               :: PI
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_LOCALSPECTRALFUNCTION')
      IF(.NOT.TON) RETURN
      PI=4.D0*ATAN(1.D0)
      DO IAT=1,NAT
        NORB=ATOMSET(IAT)%ICHI2-ATOMSET(IAT)%ICHI1+1
PRINT*,'NORB',NORB,IAT
        ALLOCATE(SPECTR(IAT)%MAT(NORB,NORB,NDIMD,NE))
        SPECTR(IAT)%MAT=(0.D0,0.D0)
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
      DO IE=1,NE
        IAT=2
        NORB=ATOMSET(IAT)%ICHI2-ATOMSET(IAT)%ICHI1+1
        WRITE(*,FMT='(200F15.5)')E(IE)*27.211D0 &
     &                                ,(SPECTR(IAT)%MAT(I,I,1,IE),I=1,NORB)
      ENDDO
!STOP 'FORCED AFTER WRITING SPECTRAL FUNCTION'
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
!IF(IAT.EQ.2)CALL DMFT_SOLVERTEST(IAT)
CALL TESTRHOOFG(NLOC,NDIMD,NOMEGA,NLAU,KBT,ATOMSET(IAT)%GLOC,ATOMSET(IAT)%GLOCLAUR)
!       __SELF ENERGY IS RETURNED AS ATOMSET%DPHIDG_____________________________
        CALL DMFT_DYNAMICSOLVER(NLOC,NDIMD,NOMEGA,NLAU,KBT &
     &                         ,ATOMSET(IAT)%GLOC,ATOMSET(IAT)%GLOCLAUR &
     &                         ,ATOMSET(IAT)%U,PHILW &
     &                         ,ATOMSET(IAT)%DPHIDG,ATOMSET(IAT)%DPHIDGLAUR &
     &                         ,ATOMSET(IAT)%DEDU &
     &                         ,ATOMSET(IAT)%NATORB%PIPHI &
     &                         ,ATOMSET(IAT)%NATORB%CHIPHI)
       CALL SPINOR_PRINTMATRIX(6,'DPHIDG(NU=1) IN SOLVER ' &
    &                      ,1,NLOC,NDIMD,NLOC,ATOMSET(IAT)%DPHIDG(:,:,:,1))
!PRINT*,'DPHIDG ',ATOMSET(IAT)%DPHIDG(1,1,1,:)
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
      DGLOC(1,1,1,1)=(1.D0,2.D0)

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
CALL SPINOR_PRINTMATRIX(6,'SLOC',1,NLOC,NDIMD,NLOC,SLOC(:,:,:,1))

!!$IF(SUM(ABS(SLOC))+SUM(ABS(SLOCLAUR)).GT.1.D-8) THEN
!!$  WRITE(*,FMT='(82("="),T10," AFTER SHRINKDOWN  ")')
!!$  DO IDIMD=1,NDIMD
!!$    WRITE(*,FMT='(82("="),T10," REAL(SIGMA)  ")')
!!$    DO I=1,NLOC
!!$      WRITE(*,FMT='(100F10.5)')REAL(SLOC(I,:,IDIMD,1),8)
!!$    ENDDO
!!$    WRITE(*,FMT='(82("="),T10," IMAG(SIGMA)  ")')
!!$    DO I=1,NLOC
!!$      WRITE(*,FMT='(100F10.5)')AIMAG(SLOC(I,:,IDIMD,1))
!!$    ENDDO
!!$    DO ILAU=1,NLAU
!!$      WRITE(*,FMT='(82("="),T10," REAL(SLAUR(",I5,"))  ")')ILAU
!!$      DO I=1,NLOC
!!$        WRITE(*,FMT='(100F10.5)')REAL(SLOCLAUR(I,:,IDIMD,ILAU),8)
!!$      ENDDO
!!$      WRITE(*,FMT='(82("="),T10," IMAG(SLAUR(",I5,"))  ")')ILAU
!!$      DO I=1,NLOC
!!$        WRITE(*,FMT='(100F10.5)')AIMAG(SLOCLAUR(I,:,IDIMD,ILAU))
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO
!!$ENDIF
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
      LOGICAL(4),SAVE       :: TINI=.FALSE.
      CHARACTER(32),PARAMETER :: IDOUT='DMFT_OUT_GREEN'
      CHARACTER(32),PARAMETER :: IDIN='DMFT_IN_SIGMA'
      INTEGER(4)            :: NORB_,NOMEGA_
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: NU,I
!      CHARACTER(32),PARAMETER :: TYPE='LINEAR' ! 'NONE','LINEAR','STATIC'
!      CHARACTER(32),PARAMETER :: TYPE='NONE' ! 'NONE','LINEAR',
!      CHARACTER(32),PARAMETER :: TYPE='STATIC' ! 'NONE','LINEAR',
      CHARACTER(32),PARAMETER :: TYPE='2NDORDER' ! 'NONE','LINEAR',
      REAL(8)                 :: PI
!     **************************************************************************
                                     CALL TRACE$PUSH('DMFT_SOLVERIO')
      PI=4.D0*ATAN(1.D0)
!CALL TESTG(NORB,NOMEGA,NLAU,KBT,G,GLAUR)
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        CALL FILEHANDLER$SETFILE(IDOUT,.TRUE.,-'.GLOC')
        CALL FILEHANDLER$SETSPECIFICATION(IDOUT,'STATUS','REPLACE')
        CALL FILEHANDLER$SETSPECIFICATION(IDOUT,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(IDOUT,'ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION(IDOUT,'FORM','FORMATTED')
        CALL FILEHANDLER$SETFILE(IDIN,.TRUE.,-'.SIGMA')
        CALL FILEHANDLER$SETSPECIFICATION(IDIN,'STATUS','OLD')
        CALL FILEHANDLER$SETSPECIFICATION(IDIN,'POSITION','REWIND')
        CALL FILEHANDLER$SETSPECIFICATION(IDIN,'ACTION','READ')
        CALL FILEHANDLER$SETSPECIFICATION(IDIN,'FORM','FORMATTED')
      END IF
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
      ELSE IF(TYPE.EQ.'STATIC') THEN
        PRINT*,'USING HF AS LUTTINGER-WARD FUNCTIONAL'
        PRINT*,'(USE FOR TESTING PURPOSES ONLY)'
        CALL DMFT_SOLVER_STATIC(NORB,NOMEGA,KBT,G,U,ETOT,S)  
        CALL TRACE$POP()
        RETURN
      ELSE IF(TYPE.EQ.'2NDORDER') THEN
!IF(NORB.NE.2) THEN
IF(NORB.NE.12) THEN
 S=(0.D0,0.D0)  
 SLAUR=(0.D0,0.D0)
 RETURN
END IF
        CALL MSIGMA$2NDORDER(NORB,NOMEGA,NLAU,KBT,G,GLAUR,U,ETOT,S,SLAUR)
!!$WRITE(*,FMT='(82("="),T10,"  RIGHT AFTER MSIGMA$2NDORDER  ")')
!!$WRITE(*,FMT='(82("="),T10," REAL(SIGMA)  ")')
!!$DO I=1,NORB
!!$  WRITE(*,FMT='(100F10.5)')REAL(S(I,:,1),8)
!!$ENDDO
!!$WRITE(*,FMT='(82("="),T10," IMAG(SIGMA)  ")')
!!$DO I=1,NORB
!!$  WRITE(*,FMT='(100F10.5)')AIMAG(S(I,:,1))
!!$ENDDO
!!$DO NU=1,NLAU
!!$  WRITE(*,FMT='(82("="),T10," REAL(SLAUR(",I5,"))  ")')NU
!!$  DO I=1,NORB
!!$    WRITE(*,FMT='(100F10.5)')REAL(SLAUR(I,:,NU),8)
!!$  ENDDO
!!$  WRITE(*,FMT='(82("="),T10," IMAG(SLAUR(",I5,"))  ")')NU
!!$  DO I=1,NORB
!!$    WRITE(*,FMT='(100F10.5)')AIMAG(SLAUR(I,:,NU))
!!$  ENDDO
!!$ENDDO
        RETURN
      END IF
!
!     ==========================================================================
!     == WRITE DATA
!     ==========================================================================
!WARNING!! THE TWO FILES ARE NOT YET DEFINED
!CAREFUL: IN FORTRAN, THE FIRST INDEX RUNS FIRST
      CALL FILEHANDLER$UNIT(IDOUT,NFIL)
      REWIND NFIL
      DO NU=1,NOMEGA
        WRITE(NFIL,FMT='(200F20.5)')REAL(2*NU-1,KIND=8)*PI*KBT &
     &                             ,(G(I,I,NU),I=1,NORB)
      ENDDO

!!$      WRITE(NFIL,*)NORB,NOMEGA,NLAU,KBT
!!$      WRITE(NFIL,FMT='(2F20.10)')G
!!$      WRITE(NFIL,FMT='(2F20.10)')GLAUR
!!$      WRITE(NFIL,FMT='(2F20.10)')U
      CALL FILEHANDLER$CLOSE(IDOUT)
STOP 'STOP AFTER WRITEING DMFTOSOLVER'
!
!     ==========================================================================
!     == READ DATA DATA
!     ==========================================================================
      CALL FILEHANDLER$UNIT(IDIN,NFIL)
      REWIND NFIL
      READ(NFIL,*)NORB_,NOMEGA_,ETOT
      READ(NFIL,FMT='(2F20.10)')S
      READ(NFIL,FMT='(2F20.10)')SLAUR
      READ(NFIL,FMT='(2F20.10)')DEDU
      CALL FILEHANDLER$CLOSE(IDOUT)
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
      REAL(8)   ,PARAMETER  :: A=7.D-3,B=0.D0,LAMBDA=1.D-1
!      REAL(8)   ,PARAMETER  :: A=0.D0,B=1.D-1,LAMBDA=1.D-1
!      REAL(8)   ,PARAMETER  :: A=-7.D-3,B=1.D0,LAMBDA=1.D-1
      REAL(8)               :: SVAR
      COMPLEX(8)            :: CMAT(NORB,NORB)
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      S=(0.D0,0.D0)
      ETOT=0.D0
      DO NU=1,NOMEGA
        OMEGA=REAL(2*NU-1,KIND=8)*PI*KBT
!       == DEFINE SELF ENERGY ==================================================
        DO I=1,NORB
!          S(I,I,NU)=0.007D0*EXP(-0.1D0/ABS(OMEGA))/(CI*OMEGA)
          S(I,I,NU)=EXP(-LAMBDA/ABS(OMEGA))*(A/(CI*OMEGA)+B*(1.D0,0.D0))
        ENDDO
!       == CALCULATE ENERGY=====================================================
        CMAT=MATMUL(S(:,:,NU),G(:,:,NU))
        SVAR=0.D0
        DO I=1,NORB
          SVAR=SVAR+2.D0*REAL(CMAT(I,I),KIND=8) ! FACTOR 2 FROM -OMEGA_NU 
        ENDDO
        ETOT=ETOT+KBT*SVAR
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
      IMPLICIT NONE
      INTEGER(4),PARAMETER     :: NITER1=100  ! X#(CONJUGATE GRADIENT STEPS)
      INTEGER(4),PARAMETER     :: NITER2=100  ! X#(STEPS) FOR LINE SARCH
      REAL(8)   ,PARAMETER     :: TOL1=1.D-5  ! TOLERANCE FOR CG ITERATION
      REAL(8)   ,PARAMETER     :: TOL2=1.D-5  ! TOLERANCE FOR LINE SEARCH
      REAL(8)   ,PARAMETER     :: SCALEINI=1.D-5  ! INITIAL MIXING
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: DGAMMA(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: FORCE(NCHI,NCHI,NDIMD)
      COMPLEX(8)               :: DRHO(NCHI,NCHI,NDIMD)
      LOGICAL(4)               :: CONVG
      REAL(8)                  :: MAXDEV,MAXDEVOLD
      REAL(8)                  :: SCALE
      REAL(8)                  :: FF,PF,PFOLD
      INTEGER(4)               :: IKPT,IDIMD
      INTEGER(4)               :: ITER1,ITER2,ITER,I,J
      REAL(8)                  :: FN(2)
      REAL(8)                  :: SVAR,SVAR1,SVAR2
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
!     ==  DEVIATION FROM TARGET DENSITY MATRIX                                ==
!     ==========================================================================
      DO IKPT=1,NKPTL
!       
!       ========================================================================
!       ==  ADJUST H0 TO OBEY DENSITY MATRIX CONSTRAINT                       ==
!       ========================================================================
        DO ITER1=1,NITER1
          IF(ITER1.EQ.1) THEN
            CALL DMFT_DRHO(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                  ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
     &                  ,KSET(IKPT)%GAMMA,DRHO)
!
!           == DGAMMA=DMAXDEV/DGAMMA ===========================================
            CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                  ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
     &                  ,KSET(IKPT)%GAMMA,DRHO,MAXDEV,FORCE)
            FORCE=-FORCE  ! DOWNHILL DIRECTION
            DGAMMA=(0.D0,0.D0)
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
!         == CHECK CONVERGENCE                                                ==
!         ======================================================================
          CONVG=SQRT(MAXDEV).LT.TOL1
          IF(CONVG) EXIT
!
!         ======================================================================
!         == PREPARE NEXT SEARCH DIRECTION DGAMMA                             **
!         ======================================================================
          FF=SUM(ABS(FORCE)**2)
          DGAMMA=DGAMMA+FORCE/FF   ! NEW SEARCH DIRECTION
!!$DO IDIMD=1,NDIMD
!!$  DO I=1,NCHI
!!$    DO J=1,NCHI
!!$      SVAR1=MAXVAL(ABS(DRHO))/10.D0
!!$      SVAR2=MAXVAL(ABS(DGAMMA))/10.D0
!!$      IF(ABS(DRHO(I,J,IDIMD)).GT.SVAR1.OR.ABS(DGAMMA(I,J,IDIMD)).GT.SVAR2) THEN
!!$        WRITE(*,FMT='(3I4,10F10.5)')I,J,IDIMD,DRHO(I,J,IDIMD),DGAMMA(I,J,IDIMD)
!!$      END IF
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
!
!         ======================================================================
!         == LINE SEARCH                                                      ==
!         ======================================================================
          PFOLD=REAL(SUM(CONJG(FORCE)*DGAMMA),KIND=8)
          MAXDEVOLD=MAXDEV
          SCALE=SCALEINI*FF
          DO ITER2=1,NITER2
            KSET(IKPT)%GAMMA=KSET(IKPT)%GAMMA+DGAMMA*SCALE
!
            CALL DMFT_DRHO(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                    ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
     &                    ,KSET(IKPT)%GAMMA,DRHO)
            CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                  ,KSET(IKPT)%SMAT,KSET(IKPT)%SINV,KSET(IKPT)%HRHO &
     &                  ,KSET(IKPT)%GAMMA,DRHO,MAXDEV,FORCE)
            FORCE=-FORCE  ! DOWNHILL DIRECTION
            PF=REAL(SUM(CONJG(FORCE)*DGAMMA),KIND=8)
            CONVG=ABS(PF).LT.TOL2
            IF(CONVG) EXIT
            SVAR=PF/(PFOLD-PF)
!           __LIMIT UNUSUALLY LARGE STEPS_______________________________________
            IF(0.5D0*ABS(SVAR).GT.1.D0)SVAR=SVAR/(0.5D0*ABS(SVAR))
            SCALE=SVAR*SCALE
            PFOLD=PF
            MAXDEVOLD=MAXDEV
          ENDDO ! END OF LINE SEARCH
          IF(.NOT.CONVG) THEN
            CALL ERROR$MSG('LINE SEARCH NOT CONVERGED')
            CALL ERROR$R8VAL('CONV CITERIUM',PF)
            CALL ERROR$R8VAL('TOLERANCE',TOL2)
            CALL ERROR$STOP('DMFT_CONSTRAINTS')
          END IF
PRINT*,'ITER ',ITER1,SQRT(MAXDEV)
          CONVG=.FALSE.
        ENDDO ! END OF CONJUGATE GRADIENT LOOP
        IF(.NOT.CONVG) THEN
          CALL ERROR$MSG('LOOP NOT CONVERGED')
          CALL ERROR$I4VAL('IKPT',IKPT)
          CALL ERROR$I4VAL('ITER1',ITER1)
          CALL ERROR$R8VAL('MAX. DEVIATION',MAXDEV)
          CALL ERROR$R8VAL('TOLERANCE',TOL1)
          CALL ERROR$STOP('DMFT_CONSTRAINTS')
        END IF
      ENDDO   !END OF LOOP OVER K-POINTS
                              CALL TRACE$POP()
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
      COMPLEX(4),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: DRHO(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: GAMMA(NCHI,NCHI,NDIMD)
      REAL(8)               :: MAXDEV0,MAXDEVP,MAXDEVM
      COMPLEX(8)            :: FORCE0(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: FORCEP(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: FORCEM(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: DIR(NCHI,NCHI,NDIMD)
      REAL(8)               :: SVARP,SVARM
      INTEGER(4)            :: IDIS
!     **************************************************************************
      CALL DMFT_DRHO(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA0,DRHO)
      CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA0,DRHO,MAXDEV0,FORCE0)
      FORCE0=-FORCE0
      DIR=FORCE0
      DIR=DIR/SQRT(SUM(ABS(DIR)**2))
      PRINT*,'--',0.D0,-REAL(SUM(DIR*FORCE0),KIND=8)
      DO IDIS=1,5
        SVARP=+1.D-5*REAL(IDIS,8)
        SVARM=-SVARP
        GAMMA=GAMMA0+SVARP*DIR
        CALL DMFT_DRHO(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA,DRHO)
        CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA,DRHO,MAXDEVP,FORCEP)
        GAMMA=GAMMA0+SVARM*DIR
        CALL DMFT_DRHO(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA,DRHO)
        CALL DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &              ,SMAT,SINV,HRHO,GAMMA,DRHO,MAXDEVM,FORCEM)
        PRINT*,'--',SVARP,(MAXDEVP-MAXDEVM)/(SVARP-SVARM)
      ENDDO
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
      COMPLEX(4),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: GLAUR(NCHI,NCHI,NDIMD,NLAU+1)
      COMPLEX(8)            :: MATX(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: DMATX(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: G(NCHI,NCHI,NDIMD)
      INTEGER(4)            :: NU,IAT,I1,I2,IDIMD,ILAU,I
!     **************************************************************************
!
!     == LAURENT EXPANSION FOR THE GREENS FUNCTION DIFFERENCE ==============
      GLAUR=(0.D0,0.D0)
!     == MATX=-MU*S+HRHO ; DMATX=SLAUR(0)-GAMMA ============================
!      MATX=-MU*SMAT+HRHO
!     == SELF ENERGY MINUS GAMMA
      DMATX=(0.D0,0.D0)
      DO IAT=1,NAT
        I1=ATOMSET(IAT)%ICHI1
        I2=ATOMSET(IAT)%ICHI2
        DMATX(I1:I2,I1:I2,:)=DMATX(I1:I2,I1:I2,:)+ATOMSET(IAT)%SLOCLAUR(:,:,:,1)
      ENDDO
      DMATX=DMATX-GAMMA
!     == GLAUR2=SINV*MATX*SINV =============================================
      CALL SPINOR$MATMUL(NDIMD,NCHI,DMATX,SINV,MAT)
      CALL SPINOR$MATMUL(NDIMD,NCHI,SINV,MAT,GLAUR(:,:,:,2))
!     == THE NEXT TERM GLAUR3 DOES NOT CONTRIBUTE TO THE DENSITY MATRIX. ===
!     == GLAUR4 IS NOT IMPLEMENTED BECAUSE NLAU<3 ==========================
!
!     == LOOP OVER OMEGA ===================================================
      DRHO=(0.D0,0.D0)
      DO NU=1,NOMEGA
!       == CONSTRUCT LATTICE GREENS FUNCTION =============================
        MAT=(CI*OMEGA(NU)+MU)*SMAT-HRHO 
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
        DRHO=DRHO-KBT*G  ! SUBTRACT NONINTERACTING G
        MAT=MAT+GAMMA
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          MAT(I1:I2,I1:I2,:)=MAT(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
        ENDDO
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
        DRHO=DRHO+KBT*G
      ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!
!     == INCLUDE NEGATIVE FREQUENCIES ======================================
      DO IDIMD=1,NDIMD
        DRHO(:,:,IDIMD)=DRHO(:,:,IDIMD)+CONJG(TRANSPOSE(DRHO(:,:,IDIMD)))
      ENDDO
!
!     == ADD TAILS (GLAUR3 DOES NOT CONTRIBUTE) ============================
!SWITCHED OFF BECAUSE THEY ARE NOT ACCOUNTED FOR IN THE DERIVATIVES
!!$      DO ILAU=1,NLAU+1
!!$        DRHO=DRHO+FN(ILAU)*GLAUR(:,:,:,ILAU)
!!$      ENDDO


!!$DO IDIMD=1,NDIMD
!!$  WRITE(*,FMT='(82("="),T10," REAL(DRHO) IDIMD=",I5," ")')IDIMD
!!$  DO I=1,NCHI
!!$    WRITE(*,FMT='(100F10.5)')REAL(DRHO(I,:,IDIMD),8)
!!$  ENDDO
!!$  WRITE(*,FMT='(82("="),T10," IMAG(DRHO) ",I5," ")')IDIMD
!!$  DO I=1,NCHI
!!$    WRITE(*,FMT='(100F10.5)')AIMAG(DRHO(I,:,IDIMD))
!!$  ENDDO
!!$ENDDO
!!$STOP 'FORCED'

!!$SVAR=0.D0
!!$DO IDIMD=1,NDIMD
!!$  SVAR=MAX(SVAR,MAXVAL(ABS(DRHO(:,:,IDIMD)-CONJG(TRANSPOSE(DRHO(:,:,IDIMD))))))
!!$ENDDO
!!$PRINT*,'DEV FROM HERMITEAN DRHO',SVAR
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_SUMDRHO2(NCHI,NDIMD,NOMEGA,NLAU,KBT,MU,OMEGA,FN &
     &                        ,SMAT,SINV,HRHO,GAMMA,DRHO,F,DFDGAMMA)
!     **************************************************************************
!     ** F=SUM(ABS(DRHO)**2)=TR[DRHO*DRHO^\DAGGER]
!     ** DF=DFDGAMMA*DGAMMA
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
      COMPLEX(8),INTENT(IN) :: DRHO(NCHI,NCHI,NDIMD)
      REAL(8)   ,INTENT(OUT):: F
      COMPLEX(8),INTENT(OUT):: DFDGAMMA(NCHI,NCHI,NDIMD)
      COMPLEX(4),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: GLAUR(NCHI,NCHI,NDIMD,NLAU+1)
      COMPLEX(8)            :: MATX(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: DMATX(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: MAT(NCHI,NCHI,NDIMD)
      COMPLEX(8)            :: G(NCHI,NCHI,NDIMD)
      INTEGER(4)            :: NU,IAT,I1,I2,IDIMD,ILAU,I
      REAL(8)               :: SVAR
!     **************************************************************************
      F=SUM(ABS(DRHO)**2)
!
!     == LOOP OVER OMEGA =======================================================
      DFDGAMMA=(0.D0,0.D0)
      DO NU=1,NOMEGA
!       == CONSTRUCT LATTICE GREENS FUNCTION ===================================
        MAT=(CI*OMEGA(NU)+MU)*SMAT-HRHO+GAMMA
        DO IAT=1,NAT
          I1=ATOMSET(IAT)%ICHI1
          I2=ATOMSET(IAT)%ICHI2
          MAT(I1:I2,I1:I2,:)=MAT(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
        ENDDO
        CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
        DO IDIMD=1,NDIMD
          DFDGAMMA(:,:,IDIMD)=DFDGAMMA(:,:,IDIMD) &
     &            -KBT*MATMUL(G(:,:,IDIMD),MATMUL(DRHO(:,:,IDIMD),G(:,:,IDIMD)))
        ENDDO
      ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!
!     == INCLUDE NEGATIVE FREQUENCIES ==========================================
      DO IDIMD=1,NDIMD
        DFDGAMMA(:,:,IDIMD)=DFDGAMMA(:,:,IDIMD) &
     &                     +CONJG(TRANSPOSE(DFDGAMMA(:,:,IDIMD)))
      ENDDO
      DFDGAMMA=2.D0*DFDGAMMA  
!
!     == THIS IS A FUDGE FACTOR, WHICH IS NOT YET UNDERSTOOD. I BELIEVE THAT
!     == IT HAS TO DO WITH THE SPINOR REPRESENTSTION. IT WAS NUMERICALLY 
!     == REQUIRED WITH NON-SPIN POLARIZED CALCULATKION
      DFDGAMMA=0.25D0*DFDGAMMA  
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
      REAL(8)   ,PARAMETER :: MIX1=0.1D0
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
!SVAR=MIX1
          ATOMSET(IAT)%SLOC(:,:,:,NU)=ATOMSET(IAT)%SLOC(:,:,:,NU) &
     &                          +SVAR*ATOMSET(IAT)%DPHIDG(:,:,:,NU)
        ENDDO
!
!       == MIX IN WITH FREQUENCY INDEPENDENT MASS ==============================
        SVAR=1.D0
! SVAR=MIX1
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
        ETOT=ETOT+EDC
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
      REAL(8)                :: XSUM
      REAL(8)                :: ETOT1
      INTEGER(4)             :: NU,I,IAT,I1,I2,IDIMD,IKPT
!     **************************************************************************
      IF(NDIMD.EQ.4) THEN
        ALLOCATE(BMAT1(2*NCHI,2*NCHI))
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
            MAT1(I,I,1)=-MAT1(I,I,1)+(2.D0,0.D0) !(UNIT IN SPINOR REPRESENTATION)
          ENDDO
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,MAT1)
          IF(NDIMD.EQ.4) THEN
            BMAT1(:NCHI,:NCHI)    =MAT1(:,:,1)
            BMAT1(NCHI+1:,:NCHI)  =MAT1(:,:,2)
            BMAT1(:NCHI,NCHI+1:)  =MAT1(:,:,3)
            BMAT1(NCHI+1:,NCHI+1:)=MAT1(:,:,4)
            CALL LIB$EIGVALNONHERMITEANC8(2*NCHI,BMAT1,CVEC)
            XSUM=0.D0
            DO I=1,2*NCHI
              XSUM=XSUM+LOG(ABS(CVEC(I)))
            ENDDO
            ETOT1=ETOT1+XSUM
          ELSE
            DO IDIMD=1,NDIMD
              CALL LIB$EIGVALNONHERMITEANC8(NCHI,MAT1(:,:,IDIMD),CVEC)
              XSUM=0.D0
              DO I=1,NCHI
                XSUM=XSUM+2.D0*LOG(ABS(CVEC(I))) ! FACTOR TWO FROM -NU
              ENDDO
              IF(NDIMD.EQ.1) XSUM=XSUM*2.D0 !SPIN DEGENERACY
              ETOT1=ETOT1+XSUM
            ENDDO
          END IF
        ENDDO ! END OF LOOP OVER MATSUBARA FREQUENCIES
        ETOT=ETOT+KSET(IKPT)%WKPT*ETOT1
      ENDDO ! END OF LOOK OVER K-POINTS
      ETOT=-KBT*ETOT
      WRITE(*,FMT='(60("."),T1'&
     &        //',"--DMFT--ENERGY FROM MATSUBARA SUM"' &
     &        //',T60,":",F20.10)')ETOT
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
      REAL(8)               :: FN(3)
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
      INTEGER(4)               :: IKPT,NU,IAT,I1,I2,IDIMD,ILAU
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
     &                          ,GLAU_IN,UTENSOR,LWENERGY,SIGMA,SLAU_OUT)
!     **************************************************************************
!     ** THIS FUNCTION SHOULD BE CALLED FROM OUTSIDE                          **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NORB
      INTEGER(4),INTENT(IN)  :: NOMEGA
      INTEGER(4),INTENT(IN)  :: NLAU
      REAL(8)   ,INTENT(IN)  :: KBT
      COMPLEX(8),INTENT(IN)  :: GF(NORB,NORB,NOMEGA)
      COMPLEX(8),INTENT(IN)  :: GLAU_IN(NORB,NORB,NLAU+1)
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
      REAL(8)                :: EV
      REAL(8)                :: PI
!     **************************************************************************
      CALL CONSTANTS$GET('EV',EV)
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
        CALL MSIGMA_WRITEGF('G',NORB,NOMEGA,OMEGA,GF,"GF_IN.DAT")
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
      IF(TPR) THEN
!!$        PRINT*,'SUM(GFFULL) ',SUM(GFFULL)
!!$        OPEN(888,FILE='GFFFULL.DAT')
!!$        DO I=1,NMAXINDEX-NMININDEX+1
!!$          WRITE(888,FMT='(I5,200(F10.5,2X))')I,(GFFULL(J,J,I),J=1,NORB)
!!$       ENDDO
!!$       CLOSE(888)
      END IF
!
!     ==========================================================================
!     == CALCULATE SELF ENERGY                                                ==
!     ==========================================================================
      SIGMA = (0.D0,0.D0)
      DO N=1,NOMEGA
        WRITE(*,'(F6.2,A6)') (N-1)*100.D0/NOMEGA,'% DONE'
        DO M=1,NORB
!         __SIGMA(M,M,N+1) = GETSIGMA(M,N,UTENSOR)______________________________
          CALL MSIGMA_GETSIGMA(NORB,SUM_MAX,NMININDEX,NMAXINDEX &
    &                          ,M,N,BETA,UTENSOR,GFFULL,SIGMA(M,M,N))
        END DO
!WRITE(*,FMT='(I5,200F10.5)')N,(SIGMA(M,M,N),M=1,NORB)
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
        CALL MSIGMA_WRITEGF('S',NORB,NOMEGA,OMEGA,SIGMA,'SIGMA.DAT') 
      END IF
      RETURN
      END SUBROUTINE MSIGMA$2NDORDER
!   
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MSIGMA_WRITEGF(ID,NORB,NOMEGA,OMEGA,GF,FILENAME)
!     **************************************************************************
!     ** WRITE GREENS FUNCTION OR SELF ENERGY TO FILE
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: NORB
      INTEGER(4)  ,INTENT(IN) :: NOMEGA
      REAL(8)     ,INTENT(IN) :: OMEGA(NOMEGA)
      COMPLEX(8)  ,INTENT(IN) :: GF(NORB,NORB,NOMEGA)
      CHARACTER(*),INTENT(IN) :: FILENAME
      INTEGER(4)  ,PARAMETER  :: IOPARAM = 1999
      INTEGER(4)  ,PARAMETER  :: IOPARAMEV = 2000
      INTEGER(4)              :: M,N
      REAL(8)                 :: EV
      REAL(8)                 :: FACTOR
!     **************************************************************************
      CALL CONSTANTS$GET('EV',EV)
      IF(ID.EQ.'G') THEN
        FACTOR=1.D0/EV
      ELSE IF(ID.EQ.'S') THEN
        FACTOR=EV
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED (MUST BE "G" OR "S")')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('MSIGMA_WRITEGF')
      END IF
      OPEN(UNIT=IOPARAM,FILE=FILENAME)
      OPEN(UNIT=IOPARAMEV,FILE=FILENAME//"_EV")
      DO N=1,NOMEGA
        WRITE(IOPARAM,'(ES24.16E3,5X)',ADVANCE='NO') OMEGA(N)
        WRITE(IOPARAMEV,'(ES24.16E3,5X)',ADVANCE='NO')OMEGA(N)*EV
        DO M=1,NORB
          WRITE(IOPARAM,'(ES24.16E3,3X,ES24.16E3,5X)',ADVANCE='NO') &
     &                      REAL( GF(M,M,N),KIND=8), AIMAG( GF(M,M,N))
          WRITE(IOPARAMEV,'(ES24.16E3,3X,ES24.16E3,5X)',ADVANCE='NO') &
     &                      REAL( GF(M,M,N),KIND=8)*FACTOR &
     &                     ,AIMAG( GF(M,M,N))*FACTOR
        ENDDO !END LOOP OVER ORBITALS
        WRITE(IOPARAM,*)
        WRITE(IOPARAMEV,*)
      ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
      CLOSE(UNIT=IOPARAM)
      CLOSE(UNIT=IOPARAMEV)
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
          FILL=FILL+REAL(GF(ORB,ORB,NU),KIND=8)+GLAU(ORB,ORB,3)/OMEGA(NU)**2
        END DO
        FILL = FILL*2.0/BETA 
        FILL=FILL + GLAU(ORB,ORB,2)*0.5D0 - GLAU(ORB,ORB,3)*BETA/4.0
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
      INTEGER(4)            :: N,ORBITAL,ILAU,NU
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
!     ** CALCULATE SIGMA IN 2ND ORDER
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
          S = S - REAL(U(ORB,M,ORB,M),KIND=8)**2 * S_PRE_N1/(2.0*BETA*BETA)
        END IF
      END DO
      RETURN
      END SUBROUTINE MSIGMA_GETSIGMA
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE MSIGMA_GETLW(NORB,NOMEGA,BETA,GF,SIGMA,LW)
!     **************************************************************************
!     ** CALCULATE LUTTINGER-WARD ENERGY BY G*SIGMA                           **
!     **  REAL(8) FUNCTION GETLW(GF,SIGMA) RESULT(LW)                     **
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
          LW = LW +  REAL(GF(M,M,N),KIND=8)* REAL(SIGMA(M,M,N),KIND=8) & 
     &            - AIMAG(GF(M,M,N))       *AIMAG(SIGMA(M,M,N))    
         END DO
      END DO
      LW=LW/(4.D0*BETA)
      WRITE(6,'(A11,ES14.6E3)') 'LW-ENERGY: ',LW
      RETURN
      END SUBROUTINE MSIGMA_GETLW
   
   
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
