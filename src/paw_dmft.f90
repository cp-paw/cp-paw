!
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
LOGICAL(4),PARAMETER   :: TON=.true.
LOGICAL(4),SAVE        :: TINI=.FALSE.
INTEGER(4)             :: NOMEGA
INTEGER(4)             :: NCHI          ! #(CORRELATED ORBITALS)
INTEGER(4)             :: NB            ! #(BAND STATES PER K-POINT)
INTEGER(4)             :: NKPTL         ! #(KPOINTS ON THIS TASK)
INTEGER(4)             :: NSPIN         ! #(SPIN COMPONENTS)
INTEGER(4)             :: NDIM          ! #(SPINOR COMPONENTS)
REAL(8)   ,ALLOCATABLE :: OMEGA(:)      ! MATSUBARA FREQUENCIES
INTEGER(4),ALLOCATABLE :: IPROOFCHI(:)
REAL(8)                :: KBT           ! TEMPERATURE (K_B*T)
REAL(8)                :: MU            ! CHEMICAL POTENTIAL
REAL(8)                :: DELTAT        ! TIMESTEP
REAL(8)   ,ALLOCATABLE :: WKPTL(:)      ! (NKPTL) K-POINT WEIGHTS ON THIS TASK
COMPLEX(8),ALLOCATABLE :: PIPSI(:,:,:,:)    !(NCHI,NB,NKPTL,NSPIN) <PI|PSI>
COMPLEX(8),ALLOCATABLE :: DENMAT(:,:,:,:)   !(NB,NB,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: HAMILTON(:,:,:,:) !(NB,NB,NKPTL,NSPIN)<PSI|HEFF|PSI>
real(8)   ,ALLOCATABLE :: erho(:,:,:)       !(NB,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: dedrho(:,:,:,:)   !(Nchi,Nchi,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: h0(:,:,:,:)       !(Nchi,Nchi,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: hrho(:,:,:,:)     !(Nchi,Nchi,NKPTL,NSPIN)
!
!== LOCAL GREENS FUNCTION (WITH LAURENT EXPANSION COEFFICIENTS) ================
COMPLEX(8),ALLOCATABLE :: GLOC(:,:,:,:)     !(NCHI,NCHI,NOMEGA,NSPIN)
COMPLEX(8),ALLOCATABLE :: GLOCLAUR(:,:,:,:)  !(NCHI,NCHI,3,NSPIN)
COMPLEX(8),ALLOCATABLE :: GLOCLAUR1(:,:,:)  !(NCHI,NCHI,NSPIN)
COMPLEX(8),ALLOCATABLE :: GLOCLAUR2(:,:,:)  !(NCHI,NCHI,NSPIN)
COMPLEX(8),ALLOCATABLE :: GLOCLAUR3(:,:,:)  !(NCHI,NCHI,NSPIN)
!__MATRIX TO DIAGONALIZE OVERLAP MATRIX_________________________________________
COMPLEX(8),ALLOCATABLE :: DIAGSLOC(:,:,:) !(NCHI,NCHI,NSPIN)
COMPLEX(8),ALLOCATABLE :: smat(:,:,:)     !(NCHI,NCHI,NSPIN)
COMPLEX(8),ALLOCATABLE :: sinv(:,:,:,:)   !(NCHI,NCHI,nkptl,NSPIN)
!
!== SELF ENERGY ================================================================
COMPLEX(8),ALLOCATABLE :: sigma(:,:,:,:) !(NCHI,NCHI,NOMEGA,NSPIN)
COMPLEX(8),ALLOCATABLE :: sigmadc(:,:,:) !(NCHI,NCHI,NSPIN) (double counting)
COMPLEX(8),ALLOCATABLE :: sigmalaur1(:,:,:) !(NCHI,NCHI,NSPIN)
COMPLEX(8),ALLOCATABLE :: sigmalaur2(:,:,:) !(NCHI,NCHI,NSPIN)
COMPLEX(8),ALLOCATABLE :: sigmalaur3(:,:,:) !(NCHI,NCHI,NSPIN)
!__SELF ENERGY minus double counting ___________________________________________
COMPLEX(8),ALLOCATABLE :: DSIG0(:,:,:,:) !(NCHI,NCHI,NOMEGA,NSPIN)
!__LAURENT EXPANSION FOR DSIGMA_________________________________________________
COMPLEX(8),ALLOCATABLE :: DSIGLAUR0(:,:,:,:)    !(NCHI,NCHI,NSPIN,2)
!
!== LAGRANGE MULTIPLIER FOR DENSITY MATRIX CONSTRAINT ==========================
COMPLEX(8),ALLOCATABLE :: GAMMA0(:,:,:,:) !(NB,NB,NKPTL,NSPIN)
CONTAINS
END MODULE DMFT_MODULE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_INI()
!     **************************************************************************
      USE DMFT_MODULE
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR,LNX,LOX
      USE WAVES_MODULE, ONLY : KMAP,NDIM_W=>NDIM,NKPTL_W=>NKPTL,NSPIN_W=>NSPIN
      IMPLICIT NONE
      REAL(8)                :: PI
      INTEGER(4)             :: ISP0
      INTEGER(4)             :: L0
      REAL(8)   ,ALLOCATABLE :: WKPT(:)    ! K-POINT WEIGHTS (GLOBAL)
      INTEGER(4)             :: NTASKS_K,THISTASK_K
      INTEGER(4)             :: NTASKS_M,THISTASK_M
      INTEGER(4)             :: NKPT
      INTEGER(4)             :: NU,ISP,LN,IPRO,IKPTL,IKPT
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
!     == COLLECT PERMANENT DATA                                               ==
!     ==========================================================================
      NDIM=NDIM_W   ! IMPORT NDIM FROM WAVES_MODULE INTO DMFT MODULE
      IF(NDIM.NE.1) THEN
        CALL ERROR$MSG('NON-COLLINEAR OPTION NOT IMPLEMENTED')
        CALL ERROR$STOP('DMFT_INI')
      END IF
!
      NKPTL=NKPTL_W ! IMPORT NKPTL FROM WAVES_MODULE INTO DMFT MODULE
      NSPIN=NSPIN_W ! IMPORT NSPIN FROM WAVES_MODULE INTO DMFT MODULE
      CALL DYNOCC$GETI4('NB',NB)
      CALL DYNOCC$GETI4('NKPT',NKPT)
      ALLOCATE(WKPT(NKPT))
      ALLOCATE(WKPTL(NKPTL))
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
          WKPTL(IKPTL)=WKPT(IKPT)
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
!     == SELECT CORRELATED ORBITALS                                           ==
!     ==========================================================================
!THIS IS SPECIFIC FOR SRVO3
      ISP0=2
      L0=2
      NCHI=3    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IPRO=0
      DO ISP=1,ISP0-1
        IPRO=IPRO+SUM(2*LOX(:LNX(ISP),ISP)+1)
      ENDDO
      DO LN=1,LNX(ISP0)
        IF(LOX(LN,ISP0).EQ.L0) EXIT
        IPRO=IPRO+2*LOX(LN,ISP0)+1
      ENDDO
      ALLOCATE(IPROOFCHI(NCHI))
      IPROOFCHI(1)=IPRO+2   !T2G ORBITALS
      IPROOFCHI(2)=IPRO+4
      IPROOFCHI(3)=IPRO+5
!
!     ==========================================================================
!     == OTHER VARIABLES                                                      ==
!     ==========================================================================
!THIS IS HARD WIRED
      NOMEGA=200
      KBT=0.333D0*EV
      DELTAT=10.D0
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
!     ==  ALLOCATE PERMANENT ARRAYS
!     ==========================================================================
      ALLOCATE(sigma(NCHI,NCHI,NOMEGA,NSPIN))
      ALLOCATE(sigmalaur1(NCHI,NCHI,NSPIN))
      ALLOCATE(sigmalaur2(NCHI,NCHI,NSPIN))
      ALLOCATE(sigmalaur3(NCHI,NCHI,NSPIN))
      ALLOCATE(sigmadc(NCHI,NCHI,NSPIN))

      ALLOCATE(DSIG0(NCHI,NCHI,NOMEGA,NSPIN))
      ALLOCATE(DSIGLAUR0(NCHI,NCHI,NSPIN,2))
      ALLOCATE(GAMMA0(NB,NB,NKPT,NSPIN))

      ALLOCATE(GLOC(NCHI,NCHI,NOMEGA,NSPIN))
      ALLOCATE(DIAGSLOC(NCHI,NCHI,NSPIN))
      ALLOCATE(smat(NCHI,NCHI,NSPIN))
      ALLOCATE(sinv(NCHI,NCHI,nkptl,NSPIN))
      ALLOCATE(GLOCLAUR(NCHI,NCHI,3,NSPIN))
      ALLOCATE(GLOCLAUR1(NCHI,NCHI,NSPIN))
      ALLOCATE(GLOCLAUR2(NCHI,NCHI,NSPIN))
      ALLOCATE(GLOCLAUR3(NCHI,NCHI,NSPIN))
      ALLOCATE(DENMAT(NB,NB,NKPTL,NSPIN))
      ALLOCATE(HAMILTON(NB,NB,NKPTL,NSPIN))
      ALLOCATE(erho(NB,NKPTL,NSPIN))
      ALLOCATE(PIPSI(NCHI,NB,NKPTL,NSPIN))
      ALLOCATE(dedrho(NCHI,Nchi,NKPTL,NSPIN))
      ALLOCATE(h0(NCHI,Nchi,NKPTL,NSPIN))
      ALLOCATE(hrho(NCHI,Nchi,NKPTL,NSPIN))
      SIGMA=(0.D0,0.D0)
      SIGMALAUR1=(0.D0,0.D0)
      SIGMALAUR2=(0.D0,0.D0)
      SIGMALAUR3=(0.D0,0.D0)
      sigmadc=(0.D0,0.D0)
      dedrho=(0.D0,0.D0)
      hrho=(0.D0,0.D0)
      h0=(0.D0,0.D0)
      DSIG0=(0.D0,0.D0)
      DSIGLAUR0=(0.D0,0.D0)
      GAMMA0=(0.D0,0.D0)
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
      USE DMFT_MODULE, ONLY: TON,NB,NCHI,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                      ,WKPTL,h0,sigma,GLOC,gloclaur,DENMAT,hamilton &
     &                      ,GLOCLAUR1,GLOCLAUR2,GLOCLAUR3,DIAGSLOC
      USE MPE_MODULE
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8),ALLOCATABLE :: GREENINV(:,:)
      COMPLEX(8),ALLOCATABLE :: GREEN(:,:) !(NB,NB)
      COMPLEX(8),ALLOCATABLE :: HPS0(:,:)  !(NB,NB)
      COMPLEX(8),ALLOCATABLE :: L3MAT(:,:)  !(NB,NB)
      REAL(8)   ,ALLOCATABLE :: FLOC(:) 
      COMPLEX(8)             :: CSVAR
      INTEGER(4)             :: NTASKS_K,THISTASK_K
      INTEGER(4)             :: NTASKS_M,THISTASK_M
      INTEGER(4)             :: ISPIN,IMU,NU,I,IB,IBH,IKPT,ICHI,IPRO
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT$GREEN')
      CALL DMFT_INI()
      CALL MPE$QUERY('K',NTASKS_K,THISTASK_K)
      CALL MPE$QUERY('MONOMER',NTASKS_M,THISTASK_M)
!
!     ==========================================================================
!     ==  COLLECT DFT HAMILTONIAN  (hamiltonian initialized to erho)          ==
!     ==========================================================================
      CALL DMFT$COLLECTHAMILTONIAN()
!
!     ==========================================================================
!     ==  TRANSFORMATION ONTO A ORTHONORMAL CORRELATED BASIS SET              ==
!     ==    |CHIORTHO>   =|CHINONORTHO>*DIAGSLOC
!     ==========================================================================
      call DMFT_diagsloc()
!
!     ==========================================================================
!     ==  DETERMINE CHEMICAL POTENTIAL                                        ==
!     ==  just for testing. mu should be zero                                 ==
!     ==========================================================================
      CALL DMFT$CHEMPOT()
!
!     ==========================================================================
!     ==  CALCULATE LOCAL GREENS FUNCTION AND THE ONE-PARTICLE DENSITY MATRIX ==
!     ==========================================================================
      CALL DMFT_GREENANDDENMAT()
!
!     ==========================================================================
!     ==  WRITE LOCAL FILE FOR SOLVER                                         ==
!     ==========================================================================
      CALL DMFT_WRITEGLOC_old(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC,OMEGA &
     &                   ,GLOC,GLOCLAUR1,GLOCLAUR2,GLOCLAUR3)
!CALL DMFT_TESTGLOC()
!
!     ==========================================================================
!     ==  PLOT DN/DMU AS TEST OF THE SPECTRAL PROPERTIES                      ==
!     ==========================================================================
      CALL DMFT$DOS(NOMEGA,OMEGA,NCHI,NSPIN,KBT,MU,GLOC)
!STOP 'FORCED STOP IN DMFT$GREEN'
!
!     ==========================================================================
!     ==  CHECK IF THINGS MAKE SENSE                                          ==
!     ==========================================================================
      CALL DMFT$TEST1(NOMEGA,OMEGA,NB,NKPTL,NCHI,NSPIN,KBT,WKPTL,GLOC,DENMAT &
     &               ,GLOCLAUR1,GLOCLAUR2,GLOCLAUR3)
!STOP 'STOPPING AFTER TEST1'
!
!     ==========================================================================
!     == determine local greens function                                      ==
!     ==========================================================================
      mu=0.d0
do i=1,3
      CALL DMFT_GLOC(H0,SIGma,GLOC,GLOCLAUR)
!
!     ==========================================================================
!     ==  WRITE LOCAL FILE FOR SOLVER                                         ==
!     ==========================================================================
      CALL DMFT_WRITEGLOC(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC,OMEGA &
     &                   ,GLOC,GLOCLAUR)
!
!     ==========================================================================
!     ==  CALL THE SOLVER                                                     ==
!     ==========================================================================
      CALL DMFT$SOLVER()
!
!     ==========================================================================
!     ==  read self energy from file
!     ==========================================================================
      call DMFT_GETSIGMA()
!
!     ==========================================================================
!     ==  constraints
!     ==========================================================================
      call DMFT_constraints()
enddo

stop 'end of loop. stopping.'
                                       CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_DIAGSLOC()
!     **************************************************************************
!     **  DETERMINE OVERLAP MATRIX smat(nchi,nchi,ispin) OF LOCAL ORBITALS,   **
!     **  the k-dependent inverse sinv(nchi,nchi,ikpt,ispin)                  **
!     **  AND THE MATRIX THAT CONVERTS THE LOCAL ORBITALS                     **
!     **  INTO AN ORTHONORMAL BASISSET OF CORRELATED ORBITALS                 **
!     **                                                                      **
!     **   |PI_NEU>  = |PI_ALT>  DIAGSLOC                                     **
!     **   |PI_ALT>  = |PI_NEU>  INVERT(DIAGSLOC)                             **
!     **   |CHI_ALT> = |CHI_NEU> TRANSPOSE(DIAGSLOC)                          **
!     **                                                                      **
!     **  WHERE THE PI ARE PROJECTOR FUNCTIONS AND CHI ARE ORBITALS           **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NSPIN,WKPTL,PIPSI,diagsloc,smat,sinv
      implicit none
      integer(4)      :: ikpt,ispin,i
      complex(8)      :: mat(nchi,nchi)
      real(8)         :: floc(nchi)
      logical(4)      :: ttest=.false.
!     **************************************************************************
      do ispin=1,nspin
!       ========================================================================
!       == inverse overlap matrix from sum_n <pi|psi><psi|pi>                 ==
!       ========================================================================
        mat=(0.d0,0.d0)
        DO IKPT=1,NKPTL
          sinv(:,:,ikpt,ispin)=MATMUL(PIPSI(:,:,IKPT,ISPIN) &
     &                                ,CONJG(TRANSPOSE(PIPSI(:,:,IKPT,ISPIN))))
          mat=mat+WKPTL(IKPT)*sinv(:,:,ikpt,ispin)
        enddo
!
!       ========================================================================
!       == overlap matrix                                                     **
!       ========================================================================
        call lib$invertc8(nchi,mat,smat(:,:,ispin))
!
!       ========================================================================
!       == obtain matrix that make the overlap equal to the unit matrix       ==
!       ========================================================================
        CALL LIB$DIAGC8(NCHI,mat,FLOC,DIAGSLOC(:,:,ISPIN))
        DO I=1,NCHI
          DIAGSLOC(:,I,ISPIN)=DIAGSLOC(:,I,ISPIN)/SQRT(FLOC(I))
        ENDDO
!
!       ========================================================================
!       == OPTIONAL TESTS
!       ========================================================================
        IF(TTEST) THEN
          MAT=MATMUL(CONJG(TRANSPOSE(DIAGSLOC(:,:,ISPIN))) &
                                ,MATMUL(mat,DIAGSLOC(:,:,ISPIN)))
          DO I=1,NCHI
            WRITE(*,FMT='("TEST",10("(",2F10.5,")"))')MAT(I,:)
          ENDDO
          CALL ERROR$MSG('SCHEDULED STOP AFTER TEST')
          CALL ERROR$STOP('DMFT_DIAGSLOC')
        END IF
      ENDDO
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$COLLECTHAMILTONIAN()
!     **************************************************************************
!     ** COLLECTS THE HAMILTONIAN AND STORES IT ON THE MODULE                 **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,IPROOFCHI,kbt &
     &                       ,erho,hamilton,PIPSI
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: NGL
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IPRO,ib
      real(8)                :: f(nb,nkptl,nspin)
      real(8)                :: svar
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                      CALL TRACE$PUSH('DMFT$COLLECTHAMILTONIAN')
!
!     ==========================================================================
!     ==  
!     ==========================================================================
      call DMFT$COLLECTOCCUPATIONS(NKPTL,NSPIN,NB,F)
      hamilton=(0.d0,0.d0)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          do ib=1,nb
            call dmft_eoff(kbt,f(ib,ikpt,ispin),erho(ib,ikpt,ispin))
            hamilton(ib,ib,ikpt,ispin)=cmplx(erho(ib,ikpt,ispin))
          enddo
        enddo
      enddo 
print*,'bounds of spectrum [ev] ',minval(erho)*27.211d0,maxval(erho)*27.211d0
!
!     ==========================================================================
!     ==  
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
          DO ICHI=1,NCHI
            IPRO=IPROOFCHI(ICHI)
            DO IBH=1,NBH
              IF(NBH.NE.NB) THEN
                PIPSI(ICHI,2*IBH-1,IKPT,ISPIN)= REAL(THIS%TBC(1,IBH,IPRO))
                PIPSI(ICHI,2*IBH  ,IKPT,ISPIN)=AIMAG(THIS%TBC(1,IBH,IPRO))
              ELSE
                PIPSI(ICHI,IBH,IKPT,ISPIN)  =THIS%TBC(1,IBH,IPRO)
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
      REAL(8)                :: WKPT(NKPTL,NSPIN)
      integer(4)             :: ikpt,ispin
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
      DO IKPT=1,NKPTl
        DO ISPIN=1,NSPIN
          F(:,IKPT,ISPIN)=F(:,IKPT,ISPIN)/WKPT(IKPT,ISPIN)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE DMFT_EOFF(KBT,F,E)
!     **************************************************************************
!     ** CALCULATES ENERGY FROM F=1/(1+EXP(-E/KBT))                           **
!     ** SMALL AND LARGE OCCUPATIONS ARE SIMPLIFIED BY LINEAR EXTRAPOLATION   **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)  :: KBT
      REAL(8),INTENT(IN)  :: F
      REAL(8),INTENT(OUT) :: E
      REAL(8),PARAMETER   :: FMIN=1.D-6
      real(8)             :: f0,slope,val
!     **************************************************************************
      IF(F.GT.FMIN.AND.F.lT.1.D0-FMIN) THEN
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
      e=e*kbt
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$CHEMPOT()
!     **************************************************************************
!     ** CALCULATES THE LOCAL INTERACTING GREENS FUNCTION                     **
!     **                                                                      **
!     ** PIPSI <PI_A|PSI_N> IS THE PRE-FACTOR OF LOCAL ORBITAL |CHI_A> IN     **
!     **       THE LOCAL ORBITAL EXPANSION OF |PSI_N>                         **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NB,NCHI,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                      ,WKPTL,GLOC,DSIG0,DENMAT &
     &                      ,GLOCLAUR1,GLOCLAUR2,GLOCLAUR3,DIAGSLOC,DSIGLAUR0 &
     &                      ,HAMILTON,PIPSI
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:) ! OCCUPATIONS (INCLUDING K-WEIGHT)
      REAL(8)                :: NEL        ! #(ELECTRONS) 
      REAL(8)                :: NOFMU      ! PARTICLE NUMBER FOR MUGRID
      REAL(8)                :: NOFMU1     ! NOFMU CONTRIB FROM 1 KPOINT
      COMPLEX(8),ALLOCATABLE :: HPS0(:,:)  !(NB,NB)
      COMPLEX(8),ALLOCATABLE :: GREENINV(:,:)
      COMPLEX(8),ALLOCATABLE :: GREEN(:,:) !(NB,NB)
      COMPLEX(8)             :: TR1,TRDS1,TRHPS0,TRHPS0SQUARE,LAURENT1,LAURENT2
      INTEGER(4)             :: NTASKS_K,THISTASK_K
      INTEGER(4)             :: NTASKS_M,THISTASK_M
      COMPLEX(8)             :: CSVAR
      REAL(8)                :: SVAR
      INTEGER(4)             :: IMU,IKPT,ISPIN,NU,IB
      INTEGER(4)             :: ISTART
      INTEGER(4)             :: IBI
      REAL(8)                :: X0,XM,DX,Y0,YM
      REAL(8)                :: EV
      REAL(8)   ,PARAMETER   :: TOL=1.D-5
!     **************************************************************************
      IF(.NOT.TON) RETURN
PRINT*,'DETERMINING CHEMPOT....'
      CALL CONSTANTS$GET('EV',EV)
      CALL MPE$QUERY('K',NTASKS_K,THISTASK_K)
      CALL MPE$QUERY('MONOMER',NTASKS_M,THISTASK_M)
!
!     ==========================================================================
!     ==  DETERMINE TARGET NUMBER OF ELECTRONS                                ==
!     ==========================================================================
      ALLOCATE(OCC(NB,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NB*NKPTL*NSPIN,OCC)
      NEL=SUM(OCC)
      DEALLOCATE(OCC)
PRINT*,'DMFT$CHEMPOT: #(ELECTRONS)=',NEL,' #(BANDS)=',NB
!
!     ==========================================================================
!     ==  ALLOCATE ARRAYS                                                     ==
!     ==========================================================================
      ALLOCATE(GREENINV(NB,NB))
      ALLOCATE(GREEN(NB,NB))
      ALLOCATE(HPS0(NB,NB))
!
!     ==========================================================================
!     ==  DETERMINE CHEMICAL POTENTIAL
!     ==========================================================================
      X0=MU
      DX=1.D-2
      ISTART=1
      CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
1000  CONTINUE
      NOFMU=0.D0
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
!
!         ======================================================================
!         ==  PREPARE TERMS REQUIRED FOR THE LAURENT EXPANSION                ==
!         ======================================================================
          TR1=REAL(NB,KIND=8)  ! tr1=trace of 1
!         == DO THE HIGHER ORDER TERM FIRST TO BE ABLE TO USE HPS0 AS WORK ARRAY
!         == hps0=<psi|pi>S(1)<pi|psi> =========================================
          HPS0=MATMUL(TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN))) &
      &              ,MATMUL(DSIGLAUR0(:,:,ISPIN,2),PIPSI(:,:,IKPT,ISPIN)))
          TRDS1=(0.D0,0.D0) !trds=trace of ds(1)
          DO IB=1,NB
            TRDS1=TRDS1+HPS0(IB,IB)
          ENDDO
!         == HPS0= HAMILTON+<psi|pi>S(0)<pi|psi> ===============================
          HPS0=HAMILTON(:,:,IKPT,ISPIN) &
      &            +MATMUL(TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN))) &
      &                   ,MATMUL(DSIGLAUR0(:,:,ISPIN,1),PIPSI(:,:,IKPT,ISPIN)))
          TRHPS0=(0.D0,0.D0) ! trhps0=trace of h plus s0
          DO IB=1,NB
            TRHPS0=TRHPS0+HPS0(IB,IB)
          ENDDO
!
          TRHPS0SQUARE=SUM(HPS0**2)
!
!         ======================================================================
!         == LAURENT EXPANSION COEFFICIENTS REQUIRED FOR THE REGULARIZATION   ==
!         ======================================================================
          LAURENT1=TRHPS0-MU*TR1
          LAURENT2=(TRHPS0SQUARE-2.D0*TRHPS0*MU+MU**2*TR1)+TRDS1
!
!         ======================================================================
!         ==  PERFORM MATSUBARA SUMS                                          ==
!         ======================================================================
          NOFMU1=0.D0
          DO NU=1,NOMEGA
            IF(MODULO(NU,NTASKS_K).NE.0) CYCLE
!           == CONSTRUCT LATTICE GREENS FUNCTION ===============================
            GREENINV(:,:)=-HAMILTON(:,:,IKPT,ISPIN) &
     &            -MATMUL(TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN))) &
     &                   ,MATMUL(DSIG0(:,:,NU,ISPIN),PIPSI(:,:,IKPT,ISPIN)))
            CSVAR=CI*OMEGA(NU)+MU
            DO IB=1,NB
              GREENINV(IB,IB)=GREENINV(IB,IB)+CSVAR
            ENDDO
            CALL LIB$INVERTC8(NB,GREENINV,GREEN)
            DO IB=1,NB
              NOFMU1=NOFMU1+WKPTL(IKPT)*KBT*REAL(GREEN(IB,IB))
            ENDDO
!
!           == REGULARIZE THE GREENS FUNCTION BEFORE ADDING IT UP ==============
            CSVAR=TR1/(CI*OMEGA(NU))+LAURENT1/(CI*OMEGA(NU))**2 &
        &                           +LAURENT2/(CI*OMEGA(NU))**3
            NOFMU1=NOFMU1-WKPTL(IKPT)*KBT*REAL(CSVAR)
          ENDDO
          NOFMU=NOFMU+2.D0*NOFMU1   ! ADD NEGATIVE FREQUENCIES
!         == ADD SUM OF LONG-RANGE TAILS =======================================
          SVAR=REAL(0.5D0*TR1-0.25D0/KBT*LAURENT1)
          NOFMU=NOFMU+WKPTL(IKPT)*SVAR/REAL(NTASKS_K)
        ENDDO
      ENDDO
      CALL MPE$COMBINE('MONOMER','+',NOFMU)
!     == TAKE CARE OF SPIN MULTIPLICITY ========================================
      IF(NSPIN.EQ.1.AND.NDIM.EQ.1) NOFMU=NOFMU*2.D0
!
      Y0=NOFMU-NEL
      CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
      MU=X0
WRITE(*,FMT='("MU[EV]=",F15.8," N=",F15.8)')X0/EV,NOFMU
      IF(ABS(Y0).GT.TOL) GOTO 1000
!
PRINT*,'NEL  ',NEL,' KBT=',KBT
PRINT*,'KBT= ',KBT,' KBT[EV]=',KBT*27.211D0
PRINT*,'MU = ',MU ,' MU[EV] =',MU*27.211D0
PRINT*,'..... CHEMPOT DETERMINED'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_GREENANDDENMAT()
!     **************************************************************************
!     ** CALCULATES THE LOCAL INTERACTING GREENS FUNCTION                     **
!     ** AND THE DENSITY MATRIX                                               **
!     **                                                                      **
!     ** PIPSI <PI_A|PSI_N> IS THE PRE-FACTOR OF LOCAL ORBITAL |CHI_A> IN     **
!     **       THE LOCAL ORBITAL EXPANSION OF |PSI_N>                         **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NB,NCHI,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                      ,WKPTL,GLOC,DSIG0,DENMAT &
     &                      ,GLOCLAUR1,GLOCLAUR2,GLOCLAUR3,DSIGLAUR0 &
     &                      ,HAMILTON,PIPSI
      USE MPE_MODULE
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8),ALLOCATABLE :: GREENINV(:,:)
      COMPLEX(8),ALLOCATABLE :: GREEN(:,:) !(NB,NB)
      COMPLEX(8),ALLOCATABLE :: HPS0(:,:)  !(NB,NB)
      COMPLEX(8),ALLOCATABLE :: L3MAT(:,:)  !(NB,NB)
      COMPLEX(8),ALLOCATABLE :: DENMAT1(:,:) !(NB,NB)
      REAL(8)   ,ALLOCATABLE :: FLOC(:) 
      COMPLEX(8)             :: CSVAR
      INTEGER(4)             :: NTASKS_K,THISTASK_K
      INTEGER(4)             :: NTASKS_M,THISTASK_M
      INTEGER(4)             :: ISPIN,IMU,NU,I,IB,IBH,IKPT,ICHI,IPRO
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT_GREENANDDENMAT')
      CALL MPE$QUERY('K',NTASKS_K,THISTASK_K)
      CALL MPE$QUERY('MONOMER',NTASKS_M,THISTASK_M)
!
!     ==========================================================================
!     ==  ALLOCATE ARRAYS                                                     ==
!     ==========================================================================
      ALLOCATE(GREENINV(NB,NB))
      ALLOCATE(GREEN(NB,NB))
      ALLOCATE(HPS0(NB,NB))
      ALLOCATE(L3MAT(NB,NB)) ! THIRD TERM OF LAURENT EXPANSION
      ALLOCATE(DENMAT1(NB,NB)) 
!
!     ==========================================================================
!     ==  ACCUMULATE LOCAL GREENS FUNCTION AND ONE-PARTICLE DENSITY MATRIX    ==
!     ==========================================================================
      GLOC=(0.D0,0.D0)
      DENMAT=(0.D0,0.D0)
      GLOCLAUR1=(0.D0,0.D0)
      GLOCLAUR2=(0.D0,0.D0)
      GLOCLAUR3=(0.D0,0.D0)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
!
!         ======================================================================
!         ==  CALCULATE MATRICES FOR THE REGULARIZATION                       ==
!         ==  GLAT -> 1/IOMEGA + HPS0/IOMEGA^2 + L3MAT/IOMEGA^3               ==
!         ==  GLOC -> GLOCLAUR1/IOMEGA + GLOCLAUR2/IOMEGA^2+GLOCLAUR3/IOMEGA^3==
!         ======================================================================
!         == HPS0=HAMILTON+SIGMA-MU ============================================
          HPS0=HAMILTON(:,:,IKPT,ISPIN) &
       &         +MATMUL(TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN))) &
       &                  ,MATMUL(DSIGLAUR0(:,:,ISPIN,1),PIPSI(:,:,IKPT,ISPIN)))
          DO IB=1,NB
            HPS0(IB,IB)=HPS0(IB,IB)-MU
          ENDDO
!         == L3MAT=HPS0**2 =====================================================
          L3MAT=MATMUL(HPS0,HPS0) &
       &       +MATMUL(TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN))) &
       &                ,MATMUL(DSIGLAUR0(:,:,ISPIN,2),PIPSI(:,:,IKPT,ISPIN)))
!
!         == LAURENT EXPANSION TERMS FOR THE LOCAL GREENS FUNCTION =============
          GLOCLAUR1(:,:,ISPIN)=GLOCLAUR1(:,:,ISPIN)+WKPTL(IKPT) &
       &                               *MATMUL(PIPSI(:,:,IKPT,ISPIN) &
       &                               ,CONJG(TRANSPOSE(PIPSI(:,:,IKPT,ISPIN))))
          GLOCLAUR2(:,:,ISPIN)=GLOCLAUR2(:,:,ISPIN)+WKPTL(IKPT) &
       &                  *MATMUL(PIPSI(:,:,IKPT,ISPIN) &
       &                  ,MATMUL(HPS0,CONJG(TRANSPOSE(PIPSI(:,:,IKPT,ISPIN)))))
          GLOCLAUR3(:,:,ISPIN)=GLOCLAUR3(:,:,ISPIN)+WKPTL(IKPT) &
       &                 *MATMUL(PIPSI(:,:,IKPT,ISPIN) &
       &                 ,MATMUL(L3MAT,CONJG(TRANSPOSE(PIPSI(:,:,IKPT,ISPIN)))))
!
!         ======================================================================
!         ==  PERFORM MATSUBARA SUMS                                          ==
!         ======================================================================
          DENMAT1(:,:)=(0.D0,0.D0)
          DO NU=1,NOMEGA
            IF(MODULO(NU,NTASKS_K).NE.0) CYCLE
!
!           == INVERSE LATTICE GREENS FUNCTION =================================
            GREENINV(:,:)=-HAMILTON(:,:,IKPT,ISPIN) &
     &             -MATMUL(TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN))) &
     &                    ,MATMUL(DSIG0(:,:,NU,ISPIN),PIPSI(:,:,IKPT,ISPIN)))
            CSVAR=CI*OMEGA(NU)+MU
            DO IB=1,NB
              GREENINV(IB,IB)=GREENINV(IB,IB)+CSVAR
            ENDDO
!
!           == INVERT TO OBTAIN LATTICE GREENS FUNCTION ========================
            CALL LIB$INVERTC8(NB,GREENINV,GREEN)
!
!           == CONSTRUCT INTERACTING LOCAL GREENS FUNCTION =====================
            GLOC(:,:,NU,ISPIN)=GLOC(:,:,NU,ISPIN)+WKPTL(IKPT) &
     &                 *MATMUL(PIPSI(:,:,IKPT,ISPIN) &
     &                   ,MATMUL(GREEN,CONJG(TRANSPOSE(PIPSI(:,:,IKPT,ISPIN)))))
!
!           == REGULARIZE BEFORE ADDING UP 1P DENSITY MATRIX====================
            CSVAR=1.D0/(CI*OMEGA(NU))
            DO IB=1,NB
              GREEN(IB,IB)=GREEN(IB,IB)-CSVAR
            ENDDO
            GREEN=GREEN-CSVAR**2*HPS0-CSVAR**3*L3MAT
!
!           == ACCUMULATE DENSITY MATRIX =======================================
            DENMAT1(:,:)=DENMAT1(:,:)+KBT*GREEN(:,:)
          ENDDO
!         == CONTRIBUTION FROM NEGATIVE FREQUENCIES ============================
          DENMAT1=DENMAT1+TRANSPOSE(CONJG(DENMAT1))
          DENMAT(:,:,IKPT,ISPIN)=DENMAT(:,:,IKPT,ISPIN)+DENMAT1
!
!         == FINISH REGULARIZATION =============================================
          DO IB=1,NB
            IF(MODULO(IB,NTASKS_K).NE.0) CYCLE
            DENMAT(IB,IB,IKPT,ISPIN)=DENMAT(IB,IB,IKPT,ISPIN)+(0.5D0,0.D0)
          ENDDO
          DENMAT(:,:,IKPT,ISPIN)=DENMAT(:,:,IKPT,ISPIN)-0.25D0/KBT*HPS0
        ENDDO
      ENDDO
      IF(NSPIN.EQ.1.AND.NDIM.EQ.1) THEN  !SPIN MULTIPLICITY
        DENMAT=DENMAT*2.D0
      END IF
!
!     == ADD OVER ALL NODES (INCLUDES SUM IN K-GROUPS AND OVER K-GROUPS) =======
      CALL MPE$COMBINE('MONOMER','+',GLOC)
      CALL MPE$COMBINE('K','+',DENMAT)
                                       CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_WRITEGLOC(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC,OMEGA &
     &                          ,GLOC,GLOCLAUR)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NOMEGA
      INTEGER(4),INTENT(IN)  :: NSPIN
      REAL(8)   ,INTENT(IN)  :: KBT
      REAL(8)   ,INTENT(IN)  :: MU
      REAL(8)   ,INTENT(IN)  :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN)  :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(IN)  :: DIAGSLOC(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)             :: GLOCPRIME(NCHI,NCHI)
      COMPLEX(8)             :: T(NCHI,NCHI)
      COMPLEX(8)             :: TPLUS(NCHI,NCHI)
      INTEGER(4)             :: NU,ISPIN,I,J
      INTEGER(4)             :: NFIL=11
!     **************************************************************************
      OPEN(NFIL,FILE=-'GLOC.DATA')
      WRITE(NFIL,*)NCHI,NOMEGA,NSPIN,KBT,MU !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DO ISPIN=1,NSPIN
        T=DIAGSLOC(:,:,ISPIN)
        TPLUS=TRANSPOSE(CONJG(T)) 
        DO NU=1,NOMEGA
          GLOCPRIME=MATMUL(TPLUS,MATMUL(GLOC(:,:,NU,ISPIN),T))
          WRITE(NFIL,*)((REAL(GLOCPRIME(I,J)),AIMAG(GLOCPRIME(I,J)) &
     &                ,I=1,NCHI),J=1,NCHI),OMEGA(NU)
        ENDDO
        GLOCPRIME=MATMUL(TPLUS,MATMUL(GLOCLAUR(:,:,1,ISPIN),T))
        WRITE(NFIL,*)((REAL(GLOCPRIME(I,J)),AIMAG(GLOCPRIME(I,J)) &
     &                ,I=1,NCHI),J=1,NCHI)
        GLOCPRIME=MATMUL(TPLUS,MATMUL(GLOCLAUR(:,:,2,ISPIN),T))
        WRITE(NFIL,*)((REAL(GLOCPRIME(I,J)),AIMAG(GLOCPRIME(I,J)) &
     &                ,I=1,NCHI),J=1,NCHI)
        GLOCPRIME=MATMUL(TPLUS,MATMUL(GLOCLAUR(:,:,3,ISPIN),T))
        WRITE(NFIL,*)((REAL(GLOCPRIME(I,J)),AIMAG(GLOCPRIME(I,J)) &
     &                ,I=1,NCHI),J=1,NCHI)
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_WRITEGLOC_old(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC,OMEGA &
     &                          ,GLOC,GLOCLAUR1,GLOCLAUR2,GLOCLAUR3)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NOMEGA
      INTEGER(4),INTENT(IN)  :: NSPIN
      REAL(8)   ,INTENT(IN)  :: KBT
      REAL(8)   ,INTENT(IN)  :: MU
      REAL(8)   ,INTENT(IN)  :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN)  :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(IN)  :: DIAGSLOC(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR1(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR2(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR3(NCHI,NCHI,NSPIN)
      COMPLEX(8)             :: GLOCPRIME(NCHI,NCHI)
      COMPLEX(8)             :: T(NCHI,NCHI)
      COMPLEX(8)             :: TPLUS(NCHI,NCHI)
      INTEGER(4)             :: NU,ISPIN,I,J
      INTEGER(4)             :: NFIL=11
!     **************************************************************************
      OPEN(NFIL,FILE=-'GLOC.DATA')
      WRITE(NFIL,*)NCHI,NOMEGA,NSPIN,KBT,MU !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DO ISPIN=1,NSPIN
        T=DIAGSLOC(:,:,ISPIN)
        TPLUS=TRANSPOSE(CONJG(T)) 
        DO NU=1,NOMEGA
          GLOCPRIME=MATMUL(TPLUS,MATMUL(GLOC(:,:,NU,ISPIN),T))
          WRITE(NFIL,*)((REAL(GLOCPRIME(I,J)),AIMAG(GLOCPRIME(I,J)) &
     &                ,I=1,NCHI),J=1,NCHI),OMEGA(NU)
        ENDDO
        GLOCPRIME=MATMUL(TPLUS,MATMUL(GLOCLAUR1(:,:,ISPIN),T))
        WRITE(NFIL,*)((REAL(GLOCPRIME(I,J)),AIMAG(GLOCPRIME(I,J)) &
     &                ,I=1,NCHI),J=1,NCHI)
        GLOCPRIME=MATMUL(TPLUS,MATMUL(GLOCLAUR2(:,:,ISPIN),T))
        WRITE(NFIL,*)((REAL(GLOCPRIME(I,J)),AIMAG(GLOCPRIME(I,J)) &
     &                ,I=1,NCHI),J=1,NCHI)
        GLOCPRIME=MATMUL(TPLUS,MATMUL(GLOCLAUR3(:,:,ISPIN),T))
        WRITE(NFIL,*)((REAL(GLOCPRIME(I,J)),AIMAG(GLOCPRIME(I,J)) &
     &                ,I=1,NCHI),J=1,NCHI)
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_TESTGLOC()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NSPIN,NOMEGA
      USE STRINGS_MODULE
      IMPLICIT NONE
      REAL(8)    :: KBT
      REAL(8)    :: MU
      COMPLEX(8) :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8) :: GLOCLAUR1(NCHI,NCHI)
      COMPLEX(8) :: GLOCLAUR2(NCHI,NCHI)
      COMPLEX(8) :: GLOCLAUR3(NCHI,NCHI)
      COMPLEX(8) :: GLOC1(NCHI,NCHI)
      REAL(8)    :: OMEGA(NOMEGA)      
      INTEGER(4) :: NU
      INTEGER(4) :: NFIL1=11,NFIL2=12,NFIL3=13,NFIL4=14
      COMPLEX(8),PARAMETER :: CI=(0.D0,1.D0)
!     **************************************************************************
      CALL DMFT_READGLOC(NCHI,NOMEGA,NSPIN,KBT,MU,OMEGA &
     &                  ,GLOC,GLOCLAUR1,GLOCLAUR2,GLOCLAUR3)
      OPEN(NFIL1,FILE=-'READGLOC_RE.DAT')
      OPEN(NFIL2,FILE=-'READGLOC_IM.DAT')
      OPEN(NFIL3,FILE=-'READGLOC_WOT_RE.DAT')
      OPEN(NFIL4,FILE=-'READGLOC_WOT_IM.DAT')
      DO NU=1,NOMEGA
        GLOC1=GLOC(:,:,NU,1)
        WRITE(NFIL1,*)OMEGA(NU),REAL(GLOC1)
        WRITE(NFIL2,*)OMEGA(NU),AIMAG(GLOC1)
        GLOC1=GLOC1-(GLOCLAUR1/(CI*OMEGA(NU)) &
     &              +GLOCLAUR2/(CI*OMEGA(NU))**2+GLOCLAUR3/(CI*OMEGA(NU))**3)
        WRITE(NFIL3,*)OMEGA(NU),REAL(GLOC1)
        WRITE(NFIL4,*)OMEGA(NU),AIMAG(GLOC1)
      ENDDO
      CLOSE(NFIL1)
      CLOSE(NFIL2)
      CLOSE(NFIL3)
      CLOSE(NFIL4)
STOP 'FORCED IN DMFT_TESTGLOC'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_READGLOC(NCHI,NOMEGA,NSPIN,KBT,MU,OMEGA &
     &                        ,GLOC,GLOCLAUR1,GLOCLAUR2,GLOCLAUR3)
!     **************************************************************************
!     ** THIS ROUTINE IS FOR TESTING WHAT HAS BEEN WRITTEN TO THE FILE        **
!     ** IT DOES NOT UNDO THE TRANSFORMATION TO ORTHONORMAL LOCAL ORBITALS    **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NOMEGA
      INTEGER(4),INTENT(IN)  :: NSPIN
      REAL(8)   ,INTENT(OUT) :: KBT
      REAL(8)   ,INTENT(OUT) :: MU
      REAL(8)   ,INTENT(OUT) :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(OUT) :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(OUT) :: GLOCLAUR1(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(OUT) :: GLOCLAUR2(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(OUT) :: GLOCLAUR3(NCHI,NCHI,NSPIN)
      REAL(8)                :: GLOC_RE(NCHI,NCHI),GLOC_IM(NCHI,NCHI)
      INTEGER(4)             :: NCHI_,NOMEGA_,NSPIN_
      INTEGER(4)             :: NU,ISPIN,I,J
      INTEGER(4)             :: NFIL=11
!     **************************************************************************
      OPEN(NFIL,FILE=-'GLOC.DATA')
      READ(NFIL,*)NCHI_,NOMEGA_,NSPIN_,KBT,MU !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF(NCHI_.NE.NCHI) THEN
        CALL ERROR$MSG('NCHI ON FILE INCONSISTENT WITH INTERNAL VALUE')
        CALL ERROR$CHVAL('NCHI ON FILE',NCHI_)
        CALL ERROR$CHVAL('INTERNAL NCHI',NCHI)
        CALL ERROR$STOP('DMFT_READGLOC')
      END IF
      IF(NOMEGA_.NE.NOMEGA) THEN
        CALL ERROR$MSG('NOMEGA ON FILE INCONSISTENT WITH INTERNAL VALUE')
        CALL ERROR$CHVAL('NOMEGA ON FILE',NOMEGA_)
        CALL ERROR$CHVAL('INTERNAL NOMEGA',NOMEGA)
        CALL ERROR$STOP('DMFT_READGLOC')
      END IF
      IF(NSPIN_.NE.NSPIN) THEN
        CALL ERROR$MSG('NSPIN ON FILE INCONSISTENT WITH INTERNAL VALUE')
        CALL ERROR$CHVAL('NSPIN ON FILE',NSPIN_)
        CALL ERROR$CHVAL('INTERNAL NSPIN',NSPIN)
        CALL ERROR$STOP('DMFT_READGLOC')
      END IF
      DO ISPIN=1,NSPIN
        DO NU=1,NOMEGA
          READ(NFIL,*)((GLOC_RE(I,J),GLOC_IM(I,J),I=1,NCHI),J=1,NCHI),OMEGA(NU)
          GLOC(:,:,NU,ISPIN)=CMPLX(GLOC_RE,GLOC_IM)
        ENDDO
        READ(NFIL,*)((GLOC_RE(I,J),GLOC_IM(I,J),I=1,NCHI),J=1,NCHI)
        GLOCLAUR1(:,:,ISPIN)=CMPLX(GLOC_RE,GLOC_IM)
        READ(NFIL,*)((GLOC_RE(I,J),GLOC_IM(I,J),I=1,NCHI),J=1,NCHI)
        GLOCLAUR2(:,:,ISPIN)=CMPLX(GLOC_RE,GLOC_IM)
        READ(NFIL,*)((GLOC_RE(I,J),GLOC_IM(I,J),I=1,NCHI),J=1,NCHI)
        GLOCLAUR3(:,:,ISPIN)=CMPLX(GLOC_RE,GLOC_IM)
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_READSIGMA(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC &
     &                          ,DEDGLOC,DEDGLOCLAUR1,DEDGLOCLAUR2 &
     &                          ,DEDGLOCLAUR3,SDC)
!     **************************************************************************
!     ** READ DERIVATIVE OF THE LUTTINGER WARD FUNCTIONAL WITH RESPECT TO     **
!     ** THE LOCAL GREENS FUNCTION AND ITS LAURENT EXPANSION TERMS            **
!     **                                                                      **
!     **  A PROPER DEFINITION OF THE DERIVATIVE OF THE LAURENT EXPANION TERMS **
!     **  IS MISSING                                                          **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NOMEGA
      INTEGER(4),INTENT(IN)  :: NSPIN
      REAL(8)   ,INTENT(IN)  :: KBT
      REAL(8)   ,INTENT(IN)  :: MU
      COMPLEX(8),INTENT(IN)  :: DIAGSLOC(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(OUT) :: DEDGLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(OUT) :: DEDGLOCLAUR1(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(OUT) :: DEDGLOCLAUR2(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(OUT) :: DEDGLOCLAUR3(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(OUT) :: SDC(NCHI,NCHI,NSPIN) ! DOUBLE COUNTING
      COMPLEX(8)             :: SIGMAPRIME(NCHI,NCHI)
      REAL(8)                :: S_RE(NCHI,NCHI)
      REAL(8)                :: S_IM(NCHI,NCHI)
      COMPLEX(8)             :: T(NCHI,NCHI)
      COMPLEX(8)             :: TPLUS(NCHI,NCHI)
      INTEGER(4)             :: NCHI_
      INTEGER(4)             :: NOMEGA_
      INTEGER(4)             :: NSPIN_
      REAL(8)                :: KBT_
      REAL(8)                :: MU_
      INTEGER(4)             :: NU,ISPIN,I,J
      INTEGER(4)             :: NFIL=11
!     **************************************************************************
      OPEN(NFIL,FILE=-'SIGMA.DATA')
      READ(NFIL,*)NCHI_,NOMEGA_,NSPIN_,KBT_,MU_ !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF(NCHI_.NE.NCHI) THEN
        CALL ERROR$MSG('NCHI INCONSISTENT')
        CALL ERROR$STOP('DMFT_WRITEDEGLOC')
      END IF
      IF(NOMEGA_.NE.NOMEGA) THEN
        CALL ERROR$MSG('NOMEGA INCONSISTENT')
        CALL ERROR$STOP('DMFT_WRITEDEGLOC')
      END IF
      IF(NSPIN_.NE.NSPIN) THEN
        CALL ERROR$MSG('NOMEGA INCONSISTENT')
        CALL ERROR$STOP('DMFT_WRITEDEGLOC')
      END IF
      DO ISPIN=1,NSPIN
        T=TRANSPOSE(CONJG(DIAGSLOC(:,:,ISPIN))) 
        TPLUS=TRANSPOSE(CONJG(T)) 
        DO NU=1,NOMEGA
          READ(NFIL,*)((S_RE(I,J),S_IM(I,J),I=1,NCHI),J=1,NCHI)
          SIGMAPRIME=CMPLX(S_RE,S_IM)
          DEDGLOC(:,:,NU,ISPIN)=MATMUL(T,MATMUL(SIGMAPRIME,TPLUS))
        ENDDO
        READ(NFIL,*)((S_RE(I,J),S_IM(I,J),I=1,NCHI),J=1,NCHI)
        SIGMAPRIME=CMPLX(S_RE,S_IM)
        DEDGLOCLAUR1(:,:,ISPIN)=MATMUL(T,MATMUL(SIGMAPRIME,TPLUS))
        READ(NFIL,*)((S_RE(I,J),S_IM(I,J),I=1,NCHI),J=1,NCHI)
        SIGMAPRIME=CMPLX(S_RE,S_IM)
        DEDGLOCLAUR2(:,:,ISPIN)=MATMUL(T,MATMUL(SIGMAPRIME,TPLUS))
        READ(NFIL,*)((S_RE(I,J),S_IM(I,J),I=1,NCHI),J=1,NCHI)
        SIGMAPRIME=CMPLX(S_RE,S_IM)
        DEDGLOCLAUR3(:,:,ISPIN)=MATMUL(T,MATMUL(SIGMAPRIME,TPLUS))
        READ(NFIL,*)((S_RE(I,J),S_IM(I,J),I=1,NCHI),J=1,NCHI)
        SIGMAPRIME=CMPLX(S_RE,S_IM)
        SDC(:,:,ISPIN)=MATMUL(T,MATMUL(SIGMAPRIME,TPLUS))
      ENDDO
      CLOSE(NFIL)
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_WRITESIGMA(NCHI,NOMEGA,NSPIN,KBT,MU &
     &                          ,DEDGLOC,DEDGLOCLAUR1,DEDGLOCLAUR2,DEDGLOCLAUR3,SDC)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NOMEGA
      INTEGER(4),INTENT(IN)  :: NSPIN
      REAL(8)   ,INTENT(IN)  :: KBT
      REAL(8)   ,INTENT(IN)  :: MU
      COMPLEX(8),INTENT(OUT) :: DEDGLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(OUT) :: DEDGLOCLAUR1(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(OUT) :: DEDGLOCLAUR2(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(OUT) :: DEDGLOCLAUR3(NCHI,NCHI,NSPIN)
      COMPLEX(8),INTENT(OUT) :: SDC(NCHI,NCHI,NSPIN) ! DOUBLE COUNTING
      COMPLEX(8)             :: SIGMAPRIME(NCHI,NCHI)
      INTEGER(4)             :: NU,ISPIN,I,J
      INTEGER(4)             :: NFIL=11
!     **************************************************************************
      OPEN(NFIL,FILE=-'SIGMA.DATA')
      WRITE(NFIL,*)NCHI,NOMEGA,NSPIN,KBT,MU !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DO ISPIN=1,NSPIN
        DO NU=1,NOMEGA
          SIGMAPRIME=DEDGLOC(:,:,NU,ISPIN)
          WRITE(NFIL,*)((REAL(SIGMAPRIME(I,J)),AIMAG(SIGMAPRIME(I,J)) &
     &                                                    ,I=1,NCHI),J=1,NCHI)
        ENDDO
        SIGMAPRIME=DEDGLOCLAUR1(:,:,ISPIN)
        WRITE(NFIL,*)((REAL(SIGMAPRIME(I,J)),AIMAG(SIGMAPRIME(I,J)) &
     &                                                    ,I=1,NCHI),J=1,NCHI)
        SIGMAPRIME=DEDGLOCLAUR2(:,:,ISPIN)
        WRITE(NFIL,*)((REAL(SIGMAPRIME(I,J)),AIMAG(SIGMAPRIME(I,J)) &
     &                                                    ,I=1,NCHI),J=1,NCHI)
        SIGMAPRIME=DEDGLOCLAUR3(:,:,ISPIN)
        WRITE(NFIL,*)((REAL(SIGMAPRIME(I,J)),AIMAG(SIGMAPRIME(I,J)) &
     &                                                    ,I=1,NCHI),J=1,NCHI)
        SIGMAPRIME=SDC(:,:,ISPIN)
        WRITE(NFIL,*)((REAL(SIGMAPRIME(I,J)),AIMAG(SIGMAPRIME(I,J)) &
     &                                                    ,I=1,NCHI),J=1,NCHI)
      ENDDO
      CLOSE(NFIL)
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$SOLVER()
!     **************************************************************************
!     ** mimicks a dmft solver                                                **
!     **  e=1/beta * sum_nu e(iomega_nu0+) Tr[ dv * gloc^\dagger]             **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NB,NCHI,NKPTL,NSPIN,NDIM,NOMEGA
      IMPLICIT NONE
      REAL(8)          :: KBT
      REAL(8)          :: MU
      REAL(8)          :: OMEGA(NOMEGA)
      COMPLEX(8)       :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: GLOCLAUR1(NCHI,NCHI,NSPIN)
      COMPLEX(8)       :: GLOCLAUR2(NCHI,NCHI,NSPIN)
      COMPLEX(8)       :: GLOCLAUR3(NCHI,NCHI,NSPIN)
      COMPLEX(8)       :: DH(NCHI,NCHI,NSPIN)
      COMPLEX(8)       :: DEDGLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: DEDGLOCLAUR1(NCHI,NCHI,NSPIN)
      COMPLEX(8)       :: DEDGLOCLAUR2(NCHI,NCHI,NSPIN)
      COMPLEX(8)       :: DEDGLOCLAUR3(NCHI,NCHI,NSPIN)
      COMPLEX(8)       :: SDC(NCHI,NCHI,NSPIN)
      REAL(8)          :: EV
      INTEGER(4)       :: ISPIN,NU
!     **************************************************************************
PRINT*,'ENTERING DMFT$SOLVER...'
      CALL CONSTANTS('EV',EV)
      DH=(0.D0,0.D0)
      DO ISPIN=1,NSPIN
        DH(1,1,ISPIN)=-1.D0*EV
        DH(2,2,ISPIN)=-1.D0*EV
        DH(3,3,ISPIN)=-1.D0*EV
      ENDDO
!
!     ==========================================================================
!     == READ GREENS FUNCTION                                                 ==
!     ==========================================================================
PRINT*,'READING GREENS FUNCTION...'
      CALL DMFT_READGLOC(NCHI,NOMEGA,NSPIN,KBT,MU,OMEGA &
     &                        ,GLOC,GLOCLAUR1,GLOCLAUR2,GLOCLAUR3)
!
!     ==========================================================================
!     == CONSTRUCT SELF ENERGY                                                ==
!     ==========================================================================
PRINT*,'CONSTRUCTING SELF ENERGY...'
      DO ISPIN=1,NSPIN
        DO NU=1,NOMEGA
          DEDGLOC(:,:,NU,ISPIN)=DH(:,:,ISPIN)
        ENDDO
        DEDGLOCLAUR1(:,:,ISPIN)=0.D0
        DEDGLOCLAUR2(:,:,ISPIN)=DH(:,:,ISPIN)
        DEDGLOCLAUR3(:,:,ISPIN)=0.D0
        SDC(:,:,ISPIN)=0.D0
      ENDDO
!
!     ==========================================================================
!     == WRITE SELF ENERGY
!     ==========================================================================
PRINT*,'WRITING SELF ENERGY   '
      CALL DMFT_WRITESIGMA(NCHI,NOMEGA,NSPIN,KBT,MU &
     &                            ,DEDGLOC,DEDGLOCLAUR1,DEDGLOCLAUR2,DEDGLOCLAUR3,SDC)
PRINT*,'... DMFT$SOLVER DONE.'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$HFsolver(nchi,nspin,nomega,kbt,mu,omega &
     &        ,gloc,gloclaur1,gloclaur2,gloclaur3,philw &
     &        ,sigma,siglaur1,siglaur2,siglaur3)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      implicit none
      integer(4),intent(in)  :: nchi
      integer(4),intent(in)  :: nspin
      integer(4),intent(in)  :: nomega
      real(8)   ,intent(in)  :: kbt
      real(8)   ,intent(in)  :: mu
      real(8)   ,intent(in)  :: omega(nomega)
      complex(8),intent(in)  :: gloc(nchi,nchi,nomega,nspin)
      complex(8),intent(in)  :: gloclaur1(nchi,nchi,nspin)
      complex(8),intent(in)  :: gloclaur2(nchi,nchi,nspin)
      complex(8),intent(in)  :: gloclaur3(nchi,nchi,nspin)
      real(8)   ,intent(out) :: philw
      complex(8),intent(out) :: sigma(nchi,nchi,nomega,nspin)
      complex(8),intent(out) :: siglaur1(nchi,nchi,nspin)
      complex(8),intent(out) :: siglaur2(nchi,nchi,nspin)
      complex(8),intent(out) :: siglaur3(nchi,nchi,nspin)
      real(8)   ,parameter   :: upar=1.d0
      real(8)   ,parameter   :: jpar=1.d0
      real(8)                :: u(nchi,nchi,nchi,nchi)
      complex(8)             :: rho(nchi,nchi,nspin)
      complex(8)             :: rho1(nchi,nchi)
      complex(8)             :: h(nchi,nchi)
      real(8)                :: estat
      integer(4)             :: nu,ispin,i,j,k,l
      real(8)                :: svar
      complex(8)             :: csvar
!     **************************************************************************
!
!     ==========================================================================
!     == set up u-tensor
!     ==========================================================================
      u=0.d0
      do i=1,nchi
        do j=1,nchi
          u(i,j,i,j)=upar
          u(i,j,j,i)=-jpar
        enddo
      enddo
!
!     ==========================================================================
!     == calculate one-particle density matrix                                ==
!     ==========================================================================
      do ispin=1,nspin
        rho(:,:,ispin)=(0.d0,0.d0)
        do nu=1,nomega
          csvar=1.d0/cmplx(0.d0,omega(nu))
          rho(:,:,ispin)=rho(:,:,ispin)+kbt*(gloc(:,:,nu,ispin) &
    &            -csvar*(gloclaur1(:,:,ispin) &
    &                   +csvar*(gloclaur2(:,:,ispin) &
    &                           +csvar*gloclaur3(:,:,ispin))))
        enddo
        rho(:,:,ispin)=rho(:,:,ispin)+0.5d0*gloclaur1(:,:,ispin)  &
    &                               -0.25d0*kbt*gloclaur2(:,:,ispin)  &
    &                               -0.125d0*kbt**2*gloclaur3(:,:,ispin)  
      enddo
!
!     ==========================================================================
!     == calculate hartree fock potential                                     ==
!     ==========================================================================
      DO ISPIN=1,NSPIN
        RHO1=RHO(:,:,ISPIN)
        IF(NSPIN.EQ.1) RHO1=0.5D0*RHO1
        H=(0.D0,0.D0)
        DO I=1,NCHI
          DO J=1,NCHI
            DO K=1,NCHI
              DO L=1,NCHI
                SVAR=0.5D0*U(I,J,K,L) 
                ESTAT=ESTAT+SVAR*real(RHO1(K,I)*RHO1(L,J)-RHO1(K,J)*RHO1(L,I))
                H(I,K)=H(I,K)+SVAR*RHO1(L,J)
                H(J,L)=H(J,L)+SVAR*RHO1(K,I)
                H(J,K)=H(J,K)-SVAR*RHO1(L,I)
                H(I,L)=H(I,L)-SVAR*RHO1(K,J)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO NU=1,NOMEGA
          SIGMA(:,:,NU,ISPIN)=KBT*H(:,:)
        ENDDO
        IF(NSPIN.EQ.1) ESTAT=2.D0*ESTAT
      ENDdo
      siglaur1=0.d0
      siglaur2=0.d0
      siglaur3=0.d0
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_GETSIGMA()
!     **************************************************************************
!     ** PIPSI <PI_A|PSI_N> IS THE PRE-FACTOR OF LOCAL ORBITAL |CHI_A> IN     **
!     **       THE LOCAL ORBITAL EXPANSION OF |PSI_N>                         **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                      ,sigma,sigmalaur1,sigmalaur2,sigmalaur3,sigmadc &
     &                      ,diagsloc
      IMPLICIT NONE
      INTEGER(4)             :: NU
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT_getsigma')
PRINT*,'ENTERING dmft_getsigma'
      CALL DMFT_INI()
!
!     ==========================================================================
!     ==  obtain self energy
!     ==========================================================================
      CALL DMFT_READSIGMA(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC &
     &                   ,sigma,sigmaLAUR1,sigmaLAUR2,sigmaLAUR3,SigmaDC)
!
! attention: from here on delta-sigma and dc-sigma is stored
!            readsigma probably returns the full sigma and dc-sigma
!
!!$OPEN(11,FILE='DEDG.DAT')
!!$DO NU=1,NOMEGA
!!$  WRITE(11,*)OMEGA(NU),(REAL(DEDGLOC(I,I,NU,1)),AIMAG(DEDGLOC(I,I,NU,1)),I=1,NCHI)
!!$ENDDO
!!$WRITE(11,*)2.D0*OMEGA(NOMEGA)-OMEGA(NOMEGA-1) &
!!$           ,(REAL(DEDGLOCLAUR1(I,I,1)),AIMAG(DEDGLOCLAUR1(I,I,1)),I=1,NCHI)
!!$WRITE(11,*)2.D0*OMEGA(NOMEGA)-OMEGA(NOMEGA-2) &
!!$           ,(REAL(DEDGLOCLAUR2(I,I,1)),AIMAG(DEDGLOCLAUR2(I,I,1)),I=1,NCHI)
!!$WRITE(11,*)2.D0*OMEGA(NOMEGA)-OMEGA(NOMEGA-3) &
!!$           ,(REAL(DEDGLOCLAUR3(I,I,1)),AIMAG(DEDGLOCLAUR3(I,I,1)),I=1,NCHI)
!!$CLOSE(11)

!!$DEDGLOC=(0.D0,0.D0)
!!$DEDGLOCLAUR1=(0.D0,0.D0)
!!$DEDGLOCLAUR2=(0.D0,0.D0)
!!$DEDGLOCLAUR3=(0.D0,0.D0)
PRINT*,'... dmft_getsigma DONE'
                                       CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_CONSTRAINTS()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                      ,SIGMA,SIGMALAUR1,SIGMALAUR2,SIGMALAUR3,SIGMADC &
     &                      ,DIAGSLOC,PIPSI,HAMILTON,ERHO,DEDRHO,h0,hrho
      IMPLICIT NONE
      INTEGER(4)             :: IKPT,ISPIN,I
      COMPLEX(8)             :: SIG(NCHI,NCHI,NOMEGA,NSPIN)
      logical(4),parameter   :: tprint=.false.
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT_GETSIGMA')
!
!     ==========================================================================
!     ==  CONSTRUCT ONE-PARTICLE HAMILTONIANS THAT OBEY                       ==
!     ==  DENSITY MATRIX CONSTRAINT                                           ==
!     ==========================================================================
      CALL DMFT_CONSTRAINTS_TWO(H0,SIGMA)
      SIG=(0.D0,0.D0)
      CALL DMFT_CONSTRAINTS_TWO(HRHO,SIG)
      DEDRHO=HRHO-H0
!
!     ==========================================================================
!     == MAP ONTO HAMILTONIAN
!     ==========================================================================
      if(tprint) then
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            WRITE(*,FMT='(82("="),T10,"IKPT=",I5," ISPIN=",I2)')IKPT,ISPIN
            WRITE(*,FMT='("erho ",10f10.5)')erho(:,ikpt,ISPIN)
            DO I=1,NCHI
              WRITE(*,FMT='("DedRHO",80("(",2F10.5,")"))')DEDRHO(I,:,IKPT,ISPIN)
            ENDDO
            DO I=1,NCHI
              WRITE(*,FMT='("hRHO  ",80("(",2F10.5,")"))')hRHO(I,:,IKPT,ISPIN)
            ENDDO
            DO I=1,NCHI
              WRITE(*,FMT='("h0   =",80("(",2F10.5,")"))')h0(I,:,IKPT,ISPIN)
            ENDDO
          ENDDO
        ENDDO
      end if
!
!     ==========================================================================
!     == MAP ONTO HAMILTONIAN
!     ==========================================================================
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          HAMILTON(:,:,IKPT,ISPIN) &
     &                 =MATMUL(TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN))) &
     &                        ,MATMUL(H0(:,:,IKPT,ISPIN),PIPSI(:,:,IKPT,ISPIN)))
        ENDDO
      ENDDO

                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_CONSTRAINTS_TWO(H0,SIG)
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                      ,DIAGSLOC,PIPSI,ERHO
      IMPLICIT NONE
      COMPLEX(8),INTENT(INOUT) :: H0(NCHI,NCHI,NKPTL,NSPIN)
      COMPLEX(8),INTENT(IN)    :: SIG(NCHI,NCHI,NOMEGA,NSPIN)
      REAL(8)   ,PARAMETER     :: AMIX=1.D-1
      REAL(8)   ,PARAMETER     :: TOL=1.D-6
      INTEGER(4),PARAMETER     :: NITER=100
      LOGICAL(4)               :: CONVG
      REAL(8)                  :: MAXDEV
      INTEGER(4)               :: NU
      COMPLEX(8),ALLOCATABLE   :: GREENINV(:,:)
      COMPLEX(8),ALLOCATABLE   :: GREEN(:,:)
      COMPLEX(8)               :: DEVRHO(NCHI,NCHI)
      INTEGER(4)               :: NTASKS_K,THISTASK_K
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: CSVAR
      INTEGER(4)               :: IKPT,ISPIN,IB,I,ITER
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT_GETSIGMA')
!       
!     ==========================================================================
!     ==  ADJUST H0 TO OBEY DENSITY MATRIX CONSTRAINT                         ==
!     ==========================================================================
      ALLOCATE(GREENINV(NB,NB))
      ALLOCATE(GREEN(NB,NB))
      DO ITER=1,NITER
        MAXDEV=0.D0
!       
!       ========================================================================
!       ==  DEVIATION FROM TARGET DENSITY MATRIX                              ==
!       ========================================================================
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            DEVRHO=(0.D0,0.D0)
            DO NU=1,NOMEGA
        
!             == CONSTRUCT LATTICE GREENS FUNCTION =============================
              GREENINV(:,:)=-MATMUL(TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN))) &
     &                          ,MATMUL(H0(:,:,IKPT,ISPIN)+SIG(:,:,NU,ISPIN) &
     &                                   ,PIPSI(:,:,IKPT,ISPIN)))
              CSVAR=CI*OMEGA(NU)+MU
              DO IB=1,NB
                GREENINV(IB,IB)=GREENINV(IB,IB)+CSVAR
              ENDDO
              CALL LIB$INVERTC8(NB,GREENINV,GREEN)
!       
!             == SUBTRACT REFERENCE GREENS FUNCTION ============================
              DO IB=1,NB
                CSVAR=1.D0/(CI*OMEGA(NU)-ERHO(IB,IKPT,ISPIN))
                GREEN(IB,IB)=GREEN(IB,IB)-CSVAR
              ENDDO 
!       
!             == ADD UP DEVIATION OF THE DENSITY MATRICES ======================
              DEVRHO(:,:)=DEVRHO(:,:)+KBT*MATMUL(PIPSI(:,:,IKPT,ISPIN) &
        &               ,MATMUL(GREEN,TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN)))))
            ENDDO
!
!           == INCLUDE NEGATIVE FREQUENCIES ====================================
            DEVRHO(:,:)=DEVRHO(:,:)+CONJG(TRANSPOSE(DEVRHO(:,:)))
            MAXDEV=MAX(MAXDEV,MAXVAL(ABS(DEVRHO)))
!
!           == MIX INTO H0 =====================================================
            H0(:,:,IKPT,ISPIN)=H0(:,:,IKPT,ISPIN)+AMIX*DEVRHO(:,:)
          ENDDO
        ENDDO
!
!       ========================================================================
!       == CHECK CONVERGENCE
!       ========================================================================
        PRINT*,'MAXDEV',ITER,MAXDEV
        CONVG=MAXDEV.LT.TOL
        IF(CONVG) EXIT
      ENDDO
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('LOOP NOT CONVERGED')
        CALL ERROR$STOP('DMFT_CONSTRAINTS_TWO')
      END IF

                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_gloc(H0,SIG,gloc,laur)
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                       ,wkptl,smat,sinv,dsiglaur0
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN)    :: H0(NCHI,NCHI,NKPTL,NSPIN)
      COMPLEX(8),INTENT(IN)    :: SIG(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(out)   :: gloc(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(out)   :: laur(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)               :: Ginv(nchi,nchi)
      COMPLEX(8)               :: mat(nchi,nchi)
      COMPLEX(8)               :: G(nchi,nchi)
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: CSVAR
      INTEGER(4)               :: IKPT,ISPIN,nu,i,j
      INTEGER(4)               :: nfil
      logical(4),parameter     :: tprint=.false.
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT_GETSIGMA')
      DO ISPIN=1,NSPIN
!       ========================================================================
!       == LAURENT EXPANSION OF THE LOCAL GREENS FUNCTION FROM 0 TO 2
!       == REMEMBER THAT MU=0 !!!
!       ========================================================================
        LAUR(:,:,1,ISPIN)=SMAT(:,:,ISPIN)
        LAUR(:,:,2,ISPIN)=(0.d0,0.d0)
        LAUR(:,:,3,ISPIN)=(0.d0,0.d0)
        do ikpt=1,nkptl
          MAT=h0(:,:,ikpt,ispin)+DSIGLAUR0(:,:,ISPIN,1)
          mat=MATMUL(SINV(:,:,ikpt,ISPIN),MATMUL(MAT,SINV(:,:,ikpt,ISPIN)))
          LAUR(:,:,2,ISPIN)=LAUR(:,:,2,ISPIN)+wkptl(ikpt)*mat
          MAT=h0(:,:,ikpt,ispin)+DSIGLAUR0(:,:,ISPIN,1)
          mat=matmul(mat,matmul(sinv(:,:,ikpt,ispin),mat))
          mat=mat+DSIGLAUR0(:,:,ISPIN,2)
          mat=MATMUL(SINV(:,:,ikpt,ISPIN),MATMUL(MAT,SINV(:,:,ikpt,ISPIN)))
          LAUR(:,:,3,ISPIN)=LAUR(:,:,3,ISPIN)+wkptl(ikpt)*mat
        enddo
!
!       ========================================================================
!       == LOCAL GREENS FUNCTION ON THE MATSUBARA GRID
!       ========================================================================
        DO IKPT=1,NKPTL
          DO NU=1,NOMEGA
            CSVAR=CI*OMEGA(NU)+MU
            GINV(:,:)=SMAT(:,:,ISPIN)*CSVAR-H0(:,:,IKPT,ISPIN)-SIG(:,:,NU,ISPIN)
            CALL LIB$INVERTC8(NCHI,GINV,G)
            GLOC(:,:,NU,ISPIN)=GLOC(:,:,NU,ISPIN)+WKPTL(IKPT)*G
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == LOCAL GREENS FUNCTION ON THE MATSUBARA GRID                          ==
!     ==========================================================================
      IF(TPRINT) THEN
        DO ISPIN=1,NSPIN  
          DO I=1,3
            WRITE(*,FMT='(82("="),T10," LAUR(",I1,") IKPT,ISPIN",I4)'),I,ISPIN
            DO J=1,3
               WRITE(*,FMT='(10("(",2F10.5,")"))')LAUR(J,:,I,ISPIN)
            ENDDO
          ENDDO
        ENDDO
!
        CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,'DMFT_GLOC.DAT')
        CALL FILEHANDLER$UNIT('HOOK',NFIL)
        DO NU=1,NOMEGA
          WRITE(NFIL,*)OMEGA(NU),REAL(GLOC(:,:,NU,:)),AIMAG(GLOC(:,:,NU,:))
        ENDDO
        CALL FILEHANDLER$CLOSE('HOOK')
        CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOKERROR')
        CALL ERROR$MSG('SCHEDULED STOP AFTER PRINTING FILE DMFT_GLOC.DAT')
        CALL ERROR$STOP('DMFT_GLOC')
      END IF
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$ADDTOEIGVAL()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY  : TON,NKPTL,NSPIN,NB,GAMMA0
      IMPLICIT NONE
      REAL(8)        :: EIG(NB,NKPTL,NSPIN)
      INTEGER(4)     :: IKPT,ISPIN,IB
!     **************************************************************************
      IF(.NOT.TON) RETURN
      CALL WAVES_DYNOCCGETR8A('EPSILON',NB*NKPTL*NSPIN,EIG)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          DO IB=1,NB
            EIG(IB,IKPT,ISPIN)=EIG(IB,IKPT,ISPIN) &
     &                        +REAL(GAMMA0(IB,IB,IKPT,ISPIN))
          ENDDO
        ENDDO
      ENDDO
OPEN(11,FILE='EIGS.DAT')
DO IKPT=1,NKPTL
  DO IB=1,NB
    WRITE(11,*)EIG(IB,IKPT,1)*27.211D0,REAL(GAMMA0(IB,IB,IKPT,1))*27.211D0 &
  &                                   ,AIMAG(GAMMA0(IB,IB,IKPT,1))*27.211D0
ENDDO
ENDDO
CLOSE(11)
      CALL WAVES_DYNOCCSETR8A('EPSILON',NB*NKPTL*NSPIN,EIG)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$ADDTOHPSI()
!     **************************************************************************
!     ** 
!     **************************************************************************
      USE DMFT_MODULE, ONLY  : TON,NKPTL,NSPIN,NB,NDIM,DENMAT,GAMMA0
      USE WAVES_MODULE, ONLY : GSET,WAVES_SELECTWV,THIS
      IMPLICIT NONE
      REAL(8)   ,PARAMETER   :: MINOCC=1.D-1
      COMPLEX(8),ALLOCATABLE :: AMAT(:,:)
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      INTEGER(4)             :: NBH
      INTEGER(4)             :: NGL
      INTEGER(4)             :: IKPT,ISPIN,IB,IDIM
!     **************************************************************************
      IF(.NOT.TON) RETURN
PRINT*,'ENTERING DMFT$ADDTOHPSI'
                              CALL TRACE$PUSH('DMFT$ADDTOHAPSI')
      ALLOCATE(OCC(NB,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NB*NKPTL*NSPIN,OCC)
      ALLOCATE(AMAT(NB,NB))
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          IF(THIS%NB.NE.NB) THEN
            CALL ERROR$MSG('INCONSISTENT NUMBER OF STATES IN WAVES AND DYNOCC')
            CALL ERROR$I4VAL('NB IN DYNOCC',NB)
            CALL ERROR$I4VAL('NB IN WAVES ',THIS%NB)
            CALL ERROR$STOP('DMFT$GREEN')
          END IF
          NBH=THIS%NBH
          NGL=GSET%NGL
!
          AMAT=MATMUL(DENMAT(:,:,IKPT,ISPIN),GAMMA0(:,:,IKPT,ISPIN))
          DO IB=1,NB          
            AMAT(:,IB)=AMAT(:,IB)/(OCC(IB,IKPT,ISPIN)+MINOCC)
          ENDDO
!         == MAY PRODUCE A MEMORY SPIKE... =====================================
          CALL WAVES_ADDPSI(NGL,NDIM,NBH,NB,THIS%HPSI,THIS%PSI0,AMAT)
        ENDDO
      ENDDO
PRINT*,'.... DMFT$ADDTOHPSI DONE'
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$DOS(NOMEGA,OMEGA,NCHI,NSPIN,KBT,MU,GLOC)
!     **************************************************************************
!     ** CALCULATES A KIND OF DENSITY OF STATES FROM DN/DCHEMPOT              **
!     ** THIS SHOULD BE THE THERMALLY BROADENED SPECTRAL FUNCTION??           **
!     ** IT ALSO CALLS DMFT$DOS1 FOR COMPARISON                               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NOMEGA
      REAL(8)   ,INTENT(IN)  :: OMEGA(NOMEGA)
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NSPIN
      REAL(8)   ,INTENT(IN)  :: KBT
      REAL(8)   ,INTENT(IN)  :: MU
      COMPLEX(8),INTENT(IN)  :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      INTEGER(4),PARAMETER   :: NE=10000
      REAL(8)                :: EMIN
      REAL(8)                :: EMAX
      REAL(8)                :: DE
      REAL(8)                :: E
      COMPLEX(8)             :: DOS(NCHI,NCHI,NE,NSPIN)
      COMPLEX(8)             :: ONE(NCHI,NCHI)
      COMPLEX(8)             :: A(NCHI,NCHI)
      COMPLEX(8)             :: B(NCHI,NCHI)
      COMPLEX(8)             :: GINV0(NCHI,NCHI)
      INTEGER(4)             :: I,NU,ISPIN,IE
      REAL(8)                :: EV
!     **************************************************************************
      ONE=(0.D0,0.D0)
      DO I=1,NCHI
        ONE(I,I)=(1.D0,0.D0)
      ENDDO
      CALL CONSTANTS$GET('EV',EV)
!
!     ==========================================================================
!     == DEFINE ENERGY GRID                                                   ==
!     ==========================================================================
      EMIN=-20.D0*EV
      EMAX=+40.D0*EV
      DE=(EMAX-EMIN)/REAL(NE-1,KIND=8)
!
!     ==========================================================================
!     == CONSTRUCT DENSITY OF STATES                                          ==
!     ==========================================================================
      DOS=(0.D0,0.D0)
      DO ISPIN=1,NSPIN
        DO NU=1,NOMEGA
          CALL LIB$INVERTC8(NCHI,GLOC(:,:,NU,ISPIN),GINV0)
          DO IE=1,NE
            E=EMIN+DE*REAL(IE-1,KIND=8)
            A=GINV0+ONE*(E-MU)
            A=MATMUL(A,A)
            CALL LIB$INVERTC8(NCHI,A,B)
            DOS(:,:,IE,ISPIN)=DOS(:,:,IE,ISPIN)-KBT*B
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == WRITE RESULT TO FILE                                                 ==
!     ==========================================================================
      OPEN(11,FILE='DOS_RE.DAT')
      DO IE=1,NE
        E=EMIN+DE*REAL(IE-1,KIND=8)
        WRITE(11,*)E/EV,REAL(DOS(:,:,IE,:))
      ENDDO
      CLOSE(11)
      OPEN(11,FILE='DOS_IM.DAT')
      DO IE=1,NE
        E=EMIN+DE*REAL(IE-1,KIND=8)
        WRITE(11,*)E/EV,AIMAG(DOS(:,:,IE,:))
      ENDDO
      CLOSE(11)
!
!     ==========================================================================
!     == 
!     ==========================================================================
      CALL DMFT$DOS1()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$DOS1()
!     **************************************************************************
!     ** DENSITY OF STATES OF THE EFFECTIVE HAMILTONIAN PROJECTED ONTO        **
!     ** THE CORRELATED SUBSPACE. USED FOR COMPARING WITH THE FUNCTION FROM   **
!     ** DMFT$ DOS OBTAINED FROM THE GREENS FUNCTION                          **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NB,NCHI,NKPTL,NSPIN,NDIM,IPROOFCHI,WKPTL,KBT
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      INTEGER(8),PARAMETER   :: NE=1000
      REAL(8)                :: EMIN,EMAX
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: NGL
      COMPLEX(8)             :: CSVAR
      COMPLEX(8),ALLOCATABLE :: HAMILTON(:,:)
      COMPLEX(8)             :: PIPSI(NCHI,NB)
      COMPLEX(8)             :: DOS(NCHI,NCHI,NE)
      COMPLEX(8)             :: D(NCHI,NCHI)
      COMPLEX(8)             :: UMAT(NB,NB)
      REAL(8)                :: EIG(NB)
      REAL(8)                :: EV
      REAL(8)                :: E
      REAL(8)                :: DE
      INTEGER(4)             :: IKPT,ISPIN,ICHI,IPRO,IBH,IE,IB
!     **************************************************************************
      CALL CONSTANTS$GET('EV',EV)
      EMIN=-20.D0*EV
      EMAX=+40.D0*EV
      DE=(EMAX-EMIN)/REAL(NE-1,KIND=8)

      ALLOCATE(HAMILTON(NB,NB))
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)      
          IF(THIS%NB.NE.NB) THEN
            CALL ERROR$MSG('INCONSISTENT NUMBER OF STATES IN WAVES AND DYNOCC')
            CALL ERROR$I4VAL('NB IN DYNOCC',NB)
            CALL ERROR$I4VAL('NB IN WAVES ',THIS%NB)
            CALL ERROR$STOP('DMFT$GREEN')
          END IF
          NBH=THIS%NBH
          NGL=GSET%NGL
!
!         ======================================================================
!         ==  CALCULATE ONE-PARTICLE HAMILTONIAN                              ==
!         ======================================================================
          CALL WAVES_OVERLAP(.FALSE.,NGL,NDIM,NBH,NB,THIS%PSI0,THIS%HPSI &
     &                     ,HAMILTON)
!
!         ======================================================================
!         ==  DETERMINE ORBITAL PROJECTIONS                                   ==
!         ======================================================================
          DO ICHI=1,NCHI
            IPRO=IPROOFCHI(ICHI)
            DO IBH=1,NBH
              IF(NBH.NE.NB) THEN
                 PIPSI(ICHI,2*IBH-1)=REAL(THIS%TBC(1,IBH,IPRO))
                 PIPSI(ICHI,2*IBH)  =AIMAG(THIS%TBC(1,IBH,IPRO))
              ELSE
                 PIPSI(ICHI,IBH)  =THIS%TBC(1,IBH,IPRO)
              END IF
            ENDDO
          ENDDO
!
!         ======================================================================
!         ==  DIAGONALIZE  (CURRENTLY ONLY DIAGONAL ELEMENTS)                ==
!         ======================================================================
!!$          UMAT(:,:)=(0.D0,0.D0)
!!$          DO IB=1,NB
!!$            EIG(IB)=REAL(HAMILTON(IB,IB))
!!$            UMAT(IB,IB)=(1.D0,0.D0)
!!$          ENDDO
          CALL LIB$DIAGC8(NB,HAMILTON,EIG,UMAT)
          PIPSI=MATMUL(PIPSI,UMAT)
!
!         ======================================================================
!         ==                                                                  ==
!         ======================================================================
          DO IB=1,NB
            IE=1+NINT((EIG(IB)-EMIN)/DE)
            IF(IE.LT.1.OR.IE.GT.NE) CYCLE
            DO ICHI=1,NCHI
              D(:,ICHI)=PIPSI(:,IB)*CONJG(PIPSI(ICHI,IB))
            ENDDO
            DOS(:,:,IE)=DOS(:,:,IE)+WKPTL(IKPT)*D(:,:)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == THERMAL BROADENING                                                   ==
!     ==========================================================================
!!$      DO IE=-100,100
!!$        E=DE*REAL(IE)
!!$        F(IE)=1.D0/(KBT*COSH(0.5D0*E/KBT)**2)
!!$      ENDDO
!
!     ==========================================================================
!     == WRITE RESULT TO FILE                                                 ==
!     ==========================================================================
      OPEN(11,FILE='DOS1_RE.DAT')
      DO IE=1,NE
        E=EMIN+DE*REAL(IE-1,KIND=8)
        WRITE(11,*)E/EV,REAL(DOS(:,:,IE))
      ENDDO
      CLOSE(11)
      OPEN(11,FILE='DOS1_IM.DAT')
      DO IE=1,NE
        E=EMIN+DE*REAL(IE-1,KIND=8)
        WRITE(11,*)E/EV,AIMAG(DOS(:,:,IE))
      ENDDO
      CLOSE(11)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$TEST1(NOMEGA,OMEGA,NB,NKPT,NCHI,NSPIN,KBT,WKPTL &
     &                     ,GLOC,DENMAT,GLOCLAUR1,GLOCLAUR2,GLOCLAUR3)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NOMEGA
      REAL(8)   ,INTENT(IN)  :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(IN)  :: KBT
      INTEGER(4),INTENT(IN)  :: NKPT
      INTEGER(4),INTENT(IN)  :: NSPIN
      INTEGER(4),INTENT(IN)  :: NB
      INTEGER(4),INTENT(IN)  :: NCHI
      REAL(8)   ,INTENT(IN)  :: WKPTL(NKPT)
      COMPLEX(8),INTENT(IN)  :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(IN)  :: DENMAT(NB,NB,NKPT,NSPIN)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR1(NCHI,NCHI)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR2(NCHI,NCHI)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR3(NCHI,NCHI)
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      REAL(8)                :: NEL
      INTEGER(4)             :: IB,IKPT,ISPIN,I,J,NU
      REAL(8)                :: WEIGHT(NOMEGA)
      COMPLEX(8),ALLOCATABLE :: B0(:,:),B1(:,:),S(:,:),H(:,:),MAT(:,:)
      COMPLEX(8)             :: A0,A1,A2,CSVAR
      REAL(8)                :: SVAR
      INTEGER(4)             :: ISVAR
!     **************************************************************************
      PRINT*,'ENTERING DMFT$TEST1....'
!
!     ==========================================================================
!     ==  CALCULATE NUMBER OF ELECTRONS FROM LATTICE DENSITY MATRIX           ==
!     ==========================================================================
      NEL=0.D0
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          DO IB=1,NB
            NEL=NEL+WKPTL(IKPT)*REAL(DENMAT(IB,IB,IKPT,ISPIN))
          ENDDO
        ENDDO
      ENDDO
      PRINT*,'NUMBER OF ELECTRONS FROM LATTICE DENSITY MATRIX:',NEL
!
!     ==========================================================================
!     ==  CALCULATE NUMBER OF ELECTRONS FROM LATTICE DENSITY MATRIX           ==
!     ==========================================================================
      PRINT*,'BAND OCCUPATIONS'
      PRINT*,'SHOULD BEGIN WITH TWO AND END WITH ZERO'
      DO IB=1,NB
        SVAR=0.D0
        DO IKPT=1,NKPT
          DO ISPIN=1,NSPIN
            SVAR=SVAR+WKPTL(IKPT)*REAL(DENMAT(IB,IB,IKPT,ISPIN))
          ENDDO
        ENDDO
        WRITE(*,FMT='(10("(",2F10.5,")"))')SVAR
      ENDDO
!
!     ==========================================================================
!     ==  ONSITE DENSITY MATRIX                                               ==
!     ==========================================================================
      WEIGHT(:)=(OMEGA(:)/OMEGA(NOMEGA))**4
      ISVAR=NOMEGA/2
      WEIGHT(1+ISVAR:NOMEGA-ISVAR)=0.D0
      ALLOCATE(B0(NCHI,NCHI))
      ALLOCATE(B1(NCHI,NCHI))
      ALLOCATE(S(NCHI,NCHI))
      ALLOCATE(H(NCHI,NCHI))
      ALLOCATE(MAT(NCHI,NCHI))
      A2=SUM(WEIGHT(:)*(-CI*OMEGA(:))**2)
      A1=SUM(WEIGHT(:)*(-CI*OMEGA(:)))
      A0=SUM(WEIGHT(:))
      B0=(0.D0,0.D0)
      B1=(0.D0,0.D0)
      DO NU=1,NOMEGA
        IF(WEIGHT(NU).EQ.0.D0) CYCLE
        CALL LIB$INVERTC8(NCHI,GLOC(:,:,NU,1),MAT)
        B0=B0+WEIGHT(NU)*MAT
        B1=B1+WEIGHT(NU)*MAT*(-CI*OMEGA(NU))
      ENDDO
      CSVAR=1.D0/(A2*A0-A1*A1)
      S=(-B1*A0+B0*A1)*CSVAR
      H=(+B1*A1-B0*A2)*CSVAR
!
      PRINT*,'OVERLAP FROM TAILS OF G^-1'
      DO I=1,NCHI
        WRITE(*,FMT='("S",10("(",2F10.5,")"))')S(I,:)
      ENDDO
      PRINT*,'ONSITE HAMILTONIAN FROM TAILS OF G^-1'
      DO I=1,NCHI
        WRITE(*,FMT='("H",10("(",2F10.5,")"))')H(I,:)
      ENDDO
      PRINT*,'REAL(<PI|PSI><PSI|PI>)'
      PRINT*,'MULTIPLIED WITH THE OVERLAP MATRIX THE UNIT MATRIX SHOULD RESULT'
      DO I=1,NCHI
        WRITE(*,FMT='("RE(1/S)",10F10.5)')REAL(GLOCLAUR1(I,:))
      ENDDO
      PRINT*,'AIMAG(<PI|PSI><PSI|PI>)'
      DO I=1,NCHI
        WRITE(*,FMT='("IM(1/S)",10F10.5)')AIMAG(GLOCLAUR1(I,:))
      ENDDO

      DO I=1,NCHI
        WRITE(*,FMT='(10("(",2F10.5,")"))') &
     &                  (KBT*SUM(GLOC(I,J,:,1))*2.D0/REAL(NSPIN),J=1,NCHI)
      ENDDO

      CALL DMFT$TEST1A()
!
!     ==========================================================================
!     ==  WRITE HYBRIDIZATION FUNCTION TO FILE
!     ==========================================================================
      PRINT*,'WRITING DELTA.DAT'
      OPEN(UNIT=11,FILE='DELTA.DAT')
      DO NU=1,NOMEGA
        CALL LIB$INVERTC8(NCHI,GLOC(:,:,NU,1),MAT)
        MAT=CI*OMEGA(NU)*S-H-MAT
        WRITE(11,*)OMEGA(NU),(REAL(MAT(I,I)),AIMAG(MAT(I,I)),I=1,3)
      ENDDO
      CLOSE(11)
      PRINT*,'... DELTA.DAT WRITTEN'
!
!     ==========================================================================
!     ==  WRITE LOCAL GREENS FUNCTION TO FILE                                 ==
!     ==========================================================================
      PRINT*,'WRITING GLOC.DAT'
      OPEN(UNIT=11,FILE='GLOC.DAT')
      DO NU=1,NOMEGA
        B0=CI*OMEGA(NU)*S-H
        CALL LIB$INVERTC8(NCHI,B0,MAT)
!       == EXTRAPOLATION FROM CALCULATED EXPANSION COEFFICIENTS
        CSVAR=1.D0/(CI*OMEGA(NU))
        MAT=GLOCLAUR1*CSVAR+GLOCLAUR2*CSVAR**2+GLOCLAUR3*CSVAR**3
        WRITE(11,*)OMEGA(NU) &
    &            ,(REAL(GLOC(I,I,NU,1)),AIMAG(GLOC(I,I,NU,1)) &
    &             ,REAL(MAT(I,I)),AIMAG(MAT(I,I)) &
    &                                            ,I=1,3)
      ENDDO
      CLOSE(11)
      PRINT*,'... GLOC.DAT WRITTEN'
!
!     ==========================================================================
!     ==  FINISHED                                                            ==
!     ==========================================================================
      PRINT*,'.... DMFT$TEST1 DONE'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$TEST1A()
!     **************************************************************************
!     ** CALCULATE OVERLAP MATRIX FROM TAILED ORBITALS                        **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR,LNX,LOX
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: IAT=2
      INTEGER(4)           :: ISP
      INTEGER(4)           :: I
      INTEGER(4)           :: LMNXT
      INTEGER(4)           :: LMNX
      REAL(8)   ,ALLOCATABLE:: OVERLAP(:,:)
!     *************************************************************************
      ISP=ISPECIES(IAT)
      LMNX=SUM(2*LOX(:LNX(ISP),ISP)+1)
      LMNXT=POTPAR(ISP)%TAILED%LMNX
      ALLOCATE(OVERLAP(LMNX,LMNX))
      CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,1 &
   &               ,LMNXT,LMNXT,POTPAR(ISP)%TAILED%OVERLAP,LMNX,LMNX,OVERLAP)
      PRINT*,'OVERLAP FOR ATOM ',IAT, ' FROM TAILED ORBITALS'
      DO I=1,LMNX
        WRITE(*,FMT='(20F10.5)')OVERLAP(I,:)
      ENDDO
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE DMFT$PLOTSIGMA(FILE)
      USE DMFT_MODULE, ONLY : NOMEGA,OMEGA,DSIG0
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)              :: NU
      INTEGER(4)  ,PARAMETER  :: NFIL=11
!     **************************************************************************
      OPEN(NFIL,FILE=FILE)
      REWIND NFIL
      DO NU=1,NOMEGA
        WRITE(NFIL,*)OMEGA(NU),DSIG0(:,:,NU,:)
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE DMFT$PLOTGLOC(FILE,GLOC)
      USE DMFT_MODULE, ONLY : NOMEGA,OMEGA,NCHI,NSPIN
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      COMPLEX(8)  ,INTENT(IN) :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      INTEGER(4)              :: NU
      INTEGER(4)  ,PARAMETER  :: NFIL=11
!     **************************************************************************
      OPEN(NFIL,FILE=FILE)
      REWIND NFIL
      DO NU=1,NOMEGA
        WRITE(NFIL,*)OMEGA(NU),REAL(GLOC(:,:,NU,:)),AIMAG(GLOC(:,:,NU,:))
      ENDDO
      CLOSE(NFIL)
      RETURN
      END

