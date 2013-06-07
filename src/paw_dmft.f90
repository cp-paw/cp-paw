
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
LOGICAL(4),PARAMETER   :: TON=.TRUE.
LOGICAL(4),SAVE        :: TINI=.FALSE.
REAL(8)   ,PARAMETER   :: AMIX=1.D-2
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
REAL(8)   ,ALLOCATABLE :: ERHO(:,:,:)       !(NB,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: H0(:,:,:,:)       !(NCHI,NCHI,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: HRHO(:,:,:,:)     !(NCHI,NCHI,NKPTL,NSPIN)
!__MATRIX TO DIAGONALIZE OVERLAP MATRIX_________________________________________
COMPLEX(8),ALLOCATABLE :: DIAGSLOC(:,:,:) !(NCHI,NCHI,NSPIN)
COMPLEX(8),ALLOCATABLE :: SMAT(:,:,:,:)   !(NCHI,NCHI,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: SINV(:,:,:,:)   !(NCHI,NCHI,NKPTL,NSPIN)
!== LOCAL GREENS FUNCTION (WITH LAURENT EXPANSION COEFFICIENTS) ================
COMPLEX(8),ALLOCATABLE :: GLOC(:,:,:,:)     !(NCHI,NCHI,NOMEGA,NSPIN)
COMPLEX(8),ALLOCATABLE :: GLOCLAUR(:,:,:,:)  !(NCHI,NCHI,3,NSPIN)
!== SELF ENERGY ================================================================
COMPLEX(8),ALLOCATABLE :: SIGMA(:,:,:,:)  !(NCHI,NCHI,NOMEGA,NSPIN)
COMPLEX(8),ALLOCATABLE :: SIGLAUR(:,:,:,:) !(NCHI,NCHI,3,NSPIN)
COMPLEX(8),ALLOCATABLE :: SIGMADC(:,:,:) !(NCHI,NCHI,NSPIN) (DOUBLE COUNTING)
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
      ALLOCATE(ERHO(NB,NKPTL,NSPIN))
      ALLOCATE(PIPSI(NCHI,NB,NKPTL,NSPIN))
      ALLOCATE(HRHO(NCHI,NCHI,NKPTL,NSPIN))
!     == H0 IS PURPOSELY NOT ALLOCATED
      ALLOCATE(SIGMADC(NCHI,NCHI,NSPIN))
      ALLOCATE(SMAT(NCHI,NCHI,NKPTL,NSPIN))
      ALLOCATE(SINV(NCHI,NCHI,NKPTL,NSPIN))
      ALLOCATE(GLOC(NCHI,NCHI,NOMEGA,NSPIN))
      ALLOCATE(GLOCLAUR(NCHI,NCHI,3,NSPIN))
      ALLOCATE(DIAGSLOC(NCHI,NCHI,NSPIN))
      ALLOCATE(SIGMA(NCHI,NCHI,NOMEGA,NSPIN))
      ALLOCATE(SIGLAUR(NCHI,NCHI,3,NSPIN))
      SIGMA=(0.D0,0.D0)
      SIGLAUR=(0.D0,0.D0)
      SIGMADC=(0.D0,0.D0)
      HRHO=(0.D0,0.D0)
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
     &              ,WKPTL,H0,HRHO,SIGMA,SIGLAUR,GLOC,GLOCLAUR,DIAGSLOC
      USE MPE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: I
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT$GREEN')
      CALL DMFT_INI()
!
!     ==========================================================================
!     ==  COLLECT DFT HAMILTONIAN  (HAMILTONIAN INITIALIZED TO ERHO)          ==
!     ==========================================================================
      CALL DMFT$COLLECTHAMILTONIAN()
!
!     ==========================================================================
!     ==  CONSTRUCT NON-INTERTACTING HAMILTONIAN THAT PRODUCES THE CORRECT    ==
!     == ONE-PARTICLE DENSITY MATRIX
!     ==========================================================================
      CALL DMFT_HRHO()
!
      IF(.NOT.ALLOCATED(H0)) THEN
        ALLOCATE(H0(NCHI,NCHI,NKPTL,NSPIN))
        H0=HRHO
      END IF
!
!     ==========================================================================
!     ==  TRANSFORMATION ONTO A ORTHONORMAL CORRELATED BASIS SET              ==
!     ==    |CHIORTHO>   =|CHINONORTHO>*DIAGSLOC                              ==
!     ==========================================================================
      CALL DMFT_DIAGSLOC()
!
!     ==========================================================================
!     == DETERMINE LOCAL GREENS FUNCTION                                      ==
!     ==========================================================================
      MU=0.D0
DO I=1,3
      CALL DMFT_GLOC(H0,SIGMA,SIGLAUR,GLOC,GLOCLAUR)
      CALL DMFT$PLOTGLOC(-'GLOC.DAT')
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
!     ==  READ SELF ENERGY FROM FILE
!     ==========================================================================
      CALL DMFT_GETSIGMA()
      CALL DMFT$PLOTSIGMA(-'SIG.DAT')
!
!     ==========================================================================
!     ==  CONSTRAINTS
!     ==========================================================================
      CALL DMFT_CONSTRAINTS()
ENDDO

STOP 'END OF LOOP. STOPPING.'
                                       CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_DIAGSLOC()
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
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NSPIN,WKPTL,PIPSI,DIAGSLOC,SMAT,SINV
      IMPLICIT NONE
      INTEGER(4)      :: IKPT,ISPIN,I
      COMPLEX(8)      :: MAT(NCHI,NCHI)
      REAL(8)         :: FLOC(NCHI)
      LOGICAL(4)      :: TTEST=.FALSE.
!     **************************************************************************
      DO ISPIN=1,NSPIN
!       ========================================================================
!       == INVERSE OVERLAP MATRIX FROM SUM_N <PI|PSI><PSI|PI>                 ==
!       ========================================================================
        MAT=(0.D0,0.D0)
        DO IKPT=1,NKPTL
          SINV(:,:,IKPT,ISPIN)=MATMUL(PIPSI(:,:,IKPT,ISPIN) &
     &                                ,CONJG(TRANSPOSE(PIPSI(:,:,IKPT,ISPIN))))
          MAT=MAT+WKPTL(IKPT)*SINV(:,:,IKPT,ISPIN)
        ENDDO
!
!       ========================================================================
!       == OVERLAP MATRIX                                                     **
!       ========================================================================
        DO IKPT=1,NKPTL
          CALL LIB$INVERTC8(NCHI,SINV(:,:,IKPT,ISPIN),SMAT(:,:,IKPT,ISPIN))
        ENDDO
!
!       ========================================================================
!       == OBTAIN MATRIX THAT MAKE THE OVERLAP EQUAL TO THE UNIT MATRIX       ==
!       ========================================================================
        CALL LIB$DIAGC8(NCHI,MAT,FLOC,DIAGSLOC(:,:,ISPIN))
        DO I=1,NCHI
          DIAGSLOC(:,I,ISPIN)=DIAGSLOC(:,I,ISPIN)/SQRT(FLOC(I))
        ENDDO
!
!       ========================================================================
!       == OPTIONAL TESTS
!       ========================================================================
        IF(TTEST) THEN
          MAT=MATMUL(CONJG(TRANSPOSE(DIAGSLOC(:,:,ISPIN))) &
                                ,MATMUL(MAT,DIAGSLOC(:,:,ISPIN)))
          DO I=1,NCHI
            WRITE(*,FMT='("TEST",10("(",2F10.5,")"))')MAT(I,:)
          ENDDO
          CALL ERROR$MSG('SCHEDULED STOP AFTER TEST')
          CALL ERROR$STOP('DMFT_DIAGSLOC')
        END IF
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$COLLECTHAMILTONIAN()
!     **************************************************************************
!     ** COLLECTS THE HAMILTONIAN AND STORES IT ON THE MODULE                 **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,IPROOFCHI,KBT &
     &                       ,ERHO,PIPSI
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: NGL
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IPRO,IB
      REAL(8)                :: F(NB,NKPTL,NSPIN)
      REAL(8)                :: SVAR
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                      CALL TRACE$PUSH('DMFT$COLLECTHAMILTONIAN')
!
!     ==========================================================================
!     ==  
!     ==========================================================================
      CALL DMFT$COLLECTOCCUPATIONS(NKPTL,NSPIN,NB,F)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          DO IB=1,NB
            CALL DMFT_EOFF(KBT,F(IB,IKPT,ISPIN),ERHO(IB,IKPT,ISPIN))
          ENDDO
        ENDDO
      ENDDO 
PRINT*,'BOUNDS OF SPECTRUM [EV] ',MINVAL(ERHO)*27.211D0,MAXVAL(ERHO)*27.211D0
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
      INTEGER(4)             :: NU,ISPIN,I,J,N
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
        DO N=1,3
          GLOCPRIME=MATMUL(TPLUS,MATMUL(GLOCLAUR(:,:,N,ISPIN),T))
          WRITE(NFIL,*)((REAL(GLOCPRIME(I,J)),AIMAG(GLOCPRIME(I,J)) &
     &                  ,I=1,NCHI),J=1,NCHI)
        ENDDO
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_READGLOC(NCHI,NOMEGA,NSPIN,KBT,MU,OMEGA &
     &                        ,GLOC,GLOCLAUR)
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
      COMPLEX(8),INTENT(OUT) :: GLOCLAUR(NCHI,NCHI,3,NSPIN)
      REAL(8)                :: GLOC_RE(NCHI,NCHI),GLOC_IM(NCHI,NCHI)
      INTEGER(4)             :: NCHI_,NOMEGA_,NSPIN_
      INTEGER(4)             :: NU,ISPIN,I,J,N
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
        DO N=1,3
          READ(NFIL,*)((GLOC_RE(I,J),GLOC_IM(I,J),I=1,NCHI),J=1,NCHI)
          GLOCLAUR(:,:,N,ISPIN)=CMPLX(GLOC_RE,GLOC_IM)
        ENDDO
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_READSIGMA(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC &
     &                          ,SIGMA,SIGLAUR,SDC)
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
      COMPLEX(8),INTENT(OUT) :: SIGMA(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(OUT) :: SIGLAUR(NCHI,NCHI,3,NSPIN)
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
      INTEGER(4)             :: NU,ISPIN,I,J,N
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
          SIGMA(:,:,NU,ISPIN)=MATMUL(T,MATMUL(SIGMAPRIME,TPLUS))
        ENDDO
        DO N=1,3
          READ(NFIL,*)((S_RE(I,J),S_IM(I,J),I=1,NCHI),J=1,NCHI)
          SIGMAPRIME=CMPLX(S_RE,S_IM)
          SIGLAUR(:,:,N,ISPIN)=MATMUL(T,MATMUL(SIGMAPRIME,TPLUS))
        ENDDO
!
!       == DOUBLE COUNTING TERM ================================================
        SDC(:,:,ISPIN)=MATMUL(T,MATMUL(SIGMAPRIME,TPLUS))
      ENDDO
      CLOSE(NFIL)
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_WRITESIGMA(NCHI,NOMEGA,NSPIN,KBT,MU &
     &                      ,SIGMA,SIGLAUR,SDC)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NOMEGA
      INTEGER(4),INTENT(IN)  :: NSPIN
      REAL(8)   ,INTENT(IN)  :: KBT
      REAL(8)   ,INTENT(IN)  :: MU
      COMPLEX(8),INTENT(OUT) :: SIGMA(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(OUT) :: SIGLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8),INTENT(OUT) :: SDC(NCHI,NCHI,NSPIN) ! DOUBLE COUNTING
      COMPLEX(8)             :: SIGMAPRIME(NCHI,NCHI)
      INTEGER(4)             :: NU,ISPIN,I,J,N
      INTEGER(4)             :: NFIL=11
!     **************************************************************************
      OPEN(NFIL,FILE=-'SIGMA.DATA')
      WRITE(NFIL,*)NCHI,NOMEGA,NSPIN,KBT,MU !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DO ISPIN=1,NSPIN
!       == SELF ENERGY =========================================================
        DO NU=1,NOMEGA
          SIGMAPRIME=SIGMA(:,:,NU,ISPIN)
          WRITE(NFIL,*)((REAL(SIGMAPRIME(I,J)),AIMAG(SIGMAPRIME(I,J)) &
     &                                                    ,I=1,NCHI),J=1,NCHI)
        ENDDO
!
!       == LAURENT EXPANSION ===================================================
        DO N=1,3
          SIGMAPRIME=SIGLAUR(:,:,N,ISPIN)
          WRITE(NFIL,*)((REAL(SIGMAPRIME(I,J)),AIMAG(SIGMAPRIME(I,J)) &
     &                                                    ,I=1,NCHI),J=1,NCHI)
        ENDDO
!
!       == DOUBLE COUNTING =====================================================
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
!     ** MIMICKS A DMFT SOLVER                                                **
!     **  E=1/BETA * SUM_NU E(IOMEGA_NU0+) TR[ DV * GLOC^\DAGGER]             **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NB,NCHI,NKPTL,NSPIN,NDIM,NOMEGA
      IMPLICIT NONE
      REAL(8)          :: KBT
      REAL(8)          :: MU
      REAL(8)          :: OMEGA(NOMEGA)
      COMPLEX(8)       :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: GLOCLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)       :: DH(NCHI,NCHI,NSPIN)
      COMPLEX(8)       :: SIGMA(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: SIGLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)       :: SDC(NCHI,NCHI,NSPIN)
      REAL(8)          :: EV
      REAL(8)          :: PHILW
      CHARACTER(2)     :: CHOICE
      INTEGER(4)       :: ISPIN,NU
!     **************************************************************************
PRINT*,'ENTERING DMFT$SOLVER...'
      CALL CONSTANTS('EV',EV)
!
!     ==========================================================================
!     == READ GREENS FUNCTION                                                 ==
!     ==========================================================================
PRINT*,'READING GREENS FUNCTION...'
      CALL DMFT_READGLOC(NCHI,NOMEGA,NSPIN,KBT,MU,OMEGA,GLOC,GLOCLAUR)
!
!     ==========================================================================
!     == CONSTRUCT SELF ENERGY                                                ==
!     ==========================================================================
PRINT*,'CONSTRUCTING SELF ENERGY...'
      CHOICE='FP'
      CHOICE='HF'
!
      IF(CHOICE.EQ.'FP') THEN
        DH=(0.D0,0.D0)
        DO ISPIN=1,NSPIN
          DH(1,1,ISPIN)=-1.D0*EV
          DH(2,2,ISPIN)=-1.D0*EV
          DH(3,3,ISPIN)=+1.D0*EV
        ENDDO
        DO ISPIN=1,NSPIN
          DO NU=1,NOMEGA
            SIGMA(:,:,NU,ISPIN)=DH(:,:,ISPIN)
          ENDDO
          SIGLAUR(:,:,1,ISPIN)=DH(:,:,ISPIN)
          SIGLAUR(:,:,2,ISPIN)=0.D0
          SIGLAUR(:,:,3,ISPIN)=0.D0
          SDC(:,:,ISPIN)=0.D0
        ENDDO
!
      ELSE IF(CHOICE.EQ.'HF') THEN
        CALL DMFT$HFSOLVER(NCHI,NSPIN,NOMEGA,KBT,MU,OMEGA &
     &                        ,GLOC,GLOCLAUR,PHILW,SIGMA,SIGLAUR)
!
      ELSE
        CALL ERROR$MSG('CHOICE NOT RECOGNIZED')
        CALL ERROR$STOP('DMFT$SOLVER')
      END IF
!
!     ==========================================================================
!     == WRITE SELF ENERGY
!     ==========================================================================
PRINT*,'WRITING SELF ENERGY   '
      CALL DMFT_WRITESIGMA(NCHI,NOMEGA,NSPIN,KBT,MU,SIGMA,SIGLAUR,SDC) 
PRINT*,'... DMFT$SOLVER DONE.'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$HFSOLVER(NCHI,NSPIN,NOMEGA,KBT,MU,OMEGA &
     &                        ,GLOC,GLOCLAUR,PHILW,SIGMA,SIGLAUR)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NCHI
      INTEGER(4),INTENT(IN)  :: NSPIN
      INTEGER(4),INTENT(IN)  :: NOMEGA
      REAL(8)   ,INTENT(IN)  :: KBT
      REAL(8)   ,INTENT(IN)  :: MU
      REAL(8)   ,INTENT(IN)  :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN)  :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR(NCHI,NCHI,3,NSPIN)
      REAL(8)   ,INTENT(OUT) :: PHILW
      COMPLEX(8),INTENT(OUT) :: SIGMA(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(OUT) :: SIGLAUR(NCHI,NCHI,3,NSPIN)
      REAL(8)   ,PARAMETER   :: UPAR=1.D-2
      REAL(8)   ,PARAMETER   :: JPAR=1.D-3
      REAL(8)                :: U(NCHI,NCHI,NCHI,NCHI)
      COMPLEX(8)             :: RHO(NCHI,NCHI,NSPIN)
      COMPLEX(8)             :: RHO1(NCHI,NCHI)
      COMPLEX(8)             :: H(NCHI,NCHI)
      REAL(8)                :: ESTAT
      INTEGER(4)             :: NU,ISPIN,I,J,K,L
      REAL(8)                :: SVAR
      COMPLEX(8)             :: CSVAR
!     **************************************************************************
!
!     ==========================================================================
!     == SET UP U-TENSOR
!     ==========================================================================
      U=0.D0
      DO I=1,NCHI
        DO J=1,NCHI
          U(I,J,I,J)=UPAR
          U(I,J,J,I)=-JPAR
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == CALCULATE ONE-PARTICLE DENSITY MATRIX                                ==
!     ==========================================================================
      DO ISPIN=1,NSPIN
        RHO(:,:,ISPIN)=(0.D0,0.D0)
        DO NU=1,NOMEGA
          CSVAR=1.D0/CMPLX(0.D0,OMEGA(NU))
          RHO(:,:,ISPIN)=RHO(:,:,ISPIN)+KBT*(GLOC(:,:,NU,ISPIN) &
    &            -CSVAR*(GLOCLAUR(:,:,1,ISPIN) &
    &                   +CSVAR*(GLOCLAUR(:,:,1,ISPIN) &
    &                           +CSVAR*GLOCLAUR(:,:,3,ISPIN))))
        ENDDO
        RHO(:,:,ISPIN)=RHO(:,:,ISPIN)+0.5D0*GLOCLAUR(:,:,1,ISPIN)  &
    &                               -0.25D0*KBT*GLOCLAUR(:,:,2,ISPIN)  &
    &                               -0.125D0*KBT**2*GLOCLAUR(:,:,3,ISPIN)  
!
PRINT*,'IN HFSOLVER',ISPIN
DO I=1,NCHI
  WRITE(*,FMT='("RHO",10("(",2F10.5,")"))')RHO(I,:,ISPIN)
ENDDO
!
      ENDDO
!
!     ==========================================================================
!     == CALCULATE HARTREE FOCK POTENTIAL                                     ==
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
                ESTAT=ESTAT+SVAR*REAL(RHO1(K,I)*RHO1(L,J)-RHO1(K,J)*RHO1(L,I))
                H(I,K)=H(I,K)+SVAR*RHO1(L,J)
                H(J,L)=H(J,L)+SVAR*RHO1(K,I)
                H(J,K)=H(J,K)-SVAR*RHO1(L,I)
                H(I,L)=H(I,L)-SVAR*RHO1(K,J)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO NU=1,NOMEGA
          SIGMA(:,:,NU,ISPIN)=H(:,:)
        ENDDO
        SIGLAUR(:,:,1,ISPIN)=H
        SIGLAUR(:,:,2,ISPIN)=(0.D0,0.D0)
        SIGLAUR(:,:,3,ISPIN)=(0.D0,0.D0)
        IF(NSPIN.EQ.1) ESTAT=2.D0*ESTAT
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_GETSIGMA()
!     **************************************************************************
!     ** PIPSI <PI_A|PSI_N> IS THE PRE-FACTOR OF LOCAL ORBITAL |CHI_A> IN     **
!     **       THE LOCAL ORBITAL EXPANSION OF |PSI_N>                         **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                      ,SIGMA,SIGLAUR,SIGMADC,DIAGSLOC
      IMPLICIT NONE
      INTEGER(4)             :: NU
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT_GETSIGMA')
PRINT*,'ENTERING DMFT_GETSIGMA'
      CALL DMFT_INI()
!
!     ==========================================================================
!     ==  OBTAIN SELF ENERGY
!     ==========================================================================
      CALL DMFT_READSIGMA(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC &
     &                   ,SIGMA,SIGLAUR,SIGMADC)
!
! ATTENTION: FROM HERE ON DELTA-SIGMA AND DC-SIGMA IS STORED
!            READSIGMA PROBABLY RETURNS THE FULL SIGMA AND DC-SIGMA
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

PRINT*,'... DMFT_GETSIGMA DONE'
                                       CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_HRHO()
!     **************************************************************************
!     ** CONSTRUCT HAMILTONIAN IN THE SPACE OF CORRELATED ORBITALS            **
!     ** THAT PRODUCES THE CORRECT DENSITY MATRIX
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NSPIN,NOMEGA,HRHO
      IMPLICIT NONE
      COMPLEX(8)             :: SIG(NCHI,NCHI,NOMEGA,NSPIN)
!     **************************************************************************
      IF(.NOT.TON) RETURN
      SIG=(0.D0,0.D0)
      CALL DMFT_CONSTRAINTS_TWO(HRHO,SIG)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_CONSTRAINTS()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NSPIN,NDIM,NOMEGA,KBT,MU &
     &                      ,SIGMA,ERHO,H0,HRHO
      IMPLICIT NONE
      INTEGER(4)             :: IKPT,ISPIN,I
      COMPLEX(8)             :: SIG(NCHI,NCHI,NOMEGA,NSPIN)
      LOGICAL(4),PARAMETER   :: TPRINT=.TRUE.
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT_GETSIGMA')
!
!     ==========================================================================
!     ==  CONSTRUCT ONE-PARTICLE HAMILTONIANS THAT OBEY                       ==
!     ==  DENSITY MATRIX CONSTRAINT                                           ==
!     ==========================================================================
      CALL DMFT_CONSTRAINTS_TWO(H0,SIGMA)
!
!     ==========================================================================
!     == MAP ONTO HAMILTONIAN
!     ==========================================================================
      IF(TPRINT) THEN
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            WRITE(*,FMT='(82("="),T10,"IKPT=",I5," ISPIN=",I2)')IKPT,ISPIN
            WRITE(*,FMT='("ERHO ",10F10.5)')ERHO(:,IKPT,ISPIN)
            DO I=1,NCHI
              WRITE(*,FMT='("DEDRHO",80("(",2F10.5,")"))') &
    &                 HRHO(I,:,IKPT,ISPIN)-H0(I,:,IKPT,ISPIN)
            ENDDO
            DO I=1,NCHI
              WRITE(*,FMT='("HRHO  ",80("(",2F10.5,")"))')HRHO(I,:,IKPT,ISPIN)
            ENDDO
            DO I=1,NCHI
              WRITE(*,FMT='("H0   =",80("(",2F10.5,")"))')H0(I,:,IKPT,ISPIN)
            ENDDO
          ENDDO
        ENDDO
      END IF

                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_CONSTRAINTS_TWO(H0,SIG)
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                      ,PIPSI,ERHO,AMIX
      IMPLICIT NONE
      COMPLEX(8),INTENT(INOUT) :: H0(NCHI,NCHI,NKPTL,NSPIN)
      COMPLEX(8),INTENT(IN)    :: SIG(NCHI,NCHI,NOMEGA,NSPIN)
      REAL(8)   ,PARAMETER     :: TOL=1.D-6
      INTEGER(4),PARAMETER     :: NITER=1000
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
      SUBROUTINE DMFT_GLOC(H0,SIGMA,SIGLAUR,GLOC,GLAUR)
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                       ,WKPTL,SMAT,SINV
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN)    :: H0(NCHI,NCHI,NKPTL,NSPIN)
      COMPLEX(8),INTENT(IN)    :: SIGMA(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(IN)    :: SIGLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8),INTENT(OUT)   :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(OUT)   :: GLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)               :: GINV(NCHI,NCHI)
      COMPLEX(8)               :: MAT(NCHI,NCHI)
      COMPLEX(8)               :: G(NCHI,NCHI)
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: CSVAR
      INTEGER(4)               :: IKPT,ISPIN,NU,I,J
      INTEGER(4)               :: NFIL
      LOGICAL(4),PARAMETER     :: TPRINT=.FALSE.
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT_GETSIGMA')
      DO ISPIN=1,NSPIN
!       ========================================================================
!       == LAURENT EXPANSION OF THE LOCAL GREENS FUNCTION FROM 0 TO 2
!       == REMEMBER THAT MU=0 !!!
!       ========================================================================
        GLAUR(:,:,1,ISPIN)=(0.D0,0.D0)
        GLAUR(:,:,2,ISPIN)=(0.D0,0.D0)
        GLAUR(:,:,3,ISPIN)=(0.D0,0.D0)
        DO IKPT=1,NKPTL
          MAT=SINV(:,:,IKPT,ISPIN)
          GLAUR(:,:,1,ISPIN)=GLAUR(:,:,1,ISPIN)+WKPTL(IKPT)*MAT
          MAT=H0(:,:,IKPT,ISPIN)+SIGLAUR(:,:,1,ISPIN)-MU*SMAT(:,:,IKPT,ISPIN)
          MAT=MATMUL(SINV(:,:,IKPT,ISPIN),MATMUL(MAT,SINV(:,:,IKPT,ISPIN)))
          GLAUR(:,:,2,ISPIN)=GLAUR(:,:,2,ISPIN)+WKPTL(IKPT)*MAT
          MAT=H0(:,:,IKPT,ISPIN)+SIGLAUR(:,:,1,ISPIN)-MU*SMAT(:,:,IKPT,ISPIN)
          MAT=MATMUL(MAT,MATMUL(SINV(:,:,IKPT,ISPIN),MAT))
          MAT=MAT+SIGLAUR(:,:,2,ISPIN)
          MAT=MATMUL(SINV(:,:,IKPT,ISPIN),MATMUL(MAT,SINV(:,:,IKPT,ISPIN)))
          GLAUR(:,:,3,ISPIN)=GLAUR(:,:,3,ISPIN)+WKPTL(IKPT)*MAT
        ENDDO
!
!       ========================================================================
!       == LOCAL GREENS FUNCTION ON THE MATSUBARA GRID
!       ========================================================================
        GLOC(:,:,:,ISPIN)=(0.D0,0.D0)
        DO IKPT=1,NKPTL
          DO NU=1,NOMEGA
            CSVAR=CI*OMEGA(NU)+MU
            GINV(:,:)=SMAT(:,:,IKPT,ISPIN)*CSVAR &
     &               -H0(:,:,IKPT,ISPIN)-SIGMA(:,:,NU,ISPIN)
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
               WRITE(*,FMT='(10("(",2F10.5,")"))')GLAUR(J,:,I,ISPIN)
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
      USE DMFT_MODULE, ONLY  : TON,NKPTL,NSPIN,NB
      IMPLICIT NONE
      REAL(8)        :: EIG(NB,NKPTL,NSPIN)
      INTEGER(4)     :: IKPT,ISPIN,IB
!     **************************************************************************
      IF(.NOT.TON) RETURN
CALL ERROR$MSG('THIS ROUTINE IS CORRUPTED')
CALL ERROR$STOP('DMFT$ADDTOEIGVAL')
      CALL WAVES_DYNOCCGETR8A('EPSILON',NB*NKPTL*NSPIN,EIG)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          DO IB=1,NB
!            EIG(IB,IKPT,ISPIN)=EIG(IB,IKPT,ISPIN) &
!     &                        +REAL(GAMMA0(IB,IB,IKPT,ISPIN))
          ENDDO
        ENDDO
      ENDDO
OPEN(11,FILE='EIGS.DAT')
DO IKPT=1,NKPTL
  DO IB=1,NB
!    WRITE(11,*)EIG(IB,IKPT,1)*27.211D0,REAL(GAMMA0(IB,IB,IKPT,1))*27.211D0 &
!  &                                   ,AIMAG(GAMMA0(IB,IB,IKPT,1))*27.211D0
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
      USE DMFT_MODULE, ONLY  : TON,NKPTL,NSPIN,NB,NDIM
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
CALL ERROR$MSG('THIS ROUTINE IS CORRUPTED')
CALL ERROR$STOP('DMFT$ADDTOHPSI')
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
!          AMAT=MATMUL(DENMAT(:,:,IKPT,ISPIN),GAMMA0(:,:,IKPT,ISPIN))
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
!     ..........................................................................
      SUBROUTINE DMFT$PLOTSIGMA(FILE)
      USE DMFT_MODULE, ONLY : NCHI,NSPIN,NOMEGA,OMEGA,SIGMA,SIGLAUR
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)              :: NU
      INTEGER(4)  ,PARAMETER  :: NFIL=11
      COMPLEX(8)              :: CARG,CMAT(NCHI,NCHI,NSPIN)
!     **************************************************************************
      OPEN(NFIL,FILE=FILE)
      REWIND NFIL
      DO NU=2,NOMEGA
        CARG=1.D0/CMPLX(0.D0,OMEGA(NU))
        CMAT=SIGLAUR(:,:,1,:) &
     &      +CARG*(SIGLAUR(:,:,2,:) &
     &            +CARG*(SIGLAUR(:,:,3,:)))
        WRITE(NFIL,*)OMEGA(NU),REAL(SIGMA(:,:,NU,:)),AIMAG(SIGMA(:,:,NU,:)) &
    &               ,REAL(CMAT),AIMAG(CMAT)
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE DMFT$PLOTGLOC(FILE)
      USE DMFT_MODULE, ONLY : NOMEGA,OMEGA,NCHI,NSPIN,GLOC,GLOCLAUR
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)              :: NU
      INTEGER(4)  ,PARAMETER  :: NFIL=11
      COMPLEX(8)              :: CARG,CMAT(NCHI,NCHI,NSPIN)
!     **************************************************************************
      OPEN(NFIL,FILE=FILE)
      REWIND NFIL
      DO NU=2,NOMEGA
        CARG=1.D0/CMPLX(0.D0,OMEGA(NU))
        CMAT=CARG*(GLOCLAUR(:,:,1,:) &
     &              +CARG*(GLOCLAUR(:,:,2,:) &
     &                     +CARG*(GLOCLAUR(:,:,3,:))))
        WRITE(NFIL,*)OMEGA(NU),REAL(GLOC(:,:,NU,:)),AIMAG(GLOC(:,:,NU,:)) &
    &               ,REAL(CMAT),AIMAG(CMAT)
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
