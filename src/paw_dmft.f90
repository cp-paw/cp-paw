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
REAL(8)   ,PARAMETER   :: AMIX=1.D-1
INTEGER(4)             :: NOMEGA
INTEGER(4)             :: NCHI          ! #(CORRELATED ORBITALS)
INTEGER(4)             :: NB            ! #(BAND STATES PER K-POINT)
INTEGER(4)             :: NKPTL         ! #(KPOINTS ON THIS TASK)
INTEGER(4)             :: NSPIN         ! #(SPIN COMPONENTS)
INTEGER(4)             :: NDIM          ! #(SPINOR COMPONENTS)
INTEGER(4)             :: Nat           ! #(atoms)
REAL(8)   ,ALLOCATABLE :: OMEGA(:)      ! MATSUBARA FREQUENCIES
INTEGER(4),ALLOCATABLE :: IPROOFCHI(:)  !(nchi) map ichi to ipro
INTEGER(4),ALLOCATABLE :: Ichibnd(:,:)  !(2,nat) map iat to ichistart,ichiend
REAL(8)                :: KBT           ! TEMPERATURE (K_B*T)
REAL(8)                :: MU            ! CHEMICAL POTENTIAL
REAL(8)                :: DELTAT        ! TIMESTEP
REAL(8)   ,ALLOCATABLE :: WKPTL(:)      ! (NKPTL) K-POINT WEIGHTS ON THIS TASK
COMPLEX(8),ALLOCATABLE :: PIPSI(:,:,:,:)    !(NCHI,NB,NKPTL,NSPIN) <PI|PSI>
REAL(8)   ,ALLOCATABLE :: ERHO(:,:,:)       !(NB,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: H0(:,:,:,:)       !(NCHI,NCHI,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: HRHO(:,:,:,:)     !(NCHI,NCHI,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: RHOOFK(:,:,:,:)   !(NCHI,NCHI,NKPTL,NSPIN)
!__MATRIX TO DIAGONALIZE OVERLAP MATRIX_________________________________________
COMPLEX(8),ALLOCATABLE :: DIAGSLOC(:,:,:) !(NCHI,NCHI,NSPIN)
COMPLEX(8),ALLOCATABLE :: SMAT(:,:,:,:)   !(NCHI,NCHI,NKPTL,NSPIN)
COMPLEX(8),ALLOCATABLE :: SINV(:,:,:,:)   !(NCHI,NCHI,NKPTL,NSPIN)
!== LOCAL GREENS FUNCTION (WITH LAURENT EXPANSION COEFFICIENTS) ================
COMPLEX(8),ALLOCATABLE :: GLOC(:,:,:,:)     !(NCHI,NCHI,NOMEGA,NSPIN)
COMPLEX(8),ALLOCATABLE :: GLOCLAUR(:,:,:,:) !(NCHI,NCHI,3,NSPIN)
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
      USE WAVES_MODULE, ONLY : KMAP,NDIM_W=>NDIM,NKPTL_W=>NKPTL,NSPIN_W=>NSPIN
      IMPLICIT NONE
      REAL(8)                :: PI
      INTEGER(4)             :: ISP0
      INTEGER(4)             :: L0
      REAL(8)   ,ALLOCATABLE :: WKPT(:)    ! K-POINT WEIGHTS (GLOBAL)
      INTEGER(4)             :: NTASKS_K,THISTASK_K
      INTEGER(4)             :: NTASKS_M,THISTASK_M
      INTEGER(4)             :: NKPT
      INTEGER(4)             :: nsp
      INTEGER(4)             :: lnx
      INTEGER(4)             :: l
      INTEGER(4)             :: NU,ISP,iat,LN,im,IKPTL,IKPT,ichi,ipro
      integer(4),allocatable :: lox(:) !(lnx)
      logical(4),allocatable :: torb(:) !(lnx)
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
      call atomlist$natom(nat)
      allocate(ichibnd(2,nat))
      CALL SETUP$GETI4('NSP',NSP) 
!     == get nchi
      nchi=0
      ichibnd(1,:)=1
      ichibnd(2,:)=0
      do iat=1,nat
        CALL ATOMLIST$GETI4('ISPECIES',iat,isp)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(TORB(LNx))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETL4A('TORB',LNX,TORB)
print*,'torb ',iat,torb
        ichibnd(1,iat)=nchi+1
        DO LN=1,LNX
          l=lox(ln)
          if(torb(ln)) then
            nchi=nchi+2*l+1
          end if
        enddo
        ichibnd(2,iat)=nchi
        deallocate(lox)
        deallocate(torb)
        CALL SETUP$unSELECT()
      ENDDO
print*,'nchi ',nchi
!
!     == accumulate iproofchi ==================================================
      ALLOCATE(IPROOFCHI(NCHI))
      ichi=0
      ipro=0
      do iat=1,nat
        CALL ATOMLIST$GETI4('ISPECIES',iat,isp)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(TORB(LNx))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETL4A('TORB',LNX,TORB)
        DO LN=1,LNX
          l=lox(ln)
          if(torb(ln)) then
            do im=1,2*l+1
              ichi=ichi+1
              ipro=ipro+1
              iproofchi(ichi)=ipro
            enddo
          else
            ipro=ipro+2*l+1
          end if
        enddo
        deallocate(lox)
        deallocate(torb)
        CALL SETUP$unSELECT()
      ENDDO
print*,'iproofchi ',iproofchi
!
!!$!THIS IS SPECIFIC FOR SRVO3
!!$      ISP0=2
!!$      L0=2
!!$      NCHI=3    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$      IPRO=0
!!$      DO ISP=1,ISP0-1
!!$        IPRO=IPRO+SUM(2*LOX(:LNX(ISP),ISP)+1)
!!$      ENDDO
!!$      DO LN=1,LNX(ISP0)
!!$        IF(LOX(LN,ISP0).EQ.L0) EXIT
!!$        IPRO=IPRO+2*LOX(LN,ISP0)+1
!!$      ENDDO
!!$      ALLOCATE(IPROOFCHI(NCHI))
!!$      IPROOFCHI(1)=IPRO+2   !T2G ORBITALS
!!$      IPROOFCHI(2)=IPRO+4
!!$      IPROOFCHI(3)=IPRO+5
!
!     ==========================================================================
!     == OTHER VARIABLES                                                      ==
!     ==========================================================================
!THIS IS HARD WIRED
      NOMEGA=400
      KBT=0.333D0*EV
      DELTAT=10.D0
      MU=0.D0
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
      ALLOCATE(RHOOFK(NCHI,NCHI,NKPTL,NSPIN))
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
      INTEGER(4)             :: I,ITER,IKPT,ISPIN
      LOGICAL(4)             :: TPRINT=.FALSE.
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
!     ==  TRANSFORMATION ONTO A ORTHONORMAL CORRELATED BASIS SET              ==
!     ==    |CHIORTHO>   =|CHINONORTHO>*DIAGSLOC                              ==
!     ==========================================================================
      CALL DMFT_DIAGSLOC()
!
      CALL DMFT_GRHO() ! TEST ONLY: DENSITY FROM REFERENCE GREENS FUNCTION
!
!     ==========================================================================
!     ==  CONSTRUCT NON-INTERACTING HAMILTONIAN THAT PRODUCES THE CORRECT     ==
!     ==  ONE-PARTICLE DENSITY MATRIX                                         ==
!     ==========================================================================
      HRHO=(0.D0,0.D0)
      SIGMA=(0.D0,0.D0)
      SIGLAUR=(0.D0,0.D0)
      CALL DMFT_CONSTRAINTS(HRHO,SIGMA,SIGLAUR)
!
      IF(.NOT.ALLOCATED(H0)) THEN
        ALLOCATE(H0(NCHI,NCHI,NKPTL,NSPIN))
        H0=HRHO
      END IF
!
!      CALL DMFT_TEST()  !TEST ONLY
!
!     ==========================================================================
!     == DETERMINE LOCAL GREENS FUNCTION                                      ==
!     ==========================================================================
      MU=0.D0
DO ITER=1,10
WRITE(*,FMT='(82("="),T20," ITERATION ",I5)')ITER
      CALL DMFT_GLOC(H0,SIGMA,SIGLAUR,GLOC,GLOCLAUR)
!      CALL DMFT$PLOTGLOC(-'GLOC.DAT')
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
!      CALL DMFT$PLOTSIGMA(-'SIG.DAT')
!
!     ==========================================================================
!     ==  CONSTRAINTS
!     ==========================================================================
print*,'marke 13'
      CALL DMFT_CONSTRAINTS(H0,SIGMA,SIGLAUR)
print*,'marke 14'
!
!     ==========================================================================
!     == MAP ONTO HAMILTONIAN
!     ==========================================================================
      IF(TPRINT) THEN
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            WRITE(*,FMT='(82("="),T10,"IKPT=",I5," ISPIN=",I2)')IKPT,ISPIN
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
print*,'marke 15'
print*,'iteration completed ',iter
ENDDO
      call DMFT$ADDTOHPSI()

!STOP 'END OF LOOP. STOPPING.'
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
     &                       ,ERHO,PIPSI,RHOOFK,WKPTL,ichibnd,nat
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: NGL
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IPRO,IB,I,J,iat
      REAL(8)                :: F(NB,NKPTL,NSPIN)
      REAL(8)                :: SVAR
      COMPLEX(8)             :: RHO(NCHI,NCHI,NSPIN)
      COMPLEX(8)             :: CSVAR
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                      CALL TRACE$PUSH('DMFT$COLLECTHAMILTONIAN')
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
!
!     ==========================================================================
!     ==  DETERMINE ERHO 
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
!     ==  DETERMINE CONSTRAINT RHOOFK 
!     ==========================================================================
      RHOOFK=(0.D0,0.D0)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          DO IB=1,NB
            DO J=1,NCHI
              CSVAR=F(IB,IKPT,ISPIN)*CONJG(PIPSI(J,IB,IKPT,ISPIN))
              RHOOFK(:,J,IKPT,ISPIN)=RHOOFK(:,J,IKPT,ISPIN) &
     &                               +PIPSI(:,IB,IKPT,ISPIN)*CSVAR
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE LOCAL DENSITY MATRIX FOR TESTING 
!     ==========================================================================
      IF(TTEST) THEN
!       == CALCULATE DENSITY MATRIX FOR TESTING ================================
        RHO(:,:,:)=(0.D0,0.D0)
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            RHO(:,:,ISPIN)=RHO(:,:,ISPIN)+RHOOFK(:,:,IKPT,ISPIN)*WKPTL(IKPT)
          ENDDO
        ENDDO
        IF(NSPIN.EQ.1)RHO=2.D0*RHO
        DO IAT=1,NAT
          IF(ICHIBND(2,IAT).LT.ICHIBND(1,IAT)) CYCLE
          WRITE(*,FMT='(82("="),T10," DENSITY MATRIX (1) FOR ATOM ",I5," FROM BAND OCCUPATIONS ")')IAT
          DO ISPIN=1,NSPIN
            DO I=ICHIBND(1,IAT),ICHIBND(2,IAT)
              WRITE(*,FMT='("RHO(IS=",I1,"):",100("(",2F10.5,")"))') &
    &                      ISPIN,RHO(I,ICHIBND(1,IAT):ICHIBND(2,IAT),ISPIN)
            ENDDO
          ENDDO
        ENDDO
      END IF
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
      SUBROUTINE DMFT_READGLOC(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC,OMEGA &
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
      COMPLEX(8),INTENT(IN)  :: DIAGSLOC(NCHI,NCHI,NSPIN)
      REAL(8)   ,INTENT(OUT) :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(OUT) :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(OUT) :: GLOCLAUR(NCHI,NCHI,3,NSPIN)
      REAL(8)                :: GLOC_RE(NCHI,NCHI),GLOC_IM(NCHI,NCHI)
      COMPLEX(8)             :: T(NCHI,NCHI)
      COMPLEX(8)             :: TPLUS(NCHI,NCHI)
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
!
!     ==========================================================================
!     ==  TRANSFORM BACK                                                      ==
!     ==========================================================================
      DO ISPIN=1,NSPIN
        CALL LIB$INVERTC8(NCHI,DIAGSLOC(:,:,ISPIN),TPLUS)
        T=TRANSPOSE(CONJG(TPLUS)) 
        DO NU=1,NOMEGA
          GLOC(:,:,NU,ISPIN)=MATMUL(T,MATMUL(GLOC(:,:,NU,ISPIN),TPLUS))
        ENDDO
        DO NU=1,3
          GLOCLAUR(:,:,NU,ISPIN)=MATMUL(T,MATMUL(GLOCLAUR(:,:,NU,ISPIN),TPLUS))
        ENDDO
     ENDDO
 
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
      USE DMFT_MODULE, ONLY: NB,NCHI,NKPTL,NSPIN,NDIM,NOMEGA,GLOC,GLOCLAUR,DIAGSLOC
      IMPLICIT NONE
      REAL(8)          :: KBT
      REAL(8)          :: MU
      REAL(8)          :: OMEGA(NOMEGA)
      COMPLEX(8)       :: GLOC1(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: GLOCLAUR1(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)       :: DH(NCHI,NCHI,NSPIN)
      COMPLEX(8)       :: SIGMA(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: SIGLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)       :: SDC(NCHI,NCHI,NSPIN)
      REAL(8)          :: EV
      REAL(8)          :: PHILW
      REAL(8)          :: svar
      CHARACTER(2)     :: CHOICE
      INTEGER(4)       :: ISPIN,NU,ichi
!     **************************************************************************
PRINT*,'ENTERING DMFT$SOLVER...'
      CALL CONSTANTS('EV',EV)
!
!     ==========================================================================
!     == READ GREENS FUNCTION                                                 ==
!     ==========================================================================
PRINT*,'READING GREENS FUNCTION...'
      CALL DMFT_READGLOC(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC,OMEGA &
     &        ,GLOC,GLOCLAUR)
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
          svar=ev
          do ichi=1,nchi
            DH(ichi,ichi,ISPIN)=svar
            svar=-svar
          enddo
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
      SUBROUTINE DMFT_ULOCAL(IAT,LMNX,U)
!     **************************************************************************
!     **  FOLD DOWN THE U-TENSOR FOR THE TAILED VALENCE AND SCATTERING WAVES  **
!     **  ONTO THE TAILED LOCAL ORBITALS FOR THE ACTUAL STRUCTURE CONSTANTS   **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR,SBAR,sbarli1,LNX,LOX
      INTEGER(4),INTENT(IN)  :: IAT
      INTEGER(4),INTENT(IN)  :: LMNX          !  #(LOCAL ORBITALS ON THIS SITE)
      REAL(8)   ,INTENT(OUT) :: U(LMNX,LMNX,LMNX,LMNX)
      logical(4),parameter   :: ttest=.true.
      REAL(8)   ,ALLOCATABLE :: UT(:,:,:,:)
      INTEGER(4)             :: ISP     ! ATOM TYPE
      INTEGER(4)             :: LMNXT   ! #(VALENCE+SCATTERING WAVES)
      INTEGER(4)             :: LMX     ! #(SCATTERING WAVES)
      INTEGER(4)             :: NNS     ! #(PAIRS FOR STRUCTURE CONSTANTS)
      REAL(8)                :: SVAR
      INTEGER(4)             :: LMN,L,LM,LN,IM,NN,NN0
!     **************************************************************************
      ISP=ISPECIES(IAT)
      LMNXT=POTPAR(ISP)%TAILED%LMNX  ! SIZE OF U-TENSOR ON POTPAR
!
!     ==========================================================================
!     == COLLECT U-TENSOR IN EXPANDED BASIS                                   ==
!     ==========================================================================
      ALLOCATE(Ut(LMNXT,LMNXT,LMNXT,LMNXT))
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
            UT(:,:,LMN,:lmnx)=UT(:,:,LMN,:lmnx)-SVAR*UT(:,:,LMNX+J,:lmnx)
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
            UT(:,LMN,:lmnx,:lmnx)=UT(:,LMN,:lmnx,:lmnx) &
      &                         -SVAR*UT(:,LMNX+J,:lmnx,:lmnx)
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
            UT(lmn,:LMNX,:LMNX,:LMNX)=UT(lmn,:LMNX,:LMNX,:LMNX) &
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
      subroutine dmft$EX(iat,lmnx,ndimd,d,e,h)
      use lmto_module, only : potpar
      integer(4),intent(in)  :: iat
      integer(4),intent(in)  :: lmnx
      integer(4),intent(in)  :: ndimd
      real(8)   ,intent(in)  :: d(lmnx,lmnx,ndimd)
      real(8)   ,intent(out) :: e
      real(8)   ,intent(out) :: h(lmnx,lmnx,ndimd)
      integer(4)             :: lnxt
      integer(4)             :: lmnxt
      REAL(8)   ,ALLOCATABLE :: DT(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HT(:,:,:)
      REAL(8)   ,ALLOCATABLE :: U(:,:,:,:)
      real(8)                :: svar
      real(8)                :: eh,ex
      integer(4)             :: i,j,k,l
!     **************************************************************************
      LNXT=POTPAR(ISP)%TAILED%LNX
      LMNXT=POTPAR(ISP)%TAILED%LMNX
      ALLOCATE(DT(LMNXT,LMNXT,NDIMD))
      ALLOCATE(HT(LMNXT,LMNXT,NDIMD))
      CALL LMTO_BLOWUPDENMATNL(IAT,IAT,NDIMD,LMNX,LMNX,D,LMNXT,LMNXT,DT)
!
!     ========================================================================
!     == CALCULATE U-TENSOR                                                 ==
!     ========================================================================
      ALLOCATE(U(LMNXT,LMNXT,LMNXT,LMNXT))
      U=POTPAR(ISP)%TAILED%U
!
!     ========================================================================
!     == WORK OUT TOTAL ENERGY AND HAMILTONIAN                              ==
!     ========================================================================
      EH=0.D0
      EX=0.D0
      HT(:,:,:)=0.D0
      DO I=1,LMNXT
        DO J=1,LMNXT
          DO K=1,LMNXT
            DO L=1,LMNXT
!             ================================================================
!             == HARTREE TERM (NOT CONSIDERED)                              ==
!             == EH=EH+0.5D0*U(I,J,K,L)*D(K,I,1)*D(L,J,1)                   ==
!             == H(K,I,1)=H(K,I,1)+0.5D0*U(I,J,K,L)*D(L,J,1) !DE/DRHO(K,I,1)==
!             == H(L,J,1)=H(L,J,1)+0.5D0*U(I,J,K,L)*D(K,I,1) !DE/DRHO(L,J,1)==
!             ================================================================
!             ================================================================
!             == EXCHANGE ENERGY =============================================
!             ================================================================
!             == AN ADDITIONAL FACTOR COMES FROM THE REPRESENTATION INTO TOTAL AND SPIN
              SVAR=-0.25D0*U(I,J,K,L)
              EX=EX+SVAR*SUM(DT(K,J,:)*DT(I,L,:))
              HT(K,J,:)=HT(K,J,:)+SVAR*DT(I,L,:) 
              HT(I,L,:)=HT(I,L,:)+SVAR*DT(K,J,:) 
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,NDIMD,LMNXT,LMNXT,HT,LMNX,LMNX,H)
      e=ex
      return
      end
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
      CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,GLOC,GLOCLAUR,RHO)
!
PRINT*,'RHO FROM GLOC IN HFSOLVER',ISPIN
DO ISPIN=1,NSPIN
  DO I=1,NCHI
    WRITE(*,FMT='("RHO",10("(",2F10.5,")"))')RHO(I,:,ISPIN)
  ENDDO
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
      SUBROUTINE DMFT_CONSTRAINTS(H0,SIG,SIGLAUR)
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
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                      ,AMIX,SMAT,HRHO,RHOOFK
      IMPLICIT NONE
      COMPLEX(8),INTENT(INOUT) :: H0(NCHI,NCHI,NKPTL,NSPIN)
      COMPLEX(8),INTENT(IN)    :: SIG(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(IN)    :: SIGLAUR(NCHI,NCHI,3,NSPIN)
      REAL(8)   ,PARAMETER     :: TOL=1.D-6
      INTEGER(4),PARAMETER     :: NITER=1000
      LOGICAL(4),PARAMETER     :: TMIX=.FALSE.
      LOGICAL(4)               :: CONVG
      REAL(8)                  :: MAXDEV
      INTEGER(4)               :: NU
      COMPLEX(8)               :: DEVRHO(NCHI,NCHI)
      INTEGER(4)               :: NTASKS_K,THISTASK_K
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: CSVAR
      COMPLEX(8)               :: X4(NCHI*NCHI,NCHI*NCHI)
      COMPLEX(8)               :: MAT(NCHI,NCHI)
      COMPLEX(8)               :: SINV(NCHI,NCHI)
      COMPLEX(8)               :: G(NCHI,NCHI)
      COMPLEX(8)               :: GLAUR1(NCHI,NCHI)
      COMPLEX(8)               :: GLAUR2(NCHI,NCHI)
      COMPLEX(8)               :: GLAUR3(NCHI,NCHI)
      COMPLEX(8)               :: DH0(NCHI,NCHI)
      COMPLEX(8)               :: X4INV(NCHI*NCHI,NCHI*NCHI)
      INTEGER(4)               :: IKPT,ISPIN,IB,ITER,I,J,K,L,IND1,IND2
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_GETSIGMA')
print*,'h0(up) ',maxval(abs(h0(:,:,:,1)))
print*,'h0(dn) ',maxval(abs(h0(:,:,:,2)))
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
          DO ISPIN=1,NSPIN
!           == LAURENT EXPANSION FOR THE GREENS FUNCTION =======================
            CALL LIB$INVERTC8(NCHI,SMAT(:,:,IKPT,ISPIN),SINV)
            GLAUR1=SINV
            MAT=H0(:,:,IKPT,ISPIN)+SIGLAUR(:,:,1,ISPIN)-MU*SMAT(:,:,IKPT,ISPIN)
            GLAUR2=MATMUL(SINV,MATMUL(MAT,SINV))
            MAT=H0(:,:,IKPT,ISPIN)+SIGLAUR(:,:,1,ISPIN)-MU*SMAT(:,:,IKPT,ISPIN)
            MAT=MATMUL(MAT,MATMUL(SINV,MAT))
            MAT=MAT+SIGLAUR(:,:,2,ISPIN)
            GLAUR3=MATMUL(SINV,MATMUL(MAT,SINV))
!
            X4=(0.D0,0.D0)        
            DEVRHO=(0.D0,0.D0)
            DO NU=1,NOMEGA
        
!             == CONSTRUCT LATTICE GREENS FUNCTION =============================
              MAT=(CI*OMEGA(NU)+MU)*SMAT(:,:,IKPT,ISPIN) &
                                   -H0(:,:,IKPT,ISPIN)-SIG(:,:,NU,ISPIN)
              CALL LIB$INVERTC8(NCHI,MAT,G)
              CSVAR=1.D0/(CI*OMEGA(NU))
              DEVRHO=DEVRHO+KBT*(G-CSVAR*(GLAUR1+CSVAR*(GLAUR2+CSVAR*GLAUR3)))
!
!             == ACCUMULATE SECOND ORDER TERM ==================================
              DO I=1,NCHI
                DO J=1,NCHI
                  IND1=I+(J-1)*NCHI  !(I,J)
                  DO K=1,NCHI
                    DO L=1,NCHI
                      IND2=K+(L-1)*NCHI  !(K,L)
                      X4(IND1,IND2)=X4(IND1,IND2)+KBT*(G(I,K)*G(L,J) &
    &                                           +CONJG(G(K,I)*G(J,L)))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
!
!           == INCLUDE NEGATIVE FREQUENCIES ====================================
!           == ALREADY INCLUDED FOR X4
            DEVRHO=DEVRHO+CONJG(TRANSPOSE(DEVRHO))
!           == ADD TAILS (GLAUR3 DOES NOT CONTRIBUTE) ==========================
            DEVRHO=DEVRHO+0.5D0*GLAUR1
            DEVRHO=DEVRHO-0.25D0/KBT*GLAUR2
!
!           == SUBTRACT TARGET DENSITY =========================================
            DEVRHO=DEVRHO-RHOOFK(:,:,IKPT,ISPIN)
!
            MAXDEV=MAX(MAXDEV,MAXVAL(ABS(DEVRHO)))

PRINT*,'ITER=',ITER,' IKPT=',IKPT,' ispin=',ispin &
&     ,'MAXVAL(X4) ',MAXVAL(ABS(X4)),' MAX(DEVRHO) ',MAXVAL(ABS(DEVRHO))
!

!
            IF(TMIX) THEN
              DH0(:,:)=AMIX*DEVRHO(:,:)
            ELSE
              CALL LIB$INVERTC8(NCHI**2,X4,X4INV)
              DH0=(0.D0,0.D0)
              DO I=1,NCHI
                DO J=1,NCHI
                  IND1=I+NCHI*(J-1)
                  DO K=1,NCHI
                    DO L=1,NCHI
                      IND2=K+NCHI*(L-1)
                      DH0(I,J)=DH0(I,J)-X4INV(IND1,IND2)*DEVRHO(K,L)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            END IF
!
!           =================================================================== 
!           == ADD DH
!           ====================================================================
            H0(:,:,IKPT,ISPIN)=H0(:,:,IKPT,ISPIN)+DH0(:,:)
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
        CALL ERROR$STOP('DMFT_CONSTRAINTS')
      END IF

                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,G,GLAUR,RHO)
!     **************************************************************************
!     **  CALCULATES THE DENSITY MATRIX FROM A GREENS FUNCTION                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: NSPIN
      REAL(8)   ,INTENT(IN) :: KBT
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN) :: G(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(IN) :: GLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8),INTENT(OUT):: RHO(NCHI,NCHI,NSPIN)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)  ! SQRT(-1)
      INTEGER(4)            :: ISPIN,NU
      COMPLEX(8)            :: CSVAR
!     **************************************************************************
      DO ISPIN=1,NSPIN
        RHO(:,:,ISPIN)=(0.D0,0.D0)
        DO NU=1,NOMEGA
          CSVAR=1.D0/(CI*OMEGA(NU))
          RHO(:,:,ISPIN)=RHO(:,:,ISPIN)+G(:,:,NU,ISPIN) &
     &                                 -GLAUR(:,:,1,ISPIN)*CSVAR &
     &                                 -GLAUR(:,:,2,ISPIN)*CSVAR**2 &
     &                                 -GLAUR(:,:,3,ISPIN)*CSVAR**3
        ENDDO
        RHO(:,:,ISPIN)=RHO(:,:,ISPIN)+CONJG(TRANSPOSE(RHO(:,:,ISPIN)))
        RHO(:,:,ISPIN)=RHO(:,:,ISPIN)*KBT
        RHO(:,:,ISPIN)=RHO(:,:,ISPIN)+0.5D0*GLAUR(:,:,1,ISPIN)
        RHO(:,:,ISPIN)=RHO(:,:,ISPIN)-0.25D0/KBT*GLAUR(:,:,2,ISPIN)
        IF(NSPIN.EQ.1) RHO(:,:,ISPIN)=RHO(:,:,ISPIN)*2.D0 ! SPIN MULTIPLICITY
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_GRHO()
!     **************************************************************************
!     **  CALCULATES GRHO AND FROM IT THE DENSITY MATRIX FOR TESTING PURPOSES **
!     **  REMEMBER THAT WE NEED THE K-DEPENDENT GRHO, WHILE THIS ROUTINE ONLY **
!     **  CALCULATES THE K-SUMMED GRHO                                        **
!     **                                                                      **
!     **  GRHO(K,NU)=<PI(K)|PSI(K)>/[I*OMEGA_NU+MU-HRHO(K)]<PI(K)|PSI(K)>     **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                      ,PIPSI,ERHO,WKPTL
      IMPLICIT NONE
      LOGICAL(4),PARAMETER     :: TTEST=.TRUE.
      INTEGER(4)               :: NU
      COMPLEX(8)               :: GRHO(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)               :: GRHOLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)               :: RHO(NCHI,NCHI,NSPIN)
      COMPLEX(8)               :: RHO1(NCHI,NCHI)
      COMPLEX(8)               :: MAT(NCHI,NCHI)
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: CSVAR
      REAL(8)                  :: SVAR
      INTEGER(4)               :: IKPT,ISPIN,IB,I,J
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_GRHO')
!       
!     ==========================================================================
!     ==  GRHO                                                                ==
!     ==========================================================================
      GRHO(:,:,:,:)=(0.D0,0.D0)
      GRHOLAUR(:,:,:,:)=(0.D0,0.D0)
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPTL
          DO IB=1,NB
            DO J=1,NCHI
              MAT(:,J)=PIPSI(:,IB,IKPT,ISPIN)*CONJG(PIPSI(J,IB,IKPT,ISPIN)) &
   &                                         *WKPTL(IKPT)      
            ENDDO
            DO NU=1,NOMEGA
              CSVAR=1.D0/(CI*OMEGA(NU)-ERHO(IB,IKPT,ISPIN))
              GRHO(:,:,NU,ISPIN)=GRHO(:,:,NU,ISPIN)+MAT*CSVAR
            ENDDO
            GRHOLAUR(:,:,1,ISPIN)=GRHOLAUR(:,:,1,ISPIN)+MAT
            SVAR=ERHO(IB,IKPT,ISPIN)
            GRHOLAUR(:,:,2,ISPIN)=GRHOLAUR(:,:,2,ISPIN)+MAT*SVAR
            SVAR=SVAR*ERHO(IB,IKPT,ISPIN)
            GRHOLAUR(:,:,3,ISPIN)=GRHOLAUR(:,:,3,ISPIN)+MAT*SVAR
          ENDDO
        ENDDO
      ENDDO
!       
!     ==========================================================================
!     ==  RHO                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,GRHO,GRHOLAUR,RHO)
        DO ISPIN=1,NSPIN
          WRITE(*,FMT='(82("="),T10," DENSITY MATRIX FROM GRHO ")')
          DO I=1,NCHI
            WRITE(*,FMT='("RHO",100("(",2F12.5,")"))')RHO(I,:,ISPIN)
          ENDDO
        ENDDO
      END IF
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_TEST()
!     **************************************************************************
!     ** COMPARES DIFFERENT WAYS TO CALCULATE THE GREENS FUNCTION AND         **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NSPIN,NDIM,NOMEGA,OMEGA,KBT,MU &
     &                       ,WKPTL,SMAT,SINV,HRHO,ERHO,NB,PIPSI
      IMPLICIT NONE
      COMPLEX(8)       :: SIGMA(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: SIGLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)       :: G(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: GLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)       :: GS(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: GSLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)       :: RHO(NCHI,NCHI,NSPIN)
      INTEGER(4)       :: IKPT,ISPIN,I
!
      COMPLEX(8)       :: CI=(0.D0,1.D0)
      COMPLEX(8)       :: CSVAR
      INTEGER(4)       :: NU,IB,J
      COMPLEX(8)       :: G1(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: G2(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: G3(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8)       :: G1LAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)       :: G2LAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)       :: G3LAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8)       :: RHOK0(NCHI,NCHI,NKPTL)
      COMPLEX(8)       :: RHOK1(NCHI,NCHI,NKPTL)
      COMPLEX(8)       :: RHOK2(NCHI,NCHI,NKPTL)
      COMPLEX(8)       :: RHOK3(NCHI,NCHI,NKPTL)
      COMPLEX(8)       :: RHO0(NCHI,NCHI,Nspin)
      COMPLEX(8)       :: RHO1(NCHI,NCHI,Nspin)
      COMPLEX(8)       :: RHO2(NCHI,NCHI,Nspin)
      COMPLEX(8)       :: RHO3(NCHI,NCHI,Nspin)
      COMPLEX(8)       :: MAT1(NCHI,NCHI)
      COMPLEX(8)       :: MAT2(NCHI,NCHI)
      COMPLEX(8)       :: BMAT1(NB,NB)
      COMPLEX(8)       :: BMAT2(NB,NB)
      COMPLEX(8)       :: BMAT3(NB,NB)
!     **************************************************************************
      SIGMA=(0.D0,0.D0)
      SIGLAUR=(0.D0,0.D0)
      
      GS=(0.D0,0.D0)
      GSLAUR=(0.D0,0.D0)
      DO IKPT=1,NKPTL
        CALL DMFT_GOFK(NCHI,NKPTL,NSPIN,NOMEGA,OMEGA,KBT,MU &
     &                    ,SMAT,HRHO,SIGMA,SIGLAUR,IKPT,G,GLAUR)
! SET TO ZERO BECAUSE THESE TERMS ARE NOT IMPLEMENTED IN THE COMPARISONS
GLAUR(:,:,2:,:)=(0.D0,0D0)
        CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,G,GLAUR,RHOK0(:,:,IKPT))
        GS=GS+WKPTL(IKPT)*G
        GSLAUR=GSLAUR+WKPTL(IKPT)*GLAUR
      ENDDO
      CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,GS,GSLAUR,RHO0)
!
DO ISPIN=1,NSPIN
PRINT*,'RHO FROM GLOC IN DMFT_TEST',ISPIN
  DO I=1,NCHI
    WRITE(*,FMT='("RHO",10("(",2F10.5,")"))')RHO0(I,:,ISPIN)
  ENDDO
ENDDO
!
!     ==========================================================================
!     == GREENS FUNCTION CALCULATED IN THE SPACE OF CORRELATED ORBITALS ========
!     ==========================================================================
      G1=(0.D0,0.D0)
      G1LAUR=(0.D0,0.D0)
      DO IKPT=1,NKPTL
        G=(0.D0,0.D0)
        GLAUR=(0.D0,0.D0)
        DO ISPIN=1,NSPIN
          DO NU=1,NOMEGA
            MAT1=CI*OMEGA(NU)*SMAT(:,:,IKPT,ISPIN)-HRHO(:,:,IKPT,ISPIN)
            CALL LIB$INVERTC8(NCHI,MAT1,MAT2)
            G(:,:,NU,ISPIN)=MAT2
          ENDDO
          GLAUR(:,:,1,ISPIN)=SINV(:,:,IKPT,ISPIN)
        ENDDO
        CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,G,GLAUR,RHOK1(:,:,IKPT))
!
        G1=G1+WKPTL(IKPT)*G
        G1LAUR=G1LAUR+WKPTL(IKPT)*GLAUR
      ENDDO
      CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,G1,G1LAUR,RHO1)
!
DO ISPIN=1,NSPIN
PRINT*,'RHO FROM GLOC IN DMFT_TEST (G1)',ISPIN
  DO I=1,NCHI
    WRITE(*,FMT='("RHO",10("(",2F10.5,")"))')RHO1(I,:,ISPIN)
  ENDDO
ENDDO
!
!     ==========================================================================
!     == GREENS FUNCTION IN NAT. ORBITALS FROM HRHO ============================
!     ==========================================================================
      G2=(0.D0,0.D0)
      G2LAUR=(0.D0,0.D0)
      DO IKPT=1,NKPTL
        G=(0.D0,0.D0)
        GLAUR=(0.D0,0.D0)
        DO ISPIN=1,NSPIN
          BMAT1(:,:)=MATMUL(CONJG(TRANSPOSE(PIPSI(:,:,IKPT,ISPIN))) &
    &                       ,MATMUL(HRHO(:,:,IKPT,ISPIN),PIPSI(:,:,IKPT,ISPIN)))
          DO NU=1,NOMEGA
            BMAT2(:,:)=-BMAT1(:,:)
            DO IB=1,NB
              BMAT2(IB,IB)=BMAT2(IB,IB)+CI*OMEGA(NU)
            ENDDO
            CALL LIB$INVERTC8(NB,BMAT2,BMAT3)
            MAT2=MATMUL(PIPSI(:,:,IKPT,ISPIN) &
    &                   ,MATMUL(BMAT3,TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN)))))
            G(:,:,NU,ISPIN)=MAT2
          ENDDO
          MAT2=MATMUL(PIPSI(:,:,IKPT,ISPIN) &
    &                                  ,TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN))))
          GLAUR(:,:,1,ISPIN)=MAT2
        ENDDO
        CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,G,GLAUR,RHOK2(:,:,IKPT))

        G2=G2+WKPTL(IKPT)*G
        G2LAUR=G2LAUR+WKPTL(IKPT)*GLAUR
      ENDDO
      CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,G2,G2LAUR,RHO2)
!
DO ISPIN=1,NSPIN
PRINT*,'RHO FROM GLOC IN DMFT_TEST (G2)',ISPIN
  DO I=1,NCHI
    WRITE(*,FMT='("RHO",10("(",2F10.5,")"))')RHO2(I,:,ISPIN)
  ENDDO
ENDDO
!
!     ==========================================================================
!     == GREENS FUNCTION IN NAT. ORBITALS FROM ERHO ============================
!     ==========================================================================
      G3=(0.D0,0.D0)
      G3LAUR=(0.D0,0.D0)
      DO IKPT=1,NKPTL
        G=(0.D0,0.D0)
        GLAUR=(0.D0,0.D0)
        DO ISPIN=1,NSPIN
          DO NU=1,NOMEGA
            MAT2(:,:)=(0.D0,0.D0)
            DO IB=1,NB
              CSVAR=1.D0/(CI*OMEGA(NU)-ERHO(IB,IKPT,ISPIN))
              DO J=1,NCHI
                MAT2(:,J)=MAT2(:,J)+PIPSI(:,IB,IKPT,ISPIN)*CSVAR &
      &                      *CONJG(PIPSI(J,IB,IKPT,ISPIN))
              ENDDO
            ENDDO
            G(:,:,NU,ISPIN)=MAT2
          ENDDO
          MAT2=MATMUL(PIPSI(:,:,IKPT,ISPIN) &
    &                                  ,TRANSPOSE(CONJG(PIPSI(:,:,IKPT,ISPIN))))
          GLAUR(:,:,1,ISPIN)=MAT2
        ENDDO
        CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,G,GLAUR,RHOK3(:,:,IKPT))
        G3=G3+WKPTL(IKPT)*G
        G3LAUR=G3LAUR+WKPTL(IKPT)*GLAUR
      ENDDO
      CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,G3,G3LAUR,RHO3)
!
DO ISPIN=1,NSPIN
PRINT*,'RHO FROM GLOC IN DMFT_TEST (G3)',ISPIN
  DO I=1,NCHI
    WRITE(*,FMT='("RHO",10("(",2F10.5,")"))')RHO3(I,:,ISPIN)
  ENDDO
ENDDO
print*,'================================================================='
print*,'================================================================='
!
     WRITE(*,FMT='("K-inDEPENDENT DENSITY MATRICES SHOULD BE IDENTICAL")') 
     WRITE(*,FMT='("DEV(0,1)= ",E10.5)')MAXVAL(ABS(RHO0-RHO1))
     WRITE(*,FMT='("DEV(1,2)= ",E10.5)')MAXVAL(ABS(RHO1-RHO2))
     WRITE(*,FMT='("DEV(1,3)= ",E10.5)')MAXVAL(ABS(RHO1-RHO3))
     WRITE(*,FMT='("DEV(2,3)= ",E10.5)')MAXVAL(ABS(RHO2-RHO3))
     WRITE(*,FMT='("K-DEPENDENT DENSITY MATRICES SHOULD BE IDENTICAL")') 
     WRITE(*,FMT='("DEV(0,1)= ",E10.5)')MAXVAL(ABS(RHOK0-RHOK1))
     WRITE(*,FMT='("DEV(1,2)= ",E10.5)')MAXVAL(ABS(RHOK1-RHOK2))
     WRITE(*,FMT='("DEV(1,3)= ",E10.5)')MAXVAL(ABS(RHOK1-RHOK3))
     WRITE(*,FMT='("DEV(2,3)= ",E10.5)')MAXVAL(ABS(RHOK2-RHOK3))
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_GOFK(NCHI,NKPTL,NSPIN,NOMEGA,OMEGA,KBT,MU &
     &                    ,SMAT,H0,SIGMA,SIGLAUR,IKPT,GLOC,GLAUR)
!     **************************************************************************
!     ** CALCULATES THE GREENS FUNCTION FOR A SPECIFIC K-POINT                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NCHI
      INTEGER(4),INTENT(IN)    :: NOMEGA
      INTEGER(4),INTENT(IN)    :: NSPIN
      INTEGER(4),INTENT(IN)    :: NKPTL
      REAL(8)   ,INTENT(IN)    :: OMEGA(NOMEGA)
      REAL(8)   ,INTENT(IN)    :: KBT
      REAL(8)   ,INTENT(IN)    :: MU
      INTEGER(4),INTENT(IN)    :: IKPT
      COMPLEX(8),INTENT(IN)    :: SMAT(NCHI,NCHI,NKPTL,NSPIN)
      COMPLEX(8),INTENT(IN)    :: H0(NCHI,NCHI,NKPTL,NSPIN)
      COMPLEX(8),INTENT(IN)    :: SIGMA(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(IN)    :: SIGLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8),INTENT(OUT)   :: GLOC(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(OUT)   :: GLAUR(NCHI,NCHI,3,NSPIN)
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER     :: TPRINT=.FALSE.
      LOGICAL(4),PARAMETER     :: TTEST=.FALSE.
      COMPLEX(8)               :: GINV(NCHI,NCHI)
      COMPLEX(8)               :: SINV(NCHI,NCHI)
      COMPLEX(8)               :: MAT(NCHI,NCHI)
      COMPLEX(8)               :: G(NCHI,NCHI)
      COMPLEX(8)               :: CSVAR
      COMPLEX(8)               :: RHO(NCHI,NCHI,NSPIN)
      INTEGER(4)               :: ISPIN,NU,I,J
      INTEGER(4)               :: NFIL
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_GOFK')
      DO ISPIN=1,NSPIN
!       ========================================================================
!       ========================================================================
        CALL LIB$INVERTC8(NCHI,SMAT(:,:,IKPT,ISPIN),SINV)

!       ========================================================================
!       == LAURENT EXPANSION OF THE LOCAL GREENS FUNCTION FROM 0 TO 2
!       == REMEMBER THAT MU=0 !!!
!       ========================================================================
        GLAUR(:,:,1,ISPIN)=(0.D0,0.D0)
        GLAUR(:,:,2,ISPIN)=(0.D0,0.D0)
        GLAUR(:,:,3,ISPIN)=(0.D0,0.D0)
        MAT=SINV
        GLAUR(:,:,1,ISPIN)=GLAUR(:,:,1,ISPIN)+MAT
        MAT=H0(:,:,IKPT,ISPIN)+SIGLAUR(:,:,1,ISPIN)-MU*SMAT(:,:,IKPT,ISPIN)
        MAT=MATMUL(SINV,MATMUL(MAT,SINV))
        GLAUR(:,:,2,ISPIN)=GLAUR(:,:,2,ISPIN)+MAT
        MAT=H0(:,:,IKPT,ISPIN)+SIGLAUR(:,:,1,ISPIN)-MU*SMAT(:,:,IKPT,ISPIN)
        MAT=MATMUL(MAT,MATMUL(SINV,MAT))
        MAT=MAT+SIGLAUR(:,:,2,ISPIN)
        MAT=MATMUL(SINV,MATMUL(MAT,SINV))
        GLAUR(:,:,3,ISPIN)=GLAUR(:,:,3,ISPIN)+MAT
!
!       ========================================================================
!       == LOCAL GREENS FUNCTION ON THE MATSUBARA GRID
!       ========================================================================
        GLOC(:,:,:,ISPIN)=(0.D0,0.D0)
        DO NU=1,NOMEGA
          CSVAR=CI*OMEGA(NU)+MU
          GINV(:,:)=SMAT(:,:,IKPT,ISPIN)*CSVAR &
     &             -H0(:,:,IKPT,ISPIN)-SIGMA(:,:,NU,ISPIN)
          CALL LIB$INVERTC8(NCHI,GINV,G)
          GLOC(:,:,NU,ISPIN)=GLOC(:,:,NU,ISPIN)+G
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == CALCULATE DENSITY MATRIX FROM THE GREENS FUNCTION FOR TESTING        ==
!     ==========================================================================
      IF(TTEST) THEN
        CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,GLOC,GLAUR,RHO)
        DO ISPIN=1,NSPIN
          WRITE(*,FMT='(82("="),T10," DENSITY MATRIX FROM GOFK ")')
          DO I=1,NCHI
            WRITE(*,FMT='("RHO",100("(",2F12.5,")"))')RHO(I,:,ISPIN)
          ENDDO
        ENDDO
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
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER     :: TPRINT=.FALSE.
      LOGICAL(4),PARAMETER     :: TTEST=.false.
      COMPLEX(8)               :: GINV(NCHI,NCHI)
      COMPLEX(8)               :: MAT(NCHI,NCHI)
      COMPLEX(8)               :: G(NCHI,NCHI)
      COMPLEX(8)               :: CSVAR
      COMPLEX(8)               :: RHO(NCHI,NCHI,NSPIN)
      INTEGER(4)               :: IKPT,ISPIN,NU,I,J
      INTEGER(4)               :: NFIL
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT_GLOC')
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
!
!     ==========================================================================
!     == CALCULATE DENSITY MATRIX FROM THE GREENS FUNCTION FOR TESTING        ==
!     ==========================================================================
      IF(TTEST) THEN
        CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,GLOC,GLAUR,RHO)
        DO ISPIN=1,NSPIN
          WRITE(*,FMT='(82("="),T10," DENSITY MATRIX FROM GLOC ")')
          DO I=1,NCHI
            WRITE(*,FMT='("RHO",100("(",2F12.5,")"))')RHO(I,:,ISPIN)
          ENDDO
        ENDDO
!STOP 'FORCED'
      END IF
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$ADDTOHPSI()
!     **************************************************************************
!     ** 
!     **************************************************************************
      USE DMFT_MODULE, ONLY  : TON,NKPTL,NSPIN,NB,h0,hrho,pipsi,nchi &
     &                        ,iproofchi
      USE WAVES_MODULE, ONLY : GSET,WAVES_SELECTWV,THIS,ndim,map
      IMPLICIT NONE
      REAL(8)   ,PARAMETER   :: MINOCC=1.D-1
      complex(8),parameter   :: ci=(0.d0,1.d0)
      COMPLEX(8),ALLOCATABLE :: dhpipsi(:,:)
      logical(4)             :: tchk,treset
      INTEGER(4)             :: NBH
      INTEGER(4)             :: npro
      INTEGER(4)             :: IKPT,ISPIN,IB,IDIM,ibh,ichi,ipro
!     **************************************************************************
      IF(.NOT.TON) RETURN
PRINT*,'ENTERING DMFT$ADDTOHPSI'
                              CALL TRACE$PUSH('DMFT$ADDTOHPSI')
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
      ALLOCATE(DHPIPSI(NCHI,NB))
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
            CALL ERROR$STOP('DMFT$ADDTOHPSI')
          END IF
          NBH=THIS%NBH
!
!         ======================================================================
!         == DF/DRHO=H0-HRHO
!         == DHDPIPSI= (H0-HRHO)<PI|PSI>; 
!         ======================================================================
          DHPIPSI(:,:)=MATMUL(H0(:,:,ikpt,ispin)-HRHO(:,:,ikpt,ispin) &
      &                     ,PIPSI(:,:,IKPT,ISPIN))
!
!         ======================================================================
!         ==  DETERMINE CONTRIBUTION TO PROJECTOR PART OF H|PSI>
!         ======================================================================
          IF(TRESET) THEN
            IF(.NOT.ASSOCIATED(THIS%HTBC))ALLOCATE(THIS%HTBC(NDIM,NBH,NPRO))
            THIS%HTBC=(0.D0,0.D0)
          END IF
!
          DO ICHI=1,NCHI
            IPRO=IPROOFCHI(ICHI)
            DO IBH=1,NBH
              IF(NBH.NE.NB) THEN
!               ==  PIPSI(ICHI,2*IBH-1,IKPT,ISPIN)= REAL(THIS%TBC(1,IBH,IPRO)) =
!               ==  PIPSI(ICHI,2*IBH  ,IKPT,ISPIN)=AIMAG(THIS%TBC(1,IBH,IPRO)) =
                THIS%HTBC(1,IBH,IPRO)=THIS%HTBC(1,IBH,IPRO) &
         &                           +   REAL(DHPIPSI(ICHI,2*IBH-1)) &
         &                           -CI*REAL(DHPIPSI(ICHI,2*IBH))
              ELSE
!               == PIPSI(ICHI,IBH,IKPT,ISPIN)  =THIS%TBC(1,IBH,IPRO) ===========
                THIS%HTBC(1,IBH,IPRO)=THIS%HTBC(1,IBH,IPRO)+DHPIPSI(ICHI,IBH)
              END IF
            ENDDO
          ENDDO
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
      SUBROUTINE testbessel()
      integer(4)           :: l=2
      integer(4),parameter :: nkappa=5
      integer(4)           :: i
      real(8)              :: rad
      real(8)              :: k2(nkappa)=(/-0.2d0,-0.1d0,0.d0,1.d0,2.d0/)
      real(8)              :: jval(nkappa),jder(nkappa)
      real(8)              :: kval(nkappa),kder(nkappa)
      integer(4)           :: nfil1=12,nfil2=14
      open(nfil1,file='jval.dat')
      open(nfil2,file='kval.dat')
      do i=1,1000
        rad=1.d-2*real(i)
        do ikappa=1,nkappa
          CALL LMTO$SOLIDBESSELRAD(L,RAD,K2(ikappa),JVAL(ikappa),JDER(ikappa))
          CALL LMTO$SOLIDHANKELRAD(L,RAD,K2(ikappa),KVAL(ikappa),KDER(ikappa))
        enddo
        write(nfil1,*)rad,jval,jder
        write(nfil2,*)rad,kval,kder
      enddo
      close(nfil1)
      close(nfil2)
print*,'.... besseltest done'
stop 'forced after besseltest'
      return
      end
