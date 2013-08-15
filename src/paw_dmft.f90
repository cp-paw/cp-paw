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
TYPE denmat_TYPE
INTEGER(4)          :: LNx
INTEGER(4),POINTER  :: LOX(:) !(lnx)
INTEGER(4),POINTER  :: lmn(:) !(nloc)
integer(4)          :: lmnx
REAL(8)   ,POINTER  :: MAT(:,:,:)  ! density matrix
REAL(8)   ,POINTER  :: H(:,:,:)    ! Hamiltonian from double counting
END TYPE DENMAT_TYPE
TYPE ATOMSET_TYPE
  INTEGER(4)           :: NLOC
  REAL(8)   ,POINTER   :: U(:,:,:,:) => NULL() !(NLOC,NLOC,NLOC,NLOC) 
  COMPLEX(8),POINTER   :: GLOC(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: GLOCLAUR(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,3)
  COMPLEX(8),POINTER   :: SLOC(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,NOMEGA)
  COMPLEX(8),POINTER   :: SLOCLAUR(:,:,:,:) => NULL() !(NLOC,NLOC,NDIMD,3)
  type(denmat_type)    :: denmat
END TYPE ATOMSET_TYPE
TYPE KSET_TYPE
  REAL(8)            :: WKPTL
  logical(4)         :: taddminusk    !add  the term for -k
  REAL(8)   ,POINTER :: ERHO(:,:)     !(NB,NSPIN)
  COMPLEX(8),POINTER :: PIPSI(:,:,:,:)!(NDIM,NCHI,NB,NSPIN)
  COMPLEX(8),POINTER :: RHO(:,:,:)    !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: H0(:,:,:)     !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: HRHO(:,:,:)   !(NCHI,NCHI,NDIMD)
  COMPLEX(8),POINTER :: SINV(:,:,:)   !(NCHI,NCHI,ndimd)
  COMPLEX(8),POINTER :: SMAT(:,:,:)   !(NCHI,NCHI,ndimd)
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
INTEGER(4)             :: NDIMd         ! can be 1,2,4
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
!== U-TENSOR ===================================================================
real(8)   ,allocatable :: uchi(:,:,:,:)
!== kset =======================================================================
type(kset_type)   ,allocatable :: kset(:)  !(nkptl)
type(atomset_type),allocatable :: atomset(:)  !(nat)
END MODULE DMFT_MODULE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_INI()
!     **************************************************************************
      USE DMFT_MODULE, only: tini,ndim,ndimd,nspin,nkptl,nb,nat,nchi &
     &                       ,nomega,kbt,deltat,mu &
     &                       ,wkptl,omega,ichibnd,iproofchi &
     &                       ,erho,pipsi,hrho,rhoofk,sigmadc &
     &                       ,smat,sinv,diagsloc &
     &                       ,gloc,gloclaur,sigma,siglaur &
     &                       ,kset,atomset
      USE WAVES_MODULE, ONLY : KMAP,NDIM_W=>NDIM,NKPTL_W=>NKPTL,NSPIN_W=>NSPIN
      IMPLICIT NONE
      REAL(8)                :: PI
      INTEGER(4)             :: ISP0
      INTEGER(4)             :: L0
      REAL(8)   ,ALLOCATABLE :: WKPT(:)    ! K-POINT WEIGHTS (GLOBAL)
      INTEGER(4)             :: NTASKS_K,THISTASK_K
      INTEGER(4)             :: NTASKS_M,THISTASK_M
      INTEGER(4)             :: NKPT
      INTEGER(4)             :: lnx
      INTEGER(4)             :: lmnx
      INTEGER(4)             :: nloc
      INTEGER(4)             :: l
      INTEGER(4)             :: NU,ISP,iat,LN,im,IKPTL,IKPT,ichi,ipro,ispin,i
      INTEGER(4)             :: lmn
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
      NOMEGA=40
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
!
!     == kset ==================================================================
      allocate(kset(nkptl))
      do ikptl=1,nkptl
        kset(ikptl)%wkptl=wkptl(ikptl)
        allocate(kset(ikptl)%erho(nb,nspin))
        allocate(kset(ikptl)%pipsi(ndim,nchi,nb,nspin))
        allocate(kset(ikptl)%rho(nchi,nchi,ndimd))
        allocate(kset(ikptl)%h0(nchi,nchi,ndimd))
        allocate(kset(ikptl)%hrho(nchi,nchi,ndimd))
        allocate(kset(ikptl)%sinv(nchi,nchi,ndimd))
        allocate(kset(ikptl)%smat(nchi,nchi,ndimd))
        kset(ikptl)%erho=0.d0
        kset(ikptl)%pipsi=(0.d0,0.d0)
        kset(ikptl)%rho=(0.d0,0.d0)
        kset(ikptl)%hrho=(0.d0,0.d0)
        kset(ikptl)%h0=(0.d0,0.d0)
        kset(ikptl)%sinv=(0.d0,0.d0)
        kset(ikptl)%smat=(0.d0,0.d0)
      enddo
!
!     == atomset ===============================================================
      allocate(atomset(nat))
      do iat=1,nat
        nloc=ichibnd(2,iat)-ichibnd(1,iat)+1
        atomset(iat)%nloc=nloc
        allocate(atomset(iat)%u(nloc,nloc,nloc,nloc))
        allocate(atomset(iat)%gloc(nloc,nloc,ndimd,nomega))
        allocate(atomset(iat)%gloclaur(nloc,nloc,ndimd,3))
        allocate(atomset(iat)%sloc(nloc,nloc,ndimd,nomega))
        allocate(atomset(iat)%sloclaur(nloc,nloc,ndimd,3))
        atomset(iat)%u=0.d0
        atomset(iat)%gloc=(0.d0,0.d0)
        atomset(iat)%gloclaur=(0.d0,0.d0)
        atomset(iat)%sloc=(0.d0,0.d0)
        atomset(iat)%sloclaur=(0.d0,0.d0)
      enddo
!
!     == atomset%denmat%... for double counting correction =====================
      DO IAT=1,NAT
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        ALLOCATE(torb(LNX))
        CALL SETUP$GETL4A('TORB',LNX,TORB)
        CALL SETUP$unSELECT()
!
        ATOMSET(IAT)%DENMAT%LNX=LNX
        ALLOCATE(ATOMSET(IAT)%DENMAT%LOX(LNX))
        ATOMSET(IAT)%DENMAT%LOX=LOX
        LMNX=SUM(2*LOX+1)
        ATOMSET(IAT)%DENMAT%LMNX=LMNX
        nloc=atomset(iat)%nloc
        ALLOCATE(ATOMSET(IAT)%DENMAT%LMN(NLOC))
        LMN=0
        i=0
        do ln=1,lnx
          l=lox(ln)
          if(torb(ln)) then
            do im=1,2*l+1
              i=i+1
              lmn=lmn+1
              ATOMSET(IAT)%DENMAT%lmn(i)=lmn
            enddo
          else
            lmn=lmn+2*l+1
          end if
        enddo
        ALLOCATE(ATOMSET(IAT)%DENMAT%MAT(LMNX,LMNX,NDIMD))
        ALLOCATE(ATOMSET(IAT)%DENMAT%h(LMNX,LMNX,NDIMD))
        ATOMSET(IAT)%DENMAT%MAT=0.d0
        ATOMSET(IAT)%DENMAT%h=0.d0
        deallocate(lox)
        deallocate(torb)
      ENDDO        
!
!     == now the other allocations =============================================
      ALLOCATE(ERHO(NB,NKPTL,NSPIN))
      ALLOCATE(PIPSI(NCHI,NB,NKPTL,NSPIN))
      ALLOCATE(HRHO(NCHI,NCHI,NKPTL,Ndimd))
      ALLOCATE(RHOOFK(NCHI,NCHI,NKPTL,Ndimd))
!     == H0 IS PURPOSELY NOT ALLOCATED
      ALLOCATE(SIGMADC(NCHI,NCHI,Ndimd))
      ALLOCATE(SMAT(NCHI,NCHI,NKPTL,NSPIN))
      ALLOCATE(SINV(NCHI,NCHI,NKPTL,NSPIN))
      ALLOCATE(GLOC(NCHI,NCHI,NOMEGA,Ndimd))
      ALLOCATE(GLOCLAUR(NCHI,NCHI,3,Ndimd))
      ALLOCATE(DIAGSLOC(NCHI,NCHI,NSPIN))
      ALLOCATE(SIGMA(NCHI,NCHI,NOMEGA,Ndimd))
      ALLOCATE(SIGLAUR(NCHI,NCHI,3,Ndimd))
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
      USE DMFT_MODULE, ONLY: TON,NB,NCHI,NKPTL,NSPIN,NDIM,ndimd,nat,NOMEGA &
     &                     ,OMEGA,KBT,MU &
     &            ,WKPTL,H0,HRHO,SIGMA,SIGLAUR,GLOC,GLOCLAUR,DIAGSLOC &
     &            ,uchi,kset,ichibnd,smat
      USE MPE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: I,ITER,IKPT,ISPIN
      LOGICAL(4)             :: TPRINT=.FALSE.
!     **************************************************************************
      IF(.NOT.TON) RETURN
                              CALL TRACE$PUSH('DMFT$GREEN')
print*,'before dmft_ini'
      CALL DMFT_INI()
!
!     ==========================================================================
!     ==  COLLECT DFT HAMILTONIAN  (HAMILTONIAN INITIALIZED TO ERHO)          ==
!     ==========================================================================
print*,'old version'
      CALL DMFT$COLLECTHAMILTONIAN()   !old
print*,'new version'
      CALL DMFT$COLLECTHAMILTONIAN_WITHKSET()  !new
      call DMFT$COLLECTfulldenmat()  !new
!
!     ==========================================================================
!     ==  OBTAIN BARE U-TENSOR FROM LMTO OBJECT.                              ==
!     ==     SHAPE OF ORBITALS ARE DEFINED BY NPRO AND ATOMIC STRUCTRE        ==
!     ==     ORBITALS ARE SELECTED BY TORB.                                   ==
!     ==========================================================================
print*,'setting up u-tensor...'
      allocate(uchi(nchi,nchi,nchi,nchi))
      CALL DMFT_UCHI() !old
      call DMFT_Uchi_withatomset() !new
print*,'... u-tensor done'
!
!     ==========================================================================
!     ==  TRANSFORMATION ONTO A ORTHONORMAL CORRELATED BASIS SET              ==
!     ==    |CHIORTHO>   =|CHINONORTHO>*DIAGSLOC                              ==
!     ==========================================================================
      CALL DMFT_DIAGSLOC() !old
      CALL DMFT_PRINTspinorMATRIX(6,'DIAGSLOC(1)','DIRECT',nat,ichibnd &
     &                           ,NDIMD,NCHI,DIAGSLOC)
!
      CALL DMFT_DIAGSLOC_WITHKSET() !new
      CALL DMFT_PRINTSPINORMATRIX(6,'DIAGSLOC(2)','UPDOWN',NAT,ICHIBND &
     &                           ,NDIMD,NCHI,DIAGSLOC)
!
!      CALL DMFT_GRHO() ! TEST ONLY: DENSITY FROM REFERENCE GREENS FUNCTION
!
!     ==========================================================================
!     ==  CONSTRUCT NON-INTERACTING HAMILTONIAN THAT PRODUCES THE CORRECT     ==
!     ==  ONE-PARTICLE DENSITY MATRIX                                         ==
!     ==========================================================================
print*,'setting up non-interacting hamiltonian hrho...'
      HRHO=(0.D0,0.D0)
      SIGMA=(0.D0,0.D0)
      SIGLAUR=(0.D0,0.D0)
      CALL DMFT_CONSTRAINTS(HRHO,SIGMA,SIGLAUR) !old
print*,'now with kset'
      CALL DMFT_CONSTRAINTS_WITHKSET('HRHO') !new
PRINT*,'... HRHO DONE'
!
      IF(.NOT.ALLOCATED(H0)) THEN
        ALLOCATE(H0(NCHI,NCHI,NKPTL,NSPIN))
        H0=HRHO
      END IF
      CALL DMFT_PRINTSPINORMATRIX(6,'h0','UPDOWN',NAT,ICHIBND &
     &                           ,NDIMD,NCHI,kset(1)%h0)
!
!      CALL DMFT_TEST()  !TEST ONLY
!
!     ==========================================================================
!     == ITERATION TO ENFORCE CONSTRAINTS                                     ==
!     ==========================================================================
      MU=0.D0
      DO ITER=1,10
WRITE(*,FMT='(82("="),T20," ITERATION ",I5)')ITER
        CALL DMFT_GLOC(H0,SIGMA,SIGLAUR,GLOC,GLOCLAUR) !old
        call DMFT_GLOC_withatomset() !new

!       CALL DMFT$PLOTGLOC(-'GLOC.DAT')
!
!       ========================================================================
!       ==  WRITE LOCAL FILE FOR SOLVER                                       ==
!       ========================================================================
        CALL DMFT_WRITEGLOC(NCHI,NOMEGA,NSPIN,KBT,MU,DIAGSLOC,OMEGA &
     &                   ,GLOC,GLOCLAUR)
!
!       ========================================================================
!       ==  CALL THE SOLVER                                                   ==
!       ========================================================================
        call DMFT$SOLVER_withkset() !new
!!$        CALL DMFT$SOLVER() !old
!!$print*,'stopping'
!!$stop 'forced'
!
!       ========================================================================
!       ==  READ SELF ENERGY FROM FILE                                        ==
!       ========================================================================
!!$        CALL DMFT_GETSIGMA()
!!$!      CALL DMFT$PLOTSIGMA(-'SIG.DAT')
!!$!
!!$!       ========================================================================
!!$!       ==  CONSTRAINTS                                                       ==
!!$!       ========================================================================
!!$        CALL DMFT_CONSTRAINTS(H0,SIGMA,SIGLAUR)
        CALL DMFT_CONSTRAINTS_WITHKSET('H0')
print*,'forced stop after constraints_withkset("h0")'
stop
!
!       ========================================================================
!       == MAP ONTO HAMILTONIAN                                               ==
!       ========================================================================
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
      ENDDO ! end of loop over iterations to enforce constraint
      deallocate(uchi)
!
      call DMFT$ADDTOHPSI_withkset()

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
!     **                                                                      **
!     **   sinv(k) =sum_n <pi|psi(n,k)><psi(n,k)|pi>                          **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NSPIN,ndim &
     &                      ,WKPTL,PIPSI,DIAGSLOC,SMAT,SINV
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
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
      SUBROUTINE DMFT_DIAGSLOC_withkset()
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
!     **   sinv(k) =sum_n <pi|psi(n,k)><psi(n,k)|pi>                          **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NSPIN,ndim,ndimd,kset,diagsloc &
& ,nat,ichibnd
      IMPLICIT NONE
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)  ! SQRT(-1)
      INTEGER(4)            :: IKPT,ISPIN,idim1,idim2,idimd,I
      COMPLEX(8)            :: MAT(NCHI,NCHI)
      COMPLEX(8),allocatable:: bmat1(:,:) !(2*NCHI,2*NCHI)
      COMPLEX(8),allocatable:: bmat2(:,:) !(2*NCHI,2*NCHI)
      COMPLEX(8),allocatable:: floc2(:)
      REAL(8)               :: FLOC(NCHI)
      REAL(8)               :: svar
      LOGICAL(4)            :: TTEST=.FALSE.
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
        kset(ikpt)%sinv(:,:,:)=(0.d0,0.d0)
        DO Ispin=1,nspin
          do idim2=1,ndim
            do idim1=1,ndim
!             == non-collinear               idimd=(1)
!             == collinear spin-polaried     idimd=(1,2)
!             == non-collinear spin-polaried idimd=(1,2,3,4)
              idimd=idim1+ndim*(idim2-1)+ispin-1
              mat=matmul(kset(ikpt)%pipsi(idim1,:,:,ispin) &
     &                  ,conjg(transpose(kset(ikpt)%pipsi(idim2,:,:,ispin))))
              kset(ikpt)%sinv(:,:,idimd)=mat
            enddo
          enddo
        ENDDO
        CALL SPINOR$CONVERT('FWRD',nchi,NDIMD,KSET(IKPT)%SINV)
      enddo
!
!     ==========================================================================
!     == OVERLAP MATRIX                                                       ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        call spinor$INVERT(NDIMD,Nchi,kset(ikpt)%sinv,kset(ikpt)%smat)
      ENDDO
!
!     ==========================================================================
!     == OBTAIN MATRIX THAT MAKES THE OVERLAP EQUAL TO THE UNIT MATRIX        ==
!     ==========================================================================
      DIAGSLOC=(0.D0,0.D0)
      DO IKPT=1,NKPTL
        DIAGSLOC=DIAGSLOC+KSET(IKPT)%WKPTL*KSET(IKPT)%SINV
      ENDDO
!
!     == make hermitean (broken, when only part of the k-points are included) ==
      do idimd=1,ndimd
        diagsloc(:,:,idimd)=0.5d0*(diagsloc(:,:,idimd) &
     &                            +conjg(transpose(diagsloc(:,:,idimd))))
      enddo
!
!     == convert from (txyz) to (up-down representation) =======================
      CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,DIAGSLOC)

      IF(NDIMD.EQ.1.OR.NDIMD.EQ.2) THEN
        DO IDIMD=1,NDIMD
          mat=diagsloc(:,:,idimd)
          CALL LIB$DIAGC8(NCHI,MAT,FLOC,DIAGSLOC(:,:,IDIMD))
          DO I=1,NCHI
            DIAGSLOC(:,I,IDIMD)=DIAGSLOC(:,I,IDIMD)/SQRT(FLOC(I))
          ENDDO
        ENDDO
      ELSE if(ndimd.eq.4) then  ! NON-COLLINEAR CASE
        ALLOCATE(BMAT1(2*NCHI,2*NCHI))
        ALLOCATE(BMAT2(2*NCHI,2*NCHI))
        ALLOCATE(FLOC2(2*NCHI))
        BMAT1(1:NCHI ,1:NCHI )=DIAGSLOC(:,:,1)      
        BMAT1(1:NCHI ,NCHI+1:)=DIAGSLOC(:,:,2)      
        BMAT1(NCHI+1:,1:NCHI )=DIAGSLOC(:,:,3)      
        BMAT1(NCHI+1:,NCHI+1:)=DIAGSLOC(:,:,4)      
        CALL LIB$DIAGC8(2*NCHI,BMAT1,FLOC2,BMAT2)
        DO I=1,2*NCHI
          BMAT2(:,I)=BMAT2(:,I)/SQRT(FLOC2(I))
        ENDDO
        DIAGSLOC(:,:,1)=BMAT2(1:NCHI ,1:NCHI )
        DIAGSLOC(:,:,2)=BMAT2(1:NCHI ,NCHI+1:)      
        DIAGSLOC(:,:,3)=BMAT2(NCHI+1:,1:NCHI )      
        DIAGSLOC(:,:,4)=BMAT2(NCHI+1:,NCHI+1:)
        deallocate(bmat1)
        deallocate(bmat2)
        deallocate(floc2)
      else
        CALL ERROR$MSG('illegal value of ndimd')
        CALL ERROR$STOP('DMFT_DIAGSLOC')
      END IF
!
!     == convert from (up-down representation to (txyz) ========================
      CALL spinor$convert('FWRD',NCHI,NDIMD,DIAGSLOC)
!
!     ==========================================================================
!     == OPTIONAL TESTS
!     ==========================================================================
      IF(TTEST) THEN
        DO IDIMD=1,NDIMD
          MAT=(0.D0,0.D0)
          DO IKPT=1,NKPTL
            MAT=MAT+KSET(IKPT)%WKPTL*KSET(IKPT)%SINV(:,:,idimd)
          END do
          mat=0.5d0*(mat+conjg(transpose(mat)))
!
          MAT=MATMUL(CONJG(TRANSPOSE(DIAGSLOC(:,:,IDIMD))) &
                                    ,MATMUL(MAT,DIAGSLOC(:,:,IDIMD)))
          WRITE(*,FMT='(82("="),T10,"  TEST DIAGSLOC FOR IDIMD=",I2,"  ")')IDIMD
          DO I=1,NCHI
            WRITE(*,FMT='("TEST",10("(",2F10.5,")"))')MAT(I,:)
          ENDDO
        ENDDO
        CALL ERROR$MSG('SCHEDULED STOP AFTER TEST')
        CALL ERROR$STOP('DMFT_DIAGSLOC')
      END IF
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
!     ==========================================================================
!     == CALCULATE INVERSE OVERLAP MATRIX                                     ==
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
          WRITE(*,FMT='(82("="),T10," DENSITY MATRIX FOR ATOM ",I5," FROM BAND OCCUPATIONS ")')IAT
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
      SUBROUTINE DMFT$COLLECTHAMILTONIAN_withkset()
!     **************************************************************************
!     ** COLLECTS THE HAMILTONIAN AND STORES IT ON THE MODULE                 **
!     **                                                                      **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NDIM,ndimd,IPROOFCHI,KBT &
     &                      ,ichibnd,nat,kset
      USE MPE_MODULE
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      COMPLEX(8),allocatable :: RHO(:,:,:) !(NCHI,NCHI,Ndimd)
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: NGL
      INTEGER(4)             :: IKPT,ISPIN,IBH,ICHI,IPRO,IB,I,J,iat,idim1,idim2
      INTEGER(4)             :: Idimd
      REAL(8)                :: F(NB,NKPTL,NSPIN)
      REAL(8)                :: SVAR
      REAL(8)                :: emin,emax
      COMPLEX(8)             :: mat(NCHI,NCHI)
      COMPLEX(8)             :: CSVAR
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                      CALL TRACE$PUSH('DMFT$COLLECTHAMILTONIAN')
!
!     ==========================================================================
!     ==  DETERMINE ERHO 
!     ==========================================================================
      CALL DMFT$COLLECTOCCUPATIONS(NKPTL,NSPIN,NB,F)
      emin=+huge(emin)
      emax=-huge(emax)
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          DO IB=1,NB
            CALL DMFT_EOFF(KBT,F(IB,IKPT,ISPIN),kset(ikpt)%erho(ib,ispin))
            emin=min(emin,kset(ikpt)%erho(ib,ispin))
            emax=max(emax,kset(ikpt)%erho(ib,ispin))
          ENDDO
        ENDDO
      ENDDO 
PRINT*,'BOUNDS OF SPECTRUM [EV] ',emin*27.211D0,eMAX*27.211D0
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
!         == specify if second k-point is obtained from time-inversion.
!         == non-collinear: no, because there is no time inversion symmetry
!         == tinv: no, k-point falls onto itself under time-inversion symmetry.
          kset(ikpt)%taddminusk=(ndimd.ne.4).and.(.not.gset%tinv)
!
          DO ICHI=1,NCHI
            IPRO=IPROOFCHI(ICHI)
            DO IBH=1,NBH
              if(nbh.ne.nb) then
                kset(ikpt)%pipsi(:,ichi,2*ibh-1,ispin) &
     &                                              = REAL(THIS%TBC(:,IBH,IPRO))
                kset(ikpt)%pipsi(:,ichi,2*ibh  ,ispin)& 
     &                                              =aimag(THIS%TBC(:,IBH,IPRO))
              else
                kset(ikpt)%pipsi(:,ichi,ibh,ispin)=THIS%TBC(:,IBH,IPRO)
              end if
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==   rho(A,S1,B,S2,K)=SUM_N <PI(A,S1)|PSI(N,K)>f_n <PSI(N,K)|PI(B,S2)>  ==
!     ==                                                                      ==
!     ==   NON-SPIN POLARIZED:                                                ==
!     ==     rho(S1,S2)=1/2*[ rho(1)*(1,0/0,1) ]                              ==
!     ==   COLLINEAR SPIN-POLARIZED                                           ==
!     ==     RHO(S1,S2)=1/2*[ RHO(1)*(1,0/0,1)+RHO(2)*(1,0/0,-1) ]            ==
!     ==   NON-COLLINEAR SPIN-POLARIZED                                       ==
!     ==     RHO(S1,S2)=1/2*[ RHO(1)*(1,0/0,1)+RHO(2)*(0,1/1,0)               ==
!     ==                      +RHO(3)*(0,-I/I,0)+RHO(2)*(1,0/0,-1) ]          ==
!     ==   WHERE                                                              ==
!     ==     RHO(1)=   RHO(1,1)+RHO(2,2)                                      ==
!     ==     RHO(2)=   RHO(2,1)+RHO(1,2)                                      ==
!     ==     RHO(3)=-i(RHO(2,1)-RHO(1,2))                                     ==
!     ==     RHO(4)=   RHO(1,1)-RHO(2,2)  !FOR NSPIN=2 THIS IS SIMV(2)        ==
!     ==                                                                      ==
!     ==  Caution! matrices are arranged in a fortran-style notation with the ==
!     ==     first index running first. The notation (a,b/c,d) however is a   ==
!     ==     latex style notation with the second index running first!        ==
!     ==                                                                      ==
!     ==========================================================================
!     ==  ndim=1,nspin=1  ndimd=1: (t)                                        ==
!     ==  ndim=1,nspin=2  ndimd=2: (t,z)                                      ==
!     ==  ndim=2,nspin=1  ndimd=4: (t,x,y,z)                                  ==
!     ==========================================================================
      DO IKPT=1,NKPTL
        kset(ikpt)%rho=(0.d0,0.d0)
        DO ISPIN=1,NSPIN
          do idim1=1,ndim
            do idim2=1,ndim
              mat(:,:)=(0.d0,0.d0)
              DO IB=1,NB
                DO J=1,NCHI
                  CSVAR=F(IB,IKPT,ISPIN) &
       &               *CONJG(kset(ikpt)%PIPSI(idim2,j,IB,ispin))
                  mat(:,j)=mat(:,j) &
       &                  +kset(ikpt)%pipsi(idim1,:,ib,ispin)*csvar
                enddo
              ENDDO
!
!             ==================================================================
!             == DISTRIBUTE DENSITY MATRIX ONTO VARIOUS SPIN COMPONENTS       ==
!             == using the up-down notation (1=uu,2=du,3=ud,4=dd)             ==
!             == nspin=1,2: idim1=idim2=ndim=1 =>idimd=ispin                  ==
!             == nspin=1,ndim=2,ispin=1: (1,2,3,4)                            ==
!             ==================================================================
              idimd=idim1+ndim*(idim2-1)+ispin-1
              KSET(IKPT)%RHO(:,:,idimd)=MAT
            ENDDO
          ENDDO
        ENDDO
!
!       ========================================================================
!       == convert from up-down representation into (txyz)                    ==
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
!       == make hermitean ======================================================
        do idimd=1,ndimd
          rho(:,:,idimd)=0.5d0*(rho(:,:,idimd)+conjg(transpose(rho(:,:,idimd))))
        enddo
!       == print ===============================================================
        print*,' in DMFT$COLLECTHAMILTONIAN_withkset'
        CALL DMFT_PRINTSPINORMATRIX(6,'DENSITY MATRIX(UPDOWN)','UPDOWN' &
     &                             ,NAT,ICHIBND,NDIMD,NCHI,RHO)
        DEALLOCATE(RHO)
      END IF
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$COLLECTfulldenmat()
!     **************************************************************************
!     ==  EXTRACT onsite density matrices with all projector functions        ==
!     ==  for double counting                                                 ==
!     **************************************************************************
      USE DMFT_MODULE, ONLY: nkptl,nspin,nat,ndim,ndimd,nb,atomset
      USE WAVES_MODULE, ONLY : GSET,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      INTEGER(4)             :: NBH     !#(SUPER STATES)
      INTEGER(4)             :: lnx
      INTEGER(4)             :: lmnx
      INTEGER(4)             :: i1,i2
      COMPLEX(8)             :: CSVAR
      COMPLEX(8),allocatable :: mat(:,:,:)
      INTEGER(4)             :: nb_,nspin_
      INTEGER(4)             :: nloc
      INTEGER(4)             :: iat,isp,ikpt,ispin,ipro,ibh,idim1,idim2,idimd,i
      INTEGER(4)             :: lmn,im,l,ln
      real(8)                :: occ(nb,nkptl,nspin)
!     **************************************************************************
!
!     ==========================================================================
!     == collect occupations                                                  ==
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
      CALL WAVES_DYNOCCGETR8A('OCC',NB*NKPTL*NSPIN,occ)
!
!     ==========================================================================
!     == set density matrix to zero
!     ==========================================================================
      DO IAT=1,NAT
        atomset(iat)%denmat%mat=(0.d0,0.d0)
      ENDDO        
!
!     ==========================================================================
!     == add up density matrix                                                ==
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
            lmnx=ATOMset(IAT)%denmat%LMNX
            IF(ATOMset(IAT)%NLOC.LE.0) THEN
              IPRO=IPRO+lmnx
              CYCLE
            END IF   
            I1=IPRO+1         
            I2=IPRO+lmnx
            DO IBH=1,NBH
              IF(NBH.NE.NB) THEN
!               == this is for general special k-points ========================
                DO IDIM2=1,NDIM
                  DO IDIM1=1,NDIM
                    IDIMD=IDIM1+NDIM*(IDIM2-1)+(ISPIN-1)
                    DO lmn=1,lmnx
                      i=ipro+lmn
                      CSVAR=REAL(THIS%TBC(IDIM1,IBH,I))*OCC(2*IBH-1,IKPT,ISPIN)
                      ATOMSET(IAT)%DENMAT%MAT(:,lmn,IDIMD) &
     &                               =ATOMSET(IAT)%DENMAT%MAT(:,lmn,IDIMD) &
     &                               +REAL(THIS%TBC(IDIM2,IBH,I1:I2)*CSVAR)
!
                      CSVAR=AIMAG(THIS%TBC(IDIM1,IBH,I))*OCC(2*IBH,IKPT,ISPIN)
                      ATOMSET(IAT)%DENMAT%MAT(:,lmn,IDIMD) &
     &                             =ATOMSET(IAT)%DENMAT%MAT(:,lmn,IDIMD) &
     &                             +AIMAG(THIS%TBC(IDIM2,IBH,I1:I2)*CSVAR)
                    ENDDO
                  ENDDO
                ENDDO
              ELSE
!               == this is for general k-points ================================
                DO IDIM2=1,NDIM
                  DO IDIM1=1,NDIM
                    IDIMD=IDIM1+NDIM*(IDIM2-1)+(ISPIN-1)
                    DO lmn=1,lmnx
                      i=ipro+lmn
                      CSVAR=CONJG(THIS%TBC(IDIM1,IBH,I))*OCC(IBH,IKPT,ISPIN)
                      ATOMSET(IAT)%DENMAT%MAT(:,lmn,IDIMD) &
     &                             =ATOMSET(IAT)%DENMAT%MAT(:,lmn,IDIMD) &
     &                             +real(THIS%TBC(IDIM2,IBH,I1:I2)*CSVAR)
                    ENDDO
                  ENDDO
                ENDDO
              END IF
            ENDDO  ! end of loop over bands
            IPRO=IPRO+lmnx
          ENDDO ! end of loop over atoms
        ENDDO ! end of loop over spin
      ENDDO ! end of loop over k-points
!
!     ==========================================================================
!     ==  transform to (t,x,y,z) representation and make hermitean            ==
!     ==========================================================================
      DO IAT=1,NAT
        IF(ATOMSET(IAT)%NLOC.LE.0) CYCLE
        lmnx=ATOMSET(IAT)%DENMAT%LMNX
        ALLOCATE(MAT(LMNX,LMNX,NDIMD))
        MAT=CMPLX(ATOMSET(IAT)%DENMAT%MAT)
        CALL SPINOR$CONVERT('FWRD',LMNX,NDIMD,MAT)
        DO IDIMD=1,NDIMD
          MAT(:,:,IDIMD)=0.5D0*(MAT(:,:,IDIMD)+CONJG(TRANSPOSE(MAT(:,:,IDIMD))))
        ENDDO
        ATOMSET(IAT)%DENMAT%MAT=REAL(MAT)
        DEALLOCATE(MAT)
      ENDDO
!
!     ==========================================================================
!     ==  print
!     ==========================================================================
!!$      DO IAT=1,NAT
!!$        IF(ATOMSET(IAT)%NLOC.LE.0) CYCLE
!!$        WRITE(*,FMT='(82("="),T10," DENSITY MATRIX FOR ATOM ",I3," ")')IAT
!!$        DO IDIMD=1,NDIMD
!!$          DO LMN=1,LMNX
!!$            WRITE(*,FMT='("IDIMD=",I1,":",100F10.5)')IDIMD &
!!$     &              ,ATOMSET(IAT)%DENMAT%MAT(LMN,:,IDIMD)
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
      return
      end
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
      USE DMFT_MODULE, ONLY: nat,NB,NCHI,NKPTL,NSPIN,NDIM,NOMEGA &
     &                      ,GLOC,GLOCLAUR,DIAGSLOC,uchi,ichibnd
      IMPLICIT NONE
      REAL(8)                :: KBT
      REAL(8)                :: MU
      REAL(8)                :: OMEGA(NOMEGA)
      real(8)   ,allocatable :: uloc1(:,:,:,:)     !(NCHI1,NCHI1,nchi1,nchi1)
      COMPLEX(8),allocatable :: GLOC1(:,:,:,:)     !(NCHI1,NCHI1,NOMEGA,NSPIN)
      COMPLEX(8),allocatable :: GLOCLAUR1(:,:,:,:) !(NCHI1,NCHI1,3,NSPIN)
      COMPLEX(8),allocatable :: sLOC1(:,:,:,:)     !(NCHI1,NCHI1,NOMEGA,NSPIN)
      COMPLEX(8),allocatable :: sLOCLAUR1(:,:,:,:) !(NCHI1,NCHI1,3,NSPIN)
      REAL(8)                :: EV
      REAL(8)                :: PHILW,sdc
      REAL(8)                :: etot
      REAL(8)                :: svar
      INTEGER(4)             :: ichi1,ichi2
      INTEGER(4)             :: nchi1
      INTEGER(4)             :: ISPIN,NU,ichi,iat
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
PRINT*,'...  READING GREENS FUNCTION done'
!
!     ==========================================================================
!     == loop over atoms                                                      ==
!     ==========================================================================
      etot=0.d0
      do iat=1,nat
        ichi1=ichibnd(1,iat)
        ichi2=ichibnd(2,iat)
        nchi1=ichi2-ichi1+1
        if(nchi1.le.0) cycle
!
!       ========================================================================
!       == collect local arrays gloc1 and uloc1                               ==
!       ========================================================================
PRINT*,'CONSTRUCTING SELF ENERGY for atom ',iat,' ...'
        allocate(uloc1(nchi1,nchi1,nchi1,nchi1))
        uloc1=uchi(ichi1:ichi2,ichi1:ichi2,ichi1:ichi2,ichi1:ichi2)
!
        allocate(gloc1(nchi1,nchi1,nomega,nspin))
        allocate(gloclaur1(nchi1,nchi1,3,nspin))
        gloc1=gloc(ichi1:ichi2,ichi1:ichi2,:,:)
        gloclaur1=gloclaur(ichi1:ichi2,ichi1:ichi2,:,:)
!
        allocate(sloc1(nchi1,nchi1,nomega,nspin))
        allocate(sloclaur1(nchi1,nchi1,3,nspin))
        CALL DMFT$HFSOLVER(NCHI1,NSPIN,NOMEGA,KBT,MU,OMEGA &
     &                    ,GLOC1,GLOCLAUR1,uloc1,PHILW,Sloc1,SloclaUR1)
        etot=etot+philw
!
!       ========================================================================
!       == map self energy onto nchi-array                                    ==
!       ========================================================================
!
!       ========================================================================
!       == clean up                                                           ==
!       ========================================================================
        deallocate(uloc1)
        deallocate(gloc1)
        deallocate(gloclaur1)
        deallocate(sloc1)
        deallocate(sloclaur1)
      enddo
!
!     ==========================================================================
!     == WRITE SELF ENERGY                                                    ==
!     ==========================================================================
!!$PRINT*,'WRITING SELF ENERGY   '
!!$      CALL DMFT_WRITESIGMA(NCHI,NOMEGA,NSPIN,KBT,MU,SIGMA,SIGLAUR,SDC) 
!!$PRINT*,'... DMFT$SOLVER DONE.'
stop 'forced: implement writing self energy before continuing!'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$SOLVER_withkset()
!     **************************************************************************
!     ** MIMICKS A DMFT SOLVER                                                **
!     **  E=1/BETA * SUM_NU E(IOMEGA_NU0+) TR[ DV * GLOC^\DAGGER]             **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: nat,Ndimd,NOMEGA,omega,kbt,mu,atomset
      IMPLICIT NONE
      real(8)    ,parameter  :: dmfthfweight=0.25d0
      REAL(8)                :: EV
      REAL(8)                :: PHILW,sdc
      REAL(8)                :: etot
      REAL(8)                :: edc
      INTEGER(4)             :: nloc  !#(local orbitals on this site)
      INTEGER(4)             :: lmnx  !#(local orbitals on this site)
      INTEGER(4)             :: iat
      INTEGER(4)             :: idimd,i,j,lmn1,lmn2,nu
      complex(8),allocatable :: h(:,:,:)
!     **************************************************************************
PRINT*,'ENTERING DMFT$SOLVER...'
      CALL CONSTANTS('EV',EV)
!
!     ==========================================================================
!     == loop over atoms                                                      ==
!     ==========================================================================
      etot=0.d0
      do iat=1,nat
        nloc=atomset(iat)%nloc
        if(nloc.le.0) cycle
!
!       ========================================================================
!       == collect local arrays gloc1 and uloc1                               ==
!       ========================================================================
PRINT*,'CONSTRUCTING SELF ENERGY for atom ',iat,' ...'
        CALL DMFT$HFSOLVER_withset(Nloc,Ndimd,NOMEGA,KBT,MU,OMEGA &
     &                            ,atomset(iat)%GLOC,atomset(iat)%GLOCLAUR &
     &                            ,atomset(iat)%u,philw &
     &                            ,atomset(iat)%sLOC,atomset(iat)%sLOCLAUR)
        etot=etot+philw
!!$write(*,fmt='(82("="),t10," sloc from hfsolver_withkset ",i4,"  ")'),iat
!!$do nu=1,nomega
!!$  do idimd=1,ndimd
!!$    do i=1,nloc
!!$      write(*,fmt='(4i5,100f10.5)')iat,nu,idimd,i &
!!$                              ,atomset(iat)%sloc(i,:,idimd,nu)
!!$    enddo
!!$  enddo
!!$enddo
!!$write(*,fmt='(82("="),t10," sloclaur from hfsolver_withkset ",i4,"  ")'),iat
!!$do nu=1,3
!!$  do idimd=1,ndimd
!!$    do i=1,nloc
!!$      write(*,fmt='(4i5,100f10.5)')iat,nu,idimd,i &
!!$                              ,atomset(iat)%sloclaur(i,:,idimd,nu)
!!$    enddo
!!$  enddo
!!$enddo
!
!       ========================================================================
!       == subtract double counting term                                      ==
!       ========================================================================
        lmnx=atomset(iat)%denmat%lmnx
        call DMFT_DC(IAT,lmnx,NDIMD,atomset(iat)%denmat%mat &
       &                       ,EDC,atomset(iat)%denmat%H)
        etot=etot+edc
!
!       ========================================================================
!       == shift static self energy to the dc term                            ==
!       ========================================================================
        nloc=atomset(iat)%nloc
        allocate(h(nloc,nloc,ndimd))
!       == SHIFT ONLY REAL PART BECAUSE ATOMSET%DENMAT%H IS REAL ===============
        H=REAL(ATOMSET(IAT)%SLOCLAUR(:,:,:,1))
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
   &                     =ATOMSET(IAT)%DENMAT%H(LMN1,LMN2,:)+REAL(H(I,J,:))
          ENDDO
        ENDDO
        DEALLOCATE(H)
!
!       ========================================================================
!       == Scale                                                              ==
!       ========================================================================
        atomset(iat)%sloc    =dmfthfweight*atomset(iat)%sloc
        atomset(iat)%sloclaur=dmfthfweight*atomset(iat)%sloclaur
        atomset(iat)%denmat%h=dmfthfweight*atomset(iat)%denmat%h
!
write(*,fmt='(82("="),t10," hamiltonian from hfsolver_withkset ",i4,"  ")'),iat
do nu=1,nomega
  do idimd=1,ndimd
    do i=1,atomset(iat)%nloc
      write(*,fmt='(4i5,100f10.5)')iat,nu,idimd,i &
                              ,atomset(iat)%sloc(i,:,idimd,nu)
    enddo
  enddo
enddo
!
      enddo ! end of loop over atoms (iat)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_DC(IAT,LMNX,NDIMD,rho,EDC,ham)
!     **************************************************************************
!     ** calculates the DFT exchange energy for the correlated orbitals.      **
!     ** edc is to be added! to the total energy (minus sign is included)     **
!     **                                                                      **
!     ** remark: some routines only work in the non-collinear mode            **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR,SBAR,SBARLI1,LNX,LOX
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      INTEGER(4),INTENT(IN)  :: LMNX          !  #(LOCAL ORBITALS ON THIS SITE)
      INTEGER(4),INTENT(IN)  :: NDIMD
      REAL(8)   ,INTENT(OUT) :: EDC           ! -e_(dft-exchange)
      real(8)   ,INTENT(IN)  :: rho(LMNX,LMNX,NDIMD) ! density matrix
      real(8)   ,INTENT(OUT) :: Ham(LMNX,LMNX,NDIMD) ! hamiltonian contrib.
      real(8)                :: d(LMNX,LMNX,4)
      real(8)                :: H(LMNX,LMNX,4)
      real(8)                :: Hall(LMNX,LMNX,4)
      REAL(8)   ,ALLOCATABLE :: DT(:,:,:)    !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: DTALL(:,:,:) !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: HT(:,:,:)    !(LMNXT,LMNXT,NDIMD)
      REAL(8)   ,ALLOCATABLE :: HTALL(:,:,:) !(LMNXT,LMNXT,NDIMD)
      REAL(8)                :: EX
      INTEGER(4)             :: IDFTTYPE
      INTEGER(4)             :: LMNXT   ! #(VALENCE+SCATTERING WAVES)
      INTEGER(4)             :: LNXT    ! #(VALENCE+SCATTERING WAVES w/o m)
      INTEGER(4),allocatable :: LoXT(:)  
      INTEGER(4)             :: lrx
      INTEGER(4)             :: LMRX    ! #(SPHERICAL HARMONICS IN THE DENSITY)
      INTEGER(4)             :: GID     ! GRID ID
      INTEGER(4)             :: NR      ! #(GRID POINTS)
      REAL(8)   ,ALLOCATABLE :: AECORE(:) !(NR) CORE DENSITY
      INTEGER(4)             :: isp,idimd,lmn
!     **************************************************************************
      ham=0.d0
      CALL DFT$GETI4('TYPE',IDFTTYPE)
      PRINT*,'IDFTTYPE ',IDFTTYPE
      IF(IDFTTYPE.eq.5002) return
!
!     ==========================================================================
!     == collect data                                                         ==
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
      allocate(loxt(lnxt))
      LoXT=POTPAR(ISP)%TAILED%Lox  
!
!     ==========================================================================
!     ==  MAP INTO NON-COLLINEAR ARRAYS
!     ==========================================================================
      D(:,:,:)=0.D0
      D(:,:,1)=RHO(:,:,1)
      IF(NDIMD.EQ.2) THEN 
        D(:,:,4)=RHO(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        D(:,:,2:4)=RHO(:,:,2:4)
      END IF
!
!     ==========================================================================
!     ==  MAP ONTO REAL ARRAY. (T,x,y,z) representation                       ==
!     ==========================================================================
      ALLOCATE(DTALL(LMNXT,LMNXT,4))
      ALLOCATE(DT(LMNXT,LMNXT,4))
      ALLOCATE(HTALL(LMNXT,LMNXT,4))
      ALLOCATE(HT(LMNXT,LMNXT,4))
      CALL LMTO_BLOWUPDENMATNL(IAT,IAT,4,LMNX,LMNX,D,LMNXT,LMNXT,DT)
      POTPAR(ISP)%TALLORB=.TRUE.
      CALL LMTO_BLOWUPDENMATNL(IAT,IAT,4,LMNX,LMNX,D,LMNXT,LMNXT,DTALL)
      CALL LMTO_SIMPLEDC(GID,NR,LMNXT,LNXT,LOXT,POTPAR(ISP)%TAILED%AEF &
     &                  ,LRX,AECORE,DT,DTALL,ex,HT,HTALL)
      CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,4,LMNXT,LMNXT,HTALL,LMNX,LMNX,Hall)
      POTPAR(ISP)%TALLORB=.FALSE.  ! DO NOT FORGET THIS!!!!!
      CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,4,LMNXT,LMNXT,HT,LMNX,LMNX,H)
      EDC=-EX
      h=-(h+hall)
PRINT*,'DOUBLE COUNTING CORRECTION ENERGY FOR ATOM=',IAT,-EX
!
!     ==========================================================================
!     ==  MAP back from  NON-COLLINEAR ARRAYS                                 ==
!     ==========================================================================
      HAM=0.D0
      HAM(:,:,1)=H(:,:,1)
      IF(NDIMD.EQ.2) THEN
         HAM(:,:,2)=H(:,:,4)
      ELSE IF(NDIMD.EQ.4) THEN
         HAM(:,:,2:4)=H(:,:,2:4)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_Uchi()
!     **************************************************************************
!     ** calculate the u-tensor of the correlated orbitals in the selected set
!     **************************************************************************
      use lmto_module, only: ispecies,lnx,lox,potpar
      use dmft_module, only: nat,nchi,uchi
      implicit none
      real(8)  ,allocatable :: ub(:,:,:,:)
      integer(4)            :: isp ! atom type
      integer(4)            :: l,ln,lmn,ichi,ichi1,ichi2
      integer(4)            :: lmnx
      integer(4)            :: nchi1        ! #(selected orbitals on this atom)
      integer(4)            :: iat
!     **************************************************************************
                                          call trace$push('dmft_uchi')
      IF(.NOT.ALLOCATED(UCHI)) THEN
        CALL ERROR$MSG('UCHI IN DMFT_MODULE NOT ALLOCATED')
        CALL ERROR$STOP('DMFT_UCHI')
      END IF
      uchi(:,:,:,:)=0.d0
!
      ichi2=0
      ichi1=1
      do iat=1,nat
        isp=ispecies(iat)
!
!       == calculate local dimensions lmnx and nchi1 ===========================
        lmn=0
        ichi=0
        do ln=1,lnx(isp)
          l=lox(ln,isp)
          if(potpar(isp)%torb(ln)) ichi=ichi+2*l+1
          lmn=lmn+2*l+1
        enddo
        lmnx=lmn
        nchi1=ichi
        if(nchi1.eq.0) cycle
!
!       == matrix will be inserted in the section ichi1:ichi2 of 1:nchi ========
        ichi1=ichi2+1
        ichi2=ichi2+nchi1
        IF(ICHI2.GT.NCHI) THEN
          CALL ERROR$MSG('ARRAY SIZE EXCEEDED')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$I4VAL('NCHI',NCHI)
          CALL ERROR$I4VAL('ICHI2',ICHI2)
          CALL ERROR$MSG('NEEDS EXPERT HELP')
          CALL ERROR$STOP('DMFT_UCHI')
        END IF
!
!       ========================================================================
!       ==  calculate U-tensor in the basis of local orbitals ignoring torb   ==
!       ========================================================================
        allocate(ub(lmnx,lmnx,lmnx,lmnx))
        call DMFT_ULOCAL(IAT,LMNX,Ub)
!
!       ========================================================================
!       ==  throw out elements that are not used                              ==
!       ========================================================================
        lmn=0
        ichi=0
        do ln=1,lnx(isp)
          l=lox(ln,isp)
          if(potpar(isp)%torb(ln)) then
            ub(:,:,:,ichi+1:ichi+2*l+1)=ub(:,:,:,lmn+1:lmn+2*l+1)
            ichi=ichi+2*l+1
          end if
          lmn=lmn+2*l+1
        enddo
!  
!       ========================================================================
        lmn=0
        ichi=0
        do ln=1,lnx(isp)
          l=lox(ln,isp)
          if(potpar(isp)%torb(ln)) then
            ub(:,:,ichi+1:ichi+2*l+1,:nchi1)=ub(:,:,lmn+1:lmn+2*l+1,:nchi1)
            ichi=ichi+2*l+1
          end if
          lmn=lmn+2*l+1
        enddo
!  
!       ========================================================================
        lmn=0
        ichi=0
        do ln=1,lnx(isp)
          l=lox(ln,isp)
          if(potpar(isp)%torb(ln)) then
            ub(:,ichi+1:ichi+2*l+1,:nchi1,:nchi1) &
      &                                     =ub(:,lmn+1:lmn+2*l+1,:nchi1,:nchi1)
            ichi=ichi+2*l+1
          end if
          lmn=lmn+2*l+1
        enddo
!  
!       ========================================================================
        lmn=0
        ichi=0
        do ln=1,lnx(isp)
          l=lox(ln,isp)
          if(potpar(isp)%torb(ln)) then
            ub(ichi+1:ichi+2*l+1,:nchi1,:nchi1,:nchi1) &
       &                               =ub(lmn+1:lmn+2*l+1,:nchi1,:nchi1,:nchi1)
            ichi=ichi+2*l+1
          end if
          lmn=lmn+2*l+1
        enddo  
!
!       ========================================================================
!       ==  throw out elements that are not used                              ==
!       ========================================================================
        uchi(ichi1:ichi2,ichi1:ichi2,ichi1:ichi2,ichi1:ichi2) &
     &                                          =ub(:nchi1,:nchi1,:nchi1,:nchi1)
        deallocate(ub)
      enddo  ! end loop over atoms
                                          call trace$pop()
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_Uchi_withatomset()
!     **************************************************************************
!     ** calculate the u-tensor of the correlated orbitals in the selected set
!     **************************************************************************
      use lmto_module, only: ispecies,lnx,lox,potpar
      use dmft_module, only: nat,nchi,uchi,atomset
      implicit none
      real(8)  ,allocatable :: ub(:,:,:,:)
      integer(4)            :: isp ! atom type
      integer(4)            :: l,ln,lmn,ichi
      integer(4)            :: lmnx
      integer(4)            :: nchi1        ! #(selected orbitals on this atom)
      integer(4)            :: iat
!     **************************************************************************
                                          call trace$push('dmft_uchi')
      do iat=1,nat
        isp=ispecies(iat)
!
!       == calculate local dimensions lmnx and nchi1 ===========================
        lmn=0
        ichi=0
        do ln=1,lnx(isp)
          l=lox(ln,isp)
          if(potpar(isp)%torb(ln)) ichi=ichi+2*l+1
          lmn=lmn+2*l+1
        enddo
        lmnx=lmn
        nchi1=ichi
        if(nchi1.eq.0) cycle
        IF(NCHI1.NE.ATOMset(IAT)%NLOC) THEN
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
!       ==  throw out elements that are not used                              ==
!       ========================================================================
        lmn=0
        ichi=0
        do ln=1,lnx(isp)
          l=lox(ln,isp)
          if(potpar(isp)%torb(ln)) then
            ub(:,:,:,ichi+1:ichi+2*l+1)=ub(:,:,:,lmn+1:lmn+2*l+1)
            ichi=ichi+2*l+1
          end if
          lmn=lmn+2*l+1
        enddo
!  
!       ========================================================================
        lmn=0
        ichi=0
        do ln=1,lnx(isp)
          l=lox(ln,isp)
          if(potpar(isp)%torb(ln)) then
            ub(:,:,ichi+1:ichi+2*l+1,:nchi1)=ub(:,:,lmn+1:lmn+2*l+1,:nchi1)
            ichi=ichi+2*l+1
          end if
          lmn=lmn+2*l+1
        enddo
!  
!       ========================================================================
        lmn=0
        ichi=0
        do ln=1,lnx(isp)
          l=lox(ln,isp)
          if(potpar(isp)%torb(ln)) then
            ub(:,ichi+1:ichi+2*l+1,:nchi1,:nchi1) &
      &                                     =ub(:,lmn+1:lmn+2*l+1,:nchi1,:nchi1)
            ichi=ichi+2*l+1
          end if
          lmn=lmn+2*l+1
        enddo
!  
!       ========================================================================
        lmn=0
        ichi=0
        do ln=1,lnx(isp)
          l=lox(ln,isp)
          if(potpar(isp)%torb(ln)) then
            ub(ichi+1:ichi+2*l+1,:nchi1,:nchi1,:nchi1) &
       &                               =ub(lmn+1:lmn+2*l+1,:nchi1,:nchi1,:nchi1)
            ichi=ichi+2*l+1
          end if
          lmn=lmn+2*l+1
        enddo  
!
!       ========================================================================
!       ==  put onto atomset                                                  ==
!       ========================================================================
        atomset(iat)%u=ub(:nchi1,:nchi1,:nchi1,:nchi1)
        deallocate(ub)
      enddo  ! end loop over atoms
                                          call trace$pop()
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_ULOCAL(IAT,LMNX,U)
!     **************************************************************************
!     **  calculates the U-tensor of the natural tight-binding orbitals       **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR,SBAR,sbarli1,LNX,LOX
      implicit none
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
      INTEGER(4)             :: LMN,L,LM,LN,IM,NN,NN0,j
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
      subroutine dmft_EX(iat,lmnx,ndimd,d,e,h)
      use lmto_module, only : potpar,ispecies
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
      isp=ispecies(iat)
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
     &                        ,GLOC,GLOCLAUR,uchi,PHILW,SIGMA,SIGLAUR)
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
      REAL(8)   ,INTENT(OUT) :: uchi(nchi,nchi,nchi,nchi)
      REAL(8)   ,INTENT(OUT) :: PHILW
      COMPLEX(8),INTENT(OUT) :: SIGMA(NCHI,NCHI,NOMEGA,NSPIN)
      COMPLEX(8),INTENT(OUT) :: SIGLAUR(NCHI,NCHI,3,NSPIN)
      logical(4),parameter   :: tprint=.true.
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
                                     call trace$push('DMFT$HFSOLVER')
!
!     ==========================================================================
!     == CALCULATE ONE-PARTICLE DENSITY MATRIX                                ==
!     ==========================================================================
print*,'before dmft_rho'
      CALL DMFT_RHO(NCHI,NOMEGA,NSPIN,KBT,OMEGA,GLOC,GLOCLAUR,RHO)
print*,'after dmft_rho'
!
       if(tprint) then
          WRITE(*,FMT='(82("="),T10," DENSITY MATRIX from gloc IN HFSOLVER")')
          DO ISPIN=1,NSPIN
            DO I=1,nchi
              WRITE(*,FMT='("RHO(IS=",I1,"):",100("(",2F10.5,")"))') &
    &                      ISPIN,RHO(I,:,ISPIN)
            ENDDO
          ENDDO
       end if
!
!     ==========================================================================
!     == CALCULATE HARTREE FOCK POTENTIAL                                     ==
!     ==========================================================================
IF(NSPIN.NE.2) THEN
 CALL ERROR$MSG('IMPLEMENTATION VAID ONLY FOR NSPIN=2')
 CALL ERROR$I4VAL('NSPIN',NSPIN)
 CALL ERROR$STOP('DMFT$HFSOLVER')
END IF
      estat=0.d0
      DO ISPIN=1,NSPIN
        RHO1=RHO(:,:,ISPIN)
        IF(NSPIN.EQ.1) RHO1=0.5D0*RHO1
        H=(0.D0,0.D0)
        DO I=1,NCHI
          DO J=1,NCHI
            DO K=1,NCHI
              DO L=1,NCHI
                SVAR=0.5D0*Uchi(I,J,K,L) 
                if(svar.eq.0.d0) cycle
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
!
        IF(NSPIN.EQ.1) ESTAT=2.D0*ESTAT
      ENDDO
print*,'estat ',estat
!!$print*,'stopping now in hfsolver'
!!$stop 'forced'
!
                                     call trace$pop()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$HFSOLVER_withset(Nloc,ndimd,NOMEGA,KBT,MU,OMEGA &
     &                        ,GLOC,GLOCLAUR,uchi,PHILW,Sloc,SlocLAUR)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NLOC
      INTEGER(4),INTENT(IN)  :: Ndimd
      INTEGER(4),INTENT(IN)  :: NOMEGA
      REAL(8)   ,INTENT(IN)  :: KBT
      REAL(8)   ,INTENT(IN)  :: MU
      REAL(8)   ,INTENT(IN)  :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN)  :: GLOC(NLOC,NLOC,Ndimd,nomega)
      COMPLEX(8),INTENT(IN)  :: GLOCLAUR(NLOC,NLOC,Ndimd,3)
      REAL(8)   ,INTENT(OUT) :: uchi(nloc,nloc,nloc,nloc)
      REAL(8)   ,INTENT(OUT) :: PHILW
      COMPLEX(8),INTENT(OUT) :: Sloc(NLOC,NLOC,Ndimd,nomega)
      COMPLEX(8),INTENT(OUT) :: SlocLAUR(NLOC,NLOC,Ndimd,3)
      logical(4),parameter   :: tprint=.true.
      COMPLEX(8)             :: RHO(NLOC,NLOC,Ndimd)
      COMPLEX(8)             :: H(NLOC,NLOC,ndimd)
      REAL(8)                :: ESTAT
      INTEGER(4)             :: NU,Idimd,I,J,K,L
      REAL(8)                :: SVAR
!     **************************************************************************
                                     call trace$push('DMFT$HFSOLVER_withset')
!
!     ==========================================================================
!     == CALCULATE ONE-PARTICLE DENSITY MATRIX                                ==
!     ==========================================================================
      call DMFT_RHO_withset(Nloc,ndimd,NOMEGA,KBT,OMEGA,Gloc,GlocLAUR,RHO)
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
      estat=0.d0
      H=(0.D0,0.D0)
      DO I=1,NLOC
        DO J=1,NLOC
          DO K=1,NLOC
            DO L=1,NLOC
              SVAR=-0.25D0*Uchi(I,J,K,L) 
              if(svar.eq.0.d0) cycle
              ESTAT=ESTAT+SVAR*REAL(sum(RHO(K,J,:)*RHO(l,i,:)))
              H(k,j,:)=H(k,j,:)+SVAR*RHO(i,l,:)
              H(I,L,:)=H(I,L,:)+SVAR*RHO(K,J,:)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO NU=1,NOMEGA
        Sloc(:,:,:,nu)=H(:,:,:)
      ENDDO
      SlocLAUR(:,:,:,1)=H
      SlocLAUR(:,:,:,2)=(0.D0,0.D0)
      SlocLAUR(:,:,:,3)=(0.D0,0.D0)
print*,'estat ',estat
!
      IF(TPRINT) THEN
        CALL SPINOR$PRINTMATRIX(6 &
     &    ,'h_xc (UPDOWN) IN DMFT$HFSOLVER_WITHSET','UPDOWN',1,NLOC &
     &                       ,NDIMD,NLOC,h)
      END IF
!
                                     call trace$pop()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_GETSIGMA()
!     **************************************************************************
!     ** PIPSI <PI_A|PSI_N> IS THE PRE-FACTOR OF LOCAL ORBITAL |CHI_A> IN     **
!     **       THE LOCAL ORBITAL EXPANSION OF |PSI_N>                         **
!     **************************************************************************
      USE DMFT_MODULE, ONLY: TON,NCHI,NSPIN,NOMEGA,OMEGA,KBT,MU &
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
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NSPIN,NOMEGA,OMEGA,KBT,MU &
     &                      ,AMIX,SMAT,RHOOFK
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
      INTEGER(4)               :: IKPT,ISPIN,ITER,I,J,K,L,IND1,IND2
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_CONSTRAINTS')
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
PRINT*,'minMAX(DH0) old ',ispin,MAXVAL(abs(DH0(:,:))),MINVAL(abs(DH0(:,:)))
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
      SUBROUTINE DMFT_CONSTRAINTS_withkset(type)
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
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NSPIN,ndimd,nat,NOMEGA &
     &                      ,OMEGA,KBT,MU &
     &                      ,AMIX,kset,atomset,ichibnd
      IMPLICIT NONE
      character(*),intent(in)  :: type  ! can be 'hrho','h0'
      REAL(8)   ,PARAMETER     :: TOL=1.D-6
      INTEGER(4),PARAMETER     :: NITER=1000
      LOGICAL(4),PARAMETER     :: TMIX=.FALSE.
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      COMPLEX(8)               :: MAT(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: MAT3(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: G(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: GLAUR1(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: GLAUR2(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: GLAUR3(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: s(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: sLAUR1(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: sLAUR2(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: sLAUR3(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: DEVRHO(NCHI,NCHI,ndimd)
      COMPLEX(8),allocatable   :: X4(:,:,:)
      COMPLEX(8)               :: DH0(NCHI,NCHI,ndimd)
      LOGICAL(4)               :: CONVG
      REAL(8)                  :: MAXDEV
      COMPLEX(8)               :: CSVAR
      INTEGER(4)               :: IKPT,ISPIN,ITER,idimd,nu,iat
      INTEGER(4)               :: i1,i2
      logical(4)               :: th0
      integer(4)               :: lx4 ! first dimension of x4
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_CONSTRAINTS_withkset')
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
!     == allocate x4 ===========================================================
      if(ndimd.le.2) then
        lx4=nchi**2
      else if(ndimd.eq.4) then
        lx4=(2*nchi)**2
      endif
      allocate(x4(lx4,lx4,nspin))
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
          slaur1=(0.d0,0.d0)
          slaur2=(0.d0,0.d0)
          slaur3=(0.d0,0.d0)
          if(th0) then
            do iat=1,nat
              i1=ichibnd(1,iat)
              i2=ichibnd(2,iat)
              slaur1(i1:i2,i1:i2,:)=atomset(iat)%sloclaur(:,:,:,1)
              slaur2(i1:i2,i1:i2,:)=atomset(iat)%sloclaur(:,:,:,2)
              slaur3(i1:i2,i1:i2,:)=atomset(iat)%sloclaur(:,:,:,3)
            enddo
          end if
          GLAUR1=KSET(IKPT)%SINV
!
          MAT=-MU*KSET(IKPT)%SMAT
          MAT=MAT+KSET(IKPT)%H0
          IF(TH0) MAT=MAT+Slaur1(:,:,:)
!         == glaur2=sinv*mat*sinv ==============================================
          call spinor$matmul(NDIMD,Nchi,mat,KSET(IKPT)%SINV,mat2)
          call spinor$matmul(NDIMD,Nchi,KSET(IKPT)%SINV,mat2,glaur2)
!
          MAT=-MU*KSET(IKPT)%SMAT
          MAT=mat+KSET(IKPT)%H0
          if(th0) MAT=mat+SLAUR1(:,:,:)
!         == mat3=mat*(sinv*mat) ===============================================
          call spinor$matmul(NDIMD,Nchi,mat,kset(ikpt)%sinv,mat2)
          call spinor$matmul(NDIMD,Nchi,mat,mat2,mat3)
          if(th0) MAT3=MAT3+SLAUR2(:,:,:)
!         == glaur3=sinv*(mat3*sinv) ===========================================
          call spinor$matmul(NDIMD,Nchi,mat3,kset(ikpt)%sinv,mat)
          call spinor$matmul(NDIMD,Nchi,kset(ikpt)%sinv,mat,glaur3)
!
          X4=(0.D0,0.D0)        
          DEVRHO=(0.D0,0.D0)
          DO NU=1,NOMEGA
        
!           == CONSTRUCT LATTICE GREENS FUNCTION =============================
            MAT=(CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT
            mat=mat-kset(ikpt)%H0
            if(th0) then
              do iat=1,nat
                i1=ichibnd(1,iat)
                i2=ichibnd(2,iat)
                mat(i1:i2,i1:i2,:)=mat(i1:i2,i1:i2,:) &
     &                            -atomset(iat)%sloc(:,:,:,nu)
              enddo
            end if
            call SPINOR$INVERT(NDIMD,NCHI,mat,g)
!           == subtract Laurent expansion to improve omega-convergence =========
            CSVAR=1.D0/(CI*OMEGA(NU))
            DEVRHO=DEVRHO+KBT*(G-CSVAR*(GLAUR1+CSVAR*(GLAUR2+CSVAR*GLAUR3)))
!
!           == ACCUMULATE SECOND ORDER TERM ====================================
            CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,g) ! convert to updown
            call DMFT_addx4(ndimd,nchi,g,lx4,nspin,x4)
          ENDDO ! END LOOP OVER MATSUBARA FREQUENCIES
          x4=x4*kbt
!
!         == INCLUDE NEGATIVE FREQUENCIES ======================================
!         == ALREADY INCLUDED FOR X4
          do idimd=1,ndimd
            DEVRHO(:,:,idimd)=DEVRHO(:,:,idimd) &
     &                       +CONJG(TRANSPOSE(DEVRHO(:,:,idimd)))
          enddo
!         == ADD TAILS (GLAUR3 DOES NOT CONTRIBUTE) ============================
          DEVRHO=DEVRHO+0.5D0*GLAUR1-0.25D0/KBT*GLAUR2
!
!         == SUBTRACT TARGET DENSITY =========================================
          DEVRHO=DEVRHO-kset(ikpt)%rho
!
!!$          MAXDEV=MAX(MAXDEV,MAXVAL(ABS(DEVRHO)))

CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,DEVRHO) ! CONVERT TO UPDOWN
MAXDEV=MAX(MAXDEV,MAXVAL(ABS(DEVRHO)))
DO IDIMD=1,NDIMD
    PRINT*,'ITER=',ITER,' IKPT=',IKPT,' IDIMD=',IDIMD &
  &       ,' MAXVAL(X4)  ',MAXVAL(ABS(X4(:,:,Idimd))) &
  &       ,' MAX(DEVRHO) ',MAXVAL(ABS(DEVRHO(:,:,IDIMD)))
ENDDO
CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,DEVRHO) ! CONVERT TO txyz
!
          IF(TMIX) THEN
            DH0=AMIX*DEVRHO
          ELSE
            CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,DEVRHO) ! CONVERT TO UPDOWN
!           == devrho+x4*dh0=0 -> dh0
            call DMFT_solvex4(ndimd,nchi,devrho,lx4,nspin,x4,dh0)
            CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,dh0)    ! CONVERT TO txyz
          END IF
!
!         =================================================================== 
!         == ADD DH
!         ====================================================================
          kset(ikpt)%H0=kset(ikpt)%H0+DH0
        ENDDO ! end of loop over k-points
!
!       ========================================================================
!       == CHECK CONVERGENCE
!       ========================================================================
        PRINT*,'MAXDEV',ITER,MAXDEV
        CONVG=MAXDEV.LT.TOL
        IF(CONVG) EXIT
      ENDDO   !end of loop over iterations
      IF(.NOT.CONVG) THEN
        CALL ERROR$MSG('LOOP NOT CONVERGED')
        CALL ERROR$STOP('DMFT_CONSTRAINTS')
      END IF
!
!     ==========================================================================
!     == map onto hrho                                                        ==
!     ==========================================================================
      if(.not.th0) then
        do ikpt=1,nkptl
          kset(ikpt)%hrho=kset(ikpt)%h0
        enddo
      end if
print*,'th0 ',th0

                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_solvex4(ndimd,nchi,drho,lx4,nspin,x4,dh)
!     **************************************************************************
!     ** drho+x4*dh=0
!     **************************************************************************
      implicit none
      integer(4),intent(in)    :: ndimd
      integer(4),intent(in)    :: nchi
      integer(4),intent(in)    :: lx4
      integer(4),intent(in)    :: nspin
      complex(8),intent(in)    :: drho(nchi,nchi,ndimd)
      complex(8),intent(in)    :: x4(lx4,lx4,nspin)
      complex(8),intent(out)   :: dh(nchi,nchi,ndimd)
      complex(8)               :: bvec(lx4)
      complex(8)               :: xvec(lx4)
      integer(4)               :: ia,ib
      integer(4)               :: ida,idb
      integer(4)               :: ichia,ichib
      integer(4)               :: iab
      integer(4)               :: idimab
      integer(4)               :: ispin 
!     **************************************************************************
!
!     ==========================================================================
!     == collinear case: spin indices are decoupled                           ==
!     ==========================================================================
      if(ndimd.le.2) then
        if(lx4.ne.nchi**2.or.nspin.ne.ndimd) then
          call error$msg('inconsistent array dimensions')
          call error$i4val('ndimd',ndimd)
          call error$i4val('nspin',nspin)
          call error$i4val('nchi',nchi)
          call error$i4val('lx4',lx4)
          call error$stop('dmft_addx4')
        end if
        do ispin=1,nspin
          do ia=1,nchi
            do ib=1,nchi
              iab=ia+nchi*(ib-1)
              bvec(iab)=drho(ia,ib,ispin)
            enddo
          enddo              
!         == a*x=b -> x (Therefore minus sign needed)===========================
          call lib$matrixsolvec8(lx4,lx4,1,x4(:,:,ispin),xvec,bvec)
          do ia=1,nchi
            do ib=1,nchi
              iab=ia+nchi*(ib-1)
              dh(ia,ib,ispin)=-xvec(iab)
            enddo
          enddo              
        enddo
!
!     ==========================================================================
!     == non-collinear case:                                                  ==
!     ==========================================================================
      else if(ndimd.eq.4) then
        if(lx4.ne.(2*nchi)**2.or.nspin.ne.1) then
          call error$msg('inconsistent array dimensions')
          call error$i4val('ndimd',ndimd)
          call error$i4val('nspin',nspin)
          call error$i4val('nchi',nchi)
          call error$i4val('lx4',lx4)
          call error$stop('dmft_addx4')
        end if
        do ida=1,2
          do ichia=1,nchi
            ia=ichia+nchi*(ida-1)
            do idb=1,2
              do ichib=1,nchi
                ib=ichib+nchi*(idb-1)
                iab=ia+2*nchi*(ib-1)
                idimab=ida+2*(idb-1)
                bvec(iab)=drho(ichia,ichib,idimab)
              enddo
            enddo              
          enddo              
        enddo              
!         == a*x=b -> x (Therefore minus sign needed)===========================
        call lib$matrixsolvec8(lx4,lx4,1,x4,xvec,bvec)
        do ida=1,2
          do ichia=1,nchi
            ia=ichia+nchi*(ida-1)
            do idb=1,2
              do ichib=1,nchi
                ib=ichib+nchi*(idb-1)
                iab=ia+2*nchi*(ib-1)
                idimab=ida+2*(idb-1)
                dh(ichia,ichib,idimab)=-xvec(iab)
              enddo
            enddo              
          enddo              
        enddo              
      end if
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_addx4(ndimd,nchi,g,lx4,nspin,x4)
!     **************************************************************************
      implicit none
      integer(4),intent(in)    :: ndimd
      integer(4),intent(in)    :: nchi
      integer(4),intent(in)    :: lx4
      integer(4),intent(in)    :: nspin
      complex(8),intent(in)    :: g(nchi,nchi,ndimd)
      complex(8),intent(inout) :: x4(lx4,lx4,nspin)
      integer(4)               :: ia,ib,ic,id
      integer(4)               :: ida,idb,idc,idd
      integer(4)               :: ichia,ichib,ichic,ichid
      integer(4)               :: iab,icd
      integer(4)               :: idimac,idimdb,idimca,idimbd
      integer(4)               :: ispin 
!     **************************************************************************
!
!     ==========================================================================
!     == collinear case: spin indices are decoupled                           ==
!     ==========================================================================
      if(ndimd.le.2) then
        if(lx4.ne.nchi**2.or.nspin.ne.ndimd) then
          call error$msg('inconsistent array dimensions')
          call error$i4val('ndimd',ndimd)
          call error$i4val('nspin',nspin)
          call error$i4val('nchi',nchi)
          call error$i4val('lx4',lx4)
          call error$stop('dmft_addx4')
        end if
        do ispin=1,ndimd
          iab=0
          do ia=1,nchi
            do ib=1,nchi
              iab=iab+1
              icd=0
              do ic=1,nchi
                do id=1,nchi
                  icd=icd+1
!                 == the second term is for -i*omega_\nu =======================
!                 == g(-i*omega)=g^\dagger(i*omega) ============================
                  x4(iab,icd,ispin)=x4(iab,icd,ispin) &
     &                             +g(ia,ic,ispin)*g(id,ib,ispin) &
     &                             +conjg(g(ic,ia,ispin)*g(ib,id,ispin))
                enddo
              enddo
            enddo
          enddo
        enddo
!
!     ==========================================================================
!     == non-collinear case:                                                  ==
!     ==========================================================================
      else
        if(lx4.ne.(2*nchi)**2.or.nspin.ne.1) then
          call error$msg('inconsistent array dimensions')
          call error$i4val('ndimd',ndimd)
          call error$i4val('nspin',nspin)
          call error$i4val('nchi',nchi)
          call error$i4val('lx4',lx4)
          call error$stop('dmft_addx4')
        end if
        iab=0
        do ida=1,2
        do ichia=1,nchi
          do idb=1,2
          do ichib=1,nchi
            iab=iab+1
            icd=0
            do idc=1,2
            do ichic=1,nchi
              do idd=1,2
              do ichid=1,nchi
                icd=icd+1
!
                idimac=ida+2*(idc-1)
                idimca=idc+2*(ida-1)
                idimdb=idd+2*(idb-1)
                idimbd=idb+2*(idd-1)
                x4(iab,icd,1)=x4(iab,icd,1) &
     &                       +g(ichia,ichic,idimac)*g(ichid,ichib,idimdb) &
     &                 +conjg(g(ichic,ichia,idimca)*g(ichib,ichid,idimbd))
              enddo
              enddo
            enddo
            enddo
          enddo
          enddo
        enddo
        enddo
      end if      
      return
      end
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
      SUBROUTINE DMFT_RHO_withset(NCHI,ndimd,NOMEGA,KBT,OMEGA,G,GLAUR,RHO)
!     **************************************************************************
!     **  CALCULATES THE DENSITY MATRIX FROM A GREENS FUNCTION                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: NOMEGA
      INTEGER(4),INTENT(IN) :: Ndimd
      REAL(8)   ,INTENT(IN) :: KBT
      REAL(8)   ,INTENT(IN) :: OMEGA(NOMEGA)
      COMPLEX(8),INTENT(IN) :: G(NCHI,NCHI,Ndimd,nomega)
      COMPLEX(8),INTENT(IN) :: GLAUR(NCHI,NCHI,Ndimd,3)
      COMPLEX(8),INTENT(OUT):: RHO(NCHI,NCHI,Ndimd)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)  ! SQRT(-1)
      INTEGER(4)            :: Idimd,NU
      COMPLEX(8)            :: CSVAR
!     **************************************************************************
      RHO=(0.D0,0.D0)
      DO NU=1,NOMEGA
        CSVAR=1.D0/(CI*OMEGA(NU))
        RHO=RHO+G(:,:,:,nu)-GLAUR(:,:,:,1)*CSVAR &
     &                     -GLAUR(:,:,:,2)*CSVAR**2 &
     &                     -GLAUR(:,:,:,3)*CSVAR**3
      ENDDO
      do idimd=1,ndimd
        RHO(:,:,Idimd)=RHO(:,:,Idimd)+CONJG(TRANSPOSE(RHO(:,:,Idimd)))
      enddo
      RHO=RHO*KBT
      RHO=RHO+0.5D0*GLAUR(:,:,:,1)-0.25D0/KBT*GLAUR(:,:,:,2)
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
      USE DMFT_MODULE, ONLY: TON,NCHI,NB,NKPTL,NSPIN,NOMEGA,OMEGA,KBT,MU &
     &                      ,PIPSI,ERHO,WKPTL
      IMPLICIT NONE
      LOGICAL(4),PARAMETER     :: TTEST=.false.
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
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NSPIN,NOMEGA,OMEGA,KBT,MU &
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
      INTEGER(4)               :: ISPIN,NU,I
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
      USE DMFT_MODULE, ONLY: TON,NCHI,NKPTL,NSPIN,NOMEGA,OMEGA,KBT,MU &
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
      SUBROUTINE DMFT_GLOC_withatomset()
!     **************************************************************************
!     **************************************************************************
      USE DMFT_MODULE, ONLY: NCHI,NKPTL,NSPIN,ndimd,nat,NOMEGA,OMEGA,KBT,MU &
     &                       ,SMAT,SINV,kset,atomset,ichibnd
      IMPLICIT NONE
      COMPLEX(8),PARAMETER     :: CI=(0.D0,1.D0)  ! SQRT(-1)
      LOGICAL(4),PARAMETER     :: TPRINT=.FALSE.
      LOGICAL(4),PARAMETER     :: TTEST=.false.
      COMPLEX(8)               :: MAT(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: MAT2(NCHI,NCHI,ndimd)
      COMPLEX(8)               :: MAT3(NCHI,NCHI,ndimd)
      complex(8)               :: glaur1(nchi,nchi,ndimd)
      complex(8)               :: glaur2(nchi,nchi,ndimd)
      complex(8)               :: glaur3(nchi,nchi,ndimd)
      complex(8)               :: slaur1(nchi,nchi,ndimd)
      complex(8)               :: slaur2(nchi,nchi,ndimd)
      complex(8)               :: slaur3(nchi,nchi,ndimd)
      COMPLEX(8)               :: G(NCHI,NCHI,ndimd)
      real(8)                  :: wkptl
      INTEGER(4)               :: IKPT,NU,Iat,i1,i2,i,idimd
!     **************************************************************************
                              CALL TRACE$PUSH('DMFT_GLOC')
      do iat=1,nat
        atomset(iat)%gloc=(0.d0,0.d0)
        atomset(iat)%gloclaur=(0.d0,0.d0)
      enddo
!
      DO IKPT=1,NKPTL
        wkptl=kset(ikpt)%wkptl
!       == LAURENT EXPANSION FOR THE GREENS FUNCTION =========================
        slaur1=(0.d0,0.d0)
        slaur2=(0.d0,0.d0)
        slaur3=(0.d0,0.d0)
        do iat=1,nat
          i1=ichibnd(1,iat)
          i2=ichibnd(2,iat)
          slaur1(i1:i2,i1:i2,:)=atomset(iat)%sloclaur(:,:,:,1)
          slaur2(i1:i2,i1:i2,:)=atomset(iat)%sloclaur(:,:,:,2)
          slaur3(i1:i2,i1:i2,:)=atomset(iat)%sloclaur(:,:,:,3)
        enddo
        GLAUR1=KSET(IKPT)%SINV
!
        MAT=-MU*KSET(IKPT)%SMAT
        MAT=MAT+KSET(IKPT)%H0
        MAT=MAT+SLAUR1
!       == glaur2=sinv*mat*sinv ==============================================
        call spinor$matmul(NDIMD,Nchi,mat,KSET(IKPT)%SINV,mat2)
        call spinor$matmul(NDIMD,Nchi,KSET(IKPT)%SINV,mat2,glaur2)
!
        MAT=-MU*KSET(IKPT)%SMAT
        MAT=mat+KSET(IKPT)%H0
        MAT=mat+Slaur1
!       == mat3=mat*(sinv*mat) ===============================================
        call spinor$matmul(NDIMD,Nchi,mat,kset(ikpt)%sinv,mat2)
        call spinor$matmul(NDIMD,Nchi,mat,mat2,mat3)
        MAT3=MAT3+Slaur2
!       == glaur3=sinv*(mat3*sinv) ===========================================
        call spinor$matmul(NDIMD,Nchi,mat3,kset(ikpt)%sinv,mat)
        call spinor$matmul(NDIMD,Nchi,kset(ikpt)%sinv,mat,glaur3)
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
        do iat=1,nat
          i1=ichibnd(1,iat)
          i2=ichibnd(2,iat)
          atomset(iat)%gloclaur(:,:,:,1)=atomset(iat)%gloclaur(:,:,:,1) &
     &                                  +wkptl*glaur1(i1:i2,i1:i2,:)
          atomset(iat)%gloclaur(:,:,:,2)=atomset(iat)%gloclaur(:,:,:,2) &
     &                                  +wkptl*glaur2(i1:i2,i1:i2,:)
          atomset(iat)%gloclaur(:,:,:,3)=atomset(iat)%gloclaur(:,:,:,3) &
     &                                  +wkptl*glaur3(i1:i2,i1:i2,:)
        enddo
!
        DO NU=1,NOMEGA
!         == CONSTRUCT LATTICE GREENS FUNCTION =================================
          MAT=(CI*OMEGA(NU)+MU)*KSET(IKPT)%SMAT
          MAT=MAT-KSET(IKPT)%H0
          DO IAT=1,NAT
            I1=ICHIBND(1,IAT)
            I2=ICHIBND(2,IAT)
            MAT(I1:I2,I1:I2,:)=MAT(I1:I2,I1:I2,:)-ATOMSET(IAT)%SLOC(:,:,:,NU)
          ENDDO
          CALL SPINOR$INVERT(NDIMD,NCHI,MAT,G)
!         == ACCOUNT FOR K-POINT (-K) ==========================================
          IF(KSET(IKPT)%TADDMINUSK) THEN
            do idimd=1,ndimd
              MAT(:,:,idimd)=2.D0*CI*OMEGA(NU) &
      &                          *TRANSPOSE(CONJG(KSET(IKPT)%SMAT(:,:,idimd))) &
      &                          -TRANSPOSE(CONJG(MAT(:,:,idimd)))
            enddo
            CALL SPINOR$INVERT(NDIMD,NCHI,MAT,MAT2)
            G=0.5D0*(G+MAT2)
          END IF
          DO IAT=1,NAT
            I1=ICHIBND(1,IAT)
            I2=ICHIBND(2,IAT)
            ATOMSET(IAT)%GLOC(:,:,:,NU)=ATOMSET(IAT)%GLOC(:,:,:,NU) &
     &                                 +WKPTL*G(I1:I2,I1:I2,:)
          ENDDO
!
        ENDDO ! END OF LOOP OVER MATSUBARA FREQUENCIES
      ENDDO   ! END OF LOOP OVER KPOINTS
!
!     ==========================================================================
!     == take k-points into account that have been left out due to            ==
!     == time inversion symmetry                                              ==
!     ==========================================================================
!     == all Laurent expansion coefficients are hermitian ======================
      do iat=1,nat
        I1=ICHIBND(1,IAT)
        I2=ICHIBND(2,IAT)
        do i=1,3
          do idimd=1,ndimd
            atomset(iat)%gloclaur(:,:,idimd,i)=0.5d0 &
     &                   *(atomset(iat)%gloclaur(:,:,idimd,i) &
     &                    +conjg(transpose(atomset(iat)%gloclaur(:,:,idimd,i))))
          enddo
        enddo
      enddo

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
      INTEGER(4)             :: IKPT,ISPIN,ibh,ichi,ipro
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT$ADDTOHPSI_withkset()
!     **************************************************************************
!     ** 
!     ** 
!     ** 
!     ** 
!     ** 
!     **************************************************************************
      USE DMFT_MODULE, ONLY  : TON,NKPTL,NSPIN,ndim,ndimd,NB,nchi,nat &
     &                        ,iproofchi,ichibnd,atomset,kset
      USE WAVES_MODULE, ONLY : GSET,WAVES_SELECTWV,THIS,map
      IMPLICIT NONE
      complex(8),parameter   :: ci=(0.d0,1.d0)
      COMPLEX(8),ALLOCATABLE :: dhpipsi(:,:,:) !(ndim,nchi,nb)
      COMPLEX(8),ALLOCATABLE :: mat(:,:,:)     !(nchi,nchi,ndimd)
      logical(4)             :: tchk,treset
      INTEGER(4)             :: NBH
      INTEGER(4)             :: lmnx
      INTEGER(4)             :: npro
      INTEGER(4)             :: IKPT,ISPIN,ibh,ichi,ipro,idim1,idim2,idimd,iat
      INTEGER(4)             :: i1,i2
!     **************************************************************************
      IF(.NOT.TON) RETURN
PRINT*,'ENTERING DMFT$ADDTOHPSI_withkset'
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
      ALLOCATE(DHPIPSI(ndim,NCHI,NB))
      allocate(mat(nchi,nchi,ndimd))
!
!     ==========================================================================
!     ==  convert on-site terms from (txyz) to updown                         ==
!     ==========================================================================
      do iat=1,nat
        lmnx=atomset(iat)%denmat%lmnx
        CALL SPINOR$CONVERT('BACK',lmnx,NDIMD,atomset(iat)%denmat%h)
        atomset(iat)%denmat%h=2.d0*atomset(iat)%denmat%h
      enddo
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
!         == allocate this%htbc if necessary                                  ==
!         ======================================================================
          IF(TRESET) THEN
            IF(.NOT.ASSOCIATED(THIS%HTBC))ALLOCATE(THIS%HTBC(NDIM,NBH,NPRO))
            THIS%HTBC=(0.D0,0.D0)
          END IF
!
!         ======================================================================
!         == DF/DRHO=H0-HRHO; DHDPIPSI= (H0-HRHO)<PI|PSI>                     ==
!         ======================================================================
          mat=kset(ikpt)%H0(:,:,:)-kset(ikpt)%hrho(:,:,:)
          CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,mat)
          mat=2.d0*mat   
          IDIMD=0
          DO IDIM2=1,NDIM
            DO IDIM1=1,NDIM
              IDIMD=IDIM1+NDIM*(IDIM2-1)+ispin-1
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
!         ==  expand to all projector functions                               ==
!         ======================================================================
          DO ICHI=1,NCHI
            IPRO=IPROOFCHI(ICHI)
!           == PIPSI(ICHI,IBH,IKPT,ISPIN)  =THIS%TBC(1,IBH,IPRO) ===============
            THIS%HTBC(:,:,IPRO)=THIS%HTBC(:,:,IPRO)+DHPIPSI(:,ICHI,:nbh)
          ENDDO
!
!         ======================================================================
!         ==  NOW THE SITE-LOCAL TERM FROM THE DOUBLE COUNTING                ==
!         ======================================================================
          DO IAT=1,NAT
            I1=Ichibnd(1,IAT)
            I2=Ichibnd(2,IAT)
            DO IDIM2=1,NDIMD
              DO IDIM1=1,NDIMD
                IDIMD=IDIM1+NDIM*(IDIM2-1)+ISPIN-1
!               == HTBC=HAM*TBC ================================================
                THIS%HTBC(IDIM1,:,I1:I2)=THIS%HTBC(IDIM1,:,I1:I2) &
         &                +MATMUL(THIS%TBC(IDIM1,:,I1:I2) &
         &                       ,TRANSPOSE(ATOMSET(IAT)%denmat%H(:,:,IDIMD)))
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPINOR$TEST()
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NCHI=3
      INTEGER(4),PARAMETER :: NDIMD=4
      COMPLEX(8)           :: A(NCHI,NCHI,NDIMd)
      COMPLEX(8)           :: B(NCHI,NCHI,NDIMd)
      COMPLEX(8)           :: C(NCHI,NCHI,NDIMd)
      COMPLEX(8)           :: d(NCHI,NCHI,NDIMd)
      COMPLEX(8)           :: mat1(2*NCHI,2*NCHI)
      COMPLEX(8)           :: mat2(2*NCHI,2*NCHI)
      INTEGER(4),parameter :: nfil=6
      real(8)              :: ran1,ran2
      INTEGER(4)           :: IDIMD,I,J
      real(8)  ,parameter  :: tol=1.d-7
      real(8)              :: dev
!     **************************************************************************
                                             call trace$push('spinor$test')
print*,' starting spinor test'
      DO IDIMD=1,NDIMD
        DO I=1,NCHI
          DO J=1,NCHI
            CALL RANDOM_NUMBER(RAN1)
            CALL RANDOM_NUMBER(RAN2)
            A(I,J,idimd)=CMPLX(RAN1,RAN2)
          ENDDO
        ENDDO
      ENDDO
      DO IDIMD=1,NDIMD
        DO I=1,NCHI
          DO J=1,NCHI
            CALL RANDOM_NUMBER(RAN1)
            CALL RANDOM_NUMBER(RAN2)
            b(I,J,idimd)=CMPLX(RAN1,RAN2)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  test forward and back transformation                                ==
!     ==========================================================================
      c=a
      CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,A)
      CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,A)
      dev=MAXVAL(ABS(c-A))
      IF(dev.GT.tol) THEN 
        CALL SPINOR$PRINTMATRIX(NFIL,'CONVERSION TEST (ZERO)','DIRECT',1,NCHI &
     &                       ,NDIMD,NCHI,C-a)
        CALL ERROR$MSG('CONVERSION TEST FAILED')
        CALL ERROR$STOP('SPINOR$TEST')
      end if
!
!     ==========================================================================
!     ==  test multiplication                                                 ==
!     ==========================================================================
      C=A
      D=B
      CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,C)
      CALL SPINOR$CONVERT('BACK',NCHI,NDIMD,D)
      MAT1=(0.D0,0.D0)
      mat2=(0.d0,0.d0)
      if(ndimd.eq.1) then
        mat1(:nchi  ,:nchi)  =c(:,:,1)
        mat1(nchi+1:,nchi+1:)=c(:,:,1)
        mat2(:nchi  ,:nchi)  =d(:,:,1)
        mat2(nchi+1:,nchi+1:)=d(:,:,1)
      else if(ndimd.eq.2) then
        mat1(:nchi  ,:nchi)  =c(:,:,1)
        mat1(nchi+1:,nchi+1:)=c(:,:,2)
        mat2(:nchi  ,:nchi)  =d(:,:,1)
        mat2(nchi+1:,nchi+1:)=d(:,:,2)
      else if(ndimd.eq.4) then
        mat1(:nchi  ,:nchi)  =c(:,:,1)
        mat1(nchi+1:,:nchi)  =c(:,:,2)
        mat1(:nchi  ,nchi+1:)=c(:,:,3)
        mat1(nchi+1:,nchi+1:)=c(:,:,4)
        mat2(:nchi  ,:nchi)  =d(:,:,1)
        mat2(nchi+1:,:nchi)  =d(:,:,2)
        mat2(:nchi  ,nchi+1:)=d(:,:,3)
        mat2(nchi+1:,nchi+1:)=d(:,:,4)
      end if
      mat1=matmul(mat1,mat2)
      if(ndimd.eq.1) then
        c(:,:,1)=mat1(:nchi  ,:nchi) 
      else if(ndimd.eq.2) then
        c(:,:,1)=mat1(:nchi  ,:nchi)
        c(:,:,2)=mat1(nchi+1:,nchi+1:)
      else if(ndimd.eq.4) then
        c(:,:,1)=mat1(:nchi  ,:nchi)
        c(:,:,2)=mat1(nchi+1:,:nchi)
        c(:,:,3)=mat1(:nchi  ,nchi+1:)
        c(:,:,4)=mat1(nchi+1:,nchi+1:)
      end if
      CALL SPINOR$CONVERT('FWRD',NCHI,NDIMD,C)
!
      CALL SPINOR$MATMUL(NDIMD,NCHI,A,B,d)
      dev=MAXVAL(ABS(d-c))
      IF(dev.GT.tol) THEN 
        CALL SPINOR$PRINTMATRIX(NFIL,'matmul TEST (ZERO)','DIRECT',1,NCHI &
     &                       ,NDIMD,NCHI,C-d)
        CALL ERROR$MSG('multiplication  TEST FAILED')
        CALL ERROR$STOP('SPINOR$TEST')
      end if
!
!     ==========================================================================
!     ==  test inversion and multiplication                                   ==
!     ==========================================================================
print*,' before spinor$invert'
      CALL SPINOR$INVERT(NDIMD,Nchi,A,B)
!
print*,' before spinor$matmul'
      CALL SPINOR$MATMUL(NDIMD,Nchi,A,B,C)
!
print*,' before spinor$printl'
!      CALL SPINOR$PRINTMATRIX(NFIL,'IDENTITY TEST','TXYZ',1,NCHI,NDIMD,NCHI,C)
      CALL SPINOR$PRINTMATRIX(NFIL,'IDENTITY TEST','UPDOWN',1,NCHI,NDIMD,NCHI,C)
                                             call trace$pop()
      stop
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE spinor$matmul(NDIMD,Nchi,A,B,C)
!     **************************************************************************
!     **  C=A*B                                                               **
!     **                                                                      **
!     **  DERIVATION IN METHODS: SECTION "SECOND QUANTIZATION WITH            **
!     **  NON-ORTHONORMAL BVAIS SETS, SUBSECTION "INVERSION OF A MATRIX IN    **
!     **  SPINOR REPRESENTATION                                               **
!     **************************************PETER BLOECHL GOSLAR 2013***********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: Nchi
      INTEGER(4),INTENT(IN) :: ndimd
      complex(8),intent(in) :: a(nchi,nchi,ndimd)
      complex(8),intent(in) :: b(nchi,nchi,ndimd)
      complex(8),intent(out) :: c(nchi,nchi,ndimd)
      complex(8),parameter   :: ci=(0.d0,1.d0)
      integer(4)             :: idimd,i,ip,ipp
!     **************************************************************************
      c(:,:,:)=0.d0
      c(:,:,1)=matmul(a(:,:,1),b(:,:,1))
      do idimd=2,ndimd
        c(:,:,1)=c(:,:,1)+matmul(a(:,:,idimd),b(:,:,idimd))
        c(:,:,idimd)=c(:,:,idimd)+matmul(a(:,:,1),b(:,:,idimd)) &
     &                           +matmul(a(:,:,idimd),b(:,:,1)) 
      enddo
!
!     == this term is only present in the non-collinear case ===================
      if(ndimd.eq.4) then
        do i=1,3
!         == (x,y,z)->(2,3,4)
          ip=2+mod(i,3)
          ipp=2+mod(i+1,3)
          c(:,:,1+i)=c(:,:,1+i)+ci*(matmul(a(:,:,ip),b(:,:,ipp)) &
     &                             -matmul(a(:,:,ipp),b(:,:,ip)))
        enddo
      end if
!
!     == finally apply factor 1/2 ==============================================
      c(:,:,:)=0.5d0*c(:,:,:)
      return
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
!     **  (uu,du,ud,dd) TO THE (TOTAL,X,Y,Z) REPRESENTATION AND BACK.         **
!     **                                                                      **
!     **  caution! the up-down representation uses the index order of fortran **
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
      INTEGER(4)              :: I,IDIMD,IAT
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
      SUBROUTINE DMFT_PRINTSPINORMATRIX(NFIL,NAME,TYPE,NAT,ICHIBND,NDIMD,NCHI,A)
!     **************************************************************************
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: NAME
      CHARACTER(*),INTENT(IN) :: TYPE
      INTEGER(4)  ,INTENT(IN) :: NDIMD
      INTEGER(4)  ,INTENT(IN) :: NCHI
      INTEGER(4)  ,INTENT(IN) :: NAT
      INTEGER(4)  ,INTENT(IN) :: ICHIBND(2,NAT)
      COMPLEX(8)  ,INTENT(IN) :: A(NCHI,NCHI,NDIMD)
      INTEGER(4)              :: I,IDIMD,IAT
!     **************************************************************************
      DO IAT=1,NAT
        IF(ICHIBND(2,IAT).LT.ICHIBND(1,IAT)) CYCLE
        CALL SPINOR$PRINTMATRIX(NFIL,NAME,TYPE,ICHIBND(1,IAT),ICHIBND(2,IAT) &
     &                         ,NDIMD,NCHI,A)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DMFT_PRINTh0()
!     **************************************************************************
!     ** 
!     **************************************************************************
      use dmft_module
      IMPLICIT NONE
      INTEGER(4)   :: NFIL=6
      INTEGER(4)              :: iat,i1,i2,ispin,i
!     **************************************************************************
      DO IAT=1,NAT
        i1=ichibnd(1,iat)
        i2=ichibnd(2,iat)
        IF(ICHIBND(2,IAT).LT.ICHIBND(1,IAT)) CYCLE
        WRITE(NFIL,FMT='(82("="),T10," ",A," ")')'h0'
        DO Ispin=1,Nspin
          DO I=I1,I2
            WRITE(NFIL,FMT='("Ispin=",I1,":",100("(",2F10.5,")"))') &
    &                    Ispin,h0(I,I1:I2,1,ispin)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
