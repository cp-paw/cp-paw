!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE SESM_MODULE
TYPE ION_type
  REAL(8)            :: AEZ
  INTEGER(4)         :: LMX
  INTEGER(4)         :: LMRX
  INTEGER(4)         :: LPhiX
  INTEGER(4)         :: LRhoX
  INTEGER(4)         :: LREPX
  INTEGER(4)         :: NB
  INTEGER(4),pointer :: lofi(:)
  INTEGER(4),pointer :: nnofi(:)
  INTEGER(4)         :: NBG
  INTEGER(4)         :: NcG
  REAL(8)   ,POINTER :: EBG(:)
  REAL(8)   ,POINTER :: FBG(:)
  REAL(8)   ,POINTER :: aePOT(:,:)   !(NR,LMRX)
  REAL(8)   ,POINTER :: vemb(:,:)  !(NR,LMRX)
  REAL(8)   ,POINTER :: aerho(:,:)   !(NR,LMRX)
  REAL(8)   ,POINTER :: xcpot(:,:)   !(NR,LMRX)
  REAL(8)   ,POINTER :: xceden(:)    !(NR)
  REAL(8)   ,POINTER :: PHI(:,:,:)   !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: TPHI(:,:,:)  !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: HPHI(:,:,:)  !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: testHkinPHI(:,:,:)  !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: testHhartreePHI(:,:,:)  !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: testHxcPHI(:,:,:)  !(NR,LMX,NBG)
  real(8)            :: etot
  real(8)            :: ekin
  real(8)            :: eh
  real(8)            :: exc
END TYPE ION_type
END MODULE SESM_MODULE
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE ONEATOM_MODULE
TYPE THISTYPE
  INTEGER(4)      :: GID
  INTEGER(4)      :: NR
  REAL(8)         :: AEZ=-1.D0
  REAL(8)         :: EREF
  INTEGER(4)      :: NB
  INTEGER(4)      :: NC
  INTEGER(4)      :: LOFI(19)
  INTEGER(4)      :: NNOFI(19)
  REAL(8)         :: FOFI(19)
  REAL(8)         :: EOFI(19)
  REAL(8)         :: EKINC
  REAL(8),POINTER :: atomPOT(:)
  REAL(8),POINTER :: POTh(:)
  REAL(8),POINTER :: POTxc(:)
  REAL(8),POINTER :: AEPOT(:)
  REAL(8),POINTER :: RHOCORE(:)
  REAL(8),POINTER :: DREL(:)
  REAL(8),POINTER :: ETAPAULI(:,:)   ! VALENCE REPULSION POTENTIAL
  REAL(8),POINTER :: ETACPAULI(:,:)  ! CORE REPULSION POTENTIAL
  REAL(8),POINTER :: VCPAULI(:,:)  ! CORE REPULSION POTENTIAL
  REAL(8),POINTER :: OCPAULI(:,:)  ! CORE REPULSION POTENTIAL
  REAL(8),POINTER :: UDOT(:,:)     ! HIGHEST NODLESS PHIDOT FUNCTION
  REAL(8),POINTER :: UN(:,:)       ! HIGHEST NODLESS VALENCE STATE PER L
  REAL(8),POINTER :: RHO(:)     
  REAL(8)         :: Q
  REAL(8)         :: RAD
  REAL(8)         :: DEDRAD
  REAL(8)         :: DEDQ
END TYPE THISTYPE
INTEGER(4)    ,PARAMETER   :: NATX=10
TYPE(THISTYPE),TARGET,SAVE :: THISARR(NATX)
TYPE(THISTYPE),POINTER     :: THIS
LOGICAL(4)    ,PARAMETER   :: TREL=.false.
LOGICAL(4)    ,SAVE        :: TSPECIAL=.FALSE.
END MODULE ONEATOM_MODULE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM MAIN
      CALL DEBUG_PETERS_NONSPH()
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DEBUG_PETERS_NONSPH()
!     **************************************************************************
!     **************************************************************************
      USE ONEATOM_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)                :: GID,GID1
      INTEGER(4)                :: NR
      INTEGER(4)                :: NB        ! #(SHELLS)
      INTEGER(4)                :: LMRX      ! X#(POTENTIAL ANGULAR MOMENTA)
      INTEGER(4) ,ALLOCATABLE   :: LOFI(:)   !(NB)ANGULAR MOMENTUM
      INTEGER(4) ,ALLOCATABLE   :: NN(:)     !(NB)#(NODES)
      REAL(8)    ,ALLOCATABLE   :: POT(:,:)  !(NR,LMRX) POTENTIAL
      REAL(8)    ,ALLOCATABLE   :: E(:)      !(NB) ONE-PARTICLE ENERGIES
      INTEGER(4)                :: LMX       ! #(WAVE FUNCTION ANGULAR MOMENTA)
      INTEGER(4)                :: NBG       
      REAL(8)    ,ALLOCATABLE   :: ebg(:)
      REAL(8)    ,ALLOCATABLE   :: fbg(:)
      REAL(8)    ,ALLOCATABLE   :: phi(:,:,:)
      REAL(8)    ,ALLOCATABLE   :: tphi(:,:,:)
      REAL(8)    ,ALLOCATABLE   :: hphi(:,:)
      REAL(8)    ,ALLOCATABLE   :: DREL(:)
      REAL(8)    ,ALLOCATABLE   :: AUX(:)
      REAL(8)    ,ALLOCATABLE   :: R(:)
      INTEGER(4)                :: LM,IR,IB,ILM,IDIS
      INTEGER(4)                :: NFIL,ISVAR
      REAL(8)                   :: R1,DEX
      REAL(8)                   :: SVAR,DETOT,DEDQ,DEDRAD
      CHARACTER(250)            :: STRING
      CHARACTER(250)            :: STRINGS(7)
      REAL(8)                   :: DIS(7)
      CHARACTER(10)             :: ROOTNAME='ROOT'
      REAL(8)    ,ALLOCATABLE   :: POTREP(:,:)
      REAL(8)                   :: AEZ
      REAL(8)                   :: DR(3),val,dr2(3)
      REAL(8)                   :: cg,pi,y0
      REAL(8)    ,ALLOCATABLE   :: FIN(:,:),FOUT(:,:)
      integer(4)                :: lx1,lm1,lm2,lm3,ibg,l,lx
!     **************************************************************************
      CALL TRACE$SETL4('ON',.FALSE.)
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
!
!     ==========================================================================        
!     ==  connect protocoll file                                              ==        
!     ==========================================================================        
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL FILEHANDLER$SETFILE('PROT',.FALSE.,TRIM(ROOTNAME)//-'.PROT')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','FORM','FORMATTED')
!
!     ==========================================================================        
!     ==  test spherical                                                      ==        
!     ==========================================================================        
!      call PLOTPLGNDR()
!      call testtransform()
!      call comparetransform()
!stop
!
!     ==========================================================================        
!     ==                                                                      ==        
!     ==========================================================================        
      call bigone()
      RETURN
      END SUBROUTINE DEBUG_PETERS_NONSPH


!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE bigone()
!     **************************************************************************
!     **************************************************************************
      use sesm_module
      IMPLICIT NONE
      integer(4),parameter :: nr=250
      integer(4),parameter :: nat=2
      integer(4),parameter :: nbx=19
      integer(4)           :: gid
      type(ion_type)       :: ion(nat)
      real(8)              :: aez(nat)
      real(8)              :: dis
      real(8)              :: pos(3,nat)
      integer(4)           :: nb
      integer(4)           :: iat
!     **************************************************************************
!
!     ==========================================================================
!     ==  initializations                                                     ==
!     ==========================================================================
      CALL RADIAL$NEW('LOG',GID)
      CALL RADIAL$SETI4(GID,'NR',nr)
      CALL RADIAL$SETR8(GID,'DEX',0.05d0)
      CALL RADIAL$SETR8(GID,'R1',1.056d-4)
      CALL DFT$SETI4('TYPE',1)
      CALL DFT$SETL4('SPIN',.FALSE.)
!
!     ==========================================================================
!     ==  system specific                                                     ==
!     ==========================================================================
      ion(:)%aez=7.d0
      ion(:)%lphix=2
      ion(:)%lrhox=2
      ion(:)%lrepx=2*ion(:)%lphix
      ion(:)%nb=3
      nb=3
      do iat=1,nat
        allocate(ion(iat)%lofi(nb)) 
        ion(iat)%lofi(:)=(/0,0,1/)
        allocate(ion(iat)%nnofi(nb))
        ion(iat)%nnofi(:)=(/0,1,0/)
      enddo
      dis=1.1d0/0.529177d0
! dis=1.5d0/0.529177d0
      pos(:,:)=0.d0
      pos(3,2)=dis
!
!     ==========================================================================
!     ==  call sesm                                                           ==
!     ==========================================================================
      call SESM_peter(GID,NR,nat,pos,ion)
      return
      end

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SESM_peter(GID,NR,nat,pos,ion)
!     **************************************************************************
!     **************************************************************************
      USE ONEATOM_MODULE
      use sesm_module
      IMPLICIT NONE
      integer(4),intent(in) :: gid
      integer(4),intent(in) :: nr
      integer(4),intent(in) :: nat
      real(8)   ,intent(in) :: pos(3,nat)
      type(ion_type),intent(inout) :: ion(nat)
      integer(4)            :: nb
      integer(4)            :: nc
      integer(4),allocatable:: lofi(:)
      integer(4),allocatable:: nnofi(:)
      integer(4)            :: lmrx,lmrx1,lmrx2,lmx,lmx1,lmx2
      integer(4)            :: lmax
      real(8)               :: rad(nat)
      real(8)               :: detot,etot,EHARTREE
      real(8)               :: dedq
      real(8)               :: dr(3)
      real(8)               :: dedrad
      real(8)   ,allocatable:: testpotxc(:,:)
      real(8)   ,allocatable:: testpothartree(:,:)
      real(8)   ,allocatable:: pot(:,:)
      real(8)   ,allocatable:: rho(:,:)
      real(8)   ,allocatable:: rhos(:,:)
      real(8)   ,allocatable:: rhoext(:,:)
      real(8)   ,allocatable:: vemb(:,:)
      real(8)   ,allocatable:: FOUT(:,:)
      real(8)   ,allocatable:: FOUT2(:,:)
      real(8)   ,allocatable:: xceden(:,:)
      real(8)   ,allocatable:: aerho(:,:)
      real(8)               :: weight(nr)
      real(8)               :: aux(nr)
      real(8)               :: eden(nr)
      integer(4)            :: nbg
      real(8)   ,allocatable:: ebg(:)
      real(8)   ,allocatable:: fbg(:)
      real(8)   ,allocatable:: phi(:,:,:)
      real(8)   ,allocatable:: tphi(:,:,:)
      real(8)   ,allocatable:: hphi(:,:,:)
      real(8)   ,allocatable:: testhkinphi(:,:,:)
      real(8)   ,allocatable:: testhhartreephi(:,:,:)
      real(8)   ,allocatable:: testhxcphi(:,:,:)
      real(8)               :: drel(nr)
      real(8)               :: ekin,exc,eh,eb
      real(8)               :: sval,hval
      real(8)               :: testhkinval,testhhartreeval,testhxcval
      integer(4)            :: iat,iat1,iat2,ibg,ibg1,ibg2,ichi1,ichi2,ichi
      integer(4)            :: lm1,lm2,lm3,ib
      integer(4)            :: nchi
      character(32)         :: string
      real(8)  ,allocatable :: s(:,:)
      real(8)  ,allocatable :: H(:,:)
      real(8)  ,allocatable :: testHkin(:,:)
      real(8)  ,allocatable :: testHhartree(:,:)
      real(8)  ,allocatable :: testHxc(:,:)
      real(8)  ,allocatable :: eig(:)
      real(8)  ,allocatable :: u(:,:)
      real(8)               :: cg,pi,y0,fourpi
      real(8)               :: svar
      real(8)               :: eband
      real(8)               :: occ
      real(8)               :: r(nr)
      integer(4)            :: lx1
      integer(4)            :: lx,lrx,lmrepx,ir
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      fourpi=4.d0*pi
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
!
!     ==========================================================================
!     ==  CALCULATE ATOM IN A BOX TO OBTAIN THE REPULSIVE POTENTIAL           ==
!     ==========================================================================
      etot=0.d0
      DO IAT1=1,NAT
        CALL ONEATOM$NEW(IAT1,ion(iat1)%AEZ)
        RAD(iat1)=25.D0
        DO IAT2=1,NAT
          IF(IAT2.EQ.IAT1) CYCLE
          RAD(iat1)=MIN(RAD(iat1),sqRT(SUM((POS(:,IAT2)-POS(:,IAT1))**2)))
        ENDDO
        rad(iat1)=0.7d0*rad(iat1)
rad(iat1)=1.483d0
print*,'rad(iat1)',rad(iat1)
        CALL ONEATOM(IAT1,0.D0,RAD(iat1),DETOT,DEDQ,DEDRAD)
        etot=etot+detot
PRINT*,'IAT,ETOT ',IAT1,DETOT
!       == select number of core states
        nc=thisarr(iat1)%nc
        ion(iat1)%ncg=sum(2*thisarr(iat1)%lofi(1:nc)+1)
      ENDDO
      CALL WRITEPHI('ETAPAULI.DAT',GID,NR,2,THISARR(1)%ETAPAULI(:,:))
      CALL WRITEPHI('ETACPAULI.DAT',GID,NR,2,THISARR(1)%ETACPAULI(:,:))
      CALL WRITEPHI('VCPAULI.DAT',GID,NR,2,THISARR(1)%VCPAULI(:,:))
      CALL WRITEPHI('OCPAULI.DAT',GID,NR,2,THISARR(1)%OCPAULI(:,:))
      CALL WRITEPHI('UN.DAT',GID,NR,2,THISARR(1)%UN(:,:))
      CALL WRITEPHI('UDOT.DAT',GID,NR,2,THISARR(1)%UDOT(:,:))
open(991,file='etapaulis.dat')
do iat1=1,nr
  write(991,*)r(iat1),THISARR(1)%ETAPAULI(iat1,1)
enddo
close(991)
open(991,file='setup_peter_.eta_and_vae_hb_offsite_l1.nxy')
do iat1=1,nr
  read(991,*)svar,THISARR(1)%ETAPAULI(iat1,1)
enddo
close(991)
THISARR(2)%ETAPAULI(:,1)=THISARR(1)%ETAPAULI(:,1)

!
!     ==========================================================================
!     == CONSTRUCT EMBEDDING POTENTIAL                                        ==
!     ==========================================================================
print*,'constructing embedding potential..........................'
      DO IAT1=1,NAT
        LMX =(ION(IAT1)%LPHIX+1)**2
        LMRX=(ION(IAT1)%LRHOX+1)**2
        LMREPX=(ION(IAT1)%LREPX+1)**2
        ALLOCATE(VEMB(NR,LMREPX))
        ALLOCATE(RHOEXT(NR,LMREPX))
        ALLOCATE(FOUT(NR,LMREPX))
        ALLOCATE(FOUT2(NR,LMREPX))
        VEMB(:,:)=0.D0
        RHOEXT(:,:)=0.D0
        DO IAT2=1,NAT
          IF(IAT2.EQ.IAT1) CYCLE
          DR(:)=POS(:,IAT1)-POS(:,IAT2)
!         == ADD UP EFFECTIVE POTENTIALS OF IONS IN A BOX ========================
          AUX(:)=THISARR(IAT2)%ETAPAULI(:,1)+THISARR(IAT2)%POTH(:)+THISARR(IAT2)%POTXC(:)
          CALL SPHERICAL$SHIFTCENTER(GID,NR,DR,1,AUX,LMREPX,FOUT)
          VEMB(:,:)=VEMB(:,:)+FOUT(:,:)
!         == ADD UP DENSITIES ====================================================
          AUX(:)=THISARR(IAT2)%RHO(:)
          CALL SPHERICAL$SHIFTCENTER(GID,NR,DR,1,AUX,LMREPX,FOUT)
          RHOEXT(:,:)=RHOEXT(:,:)+FOUT(:,:)
!         == SUBTRACT XC-POTENTIAL OF EXTERNAL DENSITIES FROM VEMB================
          CALL AUGMENTATION_XC(GID,NR,LMREPX,1,FOUT,SVAR,FOUT2,AUX)
          DO IR=1,NR
            IF(FOUT(IR,1)*Y0.LT.1.D-6) FOUT2(IR,:)=0.D0
          ENDDO
          VEMB(:,:)=VEMB(:,:)-FOUT2(:,:)
        ENDDO
!       ==  SUBTRACT ONSITE XC-POTENTIAL ========================================
        VEMB(:,1)=VEMB(:,1)-THISARR(IAT1)%POTXC(:)
!       == CONSTRUCT XC-POTENTIAL OF SUPERIMPOSED DENSITIES
        FOUT(:,:)=RHOEXT(:,:)
        FOUT(:,1)=FOUT(:,1)+THISARR(IAT1)%RHO(:)
        CALL AUGMENTATION_XC(GID,NR,LMREPX,1,FOUT,SVAR,FOUT2,AUX)
        VEMB(:,:)=VEMB(:,:)+FOUT2(:,:)
!
!       == MAP EMBEDDING POTENTIAL ON ION STRUCTURE =============================
        ALLOCATE(ION(IAT1)%VEMB(NR,LMREPX))
        ION(IAT1)%VEMB(:,:)=VEMB(:,:)
!
!       == WRAP UP ==============================================================
        DEALLOCATE(FOUT)
        DEALLOCATE(FOUT2)
        DEALLOCATE(VEMB)
        DEALLOCATE(RHOEXT)
      ENDDO
      CALL WRITEPHI('VEMBEDDING.DAT',GID,NR,LMREPX,ION(1)%VEMB)
!
!     ============================================================================
!     == self-consistent calculation of the deformed ion                        ==
!     ============================================================================
print*,'scf of deformed ion...........................................'
      ETOT=0.D0
      DO IAT1=1,NAT
        LMX =(ION(IAT1)%LPHIX+1)**2
        LMRX=(ION(IAT1)%LRHOX+1)**2
        LMREPX=(ION(IAT1)%LREPX+1)**2
        LX=SQRT(REAL(LMX-1,KIND=8)+1.D-10)
        nc=thisarr(iat1)%nc
        lmax=maxval(thisarr(iat1)%lofi(1:nc))
!
!       == collect data ==========================================================
        NB=ION(IAT1)%NB
        ALLOCATE(LOFI(NB))
        ALLOCATE(NNOFI(NB))
        LOFI(:)=ION(IAT1)%LOFI(:)
        NNOFI(:)=ION(IAT1)%NNOFI(:)
!
!       ==  allocate arrays ======================================================
        ALLOCATE(AERHO(NR,LMRX))
        NBG=SUM(2*LOFI(:)+1)
        ALLOCATE(EBG(NBG))
        ALLOCATE(FBG(NBG))
        ALLOCATE(PHI(NR,LMX,NBG))
        ALLOCATE(TPHI(NR,LMX,NBG))
!
!       == START POTENTIAL =======================================================
        ALLOCATE(POT(NR,LMRX))
        POT(:,:)=0.D0
        POT(:,1)=THISARR(IAT1)%ATOMPOT(:)
!
!       ==========================================================================
!       == perform scf calculation of deformed ion                              ==
!       ==========================================================================
        DREL(:)=0.D0
        CALL SCF(GID,NR,LMRX,LMX,lmrepx,ion(iat1)%AEZ,DREL,ion(iat1)%vemb &
       &        ,lmax,ion(iat1)%vcpauli(1:lmax),ion(iat1)%ocpauli(1:lmax) &
       &        ,NB,LOFI,NNofi,POT,aerho,NBG,EBG,FBG,PHI,TPHI)
!
!       ==========================================================================
!       == calculate total energy                                               ==
!       ==========================================================================
        call deform$etot(gid,nr,lmrx,lmx,lmrepx,nbg,ion(iat1)%aez &
     &                  ,fbg,ebg,phi,tphi,detot,eb,ekin,eh,exc)
        etot=etot+detot-thisarr(iat1)%eref
!
!       == report energies =======================================================
        write(*,fmt='(72("="))')
        write(*,fmt='(72("="),t20," energies of the deformed ion ")')
        write(*,fmt='("eion    =",f10.5)')detot
        write(*,fmt='("edeform =",f10.5)')detot-thisarr(iat1)%eref
        write(*,fmt='("ekin    =",f10.5)')ekin
        write(*,fmt='("ehartree=",f10.5)')eh
        write(*,fmt='("exc     =",f10.5)')exc
        do ibg=1,nbg
          write(*,*)ibg,fbg(ibg),ebg(ibg),ebg(ibg)/27.211d0
        enddo
        write(*,fmt='(72("="))')
!
!       == map result onto ion structure =========================================
        CALL WRITEPHI('POTIN.DAT',GID,NR,LMRX,POT)
        allocate(ion(iat1)%ebg(nbg))
        allocate(ion(iat1)%fbg(nbg))
        allocate(ion(iat1)%phi(nr,lmx,nbg))
        allocate(ion(iat1)%tphi(nr,lmx,nbg))
        allocate(ion(iat1)%aepot(nr,lmrx))
        allocate(ion(iat1)%aerho(nr,lmrx))
        ion(iat1)%nbg=nbg
        ion(iat1)%ebg(:)=ebg(:)
        ion(iat1)%fbg(:)=fbg(:)
        ion(iat1)%phi(:,:,:)=phi(:,:,:)
        ion(iat1)%tphi(:,:,:)=tphi(:,:,:)
        ion(iat1)%aepot(:,:)=pot(:,:)
        ion(iat1)%aerho(:,:)=aerho(:,:)
!       == wrap up ===============================================================
        DEALLOCATE(EBG)
        DEALLOCATE(FBG)
        DEALLOCATE(PHI)
        DEALLOCATE(TPHI)
        DEALLOCATE(lofi)
        DEALLOCATE(nnofi)
        deallocate(pot)
        deallocate(aerho)
      ENDDO
!
!     ============================================================================
!     == inspect results                                                        ==
!     ============================================================================
      iat1=1
      lmrx=(ion(iat1)%lrhox+1)**2
      do ibg=1,nbg
        write(string,*)ibg
        string='phi'//trim(adjustl(string))//'.dat'
        CALL WRITEPHI(string,GID,NR,LMX,ION(IAT1)%phi(:,:,ibg))
      enddo
      CALL WRITEPHI('AEPOT.DAT',GID,NR,LMRX,ion(iat1)%aePOT)
      CALL WRITEPHI('AERHO.DAT',GID,NR,LMRX,ion(iat1)%aeRHO)
!
!     ==========================================================================
!     == determine hartree energy                                             ==
!     ==========================================================================
print*,'calculating Hartree energy ...'
      ehartree=0.d0
      do iat1=1,nat
         lmrx1=(ion(iat1)%lrhox+1)**2
        do iat2=iat1+1,nat
          lmrx2=(ion(iat2)%lrhox+1)**2
          call hartreepair(GID,NR,pos(:,iat2)-pos(:,iat1) &
     &                    ,lmrx1,ion(iat1)%aez,ion(iat1)%aerho &
     &                    ,lmrx2,ion(iat2)%aez,ion(iat2)%aerho,detot)
          ehartree=ehartree+detot
        enddo
      enddo
      print*,'ehartree ',ehartree
      etot=etot+ehartree
      print*,'etot after ehartree ',etot
!
!     ==========================================================================
!     == CALCULATE xc-correction                                              ==
!     ==========================================================================
      do iat1=1,nat
        LMRX=(ION(IAT1)%Lrhox+1)**2
        allocate(ion(iat1)%xcpot(nr,lmrx))
        allocate(ion(iat1)%xceden(nr))
        CALL AUGMENTATION_XC(GID,NR,LMRX,1,ion(iat1)%aerho,svar &
     &                      ,ion(iat1)%xcpot,ion(iat1)%xceden)
!       == cut out low density region ==========================================
        do ir=1,nr
          if(ion(iat1)%aerho(ir,1)*y0.lt.1.d-6) ion(iat1)%xcpot(ir,:)=0.d0
        enddo
       enddo
!
!     ==========================================================================
!     == CALCULATE H|PHI>
!     ==========================================================================
print*,'calculating h|\phi>....'
      exc=0.d0
      DO IAT1=1,NAT
        NBG=ION(IAT1)%NBG
        LMX =(ION(IAT1)%Lphix+1)**2
        LMRX=(ION(IAT1)%Lrhox+1)**2
        ALLOCATE(PHI(NR,LMX,NBG))
        ALLOCATE(HPHI(NR,LMX,NBG))
        ALLOCATE(testHkinPHI(NR,LMX,NBG))
        ALLOCATE(testHhartreePHI(NR,LMX,NBG))
        ALLOCATE(testHxcPHI(NR,LMX,NBG))
        PHI(:,:,:)=ION(IAT1)%PHI(:,:,:)
!
!       ========================================================================
!       == construct sum of potentials and sum of densities                   ==
!       ========================================================================
        ALLOCATE(POT(NR,LMRX))
        ALLOCATE(testPOThartree(NR,LMRX))
        ALLOCATE(testPOTxc(NR,LMRX))
        ALLOCATE(rho(NR,LMRX))
        ALLOCATE(rhos(NR,LMRX))
        ALLOCATE(xceden(NR,lmrx))
        ALLOCATE(fout(NR,LMRX))
        pot(:,:)=0.d0
        rho(:,:)=0.d0
        rhos(:,:)=0.d0
        xceden(:,:)=0.d0
        do iat2=1,nat
          if(iat2.eq.iat1) cycle
          LMRX2=(ION(IAT2)%Lrhox+1)**2
!         == total potential except xc potential ===============================
          CALL SPHERICAL$shiftcenter(GID,NR,pos(:,iat1)-pos(:,iat2) &
    &                  ,LMrx2,ion(iat2)%aepot-ion(iat2)%xcpot,Lmrx,fout(:,:))
          pot(:,:)=pot(:,:)+fout(:,:)
!
!         ==  total density ====================================================
          CALL SPHERICAL$shiftcenter(GID,NR,pos(:,iat1)-pos(:,iat2) &
    &                           ,lmrx2,ion(iat2)%aerho,LMrX,fout(:,:))
          rho(:,:)=rho(:,:)+fout(:,:)
!
!         == weighting function for xc energy density ==========================
          aux(:)=exp(-r(:)**2)
          CALL SPHERICAL$shiftcenter(GID,NR,pos(:,iat1)-pos(:,iat2) &
    &                           ,1,aux,LMrX,fout(:,:))
          rhos(:,:)=rhos(:,:)+fout(:,:)
!
!         == xc energy density                       ========================
          CALL SPHERICAL$shiftcenter(GID,NR,pos(:,iat1)-pos(:,iat2) &
    &                           ,1,ion(iat2)%xceden,lmrx,fout)
          xceden(:,:)=xceden(:,:)-fout(:,:)
        enddo
        pot(:,:)=ION(IAT1)%AEPOT(:,:)-ION(IAT1)%xcPOT(:,:)+pot(:,:)
        testpothartree(:,:)=pot(:,:)
        rho(:,:)=ION(IAT1)%AErho(:,:)+rho(:,:)
        CALL AUGMENTATION_XC(GID,NR,LMRX,1,RHO,svar,fout,aux)
!       == cut out low density region ==========================================
        do ir=1,nr
          if(rho(ir,1)*y0.lt.1.d-6) fout(ir,:)=0.d0
        enddo
 CALL WRITEPHI('xcpot.DAT',GID,NR,lmrx,fout)
        pot=pot+fout
        testpotxc(:,:)=fout(:,:)
        xceden(:,1)=xceden(:,1)-ION(IAT1)%xceden(:)+aux
 CALL WRITEPHI('xceden.DAT',GID,NR,lmrx,xceden)
 CALL WRITEPHI('xceden1.DAT',GID,NR,1,ion(1)%xceden)
 CALL WRITEPHI('xceden2.DAT',GID,NR,1,ion(2)%xceden)
!
!       ========================================================================
!       == correct xc energy                                                  ==
!       ========================================================================
!       == weighting
        aux(:)=exp(-r(:)**2)
        rhos(:,1)=aux(:)+rhos(:,1)
        aux(:)=aux(:)/(rhos(:,1)*y0)
        do lm1=2,lmrx
          rhos(:,lm1)=rhos(:,lm1)/(rhos(:,1)*y0)
          rhos(:,lm1)=rhos(:,lm1)*aux(:)
        enddo
        rhos(:,1)=aux(:)
        aux(:)=0.d0
        do lm1=1,lmrx
          aux(:)=aux(:)+rhos(:,lm1)*xceden(:,lm1)
        enddo
 CALL WRITEPHI('wxceden.DAT',GID,NR,1,aux)
 CALL WRITEPHI('wght.DAT',GID,NR,lmrx,rhos)
        call radial$integral(gid,nr,aux*r(:)**2,svar)
        exc=exc+fourpi*svar*y0
!
print*,' exc?? ',fourpi*svar*y0
if(iat1.eq.1) then
  CALL WRITEPHI('vtotc1.DAT',GID,NR,LMRX,pot)
  CALL WRITEPHI('rhototc1.DAT',GID,NR,LMRX,rho)
  CALL WRITEPHI('rhosc1.DAT',GID,NR,LMRX,rhos)
  CALL WRITEPHI('rhos1.DAT',GID,NR,1,ion(iat1)%aerho(:,1))
else if(iat1.eq.2) then
  CALL WRITEPHI('vtotc2.DAT',GID,NR,LMRX,pot)
  CALL WRITEPHI('rhototc2.DAT',GID,NR,LMRX,rho)
  CALL WRITEPHI('rhosc2.DAT',GID,NR,LMRX,rhos)
end if
        HPHI(:,:,:)=0.D0
        testHhartreePHI(:,:,:)=0.D0
        testHxcPHI(:,:,:)=0.D0
        DO IBG=1,NBG
          DO LM1=1,LMX
            DO LM2=1,LMRX
              DO LM3=1,LMX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                HPHI(:,LM1,IBG)=HPHI(:,LM1,IBG)+CG*POT(:,LM2)*PHI(:,LM3,IBG)
                testHhartreePHI(:,LM1,IBG)=testHhartreePHI(:,LM1,IBG)+CG*testPOThartree(:,LM2)*PHI(:,LM3,IBG)
                testHxcPHI(:,LM1,IBG)=testHxcPHI(:,LM1,IBG)+CG*testPOTxc(:,LM2)*PHI(:,LM3,IBG)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        deallocate(phi)
        deallocate(pot)
        deallocate(testpothartree)
        deallocate(testpotxc)
        deallocate(rho)
        deallocate(rhos)
        deallocate(xceden)
        deallocate(fout)
        ALLOCATE(ION(IAT1)%HPHI(NR,LMX,NBG))
        ALLOCATE(ION(IAT1)%testHkinPHI(NR,LMX,NBG))
        ALLOCATE(ION(IAT1)%testHhartreePHI(NR,LMX,NBG))
        ALLOCATE(ION(IAT1)%testHxcPHI(NR,LMX,NBG))
        ION(IAT1)%HPHI(:,:,:)=ION(IAT1)%TPHI(:,:,:)+HPHI(:,:,:)
        ION(IAT1)%testHkinPHI(:,:,:)=ION(IAT1)%TPHI(:,:,:)
        ION(IAT1)%testHhartreePHI(:,:,:)=testhhartreephi(:,:,:)
        ION(IAT1)%testHxcPHI(:,:,:)=testhxcphi(:,:,:)
        deallocate(hphi)
        deallocate(testhhartreephi)
        deallocate(testhxcphi)
        deallocate(testhkinphi)
      ENDDO
print*,'delta exc ',exc
!      etot=etot+exc
      print*,'etot after exchange-correlation energy ',etot
!
!     ==========================================================================
!     == CALCULATE HAMILTON AND OVERLAP MATRIX                                ==
!     ==========================================================================
print*,'calculating Hamilton and overlapmatrix ...'
      NCHI=0
      DO IAT1=1,NAT
        NCHI=NCHI+ION(IAT1)%NBG
      ENDDO
      allocate(s(nchi,nchi))
      allocate(H(nchi,nchi))
      allocate(testHkin(nchi,nchi))
      allocate(testHxc(nchi,nchi))
      allocate(testHhartree(nchi,nchi))
      S(:,:)=0.d0
      H(:,:)=0.d0
      ICHI1=0
      DO IAT1=1,NAT
        lmx1=(ion(iat1)%lphix+1)**2
        do ibg1=1,ion(iat1)%nbg
          H(ichi1+ibg1,ichi1+ibg1)=ion(iat1)%ebg(ibg1)
          s(ichi1+ibg1,ichi1+ibg1)=1.d0
        enddo
        ICHI2=0
        DO IAT2=1,NAT
!          IF(Iat1.Ge.Iat2) then
          lmx2=(ion(iat2)%lphix+1)**2
          IF(Iat1.Ge.0) then
            DO IBG1=1,ION(IAT1)%NBG
              DO IBG2=1,ION(IAT2)%NBG 
                CALL OVERLAP(GID,NR,POS(:,IAT2)-POS(:,IAT1) &
     &            ,LMX1,ION(IAT1)%PHI(:,:,IBG1) &
     &            ,LMX2,ion(iat2)%hphi(:,:,ibg2) &
     &                 ,ION(IAT2)%testhkinPHI(:,:,IBG2) &
     &                 ,ION(IAT2)%testhhartreePHI(:,:,IBG2) &
     &                 ,ION(IAT2)%testhxcPHI(:,:,IBG2) &
     &                 ,ION(IAT2)%PHI(:,:,IBG2) &
     &                 ,hval,testhkinval,testhhartreeval,testhxcval,sval)
                h(ichi1+ibg1,ichi2+ibg2)=hval
                testhkin(ichi1+ibg1,ichi2+ibg2)=testhkinval
                testhhartree(ichi1+ibg1,ichi2+ibg2)=testhhartreeval
                testhxc(ichi1+ibg1,ichi2+ibg2)=testhxcval
!                h(ichi2+ibg2,ichi1+ibg1)=hval
                s(ichi1+ibg1,ichi2+ibg2)=sval
!                s(ichi2+ibg2,ichi1+ibg1)=sval
              ENDDO
            ENDDO
          END IF
          ICHI2=ICHI2+ION(IAT2)%NBG
        ENDDO
        ICHI1=ICHI1+ION(IAT1)%NBG
      ENDDO
      WRITE(*,FMT='(72("="),T20," HAMILTON MATRIX ")')
      DO ICHI1=1,NCHI
        WRITE(*,FMT='(20F10.5)')H(ICHI1,:)
      ENDDO
      WRITE(*,FMT='(72("="),T20," KINETIC ENERGY  MATRIX ")')
      DO ICHI1=1,NCHI
        WRITE(*,FMT='(20F10.5)')TESTHKIN(ICHI1,:)
      ENDDO
      WRITE(*,FMT='(72("="),T20," HARTREE ENERGY  MATRIX ")')
      DO ICHI1=1,NCHI
        WRITE(*,FMT='(20F10.5)')TESTHHARTREE(ICHI1,:)
      ENDDO
      WRITE(*,FMT='(72("="),T20," XC ENERGY  MATRIX ")')
      DO ICHI1=1,NCHI
        WRITE(*,FMT='(20F10.5)')TESTHXC(ICHI1,:)
      ENDDO
      WRITE(*,FMT='(72("="),T20," OVERLAP MATRIX ")')
      DO ICHI1=1,NCHI
        WRITE(*,FMT='(20F10.5)')S(ICHI1,:)
      enddo
      deallocate(testHkin)
      deallocate(testHxc)
      deallocate(testHhartree)
!
!     ==========================================================================
!     == solve schroedinegr equation                                          ==
!     ==========================================================================
      s=0.5d0*(s+transpose(s))
      h=0.5d0*(h+transpose(h))
      allocate(u(nchi,nchi))
      allocate(eig(nchi))
      call lib$generaleigenvaluer8(nchi,h,s,eig,u)
!
!     ==========================================================================
!     == calculate difference in bandstucture energy                          ==
!     ==========================================================================
      eband=0.d0
      svar=0.d0
      ichi=0
      do iat=1,nat
        svar=svar+ion(iat)%aez
        nbg=ion(iat)%nbg
        do ib=1,nbg
          ichi=ichi+1
          eband=eband-ion(iat)%fbg(ib)*H(ichi,ichi)
        enddo 
      enddo
      write(*,fmt='("bandstructure energy(isolated) ",f10.5)')-eband
      do ib=1,nchi
        occ=min(2.d0,svar)
        write(*,fmt='(i3,"occ=",f5.2," e[h]=",f10.5," e[ev]=",f10.5)')ib,occ,eig(ib),eig(ib)*27.211d0
        svar=svar-occ
        eband=eband+eig(ib)*occ
      enddo
      deallocate(u)
      write(*,fmt='("bandstructure energy(difference) ",f10.5)')eband
      etot=etot+eband
      print*,'etot after bandstructure energy ',etot
!
!     ==========================================================================
!     == wrap up                                                              ==
!     ==========================================================================
      deallocate(s)
      deallocate(H)
      deallocate(eig)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE overlap(GID,NR,dr21,LMX1,phi1,lmx2,hphi2 &
     &                  ,testhkinphi2,testhhartreephi2,testhxcphi2 &
     &                  ,phi2,h,testhkin,testhhartree,testhxc,s)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: GID
      INTEGER(4) ,INTENT(IN)    :: NR
      REAL(8)    ,INTENT(IN)    :: dr21(3)
      INTEGER(4) ,INTENT(IN)    :: LMX1       ! X#(POTENTIAL ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)    :: LMX2       ! X#(POTENTIAL ANGULAR MOMENTA)
      REAL(8)    ,INTENT(IN)    :: phi1(NR,lmx1)
      REAL(8)    ,INTENT(IN)    :: hphi2(NR,lmx2)
      REAL(8)    ,INTENT(IN)    :: testhkinphi2(NR,lmx2)
      REAL(8)    ,INTENT(IN)    :: testhhartreephi2(NR,lmx2)
      REAL(8)    ,INTENT(IN)    :: testhxcphi2(NR,lmx2)
      REAL(8)    ,INTENT(IN)    :: phi2(NR,lmx2)
      REAL(8)    ,INTENT(OUT)   :: h
      REAL(8)    ,INTENT(OUT)   :: testhkin,testhhartree,testhxc
      REAL(8)    ,INTENT(OUT)   :: s
      real(8)                   :: shiftedphi(nr,lmx2)
      real(8)                   :: aux(nr)
      integer(4)                :: lm      
      real(8)                   :: r(nr)
      real(8)                   :: dis
!     **************************************************************************
      call radial$r(gid,nr,r)
!
!     ==========================================================================
!     ==  shift density rho2 to position 1                                    ==
!     ==========================================================================
      dis=sqrt(sum(dr21**2))
      if(dis.gt.1.d-5) then
        CALL SPHERICAL$shiftcenter(GID,NR,DR21,LMx1,PHI1,LMX2,SHIFTEDphi)
      else
        lm=min(lmx1,lmx2)
        shiftedphi=0.d0
        shiftedphi(:,:lm)=phi1(:,:lm)
      endif
!
!     ==========================================================================
!     ==  CALCULATE energy                                                    ==
!     ==========================================================================
      AUX(:)=0.D0
      DO LM=1,LMX2
        AUX(:)=AUX(:)+shiftedphi(:,LM)*Hphi2(:,LM)
      ENDDO
      aux(:)=aux(:)*r(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,h)
!
!     ==test=========================================================
      AUX(:)=0.D0
      DO LM=1,LMX2
        AUX(:)=AUX(:)+shiftedphi(:,LM)*testHkinphi2(:,LM)
      ENDDO
      aux(:)=aux(:)*r(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,testhkin)
!
!     ==test=========================================================
      AUX(:)=0.D0
      DO LM=1,LMX2
        AUX(:)=AUX(:)+shiftedphi(:,LM)*testHhartreephi2(:,LM)
      ENDDO
      aux(:)=aux(:)*r(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,testhhartree)
!
!     ==test=========================================================
      AUX(:)=0.D0
      DO LM=1,LMX2
        AUX(:)=AUX(:)+shiftedphi(:,LM)*testHxcphi2(:,LM)
      ENDDO
      aux(:)=aux(:)*r(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,testhxc)
!
!     == overlap =========================================================
      AUX(:)=0.D0
      DO LM=1,LMX2
        AUX(:)=AUX(:)+shiftedphi(:,LM)*phi2(:,LM)
      ENDDO
      aux(:)=aux(:)*r(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,s)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE hartreepair(GID,NR,dr21,LMRX1,aez1,rho1,lmrx2,aez2,rho2,e)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: GID
      INTEGER(4) ,INTENT(IN)    :: NR
      REAL(8)    ,INTENT(IN)    :: dr21(3)
      INTEGER(4) ,INTENT(IN)    :: LMRX1       ! X#(POTENTIAL ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)    :: LMRX2       ! X#(POTENTIAL ANGULAR MOMENTA)
      REAL(8)    ,INTENT(IN)    :: AEZ1
      REAL(8)    ,INTENT(IN)    :: AEZ2
      REAL(8)    ,INTENT(IN)    :: RHO1(NR,lmrx1)
      REAL(8)    ,INTENT(IN)    :: RHO2(NR,lmrx2)
      REAL(8)    ,INTENT(OUT)   :: E
      real(8)                   :: pot(nr,lmrx1)
      real(8)                   :: shiftedpot(nr,lmrx2)
      real(8)                   :: aux(nr)
      integer(4)                :: l,lm      
      integer(4)                :: lx1      
      real(8)                   :: pi,y0
      real(8)                   :: dis
      real(8)                   :: r(nr)
      real(8)                   :: val
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      call radial$r(gid,nr,r)
!
!     ==========================================================================
!     ==  calculate hartree potential                                         ==
!     ==========================================================================
      do lm=1,lmrx1
        l=sqrt(real(lm-1,kind=8)+1.d-10)
        call radial$poisson(gid,nr,l,rho1(:,lm),pot(:,lm))
      enddo
!     == add nuclear charge ====================================================
      pot(:,1)=pot(:,1)-aez1/r(:)/y0
!
!     ==========================================================================
!     ==  shift density rho2 to position 1                                    ==
!     ==========================================================================
      CALL SPHERICAL$shiftcenter(GID,NR,DR21,LMRX1,pot,LMRX2,SHIFTEDpot)
!
!     ==========================================================================
!     ==  CALCULATE energy                                                    ==
!     ==========================================================================
      AUX(:)=0.D0
      DO LM=1,LMRX2
        AUX(:)=AUX(:)+shiftedPOT(:,LM)*rho2(:,LM)
      ENDDO
      aux(:)=aux(:)*r(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,E)
      call radial$value(gid,nr,shiftedpot(:,1),0.d0,val)
      e=e-aez2*val*y0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCF(GID,NR,LMRX,LMX,lmrepx,AEZ,DREL,vemb,lmax,vcpauli,ocpauli &
     &              ,NB,LOFI,NN,POT,rho,nbg,ebg,fbg,phi,tphi)
!     **************************************************************************
!     **  the deformed ion fulfilles the schroedinger equation with           **
!     **  potential "pot+potrep", where                                       **
!     **  pot contains                                                        **
!     **     1) the hartree potential of this atom                            **
!     **     2) the xc potential of this and the neighboring atoms            **
!     **  potrep contains:                                                    **
!     **     1) the hartree potential of the neighboring atoms                **
!     **     2) the Pauli repulsion of the neighboring atoms                  **
!     **  the embedding potential is obtained as                              **
!     **     vemb= potrep
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: GID
      INTEGER(4) ,INTENT(IN)    :: NR
      INTEGER(4) ,INTENT(IN)    :: LMRX       ! X#(POTENTIAL ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)    :: LMX        ! X#(WAVEF. ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)    :: LMrepX     ! X#(WAVEF. ANGULAR MOMENTA)
      REAL(8)    ,INTENT(IN)    :: AEZ
      REAL(8)    ,INTENT(IN)    :: DREL(NR)
      REAL(8)    ,INTENT(IN)    :: vemb(NR,LMrepX)
      INTEGER(4) ,INTENT(IN)    :: LMax
      REAL(8)    ,INTENT(IN)    :: vcpauli(NR,LMax)
      REAL(8)    ,INTENT(IN)    :: ocpauli(NR,LMax)
      REAL(8)    ,INTENT(INOUT) :: POT(NR,LMRX)  ! POTENTIAL
      REAL(8)    ,intent(out)   :: rho(nr,lmrx)    ! all-electron density
      INTEGER(4) ,INTENT(IN)    :: NB            ! #(SHELLS)
      INTEGER(4) ,INTENT(IN)    :: LOFI(NB)      ! ANGULAR MOMENTA
      INTEGER(4) ,INTENT(IN)    :: NN(NB)        ! #(NODES)
      INTEGER(4) ,intent(in)    :: NBG       
      REAL(8)    ,intent(out)   :: EBG(nbg)        ! energy eigenvalues 
      REAL(8)    ,intent(out)   :: fBG(nbg)        ! occupations
      REAL(8)    ,intent(out)   :: PHI(nr,lmx,nbg) ! wave functions
      REAL(8)    ,intent(out)   :: TPHI(nr,lmx,nbg)! Kinetic energy times phi
      REAL(8)    ,allocatable   :: POTin(:,:)
      REAL(8)    ,allocatable   :: rhotot(:,:)
      REAL(8)                   :: POTOUT(NR,LMRX)
      REAL(8)                   :: POTbig(NR,LMRepX)
      REAL(8)                   :: r(NR)
      REAL(8)                   :: AUX(NR),aux1(nr)
      LOGICAL(4)                :: TBROYDEN=.FALSE.
      INTEGER(4) ,PARAMETER     :: NITER=50
      INTEGER(4)                :: ITER,LM1,LM2,LM3,IBG,IB,L,LM,ir
      INTEGER(4)                :: LX
      INTEGER(4)                :: Lmpotinx
      CHARACTER(64)             :: STRING
      REAL(8)                   :: MAXDEV
      REAL(8)    ,PARAMETER     :: TOL=1.D-6
      REAL(8)                   :: etot,eb,eh,EXC,epot
      REAL(8)    ,ALLOCATABLE   :: FOFBG(:)
      LOGICAL(4)                :: TCONV
      REAL(8)                   :: SVAR
      REAL(8)                   :: CG
      REAL(8)                   :: pi,y0
!     **************************************************************************
      TBROYDEN=.FALSE.
      call radial$r(gid,nr,r)
      LX=INT(SQRT(REAL(LMX-1)+1.D-8))
      if(NBG.ne.SUM(2*LOFI(:)+1)) then
        call error$stop('scf')
      end if
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
!  
!     ==========================================================================        
!     == allocate arrays                                                      ==
!     ==========================================================================        
      lmpotinx=max(lmrx,lmrepx)
      allocate(potin(nr,lmpotinx))
      allocate(rhotot(nr,lmpotinx))
!  
!     ==========================================================================        
!     == determine occupations                                                ==
!     ==========================================================================        
      SVAR=AEZ
      FBG(:)=0.D0
      IBG=0
      DO IB=1,NB
        L=LOFI(IB)
        IF(SVAR.GT.REAL(2*(2*L+1))) THEN        
          FBG(IBG+1:IBG+2*L+1)=2.D0
          SVAR=SVAR-2.D0*REAL(2*L+1,KIND=8)
          IBG=IBG+2*L+1
        ELSE 
          FBG(IBG+1:IBG+2*L+1)=SVAR/REAL(2*L+1,KIND=8)
          EXIT
        END IF
      ENDDO
!  
!     ==========================================================================        
!     == start self-consistency loop                                          ==
!     ==========================================================================        
      DO ITER=1,NITER
!  
!       ========================================================================        
!       == NONSCF CALCULATION                                                 ==
!       ========================================================================        
        potin(:,:)=vemb(:,:)
        potin(:,:lmrx)=potin(:,:lmrx)+pot(:,:)
        CALL NONSCF(GID,NR,LMPOTINX,POTIN,DREL,NB,LOFI,NN,NBG,EBG,LX,PHI,TPHI)
!  
!       ========================================================================        
!       == CALCULATE DENSITY                                                  ==
!       ========================================================================        
        RHO(:,:)=0.D0
        DO IBG=1,NBG
          DO LM1=1,LMX
            DO LM2=1,LMX
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                RHO(:,LM3)=RHO(:,LM3)+FBG(IBG)*CG*PHI(:,LM1,IBG)*PHI(:,LM2,IBG)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!  
!       ========================================================================        
!       == CALCULATE POTENTIAL                                                ==
!       ========================================================================
        POTOUT(:,:)=0.D0
        CALL AUGMENTATION_XC(GID,NR,LMRX,1,RHO,EXC,POTOUT,AUX)
!       == CUT OUT LOW DENSITY REGION ==========================================
        DO IR=1,NR
          IF(RHO(IR,1)*Y0.LT.1.D-6) POTOUT(IR,:)=0.D0
        ENDDO
!       
        CALL RADIAL$NUCPOT(GID,NR,AEZ,AUX)
        POTOUT(:,1)=POTOUT(:,1)+AUX(:)
        DO LM=1,LMRX
          L=INT(SQRT(REAL(LM-1)+1.D-10))
          CALL RADIAL$POISSON(GID,NR,L,RHO(:,LM),AUX)
          POTOUT(:,LM)=POTOUT(:,LM)+AUX(:)
        ENDDO
!  
!       ========================================================================        
!       == PRINT WAVE FUNCTIONS                                               ==
!       ========================================================================        
        WRITE(*,FMT='("iter",i5," EBG: ",10F10.5)')iter,EBG
!  
!       ========================================================================        
!       == MIX POTENTIALS                                                     ==
!       ========================================================================        
        MAXDEV=MAXVAL(ABS(POTOUT-POT))
        TCONV=MAXDEV.LT.TOL
        IF(TCONV) EXIT        
        IF(ITER.EQ.2) THEN
          TBROYDEN=.TRUE.
          CALL BROYDEN$NEW(NR*LMRX,10,1.D-1)
        END IF
        IF(TBROYDEN) THEN
          CALL BROYDEN$STEP(NR*LMRX,POT,POTOUT-POT)
        ELSE
          POT=POT+1.D0*(POTOUT-POT)
          POT(1,1)=POTOUT(1,1)
        END IF
      ENDDO      
      IF(TBROYDEN)CALL BROYDEN$CLEAR()
      IF(.NOT.TCONV) THEN
        CALL ERROR$MSG('LOOP NOT CONVERGED')
        CALL ERROR$STOP('SCF')
      END IF
!
!     ==========================================================================
!     == wrap up                                                              ==
!     ==========================================================================
      RETURN
      END SUBROUTINE SCF
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DEFORM$ETOT(GID,NR,LMRX,LMX,lmrepx,NBG,AEZ &
     &                         ,FBG,EBG,PHI,TPHI,ETOT,EB,EKIN,EH,EXC)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: LMX
      INTEGER(4),INTENT(IN) :: LMrepX
      INTEGER(4),INTENT(IN) :: NBG
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: EBG(NBG)
      REAL(8)   ,INTENT(IN) :: FBG(NBG)
      REAL(8)   ,INTENT(IN) :: PHI(NR,LMX,NBG)
      REAL(8)   ,INTENT(IN) :: TPHI(NR,LMX,NBG)
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,INTENT(OUT):: EB
      REAL(8)   ,INTENT(OUT):: EKIN
      REAL(8)   ,INTENT(OUT):: EH
      REAL(8)   ,INTENT(OUT):: EXC
      REAL(8)               :: POTOUT(NR,LMRX)
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)               :: RHO1(NR,LMRX)
      REAL(8)               :: R(NR)
      REAL(8)               :: CG
      INTEGER(4)            :: LM,LM1,LM2,LM3,IBG,ir
      INTEGER(4)            :: L
      REAL(8)               :: pi,y0
!     ****************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      CALL RADIAL$R(GID,NR,R)
!
!     ============================================================================
!     == BANDSTRUCTURE ENERGY: EB
!     ============================================================================
      EB=SUM(EBG*FBG)
!
!     ============================================================================
!     == CALCULATE DENSITY            
!     ============================================================================
      RHO1(:,:)=0.D0
      DO LM1=1,LMX
        DO LM2=1,LMX
          DO LM3=1,LMRX
            CALL CLEBSCH(LM1,LM2,LM3,CG)
            DO IBG=1,NBG
              RHO1(:,LM3)=RHO1(:,LM3)+CG*FBG(IBG)*PHI(:,LM1,IBG)*PHI(:,LM2,IBG)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ============================================================================
!     == CALCULATE KINETIC ENERGY            
!     ============================================================================
      AUX(:)=0.D0
      DO IBG=1,NBG
        DO LM=1,LMX
          AUX(:)=AUX(:)+FBG(IBG)*PHI(:,LM,IBG)*TPHI(:,LM,IBG)
        ENDDO
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAl(GID,NR,AUX,ekin)
!
!     ============================================================================
!     == EXCHANGE AND CORRELATION ENERGY: EXC
!     ============================================================================
      CALL AUGMENTATION_XC(GID,NR,LMRX,1,RHO1,EXC,POTOUT,aux)
!     == cut out low density region ==========================================
      do ir=1,nr
        if(rho1(ir,1)*y0.lt.1.d-6) potout(ir,:)=0.d0
      enddo
!
!     ============================================================================
!     == COULOMB ENERGY: EH
!     ============================================================================
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AUX)
      POTOUT(:,1)=POTOUT(:,1)+AUX(:)
      AUX1(:)=AUX(:)*RHO1(:,1)
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1)+1.D-10))
        CALL RADIAL$POISSON(GID,NR,L,RHO1(:,LM),AUX)
        POTOUT(:,LM)=POTOUT(:,LM)+AUX(:)
        AUX1(:)=AUX1+0.5D0*AUX(:)*RHO1(:,LM)
      ENDDO
      AUX1(:)=AUX1(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX1,EH)
!
!     ============================================================================
!     == TOTAL ENERGY
!     ============================================================================
      ETOT=ekin+EH+EXC
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NONSCF(GID,NR,LMRX,POT,DREL,NB,LOFI,NNOFI,NBG,EBG,LX,PHI,TPHI)
!     **************************************************************************
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID 
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      REAL(8)   ,INTENT(IN) :: POT(NR,LMRX)
      REAL(8)   ,INTENT(IN) :: DREL(NR)
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      INTEGER(4),INTENT(IN) :: NNOFI(NB)
      INTEGER(4),INTENT(IN) :: NBG
      INTEGER(4),INTENT(IN) :: LX
      REAL(8)   ,INTENT(OUT):: EBG(NBG)
      REAL(8)   ,INTENT(OUT):: PHI(NR,(LX+1)**2,NBG)
      REAL(8)   ,INTENT(OUT):: TPHI(NR,(LX+1)**2,NBG)
      REAL(8)   ,PARAMETER  :: RAD=5.D0      
      INTEGER(4),PARAMETER  :: NITER=5
      REAL(8)               :: EBGPREV(NBG)
      REAL(8)               :: R(NR)
      REAL(8)               :: G(NR,(LX+1)**2)
      INTEGER(4)            :: LMX
      REAL(8)               :: EOFI(NB)
      REAL(8)               :: ENU
      REAL(8)               :: VAL
      LOGICAL(4)            :: TMAINSH(LX+1)
      INTEGER(4)            :: NPHI
      REAL(8)               :: PHI0(NR),AUX(NR)
      LOGICAL(4)            :: TOK
      INTEGER(4)            :: IB,ITER,L,I1,I2,IBG,LM,M,I,J,J1,J2
      REAL(8)               :: EMIN1,EMAX1,EMIN2,EMAX2
      INTEGER(4)            :: NCOLOR
      INTEGER(4)            :: ICOLOR(NBG)
      INTEGER(4)            :: LOFBG(NBG)
      REAL(8)               :: SPACING=0.5D0
      REAL(8)               :: ESPACE=0.1D0
      REAL(8)               :: ENUMAX=-0.3D0
      CHARACTER(64)         :: STRING
      REAL(8)  ,ALLOCATABLE :: EBG1(:)
      REAL(8)  ,ALLOCATABLE :: PHI1(:,:,:)
      REAL(8)  ,ALLOCATABLE :: TPHI1(:,:,:)
      REAL(8)  ,ALLOCATABLE :: OVER(:,:)
      LOGICAL(4),ALLOCATABLE :: TMAP(:)
      LOGICAL(4)            :: TDO
!     **************************************************************************
      IF(NBG.NE.SUM(2*LOFI(:)+1)) THEN 
        CALL ERROR$STOP('NBG AND LOFI ARE INCONSISTENT')
        CALL ERROR$STOP('NONSCF')
      END IF
      LMX=(LX+1)**2
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == DETERMINE STARTING ENERGIES FROM SPHERICAL POTENTIAL                 ==
!     ==========================================================================
      EOFI(:)=0.D0
      DO IB=1,NB
        G(:,:)=0.D0
        CALL BOUNDSTATE(GID,NR,LOFI(IB),0,RAD,DREL,G,NNOFI(IB),POT(:,1),EOFI(IB),PHI0)
        IF(IB.LT.NB)EOFI(IB+1)=EOFI(IB)
      ENDDO
!
!     ==========================================================================
!     == DETERMINE STARTING ENERGIES FROM SPHERICAL POTENTIAL                 ==
!     ==========================================================================
      IBG=0
      DO IB=1,NB
        EBG(IBG+1:IBG+2*LOFI(IB)+1)=EOFI(IB)
        IBG=IBG+2*LOFI(IB)+1
      ENDDO
      DO ITER=1,NITER
!       == START COLORS ========================================================
!       == COLORS IDENTIFY GROUPS OF STATES ====================================
        NCOLOR=0
        ICOLOR(:)=0
        IBG=0
        DO IB=1,NB
          L=LOFI(IB)
          NCOLOR=NCOLOR+1
          DO M=1,2*L+1
            IBG=IBG+1
            ICOLOR(IBG)=NCOLOR
            LOFBG(IBG)=L
          ENDDO
        ENDDO
!       == JOIN GROUPS =========================================================
        DO I=1,NCOLOR
          EMIN1=MINVAL(EBG,ICOLOR(:).EQ.I)
          EMAX1=MAXVAL(EBG,ICOLOR(:).EQ.I)
          DO J=I+1,NCOLOR
            EMIN2=MINVAL(EBG,ICOLOR(:).EQ.J)
            EMAX2=MAXVAL(EBG,ICOLOR(:).EQ.J)
            IF(EMIN2.LT.EMAX1+SPACING) THEN
              DO IBG=1,NBG
                IF(ICOLOR(IBG).EQ.J)ICOLOR(IBG)=I
              ENDDO
            END IF
          ENDDO
        ENDDO
!       =================================================================
        EBGPREV(:)=EBG(:)
        IBG=0
        DO I=1,NCOLOR
          TMAINSH(:)=.FALSE.
          NPHI=0
          DO J=1,NBG
            IF(ICOLOR(J).NE.I) CYCLE
            NPHI=NPHI+1
            TMAINSH(LOFBG(J)+1)=.TRUE.
          ENDDO
          IF(NPHI.EQ.0) CYCLE
          ALLOCATE(EBG1(NPHI))
          ALLOCATE(PHI1(NR,LMX,NPHI))
          ALLOCATE(TPHI1(NR,LMX,NPHI))
          ALLOCATE(OVER(NPHI,NPHI))
          ALLOCATE(TMAP(NPHI))
          I1=IBG+1
          DO I2=IBG+1,IBG+NPHI
            TDO=I2.EQ.IBG+NPHI
            IF(I2.LT.IBG+NPHI)TDO=TDO.OR.ABS(EBG(I2+1)-EBG(I2)).GT.ESPACE
            IF(TDO) THEN
              ENU=SUM(EBG(I1:I2))/REAL(I2-I1+1)
              ENU=MIN(ENU,ENUMAX)
              G(:,:)=0.D0
              CALL RADIAL$NONSPHBOUND_NONSO(GID,NR,LX,LMX,LMRX,TMAINSH,POT,DREL,G,ENU &
      &                  ,NPHI,EBG1,PHI1,TPHI1,TOK)
!
!             == DETERMINE OVERLAP OF PREVOUSLY CALCULATED STATES WITH THIS SET= 
              OVER(:,:)=0.D0
              DO J1=IBG+1,I1-1
                DO J2=1,NPHI
                  AUX(:)=0.D0
                  DO LM=1,LMX
                    AUX(:)=AUX(:)+PHI(:,LM,J1)*PHI1(:,LM,J2)
                  ENDDO
                  AUX(:)=AUX(:)*R(:)**2
                  CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
                  OVER(J1-IBG,J2)=VAL
                ENDDO
              ENDDO    
!
!             == TMAP IDENTIFIES STATE THAT HAVE BEEN CALCULATED PREVIOUSLY   ==
              TMAP(:)=.FALSE.
              DO J1=IBG+1,I1-1
                TMAP(MAXLOC(ABS(OVER(J1-IBG,:))))=.TRUE.
              ENDDO
!
!             == MAP STATES OF THIS SUBSET INTO GLOBAL ARRAYS ==================
              J1=I1
              DO J2=1,NPHI
                IF(TMAP(J2)) CYCLE
                EBG(J1)=EBG1(J2)
                PHI(:,:,J1)=PHI1(:,:,J2)
                TPHI(:,:,J1)=TPHI1(:,:,J2)
                IF(J1.GE.I2) EXIT
                J1=J1+1
              ENDDO
              I1=I2+1
            END IF
          ENDDO
          DEALLOCATE(EBG1)
          DEALLOCATE(PHI1)
          DEALLOCATE(TPHI1)
          DEALLOCATE(OVER)
          DEALLOCATE(TMAP)
          IBG=IBG+NPHI
        ENDDO
      ENDDO
!
!     ==========================================================================        
!     == CHECK ORTHOGONALITY                                                  ==        
!     ==========================================================================        
      DO I1=1,NBG
        DO I2=I1,NBG
          AUX(:)=0.D0
          DO LM=1,LMX
            AUX(:)=AUX(:)+PHI(:,LM,I1)*PHI(:,LM,I2)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
          IF(I1.EQ.I2)VAL=VAL-1.D0
IF(ABS(VAL).GT.1.D-3)PRINT*,'0VERLAP ',I1,I2,VAL
          IF(ABS(VAL).GT.3.D-1) THEN
            WRITE(*,FMT='("OVERLAP ",2I5,F15.10)')I1,I2,VAL
            WRITE(*,FMT='("EBG   :",10F10.5)')EBG
            WRITE(*,FMT='("ICOLOR:",10I10)')ICOLOR
            DO IBG=1,NBG
              WRITE(STRING,*)IBG
              STRING='PHI'//TRIM(ADJUSTL(STRING))//'.DAT'
              CALL WRITEPHI(STRING,GID,NR,LMX,PHI(:,:,IBG))
            ENDDO
            CALL ERROR$MSG('DEVIATION FROM ORTHONORMALITY LARGER THAN PERMITTED')
            CALL ERROR$STOP('NONSCF')
          END IF
        ENDDO
      ENDDO
      RETURN
      END
!!$!..................................................................................
!!$MODULE BROYDEN_MODULE
!!$LOGICAL(4)         :: TON=.FALSE.
!!$INTEGER(4)         :: NSTEPX=0
!!$INTEGER(4)         :: NSTEP=0
!!$INTEGER(4)         :: NX=0
!!$REAL(8)            :: ALPHA
!!$REAL(8),ALLOCATABLE :: XPREV(:,:)
!!$REAL(8),ALLOCATABLE :: YPREV(:,:)
!!$END MODULE BROYDEN_MODULE
!!$!      .............................................................................
!!$       SUBROUTINE BROYDEN$NEW(NX_,NSTEPX_,ALPHA_)
!!$       USE BROYDEN_MODULE
!!$       IMPLICIT NONE
!!$       INTEGER(4),INTENT(IN)    :: NX_
!!$       INTEGER(4),INTENT(IN)    :: NSTEPX_
!!$       REAL(8)   ,INTENT(IN)    :: ALPHA_
!!$!      *****************************************************************************
!!$       IF(TON) THEN
!!$         CALL ERROR$MSG('BROYDEN OBJECT ALREADY IN USE')
!!$         CALL ERROR$STOP('BROYDEN$NEW')
!!$       END IF
!!$       TON=.TRUE.
!!$       NSTEP=0
!!$       NX=NX_
!!$       NSTEPX=NSTEPX_
!!$       ALPHA=ALPHA_
!!$       ALLOCATE(XPREV(NX,NSTEPX))
!!$       ALLOCATE(YPREV(NX,NSTEPX))
!!$       RETURN
!!$       END
!!$!      .............................................................................
!!$       SUBROUTINE BROYDEN$CLEAR
!!$       USE BROYDEN_MODULE
!!$       IMPLICIT NONE
!!$!      *****************************************************************************
!!$       IF(.NOT.TON) THEN
!!$         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
!!$         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
!!$         CALL ERROR$STOP('BROYDEN$CLEAR')
!!$       END IF
!!$       TON=.FALSE.
!!$       NSTEPX=0
!!$       NSTEP=0
!!$       NX=0
!!$       ALPHA=0.D0
!!$       DEALLOCATE(XPREV)
!!$       DEALLOCATE(YPREV)
!!$       RETURN
!!$       END
!!$!      .............................................................................
!!$       SUBROUTINE BROYDEN$STEP(NX_,X,Y)
!!$       USE BROYDEN_MODULE
!!$       IMPLICIT NONE
!!$       INTEGER(4),INTENT(IN)    :: NX_
!!$       REAL(8)   ,INTENT(INOUT) :: X(NX_)
!!$       REAL(8)   ,INTENT(IN)    :: Y(NX_)
!!$       REAL(8)   ,ALLOCATABLE   :: DX(:,:)
!!$       REAL(8)   ,ALLOCATABLE   :: DY(:,:)
!!$       REAL(8)   ,ALLOCATABLE   :: B(:,:)
!!$       REAL(8)   ,ALLOCATABLE   :: BINV(:,:)
!!$ REAL(8)   ,ALLOCATABLE   :: W(:,:)
!!$       INTEGER(4)               :: I
!!$!      *****************************************************************************
!!$       IF(.NOT.TON) THEN
!!$         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
!!$         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
!!$         CALL ERROR$STOP('BROYDEN$STEP')
!!$       END IF
!!$       IF(NX_.NE.NX) THEN
!!$         CALL ERROR$MSG('SIZE INCONSISTENT')
!!$         CALL ERROR$STOP('BROYDEN$STEP')
!!$       END IF
!!$!PRINT*,'NSTEP',NSTEP
!!$!
!!$!      =================================================================
!!$!      == SIMPLE MIXING IN THE FIRST STEP                             ==
!!$!      =================================================================
!!$       IF(NSTEP.EQ.0) THEN
!!$         IF(NSTEPX.GT.0)THEN
!!$           NSTEP=1
!!$           XPREV(:,1)=X(:)     
!!$           YPREV(:,1)=Y(:)     
!!$         END IF
!!$         X=X+ALPHA*Y
!!$         RETURN
!!$       END IF
!!$!
!!$!      =================================================================
!!$!      == DETERMINE INVERSE HESSIAN ALPHA+DX OTIMES DY                ==
!!$!      =================================================================
!!$       ALLOCATE(DX(NX,NSTEP))
!!$       ALLOCATE(DY(NX,NSTEP))
!!$       DO I=1,NSTEP
!!$         DY(:,I)=YPREV(:,I)-Y(:)  
!!$         DX(:,I)=XPREV(:,I)-X(:)+ALPHA*DY(:,I)
!!$       ENDDO
!!$       ALLOCATE(B(NSTEP,NSTEP))
!!$       ALLOCATE(BINV(NSTEP,NSTEP))
!!$       B=MATMUL(TRANSPOSE(DY),DY)   !OVERLAP MATRIX OF DY
!!$!PRINT*,'B',B
!!$       CALL LIB$INVERTR8(NSTEP,B,BINV)           
!!$ALLOCATE(W(NX,NSTEP))
!!$W=MATMUL(DY,BINV)            !NEW DY IS BIORTHONORMAL TO OLD DY
!!$!PRINT*,'W ',MATMUL(TRANSPOSE(W),DY)
!!$DEALLOCATE(W)
!!$       DY=MATMUL(DY,BINV)            !NEW DY IS BIORTHONORMAL TO OLD DY
!!$       DEALLOCATE(B)
!!$       DEALLOCATE(BINV)
!!$!
!!$!      =================================================================
!!$!      == STORE HISTORY                                               ==
!!$!      =================================================================
!!$       IF(NSTEP.LT.NSTEPX)NSTEP=NSTEP+1
!!$       DO I=NSTEP,2,-1
!!$         YPREV(:,I)=YPREV(:,I-1)
!!$         XPREV(:,I)=XPREV(:,I-1)
!!$       ENDDO
!!$       XPREV(:,1)=X(:)     
!!$       YPREV(:,1)=Y(:)     
!!$!
!!$!      =================================================================
!!$!      == PREDICT NEW VECTOR                                          ==
!!$!      =================================================================
!!$       X=X+ALPHA*Y-MATMUL(DX,MATMUL(TRANSPOSE(DY),Y))
!!$       DEALLOCATE(DX)
!!$       DEALLOCATE(DY)
!!$       RETURN
!!$       END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONEATOM$NEW(IAT,AEZ)
!     **************************************************************************
!     ** INITIALIZES AN ATOM                                                  **
!     **************************************************************************
      USE ONEATOM_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      REAL(8)   ,INTENT(IN)  :: AEZ          ! ATOMIC NUMBER
      REAL(8)   ,PARAMETER   :: R1=1.056D-4  ! INNERMOST POINT OF RADIAL GRID
      REAL(8)   ,PARAMETER   :: DEX=0.05D0   ! LOG SPACING OF RADIAL GRID
      INTEGER(4),PARAMETER   :: NR=250       ! #(RADIAL GRID POINTS)
      INTEGER(4)             :: GID          ! GRID ID
      INTEGER(4)             :: NS(10),NP(10),ND(10),NF(10)
      INTEGER(4)             :: SOFI(19)
      INTEGER(4)             :: LOFI(19)
      INTEGER(4)             :: NNOFI(19)
      REAL(8)                :: FOFI(19)
      REAL(8)                :: EOFI(19)
      REAL(8)                :: RBOX=20.D0
      REAL(8)                :: EV
      REAL(8)                :: SVAR
      INTEGER(4)             :: NB       ! #(STATES IN OCCUPIED SHELLS)
      INTEGER(4)             :: NC       ! #(CORE STATES)
      INTEGER(4)             :: L,IB,ISVAR
      REAL(8)   ,ALLOCATABLE :: AEPOT(:)
      REAL(8)   ,ALLOCATABLE :: DREL(:)
      REAL(8)   ,ALLOCATABLE :: RHOADD(:)
      REAL(8)   ,ALLOCATABLE :: RHOCORE(:)
      REAL(8)   ,ALLOCATABLE :: RHOVALENCE(:)
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)   ,ALLOCATABLE :: PHI(:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)                :: ETOT,DEDQ,DEDRAD,EKINC
      REAL(8)                :: eband,epot,ekin,exc,eh
      REAL(8)   ,PARAMETER   :: Q=0.D0
      REAL(8)                :: RAD
      REAL(8)                :: PI,C0LL
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL CONSTANTS('EV',EV)
      IF(IAT.GT.NATX) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$STOP('ONEATOM$NEW')
      END IF
      THIS=>THISARR(IAT)
      THIS%AEZ=AEZ
!
!     ==========================================================================
!     ==  DEFINE RADIAL GRID                                                  ==
!     ==========================================================================
      CALL RADIAL$NEW('SHLOG',GID)
      THIS%GID=GID
      CALL RADIAL$SETR8(GID,'R1',R1)
      CALL RADIAL$SETR8(GID,'DEX',DEX)
      CALL RADIAL$SETI4(GID,'NR',NR)
      THIS%NR=NR
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  DEFINE BANDS                                                        ==
!     ==========================================================================
      LOFI( 1)=0   ! 1S
      LOFI( 2)=0   ! 2S
      LOFI( 3)=1   ! 2P
      LOFI( 4)=0   ! 3S
      LOFI( 5)=1   ! 3P
      LOFI( 6)=0   ! 4S
      LOFI( 7)=2   ! 3D
      LOFI( 8)=1   ! 4P
      LOFI( 9)=0   ! 5S
      LOFI(10)=2   ! 4D
      LOFI(11)=1   ! 5P
      LOFI(12)=0   ! 6S
      LOFI(13)=3   ! 4F
      LOFI(14)=2   ! 5D
      LOFI(15)=1   ! 6P
      LOFI(16)=0   ! 7S
      LOFI(17)=3   ! 5F
      LOFI(18)=2   ! 6D
      LOFI(19)=1   ! 7P
!
!     == DETERMINE NUMBER OF NODES ===================================
      DO L=0,3
        ISVAR=0
        DO IB=1,19
          IF(LOFI(IB).NE.L) CYCLE
          NNOFI(IB)=ISVAR
          ISVAR=ISVAR+1
        ENDDO
      ENDDO
!
!     == OCCUPATIONS  AND NUMBER OF BANDS===============================
      FOFI(:)=0.D0
      SVAR=AEZ
      DO IB=1,19
        FOFI(IB)=MIN(SVAR,REAL(2*(2*LOFI(IB)+1),KIND=8))
        SVAR=SVAR-FOFI(IB)
      ENDDO
!
!     == DETERMINE THE NUMBER OF CORE STATES =============================
      NB=19
      DO IB=1,19
        IF(LOFI(IB).NE.0) CYCLE
        IF(FOFI(IB).NE.0.D0) NC=IB-1
        IF(FOFI(IB).EQ.0.D0) THEN
          NB=IB-1   ! INCLUDE AT LEAST ONE EMPTY ORBITAL TO FACILITATE FERMI LEVEL SEARCH
          EXIT
        END IF
      ENDDO
!
!     == MAP ON THIS ==================================================
      THIS%NB=NB
      THIS%NC=NC
      THIS%LOFI=LOFI
      THIS%NNOFI=NNOFI
      THIS%FOFI=FOFI
!
!     ==========================================================================
!     ==  DETERMINE ATOMIC POTENTIAL                                          ==
!     ==========================================================================
PRINT*,'AEZ',AEZ
PRINT*,'CORE'
DO IB=1,NC
  PRINT*,IB,LOFI(IB),NNOFI(IB),FOFI(IB)
ENDDO
PRINT*,'VALENCE'
DO IB=NC+1,NB
  PRINT*,IB,LOFI(IB),NNOFI(IB),FOFI(IB)
ENDDO
      ALLOCATE(AEPOT(NR))
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AEPOT)
      ALLOCATE(DREL(NR))
      DREL(:)=0.D0
      ALLOCATE(RHOADD(NR))
      RHOADD=0.D0
      SOFI(:)=0
      EOFI(:)=0.D0
      Rbox=R(NR-3)
      CALL AESCF(GID,NR,NB,TREL,LOFI,SOFI,FOFI,NNOFI,AEZ,RHOADD,RBOX,DREL,AEPOT,EOFI)
!
!     ==========================================================================
!     == CORE DENSITY
!     ==========================================================================
      ALLOCATE(PHI(NR))
      ALLOCATE(AUX(NR))
      ALLOCATE(RHOCORE(NR))
      RHOCORE(:)=0.D0
      DO IB=1,NC
        AUX(:)=0.D0 !INHOMOGENITY
        CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RBOX,DREL,AUX,NNOFI(IB),AEPOT,EOFI(IB),PHI)
        AUX(:)=(R(:)*PHI(:))**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        PHI(:)=PHI(:)/SQRT(SVAR)
        RHOCORE(:)=RHOCORE(:)+FOFI(IB)*C0LL*PHI(:)**2
      ENDDO
!
!     ==========================================================================
!     == valence density DENSITY
!     ==========================================================================
      ALLOCATE(RHOVALENCE(NR))
      rhovalence(:)=0.d0
      DO IB=NC+1,NB
        AUX(:)=0.D0 !INHOMOGENITY
        CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RBOX,DREL,AUX,NNOFI(IB),AEPOT,EOFI(IB),PHI)
        AUX(:)=(R(:)*PHI(:))**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        PHI(:)=PHI(:)/SQRT(SVAR)
        RHOVALENCE(:)=RHOVALENCE(:)+FOFI(IB)*C0LL*PHI(:)**2
      ENDDO
      DEALLOCATE(PHI)
      DEALLOCATE(AUX)
!
!     ==========================================================================
!     == CORE KINETIC ENERGY                                                  ==
!     ==========================================================================
      ALLOCATE(AUX(NR))
      AUX(:)=R(:)**2*AEPOT(:)*RHOCORE(:)
      CALL RADIAL$INTEGRAL(GID,NR,AUX,EKINC)
      EKINC=-EKINC
      DO IB=1,NC
        EKINC=EKINC+FOFI(IB)*EOFI(IB)
      ENDDO
      DEALLOCATE(AUX)
!
!     ==========================================================================
!     == determine total energy                                               ==
!     ==========================================================================
      ALLOCATE(AUX(NR))
      CALL MYVOFRHO(GID,NR,aeZ,RHOcore+rhovalence,aux,EH,EXC)
      DEALLOCATE(AUX)
      eband=sum(eofi*fofi)
      call radial$integral(gid,nr,r(:)**2*aepot(:)*(rhocore(:)+rhovalence(:)),epot)
      ekin=eband-epot
      etot=ekin+eh+exc
!
!     ==========================================================================
!     == map data onto "this"                                                 ==
!     ==========================================================================
      THIS%EREF=ETOT
      THIS%EKINC=EKINC
      THIS%EOFI=1.D0
      THIS%EOFI(1:NB)=EOFI(1:NB)
      ALLOCATE(THIS%atomPOT(NR))
      THIS%atomPOT(:)=AEPOT(:)
      DEALLOCATE(AEPOT)
      ALLOCATE(THIS%RHOCORE(NR))
      THIS%RHOCORE(:)=RHOCORE(:)
      DEALLOCATE(RHOCORE)
      ALLOCATE(THIS%DREL(NR))
      THIS%DREL=DREL
      DEALLOCATE(DREL)
!
!     ==========================================================================
!     == REPORT ================================================================
!     ==========================================================================
      WRITE(*,FMT='(72("="))')
      WRITE(*,FMT='("ISOLATED ATOM IN A BOX WITH RADIUS ",F5.2," A.U")')RBOX
      WRITE(*,FMT='("EIGENSTATES AFTER NON-SCF CALCULATION IN SCF POT")')
      DO IB=1,NB
        WRITE(*,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
     &                ," E[H]=",F20.9," E[EV]=",F20.9)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/EV
      ENDDO
      write(*,fmt='("aescf  eband ",f10.5)')eband
      write(*,fmt='("aescf  epot  ",f10.5)')epot
      write(*,fmt='("aescf  ekin  ",f10.5)')ekin
      write(*,fmt='("aescf  eh    ",f10.5)')eh
      write(*,fmt='("aescf  exc   ",f10.5)')exc
      write(*,fmt='("aescf  etot  ",f10.5)')etot
      WRITE(*,FMT='(72("="))')
!
!     ==========================================================================
!     ==  calculate atom in a large box as reference                          ==
!     ==========================================================================
!!$      RAD=R(NR-3)
!!$      CALL ONEATOM(IAT,Q,RAD,ETOT,DEDQ,DEDRAD)
!!$!
!!$!     ==========================================================================
!!$!     == REPORT                                                               ==
!!$!     ==========================================================================
!!$      WRITE(*,FMT='(72("="))')
!!$      WRITE(*,FMT='("ISOLATED ATOM IN A BOX WITH RADIUS ",F5.2," A.U")')THIS%RAD
!!$      WRITE(*,FMT='("EIGENSTATES AT TEMPERATURE 1.D-3 H")')
!!$      EOFI(:)=THIS%EOFI(:)
!!$      FOFI(:)=THIS%FOFI(:)
!!$      NNOFI(:)=THIS%NNOFI(:)
!!$      WRITE(*,FMT='("TOTAL ENERGY ",T50,F10.5," H")')etot
!!$      DO IB=1,NB
!!$        WRITE(*,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
!!$     &                ," E[H]=",F20.9," E[EV]=",F20.9)') &
!!$     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/EV
!!$      ENDDO
!!$      WRITE(*,FMT='(72("="))')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONEATOM(IAT,Q,RAD,ETOT,DEDQ,DEDRAD)
!     **************************************************************************
!     **  THE TEMPERATURE MUST NOT BE TOO SMALL. OTHERWISE DEDQ CANNOT BE     **
!     **  ACCURATELY DETERMINED                                               **
!     **************************************************************************
      USE ONEATOM_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      REAL(8)   ,INTENT(IN) :: Q
      REAL(8)   ,INTENT(IN) :: RAD
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,INTENT(OUT):: DEDQ
      REAL(8)   ,INTENT(OUT):: DEDRAD
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: NB
      INTEGER(4)            :: NC
      REAL(8)               :: AEZ
      REAL(8)               :: QVALENCE
      INTEGER(4)            :: LOFI(19)
      INTEGER(4)            :: SOFI(19)
      INTEGER(4)            :: NNOFI(19)
      REAL(8)               :: EOFI(19)
      REAL(8)               :: FOFI(19)
      REAL(8)   ,ALLOCATABLE:: AEPOT(:)
      REAL(8)   ,ALLOCATABLE:: POTIN(:)
      REAL(8)   ,ALLOCATABLE:: POTh(:)
      REAL(8)   ,ALLOCATABLE:: POTxc(:)
      REAL(8)   ,ALLOCATABLE:: RHOADD(:)
      REAL(8)   ,ALLOCATABLE:: RHO(:)
      REAL(8)   ,ALLOCATABLE:: DREL(:)
      REAL(8)   ,ALLOCATABLE:: PHI(:,:)
      REAL(8)   ,ALLOCATABLE:: PHIDOT(:)
      REAL(8)   ,ALLOCATABLE:: G(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:),AUX1(:)
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)   ,ALLOCATABLE:: ETAPAULI(:,:)
      REAL(8)   ,ALLOCATABLE:: Etacpauli(:,:)
      REAL(8)   ,ALLOCATABLE:: vcpauli(:,:)
      REAL(8)   ,ALLOCATABLE:: ocpauli(:,:)
      REAL(8)   ,ALLOCATABLE:: UDOT(:,:)
      REAL(8)   ,ALLOCATABLE:: UdDOT(:,:)
      REAL(8)   ,ALLOCATABLE:: UN(:,:)
      REAL(8)               :: SVAR,DER,VAL,DERDOT,VALDOT
      REAL(8)               :: PI,C0LL,y0
      REAL(8)               :: EKIN,EH,EXC
      integer(4)            :: ibi,istart
      REAL(8)               :: x0,z0,dx,xm,zm,e
      REAL(8)               :: XMAX
      LOGICAL(4)            :: CONVG
      INTEGER(4)            :: IB,IR,ITER,IRBOX,IL,L,LMAX
      INTEGER(4),PARAMETER  :: NITER=100
      REAL(8)   ,PARAMETER  :: TOL=1.D-5
      REAL(8)   ,PARAMETER  :: KBT=1.D-3
      REAL(8)               :: TS,F
INTEGER(4) :: I,J
REAL(8)    :: DEG,SVAR1
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      C0LL=1.D0/SQRT(4.D0*PI)
      y0=1.D0/SQRT(4.D0*PI)
      THIS=>THISARR(IAT)
      GID=THIS%GID
      NR=THIS%NR
      AEZ=THIS%AEZ      
      NB=THIS%NB
      NC=THIS%NC
      LOFI(:)=THIS%LOFI(:)
      NNOFI(:)=THIS%NNOFI(:)
      EOFI(:)=THIS%EOFI(:)
      FOFI(:)=THIS%FOFI(:)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GE.RAD) EXIT
      ENDDO
!
!     ==========================================================================
!     ==  OBTAIN SELF-CONSISTENT POTENTIAL                                    ==
!     ==========================================================================
      SOFI(:)=0.D0
      ALLOCATE(AEPOT(NR));       AEPOT(:)=THIS%atomPOT(:)
      ALLOCATE(POTIN(NR));       POTIN(:)=AEPOT(:)
      ALLOCATE(RHOADD(NR));      RHOADD=0.D0
      ALLOCATE(DREL(NR));        DREL=THIS%DREL
      ALLOCATE(RHO(NR))
      ALLOCATE(PHI(NR,19))
      ALLOCATE(PHIDOT(NR))
      ALLOCATE(G(NR))
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
      ALLOCATE(poth(NR))
      ALLOCATE(potxc(NR))
      XMAX=0.D0
      CONVG=.FALSE.
      CALL BROYDEN$NEW(NR,4,1.D0)
      POTIN=AEPOT
      DO ITER=1,NITER
        DO IB=NC+1,NB
          G(:)=0.D0
          CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RAD,DREL,G,NNOFI(IB),POTIN,EOFI(IB),PHI(:,IB))
          AUX(:)=(R(:)*PHI(:,IB))**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1(:),RAD,SVAR)
          PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
        ENDDO
!
!       == DETERMINE OCCUPATIONS ===============================================
        QVALENCE=AEZ+Q-REAL(SUM(2*(2*LOFI(1:NC)+1)),KIND=8)
        CALL ONEATOM_OPTF(NB-NC,QVALENCE,KBT,LOFI(NC+1:NB),EOFI(NC+1:NB),FOFI(NC+1:NB),DEDQ)
!
!       == DETERMINE DENSITY ===================================================
        RHO(:)=THIS%RHOCORE(:)
        DO IB=NC+1,NB
          RHO(:)=RHO(:)+FOFI(IB)*C0LL*PHI(:,IB)**2
        ENDDO
!
!       == DETERMINE OUTPUT POTENTIAL ==========================================
        CALL MYVOFRHOWITHRAD(GID,NR,RAD,AEZ,RHO,POTh,potxc,EH,EXC)
        aepot=poth+potxc
!
!       == EXIT IF CONVERGED ===================================================
        IF(CONVG) THEN
          CALL BROYDEN$CLEAR
          EXIT
        END IF
!
!       == POTENTIAL MIXING ====================================================
!       == THE POTENTIAL SIS SET EQUAL BEYOND THE BOX RADIUS, BECAUSE OTHERWISE
!       == THIS REGION DETERMINES THE CONVERGENCE ==============================
        POTIN(IRBOX:)=AEPOT(IRBOX:)
        XMAX=MAXVAL(ABS(AEPOT-POTIN)) 
        CALL BROYDEN$STEP(NR,POTIN,(AEPOT-POTIN))
        CONVG=(XMAX.LT.TOL)
      ENDDO
      IF(.NOT.CONVG) THEN
        CALL BROYDEN$CLEAR
        CALL ERROR$MSG('SELFCONSISTENCY LOOP NOT CONVERGED')
        CALL ERROR$STOP('ONEATOM')
      END IF
!
!     ==========================================================================
!     == CALCULATE PRESSURE                                                   ==
!     ==========================================================================
      DEDRAD=0.D0
      DO IB=NC+1,NB
        IF(FOFI(IB).EQ.0.D0) CYCLE
        CALL RADIAL$DERIVATIVE(GID,NR,PHI(:,IB),RAD,DER)
        IF(DER.EQ.0.D0) CYCLE  ! AVOID CORE WAVE FUNCTIONS
        G(:)=0.D0
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,SOFI(IB),G,LOFI(IB),EOFI(IB),1,PHI(:,IB))
! NORMALIZATION NOT REQUIRED
 AUX=R (:)**2*PHI(:,IB)**2
 CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
 CALL RADIAL$VALUE(GID,NR,AUX1(:),RAD,SVAR)
 PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
!
        CALL RADIAL$VALUE(GID,NR,PHI(:,IB),RAD,VAL)
        CALL RADIAL$DERIVATIVE(GID,NR,PHI(:,IB),RAD,DER)
        G(:)=PHI(:,IB)         
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,SOFI(IB),G,LOFI(IB),EOFI(IB),1,PHIDOT)
! == ORTHONORMALIZATION IS NOT NECESSARY =======
  AUX=R (:)**2*PHI(:,IB)*PHIDOT(:)
  CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
  CALL RADIAL$VALUE(GID,NR,AUX1(:),RAD,SVAR)
  PHIDOT(:)=PHIDOT(:)-PHI(:,IB)*SVAR
        CALL RADIAL$VALUE(GID,NR,PHIDOT,RAD,VALDOT)
        CALL RADIAL$DERIVATIVE(GID,NR,PHIDOT,RAD,DERDOT)
        DEDRAD=DEDRAD-FOFI(IB)*DER/VALDOT
      ENDDO
!
!     ==========================================================================
!     == CALCULATE TOTAL ENERGY                                               ==
!     ==========================================================================
print*,'========================  oneatom ======================='
print*,'nc ',nc,nb
print*,'fofi ',fofi(:nb)
print*,'eofi ',eofi(:nb)
      AUX(:)=(RHO(:)-THIS%RHOCORE(:))*AEPOT(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,SVAR)
      EKIN=THIS%EKINC+SUM(EOFI(nc+1:nb)*FOFI(nc+1:nb))-SVAR
      CALL MYVOFRHOWITHRAD(GID,NR,RAD,AEZ,RHO,AUX,aux1,EH,EXC)
      ETOT=EKIN+EH+EXC-THIS%EREF
print*,'ekin ',ekin
print*,'eh ',eh
print*,'exc ',exc
!     == ENTROPY TERM OF THE ELECTRONS
      TS=0.D0
      DO IB=NC+1,NB
        SVAR=REAL(2*(2*LOFI(IB)+1),KIND=8)
        F=FOFI(IB)/SVAR
        IF(F.LT.1.D-50) CYCLE
        IF(1.D0-F.LT.1.D-50) CYCLE
        TS=TS+KBT*SVAR*(F*LOG(F)+(1.D0-F)*LOG(1.D0-F))
      ENDDO   
      ETOT=ETOT+TS   
print*,'ets ',ts
print*,'etot ',etot,THIS%EREF
!
!     ==========================================================================
!     == CALCULATE REPULSIVE POTENTIAL FOR EACH ANGULAR MOMENTUM              ==
!     ==========================================================================
      LMAX=MAXVAL(LOFI(1:NB))
      ALLOCATE(ETAPAULI(NR,LMAX+1))
      ALLOCATE(Etacpauli(NR,LMAX+1))
      ALLOCATE(ocpauli(NR,LMAX+1))
      ALLOCATE(vcpauli(NR,LMAX+1))
      ALLOCATE(UDOT(NR,LMAX+1))
      ALLOCATE(UDdOT(NR,LMAX+1))
      ALLOCATE(UN(NR,LMAX+1))
      DO IL=1,LMAX+1
        L=IL-1
!
!       =========================================================================
!       == construct highes nodeless core state                                ==
!       =========================================================================
        un(:,il)=0.d0
        DO IB=1,NC
          IF(LOFI(IB).NE.L) CYCLE
          G(:)=UN(:,IL)
          CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RAD,DREL,G,0 &
     &                   ,POTIN,EOFI(IB),UN(:,IL))
print*,'eofi(ib)',eofi(ib)
        ENDDO
!
!       =========================================================================
!       == find binding energy e                                               ==
!       =========================================================================
        do i=nc+1,nb
          IF(LOFI(I).NE.L) CYCLE
          ib=i
          exit
        enddo
        ISTART=1
        X0=Eofi(ib)
        DX=-1.D-3
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
        DO I=1,NITER
          E=X0
          g(:)=un(:,il)
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POTin,DREL,SOfi(ib),G,L,E,1,PHI)
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,rad,Z0)
          Z0=Z0-0.5d0
print*,'e,z0 ',e,z0,eofi(ib)
          IF(ABS(2.D0*DX).LE.TOL) EXIT
          CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
        enddo
        IF(ABS(DX).GT.TOL) THEN
          CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
          CALL ERROR$MSG('BOUND STATE NOT FOUND')
          CALL ERROR$STOP('oneatom')
        END IF
!
!       ========================================================================
!       == determine core repulsion potential                                 ==
!       ========================================================================
        g=un(:,il)
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,SOFI(IB),G,LOFI(IB) &
     &                                                        ,E,1,UDOT(:,IL))
        g=udot(:,il)
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,SOFI(IB),G,LOFI(IB) &
     &                                                        ,E,1,UdDOT(:,IL))
!       ==  AVOID DIVIDE BY ZERO IF UDOT=0 NEAR THE NUCLEUS
        I=1
        DO IR=1,NR
          IF(UDOT(IR,IL).NE.0.D0) THEN
            I=IR
            EXIT
          END IF
        ENDDO
!       == NOW CONSTRUCT PAULI POTENTIAL
        Etacpauli(I:,IL)=-UN(I:,IL)/UDOT(I:,IL)/y0
        Etacpauli(1:I,IL)=Etacpauli(I,IL)
!       ==  AVOID DIVIDE BY ZERO IF UDOT=0 NEAR THE NUCLEUS
        I=1
        DO IR=1,NR
          IF(UDdOT(IR,IL).NE.0.D0) THEN
            I=IR
            EXIT
          END IF
        ENDDO
!       == construct ocpauli=-deta/de ==========================================
        ocpauli(I:,IL)=-etacpauli(i:,il)*Uddot(I:,IL)/UdOT(I:,IL)
        ocpauli(1:I,IL)=ocpauli(I,IL)
!       == construct vcpauli=eta+e*ocpauli ====================================
        vcpauli(:,il)=etacpauli(:,il)+e*ocpauli(:,il)
!
!       ========================================================================
!       == determine valence bound state in the box (antibonding)             ==
!       ========================================================================
        g=un(:,il)
        CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RAD,DREL,G,0 &
     &                                             ,POTIN,EOFI(IB),UN(:,IL))
print*,'energy of valence state ',l,eofi(ib),e
!
!       ========================================================================
!       == determine valence repulsion potential                              ==
!       ========================================================================
        g=un(:,il)
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,SOFI(IB),G,LOFI(IB) &
     &                                                        ,E,1,UDOT(:,IL))
!       ==  AVOID DIVIDE BY ZERO IF UDOT=0 NEAR THE NUCLEUS
        I=1
        DO IR=1,NR
          IF(UDOT(IR,IL).NE.0.D0) THEN
            I=IR
            EXIT
          END IF
        ENDDO
!       == NOW CONSTRUCT PAULI POTENTIAL
        ETAPAULI(I:,IL)=-UN(I:,IL)/UDOT(I:,IL)/y0
        ETAPAULI(1:I,IL)=ETAPAULI(I,IL)
      ENDDO
!
!     ==========================================================================
!     == MAP DATA BACK TO THIS                                                ==
!     ==========================================================================
      IF(.NOT.ASSOCIATED(THIS%ETAPAULI))ALLOCATE(THIS%ETAPAULI(NR,LMAX+1))
      ETAPAULI(IRBOX:,:)=0.D0
      THIS%ETAPAULI(:,:)=ETAPAULI(:,:)
      IF(.NOT.ASSOCIATED(THIS%ETACPAULI))ALLOCATE(THIS%ETACPAULI(NR,LMAX+1))
      ETACPAULI(IRBOX:,:)=0.D0
      THIS%ETACPAULI(:,:)=ETACPAULI(:,:)
      IF(.NOT.ASSOCIATED(THIS%VCPAULI))ALLOCATE(THIS%VCPAULI(NR,LMAX+1))
      VCPAULI(IRBOX:,:)=0.D0
      THIS%VCPAULI(:,:)=VCPAULI(:,:)
      IF(.NOT.ASSOCIATED(THIS%OCPAULI))ALLOCATE(THIS%OCPAULI(NR,LMAX+1))
      OCPAULI(IRBOX:,:)=0.D0
      THIS%OCPAULI(:,:)=OCPAULI(:,:)
      IF(.NOT.ASSOCIATED(THIS%UDOT))ALLOCATE(THIS%UDOT(NR,LMAX+1))
      THIS%UDOT=UDOT
      IF(.NOT.ASSOCIATED(THIS%UN))ALLOCATE(THIS%UN(NR,LMAX+1))
      THIS%UN=UN
      IF(.NOT.ASSOCIATED(THIS%RHO))ALLOCATE(THIS%RHO(NR))
      THIS%RHO=0.D0
      THIS%RHO(1:IRBOX-1)=RHO(1:IRBOX-1)
      IF(.NOT.ASSOCIATED(THIS%aepot))ALLOCATE(this%AEPOT(NR))
      THIS%AEPOT(:)=0.D0
      THIS%AEPOT(:)=POTIN(1:IRBOX-1)
      if(.not.associated(this%poth))allocate(this%poth(nr))
      THIS%POTH(:)=poth(:)
      if(.not.associated(this%potxc))allocate(this%potxc(nr))
      THIS%POTxc(:)=potxc(:)

      THIS%EOFI(:)=EOFI(:)
      THIS%FOFI(:)=FOFI(:)
      THIS%RAD=RAD
      THIS%Q=Q
      THIS%DEDRAD=DEDRAD
      THIS%DEDQ=DEDQ
!
!     ==========================================================================
!     == wrap up                                                              ==
!     ==========================================================================
      DEALLOCATE(ETAPAULI)
      DEALLOCATE(Etacpauli)
      DEALLOCATE(vcpauli)
      DEALLOCATE(ocpauli)
      DEALLOCATE(UN)
      DEALLOCATE(UDOT)
      DEALLOCATE(UDdOT)
      DEALLOCATE(poth)
      DEALLOCATE(potxc)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONEATOM_OPTF(N,Q,KBT,L,E,F,MU)
!     **************************************************************************
!     **  DETERMINE OCCUPATIONS FOR A GIVEN CHARGE AND TEMPERATURE            **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N       ! #(ANGULAR MOMENTUM SHELLS)
      REAL(8)   ,INTENT(IN) :: Q       ! TOTAL NUMBER OF ELECTRONS
      REAL(8)   ,INTENT(IN) :: KBT     ! TEMPERATURE AS K_B*T
      INTEGER(4),INTENT(IN) :: L(N)    ! MAIN ANGULAR MOMENTUM QUANTUM NUMBER
      REAL(8)   ,INTENT(IN) :: E(N)    ! ENERGY EIGENVALUE
      REAL(8)   ,INTENT(OUT):: F(N)    ! NUMBER OF ELECTRONS PER SHELL
      REAL(8)   ,INTENT(OUT):: MU      ! CHEMICAL POTENTIAL
      REAL(8)               :: MUMIN,MUMAX
      REAL(8)               :: QMAX
      REAL(8)               :: DQ
      REAL(8)               :: X
      INTEGER(4),PARAMETER  :: NITER=1000
      REAL(8)   ,PARAMETER  :: QTOL=1.D-10
      INTEGER(4)            :: ITER,I
!     **************************************************************************
      MUMIN=MINVAL(E)-10.D0*KBT
      MUMAX=MAXVAL(E)+10.D0*KBT
      QMAX=2.D0*REAL(SUM(2*L(:)+1),KIND=8)
      IF(Q.LT.0.D0.OR.Q.GT.QMAX) THEN
        CALL ERROR$MSG('CHARGE OUT OF RANGE')
        CALL ERROR$R8VAL('Q',Q)
        CALL ERROR$R8VAL('QMIN',0.D0)
        CALL ERROR$R8VAL('QMAX',QMAX)
        CALL ERROR$STOP('ONEATOM_OPTF')
      END IF
      DO ITER=1,NITER
        MU=0.5D0*(MUMIN+MUMAX)
        DO I=1,N
          X=(E(I)-MU)/KBT
          IF(X.GT.50.D0) THEN
            F(I)=0.D0
          ELSE IF(X.LT.-50.D0) THEN
            F(I)=REAL(2*(2*L(I)+1),KIND=8)
          ELSE
            F(I)=REAL(2*(2*L(I)+1),KIND=8)/(1.D0+EXP(X))
          END IF
        ENDDO
        DQ=SUM(F)-Q
        IF(ABS(DQ).LT.QTOL) EXIT
        IF(DQ.GT.0.D0) THEN
          MUMAX=MU
        ELSE
          MUMIN=MU
        END IF
      ENDDO
      IF(ABS(DQ).GT.QTOL) THEN
        CALL ERROR$MSG('LOOP NOT CONVERGED')
        CALL ERROR$MSG('ONEATOM_OPTF')
      END IF
      RETURN
      END
!
!     .....................................................VOUT.........
      SUBROUTINE MYVOFRHOWITHRAD(GID,NR,RAD,AEZ,RHO,POTH,POTXC,EH,EXC)
!     ******************************************************************
!     **                                                              **
!     **  ELECTROSTATIC AND EXCHANGE-CORRELATION POTENTIAL            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: AEZ
      REAL(8),   INTENT(IN) :: RAD
      REAL(8),   INTENT(IN) :: RHO(NR)
      REAL(8),   INTENT(OUT):: POTh(NR)
      REAL(8),   INTENT(OUT):: POTxc(NR)
      REAL(8),   INTENT(OUT):: EH
      REAL(8),   INTENT(OUT):: EXC
      REAL(8)               :: ALPHA
      REAL(8),   PARAMETER  :: RHOMIN=1.D-2
      REAL(8)               :: PI
      REAL(8)               :: FOURPI
      REAL(8)               :: Y0
      REAL(8)               :: RHO1(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: RHOPLUS(NR)
      REAL(8)               :: DRHOPLUSDRHO(NR)
      REAL(8)               :: EDEN(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      REAL(8)               :: EVAL,MUVAL
      REAL(8)               :: SVAR
      REAL(8)               :: F,X,DF,YP,YM,NBYF,DFDN
      INTEGER(4)            :: IR,IRBOX
      REAL(8)               :: Q
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).Ge.RAD) EXIT
      ENDDO
!
!     ==================================================================
!     ==  TOTAL CHARGE                                                ==
!     ==================================================================
      AUX(:)=4.D0*PI*RHO(:)*R(:)**2*Y0
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,Q)
      Q=Q-AEZ
!
!     ==================================================================
!     ==  TOTAL POTENTIAL                                             ==
!     ==================================================================
      RHO1=RHO
      IR=MIN(IRBOX+3,NR)
      RHO1(IR:)=0.D0
      CALL RADIAL$NUCPOT(GID,NR,AEZ,POTH)
      EDEN(:)=0.5D0*RHO1(:)*POTH(:)
      CALL RADIAL$POISSON(GID,NR,0,RHO1,AUX)
      POTH(:)=POTH(:)+AUX(:)
      CALL RADIAL$VALUE(GID,NR,POTH,RAD,SVAR)
      SVAR=Q/RAD/Y0-SVAR
      POTH(1:IRBOX-1)=POTH(1:IRBOX-1)+SVAR
      POTH(IRBOX:NR)=Q/R(IRBOX:NR)/Y0
      EDEN(:)=EDEN(:)+0.5D0*RHO1(:)*POTH(:)
!
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,EDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EH)
!
!     ==================================================================
!     ==  EXCHANGE CORRELATION
!     ==================================================================
      RHO1(:IRBOX)=RHO(:)
      RHO1(IRBOX:)=0.D0     
      ALPHA=LOG(2.D0)/RHOMIN
      DO IR=1,NR
        X=RHO1(IR)*Y0
        IF(R(IR).GT.RAD) X=0.D0
        IF(ALPHA*X.GT.30.D0) THEN
          F=X
          DF=1.D0
        ELSE
          YP=EXP(ALPHA*X)
          YM=1.D0/YP
          F=LOG(YP+YM)/ALPHA
          DF=(YP-YM)/(YP+YM)
        END IF
        RHOPLUS(IR)=F/Y0
        DRHOPLUSDRHO(IR)=DF/Y0
      ENDDO
!
      CALL RADIAL$DERIVE(GID,NR,RHOPLUS(:),GRHO)
      DO IR=1,NR
        RH=RHOPLUS(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC1,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        EDEN(IR)=4.D0*PI*EXC1   ! ANGULAR INTEGRATION ALREADY INCLUDED
        POTXC(IR)=VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      CALL RADIAL$DERIVE(GID,NR,GRHO(:),AUX)
      POTXC(:)=POTXC(:)-AUX
      IF(R(1).GT.1.D+10) THEN
         POTXC(:)=POTXC(:)-2.D0/R(:)*GRHO(:)
      ELSE
        POTXC(2:)=POTXC(2:)-2.D0/R(2:)*GRHO(2:)
        POTXC(1)=POTXC(1)-2.D0/R(2)*GRHO(2)
        POTXC(1)=POTXC(2)
      END IF
!
      DO IR=1,NR
        EVAL=EDEN(IR)/(4.D0*PI)
        MUVAL=POTXC(IR)*Y0
!       ==  N/F(N) =============================
        NBYF=RHO1(IR)/RHOPLUS(IR)
        DFDN=DRHOPLUSDRHO(IR)*Y0
        SVAR=(1.D0-NBYF*DFDN)/RHOPLUS(IR)
        EDEN(IR)=NBYF*EVAL
        POTXC(IR)=NBYF*MUVAL+SVAR*EVAL
        EDEN(IR)=EDEN(IR)*4.D0*PI
        POTXC(IR)=POTXC(IR)/Y0
      ENDDO
!
!     == exchange correlation potential is set to zero outside the box. 
!     == otherwise it will have a constant value in the zero-density region
      potxc(irbox:)=0.d0
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,EDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EXC)
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE STRUCTURE_MODULE
LOGICAL(4),SAVE      :: TINI=.FALSE.
INTEGER(4),PARAMETER :: NAT=2
REAL(8)              :: RBAS(3,3)
REAL(8)              :: RPOS(3,NAT)
REAL(8)              :: BAREVOL
REAL(8)              :: MADMAT(NAT,NAT)
END MODULE STRUCTURE_MODULE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ESTRUCTURE_INITIALIZE()
!     **************************************************************************
!     **************************************************************************
      USE STRUCTURE_MODULE
      IMPLICIT NONE
      REAL(8)         ::GBAS(3,3)
      INTEGER(4)      :: IAT
      REAL(8)         :: Q(NAT)
      REAL(8)         :: ETOT
      REAL(8)         :: POT(NAT)
      REAL(8)         :: FORCE(3,NAT)
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
!     == FCC ===================================================================
      RBAS(:,1)=(/0.0D0,0.5D0,0.5D0/)
      RBAS(:,2)=(/0.5D0,0.0D0,0.5D0/)
      RBAS(:,3)=(/0.5D0,0.5D0,0.0D0/)
      RPOS(:,1)=(/0.0D0,0.0D0,0.0D0/)
      RPOS(:,2)=(/0.5D0,0.5D0,0.5D0/)
!
!     == UNIT CELL VOLUME ======================================================
      CALL GBASS(RBAS,GBAS,BAREVOL)
!     == DETERMINE MADELUNG MATRIX =============================================
      DO IAT=1,NAT
        Q(:)=0.D0
        Q(IAT)=1.D0
        CALL MADELUNG(NAT,RBAS,RPOS,Q,ETOT,POT,FORCE)
        MADMAT(:,IAT)=POT(:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ESTRUCTURE(NAT_,Q,VOL,ETOT,POT,PRESSURE)
!     **************************************************************************
!     **************************************************************************
      USE STRUCTURE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT_
      REAL(8)   ,INTENT(IN) :: Q(NAT_)
      REAL(8)   ,INTENT(IN) :: VOL
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,INTENT(OUT):: POT(NAT_)
      REAL(8)   ,INTENT(OUT):: PRESSURE
      REAL(8)               :: FORCE(3,NAT)
      REAL(8)               :: SCALE
!     **************************************************************************
      CALL ESTRUCTURE_INITIALIZE()
      IF(NAT_.NE.NAT) THEN
        CALL ERROR$MSG('NAT INCONSISTENT')
        CALL ERROR$STOP('ESTRUCTURE')
      END IF
      SCALE=(VOL/BAREVOL)**(1.D0/3.D0)
      CALL MADELUNG(NAT,RBAS*SCALE,RPOS*SCALE,Q,ETOT,POT,FORCE)
! POT(:)=MATMUL(MADMAT,Q)/SCALE
! ETOT(:)=0.5D0*DOT_PRODUCT(POT,Q)
      PRESSURE=ETOT/(3.D0*VOL)
      RETURN
      END

!
!     ......................................................................
      SUBROUTINE AESCF(GID,NR,NB,TREL,LOFI,SO,F,NN,Z,RHOADD,RBOX,DREL,POT,EOFI)
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NB        ! #(STATES)
      LOGICAL(4) ,INTENT(IN)     :: TREL
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)  !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO(NB)    !SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: F(NB)     !OCCUPATION
      INTEGER(4) ,INTENT(IN)     :: NN(NB)    !#(NODES)
      REAL(8)    ,INTENT(IN)     :: Z         !ATOMIC NUMBER
      REAL(8)    ,INTENT(IN)     :: RHOADD(NR)! FIXED CONTRIBUTION OF DENSITY
      REAL(8)    ,INTENT(IN)     :: RBOX      !ATOM ENCLOSED IN A BOX 
      REAL(8)    ,INTENT(OUT)    :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(INOUT)  :: POT(NR)   !POTENTIAL
      REAL(8)    ,INTENT(OUT)    :: EOFI(NB)  !ONE-PARTICLE ENERGIES
      REAL(8)                    :: R(NR)
      REAL(8)                    :: RHO(NR)
      REAL(8)                    :: EREF
      INTEGER(4)                 :: ITER
      INTEGER(4)                 :: NITER=50
      REAL(8)                    :: XAV,XMAX
      LOGICAL(4)                 :: CONVG
      REAL(8)   ,PARAMETER       :: TOL=1.D-5
      INTEGER(4)                  :: NFILO
      REAL(8)                    :: EH,EXC
      LOGICAL(4),PARAMETER       :: TBROYDEN=.TRUE.
      REAL(8)                    :: POTIN(NR)
      REAL(8)                    :: svar
      INTEGER(4)                 :: IRBOX
      INTEGER(4)                 :: IR
      CHARACTER(32)              :: ID
REAL(8)                    :: POTS(NR,100)
!     ***********************************************************************
      CALL RADIAL$R(GID,NR,R)
      IRBOX=0
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GT.RBOX) EXIT
      ENDDO
      EOFI=0.D0
      EREF=0.D0
      XAV=0.D0
      XMAX=0.D0
      CONVG=.FALSE.
      IF(TBROYDEN) THEN
        CALL BROYDEN$NEW(NR,4,1.D-1)
      ELSE
        CALL MYMIXPOT('ON',GID,NR,POT,XMAX,XAV)
      END IF
      DO ITER=1,NITER
        IF(TREL) THEN
          CALL SCHROEDINGER$DREL(GID,NR,POT,EREF,DREL)
          ID='FULL'
        ELSE
          DREL=0.D0
          ID='NONREL'
        END IF
POTS(:,ITER)=POT
        CALL AERHO(ID,GID,NR,NB,LOFI,SO,F,NN,RBOX,DREL,POT,RHO,EOFI)
        RHO(:)=RHO(:)+RHOADD(:)
        EREF=EOFI(NB)
!
!       ================================================================
!       ==  EXIT IF CONVERGED                                         ==
!       ================================================================
        IF(CONVG) THEN
          CALL FILEHANDLER$UNIT('PROT',NFILO)
          CALL REPORT$I4VAL(NFILO,"SCF LOOP CONVERGED AFTER ",ITER," ITERATIONS")
          CALL REPORT$R8VAL(NFILO,'AV. DIFF BETWEEN IN- AND OUTPUT POTENTIAL',XAV,'H')
          CALL REPORT$R8VAL(NFILO,'MAXIMUM  OF (VIN-VOUT)*R^2',XMAX,'H')
          IF(TBROYDEN) THEN
            CALL BROYDEN$CLEAR
          ELSE
            CALL MYMIXPOT('OFF',GID,NR,POT,XMAX,XAV)
          END IF
          RETURN
        END IF

!       ====================================================================
!       == CALCULATE OUTPUT POTENTIAL                                     ==
!       ====================================================================
        IF(TBROYDEN) POTIN=POT
        CALL MYVOFRHO(GID,NR,Z,RHO,POT,EH,EXC)
!
!       ================================================================
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSON'S METHOD     ==
!       ================================================================
        IF(TBROYDEN) THEN
          XAV=SQRT(DOT_PRODUCT(POT-POTIN,POT-POTIN)/REAL(NR,KIND=8))
          XMAX=MAXVAL(ABS(POT-POTIN)) 
          CALL BROYDEN$STEP(NR,POTIN,(POT-POTIN))
          POT=POTIN
       ELSE
          CALL MYMIXPOT('GO',GID,NR,POT,XMAX,XAV)
        END IF
PRINT*,'XMAX ',XMAX,XAV
        CONVG=(XMAX.LT.TOL)
      ENDDO
      CALL WRITEPHI('POTS.DAT',GID,NR,10,POTS(:,1:10))
      CALL ERROR$MSG('SELFCONSISTENCY LOOP NOT CONVERGED')
      CALL ERROR$STOP('AESCF')
      RETURN
      END
!
!     .....................................................VOUT.........
      SUBROUTINE MYVOFRHO(GID,NR,AEZ,RHO,POT,EH,EXC)
!     ******************************************************************
!     **                                                              **
!     **  ELECTROSTATIC AND EXCHANGE-CORRELATION POTENTIAL            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: AEZ
      REAL(8),   INTENT(IN) :: RHO(NR)
      REAL(8),   INTENT(OUT):: POT(NR)
      REAL(8),   INTENT(OUT):: EH
      REAL(8),   INTENT(OUT):: EXC
      REAL(8)               :: PI
      REAL(8)               :: FOURPI
      REAL(8)               :: Y0
      REAL(8)               :: AUX(NR)
      REAL(8)               :: EDEN(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      INTEGER(4)            :: IR
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  TOTAL POTENTIAL                                             ==
!     ==================================================================
      EDEN=0.D0
      CALL RADIAL$POISSON(GID,NR,0,RHO,POT)
      EDEN(:)=EDEN(:)+0.5D0*RHO(:)*POT(:)
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AUX)
      POT(:)=POT(:)+AUX(:)
      EDEN(:)=EDEN(:)+RHO(:)*AUX(:)
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,EDEN,EH)
!
!     ==================================================================
!     ==  EXCHANGE CORRELATION
!     ==================================================================
      CALL RADIAL$DERIVE(GID,NR,RHO(:),GRHO)
      DO IR=1,NR
        RH=RHO(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC1,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        EDEN(IR)=4.D0*PI*EXC1   ! ANGULAR INTEGRATION ALREADY INCLUDED
        POT(IR)=POT(IR)+VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      CALL RADIAL$DERIVE(GID,NR,GRHO(:),AUX)
      POT(:)=POT(:)-AUX
      IF(R(1).GT.1.D+10) THEN
         POT(:)=POT(:)-2.D0/R(:)*GRHO(:)
      ELSE
        POT(2:)=POT(2:)-2.D0/R(2:)*GRHO(2:)
        POT(1)=POT(1)-2.D0/R(2)*GRHO(2)
        POT(1)=POT(2)
      END IF
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,EDEN,EXC)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MYMIXPOT(SWITCH,GID,NR,POT,XMAX,XAV)
!     ******************************************************************
!     **                                                              **
!     **  MIX POTENTIAL USING D.G.ANDERSEN'S METHOD                   **
!     **                                                              **
!     **  1) INITIALIZE WITH SWITCH='ON'                              **
!     **     STORES THE FIRST INPUT POTENTIAL <-POT                   **
!     **     ALLOCATES  INTERNAL ARRAYS                               **
!     **                                                              **
!     **  2) ITERATE WITH WITH SWITCH='GO'                            **
!     **     RECEIVES THE OUTPUT POTENTIAL    <-POT           E       **
!     **     CALCULATES NEW INPUT POTENTIAL   ->POT                   **
!     **     STORES THE NEW INPUT POTENTIAL                           **
!     **                                                              **
!     **  3) CLEAR MEMORY WITH SWITCH='OFF'                           **
!     **                                                              **
!     **  WARNING! DO NOT USE SIMULTANEOUSLY FOR TOW DIFFERENT SCHEMES**
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: SWITCH !CAN BE 'ON', 'GO' OR 'OFF'
      INTEGER(4)  ,INTENT(IN)   :: GID
      INTEGER(4)  ,INTENT(IN)   :: NR
      REAL(8)     ,INTENT(INOUT):: POT(NR) ! 
      REAL(8)     ,INTENT(OUT)  :: XMAX    ! MAX.|R**2*(VOUT-VIN)|
      REAL(8)     ,INTENT(OUT)  :: XAV     ! <R**2(VOUT-VIN)**2>
      REAL(8)     ,PARAMETER    :: ALPHA=1.D-1 !MIXING PARAMETER
      CHARACTER(8)       ,SAVE :: STATUS='OFF'
      LOGICAL(4)         ,SAVE :: TSTART
      INTEGER(4)         ,SAVE :: NRSAVE
      REAL(8),ALLOCATABLE,SAVE :: OLDPOTIN(:)
      REAL(8),ALLOCATABLE,SAVE :: OLDPOTOUT(:)
      REAL(8),ALLOCATABLE,SAVE :: NEWPOTOUT(:)
      REAL(8),ALLOCATABLE,SAVE :: NEWPOTIN(:)
      REAL(8)                  :: BETA
      REAL(8)                  :: SVAR1,SVAR2
      REAL(8)                  :: R(NR)
      REAL(8)                  :: AUX(NR),AUX1(NR),AUX2(NR)
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
      IF(SWITCH.EQ.'GO') THEN
        IF(NR.NE.NRSAVE) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF GRIDPOINTS')
          CALL ERROR$STOP('MIXPOT')
        END IF
!       ================================================================
!       == COPY POT INTO NEWPOTIN                                     ==
!       ================================================================
        NEWPOTOUT(:)=POT(:)
!       ================================================================
!       == CALCULATE MAX AND VARIANCE OF POTOUT-POTIN                 ==
!       ================================================================
        AUX(:)=NEWPOTOUT(:)-NEWPOTIN(:)
        XMAX=MAXVAL(ABS(AUX(:)))
        XAV=SUM(AUX(:)**2)
        SVAR1=REAL(NR,KIND=8)
        XAV=SQRT(XAV/SVAR1)
!       ================================================================
!       ==  CALCULATE MIXING FACTOR BETA                              ==
!       ================================================================
        IF(TSTART) THEN
          TSTART=.FALSE.
          BETA=0.D0
        ELSE
          AUX1(:)=NEWPOTOUT(:)-NEWPOTIN(:)
          AUX2(:)=OLDPOTOUT(:)-OLDPOTIN(:)
          AUX(:)=AUX1(:)-AUX2(:)
          SVAR1=SUM(AUX1(:)*AUX(:))
          SVAR2=SUM(AUX(:)**2)
          BETA=SVAR1/SVAR2
        END IF
!       ================================================================
!       == MIX POTENTIALS                                             ==
!       ================================================================
        AUX1(:)=(1.D0-BETA)*NEWPOTIN(:) + BETA*OLDPOTIN(:)
        AUX2(:)=(1.D0-BETA)*NEWPOTOUT(:)+ BETA*OLDPOTOUT(:)
        POT(:)=AUX1(:) + ALPHA*(AUX2(:)-AUX1(:))
        OLDPOTIN(:) =NEWPOTIN(:)
        OLDPOTOUT(:)=NEWPOTOUT(:)
        NEWPOTIN(:) =POT(:)
!
!     ==================================================================
!     == INITIALIZE MIXING                                            ==
!     ==================================================================
      ELSE IF(SWITCH.EQ.'ON') THEN
        IF(TRIM(STATUS).NE.'OFF') THEN
          DEALLOCATE(OLDPOTIN)
          DEALLOCATE(OLDPOTOUT)
          DEALLOCATE(NEWPOTIN)
          DEALLOCATE(NEWPOTOUT)
        END IF
        NRSAVE=NR
        ALLOCATE(OLDPOTIN(NRSAVE))
        ALLOCATE(OLDPOTOUT(NRSAVE))
        ALLOCATE(NEWPOTIN(NRSAVE))
        ALLOCATE(NEWPOTOUT(NRSAVE))
        TSTART=.TRUE.
        STATUS='ON'
        OLDPOTOUT(:)=0.D0
        OLDPOTIN(:) =0.D0
        NEWPOTOUT(:)=0.D0
        NEWPOTIN(:) =POT(:)
!
!     ==================================================================
!     == CLEAR ARRAYS                                                 ==
!     ==================================================================
      ELSE IF(SWITCH.EQ.'OFF') THEN
        DEALLOCATE(OLDPOTIN)
        DEALLOCATE(OLDPOTOUT)
        DEALLOCATE(NEWPOTIN)
        DEALLOCATE(NEWPOTOUT)
        TSTART=.FALSE.
        STATUS='OFF'
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AERHO(ID,GID,NR,NB,LOFI,SO,F,NN,RBOX,DREL,POT,RHO,E)
!     **                                                                      **
!     **  DETERMINES THE DENSITY FOR A GIVEN POTENTIAL AND A SPECIFIED SET    **
!     **  OF WAVE FUNCTIONS.                                                  **
!     **  UPDATES THE ONE-PARTICLE ENERGIES AND RELATIVISTIC CORRECTION "DREL"**
!     **                                                                      **
!     **  THE ORBITALS ARE SPECIFIED BY ANGULAR MOMENTUM "LOFI",              **
!     **                             THE NUMBER OF NODES "NN",                **
!     **                             AND THE SPIN ORBIT PARAMETER "SO"        **
!     **                                                                      **
!     **  VALUES OF ID:                                                       **
!     **    NONREL: NONRELATIVISTIC CALCULATION, DREL=0 IS RETURNED           **
!     **    EFFZORA: THE RELATIVISTIC CORRECTIONS ARE TAKEN FROM INPUT        **
!     **            AND REMAIN UNCHANED                                       **
!     **    FULL: PERFORMS RELATIVISTIC CALCULATION USING DREL FROM INPUT     **
!     **          AND RECALCULATES DREL FROM THE MEAN KINETIC ENERGY DENSITY  **
!     **          DIVIDED BY THE DENSITY.                                     **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4) ,INTENT(IN)     :: GID       ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR        ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: NB        ! #(STATES)
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)  !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO(NB)    !SWITCH FOR SPIN-ORBIT COUP.
      INTEGER(4) ,INTENT(IN)     :: NN(NB)    !#(NODES)
      REAL(8)    ,INTENT(IN)     :: F(NB)     !OCCUPATION
      REAL(8)    ,INTENT(IN)     :: RBOX      !ATOM ENCLOSED IN A BOX
      REAL(8)    ,INTENT(INOUT)  :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: POT(NR)   !POTENTIAL
      REAL(8)    ,INTENT(OUT)    :: RHO(NR)   !DENSITY
      REAL(8)    ,INTENT(INOUT)  :: E(NB)     !ONE-PARTICLE ENERGIES
      REAL(8)                    :: R(NR)     !RADIAL GRID POINTS
      REAL(8)                    :: G(NR)     !INHOMOGENEITY (NOT USED)
      REAL(8)                    :: AUX(NR),AUX1(NR)
      REAL(8)                    :: PHI(NR)
      REAL(8)                    :: EKIN(NR)
      REAL(8)                    :: EKIN2(NR)
      REAL(8)                    :: ARRAY(NR,5)
      INTEGER(4)                 :: IB,IR
      REAL(8)                    :: SVAR
      REAL(8)                    :: C0LL,PI,Y0
!     **************************************************************************
      IF(ID.NE.'FULL'.AND.ID.NE.'EFFZORA'.AND.ID.NE.'NONREL') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$MSG('ALLOWED VALUES ARE "FULL",EFFZORA","NONREL"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AESCF')
      ENDIF
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == DETERMINE BOUND STATES FOR A GIVEN POTENTIAL AND ADD TO DENSITY      ==
!     ==========================================================================
      RHO(:)=0.D0
      EKIN(:)=0.D0
      EKIN2(:)=0.D0
      DO IB=1,NB
         G(:)=0.D0
         IF(ID.EQ.'FULL') THEN
           CALL SCHROEDINGER$DREL(GID,NR,POT,E(IB),DREL)
         ELSE IF(ID.EQ.'NONREL') THEN
           DREL(:)=0.D0 
         END IF
         CALL BOUNDSTATE(GID,NR,LOFI(IB),SO(IB),RBOX,DREL,G,NN(IB),POT,E(IB),PHI)
         AUX(:)=(R(:)*PHI(:))**2
         CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
         CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
         PHI(:)=PHI(:)/SQRT(SVAR)
         RHO(:)  =RHO(:)  +F(IB)*C0LL         *PHI(:)**2
         EKIN(:) =EKIN(:) +F(IB)*C0LL*E(IB)   *PHI(:)**2
         EKIN2(:)=EKIN2(:)+F(IB)*C0LL*E(IB)**2*PHI(:)**2
      ENDDO
      DO IR=1,NR
        IF(R(IR).GT.RBOX) RHO(IR)=0.D0
      ENDDO
!
!     ==========================================================================
!     == DETERMINE RELATIVISTIC CORRECTION                                    ==
!     ==========================================================================
      IF(ID.EQ.'FULL') THEN
        AUX(:)=MAX(RHO(:),1.D-6)
        AUX=1.D0/AUX
        EKIN(:)=EKIN(:)*AUX(:)
        EKIN2(:)=EKIN2(:)*AUX(:)
        EKIN2(:)=MAX(EKIN2(:)-EKIN(:)**2,0.D0)
        EKIN2(:)=SQRT(EKIN2(:))
        DO IR=1,NR
          IF(RHO(IR).LT.1.D-5) THEN
            EKIN(IR)=0.D0
            EKIN2(IR)=0.D0
          END IF
        ENDDO
        EKIN(:)=EKIN(:)/Y0-POT(:)
        CALL SCHROEDINGER$DREL(GID,NR,-EKIN,0.D0,DREL)
      END IF
      IF(ID.EQ.'FULL') THEN
        ARRAY(:,1)=POT(:)*Y0
        ARRAY(:,2)=(EKIN+POT(:))*Y0
        ARRAY(:,3)=(EKIN+POT(:))*Y0+EKIN2(:)
        ARRAY(:,4)=(EKIN+POT(:))*Y0-EKIN2(:)
        CALL WRITEPHI('EKIN.DAT',GID,NR,4,ARRAY)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BOUNDSTATE(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,PHI)
!     **************************************************************************
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND         **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                     **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G                **
!     **                                                                      **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE                 **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.         **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.         **
!     **                                                                      **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND           **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT           **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: RBOX    ! BOX RADIUS
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-12
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: DER,DERO
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      REAL(8)                    :: PHI1(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20 ! MAXIMUM FACTOR IN THE WAVE FUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      REAL(8)                    :: ROUT
      REAL(8)                    :: VAL,VAL1,VAL2,R1,R2
      INTEGER(4)                 :: IROUT,IRCL,IRBOX,IREND
!     *********************************************************************
                                 CALL TRACE$PUSH('BOUNDSTATE')
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST IOUTSIDE THE BOX
      IRBOX=1
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).Ge.RBOX) EXIT
      ENDDO
!          
      ISTART=1
      X0=E
      DX=1.D-2
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO I=1,NITER
        E=X0
!        E=MIN(X0,EMAX)
!
!       =======================================================================
!       ==  CUT OFF THE POTENTIAL                                            ==
!       =======================================================================
        CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)
!       == BOUNDARY CONDITION PHI(ROUT)=0 =====================================
        IF(R(IROUT).LT.RBOX) THEN
          ROUT=R(IROUT)-1.D-5    !ENSURE THAT ROUT<R(IROUT)
        ELSE
          ROUT=RBOX
          IROUT=IRBOX
        END IF
!       ==  SET KINETIC ENERGY TO ZERO BEYOND ROUT
        POT1(:)=POT(:)
        POT1(IROUT:)=POT(IROUT)
!
!       =======================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       =======================================================================
        IDIR=1
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(IROUT).GT.0.OR.PHI(IROUT).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NN',NN)
          CALL ERROR$STOP('BOUNDSTATE')
        END IF
!
!       =======================================================================
!       == ESTIMATE PHASE SHIFT                                              ==
!       =======================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,ROUT,Z0)
        Z0=Z0-REAL(NN+1)
        IF(ABS(2.D0*DX).LE.TOL) EXIT
!       =====================================================================
!       ==  BISECTION                                                      ==
!       =====================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('BOUNDSTATE')
      END IF
!
!     =======================================================================
!     ==  DETERMINE MATCHING POINT                                         ==
!     =======================================================================
      IRMATCH=IRCL
      IF(R(IRCL).GT.5.D0) THEN
        CALL RADIAL$XOFR(GID,5.D0,SVAR)
        IRMATCH=INT(SVAR)
      END IF
!
!     =======================================================================
!     ==  INTEGRATE INWARD                                                 ==
!     =======================================================================
      IF(IRMATCH.LT.IROUT) THEN
        THOM=MAXVAL(ABS(G(:))).EQ.0.D0
        IDIR=-1
!       ==  HOMOGENEOUS SOLUTION THAT FULFILLS THE OUTER BOUNDARY CONDITION
        IREND=MIN(IROUT,NR-3)
        GHOM(:)=0.D0
        GHOM(IREND+1)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
        GHOM(:)=0.D0
        GHOM(IREND+2)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHI1)
        PHI1(:)=PHI1(:)/PHI1(IRMATCH)
!       == EXTRAPOLATE =======================
        IF(IREND.LT.IROUT) THEN
          R1=R(IREND)
          R2=R(IREND+1)
          VAL1=PHIHOM(IREND)
          VAL2=PHIHOM(IREND+1)
          PHIHOM(IREND:)=VAL1+(VAL2-VAL1)/(R2-R1)*(R(IREND:)-R1)
          VAL1=PHI1(IREND)
          VAL2=PHI1(IREND+1)
          PHI1(IREND:)=VAL1+(VAL2-VAL1)/(R2-R1)*(R(IREND:)-R1)
        END IF
!       == FULFILL OUTER BOUNDARY CONDITION =====================================
        CALL RADIAL$VALUE(GID,NR,PHIHOM,ROUT,VAL1)
        CALL RADIAL$VALUE(GID,NR,PHI1,ROUT,VAL2)
        SVAR=VAL1+VAL2
        VAL1=VAL1/SVAR
        VAL2=VAL2/SVAR
        PHIHOM(:)=VAL2*PHIHOM(:)-VAL1*PHI1(:)
        PHIHOM(:)=PHIHOM(:)/PHIHOM(IRMATCH)
!       == INHOMOGENEOUS SOLUTION WITH CORRECT BOUNDARY CONDITIONS
        IF(.NOT.THOM) THEN     
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
          CALL RADIAL$VALUE(GID,NR,PHIINHOM,ROUT,VAL1)
          CALL RADIAL$VALUE(GID,NR,PHI1,ROUT,VAL2)
          PHIINHOM(:)=PHIINHOM(:)-VAL1/VAL2*PHI1(:)
         ELSE
          PHIINHOM(:)=0.D0
        END IF
!
!       =======================================================================
!       ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!       =======================================================================
        SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
        PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
        CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
        CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
        SVAR=(DERO-DER)/PHI(IRMATCH)
        PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
      END IF
!
!     =======================================================================
!     ==  SET WAVE FUNCTION TO ZERO BEYOND RBOX                            ==
!     =======================================================================
DO IR=1,NR
  IF(.NOT.(PHI(IR).GT.0.D0.OR.PHI(IR).LE.0.D0)) THEN
    PRINT*,'ERROR'
!    PRINT*,'PHIIN',PHI(:IRMATCH-1)
!    PRINT*,'PHIOUT',PHI(IRMATCH:)
    PRINT*,'SVAR ',SVAR
    PRINT*,'IROUT,IRCL ',IROUT,IRCL,NR
    PRINT*,'R,PHI ',R(IR),PHI(IR)
!!$    OPEN(UNIT=8,FILE='XXX.DAT')
!!$    DO I=1,NR
!!$      WRITE(8,*)I,R(I),PHIHOM(I),POT1(I)*Y0
!!$    ENDDO
!!$    CLOSE(8)      
    CALL ERROR$MSG('PHI CONTAINS NANS')
    CALL ERROR$STOP('BOUNDSTATE')
  END IF
ENDDO
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEPHI(FILE,GID,NR,NPHI,PHI)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI(NR,NPHI)
      INTEGER(4)              :: IR
      REAL(8)                 :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        WRITE(100,FMT='(F15.10,2X,100(F25.10,2X))')R(IR),PHI(IR,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE spherical$shiftcenter(GID,NR,centernew,lmx1,fin,lmx2,fout)
!     **                                                                      **
!     **  transforms a function given in true spherical harmonics             **
!     **  (angular momentum eigenstates) to a new center                      **
!     **  ths center is displaced in z-direction                              **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: lmx1
      INTEGER(4)  ,INTENT(IN) :: lmx2
      real(8)     ,intent(in) :: centernew(3)
      real(8)     ,intent(in) :: fin(nr,lmx1)
      real(8)     ,intent(out):: fout(nr,lmx2)
      real(8)                 :: finr(nr,lmx1)
      real(8)                 :: foutr(nr,lmx2)
      real(8)                 :: fint(nr,lmx2)
      real(8)                 :: dis
      real(8)                 :: r(nr)
      real(8)                 :: rot(3,3)
      real(8)                 :: ylmrot1(lmx1,lmx1),ylmrot2(lmx2,lmx2)
      integer(4)              :: lx1,lx2
      real(8)                 :: plm1(lmx1),plm2(lmx2)
      real(8)                 :: fone(lmx1)
      logical(4)              :: tchk
      real(8)                 :: cosalpha,costheta
      real(8)                 :: sinalpha,sintheta
      real(8)                 :: dxdalpha,dxdtheta
      real(8)                 :: aux(nr)
      integer(4)              :: ir1,ir2
      real(8)                 :: p,x,dx             
      integer(4)              :: lm1,lm2,l1,l2,m1,m2,l1start,itheta
      real(8)                 :: fac,svar
      integer(4) ,parameter   :: ntheta=300
      real(8)                 :: theta,dtheta
      real(8)                 :: pi
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      lx1=sqrt(real(lmx1-1,kind=8)+1.d-10)
      lx2=sqrt(real(lmx2-1,kind=8)+1.d-10)
      call radial$r(gid,nr,r)
      dis=sqrt(sum(centernew(:)**2))
!
!     ==========================================================================
!     == Rotate coordinate system such that the new center is in positive z dir=
!     ==========================================================================
      call SPHERICAL$align(centernew,rot)
      call ROTATEYLM(lmx1,rot,ylmrot1(:,:))

!     ==  finr is the rotated  function
      finr(:,:)=0.d0
      do lm1=1,lmx1    !loop over the components of the rotated vector
        do lm2=1,lmx1 !the ilm component of the lmx vector for a given r after rotation
          svar=ylmrot1(lm1,lm2)
          if(svar.eq.0.d0) cycle
          finr(:,lm1)=finr(:,lm1)+fin(:,lm2)*svar
        enddo
      enddo
!
!     ==========================================================================
!     ==  transform to new center                                             ==
!     ==========================================================================
      dtheta=pi/real(ntheta,kind=8)
      fout(:,:)=0.d0
      DO IR2=1,NR
        P=R(IR2)
        DO itheta=1,ntheta
          theta=dtheta*(-0.5d0+real(itheta,kind=8))
          costheta=cos(theta)
          x=sqrt(dis**2+p**2+2.d0*dis*p*costheta)
          COSALPHA=(dis+P*costheta)/x
          SINTHETA=SQRT(1.D0-COSTHETA**2)
!         == START RECURSION FOR ASSOCIATED LEGENDRE POLYNOMIALS WITH P_LL
          CALL PLGNDR(LMX1,LX1,COSALPHA,PLM1)
          CALL PLGNDR(LMX2,LX2,COSTHETA,PLM2)
          do lm1=1,lmx1
            call radial$value(gid,nr,finr(:,lm1),x,fone(lm1))
          enddo
          DO LM2=1,LMX2
            L2=SQRT(REAL(LM2-1,KIND=8)+1.D-10)
            m2=L2*(l2+1)-LM2+1
            M1=M2
            L1START=ABS(M2)
            svar=0.d0
            DO L1=L1START,LX1
              LM1=L1*(l1+1)-M2+1
              FAC=0.5D0*SQRT(REAL((2*L1+1)*(2*L2+1),KIND=8))
              svar=svar+FAC*fone(LM1)*PLM1(LM1)*PLM2(LM2)
            ENDDO
            fout(ir2,lm2)=fout(ir2,lm2)+svar*sintheta*dtheta
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  rotate back to original coordinate system                           ==
!     ==========================================================================
      rot=transpose(rot)
      call ROTATEYLM(lmx2,rot,ylmrot2(:,:))
!     ==  finr is the rotated  function
      fint(:,:)=fout(:,:)
      fout(:,:)=0.d0
      do lm1=1,lmx2    !loop over the components of the rotated vector
        do lm2=1,lmx2 !the ilm component of the lmx vector for a given r after rotation
          svar=ylmrot2(lm1,lm2)
          if(svar.eq.0.d0) cycle
          fout(:,lm1)=fout(:,lm1)+fint(:,lm2)*svar
        enddo
      enddo
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine SPHERICAL$align(vec,rot)
!     **************************************************************************    
!     ** construct a rotation matrix rot which transforms the vector vec      **
!     ** onto the positive z-axis                                             **
!     **************************************************************************    
      implicit none
      real(8),intent(in)         :: vec(3)
      real(8),intent(out)        :: rot(3,3)
      integer(4)                 :: i,ivec(1)
      real(8)                    :: dis,dr1(3),DR2(3),DR3(3)
      real(8),parameter          :: ex(3)=(/1.d0,0.d0,0.d0/)
      real(8),parameter          :: ey(3)=(/0.d0,1.d0,0.d0/)
      real(8),parameter          :: ez(3)=(/0.d0,0.d0,1.d0/)
!     **************************************************************************    
      dis=sqrt(sum(vec(:)**2))
      if(dis.eq.0.d0) then
        rot(:,:)=0.d0
        do i=1,3
          rot(i,i)=1.d0
        enddo
      end if
      dr1(:)=vec(:)/dis
      ivec=minloc(abs(dr1))
      i=ivec(1)
      if(i.eq.1) Then
        dr2(1)=dr1(2)*ex(3)-dr1(3)*ex(2)
        dr2(2)=dr1(3)*ex(1)-dr1(1)*ex(3)
        dr2(3)=dr1(1)*ex(2)-dr1(2)*ex(1)
      else if(i.eq.2) Then
        dr2(1)=dr1(2)*ey(3)-dr1(3)*ey(2)
        dr2(2)=dr1(3)*ey(1)-dr1(1)*ey(3)
        dr2(3)=dr1(1)*ey(2)-dr1(2)*ey(1)
      else if(i.eq.3) Then
        dr2(1)=dr1(2)*ez(3)-dr1(3)*ez(2)
        dr2(2)=dr1(3)*ez(1)-dr1(1)*ez(3)
        dr2(3)=dr1(1)*ez(2)-dr1(2)*ez(1)
      end if
      dr2(:)=dr2(:)/sqrt(sum(dr2(:)**2))
      dr3(1)=dr1(2)*dr2(3)-dr1(3)*dr2(2)
      dr3(2)=dr1(3)*dr2(1)-dr1(1)*dr2(3)
      dr3(3)=dr1(1)*dr2(2)-dr1(2)*dr2(1)
!     ==  the strange mapping is for consistency reasons
      ROT(3,:)=DR1(:)
      ROT(2,:)=-DR2(:)
      ROT(1,:)=DR3(:)
      return
      end subroutine SPHERICAL$align
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PLOTPLGNDR()
      implicit none
      INTEGER(4),PARAMETER :: N=200
      INTEGER(4),PARAMETER :: LMX=9
      INTEGER(4),PARAMETER :: NFIL=888
      INTEGER(4)           :: I
      INTEGER(4)           :: LX
      REAL(8)              :: PLM(LMX),x
      real(8)              :: pi
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      LX=SQRT(REAL(LMX-1,KIND=8)+1.D-10)
      OPEN(NFIL,FILE='PLGNDR.DAT')
      DO I=1,N
!        X=-1.D0+2.D0*REAL(I-1)/REAL(N-1)
!        CALL PLGNDR(LMX,LX,X,PLM)
 X=2.d0*pi*REAL(I-1)/REAL(N-1)
 CALL PLGNDR(LMX,LX,cos(x),PLM)
       WRITE(NFIL,FMT='(50F10.5)')X,PLM
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE testtransform()
      implicit none
      integer(4),parameter  :: nr=250
      integer(4),parameter  :: lx2=12
      integer(4),parameter  :: nx=100
      integer(4),parameter  :: nfil=192
      integer(4)            :: gid
      real(8)               :: dr(3)=(/0.d0,2.d0,3.d0/)
      real(8)               :: r(nr)
      integer(4)            :: lmx1
      integer(4)            :: lmx2
      integer(4)            :: l,i,lm2
      real(8)   ,allocatable:: ylm(:)
      real(8)   ,allocatable:: fin(:,:)
      real(8)   ,allocatable:: fout(:,:)
      real(8)               :: auxi(nr,lx2+1)
      real(8)               :: dr2(3)
      real(8)               :: pi,y0
      real(8)               :: svar,val
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      lmx1=1
      CALL RADIAL$NEW('SHLOG',GID)
      CALL RADIAL$SETI4(GID,'NR',NR)
      CALL RADIAL$SETR8(GID,'DEX',0.05d0)
      CALL RADIAL$SETR8(GID,'R1',1.056d-4)
!
!     ==========================================================================
!     == declare raw function ==================================================
!     ==========================================================================
      CALL RADIAL$r(GID,nr,r)
      ALLOCATE(FIN(NR,lmx1))
      FIN(:,:)=0.D0
      fin(:,1)=exp(-0.25d0*r(:)**2)/y0
!
!     ==========================================================================
!     == loop to check convergence
!     ==========================================================================
      do l=0,lx2
        LMx2=(l+1)**2
        ALLOCATE(FOUT(NR,LMx2))
!       ========================================================================
!       == perform transformation                                             ==
!       ========================================================================
        call spherical$shiftcenter(GID,NR,dr,lmx1,fin,lmx2,fout)
!
!       ========================================================================
!       == reconstruct function                                               ==
!       ========================================================================
        allocate(ylm(lmx2))
        do i=1,nx
          dr2(:)=dr(:)*(-1.5d0+3.d0*real(i-1)/real(nx-1))
          call getylm(lmx2,dr2,ylm)
          svar=0.d0
          do lm2=1,lmx2
            call radial$value(gid,nr,fout(:,lm2),sqrt(sum(dr2**2)),val)
            svar=svar+val*ylm(lm2)
          enddo
          auxi(i,l+1)=svar
        ENDDO
        deallocate(ylm)
        deallocate(fout)
      enddo
!
!     ==========================================================================
!     == write result to file                                                 ==
!     ==========================================================================
      OPEN(nfil,FILE='TESTTRANSFORM.DAT')
      do i=1,nx
        svar=(-1.5d0+3.d0*real(i-1)/real(nx-1))
        write(nfil,fmt='(30f10.5)')svar,auxi(i,:)
      ENDDO
      close(222)
      return
      end

     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine comparetransform()
      implicit none
      integer(4),parameter  :: nr=250
      integer(4),parameter  :: lx1=2
      integer(4),parameter  :: lx2=3
      integer(4)            :: gid
      integer(4)            :: lmx1
      integer(4)            :: lmx2
      real(8)               :: dr(3)=(/0.d0,2.d0,3.d0/)
      real(8)               :: r(nr)
      real(8)   ,allocatable:: fin(:,:),fout1(:,:),fout2(:,:)
      real(8)               :: pi,y0
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      lmx1=(lx1+1)**2
      lmx2=(lx2+1)**2
      CALL RADIAL$NEW('SHLOG',GID)
      CALL RADIAL$SETI4(GID,'NR',NR)
      CALL RADIAL$SETR8(GID,'DEX',0.05d0)
      CALL RADIAL$SETR8(GID,'R1',1.056d-4)
      CALL RADIAL$r(GID,nr,r)
!
!     ==========================================================================
!     == declare raw function ==================================================
!     ==========================================================================
      ALLOCATE(FIN(NR,lmx1))
      FIN(:,:)=0.D0
      fin(:,1)=exp(-0.1d0*r(:)**2)/y0
      if(lmx1.ge.2)fin(:,2)=r(:)*exp(-0.3d0*r(:)**2)/y0
      if(lmx1.ge.8)fin(:,8)=r(:)**2*exp(-0.5d0*r(:)**2)/y0
!
!     ==========================================================================
!     == declare raw function ==================================================
!     ==========================================================================
      ALLOCATE(FOUT1(NR,lmx2))
      ALLOCATE(FOUT2(NR,lmx2))
print*,'peter'
      call spherical$shiftcenter(GID,NR,-dr,lmx1,fin,lmx2,fout1)
      CALL WRITEPHI('FOUT1.DAT',GID,NR,LMx2,FOUT1)
print*,'alex'
      CALL SPHERICAL$TRANSFORM(GID,NR,DR,lmx1,LMx2,FIN,FOUT2)
      CALL WRITEPHI('FOUT2.DAT',GID,NR,LMx2,FOUT2)
      CALL WRITEPHI('FOUTDIFF.DAT',GID,NR,LMx2,FOUT2-FOUT1)
      DEALLOCATE(FIN)     
      DEALLOCATE(FOUT1)     
      DEALLOCATE(FOUT2)
      return
      end
!......................................................................................
subroutine pbscfwrite(nr,lmrx,pot,potrep)
integer(4),intent(in) :: nr
integer(4),intent(in) :: lmrx
real(8)   ,intent(in) :: pot(nr,lmrx)
real(8)   ,intent(in) :: potrep(nr,lmrx)
integer(4)            :: nfil=888
integer(4)            :: lm,ir
open(nfil,file='scfiofile')
rewind(nfil)
write(nfil,*)nr,lmrx
do lm=1,lmrx
do ir=1,nr
write(nfil,*)potrep(ir,lm),pot(ir,lm)
enddo
enddo
close(nfil)
return
end
subroutine pbscfread(nr,lmrx,pot,potrep)
integer(4),intent(in) :: nr
integer(4),intent(in) :: lmrx
real(8)   ,intent(out) :: pot(nr,lmrx)
real(8)   ,intent(out) :: potrep(nr,lmrx)
integer(4)            :: nfil=888
integer(4)            :: lm,ir
integer(4)            :: lmrx1,nr1
open(nfil,file='Poddey080313/scfiofile_alex')
rewind(nfil)
read(nfil,*)nr1,lmrx1
if(nr1.ne.nr.or.lmrx1.ne.lmrx) stop 'stop in pbscfread'
do lm=1,lmrx
do ir=1,nr
read(nfil,*)potrep(ir,lm),pot(ir,lm)
enddo
enddo
close(nfil)
print*,'done'
return
end

