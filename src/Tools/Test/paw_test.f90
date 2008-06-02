!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE SESM_MODULE
TYPE ION_TYPE
  REAL(8)            :: AEZ
  INTEGER(4)         :: LMX
  INTEGER(4)         :: LMRX
  INTEGER(4)         :: LPHIX
  INTEGER(4)         :: LRHOX
  INTEGER(4)         :: LREPX
  INTEGER(4)         :: NB
  INTEGER(4),POINTER :: LOFI(:)
  INTEGER(4),POINTER :: NNOFI(:)
  INTEGER(4)         :: NBG
  INTEGER(4)         :: NCG
  REAL(8)   ,POINTER :: EBG(:)
  REAL(8)   ,POINTER :: FBG(:)
  REAL(8)   ,POINTER :: AEPOT(:,:)   !(NR,LMRX)
  REAL(8)   ,POINTER :: VEMB(:,:)    !(NR,LMRX)
  REAL(8)   ,POINTER :: AERHO(:,:)   !(NR,LMRX)
  REAL(8)   ,POINTER :: XCPOT(:,:)   !(NR,LMRX)
  REAL(8)   ,POINTER :: XCEDEN(:)    !(NR)
  REAL(8)   ,POINTER :: PHINLS(:,:,:)   !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: TPHINLS(:,:,:)   !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: PHI(:,:,:)   !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: TPHI(:,:,:)  !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: HPHI(:,:,:)  !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: TESTHKINPHI(:,:,:)  !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: TESTHHARTREEPHI(:,:,:)  !(NR,LMX,NBG)
  REAL(8)   ,POINTER :: TESTHXCPHI(:,:,:)  !(NR,LMX,NBG)
  REAL(8)            :: ETOT
  REAL(8)            :: EKIN
  REAL(8)            :: EH
  REAL(8)            :: EXC
END TYPE ION_TYPE
END MODULE SESM_MODULE
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE ONEATOM_MODULE
TYPE THISTYPE
  INTEGER(4)         :: GID
  INTEGER(4)         :: NR
  REAL(8)            :: AEZ=-1.D0
  REAL(8)            :: EREF
  INTEGER(4)         :: NB
  INTEGER(4)         :: NC
  REAL(8)   ,POINTER :: ATOMPOT(:)
  REAL(8)            :: EKINC
  REAL(8)   ,POINTER :: RHOCORE(:)
  INTEGER(4)         :: NLC             ! NUMBER CORE ANGULAR MOMENTA
  INTEGER(4),POINTER :: IBBOND(:)       ! BAND INDEX OF LOWEST VALENCE STATE PER L
  REAL(8)   ,POINTER :: ETAPAULI(:,:)   ! VALENCE REPULSION POTENTIAL
  REAL(8)   ,POINTER :: ETACPAULI(:,:)  ! CORE REPULSION POTENTIAL
  REAL(8)   ,POINTER :: VCPAULI(:,:)    ! CORE REPULSION POTENTIAL
  REAL(8)   ,POINTER :: OCPAULI(:,:)    ! CORE REPULSION POTENTIAL
  REAL(8)   ,POINTER :: PHIAT(:,:)      ! ATOMIC WAVE FUNCTIONS
  REAL(8)   ,POINTER :: TPHIAT(:,:)     ! ATOMIC WAVE FUNCTIONS
  REAL(8)   ,POINTER :: UAT(:,:)        ! NODELESS ATOMIC WAVE FUNCTIONS
! == FROM HERE INFORMATION DEPENDING ON CHARGE AND RADIUS
  REAL(8)            :: Q
  REAL(8)            :: RAD
  REAL(8)            :: DEDRAD
  REAL(8)            :: DEDQ
  INTEGER(4),POINTER :: LOFI(:)
  INTEGER(4),POINTER :: NNOFI(:)
  REAL(8)   ,POINTER :: FOFI(:)
  REAL(8)   ,POINTER :: EOFI(:)
  REAL(8)   ,POINTER :: DREL(:)
  REAL(8)   ,POINTER :: RHO(:)     
  REAL(8)   ,POINTER :: AEPOT(:)
  REAL(8)   ,POINTER :: POTH(:)
  REAL(8)   ,POINTER :: POTXC(:)
  REAL(8)   ,POINTER :: UDOT(:,:)     ! HIGHEST NODLESS PHIDOT FUNCTION
  REAL(8)   ,POINTER :: UN(:,:)       ! HIGHEST NODLESS VALENCE STATE PER L
  REAL(8)   ,POINTER :: Uc(:,:)       ! HIGHEST NODLESS core state
  REAL(8)   ,POINTER :: UBOX(:,:)
  REAL(8)   ,POINTER :: TUBOX(:,:)
  REAL(8)   ,POINTER :: PHIBOX(:,:)
  REAL(8)   ,POINTER :: TPHIBOX(:,:)
END TYPE THISTYPE
INTEGER(4)    ,PARAMETER   :: NATX=10
TYPE(THISTYPE),TARGET,SAVE :: THISARR(NATX)
TYPE(THISTYPE),POINTER     :: THIS
LOGICAL(4)    ,PARAMETER   :: TREL=.FALSE.
LOGICAL(4)    ,SAVE        :: TSPECIAL=.FALSE.
integer(4)    ,SAVE        :: gidoneatom=0
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
      REAL(8)    ,ALLOCATABLE   :: EBG(:)
      REAL(8)    ,ALLOCATABLE   :: FBG(:)
      REAL(8)    ,ALLOCATABLE   :: PHI(:,:,:)
      REAL(8)    ,ALLOCATABLE   :: TPHI(:,:,:)
      REAL(8)    ,ALLOCATABLE   :: HPHI(:,:)
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
      CHARACTER(10)             :: ROOTNAME='SESM'
      REAL(8)    ,ALLOCATABLE   :: POTREP(:,:)
      REAL(8)                   :: AEZ
      REAL(8)                   :: DR(3),VAL,DR2(3)
      REAL(8)                   :: CG,PI,Y0
      REAL(8)    ,ALLOCATABLE   :: FIN(:,:),FOUT(:,:)
      INTEGER(4)                :: LX1,LM1,LM2,LM3,IBG,L,LX
!     **************************************************************************
      CALL TRACE$SETL4('ON',.FALSE.)
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
!
!     ==========================================================================        
!     ==  CONNECT PROTOCOLL FILE                                              ==        
!     ==========================================================================        
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL FILEHANDLER$SETFILE('PROT',.FALSE.,TRIM(ROOTNAME)//-'.PROT')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','FORM','FORMATTED')
!
!     ==========================================================================        
!     ==  TEST SPHERICAL                                                      ==        
!     ==========================================================================        
!      CALL PLOTPLGNDR()
!      CALL TESTTRANSFORM()
!      CALL COMPARETRANSFORM()
!STOP
!
!     ==========================================================================        
!     ==  THIS IS FOR SESM                                                    ==        
!     ==========================================================================        
!      CALL BIGONE()
!
!     ==========================================================================
!     ==  THIS IS FOR SIMPLE SESM                                             ==
!     ==========================================================================
      CALL LATTICE()
      RETURN
      END SUBROUTINE DEBUG_PETERS_NONSPH
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BIGONE()
!     **************************************************************************
!     **************************************************************************
      USE SESM_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NR=250
      INTEGER(4),PARAMETER :: NAT=2
      INTEGER(4),PARAMETER :: NBX=19
      INTEGER(4)           :: GID
      TYPE(ION_TYPE)       :: ION(NAT)
      REAL(8)              :: AEZ(NAT)
      REAL(8)              :: DIS
      REAL(8)              :: POS(3,NAT)
      INTEGER(4)           :: NB
      INTEGER(4)           :: IAT
!     **************************************************************************
!
!     ==========================================================================
!     ==  INITIALIZATIONS                                                     ==
!     ==========================================================================
      CALL RADIAL$NEW('LOG',GID)
      CALL RADIAL$SETI4(GID,'NR',NR)
      CALL RADIAL$SETR8(GID,'DEX',0.05D0)
      CALL RADIAL$SETR8(GID,'R1',1.056D-4)
      CALL DFT$SETI4('TYPE',1)
      CALL DFT$SETL4('SPIN',.FALSE.)
!
!     ==========================================================================
!     ==  SYSTEM SPECIFIC                                                     ==
!     ==========================================================================
      ION(:)%AEZ=7.D0
      ION(:)%LPHIX=2
      ION(:)%LRHOX=2
      ION(:)%LREPX=2*ION(:)%LPHIX
      ION(:)%NB=3
      NB=3
      DO IAT=1,NAT
        ALLOCATE(ION(IAT)%LOFI(NB)) 
        ION(IAT)%LOFI(:)=(/0,0,1/)
        ALLOCATE(ION(IAT)%NNOFI(NB))
        ION(IAT)%NNOFI(:)=(/0,1,0/)
      ENDDO
      DIS=1.1D0/0.529177D0
! DIS=1.5D0/0.529177D0
      POS(:,:)=0.D0
      POS(3,2)=DIS
!
!     ==========================================================================
!     ==  CALL SESM                                                           ==
!     ==========================================================================
      CALL SESM_PETER(GID,NR,NAT,POS,ION)
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SESM_PETER(GID,NR,NAT,POS,ION)
!     **************************************************************************
!     **************************************************************************
      USE ONEATOM_MODULE
      USE SESM_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: POS(3,NAT)
      TYPE(ION_TYPE),INTENT(INOUT) :: ION(NAT)
      INTEGER(4)            :: NB
      INTEGER(4)            :: NC
      INTEGER(4),ALLOCATABLE:: LOFI(:)
      INTEGER(4),ALLOCATABLE:: NNOFI(:)
      INTEGER(4)            :: LMRX,LMRX1,LMRX2,LMX,LMX1,LMX2
      INTEGER(4)            :: LMAX
      REAL(8)               :: RAD(NAT)
      REAL(8)               :: DETOT,ETOT,EHARTREE
      REAL(8)               :: DEDQ
      REAL(8)               :: DR(3)
      REAL(8)               :: DEDRAD
      REAL(8)   ,ALLOCATABLE:: TESTPOTXC(:,:)
      REAL(8)   ,ALLOCATABLE:: TESTPOTHARTREE(:,:)
      REAL(8)   ,ALLOCATABLE:: POT(:,:)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:)
      REAL(8)   ,ALLOCATABLE:: RHOS(:,:)
      REAL(8)   ,ALLOCATABLE:: RHOEXT(:,:)
      REAL(8)   ,ALLOCATABLE:: VEMB(:,:)
      REAL(8)   ,ALLOCATABLE:: FOUT(:,:)
      REAL(8)   ,ALLOCATABLE:: FOUT2(:,:)
      REAL(8)   ,ALLOCATABLE:: XCEDEN(:,:)
      REAL(8)   ,ALLOCATABLE:: AERHO(:,:)
      REAL(8)               :: WEIGHT(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: EDEN(NR)
      INTEGER(4)            :: NCG
      INTEGER(4)            :: NBG
      REAL(8)   ,ALLOCATABLE:: EBG(:)
      REAL(8)   ,ALLOCATABLE:: FBG(:)
      REAL(8)   ,ALLOCATABLE:: PHI(:,:,:)
      REAL(8)   ,ALLOCATABLE:: TPHI(:,:,:)
      REAL(8)   ,ALLOCATABLE:: HPHI(:,:,:)
      REAL(8)   ,ALLOCATABLE:: TESTHKINPHI(:,:,:)
      REAL(8)   ,ALLOCATABLE:: TESTHHARTREEPHI(:,:,:)
      REAL(8)   ,ALLOCATABLE:: TESTHXCPHI(:,:,:)
      REAL(8)               :: DREL(NR)
      REAL(8)               :: EKIN,EXC,EH,EB
      REAL(8)               :: SVAL,HVAL
      REAL(8)               :: TESTHKINVAL,TESTHHARTREEVAL,TESTHXCVAL
      INTEGER(4)            :: IAT,IAT1,IAT2,IBG,IBG1,IBG2,ICHI1,ICHI2,ICHI
      INTEGER(4)            :: LM1,LM2,LM3,IB,I
      INTEGER(4)            :: NCHI
      CHARACTER(32)         :: STRING
      REAL(8)  ,ALLOCATABLE :: S(:,:)
      REAL(8)  ,ALLOCATABLE :: H(:,:)
      REAL(8)  ,ALLOCATABLE :: TESTHKIN(:,:)
      REAL(8)  ,ALLOCATABLE :: TESTHHARTREE(:,:)
      REAL(8)  ,ALLOCATABLE :: TESTHXC(:,:)
      REAL(8)  ,ALLOCATABLE :: EIG(:)
      REAL(8)  ,ALLOCATABLE :: U(:,:)
      REAL(8)               :: CG,PI,Y0,FOURPI
      REAL(8)               :: SVAR
      REAL(8)               :: EBAND,EBANDISO
      REAL(8)               :: OCC
      REAL(8)               :: R(NR)
      INTEGER(4)            :: LX1
      INTEGER(4)            :: LX,LRX,LMREPX,IR
      INTEGER(4)            :: NFILO
      REAL(8)               :: EV
      LOGICAL(4),PARAMETER  :: TSECONDLINE=.FALSE.
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(4.D0*PI)
      CALL CONSTANTS('EV',EV)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  CALCULATE ATOM IN A BOX TO OBTAIN THE REPULSIVE POTENTIAL           ==
!     ==========================================================================
      ETOT=0.D0
      DO IAT1=1,NAT
        CALL ONEATOM$NEW(IAT1,ION(IAT1)%AEZ)
        RAD(IAT1)=25.D0
        DO IAT2=1,NAT
          IF(IAT2.EQ.IAT1) CYCLE
          RAD(IAT1)=MIN(RAD(IAT1),SQRT(SUM((POS(:,IAT2)-POS(:,IAT1))**2)))
        ENDDO
        RAD(IAT1)=0.7D0*RAD(IAT1)
PRINT*,'RAD(IAT1)',RAD(IAT1)
        CALL ONEATOM(IAT1,0.D0,RAD(IAT1),DETOT,DEDQ,DEDRAD)
        ETOT=ETOT+DETOT
PRINT*,'IAT,ETOT ',IAT1,DETOT
!       == SELECT NUMBER OF CORE STATES
        NC=THISARR(IAT1)%NC
        ION(IAT1)%NCG=SUM(2*THISARR(IAT1)%LOFI(1:NC)+1)
      ENDDO
!!$      CALL WRITEPHI('ETAPAULI.DAT',GID,NR,2,THISARR(1)%ETAPAULI(:,:))
!!$      CALL WRITEPHI('ETACPAULI.DAT',GID,NR,2,THISARR(1)%ETACPAULI(:,:))
!!$      CALL WRITEPHI('VCPAULI.DAT',GID,NR,2,THISARR(1)%VCPAULI(:,:))
!!$      CALL WRITEPHI('OCPAULI.DAT',GID,NR,2,THISARR(1)%OCPAULI(:,:))
!!$      CALL WRITEPHI('UN.DAT',GID,NR,2,THISARR(1)%UN(:,:))
!!$      CALL WRITEPHI('UDOT.DAT',GID,NR,2,THISARR(1)%UDOT(:,:))
!!$OPEN(991,FILE='ETAPAULIS.DAT')
!!$DO IAT1=1,NR
!!$  WRITE(991,*)R(IAT1),THISARR(1)%ETAPAULI(IAT1,1)
!!$ENDDO
!!$CLOSE(991)
!!$OPEN(991,FILE='SETUP_PETER_.ETA_AND_VAE_HB_OFFSITE_L1.NXY')
!!$DO IAT1=1,NR
!!$!  READ(991,*)SVAR,THISARR(1)%ETAPAULI(IAT1,1)
!!$ENDDO
!!$CLOSE(991)
!!$!THISARR(2)%ETAPAULI(:,1)=THISARR(1)%ETAPAULI(:,1)
!!$
!
!     ==========================================================================
!     == CONSTRUCT EMBEDDING POTENTIAL                                        ==
!     ==========================================================================
PRINT*,'CONSTRUCTING EMBEDDING POTENTIAL..........................'
      DO IAT1=1,NAT
        LMX =(ION(IAT1)%LPHIX+1)**2
        LMRX=(ION(IAT1)%LRHOX+1)**2
        LMREPX=(ION(IAT1)%LREPX+1)**2
        ALLOCATE(VEMB(NR,LMREPX))
        ALLOCATE(FOUT(NR,LMREPX))
        ALLOCATE(FOUT2(NR,LMREPX))
        VEMB(:,:)=0.D0
        ALLOCATE(RHOEXT(NR,LMREPX))
        RHOEXT(:,:)=0.D0
!
        DO IAT2=1,NAT
          IF(IAT2.EQ.IAT1) CYCLE
          DR(:)=POS(:,IAT1)-POS(:,IAT2)
!         == ADD UP EFFECTIVE POTENTIALS OF IONS IN A BOX ========================
          AUX(:)=THISARR(IAT2)%ETAPAULI(:,1)+THISARR(IAT2)%POTH(:)+THISARR(IAT2)%POTXC(:)
          CALL SPHERICAL$SHIFTCENTER(GID,NR,DR,1,AUX,LMREPX,FOUT)
          VEMB(:,:)=VEMB(:,:)+FOUT(:,:)
        ENDDO
!
        IF(TSECONDLINE) THEN
          DO IAT2=1,NAT
            IF(IAT2.EQ.IAT1) CYCLE
            DR(:)=POS(:,IAT1)-POS(:,IAT2)
!           == ADD UP DENSITIES ====================================================
            AUX(:)=THISARR(IAT2)%RHO(:)
            CALL SPHERICAL$SHIFTCENTER(GID,NR,DR,1,AUX,LMREPX,FOUT)
            RHOEXT(:,:)=RHOEXT(:,:)+FOUT(:,:)
!           == SUBTRACT XC-POTENTIAL OF EXTERNAL DENSITIES FROM VEMB================
            CALL AUGMENTATION_XC(GID,NR,LMREPX,1,FOUT,SVAR,FOUT2,AUX)
            DO IR=1,NR
              IF(FOUT(IR,1)*Y0.LT.1.D-6) FOUT2(IR,:)=0.D0
            ENDDO
            VEMB(:,:)=VEMB(:,:)-FOUT2(:,:)
          ENDDO
!         ==  SUBTRACT ONSITE XC-POTENTIAL ========================================
          VEMB(:,1)=VEMB(:,1)-THISARR(IAT1)%POTXC(:)
!         == CONSTRUCT XC-POTENTIAL OF SUPERIMPOSED DENSITIES
          FOUT(:,:)=RHOEXT(:,:)
          FOUT(:,1)=FOUT(:,1)+THISARR(IAT1)%RHO(:)
          CALL AUGMENTATION_XC(GID,NR,LMREPX,1,FOUT,SVAR,FOUT2,AUX)
          VEMB(:,:)=VEMB(:,:)+FOUT2(:,:)
        END IF
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
!      CALL WRITEPHI('VEMBEDDING.DAT',GID,NR,LMREPX,ION(1)%VEMB)
!
!     ============================================================================
!     == SELF-CONSISTENT CALCULATION OF THE DEFORMED ION                        ==
!     ============================================================================
PRINT*,'SCF OF DEFORMED ION...........................................'
      ETOT=0.D0
      DO IAT1=1,NAT
        LMX =(ION(IAT1)%LPHIX+1)**2
        LMRX=(ION(IAT1)%LRHOX+1)**2
        LMREPX=(ION(IAT1)%LREPX+1)**2
        LX=SQRT(REAL(LMX-1,KIND=8)+1.D-10)
        NC=THISARR(IAT1)%NC
        LMAX=THISARR(IAT1)%NLC
!
!       == COLLECT DATA ==========================================================
        NB=ION(IAT1)%NB
        ALLOCATE(LOFI(NB))
        ALLOCATE(NNOFI(NB))
        LOFI(:)=ION(IAT1)%LOFI(:)
        NNOFI(:)=ION(IAT1)%NNOFI(:)
!
!       ==  ALLOCATE ARRAYS ======================================================
        ALLOCATE(AERHO(NR,LMRX))
        NCG=ION(IAT1)%NCG
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
!POT(:,1)=THISARR(IAT1)%AEPOT(:)
        DO IB=NC+1,NB
          DO I=1,NC
            IF(LOFI(I).NE.LOFI(IB)) CYCLE
            NNOFI(IB)=NNOFI(IB)-1
          ENDDO
        ENDDO
!
!       ==========================================================================
!       == PERFORM SCF CALCULATION OF DEFORMED ION                              ==
!       ==========================================================================
        DREL(:)=0.D0
        CALL SCF(GID,NR,LMRX,LMX,LMREPX,ION(IAT1)%AEZ,DREL,ION(IAT1)%VEMB &
       &        ,LMAX,THISARR(IAT1)%VCPAULI(:,1:LMAX),THISARR(IAT1)%OCPAULI(:,1:LMAX) &
       &        ,NC,THISARR(IAT1)%PHIAT(:,1:NC),THISARR(IAT1)%TPHIAT(:,1:NC) &
       &        ,NB,LOFI,NNOFI,POT,AERHO &
       &        ,NCG,NBG,EBG,FBG,PHI,TPHI)
!
        ALLOCATE(ION(IAT1)%PHINLS(NR,LMX,NBG))
        ALLOCATE(ION(IAT1)%TPHINLS(NR,LMX,NBG))
        ION(IAT1)%PHINLS=PHI
        ION(IAT1)%TPHINLS=TPHI
        CALL SCF_COREORTHO(GID,NR,NC,LOFI &
       &             ,THISARR(IAT1)%PHIAT(:,1:NC),THISARR(IAT1)%TPHIAT(:,1:NC) &
       &             ,NBG,NCG,LMX,PHI,TPHI)
!
!       ==========================================================================
!       == CALCULATE TOTAL ENERGY                                               ==
!       ==========================================================================
        CALL DEFORM$ETOT(GID,NR,LMRX,LMX,LMREPX,NBG,ION(IAT1)%AEZ &
     &                  ,FBG,EBG,PHI,TPHI,DETOT,EB,EKIN,EH,EXC)
        ETOT=ETOT+DETOT-THISARR(IAT1)%EREF
!
!       == REPORT ENERGIES =======================================================
        WRITE(NFILO,FMT='(72("="))')
        WRITE(NFILO,FMT='(72("="),T20," ENERGIES OF THE DEFORMED ION ")')
        WRITE(NFILO,FMT='("EDEFORM =",F10.5," H")')DETOT-THISARR(IAT1)%EREF
        WRITE(NFILO,FMT='("EDEFORM =",F10.5," EV")')(DETOT-THISARR(IAT1)%EREF)/EV
        WRITE(NFILO,FMT='("EION    =",F10.5," H")')DETOT
        WRITE(NFILO,FMT='("EKIN    =",F10.5," H")')EKIN
        WRITE(NFILO,FMT='("EHARTREE=",F10.5," H")')EH
        WRITE(NFILO,FMT='("EXC     =",F10.5," H")')EXC
        DO IBG=1,NBG
          WRITE(NFILO,FMT='(I3,"OCC=",F5.2," E[H]=",F10.5," E[EV]=",F10.5)') &
       &                  IBG,FBG(IBG),EBG(IBG),EBG(IBG)/EV
        ENDDO
        WRITE(NFILO,FMT='(72("="))')
!
!       == MAP RESULT ONTO ION STRUCTURE =========================================
!        CALL WRITEPHI('POTIN.DAT',GID,NR,LMRX,POT)
        ALLOCATE(ION(IAT1)%EBG(NBG))
        ALLOCATE(ION(IAT1)%FBG(NBG))
        ALLOCATE(ION(IAT1)%PHI(NR,LMX,NBG))
        ALLOCATE(ION(IAT1)%TPHI(NR,LMX,NBG))
        ALLOCATE(ION(IAT1)%AEPOT(NR,LMRX))
        ALLOCATE(ION(IAT1)%AERHO(NR,LMRX))
        ION(IAT1)%NBG=NBG
        ION(IAT1)%EBG(:)=EBG(:)
        ION(IAT1)%FBG(:)=FBG(:)
        ION(IAT1)%PHI(:,:,:)=PHI(:,:,:)
        ION(IAT1)%TPHI(:,:,:)=TPHI(:,:,:)
        ION(IAT1)%AEPOT(:,:)=POT(:,:)
        ION(IAT1)%AERHO(:,:)=AERHO(:,:)
!       == WRAP UP ===============================================================
        DEALLOCATE(EBG)
        DEALLOCATE(FBG)
        DEALLOCATE(PHI)
        DEALLOCATE(TPHI)
        DEALLOCATE(LOFI)
        DEALLOCATE(NNOFI)
        DEALLOCATE(POT)
        DEALLOCATE(AERHO)
      ENDDO
!
!     ============================================================================
!     == INSPECT RESULTS                                                        ==
!     ============================================================================
!!$      IAT1=1
!!$      LMRX=(ION(IAT1)%LRHOX+1)**2
!!$      DO IBG=1,NBG
!!$        WRITE(STRING,*)IBG
!!$        STRING='PHI'//TRIM(ADJUSTL(STRING))//'.DAT'
!!$        CALL WRITEPHI(STRING,GID,NR,LMX,ION(IAT1)%PHI(:,:,IBG))
!!$        WRITE(STRING,*)IBG
!!$        STRING='PHINLS'//TRIM(ADJUSTL(STRING))//'.DAT'
!!$        CALL WRITEPHI(STRING,GID,NR,LMX,ION(IAT1)%PHINLS(:,:,IBG))
!!$      ENDDO
!!$      CALL WRITEPHI('AEPOT.DAT',GID,NR,LMRX,ION(IAT1)%AEPOT)
!!$      CALL WRITEPHI('AERHO.DAT',GID,NR,LMRX,ION(IAT1)%AERHO)
!
!     ==========================================================================
!     == DETERMINE HARTREE ENERGY                                             ==
!     ==========================================================================
PRINT*,'CALCULATING HARTREE ENERGY ...'
      EHARTREE=0.D0
      DO IAT1=1,NAT
         LMRX1=(ION(IAT1)%LRHOX+1)**2
        DO IAT2=IAT1+1,NAT
          LMRX2=(ION(IAT2)%LRHOX+1)**2
          CALL HARTREEPAIR(GID,NR,POS(:,IAT2)-POS(:,IAT1) &
     &                    ,LMRX1,ION(IAT1)%AEZ,ION(IAT1)%AERHO &
     &                    ,LMRX2,ION(IAT2)%AEZ,ION(IAT2)%AERHO,DETOT)
          EHARTREE=EHARTREE+DETOT
        ENDDO
      ENDDO
      PRINT*,'EHARTREE ',EHARTREE
      WRITE(NFILO,FMT='("HARTREE ENERGY                  ",F10.5," H")')EHARTREE
      ETOT=ETOT+EHARTREE
      PRINT*,'ETOT AFTER EHARTREE ',ETOT
      WRITE(NFILO,FMT='("ETOT AFTER CALCULATING EHARTREE ",F10.5," H")')ETOT
      WRITE(NFILO,FMT='("ETOT AFTER CALCULATING EHARTREE ",F10.5," EV")')ETOT/EV
!
!     ==========================================================================
!     == CALCULATE XC-CORRECTION                                              ==
!     ==========================================================================
      DO IAT1=1,NAT
        LMRX=(ION(IAT1)%LRHOX+1)**2
        ALLOCATE(ION(IAT1)%XCPOT(NR,LMRX))
        ALLOCATE(ION(IAT1)%XCEDEN(NR))
        CALL AUGMENTATION_XC(GID,NR,LMRX,1,ION(IAT1)%AERHO,SVAR &
     &                      ,ION(IAT1)%XCPOT,ION(IAT1)%XCEDEN)
!       == CUT OUT LOW DENSITY REGION ==========================================
        DO IR=1,NR
          IF(ION(IAT1)%AERHO(IR,1)*Y0.LT.1.D-6) ION(IAT1)%XCPOT(IR,:)=0.D0
        ENDDO
       ENDDO
!
!     ==========================================================================
!     == CALCULATE H|PHI>
!     ==========================================================================
PRINT*,'CALCULATING H|\PHI>....'
      EXC=0.D0
      DO IAT1=1,NAT
        NBG=ION(IAT1)%NBG
        LMX =(ION(IAT1)%LPHIX+1)**2
        LMRX=(ION(IAT1)%LRHOX+1)**2
        ALLOCATE(PHI(NR,LMX,NBG))
        ALLOCATE(HPHI(NR,LMX,NBG))
        ALLOCATE(TESTHKINPHI(NR,LMX,NBG))
        ALLOCATE(TESTHHARTREEPHI(NR,LMX,NBG))
        ALLOCATE(TESTHXCPHI(NR,LMX,NBG))
        PHI(:,:,:)=ION(IAT1)%PHI(:,:,:)
!
!       ========================================================================
!       == CONSTRUCT SUM OF POTENTIALS AND SUM OF DENSITIES                   ==
!       ========================================================================
        ALLOCATE(POT(NR,LMRX))
        ALLOCATE(TESTPOTHARTREE(NR,LMRX))
        ALLOCATE(TESTPOTXC(NR,LMRX))
        ALLOCATE(RHO(NR,LMRX))
        ALLOCATE(RHOS(NR,LMRX))
        ALLOCATE(XCEDEN(NR,LMRX))
        ALLOCATE(FOUT(NR,LMRX))
        POT(:,:)=0.D0
        RHO(:,:)=0.D0
        RHOS(:,:)=0.D0
        XCEDEN(:,:)=0.D0
        DO IAT2=1,NAT
          IF(IAT2.EQ.IAT1) CYCLE
          LMRX2=(ION(IAT2)%LRHOX+1)**2
!         == TOTAL POTENTIAL EXCEPT XC POTENTIAL ===============================
          CALL SPHERICAL$SHIFTCENTER(GID,NR,POS(:,IAT1)-POS(:,IAT2) &
    &                  ,LMRX2,ION(IAT2)%AEPOT-ION(IAT2)%XCPOT,LMRX,FOUT(:,:))
          POT(:,:)=POT(:,:)+FOUT(:,:)
!
!         ==  TOTAL DENSITY ====================================================
          CALL SPHERICAL$SHIFTCENTER(GID,NR,POS(:,IAT1)-POS(:,IAT2) &
    &                           ,LMRX2,ION(IAT2)%AERHO,LMRX,FOUT(:,:))
          RHO(:,:)=RHO(:,:)+FOUT(:,:)
!
!         == WEIGHTING FUNCTION FOR XC ENERGY DENSITY ==========================
          AUX(:)=EXP(-R(:)**2)
          CALL SPHERICAL$SHIFTCENTER(GID,NR,POS(:,IAT1)-POS(:,IAT2) &
    &                           ,1,AUX,LMRX,FOUT(:,:))
          RHOS(:,:)=RHOS(:,:)+FOUT(:,:)
!
!         == XC ENERGY DENSITY                       ========================
          CALL SPHERICAL$SHIFTCENTER(GID,NR,POS(:,IAT1)-POS(:,IAT2) &
    &                           ,1,ION(IAT2)%XCEDEN,LMRX,FOUT)
          XCEDEN(:,:)=XCEDEN(:,:)-FOUT(:,:)
        ENDDO
        POT(:,:)=ION(IAT1)%AEPOT(:,:)-ION(IAT1)%XCPOT(:,:)+POT(:,:)
        TESTPOTHARTREE(:,:)=POT(:,:)
        RHO(:,:)=ION(IAT1)%AERHO(:,:)+RHO(:,:)
        CALL AUGMENTATION_XC(GID,NR,LMRX,1,RHO,SVAR,FOUT,AUX)
!       == CUT OUT LOW DENSITY REGION ==========================================
        DO IR=1,NR
          IF(RHO(IR,1)*Y0.LT.1.D-6) FOUT(IR,:)=0.D0
        ENDDO
! CALL WRITEPHI('XCPOT.DAT',GID,NR,LMRX,FOUT)
        POT=POT+FOUT
        TESTPOTXC(:,:)=FOUT(:,:)
        XCEDEN(:,1)=XCEDEN(:,1)-ION(IAT1)%XCEDEN(:)+AUX
! CALL WRITEPHI('XCEDEN.DAT',GID,NR,LMRX,XCEDEN)
! CALL WRITEPHI('XCEDEN1.DAT',GID,NR,1,ION(1)%XCEDEN)
! CALL WRITEPHI('XCEDEN2.DAT',GID,NR,1,ION(2)%XCEDEN)
!
!       ========================================================================
!       == CORRECT XC ENERGY                                                  ==
!       ========================================================================
!       == WEIGHTING
        AUX(:)=EXP(-R(:)**2)
        RHOS(:,1)=AUX(:)+RHOS(:,1)
        AUX(:)=AUX(:)/(RHOS(:,1)*Y0)
        DO LM1=2,LMRX
          RHOS(:,LM1)=RHOS(:,LM1)/(RHOS(:,1)*Y0)
          RHOS(:,LM1)=RHOS(:,LM1)*AUX(:)
        ENDDO
        RHOS(:,1)=AUX(:)
        AUX(:)=0.D0
        DO LM1=1,LMRX
          AUX(:)=AUX(:)+RHOS(:,LM1)*XCEDEN(:,LM1)
        ENDDO
! CALL WRITEPHI('WXCEDEN.DAT',GID,NR,1,AUX)
! CALL WRITEPHI('WGHT.DAT',GID,NR,LMRX,RHOS)
        CALL RADIAL$INTEGRAL(GID,NR,AUX*R(:)**2,SVAR)
        EXC=EXC+FOURPI*SVAR*Y0
!
PRINT*,' EXC?? ',FOURPI*SVAR*Y0
!!$IF(IAT1.EQ.1) THEN
!!$  CALL WRITEPHI('VTOTC1.DAT',GID,NR,LMRX,POT)
!!$  CALL WRITEPHI('RHOTOTC1.DAT',GID,NR,LMRX,RHO)
!!$  CALL WRITEPHI('RHOSC1.DAT',GID,NR,LMRX,RHOS)
!!$  CALL WRITEPHI('RHOS1.DAT',GID,NR,1,ION(IAT1)%AERHO(:,1))
!!$ELSE IF(IAT1.EQ.2) THEN
!!$  CALL WRITEPHI('VTOTC2.DAT',GID,NR,LMRX,POT)
!!$  CALL WRITEPHI('RHOTOTC2.DAT',GID,NR,LMRX,RHO)
!!$  CALL WRITEPHI('RHOSC2.DAT',GID,NR,LMRX,RHOS)
!!$END IF
        HPHI(:,:,:)=0.D0
        TESTHHARTREEPHI(:,:,:)=0.D0
        TESTHXCPHI(:,:,:)=0.D0
        DO IBG=1,NBG
          DO LM1=1,LMX
            DO LM2=1,LMRX
              DO LM3=1,LMX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.EQ.0.D0) CYCLE
                HPHI(:,LM1,IBG)=HPHI(:,LM1,IBG)+CG*POT(:,LM2)*PHI(:,LM3,IBG)
                TESTHHARTREEPHI(:,LM1,IBG)=TESTHHARTREEPHI(:,LM1,IBG)+CG*TESTPOTHARTREE(:,LM2)*PHI(:,LM3,IBG)
                TESTHXCPHI(:,LM1,IBG)=TESTHXCPHI(:,LM1,IBG)+CG*TESTPOTXC(:,LM2)*PHI(:,LM3,IBG)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(PHI)
        DEALLOCATE(POT)
        DEALLOCATE(TESTPOTHARTREE)
        DEALLOCATE(TESTPOTXC)
        DEALLOCATE(RHO)
        DEALLOCATE(RHOS)
        DEALLOCATE(XCEDEN)
        DEALLOCATE(FOUT)
        ALLOCATE(ION(IAT1)%HPHI(NR,LMX,NBG))
        ALLOCATE(ION(IAT1)%TESTHKINPHI(NR,LMX,NBG))
        ALLOCATE(ION(IAT1)%TESTHHARTREEPHI(NR,LMX,NBG))
        ALLOCATE(ION(IAT1)%TESTHXCPHI(NR,LMX,NBG))
        ION(IAT1)%HPHI(:,:,:)=ION(IAT1)%TPHI(:,:,:)+HPHI(:,:,:)
        ION(IAT1)%TESTHKINPHI(:,:,:)=ION(IAT1)%TPHI(:,:,:)
        ION(IAT1)%TESTHHARTREEPHI(:,:,:)=TESTHHARTREEPHI(:,:,:)
        ION(IAT1)%TESTHXCPHI(:,:,:)=TESTHXCPHI(:,:,:)
        DEALLOCATE(HPHI)
        DEALLOCATE(TESTHHARTREEPHI)
        DEALLOCATE(TESTHXCPHI)
        DEALLOCATE(TESTHKINPHI)
      ENDDO
PRINT*,'DELTA EXC ',EXC
!      ETOT=ETOT+EXC
      PRINT*,'ETOT AFTER EXCHANGE-CORRELATION ENERGY ',ETOT
!
!     ==========================================================================
!     == CALCULATE HAMILTON AND OVERLAP MATRIX                                ==
!     ==========================================================================
PRINT*,'CALCULATING HAMILTON AND OVERLAPMATRIX ...'
      NCHI=0
      DO IAT1=1,NAT
        NCHI=NCHI+ION(IAT1)%NBG
      ENDDO
      ALLOCATE(S(NCHI,NCHI))
      ALLOCATE(H(NCHI,NCHI))
      ALLOCATE(TESTHKIN(NCHI,NCHI))
      ALLOCATE(TESTHXC(NCHI,NCHI))
      ALLOCATE(TESTHHARTREE(NCHI,NCHI))
      S(:,:)=0.D0
      H(:,:)=0.D0
      ICHI1=0
      DO IAT1=1,NAT
        LMX1=(ION(IAT1)%LPHIX+1)**2
        DO IBG1=1,ION(IAT1)%NBG
          H(ICHI1+IBG1,ICHI1+IBG1)=ION(IAT1)%EBG(IBG1)
          S(ICHI1+IBG1,ICHI1+IBG1)=1.D0
        ENDDO
        ICHI2=0
        DO IAT2=1,NAT
!          IF(IAT1.GE.IAT2) THEN
          LMX2=(ION(IAT2)%LPHIX+1)**2
          IF(IAT1.GE.0) THEN
            DO IBG1=1,ION(IAT1)%NBG
              DO IBG2=1,ION(IAT2)%NBG 
                CALL OVERLAP(GID,NR,POS(:,IAT2)-POS(:,IAT1) &
     &            ,LMX1,ION(IAT1)%PHI(:,:,IBG1) &
     &            ,LMX2,ION(IAT2)%HPHI(:,:,IBG2) &
     &                 ,ION(IAT2)%TESTHKINPHI(:,:,IBG2) &
     &                 ,ION(IAT2)%TESTHHARTREEPHI(:,:,IBG2) &
     &                 ,ION(IAT2)%TESTHXCPHI(:,:,IBG2) &
     &                 ,ION(IAT2)%PHI(:,:,IBG2) &
     &                 ,HVAL,TESTHKINVAL,TESTHHARTREEVAL,TESTHXCVAL,SVAL)
                H(ICHI1+IBG1,ICHI2+IBG2)=HVAL
                TESTHKIN(ICHI1+IBG1,ICHI2+IBG2)=TESTHKINVAL
                TESTHHARTREE(ICHI1+IBG1,ICHI2+IBG2)=TESTHHARTREEVAL
                TESTHXC(ICHI1+IBG1,ICHI2+IBG2)=TESTHXCVAL
!                H(ICHI2+IBG2,ICHI1+IBG1)=HVAL
                S(ICHI1+IBG1,ICHI2+IBG2)=SVAL
!                S(ICHI2+IBG2,ICHI1+IBG1)=SVAL
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
      ENDDO
      DEALLOCATE(TESTHKIN)
      DEALLOCATE(TESTHXC)
      DEALLOCATE(TESTHHARTREE)
!
!     ==========================================================================
!     == SOLVE SCHROEDINEGR EQUATION                                          ==
!     ==========================================================================
      S=0.5D0*(S+TRANSPOSE(S))
      H=0.5D0*(H+TRANSPOSE(H))
      ALLOCATE(U(NCHI,NCHI))
      ALLOCATE(EIG(NCHI))
      CALL LIB$GENERALEIGENVALUER8(NCHI,H,S,EIG,U)
!
!     ==========================================================================
!     == CALCULATE DIFFERENCE IN BANDSTUCTURE ENERGY                          ==
!     ==========================================================================
      EBANDISO=0.D0
      SVAR=0.D0
      ICHI=0
      DO IAT=1,NAT
        SVAR=SVAR+ION(IAT)%AEZ
        NBG=ION(IAT)%NBG
        DO IB=1,NBG
          ICHI=ICHI+1
          EBANDISO=EBANDISO+ION(IAT)%FBG(IB)*H(ICHI,ICHI)
        ENDDO 
      ENDDO
      WRITE(NFILO,FMT='(72("="),T30," TIGHT-BINDING ")')
      EBAND=0.D0
      DO IB=1,NCHI
        OCC=MIN(2.D0,SVAR)
        WRITE(NFILO,FMT='(I3,"OCC=",F5.2," E[H]=",F10.5," E[EV]=",F10.5)')IB,OCC,EIG(IB),EIG(IB)/EV
        SVAR=SVAR-OCC
        EBAND=EBAND+EIG(IB)*OCC
      ENDDO
      DEALLOCATE(U)
      WRITE(NFILO,FMT='("BANDSTRUCTURE ENERGY(ISOLATED)   ",F10.5," H")')EBANDISO
      WRITE(NFILO,FMT='("BANDSTRUCTURE ENERGY             ",F10.5," H")')EBAND
      WRITE(NFILO,FMT='("BANDSTRUCTURE ENERGY(DIFFERENCE) ",F10.5," H")')EBAND-EBANDISO
      WRITE(NFILO,FMT='("BANDSTRUCTURE ENERGY(DIFFERENCE) ",F10.5," EV")')(EBAND-EBANDISO)/EV
      ETOT=ETOT+EBAND-EBANDISO
      WRITE(NFILO,FMT='("ETOT AFTER BANDSTRUCTURE ENERGY  ",F10.5," H")')ETOT
      WRITE(NFILO,FMT='("ETOT AFTER BANDSTRUCTURE ENERGY  ",F10.5," EV")')ETOT/EV
!
!     ==========================================================================
!     == WRAP UP                                                              ==
!     ==========================================================================
      DEALLOCATE(S)
      DEALLOCATE(H)
      DEALLOCATE(EIG)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE OVERLAP(GID,NR,DR21,LMX1,PHI1,LMX2,HPHI2 &
     &                  ,TESTHKINPHI2,TESTHHARTREEPHI2,TESTHXCPHI2 &
     &                  ,PHI2,H,TESTHKIN,TESTHHARTREE,TESTHXC,S)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: GID
      INTEGER(4) ,INTENT(IN)    :: NR
      REAL(8)    ,INTENT(IN)    :: DR21(3)
      INTEGER(4) ,INTENT(IN)    :: LMX1       ! X#(POTENTIAL ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)    :: LMX2       ! X#(POTENTIAL ANGULAR MOMENTA)
      REAL(8)    ,INTENT(IN)    :: PHI1(NR,LMX1)
      REAL(8)    ,INTENT(IN)    :: HPHI2(NR,LMX2)
      REAL(8)    ,INTENT(IN)    :: TESTHKINPHI2(NR,LMX2)
      REAL(8)    ,INTENT(IN)    :: TESTHHARTREEPHI2(NR,LMX2)
      REAL(8)    ,INTENT(IN)    :: TESTHXCPHI2(NR,LMX2)
      REAL(8)    ,INTENT(IN)    :: PHI2(NR,LMX2)
      REAL(8)    ,INTENT(OUT)   :: H
      REAL(8)    ,INTENT(OUT)   :: TESTHKIN,TESTHHARTREE,TESTHXC
      REAL(8)    ,INTENT(OUT)   :: S
      REAL(8)                   :: SHIFTEDPHI(NR,LMX2)
      REAL(8)                   :: AUX(NR)
      INTEGER(4)                :: LM      
      REAL(8)                   :: R(NR)
      REAL(8)                   :: DIS
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  SHIFT DENSITY RHO2 TO POSITION 1                                    ==
!     ==========================================================================
      DIS=SQRT(SUM(DR21**2))
      IF(DIS.GT.1.D-5) THEN
        CALL SPHERICAL$SHIFTCENTER(GID,NR,DR21,LMX1,PHI1,LMX2,SHIFTEDPHI)
      ELSE
        LM=MIN(LMX1,LMX2)
        SHIFTEDPHI=0.D0
        SHIFTEDPHI(:,:LM)=PHI1(:,:LM)
      ENDIF
!
!     ==========================================================================
!     ==  CALCULATE ENERGY                                                    ==
!     ==========================================================================
      AUX(:)=0.D0
      DO LM=1,LMX2
        AUX(:)=AUX(:)+SHIFTEDPHI(:,LM)*HPHI2(:,LM)
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,H)
!
!     ==TEST=========================================================
      AUX(:)=0.D0
      DO LM=1,LMX2
        AUX(:)=AUX(:)+SHIFTEDPHI(:,LM)*TESTHKINPHI2(:,LM)
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,TESTHKIN)
!
!     ==TEST=========================================================
      AUX(:)=0.D0
      DO LM=1,LMX2
        AUX(:)=AUX(:)+SHIFTEDPHI(:,LM)*TESTHHARTREEPHI2(:,LM)
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,TESTHHARTREE)
!
!     ==TEST=========================================================
      AUX(:)=0.D0
      DO LM=1,LMX2
        AUX(:)=AUX(:)+SHIFTEDPHI(:,LM)*TESTHXCPHI2(:,LM)
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,TESTHXC)
!
!     == OVERLAP =========================================================
      AUX(:)=0.D0
      DO LM=1,LMX2
        AUX(:)=AUX(:)+SHIFTEDPHI(:,LM)*PHI2(:,LM)
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,S)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE HARTREEPAIR(GID,NR,DR21,LMRX1,AEZ1,RHO1,LMRX2,AEZ2,RHO2,E)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: GID
      INTEGER(4) ,INTENT(IN)    :: NR
      REAL(8)    ,INTENT(IN)    :: DR21(3)
      INTEGER(4) ,INTENT(IN)    :: LMRX1       ! X#(POTENTIAL ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)    :: LMRX2       ! X#(POTENTIAL ANGULAR MOMENTA)
      REAL(8)    ,INTENT(IN)    :: AEZ1
      REAL(8)    ,INTENT(IN)    :: AEZ2
      REAL(8)    ,INTENT(IN)    :: RHO1(NR,LMRX1)
      REAL(8)    ,INTENT(IN)    :: RHO2(NR,LMRX2)
      REAL(8)    ,INTENT(OUT)   :: E
      REAL(8)                   :: POT(NR,LMRX1)
      REAL(8)                   :: SHIFTEDPOT(NR,LMRX2)
      REAL(8)                   :: AUX(NR)
      INTEGER(4)                :: L,LM      
      INTEGER(4)                :: LX1      
      REAL(8)                   :: PI,Y0
      REAL(8)                   :: DIS
      REAL(8)                   :: R(NR)
      REAL(8)                   :: VAL
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  CALCULATE HARTREE POTENTIAL                                         ==
!     ==========================================================================
      DO LM=1,LMRX1
        L=SQRT(REAL(LM-1,KIND=8)+1.D-10)
        CALL RADIAL$POISSON(GID,NR,L,RHO1(:,LM),POT(:,LM))
      ENDDO
!     == ADD NUCLEAR CHARGE ====================================================
      POT(:,1)=POT(:,1)-AEZ1/R(:)/Y0
!
!     ==========================================================================
!     ==  SHIFT DENSITY RHO2 TO POSITION 1                                    ==
!     ==========================================================================
      CALL SPHERICAL$SHIFTCENTER(GID,NR,DR21,LMRX1,POT,LMRX2,SHIFTEDPOT)
!
!     ==========================================================================
!     ==  CALCULATE ENERGY                                                    ==
!     ==========================================================================
      AUX(:)=0.D0
      DO LM=1,LMRX2
        AUX(:)=AUX(:)+SHIFTEDPOT(:,LM)*RHO2(:,LM)
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,E)
      CALL RADIAL$VALUE(GID,NR,SHIFTEDPOT(:,1),0.D0,VAL)
      E=E-AEZ2*VAL*Y0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCF(GID,NR,LMRX,LMX,LMREPX,AEZ,DREL,VEMB,LMAX,VCPAULI,OCPAULI &
     &              ,NC,PHIC,TPHIC,NB,LOFI,NN,POT,RHO,NCG,NBG,EBG,FBG,PHI,TPHI)
!     **************************************************************************
!     **  THE DEFORMED ION FULFILLES THE SCHROEDINGER EQUATION WITH           **
!     **  POTENTIAL "POT+POTREP", WHERE                                       **
!     **  POT CONTAINS                                                        **
!     **     1) THE HARTREE POTENTIAL OF THIS ATOM                            **
!     **     2) THE XC POTENTIAL OF THIS AND THE NEIGHBORING ATOMS            **
!     **  POTREP CONTAINS:                                                    **
!     **     1) THE HARTREE POTENTIAL OF THE NEIGHBORING ATOMS                **
!     **     2) THE PAULI REPULSION OF THE NEIGHBORING ATOMS                  **
!     **  THE EMBEDDING POTENTIAL IS OBTAINED AS                              **
!     **     VEMB= POTREP
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)    :: GID
      INTEGER(4) ,INTENT(IN)    :: NR
      INTEGER(4) ,INTENT(IN)    :: LMRX       ! X#(POTENTIAL ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)    :: LMX        ! X#(WAVEF. ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)    :: LMREPX     ! X#(WAVEF. ANGULAR MOMENTA)
      REAL(8)    ,INTENT(IN)    :: AEZ
      REAL(8)    ,INTENT(IN)    :: DREL(NR)
      REAL(8)    ,INTENT(IN)    :: VEMB(NR,LMREPX)
      INTEGER(4) ,INTENT(IN)    :: LMAX
      REAL(8)    ,INTENT(IN)    :: VCPAULI(NR,LMAX)
      REAL(8)    ,INTENT(IN)    :: OCPAULI(NR,LMAX)
      REAL(8)    ,INTENT(INOUT) :: POT(NR,LMRX)  ! POTENTIAL
      REAL(8)    ,INTENT(OUT)   :: RHO(NR,LMRX)  ! ALL-ELECTRON DENSITY
      INTEGER(4) ,INTENT(IN)    :: NC            ! #(SHELLS)
      REAL(8)    ,INTENT(IN)    :: PHIC(NR,NC)   ! CORE WAVE FUNCTIONS        
      REAL(8)    ,INTENT(IN)    :: TPHIC(NR,NC)  ! KINETIC ENERGY OF CORE WF.
      INTEGER(4) ,INTENT(IN)    :: NB            ! #(SHELLS)
      INTEGER(4) ,INTENT(IN)    :: LOFI(NB)      ! ANGULAR MOMENTA
      INTEGER(4) ,INTENT(IN)    :: NN(NB)        ! #(NODES)
      INTEGER(4) ,INTENT(IN)    :: NCG       
      INTEGER(4) ,INTENT(IN)    :: NBG       
      REAL(8)    ,INTENT(OUT)   :: EBG(NBG)        ! ENERGY EIGENVALUES 
      REAL(8)    ,INTENT(OUT)   :: FBG(NBG)        ! OCCUPATIONS
      REAL(8)    ,INTENT(OUT)   :: PHI(NR,LMX,NBG) ! WAVE FUNCTIONS
      REAL(8)    ,INTENT(OUT)   :: TPHI(NR,LMX,NBG)! KINETIC ENERGY TIMES PHI
      REAL(8)                   :: AEPHI(NR,LMX,NBG)
      REAL(8)                   :: TAEPHI(NR,LMX,NBG)
      REAL(8)    ,ALLOCATABLE   :: POTIN(:,:)
      REAL(8)    ,ALLOCATABLE   :: RHOTOT(:,:)
      REAL(8)                   :: POTOUT(NR,LMRX)
      REAL(8)                   :: POTBIG(NR,LMREPX)
      REAL(8)                   :: R(NR)
      REAL(8)                   :: AUX(NR),AUX1(NR)
      LOGICAL(4)                :: TBROYDEN=.FALSE.
      INTEGER(4) ,PARAMETER     :: NITER=50
      REAL(8)    ,PARAMETER     :: TOL=1.D-6
      REAL(8)    ,PARAMETER     :: KBT=1.D-1
      INTEGER(4)                :: ITER,LM1,LM2,LM3,IBG,IB,L,LM,IR,M
      INTEGER(4)                :: LX
      INTEGER(4)                :: LMPOTINX
      CHARACTER(64)             :: STRING
      REAL(8)                   :: MAXDEV
      REAL(8)                   :: ETOT,EB,EH,EXC,EPOT
      REAL(8)                   :: EBGCORR(NBG)
      LOGICAL(4)                :: TCONV
      REAL(8)                   :: SVAR
      REAL(8)                   :: CG
      REAL(8)                   :: PI,Y0
      REAL(8)                   :: MU
!     **************************************************************************
      TBROYDEN=.FALSE.
      CALL RADIAL$R(GID,NR,R)
      LX=INT(SQRT(REAL(LMX-1)+1.D-8))
      IF(NBG.NE.SUM(2*LOFI(:)+1)) THEN
        CALL ERROR$STOP('SCF')
      END IF
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
!  
!     ==========================================================================        
!     == ALLOCATE ARRAYS                                                      ==
!     ==========================================================================        
      LMPOTINX=MAX(LMRX,LMREPX)
      ALLOCATE(POTIN(NR,LMPOTINX))
      ALLOCATE(RHOTOT(NR,LMPOTINX))
!  
!     ==========================================================================        
!     == DETERMINE OCCUPATIONS                                                ==
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
!     == START SELF-CONSISTENCY LOOP                                          ==
!     ==========================================================================        
      DO ITER=1,NITER
!  
!       ========================================================================        
!       == NONSCF CALCULATION                                                 ==
!       ========================================================================        
        POTIN(:,:)=VEMB(:,:)
        POTIN(:,:LMRX)=POTIN(:,:LMRX)+POT(:,:) 
!CALL WRITEPHI('POTIN',GID,NR,LMRX,POTIN)
       CALL NONSCF(GID,NR,LMPOTINX,POTIN,DREL,LMAX,VCPAULI,OCPAULI &
       &                               ,NC,NB,LOFI,NN,NCG,NBG,EBG,LX,PHI,TPHI)
!!$PRINT*,'EBG ',EBG
!!$PRINT*,'NCG',NCG
!!$DO IBG=NCG+1,NBG
!!$WRITE(STRING,*)IBG
!!$STRING='NODELESS_PHI'//TRIM(ADJUSTL(STRING))//'.DAT'
!!$CALL WRITEPHI(STRING,GID,NR,LMX,PHI(:,:,IBG))
!!$ENDDO
!  
!       ========================================================================        
!       == FILL IN CORE STATES                                                ==
!       ========================================================================        
        AEPHI=PHI
        TAEPHI=PHI
        CALL SCF_COREORTHO(GID,NR,NC,LOFI,PHIC,TPHIC,NBG,NCG,LMX,AEPHI,TAEPHI)
!!$
!!$
!!$
!!$        IBG=0
!!$        DO IB=1,NC
!!$          L=LOFI(IB)
!!$          DO M=1,2*L+1
!!$            IBG=IBG+1
!!$            PHI(:,:,IBG)=0.D0
!!$            TPHI(:,:,IBG)=0.D0
!!$            LM=L**2+M
!!$            PHI(:,LM,IBG)=PHIC(:,IB)
!!$            TPHI(:,LM,IBG)=TPHIC(:,IB)
!!$          ENDDO
!!$        ENDDO
!!$!  
!!$!       ========================================================================        
!!$!       == ORTHOGONALIZE VALENCE STATES TO CORE                               ==
!!$!       ========================================================================        
!!$        DO IBG=NCG+1,NBG
!!$          DO IB=NC,1,-1
!!$            L=LOFI(IB)
!!$            DO M=1,2*L+1
!!$              LM=L**2+M
!!$              AUX(:)=R(:)**2*PHI(:,LM,IBG)*PHIC(:,IB)
!!$              CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$              PHI(:,LM,IBG)=PHI(:,LM,IBG)-PHIC(:,IB)*SVAR
!!$              TPHI(:,LM,IBG)=TPHI(:,LM,IBG)-TPHIC(:,IB)*SVAR
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$DO IBG=1,NBG
!!$WRITE(STRING,*)IBG
!!$STRING='PHI'//TRIM(ADJUSTL(STRING))//'.DAT'
!!$CALL WRITEPHI(STRING,GID,NR,LMX,PHI(:,:,IBG))
!!$ENDDO
!!$!  
!!$!       ========================================================================        
!!$!       == NORMALIZE STATES TO ACCOUNT FOR POOR NUMERICS                      ==
!!$!       ========================================================================        
!!$        DO IBG=1,NBG
!!$          AUX(:)=0.D0
!!$          DO LM=1,LMX
!!$            AUX(:)=AUX(:)+PHI(:,LM,IBG)**2
!!$          ENDDO
!!$          AUX(:)=AUX(:)*R(:)**2
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$          SVAR=1.D0/SQRT(SVAR)
!!$          PHI(:,:,IBG)=PHI(:,:,IBG)*SVAR
!!$          TPHI(:,:,IBG)=TPHI(:,:,IBG)*SVAR
!!$        ENDDO
!  
!       ========================================================================        
!       == DETERMINE OCCUPATIONS                                              ==
!       ========================================================================        
        DO IBG=1,NBG
          AUX(:)=0.D0
          DO LM1=1,LMX
            DO LM2=1,LMX
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.EQ.0.D0) CYCLE
                AUX(:)=AUX(:)+CG*VEMB(:,LM3)*PHI(:,LM1,IBG)*PHI(:,LM2,IBG)
              ENDDO
            ENDDO
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          EBGCORR(IBG)=EBG(IBG)-SVAR
        ENDDO
!PRINT*,'EBGCORR ',EBGCORR
        CALL SCF_OPTF(NBG,AEZ,KBT,EBGCORR,FBG,MU)
!PRINT*,'FBG ',FBG
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
                IF(CG.EQ.0.D0) CYCLE
                RHO(:,LM3)=RHO(:,LM3)+FBG(IBG)*CG*AEPHI(:,LM1,IBG)*AEPHI(:,LM2,IBG)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!!$! TEST TOTAL CHARGE
!!$AUX(:)=R(:)**2*RHO(:,1)*Y0*4.D0*PI
!!$CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$SVAR=SVAR-AEZ
!!$PRINT*,'TOTAL CHARGE ',SVAR
!CALL WRITEPHI('RHO.DAT',GID,NR,LMRX,RHO)
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
!CALL WRITEPHI('POTOUT.DAT',GID,NR,LMRX,POTOUT)
!  
!       ========================================================================        
!       == PRINT WAVE FUNCTIONS                                               ==
!       ========================================================================        
        WRITE(*,FMT='("ITER",I5," EBG: ",20F10.5)')ITER,EBG,FBG
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
!CALL WRITEPHI('NEWPOT.DAT',GID,NR,LMRX,POT)
!STOP 'KIKI'
      ENDDO      
      IF(TBROYDEN)CALL BROYDEN$CLEAR()
      IF(.NOT.TCONV) THEN
        CALL ERROR$MSG('LOOP NOT CONVERGED')
        CALL ERROR$STOP('SCF')
      END IF
!
!     ==========================================================================
!     == WRAP UP                                                              ==
!     ==========================================================================
      RETURN
      END SUBROUTINE SCF
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCF_COREORTHO(GID,NR,NC,LOFI,PHIC,TPHIC,NBG,NCG,LMX,PHI,TPHI)
!     ****************************************************************************
!     **  COMPLETE MISSING CORE STATES                                          **
!     **  ORTHOGONALIZE TO THE CORE STATES                                      **
!     **  NORMALIZE WAVE FUNCTIONS                                              **
!     ****************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NC
      INTEGER(4),INTENT(IN) :: LOFI(NC)
      REAL(8)   ,INTENT(IN) :: PHIC(NR,NC)
      REAL(8)   ,INTENT(IN) :: TPHIC(NR,NC)
      INTEGER(4),INTENT(IN) :: NBG
      INTEGER(4),INTENT(IN) :: NCG
      INTEGER(4),INTENT(IN) :: LMX
      REAL(8)   ,INTENT(INOUT) :: PHI(NR,LMX,NBG)
      REAL(8)   ,INTENT(INOUT) :: TPHI(NR,LMX,NBG)
      INTEGER(4)            :: IBG,IB,L,M,LM
      REAL(8)               :: AUX(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: SVAR
!     ****************************************************************************
      CALL RADIAL$R(GID,NR,R)
!  
!     ========================================================================        
!     == FILL IN CORE STATES                                                ==
!     ========================================================================        
      IBG=0
      DO IB=1,NC
        L=LOFI(IB)
        DO M=1,2*L+1
          IBG=IBG+1
          PHI(:,:,IBG)=0.D0
          TPHI(:,:,IBG)=0.D0
          LM=L**2+M
          PHI(:,LM,IBG)=PHIC(:,IB)
          TPHI(:,LM,IBG)=TPHIC(:,IB)
        ENDDO
      ENDDO
!  
!     ========================================================================        
!     == ORTHOGONALIZE VALENCE STATES TO CORE                               ==
!     ========================================================================        
      DO IBG=NCG+1,NBG
        DO IB=NC,1,-1
          L=LOFI(IB)
          DO M=1,2*L+1
            LM=L**2+M
            AUX(:)=R(:)**2*PHI(:,LM,IBG)*PHIC(:,IB)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            PHI(:,LM,IBG)=PHI(:,LM,IBG)-PHIC(:,IB)*SVAR
            TPHI(:,LM,IBG)=TPHI(:,LM,IBG)-TPHIC(:,IB)*SVAR
          ENDDO
        ENDDO
      ENDDO
!  
!     ========================================================================        
!     == NORMALIZE STATES TO ACCOUNT FOR POOR NUMERICS                      ==
!     ========================================================================        
      DO IBG=1,NBG
        AUX(:)=0.D0
        DO LM=1,LMX
          AUX(:)=AUX(:)+PHI(:,LM,IBG)**2
        ENDDO
        AUX(:)=AUX(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        SVAR=1.D0/SQRT(SVAR)
        PHI(:,:,IBG)=PHI(:,:,IBG)*SVAR
        TPHI(:,:,IBG)=TPHI(:,:,IBG)*SVAR
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DEFORM$ETOT(GID,NR,LMRX,LMX,LMREPX,NBG,AEZ &
     &                         ,FBG,EBG,PHI,TPHI,ETOT,EB,EKIN,EH,EXC)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: LMX
      INTEGER(4),INTENT(IN) :: LMREPX
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
      INTEGER(4)            :: LM,LM1,LM2,LM3,IBG,IR
      INTEGER(4)            :: L
      REAL(8)               :: PI,Y0,SVAR
!     ****************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!
!     ============================================================================
!     == BANDSTRUCTURE ENERGY: EB
!     ============================================================================
!PRINT*,'FBG ',FBG
!PRINT*,'EBG ',EBG
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
            IF(CG.EQ.0.D0) CYCLE
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
!CALL RADIAL$INTEGRAL(GID,NR,AUX*R**2,SVAR)
!PRINT*,'IBG ',IBG,SVAR,FBG(IBG)
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,EKIN)
!
!     ============================================================================
!     == EXCHANGE AND CORRELATION ENERGY: EXC
!     ============================================================================
      CALL AUGMENTATION_XC(GID,NR,LMRX,1,RHO1,EXC,POTOUT,AUX)
!     == CUT OUT LOW DENSITY REGION ==========================================
      DO IR=1,NR
        IF(RHO1(IR,1)*Y0.LT.1.D-6) POTOUT(IR,:)=0.D0
      ENDDO
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
      ETOT=EKIN+EH+EXC
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NONSCF(GID,NR,LMRX,POT,DREL,LMAX,VCPAULI,OCPAULI &
     &                 ,NC,NB,LOFI,NNOFI,NCG,NBG,EBG,LX,PHI,TPHI)
!     **************************************************************************
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID 
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      REAL(8)   ,INTENT(IN) :: POT(NR,LMRX)
      REAL(8)   ,INTENT(IN) :: DREL(NR)
      INTEGER(4),INTENT(IN) :: LMAX
      REAL(8)   ,INTENT(IN) :: VCPAULI(NR,LMAX)
      REAL(8)   ,INTENT(IN) :: OCPAULI(NR,LMAX)
      INTEGER(4),INTENT(IN) :: NC
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      INTEGER(4),INTENT(IN) :: NNOFI(NB)
      INTEGER(4),INTENT(IN) :: NBG
      INTEGER(4),INTENT(IN) :: NCG
      INTEGER(4),INTENT(IN) :: LX
      REAL(8)   ,INTENT(OUT):: EBG(NBG)
      REAL(8)   ,INTENT(OUT):: PHI(NR,(LX+1)**2,NBG)
      REAL(8)   ,INTENT(OUT):: TPHI(NR,(LX+1)**2,NBG)
      REAL(8)   ,PARAMETER  :: RAD=5.D0      
      INTEGER(4),PARAMETER  :: NITER=5
      REAL(8)               :: R(NR)
      REAL(8)               :: G(NR,(LX+1)**2)
      REAL(8)               :: POT1(NR)
      REAL(8)               :: OV1(NR)
      INTEGER(4)            :: LMX
      REAL(8)               :: EOFI(NB)
      REAL(8)               :: ENU
      REAL(8)               :: VAL
      LOGICAL(4)            :: TMAINSH(LX+1)
      INTEGER(4)            :: NPHI
      REAL(8)               :: PHI0(NR),AUX(NR)
      LOGICAL(4)            :: TOK
      INTEGER(4)            :: IB,ITER,L,I1,I2,IBG,LM,M,I,J,J1,J2,IL,IBG1,IBG2
      REAL(8)               :: EMIN1,EMAX1,EMIN2,EMAX2
      INTEGER(4)            :: NCOLOR
      INTEGER(4)            :: ICOLOR(NBG)
      INTEGER(4)            :: LOFBG(NBG)
      REAL(8)               :: SPACING=0.5D0
      REAL(8)               :: ESPACE=0.1D0
      REAL(8)               :: ENUMAX=-0.3D0
      REAL(8)               :: SVAR
      CHARACTER(64)         :: STRING
      REAL(8)  ,ALLOCATABLE :: EBG1(:)
      REAL(8)  ,ALLOCATABLE :: PHI1(:,:,:)
      REAL(8)  ,ALLOCATABLE :: OPHI1(:,:,:)
      REAL(8)  ,ALLOCATABLE :: TPHI1(:,:,:)
      REAL(8)  ,ALLOCATABLE :: OVER(:,:)
      LOGICAL(4),ALLOCATABLE :: TMAP(:)
      LOGICAL(4)            :: TDO
      REAL(8)               :: PI,Y0
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      IF(NBG.NE.SUM(2*LOFI(:)+1)) THEN 
        CALL ERROR$MSG('NBG AND LOFI ARE INCONSISTENT')
        CALL ERROR$STOP('NONSCF')
      END IF
      EBG(:)=0.D0
      PHI(:,:,:)=0.D0
      TPHI(:,:,:)=0.D0
      LMX=(LX+1)**2
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == DETERMINE STARTING ENERGIES FROM SPHERICAL POTENTIAL                 ==
!     ==========================================================================
      EOFI(:)=0.D0
      DO IB=NC+1,NB
        G(:,:)=0.D0
        POT1(:)=POT(:,1)
        OV1(:)=0.D0
        IF(LOFI(IB)+1.LE.LMAX) THEN
          POT1(:)=POT1(:)+VCPAULI(:,LOFI(IB)+1)
          OV1(:)=OV1(:)+OCPAULI(:,LOFI(IB)+1)
        END IF
        CALL BOUNDSTATE_OV(GID,NR,LOFI(IB),0,RAD,DREL,G,NNOFI(IB),POT1(:),OV1(:),EOFI(IB),PHI0)
        IF(IB.LT.NB)EOFI(IB+1)=EOFI(IB)
      ENDDO
      EOFI(1:NC)=-1.D+10
!
!     ==========================================================================
!     == DETERMINE STARTING ENERGIES FROM SPHERICAL POTENTIAL                 ==
!     ==========================================================================
      IBG=0
      DO IB=1,NB
        EBG(IBG+1:IBG+2*LOFI(IB)+1)=EOFI(IB)
        IBG=IBG+2*LOFI(IB)+1
      ENDDO
!!$PRINT*,'====='
!!$PRINT*,'====='
!!$PRINT*,'====='
!!$PRINT*,'====='
!!$PRINT*,'EOFI ',EOFI
!!$PRINT*,'EBG ',EBG
      DO ITER=1,NITER
!       == ORDER ENERGIES ======================================================
        DO IBG1=1,NBG
          DO IBG2=IBG1+1,NBG
            IF(EBG(IBG2).LT.EBG(IBG1)) THEN
              SVAR=EBG(IBG1)
              EBG(IBG1)=EBG(IBG2)
              EBG(IBG2)=SVAR
            END IF
          ENDDO
        ENDDO
!
!       == START COLORS ========================================================
!       == COLORS IDENTIFY GROUPS OF STATES ====================================
        NCOLOR=0
        ICOLOR(:)=0
        IBG=0
        DO IB=1,NB
          L=LOFI(IB)
          IF(IB.GT.NC)NCOLOR=NCOLOR+1
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
!       ========================================================================
        IBG=NCG
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
          ALLOCATE(OPHI1(NR,LMX,NPHI))
          ALLOCATE(TPHI1(NR,LMX,NPHI))
          ALLOCATE(OVER(NPHI,NPHI))
          ALLOCATE(TMAP(NPHI))
          I1=IBG+1
          DO I2=IBG+1,IBG+NPHI
            TDO=I2.EQ.IBG+NPHI
            IF(I2.LT.IBG+NPHI) TDO=TDO.OR.ABS(EBG(I2+1)-EBG(I2)).GT.ESPACE
            IF(TDO) THEN
              ENU=SUM(EBG(I1:I2))/REAL(I2-I1+1)
              ENU=MIN(ENU,ENUMAX)
              G(:,:)=0.D0
!              CALL RADIAL$NONSPHBOUND_NONSO(GID,NR,LX,LMX,LMRX,TMAINSH,POT,DREL,G,ENU &
!      &                  ,NPHI,EBG1,PHI1,TPHI1,TOK)
              CALL RADIAL$NONSPHBOUND_NSO_SLOC_OV(GID,NR,LX,LMX,LMRX,TMAINSH &
       &               ,POT,LMAX,VCPAULI*Y0,OCPAULI*Y0,G,ENU,NPHI,EBG1,PHI1,TPHI1,TOK)
!
!             == DETERMINE OVERLAP OF PREVOUSLY CALCULATED STATES WITH THIS SET= 
              OPHI1=PHI1
              DO IL=1,LMAX
                L=IL-1
                DO M=1,2*L+1
                  LM=L**2+M
                  DO J1=1,NPHI
                    OPHI1(:,LM,J1)=OPHI1(:,LM,J1)*(1.D0+OCPAULI(:,IL))
                  ENDDO
                ENDDO
              ENDDO
!
              OVER(:,:)=0.D0
              DO J1=IBG+1,I1-1
                DO J2=1,NPHI
                  AUX(:)=0.D0
                  DO LM=1,LMX
                    AUX(:)=AUX(:)+PHI(:,LM,J1)*OPHI1(:,LM,J2)
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
!!$PRINT*,'THIS ICOLOR ',I,'I1,I2 ',I1,I2
!!$PRINT*,'EBG1',EBG1
!!$PRINT*,'TMAP',TMAP
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
          DEALLOCATE(OPHI1)
          DEALLOCATE(OVER)
          DEALLOCATE(TMAP)
          IBG=IBG+NPHI
        ENDDO
!!$PRINT*,'ICOLOR ',ICOLOR
!!$PRINT*,'EBG ',EBG
!!$PRINT*,'============================================================='
      ENDDO
!
!     ==========================================================================        
!     == CHECK ORTHOGONALITY                                                  ==        
!     ==========================================================================        
      DO I1=NCG+1,NBG
        DO I2=I1,NBG
          AUX(:)=0.D0
          DO LM=1,LMX
            AUX(:)=AUX(:)+PHI(:,LM,I1)*PHI(:,LM,I2)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
          IF(I1.EQ.I2)VAL=VAL-1.D0
!IF(ABS(VAL).GT.1.D-3)PRINT*,'0VERLAP ',I1,I2,VAL
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONEATOM$NEW(IAT,AEZ)
!     **************************************************************************
!     ** INITIALIZES AN ATOM                                                  **
!     **************************************************************************
      USE ONEATOM_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      REAL(8)   ,INTENT(IN)  :: AEZ          ! ATOMIC NUMBER
      REAL(8)   ,PARAMETER   :: R1=1.056D-4  ! INNERMOST POINT OF RADIAL GRID
      REAL(8)   ,PARAMETER   :: DEX=0.05D0   ! LOG SPACING OF RADIAL GRID
      INTEGER(4),PARAMETER   :: NR=250       ! #(RADIAL GRID POINTS)
      LOGICAL(4),PARAMETER   :: TFROZENCORE=.TRUE.
      INTEGER(4)             :: GID          ! GRID ID
      INTEGER(4)             :: NS(10),NP(10),ND(10),NF(10)
      INTEGER(4)             :: SOFI(19)
      INTEGER(4)             :: LOFI(19)
      INTEGER(4)             :: NNOFI(19)
      REAL(8)                :: FOFI(19)
      REAL(8)                :: EOFI(19)
      REAL(8)                :: RBOX
      REAL(8)                :: EV
      REAL(8)                :: SVAR
      INTEGER(4)             :: NB       ! #(STATES IN OCCUPIED SHELLS)
      INTEGER(4)             :: NC       ! #(CORE STATES)
      REAL(8)   ,ALLOCATABLE :: AEPOT(:)
      REAL(8)   ,ALLOCATABLE :: DREL(:)
      REAL(8)   ,ALLOCATABLE :: RHOADD(:)
      REAL(8)   ,ALLOCATABLE :: RHOCORE(:)
      REAL(8)   ,ALLOCATABLE :: RHOVALENCE(:)
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)   ,ALLOCATABLE :: PHI(:)
      REAL(8)   ,ALLOCATABLE :: PHIAT(:,:)
      REAL(8)   ,ALLOCATABLE :: TPHIAT(:,:)
      REAL(8)   ,ALLOCATABLE :: UAT(:,:)
      REAL(8)   ,ALLOCATABLE :: TUAT(:,:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)                :: ETOT,DEDQ,DEDRAD,EKINC
      REAL(8)                :: EBAND,EPOT,EKIN,EXC,EH
      REAL(8)   ,PARAMETER   :: Q=0.D0
      REAL(8)                :: RAD
      REAL(8)                :: PI,C0LL,Y0
      INTEGER(4)             :: L,IB,ISVAR,IC,IB2
      INTEGER(4)             :: NFILO
!     **************************************************************************
PRINT*,'ONEATOM_NEW ..............................',IAT
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL CONSTANTS('EV',EV)
      IF(IAT.GT.NATX) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$STOP('ONEATOM$NEW')
      END IF
      CALL FILEHANDLER$UNIT('PROT',NFILO)

      THIS=>THISARR(IAT)
      THIS%AEZ=AEZ
!
!     ==========================================================================
!     ==  DEFINE RADIAL GRID                                                  ==
!     ==========================================================================
      IF(GIDONEATOM.EQ.0) then
        CALL RADIAL$NEW('SHLOG',GIDONEATOM)
      END IF
      GID=GIDONEATOM
      THIS%GID=GID
      CALL RADIAL$SETR8(GID,'R1',R1)
      CALL RADIAL$SETR8(GID,'DEX',DEX)
      CALL RADIAL$SETI4(GID,'NR',NR)
      THIS%NR=NR
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      RBOX=R(NR-3)
      WRITE(NFILO,FMT='(72("="),T30,"  ATOM NR. ",I5,"  ")')IAT
      WRITE(NFILO,FMT='("ISOLATED ATOM IN A BOX WITH RADIUS ",F5.2," A.U")')RBOX
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
!     == SET POINTER TO LOWEST VALENCE STATE PER L ======================
      L=MAXVAL(LOFI(1:NB))
      ALLOCATE(THIS%IBBOND(L+1))
      THIS%IBBOND(:)=0
      DO IB=NC+1,NB
        L=LOFI(IB)
        IF(THIS%IBBOND(L+1).EQ.0) THIS%IBBOND(L+1)=IB
      ENDDO
PRINT*,'IBBOND',THIS%IBBOND
!
!     == SWITCH FROM FROZEN CORE TO ALL-ELECTRON ==============================
      IF(.NOT.TFROZENCORE)NC=0
!
!     == MAP ON THIS ==================================================
      THIS%NB=NB
      THIS%NC=NC
      ALLOCATE(THIS%LOFI(NB))
      ALLOCATE(THIS%NNOFI(NB))
      ALLOCATE(THIS%FOFI(NB))
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
      CALL AESCF(GID,NR,NB,TREL,LOFI,SOFI,FOFI,NNOFI,AEZ,RHOADD,RBOX,DREL,AEPOT,EOFI)
PRINT*,'EOFI FROM AESCF ',EOFI(1:NB)
!
!     ==========================================================================
!     == CONSTRUCT ATOMIC WAVE FUNCTIONS
!     ==========================================================================
      ALLOCATE(PHIAT(NR,NB))
      ALLOCATE(TPHIAT(NR,NB))
      ALLOCATE(AUX(NR))
      DO IB=1,NB
        AUX(:)=0.D0 !INHOMOGENITY
        CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RBOX,DREL,AUX,NNOFI(IB),AEPOT,EOFI(IB),PHIAT(:,IB))
        TPHIAT(:,IB)=-(AEPOT(:)*Y0-EOFI(IB))*PHIAT(:,IB)        
        PHIAT(NR-2:,IB)=0.D0
        AUX(:)=(R(:)*PHIAT(:,IB))**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        PHIAT(:,IB)=PHIAT(:,IB)/SQRT(SVAR)
      ENDDO
PRINT*,'EOFI FROM AE ',EOFI(1:NB)
!
!     ==========================================================================
!     == CONSTRUCT NODELESS ATOMIC WAVE FUNCTIONS
!     ==========================================================================
      ALLOCATE(UAT(NR,NB))
      DO IB=1,NB
        IC=0
        DO IB2=1,IB-1
          IF(LOFI(IB2).NE.LOFI(IB)) CYCLE
          IC=IB2
          EXIT
        ENDDO
        IF(IC.EQ.0) THEN
          UAT(:,IB)=PHIAT(:,IB)
        ELSE
          AUX(:)=UAT(:,IC) !INHOMOGENITY
          CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RBOX,DREL,AUX,0,AEPOT,EOFI(IB),UAT(:,IB))
          UAT(NR-2:,IB)=0.D0
        END IF
      ENDDO
!      CALL WRITEPHI('UAT1.DAT',GID,NR,NB,UAT)
PRINT*,'EOFI FROM NODELESS ',EOFI(1:NB)
      DEALLOCATE(AUX)
!
!     ==========================================================================
!     == CORE DENSITY
!     ==========================================================================
      ALLOCATE(RHOCORE(NR))
      RHOCORE(:)=0.D0
      DO IB=1,NC
        RHOCORE(:)=RHOCORE(:)+FOFI(IB)*C0LL*PHIAT(:,IB)**2
      ENDDO
!
!     ==========================================================================
!     == VALENCE DENSITY DENSITY
!     ==========================================================================
      ALLOCATE(RHOVALENCE(NR))
      RHOVALENCE(:)=0.D0
      DO IB=NC+1,NB
        RHOVALENCE(:)=RHOVALENCE(:)+FOFI(IB)*C0LL*PHIAT(:,IB)**2
      ENDDO
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
!     == DETERMINE TOTAL ENERGY                                               ==
!     ==========================================================================
      ALLOCATE(AUX(NR))
      CALL MYVOFRHO(GID,NR,AEZ,RHOCORE+RHOVALENCE,AUX,EH,EXC)
      DEALLOCATE(AUX)
      EBAND=SUM(EOFI*FOFI)
      CALL RADIAL$INTEGRAL(GID,NR,R(:)**2*AEPOT(:)*(RHOCORE(:)+RHOVALENCE(:)),EPOT)
      EKIN=EBAND-EPOT
      ETOT=EKIN+EH+EXC
!
!     ==========================================================================
!     == MAP DATA ONTO "THIS"                                                 ==
!     ==========================================================================
      THIS%EREF=ETOT
      THIS%EKINC=EKINC
      ALLOCATE(THIS%EOFI(NB))
      THIS%EOFI(1:NB)=EOFI(1:NB)
      ALLOCATE(THIS%ATOMPOT(NR))
      THIS%ATOMPOT(:)=AEPOT(:)
      DEALLOCATE(AEPOT)
      ALLOCATE(THIS%RHOCORE(NR))
      THIS%RHOCORE(:)=RHOCORE(:)
      DEALLOCATE(RHOCORE)
      ALLOCATE(THIS%DREL(NR))
      THIS%DREL=DREL
      DEALLOCATE(DREL)
      ALLOCATE(THIS%PHIAT(NR,NB))
      THIS%PHIAT(:,:)=PHIAT(:,:)
      DEALLOCATE(PHIAT)
      ALLOCATE(THIS%TPHIAT(NR,NB))
      THIS%TPHIAT(:,:)=TPHIAT(:,:)
      DEALLOCATE(TPHIAT)
      ALLOCATE(THIS%UAT(NR,NB))
      THIS%UAT(:,:)=UAT(:,:)
      DEALLOCATE(UAT)
!
!     ==========================================================================
!     == REPORT ================================================================
!     ==========================================================================
      WRITE(NFILO,FMT='("EIGENSTATES AFTER NON-SCF CALCULATION IN SCF POT")')
      DO IB=1,NB
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
     &                ," E[H]=",F20.9," E[EV]=",F20.9)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/EV
      ENDDO
      WRITE(NFILO,FMT='("AESCF  EBAND ",F10.5)')EBAND
      WRITE(NFILO,FMT='("AESCF  EPOT  ",F10.5)')EPOT
      WRITE(NFILO,FMT='("AESCF  EKIN  ",F10.5)')EKIN
      WRITE(NFILO,FMT='("AESCF  EH    ",F10.5)')EH
      WRITE(NFILO,FMT='("AESCF  EXC   ",F10.5)')EXC
      WRITE(NFILO,FMT='("AESCF  ETOT  ",F10.5)')ETOT
      WRITE(NFILO,FMT='(72("="))')
!
!     ==========================================================================
!     == DETERMINE CORE REPULSION POTENTIAL                                   ==
!     ==========================================================================
!      CALL ONEATOM_CORE(IAT)
!
!     ==========================================================================
!     ==  WRITE ANALYSIS FILES                                                ==
!     ==========================================================================
!!$      CALL WRITEPHI('PHIAT.DAT',GID,NR,THIS%NB,THIS%PHIAT)
!!$      CALL WRITEPHI('UAT.DAT',GID,NR,THIS%NB,THIS%UAT)
!!$      CALL WRITEPHI('ETAC.DAT',GID,NR,THIS%NLC,THIS%ETACPAULI)
!!$      CALL WRITEPHI('VCPAULI.DAT',GID,NR,THIS%NLC,THIS%VCPAULI)
!!$      CALL WRITEPHI('OCPAULI.DAT',GID,NR,THIS%NLC,THIS%OCPAULI)
PRINT*,'ONEATOM_NEW DONE..........................',IAT
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONEATOM_CORE(IAT)
!     **************************************************************************
!     **************************************************************************
      USE ONEATOM_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      REAL(8)               :: AEZ
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: NB
      INTEGER(4)            :: NC
      REAL(8)               :: RCOV
      REAL(8)   ,ALLOCATABLE:: R(:)
      INTEGER(4)            :: LMAX
      INTEGER(4),ALLOCATABLE:: LOFI(:)
      REAL(8)   ,ALLOCATABLE:: ETACPAULI(:)
      REAL(8)   ,ALLOCATABLE:: OCPAULI(:)
      REAL(8)   ,ALLOCATABLE:: VCPAULI(:)
      REAL(8)   ,ALLOCATABLE:: UN(:)
      REAL(8)   ,ALLOCATABLE:: UDOT(:)
      REAL(8)   ,ALLOCATABLE:: UDDOT(:)
      REAL(8)   ,ALLOCATABLE:: PHI(:)
      REAL(8)   ,ALLOCATABLE:: DREL(:)
      REAL(8)   ,ALLOCATABLE:: G(:)
      REAL(8)               :: E
      REAL(8)               :: X0,DX,XM,Z0,ZM
      INTEGER(4)            :: ISTART,IBI
      INTEGER(4)            :: NLC
      INTEGER(4),PARAMETER  :: NITER=100
      REAL(8)   ,PARAMETER  :: TOL=1.D-6
      INTEGER(4)            :: I,IL,L,IC,IB,IR
      REAL(8)               :: PI,Y0
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      THIS=>THISARR(IAT)
      AEZ=THIS%AEZ
      GID=THIS%GID
      NR=THIS%NR
      NB=THIS%NB
      NC=THIS%NC
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
      ALLOCATE(DREL(NR))
      DREL(:)=THIS%DREL(:)
      ALLOCATE(LOFI(NB))
      LOFI(:)=THIS%LOFI(:)
!
!     ==========================================================================
!     ==  DEFINE RADIAL GRID                                                  ==
!     ==========================================================================
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  ALLOCATE POINTERS                                                   ==
!     ==========================================================================
      NLC=MAXVAL(LOFI(1:NC))+1
      IF(NC.EQ.0) NLC=0
      THIS%NLC=NLC
      IF(.NOT.ASSOCIATED(THIS%ETACPAULI))ALLOCATE(THIS%ETACPAULI(NR,NLC))
      IF(.NOT.ASSOCIATED(THIS%VCPAULI))ALLOCATE(THIS%VCPAULI(NR,NLC))
      IF(.NOT.ASSOCIATED(THIS%OCPAULI))ALLOCATE(THIS%OCPAULI(NR,NLC))
      THIS%ETACPAULI=0.D0
      THIS%VCPAULI=0.D0
      THIS%OCPAULI=0.D0
!
!     ==========================================================================
!     ==  ALLOCATE WORK ARRAYS                                                ==
!     ==========================================================================
      ALLOCATE(ETACPAULI(NR))
      ALLOCATE(OCPAULI(NR))
      ALLOCATE(VCPAULI(NR))
      ALLOCATE(UDOT(NR))
      ALLOCATE(UDDOT(NR))
      ALLOCATE(UN(NR))
      ALLOCATE(PHI(NR))
      ALLOCATE(G(NR))
      DO IL=1,NLC
        L=IL-1
!
!       =========================================================================
!       == FIND ENERGY E OF BONDING STATE                                      ==
!       =========================================================================
        IB=THIS%IBBOND(L+1)
        do i=1,ib-1
          if(lofi(i).ne.l) cycle
          ic=i
        enddo
        if(ic.eq.0) then
!         == no core state -> no core repulsion.... =============================
          THIS%VCPAULI(:,IL)=0.d0
          THIS%OCPAULI(:,IL)=0.d0
          THIS%ETACPAULI(:,IL)=0.d0
          cycle
        end if
        UN(:)=THIS%Ubox(:,IC)
        IF(IB.EQ.0) THEN
          E=0.D0
          GOTO 1000
        END IF
        ISTART=1
        X0=THIS%EOFI(IB)
        DX=-1.D-3
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
        DO I=1,NITER
          E=X0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,THIS%ATOMPOT,THIS%DREL,0,UN,L,E,1,PHI)
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RCOV,Z0)
          Z0=Z0-0.5D0
          IF(ABS(2.D0*DX).LE.TOL) EXIT
          CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
        ENDDO
        IF(ABS(DX).GT.TOL) THEN
          CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
          CALL ERROR$MSG('BOUND STATE NOT FOUND')
          CALL ERROR$STOP('ONEATOM')
        END IF
PRINT*,'ENERGY OF THE BINDING STATE FOR L= ',L,E
1000 CONTINUE
!
!       ========================================================================
!       == DETERMINE CORE REPULSION POTENTIAL                                 ==
!       ========================================================================
        CALL SCHROEDINGER$SPHERICAL(GID,NR,THIS%ATOMPOT,THIS%DREL,0,UN,LOFI(IB) &
     &                                                        ,E,1,UDOT)
        CALL SCHROEDINGER$SPHERICAL(GID,NR,THIS%ATOMPOT,THIS%DREL,0,UDOT,LOFI(IB) &
     &                                                        ,E,1,UDDOT)
!       ==  AVOID DIVIDE BY ZERO IF UDOT=0 NEAR THE NUCLEUS
        I=1
        DO IR=1,NR
          IF(UDOT(IR).NE.0.D0) THEN
            I=IR
            EXIT
          END IF
        ENDDO
!       == NOW CONSTRUCT PAULI POTENTIAL
        ETACPAULI(I:)=-UN(I:)/UDOT(I:)/Y0
        ETACPAULI(1:I)=ETACPAULI(I)
!       ==  AVOID DIVIDE BY ZERO IF UDOT=0 NEAR THE NUCLEUS
        I=1
        DO IR=1,NR
          IF(UDDOT(IR).NE.0.D0) THEN
            I=IR
            EXIT
          END IF
        ENDDO
!       == CONSTRUCT OCPAULI=-DETA/DE ==========================================
        OCPAULI(I:)=-ETACPAULI(I:)*UDDOT(I:)/UDOT(I:)
        OCPAULI(1:I)=OCPAULI(I)
!       == CONSTRUCT VCPAULI=ETA+E*OCPAULI ====================================
        VCPAULI(:)=ETACPAULI(:)+E*OCPAULI(:)
        THIS%VCPAULI(:,IL)=VCPAULI
        THIS%OCPAULI(:,IL)=OCPAULI
        THIS%ETACPAULI(:,IL)=ETACPAULI
      ENDDO
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
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      REAL(8)   ,INTENT(IN) :: Q    ! charge in electron charges
      REAL(8)   ,INTENT(IN) :: RAD  ! radius of the spherical box
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,INTENT(OUT):: DEDQ
      REAL(8)   ,INTENT(OUT):: DEDRAD
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: NB
      INTEGER(4)            :: NC
      REAL(8)               :: AEZ
      REAL(8)               :: QVALENCE
      INTEGER(4),ALLOCATABLE:: LOFI(:)
      INTEGER(4),ALLOCATABLE:: SOFI(:)
      INTEGER(4),ALLOCATABLE:: NNOFI(:)
      REAL(8)   ,ALLOCATABLE:: EOFI(:)
      REAL(8)   ,ALLOCATABLE:: FOFI(:)
      REAL(8)   ,ALLOCATABLE:: AEPOT(:)
      REAL(8)   ,ALLOCATABLE:: POTIN(:)
      REAL(8)   ,ALLOCATABLE:: POTH(:)
      REAL(8)   ,ALLOCATABLE:: POTXC(:)
      REAL(8)   ,ALLOCATABLE:: RHOADD(:)
      REAL(8)   ,ALLOCATABLE:: RHO(:)
      REAL(8)   ,ALLOCATABLE:: DREL(:)
      REAL(8)   ,ALLOCATABLE:: PHI(:,:)
      REAL(8)   ,ALLOCATABLE:: TPHI(:,:)
      REAL(8)   ,ALLOCATABLE:: PHI1(:,:)
      REAL(8)   ,ALLOCATABLE:: UBOX(:,:)
      REAL(8)   ,ALLOCATABLE:: TUBOX(:,:)
      REAL(8)   ,ALLOCATABLE:: PHIDOT(:)
      REAL(8)   ,ALLOCATABLE:: G(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:),AUX1(:)
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)   ,ALLOCATABLE:: ETAPAULI(:,:)
      REAL(8)   ,ALLOCATABLE:: UDOT(:,:)
      REAL(8)   ,ALLOCATABLE:: UN(:,:)
      REAL(8)               :: SVAR,DER,VAL,DERDOT,VALDOT
      REAL(8)               :: PI,C0LL,Y0
      REAL(8)               :: EKIN,EH,EXC
      INTEGER(4)            :: IBI,ISTART
      REAL(8)               :: X0,Z0,DX,XM,ZM,E
      REAL(8)               :: XMAX
      LOGICAL(4)            :: CONVG
      INTEGER(4)            :: IB,IR,ITER,IRBOX,IL,L,LMAX,ncount
      INTEGER(4)            :: IC
      INTEGER(4),PARAMETER  :: NITER=100
      REAL(8)   ,PARAMETER  :: TOL=1.D-5
      REAL(8)   ,PARAMETER  :: KBT=1.D-3
      REAL(8)               :: TS,F
      REAL(8)               :: RCOV
      REAL(8)               :: EV
      INTEGER(4)            :: NFILO
INTEGER(4) :: I,J
REAL(8)    :: DEG,SVAR1
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      C0LL=1.D0/SQRT(4.D0*PI)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL CONSTANTS('EV',EV)
      THIS=>THISARR(IAT)
      GID=THIS%GID
      NR=THIS%NR
      AEZ=THIS%AEZ      
      NB=THIS%NB
      NC=THIS%NC
      ALLOCATE(LOFI(NB))
      ALLOCATE(NNOFI(NB))
      ALLOCATE(EOFI(NB))
      ALLOCATE(FOFI(NB))
      ALLOCATE(SOFI(NB))
      LOFI(:)=THIS%LOFI(:)
      NNOFI(:)=THIS%NNOFI(:)
      EOFI(:)=THIS%EOFI(:)
      FOFI(:)=THIS%FOFI(:)
      SOFI(:)=0.D0
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GE.RAD) EXIT
      ENDDO
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)      
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      WRITE(NFILO,FMT='(72("="),T30,"  ATOM NR. ",I5,"  ")')IAT
      WRITE(NFILO,FMT='("ISOLATED ATOM IN A BOX WITH RADIUS ",F5.2," A.U")')RAD
!
!     ==========================================================================
!     ==  OBTAIN SELF-CONSISTENT POTENTIAL                                    ==
!     ==========================================================================
      ALLOCATE(AEPOT(NR));       AEPOT(:)=THIS%ATOMPOT(:)
      ALLOCATE(POTIN(NR));       POTIN(:)=AEPOT(:)
      ALLOCATE(RHOADD(NR));      RHOADD=0.D0
      ALLOCATE(DREL(NR));        DREL=THIS%DREL
      ALLOCATE(RHO(NR))
      ALLOCATE(PHI(NR,NB));      PHI(:,:)=0.D0
      ALLOCATE(tPHI(NR,NB));     tPHI(:,:)=0.D0
      ALLOCATE(PHIDOT(NR))
      ALLOCATE(G(NR))
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
      ALLOCATE(POTH(NR))
      ALLOCATE(POTXC(NR))
      XMAX=0.D0
      CONVG=.FALSE.
      CALL BROYDEN$NEW(NR,4,1.D0)
      POTIN=AEPOT   ! start with atomic potential
      DO ITER=1,NITER
        DO IB=1,NB
          G(:)=0.D0
          CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RAD,DREL,G,NNOFI(IB),POTIN,EOFI(IB),PHI(:,IB))
!         == normalize wave functions ==========================================
          AUX(:)=(R(:)*PHI(:,IB))**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1(:),RAD,SVAR)
          PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
!         == kinetic energy times wave function ================================
          TPHI(:,IB)=-(potin(:)*Y0-EOFI(IB))*PHI(:,IB)        
        ENDDO
!
!       == DETERMINE OCCUPATIONS ===============================================
        CALL ONEATOM_OPTF(NB,aez+q,KBT,LOFI,EOFI,FOFI,DEDQ)
!
!       == DETERMINE DENSITY ===================================================
        RHO(:)=0.d0
        DO IB=1,NB
          RHO(:)=RHO(:)+FOFI(IB)*C0LL*PHI(:,IB)**2
        ENDDO
!
!       == DETERMINE OUTPUT POTENTIAL ==========================================
        CALL MYVOFRHOWITHRAD(GID,NR,RAD,AEZ,RHO,POTH,POTXC,EH,EXC)
        AEPOT=POTH+POTXC
!
!       == EXIT IF CONVERGED ===================================================
        IF(CONVG) THEN
          CALL BROYDEN$CLEAR
          EXIT
        END IF
!
!       == POTENTIAL MIXING ====================================================
!       == THE POTENTIAL IS SET EQUAL BEYOND THE BOX RADIUS, BECAUSE OTHERWISE
!       == THIS REGION DETERMINES THE CONVERGENCE ==============================
!        POTIN(IRBOX:)=AEPOT(IRBOX:)
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
!     == map wave functions onto this ==========================================
      IF(.NOT.ASSOCIATED(THIS%PHIBOX))ALLOCATE(THIS%PHIBOX(NR,NB))
      IF(.NOT.ASSOCIATED(THIS%TPHIBOX))ALLOCATE(THIS%TPHIBOX(NR,NB))
      THIS%PHIBOX(:,:)=PHI(:,:)
      THIS%TPHIBOX(:,:)=TPHI(:,:)
!
!     ==========================================================================
!     == CALCULATE NODELESS WAVE FUNCTIONS                                    ==
!     ==========================================================================
      ALLOCATE(UBOX(NR,NB))
      ALLOCATE(TUBOX(NR,NB))
      DO IB=1,NB
        IC=0
        DO I=1,IB-1
          IF(LOFI(I).NE.LOFI(IB))CYCLE
          IC=I
        ENDDO
        IF(IC.EQ.0) THEN
          UBOX(:,IB)=PHI(:,IB)
          TUBOX(:,IB)=TPHI(:,IB)
        ELSE 
          G=UBOX(:,IC)
          AUX(:)=0.D0
          CALL BOUNDSTATE(GID,NR,LOFI(IB),0,RAD,AUX,G,0,AEPOT,EOFI(IB),UBOX(:,IB))
          TUBOX(:,IB)=G(:)+(EOFI(IB)-AEPOT(:)*Y0)*UBOX(:,IB)        
        END IF
      ENDDO
!
!     == map wave functions onto this ==========================================
      IF(.NOT.ASSOCIATED(THIS%UBOX))ALLOCATE(THIS%UBOX(NR,NB))
      IF(.NOT.ASSOCIATED(THIS%TUBOX))ALLOCATE(THIS%TUBOX(NR,NB))
      THIS%UBOX(:,:)=UBOX(:,:)
      THIS%TUBOX(:,:)=TUBOX(:,:)
!
!     ==========================================================================
!     == CALCULATE PRESSURE                                                   ==
!     ==========================================================================
      DEDRAD=0.D0
      DO IB=1,NB
        IF(FOFI(IB).EQ.0.D0) CYCLE
        CALL RADIAL$DERIVATIVE(GID,NR,PHI(:,IB),RAD,DER)
        IF(abs(DER).lt.1.d-8) CYCLE  ! AVOID CORE WAVE FUNCTIONS
!!$        G(:)=0.D0
!!$        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,SOFI(IB),G,LOFI(IB),EOFI(IB),1,PHI(:,IB))
!!$! NORMALIZATION NOT REQUIRED
!!$ AUX=R (:)**2*PHI(:,IB)**2
!!$ CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$ CALL RADIAL$VALUE(GID,NR,AUX1(:),RAD,SVAR)
!!$ PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
!!$!
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
      AUX(:)=RHO(:)*AEPOT(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,SVAR)
      EKIN=SUM(EOFI*FOFI)-SVAR
      CALL MYVOFRHOWITHRAD(GID,NR,RAD,AEZ,RHO,AUX,AUX1,EH,EXC)
      ETOT=EKIN+EH+EXC-THIS%EREF
!     == ENTROPY TERM OF THE ELECTRONS
      TS=0.D0
      DO IB=1,NB
        SVAR=REAL(2*(2*LOFI(IB)+1),KIND=8)
        F=FOFI(IB)/SVAR
        IF(F.LT.1.D-50) CYCLE
        IF(1.D0-F.LT.1.D-50) CYCLE
        TS=TS+KBT*SVAR*(F*LOG(F)+(1.D0-F)*LOG(1.D0-F))
      ENDDO   
      ETOT=ETOT+TS   
!
!     ==========================================================================
!     == REPORT ================================================================
!     ==========================================================================
      WRITE(NFILO,FMT='("EIGENSTATES AFTER NON-SCF CALCULATION IN SCF POT")')
      DO IB=NC+1,NB
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
     &                ," E[H]=",F20.9," E[EV]=",F20.9)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/EV
      ENDDO
      WRITE(NFILO,FMT='("ONEATOM  ETOT-AT  ",F10.5)')ETOT
      WRITE(NFILO,FMT='("ONEATOM  ETOT     ",F10.5)')ETOT+THIS%EREF
      WRITE(NFILO,FMT='("ONEATOM  EKIN     ",F10.5)')EKIN
      WRITE(NFILO,FMT='("ONEATOM  EHARTREE ",F10.5)')EH
      WRITE(NFILO,FMT='("ONEATOM  EXC      ",F10.5)')EXC
      WRITE(NFILO,FMT='("ONEATOM  ETS      ",F10.5)')TS
      WRITE(NFILO,FMT='("ONEATOM  DE/DRAD  ",F10.5)')DEDRAD
      WRITE(NFILO,FMT='("ONEATOM  DE/DQ    ",F10.5)')DEDQ
      WRITE(NFILO,FMT='(72("="))')
!
!     ==========================================================================
!     == DETERMINE CORE REPULSION POTENTIAL                                   ==
!     ==========================================================================
      CALL ONEATOM_CORE(IAT)
!
!     ==========================================================================
!     == CALCULATE REPULSIVE POTENTIAL FOR EACH ANGULAR MOMENTUM              ==
!     ==========================================================================
PRINT*,'IN ONEATOM: CONSTRUCTING NODELESS WAVE FUNCTIONS .............'
      LMAX=MAXVAL(LOFI(1:NB))
      ALLOCATE(ETAPAULI(NR,LMAX+1))
      ALLOCATE(UN(NR,LMAX+1))
      ALLOCATE(UDOT(NR,LMAX+1))
      DO IL=1,LMAX+1
        L=IL-1
        IB=THIS%IBBOND(L+1)   ! band index of lowest valence state per l
        IC=0
        DO I=1,IB-1
          IF(LOFI(I).EQ.L)IC=I
        ENDDO
        IF(IC.GT.0) THEN
          UN(:,IL)=THIS%Ubox(:,IC)
        ELSE
          UN(:,IL)=0.D0
        END IF
PRINT*,'IB,IC ',IB,IC
!
!       ========================================================================
!       == DETERMINE VALENCE BOUND STATE IN THE BOX (ANTIBONDING)             ==
!       ========================================================================
        IF(IB.EQ.0) THEN
          E=0.D0
          GOTO 1000
        END IF
        ISTART=1
        X0=EOFI(IB)
        DX=-1.D-3
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
        DO I=1,NITER
          E=X0
          G(:)=UN(:,IL)
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,0,G,L,E,1,PHI)
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RCOV,Z0)
          Z0=Z0-0.5D0
          IF(ABS(2.D0*DX).LE.TOL) EXIT
          CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
        ENDDO
        IF(ABS(DX).GT.TOL) THEN
          CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
          CALL ERROR$MSG('BOUND STATE NOT FOUND')
          CALL ERROR$STOP('ONEATOM')
        END IF
1000 CONTINUE
PRINT*,'ENERGY OF THE BINDING STATE FOR L= ',L,E
!
!       ========================================================================
!       == DETERMINE VALENCE BOUND STATE IN THE BOX (ANTIBONDING)             ==
!       ========================================================================
        G=UN(:,IL)
        CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RAD,DREL,G,0 &
     &                                             ,POTIN,EOFI(IB),UN(:,IL))
PRINT*,'ENERGY OF VALENCE STATE ',L,EOFI(IB),E
!
!       ========================================================================
!       == DETERMINE VALENCE REPULSION POTENTIAL                              ==
!       ========================================================================
        G=UN(:,IL)
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
        ETAPAULI(I:,IL)=-UN(I:,IL)/UDOT(I:,IL)/Y0
        ETAPAULI(1:I,IL)=ETAPAULI(I,IL)
      ENDDO
!
!     ==========================================================================
!     == MAP DATA BACK TO THIS                                                ==
!     ==========================================================================
      IF(.NOT.ASSOCIATED(THIS%ETAPAULI))ALLOCATE(THIS%ETAPAULI(NR,LMAX+1))
      ETAPAULI(IRBOX:,:)=0.D0
      THIS%ETAPAULI(:,:)=ETAPAULI(:,:)
      IF(.NOT.ASSOCIATED(THIS%UDOT))ALLOCATE(THIS%UDOT(NR,LMAX+1))
      THIS%UDOT=UDOT
      IF(.NOT.ASSOCIATED(THIS%UN))ALLOCATE(THIS%UN(NR,LMAX+1))
      THIS%UN=UN
      IF(.NOT.ASSOCIATED(THIS%RHO))ALLOCATE(THIS%RHO(NR))
      THIS%RHO=0.D0
      THIS%RHO(1:IRBOX-1)=RHO(1:IRBOX-1)
      IF(.NOT.ASSOCIATED(THIS%AEPOT))ALLOCATE(THIS%AEPOT(NR))
      THIS%AEPOT(:)=0.D0
      THIS%AEPOT(:)=POTIN(1:IRBOX-1)
      IF(.NOT.ASSOCIATED(THIS%POTH))ALLOCATE(THIS%POTH(NR))
      THIS%POTH(:)=POTH(:)
      IF(.NOT.ASSOCIATED(THIS%POTXC))ALLOCATE(THIS%POTXC(NR))
      THIS%POTXC(:)=POTXC(:)

      THIS%EOFI(:)=EOFI(:)
      THIS%FOFI(:)=FOFI(:)
      THIS%RAD=RAD
      THIS%Q=Q
      THIS%DEDRAD=DEDRAD
      THIS%DEDQ=DEDQ
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
PRINT*,'TESTING ONEATOM.................................'
!      CALL WRITEPHI('PHIBEFORE.DAT',GID,NR,NB-NC,PHI(:,NC+1:NB))
      ALLOCATE(PHI1(NR,NB));PHI1(:,:)=0.D0
      DO IB=NC+1,NB
        L=THIS%LOFI(IB)
        E=THIS%EOFI(IB)
!       == PERFORM CALCULATION WITH PAULI POTENTIAL ============================
        POTIN=THIS%AEPOT
        G=0.D0
        IF(L+1.LE.THIS%NLC) THEN
          POTIN=POTIN+THIS%VCPAULI(:,L+1) -E*THIS%OCPAULI(:,L+1)
!          G=THIS%UAT(:,L+1)
        END IF
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,0,G,L,E,1,PHI(:,IB))
        AUX(:)=R(:)**2*PHI(:,IB)**2
!        IF(L+1.LE.THIS%NLC) AUX(:)=AUX(:)*(1.D0+THIS%OCPAULI(:,L+1))
        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
! CALL WRITEPHI('AUX',GID,NR,1,AUX)
! CALL WRITEPHI('PHI',GID,NR,1,PHI(:,IB))
        CALL RADIAL$VALUE(GID,NR,AUX1(:),RAD,SVAR)
! CALL WRITEPHI('AUX1',GID,NR,1,AUX1)
        PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
!
!       == PERFORM CALCULATION WITH INHOMOGENEITY   ============================
        POTIN=THIS%AEPOT
        G=0.D0
        IF(L+1.LE.THIS%NLC) THEN
          G=THIS%UAT(:,L+1)
        END IF
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,0,G,L,E,1,PHI1(:,IB))
        AUX(:)=R(:)**2*PHI1(:,IB)**2
        IF(L+1.LE.THIS%NLC) AUX(:)=AUX(:)*(1.D0+THIS%OCPAULI(:,L+1))
        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
        CALL RADIAL$VALUE(GID,NR,AUX1(:),RAD,SVAR)
        PHI1(:,IB)=PHI1(:,IB)/SQRT(SVAR)

      ENDDO
!
!     ==========================================================================
!     == WRAP UP                                                              ==
!     ==========================================================================
      DEALLOCATE(ETAPAULI)
      DEALLOCATE(UN)
      DEALLOCATE(UDOT)
      DEALLOCATE(POTH)
      DEALLOCATE(POTXC)
! call ONEATOM_WRITE(THIS)
PRINT*,'ONEATOM done.................................'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONEATOM_WRITE(THIS)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE ONEATOM_MODULE, ONLY : THISTYPE
      IMPLICIT NONE
      type(thistype),intent(in) :: this
      INTEGER(4) :: GID
      INTEGER(4) :: NR
      INTEGER(4) :: NLC
      INTEGER(4) :: NB
!     **************************************************************************
      GID=THIS%GID
      NR=THIS%NR
      NLC=THIS%NLC
      NB=THIS%NB
      IF(ASSOCIATED(THIS%ATOMPOT))CALL WRITEPHI('ATOMPOT',GID,NR,1,THIS%ATOMPOT)
      IF(ASSOCIATED(THIS%ETAPAULI))CALL WRITEPHI('ETAPAULI',GID,NR,NLC,THIS%ETAPAULI)
      IF(ASSOCIATED(THIS%ETACPAULI))CALL WRITEPHI('ETACPAULI',GID,NR,NLC,THIS%ETACPAULI)
      IF(ASSOCIATED(THIS%VCPAULI))CALL WRITEPHI('VCPAULI',GID,NR,NLC,THIS%VCPAULI)
      IF(ASSOCIATED(THIS%OCPAULI))CALL WRITEPHI('OCPAULI',GID,NR,NLC,THIS%OCPAULI)
      IF(ASSOCIATED(THIS%PHIAT))CALL WRITEPHI('PHIAT',GID,NR,NB,THIS%PHIAT)
      IF(ASSOCIATED(THIS%UAT))CALL WRITEPHI('UAT',GID,NR,NB,THIS%UAT)
      IF(ASSOCIATED(THIS%RHO))CALL WRITEPHI('RHO',GID,NR,1,THIS%RHO)
      IF(ASSOCIATED(THIS%AEPOT))CALL WRITEPHI('AEPOT',GID,NR,1,THIS%AEPOT)
      IF(ASSOCIATED(THIS%POTH))CALL WRITEPHI('POTH',GID,NR,1,THIS%POTH)
      IF(ASSOCIATED(THIS%POTXC))CALL WRITEPHI('POTXC',GID,NR,1,THIS%POTXC)
      IF(ASSOCIATED(THIS%UN))CALL WRITEPHI('UN',GID,NR,NLC,THIS%UN)
      IF(ASSOCIATED(THIS%UDOT))CALL WRITEPHI('UDOT',GID,NR,NLC,THIS%UDOT)
      IF(ASSOCIATED(THIS%UC))CALL WRITEPHI('UC',GID,NR,NLC,THIS%UC)
      IF(ASSOCIATED(THIS%UBOX))CALL WRITEPHI('UBOX',GID,NR,Nb,THIS%UBOX)
      IF(ASSOCIATED(THIS%PHIBOX))CALL WRITEPHI('PHIBOX',GID,NR,Nb,THIS%PHIBOX)
      STOP 'stop IN ONEATOM_WRITE'
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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SCF_OPTF(N,Q,KBT,E,F,MU)
!     **************************************************************************
!     **  DETERMINE OCCUPATIONS FOR A GIVEN CHARGE AND TEMPERATURE            **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N       ! ENERGY LEVELS
      REAL(8)   ,INTENT(IN) :: Q       ! TOTAL CHARGE ASSUMING 2E PER STATE
      REAL(8)   ,INTENT(IN) :: KBT     ! TEMPERATURE AS K_B*T
      REAL(8)   ,INTENT(IN) :: E(N)    ! ENERGY EIGENVALUE
      REAL(8)   ,INTENT(OUT):: F(N)    ! NUMBER OF ELECTRONS PER SHELL
      REAL(8)   ,INTENT(OUT):: MU      ! CHEMICAL POTENTIAL
      REAL(8)               :: MUMIN,MUMAX
      REAL(8)               :: DQ
      REAL(8)               :: X
      INTEGER(4),PARAMETER  :: NITER=1000
      REAL(8)   ,PARAMETER  :: QTOL=1.D-10
      INTEGER(4)            :: ITER,I
!     **************************************************************************
      MUMIN=MINVAL(E)-10.D0*KBT
      MUMAX=MAXVAL(E)+10.D0*KBT
      IF(Q.LT.0.D0.OR.Q.GT.2.D0*REAL(N)) THEN
        CALL ERROR$MSG('CHARGE OUT OF RANGE')
        CALL ERROR$R8VAL('Q',Q)
        CALL ERROR$R8VAL('QMIN',0.D0)
        CALL ERROR$I4VAL('QMAX',2*N)
        CALL ERROR$STOP('SCF_OPTF')
      END IF
      DO ITER=1,NITER
        MU=0.5D0*(MUMIN+MUMAX)
        DO I=1,N
          X=(E(I)-MU)/KBT
          IF(X.GT.50.D0) THEN
            F(I)=0.D0
          ELSE IF(X.LT.-50.D0) THEN
            F(I)=2.D0
          ELSE
            F(I)=2.D0/(1.D0+EXP(X))
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
        CALL ERROR$MSG('SCF_OPTF')
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
      REAL(8),   INTENT(OUT):: POTH(NR)
      REAL(8),   INTENT(OUT):: POTXC(NR)
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
        IF(R(IR).GE.RAD) EXIT
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
!!$      RHO1(:IRBOX)=RHO(:)
!!$      RHO1(IRBOX:)=0.D0     
!!$      ALPHA=LOG(2.D0)/RHOMIN
!!$      DO IR=1,NR
!!$        X=RHO1(IR)*Y0
!!$        IF(R(IR).GT.RAD) X=0.D0
!!$        IF(ALPHA*X.GT.30.D0) THEN
!!$          F=X
!!$          DF=1.D0
!!$        ELSE
!!$          YP=EXP(ALPHA*X)
!!$          YM=1.D0/YP
!!$          F=LOG(YP+YM)/ALPHA
!!$          DF=(YP-YM)/(YP+YM)
!!$        END IF
!!$        RHOPLUS(IR)=F/Y0
!!$        DRHOPLUSDRHO(IR)=DF/Y0
!!$      ENDDO
!!$
!!$      CALL RADIAL$DERIVE(GID,NR,RHOPLUS(:),GRHO)
!
      CALL RADIAL$DERIVE(GID,NR,RHO,GRHO)
!
      DO IR=1,NR
        RH=RHO(IR)*Y0
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
!!$      DO IR=1,NR
!!$        EVAL=EDEN(IR)/(4.D0*PI)
!!$        MUVAL=POTXC(IR)*Y0
!!$!       ==  N/F(N) =============================
!!$        NBYF=RHO1(IR)/RHOPLUS(IR)
!!$        DFDN=DRHOPLUSDRHO(IR)*Y0
!!$        SVAR=(1.D0-NBYF*DFDN)/RHOPLUS(IR)
!!$        EDEN(IR)=NBYF*EVAL
!!$        POTXC(IR)=NBYF*MUVAL+SVAR*EVAL
!!$        EDEN(IR)=EDEN(IR)*4.D0*PI
!!$        POTXC(IR)=POTXC(IR)/Y0
!!$      ENDDO
!
!     == EXCHANGE CORRELATION POTENTIAL IS SET TO ZERO OUTSIDE THE BOX. 
!     == OTHERWISE IT WILL HAVE A CONSTANT VALUE IN THE ZERO-DENSITY REGION
      POTXC(IRBOX:)=0.D0
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,EDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EXC)
!
!     ==========================================================================
!     ==  CUT OF POTENTIAL FOTR LOW DENSITIES                                 ==
!     ==========================================================================
      DO IR=1,NR
        IF(RHO(IR)*Y0.LT.1.D-6) POTXC(IR)=0.D0
      ENDDO
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
      REAL(8)                    :: SVAR
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
!PRINT*,'XMAX ',XMAX,XAV
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
      REAL(8)               :: POTXC(NR)
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
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,EDEN,EXC)
!
!     ==========================================================================
!     ==  CUT OF POTENTIAL FOR LOW DENSITIES                                  ==
!     ==========================================================================
      DO IR=1,NR
        IF(RHO(IR)*Y0.LT.1.D-6) POTXC(IR)=0.D0
      ENDDO
!
!     ==========================================================================
!     ==  NOW ADD MASSAGED XCHANGE POTENTIAL TO HARTREE POTENTIAL             ==
!     ==========================================================================
      POT=POT+POTXC
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
!        CALL WRITEPHI('EKIN.DAT',GID,NR,4,ARRAY)
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
      REAL(8)                    :: AUX(NR),AUX1(NR)
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
        IF(R(IR).GE.RBOX) EXIT
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
!     ==  NORMALIZE WAVE FUNCTION                                          ==
!     =======================================================================
      AUX(:)=R(:)**2*PHI(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
      PHI(:)=PHI(:)/SQRT(SVAR)
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
      SUBROUTINE BOUNDSTATE_OV(GID,NR,L,SO,RBOX,DREL,G,NN,POT,OV,E,PHI)
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
      REAL(8)    ,INTENT(IN)     :: OV(NR)  !OVERLAP OPERATOR
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
      REAL(8)                    :: AUX(NR),AUX1(NR)
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
        IF(R(IR).GE.RBOX) EXIT
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
        POT1(:)=POT(:)-E*OV(:)
        CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT1,E,IRCL,IROUT)
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
!     ==  NORMALIZE WAVE FUNCTION                                          ==
!     =======================================================================
      AUX(:)=R(:)**2*(1.D0+OV(:)*Y0)*PHI(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
      PHI(:)=PHI(:)/SQRT(SVAR)
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
      SUBROUTINE SPHERICAL$SHIFTCENTER(GID,NR,CENTERNEW,LMX1,FIN,LMX2,FOUT)
!     **                                                                      **
!     **  TRANSFORMS A FUNCTION GIVEN IN TRUE SPHERICAL HARMONICS             **
!     **  (ANGULAR MOMENTUM EIGENSTATES) TO A NEW CENTER                      **
!     **  THS CENTER IS DISPLACED IN Z-DIRECTION                              **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: LMX1
      INTEGER(4)  ,INTENT(IN) :: LMX2
      REAL(8)     ,INTENT(IN) :: CENTERNEW(3)
      REAL(8)     ,INTENT(IN) :: FIN(NR,LMX1)
      REAL(8)     ,INTENT(OUT):: FOUT(NR,LMX2)
      REAL(8)                 :: FINR(NR,LMX1)
      REAL(8)                 :: FOUTR(NR,LMX2)
      REAL(8)                 :: FINT(NR,LMX2)
      REAL(8)                 :: DIS
      REAL(8)                 :: R(NR)
      REAL(8)                 :: ROT(3,3)
      REAL(8)                 :: YLMROT1(LMX1,LMX1),YLMROT2(LMX2,LMX2)
      INTEGER(4)              :: LX1,LX2
      REAL(8)                 :: PLM1(LMX1),PLM2(LMX2)
      REAL(8)                 :: FONE(LMX1)
      LOGICAL(4)              :: TCHK
      REAL(8)                 :: COSALPHA,COSTHETA
      REAL(8)                 :: SINALPHA,SINTHETA
      REAL(8)                 :: DXDALPHA,DXDTHETA
      REAL(8)                 :: AUX(NR)
      INTEGER(4)              :: IR1,IR2
      REAL(8)                 :: P,X,DX             
      INTEGER(4)              :: LM1,LM2,L1,L2,M1,M2,L1START,ITHETA
      REAL(8)                 :: FAC,SVAR
      INTEGER(4) ,PARAMETER   :: NTHETA=300
      REAL(8)                 :: THETA,DTHETA
      REAL(8)                 :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      LX1=SQRT(REAL(LMX1-1,KIND=8)+1.D-10)
      LX2=SQRT(REAL(LMX2-1,KIND=8)+1.D-10)
      CALL RADIAL$R(GID,NR,R)
      DIS=SQRT(SUM(CENTERNEW(:)**2))
!
!     ==========================================================================
!     == ROTATE COORDINATE SYSTEM SUCH THAT THE NEW CENTER IS IN POSITIVE Z DIR=
!     ==========================================================================
      CALL SPHERICAL$ALIGN(CENTERNEW,ROT)
      CALL ROTATEYLM(LMX1,ROT,YLMROT1(:,:))

!     ==  FINR IS THE ROTATED  FUNCTION
      FINR(:,:)=0.D0
      DO LM1=1,LMX1    !LOOP OVER THE COMPONENTS OF THE ROTATED VECTOR
        DO LM2=1,LMX1 !THE ILM COMPONENT OF THE LMX VECTOR FOR A GIVEN R AFTER ROTATION
          SVAR=YLMROT1(LM1,LM2)
          IF(SVAR.EQ.0.D0) CYCLE
          FINR(:,LM1)=FINR(:,LM1)+FIN(:,LM2)*SVAR
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  TRANSFORM TO NEW CENTER                                             ==
!     ==========================================================================
      DTHETA=PI/REAL(NTHETA,KIND=8)
      FOUT(:,:)=0.D0
      DO IR2=1,NR
        P=R(IR2)
        DO ITHETA=1,NTHETA
          THETA=DTHETA*(-0.5D0+REAL(ITHETA,KIND=8))
          COSTHETA=COS(THETA)
          X=SQRT(DIS**2+P**2+2.D0*DIS*P*COSTHETA)
          COSALPHA=(DIS+P*COSTHETA)/X
          SINTHETA=SQRT(1.D0-COSTHETA**2)
!         == START RECURSION FOR ASSOCIATED LEGENDRE POLYNOMIALS WITH P_LL
          CALL PLGNDR(LMX1,LX1,COSALPHA,PLM1) 
          CALL PLGNDR(LMX2,LX2,COSTHETA,PLM2)
          DO LM1=1,LMX1
            CALL RADIAL$VALUE(GID,NR,FINR(:,LM1),X,FONE(LM1))
          ENDDO
          DO LM2=1,LMX2
            L2=SQRT(REAL(LM2-1,KIND=8)+1.D-10)
            M2=L2*(L2+1)-LM2+1
            M1=M2
            L1START=ABS(M2)
            SVAR=0.D0
            DO L1=L1START,LX1
              LM1=L1*(L1+1)-M2+1
              FAC=0.5D0*SQRT(REAL((2*L1+1)*(2*L2+1),KIND=8))
              SVAR=SVAR+FAC*FONE(LM1)*PLM1(LM1)*PLM2(LM2)
            ENDDO
            FOUT(IR2,LM2)=FOUT(IR2,LM2)+SVAR*SINTHETA*DTHETA
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ROTATE BACK TO ORIGINAL COORDINATE SYSTEM                           ==
!     ==========================================================================
      ROT=TRANSPOSE(ROT)
      CALL ROTATEYLM(LMX2,ROT,YLMROT2(:,:))
!     ==  FINR IS THE ROTATED  FUNCTION
      FINT(:,:)=FOUT(:,:)
      FOUT(:,:)=0.D0
      DO LM1=1,LMX2    !LOOP OVER THE COMPONENTS OF THE ROTATED VECTOR
        DO LM2=1,LMX2 !THE ILM COMPONENT OF THE LMX VECTOR FOR A GIVEN R AFTER ROTATION
          SVAR=YLMROT2(LM1,LM2)
          IF(SVAR.EQ.0.D0) CYCLE
          FOUT(:,LM1)=FOUT(:,LM1)+FINT(:,LM2)*SVAR
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPHERICAL$ALIGN(VEC,ROT)
!     **************************************************************************    
!     ** CONSTRUCT A ROTATION MATRIX ROT WHICH TRANSFORMS THE VECTOR VEC      **
!     ** ONTO THE POSITIVE Z-AXIS                                             **
!     **************************************************************************    
      IMPLICIT NONE
      REAL(8),INTENT(IN)         :: VEC(3)
      REAL(8),INTENT(OUT)        :: ROT(3,3)
      INTEGER(4)                 :: I,IVEC(1)
      REAL(8)                    :: DIS,DR1(3),DR2(3),DR3(3)
      REAL(8),PARAMETER          :: EX(3)=(/1.D0,0.D0,0.D0/)
      REAL(8),PARAMETER          :: EY(3)=(/0.D0,1.D0,0.D0/)
      REAL(8),PARAMETER          :: EZ(3)=(/0.D0,0.D0,1.D0/)
!     **************************************************************************    
      DIS=SQRT(SUM(VEC(:)**2))
      IF(DIS.EQ.0.D0) THEN
        ROT(:,:)=0.D0
        DO I=1,3
          ROT(I,I)=1.D0
        ENDDO
      END IF
      DR1(:)=VEC(:)/DIS
      IVEC=MINLOC(ABS(DR1))
      I=IVEC(1)
      IF(I.EQ.1) THEN
        DR2(1)=DR1(2)*EX(3)-DR1(3)*EX(2)
        DR2(2)=DR1(3)*EX(1)-DR1(1)*EX(3)
        DR2(3)=DR1(1)*EX(2)-DR1(2)*EX(1)
      ELSE IF(I.EQ.2) THEN
        DR2(1)=DR1(2)*EY(3)-DR1(3)*EY(2)
        DR2(2)=DR1(3)*EY(1)-DR1(1)*EY(3)
        DR2(3)=DR1(1)*EY(2)-DR1(2)*EY(1)
      ELSE IF(I.EQ.3) THEN
        DR2(1)=DR1(2)*EZ(3)-DR1(3)*EZ(2)
        DR2(2)=DR1(3)*EZ(1)-DR1(1)*EZ(3)
        DR2(3)=DR1(1)*EZ(2)-DR1(2)*EZ(1)
      END IF
      DR2(:)=DR2(:)/SQRT(SUM(DR2(:)**2))
      DR3(1)=DR1(2)*DR2(3)-DR1(3)*DR2(2)
      DR3(2)=DR1(3)*DR2(1)-DR1(1)*DR2(3)
      DR3(3)=DR1(1)*DR2(2)-DR1(2)*DR2(1)
!     ==  THE STRANGE MAPPING IS FOR CONSISTENCY REASONS
      ROT(3,:)=DR1(:)
      ROT(2,:)=-DR2(:)
      ROT(1,:)=DR3(:)
      RETURN
      END SUBROUTINE SPHERICAL$ALIGN
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PLOTPLGNDR()
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: N=200
      INTEGER(4),PARAMETER :: LMX=9
      INTEGER(4),PARAMETER :: NFIL=888
      INTEGER(4)           :: I
      INTEGER(4)           :: LX
      REAL(8)              :: PLM(LMX),X
      REAL(8)              :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      LX=SQRT(REAL(LMX-1,KIND=8)+1.D-10)
      OPEN(NFIL,FILE='PLGNDR.DAT')
      DO I=1,N
!        X=-1.D0+2.D0*REAL(I-1)/REAL(N-1)
!        CALL PLGNDR(LMX,LX,X,PLM)
 X=2.D0*PI*REAL(I-1)/REAL(N-1)
 CALL PLGNDR(LMX,LX,COS(X),PLM)
       WRITE(NFIL,FMT='(50F10.5)')X,PLM
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TESTTRANSFORM()
      IMPLICIT NONE
      INTEGER(4),PARAMETER  :: NR=250
      INTEGER(4),PARAMETER  :: LX2=12
      INTEGER(4),PARAMETER  :: NX=100
      INTEGER(4),PARAMETER  :: NFIL=192
      INTEGER(4)            :: GID
      REAL(8)               :: DR(3)=(/0.D0,2.D0,3.D0/)
      REAL(8)               :: R(NR)
      INTEGER(4)            :: LMX1
      INTEGER(4)            :: LMX2
      INTEGER(4)            :: L,I,LM2
      REAL(8)   ,ALLOCATABLE:: YLM(:)
      REAL(8)   ,ALLOCATABLE:: FIN(:,:)
      REAL(8)   ,ALLOCATABLE:: FOUT(:,:)
      REAL(8)               :: AUXI(NR,LX2+1)
      REAL(8)               :: DR2(3)
      REAL(8)               :: PI,Y0
      REAL(8)               :: SVAR,VAL
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      LMX1=1
      CALL RADIAL$NEW('SHLOG',GID)
      CALL RADIAL$SETI4(GID,'NR',NR)
      CALL RADIAL$SETR8(GID,'DEX',0.05D0)
      CALL RADIAL$SETR8(GID,'R1',1.056D-4)
!
!     ==========================================================================
!     == DECLARE RAW FUNCTION ==================================================
!     ==========================================================================
      CALL RADIAL$R(GID,NR,R)
      ALLOCATE(FIN(NR,LMX1))
      FIN(:,:)=0.D0
      FIN(:,1)=EXP(-0.25D0*R(:)**2)/Y0
!
!     ==========================================================================
!     == LOOP TO CHECK CONVERGENCE
!     ==========================================================================
      DO L=0,LX2
        LMX2=(L+1)**2
        ALLOCATE(FOUT(NR,LMX2))
!       ========================================================================
!       == PERFORM TRANSFORMATION                                             ==
!       ========================================================================
        CALL SPHERICAL$SHIFTCENTER(GID,NR,DR,LMX1,FIN,LMX2,FOUT)
!
!       ========================================================================
!       == RECONSTRUCT FUNCTION                                               ==
!       ========================================================================
        ALLOCATE(YLM(LMX2))
        DO I=1,NX
          DR2(:)=DR(:)*(-1.5D0+3.D0*REAL(I-1)/REAL(NX-1))
          CALL GETYLM(LMX2,DR2,YLM)
          SVAR=0.D0
          DO LM2=1,LMX2
            CALL RADIAL$VALUE(GID,NR,FOUT(:,LM2),SQRT(SUM(DR2**2)),VAL)
            SVAR=SVAR+VAL*YLM(LM2)
          ENDDO
          AUXI(I,L+1)=SVAR
        ENDDO
        DEALLOCATE(YLM)
        DEALLOCATE(FOUT)
      ENDDO
!
!     ==========================================================================
!     == WRITE RESULT TO FILE                                                 ==
!     ==========================================================================
      OPEN(NFIL,FILE='TESTTRANSFORM.DAT')
      DO I=1,NX
        SVAR=(-1.5D0+3.D0*REAL(I-1)/REAL(NX-1))
        WRITE(NFIL,FMT='(30F10.5)')SVAR,AUXI(I,:)
      ENDDO
      CLOSE(222)
      RETURN
      END

     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE COMPARETRANSFORM()
      IMPLICIT NONE
      INTEGER(4),PARAMETER  :: NR=250
      INTEGER(4),PARAMETER  :: LX1=2
      INTEGER(4),PARAMETER  :: LX2=3
      INTEGER(4)            :: GID
      INTEGER(4)            :: LMX1
      INTEGER(4)            :: LMX2
      REAL(8)               :: DR(3)=(/0.D0,2.D0,3.D0/)
      REAL(8)               :: R(NR)
      REAL(8)   ,ALLOCATABLE:: FIN(:,:),FOUT1(:,:),FOUT2(:,:)
      REAL(8)               :: PI,Y0
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      LMX1=(LX1+1)**2
      LMX2=(LX2+1)**2
      CALL RADIAL$NEW('SHLOG',GID)
      CALL RADIAL$SETI4(GID,'NR',NR)
      CALL RADIAL$SETR8(GID,'DEX',0.05D0)
      CALL RADIAL$SETR8(GID,'R1',1.056D-4)
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == DECLARE RAW FUNCTION ==================================================
!     ==========================================================================
      ALLOCATE(FIN(NR,LMX1))
      FIN(:,:)=0.D0
      FIN(:,1)=EXP(-0.1D0*R(:)**2)/Y0
      IF(LMX1.GE.2)FIN(:,2)=R(:)*EXP(-0.3D0*R(:)**2)/Y0
      IF(LMX1.GE.8)FIN(:,8)=R(:)**2*EXP(-0.5D0*R(:)**2)/Y0
!
!     ==========================================================================
!     == DECLARE RAW FUNCTION ==================================================
!     ==========================================================================
      ALLOCATE(FOUT1(NR,LMX2))
      ALLOCATE(FOUT2(NR,LMX2))
PRINT*,'PETER'
      CALL SPHERICAL$SHIFTCENTER(GID,NR,-DR,LMX1,FIN,LMX2,FOUT1)
!      CALL WRITEPHI('FOUT1.DAT',GID,NR,LMX2,FOUT1)
PRINT*,'ALEX'
      CALL SPHERICAL$TRANSFORM(GID,NR,DR,LMX1,LMX2,FIN,FOUT2)
!      CALL WRITEPHI('FOUT2.DAT',GID,NR,LMX2,FOUT2)
!      CALL WRITEPHI('FOUTDIFF.DAT',GID,NR,LMX2,FOUT2-FOUT1)
      DEALLOCATE(FIN)     
      DEALLOCATE(FOUT1)     
      DEALLOCATE(FOUT2)
      RETURN
      END
!......................................................................................
SUBROUTINE PBSCFWRITE(NR,LMRX,POT,POTREP)
INTEGER(4),INTENT(IN) :: NR
INTEGER(4),INTENT(IN) :: LMRX
REAL(8)   ,INTENT(IN) :: POT(NR,LMRX)
REAL(8)   ,INTENT(IN) :: POTREP(NR,LMRX)
INTEGER(4)            :: NFIL=888
INTEGER(4)            :: LM,IR
OPEN(NFIL,FILE='SCFIOFILE')
REWIND(NFIL)
WRITE(NFIL,*)NR,LMRX
DO LM=1,LMRX
DO IR=1,NR
WRITE(NFIL,*)POTREP(IR,LM),POT(IR,LM)
ENDDO
ENDDO
CLOSE(NFIL)
RETURN
END
SUBROUTINE PBSCFREAD(NR,LMRX,POT,POTREP)
INTEGER(4),INTENT(IN) :: NR
INTEGER(4),INTENT(IN) :: LMRX
REAL(8)   ,INTENT(OUT) :: POT(NR,LMRX)
REAL(8)   ,INTENT(OUT) :: POTREP(NR,LMRX)
INTEGER(4)            :: NFIL=888
INTEGER(4)            :: LM,IR
INTEGER(4)            :: LMRX1,NR1
OPEN(NFIL,FILE='PODDEY080313/SCFIOFILE_ALEX')
REWIND(NFIL)
READ(NFIL,*)NR1,LMRX1
IF(NR1.NE.NR.OR.LMRX1.NE.LMRX) STOP 'STOP IN PBSCFREAD'
DO LM=1,LMRX
DO IR=1,NR
READ(NFIL,*)POTREP(IR,LM),POT(IR,LM)
ENDDO
ENDDO
CLOSE(NFIL)
PRINT*,'DONE'
RETURN
END

!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LATTICE()
      USE ONEATOM_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NAT=2
      REAL(8)              :: AEZARR(NAT)
      REAL(8)              :: RBASFCC(3,3)
      REAL(8)              :: RBAS(3,3)
      REAL(8)              :: POS(3,NAT)
      REAL(8)              :: ALAT
      REAL(8)              :: RCOV(NAT)
      REAL(8)              :: Q(NAT)
      REAL(8)              :: RAD0(NAT)
      REAL(8)              :: RAD(NAT)
      REAL(8)              :: ETOT
      INTEGER(4)           :: IAT,ip,if
      integer(4),parameter :: np=30
      integer(4),parameter :: nf=1
      real(8)              :: x(np)
      real(8)              :: y(np,nf)
!     **************************************************************************
      AEZARR(:)=(/12.D0,8.D0/)
      ALAT=4.214D0/0.529177d0
      RBASFCC(:,1)=(/0.0D0,0.5D0,0.5D0/)
      RBASFCC(:,2)=(/0.5D0,0.0D0,0.5D0/)
      RBASFCC(:,3)=(/0.5D0,0.5D0,0.0D0/)
      POS(:,:)=0.D0
      POS(1,2)=0.5D0
      RBAS=RBASFCC*ALAT
      POS=POS*ALAT
      Q(:)=(/-2.d0,2.d0/)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      CALL DFT$SETI4('TYPE',2)
      CALL DFT$SETL4('SPIN',.FALSE.)
      DO IAT=1,NAT
        CALL ONEATOM$NEW(IAT,AEZARR(IAT))
        CALL PERIODICTABLE$GET(NINT(AEZARR(IAT)),'R(COV)',RCOV(IAT))
      ENDDO
      rad0(:)=rcov(:)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      open(unit=8,file='lattice_gnu.dat')
      do ip=1,np
        x(ip)=7.8d0+0.2d0*real(ip-1)
        do if=1,nf
          rad(2)=x(ip)
          rad(1)=2.0d0+0.2d0*real(if-1)
          CALL LATTICE$ENERGY(NAT,RBAS,POS,Q,RAD,ETOT)
          y(ip,if)=etot          
          write(8,FMT='(30f10.5)')rad(1),rad(2),etot
        enddo
      enddo
      close(8)
      open(unit=8,file='lattice.dat')
      do ip=1,np
        write(8,FMT='(30f10.5)')x(ip),y(ip,:)
      enddo
      close(8)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LATTICE$ENERGY(NAT,RBAS,POS,Q,RAD,ETOT)
      USE ONEATOM_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RBAS(3,3)
      REAL(8)   ,INTENT(IN) :: POS(3,NAT)
      REAL(8)   ,INTENT(IN) :: Q(NAT)
      REAL(8)   ,INTENT(IN) :: RAD(NAT)
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)               :: DEDQ(NAT),DEDQ1(NAT),DEDRAD(NAT),FORCE(3,NAT)
      INTEGER(4),parameter  :: ndiv=2
      INTEGER(4)            :: IAT,IAT1,IAT2,it1,it2,it3
      INTEGER(4)            :: NR
      INTEGER(4)            :: GID
      REAL(8)               :: DETOT
      REAL(8)               :: DIS,DR(3)
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:),AUX1(:),RHO(:),POT(:)
      REAL(8)               :: VAL,VAL1,VAL2
      REAL(8)               :: PI,y0
      LOGICAL(4),PARAMETER  :: TPR=.true.
      INTEGER(4)            :: nfilo
!     **************************************************************************
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      PI=4.D0*ATAN(1.D0)
      y0=1.d0/sqrt(4.d0*pi)
      GID=GIDONEATOM
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      DO IAT=2,NAT
        IF(GID.NE.THISARR(IAT)%GID) THEN
          CALL ERROR$MSG('GRID ID MUST BE UNIQUE')
          CALL ERROR$I4VAL('GID1',GID)
          CALL ERROR$I4VAL('GID2',THISARR(IAT)%GID)
          CALL ERROR$STOP('LATTICE$ENERGY')
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  ENERGY OF COMPRESSED IONS                                           ==
!     ==========================================================================
      ETOT=0.D0
      DO IAT=1,NAT
        CALL ONEATOM(IAT,Q(IAT),RAD(IAT),DETOT,DEDQ(IAT),DEDRAD(IAT))
      IF(TPR)WRITE(NFILO,FMT='(50("."),T1,"INDIVIDUAL ION-ENERGY",T50,F10.5)')DETOT
        ETOT=ETOT+DETOT
      ENDDO
      IF(TPR)WRITE(NFILO,FMT='(50("."),T1,"TOTAL ION-ENERGY",T50,F10.5)')ETOT
!
!     ==========================================================================
!     == MADELUNG ENERGY                                                      ==
!     ==========================================================================
      CALL MADELUNG(NAT,RBAS,POS,Q,DETOT,DEDQ1,FORCE)
! POT(:)=MATMUL(MADMAT,Q)/SCALE
! ETOT(:)=0.5D0*DOT_PRODUCT(POT,Q)
      DEDQ=DEDQ+DEDQ1
!      PRESSURE=ETOT/(3.D0*VOL)
      ETOT=ETOT+DETOT
      IF(TPR)WRITE(NFILO,FMT='(50("."),T1,"MADELUNG ENERGY",T50,F10.5)')DETOT
!
!     ==========================================================================
!     == ELECTROSTATIC OVERLAP ENERGY AND PAULI REPULSION                     ==
!     ==========================================================================
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
      ALLOCATE(RHO(NR))
      ALLOCATE(POT(NR))
      DETOT=0.D0
      DO IAT1=1,NAT
        DO IAT2=IAT1,NAT
          DO IT1=-NDIV,NDIV
            DO IT2=-NDIV,NDIV
              DO IT3=-NDIV,NDIV
                IF(IAT1.EQ.IAT2) THEN
                  IF(IT1.LT.0) CYCLE
                  IF(IT1.EQ.0) THEN
                    IF(IT2.LT.0) CYCLE
                    IF(IT2.EQ.0) THEN
                      IF(IT3.LE.0) CYCLE
                    END IF
                  END IF
                END IF
                DR(:)=RBAS(:,1)*REAL(IT1)+RBAS(:,2)*REAL(IT2)+RBAS(:,3)*REAL(IT3)
                DR(:)=POS(:,IAT2)-POS(:,IAT1)+DR(:)
                DIS=SQRT(SUM(DR(:)**2))
                IF(DIS.GT.RAD(IAT1)+RAD(IAT2)) CYCLE
!               == TRANSFORM DENSITY OF ATOM 1 TO SITE OF ATOM 2 ===========================
                AUX(:)=THISARR(IAT1)%RHO(:)
                CALL SPHERICAL$SHIFTCENTER(GID,NR,DR,1,AUX,1,RHO)
                POT(:)=THISARR(IAT2)%ETAPAULI(:,1)+0.5D0*THISARR(IAT2)%POTXC(:)+0.5D0*(THISARR(IAT2)%POTH(:)-Q(IAT2)/Y0/R(:))
                AUX(:)=RHO(:)*POT(:)*R(:)**2
                CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
                CALL RADIAL$VALUE(GID,NR,AUX1,RAD(IAT2),VAL)
                DETOT=DETOT+VAL
!               == TAKE CARE OF NUCLEUS ======================================================
                POT(:)=THISARR(IAT2)%POTH(:)-Q(IAT2)/Y0/R(:)
                CALL RADIAL$VALUE(GID,NR,POT,DIS,VAL)
                DETOT=DETOT-0.5D0*VAL*Y0*THISARR(IAT1)%AEZ/DIS
!               == TAKE CARE OF OVERLAP WITH NEIGBORING NUCLEUS ===============================
                IF(RAD(IAT1).GT.DIS) THEN
                   AUX(:)=THISARR(IAT1)%RHO(:)*Q(IAT2)/Y0*(1.D0/R(:)-1.D0/DIS)*R**2
                   CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
                   CALL RADIAL$VALUE(GID,NR,AUX1,DIS,VAL1)
                   CALL RADIAL$VALUE(GID,NR,AUX1,RAD(IAT1),VAL2)
                   DETOT=DETOT+0.5D0*(VAL2-VAL1)
                END IF
! 
!               == TRANSFORM DENSITY OF ATOM 2 TO SITE OF ATOM 1 ===========================
                AUX(:)=THISARR(IAT2)%RHO(:)
                CALL SPHERICAL$SHIFTCENTER(GID,NR,-DR,1,AUX,1,RHO)
                POT(:)=THISARR(IAT1)%ETAPAULI(:,1)+0.5D0*THISARR(IAT1)%POTXC(:)+0.5D0*(THISARR(IAT1)%POTH(:)-Q(IAT1)/Y0/R(:))
                AUX(:)=RHO(:)*POT(:)*R(:)**2
                CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
                CALL RADIAL$VALUE(GID,NR,AUX1,RAD(IAT1),VAL)
                DETOT=DETOT+VAL
!               == TAKE CARE OF NUCLEUS ======================================================
                POT(:)=THISARR(IAT1)%POTH(:)-Q(IAT1)/Y0/R(:)
                CALL RADIAL$VALUE(GID,NR,POT,DIS,VAL)
                DETOT=DETOT-0.5D0*VAL*Y0*THISARR(IAT2)%AEZ/DIS
!               == TAKE CARE OF OVERLAP WITH NEIGBORING NUCLEUS ===============================
                IF(RAD(IAT2).GT.DIS) THEN
                   AUX(:)=THISARR(IAT2)%RHO(:)*Q(IAT1)/Y0*(1.D0/R(:)-1.D0/DIS)*R**2
                   CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
                   CALL RADIAL$VALUE(GID,NR,AUX1,DIS,VAL1)
                   CALL RADIAL$VALUE(GID,NR,AUX1,RAD(IAT2),VAL2)
                   DETOT=DETOT+0.5D0*(VAL2-VAL1)
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF(TPR)WRITE(NFILO,FMT='(50("."),T1,"ELECTROSTATIC OVERLAP",T50,F10.5)')DETOT
      etot=etot+detot
      DEALLOCATE(AUX)
      DEALLOCATE(AUX1)
      DEALLOCATE(RHO)
      DEALLOCATE(POT)
!
!     ==========================================================================
!     == TRANSFORM RADIUS DERIVATIVE IN VOLUME DERIVATIVE                     ==
!     ==========================================================================
      IF(TPR)WRITE(NFILO,FMT='(50("."),T1,"TOTAL ENERGY",T50,F10.5)')ETOT
      DEDRAD=DEDRAD/(4.D0*PI*RAD**2)
      RETURN
      END
