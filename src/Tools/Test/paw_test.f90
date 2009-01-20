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
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM MAIN
      CALL TRACE$SETL4('ON',.FALSE.)
CALL BUILD()
STOP
CALL TESTATOMLIB()
STOP
CALL TESTSTRUCTURECONSTANTS
! CALL TEST
STOP
      CALL DEBUG_PETERS_NONSPH()
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BUILD()
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NFIL=1001
      CHARACTER(2)     :: EL
      CHARACTER(12)    :: ID
      CHARACTER(14)    :: SY
      CHARACTER(5)    :: TYPE="'HBS'"
      INTEGER(4)       :: I
      REAL(8)          :: AEZ,ZCORE,ZV
!     **************************************************************************
      OPEN(NFIL,FILE='STPS')
      WRITE(NFIL,FMT='("!SCNTL")')
      DO I=1,105
        CALL PERIODICTABLE$GET(I,'SYMBOL',EL)
        CALL PERIODICTABLE$GET(I,'Z',AEZ)
        CALL PERIODICTABLE$GET(I,'ZCORE',ZCORE)
        ZV=AEZ-ZCORE
        ID="'"//TRIM(EL)//'_HBS'//"'"
        SY="'"//EL//"'"
        IF(SY(3:3).EQ.' ')SY(3:3)='_'
        WRITE(NFIL,FMT='(T2,"!SETUP  ID=",A," EL=",A," ZV=",F3.0)') &
     &     TRIM(ID),TRIM(SY),ZV
        WRITE(NFIL,FMT='(T10,"TYPE=",A," RBOX/RCOV=2.0 RCSM/RCOV=0.25")') &
           TYPE
        WRITE(NFIL,FMT='(T10,"RCL/RCOV= 0.75 0.75 0.75 0.75")') 
        WRITE(NFIL,FMT='(T10,"LAMBDA= 6. 6. 6. 6.")') 
        WRITE(NFIL,FMT='(T4,"!GRID DMIN=5.E-6 DMAX=0.1 RMAX=7.0 !END")') 
        WRITE(NFIL,FMT='(T4,"!POT  POW=3. VAL=-1.2 RC/RCOV=0.67 !END")') 
        WRITE(NFIL,FMT='(T4,"!CORE POW=3. VAL=0.1  RC/RCOV=0.67 !END")') 
        WRITE(NFIL,FMT='(T2,"!END")') 
        WRITE(NFIL,*) 
!       == NOW THE SEMICORE SETUP
        IF(    (AEZ.GE.11.AND.AEZ.LE.12) &
    &      .OR.(AEZ.GE.19.AND.AEZ.LE.28) &
    &      .OR.(AEZ.GE.37.AND.AEZ.LE.46) &
    &      .OR.(AEZ.GE.55.AND.AEZ.LE.78) &
    &      .OR.(AEZ.GE.87.AND.AEZ.LE.110)) THEN
          ID="'"//TRIM(EL)//'_HBS_SC'//"'"
          SY="'"//EL//"'"
          IF(SY(3:3).EQ.' ')SY(3:3)='_'
          WRITE(NFIL,FMT='(T2,"!SETUP  ID=",A," EL=",A," ZV=",F3.0)') &
       &     TRIM(ID),TRIM(SY),ZV+8.d0
          WRITE(NFIL,FMT='(T10,"TYPE=",A," RBOX/RCOV=1.2 RCSM/RCOV=0.25")') &
           TYPE
          WRITE(NFIL,FMT='(T10,"RCL/RCOV= 0.5 0.5 0.5 0.5")') 
          WRITE(NFIL,FMT='(T10,"LAMBDA= 6. 6. 6. 6.")') 
          WRITE(NFIL,FMT='(T4,"!GRID DMIN=5.E-6 DMAX=0.1 RMAX=7.0 !END")') 
          WRITE(NFIL,FMT='(T4,"!POT  POW=3. VAL=-2.2 RC/RCOV=0.5  !END")') 
          WRITE(NFIL,FMT='(T4,"!CORE POW=3. VAL= 0.1 RC/RCOV=0.5  !END")') 
          WRITE(NFIL,FMT='(T2,"!END")') 
          WRITE(NFIL,*) 
        END IF
      ENDDO
      WRITE(NFIL,FMT='("!END")') 
      WRITE(NFIL,FMT='("!EOB")') 
      WRITE(NFIL,*)
      CLOSE(NFIL)

      RETURN
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
      REAL(8)                   :: RAD
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
CALL DFT$SETI4('TYPE',2)
CALL DFT$SETL4('SPIN',.FALSE.)
RAD=2.D0
PRINT*,'=========================== FADER ==================================='
CALL X_ONEATOM$NEW(1,8.D0)
CALL X_ONEATOM(1,0.D0,RAD,DETOT,DEDQ,DEDRAD)
PRINT*,'RAD ',RAD
PRINT*,'DETOT ',DETOT
PRINT*,'=========================== PETER ==================================='
CALL ONEATOM$NEW(1,8.D0)
CALL ONEATOM(1,0.D0,RAD,DETOT,DEDQ,DEDRAD)
PRINT*,'RAD ',RAD
PRINT*,'DETOT ',DETOT
STOP

!
!     ==========================================================================        
!     ==  THIS IS FOR SESM                                                    ==        
!     ==========================================================================        
      CALL BIGONE()
!
!     ==========================================================================
!     ==  THIS IS FOR SIMPLE SESM (ROLF'S PROJECT)                            ==
!     ==========================================================================
!      CALL LATTICE()
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
PRINT*,'NC ',NC
PRINT*,'LOFI ',THISARR(IAT1)%LOFI
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
        WRITE(*,FMT='("ITER",I5," EBG: ",20F10.5)')ITER,EBG(NC+1:)
        WRITE(*,FMT='("ITER",I5," FBG: ",20F10.5)')ITER,FBG(NC+1:)
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
      INTEGER(4)           :: IAT,IP,IF
      INTEGER(4),PARAMETER :: NP=30
      INTEGER(4),PARAMETER :: NF=1
      REAL(8)              :: X(NP)
      REAL(8)              :: Y(NP,NF)
!     **************************************************************************
      AEZARR(:)=(/12.D0,8.D0/)
      ALAT=4.214D0/0.529177D0
      RBASFCC(:,1)=(/0.0D0,0.5D0,0.5D0/)
      RBASFCC(:,2)=(/0.5D0,0.0D0,0.5D0/)
      RBASFCC(:,3)=(/0.5D0,0.5D0,0.0D0/)
      POS(:,:)=0.D0
      POS(1,2)=0.5D0
      RBAS=RBASFCC*ALAT
      POS=POS*ALAT
      Q(:)=(/-2.D0,2.D0/)
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
      RAD0(:)=RCOV(:)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      OPEN(UNIT=8,FILE='LATTICE_GNU.DAT')
      DO IP=1,NP
        X(IP)=7.8D0+0.2D0*REAL(IP-1)
        DO IF=1,NF
          RAD(2)=X(IP)
          RAD(1)=2.0D0+0.2D0*REAL(IF-1)
          CALL LATTICE$ENERGY(NAT,RBAS,POS,Q,RAD,ETOT)
          Y(IP,IF)=ETOT          
          WRITE(8,FMT='(30F10.5)')RAD(1),RAD(2),ETOT
        ENDDO
      ENDDO
      CLOSE(8)
      OPEN(UNIT=8,FILE='LATTICE.DAT')
      DO IP=1,NP
        WRITE(8,FMT='(30F10.5)')X(IP),Y(IP,:)
      ENDDO
      CLOSE(8)
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
      INTEGER(4),PARAMETER  :: NDIV=2
      INTEGER(4)            :: IAT,IAT1,IAT2,IT1,IT2,IT3
      INTEGER(4)            :: NR
      INTEGER(4)            :: GID
      REAL(8)               :: DETOT
      REAL(8)               :: DIS,DR(3)
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:),AUX1(:),RHO(:),POT(:)
      REAL(8)               :: VAL,VAL1,VAL2
      REAL(8)               :: PI,Y0
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      INTEGER(4)            :: NFILO
!     **************************************************************************
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
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
      ETOT=ETOT+DETOT
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
!==================================================================================
!===================================================================================
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TEST
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: L1X=2
      INTEGER(4),PARAMETER :: L2X=3
      INTEGER(4),PARAMETER :: LM1X=(L1X+1)**2
      INTEGER(4),PARAMETER :: LM2X=(L2X+1)**2
      REAL(8)   ,PARAMETER :: R0(3)=(/-1.D0,-1.D0,0.D0/)
      REAL(8)   ,PARAMETER :: DR1(3)=(/1.D-1,0.D0,0.D0/)
      REAL(8)   ,PARAMETER :: DR2(3)=(/0.D0,1.D-1,0.D0/)
      INTEGER(4),PARAMETER :: N1=31
      INTEGER(4),PARAMETER :: N2=21
      REAL(8)   ,PARAMETER :: RCENTER(3)=(/0.D0,0.D0,0.D0/)
      REAL(8)   ,PARAMETER :: ROFFCENTER(3)=(/1.4D0,0.3D0,0.2D0/)
      REAL(8)   ,PARAMETER :: K2=0.D0
      REAL(8)   ,PARAMETER :: VALX=2.D0
      INTEGER(4),PARAMETER :: NFIL1=10
      INTEGER(4),PARAMETER :: NFIL2=11
      INTEGER(4),PARAMETER :: NFIL3=12
      INTEGER(4)           :: I1,I2
      REAL(8)              :: H(LM1X)
      REAL(8)              :: H2(LM1X)
      REAL(8)              :: S(LM2X,LM1X)
      REAL(8)              :: J(LM2X)
      REAL(8)              :: R(3)
      REAL(8)              :: VAL,SVAR
      INTEGER(4)           :: LMSELECT
      REAL(8)              :: PI
      INTEGER(4)           :: LM1,LM2,L1,L2,M1,M2
      COMPLEX(8)           :: KK
!     **************************************************************************      
      PI=4.D0*ATAN(1.D0)
      IF(K2.GT.0.D0) THEN
        KK=CMPLX(1.D0,0.D0)*SQRT(K2)
      ELSE IF(K2.LT.0.D0) THEN
        KK=CMPLX(0.D0,1.D0)*SQRT(-K2)
      ELSE
        KK=(0.D0,0.D0)
      END IF
      OPEN(UNIT=NFIL1,FILE='H.DAT')
      OPEN(UNIT=NFIL2,FILE='JS.DAT')
      OPEN(UNIT=NFIL3,FILE='H-JS.DAT')
      LMSELECT=1
      CALL LMTO$STRUCTURECONSTANTS(RCENTER,L1X,ROFFCENTER,L2X,K2,S)
!
!     == COMPARE STRUCTURE CONSTANTS WITH DANIEL
!!$      DO LM1=1,LM1X
!!$        L1=SQRT(REAL(LM1-1)+1.D-7)
!!$        M1=1+L1**2+L1-LM1
!!$        DO LM2=1,LM2X
!!$          L2=SQRT(REAL(LM2-1)+1.D-7)
!!$          M2=1+L2**2+L2-LM2
!!$          CALL STRUCTURECONSTANTS(ROFFCENTER,L2,M2,RCENTER,L1,M1,KK,VAL)
!!$          IF(K2.GT.0.D0) THEN
!!$            VAL=-SQRT(K2)**(L1+L2+1)*VAL
!!$          ELSE IF(K2.LT.0.D0) THEN
!!$            VAL=2.D0/PI*SQRT(-K2)**(L1+L2+1)*VAL
!!$          ELSE
!!$            VAL=-VAL
!!$          END IF
!!$          WRITE(*,FMT='(6I5,3F10.5)')LM1,L1,M1,LM2,L2,M2,VAL,S(LM2,LM1),S(LM2,LM1)-VAL
!!$        ENDDO
!!$      ENDDO
!
!     == PUT FUNCTION AND EXPANSION ONI GRID
      DO I1=1,N1
        DO I2=1,N2
          R(:)=R0(:)+DR1(:)*REAL(I1-1,KIND=8)+DR2(:)*REAL(I2-1,KIND=8)
          CALL LMTO$SOLIDHANKEL(R-RCENTER,1.D-3,K2,LM1X,H)
          VAL=H(LMSELECT)
          VAL=MAX(-VALX,MIN(VALX,VAL))
          WRITE(NFIL1,FMT='(3F20.5)')R(1),R(2),VAL
          SVAR=VAL
          CALL LMTO$SOLIDBESSEL(R-ROFFCENTER,K2,LM2X,J)
          H2(:)=MATMUL(J,S)
          VAL=H2(LMSELECT)
          VAL=MAX(-VALX,MIN(VALX,VAL))
          WRITE(NFIL2,FMT='(3F20.5)')R(1),R(2),VAL
          WRITE(NFIL3,FMT='(3F20.5)')R(1),R(2),SVAR-VAL
        ENDDO
      ENDDO
      CLOSE(NFIL1)
      CLOSE(NFIL2)
      CLOSE(NFIL3)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TESTSTRUCTURECONSTANTS
!     **                                                                      **
!     **************************************************************************      
      USE PERIODICTABLE_MODULE
      USE ONEATOM_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NFIL=6
      INTEGER(4),PARAMETER :: NAT=7
      INTEGER(4),PARAMETER :: LXX=2
      REAL(8)   ,PARAMETER :: K2=0.D0
      REAL(8)   ,PARAMETER :: LUNIT=5.D0
      CHARACTER(4),PARAMETER:: ATOMID='FE'
      REAL(8)              :: RPOS(3,NAT)
      INTEGER(4)           :: LX(NAT)
      INTEGER(4)           :: N
      INTEGER(4)           :: NORB
      REAL(8)   ,ALLOCATABLE :: SBAR(:,:)
      REAL(8)              :: QBAR(LXX+1,NAT)
      INTEGER(4)           :: II,I
      INTEGER(4)           :: L
      INTEGER(4)           :: IAT
      REAL(8)              :: AEZ
      REAL(8)              :: RAD(NAT)
      REAL(8)              :: DETOT,DEDQ,DEDRAD
      REAL(8)              :: VAL,DER
      CHARACTER(32)        :: ROOTNAME='XX'
      REAL(8)      :: S1(4,1)
      INTEGER(4),PARAMETER :: LMSELECT=2
      REAL(8)   ,PARAMETER :: R0(3)=(/-1.D0,-1.D0,0.D0/)*10.D0
      REAL(8)   ,PARAMETER :: DR1(3)=(/1.D-1,0.D0,0.D0/)*10.D0
      REAL(8)   ,PARAMETER :: DR2(3)=(/0.D0,1.D-1,0.D0/)*10.D0
      INTEGER(4),PARAMETER :: N1=31
      INTEGER(4),PARAMETER :: N2=21
      REAL(8)   ,PARAMETER :: VALX=2.D-2
      INTEGER(4),PARAMETER :: NFIL1=10
      INTEGER(4)           :: I1,I2
      REAL(8)              :: R(3)
      REAL(8)              :: H((LXX+1)**2)
      REAL(8)   ,ALLOCATABLE :: C(:,:)
      REAL(8)   ,ALLOCATABLE :: QBARVEC(:)
!     **************************************************************************      
      LX(:)=LXX
      RPOS(:,1)=(/0.D0,0.D0,0.D0/)
      RPOS(:,2)=(/1.D0,0.D0,0.D0/)
      RPOS(:,3)=(/0.D0,1.D0,0.D0/)
      RPOS(:,4)=(/0.D0,0.D0,1.D0/)
      RPOS(:,5)=(/-1.D0,0.D0,0.D0/)
      RPOS(:,6)=(/0.D0,-1.D0,0.D0/)
      RPOS(:,7)=(/0.D0,0.D0,-1.D0/)
      RPOS(:,:)=RPOS(:,:)*LUNIT
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
!     == DETERMINE SCREENED STRUCTURE CONSTANTS                               ==
!     ==========================================================================
      CALL DFT$SETI4('TYPE',2)
      CALL DFT$SETL4('SPIN',.FALSE.)
      CALL PERIODICTABLE$GET(ATOMID,'Z',AEZ)
      CALL ONEATOM$NEW(1,AEZ)
      CALL ONEATOM(1,0.D0,20.D0,DETOT,DEDQ,DEDRAD)
      CALL PERIODICTABLE$GET(ATOMID,'R(COV)',RAD(1))
      RAD(2:)=RAD(1)
      THIS=>THISARR(1)
      QBAR(:,:)=0.D0
      DO L=0,LX(1)
        CALL RADIAL$VALUE(THIS%GID,THIS%NR,THIS%UDOT(:,L+1),RAD,VAL)
        CALL RADIAL$VALUE(THIS%GID,THIS%NR,THIS%UDOT(:,L+1),RAD,DER)
        CALL LMTO$Q(L,RAD(1),VAL,DER,K2,QBAR(L+1,1))
PRINT*,'-- ',L,RAD(1),VAL,DER,QBAR(L+1,1)
        QBAR(L+1,2:)=QBAR(L+1,1)
      ENDDO
!     
!     ==========================================================================
!     == DETERMINE SCREENED STRUCTURE CONSTANTS                               ==
!     ==========================================================================
      N=SUM((LX(:)+1)**2)
      NORB=(LX(1)+1)**2
      ALLOCATE(SBAR(N,NORB))
      ALLOCATE(QBARVEC(N))
      ALLOCATE(C(N,NORB))
      CALL LMTO$EXPANDQBAR(NAT,LXX,LX,QBAR,N,QBARVEC)
      CALL LMTO$CLUSTERSTRUCTURECONSTANTS(K2,NAT,RPOS,LX,QBARVEC,N,NORB,SBAR)
      CALL LMTO$KBARMULTICENTER(N,NORB,QBARVEC,SBAR,C)
!
!     ==========================================================================
!     == PRINT STRUCTURE CONSTANTS                                            ==
!     ==========================================================================
      II=0
      DO IAT=1,NAT
        WRITE(NFIL,FMT='(82("="),T20," STRUCTURECONSTANTS FOR ATOM ",I2)')IAT
        DO I=1,(LX(IAT)+1)**2
          II=II+1
          WRITE(NFIL,FMT='(20F10.5)')SBAR(II,:)
        ENDDO
      ENDDO
 !
!     ==========================================================================
!     == PLOT ENVELOPE FUNCTION                                               ==
!     ==========================================================================
      OPEN(UNIT=NFIL1,FILE='KBAR.DAT')
      DO I1=1,N1
        DO I2=1,N2
          R(:)=R0(:)+DR1(:)*REAL(I1-1,KIND=8)+DR2(:)*REAL(I2-1,KIND=8)
          CALL LMTO$SCREENEDSOLIDHANKEL(K2,NAT,RPOS,RAD,LX,N,NORB,C,R,H)
          VAL=H(LMSELECT)
          VAL=MAX(-VALX,MIN(VALX,VAL))
          WRITE(NFIL1,FMT='(3F20.5)')R(1),R(2),VAL
        ENDDO
      ENDDO
      CLOSE(NFIL1)
      DEALLOCATE(SBAR)
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE STRUCTURECONSTANTS_MODULE
TYPE NMAT_TYPE
  INTEGER(4) :: IAT1=0
  INTEGER(4) :: IAT2=0
  INTEGER(4) :: IT(3)=0
  INTEGER(4) :: LMN1=0
  INTEGER(4) :: LMN2=0
  REAL(8)  ,POINTER :: S(:,:)
  TYPE(NMAT_TYPE),POINTER :: NEXT
END TYPE NMAT_TYPE
TYPE(NMAT_TYPE),TARGET,SAVE  :: FIRST
TYPE(NMAT_TYPE),POINTER,SAVE :: THIS
TYPE(NMAT_TYPE),POINTER,SAVE :: PREVIOUS
CONTAINS
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NMAT_SELECT(IAT1,IAT2,IT)
!     ** POINTS THIS ONTO THE DESIRED ELEMENT OR DISASSOCIATES IT          **
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IT(3)
!     **************************************************************************
      THIS=>FIRST
      NULLIFY(PREVIOUS)
      DO
        IF(THIS%IAT1.EQ.IAT1) THEN
          IF(THIS%IAT2.NE.IAT2) THEN
            IF(THIS%IT(1).EQ.IT(1))THEN
              IF(THIS%IT(2).EQ.IT(2)) THEN
                IF(THIS%IT(3).EQ.IT(3)) THEN
                  EXIT     ! SUCCESS: CORRECT MEMBER FOUND
                ELSE IF(THIS%IT(3).GT.IT(3)) THEN
                  NULLIFY(THIS)
                  EXIT
                END IF
              ELSE IF(THIS%IT(2).GT.IT(2)) THEN
                NULLIFY(THIS)
                EXIT
              END IF
            ELSE IF(THIS%IT(1).GT.IT(1)) THEN
              NULLIFY(THIS)
              EXIT
            END IF
          ELSE IF(THIS%IAT2.GT.IAT2) THEN
            NULLIFY(THIS)
            EXIT
          END IF
        ELSE IF(THIS%IAT1.GT.IAT1) THEN
          NULLIFY(THIS)
          EXIT
        END IF
!       == SEARCH FOR THE NEXT MEMBER IN THE LIST
        IF(ASSOCIATED(THIS%NEXT)) THEN
          PREVIOUS=>THIS
          THIS=>THIS%NEXT
        ELSE
          NULLIFY(THIS)
          EXIT
        END IF
      ENDDO
      RETURN
      END SUBROUTINE NMAT_SELECT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE NMAT_INSERT()
!     ** INSERT A NEW NEIGHBOR PAIR AFTER PREVIOUS                            **
      IMPLICIT NONE
!     **************************************************************************
      IF(ASSOCIATED(PREVIOUS)) THEN
        ALLOCATE(THIS)
        THIS%NEXT=>PREVIOUS%NEXT
        PREVIOUS%NEXT=>THIS
      ELSE
        IF(FIRST%IAT1.NE.0) THEN 
!         == FIRST ELEMENT HAS ALREADY BEEN DEFINED. CREATE NEW ONE ON THE    ==
!         == FIRST POSITION AND SHIFT THE FIRST TO THE SECOND PLACE           ==
          ALLOCATE(THIS)
          THIS%NEXT=>FIRST%NEXT
          FIRST%NEXT=>THIS
          THIS%IAT1=FIRST%IAT1
          THIS%IAT2=FIRST%IAT2
          THIS%IT=FIRST%IT
          THIS%LMN1=FIRST%LMN1
          THIS%LMN2=FIRST%LMN2
          THIS%S=>FIRST%S
          THIS=>FIRST
        END IF
      END IF
!     == BLANK OUT THE NEW ELEMENT
      THIS%IAT1=0
      THIS%IAT2=0
      THIS%IT=0
      THIS%LMN1=0
      THIS%LMN2=0
      NULLIFY(THIS%S)
      RETURN
      END SUBROUTINE NMAT_INSERT
END MODULE STRUCTURECONSTANTS_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SMAT$GET(IAT1,IAT2,IT,LMN1,LMN2,S,TCHK)
      USE STRUCTURECONSTANTS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IT(3)
      INTEGER(4),INTENT(IN) :: LMN1
      INTEGER(4),INTENT(IN) :: LMN2
      REAL(8)   ,INTENT(OUT):: S(LMN1,LMN2)
      LOGICAL(4),INTENT(OUT):: TCHK
!     **************************************************************************      
      CALL NMAT_SELECT(IAT1,IAT2,IT)
      TCHK=ASSOCIATED(THIS)
      IF(TCHK) THEN
        IF(LMN1.NE.THIS%LMN1.OR.LMN2.NE.THIS%LMN2) THEN
          CALL ERROR$STOP('SMAT$GET')
        END IF
        S=THIS%S
      ELSE
        S(:,:)=0.D0
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SMAT$SET(IAT1,IAT2,IT,LMN1,LMN2,S)
      USE STRUCTURECONSTANTS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT1
      INTEGER(4),INTENT(IN) :: IAT2
      INTEGER(4),INTENT(IN) :: IT(3)
      INTEGER(4),INTENT(IN) :: LMN1
      INTEGER(4),INTENT(IN) :: LMN2
      REAL(8)   ,INTENT(IN) :: S(LMN1,LMN2)
      LOGICAL(4)            :: TCHK
!     **************************************************************************      
      CALL NMAT_SELECT(IAT1,IAT2,IT)
      TCHK=ASSOCIATED(THIS)
      IF(TCHK) THEN
        IF(LMN1.NE.THIS%LMN1.OR.LMN2.NE.THIS%LMN2) THEN
          DEALLOCATE(THIS%S)
          THIS%LMN1=LMN1
          THIS%LMN2=LMN2
          ALLOCATE(THIS%S(LMN1,LMN2))
        END IF
        THIS%S=S
      ELSE
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SMAT$REPORT(NFIL,TITLE)
      USE STRUCTURECONSTANTS_MODULE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: NFIL
      CHARACTER(*)   ,INTENT(IN) :: TITLE
      TYPE(NMAT_TYPE),POINTER:: CURRENT
!     **************************************************************************      
      WRITE(NFIL,FMT='(82("="),T20,A)')TITLE
      IF(FIRST%IAT1.EQ.0) THEN
        WRITE(NFIL,FMT='("LIST IS EMPTY")')
      END IF
      CURRENT=>FIRST
      DO 
        WRITE(NFIL,FMT='("IAT1=",I5," IAT2=",I5," IT=",3I2)')CURRENT%IAT1,CURRENT%IAT2,CURRENT%IT(:)
        IF(.NOT.ASSOCIATED(CURRENT%NEXT)) EXIT
        CURRENT=>CURRENT%NEXT
      ENDDO
      WRITE(NFIL,FMT='("END OF LIST")')
      RETURN
      END
      

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TESTATOMLIB()
      REAL(8)   ,PARAMETER :: R1=1.056D-4
      REAL(8)   ,PARAMETER :: DEX=0.05
      INTEGER(4),PARAMETER :: NR=250
      INTEGER(4),PARAMETER :: NP=20
!      REAL(8)   ,PARAMETER :: XMIN=-7.D0,XMAX=+3.D0
      REAL(8)   ,PARAMETER :: XMIN=0.D0,XMAX=+20.D0
      REAL(8)              :: POTPAR(NP,108)
      INTEGER(4)           :: GID
      CHARACTER(128)       :: KEY='START,REL,NONSO,ROBUST'
      REAL(8)              :: RBOX
      INTEGER(4),PARAMETER :: NX=40
      INTEGER(4)           :: NB
      REAL(8)              :: AEZ
      INTEGER(4)           :: LOFI(NX)
      INTEGER(4)           :: SO(NX)
      INTEGER(4)           :: NN(NX)
      REAL(8)              :: F(NX)
      REAL(8)              :: EOFI(NX)
      REAL(8)              :: PHI(NR,NX)
      REAL(8)              :: ETOT
      REAL(8)              :: POT(NR),POT1(NR),POT2(NR)
      REAL(8)              :: R(NR)
      REAL(8)              :: Y,RI,Y1
      INTEGER(4)           :: NFILDAT=10
      INTEGER(4)           :: NFILPOT=11
      INTEGER(4)           :: I,IR,J
      CHARACTER(64)        :: STRING
      REAL(8)              :: PI,Y0
      REAL(8)              :: AUX(NR),NUCPOT(NR)
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$NEW('SHLOG',GID)
      CALL RADIAL$SETR8(GID,'R1',R1)
      CALL RADIAL$SETR8(GID,'DEX',DEX)
      CALL RADIAL$SETI4(GID,'NR',NR)
      CALL RADIAL$R(GID,NR,R)
      RBOX=R(NR-2)
      OPEN(UNIT=NFILDAT,FILE='ATOMLIB_E.DAT')
      DO I=1,108,1
        AEZ=REAL(I,KIND=8)
        CALL ATOMLIB$AESCF(GID,NR,KEY,RBOX,AEZ,NX,NB,LOFI,SO,F,NN,ETOT,POT,EOFI,PHI)
        PRINT*,'ETOT',AEZ,ETOT
        WRITE(NFILDAT,FMT='(41F12.2)')AEZ,EOFI

   STRING=TRIM('APOT'//TRIM(ADJUSTL(STRING)))//'.DAT'
   OPEN(UNIT=NFILPOT,FILE=STRING)
    DO IR=1,NR
      WRITE(NFILPOT,*)R(IR),-POT(IR)*R(IR)
    ENDDO      
   CLOSE(NFILPOT)
CYCLE

!
!       == DETERMINE DATA FOR INTERPOLATION ======================================
        CALL RADIAL$NUCPOT(GID,NR,AEZ,NUCPOT)
        AUX=POT*R
        DO J=1,NP
          RI=XMIN+REAL(J-1,KIND=8)*(XMAX-XMIN)/REAL(NP-1,KIND=8)
!PRINT*,'RI ',RI,R(NR)
          CALL RADIAL$VALUE(GID,NR,AUX,RI,POTPAR(J,I))
        ENDDO
        CALL ATOMLIB_UNWRAPPOT(GID,NR,AEZ,NP,XMIN,XMAX,POTPAR(:,I),POT1)
!PRINT*,'MARKE 1'       
POT2=POT1
!CALL ATOMLIB$AESCF(GID,NR,'REL,NONSO',RBOX,AEZ,NX,NB,LOFI,SO,F,NN,ETOT,POT2,EOFI,PHI)
PRINT*,'========================= MARKE 2=============================='       
CYCLE
!
!       ====
        WRITE(STRING,*)I
        STRING=TRIM('POT'//TRIM(ADJUSTL(STRING)))//'.DAT'
        OPEN(UNIT=NFILPOT,FILE=STRING)
        DO IR=1,NR
!          Y=LOG(POT(IR)-NUCPOT(IR))
!          Y1=LOG(POT1(IR)-NUCPOT(IR))
IF(IR.LT.EXP(XMIN)) CYCLE
Y=POT(IR)*R(IR)**3
Y1=POT1(IR)
          WRITE(NFILPOT,*)R(IR),Y,Y1
        ENDDO
        CLOSE(NFILPOT)
      ENDDO
      CLOSE(NFILDAT)
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB_UNWRAPPOT(GID,NR,AEZ,NX,XMIN,XMAX,POTPAR,POT)
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      REAL(8)   ,INTENT(IN) :: XMIN
      REAL(8)   ,INTENT(IN) :: XMAX
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: POTPAR(NX)
      REAL(8)   ,INTENT(OUT):: POT(NR)
      REAL(8)               :: XI(NX)
      REAL(8)               :: R(NR)
      REAL(8)               :: PI,Y0
      REAL(8)               :: A(NX-3),B(NX-3),C(NX-3),D(NX-3)
      INTEGER(4)            :: I,IR,J
      REAL(8)               :: X(4),F(4)
      REAL(8)               :: XC
      REAL(8)               :: XR
      REAL(8)               :: NUCPOT(NR)
!     *************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      DO I=1,NX
        XI(I)=XMIN+REAL(I-1,KIND=8)*(XMAX-XMIN)/REAL(NX-1,KIND=8)
      ENDDO
      DO I=1,NX-3
        F(:)=POTPAR(I:I+3)
        X(:)=XI(I:I+3)
        XC=SUM(X)/4.D0
        X(:)=X(:)-XC
        C(I)=0.5D0*(F(1)-F(2)-F(3)+F(4))/(X(4)**2-X(3)**2)
        F(:)=F(:)-C(I)*X(:)**2
        A(I)=SUM(F)/4.D0
        F(:)=F(:)-A(I)
        F(:)=F(:)/X(:)
        D(I)=0.5D0*(F(1)-F(2)-F(3)+F(4))/(X(4)**2-X(3)**2)
        F(:)=F(:)-D(I)*X(:)**2
        B(I)=SUM(F)/4.D0
!       == SHIFT ORIGIN
        A(I)=A(I)-B(I)*XC+C(I)*XC**2-D(I)*XC**3
        B(I)=B(I)-2.D0*C(I)*XC+3.D0*D(I)*XC**2
        C(I)=C(I)-3.D0*D(I)*XC
!       == TEST
!        X(:)=XI(I:I+3)
!        F(:)=POTPAR(I:I+3)-A(I)-B(I)*X(:)-C(I)*X(:)**2-D(I)*X(:)**3
!        PRINT*,'DEV',MAXVAL(ABS(F))
      ENDDO
      DO IR=1,NR
        XR=LOG(R(IR))
   XR=R(IR)
        IF(XR.LT.XI(1)) THEN
          POT(IR)=POTPAR(1)
        ELSE IF(XR.LT.XI(2)) THEN
          POT(IR)=A(1)+B(1)*XR+C(1)*XR**2+D(1)*XR**3
        ELSE IF(XR.LT.XI(NX-1)) THEN
          DO J=1,NX-3
            IF(XR.GT.XI(J+1).AND.XR.LE.XI(J+3)) THEN
              POT(IR)=A(J)+B(J)*XR+C(J)*XR**2+D(J)*XR**3
              EXIT
            END IF
          ENDDO
        ELSE 
           POT(IR)=A(NX-3)+B(NX-3)*XR+C(NX-3)*XR**2+D(NX-3)*XR**3
        END IF
      ENDDO
!      POT=EXP(POT)
!      CALL RADIAL$NUCPOT(GID,NR,AEZ,NUCPOT)
!      POT=POT+NUCPOT
      RETURN 
      END
