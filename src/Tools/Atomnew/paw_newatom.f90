! MIXPOT COULD USE SOME PRECONDITIONING OR BROYDEN
! IT IS NOT CLEAR WHICH REFERENCE ENERGY SHOULD BE TAKEN FOR DREL
! DERIVATIVE DOES NOT HANDLE THE FIRST POINT RIGHT?
! INCLUDE SMALL COMPONENT FOR NORMALIZATION AND DENSITY
      PROGRAM TEST
      USE PERIODICTABLE_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTERFACE 
        SUBROUTINE LINKEDLIST$WRITE(LL_,NFIL,CID_)
        USE LINKEDLIST_MODULE, ONLY: LL_TYPE 
        TYPE(LL_TYPE),INTENT(IN) :: LL_
        INTEGER(4)   ,INTENT(IN) :: NFIL
        CHARACTER(*) ,INTENT(IN),OPTIONAL :: CID_ ! RELEVANT PROCESSOR GROUP (SEE MPE OBECT)
        END SUBROUTINE LINKEDLIST$WRITE
      END INTERFACE
      INTERFACE 
        SUBROUTINE LINKEDLIST$read(LL_,NFIL,CID_)
        USE LINKEDLIST_MODULE, ONLY: LL_TYPE 
        TYPE(LL_TYPE),INTENT(IN) :: LL_
        INTEGER(4)   ,INTENT(IN) :: NFIL
        CHARACTER(*) ,INTENT(IN),OPTIONAL :: CID_ ! RELEVANT PROCESSOR GROUP (SEE MPE OBECT)
        END SUBROUTINE LINKEDLIST$read
      END INTERFACE
      TYPE(LL_TYPE) :: LL_STP
      TYPE(LL_TYPE) :: LL_CNTL
      INTEGER(4)             :: NR
      INTEGER(4)             :: GID
      INTEGER(4)             :: GID0
      REAL(8)                :: AEZ
      character(2)           :: elementsymbol
      REAL(8)                :: PSZ
      REAL(8)   ,ALLOCATABLE :: AEPOT(:)
      REAL(8)   ,ALLOCATABLE :: DREL(:)
      REAL(8)   ,ALLOCATABLE :: DRELZORA(:)
      INTEGER(4)             :: NB
      INTEGER(4),ALLOCATABLE :: LOFI(:)
      INTEGER(4),ALLOCATABLE :: SOFI(:)
      INTEGER(4),ALLOCATABLE :: NNOFI(:)
      REAL(8)   ,ALLOCATABLE :: FOFI(:)
      REAL(8)   ,ALLOCATABLE :: EOFI(:)
      LOGICAL   ,ALLOCATABLE :: TCORE(:)
      INTEGER(4)             :: IB,JB,IR,IB1,IB2,L1,L2,I1,I2
      REAL(8)   ,ALLOCATABLE :: RHOADD(:)
      REAL(8)   ,ALLOCATABLE :: R(:)    ! radial grod
      REAL(8)   ,ALLOCATABLE :: RCL(:)
      real(8)                :: rcut    ! partial waves truncated at this radius
      integer(4)             :: ircut
      REAL(8)                :: RCOV
      REAL(8)                :: Y0,PI,C0LL
      INTEGER(4)             :: NC      ! #(core states)
      INTEGER(4)             :: NN      ! #(nodes(
      REAL(8)   ,ALLOCATABLE :: AERHOC(:),AERHOV(:)
      REAL(8)   ,ALLOCATABLE :: PSRHOC(:),PSRHOV(:)
      REAL(8)   ,ALLOCATABLE :: UOFI(:,:),TUOFI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:,:),TAEPHI(:,:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:,:),TPSPHI(:,:,:)
      REAL(8)   ,ALLOCATABLE :: UPHI(:,:,:),TUPHI(:,:,:)
      REAL(8)   ,ALLOCATABLE :: PRO(:,:,:)
      REAL(8)   ,ALLOCATABLE :: VADD(:)
      REAL(8)                :: RCSM
      INTEGER(4)             :: POW_CORE,POW_POT
      INTEGER(4),ALLOCATABLE :: NPRO(:)
      REAL(8)                :: VAL0_CORE,VAL0_POT
      REAL(8)                :: RC_CORE,RC_POT
      REAL(8)                :: RC
      REAL(8)   ,ALLOCATABLE :: PSPOT(:)
      INTEGER(4)             :: IBC
      INTEGER(4)             :: NNSTART,NNEND
      INTEGER(4)             :: LMAX
      INTEGER(4)             :: NAUG
      INTEGER(4)             :: L,SO
      INTEGER(4)             :: IC ! STATE INDEX OF HIGHEST CORDE STATE
      REAL(8)   ,ALLOCATABLE :: UC(:)   ! HIGHEST CORE STATE FOR THIS L
      REAL(8)                :: E
      REAL(8)   ,ALLOCATABLE :: EGRID(:,:)
      REAL(8)   ,ALLOCATABLE :: fGRID(:,:)
      REAL(8)   ,ALLOCATABLE :: DT(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DO(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DH(:,:,:)
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)   ,ALLOCATABLE :: WORK1(:),WORK2(:),WORK3(:),work2d(:,:)
      REAL(8)                :: RMAX,R1,DEX
      REAL(8)                :: DMIN,DMAX,RX
      INTEGER(4)             :: NR1
      REAL(8)   ,ALLOCATABLE :: G(:)
      REAL(8)   ,ALLOCATABLE :: PHI2(:,:)
      REAL(8)   ,ALLOCATABLE :: PHI3(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHISAVE(:,:,:)
      REAL(8)   ,ALLOCATABLE :: TRANSPRO(:,:,:)
      REAL(8)   ,ALLOCATABLE :: TRANSPHI(:,:,:)
      REAL(8)   ,ALLOCATABLE :: TRANSPROINV(:,:,:)
      REAL(8)   ,ALLOCATABLE :: TRANSPHIINV(:,:,:)
      REAL(8)   ,ALLOCATABLE :: PROPSIBAR(:,:,:)
      REAL(8)   ,ALLOCATABLE :: SMAT(:,:)
      REAL(8)   ,ALLOCATABLE :: SVEC(:)
      REAL(8)                :: FAC,SVAR
      LOGICAL(4)             :: TCHK
      CHARACTER(32)          :: EXT
      INTEGER(4)             :: ISVAR
      INTEGER(4)             :: NFIL
      INTEGER(4)             :: IDFT
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      CHARACTER(128)         :: CNTLNAME
      CHARACTER(128)         :: ROOTNAME
      INTEGER(4)             :: I
      INTEGER(4)             :: NFILO
      LOGICAL(4)             :: TLOG=.TRUE.
      INTEGER(4)             :: GID2
      INTEGER(4)             :: NR2,ll,iaug
      real(8)                :: ev
      real(8)                :: svar1,svar2,svar3
      integer(4)             ::ISVAR1ARR(1)  ! USED TO RESHAPE AN ARRAY OF LENGTH 1
!     ****************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=Y0
      CALL CONSTANTS('EV',EV)
!     =================================================================
!     ==  DEFINE FILES                                               ==
!     =================================================================
      CALL GETARG(1,CNTLNAME)
      I=INDEX(CNTLNAME,'.',BACK=.TRUE.)
      ROOTNAME=CNTLNAME(1:I-1)
      CALL FILEHANDLER$SETROOT(ROOTNAME)
      CALL FILEHANDLER$SETFILE('PROT',.FALSE.,TRIM(ROOTNAME)//-'.PROT')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','FORM','FORMATTED')
!
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==================================================================
!     == DEFINE INPUT FILE                                            ==
!     ==================================================================
      CALL FILEHANDLER$SETFILE('INPUT',.FALSE.,TRIM(ROOTNAME)//-'.ACNTL')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','FORM','FORMATTED')
!
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL FILEHANDLER$UNIT('INPUT',NFIL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL)
!
!CALL LINKEDLIST$SELECT(LL_CNTL,'~')
!CALL LINKEDLIST$REPORT(LL_CNTL,NFILO)
!
!     ==================================================================
!     == DEFINE SETUP FILE                                            ==
!     ==================================================================
      CALL FILEHANDLER$SETFILE('SETUP',.FALSE.,TRIM(ROOTNAME)//-'.STP')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','FORM','FORMATTED')
!
!     ================================================================
!     == IDENTIFY ATOM                                              ==
!     ================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')

      CALL LINKEDLIST$EXISTD(LL_CNTL,'AEZ',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'AEZ',1,AEZ)
        CALL PERIODICTABLE$GET(NINT(AEZ),'SYMBOL',ELEMENTSYMBOL)
      ELSE
        CALL LINKEDLIST$EXISTD(LL_CNTL,'EL',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('EITHER AEZ OR EL NEEDS TO BE SPECIFIED')
          CALL ERROR$STOP('MAIN')
        END IF
        CALL LINKEDLIST$GET(LL_CNTL,'EL',1,ELEMENTSYMBOL)
        CALL PERIODICTABLE$GET(ELEMENTSYMBOL,'Z',AEZ)
      END IF
!
      CALL PERIODICTABLE$GET(NINT(AEZ),'R(COV)',RCOV)
!
      CALL REPORT$TITLE(NFILO,'ELEMENT')
      CALL REPORT$chVAL(NFILO,'ELEMENT',ELEMENTSYMBOL)
      CALL REPORT$R8VAL(NFILO,'Z',AEZ,'')
!
!     ================================================================
!     == DEFINE RADIAL GRID
!     ================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'GRID')
      DMIN=5.D-6
      DMAX=5.D-1
      RX=1.056D-4*EXP(5.D-2*249.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'DMIN',1,DMIN)
      CALL LINKEDLIST$GET(LL_CNTL,'DMAX',1,DMAX)
      CALL LINKEDLIST$GET(LL_CNTL,'RMAX',1,RX)
      CALL LINKEDLIST$GET(LL_CNTL,'RCUT/RCOV',1,RCUT)
      RCUT=RCUT*RCOV   ! projectors and partial waves truncated beyond this point
!
!     == DEFINE GRID ==================================================
      CALL GRIDPARAMETERS(DMIN,DMAX,RX,R1,DEX,NR)
      CALL RADIAL$NEW('SHLOG',GID)
      CALL RADIAL$SETR8(GID,'R1',R1)
      CALL RADIAL$SETR8(GID,'DEX',DEX)
      CALL RADIAL$SETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      do ir=1,nr
        ircut=ir
        if(r(ir).gt.rcut) exit
      enddo
!
!     == REPORT GRID ==================================================
      CALL REPORT$TITLE(NFILO,'RADIAL GRID')
      CALL REPORT$R8VAL(NFILO,'R1',R1,'A0')
      CALL REPORT$R8VAL(NFILO,'DEX',DEX,'')
      CALL REPORT$I4VAL(NFILO,'#(RADIAL GRID POINTS)',NR,'')
      CALL REPORT$R8VAL(NFILO,'INNERMOST SPACING',DMIN,'A0')
      CALL REPORT$R8VAL(NFILO,'OUTERMOST SPACING',DMAX,'A0')
      CALL REPORT$R8VAL(NFILO,'OUTERMOST GRID POINT',RX,'A0')
      CALL REPORT$R8VAL(NFILO,'TRANCATION RADIUS',RCUT,'A0')
      CALL REPORT$R8VAL(NFILO,'COVALENT RADIUS',RCOV,'A0')
!
!     ================================================================
!     == DEFINE DFT FUNCTIONAL
!     ================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
      CALL LINKEDLIST$GET(LL_CNTL,'DFT',1,IDFT)
!
      CALL DFT$SETI4('TYPE',IDFT)
      CALL DFT$SETL4('SPIN',.FALSE.)
!
!     == REPORT =====================================================
      CALL DFT$REPORT(NFILO)
!
!     ================================================================
!     == DEFINE STATES AND THEIR OCCUPATIONS
!     ================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$NLISTS(LL_CNTL,'STATE',NB)
      ALLOCATE(LOFI(NB))
      ALLOCATE(SOFI(NB))
      ALLOCATE(NNOFI(NB))
      ALLOCATE(FOFI(NB))
      ALLOCATE(TCORE(NB))
      NC=0
      DO IB=1,NB
        CALL LINKEDLIST$SELECT(LL_CNTL,'STATE',IB)
!       == MAIN ANGULAR MOMENTUM QUANTUM NUMBER
        CALL LINKEDLIST$GET(LL_CNTL,'L',1,LOFI(IB))
!       == #(NODES) PRINCIPAL QUANTUM NUMBER IS READ AND CONVERTED
        CALL LINKEDLIST$GET(LL_CNTL,'N',1,NNOFI(IB))
        NNOFI(IB)=NNOFI(IB)-LOFI(IB)-1
!       == SOFI DESCRIBES SPIN-ORBIT COUPLING:
!       ==  0: NO SPIN ORBIT COUPLING
!       ==  1: PARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM
!       == -1: ANTIPARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM
        SOFI(IB)=0
!       == OCCUPATION
        FOFI(IB)=2.D0*REAL(2*LOFI(IB)+1)
        CALL LINKEDLIST$EXISTD(LL_CNTL,'F',1,TCHK)
!       == SPECIFY IF CORE OR VALENCE
        IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'F',1,FOFI(IB))
        CALL LINKEDLIST$GET(LL_CNTL,'CORE',1,TCORE(IB))
        IF(TCORE(IB)) NC=NC+1
        CALL LINKEDLIST$SELECT(LL_CNTL,'..')
      ENDDO
      DO IB=1,NC
        IF(.NOT.TCORE(IB)) THEN
          CALL ERROR$MSG('CORE STATES MUST PRECEED VALENCE STATES')
          CALL ERROR$I4VAL('NC',NC)
          CALL ERROR$I4VAL('IB',IB)
          CALL ERROR$STOP('MAIN')
        END IF
      ENDDO
      CALL REPORT$TITLE(NFILO,'INPUT ON STATES')
      DO IB=1,NB
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB)
      ENDDO
!
!     ================================================================
!     == PERFORM SELF CONSISTENT ALL-ELECTRON CALCULATION           ==
!     ================================================================
      ALLOCATE(EOFI(NB)) 
      ALLOCATE(DREL(NR))
      ALLOCATE(DRELZORA(NR))
      ALLOCATE(AEPOT(NR))      
      ALLOCATE(RHOADD(NR))
!
      CALL NUCPOT(AEZ,NR,R,AEPOT)
      AEPOT(:)=AEPOT(:)/Y0
      RHOADD=0.D0
      CALL AESCF(GID,NR,NB,LOFI,SOFI,FOFI,NNOFI,AEZ,RHOADD,DREL,AEPOT,EOFI)
      DRELZORA(:)=DREL(:)
!
!     == REPORT ======================================================== 
      CALL REPORT$TITLE(NFILO,'EIGENSTATES AFTER SCF CALCULATION')
      DO IB=1,NB
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
     &                ," E[H]=",F20.9," E[ev]=",F20.9)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/ev
      ENDDO
      DEALLOCATE(RHOADD)
!
!     ================================================================
!     == CALCULATE CORE AND VALENCE DENSITY           
!     ================================================================
      ALLOCATE(AERHOC(NR))
      ALLOCATE(AERHOV(NR))
!
      CALL AERHO('FULL',GID,NR,NC,LOFI(1:NC),SOFI(1:NC),FOFI(1:NC),NNOFI(1:NC) &
     &          ,DREL,AEPOT,AERHOC,EOFI(1:NC))
!
!     == RELATIVISTIC CORRECTION WILL BE USED FROM THE VALENCE STATES
      CALL AERHO('FULL',GID,NR,NB-NC,LOFI(NC+1:NB),SOFI(NC+1:NB),FOFI(NC+1:NB) &
     &          ,NNOFI(NC+1:NB),DREL,AEPOT,AERHOV,EOFI(NC+1:NB))
!
!     == REPORT ======================================================
      WRITE(NFILO,FMT='("EIGENSTATES AFTER NON-SCF CALCULATION IN SCF POT")')
      DO IB=1,NB
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
     &                ," E[H]=",F20.9," E[ev]=",F20.9)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/ev
      ENDDO
!
!     ================================================================
!     == CALCULATE NODELESS WAVE FUNCTIONS
!     ================================================================
      ALLOCATE(UOFI(NR,NB))  ; uofi=0.d0
      ALLOCATE(TUOFI(NR,NB)) ; tuofi=0.d0
      CALL NODELESS('EFFZORA',GID,NR,NB,LOFI,SOFI,NNOFI,DREL,AEPOT &
     &             ,UOFI,TUOFI,EOFI)
!
!     == REPORT =======================================================
      WRITE(NFILO,FMT='("EIGENSTATES FROM NODELESS CALCULATION")')
      DO IB=1,NB
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
     &                ," E[H]=",F20.9," E[ev]=",F20.9)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/ev
      ENDDO
!     CALL HAMTEST(GID,NR,AEPOT,NB,LOFI,SOFI,UOFI,TUOFI)
!
!     ================================================================
!     == SELF-CONSISTENT ATOM FINISHED; WRITE RESULT                ==
!     ================================================================
      CALL WRITEPHI('ALLNODELESS.DAT',GID,NR,NB,UOFI)
      CALL WRITEPHI('ALLNODELESS_T.DAT',GID,NR,NB,TUOFI)
!
!     ================================================================
!     == TOTAL ENERGY                                               ==
!     ================================================================
      CALL AEETOT(GID,NR,NB,AEZ,NC,EOFI,FOFI,AERHOC,AERHOV,AEPOT)
!
!     ================================================================
!     ================================================================
!     == DEFINE AUGMENTATION                     
!     ================================================================
!     ================================================================
!
!     ================================================================
!     ==  PSEUDO CORE                                               ==
!     ================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'PSCORE')
      CALL LINKEDLIST$GET(LL_CNTL,'RC/RCOV',1,RC_CORE)!CUTOFF RADIUS FOR CORE PSEUDIZATION
      RC_CORE=RC_CORE*RCOV
      CALL LINKEDLIST$GET(LL_CNTL,'POWER',1,POW_CORE)  !POWER FOR CORE PSEUDIZATION
      CALL LINKEDLIST$GET(LL_CNTL,'VAL0',1,VAL0_CORE)  !VALUE AT THE ORIGIN FOR CORE PSEUDIZATION
      val0_core=val0_core/y0

      ALLOCATE(PSRHOC(NR))
      CALL PSEUDIZE(GID,NR,POW_CORE,VAL0_CORE,RC_CORE,AERHOC,PSRHOC)
!
!     == REPORT ======================================================
      CALL REPORT$TITLE(NFILO,'PSEUDO CORE DENSITY')
      CALL REPORT$R8VAL(NFILO,'RC',RC_CORE,'')
      CALL REPORT$I4VAL(NFILO,'POWER',POW_CORE,'')
      CALL REPORT$R8VAL(NFILO,'VAL(0)',VAL0_CORE*y0,'')
!
!     ================================================================
!     ==  PSEUDO POTENTIAL                                          ==
!     ================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'VTILDE')
      CALL LINKEDLIST$GET(LL_CNTL,'RC/RCOV',1,RC_POT)!CUTOFF RADIUS FOR POTENTIAL PSEUDIZATION
      RC_POT=RC_POT*RCOV
      CALL LINKEDLIST$GET(LL_CNTL,'POWER',1,POW_POT)  !POWER FOR POTENTIAL PSEUDIZATION 
      CALL LINKEDLIST$GET(LL_CNTL,'VAL0',1,VAL0_POT)  !VALUE AT THE ORIGIN FOR POTENTIAL PSEUDIZATION
      val0_pot=val0_pot/y0
!
      ALLOCATE(PSPOT(NR))
      CALL PSEUDIZE(GID,NR,POW_POT,VAL0_POT,RC_POT,AEPOT,PSPOT)
!
!     == REPORT ======================================================
      CALL REPORT$TITLE(NFILO,'PSEUDO POTENTIAL')
      CALL REPORT$R8VAL(NFILO,'RC',RC_POT,'')
      CALL REPORT$I4VAL(NFILO,'POWER',POW_POT,'')
      CALL REPORT$R8VAL(NFILO,'VAL(0)',VAL0_POT*y0,'')
!
!     ================================================================
!     ==  COMPENSATION GAUSSIAN                                     ==
!     ================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'COMPENSATE')
      CALL LINKEDLIST$GET(LL_CNTL,'RC/RCOV',1,RCSM)  
      RCSM=RCSM*RCOV
!
      CALL REPORT$TITLE(NFILO,'COMPENSATION DENSITY')
      CALL REPORT$R8VAL(NFILO,'RC',RCSM,'')
      SVAR=AEZ-SUM(FOFI(1:NC))       
      SVAR=SVAR*EXP(-(RCOV/RCSM)**2)/(SQRT(PI)*RCSM)**3
      CALL REPORT$R8VAL(NFILO,'VALUE AT THE COVALENT RADIUS',SVAR,'')
!
!     ================================================================
!     ==  AUGMENTATION                                              ==
!     ================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'AUGMENTATION')
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NPRO',1,TCHK)
      LMAX=2
      IF(TCHK) THEN
        CALL LINKEDLIST$SIZE(LL_CNTL,'NPRO',1,LMAX)
        LMAX=LMAX-1
      END IF
!
      ALLOCATE(NPRO(LMAX+1))
      NPRO=2
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'NPRO',1,NPRO)
!
      ALLOCATE(RCL(LMAX+1))
      RCL(:)=RCOV
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RC/RCOV',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'RC/RCOV',1,RCL)
      RCL(:)=RCL(:)*RCOV
!
      CALL REPORT$TITLE(NFILO,'AUGMENTATION')
      CALL REPORT$I4VAL(NFILO,'LMAX',LMAX,'')
      DO L=0,LMAX
        IF(L.EQ.0) THEN
          CALL REPORT$I4VAL(NFILO,'#(PROJECTORS)FOR L=0',NPRO(L+1),'')
        ELSE IF(L.EQ.1) THEN
          CALL REPORT$I4VAL(NFILO,'#(PROJECTORS)FOR L=1',NPRO(L+1),'')
        ELSE IF(L.EQ.2) THEN
          CALL REPORT$I4VAL(NFILO,'#(PROJECTORS)FOR L=2',NPRO(L+1),'')
        ELSE IF(L.EQ.3) THEN
          CALL REPORT$I4VAL(NFILO,'#(PROJECTORS)FOR L=3',NPRO(L+1),'')
        END IF
      ENDDO
      DO L=0,LMAX
        IF(L.EQ.0) THEN
          CALL REPORT$R8VAL(NFILO,'CUTOFF FOR L=0',RCL(L+1),'')
        ELSE IF(L.EQ.1) THEN
          CALL REPORT$R8VAL(NFILO,'CUTOFF FOR L=1',RCL(L+1),'')
        ELSE IF(L.EQ.2) THEN
          CALL REPORT$R8VAL(NFILO,'CUTOFF FOR L=2',RCL(L+1),'')
        ELSE IF(L.EQ.3) THEN
          CALL REPORT$R8VAL(NFILO,'CUTOFF FOR L=3',RCL(L+1),'')
        END IF
      ENDDO
!
!     ================================================================
!     == ALLOCATE ARRAYS
!     ================================================================
      ALLOCATE(PSRHOV(NR))
      ALLOCATE(VADD(NR))
      ALLOCATE(UC(NR))
      ALLOCATE(AUX(NR))
!
!
      NAUG=MAXVAL(NPRO)
      ALLOCATE(AEPHI(NR,NAUG,LMAX+1))  ; aephi=0.d0
      ALLOCATE(TAEPHI(NR,NAUG,LMAX+1)) ; taephi=0.d0
      ALLOCATE(PSPHI(NR,NAUG,LMAX+1))  ; psphi=0.d0
      ALLOCATE(TPSPHI(NR,NAUG,LMAX+1)) ; tpsphi=0.d0
      ALLOCATE(PRO(NR,NAUG,LMAX+1))    ; pro=0.d0
      ALLOCATE(UPHI(NR,NAUG,LMAX+1))   ; uphi=0.d0
      ALLOCATE(TUPHI(NR,NAUG,LMAX+1))  ; tuphi=0.d0
      ALLOCATE(TRANSPRO(NAUG,NAUG,LMAX+1))    ; transpro=0.d0
      ALLOCATE(TRANSPHI(NAUG,NAUG,LMAX+1))    ; transphi=0.d0
      ALLOCATE(TRANSPROINV(NAUG,NAUG,LMAX+1)) ; transproinv=0.d0
      ALLOCATE(TRANSPHIINV(NAUG,NAUG,LMAX+1)) ;TRANSPHIINV=0.d0
      ALLOCATE(DT(NAUG,NAUG,LMAX+1))
      ALLOCATE(DO(NAUG,NAUG,LMAX+1))
      ALLOCATE(DH(NAUG,NAUG,LMAX+1))
      ALLOCATE(PROPSIBAR(NAUG,NAUG,LMAX+1))
      ALLOCATE(EGRID(NAUG,LMAX+1))
      ALLOCATE(fGRID(NAUG,LMAX+1))
      ALLOCATE(PSPHISAVE(NR,NAUG,LMAX+1))
!
!     ================================================================
!     == CONSTRUCT PARTIAL WAVES AND PROJECTORS FOR EACH L          ==
!     ================================================================
      PSRHOV(:)=0.D0
      DO L=0,LMAX
        WRITE(NFILO,FMT='(30("="),2X,"L=",I1,2X,30("="))')L
        NAUG=NPRO(L+1)
!       == ext is used as ending for printouts of l-dependent FUNCTIONS
        IF(L.EQ.0) THEN
          EXT='_S.DAT'
        ELSE IF (L.EQ.1) THEN
          EXT='_P.DAT'
        ELSE IF (L.EQ.2) THEN
          EXT='_D.DAT'
        ELSE IF (L.EQ.3) THEN
          EXT='_F.DAT'
        END IF
!
!ATTENTION CURRENTLY NPRO MUST BE THE SAME FOR ALL L!!!!
!
!       =================================================================
!       ==  DETERMINE HIGHEST NODELESS CORE WAVE FUNCTION UC            ==
!       =================================================================
        UC(:)=0.D0
!
!       == determine index ic of highest core state. (-1 if there is none)
!       == determine highest core state uc
!       == determine NNstart, the number of nodes for the lowest valence state
        NNSTART=0   ! #(NODES OF FIRST  VALENCE WAVE FUNCTION)
        IC=-1
        DO IB=1,NC
          IF(LOFI(IB).NE.L) CYCLE
          IF(NNOFI(IB)+1.GT.NNSTART) THEN
            UC(:)=UOFI(:,IB)
            NNSTART=NNOFI(IB)+1  ! lowest valence state has one node more 
                                 ! than highest core state
            IC=IB
          END IF
        ENDDO        
!
!       ==  determine energies for partial wave construction ==================
        egrid(:,l+1)=0.d0
        fgrid(:,l+1)=0.d0
        TCHK=.FALSE.
        NNEND=0
        DO IB=NC+1,NB
          IF(LOFI(IB).NE.L) CYCLE
          IF(NNOFI(IB).GE.NNSTART) THEN
            iaug=NNOFI(IB)-NNSTART+1
            fGRID(iaug,L+1)=fOFI(IB)
            EGRID(iaug,L+1)=EOFI(IB)
            NNEND=MAX(iaug,NNEND)
            TCHK=.TRUE.
          END IF
        ENDDO
! the lowest state of each angular momentum shall be a bound state.
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('LOWEST VALENCE STATE NOT FOUND')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('MAIN')
        END IF
!
        CALL REPORT$I4VAL(NFILO,'NN OF HIGHEST CORE STATES ',NNSTART-1,'')
!
!       =================================================================
!       ==  DETERMINE NODELESS ENERGY EXPANSION COEFFICIENTS           ==
!       =================================================================
!       use relativistic corrections from above
        SO=0
        DO IB=1,NAUG
          IF(IB.LE.NNEND) THEN
            NN=IB-1
            CALL ONENODELESS('BOUND',GID,NR,L,NN,SO,DREL,AEPOT,UC &
       &                    ,EGRID(IB,L+1),UPHI(:,IB,L+1),TUPHI(:,IB,L+1))
            CALL REPORT$R8VAL(NFILO,'ENERGY SUPPORT GRID POINT (BOUND)',EGRID(IB,L+1),'H')
          ELSE
            IF(IB.EQ.1.AND.NNEND.GT.0) THEN   
              EGRID(IB,L+1)=0.D0
            ELSE
              EGRID(IB,L+1)=EGRID(NNEND,L+1)+0.5D0*REAL(IB-NNEND,KIND=8)**2
            END IF
            CALL ONENODELESS('SCATT',GID,NR,L,NN,SO,DREL,AEPOT,UC &
       &                    ,EGRID(IB,L+1),UPHI(:,IB,L+1),TUPHI(:,IB,L+1))
            CALL REPORT$R8VAL(NFILO,'ENERGY SUPPORT GRID POINT (UNBOUND)',EGRID(IB,L+1),'H')
!
          END IF
        ENDDO
!
!       =================================================================
!       ==  ALL ELECTRON FUNCTIONS BY CORE ORTHOGONALIZATION           ==
!       =================================================================
        CALL COREORTHOGONALIZE(GID,NR,NC,LOFI(1:NC),SOFI(1:NC) &
      &                      ,UOFI(:,1:NC),TUOFI(:,1:NC) &
      &                      ,L,NAUG,UPHI(:,:,L+1),TUPHI(:,:,L+1) &
      &                      ,AEPHI(:,:,L+1),TAEPHI(:,:,L+1))
!
!       =================================================================
!       ==  CHOSE NORM OF PARTIAL WAVES                                ==
!       =================================================================
        DO IB=1,NAUG
!         __ only bound states can be normalized. _______________________
!         __ for unbound states the boundary conditions at the nucleus __
!         __ are kept identical to that of the highest bound state ______
          IF(IB.LE.NNEND) THEN
            CALL RADIAL$INTEGRAL(GID,NR,(AEPHI(:,IB,L+1)*R(:))**2,SVAR)
            SVAR=1.D0/SQRT(SVAR)
          else
!           == choose 5th point to avoid divideing by small numbers for L.neq.0
            if(ib.eq.1) then
              svar=r(5)**l/aephi(5,ib,l+1)
            else
              svar=aephi(5,ib-1,l+1)/aephi(5,ib,l+1)
            end if
          end if
          UPHI(:,IB,L+1)  =  UPHI(:,IB,L+1)*SVAR
          TUPHI(:,IB,L+1) = TUPHI(:,IB,L+1)*SVAR
          AEPHI(:,IB,L+1) = AEPHI(:,IB,L+1)*SVAR
          TAEPHI(:,IB,L+1)=TAEPHI(:,IB,L+1)*SVAR
        ENDDO
!
!       =================================================================
!       ==  CONSTRUCT PS PARTIAL WAVES                                 ==
!       =================================================================
        PSPHI(:,:,L+1)=UPHI(:,:,L+1)
        TPSPHI(:,:,L+1)=TUPHI(:,:,L+1)
        RC=RCL(L+1)
!       USE METHOD ADJUSTING PSEUDOPOTENTIAL
        CALL CONSTRUCTPSPHI1(GID,NR,RC,POW_POT,PSPOT,L,NAUG &
     &                      ,EGRID(:,L+1),PSPHI(:,:,L+1),TPSPHI(:,:,L+1))
!PRINT*,'V0',-(TPSPHI(1,:,L+1)/PSPHI(1,:,L+1)-EGRID(:,L+1))/Y0
!       == USE METHOD MATCHING BESSEL FUNCTIONS
!       CALL CONSTRUCTPSPHI(GID,NR,RC,L,NAUG,PSPHI(:,:,L+1),TPSPHI(:,:,L+1))
!
!       =================================================================
!       ==  CONSTRUCT BARE PROJECTOR FUNCTIONS                         ==
!       =================================================================
        CALL BAREPROJECTORS(GID,NR,L,EGRID(:,L+1),PSPOT,10.D0,NAUG &
       &     ,PSPHI(:,:,L+1),TPSPHI(:,:,L+1),PRO(:,:,L+1))
!
!       =================================================================
!       ==  CONSTRUCT 1C KINETIC ENERGY AND OVERLAP                    ==
!       =================================================================
        CALL DTDO(GID,NR,NAUG,AEPHI(:,:,L+1),TAEPHI(:,:,L+1) &
       &         ,PSPHI(:,:,L+1),TPSPHI(:,:,L+1),DT(:,:,L+1),DO(:,:,L+1))
        CALL MAKEDATH(GID,NR,AEPOT,PSPOT,L,NAUG,AEPHI(:,:,L+1),PSPHI(:,:,L+1) &
       &             ,DT(:,:,L+1),DH(:,:,L+1))
!
!       ======================================================================
!       == ADD TO PSEUDO DENSITY (MUST BE DONE BEFORE BIORTHOGONALIZATION)  ==
!       == ASSUMES THAT THERE IS ONE PARTIAL WAVE FOR EACH OCCUPIED STATE   ==
!       ======================================================================
        DO IB=NC+1,NB
          IF(LOFI(IB).NE.L) CYCLE
          IAUG=NNOFI(IB)-NNSTART+1
          AUX(:)=(AEPHI(:,IAUG,L+1)*R(:))**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          SVAR=FOFI(IB)/SVAR*Y0
          PSRHOV(:)=PSRHOV(:)+SVAR*PSPHI(:,IAUG,L+1)**2
        ENDDO
!
!       =================================================================
!       ==  TEST RESULT BY RECALCULATING PSPHI FROM BARE PROJECTORS
!       =================================================================
        WRITE(NFILO,FMT='(72("=")/72("="),T20,"  BARE QUANTITIES DONE: TESTS...  "/72("="))')
        ALLOCATE(PHI2(NR,NAUG))
        ALLOCATE(PHI3(NR,NAUG))
        ALLOCATE(G(NR))
        WRITE(NFILO,FMT='("TEST: SOLVE [T+PSV-E|F>=|PRO>")')
        WRITE(NFILO,FMT='("     AND COMPARE |F> WITH |PSPHI> ... (MUST VANISH)")')
        DO IB=1,NAUG
          G(:)=PRO(:,IB,L+1)
          CALL RADIAL$SCHRODINGER(GID,NR,PSPOT,DREL,0,G,L,EGRID(IB,L+1) &
       &                        ,1,PHI2(:,IB))
!         == ENSURE SAME NORM AS PSPHI BY ADDING THE HOMOGENEOUS SOLUTION
          G=0.D0
          CALL RADIAL$SCHRODINGER(GID,NR,PSPOT,DREL,0,G,L,EGRID(IB,L+1) &
       &                        ,1,PHI3(:,IB))
          isvar1arr=maxloc(abs(psphi(:,ib,l+1)))
          ir=min(isvar1arr(1),ircut-1)
          SVAR=(PSPHI(IR,IB,L+1)-PHI2(IR,IB))/PHI3(IR,IB)
          PHI2(:,IB)=PHI2(:,IB)+PHI3(:,IB)*SVAR
          write(nfilo,fmt='(i5,e20.10)')ib,maxval(abs(phi2(:,ib)-PSPHI(:,ib,L+1)))
        ENDDO
        DEALLOCATE(G)
!
        ALLOCATE(WORK2D(NR,3*NAUG))
        WORK2D(:,1:NAUG)=PHI2(:,:)
        WORK2D(:,NAUG+1:2*NAUG)=PSPHI(:,:,L+1)
        WORK2D(:,2*NAUG+1:3*NAUG)=PSPHI(:,:,L+1)-PHI2(:,:)
        CALL WRITEPHI('TESTBAREPSPHI'//TRIM(EXT),GID,NR,3*NAUG,WORK2D)
        DEALLOCATE(WORK2D)
!
        DEALLOCATE(PHI3) 
        DEALLOCATE(PHI2)
!
!       =================================================================
!       ==  test schroedinger equation for aephi and psphi             ==
!       =================================================================
        ALLOCATE(PHI2(NR,NAUG))
        ALLOCATE(PHI3(NR,NAUG))
!       == phi2 and phi3 must vanish if everything is ok
        DO IB=1,NAUG
          PHI2(:,IB)=TAEPHI(:,IB,L+1) &
      &             +(AEPOT(:)*Y0-EGRID(IB,L+1))*AEPHI(:,IB,L+1)
          PHI3(:,IB)=TPSPHI(:,IB,L+1) &
      &             +(PSPOT(:)*Y0-EGRID(IB,L+1))*PSPHI(:,IB,L+1) &
      &             -PRO(:,IB,L+1)
          DO JB=1,NAUG
            AUX(:)=PRO(:,IB,L+1)*PSPHI(:,JB,L+1)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            PROPSIBAR(IB,JB,L+1)=SVAR
          ENDDO
        ENDDO
!       == test if [t+vtilde-e+|pro>(dh-edo)<pro|]|phi>=0
        WRITE(NFILO,FMT='("TEST SCHROEDINGER EQUATION FOR BARE WAVE FUNCTIONS")')
        WRITE(NFILO,FMT='("....THE FOLLOWING VALUES MUST VANISH")')
        WRITE(NFILO,FMT='(t11,"[T+V-E]|PHI>",t31,"[T+PSV-E]|PSPHI>-|PRO>")') 
        DO IB=1,NAUG
          svar1=maxval(abs(phi2(:,ib)))
          svar2=maxval(abs(phi3(:,ib)))
          write(nfilo,fmt='(i5,3e20.10)')ib,svar1,svar2
        ENDDO
        write(nfilo,fmt='("(dH-edO<pro|pspsi>+<pro|pspsi> ... must vanish")')
        do ib=1,naug
          write(nfilo,fmt='(20f20.10)')DH(:,IB,L+1)-EGRID(IB,L+1)*DO(:,IB,L+1)+PROPSIBAR(IB,:,L+1)
        enddo
        write(nfilo,fmt='("<pro|pspsi>")')
        do ib=1,naug
          write(nfilo,fmt='(20f20.10)')propsibar(ib,:,l+1)
        enddo
        DEALLOCATE(PHI2)
        DEALLOCATE(PHI3)
        WRITE(NFILO,FMT='(72("=")/72("="),T20,"TESTS DONE"/72("="))')
!
!       ======================================================================
!       == TRUNCATE PARTIAL WAVES TO AVOID NUMERICAL ERRORS WITH EXPONENTIAL
!       == INCREASE
!       ======================================================================
        uphi(ircut:,:,l+1)=0.d0
        tuphi(ircut:,:,l+1)=0.d0
        aephi(ircut:,:,l+1)=0.d0
        taephi(ircut:,:,l+1)=0.d0
        psphi(ircut:,:,l+1)=0.d0
        tpsphi(ircut:,:,l+1)=0.d0
        pro(ircut:,:,l+1)=0.d0
!
!       ======================================================================
!       == BARE QUANTITIES FINISHED. WRITE RESULT...                        ==
!       ======================================================================
!       == bareaugment contains ae,nodeless, and ps partial waves
        ALLOCATE(WORK2D(NR,3*NAUG))
        WORK2D(:,1:NAUG)=AEPHI(:,:,L+1)
        WORK2D(:,NAUG+1:2*NAUG)=UPHI(:,:,L+1)
        WORK2D(:,2*NAUG+1:3*NAUG)=PSPHI(:,:,L+1)
        CALL WRITEPHI('BAREPARTIALWAVES'//TRIM(EXT),GID,NR,3*NAUG,WORK2D)
        DEALLOCATE(WORK2D)
        CALL WRITEPHI('BAREPRO'//TRIM(EXT),GID,NR,NAUG,PRO(:,:,L+1))
!
!       =================================================================
!       ==  ENSURE BIORTHOGONALIZATION CONDITION                       ==
!       =================================================================
PSPHISAVE(:,:,L+1)=PSPHI(:,:,L+1)
        CALL BIORTHOMATRICES(GID,NR,L,NAUG,PSPHI(:,:,L+1),PRO(:,:,L+1) &
       &                    ,TRANSPHI(:,:,L+1),TRANSPRO(:,:,L+1))
        PRO(:,:,L+1)   =MATMUL(PRO(:,:,L+1),TRANSPRO(:,:,L+1))
        PSPHI(:,:,L+1) =MATMUL(PSPHI(:,:,L+1),TRANSPHI(:,:,L+1))
        TPSPHI(:,:,L+1)=MATMUL(TPSPHI(:,:,L+1),TRANSPHI(:,:,L+1))
        AEPHI(:,:,L+1) =MATMUL(AEPHI(:,:,L+1),TRANSPHI(:,:,L+1))
        TAEPHI(:,:,L+1)=MATMUL(TAEPHI(:,:,L+1),TRANSPHI(:,:,L+1))
        UPHI(:,:,L+1) =MATMUL(UPHI(:,:,L+1),TRANSPHI(:,:,L+1))
        TUPHI(:,:,L+1)=MATMUL(TUPHI(:,:,L+1),TRANSPHI(:,:,L+1))
        DT(:,:,L+1)=MATMUL(TRANSPOSE(TRANSPHI(:,:,L+1)) &
       &                  ,MATMUL(DT(:,:,L+1),TRANSPHI(:,:,L+1)))
        DO(:,:,L+1)=MATMUL(TRANSPOSE(TRANSPHI(:,:,L+1)) &
       &                  ,MATMUL(DO(:,:,L+1),TRANSPHI(:,:,L+1)))
        DH(:,:,L+1)=MATMUL(TRANSPOSE(TRANSPHI(:,:,L+1)) &
       &                  ,MATMUL(DH(:,:,L+1),TRANSPHI(:,:,L+1)))
        CALL LIB$INVERTR8(NAUG,TRANSPHI(:,:,L+1),TRANSPHIINV(:,:,L+1))
        CALL LIB$INVERTR8(NAUG,TRANSPRO(:,:,L+1),TRANSPROINV(:,:,L+1))

WRITE(*,FMT='("1C-KINETIC ENERGY")')
DO IB=1,NAUG
  WRITE(*,FMT='(20E20.5)')DT(IB,:,L+1)
ENDDO
WRITE(*,FMT='("1C-OVERLAP")')
DO IB=1,NAUG
  WRITE(*,FMT='(20E20.5)')DO(IB,:,L+1)
ENDDO
print*,'write files'
        ALLOCATE(WORK2D(NR,3*NAUG))
        WORK2D(:,1:NAUG)=AEPHI(:,:,L+1)
        WORK2D(:,NAUG+1:2*NAUG)=PSPHI(:,:,L+1)
        WORK2D(:,2*naug+1:3*NAUG)=PRO(:,:,L+1)
        CALL WRITEPHI('AUGMENT'//TRIM(EXT),GID,NR,3*NAUG,WORK2D)
        DEALLOCATE(WORK2D)
        CALL WRITEPHI('AEPHI'//TRIM(EXT),GID,NR,NAUG,AEPHI(:,:,L+1))
        CALL WRITEPHI('AEPHI_T'//TRIM(EXT),GID,NR,NAUG,TAEPHI(:,:,L+1))
        CALL WRITEPHI('PSPHI'//TRIM(EXT),GID,NR,NAUG,PSPHI(:,:,L+1))
        CALL WRITEPHI('PSPHI_T'//TRIM(EXT),GID,NR,NAUG,TPSPHI(:,:,L+1))
        CALL WRITEPHI('PRO'//TRIM(EXT),GID,NR,NAUG,PRO(:,:,L+1))
!
!       =================================================================
!       ==  TEST SCATTERING PROPERTIES                                 ==
!       =================================================================
!!$print*,'test1'
!!$DO IB=1,NAUG
!!$!   -- <P|PHIBAR>=<P|PHI>/TPHI=1/TPHI
!!$  DO JB=1,NAUG
!!$    AUX=PRO(:,JB,L+1)*PSPHISAVE(:,IB,L+1)*R(:)**2
!!$    CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$!    PRINT*,'SVAR ',IB,JB,SVAR,TRANSPHIINV(JB,IB,L+1)
!!$  ENDDO 
!!$  PRINT*,'SVEC        ',IB,MATMUL(TRANSPRO(:,:,L+1),MATMUL(DH(:,:,L+1)-EGRID(IB,L+1)*DO(:,:,L+1),TRANSPHIINV(:,IB,L+1)))
!!$ENDDO

print*,'test2'
       IF(L.EQ.0) THEN
         E=-9.D0
         AUX(:)=0.D0
         ALLOCATE(PHI2(NR,1))
         CALL PAWDER(GID,NR,L,E,PSPOT,NAUG,PRO(:,:,L+1),DH(:,:,L+1) &
      &           ,DO(:,:,L+1),AUX,PHI2(:,1))
         CALL WRITEPHI('XXX.DAT',GID,NR,1,PHI2)
         DEALLOCATE(PHI2)
       END IF
!      STOP 'IN TEST OK'
!
!       =================================================================
!       ==  TEST SCATTERING PROPERTIES                                 ==
!       =================================================================
        CALL TESTPHASESHIFT(GID,NR,L,AEPOT,PSPOT,NAUG,PRO(:,:,L+1) &
       &                        ,DH(:,:,L+1),DO(:,:,L+1))
      ENDDO
print*,'================loop done========================'
!
!     ===================================================================
!     ==  UNSCREENING                                                  ==
!     ===================================================================
      CALL UNSCREEN(GID,NR,AEZ,AERHOC+AERHOV,PSRHOC+PSRHOV,PSPOT,RCSM,VADD)
!     == TEST
      allocate(work2d(nr,3))
      work2d(:,1)=aepot*y0
      work2d(:,2)=pspot*y0
      work2d(:,3)=(pspot-vadd)*y0
      CALL WRITEPHI('POTENTIALS.DAT',GID,NR,3,WORK2D)
      DEALLOCATE(WORK2D)
      CALL WRITEPHI('PSRHOV.DAT',GID,NR,1,PSRHOV)
      CALL WRITECOMPARE('AERHO-PSRHO.DAT',GID,NR,1,AERHOV+AERHOC,PSRHOV+PSRHOC)
!
!     ===================================================================
!     ==  CALCULATE TOTAL ENERGY                                       ==
!     ===================================================================
      CALL PSETOTSHELL(GID,NR,AEZ,PSRHOC,AERHOC,PSPOT &
     &                      ,NB,NC,LMAX,NAUG,NPRO,LOFI,FOFI,EOFI &
     &                      ,PRO,AEPHI,PSPHI,RCSM,DT,DH,DO,VADD)
!
!     ===================================================================
!     ==  WRITE TO FILE                                                ==
!     ===================================================================
      CALL LINKEDLIST$NEW(LL_STP)
!
!     == GRID  ======================================
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      CALL LINKEDLIST$SELECT(LL_STP,'GRID')
      IF(TLOG) THEN
        CALL LINKEDLIST$SET(LL_STP,'TYPE',0,'SHLOG')
        CALL LINKEDLIST$SET(LL_STP,'R1',0,R1)
        CALL LINKEDLIST$SET(LL_STP,'DEX',0,DEX)
        CALL LINKEDLIST$SET(LL_STP,'NR',0,NR)
      ELSE
        CALL RADIAL$NEW('LOG',GID2)
        CALL RADIAL$SETR8(GID2,'R1',1.056D-4)
        CALL RADIAL$SETR8(GID2,'DEX',5.D-2)
        NR2=250
        CALL RADIAL$SETI4(GID2,'NR',NR2)
        CALL LINKEDLIST$SET(LL_STP,'TYPE',0,'LOG')
        CALL LINKEDLIST$SET(LL_STP,'R1',0,1.056D-4)
        CALL LINKEDLIST$SET(LL_STP,'DEX',0,5.D-2)
        CALL LINKEDLIST$SET(LL_STP,'NR',0,NR2)
        DEALLOCATE(AUX)
        ALLOCATE(AUX(NR2))
      END IF
!
!     == GENERIC ======================================
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      CALL LINKEDLIST$SELECT(LL_STP,'GENERIC')
      CALL LINKEDLIST$SET(LL_STP,'IDFT',0,IDFT)
      CALL LINKEDLIST$SET(LL_STP,'AEZ',0,AEZ)
!
!     == AUGMENTATION ======================================
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      CALL LINKEDLIST$SELECT(LL_STP,'AESCF')
      CALL LINKEDLIST$SET(LL_STP,'NB',0,NB)
      CALL LINKEDLIST$SET(LL_STP,'NC',0,NC)
      CALL LINKEDLIST$SET(LL_STP,'L',0,LOFI)
      CALL LINKEDLIST$SET(LL_STP,'NN',0,NNOFI)
      CALL LINKEDLIST$SET(LL_STP,'OCC',0,FOFI)
      CALL LINKEDLIST$SET(LL_STP,'E',0,EOFI)
      IF(TLOG) THEN
        CALL LINKEDLIST$SET(LL_STP,'AEPOT',0,AEPOT)
!!$print*,'aepot ',aepot(1:nr:10)
!!$allocate(g(nr))
!!$svar=-18.d0
!!$g(:)=0.d0
!!$CALL BOUNDSTATE(GID,NR,0,0,DREL,G,0,aePOT,svar,aePHI(:,1,1))
!!$print*,'svar',svar
!!$stop
        DO IB=1,NB
          CALL LINKEDLIST$SET(LL_STP,'PHI-NODELESS',0,UOFI(:,IB))
          CALL LINKEDLIST$SET(LL_STP,'TPHI-NODELESS',0,TUOFI(:,IB))
        ENDDO
      ELSE
        CALL RADIAL$CHANGEGRID(GID,NR,AEPOT,GID2,NR2,AUX)
        CALL LINKEDLIST$SET(LL_STP,'AEPOT',0,AUX)
        DO IB=1,NB
          CALL RADIAL$CHANGEGRID(GID,NR,UOFI(:,IB),GID2,NR2,AUX)
          CALL LINKEDLIST$SET(LL_STP,'PHI-NODELESS',0,AUX)
          CALL RADIAL$CHANGEGRID(GID,NR,TUOFI(:,IB),GID2,NR2,AUX)
          CALL LINKEDLIST$SET(LL_STP,'TPHI-NODELESS',0,AUX)
        ENDDO
      END IF
!
!     == AUGMENTATION ======================================
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION')
      PSZ=AEZ
      DO IB=1,NC
        PSZ=PSZ-FOFI(IB)
      ENDDO
      CALL LINKEDLIST$SET(LL_STP,'PSZ',0,PSZ)
      LNX=SUM(NPRO(:))
      CALL LINKEDLIST$SET(LL_STP,'LNX',0,LNX)
      ALLOCATE(LOX(LNX))
      ISVAR=0
      DO L=0,LMAX
        DO IB=1,NPRO(L+1)
          ISVAR=ISVAR+1
          LOX(ISVAR)=L
          IF(TLOG) THEN
             AUX(:)=PSPHI(:,IB,L+1)
             CALL LINKEDLIST$SET(LL_STP,'PSPHI',-1,AUX)
             AUX(:)=AEPHI(:,IB,L+1)
             CALL LINKEDLIST$SET(LL_STP,'AEPHI',-1,AUX)
             AUX(:)=PRO(:,IB,L+1)
             CALL LINKEDLIST$SET(LL_STP,'PRO',-1,AUX)
             AUX(:)=UPHI(:,IB,L+1)
             CALL LINKEDLIST$SET(LL_STP,'NDLSPHI',-1,AUX)
             AUX(:)=TUPHI(:,IB,L+1)
             CALL LINKEDLIST$SET(LL_STP,'NDLSTPHI',-1,AUX)
!            CALL LINKEDLIST$SET(LL_STP,'PSPHI',-1,PSPHI(:,IB,LL))
!            CALL LINKEDLIST$SET(LL_STP,'AEPHI',-1,AEPHI(:,IB,LL))
!            CALL LINKEDLIST$SET(LL_STP,'PRO',-1,PRO(:,IB,LL))
!            CALL LINKEDLIST$SET(LL_STP,'NDLSPHI',-1,UPHI(:,IB,LL))
!            CALL LINKEDLIST$SET(LL_STP,'NDLSTPHI',-1,TUPHI(:,IB,LL))
          ELSE
            CALL RADIAL$CHANGEGRID(GID,NR,PSPHI(:,IB,L+1),GID2,NR2,AUX)
            CALL LINKEDLIST$SET(LL_STP,'PSPHI',-1,AUX)
            CALL RADIAL$CHANGEGRID(GID,NR,AEPHI(:,IB,L+1),GID2,NR2,AUX)
            CALL LINKEDLIST$SET(LL_STP,'AEPHI',-1,AUX)
            CALL RADIAL$CHANGEGRID(GID,NR,PRO(:,IB,L+1),GID2,NR2,AUX)
            CALL LINKEDLIST$SET(LL_STP,'PRO',-1,AUX)
            CALL RADIAL$CHANGEGRID(GID,NR,UPHI(:,IB,L+1),GID2,NR2,AUX)
            CALL LINKEDLIST$SET(LL_STP,'NDLSPHI',-1,AUX)
            CALL RADIAL$CHANGEGRID(GID,NR,TUPHI(:,IB,L+1),GID2,NR2,AUX)
            CALL LINKEDLIST$SET(LL_STP,'NDLSTPHI',-1,AUX)
          END IF
        ENDDO
      ENDDO
      CALL LINKEDLIST$SET(LL_STP,'LOX',0,LOX)
!
!     == SET MATRIX ELEMENTS ===========================================
      ALLOCATE(WORK1(LNX*LNX))
      ALLOCATE(WORK2(LNX*LNX))
      ALLOCATE(WORK3(LNX*LNX))
      ISVAR=0
      I2=0
      DO L2=0,LMAX
        DO IB2=1,NPRO(L2+1)
          I2=I2+1
          I1=0
          DO L1=0,LMAX
            DO IB1=1,NPRO(L1+1)
              I1=I1+1
              ISVAR=ISVAR+1
              IF(L1.NE.L2) THEN
                WORK1(ISVAR)=0.D0
                WORK2(ISVAR)=0.D0
                WORK3(ISVAR)=0.D0
                CYCLE
              END IF
              WORK1(ISVAR)=DT(IB1,IB2,L1+1)
              WORK2(ISVAR)=DO(IB1,IB2,L1+1)
              WORK3(ISVAR)=DH(IB1,IB2,L1+1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL LINKEDLIST$SET(LL_STP,'DT',0,WORK1)
      CALL LINKEDLIST$SET(LL_STP,'DO',0,WORK2)
      CALL LINKEDLIST$SET(LL_STP,'DH',0,WORK3)
      DEALLOCATE(WORK1)
      DEALLOCATE(WORK2)
      DEALLOCATE(WORK3)
      CALL LINKEDLIST$SET(LL_STP,'RCSM',0,RCSM)
      IF(TLOG) THEN
        CALL LINKEDLIST$SET(LL_STP,'PSCORE',0,PSRHOC)
        CALL LINKEDLIST$SET(LL_STP,'AECORE',0,AERHOC)
        CALL LINKEDLIST$SET(LL_STP,'PSPOT',0,PSPOT)
        CALL LINKEDLIST$SET(LL_STP,'VADD',0,VADD)
      ELSE
        CALL RADIAL$CHANGEGRID(GID,NR,PSRHOC,GID2,NR2,AUX)
        CALL LINKEDLIST$SET(LL_STP,'PSCORE',0,AUX)
        CALL RADIAL$CHANGEGRID(GID,NR,AERHOC,GID2,NR2,AUX)
        CALL LINKEDLIST$SET(LL_STP,'AECORE',0,AUX)
        CALL RADIAL$CHANGEGRID(GID,NR,PSPOT,GID2,NR2,AUX)
        CALL LINKEDLIST$SET(LL_STP,'PSPOT',0,AUX)
        CALL RADIAL$CHANGEGRID(GID,NR,VADD,GID2,NR2,AUX)
        CALL LINKEDLIST$SET(LL_STP,'VADD',0,AUX)
      END IF
!
!     ==  WRITE TO FILE ================================================
      CALL FILEHANDLER$UNIT('SETUP',NFIL)
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$WRITE(LL_STP,NFIL)
!     ==  CHECK CONTENTS ================================================
      CALL LINKEDLIST$SELECT(LL_STP,'~')
!     CALL LINKEDLIST$REPORT(LL_STP,NFILO)
      CALL FILEHANDLER$CLOSE('SETUP')

      CALL TESTSETUPREAD(ROOTNAME)
      STOP
      END
!
!     ......................................................................
      SUBROUTINE TESTSETUPREAD(ROOTNAME)
!     **                                                                  **
!     **                                                                  **
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ROOTNAME
      REAL(8)   ,ALLOCATABLE  :: R(:)
      INTEGER(4)              :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)              :: NFIL
      INTEGER(4)              :: NFILO
      LOGICAL(4)              :: TCHK
      INTEGER(4)              :: LNX
      INTEGER(4)              :: GID
      REAL(8)   ,ALLOCATABLE  :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE  :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE  :: PRO(:,:)
      INTEGER(4),ALLOCATABLE  :: LOX(:)
      REAL(8)                 :: AEZ,PSZ
      REAL(8)                 :: RCSM
      REAL(8)   ,ALLOCATABLE  :: PSCORE(:)
      REAL(8)   ,ALLOCATABLE  :: AECORE(:)
      REAL(8)   ,ALLOCATABLE  :: VADD(:)
!     **********************************************************************
PRINT*,' START READING SETUP FILE'
!
!     ==================================================================
!     == DEFINE INPUT FILE                                            ==
!     ==================================================================
      CALL FILEHANDLER$SETFILE('INPUT1',.FALSE.,TRIM(ROOTNAME)//-'.STP')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT1','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT1','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT1','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT1','FORM','FORMATTED')
!
      CALL FILEHANDLER$UNIT('INPUT1',NFIL)
      CALL SETUPREAD$NEW(NFIL,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$STOP('TESTREAD')
      END IF
      CALL SETUPREAD$GETI4('LNX',LNX)
      CALL SETUPREAD$GETI4('GID',GID)
      CALL SETUPREAD$GETI4('NR',NR)
      CALL SETUPREAD$GETR8('AEZ',AEZ)
      CALL SETUPREAD$GETR8('PSZ',PSZ)
      CALL SETUPREAD$GETR8('RCSM',RCSM)
      ALLOCATE(LOX(LNX))
      CALL SETUPREAD$GETI4A('LOX',LNX,LOX)
      ALLOCATE(PRO(NR,LNX))
      ALLOCATE(AEPHI(NR,LNX))
      ALLOCATE(PSPHI(NR,LNX))
      CALL SETUPREAD$GETR8A('PRO',NR*LNX,PRO)
      CALL SETUPREAD$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUPREAD$GETR8A('PSPHI',NR*LNX,PSPHI)
      ALLOCATE(PSCORE(NR))
      ALLOCATE(AECORE(NR))
      CALL SETUPREAD$GETR8A('PSCORE',NR,PSCORE)
      CALL SETUPREAD$GETR8A('AECORE',NR,AECORE)
!
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE UNSCREEN(GID,NR,AEZ,AERHO,PSRHO,PSPOT,RCSM,VADD)
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      REAL(8)     ,INTENT(IN) :: AEZ
      REAL(8)     ,INTENT(IN) :: AERHO(NR)
      REAL(8)     ,INTENT(IN) :: PSRHO(NR)
      REAL(8)     ,INTENT(IN) :: PSPOT(NR)
      REAL(8)     ,INTENT(IN) :: RCSM
      REAL(8)     ,INTENT(OUT):: VADD(NR)
      REAL(8)                 :: R(NR)
      REAL(8)                 :: AUX(NR),SVAR
      REAL(8)                 :: POT(NR)
      REAL(8)                 :: PI,Y0
      REAL(8)                 :: QLM
      REAL(8)                 :: ALPHA,CL
      REAL(8)                 :: GRHO(NR)
      INTEGER(4)              :: IR
      REAL(8)                 :: RH,GRHO2,VXC,VGXC,EXC,DUMMY1,DUMMY2,DUMMY3
!     ************************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!
!     ========================================================================
!     == MOMENT OF DIFFERENCE CHARGE DENSITY                                ==
!     ========================================================================
      AUX(:)=(AERHO(:)-PSRHO(:))*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,QLM)
      QLM=QLM-AEZ*Y0    ! CHARGE =QM/Y0
!
!     ========================================================================
!     == ADD COMPENSATION DENSITY AND DETERMINE ELECTROSTATIC POTENTIA      ==
!     ========================================================================
      ALPHA=1.D0/RCSM**2
      CALL GAUSSN(0,ALPHA,CL)
      SVAR=QLM*CL
      AUX(:)=PSRHO(:)+SVAR*EXP(-ALPHA*R(:)**2)
      CALL RADIAL$POISSON(GID,NR,0,AUX,POT)
!
!     ========================================================================
!     == EXCHANGE AND CORRELATION                                           ==
!     ========================================================================
      CALL RADIAL$DERIVE(GID,NR,PSRHO(:),GRHO)
      DO IR=1,NR
        RH=PSRHO(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        POT(IR)=POT(IR)+VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      CALL RADIAL$DERIVE(GID,NR,GRHO(:),AUX)
      POT(:)=POT(:)-AUX(:)
      IF(R(1).GT.1.D+10) THEN
         POT(:)=POT(:)-2.D0/R(:)*GRHO(:)
      ELSE
        POT(2:)=POT(2:)-2.D0/R(2:)*GRHO(2:)
        POT(1)=POT(1)-2.D0/R(2)*GRHO(2)
      END IF
!
!     ========================================================================
!     == VADD                                                               ==
!     ========================================================================
      VADD(:)=PSPOT(:)-POT(:)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE WRITEPHI(FILE,GID,NR,NPHI,PHI)
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI(NR,NPHI)
      INTEGER(4)              :: IR
      REAL(8)                 :: R(NR)
!     ***********************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        WRITE(100,FMT='(F15.10,2X,20(F25.15,2X))')R(IR),PHI(IR,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE WRITECOMPARE(FILE,GID,NR,NPHI,PHI1,PHI2)
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI1(NR,NPHI)
      REAL(8)     ,INTENT(IN) :: PHI2(NR,NPHI)
      INTEGER(4)              :: IR
      REAL(8)                 :: R(NR)
!     ***********************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        WRITE(100,FMT='(F15.10,2X,8(F25.15,2X))')R(IR),PHI1(IR,:),PHI2(IR,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE TESTPHASESHIFT(GID,NR,L,AEPOT,PSPOT,NAUG,PRO,DH,DO)
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID      ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L        ! ANGULAR MOMENTUM
      REAL(8)    ,INTENT(IN)     :: AEPOT(NR)  ! AE POTENTIAL (RADIAL PART ONLY)
      REAL(8)    ,INTENT(IN)     :: PSPOT(NR)  ! PS POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NAUG       ! 
      REAL(8)    ,INTENT(IN)     :: PRO(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: DH(NAUG,NAUG)
      REAL(8)    ,INTENT(IN)     :: DO(NAUG,NAUG)
      INTEGER(4) ,PARAMETER      :: NE=1000
      REAL(8)    ,PARAMETER      :: EMIN=-10.D0,EMAX=10.D0
      REAL(8)    ,PARAMETER      :: RPHASE=3.D0
      REAL(8)                    :: PSPHASE(NE,NAUG)
      REAL(8)                    :: AEPHASE(NE)
      REAL(8)                    :: PSPHI(NR)
      REAL(8)                    :: AEPHI(NR)
      REAL(8)                    :: G(NR)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: DREL(NR)
      INTEGER(4)                 :: SO
      INTEGER(4)                 :: IR,IE,NPRO
      REAL(8)                    :: E
      REAL(8)                    :: NN
      REAL(8)                    :: VAL,DER
      REAL(8)                    :: PI
!     ***********************************************************************
      PI=4.D0*DATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
      DO IE=1,NE
        E=EMIN+(EMAX-EMIN)/REAL(NE-1)*REAL(IE-1)
        G(:)=0.D0
        DO NPRO=1,NAUG
          CALL PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO(:,:NPRO) &
                    ,DH(:NPRO,:NPRO),DO(:NPRO,:NPRO),G,PSPHI)
          CALL PHASESHIFT(GID,NR,PSPHI,RPHASE,PSPHASE(IE,NPRO))
        ENDDO
!       ================================================================
        DREL(:)=0.D0
        SO=0.D0
        CALL RADIAL$SCHRODINGER(GID,NR,AEPOT,DREL,SO,G,L,E,1,AEPHI)
        CALL PHASESHIFT(GID,NR,AEPHI,RPHASE,AEPHASE(IE))
      ENDDO
!
!     ======================================================================
!     == WRITE RESULT                                                     ==
!     ======================================================================
      IF(L.EQ.0) THEN
        OPEN(100,FILE='PHASESHIFT_S.DAT')
      ELSE IF (L.EQ.1) THEN
        OPEN(100,FILE='PHASESHIFT_P.DAT')
      ELSE IF (L.EQ.2) THEN
        OPEN(100,FILE='PHASESHIFT_D.DAT')
      ELSE IF (L.EQ.3) THEN
        OPEN(100,FILE='PHASESHIFT_F.DAT')
      END IF
      DO IE=1,NE
        E=EMIN+(EMAX-EMIN)/REAL(NE-1)*REAL(IE-1)
        WRITE(100,FMT='(F20.10,2X,F20.10,2X,20F20.10)')E,AEPHASE(IE),PSPHASE(IE,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE CONSTRUCTPSPHI(GID,NR,RC,L,NB,PSPHI,TPSPHI)
!     **                                                                  **
!     **  PSEUDIZES A WAVE FUNCTION BY MATCHING SPHERICAL BESSEL FUNCTIONS**
!     **                                                                  **
!     **  FIRST WE FIND THE VALUES K_I SO THAT THE SPHERICAL BESSEL       **
!     **  FUNCTION J_L(K_I*R) MATCHES THE LOGARITHMIC DERIVATIVE          **
!     **  OF THE FUNCTION TO BE PSEUDIZED.                                **
!     **                                                                  **
!     **  THE FIRST TWO ARE THEN MATCHED SO THAT THE PHI AND TPHI         **
!     **  ARE CONTINUOUS. THE DERIVATIVE OF PHI IS AUTOMATICALLY          **
!     **  CONTINUOUS THROUGH THE CORRECT CHOICE OF THE K_I.               **
!     **                                                                  **
!     **  ATTENTION! SINCE WE MATCH TO TAYLOR COEFFICIENTS, ONE SHOULD    **
!     **  ALSO MATCH ENERGY DERIVATIVES OF THE SPHERICAL BESSEL FUNCTIONS.**
!     **                                                                  **
!     **  NOTE THAT THE TAYLOR EXPANSION OF THE WAVE FUNCTIONS IN ENERGY  **
!     **  BEHAVE DIFFERENT FROM THE ENERGY DEPENDENT SOLUTION.            **
!     **  THE TAYLOR EXPANSION CAN CREATE TWO NODES AT ONCE, WHICH        **
!     **  CAUSES DIFFICULTIES.                                            **
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: L
      INTEGER(4) ,INTENT(IN)     :: NB
      REAL(8)    ,INTENT(IN)     :: RC
      REAL(8)    ,INTENT(INOUT)  :: PSPHI(NR,NB)
      REAL(8)    ,INTENT(INOUT)  :: TPSPHI(NR,NB)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: SVAR
      INTEGER(4)                 :: IR,I,ITER
      INTEGER(4),PARAMETER       :: NITER=100
      REAL(8)                    :: TOL=1.D-8
      REAL(8)                    :: PI
      REAL(8)                    :: VAL,DER,VALT,DERT
      REAL(8)                    :: JL(NR),DJLDR(NR),PHASE(NR)
      REAL(8)                    :: NN
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,Y0,DX,XM,YM
      INTEGER(4),PARAMETER       :: NZEROS=3
      REAL(8)                    :: KI(NZEROS)
      REAL(8)                    :: TJL(NR)
      REAL(8)                    :: A11,A12,A21,A22,DET
      REAL(8)                    :: AINV11,AINV12,AINV21,AINV22
      REAL(8)                    :: V1,V2,V3,C1,C2,C3,C1A,C2A,C1B,C2B
      REAL(8)                    :: DTPHI,DTJ1,DTJ2,DTJ3
      REAL(8)                    :: SVAR1,SVAR2,FAC
      INTEGER(4)                 :: IB,JB
      REAL(8)                    :: SUPPHI(NR,NB)
      REAL(8)                    :: SUPTPHI(NR,NB)
      REAL(8)                    :: SHIFT
      LOGICAL(4),PARAMETER       :: TDIFFERENTIABLELAPLACIAN=.FALSE.
!     ***********************************************************************
      PI=4.D0*DATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
!     =======================================================================
!     == MAP POWERSERIES EXPANSION ONTO SUPPORT ENERGY GRID                ==
!     =======================================================================
      SUPPHI(:,:)=PSPHI(:,:)
      SUPTPHI(:,:)=TPSPHI(:,:)
!
!     =======================================================================
!     == CALCULATE BESSEL FUNCTION                                         ==
!     =======================================================================
      DO IR=1,NR
        CALL SPECIALFUNCTION$BESSEL(L,R(IR),JL(IR)) 
      ENDDO
!     == 2*TJL-JL=0 =======================================================
      TJL(:)=0.5D0*JL(:)   
!
!     =======================================================================
!     == CALCULATE PHASESHIFT OF BESSEL FUNCTION                           ==
!     =======================================================================
      CALL RADIAL$VERLETD1(GID,NR,JL,DJLDR)
      NN=0.D0
      DO IR=2,NR
        IF(JL(IR)*JL(IR-1).LT.0.D0)NN=NN+1.D0
        PHASE(IR)=0.5D0-ATAN(R(IR)/RC*DJLDR(IR)/JL(IR))/PI+NN
      ENDDO
!     == AVOID GLITCH
      PHASE(1)=PHASE(2)
      PHASE(NR)=PHASE(NR-1)
!
!     =======================================================================
!     == PSEUDIZE EVERY FUNCTION ON THE SUPPORT ENERGY GRID                ==
!     =======================================================================
      DO IB=1,NB
!       =======================================================================
!       == DETERMINE K-VALUES FOR WHICH THE LOGARITHIC DERIVATIVE MATCHES    ==
!       =======================================================================
        CALL PHASESHIFT(GID,NR,SUPPHI,RC,SHIFT)
        IF(PHASE(1)-SHIFT.GT.0.D0) THEN
PRINT*,'WARNING SHIFT WAS REQUIRED ',L,IB
           SHIFT=SHIFT+REAL(INT(PHASE(1)-SHIFT+1.D0))
        END IF
!!$OPEN(100,FILE='PHASE.DAT')
!!$REWIND(100)
!!$DO IR=1,NR
!!$  WRITE(100,FMT='(20F20.8)')R(IR),PHASE(IR)
!!$ENDDO
!!$CLOSE(100)
PRINT*,'SHIFT ',IB,SHIFT,NN
!       =======================================================================
!       == SEARCH ZERO OF PHASE USING BISECTION                              ==
!       =======================================================================
         DO I=1,NZEROS
          ISTART=1
          X0=0.4D0+1.D-3
          DX=0.1D0
          CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
          DO ITER=1,NITER
            CALL RADIAL$VALUE(GID,NR,PHASE,X0,Y0)
            Y0=Y0-SHIFT
            CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
            IF(ABS(DX).LT.TOL) EXIT
          ENDDO
          KI(I)=X0/RC
          IF(KI(I).LT.0.D0) THEN
            CALL ERROR$MSG('KI NEGATIVE')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$I4VAL('IB',IB)
            CALL ERROR$I4VAL('I',I)
            CALL ERROR$R8VAL('KI(I)',KI(I))
            CALL ERROR$STOP('CONSTRUCTPSPHI')
          END IF
          SHIFT=SHIFT+1.D0
        ENDDO
PRINT*,'KI ',KI

!       =======================================================================
!       == CONSTRUCT NEW PSEUDO WAVE FUNCTION                                ==
!       =======================================================================
        CALL RADIAL$VALUE(GID,NR,SUPTPHI(:,IB),RC,VALT)
        CALL RADIAL$DERIVATIVE(GID,NR,SUPTPHI(:,IB),RC,DERT)     
        IF(KI(3)*RC.GT.R(NR-2)) THEN
          CALL ERROR$MSG('KI(3)*RC>RX')
          CALL ERROR$R8VAL('KI(3)*RC',KI(3)*RC)
          CALL ERROR$R8VAL('KI(3)',KI(3))
          CALL ERROR$R8VAL('RC',RC)
          CALL ERROR$R8VAL('R(NR-2)',R(NR-2))
          CALL ERROR$STOP('CONSTRUCTPSPHI')
        END IF
!       == MATCH K1 AND K2
        CALL RADIAL$VALUE(GID,NR,JL,KI(1)*RC,A11)
        CALL RADIAL$VALUE(GID,NR,JL,KI(2)*RC,A12)
        CALL RADIAL$VALUE(GID,NR,TJL,KI(1)*RC,A21)
        A21=A21*KI(1)**2
        CALL RADIAL$VALUE(GID,NR,TJL,KI(2)*RC,A22)
        A22=A22*KI(2)**2
        DET=A11*A22-A12*A21
        AINV11=A22/DET
        AINV22=A11/DET
        AINV12=-A12/DET
        AINV21=-A21/DET
        V1=VAL
        V2=VALT
        C1A=AINV11*V1+AINV12*V2
        C2A=AINV21*V1+AINV22*V2
        IF(TDIFFERENTIABLELAPLACIAN) THEN
!         ==  USE A THIRD BESSELFUNCTION TO OBTAIN A DIFFERENTIABLE KINETIC 
!         ==  ENERGY DENSITY. THIS DOES NOT SEEM TO BE A GOOD IDEA
!         ==  BECAUSE IT DETERIORATES THE SCATTERING PROPERTIES.
!         == MATCH K2 AND K3
          CALL RADIAL$VALUE(GID,NR,JL,KI(2)*RC,A11)
          CALL RADIAL$VALUE(GID,NR,JL,KI(3)*RC,A12)
          CALL RADIAL$VALUE(GID,NR,TJL,KI(2)*RC,A21)
          A21=A21*KI(2)**2
          CALL RADIAL$VALUE(GID,NR,TJL,KI(3)*RC,A22)
          A22=A22*KI(3)**2
          DET=A11*A22-A12*A21
          AINV11=A22/DET
          AINV22=A11/DET
          AINV12=-A12/DET
          AINV21=-A21/DET
          V1=VAL
          V2=VALT
          C1B=AINV11*V1+AINV12*V2
          C2B=AINV21*V1+AINV22*V2
!         ======================================================
!         ==A SUPERPOSITION OF THE TWO SOLUTIONS WITH FACTORS ==
!         == THAT ADD UP TO ONE, DO NOT AFFECT VALUE OR       ==
!         == DERIVATIVE. USE THIS TO MATCH DERIVATIVE OF      ==
!         == KINETIC ENERGY                                   ==
!         ======================================================
          CALL RADIAL$DERIVATIVE(GID,NR,TJL,KI(1)*RC,DTJ1)      
          DTJ1=KI(1)**3*DTJ1
          CALL RADIAL$DERIVATIVE(GID,NR,TJL,KI(2)*RC,DTJ2)      
          DTJ2=KI(2)**3*DTJ2
          CALL RADIAL$DERIVATIVE(GID,NR,TJL,KI(3)*RC,DTJ3)      
          DTJ3=KI(3)**3*DTJ3
          SVAR1=DTJ1*C1A+DTJ2*C2A
          SVAR2=DTJ2*C1B+DTJ3*C2B
          FAC=(DERT-SVAR2)/(SVAR1-SVAR2)
          C1=C1A*FAC
          C2=C2A*FAC+C1B*(1.D0-FAC)
          C3=C2B*(1.D0-FAC)
        ELSE
          C1=C1A
          C2=C2A
          C3=0.D0
        END IF
        DO IR=1,NR
          IF(R(IR).GT.RC) EXIT
          CALL RADIAL$VALUE(GID,NR,JL,KI(1)*R(IR),V1)
          CALL RADIAL$VALUE(GID,NR,JL,KI(2)*R(IR),V2)
          CALL RADIAL$VALUE(GID,NR,JL,KI(3)*R(IR),V3)
          PSPHI(IR,IB)=V1*C1+V2*C2+V3*C3
          CALL RADIAL$VALUE(GID,NR,TJL,KI(1)*R(IR),V1)
          V1=V1*KI(1)**2
          CALL RADIAL$VALUE(GID,NR,TJL,KI(2)*R(IR),V2)
          V2=V2*KI(2)**2
          CALL RADIAL$VALUE(GID,NR,TJL,KI(3)*R(IR),V3)
          V3=V3*KI(3)**2
          TPSPHI(IR,IB)=V1*C1+V2*C2+V3*C3
        ENDDO
        PSPHI(1,IB)=PSPHI(2,IB)
        TPSPHI(1,IB)=TPSPHI(2,IB)
      ENDDO
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE CONSTRUCTPSPHI1(GID,NR,RC,POW,PSPOT,L,NB,E,PSPHI,TPSPHI)
!     **                                                                  **
!     **  PSEUDIZES A WAVE FUNCTION                                       **
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: L
      INTEGER(4) ,INTENT(IN)     :: NB
      REAL(8)    ,INTENT(IN)     :: RC
      integer(4) ,INTENT(IN)     :: POW
      REAL(8)    ,INTENT(IN)     :: PSPOT(NR)
      REAL(8)    ,INTENT(IN)     :: E(NB)
      REAL(8)    ,INTENT(INOUT)  :: PSPHI(NR,NB)
      REAL(8)    ,INTENT(INOUT)  :: TPSPHI(NR,NB)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: PHASE,PHASE1
      INTEGER(4)                 :: ITER,IB,IR
      INTEGER(4) ,PARAMETER      :: NITER=300
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,F0,XM,FM,DX
      REAL(8)                    :: VAL0
      REAL(8)                    :: PHI(NR)
      REAL(8)                    :: POT(NR)
      REAL(8)    ,PARAMETER      :: TOL=1.D-10
      REAL(8)                    :: DREL(NR),G(NR)
      INTEGER(4)                 :: IR0
      REAL(8)                    :: SVAR,SVAR1
      REAL(8)                    :: RCMATCH
!     ***********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      RCMATCH=2.D0*RC  !MATCH SUFFICIENTLY FAR OUT TO AVOID PROBLEMS WITH
                       !INHOMOGENEITY OF THE NODELESS EQUATION
      DO IR=1,NR
        IF(R(IR).GT.RCMATCH)  THEN
          IR0=IR-1
          EXIT
        END IF
      ENDDO
!     =======================================================================
!     == MAP POWERSERIES EXPANSION ONTO SUPPORT ENERGY GRID                ==
!     =======================================================================
      DO IB=1,NB
print*,'constructpsphi1',ib
!
!       =====================================================================
!       ==  DETERMINE PHASESHIFT OF PSPHI                                  ==
!       =====================================================================
        CALL PHASESHIFT(GID,NR,PSPHI(:,IB),RCMATCH,PHASE)
!
!       =====================================================================
!       ==  FIND POTENTIAL THAT PRODUCES CORRECT PHASESHIFT                ==
!       =====================================================================
        ISTART=-1
        X0=pspot(1)
        DX=0.1D0
        CALL BISEC(ISTART,IBI,X0,F0,DX,XM,FM)
        DO ITER=1,NITER
          VAL0=X0
          CALL PSEUDIZE(GID,NR,POW,VAL0,RC,PSPOT,POT)
          DREL=0.D0
          G=0.D0
          CALL RADIAL$SCHRODINGER(GID,NR,POT,DREL,0,G,L,E(IB),1,PHI(:))
          CALL PHASESHIFT(GID,NR,PHI,RCMATCH,PHASE1)
          F0=PHASE1-PHASE
          CALL BISEC(ISTART,IBI,X0,F0,DX,XM,FM)
          IF(ABS(DX).LT.TOL) EXIT
        ENDDO
        IF(DX.GT.TOL) THEN
          CALL ERROR$MSG('LOOP NOT CONVERGED')
          CALL ERROR$STOP('CONSTRUCTPSPHI1')
        END IF
!       == SCALE 
        CALL RADIAL$VALUE(GID,NR,PHI,RCMATCH,SVAR)
        CALL RADIAL$VALUE(GID,NR,PSPHI(:,IB),RCMATCH,SVAR1)
        PHI(:)=PHI(:)*SVAR1/SVAR
!       == ADJUST PSPHI
        PSPHI(1:IR0,IB)=PHI(1:IR0)
        TPSPHI(1:IR0,IB)=(E(IB)-POT(1:IR0)*Y0)*PHI(1:IR0)
      ENDDO
      RETURN
      END
!
!     .......................................................................
      SUBROUTINE PHASESHIFT(GID,NR,PHI,RC,PHASE)
!     **                                                                   **
!     **  RETURNS THE PHASE SHIFT DEFINED AS                               **
!     **        PHASE=0.5-ARCTAN(DPHI/DR/PHI)+#(NODES)                     **
!     **  PHASE IS A MONOTONICALLY INCREASING FUNCTION WITH ENERGY         **
!     **  IF PHI IS OBTAINED FROM A LOCAL OR SEMILOCAL POTETIAL            **
!     **                                                                   **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: RC
      REAL(8)   ,INTENT(IN) :: PHI(NR)
      REAL(8)   ,INTENT(OUT):: PHASE
      REAL(8)               :: PI
      REAL(8)               :: R(NR)
      REAL(8)               :: VAL,DER
      INTEGER(4)            :: IR
!     ***********************************************************************   
      PI=4.D0*DATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$VALUE(GID,NR,PHI,RC,VAL)
      CALL RADIAL$DERIVATIVE(GID,NR,PHI,RC,DER)
      PHASE=0.5D0-ATAN(DER/VAL)/PI
      DO IR=3,NR  ! LEAVE OUT FIRST POINT WHICH IS OFTEN ZERO
        IF(R(IR).GT.RC) THEN
          IF(PHI(IR-1)*VAL.LT.0.D0)PHASE=PHASE+1.D0
          EXIT
        END IF
        IF(PHI(IR)*PHI(IR-1).LT.0.D0)PHASE=PHASE+1.D0
      ENDDO
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE MAKEDATH(GID,NR,AEPOT,PSPOT,L,NAUG,AEPHI,PSPHI,DT,DH)
!     **                                                                  **
!     **  ONE-CENTER HAMILTONIAN                                          **
!     **  DH=DT+<AEPHI|V|AEPHI>-<PSPHI|VTILDE|PSPHI>                      **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID            ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR             ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L              ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: NAUG           ! #(PARTIAL WAVES)
      REAL(8)    ,INTENT(IN)     :: AEPOT(NR)      ! AE-POTENTIAL
      REAL(8)    ,INTENT(IN)     :: PSPOT(NR)      ! PS-POTENTIAL
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,NAUG) ! PS-PARTIAL WAVE
      REAL(8)    ,INTENT(IN)     :: AEPHI(NR,NAUG) ! AE PARTIAL WAVE
      REAL(8)    ,INTENT(IN)     :: DT(NAUG,NAUG)  ! 1C KINETIC ENERGY
      REAL(8)    ,INTENT(OUT)    :: DH(NAUG,NAUG)  ! 1C-HAMILTONIAN
      REAL(8)                    :: R(NR)
      REAL(8)                    :: AUX(NR)
      INTEGER(4)                 :: IB,JB
      REAL(8)                    :: SVAR
      REAL(8)                    :: PI,Y0
!     ***********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      DO IB=1,NAUG
        DO JB=1,IB
          AUX(:)=AEPHI(:,IB)*AEPOT(:)*AEPHI(:,JB)-PSPHI(:,IB)*PSPOT(:)*PSPHI(:,JB)
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          SVAR=0.5D0*(DT(IB,JB)+DT(JB,IB))+SVAR*Y0
          DH(IB,JB)=SVAR
          DH(JB,IB)=SVAR
        ENDDO
      ENDDO
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE DTDO(GID,NR,NAUG,AEPHI,TAEPHI,PSPHI,TPSPHI,DT,DO)
!     **                                                                  **
!     **  ONE-CENTER KINETIC ENERGY AND ONE-CENTER OVERLAP                **
!     **     DT=<AEPHI|T|AEPHI>-<PSPHI|T|PSPHI>                           **
!     **     DO=<AEPHI|AEPHI>-<PSPHI|PSPHI>                               **
!     **  FROM THE PARTIAL WAVES. T IS THE KINETIC ENERGY OPERATOR        **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NAUG
      REAL(8)    ,INTENT(IN)     :: AEPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: TAEPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: TPSPHI(NR,NAUG)
      REAL(8)    ,INTENT(OUT)    :: DT(NAUG,NAUG)
      REAL(8)    ,INTENT(OUT)    :: DO(NAUG,NAUG)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: R2(NR)
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: SVAR
      INTEGER(4)                 :: IB,JB,IR
      LOGICAL(4)                 :: TTEST=.TRUE.
!     *********************************************************************
      CALL RADIAL$R(GID,NR,R)
      R2(:)=R(:)**2
      DO IB=1,NAUG
        DO JB=1,NAUG
          AUX(:)=(AEPHI(:,IB)*TAEPHI(:,JB)-PSPHI(:,IB)*TPSPHI(:,JB))*R2(:)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          DT(IB,JB)=SVAR
          AUX(:)=(AEPHI(:,IB)*AEPHI(:,JB)-PSPHI(:,IB)*PSPHI(:,JB))*R2(:)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          DO(IB,JB)=SVAR
        ENDDO
      ENDDO
!     ======================================================================
!     ==  TEST
!     ======================================================================
      IF(TTEST) THEN
        WRITE(*,FMT='("1C-KINETIC ENERGY")')
        DO IB=1,NAUG
          WRITE(*,FMT='(20E20.5)')DT(IB,:)
        ENDDO
        WRITE(*,FMT='("1C-OVERLAP")')
        DO IB=1,NAUG
          WRITE(*,FMT='(20E20.5)')DO(IB,:)
        ENDDO
      END IF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE COREORTHOGONALIZE(GID,NR,NC,LOFI,SOFI,UOFI,TUOFI &
     &                            ,L,NAUG,UPHI,TUPHI,AEPHI,TAEPHI)
!     **                                                                  **
!     **  ORTHOGONALIZE NODE-LESS WAVE FUNCTIONS TO THE CORE              **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NC
      INTEGER(4) ,INTENT(IN)     :: LOFI(NC)
      INTEGER(4) ,INTENT(IN)     :: SOFI(NC)
      REAL(8)    ,INTENT(IN)     :: UOFI(NR,NC)
      REAL(8)    ,INTENT(IN)     :: TUOFI(NR,NC)
      INTEGER(4) ,INTENT(IN)     :: L
      INTEGER(4) ,INTENT(IN)     :: NAUG
      REAL(8)    ,INTENT(IN)     :: UPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: TUPHI(NR,NAUG)
      REAL(8)    ,INTENT(OUT)    :: AEPHI(NR,NAUG)
      REAL(8)    ,INTENT(OUT)    :: TAEPHI(NR,NAUG)
      INTEGER(4)                 :: IB,JB
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: SVAR
      INTEGER(4)                 :: NCL
      REAL(8)   ,ALLOCATABLE     :: PHIC(:,:)
      REAL(8)   ,ALLOCATABLE     :: TPHIC(:,:)
      LOGICAL(4),PARAMETER       :: TTEST=.FALSE.
!     *********************************************************************
      CALL RADIAL$R(GID,NR,R)
      AEPHI(:,:)=UPHI(:,:)
      TAEPHI(:,:)=TUPHI(:,:)
!     =======================================================================
!     == COLLECT CORE WAVE FUNCTIONS
!     =======================================================================
      NCL=0
      DO IB=1,NC
        IF(LOFI(IB).NE.L) CYCLE
        NCL=NCL+1
      ENDDO
      ALLOCATE(PHIC(NR,NCL))
      ALLOCATE(TPHIC(NR,NCL))
      NCL=0
      DO IB=1,NC
        IF(LOFI(IB).NE.L) CYCLE
        NCL=NCL+1
        PHIC(:,NCL) =UOFI(:,IB)
        TPHIC(:,NCL)=TUOFI(:,IB)
      ENDDO

!     =======================================================================
!     == ORTHOGONALIZE CORE WAVE FUNCTIONS                 
!     =======================================================================
      DO IB=1,NCL
        DO JB=1,IB-1
          AUX(:)=PHIC(:,IB)*PHIC(:,JB)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PHIC(:,IB) = PHIC(:,IB)- PHIC(:,JB)*SVAR
          TPHIC(:,IB)=TPHIC(:,IB)-TPHIC(:,JB)*SVAR
        ENDDO
        AUX(:)=(PHIC(:,IB)*R(:))**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        SVAR=1.D0/SQRT(SVAR)
        PHIC(:,IB) = PHIC(:,IB)*SVAR
        TPHIC(:,IB)=TPHIC(:,IB)*SVAR
      ENDDO
!
!     =======================================================================
!     == ORTHOGONALIZE NODE-LESS WAVE FUNCTIONS                            ==
!     =======================================================================
      DO IB=1,NAUG
        DO JB=1,NCL
!         == ORTHOGONALIZE PARTIAL WAVES ===================================
          AUX(:)=PHIC(:,JB)*AEPHI(:,IB)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)    
          AEPHI(:,IB) = AEPHI(:,IB)- PHIC(:,JB)*SVAR
          TAEPHI(:,IB)=TAEPHI(:,IB)-TPHIC(:,JB)*SVAR
        ENDDO
      ENDDO
!
      DEALLOCATE(PHIC)
      DEALLOCATE(TPHIC)
!     =======================================================================
!     == TEST
!     =======================================================================
      IF(TTEST) THEN
        DO IB=1,NAUG
          DO JB=1,NC
            IF(LOFI(JB).NE.L) CYCLE
!           == ORTHOGONALIZE PARTIAL WAVES ===================================
            AUX(:)=UOFI(:,JB)*AEPHI(:,IB)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)       
            PRINT*,'<UCORE|AEPHI> ',IB,JB,SVAR
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE BIORTHOMATRICES(GID,NR,L,NAUG,PSPHI,PRO,TRANSPHI,TRANSPRO)
!     **                                                                  **
!     ** DETERMINES THE MATRICES TRANSPHI AND TRANSPRO SUCH THAT          **
!     **     |PHI-BAR>:=|PHI>TRANSSPHI                                    **
!     **     |PRO-BAR>:=|PRO>TRANSSPRO                                    **
!     **  OBEY  <PHIBAR(I)|PROBAR(J)>=DELTA(I,J)    (KRONECKER DELTA)     **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NAUG
      INTEGER(4) ,INTENT(IN)     :: L
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: PRO(NR,NAUG)
      REAL(8)    ,INTENT(OUT)    :: TRANSPHI(NAUG,NAUG)
      REAL(8)    ,INTENT(OUT)    :: TRANSPRO(NAUG,NAUG)
      INTEGER(4)                 :: IB,JB
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: SVAR
      LOGICAL(4),PARAMETER       :: TTEST=.TRUE.
      REAL(8)                    :: PROPSI(NAUG,NAUG)
      REAL(8)                    :: TRANSPROPSI(NAUG,NAUG)
!     *********************************************************************
      CALL RADIAL$R(GID,NR,R)
!     =====================================================================
!     == CALCULATE INITIAL VIOLATION OF BIORTHOGONALITY                  ==
!     =====================================================================
      DO IB=1,NAUG
        DO JB=1,NAUG
          AUX(:)=PRO(:,IB)*PSPHI(:,JB)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)       
          PROPSI(IB,JB)=SVAR
        ENDDO
      ENDDO
!
!     =====================================================================
!     == COLLECT TRANSFORMATION MATRIX BETWEEN NEW AND OLD               ==
!     =====================================================================
      TRANSPRO(:,:)=0.D0
      TRANSPHI(:,:)=0.D0
      DO IB=1,NAUG
        TRANSPRO(IB,IB)=1.D0
        TRANSPHI(IB,IB)=1.D0
      ENDDO
      TRANSPROPSI(:,:)=PROPSI(:,:)
!
      DO IB=1,NAUG
        DO JB=1,IB-1
!         == ORTHOGONALIZE PARTIAL WAVES ===================================
          SVAR=TRANSPROPSI(JB,IB)
          TRANSPHI(:,IB)   =TRANSPHI(:,IB)-TRANSPHI(:,JB)*SVAR
          TRANSPROPSI(:,IB)=TRANSPROPSI(:,IB)-TRANSPROPSI(:,JB)*SVAR          
!         == ORTHOGONALIZE PROJECTOR=====================================
          SVAR=TRANSPROPSI(IB,JB)
          TRANSPRO(:,IB)   =TRANSPRO(:,IB)-TRANSPRO(:,JB)*SVAR
          TRANSPROPSI(IB,:)=TRANSPROPSI(IB,:)-TRANSPROPSI(JB,:)*SVAR          
        ENDDO
        SVAR=TRANSPROPSI(IB,IB)
        TRANSPRO(:,IB)=TRANSPRO(:,IB)/SVAR
        TRANSPROPSI(IB,:)=TRANSPROPSI(IB,:)/SVAR
      ENDDO
!
!     =====================================================================
!     == CHECK RESULT                                                    ==
!     =====================================================================
      TRANSPROPSI(:,:)=MATMUL(TRANSPOSE(TRANSPRO),MATMUL(PROPSI,TRANSPHI))
      DO IB=1,NAUG
        TRANSPROPSI(IB,IB)=TRANSPROPSI(IB,IB)-1.D0
      ENDDO
      SVAR=MAXVAL(TRANSPROPSI)
      IF(SVAR.GT.1.D-5) THEN
do jb=1,naug
write(*,fmt='(20e10.2)')transpropsi(:,jb)
enddo
        CALL ERROR$MSG('BIORTHOGONALIZATION INACCURATE')
        CALL ERROR$R8VAL('MAX. DEV.',SVAR)
        CALL ERROR$STOP('BIORTHOMATRICES')
      END IF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE BAREPROJECTORS(GID,NR,L,E,PSPOT,RC &
     &                         ,NAUG,PSPHI,TPSPHI,PRO)
!     ** THE IDEA IS                                                      **
!     **    (T+VTILDE-E)|PSPHI(E)>=|P(E)>                                 **
!     **  WITH                                                            **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: L
      REAL(8)    ,INTENT(IN)     :: E(NAUG)
      REAL(8)    ,INTENT(IN)     :: PSPOT(NR)
      REAL(8)    ,INTENT(IN)     :: RC
      INTEGER(4) ,INTENT(IN)     :: NAUG
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: TPSPHI(NR,NAUG)
      REAL(8)    ,INTENT(OUT)    :: PRO(NR,NAUG)
      INTEGER(4)                 :: IB,IR
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: R(NR)
      LOGICAL(4)                 :: TTEST=.TRUE.
!     *********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      DO IB=1,NAUG
        PRO(:,IB)=TPSPHI(:,IB)+(PSPOT(:)*Y0-E(IB))*PSPHI(:,IB)
        PRO(NR,IB)=2.D0*PRO(NR-1,IB)-PRO(NR-2,IB)
      ENDDO
!     =====================================================================
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IF(R(IR).GT.RC) THEN
          PRO(IR:,:)=0.D0
          EXIT
        END IF
      ENDDO
!     =====================================================================
      IF(TTEST) THEN
        IF(L.EQ.0) THEN
          OPEN(100,FILE='BAREPROJECTORS_S.DAT')
        ELSE IF(L.EQ.1) THEN
          OPEN(100,FILE='BAREPROJECTORS_P.DAT')
        ELSE IF(L.EQ.2) THEN
          OPEN(100,FILE='BAREPROJECTORS_D.DAT')
        ELSE IF(L.EQ.3) THEN
          OPEN(100,FILE='BAREPROJECTORS_F.DAT')
        END IF
        DO IR=1,NR
          WRITE(100,FMT='(22F30.10)')R(IR),PRO(IR,:)
        ENDDO 
        CLOSE(10)
      END IF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE PSEUDIZE(GID,NR,POW,VAL0,RC,AEF,PSF)
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: POW   !POWER
      REAL(8)    ,INTENT(IN)     :: VAL0
      REAL(8)    ,INTENT(IN)     :: RC
      REAL(8)    ,INTENT(IN)     :: AEF(NR)
      REAL(8)    ,INTENT(OUT)    :: PSF(NR)
      REAL(8)                    :: VAL,DER
      REAL(8)                    :: V1,V2
      REAL(8)                    :: A11,A12,A21,A22
      REAL(8)                    :: DET
      REAL(8)                    :: AINV11,AINV12,AINV21,AINV22
      REAL(8)                    :: C0,C1,C2
      INTEGER(4)                 :: IR,IR0
      REAL(8)                    :: R(NR)
      LOGICAL(4),PARAMETER       :: TTEST=.FALSE.
      REAL(8)                    :: SVAR
      REAL(8)                    :: AUX(NR)
!     *********************************************************************
      CALL RADIAL$VALUE(GID,NR,AEF,RC,VAL)
      CALL RADIAL$DERIVATIVE(GID,NR,AEF,RC,DER)
!     == EQUATION SYSTEM A*C=V      
      V1=VAL-VAL0
      V2=DER
      A11=RC**POW
      A12=RC**(POW+2)
      A21=REAL(POW,KIND=8)*RC**(POW-1)
      A22=REAL(POW+2,KIND=8)*RC**(POW+1)
!     == INVERT MATRIX
      DET=A11*A22-A12*A21
      AINV11=A22/DET
      AINV22=A11/DET
      AINV12=-A12/DET
      AINV21=-A21/DET
!     == SOLVE EQUATION SYSTEM
      C0=VAL0
      C1=AINV11*V1+AINV12*V2
      C2=AINV21*V1+AINV22*V2
!     == CALCULATE PS FUNCTION
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IF(R(IR).GT.RC) THEN
          IR0=IR
          EXIT
        END IF
      ENDDO
      PSF(:IR0-1)=C0+R(:IR0-1)**POW*(C1+C2*R(:IR0-1)**2)
      PSF(IR0:)=AEF(IR0:)
      IF(TTEST) THEN
        SVAR=C0+C1*RC**POW+C2*RC**(POW+2)
        PRINT*,'TEST PESUDIZE 1 ',SVAR,VAL,SVAR-VAL
        SVAR=REAL(POW,8)*C1*RC**(POW-1)+REAL(POW+2,8)*C2*RC**(POW+1)
        PRINT*,'TEST PESUDIZE 2 ',SVAR,DER,SVAR-DER
        SVAR=C0
        PRINT*,'TEST PESUDIZE 3 ',SVAR,VAL0,SVAR-VAL0
        AUX(:)=C0+R(:)**POW*(C1+C2*R(:)**2)
!CALL WRITECOMPARE('TEST1.DAT',GID,NR,1,AEF,AUX)
!STOP
      END IF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE AESCF(GID,NR,NB,LOFI,SO,F,NN,Z,RHOADD,DREL,POT,EOFI)
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NB        ! #(STATES)
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)  !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO(NB)    !SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: F(NB)     !OCCUPATION
      INTEGER(4) ,INTENT(IN)     :: NN(NB)    !#(NODES)
      REAL(8)    ,INTENT(IN)     :: Z         !ATOMIC NUMBER
      REAL(8)    ,INTENT(IN)     :: RHOADD(NR)! FIXED CONTRIBUTION OF DENSITY
      REAL(8)    ,INTENT(OUT)    :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(INOUT)  :: POT(NR)   !POTENTIAL
      REAL(8)    ,INTENT(OUT)    :: EOFI(NB)  !ONE-PARTICLE ENERGIES
      REAL(8)                    :: R(NR)
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: RHO(NR)
      REAL(8)                    :: MASS
      REAL(8)                    :: U
      REAL(8)                    :: EREF
      INTEGER(4)                 :: ITER
      INTEGER(4)                 :: NITER=5000
      REAL(8)                    :: XAV,XMAX
      LOGICAL(4)                 :: CONVG
      REAL(8)   ,PARAMETER       :: TOL=1.D-5
      INTEGER(4)                  :: NFILO
      REAL(8)                    :: EH,EXC
      logical(4),parameter       :: tbroyden=.true.
      real(8)                    :: potin(nr)
!     ***********************************************************************
      call radial$r(gid,nr,r)
      eofi=0.d0
      EREF=0.D0
      XAV=0.D0
      XMAX=0.D0
      CONVG=.FALSE.
      IF(TBROYDEN) THEN
        CALL BROYDEN$NEW(NR,4,1.D0)
      ELSE
        CALL MYMIXPOT('ON',GID,NR,POT,XMAX,XAV)
      END IF
      DO ITER=1,NITER
        CALL RELATIVISTICCORRECTION(GID,NR,POT,EREF,DREL)
        CALL AERHO('FULL',GID,NR,NB,LOFI,SO,F,NN,DREL,POT,RHO,EOFI)
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
          if(tbroyden) then
            call broyden$clear
          else
            CALL MYMIXPOT('OFF',GID,NR,POT,XMAX,XAV)
          end if
          RETURN
        END IF

!       ====================================================================
!       == CALCULATE OUTPUT POTENTIAL                                     ==
!       ====================================================================
        if(tbroyden) potin=pot
        CALL MYVOFRHO(GID,NR,Z,RHO,POT,EH,EXC)
!
!       ================================================================
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSON'S METHOD     ==
!       ================================================================
        if(tbroyden) then
          xav=sqrt(dot_product(pot-potin,pot-potin)/real(nr,kind=8))
          xmax=maxval(abs(pot-potin)) 
          call broyden$step(nr,potin,(pot-potin))
          pot=potin
       else
          CALL MYMIXPOT('GO',GID,NR,POT,XMAX,XAV)
        end if
print*,'xmax ',xmax,xav
        CONVG=(XMAX.LT.TOL)
      ENDDO
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
      REAL(8)               :: XEXP
      REAL(8)               :: AUX(NR)
      REAL(8)               :: EDEN(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR),RI
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      INTEGER(4)            :: IR
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/DSQRT(FOURPI)
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  TOTAL POTENTIAL                                             ==
!     ==================================================================
      EDEN=0.D0
      CALL RADIAL$POISSON(GID,NR,0,RHO,POT)
      EDEN(:)=EDEN(:)+0.5D0*RHO(:)*POT(:)
      CALL NUCPOT(AEZ,NR,R,AUX)
      POT(:)=POT(:)+AUX(:)/Y0
      EDEN(:)=EDEN(:)+RHO(:)*AUX(:)/Y0
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
        XAV=DSQRT(XAV/SVAR1)
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
!     ......................................................................
      SUBROUTINE NODELESS(ID,GID,NR,NB,LOFI,SO,NN,DREL,POT,PHI,TPHI,E)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NB        ! #(STATES)
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)  !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO(NB)    !SWITCH FOR SPIN-ORBIT COUP.
      INTEGER(4) ,INTENT(IN)     :: NN(NB)    !#(NODES)
      REAL(8)    ,INTENT(IN)     :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: POT(NR)   !POTENTIAL
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,NB)   !WAVE FUNCTIONS
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR,NB)  !KINETIC ENERGY*WAVE FUNCTIONS
      REAL(8)    ,INTENT(inOUT)    :: E(NB)     !ONE-PARTICLE ENERGIES
      REAL(8)                    :: R(NR)
      REAL(8)                    :: G(NR)
      REAL(8)                    :: DREL1(NR)
      REAL(8)                    :: AUX(NR)
      INTEGER(4)                 :: IB,IR,JB,IBLAST
      REAL(8)                    :: SVAR
      REAL(8)                    :: C0LL,PI,Y0
      LOGICAL(4)                 :: TCHK
!     ***********************************************************************
      IF(ID.NE.'FULL'.AND.ID.NE.'EFFZORA'.AND.ID.NE.'NONREL') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$MSG('ALLOWED VALUES ARE "FULL",EFFZORA","NONREL"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('nodeless')
      ENDIF
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      DO IB=1,NB
        IF(NN(IB).EQ.0) THEN
          G(:)=0.D0
        ELSE
          TCHK=.FALSE.
          DO JB=1,IB-1
            IF(LOFI(JB).NE.LOFI(IB)) CYCLE
            IF(SO(JB).NE.SO(IB)) CYCLE
            IF(NN(JB).NE.NN(IB)-1) CYCLE
            G(:)=PHI(:,JB)
            TCHK=.TRUE.
            EXIT
          ENDDO
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('PREVIOUS STATE NOT FOUND')
            CALL ERROR$STOP('NODELESS')
          END IF
        END IF
        IF(ID.EQ.'FULL') THEN
          CALL RELATIVISTICCORRECTION(GID,NR,POT,E(IB),DREL1)
        ELSE IF(ID.EQ.'NONREL') THEN
          DREL1(:)=0.D0 
        ELSE IF(ID.EQ.'EFFZORA') THEN
          DREL1(:)=DREL(:)
        END IF
        CALL BOUNDSTATE(GID,NR,LOFI(IB),SO(IB),DREL1,G,0,POT,E(IB),PHI(:,IB))
        TPHI(:,IB)=G(:)+(E(IB)-POT(:)*Y0)*PHI(:,IB)
!
!       =================================================================
!       == NORMALIZE LOWEST WAVE FUNCTION                              ==
!       =================================================================
        IF(NN(IB).EQ.0) THEN
          AUX(:)=(R(:)*PHI(:,IB))**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          SVAR=1.D0/SQRT(SVAR)
!         == CHANGE SIGN SO THAT NODELESS FUNCTIONS ARE POSITIVE
          DO IR=1,NR
            IF(PHI(IR,IB).EQ.0.D0) CYCLE
            IF(PHI(IR,IB).LT.0.D0) SVAR=-SVAR
            EXIT
          ENDDO
          PHI(:,IB) =PHI(:,IB) *SVAR
          TPHI(:,IB)=TPHI(:,IB)*SVAR
        END IF
      ENDDO
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE ONENODELESS(ID,GID,NR,L,NN,SO,DREL,POT,PHIC,E,PHI,TPHI)
!     **                                                                  **
!     **  CALCULATES THE NODELESS FUNCTIONS FOR A GIVEN NODE-LESS         **
!     **  CORE FUNCTION PHIC                                              **
!     **     [T+POT-E]|PHI>=|PHIC>                                        **
!     **  FOR ID='BOUND' THE BOUND STATE WITH NN NODES IS DETERMINED      **
!     **      AND THE NEW ENERGY IS RETURNED                              **
!     **  FOR ID='SCATT' THE SOLUTION FOR THE ENERGY GIVEN IS PROVIDED    **
!     **                                                                  **
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: L         !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: NN        !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO        !SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: POT(NR)   !POTENTIAL
      REAL(8)    ,INTENT(IN)     :: PHIC(NR)  !HIGHEST NODELESS CORE FUNCTION
      REAL(8)    ,INTENT(INOUT)  :: E         !EXPANSION ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE FUNCTIONS
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR)!KINETIC ENERGY*WAVE FUNCTIONS
      INTEGER(4)                 :: IB
      REAL(8)                    :: PI,Y0
!     ***********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      IF(ID.EQ.'BOUND') THEN
        CALL BOUNDSTATE(GID,NR,L,SO,DREL,PHIC,NN,POT,E,PHI)
      ELSE IF(ID.EQ.'SCATT') THEN
        CALL RADIAL$SCHRODINGER(GID,NR,POT,DREL,SO,PHIC,L,E,1,PHI)
      END IF
      TPHI(:)=PHIC(:)+(E-POT(:)*Y0)*PHI(:)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE AERHO(ID,GID,NR,NB,LOFI,SO,F,NN,DREL,POT,RHO,E)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NB        ! #(STATES)
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)  !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO(NB)    !SWITCH FOR SPIN-ORBIT COUP.
      INTEGER(4) ,INTENT(IN)     :: NN(NB)    !#(NODES)
      REAL(8)    ,INTENT(IN)     :: F(NB)     !OCCUPATION
      REAL(8)    ,INTENT(INOUT)  :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: POT(NR)   !POTENTIAL
      REAL(8)    ,INTENT(OUT)    :: RHO(NR)   !DENSITY
      REAL(8)    ,INTENT(INOUT)  :: E(NB)     !ONE-PARTICLE ENERGIES
      REAL(8)                    :: R(NR)
      REAL(8)                    :: G(NR)
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: PHI(NR)
      REAL(8)                    :: EKIN(NR)
      REAL(8)                    :: EKIN2(NR)
      REAL(8)                    :: ARRAY(NR,5)
      INTEGER(4)                 :: IB,IR
      REAL(8)                    :: SVAR
      REAL(8)                    :: C0LL,PI,Y0
!     ***********************************************************************
      IF(ID.NE.'FULL'.AND.ID.NE.'EFFZORA'.AND.ID.NE.'NONREL') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$MSG('ALLOWED VALUES ARE "FULL",EFFZORA","NONREL"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AESCF')
      ENDIF
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      RHO(:)=0.D0
      EKIN(:)=0.D0
      EKIN2(:)=0.D0
      DO IB=1,NB
         G(:)=0.D0
         IF(ID.EQ.'FULL') THEN
           CALL RELATIVISTICCORRECTION(GID,NR,POT,E(IB),DREL)
         ELSE IF(ID.EQ.'NONREL') THEN
           DREL(:)=0.D0 
         END IF
         CALL BOUNDSTATE(GID,NR,LOFI(IB),SO(IB),DREL,G,NN(IB),POT,E(IB),PHI)
         AUX(:)=(R(:)*PHI(:))**2
         CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
         PHI(:)=PHI(:)/SQRT(SVAR)
         RHO(:)  =RHO(:)  +F(IB)*C0LL         *PHI(:)**2
         EKIN(:) =EKIN(:) +F(IB)*C0LL*E(IB)   *PHI(:)**2
         EKIN2(:)=EKIN2(:)+F(IB)*C0LL*E(IB)**2*PHI(:)**2
      ENDDO
      IF(ID.EQ.'FULL') THEN
        EKIN(:)=EKIN(:)/RHO(:)
        EKIN2(:)=EKIN2(:)/RHO(:)
        EKIN2(:)=SQRT(EKIN2(:)-EKIN(:)**2)
        EKIN(:)=EKIN(:)/Y0-POT(:)
        CALL RELATIVISTICCORRECTION(GID,NR,-EKIN,0.D0,DREL)
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
!     ......................................................................
      SUBROUTINE BOUNDSTATE(GID,NR,L,SO,DREL,G,NN,POT,E,PHI)
!     **                                                                  **
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND     **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                 **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G            **
!     **                                                                  **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE             **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.     **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.     **
!     **                                                                  **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND       **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT       **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                     **
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,Z0,DX,XM,ZM
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: VAL,DER,DERO,RTEST
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      INTEGER(4)                 :: NN1
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+100 ! MAXIMUM FACTOR IN THE WAVE FUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      integer(4)                 :: nn1m,nn10,nn1p
      real(8)                    :: phip(nr),phim(nr)
      real(8)                    :: rcl
      real(8)                    :: swkb(nr)  ! phi=e^S
      real(8)                    :: aux(nr),aux1(nr)
      integer(4)                 :: irout,ircl
!     *********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      RTEST=R(NR)
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
        POT1(:)=POT(:)
!       == FIND CLASSICAL TURNING POINT
        RCL=R(NR)
        IF(E.LT.POT(NR)*Y0) THEN
          DO IR=NR-1,1,-1
            IF(E.GT.POT(IR)*Y0) THEN
              IRCL=IR
              RCL=R(IR)-(POT(IR)-E/Y0)/(POT(IR+1)-POT(IR))*(R(IR+1)-R(IR))
! RCL=MIN(RTEST,R(NR))
              EXIT
            END IF
          ENDDO
          RTEST=RCL
        END IF
!
!       == USE WKB SOLUTION FOR THE SCHR.GL. FOR A CONSTANT POTENTIAL AND L=0
!       == TO ESTIMATE FACTOR FROM RTEST TO OUTERMOST POINT
        irout=nr
        IF(RTEST.LT.R(NR)) THEN
          AUX(:IRCL)=0.D0
          AUX(IRCL+1:)=SQRT(REAL(L*(L+1),kind=8)/R(ircl+1:)**2+2.D0*(POT(IRCL+1:)*Y0-E))
          CALL RADIAL$DERIVe(GID,NR,AUX,AUX1)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,Swkb)
!          Swkb(:)=Swkb(:)-0.5D0*LOG(AUX1(:))-LOG(R(:))
          Swkb(:)=Swkb(:)-LOG(R(:))
!         == DETERMINE IROUT WHERE THE WAVE FUNCTION CAN GROW BY A FACTOR 
!         == OF XMAX FROM THE CLASSICAL TURNING POINT
          SVAR=LOG(XMAX)
          DO IR=1,NR
            IF(Swkb(IR).GT.SVAR) THEN
              IROUT=IR-1
              EXIT
            END IF
          ENDDO
        END IF
        svar=pot(irout)
        pot1(:)=pot(:)
        pot1(irout:)=pot(irout)
!print*,'r(irout)',rcl,r(irout),pot1(nr)
!
!       =======================================================================
!       == INTEGRATE RADIAL SCHRODINGER EQUATION OUTWARD                     ==
!       =======================================================================
        IDIR=1
        CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(irout).GT.0.OR.PHI(irout).Le.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NN',NN)
          CALL ERROR$STOP('BOUNDSTATE')
        END IF
!
!       =======================================================================
!       == CALCULATE PHASE SHIFT =========================================
!       =======================================================================
        NN1=0
        DO IR=3,irout
          IF(PHI(IR)*PHI(IR-1).LT.0.D0) NN1=NN1+1
        ENDDO
        Z0=-1.D0
        IF(NN1.Gt.NN) Z0=1.D0
!write(*,fmt='("iter",4i5,4e20.10)')I,L,NN1,NN,E,2.d0*DX,phi(irout),z0
        IF(ABS(2.d0*DX).LE.TOL) EXIT
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
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E+2.d0*dx,IDIR,PHIp)
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E-2.d0*dx,IDIR,PHIm)
      NN1p=-nn
      NN1m=-nn
      NN10=-nn
      DO IR=3,irout
        IF(PHI(IR)*PHI(IR-1).LT.0.D0)   NN10=NN10+1
        IF(PHIm(IR)*PHIm(IR-1).LT.0.D0) NN1m=NN1m+1
        IF(PHIp(IR)*PHIp(IR-1).LT.0.D0) NN1p=NN1p+1
      ENDDO
!Print*,'nn1-nn ',nn1m,nn10,nn1p
!Print*,'e ',e-2.d0*dx,e,e+2.d0*dx
!
!     =======================================================================
!     ==  DETERMINE MATCHING POINT                                         ==
!     =======================================================================
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
      DO IR=1,NR
         IRMATCH=IR
         IF(POT(IR)-E/y0.GT.0.D0.OR.R(IR).GT.5.D0) EXIT
      ENDDO
!
!     =======================================================================
!     ==  INTEGRATE INWARD                                                 ==
!     =======================================================================
      IDIR=-1
      GHOM(:)=0.D0
      if(irout.lt.nr) GHOM(IROUT)=1.D-5
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
      IF(.NOT.THOM) THEN     
        CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
      ELSE
        PHIINHOM(:)=0.D0
      END IF
!
!     =======================================================================
!     ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!     =======================================================================
      SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
      PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
      CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
      CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
      SVAR=(DERO-DER)/PHI(IRMATCH)
      PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
do ir=1,nr
  if(.not.(phi(ir).gt.0.d0.or.phi(ir).le.0.d0)) then
    print*,'error'
    print*,'phiin',phi(:irmatch-1)
    print*,'phiout',phi(irmatch:)
    print*,'svar ',svar
    call error$stop('boundstate')
  end if
enddo
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE old1BOUNDSTATE(GID,NR,L,SO,DREL,G,NN,POT,E,PHI)
!     **                                                                  **
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND     **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                 **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G            **
!     **                                                                  **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE             **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.     **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.     **
!     **                                                                  **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND       **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT       **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                     **
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,Z0,DX,XM,ZM
      REAL(8)    ,PARAMETER      :: TOL=1.D-8
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: VAL,DER,DERO,RTEST
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      INTEGER(4)                 :: NN1
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+100 ! MAXIMUM FACTOR IN THE WAVE FUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      integer(4)                 :: nn1m,nn10,nn1p
      real(8)                    :: phip(nr),phim(nr)
      real(8)                    :: rcl,ircl
!     *********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      RTEST=R(NR)
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
        pot1(:)=pot(:)
!       == FIND CLASSICAL TURNING POINT
        Rcl=R(NR)
        IF(E.lt.POT(NR)*y0) THEN
          DO IR=NR-1,1,-1
            IF(E.gt.POT(IR)*y0) THEN
              rcl=R(IR)-(POT(IR)-E/y0)/(POT(IR+1)-POT(IR))*(R(IR+1)-R(IR))
! rcl=MIN(RTEST,R(NR))
              EXIT
            END IF
          ENDDO
          rtest=rcl
        END IF
!
!       == USE WKB SOLUTION FOR THE SCHR.GL. FOR A CONSTANT POTENTIAL AND L=0
!       == TO ESTIMATE FACTOR FROM RTEST TO OUTERMOST POINT
        IF(RTEST.LT.R(NR)) THEN
          SVAR=0.5D0*(LOG(XMAX*R(NR)/RTEST)/(R(NR)-RTEST))**2
          do ir=nr-1,1,-1
            if(pot(ir)*y0.gt.E+SVAR) cycle
write(*,fmt='("pot changed",8f20.5)')svar,e,pot(ir)*y0,pot(ir+1)*y0,rtest,r(ir)
            exit
          enddo
          POT1(:)=MIN(POT(:),(E+SVAR)/y0)
        ELSE
          POT1(:)=POT(:)
        END IF
!
!       =======================================================================
!       == INTEGRATE RADIAL SCHRODINGER EQUATION OUTWARD                     ==
!       =======================================================================
        IDIR=1
        CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(NR).GT.0.OR.PHI(NR).LT.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NN',NN)
          CALL ERROR$STOP('BOUNDSTATE')
        END IF
!
!       =======================================================================
!       == CALCULATE PHASE SHIFT =========================================
!       =======================================================================
        NN1=0
        DO IR=3,NR
          IF(PHI(IR)*PHI(IR-1).LT.0.D0) NN1=NN1+1
        ENDDO
        Z0=-1.D0
        IF(NN1.Gt.NN) Z0=1.D0
write(*,fmt='("iter",4i5,4e20.10)')I,L,NN1,NN,E,2.d0*DX,phi(nr),z0
        IF(ABS(2.d0*DX).LE.TOL) EXIT
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
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E+2.d0*dx,IDIR,PHIp)
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E-2.d0*dx,IDIR,PHIm)
      NN1p=-nn
      NN1m=-nn
      NN10=-nn
      DO IR=3,NR-1
        IF(PHI(IR)*PHI(IR-1).LT.0.D0)   NN10=NN10+1
        IF(PHIm(IR)*PHIm(IR-1).LT.0.D0) NN1m=NN1m+1
        IF(PHIp(IR)*PHIp(IR-1).LT.0.D0) NN1p=NN1p+1
      ENDDO
!Print*,'nn1-nn ',nn1m,nn10,nn1p
!Print*,'e ',e-2.d0*dx,e,e+2.d0*dx
!
!     =======================================================================
!     ==  DETERMINE MATCHING POINT                                         ==
!     =======================================================================
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
      DO IR=1,NR
         IRMATCH=IR
         IF(POT(IR)-E/y0.GT.0.D0.OR.R(IR).GT.5.D0) EXIT
      ENDDO
!
!     =======================================================================
!     ==  INTEGRATE INWARD                                                 ==
!     =======================================================================
      IDIR=-1
      GHOM(:)=0.D0
      CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
      IF(.NOT.THOM) THEN     
        CALL RADIAL$SCHRODINGER(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
      ELSE
        PHIINHOM(:)=0.D0
      END IF
!
!     =======================================================================
!     ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!     =======================================================================
      SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
      PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
      CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
      CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
      SVAR=(DERO-DER)/PHI(IRMATCH)
!PRINT*,'SVAR',SVAR,R(IRMATCH),DER/PHI(IRMATCH),DERO/PHI(IRMATCH)
!do ir=1,nr
!  print*,r(ir),phi(ir),phiinhom(ir)
!enddo
!!$      IF(ABS(SVAR).GT.1.D-5) THEN
!!$        CALL ERROR$MSG('DERIVATIVES DO NOT MATCH')
!!$        CALL ERROR$R8VAL('STEP IN LOG. DERIVATIVE',SVAR)
!!$        CALL ERROR$R8VAL('RC',R(IRMATCH))
!!$        CALL ERROR$STOP('BOUNDSTATE')
!!$      END IF
!
      PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE NUCPOT(Z,NR,R,POT)
!     **                                                                  **
!     **  ELECTROSTATIC POTENTIAL OF A NUCLEUS WITH FINITE RADIUS         **
!     **  THE NUCLEUS IS CONSIDERED AS A HOMOGENEOUSLY CHARGED SPHERE.    **
!     **  THE RADIUS IS RELATED  TO THE TOTAL MASS (NUMBER OF NUCLEONS),  **
!     **  WHICH CAN BE LOOKED UP KNOWING THE ATOMIC NUMBER Z.             **
!     **                                                                  **
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: Z        ! ATOMIC NUMBER
      INTEGER(4),INTENT(IN) :: NR       ! #(RADIAL GRID POINTS)
      REAL(8)   ,INTENT(IN) :: R(NR)    ! RADIAL GRID
      REAL(8)   ,INTENT(OUT):: POT(NR)  ! POTENTIAL
      REAL(8)               :: RNUC     ! NUCLEAR RADIUS
      INTEGER(4)            :: IR
!     ***********************************************************************
      CALL PERIODICTABLE$GET(NINT(Z),'RNUC',RNUC)
      DO IR=1,NR
        IF(R(IR).GT.RNUC) THEN
           POT(IR)=-Z/R(IR)
        ELSE
          POT(IR)=-Z/RNUC*(1.5D0-0.5D0*(R(IR)/RNUC)**2)
        END IF
      ENDDO
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE RELATIVISTICCORRECTION(GID,NR,POT,E,DREL)
!     **                                                                  **
!     **  DREL IS A MEASURE OF THE RELATIVISTIC CORRECTIONS               **
!     **     D:=1/MREL-1/M0                                               **
!     **  WHERE MREL IS THE RELATIVISTIC MASS  MREL=M0+(E-POT)/(2C**2)    **
!     **  AND  M0 IS THE REST MASS                                        **
!     **                                                                  **
!     **  REMARKS:                                                        **
!     **  -  RELATIVISTIC CORRECTIONS FOR EKIN<0 ARE SWITCHED OFF!        **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID         ! GRID ID
      INTEGER(4),INTENT(IN) :: NR          ! #(RADIAL GRID POINTS
      REAL(8)   ,INTENT(IN) :: POT(NR)     ! POTENTIAL (MULTIPLY WITH Y0!)
      REAL(8)   ,INTENT(IN) :: E           ! ONE-PARTICLE ENERGY
      REAL(8)   ,INTENT(OUT):: DREL(NR)    ! RELATIVISTIC CORRECTION 
      INTEGER(4)            :: IR
      REAL(8)               :: C           ! SPEED OF LIGHT
      REAL(8)               :: PI,Y0      
      REAL(8)               :: DRELDOT(NR)
LOGICAL(4),SAVE :: TFIRST=.TRUE.
!     ***********************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL CONSTANTS$GET('C',C)
      DREL(:)=-1.D0/(1.D0+2.D0*C**2/MAX(E-POT(:)*Y0,0.D0))
      DRELDOT(:)=DREL(:)**2*2.D0*C**2/MAX(E-POT(:)*Y0,0.D0)**2
      DO IR=1,NR
        IF(E.LT.POT(IR)*Y0) THEN 
          DRELDOT(IR)=0.D0
        END IF
      ENDDO


!!$IF(TFIRST)PRINT*,'WARNING! RELATIVISTIC CORRECTION SWITCHED OFF',C
!!$TFIRST=.FALSE.
!!$DREL(:)=0.D0
      RETURN
      END

!
!     ..................................................................
      SUBROUTINE PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO,DH,DO,G,PHI)
!     **                                                              **
!     **  SOLVES THE RADIAL PAW -SCHROEDINGER EQUATION.               **
!     **    (T+VTILDE-E+|P>(DH-E*DO<P|]|PHI>=|G>                      **
!     **  WHERE T IS THE NONRELATIVISTIC KINETIC ENERGY.              **
!     **                                                              **
!     **    DH=<AEPHI|T+V|AEPHI>-<PSPHI|T+VTILDE|PSPHI>               **
!     **    DO=<AEPHI|AEPHI>-<PSPHI|PSPHI>                            **
!     **                                                              **
!     **                                                              **
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID     !GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR      ! #(RADIAL GRID POINTS
      INTEGER(4)  ,INTENT(IN) :: L       ! ANGULAR MOMENTUM
      REAL(8)     ,INTENT(IN) :: E       ! ENERGY
      REAL(8)     ,INTENT(IN) :: PSPOT(NR)  !VTILDE=PSPOT*Y0
      INTEGER(4)  ,INTENT(IN) :: NPRO       ! #(PROJECTOR FUNCTIONS)
      REAL(8)     ,INTENT(IN) :: PRO(NR,NPRO) ! PROJECTOR FUNCTIONS
      REAL(8)     ,INTENT(IN) :: DH(NPRO,NPRO)     
      REAL(8)     ,INTENT(IN) :: DO(NPRO,NPRO)     
      REAL(8)     ,INTENT(IN) :: G(NR)        !INHOMOGENEITY
      REAL(8)     ,INTENT(OUT):: PHI(NR)      !PAW RADIAL FUNCTION
      REAL(8)                 :: U(NR)
      REAL(8)                 :: V(NR,NPRO)
      REAL(8)                 :: AMAT(NPRO,NPRO)
      REAL(8)                 :: BMAT(NPRO,NPRO)
      REAL(8)                 :: BMATINV(NPRO,NPRO)
      REAL(8)                 :: CMAT(NPRO,NPRO)
      REAL(8)                 :: CVEC(NPRO)
      REAL(8)                 :: DVEC(NPRO)
      INTEGER(4)              :: IB,JB
      INTEGER(4)              :: I1,I2,I3
      REAL(8)                 :: SVAR
      REAL(8)                 :: AUX(NR)
      REAL(8)                 :: R(NR)
      REAL(8)                 :: DREL(NR)
      INTEGER(4)              :: SO
INTEGER(4)              :: IR
REAL(8)                 :: PI,Y0
!     ******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  -1/2NABLA^2+POT-E|U>=|G>                                    ==
!     ==================================================================
      SO=0
      DREL(:)=0.D0
      CALL RADIAL$SCHRODINGER(GID,NR,PSPOT,DREL,SO,G,L,E,1,U)
!
!     ==================================================================
!     ==  -1/2NABLA^2+POT-E|V>=|PRO>                                 ==
!     ==================================================================
      DO I1=1,NPRO
        CALL RADIAL$SCHRODINGER(GID,NR,PSPOT,DREL,SO,PRO(:,I1),L,E,1,V(:,I1))
      ENDDO
!!$PRINT*,'E ',E,SO,L,NR,GID
!!$OPEN(100,FILE='XXX.DAT')
!!$ DO IR=1,NR
!!$    WRITE(100,FMT='(22F30.10)')R(IR),U(IR),V(IR,:),PRO(IR,:)
!!$ ENDDO
!!$CLOSE(10)
!!$STOP 'IN PAWDER'
!
!     ==================================================================
!     ==  AMAT=<PRO|V>  CVEC=<PRO|U>                                  ==
!     ==================================================================
      DO I1=1,NPRO
        DO I2=1,NPRO
          AUX(:)=PRO(:,I1)*V(:,I2)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,AMAT(I1,I2))
        ENDDO
      ENDDO
      DO I1=1,NPRO
        AUX(:)=PRO(:,I1)*U(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,CVEC(I1))
      ENDDO  
!
!     ==================================================================
!     ==  BMAT=1+(DATH-EDO)<PRO|V>                                    ==
!     ==================================================================
!     BMAT(:,:)=MATMUL(DH(:,:)-E*DO(:,:),AMAT(:,:))
!     DO I1=1,NPRO
!       BMAT(I1,I1)=BMAT(I1,I1)+1.D0
!     ENDDO
!
      DO I1=1,NPRO
        DO I2=1,NPRO
          BMAT(I1,I2)=0.D0
          DO I3=1,NPRO
            BMAT(I1,I2)=BMAT(I1,I2)+(DH(I1,I3)-E*DO(I1,I3))*AMAT(I3,I2)     
          ENDDO
        ENDDO
        BMAT(I1,I1)=BMAT(I1,I1)+1.D0
      ENDDO
!
!     ==================================================================
!     ==  BMAT = BMAT^-1 = [1+(DATH-EDO)<PRO|V>]^-1                   ==
!     ==================================================================
      IF(NPRO.EQ.0) THEN
        CALL ERROR$STOP('NPRO=0 NOT ALLOWED')
      END IF
      IF(NPRO.EQ.1) THEN 
        BMAT(1,1)=1.D0/BMAT(1,1)
      ELSE 
        CALL LIB$INVERTR8(NPRO,BMAT,BMATINV)
        BMAT=BMATINV
      END IF
!
!     ==================================================================
!     ==  CMAT = -BMAT*(DATH-EDO)                                     ==
!     ==       = -[1+(DATH-EDO)*<PRO|V>]^-1 (DATH-EDO)                ==
!     ==================================================================
!     CMAT(:,:)=MATMUL(BMAT(:,:),DH(:,:)-E*DO(:,:))
      DO I1=1,NPRO
        DO I2=1,NPRO
          CMAT(I1,I2)=0.D0
          DO I3=1,NPRO
            CMAT(I1,I2)=CMAT(I1,I2)-BMAT(I1,I3)*(DH(I3,I2)-E*DO(I3,I2))
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DVEC = CMAT*CVEC                                           ==
!     ==       = -[1+(DATH-EDO)*<PRO|V>]^-1 (DATH-EDO) <PRO|U>        ==
!     ==================================================================
!     DVEC(:)=MATMUL(CMAT(:,:),CVEC(:))
      DO I1=1,NPRO
        DVEC(I1)=0.D0
        DO I2=1,NPRO
          DVEC(I1)=DVEC(I1)+CMAT(I1,I2)*CVEC(I2)
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  |PHI> = |U>+|V>DVEC                                         ==
!     ==  = [1-|V>[1+(DATH-EDO)*<PRO|V>]^-1(DATH-EDO)<PRO|] |U>       ==
!     ==================================================================
!     PHI(:)=U(:)+MATMUL(V(:,:),DVEC(:))
      PHI(:)=U(:)
      DO I1=1,NPRO
        PHI(:)=PHI(:)+V(:,I1)*DVEC(I1)
      ENDDO
      RETURN
      STOP
      END
!
!     ...................................................................
      SUBROUTINE GRIDPARAMETERS(DMIN,DMAX,RX,R1,DEX,NR)
!     **                                                               **
!     **  DETERMINES THE GRID PARAMETERS FOR THE SHIFTED LOGARITHMIC   **
!     **  GRID FROM A SPECIFIED MINIMUM SPACING DMIN, A MAXIMUM        **
!     **  SPACING DMAX AND A MAXIMUM RADIUS RX                         **
!     **                                                               **
      REAL(8),   INTENT(IN) :: DMIN
      REAL(8),   INTENT(IN) :: DMAX
      REAL(8),   INTENT(IN) :: RX
      REAL(8),   INTENT(OUT):: R1
      REAL(8),   INTENT(OUT):: DEX
      INTEGER(4),INTENT(OUT):: NR
      REAL(8)               :: RN
      REAL(8)               :: Q   ! EXP(DEX)
!     *******************************************************************
      RN=2.D0+LOG(DMAX/DMIN)/LOG((RX-DMIN)/(RX-DMAX))
      Q=(DMAX/DMIN)**(1.D0/(RN-2.D0))
      DEX=LOG(Q)
      R1=DMIN/(Q-1)
!
!      PRINT*,'DMIN ',DMIN,R1*(EXP(DEX)-1.D0)            
!      PRINT*,'DMAX ',DMAX,R1*(EXP(DEX*REAL(RN-1))-EXP(DEX*REAL(RN-2)))
!      PRINT*,'RX   ',RX,R1*(EXP(DEX*REAL(RN-1))-1.D0)
      NR=NINT(RN)
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE AEETOT(GID,NR,NB,Z,NC,EOFI,FOFI,AERHOC,AERHOV,AEPOT)
!     **                                                               **
!     **                                                               **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NC
      REAL(8),   INTENT(IN) :: Z
      REAL(8),   INTENT(IN) :: AERHOC(NR)
      REAL(8),   INTENT(IN) :: AERHOV(NR)
      REAL(8)   ,INTENT(IN) :: EOFI(NB)
      REAL(8)   ,INTENT(IN) :: FOFI(NB)
      REAL(8)   ,INTENT(IN) :: AEPOT(NR)
      REAL(8)               :: EKINC,EKINV
      REAL(8)               :: EHC,EHV
      REAL(8)               :: EXCC,EXCV
      REAL(8)               :: R(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR
      INTEGER(4)            :: NFILO
!     *******************************************************************
      CALL RADIAL$R(GID,NR,R)
!     == CORE KINETIC ENERGY ============================================
      AUX(:)=AERHOC(:)*AEPOT(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      EKINC=SUM(EOFI(1:NC)*FOFI(1:NC))-SVAR
!     == VALENCE KINETCI ENERGY ==========================================
      AUX=AERHOV(:)*AEPOT(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      EKINV=SUM(EOFI(NC+1:)*FOFI(NC+1:))-SVAR
!     == CORE POTENTIAL ENERGY  ==========================================
      CALL MYVOFRHO(GID,NR,Z,AERHOC,AUX,EHC,EXCC)
!     == CORE POTENTIAL ENERGY  ==========================================
      CALL MYVOFRHO(GID,NR,Z,AERHOC+AERHOV,AUX,EHV,EXCV)
      EHV=EHV-EHC
      EXCV=EXCV-EXCC
!     == WRITE ===========================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL REPORT$TITLE(NFILO,'AE-ENERGIES')
      CALL REPORT$R8VAL(NFILO,'ETOT(VALENCE ONLY)',EKINV+EHV+EXCV,'H')
      CALL REPORT$R8VAL(NFILO,'ETOT',EKINV+EKINC+EHV+EHC+EXCV+EXCC,'H')
      CALL REPORT$R8VAL(NFILO,'E-BAND',SUM(EOFI(NC+1:)*FOFI(NC+1:)),'H')
      CALL REPORT$R8VAL(NFILO,'E-KIN',EKINV,'H')
      CALL REPORT$R8VAL(NFILO,'E-HARTREE',EHV,'H')
      CALL REPORT$R8VAL(NFILO,'E-XC',EXCV,'H')
      CALL REPORT$R8VAL(NFILO,'E-XC-CORE',EXCC,'H')
      CALL REPORT$R8VAL(NFILO,'E-XC-VALENCE',EXCV,'H')
      CALL REPORT$R8VAL(NFILO,'E-CORE',EKINC+EHC+EXCC,'H')
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE PSETOTSHELL(GID,NR,AEZ,PSRHOC,AERHOC,PSPOT &
     &                      ,NB,NC,LMAX,NAUG,NPRO,LOFI,FOFI,EOFI &
     &                      ,PRO1,AEPHI1,PSPHI1,RCSM,DTKIN1,DH1,DO1,VADD)
!     **                                                               **
!     **                                                               **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: AEZ
      REAL(8),   INTENT(IN) :: PSRHOC(NR)
      REAL(8),   INTENT(IN) :: AERHOC(NR)
      REAL(8),   INTENT(IN) :: PSPOT(NR)
      INTEGER(4),INTENT(IN) :: NC
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LMAX
      INTEGER(4),INTENT(IN) :: NAUG
      INTEGER(4),INTENT(IN) :: NPRO(LMAX+1)
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      REAL(8),   INTENT(IN) :: FOFI(NB)
      REAL(8),   INTENT(IN) :: EOFI(NB)
      REAL(8)   ,INTENT(IN) :: PRO1(NR,NAUG,LMAX+1)
      REAL(8)   ,INTENT(IN) :: AEPHI1(NR,NAUG,LMAX+1)
      REAL(8)   ,INTENT(IN) :: PSPHI1(NR,NAUG,LMAX+1)
      REAL(8)   ,INTENT(IN) :: RCSM
      REAL(8)   ,INTENT(IN) :: DTKIN1(NAUG,NAUG,LMAX+1)
      REAL(8)   ,INTENT(IN) :: DH1(NAUG,NAUG,LMAX+1)
      REAL(8)   ,INTENT(IN) :: DO1(NAUG,NAUG,LMAX+1)
      REAL(8)   ,INTENT(IN) :: VADD(NR)
      REAL(8)               :: PI,Y0
      REAL(8)               :: R(NR)
      INTEGER(4)            :: ISTART,IBI
      INTEGER(4),PARAMETER  :: NITER=100
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: X0,F0,XM,YM,DX
      REAL(8)               :: DF0
      REAL(8)               :: E
      REAL(8)               :: G(NR)
      REAL(8)   ,ALLOCATABLE:: PSPSI(:,:)
      REAL(8)   ,ALLOCATABLE:: TPSPSI(:,:)
      INTEGER(4)            :: LNX
      INTEGER(4)            :: L1,L2,LN1,LN2,IPRO,IPRO1,IPRO2,IB,I,ir
      INTEGER(4)            :: L,NN
      INTEGER(4)            :: IRMATCH
      REAL(8)               :: DREL(NR)
      REAL(8)               :: AUX(NR)
      INTEGER(4),ALLOCATABLE:: LOX(:)
      REAL(8)               :: SVAR
      REAL(8)               :: PROJ(NAUG)
      REAL(8)   ,ALLOCATABLE:: DTKIN(:,:)
      REAL(8)   ,ALLOCATABLE:: PRO(:,:)
      REAL(8)   ,ALLOCATABLE:: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE:: PSPHI(:,:)
      integer(4)            :: nfilo
!     *******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL REPORT$TITLE(NFILO,'EIGENSTATES FROM PAW CALCULATION')
!
!     ===================================================================================
!     == find paw bound state by bisection                                             ==
!     ===================================================================================
      ALLOCATE(PSPSI(NR,NB-NC))
      ALLOCATE(TPSPSI(NR,NB-NC))
      DO IB=NC+1,NB
        L=LOFI(IB)
        DF0=1.D0
        DO I=NC+1,IB-1
          IF(LOFI(I).EQ.L) DF0=DF0+1.D0
        ENDDO
!
!       =================================================================================
!       == DETERMINE IRMATCH, WHERE THE OUTWARD AND INWARD SOLUTIONS ARE MATCHED ========
!       =================================================================================
        irmatch=10.d0
        DO IR=1,NR
          DO I=1,NPRO(L+1)
            IF(ABS(PRO1(IR,I,L+1)).gt.1.D-12) IRMATCH=IR
          ENDDO
        ENDDO
        IF(R(IRMATCH).GT.5.D0) THEN
          CALL ERROR$MSG('MATCHING RADII NOT FOUND')
          CALL ERROR$STOP('PSETOTSHELL')
        END IF
!
!       =================================================================================
!       == find paw bound state by bisection                                           ==
!       =================================================================================
        ISTART=1
        X0=EOFI(IB)
        DX=0.1D0
        CALL BISEC(ISTART,IBI,X0,F0,DX,XM,YM)
        DO I=1,NITER
          E=X0
          G(:)=0.D0
          CALL PAWDER(GID,NR,L,E,PSPOT,NPRO(L+1),PRO1(:,1:NPRO(L+1),L+1) &
                    ,DH1(1:NPRO(L+1),1:NPRO(L+1),L+1),DO1(1:NPRO(L+1),1:NPRO(L+1),L+1) &
       &            ,G,PSPSI(:,IB-NC))
          CALL PHASESHIFT(GID,NR,PSPSI(:,IB-NC),R(nr),F0)
          F0=F0-DF0
          CALL BISEC(ISTART,IBI,X0,F0,DX,XM,YM)
          IF(ABS(DX).LE.TOL) EXIT
        ENDDO
        IF(ABS(DX).GT.TOL) THEN
          CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
          CALL ERROR$MSG('BOUND STATE NOT FOUND')
          CALL ERROR$STOP('BOUNDSTATE')
        END IF
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2," E[H]=",F15.9)') &
     &              L,0,0,FOFI(IB),E
!
!       == INTEGRATE INWARD TO AVOID EXPONENTIALLY INCREASING FUNCTION
        G=0.D0
        DREL=0.D0
        CALL RADIAL$SCHRODINGER(GID,NR,PSPOT,DREL,0,G,L,E,-1,AUX)
        PSPSI(IRMATCH:,IB-NC)=AUX(IRMATCH:)*PSPSI(IRMATCH,IB-NC)/AUX(IRMATCH)
!
!       == NORMALIZE WAVE FUNCTION
        DO IPRO=1,NPRO(L+1)
          AUX(:)=PSPSI(:,IB-NC)*PRO1(:,IPRO,L+1)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,PROJ(IPRO))
        ENDDO
        AUX(:)=PSPSI(:,IB-NC)**2*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        SVAR=SVAR+DOT_PRODUCT(PROJ(:NPRO(L+1)) &
     &           ,MATMUL(DO1(:NPRO(L+1),:NPRO(L+1),L+1),PROJ(1:NPRO(L+1))))
        PSPSI(:,IB-NC)=PSPSI(:,IB-NC)/SQRT(SVAR)
        PROJ(:)=PROJ(:)/SQRT(SVAR)
        PROJ(:)=MATMUL(DH1(:NPRO(L+1),:NPRO(L+1),L+1)-E*DO1(:NPRO(L+1),:NPRO(L+1),L+1) &
     &                ,PROJ(1:NPRO(L+1)))
        TPSPSI(:,IB-NC)=(E-PSPOT(:)*Y0)*PSPSI(:,IB-NC) &
     &                 -MATMUL(PRO1(:,:NPRO(L+1),L+1),PROJ(:NPRO(L+1)))
      ENDDO
!
!     ====================================================================
!     ==  MAP AUGMENTATION INFORMATION                                  ==
!     ====================================================================
      LNX=SUM(NPRO)
      ALLOCATE(AEPHI(NR,LNX))
      ALLOCATE(PSPHI(NR,LNX))
      ALLOCATE(PRO(NR,LNX))
      ALLOCATE(DTKIN(LNX,LNX))
      ALLOCATE(LOX(LNX))
      DTKIN(:,:)=0.D0
      LN1=0
      DO L1=0,LMAX
        DO IPRO1=1,NPRO(L1+1)
          LN1=LN1+1
          LOX(LN1)=L1
          AEPHI(:,LN1)=AEPHI1(:,IPRO1,L1+1)
          PSPHI(:,LN1)=PSPHI1(:,IPRO1,L1+1)
          PRO(:,LN1)  =PRO1(:,IPRO1,L1+1)
          LN2=0
          DO L2=0,LMAX
            DO IPRO2=1,NPRO(L2+1)
              LN2=LN2+1
              IF(L1.NE.L2) CYCLE
              DTKIN(LN1,LN2)=DTKIN1(IPRO1,IPRO2,L2+1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ====================================================================
!     ==  TOTAL ENERGY                                                  ==
!     ====================================================================
      CALL PSETOT(GID,NR,AEZ,PSRHOC,AERHOC &
     &                 ,NB-NC,LOFI(NC+1:NB),FOFI(NC+1:NB),PSPSI,TPSPSI &
     &                 ,LNX,LOX,PRO,AEPHI,PSPHI,RCSM,DTKIN,VADD)
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE PSETOT(GID,NR,AEZ,PSRHOC,AERHOC,NB,LOFI,FOFI,PSPSI,TPSPSI &
     &                 ,LNX,LOX,PRO,AEPHI,PSPHI,RCSM,DTKIN,VADD)
!     **                                                               **
!     **                                                               **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: AEZ
      REAL(8),   INTENT(IN) :: PSRHOC(NR)
      REAL(8),   INTENT(IN) :: AERHOC(NR)
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      REAL(8),   INTENT(IN) :: FOFI(NB)
      REAL(8)   ,INTENT(IN) :: PSPSI(NR,NB)
      REAL(8)   ,INTENT(IN) :: TPSPSI(NR,NB)
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: PRO(NR,LNX)
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: RCSM
      REAL(8)   ,INTENT(IN) :: DTKIN(LNX,LNX)
      REAL(8)   ,INTENT(IN) :: VADD(NR)
      REAL(8)               :: PI,Y0,C0LL
      REAL(8)               :: R(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR,FAC
      INTEGER(4)            :: NFILO
      REAL(8)               :: DENMAT(LNX,LNX)  ! 1-C DENSITY MATRIX
      REAL(8)               :: PROJ(LNX)    ! PROJECTIONS
      REAL(8)               :: QLM          ! MULTIPOLE MOMENT
      INTEGER(4)            :: IB,LN,LN1,LN2,IR
      REAL(8)               :: PSRHO(NR)
      REAL(8)               :: PSRHO1(NR)
      REAL(8)               :: AERHO1(NR)
      REAL(8)               :: EH1,EXC1
      REAL(8)               :: ALPHA,CL
      REAL(8)               :: PSEKIN,PSEXC,PSEH
      REAL(8)               :: DEKIN,PSEXC1,AEEXC1,PSEH1,AEEH1
      REAL(8)               :: PSEADD,PSEADD1
      REAL(8)               :: ETOT
      REAL(8)               :: POT(NR)
!     *******************************************************************
      PI=4.D0*DATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ===================================================================
!     == DENSITY MATRIX ==================================================
!     ===================================================================
      DENMAT(:,:)=0.D0
      DO IB=1,NB
        PROJ(:)=0.D0
        DO LN=1,LNX
          IF(LOX(LN).NE.LOFI(IB)) CYCLE
          AUX(:)=PSPSI(:,IB)*PRO(:,LN)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PROJ(LN)=SVAR
        ENDDO
        DO LN1=1,LNX
          DO LN2=1,LNX
            DENMAT(LN1,LN2)=DENMAT(LN1,LN2)+FOFI(IB)*PROJ(LN1)*PROJ(LN2)
          ENDDO
        ENDDO
      ENDDO
DO LN=1,LNX
   WRITE(*,FMT='(10F10.5)')DENMAT(LN,:)
ENDDO
!
!     ===================================================================
!     == MULTIPOLE MOMENT=================================================
!     ===================================================================
      AUX(:)=(AERHOC(:)-PSRHOC(:))/Y0
      DO LN1=1,LNX
        DO LN2=1,LNX
          AUX(:)=AUX(:)+DENMAT(LN1,LN2) &
     &            *(AEPHI(:,LN1)*AEPHI(:,LN2)-PSPHI(:,LN1)*PSPHI(:,LN2))
        ENDDO
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      QLM=(SVAR-AEZ)/Y0
!TAKE CARE OF FACTOR YLM!!!!
PRINT*,'COMPENSATION CHARGE ',QLM*Y0
!
!     ===================================================================
!     == KINETIC ENERGY AND CHARGE DENSITY ==============================
!     ===================================================================
      PSRHO(:)=PSRHOC
      AUX(:)=0.D0
      DO IB=1,NB
        PSRHO(:)=PSRHO(:)+FOFI(IB)*PSPSI(:,IB)**2*C0LL
        AUX(:)=AUX(:)+FOFI(IB)*PSPSI(:,IB)*TPSPSI(:,IB)*R(:)**2
      ENDDO
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      PSEKIN=SVAR
PRINT*,'PS-EKIN ',PSEKIN
!
!     == CALCULATE PSEUDO CHARGE FOR TEST
      AUX(:)=4.D0*PI*PSRHO(:)*Y0*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
PRINT*,'PS-CHARGE ',SVAR
!
!     ===================================================================
!     == EXCHANGE AND CORRELATION ENERGY=================================
!     ===================================================================
      CALL MYVOFRHO(GID,NR,0.D0,PSRHO,AUX,EH1,EXC1)
      PSEXC=EXC1
PRINT*,'PS-EXC ',PSEXC
!
!     ===================================================================
!     == HARTREE ENERGY   ===============================================
!     ===================================================================
      ALPHA=1.D0/RCSM**2
      CALL GAUSSN(0,ALPHA,CL)
      FAC=QLM*CL/(4.D0*PI)
CALL RADIAL$INTEGRAL(GID,NR,4.D0*PI*PSRHO*Y0*R**2,SVAR)
PRINT*,'PSEUDO CHARGE ',SVAR
AUX(:)=FAC*EXP(-ALPHA*R(:)**2)
CALL RADIAL$INTEGRAL(GID,NR,4.D0*PI*AUX*Y0*R**2,SVAR)
PRINT*,'COMPENSATION  CHARGE ',SVAR
      AUX(:)=PSRHO(:)+FAC*EXP(-ALPHA*R(:)**2)
CALL RADIAL$INTEGRAL(GID,NR,4.D0*PI*AUX*Y0*R**2,SVAR)
PRINT*,'TOTAL CHARGE ',SVAR
      CALL RADIAL$POISSON(GID,NR,0,AUX,POT)
      AUX(:)=0.5D0*AUX(:)*POT(:)*R(:)**2      
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      PSEH=SVAR
PRINT*,'PS-HARTREE ',PSEH
!
!     ===================================================================
!     == EADD             ===============================================
!     ===================================================================
      AUX(:)=VADD(:)*PSRHO(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      PSEADD=SVAR      
!
!     ===================================================================
!     == ONE-CENTER DENSITIES ===========================================
!     ===================================================================
      AERHO1(:)=AERHOC
      PSRHO1(:)=PSRHOC
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AERHO1(:)=AERHO1(:)+DENMAT(LN1,LN2)*AEPHI(:,LN1)*AEPHI(:,LN2)*C0LL
          PSRHO1(:)=PSRHO1(:)+DENMAT(LN1,LN2)*PSPHI(:,LN1)*PSPHI(:,LN2)*C0LL
        ENDDO
      ENDDO
CALL RADIAL$INTEGRAL(GID,NR,4.D0*PI*PSRHO1*Y0*R**2,SVAR)
PRINT*,'PSEUDO-1 CHARGE ',SVAR
!
!     ===================================================================
!     == KINETIC ENERGY =================================================
!     ===================================================================
      DEKIN=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          DEKIN=DEKIN+DENMAT(LN1,LN2)*DTKIN(LN1,LN2)
        ENDDO
      ENDDO
PRINT*,'AEEKIN1-PSEKIN1 ',DEKIN
!
!     ===================================================================
!     == ONE CENTER EXCHANGE ENERGY    ==================================
!     ===================================================================
      CALL MYVOFRHO(GID,NR,0.D0,PSRHO1,AUX,EH1,EXC1)
      PSEXC1=EXC1
      CALL MYVOFRHO(GID,NR,AEZ,AERHO1,AUX,EH1,EXC1)
      AEEXC1=EXC1
      AEEH1=EH1
PRINT*,'FULL 1-C HARTREE ENERGY ',AEEH1,aez
      CALL MYVOFRHO(GID,NR,AEZ,AERHOC,AUX,EH1,EXC1)
      AEEXC1=AEEXC1-EXC1
      AEEH1=AEEH1-EH1
PRINT*,'CORE 1-C HARTREE ENERGY ',EH1
PRINT*,'VALENCE 1-C HARTREE ENERGY ',AEEH1
PRINT*,'AEEXC1,PSEXC1 ',AEEXC1,PSEXC1
!
!     ===================================================================
!     == HARTREE ENERGIES  ==============================================
!     ===================================================================
      ALPHA=1.D0/RCSM**2
      CALL GAUSSN(0,ALPHA,CL)
      FAC=QLM*CL/(4.D0*PI)
CALL RADIAL$INTEGRAL(GID,NR,4.D0*PI*PSRHO1*Y0*R**2,SVAR)
PRINT*,'PSEUDO-1 CHARGE ',SVAR
      AUX(:)=PSRHO1(:)+FAC*EXP(-ALPHA*R(:)**2)
CALL RADIAL$INTEGRAL(GID,NR,4.D0*PI*AUX*Y0*R**2,SVAR)
PRINT*,'PSEUDO-1+COMP CHARGE ',SVAR
      CALL RADIAL$POISSON(GID,NR,0,AUX,POT)
      AUX(:)=0.5D0*AUX(:)*POT(:)*R(:)**2      
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      PSEH1=SVAR
PRINT*,'AEEH1,PSEH1 ',AEEH1,PSEH1
!
!     ===================================================================
!     == EADD             ===============================================
!     ===================================================================
      AUX(:)=VADD(:)*PSRHO1(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      PSEADD1=SVAR      
PRINT*,'PSEADD,PSEADD1 ',PSEADD,PSEADD1
!
!     ===================================================================
!     == WRITE ===========================================================
!     ===================================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL REPORT$TITLE(NFILO,'PAW-ENERGIES')
      ETOT=PSEKIN+DEKIN &
     &    +PSEH+AEEH1-PSEH1 &
     &    +PSEADD-PSEADD1 &
     &    +PSEXC+AEEXC1-PSEXC1
      CALL REPORT$R8VAL(NFILO,'ETOT',ETOT,'H')
!
      CALL REPORT$TITLE(NFILO,'ALL-ELECTRON ENERGIES')
      CALL REPORT$R8VAL(NFILO,'AE-EKIN',PSEKIN+DEKIN,'H')
      CALL REPORT$R8VAL(NFILO,'AE-EHARTREE',PSEH+AEEH1-PSEH1,'H')
      CALL REPORT$R8VAL(NFILO,'AE-EXC',PSEXC+AEEXC1-PSEXC1,'H')
!
      CALL REPORT$TITLE(NFILO,'PSEUDO ENERGIES')
      CALL REPORT$R8VAL(NFILO,'PS-EKIN',PSEKIN,'H')
      CALL REPORT$R8VAL(NFILO,'PS-EHARTREE',PSEH,'H')
      CALL REPORT$R8VAL(NFILO,'PS-EXC',PSEXC,'H')
      CALL REPORT$R8VAL(NFILO,'PS-EADD',PSEADD,'H')
!
      CALL REPORT$TITLE(NFILO,'PS-PS1 ENERGIES')
      CALL REPORT$R8VAL(NFILO,'PS-PS1-EADD',PSEADD-PSEADD1,'H')
      CALL REPORT$R8VAL(NFILO,'PS-PS1-EHARTREE',PSEH-PSEH1,'H')
      CALL REPORT$R8VAL(NFILO,'PS-PS1-EXC',PSEXC-PSEXC1,'H')
!
      CALL REPORT$TITLE(NFILO,'1C ENERGIES')
      CALL REPORT$R8VAL(NFILO,'PS1-EADD',PSEADD1,'H')
      CALL REPORT$R8VAL(NFILO,'AE1-EHARTREE',AEEH1,'H')
      CALL REPORT$R8VAL(NFILO,'PS1-EHARTREE',PSEH1,'H')
      RETURN
      END
!..................................................................................
MODULE BROYDEN_MODULE
LOGICAL(4)         :: TON=.FALSE.
INTEGER(4)         :: NSTEPX=0
INTEGER(4)         :: NSTEP=0
INTEGER(4)         :: Nx=0
REAL(8)            :: ALPHA
REAL(8),ALLOCATaBLE :: XPREV(:,:)
REAL(8),ALLOCATaBLE :: YPREV(:,:)
END MODULE BROYDEN_MODULE
!      .............................................................................
       subroutine BROYDEN$NEW(NX_,NSTEPX_,ALPHA_)
       USE BROYDEN_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)    :: NX_
       INTEGER(4),INTENT(IN)    :: NSTEPX_
       REAL(8)   ,INTENT(IN)    :: ALPHA_
!      *****************************************************************************
       IF(TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT ALREADY IN USE')
         CALL ERROR$STOP('BROYDEN$NEW')
       END IF
       TON=.TRUE.
       NSTEP=0
       NX=NX_
       NSTEPX=NSTEPX_
       ALPHA=ALPHA_
       ALLOCATE(XPREV(NX,NSTEPX))
       ALLOCATE(YPREV(NX,NSTEPX))
       RETURN
       END
!      .............................................................................
       subroutine BROYDEN$CLEAR
       USE BROYDEN_MODULE
       IMPLICIT NONE
!      *****************************************************************************
       IF(.NOT.TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
         CALL ERROR$STOP('BROYDEN$CLEAR')
       END IF
       TON=.FALSE.
       NSTEPX=0
       NSTEP=0
       NX=0
       ALPHA=0.D0
       DEALLOCATE(XPREV)
       DEALLOCATE(YPREV)
       RETURN
       END
!      .............................................................................
       subroutine BROYDEN$STEP(NX_,X,Y)
       USE BROYDEN_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)    :: NX_
       REAL(8)   ,INTENT(INOUT) :: X(NX_)
       REAL(8)   ,INTENT(IN)    :: Y(NX_)
       REAL(8)   ,ALLOCATABLE   :: DX(:,:)
       REAL(8)   ,ALLOCATABLE   :: DY(:,:)
       REAL(8)   ,ALLOCATABLE   :: B(:,:)
       REAL(8)   ,ALLOCATABLE   :: BINV(:,:)
 REAL(8)   ,ALLOCATABLE   :: w(:,:)
       integer(4)               :: i
!      *****************************************************************************
       IF(.NOT.TON) THEN
         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
         CALL ERROR$STOP('BROYDEN$STEP')
       END IF
       if(nx_.ne.nx) then
         CALL ERROR$MSG('size inconsistent')
         CALL ERROR$STOP('BROYDEN$STEP')
       END IF
!print*,'nstep',nstep
!
!      =================================================================
!      == simple mixing in the first step                             ==
!      =================================================================
       if(nstep.eq.0) then
         if(nstepx.gt.0)then
           nstep=1
           XPREV(:,1)=X(:)     
           YPREV(:,1)=Y(:)     
         END IF
         X=X+ALPHA*Y
         return
       end if
!
!      =================================================================
!      == determine inverse hessian alpha+dx otimes dy                ==
!      =================================================================
       ALLOCATE(DX(NX,NSTEP))
       ALLOCATE(DY(NX,NSTEP))
       DO I=1,NSTEP
         DY(:,I)=YPREv(:,I)-Y(:)  
         DX(:,I)=XPREv(:,I)-X(:)+ALPHA*DY(:,I)
       ENDDO
       ALLOCATE(B(NSTEP,NSTEP))
       ALLOCATE(BINV(NSTEP,NSTEP))
       B=MATMUL(TRANSPOSE(DY),DY)   !OVERLAP MATRIX OF DY
!print*,'b',b
       CALL LIB$INVERTR8(NSTEP,B,BINV)           
allocate(w(nx,nstep))
w=MATMUL(DY,BINV)            !NEW DY IS BIORTHOnormal TO OLD DY
!print*,'w ',matmul(transpose(w),dy)
deallocate(w)
       DY=MATMUL(DY,BINV)            !NEW DY IS BIORTHOnormal TO OLD DY
       DEALLOCATE(B)
       DEALLOCATE(BINV)
!
!      =================================================================
!      == store history                                               ==
!      =================================================================
       IF(NSTEP.LT.NSTEPX)NSTEP=NSTEP+1
       DO I=NSTEP,2,-1
         YPREV(:,I)=YPREV(:,I-1)
         XPREV(:,I)=XPREV(:,I-1)
       ENDDO
       XPREV(:,1)=X(:)     
       YPREV(:,1)=Y(:)     
!
!      =================================================================
!      == predict new vector                                          ==
!      =================================================================
       X=X+ALPHA*Y-MATMUL(DX,MATMUL(TRANSPOSE(DY),Y))
       DEALLOCATE(DX)
       DEALLOCATE(DY)
       RETURN
       END


