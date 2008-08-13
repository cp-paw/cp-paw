! For heavy atoms the code crashes if the grid reached too far out
!
! MIXPOT COULD USE SOME PRECONDITIONING OR BROYDEN
! IT IS NOT CLEAR WHICH REFERENCE ENERGY SHOULD BE TAKEN FOR DREL
! DERIVATIVE DOES NOT HANDLE THE FIRST POINT RIGHT?
! INCLUDE SMALL COMPONENT FOR NORMALIZATION AND DENSITY
!
!   RCUT/RCOV MUST NOT BE TOO SMALL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM TEST
!     **************************************************************************
!     ** NEW ATOMIC PROGRAM FOR THE CALCULATIONS OF ATOMIC SETUPS FOR PAW     **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      logical   ,parameter   :: tphidots=.false. !(does not yet work!)
      TYPE(LL_TYPE)          :: LL_CNTL   ! INPUT LINKED-LIST
      INTEGER(4)             :: NR        ! #(RADIAL GRID POINTS)
      INTEGER(4)             :: GID       ! RADIAL GRID ID
      REAL(8)                :: AEZ       ! ATOMIC NUMBER
      REAL(8)   ,ALLOCATABLE :: AEPOT(:)  ! ALL-ELECTRON ATOM
      REAL(8)   ,ALLOCATABLE :: DREL(:)   ! RELATIVISTIC CORRECTION
      INTEGER(4)             :: NB        ! #(ATOMIC STATES)
      INTEGER(4),ALLOCATABLE :: LOFI(:)
      INTEGER(4),ALLOCATABLE :: SOFI(:)
      INTEGER(4),ALLOCATABLE :: NNOFI(:)
      REAL(8)   ,ALLOCATABLE :: FOFI(:)
      REAL(8)   ,ALLOCATABLE :: EOFI(:)
      LOGICAL   ,ALLOCATABLE :: TCORE(:)
      INTEGER(4)             :: IB,JB,IR
      REAL(8)   ,ALLOCATABLE :: RHOADD(:)
      REAL(8)   ,ALLOCATABLE :: R(:)    ! RADIAL GROD
      REAL(8)   ,ALLOCATABLE :: RCL(:)
      INTEGER(4)             :: IRCUT ! PARTIAL WAVES TRUNCATED AT THIS RADIUS
      REAL(8)                :: RBOX    ! atom is calculated in a box with radius rbox
      REAL(8)                :: Rcut   
      integer(4)             :: irbox
      REAL(8)                :: RCOV
      REAL(8)                :: Y0,PI
      INTEGER(4)             :: NC      ! #(CORE STATES)
      INTEGER(4)             :: NN      ! #(NODES(
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
      INTEGER(4)             :: NNSTART,NNEND
      INTEGER(4)             :: LMAX
      INTEGER(4)             :: NAUG
      INTEGER(4)             :: NAUGX
      INTEGER(4)             :: L,SO
      REAL(8)   ,ALLOCATABLE :: UC(:)   ! HIGHEST CORE STATE FOR THIS L
      REAL(8)                :: E
      REAL(8)   ,ALLOCATABLE :: EGRID(:,:)
      REAL(8)   ,ALLOCATABLE :: DT(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DO(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DH(:,:,:)
      REAL(8)   ,ALLOCATABLE :: AUX(:),aux1(:)
      REAL(8)   ,ALLOCATABLE :: WORK2D(:,:)
      REAL(8)   ,ALLOCATABLE :: G(:)
      REAL(8)   ,ALLOCATABLE :: PHI1(:,:)
      REAL(8)   ,ALLOCATABLE :: PHI2(:,:)
      REAL(8)   ,ALLOCATABLE :: PHI3(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHISAVE(:,:,:)  ! BARE PS PARTIAL WAVE
      REAL(8)   ,ALLOCATABLE :: TPSPHISAVE(:,:,:) 
      REAL(8)   ,ALLOCATABLE :: PROSAVE(:,:,:)    ! BARE PROJECTOR
      REAL(8)   ,ALLOCATABLE :: TRANSPRO(:,:,:)
      REAL(8)   ,ALLOCATABLE :: TRANSPHI(:,:,:)
      REAL(8)   ,ALLOCATABLE :: TRANSPROINV(:,:,:)
      REAL(8)   ,ALLOCATABLE :: TRANSPHIINV(:,:,:)
      REAL(8)   ,ALLOCATABLE :: TRANSPRO1(:,:)
      REAL(8)   ,ALLOCATABLE :: TRANSPHI1(:,:)
      REAL(8)   ,ALLOCATABLE :: PROPSIBAR(:,:,:)
      REAL(8)   ,ALLOCATABLE :: SVEC(:),MAT1(:,:),MAT2(:,:)
      REAL(8)                :: SVAR
      LOGICAL(4)             :: TCHK
      CHARACTER(32)          :: EXT
      INTEGER(4)             :: NFIL
      INTEGER(4)             :: NFILO
      INTEGER(4)             :: IAUG
      REAL(8)                :: EV
      REAL(8)                :: SVAR1,SVAR2
      INTEGER(4)             ::ISVAR1ARR(1)  ! USED TO RESHAPE AN ARRAY OF LENGTH 1
      CHARACTER(16)          :: PSEUDOPARTIALWAVEMETHOD
      character(16)          :: id
      logical(4) ,parameter  :: trel=.false.  !switch between relativistic and non-relativistic
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL CONSTANTS('EV',EV)
!     ==========================================================================
!     == connect standard files to the filehandler (reads input argument)     ==
!     ==========================================================================
      CALL ATTACHFILES()
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL TRACE$SETL4('ON',.false.)
!
!     ==========================================================================
!     == READ INPUT FILE INTO LINKEDLIST LL_CNTL                              ==
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_CNTL)
      CALL FILEHANDLER$UNIT('INPUT',NFIL)
      CALL LINKEDLIST$READ(LL_CNTL,NFIL,'~')
!
!     ==========================================================================
!     == READ CONTROL FILE                                                    ==
!     ==    rbox: atom is calculated in a hard box with radius rbox ,         ==
!     ==          i.e.phi(rbox)=0                                             ==
!     ==    ircut: partial waves are set to zero beyond ircut                 ==
!     ==           default defined by rbox                                    ==
!     ==    rcov : scale factor for input values                              ==
!     ==           used for test of localizatiuon of compensation charge      ==
!     ==           default value for rcl (cutoff for partial wave constr.)    ==
!     ==========================================================================
      IF(TREL) THEN
        CALL REPORT$TITLE(NFILO,'NONRELATIVISTIC CALCULATION')
      ELSE
        CALL REPORT$TITLE(NFILO,'SCALAR RELATIVISTIC CALCULATION')
      END IF
      CALL READCNTL_DIMENSIONS(LL_CNTL,NB,LMAX,NAUGX)
      ALLOCATE(LOFI(NB))
      ALLOCATE(SOFI(NB))
      ALLOCATE(NNOFI(NB))
      ALLOCATE(FOFI(NB))
      ALLOCATE(TCORE(NB))
      CALL READCNTL1(LL_CNTL,NB,AEZ,RCOV,RBOX,rcut,GID,NR,IRCUT,LOFI,SOFI,NNOFI,FOFI,TCORE)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
!
!     == irbox is the grid point just outside rbox =============================
      do ir=1,nr
        irbox=ir
        if(r(ir).gt.rbox) exit
      enddo
!
!     == count #(core states) ==================================================
!     == core states preceed valence states (see readcntl1) ====================
      NC=0
      DO IB=1,NB
        IF(TCORE(IB)) NC=NC+1
      ENDDO
!
!     ==========================================================================
!     == PERFORM SELF CONSISTENT ALL-ELECTRON CALCULATION                     ==
!     ==========================================================================
      ALLOCATE(EOFI(NB)) 
      ALLOCATE(DREL(NR))
      ALLOCATE(AEPOT(NR))      
      ALLOCATE(RHOADD(NR))
!
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AEPOT)
      RHOADD=0.D0
call trace$pass('all-electron scf calculation')
      CALL AESCF(GID,NR,NB,trel,LOFI,SOFI,FOFI,NNOFI,AEZ,RHOADD,rbox,DREL,AEPOT,EOFI)
!
!     == REPORT ================================================================
      CALL REPORT$TITLE(NFILO,'EIGENSTATES AFTER SCF CALCULATION')
      DO IB=1,NB
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
     &                ," E[H]=",F20.9," E[EV]=",F20.9)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/EV
      ENDDO
      DEALLOCATE(RHOADD)
!
!     ==========================================================================
!     == CALCULATE CORE AND VALENCE DENSITY                                   ==
!     ==========================================================================
      ALLOCATE(AERHOC(NR))
      ALLOCATE(AERHOV(NR))
      IF(TREL) THEN
        ID='FULL'
      ELSE
        ID='NONREL'
        drel=0.d0
      END IF
      CALL AERHO(id,GID,NR,NC,LOFI(1:NC),SOFI(1:NC),FOFI(1:NC),NNOFI(1:NC) &
     &          ,rbox,DREL,AEPOT,AERHOC,EOFI(1:NC))
!
!     == RELATIVISTIC CORRECTION WILL BE USED FROM THE VALENCE STATES
call trace$pass('determine valence density')
      CALL AERHO(id,GID,NR,NB-NC,LOFI(NC+1:NB),SOFI(NC+1:NB),FOFI(NC+1:NB) &
     &          ,NNOFI(NC+1:NB),rbox,DREL,AEPOT,AERHOV,EOFI(NC+1:NB))
!
!     == REPORT ================================================================
      WRITE(NFILO,FMT='("EIGENSTATES AFTER NON-SCF CALCULATION IN SCF POT")')
      DO IB=1,NB
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
     &                ," E[H]=",F20.9," E[EV]=",F20.9)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/EV
      ENDDO
!
!     ==========================================================================
!     == CALCULATE NODELESS WAVE FUNCTIONS                                    ==
!     ==========================================================================
      ALLOCATE(UOFI(NR,NB))  ; UOFI=0.D0
      ALLOCATE(TUOFI(NR,NB)) ; TUOFI=0.D0
      IF(TREL) THEN
        ID='EFFZORA'
      ELSE
        ID='NONREL'
      END IF
      CALL NODELESS(ID,GID,NR,NB,LOFI,SOFI,NNOFI,RBOX,DREL,AEPOT &
     &             ,UOFI,TUOFI,EOFI)
!
!     == REPORT ================================================================
!     == SMALL DEVIATIONS FOR THE HIGHEST STATES ARE PRESENT IF THE BOX RADIUS =
!     == RBOX IS CHOSEN TOO SMALL. =============================================
      WRITE(NFILO,FMT='("EIGENSTATES FROM NODELESS CALCULATION")')
      DO IB=1,NB
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
     &                ," E[H]=",F20.9," E[EV]=",F20.9)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/EV
      ENDDO
!     CALL HAMTEST(GID,NR,AEPOT,NB,LOFI,SOFI,UOFI,TUOFI)
!
!     ==========================================================================
!     == SELF-CONSISTENT ATOM FINISHED; WRITE RESULT                          ==
!     ==========================================================================
      CALL WRITEPHI('ALLNODELESS.DAT',GID,NR,NB,UOFI)
      CALL WRITEPHI('ALLNODELESS_T.DAT',GID,NR,NB,TUOFI)
!
!     ==========================================================================
!     == TOTAL ENERGY                                                         ==
!     ==========================================================================
      CALL AEETOT(GID,NR,NB,AEZ,NC,EOFI,FOFI,AERHOC,AERHOV,AEPOT)
!
!     ==========================================================================
!     ==========================================================================
!     == DEFINE AUGMENTATION                                                  ==
!     ==========================================================================
!     ==========================================================================
!
!     ==========================================================================
!     ==  read information from control input file                            ==
!     ==========================================================================
      ALLOCATE(NPRO(LMAX+1))
      ALLOCATE(RCL(LMAX+1))
      CALL READCNTL2(LL_CNTL,LMAX,RCOV,RC_CORE,POW_CORE,VAL0_CORE &
     &              ,RC_POT,POW_POT,VAL0_POT,RCSM,PSEUDOPARTIALWAVEMETHOD,NPRO,RCL)
!
!     ==========================================================================
!     == CONSTRUCT PSEUDO CORE                                                ==
!     ==========================================================================
      ALLOCATE(PSRHOC(NR))
      CALL PSEUDIZE(GID,NR,POW_CORE,VAL0_CORE,RC_CORE,AERHOC,PSRHOC)
!
!     ==========================================================================
!     ==  CONSTRUCT PSEUDO POTENTIAL                                          ==
!     ==========================================================================
      ALLOCATE(PSPOT(NR))
      CALL PSEUDIZE(GID,NR,POW_POT,VAL0_POT,RC_POT,AEPOT,PSPOT)
!
!     ==========================================================================
!     == ALLOCATE ARRAYS                                                      ==
!     ==========================================================================
      ALLOCATE(PSRHOV(NR))
      ALLOCATE(VADD(NR))
      ALLOCATE(UC(NR))
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
!
      NAUGX=MAXVAL(NPRO)
      ALLOCATE(AEPHI(NR,NAUGX,LMAX+1))  ; AEPHI=0.D0
      ALLOCATE(TAEPHI(NR,NAUGX,LMAX+1)) ; TAEPHI=0.D0
      ALLOCATE(PSPHI(NR,NAUGX,LMAX+1))  ; PSPHI=0.D0
      ALLOCATE(TPSPHI(NR,NAUGX,LMAX+1)) ; TPSPHI=0.D0
      ALLOCATE(PRO(NR,NAUGX,LMAX+1))    ; PRO=0.D0
      ALLOCATE(UPHI(NR,NAUGX,LMAX+1))   ; UPHI=0.D0
      ALLOCATE(TUPHI(NR,NAUGX,LMAX+1))  ; TUPHI=0.D0
      ALLOCATE(TRANSPRO(NAUGX,NAUGX,LMAX+1))    ; TRANSPRO=0.D0
      ALLOCATE(TRANSPHI(NAUGX,NAUGX,LMAX+1))    ; TRANSPHI=0.D0
      ALLOCATE(TRANSPROINV(NAUGX,NAUGX,LMAX+1)) ; TRANSPROINV=0.D0
      ALLOCATE(TRANSPHIINV(NAUGX,NAUGX,LMAX+1)) ;TRANSPHIINV=0.D0
      ALLOCATE(DT(NAUGX,NAUGX,LMAX+1))
      ALLOCATE(DO(NAUGX,NAUGX,LMAX+1))
      ALLOCATE(DH(NAUGX,NAUGX,LMAX+1))
      ALLOCATE(PROPSIBAR(NAUGX,NAUGX,LMAX+1))
      ALLOCATE(EGRID(NAUGX,LMAX+1))
      ALLOCATE(PSPHISAVE(NR,NAUGX,LMAX+1))
      ALLOCATE(TPSPHISAVE(NR,NAUGX,LMAX+1))
      ALLOCATE(PROSAVE(NR,NAUGX,LMAX+1))   
!
!     ==========================================================================
!     == CONSTRUCT PARTIAL WAVES AND PROJECTORS FOR EACH L                    ==
!     ==========================================================================
      PSRHOV(:)=0.D0
      DO L=0,LMAX
        WRITE(NFILO,FMT='(30("="),2X,"L=",I1,2X,30("="))')L
        NAUG=NPRO(L+1)
!       == EXT IS USED AS ENDING FOR PRINTOUTS OF L-DEPENDENT FUNCTIONS ========
        IF(L.EQ.0) THEN
          EXT='_S.DAT'
        ELSE IF (L.EQ.1) THEN
          EXT='_P.DAT'
        ELSE IF (L.EQ.2) THEN
          EXT='_D.DAT'
        ELSE IF (L.EQ.3) THEN
          EXT='_F.DAT'
        ELSE
          CALL ERROR$MSG('L>3 IS NOT ALLOWED YET')
          CALL ERROR$STOP('TEST')
        END IF
!
!ATTENTION CURRENTLY NPRO MUST BE THE SAME FOR ALL L!!!!
!
!       ========================================================================
!       ==  DETERMINE HIGHEST NODELESS CORE WAVE FUNCTION UC                  ==
!       ========================================================================
        UC(:)=0.D0
!
!       == DETERMINE HIGHEST CORE STATE UC
!       == DETERMINE NNSTART, THE NUMBER OF NODES FOR THE LOWEST VALENCE STATE
        NNSTART=0   ! #(NODES OF FIRST  VALENCE WAVE FUNCTION)
        DO IB=1,NC
          IF(LOFI(IB).NE.L) CYCLE
          IF(NNOFI(IB)+1.GT.NNSTART) THEN
            UC(:)=UOFI(:,IB)
            NNSTART=NNOFI(IB)+1  ! LOWEST VALENCE STATE HAS ONE NODE MORE 
                                 ! THAN HIGHEST CORE STATE
          END IF
        ENDDO        
!
!       ==  DETERMINE ENERGIES FOR PARTIAL WAVE CONSTRUCTION ===================
        EGRID(:,L+1)=0.D0
        TCHK=.FALSE.
        NNEND=0
        DO IB=NC+1,NB
          IF(LOFI(IB).NE.L) CYCLE
          IF(NNOFI(IB).GE.NNSTART) THEN
            IAUG=NNOFI(IB)-NNSTART+1   !IDENTIFY IAUG WITH NUMBER OF NODES +1
            IF(IAUG.GT.NAUG) THEN
              CALL ERROR$MSG('IAUG OUT OF RANGE')
              CALL ERROR$I4VAL('IB',IB)
              CALL ERROR$I4VAL('L',L)
              CALL ERROR$I4VAL('NN',NNOFI(IB))
              CALL ERROR$I4VAL('IAUG',IAUG)
              CALL ERROR$I4VAL('NAUG',NAUGX)
              CALL ERROR$STOP('MAIN')
            END IF 
            EGRID(IAUG,L+1)=EOFI(IB)
            NNEND=MAX(IAUG,NNEND)   ! egrid(1:nnend) corresponds to bound states
            TCHK=.TRUE.
          END IF
        ENDDO
!
!       == CHECK REQUIREMENT FROM BIORTHOMATRICES_TAYLOR ========================
        IF(TPHIDOTS.AND.NNEND+1.GT.NAUG) THEN
          CALL ERROR$MSG('NNEND+1>NAUG')
          CALL ERROR$MSG('AT LEAST ONE SCATTERING STATE IS REQUIRED IF TPHIDOTS=.TRUE.')
          CALL ERROR$L4VAL('TPHIDOTS',TPHIDOTS)
          CALL ERROR$I4VAL('IB',IB)
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NNEND',NNEND)
          CALL ERROR$I4VAL('NAUG',NAUGX)
          CALL ERROR$STOP('MAIN')
        END IF 
!
!       ==  THE LOWEST STATE OF EACH ANGULAR MOMENTUM SHALL BE A BOUND STATE ===
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('LOWEST VALENCE STATE NOT FOUND')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('MAIN')
        END IF
!
        CALL REPORT$I4VAL(NFILO,'NN OF HIGHEST CORE STATES ',NNSTART-1,'')
!
!       ========================================================================
!       ==  DETERMINE NODELESS ENERGY EXPANSION COEFFICIENTS                  ==
!       ========================================================================
!       USE RELATIVISTIC CORRECTIONS FROM ABOVE
        SO=0
        IF(TPHIDOTS) THEN
          ALLOCATE(G(NR))
          DO IB=1,NAUG
            NN=0
            IF(IB.EQ.1) THEN
              G(:)=UC
            ELSE
              G(:)=UPHI(:,IB-1,L+1)
            END IF
            IF(IB.LE.NNEND) THEN
              EGRID(IB,L+1)=EGRID(IB,L+1)
              CALL ONENODELESS('BOUND',GID,NR,L,NN,SO,RBOX,DREL,AEPOT,G &
       &                      ,EGRID(IB,L+1),UPHI(:,IB,L+1),TUPHI(:,IB,L+1))
              CALL REPORT$R8VAL(NFILO,'ENERGY SUPPORT GRID POINT (BOUND)' & 
       &                           ,EGRID(IB,L+1),'H')
            ELSE
              EGRID(IB,L+1)=EGRID(NNEND,L+1)
              CALL ONENODELESS('SCATT',GID,NR,L,NN,SO,RBOX,DREL,AEPOT,G &
       &                      ,EGRID(NNEND,L+1),UPHI(:,IB,L+1),TUPHI(:,IB,L+1))
            END IF
          ENDDO  
          DEALLOCATE(G)
          CALL WRITEPHI('XXX'//TRIM(EXT),GID,NR,NAUG,UPHI(:,:,L+1))
        ELSE if(1.eq.0) then
          DO IB=1,NAUG
            IF(IB.LE.NNEND) THEN
              NN=IB-1
              CALL ONENODELESS('BOUND',GID,NR,L,NN,SO,RBOX,DREL,AEPOT,UC &
       &                    ,EGRID(IB,L+1),UPHI(:,IB,L+1),TUPHI(:,IB,L+1))
              CALL REPORT$R8VAL(NFILO,'ENERGY SUPPORT GRID POINT (BOUND)' & 
       &                           ,EGRID(IB,L+1),'H')
            ELSE
              IF(IB.EQ.1.AND.NNEND.GT.0) THEN   
                EGRID(IB,L+1)=0.D0
              ELSE
                EGRID(IB,L+1)=EGRID(NNEND,L+1)+0.5D0*REAL(IB-NNEND,KIND=8)**2
              END IF
              CALL ONENODELESS('SCATT',GID,NR,L,NN,SO,RBOX,DREL,AEPOT,UC &
       &                      ,EGRID(IB,L+1),UPHI(:,IB,L+1),TUPHI(:,IB,L+1))
              CALL REPORT$R8VAL(NFILO,'ENERGY SUPPORT GRID POINT (UNBOUND)' &
       &                           ,EGRID(IB,L+1),'H')
            END IF
          ENDDO
        ELSE 
          DO IB=1,NAUG
            NN=IB-1
            CALL ONENODELESS('BOUND',GID,NR,L,NN,SO,RBOX,DREL,AEPOT,UC &
       &                    ,EGRID(IB,L+1),UPHI(:,IB,L+1),TUPHI(:,IB,L+1))
            CALL REPORT$R8VAL(NFILO,'ENERGY SUPPORT GRID POINT (BOUND)' & 
       &                           ,EGRID(IB,L+1),'H')
          ENDDO
        END IF
!
!       ========================================================================
!       ==  ALL ELECTRON FUNCTIONS BY CORE ORTHOGONALIZATION                  ==
!       ========================================================================
        CALL COREORTHOGONALIZE(GID,NR,rbox,NC,LOFI(1:NC),SOFI(1:NC) &
      &                      ,UOFI(:,1:NC),TUOFI(:,1:NC) &
      &                      ,L,NAUG,UPHI(:,:naug,L+1),TUPHI(:,:naug,L+1) &
      &                      ,AEPHI(:,:naug,L+1),TAEPHI(:,:naug,L+1))
!
!       ========================================================================
!       ==  CHOSE NORM OF PARTIAL WAVES                                       ==
!       ========================================================================
        IF(TPHIDOTS) THEN
          CALL RADIAL$INTEGRAL(GID,NR,(AEPHI(:,1,L+1)*R(:))**2,SVAR)
          SVAR=1.D0/SQRT(SVAR)
          UPHI(:,:,L+1)  =  UPHI(:,:,L+1)*SVAR
          TUPHI(:,:,L+1) = TUPHI(:,:,L+1)*SVAR
          AEPHI(:,:,L+1) = AEPHI(:,:,L+1)*SVAR
          TAEPHI(:,:,L+1)=TAEPHI(:,:,L+1)*SVAR
        ELSE
          DO IB=1,NAUG
!           __ ONLY BOUND STATES CAN BE NORMALIZED. ______________________________
!           __ FOR UNBOUND STATES THE BOUNDARY CONDITIONS AT THE NUCLEUS _________
!           __ ARE KEPT IDENTICAL TO THAT OF THE HIGHEST BOUND STATE _____________
            IF(IB.LE.NNEND) THEN
              CALL RADIAL$INTEGRAL(GID,NR,(AEPHI(:,IB,L+1)*R(:))**2,SVAR)
              SVAR=1.D0/SQRT(SVAR)
            ELSE
!             == FOR ENERGY-TAYLOR EXPANSION THE SCALING MUST AGREE WITH THE LAST
!             == NODELESS BOUND STATE
!             == CHOOSE 5TH POINT TO AVOID DIVIDING BY SMALL NUMBERS FOR L.NEQ.0 =
              IF(IB.EQ.1) THEN
                SVAR=R(5)**L/AEPHI(5,IB,L+1)
              ELSE
                SVAR=AEPHI(5,IB-1,L+1)/AEPHI(5,IB,L+1)
              END IF
            END IF
            UPHI(:,IB,L+1)  =  UPHI(:,IB,L+1)*SVAR
            TUPHI(:,IB,L+1) = TUPHI(:,IB,L+1)*SVAR
            AEPHI(:,IB,L+1) = AEPHI(:,IB,L+1)*SVAR
            TAEPHI(:,IB,L+1)=TAEPHI(:,IB,L+1)*SVAR
          ENDDO
        END IF
!
!       ========================================================================
!       ==  CONSTRUCT PS PARTIAL WAVES                                        ==
!       ========================================================================
        PSPHI(:,:,L+1)=UPHI(:,:,L+1)
        TPSPHI(:,:,L+1)=TUPHI(:,:,L+1)
 PSPHI(:,:,L+1)=aePHI(:,:,L+1)
 TPSPHI(:,:,L+1)=TaePHI(:,:,L+1)
        RC=RCL(L+1)
!       USE METHOD ADJUSTING PSEUDOPOTENTIAL
!PSEUDOPARTIALWAVEMETHOD='HBS'
PSEUDOPARTIALWAVEMETHOD='KRESSE'
        IF(PSEUDOPARTIALWAVEMETHOD.EQ.'HBS') THEN
          CALL CONSTRUCTPSPHI1(GID,NR,RC,POW_POT,PSPOT,L,NAUG &
     &                      ,EGRID(:,L+1),PSPHI(:,:,L+1),TPSPHI(:,:,L+1))
!       == USE METHOD MATCHING BESSEL FUNCTIONS
        ELSE IF(PSEUDOPARTIALWAVEMETHOD.EQ.'KRESSE') THEN
          if(tphidots) then
            CALL CONSTRUCTPSPHI(GID,NR,RC,L,1,PSPHI(:,1,L+1),TPSPHI(:,1,L+1))
          else
            CALL CONSTRUCTPSPHI(GID,NR,RC,L,NAUG,PSPHI(:,:,L+1),TPSPHI(:,:,L+1))
          end if
        ELSE
          CALL ERROR$MSG('METHOD FOR CONSTRUCTING PSEUDO PARTIAL WAVES')
          CALL ERROR$MSG('NOT RECOGNIZED')
          CALL ERROR$CHVAL('PSEUDOPARTIALWAVEMETHOD',PSEUDOPARTIALWAVEMETHOD)
          CALL ERROR$STOP('ATOM')
        END IF
!
!       ========================================================================
!       ==  WRITE PSEUDOPOTENTIAL FOR EACH PARTIAL WAVE                       ==
!       ========================================================================
        ALLOCATE(WORK2D(NR,NAUG))
        WORK2D=0.D0
        DO IR=1,IRBOX
          DO IB=1,NAUG
            IF(PSPHI(IR,IB,L+1).NE.0.D0) THEN
              WORK2D(IR,IB)=-TPSPHI(IR,IB,L+1)/PSPHI(IR,IB,L+1)+EGRID(IB,L+1)
            ELSE
             WORK2D(IR,IB)=0.D0
           END IF
          ENDDO
        ENDDO
        CALL WRITEPHI('PSPOTOFPHI'//TRIM(EXT),GID,NR,NAUG,WORK2D)
        DEALLOCATE(WORK2D)
!
!       ========================================================================
!       ==  CONSTRUCT BARE PROJECTOR FUNCTIONS                                ==
!       ========================================================================
        IF(TPHIDOTS) THEN
          CALL BAREPROJECTORS_TAYLOR(GID,NR,L,EGRID(:,L+1),PSPOT,10.D0,NAUG &
       &             ,NNEND,PSPHI(:,:,L+1),TPSPHI(:,:,L+1),PRO(:,:,L+1))
        ELSE
          CALL BAREPROJECTORS(GID,NR,L,EGRID(:,L+1),PSPOT,10.D0,NAUG &
       &                   ,PSPHI(:,:,L+1),TPSPHI(:,:,L+1),PRO(:,:,L+1))
        END IF
!
!       ========================================================================
!       ==  CONSTRUCT 1C KINETIC ENERGY AND OVERLAP                           ==
!       ========================================================================
        CALL DTDO(GID,NR,rbox,NAUG,AEPHI(:,:naug,L+1),TAEPHI(:,:naug,L+1) &
       &         ,PSPHI(:,:naug,L+1),TPSPHI(:,:naug,L+1),DT(:naug,:naug,L+1),DO(:naug,:naug,L+1))
        CALL MAKEDATH(GID,NR,rbox,AEPOT,PSPOT,L,NAUG,AEPHI(:,:naug,L+1),PSPHI(:,:naug,L+1) &
       &             ,DT(:naug,:naug,L+1),DH(:naug,:naug,L+1))
!
!       ========================================================================
!       == ADD TO PSEUDO DENSITY (MUST BE DONE BEFORE BIORTHOGONALIZATION)    ==
!       == ASSUMES THAT THERE IS ONE PARTIAL WAVE FOR EACH OCCUPIED STATE     ==
!       ========================================================================
        DO IB=NC+1,NB
          IF(LOFI(IB).NE.L) CYCLE
          IAUG=NNOFI(IB)-NNSTART+1
          AUX(:)=(AEPHI(:,IAUG,L+1)*R(:))**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
          SVAR=FOFI(IB)/SVAR*Y0
          PSRHOV(:)=PSRHOV(:)+SVAR*PSPHI(:,IAUG,L+1)**2
        ENDDO
!
!       ========================================================================
!       ==  TEST RESULT BY RECALCULATING PSPHI FROM BARE PROJECTORS           ==
!       ========================================================================
        WRITE(NFILO,FMT='(72("=")/72("="),T20,"  BARE QUANTITIES DONE:' &
       &                                              //' TESTS...  "/72("="))')
        ALLOCATE(PHI1(NR,NAUG))
        ALLOCATE(PHI2(NR,NAUG))
        ALLOCATE(PHI3(NR,NAUG))
        ALLOCATE(G(NR))
        WRITE(NFILO,FMT='("COMPARE SOLUTION OF SCHROEDINGER EQUATIONS' &
       &                                    //' WITH ORIGINAL PARTIAL WAVES:")')
        WRITE(NFILO,FMT='("LEFT: ALL-ELECTRON; ' &
       &                              //'[T+AEPOT-E|AEF>=0 -> MAX|AEF-AEPHI|")')
        WRITE(NFILO,FMT='("RIGHT: PSEUDO; ' &
       &                          //'[T+PSPOT-E|PSF>=|PRO> -> MAX|PSF-PSPHI|")')
        WRITE(NFILO,FMT='("BOTH MUST VANISH: ' &
       &            //'if not see TESTBAREAEPSI_*.DAT OR TESTBAREPSPSI_*.DAT")')
        DO IB=1,NAUG
!         == all electron partial wave =========================================
          G=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,0,G,L,EGRID(IB,L+1) &
       &                        ,1,PHI1(:,IB))
!         == NORMALIZE ALL-ELECTRON PARTIAL WAVE ===============================
          CALL RADIAL$INTEGRATE(GID,NR,R**2*AEPHI(:,IB,L+1)**2,AUX)
          CALL RADIAL$VALUE(GID,NR,AUX,RBOX,SVAR1)
          CALL RADIAL$INTEGRATE(GID,NR,R**2*PHI1(:,IB)**2,AUX)
          CALL RADIAL$VALUE(GID,NR,AUX,RBOX,SVAR2)
          PHI1(:,IB)=PHI1(:,IB)*SQRT(SVAR1/SVAR2)

!         == INHOMOGENEOUS PSEUDO PARTIAL WAVE =================================
          G(:)=PRO(:,IB,L+1)
          CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,0,G,L,EGRID(IB,L+1) &
       &                        ,1,PHI2(:,IB))
!         == HOMOGENEOUS PSEUDO PARTIAL WAVE ===================================
          G=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,0,G,L,EGRID(IB,L+1) &
       &                        ,1,PHI3(:,IB))
!
!         == NORMALIZE PSEUDO PARTIAL WAVE =====================================
          PHI2(IRBOX:,IB)=0.D0
          PHI3(IRBOX:,IB)=0.D0
          ISVAR1ARR=MAXLOC(ABS(PSPHI(:IRBOX-1,IB,L+1)-PHI2(:IRBOX-1,IB)))
          IR=MIN(ISVAR1ARR(1),IRBOX-1)
          IF(ABS(PHI3(IR,IB)).GT.1.D-20) THEN
            SVAR=(PSPHI(IR,IB,L+1)-PHI2(IR,IB))/PHI3(IR,IB)
            PHI2(:,IB)=PHI2(:,IB)+PHI3(:,IB)*SVAR
          ENDIF
          WRITE(NFILO,FMT='(I5,2E20.10)')IB &
       &         ,MAXVAL(ABS(PHI1(:IRBOX-1,IB)-AEPHI(:IRBOX-1,IB,L+1))) &
       &         ,MAXVAL(ABS(PHI2(:IRBOX-1,IB)-PSPHI(:IRBOX-1,IB,L+1)))
        ENDDO
        DEALLOCATE(G)
!
        ALLOCATE(WORK2D(NR,3*NAUG))
        WORK2D(:,1:NAUG)=PHI2(:,:naug)
        WORK2D(:,NAUG+1:2*NAUG)=PSPHI(:,:naug,L+1)
        WORK2D(:,2*NAUG+1:3*NAUG)=PSPHI(:,:naug,L+1)-PHI2(:,:naug)
        CALL WRITEPHI('TESTBAREPSPHI'//TRIM(EXT),GID,NR,3*NAUG,WORK2D)
        WORK2D(:,1:NAUG)=PHI1(:,:naug)
        WORK2D(:,NAUG+1:2*NAUG)=AEPHI(:,:naug,L+1)
        WORK2D(:,2*NAUG+1:3*NAUG)=AEPHI(:,:naug,L+1)-PHI1(:,:naug)
        CALL WRITEPHI('TESTBAREAEPHI'//TRIM(EXT),GID,NR,3*NAUG,WORK2D)
        DEALLOCATE(WORK2D)
!
        DEALLOCATE(PHI3) 
        DEALLOCATE(PHI2)
        DEALLOCATE(PHI1)
!
!       ========================================================================
!       ==  TEST SCHROEDINGER EQUATION FOR AEPHI AND PSPHI                    ==
!       ========================================================================
        ALLOCATE(PHI2(NR,NAUG))
        ALLOCATE(PHI3(NR,NAUG))
!       == PHI2 AND PHI3 MUST VANISH IF EVERYTHING IS OK =======================
        DO IB=1,NAUG
          PHI2(:,IB)=TAEPHI(:,IB,L+1) &
      &             +(AEPOT(:)*Y0-EGRID(IB,L+1))*AEPHI(:,IB,L+1)
          PHI3(:,IB)=TPSPHI(:,IB,L+1) &
      &             +(PSPOT(:)*Y0-EGRID(IB,L+1))*PSPHI(:,IB,L+1) &
      &             -PRO(:,IB,L+1)
          DO JB=1,NAUG
            AUX(:)=PRO(:,IB,L+1)*PSPHI(:,JB,L+1)*R(:)**2
            CALL RADIAL$INTEGRAte(GID,NR,AUX,aux1)
            CALL RADIAL$value(GID,NR,AUX1,rbox,svar)
            PROPSIBAR(IB,JB,L+1)=SVAR
          ENDDO
        ENDDO
!
!       == TEST IF [T+VTILDE-E+|PRO>(DH-EDO)<PRO|]|PHI>=0 ======================
        WRITE(NFILO,FMT='("TEST SCHROEDINGER EQUATION FOR BARE WAVE FUNCTIONS")')
        WRITE(NFILO,FMT='("....THE FOLLOWING VALUES MUST VANISH")')
        WRITE(NFILO,FMT='(T11,"[T+V-E]|AEPHI>",T31,"[T+PSV-E]|PSPHI>-|PRO>")') 
        DO IB=1,NAUG
          SVAR1=MAXVAL(ABS(PHI2(:,IB)))
          SVAR2=MAXVAL(ABS(PHI3(:,IB)))
          WRITE(NFILO,FMT='(I5,3E20.10)')IB,SVAR1,SVAR2
        ENDDO
        CALL WRITEPHI('TESTAEPHISCHRGL'//TRIM(EXT),GID,NR,NAUG,PHI2)
        CALL WRITEPHI('TESTPSPHISCHRGL'//TRIM(EXT),GID,NR,NAUG,PHI3)
!
        WRITE(NFILO,FMT='("(DH-EDO<PRO|PSPSI>+<PRO|PSPSI> ... MUST VANISH")')
        DO IB=1,NAUG
          WRITE(NFILO,FMT='(20F20.10)')DH(:naug,IB,L+1)-EGRID(IB,L+1)*DO(:naug,IB,L+1)+PROPSIBAR(IB,:,L+1)
        ENDDO
        WRITE(NFILO,FMT='("<PRO|PSPSI> (CAN BE ANYTHING)")')
        DO IB=1,NAUG
          WRITE(NFILO,FMT='(20F20.10)')PROPSIBAR(IB,:NAUG,L+1)
        ENDDO
        DEALLOCATE(PHI2)
        DEALLOCATE(PHI3)
        WRITE(NFILO,FMT='(72("=")/72("="),T20,"TESTS DONE"/72("="))')
!
!       ========================================================================
!       ==  keep bare partial waves for test purposes                         ==
!       ========================================================================
        PSPHISAVE(:,:,L+1)=PSPHI(:,:,L+1)
        TPSPHISAVE(:,:,L+1)=TPSPHI(:,:,L+1)
        PROSAVE(:,:,L+1)=PRO(:,:,L+1)
!
!       ========================================================================
!       == TRUNCATE PARTIAL WAVES TO AVOID NUMERICAL ERRORS WITH EXPONENTIAL  ==
!       == INCREASE                                                           ==
!       ========================================================================
        UPHI(IRCUT:,:,L+1)=0.D0
        TUPHI(IRCUT:,:,L+1)=0.D0
        AEPHI(IRCUT:,:,L+1)=0.D0
        TAEPHI(IRCUT:,:,L+1)=0.D0
        PSPHI(IRCUT:,:,L+1)=0.D0
        TPSPHI(IRCUT:,:,L+1)=0.D0
        PRO(IRCUT:,:,L+1)=0.D0
!
!       ========================================================================
!       == BARE QUANTITIES FINISHED. WRITE RESULT...                          ==
!       ========================================================================
!       == BAREAUGMENT CONTAINS AE,NODELESS, AND PS PARTIAL WAVES
        ALLOCATE(WORK2D(NR,3*NAUG))
        WORK2D(:,1:NAUG)=AEPHI(:,:naug,L+1)
        WORK2D(:,NAUG+1:2*NAUG)=UPHI(:,:naug,L+1)
        WORK2D(:,2*NAUG+1:3*NAUG)=PSPHI(:,:naug,L+1)
        CALL WRITEPHI('BAREPARTIALWAVES'//TRIM(EXT),GID,NR,3*NAUG,WORK2D)
        DEALLOCATE(WORK2D)
        CALL WRITEPHI('BAREPRO'//TRIM(EXT),GID,NR,NAUG,PRO(:,:naug,L+1))
!
!       ========================================================================
!       ==  ENSURE BIORTHOGONALIZATION CONDITION                              ==
!       ========================================================================
        allocate(transphi1(naug,naug))
        allocate(transpro1(naug,naug))
        if(tphidots) then
          CALL BIORTHOMATRICES_Taylor(GID,NR,L,NAUG,nnend,rcov &
       &                  ,PSPHI(:,:naug,L+1),PRO(:,:naug,L+1) &
       &                  ,TRANSPHI1,TRANSPRO1)
        else
          CALL BIORTHOMATRICES(GID,NR,rbox,L,NAUG,PSPHI(:,:NAUG,L+1),PRO(:,:NAUG,L+1) &
       &                   ,TRANSPHI1,TRANSPRO1)
        END IF
!
        PRO(:,:naug,L+1)   =MATMUL(PRO(:,:naug,L+1),TRANSPRO1)
        PSPHI(:,:naug,L+1) =MATMUL(PSPHI(:,:naug,L+1),TRANSPHI1)
        TPSPHI(:,:naug,L+1)=MATMUL(TPSPHI(:,:naug,L+1),TRANSPHI1)
        AEPHI(:,:naug,L+1) =MATMUL(AEPHI(:,:naug,L+1),TRANSPHI1)
        TAEPHI(:,:naug,L+1)=MATMUL(TAEPHI(:,:naug,L+1),TRANSPHI1)
        UPHI(:,:naug,L+1) =MATMUL(UPHI(:,:naug,L+1),TRANSPHI1)
        TUPHI(:,:naug,L+1)=MATMUL(TUPHI(:,:naug,L+1),TRANSPHI1)
!
        DT(:naug,:naug,L+1)=MATMUL(TRANSPOSE(TRANSPHI1) &
       &                ,MATMUL(DT(:naug,:naug,L+1),TRANSPHI1))
        DO(:naug,:naug,L+1)=MATMUL(TRANSPOSE(TRANSPHI1) &
       &                ,MATMUL(DO(:naug,:naug,L+1),TRANSPHI1))
        DH(:naug,:naug,L+1)=MATMUL(TRANSPOSE(TRANSPHI1) &
       &                ,MATMUL(DH(:naug,:naug,L+1),TRANSPHI1))
!
        transphiinv(:,:,l+1)=0.d0
        transproinv(:,:,l+1)=0.d0
        CALL LIB$INVERTR8(NAUG,TRANSPHI1,TRANSPHIINV(:naug,:naug,L+1))
        CALL LIB$INVERTR8(NAUG,TRANSPRO1,TRANSPROINV(:naug,:naug,L+1))
        transphi(:,:,l+1)=0.d0
        transpro(:,:,l+1)=0.d0
        transphi(:naug,:naug,l+1)=transphi1(:,:)
        transpro(:naug,:naug,l+1)=transpro1(:,:)
!
!       == test biorthogonality ================================================
        DO IB=1,NAUG
          DO JB=1,NAUG
            AUX(:)=PRO(:,IB,L+1)*PSPHI(:,JB,L+1)*R(:)**2
            CALL RADIAL$INTEGRAl(GID,NR,AUX,svar)
            PROPSIBAR(IB,JB,L+1)=SVAR
          ENDDO
        ENDDO
        WRITE(NFILO,FMT='("<PRO|PSPSI> AFTER BIORTHOGONALIZATION")')
        WRITE(NFILO,FMT='("   (SHOULD BE <PRO_I|PSPSI_J>=DELTA_I,J)")')
        DO IB=1,NAUG
          WRITE(NFILO,FMT='(20F20.10)')PROPSIBAR(IB,:NAUG,L+1)
        ENDDO
!
!       ========================================================================
        PROPSIBAR(:,:,L+1)=0.D0
        DO IB=1,NAUG
          DO JB=1,NAUG
            AUX(:)=PROSAVE(:,IB,L+1)*PSPHISAVE(:,JB,L+1)*R(:)**2
            CALL RADIAL$INTEGRAl(GID,NR,AUX,svar)
            PROPSIBAR(IB,JB,L+1)=SVAR
          ENDDO
        ENDDO

        ALLOCATE(MAT1(NAUG,NAUG))
        ALLOCATE(MAT2(NAUG,NAUG))
        MAT1=MATMUL(MATMUL(TRANSPRO1,MATMUL(DH(:NAUG,:NAUG,L+1),TRANSPOSE(TRANSPRO1))),PROPSIBAR(:NAUG,:NAUG,L+1))
        MAT2=MATMUL(MATMUL(TRANSPRO1,MATMUL(Do(:NAUG,:NAUG,L+1),TRANSPOSE(TRANSPRO1))),PROPSIBAR(:NAUG,:NAUG,L+1))
        DO IB=1,NAUG
          MAT1(IB,IB)=MAT1(IB,IB)+1.D0
          MAT2(:,IB)=MAT2(:,IB)*EGRID(IB,L+1)
        ENDDO
        MAT1=MAT1-MAT2
        write(nfilo,*)'TEST:  MUST VANISH'
        DO IB=1,NAUG
          WRITE(NFILO,FMT='(20F20.10)')MAT1(:,IB)
        ENDDO
        DEALLOCATE(MAT1)
        DEALLOCATE(MAT2)
        deallocate(transphi1)
        deallocate(transpro1)
!
!       ========================================================================
!       ==  CHECK IF THE BARE PSEUDO PARTIAL WAVES FULFILL PAW DGL            ==
!       ========================================================================
        WRITE(NFILO,*)'CHECK IF THE BARE PSEUDO PARTIAL WAVES FULFILL PAW DGL'
        ALLOCATE(PHI2(NR,NAUG))
        ALLOCATE(SVEC(NAUG))
        DO IB=1,NAUG
          PHI2(:,IB)=TPSPHISAVE(:,IB,L+1)+(PSPOT(:)*Y0-EGRID(IB,L+1))*PSPHISAVE(:,IB,L+1)
          DO JB=1,NAUG
            CALL RADIAL$INTEGRAL(GID,NR,R(:)**2*PRO(:,JB,L+1)*PSPHISAVE(:,IB,L+1),SVEC(JB))
          ENDDO
          WRITE(NFILO,fmt='(A,20f10.5)')'<P|PSI>: ',SVEC
          SVEC(:)=MATMUL(DH(:NAUG,:NAUG,L+1)-EGRID(IB,L+1)*DO(:NAUG,:NAUG,L+1),SVEC)
          PHI2(:,IB)=PHI2(:,IB)+MATMUL(PRO(:,:NAUG,L+1),SVEC)
          WRITE(NFILO,fmt='(a,f20.10)')'MAX: (HTILDE -E*OTILDE)|PS-PSI> ',MAXVAL(ABS(PHI2(:,IB)))
        ENDDO
        DEALLOCATE(PHI2)
        DEALLOCATE(SVEC)

WRITE(*,FMT='("1C-KINETIC ENERGY")')
DO IB=1,NAUG
  WRITE(*,FMT='(20E20.5)')DT(IB,:naug,L+1)
ENDDO
WRITE(*,FMT='("1C-OVERLAP")')
DO IB=1,NAUG
  WRITE(*,FMT='(20E20.5)')DO(IB,:naug,L+1)
ENDDO
PRINT*,'WRITE FILES'
        ALLOCATE(WORK2D(NR,4*NAUG))
        WORK2D(:,1:NAUG)=AEPHI(:,:naug,L+1)
        WORK2D(:,NAUG+1:2*NAUG)=PSPHI(:,:naug,L+1)
        WORK2D(:,2*NAUG+1:3*NAUG)=PRO(:,:naug,L+1)
        WORK2D(:,3*NAUG+1:4*NAUG)=UPHI(:,:naug,L+1)
        CALL WRITEPHI('AUGMENT'//TRIM(EXT),GID,NR,4*NAUG,WORK2D)
        DEALLOCATE(WORK2D)
        CALL WRITEPHI('AEPHI'//TRIM(EXT),GID,NR,NAUG,AEPHI(:,:naug,L+1))
        CALL WRITEPHI('AEPHI_T'//TRIM(EXT),GID,NR,NAUG,TAEPHI(:,:naug,L+1))
        CALL WRITEPHI('PSPHI'//TRIM(EXT),GID,NR,NAUG,PSPHI(:,:naug,L+1))
        CALL WRITEPHI('PSPHI_T'//TRIM(EXT),GID,NR,NAUG,TPSPHI(:,:naug,L+1))
        CALL WRITEPHI('PRO'//TRIM(EXT),GID,NR,NAUG,PRO(:,:naug,L+1))
!
!       ========================================================================
!       ==  TEST SCATTERING PROPERTIES                                        ==
!       ========================================================================
!!$PRINT*,'TEST1'
!!$DO IB=1,NAUG
!!$!   -- <P|PHIBAR>=<P|PHI>/TPHI=1/TPHI
!!$  DO JB=1,NAUG
!!$    AUX=PRO(:,JB,L+1)*PSPHISAVE(:,IB,L+1)*R(:)**2
!!$    CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$!    PRINT*,'SVAR ',IB,JB,SVAR,TRANSPHIINV(JB,IB,L+1)
!!$  ENDDO 
!!$  PRINT*,'SVEC        ',IB,MATMUL(TRANSPRO(:,:,L+1),MATMUL(DH(:,:,L+1)-EGRID(IB,L+1)*DO(:,:,L+1),TRANSPHIINV(:,IB,L+1)))
!!$ENDDO

PRINT*,'TEST2'
       IF(L.EQ.0) THEN
         E=-9.D0
         AUX(:)=0.D0
         ALLOCATE(PHI2(NR,1))
         CALL PAWDER(GID,NR,L,E,PSPOT,NAUG,PRO(:,:naug,L+1),DH(:naug,:naug,L+1) &
      &           ,DO(:naug,:naug,L+1),AUX,PHI2(:,1))
         CALL WRITEPHI('XXX.DAT',GID,NR,1,PHI2)
         DEALLOCATE(PHI2)
       END IF
!      STOP 'IN TEST OK'
!
!       ========================================================================
!       ==  TEST SCATTERING PROPERTIES                                        ==
!       ========================================================================
        CALL TESTPHASESHIFT(GID,NR,L,AEPOT,uc,PSPOT,NAUG,PRO(:,:naug,L+1) &
       &                        ,DH(:naug,:naug,L+1),DO(:naug,:naug,L+1))
      ENDDO
PRINT*,'================LOOP DONE========================'
!
!     ==========================================================================
!     ==  UNSCREENING                                                         ==
!     ==========================================================================
print*,'unscreening....'
      CALL UNSCREEN(GID,NR,AEZ,AERHOC+AERHOV,PSRHOC+PSRHOV,PSPOT,RCSM,VADD)
!     == TEST
      ALLOCATE(WORK2D(NR,3))
      WORK2D(:,1)=AEPOT*Y0
      WORK2D(:,2)=PSPOT*Y0
      WORK2D(:,3)=(PSPOT-VADD)*Y0
      CALL WRITEPHI('POTENTIALS.DAT',GID,NR,3,WORK2D)
      DEALLOCATE(WORK2D)
      CALL WRITEPHI('PSRHOV.DAT',GID,NR,1,PSRHOV)
      CALL WRITECOMPARE('AERHO-PSRHO.DAT',GID,NR,1,AERHOV+AERHOC,PSRHOV+PSRHOC)
!
!     ==========================================================================
!     ==  CALCULATE TOTAL ENERGY                                              ==
!     ==========================================================================
print*,'calculating total energy....'
      CALL PSETOTSHELL(GID,NR,AEZ,rbox,PSRHOC,AERHOC,PSPOT &
     &                ,NB,NC,LMAX,NAUGX,NPRO,LOFI,FOFI,EOFI &
     &                ,PRO,AEPHI,PSPHI,RCSM,DT,DH,DO,VADD)
!
!     ==========================================================================
!     ==  WRITE SETUP TO FILE                                                 ==
!     ==========================================================================
print*,'writing setup.....'
      CALL WRITESETUP(GID,NR,AEZ,NB,NC,LOFI,NNOFI,FOFI,EOFI,UOFI,TUOFI &
     &               ,AEPOT,PSPOT,PSRHOC,AERHOC,VADD,RCSM &
     &               ,NPRO,NAUGX,LMAX,AEPHI,PSPHI,PRO,UPHI,TUPHI,DT,DO,DH)
!
!     ==========================================================================
!     ==  READ SETUP FILE AGAIN AS CONSISTENCY CHECK                          ==
!     ==========================================================================
print*,'testing setupfile.....'
      CALL TESTSETUPREAD()
print*,'done....'
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATTACHFILES()
!     **************************************************************************
!     **                                                                      **
!     **  READS THE NAME OF THE CONTROL FILE FROM THE ARGUMENT LIST AND       **
!     **  DEFINES THE FILES IN THE FILEHANDLER OBJECT                         **
!     **                                                                      **
!     **  THE FILES DEFINED ARE                                               **
!     **    PROT     PROTOCOLL FILE                                           **
!     **    INPUT    CONTROL FILE X.ACNTL                                     **
!     **    STP      SETU OUTPUT FILE                                         **
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(128)         :: CNTLNAME
      CHARACTER(128)         :: ROOTNAME
      INTEGER(4)             :: I
!     **************************************************************************

!     ==========================================================================
!     ==  COLLECT INPUT ARGUMENT AND DEFINES ROOTNAME                         ==
!     ==========================================================================
      CALL GETARG(1,CNTLNAME)
      I=INDEX(CNTLNAME,'.',BACK=.TRUE.)
      ROOTNAME=CNTLNAME(1:I-1)
      CALL FILEHANDLER$SETROOT(ROOTNAME)
!
!     ==========================================================================
!     ==  ATTACH PROTOCOLL FILE  X.PROT                                       ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('PROT',.FALSE.,TRIM(ROOTNAME)//-'.PROT')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','FORM','FORMATTED')
!
!     ==========================================================================
!     == DEFINE AND READ INPUT FILE X.ACNTL                                   ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('INPUT',.FALSE.,TRIM(ROOTNAME)//-'.ACNTL')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','FORM','FORMATTED')
!
!     ==========================================================================
!     == DEFINE SETUP FILE                                                    ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('SETUP',.FALSE.,TRIM(ROOTNAME)//-'.STP')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','FORM','FORMATTED')
!
!     ==========================================================================
!     == ATTACH SETUP FILE AS INPUT FILE FOR TESTING                          ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('SETUP_IN',.FALSE.,TRIM(ROOTNAME)//-'.STP')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP_IN','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP_IN','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP_IN','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP_IN','FORM','FORMATTED')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL_DIMENSIONS(LL_CNTL,NB,LMAX,NAUGX)
!     **************************************************************************
!     ** EXTRACT DIMENSION FROM INPUT CONTROL FILE                            **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_CNTL
      INTEGER(4)   ,INTENT(OUT)   :: NB
      INTEGER(4)   ,INTENT(OUT)   :: LMAX
      INTEGER(4)   ,INTENT(OUT)   :: NAUGX
      LOGICAL(4)                  :: TCHK
      INTEGER(4)  ,ALLOCATABLE    :: NPRO(:)
!     **************************************************************************
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$NLISTS(LL_CNTL,'STATE',NB)
!
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'AUGMENTATION')
      CALL LINKEDLIST$EXISTD(LL_CNTL,'NPRO',1,TCHK)
      LMAX=2
      IF(TCHK) THEN
        CALL LINKEDLIST$SIZE(LL_CNTL,'NPRO',1,LMAX)
        LMAX=LMAX-1
      END IF
      ALLOCATE(NPRO(LMAX+1))
      NPRO=2
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'NPRO',1,NPRO)
      NAUGX=MAXVAL(NPRO)
      DEALLOCATE(NPRO)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL1(LL_CNTL,NB,AEZ,RCOV,RBOX,rcut,GID,NR,IRCUT,LOFI,SOFI &
     &                   ,NNOFI,FOFI,TCORE)
!     **************************************************************************
!     ** READ INPUT FILE FOR ATOMIC CALCULATION                               **
!     **                                                                      **
!     **   RMAX: MAXIMUM GRID POINT                                           **
!     **   RBOX: ATOM IS CONFINED IN A BOX WITH RADIUS RBOX                   **
!     **   RCUT: PARTIAL WAVES ARE SET TO ZERO BEYOND RCUT                    **
!     **                                                                      **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_CNTL
      INTEGER(4)   ,INTENT(IN)    :: NB
      REAL(8)     ,INTENT(OUT)    :: AEZ
      INTEGER(4)  ,INTENT(OUT)    :: IRCUT
      REAL(8)     ,INTENT(OUT)    :: RBOX     ! ATOM IS ENCLOSED IN A HARD BOX
      REAL(8)     ,intent(out)    :: RCUT
      INTEGER(4)  ,INTENT(OUT)    :: GID
      INTEGER(4)  ,INTENT(OUT)    :: NR
      INTEGER(4)  ,INTENT(OUT)    :: LOFI(NB)
      INTEGER(4)  ,INTENT(OUT)    :: SOFI(NB)
      INTEGER(4)  ,INTENT(OUT)    :: NNOFI(NB)
      REAL(8)    ,INTENT(OUT)     :: FOFI(NB)
      LOGICAL(4) ,INTENT(OUT)     :: TCORE(NB)
      CHARACTER(2)                :: ELEMENTSYMBOL
      REAL(8)                     :: RCOV
      LOGICAL(4)                  :: TCHK
      REAL(8)                     :: DMIN,DMAX,RX
      INTEGER(4)                  :: NFILO
      INTEGER(4)                  :: IDFT
      INTEGER(4)                  :: NC
      INTEGER(4)                  :: IB,IR
      REAL(8)                     :: R1,DEX
      REAL(8)     ,ALLOCATABLE    :: R(:)
!     **************************************************************************
      CALL FILEHANDLER$UNIT('PROT',NFILO)
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
      CALL REPORT$CHVAL(NFILO,'ELEMENT',ELEMENTSYMBOL)
      CALL REPORT$R8VAL(NFILO,'Z',AEZ,'')
!
!     ==========================================================================
!     == DEFINE RADIAL GRID                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'GRID')
      DMIN=5.D-6
      DMAX=5.D-1
      RX=1.056D-4*EXP(5.D-2*249.D0)
      CALL LINKEDLIST$GET(LL_CNTL,'DMIN',1,DMIN)
      CALL LINKEDLIST$GET(LL_CNTL,'DMAX',1,DMAX)
      CALL LINKEDLIST$GET(LL_CNTL,'RMAX',1,RX)
!
!     == DEFINE GRID ===========================================================
      CALL GRIDPARAMETERS(DMIN,DMAX,RX,R1,DEX,NR)
      CALL RADIAL$NEW('SHLOG',GID)
      CALL RADIAL$SETR8(GID,'R1',R1)
      CALL RADIAL$SETR8(GID,'DEX',DEX)
      CALL RADIAL$SETI4(GID,'NR',NR)
!
!     == REPORT GRID ===========================================================
      CALL REPORT$TITLE(NFILO,'RADIAL GRID')
      CALL REPORT$R8VAL(NFILO,'R1',R1,'A0')
      CALL REPORT$R8VAL(NFILO,'DEX',DEX,'')
      CALL REPORT$I4VAL(NFILO,'#(RADIAL GRID POINTS)',NR,'')
      CALL REPORT$R8VAL(NFILO,'INNERMOST SPACING',DMIN,'A0')
      CALL REPORT$R8VAL(NFILO,'OUTERMOST SPACING',DMAX,'A0')
      CALL REPORT$R8VAL(NFILO,'OUTERMOST GRID POINT',RX,'A0')
      CALL REPORT$R8VAL(NFILO,'COVALENT RADIUS',RCOV,'A0')
!
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == atom is inclosed in a box with radius rbox                           ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RBOX/RCOV',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'RBOX/RCOV',1,RBOX)
        RBOX=RBOX*RCOV   ! PROJECTORS AND PARTIAL WAVES TRUNCATED BEYOND THIS POINT
      ELSE
        CALL LINKEDLIST$EXISTD(LL_CNTL,'RBOX',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_CNTL,'RBOX',1,RBOX)
        ELSE
          RBOX=R(NR)
        END IF
      END IF
      CALL REPORT$R8VAL(NFILO,'BOX RADIUS',RBOX,'A0')
!
!     ==========================================================================
!     == DEFINE RADIUS UP TO WHICH THE DATA ARE STORED                        ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RCUT/RCOV',1,TCHK)
      IF(TCHK) THEN
        CALL LINKEDLIST$GET(LL_CNTL,'RCUT/RCOV',1,RCUT)
        RCUT=RCUT*RCOV   ! PROJECTORS AND PARTIAL WAVES TRUNCATED BEYOND THIS POINT
      ELSE
        RCUT=RBOX
      end if
      CALL REPORT$R8VAL(NFILO,'TRUNCATION RADIUS',RCUT,'A0')
      DO IR=1,NR
        IRCUT=IR
        IF(R(IR).GT.RCUT) EXIT
      ENDDO
!
!     ==========================================================================
!     == DEFINE DFT FUNCTIONAL                                                ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'GENERIC')
      CALL LINKEDLIST$GET(LL_CNTL,'DFT',1,IDFT)
!
      CALL DFT$SETI4('TYPE',IDFT)
      CALL DFT$SETL4('SPIN',.FALSE.)
!
!     == REPORT ================================================================
      CALL DFT$REPORT(NFILO)
!
!     ==========================================================================
!     == DEFINE STATES AND THEIR OCCUPATIONS                                  ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
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
      DO IB=1,Nb
        IF((IB.LE.NC.AND..NOT.TCORE(IB)).OR.(IB.GT.NC.AND.TCORE(IB))) THEN
          CALL ERROR$MSG('CORE STATES MUST PRECEED VALENCE STATES')
          CALL ERROR$I4VAL('NC',NC)
          CALL ERROR$I4VAL('IB',IB)
          CALL ERROR$L4VAL('TCORE',TCORE(IB))
          CALL ERROR$STOP('READCNTL1')
        END IF
      ENDDO
      
      CALL REPORT$TITLE(NFILO,'INPUT ON STATES')
      DO IB=1,NB
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE READCNTL2(LL_CNTL,LMAX,RCOV,RC_CORE,POW_CORE,VAL0_CORE &
     &                    ,RC_POT,POW_POT,VAL0_POT,RCSM,PSEUDOPARTIALWAVEMETHOD &
     &                    ,NPRO,RCL)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_CNTL
      INTEGER(4)   ,INTENT(IN)    :: LMAX
      REAL(8)      ,INTENT(IN)    :: RCOV
      REAL(8)      ,INTENT(OUT)   :: RC_CORE,VAL0_CORE
      INTEGER(4)   ,INTENT(OUT)   :: POW_CORE
      REAL(8)      ,INTENT(OUT)   :: RC_POT,VAL0_POT
      INTEGER(4)   ,INTENT(OUT)   :: POW_POT
      REAL(8)      ,INTENT(OUT)   :: RCSM
      INTEGER(4)   ,INTENT(OUT)   :: NPRO(LMAX+1)
      REAL(8)      ,INTENT(OUT)   :: RCL(LMAX+1)
      CHARACTER(*) ,INTENT(OUT)   :: PSEUDOPARTIALWAVEMETHOD
      REAL(8)                     :: PI,Y0
      REAL(8)                     :: SVAR
      INTEGER(4)                  :: NFILO
      INTEGER(4)                  :: L
      LOGICAL(4)                  :: TCHK
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
!
!     ==========================================================================
!     == COLLECT PARAMETERS FOR CORE PSUEDIZATION                             ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'PSCORE')
      CALL LINKEDLIST$GET(LL_CNTL,'RC/RCOV',1,RC_CORE)!CUTOFF RADIUS FOR CORE PSEUDIZATION
      RC_CORE=RC_CORE*RCOV
      CALL LINKEDLIST$GET(LL_CNTL,'POWER',1,POW_CORE)  !POWER FOR CORE PSEUDIZATION
      CALL LINKEDLIST$GET(LL_CNTL,'VAL0',1,VAL0_CORE)  !VALUE AT THE ORIGIN FOR CORE PSEUDIZATION
      VAL0_CORE=VAL0_CORE/Y0
!
!     == REPORT ================================================================
      CALL REPORT$TITLE(NFILO,'PSEUDO CORE DENSITY')
      CALL REPORT$R8VAL(NFILO,'RC',RC_CORE,'')
      CALL REPORT$I4VAL(NFILO,'POWER',POW_CORE,'')
      CALL REPORT$R8VAL(NFILO,'VAL(0)',VAL0_CORE*Y0,'')
!
!     ==========================================================================
!     == COLLECT PARAMETERS FOR POTENTIAL PSEUDIZATION                        ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'VTILDE')
      CALL LINKEDLIST$GET(LL_CNTL,'RC/RCOV',1,RC_POT)!CUTOFF RADIUS FOR POTENTIAL PSEUDIZATION
      RC_POT=RC_POT*RCOV
      CALL LINKEDLIST$GET(LL_CNTL,'POWER',1,POW_POT)  !POWER FOR POTENTIAL PSEUDIZATION 
      CALL LINKEDLIST$GET(LL_CNTL,'VAL0',1,VAL0_POT)  !VALUE AT THE ORIGIN FOR POTENTIAL PSEUDIZATION
      VAL0_POT=VAL0_POT/Y0
!
!     == REPORT ================================================================
      CALL REPORT$TITLE(NFILO,'PSEUDO POTENTIAL')
      CALL REPORT$R8VAL(NFILO,'RC',RC_POT,'')
      CALL REPORT$I4VAL(NFILO,'POWER',POW_POT,'')
      CALL REPORT$R8VAL(NFILO,'VAL(0)',VAL0_POT*Y0,'')
!
!     ==========================================================================
!     ==  COMPENSATION GAUSSIAN                                               ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'COMPENSATE')
      CALL LINKEDLIST$GET(LL_CNTL,'RC/RCOV',1,RCSM)  
      RCSM=RCSM*RCOV
!
!     == REPORT ================================================================
      CALL REPORT$TITLE(NFILO,'COMPENSATION DENSITY')
      CALL REPORT$R8VAL(NFILO,'RC',RCSM,'')
      SVAR=EXP(-(RCOV/RCSM)**2)/(SQRT(PI)*RCSM)**3
      CALL REPORT$R8VAL(NFILO,'VALUE AT THE COVALENT RADIUS',SVAR,'')
!
!     ==========================================================================
!     ==  PSEUDIZATION ALGORITHM FOR PSEUDO PARTIAL WAVES                     ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_CNTL,'~')
      CALL LINKEDLIST$SELECT(LL_CNTL,'ACNTL')  
      CALL LINKEDLIST$SELECT(LL_CNTL,'AUGMENTATION')
      PSEUDOPARTIALWAVEMETHOD='HBS'
      CALL LINKEDLIST$EXISTD(LL_CNTL,'PSEUDOMETHOD',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'PSEUDOMETHOD' &
     &                                   ,1,PSEUDOPARTIALWAVEMETHOD)
      IF(PSEUDOPARTIALWAVEMETHOD.NE.'HBS'.AND. &
     &   PSEUDOPARTIALWAVEMETHOD.NE.'KRESSE') THEN
        CALL ERROR$MSG('ILLEGAL VALUE FOR PSEUDIZATION METHOD')
        CALL ERROR$CHVAL('PSEUDOPARTIALWAVEMETHOD',PSEUDOPARTIALWAVEMETHOD)
        CALL ERROR$STOP('MAIN')
      END IF
!
!     ==========================================================================
!     ==  NORB                                                                ==
!     ==========================================================================
      CALL REPORT$TITLE(NFILO,'AUGMENTATION')
      CALL REPORT$I4VAL(NFILO,'LMAX',LMAX,'')
!
      NPRO=2
      CALL LINKEDLIST$existd(LL_CNTL,'NPRO',1,tchk)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'NPRO',1,NPRO)
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
!
!     ==========================================================================
!     ==  RCL                                                                 ==
!     ==========================================================================
      RCL(:)=RCOV
      CALL LINKEDLIST$EXISTD(LL_CNTL,'RC/RCOV',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_CNTL,'RC/RCOV',1,RCL)
      RCL(:)=RCL(:)*RCOV
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
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITESETUP(GID,NR,AEZ,NB,NC,LOFI,NNOFI,FOFI,EOFI,UOFI,TUOFI &
     &                     ,AEPOT,PSPOT,PSRHOC,AERHOC,VADD,RCSM &
     &                     ,NPRO,NAUGX,LMAX,AEPHI,PSPHI,PRO,UPHI,TUPHI,DT,DO,DH)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)    :: GID
      INTEGER(4)   ,INTENT(IN)    :: NR
      REAL(8)      ,INTENT(IN)    :: AEZ
      REAL(8)                     :: PSZ
      INTEGER(4)   ,INTENT(IN)    :: NB
      INTEGER(4)   ,INTENT(IN)    :: NC
      INTEGER(4)   ,INTENT(IN)    :: LOFI(NB)
      INTEGER(4)   ,INTENT(IN)    :: NNOFI(NB)
      REAL(8)      ,INTENT(IN)    :: FOFI(NB)
      REAL(8)      ,INTENT(IN)    :: EOFI(NB)
      REAL(8)      ,INTENT(IN)    :: UOFI(NR,NB)
      REAL(8)      ,INTENT(IN)    :: TUOFI(NR,NB)
      REAL(8)      ,INTENT(IN)    :: AEPOT(NR)
      REAL(8)      ,INTENT(IN)    :: PSPOT(NR)
      REAL(8)      ,INTENT(IN)    :: PSRHOC(NR)
      REAL(8)      ,INTENT(IN)    :: AERHOC(NR)
      REAL(8)      ,INTENT(IN)    :: VADD(NR)
      REAL(8)      ,INTENT(IN)    :: RCSM
      INTEGER(4)   ,INTENT(IN)    :: NPRO(LMAX+1)
      INTEGER(4)   ,INTENT(IN)    :: NAUGX
      INTEGER(4)   ,INTENT(IN)    :: LMAX
      REAL(8)      ,INTENT(IN)    :: AEPHI(NR,NAUGX,LMAX+1)
      REAL(8)      ,INTENT(IN)    :: PSPHI(NR,NAUGX,LMAX+1)
      REAL(8)      ,INTENT(IN)    :: PRO(NR,NAUGX,LMAX+1)
      REAL(8)      ,INTENT(IN)    :: UPHI(NR,NAUGX,LMAX+1)
      REAL(8)      ,INTENT(IN)    :: TUPHI(NR,NAUGX,LMAX+1)
      REAL(8)      ,INTENT(IN)    :: DT(NAUGX,NAUGX,LMAX+1)
      REAL(8)      ,INTENT(IN)    :: DO(NAUGX,NAUGX,LMAX+1)
      REAL(8)      ,INTENT(IN)    :: DH(NAUGX,NAUGX,LMAX+1)
      TYPE(LL_TYPE)               :: LL_STP
      REAL(8)                     :: R1
      REAL(8)                     :: DEX
      LOGICAL(4)  ,PARAMETER      :: TLOG=.false.
      INTEGER(4)                  :: NR2
      INTEGER(4)                  :: GID2
      INTEGER(4)                  :: IDFT
      INTEGER(4)                  :: LNX
      INTEGER(4)                  :: L,IB,IB1,IB2,L1,L2,I1,I2,ir
      INTEGER(4)                  :: ISVAR
      INTEGER(4) ,ALLOCATABLE     :: LOX(:)
      REAL(8)    ,ALLOCATABLE     :: AUX(:)
      REAL(8)    ,ALLOCATABLE     :: WORK1(:),WORK2(:),WORK3(:)
      INTEGER(4)                  :: NFIL
      real(8)                     :: pi,y0,val
      REAL(8)    ,ALLOCATABLE     :: r(:),auxa(:)
      real(8)                     :: rcut
!     **************************************************************************
      pi=4.d0*atan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      allocate(r(nr))
      call radial$r(gid,nr,r)
      rcut=r(nr)
      deallocate(r)

      CALL LINKEDLIST$NEW(LL_STP)
!
!     == GRID  =================================================================
      CALL RADIAL$GETR8(GID,'R1',R1)
      CALL RADIAL$GETR8(GID,'DEX',DEX)
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      CALL LINKEDLIST$SELECT(LL_STP,'GRID')
      IF(TLOG) THEN
        CALL LINKEDLIST$SET(LL_STP,'TYPE',0,'SHLOG')
        CALL LINKEDLIST$SET(LL_STP,'R1',0,R1)
        CALL LINKEDLIST$SET(LL_STP,'DEX',0,DEX)
        CALL LINKEDLIST$SET(LL_STP,'NR',0,NR)
        ALLOCATE(AUX(NR))
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
        ALLOCATE(AUX(NR2))
      END IF
!
!     == GENERIC ===============================================================
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$SELECT(LL_STP,'SETUP')
      CALL LINKEDLIST$SELECT(LL_STP,'GENERIC')
      CALL DFT$GETI4('TYPE',IDFT)
      CALL LINKEDLIST$SET(LL_STP,'IDFT',0,IDFT)
      CALL LINKEDLIST$SET(LL_STP,'AEZ',0,AEZ)
!
!     == AUGMENTATION ==========================================================
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
!!$PRINT*,'AEPOT ',AEPOT(1:NR:10)
!!$ALLOCATE(G(NR))
!!$SVAR=-18.D0
!!$G(:)=0.D0
!!$CALL BOUNDSTATE(GID,NR,0,0,DREL,G,0,AEPOT,SVAR,AEPHI(:,1,1))
!!$PRINT*,'SVAR',SVAR
!!$STOP
        DO IB=1,NB
          CALL LINKEDLIST$SET(LL_STP,'PHI-NODELESS',ib,UOFI(:,IB))
          CALL LINKEDLIST$SET(LL_STP,'TPHI-NODELESS',ib,TUOFI(:,IB))
        ENDDO
      ELSE
        CALL RADIAL$CHANGEGRID(GID,NR,AEPOT,GID2,NR2,AUX)
        CALL LINKEDLIST$SET(LL_STP,'AEPOT',0,AUX)
        DO IB=1,NB
          CALL RADIAL$CHANGEGRID(GID,NR,UOFI(:,IB),GID2,NR2,AUX)
          CALL LINKEDLIST$SET(LL_STP,'PHI-NODELESS',-1,AUX)
          CALL RADIAL$CHANGEGRID(GID,NR,TUOFI(:,IB),GID2,NR2,AUX)
          CALL LINKEDLIST$SET(LL_STP,'TPHI-NODELESS',-1,AUX)
        ENDDO
      END IF
!
!     == AUGMENTATION ==========================================================
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
!     == SET MATRIX ELEMENTS ===================================================
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
      DEALLOCATE(AUX)
!
!     ==========================================================================
!     ==  WRITE SETUP TO FILE                                                 ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('SETUP',NFIL)
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$WRITE(LL_STP,NFIL,'~')
!
!     ==  CHECK CONTENTS =======================================================
!     CALL LINKEDLIST$SELECT(LL_STP,'~')
!     CALL LINKEDLIST$REPORT(LL_STP,NFILO)
!
      CALL FILEHANDLER$CLOSE('SETUP')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TESTSETUPREAD()
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)              :: NFIL
      INTEGER(4)              :: NFILO
      LOGICAL(4)              :: TCHK
      INTEGER(4)              :: LNX
      INTEGER(4)              :: NC
      INTEGER(4)              :: GID
      REAL(8)                 :: R1
      REAL(8)                 :: DEX
      INTEGER(4)              :: NR       ! #(RDIAL GRID POINTS)
      REAL(8)   ,ALLOCATABLE  :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE  :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE  :: PRO(:,:)
      INTEGER(4),ALLOCATABLE  :: LOX(:)
      REAL(8)                 :: AEZ,PSZ
      REAL(8)                 :: RCSM
      REAL(8)   ,ALLOCATABLE  :: PSCORE(:)
      REAL(8)   ,ALLOCATABLE  :: AECORE(:)
      REAL(8)   ,ALLOCATABLE  :: AEPOT(:)
      REAL(8)   ,ALLOCATABLE  :: DTKIN(:,:)
      REAL(8)   ,ALLOCATABLE  :: DO(:,:)
      REAL(8)   ,ALLOCATABLE  :: VHAT(:)
      REAL(8)   ,ALLOCATABLE  :: AEPSICORE(:,:)
      INTEGER(4),ALLOCATABLE  :: LB(:)
      REAL(8)   ,ALLOCATABLE  :: FB(:)
      REAL(8)   ,ALLOCATABLE  :: EB(:)
!     **************************************************************************
      CALL FILEHANDLER$UNIT('SETUP_IN',NFIL)
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
      ALLOCATE(DTKIN(LNX,LNX))
      ALLOCATE(DO(LNX,LNX))
      CALL SETUPREAD$GETR8A('DT',LNX*LNX,DTKIN)
      CALL SETUPREAD$GETR8A('DO',LNX*LNX,DO)
      ALLOCATE(AEPOT(NR))
      CALL SETUPREAD$GETR8A('AEPOT',NR,AEPOT)
      ALLOCATE(VHAT(NR))
      CALL SETUPREAD$GETR8A('VADD',NR,VHAT)
      CALL SETUPREAD$GETI4('NC',NC)
      ALLOCATE(LB(NC))
      ALLOCATE(FB(NC))
      ALLOCATE(EB(NC))
      CALL SETUPREAD$GETI4A('LOFC',NC,LB)
      CALL SETUPREAD$GETR8A('FOFC',NC,FB)
      CALL SETUPREAD$GETR8A('EOFC',NC,EB)
      ALLOCATE(AEPSICORE(NR,NC))
      CALL SETUPREAD$GETR8A('AEPSICORE',NR*NC,AEPSICORE)
!
!     === REPORT INFORMATION ON THE SETP FILE ==============================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL RADIAL$GETR8(GID,'R1',R1)
      CALL RADIAL$GETR8(GID,'DEX',DEX)
!
! WRITE TRADITIONAL OUT FILE ===============================================
      CALL FILEHANDLER$SETFILE('SETUPO',.TRUE.,-'.OUT')
      CALL FILEHANDLER$SETSPECIFICATION('SETUPO','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('SETUPO','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('SETUPO','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('SETUPO','FORM','FORMATTED')
      CALL FILEHANDLER$UNIT('SETUPO',NFIL)
      CALL WRITEOLDSETUP(NFIL,R1,DEX,NR,LNX,LOX,NC,PSZ,AEZ,RCSM &
     &              ,VHAT,AECORE,PSCORE,DTKIN,DO,AEPHI,PSPHI,PRO,AEPOT &
     &              ,LB,FB,EB,AEPSICORE)
      CALL FILEHANDLER$CLOSE('SETUPO')
      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEOLDSETUP(NFIL,R1,DEX,NR,LNX,LOX,NC,PSZ,AEZ,RCSMALL &
     &              ,VHAT,AECORE,PSCORE,DTKIN,DO,AEPHI,PSPHI,PRO,AEPOT &
     &              ,LB,FB,EB,AEPSI)
      IMPLICIT NONE      
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX) 
      INTEGER(4),INTENT(IN) :: NC
      REAL(8)   ,INTENT(IN) :: PSZ
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: RCSMALL
      REAL(8)   ,INTENT(IN) :: VHAT(NR)
      REAL(8)   ,INTENT(IN) :: AECORE(NR)
      REAL(8)   ,INTENT(IN) :: PSCORE(NR)
      REAL(8)   ,INTENT(IN) :: DTKIN(LNX,LNX)
      REAL(8)   ,INTENT(IN) :: DO(LNX,LNX)
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: PRO(NR,LNX)
      REAL(8)   ,INTENT(IN) :: AEPOT(NR)
      INTEGER(4),INTENT(IN) :: LB(NC)
      REAL(8)   ,INTENT(IN) :: FB(NC)
      REAL(8)   ,INTENT(IN) :: EB(NC)
      REAL(8)   ,INTENT(IN) :: AEPSI(NR,NC)
      INTEGER(4)            :: LN,IC
!     **************************************************************************
!     ==========================================================================
      WRITE(NFIL,FMT='(F15.10,F10.5,2I4,2F5.2,F15.12,I5)') &
     &                 R1,DEX,NR,LNX,PSZ,AEZ,RCSMALL,NC
      WRITE(NFIL,FMT='(14I5)')LOX(:)
      WRITE(NFIL,FMT='(SP,5E14.8)')VHAT(:)
!     ====  AECORE = CORE CHARGE DENSITY  ==============================
      WRITE(NFIL,FMT='(SP,5E14.8)')AECORE(:)
!     ====  PSCORE = PSEUDIZED CHARGE DENSITY===========================
      WRITE(NFIL,FMT='(SP,5E14.8)')PSCORE(:)
!     ====  DTKIN = <AEPHI|-DELTA/2|AEPHI> - <PSPHI|-DELTA/2|PSPHI> ====
      WRITE(NFIL,FMT='(SP,5E14.8)')DTKIN(:,:)
!     ====  DOVER = <AEPHI|AEPHI> - <PSPHI|PSPHI> ======================
      WRITE(NFIL,FMT='(SP,5E14.8)')DO(:,:)
      DO LN=1,LNX
        WRITE(NFIL,FMT='(SP,5E14.8)')PRO(:,LN)
        WRITE(NFIL,FMT='(SP,5E14.8)')AEPHI(:,LN)
        WRITE(NFIL,FMT='(SP,5E14.8)')PSPHI(:,LN)
      ENDDO
!     ==================================================================
!     == ALL ELECTRON POTENTIAL FOR 
!     ==================================================================
      WRITE(NFIL,FMT='(SP,5E14.8)')AEPOT(:)

!     ==================================================================
!     == AE CORE PARTIAL WAVES
!     ==================================================================
! SANTOS040616 BEGIN
      WRITE(NFIL,FMT='(14I5)')LB(1:NC)      ! MAIN ANGULAR MOMENTUM
      WRITE(NFIL,FMT='(SP,5E14.8)')FB(1:NC) ! OCCUPATIONS
      WRITE(NFIL,FMT='(SP,5E14.8)')EB(1:NC) ! ONE-PARTICLE ENERGIES
      DO IC=1,NC
        WRITE(NFIL,FMT='(SP,5E14.8)')AEPSI(:,IC)
      ENDDO
! SANTOS040616 END
      
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
      PI=4.D0*ATAN(1.D0)
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
        WRITE(100,FMT='(F15.10,2X,20(F25.15,2X))')R(IR),PHI(IR,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITECOMPARE(FILE,GID,NR,NPHI,PHI1,PHI2)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI1(NR,NPHI)
      REAL(8)     ,INTENT(IN) :: PHI2(NR,NPHI)
      INTEGER(4)              :: IR
      REAL(8)                 :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        WRITE(100,FMT='(F15.10,2X,8(F25.15,2X))')R(IR),PHI1(IR,:),PHI2(IR,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TESTPHASESHIFT(GID,NR,L,AEPOT,UC,PSPOT,NAUG,PRO,DH,DO)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID      ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L        ! ANGULAR MOMENTUM
      REAL(8)    ,INTENT(IN)     :: AEPOT(NR)  ! AE POTENTIAL (RADIAL PART ONLY)
      REAL(8)    ,INTENT(IN)     :: PSPOT(NR)  ! PS POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NAUG       ! 
      REAL(8)    ,INTENT(IN)     :: UC(NR)     ! NODELESS CORE STATE OR ZERO
      REAL(8)    ,INTENT(IN)     :: PRO(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: DH(NAUG,NAUG)
      REAL(8)    ,INTENT(IN)     :: DO(NAUG,NAUG)
      INTEGER(4) ,PARAMETER      :: NE=100
      REAL(8)    ,PARAMETER      :: EMIN=-3.D0,EMAX=2.D0
      REAL(8)    ,PARAMETER      :: RPHASE=3.D0
      REAL(8)                    :: PSPHASE(NE,NAUG)
      REAL(8)                    :: AEPHASE(NE)
      REAL(8)                    :: NLPHASE(NE)
      REAL(8)                    :: PSPHI(NR)
      REAL(8)                    :: AEPHI(NR)
      REAL(8)                    :: G(NR)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: DREL(NR)
      INTEGER(4)                 :: SO
      INTEGER(4)                 :: IE,NPRO,IR
      REAL(8)                    :: E
      REAL(8)                    :: SVAR
!     ***************************************************************************
      CALL RADIAL$R(GID,NR,R)
      DO IE=1,NE
        E=EMIN+(EMAX-EMIN)/REAL(NE-1)*REAL(IE-1)
!
!       ==  PHASESHIFT OF ALL-ELECTRON SOLUTION  ===============================
        DREL(:)=0.D0
        SO=0.D0
        CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,SO,G,L,E,1,AEPHI)
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,AEPHI,RPHASE,AEPHASE(IE))
!
!       ==  PHASESHIFT OF NODELESS WAVE FUNCTION   =============================
        DREL(:)=0.D0
        SO=0.D0
        CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,SO,UC,L,E,1,AEPHI)
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,AEPHI,RPHASE,NLPHASE(IE))
!
!       ==  PHASESHIFT FROM PAW PSEUDO WAVE FUNCTIONS  =========================
        G(:)=0.D0
        DO NPRO=1,NAUG
          CALL PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO(:,:NPRO) &
                    ,DH(:NPRO,:NPRO),DO(:NPRO,:NPRO),G,PSPHI)
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSPHI,RPHASE,PSPHASE(IE,NPRO))
!
!         == CLEAN UP STEPS IN THE PHASE SHIFT
          IF(IE.GT.1) THEN
            SVAR=PSPHASE(IE,NPRO)-PSPHASE(IE-1,NPRO)
            IF(ABS(SVAR).GT.0.9D0) THEN
              PSPHASE(IE,NPRO)=PSPHASE(IE,NPRO)-NINT(SVAR)
            END IF
          ELSE
            SVAR=PSPHASE(1,NPRO)
            IF(ABS(SVAR).GT.0.9D0) THEN
              PSPHASE(IE,NPRO)=PSPHASE(IE,NPRO)-NINT(SVAR)
            END IF
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == WRITE RESULT                                                         ==
!     ==========================================================================
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
        WRITE(100,FMT='(F20.10,2X,F20.10,2X,20F20.10)')E,NLPHASE(IE),AEPHASE(IE),PSPHASE(IE,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CONSTRUCTPSPHI(GID,NR,RC,L,NB,PSPHI,TPSPHI)
!     **                                                                      **
!     **  PSEUDIZES A WAVE FUNCTION BY MATCHING SPHERICAL BESSEL FUNCTIONS    **
!     **                                                                      **
!     **  FIRST WE FIND THE VALUES K_I SO THAT THE SPHERICAL BESSEL           **
!     **  FUNCTION J_L(K_I*R) MATCHES THE LOGARITHMIC DERIVATIVE              **
!     **  OF THE FUNCTION TO BE PSEUDIZED.                                    **
!     **                                                                      **
!     **  THE FIRST TWO ARE THEN MATCHED SO THAT THE PHI AND TPHI             **
!     **  ARE CONTINUOUS. THE DERIVATIVE OF PHI IS AUTOMATICALLY              **
!     **  CONTINUOUS THROUGH THE CORRECT CHOICE OF THE K_I.                   **
!     **                                                                      **
!     **  ATTENTION! SINCE WE MATCH TO TAYLOR COEFFICIENTS, ONE SHOULD        **
!     **  ALSO MATCH ENERGY DERIVATIVES OF THE SPHERICAL BESSEL FUNCTIONS.    **
!     **                                                                      **
!     **  NOTE THAT THE TAYLOR EXPANSION OF THE WAVE FUNCTIONS IN ENERGY      **
!     **  BEHAVE DIFFERENT FROM THE ENERGY DEPENDENT SOLUTION.                **
!     **  THE TAYLOR EXPANSION CAN CREATE TWO NODES AT ONCE, WHICH            **
!     **  CAUSES DIFFICULTIES.                                                **
!     **                                                                      **
!     **                                                                      ** 
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)     :: GID
      INTEGER(4),INTENT(IN)     :: NR
      INTEGER(4),INTENT(IN)     :: L
      INTEGER(4),INTENT(IN)     :: NB
      REAL(8)   ,INTENT(IN)     :: RC
      REAL(8)   ,INTENT(INOUT)  :: PSPHI(NR,NB)
      REAL(8)   ,INTENT(INOUT)  :: TPSPHI(NR,NB)
      INTEGER(4),PARAMETER      :: NZEROS=2
      INTEGER(4),PARAMETER      :: NITER=100
      INTEGER(4),PARAMETER      :: NX=1000
      LOGICAL(4),PARAMETER      :: TDIFFERENTIABLELAPLACIAN=.FALSE.
      REAL(8)   ,PARAMETER      :: XMAX=30.D0
      REAL(8)                   :: R(NR)
      INTEGER(4)                :: IR,I,J,ITER
      REAL(8)                   :: TOL=1.D-12
      REAL(8)                   :: PI
      REAL(8)                   :: VAL,VALT,DERT,VALJ,VALTJ
      REAL(8)                   :: X(NX)
      REAL(8)                   :: JL(NX),DJLDR(NX),D2JLDR2(NX),PHASE(NX)
      REAL(8)                   :: TJL(NX)
      REAL(8)                   :: JLOFKR(NR,NZEROS),TJLOFKR(NR,NZEROS)
      REAL(8)                   :: DJLOFKR(NR,NZEROS),D2JLOFKR(NR,NZEROS),RDUMMY(NR,NZEROS)
      REAL(8)                   :: AUX1(NR),AUX2(NR)
      REAL(8)                   :: NN
      INTEGER(4)                :: ISTART,IBI
      REAL(8)                   :: X0,Y0,DX,XM,YM
      REAL(8)                   :: KI(NZEROS)
      REAL(8)                   :: A11,A12,A21,A22,DET
      REAL(8)                   :: AINV11,AINV12,AINV21,AINV22
      REAL(8)                   :: V1,V2,C1,C2,C3,C1A,C2A,C1B,C2B
      REAL(8)                   :: DTJ1,DTJ2,DTJ3
      REAL(8)                   :: SVAR1,SVAR2,FAC,SVAR
      INTEGER(4)                :: IB
      REAL(8)                   :: SUPPHI(NR,NB)
      REAL(8)                   :: SUPTPHI(NR,NB)
      REAL(8)                   :: SHIFT
      REAL(8)                   :: X1,XDEX
      REAL(8)                   :: OFFSETNN
      INTEGER(4)                :: GIDX
      INTEGER(4)                :: IRMATCH
!     **************************************************************************
!CALL WRITEPHI('TESTCONSTRUCTPSPHI1',GID,NR,NB,PSPHI)
      PI=4.D0*ATAN(1.D0)
      CALL RADIAL$NEW('SHLOG',GIDX)
      XDEX=1.D-5
      X1=XMAX/XDEX/REAL(NX)
      CALL RADIAL$SETR8(GIDX,'R1',X1)
      CALL RADIAL$SETR8(GIDX,'DEX',XDEX)
      CALL RADIAL$SETI4(GIDX,'NR',NX)
      CALL RADIAL$R(GIDX,NX,X)
!
!     ==========================================================================
!     ==  MAKE A COPY OF THE INPUT WAVE FUNCTIONS                             ==
!     ==========================================================================
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$XOFR(GID,RC,VAL)
      IRMATCH=INT(VAL)
      SUPPHI(:,:)=PSPHI(:,:)
      SUPTPHI(:,:)=TPSPHI(:,:)
!
!     =======================================================================
!     == CALCULATE BESSEL FUNCTION                                         ==
!     =======================================================================
      DO IR=1,NX
        CALL SPECIALFUNCTION$BESSEL(L,X(IR),JL(IR)) 
      ENDDO
!
!     ==========================================================================
!     == CALCULATE PHASESHIFT OF BESSEL FUNCTION                              ==
!     ==========================================================================
      CALL RADIAL$VERLETD1(GIDX,NX,JL,DJLDR)
      NN=0.D0
      DO IR=2,NX
        IF(JL(IR)*JL(IR-1).LT.0.D0)NN=NN+1.D0
        PHASE(IR)=0.5D0-ATAN(X(IR)/RC*DJLDR(IR)/JL(IR))/PI+NN
      ENDDO
!     == AVOID GLITCH
      PHASE(1)=PHASE(2)
      PHASE(NX)=PHASE(NX-1)
!
!     ==========================================================================
!     == DETERMINE OFFSET IN THE NUMBER OF NODES                              ==
!     ==========================================================================
      CALL SCHROEDINGER$PHASESHIFT(GID,NR,SUPPHI(:,1),RC,SHIFT)
      OFFSETNN=REAL(INT(SHIFT),KIND=8)
!     
!     ==========================================================================
!     == PSEUDIZE EVERY FUNCTION ON THE SUPPORT ENERGY GRID                   ==
!     ==========================================================================
      DO IB=1,NB
!       ========================================================================
!       == DETERMINE K-VALUES FOR WHICH THE LOGARITHIC DERIVATIVE MATCHES     ==
!       ========================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,SUPPHI(:,IB),RC,SHIFT)
        SHIFT=SHIFT-OFFSETNN
        IF(SHIFT.LT.PHASE(1)) THEN
          CALL ERROR$MSG('MATCHING INCLUDING FIRST BESSEL FUNCTION NOT POSSIBLE')
          CALL ERROR$MSG('TRY TO INCREASE MATCHING RADIUS')
          CALL ERROR$I4VAL('IB',IB)
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('CONSTRUCTPSPHI')
PRINT*,'WARNING SHIFT WAS REQUIRED ',L,IB,PHASE(1),SHIFT
           SHIFT=SHIFT+REAL(INT(PHASE(1)-SHIFT+1.D0))
        END IF
!
!       ========================================================================
!       == SEARCH BESSEL FUNCTIONS JL(K*R) WHICH HAVE IDENTICAL LOGARITHMIC   ==
!       == DERIVATIVES TO THE NODELESS FUNCTION AT THE RADIUS RCS             ==
!       ========================================================================
        DO I=1,NZEROS
          ISTART=1
          X0=0.4D0+1.D-3
          DX=0.1D0
          CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
          DO ITER=1,NITER
            CALL RADIAL$VALUE(GIDX,NX,PHASE,X0,Y0)
            Y0=Y0-SHIFT
            CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
            IF(X0.GT.XMAX) THEN
              CALL ERROR$MSG('RESULTING K OUT OF RANGE')
              CALL ERROR$STOP('CONSTRUCTPSPHI')
            END IF
            IF(ABS(DX).LT.TOL) EXIT
          ENDDO
          IF(ABS(DX).GT.TOL) THEN
            CALL ERROR$MSG('LOOP NOT CONVERGED')
            CALL ERROR$STOP('CONSTRUCTPSPHI')
          END IF
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
!
!       ========================================================================
!       == PLACE BESSELFUNCTIONS ON RADIAL GRID                               ==
!       ========================================================================
        DO I=1,NZEROS
          DO IR=1,NR
            CALL SPECIALFUNCTION$BESSEL(L,KI(I)*R(IR),JLOFKR(IR,I)) 
          ENDDO
          TJLOFKR(:,I)=0.5D0*KI(I)**2*JLOFKR(:,I)
        ENDDO
!
!       ========================================================================
!       == CONSTRUCT NEW PSEUDO WAVE FUNCTION                                 ==
!       ========================================================================
        DO I=1,NZEROS
!         == SELECT AND NORMALIZE THIS AND HIGHER BESSEL FUNCTIONS WITH ========
!         == RESPECT TO PROPERTY ===============================================
          DO J=I,NZEROS
            IF(I.EQ.1) THEN
              CALL RADIAL$VALUE(GID,NR,JLOFKR(:,J),RC,SVAR)
            ELSE IF(I.EQ.2) THEN
              CALL RADIAL$VALUE(GID,NR,TJLOFKR(:,J),RC,SVAR)
            ELSE IF(I.EQ.3) THEN
              CALL RADIAL$DERIVATIVE(GID,NR,TJLOFKR(:,J),RC,SVAR)     
              SVAR=SVAR*1.D+3
            END IF
            IF(ABS(SVAR).LT.1.D-7) THEN
              CALL ERROR$MSG('PROBLEM WITH DIVIDE-BY ZERO')
              CALL ERROR$MSG('CHANGE MATCHING RADIUS FOR PS-PARTIALWAVE CONSTRUCTION')
              CALL ERROR$I4VAL('L',L)
              CALL ERROR$I4VAL('IB',IB)
              CALL ERROR$R8VAL('RC',RC)
              CALL ERROR$I4VAL('I',I)
              CALL ERROR$I4VAL('J',J)
              CALL ERROR$STOP('CONSTRUCTPSPHI')
            END IF
            JLOFKR(:,J) =JLOFKR(:,J)/SVAR
            TJLOFKR(:,J)=TJLOFKR(:,J)/SVAR
          ENDDO
!
!         == PROJECT PROPERTY FROM HIGHER BESSEL FUNCTIONS =====================
          DO J=I+1,NZEROS
            JLOFKR(:,J) =JLOFKR(:,J)-JLOFKR(:,I)
            TJLOFKR(:,J)=TJLOFKR(:,J)-TJLOFKR(:,I)
          ENDDO
!
!         == EXTRACT MISSING PART OF THIS PROPERTY IN THE PARTIAL WAVE =========
          AUX1(:)=SUPPHI(:,IB)
          AUX2(:)=SUPTPHI(:,IB)
          DO J=1,I-1
            AUX1(:)=AUX1(:)-JLOFKR(:,J)
            AUX2(:)=AUX2(:)-TJLOFKR(:,J)
          ENDDO
          IF(I.EQ.1) THEN
            CALL RADIAL$VALUE(GID,NR,AUX1,RC,SVAR)
          ELSE IF(I.EQ.2) THEN
            CALL RADIAL$VALUE(GID,NR,AUX2,RC,SVAR)
          ELSE IF(I.EQ.3) THEN
            CALL RADIAL$DERIVATIVE(GID,NR,AUX2,RC,SVAR)     
            SVAR=SVAR*1.D+3
          END IF
!         == SCALE THIS BESSEL FUNCTION TO MATCH THE MISSING PART ==============
          JLOFKR(:,I) =JLOFKR(:,I)*SVAR
          TJLOFKR(:,I)=TJLOFKR(:,I)*SVAR
        ENDDO
!
!       ========================================================================
!       ==                                                                    ==
!       ========================================================================
        DO I=1,NZEROS
          DO J=1,I-1
            JLOFKR(:,I)=JLOFKR(:,I)+JLOFKR(:,J)
            TJLOFKR(:,I)=TJLOFKR(:,I)+TJLOFKR(:,J)
          ENDDO
        ENDDO
        SUPPHI(:IRMATCH,IB) =JLOFKR(:IRMATCH,NZEROS)
        SUPTPHI(:IRMATCH,IB)=TJLOFKR(:IRMATCH,NZEROS)
!
!!$IF(IB.EQ.1)OPEN(100,FILE='JLSD1.DAT')
!!$IF(IB.EQ.2)OPEN(100,FILE='JLSD2.DAT')
!!$IF(IB.EQ.3)OPEN(100,FILE='JLSD3.DAT')
!!$REWIND(100)
!!$DO IR=1,NR
!!$  WRITE(100,FMT='(20F20.8)')R(IR),PSPHI(IR,IB),JLOFKR(IR,:) &
!!$ &                               ,TPSPHI(IR,IB),TJLOFKR(IR,:)
!!$ENDDO
!!$CLOSE(100)
      ENDDO
!CALL WRITEPHI('TESTCONSTRUCTPSPHI2',GID,NR,NB,SUPPHI)
!STOP 'FORCED IN CONSTRUCTPSPHI'
      PSPHI(:IRMATCH,:)=SUPPHI(:IRMATCH,:)
      TPSPHI(:IRMATCH,:)=SUPTPHI(:IRMATCH,:)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE CONSTRUCTPSPHIOLD(GID,NR,RC,L,NB,PSPHI,TPSPHI)
!     **                                                                      **
!     **  PSEUDIZES A WAVE FUNCTION BY MATCHING SPHERICAL BESSEL FUNCTIONS    **
!     **                                                                      **
!     **  FIRST WE FIND THE VALUES K_I SO THAT THE SPHERICAL BESSEL           **
!     **  FUNCTION J_L(K_I*R) MATCHES THE LOGARITHMIC DERIVATIVE              **
!     **  OF THE FUNCTION TO BE PSEUDIZED.                                    **
!     **                                                                      **
!     **  THE FIRST TWO ARE THEN MATCHED SO THAT THE PHI AND TPHI             **
!     **  ARE CONTINUOUS. THE DERIVATIVE OF PHI IS AUTOMATICALLY              **
!     **  CONTINUOUS THROUGH THE CORRECT CHOICE OF THE K_I.                   **
!     **                                                                      **
!     **  ATTENTION! SINCE WE MATCH TO TAYLOR COEFFICIENTS, ONE SHOULD        **
!     **  ALSO MATCH ENERGY DERIVATIVES OF THE SPHERICAL BESSEL FUNCTIONS.    **
!     **                                                                      **
!     **  NOTE THAT THE TAYLOR EXPANSION OF THE WAVE FUNCTIONS IN ENERGY      **
!     **  BEHAVE DIFFERENT FROM THE ENERGY DEPENDENT SOLUTION.                **
!     **  THE TAYLOR EXPANSION CAN CREATE TWO NODES AT ONCE, WHICH            **
!     **  CAUSES DIFFICULTIES.                                                **
!     **                                                                      **
!     **                                                                      ** 
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)     :: GID
      INTEGER(4),INTENT(IN)     :: NR
      INTEGER(4),INTENT(IN)     :: L
      INTEGER(4),INTENT(IN)     :: NB
      REAL(8)   ,INTENT(IN)     :: RC
      REAL(8)   ,INTENT(INOUT)  :: PSPHI(NR,NB)
      REAL(8)   ,INTENT(INOUT)  :: TPSPHI(NR,NB)
      INTEGER(4),PARAMETER      :: NZEROS=3
      INTEGER(4),PARAMETER      :: NITER=100
      INTEGER(4),PARAMETER      :: NX=1000
      LOGICAL(4),PARAMETER      :: TDIFFERENTIABLELAPLACIAN=.FALSE.
      REAL(8)   ,PARAMETER      :: XMAX=30.D0
      REAL(8)                   :: R(NR)
      INTEGER(4)                :: IR,I,ITER
      REAL(8)                   :: TOL=1.D-12
      REAL(8)                   :: PI
      REAL(8)                   :: VAL,VALT,DERT
      REAL(8)                   :: JL(NX),DJLDR(NX),PHASE(NX)
      REAL(8)                   :: TJL(NX)
      REAL(8)                   :: JLOFKR(NR,NZEROS),TJLOFKR(NR,NZEROS)
      REAL(8)                   :: NN
      INTEGER(4)                :: ISTART,IBI
      REAL(8)                   :: X0,Y0,DX,XM,YM
      REAL(8)                   :: KI(NZEROS)
      REAL(8)                   :: A11,A12,A21,A22,DET
      REAL(8)                   :: AINV11,AINV12,AINV21,AINV22
      REAL(8)                   :: V1,V2,C1,C2,C3,C1A,C2A,C1B,C2B
      REAL(8)                   :: DTJ1,DTJ2,DTJ3
      REAL(8)                   :: SVAR1,SVAR2,FAC
      INTEGER(4)                :: IB
      REAL(8)                   :: SUPPHI(NR,NB)
      REAL(8)                   :: SUPTPHI(NR,NB)
      REAL(8)                   :: SHIFT
      REAL(8)                   :: X1,XDEX
      INTEGER(4)                :: GIDX
      INTEGER(4)                :: IRMATCH
      REAL(8)                   :: X(NX)
REAL(8) :: D2JLDR2(NR)
!     **************************************************************************
CALL WRITEPHI('TESTCONSTRUCTPSPHI1',GID,NR,NB,PSPHI)
      PI=4.D0*ATAN(1.D0)
      CALL RADIAL$NEW('SHLOG',GIDX)
      XDEX=1.D-5
      X1=XMAX/XDEX/REAL(NX)
      CALL RADIAL$SETR8(GIDX,'R1',X1)
      CALL RADIAL$SETR8(GIDX,'DEX',XDEX)
      CALL RADIAL$SETI4(GIDX,'NR',NX)
      CALL RADIAL$R(GIDX,NX,X)
!
!     ==========================================================================
!     ==  MAKE A COPY OF THE INPUT WAVE FUNCTIONS                             ==
!     ==========================================================================
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$XOFR(GID,RC,VAL)
      IRMATCH=INT(VAL)
      SUPPHI(:,:)=PSPHI(:,:)
      SUPTPHI(:,:)=TPSPHI(:,:)
!
!     =======================================================================
!     == CALCULATE BESSEL FUNCTION                                         ==
!     =======================================================================
      DO IR=1,NX
        CALL SPECIALFUNCTION$BESSEL(L,X(IR),JL(IR)) 
      ENDDO
!     == 2*TJL-JL=0 =======================================================
      TJL(:)=0.5D0*JL(:)   !KINETIC ENERGY TIMES BESSEL FUNCTION
!CALL RADIAL$VERLETD1(GID,NR,JL,DJLDR)
!CALL RADIAL$VERLETD2(GID,NR,JL,D2JLDR2)
!PRINT*,'JL ',L,(2.D0*DJLDR/R+D2JLDR2-REAL(L*(L+1),KIND=8)/R**2*JL)/JL

CALL WRITEPHI('JL',GIDX,NX,1,JL)
!
!     ==========================================================================
!     == CALCULATE PHASESHIFT OF BESSEL FUNCTION                              ==
!     ==========================================================================
      CALL RADIAL$VERLETD1(GIDX,NX,JL,DJLDR)
      NN=0.D0
      DO IR=2,NX
        IF(JL(IR)*JL(IR-1).LT.0.D0)NN=NN+1.D0
        PHASE(IR)=0.5D0-ATAN(X(IR)/RC*DJLDR(IR)/JL(IR))/PI+NN
      ENDDO
!     == AVOID GLITCH
      PHASE(1)=PHASE(2)
      PHASE(NX)=PHASE(NX-1)
CALL WRITEPHI('PHASE',GIDX,NX,1,PHASE)
!
!     ==========================================================================
!     == PSEUDIZE EVERY FUNCTION ON THE SUPPORT ENERGY GRID                   ==
!     ==========================================================================
      DO IB=1,NB
!       ========================================================================
!       == DETERMINE K-VALUES FOR WHICH THE LOGARITHIC DERIVATIVE MATCHES     ==
!       ========================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,SUPPHI(:,IB),RC,SHIFT)
        IF(SHIFT.LT.PHASE(1)) THEN
          CALL ERROR$MSG('MATCHING INCLUDING FIRST BESSEL FUNCTION NOT POSSIBLE')
          CALL ERROR$MSG('TRY TO INCREASE MATCHING RADIUS')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$STOP('CONSTRUCTPSPHI')
PRINT*,'WARNING SHIFT WAS REQUIRED ',L,IB,PHASE(1),SHIFT
           SHIFT=SHIFT+REAL(INT(PHASE(1)-SHIFT+1.D0))
        END IF
PRINT*,'SHIFT ',SHIFT
!
!       ========================================================================
!       == SEARCH BESSEL FUNCTIONS JL(K*R) WHICH HAVE IDENTICAL LOGARITHMIC   ==
!       == DERIVATIVES TO THE NODELESS FUNCTION AT THE RADIUS RCS             ==
!       ========================================================================
        DO I=1,NZEROS
          ISTART=1
          X0=0.4D0+1.D-3
          DX=0.1D0
          CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
          DO ITER=1,NITER
            CALL RADIAL$VALUE(GIDX,NX,PHASE,X0,Y0)
            Y0=Y0-SHIFT
            CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
            IF(X0.GT.XMAX) THEN
              CALL ERROR$MSG('RESULTING K OUT OF RANGE')
              CALL ERROR$STOP('CONSTRUCTPSPHI')
            END IF
            IF(ABS(DX).LT.TOL) EXIT
          ENDDO
          IF(ABS(DX).GT.TOL) THEN
            CALL ERROR$MSG('LOOP NOT CONVERGED')
            CALL ERROR$STOP('CONSTRUCTPSPHI')
          END IF
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
!
!       ========================================================================
!       == PLACE BESSELFUNCTIONS ON RADIAL GRID                               ==
!       ========================================================================
        DO I=1,NZEROS
          DO IR=1,NR
            CALL SPECIALFUNCTION$BESSEL(L,KI(I)*R(IR),JLOFKR(IR,I)) 
          ENDDO
          TJLOFKR(:,I)=0.5D0*KI(I)**2*JLOFKR(:,I)
          CALL RADIAL$VALUE(GID,NR,JLOFKR(:,I),RC,VAL)
          JLOFKR(:,I)=JLOFKR(:,I)/VAL
          TJLOFKR(:,I)=TJLOFKR(:,I)/VAL
        ENDDO
        CALL RADIAL$VALUE(GID,NR,SUPPHI(:,IB),RC,VAL)
        JLOFKR(:,:) =JLOFKR(:,:)*VAL
        TJLOFKR(:,:)=TJLOFKR(:,:)*VAL
IF(IB.EQ.1)OPEN(100,FILE='JLS1.DAT')
IF(IB.EQ.2)OPEN(100,FILE='JLS2.DAT')
IF(IB.EQ.3)OPEN(100,FILE='JLS3.DAT')
REWIND(100)
DO IR=1,NR
  WRITE(100,FMT='(20F20.8)')R(IR),SUPPHI(IR,IB),JLOFKR(IR,:)
ENDDO
CLOSE(100)
!
!       ========================================================================
!       == CONSTRUCT NEW PSEUDO WAVE FUNCTION                                 ==
!       ========================================================================
        CALL RADIAL$VALUE(GID,NR,SUPPHI(:,IB),RC,VAL)
        CALL RADIAL$VALUE(GID,NR,SUPTPHI(:,IB),RC,VALT)
        CALL RADIAL$DERIVATIVE(GID,NR,SUPTPHI(:,IB),RC,DERT)     
!       == MATCH K1 AND K2 =====================================================
        CALL RADIAL$VALUE(GID,NR,JLOFKR(:,1),RC,A11)
        CALL RADIAL$VALUE(GID,NR,JLOFKR(:,2),RC,A12)
        CALL RADIAL$VALUE(GID,NR,TJLOFKR(:,1),RC,A21)
        CALL RADIAL$VALUE(GID,NR,TJLOFKR(:,2),RC,A22)
        DET=A11*A22-A12*A21
        AINV11=+A22/DET
        AINV22=+A11/DET
        AINV12=-A12/DET
        AINV21=-A21/DET
        V1=VAL
        V2=VALT
        C1A=AINV11*V1+AINV12*V2
        C2A=AINV21*V1+AINV22*V2
        IF(TDIFFERENTIABLELAPLACIAN) THEN
!         ======================================================================
!         ==  USE A THIRD BESSELFUNCTION TO OBTAIN A DIFFERENTIABLE KINETIC   ==
!         ==  ENERGY DENSITY. THIS DOES NOT SEEM TO BE A GOOD IDEA            ==
!         ==  BECAUSE IT DETERIORATES THE SCATTERING PROPERTIES.              ==
!         == MATCH K2 AND K3                                                  ==
!         ======================================================================
          IF(KI(3)*RC.GT.R(NR-2)) THEN
            CALL ERROR$MSG('KI(3)*RC>RX')
            CALL ERROR$R8VAL('KI(3)*RC',KI(3)*RC)
            CALL ERROR$R8VAL('KI(3)',KI(3))
            CALL ERROR$R8VAL('RC',RC)
            CALL ERROR$R8VAL('R(NR-2)',R(NR-2))
            CALL ERROR$STOP('CONSTRUCTPSPHI')
          END IF
          CALL RADIAL$VALUE(GID,NR,JLOFKR(:,2),RC,A11)
          CALL RADIAL$VALUE(GID,NR,JLOFKR(:,3),RC,A12)
          CALL RADIAL$VALUE(GID,NR,TJLOFKR(:,2),RC,A21)
          CALL RADIAL$VALUE(GID,NR,TJLOFKR(:,3),RC,A22)
          DET=A11*A22-A12*A21
          AINV11=A22/DET
          AINV22=A11/DET
          AINV12=-A12/DET
          AINV21=-A21/DET
          V1=VAL
          V2=VALT
          C1B=AINV11*V1+AINV12*V2
          C2B=AINV21*V1+AINV22*V2
!         ======================================================================
!         == A SUPERPOSITION OF THE TWO SOLUTIONS WITH FACTORS                ==
!         == THAT ADD UP TO ONE, DO NOT AFFECT VALUE OR                       ==
!         == DERIVATIVE. USE THIS TO MATCH DERIVATIVE OF                      ==
!         == KINETIC ENERGY                                                   ==
!         ======================================================================
          CALL RADIAL$DERIVATIVE(GID,NR,TJLOFKR(:,1),RC,DTJ1)      
          CALL RADIAL$DERIVATIVE(GID,NR,TJLOFKR(:,2),RC,DTJ2)      
          CALL RADIAL$DERIVATIVE(GID,NR,TJLOFKR(:,3),RC,DTJ3)      
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
PRINT*,'C1,C2,C3 ',C1,C2,C3 
        SUPPHI(:,IB)=JLOFKR(:,1)*C1+JLOFKR(:,2)*C2+JLOFKR(:,3)*C3
        SUPTPHI(:,IB)=TJLOFKR(:,1)*C1+TJLOFKR(:,2)*C2+TJLOFKR(:,3)*C3
      ENDDO
!!$OPEN(100,FILE='JLS.DAT')
!!$REWIND(100)
!!$DO IR=1,NR
!!$  WRITE(100,FMT='(40F20.8)')R(IR),PSPHI(IR,:),SUPPHI(IR,:)
!!$ENDDO
!!$CLOSE(100)
!!$STOP
CALL WRITEPHI('TESTCONSTRUCTPSPHI2',GID,NR,NB,SUPPHI)
STOP 'FORCED IN CONSTRUCTPSPHI'
      PSPHI(:IRMATCH,:)=SUPPHI(:IRMATCH,:)
      TPSPHI(:IRMATCH,:)=SUPTPHI(:IRMATCH,:)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
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
      INTEGER(4) ,INTENT(IN)     :: POW
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
      PI=4.D0*ATAN(1.D0)
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
PRINT*,'CONSTRUCTPSPHI1',IB
!
!       =====================================================================
!       ==  DETERMINE PHASESHIFT OF PSPHI                                  ==
!       =====================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSPHI(:,IB),RCMATCH,PHASE)
!
!       =====================================================================
!       ==  FIND POTENTIAL THAT PRODUCES CORRECT PHASESHIFT                ==
!       =====================================================================
        ISTART=-1
        X0=PSPOT(1)
        DX=0.1D0
        CALL BISEC(ISTART,IBI,X0,F0,DX,XM,FM)
        DO ITER=1,NITER
          VAL0=X0
          CALL PSEUDIZE(GID,NR,POW,VAL0,RC,PSPOT,POT)
          DREL=0.D0
          G=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,0,G,L,E(IB),1,PHI(:))
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,RCMATCH,PHASE1)
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
!
!     ......................................................................
      SUBROUTINE MAKEDATH(GID,NR,RBOX,AEPOT,PSPOT,L,NAUG,AEPHI,PSPHI,DT,DH)
!     **                                                                  **
!     **  ONE-CENTER HAMILTONIAN                                          **
!     **  DH=DT+<AEPHI|V|AEPHI>-<PSPHI|VTILDE|PSPHI>                      **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID            ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR             ! #(RADIAL GRID POINTS)
      REAL(8)    ,INTENT(IN)     :: RBOX           ! BOX RADIUS
      INTEGER(4) ,INTENT(IN)     :: L              ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: NAUG           ! #(PARTIAL WAVES)
      REAL(8)    ,INTENT(IN)     :: AEPOT(NR)      ! AE-POTENTIAL
      REAL(8)    ,INTENT(IN)     :: PSPOT(NR)      ! PS-POTENTIAL
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,NAUG) ! PS-PARTIAL WAVE
      REAL(8)    ,INTENT(IN)     :: AEPHI(NR,NAUG) ! AE PARTIAL WAVE
      REAL(8)    ,INTENT(IN)     :: DT(NAUG,NAUG)  ! 1C KINETIC ENERGY
      REAL(8)    ,INTENT(OUT)    :: DH(NAUG,NAUG)  ! 1C-HAMILTONIAN
      REAL(8)                    :: R(NR)
      REAL(8)                    :: AUX(NR),AUX1(NR)
      INTEGER(4)                 :: IB,JB
      REAL(8)                    :: SVAR
      REAL(8)                    :: PI,Y0
!     ***********************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      DO IB=1,NAUG
        DO JB=1,IB
          AUX(:)=AEPHI(:,IB)*AEPOT(:)*AEPHI(:,JB)-PSPHI(:,IB)*PSPOT(:)*PSPHI(:,JB)
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
          SVAR=0.5D0*(DT(IB,JB)+DT(JB,IB))+SVAR*Y0
          DH(IB,JB)=SVAR
          DH(JB,IB)=SVAR
        ENDDO
      ENDDO
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE DTDO(GID,NR,RBOX,NAUG,AEPHI,TAEPHI,PSPHI,TPSPHI,DT,DO)
!     **                                                                  **
!     **  ONE-CENTER KINETIC ENERGY AND ONE-CENTER OVERLAP                **
!     **     DT=<AEPHI|T|AEPHI>-<PSPHI|T|PSPHI>                           **
!     **     DO=<AEPHI|AEPHI>-<PSPHI|PSPHI>                               **
!     **  FROM THE PARTIAL WAVES. T IS THE KINETIC ENERGY OPERATOR        **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      REAL(8)    ,INTENT(IN)     :: RBOX           ! BOX RADIUS
      INTEGER(4) ,INTENT(IN)     :: NAUG
      REAL(8)    ,INTENT(IN)     :: AEPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: TAEPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: TPSPHI(NR,NAUG)
      REAL(8)    ,INTENT(OUT)    :: DT(NAUG,NAUG)
      REAL(8)    ,INTENT(OUT)    :: DO(NAUG,NAUG)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: R2(NR)
      REAL(8)                    :: AUX(NR),AUX1(NR)
      REAL(8)                    :: SVAR
      INTEGER(4)                 :: IB,JB
      LOGICAL(4)                 :: TTEST=.TRUE.
!     *********************************************************************
      CALL RADIAL$R(GID,NR,R)
      R2(:)=R(:)**2
      DO IB=1,NAUG
        DO JB=1,NAUG
          AUX(:)=(AEPHI(:,IB)*TAEPHI(:,JB)-PSPHI(:,IB)*TPSPHI(:,JB))*R2(:)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
          DT(IB,JB)=SVAR
          AUX(:)=(AEPHI(:,IB)*AEPHI(:,JB)-PSPHI(:,IB)*PSPHI(:,JB))*R2(:)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
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
      SUBROUTINE COREORTHOGONALIZE(GID,NR,RBOX,NC,LOFI,SOFI,UOFI,TUOFI &
     &                            ,L,NAUG,UPHI,TUPHI,AEPHI,TAEPHI)
!     **                                                                  **
!     **  ORTHOGONALIZE NODE-LESS WAVE FUNCTIONS TO THE CORE              **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      REAL(8)    ,INTENT(IN)     :: RBOX
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
      REAL(8)                    :: AUX(NR),AUX1(NR)
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
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
          PHIC(:,IB) = PHIC(:,IB)- PHIC(:,JB)*SVAR
          TPHIC(:,IB)=TPHIC(:,IB)-TPHIC(:,JB)*SVAR
        ENDDO
        AUX(:)=(PHIC(:,IB)*R(:))**2
        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
        CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
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
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
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
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
            PRINT*,'<UCORE|AEPHI> ',IB,JB,SVAR
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE BIORTHOMATRICES(GID,NR,RBOX,L,NAUG,PSPHI,PRO,TRANSPHI,TRANSPRO)
!     **                                                                  **
!     ** DETERMINES THE MATRICES TRANSPHI AND TRANSPRO SUCH THAT          **
!     **     |PHI-BAR>:=|PHI>TRANSSPHI                                    **
!     **     |PRO-BAR>:=|PRO>TRANSSPRO                                    **
!     **  OBEY  <PHIBAR(I)|PROBAR(J)>=DELTA(I,J)    (KRONECKER DELTA)     **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      REAL(8)    ,INTENT(IN)     :: RBOX
      INTEGER(4) ,INTENT(IN)     :: NAUG
      INTEGER(4) ,INTENT(IN)     :: L
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: PRO(NR,NAUG)
      REAL(8)    ,INTENT(OUT)    :: TRANSPHI(NAUG,NAUG)
      REAL(8)    ,INTENT(OUT)    :: TRANSPRO(NAUG,NAUG)
      INTEGER(4)                 :: IB,JB
      REAL(8)                    :: AUX(NR),AUX1(NR)
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
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)       
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)       
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
DO JB=1,NAUG
WRITE(*,FMT='(20E10.2)')TRANSPROPSI(:,JB)
ENDDO
        CALL ERROR$MSG('BIORTHOGONALIZATION INACCURATE')
        CALL ERROR$R8VAL('MAX. DEV.',SVAR)
        CALL ERROR$STOP('BIORTHOMATRICES')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BIORTHOMATRICES_TAYLOR(GID,NR,L,NAUG,NNEND,RCOV,PSPHI,PRO,TRANSPHI,TRANSPRO)
!     **                                                                      **
!     ** DETERMINES THE MATRICES TRANSPHI AND TRANSPRO SUCH THAT              **
!     **     |PHI-BAR>:=|PHI>TRANSSPHI                                        **
!     **     |PRO-BAR>:=|PRO>TRANSSPRO                                        **
!     **  OBEY  <PHIBAR(I)|PROBAR(J)>=DELTA(I,J)    (KRONECKER DELTA)         **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NAUG
      INTEGER(4) ,INTENT(IN)     :: L
      INTEGER(4) ,INTENT(IN)     :: NNEND
      REAL(8)    ,INTENT(IN)     :: RCOV
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
      REAL(8)                    :: PROPSIINV(NAUG,NAUG)
      REAL(8)                    :: TRANSPROPSI(NAUG,NAUG)
      REAL(8)                    :: VAL,VALDOT,DER,DERDOT
      REAL(8)                    :: BANDWIDTH
!     *********************************************************************
      CALL RADIAL$R(GID,NR,R)
      IF(NNEND+1.GT.NAUG) THEN
        CALL ERROR$MSG('NNEND TOO LARGE (NNEND+1>NAUG)')
        CALL ERROR$I4VAL('NNEND',NNEND)
        CALL ERROR$I4VAL('NAUG',NAUG)
        CALL ERROR$STOP('BIORTHOMATRICES_TAYLOR')
      END IF
!     =====================================================================
!     == RESCALE PARTIAL WAVES
!     =====================================================================
      CALL RADIAL$VALUE(GID,NR,PSPHI(:,NNEND),RCOV,VAL)
      CALL RADIAL$VALUE(GID,NR,PSPHI(:,NNEND+1),RCOV,VALDOT)
      CALL RADIAL$DERIVATIVE(GID,NR,PSPHI(:,NNEND),RCOV,DER)
      CALL RADIAL$DERIVATIVE(GID,NR,PSPHI(:,NNEND+1),RCOV,DERDOT)
      BANDWIDTH=-VAL/VALDOT+DER/DERDOT
PRINT*,'BANDWIDTH',L,BANDWIDTH,VAL,DER,VALDOT,DERDOT
      TRANSPHI(:,:)=0.D0
      DO IB=1,NNEND
        TRANSPHI(IB,IB)=1.D0
      ENDDO
      DO IB=NNEND+1,NAUG
        CALL RADIAL$VALUE(GID,NR,PSPHI(:,IB),RCOV,SVAR)
        TRANSPHI(IB,IB)=1.D0 !/SVAR
      ENDDO
!
!     =====================================================================
!     == CALCULATE INITIAL VIOLATION OF BIORTHOGONALITY                  ==
!     == COLLECT TRANSFORMATION MATRIX BETWEEN NEW AND OLD               ==
!     =====================================================================
      DO IB=1,NAUG
        DO JB=1,NAUG
          AUX(:)=PRO(:,IB)*PSPHI(:,JB)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)       
          PROPSI(IB,JB)=SVAR
        ENDDO
      ENDDO
      CALL LIB$INVERTR8(NAUG,TRANSPOSE(MATMUL(PROPSI,TRANSPHI)),TRANSPRO)
DO IB=1,NAUG
 WRITE(*,FMT='("PROPSI",10F20.10)')PROPSI(IB,:)
ENDDO
DO IB=1,NAUG
 WRITE(*,FMT='("TRANSPRO",10F20.10)')TRANSPRO(IB,:)
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
DO JB=1,NAUG
WRITE(*,FMT='(20E10.2)')TRANSPROPSI(:,JB)
ENDDO
        CALL ERROR$MSG('BIORTHOGONALIZATION INACCURATE')
        CALL ERROR$R8VAL('MAX. DEV.',SVAR)
        CALL ERROR$STOP('BIORTHOMATRICES_TAYLOR')
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
      PI=4.D0*ATAN(1.D0)
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
      SUBROUTINE BAREPROJECTORS_TAYLOR(GID,NR,L,E,PSPOT,RC &
     &                         ,NAUG,NNEND,PSPHI,TPSPHI,PRO)
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
      INTEGER(4) ,INTENT(IN)     :: NNEND
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,NAUG)
      REAL(8)    ,INTENT(IN)     :: TPSPHI(NR,NAUG)
      REAL(8)    ,INTENT(OUT)    :: PRO(NR,NAUG)
      INTEGER(4)                 :: IB,IR
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: R(NR)
      LOGICAL(4)                 :: TTEST=.TRUE.
!     *********************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      DO IB=1,1
        PRO(:,IB)=TPSPHI(:,IB)+(PSPOT(:)*Y0-E(IB))*PSPHI(:,IB)
        PRO(NR,IB)=2.D0*PRO(NR-1,IB)-PRO(NR-2,IB)
      ENDDO
      DO IB=2,NAUG
        PRO(:,IB)=TPSPHI(:,IB)+(PSPOT(:)*Y0-E(IB))*PSPHI(:,IB)-PSPHI(:,IB-1)
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
      INTEGER(4)                 :: NITER=5000
      REAL(8)                    :: XAV,XMAX
      LOGICAL(4)                 :: CONVG
      REAL(8)   ,PARAMETER       :: TOL=1.D-5
      INTEGER(4)                  :: NFILO
      REAL(8)                    :: EH,EXC
      LOGICAL(4),PARAMETER       :: TBROYDEN=.TRUE.
      REAL(8)                    :: POTIN(NR)
      INTEGER(4)                 :: IRBOX
      INTEGER(4)                 :: IR
      CHARACTER(32)              :: ID
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
        CALL BROYDEN$NEW(NR,4,1.D0)
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
!     ......................................................................
      SUBROUTINE NODELESS(ID,GID,NR,NB,LOFI,SO,NN,RBOX,DREL,POT,PHI,TPHI,E)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NB        ! #(STATES)
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)  !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO(NB)    !SWITCH FOR SPIN-ORBIT COUP.
      INTEGER(4) ,INTENT(IN)     :: NN(NB)    !#(NODES)
      REAL(8)    ,INTENT(IN)     :: RBOX      !BOX RADIUS
      REAL(8)    ,INTENT(IN)     :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: POT(NR)   !POTENTIAL
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,NB)   !WAVE FUNCTIONS
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR,NB)  !KINETIC ENERGY*WAVE FUNCTIONS
      REAL(8)    ,INTENT(INOUT)    :: E(NB)     !ONE-PARTICLE ENERGIES
      REAL(8)                    :: R(NR)
      REAL(8)                    :: G(NR)
      REAL(8)                    :: DREL1(NR)
      REAL(8)                    :: AUX(NR)
      INTEGER(4)                 :: IB,IR,JB
      REAL(8)                    :: SVAR
      REAL(8)                    :: PI,Y0
      LOGICAL(4)                 :: TCHK
!     ***********************************************************************
      IF(ID.NE.'FULL'.AND.ID.NE.'EFFZORA'.AND.ID.NE.'NONREL') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$MSG('ALLOWED VALUES ARE "FULL",EFFZORA","NONREL"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('NODELESS')
      ENDIF
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
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
          CALL SCHROEDINGER$DREL(GID,NR,POT,E(IB),DREL1)
        ELSE IF(ID.EQ.'NONREL') THEN
          DREL1(:)=0.D0 
        ELSE IF(ID.EQ.'EFFZORA') THEN
          DREL1(:)=DREL(:)
        END IF
        CALL BOUNDSTATE(GID,NR,LOFI(IB),SO(IB),RBOX,DREL1,G,0,POT,E(IB),PHI(:,IB))
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
      SUBROUTINE ONENODELESS(ID,GID,NR,L,NN,SO,RBOX,DREL,POT,PHIC,E,PHI,TPHI)
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
      REAL(8)    ,INTENT(IN)     :: RBOX      !BOX RADIUS
      REAL(8)    ,INTENT(INOUT)  :: E         !EXPANSION ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE FUNCTIONS
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR)!KINETIC ENERGY*WAVE FUNCTIONS
      INTEGER(4)                 :: IB,IR
      INTEGER(4)                 :: IRBOX
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: R(NR)
!     ***********************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GT.RBOX) EXIT
      ENDDO
      IF(ID.EQ.'BOUND') THEN
        CALL BOUNDSTATE(GID,NR,L,SO,RBOX,DREL,PHIC,NN,POT,E,PHI)
      ELSE IF(ID.EQ.'SCATT') THEN
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,PHIC,L,E,1,PHI)
        PHI(IRBOX+2:)=0.D0 !CHOP TAILS OFF BUT LEAVE ENOUGH FOR DERIVATIVES
      END IF
      TPHI(:)=PHIC(:)+(E-POT(:)*Y0)*PHI(:)
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
      REAL(8)    ,INTENT(IN)     :: RBOX      !ATOM ENCLOSED WITHIN RADIUS RBOX
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
!     ==  R(IRBOX) IS THE FIRST GRIDPOINT JUST OUTSIDE THE BOX
      IRBOX=1
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GT.RBOX) EXIT
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
!        POT1(IROUT:)=POT(IROUT)
        POT1(IROUT:)=E
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
    OPEN(UNIT=8,FILE='XXX.DAT')
    DO I=1,NR
      WRITE(8,*)I,R(I),PHIHOM(I),POT1(I)*Y0
    ENDDO
    CLOSE(8)      
    CALL ERROR$MSG('PHI CONTAINS NANS')
    CALL ERROR$STOP('BOUNDSTATE')
  END IF
ENDDO
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE BOUNDSTATEY(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,PHI)
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
      REAL(8)                    :: VAL,VAL1,VAL2
      INTEGER(4)                 :: IROUT,IRCL,IRBOX
!     *********************************************************************
                                 CALL TRACE$PUSH('BOUNDSTATE')
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      IRBOX=1
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GT.RBOX) EXIT
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
        ROUT=MIN(R(IROUT),RBOX)
        IROUT=MIN(IROUT,IRBOX)
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
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
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
        IDIR=-1
        IF(IROUT.EQ.IRBOX.AND.IROUT+4.LT.NR) THEN
          GHOM(:)=0.D0
          GHOM(IROUT+3)=1.D-8
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
          GHOM(:)=0.D0
          GHOM(IROUT+4)=1.D-8
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIINHOM)
          CALL RADIAL$VALUE(GID,NR,PHIHOM,RBOX,VAL1)
          CALL RADIAL$VALUE(GID,NR,PHIINHOM,RBOX,VAL2)
          SVAR=VAL1+VAL2
          VAL1=VAL1/SVAR
          VAL2=VAL2/SVAR
          PHIHOM(:)=VAL2*PHIHOM(:)-VAL1*PHIINHOM(:)
        ELSE
          GHOM(:)=0.D0
          GHOM(IROUT)=1.D-8
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
        END IF
        IF(.NOT.THOM) THEN     
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
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
      PHI(IROUT:)=0.D0
!
DO IR=1,NR
  IF(.NOT.(PHI(IR).GT.0.D0.OR.PHI(IR).LE.0.D0)) THEN
    PRINT*,'ERROR'
!    PRINT*,'PHIIN',PHI(:IRMATCH-1)
!    PRINT*,'PHIOUT',PHI(IRMATCH:)
    PRINT*,'SVAR ',SVAR
    PRINT*,'IROUT,IRCL ',IROUT,IRCL,NR
    PRINT*,'R,PHI ',R(IR),PHI(IR)
    OPEN(UNIT=8,FILE='XXX.DAT')
    DO I=1,NR
      WRITE(8,*)I,R(I),PHIHOM(I),POT1(I)*Y0
    ENDDO
    CLOSE(8)      
    CALL ERROR$MSG('PHI CONTAINS NANS')
    CALL ERROR$STOP('BOUNDSTATE')
  END IF
ENDDO
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE BOUNDSTATEX(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,PHI)
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
      REAL(8)                    :: PHIHOM(NR),PHIHOM2(NR),PHIINHOM(NR),GHOM(NR)
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
      REAL(8)                    :: VAL,VAL1,VAL2
      INTEGER(4)                 :: IROUT,IRCL,IRBOX
!     *********************************************************************
                                 CALL TRACE$PUSH('BOUNDSTATE')
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      IRBOX=1
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GT.RBOX) EXIT
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
        ROUT=MIN(R(IROUT),RBOX)
        IROUT=MIN(IROUT,IRBOX)
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
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
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
        IDIR=-1
        GHOM(:)=0.D0
        IR=MIN(IROUT+3,NR)
        GHOM(IR)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
        GHOM=0.D0
        IR=MIN(IROUT+10,NR)
        GHOM(IR)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM2)
        CALL RADIAL$VALUE(GID,NR,PHIHOM,RBOX,VAL1)
        CALL RADIAL$VALUE(GID,NR,PHIHOM2,RBOX,VAL2)
        IF(VAL2.NE.0.D0) THEN
          SVAR=-VAL1/VAL2
          PHIHOM(:)=PHIHOM(:)+SVAR*PHIHOM2(:)
        END IF
        IF(.NOT.THOM) THEN     
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
          IF(VAL2.NE.0.D0) THEN
            CALL RADIAL$VALUE(GID,NR,PHIINHOM,RBOX,VAL1)
            SVAR=-VAL1/VAL2
            PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM2(:)
          END IF
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
      PHI(IROUT:)=0.D0
!
DO IR=1,NR
  IF(.NOT.(PHI(IR).GT.0.D0.OR.PHI(IR).LE.0.D0)) THEN
    PRINT*,'ERROR'
!    PRINT*,'PHIIN',PHI(:IRMATCH-1)
!    PRINT*,'PHIOUT',PHI(IRMATCH:)
    PRINT*,'SVAR ',SVAR
    PRINT*,'IROUT,IRCL ',IROUT,IRCL,NR
    PRINT*,'R,PHI ',R(IR),PHI(IR)
    OPEN(UNIT=8,FILE='XXX.DAT')
    DO I=1,NR
      WRITE(8,*)I,R(I),PHIHOM(I),POT1(I)*Y0
    ENDDO
    CLOSE(8)      
    CALL ERROR$MSG('PHI CONTAINS NANS')
    CALL ERROR$STOP('BOUNDSTATE')
  END IF
ENDDO
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO,DH,DO,G,PHI)
!     **************************************************************************
!     **                                                                      **
!     **  SOLVES THE RADIAL PAW -SCHROEDINGER EQUATION.                       **
!     **    (T+VTILDE-E+|P>(DH-E*DO<P|]|PHI>=|G>                              **
!     **  WHERE T IS THE NONRELATIVISTIC KINETIC ENERGY.                      **
!     **                                                                      **
!     **    DH=<AEPHI|T+V|AEPHI>-<PSPHI|T+VTILDE|PSPHI>                       **
!     **    DO=<AEPHI|AEPHI>-<PSPHI|PSPHI>                                    **
!     **                                                                      **
!     **************************************************************************
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
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  -1/2NABLA^2+POT-E|U>=|G>                                    ==
!     ==================================================================
      SO=0
      DREL(:)=0.D0
      CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,G,L,E,1,U)
!
!     ==================================================================
!     ==  -1/2NABLA^2+POT-E|V>=|PRO>                                 ==
!     ==================================================================
      DO I1=1,NPRO
        CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,PRO(:,I1),L,E,1,V(:,I1))
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
      SUBROUTINE PSETOTSHELL(GID,NR,AEZ,RBOX,PSRHOC,AERHOC,PSPOT &
     &                      ,NB,NC,LMAX,NAUG,NPRO,LOFI,FOFI,EOFI &
     &                      ,PRO1,AEPHI1,PSPHI1,RCSM,DTKIN1,DH1,DO1,VADD)
!     **                                                               **
!     **                                                               **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: AEZ
      REAL(8),   INTENT(IN) :: RBOX
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
      INTEGER(4)            :: L1,L2,LN1,LN2,IPRO,IPRO1,IPRO2,IB,I,IR
      INTEGER(4)            :: L,NN
      INTEGER(4)            :: IRMATCH
      INTEGER(4)            :: IRBOX
      REAL(8)               :: DREL(NR)
      REAL(8)               :: AUX(NR)
      INTEGER(4),ALLOCATABLE:: LOX(:)
      REAL(8)               :: SVAR
      REAL(8)               :: PROJ(NAUG)
      REAL(8)   ,ALLOCATABLE:: DTKIN(:,:)
      REAL(8)   ,ALLOCATABLE:: PRO(:,:)
      REAL(8)   ,ALLOCATABLE:: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE:: PSPHI(:,:)
      INTEGER(4)            :: NFILO
      CHARACTER(64)         :: STRING
!     *******************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL REPORT$TITLE(NFILO,'EIGENSTATES FROM PAW CALCULATION')
      IRBOX=1
      DO IR=1,NR-1
        IRBOX=IR
        IF(R(IR).GT.RBOX) EXIT
      ENDDO
!
!     ===================================================================================
!     == FIND PAW BOUND STATE BY BISECTION                                             ==
!     ===================================================================================
IF(NITER.EQ.1) WRITE(NFILO,*) 'NO BAND ENERGY SEARCH IN PSETOTSHELL!!!'
      ALLOCATE(PSPSI(NR,NB-NC))
      ALLOCATE(TPSPSI(NR,NB-NC))
      DO IB=NC+1,NB
        L=LOFI(IB)
!
        DF0=1.D0
        DO I=NC+1,IB-1
          IF(LOFI(I).EQ.L) DF0=DF0+1.D0
        ENDDO
!
!       =================================================================================
!       == DETERMINE IRMATCH, WHERE THE OUTWARD AND INWARD SOLUTIONS ARE MATCHED ========
!       =================================================================================
        IRMATCH=10.D0
        DO IR=1,NR
          IF(R(IR).GT.3.D0) EXIT
          DO I=1,NPRO(L+1)
            IF(ABS(PRO1(IR,I,L+1)).GT.1.D-12) IRMATCH=IR
          ENDDO
        ENDDO
!
!       =================================================================================
!       == FIND PAW BOUND STATE BY BISECTION                                           ==
!       =================================================================================
        ISTART=1
        X0=EOFI(IB)
        DX=0.1D0
        CALL BISEC(ISTART,IBI,X0,F0,DX,XM,YM)
!PRINT*,'IB ',IB,L,DF0,EOFI(IB),RBOX
        DO I=1,NITER
          E=X0
          G(:)=0.D0
          CALL PAWDER(GID,NR,L,E,PSPOT,NPRO(L+1),PRO1(:,1:NPRO(L+1),L+1) &
                    ,DH1(1:NPRO(L+1),1:NPRO(L+1),L+1),DO1(1:NPRO(L+1),1:NPRO(L+1),L+1) &
       &            ,G,PSPSI(:,IB-NC))
!WRITE(STRING,*)IB
!STRING='PAWTEST'//TRIM(ADJUSTL(STRING))//'.DAT'         
!CALL WRITEPHI('PAWTEST2.DAT',GID,NR,1,PSPSI(:,IB-NC))
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSPSI(:,IB-NC),RBOX,F0)
!PRINT*,'MARKE 1B',I,X0,F0-DF0
!IF(I.EQ.10)STOP
          F0=F0-DF0
          CALL BISEC(ISTART,IBI,X0,F0,DX,XM,YM)
          IF(ABS(DX).LE.TOL) EXIT
        ENDDO
        IF(ABS(DX).GT.TOL.AND.NITER.GT.1) THEN
          CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
          CALL ERROR$MSG('BOUND STATE NOT FOUND')
          CALL ERROR$STOP('PSETOTSHELL')
        END IF
        WRITE(NFILO,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2," E[H]=",F15.9)') &
     &              L,0,0,FOFI(IB),E
        WRITE(NFILO,*)NPRO(L+1),F0
!
!       == INTEGRATE INWARD TO AVOID EXPONENTIALLY INCREASING FUNCTION
        IF(IRMATCH.LT.IRBOX) THEN
          G=0.D0
          G(IRBOX)=1.D-8
          DREL=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,0,G,L,E,-1,AUX)
          PSPSI(IRMATCH:,IB-NC)=AUX(IRMATCH:)*PSPSI(IRMATCH,IB-NC)/AUX(IRMATCH)
        ELSE
          PSPSI(IRBOX+1:,IB-NC)=0.D0
        END IF
!
!       == NORMALIZE WAVE FUNCTION
        DO IPRO=1,NPRO(L+1)
          AUX(:)=PSPSI(:,IB-NC)*PRO1(:,IPRO,L+1)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,PROJ(IPRO))
        ENDDO
        AUX(:)=PSPSI(:,IB-NC)**2*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        SVAR=SVAR+DOT_PRODUCT(PROJ(:NPRO(L+1)) &
     &           ,MATMUL(DO1(:NPRO(L+1),:NPRO(L+1),L+1),PROJ(:NPRO(L+1))))
        IF(SVAR.LE.0.D0) THEN
          CALL ERROR$MSG('NEGATIVE NORM')
          CALL ERROR$I4VAL('IB',IB)
          CALL ERROR$R8VAL('SVAR',SVAR)
          CALL ERROR$STOP('PSETOTSHELL')
        END IF
        PSPSI(:,IB-NC)=PSPSI(:,IB-NC)/SQRT(SVAR)
        PROJ(:)=PROJ(:)/SQRT(SVAR)

IF(NPRO(L+1).GT.NAUG) THEN
  CALL ERROR$MSG('ARRAY BOUNDS VIOLATED')
  CALL ERROR$I4VAL('NPRO',NPRO(L+1))
  CALL ERROR$I4VAL('NAUG',NAUG)
  CALL ERROR$STOP('PSETOTSHELL')
END IF
IF(L+1.GT.LMAX+1) THEN
  CALL ERROR$MSG('ARRAY BOUNDS VIOLATED')
  CALL ERROR$I4VAL('L',L)
  CALL ERROR$I4VAL('LMAX',LMAX)
  CALL ERROR$STOP('PSETOTSHELL')
END IF
        PROJ(:NPRO(L+1))=MATMUL(DH1(:NPRO(L+1),:NPRO(L+1),L+1)-E*DO1(:NPRO(L+1),:NPRO(L+1),L+1) &
     &                ,PROJ(1:NPRO(L+1)))
        TPSPSI(:,IB-NC)=(E-PSPOT(:)*Y0)*PSPSI(:,IB-NC) &
     &                 -MATMUL(PRO1(:,:NPRO(L+1),L+1),PROJ(:NPRO(L+1)))
      ENDDO
CALL WRITEPHI('PSPSIFROMPSETOTSHELL.DAT',GID,NR,NB-NC,PSPSI)
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
PRINT*,'LNX ',LNX
      CALL PSETOT(GID,NR,AEZ,PSRHOC,AERHOC &
     &                 ,NB-NC,LOFI(NC+1:NB),FOFI(NC+1:NB),PSPSI,TPSPSI &
     &                 ,LNX,LOX,PRO,AEPHI,PSPHI,RCSM,DTKIN,VADD)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE PSETOT(GID,NR,AEZ,PSRHOC,AERHOC,NB,LOFI,FOFI,PSPSI,TPSPSI &
     &                 ,LNX,LOX,PRO,AEPHI,PSPHI,RCSM,DTKIN,VADD)
!     **************************************************************************
!     **                                                               **
!     **                                                               **
!     **************************************************************************
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
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL WRITEPHI('TEST.DAT',GID,NR,NB,PSPSI)
!
!     ==========================================================================
!     == DENSITY MATRIX                                                       ==
!     ==========================================================================
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
!     ==========================================================================
!     == MULTIPOLE MOMENT                                                     ==
!     ==========================================================================
      AUX(:)=(AERHOC(:)-PSRHOC(:))/Y0
      DO LN1=1,LNX
        DO LN2=1,LNX
          AUX(:)=AUX(:)+DENMAT(LN1,LN2) &
     &                 *(AEPHI(:,LN1)*AEPHI(:,LN2)-PSPHI(:,LN1)*PSPHI(:,LN2))
        ENDDO
      ENDDO
      AUX(:)=AUX(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      QLM=(SVAR-AEZ)/Y0
!TAKE CARE OF FACTOR YLM!!!!
PRINT*,'COMPENSATION CHARGE ',QLM*Y0
!
!     ==========================================================================
!     == KINETIC ENERGY AND CHARGE DENSITY                                    ==
!     ==========================================================================
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
PRINT*,'FULL 1-C HARTREE ENERGY ',AEEH1,AEZ
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
!.............................................................
subroutine testintegral(string,gid,nr,aecore)
implicit none
character(*),intent(in) :: string
integer(4)  ,intent(in) :: gid
integer(4)  ,intent(in) :: nr
real(8)     ,intent(in) :: aecore(nr)
real(8)                 :: aux(nr)
real(8)                 :: r(nr)
real(8)                 :: val
real(8)                 :: pi,y0
pi=4.d0*atan(1.d0)
y0=1.d0/sqrt(4.d0*pi)
call radial$r(gid,nr,r)
aux=r**2*aecore/y0
call radial$integral(gid,nr,aux,val)
print*,string,' integral of aecore ',val
return
end
!.............................................................i
subroutine testintegral2(string,gid,nr,phi)
implicit none
character(*),intent(in) :: string
integer(4)  ,intent(in) :: gid
integer(4)  ,intent(in) :: nr
real(8)     ,intent(in) :: phi(nr)
real(8)                 :: aux(nr)
real(8)                 :: r(nr)
real(8)                 :: val
real(8)                 :: pi,y0
pi=4.d0*atan(1.d0)
y0=1.d0/sqrt(4.d0*pi)
call radial$r(gid,nr,r)
aux=r**2*phi**2
call radial$integral(gid,nr,aux,val)
print*,string,' integral of phi**2 ',val
return
end
