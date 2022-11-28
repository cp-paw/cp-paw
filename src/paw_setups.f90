!===============================================================================
!===============================================================================
!===============================================================================
!====      FORMFACTORS                                                      ====
!===============================================================================
!===============================================================================
!===============================================================================
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE SETUP_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: SETUP                                                              **
!**                                                                           **
!**  PURPOSE:                                                                 **
!**                                                                           **
!**  READS PSEUDOPOTENTIALS AND ATOMIC PSEUDO-WAVEFUNCTIONS                   **
!**  GIVEN ON A LINEAR, RADIAL MESH                                           **
!**              AND                                                          **
!**  CALCULATES THE FORMFACTORS FOR LOCAL AND NONLOCAL                        **
!**  CONTRIBUTIONS TO THE PSEUDOPOTENTIAL                                     **
!**  IN PLANE WAVE REPRESENTATION                                             **
!**                                                                           **
!**  THE LOCAL COMPONENT CONTAINS ALSO THE POTENTIAL OF A                     **
!**  (GAUSSIAN SHAPED) CHARGEDENSITY, WHICH COMPENSATES                       **
!**  THE IONIC CHARGE.                                                        **
!**                                                                           **
!**                                                                           **
!**  RCBG (R) = -ZATOM/(SQRT(PI)*RCRHO)**3 * EXP(-(R/RCRHO)**2)               **
!**                                                                           **
!**  OUTPUT:                                                                  **
!**    RCRHO     RADIUS OF GAUSSIAN IN RHOPS                                  **
!**    RHOPS     COMPENSATION CHARGE DENSITY                                  **
!**    VLOC      LOCAL CONTRIBUTION TO PSEUDOPOTENTIAL                        **
!**    WNL       NONLOCAL CONTRIBUTION (PROJECTOR)                            **
!**                                                                           **
!**                                                                           **
!**                                                                           **
!**  LOCAL ORBITALS:                                                          **
!**    LOCAL ORBITALS ARE APPROXIMATIONS FOR WANNIER FUNCTIONS FOR TREATING   **
!**    LOCAL CORRELATION EFFECTS AND FOR CHEMICAL BOND ANALYSIS               **
!**    THEY ARE PARTIAL WAVES SUPERIMPOSED SO THAT THEY HAVE A NODE           **
!**    AT A GIVEN RADIUS.                                                     **
!**                                                                           **
!**      |CHI_I>=\SUM_J |PHI_J> AMAT(J,I)                                     **
!**      |PSI>=\SUM_J |CHI_J> BMAT(J,I) <PRO|PSITILDE>                        **
!**                                                                           **
!**    SETR4('RADCHI',RAD)                                                    **
!**    SETI4A('NOFLCHI',4,NOFL)                                               **
!**    GETI4('LNXCHI',LNXCHI)                                                 **
!**    GETI4A('LOXCHI',LNXCHI,LOXCHI)                                         **
!**    GETR8A('AMATCHI',LNXPHI*LNXCHI,AMAT)                                   **
!**    GETR8A('BMATCHI',LNXCHI*LNXPHI,BMAT)                                   **
!**                                                                           **
!**  AN UNSELECT FUNCTION HAS BEEN INCLUDED BUT NOT FULLY IMPLEMENTED.        **
!**  BY FORCING UNSELECT BEFORE SELECT, ONE CAN SAFEGUARD THAT A SUBROUTINE   **
!**  CHANGES THE SETTING OF A PARENT ROUTINE.                                 **
!**                                                                           **
!**  RBOX HAS THE FOLLOWING FUNCTIONS:                                        **
!**    PSG2 AND PSG4 IS EVALUATED WITHIN A SPHERE WITH RADIUS RBOX            **
!**                                                                           **
!**  THE RADIAL GRID DEFINES THE BOUNDARY CONDITIONS FOR THE ATOMIC           **
!**    CALCULATIONS: R(NR-3) IS THE POSITION OF THE OUTERMOST ZERO.           **
!**                                                                           **
!**                                                                           **
!**                                              P.E. BLOECHL, (1991-2010)    **
!*******************************************************************************
!== PARAMETERS DEFINING THE SETUP CONSTRUCTION =================================
TYPE SETUPPARMS_TYPE
  CHARACTER(128)      :: ID 
  REAL(8)             :: POW_POT=0.D0
  LOGICAL(4)          :: TVAL0_POT       
  REAL(8)             :: VAL0_POT
  REAL(8)             :: RC_POT
  REAL(8)             :: POW_CORE
  LOGICAL(4)          :: TVAL0_CORE       
  REAL(8)             :: VAL0_CORE
  REAL(8)             :: RC_CORE
  CHARACTER(32)       :: TYPE
  REAL(8),ALLOCATABLE :: RCL(:)    !(LX+1)
  REAL(8),ALLOCATABLE :: LAMBDA(:) !(LX+1)
END TYPE SETUPPARMS_TYPE
TYPE ATOMWAVES_TYPE
  INTEGER(4)             :: NB=-1
  INTEGER(4)             :: NC=-1
  INTEGER(4),ALLOCATABLE :: LOFI(:)
  INTEGER(4),ALLOCATABLE :: SOFI(:)
  INTEGER(4),ALLOCATABLE :: NNOFI(:)
  REAL(8)   ,ALLOCATABLE :: EOFI(:)
  REAL(8)   ,ALLOCATABLE :: FOFI(:)
  REAL(8)   ,ALLOCATABLE :: AEPSI(:,:)
  REAL(8)   ,ALLOCATABLE :: AEPSISM(:,:)
  REAL(8)   ,ALLOCATABLE :: AEPOT(:)
END TYPE ATOMWAVES_TYPE
TYPE SETTING_TYPE
  LOGICAL  :: TREL ! RELATIVISTIC OR NON-RELATIVISTIC
  LOGICAL  :: SO   ! SPIN ORBIT SPLITTING OR SCALAR RELATIVISTIC
  LOGICAL  :: ZORA ! ZEROTH ORDER RELATIVISTIC CORRECTION
  REAL(8)  :: FOCK ! ADMIXTURE OF EXACT EXCHANGE - MU_X
END TYPE SETTING_TYPE
TYPE THIS_TYPE
INTEGER(4)             :: I            ! INTEGER IDENTIFIER (ISPECIES)
CHARACTER(32)          :: ID           ! IDENTIFIER (SPECIES-NAME)
INTEGER(4)             :: GID          ! GRID ID FOR R-SPACE GRID
INTEGER(4)             :: GIDG         ! GRID ID FOR G-SPACE GRID
REAL(8)                :: AEZ          ! ATOMIC NUMBER
REAL(8)                :: RAD          ! ATOMIC SPHERE RADIUS FOR PDOS ETC.
REAL(8)                :: RCBG
REAL(8)                :: RCSM         ! GAUSSIAN DECAY FOR COMPENSATION CHARGE
INTEGER(4)             :: LX           ! HIGHEST ANGULAR MOMENTUM
INTEGER(4)             :: LNX          ! #(ORBITAL SHELLS)
INTEGER(4)             :: LMNX
INTEGER(4)             :: LMRX         ! #(ANGULAR MOMENTA FOR 1C-DENSITY)
INTEGER(4),ALLOCATABLE :: LOX(:)       !(LNX) MAIN ANGULAR MOMENTA 
INTEGER(4),ALLOCATABLE :: ISCATT(:)    !(LNX) =-1 FOR SEMI-CORE STATE
                                   !      = 0 FOR VALENCE STATE   (PHI)
                                   !   = 1 FOR 1ST SCATTERING STATE (PHIDOT)
REAL(8)   ,ALLOCATABLE :: VADD(:)      !(NR)
REAL(8)   ,ALLOCATABLE :: PSPOT(:)     !(NR)
REAL(8)   ,ALLOCATABLE :: AECORE(:)    !(NR)  CORE ELECTRON DENSITY
REAL(8)   ,ALLOCATABLE :: PSCORE(:)    !(NR)  PSEUDIZED ELECTRON DENSITY
REAL(8)   ,ALLOCATABLE :: PRO(:,:)     !(NR,LNX)  PROJECTOR FUNCTIONS
REAL(8)                :: RBOX         ! PARTIAL WAVES HAVE OUTER NODE AT RBOX
REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)   !(NR,LNX)  AE PARTIAL WAVES
REAL(8)   ,ALLOCATABLE :: AEPHISM(:,:) !(NR,LNX)  
REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)   !(NR,LNX)  PS PARTIAL WAVES
REAL(8)   ,ALLOCATABLE :: PSPHISM(:,:) !(NR,LNX)  
REAL(8)   ,ALLOCATABLE :: UPHI(:,:)    !(NR,LNX)  NODELESS PARTIAL WAVES
REAL(8)   ,ALLOCATABLE :: UPHISM(:,:)  !(NR,LNX)  
REAL(8)   ,ALLOCATABLE :: QPHI(:,:)    !(NR,LNX)  REDUCED-NODE PARTIAL WAVES
REAL(8)   ,ALLOCATABLE :: QPHISM(:,:)  !(NR,LNX)
REAL(8)   ,ALLOCATABLE :: QPHIDOT(:,:) !(NR,LNX)  R-N SCATTERING PARTIAL WAVES
REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)!(NR,LNX)  NL SCATTERING PARTIAL WAVES
REAL(8)   ,ALLOCATABLE :: PSPHIDOT(:,:)!(NR,LNX)  PS SCATTERING PARTIAL WAVES
REAL(8)   ,ALLOCATABLE :: PSPHIDOTSM(:,:)!(NR,LNX)  
REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:,:)!(NR,LNX)  AE SCATTERING PARTIAL WAVES
REAL(8)   ,ALLOCATABLE :: AEPHIDOTSM(:,:)!(NR,LNX)  
REAL(8)   ,ALLOCATABLE :: DTKIN(:,:)   !(LNX,LNX) 1C-KIN. EN. MATRIX ELEMENTS
REAL(8)   ,ALLOCATABLE :: DOVER(:,:)   !(LNX,LNX) 1C-OVERLAP MATRIX ELEMENTS
REAL(8)   ,ALLOCATABLE :: DATH(:,:)   !(LNX,LNX) 1C-OVERLAP MATRIX ELEMENTS
REAL(8)   ,ALLOCATABLE :: PROPHIDOT(:,:)  !(LNX,LNX) <PRO|PSPHIDOT>
REAL(8)   ,ALLOCATABLE :: COREVALENCEX(:,:)  !(LNX,LNX) CORE VALENCE EXCHANGE
REAL(8)   ,ALLOCATABLE :: VADDOFG(:)   !(NGX)
REAL(8)   ,ALLOCATABLE :: PSCOREOFG(:) !(NGX)
REAL(8)   ,ALLOCATABLE :: VHATOFG(:)   !(NGX)
REAL(8)   ,ALLOCATABLE :: NHATPRIMEOFG(:)  !(NGX)
REAL(8)   ,ALLOCATABLE :: PROOFG(:,:)  !(NR,LNX)  
LOGICAL(4)             :: LOCORBINI=.FALSE. !LOCAL ORBS ARE INITIALIZED IF TRUE
REAL(8)                :: LOCORBRAD(4)=5.D0  ! RADIUS OF LOCAL ORBITALS
INTEGER(4)             :: LOCORBNOFL(4)=0 ! #(LOCAL ORBITALS PER L)
INTEGER(4)             :: LOCORBLNX=0     ! #(LOCAL ORBITAL-SHELLS)
INTEGER(4),ALLOCATABLE :: LOCORBLOX(:)    ! L FOR EACH LOCAL ORBITAL-SHELL
REAL(8)   ,ALLOCATABLE :: LOCORBAMAT(:,:) ! |CHI>=|PHI>*AMAT
REAL(8)   ,ALLOCATABLE :: LOCORBBMAT(:,:) ! |PSI>=|CHI>BMAT<PTILDE|PSITILDE>
REAL(8)                :: M
REAL(8)                :: ZV
CHARACTER(64)          :: COREID
REAL(8)                :: PSG2
REAL(8)                :: PSG4
CHARACTER(32)          :: SOFTCORETYPE
CHARACTER(16)          :: FILEID
TYPE(SETTING_TYPE)     :: SETTING
TYPE(SETUPPARMS_TYPE)  :: PARMS
TYPE(ATOMWAVES_TYPE)   :: ATOM
TYPE(THIS_TYPE),POINTER:: NEXT
END TYPE THIS_TYPE
!
INTEGER(4)              :: GIDG_PROTO=0 !PROTOTYPE G-GRID
INTEGER(4)              :: NSP=0
LOGICAL,SAVE            :: SELECTED=.FALSE. ! CONTAINER FOR ACTUAL SETTING
INTEGER(4)              :: LMRXX=0
INTEGER(4)              :: LMNXX=0
INTEGER(4)              :: LNXX=0
TYPE(THIS_TYPE),POINTER :: FIRST
TYPE(THIS_TYPE),POINTER :: THIS
TYPE FASTACCESS_TYPE  
  TYPE(THIS_TYPE),POINTER :: THIS
END TYPE FASTACCESS_TYPE
TYPE(FASTACCESS_TYPE),ALLOCATABLE :: FASTACCESS(:)
END MODULE SETUP_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$ISELECT(I)
!     **************************************************************************
!     **  SELECTS A SETUP PER INTEGER INDEX                                   **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: I
      INTEGER(4)            :: J
!     **************************************************************************
!
!     ==========================================================================
!     ==  DEFINE FASTACCESS LOOKUP ARRAY                                      ==
!     ==========================================================================
      IF(.NOT.ALLOCATED(FASTACCESS)) THEN
        IF(.NOT.ASSOCIATED(FIRST)) THEN
          CALL ERROR$MSG('NO SETUPS DEFINED')
          CALL ERROR$STOP('SETUP$ISELECT')
        END IF
        THIS=>FIRST
        NSP=1
        DO WHILE(ASSOCIATED(THIS%NEXT))
          NSP=NSP+1
          THIS=>THIS%NEXT
        ENDDO
        ALLOCATE(FASTACCESS(NSP))
        FASTACCESS(1)%THIS=>FIRST
        DO J=2,NSP
          FASTACCESS(J)%THIS=>FASTACCESS(J-1)%THIS%NEXT 
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  CHECK IF I LIES OUT OF RANGE                                        ==
!     ==========================================================================
      IF(I.GT.NSP.OR.I.LT.0) THEN
        CALL ERROR$MSG('INDEX I OUT OF RANGE')
        CALL ERROR$I4VAL('I',I)
        CALL ERROR$I4VAL('NSP',NSP)
        CALL ERROR$STOP('SETUP$ISELECT')
      END IF
!
!     ==========================================================================
!     ==  SELECT OR UNSELECT                                                  ==
!     ==========================================================================
      IF(I.EQ.0) THEN
!       ==  UNSELECT ===========================================================
        IF(.NOT.SELECTED) THEN
          CALL ERROR$MSG('SAFEGUARD FUNCTION:')
          CALL ERROR$MSG('CANNOT UNSELECT A SETUP THAT IS NOT SELECTED')
          CALL ERROR$STOP('SETUP$ISELECT')
        END IF
        SELECTED=.FALSE.
        NULLIFY(THIS)
      ELSE 
!       ==  SELECT =============================================================
        IF(SELECTED) THEN
          CALL ERROR$MSG('SAFEGUARD FUNCTION:')
          CALL ERROR$MSG('ANOTHER SETUP IS ALREADY SELECTED:')
          CALL ERROR$MSG('UNSELECT BEFORE SELECT')
          CALL ERROR$STOP('SETUP$ISELECT')
        END IF
        THIS=>FASTACCESS(I)%THIS
        SELECTED=.TRUE.
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$UNSELECT()
!     **************************************************************************
!     **  UNSELECTS                                                           **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
!     **************************************************************************
      IF(.NOT.SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT UNSELECT A SETUP THAT IS NOT SELECTED')
        CALL ERROR$STOP('SETUP$UNSELECT')
      END IF
      SELECTED=.FALSE.
      NULLIFY(THIS)
      RETURN
      END
!!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$SELECT(ID)
!     **************************************************************************
!     **  SELECTS A SETUP PER ID                                              **
!     **  AND CREATES A NEW, IF IT DOES NOT EXIST                             **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,PARAMETER  :: GMAX=15.D0 ! <-> EPW OF ABOUT 200 RY
      REAL(8)     ,PARAMETER  :: G1=1.D-3
      INTEGER(4)  ,PARAMETER  :: NG=512
!!$!     ==  THE RADIAL GRID HAS BEEN CHANGED BETWEEN REVISION 857 AND 1079 
!!$!     ==  IN THE DEVEL BRANCH. THIS CAUSED SMALL DIFFERENCES IN THE RESULTS
!!$!     ==  HERE THE PREVIOUS PARAMETERS:
!!$      REAL(8)     ,PARAMETER  :: GMAX=30.D0 ! <-> EPW OF ABOUT 200 RY
!!$      REAL(8)     ,PARAMETER  :: G1=1.175316829807299D-4
!!$      INTEGER(4)  ,PARAMETER  :: NG=250
      REAL(8)                 :: DEX
      CHARACTER(6)            :: TYPEID
!     **************************************************************************
      IF(SELECTED) THEN
        CALL ERROR$MSG('SAFEGUARD FUNCTION:')
        CALL ERROR$MSG('CANNOT SELECT A SETUP WHILE ANOTHER ONE IS SELECTED')
        CALL ERROR$CHVAL('SELECTED ID',THIS%ID)
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$SELECT')
      END IF
      SELECTED=.TRUE.    ! THIS FUNCTION WILL SELECT AN INSTANCE
!
!     ==========================================================================
!     == CHECK IF ALREADY SELECTED                                            ==
!     ==========================================================================
      IF(ASSOCIATED(THIS)) THEN
        IF(THIS%ID.EQ.ID) RETURN
      END IF
!
!     ==========================================================================
!     == CHECK IF PRESENT                                                     ==
!     ==========================================================================
      IF(ASSOCIATED(FIRST)) THEN
        THIS=>FIRST
        DO 
          IF(THIS%ID.EQ.ID) RETURN
          IF(.NOT.ASSOCIATED(THIS%NEXT))EXIT
          THIS=>THIS%NEXT
        ENDDO
        ALLOCATE(THIS%NEXT)
        THIS%NEXT%I=THIS%I+1
        THIS=>THIS%NEXT
      ELSE
        ALLOCATE(FIRST)
        THIS=>FIRST
        THIS%I=1
      END IF
!
!     ==========================================================================
!     == CREATE NEW
!     ==========================================================================
      IF(GIDG_PROTO.EQ.0) THEN
        TYPEID='LOG'
        CALL RADIAL$NEW(TYPEID,GIDG_PROTO)
        DEX=LOG(GMAX/G1)/REAL(NG-1,KIND=8)
        CALL RADIAL$SETI4(GIDG_PROTO,'NR',NG)
        CALL RADIAL$SETR8(GIDG_PROTO,'R1',G1)
        CALL RADIAL$SETR8(GIDG_PROTO,'DEX',DEX)

        !WRITE GRID FOR PROJECTORS TO BANDDATA MODULE
        CALL BANDDATA$SETI4('NG_PROTO',NG)
        CALL BANDDATA$SETCH('TYPEID_PROTO',TYPEID)
        CALL BANDDATA$SETR8('GMAX_PROTO',GMAX)
        CALL BANDDATA$SETR8('G1_PROTO',G1)
        CALL BANDDATA$SETR8('DEX_PROTO',DEX)
      END IF
!
!     ==========================================================================
!     == CREATE NEW
!     ==========================================================================
      IF(ALLOCATED(FASTACCESS)) DEALLOCATE(FASTACCESS)
      THIS%ID    =ID
      THIS%GID   =0
      THIS%GIDG   =0
      THIS%AEZ   =0.D0
      THIS%RAD   =0.D0
      THIS%RCBG  =0.D0
      THIS%RCSM  =0.D0
      THIS%LX    =0
      THIS%LNX   =0
      THIS%LMNX  =0
      THIS%LMRX  =0
!!$      NULLIFY(THIS%LOX)     !(LNX)
!!$      NULLIFY(THIS%ISCATT)  !(LNX)
!!$      NULLIFY(THIS%VADD)    !(NRX)
!!$      NULLIFY(THIS%PSPOT)   !(NRX)
!!$      NULLIFY(THIS%AECORE)  !(NRX)
!!$      NULLIFY(THIS%PSCORE)  !(NRX)
!!$      NULLIFY(THIS%PRO)     !(NRX,LNX)
      THIS%RBOX  =0.D0
!!$      NULLIFY(THIS%AEPHI)   !(NRX,LNX)
!!$      NULLIFY(THIS%AEPHISM) !(NRX,LNX)
!!$      NULLIFY(THIS%PSPHI)   !(NRX,LNX)
!!$      NULLIFY(THIS%PSPHISM) !(NRX,LNX)
!!$      NULLIFY(THIS%UPHI)    !(NRX,LNX)
!!$      NULLIFY(THIS%UPHISM)  !(NRX,LNX)
!!$      NULLIFY(THIS%QPHI)    !(NRX,LNX)
!!$      NULLIFY(THIS%QPHISM)  !(NRX,LNX)
!!$      NULLIFY(THIS%QPHIDOT) !(NRX,LNX)
!!$      NULLIFY(THIS%NLPHIDOT)!(NRX,LNX)
!!$      NULLIFY(THIS%PSPHIDOT)!(NRX,LNX)
!!$      NULLIFY(THIS%PSPHIDOTSM)!(NRX,LNX)
!!$      NULLIFY(THIS%AEPHIDOT)!(NRX,LNX)
!!$      NULLIFY(THIS%AEPHIDOTSM)!(NRX,LNX)
!!$      NULLIFY(THIS%DTKIN)   !(LNXX,LNX)
!!$      NULLIFY(THIS%DOVER)   !(LNXX,LNX)
!!$      NULLIFY(THIS%PROPHIDOT)!(LNX,LNX)
!!$      NULLIFY(THIS%COREVALENCEX)!(LNX,LNX)
!!$      NULLIFY(THIS%VADDOFG) !(NGX)
!!$      NULLIFY(THIS%PSCOREOFG) !(NGX)
!!$      NULLIFY(THIS%VHATOFG) !(NGX)
!!$      NULLIFY(THIS%NHATPRIMEOFG) !(NGX)
!!$      NULLIFY(THIS%PROOFG)  !(NGX,LNX)
      THIS%LOCORBINI=.FALSE. ! INITIALIZED?
      THIS%LOCORBRAD(:)=5.D0  ! RADIUS OF LOCAL ORBITALS
      THIS%LOCORBNOFL(:)=0 ! (4) #(LOCAL ORBITALS PER L)
      THIS%LOCORBLNX=0     ! #(LOCAL ORBITAL-SHELLS)
!!$      NULLIFY(THIS%LOCORBLOX)  !(LOCORBLNX)  L FOR EACH LOCAL ORBITAL-SHELL
!!$      NULLIFY(THIS%LOCORBAMAT) !(LNX,LOCORBLNX) |CHI>=|PHI>*AMAT
!!$      NULLIFY(THIS%LOCORBBMAT) !(LOCORBLNX,LNX) |PSI>=|CHI>BMAT<PTILDE|PSITILDE>
      THIS%M     =0.D0
      THIS%ZV     =0.D0
      THIS%PSG2  =0.D0
      THIS%PSG4  =0.D0
      THIS%SOFTCORETYPE='NONE'
      THIS%SETTING%TREL=.TRUE.
      THIS%SETTING%SO=.FALSE.
      THIS%SETTING%ZORA=.FALSE.
      THIS%SETTING%FOCK=0.D0
      THIS%PARMS%ID=' '
      THIS%PARMS%POW_POT=0.D0
      THIS%PARMS%TVAL0_POT=.FALSE.
      THIS%PARMS%VAL0_POT=0.D0
      THIS%PARMS%RC_POT=0.D0
       
!      THIS%ATOM%
      WRITE(THIS%FILEID,*)THIS%I
      THIS%FILEID='ATOM'//TRIM(ADJUSTL(THIS%FILEID))
      NULLIFY(THIS%NEXT)
!
!     == INITIALIZE DEFAULT VALUES
      THIS%RCBG  =1.D0/SQRT(0.218D0)
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETCH(ID,VAL)
!     **************************************************************************
!     **  COLLECTS INTERNAL DATA                                              **
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      CHARACTER(*),INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'ID') THEN
        VAL=THIS%ID
      ELSE IF(ID.EQ.'SOFTCORETYPE') THEN
        VAL=THIS%SOFTCORETYPE
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('SETUP$GETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETI4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'NSP') THEN
        VAL=NSP
      ELSE IF(ID.EQ.'ISP') THEN
        VAL=THIS%I
      ELSE IF(ID.EQ.'LNXX') THEN
        VAL=LNXX
      ELSE IF(ID.EQ.'LMNXX') THEN
        VAL=LMNXX
      ELSE IF(ID.EQ.'LMRXX') THEN
        VAL=LMRXX
      ELSE IF(ID.EQ.'GID') THEN
        VAL=THIS%GID
      ELSE IF(ID.EQ.'GIDG') THEN
        VAL=THIS%GIDG
      ELSE IF(ID.EQ.'NR') THEN
        CALL RADIAL$GETI4(THIS%GID,'NR',VAL)
      ELSE IF(ID.EQ.'LNX') THEN
        VAL=THIS%LNX
      ELSE IF(ID.EQ.'LMNX') THEN
        VAL=THIS%LMNX
      ELSE IF(ID.EQ.'LMRX') THEN
        VAL=THIS%LMRX
      ELSE IF(ID.EQ.'NC') THEN
        VAL=THIS%ATOM%NC
      ELSE IF(ID.EQ.'NB') THEN
        VAL=THIS%ATOM%NB
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$SETI4A(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      INTEGER(4)  ,INTENT(IN)  :: VAL(LEN)
      INTEGER(4)               :: I
!     **************************************************************************
      IF(ID.EQ.'NOFLCHI') THEN
        IF(LEN.NE.4) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$SETI4A')
        END IF
!       == CHECK IF VALUE HAS CHANGED ==========================================
        DO I=1,4
          THIS%LOCORBINI=THIS%LOCORBINI.AND.(THIS%LOCORBNOFL(I).EQ.VAL(I))
        ENDDO
        THIS%LOCORBNOFL=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$SETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETI4A(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      INTEGER(4)  ,INTENT(OUT) :: VAL(LEN)
!     **************************************************************************
!     ==========================================================================
!     == ANGULAR MOMENTA OF THE PARTIAL WAVES                                 ==
!     ==========================================================================
      IF(ID.EQ.'LOX') THEN
        IF(LEN.NE.THIS%LNX) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        VAL=THIS%LOX
!
!     ==========================================================================
!     == COUNTS PARTIAL WAVES RELATIVE TO THE HIGHEST VALENCE STATE           ==
!     ==========================================================================
      ELSE IF(ID.EQ.'ISCATT') THEN
        IF((LEN.NE.THIS%LNX).OR.(THIS%LNX.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
!!$        IF(.NOT.ASSOCIATED(THIS%ISCATT)) THEN
        IF(.NOT.ALLOCATED(THIS%ISCATT)) THEN
          CALL ERROR$MSG('DATA NOT AVAILABLE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        VAL=THIS%ISCATT
!
!     ==========================================================================
!     == MAIN ANGULAR MOMENTA OF ATOMIC STATES                                ==
!     ==========================================================================
      ELSE IF(ID.EQ.'LB') THEN
        IF(LEN.NE.THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        VAL=THIS%ATOM%LOFI
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$SETL4A(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      LOGICAL(4)  ,INTENT(IN)  :: VAL(LEN)
!     **************************************************************************
      IF(ID.EQ.'TORB') THEN
        CALL ERROR$MSG('ID="TORB" IS OBSOLETE')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$SETL4A')
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$SETL4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETL4A(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      LOGICAL(4)  ,INTENT(OUT)  :: VAL(LEN)
!     **************************************************************************
      VAL(:)=.FALSE.
      IF(ID.EQ.'TORB') THEN
        CALL ERROR$MSG('ID="TORB" IS OBSOLETE')
        CALL ERROR$STOP('SETUP$GETL4A')
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETL4A')
      END IF
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE SETUP$GETR8(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'AEZ') THEN
        VAL=THIS%AEZ
      ELSE IF(ID.EQ.'RCSM') THEN
        VAL=THIS%RCSM
      ELSE IF(ID.EQ.'RCBG') THEN
        VAL=THIS%RCBG
      ELSE IF(ID.EQ.'ZV') THEN
        VAL=THIS%ZV
      ELSE IF(ID.EQ.'M') THEN
        VAL=THIS%M
      ELSE IF(ID.EQ.'PS<G2>') THEN
        VAL=THIS%PSG2
      ELSE IF(ID.EQ.'PS<G4>') THEN
        VAL=THIS%PSG4
      ELSE IF(ID.EQ.'RBOX') THEN
        VAL=THIS%RBOX  ! USED IN PAW_OPTEELS.F90
      ELSE IF(ID.EQ.'RAD') THEN
        VAL=THIS%RAD  ! ATOM-RADIUS FOR PDOS ETC.
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETR8A(ID,LEN,VAL)
!     **************************************************************************
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      REAL(8)     ,INTENT(OUT) :: VAL(LEN)
      INTEGER(4)               :: NR,NG
      INTEGER(4)               :: I
!     **************************************************************************
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
      CALL RADIAL$GETI4(THIS%GIDG,'NR',NG)
!
!     ==========================================================================
!     == PROJECTOR FUNCTIONS                                                  ==
!     ==========================================================================
      IF(ID.EQ.'PRO') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PRO,(/LEN/))
      ELSE IF(ID.EQ.'PROOFG') THEN
        IF(LEN.NE.THIS%LNX*NG) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PROOFG,(/LEN/))
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%AEPHI,(/LEN/))
!
!     ==========================================================================
!     ==  PSEUDO PARTIAL WAVES                                                ==
!     ==========================================================================
      ELSE IF(ID.EQ.'PSPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PSPHI,(/LEN/))
!
!     ==========================================================================
!     ==  NODELESS PARTIAL WAVES                                              ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NLPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%UPHI,(/LEN/))
!
!     ==========================================================================
!     ==  REDUCED NODE PARTIAL WAVES                                          ==
!     ==========================================================================
      ELSE IF(ID.EQ.'QPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%QPHI,(/LEN/))
!
!     ==========================================================================
!     ==  NODELESS SCATTERING PARTIAL WAVES                                   ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NLPHIDOT') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%NLPHIDOT,(/LEN/))
!
!     ==========================================================================
!     ==  REDUCED-NODE SCATTERING PARTIAL WAVES                               ==
!     ==========================================================================
      ELSE IF(ID.EQ.'QPHIDOT') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%QPHIDOT,(/LEN/))
!
!     ==========================================================================
!     ==  SCATTERING PSEUDO PARTIAL WAVES                                     ==
!     ==========================================================================
      ELSE IF(ID.EQ.'PSPHIDOT') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PSPHIDOT,(/LEN/))
!
!     ==========================================================================
!     ==  SCATTERING ALL-ELECTRON PARTIAL WAVES                               ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPHIDOT') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%AEPHIDOT,(/LEN/))
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AECORE') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%AECORE
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'PSCORE') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%PSCORE
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'VADD') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%VADD
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'DEKIN') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%DTKIN,(/LEN/))
!
!     ==========================================================================
!     ==  OVERLAP DIFFERENCE MATRIX ELEMENTS <AEPSI|AEPSI>-<PSPSI|PSPSI>      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'DO') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%DOVER,(/LEN/))
!
!     ==========================================================================
!     ==  MATRIX ELEMENTS FOR CORE VALENCE EXCHANGE WITH PARTIAL WAVES        ==
!     ==========================================================================
      ELSE IF(ID.EQ.'CVX') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%COREVALENCEX,(/LEN/))
!
!     ==========================================================================
!     ==  <PRO|PSPHIDOT>                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'PROPHIDOT') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PROPHIDOT,(/LEN/))
!
!     ==========================================================================
!     ==  ATOMIC WAVE FUNCTIONS                                               ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPSI') THEN
        IF(LEN.NE.NR*THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%ATOM%AEPSI,(/LEN/))
!
!     ==========================================================================
!     ==  ATOMIC OCCUPATIONS                                                  ==
!     ==========================================================================
      ELSE IF(ID.EQ.'FOFI') THEN
        IF(LEN.NE.THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%ATOM%FOFI
!
!     ==========================================================================
!     ==  ATOMIC ONE-PARTICLE ENERGIES                                        ==
!     ==========================================================================
      ELSE IF(ID.EQ.'EOFI') THEN
        IF(LEN.NE.THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%ATOM%EOFI
!
!     ==========================================================================
!     ==  OCCUPATION OF ATOMIC ORBITALS                                       ==
!     ==========================================================================
      ELSE IF(ID.EQ.'FB') THEN
        IF(LEN.NE.THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%ATOM%FOFI
!
!     ==========================================================================
!     ==  ATOMIC ENERGY LEVELS                                                ==
!     ==========================================================================
      ELSE IF(ID.EQ.'EB') THEN
        IF(LEN.NE.THIS%ATOM%NB) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%ATOM%EOFI
!
!     ==========================================================================
!     ==  ATOMIC DENSITY                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AERHO') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL(:)=0.D0
        DO I=1,THIS%ATOM%NB
          VAL(:)=VAL(:)+THIS%ATOM%FOFI(I) &
     &                 *(THIS%ATOM%AEPSI(:,I)**2+THIS%ATOM%AEPSISM(:,I)**2)
        ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AEPOT') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL(:)=THIS%ATOM%AEPOT(:)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NUCPOT') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        CALL RADIAL$NUCPOT(THIS%GID,NR,THIS%AEZ,VAL)
!
!     ==========================================================================
!     ==  WRONG ID                                                            ==
!     ==========================================================================
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETR8A')
      END IF
      RETURN
      END  
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$GETFOFG(ID,TDER,IND,NG_,G2,CELLVOL,F)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID      !IDENTIFIER
      LOGICAL(4)  ,INTENT(IN)  :: TDER    ! CALCULATE RADIAL DERVATIVE
      INTEGER(4)  ,INTENT(IN)  :: IND     ! SELECTOR (USED ONLY FOR ID=PRO)
      INTEGER(4)  ,INTENT(IN)  :: NG_     ! #(PLANE WAVES)
      REAL(8)     ,INTENT(IN)  :: G2(NG_) ! G**2
      REAL(8)     ,INTENT(OUT) :: F(NG_)  
      REAL(8)     ,ALLOCATABLE :: FOFG(:)  !(NG)
      REAL(8)     ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)               :: IG
      INTEGER(4)               :: NG
      INTEGER(4)               :: NGAMMA
      REAL(8)                  :: CELLVOL
      REAL(8)                  :: G
      INTEGER(4)               :: GIDG
!     ******************************************************************
      IF(NG_.EQ.0) RETURN
      CALL RADIAL$GETI4(THIS%GIDG,'NR',NG)
      ALLOCATE(FOFG(NG))
      IF(ID.EQ.'PRO') THEN
        IF(IND.LT.1.OR.IND.GT.THIS%LNX) THEN
          CALL ERROR$MSG('LN OUT OF RANGE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('IND',IND)
          CALL ERROR$I4VAL('LNX',THIS%LNX)
          CALL ERROR$STOP('SETUP$GETFOFG')
        END IF
        FOFG(:)=THIS%PROOFG(:,IND)
      ELSE IF(ID.EQ.'PSCORE') THEN
        FOFG(:)=THIS%PSCOREOFG(:)
      ELSE IF(ID.EQ.'VADD') THEN
        FOFG(:)=THIS%VADDOFG(:)
      ELSE IF(ID.EQ.'V0') THEN
        FOFG(:)=THIS%VHATOFG(:)
      ELSE IF(ID.EQ.'G0') THEN
        FOFG(:)=THIS%NHATPRIMEOFG(:)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETFOFG')
      END IF
!
!     ==================================================================
!     == INTERPOLATE VALUES FROM RADIAL GRID
!     ==================================================================
      GIDG=THIS%GIDG
      NGAMMA=0
      IF(TDER) THEN
        G=SQRT(G2(1))
        IF(G.LT.1.D-6) NGAMMA=1
        CALL RADIAL$DERIVATIVE(GIDG,NG,FOFG,G,F(1))
        F(1)=G*F(1)
        DO IG=2,NG_
          IF(ABS(G2(IG)-G2(IG-1)).LT.1.D-6) THEN
            F(IG) =F(IG-1)
          ELSE
            G=SQRT(G2(IG))
            IF(G.LT.1.D-6) NGAMMA=IG
            CALL RADIAL$DERIVATIVE(GIDG,NG,FOFG,G,F(IG))
            F(IG)=G*F(IG)
          END IF
        ENDDO
      ELSE
        G=SQRT(G2(1))
        IF(G.LT.1.D-6) NGAMMA=1
        CALL RADIAL$VALUE(GIDG,NG,FOFG,G,F(1))
        DO IG=2,NG_
          IF(ABS(G2(IG)-G2(IG-1)).LT.1.D-6) THEN
            F(IG) =F(IG-1)
          ELSE
            G=SQRT(G2(IG))
            IF(G.LT.1.D-6) NGAMMA=IG
            CALL RADIAL$VALUE(GIDG,NG,FOFG,G,F(IG))
          END IF
        ENDDO
      END IF
      DEALLOCATE(FOFG)
!
!     ==================================================================
!     == CORRECT EXTRAPOLATION TO THE GAMMA POINT                     ==
!     ==================================================================
      IF(NGAMMA.NE.0) THEN
        IF(TDER) THEN
          NGAMMA=0
        ELSE
          IF(ID.EQ.'G0') THEN 
            F(NGAMMA)=4.D0*PI
          ELSE IF(ID.EQ.'V0') THEN
            F(NGAMMA)=PI*(THIS%RCBG**2-THIS%RCSM**2)*4.D0*PI
          END IF
        END IF
      END IF
!
!     ==================================================================
!     == DIVIDE BY CELLVOL                                            ==
!     ==================================================================
      F=F/CELLVOL
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE SETUP_QNJOFE(NR,L,N,J,QNJOFE)
!!$!     **************************************************************************
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: NR
!!$      INTEGER(4),INTENT(IN) :: L
!!$      INTEGER(4),INTENT(IN) :: N
!!$      INTEGER(4),INTENT(IN) :: J
!!$      REAL(8)   ,INTENT(OUT):: QNJOFE
!!$!     ***********************************************************************
!!$
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$READSTRCIN(LL_STRC_)
!     **************************************************************************
!     ** ANALYZES THE !STRUCTURE!SPECIES BLOCK OF THE STRUCTURE INPUT FILE.   **
!     ** IT IS CALLED DIRECTLY FROM STRCIN_SPECIES OF PAW_IOROUTINES.F90.     **
!     **                                                                      **
!     **  SETUP$RESOLVE_SETUPID                                               **
!     **     IF(!SPECIED:ID PRESENT) THEN                                     **
!     **       SETUP_BUILDPARMSON_NDLSS/HBS/7560                              **
!     **       ATOMLIB$SCNTLLOOKUPONE                                         **
!     **     END IF                                                           **
!     **     SETUP_LOOKUPSETUP  (READ !STRUCTURE!SPECIES!AUGMENT)             **
!     **       SETUP_LOOKUPGENERIC                                            **
!     **       SETUP_LOOKUP_KERKER/HBS/NDLSS                                  **
!     **       SETUP_LOOKUP_GRID                                              **
!     **       SETUP_LOOKUP_PSPOT                                             **
!     **       SETUP_LOOKUP_PSCORE                                            **
!     **                                                                      **
!     ********************************************PETER BLOECHL, GOSLAR 2014****
      USE PERIODICTABLE_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      USE SETUP_MODULE, ONLY: NSP,THIS
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN):: LL_STRC_
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      TYPE(LL_TYPE)           :: LL_STRC
      INTEGER(4)   ,PARAMETER :: LX=3
      REAL(8)                 :: AEZ     ! ATOMIC NUMBER
      REAL(8)                 :: ZV      ! #(VALENCE ELECTRONS)
      CHARACTER(64)           :: COREID  ! IDENTIFIER FOR THE FROZEN CORE
      LOGICAL(4)              :: TSO     ! SPIN-ORBIT SWITCH
      REAL(8)                 :: RBOX
      CHARACTER(32)           :: TYPE    ! SETUP TYPE
      CHARACTER(32)           :: SPNAME  ! SPECIES NAME
      REAL(8)                 :: RCL(LX+1)
      REAL(8)                 :: LAMBDA(LX+1)
      REAL(8)                 :: RCSM
      REAL(8)                 :: POW_POT
      REAL(8)                 :: RC_POT
      LOGICAL(4)              :: TVAL0_POT
      REAL(8)                 :: VAL0_POT
      REAL(8)                 :: POW_CORE
      REAL(8)                 :: RC_CORE
      LOGICAL(4)              :: TVAL0_CORE
      REAL(8)                 :: VAL0_CORE
      REAL(8)                 :: DMIN
      REAL(8)                 :: DMAX
      REAL(8)                 :: RMAX
      REAL(8)                 :: R1
      REAL(8)                 :: DEX
      INTEGER(4)              :: NR
      INTEGER(4)              :: GID
      LOGICAL(4)              :: TCHK,TCHK1
      CHARACTER(128)          :: ID
      REAL(8)                 :: RCOV
      REAL(8)                 :: PROTONMASS
      REAL(8)                 :: EV
      REAL(8)                 :: SVAR
      INTEGER(4),ALLOCATABLE  :: NPRO(:)
      INTEGER(4)              :: IND,L,ISP
      INTEGER(4)              :: LENG
      INTEGER(4)              :: LRHOX
!     **************************************************************************
                                         CALL TRACE$PUSH('SETUP$READSTRCIN')
      LL_STRC=LL_STRC_  !AVOID REPOSITIONING OF THE POINTER
      CALL CONSTANTS('U',PROTONMASS)
      CALL CONSTANTS('EV',EV)
!
!     ==========================================================================
!     == RESOLVE SETUPID AND EXPAND SPECIES BY A SETUP BLOCK IN THE LINKEDLIST==
!     ==                                                                      ==
!     == THE KEY IS COMPARED TO THE INTERNAL SETUPS AND THEN TO THE ONES ON   ==
!     == A SPECIFIED PARAMETER FILE.                                          ==
!     ==                                                                      ==
!     ==  THE KEYWORD "ID" BECOMES MEANINGLESS AND EVERY SPECIES              ==
!     ==  WILL HAVE A SUB BLOCK !AUGMENT WHEN FINISHED                        ==
!     ==========================================================================
      CALL SETUP_RESOLVESETUPID(LL_STRC)
!
!     ==========================================================================
!     == NOW ANALYZE THE SPECIES BLOCK
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~')
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE')
      CALL LINKEDLIST$NLISTS(LL_STRC,'SPECIES',NSP)
      DO ISP=1,NSP
        CALL LINKEDLIST$SELECT(LL_STRC,'SPECIES',ISP)
!
!       ========================================================================
!       ==  SPECIES NAME:  CREATE AND SELECT NEW SETUP INSTANCE               ==
!       ========================================================================
        CALL LINKEDLIST$GET(LL_STRC,'NAME',1,SPNAME)
        CALL SETUP$SELECT(SPNAME)
        THIS%ID=SPNAME
        THIS%I=ISP
!
!       ========================================================================
!       ==  #(PARTIAL WAVES PER MAIN ANGULAR MOMENTUM)                        ==
!       ========================================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'NPRO',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('VARIABLE !STRUCTURE!SPECIES:NPRO IS MANDATORY')
          CALL ERROR$STOP('STRCIN_SPECIES')
        END IF
        CALL LINKEDLIST$SIZE(LL_STRC,'NPRO',1,LENG)
        ALLOCATE(NPRO(LENG))
        CALL LINKEDLIST$GET(LL_STRC,'NPRO',1,NPRO)
!       __ DEFINE LNX___________________________________________________________
        THIS%LNX=SUM(NPRO)
!       __ DEFINE LOX___________________________________________________________
        ALLOCATE(THIS%LOX(THIS%LNX))
        IND=0
        DO L=0,LENG-1
          THIS%LOX(IND+1:IND+NPRO(L+1))=L
          IND=IND+NPRO(L+1)
        ENDDO
        DEALLOCATE(NPRO)
!       __LMNX__________________________________________________________________
        THIS%LMNX=SUM(2*THIS%LOX(:)+1)
!
!       ========================================================================
!       ==  MAX. #(ANGULAR MOMENTA) FOR ONE-CENTER DENSITY                    ==
!       ========================================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'LRHOX',1,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$GET(LL_STRC,'LRHOX',1,LRHOX)
          LRHOX=MIN(2*MAXVAL(THIS%LOX),LRHOX)
        ELSE
          LRHOX=2*MAXVAL(THIS%LOX)
        END IF
        THIS%LMRX=(LRHOX+1)**2
!
!       ========================================================================
!       == IDENTIFY SETUP ID                                                  ==
!       ========================================================================
        CALL SETUP_LOOKUPSETUP(LL_STRC,ID,AEZ,ZV,COREID,TSO,RBOX,LX &
     &                             ,TYPE,RCL,LAMBDA &
     &                             ,RCSM,POW_POT,RC_POT,TVAL0_POT,VAL0_POT &
     &                             ,POW_CORE,RC_CORE,TVAL0_CORE,VAL0_CORE &
     &                             ,DMIN,DMAX,RMAX)
        THIS%PARMS%ID=ID       
        THIS%AEZ=AEZ        ! ATOMIC NUMBER
        THIS%RCSM=RCSM      ! DECAY RADIUS FOR COMPENSATION CHARGE
        THIS%RBOX=RBOX    
        THIS%ZV=ZV          ! #(VALENCE ELECTRONS)
        THIS%COREID=COREID  ! IDENTIFIER OF THE FROZEN CORE
        THIS%SETTING%SO=TSO ! SPIN-ORBIT SWITCH
!       __ PARTIAL WAVES________________________________________________________
        THIS%PARMS%TYPE     =TYPE         ! PARTIAL WAVE PSEUDIZATION METHOD
        ALLOCATE(THIS%PARMS%RCL(LX+1))    
        ALLOCATE(THIS%PARMS%LAMBDA(LX+1))
        THIS%PARMS%RCL      =RCL          ! PARTIAL WAVE RADIUS
        THIS%PARMS%LAMBDA   =LAMBDA
!       __PSEUDO POTENTIAL______________________________________________________
        THIS%PARMS%POW_POT  =POW_POT
        THIS%PARMS%RC_POT   =RC_POT
        THIS%PARMS%TVAL0_POT=TVAL0_POT
        THIS%PARMS%VAL0_POT =VAL0_POT
!       __PSEUDO POTENTIAL______________________________________________________
        THIS%PARMS%POW_CORE  =POW_CORE
        THIS%PARMS%RC_CORE   =RC_CORE
        THIS%PARMS%TVAL0_CORE=TVAL0_CORE
        THIS%PARMS%VAL0_CORE =VAL0_CORE
!       __ RADIAL GRID__________________________________________________________
        CALL RADIAL$GRIDPARAMETERS(DMIN,DMAX,RMAX,R1,DEX,NR)
        CALL RADIAL$NEW('SHLOG',GID)
        CALL RADIAL$SETR8(GID,'R1',R1)
        CALL RADIAL$SETR8(GID,'DEX',DEX)
        CALL RADIAL$SETI4(GID,'NR',NR)
        THIS%GID=GID
!
!       ========================================================================
!       ==  ATOMIC MASS                                                       ==
!       ========================================================================
        CALL LINKEDLIST$EXISTD(LL_STRC,'M',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL PERIODICTABLE$GET(AEZ,'MASS',SVAR)
          CALL LINKEDLIST$SET(LL_STRC,'M',0,SVAR/PROTONMASS)
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'M',1,SVAR)
        THIS%M=SVAR*PROTONMASS
!
!       ========================================================================
!       ==  ATOMIC RADIUS FOR PDOS ETC.                                       ==
!       ========================================================================
        CALL PERIODICTABLE$GET(AEZ,'R(ASA)',THIS%RAD)  !SET DEFAULT
        CALL LINKEDLIST$EXISTD(LL_STRC,'RAD/RCOV',1,TCHK)
        CALL LINKEDLIST$EXISTD(LL_STRC,'RAD',1,TCHK1)
        IF(TCHK) THEN
          CALL PERIODICTABLE$GET(AEZ,'R(COV)',SVAR)
          CALL LINKEDLIST$GET(LL_STRC,'RAD/RCOV',0,THIS%RAD)
          THIS%RAD=THIS%RAD*SVAR
        ELSE IF(TCHK1) THEN
          CALL LINKEDLIST$GET(LL_STRC,'RAD',0,THIS%RAD)
        ENDIF
!!$        IF(THIS%RAD.GT.THIS%RBOX) THEN
!!$          CALL ERROR$MSG('VARIABLE "!STRUCTURE!SPECIES:RAD" MUST BE SMALLER')
!!$          CALL ERROR$MSG('THAN VARIABLE "!STRUCTURE!SPECIES!AUGMENT:RBOX".')
!!$          CALL ERROR$R8VAL('RAD[ABOHR]',THIS%RAD)
!!$          CALL ERROR$R8VAL('RBOX[ABOHR]',THIS%RBOX)
!!$          CALL ERROR$STOP('SETUP$READSTRCIN')
!!$        END IF
!
        CALL SETUP$UNSELECT()
        CALL LINKEDLIST$SELECT(LL_STRC,'..')
      ENDDO   ! END OF LOOP OVER SPEECIES
                                         CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_RESOLVESETUPID(LL_STRC)
!     **************************************************************************
!     ** SCANS THE LINKEDLIST FROM THE STRUCTURE INPUT FILE                   **
!     ** AND INSERTS A SETUP BLOCK FOR EVERY KNOWN SETUP ID.                  **
!     ** IT DOES NOT EXPAND SETUPS FOR WHICH THE ID IS UNKNOWN.               **
!     ** IT SHOULD BE EXECUTED BEFORE THE SETUP PARAMTERS ARE USED.           **
!     **                                                                      **
!     ** THERE ARE THREE ALLOWED CHANNELS TO PROVIDE SETUP PARAMETERS:        **
!     ** (1) PARAMETER SET SPECIFIED IN THE SPECIES BLOCK OF THE STRC FILE    **
!     ** (2) INTERNAL PARAMETER SETS FROM SETUP_BUILDPARMSONE_...             **
!     ** (3) STP.CNTL FILE                                                    **
!     **                                                                      **
!     ** THIS ROUTINE IS CALLED BY STRCIN PAW_IOROUTINES.F90 BEFORE           **
!     ** STRUCTURE INPUT FILE IS SCANNED (BEFORE STRCIN_SPECIES)              **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2014******************
      USE PERIODICTABLE_MODULE
      USE CONSTANTS_MODULE
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STRC
      INTEGER(4)   ,PARAMETER  :: LX=3
      CHARACTER(64)            :: ID
      CHARACTER(64)            :: ID1
      INTEGER(4)               :: NSP
      INTEGER(4)               :: ISP
      LOGICAL(4)               :: TCHK
      INTEGER(4)               :: I
      INTEGER(4)               :: NFIL
      CHARACTER(64)            :: STPTYPE
      CHARACTER(2)             :: EL
      REAL(8)                  :: AEZ    ! ATOMIC NUMBER
      REAL(8)                  :: ZV     ! NUMBER OF VALENCE ELECTRONS
      CHARACTER(64)            :: COREID ! SPECIFIES ATOMCORE E.G. 'KR-D+F'
      LOGICAL(4)               :: TSO    ! SPIN-ORBIT SWITCH
      CHARACTER(64)            :: TYPE
      REAL(8)                  :: RBOX
      REAL(8)                  :: RCOV   !COVALENT RADIUS
      REAL(8)                  :: RCSM
      REAL(8)                  :: RCL(LX+1)
      REAL(8)                  :: LAMBDA(LX+1)
      REAL(8)                  :: POTPOW,POTVAL,POTRC
      REAL(8)                  :: CORPOW,CORVAL,CORRC
      LOGICAL(4)               :: TCORVAL,TPOTVAL
      REAL(8)                  :: DMIN,DMAX,RMAX
      LOGICAL(4)               :: TFOUND
      TYPE(LL_TYPE)            :: LL_SCNTL
      INTEGER(4)               :: NENTRY
      LOGICAL(4)               :: TSTPFILE
!     **************************************************************************
                                    CALL TRACE$PUSH('SETUP$RESOLVESETUPID')
!     ==========================================================================
!     == READ SETUP PARAMETER FILE                                            ==
!     ==========================================================================
      CALL FILEHANDLER$ATTACHED('AUGPARMS',TSTPFILE)
      IF(TSTPFILE) THEN
        CALL LINKEDLIST$NEW(LL_SCNTL)
        CALL FILEHANDLER$UNIT('AUGPARMS',NFIL)
        CALL LINKEDLIST$READ(LL_SCNTL,NFIL,'~')
PRINT*,'SETUP PARAMETER FILE READ'
      END IF

!     ==========================================================================
!     == READ SETUP PARAMETER FILE                                            ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STRC,'~',0)
      CALL LINKEDLIST$SELECT(LL_STRC,'STRUCTURE',0)
      CALL LINKEDLIST$NLISTS(LL_STRC,'SPECIES',NSP)
      DO ISP=1,NSP
        CALL LINKEDLIST$SELECT(LL_STRC,'SPECIES',ISP)
!
!       == SETUP BLOCK HAS HIGHEST PRIORITY ====================================
!       == HENCE DO NOTHING IF BLOCK !AUGMENT EXISTS ===========================
        CALL LINKEDLIST$EXISTL(LL_STRC,'AUGMENT',0,TCHK)
        IF(TCHK) THEN
          CALL LINKEDLIST$SELECT(LL_STRC,'..',0)
          CYCLE
        END IF
!
!       == RESOLVE SETUP ID FROM EITHER INTERNAL SET OF STP.CNTL FILE ==========
        CALL LINKEDLIST$EXISTD(LL_STRC,'ID',0,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('NO SETUP PARAMETER AVAILABLE')
          CALL ERROR$MSG('!SPECIES:ID AND !SPECIES!AUGMENT ARE MISSING')
          CALL ERROR$STOP('SETUP$RESOLVESETUPID')
        END IF
        CALL LINKEDLIST$GET(LL_STRC,'ID',0,ID)  
!    
!       ========================================================================
!       == ID HAS THE STRUCTURE EL_TYPE                                       ==
!       ========================================================================
        TFOUND=.TRUE.
        I=INDEX(ID,'_')
        STPTYPE=ID(I+1:)
        EL=ID(1:2)
!    
!       ========================================================================
!       == PARSE INTERNAL SETUPS                                              ==
!       ========================================================================
        TSO=.FALSE.   !SET DEFAULT FOR SPIN-ORBIT SWITCH
        IF(STPTYPE.EQ.'NDLSS_V0') THEN
          CALL SETUP_BUILDPARMSONE_NDLSS(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
     &              ,RCL &
     &              ,POTPOW,POTRC,TPOTVAL,POTVAL &
     &              ,CORPOW,CORRC,TCORVAL,CORVAL &
     &              ,DMIN,DMAX,RMAX)
          LAMBDA=0.D0   !NOT USED
        ELSE IF(STPTYPE.EQ.'NDLSS_SC_V0') THEN
          CALL SETUP_BUILDPARMSONE_NDLSS(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
     &              ,RCL &
     &              ,POTPOW,POTRC,TPOTVAL,POTVAL &
     &              ,CORPOW,CORRC,TCORVAL,CORVAL &
     &              ,DMIN,DMAX,RMAX)
          LAMBDA=0.D0   !NOT USED
        ELSE IF(STPTYPE.EQ.'HBS') THEN
          CALL SETUP_BUILDPARMSONE_HBS(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
     &              ,RCL,LAMBDA &
     &              ,POTPOW,POTRC,TPOTVAL,POTVAL &
     &              ,CORPOW,CORRC,TCORVAL,CORVAL &
     &              ,DMIN,DMAX,RMAX)
        ELSE IF(STPTYPE.EQ.'HBS_SC') THEN
          CALL SETUP_BUILDPARMSONE_HBS(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
     &              ,RCL,LAMBDA &
     &              ,POTPOW,POTRC,TPOTVAL,POTVAL &
     &              ,CORPOW,CORRC,TCORVAL,CORVAL &
     &              ,DMIN,DMAX,RMAX)
        ELSE IF(STPTYPE.EQ.'.75_6.0') THEN
          CALL SETUP_BUILDPARMSONE_7560(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
     &              ,RCL,LAMBDA &
     &              ,POTPOW,POTRC,TPOTVAL,POTVAL &
     &              ,CORPOW,CORRC,TCORVAL,CORVAL &
     &              ,DMIN,DMAX,RMAX)
        ELSE
          TFOUND=.FALSE.
          LAMBDA=0.D0   !INITIALIZE SOMEHOW
        END IF
!    
!       ========================================================================
!       == PARSE SETUP-PARAMETER FILE                                         ==
!       ========================================================================
        IF(TSTPFILE) THEN
          CALL LINKEDLIST$SELECT(LL_SCNTL,'~',0)
          CALL LINKEDLIST$SELECT(LL_SCNTL,'ACNTL',0)
          CALL LINKEDLIST$NLISTS(LL_SCNTL,'AUGMENT',NENTRY)
          DO I=1,NENTRY
            CALL LINKEDLIST$SELECT(LL_SCNTL,'AUGMENT',I)
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'ID',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!ACNTL!AUGMENT:ID NOT SPECIFIED')
              CALL ERROR$STOP('SETUP$RESOLVESETUPID')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'ID',1,ID1)
            IF(ID1.EQ.ID) THEN
!             == CHECK IF THE SAME NAME HAS ALREADY BEEN FOUND =================
              IF(TFOUND) THEN
                CALL ERROR$MSG('AUGMENTATION ID IS AMBIGOUS')
                CALL ERROR$MSG('ID HAS BEEN ENCOUNTERED TWICE')
                CALL ERROR$CHVAL('ID',ID)
                CALL ERROR$STOP('SETUP$RESOLVESETUPID')
              END IF
              TFOUND=.TRUE.
              CALL ATOMLIB$SCNTLLOOKUPONE(LL_SCNTL,AEZ,ZV,COREID,TSO &
       &                              ,RBOX,LX,TYPE &
       &                              ,RCL &
       &                              ,LAMBDA &
       &                              ,RCSM,POTPOW,POTRC,TPOTVAL,POTVAL &
       &                              ,CORPOW,CORRC,TCORVAL,CORVAL &
       &                              ,DMIN,DMAX,RMAX)
            END IF
            CALL LINKEDLIST$SELECT(LL_SCNTL,'..',0)
          ENDDO
        END IF
!
        IF(.NOT.TFOUND) THEN
          CALL ERROR$MSG('SETUP TYPE NOT RECOGNIZED')
          CALL ERROR$CHVAL('STPTYPE',STPTYPE)
          CALL ERROR$CHVAL('ID',ID)
          IF(TSTPFILE) THEN
            CALL ERROR$MSG('SETUP FILE HAS BEEN READ')
            CALL ERROR$MSG('CAUTION: CHANGED NOTATION IN SETUP FILE:')
            CALL ERROR$MSG('USE !ACNTL!AUGMENT INSTEAD OF !SCNTL!SETUP')
          ELSE
            CALL ERROR$MSG('NO SETUP FILE HAS BEEN READ')
          END IF
          CALL ERROR$STOP('SETUP_BUILDPARMS')
        END IF
!    
!       ========================================================================
!       == PLACE INFORMATION INTO LINKEDLIST                                  ==
!       ========================================================================
        CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
        CALL LINKEDLIST$SELECT(LL_STRC,'AUGMENT',0)
        CALL LINKEDLIST$SET(LL_STRC,'ID',0,ID)
        CALL LINKEDLIST$SET(LL_STRC,'EL',0,EL)
        CALL LINKEDLIST$SET(LL_STRC,'ZV',0,ZV)
        CALL LINKEDLIST$SET(LL_STRC,'COREID',0,COREID)
        CALL LINKEDLIST$SET(LL_STRC,'SO',0,TSO)
        CALL LINKEDLIST$SET(LL_STRC,'TYPE',0,TYPE)
        CALL LINKEDLIST$SET(LL_STRC,'RBOX/RCOV',0,RBOX/RCOV)
        CALL LINKEDLIST$SET(LL_STRC,'RCSM/RCOV',0,RCSM/RCOV)
        CALL LINKEDLIST$SET(LL_STRC,'RCL/RCOV',0,RCL/RCOV)
PRINT*,'PAW_SETUPS.F90 A LAMBDA=',LAMBDA
        CALL LINKEDLIST$SET(LL_STRC,'LAMBDA',0,LAMBDA)
!   
        CALL LINKEDLIST$SELECT(LL_STRC,'GRID',0)
          CALL LINKEDLIST$SET(LL_STRC,'DMIN',0,DMIN)
          CALL LINKEDLIST$SET(LL_STRC,'DMAX',0,DMAX)
          CALL LINKEDLIST$SET(LL_STRC,'RMAX',0,RMAX)
        CALL LINKEDLIST$SELECT(LL_STRC,'..',0) !RETURN FROM !GRID
!   
        CALL LINKEDLIST$SELECT(LL_STRC,'POT',0)
          CALL LINKEDLIST$SET(LL_STRC,'POW',0,POTPOW)
          IF(TPOTVAL)CALL LINKEDLIST$SET(LL_STRC,'VAL0',0,POTVAL)
          CALL LINKEDLIST$SET(LL_STRC,'RC/RCOV',0,POTRC/RCOV)
        CALL LINKEDLIST$SELECT(LL_STRC,'..',0) !RETURN FROM !POT
!   
        CALL LINKEDLIST$SELECT(LL_STRC,'CORE',0)
          CALL LINKEDLIST$SET(LL_STRC,'POW',0,CORPOW)
          IF(TCORVAL) CALL LINKEDLIST$SET(LL_STRC,'VAL0',0,CORVAL)
          CALL LINKEDLIST$SET(LL_STRC,'RC/RCOV',0,CORRC/RCOV)
        CALL LINKEDLIST$SELECT(LL_STRC,'..',0)  !RETURN FROM !CORE
!   
        CALL LINKEDLIST$SELECT(LL_STRC,'..',0)  !RETURN FROM !AUGMENT
        CALL LINKEDLIST$SELECT(LL_STRC,'..',0)  !RETURN FROM !SPECIES
      ENDDO ! END OF LOOP OVER SPECIES TYPES

!!$CALL LINKEDLIST$SELECT(LL_STRC,'~',0)  !RETURN FROM !AUGMENT
!!$CALL LINKEDLIST$REPORT(LL_STRC,6)
!CALL SETUP_BUILDPARMS()
!!$CALL SETUP_BUILDPARMS_NDLSS()
!!$!CALL SETUP_IDENTIFYRMAX()
!STOP
                                    CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_BUILDPARMSONE_NDLSS(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
     &                    ,RCL,POTPOW,POTRC,TPOTVAL,POTVAL &
     &                    ,CORPOW,CORRC,TCORVAL,CORVAL &
     &                    ,DMIN,DMAX,RMAX)
!     **************************************************************************
!     ** GENERATES AUTOMATICALLY A PARAMETER FILE FOR THE SETUP CONSTRUCTION  **
!     ** TYPE: NDLSS                                                          **
!     **                                                                      **
!     ** RADII RELATIVE TO RCOV ARE LIMITED BY A MINIMUM AND A MAXIMUM        **
!     ******************************PETER BLOECHL, GOSLAR 2014******************
      USE PERIODICTABLE_MODULE
      USE CONSTANTS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID       ! SETUP ID
      INTEGER(4)  ,INTENT(IN)  :: LX       ! RC PRODUCED UP TO ANGULAR MOM LX.
      REAL(8)     ,INTENT(OUT) :: AEZ      ! ATOMIC NUMBER
      REAL(8)     ,INTENT(OUT) :: ZV       ! #(VALENCE ELECTRONS)
      CHARACTER(*),INTENT(OUT) :: COREID   ! SPECIFIES FROZEN CORE 
      CHARACTER(*),INTENT(OUT) :: TYPE     ! PSEUDIZATION TYPE
      REAL(8)     ,INTENT(OUT) :: RBOX     ! RBOX
      REAL(8)     ,INTENT(OUT) :: RCSM     ! RCSM
      REAL(8)     ,INTENT(OUT) :: RCL(LX+1)! RC(PHI)
      REAL(8)     ,INTENT(OUT) :: POTPOW   ! LEADING POWER FOR PSEUDO POTENTIAL
      REAL(8)     ,INTENT(OUT) :: POTRC    ! RC FOR PSEUDO POTENTIAL
      LOGICAL(4)  ,INTENT(OUT) :: TPOTVAL  ! CONSIDER POTVAL?
      REAL(8)     ,INTENT(OUT) :: POTVAL   ! PSEUDO POTENTIAL VALUE AT R=0
      REAL(8)     ,INTENT(OUT) :: CORPOW   ! LEADING POWER FOR PSEUDO CORE
      REAL(8)     ,INTENT(OUT) :: CORRC    ! RC FOR PSEUDO CORE
      LOGICAL(4)  ,INTENT(OUT) :: TCORVAL  ! CONSIDER CORVAL?
      REAL(8)     ,INTENT(OUT) :: CORVAL   ! PSEUDO CORE VALUE AT R=0
      REAL(8)     ,INTENT(OUT) :: DMIN     ! SMALLEST GRID SPACING
      REAL(8)     ,INTENT(OUT) :: DMAX     ! LARGEST GRID SPACING
      REAL(8)     ,INTENT(OUT) :: RMAX     ! OUTERMOST GRID POINT 
      REAL(8)     ,PARAMETER   :: SCALE=1.D0 
      CHARACTER(32)            :: STPTYPE
      CHARACTER(2)             :: EL       ! ELEMENT SYMBOL
      INTEGER(4)               :: I,L
      INTEGER(4)               :: IZ
      LOGICAL(4)               :: TSC
      REAL(8)                  :: ZCORE
      REAL(8)                  :: RCOV
      REAL(8)                  :: ANGSTROM
      INTEGER(4)               :: NMAIN(0:LX)
!     **************************************************************************
                      CALL TRACE$PUSH('SETUP_BUILDPARMSONE_NDLSS')
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      TYPE='NDLSS'
!
!     ==========================================================================
!     == CHECK ID AND PREPARE SOME DATA                                      ==
!     ==========================================================================
      I=INDEX(ID,'_')
      STPTYPE=ID(I+1:)
      TSC=.FALSE.
      IF(STPTYPE.EQ.'NDLSS_V0') THEN
        TSC=.FALSE.
      ELSE IF(STPTYPE.EQ.'NDLSS_SC_V0') THEN
        TSC=.TRUE.
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$MSG('ID MUST END ON "NDLSS" OR "NDLSS_SC"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_BUILDPARMSONE_NDLSS')
      END IF
      EL=ID(1:2)
      CALL PERIODICTABLE$GET(EL,'Z',IZ)
      CALL PERIODICTABLE$GET(IZ,'R(COV)',RCOV)
!
!     ==========================================================================
!     == DETERMINE NUMBER OF VALENCE ELECTRONS                                ==
!     ==========================================================================
      CALL PERIODICTABLE$GET(IZ,'Z',AEZ)
      CALL PERIODICTABLE$GET(IZ,'ZCORE',ZCORE)
      ZV=AEZ-ZCORE
      CALL PERIODICTABLE$GET(NINT(ZCORE),'SYMBOL',COREID)
!
!     ==========================================================================
!     == DETERMINE MAIN QUANTUM NUMBERS ACCORDIG TO PERIODIC TABLE            ==
!     ==========================================================================
!     == START WITH LOWEST POSSIBLE VALUE
      DO L=0,LX
        NMAIN(L)=L+1
      ENDDO
!     == INCREASE ACCORDING TO NORMAL SHELLS
      IF(IZ.GE. 3) NMAIN(:)=MAX(NMAIN,2)
      IF(IZ.GE.11) NMAIN(:)=MAX(NMAIN,3)
      IF(IZ.GE.19) NMAIN(:)=MAX(NMAIN,4)
      IF(IZ.GE.37) NMAIN(:)=MAX(NMAIN,5)
      IF(IZ.GE.55) NMAIN(:)=MAX(NMAIN,6)
      IF(IZ.GE.87) NMAIN(:)=MAX(NMAIN,7)
!     == TAKE INTO ACCOUNT D AND F ANOMALIES
      IF(IZ.GE.19) NMAIN(2)=3    ! 3D SHELL
      IF(IZ.GE.37) NMAIN(2)=4    ! 4D SHELL
      IF(IZ.GE.55) NMAIN(2)=5    ! 5D SHELL
      IF(IZ.GE.87) NMAIN(2)=6    ! 6D SHELL
      IF(IZ.GT.55) NMAIN(3)=4    ! 4F SHELL
      IF(IZ.GT.87) NMAIN(3)=5    ! 4F SHELL
!
!     ==========================================================================
!     ==  SPECIFY CORE AND VALENCE SHELLS FOR NORMAL AND SEMI-CORE SETUPS     ==
!     ==========================================================================
      IF(TSC) THEN
!       == SETUPS WITH SEMICORE ================================================
        IF(IZ.GE.19.AND.IZ.LE.30) THEN      ! PUT 3SP SHELL IN VALENCE
          NMAIN(0)=NMAIN(0)-1
          NMAIN(1)=NMAIN(1)-1
          ZV=ZV+8
          COREID=TRIM(ADJUSTL(COREID))//'-S-P'
        ELSE IF(IZ.GE.37.AND.IZ.LE.48) THEN ! PUT 4SP SHELL IN VALENCE
          NMAIN(0)=NMAIN(0)-1
          NMAIN(1)=NMAIN(1)-1
          ZV=ZV+8
          COREID=TRIM(ADJUSTL(COREID))//'-S-P'
        ELSE IF(IZ.GE.55.AND.IZ.LE.80) THEN ! PUT 5SP SHELL IN VALENCE
          NMAIN(0)=NMAIN(0)-1
          NMAIN(1)=NMAIN(1)-1
          ZV=ZV+8
          COREID=TRIM(ADJUSTL(COREID))//'-S-P'
          IF(IZ.GT.72) THEN                 ! PUT 4F INTO CORE
            NMAIN(3)=NMAIN(3)+1
            ZV=ZV-14
            COREID=TRIM(ADJUSTL(COREID))//'+F'
          END IF
        ELSE IF(IZ.GE.97.AND.IZ.LE.112) THEN ! PUT 6SP SHELL IN VALENCE
          NMAIN(0)=NMAIN(0)-1
          NMAIN(1)=NMAIN(1)-1
          ZV=ZV+8
          COREID=TRIM(ADJUSTL(COREID))//'-S-P'
          IF(IZ.GT.104) THEN                 ! PUT 5F INTO CORE
            NMAIN(3)=NMAIN(3)+1
            ZV=ZV-14
            COREID=TRIM(ADJUSTL(COREID))//'+F'
          END IF
        END IF
      ELSE
!       == SETUPS WITHOUT SEMICORE =============================================
        IF(IZ.GE.31.AND.IZ.LE.36) THEN     ! PUT 3D-SHELL IN CORE
          NMAIN(2)=4
          ZV=ZV-10
          COREID=TRIM(ADJUSTL(COREID))//'+D'
        ELSE IF(IZ.GE.49.AND.IZ.LE.54) THEN  ! PUT 4D-SHELL IN CORE
          NMAIN(2)=5
          ZV=ZV-10
          COREID=TRIM(ADJUSTL(COREID))//'+D'
        ELSE IF(IZ.GE.72.AND.IZ.LE.80) THEN  ! PUT 4F SHELL IN CORE
          NMAIN(3)=5
          ZV=ZV-14
          COREID=TRIM(ADJUSTL(COREID))//'+F'
        ELSE IF(IZ.GE.81.AND.IZ.LE.86) THEN  ! PUT 4F AND 5D SHELLS IN CORE
          NMAIN(2)=6
          NMAIN(3)=5
          ZV=ZV-24
          COREID=TRIM(ADJUSTL(COREID))//'+D+F'
        ELSE IF(IZ.GE.104.AND.IZ.LE.112) THEN ! PUT 5F SHELL IN CORE
          NMAIN(3)=6
          ZV=ZV-14
          COREID=TRIM(ADJUSTL(COREID))//'+F'
        END IF  
      END IF
!
!     ==========================================================================
!     == DETERMINE RADIAL GRID
!     ==========================================================================
      DO L=0,LX
        IF(AEZ.LT.1.D0) THEN
          CALL SETUP_RADNDLSSMAX(1.D0,NMAIN(L),RCL(L+1))
        ELSE
          CALL SETUP_RADNDLSSMAX(AEZ,NMAIN(L),RCL(L+1))
        END IF
        RCL(L+1)=MIN(1.4D0*RCOV,RCL(L+1))
        RCL(L+1)=MAX(0.5D0*RCOV,RCL(L+1))
      ENDDO
      RCL(:)=RCL(:)*SCALE  ! APPLY SCALE FACTOR
RCL=RCOV
!
!     ==========================================================================
!     == OTHER PARAMETERS
!     ==========================================================================
      RBOX     =1.2D0*RCOV
      RCSM     =0.25D0*RCOV
      POTPOW   =3.D0
      TPOTVAL  =.FALSE.    ! POTVAL NOT USED
      POTVAL   =0.D0
      POTRC    =0.702D0*RCOV
      CORPOW   =3.D0
      TCORVAL  =.FALSE.    ! CORVAL NOT USED
      CORVAL   =0.D0
      CORRC    =0.702D0*RCOV
!
!     ==========================================================================
!     == DETERMINE RADIAL GRID
!     ==========================================================================
      DMIN=1.D-6
      DMAX=1.D-1
!     == RMAX MUST BE LARGER THAN THE BOX RADIUS AND SHOULD AT LEAST BE AS ===
!     == LARGE AS THE BOND DISTANCE TO OXYGEN ================================
      RMAX=20.D0
      IF(RMAX-3.D0*DMAX.LT.RBOX) THEN
        CALL ERROR$MSG('EXTEND OF RADIAL GRID TOO SHORT FOR SPECIFIED RBOX')
        CALL ERROR$CHVAL('SETUP ID',ID)
        CALL ERROR$R8VAL('AEZ',AEZ)
        CALL ERROR$R8VAL('RMAX',RMAX)
        CALL ERROR$R8VAL('RBOX',RBOX)
        CALL ERROR$STOP('SETUP_BUILDPARMS_HBS')
      END IF
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_BUILDPARMSONE_HBS(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
     &               ,RCL,LAMBDA &
     &               ,POTPOW,POTRC,TPOTVAL,POTVAL &
     &               ,CORPOW,CORRC,TCORVAL,CORVAL &
     &               ,DMIN,DMAX,RMAX)
!     **************************************************************************
!     ** GENERATES AUTOMATICALLY A PARAMETER FILE FOR THE SETUP CONSTRUCTION  **
!     ** TYPE: HBS AND HBS_SC                                                 **
!     ******************************PETER BLOECHL, GOSLAR 2014******************
      USE PERIODICTABLE_MODULE
      USE CONSTANTS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID       ! SETUP ID
      INTEGER(4)  ,INTENT(IN)  :: LX       ! RC PRODUCED UP TO ANGULAR MOM LX.
      REAL(8)     ,INTENT(OUT) :: AEZ      ! ATOMIC NUMBER
      REAL(8)     ,INTENT(OUT) :: ZV       ! #(VALENCE ELECTRONS)
      CHARACTER(*),INTENT(OUT) :: COREID   ! SPECIFIES FROZEN CORE 
      CHARACTER(*),INTENT(OUT) :: TYPE     ! PSEUDIZATION TYPE
      REAL(8)     ,INTENT(OUT) :: RBOX     ! RBOX
      REAL(8)     ,INTENT(OUT) :: RCSM     ! RCSM
      REAL(8)     ,INTENT(OUT) :: RCL(LX+1)! RC(PHI)
      REAL(8)     ,INTENT(OUT) :: LAMBDA(LX+1) ! LAMBDA
      REAL(8)     ,INTENT(OUT) :: POTPOW   ! LEADING POWER FOR PSEUDO POTENTIAL
      REAL(8)     ,INTENT(OUT) :: POTRC    ! RC FOR PSEUDO POTENTIAL
      LOGICAL(4)  ,INTENT(OUT) :: TPOTVAL  ! POTVAL RELEVANT?
      REAL(8)     ,INTENT(OUT) :: POTVAL   ! PSEUDO POTENTIAL VALUE AT R=0
      REAL(8)     ,INTENT(OUT) :: CORPOW   ! LEADING POWER FOR PSEUDO CORE
      REAL(8)     ,INTENT(OUT) :: CORRC    ! RC FOR PSEUDO CORE
      LOGICAL(4)  ,INTENT(OUT) :: TCORVAL  ! CORVAL RELEVANT?
      REAL(8)     ,INTENT(OUT) :: CORVAL   ! PSEUDO CORE VALUE AT R=0
      REAL(8)     ,INTENT(OUT) :: DMIN     ! SMALLEST GRID SPACING
      REAL(8)     ,INTENT(OUT) :: DMAX     ! LARGEST GRID SPACING
      REAL(8)     ,INTENT(OUT) :: RMAX     ! OUTERMOST GRID POINT 
      CHARACTER(32)            :: STPTYPE
      CHARACTER(2)             :: EL       ! ELEMENT SYMBOL
      INTEGER(4)               :: I
      INTEGER(4)               :: IZ
      LOGICAL(4)               :: TSC
      REAL(8)                  :: ZCORE
      REAL(8)                  :: RCOV
      REAL(8)                  :: ANGSTROM
!     **************************************************************************
                      CALL TRACE$PUSH('SETUP_BUILDPARMSONE_HBS')
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
!
!     ==========================================================================
!     == CHECK ID AND PREPARE SOME DATA                                      ==
!     ==========================================================================
      I=INDEX(ID,'_')
      STPTYPE=ID(I+1:)
      TSC=.FALSE.
      IF(STPTYPE.EQ.'HBS') THEN
        TSC=.FALSE.
      ELSE IF(STPTYPE.EQ.'HBS_SC') THEN
        TSC=.TRUE.
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$MSG('ID MUST END ON "HBS" OR "HBS_SC"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_BUILDPARMSONE_HBS')
      END IF
      EL=ID(1:2)
      CALL PERIODICTABLE$GET(EL,'Z',IZ)
      CALL PERIODICTABLE$GET(IZ,'R(COV)',RCOV)
      TYPE='HBS'
!
!     ==========================================================================
!     == DETERMINE NUMBER OF VALENCE ELECTRONS                                ==
!     ==========================================================================
      CALL PERIODICTABLE$GET(IZ,'Z',AEZ)
      CALL PERIODICTABLE$GET(IZ,'ZCORE',ZCORE)
      ZV=AEZ-ZCORE
      CALL PERIODICTABLE$GET(NINT(ZCORE),'SYMBOL',COREID)
      IF(IZ.GE.31.AND.IZ.LE.36) ZV=ZV-10.D0
      IF(IZ.GE.49.AND.IZ.LE.54) ZV=ZV-10.D0
      IF(IZ.GE.72.AND.IZ.LE.80) ZV=ZV-14.D0
      IF(IZ.GE.81.AND.IZ.LE.86) ZV=ZV-24.D0
      IF(IZ.GE.104) ZV=ZV-14.D0
      IF(TSC) THEN
        IF(AEZ.GE.11)ZV=AEZ-ZCORE+8.D0
      END IF

      IF(IZ.GE.31.AND.IZ.LE.36) COREID=TRIM(ADJUSTL(COREID))//'+D'
      IF(IZ.GE.49.AND.IZ.LE.54) COREID=TRIM(ADJUSTL(COREID))//'+D'
      IF(IZ.GE.72.AND.IZ.LE.80) COREID=TRIM(ADJUSTL(COREID))//'+F'
      IF(IZ.GE.81.AND.IZ.LE.86) COREID=TRIM(ADJUSTL(COREID))//'+D+F'
      IF(IZ.GE.104)             COREID=TRIM(ADJUSTL(COREID))//'+F'
      IF(TSC) THEN
        IF(AEZ.GE.11) THEN
          ZV=AEZ-ZCORE+8.D0
          CALL PERIODICTABLE$GET(NINT(ZCORE),'SYMBOL',COREID)
          COREID=TRIM(ADJUSTL(COREID))//'-S-P'
        END IF
      END IF
!
!     ==========================================================================
!     == DETERMINE CUTOFF RADII ETC.
!     ==========================================================================
      IF(.NOT.TSC) THEN
        RBOX     =1.2D0*RCOV
        RCSM     =0.25D0*RCOV
        RCL(:)   =0.75D0*RCOV
        LAMBDA(:)=6.D0
        POTPOW   =3.D0
        POTRC    =0.67D0*RCOV
        TPOTVAL  =.FALSE.
        POTVAL   =0.D0
        CORPOW   =3.D0
        CORRC    =0.67D0*RCOV
        TCORVAL  =.FALSE.
        CORVAL   =0.D0
      ELSE
        RBOX     =1.2D0*RCOV
        RCSM     =0.25D0*RCOV
        RCL(:)   =0.5D0*RCOV
        LAMBDA(:)=6.D0
        POTPOW   =3.D0
        POTRC    =0.5D0*RCOV
        TPOTVAL  =.FALSE.
        POTVAL   =-2.2D0  
        CORPOW   =3.D0
        CORRC    =0.5D0*RCOV
        TCORVAL  =.FALSE.
        CORVAL   =0.1D0   
      END IF
!     __ROUND POTRC/RCOV TO THREE DIGITS TO BE CONSISTENT WITH SETUP FILE_______
      POTRC=REAL(INT(1.D+3*POTRC/RCOV),KIND=8)*RCOV*1.D-3
!     __ROUND CORRC/RCOV TO THREE DIGITS TO BE CONSISTENT WITH SETUP FILE_______
      CORRC=REAL(INT(1.D+3*CORRC/RCOV),KIND=8)*RCOV*1.D-3
!
!     ==========================================================================
!     == DETERMINE RADIAL GRID
!     ==========================================================================
      DMIN=1.D-6
      DMAX=1.D-1
!     == RMAX MUST BE LARGER THAN THE BOX RADIUS AND SHOULD AT LEAST BE AS ===
!     == LARGE AS THE BOND DISTANCE TO OXYGEN ================================
      RMAX=9.D0
      IF(RMAX-3.D0*DMAX.LT.RBOX) THEN
        CALL ERROR$MSG('EXTEND OF RADIAL GRID TOO SHORT FOR SPECIFIED RBOX')
        CALL ERROR$CHVAL('SETUP ID',ID)
        CALL ERROR$R8VAL('AEZ',AEZ)
        CALL ERROR$R8VAL('RMAX',RMAX)
        CALL ERROR$R8VAL('RBOX',RBOX)
        CALL ERROR$STOP('SETUP_BUILDPARMS_HBS')
      END IF
!
!     ========================================================================
!     == SEMI-CORE SETUPS                                                   ==
!     ========================================================================
!     == INCLUDE CORE-S-P SHELL INTO THE VALENCE SO THAT ALL STATES 
!     == WITH A MINIMUM MAIN QUANTUM NUMBER ARE INCLUDED
      IF(TSC.AND.AEZ.GE.11) THEN
        ZV=AEZ-ZCORE+8.D0
!       == F-STATES ARE INCLUDED ONLY FOR F-ELEMENTS
!       IF(IZ.GE.72.AND.IZ.LE.86) ZV=ZV-14.D0
!       IF(IZ.GE.104.AND.IZ.LE.118) ZV=ZV-14.D0
      END IF
                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_BUILDPARMSONE_7560(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
     &                                   ,RCL,LAMBDA &
     &                                   ,POTPOW,POTRC,TPOTVAL,POTVAL &
     &                                   ,CORPOW,CORRC,TCORVAL,CORVAL &
     &                                   ,DMIN,DMAX,RMAX)
!     **************************************************************************
!     ** GENERATES AUTOMATICALLY A PARAMETER FILE FOR THE SETUP CONSTRUCTION  **
!     ** TYPE: .75_6.0                                                        **
!     ******************************PETER BLOECHL, GOSLAR 2014******************
      USE PERIODICTABLE_MODULE
      USE CONSTANTS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID       ! SETUP ID
      INTEGER(4)  ,INTENT(IN)  :: LX       ! RC PRODUCED UP TO ANGULAR MOM LX.
      REAL(8)     ,INTENT(OUT) :: AEZ      ! ATOMIC NUMBER
      REAL(8)     ,INTENT(OUT) :: ZV       ! #(VALENCE ELECTRONS)
      CHARACTER(*),INTENT(OUT) :: COREID   ! SPECIFIES FROZEN CORE 
      CHARACTER(*),INTENT(OUT) :: TYPE     ! PSEUDIZATION TYPE
      REAL(8)     ,INTENT(OUT) :: RBOX     ! RBOX
      REAL(8)     ,INTENT(OUT) :: RCSM     ! RCSM
      REAL(8)     ,INTENT(OUT) :: RCL(LX+1)! RC(PHI)
      REAL(8)     ,INTENT(OUT) :: LAMBDA(LX+1) ! LAMBDA
      REAL(8)     ,INTENT(OUT) :: POTPOW   ! LEADING POWER FOR PSEUDO POTENTIAL
      REAL(8)     ,INTENT(OUT) :: POTRC    ! RC FOR PSEUDO POTENTIAL
      LOGICAL(4)  ,INTENT(OUT) :: TPOTVAL  ! POTVAL RELEVANT?
      REAL(8)     ,INTENT(OUT) :: POTVAL   ! PSEUDO POTENTIAL VALUE AT R=0
      REAL(8)     ,INTENT(OUT) :: CORPOW   ! LEADING POWER FOR PSEUDO CORE
      REAL(8)     ,INTENT(OUT) :: CORRC    ! RC FOR PSEUDO CORE
      LOGICAL(4)  ,INTENT(OUT) :: TCORVAL ! CORVAL RELEVANT?
      REAL(8)     ,INTENT(OUT) :: CORVAL   ! PSEUDO CORE VALUE AT R=0
      REAL(8)     ,INTENT(OUT) :: DMIN     ! SMALLEST GRID SPACING
      REAL(8)     ,INTENT(OUT) :: DMAX     ! LARGEST GRID SPACING
      REAL(8)     ,INTENT(OUT) :: RMAX     ! OUTERMOST GRID POINT 
      CHARACTER(32)            :: STPTYPE
      CHARACTER(2)             :: EL       ! ELEMENT SYMBOL
      INTEGER(4)               :: NS,NP,ND,NF
      INTEGER(4)               :: I
      INTEGER(4)               :: IZ
      LOGICAL(4)               :: TSC
      REAL(8)                  :: ZCORE
      REAL(8)                  :: RCOV
!     **************************************************************************
                      CALL TRACE$PUSH('SETUP_BUILDPARMSONE_NDLSS')
!
!     ==========================================================================
!     == CHECK ID AND PREPARE SOME DATA                                      ==
!     ==========================================================================
      I=INDEX(ID,'_')
      STPTYPE=ID(I+1:)
      TSC=.FALSE.
      IF(STPTYPE.EQ.'.75_6.0') THEN
        TSC=.FALSE.
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$MSG('ID MUST END ON ".75_6.0"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_BUILDPARMSONE_NDLSS')
      END IF
      TYPE='HBS'
      EL=ID(1:2)
      CALL PERIODICTABLE$GET(EL,'Z',IZ)
      CALL PERIODICTABLE$GET(IZ,'R(COV)',RCOV)
      CALL PERIODICTABLE$GET(IZ,'Z',AEZ)
      CALL PERIODICTABLE$GET(IZ,'ZCORE',ZCORE)
      ZV=AEZ-ZCORE
      CALL PERIODICTABLE$GET(NINT(ZCORE),'SYMBOL',COREID)
!
!     ==========================================================================
!     == DETERMINE NUMBER OF VALENCE ELECTRONS                                ==
!     ==========================================================================
      CALL PERIODICTABLE$GET(IZ,'OCC(S)',NS)
      CALL PERIODICTABLE$GET(IZ,'OCC(P)',NP)
      CALL PERIODICTABLE$GET(IZ,'OCC(D)',ND)
      CALL PERIODICTABLE$GET(IZ,'OCC(F)',NF)
!     == DO NOT USE D OR F-STATES IN THE VALENCE IF AT LEAST TWO OTHER
!     == ELECTRONS ARE PRESENT.
      IF(ND.EQ.10.AND.NS.EQ.2.AND.NP.GT.0) THEN
        NF=0
        ND=0
        COREID=TRIM(ADJUSTL(COREID))//'+D+F'
      END IF
      ZV=REAL(NS+NP+ND+NF)
!
!     ==========================================================================
!     == DETERMINE CUTOFF RADII ETC.
!     ==========================================================================
      RBOX     =1.2D0*RCOV
      RCSM     =0.25D0*RCOV
      RCL(:)   =0.75D0*RCOV
      LAMBDA(:)=6.D0
      POTPOW   =3.D0
      POTRC    =0.75D0*RCOV-0.1D0
!     __ROUND POTRC/RCOV TO THREE DIGITS TO BE CONSISTENT WITH SETUP FILE_______
      POTRC=REAL(INT(1.D+3*POTRC/RCOV),KIND=8)*RCOV*1.D-3
      TPOTVAL  =.FALSE.
      POTVAL   =999.999D0  !APPARENTLY INCORRECT NUMBER
      CORPOW   =3.D0
      CORRC    =0.75D0*RCOV-0.1D0
!     __ROUND CORRC/RCOV TO THREE DIGITS TO BE CONSISTENT WITH SETUP FILE_______
      CORRC=REAL(INT(1.D+3*CORRC/RCOV),KIND=8)*RCOV*1.D-3
      TCORVAL  =.FALSE.
      CORVAL   =999.999D0  !APPARENTLY INCORRECT NUMBER
!
!     ==========================================================================
!     == DETERMINE RADIAL GRID
!     ==========================================================================
      DMIN=1.D-6
      DMAX=1.D-1
      RMAX=20.
      IF(RMAX-3.D0*DMAX.LT.RBOX) THEN
        CALL ERROR$MSG('EXTEND OF RADIAL GRID TOO SHORT FOR SPECIFIED RBOX')
        CALL ERROR$CHVAL('SETUP ID',ID)
        CALL ERROR$R8VAL('AEZ',AEZ)
        CALL ERROR$R8VAL('RMAX',RMAX)
        CALL ERROR$R8VAL('RBOX',RBOX)
        CALL ERROR$STOP('SETUP_BUILDPARMS_HBS')
      END IF
                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$SCNTLLOOKUPONE(LL_STP_,AEZ,ZV,COREID,TSO &
     &                              ,RBOX,LX,TYPE,RCL,LAMBDA &
     &                              ,RCSM,POW_POT,RC_POT,TVAL0_POT,VAL0_POT &
     &                              ,POW_CORE,RC_CORE,TVAL0_CORE,VAL0_CORE &
     &                              ,DMIN,DMAX,RMAX)
!     **************************************************************************
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN):: LL_STP_
      INTEGER(4)  ,INTENT(IN) :: LX
      REAL(8)     ,INTENT(OUT):: AEZ
      REAL(8)     ,INTENT(OUT):: ZV
      CHARACTER(*),INTENT(OUT):: COREID
      LOGICAL(4)  ,INTENT(OUT):: TSO
      REAL(8)     ,INTENT(OUT):: RBOX
      CHARACTER(*),INTENT(OUT):: TYPE
      REAL(8)     ,INTENT(OUT):: RCL(LX+1)
      REAL(8)     ,INTENT(OUT):: LAMBDA(LX+1)
      REAL(8)     ,INTENT(OUT):: RCSM
      REAL(8)     ,INTENT(OUT):: POW_POT
      REAL(8)     ,INTENT(OUT):: RC_POT
      LOGICAL(4)  ,INTENT(OUT):: TVAL0_POT
      REAL(8)     ,INTENT(OUT):: VAL0_POT
      REAL(8)     ,INTENT(OUT):: POW_CORE
      REAL(8)     ,INTENT(OUT):: RC_CORE
      LOGICAL(4)  ,INTENT(OUT):: TVAL0_CORE
      REAL(8)     ,INTENT(OUT):: VAL0_CORE
      REAL(8)     ,INTENT(OUT):: DMIN
      REAL(8)     ,INTENT(OUT):: DMAX
      REAL(8)     ,INTENT(OUT):: RMAX
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      TYPE(LL_TYPE)           :: LL_STP
      LOGICAL(4)              :: TCHK
      CHARACTER(128)          :: ID
      REAL(8)                 :: RCOV
!     **************************************************************************
      LL_STP=LL_STP_  !AVOID REPOSITIONING OF THE POINTER

      CALL LINKEDLIST$EXISTD(LL_STP,'ID',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:ID NOT SPECIFIED')
        CALL ERROR$STOP('ATOMLIB$SCNTLLOOKUPONE')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'ID',1,ID)
!
!     ==========================================================================
!     == IDENTIFY ATOM                                                        ==
!     ==========================================================================
      CALL SETUP_LOOKUP_GENERIC(LL_STP_,ID,AEZ,ZV,RCSM,COREID,TSO)
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
!
!     ==========================================================================
!     == PARAMETER FOR PARTIALWAVE CONSTRUCTION                               ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'TYPE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('SETUP TYPE NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('ATOMLIB$SCNTLLOOKUPONE')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'TYPE',1,TYPE)
!
      IF(TYPE.EQ.'KERKER') THEN
        CALL SETUP_LOOKUP_KERKER(LL_STP_,ID,RCOV,LX,RCL,RBOX)
      ELSE IF(TYPE.EQ.'HBS') THEN
        CALL SETUP_LOOKUP_HBS(LL_STP_,ID,RCOV,LX,RCL,LAMBDA,RBOX)
      ELSE IF(TYPE.EQ.'NDLSS') THEN
        CALL SETUP_LOOKUP_NDLSS(LL_STP_,ID,RCOV,LX,RCL,RBOX)
      ELSE
        CALL ERROR$MSG('SETUP TYPE NOT RECOGNIZED')
        CALL ERROR$CHVAL('TYPE',TYPE)
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('ATOMLIB$SCNTLLOOKUPONE')
      END IF
!
!     ==========================================================================
!     == COLLECT RADIAL GRID                                                  ==
!     ==========================================================================
      CALL SETUP_LOOKUP_GRID(LL_STP_,ID,DMIN,DMAX,RMAX)
!
!     ==========================================================================
!     == DEFINE POTENTIAL PSEUDIZATION                                        ==
!     ==========================================================================
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
      CALL SETUP_LOOKUP_PSPOT(LL_STP_,ID,RCOV,POW_POT,RC_POT,TVAL0_POT,VAL0_POT)
!
!     ==========================================================================
!     == DEFINE CORE PSEUDIZATION                                             ==
!     ==========================================================================
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
      CALL SETUP_LOOKUP_PSCORE(LL_STP_,ID,RCOV &
     &                       ,POW_CORE,RC_CORE,TVAL0_CORE,VAL0_CORE)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_LOOKUPSETUP(LL_STP_,ID,AEZ,ZV,COREID,TSO &
     &                              ,RBOX,LX,TYPE,RCL,LAMBDA &
     &                              ,RCSM,POW_POT,RC_POT,TVAL0_POT,VAL0_POT &
     &                              ,POW_CORE,RC_CORE,TVAL0_CORE,VAL0_CORE &
     &                              ,DMIN,DMAX,RMAX)
!     **************************************************************************
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN) :: LL_STP_
      INTEGER(4)   ,INTENT(IN) :: LX
      CHARACTER(*) ,INTENT(OUT):: ID
      REAL(8)      ,INTENT(OUT):: AEZ
      REAL(8)      ,INTENT(OUT):: ZV
      CHARACTER(*) ,INTENT(OUT):: COREID
      LOGICAL(4)   ,INTENT(OUT):: TSO      ! SWITCH FOR SPIN-ORBIT COUPLING
      REAL(8)      ,INTENT(OUT):: RBOX
      CHARACTER(32),INTENT(OUT):: TYPE
      REAL(8)      ,INTENT(OUT):: RCL(LX+1)
      REAL(8)      ,INTENT(OUT):: LAMBDA(LX+1)
      REAL(8)      ,INTENT(OUT):: RCSM
      REAL(8)      ,INTENT(OUT):: POW_POT
      REAL(8)      ,INTENT(OUT):: RC_POT
      LOGICAL(4)   ,INTENT(OUT):: TVAL0_POT
      REAL(8)      ,INTENT(OUT):: VAL0_POT
      REAL(8)      ,INTENT(OUT):: POW_CORE
      REAL(8)      ,INTENT(OUT):: RC_CORE
      LOGICAL(4)   ,INTENT(OUT):: TVAL0_CORE
      REAL(8)      ,INTENT(OUT):: VAL0_CORE
      REAL(8)      ,INTENT(OUT):: DMIN
      REAL(8)      ,INTENT(OUT):: DMAX
      REAL(8)      ,INTENT(OUT):: RMAX
      TYPE(LL_TYPE)            :: LL_STP
      LOGICAL(4)               :: TCHK
      CHARACTER(64)            :: STRING
      REAL(8)                  :: RCOV
      REAL(8)      ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)      ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
!     **************************************************************************
      LL_STP=LL_STP_  !AVOID REPOSITIONING OF THE POINTER

      CALL LINKEDLIST$EXISTL(LL_STP,'AUGMENT',0,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('BLOCK !AUGMENT MISSING')
        CALL ERROR$STOP('SETUP_LOOKUPSETUP')
      END IF
      CALL LINKEDLIST$SELECT(LL_STP,'AUGMENT')
      CALL LINKEDLIST$EXISTD(LL_STP,'ID',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:ID NOT SPECIFIED')
        CALL ERROR$STOP('SETUP_LOOKUPSETUP')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'ID',1,ID)
!
!     ==========================================================================
!     == IDENTIFY ATOM                                                        ==
!     ==========================================================================
      CALL SETUP_LOOKUP_GENERIC(LL_STP,ID,AEZ,ZV,RCSM,COREID,TSO)
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
!
!     ==========================================================================
!     == PARAMETER FOR PARTIALWAVE CONSTRUCTION                               ==
!     ==========================================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'TYPE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('AUGMENTATION METHOD "TYPE" NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('ATOMLIB$SCNTLLOOKUPONE')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'TYPE',1,STRING)
      TYPE=TRIM(STRING) ! COPY REQUIRED BECAUSE OF MAPPING FROM CH(64) TO CH(32)
!
      IF(TYPE.EQ.'KERKER') THEN
        LAMBDA=0.D0  !NOT USED
        CALL SETUP_LOOKUP_KERKER(LL_STP,ID,RCOV,LX,RCL,RBOX)
      ELSE IF(TYPE.EQ.'HBS') THEN
        CALL SETUP_LOOKUP_HBS(LL_STP,ID,RCOV,LX,RCL,LAMBDA,RBOX)
      ELSE IF(TYPE.EQ.'NDLSS') THEN
        LAMBDA=0.D0  !NOT USED
        CALL SETUP_LOOKUP_NDLSS(LL_STP,ID,RCOV,LX,RCL,RBOX)
      ELSE
        CALL ERROR$MSG('AUGMENTATION METHOD "TYPE" NOT RECOGNIZED')
        CALL ERROR$CHVAL('TYPE',TYPE)
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('ATOMLIB$SCNTLLOOKUPONE')
      END IF
!
!     ==========================================================================
!     == COLLECT RADIAL GRID                                                  ==
!     ==========================================================================
      CALL SETUP_LOOKUP_GRID(LL_STP,ID,DMIN,DMAX,RMAX)
!
!     ==========================================================================
!     == DEFINE POTENTIAL PSEUDIZATION                                        ==
!     ==========================================================================
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
      CALL SETUP_LOOKUP_PSPOT(LL_STP,ID,RCOV,POW_POT,RC_POT,TVAL0_POT,VAL0_POT)
!
!     ==========================================================================
!     == DEFINE CORE PSEUDIZATION                                             ==
!     ==========================================================================
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
      CALL SETUP_LOOKUP_PSCORE(LL_STP,ID,RCOV &
     &                       ,POW_CORE,RC_CORE,TVAL0_CORE,VAL0_CORE)
      CALL LINKEDLIST$SELECT(LL_STP,'..')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_LOOKUP_GENERIC(LL_STP_,ID,AEZ,ZV,RCSM,COREID,TSO)
!     **************************************************************************
!     **  COLLECT INFORMATION ON CORE PSEUDIZATION FROM CORE BLOCK            **
!     **  LINKED LIST (LL_STP_) MUST BE POSITIONED IN THE PARENT OF THE       **
!     **  GRID BLOCK                                                          **
!     **                                                                      **
!     **  CAUTION: VAL0 IS THE VALUE OF THE PSEUDO CORE NOT ITS RADIAL PART   **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN)  :: LL_STP_  !IDENTIFIER OF THE LINKEDLIST
      CHARACTER(*) ,INTENT(IN)  :: ID       ! ATOM ID
      REAL(8)      ,INTENT(OUT) :: AEZ      ! ATOMIC NUMBER
      REAL(8)      ,INTENT(OUT) :: ZV
      REAL(8)      ,INTENT(OUT) :: RCSM     ! SMALL GAUSSIAN DECAY FOR COMP.CH. 
      CHARACTER(*) ,INTENT(OUT) :: COREID   ! IDENTIFIER FOR THE FROZEN CORE
      LOGICAL(4)   ,INTENT(OUT) :: TSO      ! SWITCH FOR SPIN-ORBIT COUPLING
      TYPE(LL_TYPE)             :: LL_STP
      LOGICAL(4)                :: TCHK,TCHK1,TCHK2
      CHARACTER(2)              :: EL
      REAL(8)                   :: RCOV
      INTEGER(4)                :: ZCORE
!     **************************************************************************
      LL_STP=LL_STP_ !AVOID REPOSITIONING OF THE POINTER
!
!     ==========================================================================
!     == IDENTIFY ATOM                                                        ==
!     ==========================================================================
!     == COLLECT ATOMIC NUMBER =================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'EL',1,TCHK1)
      CALL LINKEDLIST$EXISTD(LL_STP,'Z',1,TCHK2)
      IF(TCHK1.AND.TCHK2) THEN
        CALL ERROR$MSG('!AUGMENT:AEZ AND EL="XX"')
        CALL ERROR$MSG('MUST NOT BE SPECIFIED SIMULTANEOUSLY')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUPGENERIC')
      END IF 
      IF(.NOT.(TCHK1.OR.TCHK2)) THEN
        CALL ERROR$MSG('NEITHER EL NOR AEZ SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUPGENERIC')
      END IF 
      IF(TCHK1) THEN
        CALL LINKEDLIST$GET(LL_STP,'EL',1,EL)
        CALL PERIODICTABLE$GET(EL,'Z',AEZ)
      ELSE 
        CALL LINKEDLIST$GET(LL_STP,'Z',1,AEZ)
      END IF
!
!     == COLLECT #VALENCE ELECTRONS ============================================
      CALL LINKEDLIST$EXISTD(LL_STP,'CORE',0,TCHK)
      IF(.NOT.TCHK) THEN
        IF(AEZ.GE. 71.AND.AEZ.LE. 86) COREID='+F'   ! FREEZE 4F SHELL FOR LU-RN
        IF(AEZ.GE. 89.AND.AEZ.LE.118) COREID='+F'   ! FREEZE 5F SHELL FOR LR-OG
        IF(AEZ.GE. 31.AND.AEZ.LE. 36) COREID='+D'   ! FREEZE 3D SHELL FOR GA-KR
        IF(AEZ.GE. 49.AND.AEZ.LE. 54) COREID='+D'   ! FREEZE 4D SHELL FOR IN-XE
        IF(AEZ.GE. 81.AND.AEZ.LE. 86) COREID='+D'   ! FREEZE 5D SHELL FOR TL-RN
        IF(AEZ.GE.113.AND.AEZ.LE.118) COREID='+D'   ! FREEZE 6D SHELL FOR NH-OG
!
!       ========================================================================
!       == OVERWRITE DEFAULT IF ZV IS SPECIFIED. ===============================
!       == THIS OPTION MAINTAINS BACKWARD COMPATIBILITY ========================
!       == AND IS MARKED FOR DELETION ==========================================
!       ========================================================================
        CALL LINKEDLIST$EXISTD(LL_STP,'ZV',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('NEITHER COREID NOR ZV NOT SPECIFIED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP_LOOKUPGENERIC')
        END IF
        CALL LINKEDLIST$GET(LL_STP,'ZV',1,ZV)
        COREID='0'
        IF(AEZ.GE.2) COREID='HE'
        IF(AEZ.GE.10)COREID='NE'
        IF(AEZ.GE.18)COREID='AR'
        IF(AEZ.GE.36)COREID='KR'
        IF(AEZ.GE.54)COREID='XE'
        IF(AEZ.GE.86)COREID='RN'
        ZCORE=0
        IF(AEZ.GE.2) ZCORE=2
        IF(AEZ.GE.10)ZCORE=10
        IF(AEZ.GE.18)ZCORE=18
        IF(AEZ.GE.36)ZCORE=36
        IF(AEZ.GE.54)ZCORE=54
        IF(AEZ.GE.86)ZCORE=86
        ZCORE=INT(AEZ-ZV)-ZCORE
        IF(ZCORE.EQ.0) THEN
         COREID=TRIM(ADJUSTL(COREID))
        ELSE IF(ZCORE.EQ.-2) THEN
         COREID=TRIM(ADJUSTL(COREID))//'-S'
        ELSE IF(ZCORE.EQ.-6) THEN
          COREID=TRIM(ADJUSTL(COREID))//'-P'
        ELSE IF(ZCORE.EQ.-8)  THEN
          COREID=TRIM(ADJUSTL(COREID))//'-S-P'
        ELSE IF(ZCORE.EQ.-10) THEN
          COREID=TRIM(ADJUSTL(COREID))//'-D'
        ELSE IF(ZCORE.EQ.-12) THEN
          COREID=TRIM(ADJUSTL(COREID))//'-S-D'
        ELSE IF(ZCORE.EQ.-16) THEN
          COREID=TRIM(ADJUSTL(COREID))//'-P-D'
        ELSE IF(ZCORE.EQ.-18) THEN
          COREID=TRIM(ADJUSTL(COREID))//'-S-P-D'
        ELSE IF(ZCORE.EQ.+2) THEN
          COREID=TRIM(ADJUSTL(COREID))//'+S'
        ELSE IF(ZCORE.EQ.+10) THEN
          COREID=TRIM(ADJUSTL(COREID))//'+D'
        ELSE IF(ZCORE.EQ.+14) THEN
          COREID=TRIM(ADJUSTL(COREID))//'+F'
        ELSE IF(ZCORE.EQ.+24) THEN
          COREID=TRIM(ADJUSTL(COREID))//'+D+F'
        ELSE
          CALL ERROR$MSG('VALUE OF ZV NOT RECOGNIZED')
          CALL ERROR$MSG('"ZV" IS OBSOLETE')
          CALL ERROR$MSG('SPECIFY "COREID" INSTEAD')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$R8VAL('Z',AEZ)
          CALL ERROR$R8VAL('ZV',ZV)
          CALL ERROR$I4VAL('ZCORE',ZCORE)
          CALL ERROR$STOP('SETUP_LOOKUPGENERIC')
        END IF
!
!       == DISCARD ZV TO AVOID ITS USAGE. IT WILL BE LATER RECALCULATED FROM ===
!       == COREID ==============================================================
        ZV=999999.D0
      ELSE      
        CALL LINKEDLIST$GET(LL_STP,'CORE',1,COREID)
        COREID=+COREID
      END IF
!
!     == DECAY CONSTANT FOR SHORT-RANGED COMPENSATION CHARGE DENSITY ===========
      CALL LINKEDLIST$EXISTD(LL_STP,'RCSM/RCOV',1,TCHK)
      RCSM=0.25D0
      IF(TCHK)CALL LINKEDLIST$GET(LL_STP,'RCSM/RCOV',1,RCSM)
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
      RCSM=RCSM*RCOV
!
!     == SWITCH FOR SPIN-ORBIT COUPLING ========================================
      TSO=.FALSE.
      CALL LINKEDLIST$EXISTD(LL_STP,'SO',1,TCHK)
      IF(TCHK)CALL LINKEDLIST$GET(LL_STP,'SO',1,TSO)

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_LOOKUP_HBS(LL_STP_,ID,RCOV,LX,RC,LAMBDA,RBOX)
!     **************************************************************************
!     **  COLLECT INFORMATION ON CORE PSEUDIZATION FROM CORE BLOCK            **
!     **  LINKED LIST (LL_STP_) MUST BE POSITIONED IN THE PARENT OF THE       **
!     **  GRID BLOCK                                                          **
!     **                                                                      **
!     **  CAUTION: VAL0 IS THE VALUE OF THE PSEUDO CORE NOT ITS RADIAL PART   **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN)  :: LL_STP_
      CHARACTER(*) ,INTENT(IN)  :: ID
      REAL(8)      ,INTENT(IN)  :: RCOV
      INTEGER(4)   ,INTENT(IN)  :: LX
      REAL(8)      ,INTENT(OUT) :: RC(LX+1)
      REAL(8)      ,INTENT(OUT) :: LAMBDA(LX+1)
      REAL(8)      ,INTENT(OUT) :: RBOX
      TYPE(LL_TYPE)             :: LL_STP
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: LENG
      REAL(8)      ,ALLOCATABLE :: WORK(:)
!     **************************************************************************
      LL_STP=LL_STP_  !AVOID REPOSITIONING OF THE POINTER
!
!     == COLLECT BOX RADIUS ====================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'RBOX/RCOV',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:RBOX/RCOV NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_HBS')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'RBOX/RCOV',1,RBOX)
      RBOX=RBOX*RCOV 
!
!     == ARRAY OF CUTOFF RADII RC ==============================================
      CALL LINKEDLIST$EXISTD(LL_STP,'RCL/RCOV',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:RCL NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_HBS')
      END IF
!
      CALL LINKEDLIST$SIZE(LL_STP,'RCL/RCOV',1,LENG)
      IF(LENG.LT.LX+1) THEN
        CALL ERROR$MSG('!AUGMENT:RCL ARRAY TOO SHORT')
        CALL ERROR$I4VAL('MAX L ON FILE ',LENG-1)
        CALL ERROR$I4VAL('MAX L REQUESTED',LX) 
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_HBS')
      END IF
!
      ALLOCATE(WORK(LENG))
      CALL LINKEDLIST$GET(LL_STP,'RCL/RCOV',1,WORK)
      RC=WORK(:LX+1)*RCOV
      DEALLOCATE(WORK)
!
!     == ADDITIONAL LAMBDA PARAMETER REQUIRED FOR 'HBS' TYPE CONSTRUCTION ======
      CALL LINKEDLIST$EXISTD(LL_STP,'LAMBDA',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:LAMBDA NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_HBS')
      END IF
!     ___CHECK LENGTH OF ARRAY__________________________________________________
      CALL LINKEDLIST$SIZE(LL_STP,'LAMBDA',1,LENG)
      IF(LENG.LT.LX+1) THEN
        CALL ERROR$MSG('!AUGMENT:RCL ARRAY TOO SHORT')
        CALL ERROR$I4VAL('MAX L ON FILE ',LENG-1)
        CALL ERROR$I4VAL('MAX L REQUESTED',LX) 
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_HBS')
      END IF
!
      ALLOCATE(WORK(LENG))
      CALL LINKEDLIST$GET(LL_STP,'LAMBDA',1,WORK)
      LAMBDA=WORK(:LX+1)
      DEALLOCATE(WORK)
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_LOOKUP_KERKER(LL_STP_,ID,RCOV,LX,RC,RBOX)
!     **************************************************************************
!     **  COLLECT INFORMATION ON CORE PSEUDIZATION FROM CORE BLOCK            **
!     **  LINKED LIST (LL_STP_) MUST BE POSITIONED IN THE PARENT OF THE       **
!     **  GRID BLOCK                                                          **
!     **                                                                      **
!     **  CAUTION: VAL0 IS THE VALUE OF THE PSEUDO CORE NOT ITS RADIAL PART   **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN)  :: LL_STP_
      CHARACTER(*) ,INTENT(IN)  :: ID
      REAL(8)      ,INTENT(IN)  :: RCOV
      INTEGER(4)   ,INTENT(IN)  :: LX
      REAL(8)      ,INTENT(OUT) :: RC(LX+1)
      REAL(8)      ,INTENT(OUT) :: RBOX
      TYPE(LL_TYPE)             :: LL_STP
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: LENG
      REAL(8)      ,ALLOCATABLE :: WORK(:)
!     **************************************************************************
      LL_STP=LL_STP_  !AVOID REPOSITIONING OF THE POINTER
!
!     == COLLECT BOX RADIUS ====================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'RBOX/RCOV',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:RBOX/RCOV NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_KERKER')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'RBOX/RCOV',1,RBOX)
      RBOX=RBOX*RCOV 
!
!     == ARRAY OF CUTOFF RADII RC ==============================================
      CALL LINKEDLIST$EXISTD(LL_STP,'RCL/RCOV',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:RCL NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_KERKER')
      END IF
!
      CALL LINKEDLIST$SIZE(LL_STP,'RCL/RCOV',1,LENG)
      IF(LENG.LT.LX+1) THEN
        CALL ERROR$MSG('!AUGMENT:RCL ARRAY TOO SHORT')
        CALL ERROR$I4VAL('MAX L ON FILE ',LENG-1)
        CALL ERROR$I4VAL('MAX L REQUESTED',LX) 
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_KERKER')
      END IF
!
      ALLOCATE(WORK(LENG))
      CALL LINKEDLIST$GET(LL_STP,'RCL/RCOV',1,WORK)
      RC=WORK(:LX+1)*RCOV
      DEALLOCATE(WORK)

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_LOOKUP_NDLSS(LL_STP_,ID,RCOV,LX,RC,RBOX)
!     **************************************************************************
!     **  COLLECT INFORMATION ON CORE PSEUDIZATION FROM CORE BLOCK            **
!     **  LINKED LIST (LL_STP_) MUST BE POSITIONED IN THE PARENT OF THE       **
!     **  GRID BLOCK                                                          **
!     **                                                                      **
!     **  CAUTION: VAL0 IS THE VALUE OF THE PSEUDO CORE NOT ITS RADIAL PART   **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN)  :: LL_STP_
      CHARACTER(*) ,INTENT(IN)  :: ID
      REAL(8)      ,INTENT(IN)  :: RCOV
      INTEGER(4)   ,INTENT(IN)  :: LX
      REAL(8)      ,INTENT(OUT) :: RC(LX+1)
      REAL(8)      ,INTENT(OUT) :: RBOX
      TYPE(LL_TYPE)             :: LL_STP
      LOGICAL(4)                :: TCHK
      INTEGER(4)                :: LENG
      REAL(8)      ,ALLOCATABLE :: WORK(:)
!     **************************************************************************
      LL_STP=LL_STP_  !AVOID REPOSITIONING OF THE POINTER
!
!     == COLLECT BOX RADIUS ====================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'RBOX/RCOV',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:RBOX/RCOV NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_NDLSS')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'RBOX/RCOV',1,RBOX)
      RBOX=RBOX*RCOV 
!
!     == ARRAY OF CUTOFF RADII RC ==============================================
      CALL LINKEDLIST$EXISTD(LL_STP,'RCL/RCOV',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:RCL NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_NDLSS')
      END IF
!
      CALL LINKEDLIST$SIZE(LL_STP,'RCL/RCOV',1,LENG)
      IF(LENG.LT.LX+1) THEN
        CALL ERROR$MSG('!AUGMENT:RCL ARRAY TOO SHORT')
        CALL ERROR$I4VAL('MAX L ON FILE ',LENG-1)
        CALL ERROR$I4VAL('MAX L REQUESTED',LX) 
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_NDLSS')
      END IF
!
      ALLOCATE(WORK(LENG))
      CALL LINKEDLIST$GET(LL_STP,'RCL/RCOV',1,WORK)
      RC=WORK(:LX+1)*RCOV
      DEALLOCATE(WORK)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_LOOKUP_PSPOT(LL_STP_,ID,RCOV,POW,RC,TVAL0,VAL0)
!     **************************************************************************
!     **  COLLECT INFORMATION ON CORE PSEUDIZATION FROM CORE BLOCK            **
!     **  LINKED LIST (LL_STP_) MUST BE POSITIONED IN THE PARENT OF THE       **
!     **  GRID BLOCK                                                          **
!     **                                                                      **
!     **  CAUTION: VAL0 IS THE VALUE OF THE PSEUDO CORE NOT ITS RADIAL PART   **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN)  :: LL_STP_
      CHARACTER(*) ,INTENT(IN)  :: ID
      REAL(8)      ,INTENT(IN)  :: RCOV
      REAL(8)      ,INTENT(OUT) :: POW
      REAL(8)      ,INTENT(OUT) :: RC
      LOGICAL(4)   ,INTENT(OUT) :: TVAL0
      REAL(8)      ,INTENT(OUT) :: VAL0
      REAL(8)      ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)      ,PARAMETER   :: Y0=1.D0/SQRT(4.D0*PI)
      TYPE(LL_TYPE)             :: LL_STP
      LOGICAL(4)                :: TCHK
!     **************************************************************************
      LL_STP=LL_STP_  !AVOID REPOSITIONING OF THE POINTER
!
      CALL LINKEDLIST$EXISTL(LL_STP,'POT',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!POT NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_PSPOT')
      END IF
      CALL LINKEDLIST$SELECT(LL_STP,'POT')
!  
!     == POW ===================================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'POW',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!POT:POW NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_PSPOT')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'POW',1,POW)

!     == RC ====================================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'RC/RCOV',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!POT:RC/RCOV NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_PSPOT')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'RC/RCOV',1,RC)
      RC=RC*RCOV
!
      VAL0=0.D0
      CALL LINKEDLIST$EXISTD(LL_STP,'VAL0',1,TCHK)
      TVAL0=TCHK
      IF(TCHK)CALL LINKEDLIST$GET(LL_STP,'VAL0',1,VAL0)
      VAL0=VAL0/Y0
      CALL LINKEDLIST$SELECT(LL_STP,'..')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_LOOKUP_PSCORE(LL_STP_,ID,RCOV,POW,RC,TVAL0,VAL0)
!     **************************************************************************
!     **  COLLECT INFORMATION ON CORE PSEUDIZATION FROM CORE BLOCK            **
!     **  LINKED LIST (LL_STP_) MUST BE POSITIONED IN THE PARENT OF THE       **
!     **  GRID BLOCK                                                          **
!     **                                                                      **
!     **  CAUTION: VAL0 IS THE VALUE OF THE PSEUDO CORE NOT ITS RADIAL PART   **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN)  :: LL_STP_
      CHARACTER(*) ,INTENT(IN)  :: ID
      REAL(8)      ,INTENT(IN)  :: RCOV
      REAL(8)      ,INTENT(OUT) :: POW
      REAL(8)      ,INTENT(OUT) :: RC
      LOGICAL(4)   ,INTENT(OUT) :: TVAL0
      REAL(8)      ,INTENT(OUT) :: VAL0
      REAL(8)      ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)      ,PARAMETER   :: Y0=1.D0/SQRT(4.D0*PI)
      TYPE(LL_TYPE)             :: LL_STP
      LOGICAL(4)                :: TCHK
!     **************************************************************************
      LL_STP=LL_STP_  !AVOID REPOSITIONING OF THE POINTER
!
      CALL LINKEDLIST$EXISTL(LL_STP,'CORE',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!CORE NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_PSCORE')
      END IF
      CALL LINKEDLIST$SELECT(LL_STP,'CORE')
!  
!     == POW ===================================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'POW',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!CORE:POW NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_PSCORE')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'POW',1,POW)

!     == RC ====================================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'RC/RCOV',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!CORE:RC/RCOV NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_PSCORE')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'RC/RCOV',1,RC)
      RC=RC*RCOV
!
      VAL0=0.D0
      CALL LINKEDLIST$EXISTD(LL_STP,'VAL0',1,TCHK)
      TVAL0=TCHK
      IF(TCHK)CALL LINKEDLIST$GET(LL_STP,'VAL0',1,VAL0)
      VAL0=VAL0/Y0
      CALL LINKEDLIST$SELECT(LL_STP,'..')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_LOOKUP_GRID(LL_STP_,ID,DMIN,DMAX,RMAX)
!     **************************************************************************
!     **  COLLECT GRID INFORMATION FROM GRID BLOCK                            **
!     **  LINKED LIST (LL_STP_) MUST BE POSITIONED IN THE PARENT OF THE       **
!     **  GRID BLOCK                                                          **
!     **************************************************************************
      USE LINKEDLIST_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(IN)  :: LL_STP_
      CHARACTER(*) ,INTENT(IN)  :: ID      !AUGMENT ID (USED ONLY FOR ERROR MSG)
      REAL(8)      ,INTENT(OUT) :: DMIN
      REAL(8)      ,INTENT(OUT) :: DMAX
      REAL(8)      ,INTENT(OUT) :: RMAX
      TYPE(LL_TYPE)             :: LL_STP
      LOGICAL(4)                :: TCHK
!     **************************************************************************
      LL_STP=LL_STP_  !AVOID REPOSITIONING OF THE POINTER
!
      CALL LINKEDLIST$EXISTL(LL_STP,'GRID',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!GRID NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_GRID')
      END IF
      CALL LINKEDLIST$SELECT(LL_STP,'GRID')   ! ENTER GRID BLOCK

!     == DMIN ==================================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'DMIN',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:DMIN NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_GRID')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'DMIN',1,DMIN)
!
!     == DMAX ==================================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'DMAX',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:DMAX NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_GRID')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'DMAX',1,DMAX)
!
!     == RMAX ==================================================================
      CALL LINKEDLIST$EXISTD(LL_STP,'RMAX',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('!AUGMENT:RMAX NOT SPECIFIED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP_LOOKUP_GRID')
      END IF
      CALL LINKEDLIST$GET(LL_STP,'RMAX',1,RMAX)
      CALL LINKEDLIST$SELECT(LL_STP,'..')  ! LEAVE !GRID BLOCK
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_READ
!     **************************************************************************
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      USE SETUP_MODULE
      USE LINKEDLIST_MODULE
      USE RADIALFOCK_MODULE, ONLY: VFOCK_TYPE
      IMPLICIT NONE
      TYPE(LL_TYPE)         :: LL_STP
      INTEGER(4),PARAMETER  :: NBX=40
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)   ,PARAMETER  :: FOURPI=4.D0*PI
      REAL(8)   ,PARAMETER  :: C0LL=1.D0/SQRT(FOURPI)
      INTEGER(4)            :: LOFI(NBX)  ! AZIMUTAL ANGULAR MOMENTUM
      INTEGER(4)            :: SOFI(NBX)  ! SPIN ORBIT ORIENTATION 
      REAL(8)               :: FOFI(NBX)  ! OCCUPATION OF THE SHELL
      INTEGER(4)            :: NNOFI(NBX) ! NUMBER OF NODES OF THE PARTIAL WAVE
      REAL(8)               :: EOFI(NBX)  ! ENERGY LEVEL OF THE ATOMIC CALC.
      INTEGER(4)            :: GID        ! GRID ID
      INTEGER(4)            :: GIDG       ! GRID ID FOR G-GRID
      INTEGER(4)            :: NG         ! #(RECIPROCAL GRID POINTS)
      INTEGER(4)            :: NR         ! #(REAL-SPACE GRID POINTS)
      INTEGER(4)            :: LX,LNX
      INTEGER(4)            :: NB
      INTEGER(4)            :: NC         ! #(CORE STATES)
      INTEGER(4),ALLOCATABLE:: NPRO(:)
      REAL(8)   ,ALLOCATABLE:: R(:)       ! RADIAL GRID
      INTEGER(4),ALLOCATABLE:: LOX(:)     ! AZ. ANGULAR MOMENTUM OF PARTIAL WAVE
      REAL(8)   ,ALLOCATABLE:: RC(:)
      LOGICAL(4),ALLOCATABLE:: TC(:)      ! FROZEN-CORE STATE SELECTOR
      REAL(8)   ,ALLOCATABLE:: LAMBDA(:)  ! SECOND PARAMETER FOR PSEUDIZATION
      REAL(8)               :: RBOX
      REAL(8)               :: ROUT
      REAL(8)               :: AEZ        !ATOMIC NUMBER
      REAL(8)               :: ZV         ! #(VALENCE ELECTRONS)
      REAL(8)               :: ETOT
      CHARACTER(64)         :: KEY
      REAL(8)   ,ALLOCATABLE:: PSI(:,:)
      REAL(8)   ,ALLOCATABLE:: PSISM(:,:)  !SMALL COMPONENT
      LOGICAL(4)            :: TCHK,TCHK1
      INTEGER(4)            :: IR,IB,LN,L,IB1,LENG
      INTEGER(4)            :: LRHOX
      CHARACTER(64)         :: STRING
      TYPE(VFOCK_TYPE)      :: VFOCK
      REAL(8)   ,ALLOCATABLE:: AUX(:)
      REAL(8)               :: SVAR
      INTEGER(4),ALLOCATABLE:: IWORK(:)
!     **************************************************************************
                            CALL TRACE$PUSH('SETUP_READ')
      CALL TIMING$CLOCKON('SETUP CONSTRUCTION')
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('SETUP_READ')
      END IF
      WRITE(*,'(80("="))')
      WRITE(*,'(80("="),T20,"  SETUP CONSTRUCTION FOR ",A,"  ")')TRIM(THIS%ID)
      WRITE(*,'(80("="))')
      CALL LINKEDLIST$NEW(LL_STP)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',0)
!
!     ==========================================================================
!     ==  THE FOLLOWING DATA HAVE BEEN SET IN SETUP$READSTRCIN                ==
!     ==  NSP                                                                 ==
!     ==  THIS%I,ID,GID,AEZ,ZV,M,RAD,RCSM,LNX,LOX,LMNX,LMRX                   ==
!     ==  THIS%PARMS%ID,TYPE,RCL,LAMBDA                                       ==
!     ==  THIS%PARMS%POW_POT ,RC_POT ,TVAL0_POT ,VAL0_POT                     ==
!     ==  THIS%PARMS%POW_CORE,RC_CORE,TVAL0_CORE,VAL0_CORE                    ==
!     ==========================================================================
      THIS%SOFTCORETYPE='NONE'   ! OBSOLETE
!
      AEZ=THIS%AEZ
      GID=THIS%GID
      ZV=THIS%ZV
      LNX=THIS%LNX
      ALLOCATE(LOX(LNX))
      LOX=THIS%LOX      
      LX=MAXVAL(LOX)
      RBOX=THIS%RBOX
!
!     ==========================================================================
!     == EXTRACT LNX,LOX,LX FROM NPRO AS DEFINED IN STRC INPUT FILE           ==
!     ==========================================================================
      LMNXX=MAX(LMNXX,THIS%LMNX)
      LMRXX=MAX(LMRXX,THIS%LMRX)
      LNXX=MAX(LNXX,THIS%LNX)
!
!     ==========================================================================
!     == READ SETUP INFORMATION FROM PARAMETER FILE                           ==
!     ==========================================================================
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      ROUT=R(NR-3)
!     
!     ==========================================================================
!     ==  PERFORM ALL-ELECTRON CALCULATION FOR THE ATOM IN A BOX              ==
!     ==========================================================================
      ALLOCATE(THIS%ATOM%AEPOT(NR))
      ALLOCATE(PSI(NR,NBX))
      ALLOCATE(PSISM(NR,NBX))
!
!     == SET SWITCHES FOR RELATIVISTIC/NON-RELATIVISTIC HERE ===================
      THIS%SETTING%TREL=.TRUE.
!THIS%SETTING%TREL=.FALSE.
!     __SPIN ORBIT SWITCH IS READ FROM !AUGMENT BLOCK. DEFAULT=.FALSE.__________
!THIS%SETTING%SO=.FALSE.
PRINT*,'THIS%SETTING%SO=',THIS%SETTING%SO
!THIS%SETTING%SO=.TRUE.
      THIS%SETTING%ZORA=.TRUE.
      IF(THIS%SETTING%SO)THIS%SETTING%ZORA=.FALSE.
!THIS%SETTING%ZORA=.FALSE.
!     == SELECT HARTREE FOCK ADMIXTURE =========================================
      THIS%SETTING%FOCK=0.D0
!
      CALL LMTO$GETL4('ON',TCHK)
      IF(TCHK) THEN
        CALL LMTO$SETI4('ISP',THIS%I)
        CALL LMTO$GETL4('FOCKSETUP',TCHK1)
        IF(TCHK1) THEN
          CALL LMTO$GETR8('LHFWEIGHT',THIS%SETTING%FOCK)
        ELSE
          THIS%SETTING%FOCK=0.D0
        END IF
        CALL LMTO$SETI4('ISP',0)
      END IF
!
      IF(THIS%SETTING%TREL) THEN
        IF(THIS%SETTING%SO) THEN
          KEY='START,REL,SO'
        ELSE 
          KEY='START,REL,NONSO'
        END IF
        IF(THIS%SETTING%ZORA) THEN
          KEY=TRIM(ADJUSTL(KEY))//',ZORA'
        ELSE
          KEY=TRIM(ADJUSTL(KEY))//',NONZORA'
        END IF
      ELSE 
         KEY='START,NONREL,NONSO,NONZORA'
      END IF
!
      IF(THIS%SETTING%FOCK.NE.0.D0) THEN
        WRITE(STRING,*)THIS%SETTING%FOCK
!PRINT*,'FOCK SWITCHED OFF'
        KEY=TRIM(KEY)//',FOCK='//TRIM(STRING)
      END IF
!
      CALL TRACE$PASS('BEFORE SCF-ATOM')
      CALL TIMING$CLOCKON('SCF-ATOM')
!
!     == THESE VARIABLES ARE INTENT(INOUT) BUT NOT YET DEFINED =================
      LOFI=-1111
      SOFI=-1111
      FOFI=0.D0
      NNOFI=-1111
!
!     == PERFORM ALL-ELECTRON SELF-CONSISTENT CALCULATION OF THE ATOM
      CALL ATOMLIB$AESCF(GID,NR,KEY,ROUT,AEZ,NBX,NB,LOFI,SOFI,FOFI,NNOFI &
    &                   ,ETOT,THIS%ATOM%AEPOT,VFOCK,EOFI,PSI,PSISM)
      CALL TIMING$CLOCKOFF('SCF-ATOM')
      CALL TRACE$PASS('AFTER SCF-ATOM')
!     
!     ==========================================================================
!     ==  SELECT CORE, AND REORDER STATES                                     ==
!     ==========================================================================
      ALLOCATE(TC(NB))
      CALL SETUP_CORESELECT(AEZ,THIS%COREID,NB,LOFI,SOFI,TC)
!     
!     ==========================================================================
!     ==  MAP DATA OF THE ATOMIC CALCULATION ONTO DATA STRUCTURE "THIS"       ==
!     ==  NC CORE STATES FIRST THEN THE VALENCE STATES                        ==
!     ==========================================================================
      NC=0
      DO IB=1,NB
        IF(TC(IB))NC=NC+1
      ENDDO
      THIS%ATOM%NC=NC
!     == MAP ATOMIC DATA ON GRID ===============================================
      THIS%ATOM%NB=NB
      ALLOCATE(THIS%ATOM%LOFI(NB))
      ALLOCATE(THIS%ATOM%FOFI(NB))
      ALLOCATE(THIS%ATOM%SOFI(NB))
      ALLOCATE(THIS%ATOM%NNOFI(NB))
      ALLOCATE(THIS%ATOM%EOFI(NB))
      ALLOCATE(THIS%ATOM%AEPSI(NR,NB))
      ALLOCATE(THIS%ATOM%AEPSISM(NR,NB))
      IB1=0
      DO IB=1,NB 
        IF(TC(IB)) THEN
          IB1=IB1+1  
          THIS%ATOM%LOFI(IB1)=LOFI(IB)
          THIS%ATOM%SOFI(IB1)=SOFI(IB)
          THIS%ATOM%FOFI(IB1)=FOFI(IB)
          THIS%ATOM%NNOFI(IB1)=NNOFI(IB)
          THIS%ATOM%EOFI(IB1)=EOFI(IB)
          THIS%ATOM%AEPSI(:,IB1)=PSI(:,IB)
          THIS%ATOM%AEPSISM(:,IB1)=PSISM(:,IB)
        END IF
      ENDDO
      DO IB=1,NB 
        IF(.NOT.TC(IB)) THEN
          IB1=IB1+1  
          THIS%ATOM%LOFI(IB1)=LOFI(IB)
          THIS%ATOM%SOFI(IB1)=SOFI(IB)
          THIS%ATOM%FOFI(IB1)=FOFI(IB)
          THIS%ATOM%NNOFI(IB1)=NNOFI(IB)
          THIS%ATOM%EOFI(IB1)=EOFI(IB)
          THIS%ATOM%AEPSI(:,IB1)=PSI(:,IB)
          THIS%ATOM%AEPSISM(:,IB1)=PSISM(:,IB)
        END IF
      ENDDO
      LOFI(:NB) =THIS%ATOM%LOFI
      SOFI(:NB) =THIS%ATOM%SOFI
      FOFI(:NB) =THIS%ATOM%FOFI
      NNOFI(:NB)=THIS%ATOM%NNOFI
      EOFI(:NB) =THIS%ATOM%EOFI
      PSI(:,:NB)=THIS%ATOM%AEPSI
      PSISM(:,:NB)=THIS%ATOM%AEPSISM
      DEALLOCATE(TC)
!
!     == ZV IS OVERWRITTEN BY INFORMATION FROM COREID ==========================
      ZV=AEZ-SUM(THIS%ATOM%FOFI(:NC))  
      THIS%ZV=ZV
!     
!     ==========================================================================
!     == CUT OFF TAIL OF THE CORE STATES. BOUND STATES HAVE A NODE AT ROUT,   ==
!     == BUT MAY STILL HAVE AN EXPONENTIALLY INCREASING TAIL                  ==
!     ==========================================================================
      DO IR=1,NR
        IF(R(IR).LT.ROUT) CYCLE
        PSI(IR:,:NC)=0.D0
        PSISM(IR:,:NC)=0.D0
        EXIT
      ENDDO
      ALLOCATE(AUX(NR))
      DO IB=1,NC
        AUX(:)=R(:)**2*(PSI(:,IB)**2+PSISM(:,IB)**2)
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        SVAR=1.D0/SQRT(SVAR)
        PSI(:,IB)=PSI(:,IB)*SVAR
        PSISM(:,IB)=PSISM(:,IB)*SVAR
        THIS%ATOM%AEPSI(:,IB)=PSI(:,IB)
        THIS%ATOM%AEPSISM(:,IB)=PSISM(:,IB)
      ENDDO
      DEALLOCATE(AUX)
!     
!     ==========================================================================
!     == REPORT ENERGIES                                                      ==
!     ==========================================================================
      WRITE(6,FMT='(56("="))')
      WRITE(6,FMT='(56("="),T10," ALL-ELECTRON ATOM CALCULATION ")')
      WRITE(6,FMT='(56("="))')
      WRITE(6,FMT='("TOTAL ENERGY=",F25.7)')ETOT
      WRITE(6,FMT='(4A4,A20,2A10)')"IB","N","L","SO","E","F","#(REM. EL.)"
      SVAR=0.D0
      DO IB=1,NB
        SVAR=SVAR+FOFI(IB)  
        WRITE(6,FMT='(4I4,F20.4,2F10.5)')IB,NNOFI(IB)+LOFI(IB)+1,LOFI(IB)&
     &                                  ,SOFI(IB),EOFI(IB),FOFI(IB),AEZ-SVAR
        IF(IB.EQ.NC)WRITE(6,FMT='(56("-"),T10," CORE ABOVE, VALENCE BELOW ")')
      ENDDO
      WRITE(6,FMT='(56("="))')
      WRITE(6,FMT='(56("="))')
!
!     ==========================================================================
!     == CALCULATE AND PSEUDIZE CORE DENSITY                                  ==
!     ==========================================================================
      ALLOCATE(THIS%AECORE(NR))
      ALLOCATE(THIS%PSCORE(NR))
      THIS%AECORE(:)=0.D0
      DO IB=1,NC
        THIS%AECORE(:)=THIS%AECORE(:)+FOFI(IB)*(PSI(:,IB)**2 &
     &                                         +PSISM(:,IB)**2)*C0LL
      ENDDO
      CALL ATOMIC_PSEUDIZE(GID,NR,THIS%PARMS%POW_CORE,THIS%PARMS%TVAL0_CORE &
     &         ,THIS%PARMS%VAL0_CORE,THIS%PARMS%RC_CORE,THIS%AECORE,THIS%PSCORE)
      DEALLOCATE(PSI)
!
!     ==========================================================================
!     == IDENTIFY SCATTERING CHARACTER OF THE PARTIAL WAVES                   ==
!     == (BAND INDEX RELATIVE TO HOMO FOR EACH L)                             ==
!     ==========================================================================
      ALLOCATE(THIS%ISCATT(LNX))
      CALL SETUP_MAKEISCATT(AEZ,NB,NC,LOFI(1:NB),NNOFI(1:NB),LNX,LOX &
     &                                                             ,THIS%ISCATT)
WRITE(6,*)'ISCATT=',THIS%ISCATT
WRITE(6,*)'LOX=   ',LOX

!
!     ==========================================================================
!     == CONSTRUCT PARTIAL WAVES                                              ==
!     ==========================================================================
      ALLOCATE(THIS%VADD(NR))
      ALLOCATE(THIS%PSPOT(NR))
      ALLOCATE(THIS%PRO(NR,LNX))
      ALLOCATE(THIS%AEPHI(NR,LNX))
      ALLOCATE(THIS%AEPHISM(NR,LNX))
      ALLOCATE(THIS%PSPHI(NR,LNX))
      ALLOCATE(THIS%PSPHISM(NR,LNX))
      ALLOCATE(THIS%UPHI(NR,LNX))
      ALLOCATE(THIS%UPHISM(NR,LNX))
      ALLOCATE(THIS%QPHI(NR,LNX))
      ALLOCATE(THIS%QPHISM(NR,LNX))
      ALLOCATE(THIS%DTKIN(LNX,LNX))
      ALLOCATE(THIS%DOVER(LNX,LNX))
      ALLOCATE(THIS%DATH(LNX,LNX))
      ALLOCATE(THIS%PSPHIDOT(NR,LNX))
      ALLOCATE(THIS%PSPHIDOTSM(NR,LNX))
      ALLOCATE(THIS%AEPHIDOT(NR,LNX))
      ALLOCATE(THIS%AEPHIDOTSM(NR,LNX))
      THIS%VADD=0.D0
      THIS%PRO=0.D0
      THIS%AEPHI=0.D0
      THIS%PSPHI=0.D0
!
      ALLOCATE(RC(LNX))
      ALLOCATE(LAMBDA(LNX))
      DO LN=1,LNX
        L=LOX(LN)
        RC(LN)=THIS%PARMS%RCL(L+1)
        LAMBDA(LN)=THIS%PARMS%LAMBDA(L+1)
      ENDDO
!
      CALL TIMING$CLOCKON('MAKEPARTIALWAVES')
      CALL SETUP_MAKEPARTIALWAVES(GID,NR,KEY,LL_STP,AEZ,THIS%ATOM%AEPOT &
     &          ,VFOCK,NB,NC &
     &          ,LOFI(1:NB),SOFI(1:NB),NNOFI(1:NB),EOFI(1:NB),FOFI(1:NB) &
     &          ,RBOX,ROUT,LNX,LOX,THIS%PARMS%TYPE,RC,LAMBDA &
     &          ,THIS%ISCATT &
     &          ,THIS%AEPHI,THIS%AEPHISM,THIS%PSPHI,THIS%PSPHISM &
     &          ,THIS%UPHI,THIS%UPHISM,THIS%QPHI,THIS%QPHISM &
     &          ,THIS%AEPHIDOT,THIS%AEPHIDOTSM,THIS%PSPHIDOT,THIS%PSPHIDOTSM &
     &          ,THIS%PRO,THIS%DTKIN,THIS%DOVER,THIS%DATH &
     &          ,THIS%AECORE,THIS%PSCORE &
     &          ,THIS%PSPOT,THIS%PARMS%POW_POT,THIS%PARMS%TVAL0_POT &
     &          ,THIS%PARMS%VAL0_POT,THIS%PARMS%RC_POT &
     &          ,THIS%RCSM,THIS%VADD,THIS%PSG2,THIS%PSG4)
      ALLOCATE(THIS%NLPHIDOT(NR,LNX))
      ALLOCATE(THIS%QPHIDOT(NR,LNX))
!     == FOR THE NODELESS CONSTRUCTION THE NLPHIDOT, QNPHIDOT AND PSPHIDOT  ====
!     == FUNCTIONS ARE IDENTICAL PER DEFINITION.
      THIS%NLPHIDOT=THIS%PSPHIDOT 
      THIS%QPHIDOT=THIS%PSPHIDOT
CALL TRACE$PASS('AFTER MAKEPARTIALWAVES')

      CALL TIMING$CLOCKOFF('MAKEPARTIALWAVES')
      IF(THIS%SETTING%FOCK.NE.0.D0) THEN
        CALL RADIALFOCK$CLEANVFOCK(VFOCK)
      END IF
!     
!     ==========================================================================
!     ==  DETERMINE CORE-VALENCE EXCHANGE MATRIX ELEMENTS                     ==
!     ==========================================================================
      CALL SETUP_CVXSETUP()
!     
!     ==========================================================================
!     ==  DETERMINE MATRIX ELEMENTS BETWEEN PROJECTOR FUNCTIONS               ==
!     ==  AND ENERGY DERIVATIVE PSEUDO PARTIAL WAVES
!     ==========================================================================
      ALLOCATE(THIS%PROPHIDOT(LNX,LNX))
      CALL SETUP$PROPHIDOT(GID,NR,LNX,THIS%LOX,THIS%PRO,THIS%PSPHIDOT &
     &                     ,THIS%PROPHIDOT)
!!$!     
!!$!     =======================================================================
!!$!     ==  THIS FOR READING AND WRITING THE CORE DENSITY                    ==
!!$!     =======================================================================
!!$      OPEN(UNIT=11,FILE='CORE.DATA',FORM='UNFORMATTED')
!!$      REWIND 11
!!$!     ==  USE THIS FOR WRITEING ...
!!$      WRITE(11)THIS%AECORE      
!!$      WRITE(11)THIS%PSCORE      
!!$!     == ..AND THIS FOR READING
!!$      READ(11)THIS%AECORE      
!!$      READ(11)THIS%PSCORE      
!!$      CLOSE(11)
!     
!     ==========================================================================
!     ==  CALCULATE AND PRINT SCATTERING PROPERTIES                           ==
!     ==========================================================================
      CALL TIMING$CLOCKON('TEST SCATTERING')
      CALL SETUP_TESTSCATTERING(LL_STP)
      CALL TIMING$CLOCKOFF('TEST SCATTERING')
!
      CALL TIMING$CLOCKON('TEST GHOSTS')
      CALL SETUP_TESTGHOST1()
      CALL SETUP_TESTGHOST()
      CALL TIMING$CLOCKOFF('TEST GHOSTS')
!     
!     ==========================================================================
!     ==  PERFORM BESSELTRANSFORMS                                            ==
!     ==========================================================================
      CALL TIMING$CLOCKON('BESSELTRANSFORMS')
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      GIDG=GIDG_PROTO
      THIS%GIDG=GIDG
      CALL RADIAL$GETI4(GIDG,'NR',NG)
!       
!     == VADD (VBAR) ===========================================================
      ALLOCATE(THIS%VADDOFG(NG))
      CALL RADIAL$BESSELTRANSFORM(0,GID,NR,THIS%VADD,GIDG,NG,THIS%VADDOFG)
      THIS%VADDOFG(:)=FOURPI*THIS%VADDOFG(:)
!     == PSCORE (VBAR) =========================================================
      ALLOCATE(THIS%PSCOREOFG(NG))
      CALL RADIAL$BESSELTRANSFORM(0,GID,NR,THIS%PSCORE,GIDG,NG,THIS%PSCOREOFG)
      THIS%PSCOREOFG(:)=FOURPI*THIS%PSCOREOFG(:)
!     == PROJECTORS ============================================================
      ALLOCATE(THIS%PROOFG(NG,LNX))
      DO LN=1,LNX
        L=THIS%LOX(LN)
        CALL RADIAL$BESSELTRANSFORM(L,GID,NR,THIS%PRO(:,LN),GIDG,NG &
     &                             ,THIS%PROOFG(:,LN))
        THIS%PROOFG(:,LN)=FOURPI*THIS%PROOFG(:,LN)
      ENDDO
!     == COMPENSATION GAUSSIAN =================================================
      ALLOCATE(THIS%NHATPRIMEOFG(NG))
      ALLOCATE(THIS%VHATOFG(NG))
      CALL SETUP_COMPOFG(THIS%RCBG,THIS%RCSM,GIDG,NG &
     &                  ,THIS%NHATPRIMEOFG,THIS%VHATOFG)
      CALL TIMING$CLOCKOFF('BESSELTRANSFORMS')
      CALL TIMING$CLOCKOFF('SETUP CONSTRUCTION')
!
!     ==========================================================================
!     == WRITE REPORT
!     ==========================================================================
      CALL SETUP_REPORTFILE(LL_STP)
      WRITE(*,'(80("="))')
      WRITE(*,'(80("="),T20,"  SETUP CONSTRUCTION FINISHED  ")')
      WRITE(*,'(80("="))')
                            CALL TRACE$POP
      RETURN
      END
!
!  
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_REPORTFILE(LL_STP)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE SETUP_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      CHARACTER(64)               :: STRING
      INTEGER(4)                  :: GID
      INTEGER(4)                  :: NR
      INTEGER(4)                  :: NG
      INTEGER(4)                  :: GIDG
      REAL(8)       ,ALLOCATABLE  :: R(:)
      INTEGER(4)                  :: NTASKS,THISTASK
      INTEGER(4)                  :: NFIL
      INTEGER(4)                  :: I
!     **************************************************************************
                                  CALL TRACE$PUSH('SETUP_REPORTFILE')
!     == WRITE REPORT ==========================================================
      CALL LINKEDLIST$SELECT(LL_STP,'~',1)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'GENERIC',0)
      CALL LINKEDLIST$SET(LL_STP,'ID',0,TRIM(THIS%PARMS%ID))
      CALL LINKEDLIST$SET(LL_STP,'Z',0,THIS%AEZ)
      CALL LINKEDLIST$SET(LL_STP,'ZV',0,THIS%ZV)
!
!     ==========================================================================
!     == RADIAL GRID                                                          ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STP,'~',1)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'RGRID',0)
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      CALL LINKEDLIST$SET(LL_STP,'NR',0,NR)
      CALL LINKEDLIST$SET(LL_STP,'R',0,R)
      DEALLOCATE(R)
!
!     ==========================================================================
!     == ATOM                                                                 ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STP,'~',1)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'ATOM',0)
      CALL LINKEDLIST$SET(LL_STP,'NB',0,THIS%ATOM%NB)
      CALL LINKEDLIST$SET(LL_STP,'NC',0,THIS%ATOM%NC)
      CALL LINKEDLIST$SET(LL_STP,'L',0,THIS%ATOM%LOFI)
      CALL LINKEDLIST$SET(LL_STP,'SO',0,THIS%ATOM%SOFI)
      CALL LINKEDLIST$SET(LL_STP,'E',0,THIS%ATOM%EOFI)
      CALL LINKEDLIST$SET(LL_STP,'F',0,THIS%ATOM%FOFI)
      DO I=1,THIS%ATOM%NB
        CALL LINKEDLIST$SET(LL_STP,'AEPSI',-1,THIS%ATOM%AEPSI(:,I))
        CALL LINKEDLIST$SET(LL_STP,'AEPSISM',-1,THIS%ATOM%AEPSISM(:,I))
      ENDDO
      CALL LINKEDLIST$SET(LL_STP,'AEPOT',-1,THIS%ATOM%AEPOT)
      CALL LINKEDLIST$SET(LL_STP,'PSPOT',0,THIS%PSPOT)
      CALL LINKEDLIST$SET(LL_STP,'VADD',0,THIS%VADD)
!
!     ==========================================================================
!     == SETUP PARAMETERS                                                     ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STP,'~',1)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'PARAMETERS',0)
      CALL LINKEDLIST$SET(LL_STP,'RBOX',0,THIS%RBOX)
      CALL LINKEDLIST$SET(LL_STP,'RCSM',0,THIS%RCSM)
      CALL LINKEDLIST$SELECT(LL_STP,'PSI',0)
      CALL LINKEDLIST$SET(LL_STP,'TYPE',0,THIS%PARMS%TYPE)
      CALL LINKEDLIST$SET(LL_STP,'RCL',0,THIS%PARMS%RCL)
      CALL LINKEDLIST$SET(LL_STP,'LAMBDA',0,THIS%PARMS%LAMBDA)
      CALL LINKEDLIST$SELECT(LL_STP,'..',0)
      CALL LINKEDLIST$SELECT(LL_STP,'CORE',0)
      CALL LINKEDLIST$SET(LL_STP,'RC',0,THIS%PARMS%RC_CORE)
      CALL LINKEDLIST$SET(LL_STP,'POW',0,THIS%PARMS%POW_CORE)
      IF(THIS%PARMS%TVAL0_CORE) THEN
        CALL LINKEDLIST$SET(LL_STP,'VAL0',0,THIS%PARMS%VAL0_CORE)
      END IF
      CALL LINKEDLIST$SELECT(LL_STP,'..',0)
      CALL LINKEDLIST$SELECT(LL_STP,'POT',0)
      CALL LINKEDLIST$SET(LL_STP,'RC',0,THIS%PARMS%RC_POT)
      CALL LINKEDLIST$SET(LL_STP,'POW',0,THIS%PARMS%POW_POT)
      IF(THIS%PARMS%TVAL0_POT) THEN
        CALL LINKEDLIST$SET(LL_STP,'VAL0',0,THIS%PARMS%VAL0_POT)
      END IF
!
!     ==========================================================================
!     == AUGMENTATION                                                         ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STP,'~',1)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION',0)
      CALL LINKEDLIST$SET(LL_STP,'AECORE',0,THIS%AECORE)
      CALL LINKEDLIST$SET(LL_STP,'PSCORE',0,THIS%PSCORE)
!
!     ==========================================================================
!     == BESSELTRANSFORMS                                                     ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STP,'~',1)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'BESSELTRANSFORMED',0)
      GIDG=THIS%GIDG
      CALL RADIAL$GETI4(GIDG,'NR',NG)
      ALLOCATE(R(NG))
      CALL RADIAL$R(GIDG,NG,R)
      CALL LINKEDLIST$SET(LL_STP,'NG',0,NG)
      CALL LINKEDLIST$SET(LL_STP,'G',0,R)
      DEALLOCATE(R)
      CALL LINKEDLIST$SET(LL_STP,'PRO',0,THIS%PROOFG)
      CALL LINKEDLIST$SET(LL_STP,'VADD',0,THIS%VADDOFG)
!
!      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
!      CALL LINKEDLIST$REPORT(LL_STP,6)
!
!     ==========================================================================
!     == WRITE SETUP REPORT TO FILE                                           ==
!     ==========================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.EQ.1) THEN
        WRITE(STRING,*)NINT(THIS%AEZ)
        STRING='_STPFORZ'//TRIM(ADJUSTL(STRING))//'.MYXML' 
        CALL FILEHANDLER$SETFILE('STP_REPORT',.TRUE.,-STRING)
        CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','STATUS','UNKNOWN')
        CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION('STP_REPORT','FORM','FORMATTED')
        CALL FILEHANDLER$UNIT('STP_REPORT',NFIL)
        REWIND(NFIL)
        CALL LINKEDLIST$SELECT(LL_STP,'~')
        CALL LINKEDLIST$WRITE(LL_STP,NFIL,'MONOMER')
        CALL FILEHANDLER$CLOSE('STP_REPORT')
      END IF
      CALL TRACE$PASS('SETUP REPORT FILE WRITTEN')
!
!     ==========================================================================
!     == REMOVE LINKED LIST                                                   ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STP,'~')
      CALL LINKEDLIST$RMLIST(LL_STP,'SETUPREPORT')
      CALL LINKEDLIST$SELECT(LL_STP,'~')
!      CALL LINKEDLIST$REPORT(LL_STP,6)
                            CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$PROPHIDOT(GID,NR,LNX,LOX,PRO,PSPHIDOT,PROPHIDOT)
!     **************************************************************************
!     ** DETERMINES <PRO|PSPHIDOT>                                            **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: PRO(NR,LNX)
      REAL(8)   ,INTENT(IN) :: PSPHIDOT(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PROPHIDOT(LNX,LNX)
      REAL(8)               :: R(NR)
      REAL(8)               :: AUX(NR)
      INTEGER(4)            :: LN1,LN2
      LOGICAL   ,PARAMETER  :: TPRINT=.FALSE.
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      PROPHIDOT(:,:)=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN2).NE.LOX(LN1)) CYCLE
          AUX(:)=R(:)**2*PRO(:,LN1)*PSPHIDOT(:,LN2)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,PROPHIDOT(LN1,LN2))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  PRINT FOR TESTING                                                   ==
!     ==========================================================================
      IF(TPRINT) THEN
        PRINT*,'PROPHIDOT'
        DO LN1=1,LNX
          WRITE(*,FMT='(I5,20F20.5)')LOX(LN1),PROPHIDOT(LN1,:)
        ENDDO
        STOP
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_CORESELECT(AEZ,COREID,NB,LOFI,SOFI,TC)
!     **************************************************************************
!     ** SELECT THE STATES THAT SHOULD BE KEPT FROZEN BY TC(IB)               **
!     **                                                                      **
!     ** THE CHOICE OF CORE SHELLS IS STRONGLY CONSTRAINED.                   **
!     ** 1) ZERO CORE STATES                                                  **
!     ** 2) NOBLE GAS CONFIGURATIONS                                          **
!     ** 3) NOBLE GAS CONFIGURATION + COMPLETE D-SHELL                        **
!     ** 4) NOBLE GAS CONFIGURATION + PARTIAL OR COMPLETE F-SHELL             **
!     ** 5) NOBLE GAS CONFIGURATION + COMPLETE F-SHELL + COMPLETE D-SHELL     **
!     **                                                                      **
!     ** CORE HOLES MUST BE TREATED SEPARATELY                                **
!     **                                                                      **
!     ** THE VARIALBLE LMAP(L+1) SPECIFIES THE MAIN QUANTUM NUMBER OF THE     **
!     ** HIGHEST FROZEN CORE STATE FOR EACH ANGULAR MOMENTUM L                **
!     ** (CAUTION! THE SAME VARIABLE NAME IS USED IN OTHER ROUTINES WITH A    **
!     **  DIFFERENT MEANING!)                                                 **
!     **                                                                      **
!     ***************************PETER BLOECHL, GOSLAR 2021*********************
      IMPLICIT NONE
      REAL(8)     ,INTENT(IN) :: AEZ      ! ATOMIC NUMBER
      CHARACTER(*),INTENT(IN) :: COREID
      INTEGER(4)  ,INTENT(IN) :: NB       ! #(ATOMIC SHELLS)
      INTEGER(4)  ,INTENT(IN) :: LOFI(NB) ! ANGULAR MOMENTA
      INTEGER(4)  ,INTENT(IN) :: SOFI(NB) ! SPIN ORBIT INDEX (-1,1) OR 0
      LOGICAL(4)  ,INTENT(OUT):: TC(NB)   ! FROZEN-CORE SELECTOR
      INTEGER(4)              :: LMAP(4)  ! FOR L=0,1,2,3
      INTEGER(4)              :: ISVAR,IB,L,ISO,IL,IPOS
      INTEGER(4)              :: ZCORE
      CHARACTER(1)            :: CH
!     **************************************************************************
      IF(AEZ.GT.118.D0) THEN
        CALL ERROR$MSG('IMPLEMENTATION IS LIMITED TO ATOMIC NUMBER BELOW 119')
        CALL ERROR$R8VAL('ATOMIC NUMBER AEZ',AEZ)
        CALL ERROR$STOP('SETUP_CORESELECT')
      END IF
!
!     ==========================================================================
!     == PARSE COREID                                                         ==
!     == START WITH AN OPTIONAL CORE ID: NOBLE GAS OR "0"                     ==
!     == USE ZERO OR MORE MODIFIERS -S,-P,-D,-F,+S,+P,+D,+F                   ==
!     == "-" REMOVES ONE SHELL WITH SPECIFIED ANGULAR MOMENTUM FROM CORE      ==
!     == "+" ADDS ONE SHELL WITH SPECIFIED ANGULAR MOMENTUM TO THE CORE       ==
!     ==========================================================================
!     == SET DEFAULT LMAP: NEXT LOWER NOBLE GAS CONFIGRATION
      LMAP=(/0,0,0,0/) !EMPTY 
      IF(AEZ.GT.  2.D0) LMAP=(/1,0,0,0/) !HE
      IF(AEZ.GT. 10.D0) LMAP=(/2,1,0,0/) !NE
      IF(AEZ.GT. 18.D0) LMAP=(/3,2,0,0/) !AR
      IF(AEZ.GT. 36.D0) LMAP=(/4,3,1,0/) !KR
      IF(AEZ.GT. 54.D0) LMAP=(/5,4,2,0/) !XE
      IF(AEZ.GT. 86.D0) LMAP=(/6,5,3,1/) !RN


      IPOS=1
      IF(COREID(1:1).EQ.'+'.OR.COREID(1:1).EQ.'-') THEN
!       == DO NOTHING
      ELSE IF(LEN_TRIM(COREID).EQ.0) THEN 
!       == DO NOTHING
      ELSE IF(COREID(1:1).EQ.'0') THEN
        LMAP=(/0,0,0,0/)
        IPOS=IPOS+1
      ELSE IF(COREID(1:2).EQ.'HE') THEN
        LMAP=(/1,0,0,0/)
        IPOS=IPOS+2
      ELSE IF(COREID(1:2).EQ.'NE') THEN
        LMAP=(/2,1,0,0/)
        IPOS=IPOS+2
      ELSE IF(COREID(1:2).EQ.'AR') THEN
        LMAP=(/3,2,0,0/)
        IPOS=IPOS+2
      ELSE IF(COREID(1:2).EQ.'KR') THEN
        LMAP=(/4,3,1,0/)
        IPOS=IPOS+2
      ELSE IF(COREID(1:2).EQ.'XE') THEN
        LMAP=(/5,4,2,0/)
        IPOS=IPOS+2
      ELSE IF(COREID(1:2).EQ.'RN') THEN
        LMAP=(/6,5,3,1/)
        IPOS=IPOS+2
      ELSE
        CALL ERROR$MSG('MAL-FORMED CORE ID')
        CALL ERROR$MSG('FIRST LETTERS MUST SPECIFY NOBLE GAS ELEMENT OR "0"')
        CALL ERROR$MSG('OR THE SIGN, "+" OR "-" OF A MODIFIER')
        CALL ERROR$MSG('EXAMPLE:  KR-S-D+4F')
        CALL ERROR$R8VAL('ATOMIC NUMBER AEZ',AEZ)
        CALL ERROR$CHVAL('COREID',COREID)
        CALL ERROR$STOP('SETUP_CORESELECT')
      END IF
!
!     == EXTRACT AND APPLY MODIFIERS ===========================================
      DO WHILE(LEN_TRIM(COREID(IPOS:)).GT.0)
        CH=COREID(IPOS+1:IPOS+1)
        IF(CH.EQ.'S') THEN
          IL=1
        ELSE IF(CH.EQ.'P') THEN
          IL=2
        ELSE IF(CH.EQ.'D') THEN
          IL=3
        ELSE IF(CH.EQ.'F') THEN
          IL=4
        ELSE
          CALL ERROR$MSG('MAL-FORMED CORE ID (EXAMPLE: KR+6S+4D+4F)')
          CALL ERROR$MSG('ANGULAR MOMENTUM ID NOT RECOGNIZED')
          CALL ERROR$R8VAL('ATOMIC NUMBER AEZ',AEZ)
          CALL ERROR$CHVAL('COREID',COREID)
          CALL ERROR$CHVAL('ANGULAR MOMENTUM ID',CH)
          CALL ERROR$STOP('SETUP_CORESELECT')
        END IF
        IF(COREID(IPOS:IPOS).EQ.'+') THEN
          LMAP(IL)=LMAP(IL)+1
        ELSE IF(COREID(IPOS:IPOS).EQ.'-') THEN
          LMAP(IL)=LMAP(IL)-1
        ELSE
          CALL ERROR$MSG('MAL-FORMED CORE ID (EXAMPLE: KR+6S+4D+4F)')
          CALL ERROR$MSG('SIGN NOT RECOGNIZED')
          CALL ERROR$R8VAL('ATOMIC NUMBER AEZ',AEZ)
          CALL ERROR$CHVAL('COREID',COREID)
          CALL ERROR$CHVAL('SIGN',COREID(IPOS:IPOS))
          CALL ERROR$STOP('SETUP_CORESELECT')
        END IF
        IPOS=IPOS+2
      ENDDO
!
!     == CONSISTENCY CHECK =====================================================
      DO IL=1,4
        IF(IL.LT.0) THEN
          CALL ERROR$MSG('NUMBER OF FROZEN CORE STATES MUST NOT BE NEGATIVE')
          CALL ERROR$R8VAL('ATOMIC NUMBER AEZ',AEZ)
          CALL ERROR$CHVAL('COREID',COREID)
          CALL ERROR$I4VAL('L',IL-1)
          CALL ERROR$I4VAL('#(FROZEN CORE STATES)',LMAP(IL))
          CALL ERROR$STOP('SETUP_CORESELECT')
        END IF
      ENDDO
!
!     ==========================================================================
!     == SET FROZEN-CORE SELECTOR TC                                          ==
!     ==========================================================================
      TC=.FALSE.
      DO L=0,3
        DO ISO=-1,1
          ISVAR=0
          DO IB=1,NB
            IF(LOFI(IB).NE.L) CYCLE
            IF(SOFI(IB).NE.ISO) CYCLE
            ISVAR=ISVAR+1
            TC(IB)=(ISVAR.LE.LMAP(L+1)) 
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_CORESELECT_OLD(AEZ,ZV,NB,LOFI,SOFI,TC)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: AEZ
      REAL(8)   ,INTENT(IN)  :: ZV
      INTEGER(4),INTENT(IN)  :: NB
      INTEGER(4),INTENT(IN)  :: LOFI(NB)
      INTEGER(4),INTENT(IN)  :: SOFI(NB)
      LOGICAL(4),INTENT(OUT) :: TC(NB)
      INTEGER(4)             :: ZNOBLE(0:7)=(/0,2,10,18,36,54,86,118/)
      INTEGER(4)             :: IS(0:15)=(/0,1,0,0,1,0,1,0,0,1,0,1,0,1,0,1/)
      INTEGER(4)             :: IP(0:15)=(/0,0,0,1,1,0,0,0,1,1,1,1,0,0,1,1/)
      INTEGER(4)             :: ID(0:15)=(/0,0,0,0,0,1,1,0,1,1,0,0,1,1,1,1/)
      INTEGER(4)             :: IF(0:15)=(/0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,1/)
      INTEGER(4)           :: LMAP(19)=(/0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1/)
      INTEGER(4)             :: I,IB,L,ISO,IB0,ISVAR,NEL
      INTEGER(4)             :: IZ1,IZ2
      INTEGER(4)             :: ION(0:3)
      INTEGER(4)             :: ISUMEL
      LOGICAL(4)             :: TMAP(19)
!     **************************************************************************
!     == DETERMINE ATOMIC NUMBER OF LAST NOBEL GAS ATOM BEFORE VALENCE STATES ==
      DO I=1,7
        IZ1=ZNOBLE(I-1)
        IZ2=ZNOBLE(I)
        IF(INT(AEZ-ZV).LT.IZ2) EXIT
      ENDDO
!
!     == DETERMINE OCCUPIED CORE SHELLS ========================================
      ISVAR=NINT(AEZ-ZV)-IZ1   
      ISVAR=ISVAR/2
      IF(ISVAR.LT.0.OR.ISVAR.GT.15) THEN
        CALL ERROR$MSG('LESS THAN 0 OR MORE THAN 30 CORE ELECTRONS REQUESTED')
        CALL ERROR$R8VAL('AEZ',AEZ)
        CALL ERROR$R8VAL('ZV',ZV)
        CALL ERROR$I4VAL('IZ1',IZ1)
        CALL ERROR$I4VAL('IZ2',IZ2)
        CALL ERROR$I4VAL('ISVAR',ISVAR)
        CALL ERROR$STOP('SETUP_CORESELECT')
      END IF
      ION(0)=IS(ISVAR)
      ION(1)=IP(ISVAR)
      ION(2)=ID(ISVAR)
      ION(3)=IF(ISVAR)
!      
!     ==  FILL COMPLETE CORE SHELLS ============================================
      TMAP=.FALSE.
      ISUMEL=IZ1
      IB0=1
      DO IB=1,19
        NEL=2*(2*LMAP(IB)+1)
        TMAP(IB)=ISUMEL.GE.NEL
        IF(TMAP(IB)) ISUMEL=ISUMEL-NEL
        IF(TMAP(IB)) IB0=IB+1  !IB0 POINTS TO FIRST VALENCE SHELL
      ENDDO
      IF(ISUMEL.NE.0) THEN
        CALL ERROR$MSG('INTERNAL ERROR 1')
        CALL ERROR$STOP('SETUP_CORESELECT')
      END IF
!
!     == FILL REMAINING SUB-SHELLS
      DO IB=IB0,19
        IF(ION(LMAP(IB)).EQ.1) THEN
          TMAP(IB)=.TRUE.
          ION(LMAP(IB))=0
        END IF
      ENDDO
!
!     == IDENTIFY CORE SHELLS ==================================================
      TC=.FALSE.
      DO L=0,3
        DO ISO=-1,1
!         == COUNT NUMBER OF CORE SHELLS WITH THIS L
          ISVAR=0
          DO I=1,19
            IF(LMAP(I).NE.L) CYCLE
            IF(TMAP(I))ISVAR=ISVAR+1
          ENDDO
!         == SWITCH CORE STATES ON
          DO IB=1,NB
            IF(LOFI(IB).NE.L) CYCLE
            IF(SOFI(IB).NE.ISO) CYCLE
            IF(ISVAR.EQ.0) EXIT
            TC(IB)=.TRUE.
            ISVAR=ISVAR-1
          ENDDO
        ENDDO
      ENDDO     
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$REPORT(NFIL)
!     **************************************************************************
!     **  REPORT SETUP INFORMATION                                            **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: NFIL
      REAL(8)        ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)        ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4)                 :: L,NPRO,LN,NPROSUM
      TYPE(THIS_TYPE),POINTER    :: THIS1
      CHARACTER(64)              :: STRING
      REAL(8)                    :: U
      INTEGER(4)                 :: THISTASK,NTASKS
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      CALL CONSTANTS$GET('U',U)
!
      IF(.NOT.ASSOCIATED(FIRST)) THEN
        CALL REPORT$TITLE(NFIL,'ATOMIC SETUP: '//'NO SETUPS SPECIFIED')
        RETURN
      END IF
      THIS1=>FIRST
      DO 
        WRITE(NFIL,*)
        CALL REPORT$TITLE(NFIL,'ATOMIC SETUP: '//TRIM(THIS1%ID))
        CALL REPORT$R8VAL(NFIL,'ATOMIC NUMBER',THIS1%AEZ,' ')
        CALL REPORT$R8VAL(NFIL,'ATOMIC MASS  ',THIS1%M/U,'U')
        CALL REPORT$R8VAL(NFIL,'#(VALENCE ELECTRONS)',THIS1%ZV,' ')
        L=0
        NPROSUM=0
        DO WHILE (NPROSUM.LT.THIS1%LNX)
          NPRO=0
          DO LN=1,THIS1%LNX
            IF(THIS1%LOX(LN).EQ.L) NPRO=NPRO+1
          ENDDO
          IF(NPRO.NE.0) THEN
             WRITE(STRING,*)L
             STRING='#(PROJECTORS FOR L='//TRIM(ADJUSTL(STRING))
             STRING=TRIM(STRING)//')'
             CALL REPORT$I4VAL(NFIL,TRIM(STRING),NPRO,' ')
          END IF
          L=L+1
          NPROSUM=NPROSUM+NPRO
        ENDDO        

        CALL REPORT$I4VAL(NFIL,'MAX #(ANGULAR MOMENTA (L,M) FOR 1C-DENSITY)' &
     &                        ,THIS1%LMRX,' ')
        CALL REPORT$R8VAL(NFIL,'GAUSSIAN DECAY FOR COMPENSATION DENSITY ' &
     &                        ,THIS1%RCSM,'ABOHR ')
        CALL REPORT$R8VAL(NFIL,'GAUSSIAN DECAY FOR EXTENDED COMPENSATION' &
     &                        //' DENSITY',THIS1%RCBG,'ABOHR ')
        CALL REPORT$R8VAL(NFIL,'PARAMETER PS<G2> FOR MASS RENORMALIZATION' &
     &                        ,THIS1%PSG2,'H')
        CALL REPORT$R8VAL(NFIL,'PARAMETER PS<G4> FOR MASS RENORMALIZATION' &
     &                        ,THIS1%PSG4,'H')
!
!       == REPORT SETTINGS FOR THE SETUP CONSTRUCTION ==========================
        CALL REPORT$CHVAL(NFIL,'SETUP ID',THIS1%PARMS%ID)
        IF(THIS1%PARMS%TYPE.EQ.'KERKER') THEN
          CALL REPORT$CHVAL(NFIL,'CONSTRUCTION METHOD','KERKER TYPE')
        ELSE IF(THIS1%PARMS%TYPE.EQ.'HBS') THEN
          CALL REPORT$CHVAL(NFIL,'CONSTRUCTION METHOD' &
     &                          ,'HAMANN-BACHELET-SCHLUETER TYPE')
        ELSE IF(THIS1%PARMS%TYPE.EQ.'NDLSS') THEN
          CALL REPORT$CHVAL(NFIL,'CONSTRUCTION METHOD' &
     &                          ,'NODELESS ORBITALS')
        ELSE
          CALL REPORT$CHVAL(NFIL,'CONSTRUCTION METHOD',THIS1%PARMS%TYPE)
        END IF
        DO L=0,MAXVAL(THIS1%LOX)
          WRITE(STRING,*)L
          STRING='PARTIAL WAVE PSEUDIZATION PARAMETER RC FOR L=' &
     &          //TRIM(ADJUSTL(STRING))
          CALL REPORT$R8VAL(NFIL,TRIM(STRING),THIS1%PARMS%RCL(L+1),'ABOHR')
          IF(THIS1%PARMS%TYPE.EQ.'HBS') THEN
            WRITE(STRING,*)L
            STRING='PARTIAL WAVE PSEUDIZATION PARAMETER LAMBDA FOR L=' &
     &                                           //TRIM(ADJUSTL(STRING))
            CALL REPORT$R8VAL(NFIL,TRIM(STRING),THIS1%PARMS%LAMBDA(L+1) &
     &                                         ,'ABOHR')
          END IF 
        ENDDO
        IF(THIS1%PARMS%TVAL0_POT) THEN
          CALL REPORT$R8VAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: VAL(0)' &
     &                          ,THIS1%PARMS%VAL0_POT*Y0,'H')
        ELSE
          CALL REPORT$CHVAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: VAL(0)' &
     &                          ,'NOT DETERMINED')
        END IF
        CALL REPORT$R8VAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: RC' &
     &                          ,THIS1%PARMS%RC_POT,'ABOHR')
        CALL REPORT$R8VAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: POWER' &
     &                          ,THIS1%PARMS%POW_POT,' ')
        IF(THIS1%PARMS%TVAL0_CORE) THEN
          CALL REPORT$R8VAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER:' &
     &                          //' VAL(0)',THIS1%PARMS%VAL0_CORE*Y0,'A_0^(-3)')
        ELSE
          CALL REPORT$CHVAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER:' & 
     &                          //' VAL(0)','NOT DETERMINED')
        END IF
        CALL REPORT$R8VAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER: RC' &
     &                         ,THIS1%PARMS%RC_CORE,'ABOHR')
        CALL REPORT$R8VAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER: POWER' &
     &                         ,THIS1%PARMS%POW_CORE,' ')
        IF(.NOT.ASSOCIATED(THIS1%NEXT)) EXIT
        THIS1=>THIS1%NEXT
      ENDDO
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP$READ()
!     **************************************************************************
!     **  WORK OUT SETUP INFORMATION FROM THE PARAMETERS, WHICH HAVE ALREADY  **
!     **  BEEN SET                                                            **
!     **************************************************************************
      USE SETUP_MODULE, ONLY : NSP
      IMPLICIT NONE
      INTEGER(4)              :: ISP
!     **************************************************************************
                            CALL TRACE$PUSH('SETUP$READ')
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP_READ
        CALL SETUP$UNSELECT()
      ENDDO
                            CALL TRACE$POP
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE SETUP_COMPOFG(RCBG,RCSM,GIDG,NG,G0,V0)
!     **************************************************************************
!     **                                                                      **
!     **  COMPENSATION DENSITY AND POTENTIAL ON A RADIAL GRID                 **
!     **  IN G-SPACE                                                          **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RCBG
      REAL(8)   ,INTENT(IN) :: RCSM
      INTEGER(4),INTENT(IN) :: GIDG
      INTEGER(4),INTENT(IN) :: NG
      REAL(8)   ,INTENT(OUT):: G0(NG)
      REAL(8)   ,INTENT(OUT):: V0(NG)
      REAL(8)   ,PARAMETER  :: EPSILONGAMMA=1.D-7
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: SVAR1,SVAR2,SVAR3,SVAR4
      REAL(8)               :: BGGAUSS,SMGAUSS
      INTEGER(4)            :: IG
      REAL(8)               :: GARR(NG)
!     ******************************************************************
      SVAR1=-0.25D0*RCBG**2
      SVAR2=-0.25D0*RCSM**2
      SVAR3=4.D0*PI
      SVAR4=-4.D0*PI*SVAR3
      CALL RADIAL$R(GIDG,NG,GARR)
      DO IG=1,NG
        BGGAUSS=EXP(SVAR1*GARR(IG)**2)
        SMGAUSS=EXP(SVAR2*GARR(IG)**2)
        G0(IG)=SVAR3*BGGAUSS
        V0(IG)=SVAR4*(BGGAUSS-SMGAUSS)/GARR(IG)**2
      ENDDO
      RETURN  
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_MAKEPARTIALWAVES(GID,NR,KEY,LL_STP,AEZ,AEPOT,VFOCK &
     &                  ,NB,NC,LOFI,SOFI,NNOFI,EOFI,FOFI &
     &                  ,RBOX,ROUT,LNX,LOX,TYPE,RC,LAMBDA,ISCATT &
     &                  ,AEPHI,AEPHISM,PSPHI,PSPHISM,NLPHI,NLPHISM,QN,QNSM &
     &                  ,AEPHIDOT,AEPHIDOTSM,PSPHIDOT,PSPHIDOTSM &
     &                  ,PRO,DT,DOVER,DH,AECORE,PSCORE,PSPOT &
     &                  ,POW_POT,TVAL0_POT,VAL0_POT,RC_POT,RCSM,VADD,PSG2,PSG4)
!     **************************************************************************
!     **  SETUP CONSTRUCTION: PARTIAL WAVES, POTENTIALS, ETC                  **
!     **                                                                      **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      USE PERIODICTABLE_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      USE RADIALFOCK_MODULE,ONLY: VFOCK_TYPE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: TYPE
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      CHARACTER(*),INTENT(IN) :: KEY
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: AEPOT(NR)   ! ALL ELECTRON POTENTIAL
      TYPE(VFOCK_TYPE),INTENT(INOUT) :: VFOCK
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NC
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      INTEGER(4),INTENT(IN) :: SOFI(NB)
      INTEGER(4),INTENT(IN) :: NNOFI(NB)
      REAL(8)   ,INTENT(IN) :: EOFI(NB)
      REAL(8)   ,INTENT(IN) :: ROUT    ! HARD SPHERE FOR ATOM HAS RADIUS RBOX
      REAL(8)   ,INTENT(IN) :: RBOX    ! RADIUS FOR BOUNDARY CONDITION OF PHI
      REAL(8)   ,INTENT(IN) :: FOFI(NB)
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: ISCATT(LNX)
      REAL(8)   ,INTENT(IN) :: RC(LNX)
      REAL(8)   ,INTENT(IN) :: LAMBDA(LNX)
      REAL(8)   ,INTENT(IN) :: POW_POT
      LOGICAL(4),INTENT(IN) :: TVAL0_POT
      REAL(8)   ,INTENT(IN) :: VAL0_POT
      REAL(8)   ,INTENT(IN) :: RC_POT
      REAL(8)   ,INTENT(IN) :: RCSM
      REAL(8)   ,INTENT(IN) :: AECORE(NR)
      REAL(8)   ,INTENT(IN) :: PSCORE(NR)
      REAL(8)   ,INTENT(OUT):: VADD(NR)
      REAL(8)   ,INTENT(OUT):: PSPOT(NR)
      REAL(8)   ,INTENT(OUT):: AEPHI(NR,LNX),AEPHISM(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PSPHI(NR,LNX),PSPHISM(NR,LNX)
      REAL(8)   ,INTENT(OUT):: NLPHI(NR,LNX),NLPHISM(NR,LNX)
      REAL(8)   ,INTENT(OUT):: QN(NR,LNX),QNSM(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PRO(NR,LNX)
      REAL(8)   ,INTENT(OUT):: DT(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: DOVER(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: DH(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: PSPHIDOT(NR,LNX),PSPHIDOTSM(NR,LNX)  
      REAL(8)   ,INTENT(OUT):: AEPHIDOT(NR,LNX),AEPHIDOTSM(NR,LNX)  
      REAL(8)   ,INTENT(OUT):: PSG2
      REAL(8)   ,INTENT(OUT):: PSG4
      REAL(8)   ,PARAMETER  :: TOL=1.D-7
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      LOGICAL   ,PARAMETER  :: TWRITE=.FALSE.
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)   ,PARAMETER  :: C0LL=Y0
      INTEGER(4),ALLOCATABLE:: NCL(:)
      REAL(8)               :: UOFI(NR,NB)   ! NODELESS WAVE FUNCTION
      REAL(8)               :: UOFISM(NR,NB) ! SMALL COMPONENT
      REAL(8)               :: PHISCALE(LNX)
      REAL(8)               :: PSISCALE(NB)
      REAL(8)               :: PSRHO(NR),PAWRHO(NR)
      INTEGER(4)            :: LX
      INTEGER(4)            :: L,IB,LN,IR,LN1,LN2
      REAL(8)               :: R(NR)
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)   ,ALLOCATABLE:: AUXARR(:,:)
      REAL(8)               :: VAL,VAL1,VAL2
      REAL(8)               :: SVAR
      LOGICAL(4)            :: TREL,TSO,TZORA
      REAL(8)   ,ALLOCATABLE:: A(:,:)
      REAL(8)               :: AEPSIF(NR,NB-NC)
      REAL(8)               :: PSPSIF(NR,NB-NC)
      REAL(8)               :: AUGPSIF(NR,NB-NC)
      CHARACTER(64)         :: STRING
      REAL(8)               :: RCOV    !COVALENT RADIUS
      REAL(8)               :: RNORM   !NORMALIZATIONS ARE DONE WITHIN RNORM
      REAL(8)               :: RBNDOUT
      CHARACTER(32)         :: RELTYPE
!     **************************************************************************
                                CALL TRACE$PUSH('SETUP_MAKEPARTIALWAVES')
!VFOCK%SCALE=0.D0
      LX=MAX(MAXVAL(LOX),MAXVAL(LOFI))
      CALL RADIAL$R(GID,NR,R)
      CALL PERIODICTABLE$GET(NINT(AEZ),'R(COV)',RCOV)
      RNORM=RBOX
!     == RBNDOUT HAS EARLIER BEEN SET EQUAL TO RBOX. IT SHOULD BE EQUAL TO ROUT.
!     == THIS VARIABLE IS USED TO TRACK THE CHANGES.
      RBNDOUT=ROUT  
IF(TWRITE) THEN
 PRINT*,'SETUP CONSTRUCTION TYPE ',TYPE
 PRINT*,'R(NR)   ',R(NR),' OUTERMOST GRIDPOINT'
 PRINT*,'R(NR-1) ',R(NR-1),' SECOND TO OUTERMOST GRIDPOINT'
 PRINT*,'R(NR-2) ',R(NR-2),' THIRD TOOUTERMOST GRIDPOINT'
 PRINT*,'ROUT    ',ROUT,' HARD SPHERE RADIUS FOR ATOMIC CALCULATION'
 PRINT*,'RBOX    ',RBOX,' RADIUS FOR BOUNDARY CONDITIONS FOR PARTIAL WAVES'
 PRINT*,'RNORM   ',RNORM,' NORMALIZATION WILL BE DONE UP TO RNORM'
 PRINT*,'RCOV    ',RCOV,' COVALENT RADIUS'
END IF
!
!     ==========================================================================
!     == RESOLVE KEY                                                          ==
!     ==========================================================================
      TREL=INDEX(KEY,'NONREL').EQ.0
      IF(TREL.AND.INDEX(KEY,'REL').EQ.0) THEN
        CALL ERROR$STOP('SETUP_MAKEPARTIALWAVES')
      END IF
      TSO=INDEX(KEY,'NONSO').EQ.0
      IF(TSO.AND.INDEX(KEY,'SO').EQ.0) THEN
        CALL ERROR$STOP('SETUP_MAKEPARTIALWAVES')
      END IF
      TZORA=INDEX(KEY,'NONZORA').EQ.0
      IF(TZORA.AND.INDEX(KEY,'ZORA').EQ.0) THEN
        CALL ERROR$STOP('SETUP_MAKEPARTIALWAVES')
      END IF
!
!     == DETERMINE HIGHEST CORE STATE FOR EACH ANGULAR MOMENTUM ================
      ALLOCATE(NCL(0:LX))
      NCL(:)=0
      DO IB=1,NC
        L=LOFI(IB)
        NCL(L)=MAX(NCL(L),IB)
      ENDDO
!
!     ==========================================================================
!     == CONSTRUCT PSEUDO POTENTIAL                                           ==
!     ==========================================================================
      CALL ATOMIC_PSEUDIZE(GID,NR,POW_POT,TVAL0_POT,VAL0_POT,RC_POT,AEPOT,PSPOT)
!
!     ==========================================================================
!     == ALTERNATE PROJECTOR CONSTRUCTION                                     ==
!     ==========================================================================
      IF(TYPE.EQ.'NDLSS') THEN
        IF(TREL.AND.(.NOT.TZORA).AND.TSO) THEN
          RELTYPE='SPINORBIT'
        ELSE IF(TREL.AND.TZORA.AND.(.NOT.TSO)) THEN
          RELTYPE='ZORA'
        ELSE IF((.NOT.TREL).AND.(.NOT.TZORA).AND.(.NOT.TSO)) THEN
          RELTYPE='NONREL'
        ELSE
          CALL ERROR$MSG('CAN NOT SELECT RELTYPE: UNKNOWN SELECTION')
          CALL ERROR$L4VAL('TREL',TREL)
          CALL ERROR$L4VAL('TZORA',TZORA)
          CALL ERROR$L4VAL('TSO',TSO)
          CALL ERROR$STOP('MAKEPARTIALWAVES')
        END IF
        CALL SETUP_OUTERNEWPROWRAPPER(GID,NR,ROUT,RELTYPE &
      &                   ,NB,NC,LOFI,SOFI,EOFI,LNX,LOX,RC,AEPOT,PSPOT,VFOCK &
      &                   ,QN,QNSM,AEPHI,AEPHISM,PSPHI,PSPHISM,PRO,DT,DOVER &
      &                   ,AEPHIDOT,AEPHIDOTSM,PSPHIDOT,PSPHIDOTSM)
!       == PSPHIDOT,PSPHIDOTSM ARE PER DEFINITION IDENTICAL TO QNDOT,QNDOTSM ===
        NLPHI=QN
        NLPHISM=QNSM
        IF(TWRITE) THEN
          WRITE(STRING,FMT='(F3.0)')AEZ
          STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
          CALL SETUP_WRITEPHI(-'NEW_PRO'//TRIM(STRING),GID,NR,LNX,PRO)
          CALL SETUP_WRITEPHI(-'NEW_QN'//TRIM(STRING),GID,NR,LNX,QN)
          CALL SETUP_WRITEPHI(-'NEW_AEPHI'//TRIM(STRING),GID,NR,LNX,AEPHI)
          CALL SETUP_WRITEPHI(-'NEW_PSPHI'//TRIM(STRING),GID,NR,LNX,PSPHI)
          CALL SETUP_WRITEPHI(-'NEW_AEPHIDOT'//TRIM(STRING),GID,NR,LNX,AEPHIDOT)
          CALL SETUP_WRITEPHI(-'NEW_PSPHIDOT'//TRIM(STRING),GID,NR,LNX,PSPHIDOT)
        END IF
        !
        !  MISSING VARIABLES:
        !
        PHISCALE=1.D0
        PSISCALE=1.D0
      ELSE 
        CALL SETUP_OUTEROLDPROWRAPPER(GID,NR,ROUT,RBOX,RCOV &
     &                    ,NC,NB,LOFI,SOFI,EOFI,LNX,LOX,TYPE,RC,LAMBDA,ISCATT &
     &                    ,AEPOT,PSPOT,VFOCK &
     &                    ,QN,AEPHI,PSPHI,PRO,DT,DOVER &
     &                    ,AEPHIDOT,PSPHIDOT,TREL,TZORA,PSISCALE,PHISCALE)
        NLPHI=QN
        NLPHISM=0.D0
        QNSM=0.D0
        AEPHISM=0.D0
        PSPHISM=0.D0
        AEPHIDOTSM=0.D0
        PSPHIDOTSM=0.D0
        IF(TWRITE) THEN
          WRITE(STRING,FMT='(F3.0)')AEZ
          STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
          CALL SETUP_WRITEPHI(-'OLD_PRO'//TRIM(STRING),GID,NR,LNX,PRO)
          CALL SETUP_WRITEPHI(-'OLD_QN'//TRIM(STRING),GID,NR,LNX,QN)
          CALL SETUP_WRITEPHI(-'OLD_AEPHI'//TRIM(STRING),GID,NR,LNX,AEPHI)
          CALL SETUP_WRITEPHI(-'OLD_PSPHI'//TRIM(STRING),GID,NR,LNX,PSPHI)
          CALL SETUP_WRITEPHI(-'OLD_AEPHIDOT'//TRIM(STRING),GID,NR,LNX,AEPHIDOT)
          CALL SETUP_WRITEPHI(-'OLD_PSPHIDOT'//TRIM(STRING),GID,NR,LNX,PSPHIDOT)
        END IF
      END IF
!
!     ==========================================================================
!     == CALCULATE DH                                                         ==
!     ==========================================================================
      DH=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOX(LN2),AEPHI(:,LN2),AUX1)
          AUX=AEPOT*(AEPHI(:,LN1)*AEPHI(:,LN2)+AEPHISM(:,LN1)*AEPHISM(:,LN2)) &
     &       -PSPOT*(PSPHI(:,LN1)*PSPHI(:,LN2)+PSPHISM(:,LN1)*PSPHISM(:,LN2)) &
     &       +AEPHI(:,LN1)*AUX1(:)/Y0
          AUX=Y0*R**2*AUX  !Y0 FOR THE POTENTIAL
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!         == IN THE OLDER VERSION RBOX IS USED INSTEAD OF RBOX, WHICH IS 
!         == INCONSISTENT. 
!         == ROUT IS THE SIZE OF THE BOX IN WHICH THE ATOM "LIVES", WHILE 
!         == RBOX SELECTS THE "EXCITED" PARTIAL WAVES
          CALL RADIAL$VALUE(GID,NR,AUX1,RBNDOUT,SVAR)
          DH(LN1,LN2)=DT(LN1,LN2)+SVAR
        ENDDO
      ENDDO
!!$DO LN1=1,LNX
!!$  WRITE(*,FMT='(A,100E12.5)')'DT   ',DT(LN1,:)
!!$ENDDO
!!$DO LN1=1,LNX
!!$  WRITE(*,FMT='(A,100E12.5)')'DH   ',DH(LN1,:)
!!$ENDDO
!!$DO LN1=1,LNX
!!$  WRITE(*,FMT='(A,100E12.5)')'DOVER',DOVER(LN1,:)
!!$ENDDO
!
!     ==========================================================================
!     == CALCULATE DENSITY FOR UNSCREENING                                    ==
!     ==========================================================================
                      CALL TRACE$PASS('CONSTRUCT DENSITY OR UNSCREENING')
      CALL SETUP_PAWDENSITY(GID,NR,LNX,LOX,NB,NC,LOFI,NNOFI,EOFI,FOFI &
     &                      ,AECORE,PSCORE &
     &                      ,AEPHI,AEPHISM,PSPHI,PSPHISM,PRO &
     &                      ,DH,DOVER,AEPOT,PSPOT,VFOCK &
     &                      ,ROUT,TREL,TZORA &
     &                      ,PAWRHO,PSRHO,AEPSIF,PSPSIF,AUGPSIF)
!!$PRINT*,'GID,NR(3) ',GID,NR,R(1:3),'...',R(NR)
!!$PRINT*,'ROUT=',ROUT
!!$CALL SETUP_WRITEPHI('AECORE',GID,NR,1,AECORE)
!!$CALL SETUP_WRITEPHI('PAWVALRHO',GID,NR,1,PAWRHO-AECORE)
!!$CALL SETUP_WRITEPHI('PAWRHO',GID,NR,1,PAWRHO)
!!$CALL SETUP_WRITEPHI('PSRHO',GID,NR,1,PSRHO)
!!$STOP 'FORCED IN PAW_SETUPS.F90 AFTER PAWDENSITY'
!
!     ==========================================================================
!     == UNSCREENING                                                          ==
!     ==========================================================================
                      CALL TRACE$PASS('CONSTRUCT VADD (UNSCREEN)')
      CALL ATOMIC_UNSCREEN(GID,NR,ROUT,AEZ,PAWRHO,PSRHO,PSPOT,RCSM,VADD)
      DO IR=1,NR
        IF(R(IR).GT.MAX(1.2D0*RCOV,RNORM)) THEN
          VADD(IR+1:)=0.D0
          EXIT
        END IF
      ENDDO
!
!     ==========================================================================
!     == MAKE DT AND DO SYMMETRIC                                             ==
!     ==========================================================================
      WRITE(6,FMT='(80("="),T20,"  DTKIN BEFORE SYMMETRIZATION ")')
      DO LN1=1,LNX
        WRITE(6,FMT='(20F12.3)')DT(LN1,:)
      ENDDO
      WRITE(6,FMT='(80("="),T20,"  DO BEFORE SYMMETRIZATION ")')
      DO LN1=1,LNX
        WRITE(6,FMT='(20F12.3)')DOVER(LN1,:)
      ENDDO
!
      DO LN1=1,LNX
        DO LN2=LN1+1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          SVAR=0.5D0*(DT(LN1,LN2)+DT(LN2,LN1))
          DT(LN1,LN2)=SVAR
          DT(LN2,LN1)=SVAR
          SVAR=0.5D0*(DOVER(LN1,LN2)+DOVER(LN2,LN1))
          DOVER(LN1,LN2)=SVAR
          DOVER(LN2,LN1)=SVAR
          SVAR=0.5D0*(DH(LN1,LN2)+DH(LN2,LN1))
          DH(LN1,LN2)=SVAR
          DH(LN2,LN1)=SVAR
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == WRITE INFORMATION TO FILE                                            ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(80("="),T20,"  DTKIN  ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F12.3)')DT(LN1,:)
        ENDDO
        WRITE(6,FMT='(80("="),T20,"  DOVERLAP  ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F12.3)')DOVER(LN1,:)
        ENDDO
        WRITE(6,FMT='(80("="),T20,"  DHAMILTONIAN ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F12.3)')DH(LN1,:)
        ENDDO
!
        WRITE(6,FMT='(80("="),T20,"  AE-OVERLAP ")')
        ALLOCATE(A(LNX,LNX))
        A(:,:)=0.D0
        DO LN1=1,LNX
          DO LN2=1,LNX
            IF(LOX(LN2).NE.LOX(LN1)) CYCLE
            AUX(:)=R(:)**2*AEPHI(:,LN1)*AEPHI(:,LN2)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL)
            A(LN1,LN2)=VAL
          ENDDO
          WRITE(6,FMT='(20F12.3)')A(LN1,:)
        ENDDO
        DEALLOCATE(A)
!
        WRITE(6,FMT='(80("="),T20," <PRO|PSPHI> ")')
        ALLOCATE(A(LNX,LNX))
        A(:,:)=0.D0
        DO LN1=1,LNX
          DO LN2=1,LNX
            IF(LOX(LN2).NE.LOX(LN1)) CYCLE
            AUX(:)=R(:)**2*PRO(:,LN1)*PSPHI(:,LN2)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL)
            A(LN1,LN2)=VAL
          ENDDO
          WRITE(6,FMT='(20F12.3)')A(LN1,:)
        ENDDO
        DEALLOCATE(A)
      END IF
!
!     ==========================================================================
!     == CONSTRUCT PARAMETERS FOR MASS RENORMALIZATION                        ==
!     ==========================================================================
      CALL SETUP_PARMSMASSRENORMALIZATION(GID,NR,ROUT,NB-NC,LOFI(NC+1:) &
     &                           ,FOFI(NC+1:),PSPSIF,PSG2,PSG4)
!
!     ==========================================================================
!     == WRITE DATA TO FILE                                                   ==
!     ==========================================================================
                      CALL TRACE$PASS('WRITING DATA FILES FOR DIAGNOSIS')
      IF(TWRITE) THEN
        WRITE(STRING,FMT='(F3.0)')AEZ
        STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
!       == AE PARTIAL WAVES ====================================================
        CALL SETUP_WRITEPHI(-'AEPHI'//TRIM(STRING),GID,NR,LNX,AEPHI)
!
!       == PS PARTIAL WAVES ====================================================
        CALL SETUP_WRITEPHI(-'PSPHI'//TRIM(STRING),GID,NR,LNX,PSPHI)
!
!       == NODELESS PARTIAL WAVES ==============================================
        CALL SETUP_WRITEPHI(-'NLPHI'//TRIM(STRING),GID,NR,LNX,NLPHI)
!
!       == NODELESS PARTIAL WAVES ==============================================
        CALL SETUP_WRITEPHI(-'QPHI'//TRIM(STRING),GID,NR,LNX,QN)
!
!       == PROJECTOR FUNCTIONS =================================================
        CALL SETUP_WRITEPHI(-'PRO'//TRIM(STRING),GID,NR,LNX,PRO)
!
!       == SCATTERING PSEUDO PARTIAL WAVES =====================================
        CALL SETUP_WRITEPHI(-'PSPHIDOT'//TRIM(STRING),GID,NR,LNX,PSPHIDOT)
!
!       == SCATTERING ALL-ELECTRON PARTIAL WAVES ===============================
        CALL SETUP_WRITEPHI(-'AEPHIDOT'//TRIM(STRING),GID,NR,LNX,AEPHIDOT)
!
!       == ALL-ELECTRON WAVE FUNCTIONS =========================================
        CALL SETUP_WRITEPHI(-'AEPSIF'//TRIM(STRING),GID,NR,NB-NC,AEPSIF)
!
!       == PSEUDO WAVE FUNCTIONS =========================================
        CALL SETUP_WRITEPHI(-'PSPSIF'//TRIM(STRING),GID,NR,NB-NC,PSPSIF)
!
!       == PSEUDO WAVE FUNCTIONS =========================================
        CALL SETUP_WRITEPHI(-'AUGPSIF'//TRIM(STRING),GID,NR,NB-NC,AUGPSIF)
!
!       == POTENTIALS  =========================================================
        ALLOCATE(AUXARR(NR,4))
        AUXARR(:,1)=AEPOT
        AUXARR(:,2)=PSPOT
        AUXARR(:,3)=PSPOT-VADD
        AUXARR(:,4)=VADD
        CALL SETUP_WRITEPHI(-'POT'//TRIM(STRING),GID,NR,4,AUXARR)
        DEALLOCATE(AUXARR)
!
!       == QBAR ================================================================
        DO LN=1,LNX
          WRITE(*,FMT='("LN= ",I3," L=",I1," ISCATT=",I3," QBAR ",F10.5)') &
     &                 LN,LOX(LN),ISCATT(LN)
        ENDDO
      END IF
!
!     ==========================================================================
!     == WRITE SETUP REPORT                                                   ==
!     ==========================================================================
                      CALL TRACE$PASS('WRITING SETUP FILE')
!!$!ERROR UOFI IS NOT INITIALIZED!!!! REMOVED BLOCK
!!$ IT IS NOT USED IN THE SIMULATION CODE
!!$PRINT*,'UOFI= ',UOFI(:,:)
!!$      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
!!$      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
!!$      CALL LINKEDLIST$SELECT(LL_STP,'ATOM',0)
!!$      DO IB=1,NB
!!$        CALL LINKEDLIST$SET(LL_STP,'UPSI',-1,UOFI(:,IB)*PSISCALE(IB))
!!$        CALL LINKEDLIST$SET(LL_STP,'UPSI_SMALL',-1,UOFISM(:,IB)*PSISCALE(IB))
!!$      ENDDO
!
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'AUGMENTATION',0)
      CALL LINKEDLIST$SET(LL_STP,'LNX',0,LNX)
      CALL LINKEDLIST$SET(LL_STP,'LOX',0,LOX)
      DO LN=1,LNX
        CALL LINKEDLIST$SET(LL_STP,'AEPHI',-1,AEPHI(:,LN)*PHISCALE(LN))
        CALL LINKEDLIST$SET(LL_STP,'PSPHI',-1,PSPHI(:,LN)*PHISCALE(LN))
        CALL LINKEDLIST$SET(LL_STP,'NLPHI',-1,NLPHI(:,LN)*PHISCALE(LN))
        CALL LINKEDLIST$SET(LL_STP,'QPHI',-1,QN(:,LN)*PHISCALE(LN))
        CALL LINKEDLIST$SET(LL_STP,'PRO',-1,PRO(:,LN)/PHISCALE(LN))
        CALL LINKEDLIST$SET(LL_STP,'AEPHIDOT',-1,AEPHIDOT(:,LN)*PHISCALE(LN))
        CALL LINKEDLIST$SET(LL_STP,'PSPHIDOT',-1,PSPHIDOT(:,LN)*PHISCALE(LN))
        CALL LINKEDLIST$SET(LL_STP,'NLPHIDOT',-1,PSPHIDOT(:,LN)*PHISCALE(LN))
      ENDDO
      CALL LINKEDLIST$SET(LL_STP,'NV',0,NB-NC)
      DO IB=1,NB-NC  
        CALL LINKEDLIST$SET(LL_STP,'AEPSI',-1,AEPSIF(:,IB))
        CALL LINKEDLIST$SET(LL_STP,'PSPSI',-1,PSPSIF(:,IB))
        CALL LINKEDLIST$SET(LL_STP,'AUGPSI',-1,AUGPSIF(:,IB))
      ENDDO
      CALL LINKEDLIST$SET(LL_STP,'AEPOT',0,AEPOT)
      CALL LINKEDLIST$SET(LL_STP,'PSPOT',0,PSPOT)
      CALL LINKEDLIST$SET(LL_STP,'POTOFPSRHO',0,PSPOT-VADD)
!
!!$DO LN=1,LNX
!!$  CALL RADIAL$DERIVE(GID,NR,QN(:,LN),AUX)
!!$  PSPHI(:,LN)=AUX/QN(:,LN)
!!$  CALL RADIAL$DERIVE(GID,NR,NLPHIDOT(:,LN),AUX)
!!$  PSPHIDOT(:,LN)=AUX/NLPHIDOT(:,LN)
!!$ENDDO
!!$PSPHI(:10,:)=0.D0
!!$PSPHIDOT(:10,:)=0.D0
!!$DO IR=1,NR
!!$  IF(R(IR).GT.6.D0) THEN
!!$    PSPHI(IR,:)=0.D0
!!$    PSPHIDOT(IR,:)=0.D0
!!$  END IF
!!$ENDDO
!!$WRITE(STRING,FMT='(F3.0)')AEZ
!!$STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
!!$CALL LMTO_WRITEPHI(-'DLOG'//TRIM(STRING),GID,NR,LNX,PSPHI)
!!$CALL LMTO_WRITEPHI(-'DLOGDOT'//TRIM(STRING),GID,NR,LNX,PSPHIDOT)

                                CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_PAWDENSITY(GID,NR,LNX,LOX,NB,NC,LOFI,NNOFI,EOFI,FOFI &
     &                            ,AECORE,PSCORE &
     &                            ,AEPHI,AEPHISM,PSPHI,PSPHISM,PRO &
     &                            ,DH,DOVER,AEPOT,PSPOT,VFOCK &
     &                            ,ROUT,TREL,TZORA &
     &                            ,PAWRHO,PSRHO,AEPSIF,PSPSIF,AUGPSIF)
!     **************************************************************************
!     ** PERFORM A PAW CALCULATION FOR A GIVEN SET OF ALL-ELECTRON AND PSEUDO **
!     ** POTENTIALS AND PROVIDE THE ALL-ELECTRON AND PSEUDO DENSITIES USED    **
!     ** FOR UNSCREENING (I.E. THE EXTRACTION OF VADD)                        **
!     **                                                                      **
!     **************************************************************************
      USE RADIALFOCK_MODULE,ONLY: VFOCK_TYPE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NC
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      INTEGER(4),INTENT(IN) :: NNOFI(NB)
      REAL(8)   ,INTENT(IN) :: EOFI(NB)
      REAL(8)   ,INTENT(IN) :: FOFI(NB)
      REAL(8)   ,INTENT(IN) :: ROUT
      REAL(8)   ,INTENT(IN) :: AECORE(NR)
      REAL(8)   ,INTENT(IN) :: PSCORE(NR)
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: AEPHISM(NR,LNX)
      REAL(8)   ,INTENT(IN) :: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: PSPHISM(NR,LNX)
      REAL(8)   ,INTENT(IN) :: PRO(NR,LNX)
      REAL(8)   ,INTENT(IN) :: DH(LNX,LNX)
      REAL(8)   ,INTENT(IN) :: DOVER(LNX,LNX)
      REAL(8)   ,INTENT(IN) :: AEPOT(NR)
      REAL(8)   ,INTENT(IN) :: PSPOT(NR)
      TYPE(VFOCK_TYPE),INTENT(INOUT) :: VFOCK
      LOGICAL(4),INTENT(IN) :: TREL
      LOGICAL(4),INTENT(IN) :: TZORA
      REAL(8)   ,INTENT(OUT):: PAWRHO(NR)        ! 
      REAL(8)   ,INTENT(OUT):: PSRHO(NR)
      REAL(8)   ,INTENT(OUT):: AEPSIF(NR,NB-NC)
      REAL(8)   ,INTENT(OUT):: PSPSIF(NR,NB-NC)
      REAL(8)   ,INTENT(OUT):: AUGPSIF(NR,NB-NC)
      LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
      LOGICAL(4)            :: TVARDREL
      REAL(8)   ,ALLOCATABLE:: AEPHI1(:,:)
      REAL(8)   ,ALLOCATABLE:: AEPHISM1(:,:)
      REAL(8)   ,ALLOCATABLE:: PSPHI1(:,:)
      REAL(8)   ,ALLOCATABLE:: PSPHISM1(:,:)
      REAL(8)   ,ALLOCATABLE:: PRO1(:,:)
      REAL(8)   ,ALLOCATABLE:: DH1(:,:)
      REAL(8)   ,ALLOCATABLE:: DO1(:,:)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)   ,PARAMETER  :: C0LL=Y0
      REAL(8)               :: E
      REAL(8)               :: EOFI1(NB)
      REAL(8)               :: EOFICOMP(2,NB-NC)
      REAL(8)               :: AERHO(NR)
      REAL(8)               :: AUGRHO(NR)
      REAL(8)   ,ALLOCATABLE:: PROJ(:)
      REAL(8)               :: R(NR)
      REAL(8)               :: DREL(NR)
      REAL(8)               :: AUX(NR),AUX1(NR),SVAR,SVAR1,SVAR2,VAL
      REAL(8)               :: G(NR)
      INTEGER(4),ALLOCATABLE:: NPROL(:)
      INTEGER(4)            :: LX
      INTEGER(4)            :: ISO
      INTEGER(4)            :: IVB
      INTEGER(4)            :: L,LN,LN1,LN2,IPRO,IPRO1,IPRO2,IB,IR
      INTEGER(4)            :: NN,NN0,NPRO
      CHARACTER(128)        :: STRING
!     **************************************************************************
                            CALL TRACE$PUSH('SETUP_PAWDENSITY')
      LX=MAX(MAXVAL(LOX),MAXVAL(LOFI))
      CALL RADIAL$R(GID,NR,R)
      EOFI1=EOFI
!
!     == DETERMINE NUMBER OF PROJECTORS FOR EACH ANGULAR MOMENTUM ==============
      ALLOCATE(NPROL(0:LX))
      NPROL(:)=0
      DO LN=1,LNX
        L=LOX(LN)
        NPROL(L)=NPROL(L)+1
      ENDDO
!
!     ==========================================================================
!     ==  ACCUMULATE DENSITY                                                  ==
!     ==     PSRHO IS THE DENSITY OF THE PSEUDO PARTIAL WAVES                 ==
!     ==     PAWRHO IS AT FIRST JUST THE AUGMENTATION DENSITY, BUT LATER, THE ==
!     ==            PSEUDO DENSITY PSRHO IS ADDED                             ==
!     ==     AERHO IS THE ALL-ELECTRON DENSITY OBTAINED WITHOUT PAW           ==
!     ==     AUGRHO IS USED FOR TESTING. IT IS THE DENSITY OF THE ALL-ELECTRON==
!     ==            WACE FUNCTIONS OBTAINED FROM THE PAW CONSTRUCTION         ==
!     ==                                                                      ==
!     ==  PAWRHO AND PSRHO DIFFER BY THE AUGMENTATION. BOTH DO NOT CONSIDER   ==
!     ==  THE FOCK POTENTIAL VFOCK!                                           ==
!     ==========================================================================
      AERHO(:)=AECORE(:)
      AUGRHO(:)=AECORE(:)
      PSRHO(:)=PSCORE(:)
      PAWRHO(:)=AECORE(:)-PSCORE(:)
      EOFICOMP(:,:)=0.D0
      DO L=0,LX
        ISO=0
PRINT*,'=================== L=',L,' ================================='
        NPRO=NPROL(L)
        IF(NPRO.EQ.0) CYCLE
        ALLOCATE(DH1(NPRO,NPRO))
        ALLOCATE(DO1(NPRO,NPRO))
        ALLOCATE(PRO1(NR,NPRO))
        ALLOCATE(AEPHI1(NR,NPRO))
        ALLOCATE(AEPHISM1(NR,NPRO))
        ALLOCATE(PSPHI1(NR,NPRO))
        ALLOCATE(PSPHISM1(NR,NPRO))
        ALLOCATE(PROJ(NPRO))
        IPRO1=0
        DO LN1=1,LNX
          IF(LOX(LN1).NE.L) CYCLE
          IPRO1=IPRO1+1
          PRO1(:,IPRO1)=PRO(:,LN1)
          AEPHI1(:,IPRO1)=AEPHI(:,LN1)
          AEPHISM1(:,IPRO1)=AEPHISM(:,LN1)
          PSPHI1(:,IPRO1)=PSPHI(:,LN1)
          PSPHISM1(:,IPRO1)=PSPHISM(:,LN1)
          IPRO2=0
          DO LN2=1,LNX
            IF(LOX(LN2).NE.L) CYCLE
            IPRO2=IPRO2+1
            DH1(IPRO1,IPRO2)=DH(LN1,LN2)
            DO1(IPRO1,IPRO2)=DOVER(LN1,LN2)
          ENDDO
        ENDDO
!
        DREL=0.D0
        TVARDREL=TREL.AND.(.NOT.TZORA) 
        IF(TREL.AND.TZORA)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,0.D0,DREL)
!
        NN0=-1
        G(:)=0.D0
        DO IB=NC+1,NB
          IVB=IB-NC
          IF(LOFI(IB).NE.L) CYCLE
!PRINT*,'        ---- IB=',IB,' --------------------------------'
          IF(NN0.EQ.-1)NN0=NNOFI(IB)
          E=EOFI1(IB)
!
!         ======================================================================
!         ==  CONSTRUCT ALL-ELECTRON WAVE FUNCTION AEPSIF                     ==
!         ======================================================================
          G(:)=0.D0
          CALL ATOMLIB$BOUNDSTATE(GID,NR,L,0,0.D0,ROUT,TVARDREL &
       &                         ,DREL,G,NNOFI(IB),AEPOT,E,AEPSIF(:,IVB))
          CALL ATOMLIB$UPDATESTATEWITHHF(GID,NR,L,0,DREL,G,AEPOT,VFOCK &
       &                              ,ROUT,E,AEPSIF(:,IB-NC))
!!$          IF(TREL.AND.(.NOT.TZORA)) THEN
!!$            CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,L,ISO &
!!$     &                                    ,DREL,G,AEPSIF(:,IVB),AEPSIFSM(:,IVB))
!!$          ELSE
!!$            AEPSIFSM(:,IB-NC)=0.D0
!!$          END IF
!         == LEAVE OUT THE FOCK POTENTIAL FOR CONSISTENCY ======================
          SVAR1=E
          EOFICOMP(1,IB-NC)=E
!
!         ======================================================================
!         ==  CONSTRUCT PAW PSEUDO WAVE FUNCTION                              ==
!         ======================================================================
!         == THIS DOES NOT WORK WITH THE FOCK POTENTIAL BECAUSE THE NUMBER OF 
!         == NODES DOES NOT INCREASE WITH ENERGY. HENCE THE NODE TRACING FAILS
          NN=NNOFI(IB)-NN0
          G(:)=0.D0
PRINT*,'BEFORE PAWBOUNDSTATE'
PRINT*,'L=',L,' NN=',NN,' E=',E,' NPRO=',NPRO,' ROUT=',ROUT,' IB=',IB
          CALL ATOMLIB$PAWBOUNDSTATE(GID,NR,L,NN,ROUT,PSPOT,NPRO,PRO1,DH1,DO1 &
     &                              ,G,E,PSPSIF(:,IB-NC))
PRINT*,' E=',E
          SVAR2=E
          EOFICOMP(2,IB-NC)=E
!
!         ======================================================================
!         ==  COMPARE PSEUDO AND ALL-ELECTRON RESULT                          ==
!         ======================================================================
          IF(ABS(SVAR2-EOFI1(IB)).GT.1.D-1) THEN
            CALL SETUP_WRITEPHI(-'ERROR_AEPSIF',GID,NR,1,AEPSIF(:,IB-NC))
            CALL SETUP_WRITEPHI(-'ERROR_PSPSIF',GID,NR,1,PSPSIF(:,IB-NC))
            CALL SETUP_WRITEPHI(-'ERROR_PRO',GID,NR,NPRO,PRO1)
            CALL SETUP_WRITEPHI(-'ERROR_PSPOT',GID,NR,1,PSPOT)
            CALL SETUP_WRITEPHI(-'ERROR_G',GID,NR,1,G)
            CALL ERROR$MSG('INACCURACY WHILE UNSCREENING PS POTENTIAL')
            CALL ERROR$MSG('ONE-PARTICLE ENERGIES OBTAINED FROM PAW ')
            CALL ERROR$MSG('DISAGREE WITH THOSE FROM THE AE CALCULATION')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$I4VAL('IB',IB)
            CALL ERROR$I4VAL('NN',NN)
            CALL ERROR$I4VAL('NNOFI(IB)',NNOFI(IB))
            CALL ERROR$I4VAL('NN0',NN0)
            CALL ERROR$R8VAL('ROUT',ROUT)
            CALL ERROR$R8VAL('TARGET: E[EV]',EOFI1(IB)*27.211D0)
            CALL ERROR$R8VAL('TARGET: E[EV]',EOFI(IB)*27.211D0)
            CALL ERROR$R8VAL('AE:     E[EV]',SVAR1*27.211D0)
            CALL ERROR$R8VAL('PAW:    E[EV]',SVAR2*27.211D0)
            CALL ERROR$R8VAL('( PAW-AE)[EV]',(SVAR2-SVAR1)*27.211D0)
            CALL ERROR$R8VAL('(PAW-REF)[EV]',(SVAR2-EOFI1(IB))*27.211D0)
            CALL ERROR$R8VAL('( AE-REF)[EV]',(SVAR1-EOFI1(IB))*27.211D0)
            CALL ERROR$STOP('SETUP_PAWDENSITY')
          END IF
          IF(TTEST) THEN
             WRITE(6,FMT='("DEVIATION OF THE ATOMIC ENERGY LEVELS IN EV")')
             WRITE(6,FMT='("OBTAINED ONCE WITH PAW AND THE AE CALCULATION")')
             WRITE(6,FMT='("DEVIATION PROBABLY DUE TO INCONSISTENCY..")')
             WRITE(6,FMT='("...OF RELATIVISTIC EFFECTS")')
             WRITE(6,FMT='("DEVIATION AE-REF DUE TO DIFFERENCE OF NODAL ...")')
             WRITE(6,FMT='("..AND NODELESS CONSTRUCTION")')
             WRITE(6,FMT='("L",I2," PAW-REF",F10.5,"EV; AE-REF ",F10.5," EV")')&
         &          L,(SVAR2-EOFI1(IB))*27.211D0,(SVAR1-EOFI1(IB))*27.211D0
          END IF
!
!         ==  ENSURE THAT THE TAILS OF AE AND PS WAVE FUNCTION HAVE SAME SIGN ==
!         == AND CUT OFF WAVE FUNCTION TWO GRID-POINTS OUTSIDE ROUT ============
          DO IR=1,NR-2
            IF(R(IR).LT.ROUT) CYCLE
            PSPSIF(IR+2:,IB-NC)=0.D0
            AEPSIF(IR+2:,IB-NC)=0.D0
            IF(PSPSIF(IR-2,IB-NC)*AEPSIF(IR-2,IB-NC).LT.0.D0) &
     &                                       AEPSIF(:,IB-NC)=-AEPSIF(:,IB-NC)
            EXIT
          ENDDO
!
!         == CALCULATE PROJECTIONS PROJ=<P|PS-PSI>  ============================
          DO IPRO=1,NPRO
            AUX(:)=R(:)**2*PSPSIF(:,IB-NC)*PRO1(:,IPRO)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,PROJ(IPRO))
          ENDDO
!
!         == NORMALIZE PS WAVE FUNCTION=========================================
          AUX(:)=R(:)**2*PSPSIF(:,IB-NC)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL)
          DO IPRO1=1,NPRO
            DO IPRO2=1,NPRO
              VAL=VAL+PROJ(IPRO1)*DO1(IPRO1,IPRO2)*PROJ(IPRO2)
            ENDDO
          ENDDO
          PSPSIF(:,IB-NC)=PSPSIF(:,IB-NC)/SQRT(VAL)
          PROJ=PROJ/SQRT(VAL)
!PRINT*,'PROJ ',L,PROJ
!
!         == ADD TO AUGMENTATION DENSITY =======================================
          DO IPRO1=1,NPRO
            DO IPRO2=1,NPRO
              SVAR=FOFI(IB)*PROJ(IPRO1)*PROJ(IPRO2)*C0LL
              AUX= AEPHI1(:,IPRO1)*AEPHI1(:,IPRO2) &
       &         + AEPHISM1(:,IPRO1)*AEPHISM1(:,IPRO2) &
       &         - PSPHI1(:,IPRO1)*PSPHI1(:,IPRO2) &
       &         - PSPHISM1(:,IPRO1)*PSPHISM1(:,IPRO2) 
              PAWRHO(:)=PAWRHO(:)+SVAR*AUX(:)
            ENDDO
          ENDDO
!
!         == AUGMENT PS WAVE FUNCTIONS =========================================
          AUGPSIF(:,IB-NC)=PSPSIF(:,IB-NC)
          DO IPRO=1,NPRO
             AUGPSIF(:,IB-NC)=AUGPSIF(:,IB-NC) &
       &                     +(AEPHI1(:,IPRO)-PSPHI1(:,IPRO))*PROJ(IPRO)
          ENDDO
!
!         == NORMALIZE AE WAVE FUNCTION=========================================
          AUX(:)=R(:)**2*AEPSIF(:,IB-NC)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL)
          AEPSIF(:,IB-NC)=AEPSIF(:,IB-NC)/SQRT(VAL)
!
          PSRHO(:) =PSRHO(:) +FOFI(IB)*PSPSIF(:,IB-NC)**2*Y0
          AUGRHO(:)=AUGRHO(:)+FOFI(IB)*AUGPSIF(:,IB-NC)**2*Y0
          AERHO(:) =AERHO(:) +FOFI(IB)*AEPSIF(:,IB-NC)**2*Y0
        ENDDO !END OF LOOP OVER VALENCE STATES
        DEALLOCATE(DH1)
        DEALLOCATE(DO1)
        DEALLOCATE(PRO1)
        DEALLOCATE(AEPHI1)
        DEALLOCATE(AEPHISM1)
        DEALLOCATE(PSPHI1)
        DEALLOCATE(PSPHISM1)
        DEALLOCATE(PROJ)
      ENDDO      
      PAWRHO(:)=PAWRHO(:)+PSRHO(:)
!
!     == REPORT ===============================================================
      WRITE(6,FMT='(80("="),T20,"  STRAIGHT AE- AND PAW-ENERGY LEVELS ")')
      STRING='("IB=",I2," L=",I2," AE-ENERGY:",F10.5," EV;"'
      STRING=TRIM(ADJUSTL(STRING))//'," PS-ENERGY:",F10.5," EV")' 
      DO IB=NC+1,NB
        WRITE(6,FMT=STRING)IB,LOFI(IB),EOFICOMP(:,IB-NC)*27.211D0
      ENDDO
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_PARMSMASSRENORMALIZATION(GID,NR,RBOX,NB &
     &                                 ,LOFI,FOFI,PSPSI,PSG2,PSG4)
!     **************************************************************************
!     **************************************************************************
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      REAL(8)   ,INTENT(IN) :: FOFI(NB)
      REAL(8)   ,INTENT(IN) :: RBOX
      REAL(8)   ,INTENT(IN) :: PSPSI(NR,NB)
      REAL(8)   ,INTENT(OUT):: PSG2
      REAL(8)   ,INTENT(OUT):: PSG4
      REAL(8)   ,PARAMETER  :: GMAX=15.D0   ! EPW[RY]<GMAX**2 FOR PSI AND RHO
      REAL(8)   ,PARAMETER  :: G1=1.D-3     
      INTEGER(4),PARAMETER  :: NG=250
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: CHARGE  ! PSEUDO VALENCE CHARGE
      REAL(8)               :: EKIN    ! PSEUDO KINETIC ENERGY
      REAL(8)               :: EKING2  ! SUM F<PSPSI|G**4|PSPSI>
      REAL(8)               :: G(NG)
      REAL(8)               :: AUXG(NG)
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)               :: PSPSIG(NG)
      INTEGER(4)            :: GIDG
      REAL(8)               :: DEX
      REAL(8)               :: VAL
      REAL(8)               :: R(NR)
      INTEGER(4)            :: IRBOX ! FIRST GRID POINT BEYOND BOX RADIUS
      REAL(8)               :: CHARGE1,EKIN1
      INTEGER(4)            :: L,IB,IR
!     **************************************************************************
                      CALL TRACE$PUSH('SETUP_PARMSMASSRENORMALIZATION')
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IF(R(IR).LT.RBOX) CYCLE
        IRBOX=IR
        EXIT
      ENDDO
!       
!     ==========================================================================
!     == DEFINE THE RADIAL GRID FOR THE FOURIER TRANSFORM                     ==
!     ==========================================================================
      DEX=LOG(GMAX/G1)/REAL(NG-1,KIND=8)
      CALL RADIAL$NEW('LOG',GIDG)
      CALL RADIAL$SETI4(GIDG,'NR',NG)
      CALL RADIAL$SETR8(GIDG,'R1',G1)
      CALL RADIAL$SETR8(GIDG,'DEX',DEX)
!
!     ==========================================================================
!     == DETERMINE FACTORS FROM BESSEL TRANSFORM                              ==
!     ==========================================================================
      CALL RADIAL$R(GIDG,NG,G)
      CHARGE=0.D0
      EKIN=0.D0
      EKING2=0.D0
      DO IB=1,NB
        L=LOFI(IB)
        AUX(:)=PSPSI(:,IB)
        AUX(IRBOX:)=0.D0
        CALL RADIAL$BESSELTRANSFORM(L,GID,NR,AUX,GIDG,NG,PSPSIG)
        PSPSIG(:)=PSPSIG(:)*SQRT(2.D0/PI)
        AUXG(:)=PSPSIG(:)**2*G(:)**2
        CALL RADIAL$INTEGRAL(GIDG,NG,AUXG,VAL)
        CHARGE=CHARGE+FOFI(IB)*VAL
        AUXG(:)=PSPSIG(:)**2*G(:)**4
        CALL RADIAL$INTEGRAL(GIDG,NG,AUXG,VAL)
        EKIN=EKIN+FOFI(IB)*0.5D0*VAL
        AUXG(:)=PSPSIG(:)**2*G(:)**6
        CALL RADIAL$INTEGRAL(GIDG,NG,AUXG,VAL)
        EKING2=EKING2+FOFI(IB)*0.5D0*VAL
      ENDDO
      PSG2=2.D0*EKIN
      PSG4=2.D0*EKING2
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        EKIN1=0.D0
        CHARGE1=0.D0
        DO IB=1,NB
          L=LOFI(IB)
          AUX(:)=(R(:)*PSPSI(:,IB))**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          CHARGE1=CHARGE1+FOFI(IB)*VAL
!
          CALL RADIAL$DERIVE(GID,NR,R(:)*PSPSI(:,IB),AUX)
          AUX(:)=AUX(:)**2+REAL(L*(L+1),KIND=8)*PSPSI(:,IB)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          EKIN1=EKIN1+FOFI(IB)*0.5D0*VAL
        ENDDO
        WRITE(*,FMT='("PS CHARGE IN G-SPACE ",F10.6,:" IN R-SPACE ",F10.5)') &
     &             CHARGE,CHARGE1
        WRITE(*,FMT='("PS EKIN   IN G-SPACE ",F10.6,:" IN R-SPACE ",F10.5)') &
     &             EKIN,EKIN1
        CALL ERROR$MSG('PLANNED STOP AFTER PRINTOUT OF DIAGNOSTIC INFORMATION')
        CALL ERROR$R8VAL('RBOX',RBOX)
        CALL ERROR$STOP('SETUP_PARMSMASSRENORMALIZATION')
      END IF
                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_TESTGHOST1
!     **************************************************************************
!     **  TEST FOR STATES WITH NEGATIVE NORM.                                 **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: AMAT(:,:)
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: E(:)
      REAL(8)   ,ALLOCATABLE :: U(:,:)
      INTEGER(4)             :: NR
      INTEGER(4)             :: GID
      INTEGER(4)             :: LN,LN1,LN2
      LOGICAL(4),PARAMETER   :: TFORCESTOP=.FALSE.
!     **************************************************************************
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      LNX=THIS%LNX
      ALLOCATE(LOX(LNX))
      LOX=THIS%LOX
      ALLOCATE(AMAT(LNX,LNX))
      ALLOCATE(AUX(NR))
      AMAT(:,:)=0.D0
      DO LN1=1,LNX
        DO LN2=LN1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AUX(:)=R(:)**2*THIS%PRO(:,LN1)*THIS%PRO(:,LN2)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,AMAT(LN1,LN2))
          AMAT(LN2,LN1)=AMAT(LN1,LN2)
        ENDDO
      ENDDO
      AMAT=MATMUL(THIS%DOVER,AMAT)
      DO LN=1,LNX
        AMAT(LN,LN)=1.D0+AMAT(LN,LN)
      ENDDO
!
      ALLOCATE(E(LNX))
      ALLOCATE(U(LNX,LNX))
      CALL LIB$DIAGR8(LNX,AMAT,E,U)
!
      IF(MINVAL(E).LE.0.D0) THEN
        IF(TFORCESTOP) THEN
          CALL ERROR$MSG('1+DO<P|P> IS NOT POSITIVE DEFINITE')
          CALL ERROR$MSG('SETUP WILL PRODUCE GHOST STATES')
          CALL ERROR$CHVAL('SPECIES NAME',THIS%ID)
          DO LN=1,LNX
            IF(E(LN).LE.0.D0) THEN
              DO LN2=1,LNX
                IF(ABS(U(LN2,LN)).GT.1.D-6) THEN
                  CALL ERROR$I4VAL('ANGULAR MOMENTUM ',LOX(LN2))
                  EXIT
                END IF
              ENDDO
              WRITE(*,FMT='("E=",F10.3," U=",20F10.5)')E(LN),U(:,LN)
            END IF
          ENDDO
          CALL ERROR$STOP('SETUP_TESTGHOST1')
        ELSE
          WRITE(*,FMT='(80("="))')
          WRITE(*,FMT='(80("="),T10,A)') &
     &                          ' THIS SETUP IS LIKELY TO PRODUCE GHOST STATES '
          WRITE(*,FMT='(80("="),T10," SPECIES=",A5," Z=",F6.1," ")') &
     &                               TRIM(THIS%ID),THIS%AEZ
          WRITE(*,FMT='(80("="))')
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_TESTGHOST
!     **************************************************************************
!     **  TEST FOR GHOST STATES.                                              **
!     **  METHOD IS INSPIRED BY X. GONZE ET AL, PRB 44, 8503 (1991)           **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: PROPRO(:,:)    !<P|P>
      REAL(8)   ,ALLOCATABLE :: PROHPRO(:,:)   !<P|HLOC|P>
      REAL(8)   ,ALLOCATABLE :: HMAT(:,:)      !
      REAL(8)   ,ALLOCATABLE :: OMAT(:,:)      !
      REAL(8)   ,ALLOCATABLE :: AUX(:),AUX1(:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: E(:)
      REAL(8)   ,ALLOCATABLE :: U(:,:)
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER   :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4)             :: NR
      INTEGER(4)             :: GID
      INTEGER(4)             :: L
      INTEGER(4)             :: LN,LN1,LN2
      LOGICAL(4),PARAMETER   :: TFORCESTOP=.FALSE.
      REAL(8)                :: SVAR
!     **************************************************************************
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
      CALL RADIAL$R(GID,NR,R)
      LNX=THIS%LNX
      ALLOCATE(LOX(LNX))
      LOX=THIS%LOX
!
      CALL SETUP_TESTGHOST_OUTER
      CALL SETUP_TESTGHOST_INNER
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_TESTGHOST_OUTER
!     **************************************************************************
!     **  TEST FOR GHOST STATES.                                              **
!     **  METHOD IS INSPIRED BY X. GONZE ET AL, PRB 44, 8503 (1991)           **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER   :: NBX=10
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER   :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: DREL(:)
      REAL(8)   ,ALLOCATABLE :: G(:)
      REAL(8)   ,ALLOCATABLE :: PHI0(:,:)
      REAL(8)   ,ALLOCATABLE :: PHI(:,:)
      REAL(8)   ,ALLOCATABLE :: SMAT(:,:),SMATIN(:,:)
      REAL(8)                :: E0(NBX)
      REAL(8)                :: E(NBX)
      REAL(8)   ,ALLOCATABLE :: PRO(:,:)
      REAL(8)   ,ALLOCATABLE :: PROPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AUX(:),AUX1(:)
      REAL(8)                :: QMAT(NBX,NBX)
      REAL(8)                :: HMAT(NBX,NBX)
      REAL(8)                :: OMAT(NBX,NBX)
      REAL(8)                :: U(NBX,NBX)
      REAL(8)                :: EPSILONX
      INTEGER(4)             :: NR
      INTEGER(4)             :: GID
      INTEGER(4)             :: L,NN,SO,LN,NPRO,IPRO,I,J
      LOGICAL(4),PARAMETER   :: TFORCESTOP=.FALSE.
      REAL(8)                :: RNS
      REAL(8)                :: SVAR
      REAL(8)                :: RBOX
!     **************************************************************************
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      RBOX=R(NR-3)  ! BOUNDARY CONDITIONS
      LNX=THIS%LNX
      ALLOCATE(LOX(LNX))
      LOX=THIS%LOX
!
!     ==========================================================================
!     ==  LOOP OVER ANGULAR MOMENTA                                           ==
!     ==========================================================================
      ALLOCATE(DREL(NR))
      ALLOCATE(G(NR))
      ALLOCATE(PHI(NR,NBX))
      ALLOCATE(PHI0(NR,NBX))
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
      DO L=0,MAXVAL(LOX)
!
!       ========================================================================
!       ==  COLLECT PROJECTOR FUNCTIONS FOR THIS ANGULAR MOMENTUM             ==
!       ========================================================================
        NPRO=0
        DO LN=1,LNX
          IF(LOX(LN).EQ.L)NPRO=NPRO+1
        ENDDO
        ALLOCATE(PRO(NR,NPRO))
        ALLOCATE(PROPHI(NPRO,NBX))
        IPRO=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L)CYCLE
          IPRO=IPRO+1
          PRO(:,IPRO)=THIS%PRO(:,LN)
        ENDDO
!
!       ========================================================================
!       ==  DETERMINE EIGENSTATES OF THE UNCONSTRAINED PROBLEM                ==
!       ========================================================================
        DO I=1,NBX
          SO=0
          DREL=0.D0
          G=0.D0
          RNS=1.D-3
          E0(I)=0.D0
          IF(I.GT.1)E0(I)=E0(I-1)
          CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RNS,RBOX &
     &                          ,.FALSE.,DREL,G,I-1,THIS%PSPOT,E0(I),PHI0(:,I))
!          WRITE(*,FMT='("L=",I2," IB=",I2," E(LOCAL)=",F10.5)')L,I,E0(I)
!         == NORMALIZE =========================================================
          AUX(:)=R(:)**2*PHI0(:,I)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
          PHI0(:,I)=PHI0(:,I)/SQRT(SVAR)
        ENDDO
!          
!       ========================================================================
!       ==  DETERMINE CONSTRAINT VIOLATION PROPHI = <PRO|PHI0>                ==
!       ========================================================================
        DO I=1,NBX
          DO IPRO=1,NPRO
            AUX(:)=R(:)**2*PRO(:,IPRO)*PHI0(:,I)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,PROPHI(IPRO,I))
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  QMAT = <PHI0|PRO>(<PRO|PRO)^{-1}<PRO|PHI0>                        ==
!       ========================================================================
        ALLOCATE(SMATIN(NPRO,NPRO))
        DO I=1,NPRO
          DO J=I,NPRO
            AUX(:)=R(:)**2*PRO(:,I)*PRO(:,J)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SMATIN(I,J))
            SMATIN(J,I)=SMATIN(I,J)
          ENDDO
        ENDDO
        ALLOCATE(SMAT(NPRO,NPRO))
        CALL LIB$INVERTR8(NPRO,SMATIN,SMAT)
        DEALLOCATE(SMATIN)
        QMAT=MATMUL(TRANSPOSE(PROPHI),MATMUL(SMAT,PROPHI))
        DEALLOCATE(SMAT)
        DEALLOCATE(PROPHI)
!
!       ========================================================================
!       ==  CALCULATE HAMILTONIAN AND OVERLAP OF THE PROJECTED FUNCTIONS      ==
!       ========================================================================
        EPSILONX=MAXVAL(E0)
        HMAT=QMAT*EPSILONX
        OMAT=-QMAT
        DO I=1,NBX
          OMAT(I,I)=OMAT(I,I)+1.D0
          HMAT(I,I)=HMAT(I,I)+E0(I)
          HMAT(:,I)=HMAT(:,I)-QMAT(:,I)*E0(I)
          HMAT(I,:)=HMAT(I,:)-E0(I)*QMAT(I,:)
          DO J=1,NBX
            HMAT(:,J)=HMAT(:,J)+QMAT(:,I)*(E0(I)-EPSILONX)*QMAT(I,J)
          ENDDO
        ENDDO
!       == [ HMAT - OMAT*E(I) ] U(:,I)=0 =======================================
        CALL LIB$GENERALEIGENVALUER8(NBX,HMAT,OMAT,E,U)
        PHI(:,:)=MATMUL(PHI0,U)
!
!       ========================================================================
!       ==  WRITE RESULT                                                      ==
!       ========================================================================
        DO I=1,MIN(NBX,4)
!          CALL ATOMLIB$ORTHOBOUNDSTATE(GID,NR,L,NN,RBOX,THIS%PSPOT,NPRO,PRO,G &
!     &                                ,E2,PHI)
          WRITE(*,FMT='("L=",I2," IB=",I2," E(LOCAL)=",F10.5' &
      &                               //'," E(OUTER)=",2F10.5)')L,I,E0(I),E(I)
        ENDDO
!        WRITE(*,FMT='("SUMTEST:",10F10.5)')(SUM(E(:I)-E0(:I)),I=1,NBX)
!!$CALL SETUP_WRITEPHI('TEST_PHI',GID,NR,NBX,PHI)
!!$CALL SETUP_WRITEPHI('TEST_PHI0',GID,NR,NBX,PHI0)
!!$STOP
        DEALLOCATE(PRO)
      ENDDO ! END OF LOOP OVER ANGULAR MOMENTA
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_TESTGHOST_INNER
!     **************************************************************************
!     **  TEST FOR GHOST STATES.                                              **
!     **  METHOD IS INSPIRED BY X. GONZE ET AL, PRB 44, 8503 (1991)           **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER   :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: PROPRO(:,:)    !<P|P>
      REAL(8)   ,ALLOCATABLE :: PROHPRO(:,:)   !<P|HLOC|P>
      REAL(8)   ,ALLOCATABLE :: HMAT(:,:)      !
      REAL(8)   ,ALLOCATABLE :: OMAT(:,:)      !
      REAL(8)   ,ALLOCATABLE :: AUX(:),AUX1(:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: E(:)
      REAL(8)   ,ALLOCATABLE :: U(:,:)
      INTEGER(4)             :: NR
      INTEGER(4)             :: GID
      INTEGER(4)             :: L
      INTEGER(4)             :: LN,LN1,LN2
      LOGICAL(4),PARAMETER   :: TFORCESTOP=.FALSE.
      REAL(8)                :: SVAR
!     **************************************************************************
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
      CALL RADIAL$R(GID,NR,R)
      LNX=THIS%LNX
      ALLOCATE(LOX(LNX))
      LOX=THIS%LOX
!
!     ==========================================================================
!     == CALCULATE <P|P>  AND <P|HLOC|P>                                      ==
!     ==========================================================================
      ALLOCATE(PROPRO(LNX,LNX))
      ALLOCATE(PROHPRO(LNX,LNX))
      PROPRO(:,:)=0.D0
      PROHPRO(:,:)=0.D0
      DO LN2=1,LNX
        L=LOX(LN2)
!
!       == AUX1 = KINETIC ENERGY OF PRO(LN2) TIMES R**2 ========================
        CALL RADIAL$VERLETD1(GID,NR,THIS%PRO(:,LN2),AUX)
        CALL RADIAL$VERLETD2(GID,NR,THIS%PRO(:,LN2),AUX1)
!       __ INCLUDE ALREADY THE R**2 WEIGHTING FROM THE RADIAL INTEGRATION______ 
        AUX1=AUX1*R**2+2.D0*AUX*R-REAL(L*(L+1),KIND=8)*THIS%PRO(:,LN2)
        AUX1=-0.5D0*AUX1
!
!       == LOOP OVER LN1 =======================================================
        DO LN1=LN2,LNX
          IF(LOX(LN1).NE.L) CYCLE
          AUX(:)=R(:)**2*THIS%PRO(:,LN1)*THIS%PRO(:,LN2)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,PROPRO(LN1,LN2))
          PROPRO(LN2,LN1)=PROPRO(LN1,LN2)
!
!         ==  EXPECTATION VALUE OF LOCAL HAMILTONIAN ===========================
          AUX=AUX*THIS%PSPOT(:)*Y0+THIS%PRO(:,LN1)*AUX1
          CALL RADIAL$INTEGRAL(GID,NR,AUX,PROHPRO(LN1,LN2))
          PROHPRO(LN2,LN1)=PROHPRO(LN1,LN2)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == CALCULATE EIGENVALUES OF LOCAL HAMILTONIAN                           ==
!     ==========================================================================
      ALLOCATE(U(LNX,LNX))
      ALLOCATE(E(LNX))
!
      ALLOCATE(HMAT(LNX,LNX))
      ALLOCATE(OMAT(LNX,LNX))
      HMAT=PROHPRO+MATMUL(PROPRO,MATMUL(THIS%DATH,PROPRO))
      OMAT=PROPRO+MATMUL(PROPRO,MATMUL(THIS%DOVER,PROPRO))
!     == CHECK IF OMAT IS POSITIVE DEFINITE ====================================
      CALL LIB$DIAGR8(LNX,OMAT,E,U)
      IF(MINVAL(E).LE.0.D0) THEN
        CALL ERROR$MSG('OMAT NOT POSITIVE DEFINITE')
        CALL ERROR$R8VAL('SMALLEST EIGENVALUE',MINVAL(E))
        CALL ERROR$STOP('SETUP_TESTGHOST')
      END IF

!     == HMAT*U=OMAT*U*E  WITH   UDAGGER*OMAT*U=1 ==============================
!     == OMAT MUST BE POSITIVE DEFINITE!   =====================================
      CALL LIB$GENERALEIGENVALUER8(LNX,HMAT,OMAT,E,U)
!
!     ==========================================================================
!     == WRITE INFORMATION ON LOCAL HAMILTONIAN                               ==
!     ==========================================================================
      WRITE(*,FMT='(80("="),T10," GHOSTBUSTER FOR=",A5," Z=",F6.1," ")') &
     &                               TRIM(THIS%ID),THIS%AEZ
      DO LN=1,LNX
        L=LOX(LN)
        WRITE(*,FMT='("LN",I2," L=",I2," E(INNER)=",F10.5," H")')LN,L,E(LN)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_TESTSCATTERING(LL_STP)
!     **************************************************************************
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      USE STRINGS_MODULE
      USE LINKEDLIST_MODULE
      USE SETUP_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE),INTENT(INOUT) :: LL_STP
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER   :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4)             :: NFIL
      INTEGER(4)             :: NC
      INTEGER(4)             :: NB
      REAL(8)                :: EMIN
      REAL(8)                :: EMAX
      INTEGER(4),PARAMETER   :: NE=200
      REAL(8)                :: DE
      REAL(8)                :: E
      INTEGER(4)             :: LX
      INTEGER(4)             :: GID ! GRID ID
      INTEGER(4)             :: NR  ! #(GRID POINTS)
      REAL(8)   ,ALLOCATABLE :: R(:)  !(NR) RADIAL GRID
      INTEGER(4)             :: LNX ! #(PARTIAL WAVES)
      INTEGER(4),ALLOCATABLE :: LOX(:) !(LNX) ANGULAR MOMENTA OF PARTIAL WAVES
      REAL(8)   ,ALLOCATABLE :: DREL(:)
      REAL(8)   ,ALLOCATABLE :: G(:)
      REAL(8)   ,ALLOCATABLE :: PHI(:)
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)   ,ALLOCATABLE :: AEPHASE(:,:)
      REAL(8)   ,ALLOCATABLE :: PAWPHASE(:,:)
      REAL(8)   ,ALLOCATABLE :: PRO(:,:)  !(NR,LNX) PROJECTOR FUNCTIONS
      REAL(8)   ,ALLOCATABLE :: DH(:,:)  !(LNX,LNX) HAMILTON-DIFFERENCE
      REAL(8)   ,ALLOCATABLE :: DO(:,:)  !(LNX,LNX) OVERLAP DIFFERENCE
      REAL(8)                :: VAL
      REAL(8)                :: RCOV
      REAL(8)                :: DPHASE
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: IE,LN1,LN2,IPRO1,IPRO2,L,LN,IB
      CHARACTER(64)          :: STRING
      LOGICAL(4),PARAMETER   :: TWRITE=.FALSE.
!     **************************************************************************
      NC=THIS%ATOM%NC
      NB=THIS%ATOM%NB
      EMIN=MIN(MINVAL(THIS%ATOM%EOFI(NC+1:NB))-0.2D0,-1.D0)
      EMAX=MAX(MAXVAL(THIS%ATOM%EOFI(NC+1:NB)),1.D0)
      DE=(EMAX-EMIN)/REAL(NE-1,KIND=8)
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      ALLOCATE(DREL(NR))
      ALLOCATE(G(NR))
      ALLOCATE(PHI(NR))
      ALLOCATE(AUX(NR))
      LNX=THIS%LNX
      ALLOCATE(LOX(LNX))
      LOX=THIS%LOX
      LX=MAX(MAXVAL(THIS%LOX),MAXVAL(THIS%ATOM%LOFI(:NB)))
      CALL PERIODICTABLE$GET(THIS%AEZ,'R(COV)',RCOV)
!
!     ==========================================================================
!     ==  LOOP OVER ANGULAR MOMENTA                                           ==
!     ==========================================================================
      ALLOCATE(AEPHASE(LX+1,NE))
      ALLOCATE(PAWPHASE(LX+1,NE))
      DO L=0,LX
        NPRO=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          NPRO=NPRO+1
        ENDDO
!       __ COUNT OFFSET IN THE NUMBER OF NODES BETWEEN AE AND PS PARTIAL WAVES
        DPHASE=0
        DO IB=1,NC
          IF(THIS%ATOM%LOFI(IB).NE.L) CYCLE
          DPHASE=DPHASE+1.D0
        ENDDO
!       __
        ALLOCATE(DH(NPRO,NPRO))
        ALLOCATE(DO(NPRO,NPRO))
        ALLOCATE(PRO(NR,NPRO))
        IPRO1=0
        DO LN1=1,LNX
          IF(LOX(LN1).NE.L) CYCLE
          IPRO1=IPRO1+1
          PRO(:,IPRO1)=THIS%PRO(:,LN1)
          IPRO2=0
          DO LN2=1,LNX
            IF(LOX(LN2).NE.L) CYCLE
            IPRO2=IPRO2+1
            AUX(:)=THIS%AEPHI(:,LN1)*THIS%AEPHI(:,LN2)*THIS%ATOM%AEPOT*Y0 &
     &            -THIS%PSPHI(:,LN1)*THIS%PSPHI(:,LN2)*THIS%PSPOT*Y0
            AUX(:)=R(:)**2*AUX(:)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
            DH(IPRO1,IPRO2)=THIS%DTKIN(LN1,LN2)+VAL
            DO(IPRO1,IPRO2)=THIS%DOVER(LN1,LN2)
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  LOOP OVER ENERGIES                                                ==
!       ========================================================================
        DO IE=1,NE
          E=EMIN+DE*REAL(IE-1,KIND=8)
          DREL=0.D0
          IF(THIS%SETTING%TREL) THEN
            CALL SCHROEDINGER$DREL(GID,NR,THIS%ATOM%AEPOT,E,DREL)
          END IF
          G(:)=0.D0          
          CALL SCHROEDINGER$SPHERICAL(GID,NR,THIS%ATOM%AEPOT,DREL,0,G,L,E,1,PHI)
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,0.D0,RCOV,AEPHASE(L+1,IE))
          IF(NPRO.GT.0) THEN
            G(:)=0.D0          
            CALL ATOMLIB_PAWDER(GID,NR,L,E,THIS%PSPOT,NPRO,PRO,DH,DO,G,PHI)
          ELSE      
            G(:)=0.D0          
            DREL(:)=0.D0
            CALL SCHROEDINGER$SPHERICAL(GID,NR,THIS%PSPOT,DREL,0,G,L,E,1,PHI)
          END IF
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,0.D0,RCOV,PAWPHASE(L+1,IE))
          PAWPHASE(L+1,IE)=PAWPHASE(L+1,IE)+DPHASE
        ENDDO
        DEALLOCATE(DH)
        DEALLOCATE(DO)
        DEALLOCATE(PRO)
      ENDDO
!
!     ==========================================================================
!     ==  ADD RESULT TO LINKEDLIST                                            ==
!     ==========================================================================
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
      CALL LINKEDLIST$SELECT(LL_STP,'SETUPREPORT',1)
      CALL LINKEDLIST$SELECT(LL_STP,'SCATTERING',0)
      CALL LINKEDLIST$SET(LL_STP,'EMIN',0,EMIN)
      CALL LINKEDLIST$SET(LL_STP,'EMAX',0,EMAX)
      CALL LINKEDLIST$SET(LL_STP,'NE',0,NE)
      CALL LINKEDLIST$SET(LL_STP,'LX',0,LX)
!     == NTH=-1 IN LINKEDLIST SET ADDS A NEW ITEM TO THE LIST =================
      DO L=0,LX
        CALL LINKEDLIST$SET(LL_STP,'AEPHASE',-1,AEPHASE(L+1,:))
        CALL LINKEDLIST$SET(LL_STP,'PAWPHASE',-1,PAWPHASE(L+1,:))
      ENDDO
      CALL LINKEDLIST$SELECT(LL_STP,'~',0)
!
!     ==========================================================================
!     ==  WRITE RESULT TO FILE                                                ==
!     ==========================================================================
      IF(TWRITE) THEN
        WRITE(STRING,FMT='(F3.0)')THIS%AEZ
        STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
        CALL FILEHANDLER$SETFILE('TMP',.FALSE.,-'SCATT'//TRIM(STRING))
        CALL FILEHANDLER$SETSPECIFICATION('TMP','STATUS','UNKNOWN')
        CALL FILEHANDLER$SETSPECIFICATION('TMP','ACTION','WRITE')
        CALL FILEHANDLER$SETSPECIFICATION('TMP','FORM','FORMATTED')
        CALL FILEHANDLER$UNIT('TMP',NFIL)
        REWIND(NFIL)
        DO IE=1,NE
          E=EMIN+DE*REAL(IE-1,KIND=8)
          WRITE(NFIL,FMT='(20F10.5)')E,AEPHASE(:,IE),PAWPHASE(:,IE)
        ENDDO
        CALL FILEHANDLER$CLOSE('TMP')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_MAKEISCATT(AEZ,NB,NC,LOFI,NNOFI,LNX,LOX,ISCATT)
!     **************************************************************************
!     **                                                                      **
!     ** THE PARAMETER ISCATT DETERMINES WHETHER A PARTIAL WAVE IS ATTRIBUTED **
!     ** TO A SEMI-CORE, VALENCE, OR SCATTERING STATE                         **
!     **   ISCATT=-2    STATE BELOW HIGHEST CORE STATE                        **
!     **   ISCATT=-1    SEMI CORE STATE                                       **
!     **   ISCATT= 0    VALENCE STATE (PHI)                                   **
!     **   ISCATT= 1    FIRST SCATTERING STATE (PHIDOT)                       **
!     **   ISCATT= 2    SECOND SCATTERING STATE (PHIDOTDOT)                   **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: AEZ
      INTEGER(4),INTENT(IN)  :: NB
      INTEGER(4),INTENT(IN)  :: NC
      INTEGER(4),INTENT(IN)  :: LOFI(NB)
      INTEGER(4),INTENT(IN)  :: NNOFI(NB)
      INTEGER(4),INTENT(IN)  :: LNX
      INTEGER(4),INTENT(IN)  :: LOX(LNX)
      INTEGER(4),INTENT(OUT) :: ISCATT(LNX)
      INTEGER(4),ALLOCATABLE :: NCL(:)
      INTEGER(4)             :: LX
      INTEGER(4)             :: IB,L,LN
      INTEGER(4)             :: ISVAR
!     **************************************************************************
      LX=MAX(MAXVAL(LOFI),MAXVAL(LOX))
!
!     == DETERMINE HIGHEST CORE STATE FOR EACH ANGULAR MOMENTUM ================
      ALLOCATE(NCL(0:LX))
      NCL(:)=-1
      DO IB=1,NC
        L=LOFI(IB)
        NCL(L)=MAX(NCL(L),IB)
      ENDDO
!
!     == DIVIDE PARTIAL WAVES TO VALENCE AND SCATTERING TYPES ==================
      DO L=0,LX
        IF(L.EQ.0) THEN
          CALL PERIODICTABLE$GET(NINT(AEZ),'#NODES(S)',ISVAR)
        ELSE IF(L.EQ.1) THEN
          CALL PERIODICTABLE$GET(NINT(AEZ),'#NODES(P)',ISVAR)
        ELSE IF(L.EQ.2) THEN
          CALL PERIODICTABLE$GET(NINT(AEZ),'#NODES(D)',ISVAR)
        ELSE IF(L.EQ.3) THEN
          CALL PERIODICTABLE$GET(NINT(AEZ),'#NODES(F)',ISVAR)
        ELSE
          ISVAR=0   !CALL ERROR$STOP('ATOMIC_MAKEISCATT')
        END IF
        IF(NCL(L).GE.1) THEN 
          ISVAR=ISVAR-NNOFI(NCL(L))-1  ! #(NODES) OF VALENCE PS-PARTIAL WAVE
        ELSE
          ISVAR=ISVAR
        END IF
!       == 1+ISVAR IS NOW THE NUMBER OF VALENCE SHELLS
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          ISCATT(LN)=-ISVAR
          ISVAR=ISVAR-1
        ENDDO
      ENDDO
      DEALLOCATE(NCL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_WRITEPHI(FILE,GID,NR,NPHI,PHI)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI(NR,NPHI)
      REAL(8)                 :: PHI1(NPHI)
      INTEGER(4)              :: IR,I
      REAL(8)                 :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        IF(R(IR).GT.3.D0.AND.MAXVAL(ABS(PHI(IR,:))).GT.1.D+3) EXIT
        PHI1(:)=PHI(IR,:)
        DO I=1,NPHI  ! AVOID CONFLICT WITH XMGRACE
          IF(ABS(PHI1(I)).LT.1.D-60) PHI1(I)=0.D0
        ENDDO
        WRITE(100,FMT='(F15.10,2X,100(E25.10,2X))')R(IR),PHI1
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_WRITESCALEDPHI(FILE,GID,NR,NPHI,SCALE,PHI)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: SCALE(NPHI)
      REAL(8)     ,INTENT(IN) :: PHI(NR,NPHI)
      REAL(8)                 :: PHI1(NPHI)
      INTEGER(4)              :: IR,I
      REAL(8)                 :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        IF(R(IR).GT.3.D0.AND.MAXVAL(ABS(PHI(IR,:))).GT.1.D+3) EXIT
        PHI1(:)=PHI(IR,:)*SCALE(:)
        DO I=1,NPHI  ! AVOID CONFLICT WITH XMGRACE
          IF(ABS(PHI1(I)).LT.1.D-60) PHI1(I)=0.D0
        ENDDO
        WRITE(100,FMT='(F15.10,2X,20(E25.10,2X))')R(IR),PHI1
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMIC_UNSCREEN(GID,NR,RBOX,AEZ,AERHO,PSRHO,PSPOT,RCSM,VADD)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     ** CAUTION: I BELIEVE THAT ROUT=R(NR-3) SHOULD BE CHOSEN IN PLACE OF RBOX*
!     **                                                                      **
!     **************************************************************************
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      REAL(8)     ,INTENT(IN) :: RBOX
      REAL(8)     ,INTENT(IN) :: AEZ
      REAL(8)     ,INTENT(IN) :: AERHO(NR)
      REAL(8)     ,INTENT(IN) :: PSRHO(NR)
      REAL(8)     ,INTENT(IN) :: PSPOT(NR)
      REAL(8)     ,INTENT(IN) :: RCSM
      REAL(8)     ,INTENT(OUT):: VADD(NR)
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                 :: R(NR)
      REAL(8)                 :: AUX(NR),AUX1(NR),SVAR
      REAL(8)                 :: POT(NR),POTXC(NR),POTH(NR)
      REAL(8)                 :: QLM
      REAL(8)                 :: ALPHA,CL
      REAL(8)                 :: GRHO(NR)
      INTEGER(4)              :: IR
      REAL(8)                 :: RH,GRHO2,VXC,VGXC,EXC,DUMMY1,DUMMY2,DUMMY3
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == MOMENT OF DIFFERENCE CHARGE DENSITY                                  ==
!     ==========================================================================
      AUX(:)=(AERHO(:)-PSRHO(:))*R(:)**2
!     == THE INTEGRATION MUST BE PERFORMED OUTWARD TO THE END, BECAUSE =========
!     == SMALL DEVIATIONS LEAD TO LONG-RANGE TAIL IN VADD ======================
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,QLM)
      QLM=QLM-AEZ*Y0    ! CHARGE =QM/Y0
!
!     ==========================================================================
!     == ADD COMPENSATION DENSITY AND DETERMINE ELECTROSTATIC POTENTIAL       ==
!     ==========================================================================
      ALPHA=1.D0/RCSM**2
      CALL GAUSSN(0,ALPHA,CL)
      SVAR=QLM*CL
      AUX(:)=PSRHO(:)+SVAR*EXP(-ALPHA*R(:)**2)
!== CHECK CHARGE NEUTRALITY PSRHO+AUGMENTATION =================================
CALL RADIAL$INTEGRATE(GID,NR,AUX*R(:)**2,AUX1)
CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
PRINT*,'PSEUDO+AUGMENTATION CHARGE ',SVAR*Y0*4.D0*PI,' (SHOULD BE ZERO)'
      CALL RADIAL$POISSON(GID,NR,0,AUX,POTH)
!     == SET POTENTIAL AT THE BOX RADIUS TO ZERO ==============================
      CALL RADIAL$VALUE(GID,NR,POTH,RBOX,SVAR)
      POTH=POTH-SVAR
!
!     ==========================================================================
!     == EXCHANGE AND CORRELATION                                             ==
!     ==========================================================================
      CALL RADIAL$DERIVE(GID,NR,PSRHO(:),GRHO)
      GRHO(1)=0.D0 ! R(1)=0 DENSITY GRADIENT VANISHES AT THE ORIGIN
      DO IR=1,NR
        RH=PSRHO(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        POTXC(IR)=VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO

      AUX(:)=R(:)**2*GRHO(:)
      CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
      GRHO(2:)=AUX1(2:)/R(2:)**2
      GRHO(1:5)=GRHO(5) ! AVOID ERRORS DUE TO TERMINATION OF THE GRID
                        ! 5 POINTS OFFSET FOR 5-POINT FORMULA APPLIED TWICE...
      POTXC(:)=POTXC(:)-GRHO(:)
!
!     ==========================================================================
!     ==  CUT OF POTENTIAL FOR LOW DENSITIES                                  ==
!     == ATTENTION: THIS CUTOFF MUST BE COMPLETELY CONSISTENT WITH            ==
!     ==            ATOMLIB$BOXVOFRHO                                         ==
!     ==========================================================================
      DO IR=1,NR
        IF(AERHO(IR)*Y0.LT.1.D-6) POTXC(IR)=0.D0
      ENDDO
!
!     ==========================================================================
!     == SET POTXC TO ZERO BEYOND RBOX  THE RADIUS WHERE THE WAVE FUNCTIONS   ==
!     == HAVE THEIR ZERO
!     ==========================================================================
      DO IR=1,NR
        IF(R(IR).GE.RBOX) THEN
          POTXC(IR:)=0.D0
          POTH(IR:)=0.D0
          EXIT
        END IF
      ENDDO
!
!     ==========================================================================
!     == VADD                                                                 ==
!     ==========================================================================
      POT=POTH+POTXC
      VADD(:)=PSPOT(:)-POT(:)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_BIORTHOMATRICES(GID,NR,RBOX,LNX,LOX,PSPHI,PRO &
     &                          ,TRANSPHI,TRANSPRO)
!     **************************************************************************
!     ** DETERMINES THE MATRICES TRANSPHI AND TRANSPRO SUCH THAT              **
!     **     |PHI-BAR>:=|PHI>TRANSSPHI                                        **
!     **     |PRO-BAR>:=|PRO>TRANSSPRO                                        **
!     **  OBEY  <PHIBAR(I)|PROBAR(J)>=DELTA(I,J)    (KRONECKER DELTA)         **
!     **                                                                      **
!     **************************************************************************
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      REAL(8)    ,INTENT(IN)     :: RBOX
      INTEGER(4) ,INTENT(IN)     :: LNX
      INTEGER(4) ,INTENT(IN)     :: LOX(LNX)
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,LNX)
      REAL(8)    ,INTENT(IN)     :: PRO(NR,LNX)
      REAL(8)    ,INTENT(OUT)    :: TRANSPHI(LNX,LNX)
      REAL(8)    ,INTENT(OUT)    :: TRANSPRO(LNX,LNX)
      LOGICAL(4),PARAMETER       :: TTEST=.TRUE.
      INTEGER(4)                 :: LN,LN1,LN2
      REAL(8)                    :: AUX(NR),AUX1(NR)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: SVAR
      REAL(8)                    :: PROPSI(LNX,LNX)
      REAL(8)                    :: TRANSPROPSI(LNX,LNX)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == CALCULATE INITIAL VIOLATION <PSPHI|BARE-PRO> OF BIORTHOGONALITY      ==
!     ==========================================================================
      PROPSI(:,:)=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AUX(:)=PRO(:,LN1)*PSPHI(:,LN2)*R(:)**2
!         == PROJECTOR FUNCTIONS HAVE A STRICT CUTOFF. INTEGRATION IS DONE =====
!         == TO THE END  =======================================================
!          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)       
!          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)       
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PROPSI(LN1,LN2)=SVAR
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == COLLECT TRANSFORMATION MATRIX BETWEEN NEW AND OLD                    ==
!     ==========================================================================
      TRANSPRO(:,:)=0.D0
      TRANSPHI(:,:)=0.D0
      DO LN=1,LNX
        TRANSPRO(LN,LN)=1.D0
        TRANSPHI(LN,LN)=1.D0
      ENDDO
      TRANSPROPSI(:,:)=PROPSI(:,:)
!
      DO LN1=1,LNX
        DO LN2=1,LN1-1
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
!         == ORTHOGONALIZE PARTIAL WAVES =======================================
          SVAR=TRANSPROPSI(LN2,LN1)
          TRANSPHI(:,LN1)   =TRANSPHI(:,LN1)-TRANSPHI(:,LN2)*SVAR
          TRANSPROPSI(:,LN1)=TRANSPROPSI(:,LN1)-TRANSPROPSI(:,LN2)*SVAR     
!         == ORTHOGONALIZE PROJECTOR===========================================
          SVAR=TRANSPROPSI(LN1,LN2)
          TRANSPRO(:,LN1)   =TRANSPRO(:,LN1)-TRANSPRO(:,LN2)*SVAR
          TRANSPROPSI(LN1,:)=TRANSPROPSI(LN1,:)-TRANSPROPSI(LN2,:)*SVAR 
        ENDDO
        SVAR=TRANSPROPSI(LN1,LN1)
        TRANSPRO(:,LN1)=TRANSPRO(:,LN1)/SVAR
        TRANSPROPSI(LN1,:)=TRANSPROPSI(LN1,:)/SVAR
      ENDDO
!
!     ==========================================================================
!     == CHECK RESULT                                                         ==
!     ==========================================================================
      IF(TTEST) THEN
        TRANSPROPSI(:,:)=MATMUL(TRANSPOSE(TRANSPRO),MATMUL(PROPSI,TRANSPHI))
        DO LN=1,LNX
          TRANSPROPSI(LN,LN)=TRANSPROPSI(LN,LN)-1.D0
        ENDDO
        SVAR=MAXVAL(TRANSPROPSI)
        IF(SVAR.GT.1.D-5) THEN
          DO LN=1,LNX
            WRITE(*,FMT='(20E10.2)')TRANSPROPSI(:,LN)
          ENDDO
          CALL ERROR$MSG('BIORTHOGONALIZATION INACCURATE')
          CALL ERROR$R8VAL('MAX. DEV.',SVAR)
          CALL ERROR$STOP('SETUP_BIORTHOMATRICES')
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMIC_MAKEPSPHI_HBS(GID,NR,LNX,LOX,EOFLN,RC,LAMBDA,PSPOT &
     &                               ,RBND,PSPHI,TPSPHI)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)     :: GID
      INTEGER(4),INTENT(IN)     :: NR
      INTEGER(4),INTENT(IN)     :: LNX
      INTEGER(4),INTENT(IN)     :: LOX(LNX)
      REAL(8)   ,INTENT(IN)     :: EOFLN(LNX)
      REAL(8)   ,INTENT(IN)     :: LAMBDA(LNX)
      REAL(8)   ,INTENT(IN)     :: RC(LNX)
      REAL(8)   ,INTENT(IN)     :: PSPOT(NR)
      REAL(8)   ,INTENT(IN)     :: RBND
      REAL(8)   ,INTENT(INOUT)  :: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(INOUT)  :: TPSPHI(NR,LNX)
      INTEGER(4),PARAMETER      :: NITER=100
      INTEGER(4),PARAMETER      :: ISO=0
      REAL(8)   ,PARAMETER      :: TOL=1.D-11
      REAL(8)   ,PARAMETER      :: CMIN=1.D-8
      REAL(8)   ,PARAMETER      :: RBNDX=6.D0
      REAL(8)   ,PARAMETER      :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER      :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                   :: E
      REAL(8)                   :: DREL(NR),G(NR),POT(NR),PHI(NR)
      REAL(8)                   :: C(NR)
      REAL(8)                   :: R(NR)
      REAL(8)                   :: X0,XM,Z0,ZM,DX
      REAL(8)                   :: PHIPHASE
      INTEGER(4)                :: IRBND,IRBND2
      INTEGER(4)                :: ISTART,IBI
      INTEGER(4)                :: L
      INTEGER(4)                :: LN,ITER,IR,II(1)
      LOGICAL(4)                :: CONVG
      REAL(8)                   :: SVAR,SVAR1
      REAL(8)                   :: ARR1(5),ARR2(5)
      REAL(8)                   :: RBND2
!     **************************************************************************
                                CALL TRACE$PUSH('ATOMIC_MAKEPSPHI_HBS')
      CALL RADIAL$R(GID,NR,R)
      DO LN=1,LNX
        L=LOX(LN)
        E=EOFLN(LN)
!
!       ========================================================================
!       ==  SET UP THE ADDITIONAL POTENTIAL USED TO ADJUST THE WAVE FUNCTION  ==
!       ==  TO BECOME EQUAL WITH THE ALL-ELECTRON WAVE FUNCTION               ==
!       ========================================================================
        C(:)=EXP(-(R(:)/RC(LN))**LAMBDA(LN))
!       == CUT OFF C, IF IT FALLS BELOW MINIMUM ================================
        SVAR=RC(LN)*(-LOG(CMIN))**(1.D0/LAMBDA(LN))
        IRBND=0
        DO IR=1,NR
          IRBND=IR
          IF(R(IR).GT.SVAR) EXIT
        ENDDO
        C(:)=(C(:)-C(IRBND))/(1.D0-C(IRBND))
        C(IRBND:)=0.D0
!
!       ========================================================================
!       == DETERMINE RBND TO MATCH THE PSEUDO PARTIAL WAVE. REQUIRING THE     ==
!       == THE POTENTIALS BEYOND RBND ARE EQUAL                               ==
!       ==   1. START WITH POINT WHERE IZ ZERO                                ==
!       ==   2. EXPAND IF AE AND PS POTENTIALS DEVIATE ATE RBND               ==
!       ==   3. OVERWRITE BY RBNDX IF MATCHING RADIUS IS TOO LARGE            ==
!       ========================================================================
!!$!       == FIND OUTERMOST POINT CONSIDERING THE POTENTIAL ===================
!!$        SVAR1=MAXVAL(ABS(PSPHI(IRBND:,LN)))
!!$        DO IR=NR-2,1,-1
!!$          IF(IR.LT.IRBND) EXIT
!!$          SVAR=TPSPHI(IR,LN)+(PSPOT(IR)*Y0-E)*PSPHI(IR,LN)
!!$          IF(ABS(SVAR).GT.1.D-4*SVAR1) THEN
!!$            IRBND=IR
!!$            EXIT
!!$          END IF
!!$        ENDDO
        RBND2=MAX(RBND,R(IRBND))
        DO IR=1,NR
          IF(R(IR).GE.RBND2) THEN
            IRBND2=IR
            EXIT
          END IF
        ENDDO
        IRBND=IRBND2   !! CRITICAL CHANGE!!: THE MATCHING RADIUS IS RBND2 BUT 
                       !! SCALING AND TAIL REPLACEMENT START AT IRBND
!
!       == THE CUTOFF CAN BE LARGE, IN PARTICULAR WHEN THE PARTIAL WAVES ARE  ==
!       == MATCHED TO THE NODELESS PARTIAL WAVES. IF THIS IS UNACCEPTABLE, =====
!       == ONE SHOULD MATCH TO THE ALL-ELECTRON PARTIAL WAVES OF INCLUDE THE  ==
!       == SEMI-CORE ORBITALS AS VALENCE ELECTRONS
        IF(RBND2.GT.RBNDX) THEN
          CALL ERROR$MSG('RADIUS FOR LOGARITHMIC DERIVATIVE EXCEEDS LIMIT')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('LN',LN)
          CALL ERROR$R8VAL('E',E)
          CALL ERROR$R8VAL('RBND2',RBND2)
          CALL ERROR$R8VAL('INTERNAL LIMIT',RBNDX)
          CALL ERROR$I4VAL('IRBND',IRBND)
          CALL ERROR$I4VAL('NR',NR)
          CALL SETUP_WRITEPHI('ERRORDATA1.DAT',GID,NR,1 &
      &                      ,TPSPHI(:,LN)+(PSPOT(:)*Y0-E)*PSPHI(:,LN))
          CALL SETUP_WRITEPHI('ERRORDATA2.DAT',GID,NR,1,(PSPOT(:)*Y0-E))
          CALL SETUP_WRITEPHI('ERRORDATA3.DAT',GID,NR,1,TPSPHI(:,LN))
          CALL SETUP_WRITEPHI('ERRORDATA4.DAT',GID,NR,1,PSPHI(:,LN))
          CALL ERROR$STOP('ATOMIC_MAKEPSPHI_HBS')
        END IF
!!$!
!!$!       == LIMIT RBND TO A REASONABLE MAXIMUM ===============================
!!$        IF(R(IRBND).GT.RBNDX) THEN
!!$          WRITE(*,FMT='("OVERWRITING MATCHING RADIUS ESTIMATED AS ",F10.5)')&
!!$     &                                                               R(IRBND)
!!$          DO IR=1,IRBND
!!$            IF(R(IR).LT.RBNDX) CYCLE
!!$            IRBND=IR
!!$            EXIT
!!$          ENDDO
!!$          WRITE(*,FMT='("NEW RADIUS ",F10.5)')R(IRBND)
!!$        END IF
!!$!        RBND2=RBND !OLD CHOICE
!
!       ========================================================================
!       ==  CORRECT FOR NODES LYING WITHIN 0.3. THEY CAN OCCUR                ==
!       ==  WITH A (NONLOCAL) FOCK POTENTIAL AND UPSET THE FORMALISM          ==
!       ========================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PSPHI(:,LN),0.D0,RBND2,PHIPHASE)
PRINT*,'L ',L,' E=',E,'PHIPHASE BEFORE ADJUSTMENT= ',PHIPHASE
        DO IR=1,IRBND
          IF(PSPHI(IR,LN)*PSPHI(IR+1,LN).LT.0.D0) THEN
            PHIPHASE=PHIPHASE-1.D0
            WRITE(*,FMT= &
     &       '("L=",I2,"LN",I3,"NR. OF NODES REDUCED BY ONE RELATIVE TO QN. R=",F10.5)') &
     &       L,LN,R(IR)
          END IF           
          IF(R(IR).GT.0.3D0) EXIT  !0.3D0 IS REQUIRED FOR URANIUM
        ENDDO
PRINT*,'L ',L,' E=',E,'PHIPHASE AFTER ADJUSTMENT= ',PHIPHASE
!
!       ========================================================================
!       ==  LOOP TO FIND PSEUDO PARTIAL WAVE                                  ==
!       ========================================================================
        ISTART=1
        X0=0.D0
        DX=1.D-2
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
        DO ITER=1,NITER
          POT(:)=PSPOT(:)+X0*C(:)
          DREL(:)=0.D0
          G(:)=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,ISO,G,L,E,1,PHI)
          IF(.NOT.(PHI(IRBND+2).GT.0.D0.OR.PHI(IRBND+2).LE.0.D0)) THEN
            CALL ERROR$MSG('WAVE FUNCTION IS NOT A NUMBER')
            CALL ERROR$MSG('OVERFLOW OF WAVE FUNCTION ENCOUNTERED')
            CALL ERROR$STOP('ATOMIC_MAKEPSPHI_HBS')
          END IF
          CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,0.D0,RBND2,Z0)
          Z0=PHIPHASE-Z0
          CONVG=(ABS(2.D0*DX).LE.TOL)
          IF(CONVG) EXIT
          IF(ABS(X0).GT.1.D+5) THEN
            CALL ERROR$MSG('NO SUITABLE PSEUDO-PARTIAL WAVE FOUND')
            CALL ERROR$MSG('EITHER INCREASE CUTOFF-RADIUS OR ..')
            CALL ERROR$MSG('.. DISCARD THIS PARTIAL WAVE IF IT HAS HIGH ENERGY')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$R8VAL('E',E)
            CALL ERROR$R8VAL('X0',X0)
            CALL ERROR$R8VAL('RBND2',RBND2)
            CALL ERROR$R8VAL('R(IRBND)',R(IRBND))
            PHI(IRBND+2:)=0.D0  ! MASK OVERFLOWS
            CALL SETUP_WRITEPHI('ERRORDATA1.DAT',GID,NR,1,PHI)
            CALL SETUP_WRITEPHI('ERRORDATA2.DAT',GID,NR,1,PSPHI(:,LN))
            CALL ERROR$STOP('ATOMIC_MAKEPSPHI_HBS')
          END IF
          CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
        ENDDO
        IF(.NOT.CONVG) THEN
          CALL ERROR$MSG('LOOP NOT CONVERGED')
          CALL ERROR$STOP('ATOMIC_MAKEPSPHI_HBS')
        END IF
!
!       ========================================================================
!       ==  RESCALE PSEUDO PARTIAL WAVES SO THAT THEY MATCH TO AE PARTIAL WAVES=
!       ========================================================================
!       == DO NOT RESCALE AT THE NODAL PLANE, BUT 5 POINTS INWARD....        
        ARR1(:)=PSPHI(IRBND-4:IRBND,LN)
        ARR2(:)=PHI(IRBND-4:IRBND)
        II=MAXLOC(ABS(ARR1*ARR2))
        SVAR=ARR1(II(1))/ARR2(II(1))
        PHI(:)=PHI(:)*SVAR
!
!       == OVERWRITE NODELESS INPUT PARTIAL WAVE BY PSEUDO PARTIAL WAVE
!       == THE PARTIAL WAVE OUTSIDE RBND IS SET TO THE INPUT FUNCTION
!       == TO AVOID AN EXPONENTIALLY GROWING DEVIATION. TPSPHI IS NOT 
!       == CONSTRUCTED FROM THE INPUT TPSPHI TO AVOID PICKING UP ??
!!$        PSPHI(:IRBND,LN) =PHI(:IRBND)
!!$        TPSPHI(:IRBND,LN)=(E-POT(:IRBND)*Y0)*PHI(:IRBND)
        PSPHI(:IRBND2,LN) =PHI(:IRBND2)
        TPSPHI(:,LN)=(E-POT(:)*Y0)*PSPHI(:,LN)
      ENDDO
                                CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMIC_MAKEPSPHI(GID,NR,RC,L,NB,PSPHI,TPSPHI)
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
!     ******************************PETER BLOECHL, GOSLAR 2010******************
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
      REAL(8)   ,PARAMETER      :: PI=4.D0*ATAN(1.D0)
      REAL(8)                   :: R(NR)
      INTEGER(4)                :: IR,I,J,ITER
      REAL(8)                   :: TOL=1.D-12
      REAL(8)                   :: VAL
      REAL(8)                   :: X(NX)
      REAL(8)                   :: JL(NX),DJLDR(NX),PHASE(NX)
      REAL(8)                   :: JLOFKR(NR,NZEROS),TJLOFKR(NR,NZEROS)
      REAL(8)                   :: AUX1(NR),AUX2(NR)
      REAL(8)                   :: NN
      INTEGER(4)                :: ISTART,IBI
      REAL(8)                   :: X0,Y0,DX,XM,YM
      REAL(8)                   :: KI(NZEROS)
      REAL(8)                   :: SVAR
      INTEGER(4)                :: IB
      REAL(8)                   :: SUPPHI(NR,NB)
      REAL(8)                   :: SUPTPHI(NR,NB)
      REAL(8)                   :: SHIFT
      REAL(8)                   :: X1,XDEX
      REAL(8)                   :: OFFSETNN
      INTEGER(4)                :: GIDX
      INTEGER(4)                :: IRMATCH
!     **************************************************************************
!CALL SETUP_WRITEPHI('TESTCONSTRUCTPSPHI1',GID,NR,NB,PSPHI)
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
      CALL SCHROEDINGER$PHASESHIFT(GID,NR,SUPPHI(:,1),0.D0,RC,SHIFT)
      OFFSETNN=REAL(INT(SHIFT),KIND=8)
!     
!     ==========================================================================
!     == PSEUDIZE EVERY FUNCTION ON THE SUPPORT ENERGY GRID                   ==
!     ==========================================================================
      DO IB=1,NB
!       ========================================================================
!       == DETERMINE K-VALUES FOR WHICH THE LOGARITHIC DERIVATIVE MATCHES     ==
!       ========================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,SUPPHI(:,IB),0.D0,RC,SHIFT)
        SHIFT=SHIFT-OFFSETNN
        IF(SHIFT.LT.PHASE(1)) THEN
          CALL ERROR$MSG('MATCHING INCLUDING FIRST BESSEL FUNCTION NOT POSSIBLE')
          CALL ERROR$MSG('TRY TO INCREASE MATCHING RADIUS')
          CALL ERROR$I4VAL('IB',IB)
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$R8VAL('SHIFT',SHIFT)
          CALL ERROR$R8VAL('PHASE(1)',PHASE(1))
          CALL ERROR$R8VAL('OFFSETNN',OFFSETNN)
          CALL ERROR$R8VAL('RC',RC)
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
      SUBROUTINE ATOMIC_PSEUDIZE(GID,NR,POW,TVAL,VAL0_,RC,AEF,PSF)
!     **************************************************************************
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      REAL(8)    ,INTENT(IN)     :: POW    !POWER
      LOGICAL(4) ,INTENT(IN)     :: TVAL   ! VALUE AT THE ORIGIN PRESCRIBED
      REAL(8)    ,INTENT(IN)     :: VAL0_  ! VALUE OF PSEUDIZED FUNCTION AT R=0
      REAL(8)    ,INTENT(IN)     :: RC     ! MATCHING RADIUS
      REAL(8)    ,INTENT(IN)     :: AEF(NR)
      REAL(8)    ,INTENT(OUT)    :: PSF(NR)
      REAL(8)                    :: A,B,C
      REAL(8)                    :: VAL,DER
      INTEGER(4)                 :: IR,IR0
      REAL(8)                    :: R(NR)
      LOGICAL(4),PARAMETER       :: TTEST=.FALSE.
      LOGICAL                    :: TOK
      REAL(8)                    :: TOL=1.D-10
      REAL(8)                    :: SVAR1,SVAR2,SVAR3
!     **************************************************************************
      CALL RADIAL$VALUE(GID,NR,AEF,RC,VAL)
      CALL RADIAL$DERIVATIVE(GID,NR,AEF,RC,DER)
      IF(TVAL) THEN
        A=VAL0_
        C=0.5D0*(RC*DER-POW*(VAL-A))
        B=VAL-A-C
        B=B/RC**POW
        C=C/RC**(POW+2.D0)
      ELSE
        B=RC*DER/POW
        A=VAL-B
        B=B/RC**POW
        C=0.D0
      END IF
!     == CALCULATE PS FUNCTION
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IF(R(IR).GT.RC) THEN
          IR0=IR
          EXIT
        END IF
      ENDDO
      PSF(:IR0-1)=A+R(:IR0-1)**POW*(B+C*R(:IR0-1)**2)
      PSF(IR0:)=AEF(IR0:)
!
!     ==========================================================================
!     ==  OPTIONAL TEST                                                       ==
!     ==========================================================================
      IF(TTEST) THEN
        SVAR1=A+B*RC**POW+C*RC**(POW+2.D0)
        SVAR2=POW*B*RC**(POW-1)+(POW+2.D0)*C*RC**(POW+1.D0)
        SVAR3=A
        IF(TVAL) THEN
          TOK=.TRUE.
          TOK=TOK.AND.ABS(SVAR1-VAL).LE.TOL
          TOK=TOK.AND.ABS(SVAR2-DER).LE.TOL
          IF(TVAL)TOK=TOK.AND.ABS(SVAR3-VAL0_).LE.TOL
          IF(.NOT.TOK) THEN
            CALL ERROR$MSG('INTERNAL ERROR')
            CALL ERROR$R8VAL('F(RC) TARGET    ',VAL)
            CALL ERROR$R8VAL('F(RC) ACTUAL   ',SVAR1)
            CALL ERROR$R8VAL('DF/DR|_RC TARGET',DER)
            CALL ERROR$R8VAL('DF/DR|_RC ACTUAL',SVAR2)
            IF(TVAL) THEN
              CALL ERROR$R8VAL('F(0) TARGET',VAL0_)
              CALL ERROR$R8VAL('F(0) ACTUAL',SVAR3)
            END IF
            CALL ERROR$STOP('ATOMIC$PSEUDIZE')
          END IF
        END IF
      END IF
      RETURN
      END
!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE SETUP_BUILDPARMS()
!!$!     **************************************************************************
!!$!     ** GENERATES AUTOMATICALLY A PARAMETER FILE FOR THE SETUP CONSTRUCTION  **
!!$!     ** TYPE: NDLSS                                                          **
!!$!     ******************************PETER BLOECHL, GOSLAR 2014******************
!!$      USE PERIODICTABLE_MODULE
!!$      USE CONSTANTS_MODULE
!!$      IMPLICIT NONE
!!$      INTEGER(4),PARAMETER :: NFIL=1001
!!$      INTEGER(4),PARAMETER :: NFIL2=1002
!!$      INTEGER(4),PARAMETER :: LX=4
!!$      INTEGER(4),PARAMETER :: NSTPTYPES=5
!!$      CHARACTER(20)        :: STPTYPE(NSTPTYPES)
!!$      INTEGER(4)           :: IZ,ISTPTYPE
!!$      CHARACTER(2)         :: EL
!!$      CHARACTER(128)       :: ID
!!$      REAL(8)              :: AEZ
!!$      REAL(8)              :: ZV
!!$      REAL(8)              :: RBOX
!!$      REAL(8)              :: RCSM
!!$      REAL(8)              :: RCL(LX+1)
!!$      REAL(8)              :: LAMBDA(LX+1)
!!$      REAL(8)              :: POTPOW
!!$      REAL(8)              :: POTRC
!!$      LOGICAL(4)           :: TPOTVAL
!!$      REAL(8)              :: POTVAL
!!$      REAL(8)              :: CORPOW
!!$      REAL(8)              :: CORRC
!!$      LOGICAL(4)           :: TCORVAL
!!$      REAL(8)              :: CORVAL
!!$      REAL(8)              :: DMIN
!!$      REAL(8)              :: DMAX
!!$      REAL(8)              :: RMAX
!!$      CHARACTER(12)        :: TYPE
!!$!     **************************************************************************
!!$!
!!$!     ==========================================================================
!!$!     == DEFINE STRINGS FOR DIFFERENT SETUP CONSTRUCTION TYPES                ==
!!$!     ==========================================================================
!!$      STPTYPE(1)='NDLSS_V0'
!!$      STPTYPE(2)='NDLSS_SC_V0'
!!$      STPTYPE(3)='.75_6.0'
!!$      STPTYPE(4)='HBS'
!!$      STPTYPE(5)='HBS_SC'
!!$!
!!$!     ==========================================================================
!!$!     == LOOP OVER ALL SETUPS                                                 ==
!!$!     ==========================================================================
!!$      OPEN(NFIL,FILE='STP.CNTL_AUTOMATIC_NEW')
!!$      OPEN(NFIL2,FILE='TEST.STRC_AUTOMATIC_NEW')
!!$      REWIND NFIL
!!$      REWIND NFIL2
!!$      WRITE(NFIL,FMT='("!SCNTL")')
!!$      DO IZ=0,105
!!$        DO ISTPTYPE=1,NSTPTYPES
!!$          CALL PERIODICTABLE$GET(IZ,'SYMBOL',EL)
!!$          ID=TRIM(EL)//'_'//TRIM(STPTYPE(ISTPTYPE))
!!$!
!!$          IF(STPTYPE(ISTPTYPE).EQ.'NDLSS_V0') THEN
!!$            CALL SETUP_BUILDPARMSONE_NDLSS(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
!!$     &              ,RCL &
!!$     &              ,POTPOW,POTRC,TPOTVAL,POTVAL &
!!$     &              ,CORPOW,CORRC,TCORVAL,CORVAL &
!!$     &              ,DMIN,DMAX,RMAX)
!!$            LAMBDA=0.D0   !NOT USED FOR NDLSS CONSTRUCTION
!!$          ELSE IF(STPTYPE(ISTPTYPE).EQ.'NDLSS_SC_V0') THEN
!!$            CALL SETUP_BUILDPARMSONE_NDLSS(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
!!$     &              ,RCL &
!!$     &              ,POTPOW,POTRC,TPOTVAL,POTVAL &
!!$     &              ,CORPOW,CORRC,TCORVAL,CORVAL &
!!$     &              ,DMIN,DMAX,RMAX)
!!$            LAMBDA=0.D0   !NOT USED FOR NDLSS CONSTRUCTION
!!$          ELSE IF(STPTYPE(ISTPTYPE).EQ.'HBS') THEN
!!$            CALL SETUP_BUILDPARMSONE_HBS(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
!!$     &              ,RCL,LAMBDA &
!!$     &              ,POTPOW,POTRC,TPOTVAL,POTVAL &
!!$     &              ,CORPOW,CORRC,TCORVAL,CORVAL &
!!$     &              ,DMIN,DMAX,RMAX)
!!$          ELSE IF(STPTYPE(ISTPTYPE).EQ.'HBS_SC') THEN
!!$            CALL SETUP_BUILDPARMSONE_HBS(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
!!$     &              ,RCL,LAMBDA &
!!$     &              ,POTPOW,POTRC,TPOTVAL,POTVAL &
!!$     &              ,CORPOW,CORRC,TCORVAL,CORVAL &
!!$     &              ,DMIN,DMAX,RMAX)
!!$          ELSE IF(STPTYPE(ISTPTYPE).EQ.'.75_6.0') THEN
!!$            CALL SETUP_BUILDPARMSONE_7560(ID,LX,AEZ,ZV,COREID,TYPE,RBOX,RCSM &
!!$     &              ,RCL,LAMBDA &
!!$     &              ,POTPOW,POTRC,TPOTVAL,POTVAL &
!!$     &              ,CORPOW,CORRC,TCORVAL,CORVAL &
!!$     &              ,DMIN,DMAX,RMAX)
!!$          ELSE
!!$            CALL ERROR$MSG('SETUP TYPE NOT RECOGNIZED')
!!$            CALL ERROR$CHVAL('STPTYPE',STPTYPE(ISTPTYPE))
!!$            CALL ERROR$CHVAL('ID',ID)
!!$            CALL ERROR$STOP('SETUP_BUILDPARMS')
!!$          END IF
!!$!
!!$!         ======================================================================
!!$!         == WRITE SPECIESBLOCKS FOR SETUPFILE                                ==
!!$!         ======================================================================
!!$!          CALL SETUP_WRITESPECIESPARMSET(NFIL,ID,LX,AEZ,ZV,TYPE,RBOX,RCSM &
!!$          CALL SETUP_WRITEPARMSET(NFIL,ID,LX,AEZ,ZV,TYPE,RBOX,RCSM &
!!$     &                           ,RCL,LAMBDA &
!!$     &                           ,POTPOW,POTRC,TPOTVAL,POTVAL &
!!$     &                           ,CORPOW,CORRC,TCORVAL,CORVAL &
!!$     &                           ,DMIN,DMAX,RMAX)
!!$!
!!$!         ======================================================================
!!$!         == WRITE TO MODEL STRUCTURE INPUT FILE FOR TEST                     ==
!!$!         ======================================================================
!!$          CALL SETUP_MKSPECIESINPUT(NFIL2,ID)
!!$        ENDDO
!!$      ENDDO
!!$      WRITE(NFIL,FMT='("!END")') 
!!$      WRITE(NFIL,FMT='("!EOB")') 
!!$      WRITE(NFIL,*)
!!$      CLOSE(NFIL)
!!$      CLOSE(NFIL2)
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_WRITEPARMSET(NFIL,ID,LX,AEZ,ZV,TYPE,RBOX,RCSM &
     &                             ,RCL,LAMBDA &
     &                             ,POTPOW,POTRC,TPOTVAL,POTVAL &
     &                             ,CORPOW,CORRC,TCORVAL,CORVAL &
     &                             ,DMIN,DMAX,RMAX)
!     **************************************************************************
!     ** WRITES THE PARAMETERS FOR A SPECIFIC SETUP FOR THE SETUP INPUT FILE  **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2014******************
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: ID       ! SETUP ID
      REAL(8)     ,INTENT(IN) :: AEZ      ! ELEMENT SYMBOL
      REAL(8)     ,INTENT(IN) :: ZV       ! #(VALENCE ELECTRONS)
      CHARACTER(*),INTENT(IN) :: TYPE     ! PSEUDIZATION TYPE
      REAL(8)     ,INTENT(IN) :: RBOX     ! RBOX
      REAL(8)     ,INTENT(IN) :: RCSM     ! RCSM
      INTEGER(4)  ,INTENT(IN) :: LX       ! HIGHEST ANGULAR MOMENTUM ON FILE
      REAL(8)     ,INTENT(IN) :: RCL(LX+1)!RC(PHI)
      REAL(8)     ,INTENT(IN) :: LAMBDA(LX+1) !LAMBDA
      REAL(8)     ,INTENT(IN) :: POTPOW   ! LEADING POWER FOR PSEUDO POTENTIAL
      REAL(8)     ,INTENT(IN) :: POTRC    ! RC FOR PSEUDO POTENTIAL
      LOGICAL(4)  ,INTENT(IN) :: TPOTVAL  ! CORVAL SHALL BE SET
      REAL(8)     ,INTENT(IN) :: POTVAL   ! PSEUD POT VALUE AT R=0
      REAL(8)     ,INTENT(IN) :: CORPOW   ! LEADING POWER FOR PSEUDO CORE
      REAL(8)     ,INTENT(IN) :: CORRC    ! RC FOR PSEUDO CORE
      LOGICAL(4)  ,INTENT(IN) :: TCORVAL  ! CORVAL SHALL BE SET
      REAL(8)     ,INTENT(IN) :: CORVAL   ! PSEUD CORE VALUE AT R=0
      REAL(8)     ,INTENT(IN) :: DMIN     ! SMALLEST GRID SPACING
      REAL(8)     ,INTENT(IN) :: DMAX     ! LARGEST GRID SPACING
      REAL(8)     ,INTENT(IN) :: RMAX     ! OUTERMOST GRID POINT 
      REAL(8)                 :: RCOV     ! COVALENT RADIUS
!     **************************************************************************
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
!
      WRITE(NFIL,FMT='(T2,"!AUGMENT  ID=",A," Z=",F10.5," ZV=",F3.0)') &
     &                            "'"//TRIM(ID)//"'",AEZ,ZV
      WRITE(NFIL,FMT='(T10,"TYPE=",A," RBOX/RCOV=",F7.3," RCSM/RCOV=",F7.3)') &
     &                              "'"//TRIM(TYPE)//"'",RBOX/RCOV,RCSM/RCOV
!
!     == SET PARAMETERS FOR PARTIAL WAVE PSEUDIZATION ==========================
      WRITE(NFIL,FMT='(T10,"RCL/RCOV=",20F7.3)')RCL/RCOV
      WRITE(NFIL,FMT='(T10,"LAMBDA=",20F7.3)')LAMBDA 
!
!     == SET RADIAL GRID =======================================================
      WRITE(NFIL,FMT='(T4,"!GRID DMIN=",E12.3," DMAX=",F10.3," RMAX=",F10.3' &
     &                      //'," !END")')DMIN,DMAX,RMAX 
!
!     == SET PARAMETERS FOR PSEUDO POTENTIAL ===================================
      IF(TPOTVAL) THEN
        WRITE(NFIL,FMT='(T4,"!POT  POW=",F7.3," VAL0=",F12.3' &
     &               //'," RC/RCOV=",F7.3," !END")') &
     &                             POTPOW,POTVAL,POTRC/RCOV
      ELSE
        WRITE(NFIL,FMT='(T4,"!POT  POW=",F7.3," RC/RCOV=",F7.3," !END")') &
     &                             POTPOW,POTRC/RCOV
      END IF
!
!     == SET PARAMETERS FOR PSEUDO CORE ========================================
      IF(TPOTVAL) THEN
         WRITE(NFIL,FMT='(T4,"!CORE POW=",F7.3," VAL0=",F12.3' &
     &                //'," RC/RCOV=",F7.3," !END")') &
     &                             CORPOW,CORVAL,CORRC/RCOV
      ELSE
         WRITE(NFIL,FMT='(T4,"!CORE POW=",F7.3," RC/RCOV=",F7.3," !END")') &
     &                             CORPOW,CORRC/RCOV
      END IF
      WRITE(NFIL,FMT='(T2,"!END")') 
      WRITE(NFIL,*) 
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_WRITESPECIESPARMSET(NFIL,ID,LX,AEZ,ZV,TYPE,RBOX,RCSM &
     &                             ,RCL,LAMBDA &
     &                             ,POTPOW,POTRC,TPOTVAL,POTVAL &
     &                             ,CORPOW,CORRC,TCORVAL,CORVAL &
     &                             ,DMIN,DMAX,RMAX)
!     **************************************************************************
!     ** WRITES THE PARAMETERS FOR A SPECIFIC SETUP FOR THE SETUP INPUT FILE  **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2014******************
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: ID       ! SETUP ID
      REAL(8)     ,INTENT(IN) :: AEZ      ! ELEMENT SYMBOL
      REAL(8)     ,INTENT(IN) :: ZV       ! #(VALENCE ELECTRONS)
      CHARACTER(*),INTENT(IN) :: TYPE     ! PSEUDIZATION TYPE
      REAL(8)     ,INTENT(IN) :: RBOX     ! RBOX
      REAL(8)     ,INTENT(IN) :: RCSM     ! RCSM
      INTEGER(4)  ,INTENT(IN) :: LX       ! HIGHEST ANGULAR MOMENTUM ON FILE
      REAL(8)     ,INTENT(IN) :: RCL(LX+1)!RC(PHI)
      REAL(8)     ,INTENT(IN) :: LAMBDA(LX+1) !LAMBDA
      REAL(8)     ,INTENT(IN) :: POTPOW   ! LEADING POWER FOR PSEUDO POTENTIAL
      REAL(8)     ,INTENT(IN) :: POTRC    ! RC FOR PSEUDO POTENTIAL
      LOGICAL(4)  ,INTENT(IN) :: TPOTVAL  ! CORVAL SHALL BE SET
      REAL(8)     ,INTENT(IN) :: POTVAL   ! PSEUD POT VALUE AT R=0
      REAL(8)     ,INTENT(IN) :: CORPOW   ! LEADING POWER FOR PSEUDO CORE
      REAL(8)     ,INTENT(IN) :: CORRC    ! RC FOR PSEUDO CORE
      LOGICAL(4)  ,INTENT(IN) :: TCORVAL  ! CORVAL SHALL BE SET
      REAL(8)     ,INTENT(IN) :: CORVAL   ! PSEUD CORE VALUE AT R=0
      REAL(8)     ,INTENT(IN) :: DMIN     ! SMALLEST GRID SPACING
      REAL(8)     ,INTENT(IN) :: DMAX     ! LARGEST GRID SPACING
      REAL(8)     ,INTENT(IN) :: RMAX     ! OUTERMOST GRID POINT 
      REAL(8)                 :: RCOV     ! COVALENT RADIUS
      CHARACTER(2)            :: EL
      LOGICAL(4)              :: TSC
      CHARACTER(128)          :: STRING
      CHARACTER(1)            :: ATOMTYPE
      INTEGER(4)              :: IZ
!     **************************************************************************
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
      EL=ID(1:2)
      TSC=INDEX(ID,'_SC').NE.0
      CALL PERIODICTABLE$GET(EL,'Z',IZ)
      CALL SETUP_ELEMENTTYPE(IZ,ATOMTYPE)

      IF(.NOT.TSC) THEN
        IF(IZ.LE.2) THEN
          STRING='" NPRO=1  LRHOX=4 ")'
        ELSE IF(IZ.LE.4) THEN
          STRING='" NPRO=1 1 0  LRHOX=4 ")'
        ELSE IF(IZ.LE.10) THEN
          STRING='" NPRO=1 1 1  LRHOX=4 ")'
        ELSE IF(IZ.LE.18) THEN
          STRING='" NPRO=2 2 1  LRHOX=4 ")'
        ELSE IF(IZ.LT.56) THEN
          STRING='" NPRO=1 1 1  LRHOX=4 ")'
        ELSE      
          STRING='" NPRO=1 1 1 1 LRHOX=4 ")'
        END IF
      ELSE
        IF(IZ.LE.2) THEN
          STRING='" NPRO=1  LRHOX=4 ")'
        ELSE IF(IZ.LE.18) THEN
          STRING='" NPRO=2 2 0  LRHOX=4 ")'
        ELSE IF(IZ.LT.56) THEN
          STRING='" NPRO=2 2 1  LRHOX=4 ")'
        ELSE      
          STRING='" NPRO=2 2 1 1 LRHOX=4 ")'
        END IF
      END IF
      STRING='(T3,"!SPECIES NAME=",A,'//TRIM(STRING)
      WRITE(NFIL,FMT=TRIM(STRING))"'"//EL//"'"
!
      WRITE(NFIL,FMT='(T5,"!AUGMENT  ID=",A," Z=",F10.5," ZV=",F3.0)') &
     &                            "'"//TRIM(ID)//"'",AEZ,ZV
      WRITE(NFIL,FMT='(T14,"TYPE=",A," RBOX/RCOV=",F7.3," RCSM/RCOV=",F7.3)') &
     &                              "'"//TRIM(TYPE)//"'",RBOX/RCOV,RCSM/RCOV
!
!     == SET PARAMETERS FOR PARTIAL WAVE PSEUDIZATION ==========================
      WRITE(NFIL,FMT='(T14,"RCL/RCOV=",20F7.3)')RCL/RCOV
      WRITE(NFIL,FMT='(T14,"LAMBDA=",20F7.3)')LAMBDA 
!
!     == SET RADIAL GRID =======================================================
      WRITE(NFIL,FMT='(T7,"!GRID DMIN=",E12.3," DMAX=",F10.3," RMAX=",F10.3' &
     &                      //'," !END")')DMIN,DMAX,RMAX 
!
!     == SET PARAMETERS FOR PSEUDO POTENTIAL ===================================
      IF(TPOTVAL) THEN
        WRITE(NFIL,FMT='(T7,"!POT  POW=",F7.3," VAL0=",F12.3' &
     &               //'," RC/RCOV=",F7.3," !END")') &
     &                             POTPOW,POTVAL,POTRC/RCOV
      ELSE
        WRITE(NFIL,FMT='(T7,"!POT  POW=",F7.3," RC/RCOV=",F7.3," !END")') &
     &                             POTPOW,POTRC/RCOV
      END IF
!
!     == SET PARAMETERS FOR PSEUDO CORE ========================================
      IF(TPOTVAL) THEN
         WRITE(NFIL,FMT='(T7,"!CORE POW=",F7.3," VAL0=",F12.3' &
     &                //'," RC/RCOV=",F7.3," !END")') &
     &                             CORPOW,CORVAL,CORRC/RCOV
      ELSE
         WRITE(NFIL,FMT='(T7,"!CORE POW=",F7.3," RC/RCOV=",F7.3," !END")') &
     &                             CORPOW,CORRC/RCOV
      END IF
      WRITE(NFIL,FMT='(T5,"!END")') 
      WRITE(NFIL,FMT='(T3,"!END")') 
      WRITE(NFIL,*) 
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_ELEMENTTYPE(IZ,ELEMENTTYPE)
!     **************************************************************************
!     ** CLASSIFIES THE ATOMS ACCORDING TO SPECIAL ELEMENT TYPES           **
!     **    A  : ALKALI OR EARTH ALKALI ATOM                                  **
!     **    M  : MAIN GROUP ELEMENT                                           **
!     **    T  : TRANSITION METAL ELEMENT                                     **
!     **    F  : F-ELECTRON METAL: LANTHANIDE OR ACTINIDE                     **
!     ******************************PETER BLOECHL, GOSLAR 2014******************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: IZ              ! ATOMIC NUMBER
      CHARACTER(1),INTENT(OUT) :: ELEMENTTYPE     ! ELEMENT TYPE
!     **************************************************************************
      ELEMENTTYPE=' '
      IF(IZ.GE.  1.AND.IZ.LE. 2) ELEMENTTYPE='A'   ! H,HE
      IF(IZ.GE.  3.AND.IZ.LE. 4) ELEMENTTYPE='A'   ! ALKALI AND EARTH ALKALI
      IF(IZ.GE.  5.AND.IZ.LE.10) ELEMENTTYPE='M'   ! MAIN GROUP
      IF(IZ.GE. 11.AND.IZ.LE.12) ELEMENTTYPE='A'   ! LI, BE
      IF(IZ.GE. 13.AND.IZ.LE.18) ELEMENTTYPE='M'   ! OTHER 1-ST GROUP ELEMENTS
      IF(IZ.GE. 19.AND.IZ.LE.20) ELEMENTTYPE='A'   ! K, CA
      IF(IZ.GE. 21.AND.IZ.LE.30) ELEMENTTYPE='T'   ! 3D TRANSITION METALS
      IF(IZ.GE. 31.AND.IZ.LE.36) ELEMENTTYPE='M'   
      IF(IZ.GE. 37.AND.IZ.LE.38) ELEMENTTYPE='A'   ! RB, SR
      IF(IZ.GE. 39.AND.IZ.LE.48) ELEMENTTYPE='T'   ! 4D TRANSITION METALS
      IF(IZ.GE. 49.AND.IZ.LE.54) ELEMENTTYPE='M'   
      IF(IZ.GE. 55.AND.IZ.LE.56) ELEMENTTYPE='A'   ! CS, BA
      IF(IZ.GE. 57.AND.IZ.LE.71) ELEMENTTYPE='F'   ! LANTHANIDES
      IF(IZ.GE. 72.AND.IZ.LE.80) ELEMENTTYPE='T'   ! 5D TRANSITION METALS
      IF(IZ.GE. 81.AND.IZ.LE.86) ELEMENTTYPE='M'   
      IF(IZ.GE. 87.AND.IZ.LE.88) ELEMENTTYPE='A'   ! FR, RD
      IF(IZ.GE. 89.AND.IZ.LE.103)ELEMENTTYPE='F'   ! ACTINIDES
      IF(IZ.GE.104.AND.IZ.LE.112)ELEMENTTYPE='T'   ! 6D TRANSITION METALS
      IF(IZ.GE.113.AND.IZ.LE.118)ELEMENTTYPE='M'   
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_RADNDLSSMAX(AEZ,N,RAD)
!     **************************************************************************
!     ** ESTIMATE THE POSITION OF THE MAXIMUM OF THE NODELESS WAVE FUNCTION   **
!     ** FOR A GIVEN ATOMIC NUMBER AND MAIN QUANTUM NUMBER N.                 **
!     **                                                                      **
!     ** THE POSITION OF THE OUTERMOST MAXIMUM HAS BEEN CALCULATED FOR ALL    **
!     ** ATOMS IN THE PERIODIC TABLE USING PBE0, AND SCALAR RELATIVISTIC      **
!     ** THE VALUES WITH 0.1 A.U.<RAD<5 A.U. HAVE BEEN FITTED AS LINEAR       **
!     ** FUNCTION LN(RAD)=A+B*LN(AEZ).                                        **
!     **                                                                      **
!     ** SET (1) CONSIDERS THE MAXIMUM OF |U(R)|                              **
!     ** SET (2) CONSIDERS THE MAXIMUM OF |R*U(R)|                            **
!     **                                                                      **
!     ** THE 1S FUNCTION HAS ITS MAXIMUM AT R=0 AND DOES NOT GIVE A MEANINGFUL**
!     ** FIT.                                                                 **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2013******************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: AEZ  ! ATOMIC NUMBER
      INTEGER(4),INTENT(IN)  :: N    ! MAIN QUANTUM NUMBER
      REAL(8)   ,INTENT(OUT) :: RAD  ! RADIUS OF MAXIMAL VALIE
      INTEGER(4),PARAMETER   :: NX=7
      REAL(8)   ,PARAMETER   :: A1(NX)=(/0.D0,1.4027D0,3.5632D0,5.6558D0 &
     &                                  ,6.9208D0,9.1717D0,9.9144D0/)
      REAL(8)   ,PARAMETER   :: B1(NX)=(/0.D0,-1.1004D0,-1.3410D0,-1.5993D0 &
     &                                  ,-1.6762D0,-1.9976D0,-1.9896D0/)
      REAL(8)   ,PARAMETER   :: A2(NX)=(/0.072414D0,2.3212D0,4.1957D0,6.2388D0 &
     &                                  ,7.6373D0,10.068D0,10.03D0/)
      REAL(8)   ,PARAMETER   :: B2(NX)=(/-1.0164D0,-1.2082D0,-1.4017D0 &
     &                                ,-1.6617D0,-1.7762D0,-2.1334D0,-1.9469D0/)
!     **************************************************************************
      IF(N.LT.1.OR.N.GT.7) THEN
        CALL ERROR$MSG('VALUE FOR MAIN QUANTUM NUMBER OUT OF RANGE')
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$STOP('SETUP_RADNDLSSMAX')
      END IF
      IF(AEZ.LE.0.D0) THEN
        CALL ERROR$MSG('VALUE FOR ATOMIC NUMBER OUT OF RANGE')
        CALL ERROR$R8VAL('AEZ',AEZ)
        CALL ERROR$STOP('SETUP_RADNDLSSMAX')
      END IF

!     RAD=EXP(A1(N)+B1(N)*LOG(AEZ))  ! MAX(U(R))  N=1 NOT AVAILABLE
      RAD=EXP(A2(N)+B2(N)*LOG(AEZ))  ! MAX(R*U(R))
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_MKSPECIESINPUT(NFIL,ID)
!     **************************************************************************
!     ** GENERATES AUTOMATICALLY A PAERAMETER FILE FOR THE SETUP CONSTRUCTION **
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(2)            :: EL
      CHARACTER(1)            :: ATOMTYPE
      LOGICAL(4)              :: TSC   ! SEMI-CORE?
      CHARACTER(64)           :: STRING
      INTEGER(4)              :: IZ
!     **************************************************************************
      EL=ID(1:2)
      TSC=INDEX(ID,'_SC').NE.0
      CALL PERIODICTABLE$GET(EL,'Z',IZ)
      CALL SETUP_ELEMENTTYPE(IZ,ATOMTYPE)
      WRITE(NFIL,FMT='(T2,"!SPECIES NAME=",A," ID=",A)') &
     &                               "'"//EL//"'","'"//TRIM(ID)//"'"
      IF(.NOT.TSC) THEN
        IF(IZ.LE.2) THEN
          WRITE(NFIL,FMT='(T4,"NPRO=1  LRHOX=4")')
        ELSE IF(IZ.LE.18) THEN
          WRITE(NFIL,FMT='(T4,"NPRO=1 1 0  LRHOX=4")')
        ELSE IF(IZ.LT.56) THEN
          WRITE(NFIL,FMT='(T4,"NPRO=1 1 1  LRHOX=4")')
        ELSE      
          WRITE(NFIL,FMT='(T4,"NPRO=1 1 1 1 LRHOX=4")')
        END IF
      ELSE
        IF(IZ.LE.2) THEN
          WRITE(NFIL,FMT='(T4,"NPRO=1  LRHOX=4")')
        ELSE IF(IZ.LE.18) THEN
          WRITE(NFIL,FMT='(T4,"NPRO=2 2 0  LRHOX=4")')
        ELSE IF(IZ.LT.56) THEN
          WRITE(NFIL,FMT='(T4,"NPRO=2 2 1  LRHOX=4")')
        ELSE      
          WRITE(NFIL,FMT='(T4,"NPRO=2 2 1 1 LRHOX=4")')
        END IF
      END IF
!
!     ==========================================================================
!     == ENTER THE HYBRID BLOCK 
!     ==========================================================================
      IF(.NOT.TSC) THEN
        IF(ATOMTYPE.EQ.'A') THEN
          STRING='NCORROFL= 1 0 0 0'
        ELSE IF(ATOMTYPE.EQ.'M') THEN
          STRING='NCORROFL= 1 2 0 0'
        ELSE IF(ATOMTYPE.EQ.'T') THEN
          STRING='NCORROFL= 1 0 1 0'
        ELSE IF(ATOMTYPE.EQ.'F') THEN
          STRING='NCORROFL= 1 0 1 1'
        END IF
      ELSE
        IF(ATOMTYPE.EQ.'A') THEN
          STRING='NCORROFL= 2 1 0 0'
        ELSE IF(ATOMTYPE.EQ.'M') THEN
          STRING='NCORROFL= 2 2 0 0'
        ELSE IF(ATOMTYPE.EQ.'T') THEN
          STRING='NCORROFL= 2 1 1 0'
        ELSE IF(ATOMTYPE.EQ.'F') THEN
          STRING='NCORROFL= 2 1 1 1'
        END IF
      END IF
      WRITE(NFIL,FMT='(T4,"!HYBRID_X ",A," CV=T HFWEIGHT=0.25 !END")')
!
      WRITE(NFIL,FMT='(T2,"!END")')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_CVXSETUP()
!     **************************************************************************
!     **  CORE-VALENCE EXCHANGE MATRIX ELEMENTS                               **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE SETUP_MODULE,ONLY : THIS
      IMPLICIT NONE
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: LMRX
      INTEGER(4)            :: LRHOX
      INTEGER(4)            :: NC
      INTEGER(4)            :: LNX
!     **************************************************************************
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      LMRX=THIS%LMRX
      LRHOX=INT(SQRT(REAL(LMRX-1)+1.D-8))
      NC=THIS%ATOM%NC
      LNX=THIS%LNX
      ALLOCATE(THIS%COREVALENCEX(LNX,LNX))
      CALL SETUP_CVXMAT(GID,NR &
     &              ,LNX,THIS%LOX,THIS%AEPHI &
     &              ,NC,THIS%ATOM%LOFI(:NC),THIS%ATOM%AEPSI(:,:NC) &
     &              ,LRHOX,THIS%COREVALENCEX)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_CVXMAT(GID,NR,LNX,LOX,AEPHI,NC,LOFC,PSIC,LRHOX,MAT)
!     **************************************************************************
!     **  CORE-VALENCE EXCHANGE MATRIX ELEMENTS                               **
!     **                                                                      **
!     ** NOTE THAT THERE IS A SIMILAR ROUTINE IN PAW_LDAPLUSU                 **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2008******************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID           ! GRID ID
      INTEGER(4),INTENT(IN) :: NR            ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: LNX           ! #(PARTIAL WAVES W/O M,SIGMA)
      INTEGER(4),INTENT(IN) :: LOX(LNX)      ! ANGULAR MOMENTA OF PARTIAL WAVES
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX) ! AE PARTIAL WAVES
      INTEGER(4),INTENT(IN) :: NC            ! #(CORE STATES)
      INTEGER(4),INTENT(IN) :: LOFC(NC)      ! ANGULAR MOMENTA OF CORE STATES
      REAL(8)   ,INTENT(IN) :: PSIC(NR,NC)   ! CORE STATES
      INTEGER(4),INTENT(IN) :: LRHOX         ! MAX ANGULAR MOMENTUM OF DENSITY
      REAL(8)   ,INTENT(OUT):: MAT(LNX,LNX)  ! CORE VALENCE X MATRIX ELEMENTS
      LOGICAL(4),PARAMETER  :: TPRINT=.FALSE.
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)            :: LX  ! MAX ANGULAR MOMENTUM PARTIAL WAVES
      INTEGER(4)            :: LMX ! MAX #(ANGULAR MOMENTA) OF PARTIAL WAVES
      INTEGER(4)            :: LCX ! HIGHEST ANGULAR MOMENTUM OF CORE STATES
      REAL(8)               :: CG
      REAL(8)               :: RHO1(NR)
      REAL(8)               :: AUX(NR),POT(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VAL
      REAL(8)   ,ALLOCATABLE:: FACTOR(:,:,:)
      INTEGER(4)            :: LMCA
      INTEGER(4)            :: LM1,LC,LRHO,IMC,LMC,LMRHOA,IMRHO,LMRHO
      INTEGER(4)            :: LN1,L1,IC,LN2,L2
      REAL(8)               :: CNORM
!     **************************************************************************
      LX=MAXVAL(LOX)
      LMX=(LX+1)**2
      LCX=MAXVAL(LOFC)
      ALLOCATE(FACTOR(LX+1,LRHOX+1,LCX+1))
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == INCLUDE ANGULAR INTEGRATIONS                                         ==
!     ==========================================================================
!     == SPHERICAL SYMMETRY EXPLOITED: 
      FACTOR=0.D0
      DO L1=0,LX           ! ANGULAR MOMENTUM OF LOCAL ORBITAL
        LM1=L1**2+1
        DO LC=0,LCX        ! ANGULAR MOMENTUM OF CORE STATE
          DO LRHO=0,LRHOX  ! ANGULAR MOMENTUM OF DENSITY
            LMCA=LC**2
            DO IMC=1,2*LC+1
              LMC=LMCA+IMC
              LMRHOA=LRHO**2
              DO IMRHO=1,2*LRHO+1
                LMRHO=LMRHOA+IMRHO
                CALL SPHERICAL$GAUNT(LM1,LMC,LMRHO,CG)
                FACTOR(L1+1,LRHO+1,LC+1)=FACTOR(L1+1,LRHO+1,LC+1)+CG**2
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     = THIS FACTOR IS INCLUDED TO HAVE A PROPER DEFINITION OF SLATER INTEGRALS
      DO LRHO=0,LRHOX
        FACTOR(:,LRHO+1,:)=FACTOR(:,LRHO+1,:)*4.D0*PI/REAL(2*LRHO+1,KIND=8)
      ENDDO
!
!     ==========================================================================
!     == WORK OUT RADIAL INTEGRATIONS                                         ==
!     ==========================================================================
      MAT(:,:)=0.D0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IC=1,NC
          LC=LOFC(IC)
          RHO1(:)=AEPHI(:,LN1)*PSIC(:,IC)
          DO LRHO=0,LRHOX
            CALL RADIAL$POISSON(GID,NR,LRHO,RHO1,POT)
            DO LN2=1,LN1  ! MATRIX IS SYMMETRIC
              L2=LOX(LN2)
              IF(L2.NE.L1) CYCLE
              AUX(:)=R(:)**2*POT(:)*PSIC(:,IC)*AEPHI(:,LN2)
              CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!
!!$PRINT*,LN1,IC,LC,LRHO,LN2,'VAL ',VAL
!!$CALL ATOMLIB_WRITEPHI('X0.DAT',GID,NR,1,AUX)
              VAL=VAL*REAL(2*LRHO+1,KIND=8)/(4.D0*PI)  !SLATER INTEGRAL
              MAT(LN1,LN2)=MAT(LN1,LN2)-FACTOR(L1+1,LRHO+1,LC+1)*VAL
              MAT(LN2,LN1)=MAT(LN1,LN2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(FACTOR)
!
!     ==========================================================================
!     ==  WRITE FOR TEST PURPOSES                                             ==
!     ==========================================================================
      IF(TPRINT) THEN
        PRINT*,'LOFC ',LOFC
        PRINT*,'NOW THE MATRIX FOR THE CORE VALENCE EXCHANGE INTERACTION'
        WRITE(*,FMT='(4A3,10A18)')'LN1','L1','LN2','L2','MAT(LN1,LN2)'
        DO LN1=1,LNX
          L1=LOX(LN1)
          DO LN2=1,LNX
            L2=LOX(LN2)
            WRITE(*,FMT='(4I3,10F18.10)')LN1,L1,LN2,L2,MAT(LN1,LN2)
          ENDDO
        ENDDO
!        STOP
!        CALL TESTPOTTEST(GID,NR,LNX,LOX,LOFC,NC,LRHOX,AEPHI,PSIC)
      END IF

      RETURN
      END      
! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TESTPOTTEST(GID,NR,LNX,LOX,LOFC,NC,LRHOX,AEPHI,PSIC)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID           ! GRID ID
      INTEGER(4),INTENT(IN) :: NR            ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: LNX           ! #(PARTIAL WAVES W/O M,SIGMA)
      INTEGER(4),INTENT(IN) :: LOX(LNX)      ! ANGULAR MOMENTA OF PARTIAL WAVES
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX) ! AE PARTIAL WAVES
      INTEGER(4),INTENT(IN) :: NC            ! #(CORE STATES)
      INTEGER(4),INTENT(IN) :: LOFC(NC)      ! ANGULAR MOMENTA OF CORE STATES
      REAL(8)   ,INTENT(IN) :: PSIC(NR,NC)   ! CORE STATES
      INTEGER(4),INTENT(IN) :: LRHOX         ! X(ANGULAR MOMENMTUM OF THE DENISTY)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,ALLOCATABLE:: RHO(:)
      REAL(8)   ,ALLOCATABLE:: POT(:,:)
      INTEGER(4),SAVE       :: COUNTER=0
      CHARACTER(128)        :: NAME
      INTEGER(4)            :: LN1,L1,IC,LC,LRHO,IND
      REAL(8)               :: R(NR)
      REAL(8)               :: RHO1(NR)
      CHARACTER(128)        :: STRING
      REAL(8)               :: VAL
!     **************************************************************************
      COUNTER=COUNTER+1
      CALL RADIAL$R(GID,NR,R)
      ALLOCATE(POT(NR,LRHOX+1))
      ALLOCATE(RHO(NR))
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IC=1,NC
          LC=LOFC(IC)
          RHO(:)=AEPHI(:,LN1)*PSIC(:,IC)
IF(LC.EQ.L1) THEN
  CALL RADIAL$INTEGRAL(GID,NR,R**2*RHO,VAL)
 PRINT*,'OVERLAP  IAT=',COUNTER,' LN=',LN1,' O=',VAL
END IF
          DO LRHO=0,LRHOX
            CALL RADIAL$POISSON(GID,NR,LRHO,RHO,POT(:,LRHO+1))
          ENDDO
          NAME='POT_'
          WRITE(STRING,*)COUNTER
          NAME=TRIM(ADJUSTL(NAME))//"_"//TRIM(ADJUSTL(STRING))
          WRITE(STRING,*)LN1
          NAME=TRIM(ADJUSTL(NAME))//"_"//TRIM(ADJUSTL(STRING))
          WRITE(STRING,*)IC
          NAME=TRIM(ADJUSTL(NAME))//"_"//TRIM(ADJUSTL(STRING))
          CALL SETUP_WRITEPHI(TRIM(NAME),GID,NR,LRHOX+1,POT)
        ENDDO
      ENDDO
STOP
      RETURN
      END
    

! 
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_OUTEROLDPROWRAPPER(GID,NR,ROUT,RBOX,RCOV &
     &            ,NC,NB,LOFI,SOFI,EOFI,LNX,LOX,TYPE,RC,LAMBDA,ISCATT &
     &            ,AEPOT,PSPOT,VFOCK &
     &            ,QN,AEPHI,PSPHI,PRO,DTKIN,DOVER,AEPHIDOT,PSPHIDOT,TREL,TZORA &
     &            ,PSISCALE,PHISCALE)
      USE RADIALFOCK_MODULE, ONLY: VFOCK_TYPE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: ROUT
      REAL(8)   ,INTENT(IN) :: RBOX
      REAL(8)   ,INTENT(IN) :: RCOV
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NC 
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      INTEGER(4),INTENT(IN) :: SOFI(NB)
      REAL(8)   ,INTENT(IN) :: EOFI(NB)
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      CHARACTER(*),INTENT(IN) :: TYPE ! PSEUDIZATION METHOD
      REAL(8)   ,INTENT(IN) :: RC(LNX)
      REAL(8)   ,INTENT(IN) :: LAMBDA(LNX)
      INTEGER(4),INTENT(IN) :: ISCATT(LNX)
      REAL(8)   ,INTENT(IN) :: AEPOT(NR)
      TYPE(VFOCK_TYPE),INTENT(IN) :: VFOCK
      REAL(8)   ,INTENT(IN) :: PSPOT(NR)
      REAL(8)   ,INTENT(OUT):: QN(NR,LNX)
      REAL(8)   ,INTENT(OUT):: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PRO(NR,LNX)
      REAL(8)   ,INTENT(OUT):: DTKIN(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: DOVER(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: AEPHIDOT(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PSPHIDOT(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PHISCALE(LNX)
      REAL(8)   ,INTENT(OUT):: PSISCALE(NB)
      LOGICAL(4),INTENT(IN) :: TREL
      LOGICAL(4),INTENT(IN) :: TZORA
      LOGICAL(4),PARAMETER  :: TSMALLBOX=.FALSE.
      LOGICAL(4),PARAMETER  :: TMATCHTOALLELECTRON=.FALSE.
      LOGICAL(4),PARAMETER  :: TCUTTAIL=.TRUE.
      LOGICAL   ,PARAMETER  :: TTEST=.TRUE.
      LOGICAL   ,PARAMETER  :: TWRITE=.FALSE.
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      LOGICAL(4)            :: TSEQUENTIALAUGMENT=.FALSE.
      REAL(8)               :: R(NR)
      REAL(8)               :: DREL(NR)
      REAL(8)               :: EOFI1(NB)
      REAL(8)               :: G(NR)
      REAL(8)               :: GS(NR)
      REAL(8)               :: UOFI(NR,NB)
      REAL(8)               :: TUOFI(NR,NB)
      REAL(8)               :: UOFISM(NR,NB)
      REAL(8)               :: NLPHI(NR,LNX)
      REAL(8)               :: TNLPHI(NR,LNX)
      REAL(8)               :: NLPHIDOT(NR,LNX)
      REAL(8)               :: TAEPHI(NR,LNX)
      REAL(8)               :: TPSPHI(NR,LNX)
      REAL(8)               :: TQN(NR,LNX)
      REAL(8)               :: QNDOT(NR,LNX)
      REAL(8)               :: BAREPRO(NR,LNX)
      REAL(8)               :: PHI(NR)
      REAL(8)   ,ALLOCATABLE::  PHITEST(:,:),TPHITEST(:,:)
      REAL(8)               :: PHITEST1(NR,LNX)
      REAL(8)               :: QNP(NR,LNX),PSPHIP(NR,LNX)
      REAL(8)   ,ALLOCATABLE:: DH1(:,:),DO1(:,:),PRO1(:,:)
      REAL(8)               :: E
      REAL(8)               :: ESCATT(LNX)
      REAL(8)               :: RBND,RBNDOUT
      REAL(8)               :: PHIPHASE
      REAL(8)               :: AUX(NR),AUX1(NR),VAL,SVAR,SVAR1,SVAR2,RC1
      REAL(8)               :: SPEEDOFLIGHT
      REAL(8)               :: EOFLN(LNX)
      REAL(8)               :: TRANSU(LNX,LNX),TRANSUINV(LNX,LNX)
      REAL(8)               :: TRANSPHI(LNX,LNX),TRANSPHIINV(LNX,LNX)
      REAL(8)               :: TRANSPRO(LNX,LNX)
      REAL(8)   ,ALLOCATABLE:: A(:,:)
      REAL(8)               :: UBYQ(LNX)
      REAL(8)               :: DH(LNX,LNX)
      REAL(8)   ,ALLOCATABLE:: PROJ(:)
      LOGICAL(4)            :: TFIRST,TVARDREL
      INTEGER(4)            :: L,LN,LN1,LN2,ISO,IB,IB1,IR,IR1,IPRO,IPRO1,IPRO2
      INTEGER(4)            :: NPRO
      INTEGER(4)            :: LNLAST
      INTEGER(4)            :: LX ! X(ANGULAR MOMENTUM)
      INTEGER(4),ALLOCATABLE :: NCL(:)
      INTEGER(4),ALLOCATABLE :: NPROL(:)
      CHARACTER(128)         :: STRING
!     **************************************************************************
                                     CALL TRACE$PUSH('SETUP_OUTEROLDPROWRAPPER')
      CALL CONSTANTS$GET('C',SPEEDOFLIGHT)
      LX=MAX(MAXVAL(LOFI),MAXVAL(LOX))
      CALL RADIAL$R(GID,NR,R)
      RBND=RBOX
!     == RBNDOUT HAS EARLIER BEEN SET EQUAL TO RBOX. IT SHOULD BE EQUAL TO ROUT.
!     == THIS VARIABLE IS USED TO TRACK THE CHANGES.
      RBNDOUT=ROUT
!
      ALLOCATE(NCL(0:LX))
      NCL(:)=0
      DO IB=1,NC
        L=LOFI(IB)
        NCL(L)=MAX(NCL(L),IB)
      ENDDO

!     == DETERMINE NUMBER OF PROJECTORS FOR EACH ANGULAR MOMENTUM ==============
      ALLOCATE(NPROL(0:LX))
      NPROL(:)=0
      DO LN=1,LNX
        L=LOX(LN)
        NPROL(L)=NPROL(L)+1
      ENDDO
!
!     ==========================================================================
!     == CONSTRUCT NODELESS WAVE FUNCTIONS (LOCAL POTENTIAL ONLY)             ==
!     ==========================================================================
                           CALL TRACE$PASS('CONSTRUCT NODELESS WAVE FUNCTIONS')
      DREL(:)=0.D0
      TVARDREL=TREL.AND.(.NOT.TZORA)
      IF(TREL.AND.TZORA) CALL SCHROEDINGER$DREL(GID,NR,AEPOT,0.D0,DREL)
!
      EOFI1(:)=EOFI(:)
      DO L=0,LX
        DO ISO=-1,1
          G(:)=0.D0
          GS(:)=0.D0
          DO IB=1,NB
            IF(LOFI(IB).NE.L) CYCLE
            IF(SOFI(IB).NE.ISO) CYCLE
            E=EOFI1(IB)
            CALL ATOMLIB$BOUNDSTATE(GID,NR,L,ISO,0.D0,ROUT,TVARDREL &
     &                             ,DREL,G,0,AEPOT,E,UOFI(:,IB))
            IF(TREL.AND.(.NOT.TZORA)) THEN
              CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,L,ISO &
     &                                         ,DREL,GS,UOFI(:,IB),UOFISM(:,IB))
            ELSE
              UOFISM(:,IB)=0.D0
            END IF
            EOFI1(IB)=E
            TUOFI(:,IB)=G+(E-AEPOT(:)*Y0)*UOFI(:,IB)
            G(:)=UOFI(:,IB) 
            GS(:)=0.D0      
!
!           == DETERMINE THE INHOMOGENEITY CONSIDERING THE SMALL COMPONENT
!           == INCLUDING THE SMALL COMPONENT CAN FAIL BECAUSE OF NON-REGULAR
!           == NUMBER OF NODES
!!$            AUX=(1.D0+DREL)*UOFISM(:,IB)
!!$            CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
!!$            IF(ISO.EQ.1) THEN
!!$              AUX(2:)=AUX1(2:)+REAL(L+2,KIND=8)/R(2:)*AUX(2:)
!!$              AUX(1)=AUX(2)
!!$              AUX(:)=-0.5/SPEEDOFLIGHT*AUX(:)
!!$            ELSE IF(ISO.EQ.-1) THEN
!!$              AUX(2:)=AUX1(2:)-REAL(L-1,KIND=8)/R(2:)*AUX(2:)
!!$              AUX(1)=AUX(2)
!!$              AUX=+0.5/SPEEDOFLIGHT*AUX(:)
!!$            ELSE
!!$              AUX=+0.5/SPEEDOFLIGHT*AUX1
!!$            END IF
!!$            G(:)=UOFI(:,IB)+AUX(:)
!!$            GS(:)=UOFISM(:,IB)
          ENDDO
        ENDDO
      ENDDO
!!$PRINT*,'EOFI ',EOFI
!!$PRINT*,'EOFI1 ',EOFI1
      IF(TTEST.AND.TWRITE)CALL SETUP_WRITEPHI('UOFI1.DAT',GID,NR,NB,UOFI)

!
!     ==========================================================================
!     == CONSTRUCT NODELESS PARTIAL WAVES  (LOCAL POTENTIAL ONLY)             ==
!     ==========================================================================
                           CALL TRACE$PASS('CONSTRUCT NODELESS PARTIAL WAVES')
!     == FIRST USE ONLY THE LOCAL POTENTIAL... =================================
      DO L=0,LX
        ISO=0
!       == GET ENERGY OF THE LOWEST VALENCE STATE WITH THIS L ==================
        E=0.D0
!!$PRINT*,'NB,NC,EOFI1',NB,NC,EOFI1

        DO IB=NC+1,NB
          IF(LOFI(IB).NE.L) CYCLE
          E=EOFI1(IB)
          EXIT
        ENDDO
!       == USE HIGHEST CORE STATE WITH THIS L AS INHOMOGENEITY =================
        G(:)=0.D0
        IF(NCL(L).NE.0)G(:)=UOFI(:,NCL(L))
        TFIRST=.TRUE.   ! SWITCH BACK TO BOX WITH TFIRST=.FALSE.
        PHIPHASE=1.D0   ! NODE 
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          IF(TREL)THEN
            IF(TZORA) THEN
              CALL SCHROEDINGER$DREL(GID,NR,AEPOT,0.D0,DREL)
            ELSE
              CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
            END IF
          END IF
          IF(TFIRST.AND.(.NOT.TSMALLBOX)) THEN
            PHIPHASE=1.D0   ! NODE AT ROUT
            CALL ATOMLIB$PHASESHIFTSTATE(GID,NR,L,ISO,DREL,G,AEPOT &
     &                                ,ROUT,PHIPHASE,E,PHI)
            CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,0.D0,RBND,PHIPHASE)
!           == FIX SIZE OF FIRST PARTIAL WAVE ABOUT EQUAL TO UOFI   ============
            IF(NCL(L).EQ.0) THEN
              IR=1
              DO IR1=1,NR
                IF(R(IR1).GT.RBND) EXIT
                IR=IR1
              ENDDO
              DO IB1=1,NB
                IB=IB1   !BAND INDEX OF LOWEST STATE WITH THIS L
                IF(LOFI(IB1).EQ.L) EXIT
              ENDDO
              PHI(:)=PHI(:)*UOFI(IR,IB)/PHI(IR)
            END IF
            TFIRST=.FALSE.
          ELSE
            CALL ATOMLIB$PHASESHIFTSTATE(GID,NR,L,ISO,DREL,G,AEPOT &
     &                                ,RBND,PHIPHASE,E,PHI)
          END IF
          EOFLN(LN)=E
          NLPHI(:,LN)=PHI(:)
          TNLPHI(:,LN)=G(:)+(E-AEPOT(:)*Y0)*PHI(:)
          G(:)=NLPHI(:,LN)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == UPDATE WAVE FUNCTIONS WITH FOCK POTENTIAL                            ==
!     == THE FOCK CORRECTION MUST
!     ==========================================================================
      IF(VFOCK%TON) THEN
                      CALL TRACE$PASS('APPLY FOCK CORRECTION TO WAVE FUNCTIONS')
        DO L=0,LX
          DO ISO=-1,1
            G(:)=0.D0
            GS(:)=0.D0
            DO IB=1,NB
              IF(LOFI(IB).NE.L) CYCLE
              IF(SOFI(IB).NE.ISO) CYCLE
              E=EOFI1(IB)
              IF(TREL)THEN
                IF(TZORA) THEN
                  CALL SCHROEDINGER$DREL(GID,NR,AEPOT,0.D0,DREL)
                ELSE
                  CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
                END IF
              ELSE
                DREL(:)=0.D0
              END IF
!             == REMARK: IF THIS CRASHES PROCEED LIKE FOR NODELESS PARTIAL WAVES
!             == WORK WITH LOCAL POTENTIAL FIRST AND THEN UPDATE WITH FOCK POT
!!$              IF(.NOT.TSMALLBOX) THEN
!!$                CALL ATOMLIB$UPDATESTATEWITHHF(GID,NR,L,ISO,DREL,G,AEPOT,VFOCK &
!!$     &                                    ,RBND,E,UOFI(:,IB))
!!$              ELSE
                CALL ATOMLIB$UPDATESTATEWITHHF(GID,NR,L,ISO,DREL,G,AEPOT,VFOCK &
     &                                        ,ROUT,E,UOFI(:,IB))
!!$              END IF
              IF(TREL) THEN
                CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,L,ISO &
     &                                         ,DREL,GS,UOFI(:,IB),UOFISM(:,IB))
              ELSE
                UOFISM(:,IB)=0.D0
              END IF
              EOFI1(IB)=E
!             == DETERMINE KINETIC ENERGY OF UOFI: (E-POT)|UOFI> ===============
              CALL RADIALFOCK$VPSI(GID,NR,VFOCK,L,UOFI(:,IB),AUX)
              TUOFI(:,IB)=G+(E-AEPOT(:)*Y0)*UOFI(:,IB)-AUX(:)
!
!             ==  DETERMINE INHOMOGENEITY FOR THE NEXT BAND ====================
              AUX=(1.D0+DREL)*UOFISM(:,IB)
              CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
              IF(ISO.EQ.1) THEN
                AUX(2:)=AUX1(2:)+REAL(L+2,KIND=8)/R(2:)*AUX(2:)
                AUX(1)=AUX(2)
                AUX(:)=-0.5/SPEEDOFLIGHT*AUX(:)
              ELSE IF(ISO.EQ.-1) THEN
                AUX(2:)=AUX1(2:)-REAL(L-1,KIND=8)/R(2:)*AUX(2:)
                AUX(1)=AUX(2)
                AUX=+0.5/SPEEDOFLIGHT*AUX(:)
              ELSE
                AUX=+0.5/SPEEDOFLIGHT*AUX1
              END IF
              G(:) =UOFI(:,IB) !+AUX(:)
              GS(:)=UOFISM(:,IB) ! REST IS DONE IN SCHR...$SPHSMALLCOMPONENT
           ENDDO
          ENDDO
        ENDDO
PRINT*,'EOFI1 A ',EOFI1
      END IF
!
!     ==========================================================================
!     == UPDATE NODELESS PARTIAL WAVES WITH FOCK POTENTIAL                    ==
!     ==========================================================================
PRINT*,'VFOCK%TON ',VFOCK%TON
      IF(VFOCK%TON) THEN
                      CALL TRACE$PASS('APPLY FOCK CORRECTION TO PARTIAL WAVES')
PRINT*,'EOFLN BEFORE VFOCK',EOFLN
        DO L=0,LX
          TFIRST=.TRUE.
          ISO=0
          G(:)=0.D0
          IF(NCL(L).NE.0)G(:)=UOFI(:,NCL(L))
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            IF(TREL)THEN
              IF(TZORA) THEN
                CALL SCHROEDINGER$DREL(GID,NR,AEPOT,0.D0,DREL)
              ELSE
                CALL SCHROEDINGER$DREL(GID,NR,AEPOT,EOFLN(LN),DREL)
              END IF
            ELSE
              DREL(:)=0.D0
            END IF
            IF(TFIRST.AND.(.NOT.TSMALLBOX)) THEN
              CALL ATOMLIB$UPDATESTATEWITHHF(GID,NR,L,ISO,DREL,G,AEPOT,VFOCK &
    &                                    ,ROUT,EOFLN(LN),NLPHI(:,LN))
              TFIRST=.FALSE. 
            ELSE
              CALL ATOMLIB$UPDATESTATEWITHHF(GID,NR,L,ISO,DREL,G,AEPOT,VFOCK &
    &                                    ,RBNDOUT,EOFLN(LN),NLPHI(:,LN))
            END IF
            CALL RADIALFOCK$VPSI(GID,NR,VFOCK,L,NLPHI(:,LN),AUX)
            TNLPHI(:,LN)=G(:)+(EOFLN(LN)-AEPOT(:)*Y0)*NLPHI(:,LN)-AUX(:)
            G(:)=NLPHI(:,LN)
          ENDDO
        ENDDO
PRINT*,'EOFLN AFTER VFOCK',EOFLN
      END IF
!
!     ==========================================================================
!     == REPORT SETTINGS ON WAVE FUNCTIONS                                    ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(80("="),T20," ENERGIES FOR ATOMIC WAVE FUNCTIONS ")')
        WRITE(6,FMT='(80("="),T20," OLD: AE SCHRODINGER EQUATION           ")')
        WRITE(6,FMT='(80("="),T20," NEW: NODELESS EQUATION                 ")')
        WRITE(6,FMT='(80("="),T20," DIFFERENCE DUE TO RELATIVISTIC EFFECTS ")')
        DO IB=1,NB
          WRITE(6,FMT='("IB=",I3," L=",I2," E[NEW]=",F15.5," E[OLD]=",F15.5)') &
     &                  IB,LOFI(IB),EOFI1(IB),EOFI(IB)
        ENDDO
        IF(TWRITE)CALL SETUP_WRITEPHI('UOFI.DAT',GID,NR,NB,UOFI)
      END IF
!
!     ==========================================================================
!     == NORMALIZE EACH ANGULAR MOMENTUM SO THAT FIRST PARTIAL WAVE IS NORMAL ==
!     ==========================================================================
      DO L=0,LX
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          PHI(:)=NLPHI(:,LN)
!         == ORTHOGONALIZE TO CORE STATES  =====================================
          DO IB=NC,1,-1
            IF(LOFI(IB).NE.L) CYCLE
            AUX(:)=R(:)**2*UOFI(:,IB)*PHI(:)
!PB CHANGED RBOX TO ROUT 140121
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR1)
            AUX(:)=R(:)**2*UOFI(:,IB)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!PB CHANGED RBOC TO ROUT 140121
            CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR2)
            VAL=SVAR1/SVAR2
            PHI(:)=PHI(:)-UOFI(:,IB)*VAL
          ENDDO
!         == NORMALIZATION FACTOR  =============================================
          AUX(:)=R(:)**2*PHI(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!PB CHANGED RADIUS FOR NORMALIZATION CONSISTENTLY FROM RBNDOUT TO 3 A.U. 
!PB THIS RADIUS IS PHYSICALLY IRRELEVANT BUT ALLOWS BETTER COMPARISON
          CALL RADIAL$VALUE(GID,NR,AUX1,3.D0,VAL)
          VAL=1.D0/SQRT(VAL)
          CALL RADIAL$VALUE(GID,NR,PHI,3.D0,SVAR)
          VAL=SIGN(VAL,SVAR)
!         == SCALE ALL WAVE FUNCTIONS ==========================================
          DO IB=1,NB
            IF(LOFI(IB).NE.L) CYCLE
            UOFI(:,IB) = UOFI(:,IB)*VAL
            TUOFI(:,IB)=TUOFI(:,IB)*VAL
          ENDDO
!         == SCALE ALL PARTIAL WAVES  ==========================================
          DO LN1=1,LNX
            IF(LOX(LN1).NE.L) CYCLE
            NLPHI(:,LN1)   =NLPHI(:,LN1)*VAL
            TNLPHI(:,LN1)  =TNLPHI(:,LN1)*VAL
          ENDDO
          EXIT ! SCALE ONLY ONCE PER L
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DEFINE PHISCALE AND PSISCALE                                         ==
!     == NEEDED TO COMPENSATE FOR THE HUGE SIZE DIFFERENCE BETWEEN STATES     ==
!     == FROM DIFFERENT SHELLS IN THE NODELESS CONSTRUCTION                   ==
!     ==========================================================================
      PHISCALE(:)=1.D0
      DO LN=1,LNX
        DO LN1=LN+1,LNX
          IF(LOX(LN1).NE.LOX(LN)) CYCLE
          PHISCALE(LN1)=PHISCALE(LN1)*(EOFLN(LN)-EOFLN(LN1))
        ENDDO
      ENDDO
!
      PSISCALE(:)=1.D0
      DO IB=1,NB
        DO IB1=IB+1,NB
          IF(LOFI(IB1).NE.LOFI(IB)) CYCLE
          PSISCALE(IB1)=PSISCALE(IB1)*(EOFI(IB)-EOFI(IB1))
        ENDDO
      ENDDO
!
      DO L=0,LX
        SVAR=1.D0
        DO IB=NC+1,NB
          IF(LOFI(IB).EQ.L) THEN
            SVAR=1.D0/PSISCALE(IB)
            EXIT
          END IF
        ENDDO
        DO IB=1,NB
          IF(LOFI(IB).EQ.L)PSISCALE(IB)=PSISCALE(IB)*SVAR
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == REPORT SETTINGS ON PARTIAL WAVES                                     ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(80("="),T20," ENERGIES FOR PARTIAL-WAVE CONSTRUCTION")')
        WRITE(6,FMT='("RBOX=",F9.5)')RBOX
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," L=",I2," E=",F10.5," RC=",F6.3)') &
     &                      LN,LOX(LN),EOFLN(LN),RC(LN)
        ENDDO
        IF(TWRITE) CALL SETUP_WRITESCALEDPHI(-'NLPHISCALED.DAT',GID,NR,LNX &
     &                                      ,PHISCALE,NLPHI)
        IF(TWRITE) CALL SETUP_WRITESCALEDPHI(-'UOFISCALED.DAT',GID,NR,NB &
     &                                      ,PSISCALE,UOFI)
      END IF
!
!     ==========================================================================
!     == CONSTRUCT QN FUNCTIONS        (H-E)|QN>=|UC>                         ==
!     == THE QN FUNCTIONS MAY HAVE NODES, BUT THE NUMBER OF NODES IS REDUCED  ==
!     == BY THE NUMBER OF CORE STATES                                         ==
!     ==========================================================================
                           CALL TRACE$PASS('CONSTRUCT QN FUNCTIONS')
      TRANSU(:,:)=0.D0
      DO L=0,LX
        IPRO=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          IPRO=IPRO+1
          SVAR=1.D0
          DO LN1=1,LN
            IF(LOX(LN1).NE.L) CYCLE
            IF(TRANSU(LN1,LN).NE.0.D0) THEN
!              === THIS CHECK SHALL BE REMOVED AFTER SOME TESTING PERIOD
               CALL ERROR$MSG('CHECKING ASSUMPTIONS UNDERLYING THE CODE')
               CALL ERROR$MSG('VARIABLE TRANSU IS NOT ZERO AS ASSUMED')
               CALL ERROR$I4VAL('LN',LN)
               CALL ERROR$I4VAL('LN1',LN1)
               CALL ERROR$R8VAL('TRANSU',TRANSU(LN1,LN))
               CALL ERROR$STOP('SETUP_MAKEPARTIALWAVES')
               TRANSU(LN1,LN)=TRANSU(LN1,LN)+SVAR  !OLD STATEMENT
            END IF
            TRANSU(LN1,LN)=SVAR
            SVAR=SVAR*(EOFLN(LN)-EOFLN(LN1))
          ENDDO
          UBYQ(LN)=1.D0/TRANSU(LN,LN) ! MATCHFACTOR |UN> <-> QN*UBYQ
        ENDDO
      ENDDO
      QN=MATMUL(NLPHI,TRANSU)
      TQN=MATMUL(TNLPHI,TRANSU)
!
      CALL LIB$INVERTR8(LNX,TRANSU,TRANSUINV)
!
      IF(TTEST) THEN
        WRITE(6,FMT='(80("="),T20,"  TRANSU ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F15.10)')TRANSU(LN1,:)
        ENDDO
        WRITE(6,FMT='(80("="),T20,"  TRANSU^(-1) ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F15.10)')TRANSUINV(LN1,:)
        ENDDO
        IF(TWRITE)CALL SETUP_WRITEPHI(-'QN.DAT',GID,NR,LNX,QN)
      END IF
!
!     ==========================================================================
!     == ADJUST SCALING OF NODELESS PARTIAL WAVES TO THE QN                   ==
!     == SO THAT LONG-RANGE TALE OF NLPHI AND QPHI ARE IDENTICAL              ==
!     ==========================================================================
      DO LN=1,LNX
        NLPHI(:,LN)=NLPHI(:,LN)/UBYQ(LN)
        TNLPHI(:,LN)=TNLPHI(:,LN)/UBYQ(LN)
      ENDDO
!
!     ==========================================================================
!     == TEST EQUATION FOR QN                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(80("="),T20,"  TEST QN EQ.  ")')
        DO L=0,LX
          IPRO=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            IPRO=IPRO+1
            CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOX(LN),QN(:,LN),AUX)
            PRO(:,LN)=TQN(:,LN)+(AEPOT*Y0-EOFLN(LN))*QN(:,LN)+AUX(:)
            IF(NCL(L).NE.0) PRO(:,LN)=PRO(:,LN)-UOFI(:,NCL(L))
            WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|QN>-|UC>=",F20.15)') &
     &                     LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
          ENDDO
        ENDDO
        IF(TWRITE)CALL SETUP_WRITEPHI(-'TEST.DAT',GID,NR,LNX,PRO)
      END IF
!
!     ==========================================================================
!     == CORE-ORTHOGONALIZE QN TO OBTAIN NODAL AE PARTIAL WAVES               ==
!     == USE LADDER OF NODELESS WAVE FUNCTIONS                                ==
!     == DUE TO NEGLECT OF THE SMALL COMPONENT THE ORTHOGONALIZATION IS NOT   ==
!     == EXACT. HOWEVER THIS CHOICE ENSURES THAT AEPHI FULFILLS THE           ==
!     == SCHRODINGER EQUATION                                                 ==
!     ==========================================================================
      AEPHI(:,:) =QN(:,:)   
      TAEPHI(:,:)=TQN(:,:)  
      DO LN=1,LNX
        SVAR=1.D0
        DO IB=NC,1,-1
          IF(LOFI(IB).NE.LOX(LN)) CYCLE
          SVAR=SVAR/(EOFLN(LN)-EOFI1(IB))
          AEPHI(:,LN) = AEPHI(:,LN)+ UOFI(:,IB)*SVAR
          TAEPHI(:,LN)=TAEPHI(:,LN)+TUOFI(:,IB)*SVAR
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == TEST ORTHOGONALITY OF AEPHI TO CORE STATES                           ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(A(LNX,NC))
        A(:,:)=0.D0
        DO LN=1,LNX
          DO IB=1,NC
            IF(LOX(LN).NE.LOFI(IB)) CYCLE
            AUX(:)=R(:)**2*UOFI(:,IB)*AEPHI(:,LN)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
            A(LN,IB)=VAL
            AUX(:)=R(:)**2*UOFI(:,IB)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
            A(LN,IB)=A(LN,IB)/SQRT(VAL)
          ENDDO
        ENDDO
        WRITE(6,FMT='(80("="),T20," <UC|AEPHI>/SQRT(<UC|UC>  ")')
        WRITE(6,FMT='(80("="),T10," DEVIATION DUE NEGLECT OF SMALL COMPONENT")')
        WRITE(6,FMT='(80("="),T10," ORTHOGONALIZATION DONE BASED ON ENERGIES")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F10.5)')A(LN1,:)
        ENDDO
        DEALLOCATE(A)
      END IF
!
!     ==========================================================================
!     == TEST EQUATION FOR ALL-ELECTRON PARTIAL WAVES AEPHI                   ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(80("="),T20,"  TEST AEPHI EQ.  ")')
        DO LN=1,LNX
          CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOX(LN),AEPHI(:,LN),AUX)
          PRO(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)+AUX(:)
          WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|AEPHI_N> =",F20.15)') &
      &                 LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
        ENDDO
!       CALL SETUP_WRITEPHI(+'AEPHI.DAT',GID,NR,LNX,AEPHI)
      END IF
!
!     ==========================================================================
!     == CONSTRUCT PSEUDO PARTIAL WAVES                                       ==
!     ==========================================================================
                           CALL TRACE$PASS('CONSTRUCT PSEUDO PARTIAL WAVES')
      IF(TMATCHTOALLELECTRON) THEN
        PSPHI=AEPHI
        TPSPHI=TAEPHI
      ELSE
        PSPHI=QN
        TPSPHI=TQN
      END IF
      IF(TTEST.AND.TWRITE)CALL SETUP_WRITEPHI('XX1.DAT',GID,NR,LNX,PSPHI)
      IF(TYPE.EQ.'KERKER') THEN
        DO L=0,LX
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            RC1=RC(LN)
          ENDDO
          NPRO=NPROL(L)
          IF(NPRO.EQ.0) CYCLE
          ALLOCATE(PHITEST(NR,NPRO))
          ALLOCATE(TPHITEST(NR,NPRO))
          IPRO=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            IPRO=IPRO+1
            PHITEST(:,IPRO)=PSPHI(:,LN)
            TPHITEST(:,IPRO)=TPSPHI(:,LN)
          ENDDO
          CALL ATOMIC_MAKEPSPHI(GID,NR,RC1,L,NPRO,PHITEST,TPHITEST)
          IPRO=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            IPRO=IPRO+1
            PSPHI(:,LN)=PHITEST(:,IPRO)
            TPSPHI(:,LN)=TPHITEST(:,IPRO)
          ENDDO
          DEALLOCATE(PHITEST)
          DEALLOCATE(TPHITEST)
        ENDDO
      ELSE IF(TYPE.EQ.'HBS') THEN
        CALL ATOMIC_MAKEPSPHI_HBS(GID,NR,LNX,LOX,EOFLN,RC,LAMBDA,PSPOT &
     &                               ,RBND,PSPHI,TPSPHI)
      ELSE
        CALL ERROR$MSG('PSEUDIZATION TYPE IS UNKNOWN')
        CALL ERROR$MSG('CAN BE "BESSEL" OR "HBS"')
        CALL ERROR$STOP('ATOMIC_MAKEPARTIALWAVES')
      END IF
      IF(TTEST.AND.TWRITE)CALL SETUP_WRITEPHI('XX2.DAT',GID,NR,LNX,PSPHI)
!
!     ==========================================================================
!     == CONSTRUCT BARE PROJECTOR FUNCTIONS                                   ==
!     ==========================================================================
                      CALL TRACE$PASS('CONSTRUCT PROJECTOR FUNCTIONS')
      DO LN=1,LNX
        BAREPRO(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN)
!PRO(:,LN)=TQN(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*QN(:,LN)
      ENDDO
!!$!     ==  CLEANUP OF NUMERICAL ERRORS INDUCED BY A TOO COARSE GRID AT LARGE ====
!!$!     ==  RADII. 2.D0*RCOV IS CHOSEN ARBITRARILY ===============================
!!$      DO IR=1,NR
!!$!        IF(R(IR).LT.RNORM) CYCLE !RNORM IS CHOSEN =RBND WHICH MAY BE TOO SMALL
!!$        IF(R(IR).LT.1.5D0*RCOV) CYCLE
!!$        BAREPRO(IR:,:)=0.D0
!!$        EXIT
!!$      ENDDO
!!$CALL SETUP_WRITEPHI(-'QN.DAT',GID,NR,LNX,QN)
!!$CALL SETUP_WRITEPHI(-'PSPHI.DAT',GID,NR,LNX,PSPHI)
!!$CALL SETUP_WRITEPHI(-'AEPHI.DAT',GID,NR,LNX,AEPHI)
!!$CALL SETUP_WRITEPHI(-'PRO-1.DAT',GID,NR,LNX,BAREPRO)
!!$CALL SETUP_WRITEPHI(-'PRO-1A.DAT',GID,NR,LNX,PRO)
!!$CALL SETUP_WRITEPHI(-'PSPOT.DAT',GID,NR,1,PSPOT)
!!$CALL SETUP_WRITEPHI(-'AEPOT.DAT',GID,NR,1,AEPOT)
!!$CALL SETUP_WRITEPHI(-'DPOT.DAT',GID,NR,1,AEPOT-PSPOT)
!!$STOP 'FORCED'
      IF(TTEST.AND.TWRITE) THEN
        CALL SETUP_WRITEPHI('PRO-BARE.DAT',GID,NR,LNX,BAREPRO)
      END IF
!
!     ==========================================================================
!     == CHECK PAW EQUATION FOR PSEUDO PARTIALWAVES                           ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(PHITEST(NR,LNX))
        ALLOCATE(TPHITEST(NR,LNX))
        DO LN=1,LNX
          TPHITEST(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN) &
     &                  -BAREPRO(:,LN)
          CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOX(LN),AEPHI(:,LN),AUX)
          PHITEST(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)+AUX(:)

        ENDDO
!
        WRITE(6,FMT='(80("="),T20,"  TEST RAW PAW EQUATION  ")')
        STRING='("LN=",I2," L=",I2," RAW PAW EQ.",F10.5'
        STRING=TRIM(ADJUSTL(STRING))//'" SCHR. EQ.",F10.5)'
        DO LN=1,LNX
          WRITE(6,FMT=STRING)LN,LOX(LN),MAXVAL(ABS(TPHITEST(:,LN))) &
    &                                  ,MAXVAL(ABS(PHITEST(:,LN)))
        ENDDO
        DEALLOCATE(PHITEST)
        DEALLOCATE(TPHITEST)
      ENDIF
!
!     ==========================================================================
!     == ENFORCE BIORTHOGONALIZATION                                          ==
!     ==========================================================================
      ALLOCATE(A(LNX,LNX))
      CALL SETUP_BIORTHOMATRICES(GID,NR,RBOX,LNX,LOX,PSPHI,BAREPRO &
     &                          ,TRANSPHI,TRANSPRO)
      A=MATMUL(TRANSPRO,TRANSPOSE(TRANSPHI))
      PRO=MATMUL(BAREPRO,A)  ! ENFORCE BIORTHOGONALITY ON THE PROJECTORS ONLY
      CALL LIB$INVERTR8(LNX,TRANSPHI,TRANSPHIINV)
      DEALLOCATE(A)
!
!     ==========================================================================
!     == CHECK BIORTHOGONALIZATION                                            ==
!     ==========================================================================
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
!         == PROJECTOR FUNCTIONS HAVE A STRICT CUTOFF. INTEGRATION IS DONE =====
!         == TO THE END  =======================================================
          AUX(:)=R(:)**2*PSPHI(:,LN1)*PRO(:,LN2)
!          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
          IF(LN1.EQ.LN2)VAL=VAL-1.D0
          IF(ABS(VAL).GT.1.D-5) THEN
            CALL ERROR$MSG('BIORTHOGONALIZATION FAILED')
            CALL ERROR$I4VAL('L',LOX(LN1))
            CALL ERROR$I4VAL('LN1',LN1)
            CALL ERROR$I4VAL('LN2',LN2)
            CALL ERROR$R8VAL('DEVIATION',VAL)
            CALL ERROR$STOP('ATOMIC_MAKEPARTIALWAVES')
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DTKIN,DOVER                                                          ==
!     == ATTENTION!! DO NOT MAKE DTKIN AND DH SYMMETRIC. THEY ARE NOT HERMITEAN!
!     ==                                                                      ==
!     == THE REASON FOR THE WARNING NOT TO SYMMETRIZE IS NO MORE UNDERSTOOD.  ==
!     == IN ORDER TO OBTAIN A HERMITEAN HAMILTONIAN LATERON, DTKIN MUST BE    ==
!     == SYMMETRIC. THEREFORE IT IS SYMMETRIZED AT THE END OF THE CONSTRUCTION.=
!     ==                                                                      ==
!     == THE ASYMMETRY DECREASES WHEN THE CUTOFF RADIUS FOR THE CONSTRUCTION  ==
!     == OF PSEUDO PARTIAL WAVES INCREASED. IT DECREASES WHE THE RADIUS       ==
!     == RBOX/RCOV AS DEFINED IN THE SETUP-PARAMETER FILE IS INCREASED        ==
!     ==                                                                      ==
!     == BOTH, AEPHI AND PSPHI ARE CONSTRUCTED FROM QN                        ==
!     == IT SEEMS THAT REPLACING AEPHI BY QN REDUCES MOST OF THE ASYMMETRY    ==
!     ==                                                                      ==
!     == THE ASYMMETRY ALSO OCCURS IN NON-RELATIVISTIC CALCULATIONS           ==
!     ==========================================================================
      DTKIN=0.D0
      DOVER=0.D0
      DH=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
! IF EVERYTHING IS CORRECT, THE INTEGRATION SHOULD BE DONE TOWARDS THE
! END OF THE GRID. HOWEVER THE EXPONENTIALLY GROWING TAIL OF THE PARTIAL 
! WAVES MAY INTRODUCE NUMERICAL ERRORS
          AUX(:)=R(:)**2*(AEPHI(:,LN1)*TAEPHI(:,LN2)-PSPHI(:,LN1)*TPSPHI(:,LN2))
!AUX(:)=R(:)**2*(AEPHI(:,LN1)*TAEPHI(:,LN2)-PSPHI(:,LN1)*TPSPHI(:,LN2))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBNDOUT,VAL)
!CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL)
          DTKIN(LN1,LN2)=VAL
          AUX(:)=R(:)**2*(AEPHI(:,LN1)*AEPHI(:,LN2)-PSPHI(:,LN1)*PSPHI(:,LN2))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBNDOUT,VAL)
!CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL)
          DOVER(LN1,LN2)=VAL
          CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOX(LN2),AEPHI(:,LN2),AUX1)
          AUX(:)=R(:)**2*(AEPHI(:,LN1)*(AEPOT(:)*Y0*AEPHI(:,LN2)+AUX1(:)) &
      &                  -PSPHI(:,LN1)*PSPOT(:)*Y0*PSPHI(:,LN2))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBNDOUT,VAL)
!CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,VAL)
          DH(LN1,LN2)=DTKIN(LN1,LN2)+VAL
        ENDDO
      ENDDO
PRINT*,'VFOCK%TON  ',VFOCK%TON
PRINT*,'RBND       ',RBND
DO LN1=1,LNX
  WRITE(*,FMT='(A,100E12.5)')'DT   ',DTKIN(LN1,:)
ENDDO
DO LN1=1,LNX
  WRITE(*,FMT='(A,100E12.5)')'DH   ',DH(LN1,:)
ENDDO
DO LN1=1,LNX
  WRITE(*,FMT='(A,100E12.5)')'DO   ',DOVER(LN1,:)
ENDDO
      IF(TTEST) THEN
!       == THE NON-HERMIEANITY COMES (FOR SILICON) TO 85 PERCENT FROM         ==
!       == THE ADMIXTURE OF THE CORE WAVE FUNCTIONS, WHEN AEPHI IS OBTAINED.  ==
!       == ABOUT 15 PERCENT CAN BE ATTRIBUTED TO THE CONSTRUCTION OF THE      ==
!       == PSEUDO PARTIAL WAVES FROM THE QN                                   ==
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," DTKIN-TRANSPOSE(DTKIN)=",10F10.5)') &
     &                                LN,(DTKIN(LN,LN2)-DTKIN(LN2,LN),LN2=1,LNX)
        ENDDO
      END IF
!
!     ==========================================================================
!     == CHECK PAW EQUATION FOR PSEUDO PARTIALWAVES                           ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(TPHITEST(NR,LNX))    ! HOLDS TEST FOR PSEUDO
        ALLOCATE(PHITEST(NR,LNX))     ! HOLDS TEST FOR ALL-ELECTRON 
        ALLOCATE(PROJ(LNX))
        DO LN=1,LNX
!         == DETERMINE PROJECTIONS
          PROJ(:)=0.D0
          DO LN1=1,LNX
            IF(LOX(LN1).NE.LOX(LN)) CYCLE
            AUX(:)=R(:)**2*PRO(:,LN1)*PSPHI(:,LN)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBNDOUT,VAL)
            PROJ(LN1)=VAL
          ENDDO
          WRITE(6,FMT='("LN=",I2," <P|PSPHI>=",10F10.5)')LN,PROJ
!
          CALL RADIALFOCK$VPSI(GID,NR,VFOCK,LOX(LN),AEPHI(:,LN),AUX)
          PHITEST(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)+AUX(:)
          TPHITEST(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN)
!TPHITEST(:,LN)=TPHITEST(:,LN)-BAREPRO(:,LN)
          PHITEST1(:,LN)=0.D0          
          DO LN1=1,LNX
            IF(LOX(LN1).NE.LOX(LN)) CYCLE
            SVAR=0.D0
            DO LN2=1,LNX
              IF(LOX(LN2).NE.LOX(LN)) CYCLE
              SVAR=SVAR+(DH(LN1,LN2)-EOFLN(LN)*DOVER(LN1,LN2))*PROJ(LN2)
            ENDDO
            TPHITEST(:,LN)=TPHITEST(:,LN)+PRO(:,LN1)*SVAR
            PHITEST1(:,LN)=PHITEST1(:,LN)+PRO(:,LN1)*SVAR
          ENDDO
        ENDDO
        WRITE(6,FMT='(80("="),T20,"  TEST PAW EQUATION  ")')
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," L=",I2," PAW EQ.",F20.5' &
     &                                //'," SCHR. EQ.",F20.5," DPRO ",F20.5)') &
     &          LN,LOX(LN),MAXVAL(ABS(TPHITEST(:,LN))) &  
     &                    ,MAXVAL(ABS(PHITEST(:,LN))) &   ! ALL-ELECTRON EQ
     &                    ,MAXVAL(ABS(PHITEST1(:,LN)+BAREPRO(:,LN)))
        ENDDO
        DEALLOCATE(PROJ)
        DEALLOCATE(PHITEST)
        DEALLOCATE(TPHITEST)
      END IF
!
!     ==========================================================================
!     == CONSTRUCT PHIDOT FUNCTIONS                                           ==
!     ==========================================================================
!     == TRANSFORMATION TO NODELESS REPRESENTATION =============================
      QNP(:,:)   =MATMUL(QN,TRANSUINV)
      PSPHIP(:,:)=MATMUL(PSPHI,TRANSUINV)
      DO L=0,LX
        NPRO=NPROL(L)
        IF(NPRO.EQ.0) CYCLE
        ALLOCATE(DH1(NPRO,NPRO))
        ALLOCATE(DO1(NPRO,NPRO))
        ALLOCATE(PRO1(NR,NPRO))
        ALLOCATE(PROJ(NPRO))
        IPRO1=0
        DO LN1=1,LNX
          IF(LOX(LN1).NE.L) CYCLE
          IPRO1=IPRO1+1
          PRO1(:,IPRO1)=PRO(:,LN1)
          IPRO2=0
          DO LN2=1,LNX
            IF(LOX(LN2).NE.L) CYCLE
            IPRO2=IPRO2+1
            DH1(IPRO1,IPRO2)=DH(LN1,LN2)
            DO1(IPRO1,IPRO2)=DOVER(LN1,LN2)
          ENDDO
        ENDDO
        G(:)=0.D0
        LNLAST=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          ESCATT(LN)=MIN(-0.0D0,EOFLN(LN))
          E=ESCATT(LN)
          IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
!
!         == CALCULATE NLPHIDOT    =============================================
          G(:)=NLPHI(:,LN)
          CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,0,G,L,E,1 &
    &                                                           ,NLPHIDOT(:,LN))
          CALL ATOMLIB$UPDATESTATEWITHHF(GID,NR,L,0,DREL,G,AEPOT,VFOCK &
    &                                    ,-1.D0,E,NLPHIDOT(:,LN))
!THIS CONSTRUCTION HAS BEEN REPLACED, BECAUSE THE DOT FUNCTIONS DID NOT HAVE 
! THE SAME TAIL BEHAVIOR.
!!$!
!!$!         == CALCULATE QNDOT  ==================================================
!!$          G(:)=QN(:,LN)
!!$          CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,0,G,L,E,1,QNDOT(:,LN))
!!$          CALL ATOMLIB$UPDATESTATEWITHHF(GID,NR,L,0,DREL,G,AEPOT,VFOCK &
!!$    &                                    ,-1.D0,E,QNDOT(:,LN))
!!$!
!!$!         == CALCULATE PSEUDO WAVE FUNCTIONS ===================================
!!$          G(:)=PSPHI(:,LN)
!!$          DO IPRO=1,NPRO
!!$            AUX(:)=R(:)**2*PRO(:,IPRO)*PSPHI(:,LN)
!!$!           == INTEGRAL IS OK BECAUSE PRO IS EXACTLY ZERO BEYOND A RADIUS ======
!!$            CALL RADIAL$INTEGRAL(GID,NR,AUX,PROJ(IPRO))
!!$          ENDDO
!!$          PROJ(:)=MATMUL(DO1,PROJ)
!!$          DO IPRO=1,NPRO
!!$            G(:)=G(:)+PRO(:,IPRO)*PROJ(IPRO)
!!$          ENDDO
!!$          CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO1,DH1,DO1,G &
!!$    &                        ,PSPHIDOT(:,LN))
!!$!
!!$!         == ADD HOMOGENEOUS SOLUTION TO MATCH OUTER BOUNDARY CONDITIONS =======
!!$          G(:)=0.D0
!!$          CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO1,DH1,DO1,G,PHI)
!!$          CALL RADIAL$VALUE(GID,NR,PHI,MAXVAL(RC),VAL)
!!$          PHI=PHI/VAL
!!$          CALL RADIAL$VALUE(GID,NR,NLPHIDOT(:,LN)-PSPHIDOT(:,LN),MAXVAL(RC),VAL)
!!$          PSPHIDOT(:,LN)=PSPHIDOT(:,LN)+PHI(:)*VAL
!!$!         == REPLACE TAILS TO AVOID NUMERICAL ERRORS ===========================
!!$          SVAR=MAXVAL(RC)
!!$          DO IR=1,NR
!!$            IF(R(IR).LE.SVAR) CYCLE
!!$            PSPHIDOT(IR:,LN)=NLPHIDOT(IR:,LN)
!!$            EXIT
!!$          ENDDO
!
!         ======================================================================
!         == CONSTRUCT 
!         ======================================================================
          QNDOT(:,LN)=NLPHIDOT(:,LN)
          PSPHIDOT(:,LN)=NLPHIDOT(:,LN)
!
!         ======================================================================
!         == CONSTRUCT AEPHIDOT BY CORE-ORTHOGONALIZATION                     ==
!         == AEPHIDOT DOES NOT OBEY (H-E)|AEPHIDOT>=|AEPHI> !!                ==
!         ======================================================================
          AEPHIDOT(:,LN)=NLPHIDOT(:,LN)
          DO IB=NC,1,-1
            IF(LOFI(IB).NE.L) CYCLE
            AUX=R(:)**2*UOFI(:,IB)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBNDOUT,SVAR1)
            AUX=R(:)**2*UOFI(:,IB)*AEPHIDOT(:,LN)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBNDOUT,SVAR2)
            AEPHIDOT(:,LN)=AEPHIDOT(:,LN)-UOFI(:,IB)*SVAR2/SVAR1
          ENDDO
          LNLAST=LN
        ENDDO
        DEALLOCATE(DH1)
        DEALLOCATE(DO1)
        DEALLOCATE(PRO1)
        DEALLOCATE(PROJ)
      ENDDO
!     == QNDOT DOES NOT HAVE THE SAME RADIAL LONG-RANGE BEHAVIOR AS THE OTHER
!     == DOT-FUNCTIONS!
!!$CALL SETUP_WRITEPHI('NL.DAT',GID,NR,LNX,NLPHI)
!!$CALL SETUP_WRITEPHI('QN.DAT',GID,NR,LNX,QN)
!!$CALL SETUP_WRITEPHI('PS.DAT',GID,NR,LNX,PSPHI)
!!$CALL SETUP_WRITEPHI('AE.DAT',GID,NR,LNX,AEPHI)
!!$CALL SETUP_WRITEPHI(+'AEPHIDOT.DAT',GID,NR,LNX,AEPHIDOT)
!!$CALL SETUP_WRITEPHI(+'PSPHIDOT.DAT',GID,NR,LNX,PSPHIDOT)
!!$CALL SETUP_WRITEPHI(+'NLPHIDOT.DAT',GID,NR,LNX,NLPHIDOT)
!!$CALL SETUP_WRITEPHI(+'QNDOT.DAT',GID,NR,LNX,QNDOT)
!!$STOP 'OK HERE'
!
!     ==========================================================================
!     == BACK TRANSFORM                                                       ==
!     ==========================================================================
GOTO 100
      QN=MATMUL(QN,TRANSUINV)
      TQN=MATMUL(TQN,TRANSUINV)
      PSPHI=MATMUL(PSPHI,TRANSUINV)
      TPSPHI=MATMUL(TPSPHI,TRANSUINV)
      AEPHI=MATMUL(AEPHI,TRANSUINV)
      TAEPHI=MATMUL(TAEPHI,TRANSUINV)
      PRO=MATMUL(PRO,TRANSPOSE(TRANSU))
      DTKIN=MATMUL(TRANSPOSE(TRANSUINV),MATMUL(DTKIN,TRANSUINV))
      DOVER=MATMUL(TRANSPOSE(TRANSUINV),MATMUL(DOVER,TRANSUINV))
      DH=MATMUL(TRANSPOSE(TRANSUINV),MATMUL(DH,TRANSUINV))
!
!     == TEST IF BACK TRANSFORM WAS SUCCESSFUL ================================
      IF(TTEST) THEN
        WRITE(6,FMT='(80("="),T20,"  TEST BACK TRANSFORM  ")')
        STRING='("LN=",I2," L=",I2," DIFF. NDLSS PHI",F10.5'
        STRING=TRIM(ADJUSTL(STRING))//'" DIFF. KIN.OP NDLSS. PHI ",F10.5)'
        DO LN=1,LNX
          WRITE(6,FMT=STRING)LN,LOX(LN),MAXVAL(ABS(QN(:,LN)-NLPHI(:,LN))) &
     &                                 ,MAXVAL(ABS(TQN(:,LN)-TNLPHI(:,LN)))
        ENDDO
      END IF
100 CONTINUE
!
!     ==========================================================================
!     == RENORMALIZE WAVE FUNCTIONS AND PROJECTOR FUNCTIONS                   ==
!     ==========================================================================
GOTO 10001
      DO L=0,LX
        IPRO=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          IPRO=IPRO+1
!         == NORMALIZE PS PARTIAL WAVE =========================================
!          IF(IPRO.EQ.1) THEN
          IF(ISCATT(LN).LE.0) THEN    ! NORMALIZE VALENCE AND SEMI-CORE STATES
            AUX(:)=R(:)**2*PSPHI(:,LN)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBNDOUT,VAL)
            VAL=VAL+DOVER(LN,LN)
            VAL=1.D0/SQRT(VAL)
          END IF
          DO LN2=LN,LNX
            IF(LOX(LN2).NE.L) CYCLE
            PSPHI(:,LN2)   =   PSPHI(:,LN2)*VAL
            TPSPHI(:,LN2)  =  TPSPHI(:,LN2)*VAL
            PSPHIDOT(:,LN2)=PSPHIDOT(:,LN2)*VAL
            AEPHI(:,LN2)   =   AEPHI(:,LN2)*VAL
            TAEPHI(:,LN2)  =  TAEPHI(:,LN2)*VAL
            AEPHIDOT(:,LN2)=AEPHIDOT(:,LN2)*VAL
            NLPHI(:,LN2)   =   NLPHI(:,LN2)*VAL
            TNLPHI(:,LN2)  =  TNLPHI(:,LN2)*VAL
            NLPHIDOT(:,LN2)=NLPHIDOT(:,LN2)*VAL
            QN(:,LN2)      =      QN(:,LN2)*VAL
            TQN(:,LN2)     =     TQN(:,LN2)*VAL
            QNDOT(:,LN2)   =   QNDOT(:,LN2)*VAL
            PRO(:,LN2)     =     PRO(:,LN2)/VAL
            DH(LN2,:)      =      DH(LN2,:)*VAL
            DTKIN(LN2,:)   =      DTKIN(LN2,:)*VAL
            DOVER(LN2,:)   =   DOVER(LN2,:)*VAL
            DH(:,LN2)      =      DH(:,LN2)*VAL
            DTKIN(:,LN2)   =      DTKIN(:,LN2)*VAL
            DOVER(:,LN2)   =   DOVER(:,LN2)*VAL
          ENDDO
        ENDDO
      ENDDO
10001 CONTINUE
!
!     ==========================================================================
!     == CUT OFF THE EXPONENTIALLY GROWING TAIL OF THE PARTIALWAVES
!     ==========================================================================
      IF(TCUTTAIL) THEN
        DO LN=1,LNX
          SVAR=0.D0
          DO IR=1,NR
            IF(R(IR).LE.RCOV) THEN
               SVAR=MAX(SVAR,ABS(PSPHI(IR,LN)))
            ELSE
              IF(ABS(PSPHI(IR,LN)).GT.10.D0*SVAR) THEN
PRINT*,'CUT PARTIAL WAVE TAIL FOR LN=',LN,' AT R=',R(IR),' AND BEYOND'
                AEPHI(IR:,LN)=0.D0
                PSPHI(IR:,LN)=0.D0
                NLPHI(IR:,LN)=0.D0
                QN(IR:,LN)=0.D0
                AEPHIDOT(IR:,LN)=0.D0
                NLPHIDOT(IR:,LN)=0.D0
                QNDOT(IR:,LN)=0.D0
                PSPHIDOT(IR:,LN)=0.D0
                PRO(IR:,LN)=0.D0
                EXIT
              END IF
            END IF
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == TRANSFORM ONTO SEQUENTIAL RESPRESENTATION                            ==
!     ==========================================================================
!     == CHANGES INDIVIDUAL ENERGY CONTRIBUTIONS LESS THAN 1.E-7 HARTREE
!     == THE SEQUENTIAL REPRESENTATION ALLOWS TO TRUNCATE THE SET OF 
!     ==  PARTIAL WAVES AND PROJECTOR FUNCTIONS AT EVERY NUMBER OF PROJECTORS 
!     == WHICH WAS HISTORICALLY RELEVANT WHEN THE SETUPS WERE STORED ON FILE.
!     == NOW IT MAY BE RELEVANT FOR THE MAPPING ON LOCAL ORBITALS. PLEASE CHECK!
!     ==========================================================================
      IF(TSEQUENTIALAUGMENT) THEN
!!$CALL SETUP_WRITEPHI('NLPHI_1.DAT',GID,NR,LNX,NLPHI)
!!$CALL SETUP_WRITEPHI('AEPHI_1.DAT',GID,NR,LNX,AEPHI)
!!$CALL SETUP_WRITEPHI('PSPHI_1.DAT',GID,NR,LNX,PSPHI)
        PSPHI=MATMUL(PSPHI,TRANSPHI)
        AEPHI=MATMUL(AEPHI,TRANSPHI)
        NLPHI=MATMUL(NLPHI,TRANSPHI)
        QN=MATMUL(QN,TRANSPHI)
        PRO=MATMUL(PRO,TRANSPOSE(TRANSPHIINV))
        DTKIN=MATMUL(TRANSPOSE(TRANSPHI),MATMUL(DTKIN,TRANSPHI))
        DOVER=MATMUL(TRANSPOSE(TRANSPHI),MATMUL(DOVER,TRANSPHI))
!!$CALL SETUP_WRITEPHI('NLPHI_2.DAT',GID,NR,LNX,NLPHI)
!!$CALL SETUP_WRITEPHI('AEPHI_2.DAT',GID,NR,LNX,AEPHI)
!!$CALL SETUP_WRITEPHI('PSPHI_2.DAT',GID,NR,LNX,PSPHI)
!!$CALL ERROR$STOP('FORCED IN PAW_SETUP')
      END IF
                                     CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_OUTERNEWPROWRAPPER(GID,NR,ROUT,RELTYPE &
     &                  ,NB,NC,LOFI,SOFI,EOFI,LNX,LOX,RC,AEPOT,PSPOT,VFOCK &
     &                  ,G_QNPHI,G_QNPHISM,G_AEPHI,G_AEPHISM,G_PSPHI,G_PSPHISM &
     &                  ,G_PRO,G_DTKIN,G_DOVER &
     &                  ,G_AEPHIDOT,G_AEPHIDOTSM,G_QNPHIDOT,G_QNPHIDOTSM)
!     **************************************************************************
!     **  CONSTRUCT ALL-ELECTRON AND PSEUDO WAVE FUNCTIONS USING THE NODE-LESS**
!     **  CONSTRUCTION.                                                       **
!     **                                                                      **
!     **  ONLY THE RADIUS RC FOR THE FIRST LN PER L IS USED TO DEFINE THE     **
!     **  PSEUDO PARTIAL WAVES.                                               **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2013******************
      USE RADIALFOCK_MODULE, ONLY: VFOCK_TYPE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: GID    ! GRID ID
      INTEGER(4),INTENT(IN)  :: NR
      REAL(8)   ,INTENT(IN)  :: ROUT   ! RADIUS OF ENCLOSING SPHERE
      CHARACTER(*),INTENT(IN):: RELTYPE
      INTEGER(4),INTENT(IN)  :: NB
      INTEGER(4),INTENT(IN)  :: NC
      INTEGER(4),INTENT(IN)  :: LOFI(NB)
      INTEGER(4),INTENT(IN)  :: SOFI(NB)
      REAL(8)   ,INTENT(IN)  :: EOFI(NB)
      INTEGER(4),INTENT(IN)  :: LNX
      INTEGER(4),INTENT(IN)  :: LOX(LNX)
      REAL(8)   ,INTENT(IN)  :: RC(LNX)  !ONLY THE FIRST PER L IS USED!!
      REAL(8)   ,INTENT(IN)  :: AEPOT(NR)
      REAL(8)   ,INTENT(IN)  :: PSPOT(NR)
      TYPE(VFOCK_TYPE),INTENT(IN) :: VFOCK

      REAL(8)   ,INTENT(OUT) :: G_QNPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT) :: G_AEPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT) :: G_PSPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT) :: G_PRO(NR,LNX)
      REAL(8)   ,INTENT(OUT) :: G_DTKIN(LNX,LNX)
      REAL(8)   ,INTENT(OUT) :: G_DOVER(LNX,LNX)
      REAL(8)   ,INTENT(OUT) :: G_QNPHISM(NR,LNX)
      REAL(8)   ,INTENT(OUT) :: G_AEPHISM(NR,LNX)
      REAL(8)   ,INTENT(OUT) :: G_PSPHISM(NR,LNX)
      REAL(8)   ,INTENT(OUT) :: G_QNPHIDOT(NR,LNX)
      REAL(8)   ,INTENT(OUT) :: G_QNPHIDOTSM(NR,LNX)
      REAL(8)   ,INTENT(OUT) :: G_AEPHIDOT(NR,LNX)
      REAL(8)   ,INTENT(OUT) :: G_AEPHIDOTSM(NR,LNX)

      REAL(8)   ,ALLOCATABLE :: EC(:)
      REAL(8)   ,ALLOCATABLE :: UCORE(:,:)
      REAL(8)   ,ALLOCATABLE :: UCORESM(:,:)
      REAL(8)   ,ALLOCATABLE :: AECORE(:,:)
      REAL(8)   ,ALLOCATABLE :: AECORESM(:,:)
      REAL(8)   ,ALLOCATABLE :: PSCORE(:,:)
      REAL(8)   ,ALLOCATABLE :: PSCORESM(:,:)
      REAL(8)   ,ALLOCATABLE :: QN(:,:)
      REAL(8)   ,ALLOCATABLE :: QNSM(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHISM(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHISM(:,:)
      REAL(8)   ,ALLOCATABLE :: QNDOT(:)
      REAL(8)   ,ALLOCATABLE :: QNDOTSM(:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOTSM(:)
      REAL(8)   ,ALLOCATABLE :: PRO(:,:)
      REAL(8)   ,ALLOCATABLE :: DTKIN(:,:)
      REAL(8)   ,ALLOCATABLE :: DOVER(:,:)
      REAL(8)   ,ALLOCATABLE :: EOFPHI(:)
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      LOGICAL(4),PARAMETER   :: TWRITE=.TRUE.
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER   :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4)             :: NBL
      INTEGER(4)             :: NCL  
      INTEGER(4)             :: NPHI !#(PARTIAL WAVES FOR THIS L)
      INTEGER(4)             :: LX
      INTEGER(4)             :: L,SO,IB,IBL,IPHI,JPHI,LN,LN1,LN2
      REAL(8)                :: RCL  !PSEUDIZATION RADIUS FOR THIS L
      REAL(8)                :: ENU  !EXPANSION ENERGY FOR TAYLOR EXPANSION
!     **************************************************************************
                                     CALL TRACE$PUSH('SETUP_OUTERNEWPROWRAPPER')
      G_QNPHI=0.D0
      G_AEPHI=0.D0
      G_PSPHI=0.D0
      G_QNPHISM=0.D0
      G_AEPHISM=0.D0
      G_PSPHISM=0.D0
      G_QNPHIDOT=0.D0
      G_AEPHIDOT=0.D0
      G_QNPHIDOTSM=0.D0
      G_AEPHIDOTSM=0.D0
      G_PRO=0.D0
      G_DTKIN=0.D0
      G_DOVER=0.D0

      IF(RELTYPE.NE.'NONREL'.AND.RELTYPE.NE.'ZORA' &
                            .AND.RELTYPE.NE.'SPINORBIT') THEN
        CALL ERROR$MSG('ILLEGAL VALUE OF IDENTIFIER RELTYPE')
        CALL ERROR$CHVAL('RELTYPE',RELTYPE)
        CALL ERROR$STOP('SETUP_OUTERNEWPROWRAPPER')
      END IF
!
      LX=MAX(MAXVAL(LOFI),MAXVAL(LOX))
      DO L=0,LX
        DO SO=MINVAL(SOFI),1,2  
          IF(L.EQ.0.AND.SO.EQ.-1) CYCLE
!
WRITE(*,FMT='(80("-"),T20," L=",I2," SO=",I2)')L,SO
!
!         ======================================================================
!         == COUNT NUMBER OF BANDS AND NUMBER OF CORE STATES                  ==
!         ======================================================================
          NBL=0
          NCL=0
          DO IB=1,NB
            IF(LOFI(IB).NE.L) CYCLE
            IF(SOFI(IB).NE.SO) CYCLE
            NBL=NBL+1
            IF(IB.LE.NC)NCL=NCL+1
          ENDDO
!
!         ======================================================================
!         == COUNT NUMBER OF PARTIAL WAVES
!         ======================================================================
          NPHI=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            NPHI=NPHI+1
          ENDDO
!
!         == COLLECT PSEUDIZATION RADIUS =======================================
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            RCL=RC(LN)
            EXIT
          ENDDO
!
!         ======================================================================
!         == COLLECT CORE ENERGIES AND ENU                                    ==
!         ======================================================================
          ALLOCATE(EC(NCL))
          IBL=0
          DO IB=1,NC
            IF(LOFI(IB).NE.L) CYCLE
            IF(SOFI(IB).NE.SO) CYCLE
            IBL=IBL+1
            EC(IBL)=EOFI(IB)
          ENDDO
!
!         == COLLECT ENU =======================================================
          ALLOCATE(EOFPHI(NPHI))
          EOFPHI=MAXVAL(AEPOT*Y0)-1.D-2  ! ONLY BOUND STATES ARE ALLOWED
          IBL=0
          DO IB=NC+1,NB
            IF(LOFI(IB).NE.L) CYCLE
            IF(SOFI(IB).NE.SO) CYCLE
            IBL=IBL+1
            EOFPHI(IBL:)=EOFI(IB)
          ENDDO
!
!         ======================================================================
!         == CALCULATE PARTIAL WAVES ETC
!         ======================================================================
          ALLOCATE(UCORE(NR,NCL))
          ALLOCATE(UCORESM(NR,NCL))
          ALLOCATE(AECORE(NR,NCL))
          ALLOCATE(AECORESM(NR,NCL))
          ALLOCATE(PSCORE(NR,NCL))
          ALLOCATE(PSCORESM(NR,NCL))
          ALLOCATE(QN(NR,NPHI))
          ALLOCATE(QNSM(NR,NPHI))
          ALLOCATE(AEPHI(NR,NPHI))
          ALLOCATE(AEPHISM(NR,NPHI))
          ALLOCATE(PSPHI(NR,NPHI))
          ALLOCATE(PSPHISM(NR,NPHI))
          ALLOCATE(QNDOT(NR))
          ALLOCATE(QNDOTSM(NR))
          ALLOCATE(AEPHIDOT(NR))
          ALLOCATE(AEPHIDOTSM(NR))
          ALLOCATE(PRO(NR,NPHI))
          ALLOCATE(DTKIN(NPHI,NPHI))
          ALLOCATE(DOVER(NPHI,NPHI))
IF(TWRITE) THEN
PRINT*,'------------- INPUT DATA FOR SETUP_NEWPRO-----------'
PRINT*,'L        ',L
PRINT*,'SO        ',SO
PRINT*,'ROUT     ',ROUT
PRINT*,'RCL      ',RCL
PRINT*,'NCL,NPHI ',NCL,NPHI
PRINT*,'EC       ',EC
PRINT*,'EOFPHI   ',EOFPHI
PRINT*,'-----------------------------------------------------'
END IF
          CALL SETUP_NEWPRO(RELTYPE,GID,NR,ROUT,L,SO,NCL,NPHI,RCL,EOFPHI,EC &
     &                  ,AEPOT,PSPOT,VFOCK &
     &                  ,UCORE,AECORE,PSCORE,QN,AEPHI,PSPHI &
     &                  ,QNDOT,QNDOTSM,AEPHIDOT,AEPHIDOTSM &
     &                  ,UCORESM,AECORESM,PSCORESM,QNSM,AEPHISM,PSPHISM &
     &                  ,PRO,DTKIN,DOVER)
!
!         ======================================================================
!         == TESTS                                                            ==
!         ======================================================================
          IF(TTEST) THEN
            CALL SETUP_TEST_NEWPRO(RELTYPE,GID,NR,ROUT,L,SO,NCL,NPHI,RCL &
     &                  ,EOFPHI,EC &
     &                  ,AEPOT,PSPOT &
     &                  ,UCORE,AECORE,PSCORE,QN,AEPHI,PSPHI,QNDOT &
     &                  ,UCORESM,AECORESM,PSCORESM,QNSM,AEPHISM,PSPHISM &
     &                  ,PRO,DTKIN,DOVER)
          END IF
!IF(L.EQ.3) STOP 'FORCED'
!
!         ======================================================================
!         == SORT INTO THE ARRAY                                              ==
!         ======================================================================
          IPHI=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
IF(SO.NE.0) THEN
  CALL ERROR$MSG('SPIN ORBIT IMPLEMENTATION NOT COMPLETE')
  CALL ERROR$STOP('SETUP_OUTERNEWPROWRAPPER')
END IF         
            IPHI=IPHI+1
            G_AEPHI(:,LN)=AEPHI(:,IPHI)
            G_PSPHI(:,LN)=PSPHI(:,IPHI)
            G_QNPHI(:,LN)=QN(:,IPHI)
            G_PRO(:,LN) =PRO(:,IPHI)
            G_AEPHISM(:,LN)=AEPHISM(:,IPHI)
            G_PSPHISM(:,LN)=PSPHISM(:,IPHI)
            G_QNPHISM(:,LN)=QNSM(:,IPHI)
            IF(IPHI.LT.NPHI) THEN
              G_QNPHIDOT(:,LN)  =QN(:,IPHI+1)
              G_QNPHIDOTSM(:,LN)=QNSM(:,IPHI+1)
              G_AEPHIDOT(:,LN)  =AEPHI(:,IPHI+1)
              G_AEPHIDOTSM(:,LN)=AEPHISM(:,IPHI+1)
            ELSE
              G_QNPHIDOT(:,LN)  =QNDOT(:)
              G_QNPHIDOTSM(:,LN)=QNDOTSM(:)
              G_AEPHIDOT(:,LN)  =AEPHIDOT(:)
              G_AEPHIDOTSM(:,LN)=AEPHIDOTSM(:)
            END IF
          ENDDO            
!
!         == MATRIX ELEMENTS ===================================================
          IPHI=0
          DO LN1=1,LNX
            IF(LOX(LN1).NE.L) CYCLE
            IPHI=IPHI+1
            JPHI=0
            DO LN2=1,LNX
              IF(LOX(LN2).NE.L) CYCLE
              JPHI=JPHI+1
              G_DTKIN(LN1,LN2)=DTKIN(IPHI,JPHI)
              G_DOVER(LN1,LN2)=DOVER(IPHI,JPHI)
            ENDDO
          ENDDO
!
!         ======================================================================
!         == SORT INTO THE ARRAY                                              ==
!         ======================================================================
          DEALLOCATE(EC)
          DEALLOCATE(EOFPHI)
          DEALLOCATE(UCORE)
          DEALLOCATE(UCORESM)
          DEALLOCATE(AECORE)
          DEALLOCATE(AECORESM)
          DEALLOCATE(PSCORE)
          DEALLOCATE(PSCORESM)
          DEALLOCATE(QN)
          DEALLOCATE(QNSM)
          DEALLOCATE(AEPHI)
          DEALLOCATE(AEPHISM)
          DEALLOCATE(PSPHI)
          DEALLOCATE(PSPHISM)
          DEALLOCATE(QNDOT)
          DEALLOCATE(QNDOTSM)
          DEALLOCATE(AEPHIDOT)
          DEALLOCATE(AEPHIDOTSM)
          DEALLOCATE(PRO)
          DEALLOCATE(DTKIN)
          DEALLOCATE(DOVER)
        ENDDO
      ENDDO
                                     CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_NEWPRO(RELTYPE,GID,NR,ROUT,L,SO,NC,NJ,RC,EOFPHI &
     &                        ,ECOREIN,AEPOT,PSPOT,VFOCK &
     &                        ,UCORE,AECORE,PSCORE,QN,AEPHI,PSPHI &
     &                        ,QNDOT,QNDOTSM,AEPHIDOT,AEPHIDOTSM &
     &                        ,UCORESM,AECORESM,PSCORESM,QNSM,AEPHISM,PSPHISM &
     &                        ,PRO,DTKIN,DOVER)
!     **************************************************************************
!     ** CALCULATES PARTIAL WAVES AND PROJECTOR FUNCTIONS FOR A GIVEN L AND SO**
!     **                                                                      **
!     **  THE CORE STATES SOLVE                                               **
!     **     (H-E(I))|U_I>=|U_{I-1}>   WITH |U_0>=|0>                         **
!     **  THE NODE-REDUCED VALENCE STATES OBEY                                **
!     **     (H-ENU)|QN_J>=|QN_{J-1}>   WITH |QN_0>=|U_{NC}>                  **
!     **                                                                      **
!     **  RELTYPE='NONREL', 'ZORA', 'SCALAR', 'SPINORBIT'                     **
!     **                                                                      **
!     **  FOR THE SMALL COMPONENT, WE USE THE MODEL THAT ALSO THE PSEUDO-     **
!     **  PARTIAL WAVES CARRY A SMALL COMPONENT, WHICH IS IDENTICAL TO THAT   **
!     **  OF THE NODE-REDUCED PARTIAL WAVES QN. DTKIN AND DOVER ARE DONE      **
!     **  ON THIS BASIS, WHICH ENSURES THAT THERE ARE NO EXPONENTIALLY GROWING**
!     **  CONTRIBUTIONS IN THE INTEGRALS.                                     **
!     **                                                                      **
!     **  THE QNDOT, QNDOTSM FUNCTIONS ARE, PER DEFINITION, IDENTICAL         **
!     **  TO THE PSDOT,PSDOTSM                                                **
!     **                                                                      **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2013******************
      USE STRINGS_MODULE
      USE RADIALFOCK_MODULE, ONLY: VFOCK_TYPE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: RELTYPE ! SELECTOR FOR RELATIVISTIC TREATMENT
      INTEGER(4),INTENT(IN) :: GID       ! GRID ID
      INTEGER(4),INTENT(IN) :: NR        ! #(GRID POINTS)
      REAL(8)   ,INTENT(IN) :: ROUT      ! RADIUS OF ATOM IN A BOX
      INTEGER(4),INTENT(IN) :: L         ! MAIN ANGULAR MOMENTUM
      INTEGER(4),INTENT(IN) :: SO        ! SPIN ORBIT ALLIGNMENT (-1,0,1)
      INTEGER(4),INTENT(IN) :: NC        ! #(CORE STATES)
      INTEGER(4),INTENT(IN) :: NJ        ! #(PARTIAL WAVES)
      REAL(8)   ,INTENT(IN) :: RC        ! PSEUDIZATION RADIUS
      REAL(8)   ,INTENT(INOUT) :: EOFPHI(NJ)! EXPANSION ENERGY FOR PARTIAL WAVES
      REAL(8)   ,INTENT(IN) :: ECOREIN(NC) ! CORE LEVEL ENERGIES
      REAL(8)   ,INTENT(IN) :: AEPOT(NR) ! ALL-ELECTRON POTENTIAL
      REAL(8)   ,INTENT(IN) :: PSPOT(NR) ! PSEUDO POTENTIAL
      TYPE(VFOCK_TYPE),INTENT(IN) :: VFOCK
!
      REAL(8)   ,INTENT(OUT):: UCORE(NR,NC)  ! NODELESS CORE WAVE FUNCTIONS
      REAL(8)   ,INTENT(OUT):: AECORE(NR,NC) ! ALL-ELECTRON CORE WAVE FUNCTIONS
      REAL(8)   ,INTENT(OUT):: PSCORE(NR,NC) ! PSEUDO CORE WAVE FUNCTIONS
      REAL(8)   ,INTENT(OUT):: UCORESM(NR,NC) ! SMALL COMPONENT OF CORE STATES
      REAL(8)   ,INTENT(OUT):: AECORESM(NR,NC) ! SMALL AE CORE WAVE FUNCTIONS
      REAL(8)   ,INTENT(OUT):: PSCORESM(NR,NC) ! SMALL PS CORE WAVE FUNCTIONS
!
      REAL(8)   ,INTENT(OUT):: QN(NR,NJ)      ! NODE REDUCED PARTIAL WAVES
      REAL(8)   ,INTENT(OUT):: QNSM(NR,NJ)    ! SMALL QN
      REAL(8)   ,INTENT(OUT):: AEPHI(NR,NJ)
      REAL(8)   ,INTENT(OUT):: AEPHISM(NR,NJ)
      REAL(8)   ,INTENT(OUT):: PSPHI(NR,NJ)
      REAL(8)   ,INTENT(OUT):: PSPHISM(NR,NJ)
!
      REAL(8)   ,INTENT(OUT):: QNDOT(NR)      
      REAL(8)   ,INTENT(OUT):: QNDOTSM(NR)      
      REAL(8)   ,INTENT(OUT):: AEPHIDOT(NR)      
      REAL(8)   ,INTENT(OUT):: AEPHIDOTSM(NR)      
      REAL(8)   ,INTENT(OUT):: PRO(NR,NJ)
      REAL(8)   ,INTENT(OUT):: DOVER(NJ,NJ)
      REAL(8)   ,INTENT(OUT):: DTKIN(NJ,NJ)
      REAL(8)   ,PARAMETER  :: RNSCORE=0.07D0 !SEE MASTERS THESIS ROBERT SCHADE
      REAL(8)   ,PARAMETER  :: RNSPHI=0.09D0  !SEE MASTERS THESIS ROBERT SCHADE
      INTEGER(4),PARAMETER  :: SWITCH=2 ! BIORTHOGONALIZATION METHOD
      INTEGER(4),PARAMETER  :: IDIR=1
      LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
      LOGICAL(4),PARAMETER  :: TPRODOT=.FALSE. ! EXPERIMENTAL OPTION
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      LOGICAL(4)            :: TOLD
      LOGICAL(4)            :: TREL     ! RELATIVISTIC EFFECTS INCLUDED
      LOGICAL(4)            :: TZORA    ! ZEROTH-ORDER RELATIVISTIC EFFECTS
      LOGICAL(4)            :: TSMALL   ! SMALL COMPONENTS CALCULATED
      LOGICAL               :: TVARDREL ! ADJUST DREL IN BOUND-STATE SEARCH
      REAL(8)               :: KAPPA
      REAL(8)               :: SPEEDOFLIGHT
      REAL(8)               :: ALPHA    ! 1/C (=FINE-STRUCTURE CONSTANT IN A.U.)
      REAL(8)               :: ECORE(NC)
      REAL(8)               :: TUCORE(NR,NC)
      REAL(8)               :: TAECORE(NR,NC)
      REAL(8)               :: TPSCORE(NR,NC)
      REAL(8)               :: TQN(NR,NJ)
      REAL(8)               :: TQNDOT(NR)
      REAL(8)               :: PSPHIDOT(NR)
      REAL(8)               :: TPSPHIDOT(NR)
      REAL(8)               :: PRODOT(NR)
      REAL(8)               :: TAEPHI(NR,NJ)
      REAL(8)               :: TPSPHI(NR,NJ)
      REAL(8)   ,ALLOCATABLE:: PHITEST(:,:),TPHITEST(:,:),DEV(:,:)
      REAL(8)               :: DREL(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)               :: G(NR),GSM(NR)
      REAL(8)               :: MAT(NJ,NJ),MATINV(NJ,NJ)
      REAL(8)               :: A(NC,NJ),B(NC,NJ)
      REAL(8)               :: SVAR,SVAR1,SVAR2
      REAL(8)               :: RNS
      INTEGER(4)            :: IPHISCALE  
      REAL(8)               :: SCALE  
      REAL(8)               :: E
      REAL(8)               :: ENU   !  MAX(EOFPHI)
      INTEGER(4)            :: IRCL ! GRID POINT BEYOND CLASSICAL TURNING POINT
      INTEGER(4)            :: JP,J,IR,K,M,I,IB
      REAL(8)               :: FACTOR
!     **************************************************************************
                                       CALL TRACE$PUSH('SETUP_NEWPRO')
      ECORE=ECOREIN
!     == KAPPA=-L-1 FOR L*S.GE.0; KAPPA=L FOR L*S<0; KAPPA=-1 FOR SO=0 =========
      KAPPA=REAL( -1+SO*(-L+(SO-1)/2) ,KIND=8)
      CALL RADIAL$R(GID,NR,R)
      CALL CONSTANTS$GET('C',SPEEDOFLIGHT)
      ALPHA=1.D0/SPEEDOFLIGHT ! FINE STRUCTURE CONSTANT IN A.U.
      ENU=0.D0
      IF(NJ.GT.0) ENU=MAXVAL(EOFPHI)  ! USED TO CUT RELATIVISTIC EFFECTS OFF
                                   ! AND AS REFERENCE ENERGY FOR QNDOT FUNCTIONS
      IPHISCALE=0
!
!     ==========================================================================
!     == RESOLVE RELTYPE                                                      ==
!     == TSMALL=TVARDREL=(TREL.AND.(NOT.TZORA))                               ==
!     ==========================================================================
      IF(RELTYPE.EQ.'NONREL') THEN
        TREL=.FALSE.
        TZORA=.FALSE.
        TSMALL=.FALSE.   !SMALL COMPONENT IS ZERO
        TVARDREL=.FALSE.
      ELSE IF(RELTYPE.EQ.'ZORA') THEN
        TREL=.TRUE.
        TZORA=.TRUE.
        TSMALL=.FALSE.
        TVARDREL=.FALSE.
      ELSE IF(RELTYPE.EQ.'SPINORBIT') THEN
        TREL=.TRUE.
        TZORA=.FALSE.
        TSMALL=.TRUE.
        TVARDREL=.TRUE.
      ELSE
        CALL ERROR$MSG('RELTYPE NOT RECOGNIZED')
        CALL ERROR$MSG('MUST BE "NONREL", "ZORA", OR "SPINORBIT"')
        CALL ERROR$CHVAL('RELTYPE',RELTYPE)
        CALL ERROR$STOP('SETUP_NEWPRO')
      END IF
      IF(.NOT.TSMALL.AND.SO.NE.0) THEN
        CALL ERROR$MSG('SPIN ORBIT REQUIRES SMALL COMPONENTS')
        CALL ERROR$CHVAL('RELTYPE',RELTYPE)
        CALL ERROR$I4VAL('SO',SO)
        CALL ERROR$STOP('SETUP_NEWPRO')
      END IF
!
!     ==========================================================================
!     == FIND CLASSICAL TURNING POINT TO SWITCH OFF RELATIVISTIC EFFECTS      ==
!     ==========================================================================
      IRCL=0
      DO IR=1,NR
        IF(AEPOT(IR).GT.ENU/Y0) THEN
          IRCL=IR
          EXIT
        END IF
      ENDDO
!
!     == IF NO CLASSICAL TURNING POINT IS FOUND, THE KINETIC ENERGY IS ALWAYS ==
!     == POSITIVE. IN THIS CASE, SET RC TO 1 ABOHR =============================
      IF(IRCL.EQ.0) THEN
        DO IR=1,NR
          IRCL=IR
          IF(R(IR).GT.1.D0) EXIT
        ENDDO
      END IF
!
!     ==========================================================================
!     == CONSTRUCT RELATIVISTIC FACTOR                                        ==
!     == KEPT CONSTANT EXCEPT FOR FULLY RELATIVISTIC CALCULATION              ==
!     ==========================================================================
      DREL=0.D0   !NON-RELATIVISTIC CASE
      IF(TREL) THEN
!       == THE ENERGY (0.D0) MUST BE CONSISTENT WITH ATOMLIB$AESCF =============
        CALL SCHROEDINGER$DREL(GID,NR,AEPOT,0.D0,DREL)
        DREL(IRCL:)=0.D0
      END IF
!
!     ==========================================================================
!     == CONSTRUCT NODELESS CORE STATES                                       ==
!     == CAUTION!                                                             ==
!     == FOR THE FULLY RELATIVISTIC CALCULATION WITH THE SMALL COMPONENT,     ==
!     == THE ENERGIES FROM THE NODELESS CONSTRUCTION DEVIATE FROM THOSE       ==
!     == OBTAINED WITH THE NODAL EQUATION . THE DEVIATION FOR AU-S STATES     ==
!     == IS IN THE RANGE 10^-4 H, WHILE 10^-9 H IS ACHIEVED FOR ZORA.         ==
!     ==========================================================================
      G(:)=0.D0
      GSM(:)=0.D0
      DO IB=1,NC
        E=ECORE(IB)
        IF(TREL.AND.(.NOT.TZORA)) THEN
          CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
          DREL(IRCL:)=0.D0
        END IF
!
!       ========================================================================
!       == PREPARE INHOMOGENEITY ===============================================
!       ========================================================================
        IF(TSMALL) THEN
          AUX=0.5D0*ALPHA*(1.D0+DREL)*GSM   !GSM=-UCORE(IB-1)
          CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
          AUX=AUX1+(1.D0-KAPPA)/R*AUX    !CAUTION: PROBABLY SIGN ERROR
          AUX(1)=AUX(2)
          G=G-AUX !GSM=-FSM(IB-1)   !CAUTION: PROBABLY SIGN ERROR
!
          IF(IB.GT.1) THEN ! THIS IS A SMALL NUMERICAL CORRECTION ==============
            CALL RADIAL__DERIVE(GID,NR,UCORE(:,IB-1),AUX1)
            CALL RADIAL__DERIVE(GID,NR,AUX1,AUX)
            CALL RADIAL__VERLETD2(GID,NR,UCORE(:,IB-1),AUX1)
            G=G-(0.5D0*ALPHA*(1.D0+DREL))**2*(AUX1-AUX)
          END IF  
        END IF
!
!       == OBTAIN LARGE COMPONENT ==============================================
        RNS=0.D0
        IF(IB.GT.1.AND.TSMALL) RNS=RNSCORE ! AVOID SPURIOUS ZEROS NEAR ORIGIN
        CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RNS,ROUT,TVARDREL &
     &                         ,DREL,G,0,AEPOT,E,UCORE(:,IB))
        TUCORE(:,IB)=G+(E-AEPOT(:)*Y0)*UCORE(:,IB)
        ECORE(IB)=E
!
!       == CONSTRUCT SMALL COMPONENT ===========================================
        IF(TSMALL) THEN
          CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
          DREL(IRCL:)=0.D0
          CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,L,SO &
     &                                      ,DREL,GSM,UCORE(:,IB),UCORESM(:,IB))
          UCORESM(IRCL:,IB)=0.D0
        ELSE
          UCORESM(:,IB)=0.D0
        END IF
!
!       ========================================================================
        PRINT*,'CORE STATE ENERGY ',L,IB,ECORE(IB) &
                               ,' OLD ',ECOREIN(IB),ECORE(IB)-ECOREIN(IB)
!
!       == PROVIDE WAVE FUNCTIONS FOR THE NEXT NODELESS LEVEL ==================
        G  =-UCORE(:,IB)
        GSM=-UCORESM(:,IB)
      ENDDO
!
!     ==========================================================================
!     == NODE-LESS PARTIAL WAVES                                              ==
!     ==========================================================================
!     == GO BACK TO THE DIRAC EQUATION TO SEE HOW THE INHOMGENEITY TRANSLATES
      IF(TREL) THEN
!       == THE ENERGY (0.D0) MUST BE CONSISTENT WITH ATOMLIB$AESCF =============
        CALL SCHROEDINGER$DREL(GID,NR,AEPOT,0.D0,DREL)
        DREL(IRCL:)=0.D0
      END IF
      G=0.D0
      GSM=0.D0
      IF(NC.NE.0) THEN
        G=-UCORE(:,NC)
        GSM=-UCORESM(:,NC)
      END IF
      DO JP=1,NJ
        J=JP-1
!
!       == PREPARE INHOMOGENEITY ===============================================
        IF(TSMALL) THEN
          AUX=0.5D0*ALPHA*(1.D0+DREL)*GSM
          CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
          G=G-AUX1-(1.D0-KAPPA)/R*AUX
        END IF
!
!       == OBTAIN LARGE COMPONENT ==============================================
!        IF(JP.EQ.1) THEN
        IF(JP.GT.0) THEN
          RNS=RNSPHI ! AVOID SPURIOUS ZEROS NEAR ORIGIN
          E=EOFPHI(JP)
          CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RNS,ROUT,.FALSE. &
     &                           ,DREL,G,0,AEPOT,E,QN(:,JP))
          EOFPHI(JP)=E
        ELSE
          E=EOFPHI(JP)
          CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,SO,G,L,E,IDIR,QN(:,JP))
        END IF
        TQN(:,JP)=G-(AEPOT*Y0-E)*QN(:,JP)
!
!       == CONSTRUCT SMALL COMPONENT ===========================================
        IF(TSMALL) THEN
          CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,L,SO &
     &                                      ,DREL,GSM,QN(:,JP),QNSM(:,JP))
          QNSM(IRCL:,JP)=0.D0
        ELSE
          QNSM(:,JP)=0.D0
        END IF
!  
!       ========================================================================
        SCALE=1.D0/MAXVAL(ABS(QN(:,JP)))
        IPHISCALE=JP
!
!       == PROVIDE WAVE FUNCTIONS FOR THE NEXT NODELESS LEVEL ==================
        G=-QN(:,JP)
        GSM=-QNSM(:,JP)
      ENDDO
!
!     ==========================================================================
!     == UPDATE WAVE FUNCTIONS WITH FOCK POTENTIAL                            ==
!     == THE FOCK CORRECTION MUST
!     ==========================================================================
      IF(VFOCK%TON) THEN
        G(:)=0.D0
        GSM(:)=0.D0
        DO IB=1,NC
          E=ECORE(IB)
          IF(TREL.AND.(.NOT.TZORA)) THEN
            CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
            DREL(IRCL:)=0.D0
          END IF
!
!         == PREPARE INHOMOGENEITY =============================================
          IF(TSMALL) THEN
            AUX=0.5D0*ALPHA*(1.D0+DREL)*GSM   
            CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
            AUX=AUX1+(1.D0-KAPPA)/R*AUX  
            AUX(1)=AUX(2)
            G=G-AUX !GSM=-FSM(IB-1)
          END IF
!
!         == OBTAIN LARGE COMPONENT ============================================
          CALL ATOMLIB$UPDATESTATEWITHHF(GID,NR,L,SO,DREL,G,AEPOT,VFOCK &
     &                                    ,ROUT,E,UCORE(:,IB))
          CALL RADIALFOCK$VPSI(GID,NR,VFOCK,L,UCORE(:,IB),AUX)
          TUCORE(:,IB)=G+(E-AEPOT(:)*Y0)*UCORE(:,IB)-AUX(:)
          ECORE(IB)=E
!
!         == CONSTRUCT SMALL COMPONENT =========================================
          IF(TSMALL) THEN
            CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,L,SO &
     &                                      ,DREL,GSM,UCORE(:,IB),UCORESM(:,IB))
            UCORESM(IRCL:,IB)=0.D0
          ELSE
            UCORESM(:,IB)=0.D0
          END IF
!
!         == PROVIDE WAVE FUNCTIONS FOR THE NEXT NODELESS LEVEL ================
          G(:) =-UCORE(:,IB) 
          GSM(:)=-UCORESM(:,IB) 
        ENDDO
      END IF
!
!     ==========================================================================
!     == UPDATE NODELESS PARTIAL WAVES WITH FOCK POTENTIAL                    ==
!     ==========================================================================
      IF(VFOCK%TON) THEN
        IF(TREL) THEN
!         == THE ENERGY (0.D0) MUST BE CONSISTENT WITH ATOMLIB$AESCF ===========
          CALL SCHROEDINGER$DREL(GID,NR,AEPOT,0.D0,DREL)
          DREL(IRCL:)=0.D0
        END IF
        G(:)=0.D0
        GSM(:)=0.D0
        IF(NC.NE.0) THEN
          G(:)=UCORE(:,NC)
          GSM(:)=UCORESM(:,NC)
        END IF
        DO JP=1,NJ
!
!         == PREPARE INHOMOGENEITY =============================================
          IF(TSMALL) THEN
            AUX=0.5D0*ALPHA*(1.D0+DREL)*GSM
            CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
            G=G-AUX1-(1.D0-KAPPA)/R*AUX
          END IF
!
!         == OBTAIN LARGE COMPONENT ============================================
          RNS=RNSPHI ! AVOID SPURIOUS ZEROS NEAR ORIGIN
          E=EOFPHI(JP)
          CALL ATOMLIB$UPDATESTATEWITHHF(GID,NR,L,SO,DREL,G,AEPOT,VFOCK &
    &                                    ,ROUT,E,QN(:,JP))
          CALL RADIALFOCK$VPSI(GID,NR,VFOCK,L,QN(:,JP),AUX)
          TQN(:,JP)=G(:)+(E-AEPOT(:)*Y0)*QN(:,JP)-AUX(:)
          EOFPHI(JP)=E
!
!         == CONSTRUCT SMALL COMPONENT =========================================
          IF(TSMALL) THEN
            CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,L,SO &
     &                                      ,DREL,GSM,QN(:,JP),QNSM(:,JP))
            QNSM(IRCL:,JP)=0.D0
          ELSE
            QNSM(:,JP)=0.D0
          END IF
!
!         == PROVIDE WAVE FUNCTIONS FOR THE NEXT NODELESS LEVEL ================
          G(:)=-QN(:,JP)
          GSM(:)=-QNSM(:,JP)
        ENDDO
      END IF
!
!     ==========================================================================
!     == RESCALE CONSISTENTLY                                                 ==
!     ==========================================================================
!     == MAX ABSOLUTE VALUE OF LOWEST VALENCE FUNCTION =1 ======================
      IF(NJ.GT.0) THEN
        SCALE=1.D0/MAXVAL(ABS(QN(:,1)))
      ELSE
        IF(NC.GT.0) THEN
          SCALE=1.D0/MAXVAL(ABS(UCORE(:,NC)))
        ELSE
          SCALE=1.D0
        END IF
      END IF
      QN=QN*SCALE
      TQN=TQN*SCALE
      QNSM=QNSM*SCALE
      UCORE=UCORE*SCALE
      TUCORE=TUCORE*SCALE
      UCORESM=UCORESM*SCALE
      SCALE=1.D0
!
!     ==========================================================================
!     == CONSTRUCT QNDOT FUNCTION                                             ==
!     ==========================================================================
!!$      IF(NJ.EQ.0) THEN
!!$     CALL ERROR$MSG('ZERO NUMBER OF PARTIAL WAVES FOR THIS ANGULAR MOMENTUM')
!!$        CALL ERROR$MSG('THE CODE IS NOT YET ABLE TO DEAL WITH THIS')
!!$        CALL ERROR$MSG('BECAUSE OF THE STATEMENT: G=-QN(:,NJ)')
!!$        CALL ERROR$I4VAL('L',L)
!!$        CALL ERROR$I4VAL('SO',SO)
!!$        CALL ERROR$STOP('SETUP_NEWPRO')
!!$      END IF
      IF(NJ.GT.0) THEN
        G=-QN(:,NJ)
        IF(TSMALL) GSM=-QNSM(:,NJ)
      ELSE
        IF(NC.GT.0) THEN
          G=-UCORE(:,NC)
          IF(TSMALL) GSM=-UCORESM(:,NC)
        ELSE
          G=0.D0
          IF(TSMALL) GSM=0.D0
          IF(NC.EQ.0) THEN
            CALL ERROR$MSG('ZERO NUMBER OF CORE STATES AND')
            CALL ERROR$MSG('ZERO NUMBER OF VALENCE STATES FOR THIS L,SO')
            CALL ERROR$MSG('THE CODE IS NOT YET ABLE TO DEAL WITH THIS')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$I4VAL('SO',SO)
            CALL ERROR$I4VAL('#(CORE STATES)',NC)
            CALL ERROR$I4VAL('#(VALENCE STATES)',NJ)
            CALL ERROR$STOP('SETUP_NEWPRO')
          END IF
        END IF
      END IF
      IF(TSMALL) THEN
        AUX=0.5D0*ALPHA*(1.D0+DREL)*GSM
        CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
        G=G-AUX1-(1.D0-KAPPA)/R*AUX
      END IF
      CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,SO,G,L,ENU,IDIR,QNDOT)
      TQNDOT=G(:)+(ENU-AEPOT(:)*Y0)*QNDOT
      IF(TSMALL) THEN
        CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,L,SO,DREL,GSM,QNDOT,QNDOTSM)
        QNDOTSM(IRCL:)=0.D0
      ELSE
        QNDOTSM(:)=0.D0
      END IF
!
!     ==========================================================================
!     == PSEUDO CORE WAVE FUNCTIONS                                           ==
!     ==========================================================================
      PSCORE=UCORE
      TPSCORE=TUCORE
      PSCORESM=UCORESM
      DO IB=1,NC
        CALL SETUP_MAKEPSCORE(GID,NR,RC,L,IB-1,PSCORE(:,IB),TPSCORE(:,IB))
        AUX=0.D0
        CALL SETUP_MAKEPSCORE(GID,NR,RC,L,IB-1,PSCORESM(:,IB),AUX)
      ENDDO
!
!     ==========================================================================
!     == ORTHOGONALIZE ALL-ELECTRON CORE STATES                               ==
!     ==========================================================================
      AECORE(:,:)  =UCORE(:,:)
      TAECORE(:,:) =TUCORE(:,:)
      AECORESM(:,:)=UCORESM(:,:)
      IF(1.EQ.1) THEN  !ALTERNATIVE 1
        DO IB=2,NC
          SVAR=1.D0
          DO M=IB-1,1,-1
            SVAR=SVAR/(ECORE(M)-ECORE(IB))
            AECORE(:,IB)  =  AECORE(:,IB)+   UCORE(:,M)*SVAR
            TAECORE(:,IB) = TAECORE(:,IB)+  TUCORE(:,M)*SVAR
            AECORESM(:,IB)=AECORESM(:,IB)+ UCORESM(:,M)*SVAR
            PSCORE(:,IB)  =  PSCORE(:,IB)+  PSCORE(:,M)*SVAR
            TPSCORE(:,IB) = TPSCORE(:,IB)+ TPSCORE(:,M)*SVAR
            PSCORESM(:,IB)=PSCORESM(:,IB)+PSCORESM(:,M)*SVAR
          ENDDO
        ENDDO
      ELSE    !ALTERNATIVE 2
        DO I=1,NC
          AUX(:)=R**2*(AECORE(:,I)*AECORE(:,I)+AECORESM(:,I)*AECORESM(:,I))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR1)
          DO J=I+1,NC
            AUX(:)=R**2*(AECORE(:,I)*AECORE(:,J)+AECORESM(:,I)*AECORESM(:,J))
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR2)
            SVAR=SVAR2/SVAR1
            AECORE(:,J)  =  AECORE(:,J)-  AECORE(:,I)*SVAR
            TAECORE(:,J) = TAECORE(:,J)- TAECORE(:,I)*SVAR
            AECORESM(:,J)=AECORESM(:,J)-AECORESM(:,I)*SVAR
            PSCORE(:,J)  =  PSCORE(:,J)-  PSCORE(:,I)*SVAR
            TPSCORE(:,J) = TPSCORE(:,J)- TPSCORE(:,I)*SVAR
            PSCORESM(:,J)=PSCORESM(:,J)-PSCORESM(:,I)*SVAR
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == CORE-VALENCE ORTHOGONALIZATION
!     ==========================================================================
      AEPHI=QN
      TAEPHI=TQN
      AEPHISM=QNSM
      AEPHIDOT=QNDOT
      AEPHIDOTSM=QNDOTSM
      IF(1.EQ.1) THEN
        DO JP=1,NJ
          FACTOR=1.D0
          DO I=JP-1,1,-1
            FACTOR=FACTOR/(EOFPHI(I)-EOFPHI(JP))
            AEPHI(:,JP)  =AEPHI(:,JP)  +(QN(:,I)-AEPHI(:,I))*FACTOR
            TAEPHI(:,JP) =TAEPHI(:,JP) +(TQN(:,I)-TAEPHI(:,I))*FACTOR
            AEPHISM(:,JP)=AEPHISM(:,JP)+(QNSM(:,I)-AEPHISM(:,I))*FACTOR
          ENDDO
          DO I=NC,1,-1
            FACTOR=FACTOR/(ECORE(I)-EOFPHI(JP))
            AEPHI(:,JP)  =AEPHI(:,JP)  +UCORE(:,I)  *FACTOR
            TAEPHI(:,JP) =TAEPHI(:,JP) +TUCORE(:,I) *FACTOR
            AEPHISM(:,JP)=AEPHISM(:,JP)+UCORESM(:,I)*FACTOR
          ENDDO
        ENDDO
!       == CORE ORTHOGONALIZE SCATTERING FUNCTIONS TRADITIONALLY ===============
        DO I=NC,1,-1
          AUX(:)=R**2*(AECORE(:,I)*AECORE(:,I)+AECORESM(:,I)*AECORESM(:,I))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR1)
          AUX(:)=R**2*(AECORE(:,I)*AEPHIDOT+AECORESM(:,I)*AEPHIDOTSM(:))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR2)
          SVAR=SVAR2/SVAR1
          AEPHIDOT(:)=AEPHIDOT(:)    -AECORE(:,I)  *SVAR
          AEPHIDOTSM(:)=AEPHIDOTSM(:)-AECORESM(:,I)*SVAR
        ENDDO
      ELSE
        DO I=NC,1,-1
          AUX(:)=R**2*(AECORE(:,I)*AECORE(:,I)+AECORESM(:,I)*AECORESM(:,I))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR1)
          DO JP=1,NJ
            AUX(:)=R**2*(AECORE(:,I)*AEPHI(:,JP)+AECORESM(:,I)*AEPHISM(:,JP))
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR2)
            SVAR2=SVAR2/SVAR1
            AEPHI(:,JP)  =AEPHI(:,JP)  -AECORE(:,I)  *SVAR2
            TAEPHI(:,JP) =TAEPHI(:,JP) -TAECORE(:,I) *SVAR2
            AEPHISM(:,JP)=AEPHISM(:,JP)-AECORESM(:,I)*SVAR2
          ENDDO
    !PB CAUTION: AEPHIDOT INCREASES EXPONENTIALLY WITH DISTANCE!!
    !PB THIS MAY BE NUMERICALLY PROBLEMATIC.
          AUX(:)=R**2*(AECORE(:,I)*AEPHIDOT+AECORESM(:,I)*AEPHIDOTSM(:))/SVAR1
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR2)
          AEPHIDOT(:)=AEPHIDOT(:)    -AECORE(:,I)  *SVAR2
          AEPHIDOTSM(:)=AEPHIDOTSM(:)-AECORESM(:,I)*SVAR2
        ENDDO
      END IF
!
!     ==========================================================================
!     == MAKE PSEUDO PARTIAL WAVE                                             ==
!     == (ONLY THE FIRST PSEUDO PARTIAL WAVE DIFFERS FROM THE QN)             ==
!     ==========================================================================
      PSPHI=QN
      TPSPHI=TQN
      PSPHISM=QNSM
!     == ONLY THE FIRST PARTIAL WAVE IS PSEUDIZED ==============================
!     == TAILS OF THE PSEUDO CORE FUNCTIONS WILL BE ADDED LATER  ===============
      IF(NJ.GT.0) CALL SETUP_MAKEPSPHI_FLAT(GID,NR,RC,L,PSPHI,TPSPHI)
!
!     == CALCULATE NON-RELATIVISTIC KINETIC ENERGY =============================
!     == FOR SMALL DISTANCES I TRY TO AVOID THE SINGULARITIES 
!     == AUX = -1/2*[1/R PARTIAL_R^2 R- L(L+1)/R^2]|Q> 
!     ==     = -1/2*R^L [ 2(L+1)/R+PARTIAL_R ] PARTIAL_R ( QN/R^L )
      DO JP=2,NJ
TOLD=.FALSE.
IF(TOLD) THEN
!!$! CHANGED ON JAN. 1, 2015 TO BE CONSISTENT WITH SCHROEDINGER EQUATION ======
        AUX(:)=PSPHI(:,JP)/R(:)**L
        CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
        CALL RADIAL$DERIVE(GID,NR,AUX1,AUX)
        AUX(:)=AUX(:)+2.D0*REAL(L+1,KIND=8)*AUX1(:)/R(:)
        TPSPHI(:,JP)=-0.5D0*R(:)**L*AUX(:)
ELSE
!
!       == THIS SHOULD BE CONSISTENT WITH THE SOLVER FOR THE SCHR. EQ.
        CALL RADIAL$VERLETD1(GID,NR,PSPHI(:,JP),AUX)
        CALL RADIAL$VERLETD2(GID,NR,PSPHI(:,JP),AUX1)
        TPSPHI(:,JP)=-0.5D0*(AUX1+2.D0*AUX/R &
     &                           -REAL(L*(L+1),KIND=8)*PSPHI(:,JP)/R**2)
        TPSPHI(1:2,JP)=TPSPHI(3,JP)
END IF
      ENDDO
      IF(TPRODOT) THEN
        PSPHIDOT=QNDOT
        TPSPHIDOT=TQNDOT
        CALL RADIAL$VERLETD1(GID,NR,PSPHIDOT,AUX)
        CALL RADIAL$VERLETD2(GID,NR,PSPHIDOT,AUX1)
        TPSPHIDOT=-0.5D0*(AUX1+2.D0*AUX/R-REAL(L*(L+1),KIND=8)*PSPHIDOT/R**2)
        TPSPHIDOT(1:2)=TPSPHIDOT(3)
      END IF
!
!     == FOR LARGE DISTANCES (R>RC), I USE THE RELATIVISTIC KINETIC ENERGY =====
!     == TO ENSURE CANCELLATION OF TAIL CONTRIBUTIONS. =========================
      DO IR=1,NR
        IF(R(IR).GT.RC) THEN
          TPSPHI(IR:,2:)=TQN(IR:,2:)
          IF(TPRODOT) TPSPHIDOT(IR:)=TQNDOT(IR:)
          EXIT
        END IF
      ENDDO
!
!     ==========================================================================
!     == RESCALE CONSISTENTLY: WAVE FUNCTION PHI(IPHISCALE) WILL HAVE NORM 1  ==
!     ==========================================================================
!!$PRINT*,'IPHISCALE ',IPHISCALE
      IPHISCALE=1
      IF(NJ.NE.0) THEN
        AUX=R**2*(AEPHI(:,IPHISCALE)**2+AEPHISM(:,IPHISCALE)**2)
        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
        CALL RADIAL$VALUE(GID,NR,AUX1,3.D0,SCALE)
!PB THIS NORMALIZATION RADIUS IS PHYSICALLY INSIGNIFICANT. THE RADIUS OF 3 A.U
!PB ALLOWS BETTER COMPARISON.
        SCALE=1.D0/SQRT(SCALE)
!       == PARTIAL WAVES =======================================================
        QN=QN*SCALE
        TQN=TQN*SCALE
        QNSM=QNSM*SCALE
        AEPHI=AEPHI*SCALE
        TAEPHI=TAEPHI*SCALE
        AEPHISM=AEPHISM*SCALE
        PSPHI=PSPHI*SCALE
        TPSPHI=TPSPHI*SCALE
        PSPHISM=PSPHISM*SCALE
!       == CORE STATES =========================================================
        UCORE=UCORE*SCALE
        TUCORE=TUCORE*SCALE
        UCORESM=UCORESM*SCALE
        AECORE=AECORE*SCALE
        TAECORE=TAECORE*SCALE
        AECORESM=AECORESM*SCALE
        PSCORE=PSCORE*SCALE
        TPSCORE=TPSCORE*SCALE
        PSCORESM=PSCORESM*SCALE
!       == SCATTERING WAVE FUNCTIONS ===========================================
        QNDOT=QNDOT*SCALE
        TQNDOT=TQNDOT*SCALE
        QNDOTSM=QNDOTSM*SCALE
        AEPHIDOT=AEPHIDOT*SCALE
        AEPHIDOTSM=AEPHIDOTSM*SCALE
        IF(TPRODOT) THEN
          PSPHIDOT=PSPHIDOT*SCALE
          TPSPHIDOT=TPSPHIDOT*SCALE
        END IF
      END IF
!
!     ==========================================================================
!     == CONSTRUCT BARE PROJECTOR FUNCTIONS                                   ==
!     ==========================================================================
      IF(NJ.GT.0) PRO(:,1)=TPSPHI(:,1)+(PSPOT*Y0-EOFPHI(1))*PSPHI(:,1)
      DO JP=2,NJ
        PRO(:,JP)=TPSPHI(:,JP)-TQN(:,JP)+(PSPOT-AEPOT)*Y0*QN(:,JP)
!        PRO(:,JP)=TPSPHI(:,JP)+(PSPOT*Y0-EOFPHI(JP))*PSPHI(:,JP)
      ENDDO
      IF(NJ.GE.2) PRO(:,2)=PRO(:,2)+(PSPHI(:,1)-QN(:,1))  ! |K>
!
      IF(TPRODOT) THEN
        PRODOT(:)=TPSPHIDOT(:)-TQNDOT(:)+(PSPOT-AEPOT)*Y0*QNDOT
        IF(NJ.EQ.1) PRODOT(:)=PRODOT(:)+(PSPHI(:,1)-QN(:,1))  ! |K>
      END IF
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST.AND.L.EQ.0) THEN
        ALLOCATE(DEV(NR,NJ))
        ALLOCATE(PHITEST(NR,NJ))
        ALLOCATE(TPHITEST(NR,NJ))
!
!       ==  TEST NODELESS DIFFERENTIAL EQUATION ================================
        DEV=0.D0
        DO J=1,NJ
          DEV(:,J)=TQN(:,J)+(AEPOT*Y0-EOFPHI(J))*QN(:,J)
          IF(J.EQ.1) THEN
            DEV(:,J)=DEV(:,J)+UCORE(:,NC)
          ELSE
            DEV(:,J)=DEV(:,J)+QN(:,J-1)
          END IF
        ENDDO
        IF(MAXVAL(ABS(DEV)).LT.1.D-6) THEN
          PRINT*,'OK: TESTING NODELESS EQUATION PASSED'
        ELSE
          PRINT*,'WARNING: TESTING NODELESS EQUATION FAILED!'
        END IF
        CALL SETUP_WRITEPHI(-'TEST_QNDEV.DAT',GID,NR,NJ,DEV)
        CALL SETUP_WRITEPHI(-'TEST_QN.DAT',GID,NR,NJ,QN)
!
!       == TEST ALL-ELECTRON SCHROEDINGER EQ ===================================
        DEV=0.D0
        DO J=1,NJ
          CALL SETUP_PHITPHIOFE(GID,NR,NJ,EOFPHI,AEPHI,TAEPHI,EOFPHI(J) &
      &                     ,PHITEST(:,J),TPHITEST(:,J))
          DEV(:,J)=TPHITEST(:,J)+(AEPOT*Y0-EOFPHI(J))*PHITEST(:,J)
        ENDDO
        IF(MAXVAL(ABS(DEV)).LT.1.D-6) THEN
          PRINT*,'OK: TESTING ALL-ELECTRON EQUATION PASSED'
        ELSE
          PRINT*,'WARNING: TESTING ALL-ELECTRON EQUATION FAILED!'
        END IF
        CALL SETUP_WRITEPHI(-'TEST_AEPHIDEV.DAT',GID,NR,NJ,DEV)
        CALL SETUP_WRITEPHI(-'TEST_AEPHI.DAT',GID,NR,NJ,PHITEST)
!
!       == TEST PSEUDO-SCHROEDINGER EQ ===================================
        DEV=0.D0
        DO J=1,NJ
          CALL SETUP_PHITPHIOFE(GID,NR,NJ,EOFPHI,PRO,TPSPHI,EOFPHI(J) &
      &                     ,DEV(:,J),TPHITEST(:,J))
        ENDDO
        DO J=1,NJ
          CALL SETUP_PHITPHIOFE(GID,NR,NJ,EOFPHI,PSPHI,TPSPHI,EOFPHI(J) &
      &                     ,PHITEST(:,J),TPHITEST(:,J))
          DEV(:,J)=TPHITEST(:,J)+(PSPOT*Y0-EOFPHI(J))*PHITEST(:,J)-DEV(:,J)
        ENDDO
        IF(MAXVAL(ABS(DEV)).LT.1.D-6) THEN
          PRINT*,'OK: TESTING PSEUDO EQUATION PASSED'
        ELSE
          PRINT*,'WARNING: TESTING PSEUDO EQUATION FAILED!'
        END IF
        CALL SETUP_WRITEPHI(-'TEST_PSPHIDEV.DAT',GID,NR,NJ,DEV)
        CALL SETUP_WRITEPHI(-'TEST_PSPHI.DAT',GID,NR,NJ,PHITEST)
        CALL SETUP_WRITEPHI(-'TEST_PRO.DAT',GID,NR,NJ,PRO)
!
!       == TEST CORE-VALENCE ORTHOGONALITY <AECORE|AEPHI> ======================
        WRITE(*,FMT='("CORE ORTHOGONALIZATION TEST")')
        DO I=1,NC
          AUX=R**2*AECORE(:,I)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR1)
          DO J=1,NJ
            AUX=R**2*AECORE(:,I)*AEPHI(:,J)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR)
            SVAR=SVAR/SQRT(SVAR1)
            WRITE(*,FMT='("IC=",I2," IB=",I2," <C|AEPHI>=",F20.10)')I,J,SVAR
          ENDDO
        ENDDO
        CALL SETUP_WRITEPHI(-'TEST_AECORE.DAT',GID,NR,NC,AECORE)
!
!       == TEST CORE-CORE ORTHONORMALITY =======================================
        WRITE(*,FMT='("TEST ORTHOGONALITY OF CORE STATES")')
        DO I=1,NC
          AUX=R**2*AECORE(:,I)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR1)
          DO J=I+1,NC
            AUX=R**2*AECORE(:,J)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR2)
            AUX=R**2*AECORE(:,I)*AECORE(:,J)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR)
            SVAR=SVAR/SQRT(SVAR1*SVAR2)
            WRITE(*,FMT='("IC=",I2," IB=",I2," <C|C>=",F20.10)')I,J,SVAR
          ENDDO
        ENDDO
!
!       ==  RECALCULATE ENERGY-DEPENDENT PARTIAL WAVES FROM SCHROEDINGER EQ. ===
!       == PUT PROJECTOR FUNCTIONS IN DEV AND PSEUDO PARTIAL WAVES IN TPHITEST
        DEV=0.D0
        DO J=1,NJ
          CALL SETUP_PHITPHIOFE(GID,NR,NJ,EOFPHI,PRO,TPSPHI,EOFPHI(J) &
      &                     ,DEV(:,J),TPHITEST(:,J))
          CALL SETUP_PHITPHIOFE(GID,NR,NJ,EOFPHI,PSPHI,TPSPHI,EOFPHI(J) &
      &                     ,PHITEST(:,J),TPHITEST(:,J))
          TPHITEST(:,J)=PHITEST(:,J)
        ENDDO
        CALL SETUP_WRITEPHI(-'TEST_PROOFE.DAT',GID,NR,NJ,DEV)
        DO J=1,NJ
          DREL=0.D0
          E=EOFPHI(J)
!         == DETERMINE INHOMOGENEOUS SOLUTION ==================================
          G=DEV(:,J)
          CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,G,L,E,IDIR &
       &                             ,PHITEST(:,J))
          CALL SETUP_WRITEPHI(-'TEST_A.DAT',GID,NR,1,PHITEST(:,J))
!         == DETERMINE HOMOGENEOUS SOLUTION ====================================
          G=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,PSPOT,DREL,SO,G,L,E,IDIR,AUX)
          CALL SETUP_WRITEPHI(-'TEST_B.DAT',GID,NR,1,AUX)
!         == ADD INHOMOGENEOUS SOLUTION TO SATISFY BOUNDARY CONDITIONS =========
          CALL RADIAL$VALUE(GID,NR,PHITEST(:,J),ROUT,SVAR1)
          CALL RADIAL$VALUE(GID,NR,AUX,ROUT,SVAR2)
          PHITEST(:,J)=PHITEST(:,J)-AUX(:)*SVAR1/SVAR2
        ENDDO
!       == REMOVE IRRELEVANT DEVIATIONS
        DO IR=1,NR
          IF(R(IR).LE.ROUT) CYCLE
          DEV(IR:,:)=0.D0
        END DO
        IF(MAXVAL(ABS(DEV)).LT.1.D-6) THEN
          PRINT*,'OK: TESTING INHOM PAW_EQ FOR PSPHI PASSED'
        ELSE
          PRINT*,'WARNING: TESTING INHOM PAW_EQ FOR PSPHI FAILED!'
!         == THE SECOND WAVE FUNCTION DIFFERS BY A CONSTANT FACTOR.       
        END IF
        CALL SETUP_WRITEPHI(-'TEST_PHI1.DAT',GID,NR,NJ,PHITEST)
        CALL SETUP_WRITEPHI(-'TEST_PHI2.DAT',GID,NR,NJ,TPHITEST)
        CALL SETUP_WRITEPHI(-'TEST_PHI12DEV.DAT',GID,NR,NJ,DEV)
        DEALLOCATE(PHITEST)
        DEALLOCATE(TPHITEST)
        DEALLOCATE(DEV)
        CALL ERROR$MSG('FORCED STOP AFTER TESTING')
        CALL ERROR$STOP('SETUP_NEWPRO')
      END IF
!
!     ==========================================================================
!     == ORTHOGONALIZE BARE PROJECTOR FUNCTION TO HIGHER PARTIAL WAVE         ==
!     ==========================================================================
      IF(TPRODOT) THEN
!!$CALL SETUP_WRITEPHI(-'TEST_PRO_OLD.DAT',GID,NR,NJ,PRO)
!!$CALL SETUP_WRITEPHI(-'TEST_PSPHI.DAT',GID,NR,NJ,PSPHI)
!!$CALL SETUP_WRITEPHI(-'TEST_PRODOT_OLD.DAT',GID,NR,1,PRODOT)
!!$CALL SETUP_WRITEPHI(-'TEST_PSPHIDOT.DAT',GID,NR,1,PSPHIDOT)
!!$CALL SETUP_WRITEPHI(-'TEST_QNDOT.DAT',GID,NR,1,QNDOT)
        CALL RADIAL$INTEGRAL(GID,NR,R**2*PSPHIDOT*PRODOT,SCALE)
        SCALE=1.D0/SCALE
        DO JP=1,NJ
          CALL RADIAL$INTEGRAL(GID,NR,R**2*PSPHIDOT*PRO(:,JP),SVAR)
          PRO(:,JP)=PRO(:,JP)-PRODOT*SVAR*SCALE
        ENDDO
      END IF
!
!     ==========================================================================
!     == UNBARE PROJECTOR FUNCTIONS                                           ==
!     ==========================================================================
      IF(NJ.GT.0) THEN
        DO JP=1,NJ
          DO J=1,NJ
            AUX(:)=R(:)**2*PSPHI(:,J)*PRO(:,JP)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,MAT(J,JP))
          ENDDO
        ENDDO        
        CALL LIB$INVERTR8(NJ,MAT,MATINV)
        PRO(:,:)=MATMUL(PRO,MATINV)
      END IF
!
!     ==========================================================================
!     == ADJUST PRODOT (PRODOT IS ONLY FOR INTERNAL CONSISTENCY CHECKS)
!     ==========================================================================
      IF(TPRODOT) THEN
        DO JP=1,NJ
          CALL RADIAL$INTEGRAL(GID,NR,R**2*PSPHI(:,JP)*PRODOT,SVAR)
          PRODOT(:)=PRODOT(:)-PRO(:,JP)*SVAR
        ENDDO
        CALL RADIAL$INTEGRAL(GID,NR,R**2*PSPHIDOT*PRODOT,SVAR)
        PRODOT=PRODOT/SVAR
!TEST
!!$      DO J=1,NJ
!!$        DO JP=1,NJ
!!$          AUX(:)=R(:)**2*PRO(:,J)*PSPHI(:,JP)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$          PRINT*,'J,JP, BIORTHO ',J,JP,SVAR
!!$        ENDDO
!!$        AUX(:)=R(:)**2*PRO(:,JP)*PSPHIDOT
!!$        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$        PRINT*,'J,JP, BIORTHO ',J,'DOT',SVAR
!!$      ENDDO
!!$      DO JP=1,NJ
!!$        AUX(:)=R(:)**2*PRODOT*PSPHI(:,JP)
!!$        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$        PRINT*,'J,JP, BIORTHO ','DOT',JP,SVAR
!!$      ENDDO
!!$      AUX(:)=R(:)**2*PRODOT*PSPHIDOT
!!$      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$      PRINT*,'J,JP, BIORTHO ','DOT','DOT',SVAR
!!$CALL SETUP_WRITEPHI(-'TEST_PRO_NEW.DAT',GID,NR,NJ,PRO)
!!$CALL SETUP_WRITEPHI(-'TEST_PRODOT_NEW.DAT',GID,NR,1,PRODOT)
!!$STOP
      END IF
!
!     ==========================================================================
!     == CALCULATE DOVER AND DTKIN                                            ==
!     ==========================================================================
      DO JP=1,NJ
        DO J=1,NJ
          AUX(:)=AEPHI(:,J)*AEPHI(:,JP)-PSPHI(:,J)*PSPHI(:,JP) &
     &          +AEPHISM(:,J)*AEPHISM(:,JP)-PSPHISM(:,J)*PSPHISM(:,JP)
          CALL RADIAL$INTEGRAL(GID,NR,AUX*R**2,DOVER(J,JP))
          AUX(:)=AEPHI(:,J)*TAEPHI(:,JP)-PSPHI(:,J)*TPSPHI(:,JP)
          IF(TSMALL) THEN
             AUX(:)=AUX(:)+TAEPHI(:,J)*AEPHI(:,JP) &
     &                    -2.D0*SPEEDOFLIGHT**2*AEPHISM(:,J)*AEPHISM(:,JP) &
     &                    -TPSPHI(:,J)*AEPHI(:,JP) &
     &                    +2.D0*SPEEDOFLIGHT**2*PSPHISM(:,J)*PSPHISM(:,JP)
          END IF
          CALL RADIAL$INTEGRAL(GID,NR,AUX*R**2,DTKIN(J,JP))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == SYMMETRIZE DTKIN                                                     ==
!     == DTKIN IS NOT EXACTLY SYMMETRIC. THE CULPRIT IS THE PSEUDO PARTIAL    ==
!     == WAVE.                                                                ==
!     ==========================================================================
      DTKIN=0.5D0*(DTKIN+TRANSPOSE(DTKIN))
                                       CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_MAKEPSPHI_FLAT(GID,NR,RC,L,PHI,TPHI)
!     **************************************************************************
!     **  REPLACES PHI BY ITS PSEUDIZED VERSION                               **
!     **                                                                      **
!     **  THE PSEUDIZATION RADIUS IS THE GRID POINT JUST OUTSIDE RC.          **
!     **  REPLACE THE POTENTIAL INSIDE RC BY A CONSTANT WITH POT(IRC) AND     **
!     **  FIND THE PARTIAL WAVE AND ITS FIRST TWO ENERGY DERIVATIVES.         **
!     **  MATCH THEM WITH VALUE AND DERIVATIVE AND (FOR T3PAR=TRUE) ALSO      **
!     **  WITH KINETIC ENERGY                                                 **
!     **                                                                      **
!     **  NOTE, THAT FDDOT IS 1/2 OF THE SECOND DERIVATIVE                    **
!     **       (H-E)|F^J>=J|F^{J-1}>.  WE DROP THE FACTOR J                   **
!     ******************************PETER BLOECHL, GOSLAR 2013 *****************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: RC
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(INOUT) :: PHI(NR)
      REAL(8)   ,INTENT(INOUT) :: TPHI(NR)
      REAL(8)               :: DREL(NR)
      INTEGER(4),PARAMETER  :: IDIR=1
      INTEGER(4),PARAMETER  :: SO=0
      LOGICAL(4),PARAMETER  :: T3PAR=.TRUE.
      LOGICAL(4),PARAMETER  :: T4PAR=.TRUE.
      LOGICAL(4),PARAMETER  :: TWRITE=.FALSE.
!     == A LARGER POTPOW MAKES THE POTENTIAL FLATTER WITH A MORE ABRUPT UPTURN
      REAL(8)   ,PARAMETER  :: POTPOW=3.D0
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4)            :: IRC
      REAL(8)               :: R(NR)
      REAL(8)               :: PHIIN(NR),TPHIIN(NR)
      REAL(8)               :: POT(NR)
      REAL(8)               :: G(NR)
      REAL(8)               :: F(NR)
      REAL(8)               :: TF(NR)
      REAL(8)               :: FDOT(NR)
      REAL(8)               :: TFDOT(NR)
      REAL(8)               :: FDDOT(NR)
      REAL(8)               :: TFDDOT(NR)
      REAL(8)               :: FDDDOT(NR)
      REAL(8)               :: TFDDDOT(NR)
      REAL(8)               :: ENU
      REAL(8)               :: VALPHI,VALF,VALFDOT,VALFDDOT,VALFDDDOT,VAL
      REAL(8)               :: DERPHI,DERF,DERFDOT,DERFDDOT,DERFDDDOT
      REAL(8)               :: VALTPHI,VALTF,VALTFDOT,VALTFDDOT,VALTFDDDOT
      REAL(8)               :: DERTPHI,DERTF,DERTFDOT,DERTFDDOT,DERTFDDDOT
      REAL(8)               :: C1,C2,C3,C4,DET
      INTEGER(4)            :: IR
      REAL(8)               :: SVAR
      CHARACTER(16) :: LSTRING
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IF(R(IR).GT.RC) EXIT
        IRC=IR         !R(IRC)>RC
      ENDDO
!     == MAKE COPIES OF INPUT FUNCTIONS
      PHIIN=PHI
      TPHIIN=TPHI

!     ==========================================================================
!     == INITIAL REPORTING                                                    ==
!     ==========================================================================
      IF(TWRITE) THEN
        WRITE(LSTRING,*)L
        LSTRING=ADJUSTL(LSTRING)
        CALL SETUP_WRITEPHI('PHIBEFORE'//TRIM(LSTRING),GID,NR,1,PHI)
        CALL SETUP_WRITEPHI('TPHIBEFORE'//TRIM(LSTRING),GID,NR,1,TPHI)
      END IF
!
!     ==========================================================================
!     == EXTRACT POTENTIAL FROM INPUT PARTIAL WAVE                            ==
!     ==========================================================================
      DREL(:)=0.D0
      POT(:)=0.D0
      DO IR=IRC-5,NR-2  ! SUFFICIENT POINTS TO EXTRACT VALUE AND DERIVATIVE
                        ! LAST TWO POINTS ARE NAN
        POT(IR)=-TPHI(IR)/PHI(IR)/Y0
      ENDDO
      ENU=0.D0
!
!     == REPLACE INSIDE RC BY C1+C2*R^POTPOW ===================================
!     == POTPOW=3 IS INSPIRED BY TROULIER MARTINS
      F(:)   =1.D0
      FDOT(:)=R(:)**POTPOW
      CALL RADIAL$VALUE(GID,NR,POT,RC,VALPHI)
      CALL RADIAL$DERIVATIVE(GID,NR,POT,RC,DERPHI)
      CALL RADIAL$VALUE(GID,NR,F,RC,VALF)
      CALL RADIAL$DERIVATIVE(GID,NR,F,RC,DERF)
      CALL RADIAL$VALUE(GID,NR,FDOT,RC,VALFDOT)
      CALL RADIAL$DERIVATIVE(GID,NR,FDOT,RC,DERFDOT)
      DET=VALF*DERFDOT-VALFDOT*DERF
      C1=( DERFDOT*VALPHI-VALFDOT*DERPHI)/DET
      C2=(-DERF   *VALPHI+VALF   *DERPHI)/DET
      POT(:IRC)=F(:IRC)*C1+FDOT(:IRC)*C2
!
!     ==========================================================================
!     == DETERMINE RADIAL FUNCTIONS  FOR PSEUDO PARTIAL WAVES                 ==
!     ==========================================================================
!     == SOLUTIONS WILL ONLY BE USED UP TO MATCHING RADIUS. CUTTING THE ========
!     == THE POTENTIAL OFF EARLY, AVOIDS OVERFLOW ==============================
      POT(IRC+5:)=0.D0
      G=0.D0
      CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,ENU,IDIR,F)
      TF=G-(POT*Y0-ENU)*F
      G=F
!      G(IRC+1:)=0.D0
      CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,ENU,IDIR,FDOT)
      TFDOT=G-(POT*Y0-ENU)*FDOT
      G=2.D0*FDOT
!      G(IRC+1:)=0.D0
      CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,ENU,IDIR,FDDOT)
      TFDDOT=G-(POT*Y0-ENU)*FDDOT
      G=3.D0*FDDOT
!      G(IRC+1:)=0.D0
      CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,SO,G,L,ENU,IDIR,FDDDOT)
      TFDDDOT=G-(POT*Y0-ENU)*FDDDOT
!
!     ==========================================================================
!     == SUPERIMPOSE SO THAT IT MATCHES TO INPUT PARTIAL WAVE                 ==
!     ==========================================================================
      CALL RADIAL$VALUE(GID,NR,PHI,RC,VALPHI)
      CALL RADIAL$DERIVATIVE(GID,NR,PHI,RC,DERPHI)
      CALL RADIAL$VALUE(GID,NR,TPHI,RC,VALTPHI)
      CALL RADIAL$DERIVATIVE(GID,NR,TPHI,RC,DERTPHI)

      CALL RADIAL$VALUE(GID,NR,F,RC,VALF)
      CALL RADIAL$DERIVATIVE(GID,NR,F,RC,DERF)
      CALL RADIAL$VALUE(GID,NR,TF,RC,VALTF)
      CALL RADIAL$DERIVATIVE(GID,NR,TF,RC,DERTF)

      CALL RADIAL$VALUE(GID,NR,FDOT,RC,VALFDOT)
      CALL RADIAL$DERIVATIVE(GID,NR,FDOT,RC,DERFDOT)
      CALL RADIAL$VALUE(GID,NR,TFDOT,RC,VALTFDOT)
      CALL RADIAL$DERIVATIVE(GID,NR,TFDOT,RC,DERTFDOT)

      CALL RADIAL$VALUE(GID,NR,FDDOT,RC,VALFDDOT)
      CALL RADIAL$DERIVATIVE(GID,NR,FDDOT,RC,DERFDDOT)
      CALL RADIAL$VALUE(GID,NR,TFDDOT,RC,VALTFDDOT)
      CALL RADIAL$DERIVATIVE(GID,NR,TFDDOT,RC,DERTFDDOT)

      CALL RADIAL$VALUE(GID,NR,FDDDOT,RC,VALFDDDOT)
      CALL RADIAL$DERIVATIVE(GID,NR,FDDDOT,RC,DERFDDDOT)
      CALL RADIAL$VALUE(GID,NR,TFDDDOT,RC,VALTFDDDOT)
      CALL RADIAL$DERIVATIVE(GID,NR,TFDDDOT,RC,DERTFDDDOT)
      DET=VALF*DERFDOT-VALFDOT*DERF

!
!     ==  FIRST REMOVE VALUE AND DERIVATIVE FROM THE FDDOT FUNCTION ============
      C1=( DERFDOT*VALFDDOT-VALFDOT*DERFDDOT)/DET
      C2=(-DERF   *VALFDDOT+VALF   *DERFDDOT)/DET
      FDDOT(:) = FDDOT(:)- F(:)*C1- FDOT(:)*C2
      TFDDOT(:)=TFDDOT(:)-TF(:)*C1-TFDOT(:)*C2
      VALTFDDOT=VALTFDDOT-VALTF*C1-VALTFDOT*C2
      DERTFDDOT=DERTFDDOT-DERTF*C1-DERTFDOT*C2
!
!     ==  ... AND FROM THE FDDDOT FUNCTION =====================================
      C1=( DERFDOT*VALFDDDOT-VALFDOT*DERFDDDOT)/DET
      C2=(-DERF   *VALFDDDOT+VALF   *DERFDDDOT)/DET
      FDDDOT(:) = FDDDOT(:)- F(:)*C1- FDOT(:)*C2
      TFDDDOT(:)=TFDDDOT(:)-TF(:)*C1-TFDOT(:)*C2
      VALTFDDDOT=VALTFDDDOT-VALTF*C1-VALTFDOT*C2
      DERTFDDDOT=DERTFDDDOT-DERTF*C1-DERTFDOT*C2
!
!     ==  NOW MATCH DIFFERENTIABLY ONTO INPUT PARTIAL WAVE  ====================
      C1=( DERFDOT*VALPHI-VALFDOT*DERPHI)/DET
      C2=(-DERF   *VALPHI+VALF   *DERPHI)/DET
      PHI(:IRC) = F(:IRC)*C1+ FDOT(:IRC)*C2
      TPHI(:IRC)=TF(:IRC)*C1+TFDOT(:IRC)*C2
!     == DEFINE VALTPHI AND DERTPHI RELATIVE TO THE OUTER SOLUTION 
      VALTPHI   =VALTF*C1+VALTFDOT*C2  - VALTPHI
      DERTPHI   =DERTF*C1+DERTFDOT*C2  - DERTPHI
!
!     ==  NOW ENFORCE DIFFERENTIABILTY OF THE KINETIC ENERGY====================
      IF(T4PAR) THEN
        DET=VALTFDDOT*DERTFDDDOT-VALTFDDDOT*DERTFDDOT
        C3=( DERTFDDDOT*VALTPHI-VALTFDDDOT*DERTPHI)/DET
        C4=(-DERTFDDOT *VALTPHI+VALTFDDOT *DERTPHI)/DET
        PHI(:IRC) = PHI(:IRC)- FDDOT(:IRC)*C3- FDDDOT(:IRC)*C4
        TPHI(:IRC)=TPHI(:IRC)-TFDDOT(:IRC)*C3-TFDDDOT(:IRC)*C4
      ELSE 
        IF(T3PAR) THEN
          C3=TPHI(IRC)/TFDDOT(IRC)
          PHI(:IRC) = PHI(:IRC)- FDDOT(:IRC)*C3
          TPHI(:IRC)= PHI(:IRC)-TFDDOT(:IRC)*C3
        END IF
      END IF
!
!     ==========================================================================
!     == FINAL REPORTING                                                      ==
!     ==========================================================================
      IF(TWRITE) THEN
        WRITE(LSTRING,*)L
        LSTRING=ADJUSTL(LSTRING)
        CALL SETUP_WRITEPHI('PHIAFTER'//TRIM(LSTRING),GID,NR,1,PHI)
        CALL SETUP_WRITEPHI('TPHIAFTER'//TRIM(LSTRING),GID,NR,1,TPHI)
        STOP 'FORCED IN FLAT'
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_MAKEPSCORE(GID,NR,RC,L,N,PHI,TPHI)
!     **************************************************************************
!     **  REPLACES PHI BY ITS PSEUDIZED VERSION                               **
!     **                                                                      **
!     **  NOTE, THAT FDDOT IS 1/2 OF THE SECOND DERIVATIVE                    **
!     **       (H-E)|F^J>=J|F^{J-1}>.  WE DROP THE FACTOR J                   **
!     ******************************PETER BLOECHL, GOSLAR 2013 *****************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: RC
      INTEGER(4),INTENT(IN) :: L
      INTEGER(4),INTENT(IN) :: N  !WAVE FUNCTION STARTS WITH R**(L+2N)
      REAL(8)   ,INTENT(INOUT) :: PHI(NR)
      REAL(8)   ,INTENT(INOUT) :: TPHI(NR)
      LOGICAL(4),PARAMETER  :: T3PAR=.TRUE.
      LOGICAL(4),PARAMETER  :: TWRITE=.FALSE.
      INTEGER(4)            :: IRC
      REAL(8)               :: R(NR)
      REAL(8)               :: F(NR)
      REAL(8)               :: TF(NR)
      REAL(8)               :: FDOT(NR)
      REAL(8)               :: TFDOT(NR)
      REAL(8)               :: FDDOT(NR)
      REAL(8)               :: TFDDOT(NR)
      REAL(8)               :: VALPHI,VALF,VALFDOT,VALFDDOT
      REAL(8)               :: DERPHI,DERF,DERFDOT,DERFDDOT
      REAL(8)               :: C1,C2,C3,DET
      INTEGER(4)            :: IR
      INTEGER(4)            :: J
      CHARACTER(16) :: LSTRING
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRC=IR
        IF(R(IR).GT.RC) EXIT
      ENDDO
!
!     ==========================================================================
!     == INITIAL REPORTING                                                    ==
!     ==========================================================================
      IF(TWRITE) THEN
        WRITE(LSTRING,*)L
        LSTRING=ADJUSTL(LSTRING)
        CALL SETUP_WRITEPHI('PHIBEFORE'//TRIM(LSTRING),GID,NR,1,PHI)
        CALL SETUP_WRITEPHI('TPHIBEFORE'//TRIM(LSTRING),GID,NR,1,TPHI)
      END IF
!
!     ==========================================================================
!     == DETERMINE RADIAL FUNCTIONS  FOR PSEUDO PARTIAL WAVES                 ==
!     ==========================================================================
      J=L+2*N
      F=R(:)**J
      TF=-0.5D0*REAL(J*(J+1)-L*(L+1),KIND=8)*R(:)**(J-2)
      J=L+2*N+2
      FDOT=R(:)**J
      TFDOT=-0.5D0*REAL(J*(J+1)-L*(L+1),KIND=8)*R(:)**(J-2)
      J=L+2*N+4
      FDDOT=R(:)**J
      TFDDOT=-0.5D0*REAL(J*(J+1)-L*(L+1),KIND=8)*R(:)**(J-2)
!
!     ==========================================================================
!     == SUPERIMPOSE SO THAT IT MATCHES TO INPUT PARTIAL WAVE                 ==
!     ==========================================================================
      CALL RADIAL$VALUE(GID,NR,PHI,RC,VALPHI)
      CALL RADIAL$DERIVATIVE(GID,NR,PHI,RC,DERPHI)
      CALL RADIAL$VALUE(GID,NR,F,RC,VALF)
      CALL RADIAL$DERIVATIVE(GID,NR,F,RC,DERF)
      CALL RADIAL$VALUE(GID,NR,FDOT,RC,VALFDOT)
      CALL RADIAL$DERIVATIVE(GID,NR,FDOT,RC,DERFDOT)
      CALL RADIAL$VALUE(GID,NR,FDDOT,RC,VALFDDOT)
      CALL RADIAL$DERIVATIVE(GID,NR,FDDOT,RC,DERFDDOT)
      DET=VALF*DERFDOT-VALFDOT*DERF
!
!     ==  FIRST REMOVE VALUE AND DERIVATIVE FROM THE FDDOT FUNCTION ============
      C1=(DERFDOT*VALFDDOT-VALFDOT*DERFDDOT)/DET
      C2=(-DERF*VALFDDOT+VALF*DERFDDOT)/DET
      FDDOT(:) = FDDOT(:)- F(:)*C1- FDOT(:)*C2
      TFDDOT(:)=TFDDOT(:)-TF(:)*C1-TFDOT(:)*C2
!
!     ==  NOW MATCH ONTO INPUT PARTIAL WAVE  ===================================
      C1=( DERFDOT*VALPHI-VALFDOT*DERPHI)/DET
      C2=(-DERF   *VALPHI+VALF   *DERPHI)/DET
      IF(T3PAR) THEN
        C3=(TPHI(IRC)-TF(IRC)*C1-TFDOT(IRC)*C2)/TFDDOT(IRC)
      ELSE
        C3=0.D0
      END IF
      PHI(:IRC) = F(:IRC)*C1+ FDOT(:IRC)*C2+ FDDOT(:IRC)*C3
      TPHI(:IRC)=TF(:IRC)*C1+TFDOT(:IRC)*C2+TFDDOT(:IRC)*C3
!
!     ==========================================================================
!     == FINAL REPORTING                                                      ==
!     ==========================================================================
      IF(TWRITE) THEN
        WRITE(LSTRING,*)L
        LSTRING=ADJUSTL(LSTRING)
        CALL SETUP_WRITEPHI('PHIAFTER'//TRIM(LSTRING),GID,NR,1,PHI)
        CALL SETUP_WRITEPHI('TPHIAFTER'//TRIM(LSTRING),GID,NR,1,TPHI)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_TEST_NEWPRO(RELTYPE,GID,NR,ROUT,L,SO,NC,NJ,RC,EOFPHI &
     &                        ,ECORE,AEPOT,PSPOT &
     &                        ,UCORE,AECORE,PSCORE,QN,AEPHI,PSPHI,QNDOT &
     &                        ,UCORESM,AECORESM,PSCORESM,QNSM,AEPHISM,PSPHISM &
     &                        ,PRO,DTKIN,DOVER)
!     **************************************************************************
!     **  THE CORE STATES SOLVE                                               **
!     **     (H-E(I))|U_I>=|U_{I-1}>   WITH |U_0>=|0>                         **
!     **  THE NODE-REDUCED VALENCE STATES OBEY                                **
!     **     (H-ENU)|QN_J>=|QN_{J-1}>   WITH |QN_0>=|U_{NC}>                  **
!     **                                                                      **
!     **  RELTYPE='NONREL', 'ZORA', 'SCALAR', 'SPINORBIT'                     **
!     **                                                                      **
!     **  FOR THE SMALL COMPONENT, WE USE THE MODEL THAT ALSO THE PSEUDO-     **
!     **  PARTIAL WAVES CARRY A SMALL COMPOMENT, WHICH IS IDENTICAL TO THAT   **
!     **  OF THE NODE-REDUCED PARTIAL WAVES QN. DTKIN AND DOVER ARE DONE      **
!     **  ON THIS BASIS, WHICH ENSURES THAT THERE ARE NO EXPOENTIALLY GROWING **
!     **  CONTRIBUTIONS IN THE INTEGRALS.                                     **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2013******************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: RELTYPE ! SELECTOR FOR RELATIVISTIC TREATMENT
      INTEGER(4),INTENT(IN) :: GID       ! GRID ID
      INTEGER(4),INTENT(IN) :: NR        ! #(GRID POINTS)
      REAL(8)   ,INTENT(IN) :: ROUT      ! RADIUS OF ATOM IN A BOX
      INTEGER(4),INTENT(IN) :: L         ! MAIN ANGULAR MOMENTUM
      INTEGER(4),INTENT(IN) :: SO        ! SPIN ORBIT ALLIGNMENT (-1,0,1)
      INTEGER(4),INTENT(IN) :: NC        ! #(CORE STATES)
      INTEGER(4),INTENT(IN) :: NJ        ! #(PARTIAL WAVES)
      REAL(8)   ,INTENT(IN) :: RC        ! PSEUDIZATION RADIUS
      REAL(8)   ,INTENT(IN) :: EOFPHI(NJ)! EXPANSION ENERGY FOR PARTIAL WAVES
      REAL(8)   ,INTENT(IN) :: ECORE(NC) ! CORE LEVEL ENERGIES
      REAL(8)   ,INTENT(IN) :: AEPOT(NR) ! ALL-ELECTRON POTENTIAL
      REAL(8)   ,INTENT(IN) :: PSPOT(NR) ! PSEUDO POTENTIAL
!
      REAL(8)   ,INTENT(IN) :: UCORE(NR,NC)  ! NODELESS CORE WAVE FUNCTIONS
      REAL(8)   ,INTENT(IN) :: AECORE(NR,NC) ! ALL-ELECTRON CORE WAVE FUNCTIONS
      REAL(8)   ,INTENT(IN) :: PSCORE(NR,NC) ! PSEUDO CORE WAVE FUNCTIONS
      REAL(8)   ,INTENT(IN) :: UCORESM(NR,NC) ! SMALL COMPONENT OF CORE STATES
      REAL(8)   ,INTENT(IN) :: AECORESM(NR,NC) ! SMALL AE CORE WAVE FUNCTIONS
      REAL(8)   ,INTENT(IN) :: PSCORESM(NR,NC) ! SMALL PS CORE WAVE FUNCTIONS
!
      REAL(8)   ,INTENT(IN) :: QN(NR,NJ)      ! NODE REDUCED PARTIAL WAVES
      REAL(8)   ,INTENT(IN) :: QNSM(NR,NJ)    ! SMALL QN
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,NJ)
      REAL(8)   ,INTENT(IN) :: AEPHISM(NR,NJ)
      REAL(8)   ,INTENT(IN) :: PSPHI(NR,NJ)
      REAL(8)   ,INTENT(IN) :: PSPHISM(NR,NJ)
!
      REAL(8)   ,INTENT(IN) :: QNDOT(NR)      
      REAL(8)   ,INTENT(IN) :: PRO(NR,NJ)
      REAL(8)   ,INTENT(IN) :: DOVER(NJ,NJ)
      REAL(8)   ,INTENT(IN) :: DTKIN(NJ,NJ)
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)               :: R(NR)
      REAL(8)               :: G(NR)
      REAL(8)               :: PHI(NR)
      REAL(8)               :: AUX(NR),AUX1(NR),SVAR,SVAR1,SVAR2,SVAR3
      REAL(8)               :: MAT(NJ,NJ)
      REAL(8)               :: E
      INTEGER(4)            :: NFIL
      CHARACTER(16)         :: LSTRING
      INTEGER(4)            :: JP,I,J
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      NFIL=11
      WRITE(LSTRING,*)L
      LSTRING=ADJUSTL(LSTRING)
      WRITE(*,FMT='(80("="),T20," SETUP_TEST_NEWPRO START  ")')
!
!     ==========================================================================
!     == WRITE DTKIN AND DOVER                                                ==
!     ==========================================================================
      WRITE(*,FMT='(80("="),T20," DTKIN  ")')
      DO JP=1,NJ
        WRITE(*,FMT='(100F25.5)')DTKIN(JP,:)
      ENDDO
!
      WRITE(*,FMT='(80("="),T20," DOVER  ")')
      DO JP=1,NJ
        WRITE(*,FMT='(100F25.5)')DOVER(JP,:)
      ENDDO
!
!     ==========================================================================
!     == REPORT <PRO|PSPHI>                                                   ==
!     ==========================================================================
      DO JP=1,NJ
        DO J=1,NJ
          AUX(:)=R(:)**2*PRO(:,J)*PSPHI(:,JP)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,MAT(J,JP))
        ENDDO
      ENDDO
      WRITE(*,FMT='(80("="),T20," <PTILDE|PSPHI>  ")')
      DO JP=1,NJ
        WRITE(*,FMT='(100F15.10)')MAT(JP,:)
      ENDDO
!
!     ==========================================================================
!     == TEST ORTHOGONALITY OF CORE WAVE FUNCTIONS
!     ==========================================================================
      WRITE(*,FMT='(80("="),T20," CHECK ORTHOGONALITY OF CORE STATES  ")')
      DO I=1,NC
        DO J=I+1,NC
          AUX(:)=R(:)**2*(AECORE(:,I)*AECORE(:,J)+AECORESM(:,I)*AECORESM(:,J))
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          AUX(:)=R(:)**2*(AECORE(:,I)**2+AECORESM(:,I)**2)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
          AUX(:)=R(:)**2*(AECORE(:,J)**2+AECORESM(:,J)**2)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR3)
          PRINT*,'OVERLAP :',I,J,SVAR1/SQRT(SVAR2*SVAR3)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == TEST ORTHOGONALITY OF VALENCE AND CORE WAVE FUNCTIONS
!     ==========================================================================
      WRITE(*,FMT='(80("="),T20," CHECK ORTHOGONALITY BETWEEN CORE STATES AND PARTIAL WAVES ")')
      DO I=1,NC
        AUX(:)=R(:)**2*(AECORE(:,I)**2+AECORESM(:,I)**2)
        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
        CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR2)
        DO J=1,NJ
          AUX(:)=R(:)**2*(AECORE(:,I)*AEPHI(:,J)+AECORESM(:,I)*AEPHISM(:,J))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR1)
          AUX(:)=R(:)**2*(AEPHI(:,J)**2+AEPHISM(:,J)**2)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,ROUT,SVAR3)
          PRINT*,'OVERLAP :',I,J,SVAR1/SQRT(SVAR2*SVAR3)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == EXLORE SIZE OF NODELESS WAVE FUNCTIONS
!     ==========================================================================
!!$      DO I=1,NC
!!$        PRINT*,'(MAX(UCORE) IB=',I,MAXVAL(UCORE(:,I)),MAXVAL(UCORESM(:,I)) &
!!$     &                      ,ECORE(I)
!!$      ENDDO
!
!     ==========================================================================
!     == SOLVE RADIAL PAW EQUATION                                            ==
!     ==========================================================================
      CALL SETUP_PAWTEST(GID,NR,L,EOFPHI,AEPOT,PSPOT &
                        ,NJ,AEPHI,PSPHI,PRO,DTKIN,DOVER)
!
!     ==========================================================================
!     == WRITE RESULT                                                         ==
!     ==========================================================================
!
!     == CUTOFF TOO LARGE AND TOO SMALL VALUES FOR PLOTTING ====================
!!$      SVAR=1.D0
!!$      QN=MAX(-SVAR,MIN(SVAR,QN))
!!$      PSPHI=MAX(-SVAR,MIN(SVAR,PSPHI))
!!$      AEPHI=MAX(-SVAR,MIN(SVAR,AEPHI))
!
!     == WRITE PARTIAL WAVES AN PROJECTORS =====================================
!!$      CALL SETUP_WRITEPHI(-'MYQN'//TRIM(LSTRING),GID,NR,NJ,QN)
!!$      CALL SETUP_WRITEPHI(-'MYQNSM'//TRIM(LSTRING),GID,NR,NJ,QNSM)
!!$      CALL SETUP_WRITEPHI(-'MYPS'//TRIM(LSTRING),GID,NR,NJ,PSPHI)
!!$      CALL SETUP_WRITEPHI(-'MYPSSM'//TRIM(LSTRING),GID,NR,NJ,PSPHISM)
!!$      CALL SETUP_WRITEPHI(-'MYAE'//TRIM(LSTRING),GID,NR,NJ,AEPHI)
!!$      CALL SETUP_WRITEPHI(-'MYAESM'//TRIM(LSTRING),GID,NR,NJ,AEPHISM)
!!$      CALL SETUP_WRITEPHI(-'MYUCORE'//TRIM(LSTRING),GID,NR,NC,UCORE)
!!$      CALL SETUP_WRITEPHI(-'MYUCORESM'//TRIM(LSTRING),GID,NR,NC,UCORESM)
!!$      CALL SETUP_WRITEPHI(-'MYPSCORE'//TRIM(LSTRING),GID,NR,NC,PSCORE)
!!$      CALL SETUP_WRITEPHI(-'MYPSCORESM'//TRIM(LSTRING),GID,NR,NC,PSCORESM)
!!$      CALL SETUP_WRITEPHI(-'MYAECORE'//TRIM(LSTRING),GID,NR,NC,AECORE)
!!$      CALL SETUP_WRITEPHI(-'MYAECORESM'//TRIM(LSTRING),GID,NR,NC,AECORESM)
!!$      CALL SETUP_WRITEPHI(-'MYPRO'//TRIM(LSTRING),GID,NR,NJ,PRO)
!!$      CALL SETUP_WRITEPHI(-'DEPOT.DAT',GID,NR,1,(PSPOT-AEPOT)*Y0)
!
!     ==  
!PRINT*,'BEFORE NEWPROANALYZE1'
!      CALL SETUP_NEWPROANALYZE1(GID,NR,L,UCORE(:,NC),AEPOT,PSPOT,ENU,NJ,QN)
!PRINT*,'BEFORE NEWPROANALYZE2'
!      CALL SETUP_NEWPROANALYZE2(GID,NR,L,NC,ECORE,UCORE)
!
      WRITE(*,FMT='(80("="),T20," SETUP_TEST_NEWPRO END  ")')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_PAWTEST(GID,NR,L,EOFPHI,AEPOT,PSPOT &
     &                         ,NPRO,AEPHI,PSPHI,PRO,DTKIN,DOVER)
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: AEPOT(NR)
      REAL(8)   ,INTENT(IN) :: PSPOT(NR)
      INTEGER(4),INTENT(IN) :: NPRO
      REAL(8)   ,INTENT(IN) :: EOFPHI(NPRO)
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,NPRO)
      REAL(8)   ,INTENT(IN) :: PSPHI(NR,NPRO)
      REAL(8)   ,INTENT(IN) :: PRO(NR,NPRO)
      REAL(8)   ,INTENT(IN) :: DTKIN(NPRO,NPRO)
      REAL(8)   ,INTENT(IN) :: DOVER(NPRO,NPRO)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)               :: DATH(NPRO,NPRO)
      REAL(8)               :: R(NR)
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)               :: G(NR)
      REAL(8)               :: PHI(NR,NPRO)
      REAL(8)               :: ARR(NR,NPRO*NPRO)
      REAL(8)               :: SVAR
      REAL(8)               :: ROUT
      REAL(8)               :: E,EIN
      CHARACTER(16)         :: LSTRING
      REAL(8)  ,PARAMETER   :: RBNDOUT=5.D0
      INTEGER(4)            :: I,J,IND
      INTEGER(4)            :: NN
      LOGICAL(4)            :: TFAIL
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      WRITE(LSTRING,*)L
      LSTRING=ADJUSTL(LSTRING)
      WRITE(*,FMT='(80("-"),T10," PAWTEST FOR L=",I3," ")')L
      WRITE(*,FMT='(A)')"TESTS IF THE PARTIAL WAVES ARE SOLUTIONS " &
     &                  //" OF THE PAW EQUATIONS"
!
!     ==========================================================================
!     ==
!     ==========================================================================
      IND=0
      DO I=1,NPRO
        DO J=1,NPRO
          IND=IND+1
          AUX=AEPOT(:)*AEPHI(:,I)*AEPHI(:,J) &
     &       -PSPOT(:)*PSPHI(:,I)*PSPHI(:,J)
          AUX=AUX*Y0*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          ARR(:,IND)=AUX1
          CALL RADIAL$VALUE(GID,NR,AUX1,RBNDOUT,SVAR)
          DATH(I,J)=DTKIN(I,J)+SVAR 
        ENDDO
      ENDDO
!!$PRINT*,'DATH ',DATH
!!$PRINT*,'DTKIN ',DTKIN
!!$      CALL SETUP_WRITEPHI(-'ARRTEST',GID,NR,NPRO**2,ARR)
!
!     ==========================================================================
!     == TEST DEVIATION OF PAW EQUATION
!     ==========================================================================
      CALL SETUP_WRITEPHI(-'PSPHI.DAT',GID,NR,NPRO,PSPHI)
      DO I=1,NPRO
!       == CONSTRUCT PSEUDO PARTIAL WAVE AUX FROM TAYLOR EXPANSION =============
        CALL SETUP_PHIOFE(GID,NR,NPRO,EOFPHI,PSPHI,EOFPHI(I),AUX)
!       == CONSTRUCT RESIDUAL OF THE PAW EQUATION AS PHI =======================
        CALL SETUP_TESTPAWEQ(GID,NR,L,PSPOT,NPRO,PRO,DATH,DOVER,EOFPHI(I) &
     &                      ,AUX,PHI(:,I))
PRINT*,'DEV FROM PAW. EQ. FOR IPRO=',I,' MAXDEV=',MAXVAL(ABS(PHI(:NR-5,I)))
      ENDDO
!      CALL SETUP_WRITEPHI(-'DEV',GID,NR,NPRO,PHI)
!
!     ==========================================================================
!     ==
!     ==========================================================================
PRINT*,'EOFPHI ',EOFPHI
      TFAIL=.FALSE.
      G(:)=0.D0
      DO I=1,NPRO
        E=EOFPHI(I)
        NN=I-1
        ROUT=R(NR-3)
        EIN=E
!        CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO,DATH,DOVER,G,PHI(:,I))
!!$CALL SETUP_WRITEPHI(-'T4.DAT',GID,NR,1,PHI(:,I))
!!$CALL SETUP_WRITEPHI(-'T5.DAT',GID,NR,1,PSPHI(:,I))
!!$PRINT*,'BEFORE ATOMLIB$PAWBOUNDSTATE'
         CALL ATOMLIB$PAWBOUNDSTATE(GID,NR,L,NN,ROUT,PSPOT,NPRO,PRO,DATH,DOVER &
     &                              ,G,E,PHI(:,I))
!!$CALL SETUP_WRITEPHI(-'T6.DAT',GID,NR,1,PHI(:,I))
!!$PRINT*,'AFTER ATOMLIB$PAWBOUNDSTATE',E,PHI(NR,I)
        WRITE(*,FMT='("IPRO=",I2," E BEFORE:",F10.5," E AFTER: ",F10.5)')I,EIN,E
        TFAIL=TFAIL.AND.(ABS(E-EIN).LT.1.D-2)
      ENDDO
!!$      CALL SETUP_WRITEPHI(-'MYPHI1OFENU'//TRIM(LSTRING),GID,NR,NPRO,PHI)
!!$      CALL SETUP_WRITEPHI(-'MYPHI2OFENU'//TRIM(LSTRING),GID,NR,NPRO,PSPHI)
      IF(TFAIL) THEN
        CALL ERROR$MSG('ENERGIES NOT CONSISTENT')
        CALL ERROR$STOP('SETUP_PAWTEST')
      END IF
      WRITE(*,FMT='(80("-"),T10," PAWTEST FOR L=",I3," DONE  ")')L
!STOP 'FORCED'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_PHITPHIOFE(GID,NR,NPRO,EOFPHI,PSPHI,TPSPHI,E,PHI,TPHI)
!     **************************************************************************
!     ** CONSTRUCT THE PARTIAL WAVES FROM ITS TAYLOR EXPANSION COEFFICIENTS   **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2013******************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID           ! GRID ID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NPRO
      REAL(8)   ,INTENT(IN) :: EOFPHI(NPRO)
      REAL(8)   ,INTENT(IN) :: PSPHI(NR,NPRO)
      REAL(8)   ,INTENT(IN) :: TPSPHI(NR,NPRO)
      REAL(8)   ,INTENT(IN) :: E
      REAL(8)   ,INTENT(OUT):: PHI(NR)
      REAL(8)   ,INTENT(OUT):: TPHI(NR)
      INTEGER(4)            :: K
      REAL(8)               :: FACTOR
!     **************************************************************************
      FACTOR=1.D0
      PHI=0.D0
      TPHI=0.D0
      DO K=1,NPRO
        PHI=PHI+PSPHI(:,K)*FACTOR
        TPHI=TPHI+TPSPHI(:,K)*FACTOR
        FACTOR=FACTOR*(EOFPHI(K)-E)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_PHIOFE(GID,NR,NPRO,EOFPHI,PSPHI,E,PHI)
!     **************************************************************************
!     ** CONSTRUCT THE PARTIAL WAVES FROM ITS TAYLOR EXPANSION COEFFICIENTS   **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2013******************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID           ! GRID ID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NPRO
      REAL(8)   ,INTENT(IN) :: EOFPHI(NPRO)
      REAL(8)   ,INTENT(IN) :: PSPHI(NR,NPRO)
      REAL(8)   ,INTENT(IN) :: E
      REAL(8)   ,INTENT(OUT):: PHI(NR)
      INTEGER(4)            :: K
      REAL(8)               :: FACTOR
!     **************************************************************************
      FACTOR=1.D0
      PHI=0.D0
      DO K=1,NPRO
        PHI=PHI+PSPHI(:,K)*FACTOR
        FACTOR=FACTOR*(EOFPHI(K)-E)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_TESTPAWEQ(GID,NR,L,PSPOT,NPRO,PRO,DH,DO,E,PSPHI,DEV)
!     **************************************************************************
!     ** CALCULATES THE RESIDUAL DEV OF THE PAW EQUATION.                     **
!     **    DEV=HTILDE-E*OTILDE|PSPHI>                                        **
!     **                                                                      **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2013******************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID           ! GRID ID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: E
      REAL(8)   ,INTENT(IN) :: PSPOT(NR)
      INTEGER(4),INTENT(IN) :: NPRO
      REAL(8)   ,INTENT(IN) :: PRO(NR,NPRO)
      REAL(8)   ,INTENT(IN) :: DH(NPRO,NPRO)
      REAL(8)   ,INTENT(IN) :: DO(NPRO,NPRO)
      REAL(8)   ,INTENT(IN) :: PSPHI(NR)
      REAL(8)   ,INTENT(OUT):: DEV(NR)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)               :: PROJ(NPRO)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: R(NR)
      INTEGER(4)            :: I,J
      REAL(8)               :: SVAR
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  CALCULATE PROJECTION <PRO|PSPHI>                                    ==
!     ==========================================================================
      DO I=1,NPRO !LOOP OVER PARTIAL WAVES
        AUX=R**2*PRO(:,I)*PSPHI(:)
        CALL RADIAL$INTEGRAL(GID,NR,AUX,PROJ(I))
      END DO
!
!     ==========================================================================
!     ==  CALCULATE KINETIC ENERGY DENSITY                                    ==
!     ==  USE DERIVATIVES CONSISTENT WITH VERLET ALGORITHM, WHICH HAS BEEN    ==
!     ==  USED FOR SOLVING THE DIFFERENTIAL EQUATION FOR PSPHI                ==
!     ==========================================================================
      CALL RADIAL$VERLETD1(GID,NR,PSPHI,AUX)
      CALL RADIAL$VERLETD2(GID,NR,PSPHI,DEV)
      DEV=DEV+2.D0*REAL(L,KIND=8)*AUX/R-REAL(L*(L+1),KIND=8)*PSPHI/R**2
      DEV(1:2)=DEV(3)
      DEV=-0.5D0*DEV   ! TURN LAPLACIAN INTO KINETIC ENERGY

!!$! THIS WAS APPARENTLY WRONG
!!$!     == 1/R PARTIAL^2_R R- L(L+1)/R^2  
!!$!     ==  =  R^L [2(L+1)/R +PARTIAL_R] PARTIAL_R R^(-L)
!!$      DEV=PSPHI(:)/R(:)**L
!!$      CALL RADIAL$DERIVE(GID,NR,DEV,AUX)
!!$      CALL RADIAL$DERIVE(GID,NR,AUX,DEV)
!!$      DEV=R**L*(DEV+REAL(2*(L+1),KIND=8)/R(:)*AUX)
!!$      DEV(1:4)=DEV(5)
!!$!     == MAKE KINETIC ENERGY OUT OF LAPLACIAN
!!$      DEV=-0.5D0*DEV

!     == ADD POTENTIAL ENERGY ==================================================
      DEV=DEV+(PSPOT(:)*Y0-E)*PSPHI(:)
!     == ADD PROJECTOR CONTRIBUTION ============================================
      DO I=1,NPRO
        SVAR=0.D0
        DO J=1,NPRO
          SVAR=SVAR+(DH(I,J)-E*DO(I,J))*PROJ(J)
        ENDDO
        DEV=DEV+PRO(:,I)*SVAR
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3 .........4.........5.........6.........7........8
      SUBROUTINE SETUP_NEWPROANALYZE1(GID,NR,L,UCORE,AEPOT,PSPOT,ENU,NJ,QN)
!     **************************************************************************
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2013******************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID           ! GRID ID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: ENU
      INTEGER(4),INTENT(IN) :: NJ
      REAL(8)   ,INTENT(IN) :: QN(NR,NJ)
      REAL(8)   ,INTENT(IN) :: AEPOT(NR)
      REAL(8)   ,INTENT(IN) :: PSPOT(NR)
      REAL(8)   ,INTENT(IN) :: UCORE(NR)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)               :: R(NR)
      REAL(8)               :: G(NR)
      REAL(8)               :: DREL(NR)
      REAL(8)               :: AUX(NR)
      INTEGER(4),PARAMETER  :: NE=100
      REAL(8)               :: EMIN=-10.D0,EMAX=10.D0,E
      REAL(8)               :: PHIOFE0(NR,NE)
      REAL(8)               :: PHIOFE1(NR,NE)
      REAL(8)               :: RDEV(NJ,NE)
      INTEGER(4)            :: IE,J,JP,IR
      INTEGER(4),PARAMETER  :: SO=0
      INTEGER(4),PARAMETER  :: IDIR=1
      INTEGER(4),PARAMETER  :: NFIL=11
      REAL(8)               :: SVAR
      REAL(8)               :: TOL
      CHARACTER(16)         :: LSTRING
      CHARACTER(128)        :: FILENAME
!     **************************************************************************
      CALL SCHROEDINGER$DREL(GID,NR,AEPOT,ENU,DREL)
      CALL RADIAL$R(GID,NR,R)
      WRITE(LSTRING,*)L
      LSTRING=ADJUSTL(LSTRING)
!
      TOL=0
      DO IR=1,NR
        SVAR=ABS(QN(IR,1))
        IF(SVAR.LT.TOL)EXIT
        IF(R(IR).GT.2.D0) EXIT
        TOL=SVAR
      ENDDO
      TOL=TOL/10.D0
!
      RDEV=0.D0
      DO IE=1,NE
        E=EMIN+(EMAX-EMIN)*REAL(IE-1,KIND=8)/REAL(NE-1,KIND=8)
!
        G=-UCORE
        CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,SO,G,L,E,IDIR &
     &                             ,PHIOFE0(:,IE))
        PHIOFE1(:,IE)=QN(:,1)
        AUX=ABS(PHIOFE1(:,IE)-PHIOFE0(:,IE))
        DO IR=1,NR
          RDEV(1,IE)=R(IR)
          IF(AUX(IR).GT.TOL) EXIT
        ENDDO
        SVAR=1.D0
        DO JP=2,NJ
          J=JP-1
          SVAR=SVAR*(E-ENU)/REAL(J,KIND=8)
          PHIOFE1(:,IE)=PHIOFE1(:,IE)+QN(:,JP)*SVAR
!
          AUX=ABS(PHIOFE1(:,IE)-PHIOFE0(:,IE))
          DO IR=1,NR
            RDEV(JP,IE)=R(IR)
            IF(AUX(IR).GT.TOL) EXIT
          ENDDO
        ENDDO
      ENDDO
!
!     == CUTT OFF TOO LARGE AND TOO SMALL VALUES FOR PLOTTING ==================
      SVAR=5.D-4
      PHIOFE0=MAX(-SVAR,MIN(SVAR,PHIOFE0))
      PHIOFE1=MAX(-SVAR,MIN(SVAR,PHIOFE1))
!
      FILENAME='PHIOFE0_'//TRIM(LSTRING)
      OPEN(NFIL,FILE=FILENAME)
      DO IR=1,NR
        WRITE(NFIL,*)R(IR),PHIOFE0(IR,:)
      ENDDO
      CLOSE(NFIL)
      FILENAME='PHIOFE1_'//TRIM(LSTRING)
      OPEN(NFIL,FILE=FILENAME)
      DO IR=1,NR
        WRITE(NFIL,*)R(IR),PHIOFE1(IR,:)
      ENDDO
      CLOSE(NFIL)

      FILENAME='RDEV_'//TRIM(LSTRING)
      OPEN(NFIL,FILE=FILENAME)
      DO IE=1,NE
        E=EMIN+(EMAX-EMIN)*REAL(IE-1,KIND=8)/REAL(NE-1,KIND=8)
        WRITE(NFIL,*)E,RDEV(:,IE)
      ENDDO
      CLOSE(NFIL)
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_NEWPROANALYZE2(GID,NR,L,NC,ECORE,UCORE)
!     **************************************************************************
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2013******************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID           ! GRID ID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: L
      INTEGER(4),INTENT(IN) :: NC
      REAL(8)   ,INTENT(IN) :: ECORE(NC)
      REAL(8)   ,INTENT(IN) :: UCORE(NR,NC)
      INTEGER(4),PARAMETER  :: NE=20
      REAL(8)               :: EMIN=-10.D0,EMAX=10.D0
      REAL(8)               :: R(NR)
      REAL(8)               :: E
      REAL(8)               :: Y(NC,NE)
      REAL(8)               :: FCORE(NR,NE)
      INTEGER(4)            :: IE,M,J
      INTEGER(4)            :: NFIL=11
      CHARACTER(32)         :: LSTRING
      CHARACTER(128)        :: FILENAME
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      WRITE(LSTRING,*)L
      LSTRING=ADJUSTL(LSTRING)
!
      DO IE=1,NE
        E=EMIN+(EMAX-EMIN)*REAL(IE-1,KIND=8)/REAL(NE-1,KIND=8)
        FCORE(:,IE)=0.D0
        DO M=1,NC
          Y(M,IE)=1.D0
          DO J=M,NC
            Y(M,IE)=Y(M,IE)/(ECORE(J)-E)
          ENDDO
          FCORE(:,IE)=FCORE(:,IE)+UCORE(:,M)*Y(M,IE)
        ENDDO
      ENDDO
!
      FILENAME='CORECOEFF_'//TRIM(LSTRING)
      OPEN(NFIL,FILE=FILENAME)
      DO IE=1,NE
        E=EMIN+(EMAX-EMIN)*REAL(IE-1,KIND=8)/REAL(NE-1,KIND=8)
        WRITE(NFIL,*)E,Y(:,IE)
      ENDDO
      CLOSE(NFIL)
!
      CALL SETUP_WRITEPHI(-'CORECONTRIB'//TRIM(LSTRING),GID,NR,NE,FCORE)
      RETURN
      END


!
!     ..........................................................................
      SUBROUTINE SETUP_IDENTIFYRMAX()
      USE PERIODICTABLE_MODULE
      USE RADIALFOCK_MODULE, ONLY: VFOCK_TYPE
      IMPLICIT NONE
      INTEGER(4),PARAMETER  :: NR=1300
      REAL(8)   ,PARAMETER  :: R1=1.D-4
      REAL(8)   ,PARAMETER  :: DEX=1.D-2
      INTEGER(4),PARAMETER  :: NBX=20
      REAL(8)   ,PARAMETER  :: RNSCORE=0.07D0 !SEE MASTERS THESIS ROBERT SCHADE
      REAL(8)   ,PARAMETER  :: RNSPHI=0.09D0  !SEE MASTERS THESIS ROBERT SCHADE
      CHARACTER(32),PARAMETER  :: KEY='START,REL,NONSO,NONZORA'
      LOGICAL(4),PARAMETER  :: TSMALL=.TRUE.
      LOGICAL(4),PARAMETER  :: TVARDREL=.TRUE.
      INTEGER(4),PARAMETER  :: SO=0
!     == KAPPA=-L-1 FOR L*S.GE.0; KAPPA=L FOR L*S<0; KAPPA=-1 FOR SO=0 =====
      REAL(8)   ,PARAMETER  :: KAPPA=-1  ! NO SPIN ORBIT
      INTEGER(4),PARAMETER  :: LX=3
      INTEGER(4),PARAMETER  :: NFILX=7
      INTEGER(4)            :: NB
      INTEGER(4)            :: LOFI(NBX)
      INTEGER(4)            :: SOFI(NBX)
      REAL(8)               :: FOFI(NBX)
      INTEGER(4)            :: NNOFI(NBX)
      REAL(8)               :: EOFI(NBX),ECORE(NBX)
      TYPE(VFOCK_TYPE)      :: VFOCK
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: IZ
      INTEGER(4)            :: GID
      INTEGER(4)            :: IB,IR,L,IND,IRARR(1),IRCL,I
      REAL(8)               :: R(NR)
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)               :: DREL(NR)
      REAL(8)               :: G(NR)
      REAL(8)               :: GSM(NR)
      REAL(8)               :: ROUT
      REAL(8)               :: RAD(NBX)
      REAL(8)               :: AEPOT(NR)
      REAL(8)  ,ALLOCATABLE :: PSI(:,:)     !(NR,NBX)
      REAL(8)  ,ALLOCATABLE :: PSISM(:,:)   !(NR,NBX)
      REAL(8)  ,ALLOCATABLE :: UCORE(:,:)   !(NR,NBX)
      REAL(8)  ,ALLOCATABLE :: UCORESM(:,:) !(NR,NBX)
      REAL(8)               :: ETOT,AEZ
      REAL(8)               :: RNS
      REAL(8)               :: E,X,Y
      CHARACTER(20)         :: NAME(NFILX)
      REAL(8)               :: SPEEDOFLIGHT,ALPHA
      REAL(8)               :: RCOV
!     **************************************************************************
      ALLOCATE(PSI(NR,NBX))
      ALLOCATE(PSISM(NR,NBX))
      ALLOCATE(UCORE(NR,NBX))
      ALLOCATE(UCORESM(NR,NBX))

      CALL CONSTANTS$GET('C',SPEEDOFLIGHT)
      ALPHA=1.D0/SPEEDOFLIGHT ! FINE STRUCTURE CONSTANT IN A.U.
      DO I=1,NFILX
        WRITE(NAME(I),*)I
        NAME(I)='MYTESY'//TRIM(ADJUSTL(NAME(I)))
        CALL FILEHANDLER$SETFILE(NAME(I),.FALSE.,TRIM(NAME(I))//'.DAT')
        CALL FILEHANDLER$SETSPECIFICATION(NAME(I),'FORM','FORMATTED')
        CALL FILEHANDLER$UNIT(NAME(I),NFIL)
        REWIND(NFIL)
      ENDDO
!
!     ==========================================================================
!     == DEFINE RADIAL GRID                                                   ==
!     ==========================================================================
      CALL RADIAL$NEW('SHLOG',GID)
      CALL RADIAL$SETR8(GID,'R1',R1)
      CALL RADIAL$SETR8(GID,'DEX',DEX)
      CALL RADIAL$SETI4(GID,'NR',NR)
      CALL RADIAL$R(GID,NR,R)
      ROUT=R(NR-3)
      DO IR=1,NR
        IRCL=IR
        IF(R(IR).GT.3.D0) EXIT
      ENDDO
!
!     ==========================================================================
!     == LOOP OVER ALL ELEMENTS                                               ==
!     ==========================================================================
      DO IZ=1,105
        CALL PERIODICTABLE$GET(IZ,'Z',AEZ)
        CALL ATOMLIB$AESCF(GID,NR,KEY,ROUT,AEZ,NBX,NB,LOFI,SOFI,FOFI,NNOFI &
    &                   ,ETOT,AEPOT,VFOCK,EOFI,PSI,PSISM)
!
!       ========================================================================
!       ==  COLLECT RADII                                                     ==
!       ========================================================================
        DO L=0,LX
          G(:)=0.D0
          GSM(:)=0.D0
          DO IB=1,NB
            IF(LOFI(IB).NE.L) CYCLE
            E=EOFI(IB)
            CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
!
!           == PREPARE INHOMOGENEITY ===========================================
            AUX=0.5D0*ALPHA*(1.D0+DREL)*GSM   
            CALL RADIAL$DERIVE(GID,NR,AUX,AUX1)
            AUX=AUX1+(1.D0-KAPPA)/R*AUX  
            AUX(1)=AUX(2)
            G=G-AUX !GSM=-FSM(IB-1)
!
!           == OBTAIN LARGE COMPONENT ==========================================
            RNS=0.D0
            IF(IB.GT.1.AND.TSMALL) RNS=RNSCORE !AVOID SPURIOUS ZEROS NEAR ORIGIN
            CALL ATOMLIB$BOUNDSTATE(GID,NR,L,SO,RNS,ROUT,TVARDREL &
     &                         ,DREL,G,0,AEPOT,E,UCORE(:,IB))
            ECORE(IB)=E
!
!           == CONSTRUCT SMALL COMPONENT =======================================
            DREL(IRCL:)=0.D0  ! MAY HAVE BEEN RESET IN ATOMLIB$BOUNDSTATE
            CALL SCHROEDINGER$SPHSMALLCOMPONENT(GID,NR,L,SO &
     &                                      ,DREL,GSM,UCORE(:,IB),UCORESM(:,IB))
            UCORESM(IRCL:,IB)=0.D0
!
!           == PROVIDE WAVE FUNCTIONS FOR THE NEXT NODELESS LEVEL ==============
            G=-UCORE(:,IB)
            GSM=-UCORESM(:,IB)
          ENDDO
        ENDDO

!       ========================================================================
!       ==  COLLECT RADII                                                     ==
!       ========================================================================
        CALL PERIODICTABLE$GET(IZ,'R(COV)',RCOV)
        DO L=0,LX
          RAD=5.1D0
          IND=L
          DO IB=1,NB
            IF(LOFI(IB).NE.L) CYCLE
            IND=IND+1
            IRARR=MAXLOC(ABS(R*UCORE(:,IB)))
            RAD(IND)=MIN(5.D0,R(IRARR(1)))
            X=AEZ
            Y=R(IRARR(1))
            IF(Y.GT.5.D0) CYCLE
            IF(Y.LT.1.D-1) CYCLE
            Y=LOG(Y)
            X=LOG(X)
            CALL FILEHANDLER$UNIT(NAME(IND),NFIL)
            WRITE(NFIL,FMT='(2F10.5)')X,Y
          ENDDO
        ENDDO
      ENDDO
      DO I=1,NFILX
        CALL FILEHANDLER$CLOSE(NAME(I))
      ENDDO
      RETURN
      END

!!$POSITION OF THE MAXIMUM RAD(Z)=EXP[A+B*LN(Z)]  0.1<Y<5.
!!$ N=1  ---     ---     0.072414   -1.0164
!!$ N=2 1.4027 -1.1004   2.3212     -1.2082
!!$ N=3 3.5632 -1.3410   4.1957     -1.4017
!!$ N=4 5.6558 -1.5993   6.2388     -1.6617
!!$ N=5 6.9208 -1.6762   7.6373     -1.7762
!!$ N=6 9.1717 -1.9976  10.068      -2.1334 
!!$ N=7 9.9144 -1.9896  10.03       -1.9469


