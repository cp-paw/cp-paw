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
!**                                                                           **
!**  AN UNSELECT FUNCTION HAS BEEN INCLUDED BUT NOT FULLY IMPLEMENTED.        **
!**  BY FORCING UNSELECT BEFORE SELECT, ONE CAN SAFEGUARD THAT A SUBROUTINE   **
!**  CHANGES THE SETTING OF A PARENT ROUTINE.                                 **
!**                                                                           **
!**                                                                           **
!**                                              P.E. BLOECHL, (1991-2008)    **
!*******************************************************************************
TYPE SETUPPARMS_TYPE
  CHARACTER(128)  :: ID 
  REAL(8)         :: POW_POT=0.D0
  REAL(8)         :: VAL0_POT
  REAL(8)         :: RC_POT
  REAL(8)         :: POW_CORE
  REAL(8)         :: VAL0_CORE
  REAL(8)         :: RC_CORE
  REAL(8),POINTER :: RCL(:)
END TYPE SETUPPARMS_TYPE
TYPE ATOMWAVES_TYPE
  INTEGER(4)         :: NB=-1
  INTEGER(4)         :: NC=-1
  INTEGER(4),POINTER :: LOFI(:)
  INTEGER(4),POINTER :: EOFI(:)
  INTEGER(4),POINTER :: FOFI(:)
  INTEGER(4),POINTER :: AEPSI(:,:)
  INTEGER(4),POINTER :: AEPOT(:)
END TYPE ATOMWAVES_TYPE
TYPE THIS_TYPE
INTEGER(4)             :: I            ! INTEGER IDENTIFIER (ISPECIES)
CHARACTER(32)          :: ID           ! IDENTIFIER (SPECIES-NAME)
INTEGER(4)             :: GID          ! GRID ID FOR R-SPACE GRID
INTEGER(4)             :: GIDG         ! GRID ID FOR G-SPACE GRID
REAL(8)                :: AEZ          ! ATOMIC NUMBER
REAL(8)                :: RCBG
REAL(8)                :: RCSM         ! GAUSSIAN DECAY FOR COMPENSATION CHARGE
INTEGER(4)             :: LX           ! HIGHEST ANGULAR MOMENTUM
INTEGER(4)             :: LNX          ! #(ORBITAL SHELLS)
INTEGER(4)             :: LMNX
INTEGER(4)             :: LMRX         ! #(ANGULAR MOMENTA FOR 1C-DENSITY)
INTEGER(4),POINTER     :: LOX(:)       !(LNX) MAIN ANGULAR MOMENTA 
INTEGER(4),POINTER     :: ISCATT(:)    !(LNX) =-1 FOR SEMI-CORE STATE
                                       !      = 0 FOR VALENCE STATE   (PHI)
                                       !      = 1 FOR 1ST SCATTERING STATE (PHIDOT)
REAL(8)   ,POINTER     :: VADD(:)      !(NR)
REAL(8)   ,POINTER     :: AECORE(:)    !(NR)  CORE ELECTRON DENSITY
REAL(8)   ,POINTER     :: PSCORE(:)    !(NR)  PSEUDIZED ELECTRON DENSITY
REAL(8)   ,POINTER     :: PRO(:,:)     !(NR,LNX)  PROJECTOR FUNCTIONS
REAL(8)   ,POINTER     :: AEPHI(:,:)   !(NR,LNX)  AE PARTIAL WAVES
REAL(8)   ,POINTER     :: PSPHI(:,:)   !(NR,LNX)  PS PARTIAL WAVES
REAL(8)   ,POINTER     :: UPHI(:,:)    !(NR,LNX)   NODELESS PARTIAL WAVES
REAL(8)   ,POINTER     :: TUPHI(:,:)   !(NR,LNX)  KINETIC ENERGY OF ABOVE
REAL(8)   ,POINTER     :: DTKIN(:,:)   !(LNX,LNX) 1C-KIN. EN. MATRIX ELEMENTS
REAL(8)   ,POINTER     :: DOVER(:,:)   !(LNX,LNX) 1C-OVERLAP MATRIX ELEMENTS
REAL(8)   ,POINTER     :: VADDOFG(:)   !(NGX)
REAL(8)   ,POINTER     :: PSCOREOFG(:) !(NGX)
REAL(8)   ,POINTER     :: VHATOFG(:)   !(NGX)
REAL(8)   ,POINTER     :: NHATPRIMEOFG(:)  !(NGX)
REAL(8)   ,POINTER     :: PROOFG(:,:)  !(NR,LNX)  
LOGICAL(4)             :: LOCORBINI=.FALSE.       ! LOCAL ORBITALS ARE INITIALIZED IF TRUE
REAL(8)                :: LOCORBRAD(4)=5.D0  ! RADIUS OF LOCAL ORBITALS
INTEGER(4)             :: LOCORBNOFL(4)=0 ! #(LOCAL ORBITALS PER L)
INTEGER(4)             :: LOCORBLNX=0     ! #(LOCAL ORBITAL-SHELLS)
INTEGER(4),POINTER     :: LOCORBLOX(:)    ! L FOR EACH LOCAL ORBITAL-SHELL
REAL(8)   ,POINTER     :: LOCORBAMAT(:,:) ! |CHI>=|PHI>*AMAT
REAL(8)   ,POINTER     :: LOCORBBMAT(:,:) ! |PSI>=|CHI>BMAT<PTILDE|PSITILDE>
! SANTOS040616 BEGIN
INTEGER(4)             :: NC     !SANTOS040616
LOGICAL(4)             :: TCORE=.FALSE. ! CORE INFORMATION PRESENT?
REAL(8)   ,POINTER     :: AEPOT(:)     !(NRX)
INTEGER(4),POINTER     :: LB(:)        !(NC)
REAL(8)   ,POINTER     :: FB(:)        !(NC)
REAL(8)   ,POINTER     :: EB(:)        !(NC)
REAL(8)   ,POINTER     :: AEPSI(:,:)   !(NRX,NC)
! SANTOS040616 END
!
REAL(8)                :: M
REAL(8)                :: ZV
REAL(8)                :: PSG2
REAL(8)                :: PSG4
CHARACTER(32)          :: SOFTCORETYPE
CHARACTER(16)          :: FILEID
TYPE(SETUPPARMS_TYPE)  :: PARMS
TYPE(ATOMWAVES_TYPE)   :: ATOM
TYPE(THIS_TYPE),POINTER:: NEXT
END TYPE THIS_TYPE

!
!INTEGER(4)             :: NR    !=250
!INTEGER(4)             :: NRX=250
!REAL(8)   ,PARAMETER   :: GMAX=30     ! EPW[RY]<GMAX**2 FOR PSI AND RHO
!INTEGER(4),PARAMETER   :: NG=256
!REAL(8)                :: G1          ! FIRST POINT ON THE RADIAL G-GRID
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
LOGICAL(4)             :: TINTERNALSETUPS=.TRUE.
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
!        IF(SELECTED) THEN
!          CALL ERROR$MSG('SAFEGUARD FUNCTION:')
!          CALL ERROR$MSG('ANOTHER SETUP IS ALREADY SELECTED:')
!          CALL ERROR$MSG('UNSELECT BEFORE SELECT')
!          CALL ERROR$STOP('SETUP$ISELECT')
!        END IF
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
!     **************************************************************************
!
!     == CHECK IF ALREADY SELECTED
      IF(ASSOCIATED(THIS)) THEN
        IF(THIS%ID.EQ.ID) RETURN
      END IF
!
!     == CHECK IF PRESENT ===
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
!     == CREATE NEW
      IF(ALLOCATED(FASTACCESS)) DEALLOCATE(FASTACCESS)
      THIS%ID    =ID
      THIS%AEZ   =0.D0
      THIS%RCBG  =0.D0
      THIS%RCSM  =0.D0
      THIS%LNX   =0
      THIS%LMNX  =0
      THIS%LMRX  =0
      THIS%PSG2  =0.D0
      THIS%PSG4  =0.D0
      THIS%SOFTCORETYPE='NONE'
      NULLIFY(THIS%LOX)     !(LNX)
      NULLIFY(THIS%VADD)    !(NRX)
      NULLIFY(THIS%AECORE)  !(NRX)
      NULLIFY(THIS%PSCORE)  !(NRX)
      NULLIFY(THIS%PRO)     !(NRX,LNX)
      NULLIFY(THIS%AEPHI)   !(NRX,LNX)
      NULLIFY(THIS%PSPHI)   !(NRX,LNX)
      NULLIFY(THIS%DTKIN)   !(LNXX,LNX)
      NULLIFY(THIS%DOVER)   !(LNXX,LNX)
      NULLIFY(THIS%VADDOFG) !(NGX)
      NULLIFY(THIS%PSCOREOFG) !(NGX)
      NULLIFY(THIS%VHATOFG) !(NGX)
      NULLIFY(THIS%NHATPRIMEOFG) !(NGX)
      NULLIFY(THIS%PROOFG)  !(NGX,LNX)
! SANTOS040616 BEGIN
      THIS%NC    =0
      THIS%TCORE=.FALSE.
      NULLIFY(THIS%AEPOT)   !(NRX)
      NULLIFY(THIS%LB)      !(NC)
      NULLIFY(THIS%FB)      !(NC)
      NULLIFY(THIS%EB)      !(NC)
      NULLIFY(THIS%AEPSI)   !(NRX,NC)
! SANTOS040616 END
      THIS%LOCORBINI=.FALSE. ! INITIALIZED?
      THIS%LOCORBRAD(:)=5.D0  ! RADIUS OF LOCAL ORBITALS
      THIS%LOCORBNOFL(:)=0 ! (4) #(LOCAL ORBITALS PER L)
      THIS%LOCORBLNX=0     ! #(LOCAL ORBITAL-SHELLS)
      NULLIFY(THIS%LOCORBLOX)  !(LOCORBLNX)  L FOR EACH LOCAL ORBITAL-SHELL
      NULLIFY(THIS%LOCORBAMAT) !(LNX,LOCORBLNX) |CHI>=|PHI>*AMAT
      NULLIFY(THIS%LOCORBBMAT) !(LOCORBLNX,LNX) |PSI>=|CHI>BMAT<PTILDE|PSITILDE>

      NULLIFY(THIS%NEXT)
      WRITE(THIS%FILEID,*)THIS%I
      THIS%FILEID='ATOM'//ADJUSTL(THIS%FILEID)
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
      SUBROUTINE SETUP$GETL4(ID,VAL)
!     **************************************************************************
!     **  COLLECTS INTERNAL DATA                                              **
!     **                                                                      **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                        **
!     **                                                                      **
!     **************************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      LOGICAL(4)  ,INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'INTERNALSETUPS') THEN
        VAL=tinternalsetups
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('SETUP$GETL4')
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
      ELSE IF(ID.EQ.'LNXCHI') THEN
        CALL ERROR$MSG('ID LNXCHI IS MARKED FOR DELETION')
        CALL ERROR$STOP('SETUP$GETI4')
        CALL SETUP_RENEWLOCORB()
        VAL=THIS%LOCORBLNX    ! #(LOCAL ORBITALS)
! SANTOS040616 BEGIN
      ELSE IF(ID.EQ.'NC') THEN
        VAL=THIS%NC
! SANTOS040616 END
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$SETI4A(ID,LEN,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      INTEGER(4)  ,INTENT(IN)  :: VAL(LEN)
      INTEGER(4)               :: I
!     ******************************************************************
      IF(ID.EQ.'NOFLCHI') THEN
        IF(LEN.NE.4) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$SETI4A')
        END IF
!       == CHECK IF VALUE HAS CHANGED ==================================        
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
!     ..................................................................
      SUBROUTINE SETUP$GETI4A(ID,LEN,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      INTEGER(4)  ,INTENT(OUT) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'LOX') THEN
        IF(LEN.NE.THIS%LNX) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        VAL=THIS%LOX
      ELSE IF(ID.EQ.'ISCATT') THEN
        IF((LEN.NE.THIS%LNX).OR.(THIS%LNX.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        IF(.NOT.ASSOCIATED(THIS%ISCATT)) THEN
          CALL ERROR$MSG('DATA NOT AVAILABLE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        VAL=THIS%ISCATT
! SANTOS040616 BEGIN
      ELSE IF(ID.EQ.'LB') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        VAL=THIS%LB
! SANTOS040616 END
      ELSE IF(ID.EQ.'LOFC') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETI4A')
        END IF
        VAL=THIS%LB
      ELSE IF(ID.EQ.'LOXCHI') THEN
        CALL ERROR$MSG('ID LoxCHI IS MARKED FOR DELETION')
        CALL ERROR$STOP('SETUP$GETI4a')
!!$!       == ANULAR MOMENTA OF LOCAL ORBITALS ==========================
!!$        CALL SETUP_RENEWLOCORB()
!!$        IF((LEN.NE.THIS%LOCORBLNX).OR.(THIS%LOCORBLNX.EQ.0)) THEN
!!$          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
!!$          CALL ERROR$CHVAL('ID',ID)
!!$          CALL ERROR$STOP('SETUP$GETRI4A')
!!$        END IF
!!$        VAL=THIS%LOCORBLOX
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETI4A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$SETR8(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(IN)  :: VAL
!     ******************************************************************
      IF(ID.EQ.'RADCHI') THEN
        CALL ERROR$MSG('ID radchi IS MARKED FOR DELETION')
        CALL ERROR$STOP('SETUP$sETr8')
!!$!       == EXTENT OF LOCAL ORBITALS
!!$        THIS%LOCORBINI=.FALSE.
!!$        THIS%LOCORBRAD(:)=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$SETR8A(ID,LENG,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LENG
      REAL(8)     ,INTENT(IN)  :: VAL(LENG)
!     ******************************************************************
      IF(ID.EQ.'RADCHI') THEN
        CALL ERROR$MSG('ID radchi IS MARKED FOR DELETION')
        CALL ERROR$STOP('SETUP$sETr8a')
!!$!       == EXTENT OF LOCAL ORBITALS
!!$        THIS%LOCORBINI=.FALSE.
!!$        IF(LENG.GT.4) THEN
!!$          CALL ERROR$MSG('INCONSISTENT SIZE')
!!$          CALL ERROR$MSG('DIMENSION OF ARRAY MUST NOT BE LARGER THAN 4')
!!$          CALL ERROR$I4VAL('LENG',LENG)
!!$          CALL ERROR$CHVAL('ID',ID)
!!$          CALL ERROR$STOP('SETUP$SETR8')
!!$        END IF
!!$        THIS%LOCORBRAD(:)=0.D0
!!$        THIS%LOCORBRAD(:)=VAL(1:LENG)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$GETR8(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(OUT) :: VAL
!     ******************************************************************
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
      ELSE IF(ID.EQ.'<G2>') THEN
        VAL=THIS%PSG2
      ELSE IF(ID.EQ.'<G4>') THEN
        VAL=THIS%PSG4
      ELSE IF(ID.EQ.'R1') THEN
        CALL ERROR$MSG('INTERFACE MARKED FOR DELETION!')
        CALL ERROR$MSG('DEX IS OWNED BY RADIAL OBJECT')
        CALL ERROR$STOP('SETUP$GETR8')
!       VAL=R1
      ELSE IF(ID.EQ.'DEX') THEN
        CALL ERROR$MSG('INTERFACE MARKED FOR DELETION!')
        CALL ERROR$MSG('DEX IS OWNED BY RADIAL OBJECT')
        CALL ERROR$STOP('SETUP$GETR8')
!       VAL=DEX
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
      INTEGER(4)               :: NR
      INTEGER(4)               :: LRHOX
!     **************************************************************************
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
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
!     ==                                                                      ==
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
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NDLSPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%UPHI,(/LEN/))
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'NDLSTPHI') THEN
!       == DO NOT FORGET TO CLEAN UP THIS%TUPHI!!!
        CALL ERROR$MSG('OPTION IS MARKED FOR DELETION')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETR8A')
!
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%TUPHI,(/LEN/))
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
!     ==  CORE VALENCE EXCHANGE MATRIX ELEMENTS                               ==
!     ==========================================================================
      ELSE IF(ID.EQ.'CVX') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        IF(.NOT.THIS%TCORE) THEN
          CALL ERROR$MSG('CORE STATES NOT AVAILABLE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
!
!       == CALCULATE CORE-VALENCE EXCHANGE MATRIX ELEMENTS =====================
        LRHOX=INT(SQRT(REAL(THIS%LMRX-1)+1.D-8))
        CALL SETUP_CVXMAT(THIS%GID,NR,THIS%LNX,THIS%LOX,THIS%AEPHI &
     &                   ,THIS%NC,THIS%LB,THIS%AEPSI,LRHOX,VAL)
!
! SANTOS040616 BEGIN
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'FB') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%FB,(/LEN/))
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'EB') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%EB,(/LEN/))
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AECOREPSI') THEN
        IF(LEN.NE.THIS%NC*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%AEPSI,(/LEN/))
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'ATOMICAEPOT') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL(:)=THIS%AEPOT(:)
! SANTOS040616 END
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'FOFC') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('THIS%NC',THIS%NC)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%FB,(/LEN/))
!     == ATOMIC ENERGIES OF THE CORE STATES ===========================
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'EOFC') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%EB,(/LEN/))
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
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AMATCHI') THEN
        CALL ERROR$MSG('ID amatchi IS MARKED FOR DELETION')
        CALL ERROR$STOP('SETUP$GETr8a')
!!$        CALL SETUP_RENEWLOCORB()
!!$        IF(LEN.NE.THIS%LOCORBLNX*THIS%LNX) THEN
!!$          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
!!$          CALL ERROR$CHVAL('ID',ID)
!!$          CALL ERROR$STOP('SETUP$GETR8A')
!!$        END IF
!!$        VAL=RESHAPE(THIS%LOCORBAMAT,(/LEN/))
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ELSE IF(ID.EQ.'BMATCHI') THEN
        CALL ERROR$MSG('ID Bmatchi IS MARKED FOR DELETION')
        CALL ERROR$STOP('SETUP$GETr8a')
!!$        CALL SETUP_RENEWLOCORB()
!!$        IF(LEN.NE.THIS%LOCORBLNX*THIS%LNX) THEN
!!$          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
!!$          CALL ERROR$CHVAL('ID',ID)
!!$          CALL ERROR$STOP('SETUP$GETR8A')
!!$        END IF
!!$        VAL=RESHAPE(THIS%LOCORBBMAT,(/LEN/))
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
      REAL(8)                  :: PI
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
        PI=4.D0*ATAN(1.D0)
        IF(TDER) THEN
          NGAMMA=0.D0
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_READ
!     ******************************************************************
!     **  READ SELECTED SETUP                                         **
!     **  REQUIRES INFORMATION FROM ATOMTYPELIST                      **
!     **    NAME; LRHOX                                               **
!     **  REQUIRES THE FILEHANDLER TO KNOW THE SETUP FILE             **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      REAL(8)   ,PARAMETER  :: TOL=1.D-6
      INTEGER(4)            :: LMRXCUT
      INTEGER(4)            :: GID
      INTEGER(4)            :: GIDG
      INTEGER(4)            :: NFIL
      REAL(8)               :: G1
      REAL(8)   ,PARAMETER  :: GMAX=30     ! EPW[RY]<GMAX**2 FOR PSI AND RHO
      INTEGER(4),PARAMETER  :: NG=250
      REAL(8)               :: R1,DEX
      INTEGER(4)            :: NR,NRX
      INTEGER(4)            :: IR
      INTEGER(4)            :: LN
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: IRMAX
      INTEGER(4)            :: L,LX,ISVAR,LNOLD,LNX
      INTEGER(4)            :: NC !SANTOS040616
      INTEGER(4)            :: LN1,LN2,LN1A,LN2A
      INTEGER(4),ALLOCATABLE:: NPRO(:)
      INTEGER(4),ALLOCATABLE:: IWORK(:)
      REAL(8)   ,ALLOCATABLE:: DWORK(:,:,:) 
      REAL(8)               :: PI,FOURPI,Y0
      LOGICAL(4)            :: TNEWFORMAT
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)               :: SVAR
      REAL(8)               :: PSZ   ! LEGACY ONLY
!     ******************************************************************
                            CALL TRACE$PUSH('SETUP_READ')
!     == certain options only work if all setups are calculated  =======
!     == internally using setup_read_new. This switch allows to ========
!     == lock options that do not work =================================
      tinternalsetups=.false.
!
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('SETUP_READ')
      END IF
!
      CALL ATOMTYPELIST$NAME(THIS%I,THIS%ID)
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETR8('M',THIS%M)
      CALL ATOMTYPELIST$GETR8('ZV',THIS%ZV)
      CALL ATOMTYPELIST$GETR8('PS<G2>',THIS%PSG2)
      CALL ATOMTYPELIST$GETR8('PS<G4>',THIS%PSG4)
      CALL ATOMTYPELIST$GETCH('SOFTCORETYPE',THIS%SOFTCORETYPE)
!
      CALL FILEHANDLER$UNIT(THIS%FILEID,NFIL)
      CALL SETUPREAD$NEW(NFIL,TNEWFORMAT)
!
      IF(TNEWFORMAT) THEN
        CALL SETUPREAD$GETI4('LNX',LNX)
        CALL SETUPREAD$GETI4('GID',THIS%GID)
        GID=THIS%GID
        CALL SETUPREAD$GETI4('NR',NR)
        NRX=NR
        THIS%LNX=LNX
      ELSE
        CALL INPOT$GRID(NFIL,R1,DEX,NR)
        CALL RADIAL$NEW('LOG',GID)
        THIS%GID=GID
        NRX=NR
        CALL RADIAL$SETI4(GID,'NR',NR)
        CALL RADIAL$SETR8(GID,'DEX',DEX)
        CALL RADIAL$SETR8(GID,'R1',R1)
        CALL INPOT$LNX(NFIL,LNX)
        THIS%LNX=LNX
      END IF
      ALLOCATE(THIS%LOX(LNX))
      ALLOCATE(THIS%VADD(NRX))
      ALLOCATE(THIS%AECORE(NRX))
      ALLOCATE(THIS%PSCORE(NRX))
      ALLOCATE(THIS%PRO(NRX,LNX))
      ALLOCATE(THIS%AEPHI(NRX,LNX))
      ALLOCATE(THIS%PSPHI(NRX,LNX))
      ALLOCATE(THIS%UPHI(NRX,LNX))
      ALLOCATE(THIS%TUPHI(NRX,LNX))
      ALLOCATE(THIS%DTKIN(LNX,LNX))
      ALLOCATE(THIS%DOVER(LNX,LNX))
      THIS%VADD=0.D0
      THIS%AECORE=0.D0
      THIS%PSCORE=0.D0
      THIS%PRO=0.D0
      THIS%AEPHI=0.D0
      THIS%PSPHI=0.D0
! SANTOS040616 BEGIN
      IF(TNEWFORMAT) THEN
        CALL SETUPREAD$GETI4('NC',NC)
      ELSE
        CALL INPOT$NC(NFIL,NC)
      END IF
      THIS%NC=NC
      ALLOCATE(THIS%AEPOT(NRX))
      ALLOCATE(THIS%LB(NC))
      ALLOCATE(THIS%FB(NC))
      ALLOCATE(THIS%EB(NC))
      ALLOCATE(THIS%AEPSI(NRX,NC))
      THIS%AEPOT=0.D0
      THIS%LB=0.D0
      THIS%FB=0.D0
      THIS%EB=0.D0
      THIS%AEPSI=0.D0
! SANTOS040616 END
!     
!     ==================================================================
!     ==  READ PSEUDOPOTENTIALS AND PSEUDO WAVE FUNCTIONS             ==
!     ==================================================================
                            CALL TRACE$PASS('READ SETUP FILES')
      THIS%RCBG=1.D0/SQRT(0.218D0)
      
!      CALL INPOT$READALL(NFIL,NRX,R1,DEX,NR,THIS%LNX,THIS%LOX &
!     &         ,THIS%AEZ,THIS%PSZ,THIS%PSPHI,THIS%AEPHI &
!     &         ,THIS%VADD,THIS%RCSM,THIS%DTKIN,THIS%DOVER &
!     &         ,IRCCOR,THIS%AECORE,THIS%PSCORE,THIS%PRO)
      IF(TNEWFORMAT) THEN
        CALL SETUPREAD$GETI4A('LOX',LNX,THIS%LOX)
        CALL SETUPREAD$GETR8('AEZ',THIS%AEZ)
        CALL SETUPREAD$GETR8A('PSPHI',NR*LNX,THIS%PSPHI)
        CALL SETUPREAD$GETR8A('AEPHI',NR*LNX,THIS%AEPHI)
        CALL SETUPREAD$GETR8A('PRO',NR*LNX,THIS%PRO)
        CALL SETUPREAD$GETR8A('NDLSPHI',NR*LNX,THIS%UPHI)
        CALL SETUPREAD$GETR8A('NDLSTPHI',NR*LNX,THIS%TUPHI)
        CALL SETUPREAD$GETR8A('VADD',NR,THIS%VADD)
        CALL SETUPREAD$GETR8('RCSM',THIS%RCSM)
        CALL SETUPREAD$GETR8A('DT',LNX*LNX,THIS%DTKIN)
        CALL SETUPREAD$GETR8A('DO',LNX*LNX,THIS%DOVER)
        CALL SETUPREAD$GETR8A('PSCORE',NR,THIS%PSCORE)
        CALL SETUPREAD$GETR8A('AECORE',NR,THIS%AECORE)
        CALL SETUPREAD$GETR8A('AEPOT',NR,THIS%AEPOT)
!ADD HERE THE CORE WAVE FUNCTIONS FOR SANTOS
        THIS%TCORE=.TRUE.
        CALL SETUPREAD$GETI4A('LOFC',NC,THIS%LB)
        CALL SETUPREAD$GETR8A('FOFC',NC,THIS%FB)
        CALL SETUPREAD$GETR8A('EOFC',NC,THIS%EB)
        CALL SETUPREAD$GETR8A('AEPSICORE',NR*NC,THIS%AEPSI)
      ELSE
        CALL INPOT$READALL(NFIL,NRX,R1,DEX,NR,THIS%LNX,THIS%LOX &
     &         ,THIS%AEZ,PSZ,THIS%PSPHI,THIS%AEPHI &
     &         ,THIS%VADD,THIS%RCSM,THIS%DTKIN,THIS%DOVER &
     &         ,THIS%AECORE,THIS%PSCORE,THIS%PRO &
     &         ,THIS%TCORE,THIS%AEPOT,THIS%NC,THIS%LB,THIS%FB,THIS%EB,THIS%AEPSI) !SANTOS040616
      END IF
      CALL FILEHANDLER$CLOSE(THIS%FILEID)
PRINT*,'NEW FORMAT?',TNEWFORMAT
ALLOCATE(R(NR))
CALL RADIAL$R(GID,NR,R)
CALL RADIAL$INTEGRAL(GID,NR,4.D0*PI*THIS%AECORE*Y0*R**2,SVAR)
PRINT*,'INT AECORE ',SVAR
CALL RADIAL$INTEGRAL(GID,NR,4.D0*PI*THIS%PSCORE*Y0*R**2,SVAR)
PRINT*,'INT PSCORE ',SVAR
PRINT*,'AEZ ',THIS%AEZ
PRINT*,'RCSM ',THIS%RCSM
!STOP
!     
!     ==================================================================
!     == LIMIT NUMBER OF PROJECTORS FOR EACH L                        ==
!     ==================================================================
      LX=0
      DO LN=1,THIS%LNX
        LX=MAX(LX,THIS%LOX(LN))
      ENDDO
      ALLOCATE(NPRO(LX+1)) 
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETI4A('NPRO',LX+1,NPRO)
      DO L=0,LX
        ISVAR=0
        DO LN=1,THIS%LNX
          IF(THIS%LOX(LN).NE.L)CYCLE
          ISVAR=ISVAR+1
          IF(ISVAR.GT.NPRO(L+1)) THEN
            THIS%LOX(LN)=-1     ! MARK PROJECTORS TO BE DELETED BY LOX=-1
            ISVAR=ISVAR-1
          END IF
        ENDDO
        NPRO(L+1)=ISVAR
      ENDDO
      CALL ATOMTYPELIST$SETI4A('NPRO',LX+1,NPRO)
      DEALLOCATE(NPRO)
!
      LNOLD=THIS%LNX
      LNX=0
      DO LN=1,LNOLD
        IF(THIS%LOX(LN).NE.-1) LNX=LNX+1
      ENDDO
      THIS%LNX=LNX
!
!     == FOLD DOWN ARRAYS FOR PROJECTORS AND PARTIALWAVES, LOX =========
      ALLOCATE(DWORK(NRX,LNOLD,5))
      ALLOCATE(IWORK(LNOLD))
      DWORK(:,:,1)=THIS%PRO(:,:)
      DWORK(:,:,2)=THIS%AEPHI(:,:)
      DWORK(:,:,3)=THIS%PSPHI(:,:)
      DWORK(:,:,4)=THIS%UPHI(:,:)
      DWORK(:,:,5)=THIS%TUPHI(:,:)
      IWORK(:)=THIS%LOX(:)
      DEALLOCATE(THIS%PRO)
      DEALLOCATE(THIS%AEPHI)
      DEALLOCATE(THIS%PSPHI)
      DEALLOCATE(THIS%UPHI)
      DEALLOCATE(THIS%TUPHI)
      DEALLOCATE(THIS%LOX)
      ALLOCATE(THIS%PRO(NRX,THIS%LNX))
      ALLOCATE(THIS%AEPHI(NRX,THIS%LNX))
      ALLOCATE(THIS%PSPHI(NRX,THIS%LNX))
      ALLOCATE(THIS%UPHI(NRX,THIS%LNX))
      ALLOCATE(THIS%TUPHI(NRX,THIS%LNX))
      ALLOCATE(THIS%LOX(THIS%LNX))
      ISVAR=0
      DO LN=1,LNOLD
        IF(IWORK(LN).EQ.-1) CYCLE
        ISVAR=ISVAR+1
        THIS%PRO(:,ISVAR)=DWORK(:,LN,1)
        THIS%AEPHI(:,ISVAR)=DWORK(:,LN,2)
        THIS%PSPHI(:,ISVAR)=DWORK(:,LN,3)
        THIS%UPHI(:,ISVAR)=DWORK(:,LN,4)
        THIS%TUPHI(:,ISVAR)=DWORK(:,LN,5)
        THIS%LOX(ISVAR)=IWORK(LN)
      ENDDO
      DEALLOCATE(DWORK)
!
!     == FOLD DOWN ARRAYS FOR DTKIN AND DOVER ==========================
      ALLOCATE(DWORK(LNOLD,LNOLD,2))
      DWORK(:,:,1)=THIS%DTKIN(:,:)
      DWORK(:,:,2)=THIS%DOVER(:,:)
      DEALLOCATE(THIS%DTKIN)
      DEALLOCATE(THIS%DOVER)
      ALLOCATE(THIS%DTKIN(LNX,LNX))
      ALLOCATE(THIS%DOVER(LNX,LNX))
      LN1A=0
      DO LN1=1,LNOLD
        IF(IWORK(LN1).EQ.-1) CYCLE
        LN1A=LN1A+1
        LN2A=0
        DO LN2=1,LNOLD
          IF(IWORK(LN2).EQ.-1) CYCLE
          LN2A=LN2A+1
          THIS%DTKIN(LN1A,LN2A)=DWORK(LN1,LN2,1)
          THIS%DOVER(LN1A,LN2A)=DWORK(LN1,LN2,2)
        ENDDO
      ENDDO
      DEALLOCATE(DWORK)
!
      DEALLOCATE(IWORK)

!!$PRINT*,'=============================================================='
!!$PRINT*,'AEZ  ',THIS%AEZ
!!$PRINT*,'RCSM ',THIS%RCSM
!!$PRINT*,'LNX  ',THIS%LNX
!!$PRINT*,'LOX  ',THIS%LOX
!!$PRINT*,'PRO  ',THIS%PRO
!!$PRINT*,'AEPHI',THIS%AEPHI
!!$PRINT*,'PSPHI',THIS%PSPHI
!!$PRINT*,'VADD ',THIS%VADD
!!$PRINT*,'PSCORE',THIS%PSCORE
!!$PRINT*,'AECORE',THIS%AECORE
!!$PRINT*,'DOVER',THIS%DOVER
!!$PRINT*,'DTKIN',THIS%DTKIN
!CALL ERROR$STOP('FORCED STOP IN SETUP')
!     
!     ==================================================================
!     == SET VALUES BEYOND A CERTAIN RADIUS EXACTLY TO ZERO           ==
!     ==================================================================
                            CALL TRACE$PASS('CHECK MAX. RADIUS')
      CALL RADIAL$R(GID,NR,R)
      IRMAX=0
      DO IR=1,NR
        TCHK=(ABS(THIS%VADD(IR)).LT.TOL)
        TCHK=TCHK.AND.(ABS(THIS%PSCORE(IR)-THIS%AECORE(IR)).LT.TOL)
        DO LN=1,THIS%LNX
          TCHK=TCHK.AND. &
     &           (ABS(THIS%AEPHI(IR,LN)-THIS%PSPHI(IR,LN)).LT.TOL)
        ENDDO
!       == LDAPLUSU REQUIRES A SOMEWHAT LARGER RADIUS ==================
        TCHK=TCHK.AND.(R(IR).GE.6.D0)  
        IF(.NOT.TCHK) IRMAX=IR
      ENDDO
      DO IR=IRMAX+1,NR
        THIS%VADD(IR)=0.D0
        DO LN=1,THIS%LNX
          THIS%AEPHI(IR,LN)=0.D0
          THIS%PSPHI(IR,LN)=0.D0
          THIS%UPHI(IR,LN)=0.D0
          THIS%TUPHI(IR,LN)=0.D0
        ENDDO
      ENDDO
!     
!     ================================================================
!     ==  DEFINE ARRAYS                                             ==
!     ================================================================
                            CALL TRACE$PASS('DEFINE ARRAYS')
!
!     == SELECT NATURAL VALUES =======================================
      THIS%LMNX=0
      THIS%LMRX=0
      DO LN=1,THIS%LNX
        THIS%LMNX=THIS%LMNX+2*THIS%LOX(LN)+1
        THIS%LMRX=MAX(THIS%LMRX,(2*THIS%LOX(LN)+1)**2)
      ENDDO
!
!     == LIMIT MAX ANGULAR MOMENTUM FOR THE DENSITY TO MAX VALUE =======
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETI4('LRHOX',LMRXCUT)
      LMRXCUT=(LMRXCUT+1)**2
      THIS%LMRX=MIN(THIS%LMRX,LMRXCUT)
      CALL ATOMTYPELIST$UNSELECT
!     
!     ==================================================================
!     ==  UPDATE GLOBAL VARIABLES                                     ==
!     ==================================================================
      LMNXX=MAX(LMNXX,THIS%LMNX)
      LMRXX=MAX(LMRXX,THIS%LMRX)
      LNXX=MAX(LNXX,THIS%LNX)
!     
!     ==================================================================
!     ==  PERFORM BESSELTRANSFORMS                                    ==
!     ==================================================================
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      CALL RADIAL$GETR8(GID,'DEX',DEX)
      G1=GMAX*EXP(-DEX*DBLE(NG-1))
      CALL RADIAL$NEW('LOG',GIDG)
      THIS%GIDG=GIDG
      CALL RADIAL$SETI4(GIDG,'NR',NG)
      CALL RADIAL$SETR8(GIDG,'R1',G1)
      CALL RADIAL$SETR8(GIDG,'DEX',DEX)
PRINT*,'GIDG ',GIDG,G1,DEX,NG
!       
!     == VADD (VBAR) ===================================================
      ALLOCATE(THIS%VADDOFG(NG))
      CALL RADIAL$BESSELTRANSFORM(0,GID,NR,THIS%VADD,GIDG,NG,THIS%VADDOFG)
      THIS%VADDOFG(:)=FOURPI*THIS%VADDOFG(:)
!     == PSCORE (VBAR) =================================================
      ALLOCATE(THIS%PSCOREOFG(NG))
      CALL RADIAL$BESSELTRANSFORM(0,GID,NR,THIS%PSCORE,GIDG,NG,THIS%PSCOREOFG)
      THIS%PSCOREOFG(:)=FOURPI*THIS%PSCOREOFG(:)
!     == PROJECTORS ====================================================
      ALLOCATE(THIS%PROOFG(NG,LNX))
      DO LN=1,LNX
        L=THIS%LOX(LN)
        CALL RADIAL$BESSELTRANSFORM(L,GID,NR,THIS%PRO(:,LN),GIDG,NG,THIS%PROOFG(:,LN))
        THIS%PROOFG(:,LN)=FOURPI*THIS%PROOFG(:,LN)
      ENDDO
!     == COMPENSATION GAUSSIAN =========================================
      ALLOCATE(THIS%NHATPRIMEOFG(NG))
      ALLOCATE(THIS%VHATOFG(NG))
      CALL SETUP_COMPOFG(THIS%RCBG,THIS%RCSM,GIDG,NG &
     &                  ,THIS%NHATPRIMEOFG,THIS%VHATOFG)
!      
                            CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_READ_NEW
!     **************************************************************************
!     **  READ SELECTED SETUP                                                 **
!     **  REQUIRES INFORMATION FROM ATOMTYPELIST                              **
!     **    NAME; LRHOX                                                       **
!     **  REQUIRES THE FILEHANDLER TO KNOW THE SETUP FILE                     **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER  :: NBX=19
      INTEGER(4)            :: LOFI(NBX)
      INTEGER(4)            :: SOFI(NBX)
      REAL(8)               :: FOFI(NBX)
      INTEGER(4)            :: NNOFI(NBX)
      REAL(8)               :: EOFI(NBX)
      INTEGER(4)            :: GID
      INTEGER(4)            :: GIDG
      REAL(8)               :: DEX
      REAL(8)               :: G1
      REAL(8)   ,PARAMETER  :: GMAX=30     ! EPW[RY]<GMAX**2 FOR PSI AND RHO
      INTEGER(4),PARAMETER  :: NG=250
      INTEGER(4)            :: NR
      INTEGER(4)            :: LX,LNX
      INTEGER(4)            :: NB
      INTEGER(4)            :: NC 
      INTEGER(4),ALLOCATABLE:: NPRO(:)
      REAL(8)   ,ALLOCATABLE:: R(:)
      INTEGER(4),ALLOCATABLE:: LOX(:)
      REAL(8)   ,ALLOCATABLE:: RC(:)
      REAL(8)               :: RBOX
      REAL(8)               :: AEZ
      REAL(8)               :: ZV
      REAL(8)               :: ETOT
      CHARACTER(32)         :: KEY
      REAL(8)   ,ALLOCATABLE:: PSI(:,:)
      LOGICAL(4)            :: TCHK
      REAL(8)               :: PI,FOURPI,Y0,C0LL
      REAL(8)               :: SVAR
      INTEGER(4)            :: IR,IB,LN,L
!     **************************************************************************
                            CALL TRACE$PUSH('SETUP_READ_NEW')
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      C0LL=1.D0/SQRT(FOURPI)
      IF(.NOT.ASSOCIATED(THIS)) THEN
        CALL ERROR$MSG('NO SETUP SELECTED')
        CALL ERROR$STOP('SETUP_READ')
      END IF
!
      CALL ATOMTYPELIST$NAME(THIS%I,THIS%ID)
      CALL ATOMTYPELIST$SELECT(THIS%ID)
      CALL ATOMTYPELIST$GETR8('M',THIS%M)
      CALL ATOMTYPELIST$GETR8('PS<G2>',THIS%PSG2)
      CALL ATOMTYPELIST$GETR8('PS<G4>',THIS%PSG4)
      CALL ATOMTYPELIST$GETCH('SOFTCORETYPE',THIS%SOFTCORETYPE)
!
!     ==========================================================================
!     == EXTRACT LNX,LOX,LX FROM NPRO AS DEFINED IN STRC INPUT FILE           ==
!     ==========================================================================
      ALLOCATE(NPRO(10)) 
      CALL ATOMTYPELIST$GETI4A('NPRO',10,NPRO)
      LNX=SUM(NPRO)
      ALLOCATE(LOX(LNX))
      LN=0
      DO L=0,9
        DO WHILE (NPRO(L+1).GT.0)
          LN=LN+1
          NPRO(L+1)=NPRO(L+1)-1
          LOX(LN)=L
        ENDDO
      ENDDO
      LX=MAXVAL(LOX)
      DEALLOCATE(NPRO)
      THIS%LMNX=SUM(2*LOX(:)+1)
      THIS%LNX=LNX
      ALLOCATE(THIS%LOX(LNX))
      THIS%LOX(:)=LOX(:LNX)
!
!     == LIMIT LRHOX TO THE MAXIMUM CONSISTENT WITH WAVE FUNCTION CUTOFF
      CALL ATOMTYPELIST$GETI4('LRHOX',THIS%LMRX)
      THIS%LMRX=(THIS%LMRX+1)**2
      THIS%LMRX=MIN(THIS%LMRX,(2*LX+1)**2)
      CALL ATOMTYPELIST$SETI4('LRHOX',THIS%LMRX)
!     
!     ==  UPDATE GLOBAL VARIABLES ==============================================
      LMNXX=MAX(LMNXX,THIS%LMNX)
      LMRXX=MAX(LMRXX,THIS%LMRX)
      LNXX=MAX(LNXX,THIS%LNX)
!
!     ==========================================================================
!     == READ SETUP INFORMATION FROM PARAMETER FILE                           ==
!     ==========================================================================
      ALLOCATE(THIS%PARMS%RCL(LX+1))
      CALL ATOMTYPELIST$GETCH('ID',THIS%PARMS%ID)
      CALL ATOMLIB$SCNTLLOOKUP(THIS%PARMS%ID,GID,AEZ,ZV,RBOX,LX,THIS%PARMS%RCL,THIS%RCSM &
     &                     ,THIS%PARMS%POW_POT,THIS%PARMS%RC_POT,THIS%PARMS%VAL0_POT &
     &                     ,THIS%PARMS%POW_CORE,THIS%PARMS%RC_CORE,THIS%PARMS%VAL0_CORE)
      THIS%GID=GID
      THIS%AEZ=AEZ
      THIS%ZV=ZV
      THIS%RCBG=1.D0/SQRT(0.218D0)
      CALL RADIAL$GETI4(GID,'NR',NR)
!     
!     ==========================================================================
!     ==  PERFORM ALL-ELECTRON CALCULATION FOR THE ATOM IN A BOX              ==
!     ==========================================================================
      ALLOCATE(THIS%AEPOT(NR))
      ALLOCATE(PSI(NR,NBX))
      KEY='START,REL,NONSO'
      CALL ATOMLIB$AESCF(GID,NR,KEY,RBOX,AEZ,NBX,NB,LOFI,SOFI,FOFI,NNOFI &
    &                   ,ETOT,THIS%AEPOT,EOFI,PSI)
!
!     == DETERMINE NUMBER OF CORE SHELLS =======================================
      SVAR=AEZ-ZV
      NC=0
      if(svar.gt.1.d-3) then
        DO IB=1,NB
          SVAR=SVAR-FOFI(IB)
          NC=IB
          IF(ABS(SVAR).LT.1.D-3) EXIT
          IF(SVAR.LT.-1.D-3) THEN
            CALL ERROR$MSG('ZV INCONSISTENT WITH COMPLETE ANGULAR MOMENTUM SHELLS')
            CALL ERROR$R8VAL('AEZ',AEZ)
            CALL ERROR$R8VAL('ZV',ZV)
            CALL ERROR$R8VAL('svar',svar)
            CALL ERROR$STOP('SETUP_READ_NEW')
          END IF
        ENDDO
      end if
!
!     ==========================================================================
      ALLOCATE(THIS%ISCATT(LNX))
      CALL SETUP_MAKEISCATT(AEZ,NB,NC,LOFI(1:NB),NNOFI(1:NB),LNX,LOX,THIS%ISCATT)
PRINT*,'AEZ',THIS%AEZ 
PRINT*,'LOX',THIS%LOX 
PRINT*,'ISCATT',THIS%ISCATT 
!
!     == MAP CORE ORBITALS ON GRID =============================================
!     == THIS IS REDUNDANT AND SHALL BE REMOVED LATER (SEE THIS%ATOM BELOW)      
      THIS%NC=NC
      THIS%TCORE=.TRUE.
      ALLOCATE(THIS%LB(NC))
      ALLOCATE(THIS%FB(NC))
      ALLOCATE(THIS%EB(NC))
      ALLOCATE(THIS%AEPSI(NR,NC))
      THIS%LB=LOFI(1:NC)
      THIS%FB=FOFI(1:NC)
      THIS%EB=EOFI(1:NC)
      THIS%AEPSI=PSI(:,1:NC)
!
!     == MAP ATOMIC DATA ON GRID ================================================
      THIS%ATOM%NB=NB
      THIS%ATOM%NC=NC
      ALLOCATE(THIS%ATOM%LOFI(NB))
      ALLOCATE(THIS%ATOM%FOFI(NB))
      ALLOCATE(THIS%ATOM%EOFI(NB))
      ALLOCATE(THIS%ATOM%AEPSI(NR,NB))
      ALLOCATE(THIS%ATOM%AEPOT(NR))
      THIS%ATOM%LOFI=LOFI(1:NB)
      THIS%ATOM%LOFI=FOFI(1:NB)
      THIS%ATOM%EOFI=EOFI(1:NB)
      THIS%ATOM%AEPSI=PSI(:,1:NB)
      THIS%ATOM%AEPOT=THIS%AEPOT
!
!     ==========================================================================
!     == CALCULATE AND PSEUDIZE CORE DENSITY                                  ==
!     ==========================================================================
      ALLOCATE(THIS%AECORE(NR))
      ALLOCATE(THIS%PSCORE(NR))
      THIS%AECORE(:)=0.D0
      DO IB=1,NC
        THIS%AECORE(:)=THIS%AECORE(:)+FOFI(IB)*PSI(:,IB)**2*C0LL
      ENDDO
      CALL ATOMIC_PSEUDIZE(GID,NR &
     &         ,THIS%PARMS%POW_CORE,THIS%PARMS%VAL0_CORE,THIS%PARMS%RC_CORE &
     &         ,THIS%AECORE,THIS%PSCORE)
      DEALLOCATE(PSI)
!
!     ==========================================================================
!     == CONSTRUCT PARTIAL WAVES                                              ==
!     ==========================================================================
      ALLOCATE(THIS%VADD(NR))
      ALLOCATE(THIS%PRO(NR,LNX))
      ALLOCATE(THIS%AEPHI(NR,LNX))
      ALLOCATE(THIS%PSPHI(NR,LNX))
      ALLOCATE(THIS%UPHI(NR,LNX))
      ALLOCATE(THIS%TUPHI(NR,LNX))
      ALLOCATE(THIS%DTKIN(LNX,LNX))
      ALLOCATE(THIS%DOVER(LNX,LNX))
      THIS%VADD=0.D0
      THIS%PRO=0.D0
      THIS%AEPHI=0.D0
      THIS%PSPHI=0.D0
!
      ALLOCATE(RC(LNX))
      DO LN=1,LNX
        L=LOX(LN)
        RC(LN)=THIS%PARMS%RCL(L+1)
      ENDDO
!
      CALL ATOMIC_MAKEPARTIALWAVES(GID,NR,KEY,AEZ,THIS%AEPOT,NB,NC &
     &           ,LOFI(1:NB),SOFI(1:NB),NNOFI(1:NB),EOFI(1:NB),FOFI(1:NB) &
     &           ,RBOX,LNX,LOX,RC,THIS%AEPHI,THIS%PSPHI,THIS%UPHI,THIS%PRO &
     &           ,THIS%DTKIN,THIS%DOVER &
     &           ,THIS%PARMS%POW_POT,THIS%PARMS%VAL0_POT,THIS%PARMS%RC_POT &
     &           ,THIS%RCSM,THIS%VADD)
!     
!     ==================================================================
!     ==  PERFORM BESSELTRANSFORMS                                    ==
!     ==================================================================
      GID=THIS%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      CALL RADIAL$GETR8(GID,'DEX',DEX)
      G1=GMAX*EXP(-DEX*DBLE(NG-1))
G1=1.175D-4     
DEX=0.05D0
      CALL RADIAL$NEW('LOG',GIDG)
      THIS%GIDG=GIDG
      CALL RADIAL$SETI4(GIDG,'NR',NG)
      CALL RADIAL$SETR8(GIDG,'R1',G1)
      CALL RADIAL$SETR8(GIDG,'DEX',DEX)
PRINT*,'GIDG ',GIDG,G1,DEX,NG
!       
!     == VADD (VBAR) ===================================================
      ALLOCATE(THIS%VADDOFG(NG))

      CALL RADIAL$BESSELTRANSFORM(0,GID,NR,THIS%VADD,GIDG,NG,THIS%VADDOFG)
      THIS%VADDOFG(:)=FOURPI*THIS%VADDOFG(:)
!     == PSCORE (VBAR) =================================================
      ALLOCATE(THIS%PSCOREOFG(NG))
      CALL RADIAL$BESSELTRANSFORM(0,GID,NR,THIS%PSCORE,GIDG,NG,THIS%PSCOREOFG)
      THIS%PSCOREOFG(:)=FOURPI*THIS%PSCOREOFG(:)
!     == PROJECTORS ====================================================
      ALLOCATE(THIS%PROOFG(NG,LNX))
      DO LN=1,LNX
        L=THIS%LOX(LN)
        CALL RADIAL$BESSELTRANSFORM(L,GID,NR,THIS%PRO(:,LN),GIDG,NG,THIS%PROOFG(:,LN))
        THIS%PROOFG(:,LN)=FOURPI*THIS%PROOFG(:,LN)
      ENDDO
!     == COMPENSATION GAUSSIAN =========================================
      ALLOCATE(THIS%NHATPRIMEOFG(NG))
      ALLOCATE(THIS%VHATOFG(NG))
      CALL SETUP_COMPOFG(THIS%RCBG,THIS%RCSM,GIDG,NG &
     &                  ,THIS%NHATPRIMEOFG,THIS%VHATOFG)
!      
                            CALL TRACE$POP
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
      INTEGER(4)                 :: L,NPRO,LN,NPROSUM
      TYPE(THIS_TYPE),POINTER    :: THIS1
      CHARACTER(64)              :: STRING
      REAL(8)                    :: U,PI,Y0
      INTEGER(4)                 :: THISTASK,NTASKS
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL CONSTANTS$GET('U',U)
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
             STRING='#(PROJECTORS FOR L='//ADJUSTL(STRING)
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
        CALL REPORT$R8VAL(NFIL,'GAUSSIAN DECAY FOR EXTENDED COMPENSATION DENSITY' &
     &                        ,THIS1%RCBG,'ABOHR ')
        CALL REPORT$R8VAL(NFIL,'PARAMETER PS<G2> FOR MASS RENORMALIZATION' &
     &                        ,THIS1%PSG2,'H')
        CALL REPORT$R8VAL(NFIL,'PARAMETER PS<G4> FOR MASS RENORMALIZATION' &
     &                        ,THIS1%PSG4,'H')
!
!       == REPORT SETTINGS FOR THE SETUP CONSTRUCTION ==========================
        IF(THIS1%PARMS%POW_POT.NE.0.D0) THEN
          CALL REPORT$CHVAL(NFIL,'SETUP ID',THIS1%PARMS%ID)
          DO L=0,MAXVAL(THIS1%LOX)
            WRITE(STRING,*)L
            STRING='PARTIAL WAVE PSEUDIZATION PARAMETER RC FOR L='//ADJUSTL(STRING)
            CALL REPORT$R8VAL(NFIL,TRIM(STRING),THIS1%PARMS%RCL(L+1),'ABOHR')
          ENDDO
          CALL REPORT$R8VAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: VAL(0)' &
    &                            ,THIS1%PARMS%VAL0_POT*Y0,'H')
          CALL REPORT$R8VAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: RC' &
    &                            ,THIS1%PARMS%RC_POT,'ABOHR')
          CALL REPORT$R8VAL(NFIL,'POTENTIAL PSEUDIZATION PARAMETER: POWER' &
    &                            ,THIS1%PARMS%POW_POT,'')
          CALL REPORT$R8VAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER: VAL(0)' &
    &                            ,THIS1%PARMS%VAL0_CORE*Y0,'H')
          CALL REPORT$R8VAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER: RC' &
    &                            ,THIS1%PARMS%RC_CORE,'ABOHR')
          CALL REPORT$R8VAL(NFIL,'CORE DENSITY PSEUDIZATION PARAMETER: POWER' &
    &                            ,THIS1%PARMS%POW_CORE,'')
        END IF
        IF(.NOT.ASSOCIATED(THIS1%NEXT)) EXIT
        THIS1=>THIS1%NEXT
      ENDDO
!
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$READ()
!     ******************************************************************
!     **  READ SETUP                                                  **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)               :: ISP,NSP1
      CHARACTER(32)            :: NAME
      LOGICAL     ,ALLOCATABLE :: TNEW(:) ! SWITCH BETWEEN INTERNAL ATOM 
                                          ! PROGRAM AND READING THE FILE
!     ******************************************************************
                            CALL TRACE$PUSH('SETUP$READ')
!
!     ==================================================================
!     ==  CREATE SETUPS                                               ==
!     ==================================================================
      CALL ATOMTYPELIST$LENGTH(NSP1)
      ALLOCATE(TNEW(NSP1))
      DO ISP=1,NSP1
        CALL ATOMTYPELIST$NAME(ISP,NAME)
        CALL ATOMTYPELIST$SELECT(NAME)
        CALL ATOMTYPELIST$GETL4('TID',TNEW(ISP))
        CALL ATOMTYPELIST$UNSELECT
        CALL SETUP$SELECT(NAME)
      ENDDO
!
!     ==================================================================
!     ==  READ SETUP FILES                                            ==
!     ==================================================================
      DO ISP=1,NSP1
        CALL SETUP$ISELECT(ISP)
        IF(TNEW(ISP)) THEN
          CALL SETUP_READ_NEW
        ELSE
          CALL SETUP_READ
        END IF
      ENDDO
      DEALLOCATE(TNEW)
                            CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$AEPARTIALWAVES(ISP_,NRX_,LNX_,AEPHI_)
!LEGACY CODE! -> SETUP$GETR8A('AEPHI'
!     ******************************************************************
!     **  RETURN AE PARTIAL WAVES ON THE RADIAL GRID                  **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      INTEGER(4),INTENT(IN) :: LNX_
      REAL(8)   ,INTENT(OUT):: AEPHI_(NRX_,LNX_)
      INTEGER(4)            :: NR
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
      IF(NRX_.NE.NR) THEN
        CALL ERROR$MSG('INCONSISTENT GRID SIZE')
        CALL ERROR$STOP('SETUP$AEPARTIALWAVES')
      END IF
      IF(LNX_.NE.THIS%LNX) THEN
        CALL ERROR$MSG('INCONSISTENT #(PARTIAL WAVES)')
        CALL ERROR$STOP('SETUP$AEPARTIALWAVES')
      END IF
      AEPHI_(:,:)=THIS%AEPHI(:,:)
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$PSPARTIALWAVES(ISP_,NRX_,LNX_,PSPHI_)
!LEGACY CODE! -> SETUP$GETR8A('PSPHI'
!     ******************************************************************
!     **  RETURN PS PARTIAL WAVE ON A RADIAL GRID                     **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      INTEGER(4),INTENT(IN) :: LNX_
      REAL(8)   ,INTENT(OUT):: PSPHI_(NRX_,LNX_)
      INTEGER(4)            :: NR
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
      IF(NRX_.NE.NR) THEN
        CALL ERROR$MSG('INCONSISTENT GRID SIZE')
        CALL ERROR$STOP('SETUP$AEPARTIALWAVES')
      END IF
      IF(LNX_.NE.THIS%LNX) THEN
        CALL ERROR$MSG('INCONSISTENT #(PARTIAL WAVES)')
        CALL ERROR$STOP('SETUP$AEPARTIALWAVES')
      END IF
      PSPHI_(:,:)=THIS%PSPHI(:,:)
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$1COVERLAP(ISP_,LNXX_,DOVER_)
!LEGACY CODE! -> SETUP$GETR8A('DO'
!     ******************************************************************
!     **  RETURN 1-C- OVERLAP OF THE PARTIAL WAVES                    **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: LNXX_
      REAL(8)   ,INTENT(OUT):: DOVER_(LNXX_,LNXX_)
      INTEGER(4)            :: LN1,LN2
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      DO LN1=1,THIS%LNX
        DO LN2=1,THIS%LNX
          DOVER_(LN1,LN2)=THIS%DOVER(LN1,LN2)
        ENDDO
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$1CKINETIC(ISP_,LNXX_,DTKIN_)
!LEGACY CODE! -> SETUP$GETR8A('DEKIN'
!     ******************************************************************
!     **  RETURN 1-C- KINETIC ENERGY OVERLAP OF THE PARTIAL WAVES     **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: LNXX_
      REAL(8)   ,INTENT(OUT):: DTKIN_(LNXX_,LNXX_)
      INTEGER(4)            :: LN1,LN2
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      DO LN1=1,THIS%LNX
        DO LN2=1,THIS%LNX
          DTKIN_(LN1,LN2)=THIS%DTKIN(LN1,LN2)
        ENDDO
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$RCSM(ISP_,RCSM_)
!LEGACY CODE! -> SETUP$GETR8('
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(OUT):: RCSM_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      RCSM_=THIS%RCSM
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$RCBG(ISP_,RCBG_)
!LEGACY CODE! -> SETUP$GETR8('
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(OUT):: RCBG_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      RCBG_=THIS%RCBG
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LNX(ISP_,LNX_)
!LEGACY CODE! -> SETUP$GETI4('
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(OUT):: LNX_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      LNX_=THIS%LNX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMNX(ISP_,LMNX_)
!LEGACY CODE! -> SETUP$GETI4('
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(OUT):: LMNX_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      LMNX_=THIS%LMNX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMNXX(LMNXX_)
!LEGACY CODE! -> SETUP$GETI4('
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT):: LMNXX_
!     ******************************************************************
      LMNXX_=LMNXX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LOFLN(ISP_,LNX_,LOX_)
!LEGACY CODE! -> SETUP$GETI4A('LOX'
!     ******************************************************************
!     **  RETURN NUMBER MAIN ANGULAR MOMENTUM OF PARTIAL WAVES        **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: LNX_
      INTEGER(4),INTENT(OUT):: LOX_(LNX_)
      INTEGER(4)            :: LN
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      IF(LNX_.GT.THIS%LNX) THEN
        CALL ERROR$MSG('LNX ON INPUT TOO SMALL')
        CALL ERROR$OVERFLOW('LNX_',LNX_,THIS%LNX)
        CALL ERROR$STOP('SETUP$LOFLN')
      END IF
      DO LN=1,THIS%LNX
        LOX_(LN)=THIS%LOX(LN)
      ENDDO
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$NSPECIES(NSP_)
!LEGACY CODE! -> SETUP$GETI4
!     ******************************************************************
!     **  RETURN NUMBER OF PARTIAL WAVES                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: NSP_
!     ******************************************************************
      NSP_=NSP
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMRXX(LMRXX_)
!LEGACY CODE! -> SETUP$GETI4
!     ******************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY          **
!     ****************************************************************** 
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: LMRXX_
!     ******************************************************************
      LMRXX_=LMRXX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LNXX(LNXX_)
!LEGACY CODE! -> SETUP$GETI4
!     ******************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY          **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: LNXX_
!     ******************************************************************
      LNXX_=LNXX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$LMRX(ISP_,LMRX_)
!LEGACY CODE! -> SETUP$GETI4
!     ******************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY          **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(OUT):: LMRX_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      LMRX_=THIS%LMRX
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$AEZ(ISP_,AEZ_)
!LEGACY CODE! ->SETUP$GETR8
!     ******************************************************************
!     **  RETURN MAXIMUM ANGULAR MOMENTUM FOR THE 1C-DENSITY          **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      REAL(8)   ,INTENT(OUT):: AEZ_
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      AEZ_=THIS%AEZ
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP_COMPOFG(RCBG,RCSM,GIDG,NG,G0,V0)
!     ******************************************************************
!     **                                                              **
!     **  COMPENSATION DENSITY AND POTENTIAL ON A RADIAL GRID         **
!     **  IN G-SPACE                                                  **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: RCBG
      REAL(8)   ,INTENT(IN) :: RCSM
      INTEGER(4),INTENT(IN) :: GIDG
      INTEGER(4),INTENT(IN) :: NG
      REAL(8)   ,INTENT(OUT):: G0(NG)
      REAL(8)   ,INTENT(OUT):: V0(NG)
      REAL(8)   ,PARAMETER  :: EPSILONGAMMA=1.D-7
      REAL(8)               :: PI
      REAL(8)               :: SVAR1,SVAR2,SVAR3,SVAR4
      REAL(8)               :: BGGAUSS,SMGAUSS
      INTEGER(4)            :: IG
      REAL(8)               :: GARR(NG)
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
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
!     ..................................................................
      SUBROUTINE SETUP_CHECKGAUSS(CELLVOL,RC,TOL,GMAX,RMAX,TCHKR,TCHKG)
!     **                                                              **
!     **                                                              **
      IMPLICIT NONE
      REAL(8)    ,INTENT(IN) :: CELLVOL
      REAL(8)    ,INTENT(IN) :: RC
      REAL(8)    ,INTENT(IN) :: TOL
      REAL(8)    ,INTENT(IN) :: GMAX
      REAL(8)    ,INTENT(IN) :: RMAX
      LOGICAL(4) ,INTENT(OUT):: TCHKR
      LOGICAL(4) ,INTENT(OUT):: TCHKG
      REAL(8)                :: CHECK
!     ******************************************************************
      TCHKG=.TRUE.
      TCHKR=.TRUE.
      CHECK=-1.D0/CELLVOL*EXP(-0.25D0*(RC*GMAX)**2)
      IF(ABS(CHECK).GT.TOL) TCHKG=.FALSE. 
      CHECK=-EXP(-(RMAX/RC)**2)
      IF(ABS(CHECK).GT.TOL) TCHKR=.FALSE. 
      RETURN
      END
!
!     ...........................................INPOT..................
      SUBROUTINE INPOT$LNX(NFIL,LNX)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4),INTENT(OUT) :: LNX
      REAL(8)                :: R1,DEX
      INTEGER(4)             :: NR
!     ******************************************************************
                              CALL TRACE$PUSH('INPOT$LNX')
      REWIND NFIL
      READ(NFIL,FMT='(F15.10,F10.5,2I4)')R1,DEX,NR,LNX
                              CALL TRACE$POP
      RETURN
      END
!
!     ...........................................INPOT..................
      SUBROUTINE INPOT$GRID(NFIL,R1,DEX,NR)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4),INTENT(OUT) :: NR
      REAL(8)   ,INTENT(OUT) :: R1
      REAL(8)   ,INTENT(OUT) :: DEX
!     ******************************************************************
                              CALL TRACE$PUSH('INPOT$GRID')
      REWIND NFIL
      READ(NFIL,FMT='(F15.10,F10.5,2I4)')R1,DEX,NR
                              CALL TRACE$POP
      RETURN
      END
!
!     ...........................................INPOT..................
      SUBROUTINE INPOT$NC(NFIL,NC)
!     ******************************************************************
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4),INTENT(OUT) :: NC
      REAL(8)                :: R1,DEX,PSZ,AEZ,RCSM
      INTEGER(4)             :: NR,LNX
!     ******************************************************************
                              CALL TRACE$PUSH('INPOT$LNX')
      REWIND NFIL
      NC=0
      READ(NFIL,ERR=1000,FMT='(F15.10,F10.5,2I4,2F5.2,F15.12,I5)') &
     &               R1,DEX,NR,LNX,PSZ,AEZ,RCSM,NC
1000 CONTINUE
                              CALL TRACE$POP
      RETURN
      END
!
!     ...........................................INPOT..................
      SUBROUTINE INPOT$READALL(NFIL,NRX,R1,DEX,NR,LNX,LOX,AEZ,PSZ &
     &         ,PSPHI,AEPHI,VADD,RCSM &
     &         ,DTKIN,DOVER,RHOCOR,PSCORR,PRO & !SANTOS040616/BLO
     &         ,TCORE,AEPOT,NC,LB,FB,EB,AEPSI)           !SANTOS040616/BLO
!     &         ,DTKIN,DOVER,IRCCOR,RHOCOR,PSCORR,PRO)
!     ******************************************************************
!     **                                                              **
!     ******************************************* P.E. BLOECHL, 1991 ***
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)  :: NFIL
      INTEGER(4) ,INTENT(IN)  :: NRX
      REAL(8)    ,INTENT(IN)  :: R1
      REAL(8)    ,INTENT(IN)  :: DEX
      INTEGER(4) ,INTENT(IN)  :: NR
      INTEGER(4) ,INTENT(IN)  :: LNX
      INTEGER(4) ,INTENT(IN)  :: NC   !SANTOS
      REAL(8)    ,INTENT(OUT) :: AEZ
      REAL(8)    ,INTENT(OUT) :: PSZ
      REAL(8)    ,INTENT(OUT) :: RCSM
      REAL(8)    ,INTENT(OUT) :: VADD(NRX)
      REAL(8)    ,INTENT(OUT) :: PRO(NRX,LNX)
      INTEGER(4) ,INTENT(OUT) :: LOX(LNX)
!INTEGER(4) ,INTENT(OUT) :: IRCCOR
      REAL(8)    ,INTENT(OUT) :: DTKIN(LNX,LNX)
      REAL(8)    ,INTENT(OUT) :: DOVER(LNX,LNX)
      REAL(8)    ,INTENT(OUT) :: AEPHI(NRX,LNX)
      REAL(8)    ,INTENT(OUT) :: PSPHI(NRX,LNX)
      REAL(8)    ,INTENT(OUT) :: RHOCOR(NRX)
      REAL(8)    ,INTENT(OUT) :: PSCORR(NRX)
! SANTOS040616 BEGIN
      LOGICAL(4) ,INTENT(OUT) :: TCORE
      REAL(8)    ,INTENT(OUT) :: AEPOT(NRX) ! ATOMIC AE POTENTIAL
      INTEGER(4) ,INTENT(OUT) :: LB(NC) ! MAIN ANGULAR MOMENTUM
      REAL(8)    ,INTENT(OUT) :: FB(NC) ! OCCUPATIONS
      REAL(8)    ,INTENT(OUT) :: EB(NC) ! ONE-PARTICLE ENERGIES
      REAL(8)    ,INTENT(OUT) :: AEPSI(NRX,NC) ! CORE STATES
! SANTOS040616 END
      REAL(8)                 :: R11,DEX1
      INTEGER(4)              :: NR1,LNX1,I,IR,LN1,LN2,LN
!     ******************************************************************
                              CALL TRACE$PUSH('INPOT$READALL')
      REWIND NFIL
!     == OPTIONAL VALUE NC FOLLOWING RCSM IS NOT READ 
      READ(NFIL,6000)R11,DEX1,NR1,LNX1,PSZ,AEZ,RCSM
!      READ(NFIL,6000)R11,DEX1,NR1,LNX1,PSZ,AEZ,RCSM,IRCCOR
6000  FORMAT(F15.10,F10.5,2I4,2F5.2,F15.12,I5)
      IF(R11.NE.R1.OR.DEX1.NE.DEX.OR.NR1.NE.NR) THEN
        CALL ERROR$MSG('ONLY ONE TYPE OF RADIAL GRID ALLOWED')
        CALL ERROR$STOP('INPOT')
      END IF


!      IF(IRCCOR.LE.0.OR.IRCCOR.GT.NR) THEN
!      PRINT*,'WARNING! NO MT-RADIUS SPECIFIED FOR ATOM WITH Z=',AEZ
!       IRCCOR=NR-2
!      END IF
      IF(LNX.NE.LNX1) THEN
        CALL ERROR$MSG('LNX OUT OF RANGE')
        CALL ERROR$I4VAL('LNX ON INPUT',LNX)
        CALL ERROR$I4VAL('LNX ON FILE',LNX1)
        CALL ERROR$STOP('INPOT')
      END IF
                              CALL TRACE$PASS('BEFORE LOX')
      READ(NFIL,6020)(LOX(I),I=1,LNX)
6020  FORMAT(14I5)
                              CALL TRACE$PASS('BEFORE VADD')
      READ(NFIL,6100)(VADD(IR),IR=1,NR)
!     ====  RHOCOR = CORE CHARGE DENSITY  ==============================
                              CALL TRACE$PASS('BEFORE RHOCOR')
      READ(NFIL,6100)(RHOCOR(IR),IR=1,NR)
!     ====  PSCORR = PSEUDO CORE CHARGE DENSITY ========================
                              CALL TRACE$PASS('BEFORE PSCORR')
      READ(NFIL,6100)(PSCORR(IR),IR=1,NR)
!     ====  DTKIN = <AEPHI|-DELTA/2|AEPHI> - <PSPHI|-DELTA/2|PSPHI> ====
                              CALL TRACE$PASS('BEFORE DTKIN')
      READ(NFIL,6100)((DTKIN(LN1,LN2),LN1=1,LNX),LN2=1,LNX)
!     PRINT*,'DTKIN ',(DTKIN(LN,LN),LN=1,LNX)
!     ====  DOVER = <AEPHI|AEPHI> - <PSPHI|PSPHI> ======================
                              CALL TRACE$PASS('BEFORE DOVER')
      READ(NFIL,6100)((DOVER(LN1,LN2),LN1=1,LNX),LN2=1,LNX)
                              CALL TRACE$PASS('BEFORE PROJECTORS')
      DO LN=1,LNX
        READ(NFIL,6100)(PRO(IR,LN),IR=1,NR)
        READ(NFIL,6100)(AEPHI(IR,LN),IR=1,NR)
        READ(NFIL,6100)(PSPHI(IR,LN),IR=1,NR)
6100    FORMAT(SP,5E14.8)
      ENDDO
! SANTOS040616 BEGIN
      TCORE=.FALSE.
      AEPOT(:)=0.D0
      LB(:)=0
      FB(:)=0.D0
      EB(:)=0.D0
      AEPSI(:,:)=0.D0
      READ(NFIL,END=1000,FMT='(SP,5E14.8)')(AEPOT(IR),IR=1,NR)
      READ(NFIL,END=1000,FMT='(14I5)')(LB(I),I=1,NC)
      READ(NFIL,END=1000,FMT='(SP,5E14.8)')(FB(I),I=1,NC)
      READ(NFIL,END=1000,FMT='(SP,5E14.8)')(EB(I),I=1,NC)
      DO I=1,NC
        READ(NFIL,END=1000,FMT='(SP,5E14.8)')(AEPSI(IR,I),IR=1,NR)
      ENDDO
      TCORE=.TRUE.
1000  CONTINUE
! SANTOS040616 END
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_CVXMAT(GID,NR,LNX,LOX,AEPHI,NC,LOFC,PSIC,LRHOX,MAT)
!     **************************************************************************
!     **  CORE-VALENCE EXCHANGE MATRIX ELEMENTS                               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)      ! ANGULAR MOMENTA OF PARTIAL WAVES
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX) ! AE PARTIAL WAVES
      INTEGER(4),INTENT(IN) :: NC            ! #(CORE STATES)
      INTEGER(4),INTENT(IN) :: LOFC(NC)      ! ANGULAR MOMENTA OF CORE STATES
      REAL(8)   ,INTENT(IN) :: PSIC(NR,NC)   ! CORE STATES
      INTEGER(4),INTENT(IN) :: LRHOX         ! MAX ANGULAR MOMENTUM OF DENSITY
      REAL(8)   ,INTENT(OUT):: MAT(LNX,LNX)  ! CORE VALENCE X MATRIX ELEMENTS
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
      INTEGER(4)            :: LM1,LC,LRHO,LMC1A,IMC,LMC,LMRHOA,IMRHO,LMRHO
      INTEGER(4)            :: LN1,L1,IC,LN2,L2
      LOGICAL(4),PARAMETER  :: TPRINT=.FALSE.
      REAL(8)               :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      LX=MAXVAL(LOX)
      LMX=(LX+1)**2
      LCX=MAXVAL(LOFC)
      ALLOCATE(FACTOR(LX+1,LRHOX+1,LCX+1))
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == INCLUDE ANGULAR INTEGRATIONS                                         ==
!     ==========================================================================
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
        PRINT*,'NOW THE MATRIX FOR THE CORE VALENCE EXCHANGE INTERACTION'
        DO LN1=1,LNX
          L1=LOX(LN1)
          DO LN2=1,LNX
            L2=LOX(LN2)
            WRITE(*,FMT='(4I3,10F18.10)')LN1,L1,LN2,L2,MAT(LN1,LN2)
          ENDDO
        ENDDO
        STOP
      END IF
      RETURN
      END      
!
!     ..........................................................................
      SUBROUTINE ATOMLIB$READSCNTL(LL_SCNTL)
!     **************************************************************************
!     **************************************************************************
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE(LL_TYPE) ,INTENT(OUT) :: LL_SCNTL
      CHARACTER(128)             :: SCNTLNAME='STP.CNTL'
      TYPE(LL_TYPE) ,SAVE        :: LL_SCNTL_SAVE
      INTEGER(4)                 :: NFIL
      LOGICAL(4)    ,SAVE        :: TINI=.FALSE.
!     **************************************************************************
      IF(TINI) THEN
        LL_SCNTL=LL_SCNTL_SAVE
        RETURN
      END IF
      TINI=.TRUE.
!!$      CALL FILEHANDLER$SETFILE(+'SCNTL',.FALSE.,-SCNTLNAME)
!!$      CALL FILEHANDLER$SETSPECIFICATION(+'SCNTL','STATUS','OLD')
!!$      CALL FILEHANDLER$SETSPECIFICATION(+'SCNTL','POSITION','REWIND')
!!$      CALL FILEHANDLER$SETSPECIFICATION(+'SCNTL','ACTION','READ')
!!$      CALL FILEHANDLER$SETSPECIFICATION(+'SCNTL','FORM','FORMATTED')
!!$      CALL FILEHANDLER$UNIT('SCNTL',NFIL)
!     ==========================================================================
      CALL LINKEDLIST$NEW(LL_SCNTL)
      CALL FILEHANDLER$UNIT('PARMS_STP',NFIL)
      CALL LINKEDLIST$READ(LL_SCNTL,NFIL,'~')
      LL_SCNTL_SAVE=LL_SCNTL
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMLIB$SCNTLLOOKUP(SETUPID,GID,AEZ,ZV,RBOX,LX,RCL &
     &           ,RCSM,POW_POT,RC_POT,VAL0_POT,POW_CORE,RC_CORE,VAL0_CORE)
!     **************************************************************************
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE LINKEDLIST_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: SETUPID
      INTEGER(4)  ,INTENT(OUT):: GID
      REAL(8)     ,INTENT(OUT):: AEZ
      REAL(8)     ,INTENT(OUT):: ZV
      REAL(8)     ,INTENT(OUT):: RBOX
      INTEGER(4)  ,INTENT(IN) :: LX
      REAL(8)     ,INTENT(OUT):: RCL(LX+1)
      REAL(8)     ,INTENT(OUT):: RCSM
      REAL(8)     ,INTENT(OUT):: POW_POT
      REAL(8)     ,INTENT(OUT):: RC_POT
      REAL(8)     ,INTENT(OUT):: VAL0_POT
      REAL(8)     ,INTENT(OUT):: POW_CORE
      REAL(8)     ,INTENT(OUT):: RC_CORE
      REAL(8)     ,INTENT(OUT):: VAL0_CORE
      INTEGER(4)              :: NFIL
      TYPE(LL_TYPE)           :: LL_SCNTL
      INTEGER(4)              :: ITH
      LOGICAL(4)              :: TCHK,TCHK1,TCHK2,TCHK3
      CHARACTER(128)          :: ID
      CHARACTER(2)            :: EL
      REAL(8)                 :: R1
      REAL(8)                 :: DEX
      REAL(8)                 :: SVAR,SVAR1,SVAR2
      INTEGER(4)              :: NR
      INTEGER(4)              :: LENG
      INTEGER(4)              :: L,LN
      INTEGER(4)              :: IZ
      REAL(8)     ,ALLOCATABLE:: R(:)
      REAL(8)     ,ALLOCATABLE:: RCL1(:)
      REAL(8)                 :: RCOV
      REAL(8)                 :: DMIN,DMAX,RX
      REAL(8)                 :: PI,Y0
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL ATOMLIB$READSCNTL(LL_SCNTL)
      CALL LINKEDLIST$SELECT(LL_SCNTL,'~')
      CALL LINKEDLIST$SELECT(LL_SCNTL,'SCNTL')
      ITH=0
      DO 
        ITH=ITH+1
        CALL LINKEDLIST$EXISTL(LL_SCNTL,'SETUP',ITH,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('SETUP DATA NOT FOUND')
          CALL ERROR$STOP('ATOMLIB$READCNTL')
        END IF
        CALL LINKEDLIST$SELECT(LL_SCNTL,'SETUP',ITH)
        CALL LINKEDLIST$EXISTD(LL_SCNTL,'ID',1,TCHK)
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('!SCNTL!SETUP:ID NOT SPECIFIED')
          CALL ERROR$STOP('ATOMLIB$READCNTL')
        END IF
        CALL LINKEDLIST$GET(LL_SCNTL,'ID',1,ID)
        IF(ID.EQ.SETUPID) THEN
!
!         ======================================================================
!         == COLLECT RADIAL GRID                                              ==
!         ======================================================================
          CALL LINKEDLIST$EXISTL(LL_SCNTL,'GRID',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SCNTL!SETUP!GRID NOT SPECIFIED')
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL LINKEDLIST$SELECT(LL_SCNTL,'GRID')

          CALL LINKEDLIST$EXISTD(LL_SCNTL,'DMIN',1,TCHK)
          IF(TCHK) THEN
!           == R1 ================================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'DMIN',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:DMIN NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'DMIN',1,DMIN)
!           == R1 ================================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'DMAX',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:DMAX NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'DMAX',1,DMAX)
!           == R1 ================================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'RMAX',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:RMAX NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'RMAX',1,RX)
!
            CALL RADIAL$GRIDPARAMETERS(DMIN,DMAX,RX,R1,DEX,NR)
!
!           == CHECK FOR CONFLICT ==============================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'R1',1,TCHK1)
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'DEX',1,TCHK2)
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'NR',1,TCHK3)
            IF(TCHK1.OR.TCHK2.OR.TCHK3) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:R1,DEX,NR AND DMIN,DMAX,RX')
              CALL ERROR$MSG('ARE NOT COMPATIBLE')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
          ELSE
!           == R1 ================================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'R1',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:R1 NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'R1',1,R1)
!           == DEX ===============================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'DEX',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:DEX NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'DEX',1,DEX)
!           === NR ===============================================================
            CALL LINKEDLIST$EXISTD(LL_SCNTL,'NR',1,TCHK)
            IF(.NOT.TCHK) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:R NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'NR',1,NR)
          END IF
!
!         == CREATE GRID =======================================================
          CALL RADIAL$NEW('SHLOG',GID)
          CALL RADIAL$SETR8(GID,'R1',R1)
          CALL RADIAL$SETR8(GID,'DEX',DEX)
          CALL RADIAL$SETI4(GID,'NR',NR)
!
          CALL LINKEDLIST$SELECT(LL_SCNTL,'..')
!
!         ======================================================================
!         == IDENTIFY ATOM                                                    ==
!         ======================================================================
!         == COLLECT ATOMIC NUMBER =============================================
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'EL',1,TCHK)
          EL='XX'
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'EL',1,EL)
!
!         == COLLECT ATOMIC NUMBER =============================================
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'Z',1,TCHK)
          IF(TCHK) THEN
            IF(EL.EQ.'XX') THEN 
              CALL ERROR$MSG('!SCNTL!SETUP:AEZ AND EL="XX"')
              CALL ERROR$MSG('MUST NOT BE SPECIFIED SIMULANEOUSLY')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
            CALL LINKEDLIST$GET(LL_SCNTL,'Z',1,AEZ)
          ELSE
            IF(EL.EQ.'XX') THEN
              CALL ERROR$MSG('!SCNTL!SETUP:AEZ NOT SPECIFIED')
              CALL ERROR$MSG('AND ELEMENT NOT SPECIFIED')
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            ELSE
              CALL PERIODICTABLE$GET(EL,'Z',AEZ)
            END IF
          END IF
!
!         == COLLECT #VALENCE ELECTRONS ========================================
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'ZV',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SCNTL!SETUP:ZV NOT SPECIFIED')
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL LINKEDLIST$GET(LL_SCNTL,'ZV',1,ZV)
!
!         == COLLECT BOX RADIUS ================================================
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'RBOX/RCOV',1,TCHK)
          ALLOCATE(R(NR))
          CALL RADIAL$R(GID,NR,R)
          IF(.NOT.TCHK) THEN
            RBOX=R(NR-3)
          ELSE
            CALL LINKEDLIST$GET(LL_SCNTL,'RBOX/RCOV',1,RBOX)
            CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
            RBOX=RBOX*RCOV 
            IF(RBOX.GT.R(NR-3)) THEN
              CALL ERROR$MSG('!SCNTL!SETUP:RBOX  EXCEEDS GRID DIMENSIONS')
              CALL ERROR$R8VAL('RBOX ',RBOX)
              CALL ERROR$R8VAL('R(NR-3)',R(NR-3))
              CALL ERROR$STOP('ATOMLIB$READCNTL')
            END IF
          END IF
          DEALLOCATE(R)
!
!         ======================================================================
!         == CUTOFF RADII FOR PARTIALWAVE CONSTRUCTION                        ==
!         ======================================================================
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'RCL/RCOV',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SCNTL!SETUP:RCL NOT SPECIFIED')
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL LINKEDLIST$SIZE(LL_SCNTL,'RCL/RCOV',1,LENG)
          ALLOCATE(RCL1(LENG))
          CALL LINKEDLIST$GET(LL_SCNTL,'RCL/RCOV',1,RCL1)
          IF(LENG.LT.LX+1) THEN
            CALL ERROR$MSG('!SCNTL!SETUP:RCL ARRAY TOO SHORT')
            CALL ERROR$I4VAL('MAX L ON FILE ',LENG-1)
            CALL ERROR$I4VAL('MAX L REQUESTED',LX)
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
          leng=min(size(rcl),leng)
          RCL=RCL1(1:LENG)*RCOV
!      
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'RCSM/RCOV',1,TCHK)
          IF(TCHK) THEN
            CALL LINKEDLIST$GET(LL_SCNTL,'RCSM/RCOV',1,RCSM)
            RCSM=RCSM*RCOV
          ELSE
            CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
            RCSM=0.25D0*RCOV
          END IF
!
!         ======================================================================
!         == DEFINE POTENTIAL PSEUDIZATION                                    ==
!         ======================================================================
          CALL LINKEDLIST$EXISTL(LL_SCNTL,'POT',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SCNTL!SETUP!POT NOT SPECIFIED')
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL LINKEDLIST$SELECT(LL_SCNTL,'POT')
!  
          POW_POT=2.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'POW',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'POW',1,POW_POT)
!
          VAL0_POT=0.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'VAL0',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'VAL0',1,VAL0_POT)
          VAL0_POT=VAL0_POT/Y0
!
          RC_POT=1.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'RC/RCOV',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'RC/RCOV',1,RC_POT)
          CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
          RC_POT=RC_POT*RCOV

          CALL LINKEDLIST$SELECT(LL_SCNTL,'..')
!
!         ======================================================================
!         == DEFINE CORE PSEUDIZATION                                         ==
!         ======================================================================
          CALL LINKEDLIST$EXISTL(LL_SCNTL,'CORE',1,TCHK)
          IF(.NOT.TCHK) THEN
            CALL ERROR$MSG('!SCNTL!SETUP!CORE NOT SPECIFIED')
            CALL ERROR$STOP('ATOMLIB$READCNTL')
          END IF
          CALL LINKEDLIST$SELECT(LL_SCNTL,'CORE')
!  
          POW_CORE=2.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'POW',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'POW',1,POW_CORE)
!
          VAL0_CORE=0.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'VAL0',1,TCHK)
          IF(TCHK)CALL LINKEDLIST$GET(LL_SCNTL,'VAL0',1,VAL0_CORE)
          VAL0_CORE=VAL0_CORE/Y0
!
          RC_CORE=1.D0
          CALL LINKEDLIST$EXISTD(LL_SCNTL,'RC/RCOV',1,TCHK)
          IF(TCHK) CALL LINKEDLIST$GET(LL_SCNTL,'RC/RCOV',1,RC_CORE)
          CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
          RC_CORE=RC_CORE*RCOV
!
          CALL LINKEDLIST$SELECT(LL_SCNTL,'..')
!
!         ======================================================================        
!         == DONE                                                             ==
!         ======================================================================
          EXIT        
        END IF
        CALL LINKEDLIST$SELECT(LL_SCNTL,'..')
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATOMIC_MAKEPARTIALWAVES(GID,NR,KEY,AEZ,AEPOT &
     &                    ,NB,NC,LOFI,SOFI,NNOFI,EOFI,FOFI &
     &                    ,RBOX,LNX,LOX,RC,AEPHI,PSPHI,NLPHI,PRO,DT,DOVER &
     &                    ,POW_POT,VAL0_POT,RC_POT,RCSM,VADD)
!     **************************************************************************
!     **************************************************************************
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      CHARACTER(*),INTENT(IN) :: KEY
      REAL(8)   ,INTENT(IN) :: AEZ
      REAL(8)   ,INTENT(IN) :: AEPOT(NR)
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NC
      INTEGER(4),INTENT(IN) :: LOFI(NB)
      INTEGER(4),INTENT(IN) :: SOFI(NB)
      INTEGER(4),INTENT(IN) :: NNOFI(NB)
      REAL(8)   ,INTENT(IN) :: EOFI(NB)
      REAL(8)   ,INTENT(IN) :: RBOX
      REAL(8)   ,INTENT(IN) :: FOFI(NB)
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: RC(LNX)
      REAL(8)   ,INTENT(IN) :: POW_POT
      REAL(8)   ,INTENT(IN) :: VAL0_POT
      REAL(8)   ,INTENT(IN) :: RC_POT
      REAL(8)   ,INTENT(IN) :: RCSM
      REAL(8)   ,INTENT(OUT):: VADD(NR)
      REAL(8)   ,INTENT(OUT):: AEPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PSPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: NLPHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: PRO(NR,LNX)
      REAL(8)   ,INTENT(OUT):: DT(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: DOVER(LNX,LNX)
      INTEGER(4),ALLOCATABLE:: NPROL(:)
      INTEGER(4),ALLOCATABLE:: NCL(:)
      REAL(8)               :: DH(LNX,LNX)
      REAL(8)               :: PSPHIPROBARE(LNX,LNX)
      REAL(8)               :: TRANSPHI(LNX,LNX)
      REAL(8)               :: TRANSPRO(LNX,LNX)
      REAL(8)               :: TRANSU(LNX,LNX)
      REAL(8)               :: TRANSUINV(LNX,LNX)
      REAL(8)               :: PSPOT(NR)
      REAL(8)               :: EOFI1(NB)
      REAL(8)               :: EOFLN(LNX)
      REAL(8)               :: UOFI(NR,NB)
      REAL(8)               :: TUOFI(NR,NB)
      REAL(8)               :: AEPSI(NR,NB)
      REAL(8)               :: TNLPHI(NR,LNX)
      REAL(8)               :: TAEPHI(NR,LNX)
      REAL(8)               :: TPSPHI(NR,LNX)
      REAL(8)               :: QN(NR,LNX)
      REAL(8)               :: TQN(NR,LNX)
      REAL(8)               :: BAREPRO(NR,LNX)
REAL(8)               :: PHITEST1(NR,LNX)
      REAL(8)   ,ALLOCATABLE:: PHITEST(:,:)
      REAL(8)   ,ALLOCATABLE:: TPHITEST(:,:)
      REAL(8)   ,ALLOCATABLE:: PRO1(:,:)
      REAL(8)   ,ALLOCATABLE:: DH1(:,:)
      REAL(8)   ,ALLOCATABLE:: DO1(:,:)
      REAL(8)               :: AERHO(NR),PSRHO(NR)
      REAL(8)               :: G(NR),DREL(NR),G1(NR),PHI(NR)
      REAL(8)               :: E
      REAL(8)               :: RC1
      REAL(8)               :: EHOMO
      INTEGER(4)            :: LX
      INTEGER(4)            :: L,IB,LN,IR,IB1,IB2,LN1,LN2,ITER
      INTEGER(4)            :: NV,NPRO,IV,IPRO,IPRO1,IPRO2
      REAL(8)               :: PI,Y0
      REAL(8)               :: R(NR)
      REAL(8)               :: AUX(NR),AUX1(NR)
      REAL(8)   ,allocatable:: AUXarr(:,:)
      REAL(8)               :: VAL
      REAL(8)               :: SVAR,SVAR1,SVAR2
      LOGICAL(4)            :: TREL,TSO,TCHK
      REAL(8)   ,ALLOCATABLE:: A(:,:),AINV(:,:)
      REAL(8)   ,ALLOCATABLE:: PROJ(:)
      LOGICAL   ,PARAMETER  :: TTEST=.TRUE.
      REAL(8)               :: AEPSIF(NR,NB-NC)
      REAL(8)               :: PSPSIF(NR,NB-NC)
      REAL(8)               :: EH,EXC
      INTEGER(4)            :: NN,NN0
      INTEGER(4)            :: NFIL
      CHARACTER(64)         :: STRING
      LOGICAL   ,PARAMETER  :: TWRITE=.TRUE.
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      LX=MAX(MAXVAL(LOX),MAXVAL(LOFI))
      CALL RADIAL$R(GID,NR,R)
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
!
!     == DETERMINE HIGHEST CORE STATE FOR EACH ANGULAR MOMENTUM ================
      ALLOCATE(NCL(0:LX))
      NCL(:)=0
      DO IB=1,NC
        L=LOFI(IB)
        NCL(L)=MAX(NCL(L),IB)
      ENDDO
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
!     == CONSTRUCT PSEUDO POTENTIAL                                           ==
!     ==========================================================================
      CALL ATOMIC_PSEUDIZE(GID,NR,POW_POT,VAL0_POT,RC_POT,AEPOT,PSPOT)
!
!     ==========================================================================
!     == CONSTRUCT NODELESS WAVE FUNCTIONS                                    ==
!     ==========================================================================
      DO L=0,LX
        G(:)=0.D0
        DO IB=1,NB
          IF(LOFI(IB).NE.L) CYCLE
          E=EOFI(IB)
          DREL(:)=0.D0
          IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
          CALL ATOMLIB$BOUNDSTATE(GID,NR,L,0,RBOX,DREL,G,0,AEPOT,E,UOFI(:,IB))
          EOFI1(IB)=E
          TUOFI(:,IB)=G+(E-AEPOT(:)*Y0)*UOFI(:,IB)
          G(:)=UOFI(:,IB)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == REPORT SETTINGS ON WAVE FUNCTIONS                                    ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20," Z=",F5.0," ")')AEZ
        WRITE(6,FMT='(82("="),T20," ENERGIES FOR ATOMIC WAVE FUNCTIONS")')
        DO IB=1,NB
          WRITE(6,FMT='("IB=",I3," L=",I2," F=",F10.5," E[NEW]=",F15.5 &
     &                                               ," E[OLD]=",F15.5)') &
     &                  IB,LOFI(IB),FOFI(IB),EOFI1(IB),EOFI(IB)
        ENDDO
!       CALL SETUP_WRITEPHI('UOFI.DAT',GID,NR,nb,UOFI)
      END IF
!
!     ==========================================================================
!     == DETERMINE ENERGY OF HIGHEST OCCUPIED STATE                           ==
!     ==========================================================================
      EHOMO=-0.5D0
      DO IB=1,NB
        IF(FOFI(IB).GT.0.D0) EHOMO=MAX(EHOMO,EOFI1(IB))
      ENDDO
      IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,EHOMO,DREL)
!
!     ==========================================================================
!     == CONSTRUCT NODELESS PARTIAL WAVES                                     ==
!     ==========================================================================
      DO L=0,LX
        E=EHOMO
        G=0.D0
        IF(NCL(L).NE.0)G(:)=UOFI(:,NCL(L))
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          CALL ATOMLIB$BOUNDSTATE(GID,NR,L,0,RBOX,DREL,G,0,AEPOT,E,NLPHI(:,LN))
          EOFLN(LN)=E
          TNLPHI(:,LN)=G(:)+(E-AEPOT(:)*Y0)*NLPHI(:,LN)
          G(:)=NLPHI(:,LN)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == REPORT SETTINGS ON PARTIAL WAVES                                     ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20," ENERGIES FOR PARTIAL-WAVE CONSTRUCTION")')
        WRITE(6,FMT='("RBOX=",F9.5)')RBOX
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," L=",I2," E=",F10.5," RC=",F6.3)') &
     &                      LN,LOX(LN),EOFLN(LN),RC(LN)
        ENDDO
      END IF
!
!     ==========================================================================
!     == NORMALIZE EACH ANGULAR MOMENTUM SO THAT FIRST PARTIAL WAVE IS NORMAL ==
!     ==========================================================================
      DO L=0,LX
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          PHI(:)=NLPHI(:,LN)
          DO IB=NC,1,-1
            IF(LOFI(IB).NE.L) CYCLE
            AUX(:)=R(:)**2*UOFI(:,IB)*PHI(:)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR1)
            AUX(:)=R(:)**2*UOFI(:,IB)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR2)
            VAL=SVAR1/SVAR2
            PHI(:)=PHI(:)-UOFI(:,IB)*VAL
          ENDDO
          AUX(:)=R(:)**2*PHI(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          VAL=1.D0/SQRT(VAL)
          DO IB=1,NB
            IF(LOFI(IB).NE.L) CYCLE
            UOFI(:,IB) = UOFI(:,IB)*VAL
            TUOFI(:,IB)=TUOFI(:,IB)*VAL
          ENDDO
          DO LN1=1,LNX
            IF(LOX(LN1).NE.L) CYCLE
            NLPHI(:,LN1) = NLPHI(:,LN1)*VAL
            TNLPHI(:,LN1)=TNLPHI(:,LN1)*VAL
          ENDDO
          EXIT
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == CONSTRUCT QN FUNCTIONS        (H-E)|QN>=|UC>                         ==
!     ==========================================================================
      TRANSU(:,:)=0.D0
      DO L=0,LX
        IPRO=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          IPRO=IPRO+1
          SVAR=1.D0
          DO LN1=1,LN
            IF(LOX(LN1).NE.L) CYCLE
            TRANSU(LN1,LN)=TRANSU(LN1,LN)+SVAR
            SVAR=SVAR*(EOFLN(LN)-EOFLN(LN1))
          ENDDO
        ENDDO
      ENDDO
      CALL LIB$INVERTR8(LNX,TRANSU,TRANSUINV)
!
      QN=MATMUL(NLPHI,TRANSU)
      TQN=MATMUL(TNLPHI,TRANSU)

      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20,"  TRANSU ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F15.10)')TRANSU(LN1,:)
        ENDDO
        WRITE(6,FMT='(82("="),T20,"  TRANSU^(-1) ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F15.10)')TRANSUINV(LN1,:)
        ENDDO
      END IF
!
!     ==========================================================================
!     == TEST EQUATION FOR QN                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20,"  TEST QN EQ.  ")')
        DO L=0,LX
          IPRO=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            IPRO=IPRO+1
            PRO(:,LN)=TQN(:,LN)+(AEPOT*Y0-EOFLN(LN))*QN(:,LN)
            IF(NCL(L).NE.0) PRO(:,LN)=PRO(:,LN)-UOFI(:,NCL(L))
            WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|QN>-|UC>     =",F20.15)') &
     &                     LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == CORE-ORTHOGONALIZE TO OBTAIN NODAL PARTIAL WAVES                     ==
!     == NOTE THAT THESE ARE NOT EIGENSTATES AS THEY ARE ONLY ORTHOGONALIZED  ==
!     == TO THE CORE STATES                                                   ==
!     ==                                                                      ==
!     == REMARK: DIRECT ORTHOGONALIZATION IMPROVES THE CORE ORTHOGONALIZATION ==
!     ==   BUT THE RESULT DOES NOT FULFILL THE SCHROEDINGER EQUATION WELL.    ==
!     ==   THE LADDER OF NODELESS WAVE FUNCTION LEADS TO A BETTER SOLUTION    ==
!     ==   OF THE SCHROEDINGER EQUATION, BUT THE CORE-ORTHOGONALITY IS WORSE  ==
!     ==========================================================================
!!$      ALLOCATE(A(NC,NC))
!!$      ALLOCATE(A(NC,NC))
!!$      ALLOCATE(A(NC,NC))
!!$!     == <UC|UC> ===============================================================
!!$      A(:,:)=0.D0
!!$      DO IB1=1,NC
!!$        DO IB2=IB1,NC
!!$          IF(LOFI(IB1).NE.LOFI(IB2)) CYCLE
!!$          AUX(:)=R(:)**2*UOFI(:,IB1)*UOFI(:,IB2)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          A(IB1,IB2)=VAL
!!$          A(IB2,IB1)=VAL
!!$        ENDDO
!!$      ENDDO
!!$!     == INVERT  <UC|UC> =======================================================
!!$      ALLOCATE(AINV(NC,NC))
!!$      CALL LIB$INVERTR8(NC,A,AINV)
!!$      A=MATMUL(A,AINV)
!!$      WRITE(6,FMT='(82("="),T20,"  INVERSION TEST  ")')
!!$      DO LN1=1,LNX
!!$        WRITE(6,FMT='(20F20.15)')A(LN1,:)
!!$      ENDDO
!!$      DEALLOCATE(A)
!!$!     == <UC|QN> ===============================================================
!!$      ALLOCATE(A(NC,LNX))
!!$      A(:,:)=0.D0
!!$      DO IB=1,NC
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.LOFI(IB)) CYCLE
!!$          AUX(:)=R(:)**2*UOFI(:,IB)*QN(:,LN)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          A(IB,LN)=VAL
!!$        ENDDO
!!$      ENDDO
!!$!     == |AEPHI>=|QN> - |UC>(<UC|UC>)^{-1}<UC|QN> ==============================
!!$      AEPHI(:,:)=QN(:,:)-MATMUL(UOFI,MATMUL(AINV,A))
!!$      TAEPHI(:,:)=TQN(:,:)-MATMUL(TUOFI,MATMUL(AINV,A))
!!$!     == AGAIN TO AVOID NUMERICAL ERRORS =======================================
!!$!     == <UC|QN> ===============================================================
!!$      DO IB=1,NC
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.LOFI(IB)) CYCLE
!!$          AUX(:)=R(:)**2*UOFI(:,IB)*AEPHI(:,LN)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          A(IB,LN)=VAL
!!$        ENDDO
!!$      ENDDO
!!$!     == |AEPHI>=|QN> - |UC>(<UC|UC>)^{-1}<UC|QN> ==============================
!!$      AEPHI(:,:)=AEPHI(:,:)-MATMUL(UOFI,MATMUL(AINV,A))
!!$      TAEPHI(:,:)=TAEPHI(:,:)-MATMUL(TUOFI,MATMUL(AINV,A))
!!$      DEALLOCATE(AINV)
!!$      DEALLOCATE(A)
!
!     == USE LADDER OF NODELESS WAVE FUNCTIONS =================================
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
          ENDDO
        ENDDO
        WRITE(6,FMT='(82("="),T20," <UC|AEPHI>  ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F20.10)')A(LN1,:)
        ENDDO
        DEALLOCATE(A)
      END IF
!
!     ==========================================================================
!     == TEST EQUATION FOR ALL-ELECTRON PARTIAL WAVES AEPHI                   ==
!     ==========================================================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20,"  TEST AEPHI EQ.  ")')
        DO LN=1,LNX
          PRO(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
          WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|AEPHI_N> =",F20.15)') &
      &                 LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
        ENDDO
      END IF
!
!     ==========================================================================
!     == CONSTRUCT PSEUDO PARTIAL WAVES                                       ==
!     ==========================================================================
      PSPHI=QN
      TPSPHI=TQN
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
          PHITEST(:,IPRO)=QN(:,LN)
          TPHITEST(:,IPRO)=TQN(:,LN)
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
!
!     ==========================================================================
!     == CONSTRUCT PROJECTOR FUNCTIONS                                        ==
!     ==========================================================================
      DO LN=1,LNX
        BAREPRO(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN)
      ENDDO

      IF(TTEST) THEN
!        CALL SETUP_WRITEPHI('PRO-BARE.DAT',GID,NR,LNX,BAREPRO)
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
          PHITEST(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
        ENDDO
!
        WRITE(6,FMT='(82("="),T20,"  TEST RAW PAW EQUATION  ")')
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," L=",I2," RAW PAW EQ.",F10.5 &
     &                                 ," SCHR. EQ.",F10.5)') &
     &          LN,LOX(LN),MAXVAL(ABS(TPHITEST(:,LN))),MAXVAL(ABS(PHITEST(:,LN)))
        ENDDO
        DEALLOCATE(PHITEST)
        DEALLOCATE(TPHITEST)
      ENDIF
!
!     ==========================================================================
!     == ENFORCE BIORTHOGONALIZATION                                          ==
!     ==========================================================================
      PRO(:,:)=BAREPRO(:,:)
!
!     ==  DETERMINE <PSPHI|PROBARE> ============================================
      PSPHIPROBARE(:,:)=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AUX(:)=R(:)**2*PSPHI(:,LN1)*PRO(:,LN2)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          PSPHIPROBARE(LN1,LN2)=VAL
        ENDDO
      ENDDO
!
      DO L=0,LX
        NPRO=NPROL(L)
        IF(NPRO.EQ.0) CYCLE
        ALLOCATE(A(NPRO,NPRO))
        ALLOCATE(AINV(NPRO,NPRO))
        ALLOCATE(PRO1(NR,NPRO))
        IPRO1=0
        DO LN1=1,LNX
          IF(LOX(LN1).NE.L) CYCLE
          IPRO1=IPRO1+1
          PRO1(:,IPRO1)=PRO(:,LN1)
          IPRO2=0
          DO LN2=1,LNX
            IF(LOX(LN2).NE.L) CYCLE
            IPRO2=IPRO2+1
            A(IPRO1,IPRO2)=PSPHIPROBARE(LN1,LN2)
          ENDDO
        ENDDO
        CALL LIB$INVERTR8(NPRO,A,AINV)
        PRO1=MATMUL(PRO1,AINV)
        DEALLOCATE(A)
        DEALLOCATE(AINV)
        IPRO=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          IPRO=IPRO+1
          PRO(:,LN)=PRO1(:,IPRO)
        ENDDO
        DEALLOCATE(PRO1)
      ENDDO
!
!     ==========================================================================
!     == CHECK BIORTHOGONALIZATION                                            ==
!     ==========================================================================
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AUX(:)=R(:)**2*PSPHI(:,LN1)*PRO(:,LN2)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          IF(LN1.EQ.LN2)VAL=VAL-1.D0
          IF(ABS(VAL).GT.1.D-5) THEN
            CALL ERROR$MSG('BIORTHOGONALIZATION FAILED')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$I4VAL('LN1',LN1)
            CALL ERROR$I4VAL('LN2',LN2)
            CALL ERROR$R8VAL('DEVIATION',VAL)
            CALL ERROR$STOP('XXXX')
          END IF
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DT,DO                                                                ==
!     ==========================================================================
      DT=0.D0
      DOVER=0.D0
      DH=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AUX(:)=R(:)**2*(AEPHI(:,LN1)*TAEPHI(:,LN2)-PSPHI(:,LN1)*TPSPHI(:,LN2))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          DT(LN1,LN2)=VAL
          AUX(:)=R(:)**2*(AEPHI(:,LN1)*AEPHI(:,LN2)-PSPHI(:,LN1)*PSPHI(:,LN2))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          DOVER(LN1,LN2)=VAL
          AUX(:)=R(:)**2*(AEPHI(:,LN1)*AEPOT(:)*Y0*AEPHI(:,LN2) &
      &                  -PSPHI(:,LN1)*PSPOT(:)*Y0*PSPHI(:,LN2))
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          DH(LN1,LN2)=DT(LN1,LN2)+VAL
        ENDDO
      ENDDO

!     == SYMMETRIZE
      DT=0.5D0*(DT+TRANSPOSE(DT))
      DOVER=0.5D0*(DOVER+TRANSPOSE(DOVER))
      DH=0.5D0*(DH+TRANSPOSE(DH))
!
!     ==========================================================================
!     == CHECK2 PAW EQUATION FOR PSEUDO PARTIALWAVES                          ==
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
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
            PROJ(LN1)=VAL
          ENDDO
          WRITE(6,FMT='("LN=",I2," <P|PSPHI>=",10F10.5)')LN,PROJ
!
          PHITEST(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
          TPHITEST(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN)
!TPHITEST(:,LN)=TPHITEST(:,LN)-BAREPRO(:,LN)
          PHITEST1(:,LN)=0.D0          
          DO LN1=1,LNX
            IF(LOX(LN1).NE.LOX(LN)) CYCLE
            SVAR=0.D0
            DO LN2=1,LNX
              IF(LOX(LN2).NE.LOX(LN)) CYCLE
              SVAR=SVAR+(DH(LN1,LN2)-EOFLN(LN2)*DOVER(LN1,LN2))*PROJ(LN2)
            ENDDO
            TPHITEST(:,LN)=TPHITEST(:,LN)+PRO(:,LN1)*SVAR
            PHITEST1(:,LN)=PHITEST1(:,LN)+PRO(:,LN1)*SVAR
          ENDDO
        ENDDO
        WRITE(6,FMT='(82("="),T20,"  TEST PAW EQUATION  ")')
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," L=",I2," PAW EQ.",F20.5 &
     &                                 ," SCHR. EQ.",F20.5," DPRO ",F20.5)') &
     &          LN,LOX(LN),MAXVAL(ABS(TPHITEST(:,LN))),MAXVAL(ABS(PHITEST(:,LN))),MAXVAL(ABS(PHITEST1(:,LN)+BAREPRO(:,LN)))
        ENDDO
!        CALL SETUP_WRITEPHI('pro-test.DAT',GID,NR,LNX,-barepro)
        DEALLOCATE(PROJ)
        DEALLOCATE(PHITEST)
        DEALLOCATE(TPHITEST)
      END IF
!
!     ==========================================================================
!     == BACK TRANSFORM                                                       ==
!     ==========================================================================
      QN=MATMUL(QN,TRANSUINV)
      TQN=MATMUL(TQN,TRANSUINV)
      PSPHI=MATMUL(PSPHI,TRANSUINV)
      TPSPHI=MATMUL(TPSPHI,TRANSUINV)
      AEPHI=MATMUL(AEPHI,TRANSUINV)
      TAEPHI=MATMUL(TAEPHI,TRANSUINV)
      PRO=MATMUL(PRO,TRANSPOSE(TRANSU))
      DT=MATMUL(TRANSPOSE(TRANSUINV),MATMUL(DT,TRANSUINV))
      DOVER=MATMUL(TRANSPOSE(TRANSUINV),MATMUL(DOVER,TRANSUINV))
      DH=MATMUL(TRANSPOSE(TRANSUINV),MATMUL(DH,TRANSUINV))
!
!     == TEST IF BACK TRANSFORM WAS SUCCESSFUL ================================
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20,"  TEST BACK TRANSFORM  ")')
        DO LN=1,LNX
          WRITE(6,FMT='("LN=",I2," L=",I2," DIFF. NDLSS PHI",F10.5 &
     &                                 ," DIFF. KIN.OP NDLSS. PHI ",F10.5)') &
     &          LN,LOX(LN),MAXVAL(ABS(QN(:,LN)-NLPHI(:,LN))) &
     &                    ,MAXVAL(ABS(TQN(:,LN)-TNLPHI(:,LN)))
        ENDDO
      END IF
!
      IF(TTEST) THEN
        WRITE(6,FMT='(82("="),T20,"  DTKIN  ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F12.3)')DT(LN1,:)
        ENDDO
        WRITE(6,FMT='(82("="),T20,"  DOVERLAP  ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F12.3)')DOVER(LN1,:)
        ENDDO
        WRITE(6,FMT='(82("="),T20,"  DOHAMILTONIAN ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F12.3)')DH(LN1,:)
        ENDDO
        WRITE(6,FMT='(82("="),T20,"  <PSPHI|PRO-BARE> ")')
        DO LN1=1,LNX
          WRITE(6,FMT='(20F12.3)')PSPHIPROBARE(LN1,:)
        ENDDO
      END IF
!
!     ==========================================================================
!     == RENORMALIZE WAVE FUNCTIONS AND PROJECTOR FUNCTIONS                   ==
!     ==========================================================================
      DO L=0,LX
        IPRO=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          IPRO=IPRO+1
!         == NORMALIZE PS PARTIAL WAVE =========================================
          IF(IPRO.EQ.1) THEN
            AUX(:)=R(:)**2*PSPHI(:,LN)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
            VAL=VAL+DOVER(LN,LN)
            VAL=1.D0/SQRT(VAL)
          END IF
          PSPHI(:,LN) = PSPHI(:,LN)*VAL
          TPSPHI(:,LN)=TPSPHI(:,LN)*VAL
          AEPHI(:,LN) = AEPHI(:,LN)*VAL
          TAEPHI(:,LN)=TAEPHI(:,LN)*VAL
          NLPHI(:,LN) = NLPHI(:,LN)*VAL
          TNLPHI(:,LN)=TNLPHI(:,LN)*VAL
          QN(:,LN)    =    QN(:,LN)*VAL
          TQN(:,LN)   =   TQN(:,LN)*VAL
          PRO(:,LN)   =   PRO(:,LN)/VAL
          DH(LN,:)=DH(LN,:)*VAL
          DT(LN,:)=DT(LN,:)*VAL
          DOVER(LN,:)=DOVER(LN,:)*VAL
          DH(:,LN)=DH(:,LN)*VAL
          DT(:,LN)=DT(:,LN)*VAL
          DOVER(:,LN)=DOVER(:,LN)*VAL
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == CALCULATE DENSITY FOR UNSCREENING                                    ==
!     ==========================================================================
      AERHO(:)=0.D0
      PSRHO(:)=0.D0
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
        NN0=-1
        G(:)=0.D0
        DO IB=NC+1,NB
          IF(LOFI(IB).NE.L) CYCLE
          IF(NN0.EQ.-1)NN0=NNOFI(IB)
          E=EOFI1(IB)
!
          G(:)=0.D0
          IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
          CALL ATOMLIB$BOUNDSTATE(GID,NR,L,0,RBOX,DREL,G,NNOFI(IB),AEPOT &
       &                             ,E,AEPSIF(:,IB-NC))
          SVAR1=E
          NN=NNOFI(IB)-NN0
          G(:)=0.D0
          CALL ATOMLIB$PAWBOUNDSTATE(GID,NR,L,NN,RBOX,PSPOT,NPRO,PRO1,DH1,DO1,G &
     &                                ,E,PSPSIF(:,IB-NC))
          SVAR2=E
          IF(ABS(SVAR2-EOFI1(IB)).GT.1.D-2) THEN
            CALL ERROR$MSG('INACCURATE BEHAVIOR DURING UNSCREENING PS POTENTIAL')
            CALL ERROR$R8VAL('( PAW-AE)[EV]',(SVAR2-SVAR1)*27.211D0)
            CALL ERROR$R8VAL('(PAW-REF)[EV]',(SVAR2-EOFI1(IB))*27.211D0)
            CALL ERROR$R8VAL('( AE-REF)[EV]',(SVAR1-EOFI1(IB))*27.211D0)
            CALL ERROR$STOP('ATOMLIB_MAKEPARTIALWAVES')
          END IF
          IF(TTEST) THEN
             WRITE(6,FMT='("DEV UNSCREENING LEVELS IN EV ",I5,2F10.5)') &
         &              L,(SVAR1-EOFI1(IB))*27.211D0,(SVAR2-EOFI1(IB))*27.211D0
          END IF
!
          DO IR=1,NR-2
            IF(R(IR).LT.RBOX) CYCLE
            PSPSIF(IR+2:,IB-NC)=0.D0
            AEPSIF(IR+2:,IB-NC)=0.D0
            IF(PSPSIF(IR-2,IB-NC)*AEPSIF(IR-2,IB-NC).LT.0.D0) &
     &                                          AEPSIF(:,IB-NC)=-AEPSIF(:,IB-NC)
            EXIT
          ENDDO
!
!         == CALCULATE PROJECTIONS PROJ=<P|PS-PSI>  ============================
          DO IPRO=1,NPRO
            AUX(:)=R(:)**2*PSPSIF(:,IB-NC)*PRO1(:,IPRO)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,PROJ(IPRO))
          ENDDO
!
!         == NORMALIZE PS WAVE FUNCTION=========================================
          AUX(:)=R(:)**2*PSPSIF(:,IB-NC)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          DO IPRO1=1,NPRO
            DO IPRO2=1,NPRO
              VAL=VAL+PROJ(IPRO1)*DO1(IPRO1,IPRO2)*PROJ(IPRO2)
            ENDDO
          ENDDO
          PSPSIF(:,IB-NC)=PSPSIF(:,IB-NC)/SQRT(VAL)
          PROJ=PROJ/SQRT(VAL)
!
    PRINT*,'PROJ ',L,PROJ
!
!         == NORMALIZE AE WAVE FUNCTION=========================================
          AUX(:)=R(:)**2*AEPSIF(:,IB-NC)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
          AEPSIF(:,IB-NC)=AEPSIF(:,IB-NC)/SQRT(VAL)
!
          PSRHO(:)=PSRHO(:)+FOFI(IB)*PSPSIF(:,IB-NC)**2*Y0
          AERHO(:)=AERHO(:)+FOFI(IB)*AEPSIF(:,IB-NC)**2*Y0
        ENDDO
        DEALLOCATE(DH1)
        DEALLOCATE(DO1)
        DEALLOCATE(PRO1)
        DEALLOCATE(PROJ)
      ENDDO      
!
!!$      IF(TTEST) THEN
!!$        CALL SETUP_WRITEPHI('PSPSIF.DAT',GID,NR,NB-NC,PSPSIF)
!!$        CALL SETUP_WRITEPHI('AEPSIF.DAT',GID,NR,NB-NC,AEPSIF)
!!$        ALLOCATE(AUXARR(NR,2))
!!$        AUXARR(:,1)=AERHO
!!$        AUXARR(:,2)=PSRHO
!!$        CALL SETUP_WRITEPHI('RHO.DAT',GID,NR,2,AUXARR)
!!$        DEALLOCATE(AUXARR)
!!$      END IF
!
!     ==========================================================================
!     == UNSCREENING                                                          ==
!     ==========================================================================
      SVAR=AEZ
      DO IB=1,NC
        SVAR=SVAR-FOFI(IB)
      ENDDO
      CALL ATOMIC_UNSCREEN(GID,NR,RBOX,SVAR,AERHO,PSRHO,PSPOT,RCSM,VADD)
!!$      IF(TTEST) THEN
!!$        ALLOCATE(AUXARR(NR,2))
!!$        AUXARR(:,1)=vadd
!!$        AUXARR(:,2)=pspot
!!$        CALL SETUP_WRITEPHI('vadd.DAT',GID,NR,2,AUXARR)
!!$        DEALLOCATE(AUXARR)
!!$      END IF
!
!     ==========================================================================
!     == WRITE INFORMATION TO FILE                                            ==
!     ==========================================================================
      IF(TWRITE) THEN
        NFIL=8   
        WRITE(STRING,FMT='(F3.0)')AEZ
        STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
!       == AE PARTIAL WAVES
        CALL SETUP_WRITEPHI(-'AEPHI'//TRIM(STRING),GID,NR,LNX,AEPHI)
!
!       == PS PARTIAL WAVES ====================================================
        CALL SETUP_WRITEPHI(-'PSPHI'//TRIM(STRING),GID,NR,LNX,PSPHI)
!
!       == NODELESS PARTIAL WAVES ==============================================
        CALL SETUP_WRITEPHI(-'NLPHI'//TRIM(STRING),GID,NR,LNX,NLPHI)
!
!       == PROJECTOR FUNCTIONS =================================================
        CALL SETUP_WRITEPHI(-'PRO'//TRIM(STRING),GID,NR,LNX,PRO)
!
!       == PROJECTOR FUNCTIONS =================================================
        ALLOCATE(AUXARR(NR,3))
        AUXARR(:,1)=AEPOT
        AUXARR(:,2)=PSPOT
        AUXARR(:,3)=PSPOT-VADD
        CALL SETUP_WRITEPHI(-'POT'//TRIM(STRING),GID,NR,3,AUXARR)
        DEALLOCATE(AUXARR)
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
!     **************************************************************************
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
      LX=MAXVAL(LOX)
!
!     == DETERMINE HIGHEST CORE STATE FOR EACH ANGULAR MOMENTUM ================
      ALLOCATE(NCL(0:LX))
      NCL(:)=0
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
          CALL ERROR$STOP('ATOMIC_MAKEISCATT')
        END IF
        IF(NCL(L).GT.0) THEN
          ISVAR=ISVAR-NNOFI(NCL(L))
        ELSE
          ISVAR=ISVAR+1
        END IF
!       == ISVAR IS NOW THE NUMBER OF VALENCE SHELLS
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          ISCATT(LN)=-ISVAR+1
          ISVAR=ISVAR-1
        ENDDO
      ENDDO
      DEALLOCATE(NCL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_WRITEPHI(FILE,GID,NR,NPHI,PHI)
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
!     ......................................................................
      SUBROUTINE ATOMIC_UNSCREEN(GID,NR,RBOX,AEZ,AERHO,PSRHO,PSPOT,RCSM,VADD)
!     **                                                                  **
!     **                                                                  **
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
      REAL(8)                 :: R(NR)
      REAL(8)                 :: AUX(NR),AUX1(NR),SVAR
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
AUX(:)=PSRHO(:)*R(:)**2
CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,QLM)
PRINT*,'PSRHO CHARGE ',QLM*Y0*4.D0*PI
 AUX(:)=AERHO(:)*R(:)**2
CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,QLM)
PRINT*,'AERHO CHARGE ',QLM*Y0*4.D0*PI
!
!     ========================================================================
!     == MOMENT OF DIFFERENCE CHARGE DENSITY                                ==
!     ========================================================================
      AUX(:)=(AERHO(:)-PSRHO(:))*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,QLM)
      QLM=QLM-AEZ*Y0    ! CHARGE =QM/Y0
PRINT*,'CHARGE ',QLM/Y0,AEZ
!
!     ========================================================================
!     == ADD COMPENSATION DENSITY AND DETERMINE ELECTROSTATIC POTENTIA      ==
!     ========================================================================
      ALPHA=1.D0/RCSM**2
      CALL GAUSSN(0,ALPHA,CL)
      SVAR=QLM*CL
      AUX(:)=PSRHO(:)+SVAR*EXP(-ALPHA*R(:)**2)
 CALL RADIAL$INTEGRATE(GID,NR,AUX*R(:)**2,AUX1)
 CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,QLM)
PRINT*,'QLM ',QLM*Y0*4.D0*PI
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
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE ATOMIC_MAKEPARTIALWAVES_SECOND(GID,NR,KEY,AEPOT,NB,NC,LOFI,SOFI,NNOFI,EOFI,FOFI,PSI &
!!$     &                    ,RBOX,LNX,LOX,RC,AEPHI,PSPHI,PRO,DT,DO &
!!$     &                    ,POW_POT,VAL0_POT,RC_POT,PSPOT)
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: GID
!!$      INTEGER(4),INTENT(IN) :: NR
!!$      CHARACTER(*),INTENT(IN) :: KEY
!!$      REAL(8)   ,INTENT(IN) :: AEPOT(NR)
!!$      INTEGER(4),INTENT(IN) :: NB
!!$      INTEGER(4),INTENT(IN) :: NC
!!$      INTEGER(4),INTENT(IN) :: LOFI(NB)
!!$      INTEGER(4),INTENT(IN) :: SOFI(NB)
!!$      INTEGER(4),INTENT(IN) :: NNOFI(NB)
!!$      REAL(8)   ,INTENT(IN) :: EOFI(NB)
!!$      REAL(8)   ,INTENT(IN) :: RBOX
!!$      REAL(8)   ,INTENT(IN) :: FOFI(NB)
!!$      REAL(8)   ,INTENT(IN) :: PSI(NR,NB)
!!$      INTEGER(4),INTENT(IN) :: LNX
!!$      INTEGER(4),INTENT(IN) :: LOX(LNX)
!!$      REAL(8)   ,INTENT(IN) :: RC(LNX)
!!$      REAL(8)   ,INTENT(IN) :: POW_POT
!!$      REAL(8)   ,INTENT(IN) :: VAL0_POT
!!$      REAL(8)   ,INTENT(IN) :: RC_POT
!!$      REAL(8)   ,INTENT(OUT):: PSPOT(NR)
!!$      REAL(8)   ,INTENT(OUT):: AEPHI(NR,LNX)
!!$      REAL(8)   ,INTENT(OUT):: PSPHI(NR,LNX)
!!$      REAL(8)   ,INTENT(OUT):: PRO(NR,LNX)
!!$      REAL(8)   ,INTENT(OUT):: DT(LNX,LNX)
!!$      REAL(8)   ,INTENT(OUT):: DO(LNX,LNX)
!!$      INTEGER(4),ALLOCATABLE:: NPROL(:)
!!$      INTEGER(4),ALLOCATABLE:: NCL(:)
!!$      REAL(8)               :: DH(LNX,LNX)
!!$      REAL(8)               :: PSPHIPROBARE(LNX,LNX)
!!$      REAL(8)               :: TRANSPHI(LNX,LNX)
!!$      REAL(8)               :: TRANSPRO(LNX,LNX)
!!$      REAL(8)               :: TRANSU(LNX,LNX)
!!$      REAL(8)               :: TRANSUINV(LNX,LNX)
!!$      REAL(8)               :: EOFLN(LNX)
!!$      REAL(8)               :: UOFI(NR,NB)
!!$      REAL(8)               :: TUOFI(NR,NB)
!!$      REAL(8)               :: AEPSI(NR,NB)
!!$      REAL(8)               :: NLPHI(NR,LNX)
!!$      REAL(8)               :: TNLPHI(NR,LNX)
!!$      REAL(8)               :: TAEPHI(NR,LNX)
!!$      REAL(8)               :: TPSPHI(NR,LNX)
!!$      REAL(8)               :: QN(NR,LNX)
!!$      REAL(8)               :: TQN(NR,LNX)
!!$      REAL(8)               :: BAREPRO(NR,LNX)
!!$REAL(8)               :: PHITEST1(NR,LNX)
!!$      REAL(8)   ,ALLOCATABLE:: PHITEST(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: TPHITEST(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: PRO1(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: DH1(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: DO1(:,:)
!!$      REAL(8)               :: AERHO(NR),PSRHO(NR)
!!$      REAL(8)               :: G(NR),DREL(NR),G1(NR),PHI(NR)
!!$      REAL(8)               :: E
!!$      REAL(8)               :: RC1
!!$      REAL(8)               :: EHOMO
!!$      INTEGER(4)            :: LX
!!$      INTEGER(4)            :: L,IB,LN,IR,IB1,LN1,LN2
!!$      INTEGER(4)            :: NV,NPRO,IV,IPRO,IPRO1,IPRO2
!!$      REAL(8)               :: PI,Y0
!!$      REAL(8)               :: R(NR)
!!$      REAL(8)               :: AUX(NR),AUX1(NR)
!!$      REAL(8)               :: VAL
!!$      REAL(8)               :: SVAR,SVAR1,SVAR2
!!$      LOGICAL(4)            ::TREL,TSO,TDOT(LNX),TCHK
!!$      REAL(8)   ,ALLOCATABLE:: A(:,:),AINV(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: PROJ(:)
!!$      LOGICAL(4),PARAMETER  :: TTEST=.TRUE.
!!$      REAL(8)   ,PARAMETER  :: ESTEP=1.D-2
!!$!     **************************************************************************
!!$      PI=4.D0*ATAN(1.D0)
!!$      Y0=1.D0/SQRT(4.D0*PI)
!!$      LX=MAX(MAXVAL(LOX),MAXVAL(LOFI))
!!$      CALL RADIAL$R(GID,NR,R)
!!$!
!!$!     ==========================================================================
!!$!     == RESOLVE KEY                                                          ==
!!$!     ==========================================================================
!!$      TREL=INDEX(KEY,'NONREL').EQ.0
!!$      IF(TREL.AND.INDEX(KEY,'REL').EQ.0) THEN
!!$        CALL ERROR$STOP('SETUP_MAKEPARTIALWAVES')
!!$      END IF
!!$      TSO=INDEX(KEY,'NONSO').EQ.0
!!$      IF(TSO.AND.INDEX(KEY,'SO').EQ.0) THEN
!!$        CALL ERROR$STOP('SETUP_MAKEPARTIALWAVES')
!!$      END IF
!!$!
!!$!     == DETERMINE HIGHEST CORE STATE FOR THIS ANGULAR MOMENTUM ================
!!$      ALLOCATE(NCL(0:LX))
!!$      NCL(:)=0
!!$      DO IB=1,NC
!!$        L=LOFI(IB)
!!$        NCL(L)=MAX(NCL(L),IB)
!!$      ENDDO
!!$!
!!$!     == DETERMINE NUMBER OF PROJECTORS FOR EACH ANGULAR MOMENTUM ==============
!!$      ALLOCATE(NPROL(0:LX))
!!$      NPROL(:)=0
!!$      DO LN=1,LNX
!!$        L=LOX(LN)
!!$        NPROL(L)=NPROL(L)+1
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT PSEUDO POTENTIAL                                           ==
!!$!     ==========================================================================
!!$      CALL PSEUDIZE(GID,NR,POW_POT,VAL0_POT,RC_POT,AEPOT,PSPOT)
!!$OPEN(UNIT=8,FILE='PSPOT.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PSPOT(IR),AEPOT(IR)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT NODELESS WAVE FUNCTIONS PSEUDO POTENTIAL                   ==
!!$!     ==========================================================================
!!$      DO L=0,LX
!!$        G(:)=0.D0
!!$        DO IB=1,NB
!!$          IF(LOFI(IB).NE.L) CYCLE
!!$          E=EOFI(IB)
!!$          DREL(:)=0.D0
!!$          IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
!!$          CALL ATOMLIB$BOUNDSTATE(GID,NR,L,0,RBOX,DREL,G,0,AEPOT,E,UOFI(:,IB))
!!$WRITE(6,FMT='("==",I3,4F10.5)')L,E,EOFI(IB),FOFI(IB)
!!$          TUOFI(:,IB)=G+(E-AEPOT(:)*Y0)*UOFI(:,IB)
!!$          G(:)=UOFI(:,IB)
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == NORMALIZE NODELESS WAVE FUNCTIONS                                    ==
!!$!     ==========================================================================
!!$      AEPSI(:,:)=UOFI(:,:)
!!$      DO IB=NB,1,-1
!!$        L=LOFI(IB)
!!$        DO IB1=IB-1,1,-1
!!$          IF(LOFI(IB1).NE.L) CYCLE
!!$          AUX(:)=R(:)**2*UOFI(:,IB1)*AEPSI(:,IB)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          AEPSI(:,IB)=AEPSI(:,IB)-UOFI(:,IB1)*VAL
!!$        ENDDO
!!$        AUX(:)=R(:)**2*AEPSI(:,IB)**2
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$        CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$        VAL=1.D0/SQRT(VAL)
!!$        AEPSI(:,IB)=AEPSI(:,IB)*VAL
!!$        UOFI(:,IB) =UOFI(:,IB)*VAL
!!$        TUOFI(:,IB)=TUOFI(:,IB)*VAL
!!$      ENDDO
!!$OPEN(UNIT=8,FILE='UOFI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),UOFI(IR,:),AEPSI(IR,:),PSI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='PSI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PSI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == DETERMINE ENERGY OF HIGHEST OCCUPIED STATE                           ==
!!$!     ==========================================================================
!!$      EHOMO=-0.5D0
!!$      DO IB=1,NB
!!$        IF(FOFI(IB).GT.0.D0) EHOMO=MAX(EHOMO,EOFI(IB))
!!$      ENDDO
!!$      IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,EHOMO,DREL)
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT PARTIAL WAVES                                              ==
!!$!     ==========================================================================
!!$      TDOT(:)=.FALSE.
!!$      AEPHI=0.D0
!!$      DO L=0,LX
!!$!
!!$!       == DETERMINE NUMBER OF VALENCE STATES AND INDEX OF LOWEST VALENCE STATE
!!$!       == FOR SEMI-CORE STATES NV=2. OTHERWISE IT IS ZERO OR ONE
!!$        IV=0
!!$        NV=0
!!$        DO IB=NC+1,NB
!!$          IF(LOFI(IB).NE.L) CYCLE
!!$          IF(FOFI(IB).NE.0.D0) NV=NV+1
!!$          IF(IV.EQ.0) IV=IB
!!$        ENDDO
!!$!       == SELECT PHIDOT FUNCTIONS, I.E. TDOT
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$!         TDOT(LN)=IPRO.EQ.NV+1
!!$        ENDDO
!!$!
!!$!       == DETERMINE NODELESS PARTIAL WAVES  ===================================
!!$        E=EHOMO
!!$        G=0.D0
!!$        IF(NCL(L).NE.0)G(:)=UOFI(:,NCL(L))
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IF(TDOT(LN)) THEN   !PHIDOT FUNCTION Q_(N+1)
!!$            CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,0,G,L,E,1,NLPHI(:,LN))
!!$            EOFLN(LN)=E
!!$            TNLPHI(:,LN)=G(:)+(E-AEPOT(:)*Y0)*NLPHI(:,LN)
!!$          ELSE
!!$            CALL ATOMLIB$BOUNDSTATE(GID,NR,L,0,RBOX,DREL,G,0,AEPOT,E,NLPHI(:,LN))
!!$            EOFLN(LN)=E
!!$            TNLPHI(:,LN)=G(:)+(E-AEPOT(:)*Y0)*NLPHI(:,LN)
!!$            G(:)=NLPHI(:,LN)
!!$          END IF
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == REPORT SETTINGS ON PARTIAL WAVES                                     ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20," ENERGIES FOR PARTIAL-WAVE CONSTRUCTION")')
!!$        WRITE(6,FMT='("RBOX=",F9.5)')RBOX
!!$        DO LN=1,LNX
!!$          WRITE(6,FMT='("LN=",I2," L=",I2," E=",F10.5," TDOT=",L1," RC=",F6.3)') &
!!$     &                      LN,LOX(LN),EOFLN(LN),TDOT(LN),RC(LN)
!!$        ENDDO
!!$      END IF
!!$OPEN(UNIT=8,FILE='NLPHI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),NLPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT QN FUNCTIONS        (H-E)|QN>=|UC>                         ==
!!$!     ==========================================================================
!!$      TRANSU(:,:)=0.D0
!!$      DO L=0,LX
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$          IF(TDOT(LN)) THEN
!!$            SVAR1=1.D0
!!$            SVAR2=0.D0
!!$            LN2=0
!!$            IPRO1=0
!!$            DO LN1=1,LN-1
!!$              IF(LOX(LN1).NE.L) CYCLE
!!$              IPRO1=IPRO1+1
!!$              IF(TDOT(LN1)) CYCLE
!!$              IF(LN2.EQ.0) THEN
!!$                LN2=LN1
!!$                CYCLE
!!$              END IF
!!$              SVAR1=SVAR1*(EOFLN(LN)-EOFLN(LN2))
!!$              SVAR2=SVAR2+1.D0/(EOFLN(LN)-EOFLN(LN2))
!!$              SVAR=SVAR1*SVAR2
!!$              TRANSU(LN1,LN)=TRANSU(LN1,LN)+SVAR
!!$              LN2=LN1
!!$            ENDDO
!!$            TRANSU(LN,LN)=TRANSU(LN,LN)+SVAR
!!$          ELSE
!!$            SVAR=1.D0
!!$            DO LN1=1,LN
!!$              IF(LOX(LN1).NE.L) CYCLE
!!$              IF(TDOT(LN1)) CYCLE
!!$              TRANSU(LN1,LN)=TRANSU(LN1,LN)+SVAR
!!$              SVAR=SVAR*(EOFLN(LN)-EOFLN(LN1))
!!$            ENDDO
!!$          END IF
!!$        ENDDO
!!$      ENDDO
!!$      CALL LIB$INVERTR8(LNX,TRANSU,TRANSUINV)
!!$      QN=MATMUL(NLPHI,TRANSU)
!!$      TQN=MATMUL(TNLPHI,TRANSU)
!!$
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20,"  TRANSU ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')TRANSU(LN1,:)
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  TRANSU^(-1) ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')TRANSUINV(LN1,:)
!!$        ENDDO
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     == TEST EQUATION FOR QN                                                 ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20,"  TEST QN EQ.  ")')
!!$        DO L=0,LX
!!$          IPRO=0
!!$          DO LN=1,LNX
!!$            IF(LOX(LN).NE.L) CYCLE
!!$            IPRO=IPRO+1
!!$            PRO(:,LN)=TQN(:,LN)+(AEPOT*Y0-EOFLN(LN))*QN(:,LN)
!!$            IF(TDOT(LN)) THEN
!!$              IF(IPRO.EQ.0) THEN
!!$                IF(NCL(L).NE.0) PRO(:,LN)=PRO(:,LN)-UOFI(:,NCL(L))
!!$                WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|QN>-|UC>     =",F20.15)') &
!!$     &                       LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
!!$              ELSE
!!$                PRO(:,LN)=PRO(:,LN)-QN(:,LN-1)
!!$                WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|QN>-|Q_{N-1}>=",F20.15)') &
!!$     &                       LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
!!$              END IF
!!$            ELSE
!!$              IF(NCL(L).NE.0) PRO(:,LN)=PRO(:,LN)-UOFI(:,NCL(L))
!!$              WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|QN>-|UC>     =",F20.15)') &
!!$     &                     LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
!!$            END IF
!!$          ENDDO
!!$        ENDDO
!!$      END IF
!!$!
!!$OPEN(UNIT=8,FILE='QN.DAT')
!!$DO IR=1,NR
!!$  IF(R(IR).GT.10.D0) EXIT
!!$  WRITE(8,'(30F20.5)')R(IR),MAX(-10.D0,MIN(10.D0,QN(IR,:))),MAX(-10.D0,MIN(10.D0,NLPHI(IR,:)))
!!$ENDDO
!!$CLOSE(8)
!!$
!!$OPEN(UNIT=8,FILE='BAREPRO-1.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PRO(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CORE-ORTHOGONALIZE TO OBTAIN NODAL PARTIAL WAVES                     ==
!!$!     == NOTE THAT THESE ARE NOT EIGENSTATES AS THEY ARE ONLY ORTHOGONALIZED  ==
!!$!     == TO THE CORE STATES                                                   ==
!!$!     ==========================================================================
!!$      AEPHI(:,:)=QN(:,:)    !NLPHI(:,:)
!!$      TAEPHI(:,:)=TQN(:,:)  !TNLPHI(:,:)
!!$      DO L=0,LX
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          DO IB=NC,1,-1
!!$            IF(LOFI(IB).NE.L) CYCLE
!!$            AUX(:)=R(:)**2*UOFI(:,IB)*AEPHI(:,LN)
!!$            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR1)
!!$            AUX(:)=R(:)**2*UOFI(:,IB)**2
!!$            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR2)
!!$            VAL=SVAR1/SVAR2
!!$            AEPHI(:,LN)=AEPHI(:,LN)-UOFI(:,IB)*VAL
!!$            TAEPHI(:,LN)=TAEPHI(:,LN)-TUOFI(:,IB)*VAL
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$OPEN(UNIT=8,FILE='AEPHI-1.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),AEPHI(IR,:),QN(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == TEST EQUATION FOR ALL-ELECTRON PARTIAL WAVES AEPHI                   ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20,"  TEST AEPHI EQ.  ")')
!!$        DO LN=1,LNX
!!$          PRO(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
!!$          IF(TDOT(LN)) THEN
!!$            DO LN1=LN-1,1,-1
!!$              IF(LOX(LN1).NE.LOX(LN)) CYCLE
!!$              PRO(:,LN)=PRO(:,LN)-AEPHI(:,LN1) 
!!$              EXIT
!!$            ENDDO
!!$            WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|AEPHI_N>-|AEPHI_N-1>=",F20.15)') &
!!$      &          LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
!!$          ELSE
!!$            WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|AEPHI_N>            =",F20.15)') &
!!$      &          LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
!!$          END IF
!!$        ENDDO
!!$      END IF
!!$
!!$OPEN(UNIT=8,FILE='BAREPRO-2.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PRO(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT PSEUDO PARTIAL WAVES                                       ==
!!$!     ==========================================================================
!!$      PSPHI=QN
!!$      TPSPHI=TQN
!!$      DO L=0,LX
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          RC1=RC(LN)
!!$        ENDDO
!!$        NPRO=NPROL(L)
!!$        IF(NPRO.EQ.0) CYCLE
!!$        ALLOCATE(PHITEST(NR,NPRO))
!!$        ALLOCATE(TPHITEST(NR,NPRO))
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$          PHITEST(:,IPRO)=QN(:,LN)
!!$          TPHITEST(:,IPRO)=TQN(:,LN)
!!$          IF(TDOT(LN).AND.IPRO.GT.1) THEN
!!$            PHITEST(:,IPRO) =PHITEST(:,IPRO-1) +ESTEP*PHITEST(:,IPRO)
!!$            TPHITEST(:,IPRO)=TPHITEST(:,IPRO-1)+ESTEP*TPHITEST(:,IPRO)
!!$          END IF                 
!!$        ENDDO
!!$        CALL ATOMIC_MAKEPSPHI(GID,NR,RC1,L,NPRO,PHITEST,TPHITEST)
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$          IF(TDOT(LN).AND.IPRO.GT.1) THEN
!!$            PHITEST(:,IPRO) =(PHITEST(:,IPRO)-PHITEST(:,IPRO-1))/ESTEP
!!$            TPHITEST(:,IPRO)=(TPHITEST(:,IPRO)-TPHITEST(:,IPRO-1))/ESTEP
!!$          END IF
!!$          PSPHI(:,LN)=PHITEST(:,IPRO)
!!$          TPSPHI(:,LN)=TPHITEST(:,IPRO)
!!$        ENDDO
!!$        DEALLOCATE(PHITEST)
!!$        DEALLOCATE(TPHITEST)
!!$      ENDDO
!!$
!!$OPEN(UNIT=8,FILE='PSPHI-BARE.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),QN(IR,:),PSPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='TPSPHI-BARE.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),TQN(IR,:),TPSPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == BACK TRANSFORM                                                       ==
!!$!     ==========================================================================
!!$      QN=MATMUL(QN,TRANSUINV)
!!$      TQN=MATMUL(TQN,TRANSUINV)
!!$      PSPHI=MATMUL(PSPHI,TRANSUINV)
!!$      TPSPHI=MATMUL(TPSPHI,TRANSUINV)
!!$      AEPHI=MATMUL(AEPHI,TRANSUINV)
!!$      TAEPHI=MATMUL(TAEPHI,TRANSUINV)
!!$!
!!$!     == TEST IF BACK TRANSFORM WAS SUCCESSFUL ================================
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20,"  TEST BACK TRANSFORM  ")')
!!$        DO LN=1,LNX
!!$          WRITE(6,FMT='("LN=",I2," L=",I2," DIFF. NDLSS PHI",F10.5 &
!!$     &                                 ," DIFF. KIN.OP NDLSS. PHI ",F10.5)') &
!!$     &          LN,LOX(LN),MAXVAL(ABS(QN(:,LN)-NLPHI(:,LN))) &
!!$     &                    ,MAXVAL(ABS(TQN(:,LN)-TNLPHI(:,LN)))
!!$        ENDDO
!!$      END IF
!!$ 
!!$OPEN(UNIT=8,FILE='PSPHI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PSPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='AEPHI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),AEPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='TPSPHI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),TNLPHI(IR,:),TPSPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT PROJECTOR FUNCTIONS                                        ==
!!$!     ==========================================================================
!!$      DO LN=1,LNX
!!$!       == PRO(:,LN)=TNLPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*NLPHI(:,LN) ==========
!!$        BAREPRO(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN)
!!$        TCHK=.FALSE.
!!$        DO LN1=LN-1,1,-1
!!$          IF(LOX(LN1).NE.LOX(LN)) CYCLE
!!$          IF(TDOT(LN1)) CYCLE
!!$!         ==  PRO(:,LN)=PRO(:,LN)-NLPHI(:,LN1) =================================
!!$          BAREPRO(:,LN)=BAREPRO(:,LN)-PSPHI(:,LN1) 
!!$          TCHK=.TRUE.
!!$          EXIT
!!$        ENDDO
!!$        IF(.NOT.TCHK) THEN
!!$          DO IB=NC,1,-1
!!$            IF(LOFI(IB).NE.LOX(LN)) CYCLE
!!$!           == PRO(:,LN)=PRO(:,LN)-UOFI(:,IB) ==================================
!!$            EXIT
!!$          ENDDO
!!$        END IF
!!$      ENDDO
!!$
!!$
!!$OPEN(UNIT=8,FILE='BAREPRO.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),BAREPRO(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CHECK PAW EQUATION FOR PSEUDO PARTIALWAVES                           ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        ALLOCATE(PHITEST(NR,LNX))
!!$        ALLOCATE(TPHITEST(NR,LNX))
!!$        DO LN=1,LNX
!!$          TPHITEST(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN) &
!!$     &                  -BAREPRO(:,LN)
!!$          PHITEST(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
!!$          DO LN1=LN-1,1,-1
!!$            IF(LOX(LN1).NE.LOX(LN)) CYCLE
!!$            IF(TDOT(LN1)) CYCLE
!!$            PHITEST(:,LN)=PHITEST(:,LN)-AEPHI(:,LN1) 
!!$            TPHITEST(:,LN)=TPHITEST(:,LN)-PSPHI(:,LN1) 
!!$            EXIT
!!$          ENDDO
!!$        ENDDO
!!$!
!!$        WRITE(6,FMT='(82("="),T20,"  TEST RAW PAW EQUATION  ")')
!!$        DO LN=1,LNX
!!$          WRITE(6,FMT='("LN=",I2," L=",I2," RAW PAW EQ.",F10.5 &
!!$     &                                 ," SCHR. EQ.",F10.5)') &
!!$     &          LN,LOX(LN),MAXVAL(ABS(TPHITEST(:,LN))),MAXVAL(ABS(PHITEST(:,LN)))
!!$        ENDDO
!!$OPEN(UNIT=8,FILE='AUX.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),TPHITEST(IR,:),PHITEST(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$STOP
!!$        DEALLOCATE(PHITEST)
!!$        DEALLOCATE(TPHITEST)
!!$      ENDIF
!!$
!!$!
!!$!     ==========================================================================
!!$!     == CUT OFF THE TAILS OF THE PROJECTOR FUNCTIONS                         ==
!!$!     ==========================================================================
!!$      PRO(:,:)=BAREPRO(:,:)
!!$      DO LN=1,LNX
!!$        DO IR=1,NR
!!$          IF(R(IR).LT.RC(LN)) CYCLE
!!$          PRO(IR:,LN)=0.D0
!!$          EXIT
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == ENFORCE BIORTHOGONALIZATION                                          ==
!!$!     ==========================================================================
!!$      CALL BIORTHOMATRICES(GID,NR,RBOX,LNX,LOX,PSPHI,PRO,TRANSPHI,TRANSPRO)
!!$        WRITE(6,FMT='(82("="),T20,"  TRANSPRO ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F20.5)')TRANSPRO(LN1,:)
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  TRANSPHI ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F20.5)')TRANSPHI(LN1,:)
!!$        ENDDO
!!$
!!$!     ==  DETERMINE <PSPHI|PROBARE> ============================================
!!$      PSPHIPROBARE(:,:)=0.D0
!!$      DO LN1=1,LNX
!!$        DO LN2=1,LNX
!!$          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
!!$          AUX(:)=R(:)**2*PSPHI(:,LN1)*PRO(:,LN2)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          PSPHIPROBARE(LN1,LN2)=VAL
!!$        ENDDO
!!$      ENDDO
!!$!
!!$      DO L=0,LX
!!$        NPRO=NPROL(L)
!!$        IF(NPRO.EQ.0) CYCLE
!!$        ALLOCATE(A(NPRO,NPRO))
!!$        ALLOCATE(AINV(NPRO,NPRO))
!!$        ALLOCATE(PRO1(NR,NPRO))
!!$        IPRO1=0
!!$        DO LN1=1,LNX
!!$          IF(LOX(LN1).NE.L) CYCLE
!!$          IPRO1=IPRO1+1
!!$          PRO1(:,IPRO1)=PRO(:,LN1)
!!$          IPRO2=0
!!$          DO LN2=1,LNX
!!$            IF(LOX(LN2).NE.L) CYCLE
!!$            IPRO2=IPRO2+1
!!$            A(IPRO1,IPRO2)=PSPHIPROBARE(LN1,LN2)
!!$          ENDDO
!!$        ENDDO
!!$        CALL LIB$INVERTR8(NPRO,A,AINV)
!!$DO IPRO1=1,NPRO
!!$  WRITE(*,FMT='("A",20F15.10)')A(IPRO1,:)
!!$ENDDO
!!$DO IPRO1=1,NPRO
!!$  WRITE(*,FMT='("AINV",20E15.3)')AINV(IPRO1,:)
!!$ENDDO
!!$A=MATMUL(AINV,A)
!!$DO IPRO1=1,NPRO
!!$  WRITE(*,FMT='("AINV*A",20F15.10)')A(IPRO1,:)
!!$ENDDO
!!$
!!$        A=MATMUL(A,AINV)
!!$        PRO1=MATMUL(PRO1,AINV)
!!$        DEALLOCATE(A)
!!$        DEALLOCATE(AINV)
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$          PRO(:,LN)=PRO1(:,IPRO)
!!$        ENDDO
!!$        DEALLOCATE(PRO1)
!!$      ENDDO
!!$OPEN(UNIT=8,FILE='PRO.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PRO(IR,:) !*R(IR)**2
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='AEPSPHI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),AEPHI(IR,:),PSPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='AEPSDIFF.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),(AEPHI(IR,:)-PSPHI(IR,:))*R(IR)**2
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CHECK BIORTHOGONALIZATION                                            ==
!!$!     ==========================================================================
!!$      DO LN1=1,LNX
!!$        DO LN2=1,LNX
!!$          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
!!$          AUX(:)=R(:)**2*PSPHI(:,LN1)*PRO(:,LN2)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          IF(LN1.EQ.LN2)VAL=VAL-1.D0
!!$          IF(ABS(VAL).GT.1.D-5) THEN
!!$            CALL ERROR$MSG('BIORTHOGONALIZATION FAILED')
!!$            CALL ERROR$I4VAL('L',L)
!!$            CALL ERROR$I4VAL('LN1',LN1)
!!$            CALL ERROR$I4VAL('LN2',LN2)
!!$            CALL ERROR$R8VAL('DEVIATION',VAL)
!!$            CALL ERROR$STOP('XXXX')
!!$          END IF
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == DT,DO                                                                ==
!!$!     ==========================================================================
!!$      DT=0.D0
!!$      DO=0.D0
!!$      DH=0.D0
!!$      DO LN1=1,LNX
!!$        DO LN2=1,LNX
!!$          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
!!$          AUX(:)=R(:)**2*(AEPHI(:,LN1)*TAEPHI(:,LN2)-PSPHI(:,LN1)*TPSPHI(:,LN2))
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          DT(LN1,LN2)=VAL
!!$          AUX(:)=R(:)**2*(AEPHI(:,LN1)*AEPHI(:,LN2)-PSPHI(:,LN1)*PSPHI(:,LN2))
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          DO(LN1,LN2)=VAL
!!$          AUX(:)=R(:)**2*(AEPHI(:,LN1)*AEPOT(:)*Y0*AEPHI(:,LN2) &
!!$      &                  -PSPHI(:,LN1)*PSPOT(:)*Y0*PSPHI(:,LN2))
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          DH(LN1,LN2)=DT(LN1,LN2)+VAL
!!$        ENDDO
!!$      ENDDO
!!$!
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20,"  DTKIN  ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')DT(LN1,:)
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  DOVERLAP  ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')DO(LN1,:)
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  DOHAMILTONIAN ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')DH(LN1,:)
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  <PSPHI|PRO-BARE> ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')PSPHIPROBARE(LN1,:)
!!$        ENDDO
!!$!
!!$        WRITE(6,FMT='(82("="),T20,"  TEST <PSPHI|PRO-BARE>+DH-DO-DO_N-1=0 ")')
!!$        ALLOCATE(A(LNX,LNX))        
!!$        A=PSPHIPROBARE+DH
!!$        DO LN=1,LNX
!!$          A(:,LN)=A(:,LN)-DO(:,LN)*EOFLN(LN)
!!$        ENDDO
!!$        DO L=0,LX
!!$          LN1=0
!!$          DO LN=1,LNX
!!$            IF(LOX(LN).NE.L) CYCLE
!!$            IF(LN1.NE.0)A(:,LN)=A(:,LN)-DO(:,LN1)
!!$            LN1=LN
!!$          ENDDO
!!$        ENDDO
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')A(LN1,:)
!!$        ENDDO
!!$OPEN(UNIT=8,FILE='AUX1.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),(SUM(PRO(IR,:)*A(:,LN)),LN=1,LNX)
!!$ENDDO
!!$CLOSE(8)
!!$!        DEALLOCATE(A)
!!$      END IF
!!$
!!$!     == SYMMETRIZE
!!$      DT=0.5D0*(DT+TRANSPOSE(DT))
!!$      DO=0.5D0*(DO+TRANSPOSE(DO))
!!$      DH=0.5D0*(DH+TRANSPOSE(DH))
!!$!
!!$!     ==========================================================================
!!$!     == CHECK2 PAW EQUATION FOR PSEUDO PARTIALWAVES                          ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        ALLOCATE(TPHITEST(NR,LNX))    ! HOLDS TEST FOR PSEUDO
!!$        ALLOCATE(PHITEST(NR,LNX))     ! HOLDS TEST FOR ALL-ELECTRON 
!!$        DO LN=1,LNX
!!$          TPHITEST(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN)
!!$! TPHITEST(:,LN)=TPHITEST(:,LN)-BAREPRO(:,LN)
!!$          PHITEST1(:,LN)=0.D0          
!!$          DO LN1=1,LNX
!!$            IF(LOX(LN1).NE.LOX(LN)) CYCLE
!!$             TPHITEST(:,LN)=TPHITEST(:,LN)+PRO(:,LN1)*(DH(LN1,LN)-EOFLN(LN)*DO(LN1,LN))
!!$             PHITEST1(:,LN)=PHITEST1(:,LN)+PRO(:,LN1)*(DH(LN1,LN)-EOFLN(LN)*DO(LN1,LN))
!!$          ENDDO
!!$!
!!$          PHITEST(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
!!$          DO LN1=LN-1,1,-1
!!$            IF(LOX(LN1).NE.LOX(LN)) CYCLE
!!$            IF(TDOT(LN1)) CYCLE
!!$!           ==  PRO(:,LN)=PRO(:,LN)-NLPHI(:,LN1) =================================
!!$            PHITEST(:,LN)=PHITEST(:,LN)-AEPHI(:,LN1) 
!!$            TPHITEST(:,LN)=TPHITEST(:,LN)-PSPHI(:,LN1) 
!!$            DO LN2=1,LNX
!!$              IF(LOX(LN2).NE.LOX(LN)) CYCLE
!!$              TPHITEST(:,LN)=TPHITEST(:,LN)-PRO(:,LN2)*DO(LN2,LN1)
!!$              PHITEST1(:,LN)=PHITEST1(:,LN)-PRO(:,LN2)*DO(LN2,LN1)
!!$            ENDDO
!!$            EXIT
!!$          ENDDO
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  TEST PAW EQUATION  ")')
!!$        DO LN=1,LNX
!!$          WRITE(6,FMT='("LN=",I2," L=",I2," PAW EQ.",F10.5 &
!!$     &                                 ," SCHR. EQ.",F10.5)') &
!!$     &          LN,LOX(LN),MAXVAL(ABS(TPHITEST(:,LN))),MAXVAL(ABS(PHITEST(:,LN)))
!!$        ENDDO
!!$OPEN(UNIT=8,FILE='AUX.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PHITEST1(IR,:),-BAREPRO(IR,:),(SUM(PRO(IR,:)*A(:,LN)),LN=1,LNX)
!!$ENDDO
!!$CLOSE(8)
!!$STOP
!!$        DEALLOCATE(PHITEST)
!!$        DEALLOCATE(TPHITEST)
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     == UNSCREENING                                                          ==
!!$!     ==========================================================================
!!$      AERHO(:)=0.D0
!!$      PSRHO(:)=0.D0
!!$      DO L=0,LX
!!$        NPRO=NPROL(L)
!!$        ALLOCATE(DH1(NPRO,NPRO))
!!$        ALLOCATE(DO1(NPRO,NPRO))
!!$        ALLOCATE(PRO1(NR,NPRO))
!!$        ALLOCATE(PROJ(NPRO))
!!$        IPRO1=0
!!$        DO LN1=1,LNX
!!$          IF(LOX(LN1).NE.L) CYCLE
!!$          IPRO1=IPRO1+1
!!$          PRO1(:,IPRO1)=PRO(:,LN1)
!!$          IPRO2=0
!!$          DO LN2=1,LNX
!!$            IF(LOX(LN2).NE.L) CYCLE
!!$            IPRO2=IPRO2+1
!!$            DH1(IPRO1,IPRO2)=DH(LN1,LN2)
!!$            DO1(IPRO1,IPRO2)=DO(LN1,LN2)
!!$          ENDDO
!!$        ENDDO
!!$        G(:)=0.D0
!!$        DO IB=NC+1,NB
!!$          IF(LOFI(IB).NE.L) CYCLE
!!$          E=EOFI(IB)
!!$          G(:)=0.D0
!!$          CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO1,DH1,DO1,G,PHI)
!!$
!!$          CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,0,G,L,E,1,AUX)
!!$          DO IR=1,NR
!!$            IF(R(IR).LT.RBOX) CYCLE
!!$            PHI(IR:)=0.D0
!!$          ENDDO
!!$PRINT*,'RBOX ',RBOX
!!$OPEN(UNIT=8,FILE='TEST.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PHI(IR),AUX(IR)
!!$ENDDO
!!$CLOSE(8)
!!$STOP
!!$          DO IPRO=1,NPRO
!!$            AUX(:)=R(:)**2*PHI(:)*PRO1(:,IPRO)
!!$            CALL RADIAL$INTEGRAL(GID,NR,AUX,PROJ(NPRO))
!!$          ENDDO
!!$PRINT*,'PROJ ',L,PROJ
!!$          AUX(:)=R(:)**2*PHI(:)**2
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!!$PRINT*,'VAL1 ',L,VAL
!!$          DO IPRO1=1,NPRO
!!$            DO IPRO2=1,NPRO
!!$              VAL=VAL+PROJ(IPRO1)*DO1(IPRO1,IPRO2)*PROJ(IPRO2)
!!$            ENDDO
!!$          ENDDO
!!$PRINT*,'VAL2 ',L,VAL
!!$          PHI(:)=PHI(:)/SQRT(VAL)
!!$          PSRHO(:)=PSRHO(:)+FOFI(IB)*PHI(:)**2
!!$          AERHO(:)=AERHO(:)+FOFI(IB)*AEPSI(:,IB)**2
!!$        ENDDO
!!$        DEALLOCATE(DH1)
!!$        DEALLOCATE(DO1)
!!$        DEALLOCATE(PRO1)
!!$        DEALLOCATE(PROJ)
!!$      ENDDO      
!!$OPEN(UNIT=8,FILE='RHO.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),AERHO(IR),PSRHO(IR)
!!$ENDDO
!!$CLOSE(8)
!!$STOP
!!$
!!$!      CALL ATOMLIB$BOXVOFRHO(GID,NR,RAD,AEZ,RHO,POT,EH,EXC)
!!$
!!$      RETURN
!!$      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE ATOMIC_MAKEPARTIALWAVES_FIRST(GID,NR,KEY,AEPOT,NB,NC,LOFI,SOFI,NNOFI,EOFI,FOFI,PSI &
!!$     &                    ,RBOX,LNX,LOX,RC,AEPHI,PSPHI,PRO,DT,DO &
!!$     &                    ,POW_POT,VAL0_POT,RC_POT,PSPOT)
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: GID
!!$      INTEGER(4),INTENT(IN) :: NR
!!$      CHARACTER(*),INTENT(IN) :: KEY
!!$      REAL(8)   ,INTENT(IN) :: AEPOT(NR)
!!$      INTEGER(4),INTENT(IN) :: NB
!!$      INTEGER(4),INTENT(IN) :: NC
!!$      INTEGER(4),INTENT(IN) :: LOFI(NB)
!!$      INTEGER(4),INTENT(IN) :: SOFI(NB)
!!$      INTEGER(4),INTENT(IN) :: NNOFI(NB)
!!$      REAL(8)   ,INTENT(IN) :: EOFI(NB)
!!$      REAL(8)   ,INTENT(IN) :: RBOX
!!$      REAL(8)   ,INTENT(IN) :: FOFI(NB)
!!$      REAL(8)   ,INTENT(IN) :: PSI(NR,NB)
!!$      INTEGER(4),INTENT(IN) :: LNX
!!$      INTEGER(4),INTENT(IN) :: LOX(LNX)
!!$      REAL(8)   ,INTENT(IN) :: RC(LNX)
!!$      REAL(8)   ,INTENT(IN) :: POW_POT
!!$      REAL(8)   ,INTENT(IN) :: VAL0_POT
!!$      REAL(8)   ,INTENT(IN) :: RC_POT
!!$      REAL(8)   ,INTENT(OUT):: PSPOT(NR)
!!$      REAL(8)   ,INTENT(OUT):: AEPHI(NR,LNX)
!!$      REAL(8)   ,INTENT(OUT):: PSPHI(NR,LNX)
!!$      REAL(8)   ,INTENT(OUT):: PRO(NR,LNX)
!!$      REAL(8)   ,INTENT(OUT):: DT(LNX,LNX)
!!$      REAL(8)   ,INTENT(OUT):: DO(LNX,LNX)
!!$      REAL(8)               :: DH(LNX,LNX)
!!$      REAL(8)               :: PSPHIPROBARE(LNX,LNX)
!!$      REAL(8)               :: TRANSPHI(LNX,LNX)
!!$      REAL(8)               :: TRANSPRO(LNX,LNX)
!!$      REAL(8)               :: EOFLN(LNX)
!!$      REAL(8)               :: UOFI(NR,NB)
!!$      REAL(8)               :: TUOFI(NR,NB)
!!$      REAL(8)               :: AEPSI(NR,NB)
!!$      REAL(8)               :: NLPHI(NR,LNX)
!!$      REAL(8)               :: TNLPHI(NR,LNX)
!!$      REAL(8)               :: TAEPHI(NR,LNX)
!!$      REAL(8)               :: TPSPHI(NR,LNX)
!!$      REAL(8)               :: QN(NR,LNX)
!!$      REAL(8)               :: TQN(NR,LNX)
!!$      REAL(8)               :: BAREPRO(NR,LNX)
!!$REAL(8)               :: PHITEST1(NR,LNX)
!!$      REAL(8)   ,ALLOCATABLE:: PHITEST(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: TPHITEST(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: PRO1(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: DH1(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: DO1(:,:)
!!$      REAL(8)               :: AERHO(NR),PSRHO(NR)
!!$      REAL(8)               :: G(NR),DREL(NR),G1(NR),PHI(NR)
!!$      REAL(8)               :: E
!!$      REAL(8)               :: RC1
!!$      REAL(8)               :: EHOMO
!!$      INTEGER(4)            :: LX
!!$      INTEGER(4)            :: L,IB,LN,IR,IB1,LN1,LN2
!!$      INTEGER(4)            :: NC1,NV,NPRO,IV,IPRO,IPRO1,IPRO2
!!$      REAL(8)               :: PI,Y0
!!$      REAL(8)               :: R(NR)
!!$      REAL(8)               :: AUX(NR),AUX1(NR)
!!$      REAL(8)               :: VAL
!!$      REAL(8)               :: SVAR,SVAR1,SVAR2
!!$      LOGICAL(4)            ::TREL,TSO,TDOT(LNX),TCHK
!!$      REAL(8)   ,ALLOCATABLE:: A(:,:),AINV(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: PROJ(:)
!!$      LOGICAL(4),PARAMETER  :: TTEST=.TRUE.
!!$!     **************************************************************************
!!$      PI=4.D0*ATAN(1.D0)
!!$      Y0=1.D0/SQRT(4.D0*PI)
!!$      LX=MAX(MAXVAL(LOX),MAXVAL(LOFI))
!!$      CALL RADIAL$R(GID,NR,R)
!!$!
!!$!     ==========================================================================
!!$!     == RESOLVE KEY                                                          ==
!!$!     ==========================================================================
!!$      TREL=INDEX(KEY,'NONREL').EQ.0
!!$      IF(TREL.AND.INDEX(KEY,'REL').EQ.0) THEN
!!$        CALL ERROR$STOP('SETUP_MAKEPARTIALWAVES')
!!$      END IF
!!$      TSO=INDEX(KEY,'NONSO').EQ.0
!!$      IF(TSO.AND.INDEX(KEY,'SO').EQ.0) THEN
!!$        CALL ERROR$STOP('SETUP_MAKEPARTIALWAVES')
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT PSEUDO POTENTIAL                                           ==
!!$!     ==========================================================================
!!$      CALL PSEUDIZE(GID,NR,POW_POT,VAL0_POT,RC_POT,AEPOT,PSPOT)
!!$OPEN(UNIT=8,FILE='PSPOT.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PSPOT(IR),AEPOT(IR)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT NODELESS WAVE FUNCTIONS PSEUDO POTENTIAL                   ==
!!$!     ==========================================================================
!!$      DO L=0,LX
!!$        G(:)=0.D0
!!$        DO IB=1,NB
!!$          IF(LOFI(IB).NE.L) CYCLE
!!$          E=EOFI(IB)
!!$          DREL(:)=0.D0
!!$          IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,E,DREL)
!!$          CALL ATOMLIB$BOUNDSTATE(GID,NR,L,0,RBOX,DREL,G,0,AEPOT,E,UOFI(:,IB))
!!$WRITE(6,FMT='("==",I3,4F10.5)')L,E,EOFI(IB),FOFI(IB)
!!$          TUOFI(:,IB)=G+(E-AEPOT(:)*Y0)*UOFI(:,IB)
!!$          G(:)=UOFI(:,IB)
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == NORMALIZE NODELESS WAVE FUNCTIONS                                    ==
!!$!     ==========================================================================
!!$      AEPSI(:,:)=UOFI(:,:)
!!$      DO IB=NB,1,-1
!!$        L=LOFI(IB)
!!$        DO IB1=IB-1,1,-1
!!$          IF(LOFI(IB1).NE.L) CYCLE
!!$          AUX(:)=R(:)**2*UOFI(:,IB1)*AEPSI(:,IB)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          AEPSI(:,IB)=AEPSI(:,IB)-UOFI(:,IB1)*VAL
!!$        ENDDO
!!$        AUX(:)=R(:)**2*AEPSI(:,IB)**2
!!$        CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$        CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$        VAL=1.D0/SQRT(VAL)
!!$        AEPSI(:,IB)=AEPSI(:,IB)*VAL
!!$        UOFI(:,IB) =UOFI(:,IB)*VAL
!!$        TUOFI(:,IB)=TUOFI(:,IB)*VAL
!!$      ENDDO
!!$OPEN(UNIT=8,FILE='UOFI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),UOFI(IR,:),AEPSI(IR,:),PSI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='PSI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PSI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == DETERMINE ENERGY OF HIGHEST OCCUPIED STATE                           ==
!!$!     ==========================================================================
!!$      EHOMO=-0.5D0
!!$      DO IB=1,NB
!!$        IF(FOFI(IB).GT.0.D0) EHOMO=MAX(EHOMO,EOFI(IB))
!!$      ENDDO
!!$     IF(TREL)CALL SCHROEDINGER$DREL(GID,NR,AEPOT,EHOMO,DREL)
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT PARTIAL WAVES                                              ==
!!$!     ==========================================================================
!!$      TDOT(:)=.FALSE.
!!$      AEPHI=0.D0
!!$      DO L=0,LX
!!$!
!!$!       == DETERMINE HIGHEST CORE STATE FOR THIS ANGULAR MOMENTUM ==============
!!$        NC1=0
!!$        DO IB=NC,1,-1
!!$          IF(LOFI(IB).NE.L) CYCLE
!!$          NC1=IB
!!$          EXIT
!!$        ENDDO
!!$!       == DETERMINE NUMBER OF VALENCE STATES AND INDEX OF LOWEST VALENCE STATE
!!$!       == FOR SEMI-CORE STATES NV=2. OTHERWISE IT IS ZERO OR ONE
!!$        IV=0
!!$        NV=0
!!$        DO IB=NC+1,NB
!!$          IF(LOFI(IB).NE.L) CYCLE
!!$          IF(FOFI(IB).NE.0.D0) NV=NV+1
!!$          IF(IV.EQ.0) IV=IB
!!$        ENDDO
!!$!       == SELECT PHIDOT FUNCTIONS, I.E. TDOT
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$          TDOT(LN)=IPRO.EQ.NV+1
!!$        ENDDO
!!$TDOT(2)=.FALSE.
!!$TDOT(3)=.TRUE.
!!$TDOT=.FALSE.
!!$!
!!$!       == DETERMINE NODELESS PARTIAL WAVES  ===================================
!!$        E=EHOMO
!!$        G=0.D0
!!$        IF(NC1.NE.0)G(:)=UOFI(:,NC1)
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IF(TDOT(LN)) THEN   !PHIDOT FUNCTION Q_(N+1)
!!$            CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,0,G,L,E,1,NLPHI(:,LN))
!!$            EOFLN(LN)=E
!!$            TNLPHI(:,LN)=G(:)+(E-AEPOT(:)*Y0)*NLPHI(:,LN)
!!$          ELSE
!!$            CALL ATOMLIB$BOUNDSTATE(GID,NR,L,0,RBOX,DREL,G,0,AEPOT,E,NLPHI(:,LN))
!!$            EOFLN(LN)=E
!!$            TNLPHI(:,LN)=G(:)+(E-AEPOT(:)*Y0)*NLPHI(:,LN)
!!$            G(:)=NLPHI(:,LN)
!!$          END IF
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == REPORT SETTINGS ON PARTIAL WAVES                                     ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20," ENERGIES FOR PARTIAL-WAVE CONSTRUCTION")')
!!$        WRITE(6,FMT='("RBOX=",F9.5)')RBOX
!!$        DO LN=1,LNX
!!$          WRITE(6,FMT='("LN=",I2," L=",I2," E=",F10.5," TDOT=",L1," RC=",F6.3)') &
!!$     &                      LN,LOX(LN),EOFLN(LN),TDOT(LN),RC(LN)
!!$        ENDDO
!!$      END IF
!!$OPEN(UNIT=8,FILE='NLPHI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),NLPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT QN FUNCTIONS        (H-E)|QN>=|UC>                         ==
!!$!     ==========================================================================
!!$      DO L=0,LX
!!$!
!!$!       == DETERMINE HIGHEST CORE STATE FOR THIS ANGULAR MOMENTUM ==============
!!$        NC1=0
!!$        DO IB=NC,1,-1
!!$          IF(LOFI(IB).NE.L) CYCLE
!!$          NC1=IB
!!$          EXIT
!!$        ENDDO
!!$!
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$          IF(TDOT(LN)) THEN
!!$            QN(:,LN)=0.D0
!!$            TQN(:,LN)=0.D0
!!$            SVAR1=1.D0
!!$            SVAR2=0.D0
!!$            LN2=0
!!$!
!!$            IPRO1=0
!!$            DO LN1=1,LN-1
!!$              IF(LOX(LN1).NE.L) CYCLE
!!$              IPRO1=IPRO1+1
!!$              IF(TDOT(LN1)) CYCLE
!!$              IF(LN2.EQ.0) THEN
!!$                LN2=LN1
!!$                CYCLE
!!$              END IF
!!$              SVAR1=SVAR1*(EOFLN(LN)-EOFLN(LN2))
!!$              SVAR2=SVAR2+1.D0/(EOFLN(LN)-EOFLN(LN2))
!!$              SVAR=SVAR1*SVAR2
!!$              QN(:,LN) =QN(:,LN) +NLPHI(:,LN1)*SVAR
!!$              TQN(:,LN)=TQN(:,LN)+TNLPHI(:,LN1)*SVAR
!!$              LN2=LN1
!!$            ENDDO
!!$            QN(:,LN) = QN(:,LN)+ NLPHI(:,LN)*SVAR1
!!$            TQN(:,LN)=TQN(:,LN)+TNLPHI(:,LN)*SVAR1
!!$          ELSE
!!$            QN(:,LN)=0.D0
!!$            TQN(:,LN)=0.D0
!!$            SVAR=1.D0
!!$            DO LN1=1,LN
!!$              IF(LOX(LN1).NE.L) CYCLE
!!$              IF(TDOT(LN1)) CYCLE
!!$              QN(:,LN) =QN(:,LN) + NLPHI(:,LN1)*SVAR
!!$              TQN(:,LN)=TQN(:,LN)+TNLPHI(:,LN1)*SVAR
!!$              SVAR=SVAR*(EOFLN(LN)-EOFLN(LN1))
!!$            ENDDO
!!$          END IF
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == TEST EQUATION FOR QN                                                 ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20,"  TEST QN EQ.  ")')
!!$        DO L=0,LX
!!$          NC1=0
!!$          DO IB=NC,1,-1
!!$            IF(LOFI(IB).NE.L) CYCLE
!!$            NC1=IB
!!$            EXIT
!!$          ENDDO
!!$          IPRO=0
!!$          DO LN=1,LNX
!!$            IF(LOX(LN).NE.L) CYCLE
!!$            IPRO=IPRO+1
!!$            PRO(:,LN)=TQN(:,LN)+(AEPOT*Y0-EOFLN(LN))*QN(:,LN)
!!$            IF(TDOT(LN)) THEN
!!$              IF(IPRO.EQ.0) THEN
!!$                IF(NC1.NE.0) PRO(:,LN)=PRO(:,LN)-UOFI(:,NC1)
!!$                WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|QN>-|UC>     =",F20.15)') &
!!$     &                       LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
!!$              ELSE
!!$                PRO(:,LN)=PRO(:,LN)-QN(:,LN-1)
!!$                WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|QN>-|Q_{N-1}>=",F20.15)') &
!!$     &                       LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
!!$              END IF
!!$            ELSE
!!$              IF(NC1.NE.0) PRO(:,LN)=PRO(:,LN)-UOFI(:,NC1)
!!$              WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|QN>-|UC>     =",F20.15)') &
!!$     &                     LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
!!$            END IF
!!$          ENDDO
!!$        ENDDO
!!$      END IF
!!$!
!!$OPEN(UNIT=8,FILE='QN.DAT')
!!$DO IR=1,NR
!!$  IF(R(IR).GT.10.D0) EXIT
!!$  WRITE(8,'(30F20.5)')R(IR),MAX(-10.D0,MIN(10.D0,QN(IR,:))),MAX(-10.D0,MIN(10.D0,NLPHI(IR,:)))
!!$ENDDO
!!$CLOSE(8)
!!$
!!$OPEN(UNIT=8,FILE='BAREPRO-1.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PRO(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CORE-ORTHOGONALIZE TO OBTAIN NODAL PARTIAL WAVES                     ==
!!$!     == NOTE THAT THESE ARE NOT EIGENSTATES AS THEY ARE ONLY ORTHOGONALIZED  ==
!!$!     == TO THE CORE STATES                                                   ==
!!$!     ==========================================================================
!!$      AEPHI(:,:)=QN(:,:)    !NLPHI(:,:)
!!$      TAEPHI(:,:)=TQN(:,:)  !TNLPHI(:,:)
!!$      DO L=0,LX
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          DO IB=NC,1,-1
!!$            IF(LOFI(IB).NE.L) CYCLE
!!$            AUX(:)=R(:)**2*UOFI(:,IB)*AEPHI(:,LN)
!!$            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR1)
!!$            AUX(:)=R(:)**2*UOFI(:,IB)**2
!!$            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$            CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR2)
!!$            VAL=SVAR1/SVAR2
!!$            AEPHI(:,LN)=AEPHI(:,LN)-UOFI(:,IB)*VAL
!!$            TAEPHI(:,LN)=TAEPHI(:,LN)-TUOFI(:,IB)*VAL
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$OPEN(UNIT=8,FILE='AEPHI-1.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),AEPHI(IR,:),QN(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == TEST EQUATION FOR ALL-ELECTRON PARTIAL WAVES AEPHI                   ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20,"  TEST AEPHI EQ.  ")')
!!$        DO LN=1,LNX
!!$          PRO(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
!!$          IF(TDOT(LN)) THEN
!!$            DO LN1=LN-1,1,-1
!!$              IF(LOX(LN1).NE.LOX(LN)) CYCLE
!!$              PRO(:,LN)=PRO(:,LN)-AEPHI(:,LN1) 
!!$              EXIT
!!$            ENDDO
!!$            WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|AEPHI_N>-|AEPHI_N-1>=",F20.15)') &
!!$      &          LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
!!$          ELSE
!!$            WRITE(6,FMT='("LN=",I2," L=",I2,"  [T+V-E_N]|AEPHI_N>            =",F20.15)') &
!!$      &          LN,LOX(LN),MAXVAL(ABS(PRO(:,LN)))
!!$          END IF
!!$        ENDDO
!!$      END IF
!!$
!!$OPEN(UNIT=8,FILE='BAREPRO-2.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PRO(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT PSEUDO PARTIAL WAVES                                       ==
!!$!     ==========================================================================
!!$      PSPHI=QN
!!$      TPSPHI=TQN
!!$      DO L=0,LX
!!$        NPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          NPRO=NPRO+1
!!$          RC1=RC(LN)
!!$        ENDDO
!!$        IF(NPRO.EQ.0) CYCLE
!!$        ALLOCATE(PHITEST(NR,NPRO))
!!$        ALLOCATE(TPHITEST(NR,NPRO))
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$          PHITEST(:,IPRO)=QN(:,LN)
!!$          TPHITEST(:,IPRO)=TQN(:,LN)
!!$          IF(TDOT(LN).AND.IPRO.GT.1) THEN
!!$            PHITEST(:,IPRO) =PHITEST(:,IPRO-1) +0.1D0*PHITEST(:,IPRO)
!!$            TPHITEST(:,IPRO)=TPHITEST(:,IPRO-1)+0.1D0*TPHITEST(:,IPRO)
!!$          END IF                 
!!$        ENDDO
!!$        CALL ATOMIC_MAKEPSPHI(GID,NR,RC1,L,NPRO,PHITEST,TPHITEST)
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$          IF(TDOT(LN).AND.IPRO.GT.1) THEN
!!$            PHITEST(:,IPRO) =10.D0*(PHITEST(:,IPRO)-PHITEST(:,IPRO-1))
!!$            TPHITEST(:,IPRO)=10.D0*(TPHITEST(:,IPRO)-TPHITEST(:,IPRO-1))
!!$          END IF
!!$          PSPHI(:,LN)=PHITEST(:,IPRO)
!!$          TPSPHI(:,LN)=TPHITEST(:,IPRO)
!!$        ENDDO
!!$        DEALLOCATE(PHITEST)
!!$        DEALLOCATE(TPHITEST)
!!$      ENDDO
!!$
!!$OPEN(UNIT=8,FILE='PSPHI-BARE.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),QN(IR,:),PSPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='TPSPHI-BARE.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),TQN(IR,:),TPSPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == BACK TRANSFORM                                                       ==
!!$!     ==========================================================================
!!$      DO L=0,LX
!!$!       == DETERMINE NPRO ======================================================
!!$        NPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).EQ.L) NPRO=NPRO+1
!!$        ENDDO
!!$!
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$          IF(TDOT(LN)) CYCLE
!!$PRINT*,'MAKE ',LN,' TO BASE'
!!$          IPRO1=IPRO
!!$          DO LN1=LN+1,LNX
!!$            IF(LOX(LN1).NE.L) CYCLE
!!$            IPRO1=IPRO1+1
!!$            IF(IPRO1.LE.IPRO) CYCLE
!!$            IF(TDOT(LN1)) THEN
!!$              IF(IPRO1.EQ.IPRO+1) CYCLE
!!$!             == LN1 -> PHIDOT
!!$!             == LN2 -> PHI
!!$              DO LN2=LN1-1,1,-1
!!$                IF(LOX(LN2).NE.L) CYCLE
!!$                SVAR=1.D0/(EOFLN(LN1)-EOFLN(LN))
!!$                QN(:,LN1) =( QN(:,LN1) -QN(:,LN2))*SVAR
!!$                TQN(:,LN1)=(TQN(:,LN1)-TQN(:,LN2))*SVAR
!!$                PSPHI(:,LN1) =( PSPHI(:,LN1) -PSPHI(:,LN2))*SVAR
!!$                TPSPHI(:,LN1)=(TPSPHI(:,LN1)-TPSPHI(:,LN2))*SVAR
!!$                AEPHI(:,LN1) =( AEPHI(:,LN1) -AEPHI(:,LN2))*SVAR
!!$                TAEPHI(:,LN1)=(TAEPHI(:,LN1)-TAEPHI(:,LN2))*SVAR
!!$                EXIT
!!$              ENDDO
!!$            ELSE
!!$              SVAR=1.D0/(EOFLN(LN1)-EOFLN(LN))
!!$              QN(:,LN1)    =(    QN(:,LN1)-    QN(:,LN))*SVAR
!!$              TQN(:,LN1)   =(   TQN(:,LN1)-   TQN(:,LN))*SVAR
!!$              PSPHI(:,LN1) =( PSPHI(:,LN1)- PSPHI(:,LN))*SVAR
!!$              TPSPHI(:,LN1)=(TPSPHI(:,LN1)-TPSPHI(:,LN))*SVAR
!!$              AEPHI(:,LN1) =( AEPHI(:,LN1)- AEPHI(:,LN))*SVAR
!!$              TAEPHI(:,LN1)=(TAEPHI(:,LN1)-TAEPHI(:,LN))*SVAR
!!$            END IF
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     == TEST IF BACK TRANSFORM WAS SUCCESSFUL ================================
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20,"  TEST BACK TRANSFORM  ")')
!!$        DO LN=1,LNX
!!$          WRITE(6,FMT='("LN=",I2," L=",I2," DIFF. NDLSS PHI",F10.5 &
!!$     &                                 ," DIFF. KIN.OP NDLSS. PHI ",F10.5)') &
!!$     &          LN,LOX(LN),MAXVAL(ABS(QN(:,LN)-NLPHI(:,LN))) &
!!$     &                    ,MAXVAL(ABS(TQN(:,LN)-TNLPHI(:,LN)))
!!$        ENDDO
!!$      END IF
!!$ 
!!$OPEN(UNIT=8,FILE='PSPHI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PSPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='AEPHI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),AEPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='TPSPHI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),TNLPHI(IR,:),TPSPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CONSTRUCT PROJECTOR FUNCTIONS                                        ==
!!$!     ==========================================================================
!!$      DO LN=1,LNX
!!$!       == PRO(:,LN)=TNLPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*NLPHI(:,LN) ==========
!!$        BAREPRO(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN)
!!$        TCHK=.FALSE.
!!$        DO LN1=LN-1,1,-1
!!$          IF(LOX(LN1).NE.LOX(LN)) CYCLE
!!$          IF(TDOT(LN1)) CYCLE
!!$!         ==  PRO(:,LN)=PRO(:,LN)-NLPHI(:,LN1) =================================
!!$          BAREPRO(:,LN)=BAREPRO(:,LN)-PSPHI(:,LN1) 
!!$          TCHK=.TRUE.
!!$          EXIT
!!$        ENDDO
!!$        IF(.NOT.TCHK) THEN
!!$          DO IB=NC,1,-1
!!$            IF(LOFI(IB).NE.LOX(LN)) CYCLE
!!$!           == PRO(:,LN)=PRO(:,LN)-UOFI(:,IB) ==================================
!!$            EXIT
!!$          ENDDO
!!$        END IF
!!$      ENDDO
!!$
!!$
!!$OPEN(UNIT=8,FILE='BAREPRO.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),BAREPRO(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CHECK PAW EQUATION FOR PSEUDO PARTIALWAVES                           ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        ALLOCATE(PHITEST(NR,LNX))
!!$        ALLOCATE(TPHITEST(NR,LNX))
!!$        DO LN=1,LNX
!!$          TPHITEST(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN) &
!!$     &                  -BAREPRO(:,LN)
!!$          PHITEST(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
!!$          DO LN1=LN-1,1,-1
!!$            IF(LOX(LN1).NE.LOX(LN)) CYCLE
!!$            IF(TDOT(LN1)) CYCLE
!!$            PHITEST(:,LN)=PHITEST(:,LN)-AEPHI(:,LN1) 
!!$            TPHITEST(:,LN)=TPHITEST(:,LN)-PSPHI(:,LN1) 
!!$            EXIT
!!$          ENDDO
!!$        ENDDO
!!$!
!!$        WRITE(6,FMT='(82("="),T20,"  TEST RAW PAW EQUATION  ")')
!!$        DO LN=1,LNX
!!$          WRITE(6,FMT='("LN=",I2," L=",I2," RAW PAW EQ.",F10.5 &
!!$     &                                 ," SCHR. EQ.",F10.5)') &
!!$     &          LN,LOX(LN),MAXVAL(ABS(TPHITEST(:,LN))),MAXVAL(ABS(PHITEST(:,LN)))
!!$        ENDDO
!!$OPEN(UNIT=8,FILE='AUX.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),TPHITEST(IR,:),PHITEST(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$STOP
!!$        DEALLOCATE(PHITEST)
!!$        DEALLOCATE(TPHITEST)
!!$      ENDIF
!!$
!!$!
!!$!     ==========================================================================
!!$!     == CUT OFF THE TAILS OF THE PROJECTOR FUNCTIONS                         ==
!!$!     ==========================================================================
!!$      PRO(:,:)=BAREPRO(:,:)
!!$      DO LN=1,LNX
!!$        DO IR=1,NR
!!$          IF(R(IR).LT.RC(LN)) CYCLE
!!$          PRO(IR:,LN)=0.D0
!!$          EXIT
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == ENFORCE BIORTHOGONALIZATION                                          ==
!!$!     ==========================================================================
!!$      CALL BIORTHOMATRICES(GID,NR,RBOX,LNX,LOX,PSPHI,PRO,TRANSPHI,TRANSPRO)
!!$        WRITE(6,FMT='(82("="),T20,"  TRANSPRO ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F20.5)')TRANSPRO(LN1,:)
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  TRANSPHI ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F20.5)')TRANSPHI(LN1,:)
!!$        ENDDO
!!$
!!$!     ==  DETERMINE <PSPHI|PROBARE> ============================================
!!$      PSPHIPROBARE(:,:)=0.D0
!!$      DO LN1=1,LNX
!!$        DO LN2=1,LNX
!!$          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
!!$          AUX(:)=R(:)**2*PSPHI(:,LN1)*PRO(:,LN2)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          PSPHIPROBARE(LN1,LN2)=VAL
!!$        ENDDO
!!$      ENDDO
!!$!
!!$      DO L=0,LX
!!$        NPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          NPRO=NPRO+1
!!$        ENDDO
!!$        IF(NPRO.EQ.0) CYCLE
!!$        ALLOCATE(A(NPRO,NPRO))
!!$        ALLOCATE(AINV(NPRO,NPRO))
!!$        ALLOCATE(PRO1(NR,NPRO))
!!$        IPRO1=0
!!$        DO LN1=1,LNX
!!$          IF(LOX(LN1).NE.L) CYCLE
!!$          IPRO1=IPRO1+1
!!$          PRO1(:,IPRO1)=PRO(:,LN1)
!!$          IPRO2=0
!!$          DO LN2=1,LNX
!!$            IF(LOX(LN2).NE.L) CYCLE
!!$            IPRO2=IPRO2+1
!!$            A(IPRO1,IPRO2)=PSPHIPROBARE(LN1,LN2)
!!$          ENDDO
!!$        ENDDO
!!$        CALL LIB$INVERTR8(NPRO,A,AINV)
!!$DO IPRO1=1,NPRO
!!$  WRITE(*,FMT='("A",20F15.10)')A(IPRO1,:)
!!$ENDDO
!!$DO IPRO1=1,NPRO
!!$  WRITE(*,FMT='("AINV",20E15.3)')AINV(IPRO1,:)
!!$ENDDO
!!$A=MATMUL(AINV,A)
!!$DO IPRO1=1,NPRO
!!$  WRITE(*,FMT='("AINV*A",20F15.10)')A(IPRO1,:)
!!$ENDDO
!!$
!!$        A=MATMUL(A,AINV)
!!$        PRO1=MATMUL(PRO1,AINV)
!!$        DEALLOCATE(A)
!!$        DEALLOCATE(AINV)
!!$        IPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).NE.L) CYCLE
!!$          IPRO=IPRO+1
!!$          PRO(:,LN)=PRO1(:,IPRO)
!!$        ENDDO
!!$        DEALLOCATE(PRO1)
!!$      ENDDO
!!$OPEN(UNIT=8,FILE='PRO.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PRO(IR,:) !*R(IR)**2
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='AEPSPHI.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),AEPHI(IR,:),PSPHI(IR,:)
!!$ENDDO
!!$CLOSE(8)
!!$OPEN(UNIT=8,FILE='AEPSDIFF.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),(AEPHI(IR,:)-PSPHI(IR,:))*R(IR)**2
!!$ENDDO
!!$CLOSE(8)
!!$!
!!$!     ==========================================================================
!!$!     == CHECK BIORTHOGONALIZATION                                            ==
!!$!     ==========================================================================
!!$      DO LN1=1,LNX
!!$        DO LN2=1,LNX
!!$          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
!!$          AUX(:)=R(:)**2*PSPHI(:,LN1)*PRO(:,LN2)
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          IF(LN1.EQ.LN2)VAL=VAL-1.D0
!!$          IF(ABS(VAL).GT.1.D-5) THEN
!!$            CALL ERROR$MSG('BIORTHOGONALIZATION FAILED')
!!$            CALL ERROR$I4VAL('L',L)
!!$            CALL ERROR$I4VAL('LN1',LN1)
!!$            CALL ERROR$I4VAL('LN2',LN2)
!!$            CALL ERROR$R8VAL('DEVIATION',VAL)
!!$            CALL ERROR$STOP('XXXX')
!!$          END IF
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == DT,DO                                                                ==
!!$!     ==========================================================================
!!$      DT=0.D0
!!$      DO=0.D0
!!$      DH=0.D0
!!$      DO LN1=1,LNX
!!$        DO LN2=1,LNX
!!$          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
!!$          AUX(:)=R(:)**2*(AEPHI(:,LN1)*TAEPHI(:,LN2)-PSPHI(:,LN1)*TPSPHI(:,LN2))
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          DT(LN1,LN2)=VAL
!!$          AUX(:)=R(:)**2*(AEPHI(:,LN1)*AEPHI(:,LN2)-PSPHI(:,LN1)*PSPHI(:,LN2))
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          DO(LN1,LN2)=VAL
!!$          AUX(:)=R(:)**2*(AEPHI(:,LN1)*AEPOT(:)*Y0*AEPHI(:,LN2) &
!!$      &                  -PSPHI(:,LN1)*PSPOT(:)*Y0*PSPHI(:,LN2))
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,VAL)
!!$          DH(LN1,LN2)=DT(LN1,LN2)+VAL
!!$        ENDDO
!!$      ENDDO
!!$!
!!$      IF(TTEST) THEN
!!$        WRITE(6,FMT='(82("="),T20,"  DTKIN  ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')DT(LN1,:)
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  DOVERLAP  ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')DO(LN1,:)
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  DOHAMILTONIAN ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')DH(LN1,:)
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  <PSPHI|PRO-BARE> ")')
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')PSPHIPROBARE(LN1,:)
!!$        ENDDO
!!$!
!!$        WRITE(6,FMT='(82("="),T20,"  TEST <PSPHI|PRO-BARE>+DH-DO-DO_N-1=0 ")')
!!$        ALLOCATE(A(LNX,LNX))        
!!$        A=PSPHIPROBARE+DH
!!$        DO LN=1,LNX
!!$          A(:,LN)=A(:,LN)-DO(:,LN)*EOFLN(LN)
!!$        ENDDO
!!$        DO L=0,LX
!!$          LN1=0
!!$          DO LN=1,LNX
!!$            IF(LOX(LN).NE.L) CYCLE
!!$            IF(LN1.NE.0)A(:,LN)=A(:,LN)-DO(:,LN1)
!!$            LN1=LN
!!$          ENDDO
!!$        ENDDO
!!$        DO LN1=1,LNX
!!$          WRITE(6,FMT='(20F10.5)')A(LN1,:)
!!$        ENDDO
!!$OPEN(UNIT=8,FILE='AUX1.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),(SUM(PRO(IR,:)*A(:,LN)),LN=1,LNX)
!!$ENDDO
!!$CLOSE(8)
!!$!        DEALLOCATE(A)
!!$      END IF
!!$
!!$!     == SYMMETRIZE
!!$      DT=0.5D0*(DT+TRANSPOSE(DT))
!!$      DO=0.5D0*(DO+TRANSPOSE(DO))
!!$      DH=0.5D0*(DH+TRANSPOSE(DH))
!!$!
!!$!     ==========================================================================
!!$!     == CHECK2 PAW EQUATION FOR PSEUDO PARTIALWAVES                          ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        ALLOCATE(TPHITEST(NR,LNX))    ! HOLDS TEST FOR PSEUDO
!!$        ALLOCATE(PHITEST(NR,LNX))     ! HOLDS TEST FOR ALL-ELECTRON 
!!$        DO LN=1,LNX
!!$          TPHITEST(:,LN)=TPSPHI(:,LN)+(PSPOT(:)*Y0-EOFLN(LN))*PSPHI(:,LN)
!!$! TPHITEST(:,LN)=TPHITEST(:,LN)-BAREPRO(:,LN)
!!$          PHITEST1(:,LN)=0.D0          
!!$          DO LN1=1,LNX
!!$            IF(LOX(LN1).NE.LOX(LN)) CYCLE
!!$             TPHITEST(:,LN)=TPHITEST(:,LN)+PRO(:,LN1)*(DH(LN1,LN)-EOFLN(LN)*DO(LN1,LN))
!!$             PHITEST1(:,LN)=PHITEST1(:,LN)+PRO(:,LN1)*(DH(LN1,LN)-EOFLN(LN)*DO(LN1,LN))
!!$          ENDDO
!!$!
!!$          PHITEST(:,LN)=TAEPHI(:,LN)+(AEPOT(:)*Y0-EOFLN(LN))*AEPHI(:,LN)
!!$          DO LN1=LN-1,1,-1
!!$            IF(LOX(LN1).NE.LOX(LN)) CYCLE
!!$            IF(TDOT(LN1)) CYCLE
!!$!           ==  PRO(:,LN)=PRO(:,LN)-NLPHI(:,LN1) =================================
!!$            PHITEST(:,LN)=PHITEST(:,LN)-AEPHI(:,LN1) 
!!$            TPHITEST(:,LN)=TPHITEST(:,LN)-PSPHI(:,LN1) 
!!$            DO LN2=1,LNX
!!$              IF(LOX(LN2).NE.LOX(LN)) CYCLE
!!$              TPHITEST(:,LN)=TPHITEST(:,LN)-PRO(:,LN2)*DO(LN2,LN1)
!!$              PHITEST1(:,LN)=PHITEST1(:,LN)-PRO(:,LN2)*DO(LN2,LN1)
!!$            ENDDO
!!$            EXIT
!!$          ENDDO
!!$        ENDDO
!!$        WRITE(6,FMT='(82("="),T20,"  TEST PAW EQUATION  ")')
!!$        DO LN=1,LNX
!!$          WRITE(6,FMT='("LN=",I2," L=",I2," PAW EQ.",F10.5 &
!!$     &                                 ," SCHR. EQ.",F10.5)') &
!!$     &          LN,LOX(LN),MAXVAL(ABS(TPHITEST(:,LN))),MAXVAL(ABS(PHITEST(:,LN)))
!!$        ENDDO
!!$OPEN(UNIT=8,FILE='AUX.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PHITEST1(IR,:),-BAREPRO(IR,:),(SUM(PRO(IR,:)*A(:,LN)),LN=1,LNX)
!!$ENDDO
!!$CLOSE(8)
!!$STOP
!!$        DEALLOCATE(PHITEST)
!!$        DEALLOCATE(TPHITEST)
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     == UNSCREENING                                                          ==
!!$!     ==========================================================================
!!$      AERHO(:)=0.D0
!!$      PSRHO(:)=0.D0
!!$      DO L=0,LX
!!$        NPRO=0
!!$        DO LN=1,LNX
!!$          IF(LOX(LN).EQ.L) NPRO=NPRO+1
!!$        ENDDO
!!$        ALLOCATE(DH1(NPRO,NPRO))
!!$        ALLOCATE(DO1(NPRO,NPRO))
!!$        ALLOCATE(PRO1(NR,NPRO))
!!$        ALLOCATE(PROJ(NPRO))
!!$        IPRO1=0
!!$        DO LN1=1,LNX
!!$          IF(LOX(LN1).NE.L) CYCLE
!!$          IPRO1=IPRO1+1
!!$          PRO1(:,IPRO1)=PRO(:,LN1)
!!$          IPRO2=0
!!$          DO LN2=1,LNX
!!$            IF(LOX(LN2).NE.L) CYCLE
!!$            IPRO2=IPRO2+1
!!$            DH1(IPRO1,IPRO2)=DH(LN1,LN2)
!!$            DO1(IPRO1,IPRO2)=DO(LN1,LN2)
!!$          ENDDO
!!$        ENDDO
!!$        G(:)=0.D0
!!$        DO IB=NC+1,NB
!!$          IF(LOFI(IB).NE.L) CYCLE
!!$          E=EOFI(IB)
!!$          G(:)=0.D0
!!$          CALL ATOMLIB_PAWDER(GID,NR,L,E,PSPOT,NPRO,PRO1,DH1,DO1,G,PHI)
!!$
!!$          CALL SCHROEDINGER$SPHERICAL(GID,NR,AEPOT,DREL,0,G,L,E,1,AUX)
!!$          DO IR=1,NR
!!$            IF(R(IR).LT.RBOX) CYCLE
!!$            PHI(IR:)=0.D0
!!$          ENDDO
!!$PRINT*,'RBOX ',RBOX
!!$OPEN(UNIT=8,FILE='TEST.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),PHI(IR),AUX(IR)
!!$ENDDO
!!$CLOSE(8)
!!$STOP
!!$          DO IPRO=1,NPRO
!!$            AUX(:)=R(:)**2*PHI(:)*PRO1(:,IPRO)
!!$            CALL RADIAL$INTEGRAL(GID,NR,AUX,PROJ(NPRO))
!!$          ENDDO
!!$PRINT*,'PROJ ',L,PROJ
!!$          AUX(:)=R(:)**2*PHI(:)**2
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!!$PRINT*,'VAL1 ',L,VAL
!!$          DO IPRO1=1,NPRO
!!$            DO IPRO2=1,NPRO
!!$              VAL=VAL+PROJ(IPRO1)*DO1(IPRO1,IPRO2)*PROJ(IPRO2)
!!$            ENDDO
!!$          ENDDO
!!$PRINT*,'VAL2 ',L,VAL
!!$          PHI(:)=PHI(:)/SQRT(VAL)
!!$          PSRHO(:)=PSRHO(:)+FOFI(IB)*PHI(:)**2
!!$          AERHO(:)=AERHO(:)+FOFI(IB)*AEPSI(:,IB)**2
!!$        ENDDO
!!$        DEALLOCATE(DH1)
!!$        DEALLOCATE(DO1)
!!$        DEALLOCATE(PRO1)
!!$        DEALLOCATE(PROJ)
!!$      ENDDO      
!!$OPEN(UNIT=8,FILE='RHO.DAT')
!!$DO IR=1,NR
!!$  WRITE(8,'(30F20.5)')R(IR),AERHO(IR),PSRHO(IR)
!!$ENDDO
!!$CLOSE(8)
!!$STOP
!!$
!!$!      CALL ATOMLIB$BOXVOFRHO(GID,NR,RAD,AEZ,RHO,POT,EH,EXC)
!!$
!!$      RETURN
!!$      END
!
!     ......................................................................
      SUBROUTINE BIORTHOMATRICES(GID,NR,RBOX,LNX,LOX,PSPHI,PRO,TRANSPHI,TRANSPRO)
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
      INTEGER(4) ,INTENT(IN)     :: LNX
      INTEGER(4) ,INTENT(IN)     :: LOX(LNX)
      REAL(8)    ,INTENT(IN)     :: PSPHI(NR,LNX)
      REAL(8)    ,INTENT(IN)     :: PRO(NR,LNX)
      REAL(8)    ,INTENT(OUT)    :: TRANSPHI(LNX,LNX)
      REAL(8)    ,INTENT(OUT)    :: TRANSPRO(LNX,LNX)
      INTEGER(4)                 :: LN,LN1,LN2
      REAL(8)                    :: AUX(NR),AUX1(NR)
      REAL(8)                    :: R(NR)
      REAL(8)                    :: SVAR
      LOGICAL(4),PARAMETER       :: TTEST=.TRUE.
      REAL(8)                    :: PROPSI(LNX,LNX)
      REAL(8)                    :: TRANSPROPSI(LNX,LNX)
!     *********************************************************************
      CALL RADIAL$R(GID,NR,R)
!
!     =====================================================================
!     == CALCULATE INITIAL VIOLATION OF BIORTHOGONALITY                  ==
!     =====================================================================
      PROPSI(:,:)=0.D0
      DO LN1=1,LNX
        DO LN2=1,LNX
          IF(LOX(LN1).NE.LOX(LN2)) CYCLE
          AUX(:)=PRO(:,LN1)*PSPHI(:,LN2)*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)       
          CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)       
          PROPSI(LN1,LN2)=SVAR
        ENDDO
      ENDDO
!
!     =====================================================================
!     == COLLECT TRANSFORMATION MATRIX BETWEEN NEW AND OLD               ==
!     =====================================================================
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
!         == ORTHOGONALIZE PARTIAL WAVES ===================================
          SVAR=TRANSPROPSI(LN2,LN1)
          TRANSPHI(:,LN1)   =TRANSPHI(:,LN1)-TRANSPHI(:,LN2)*SVAR
          TRANSPROPSI(:,LN1)=TRANSPROPSI(:,LN1)-TRANSPROPSI(:,LN2)*SVAR          
!         == ORTHOGONALIZE PROJECTOR=====================================
          SVAR=TRANSPROPSI(LN1,LN2)
          TRANSPRO(:,LN1)   =TRANSPRO(:,LN1)-TRANSPRO(:,LN2)*SVAR
          TRANSPROPSI(LN1,:)=TRANSPROPSI(LN1,:)-TRANSPROPSI(LN2,:)*SVAR          
        ENDDO
        SVAR=TRANSPROPSI(LN1,LN1)
        TRANSPRO(:,LN1)=TRANSPRO(:,LN1)/SVAR
        TRANSPROPSI(LN1,:)=TRANSPROPSI(LN1,:)/SVAR
      ENDDO
!
!     =====================================================================
!     == CHECK RESULT                                                    ==
!     =====================================================================
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
        CALL ERROR$STOP('BIORTHOMATRICES')
      END IF
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
!CALL setup_WRITEPHI('TESTCONSTRUCTPSPHI1',GID,NR,NB,PSPHI)
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
          CALL ERROR$r8val('shift',shift)
          CALL ERROR$r8val('phase(1)',phase(1))
          CALL ERROR$r8val('offsetnn',offsetnn)
          CALL ERROR$r8val('rc',rc)
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
!     ......................................................................
      SUBROUTINE ATOMIC_PSEUDIZE(GID,NR,POW,VAL0,RC,AEF,PSF)
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      REAL(8)    ,INTENT(IN)     :: POW   !POWER
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
      A12=RC**(POW+2.D0)
      A21=POW       *RC**(POW-1.D0)
      A22=(POW+2.D0)*RC**(POW+1.D0)
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
        SVAR=C0+C1*RC**POW+C2*RC**(POW+2.D0)
        PRINT*,'TEST PESUDIZE 1 ',SVAR,VAL,SVAR-VAL
        SVAR=REAL(POW,8)*C1*RC**(POW-1.D0)+REAL(POW+2,8)*C2*RC**(POW+1.D0)
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
!     ..................................................................
      SUBROUTINE SETUP$AECORE(ISP_,NRX_,AECORE_)
!LEGACY CODE! -> SETUP$GETR8A('AECORE'
!     ******************************************************************
!     **  RETURN AE CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: AECORE_(NRX_)
!     ******************************************************************
      CALL ERROR$MSG('CODE IS MARKED FOR DELETION')
      CALL ERROR$MSG('USE SETUP$GETR8A INSTEAD')
      CALL ERROR$STOP('SETUP$AECORE')
      CALL SETUP$ISELECT(ISP_)
      AECORE_(:)=THIS%AECORE(:)
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$PSCORE(ISP_,NRX_,PSCORE_)
!LEGACY CODE! -> SETUP$GETR8A('PSCORE'
!     ******************************************************************
!     **  RETURN PS CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: PSCORE_(NRX_)
!     ******************************************************************
      CALL ERROR$MSG('CODE IS MARKED FOR DELETION')
      CALL ERROR$MSG('USE SETUP$GETR8A INSTEAD')
      CALL ERROR$STOP('SETUP$PSCORE')
      CALL SETUP$ISELECT(ISP_)
      PSCORE_(:)=THIS%PSCORE(:)
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$VBAR(ISP_,NRX_,VBAR_)
!LEGACY CODE! -> SETUP$GETR8A('VADD'
!     ******************************************************************
!     **  RETURN PS CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: VBAR_(NRX_)
!     ******************************************************************
      CALL ERROR$MSG('CODE IS MARKED FOR DELETION')
      CALL ERROR$MSG('USE SETUP$GETR8A INSTEAD')
      CALL ERROR$STOP('SETUP$VBAR')
      CALL SETUP$ISELECT(ISP_)
      VBAR_(:)=THIS%VADD(:)
      RETURN  
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_RENEWLOCORB()
!     **************************************************************************
!     **                                                                      **
!     **  DRIVER ROUTINE FOR SETUP_CHIFROMPHI                                 **
!     **                                                                      **
!     **  CHECKS WHETHER AN UPDATE IS NECESSARY; REALLOCATES ARRAYS AND       **
!     **  COMPUTES NEW TRANSFORMATION MATRICES                                **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ORBRAD
      USE PERIODICTABLE_MODULE
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: LN,L,I,NR,LXX
      REAL(8)   ,ALLOCATABLE :: RAD(:)
      CHARACTER(3)           :: ID ! CAN BE OLD OR NEW
      REAL(8)                :: RASA
!     **************************************************************************
      call error$msg('routine is marked for deletion')
      call error$stop('SETUP_renewlocorb')
      IF(THIS%LOCORBINI) RETURN
      THIS%LOCORBINI=.TRUE.
                            CALL TRACE$PUSH('SETUP_RENEWLOCORB')
!     == LOCORBINI IS RESET TO FALSE, WHEVEVER THE DATA ARE CHANGED VIA SETUP$SET
PRINT*,'A1',THIS%LOCORBNOFL
PRINT*,'A2',THIS%LOCORBLNX
!
!     == REALLOCATE DATA =======================================================
      IF(SUM(THIS%LOCORBNOFL).NE.THIS%LOCORBLNX) THEN
PRINT*,'B'
        THIS%LOCORBLNX=SUM(THIS%LOCORBNOFL)
        IF(ASSOCIATED(THIS%LOCORBLOX)) DEALLOCATE(THIS%LOCORBLOX)
        IF(ASSOCIATED(THIS%LOCORBAMAT)) DEALLOCATE(THIS%LOCORBAMAT)
        IF(ASSOCIATED(THIS%LOCORBBMAT)) DEALLOCATE(THIS%LOCORBBMAT)
        ALLOCATE(THIS%LOCORBLOX(THIS%LOCORBLNX))
        ALLOCATE(THIS%LOCORBAMAT(THIS%LNX,THIS%LOCORBLNX))
        ALLOCATE(THIS%LOCORBBMAT(THIS%LOCORBLNX,THIS%LNX))
        LN=0
        DO L=0,3
          DO I=1,THIS%LOCORBNOFL(L+1)
            LN=LN+1
            THIS%LOCORBLOX(LN)=L
          ENDDO
        ENDDO
      END IF
PRINT*,'C'
!
!     == RECALCULATE LOCAL ORBITALS ============================================
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
      ALLOCATE(RAD(THIS%LOCORBLNX))
      DO LN=1,THIS%LOCORBLNX
        L=THIS%LOCORBLOX(LN)
        RAD(LN)=THIS%LOCORBRAD(L+1)
      ENDDO
      CALL SETUP_CHIFROMPHI(THIS%GID,NR,THIS%LNX,THIS%LOX,THIS%AEPHI &
     &               ,THIS%LOCORBNOFL,RAD,THIS%LOCORBLNX &
     &               ,THIS%LNX,THIS%LOCORBLOX,THIS%LOCORBAMAT,THIS%LOCORBBMAT)
      DEALLOCATE(RAD)
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_CHIFROMPHI(GID,NR,LNX,LOX,PHI,NORB,RCUT &
     &                           ,LNXCHI,LNXPHI,LOXCHI,AMAT,BMAT)
!     **                                                                      **
!     **  DEFINES PROJECTOR FUNCTIONS FOR LOCAL ORBITALS IN TERMS OF          **
!     **  THE CONVENTIONAL PARTIAL WAVE PROJECTOR FUNCTIONS                   **
!     **                                                                      **
!     **     <P_CHI_I|= SUM_J BMAT(I,J) <P_PHI_J|                             **
!     **     |CHI_I>  = SUM_J |PHI_J> AMAT(J,I)                               **
!     **                                                                      **
!     **  THE LOCAL ORBITALS ARE SELECTED BY NORB WHICH SPECIFIES THE NUMBER  **
!     **  OF LOCAL ORBITALS PER ANGULAR MOMENTUM CHANNEL.                     **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID     ! GRID ID
      INTEGER(4),INTENT(IN) :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: LNX     ! #(PARTIAL WAVES) 
      INTEGER(4),INTENT(IN) :: LOX(LNX)! ANGULAR MOMENTUM PER PARTIAL WAVE
      REAL(8)   ,INTENT(IN) :: PHI(NR,LNX)    ! PARTIAL WAVES
      INTEGER(4),INTENT(IN) :: NORB(4) ! X#(HEAD ORBITALS PER L)
      REAL(8)   ,INTENT(IN) :: RCUT(LNX)  ! EXTENT OF HEAD FUNCTIONS
      INTEGER(4),INTENT(IN) :: LNXCHI  ! #(LOCAL ORBITALS)
      INTEGER(4),INTENT(IN) :: LNXPHI  ! #(PARTIAL WAVES)
      INTEGER(4),INTENT(OUT):: LOXCHI(LNXCHI)  ! ANGULAR MOMENTUM PER LOCAL ORB.
      REAL(8)   ,INTENT(OUT):: BMAT(LNXCHI,LNXPHI) ! MAPPING OF PROJECTORS
      REAL(8)   ,INTENT(OUT):: AMAT(LNXPHI,LNXCHI)     ! MAPPING OF WAVE FUNCTIONS
      REAL(8)   ,ALLOCATABLE:: CHI(:,:)
      REAL(8)   ,ALLOCATABLE:: CHI1(:,:)
      REAL(8)   ,ALLOCATABLE:: A(:,:)
      REAL(8)   ,ALLOCATABLE:: B(:,:)
      REAL(8)   ,ALLOCATABLE:: A1(:,:)
      REAL(8)   ,ALLOCATABLE:: MAT(:,:)
      REAL(8)   ,ALLOCATABLE:: MATINV(:,:)
      REAL(8)   ,ALLOCATABLE:: R(:)        
      REAL(8)   ,ALLOCATABLE:: G(:)        
      REAL(8)   ,PARAMETER  :: RCG=1.D-2
      REAL(8)   ,ALLOCATABLE:: AUX(:),AUX1(:)
      REAL(8)               :: SVAR1,SVAR2
      INTEGER(4)            :: IR
      INTEGER(4)            :: NX,N,LX,L,LN,LNCHI,NOFL,NOFL2,ISVAR
      INTEGER(4)            :: N1,N2,LN1,LN2,L1,L2
      CHARACTER(64)         :: FILE
!     **************************************************************************
                            CALL TRACE$PUSH('SETUP_CHIFROMPHI')
      call error$msg('routine is marked for deletion')
      call error$stop('SETUP_CHIFROMPHI')
!
!     ==========================================================================
!     ==  COLLECT DATA FROM SETUP OBJECT                                      ==
!     ==========================================================================
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  CONSISTENCY CHECKS                                                  ==
!     ==========================================================================
!     == CHECK IF #(LOCAL ORBITALS) IS CONSISTENT WITH ARRAY DIMENSIONS ========
      IF(SUM(NORB).NE.LNXCHI) THEN
        CALL ERROR$MSG('NORB AND LNXCHI INCONSISTENT')
        CALL ERROR$I4VAL('NORB',LNX)
        CALL ERROR$I4VAL('LNXCHI',LNXCHI)
        CALL ERROR$STOP('SETUP_CHIFROMPHI')
      END IF
!
!     == CHECK IF #(PARTIAL WAVES) FROM SETUP OBJECT IS CONSISTENT WITH ========
!     == ARRAY DIMENSIONS ======================================================
      IF(LNX.NE.LNXPHI) THEN
        CALL ERROR$MSG('LNX AND LNXPHI INCONSISTENT')
        CALL ERROR$I4VAL('LNX',LNX)
        CALL ERROR$I4VAL('LNXPHI',LNXPHI)
        CALL ERROR$STOP('SETUP_CHIFROMPHI')
      END IF
!
!     == CHECK IF NUMBER OF LOCAL ORBITALS DOES NOT EXCEED NUMBER OF PARTIAL WAVES
      LX=MAXVAL(LOX)
      DO L=0,LX
        NOFL=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE !CONSIDER ONLY PARTIAL WAVES WITH THE CORRECT L
          NOFL=NOFL+1
        ENDDO
        IF(NORB(L+1).GT.NOFL) THEN
          CALL ERROR$MSG('#(CORRELATED ORBITALS) EXCEEDS #(PARTIAL WAVES)')
          CALL ERROR$STOP('SETUP_CHIFROMPHI')
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  START WITH PARTIAL WAVES AS LOCAL ORBITALS                          ==
!     == |CHI_I>=SUM_I |PHI_J>A(J,I)
!     ==========================================================================
      ALLOCATE(CHI(NR,LNX))        
      CHI(:,:)=0.D0
      ALLOCATE(A(LNX,LNX))        
      A(:,:)=0.D0
      DO LN=1,LNX
        L=LOX(LN)
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
        A(LN,LN)=1.D0
        CHI(:,LN)=PHI(:,LN)
      ENDDO
!
!     ==========================================================================
!     == MAKE HEAD FUNCTION ANTIBONDING WITH NODE AT RCUT ======================
!     ==========================================================================
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
      DO L=0,LX
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
        NOFL=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          NOFL=NOFL+1
          IF(NOFL.GT.NORB(L+1)) CYCLE  ! CONSIDER ONLY HEAD FUNCTIONS
!         == SEARCH FOR NEXT PARTIAL WAVE WITH THIS L =========================
          LN1=0
          DO LN2=LN+1,LNX
            IF(LOX(LN2).NE.L) CYCLE  
            LN1=LN2
            EXIT
          ENDDO
          IF(LN1.EQ.0) THEN
            CALL ERROR$MSG('CAN NOT LOCALIZE')
            CALL ERROR$MSG('NO TAIL FUNCTION AVAILABLE')
            CALL ERROR$STOP('SETUP_CHIFROMPHI')
          END IF
!
!         == NOW CONTINUE; LN1 POINTS TO THE NEXT PHI WITH THIS L ==============
!         == IMPOSE NODE CONDITION==============================================
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN),RCUT(LN),SVAR1)
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN1),RCUT(LN),SVAR2)
          IF(SVAR1.EQ.0.D0.AND.SVAR2.EQ.0.D0) THEN
            CALL ERROR$MSG('PARTIAL WAVES ARE TRUNCATED INSIDE OF RCUT')
            CALL ERROR$MSG('THIS IS A FLAW OF THE IMPLEMENTATION')
            CALL ERROR$MSG('CHOOSE SMALLER RCUT')
            CALL ERROR$STOP('LDAPLUSU_CHIFROMPHI')
          END IF
          CHI(:,LN)=CHI(:,LN)*SVAR2-CHI(:,LN1)*SVAR1
          A(:,LN)=A(:,LN)*SVAR2-A(:,LN1)*SVAR1
!         == ORTHOGONALIZE TO THE LOWER HEAD FUNCTIONS =========================
          DO LN1=1,LN-1
            IF(LOX(LN1).NE.L) CYCLE  
            AUX(:)=CHI(:,LN)*CHI(:,LN1)*R(:)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,MIN(RCUT(LN),RCUT(LN1)),SVAR1)
            CHI(:,LN)=CHI(:,LN)-CHI(:,LN1)*SVAR1
            A(:,LN)=A(:,LN)-A(:,LN1)*SVAR1
          ENDDO
!         == NORMALIZE HEAD FUNCTION =============================================
          AUX(:)=CHI(:,LN)**2*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RCUT(LN),SVAR1)
          SVAR1=1.D0/SQRT(SVAR1)
          CHI(:,LN)=CHI(:,LN)*SVAR1
          A(:,LN)=A(:,LN)*SVAR1
        ENDDO
      END DO          
!
!     ==========================================================================
!     == DELOCALIZE TAIL FUNCTIONS                                            ==
!     ==========================================================================
!     == DEFINE WEIGHT FUNCTION       
      ALLOCATE(G(NR))
      G(:)=EXP(-(R(:)/RCG)**2)
!     ==  MINIMIZE SUCCESSIVELY THE WEIGHT OF <CHI2+X*CHI1|W|CHI2+X*CHI1> ======
!     ==  I.E. <CHI1|W|CHI2>+X<CHI1|W|CHI1> ====================================
      ALLOCATE(CHI1(NR,LNX))
      CHI1(:,:)=CHI(:,:)
      ALLOCATE(A1(LNX,LNX))        
      A1(:,:)=A(:,:)
      DO L=0,LX
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
        NOFL=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          NOFL=NOFL+1
!         == MINIMIZE CONTRIBUTION NEAR THE CENTER 
          NOFL2=0
          DO LN2=1,LN-1
            IF(LOX(LN2).NE.L) CYCLE
            NOFL2=NOFL2+1
            AUX(:)=G(:)*CHI1(:,LN)*CHI1(:,LN2)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
            AUX(:)=G(:)*CHI1(:,LN2)**2*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
            SVAR1=SVAR1/SVAR2
            CHI1(:,LN)=CHI1(:,LN)-CHI1(:,LN2)*SVAR1
            A1(:,LN)=A1(:,LN)-A1(:,LN2)*SVAR1
          ENDDO
        ENDDO
      ENDDO
!
!     = NOW MAP ONLY TAIL FUNCTIONS BACK AND LEAVE HEAD FUNCTIONS UNTOUCHED        
      DO L=0,LX
        NOFL=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          NOFL=NOFL+1
          IF(NOFL.LE.NORB(L+1)) CYCLE ! LEAVE HEAD FUNCTIONS ALONE
          CHI(:,LN)=CHI1(:,LN)
          A(:,LN)=A1(:,LN)
        ENDDO
      ENDDO
      DEALLOCATE(G)
      DEALLOCATE(CHI1)
      DEALLOCATE(A1)
      DEALLOCATE(AUX)
      DEALLOCATE(AUX1)
!
!     ==========================================================================
!     == NOW INVERT THE MATRIX A: B:=A^{-1}                                   ==
!     == |CHI_I>=SUM_J |PHI_J>A(J,I)                                          ==
!     == 1=|PHI><P|=|PHI>A B <P|=|CHI> ( B<P| )                               ==
!     == <P_CHI_I|=SUM_J B(I,J)<P_PHI_J|                                      ==
!     ==========================================================================
      ALLOCATE(B(LNX,LNX))
      B(:,:)=0.D0
      LX=MAXVAL(LOX)
      DO L=0,LX
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
!
        NX=0
        DO LN=1,LNX
          IF(LOX(LN).EQ.L) NX=NX+1
        ENDDO
        IF(NX.EQ.0) CYCLE   

        ALLOCATE(MAT(NX,NX))
        ALLOCATE(MATINV(NX,NX))
!        
!       == MAP A INTO MATRIX MAT ===============================================
        N1=0
        DO LN1=1,LNX
          L1=LOX(LN1)
          IF(L1.NE.L) CYCLE
          N1=N1+1
          N2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            IF(L2.NE.L) CYCLE
            N2=N2+1
            MAT(N2,N1)=A(LN2,LN1)
          ENDDO
        ENDDO
!
!       == INVERT MATRIX =======================================================
        CALL LIB$INVERTR8(NX,MAT,MATINV)
!
!       == MAP MATINV INTO B ===================================================
        N1=0
        DO LN1=1,LNX
          L1=LOX(LN1)
          IF(L1.NE.L) CYCLE
          N1=N1+1
          N2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            IF(L2.NE.L) CYCLE
            N2=N2+1
            B(LN1,LN2)=MATINV(N1,N2) 
          ENDDO
        ENDDO
        DEALLOCATE(MAT)
        DEALLOCATE(MATINV)
      ENDDO
!
!     ==========================================================================
!     === REMOVE TAIL FUNCTIONS                                               ==
!     ==========================================================================
      LNCHI=0
      DO L=0,LX
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
        NOFL=0
        DO LN=1,LNX
          IF(L.NE.LOX(LN)) CYCLE
          NOFL=NOFL+1
          IF(NOFL.GT.NORB(L+1)) EXIT
          LNCHI=LNCHI+1
          IF(LNCHI.GT.LNXCHI) THEN
            CALL ERROR$MSG('LNCHI OUT OF RANGE')
            CALL ERROR$I4VAL('LNCHI',LNCHI)
            CALL ERROR$I4VAL('LNXCHI',LNXCHI)
            CALL ERROR$STOP('SETUP_CHIFROMPHI')
          END IF
          LOXCHI(LNCHI)=L
          BMAT(LNCHI,:)=B(LN,:)
          AMAT(:,LNCHI)=A(:,LN)
        END DO
      ENDDO
      IF(LNCHI.NE.LNXCHI) THEN
        CALL ERROR$MSG('LNCHI AND LNXCHI ARE INCONSISTENT')
        CALL ERROR$I4VAL('LNCHI',LNCHI)
        CALL ERROR$I4VAL('LNXCHI',LNXCHI)
        CALL ERROR$STOP('SETUP_CHIFROMPHI')
      END IF
!
!     ==========================================================================
!     ==  CLEAN UP                                                            ==
!     ==========================================================================
      DEALLOCATE(CHI)
      DEALLOCATE(A)
      DEALLOCATE(B)
      DEALLOCATE(R)
                            CALL TRACE$POP()
      RETURN
      END
