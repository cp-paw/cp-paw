!=======================================================================
!=======================================================================
!=======================================================================
!====      FORMFACTORS                                              ====
!=======================================================================
!=======================================================================
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
!**    GETr8A('AMATCHI',LNXPHI*LNXCHI,AMAT)                                   **
!**    GETr8A('BMATCHI',LNXCHI*LNXPHI,BMAT)                                   **
!**                                                                           **
!**                                              P.E. BLOECHL, (1991-2008)    **
!*******************************************************************************
TYPE THIS_TYPE
INTEGER(4)             :: I            ! INTEGER IDENTIFIER (ISPECIES)
CHARACTER(32)          :: ID           ! IDENTIFIER (SPECIES-NAME)
INTEGER(4)             :: GID          ! GRID ID FOR R-SPACE GRID
INTEGER(4)             :: GIDG         ! GRID ID FOR G-SPACE GRID
REAL(8)                :: AEZ          ! ATOMIC NUMBER
REAL(8)                :: PSZ
REAL(8)                :: RCBG
REAL(8)                :: RCSM         ! GAUSSIAN DECAY FOR COMPENSATION CHARGE
INTEGER(4)             :: LNX          ! #(ORBITAL SHELLS)
INTEGER(4)             :: LMNX
INTEGER(4)             :: LMRX         ! #(ANGULAR MOMENTA FOR 1C-DENSITY)
INTEGER(4)             :: NC     !SANTOS040616
INTEGER(4),POINTER     :: LOX(:)       !(LNXX) MAIN ANGULAR MOMENTA 
REAL(8)   ,POINTER     :: VADD(:)      !(NRX)
REAL(8)   ,POINTER     :: AECORE(:)    !(NRX)  CORE ELECTRON DENSITY
REAL(8)   ,POINTER     :: PSCORE(:)    !(NRX)  PSEUDIZED ELECTRON DENSITY
REAL(8)   ,POINTER     :: PRO(:,:)     !(NRX,LNXX)  PROJECTOR FUNCTIONS
REAL(8)   ,POINTER     :: AEPHI(:,:)   !(NRX,LNXX)  AE PARTIAL WAVES
REAL(8)   ,POINTER     :: PSPHI(:,:)   !(NRX,LNXX)  PS PARTIAL WAVES
REAL(8)   ,POINTER     :: UPHI(:,:)    !(NRX,LNXX)   NODELESS PARTIAL WAVES
REAL(8)   ,POINTER     :: TUPHI(:,:)   !(NRX,LNXX)  KINETIC ENERGY OF ABOVE
REAL(8)   ,POINTER     :: DTKIN(:,:)   !(LNXX,LNXX) 1C-KIN. EN. MATRIX ELEMENTS
REAL(8)   ,POINTER     :: DOVER(:,:)   !(LNXX,LNXX) 1C-OVERLAP MATRIX ELEMENTS
REAL(8)   ,POINTER     :: VADDOFG(:)   !(NGX)
REAL(8)   ,POINTER     :: PSCOREOFG(:) !(NGX)
REAL(8)   ,POINTER     :: VHATOFG(:)   !(NGX)
REAL(8)   ,POINTER     :: NHATPRIMEOFG(:)  !(NGX)
REAL(8)   ,POINTER     :: PROOFG(:,:)  !(NRX,LNXX)  
LOGICAL(4)             :: LOCORBINI       ! LOCAL ORBITALS ARE INITIALIZED IF TRUE
REAL(8)                :: LOCORBRAD=5.D0  ! RADIUS OF LOCAL ORBITALS
INTEGER(4)             :: LOCORBNOFL(4)=0 ! #(LOCAL ORBITALS PER L)
INTEGER(4)             :: LOCORBLNX=0     ! #(LOCAL ORBITAL-SHELLS)
INTEGER(4),POINTER     :: LOCORBLOX(:)    ! L FOR EACH LOCAL ORBITAL-SHELL
REAL(8)   ,POINTER     :: LOCORBAMAT(:,:) ! |CHI>=|PHI>*AMAT
REAL(8)   ,POINTER     :: LOCORBBMAT(:,:) ! |PSI>=|CHI>BMAT<PTILDE|PSITILDE>
! SANTOS040616 BEGIN
REAL(8)   ,POINTER     :: AEPOT(:)     !(NRX)
INTEGER(4),POINTER     :: LB(:)        !(NC)
REAL(8)   ,POINTER     :: FB(:)        !(NC)
REAL(8)   ,POINTER     :: EB(:)        !(NC)
REAL(8)   ,POINTER     :: AEPSI(:,:)   !(NRX,NC)
! SANTOS040616 END
REAL(8)                :: M
REAL(8)                :: ZV
REAL(8)                :: PSG2
REAL(8)                :: PSG4
CHARACTER(32)          :: SOFTCORETYPE
CHARACTER(16)          :: FILEID
TYPE(THIS_TYPE),POINTER:: NEXT
END TYPE THIS_TYPE
!
!INTEGER(4)             :: NR    !=250
!INTEGER(4)             :: NRX=250
!REAL(8)   ,PARAMETER   :: GMAX=30     ! EPW[RY]<GMAX**2 FOR PSI AND RHO
!INTEGER(4),PARAMETER   :: NG=256
!REAL(8)               :: G1          ! FIRST POINT ON THE RADIAL G-GRID
INTEGER(4)             :: NSP=0
INTEGER(4)             :: LMRXX=0
INTEGER(4)             :: LMNXX=0
INTEGER(4)             :: LNXX=0
TYPE(THIS_TYPE),POINTER :: FIRST
TYPE(THIS_TYPE),POINTER :: THIS
TYPE FASTACCESS_TYPE  
  TYPE(THIS_TYPE),POINTER :: THIS
END TYPE FASTACCESS_TYPE
TYPE(FASTACCESS_TYPE),ALLOCATABLE :: FASTACCESS(:)
END MODULE SETUP_MODULE
!
!     ..................................................................
      SUBROUTINE SETUP$ISELECT(I)
!     ******************************************************************
!     **  SELECTS A SETUP PER INTEGER INDEX                           **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: I
      INTEGER(4)            :: J
!     ******************************************************************
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
      IF(I.GT.NSP) THEN
        CALL ERROR$MSG('INDEX I OUT OF RANGE')
        CALL ERROR$I4VAL('I',I)
        CALL ERROR$I4VAL('NSP',NSP)
        CALL ERROR$STOP('SETUP$ISELECT')
      END IF
!
      THIS=>FASTACCESS(I)%THIS
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$SELECT(ID)
!     ******************************************************************
!     **  SELECTS A SETUP PER ID                                      **
!     **  AND CREATES A NEW, IF IT DOES NOT EXIST                     **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
!     ******************************************************************
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
      THIS%PSZ   =0.D0
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
      NULLIFY(THIS%AEPOT)   !(NRX)
      NULLIFY(THIS%LB)      !(NC)
      NULLIFY(THIS%FB)      !(NC)
      NULLIFY(THIS%EB)      !(NC)
      NULLIFY(THIS%AEPSI)   !(NRX,NC)
! SANTOS040616 END
      THIS%LOCORBINI=.FALSE. ! INITIALIZED?
      THIS%LOCORBRAD=5.D0  ! RADIUS OF LOCAL ORBITALS
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
!     ..................................................................
      SUBROUTINE SETUP$GETCH(ID,VAL)
!     ******************************************************************
!     **  COLLECTS INTERNAL DATA                                      **
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      CHARACTER(*),INTENT(OUT) :: VAL
!     ******************************************************************
      IF(ID.EQ.'ID') THEN
        VAL=THIS%ID
PRINT*,'FROM SETUP$GETCH ',THIS%ID
PRINT*,'FROM SETUP$GETCH ',VAL
      ELSE IF(ID.EQ.'SOFTCORETYPE') THEN
        VAL=THIS%SOFTCORETYPE
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('SETUP$GETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$GETI4(ID,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(OUT) :: VAL
!     ******************************************************************
      IF(ID.EQ.'GID') THEN
        VAL=THIS%GID
      ELSE IF(ID.EQ.'GIDG') THEN
        VAL=THIS%GIDG
      ELSE IF(ID.EQ.'NR') THEN
!       == THIS SHALL BE MARKED FOR DELETION'
        CALL RADIAL$GETI4(THIS%GID,'NR',VAL)
!       VAL=NR
      ELSE IF(ID.EQ.'NRX') THEN
!       == THIS SHALL BE MARKED FOR DELETION'
        CALL RADIAL$GETI4(THIS%GID,'NR',VAL)
!       VAL=NRX
      ELSE IF(ID.EQ.'LNX') THEN
        VAL=THIS%LNX
      ELSE IF(ID.EQ.'LMNX') THEN
        VAL=THIS%LMNX
      ELSE IF(ID.EQ.'LMRX') THEN
        VAL=THIS%LMRX
      ELSE IF(ID.EQ.'LMRXX') THEN
        VAL=LMRXX
      ELSE IF(ID.EQ.'LNXCHI') THEN
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
          CALL ERROR$STOP('SETUP$GETI4A')
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
! SANTOS040616 BEGIN
      ELSE IF(ID.EQ.'LB') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETRI4A')
        END IF
        VAL=THIS%LB
! SANTOS040616 END
      ELSE IF(ID.EQ.'LOFC') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETRI4A')
        END IF
        VAL=THIS%LB
      ELSE IF(ID.EQ.'LOXCHI') THEN
!       == ANULAR MOMENTA OF LOCAL ORBITALS ==========================
        CALL SETUP_RENEWLOCORB()
        IF((LEN.NE.THIS%LOCORBLNX).OR.(THIS%LOCORBLNX.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETRI4A')
        END IF
        VAL=THIS%LOCORBLOX
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
!       == EXTENT OF LOCAL ORBITALS
        THIS%LOCORBINI=THIS%LOCORBINI.AND.(THIS%LOCORBRAD.EQ.VAL)
        THIS%LOCORBRAD=VAL
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
!     ..................................................................
      SUBROUTINE SETUP$GETR8A(ID,LEN,VAL)
!     ******************************************************************
!     **                                                              **
!     **  REMARK: REQUIRES PROPER SETUP TO BE SELECTED                **
!     **                                                              **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      REAL(8)     ,INTENT(OUT) :: VAL(LEN)
      INTEGER(4)               :: NR
!     ******************************************************************
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
      IF(ID.EQ.'PRO') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PRO,(/LEN/))
      ELSE IF(ID.EQ.'AEPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%AEPHI,(/LEN/))
      ELSE IF(ID.EQ.'PSPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%PSPHI,(/LEN/))
      ELSE IF(ID.EQ.'NDLSPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%UPHI,(/LEN/))
      ELSE IF(ID.EQ.'NDLSTPHI') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%TUPHI,(/LEN/))
      ELSE IF(ID.EQ.'AECORE') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%AECORE
      ELSE IF(ID.EQ.'PSCORE') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%PSCORE
      ELSE IF(ID.EQ.'VADD') THEN
        IF(LEN.NE.THIS%LNX*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=THIS%VADD
      ELSE IF(ID.EQ.'DEKIN') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%DTKIN,(/LEN/))
      ELSE IF(ID.EQ.'DO') THEN
        IF(LEN.NE.THIS%LNX**2) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%DOVER,(/LEN/))
! SANTOS040616 BEGIN
      ELSE IF(ID.EQ.'FB') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%FB,(/LEN/))
      ELSE IF(ID.EQ.'EB') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%EB,(/LEN/))
      ELSE IF(ID.EQ.'AECOREPSI') THEN
        IF(LEN.NE.THIS%NC*NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%AEPSI,(/LEN/))
      ELSE IF(ID.EQ.'ATOMICAEPOT') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL(:)=THIS%AEPOT(:)
! SANTOS040616 END
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
      ELSE IF(ID.EQ.'EOFC') THEN
        IF((LEN.NE.THIS%NC).OR.(THIS%NC.EQ.0)) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%EB,(/LEN/))
      ELSE IF(ID.EQ.'NUCPOT') THEN
        IF(LEN.NE.NR) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        CALL SETUP_NUCPOT(THIS%GID,NR,THIS%AEZ,VAL)
      ELSE IF(ID.EQ.'AMATCHI') THEN
        CALL SETUP_RENEWLOCORB()
        IF(LEN.NE.THIS%LOCORBLNX*THIS%LNX) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%LOCORBAMAT,(/LEN/))
      ELSE IF(ID.EQ.'BMATCHI') THEN
        CALL SETUP_RENEWLOCORB()
        IF(LEN.NE.THIS%LOCORBLNX*THIS%LNX) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SETUP$GETR8A')
        END IF
        VAL=RESHAPE(THIS%LOCORBBMAT,(/LEN/))
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SETUP$GETR8A')
      END IF
      RETURN
      END  
!
!     ..................................................................
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
!     ..................................................................
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
!     ******************************************************************
                            CALL TRACE$PUSH('SETUP_READ')
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
        CALL SETUPREAD$GETR8('PSZ',THIS%PSZ)
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
        CALL SETUPREAD$GETI4A('LOFC',NC,THIS%LB)
        CALL SETUPREAD$GETR8A('FOFC',NC,THIS%FB)
        CALL SETUPREAD$GETR8A('EOFC',NC,THIS%EB)
        THIS%AEPSI=0.D0
      ELSE
        CALL INPOT$READALL(NFIL,NRX,R1,DEX,NR,THIS%LNX,THIS%LOX &
     &         ,THIS%AEZ,THIS%PSZ,THIS%PSPHI,THIS%AEPHI &
     &         ,THIS%VADD,THIS%RCSM,THIS%DTKIN,THIS%DOVER &
     &         ,THIS%AECORE,THIS%PSCORE,THIS%PRO &
     &         ,THIS%AEPOT,THIS%NC,THIS%LB,THIS%FB,THIS%EB,THIS%AEPSI) !SANTOS040616
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
PRINT*,'PSZ ',THIS%PSZ
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
!!$PRINT*,'PSZ  ',THIS%PSZ
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
!     ..................................................................
      SUBROUTINE SETUP$REPORT(NFIL)
!     ******************************************************************
!     **  REPORT SETUP INFORMATION                                    **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4)     ,INTENT(IN) :: NFIL
      INTEGER(4)                 :: L,NPRO,LN,NPROSUM
      TYPE(THIS_TYPE),POINTER    :: THIS1
      CHARACTER(32)              :: STRING
      REAL(8)                    :: U
      INTEGER(4)                 :: THISTASK,NTASKS
!     ******************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
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
        IF(.NOT.ASSOCIATED(THIS1%NEXT)) EXIT
        THIS1=>THIS1%NEXT
      ENDDO
!
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE SETUP_NUCPOT(GID,NR,Z,POT)
!     **                                                                  **
!     **  ELECTROSTATIC POTENTIAL OF A NUCLEUS WITH FINITE RADIUS         **
!     **  THE NUCLEUS IS CONSIDERED AS A HOMOGENEOUSLY CHARGED SPHERE.    **
!     **  THE RADIUS IS RELATED  TO THE TOTAL MASS (NUMBER OF NUCLEONS),  **
!     **  WHICH CAN BE LOOKED UP KNOWING THE ATOMIC NUMBER Z.             **
!     **                                                                  **
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID      ! GRID ID
      INTEGER(4),INTENT(IN) :: NR       ! #(RADIAL GRID POINTS)
      REAL(8)   ,INTENT(IN) :: Z        ! ATOMIC NUMBER
      REAL(8)   ,INTENT(OUT):: POT(NR)  ! POTENTIAL
      REAL(8)               :: RNUC     ! NUCLEAR RADIUS
      REAL(8)               :: R(NR)    ! RADIAL GRID
      REAL(8)               :: PI,Y0
      INTEGER(4)            :: IR
!     ***********************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL PERIODICTABLE$GET(Z,'RNUC',RNUC)
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IF(R(IR).GT.RNUC) THEN
           POT(IR)=-Z/R(IR)
        ELSE
          POT(IR)=-Z/RNUC*(1.5D0-0.5D0*(R(IR)/RNUC)**2)
        END IF
      ENDDO
      POT(:)=POT(:)/Y0
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE SETUP$READ()
!     ******************************************************************
!     **  READ SETUP                                                  **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4)           :: ISP,NSP1
      CHARACTER(32)        :: NAME
!     ******************************************************************
                            CALL TRACE$PUSH('SETUP$READ')
!
!     ==================================================================
!     ==  CREATE SETUPS                                               ==
!     ==================================================================
      CALL ATOMTYPELIST$LENGTH(NSP1)
      DO ISP=1,NSP1
        CALL ATOMTYPELIST$NAME(ISP,NAME)
        CALL ATOMTYPELIST$SELECT(NAME)
        CALL ATOMTYPELIST$UNSELECT
        CALL SETUP$SELECT(NAME)
      ENDDO
!
!     ==================================================================
!     ==  READ SETUP FILES                                            ==
!     ==================================================================
      DO ISP=1,NSP1
        CALL SETUP$ISELECT(ISP)
        CALL SETUP_READ
      ENDDO
                            CALL TRACE$POP
      RETURN
      END

!
!     ..................................................................
      SUBROUTINE SETUP$AEPARTIALWAVES(ISP_,NRX_,LNX_,AEPHI_)
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
      SUBROUTINE SETUP$AECORE(ISP_,NRX_,AECORE_)
!     ******************************************************************
!     **  RETURN AE CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: AECORE_(NRX_)
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      AECORE_(:)=THIS%AECORE(:)
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$PSCORE(ISP_,NRX_,PSCORE_)
!     ******************************************************************
!     **  RETURN PS CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: PSCORE_(NRX_)
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      PSCORE_(:)=THIS%PSCORE(:)
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$VBAR(ISP_,NRX_,VBAR_)
!     ******************************************************************
!     **  RETURN PS CORE DNSITY                                       **
!     ******************************************************************
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_
      INTEGER(4),INTENT(IN) :: NRX_
      REAL(8)   ,INTENT(OUT):: VBAR_(NRX_)
!     ******************************************************************
      CALL SETUP$ISELECT(ISP_)
      VBAR_(:)=THIS%VADD(:)
      RETURN  
      END
!
!     ..................................................................
      SUBROUTINE SETUP$1COVERLAP(ISP_,LNXX_,DOVER_)
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
     &         ,AEPOT,NC,LB,FB,EB,AEPSI)           !SANTOS040616/BLO
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
1000  CONTINUE
! SANTOS040616 END
                              CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SETUP_RENEWLOCORB()
!     **                                                                      **
!     **  DRIVER ROUTINE FOR SETUP_CHIFROMPHI                                 **
!     **                                                                      **
!     **  CHECKS WHETHER AN UPDATE IS NECESSARY; REALLOCATES ARRAYS AND       **
!     **  COMPUTES NEW TRANSFORMATION MATRICES                                **
!     **                                                                      **
      USE SETUP_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: LN,L,I,NR
!     **************************************************************************
      IF(THIS%LOCORBINI) RETURN
print*,'a1',THIS%LOCORBNOFL
print*,'a2',THIS%LOCORBlnx
!
!     == REALLOCATE DATA =======================================================
      IF(SUM(THIS%LOCORBNOFL).NE.THIS%LOCORBLNX) THEN
print*,'b'
        THIS%LOCORBLNX=SUM(THIS%LOCORBNOFL)
        IF(ASSOCIATED(THIS%LOCORBLOX)) DEALLOCATE(THIS%LOCORBLOX)
        IF(ASSOCIATED(THIS%LOCORBAMAT)) DEALLOCATE(THIS%LOCORBAMAT)
        IF(ASSOCIATED(THIS%LOCORBBMAT)) DEALLOCATE(THIS%LOCORBBMAT)
        ALLOCATE(THIS%LOCORBLOX(THIS%LOCORBLNX))
        ALLOCATE(THIS%LOCORBAMAT(THIS%LNX,THIS%LOCORBLNX))
        ALLOCATE(THIS%LOCORBBMAT(THIS%LOCORBLNX,THIS%LNX))
      END IF
print*,'c'
!
!     == RECALCULATE LOCAL ORBITALS ============================================
      CALL RADIAL$GETI4(THIS%GID,'NR',NR)
      CALL SETUP_CHIFROMPHI(THIS%GID,NR,this%lnx,THIS%LOX,THIS%AEPHI &
     &                 ,THIS%LOCORBNOFL,THIS%LOCORBRAD,THIS%LOCORBLNX &
     &                 ,THIS%LNX,THIS%LOCORBLOX,THIS%LOCORBAMAT,THIS%LOCORBBMAT)
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
      REAL(8)   ,INTENT(IN) :: RCUT    ! EXTENT OF HEAD FUNCTIONS
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
                            CALL TRACE$PUSH('setup_CHIFROMPHI')
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
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN),RCUT,SVAR1)
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN1),RCUT,SVAR2)
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
            CALL RADIAL$VALUE(GID,NR,AUX1,RCUT,SVAR1)
            CHI(:,LN)=CHI(:,LN)-CHI(:,LN1)*SVAR1
            A(:,LN)=A(:,LN)-A(:,LN1)*SVAR1
          ENDDO
!         == NORMALIZE HEAD FUNCTION =============================================
          AUX(:)=CHI(:,LN)**2*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RCUT,SVAR1)
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
