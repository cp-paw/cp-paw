MODULE LMTO_MODULE
TYPE HYBRIDSETTING_TYPE
  !== ATOM-SPECIFIC SETTINGS ==
  LOGICAL(4)  :: ACTIVE=.FALSE.    ! CONSIDER HYBRID CONTRIBUTION ON THIS ATOM
  LOGICAL(4)  :: TCV        ! INCLUDE CORE VALENCE EXCHANGE
  LOGICAL(4)  :: TNDDO      ! INCLUDE NDDO OFFSITE EXCHANGE TERMS
  LOGICAL(4)  :: T31        ! INCLUDE 31 OFFSITE EXCHANGE TERMS 
  LOGICAL(4)  :: TBONDX     ! INCLUDE BOND EXCHANGE TERMS 
  LOGICAL(4)  :: TFOCKSETUP ! CALCULATE ATOM WITH FOCK TERM
  REAL(8)     :: LHFWEIGHT  ! LOCAL EXCHANGE TERMS ARE TREATED WITH
                            ! LHFWEIGHT, EXCEPT WHEN IT IS NEGATIVE
!  REAL(8)     :: K2
  REAL(8)     :: TAILEDLAMBDA1  ! LARGER DECAY CONSTANT FOR THE NTBO TAILS
  REAL(8)     :: TAILEDLAMBDA2  ! SMALLER DECAY CONSTANT FOR THE NTBO TAILS
!  REAL(8)     :: RANGESCALE ! DETERMINES RANGE OF NEAREST NEIGHBOR LIST
END TYPE HYBRIDSETTING_TYPE
!
TYPE ORBITALGAUSSCOEFF_TYPE
  INTEGER(4)         :: NIJK     !CAN ALSO BE NPOW FOR RADIAL FUNCTION
  INTEGER(4)         :: NE
  INTEGER(4)         :: NORB
  REAL(8)   ,POINTER :: E(:)     !(NE)
  REAL(8)   ,POINTER :: C(:,:,:) !(NIJK,NE,NORB)
END TYPE ORBITALGAUSSCOEFF_TYPE
!
TYPE TAILED_TYPE
  INTEGER(4)         :: GID
  INTEGER(4)         :: LNX
  INTEGER(4)         :: LMNX
  INTEGER(4),POINTER :: LOX(:)           ! (LNX)
  INTEGER(4),POINTER :: LNDOT(:)         ! (LNX)
  INTEGER(4),POINTER :: LMNDOT(:)        ! (LMNX)
  REAL(8)   ,POINTER :: AEF(:,:)         ! (NR,LNX)
  REAL(8)   ,POINTER :: PSF(:,:)         ! (NR,LNX)
  REAL(8)   ,POINTER :: NLF(:,:)         ! (NR,LNX)
  REAL(8)   ,POINTER :: U(:,:,:,:)       ! (LMNX,LMNX,LMNX,LMNX)
  REAL(8)   ,POINTER :: OVERLAP(:,:) !(LMNX,LMNX) OVERLAP MATRIX ELEMENTS
  REAL(8)   ,POINTER :: QLN(:,:,:)   !(2,LNX,LNX) MONO- AND DIPOLE MATRIX ELEMENTS
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSNLF
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: PRODRHO
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: PRODPOT
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: TRIPLE
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: SINGLE
END TYPE TAILED_TYPE

TYPE ORBITALSPHHARM_TYPE
  INTEGER(4)         :: GID
  INTEGER(4)         :: NR
  INTEGER(4)         :: LMX
  INTEGER(4)         :: NORB
  REAL(8)   ,POINTER :: F(:,:,:) !(NR,LM,IORB)
END TYPE ORBITALSPHHARM_TYPE
!
!== POTPARRED CONSIDERS ONLY ONE ANGULAR MOMENTUM CHANNEL PER LM ===============
!== CONSISTENT WITH THE SCREENED STRUCTURE CONSTANTS ===========================
!!$TYPE POTPARRED_TYPE
!!$  REAL(8)   ,POINTER :: DOVERLAPKK(:,:)
!!$  REAL(8)   ,POINTER :: DOVERLAPKJ(:,:)
!!$  REAL(8)   ,POINTER :: DOVERLAPJJ(:,:)
!!$END TYPE POTPARRED_TYPE
!
!== HOLDS THE POTENTIAL PARAMETER FOR ONE ATOM TYPE ============================
TYPE POTPAR_TYPE
  REAL(8)            :: RAD
  REAL(8)   ,POINTER :: QBAR(:)
  REAL(8)   ,POINTER :: PHIDOTPROJ(:)   ! <P|PHIDOT>
  REAL(8)   ,POINTER :: KTOPHI(:)       ! K -> |PHI>KTOPHI+|PHIBARDOT>KTOPHIDOT 
  REAL(8)   ,POINTER :: KTOPHIDOT(:)    ! K -> |PHI>KTOPHI+|PHIBARDOT>KTOPHIDOT 
  REAL(8)   ,POINTER :: JBARTOPHIDOT(:) ! JBAR ->  |PHIBARDOT> JBARTOPHIDOT ====
  INTEGER(4),POINTER :: LNSCATT(:)      ! LN OF CORRESPONDING SCATTERING CHANNEL
!!$  REAL(8)   ,POINTER :: DOVERLAPKK(:,:) !(LNX,LNX)
!!$  REAL(8)   ,POINTER :: DOVERLAPKJ(:,:) !(LNX,LNX)
!!$  REAL(8)   ,POINTER :: DOVERLAPJJ(:,:) !(LNX,LNX)
  LOGICAL(4)         :: TALLORB=.FALSE.  ! LIKE TORB(:)=TRUE, TEMPORARY SWITCH
  LOGICAL(4),POINTER :: TORB(:)         !(LNX)
  TYPE(TAILED_TYPE)            :: TAILED
!!$  TYPE(POTPARRED_TYPE)         :: SMALL
!!$  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSKPRIME
!!$  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSTAILEDK
!!$  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSTAILEDJBAR
!!$  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSKAUGMENT
!!$  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSJAUGMENT
END TYPE POTPAR_TYPE
!
TYPE PERIODICMAT_TYPE
  INTEGER(4)      :: IAT1 ! FIRST ATOM (LINKED TO THE RIGHT INDEX OF MAT)
  INTEGER(4)      :: IAT2 ! SECOND ATOM (LINKED TO THE LEFT INDEX OF MAT)
  INTEGER(4)      :: IT(3) ! LATTICE TRANSLATIONS TO BE ADDED TO ATOM 2
  INTEGER(4)      :: N1    ! RIGHT DIMENSION OF MAT
  INTEGER(4)      :: N2    ! LEFT DIMENSION OF MAT
  REAL(8),POINTER :: MAT(:,:)  !(N2,N1)
END TYPE PERIODICMAT_TYPE
!
TYPE PERIODICMAT2_TYPE
  INTEGER(4)      :: IAT1 ! FIRST ATOM (LINKED TO THE RIGHT INDEX OF MAT)
  INTEGER(4)      :: IAT2 ! SECOND ATOM (LINKED TO THE LEFT INDEX OF MAT)
  INTEGER(4)      :: IT(3) ! LATTICE TRANSLATIONS TO BE ADDED TO ATOM 2
  INTEGER(4)      :: N1    ! RIGHT DIMENSION OF MAT
  INTEGER(4)      :: N2    ! LEFT DIMENSION OF MAT
  INTEGER(4)      :: N3    ! #(MATRICES STORED)
  REAL(8),POINTER :: MAT(:,:,:)  !(N1,N2,N3)
END TYPE PERIODICMAT2_TYPE
!
TYPE UMAT_TYPE
  INTEGER(4)      :: NN1   ! FIRST ATOM PAIR REFERRING TO SBAR
  INTEGER(4)      :: NN2   ! SECOND ATOM PAIR REFERRING TO SBAR
  INTEGER(4)      :: IT(3) ! LATTICE TRANSLATIONS TO BE ADDED TO SECOND BOND
  INTEGER(4)      :: NA    ! ->IAT1(NN1)
  INTEGER(4)      :: NB    ! ->IAT1(NN2)
  INTEGER(4)      :: NC    ! ->IAT2(NN2)
  INTEGER(4)      :: ND    ! ->IAT2(NN1)
  REAL(8),POINTER :: UABCD(:,:,:,:)  !(NA,NB,NC,ND)
END TYPE UMAT_TYPE
!
TYPE UTENSOR_TYPE
  INTEGER(4)      :: IAT1  ! FIRST ATOM 
  INTEGER(4)      :: IAT2  ! FIRST ATOM 
  INTEGER(4)      :: IAT3  ! FIRST ATOM 
  INTEGER(4)      :: IAT4  ! FIRST ATOM 
  INTEGER(4)      :: IT2(3) ! LATTICE TRANSLATIONS OF 2. ATOM
  INTEGER(4)      :: IT3(3) ! LATTICE TRANSLATIONS OF 2. ATOM
  INTEGER(4)      :: IT4(3) ! LATTICE TRANSLATIONS OF 2. ATOM
  INTEGER(4)      :: N1    ! ->IAT1
  INTEGER(4)      :: N2    ! ->IAT2
  INTEGER(4)      :: N3    ! ->IAT3
  INTEGER(4)      :: N4    ! ->IAT4
  REAL(8),POINTER :: U(:,:,:,:)  !(N1,N2,N3,N4)
END TYPE UTENSOR_TYPE

TYPE OFFSITEX_TYPE
 INTEGER(4)         :: NDIS
 INTEGER(4)         :: NF
 REAL(8)   ,POINTER :: OVERLAP(:,:)  ! OVERLAP MATRIX ELEMENTS
 REAL(8)   ,POINTER :: X22(:,:)      !
 REAL(8)   ,POINTER :: X31(:,:)
 REAL(8)   ,POINTER :: BONDU(:,:)
 REAL(8)   ,POINTER :: DIS(:)
 REAL(8)   ,POINTER :: LAMBDA(:)
END TYPE OFFSITEX_TYPE
!===============================================================================
!== PARAMETER SECTION                                                         ==
!===============================================================================
LOGICAL(4)            :: TON=.FALSE.       
LOGICAL(4)            :: TOFFSITE=.TRUE.  !INCLUDE OFFSITE EXCHANGE
LOGICAL(4)            :: TDROP=.FALSE. ! WRITE THE WAVE FUNCTIONS TO FILE
LOGICAL(4)            :: TPICK=.FALSE. ! REAL HAMILTON CORRECTION FROM FILE

REAL(8)               :: K2=-0.25D0    ! 0.5*K2 IS THE KINETIC ENERGY
REAL(8)               :: RCSCALE=1.2D0  !RADIUS SCALE FACTOR FOR NEIGHBORLIST
!         RCSCALE TIMES THE SUM OF COVALENT RADII DEFINES CUTOFF FOR NEIGBORLIST
REAL(8)               :: HFWEIGHT=0.25D0
!
!== PARAMETER SET TO DEFINE THE GAUSSIANS FOR THE UNSCREENED HANKEL FUNCTIONS ==
!== IS USED IN LMTO_GAUSSFITKPRIME  ============================================
!!$INTEGER(4),PARAMETER  :: GAUSSFITKPRIME_NPOW=2   ! -1 INDICATES LX
!!$INTEGER(4),PARAMETER  :: GAUSSFITKPRIME_NE=12
!!$REAL(8)   ,PARAMETER  :: GAUSSFITKPRIME_R1=0.6667D0
!!$REAL(8)   ,PARAMETER  :: GAUSSFITKPRIME_SCALER=1.25D0
!!$!== PARAMETER SET TO DEFINE THE GAUSSIANS FOR THE SCREENED HANKEL FUNCTIONS   ==
!!$INTEGER(4),PARAMETER  :: GAUSSFITKBARPRIME_NPOW=4
!!$INTEGER(4),PARAMETER  :: GAUSSFITKBARPRIME_NE=4
!!$REAL(8)   ,PARAMETER  :: GAUSSFITKBARPRIME_R1=1.D0
!!$REAL(8)   ,PARAMETER  :: GAUSSFITKBARPRIME_SCALER=1.5D0
!!$!== PARAMETER SET TO DEFINE THE GAUSSIANS FOR THE AUGMENTATION
!!$INTEGER(4),PARAMETER  :: GAUSSFITKAUGMENT_NPOW=4
!!$INTEGER(4),PARAMETER  :: GAUSSFITKAUGMENT_NE=6
!!$REAL(8)   ,PARAMETER  :: GAUSSFITKAUGMENT_R1=3.D-2
!!$REAL(8)   ,PARAMETER  :: GAUSSFITKAUGMENT_SCALER=2.D0
!
!===============================================================================
!== VARIABLE SECTION                                                          ==
!===============================================================================
LOGICAL(4)              :: TINI=.FALSE.
LOGICAL(4)              :: TINISTRUC=.FALSE.
LOGICAL(4)              :: THTBC=.FALSE. ! HTBC CALCULATED
INTEGER(4)              :: NSP=-1
INTEGER(4)              :: ISPSELECTOR=-1 ! USED ONLY FOR HYBRIDSETTING
TYPE(HYBRIDSETTING_TYPE),ALLOCATABLE :: HYBRIDSETTING(:)
INTEGER(4)              :: NRL   ! DIMENSION OF STRUCTURE CONSTANTS IN K-SPACE
INTEGER(4),ALLOCATABLE  :: LNX(:)      !(NSP)
INTEGER(4),ALLOCATABLE  :: LOX(:,:)    !(LNXX,NSP)
INTEGER(4),ALLOCATABLE  :: ISPECIES(:) !(NAT)
REAL(8)   ,ALLOCATABLE  :: ORBRAD(:,:) !(LXX+1,NAT) NODE-POSITION OF THE ORBITAL
TYPE(POTPAR_TYPE)     ,ALLOCATABLE :: POTPAR(:) !POTENTIAL PARAMETERS
INTEGER(4)            ,ALLOCATABLE :: SBARATOMI1(:)
INTEGER(4)            ,ALLOCATABLE :: SBARATOMI2(:)
INTEGER(4)            ,ALLOCATABLE :: SBARLI1(:,:)
!== GAUSSIAN PART OF NTBOS FROM SUPERPOSITION OF HANKEL FUNCTIONS ==============
TYPE(ORBITALGAUSSCOEFF_TYPE),ALLOCATABLE :: GAUSSORB(:) !(NAT)
!== GAUSSIAN PART OF NTBOS FROM TAILED HANKEL AND BESSEL FUNCTIONS =============
TYPE(ORBITALGAUSSCOEFF_TYPE),ALLOCATABLE :: GAUSSORB_T(:) !(NAT)
!== AUGMENTED NTBOS IN TERMS OF GAUSSIANS ======================================
TYPE(ORBITALGAUSSCOEFF_TYPE),ALLOCATABLE :: GAUSSORBAUG(:) !(NAT)
TYPE(ORBITALSPHHARM_TYPE)   ,ALLOCATABLE :: LMORB(:)
TYPE(UTENSOR_TYPE)          ,ALLOCATABLE :: UTENSOR(:)
TYPE(OFFSITEX_TYPE)         ,ALLOCATABLE :: OFFSITEX(:,:)
LOGICAL(4)                  ,PARAMETER :: TSPHERICAL=.FALSE.
!===============================================================================
!=====  STRUCTURE DEPENDENT DATA  ==============================================
!===============================================================================
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: SBAR(:)   !(NNS) SCREENED STRUCTURE CONST.
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DENMAT(:) !(NND) DENSITY MATRIX
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: HAMIL(:)  !(NND) DERIVATIVE OF ENERGY
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: OVERLAP(:)!(NNS) OVERLAP MATRIX ONLY MAIN
!INTEGER(4)                         :: NNU       !#(ELEMENTS OF UMAT)
!TYPE(UMAT_TYPE)       ,ALLOCATABLE :: UMAT(:)   !(NNUX/NNU) U-MATRIX ELEMENTS
!!$TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DENMAT_T(:) !(NNS) DENSITY MATRIX
!!$TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: HAMIL_T(:)  !(NNS) DERIVATIVE OF ENERGY
!!$INTEGER(4)                         :: NNUX      !DIMENSION OF UMAT
END MODULE LMTO_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SETL4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'ON') THEN
        TON=VAL
      ELSE IF(ID.EQ.'OFFSITE') THEN
        TOFFSITE=VAL
      ELSE IF(ID.EQ.'DROP') THEN
        TDROP=VAL
      ELSE IF(ID.EQ.'PICK') THEN
        TPICK=VAL
      ELSE IF(ID.EQ.'DHOFK') THEN
        CALL LMTO_DROPPICK$SETL4('DHOFK',VAL)
!
!     ==========================================================================
!     == IF ACTIVE THE HYBRID CONTRIBTIONS ON THIS ATOM ARE CONSIDERED        ==
!     ==========================================================================
      ELSE IF(ID.EQ.'ACTIVE') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%ACTIVE=VAL
!
      ELSE IF(ID.EQ.'COREVALENCE') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TCV=VAL
!
      ELSE IF(ID.EQ.'FOCKSETUP') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TFOCKSETUP=VAL
!
      ELSE IF(ID.EQ.'NDDO') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TNDDO=VAL
!
      ELSE IF(ID.EQ.'31') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%T31=VAL
!
      ELSE IF(ID.EQ.'BONDX') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TBONDX=VAL
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$GETL4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'ON') THEN
        VAL=TON
!
      ELSE IF(ID.EQ.'THTBC') THEN
        VAL=THTBC
!
      ELSE IF(ID.EQ.'FOCKSETUP') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$GETL4')
        END IF
        VAL=HYBRIDSETTING(ISPSELECTOR)%TFOCKSETUP
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO$GETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SETI4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'ISP') THEN
!       ========================================================================
!       == THE ATOM TYPE SELECTOR ISP IS USED ONLY FOR SETTINGS               ==
!       ========================================================================
        ISPSELECTOR=VAL
        CALL ATOMTYPELIST$LENGTH(NSP)
        IF(NSP.LE.0) THEN
          CALL ERROR$MSG('NUMBER OF ATOM TYPES UNKNOWN BY ATOMTYPELIST')
          CALL ERROR$I4VAL('NSP',NSP)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETI4')
        END IF
        IF(VAL.GT.NSP) THEN
          CALL ERROR$MSG('ATOM TYPE SPECIFIER OUT OF RANGE')
          CALL ERROR$I4VAL('ISP',ISPSELECTOR)
          CALL ERROR$I4VAL('NSP',NSP)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETI4')
        END IF
        IF(.NOT.ALLOCATED(HYBRIDSETTING)) THEN
!         == ALLOCATE AND SET DEFAULT VALUES ===================================
          ALLOCATE(HYBRIDSETTING(NSP))
          HYBRIDSETTING(:)%ACTIVE=.FALSE.
          HYBRIDSETTING(:)%TCV=.TRUE.
          HYBRIDSETTING(:)%TFOCKSETUP=.TRUE.
          HYBRIDSETTING(:)%LHFWEIGHT=-1.D0
          HYBRIDSETTING(:)%TNDDO=.TRUE.
          HYBRIDSETTING(:)%T31=.TRUE.
          HYBRIDSETTING(:)%TBONDX=.TRUE.
          HYBRIDSETTING(:)%TAILEDLAMBDA1=4.D0
          HYBRIDSETTING(:)%TAILEDLAMBDA2=2.D0
!          HYBRIDSETTING(:)%HFWEIGHT=0.25D0
!          HYBRIDSETTING(:)%K2=-0.25D0
!          HYBRIDSETTING(:)%RANGESCALE=1.2D0
        END IF
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SETR8(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
      INTEGER(4)              :: I
!     **************************************************************************
      IF(ID.EQ.'TAILLAMBDA1') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TAILEDLAMBDA1=VAL
!
!
      ELSE IF(ID.EQ.'LHFWEIGHT') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%LHFWEIGHT=VAL
!
      ELSE IF(ID.EQ.'TAILLAMBDA2') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TAILEDLAMBDA2=VAL
!
      ELSE IF(ID.EQ.'SCALERCUT') THEN
        RCSCALE=VAL
!
      ELSE IF(ID.EQ.'HFWEIGHT') THEN
        HFWEIGHT=VAL
!
      ELSE IF(ID.EQ.'K2') THEN
        K2=VAL

      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO$SETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$GETR8(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : HFWEIGHT,NSP,HYBRIDSETTING,ISPSELECTOR
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)    ,INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'HFWEIGHT') THEN
        VAL=HFWEIGHT
!
      ELSE IF(ID.EQ.'LHFWEIGHT') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$GETR8')
        END IF
        VAL=HYBRIDSETTING(ISPSELECTOR)%LHFWEIGHT
        IF(VAL.LT.0.D0) VAL=HFWEIGHT

      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO$GETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$REPORT(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : TON,NSP,HYBRIDSETTING,RCSCALE,HFWEIGHT,K2 &
     &                       ,TOFFSITE,POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: THISTASK,NTASKS
      CHARACTER(32)         :: ID
      INTEGER(4)            :: ISP,I
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: LNX1
      INTEGER(4),ALLOCATABLE :: LOX1(:)
!     **************************************************************************
      IF(.NOT.TON) RETURN
      CALL LMTO_INITIALIZE()
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      CALL REPORT$TITLE(NFIL,'LMTO OBJECT:  GENERIC VARIABLES')
      CALL REPORT$R8VAL(NFIL,'EXCHANGE ADMIXTURE ',HFWEIGHT,'')
      CALL REPORT$R8VAL(NFIL,'RANGESCALE ',RCSCALE,'*(RCOV(A)+RCOV(B))')
      CALL REPORT$R8VAL(NFIL,'K2 ',K2,'A.U.')
!
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETCH('ID',ID)
        CALL SETUP$GETI4('LNX',LNX1)
        ALLOCATE(LOX1(LNX1))
        CALL SETUP$GETI4A('LOX',LNX1,LOX1)
        WRITE(NFIL,*)
        CALL REPORT$TITLE(NFIL,'LMTO OBJECT: '//TRIM(ID))
!!$        CALL REPORT$R8VAL(NFIL,'EXCHANGE ADMIXTURE ' &
!!$     &                         ,HYBRIDSETTING(ISP)%HFWEIGHT,'')
        CALL REPORT$STRING(NFIL,'ACTIVE ORBITALS (L=SWITCH):')
        WRITE(NFIL,FMT='(100("   ",I1,"=",L1))') &
       &                                  (LOX1(I),POTPAR(ISP)%TORB(I),I=1,LNX1)
        IF(HYBRIDSETTING(ISP)%LHFWEIGHT.GE.0.D0) THEN
          CALL REPORT$R8VAL(NFIL,'LOCAL EXCHANGE ADMIXTURE ' &
     &                        ,HYBRIDSETTING(ISP)%LHFWEIGHT,'')
        END IF
        CALL REPORT$L4VAL(NFIL,'SETUP WITH FOCK TERM ' &
     &                        ,HYBRIDSETTING(ISP)%TFOCKSETUP)
        CALL REPORT$L4VAL(NFIL,'CORE VALENCE EXCHANGE ' &
     &                        ,HYBRIDSETTING(ISP)%TCV)
        TCHK=TOFFSITE.AND.HYBRIDSETTING(ISP)%TNDDO
        CALL REPORT$L4VAL(NFIL,'NDDO OFFSITE TERMS ',TCHK)
        TCHK=TOFFSITE.AND.HYBRIDSETTING(ISP)%T31
        CALL REPORT$L4VAL(NFIL,'31 OFFSITE TERMS ',TCHK)
        TCHK=TOFFSITE.AND.HYBRIDSETTING(ISP)%TBONDX
        CALL REPORT$L4VAL(NFIL,'BONDX OFFSITE TERMS ',TCHK)
!!$        CALL REPORT$R8VAL(NFIL,'K2 ' &
!!$     &                        ,HYBRIDSETTING(ISP)%K2,' ')
!!$        CALL REPORT$R8VAL(NFIL,'RANGESCALE ' &
!!$     &                        ,HYBRIDSETTING(ISP)%RANGESCALE,' ')
        CALL REPORT$R8VAL(NFIL,'LAMDA1 FOR NTBO TAILS ' &
     &                        ,HYBRIDSETTING(ISP)%TAILEDLAMBDA1,'1/ABOHR ')
        CALL REPORT$R8VAL(NFIL,'LAMDA2 FOR NTBO TAILS ' &
     &                        ,HYBRIDSETTING(ISP)%TAILEDLAMBDA2,'1/ABOHR ')
        DEALLOCATE(LOX1)
        CALL SETUP$UNSELECT()
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$MAKESTRUCTURECONSTANTS()
!     **************************************************************************
!     **  PRODUCES THE SCREENED STRUCTURE CONSTANTS SBAR                      **
!     **                                                                      **
!     **  IS NEEDED ALSO BY LDAPLUSU VIA LMTO$DOLOCORB                        **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : K2,RCSCALE,SBAR,TINISTRUC,POTPAR &
     &                      ,ISPECIES,NSP,LOX,LNX,SBARLI1
      USE PERIODICTABLE_MODULE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER   :: NNXPERATOM=100
      INTEGER(4)             :: NAT       !#(ATOMS)
      REAL(8)                :: RBAS(3,3) !LATTICE VECTORS
      REAL(8)   ,ALLOCATABLE :: R0(:,:)   !(3,NAT) ATOMIC POSITIONS
      INTEGER(4)             :: NNX
      INTEGER(4),ALLOCATABLE :: NNLIST(:,:) !(5,NNX)
      INTEGER(4)             :: NAT1
      INTEGER(4)             :: NNS
      INTEGER(4)             :: NORB
      INTEGER(4)             :: N
      INTEGER(4),ALLOCATABLE :: LX1(:)     !(NNS) MAX(ANGULAR MOMENTUM)
      REAL(8)   ,ALLOCATABLE :: RPOS(:,:)  !(3,NNS(IAT)) ATOMIC POSITIONS
      REAL(8)   ,ALLOCATABLE :: QBAR1(:,:) !(LXX+1,NSP) SCREENING PARAMETER
      REAL(8)   ,ALLOCATABLE :: QBAR(:)    !(N) SCREENING PARAMETER
      REAL(8)   ,ALLOCATABLE :: SBAR1(:,:) !
      REAL(8)                :: SVAR
      INTEGER(4)             :: IAT,IAT1,IAT2,ISP,ISP1,ISP2,LMX1,LMX2 
      INTEGER(4)             :: L,NN,NN1,NN2,NN0,I,LN,I1,I2,LX,L1,L2
      INTEGER(4)             :: I11,I12,I21,I22
      INTEGER(4)             :: J11,J12,J21,J22
      LOGICAL(4)             :: TCHK
      REAL(8)   ,ALLOCATABLE :: RC(:)
      INTEGER(4)             :: LMX
      INTEGER(4)             :: NTASKS,THISTASK,FROMTASK
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
PRINT*,'MARKE 1'
CALL SETUP$ISELECT(1)
CALL SETUP$ISELECT(0)
PRINT*,'MARKE 2'
      CALL SETUP$GETL4('INTERNALSETUPS',TCHK)
PRINT*,'MARKE 3'
      IF(.NOT.TCHK) RETURN
                              CALL TRACE$PUSH('LMTO$MAKESTRUCTURECONSTANTS')
                              CALL TIMING$CLOCKON('LMTO STRUCTURECONSTANTS')
      TINISTRUC=.TRUE.
!
!     ==========================================================================
!     ==  INITIALIZE LMTO OBJECT                                              ==
!     ==========================================================================
      CALL LMTO_INITIALIZE()
!
!     ==========================================================================
!     == OBTAIN ATOMIC STRUCTURE                                              ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     == DETERMINE SCREENING PARAMETER QBAR                                   ==
!     ==========================================================================
      LX=MAXVAL(LOX)
      ALLOCATE(QBAR1((LX+1)**2,NSP))
      QBAR1(:,:)=0.D0
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I1=L**2+1
          I2=(L+1)**2
          QBAR1(I1:I2,ISP)=POTPAR(ISP)%QBAR(POTPAR(ISP)%LNSCATT(LN))
        ENDDO
        CALL SETUP$UNSELECT()
      ENDDO
!
!     ==========================================================================
!     == NEIGHBORLIST   NNLIST(:,NN)=(IAT1,IAT2,IT1,IT2,IT3)                  ==
!     ==                     IT1,IT2,IT3 ARE THE LATTICE TRANSLATIONS OF IAT2 ==
!     ==========================================================================
      ALLOCATE(RC(NAT))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',SVAR)
        CALL SETUP$UNSELECT()
        CALL PERIODICTABLE$GET(SVAR,'R(COV)',RC(IAT))
      ENDDO
      RC(:)=RC(:)*RCSCALE
      NNX=NNXPERATOM*NAT
      ALLOCATE(NNLIST(5,NNX))
      CALL LMTO$NEIGHBORLIST(RBAS,NAT,R0,RC,NNX,NNS,NNLIST)
      DEALLOCATE(RC)
!
!     ==========================================================================
!     == CLEAN STUCTURE DEPENDENT ARRAYS                                      ==
!     ==========================================================================
      CALL LMTO_STRUCRESET()
!
!     ==========================================================================
!     == STRUCTURE CONSTANTS                                                  ==
!     ==========================================================================
      ALLOCATE(SBAR(NNS))
      DO IAT1=1,NAT
!       == MEMBERS NN1:NN2 ON THE NEIGHBOLIST BUILD THE CLUSTER AROUND ATOM 1 ==
!       == MEMBER NN0 IS THE ONSITE MEMBER                                    ==
        NN1=1
        NN0=0
        DO NN=1,NNS
          IF(NNLIST(1,NN).GT.IAT1)EXIT
          NN2=NN   ! NN2 IS THE LAST MEMBER WITH IAT1 AS FIRST ATOM 
          IF(NNLIST(1,NN).LT.IAT1)NN1=NN+1
          IF(NNLIST(1,NN).EQ.IAT1) THEN
            IF(NNLIST(2,NN).EQ.IAT1) THEN
              IF(NNLIST(3,NN).EQ.0.AND.NNLIST(4,NN).EQ.0 &
     &                            .AND.NNLIST(5,NN).EQ.0) THEN
                NN0=NN   ! NN0 IS ONSITE TERM FRO ATOM IAT1
              END IF
            END IF
          END IF
        ENDDO
!
        NAT1=NN2-NN1+1  ! #(ATOMS ON THE CLUSTER )
        ALLOCATE(LX1(NAT1))
        ALLOCATE(RPOS(3,NAT1))
        N=0  ! #ORBITALS ON THE CLUSTER
        DO NN=NN1,NN2
          IAT2=NNLIST(2,NN)
          ISP=ISPECIES(IAT2)
          LX=MAXVAL(LOX(:LNX(ISP),ISP))
          N=N+(LX+1)**2 ! #ORBITALS ON THE CLUSTER
          LX1(NN-NN1+1)=LX ! X(ANGULAR MOMENTUM FOR THIS ATOM)
          RPOS(:,NN-NN1+1)=R0(:,IAT2)+RBAS(:,1)*REAL(NNLIST(3,NN),KIND=8) &
     &                               +RBAS(:,2)*REAL(NNLIST(4,NN),KIND=8) &
     &                               +RBAS(:,3)*REAL(NNLIST(5,NN),KIND=8)
        ENDDO
        NORB=(LX1(1)+1)**2  ! #(ORBITALS ON THE CENTRAL ATOM)
        IF(NN0.NE.NN1) THEN
          CALL ERROR$MSG('ONSITE ELEMENT NOT FIRST IN NEIGHBORLIST')
          CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED')
          CALL ERROR$STOP('LMTO$MAKESTRUCTURECONSTANTS')
        END IF
!
!       ========================================================================
!       == EXPAND SCREENING PARAMETER QBAR                                    ==
!       ========================================================================
        ALLOCATE(QBAR(N))
        QBAR(:)=0.D0
        I=0
        DO NN=NN1,NN2
          IAT2=NNLIST(2,NN)
          ISP=ISPECIES(IAT2)
          LX=MAXVAL(LOX(:LNX(ISP),ISP))
          LMX=(LX+1)**2
          QBAR(I+1:I+LMX)=QBAR1(1:LMX,ISP)
          I=I+LMX
        ENDDO
!
!       ========================================================================
!       == DETERMINE STRUCTURE CONSTANTS                                      ==
!       == HERE, THE STRUCTURE CONSTANTS USE A COMPLETE SET OF ANGULAR MOMENTA==
!       == UP TO A MAXIMUM ANGULAR MOMENTUM. NOT ALL WILL BE USED LATER ON    ==
!       ========================================================================
        ALLOCATE(SBAR1(NORB,N))
!NORB=(LX1(1)+1)**2
!N=(LX(1)+2)**2?
PRINT*,'DOING LMTO$CLUSTERSTRUCTURECONSTANTS.....'
       IF(MOD(IAT1-1,NTASKS).NE.THISTASK-1) THEN
          CALL TIMING$CLOCKON('STRUCTURE CONSTANTS')
          CALL LMTO$CLUSTERSTRUCTURECONSTANTS(K2,NAT1,RPOS,LX1,QBAR &
       &                                                      ,NORB,N,SBAR1)
          CALL TIMING$CLOCKOFF('STRUCTURE CONSTANTS')
       ELSE
         SBAR1=0.D0
       END IF
PRINT*,'..... LMTO$CLUSTERSTRUCTURECONSTANTS  DONE'
!
!       ========================================================================
!       == MAP ONTO SBAR                                                      ==
!       ========================================================================
        I=0
        DO NN=NN1,NN2
          IAT2=NNLIST(2,NN)
          SBAR(NN)%IAT1=NNLIST(1,NN)
          SBAR(NN)%IAT2=NNLIST(2,NN)
          SBAR(NN)%IT(:)=NNLIST(3:5,NN)
          ISP1=ISPECIES(IAT1)
          ISP2=ISPECIES(IAT2)
          LX=MAXVAL(LOX(:,ISP1))
          LMX1=SBARLI1(LX+1,ISP1)+2*LX
          LX=MAXVAL(LOX(:,ISP2))
          LMX2=SBARLI1(LX+1,ISP2)+2*LX
          SBAR(NN)%N1=LMX1
          SBAR(NN)%N2=LMX2
          ALLOCATE(SBAR(NN)%MAT(LMX1,LMX2))
          DO L2=0,LX1(NN-NN1+1)
            I21=I+L2**2+1
            I22=I+(L2+1)**2
            J21=SBARLI1(L2+1,ISP2)
            J22=J21+2*L2
            IF(J21.LT.0) CYCLE
            DO L1=0,LX1(1)
              I11=L1**2+1
              I12=(L1+1)**2
              J11=SBARLI1(L1+1,ISP1)
              J12=J11+2*L1
              IF(J11.LT.0) CYCLE
              SBAR(NN)%MAT(J11:J12,J21:J22)=SBAR1(I11:I12,I21:I22)    !C
            ENDDO
          ENDDO
          I=I+(LX+1)**2
        ENDDO
        DEALLOCATE(LX1)
        DEALLOCATE(RPOS)
        DEALLOCATE(QBAR)
        DEALLOCATE(SBAR1)
      ENDDO
      DEALLOCATE(QBAR1)
!
!     == DISTRIBUTE STRUCTURE CONSTANTS
      DO NN=1,NNS
        IAT=SBAR(NN)%IAT1
        FROMTASK=1+MOD(IAT-1,NTASKS)
        CALL MPE$BROADCAST('MONOMER',FROMTASK,SBAR(NN)%MAT)
      ENDDO
                              CALL TIMING$CLOCKOFF('LMTO STRUCTURECONSTANTS')
                              CALL TRACE$POP()
      RETURN
      END
!     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_STRUCRESET()
!     **************************************************************************
!     **  DEALLOCATE ALL STRUCTURE DEPENDENT ARRAYS IN LMTO_MODULE            **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : SBAR,DENMAT,HAMIL,OVERLAP
      INTEGER(4)     :: NN
      INTEGER(4)     :: I
!     **************************************************************************
!
!     == CLEAN STRUCTURE CONSTANTS =============================================
      IF(ALLOCATED(SBAR)) THEN
        NN=SIZE(SBAR)
        DO I=1,NN
          DEALLOCATE(SBAR(I)%MAT)
        ENDDO
        DEALLOCATE(SBAR)
      END IF
!
!     == CLEAN DENSITY MATRIX ==================================================
      IF(ALLOCATED(DENMAT)) THEN
        NN=SIZE(DENMAT)
        DO I=1,NN
          DEALLOCATE(DENMAT(I)%MAT)
        ENDDO
        DEALLOCATE(DENMAT)
      END IF
!
!     == CLEAN HAMILTON MATRIX =================================================
      IF(ALLOCATED(HAMIL)) THEN
        NN=SIZE(HAMIL)
        DO I=1,NN
          DEALLOCATE(HAMIL(I)%MAT)
        ENDDO
        DEALLOCATE(HAMIL)
      END IF
!
!     == CLEAN OVERLAP MATRIX ==================================================
      IF(ALLOCATED(OVERLAP)) THEN
        NN=SIZE(OVERLAP)
        DO I=1,NN
          DEALLOCATE(OVERLAP(I)%MAT)
        ENDDO
        DEALLOCATE(OVERLAP)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_INITIALIZE()
!     **************************************************************************
!     **  PREPARES POTENTIAL PARAMETERS AND SIMILAR BASIC DATA.               **
!     **  IT IS CALLED BY LMTO$MAKESTRUCTURECONSTANTS                         **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : TINI,TON,HYBRIDSETTING,POTPAR,NSP,TOFFSITE
      IMPLICIT NONE
      INTEGER(4) :: NAT,ISP
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
!
!     ==========================================================================
!     == COLLECT NSP,LNX,LOX,ISPECIES                                         ==
!     ==========================================================================
      CALL LMTO_COLLECTMAPARRAYS()
!
!     ==========================================================================
!     == COLLECT SBARLI1(L+1,NSP)                                             ==
!     ==========================================================================
      CALL LMTO_SBARINDICES()   
!
!     ==========================================================================
!     == DETERMINE POTENTIAL PARAMETERS                                       ==
!     == RAD,LNSCATT,PHIDOTPROJ,QBAR,KTOPHI,KTOPHIDOT,JBARTOPHIDOT            ==
!     ==========================================================================
      CALL LMTO_MAKEPOTPAR()
!
      DO ISP=1,NSP
        IF(.NOT.HYBRIDSETTING(ISP)%ACTIVE)POTPAR(ISP)%TORB=.FALSE.
      ENDDO
!
!     ==========================================================================
!     == ATTACH EXPONENTIAL TAILS TO AUGMENTED HANKEL AND BESSEL FUNCTIONS    ==
!     ==========================================================================
      CALL LMTO_MAKETAILEDPARTIALWAVES()
      IF(.NOT.TON) RETURN
!
!     ==========================================================================
!     ==  CONSTRUCT GAUSSIAN FITS OF THE TAILES ORBITALS                      ==
!     ==========================================================================
      CALL LMTO_TAILEDGAUSSFIT()
!
!     ==========================================================================
!     ==  CONSTRUCT OFFSITE INTEGRALS OF TAILED ORBITALS                      ==
!     ==========================================================================
      IF(TOFFSITE) THEN
                            CALL TIMING$CLOCKON('OFFSITE U-TENSOR')
        CALL LMTO_OFFXINT()
                            CALL TIMING$CLOCKOFF('OFFSITE U-TENSOR')
      END IF

!      CALL LMTO_TAILEDPRODUCTS()
!!$!
!!$!     ==========================================================================
!!$!     == DETERMINE GAUSS EXPANSION OF UNSCREENED HANKEL FUNCTIONS KPRIME      ==
!!$!     ==========================================================================
!!$      CALL LMTO_GAUSSFITKPRIME()
!!$      CALL LMTO_GAUSSFITKAUGMENT()
!!$!
!!$!     ==========================================================================
!!$!     == DETERMINE KPRIME AND JBAR WITH EXPONENTIAL TAILS                     ==
!!$!     == USED TO BUILD UP APPROXIMATE NTBOS IN A ONE-CENTER EXPANSION        ==
!!$!     ==========================================================================
!!$      CALL LMTO_GAUSSFITKJTAILS()
!!$      CALL LMTO_ONSITEOVERLAP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_COLLECTMAPARRAYS()
!     **************************************************************************
!     **  STORES A LOCAL COPY OF                                              **
!     **     NSP                                                              **
!     **     LNX(NSP)                                                         **
!     **     LOX(LNXX,NSP)                                                    **
!     **     ISPECIES(NAT)                                                    **
!     **  IN THE LMTO_MODULE                                                  **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY :NSP,LNX,LOX,ISPECIES,ISPSELECTOR
      IMPLICIT NONE
      INTEGER(4)    :: ISP,NAT
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_COLLECTMAPARRAYS')
      CALL SETUP$GETI4('NSP',NSP) !FORMER CALL SETUP$NSPECIES(NSP)

      ALLOCATE(LNX(NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX(ISP))
        CALL SETUP$UNSELECT()
      ENDDO
!
      ALLOCATE(LOX(MAXVAL(LNX),NSP))
      LOX(:,:)=-1
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4A('LOX',LNX(ISP),LOX(1:LNX(ISP),ISP))
        CALL SETUP$UNSELECT()
      ENDDO

      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES) 
!
!     ==========================================================================
!     == THE FOLLOWING ENSURES THAT HYBRIDSETTING IS ALLOCATED =================
!     ==========================================================================
      IF(ISPSELECTOR.LE.0) THEN
        CALL LMTO$SETI4('ISP',1)
        CALL LMTO$SETI4('ISP',0)
      END IF
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SBARINDICES()
!     **************************************************************************
!     **  CREATES THE INDEX ARRAYS TO BE USED WITH THE SCREENED STRUCTURE     **
!     **  CONSTANTS.                                                          **
!     **     SBARATOMI1(IAT) POINTS TO THE FIRST ENTRY FOR ATOM IAT IN THE    **
!     **                     K-SPACE MATRIX                                   **
!     **     SBARATOMI2(IAT) POINTS TO THE LAST ENTRY FOR ATOM IAT IN THE     **
!     **                     K-SPACE MATRIX                                   **
!     **     SBARLI1(L+1,ISP) POINTS TO THE FIRST OF 2*L+1 ENTRIES FOR ATOM   **
!     **                     IAT AND MAIN ANGULAR MOMENTUM L                  **
!     **                     (RELATIVE TO THE FIRST ELEMENT FOR THIS ATOM)    **
!     **                                                                      **
!     **  THE ARRAYS ISPECIES1 AND LX1 HAVE STRANGE NAMES BECAUSE THE         **
!     **  ORIGINAL                                                            **
!     **  NAMES ARE ALREADY USED BY LMTO_MODULE. WE USE SEPARATE ARRAYS,      **
!     **  BECAUSE WE WANT TO REMOVE THESE ARRAYS FROM THE MODULE              **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE,ONLY : NSP,ISPECIES,LNX,LOX &
     &                      ,NRL,SBARATOMI1,SBARATOMI2,SBARLI1
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      INTEGER(4)             :: LX1
      INTEGER(4)             :: IPOS,IAT,ISP,L,LN
!     **************************************************************************
      IF(ALLOCATED(SBARATOMI1)) RETURN
                              CALL TRACE$PUSH('LMTO$SBARINDICES')
!
!     ==========================================================================
!     == POINTS TO THE FIRST STRUCTURE CONSTANT ELEMENT FOR A GIVEN ATOM      ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(SBARATOMI1(NAT))
      ALLOCATE(SBARATOMI2(NAT))
      IPOS=1
      DO IAT=1,NAT
        ISP=ISPECIES(IAT) 
        SBARATOMI1(IAT)=IPOS
        LX1=MAXVAL(LOX(:LNX(ISP),ISP))
        DO L=0,LX1
          DO LN=1,LNX(ISP)
            IF(LOX(LN,ISP).EQ.L) THEN
              IPOS=IPOS+2*L+1
              EXIT 
            END IF
          ENDDO
        ENDDO
        SBARATOMI2(IAT)=IPOS-1
      ENDDO
!
!     ==========================================================================
!     == NUMBER OF ANGULAR MOMENTA IN THE STRUCTURE CONSTANTS                 ==
!     ==========================================================================
      NRL=SBARATOMI2(NAT)  !#(TIGHT-BINDING ORBITALS)
!
!     ==========================================================================
!     == POINTER TO THE FIRST STRUCTURE CONSTANT ELEMENT OF ANGULAR MOMENTUM L==
!     ==========================================================================
      ALLOCATE(SBARLI1(MAXVAL(LOX)+1,NSP))
      SBARLI1(:,:)=-1
      DO ISP=1,NSP
        IPOS=1
        LX1=MAXVAL(LOX(:LNX(ISP),ISP))
        DO L=0,LX1
          DO LN=1,LNX(ISP)
            IF(LOX(LN,ISP).EQ.L) THEN
              SBARLI1(L+1,ISP)=IPOS
              IPOS=IPOS+2*L+1
              EXIT 
            END IF
          ENDDO
        ENDDO
      ENDDO
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_MAKEPOTPAR()
!     **************************************************************************
!     **  POTPAR(ISP)%RAD              : ASA RADIUS                           **
!     **  POTPAR(ISP)%LNSCATT(LN)      : SCATTERING WAVE FUNCTION FROM LN     **
!     **  POTPAR(ISP)%PHIDOTPROJ(LN)   : <PS-PRO|PS-PHIBARDOT>                **
!     **  POTPAR(ISP)%QBAR(LN)         : Q-BAR                                **
!     **  POTPAR(ISP)%KTOPHI(LN)       : K -> |PHI>KTOPHI                     **
!     **  POTPAR(ISP)%KTOPHIDOT(LN)    :    + |PHIBARDOT>KTOPHIDOT            **
!     **  POTPAR(ISP)%JBARTOPHIDOT(LN) : JBAR -> |PHIBARDOT>JBARTOPHIDOT      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : K2,POTPAR,NSP,LNX,LOX
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),ALLOCATABLE :: ISCATT(:)
      INTEGER(4)             :: LNX1
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: NLPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: PRO(:,:)
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      REAL(8)                :: AEZ
      REAL(8)                :: RAD
      REAL(8)                :: PHIVAL,PHIDER
      REAL(8)                :: PHIDOTVAL,PHIDOTDER
      REAL(8)                :: KVAL,KDER
      REAL(8)                :: JVAL,JDER
      REAL(8)                :: WJPHI,WJPHIDOT,WKPHI,WKPHIDOT,WJBARPHI
      REAL(8)                :: WPHIPHIDOT
      REAL(8)                :: QBAR
      REAL(8)                :: SVAR
      INTEGER(4)             :: ISP,LN,L,LN1,LN2
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO_MAKEPOTPAR')
      ALLOCATE(POTPAR(NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        LNX1=LNX(ISP)
!
!       ========================================================================
!       ==  COLLECT DATA                                                      ==
!       ========================================================================
!       == RADIAL GRID =========================================================
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
!       == MATCHING RADIUS =====================================================
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(NINT(AEZ),'R(ASA)',RAD) 
        POTPAR(ISP)%RAD=RAD
!       == SELECTION OF LOCAL ORBITALS CONSIDERED IN THE U-TENSOR ==============
        ALLOCATE(POTPAR(ISP)%TORB(LNX1))
        CALL SETUP$GETL4A('TORB',LNX1,POTPAR(ISP)%TORB)
!       == PARTIAL WAVES AND PROJECTORS ========================================
        ALLOCATE(NLPHI(NR,LNX1))
        ALLOCATE(AEPHI(NR,LNX1))
        ALLOCATE(NLPHIDOT(NR,LNX1))
        ALLOCATE(AEPHIDOT(NR,LNX1))
        ALLOCATE(PSPHI(NR,LNX1))
        ALLOCATE(PSPHIDOT(NR,LNX1))
        ALLOCATE(PRO(NR,LNX1))
        CALL SETUP$GETR8A('AEPHI',NR*LNX1,AEPHI)
        CALL SETUP$GETR8A('PSPHI',NR*LNX1,PSPHI)
        CALL SETUP$GETR8A('QPHI',NR*LNX1,NLPHI)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX1,AEPHIDOT)
        CALL SETUP$GETR8A('PSPHIDOT',NR*LNX1,PSPHIDOT)
        CALL SETUP$GETR8A('QPHIDOT',NR*LNX1,NLPHIDOT)
        CALL SETUP$GETR8A('PRO',NR*LNX1,PRO)
!
!       ========================================================================
!       ==  SELECT ONE PHIBARDOT FUNCTION PER L                               ==
!       ========================================================================
!       == GET INFO ON SCATTERING CHANNELS =====================================
!       == ISCATT=0 FOR VALENCE, ISCATT>0 FOR SCATTERING, ISCATT<0 FORSEMICORE =
        ALLOCATE(ISCATT(LNX1))
        CALL SETUP$GETI4A('ISCATT',LNX1,ISCATT)
!       == DETERIMINE SELECTION MAP ============================================
        ALLOCATE(POTPAR(ISP)%LNSCATT(LNX1))
        POTPAR(ISP)%LNSCATT(:)=-1
        DO LN=1,LNX1
          IF(POTPAR(ISP)%LNSCATT(LN).NE.-1) CYCLE  ! ALREADY DONE
!         == SELECT CHANNEL FOR SCATTERING WAVE FUNCTION =====================
          POTPAR(ISP)%LNSCATT(LN)=LN
          DO LN1=LN+1,LNX1
            IF(LOX(LN1,ISP).NE.LOX(LN,ISP)) CYCLE
            IF(ISCATT(LN1).GT.ISCATT(LN).AND.ISCATT(LN1).LE.0) THEN
              POTPAR(ISP)%LNSCATT(LN)=LN1
            END IF
          ENDDO
!         == DISTRIBUTE CHANNEL FOR SCATTERING WAVE FUNCTION =================
          DO LN1=LN,LNX1
            IF(LOX(LN1,ISP).NE.LOX(LN,ISP)) CYCLE
            POTPAR(ISP)%LNSCATT(LN1)=LN
          ENDDO
        ENDDO
        DEALLOCATE(ISCATT)
!       == MAP =================================================================
        DO LN=1,LNX1
          LN1=POTPAR(ISP)%LNSCATT(LN)
          IF(LN1.EQ.LN) CYCLE
          AEPHIDOT(:,LN)=AEPHIDOT(:,LN1)
          PSPHIDOT(:,LN)=PSPHIDOT(:,LN1)
          NLPHIDOT(:,LN)=NLPHIDOT(:,LN1)
        ENDDO
!
!       ========================================================================
!       ==  DETERMINE POTENTIAL PARAMETERS                                    ==
!       ========================================================================
        ALLOCATE(POTPAR(ISP)%QBAR(LNX1))
        ALLOCATE(POTPAR(ISP)%PHIDOTPROJ(LNX1))
        ALLOCATE(POTPAR(ISP)%KTOPHI(LNX1))
        ALLOCATE(POTPAR(ISP)%KTOPHIDOT(LNX1))
        ALLOCATE(POTPAR(ISP)%JBARTOPHIDOT(LNX1))
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
!         ====================================================================
!         == VALUE AND DERIVATIVE OF PARTIAL WAVES AND ENVELOPE FUNCTIONS   ==
!         == PHIDOT IS PHIBARDOT                                            ==
!         == THERE IS ONLY A SINGLE PHIDOT PER ANGULAR MOMENTUM             ==
!         ====================================================================
          CALL RADIAL$VALUE(GID,NR,NLPHIDOT(:,LN),RAD,PHIDOTVAL)
          CALL RADIAL$DERIVATIVE(GID,NR,NLPHIDOT(:,LN),RAD,PHIDOTDER)
          CALL RADIAL$VALUE(GID,NR,NLPHI(:,LN),RAD,PHIVAL)
          CALL RADIAL$DERIVATIVE(GID,NR,NLPHI(:,LN),RAD,PHIDER)
          CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
          CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
!
!         ====================================================================
!         == CALCULATE POTENTIAL PARAMETERS                                 ==
!         ====================================================================
          WJPHI=JVAL*PHIDER-JDER*PHIVAL
          WJPHIDOT=JVAL*PHIDOTDER-JDER*PHIDOTVAL
          WKPHI=KVAL*PHIDER-KDER*PHIVAL
          WKPHIDOT=KVAL*PHIDOTDER-KDER*PHIDOTVAL
          WPHIPHIDOT=PHIVAL*PHIDOTDER-PHIDER*PHIDOTVAL
          QBAR=WJPHIDOT/WKPHIDOT
          WJBARPHI=WJPHI-WKPHI*QBAR
!
!         ====================================================================
!         == K    -> |PHI>KTOPHI+|PHIBARDOT> KTOPHIDOT =======================
!         == JBAR ->             |PHIBARDOT> JBARTOPHIDOT ====================
!         ====================================================================
          POTPAR(ISP)%QBAR(LN)        =QBAR
          POTPAR(ISP)%KTOPHI(LN)      =WKPHIDOT/WPHIPHIDOT
          POTPAR(ISP)%KTOPHIDOT(LN)   =-WKPHI/WPHIPHIDOT
          POTPAR(ISP)%JBARTOPHIDOT(LN)=-WJBARPHI/WPHIPHIDOT
!
!         ==  <PRO|PSPHIDOT> =================================================
          CALL RADIAL$INTEGRAL(GID,NR,R**2*PRO(:,LN)*PSPHIDOT(:,LN),SVAR)
          POTPAR(ISP)%PHIDOTPROJ(LN)=SVAR
!         == CROSSCHECK BIOTHOGONALITY =======================================
          DO LN2=1,LNX(ISP)
            IF(LOX(LN2,ISP).NE.L) CYCLE
            CALL RADIAL$INTEGRAL(GID,NR,R**2*PRO(:,LN2)*PSPHI(:,LN),SVAR)
            IF(LN.EQ.LN2)SVAR=SVAR-1.D0
            IF(ABS(SVAR).GT.1.D-5) THEN
              CALL ERROR$MSG('VIOLATION OF BIORTHOGONALITY DETECTED')
              CALL ERROR$I4VAL('ISP',ISP)
              CALL ERROR$I4VAL('LN',LN)
              CALL ERROR$I4VAL('LN2',LN2)
              CALL ERROR$R8VAL('<P(LN2)|PHITILDE(LN)-DELTA(LN2,LN)',SVAR)
              CALL ERROR$STOP('LMTO_MAKEPOTPAR')
            END IF
          ENDDO
        ENDDO        
!
        DEALLOCATE(NLPHI)
        DEALLOCATE(AEPHI)
        DEALLOCATE(NLPHIDOT)
        DEALLOCATE(AEPHIDOT)
        DEALLOCATE(PSPHIDOT)
        DEALLOCATE(PSPHI)
        DEALLOCATE(PRO)
        DEALLOCATE(R)
        CALL SETUP$UNSELECT()
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_WRONSKITOPOTPAR(WKPHI,WKPHIDOT,WJPHI,WJPHIDOT &
     &                               ,WPHIPHIDOT,WJK,ENU,RAD &
     &                               ,QBAR,A,SQDELTABAR,CBAR)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: WKPHI
      REAL(8),INTENT(IN) :: WKPHIDOT
      REAL(8),INTENT(IN) :: WJPHI
      REAL(8),INTENT(IN) :: WJPHIDOT
      REAL(8),INTENT(IN) :: WPHIPHIDOT
      REAL(8),INTENT(IN) :: WJK
      REAL(8),INTENT(IN) :: ENU
      REAL(8),INTENT(IN) :: RAD
      REAL(8),INTENT(OUT):: QBAR
      REAL(8),INTENT(OUT):: CBAR
      REAL(8),INTENT(OUT):: A
      REAL(8),INTENT(OUT):: SQDELTABAR
      REAL(8)            :: WJBARPHI
      REAL(8)            :: WJBARK
      REAL(8)            :: DELTABAR
!     **************************************************************************
!     ==  SCREENING CHARGE ===============================================
      QBAR=WJPHIDOT/WKPHIDOT
      WJBARPHI=WJPHI-WKPHI*QBAR
      WJBARK=WJK ! IS INDEPENNDEN OF QBAR
!     == BAND CENTER =====================================================
      CBAR=ENU-WKPHI/WKPHIDOT
!     == BAND WIDTH ======================================================
      DELTABAR=WJBARPHI/WKPHIDOT
PRINT*,'W[JK] ',WJK,' =!=-1/RAD^2=',-1/RAD**2
PRINT*,'W[JBARK] ',WJBARK,' =!=',WJK
PRINT*,'W[PHI,PHIDOT] ',WPHIPHIDOT,' APPROX -<PHI|PHI>'
PRINT*,'CBAR-ENU ',CBAR-ENU,' CBAR ',CBAR,' ENU ',ENU
PRINT*,'DELTA (>0?) ',WJBARPHI/WKPHIDOT
PRINT*,'A^2 (>0?) ',WKPHIDOT*WJBARPHI/WPHIPHIDOT**2
PRINT*,'XX (>0?) ',WJBARPHI/WPHIPHIDOT
PRINT*,'TEST ',WJBARPHI/WKPHIDOT*WJBARK/WPHIPHIDOT
PRINT*,'W[KPHIDOT]*W[JBARPHI] ',WJBARPHI*WKPHIDOT &
       ,'=!= ',-WJBARK/WPHIPHIDOT
PRINT*,'W[KPHIDOT]/W[PHIPHIDOT] ',WKPHIDOT/WPHIPHIDOT
PRINT*,'W[JBARPHI]/W[PHIPHIDOT] ',WJBARPHI/WPHIPHIDOT

      IF(DELTABAR.LE.0.D0) THEN
        CALL ERROR$MSG('INTERNAL ERROR')
        CALL ERROR$MSG('DELTABAR MUST NOT BE NEGATIVE')
        CALL ERROR$STOP('LMTO_WRONSKITOPOTPAR')
      END IF
      SQDELTABAR=SQRT(DELTABAR)
!     == APPROXIMATE NORMALIZATION FACTOR=================================
      A=SQRT(2.D0*WJBARK/WPHIPHIDOT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_MAKETAILEDPARTIALWAVES()
!     **************************************************************************
!     ** CONSTRUCTS AUGMENTED HANKEL END BESSEL FUNCTIONS WITH                **
!     ** TWO EXPONENTIALS ATTACHED AT THE MATCHING RADIUS RAD                 **
!     **                                                                      **
!     ** HANKEL AND BESSEL FUNCTIONS ARE TREATED INDEPENDENTLY FORMING A      **
!     ** A LARGER ARRAY. AUGMENTED LMTOS ARE CONSTRUCTED AS SUPERPOSITION OF  **
!     ** HANKEL AND SCREENED BESSEL FUNCTIONS WITH ONSITE-ONLY STRUCTURE      **
!     ** CONSTANTS                                                            **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : K2,POTPAR,NSP,LNX,LOX,HYBRIDSETTING
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TCUT=.FALSE.
      REAL(8)                :: LAMBDA1
      REAL(8)                :: LAMBDA2
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)                :: RAD
      REAL(8)                :: AEZ
      REAL(8)                :: RCOV
      INTEGER(4)             :: IRAD ! FIRST POINT BEYOND RAD
      REAL(8)                :: QBAR
      REAL(8)                :: JVAL,JDER,KVAL,KDER
      REAL(8)                :: SVAR1,SVAR2,A1,A2,B1,B2
      INTEGER(4)             :: L
      INTEGER(4)             :: LRX,LMRX
      INTEGER(4)             :: LNXT
      INTEGER(4)             :: LMNXT
      INTEGER(4),ALLOCATABLE :: LOXT(:)
      INTEGER(4),ALLOCATABLE :: LNDOT(:)
      INTEGER(4),ALLOCATABLE :: LMNDOT(:)
      INTEGER(4)             :: ISP,LN,LN1,LN2,LNT,LMN,LMN1,LMN2,IM,IR
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: ULITTLE(:,:,:,:,:)
CHARACTER(128) :: STRING
REAL(8) :: PI
INTEGER(4)             :: LN3,LN4
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_MAKETAILEDPARTIALWAVES')
      DO ISP=1,NSP
        LAMBDA1=HYBRIDSETTING(ISP)%TAILEDLAMBDA1
        LAMBDA2=HYBRIDSETTING(ISP)%TAILEDLAMBDA2
!
        CALL SETUP$ISELECT(ISP)
!       == RADIAL GRID =========================================================
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        RAD=POTPAR(ISP)%RAD
        DO IR=1,NR
          IRAD=IR
          IF(R(IR).GT.RAD) EXIT
        ENDDO
        POTPAR(ISP)%TAILED%GID=GID
!       == DETERMINE MAPPING ===================================================
        LNXT=LNX(ISP)
        LMNXT=SUM(2*LOX(:LNX(ISP),ISP)+1)
        DO LN=1,LNX(ISP)
          IF(POTPAR(ISP)%LNSCATT(LN).EQ.LN) THEN
            LNXT=LNXT+1
            LMNXT=LMNXT+2*LOX(LN,ISP)+1
          END IF
        ENDDO
        POTPAR(ISP)%TAILED%LNX=LNXT
        POTPAR(ISP)%TAILED%LMNX=LMNXT
        ALLOCATE(LOXT(LNXT))
        ALLOCATE(LNDOT(LNXT))
        ALLOCATE(LMNDOT(LMNXT))
        LNT=LNX(ISP)
        DO LN=1,LNX(ISP)
          LOXT(:LN)=LOX(:LN,ISP)  ! FOR THE HANKEL FUNCTIONS
          IF(POTPAR(ISP)%LNSCATT(LN).NE.LN) CYCLE
          LNT=LNT+1
          LOXT(LNT)=LOX(LN,ISP)   ! FOR THE BESSEL FUNCTIONS
          LNDOT(LN)=LNT
          LNDOT(LNT)=LN
        ENDDO
        DO LN=1,LNX(ISP)          ! COMPLETE LNDOT ARRAY
          LNDOT(LN)=LNDOT(POTPAR(ISP)%LNSCATT(LN))
        ENDDO 
        ALLOCATE(POTPAR(ISP)%TAILED%LOX(LNXT))
        ALLOCATE(POTPAR(ISP)%TAILED%LNDOT(LNXT))
        POTPAR(ISP)%TAILED%LOX(:)=LOXT
        POTPAR(ISP)%TAILED%LNDOT=LNDOT
!
!       ========================================================================
!       == MAPPING "LMNDOT" FROM K TO JBAR FUNCTIONS AND VICE VERSA           ==
!       ==   JBAR(LMN)=F(LMNDOT(LMN)) ;                                       ==
!       ========================================================================
        LMN1=0
        LMN2=SUM(2*LOX(:LNX(ISP),ISP)+1)
        LN2=LNX(ISP)
        DO LN1=1,LNX(ISP)
          IF(POTPAR(ISP)%LNSCATT(LN1).EQ.LN1) THEN
            LN2=LN2+1
            DO IM=1,2*LOX(LN1,ISP)+1
              LMNDOT(LMN2+IM)=LMN1+IM
            ENDDO
            LMN=0
            DO LN=1,LNX(ISP)
              IF(POTPAR(ISP)%LNSCATT(LN).EQ.LN1) THEN
                DO IM=1,2*LOX(LN,ISP)+1
                  LMNDOT(LMN+IM)=LMN2+IM
                ENDDO
              END IF
              LMN=LMN+2*LOX(LN,ISP)+1
            ENDDO
            LMN2=LMN2+2*LOX(LN1,ISP)+1
          END IF
          LMN1=LMN1+2*LOX(LN1,ISP)+1
        ENDDO
        ALLOCATE(POTPAR(ISP)%TAILED%LMNDOT(LMNXT))
        POTPAR(ISP)%TAILED%LMNDOT=LMNDOT
!
!       ========================================================================
!       == AUGMENTED HANKEL AND BESSEL FUNCTIONS WITH TAILS ATTACHED          ==
!       ========================================================================
        ALLOCATE(POTPAR(ISP)%TAILED%AEF(NR,LNXT))
        ALLOCATE(POTPAR(ISP)%TAILED%PSF(NR,LNXT))
        ALLOCATE(POTPAR(ISP)%TAILED%NLF(NR,LNXT))
!
        ALLOCATE(AEPHI(NR,LNX(ISP)))
        ALLOCATE(AEPHIDOT(NR,LNX(ISP)))
        ALLOCATE(NLPHI(NR,LNX(ISP)))
        ALLOCATE(NLPHIDOT(NR,LNX(ISP)))
        ALLOCATE(PSPHI(NR,LNX(ISP)))
        ALLOCATE(PSPHIDOT(NR,LNX(ISP)))
        CALL SETUP$GETR8A('AEPHI',NR*LNX(ISP),AEPHI)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX(ISP),AEPHIDOT)
        CALL SETUP$GETR8A('QPHI',NR*LNX(ISP),NLPHI)
        CALL SETUP$GETR8A('QPHIDOT',NR*LNX(ISP),NLPHIDOT)
        CALL SETUP$GETR8A('PSPHI',NR*LNX(ISP),PSPHI)
        CALL SETUP$GETR8A('PSPHIDOT',NR*LNX(ISP),PSPHIDOT)
!
!       == TAIL PART ===========================================================
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          RAD=POTPAR(ISP)%RAD
          CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
          CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
!         -- TRANSFORM UNSCREENED BESSEL FUNCTION TO SCREENED BESSEL FUNCTION
          QBAR=POTPAR(ISP)%QBAR(LN)
          JVAL=JVAL-KVAL*QBAR
          JDER=JDER-KDER*QBAR
!         -- DETERMINE VALUE AND LOGARITHMIC DERIVATIVE OF PHI AND PHIBARDOT----
          A1=KVAL*(KDER/KVAL+LAMBDA2)/(LAMBDA2-LAMBDA1)
          A2=KVAL*(KDER/KVAL+LAMBDA1)/(LAMBDA1-LAMBDA2)
          B1=JVAL*(JDER/JVAL+LAMBDA2)/(LAMBDA2-LAMBDA1)
          B2=JVAL*(JDER/JVAL+LAMBDA1)/(LAMBDA1-LAMBDA2)
          DO IR=IRAD,NR
            SVAR1=EXP(-LAMBDA1*(R(IR)-RAD))
            SVAR2=EXP(-LAMBDA2*(R(IR)-RAD))
            POTPAR(ISP)%TAILED%NLF(IR,LN)       =A1*SVAR1+A2*SVAR2
            POTPAR(ISP)%TAILED%NLF(IR,LNDOT(LN))=B1*SVAR1+B2*SVAR2
          ENDDO
        ENDDO
!
!       == SPHERE PART =========================================================
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
          CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
          LN1=POTPAR(ISP)%LNSCATT(LN)
          A1=POTPAR(ISP)%KTOPHI(LN)
          A2=POTPAR(ISP)%KTOPHIDOT(LN)
          POTPAR(ISP)%TAILED%NLF(:IRAD-1,LN)=NLPHI(:IRAD-1,LN)*A1 &
       &                                    +NLPHIDOT(:IRAD-1,LN1)*A2
!         == THE COMPLEX ADDITION OF DIFFERENCES IN THE FOLLOWING IS ===========
!         == NECESSARY, IF THE AE, PS AND NL PARTIAL WAVES DIFFER AT THE =======
!         == MATCHING RADIUS DUE TO THE ADMIXED TAILS OF CORE STATES ===========
          POTPAR(ISP)%TAILED%AEF(:,LN)=POTPAR(ISP)%TAILED%NLF(:,LN) &
       &      +(AEPHI(:,LN)-NLPHI(:,LN))*A1+(AEPHIDOT(:,LN1)-NLPHIDOT(:,LN1))*A2
          POTPAR(ISP)%TAILED%PSF(:,LN)=POTPAR(ISP)%TAILED%NLF(:,LN) &
          &   +(PSPHI(:,LN)-NLPHI(:,LN))*A1+(PSPHIDOT(:,LN1)-NLPHIDOT(:,LN1))*A2
!
          IF(LN.EQ.LN1) THEN
            A2=POTPAR(ISP)%JBARTOPHIDOT(LN)
            POTPAR(ISP)%TAILED%NLF(:IRAD-1,LNDOT(LN))=NLPHIDOT(:IRAD-1,LN1)*A2
!           == THE COMPLEX ADDITION OF DIFFERENCES IN THE FOLLOWING IS =========
!           == NECESSARY, IF THE AE, PS AND NL PARTIAL WAVES DIFFER AT THE =====
!           == MATCHING RADIUS DUE TO THE ADMIXED TAILS OF CORE STATES =========
            POTPAR(ISP)%TAILED%AEF(:,LNDOT(LN)) &
      &                           =POTPAR(ISP)%TAILED%NLF(:,LNDOT(LN)) &
      &                           +(AEPHIDOT(:,LN1)-NLPHIDOT(:,LN1))*A2
            POTPAR(ISP)%TAILED%PSF(:,LNDOT(LN)) &
      &                           =POTPAR(ISP)%TAILED%NLF(:,LNDOT(LN)) &
      &                           +(PSPHIDOT(:,LN1)-NLPHIDOT(:,LN1))*A2
!
          END IF
        ENDDO
        DEALLOCATE(AEPHI)
        DEALLOCATE(AEPHIDOT)
        DEALLOCATE(PSPHI)
        DEALLOCATE(PSPHIDOT)
        DEALLOCATE(NLPHI)
        DEALLOCATE(NLPHIDOT)
!
!       ========================================================================
!       ==  CUT OFF TAILS FOR STABILITY                                       ==
!       ========================================================================
        IF(TCUT) THEN
          CALL SETUP$GETR8('AEZ',AEZ)
          CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
          DO IR=1,NR
            IF(R(IR).GT.2.D0*RCOV) THEN  ! 2*RCOV IS A BIT ARBITRARY
              POTPAR(ISP)%TAILED%NLF(IR:,:)=0.D0
              POTPAR(ISP)%TAILED%AEF(IR:,:)=0.D0
              POTPAR(ISP)%TAILED%PSF(IR:,:)=0.D0
              EXIT
            END IF
          ENDDO
        END IF

!!$PRINT*,'WARNING!!!! FUDGE FOR TESTING H2'
!!$PI=4.D0*ATAN(1.D0)
!!$POTPAR(ISP)%TAILED%AEF(:,:)=0.D0
!!$SVAR2=0.168856D0
!!$SVAR1=0.444635D0*(2.D0*SVAR2/PI)**0.75D0
!!$POTPAR(ISP)%TAILED%AEF(:,1)=POTPAR(ISP)%TAILED%AEF(:,1) &
!!$ &     +SVAR1*EXP(-SVAR2*R(:)**2)
!!$SVAR2=0.623913D0
!!$SVAR1=0.535328D0*(2.D0*SVAR2/PI)**0.75D0
!!$POTPAR(ISP)%TAILED%AEF(:,1)=POTPAR(ISP)%TAILED%AEF(:,1) &
!!$ &     +SVAR1*EXP(-SVAR2*R(:)**2)
!!$SVAR2=3.42525D0
!!$SVAR1=0.154329D0*(2.D0*SVAR2/PI)**0.75D0
!!$POTPAR(ISP)%TAILED%AEF(:,1)=POTPAR(ISP)%TAILED%AEF(:,1) &
!!$ &     +SVAR1*EXP(-SVAR2*R(:)**2)
!!$      ! REMOVE SPHERICAL HARMONICS
!!$POTPAR(ISP)%TAILED%AEF(:,1)=POTPAR(ISP)%TAILED%AEF(:,1)*SQRT(4.D0*PI)
!!$     ! COPY INTO OTHER ARRAYS
!!$POTPAR(ISP)%TAILED%NLF(:,:)=POTPAR(ISP)%TAILED%AEF(:,:)
!!$POTPAR(ISP)%TAILED%PSF(:,:)=POTPAR(ISP)%TAILED%AEF(:,:)
!
!
!
!!$WRITE(STRING,FMT='(I5)')ISP
!!$STRING='AETAILS_FORATOMTYPE'//TRIM(ADJUSTL(STRING))//'.DAT'
!!$CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LNXT,POTPAR(ISP)%TAILED%AEF)
!!$WRITE(STRING,FMT='(I5)')ISP
!!$STRING='NLTAILS_FORATOMTYPE'//TRIM(ADJUSTL(STRING))//'.DAT'
!!$CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LNXT,POTPAR(ISP)%TAILED%NLF)
!
        DEALLOCATE(R)
!       ========================================================================
!       == ONSITE U-TENSOR OF TAILED PARTIAL WAVES                            ==
!       ========================================================================
        CALL SETUP$GETI4('LMRX',LMRX)
        LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
        ALLOCATE(POTPAR(ISP)%TAILED%U(LMNXT,LMNXT,LMNXT,LMNXT))
        ALLOCATE(ULITTLE(LRX+1,LNXT,LNXT,LNXT,LNXT))
        CALL LMTO_ULITTLE(GID,NR,LRX,LNXT,LOXT,POTPAR(ISP)%TAILED%AEF,ULITTLE)
!!$DO L=1,LRX+1
!!$  DO LN1=1,LNXT
!!$    DO LN2=LN1,LNXT
!!$      DO LN3=1,LNXT
!!$        DO LN4=LN3,LNXT
!!$          IF(ABS(ULITTLE(L,LN1,LN2,LN3,LN4)).LT.1.D-5) CYCLE
!!$          WRITE(*,FMT='(5I4,F10.5)')L,LN1,LN2,LN3,LN4,ULITTLE(L,LN1,LN2,LN3,LN4)
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
!!$STOP 'FORCED'
        CALL LMTO_UTENSOR(LRX,LMNXT,LNXT,LOXT,ULITTLE,POTPAR(ISP)%TAILED%U)
        DEALLOCATE(ULITTLE)
!
!       ========================================================================
!       == ONSITE OVERLAP MATRIX                                              ==
!       ========================================================================
        ALLOCATE(POTPAR(ISP)%TAILED%OVERLAP(LMNXT,LMNXT))
        CALL LMTO_ONECENTEROVERLAP(GID,NR,LNXT,LOXT,POTPAR(ISP)%TAILED%AEF &
     &                            ,LMNXT,POTPAR(ISP)%TAILED%OVERLAP)
        ALLOCATE(POTPAR(ISP)%TAILED%QLN(2,LNXT,LNXT))
        CALL LMTO_ONECENTERQLN(GID,NR,LNXT,LOXT,POTPAR(ISP)%TAILED%AEF &
     &                            ,POTPAR(ISP)%TAILED%QLN)
        DEALLOCATE(LNDOT)
        DEALLOCATE(LMNDOT)
        DEALLOCATE(LOXT)
        CALL SETUP$UNSELECT()
      ENDDO
                              CALL TRACE$POP() 
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_ULITTLE(GID,NR,LRX,LNX,LOX,CHI,ULITTLE)
!     **                                                                      **
!     ** SLATER INTEGRALS.                                                    **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LRX
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: CHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: ULITTLE(LRX+1,LNX,LNX,LNX,LNX)
      INTEGER(4)            :: LN1,LN2,LN3,LN4
      INTEGER(4)            :: L
      INTEGER(4)            :: LMIN,LMAX,ISVAR1,ISVAR2
      REAL(8)               :: RHO(NR)
      REAL(8)               :: POT(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR
      REAL(8)               :: R(NR)
      REAL(8)               :: PI,FOURPI
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_ULITTLE')
      CALL RADIAL$R(GID,NR,R)
      ULITTLE=0.D0
      DO LN1=1,LNX
        DO LN2=LN1,LNX
          RHO(:)=CHI(:,LN1)*CHI(:,LN2)
!         == USE SELECTION RULE ================================================
!         == (NOT TO SAVE TIME HERE, BUT LATER FOR THE U-TENSOR) ===============
          ISVAR1=ABS(LOX(LN1)+LOX(LN2))  
          ISVAR2=ABS(LOX(LN1)-LOX(LN2))
          LMIN=MIN(ISVAR1,ISVAR2)
          LMAX=MAX(ISVAR1,ISVAR2)
          LMAX=MIN(LMAX,LRX)
          DO L=LMIN,LMAX
            CALL RADIAL$POISSON(GID,NR,L,RHO,POT)
            POT(:)=POT(:)*R(:)**2
            DO LN3=1,LNX
              DO LN4=LN3,LNX
                ISVAR1=ABS(LOX(LN3)+LOX(LN4))
                ISVAR2=ABS(LOX(LN3)-LOX(LN4))
                IF(L.LT.MIN(ISVAR1,ISVAR2)) CYCLE
                IF(L.GT.MAX(ISVAR1,ISVAR2)) CYCLE
                AUX(:)=CHI(:,LN3)*CHI(:,LN4)*POT(:)
                CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
                ULITTLE(L+1,LN1,LN2,LN3,LN4)=SVAR
                ULITTLE(L+1,LN2,LN1,LN3,LN4)=SVAR
                ULITTLE(L+1,LN1,LN2,LN4,LN3)=SVAR
                ULITTLE(L+1,LN2,LN1,LN4,LN3)=SVAR
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == ADD FACTOR CONSISTENT WITH DEFINITION OF SLATER INTEGRALS            ==
!     ==========================================================================
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      DO L=0,LRX
        ULITTLE(L+1,:,:,:,:)=ULITTLE(L+1,:,:,:,:)*REAL(2*L+1,KIND=8)/FOURPI
      ENDDO

                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_UTENSOR(LRX,NORB,LNX,LOX,ULITTLE,U)
!     **************************************************************************
!     ** EXPANDS SLATER INTEGRALS FROM LMTO_ULITTLE TO THE FULL U-TENSOR      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LRX
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: ULITTLE(LRX+1,LNX,LNX,LNX,LNX)
      INTEGER(4),INTENT(IN) :: NORB
      REAL(8)   ,INTENT(OUT):: U(NORB,NORB,NORB,NORB)
      INTEGER(4)            :: LN1,LN2,LN3,LN4
      INTEGER(4)            :: IORB1,IORB2,IORB3,IORB4
      INTEGER(4)            :: L1,L2,L3,L4
      INTEGER(4)            :: LM1,LM2,LM3,LM4
      INTEGER(4)            :: M1,M2,M3,M4
      INTEGER(4)            :: L,M,LM,LX
      REAL(8)               :: CG1,CG2
      REAL(8)               :: SVAR
      REAL(8)               :: PI,FOURPI
      REAL(8)               :: FOURPIBY2LPLUS1
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_UTENSOR')
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
!
      U(:,:,:,:)=0.D0
      IORB1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LM1=L1**2
        DO M1=1,2*L1+1
          IORB1=IORB1+1
          LM1=LM1+1
!
          IORB2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            LM2=L2**2
            DO M2=1,2*L2+1
              IORB2=IORB2+1
              LM2=LM2+1
!
              IORB3=0
              DO LN3=1,LNX
                L3=LOX(LN3)
                LM3=L3**2
                DO M3=1,2*L3+1
                  IORB3=IORB3+1
                  LM3=LM3+1
                  IF(LM3.LT.LM1) CYCLE
!
                  IORB4=0
                  DO LN4=1,LNX
                    L4=LOX(LN4)
                    LM4=L4**2
                    DO M4=1,2*L4+1
                      IORB4=IORB4+1
                      LM4=LM4+1
                      IF(LM4.LT.LM2) CYCLE
!         
                      IF(MAXVAL(ABS(ULITTLE(:,LN2,LN4,LN3,LN1))).EQ.0.D0) CYCLE
                      LX=MIN(LRX,L2+L4,L1+L3)
                      SVAR=0.D0
                      LM=0
                      DO L=0,LX
                        FOURPIBY2LPLUS1=FOURPI/REAL(2*L+1,KIND=8)
                        DO M=1,2*L+1
                          LM=LM+1
                          CALL CLEBSCH(LM2,LM4,LM,CG1)
                          CALL CLEBSCH(LM3,LM1,LM,CG2)
                          SVAR=SVAR+FOURPIBY2LPLUS1*CG1*CG2 &
    &                                          *ULITTLE(L+1,LN2,LN4,LN3,LN1)
                        ENDDO
                      ENDDO
                      U(IORB1,IORB2,IORB3,IORB4)=SVAR
                      U(IORB1,IORB4,IORB3,IORB2)=SVAR
                      U(IORB3,IORB2,IORB1,IORB4)=SVAR
                      U(IORB3,IORB4,IORB1,IORB2)=SVAR
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_ONECENTEROVERLAP(GID,NR,LNX,LOX,CHI,LMNX,OVERLAP)
!     **                                                                      **
!     ** SLATER INTEGRALS.                                                    **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LMNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: CHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: OVERLAP(LMNX,LMNX)
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: LMN10,LMN20
      INTEGER(4)            :: L1,L2,IM
      REAL(8)               :: AUX(NR),SVAR
      REAL(8)               :: R(NR)
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_ONECENTEROVERLAP')
      CALL RADIAL$R(GID,NR,R)
      OVERLAP(:,:)=0.D0
      LMN10=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LMN20=LMN10
        DO LN2=LN1,LNX
          L2=LOX(LN2)
          IF(L1.EQ.L2) THEN
            AUX(:)=R(:)**2*CHI(:,LN1)*CHI(:,LN2)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            DO IM=1,2*L1+1
              OVERLAP(LMN10+IM,LMN20+IM)=SVAR
              OVERLAP(LMN20+IM,LMN10+IM)=SVAR
            ENDDO
          END IF
          LMN20=LMN20+2*L2+1
        ENDDO
        LMN10=LMN10+2*L1+1
      ENDDO
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_ONECENTERQLN(GID,NR,LNX,LOX,CHI,QLN)
!     **************************************************************************
!     **  DETERMINES THE MATRIX ELEMENTS OF MONOPOLE AND DIPOLE BETWEEN       **
!     **  ORBITALS WITH PURE ANGULAR MOMENTUM CHARACTER                       **
!     **                                                                      **
!     **  FOR A DENSITY CHI_LN(R)Y_L(R)RHO_LN'(|R|)*Y_L'(R)                   **
!     **  THE CHARGE IS QLN(1,LN,LN')C_{L,L',S} AND                           **
!     **  THE DIPOLE IS QLN(2,LN,LN')C_{L,L',P}                               ** 
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: CHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: QLN(2,LNX,LNX)
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: L1,L2,IM
      REAL(8)               :: AUX(NR),SVAR
      REAL(8)               :: R(NR)
      REAL(8)               :: PI,SQ4PI,SQ4PITHIRD
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_ONECENTERMULTIPOLE')
      PI=4.D0*ATAN(1.D0)
      SQ4PI=SQRT(4.D0*PI)
      SQ4PITHIRD=SQRT(4.D0*PI/3.D0)
      CALL RADIAL$R(GID,NR,R)
      QLN(:,:,:)=0.D0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO LN2=LN1,LNX
          L2=LOX(LN2)
          IF(ABS(L1-L2).EQ.0.D0) THEN
            AUX(:)=R(:)**2*CHI(:,LN1)*CHI(:,LN2)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            SVAR=SQ4PI*SVAR
            QLN(1,LN1,LN2)=SVAR
            QLN(1,LN2,LN1)=SVAR
          ELSE IF(ABS(L1-L2).EQ.1.D0) THEN
            AUX(:)=R(:)**3*CHI(:,LN1)*CHI(:,LN2)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            SVAR=SQ4PITHIRD*SVAR
            QLN(2,LN1,LN2)=SVAR
            QLN(2,LN2,LN1)=SVAR
          END IF
        ENDDO
      ENDDO
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILEDGAUSSFIT()
!     **************************************************************************
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : POTPAR,NSP
      IMPLICIT NONE
      INTEGER(4),PARAMETER   :: NPOWPAR=4   !X#(POWERS), R^(L+2N)
      INTEGER(4),PARAMETER   :: NX=2*(NPOWPAR-1) !HIGHEST POWER 
      INTEGER(4),PARAMETER   :: NEPAR=2    !#(GAUSS-EXPONENTS)
      REAL(8)   ,PARAMETER   :: R1PAR=1.D0
      REAL(8)   ,PARAMETER   :: FACPAR=5.D0
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      INTEGER(4)             :: NE
      INTEGER(4)             :: NPOW
      INTEGER(4)             :: NPOW2
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:) !(LNX)
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)   ,ALLOCATABLE :: R(:) !(NR) RADIAL GRID
      REAL(8)   ,ALLOCATABLE :: W(:) !(NR)
      REAL(8)   ,ALLOCATABLE :: AUX(:) !(NR)
      REAL(8)   ,ALLOCATABLE :: E(:) !(NE)
      REAL(8)   ,ALLOCATABLE :: C(:,:,:) !(NPOW,NE,LNX)
      REAL(8)                :: SVAR,RI
      INTEGER(4)             :: ISP,IE,LN,L,IR,I
      INTEGER(4)             :: NFIL
!     **************************************************************************
                                CALL TRACE$PUSH('LMTO_TAILEDGAUSSFIT')
      DO ISP=1,NSP
        LNX=POTPAR(ISP)%TAILED%LNX
        ALLOCATE(LOX(LNX))
        LOX=POTPAR(ISP)%TAILED%LOX
!
        GID=POTPAR(ISP)%TAILED%GID
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        ALLOCATE(AUX(NR))
        ALLOCATE(W(NR))
        W(:)=R(:)**2
! 
        NE=NEPAR
        NPOW=NPOWPAR
        ALLOCATE(E(NE))
        DO IE=1,NE
          E(IE)=1.D0/(R1PAR*FACPAR**(IE-1))
        ENDDO
!
        ALLOCATE(C(NPOW,NE,LNX))
        C(:,:,:)=0.D0
        DO LN=1,LNX
          L=LOX(LN)
          AUX=POTPAR(ISP)%TAILED%NLF(:,LN)
          NPOW2=INT(0.5D0*REAL(NX-L)+1.000001D0)  !R^(L+2N), L+2N=0,NPOW-1
          IF(NPOW2.LT.1) CYCLE
          CALL GAUSSIAN_FITGAUSS(GID,NR,W,L,AUX,NE,NPOW2,E,C(:NPOW2,:,LN))
        ENDDO
        ALLOCATE(POTPAR(ISP)%TAILED%GAUSSNLF%E(NE))
        ALLOCATE(POTPAR(ISP)%TAILED%GAUSSNLF%C(NPOW,NE,LNX))
        POTPAR(ISP)%TAILED%GAUSSNLF%NIJK=NPOW
        POTPAR(ISP)%TAILED%GAUSSNLF%NE=NE
        POTPAR(ISP)%TAILED%GAUSSNLF%E=E
        POTPAR(ISP)%TAILED%GAUSSNLF%C=C
        DEALLOCATE(E)
        DEALLOCATE(C)
        DEALLOCATE(LOX)
        DEALLOCATE(R)
        DEALLOCATE(W)
        DEALLOCATE(AUX)
      ENDDO
!
!     ==========================================================================
!     ==  WRITE FIT TO FILE FOR COMPARISON                                    ==
!     ==========================================================================
      IF(TPR) THEN
        NFIL=11
        DO ISP=1,NSP
          OPEN(UNIT=NFIL,FILE='FIT.DAT')
          LNX=POTPAR(ISP)%TAILED%LNX
          ALLOCATE(LOX(LNX))
          LOX=POTPAR(ISP)%TAILED%LOX
!
          GID=POTPAR(ISP)%TAILED%GID
          CALL RADIAL$GETI4(GID,'NR',NR)
          ALLOCATE(R(NR))
          CALL RADIAL$R(GID,NR,R)
          NPOW=POTPAR(ISP)%TAILED%GAUSSNLF%NIJK
          ALLOCATE(W(LNX))
          DO IR=1,NR
            RI=R(IR)
            W(:)=0.D0
            DO IE=1,NE
              SVAR=EXP(-POTPAR(ISP)%TAILED%GAUSSNLF%E(IE)*RI**2)
              DO I=1,NPOW
                W(:)=W(:)+POTPAR(ISP)%TAILED%GAUSSNLF%C(I,IE,:)*SVAR
                SVAR=SVAR*RI**2
              ENDDO
            ENDDO
            W(:)=W(:)*RI**LOX(:)
            WRITE(NFIL,*)RI,W(:),POTPAR(ISP)%TAILED%NLF(IR,:)
          ENDDO
          CLOSE(NFIL)
          DEALLOCATE(LOX)
          DEALLOCATE(W)
          DEALLOCATE(R)
        ENDDO
        CALL ERROR$MSG('REGULAR STOP AFTER WRITING FILE FIT.DAT')
        CALL ERROR$STOP('LMTO_TAILEDGAUSSFIT')
      END IF
!
!     ==========================================================================
!     ==  MULTIPLY WITH SPHERICAL HARMONICS                                   ==
!     ==========================================================================
      CALL LMTO_TAILEDGAUSSORBTOYLM()
                                CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILEDGAUSSORBTOYLM()
!     **************************************************************************
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : POTPAR,NSP
      IMPLICIT NONE
      INTEGER(4)             :: LX         !X(ANGULAR MOMENTUM)
      INTEGER(4)             :: NPOWX      !X#(DOUBLE POWERS)
      INTEGER(4)             :: NPOW       !#(DOUBLE POWERS)
      INTEGER(4)             :: NIJKX      !
      INTEGER(4)             :: NIJK
      INTEGER(4)             :: NE
      INTEGER(4)             :: NX             !HIGHEST POWER, DETERMINES NIJK
      REAL(8)   ,ALLOCATABLE :: YLMPOL(:,:)
      REAL(8)   ,ALLOCATABLE :: POLYLM(:,:,:)
      REAL(8)   ,ALLOCATABLE :: CRAD(:,:,:)
      REAL(8)   ,ALLOCATABLE :: ORB(:,:,:)
      REAL(8)                :: SVAR1,SVAR2
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4)             :: ISP,N,I,J,K,I1,J1,K1,I2,J2,K2,IE
      INTEGER(4)             :: LN,LM,LMN,IM,IND,IND1
      INTEGER(4)             :: INDX
      INTEGER(4)             :: L
      INTEGER(4)             :: LMX
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO_TAILEDGAUSSORBTOYLM')
!
!     ==========================================================================
!     ==  DETERMINE DIMENSIONS                                                ==
!     ==========================================================================
      LX=-1
      NPOWX=-1
      NX=0
      DO ISP=1,NSP
        L=MAXVAL(POTPAR(ISP)%TAILED%LOX(:))
        NPOW=POTPAR(ISP)%TAILED%GAUSSNLF%NIJK  !#(DOUBLE POWERS)
        NPOWX=MAX(NPOW,NPOWX)
        LX=MAX(LX,L)
        NX=MAX(NX,2*(NPOW-1))   !HIGHEST POWER
      ENDDO
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKX,NX,0,0)
!
!     ==========================================================================
!     == CALCULATE YLMPOL: R^L*Y_LM = SUM_IJK X^IY^JZ^K * YLMPOL(IJK,LM)      ==
!     ==========================================================================
      INDX=(LX+1)*(LX+2)*(LX+3)/6
      LMX=(LX+1)**2
      ALLOCATE(YLMPOL(INDX,LMX))
      CALL GAUSSIAN_YLMPOL(LX,YLMPOL)
!!$PRINT*,'YLMPOL ',YLMPOL
!!$PRINT*,'INDX   ',INDX
!!$PRINT*,'NPOWX  ',NPOWX
!!$PRINT*,'LMX    ',LMX
!
!     ==========================================================================
!     == CALCULATE POLYLM(N+1,LM)=|R|**(2N)*Y(LM)      N=0,...,NPOWX-1        ==
!     ==========================================================================
      ALLOCATE(POLYLM(NIJKX,NPOWX,LMX))
      POLYLM(:,:,:)=0.D0
      DO L=0,LX
        DO IM=1,2*L+1
          LM=L**2+IM
          DO IND1=1,INDX
            IF(YLMPOL(IND1,LM).EQ.0.D0) CYCLE
            CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND1,I1,J1,K1)
            NPOW=INT(0.5D0*REAL(NX-L)+1.000001D0)
!           == MULTIPLY WITH (X^2+Y^2+Z^2)^N = SUM_{I+J+K=N}:                 ==
!           ==          :BINOM(I+J+K;J+K) * BINOM(J+K;K) * X^2I * Y^2J * Z^2K ==
            DO K2=0,NPOW-1
              K=K1+2*K2
              DO J2=0,NPOW-1
                J=J1+2*J2
                CALL BINOMIALCOEFFICIENT(J2+K2,K2,SVAR2)
                DO I2=0,NPOW-1
                  I=I1+2*I2
                  IF(I+J+K.GT.NX) CYCLE
                  CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND,I,J,K)
                  IF(IND.GT.NIJKX) CYCLE
                  CALL BINOMIALCOEFFICIENT(I2+J2+K2,J2+K2,SVAR1)
                  N=I2+J2+K2
                  POLYLM(IND,N+1,LM)=POLYLM(IND,N+1,LM) &
     &                              +SVAR1*SVAR2*YLMPOL(IND1,LM)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO ! END OF LOOP OVER ORBITALS (IM)
      ENDDO  ! END OF LOOP OVER L
!PRINT*,'POLYLM ',MAXVAL(ABS(POLYLM))
!!$DO IND1=1,NIJKX
!!$  CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND1,I1,J1,K1)
!!$  PRINT*,'POLYLM ',IND1,I1,J1,K1,POLYLM(IND1,:,:)
!!$ENDDO
!
!     ==========================================================================
!     == CALCULATE ORBITALS                                                   ==
!     ==========================================================================
      DO ISP=1,NSP
        LNX=POTPAR(ISP)%TAILED%LNX
        ALLOCATE(LOX(LNX))
        LOX=POTPAR(ISP)%TAILED%LOX
        NPOW=POTPAR(ISP)%TAILED%GAUSSNLF%NIJK
        NE=POTPAR(ISP)%TAILED%GAUSSNLF%NE
        ALLOCATE(CRAD(NPOW,NE,LNX))
        CRAD(:,:,:)=POTPAR(ISP)%TAILED%GAUSSNLF%C
        DEALLOCATE(POTPAR(ISP)%TAILED%GAUSSNLF%C)  ! WILL BE OVERWRITTEN
        CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJK,2*(NPOW-1),0,0)
        LMNX=SUM(2*LOX+1)
!
        ALLOCATE(ORB(NIJK,NE,LMNX))
        ORB(:,:,:)=0.D0
        LMN=0
        DO LN=1,LNX
          L=LOX(LN)
          LM=L**2
          DO IM=1,2*L+1
            LMN=LMN+1
            LM=LM+1
            DO IE=1,NE
              DO N=0,NPOW-1
                ORB(:,IE,LMN)=ORB(:,IE,LMN)+POLYLM(:NIJK,N+1,LM)*CRAD(N+1,IE,LN)
              ENDDO
            ENDDO
          ENDDO
        ENDDO              
        POTPAR(ISP)%TAILED%GAUSSNLF%NIJK=NIJK
        POTPAR(ISP)%TAILED%GAUSSNLF%NORB=LMNX
        ALLOCATE(POTPAR(ISP)%TAILED%GAUSSNLF%C(NIJK,NE,LMNX))
        POTPAR(ISP)%TAILED%GAUSSNLF%C=ORB
        DEALLOCATE(LOX)
        DEALLOCATE(CRAD)
        DEALLOCATE(ORB)
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILEDGAUSSOFFSITEU()
!     **************************************************************************
!     ** CALCULATE U-TENSOR MATRIX ELEMENTS OF TWO ORBITALS ON ONE ATOM       **
!     ** AND TWO ORBITALS ON THE OTHER SIDE. THE U-TENSOR IS THE INTERACTION  **
!     ** OF THE OVERLAPS IN THE BOND CENTER                                   **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : POTPAR,NSP,OFFSITEX
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: ISPA,ISPB
      INTEGER(4)             :: NIJKA,NIJKB
      INTEGER(4)             :: LMNXA,LMNXB
      INTEGER(4)             :: NEA,NEB
      REAL(8)   ,ALLOCATABLE :: EA(:),EB(:)
      REAL(8)   ,ALLOCATABLE :: ORBA(:,:,:),ORBB(:,:,:)
      REAL(8)   ,ALLOCATABLE :: UABCD(:,:,:,:)
      REAL(8)                :: DIS
      INTEGER(4)             :: LMN1,LMN2,LMN3,LMN4
      INTEGER(4)             :: IA,IB,IC,ID,IAB,ICD
      INTEGER(4)             :: IE,I,IND,IDIS,ISVAR
      REAL(8)                :: SVAR
      INTEGER(4)             :: NDIS
REAL(8) :: PI
REAL(8) :: E1,E2,E3,X1,X2,X3,EAB,SVAR0,SVAR1,SVAR2,SVAR3,SVAR4,RA(3),RB(3)
INTEGER(4) :: J,K,L
      INTEGER(4)             :: NTASKS,THISTASK,COUNT
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      COUNT=0
      DO ISPA=1,NSP
        NIJKA=POTPAR(ISPA)%TAILED%GAUSSNLF%NIJK
        LMNXA=POTPAR(ISPA)%TAILED%GAUSSNLF%NORB
        NEA  =POTPAR(ISPA)%TAILED%GAUSSNLF%NE
        ALLOCATE(EA(NEA))
        EA(:)=POTPAR(ISPA)%TAILED%GAUSSNLF%E(:)
        ALLOCATE(ORBA(NIJKA,NEA,LMNXA))
        ORBA(:,:,:)=POTPAR(ISPA)%TAILED%GAUSSNLF%C(:,:,:)
        DO ISPB=1,NSP
          NIJKB=POTPAR(ISPB)%TAILED%GAUSSNLF%NIJK
          LMNXB=POTPAR(ISPB)%TAILED%GAUSSNLF%NORB
          NEB  =POTPAR(ISPB)%TAILED%GAUSSNLF%NE
          ALLOCATE(EB(NEB))
          EB(:)=POTPAR(ISPB)%TAILED%GAUSSNLF%E(:)
          ALLOCATE(ORBB(NIJKB,NEB,LMNXB))
          ORBB(:,:,:)=POTPAR(ISPB)%TAILED%GAUSSNLF%C(:,:,:)
!
!         ======================================================================
!         == LOOP OVER DISTANCES                                              ==
!         ======================================================================
          NDIS=OFFSITEX(ISPA,ISPB)%NDIS
          ISVAR=LMNXA*LMNXB
          ISVAR=NINT(0.5D0*REAL(ISVAR*(ISVAR+1),KIND=8))
          ALLOCATE(OFFSITEX(ISPA,ISPB)%BONDU(NDIS,ISVAR))
          ALLOCATE(UABCD(LMNXA,LMNXB,LMNXA,LMNXB))
          OFFSITEX(ISPA,ISPB)%BONDU(:,:)=0.D0
          DO IDIS=1,NDIS
            COUNT=COUNT+1
            IF(MOD(COUNT-1,NTASKS).NE.THISTASK-1) CYCLE
!
!           ====================================================================
!           == FOURCENTER MATRIX ELEMENTS                                     ==
!           == THE RESULT OF GAUSSIAN$ZDIRECTION_FOURCENTER IS DEFINED AS     ==
!           == U(1,2,3,4)=INT DX IT DX': A1(X)*B2(X)][A3(X')*B4(X')]/|R-R'|   ==
!           == THE ORDER OF INDICES DEVIATES FROM THE U-TENSOR CONVENTION     ==
!           ====================================================================
            DIS=OFFSITEX(ISPA,ISPB)%DIS(IDIS)
            CALL GAUSSIAN$ZDIRECTION_FOURCENTER(NIJKA,NEA,EA,LMNXA,ORBA &
     &                                         ,NIJKB,NEB,EB,LMNXB,ORBB &
     &                                         ,DIS,UABCD)
!
!           ====================================================================
!           == MAP ONTO OFFSITEX STRUCTURE                                    ==
!           ====================================================================
            IND=0
            DO LMN2=1,LMNXB
              DO LMN1=1,LMNXA
                IAB=LMN1+LMNXA*(LMN2-1)
                DO LMN4=1,LMNXB
                  DO LMN3=1,LMNXA
                    ICD=LMN3+LMNXA*(LMN4-1)
                    IF(ICD.GT.IAB) EXIT
                    IND=IND+1
                    OFFSITEX(ISPA,ISPB)%BONDU(IDIS,IND) &
      &                                =UABCD(LMN1,LMN2,LMN3,LMN4)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO  ! END OF LOOP OVER DISTANCES
          DEALLOCATE(UABCD)
          DEALLOCATE(EB)
          DEALLOCATE(ORBB)
        ENDDO !END OF LOOP OVER SECOND ATOM TYPE ISPA
        DEALLOCATE(EA)
        DEALLOCATE(ORBA)
      ENDDO   !END OF LOOP OVER FIRST ATOM TYPE ISPA
!
!     ==========================================================================
!     ==  COMBINE RESULTS                                                     ==
!     ==========================================================================
      DO ISPA=1,NSP
        DO ISPB=1,NSP
          CALL MPE$COMBINE('MONOMER','+',OFFSITEX(ISPA,ISPB)%BONDU)
        ENDDO
      ENDDO      
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILEDPRODUCTS()
!     **************************************************************************
!     **  PRODRHO IS THE PRODUCT OF RADIAL FUNCTIONS OF TWO DIFFERENT         **
!     **  ANGULAR MOMENTA EXPANDED IN RADIAL GAUSSIANS R^(L+2N)*E(-E*^2).     **
!     **                                                                      **
!     **  - IF R1PAR IS TOO SMALL, THE LONG TAILS ARE NOT PRESENTS AND        **
!     **    GAUSS OSCILLATIONS OCCUR AT SHORTER DISTANCES                     **
!     **                                                                      **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : POTPAR,NSP
      IMPLICIT NONE
      INTEGER(4),PARAMETER   :: NEPAR=3     !#(GAUSS-EXPONENTS)
      INTEGER(4),PARAMETER   :: NPOWPAR=4   !X#(POWERS), R^(L+2N)
      INTEGER(4),PARAMETER   :: NX=2*(NPOWPAR-1) !HIGHEST POWER 
      REAL(8)   ,PARAMETER   :: R1PAR=1.D0
      REAL(8)   ,PARAMETER   :: FACPAR=2.0D0
      REAL(8)   ,PARAMETER   :: RSMOOTH=R1PAR
      INTEGER(4)             :: GID   ! GRID ID
      INTEGER(4)             :: NR    ! #(RADIAL GRID POINTS)
      INTEGER(4)             :: LNX   ! #(PARTIAL WAVES INCLUDING SCATTERING )
      INTEGER(4)             :: NS    ! #(SINGLE FUNCTIONS)
      INTEGER(4)             :: NP    ! #(PRODUCT FUNCTIONS)
      INTEGER(4)             :: NT    ! #(TRIPLE PRODUCT FUNCTIONS)
      INTEGER(4)             :: NE    ! #(EXPONENTS)
      INTEGER(4)             :: NPOW  ! #(POWERS)
      INTEGER(4)             :: NPOW2 ! #(POWERS)
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)   ,ALLOCATABLE :: AUX1(:)
      REAL(8)   ,ALLOCATABLE :: AUX2(:)
      REAL(8)   ,ALLOCATABLE :: W(:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      INTEGER(4)             :: ISP,LN1,LN2,LN3,L1,L2,L3,LR1,LR2,IS,IP,IT,IE,IR,J
      INTEGER(4)             :: IRSMOOTH
      REAL(8)                :: SVAR,SVAR1,SVAR2,A,B
CHARACTER(128) :: STRING,STRING1,STRING2
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_TAILEDPRODUCTS')
      DO ISP=1,NSP
        LNX=POTPAR(ISP)%TAILED%LNX
        ALLOCATE(LOX(LNX))
        LOX=POTPAR(ISP)%TAILED%LOX
        GID=POTPAR(ISP)%TAILED%GID
        CALL RADIAL$GETI4(GID,'NR',NR)
!
!       ========================================================================
!       == COUNT NUMBER OF PRODUCTS                                           ==
!       ========================================================================
        NT=0   !#(TRIPLES)
        NP=0   !#(PRODUCTS)
        NS=0   !#(SINGLES)
        DO LN1=1,LNX
          L1=LOX(LN1)
          NPOW2=INT(0.5D0*REAL(NX-L1)+1.000001D0)  !R^(L+2N), N=0,NPOW-1
          IF(NPOW2.LT.1) CYCLE
          NS=NS+1
          DO LN2=LN1,LNX
            L2=LOX(LN2)
            DO LR1=ABS(L1-L2),L1+L2,2  ! TRIANGLE RULE
              NPOW2=INT(0.5D0*REAL(NX-LR1)+1.000001D0)  !R^(L+2N), N=0,NPOW2-1
              IF(NPOW2.LT.1) CYCLE
              NP=NP+1
              DO LN3=1,LNX
                L3=LOX(LN3)
                DO LR2=ABS(LR1-L3),LR1+L3,2 ! TRIANGLE RULE
                  NPOW2=INT(0.5D0*REAL(NX-LR2)+1.000001D0) !R^(L+2N), N=0,NPOW-1
                  IF(NPOW2.LT.1) CYCLE
                  NT=NT+1
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
!       ========================================================================
!       == DEFINE GAUSSIANS                                                   ==
!       ========================================================================
        NE=NEPAR        
        NPOW=NPOWPAR        
        POTPAR(ISP)%TAILED%PRODRHO%NIJK=NPOW
        POTPAR(ISP)%TAILED%PRODPOT%NIJK=NPOW
        POTPAR(ISP)%TAILED%SINGLE%NIJK =NPOW
        POTPAR(ISP)%TAILED%TRIPLE%NIJK =NPOW
        POTPAR(ISP)%TAILED%PRODRHO%NE  =NE
        POTPAR(ISP)%TAILED%PRODPOT%NE  =NE
        POTPAR(ISP)%TAILED%SINGLE%NE   =NE
        POTPAR(ISP)%TAILED%TRIPLE%NE   =NE
        POTPAR(ISP)%TAILED%PRODRHO%NORB=NP
        POTPAR(ISP)%TAILED%PRODPOT%NORB=NP
        POTPAR(ISP)%TAILED%SINGLE%NORB =NS
        POTPAR(ISP)%TAILED%TRIPLE%NORB =NT
        ALLOCATE(POTPAR(ISP)%TAILED%PRODRHO%E(NE))
        ALLOCATE(POTPAR(ISP)%TAILED%PRODPOT%E(NE))
        ALLOCATE(POTPAR(ISP)%TAILED%SINGLE%E(NE))
        ALLOCATE(POTPAR(ISP)%TAILED%TRIPLE%E(NE))
        ALLOCATE(POTPAR(ISP)%TAILED%PRODRHO%C(NPOW,NE,NP))
        ALLOCATE(POTPAR(ISP)%TAILED%PRODPOT%C(NPOW,NE,NP))
        ALLOCATE(POTPAR(ISP)%TAILED%SINGLE%C(NPOW,NE,NS))
        ALLOCATE(POTPAR(ISP)%TAILED%TRIPLE%C(NPOW,NE,NT))
        DO IE=1,NE
          POTPAR(ISP)%TAILED%PRODRHO%E(IE)=1.D0/(R1PAR*FACPAR**(IE-1))
        ENDDO
        POTPAR(ISP)%TAILED%PRODPOT%E(:)=POTPAR(ISP)%TAILED%PRODRHO%E(:)
        POTPAR(ISP)%TAILED%SINGLE%E(:) =POTPAR(ISP)%TAILED%PRODRHO%E(:)
        POTPAR(ISP)%TAILED%TRIPLE%E(:) =POTPAR(ISP)%TAILED%PRODRHO%E(:)
!
!       ========================================================================
!       == DO THE FIT OF THE PRODUCTS OF TAILED ORBITALS                      ==
!       ========================================================================
        ALLOCATE(AUX(NR))
        ALLOCATE(AUX1(NR))
        ALLOCATE(AUX2(NR))
        ALLOCATE(W(NR)) ! FITTING WEIGHT FUNCTION
        ALLOCATE(R(NR)) ! FITTING WEIGHT FUNCTION
        CALL RADIAL$R(GID,NR,R)
        DO IR=1,NR
          IF(R(IR).GT.RSMOOTH) EXIT 
          IRSMOOTH=IR
        ENDDO
!       == CONSTRUCT WEIGHT FUNCTION =======================================
!       == LEAVING TAILS THAT CANNOT BE FITTED LEADS TO OSZILLATIONS
        AUX(:)=MINVAL(POTPAR(ISP)%TAILED%PRODPOT%E(:))*R(:)**2
        W(:)=1.D0
        SVAR=1.D0
        DO J=1,NX/2+2  ! IT SEEMS TO BETTER TO GO TWO ORDERS HIGHER
          SVAR=SVAR/REAL(J,KIND=8)
          W(:)=W(:)+SVAR*AUX(:)**J
        ENDDO
        W(:)=W(:)*EXP(-AUX)
        W(:)=W(:)*R(:)**2
!       == WEIGHTFUNCTION DONE =========================
        POTPAR(ISP)%TAILED%PRODRHO%C=0.D0
        POTPAR(ISP)%TAILED%PRODPOT%C=0.D0
        POTPAR(ISP)%TAILED%SINGLE%C=0.D0
        POTPAR(ISP)%TAILED%TRIPLE%C=0.D0
        IS=0
        IP=0
        IT=0
        DO LN1=1,LNX
          L1=LOX(LN1)
          AUX=POTPAR(ISP)%TAILED%AEF(:,LN1)
          NPOW2=INT(0.5D0*REAL(NX-L1)+1.000001D0)  !R^(L+2N), N=0,NPOW-1
          IF(NPOW2.LT.1) CYCLE
          IS=IS+1
          CALL GAUSSIAN_FITGAUSS(GID,NR,W,L1,AUX,NE,NPOW2 &
       &                         ,POTPAR(ISP)%TAILED%SINGLE%E &
       &                         ,POTPAR(ISP)%TAILED%SINGLE%C(:NPOW2,:,IS))
!!$WRITE(STRING1,*)ISP
!!$WRITE(STRING2,*)IS
!!$STRING='FITTEST_S_'//TRIM(ADJUSTL(STRING1))//'_'//TRIM(ADJUSTL(STRING2))
!!$CALL LMTO_TESTPLOTRADIALGAUSS(STRING,L1,GID,NR,AUX,NPOW,NE &
!!$& ,POTPAR(ISP)%TAILED%SINGLE%E,POTPAR(ISP)%TAILED%SINGLE%C(:,:,IS))
!!$STRING='F'//TRIM(ADJUSTL(STRING))
!!$CALL LMTO_TESTPLOTRADIALGAUSSALONE(STRING,L1,NPOW,NE &
!!$& ,POTPAR(ISP)%TAILED%SINGLE%E,POTPAR(ISP)%TAILED%SINGLE%C(:,:,IS))
          DO LN2=LN1,LNX
            L2=LOX(LN2)
            AUX=POTPAR(ISP)%TAILED%AEF(:,LN1)*POTPAR(ISP)%TAILED%AEF(:,LN2)
            DO LR1=ABS(L1-L2),L1+L2,2  ! TRIANGLE RULE
              NPOW2=INT(0.5D0*REAL(NX-LR1)+1.000001D0)  !R^(L+2N), N=0,NPOW2-1
              IF(NPOW2.LT.1) CYCLE
              IP=IP+1
              AUX1=AUX
!!$!             == REPLACE BY A*R^L+BR^(L+2) WITH VALUE AND MULTIPOLE MOMENT ==
!!$              AUX1(:)=AUX(:)*R(:)**(LR1+2)
!!$              CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
!!$              CALL RADIAL$VALUE(GID,NR,AUX2,RSMOOTH,SVAR1)
!!$              SVAR1=SVAR1/RSMOOTH**(LR1+3)
!!$              CALL RADIAL$VALUE(GID,NR,AUX,RSMOOTH,SVAR2)
!!$              A= 0.5D0*REAL(2*LR1+3,KIND=8) &
!!$      &              *(REAL(2*LR1+5,KIND=8)*SVAR1-SVAR2)
!!$              B=-0.5D0*REAL(2*LR1+5,KIND=8) &
!!$      &              *(REAL(2*LR1+3,KIND=8)*SVAR1-SVAR2)
!!$              AUX1=AUX
!!$              AUX1(:IRSMOOTH)=A*(R(:IRSMOOTH)/RSMOOTH)**LR1 &
!!$      &                      +B*(R(:IRSMOOTH)/RSMOOTH)**(LR1+2) 
!!$!             == REPLACEMENT DONE============================================
              CALL GAUSSIAN_FITGAUSS(GID,NR,W,LR1,AUX1,NE,NPOW2 &
       &                         ,POTPAR(ISP)%TAILED%PRODRHO%E &
       &                         ,POTPAR(ISP)%TAILED%PRODRHO%C(:NPOW2,:,IP))
!!$WRITE(STRING1,*)ISP
!!$WRITE(STRING2,*)IP
!!$STRING='FITTEST_R_'//TRIM(ADJUSTL(STRING1))//'_'//TRIM(ADJUSTL(STRING2))
!!$CALL LMTO_TESTPLOTRADIALGAUSS(STRING,LR1,GID,NR,AUX1,NPOW,NE &
!!$& ,POTPAR(ISP)%TAILED%PRODRHO%E,POTPAR(ISP)%TAILED%PRODRHO%C(:,:,IP))
!!$STRING='F'//TRIM(ADJUSTL(STRING))
!!$CALL LMTO_TESTPLOTRADIALGAUSSALONE(STRING,LR1,NPOW,NE &
!!$& ,POTPAR(ISP)%TAILED%PRODRHO%E,POTPAR(ISP)%TAILED%PRODRHO%C(:,:,IP))
!             == CONSTRUCT ELECTROSTATIC POTENTIAL =============================
              CALL RADIAL$POISSON(GID,NR,LR1,AUX1,AUX2)
              AUX1=AUX2              
              AUX1(1)=AUX1(2)  ! AVOID POTENTIAL SINGULARITY AT THE ORIGIN
              NPOW2=INT(0.5D0*REAL(NX-LR1)+1.000001D0)  !R^(L+2N), N=0,NPOW-1
              CALL GAUSSIAN_FITGAUSS(GID,NR,W,LR1,AUX1,NE,NPOW2 &
       &                         ,POTPAR(ISP)%TAILED%PRODPOT%E &
       &                         ,POTPAR(ISP)%TAILED%PRODPOT%C(:NPOW2,:,IP))
!!$WRITE(STRING1,*)ISP
!!$WRITE(STRING2,*)IP
!!$STRING='FITTEST_P_'//TRIM(ADJUSTL(STRING1))//'_'//TRIM(ADJUSTL(STRING2))
!!$CALL LMTO_TESTPLOTRADIALGAUSS(STRING,LR1,GID,NR,AUX1,NPOW,NE &
!!$& ,POTPAR(ISP)%TAILED%PRODPOT%E,POTPAR(ISP)%TAILED%PRODPOT%C(:,:,IP))
!!$STRING='F'//TRIM(ADJUSTL(STRING))
!!$CALL LMTO_TESTPLOTRADIALGAUSSALONE(STRING,LR1,NPOW,NE &
!!$& ,POTPAR(ISP)%TAILED%PRODPOT%E,POTPAR(ISP)%TAILED%PRODPOT%C(:,:,IP))
              DO LN3=1,LNX
                L3=LOX(LN3)
                DO LR2=ABS(LR1-L3),LR1+L3,2 ! TRIANGLE RULE
                  NPOW2=INT(0.5D0*REAL(NX-LR2)+1.000001D0) !R^(L+2N), N=0,NPOW-1
                  IF(NPOW2.LT.1) CYCLE
                  IT=IT+1
                  AUX2(:)=AUX1(:)*POTPAR(ISP)%TAILED%AEF(:,LN3)
                  CALL GAUSSIAN_FITGAUSS(GID,NR,W,LR2,AUX2,NE,NPOW2 &
       &                         ,POTPAR(ISP)%TAILED%TRIPLE%E &
       &                         ,POTPAR(ISP)%TAILED%TRIPLE%C(:NPOW2,:,IT))
!!$WRITE(STRING1,*)ISP
!!$WRITE(STRING2,*)IT
!!$STRING='FITTEST_T_'//TRIM(ADJUSTL(STRING1))//'_'//TRIM(ADJUSTL(STRING2))
!!$CALL LMTO_TESTPLOTRADIALGAUSS(STRING,LR2,GID,NR,AUX2,NPOW,NE &
!!$& ,POTPAR(ISP)%TAILED%TRIPLE%E,POTPAR(ISP)%TAILED%TRIPLE%C(:,:,IT))
!!$STRING='F'//TRIM(ADJUSTL(STRING))
!!$CALL LMTO_TESTPLOTRADIALGAUSSALONE(STRING,LR2,NPOW,NE &
!!$& ,POTPAR(ISP)%TAILED%TRIPLE%E,POTPAR(ISP)%TAILED%TRIPLE%C(:,:,IT))
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(AUX)
        DEALLOCATE(AUX1)
        DEALLOCATE(AUX2)
        DEALLOCATE(W)
        DEALLOCATE(R)
        DEALLOCATE(LOX)
      ENDDO
!!$CALL LMTO_TESTTAILEDP(NX)
!!$STOP 'FORCED'
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$PROJTONTBO(XK,NDIM,NBH,NPRO,PROJ)
!     **************************************************************************
!     ** TRANSFORMS THE PROJECTIONS ONTO COEFFICIENTS FOR                     **
!     ** NATURAL TIGHT-BINDING ORBITALS                                       **
!     **************************************************************************
      USE LMTO_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: XK(3)
      INTEGER(4),INTENT(IN)   :: NDIM
      INTEGER(4),INTENT(IN)   :: NBH
      INTEGER(4),INTENT(IN)   :: NPRO
      COMPLEX(8),INTENT(INOUT):: PROJ(NDIM,NBH,NPRO)
      REAL(8)                 :: A(NPRO)
      REAL(8)                 :: B(NPRO)
      REAL(8)                 :: C(NRL)
      REAL(8)                 :: D(NPRO)
      REAL(8)                 :: E(NRL)
      REAL(8)                 :: F(NRL)
      COMPLEX(8)              :: G(NRL,NRL)
      COMPLEX(8)              :: H(NRL,NRL)
      COMPLEX(8),ALLOCATABLE  :: PROJCONTR1(:,:,:)
      COMPLEX(8),ALLOCATABLE  :: PROJCONTR2(:,:,:)
      COMPLEX(8),ALLOCATABLE  :: PROJCONTR3(:,:,:)
      INTEGER(4)              :: NAT
      INTEGER(4)              :: I,L,ISP,IPRO,IAT,LN,I1,IM,J
!     **************************************************************************
      IF(.NOT.TON) RETURN
                               CALL TRACE$PUSH('LMTO$PROJTONTBO')
      CALL LMTO_PREPARE1(NPRO,NRL,A,B,C,D,E,F)
      CALL LMTO_PREPARE2(XK,NRL,C,E,F,G,H)
!
!     ==========================================================================
!     == CONTRACT PROJECTIONS SO THAT FOR EACH ANGULAR MOMENTUM ONE TERM REMAINS
!     ==========================================================================
      DO I=1,NPRO
        PROJ(:,:,I)=PROJ(:,:,I)*A(I)
      ENDDO
!
!     ==========================================================================
!     == CONTRACT PROJECTIONS SO THAT FRO EACH ANGULAR MOMENTUM ONE TERM REMAINS
!     ==========================================================================
      ALLOCATE(PROJCONTR1(NDIM,NBH,NRL))      
      ALLOCATE(PROJCONTR2(NDIM,NBH,NRL))      
      PROJCONTR1(:,:,:)=(0.D0,0.D0)
      PROJCONTR2(:,:,:)=(0.D0,0.D0)
      NAT=SIZE(ISPECIES)
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I1=SBARATOMI1(IAT)-1+SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            IPRO=IPRO+1
            I1=I1+1
            PROJCONTR1(:,:,I1)=PROJCONTR1(:,:,I1)+PROJ(:,:,IPRO)
            PROJCONTR2(:,:,I1)=PROJCONTR2(:,:,I1)+B(IPRO)*PROJ(:,:,IPRO)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == MATRIX MULTIPLICATION                                                ==
!     ==========================================================================
      ALLOCATE(PROJCONTR3(NDIM,NBH,NRL))
      PROJCONTR3(:,:,:)=(0.D0,0.D0)
      DO I=1,NRL
        DO J=1,NRL
          PROJCONTR3(:,:,I)=PROJCONTR3(:,:,I)+H(I,J)*PROJCONTR1(:,:,J) &
     &                                       +G(I,J)*PROJCONTR2(:,:,J)
        ENDDO
      ENDDO
      PROJCONTR1(:,:,:)=PROJCONTR3(:,:,:)
      DEALLOCATE(PROJCONTR3)
      DEALLOCATE(PROJCONTR2)
!
!     ==========================================================================
!     == EXPAND AGAIN
!     ==========================================================================
!     == A*PROJ IS ALREADY MAPPED ONTO PROJ !
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I1=SBARATOMI1(IAT)-1+SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            IPRO=IPRO+1
            I1=I1+1
            PROJ(:,:,IPRO)=PROJ(:,:,IPRO)-D(IPRO)*PROJCONTR1(:,:,I1)
          ENDDO
        ENDDO
      ENDDO
                               CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$NTBOTOPROJ(XK,NDIM,NBH,NPRO,PROJ)
!     **************************************************************************
!     ** TRANSFORMS THE DERIVATIVE OF THE TOTAL ENERGY WITH RESPECT TO A      **
!     ** NTBO-COEFFICIENT INTO A PREFACTOR FOR THE PROJECTOR IN THE           **
!     ** PSEUDO HAMILTONIAN                                                   **
!     **************************************************************************
      USE LMTO_MODULE,ONLY : TON,NRL,ISPECIES,LNX,LOX,SBARATOMI1,SBARLI1,THTBC
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)   :: XK(3)
      INTEGER(4),INTENT(IN)   :: NDIM
      INTEGER(4),INTENT(IN)   :: NBH
      INTEGER(4),INTENT(IN)   :: NPRO
      COMPLEX(8),INTENT(INOUT):: PROJ(NDIM,NBH,NPRO)
      REAL(8)                 :: A(NPRO)
      REAL(8)                 :: B(NPRO)
      REAL(8)                 :: C(NRL)
      REAL(8)                 :: D(NPRO)
      REAL(8)                 :: E(NRL)
      REAL(8)                 :: F(NRL)
      COMPLEX(8)              :: G(NRL,NRL)
      COMPLEX(8)              :: H(NRL,NRL)
      COMPLEX(8),ALLOCATABLE  :: PROJCONTR1(:,:,:)
      COMPLEX(8),ALLOCATABLE  :: PROJCONTR2(:,:,:)
      COMPLEX(8),ALLOCATABLE  :: PROJCONTR3(:,:,:)
      INTEGER(4)              :: NAT
      INTEGER(4)              :: I,L,ISP,IPRO,IAT,LN,I1,IM,J
!     **************************************************************************
      IF(.NOT.TON) RETURN
      IF(.NOT.THTBC) THEN
         PROJ=(0.D0,0.D0)
         RETURN
      END IF
                               CALL TRACE$PUSH('LMTO$PROJTONTBO')
      CALL LMTO_PREPARE1(NPRO,NRL,A,B,C,D,E,F)
      CALL LMTO_PREPARE2(XK,NRL,C,E,F,G,H)
      G=CONJG(G)
      H=CONJG(H)
!
!     ==========================================================================
!     == D * DE/DQ
!     ==========================================================================
      ALLOCATE(PROJCONTR1(NDIM,NBH,NRL))      
      PROJCONTR1(:,:,:)=(0.D0,0.D0)
      NAT=SIZE(ISPECIES)
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I1=SBARATOMI1(IAT)-1+SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            IPRO=IPRO+1
            I1=I1+1
            PROJCONTR1(:,:,I1)=PROJCONTR1(:,:,I1)+D(IPRO)*PROJ(:,:,IPRO)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == MATRIX MULTIPLICATION                                                ==
!     ==========================================================================
      ALLOCATE(PROJCONTR2(NDIM,NBH,NRL))      
      ALLOCATE(PROJCONTR3(NDIM,NBH,NRL))
      PROJCONTR2(:,:,:)=(0.D0,0.D0)
      PROJCONTR3(:,:,:)=(0.D0,0.D0)
      DO I=1,NRL
        DO J=1,NRL
          PROJCONTR2(:,:,I)=PROJCONTR2(:,:,I)+G(J,I)*PROJCONTR1(:,:,J)
          PROJCONTR3(:,:,I)=PROJCONTR3(:,:,I)+H(J,I)*PROJCONTR1(:,:,J)
        ENDDO
      ENDDO
      DEALLOCATE(PROJCONTR1)
!
!     ==========================================================================
!     == CONTRACT PROJECTIONS SO THAT FRO EACH ANGULAR MOMENTUM ONE TERM REMAINS
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I1=SBARATOMI1(IAT)-1+SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            IPRO=IPRO+1
            I1=I1+1
            PROJ(:,:,IPRO)=PROJ(:,:,IPRO)-B(IPRO)*PROJCONTR2(:,:,I1) &
     &                                   -PROJCONTR3(:,:,I1)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(PROJCONTR3)
      DEALLOCATE(PROJCONTR2)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO I=1,NPRO
        PROJ(:,:,I)=PROJ(:,:,I)*A(I)
      ENDDO
                               CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PREPARE1(NPRO,NRL_,A,B,C,D,E,F)
!     **************************************************************************
!     **  PREPARE ARRAYS NEEDED FOR THE TRANSFORMATION                        **
!     **  TO NATURAL TIGHT-BINDING ORBITALS                                   **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : NRL,ISPECIES,LNX,LOX,SBARATOMI1,SBARLI1,POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NPRO
      INTEGER(4),INTENT(IN)  :: NRL_
      REAL(8)   ,INTENT(OUT) :: A(NPRO)
      REAL(8)   ,INTENT(OUT) :: B(NPRO)
      REAL(8)   ,INTENT(OUT) :: C(NRL_)
      REAL(8)   ,INTENT(OUT) :: D(NPRO)
      REAL(8)   ,INTENT(OUT) :: E(NRL_)
      REAL(8)   ,INTENT(OUT) :: F(NRL_)
      INTEGER(4)             :: NAT
      INTEGER(4)             :: ISP,IAT,L,LN,IM,IPRO,IRL
!     **************************************************************************
      IF(NRL_.NE.NRL) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$I4VAL('NRL',NRL)
        CALL ERROR$I4VAL('NRL_',NRL_)
        CALL ERROR$MSG('LMTO_PREPARE1')
      END IF
      NAT=SIZE(ISPECIES)
      
!
!     == INITIALIZE OUTPUT ARRAYS WITH ZEROS ===================================
      A(:)=0.D0
      B(:)=0.D0
      C(:)=0.D0
      D(:)=0.D0
      E(:)=1.D0
      F(:)=0.D0
!      
!     ==  CONSTRUCT ARRAYS =====================================================
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IRL=SBARATOMI1(IAT)-1+SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            IPRO=IPRO+1
            IRL=IRL+1
            A(IPRO)=1.D0/POTPAR(ISP)%KTOPHI(LN)
            B(IPRO)=POTPAR(ISP)%KTOPHIDOT(LN)
            C(IRL)=C(IRL)-POTPAR(ISP)%JBARTOPHIDOT(LN)
            D(IPRO)=A(IPRO)*POTPAR(ISP)%PHIDOTPROJ(LN)
            E(IRL)=E(IRL)+B(IPRO)*D(IPRO)
            F(IRL)=F(IRL)+D(IPRO)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PREPARE2(XK,NRL,C,E,F,G,H)
!     **************************************************************************
!     **  PREPARE MATRICES NEEDED FOR THE TRANSFORMATION                      **
!     **  TO NATURAL TIGHT-BINDING ORBITALS                                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NRL
      REAL(8)   ,INTENT(IN)  :: XK(3)  ! K-POINT IN RELATIVE COORDINATES
      REAL(8)   ,INTENT(IN)  :: C(NRL)
      REAL(8)   ,INTENT(IN)  :: E(NRL)
      REAL(8)   ,INTENT(IN)  :: F(NRL)
      COMPLEX(8),INTENT(OUT) :: G(NRL,NRL)
      COMPLEX(8),INTENT(OUT) :: H(NRL,NRL)
      COMPLEX(8),ALLOCATABLE :: SBAR(:,:) ! SCREENED STRUCTURE CONSTANTS
      INTEGER(4)             :: I                  
!     **************************************************************************
      ALLOCATE(SBAR(NRL,NRL))
      CALL LMTO_SOFK(XK,NRL,SBAR)
      DO I=1,NRL
        G(:,I)=C(:)*CONJG(SBAR(I,:))*F(I)   !C
        G(I,I)=G(I,I)+E(I)
      ENDDO
      CALL LIB$INVERTC8(NRL,G,H)
      G(:,:)=H(:,:)
      DO I=1,NRL
        H(:,I)=H(:,I)*C(I)
      ENDDO
      H=MATMUL(H,TRANSPOSE(CONJG(SBAR)))    !C
      DEALLOCATE(SBAR)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SOFK(XK,N,SOFK)
!     **************************************************************************
!     **  TRANSFORMS THE SCREENED STRUCTURE CONSTANTS INTO K-SPACE            **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : SBARATOMI1,SBARATOMI2,SBAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N
      REAL(8)   ,INTENT(IN)  :: XK(3)  ! K-POINT IN FRACTIONAL COORDINATES
      COMPLEX(8),INTENT(OUT) :: SOFK(N,N)
      REAL(8)                :: KR
      COMPLEX(8)             :: EIKR
      INTEGER(4)             :: IAT1,IAT2
      INTEGER(4)             :: I1OFAT1,I2OFAT1        
      INTEGER(4)             :: I1OFAT2,I2OFAT2        
      INTEGER(4)             :: NN
      INTEGER(4)             :: NNS
      REAL(8)                :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      IF(N.NE.MAXVAL(SBARATOMI2)) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$I4VAL('MAXVAL(SBARATOMI2)',MAXVAL(SBARATOMI2))
        CALL ERROR$MSG('LMTO_SOFK')
      END IF
      SOFK(:,:)=(0.D0,0.D0)
      NNS=SIZE(SBAR)
      DO NN=1,NNS 
        IAT1=SBAR(NN)%IAT1
        IAT2=SBAR(NN)%IAT2
        KR=2.D0*PI*SUM(REAL(SBAR(NN)%IT(:))*XK(:))
        I1OFAT1=SBARATOMI1(IAT1)
        I2OFAT1=SBARATOMI2(IAT1)
        I1OFAT2=SBARATOMI1(IAT2)
        I2OFAT2=SBARATOMI2(IAT2)
        EIKR=CMPLX(COS(KR),SIN(KR)) 
        SOFK(I1OFAT1:I2OFAT1,I1OFAT2:I2OFAT2) &
     &                  =SOFK(I1OFAT1:I2OFAT1,I1OFAT2:I2OFAT2)+SBAR(NN)%MAT*EIKR
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PRPROJ(ID)
!     **************************************************************************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,THIS,WAVES_SELECTWV
      USE LMTO_MODULE, ONLY : LOX,LNX,ISPECIES
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)             :: NAT
      INTEGER(4),ALLOCATABLE :: IPRO1(:)
      INTEGER(4),ALLOCATABLE :: NPROAT(:)
      INTEGER(4)             :: IAT,ISP,IKPT,ISPIN,IB,IPRO,NB,NBH,I0
!     **************************************************************************
      NAT=SIZE(ISPECIES)
      ALLOCATE(IPRO1(NAT))
      ALLOCATE(NPROAT(NAT))
      IPRO=1
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        IPRO1(IAT)=IPRO
        NPROAT(IAT)=SUM(2*LOX(:LNX(ISP),ISP)+1)
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          NBH=THIS%NBH
          NB=THIS%NB
          PRINT*,"==============IKPT=",IKPT,' IAT=',IAT,' ID=',ID
          DO IAT=1,NAT
            I0=IPRO1(IAT)-1
            DO IB=1,NBH
              WRITE(*,FMT='(I3,40("(",2F10.5,")"))')IB,THIS%PROJ(:,IB,I0+1:I0+NPROAT(IAT))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBODENMAT()
!     **************************************************************************
!     **  CONSTRUCT DENSITY MATRIX IN A NTBO BASIS                            **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : DENMAT,LOX,LNX,ISPECIES
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      INTEGER(4)             :: N1,N2
      INTEGER(4)             :: NNS
      INTEGER(4)             :: NND
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: NB,NBH,NBX
      INTEGER(4)             :: NDIMD
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      INTEGER(4),ALLOCATABLE :: IPRO1(:)
      INTEGER(4),ALLOCATABLE :: NPROAT(:)
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      INTEGER(4)             :: IAT,NN,II,ISP,IPRO,IKPT,ISPIN,I,J,IBH,IB
      REAL(8)                :: SVAR
      REAL(8)                :: F1,F2
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IAT1,IAT2,IT(3),I0,J0,IDIM,JDIM
      COMPLEX(8)             :: EIKR,C1(NDIM),C2(NDIM),CSVAR22(NDIM,NDIM)
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      REAL(8)                :: PI
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
COMPLEX(8)  :: PHASE
      INTEGER(4)             :: NTASKS,THISTASK,ICOUNT
!     **************************************************************************
                                           CALL TRACE$PUSH('LMTO_NTBODENMAT')
      PI=4.D0*ATAN(1.D0)
      IF(.NOT.ASSOCIATED(THIS%TBC)) THEN
        CALL ERROR$MSG('THIS%TBC NOT ASSOCIATED')
        CALL ERROR$STOP('LMTO_NTBODENMAT')
      END IF
!
!     ==========================================================================
!     == ALLOCATE DENSITY MATRIX
!     ==========================================================================
      CALL LMTO_DENMATLAYOUT()
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
!
!     ==========================================================================
!     ==  CONSTRUCT INDEX ARRAYS                                              ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(IPRO1(NAT))
      ALLOCATE(NPROAT(NAT))
      IPRO=1
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        IPRO1(IAT)=IPRO
        NPROAT(IAT)=SUM(2*LOX(:LNX(ISP),ISP)+1)
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
!
!     ==========================================================================
!     ==  ADD UP DENSITY MATRIX                                               ==
!     ==========================================================================
      NND=SIZE(DENMAT)
      NPRO=MAP%NPRO
      CALL MPE$QUERY('K',NTASKS,THISTASK)
      ICOUNT=0
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          NBH=THIS%NBH
          NB=THIS%NB
!!$PRINT*,"============== IKPT=",IKPT," XK=",XK(:,IKPT),I0,J0
!!$DO IAT=1,NAT
!!$  I0=IPRO1(IAT)
!!$  J0=I0-1+NPROAT(IAT)
!!$  IF(TINV) THEN
!!$    DO IBH=1,NBH
!!$      IF(OCC(2*IBH-1,IKPT,ISPIN).LT.1.D-5) CYCLE
!!$      WRITE(*,FMT='("IAT="I3," IB=",I3,40("(",2F10.5,")"))') &
!!$&                    IAT,2*IBH-1,CMPLX(REAL(THIS%TBC(:,IBH,I0:J0)))
!!$      IF(OCC(2*IBH,IKPT,ISPIN).LT.1.D-5) CYCLE
!!$      WRITE(*,FMT='("IAT="I3," IB=",I3,40("(",2F10.5,")"))') &
!!$&                    IAT,2*IBH,CMPLX(AIMAG(THIS%TBC(:,IBH,I0:J0)))
!!$    ENDDO
!!$  ELSE
!!$    DO IB=1,NB
!!$      IF(OCC(IB,IKPT,ISPIN).LT.1.D-5) CYCLE
!!$      PHASE=THIS%TBC(1,IB,1)
!!$      PHASE=CONJG(PHASE)/SQRT((PHASE*CONJG(PHASE)))
!!$      WRITE(*,FMT='("IAT="I3," IB=",I3,40("(",2F10.5,")"))') &
!!$&                    IAT,IB,THIS%TBC(:,IB,I0:J0)*PHASE
!!$    ENDDO
!!$  END IF
!!$ENDDO
          DO NN=1,NND
            ICOUNT=ICOUNT+1
            IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) CYCLE
            IAT1=DENMAT(NN)%IAT1
            IAT2=DENMAT(NN)%IAT2
            IT=DENMAT(NN)%IT
!            SVAR=-2.D0*PI*SUM(XK(:,IKPT)*REAL(IT,KIND=8))
            SVAR=2.D0*PI*SUM(XK(:,IKPT)*REAL(IT,KIND=8))
            EIKR=EXP(CI*SVAR)  !<P_{R+T}|PSI>=<P_R|PSI>*EIKR
            I0=IPRO1(IAT1)-1
            J0=IPRO1(IAT2)-1

            DO I=1,NPROAT(IAT1)
              DO J=1,NPROAT(IAT2)
                IF(TINV) THEN
                  CSVAR22=(0.D0,0.D0)
                  DO IBH=1,NBH
                    F1=OCC(2*IBH-1,IKPT,ISPIN)
                    F2=OCC(2*IBH,IKPT,ISPIN)
                    C1(:)=THIS%TBC(:,IBH,I0+I)
                    C2(:)=THIS%TBC(:,IBH,J0+J)*EIKR   ! EXP(-I*K*T)
                    DO JDIM=1,NDIM
                      CSVAR22(:,JDIM)=CSVAR22(:,JDIM) &
        &                            +0.5D0*((F1+F2)*C1(:)*CONJG(C2(JDIM)) &
        &                                   +(F1-F2)*C1(:)*C2(JDIM))
                    ENDDO
                  ENDDO
                  CSVAR22=REAL(CSVAR22) ! IMAG(CSVAR) CONTAINS CRAP DUE TO SUPER WAVE FUNCTIONS
                ELSE
                  CSVAR22=(0.D0,0.D0)
                  DO IB=1,NB
                    F1=OCC(IB,IKPT,ISPIN)
                    C1(:)=THIS%TBC(:,IB,I0+I)
                    C2(:)=THIS%TBC(:,IB,J0+J)*EIKR ! EXP(-I*K*T)
                    DO JDIM=1,NDIM
                      CSVAR22(:,JDIM)=CSVAR22(:,JDIM)+F1*C1(:)*CONJG(C2(JDIM))
                    ENDDO
                  ENDDO
                END IF
!
!           == DISTRIBUTE ONTO DENSITY MATRIX ENTRIES ==========================
!           == D(IDIMD)=SUM_{IDIM,JDIM} D(IDIM,JDIM)*PAULI_{IDIMD}(JDIM,IDIM) ==
!           == IDIMD IN {TOTAL,X,Y,Z}; IDIM IN {UP,DOWN} =======================
!           == TRANSFORMATION MUST BE CONSISTENT WITH WAVES_DENMAT =============
                IF(NSPIN.EQ.1) THEN
                  IF(NDIM.EQ.1) THEN !NON-SPIN-POLARIZED
                    DENMAT(NN)%MAT(I,J,1)=DENMAT(NN)%MAT(I,J,1) &
          &                              +REAL(CSVAR22(1,1))
                  ELSE ! NONCOLLINEAR
                    DENMAT(NN)%MAT(I,J,1)=DENMAT(NN)%MAT(I,J,1) &
          &                              +REAL(CSVAR22(1,1)+CSVAR22(2,2))
                    DENMAT(NN)%MAT(I,J,2)=DENMAT(NN)%MAT(I,J,2) &
          &                              +REAL(CSVAR22(1,2)+CSVAR22(2,1))
                    DENMAT(NN)%MAT(I,J,3)=DENMAT(NN)%MAT(I,J,3) &
          &                              -AIMAG(CSVAR22(1,2)-CSVAR22(2,1))
                    DENMAT(NN)%MAT(I,J,4)=DENMAT(NN)%MAT(I,J,4) &
          &                              +REAL(CSVAR22(1,1)-CSVAR22(2,2))
                  END IF
                ELSE IF(NSPIN.EQ.2) THEN !SPIN POLARIZED
                  IF(ISPIN.EQ.1) THEN
                    DENMAT(NN)%MAT(I,J,1)=DENMAT(NN)%MAT(I,J,1)+REAL(CSVAR22(1,1))
                    DENMAT(NN)%MAT(I,J,4)=DENMAT(NN)%MAT(I,J,4)+REAL(CSVAR22(1,1))
                  ELSE
                    DENMAT(NN)%MAT(I,J,1)=DENMAT(NN)%MAT(I,J,1)+REAL(CSVAR22(1,1))
                    DENMAT(NN)%MAT(I,J,4)=DENMAT(NN)%MAT(I,J,4)-REAL(CSVAR22(1,1))
                  END IF
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  SUM OVER MONOMER INCLUDES ALSO THE KPOINT SUM                       ==
!     ==========================================================================
      DO NN=1,NND
        CALL MPE$COMBINE('MONOMER','+',DENMAT(NN)%MAT)
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,FMT='(82("="),T10," NON-LOCAL DENSITY MATRIX IN A NTBO BASIS ")')
        DO NN=1,NND
          IAT1=DENMAT(NN)%IAT1
          IAT2=DENMAT(NN)%IAT2
          IT=DENMAT(NN)%IT
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
       &                                                             IAT1,IAT2,IT
          N1=DENMAT(NN)%N1
          N2=DENMAT(NN)%N2
          DO I=1,1 !DENMAT(NN)%N3
            DO J=1,N1 
              WRITE(*,FMT='(I3,300F10.3)')I,DENMAT(NN)%MAT(J,:,I)
            ENDDO
            WRITE(*,FMT='(82("-"))')
          ENDDO
        ENDDO
        CALL ERROR$MSG('FORCED STOP AFTER PRINTING DENSITY MATRIX')
        CALL ERROR$STOP('LMTO_NTBODENMAT')
      END IF
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBODENMATDER()
!     **************************************************************************
!     **  CONSTRUCT HTBC FROM THE DERIVATIVE OF THE ENERGY WITH RESPECT TO    **
!     **  THE DENSITY MATRIX                                                  **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE WAVES_MODULE, ONLY: WVSET_TYPE,NKPTL,NSPIN,NDIM,THIS,MAP &
     &                       ,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : SBAR,HAMIL,LOX,LNX,ISPECIES,THTBC
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      INTEGER(4)             :: N1,N2
      INTEGER(4)             :: NND
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: NB,NBH
      INTEGER(4)             :: NDIMD
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      INTEGER(4),ALLOCATABLE :: IPRO1(:)
      INTEGER(4),ALLOCATABLE :: NPROAT(:)
      INTEGER(4)             :: IAT,NN,ISP,IPRO,IKPT,ISPIN,I,J,IBH,IB
      REAL(8)                :: SVAR
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IAT1,IAT2,IT(3),I0,J0,IDIM,JDIM
      COMPLEX(8)             :: EIKR,C1,C2,CSVAR22(NDIM,NDIM)
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      REAL(8)                :: PI
!     **************************************************************************
                                 CALL TRACE$PUSH('LMTO_NTBODENMATDER')
      PI=4.D0*ATAN(1.D0)
      THTBC=.TRUE.
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
!
!     ==========================================================================
!     ==  CONSTRUCT INDEX ARRAYS                                              ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(IPRO1(NAT))
      ALLOCATE(NPROAT(NAT))
      IPRO=1
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        IPRO1(IAT)=IPRO
        NPROAT(IAT)=SUM(2*LOX(:LNX(ISP),ISP)+1)
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
      NPRO=SUM(NPROAT(:))
!
!     ==========================================================================
!     ==  ADD UP DENSITY MATRIX                                               ==
!     ==========================================================================
      NND=SIZE(HAMIL)
      NPRO=MAP%NPRO
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          NBH=THIS%NBH
          NB=THIS%NB
          CALL PLANEWAVE$SELECT(GSET%ID)
          IF(.NOT.ASSOCIATED(THIS%HTBC))ALLOCATE(THIS%HTBC(NDIM,NBH,NPRO))
          THIS%HTBC(:,:,:)=(0.D0,0.D0)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          DO NN=1,NND
            IAT1=HAMIL(NN)%IAT1
            IAT2=HAMIL(NN)%IAT2
            IT=HAMIL(NN)%IT
            SVAR=2.D0*PI*SUM(XK(:,IKPT)*REAL(IT))
            EIKR=EXP(CI*SVAR)  !<P_{R+T}|PSI>=<P_R|PSI>*EIKR
            I0=IPRO1(IAT1)-1
            J0=IPRO1(IAT2)-1
            DO I=1,NPROAT(IAT1)
              DO J=1,NPROAT(IAT2)
!
!               == CONVERT FROM TOTAL/SPIN INTO UP-DOWN REPRESENTATION =====
                IF(NSPIN.EQ.1) THEN
                  IF(NDIM.EQ.1) THEN !NON-SPIN-POLARIZED
                    CSVAR22(1,1)=CMPLX(HAMIL(NN)%MAT(I,J,1),0.D0)
                  ELSE ! NONCOLLINEAR
                    CSVAR22(1,1)=CMPLX(HAMIL(NN)%MAT(I,J,1) &
          &                           +HAMIL(NN)%MAT(I,J,4),0.D0)
                    CSVAR22(1,2)=CMPLX(HAMIL(NN)%MAT(I,J,2) &
          &                          ,-HAMIL(NN)%MAT(I,J,3))
                    CSVAR22(2,1)=CMPLX(HAMIL(NN)%MAT(I,J,2) &
          &                          ,+HAMIL(NN)%MAT(I,J,3))
                    CSVAR22(2,2)=CMPLX(HAMIL(NN)%MAT(I,J,1) &
          &                           -HAMIL(NN)%MAT(I,J,4),0.D0)
                  END IF
                ELSE IF(NSPIN.EQ.2) THEN !SPIN POLARIZED
                  IF(ISPIN.EQ.1) THEN
                    CSVAR22(1,1)=CMPLX(HAMIL(NN)%MAT(I,J,1) &
          &                           +HAMIL(NN)%MAT(I,J,4),0.D0)
                  ELSE
                    CSVAR22(1,1)=CMPLX(HAMIL(NN)%MAT(I,J,1) &
          &                           -HAMIL(NN)%MAT(I,J,4),0.D0)
                  END IF
                END IF
                CSVAR22(:,:)=CSVAR22(:,:)*EIKR
!
                DO IB=1,NBH
                  DO IDIM=1,NDIM
                    DO JDIM=1,NDIM
! THIS IS THE OLD VERSION (THIS IS CORRECT: SEE METHODS SECTION 'SECOND QUANT..'
!!$                      THIS%HTBC(JDIM,IB,J0+J)=THIS%HTBC(JDIM,IB,J0+J) &
!!$      &                              +CSVAR22(JDIM,IDIM)*THIS%TBC(IDIM,IB,I0+I)
! THIS SHOULD BE CORRECT.(NO!)
                      THIS%HTBC(IDIM,IB,I0+I)=THIS%HTBC(IDIM,IB,I0+I) &
      &                           +CSVAR22(IDIM,JDIM)*THIS%TBC(JDIM,IB,J0+J)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SETHTBCTOZERO()
!     **************************************************************************
!     **************************************************************************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,THIS,WAVES_SELECTWV
      IMPLICIT NONE
      INTEGER(4) :: IKPT,ISPIN
!     **************************************************************************
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          IF(ASSOCIATED(THIS%HTBC))THIS%HTBC(:,:,:)=(0.D0,0.D0)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DENMATLAYOUT()
!     **************************************************************************
!     **  ALLOCATES THE DENSITY MATRIX                                        **
!     **  - ALL BONDS CONSIDERED IN SBAR ARE CONSIDERED HERE AS WELL          **
!     **  - THE SELECTION THALFBONDS SELECTS ONLY ONE OF TWO BONDS THAT       **
!     **    DIFFER ONLY IN THE DIRECTION                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : SBAR,DENMAT,ISPECIES,LNX,LOX
      IMPLICIT NONE
      INTEGER(4)      :: NNS
      INTEGER(4)      :: NND
      INTEGER(4)      :: NN
      INTEGER(4)      :: IAT1,IAT2,IT(3)
      INTEGER(4)      :: ISP,II
      INTEGER(4)      :: N1,N2
      LOGICAL(4),PARAMETER ::THALFBONDS=.FALSE.
!     **************************************************************************
                                            CALL TRACE$PUSH('LMTO_DENMATLAYOUT')
!
!     ==========================================================================
!     == COUNT NUMBER OF DENSITY-MATRIX ELEMENTS                              ==
!     ==========================================================================
      NNS=SIZE(SBAR)
      NND=0
      DO NN=1,NNS
        IAT1=SBAR(NN)%IAT1
        IAT2=SBAR(NN)%IAT2
        IT(:)=SBAR(NN)%IT(:)
        IF(THALFBONDS) THEN
          IF(IAT2.GT.IAT1) CYCLE
          IF(IAT2.EQ.IAT1) THEN
            IF(IT(1).GT.0) CYCLE
            IF(IT(1).EQ.0) THEN
              IF(IT(2).GT.0) CYCLE
              IF(IT(2).EQ.0) THEN
                IF(IT(3).GT.0) CYCLE
              END IF
            END IF
          END IF
        END IF
        NND=NND+1
      ENDDO
!
!     ==========================================================================
!     == ALLOCATE DENSITY MATRIX                                              ==
!     ==========================================================================
      IF(ALLOCATED(DENMAT)) THEN
        CALL ERROR$MSG('DENMAT MUST NOT BE ALLOCATED')
        CALL ERROR$STOP('LMTO_NTBODENMAT')
      END IF
      ALLOCATE(DENMAT(NND))
!
!     ==========================================================================
!     == INITIALIZE                                                           ==
!     ==========================================================================
      NN=0
      DO II=1,NNS
        IAT1=SBAR(II)%IAT1
        IAT2=SBAR(II)%IAT2
        IT(:)=SBAR(II)%IT(:)
        IF(THALFBONDS) THEN
          IF(IAT2.GT.IAT1) CYCLE
          IF(IAT2.EQ.IAT1) THEN
            IF(IT(1).GT.0) CYCLE
            IF(IT(1).EQ.0) THEN
              IF(IT(2).GT.0) CYCLE
              IF(IT(2).EQ.0) THEN
                IF(IT(3).GT.0) CYCLE
              END IF
            END IF
          END IF
        END IF
        NN=NN+1
        ISP=ISPECIES(IAT1)
        N1=SUM(2*LOX(:LNX(ISP),ISP)+1)
        ISP=ISPECIES(IAT2)
        N2=SUM(2*LOX(:LNX(ISP),ISP)+1)
        DENMAT(NN)%IAT1=IAT1
        DENMAT(NN)%IAT2=IAT2
        DENMAT(NN)%IT=IT
        DENMAT(NN)%N1=N1
        DENMAT(NN)%N2=N2
        DENMAT(NN)%N3=4  !(TOTAL,X,Y,Z)
        ALLOCATE(DENMAT(NN)%MAT(N1,N2,4))
        DENMAT(NN)%MAT(:,:,:)=0.D0
      ENDDO

                                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_CLEANDENMAT()
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : DENMAT,HAMIL
      IMPLICIT NONE
      INTEGER(4)  ::NN,NND
!     **************************************************************************
      NND=SIZE(DENMAT)
      DO NN=1,NND
        DEALLOCATE(DENMAT(NN)%MAT)
        DEALLOCATE(HAMIL(NN)%MAT)
      ENDDO
      DEALLOCATE(DENMAT)
      DEALLOCATE(HAMIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_ONSORTHO(IAT,LMNX,T,UNT)
!     **************************************************************************
!     **  CONSTRUCT TRANSFORMATION ONTO ONSITE ORTHONORMAL ORBITALS           **
!     **                                                                      **
!     **   |CHIPRIME_I>=SUM_J |CHI_J>T_JI;  <CHIPRIME_I|CHIPRIME_J>=DELTA_IJ  **
!     **   T*UNT=1                                                            **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR,LNX,LOX
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(OUT):: T(LMNX,LMNX)
      REAL(8)   ,INTENT(OUT):: UNT(LMNX,LMNX)
      INTEGER(4)            :: ISP
      INTEGER(4)            :: LMNXT
      REAL(8)               :: OVERLAP(LMNX,LMNX)
      REAL(8)               :: OVERLAP1(LMNX,LMNX)
      REAL(8)               :: SVAR
      INTEGER(4)            :: I,J,LMN,LN,L
      LOGICAL(4)            :: TORB(LMNX)
      LOGICAL(4),PARAMETER  :: TWRITE=.FALSE.
!     **************************************************************************
      ISP=ISPECIES(IAT)
!
!     ==========================================================================
!     == CONSTRUCT OVERLAP MATRIX                                             ==
!     == USE OVERLAP BETWEEN PUR ANGULAR MOMENTUM TAILED PARTIAL WAVES        ==
!     ==========================================================================
      LMNXT=POTPAR(ISP)%TAILED%LMNX
      CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,1 &
   &               ,LMNXT,LMNXT,POTPAR(ISP)%TAILED%OVERLAP,LMNX,LMNX,OVERLAP)
      IF(TWRITE) THEN
        PRINT*,'LMNXT ',LMNXT
        PRINT*,'LMNX  ',LMNX
        PRINT*,'TALLORB',POTPAR(ISP)%TALLORB
        PRINT*,'TORB   ',POTPAR(ISP)%TORB
        DO I=1,LMNXT
          WRITE(*,FMT='(I3," OBIG=",100E10.2)')I,POTPAR(ISP)%TAILED%OVERLAP(I,:)
        ENDDO
        DO I=1,LMNX
          WRITE(*,FMT='(I3," O=",100E10.2)')I,OVERLAP(I,:)
        ENDDO
        CALL ERROR$MSG('FORCED STOP')
        CALL ERROR$STOP('LMTO_ONSORTHO')
      END IF
!
!     ==========================================================================
!     == DETERIMINE ACTIVE ORBITALS                                           ==
!     ==========================================================================
      LMN=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        TORB(LMN+1:LMN+2*L+1)=POTPAR(ISP)%TORB(LN)
        LMN=LMN+2*L+1
      ENDDO
!
!     ==========================================================================
!     == ORTHONORMALIZE                                                       ==
!     ==========================================================================
      T(:,:)=0.D0
      DO I=1,LMNX
        T(I,I)=1.D0
      ENDDO
      DO I=1,LMNX
        IF(.NOT.TORB(I)) CYCLE
!       == ORTHOGONALIZE
        DO J=1,I-1
          IF(.NOT.(TORB(J).OR.POTPAR(ISP)%TALLORB)) CYCLE
          SVAR=OVERLAP(J,I)
          T(:,I)=T(:,I)-T(:,J)*SVAR
          OVERLAP1=OVERLAP
          OVERLAP1(:,I)=OVERLAP1(:,I)-OVERLAP(:,J)*SVAR
          OVERLAP1(I,:)=OVERLAP1(I,:)-OVERLAP(J,:)*SVAR
          OVERLAP(:,:)=OVERLAP1(:,:)
        ENDDO
!       == NORMALIZE
        SVAR=1.D0/SQRT(OVERLAP(I,I))
        T(:,I)=T(:,I)*SVAR
        OVERLAP(:,I)=OVERLAP(:,I)*SVAR
        OVERLAP(I,:)=OVERLAP(I,:)*SVAR
      ENDDO
!
      IF(TWRITE) THEN
        PRINT*,' IN ONSORTHO IAT=',IAT
        PRINT*,' TORB ',TORB
        DO I=1,LMNX
          WRITE(*,FMT='(I3," T=",100E10.2)')I,T(I,:)
        ENDDO
        DO I=1,LMNX
          WRITE(*,FMT='(I3," O=",100E10.2)')I,OVERLAP(I,:)
        ENDDO
        CALL ERROR$MSG('FORCED STOP')
        CALL ERROR$STOP('LMTO_ONSORTHO')
     END IF
!
!     ==========================================================================
!     == INVERT                                                               ==
!     ==========================================================================
      CALL LIB$INVERTR8(LMNX,T,UNT)
      DO I=1,LMNX
        IF(.NOT.(TORB(I).OR.POTPAR(ISP)%TALLORB)) THEN
          UNT(:,I)=0.D0
          T(:,I)=0.D0
          UNT(I,:)=0.D0
          T(I,:)=0.D0
        END IF
      END DO
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE LMTO_DROPPICK_MODULE
!===============================================================================
!== HARD-WIRED INPUT DATA                                                     ==
!===============================================================================
LOGICAL(4)              :: TREADDHOFK=.FALSE.
CHARACTER(32),PARAMETER :: SWITCHID='SRVO3'
!!CHARACTER(32),PARAMETER :: SWITCHID='CAFE2AS2'
!CHARACTER(32),PARAMETER :: SWITCHID='H2'
!CHARACTER(32),PARAMETER :: SWITCHID='HUBBARD'
INTEGER(4)           :: NB1        ! FIRST BAND IN W
LOGICAL(4),POINTER   :: TPRO(:)    ! SELECTOR FOR CORRELATED ORBITALS
!===============================================================================
!== DERIVED DATA                                                              ==
!===============================================================================
TYPE T_TYPE
  REAL(8),POINTER :: MAT(:,:)
  REAL(8),POINTER :: INV(:,:)
  INTEGER(4)      :: I1
  INTEGER(4)      :: I2
END TYPE T_TYPE
LOGICAL(4)           :: TINI_DROPPICK=.FALSE.
LOGICAL(4)           :: TPICKED=.FALSE.  ! READ ONLY ONCE AND KEEP
INTEGER(4)           :: NB2              ! LAST BAND IN W
INTEGER(4)           :: NBW              ! #(BANDS IN W)
INTEGER(4)           :: NCORR            ! #(CORRELATED ORBITALS)
TYPE(T_TYPE),POINTER :: T(:)             ! TRANSFORMATION TO ONSITE ORTHONORMAL
COMPLEX(8),POINTER   :: DHOFK(:,:,:,:)   ! HAMILTONIAN CORRECTION
INTEGER(4)           :: TICKET(8)=0      ! UNIQUE FILE-SET IDENTIFIER
INTEGER(4),POINTER   :: IPRO1(:)          
INTEGER(4),POINTER   :: NPROAT(:)
INTEGER(4),PARAMETER :: METHOD=2         ! LEGACY. IS TO BE REMOVED
END MODULE LMTO_DROPPICK_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPICK_INI()
!     **************************************************************************
!     **************************************************************************
      USE LMTO_DROPPICK_MODULE, ONLY: TINI_DROPPICK,TPRO,NCORR,NB1,NB2,NBW &
     &                               ,IPRO1,NPROAT,SWITCHID
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,THIS,MAP,WAVES_SELECTWV
      USE LMTO_MODULE, ONLY : TDROP,LOX,LNX,ISPECIES,POTPAR
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4)  :: IPRO,IAT,ISP,LN,L,M,IKPT,ISPIN
      INTEGER(4)  :: NPRO,NAT,NB
!     **************************************************************************
      IF(TINI_DROPPICK) RETURN
      TINI_DROPPICK=.TRUE.
!
!     ==========================================================================
!     ==  SELECT CORRELATED ORBITALS                                          ==
!     ==========================================================================
      IF(SWITCHID.EQ.'CAFE2AS2') THEN
        NB1=7
!        NB2=24
      ELSE IF(SWITCHID.EQ.'SRVO3') THEN
        NB1=17
!        NB2=19
      ELSE IF(SWITCHID.EQ.'H2') THEN
        NB1=1
      ELSE IF(SWITCHID.EQ.'HUBBARD') THEN
        NB1=1
      ELSE 
        CALL ERROR$MSG('SWITCHID NOT RECOGNIZED')
        CALL ERROR$CHVAL('SWITCHID',SWITCHID)
        CALL ERROR$STOP('LMTO_DROPICK_INI')
      END IF
!
!     ==========================================================================
!     ==  SELECT CORRELATED ORBITALS                                          ==
!     ==========================================================================
      NPRO=MAP%NPRO
      NAT=SIZE(ISPECIES)
      ALLOCATE(TPRO(NPRO))
      ALLOCATE(IPRO1(NAT))
      ALLOCATE(NPROAT(NAT))
      IPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        IPRO1(IAT)=IPRO+1
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          DO M=1,2*L+1
            IPRO=IPRO+1
            TPRO(IPRO)=POTPAR(ISP)%TORB(LN)
IF(TPRO(IPRO)) THEN
  IF(SWITCHID.EQ.'SRVO3') THEN
    IF(L.EQ.2) THEN
      IF(M.EQ.1) TPRO(IPRO)=.FALSE.   !SET EG ORBITALS OFF
      IF(M.EQ.3) TPRO(IPRO)=.FALSE.
    END IF
  ELSE IF(SWITCHID.EQ.'CAFE2AS2') THEN
  ELSE IF(SWITCHID.EQ.'H2') THEN
  ELSE IF(SWITCHID.EQ.'HUBBARD') THEN
  ELSE 
    CALL ERROR$MSG('SWITCHID NOT RECOGNIZED (2ND TEST)')
    CALL ERROR$CHVAL('SWITCHID',SWITCHID)
    CALL ERROR$STOP('LMTO_DROPICK_INI')
  END IF
END IF
          ENDDO
        ENDDO
        NPROAT(IAT)=IPRO-IPRO1(IAT)+1
      ENDDO 
!
!     ==========================================================================
!     == COUNT NUMBER OF CORRELATED ORBITALS                                  ==
!     ==========================================================================
      NCORR=0
      DO IPRO=1,NPRO
        IF(TPRO(IPRO))NCORR=NCORR+1
      ENDDO
!
!     ==========================================================================
!     == SELECT CORRELATED BANDS                                              ==
!     ==========================================================================
      NB2=HUGE(NB2)    !NB2 MUST INCLUDE ALL UNOCCUPIED BANDS
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          NB=THIS%NB
          NB2=MIN(NB2,NB)
        ENDDO
      ENDDO
      NBW=NB2-NB1+1  !#(BANDS IN THE WINDOW)

PRINT*,'INFO FROM LMTO_DROPPICK_INI SWITCHID    ',SWITCHID
PRINT*,'INFO FROM LMTO_DROPPICK_INI NB1,NB2,NBW ',NB1,NB2,NBW 
PRINT*,'INFO FROM LMTO_DROPPICK_INI NCORR       ',NCORR
!
!     ==========================================================================
!     ==  ATTACH FILE                                                         ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('DMFT2DFT',.FALSE.,-'DMFT2DFT.DAT')
      CALL FILEHANDLER$SETSPECIFICATION('DMFT2DFT','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('DMFT2DFT','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('DMFT2DFT','FORM','FORMATTED')
!
      CALL FILEHANDLER$SETFILE('DMFT2DFTDUMMY',.FALSE.,-'DMFT2DFT.DAT')
      CALL FILEHANDLER$SETSPECIFICATION('DMFT2DFTDUMMY','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('DMFT2DFTDUMMY','FORM','FORMATTED')
!
!     == THIS DOES NOT WORK, BECAUSE THE FILEHANDLER OBJECT SETS A FIXED =======
!     == RECORDLENGTH OF 1000 DIGITS FOR FORMATTED FILES TO AVOID UNINTEDED ====
!     == LINE BREAKS. FINALLY THE FILE SHOULD BE WRITTEN UNFORMATTED  ==========
      CALL FILEHANDLER$SETFILE('DFT2DMFT',.FALSE.,-'DFT2DMFT.DAT')
      CALL FILEHANDLER$SETSPECIFICATION('DFT2DMFT','FORM','FORMATTED')

      CALL FILEHANDLER$SETFILE('DHOFK_OUT',.FALSE.,-'DHOFK.DAT')
      CALL FILEHANDLER$SETSPECIFICATION('DHOFK_OUT','FORM','UNFORMATTED')
      CALL FILEHANDLER$SETFILE('DHOFK_IN',.FALSE.,-'DHOFK.DAT')
      CALL FILEHANDLER$SETSPECIFICATION('DHOFK_IN','FORM','UNFORMATTED')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPPICK$SETL4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_DROPPICK_MODULE, ONLY : TREADDHOFK
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'DHOFK') THEN
        TREADDHOFK=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO_DROPPICK$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPPICK_DROP()
!     **************************************************************************
!     **  WRITES A FILE WITH KOHN-SHAM WAVE FUNCTIONS AND ENERGIES            **
!     **  FOR USE AS DMFT INTERFACE                                           **
!     **                                                                      **
!     **    <PI_ALPHA|PSI_N> ; E_N                                            **
!     **                                                                      **
!     **  THE BASIS ARE ONSITE-ORTHOGONALIZED NATURAL TIGHT-BINDING ORBITALS  **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE STRINGS_MODULE
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : TDROP,LOX,LNX,ISPECIES
      USE LMTO_DROPPICK_MODULE, ONLY : NCORR,NB1,NB2,NBW,TPRO,DHOFK,T,TICKET &
     &                                ,IPRO1,NPROAT
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      LOGICAL(4),PARAMETER   :: TDUMMY=.FALSE. !WRITE A DUMMY DMFT2DFT FILE
      INTEGER(4)             :: NPRO
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: NFIL,NFIL2
      INTEGER(4)             :: NB,NBH
      REAL(8)                :: MU
      REAL(8)                :: KBT
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      COMPLEX(8),ALLOCATABLE :: PSI(:,:,:)
      COMPLEX(8),ALLOCATABLE :: PSI1(:,:,:)
      COMPLEX(8),ALLOCATABLE :: H0(:,:)
      COMPLEX(8),ALLOCATABLE :: H0PLUSDELTA(:,:)
      COMPLEX(8),ALLOCATABLE :: DELTAH(:,:)
      COMPLEX(8),ALLOCATABLE :: QSQ(:,:)
      COMPLEX(8),ALLOCATABLE :: QSQINV(:,:)
      REAL(8)   ,ALLOCATABLE :: WKPT(:)
      CHARACTER(16)          :: ID
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NKPT
      INTEGER(4)             :: NBX
      INTEGER(4)             :: LMNXT,LMNX
      INTEGER(4)             :: I,J,I1,I2,IDIM,IPRO,IAT,ISP,IKPT,ISPIN 
      INTEGER(4)             :: IB,JB,IBH
      INTEGER(4)             :: IWORK16(16)
      REAL(8)                :: EV
      REAL(8)                :: SVAR
      INTEGER(4)             :: L,M,LN
      COMPLEX(8)             :: CSVAR
      COMPLEX(8),ALLOCATABLE :: AMAT(:,:),UMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: PSICORR(:,:,:)
      REAL(8)   ,ALLOCATABLE :: EMAT(:)
      REAL(8)                :: NEL  !#(ELECTRONS IN SELECTED BANDS)
!     **************************************************************************
                                           CALL TRACE$PUSH('LMTO_DROPPICK_DROP')
PRINT*,'ENTERING LMTO_DROPPICK_DROP'
      IF(.NOT.TDROP) RETURN
      CALL CONSTANTS('EV',EV)
      NAT=SIZE(ISPECIES)
      CALL DYNOCC$GETR8('TEMP',KBT)
      CALL DYNOCC$GETR8('EFERMI',MU)
      IF(NDIM.NE.1) THEN
        CALL ERROR$MSG('LMTO_DROPPICK DOES NOT ALLOW NONCOLLINEAR CALCULATIONS')
        CALL ERROR$STOP('LMTO_DROPPICK_DROP')
      END IF
!
!     ==========================================================================
!     ==  DEFINE WINDOW OF BANDS AND CHOICE OF LOCAL ORBITALS                 ==
!     ==========================================================================
      CALL LMTO_DROPICK_INI()
      CALL LMTO_DROPPICK_MAKET()
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      CALL LMTO_DROPPICK_NEL(NB1,NB2,NEL)
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      NKPT=NKPTL
      IF(NKPT.NE.NKPTL) THEN
        CALL ERROR$MSG('ROUTINE NOT SUITED FOR PARALLEL CALCULATIONS')
        CALL ERROR$STOP('LMTO_DROPPICK_DROP')
      END IF
      ALLOCATE(WKPT(NKPT))
      CALL DYNOCC$GETR8A('WKPT',NKPT,WKPT)
!     == MULTIPLY WITH SPIN MULTIPLICITY =======================================
      IF(NSPIN.EQ.2) WKPT=2.D0*WKPT
!
! T IS NO MORE CALCULATED BECAUSE IT IS ALREADY CALCULATED IN DROPPICK_MAKET
!
!     ==========================================================================
!     ==  ATTACH FILE                                                         ==
!     ==========================================================================
!      CALL FILEHANDLER$UNIT('DFT2DMFT',NFIL)
NFIL=12
OPEN(NFIL,FILE=-'DFT2DMFT.DAT',FORM='FORMATTED',RECL=10000000)
      REWIND(NFIL)
!!$!FILEHANDLER DOES NOT WORK FOR FORTMATTED FILES BECAUSE OF LIMITED LINE LENGTH
!
!     ==========================================================================
!     == WRITE INFO TO FILE
!     ==========================================================================
      CALL DATE_AND_TIME(VALUES=TICKET)
      ID='INFO'
      WRITE(NFIL,*)ID,NAT,NDIM,NSPIN,NKPT,NCORR,NBW,MU,KBT,NEL,TICKET(:) !<<<<<<
PRINT*,'NBW  AAA',NBW,' NCORR ',NCORR
      IF(TDUMMY) THEN
!        CALL FILEHANDLER$UNIT('DMFT2DFTDUMMY',NFIL2)
NFIL2=14
OPEN(NFIL2,FILE=-'DMFT2DFT.DAT')
        REWIND(NFIL2)
        ID='INFO'
        WRITE(NFIL2,*)ID,NKPT,NSPIN,NBW,TICKET   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      END IF
!
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        IWORK16(:)=-1
        IF(LNX(ISP).GT.16) THEN
          CALL ERROR$MSG('INTERNAL LIMIT FOR VARIABLE LNX EXCEEDED')
          CALL ERROR$STOP('LMTO_DROP')
        END IF
        IWORK16(:LNX(ISP))=LOX(:LNX(ISP),ISP)
        ID='PSIINFO'
        WRITE(NFIL,*)ID,IAT,IPRO1(IAT),NPROAT(IAT),LNX(ISP),IWORK16 !<<<<<<<<<<
      ENDDO
!
      ID='ORBINFO'
      I=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          DO M=1,2*L+1
            I=I+1
            IF(TPRO(I))WRITE(NFIL,*)ID,I,IAT,LN,L,M !<<<<<<<<<<<<<<<<<<<<<<<<<<
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  COLLECT OCCUPATIONS                                                 ==
!     ==========================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
!
!     ==========================================================================
!     ==  WRITE WAVE FUNCTIONS                                                ==
!     ==========================================================================
      ALLOCATE(PSICORR(NDIM,NCORR,NBW))
      ALLOCATE(QSQ(NCORR,NCORR))
      ALLOCATE(QSQINV(NCORR,NCORR))
      ALLOCATE(H0(NBW,NBW))
      ALLOCATE(H0PLUSDELTA(NBW,NBW))
      ALLOCATE(DELTAH(NBW,NBW))
      NPRO=MAP%NPRO
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          NBH=THIS%NBH
          NB=THIS%NB
!
!         ==CHECKS =============================================================
          IF(.NOT.ASSOCIATED(THIS%TBC)) THEN
            CALL ERROR$MSG('THIS%TBC NOT ASSOCIATED')
            CALL ERROR$STOP('LMTO_DROPPICK_DROP')
          END IF
          IF(.NOT.ASSOCIATED(THIS%RLAM0)) THEN
            CALL ERROR$MSG('THIS%RLAM0 NOT ASSOCIATED')
            CALL ERROR$STOP('LMTO_DROPPICK_DROP')
          END IF
!
!         ==  COPY COEFFICIENTS OF LOCAL ORBITALS ONTO PSI =====================
!         ==  RESOLVE SUPER WAVE FUNCTIONS =====================================
          ALLOCATE(PSI(NDIM,NPRO,NB))
          ALLOCATE(PSI1(NDIM,NPRO,NB))
          IF(TINV) THEN
            IF(NBH.EQ.NB) THEN
              CALL ERROR$MSG('SUPER WAVE FUNCTIONS NOT PROPERLY ACCOUNTED FOR')
              CALL ERROR$STOP('LMTO_DROPPICK_DROP')
            END IF
            DO IBH=1,NBH
              PSI(:,:,2*IBH-1)=REAL(THIS%TBC(:,IBH,:))
              PSI(:,:,2*IBH)=AIMAG(THIS%TBC(:,IBH,:))
            ENDDO
          ELSE
            DO IB=1,NB
              PSI(:,:,IB)=THIS%TBC(:,IB,:)
            ENDDO
          END IF
!
!         ======================================================================
!         == TRANSFORM TO ONSITE-ORTHOGONAL ORBITALS                          ==
!         ======================================================================
          PSI1=PSI
          DO IAT=1,NAT
            I1=IPRO1(IAT)
            I2=I1-1+NPROAT(IAT)
            DO IB=1,NB
              DO IDIM=1,NDIM
                PSI(IDIM,I1:I2,IB)=MATMUL(T(IAT)%INV,PSI1(IDIM,I1:I2,IB))
              ENDDO
            ENDDO
          ENDDO
!
!         ======================================================================
!         == CONTRACT WAVE FUNCTION COEFFICIENTS ONTO CORRELATED ORBITALS     ==
!         == AND CORRELATED BANDS                                             ==
!         ======================================================================
          I=0
          DO IPRO=1,NPRO
            IF(.NOT.TPRO(IPRO)) CYCLE
            I=I+1
            PSICORR(:,I,:)=PSI(:,IPRO,NB1:NB2)
          ENDDO
!
!         ======================================================================
!         == ENFORCE SUM RULE                                                 ==
!         ======================================================================
          CALL LMTO_DROPPICK_SUMRULE(NDIM,NCORR,NBW,PSICORR,QSQ,QSQINV)
          DO IB=1,NBW
            DO IDIM=1,NDIM
              PSICORR(IDIM,:,IB)=MATMUL(QSQINV,PSICORR(IDIM,:,IB))
            ENDDO
          ENDDO
!
!         ======================================================================
!         ==  CHECK SUM RULE                                                  ==
!         ======================================================================
          IF(TTEST) THEN
            DO I=1,NCORR
              DO J=1,NCORR
                CSVAR=(0.D0,0.D0)
                IF(I.EQ.J)CSVAR=(-1.D0,0.D0)
                DO IB=1,NBW
                  DO IDIM=1,NDIM
                    CSVAR=CSVAR+PSICORR(IDIM,I,IB)*CONJG(PSICORR(IDIM,J,IB))
                  ENDDO
                ENDDO
                IF(ABS(CSVAR).GT.1.D-6) THEN
                  CALL ERROR$MSG('SUM RULE IS VIOLATED')
                  CALL ERROR$I4VAL('INDEX1',I)
                  CALL ERROR$I4VAL('INDEX2',J)
                  CALL ERROR$R8VAL('REAL(DEVIATION)',REAL(CSVAR))
                  CALL ERROR$R8VAL('IMAG(DEVIATION)',AIMAG(CSVAR))
                  CALL ERROR$STOP('LMTO_DROPPICK_DROP')
                END IF 
              ENDDO
            ENDDO
          END IF
!
!         ======================================================================
!         == COLLECT EFFECTIVE HAMILTONIAN                                    ==
!         ======================================================================
          H0PLUSDELTA(:,:)=THIS%RLAM0(NB1:NB2,NB1:NB2)
          H0PLUSDELTA(:,:)=0.5D0*(H0PLUSDELTA(:,:)+TRANSPOSE(CONJG(H0PLUSDELTA)))
!
!         ======================================================================
!         == DETERMINE NON-INTERACTING HAMILTONIAN                            ==
!         ======================================================================
          IF(ASSOCIATED(DHOFK)) THEN
            IDIM=1
            QSQ=MATMUL(QSQ,MATMUL(DHOFK(:,:,IKPT,ISPIN),QSQ))
            DELTAH=MATMUL(TRANSPOSE(CONJG(PSICORR(IDIM,:,:))) &
         &           ,MATMUL(QSQ,PSICORR(IDIM,:,:)))
          ELSE
            DELTAH=0.D0
          END IF
          H0=H0PLUSDELTA-DELTAH
!
!         ======================================================================
!         == WRITE WAVE FUNCTION TO FILE                                      ==
!         ======================================================================
          ID='H0INFO'
          WRITE(NFIL,*)ID,IKPT,ISPIN,WKPT(IKPT),H0(:,:) !<<<<<<<<<<<<<<<<<<<<<<<
          ID='H0PLUSDELTAINFO'
          WRITE(NFIL,*)ID,IKPT,ISPIN,WKPT(IKPT),H0PLUSDELTA !<<<<<<<<<<<<<<<<<<<
          ID='QSQINV'
          WRITE(NFIL,*)ID,IKPT,ISPIN,QSQINV(:,:) !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ID='PIPSI'
          DO IB=1,NBW
            WRITE(NFIL,*)ID,IB,IKPT,ISPIN,PSICORR(:,:,IB) !<<<<<<<<<<<<<<<<<<<<<
          ENDDO
          DEALLOCATE(PSI)
          DEALLOCATE(PSI1)
!
!         ======================================================================
!         == WRITE DUMMY INPUT FILE                                           ==
!         ======================================================================
          IF(TDUMMY) THEN
            H0=(0.D0,0.D0)
            DO IB=1,NBW
              SVAR=REAL(THIS%RLAM0(NB1-1+IB,NB1-1+IB),KIND=8)
              SVAR=EXP((SVAR-MU)/KBT)
              SVAR=1.D0/(1.D0+SVAR)
              H0(IB,IB)=SVAR
            ENDDO
            ID='RHO'
            WRITE(NFIL2,*)ID,IKPT,ISPIN,H0(:,:) !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          END IF
        ENDDO ! END LOOP ISPIN
      ENDDO   ! END LOOP IKPT
      DEALLOCATE(PSICORR)
      DEALLOCATE(QSQ)
      DEALLOCATE(QSQINV)
      DEALLOCATE(H0)
CLOSE(NFIL)
!      CALL FILEHANDLER$CLOSE('DFT2DMFT')
!      IF(TDUMMY) CALL FILEHANDLER$CLOSE('DMFT2DFTDUMMY')
      IF(TDUMMY) CLOSE(NFIL2)
      CALL ERROR$MSG('FORCED STOP AFTER WRITING FILE')
      CALL ERROR$STOP('LMTO_DROPPICK_DROP')
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPPICK_PICK()
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES
      USE LMTO_DROPPICK_MODULE, ONLY : TICKET,NCORR,NBW,TPRO,TPICKED &
     &                                ,IPRO1,NPROAT,DHOFK,NB1,NB2,T
      USE WAVES_MODULE, ONLY : THIS,WAVES_SELECTWV
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TREFRESH=.FALSE. ! RECALCULATE DELTA-H EACH STEP
!      CHARACTER(128),PARAMETER :: TYPE='CONSTRAINEDSEARCH'
      CHARACTER(128),PARAMETER :: TYPE='NONE'
      REAL(8)                :: MU
      REAL(8)                :: KBT
      CHARACTER(16)          :: ID
      INTEGER(4)             :: TICKET1(8)
      INTEGER(4)             :: NFIL1,NFIL2
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NKPT,NSPIN,NDIM
      INTEGER(4)             :: NKPT1,NSPIN1,NBW1,NDIM1
      INTEGER(4)             :: NCORR1
      INTEGER(4)             :: NBH,NB,NBX
      COMPLEX(8),ALLOCATABLE :: RHO(:,:)
      COMPLEX(8),ALLOCATABLE :: RHO0(:,:)
      COMPLEX(8),ALLOCATABLE :: U(:,:)
      REAL(8)   ,ALLOCATABLE :: F(:)
      COMPLEX(8),ALLOCATABLE :: H0(:,:)
      COMPLEX(8),ALLOCATABLE :: DELTAH(:,:)
      COMPLEX(8),ALLOCATABLE :: DEVRHO(:,:)
      COMPLEX(8),ALLOCATABLE :: H0PLUSDELTA(:,:)
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      COMPLEX(8),ALLOCATABLE :: PSICORR(:,:,:)
      COMPLEX(8),ALLOCATABLE :: QSQ(:,:),QSQINV(:,:)
      COMPLEX(8),ALLOCATABLE :: LAMBDAP(:,:)
      COMPLEX(8),ALLOCATABLE :: LAMBDA0(:,:)
      COMPLEX(8),ALLOCATABLE :: LAMBDAM(:,:)
      REAL(8)                :: SVAR1,SVAR2,SVAR3
      REAL(8)   ,PARAMETER   :: ALAMBDA=0.1D0
      REAL(8)   ,PARAMETER   :: MLAMBDA=1.D0
      REAL(8)   ,PARAMETER   :: DELTA=5.D0
      REAL(8)                :: SVAR
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IKPT,ISPIN,IB,IDIM,IPRO,IBH,IAT,I,J,I1,I2
      INTEGER(4)             :: IKPT1,ISPIN1,IB1
      INTEGER(4)             :: NPRO
      REAL(8)                :: WKPT1
      REAL(8)                :: NEL
!     **************************************************************************
      IF(TPICKED) RETURN
                                           CALL TRACE$PUSH('LMTO_DROPPICK_PICK')
      CALL LMTO_DROPICK_INI()
      CALL LMTO_DROPPICK_MAKET()
!
!     ==========================================================================
!     ==  PROPAGATE DFT FILE TO THE BEGINNING OF THE K-POINT LOOP             ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('DFT2DMFT',NFIL1)
      REWIND(NFIL1)
      READ(NFIL1,*)ID,NAT,NDIM,NSPIN,NKPT,NCORR1,NBW1,MU,KBT,NEL,TICKET1 !<<<<<<
!     == SET TICKET IF IT HAS NOT BEEN SET =====================================
PRINT*,'NBW1 ON DFT2DMFT ',NBW1,NBW
PRINT*,'A TICKET ',TICKET
PRINT*,'A TICKET1',TICKET1
      IF(SUM(ABS(TICKET)).EQ.0) THEN
        TICKET=TICKET1
      END IF
!     == CHECK CONSISTENCY OF THE TICKET =======================================
      IF(.NOT.ALL(TICKET.EQ.TICKET1)) THEN
        CALL ERROR$MSG('INCORRECT TICKET (1)')
PRINT*,'B TICKET ',TICKET
PRINT*,'B TICKET1',TICKET1
        CALL ERROR$STOP('LMTO_DROPPICK_PICK')
      END IF
      DO IAT=1,NAT
        READ(NFIL1,*)ID
      ENDDO
      DO I=1,NCORR
        READ(NFIL1,*)ID
      ENDDO
!
!     ==========================================================================
!     ==  READ INPUT DMFT FILE                                                ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('DMFT2DFT',NFIL2)
      REWIND(NFIL2)
      READ(NFIL2,*)ID,NKPT1,NSPIN1,NBW1,TICKET1 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
PRINT*,'NBW1 ON DMFT2DFT ',NBW1,NBW

      IF(ID.NE.'INFO') THEN
        CALL ERROR$MSG('INCORRECT ID: MUST BE "INFO"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO_DROPPICK_PICK')
      END IF
      IF(.NOT.ALL(TICKET.EQ.TICKET1)) THEN
        CALL ERROR$MSG('INCORRECT TICKET')
PRINT*,'C TICKET ',TICKET
PRINT*,'C TICKET1',TICKET1
        CALL ERROR$STOP('LMTO_DROPPICK_PICK')
      END IF
      IF(NKPT1.NE.NKPT.OR.NSPIN1.NE.NSPIN.OR.NBW1.NE.NBW) THEN
        CALL ERROR$MSG('INCONSISTENT DATA')
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$I4VAL('NSPIN1',NSPIN1)
        CALL ERROR$I4VAL('NKPT1',NKPT1)
        CALL ERROR$I4VAL('NKPTL',NKPT)
        CALL ERROR$I4VAL('NBW',NBW)
        CALL ERROR$I4VAL('NBW1',NBW1)
        CALL ERROR$STOP('LMTO_DROPPICK_PICK')
      END IF
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPT,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPT*NSPIN,OCC)
      IF(.NOT.ASSOCIATED(DHOFK)) ALLOCATE(DHOFK(NCORR,NCORR,NKPT,NSPIN))
      ALLOCATE(RHO(NBW,NBW))
      ALLOCATE(RHO0(NBW,NBW))
      ALLOCATE(H0(NBW,NBW))
      ALLOCATE(H0PLUSDELTA(NBW,NBW))
      ALLOCATE(DELTAH(NBW,NBW))
PRINT*,'SHAPE  OF DELTAH  (1) :',SHAPE(DELTAH),' NBW ',NBW
      ALLOCATE(F(NBW))
      ALLOCATE(U(NBW,NBW))
      ALLOCATE(DEVRHO(NCORR,NCORR))
      ALLOCATE(LAMBDAP(NCORR,NCORR))
      ALLOCATE(LAMBDA0(NCORR,NCORR))
      ALLOCATE(LAMBDAM(NCORR,NCORR))
      ALLOCATE(PSICORR(NDIM,NCORR,NBW))
      ALLOCATE(QSQ(NCORR,NCORR))
      ALLOCATE(QSQINV(NCORR,NCORR))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
PRINT*,'SHAPE  OF DELTAH  (1A) :',SHAPE(DELTAH),' NBW ',NBW
!
!         ======================================================================
!         == READ DENSITY MATRIX                                              ==
!         ======================================================================
          READ(NFIL2,*)ID,IKPT1,ISPIN1,RHO(:,:) !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          IF(ID.NE.'RHO') THEN
            CALL ERROR$MSG('INCORRECT ID: MUST BE "RHO"')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('LMTO_DROPPICK_PICK')
          END IF
          IF(IKPT1.NE.IKPT.OR.ISPIN1.NE.ISPIN) THEN
            CALL ERROR$MSG('INDEX CONFUSION')
            CALL ERROR$CHVAL('IKPT',IKPT)
            CALL ERROR$CHVAL('ISPIN',ISPIN)
            CALL ERROR$STOP('LMTO_DROPPICK_PICK')
          END IF
PRINT*,'SHAPE  OF DELTAH  (1B) :',SHAPE(DELTAH),' NBW ',NBW
!
!         ======================================================================
!         == READ DFT FILE                                                    ==
!         ======================================================================
PRINT*,'MARKE 1',IKPT,ISPIN,NBW
READ(NFIL1,*)ID,IKPT1,ISPIN1,WKPT1
BACKSPACE(NFIL1)
PRINT*,'MARKE 1A',TRIM(ID),IKPT1,ISPIN1,WKPT1
          READ(NFIL1,*)ID,IKPT1,ISPIN1,WKPT1,H0 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          IF(ID.NE.'H0INFO') THEN
            CALL ERROR$MSG('INCORRECT ID: MUST BE "H0INFO"')
            CALL ERROR$MSG('OLD FORMAT USED HINFO')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('LMTO_DROPPICK_PICK')
          END IF
          READ(NFIL1,*)ID,IKPT1,ISPIN1,WKPT1,H0PLUSDELTA !<<<<<<<<<<<<<<<<<<<<<<
          IF(ID.NE.'H0PLUSDELTAINFO') THEN
            CALL ERROR$MSG('INCORRECT ID: MUST BE "H0PLUSDELTAINFO"')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('LMTO_DROPPICK_PICK')
          END IF
PRINT*,'SHAPE  OF DELTAH  (1C) :',SHAPE(DELTAH),' NBW ',NBW
!
          READ(NFIL1,*)ID,IKPT1,ISPIN1,QSQINV(:,:) !<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          IF(ID.NE.'QSQINV') THEN
            CALL ERROR$MSG('INCORRECT ID: MUST BE "QSQINV"')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('LMTO_DROPPICK_PICK')
          END IF
PRINT*,'SHAPE  OF DELTAH  (1D) :',SHAPE(DELTAH),' NBW ',NBW
!
          DO IB=1,NBW
            READ(NFIL1,*)ID,IB1,IKPT1,ISPIN1,PSICORR(:,:,IB) !<<<<<<<<<<<<<<<<<<
            IF(ID.NE.'PIPSI') THEN
              CALL ERROR$MSG('INCORRECT ID: MUST BE "PIPSI"')
              CALL ERROR$CHVAL('ID',ID)
              CALL ERROR$STOP('LMTO_DROPPICK_PICK')
            END IF
          ENDDO
PRINT*,'SHAPE  OF DELTAH  (1E) :',SHAPE(DELTAH),' NBW ',NBW

!         == NON-SPIN-POLARIZED CALCULATIONS HAVE TWO ELECTRONS PER STATE
!!$          IF(NSPIN.EQ.1.AND.NDIM.EQ.1) THEN
!!$            RHO=0.5D0*RHO
!!$          END IF
!
!         ======================================================================
!         == CONVERT DENSITY MATRIX INTO FINITE TEMPERATURE HAMILTONIAN       ==
!         ======================================================================
PRINT*,'SHAPE  OF DELTAH  (1F) :',SHAPE(DELTAH),' NBW ',NBW
          IF(TYPE.EQ.'CONSTRAINEDSEARCH') THEN
            F(:)=OCC(NB1:,IKPT,ISPIN)/WKPT1
!SET UP MATRIX
!           == VIOLATION OF THE CONSTRAINT =====================================
            IDIM=1
            DEVRHO=0.D0
            DO I=1,NCORR
              DO J=1,NCORR
                DO IB=1,NBW
!HERE SHOULD BE THE PSICORR OF THE ACTUAL TIME STEP
                  DEVRHO(I,J)=DEVRHO(I,J) &
                             +PSICORR(IDIM,I,IB)*F(IB)*CONJG(PSICORR(IDIM,J,IB))
                ENDDO
              ENDDO
            ENDDO
!THE REFERENCE PSICORR SHOULD REMAIN CONSTANT WITHIN THE DFT-LOOP
            DEVRHO(:,:)=DEVRHO-MATMUL(PSICORR(IDIM,:,:), &
                               MATMUL(RHO,CONJG(TRANSPOSE(PSICORR(IDIM,:,:)))))
!
!           ==  PROPAGATE LAGRANGE MULTIPLIERS==================================
            SVAR1=2.D0/(1.D0+ALAMBDA)
            SVAR2=1.D0-SVAR1
            SVAR3=DELTA**2/MLAMBDA/(1.D0+ALAMBDA)
            LAMBDAP=LAMBDA0*SVAR1+LAMBDAM*SVAR2+DEVRHO*SVAR3
!
!           == SET DELTAH ======================================================
PRINT*,'SHAPE  OF DELTAH  (1F1) :',SHAPE(DELTAH),' NBW ',NBW
STOP 'PROGRAMM ERROR: LAMBDA0 AND DELTAH HAVE DIFFERENT SHAPE'
            DELTAH=LAMBDA0
PRINT*,'SHAPE  OF DELTAH  (1F2) :',SHAPE(DELTAH),' NBW ',NBW
!
!           == SWITCH ==========================================================
            LAMBDAM=LAMBDA0
            LAMBDA0=LAMBDAP
          ELSE 
            CALL LIB$DIAGC8(NBW,RHO,F,U)
            DO IB=1,NBW
              F(IB)=MAX(1.D-5,MIN(1.D0-1.D-5,F(IB)))
              SVAR=(1.D0-F(IB))/F(IB)
              SVAR=LOG(SVAR)
              F(IB)=MU+KBT*SVAR
            ENDDO
            DO IB=1,NBW
              RHO(IB,:)=F(IB)*U(IB,:)
            ENDDO
            RHO=MATMUL(CONJG(TRANSPOSE(U)),RHO)  !THIS IS NOW A HAMILTONIAN
!
!           ====================================================================
!           == SUBTRACT NON-INTERACTING HAMILTONIAN                           ==
!           ====================================================================
PRINT*,'SHAPE  OF DELTAH (2)  :',SHAPE(DELTAH),' NBW ',NBW
            DELTAH=RHO-H0
PRINT*,'SHAPE  OF DELTAH (3)  :',SHAPE(DELTAH),' NBW ',NBW
!!$            IF(TREFRESH) THEN
!!$CALL ERROR$STOP('LMTO_DROPPICK_PICK')
!!$              H0=DELTAH+H0PLUSDELTA
!!$              H0PLUSDELTA(:,:)=THIS%RLAM0(NB1:,NB1:)
!!$              H0PLUSDELTA(:,:)=(H0PLUSDELTA(:,:)+TRANSPOSE(CONJG(H0PLUSDELTA)))
!!$              H0=H0-H0PLUSDELTA
!!$            END IF
          END IF
PRINT*,'SHAPE  OF DELTAH  (4) :',SHAPE(DELTAH),' NBW ',NBW
!
!         ======================================================================
!         == CONVERT INTO ORBITAL BASIS                                       ==
!         ======================================================================
          IDIM=1
PRINT*,'DIMENSIONS PSICORR :',NDIM,NCORR,NBW
PRINT*,'DIMENSIONS DELTAH  :',NBW,NBW
PRINT*,'DIMENSIONS DHOFK   :',NCORR,NCORR,NKPT,NSPIN
PRINT*,'SHAPE  OF DHOFK    :',SHAPE(DHOFK)
PRINT*,'SHAPE  OF DELTAH (5)  :',SHAPE(DELTAH)
PRINT*,'SHAPE  OF PSICORR  :',SHAPE(PSICORR)
PRINT*,'SHAPE  OF TRANSPOSE(CONJG(PSICORR(IDIM,:,:))) :',SHAPE(TRANSPOSE(CONJG(PSICORR(IDIM,:,:))))
PRINT*,'SHAPE  OF MATMUL(DELTAH,TRANSPOSE(CONJG(PSICORR(IDIM,:,:)))) :' &
       ,SHAPE(MATMUL(DELTAH,TRANSPOSE(CONJG(PSICORR(IDIM,:,:)))))
PRINT*,'MARKE 5'
          DHOFK(:,:,IKPT,ISPIN)=MATMUL(PSICORR(IDIM,:,:) &
     &                      ,MATMUL(DELTAH,TRANSPOSE(CONJG(PSICORR(IDIM,:,:)))))
PRINT*,'MARKE 6'
!
!         ======================================================================
!         ==  TRANSFORM DHOFK TO THE SITE ORTHOGONALIZED REPRESENTATION       ==
!         ======================================================================
          DHOFK(:,:,IKPT,ISPIN)=MATMUL(QSQINV &
     &                                ,MATMUL(DHOFK(:,:,IKPT,ISPIN),QSQINV))
PRINT*,'MARKE 7'
        ENDDO
      ENDDO
      CALL FILEHANDLER$CLOSE('DFT2DMFT')
      CALL FILEHANDLER$CLOSE('DMFT2DFT')
!
!     ==========================================================================
!     ==  WRITE DHOFK TO FILE SO THAT ONE CAN RESTART                         ==
!     ==========================================================================
      CALL LMTO_DROPPICK_WRITEDHOFK()
      TPICKED=.TRUE.
      IF(TREFRESH) TPICKED=.FALSE.
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPPICK_HTBC()
!     **************************************************************************
!     **  READS A FILE HAMILTON CORRECTION                                    **
!     **  FOR USE AS DMFT INTERFACE                                           **
!     **                                                                      **
!     **    <PI_ALPHA|PSI_N> ; E_N                                            **
!     **                                                                      **
!     **  THE BASIS ARE ONSITE-ORTHOGONALIZED NATURAL TIGHT-BINDING ORBITALS  **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET,WVSET_TYPE
      USE LMTO_MODULE,  ONLY : TPICK,ISPECIES,LNX,LOX,THTBC
      USE LMTO_DROPPICK_MODULE, ONLY : TPRO,NCORR,T,DHOFK,TREADDHOFK,NPROAT
      IMPLICIT NONE
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      INTEGER(4)             :: NFIL
      INTEGER(4)             :: NB,NBH,NBX
      INTEGER(4)             :: NKPT,NKPT1
      INTEGER(4)             :: NGL
      INTEGER(4)             :: NAT
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: IAT,ISP
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      COMPLEX(8),ALLOCATABLE :: VEC1(:)
      COMPLEX(8),ALLOCATABLE :: VEC2(:)
      REAL(8)                :: F1,F2
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IB,IBH,IKPT,ISPIN,I,J,IDIM,IPRO
      INTEGER(4)             :: I1,I2
      INTEGER(4)             :: NPRO
      REAL(8)                :: ETOT
      COMPLEX(8)             :: CSVAR,CSVAR1,CSVAR2
      INTEGER(4)             :: NFILO
!     **************************************************************************
                                           CALL TRACE$PUSH('LMTO_DROPPICK_HTBC')
      IF(.NOT.TPICK) RETURN
      THTBC=.TRUE.
      CALL LMTO_DROPICK_INI()
      CALL LMTO_DROPPICK_MAKET()
      NAT=SIZE(ISPECIES)
      NPRO=SUM(NPROAT)
!
!     ==========================================================================
!     == READ OUTPUT OF DMFT CODE AND STORE THE AHMILTONIAN IN                ==
!     == CORE-ORTHOGONALIZED ORBITALS.                                        ==
!     ==========================================================================
      IF(TREADDHOFK) THEN
        CALL LMTO_DROPPICK_READDHOFK()
      ELSE
        CALL LMTO_DROPPICK_PICK()
      END IF
!
!     ==========================================================================
!     ==  COLLECT OCCUPATIONS                                                 ==
!     ==========================================================================
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(OCC(NBX,NKPTL,NSPIN))
      CALL WAVES_DYNOCCGETR8A('OCC',NBX*NKPTL*NSPIN,OCC)
!
!     ==========================================================================
!     ==  CALCULATE HTBC                                                      ==
!     ==========================================================================
      ETOT=0.D0
      ALLOCATE(VEC1(NPRO))
      ALLOCATE(VEC2(NCORR))
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          NBH=THIS%NBH
          NB=THIS%NB
          NGL=GSET%NGL
          IF(.NOT.ASSOCIATED(THIS%HTBC))ALLOCATE(THIS%HTBC(NDIM,NBH,NPRO))
!
!         ======================================================================
!         == ADD TO FORCES ACTING ON WAVE FUNCTIONS                           ==
!         ======================================================================
          DO IBH=1,NBH
!
!           == MAKE COPY OF THIS%TBC ===========================================
!           == VEC1 IS A SUPER WAVE FUNCTION IF TINV=TRUE ======================
            VEC1=THIS%TBC(1,IBH,:)
            IF(NDIM.NE.1) THEN
              CALL ERROR$MSG('IMPLEMENTATION ONLY FOR NDIM=1')
              CALL ERROR$STOP('LMTO_DROPPICK_HTBC')
            END IF
!
!           == TRANSFORM TO ONSITE ORTHONORMALIZED STATES ======================
            DO IAT=1,NAT
              I1=T(IAT)%I1
              I2=T(IAT)%I2
              LMNX=I2-I1+1
              VEC1(I1:I2)=MATMUL(T(IAT)%INV,VEC1(I1:I2))
            ENDDO
!
!           == CONTRACT ========================================================
            I=0
            DO IPRO=1,NPRO
              IF(.NOT.TPRO(IPRO)) CYCLE
              I=I+1
              VEC2(I)=VEC1(IPRO)
            ENDDO
!
!           == MULTIPLY WITH HAMILTONIAN =======================================
            VEC2(:)=MATMUL(DHOFK(:,:,IKPT,ISPIN),VEC2)
!
!           == EXPAND ==========================================================
            VEC1(:)=(0.D0,0.D0)
            I=0
            DO IPRO=1,NPRO
              IF(.NOT.TPRO(IPRO)) CYCLE
              I=I+1
              VEC1(IPRO)=VEC2(I)
            ENDDO
!
!           == TRANSFORM FROM ONSITE ORTHONORMALIZED STATES ====================
            DO IAT=1,NAT
              I1=T(IAT)%I1
              I2=T(IAT)%I2
              LMNX=I2-I1+1
              VEC1(I1:I2)=MATMUL(TRANSPOSE(T(IAT)%INV),VEC1(I1:I2))
            ENDDO
!
!           == MAP RESULT INTO THIS%HTBC =======================================
            THIS%HTBC(1,IBH,:)=VEC1(:)
!
!           == ADD UP TOTAL ENERGY CORRECTION ==================================
            IF(TINV) THEN
              DO IDIM=1,NDIM
                F1=OCC(2*IBH-1,IKPT,ISPIN)
                F2=OCC(2*IBH  ,IKPT,ISPIN)
                CSVAR1=DOT_PRODUCT(THIS%TBC(IDIM,IBH,:),THIS%HTBC(IDIM,IBH,:))
                CSVAR2=DOT_PRODUCT(THIS%TBC(IDIM,IBH,:) &
         &                                      ,CONJG(THIS%HTBC(IDIM,IBH,:)))
                ETOT=ETOT+0.5D0*((F1+F2)*REAL(CSVAR1) &
         &                      +(F1-F2)*REAL(CSVAR2))
              ENDDO
            ELSE
              DO IDIM=1,NDIM
                CSVAR=DOT_PRODUCT(THIS%TBC(IDIM,IBH,:),THIS%HTBC(IDIM,IBH,:))
                ETOT=ETOT+OCC(IBH,IKPT,ISPIN)*REAL(CSVAR)
              ENDDO
            END IF 
          ENDDO ! END LOOP OVER IBH
!
        ENDDO !END LOOP OVER ISPIN
      ENDDO !END LOOP OVER IKPT
!
!     ==========================================================================
!     == REPORT TO ENERGYLIST                                                 ==
!     ==========================================================================
!!$CALL FILEHANDLER$UNIT('PROT',NFILO)
!!$WRITE(NFILO,*)'LMTO INTERFACE ',ETOT
      CALL ENERGYLIST$SET('LMTO INTERFACE',ETOT)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',ETOT)
!
                                           CALL TRACE$POP()
      RETURN
      END!
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPPICK_NEL(NB1,NB2,NEL)
!     **************************************************************************
!     ** COLLECTS THE NUMBER OF ELECTRONS IN THE SELECTED BANDS (W)           **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB1
      INTEGER(4),INTENT(IN) :: NB2
      REAL(8)   ,INTENT(OUT):: NEL
      INTEGER(4)            :: NKPT
      INTEGER(4)            :: NSPIN
      INTEGER(4)            :: NB
      REAL(8)   ,ALLOCATABLE:: OCC(:,:,:)
!     **************************************************************************
      CALL DYNOCC$GETI4('NKPT',NKPT)
      CALL DYNOCC$GETI4('NSPIN',NSPIN)
      CALL DYNOCC$GETI4('NB',NB)
      ALLOCATE(OCC(NB,NKPT,NSPIN))
      CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
      NEL=SUM(OCC(NB1:NB2,:,:))
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPPICK_SUMRULE(NDIM,NORB,NB,PIPSI,QSQ,QSQINV)
!     **************************************************************************
!     ** DETERMINES THE TRANSFORMATION TO A BASIS THAT SATISFIES THE          **
!     ** SUMRULE                                                              **
!     **   QPRIME(I,J)=SUM_N <PIPRIME_I|PSI_N><PSI_N|PIPRIME_J>=DELTA_{I,J}   **
!     ** THE NEW PROJECTORS ARE                                               **
!     **   |PIPRIME_I>=SUM_J |PI_J> QSQINV(J,I)                               **
!     ** THE NEW ORBITALS ARE                                                 **
!     **   |CHIPRIME_I>=SUM_J |CHI_J> QSQ(J,I)                                **
!     **                                                                      **
!     **  Q(I,J)=SUM_N <PI_I|PSI_N><PSI_N|PI_J>                               **
!     **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NDIM                ! #(SPINOR COMPONENTS)
      INTEGER(4),INTENT(IN) :: NORB                ! #(ORBITALS)
      INTEGER(4),INTENT(IN) :: NB                  ! #(BANDS)
      COMPLEX(8),INTENT(IN) :: PIPSI(NDIM,NORB,NB) ! <PI|PSI>
      COMPLEX(8),INTENT(OUT):: QSQ(NORB,NORB)      ! SQRT(<PI|PSI><PSI|PI>
      COMPLEX(8),INTENT(OUT):: QSQINV(NORB,NORB)   ! SQRT(<PI|PSI><PSI|PI>^(-1)
      LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
      COMPLEX(8)            :: UMAT(NORB,NORB)     ! EIGENVECTORS OF QMAT
      REAL(8)               :: EIG(NORB)           ! EIGENVALUES OF QMAT
      COMPLEX(8)            :: CSVAR
      INTEGER(4)            :: I,J,IB,IDIM
      REAL(8)               :: SVAR
      COMPLEX(8)            :: AMAT(NORB,NORB)     ! WORK ARRAY
      COMPLEX(8)            :: BMAT(NORB,NORB)     ! WORK ARRAY
COMPLEX(8) :: PIPSI1(NDIM,NORB,NB)
!     **************************************************************************
!
!     ==========================================================================
!     == DETERMINE SUM RULE VIOLATION                                         ==
!     ==========================================================================
      DO I=1,NORB
        DO J=I,NORB
          CSVAR=(0.D0,0.D0)
          DO IB=1,NB
            DO IDIM=1,NDIM
              CSVAR=CSVAR+PIPSI(IDIM,I,IB)*CONJG(PIPSI(IDIM,J,IB))
            ENDDO
          ENDDO
          QSQ(I,J)=CSVAR         ! THE ARRAY QSQ IS USED INSTEAD OF Q 
          QSQ(J,I)=CONJG(CSVAR)  ! TO SAVE MEMORY
        ENDDO
      ENDDO
      IF(TTEST)AMAT=QSQ
!
!     ==========================================================================
!     == DETERMINE TRANSFORMATION THAT ENFORCES THE SUM RULE                  ==
!     ==========================================================================
      CALL LIB$DIAGC8(NORB,QSQ,EIG,UMAT)
WRITE(*,FMT='("SUMRULE",2E20.5)')MINVAL(EIG),MAXVAL(EIG)
!
!     == CHECK POSITIVE DEFINITENESS ===========================================
      DO I=1,NORB
        IF(EIG(I).LT.0.D0) THEN
          CALL ERROR$MSG('EIGENVALUES ARE NOT POSITIVE')
          CALL ERROR$I4VAL('INDEX',I)
          CALL ERROR$R8VAL('EMAT',EIG(I))
          CALL ERROR$STOP('LMTO_DROPPICK_SUMRULE')
        END IF
      ENDDO
!
!     === TEST RECONSTRUCTION 
      IF(TTEST) THEN
        DO I=1,NORB
          QSQINV(I,:)=EIG(I)*CONJG(UMAT(:,I))  ! INTERMEDIATE RESULT
        ENDDO
        QSQINV(:,:)=MATMUL(UMAT,QSQINV)
        QSQINV(:,:)=QSQINV-QSQ
        SVAR=MAXVAL(ABS(QSQINV))
        PRINT*,'RECONSTRUCTION MAX DEV ',SVAR
        IF(SVAR.GT.1.D-8) THEN 
          CALL ERROR$MSG('RECONMSTRUCTION AFTER DIAGONALIZATION FAILED')
          CALL ERROR$R8VAL('MAX DEV',SVAR)
          CALL ERROR$STOP('LMTO_DROPPICK_SUMRULE')
        END IF
      END IF
!
!     ==========================================================================
!     == DETERMINE SQRT(Q) AND ITS INVERSE (Q IS THE SUMRULE VIOLATION)       ==
!     ==========================================================================
      DO I=1,NORB
        QSQINV(I,:)=SQRT(EIG(I))*CONJG(UMAT(:,I))  ! INTERMEDIATE RESULT
      ENDDO
      QSQ(:,:)=MATMUL(UMAT,QSQINV)
!
      IF(TTEST) THEN
        SVAR=MAXVAL(ABS(MATMUL(QSQ,QSQ)-AMAT))
        PRINT*,'MAX DEV ',SVAR
        IF(SVAR.GT.1.D-8) THEN 
          CALL ERROR$MSG('CONSTRUCTION OF SQUARE ROOT FAILED')
          CALL ERROR$R8VAL('MAX DEV',SVAR)
          CALL ERROR$STOP('LMTO_DROPPICK_SUMRULE')
        END IF
      END IF
!
!     ==========================================================================
!     == DETERMINE THE INVERSE OF SQRT(Q) (Q IS THE SUMRULE VIOLATION)        ==
!     ==========================================================================
      DO I=1,NORB
        QSQINV(I,:)=1.D0/SQRT(EIG(I)) * CONJG(UMAT(:,I))  ! INTERMEDIATE RESULT
      ENDDO
      QSQINV(:,:)=MATMUL(UMAT,QSQINV)
!
      IF(TTEST) THEN
        BMAT=MATMUL(QSQINV,MATMUL(AMAT,QSQINV))
        DO I=1,NORB
          BMAT(I,I)=BMAT(I,I)-(1.D0,0.D0)
        ENDDO
        SVAR=MAXVAL(ABS(BMAT))
        PRINT*,'MAX DEV ',SVAR
        IF(SVAR.GT.1.D-8) THEN 
          CALL ERROR$MSG('CONSTRUCTION OF INVERSE SQUARE ROOT FAILED')
          CALL ERROR$R8VAL('MAX DEV',SVAR)
          CALL ERROR$STOP('LMTO_DROPPICK_SUMRULE')
        END IF
      END IF
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        DO IB=1,NB
          DO IDIM=1,NDIM
            PIPSI1(IDIM,:,IB)=MATMUL(QSQINV,PIPSI(IDIM,:,IB))
          ENDDO
        ENDDO
        DO I=1,NORB
          DO J=1,NORB
            CSVAR=(0.D0,0.D0)
            IF(I.EQ.J) CSVAR=(-1.D0,0.D0)
            DO IB=1,NB
              DO IDIM=1,NDIM
                CSVAR=CSVAR+PIPSI1(IDIM,I,IB)*CONJG(PIPSI1(IDIM,J,IB))
              ENDDO
            ENDDO
            IF(ABS(CSVAR).GT.1.D-6) THEN
              CALL ERROR$MSG('FINAL TEST FAILED')
              CALL ERROR$R8VAL('MAX DEV',ABS(CSVAR))
              CALL ERROR$I4VAL('I',I)
              CALL ERROR$I4VAL('J',J)
              CALL ERROR$STOP('LMTO_DROPPICK_SUMRULE')
           END IF
          ENDDO
        ENDDO
      END IF
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPPICK_MAKET()
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES
      USE LMTO_DROPPICK_MODULE, ONLY : IPRO1,NPROAT,T
      IMPLICIT NONE
      LOGICAL(4),PARAMETER :: TWRITE=.FALSE.
      INTEGER(4)           :: NAT
      INTEGER(4)           :: IAT
      INTEGER(4)           :: ISP
      INTEGER(4)           :: LMNX
      INTEGER(4)           :: LMN
!     **************************************************************************
      CALL LMTO_DROPICK_INI()
      NAT=SIZE(ISPECIES)
!
!     ==========================================================================
!     ==  TRANSFORMATION ONTO ONSITE ORTHOGONALIZED PARTIAL WAVES             ==
!     ==  USES PRE-CALCULATED OVERLAP OF TAILED PARTIAL WAVES                 ==
!     ==========================================================================
      IF(.NOT.ASSOCIATED(T))ALLOCATE(T(NAT))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LMNX=NPROAT(IAT)
        T(IAT)%I1=IPRO1(IAT)
        T(IAT)%I2=IPRO1(IAT)-1+NPROAT(IAT)
        ALLOCATE(T(IAT)%MAT(LMNX,LMNX))
        ALLOCATE(T(IAT)%INV(LMNX,LMNX))
        CALL LMTO_ONSORTHO(IAT,LMNX,T(IAT)%MAT,T(IAT)%INV)
      ENDDO
!
!     ==========================================================================
!     ==  REPORT                                                              ==
!     ==========================================================================
      IF(TWRITE) THEN
        DO IAT=1,NAT
          WRITE(*,FMT='(82("="),T10,"ONSITE OVERLAP FOR ATOM ",I5)')IAT
          ISP=ISPECIES(IAT)
          LMNX=NPROAT(IAT)
          DO LMN=1,LMNX      
            WRITE(*,FMT='(20F10.5)')T(IAT)%MAT(LMN,:)
          ENDDO
        ENDDO
        CALL ERROR$MSG('FORCED STOP AFTER REPORTING')
        CALL ERROR$STOP('LMTO_DROPPICK_MAKET')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPPICK_READDHOFK()
!     **************************************************************************
!     ** READ HAMILTONIAN CORRECTION FROM FILE                                **
!     **************************************************************************
      USE LMTO_DROPPICK_MODULE, ONLY: DHOFK
      IMPLICIT NONE
      INTEGER(4) :: NFIL
      INTEGER(4) :: NKPT,NSPIN,NCORR
      INTEGER(4) :: IKPT,ISPIN,I
!     **************************************************************************
                                      CALL TRACE$PUSH('LMTO_DROPPICK_READDHOFK')
      CALL FILEHANDLER$UNIT('DHOFK_IN',NFIL)
      REWIND(NFIL)
      READ(NFIL)NKPT,NSPIN,NCORR
      IF(.NOT.ASSOCIATED(DHOFK)) ALLOCATE(DHOFK(NCORR,NCORR,NKPT,NSPIN))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          READ(NFIL)I,DHOFK(:,:,IKPT,ISPIN)
        ENDDO
      ENDDO
      CALL FILEHANDLER$CLOSE('DHOFK_IN')
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPPICK_WRITEDHOFK()
!     **************************************************************************
!     ** READ HAMILTONIAN CORRECTION FROM FILE                                **
!     **************************************************************************
      USE LMTO_DROPPICK_MODULE, ONLY: DHOFK
      IMPLICIT NONE
      INTEGER(4) :: NFIL
      INTEGER(4) :: NKPT,NSPIN,NCORR
      INTEGER(4) :: IKPT,ISPIN,I
      INTEGER(4) :: XSHAPE(4)
!     **************************************************************************
                                     CALL TRACE$PUSH('LMTO_DROPPICK_WRITEDHOFK')
      XSHAPE=SHAPE(DHOFK)
      NCORR=XSHAPE(1)
      NKPT=XSHAPE(3)
      NSPIN=XSHAPE(4)
      CALL FILEHANDLER$UNIT('DHOFK_OUT',NFIL)
      REWIND(NFIL)
      WRITE(NFIL)NKPT,NSPIN,NCORR
      I=0
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          I=I+1
          WRITE(NFIL)I,DHOFK(:,:,IKPT,ISPIN)
        ENDDO
      ENDDO
      CALL FILEHANDLER$CLOSE('DHOFK_OUT')
                                      CALL TRACE$POP()
      RETURN
      END
!
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE LMTO_DROPPICK_WRITE()
!!$!     **************************************************************************
!!$!     **************************************************************************
!!$      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
!!$      USE LMTO_MODULE, ONLY : TDROP,SBAR,DENMAT,LOX,LNX,ISPECIES,POTPAR
!!$      USE LMTO_DROPPICK_MODULE, ONLY : NB1,NB2,NBW,TPRO,NCORR,DHOFK
!!$      IMPLICIT NONE
!!$      INTEGER(4)            :: NFIL
!!$      INTEGER(4)            :: NAT
!!$      INTEGER(4)            :: NKPT
!!$      INTEGER(4)            :: IAT,ISP,IKPT,ISPIN,IPRO,IORB
!!$      INTEGER(4)            :: IWORK16(16)
!!$      REAL(8)   ,ALLOCATABLE:: R0(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: XK(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: WKPT(:)
!!$      REAL(8)   ,ALLOCATABLE:: WKPT(:)
!!$      REAL(8)               :: RBAS(3,3)
!!$      REAL(8)               :: GBAS(3,3)
!!$      REAL(8)               :: SVAR
!!$!     **************************************************************************
!!$      NAT=SIZE(ISPECIES)
!!$      CALL CELL$GETR8A('T0',9,RBAS)
!!$      ALLOCATE(R0(3,NAT))
!!$      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!!$!
!!$!     ==========================================================================
!!$!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!!$!     ==========================================================================
!!$      ALLOCATE(XK(3,NKPTL))
!!$      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
!!$      NKPT=NKPTL
!!$      IF(NKPT.NE.NKPTL) THEN
!!$        CALL ERROR$MSG('ROUTINE NOT SUITED FOR PARALLEL CALCULATIONS')
!!$        CALL ERROR$STOP('LMTO_DROP')
!!$      END IF
!!$      ALLOCATE(WKPT(NKPT))
!!$      CALL DYNOCC$GETR8A('WKPT',NKPT,WKPT)
!!$!     == MULTIPLY WITH SPIN MULTIPLICITY =======================================
!!$      IF(NSPIN.EQ.2) WKPT=2.D0*WKPT
!!$!
!!$!     ==========================================================================
!!$!     == WRITE FILE                                                           ==
!!$!     ==========================================================================
!!$      CALL FILEHANDLER$UNIT('DMFTOUT',NFIL)
!!$      CALL GBASS(RBAS,GBAS,SVAR)
!!$      WRITE(NFIL)NAT,NKPTL,NSPIN,NDIM,NPRO,NCORR,KBTSTAR,RBAS,GBAS !<<<<<<<<<
!!$      DO IAT=1,NAT
!!$        ISP=ISPECIES(IAT)
!!$        IWORK16=-1
!!$        IWORK16(:LNX(ISP))=LOX(:LNX(ISP),ISP)
!!$        WRITE(NFIL,*)IAT,IZ,R0(:,IAT),LNX(ISP),IWORK16 !<<<<<<<<<<<
!!$      ENDDO   
!!$      WRITE(NFIL,*)TPRO  !<<<<<<<
!!$!
!!$      ALLOCATE(PSICORR(NDIM,NCORR,NBW))
!!$      ALLOCATE(AMAT(NCORR,NCORR))
!!$      ALLOCATE(UMAT(NCORR,NCORR))
!!$      ALLOCATE(EMAT(NCORR))
!!$      DO IKPT=1,NKPT
!!$        DO ISPIN=1,NSPIN
!!$!         == CONTRACT WAVE FUNCTION COEFFICIENTS ONTO CORRELATED ORBITALS 
!!$!         == AND CORRELATED BANDS
!!$          I=0
!!$          DO IPRO=1,NPRO
!!$            IF(.NOT.TPRO(IPRO)) CYCLE
!!$            I=I+1
!!$            PSICORR(:,I,:)=PSI(:,IPRO,NB1:NB2)
!!$          ENDDO
!!$!
!!$!         == DEFINE SUM RULE VIOLATION =======================================
!!$          AMAT(:,:)=(0.D0,0.D0)
!!$          DO I=1,NCORR
!!$            DO J=1,NCORR
!!$              CSVAR=0.D0  
!!$              DO IB=1,NBW
!!$                DO IDIM=1,NDIM
!!$                  CSVAR=CSVAR+PSICORR(IDIM,I,IB)*CONJG(PSICORR(IDIM,J,IB))
!!$                ENDDO
!!$              ENDDO
!!$              AMAT(I,J)=AMAT(I,J)+CSVAR
!!$            ENDDO
!!$          ENDDO
!!$!
!!$!         == DETERMINE TRANSFORMATION THAT ENFORCES THE SUM RULE =============
!!$          CALL LIB$DIAGC8(NCORR,AMAT,EMAT,UMAT)
!!$          DO I=1,NCORR
!!$            IF(EMAT(I).LE.0.D0) THEN
!!$              CALL ERROR$MSG('EMAT IS NOT POSITIVE ')
!!$              CALL ERROR$I4VAL('INDEX',I)
!!$              CALL ERROR$R8VAL('EMAT',EMAT(I))
!!$              CALL ERROR$STOP('LMTO_DROP')
!!$            END IF
!!$          ENDDO
!!$!
!!$          UMAT=CONJG(TRANSPOSE(UMAT))
!!$          DO IB=1,NBW
!!$            DO IDIM=1,NDIM
!!$              PSICORR(IDIM,:,IB)=MATMUL(UMAT,PSICORR(IDIM,:,IB))
!!$            ENDDO
!!$          ENDDO
!!$          DO I=1,NCORR
!!$            PSICORR(:,I,:)=PSICORR(:,I,:)/SQRT(EMAT(I))
!!$          ENDDO
!!$          UMAT=CONJG(TRANSPOSE(UMAT))
!!$          DO IB=1,NBW
!!$            DO IDIM=1,NDIM
!!$              PSICORR(IDIM,:,IB)=MATMUL(UMAT,PSICORR(IDIM,:,IB))
!!$            ENDDO
!!$          ENDDO
!!$!
!!$!         ==  CHECK SUM RULE =================================================
!!$          DO I=1,NCORR
!!$            DO J=1,NCORR
!!$              CSVAR=0.D0
!!$              IF(I.EQ.J)CSVAR=(-1.D0,0.D0)
!!$              DO IB=1,NBW
!!$                DO IDIM=1,NDIM
!!$                  CSVAR=CSVAR+PSICORR(IDIM,I,IB)*CONJG(PSICORR(IDIM,J,IB))
!!$                ENDDO
!!$              ENDDO
!!$              IF(ABS(CSVAR).GT.1.D-6) THEN
!!$                CALL ERROR$MSG('SUM RULE IS VIOLATED')
!!$                CALL ERROR$I4VAL('INDEX1',I)
!!$                CALL ERROR$I4VAL('INDEX2',J)
!!$                CALL ERROR$R8VAL('REAL(DEVIATION)',REAL(CSVAR))
!!$                CALL ERROR$R8VAL('IMAG(DEVIATION)',AIMAG(CSVAR))
!!$                CALL ERROR$STOP('LMTO_DROP')
!!$              END IF 
!!$            ENDDO
!!$          ENDDO
!!$          WRITE(NFIL)IKPT,ISPIN,NB,WKPT,XK(:,IKPT) 
!!$          WRITE(NFIL)EIGVAL(:NB),EIGVEC(:NDIM,:NCORR,:NB),DELTAH(:NB,:NB)
!!$        ENDDO
!!$      ENDDO
!!$      CALL FILEHANDLER$CLOSE('DMFTOUT')
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROP()
!     **************************************************************************
!     **  WRITES A FILE WITH KOHN-SHAM WAVE FUNCTIONS AND ENERGIES            **
!     **  FOR USE AS DMFT INTERFACE                                           **
!     **                                                                      **
!     **    <PI_ALPHA|PSI_N> ; E_N                                            **
!     **                                                                      **
!     **  THE BASIS ARE ONSITE-ORTHOGONALIZED NATURAL TIGHT-BINDING ORBITALS  **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : TDROP,SBAR,DENMAT,LOX,LNX,ISPECIES,POTPAR
      USE LMTO_DROPPICK_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: NPRO
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: NFIL,NFIL1
      INTEGER(4)             :: NB,NBH
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      COMPLEX(8),ALLOCATABLE :: PSI(:,:,:)
      COMPLEX(8),ALLOCATABLE :: PSI1(:,:,:)
      REAL(8)   ,ALLOCATABLE :: EIGVAL(:)
      COMPLEX(8),ALLOCATABLE :: EIGVEC(:,:)
      REAL(8)   ,ALLOCATABLE :: OVERLAP(:,:)
      REAL(8)   ,ALLOCATABLE :: OVERLAP1(:,:)
      REAL(8)   ,ALLOCATABLE :: TRANS(:,:)
      REAL(8)   ,ALLOCATABLE :: WKPT(:)
      CHARACTER(16)          :: ID
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NKPT
      INTEGER(4)             :: LMNXT,LMNX
      INTEGER(4)             :: I,J,I1,I2,IDIM,IPRO,IAT,ISP,IKPT,ISPIN,IB,IB1,IBH
      INTEGER(4)             :: IWORK16(16)
      REAL(8)                :: EV
      REAL(8)                :: SVAR
      INTEGER(4)             :: L,M,LN
      COMPLEX(8)             :: CSVAR
      COMPLEX(8),ALLOCATABLE :: AMAT(:,:),UMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: PSICORR(:,:,:)
      REAL(8)   ,ALLOCATABLE :: EMAT(:)
!     **************************************************************************
                                           CALL TRACE$PUSH('LMTO_DROP')
PRINT*,'ENTERING LMTO_DROP'
      IF(.NOT.TDROP) RETURN
      CALL CONSTANTS('EV',EV)
      NAT=SIZE(ISPECIES)
!
!     ==========================================================================
!     ==  DEFINE WINDOW OF BANDS AND CHOICE OF LOCAL ORBITALS                 ==
!     ==========================================================================
      CALL LMTO_DROPICK_INI()
      IF(METHOD.EQ.2) THEN
        ALLOCATE(PSICORR(NDIM,NCORR,NBW))
        ALLOCATE(AMAT(NCORR,NCORR))
        ALLOCATE(UMAT(NCORR,NCORR))
        ALLOCATE(EMAT(NCORR))
      END IF
!
!     ==========================================================================
!     ==  ATTACH FILE                                                         ==
!     ==========================================================================
!      CALL FILEHANDLER$UNIT('DMFTOUT',NFIL)
NFIL=12
OPEN(NFIL,FILE='DMFTOUT.DAT')
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      NKPT=NKPTL
      IF(NKPT.NE.NKPTL) THEN
        CALL ERROR$MSG('ROUTINE NOT SUITED FOR PARALLEL CALCULATIONS')
        CALL ERROR$STOP('LMTO_DROP')
      END IF
      ALLOCATE(WKPT(NKPT))
      CALL DYNOCC$GETR8A('WKPT',NKPT,WKPT)
!     == MULTIPLY WITH SPIN MULTIPLICITY =======================================
      IF(NSPIN.EQ.2) WKPT=2.D0*WKPT
!
!     ==========================================================================
!     ==  CONSTRUCT INDEX ARRAYS                                              ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(IPRO1(NAT))
      ALLOCATE(NPROAT(NAT))
      ID='INFO'
      WRITE(NFIL,*)ID,NAT,NDIM,NSPIN,NKPT
      IPRO=1
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        IPRO1(IAT)=IPRO
        NPROAT(IAT)=SUM(2*LOX(:LNX(ISP),ISP)+1)
!
        IWORK16(:)=-1
        IF(LNX(ISP).GT.16) THEN
          CALL ERROR$MSG('INTERNAL LIMIT FOR VARIABLE LNX EXCEEDED')
          CALL ERROR$STOP('LMTO_DROP')
        END IF
        IWORK16(:LNX(ISP))=LOX(:LNX(ISP),ISP)
        ID='PSIINFO'
        WRITE(NFIL,*)ID,IAT,IPRO1(IAT),NPROAT(IAT),LNX(ISP),IWORK16
!
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
!
      IF(METHOD.EQ.2) THEN
        ID='NORBS'
        WRITE(NFIL,*)ID,NCORR
        ID='ORBINFO'
        I=0
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          DO LN=1,LNX(ISP)
            L=LOX(LN,ISP)
            DO M=1,2*L+1
              I=I+1
              IF(TPRO(I))WRITE(NFIL,*)ID,I,IAT,LN,L,M
            ENDDO
          ENDDO
        ENDDO
      END IF    
!
!     ==========================================================================
!     ==  TRANSFORMATION ONTO ONSITE ORTHOGONALIZED PARTIAL WAVES             ==
!     ==  USES PRE-CALCULATED OVERLAP OF TAILED PARTIAL WAVES                 ==
!     ==========================================================================
      ALLOCATE(T(NAT))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LMNX=NPROAT(IAT)
        T(IAT)%I1=IPRO1(IAT)
        T(IAT)%I2=IPRO1(IAT)-1+NPROAT(IAT)
        ALLOCATE(T(IAT)%MAT(LMNX,LMNX))
        ALLOCATE(T(IAT)%INV(LMNX,LMNX))
        CALL LMTO_ONSORTHO(IAT,LMNX,T(IAT)%MAT,T(IAT)%INV)
      ENDDO
!
!     ==========================================================================
!     ==  WRITE WAVE FUNCTIONS                                                ==
!     ==========================================================================
      NPRO=MAP%NPRO
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          NBH=THIS%NBH
          NB=THIS%NB
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
!
!         ==CHECKS =============================================================
          IF(.NOT.ASSOCIATED(THIS%TBC)) THEN
            CALL ERROR$MSG('THIS%TBC NOT ASSOCIATED')
            CALL ERROR$STOP('LMTO_DROP')
          END IF
          IF(.NOT.ASSOCIATED(THIS%RLAM0)) THEN
            CALL ERROR$MSG('THIS%RLAM0 NOT ASSOCIATED')
            CALL ERROR$STOP('LMTO_DROP')
          END IF
!
!         ==  WRITE WAVE FUNCTIONS =============================================
          ALLOCATE(PSI(NDIM,NPRO,NB))
          ALLOCATE(PSI1(NDIM,NPRO,NB))
          ALLOCATE(EIGVAL(NB))
          ALLOCATE(EIGVEC(NB,NB))
          IF(TINV) THEN
            DO IBH=1,NBH
              PSI(:,:,2*IBH-1)=REAL(THIS%TBC(:,IBH,:))
              PSI(:,:,2*IBH)=AIMAG(THIS%TBC(:,IBH,:))
            ENDDO
          ELSE
            DO IB=1,NB
              PSI(:,:,IB)=THIS%TBC(:,IB,:)
            ENDDO
          END IF
!
!         =====================================================================
!         == TRANSFORM WAVE FUNCTION TO EIGENSTATES ============================
!         =====================================================================
          CALL LIB$DIAGC8(NB,THIS%RLAM0,EIGVAL,EIGVEC)
          PSI1=PSI
          PSI=(0.D0,0.D0)
          DO IB=1,NB
            DO IB1=1,NB
              PSI(:,:,IB)=PSI(:,:,IB)+PSI1(:,:,IB1)*EIGVEC(IB1,IB)
            ENDDO
          ENDDO
!
!         =====================================================================
!         == TRANSFORM TO ONSITE-ORTHOGONAL ORBITALS ===========================
!         =====================================================================
          PSI1=PSI
          DO IAT=1,NAT
            I1=IPRO1(IAT)
            I2=I1-1+NPROAT(IAT)
            DO IB=1,NB
              DO IDIM=1,NDIM
                PSI(IDIM,I1:I2,IB)=MATMUL(T(IAT)%INV,PSI1(IDIM,I1:I2,IB))
              ENDDO
            ENDDO
          ENDDO
!
!         ======================================================================
!         == SPLIT CODE FOR DIFFERENT METHODS                                 ==
!         ======================================================================
!         ======================================================================
!         == PROJECTION ONTO SITE-ORTHONORMALIZED ORBITALS                    ==
!         ======================================================================
          IF(METHOD.EQ.1) THEN
!           == WRITE WAVE FUNCTIONS TO FILE ====================================
            ID='KINFO'
            WRITE(NFIL,*)ID,NDIM,NB,NPRO,WKPT(IKPT)
            ID='PSI'
            DO IB=1,NB
              WRITE(NFIL,*)ID,IB,IKPT,ISPIN,EIGVAL(IB),PSI(:,:,IB)
!PRINT*,'NORM ',IKPT,ISPIN,IB,DOT_PRODUCT(PSI(1,:,IB),PSI(1,:,IB))
            ENDDO
!
!         ======================================================================
!         == MINIMAL METHOD TO ENSURE SUMRULE                                 ==
!         ======================================================================
          ELSE IF(METHOD.EQ.2) THEN
!           == CONTRACT WAVE FUNCTION COEFFICIENTS ONTO CORRELATED ORBITALS 
!           == AND CORRELATED BANDS
!PRINT*,'DOING IKPT=',IKPT,' AND ISPIN=',ISPIN
            I=0
            DO IPRO=1,NPRO
              IF(.NOT.TPRO(IPRO)) CYCLE
              I=I+1
              PSICORR(:,I,:)=PSI(:,IPRO,NB1:NB2)
            ENDDO
!!$DO IB=1,NBW
!!$  PRINT*,'PSICORR ',IB,PSICORR(:,:,IB)
!!$ENDDO
!
!           == DEFINE SUM RULE VIOLATION =======================================
            AMAT(:,:)=(0.D0,0.D0)
            DO I=1,NCORR
              DO J=1,NCORR
                CSVAR=0.D0  
                DO IB=1,NBW
                  DO IDIM=1,NDIM
                    CSVAR=CSVAR+PSICORR(IDIM,I,IB)*CONJG(PSICORR(IDIM,J,IB))
                  ENDDO
                ENDDO
                AMAT(I,J)=AMAT(I,J)+CSVAR
              ENDDO
            ENDDO
!
!           == DETERMINE TRANSFORMATION THAT ENFORCES THE SUM RULE =============
            CALL LIB$DIAGC8(NCORR,AMAT,EMAT,UMAT)
!!$PRINT*,'EMAT ',EMAT 
           DO I=1,NCORR
              IF(EMAT(I).LE.0.D0) THEN
                CALL ERROR$MSG('EMAT IS NOT POSITIVE ')
                CALL ERROR$I4VAL('INDEX',I)
                CALL ERROR$R8VAL('EMAT',EMAT(I))
                CALL ERROR$STOP('LMTO_DROP')
              END IF
            ENDDO
!
            UMAT=CONJG(TRANSPOSE(UMAT))
            DO IB=1,NBW
              DO IDIM=1,NDIM
                PSICORR(IDIM,:,IB)=MATMUL(UMAT,PSICORR(IDIM,:,IB))
              ENDDO
            ENDDO
            DO I=1,NCORR
              PSICORR(:,I,:)=PSICORR(:,I,:)/SQRT(EMAT(I))
            ENDDO
            UMAT=CONJG(TRANSPOSE(UMAT))
            DO IB=1,NBW
              DO IDIM=1,NDIM
                PSICORR(IDIM,:,IB)=MATMUL(UMAT,PSICORR(IDIM,:,IB))
              ENDDO
            ENDDO
!
!           ==  CHECK SUM RULE =================================================
            DO I=1,NCORR
              DO J=1,NCORR
                CSVAR=0.D0
                IF(I.EQ.J)CSVAR=(-1.D0,0.D0)
                DO IB=1,NBW
                  DO IDIM=1,NDIM
                    CSVAR=CSVAR+PSICORR(IDIM,I,IB)*CONJG(PSICORR(IDIM,J,IB))
                  ENDDO
                ENDDO
                IF(ABS(CSVAR).GT.1.D-6) THEN
                  CALL ERROR$MSG('SUM RULE IS VIOLATED')
                  CALL ERROR$I4VAL('INDEX1',I)
                  CALL ERROR$I4VAL('INDEX2',J)
                  CALL ERROR$R8VAL('REAL(DEVIATION)',REAL(CSVAR))
                  CALL ERROR$R8VAL('IMAG(DEVIATION)',AIMAG(CSVAR))
                  CALL ERROR$STOP('LMTO_DROP')
                END IF 
              ENDDO
            ENDDO
!
!           == WRITE WAVE FUNCTION TO FILE =====================================
            ID='KINFO'
            WRITE(NFIL,*)ID,NDIM,NBW,NCORR,WKPT(IKPT)
            ID='PSI'
            DO IB=1,NBW
              WRITE(NFIL,*)ID,IB,IKPT,ISPIN,EIGVAL(IB+NB1-1) &
       &                  ,((PSICORR(J,I,IB),J=1,NDIM),I=1,NCORR)
            ENDDO
          END IF
          DEALLOCATE(EIGVAL)
          DEALLOCATE(EIGVEC)
          DEALLOCATE(PSI)
          DEALLOCATE(PSI1)
        ENDDO
      ENDDO
      IF(METHOD.EQ.2) THEN
        DEALLOCATE(PSICORR)
        DEALLOCATE(EMAT)
        DEALLOCATE(AMAT)
        DEALLOCATE(UMAT)
      END IF
      DO IAT=1,NAT
        DEALLOCATE(T(IAT)%MAT)
        DEALLOCATE(T(IAT)%INV)
      ENDDO
      DEALLOCATE(T)
      CLOSE(NFIL)
!
!PRINT*,'TEST READING '
!  CALL LMTO_DROPREAD()
      CALL ERROR$MSG('FORCED STOP AFTER WRITING FILE')
      CALL ERROR$STOP('LMTO$DROP')
                                           CALL TRACE$POP()
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPREAD()
!     **************************************************************************
!     ** READS THE WAVE FUNCTIONS IN LOCAL ORBITALS, WHICH HAVE BEEN          **
!     ** ONSITE-ORTHONORMALIZED.                                              **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(16)          :: ID
      INTEGER(4)             :: NAT       ! #(ATOMS)
      INTEGER(4)             :: NDIM      ! #(SPINOR-COMPONENTS)
      INTEGER(4)             :: NSPIN     ! #(SPIN DIRECTIONS)
      INTEGER(4)             :: NKPT      ! #(K-POINTS)
      INTEGER(4)             :: NPRO      ! TOTAL #(ORBITALS)
      INTEGER(4)             :: NB        ! #(STATES)
      INTEGER(4),ALLOCATABLE :: IPRO1(:)  ! FIRST ORBITAL INDEX FOR THIS ATOM
      INTEGER(4),ALLOCATABLE :: NPROAT(:) ! #(ORBITALS) FOR THIS ATOM
      INTEGER(4),ALLOCATABLE :: LNX(:)    ! #(DIFFERENT ORBITAL SHELLS)
      INTEGER(4),ALLOCATABLE :: LOX(:,:)  ! ANGULAR MOMENTUM OF THIS ORBITAL SHELL
      REAL(8)   ,ALLOCATABLE :: E(:)      ! ENERGY EIGENVALUES
      COMPLEX(8),ALLOCATABLE :: C(:,:,:)  ! ORBITAL COEFFICIENTS
      INTEGER(4)             :: IAT,IB,IKPT,ISPIN
      INTEGER(4)             :: IAT1,IB1,IKPT1,ISPIN1
      INTEGER(4)             :: NFIL  ! FORTRAN FILE UNIT
!     **************************************************************************
      NFIL=12
      OPEN(NFIL,FILE='NTBOWV.DATA')
      REWIND(NFIL)
!
      READ(NFIL,*)ID,NAT,NDIM,NSPIN,NKPT
      IF(ID.NE.'INFO') THEN
        STOP 'INCORRECT ID_1'
      END IF
      ALLOCATE(IPRO1(NAT))
      ALLOCATE(NPROAT(NAT))
      ALLOCATE(LNX(NAT))
      ALLOCATE(LOX(16,NAT))
      DO IAT=1,NAT
        READ(NFIL,*)ID,IAT1,IPRO1(IAT),NPROAT(IAT),LNX(IAT),LOX(:,IAT)
        IF(ID.NE.'PSIINFO') THEN
          STOP 'INCORRECT ID_2'
        END IF
        IF(IAT1.NE.IAT) THEN
          STOP 'ATOM INDEX OUT OF ORDER'
        END IF
      ENDDO
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          READ(NFIL,*)ID,NDIM,NB,NPRO
          IF(ID.NE.'KINFO') THEN
            STOP 'INCORRECT ID_3'
          END IF
          ALLOCATE(E(NB))
          ALLOCATE(C(NDIM,NPRO,NB))
          DO IB=1,NB
            READ(NFIL,*)ID,IB1,IKPT1,ISPIN1,E(IB),C(:,:,IB)
            IF(ID.NE.'PSI') THEN
              PRINT*,'ID   =',TRIM(ID)
              PRINT*,'IB   =',IB1,IB
              PRINT*,'IKPT =',IKPT1,IKPT
              PRINT*,'ISPIN=',ISPIN1,ISPIN
              STOP 'INCORRECT ID_4'
            END IF
          ENDDO
!
!         ------> GRAB WAVE FUNCTIONS AND ENERGIES FROM HERE <-----
! 
          DEALLOCATE(E)
          DEALLOCATE(C)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PICKGET()
!     **************************************************************************
!     **  READS A FILE HAMILTON CORRECTION                                    **
!     **  FOR USE AS DMFT INTERFACE                                           **
!     **                                                                      **
!     **    <PI_ALPHA|PSI_N> ; E_N                                            **
!     **                                                                      **
!     **  THE BASIS ARE ONSITE-ORTHOGONALIZED NATURAL TIGHT-BINDING ORBITALS  **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : TPICK,ISPECIES,LNX,LOX
      USE LMTO_DROPPICK_MODULE, ONLY: TPICKED,NB1,NB2,NBW,NCORR,T,DHOFK
      USE STRINGS_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TWRITEDHOFK=.FALSE.
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      LOGICAL(4),SAVE        :: TCONNECTEDFILE=.FALSE.
      INTEGER(4)             :: NFIL
      INTEGER(4)             :: NB,NBH
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: NKPT1
      INTEGER(4)             :: NAT
      INTEGER(4)             :: LMNX
      REAL(8)                :: EFERMI
      COMPLEX(8),ALLOCATABLE :: DH(:,:)
      REAL(8)                :: BETA ! 1/(K_B*T)
      REAL(8)                :: EV
      COMPLEX(8)             :: A,B
      COMPLEX(8),ALLOCATABLE :: PIPSI(:,:,:)
      INTEGER(4)             :: IAT,ISP
      INTEGER(4)             :: IB1M,IB1P,IB2M,IB2P
      INTEGER(4)             :: IB,IBH,IB1,IB2,IDIM,IKPT,ISPIN,IK,I,J
      LOGICAL(4)             :: TINV
CHARACTER(128) ::STRING
!     **************************************************************************
                                           CALL TRACE$PUSH('LMTO_LMTO_PICKGET')
      IF(.NOT.TPICK) RETURN
      IF(TPICKED) RETURN
      TPICKED=.TRUE.
      CALL CONSTANTS('EV',EV)
      CALL LMTO_DROPICK_INI()
!
!     ==========================================================================
!     ==  TRANSFORMATION ONTO ONSITE ORTHOGONALIZED PARTIAL WAVES             ==
!     ==  USES PRE-CALCULATED OVERLAP OF TAILED PARTIAL WAVES                 ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(T(NAT))
      NPRO=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LMNX=SUM(2*LOX(:LNX(ISP),ISP)+1)
        T(IAT)%I1=NPRO+1
        T(IAT)%I2=NPRO+LMNX
        ALLOCATE(T(IAT)%MAT(LMNX,LMNX))
        ALLOCATE(T(IAT)%INV(LMNX,LMNX))
        CALL LMTO_ONSORTHO(IAT,LMNX,T(IAT)%MAT,T(IAT)%INV)
        NPRO=NPRO+LMNX
      ENDDO
!
!     ==========================================================================
!     ==  READ FILE HEADER                                                    ==
!     ==========================================================================
      CALL FILEHANDLER$UNIT('DMFTIN',NFIL)
!!$CALL FILEHANDLER$REPORT(6,'ALL')
!!$PRINT*,'NFIL FOR DMFTIN ',NFIL
!!$INQUIRE(NFIL,FILE=STRING)
!!$PRINT*,'FILE=',TRIM(STRING)
!!$INQUIRE(NFIL,ACTION=STRING)
!!$PRINT*,'ACTION=',TRIM(STRING)
!!$INQUIRE(NFIL,POSITION=STRING)
!!$PRINT*,'POSITION=',TRIM(STRING)
      REWIND NFIL
      READ(NFIL,*)EFERMI    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      EFERMI=EFERMI*EV
      READ(NFIL,*)BETA      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      BETA=BETA/EV
      READ(NFIL,*)NKPT1 !# (KPOINTS)   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF(NKPT1.NE.NKPTL) THEN
        CALL ERROR$MSG('INCONSISTENT NUMBER OF K-POINTS')
        CALL ERROR$I4VAL('INTERNAL NUMBER OF KPOINTS',NKPTL)
        CALL ERROR$I4VAL('NUMBER OF KPOINTS ON FILE',NKPT1)
        CALL ERROR$STOP('LMTO_PICKGET')
      END IF
      IF(NSPIN.NE.1) THEN
        CALL ERROR$MSG('DMFT INTERFACE IS LIMITED TO NSPIN=1')
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('LMTO_PICKGET')
      END IF
      IF(NDIM.NE.1) THEN
        CALL ERROR$MSG('DMFT INTERFACE IS LIMITED TO NDIM=1')
        CALL ERROR$I4VAL('NDIM',NDIM)
        CALL ERROR$STOP('LMTO_PICKGET')
      END IF
!
PRINT*,'LMTO_PICKET MARKE 4'
      IF(.NOT.ASSOCIATED(DHOFK)) ALLOCATE(DHOFK(NCORR,NCORR,NKPTL,NSPIN))
      ALLOCATE(DH(NBW,NBW))
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          NBH=THIS%NBH
          NB=THIS%NB
!
!         ======================================================================
!         ==  GET HAMILTONIAN CORRECTION FROM FILE                            ==
!         ======================================================================
          READ(NFIL,*)IK,I    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          IF(NBW.NE.I) THEN
            CALL ERROR$MSG('INCONSISTENT SIZE OF HAMILTONIAN')
            CALL ERROR$I4VAL('INTERNAL NUMBER OF CORRELATED BANDS',NBW)
            CALL ERROR$I4VAL('NUMBER OF BANDS ON FILE',I)
            CALL ERROR$STOP('LMTO_LMTO_PICKGET')
          END IF
          DO IB=1,NBW
            READ(NFIL,*)DH(IB,:)   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          ENDDO
          DH=DH*EV
!
!         ======================================================================
!         == TRANSFORM HAMILTONIAN TO ONSITE ORTHOGONALIZED ORBITALS          ==
!         ======================================================================
          ALLOCATE(PIPSI(NDIM,NPRO,NB))
          IF(TINV) THEN
            DO IBH=1,NBH
              PIPSI(:,:,2*IBH-1)=REAL(THIS%TBC(:,IBH,:))
              PIPSI(:,:,2*IBH)=AIMAG(THIS%TBC(:,IBH,:))
            ENDDO
          ELSE
            DO IB=1,NB
              PIPSI(:,:,IB)=THIS%TBC(:,IB,:)
            ENDDO
          END IF
          CALL LMTO_DROPPICK_TRANSFORM(NAT,T,NDIM,NPRO,NB,THIS%RLAM0,PIPSI &
      &                               ,DH,DHOFK(:,:,IKPT,ISPIN))
          DEALLOCATE(PIPSI)
        ENDDO
      ENDDO
      DEALLOCATE(DH)
      CALL FILEHANDLER$CLOSE('DMFTIN')
!
!     ==========================================================================
!     == WRITE DHOFK TO FILE                                                  ==
!     ==========================================================================
      IF(TWRITEDHOFK) THEN
        CALL FILEHANDLER$UNIT('DHOFK_OUT',NFIL)
        REWIND NFIL
        WRITE(NFIL)NKPTL,NSPIN,NCORR
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            WRITE(NFIL)IKPT,DHOFK(:,:,IKPT,ISPIN)
          ENDDO
        ENDDO
        CALL FILEHANDLER$CLOSE('DHOFK_OUT')
      END IF
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DROPPICK_TRANSFORM(NAT,T,NDIM,NPRO,NB,LAMBDA,PIPSI,DH,DH1)
!     **************************************************************************
!     ** CALCULATES THE DMFT CORRECTION IN A BASIS OF SITE ORTHOGONALIZED     **
!     ** ORBITALS WHICH IS CONTRACTED ON A CORRELATED SUBSET                  **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2012******************
      USE LMTO_DROPPICK_MODULE,ONLY : NCORR,NB1,NB2,NBW,TPRO,T_TYPE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NAT
      TYPE(T_TYPE),INTENT(IN) :: T(NAT)
      INTEGER(4)  ,INTENT(IN) :: NDIM
      INTEGER(4)  ,INTENT(IN) :: NPRO
      INTEGER(4)  ,INTENT(IN) :: NB
      COMPLEX(8)  ,INTENT(IN) :: LAMBDA(NB,NB)
      COMPLEX(8)  ,INTENT(IN) :: PIPSI(NDIM,NPRO,NB)
      COMPLEX(8)  ,INTENT(IN) :: DH(NBW,NBW)
      COMPLEX(8)  ,INTENT(OUT):: DH1(NCORR,NCORR)
      COMPLEX(8)              :: EIGVEC(NB,NB)
      REAL(8)                 :: EIGVAL(NB)
      COMPLEX(8)              :: PSICORR(NDIM,NCORR,NBW)
      COMPLEX(8)              :: AMAT(NCORR,NCORR)
      COMPLEX(8)              :: UMAT(NCORR,NCORR)
      REAL(8)                 :: EMAT(NCORR)
      COMPLEX(8)              :: CSVAR
      COMPLEX(8)              :: PIPSI1(NDIM,NPRO,NBW)
      INTEGER(4)              :: I,J,IPRO,IB,IDIM,I1,I2,IB1,IB2,IAT
!     **************************************************************************
!
!     ==========================================================================
!     ==  MAKE LAMBDA MATRIX SYMMETRIC
!     ==  AND TRANSFORM WAVE FUNCTION PROJECTIONS TO EIGENSTATES              ==
!     ==  AND PROJECT ON THE SUBSET OF BANDS                                  ==
!     ==========================================================================
DO IB1=1,NB
  DO IB2=IB1,NB
    CSVAR=LAMBDA(IB1,IB2)-CONJG(LAMBDA(IB2,IB1))
    CSVAR=(LAMBDA(IB2,IB1)+CONJG(LAMBDA(IB2,IB1)))/2.D0
!    IF(ABS(CSVAR).GT.1.D-6) THEN
      WRITE(*,*)IB1,IB2,LAMBDA(IB1,IB2),LAMBDA(IB2,IB1)
!    END IF
  ENDDO
ENDDO
STOP 'FORCED STOP IN LMTO_DROPPICK_TRANSFORM'
!
!     ==========================================================================
!     ==  DIAGONALIZE LAGRANGE MULTIPLIERS TO OBTAIN EIGENSTATES              ==
!     ==  AND TRANSFORM WAVE FUNCTION PROJECTIONS TO EIGENSTATES              ==
!     ==  AND PROJECT ON THE SUBSET OF BANDS                                  ==
!     ==========================================================================
      CALL LIB$DIAGC8(NB,LAMBDA,EIGVAL,EIGVEC)
      PIPSI1=PIPSI
      PIPSI1=(0.D0,0.D0)
      DO IB1=NB1,NB2
        IB=IB1-NB1+1
        DO IB2=1,NB
          PIPSI1(:,:,IB)=PIPSI1(:,:,IB)+PIPSI(:,:,IB2)*EIGVEC(IB2,IB1)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == MAP ONTO SITE ORTHOGONALIZED ORBITALS                                ==
!     ==========================================================================
      DO IAT=1,NAT
        I1=T(IAT)%I1
        I2=T(IAT)%I2
        DO IB=1,NBW
          DO IDIM=1,NDIM
            PIPSI1(IDIM,I1:I2,IB)=MATMUL(T(IAT)%INV,PIPSI1(IDIM,I1:I2,IB))
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == CONTRACT WAVE FUNCTION COEFFICIENTS ONTO CORRELATED ORBITALS         ==
!     ==========================================================================
      I=0
      DO IPRO=1,NPRO
        IF(.NOT.TPRO(IPRO)) CYCLE
        I=I+1
        PSICORR(:,I,:)=PIPSI1(:,IPRO,:)
       ENDDO
!
!     ==========================================================================
!     == ENFORCE SUM RULE                                                     ==
!     ==========================================================================
!     == DEFINE SUM RULE VIOLATION =======================================
      AMAT(:,:)=(0.D0,0.D0)
      DO I=1,NCORR
        DO J=1,NCORR
          CSVAR=0.D0  
          DO IB=1,NBW
            DO IDIM=1,NDIM
              CSVAR=CSVAR+PSICORR(IDIM,I,IB)*CONJG(PSICORR(IDIM,J,IB))
            ENDDO
          ENDDO
          AMAT(I,J)=AMAT(I,J)+CSVAR
        ENDDO
      ENDDO
!
!     == DETERMINE TRANSFORMATION THAT ENFORCES THE SUM RULE =============
      CALL LIB$DIAGC8(NCORR,AMAT,EMAT,UMAT)
      DO I=1,NCORR
        IF(EMAT(I).LE.0.D0) THEN
          CALL ERROR$MSG('EMAT IS NOT POSITIVE ')
          CALL ERROR$I4VAL('INDEX',I)
          CALL ERROR$R8VAL('EMAT',EMAT(I))
          CALL ERROR$STOP('LMTO_DROPPICK_TRANSFORM')
        END IF
      ENDDO
!
      UMAT=CONJG(TRANSPOSE(UMAT))
      DO IB=1,NBW
        DO IDIM=1,NDIM
          PSICORR(IDIM,:,IB)=MATMUL(UMAT,PSICORR(IDIM,:,IB))
        ENDDO
      ENDDO
      DO I=1,NCORR
        PSICORR(:,I,:)=PSICORR(:,I,:)/SQRT(EMAT(I))
      ENDDO
      UMAT=CONJG(TRANSPOSE(UMAT))
      DO IB=1,NBW
        DO IDIM=1,NDIM
          PSICORR(IDIM,:,IB)=MATMUL(UMAT,PSICORR(IDIM,:,IB))
        ENDDO
      ENDDO
!
!     ==  CHECK SUM RULE =================================================
      DO I=1,NCORR
        DO J=1,NCORR
          CSVAR=0.D0
          IF(I.EQ.J)CSVAR=(-1.D0,0.D0)
          DO IB=1,NBW
            DO IDIM=1,NDIM
              CSVAR=CSVAR+PSICORR(IDIM,I,IB)*CONJG(PSICORR(IDIM,J,IB))
            ENDDO
          ENDDO
          IF(ABS(CSVAR).GT.1.D-6) THEN
            CALL ERROR$MSG('SUM RULE IS VIOLATED')
            CALL ERROR$I4VAL('INDEX1',I)
            CALL ERROR$I4VAL('INDEX2',J)
            CALL ERROR$R8VAL('REAL(DEVIATION)',REAL(CSVAR))
            CALL ERROR$R8VAL('IMAG(DEVIATION)',AIMAG(CSVAR))
            CALL ERROR$STOP('LMTO_DROP')
          END IF 
        ENDDO
      ENDDO
!     === NOW WE HAVE <\PIBAR|PSI> =============================================
!
!     ==========================================================================
!     == MULTIPLICATION WITH SQRT(Q) TO GET TO ONSITE ORTHOGONALIZED ORBITALS ==
!     ==========================================================================
      UMAT=CONJG(TRANSPOSE(UMAT))
      DO IB=1,NBW
        DO IDIM=1,NDIM
          PSICORR(IDIM,:,IB)=MATMUL(UMAT,PSICORR(IDIM,:,IB))
        ENDDO
      ENDDO
      DO I=1,NCORR
        PSICORR(:,I,:)=PSICORR(:,I,:)/SQRT(EMAT(I))
      ENDDO
      UMAT=CONJG(TRANSPOSE(UMAT))
      DO IB=1,NBW
        DO IDIM=1,NDIM
          PSICORR(IDIM,:,IB)=MATMUL(UMAT,PSICORR(IDIM,:,IB))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == EXPAND WAVE FUNCTION COEFFICIENTS FROM CORRELATED ORBITALS           ==
!     ==========================================================================
      IF(NDIM.NE.1) THEN
        CALL ERROR$MSG('IMPLEMENTATION IS RESTRICTED TO NDIM=1')
        CALL ERROR$STOP('LMTO_DROPPICK_TRANSFORM')
      END IF
      DH1(:,:)=MATMUL(PSICORR(1,:,:) &
     &         ,MATMUL(DH,TRANSPOSE(CONJG(PSICORR(1,:,:)))))
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PICK2()
!     **************************************************************************
!     **  READS A FILE HAMILTON CORRECTION                                    **
!     **  FOR USE AS DMFT INTERFACE                                           **
!     **                                                                      **
!     **    <PI_ALPHA|PSI_N> ; E_N                                            **
!     **                                                                      **
!     **  THE BASIS ARE ONSITE-ORTHOGONALIZED NATURAL TIGHT-BINDING ORBITALS  **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : TPICK,SBAR,DENMAT,LOX,LNX,ISPECIES,POTPAR
      USE STRINGS_MODULE
      IMPLICIT NONE
      TYPE T_TYPE
        REAL(8),POINTER :: MAT(:,:)
        REAL(8),POINTER :: INV(:,:)
      END TYPE T_TYPE
      INTEGER(4)             :: NFIL
      INTEGER(4)             :: NB,NBH
      INTEGER(4)             :: NPRO
      INTEGER(4),ALLOCATABLE :: IPRO1(:)
      INTEGER(4),ALLOCATABLE :: NPROAT(:)
      TYPE(T_TYPE),ALLOCATABLE :: T(:)
      CHARACTER(16)          :: ID
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NKPT
      INTEGER(4)             :: LMNXT,LMNX
      INTEGER(4)             :: I,J,I1,I2,IDIM,IPRO,IAT,ISP,IKPT,ISPIN,IBH
      INTEGER(4)             :: NPRO1,NKPT1,NCHI,IK
      REAL(8)                :: EFERMI
      COMPLEX(8),ALLOCATABLE :: DH(:,:)
      REAL(8)                :: BETA ! 1/(K_B*T)
      REAL(8)                :: EV
      LOGICAL(4),SAVE        :: TCONNECTEDFILE
!     **************************************************************************
                                           CALL TRACE$PUSH('LMTO_NTBODENMAT')
      IF(.NOT.TPICK) RETURN
      CALL CONSTANTS('EV',EV)
!
!     ==========================================================================
!     ==  ATTACH FILE                                                         ==
!     ==========================================================================
      IF(.NOT.TCONNECTEDFILE) THEN
PRINT*,'MARKE 1: CONNECTING FILE'
        CALL FILEHANDLER$SETFILE('DMFTIN',.FALSE.,-'DMFTIN.DAT')
!        CALL FILEHANDLER$SETSPECIFICATION('DMFTIN','STATUS','NEW')
!        CALL FILEHANDLER$SETSPECIFICATION('DMFTIN','ACTION','READ')
        CALL FILEHANDLER$SETSPECIFICATION('DMFTIN','FORM','FORMATTED')
        TCONNECTEDFILE=.TRUE.
PRINT*,'MARKE 2:  FILE CONNECTED'
      END IF
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      NKPT=NKPTL
      IF(NKPT.NE.NKPTL) THEN
        CALL ERROR$MSG('ROUTINE NOT SUITED FOR PARALLEL CALCULATIONS')
        CALL ERROR$STOP('LMTO_DROP')
      END IF
!
!     ==========================================================================
!     ==  CONSTRUCT INDEX ARRAYS                                              ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(IPRO1(NAT))
      ALLOCATE(NPROAT(NAT))
      IPRO=1
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        IPRO1(IAT)=IPRO
        NPROAT(IAT)=SUM(2*LOX(:LNX(ISP),ISP)+1)
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
      NPRO=IPRO-1
!
!     ==========================================================================
!     ==  TRANSFORMATION ONTO ONSITE ORTHOGONALIZED PARTIAL WAVES             ==
!     ==  USES PRE-CALCULATED OVERLAP OF TAILED PARTIAL WAVES                 ==
!     ==========================================================================
      ALLOCATE(T(NAT))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LMNX=NPROAT(IAT)
        ALLOCATE(T(IAT)%MAT(LMNX,LMNX))
        ALLOCATE(T(IAT)%INV(LMNX,LMNX))
        CALL LMTO_ONSORTHO(IAT,LMNX,T(IAT)%MAT,T(IAT)%INV)
      ENDDO
!
!     ==========================================================================
!     ==  READ FILE HEADER                                                    ==
!     ==========================================================================
PRINT*,'MARKE 3:  CALL UNIT'
      CALL FILEHANDLER$UNIT('DMFTIN',NFIL)
PRINT*,'MARKE 4:  UNIT CALLED ',NFIL
      READ(NFIL,*)EFERMI
PRINT*,'MARKE 5:  EFERMI READ ',EFERMI
      EFERMI=EFERMI*EV
      READ(NFIL,*)BETA
PRINT*,'MARKE 6:  BETA READ ',BETA
      BETA=BETA/EV
      READ(NFIL,*)NKPT1 !# (KPOINTS)
PRINT*,'MARKE 6:  NKPT READ ',NKPT,NKPT1
      IF(NKPT1.NE.NKPT) THEN
        CALL ERROR$MSG('INCONSISTENT NUMBER OF K-POINTS')
        CALL ERROR$I4VAL('INTERNAL NUMBER OF KPOINTS',NKPT)
        CALL ERROR$I4VAL('NUMBER OF KPOINTS ON FILE',NKPT1)
        CALL ERROR$STOP('LMTO_PICK')
      END IF
      IF(NSPIN.NE.1) THEN
        CALL ERROR$MSG('DMFT INTERFACE IS LIMITED TO NSPIN=1')
        CALL ERROR$I4VAL('NSPIN',NSPIN)
        CALL ERROR$STOP('LMTO_PICK')
      END IF
      IF(NDIM.NE.1) THEN
        CALL ERROR$MSG('DMFT INTERFACE IS LIMITED TO NDIM=1')
        CALL ERROR$I4VAL('NDIM',NDIM)
        CALL ERROR$STOP('LMTO_PICK')
      END IF
!
      NPRO=MAP%NPRO
      ALLOCATE(DH(NPRO,NPRO))
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          NBH=THIS%NBH
          NB=THIS%NB
          CALL PLANEWAVE$SELECT(GSET%ID)
!
!         ==GET HAMILTONIAN CORRECTION FROM FILE================================
PRINT*,'READING IK,NPRO1 '
          READ(NFIL,*)IK,NB
PRINT*,'IK,NPRO1 READ: ',IK,NPRO1
          IF(NPRO1.NE.NPRO) THEN
            CALL ERROR$MSG('INCONSISTENT SIZE OF HAMILTONIAN')
            CALL ERROR$I4VAL('INTERNAL BASIS-SET SIZE',NPRO)
            CALL ERROR$I4VAL('BASIS-SET SIZE ON FILE',NPRO1)
            CALL ERROR$STOP('LMTO_READDELTAH')
          END IF
PRINT*,'READING DH ',IKPT,NPRO,NKPT,NB
          DO IPRO=1,NPRO
            READ(NFIL,*)DH(:,IPRO)
          ENDDO
PRINT*,'DH READ '
!
!         == TRANSFORM FROM ONSITE-ORTHOGONAL ORBITALS =========================
          DO IAT=1,NAT
            I1=IPRO1(IAT)
            I2=I1-1+NPROAT(IAT)
            DO IPRO=1,NPRO
              DH(I1:I2,IPRO)=MATMUL(TRANSPOSE(T(IAT)%INV),DH(I1:I2,IPRO))
              DH(IPRO,I1:I2)=MATMUL(DH(IPRO,I1:I2),T(IAT)%INV)
            ENDDO
          ENDDO
!
!         == MAP ONTO HTBC =====================================================
          IF(.NOT.ASSOCIATED(THIS%HTBC))ALLOCATE(THIS%HTBC(NDIM,NBH,NPRO))
          IF(.NOT.ASSOCIATED(THIS%TBC)) THEN
            CALL ERROR$MSG('THIS%TBC NOT ASSOCIATED')
            CALL ERROR$STOP('LMTO_PICK')
          END IF
          DO IBH=1,NBH
            THIS%HTBC(1,IBH,:)=MATMUL(DH,THIS%TBC(1,IBH,:))
          ENDDO
        ENDDO
      ENDDO
      DO IAT=1,NAT
        DEALLOCATE(T(IAT)%MAT)
        DEALLOCATE(T(IAT)%INV)
      ENDDO
      DEALLOCATE(T)
      ALLOCATE(DH(NPRO,NPRO))
      CLOSE(NFIL)
      CALL ERROR$MSG('FORCED STOP AFTER READING FILE')
      CALL ERROR$STOP('LMTO$PICK')
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$ETOT(LMNXX_,NDIMD_,NAT_,DENMAT_)
      USE LMTO_MODULE, ONLY : TON
!     **************************************************************************
!     **                                                                      **
!     **  DENMAT_ ON INPUT IS CALCULATED DIRECTLY FROM THE PROJECTIONS AND    **
!     **  IS USED IN THE AUGMENTATION                                         **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
USE LMTO_MODULE, ONLY : TDROP,TPICK,DENMAT,HAMIL,TOFFSITE,THTBC
USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMNXX_
      INTEGER(4),INTENT(IN) :: NDIMD_
      INTEGER(4),INTENT(IN) :: NAT_
      COMPLEX(8),INTENT(IN) :: DENMAT_(LMNXX_,LMNXX_,NDIMD_,NAT_)
      INTEGER(4)            :: SWITCH
INTEGER(4) ::IX,NN,IND1,IND2,IND3
REAL(8)    :: XDELTA,XSVAR,XENERGY
!     **************************************************************************
      IF(.NOT.TON) RETURN
      WRITE(*,FMT='(82("="),T30," LMTO$ENERGY START ")')
      THTBC=.FALSE.   
      CALL LMTO$SETHTBCTOZERO()
!
!     ==========================================================================
!     == WRITE DMFT INTERFACE 
!     ==========================================================================
      IF(TDROP.OR.TPICK) THEN
        IF(TPICK) THEN
          PRINT*,'CALLING DMFT INTERFACE PICK ....'
          CALL LMTO_DROPPICK_HTBC()
          PRINT*,'.... DMFT INTERFACE PICK DONE'
        END IF
        IF(TDROP) THEN
          CALL LMTO_DROPPICK_DROP()   !OLD:   CALL LMTO_DROP()
          CALL ERROR$MSG('REGULAR STOP AFTER EXECUTING LMTO_DROPPICK_DROP')
          CALL ERROR$MSG('DROP IS EXEWCUTED ONLY ONCE')
          CALL ERROR$STOP('LMTO$ETOT')
        ELSE  
          RETURN  ! IF DROP OR PICK IS TRUE INTERFACE DMFT IS USED.
        END IF
      END IF
!
!     ==========================================================================
!     ==  TEST
!     ==========================================================================
!      CALL LMTO_TESTORBITALS()
!
!     ==========================================================================
!     ==  CONSTRUCT DENSITY MATRIX
!     ==========================================================================
      CALL TIMING$CLOCKON('NTBODENMAT')
      CALL LMTO_NTBODENMAT()
      IF(TOFFSITE) CALL LMTO_OVERLAPEVAL()
      CALL TIMING$CLOCKOFF('NTBODENMAT')
!!$      CALL LMTO_TESTDENMAT_1CDENMAT(LMNXX_,NDIMD_,NAT_,DENMAT_)
!!$      CALL LMTO_TESTDENMAT()
!!$STOP 'FORCED'
!
!     ==========================================================================
!     ==  WRITE SOME INFO
!     ==========================================================================
!!$      CALL LMTO_PLOTNTBO('TAILED,CUBE',1,1)   !(TYPE,IAT,LMN)
!!$      CALL LMTO_PLOTNTBO('TAILED,STAR',1,1)   !(TYPE,IAT,LMN)
!!$      CALL LMTO_PLOTNTBO('FULL,CUBE',1,1)   !(TYPE,IAT,LMN)
!!$      CALL LMTO_PLOTNTBO('FULL,STAR',1,1)   !(TYPE,IAT,LMN)
!!$STOP 'FORCED AFTER PLOTNTBO'

!!$      CALL LMTO$REPORTPOTBAR(6)
!!$      CALL LMTO$REPORTSBAR(6)

!!$PRINT*,'FUDGE WARNING!!!!! DENSITY MATRIX OVERWRITTEN FOR H2 TEST'
!!$DO NN=1,SIZE(DENMAT)
!!$  DENMAT(NN)%MAT=0.D0
!!$  DENMAT(NN)%MAT(1,1,1)=0.60266D0
!!$ENDDO

!!$      CALL LMTO$REPORTOVERLAP(6)
!STOP 'FORCED'
!
!     ==========================================================================
!     ==  SOME INFO                                                           ==
!     ==========================================================================
!!$      CALL LMTO$REPORTORTHODENMAT(6)
!!$      CALL LMTO$REPORTDENMAT(6)
!!$STOP 'FORCED'
!
!!$NN=2
!!$IND1=1
!!$IND2=2
!!$IND3=1
!!$XSVAR=DENMAT(NN)%MAT(IND1,IND2,IND3)
!!$XDELTA=1.D-2
!!$IX=-3
!!$1000 CONTINUE
!!$IX=IX+1
!!$DENMAT(NN)%MAT(IND1,IND2,IND3)=XSVAR+XDELTA*REAL(IX,KIND=8)
!!$XENERGY=0.D0
!
!     ==========================================================================
!     ==  CALCULATE ENERGY                                                    ==
!     ==========================================================================
      CALL TIMING$CLOCKON('NTBOETOT')
      SWITCH=1
      IF(SWITCH.EQ.1) THEN
!       == TAILED, NONSPHERICAL ORBITALS =======================================
        CALL LMTO_SIMPLEENERGYTEST2()
      ELSE IF(SWITCH.EQ.2) THEN
        CALL LMTO_SIMPLEENERGYTEST()
      ELSE IF(SWITCH.EQ.3) THEN
        CALL LMTO_ENERGYTEST()
!!$      ELSE IF(SWITCH.EQ.4) THEN
!!$!       == UE MULTICENTER EXPANSIONS  ==========================================
!!$        CALL LMTO_ENERGYFULLORB()
      END IF
      CALL TIMING$CLOCKOFF('NTBOETOT')
!!$WRITE(*,*)XDELTA*REAL(IX,8),HAMIL(NN)%MAT(IND1,IND2,IND3)*4.D0 !FACTOR FOUR TO COMPENSATE HFWEIGHT
!!$IF(IX.EQ.3) STOP 'FORCED'
!!$DEALLOCATE(HAMIL)
!!$GOTO 1000
!!$STOP
!
!     ==========================================================================
!     ==  CONVERT HAMIL INTO HTBC
!     ==========================================================================
      CALL TIMING$CLOCKON('NTBODENMATDER')
      CALL LMTO_NTBODENMATDER()
      CALL TIMING$CLOCKOFF('NTBODENMATDER')
!
!!$CALL LMTO$REPORTORTHODENMAT(6)
!!$IF(TOFFSITE)CALL LMTO$REPORTOVERLAP(6)
!!$CALL LMTO$REPORTDENMAT(6)
!!$CALL LMTO$REPORTHAMIL(6)
!STOP 'FORCED'
!
!     ==========================================================================
!     ==  CLEAN DENMAT AND HAMIL                                              ==
!     ==========================================================================
      CALL LMTO_CLEANDENMAT()
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.......
!!$      SUBROUTINE LMTO_FAKEDENMAT()
!!$!     ***********************************************************************
!!$!     **  UNFINISHED TEST VERSION                                          **
!!$!     **  ENFORCE A CERTAIN DENSITY MATRIX                                 **
!!$!     ***********************************************************************
!!$      USE LMTO_MODULE, ONLY : ISPECIES,DENMAT,LNX,LOX
!!$      IMPLICIT NONE
!!$      TYPE LIST_TYPE 
!!$        CHARACTER(32) :: NAME
!!$        INTEGER(4)    :: L
!!$        INTEGER(4)    :: IM
!!$        INTEGER(4)    :: IS
!!$      END TYPE LIST_TYPE
!!$      INTEGER(4),PARAMETER :: NENTRY
!!$      TYPE(LIST_TYPE)      :: LIST(NENTRY)
!!$      INTEGER(4) :: NND
!!$      INTEGER(4) :: NN
!!$      INTEGER(4) :: IAT1,IAT2,IT(3)
!!$      INTEGER(4) :: ISP
!!$      CHARACTER(32) :: NAME  ! ATOM NAME
!!$!     ***********************************************************************
!!$      NND=SIZE(DENMAT)
!!$      DO NN=1,NND
!!$        IAT1=DENMAT(NN)%IAT1
!!$        IAT2=DENMAT(NN)%IAT2
!!$        IT(:)=DENMAT(NN)%IT(:)
!!$        IF(IAT2.NE.IAT1.OR.ABS(SUM(IT**2)).NE.0) CYCLE ! ONSITE ELEMENTS ONLY
!!$        ISP=ISPECIES(IAT1)
!!$        CALL ATOMLIST$GETCH('NAME',IAT1,NAME)
!!$        DO I=1,NENTRY
!!$          IF(NAME.NE.LIST%NAME) CYCLE
!!$        ENDDO
!!$      ENDDO
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SIMPLEENERGYTEST2()
!     **************************************************************************
!     **  WORK OUT THE ENERGY USING THE LOCAL APPROXIMATION                   **
!     **  TAILED PARTIAL WAVES                                                **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,DENMAT,HAMIL,LNX,LOX,POTPAR,TOFFSITE &
     &                       ,HYBRIDSETTING,HFWEIGHT
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      LOGICAL(4),PARAMETER  :: TPLOT=.FALSE.
      INTEGER(4)            :: NND
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IND,INH
      INTEGER(4)            :: LMNX
      INTEGER(4)            :: NDIMD
      INTEGER(4)            :: IT(3)
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: LNX1,LMRX,LRX
      INTEGER(4)            :: LNXT,LMNXT
      INTEGER(4),ALLOCATABLE:: LOXT(:)
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      REAL(8)   ,ALLOCATABLE:: DT(:,:,:)
      REAL(8)   ,ALLOCATABLE:: HT(:,:,:)
      REAL(8)   ,ALLOCATABLE:: DTALL(:,:,:)
      REAL(8)   ,ALLOCATABLE:: HTALL(:,:,:)
      REAL(8)   ,ALLOCATABLE:: AECORE(:)
      REAL(8)   ,ALLOCATABLE:: ULITTLE(:,:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: U(:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: D(:,:,:)
      REAL(8)   ,ALLOCATABLE:: H(:,:,:)
      REAL(8)               :: EH,EX,EXTOT,EHTOT,EHDC,EXDC,Q
      INTEGER(4)            :: NN,IAT,I,J,K,L,IS,IAT1,IAT2,ISP
      INTEGER(4)            :: LMN,LN,IM
      REAL(8)               :: QSPIN(4)
      REAL(8)               :: SVAR
      INTEGER(4)            :: IDFTTYPE
CHARACTER(128) :: STRING
REAL(8)   ,ALLOCATABLE:: T(:,:),UNT(:,:),MYMAT(:,:)
      REAL(8)               :: HFSCALE
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_SIMPLEENERGYTEST2')
PRINT*,'============ ENERGYTEST2 ============================='
      NAT=SIZE(ISPECIES)
      IF(ALLOCATED(HAMIL)) THEN
        CALL ERROR$MSG('HAMIL IS ALLOCATED')
        CALL ERROR$STOP('LMTO_SIMPLEENERGYTEST2')
      END IF
      NND=SIZE(DENMAT)
      ALLOCATE(HAMIL(NND))
      DO NN=1,NND
        HAMIL(NN)%IAT1=DENMAT(NN)%IAT1
        HAMIL(NN)%IAT2=DENMAT(NN)%IAT2
        HAMIL(NN)%IT=DENMAT(NN)%IT
        N1=DENMAT(NN)%N1
        N2=DENMAT(NN)%N2
        N3=DENMAT(NN)%N3
        HAMIL(NN)%N1=N1
        HAMIL(NN)%N2=N2
        HAMIL(NN)%N3=N3
        ALLOCATE(HAMIL(NN)%MAT(N1,N2,N3))
        HAMIL(NN)%MAT=0.D0
      ENDDO

      EXTOT=0.D0
      EHTOT=0.D0
      DO IAT=1,NAT
!
!       == FIND LOCAL DENSITY MATRIX ===========================================
        IND=-1
        DO NN=1,NND
          IF(DENMAT(NN)%IAT1.NE.IAT) CYCLE
          IF(DENMAT(NN)%IAT2.NE.IAT) CYCLE
          IF(SUM(DENMAT(NN)%IT**2).NE.0) CYCLE
          IND=NN
          EXIT
        ENDDO
        IF(IND.LE.0) THEN
          CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED: IND<0')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$I4VAL('NND',NND)
          CALL ERROR$STOP('LMTO_ENERGYTEST')
        END IF
!
!       == FIND LOCAL HAMILTONIAN   ===========================================
        INH=-1
        DO NN=1,NND
          IF(HAMIL(NN)%IAT1.NE.IAT) CYCLE
          IF(HAMIL(NN)%IAT2.NE.IAT) CYCLE
          IF(SUM(HAMIL(NN)%IT**2).NE.0) CYCLE
          INH=NN
          EXIT
        ENDDO
        IF(INH.LE.0) THEN
          CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED: INH<0')
          CALL ERROR$STOP('LMTO_ENERGYTEST')
        END IF
!
        ISP=ISPECIES(IAT)
        LMNX=DENMAT(IND)%N1
        NDIMD=DENMAT(IND)%N3
!
!       == ADJUSTMENT IF LOCAL HFWEIGHT IS DIFFERENT FROM GLOBAL HFWEIGHT ======
        IF(HYBRIDSETTING(ISP)%LHFWEIGHT.GE.0.D0) THEN
          HFSCALE=HYBRIDSETTING(ISP)%LHFWEIGHT/HFWEIGHT
        ELSE
          HFSCALE=1.D0
        END IF
PRINT*,'HFSCALE',IAT,HFSCALE
!
!       ========================================================================
!       == CALCULATE U-TENSOR                                                 ==
!       ========================================================================
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(AECORE(NR))
        CALL SETUP$GETR8A('AECORE',NR,AECORE)
        CALL SETUP$GETI4('LMRX',LMRX)
        CALL SETUP$UNSELECT()
        LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
        ALLOCATE(D(LMNX,LMNX,NDIMD))
        ALLOCATE(H(LMNX,LMNX,NDIMD))
        D=DENMAT(IND)%MAT
        H(:,:,:)=0.D0
!!$DO I=1,LMNX
!!$  WRITE(*,FMT='("D=",30F10.5)')D(I,:,1)
!!$ENDDO
!!$ALLOCATE(T(LMNX,LMNX))
!!$ALLOCATE(UNT(LMNX,LMNX))
!!$ALLOCATE(MYMAT(LMNX,LMNX))
!!$CALL LMTO_ONSORTHO(IAT,LMNX,T,UNT)
!!$DO I=1,NDIMD
!!$MYMAT=MATMUL(UNT,MATMUL(D(:,:,I),TRANSPOSE(UNT)))
!!$PRINT*,'============IAT=',IAT,'  IDIM=',I,'==================='
!!$DO J=1,LMNX
!!$  WRITE(*,FMT='("D(MY)=",30F10.5)')MYMAT(J,:)
!!$ENDDO
!!$ENDDO
!!$DEALLOCATE(T)
!!$DEALLOCATE(UNT)
!!$DEALLOCATE(MYMAT)
!
        LNXT=POTPAR(ISP)%TAILED%LNX
        LMNXT=POTPAR(ISP)%TAILED%LMNX
        ALLOCATE(LOXT(LNXT))
        LOXT=POTPAR(ISP)%TAILED%LOX
        ALLOCATE(DT(LMNXT,LMNXT,NDIMD))
        ALLOCATE(HT(LMNXT,LMNXT,NDIMD))
        CALL LMTO_BLOWUPDENMATNL(IAT,IAT,NDIMD,LMNX,LMNX,D,LMNXT,LMNXT,DT)
!
!       ========================================================================
!       == CALCULATE U-TENSOR                                                 ==
!       ========================================================================
        ALLOCATE(U(LMNXT,LMNXT,LMNXT,LMNXT))
        U=POTPAR(ISP)%TAILED%U
!
!       ========================================================================
!       == WORK OUT TOTAL ENERGY AND HAMILTONIAN                              ==
!       ========================================================================
        HAMIL(INH)%MAT=0.D0
        EH=0.D0
        EX=0.D0
        QSPIN=0.D0
        HT(:,:,:)=0.D0
        DO I=1,LMNXT
          DO J=1,LMNXT
            QSPIN(:NDIMD)=QSPIN(:NDIMD) &
       &                 +POTPAR(ISP)%TAILED%OVERLAP(I,J)*DT(J,I,:)
            DO K=1,LMNXT
              DO L=1,LMNXT
!               ================================================================
!               == HARTREE TERM (NOT CONSIDERED)                              ==
!               == EH=EH+0.5D0*U(I,J,K,L)*D(K,I,1)*D(L,J,1)                   ==
!               == H(K,I,1)=H(K,I,1)+0.5D0*U(I,J,K,L)*D(L,J,1) !DE/DRHO(K,I,1)==
!               == H(L,J,1)=H(L,J,1)+0.5D0*U(I,J,K,L)*D(K,I,1) !DE/DRHO(L,J,1)==
!               ================================================================
!               ================================================================
!               == EXCHANGE ENERGY =============================================
!               ================================================================
!               == AN ADDITIONAL FACTOR COMES FROM THE REPRESENTATION INTO TOTAL AND SPIN
                SVAR=-0.25D0*U(I,J,K,L)
                EX=EX+SVAR*SUM(DT(K,J,:)*DT(I,L,:))
                HT(K,J,:)=HT(K,J,:)+SVAR*DT(I,L,:) 
                HT(I,L,:)=HT(I,L,:)+SVAR*DT(K,J,:) 
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,NDIMD,LMNXT,LMNXT,HT,LMNX,LMNX,H)
        HAMIL(INH)%MAT=HAMIL(INH)%MAT+H*HFSCALE
        EXTOT=EXTOT+EX*HFSCALE
PRINT*,'TOTAL CHARGE ON ATOM=                 ',IAT,QSPIN(1)
PRINT*,'TOTAL SPIN[HBAR/2] ON ATOM=           ',IAT,QSPIN(2:NDIMD)
PRINT*,'EXACT EXCHANGE ENERGY FOR ATOM=       ',IAT,EX
!
!       ========================================================================
!       == DOUBLE COUNTING CORRECTION (EXCHANGE ONLY)                         ==
!       ========================================================================
CALL DFT$GETI4('TYPE',IDFTTYPE)
PRINT*,'IDFTTYPE ',IDFTTYPE
IF(IDFTTYPE.NE.5002) THEN
! THIS IS THE TIME CONSUMIN PART OF ENERGYTEST
CALL TIMING$CLOCKON('ENERGYTEST:DC')      
        ALLOCATE(DTALL(LMNXT,LMNXT,NDIMD))
        ALLOCATE(HTALL(LMNXT,LMNXT,NDIMD))
        POTPAR(ISP)%TALLORB=.TRUE.
        CALL LMTO_BLOWUPDENMATNL(IAT,IAT,NDIMD,LMNX,LMNX,D,LMNXT,LMNXT,DTALL)
        CALL LMTO_SIMPLEDC(GID,NR,LMNXT,LNXT,LOXT,POTPAR(ISP)%TAILED%AEF &
     &                    ,LRX,AECORE,DT,DTALL,EX,HT,HTALL)
        CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,NDIMD,LMNXT,LMNXT,HTALL,LMNX,LMNX,H)
        POTPAR(ISP)%TALLORB=.FALSE.  ! DO NOT FORGET THIS!!!!!
        HAMIL(INH)%MAT=HAMIL(INH)%MAT-H*HFSCALE
        CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,NDIMD,LMNXT,LMNXT,HT,LMNX,LMNX,H)
        HAMIL(INH)%MAT=HAMIL(INH)%MAT-H*HFSCALE
        EXTOT=EXTOT-EX*HFSCALE
        DEALLOCATE(DTALL)
        DEALLOCATE(HTALL)
PRINT*,'DOUBLE COUNTING CORRECTION ENERGY FOR ATOM=',IAT,-EX
CALL TIMING$CLOCKOFF('ENERGYTEST:DC')      
END IF
!
!       ========================================================================
!       == ADD CORE VALENCE EXCHANGE                                          ==
!       ========================================================================
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(.NOT.POTPAR(ISP)%TORB(LN)) THEN
            D(:,LMN+1:LMN+2*L+1,:)=0.D0
            D(LMN+1:LMN+2*L+1,:,:)=0.D0
          END IF
          LMN=LMN+2*L+1
        ENDDO
!
        CALL LMTO_CVX(ISP,LMNX,EX,D(:,:,1),H(:,:,1))

        IF(.NOT.HYBRIDSETTING(ISP)%TCV) THEN
          EX=0.D0
          H(:,:,:)=0.D0
        END IF
!
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(.NOT.POTPAR(ISP)%TORB(LN)) THEN
            H(:,LMN+1:LMN+2*L+1,:)=0.D0
            H(LMN+1:LMN+2*L+1,:,:)=0.D0
          END IF
          LMN=LMN+2*L+1
        ENDDO
        EXTOT=EXTOT+EX*HFSCALE
        HAMIL(INH)%MAT(:,:,1)=HAMIL(INH)%MAT(:,:,1)+H(:,:,1)*HFSCALE
PRINT*,'CORE VALENCE EXCHANGE ENERGY FOR ATOM=',IAT,EX
!
        DEALLOCATE(HT)
        DEALLOCATE(DT)
        DEALLOCATE(U)
        DEALLOCATE(H)
        DEALLOCATE(D)
        DEALLOCATE(AECORE)
        DEALLOCATE(LOXT)
      ENDDO
!
!     ==========================================================================
!     == OFFSITE EXCHANGE CONTRIBUTION                                        ==
!     ==========================================================================
      IF(TOFFSITE) THEN
        PRINT*,'NOW INTO OFFSITEX'
        IF(1.EQ.0) THEN
          CALL LMTO_OFFSITEX(EX)      ! USING GAUSS INTERPOLATION
        ELSE
          CALL LMTO_OFFSITEXEVAL(EX)  ! USING NUMERICAL MATRIX ELEMENTS
        END IF
        PRINT*,'+-+-+-+  OFFSITE EX=',EX
        EXTOT=EXTOT+EX
      END IF
!
!     ==========================================================================
!     == MAKE HAMILTONIAN HERMITESCH                                          ==
!     ==========================================================================
!      CALL LMTO$SYMMETRIZEHAMIL()
!
!     ==========================================================================
!     == RESCALE WITH HFWEIGHT                                                ==
!     ==========================================================================
      DO NN=1,NND
        HAMIL(NN)%MAT=HAMIL(NN)%MAT*HFWEIGHT
      ENDDO
      EXTOT=EXTOT*HFWEIGHT
      EHTOT=EHTOT*HFWEIGHT
!
!     ==========================================================================
!     == COMMUNICATE ENERGY TO ENERGYLIST                                     ==
!     ==========================================================================
      CALL ENERGYLIST$SET('LMTO INTERFACE',EXTOT)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EXTOT)
!
!     ==========================================================================
!     == PLOT WAVE FUNCTIONS                                                  ==
!     ==========================================================================
      IF(TPLOT) THEN
        PRINT*,'BEFORE PLOTTAILED'
        CALL LMTO_PLOTTAILED()
        PRINT*,'LMTO_GRIDPLOT_TAILED'
        CALL LMTO_GRIDPLOT_TAILED(1)
        CALL ERROR$MSG('PLANNED EXIT AFTER PLOTTING FOR ANALYSIS')
        CALL ERROR$STOP('LMTO_SIMPLEENERGYTEST2')
      END IF
!
!     ==========================================================================
!     == PRINT FOR TEST                                                       ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,FMT='(82("="),T10," DENSITY MATRIX IN A NTBO BASIS ")')
        WRITE(*,FMT='(82("="),T10," ONLY ONSITE TERMS WRITTEN ")')
        DO NN=1,NND
          IAT1=DENMAT(NN)%IAT1
          IAT2=DENMAT(NN)%IAT2
          IT=DENMAT(NN)%IT
          IF(IAT1.NE.IAT2.OR.IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
     &                                                           IAT1,IAT2,IT
          N1=DENMAT(NN)%N1
          N2=DENMAT(NN)%N2
          DO I=1,1 !DENMAT(NN)%N3
            DO J=1,N2 
              WRITE(*,FMT='(I3,30F10.3)')I,DENMAT(NN)%MAT(J,:,I)
            ENDDO
            WRITE(*,FMT='(82("-"))')
          ENDDO
        ENDDO
        WRITE(*,FMT='("XC ENERGY ",F10.5)')EXTOT
        WRITE(*,FMT='(82("="),T10," HAMILTON MATRIX IN A NTBO BASIS ")')
        WRITE(*,FMT='(82("="),T10," ONLY ONSITE TERMS WRITTEN ")')
        DO NN=1,NND
          IAT1=HAMIL(NN)%IAT1
          IAT2=HAMIL(NN)%IAT2
          IT=HAMIL(NN)%IT
          IF(IAT1.NE.IAT2.OR.IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
     &                                                            IAT1,IAT2,IT
          N1=HAMIL(NN)%N1
          N2=HAMIL(NN)%N2
          DO I=1,1 !HAMIL(NN)%N3
            DO J=1,N2 
              WRITE(*,FMT='(I3,30F10.3)')I,HAMIL(NN)%MAT(J,:,I)
            ENDDO
            WRITE(*,FMT='(82("-"))')
          ENDDO
        ENDDO
      END IF
!
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_BLOWUPDENMATNL(IAT1,IAT2,NDIMD,LMNX1,LMNX2,D &
     &                                              ,LMNXT1,LMNXT2,DT)
!     **************************************************************************
!     ** BRINGS THE DENSITY MATRIX EXPRESSED                                  **
!     **   IN LOCAL ORBITALS |CHI_I> (WITH MIXED ANGULAR MOMENTUM CHARACTER)  **
!     ** INTO A LARGER SET |Y_J> (WITH PURE ANGULAR-MOMENTUM CHARACTER)       **
!     **                                                                      **
!     **    |CHI_I>= |Y_I> - |Y_N+J>*TRANSPOSE(SBAR)_J,I                      **
!     **       (THE TRANSPOSE IS NOT YET IMPLEMENTED, FUTURE NOTATION)        **
!     **                                                                      **
!     ** DENSITY MATRIX ELEMENTS WITH TORB=F ARE SET TO ZERO                  **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY: ISPECIES,LNX,LOX,POTPAR,SBAR,SBARLI1,TSPHERICAL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT1
      INTEGER(4),INTENT(IN)  :: IAT2
      INTEGER(4),INTENT(IN)  :: NDIMD
      INTEGER(4),INTENT(IN)  :: LMNX1
      INTEGER(4),INTENT(IN)  :: LMNX2
      INTEGER(4),INTENT(IN)  :: LMNXT1
      INTEGER(4),INTENT(IN)  :: LMNXT2
      REAL(8)   ,INTENT(IN)  :: D(LMNX1,LMNX2,NDIMD)
      REAL(8)   ,INTENT(OUT) :: DT(LMNXT1,LMNXT2,NDIMD)
      REAL(8)                :: D1(LMNX1,LMNX2,NDIMD)
      INTEGER(4)             :: NNS
      REAL(8)   ,ALLOCATABLE :: SBARLOC1(:,:)
      REAL(8)   ,ALLOCATABLE :: SBARLOC2(:,:)
      INTEGER(4)             :: I,IAT,ISP,NN,LMN1,LMN2,LMNDOT1,LMNDOT2
      INTEGER(4)             :: L,LMN,LN,I1,IM,N1,N2
      REAL(8)                :: SVAR
      LOGICAL(4)             :: TLEFT
!     **************************************************************************
      D1(:,:,:)=D(:,:,:)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      DO I=1,2
        TLEFT=(I.EQ.1)
        IAT=IAT2
        IF(TLEFT) IAT=IAT1
!
!       ========================================================================
!       ==  DETERMINE SIZE OF STRUCTURE CONSTANT ARRAY                        ==
!       ========================================================================
        ISP=ISPECIES(IAT)
        N1=0
        N2=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          N2=N2+2*L+1
          IF(POTPAR(ISP)%LNSCATT(LN).EQ.LN) N1=N1+2*L+1
        ENDDO
        IF(TLEFT) THEN
          ALLOCATE(SBARLOC1(N2,N1))    !C
        ELSE
          ALLOCATE(SBARLOC2(N2,N1))    !C
        END IF
!    
!       ========================================================================
!       ==  COLLECT LOCAL STRUCTURE CONSTANTS                                 ==
!       ========================================================================
        NNS=SIZE(SBAR)
        DO NN=1,NNS
          IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
          IF(SBAR(NN)%IAT2.NE.IAT) CYCLE
          IF(SUM(SBAR(NN)%IT(:)**2).NE.0) CYCLE
          IF(SBAR(NN)%N1.NE.N1.OR.SBAR(NN)%N2.NE.N1) THEN
            CALL ERROR$MSG('INCONSISTENT ARRAY SIZES N1,N2')
            CALL ERROR$I4VAL('N1',N1)
            CALL ERROR$I4VAL('SBAR%N1',SBAR(NN)%N1)
            CALL ERROR$I4VAL('SBAR%N2',SBAR(NN)%N2)
            CALL ERROR$STOP('LMTO_BLOWUPDENMATNL')
          END IF
          LMN=0
          DO LN=1,LNX(ISP)
            L=LOX(LN,ISP)
            I1=SBARLI1(L+1,ISP)-1
            DO IM=1,2*L+1 
              IF(TLEFT) THEN
                SBARLOC1(LMN+IM,:)=SBAR(NN)%MAT(I1+IM,:)   !C
              ELSE
                SBARLOC2(LMN+IM,:)=SBAR(NN)%MAT(I1+IM,:)   !C
              END IF
            ENDDO
            LMN=LMN+2*L+1
          ENDDO
          EXIT
        ENDDO
!
!       ========================================================================
!       ==  REMOVE ORBITALS NOT IN THE SET                                    ==
!       ========================================================================
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(.NOT.(POTPAR(ISP)%TORB(LN).OR.POTPAR(ISP)%TALLORB)) THEN
            IF(TLEFT) THEN
              D1(LMN+1:LMN+2*L+1,:,:)=0.D0
            ELSE
              D1(:,LMN+1:LMN+2*L+1,:)=0.D0
            END IF
          END IF
          LMN=LMN+2*L+1
        ENDDO
!
      ENDDO
!
!     ==========================================================================
!     ==  BLOW UP DENSITY MATRIX                                              ==
!     ==========================================================================
      DO I=1,NDIMD
        DT(:,:,I)=0.D0
        DT(:LMNX1,:LMNX2,I)=D1(:,:,I)
        DT(LMNX1+1:,:LMNX2,I)=-MATMUL(TRANSPOSE(SBARLOC1),D1(:,:,I)) !C
        DT(:LMNX1,LMNX2+1:,I)=-MATMUL(D1(:,:,I),SBARLOC2)            !C
        DT(LMNX1+1:,LMNX2+1:,I)=-MATMUL(TRANSPOSE(SBARLOC1),DT(:LMNX1,LMNX2+1:,I)) !C
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_BLOWUPPSINL(IAT,LMNX,C,LMNXT,CT)
!     **************************************************************************
!     ** BRINGS THE DENSITY MATRIX EXPRESSED                                  **
!     **   IN LOCAL ORBITALS |CHI_I> (WITH MIXED ANGULAR MOMENTUM CHARACTER)  **
!     ** INTO A LARGER SET |Y_J> (WITH PURE ANGULAR-MOMENTUM CHARACTER)       **
!     **                                                                      **
!     **    |CHI_I>= |Y_I> - |Y_N+J>*TRANSPOSE(SBAR)_J,I                      **
!     **       (THE TRANSPOSE IS NOT YET IMPLEMENTED, FUTURE NOTATION)        **
!     **                                                                      **
!     ** DENSITY MATRIX ELEMENTS WITH TORB=F ARE SET TO ZERO                  **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY: ISPECIES,LNX,LOX,POTPAR,SBAR,SBARLI1,TSPHERICAL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      INTEGER(4),INTENT(IN)  :: LMNX
      INTEGER(4),INTENT(IN)  :: LMNXT
      REAL(8)   ,INTENT(IN)  :: C(LMNX)
      REAL(8)   ,INTENT(OUT) :: CT(LMNXT)
      REAL(8)                :: C1(LMNX)
      INTEGER(4)             :: NNS
      REAL(8)   ,ALLOCATABLE :: SBARLOC(:,:)
      INTEGER(4)             :: I,ISP,NN,LMN1
      INTEGER(4)             :: L,LMN,LN,I1,IM,N1,N2
      REAL(8)                :: SVAR
      LOGICAL(4)             :: TLEFT
!     **************************************************************************
      C1(:)=C(:)
!
!     ========================================================================
!     ==  DETERMINE SIZE OF STRUCTURE CONSTANT ARRAY                        ==
!     ========================================================================
      ISP=ISPECIES(IAT)
      N1=0
      N2=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        N2=N2+2*L+1
        IF(POTPAR(ISP)%LNSCATT(LN).EQ.LN) N1=N1+2*L+1
      ENDDO
      ALLOCATE(SBARLOC(N2,N1))    !C
!    
!     ========================================================================
!     ==  COLLECT LOCAL STRUCTURE CONSTANTS                                 ==
!     ========================================================================
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
        IF(SBAR(NN)%IAT2.NE.IAT) CYCLE
        IF(SUM(SBAR(NN)%IT(:)**2).NE.0) CYCLE
        IF(SBAR(NN)%N1.NE.N1.OR.SBAR(NN)%N2.NE.N1) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZES N1,N2')
          CALL ERROR$I4VAL('N1',N1)
          CALL ERROR$I4VAL('SBAR%N1',SBAR(NN)%N1)
          CALL ERROR$I4VAL('SBAR%N2',SBAR(NN)%N2)
          CALL ERROR$STOP('LMTO_BLOWUPPSINL')
        END IF
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I1=SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1 
            SBARLOC(LMN+IM,:)=SBAR(NN)%MAT(I1+IM,:)   !C
          ENDDO
          LMN=LMN+2*L+1
        ENDDO
        EXIT
      ENDDO
!
!     ========================================================================
!     ==  REMOVE ORBITALS NOT IN THE SET                                    ==
!     ========================================================================
      LMN=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        IF(.NOT.POTPAR(ISP)%TORB(LN)) THEN
          C1(LMN+1:LMN+2*L+1)=0.D0
        END IF
        LMN=LMN+2*L+1
      ENDDO
!
!     ==========================================================================
!     ==  BLOW UP DENSITY MATRIX                                              ==
!     ==========================================================================
      CT(:LMNX)=C1(:)
      CT(LMNX+1:)=-MATMUL(TRANSPOSE(SBARLOC),C1)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SHRINKDOWNHTNL(IAT1,IAT2,NDIMD,LMNXT1,LMNXT2,HT &
     &                              ,LMNX1,LMNX2,H)
!     **************************************************************************
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY: ISPECIES,LNX,LOX,POTPAR,SBAR,SBARLI1,TSPHERICAL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT1    
      INTEGER(4),INTENT(IN)  :: IAT2    
      INTEGER(4),INTENT(IN)  :: NDIMD
      INTEGER(4),INTENT(IN)  :: LMNX1
      INTEGER(4),INTENT(IN)  :: LMNX2
      INTEGER(4),INTENT(IN)  :: LMNXT1
      INTEGER(4),INTENT(IN)  :: LMNXT2
      REAL(8)   ,INTENT(IN)  :: HT(LMNXT1,LMNXT2,NDIMD)
      REAL(8)   ,INTENT(OUT) :: H(LMNX1,LMNX2,NDIMD)
      INTEGER(4)             :: NNS
      REAL(8)   ,ALLOCATABLE :: SBARLOC1(:,:)
      REAL(8)   ,ALLOCATABLE :: SBARLOC2(:,:)
      REAL(8)                :: SVAR
      INTEGER(4)             :: I,IAT,ISP,NN,LMN1,LMN2,LMN3
      INTEGER(4)             :: L,LMN,LN,I1,IM,N1,N2
      LOGICAL(4)             :: TLEFT
!     **************************************************************************
      DO I=1,2
        TLEFT=(I.EQ.1)
        IAT=IAT2
        IF(TLEFT)IAT=IAT1
!
!       ========================================================================
!       ==  DETERMINE SIZE OF STRUCTURE CONSTANT ARRAY                        ==
!       ========================================================================
        ISP=ISPECIES(IAT)
        N1=0
        N2=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          N2=N2+2*L+1
          IF(POTPAR(ISP)%LNSCATT(LN).EQ.LN) N1=N1+2*L+1
        ENDDO
        IF(TLEFT) THEN
          ALLOCATE(SBARLOC1(N2,N1))
        ELSE
          ALLOCATE(SBARLOC2(N2,N1))
        END IF
!  
!       ========================================================================
!       ==  COLLECT LOCAL STRUCTURE CONSTANTS                                 ==
!       ========================================================================
        NNS=SIZE(SBAR)
        DO NN=1,NNS
          IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
          IF(SBAR(NN)%IAT2.NE.IAT) CYCLE
          IF(SUM(SBAR(NN)%IT(:)**2).NE.0) CYCLE
          IF(SBAR(NN)%N1.NE.N1.OR.SBAR(NN)%N2.NE.N1) THEN
            CALL ERROR$MSG('INCONSISTENT ARRAY SIZES N1,N2')
            CALL ERROR$I4VAL('N1',N1)
            CALL ERROR$I4VAL('SBAR%N1',SBAR(NN)%N1)
            CALL ERROR$I4VAL('SBAR%N2',SBAR(NN)%N2)
            CALL ERROR$STOP('LMTO_SHRINKDOWNHTNL')
          END IF
          LMN=0
          DO LN=1,LNX(ISP)
            L=LOX(LN,ISP)
            I1=SBARLI1(L+1,ISP)-1
            DO IM=1,2*L+1 
              IF(TLEFT) THEN
                SBARLOC1(LMN+IM,:)=SBAR(NN)%MAT(I1+IM,:)
              ELSE
                SBARLOC2(LMN+IM,:)=SBAR(NN)%MAT(I1+IM,:)
              END IF
            ENDDO
            LMN=LMN+2*L+1
          ENDDO
          EXIT
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  SHIRNK DOWN                                                         ==
!     ==========================================================================
      DO I=1,NDIMD
        H(:,:,I)=HT(:LMNX1,:LMNX2,I) &
     &          -MATMUL(HT(:LMNX1,LMNX2+1:,I),TRANSPOSE(SBARLOC2)) &          !C
     &          -MATMUL(SBARLOC1,HT(LMNX1+1:,:LMNX2,I)) &                     !C
     &          +MATMUL(SBARLOC1 &                                            !C
     &                 ,MATMUL(HT(LMNX1+1:,LMNX2+1:,I),TRANSPOSE(SBARLOC2)))  !C
      ENDDO
!
!     ==========================================================================
!     ==  REMOVE ORBITALS NOT IN THE SET                                      ==
!     ==========================================================================
      DO I=1,2
        TLEFT=(I.EQ.1)
        IAT=IAT2
        IF(TLEFT)IAT=IAT1
        ISP=ISPECIES(IAT)
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(.NOT.(POTPAR(ISP)%TORB(LN).OR.POTPAR(ISP)%TALLORB)) THEN
            IF(TLEFT) THEN
              H(LMN+1:LMN+2*L+1,:,:)=0.D0
            ELSE
              H(:,LMN+1:LMN+2*L+1,:)=0.D0
            END IF
          END IF
          LMN=LMN+2*L+1
        ENDDO
      ENDDO

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_CVX(ISP,NORB,EX,D,H)
!     **************************************************************************
!     **  CORE VALENCE EXCHANGE ENERGY                                        **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY: POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: ISP ! ATOM-TYPE INDEX
      INTEGER(4),INTENT(IN)  :: NORB ! #(LOCAL ORBITALS)
      REAL(8)   ,INTENT(IN)  :: D(NORB,NORB)  ! TOTAL DENSITY MATRIX
      REAL(8)   ,INTENT(OUT) :: EX            ! CORE-VALENCE ENERGY
      REAL(8)   ,INTENT(OUT) :: H(NORB,NORB)  ! HAMILTONIAN CONTRIBUTION
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: CVXMAT(:,:)
      REAL(8)                :: C1,C2,SVAR
      INTEGER(4)             :: LN1,LN2,L1,L2,LMN1,LMN2,IM
!     **************************************************************************
!
!     ==========================================================================
!     == COLLECT DATA                                                         ==
!     ==========================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      ALLOCATE(CVXMAT(LNX,LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      CALL SETUP$GETR8A('CVX',LNX*LNX,CVXMAT)
      CALL SETUP$UNSELECT()
!
!     ==========================================================================
!     == CALCULATE CORE-VALENCE EXCHANGE ENERGY                               ==
!     ==========================================================================
      EX=0.D0
      H(:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        C1=POTPAR(ISP)%KTOPHI(LN1)
        LMN2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
          IF(L2.EQ.L1) THEN
            C2=POTPAR(ISP)%KTOPHI(LN2)
            SVAR=C1*CVXMAT(LN1,LN2)*C2
            DO IM=1,2*L1+1
              EX=EX+D(LMN2+IM,LMN1+IM)*SVAR
              H(LMN2+IM,LMN1+IM)=H(LMN2+IM,LMN1+IM)+SVAR
            ENDDO
          END IF
          LMN2=LMN2+2*L2+1
        ENDDO
        LMN1=LMN1+2*L1+1
      ENDDO
!
!     ==========================================================================
!     == CLOSE DOWN                                                           ==
!     ==========================================================================
      DEALLOCATE(CVXMAT)
      DEALLOCATE(LOX)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SIMPLEDC(GID,NR,LMNX,LNX,LOX,CHI,LRX,AECORE &
     &                        ,DENMAT,DENMATB,ETOT,HAM,HAMB)
!     **************************************************************************
!     **  DOUBLE COUNTING CORRECTION FOR THE HYBRID FUNCTIONAL                **
!     **                                                                      **
!     **  DETERMINES THE HARTREE AND EXCHANGE-ONLY ENERGY FROM THE            **
!     **  DFT FUNCTIONAL                                                      **
!     **  FOR THE DENSITY BUILT FROM THE LOCAL ORBITALS AND THE CORE DENSITY  **
!     **  THIS ENERGY NEEDS TO BE SUBTRACTED FROM THE TOTAL ENERGY            **
!     **                                                                      **
!     **  DENMAT DESCRIBES THE CORRELATED ORBITALS, WHILE                     **
!     **  DENMATB DESCRIBES ALL ORBITALS ON THIS SITE                         **
!     **                                                                      **
!     **   DEX=EX(RHOTOT)-EX(RHOTOT-RHOCOR)                                   **
!     **   VXTOT=MU(RHOTOT)-MU(RHOTOT-RHOCOR)                                 **
!     **   VXCOR=MU(RHOTOT-RHOCOR)                                            **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: LRX
      INTEGER(4)  ,INTENT(IN) :: LMNX       ! #(LOCAL ORBITALS)
      INTEGER(4)  ,INTENT(IN) :: LNX        ! #(RADIAL FUNCTIONS)
      INTEGER(4)  ,INTENT(IN) :: LOX(LNX)   !MAIN ANGULAR MOMENTUM OF LOCAL ORB.
      REAL(8)     ,INTENT(IN) :: CHI(NR,LNX)
      REAL(8)     ,INTENT(IN) :: AECORE(NR)
      REAL(8)     ,INTENT(IN) :: DENMAT(LMNX,LMNX,4) ! DENSITY MATRIX
      REAL(8)     ,INTENT(IN) :: DENMATB(LMNX,LMNX,4) ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT       ! DOUBLE COUNTINNG ENERGY
      REAL(8)     ,INTENT(OUT):: HAM(LMNX,LMNX,4)  ! DETOT/D(RHO^*)        
      REAL(8)     ,INTENT(OUT):: HAMB(LMNX,LMNX,4)  ! DETOT/D(RHO^*)        
      INTEGER(4)  ,PARAMETER  :: METHOD=3
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)  ,PARAMETER  :: NDIMD=4
      COMPLEX(8)              :: DENMAT1(LMNX,LMNX,NDIMD)
      COMPLEX(8)              :: HAM1(LMNX,LMNX,NDIMD)
      REAL(8)                 :: R(NR)
      REAL(8)     ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)     ,ALLOCATABLE:: RHO2(:,:,:)
      REAL(8)     ,ALLOCATABLE:: POT2(:,:,:)
      REAL(8)     ,ALLOCATABLE:: RHOWC(:,:,:)
      REAL(8)     ,ALLOCATABLE:: POT(:,:,:)
      REAL(8)                 :: EDENSITY(NR)
      REAL(8)                 :: AUX(NR),SVAR
      REAL(8)                 :: FXC(NR)
      REAL(8)                 :: PI,FOURPI,Y0
      INTEGER(4)              :: LMRX,L
      INTEGER(4)              :: IDIM,LM,LMN
      REAL(8)                 :: ETOTC,ETOTV
INTEGER(4) :: LMRX1
INTEGER(4) :: IMETHOD
 REAL(8)     ,ALLOCATABLE:: RHOTEST(:,:,:)
 REAL(8)     ,ALLOCATABLE:: POTTEST(:,:,:)
 REAL(8)     ,ALLOCATABLE:: RHOTEST2(:,:,:)
 REAL(8)     ,ALLOCATABLE:: POTTEST2(:,:,:)
 REAL(8)                 :: ETOT2
!     **************************************************************************
      LMRX=(LRX+1)**2
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(4.D0*PI)
      ETOT=0.D0
!
!     ==========================================================================
!     ==  TRANSFORM DENSITY MATRIX FROM UP/DOWN TO TOTAL/SPIN                 ==
!     ==========================================================================
      DENMAT1=CMPLX(DENMAT)
!
!     ==========================================================================
!     ==  CALCULATE DENSITY                                                   ==
!     ==========================================================================
      ALLOCATE(RHO(NR,LMRX,NDIMD))
      DO IDIM=1,NDIMD
        CALL AUGMENTATION_RHO(NR,LNX,LOX,CHI &
     &                       ,LMNX,DENMAT1(:,:,IDIM),LMRX,RHO(:,:,IDIM))
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE ENERGY AND POTENTIAL                                      ==
!     ==========================================================================
      IF(METHOD.EQ.1) THEN 
!       == THIS IS THE SAME METHOD AS 1000 JUST WITHOUT COMMENTS ===============
        ALLOCATE(RHOWC(NR,LMRX,NDIMD))  !WITH CORE
        RHOWC=RHO
        RHOWC(:,1,1)=RHO(:,1,1)+AECORE(:)
        ALLOCATE(POT(NR,LMRX,NDIMD))
        CALL DFT$SETL4('XCONLY',.TRUE.)
        CALL AUGMENTATION_XC(GID,NR,1,1,AECORE,ETOTC,POT)
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,ETOTV,POT)
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOWC,ETOT,POT)
        CALL DFT$SETL4('XCONLY',.FALSE.)
        ETOT=ETOT-ETOTC
        CALL LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
        DEALLOCATE(POT)
        HAM=REAL(HAM1)

      ELSE IF(METHOD.EQ.2) THEN 
        DENMAT1=CMPLX(DENMATB)
        ALLOCATE(RHOWC(NR,LMRX,NDIMD))  !WITH CORE
        DO IDIM=1,NDIMD
          CALL AUGMENTATION_RHO(NR,LNX,LOX,CHI &
     &                       ,LMNX,DENMAT1(:,:,IDIM),LMRX,RHOWC(:,:,IDIM))
        ENDDO
        RHOWC(:,1,1)=RHOWC(:,1,1)+AECORE(:)
!
CALL RADIAL$R(GID,NR,R)
AUX(:)=4.D0*PI*R**2*AECORE*Y0
CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
PRINT*,'INTEGRAL OF CORE',SVAR
AUX(:)=4.D0*PI*R**2*RHOWC(:,1,1)*Y0
CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
PRINT*,'INTEGRAL OF RHOTOT',SVAR
AUX(:)=4.D0*PI*R**2*RHO(:,1,1)*Y0
CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
PRINT*,'INTEGRAL OF RHOCOR',SVAR
!!$STOP
        ALLOCATE(POT(NR,LMRX,NDIMD))
        ALLOCATE(POT2(NR,LMRX,NDIMD))
        CALL DFT$SETL4('XCONLY',.TRUE.)
CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,ETOTV,POT2)
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOWC-RHO,ETOTC,POT)
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOWC,ETOT,POT2)
        CALL DFT$SETL4('XCONLY',.FALSE.)
PRINT*,'EX(TOT)=',ETOT,' EX(TOT-COR)=',ETOTC,' DEX=',ETOT-ETOTC,' EX(COR)=',ETOTV
        ETOT=ETOT-ETOTC
        CALL LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
        HAM=REAL(HAM1)
        POT=POT2-POT
        CALL LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
        HAMB=REAL(HAM1)
        DEALLOCATE(POT)
        DEALLOCATE(POT2)
!
      ELSE IF(METHOD.EQ.3) THEN
!       ==  CONSTRUCT FULL DENSITY   =====================================
        DENMAT1=CMPLX(DENMATB)
        ALLOCATE(RHOWC(NR,LMRX,NDIMD))  !WITH CORE
        DO IDIM=1,NDIMD
          CALL AUGMENTATION_RHO(NR,LNX,LOX,CHI &
     &                       ,LMNX,DENMAT1(:,:,IDIM),LMRX,RHOWC(:,:,IDIM))
        ENDDO
        RHOWC(:,1,1)=RHOWC(:,1,1)+AECORE(:)
!       ==  POTENTIAL FOR THE FULL DENSITY =====================================
        ALLOCATE(POT(NR,LMRX,NDIMD))
        CALL LMTO_RADXC(GID,NR,LMRX,NDIMD,RHOWC,FXC,POT)
        AUX(:)=(RHO(:,1,1)/(RHOWC(:,1,1)+1.D-6))**2
        DO IDIM=1,NDIMD
          DO LM=1,LMRX
            POT(:,LM,IDIM)=AUX(:)*POT(:,LM,IDIM)
          ENDDO
        ENDDO
        POT(:,1,1)=POT(:,1,1)-2.D0*FXC(:)*AUX(:)/(RHOWC(:,1,1)*Y0+1.D-6)/Y0
!       == POTENTIAL FOR THE CORRELATED DENSITY =================================
        ALLOCATE(POT2(NR,LMRX,NDIMD))
        POT2(:,:,:)=0.D0
        POT2(:,1,1)=2.D0*FXC*AUX(:)/(RHO(:,1,1)*Y0+1.D-6)/Y0
        CALL RADIAL$R(GID,NR,R)
        AUX(:)=FOURPI*R(:)**2*FXC*AUX(:)
        CALL RADIAL$INTEGRAL(GID,NR,AUX,ETOT)
PRINT*,'----EXC  ',ETOT
        CALL LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT2,CHI,HAM1)
        HAM=REAL(HAM1)
        CALL LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
        HAMB=REAL(HAM1)
        DEALLOCATE(POT)
        DEALLOCATE(POT2)
!
      ELSE IF(METHOD.EQ.1000) THEN
        ALLOCATE(RHOWC(NR,LMRX,NDIMD))  !WITH CORE
        RHOWC=RHO
        RHOWC(:,1,1)=RHO(:,1,1)+AECORE(:)
!
!     ==========================================================================
!     ==  CALCULATE ENERGY AND POTENTIAL                                      ==
!     ==========================================================================
!     == EXCHANGE ENERGY AND POTENTIAL =========================================
        CALL DFT$SETL4('XCONLY',.TRUE.)
!
!     ==========================================================================
!     == THIS FORMULATION IS BASED ON A NONCOLLINEAR FORMULATION, WHICH       ==
!     == YIELDS DIFFERENT RESULTS FROM A COLLINEAR FORMULATION EVEN FOR       ==
!     == A COLLINEAR DENSITY                                                  ==
!     ==                                                                      ==
!     == THE REASON FOR THIS DIFFERENCE IS THE TRANSFORMATION OF A            ==
!     == NON-COLLINEAR DENSITY WITHIN AUGMENTATION_NCOLLTRANS WHICH IS CALLED ==
!     == BY AUGMENTATION_XC                                                   ==
!     ==                                                                      ==
!     ==========================================================================
        ALLOCATE(POT(NR,LMRX,NDIMD))
        CALL AUGMENTATION_XC(GID,NR,1,1,AECORE,ETOTC,POT)
!!$ALLOCATE(RHO2(NR,LMRX,2))
!!$ALLOCATE(POT2(NR,LMRX,2))
!!$RHO2(:,:,1)=RHO(:,:,1)
!!$RHO2(:,:,2)=RHO(:,:,4)
!!$CALL AUGMENTATION_XC(GID,NR,LMRX,2,RHO2,ETOTV,POT2)
!!$PRINT*,'ETOTV COLLINEAR',ETOTV
!!$CALL AUGMENTATION_WRITEPHI('RHO4_Z.DAT',GID,NR,LMRX,RHO4(:,:,4))
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,ETOTV,POT)
!!$PRINT*,'ETOTV NONCOLLINEAR',ETOTV
!!$PRINT*,'GID,NR,LMRX,NDIMD ',GID,NR,LMRX,NDIMD
!!$RHO2(:,:,1)=RHOWC(:,:,1)
!!$RHO2(:,:,2)=RHOWC(:,:,4)
!!$CALL ATOMLIB_WRITEPHI('RHO2WC_T.DAT',GID,NR,LMRX,RHO2(:,:,1))
!!$CALL ATOMLIB_WRITEPHI('RHO2WC_S.DAT',GID,NR,LMRX,RHO2(:,:,2))
!!$CALL AUGMENTATION_XC(GID,NR,LMRX,2,RHO2,ETOT,POT2)
!!$CALL ATOMLIB_WRITEPHI('POT2WC_T.DAT',GID,NR,LMRX,POT2(:,:,1))
!!$CALL ATOMLIB_WRITEPHI('POT2WC_S.DAT',GID,NR,LMRX,POT2(:,:,2))
!!$PRINT*,'ETOT COLLINEAR',ETOTV
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOWC,ETOT,POT)
!!$CALL ATOMLIB_WRITEPHI('POT4WC_T.DAT',GID,NR,LMRX,POT(:,:,1))
!!$CALL ATOMLIB_WRITEPHI('POT4WC_X.DAT',GID,NR,LMRX,POT(:,:,2))
!!$CALL ATOMLIB_WRITEPHI('POT4WC_Y.DAT',GID,NR,LMRX,POT(:,:,3))
!!$CALL ATOMLIB_WRITEPHI('POT4WC_Z.DAT',GID,NR,LMRX,POT(:,:,4))
!!$PRINT*,'ETOTV NONCOLLINEAR',ETOTV
!POT(:,:,:)=0.D0
!POT(:,:,1)=POT2(:,:,1)
!POT(:,:,4)=POT2(:,:,2)
!!$DEALLOCATE(RHO2)
!!$DEALLOCATE(POT2)
!!$PRINT*,'TOTAL        EXCHANGE ENERGY (LOCAL) ',ETOT
!!$PRINT*,'VALENCE      EXCHANGE ENERGY (LOCAL) ',ETOTV
!!$PRINT*,'CORE         EXCHANGE ENERGY (LOCAL) ',ETOTC
!!$PRINT*,'CORE-VALENCE EXCHANGE ENERGY (LOCAL) ',ETOT-ETOTV-ETOTC
        ETOT=ETOT-ETOTC
!!$IF(ETOT.LT.-3.145D0) THEN
!!$  PRINT*,'FILE RHOWC.DAT WRITTEN'
!!$  CALL ATOMLIB_WRITEPHI('RHOWC1.DAT',GID,NR,LMRX,RHOWC(:,:,1))
!!$  CALL ATOMLIB_WRITEPHI('RHOWC2.DAT',GID,NR,LMRX,RHOWC(:,:,2))
!!$  CALL ATOMLIB_WRITEPHI('RHOWC3.DAT',GID,NR,LMRX,RHOWC(:,:,3))
!!$  CALL ATOMLIB_WRITEPHI('RHOWC4.DAT',GID,NR,LMRX,RHOWC(:,:,4))
!!$END IF

!!$IMETHOD=0
!!$!IMETHOD=1
!!$      IF(IMETHOD.EQ.1) THEN
!!$!       == COLLINEAR METHOD WITH COLLINEAR DENSITY
!!$        ALLOCATE(RHOTEST(NR,LMRX,2))
!!$        ALLOCATE(POTTEST(NR,LMRX,2))
!!$        POTTEST(:,:,1)=0.D0
!!$        RHOTEST(:,:,1)=RHO(:,:,1)
!!$        RHOTEST(:,:,2)=RHO(:,:,4)
!!$        CALL AUGMENTATION_XC(GID,NR,LMRX,2,RHOTEST,ETOT,POTTEST)
!!$        POT(:,:,:)=0.D0
!!$        POT(:,:,1)=POTTEST(:,:,1)
!!$        POT(:,:,4)=POTTEST(:,:,2)
!!$        DEALLOCATE(RHOTEST)
!!$        DEALLOCATE(POTTEST)
!!$!
!!$!      ELSE IF(IMETHOD.EQ.2) THEN
!!$!       == NONCOLLINEAR METHOD WITH COLLINEAR DENSITY ==========================
!!$        ALLOCATE(RHOTEST2(NR,LMRX,NDIMD))
!!$        ALLOCATE(POTTEST2(NR,LMRX,NDIMD))
!!$        RHOTEST2(:,:,:)=0.D0
!!$        POTTEST2(:,:,:)=0.D0
!!$        RHOTEST2(:,:,1)=RHO(:,:,1)
!!$        RHOTEST2(:,:,4)=RHO(:,:,4)
!!$        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOTEST2,ETOT2,POTTEST2)
!!$PRINT*,'ETOT2',ETOT2,ETOT
!!$PRINT*,'LDAPLUSUTEST',ETOT2-ETOT,MAXVAL(ABS(POTTEST2-POT)),MAXLOC(ABS(POTTEST2-POT))
!!$!        ETOT=ETOT2
!!$!        POT(:,:,:)=POTTEST2(:,:,:)
!!$        DEALLOCATE(RHOTEST2)
!!$        DEALLOCATE(POTTEST2)
!!$!
!!$      ELSE IF(IMETHOD.EQ.3) THEN
!!$!       == COMPARISON ==========================================================
!!$
!!$      END IF
        CALL DFT$SETL4('XCONLY',.FALSE.)
!!$PRINT*,'EDFT: EXC ',ETOT
!!$!
!!$!     ==========================================================================
!!$!     == HARTREE ENERGY AND POTENTIAL ==========================================
!!$!     == CORE CONTRIBUTION IS NOT INCLUDED BECAUSE IT IS NOT REPRESENTED IN   ==
!!$!     == THE U-TENSOR AND ONLY THE EXCHANGE PART OF THE CORE-VALENCE IS INCLUDED
!!$!     ==========================================================================
!!$      EDENSITY=0.D0
!!$      DO LM=1,LMRX
!!$        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
!!$        CALL RADIAL$POISSON(GID,NR,L,RHO(:,LM,1),AUX)
!!$        POT(:,LM,1)=POT(:,LM,1)+AUX(:)
!!$        EDENSITY(:)=EDENSITY(:)+0.5D0*AUX(:)*RHO(:,LM,1)
!!$      ENDDO
!!$      CALL RADIAL$R(GID,NR,R)
!!$      EDENSITY=EDENSITY*R(:)**2
!!$      CALL RADIAL$INTEGRAL(GID,NR,EDENSITY,SVAR)
!!$PRINT*,'EDFT: EH ',SVAR
!!$      ETOT=ETOT+SVAR
!
!       ========================================================================
!       ==  CALCULATE HAMILTONIAN IN TOTAL/SPIN REPRESENTATION                ==
!       ========================================================================
        CALL LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
        DEALLOCATE(POT)
!
!     ==========================================================================
!     ==  TRANSFORM HAMILTONIAN FROM TOTAL/SPIN TO UP/DOWN                    ==
!     ==========================================================================
        HAM=REAL(HAM1)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_RADXC(GID,NR,LMRX,NDIMD,RHOIN,FXC,VXC)
!     **************************************************************************
!     ** (SLIGHTLUY CHANGED AUGMENTATION_XC)                                  **
!     **                                                                      **
!     **  CALCULATES THE EXCHANGE AND CORRELATION ENERGY                      **
!     **  FOR A DENSITY GIVEN ON A RADIAL LOGARITHMIC GRID                    **
!     **  TIMES REAL SPHERICAL HARMONICS                                      **
!     **                                                                      **
!     **  THE TOTAL ENERGY IS AN EXPANSION ABOUT THE                          **
!     **  SPHERICAL CONTRIBUTION OF THE DENSITY UP TO QUADRATIC               **
!     **  ORDER IN THE NON-SPHERICAL CONTRIBUTIONS                            **
!     **                                                                      **
!     **  EXC = EXC(XVAL(L=0)*Y0)                                             **
!     **      + 0.5 * D2[EXC]/D[XVAL(L=0)*Y0]**2 * XVAL(L>0)**2               **
!     **                                                                      **
!     **  WHERE XVAL=(/RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS/)              **
!     **  IS AN SPHERICAL HARMONICS EXPANSION ON THE RADIAL GRID.             **
!     **                                                                      **
!     **  DEPENDECIES:                                                        **
!     **    DFT                                                               **
!     **    TIMING                                                            **
!     **    TRACE                                                             **
!     **                                                                      **
!     **  REMARKS: THE GRADIENTS ARE CORRECT ONLY IF DFT SUPPORTS             **
!     **    THIRD DERIVATIVES OF THE XC-ENERGY                                **
!     **   - WHEN USING SELFTEST ON THIS ROUTINE, THEN                        **
!     **     D(EXC)/DRHO(I)=POT(I)*DEX*R(I)**3                                **
!     **     AND THE VALUES AT LARGE RADII MUST BE SURPRESSED                 **
!     **                                                                      **
!     **  REMARK: FOR A COLLINEAR DENSITY THE ROUTINE GIVES DIFFERENT RESULTS **
!     **          WITH NDIMD=2 AND NDIMD=4 DUE TO THE TAYLOR EXPANSION IN     **
!     **          ANGULAR MOMENTUM EXPANSIONS                                 **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1996 ************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TNS=.TRUE. ! NON-SPHERICAL CONTRIBUTIONS ON
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD      ! CAN BE 1,2,4
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(OUT):: FXC(NR)
      REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)
      LOGICAL(4)            :: TGRA   ! SWITCH FOR GRADIENT CORRECTION
      INTEGER(4)            :: NSPIN
      REAL(8)               :: EXC1
      REAL(8)               :: R(NR)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: GRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VGRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: B(:,:)    
      REAL(8)   ,ALLOCATABLE:: C(:,:)    
      REAL(8)               :: VAL5(5),VXC5(5),V2XC5(5,5),V3XC5(5,5,5)
      REAL(8)               :: XVAL(NR,5,LMRX)
      REAL(8)               :: XDER(NR,5,LMRX)
      REAL(8)               :: PI,FOURPI
      REAL(8)               :: Y0
      INTEGER(4)            :: IR,L,II,ISPIN,ISPIN1,ISPIN2,I,J
      INTEGER(4)            :: LM
      INTEGER(4)            :: IMAX
      REAL(8)               :: FAC
      REAL(8)               :: CG0LL
      REAL(8)               :: WORK(NR)
      REAL(8)               :: WORK1(NR)
      REAL(8)               :: WORK2(NR)
!     **************************************************************************
      CALL TRACE$PUSH('AUGMENTATION_XC')
      FXC(:)=0.D0
      VXC(:,:,:)=0.D0
!
!     ==========================================================================
!     ==   CALCULATE SOME CONSTANTS NEEDED LATER                              ==
!     ==========================================================================
      CALL DFT$GETL4('GC',TGRA)
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CG0LL=Y0
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  OBTAIN SPIN DENSITY                                                 ==
!     ==========================================================================
      NSPIN=1
      IF(NDIMD.GT.1) NSPIN=2
      ALLOCATE(RHO(NR,LMRX,NSPIN))
      ALLOCATE(GRHO(NR,LMRX,NSPIN))
      ALLOCATE(VRHO(NR,LMRX,NSPIN))
      RHO(:,:,1)=RHOIN(:,:,1)
      IF(NDIMD.EQ.2) THEN
        RHO(:,:,2)=RHOIN(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
!       == HERE WE NEED TO CALCULATE THE ABSOLUTE VALUE OF THE SPIN DENSITY ====
!       == IN AN ANGULAR MOMENTUM EXPANSION. THIS IS ONLY POSSIBLE APPROXIMATELY
!       == USING A TAYLOR EXPANSION ABOUT THE SPHERICAL PART OF THE SQUARE =====
!       == OF THE SPIN DENSITY. ================================================
        VRHO(:,:,:)=0.D0
        CALL AUGMENTATION_NCOLLTRANS(GID,'RHO',NR,LMRX,RHOIN,RHO,VRHO,VXC)
      END IF
!
!     == IMAX ALLOWS TO RESTRICT SOME LOOPS (1:5) TO (1:IMAX)
      IF(TGRA) THEN
        IF(NSPIN.EQ.2) THEN; IMAX=5; ELSE; IMAX=3; END IF
      ELSE 
        IF(NSPIN.EQ.2) THEN; IMAX=2; ELSE; IMAX=1; END IF
      END IF
!
!     ==========================================================================
!     ==  CALCULATE RADIAL GRADIENT OF THE DENSITY                            ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE GRADIENTS')
      IF(TGRA) THEN
        GRHO(:,:,:)=0.D0
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            CALL RADIAL$DERIVE(GID,NR,RHO(:,LM,ISPIN),GRHO(:,LM,ISPIN))
          ENDDO
        ENDDO
      ELSE
        GRHO(:,:,:)=0.D0
      END IF
!
!     ==========================================================================
!     ==  DEFINE VECTOR (RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS)             ==
!     ==========================================================================
      XVAL(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
          XVAL(:,ISPIN,LM)=RHO(:,LM,ISPIN)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN          ! THIS LOOP PUTS T,T->3; S,S->4 ;T,S->5
          DO ISPIN2=ISPIN1,1,-1    ! AND ASSURES CONSISTENCY WITH NSPIN
            II=II+1
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
              XVAL(:,II,1)=XVAL(:,II,1) &
        &         +CG0LL*(GRHO(:,LM,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                +FAC*RHO(:,LM,ISPIN1)*RHO(:,LM,ISPIN2)/R(:)**2)
            ENDDO
            DO LM=2,LMRX 
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              XVAL(:,II,LM)=XVAL(:,II,LM) &
        &         +0.5D0*CG0LL*(GRHO(:,1,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                      +GRHO(:,LM,ISPIN1)*GRHO(:,1,ISPIN2))
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  CALCULATE EXCHANGE ENERGY FOR THE SPHERICAL DENSITY                 ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE DFT')
      WORK1(:)=0.D0
      XDER(:,:,:)=0.D0
      DO IR=1,NR
!       ==  CYCLE IF THE TOTAL DENSITY VANISHES ================================
        IF(XVAL(IR,1,1).LE.0.D0) CYCLE
!       == NOW CALL DFT ROUTINE ================================================
        VAL5(:)=XVAL(IR,:,1)*Y0
        CALL DFT3(VAL5,EXC1,VXC5,V2XC5,V3XC5)
!       == NOW CALCULATE ENERGY DENSITY AND DERIAVTIVES ========================
        FXC(IR)=EXC1
        XDER(IR,:,1)  =FOURPI*VXC5(:)*Y0
        DO LM=2,LMRX
          DO I=1,IMAX        ! IMAX=<5 
            DO J=1,IMAX
              WORK1(IR)=WORK1(IR) &
       &              +0.5D0*XVAL(IR,I,LM)*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,:,1)=XDER(IR,:,1) &
       &              +0.5D0*Y0*XVAL(IR,I,LM)*V3XC5(:,I,J)*XVAL(IR,J,LM)
              XDER(IR,I,LM)=XDER(IR,I,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,J,LM)=XDER(IR,J,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,I,LM)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  TRANSFORM POTENTIALS FOR SPHERICAL PART                             ==
!     ==========================================================================
      ALLOCATE(VGRHO(NR,LMRX,NSPIN))
      VRHO(:,:,:)=0.D0
      VGRHO(:,:,:)=0.D0
      DO ISPIN=1,NSPIN
        DO LM=1,LMRX
          VRHO(:,LM,ISPIN)=XDER(:,ISPIN,LM)
        ENDDO
      ENDDO
      IF(TGRA) THEN
        II=2
        DO ISPIN1=1,NSPIN
          DO ISPIN2=ISPIN1,1,-1
            II=II+1
!           == FIRST RESOLVE XVAL(:,II,1) ======================================
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=DBLE(L*(L+1))
!             == THE FOLLOWING LINES DIVIDE BY ZERO IF R=0. ====================
!             == THE VALUES ARE OVERWRITTEN AT THE END OF THE ROUTINE... =======
              VRHO(:,LM,ISPIN1)  =VRHO(:,LM,ISPIN1) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN2)
              VRHO(:,LM,ISPIN2)  =VRHO(:,LM,ISPIN2) &
      &                 +CG0LL*FAC/R(:)**2*XDER(:,II,1)*RHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN1) =VGRHO(:,LM,ISPIN1) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN2)
              VGRHO(:,LM,ISPIN2) =VGRHO(:,LM,ISPIN2) &
      &                 +CG0LL*XDER(:,II,1)*GRHO(:,LM,ISPIN1)
            ENDDO
!           == NOW RESOLVE XVAL(:,II,LM) =======================================
            DO LM=2,LMRX
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              VGRHO(:,1,ISPIN1) =VGRHO(:,1,ISPIN1) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN2)
              VGRHO(:,1,ISPIN2) =VGRHO(:,1,ISPIN2) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,LM,ISPIN1)
              VGRHO(:,LM,ISPIN2)=VGRHO(:,LM,ISPIN2) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN1)
              VGRHO(:,LM,ISPIN1)=VGRHO(:,LM,ISPIN1) &
      &                 +0.5D0*CG0LL*XDER(:,II,LM)*GRHO(:,1,ISPIN2)
            ENDDO
          ENDDO
        ENDDO               
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==  V = V -1/R**2 D/DR [ R**2 VGRHO ]                                   ==
!     ==  V = V -[2/R VGRHO+ D/DR VGRHO ]                                     ==
!     ==========================================================================
      IF(TGRA) THEN
        DO ISPIN=1,NSPIN
          DO LM=1,LMRX
            IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
!           == FIRST ALTERNATIVE
!           CALL RADIAL$DERIVE(GID,NR,VGRHO(:,LM,ISPIN),WORK2)   !NOT SO 
!           WORK1(:)=2.D0/R(:)*VGRHO(:,LM,ISPIN)+WORK2(:)  !GOOD
!           ==  SECOND ALTERNATIVE APPEARS TO BE MORE ACCURATE
            WORK2(:)=VGRHO(:,LM,ISPIN)*R(:)**2
            CALL RADIAL$DERIVE(GID,NR,WORK2,WORK1)
            WORK1(2:)=WORK1(2:)/R(2:)**2
            WORK(1)=WORK(2)  ! AVOID DIVIDE BY ZERO
!           == ALTERNATIVES FINISHED
            VRHO(:,LM,ISPIN)=VRHO(:,LM,ISPIN)-WORK1(:)
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(VGRHO)
!
!     ==========================================================================
!     ==   CORRECT FOR DIVERGENCE AT THE ORIGIN:                              ==
!     ==   IF A SHIFTED LOGARITHMIC GRID IS USED THE FIRST GRID POINT         ==
!     ==   IS MESSED UP BECAUSE OF FACTORS 1/R                                ==
!     ==========================================================================
      IF(R(1).LT.1.D-5) THEN
        VRHO(1,:,:)=VRHO(2,:,:)
      END IF
!
!     ==========================================================================
!     ==  TRANSFORM GRADIENT POTENTIAL BACK TO POTENTIALS                     ==
!     ==========================================================================
      VXC(:,:,1)=VRHO(:,:,1)
      IF(NDIMD.EQ.2) THEN
        VXC(:,:,2)=VRHO(:,:,2)
      ELSE IF(NDIMD.EQ.4) THEN
        CALL AUGMENTATION_NCOLLTRANS(GID,'POT',NR,LMRX,RHOIN,RHO,VRHO,VXC)
      END IF     
      DEALLOCATE(RHO)
      DEALLOCATE(GRHO)
      DEALLOCATE(VRHO)
                      CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SIMPLEDC_OLD(GID,NR,LMNX,LNX,LOX,CHI,LRX,AECORE &
     &                        ,DENMAT,ETOT,HAM)
!     **************************************************************************
!     **  DOUBLE COUNTING CORRECTION FOR THE HYBRID FUNCTIONAL                **
!     **                                                                      **
!     **  DETERMINES THE HARTREE AND EXCHANGE-ONLY ENERGY FROM THE            **
!     **  DFT FUNCTIONAL                                                      **
!     **  FOR THE DENSITY BUILT FROM THE LOCAL ORBITALS AND THE CORE DENSITY  **
!     **  THIS ENERGY NEEDS TO BE SUBTRACTED FROM THE TOTAL ENERGY            **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: LRX
      INTEGER(4)  ,INTENT(IN) :: LMNX       ! #(LOCAL ORBITALS)
      INTEGER(4)  ,INTENT(IN) :: LNX        ! #(RADIAL FUNCTIONS)
      INTEGER(4)  ,INTENT(IN) :: LOX(LNX)   !MAIN ANGULAR MOMENTUM OF LOCAL ORB.
      REAL(8)     ,INTENT(IN) :: CHI(NR,LNX)
      REAL(8)     ,INTENT(IN) :: AECORE(NR)
      REAL(8)     ,INTENT(IN) :: DENMAT(LMNX,LMNX,4) ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT       ! DOUBLE COUNTINNG ENERGY
      REAL(8)     ,INTENT(OUT):: HAM(LMNX,LMNX,4)  ! DETOT/D(RHO^*)        
      INTEGER(4)  ,PARAMETER  :: METHOD=1
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)  ,PARAMETER  :: NDIMD=4
      COMPLEX(8)              :: DENMAT1(LMNX,LMNX,NDIMD)
      COMPLEX(8)              :: HAM1(LMNX,LMNX,NDIMD)
      REAL(8)                 :: R(NR)
      REAL(8)     ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)     ,ALLOCATABLE:: RHO2(:,:,:)
      REAL(8)     ,ALLOCATABLE:: POT2(:,:,:)
      REAL(8)     ,ALLOCATABLE:: RHOWC(:,:,:)
      REAL(8)     ,ALLOCATABLE:: POT(:,:,:)
      REAL(8)                 :: EDENSITY(NR)
      REAL(8)                 :: AUX(NR),SVAR
      REAL(8)                 :: PI,Y0
      INTEGER(4)              :: LMRX,L
      INTEGER(4)              :: IDIM,LM,LMN
      REAL(8)                 :: ETOTC,ETOTV
INTEGER(4) :: LMRX1
INTEGER(4) :: IMETHOD
 REAL(8)     ,ALLOCATABLE:: RHOTEST(:,:,:)
 REAL(8)     ,ALLOCATABLE:: POTTEST(:,:,:)
 REAL(8)     ,ALLOCATABLE:: RHOTEST2(:,:,:)
 REAL(8)     ,ALLOCATABLE:: POTTEST2(:,:,:)
 REAL(8)                 :: ETOT2
!     **************************************************************************
      LMRX=(LRX+1)**2
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      ETOT=0.D0
!
!     ==========================================================================
!     ==  TRANSFORM DENSITY MATRIX FROM UP/DOWN TO TOTAL/SPIN                 ==
!     ==========================================================================
      DENMAT1=CMPLX(DENMAT)
!
!     ==========================================================================
!     ==  CALCULATE DENSITY                                                   ==
!     ==========================================================================
      ALLOCATE(RHO(NR,LMRX,NDIMD))
      DO IDIM=1,NDIMD
        CALL AUGMENTATION_RHO(NR,LNX,LOX,CHI &
     &                       ,LMNX,DENMAT1(:,:,IDIM),LMRX,RHO(:,:,IDIM))
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE ENERGY AND POTENTIAL                                      ==
!     ==========================================================================
      IF(METHOD.EQ.1) THEN 
!       == THIS IS THE SAME METHOD AS 1000 JUST WITHOUT COMMENTS ===============
        ALLOCATE(RHOWC(NR,LMRX,NDIMD))  !WITH CORE
        RHOWC=RHO
        RHOWC(:,1,1)=RHO(:,1,1)+AECORE(:)
        ALLOCATE(POT(NR,LMRX,NDIMD))
        CALL DFT$SETL4('XCONLY',.TRUE.)
        CALL AUGMENTATION_XC(GID,NR,1,1,AECORE,ETOTC,POT)
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,ETOTV,POT)
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOWC,ETOT,POT)
        CALL DFT$SETL4('XCONLY',.FALSE.)
        ETOT=ETOT-ETOTC
        CALL LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
        DEALLOCATE(POT)
        HAM=REAL(HAM1)
!
      ELSE IF(METHOD.EQ.1000) THEN
        ALLOCATE(RHOWC(NR,LMRX,NDIMD))  !WITH CORE
        RHOWC=RHO
        RHOWC(:,1,1)=RHO(:,1,1)+AECORE(:)
!
!     ==========================================================================
!     ==  CALCULATE ENERGY AND POTENTIAL                                      ==
!     ==========================================================================
!     == EXCHANGE ENERGY AND POTENTIAL =========================================
        CALL DFT$SETL4('XCONLY',.TRUE.)
!
!     ==========================================================================
!     == THIS FORMULATION IS BASED ON A NONCOLLINEAR FORMULATION, WHICH       ==
!     == YIELDS DIFFERENT RESULTS FROM A COLLINEAR FORMULATION EVEN FOR       ==
!     == A COLLINEAR DENSITY                                                  ==
!     ==                                                                      ==
!     == THE REASON FOR THIS DIFFERENCE IS THE TRANSFORMATION OF A            ==
!     == NON-COLLINEAR DENSITY WITHIN AUGMENTATION_NCOLLTRANS WHICH IS CALLED ==
!     == BY AUGMENTATION_XC                                                   ==
!     ==                                                                      ==
!     ==========================================================================
        ALLOCATE(POT(NR,LMRX,NDIMD))
        CALL AUGMENTATION_XC(GID,NR,1,1,AECORE,ETOTC,POT)
!!$ALLOCATE(RHO2(NR,LMRX,2))
!!$ALLOCATE(POT2(NR,LMRX,2))
!!$RHO2(:,:,1)=RHO(:,:,1)
!!$RHO2(:,:,2)=RHO(:,:,4)
!!$CALL AUGMENTATION_XC(GID,NR,LMRX,2,RHO2,ETOTV,POT2)
!!$PRINT*,'ETOTV COLLINEAR',ETOTV
!!$CALL AUGMENTATION_WRITEPHI('RHO4_Z.DAT',GID,NR,LMRX,RHO4(:,:,4))
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,ETOTV,POT)
!!$PRINT*,'ETOTV NONCOLLINEAR',ETOTV
!!$PRINT*,'GID,NR,LMRX,NDIMD ',GID,NR,LMRX,NDIMD
!!$RHO2(:,:,1)=RHOWC(:,:,1)
!!$RHO2(:,:,2)=RHOWC(:,:,4)
!!$CALL ATOMLIB_WRITEPHI('RHO2WC_T.DAT',GID,NR,LMRX,RHO2(:,:,1))
!!$CALL ATOMLIB_WRITEPHI('RHO2WC_S.DAT',GID,NR,LMRX,RHO2(:,:,2))
!!$CALL AUGMENTATION_XC(GID,NR,LMRX,2,RHO2,ETOT,POT2)
!!$CALL ATOMLIB_WRITEPHI('POT2WC_T.DAT',GID,NR,LMRX,POT2(:,:,1))
!!$CALL ATOMLIB_WRITEPHI('POT2WC_S.DAT',GID,NR,LMRX,POT2(:,:,2))
!!$PRINT*,'ETOT COLLINEAR',ETOTV
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOWC,ETOT,POT)
!!$CALL ATOMLIB_WRITEPHI('POT4WC_T.DAT',GID,NR,LMRX,POT(:,:,1))
!!$CALL ATOMLIB_WRITEPHI('POT4WC_X.DAT',GID,NR,LMRX,POT(:,:,2))
!!$CALL ATOMLIB_WRITEPHI('POT4WC_Y.DAT',GID,NR,LMRX,POT(:,:,3))
!!$CALL ATOMLIB_WRITEPHI('POT4WC_Z.DAT',GID,NR,LMRX,POT(:,:,4))
!!$PRINT*,'ETOTV NONCOLLINEAR',ETOTV
!POT(:,:,:)=0.D0
!POT(:,:,1)=POT2(:,:,1)
!POT(:,:,4)=POT2(:,:,2)
!!$DEALLOCATE(RHO2)
!!$DEALLOCATE(POT2)
!!$PRINT*,'TOTAL        EXCHANGE ENERGY (LOCAL) ',ETOT
!!$PRINT*,'VALENCE      EXCHANGE ENERGY (LOCAL) ',ETOTV
!!$PRINT*,'CORE         EXCHANGE ENERGY (LOCAL) ',ETOTC
!!$PRINT*,'CORE-VALENCE EXCHANGE ENERGY (LOCAL) ',ETOT-ETOTV-ETOTC
        ETOT=ETOT-ETOTC
!!$IF(ETOT.LT.-3.145D0) THEN
!!$  PRINT*,'FILE RHOWC.DAT WRITTEN'
!!$  CALL ATOMLIB_WRITEPHI('RHOWC1.DAT',GID,NR,LMRX,RHOWC(:,:,1))
!!$  CALL ATOMLIB_WRITEPHI('RHOWC2.DAT',GID,NR,LMRX,RHOWC(:,:,2))
!!$  CALL ATOMLIB_WRITEPHI('RHOWC3.DAT',GID,NR,LMRX,RHOWC(:,:,3))
!!$  CALL ATOMLIB_WRITEPHI('RHOWC4.DAT',GID,NR,LMRX,RHOWC(:,:,4))
!!$END IF

!!$IMETHOD=0
!!$!IMETHOD=1
!!$      IF(IMETHOD.EQ.1) THEN
!!$!       == COLLINEAR METHOD WITH COLLINEAR DENSITY
!!$        ALLOCATE(RHOTEST(NR,LMRX,2))
!!$        ALLOCATE(POTTEST(NR,LMRX,2))
!!$        POTTEST(:,:,1)=0.D0
!!$        RHOTEST(:,:,1)=RHO(:,:,1)
!!$        RHOTEST(:,:,2)=RHO(:,:,4)
!!$        CALL AUGMENTATION_XC(GID,NR,LMRX,2,RHOTEST,ETOT,POTTEST)
!!$        POT(:,:,:)=0.D0
!!$        POT(:,:,1)=POTTEST(:,:,1)
!!$        POT(:,:,4)=POTTEST(:,:,2)
!!$        DEALLOCATE(RHOTEST)
!!$        DEALLOCATE(POTTEST)
!!$!
!!$!      ELSE IF(IMETHOD.EQ.2) THEN
!!$!       == NONCOLLINEAR METHOD WITH COLLINEAR DENSITY ==========================
!!$        ALLOCATE(RHOTEST2(NR,LMRX,NDIMD))
!!$        ALLOCATE(POTTEST2(NR,LMRX,NDIMD))
!!$        RHOTEST2(:,:,:)=0.D0
!!$        POTTEST2(:,:,:)=0.D0
!!$        RHOTEST2(:,:,1)=RHO(:,:,1)
!!$        RHOTEST2(:,:,4)=RHO(:,:,4)
!!$        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOTEST2,ETOT2,POTTEST2)
!!$PRINT*,'ETOT2',ETOT2,ETOT
!!$PRINT*,'LDAPLUSUTEST',ETOT2-ETOT,MAXVAL(ABS(POTTEST2-POT)),MAXLOC(ABS(POTTEST2-POT))
!!$!        ETOT=ETOT2
!!$!        POT(:,:,:)=POTTEST2(:,:,:)
!!$        DEALLOCATE(RHOTEST2)
!!$        DEALLOCATE(POTTEST2)
!!$!
!!$      ELSE IF(IMETHOD.EQ.3) THEN
!!$!       == COMPARISON ==========================================================
!!$
!!$      END IF
        CALL DFT$SETL4('XCONLY',.FALSE.)
!!$PRINT*,'EDFT: EXC ',ETOT
!!$!
!!$!     ==========================================================================
!!$!     == HARTREE ENERGY AND POTENTIAL ==========================================
!!$!     == CORE CONTRIBUTION IS NOT INCLUDED BECAUSE IT IS NOT REPRESENTED IN   ==
!!$!     == THE U-TENSOR AND ONLY THE EXCHANGE PART OF THE CORE-VALENCE IS INCLUDED
!!$!     ==========================================================================
!!$      EDENSITY=0.D0
!!$      DO LM=1,LMRX
!!$        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
!!$        CALL RADIAL$POISSON(GID,NR,L,RHO(:,LM,1),AUX)
!!$        POT(:,LM,1)=POT(:,LM,1)+AUX(:)
!!$        EDENSITY(:)=EDENSITY(:)+0.5D0*AUX(:)*RHO(:,LM,1)
!!$      ENDDO
!!$      CALL RADIAL$R(GID,NR,R)
!!$      EDENSITY=EDENSITY*R(:)**2
!!$      CALL RADIAL$INTEGRAL(GID,NR,EDENSITY,SVAR)
!!$PRINT*,'EDFT: EH ',SVAR
!!$      ETOT=ETOT+SVAR
!
!       ========================================================================
!       ==  CALCULATE HAMILTONIAN IN TOTAL/SPIN REPRESENTATION                ==
!       ========================================================================
        CALL LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
        DEALLOCATE(POT)
!
!     ==========================================================================
!     ==  TRANSFORM HAMILTONIAN FROM TOTAL/SPIN TO UP/DOWN                    ==
!     ==========================================================================
        HAM=REAL(HAM1)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEX(EX)
!     ************************************************************************** 
!     **  WORK OUT THE ENERGY USING THE LOCAL APPROXIMATION                   **
!     **  TAILED PARTIAL WAVES                                                **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR,DENMAT,HAMIL,NSP
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: EX
      TYPE SPECIAL_TYPE
        INTEGER(4)      :: NIJK
        REAL(8),POINTER :: ARHO(:,:,:)
        REAL(8),POINTER :: APOT(:,:,:)
        REAL(8),POINTER :: ASINGLE(:,:)
        REAL(8),POINTER :: ATRIPLE(:,:,:,:)
      END TYPE SPECIAL_TYPE
      TYPE(SPECIAL_TYPE),ALLOCATABLE :: SPECIAL(:)
      REAL(8)                :: RBAS(3,3)
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      INTEGER(4)             :: ISP
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: NE
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4)             :: NX
      INTEGER(4)             :: NPOW
      INTEGER(4)             :: NIJK
      INTEGER(4)             :: NND
      INTEGER(4)             :: NN,NNA,NNB
      INTEGER(4)             :: NS,NP,NT !#(SINGLE, DOUBLE, TRIPLE TERMS)
      INTEGER(4)             :: IAT,IATA,IATB
      INTEGER(4)             :: ISPA,ISPB
      REAL(8)                :: RA(3),RB(3)
      INTEGER(4)             :: NIJKA,NIJKB
      INTEGER(4)             :: NEA,NEB
      INTEGER(4)             :: LMN1A,LMN2A,LMN1B,LMN2B,LMN3A,LMN3B
      INTEGER(4)             :: LMNXA,LMNXB
      INTEGER(4)             :: NA,NB
      INTEGER(4)             :: NDIMD
      REAL(8)                :: SVAR
      REAL(8)   ,ALLOCATABLE :: SAB(:,:)
REAL(8)   ,ALLOCATABLE :: ZA1(:)
      REAL(8)   ,ALLOCATABLE :: ZA(:)
      REAL(8)   ,ALLOCATABLE :: ZB(:)
      REAL(8)   ,ALLOCATABLE :: D(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DB(:,:,:)
      REAL(8)   ,ALLOCATABLE :: H(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HB(:,:,:)
      INTEGER(4),ALLOCATABLE :: NNAT(:)    !(NAT) POINTER TO ONSITE TERMS
INTEGER(4) :: I,J,K,L,I1,J1,K1
REAL(8) :: SVAR1,SVAR2
INTEGER(4) :: IX,IND1,IND2,IND3
REAL(8)    :: XDELTA,XSVAR
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_OFFSITEX')
PRINT*,'============ OFFSITEX ============================='
!!$CALL LMTO$REPORTHAMIL(6)
!
!     ==========================================================================
!     == MAP PRODUCTS TO GAUSSIANS                                            ==
!     ==========================================================================
      ALLOCATE(SPECIAL(NSP))
      DO ISP=1,NSP
        NE=POTPAR(ISP)%TAILED%PRODRHO%NE
        LNX=POTPAR(ISP)%TAILED%LNX
        LMNX=POTPAR(ISP)%TAILED%LMNX
        ALLOCATE(LOX(LNX))
        LOX(:)=POTPAR(ISP)%TAILED%LOX(:)
        NPOW=POTPAR(ISP)%TAILED%PRODRHO%NIJK  !#(POWERS)
        NX=2*(NPOW-1)                         !NIJK IS HERE THE HIGHEST POWER+1
        NIJK=(NX+1)*(NX+2)*(NX+3)/6
        NS=POTPAR(ISP)%TAILED%SINGLE%NORB
        NP=POTPAR(ISP)%TAILED%PRODRHO%NORB
        NT=POTPAR(ISP)%TAILED%TRIPLE%NORB
        ALLOCATE(SPECIAL(ISP)%ARHO(NIJK*NE,LMNX,LMNX))
        ALLOCATE(SPECIAL(ISP)%APOT(NIJK*NE,LMNX,LMNX))
        ALLOCATE(SPECIAL(ISP)%ASINGLE(NIJK*NE,LMNX))
        ALLOCATE(SPECIAL(ISP)%ATRIPLE(NIJK*NE,LMNX,LMNX,LMNX))
        SPECIAL(ISP)%NIJK=NIJK
        CALL LMTO_EXPANDPRODS(NX,NIJK,NE,LNX,LOX,LMNX,NPOW,NS,NP,NT &
     &              ,POTPAR(ISP)%TAILED%PRODRHO%C,POTPAR(ISP)%TAILED%PRODPOT%C &
     &              ,POTPAR(ISP)%TAILED%SINGLE%C,POTPAR(ISP)%TAILED%TRIPLE%C &
     &              ,SPECIAL(ISP)%ARHO,SPECIAL(ISP)%APOT &
     &              ,SPECIAL(ISP)%ASINGLE,SPECIAL(ISP)%ATRIPLE)
        DEALLOCATE(LOX)

!!$!TEST U-TENSOR ==============================================
!!$IF(ISP.EQ.1) THEN
!!$NEA=POTPAR(ISP)%TAILED%PRODRHO%NE
!!$NEB=POTPAR(ISP)%TAILED%PRODRHO%NE
!!$NIJKA=SPECIAL(ISP)%NIJK
!!$NIJKB=SPECIAL(ISP)%NIJK
!!$RA=0.D0
!!$RB=0.D0
!!$ALLOCATE(SAB(NIJKA*NEA,NIJKB*NEB))
!!$CALL GAUSSIAN$GAUSSOVERLAP(NIJKA,NEA,POTPAR(ISP)%TAILED%PRODRHO%E,RA &
!!$&                         ,NIJKB,NEB,POTPAR(ISP)%TAILED%PRODPOT%E,RB &
!!$&                         ,SAB)
!!$PRINT*,' E ',POTPAR(ISP)%TAILED%PRODRHO%E(:)
!!$DO I=1,NIJKA
!!$  CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',I,I1,J1,K1)
!!$  WRITE(*,FMT='(3I5,"  ",50E10.3)')I1,J1,K1,(SPECIAL(ISP)%ASINGLE(NIJK*(J-1)+I,1),J=1,NE)
!!$ENDDO
!!$STOP 'FORCED1'
!!$ALLOCATE(ZA(NIJKA*NEA))
!!$ALLOCATE(ZB(NIJKA*NEA))
!!$OPEN(2351,FILE='UNEU1.DAT')
!!$OPEN(2352,FILE='UNEU2.DAT')
!!$REWIND 2351
!!$REWIND 2352
!!$DO I=1,LMNX
!!$  ZB(:)=MATMUL(SAB,SPECIAL(ISP)%ASINGLE(:,I))
!!$  DO J=1,LMNX
!!$    ZA(:)=MATMUL(SAB,SPECIAL(ISP)%APOT(:,I,J))
!!$    DO K=1,LMNX
!!$      DO L=1,LMNX
!!$        SVAR1=DOT_PRODUCT(SPECIAL(ISP)%ARHO(:,K,L),ZA)
!!$        WRITE(2351,FMT='(5I5,E20.8)')ISP,I,J,K,L,SVAR1
!!$        SVAR2=DOT_PRODUCT(SPECIAL(ISP)%ATRIPLE(:,K,L,J),ZB)
!!$        WRITE(2352,FMT='(5I5,E20.8)')ISP,I,J,K,L,SVAR2
!!$WRITE(*,FMT='(5I5,3E20.10)')ISP,I,J,K,L,SVAR1,SVAR2,SVAR1-SVAR2
!!$STOP 'FORCED'
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
!!$CLOSE(2351)
!!$CLOSE(2352)
!!$DEALLOCATE(ZA)
!!$DEALLOCATE(ZB)
!!$DEALLOCATE(SAB)
!!$STOP 'FORCED'
!!$END IF
!
      ENDDO
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     == COLLECT POINTERS NNAT TO ONSITE DENSITY MATRIX ELEMENTS              ==
!     ==========================================================================
      ALLOCATE(NNAT(NAT))
      NNAT=0
      NND=SIZE(DENMAT)
      DO NN=1,NND
        IATA=DENMAT(NN)%IAT1
        IATB=DENMAT(NN)%IAT2
        IF(IATA.EQ.IATB.AND.SUM(DENMAT(NN)%IT**2).EQ.0) NNAT(IATA)=NN
      ENDDO
      DO IAT=1,NAT
        IF(NNAT(IAT).LE.0) THEN
          CALL ERROR$MSG('NO ONSITE TERMS FOUND FOR ATOM')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$STOP('LMTO_OFFSITEX')
        END IF
      ENDDO
!
!     ==========================================================================
!     == LOOP OVER PAIRS                                                      ==
!     ==========================================================================
      IF(SIZE(HAMIL).NE.NND) THEN
        CALL ERROR$MSG('SIZE OF DENSITY MATRIX AND HAMILTONIAN INONSISTENT')
        CALL ERROR$STOP('LMTO_OFFSITEX')
      END IF
      EX=0.D0
      DO NN=1,NND
        IATA=DENMAT(NN)%IAT1
        IATB=DENMAT(NN)%IAT2
!       == CONSIDER ONLY OFFSITE TERMS =========================================
        IF(IATA.EQ.IATB.AND.SUM(DENMAT(NN)%IT**2).EQ.0) CYCLE
!!$CALL LMTO$REPORTDENMAT(6)
!!$STOP
!!$IND1=1
!!$IND2=1
!!$IND3=1
!!$XSVAR=DENMAT(NN)%MAT(IND1,IND2,IND3)
!!$XDELTA=1.D-2
!!$IX=-3
!!$1000 CONTINUE
!!$IX=IX+1
!!$DENMAT(NN)%MAT(IND1,IND2,IND3)=XSVAR+XDELTA*REAL(IX,KIND=8)
!!$HAMIL(NN)%MAT(IND1,IND2,IND3)=0.D0
!!$EX=0.D0
        ISPA=ISPECIES(IATA)
        ISPB=ISPECIES(IATB)
        RA(:)=R0(:,IATA)
        RB(:)=R0(:,IATB)+MATMUL(RBAS,REAL(DENMAT(NN)%IT(:),KIND=8))
        NNA=NNAT(IATA)
        NNB=NNAT(IATB)
!
!       ========================================================================
!       == CALCULATE OVERLAP OF GAUSSIANS <G_I|G_J>                           ==
!       ========================================================================
        NEA=POTPAR(ISPA)%TAILED%PRODRHO%NE
        NEB=POTPAR(ISPB)%TAILED%PRODRHO%NE
        NIJKA=SPECIAL(ISPA)%NIJK
        NIJKB=SPECIAL(ISPB)%NIJK
        ALLOCATE(SAB(NIJKA*NEA,NIJKB*NEB))
        CALL GAUSSIAN$GAUSSOVERLAP(NIJKA,NEA,POTPAR(ISPA)%TAILED%PRODRHO%E,RA &
    &                             ,NIJKB,NEB,POTPAR(ISPB)%TAILED%PRODPOT%E,RB &
    &                             ,SAB)
!
!       ========================================================================
!       == BLOW UP DENSITY MATRIX                                             ==
!       ========================================================================
        LMNXA=POTPAR(ISPA)%TAILED%LMNX
        LMNXB=POTPAR(ISPB)%TAILED%LMNX
        NDIMD=DENMAT(NN)%N3
        NA=DENMAT(NN)%N1
        NB=DENMAT(NN)%N2
        ALLOCATE(D(LMNXA,LMNXB,NDIMD))
        CALL LMTO_BLOWUPDENMATNL(IATA,IATB,NDIMD &
     &                          ,NA,NB,DENMAT(NN)%MAT,LMNXA,LMNXB,D)
        ALLOCATE(DA(LMNXA,LMNXA,NDIMD))
        CALL LMTO_BLOWUPDENMATNL(IATA,IATA,NDIMD &
     &                          ,NA,NA,DENMAT(NNA)%MAT,LMNXA,LMNXA,DA)
        ALLOCATE(DB(LMNXB,LMNXB,NDIMD))
        CALL LMTO_BLOWUPDENMATNL(IATB,IATB,NDIMD &
     &                          ,NB,NB,DENMAT(NNB)%MAT,LMNXB,LMNXB,DB)
!
!       ========================================================================
!       == ADD UP EXCHANGE ENERGY                                             ==
!       ========================================================================
        ALLOCATE(H(LMNXA,LMNXB,NDIMD))
        ALLOCATE(HA(LMNXA,LMNXA,NDIMD))
        ALLOCATE(HB(LMNXB,LMNXB,NDIMD))
        ALLOCATE(ZA(NIJKA*NEA))
ALLOCATE(ZA1(NIJKA*NEA))
        ALLOCATE(ZB(NIJKB*NEB))
        H=0.D0
        HA=0.D0
        HB=0.D0
        DO LMN1B=1,LMNXB
          DO LMN2B=1,LMNXB
            ZA(:)=MATMUL(SAB,SPECIAL(ISPB)%APOT(:,LMN1B,LMN2B))
ZA1(:)=MATMUL(SAB,SPECIAL(ISPB)%ARHO(:,LMN1B,LMN2B))
            DO LMN1A=1,LMNXA
              DO LMN2A=1,LMNXA
                SVAR=DOT_PRODUCT(SPECIAL(ISPA)%ARHO(:,LMN1A,LMN2A),ZA)
SVAR2= DOT_PRODUCT(SPECIAL(ISPA)%APOT(:,LMN1A,LMN2A),ZA1)
!PRINT*,'TEST',SVAR,SVAR2,SVAR-SVAR2, (SVAR-SVAR2)/SQRT(0.5D0*(SVAR**2+SVAR2**2))
SVAR=0.5D0*(SVAR+SVAR2)
                EX=EX-0.25D0*SVAR*SUM(D(LMN1A,LMN1B,:)*D(LMN2A,LMN2B,:))
                H(LMN1A,LMN1B,:)=H(LMN1A,LMN1B,:)-0.25D0*SVAR*D(LMN2A,LMN2B,:)
                H(LMN2A,LMN2B,:)=H(LMN2A,LMN2B,:)-0.25D0*SVAR*D(LMN1A,LMN1B,:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!!$        DO LMN1B=1,LMNXB
!!$          ZA(:)=MATMUL(SAB,SPECIAL(ISPB)%ASINGLE(:,LMN1B))
!!$          DO LMN1A=1,LMNXA
!!$            DO LMN2A=1,LMNXA 
!!$              DO LMN3A=1,LMNXA
!!$                SVAR=DOT_PRODUCT(SPECIAL(ISPA)%ATRIPLE(:,LMN1A,LMN2A,LMN3A),ZA)
!!$                EX=EX-0.25D0*SVAR*SUM(DA(LMN2A,LMN3A,:)*D(LMN1A,LMN1B,:))
!!$                HA(LMN2A,LMN3A,:)=HA(LMN2A,LMN3A,:)-0.25D0*SVAR*D(LMN1A,LMN1B,:)
!!$                H(LMN1A,LMN1B,:) =H(LMN1A,LMN1B,:)-0.25D0*SVAR*DA(LMN2A,LMN3A,:)
!!$              ENDDO
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO               
!!$        DO LMN1A=1,LMNXA
!!$          ZB(:)=MATMUL(SPECIAL(ISPA)%ASINGLE(:,LMN1A),SAB)
!!$          DO LMN1B=1,LMNXB
!!$            DO LMN2B=1,LMNXB
!!$              DO LMN3B=1,LMNXB
!!$                SVAR=DOT_PRODUCT(ZB,SPECIAL(ISPB)%ATRIPLE(:,LMN1B,LMN2B,LMN3B))
!!$                EX=EX-0.5D0*SVAR*SUM(D(LMN1A,LMN1B,:)*DB(LMN2B,LMN3B,:))
!!$                H(LMN1A,LMN1B,:) =H(LMN1A,LMN1B,:) -0.5D0*SVAR*DB(LMN2B,LMN3B,:)
!!$                HB(LMN2B,LMN3B,:)=HB(LMN2B,LMN3B,:)-0.5D0*SVAR*D(LMN1A,LMN1B,:)
!!$              ENDDO
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO               
DEALLOCATE(ZA1)
        DEALLOCATE(ZA)
        DEALLOCATE(ZB)
        DEALLOCATE(SAB)
!
!       ========================================================================
!       == SHRINK DOWN HAMILTONIAN                                            ==
!       ========================================================================
!       __D,DA,DB IS RE-USED AS WORK ARRAY TO COLLECT HAMILTONIANS______________
        D=0.D0
        DA=0.D0
        DB=0.D0
        CALL LMTO_SHRINKDOWNHTNL(IATA,IATB,NDIMD,LMNXA,LMNXB,H &
     &                                                      ,NA,NB,D(:NA,:NB,:))
        HAMIL(NN)%MAT(:,:,:)=HAMIL(NN)%MAT(:,:,:)+D(:NA,:NB,:)
        CALL LMTO_SHRINKDOWNHTNL(IATA,IATA,NDIMD,LMNXA,LMNXA,HA &
     &                                                     ,NA,NA,DA(:NA,:NA,:))
        HAMIL(NNA)%MAT(:,:,:)=HAMIL(NNA)%MAT(:,:,:)+DA(:NA,:NA,:)
        CALL LMTO_SHRINKDOWNHTNL(IATB,IATB,NDIMD,LMNXB,LMNXB,HB &
     &                                                     ,NB,NB,DB(:NB,:NB,:))
        HAMIL(NNB)%MAT(:,:,:)=HAMIL(NNB)%MAT(:,:,:)+DB(:NB,:NB,:)
        DEALLOCATE(D)
        DEALLOCATE(DA)
        DEALLOCATE(DB)
        DEALLOCATE(H)
        DEALLOCATE(HA)
        DEALLOCATE(HB)
!!$WRITE(*,*)XDELTA*REAL(IX,8),EX,HAMIL(NN)%MAT(IND1,IND2,IND3)
!!$IF(IX.EQ.3) STOP 'FORCED'
!!$GOTO 1000
      ENDDO
!
!     ==========================================================================
!     == CLEAN UP TO AVOID MEMORY LEAK                                        ==
!     ==========================================================================
      DEALLOCATE(NNAT)
      DO ISP=1,NSP
        DEALLOCATE(SPECIAL(ISP)%ARHO)
        DEALLOCATE(SPECIAL(ISP)%APOT)
        DEALLOCATE(SPECIAL(ISP)%ASINGLE)
        DEALLOCATE(SPECIAL(ISP)%ATRIPLE)
      ENDDO
      DEALLOCATE(SPECIAL)
!!$PRINT*,'EX ',EX
!!$STOP
!!$CALL LMTO$REPORTDENMAT(6)
!!$CALL LMTO$REPORTHAMIL(6)
!!$STOP
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_EXPANDPRODS(NX,NIJK,NE,LNX,LOX,LMNX,NPOW,NS,NP,NT &
     &                           ,PRODRHO,PRODPOT,SINGLE,TRIPLE &
     &                           ,ARHO,APOT,ASINGLE,ATRIPLE)
!     **************************************************************************
!     **  GAUSSIAN REPRESENTATION OF ORBITAL PRODUCTS                         **
!     **  AND THEIR ELECTROSTATIC POTENTIALS                                  **
!     **                                                                      **
!     **    CHI_LMN1(R)*CHI_LMN2(R)                                           **
!     **                =SUM_{IJK,IE} |G_{IJK,IE}> ARHO(IJK,IE,LMN1,LMN2)     **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE        
      INTEGER(4),INTENT(IN)  :: NX
      INTEGER(4),INTENT(IN)  :: NIJK
      INTEGER(4),INTENT(IN)  :: NE
      INTEGER(4),INTENT(IN)  :: LNX
      INTEGER(4),INTENT(IN)  :: LOX(LNX)
      INTEGER(4),INTENT(IN)  :: LMNX
      INTEGER(4),INTENT(IN)  :: NPOW
      INTEGER(4),INTENT(IN)  :: NS
      INTEGER(4),INTENT(IN)  :: NP
      INTEGER(4),INTENT(IN)  :: NT
      REAL(8)   ,INTENT(IN)  :: PRODRHO(NPOW,NE,NP)
      REAL(8)   ,INTENT(IN)  :: PRODPOT(NPOW,NE,NP)
      REAL(8)   ,INTENT(IN)  :: SINGLE(NPOW,NE,NS)
      REAL(8)   ,INTENT(IN)  :: TRIPLE(NPOW,NE,NT)
      REAL(8)   ,INTENT(OUT) :: ARHO(NIJK,NE,LMNX,LMNX)
      REAL(8)   ,INTENT(OUT) :: APOT(NIJK,NE,LMNX,LMNX)
      REAL(8)   ,INTENT(OUT) :: ASINGLE(NIJK,NE,LMNX)
      REAL(8)   ,INTENT(OUT) :: ATRIPLE(NIJK,NE,LMNX,LMNX,LMNX)
      INTEGER(4)             :: IS,IP,IT
      INTEGER(4)             :: LN1,L1,IM1,LM1,LMN01,LMN1
      INTEGER(4)             :: LN2,L2,IM2,LM2,LMN02,LMN2
      INTEGER(4)             :: LN3,L3,IM3,LM3,LMN03,LMN3
      INTEGER(4)             :: LR1,IMR1,LMR1
      INTEGER(4)             :: LR2,IMR2,LMR2
      INTEGER(4)             :: J,I
      INTEGER(4)             :: NPOW2
      INTEGER(4)             :: IE
      REAL(8)                :: CG,CG1,CG2 !GAUNT COEFFICIENT
      REAL(8)                :: C(NIJK)
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_EXPANDPRODS')
      ASINGLE=0.D0
      ARHO=0.D0
      APOT=0.D0
      ATRIPLE=0.D0

      IS=0
      IP=0
      IT=0
      LMN01=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        NPOW2=INT(0.5D0*REAL(NX-L1)+1.000001D0)  !R^(L+2N), N=0,NPOW-1
        IF(NPOW2.LT.1) CYCLE
        IS=IS+1
!
!       ==  SINGLE =============================================================
        DO IM1=1,2*L1+1
          LM1=L1**2+IM1
          LMN1=LMN01+IM1
          DO J=0,(NX-L1)/2-1
            CALL GAUSSIAN$YLMTIMESRN(J,LM1,NIJK,C) ! R^(L+2*J)*YLM
            DO IE=1,NE
              ASINGLE(:,IE,LMN1)=ASINGLE(:,IE,LMN1)+C(:)*SINGLE(J+1,IE,IS)
            ENDDO
          ENDDO
        ENDDO
!       ==  SINGLE DONE ========================================================
!        
        LMN02=LMN01
        DO LN2=LN1,LNX
          L2=LOX(LN2)
          DO LR1=ABS(L1-L2),L1+L2,2  ! TRIANGLE RULE
            NPOW2=INT(0.5D0*REAL(NX-LR1)+1.000001D0)  !R^(L+2N), N=0,NPOW2-1
            IF(NPOW2.LT.1) CYCLE
            IP=IP+1
            IF(IP.GT.NP) THEN
              CALL ERROR$MSG('IP OUT OF RANGE')
              CALL ERROR$I4VAL('NP',NP)
              CALL ERROR$STOP('LMTO_EXPANDPRODS')
            END IF
!
!           ==  DOUBLE  ========================================================
            DO IM1=1,2*L1+1
              LM1=L1**2+IM1
              LMN1=LMN01+IM1
              DO IM2=1,2*L2+1
                LM2=L2**2+IM2
                LMN2=LMN02+IM2
                DO IMR1=1,2*LR1+1
                  LMR1=LR1**2+IMR1
                  CALL SPHERICAL$GAUNT(LM1,LM2,LMR1,CG)
                  IF(CG.EQ.0.D0) CYCLE
                  DO J=0,(NX-LR1)/2
                    CALL GAUSSIAN$YLMTIMESRN(J,LMR1,NIJK,C) ! R^(L+2*J)*YLM
                    C(:)=C(:)*CG
                    DO IE=1,NE
                      ARHO(:,IE,LMN1,LMN2)=ARHO(:,IE,LMN1,LMN2) &
     &                                    +C(:)*PRODRHO(J+1,IE,IP)
                      APOT(:,IE,LMN1,LMN2)=APOT(:,IE,LMN1,LMN2) &
     &                                    +C(:)*PRODPOT(J+1,IE,IP)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
!           ==  DOUBLE DONE ====================================================
!
            LMN03=0
            DO LN3=1,LNX
              L3=LOX(LN3)
              DO LR2=ABS(LR1-L3),LR1+L3,2
                NPOW2=INT(0.5D0*REAL(NX-LR2)+1.000001D0) !R^(L+2N), N=0,NPOW-1
                IF(NPOW2.LT.1) CYCLE
                IT=IT+1
                IF(IT.GT.NT) THEN
                  CALL ERROR$MSG('COUNTER FOR TRIPLE TERMS OUT OF RANGE')
                  CALL ERROR$I4VAL('NT',NT)
                  CALL ERROR$STOP('LMTO_EXPANDPRODS')
                END IF
!
!               == TRIPLE TERMS ================================================
                DO IM1=1,2*L1+1
                  LM1=L1**2+IM1
                  LMN1=LMN01+IM1
                  DO IM2=1,2*L2+1
                    LM2=L2**2+IM2
                    LMN2=LMN02+IM2
                    DO IMR1=1,2*LR1+1
                      LMR1=LR1**2+IMR1
                      CALL SPHERICAL$GAUNT(LM1,LM2,LMR1,CG1)
                      IF(CG1.EQ.0.D0) CYCLE
                      DO IM3=1,2*L3+1
                        LM3=L3**2+IM3
                        LMN3=LMN03+IM3
                        DO IMR2=1,2*LR2+1
                          LMR2=LR2**2+IMR2
                          CALL SPHERICAL$GAUNT(LMR1,LM3,LMR2,CG2)
                          IF(CG2.EQ.0.D0) CYCLE
                          DO J=0,(NX-LR2)/2
                            CALL GAUSSIAN$YLMTIMESRN(J,LMR2,NIJK,C) 
                            C(:)=C(:)*CG1*CG2
                            DO IE=1,NE
                              ATRIPLE(:,IE,LMN1,LMN2,LMN3) &
     &                                         =ATRIPLE(:,IE,LMN1,LMN2,LMN3) &
     &                                         +C(:)*TRIPLE(J+1,IE,IT)
                            ENDDO
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
!             == TRIPLE TERMS DONE =============================================
!
              LMN03=LMN03+2*L3+1
            ENDDO
          ENDDO
!         == SYMMETRIZE ========================================================
          IF(LN2.NE.LN1) THEN
            DO IM1=1,2*L1+1
              LMN1=LMN01+IM1
              DO IM2=1,2*L2+1
                LMN2=LMN02+IM2
                ARHO(:,:,LMN2,LMN1)=ARHO(:,:,LMN1,LMN2)
                APOT(:,:,LMN2,LMN1)=APOT(:,:,LMN1,LMN2)
                ATRIPLE(:,:,LMN2,LMN1,:)=ATRIPLE(:,:,LMN1,LMN2,:)
              ENDDO
            ENDDO
          END IF
          LMN02=LMN02+2*L2+1
        ENDDO
        LMN01=LMN01+2*L1+1
      ENDDO
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OVERLAPEVAL()
!     **************************************************************************
!     **  CALCULATEDS OFF-SITE OVERLAP MATRIX "OVERLAP"                       **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR,DENMAT,NSP,OFFSITEX,OVERLAP &
     &                       ,TOFFSITE
      IMPLICIT NONE
      REAL(8)                :: RBAS(3,3)
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4)             :: NND
      INTEGER(4)             :: NN
      INTEGER(4)             :: ISPA,ISPB
      INTEGER(4)             :: LMNXA,LMNXB
      INTEGER(4)             :: N1,N2
      INTEGER(4)             :: LNXA,LNXB
      INTEGER(4)             :: IAT,IATA,IATB
      REAL(8)                :: RA(3),RB(3)
      INTEGER(4)             :: NA,NB
      INTEGER(4)             :: NDIMD
      INTEGER(4)             :: I,L,LN,LM,LMN,LMX
      REAL(8)                :: SVAR
      REAL(8)                :: DR(3)
      REAL(8)                :: DIS
      REAL(8)                :: ROT(3,3)
      REAL(8)   ,ALLOCATABLE :: YLMROT(:,:)
      REAL(8)   ,ALLOCATABLE :: UROTA(:,:)
      REAL(8)   ,ALLOCATABLE :: UROTB(:,:)
      REAL(8)   ,ALLOCATABLE :: OV(:,:)
      REAL(8)   ,ALLOCATABLE :: DOV(:,:)
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_OVERLAPEVAL')
      IF(.NOT.TOFFSITE) THEN
        CALL ERROR$MSG('PARAMETER TOFFSITE MUST BE TRUE')
        CALL ERROR$L4VAL('TOFFSITE',TOFFSITE)
        CALL ERROR$STOP('LMTO_OVERLAPEVAL')
      END IF
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     == ALLOCATE OVERLAP                                                     ==
!     ==========================================================================
      NND=SIZE(DENMAT)
      IF(.NOT.ALLOCATED(OVERLAP)) THEN
        ALLOCATE(OVERLAP(NND))
PRINT*,'ALLOCATE OVERLAP ',NND
        DO NN=1,NND
          OVERLAP(NN)%IAT1=DENMAT(NN)%IAT1
          OVERLAP(NN)%IAT2=DENMAT(NN)%IAT2
          OVERLAP(NN)%IT  =DENMAT(NN)%IT
          OVERLAP(NN)%N1  =DENMAT(NN)%N1
          OVERLAP(NN)%N2  =DENMAT(NN)%N2
          N1  =DENMAT(NN)%N1
          N2  =DENMAT(NN)%N2
          ALLOCATE(OVERLAP(NN)%MAT(N1,N2))
          OVERLAP(NN)%MAT=0.D0
        ENDDO
      ELSE
PRINT*,'NOT ALLOCATE OVERLAP ',NND,SIZE(OVERLAP)
      END IF
      IF(SIZE(OVERLAP).NE.NND) THEN
        CALL ERROR$MSG('SIZE OF DENSITY MATRIX AND OVERLAP MATRIX INONSISTENT')
        CALL ERROR$I4VAL('NND',NND)
        CALL ERROR$I4VAL('SIZE(OVERLAP)',SIZE(OVERLAP))
        CALL ERROR$STOP('LMTO_OVERLAPEVAL')
      END IF
!
!     ==========================================================================
!     == LOOP OVER PAIRS                                                      ==
!     ==========================================================================
      DO NN=1,NND
        IATA=OVERLAP(NN)%IAT1
        IATB=OVERLAP(NN)%IAT2
!       == CONSIDER ONLY OFFSITE TERMS =========================================
        ISPA=ISPECIES(IATA)
        ISPB=ISPECIES(IATB)
        RA(:)=R0(:,IATA)
        RB(:)=R0(:,IATB)+MATMUL(RBAS,REAL(OVERLAP(NN)%IT(:),KIND=8))
        LMNXA=POTPAR(ISPA)%TAILED%LMNX
        LMNXB=POTPAR(ISPB)%TAILED%LMNX
        NA=OVERLAP(NN)%N1
        NB=OVERLAP(NN)%N2
!
!       == ONSITE ELEMENT
        IF(IATA.EQ.IATB.AND.SUM(OVERLAP(NN)%IT**2).EQ.0)  THEN
          ALLOCATE(OV(LMNXA,LMNXB))
          OV=POTPAR(ISPA)%TAILED%OVERLAP
          CALL LMTO_SHRINKDOWNHTNL(IATA,IATB,1,LMNXA,LMNXB,OV &
     &                          ,NA,NB,OVERLAP(NN)%MAT)
          DEALLOCATE(OV)
          CYCLE
        END IF
!
!       ========================================================================
!       == CALCULATE OVERLAP MATRIX WITH Z ATONG THE BOND AXIS                ==
!       ========================================================================
        ALLOCATE(OV(LMNXA,LMNXB))
        ALLOCATE(DOV(LMNXA,LMNXB))
        DIS=SQRT(SUM((RB-RA)**2))
        CALL LMTO_OFFSITEOVERLAP(ISPA,ISPB,DIS,LMNXA,LMNXB,OV,DOV)
!
!       ========================================================================
!       == ROTATE DENSITY MATRIX SO THAT DISTANCE VECTOR POINTS IN Z-DIRECTION=
!       ========================================================================
!       == CONSTRUCT ROTATION MATRIX ===========================================
!       == DISTANCE VECTOR WILL BE NEW Z-DIRECTION =============================
        DR(:)=RB(:)-RA(:)
        DIS=SQRT(SUM(DR**2))
        ROT(:,3)=DR(:)/DIS
!       == FIRST VECTOR IS VECTOR PRODUCT OF THE THE MOST ORTHOGNAL UNIT VECTOR
        I=MINLOC(ABS(DR),1)
        ROT(:,2)=0.D0
        ROT(I,2)=1.D0
        ROT(1,1)=ROT(2,2)*ROT(3,3)-ROT(3,2)*ROT(2,3)
        ROT(2,1)=ROT(3,2)*ROT(1,3)-ROT(1,2)*ROT(3,3)
        ROT(3,1)=ROT(1,2)*ROT(2,3)-ROT(2,2)*ROT(1,3)
        ROT(:,1)=ROT(:,1)/SQRT(SUM(ROT(:,1)**2))
        ROT(1,2)=ROT(2,3)*ROT(3,1)-ROT(3,3)*ROT(2,1)
        ROT(2,2)=ROT(3,3)*ROT(1,1)-ROT(1,3)*ROT(3,1)
        ROT(3,2)=ROT(1,3)*ROT(2,1)-ROT(2,3)*ROT(1,1)
!       == REMOVE INVERSION, IF PRESENT
        SVAR=ROT(1,1)*(ROT(2,2)*ROT(3,3)-ROT(3,2)*ROT(2,3)) &
       &    +ROT(1,2)*(ROT(2,3)*ROT(3,1)-ROT(3,3)*ROT(2,1)) &
       &    +ROT(1,3)*(ROT(2,1)*ROT(3,2)-ROT(3,1)*ROT(2,2))
        IF(SVAR.LE.0.D0) ROT(:,1)=-ROT(:,1)
!
!       == CONSTRUCT TRANSFORMATION MATRIX FOR REAL SPHERICAL HARMONICS ========
        LMX=MAX(MAXVAL(POTPAR(ISPA)%TAILED%LOX),MAXVAL(POTPAR(ISPB)%TAILED%LOX))
        LMX=(LMX+1)**2
        ALLOCATE(YLMROT(LMX,LMX))
        CALL SPHERICAL$ROTATEYLM(LMX,ROT,YLMROT)
!
!       == CONSTRUCT ROTATION MATRIX OF TAILED ORBITALS ========================
        ALLOCATE(UROTA(LMNXA,LMNXA))
        ALLOCATE(UROTB(LMNXB,LMNXB))
        LNXA=POTPAR(ISPA)%TAILED%LNX
        UROTA(:,:)=0.D0
        LMN=1
        DO LN=1,LNXA
          L=POTPAR(ISPA)%TAILED%LOX(LN)
          LM=L**2+1
          UROTA(LMN:LMN+2*L,LMN:LMN+2*L)=YLMROT(LM:LM+2*L,LM:LM+2*L)
          LMN=LMN+2*L+1
        ENDDO
        LNXB=POTPAR(ISPB)%TAILED%LNX
        UROTB(:,:)=0.D0
        LMN=1
        DO LN=1,LNXB
          L=POTPAR(ISPB)%TAILED%LOX(LN)
          LM=L**2+1
          UROTB(LMN:LMN+2*L,LMN:LMN+2*L)=YLMROT(LM:LM+2*L,LM:LM+2*L)
          LMN=LMN+2*L+1
        ENDDO
        DEALLOCATE(YLMROT)
!
!       ========================================================================
!       == ROTATE HAMILTONIAN BACK                                            ==
!       ========================================================================
        OV=MATMUL(UROTA,MATMUL(OV,TRANSPOSE(UROTB)))
        DEALLOCATE(UROTA)
        DEALLOCATE(UROTB)
!
!       ========================================================================
!       == SHRINK DOWN HAMILTONIAN                                            ==
!       ========================================================================
        CALL LMTO_SHRINKDOWNHTNL(IATA,IATB,1,LMNXA,LMNXB,OV &
     &                          ,NA,NB,OVERLAP(NN)%MAT)
        DEALLOCATE(OV)
        DEALLOCATE(DOV)
      ENDDO
!
!     ==========================================================================
!     == CLEAN UP TO AVOID MEMORY LEAK                                        ==
!     ==========================================================================
!!$CALL LMTO$REPORTDENMAT(6)
!!$CALL LMTO$REPORTOVERLAP(6)
!!$STOP 'FORCED'
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEXEVAL(EX)
!     **************************************************************************
!     **  WORK OUT THE ENERGY USING THE LOCAL APPROXIMATION                   **
!     **  TAILED PARTIAL WAVES                                                **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,POTPAR,DENMAT,HAMIL,NSP,OFFSITEX &
     &                       ,HYBRIDSETTING
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: EX
      REAL(8)                :: RBAS(3,3)
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      INTEGER(4)             :: ISP
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4)             :: NND
      INTEGER(4)             :: NN,NNA,NNB
      INTEGER(4)             :: ISPA,ISPB
      INTEGER(4)             :: LMN1A,LMN2A,LMN1B,LMN2B,LMN3A,LMN3B
      INTEGER(4)             :: LMNXA,LMNXB
      INTEGER(4)             :: LNXA,LNXB
      INTEGER(4),ALLOCATABLE :: NNAT(:)    !(NAT) POINTER TO ONSITE TERMS
      INTEGER(4)             :: IAT,IATA,IATB
      REAL(8)                :: RA(3),RB(3)
      INTEGER(4)             :: NA,NB
      INTEGER(4)             :: NDIMD
      INTEGER(4)             :: I,L,LN,LM,LMN,LMX
      REAL(8)                :: SVAR
      REAL(8)                :: DR(3)
      REAL(8)                :: DIS
      REAL(8)                :: ROT(3,3)
      REAL(8)   ,ALLOCATABLE :: YLMROT(:,:)
      REAL(8)   ,ALLOCATABLE :: UROTA(:,:)
      REAL(8)   ,ALLOCATABLE :: UROTB(:,:)
      REAL(8)   ,ALLOCATABLE :: U22(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: DU22(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: U3A1B(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: DU3A1B(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: U3B1A(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: DU3B1A(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: BONDU(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: DBONDU(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: D(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DB(:,:,:)
      REAL(8)   ,ALLOCATABLE :: H(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HB(:,:,:)
INTEGER(4) :: J,K,I1,J1,K1
REAL(8) :: SVAR1,SVAR2
INTEGER(4) :: IX,IND1,IND2,IND3
REAL(8)    :: XDELTA,XSVAR
REAL(8) :: EX1
TYPE H1_TYPE
  INTEGER(4)         :: N
  REAL(8)   ,POINTER :: MAT(:,:,:)
END TYPE H1_TYPE
TYPE(H1_TYPE),ALLOCATABLE :: H1(:)
REAL(8),ALLOCATABLE :: H1A(:,:,:)
REAL(8),ALLOCATABLE :: H1B(:,:,:)
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_OFFSITEXEVAL')
PRINT*,'============ OFFSITEXEVAL ============================='
!!$CALL LMTO$REPORTHAMIL(6)
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     == COLLECT POINTERS NNAT TO ONSITE DENSITY MATRIX ELEMENTS              ==
!     ==========================================================================
      ALLOCATE(NNAT(NAT))
      NNAT=0
      NND=SIZE(DENMAT)
      DO NN=1,NND
        IATA=DENMAT(NN)%IAT1
        IATB=DENMAT(NN)%IAT2
        IF(IATA.EQ.IATB.AND.SUM(DENMAT(NN)%IT**2).EQ.0) NNAT(IATA)=NN
      ENDDO
      DO IAT=1,NAT
        IF(NNAT(IAT).LE.0) THEN
          CALL ERROR$MSG('NO ONSITE TERMS FOUND FOR ATOM')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$STOP('LMTO_OFFSITEX')
        END IF
      ENDDO
!!$ALLOCATE(H1(NAT))
!!$EX1=0.D0
!!$DO IAT=1,NAT
!!$  J=DENMAT(NNAT(IAT))%N1
!!$  H1(IAT)%N=J
!!$  ALLOCATE(H1(IAT)%MAT(J,J,4))
!!$  H1(IAT)%MAT=0.D0
!!$ENDDO  
!
!     ==========================================================================
!     == LOOP OVER PAIRS                                                      ==
!     ==========================================================================
      IF(SIZE(HAMIL).NE.NND) THEN
        CALL ERROR$MSG('SIZE OF DENSITY MATRIX AND HAMILTONIAN INONSISTENT')
        CALL ERROR$STOP('LMTO_OFFSITEX')
      END IF
      EX=0.D0
      DO NN=1,NND
        IATA=DENMAT(NN)%IAT1
        IATB=DENMAT(NN)%IAT2
!       == CONSIDER ONLY OFFSITE TERMS =========================================
        IF(IATA.EQ.IATB.AND.SUM(DENMAT(NN)%IT**2).EQ.0) CYCLE
!!$CALL LMTO$REPORTDENMAT(6)
!!$STOP
!!$IND1=1
!!$IND2=1
!!$IND3=1
!!$XSVAR=DENMAT(NN)%MAT(IND1,IND2,IND3)
!!$XDELTA=1.D-2
!!$IX=-3
!!$1000 CONTINUE
!!$IX=IX+1
!!$DENMAT(NN)%MAT(IND1,IND2,IND3)=XSVAR+XDELTA*REAL(IX,KIND=8)
!!$HAMIL(NN)%MAT(IND1,IND2,IND3)=0.D0
!!$EX=0.D0
        ISPA=ISPECIES(IATA)
        ISPB=ISPECIES(IATB)
        RA(:)=R0(:,IATA)
        RB(:)=R0(:,IATB)+MATMUL(RBAS,REAL(DENMAT(NN)%IT(:),KIND=8))
        NNA=NNAT(IATA)
        NNB=NNAT(IATB)
!
!       ========================================================================
!       == BLOW UP DENSITY MATRIX                                             ==
!       ========================================================================
        LMNXA=POTPAR(ISPA)%TAILED%LMNX
        LMNXB=POTPAR(ISPB)%TAILED%LMNX
        NDIMD=DENMAT(NN)%N3
        NA=DENMAT(NN)%N1
        NB=DENMAT(NN)%N2
        ALLOCATE(D(LMNXA,LMNXB,NDIMD))
        CALL LMTO_BLOWUPDENMATNL(IATA,IATB,NDIMD &
     &                          ,NA,NB,DENMAT(NN)%MAT,LMNXA,LMNXB,D)
        ALLOCATE(DA(LMNXA,LMNXA,NDIMD))
        CALL LMTO_BLOWUPDENMATNL(IATA,IATA,NDIMD &
     &                          ,NA,NA,DENMAT(NNA)%MAT,LMNXA,LMNXA,DA)
        ALLOCATE(DB(LMNXB,LMNXB,NDIMD))
        CALL LMTO_BLOWUPDENMATNL(IATB,IATB,NDIMD &
     &                          ,NB,NB,DENMAT(NNB)%MAT,LMNXB,LMNXB,DB)
!
!       ========================================================================
!       == ROTATE DENSITY MATRIX SO THAT DISTANCE VECTOR POINTS IN Z-DIRECTION=
!       ========================================================================
!       == CONSTRUCT ROTATION MATRIX ===========================================
!       == DISTANCE VECTOR WILL BE NEW Z-DIRECTION =============================
        DR(:)=RB(:)-RA(:)
        DIS=SQRT(SUM(DR**2))
        ROT(:,3)=DR(:)/DIS
!       == FIRST VECTOR IS VECTOR PRODUCT OF THE THE MOST ORTHOGNAL UNIT VECTOR
        I=MINLOC(ABS(DR),1)
        ROT(:,2)=0.D0
        ROT(I,2)=1.D0
        ROT(1,1)=ROT(2,2)*ROT(3,3)-ROT(3,2)*ROT(2,3)
        ROT(2,1)=ROT(3,2)*ROT(1,3)-ROT(1,2)*ROT(3,3)
        ROT(3,1)=ROT(1,2)*ROT(2,3)-ROT(2,2)*ROT(1,3)
        ROT(:,1)=ROT(:,1)/SQRT(SUM(ROT(:,1)**2))
        ROT(1,2)=ROT(2,3)*ROT(3,1)-ROT(3,3)*ROT(2,1)
        ROT(2,2)=ROT(3,3)*ROT(1,1)-ROT(1,3)*ROT(3,1)
        ROT(3,2)=ROT(1,3)*ROT(2,1)-ROT(2,3)*ROT(1,1)
!       == REMOVE INVERSION, IF PRESENT
        SVAR=ROT(1,1)*(ROT(2,2)*ROT(3,3)-ROT(3,2)*ROT(2,3)) &
       &    +ROT(1,2)*(ROT(2,3)*ROT(3,1)-ROT(3,3)*ROT(2,1)) &
       &    +ROT(1,3)*(ROT(2,1)*ROT(3,2)-ROT(3,1)*ROT(2,2))
        IF(SVAR.LE.0.D0) ROT(:,1)=-ROT(:,1)
!
!       == CONSTRUCT TRANSFORMATION MATRIX FOR REAL SPHERICAL HARMONICS ========
        LMX=MAX(MAXVAL(POTPAR(ISPA)%TAILED%LOX),MAXVAL(POTPAR(ISPB)%TAILED%LOX))
        LMX=(LMX+1)**2
        ALLOCATE(YLMROT(LMX,LMX))
        CALL SPHERICAL$ROTATEYLM(LMX,ROT,YLMROT)
!
!       == CONSTRUCT ROTATION MATRIX OF TAILED ORBITALS ========================
        ALLOCATE(UROTA(LMNXA,LMNXA))
        ALLOCATE(UROTB(LMNXB,LMNXB))
        LNXA=POTPAR(ISPA)%TAILED%LNX
        UROTA(:,:)=0.D0
        LMN=1
        DO LN=1,LNXA
          L=POTPAR(ISPA)%TAILED%LOX(LN)
          LM=L**2+1
          UROTA(LMN:LMN+2*L,LMN:LMN+2*L)=YLMROT(LM:LM+2*L,LM:LM+2*L)
          LMN=LMN+2*L+1
        ENDDO
        LNXB=POTPAR(ISPB)%TAILED%LNX
        UROTB(:,:)=0.D0
        LMN=1
        DO LN=1,LNXB
          L=POTPAR(ISPB)%TAILED%LOX(LN)
          LM=L**2+1
          UROTB(LMN:LMN+2*L,LMN:LMN+2*L)=YLMROT(LM:LM+2*L,LM:LM+2*L)
          LMN=LMN+2*L+1
        ENDDO
        DEALLOCATE(YLMROT)
!
!       == ROTATE DENSITY MATRIX ===============================================
        !SUGGESTION: SPEED UP BY EXPLOITING THAT UROT IS SPARSE
        DO I=1,NDIMD
          DA(:,:,I)=MATMUL(TRANSPOSE(UROTA),MATMUL(DA(:,:,I),UROTA))
          D(:,:,I) =MATMUL(TRANSPOSE(UROTA),MATMUL(D(:,:,I) ,UROTB))
          DB(:,:,I)=MATMUL(TRANSPOSE(UROTB),MATMUL(DB(:,:,I),UROTB))
        ENDDO
!
!       ========================================================================
!       == CALCULATE 22-U-TENSOR                                              ==
!       ========================================================================
        ALLOCATE(U22(LMNXA,LMNXA,LMNXB,LMNXB))
        ALLOCATE(DU22(LMNXA,LMNXA,LMNXB,LMNXB))
        ALLOCATE(U3A1B(LMNXA,LMNXA,LMNXA,LMNXB))
        ALLOCATE(DU3A1B(LMNXA,LMNXA,LMNXA,LMNXB))
        ALLOCATE(U3B1A(LMNXB,LMNXB,LMNXB,LMNXA))
        ALLOCATE(DU3B1A(LMNXB,LMNXB,LMNXB,LMNXA))
        ALLOCATE(BONDU(LMNXA,LMNXB,LMNXA,LMNXB))
        ALLOCATE(DBONDU(LMNXA,LMNXB,LMNXA,LMNXB))
        DIS=SQRT(SUM((RB-RA)**2))
        CALL LMTO_OFFSITEX22U(ISPA,ISPB, DIS,LMNXA,LMNXB,U22,DU22)
        CALL LMTO_OFFSITEX31U(ISPA,ISPB, DIS,LMNXA,LMNXB,U3A1B,DU3A1B)
        CALL LMTO_OFFSITEX31U(ISPB,ISPA,-DIS,LMNXB,LMNXA,U3B1A,DU3B1A)
!       == BONDU(1,2,3,4)=INT DX INT DX': A1(X)*B2(X)][A3(X')*B4(X')]/|R-R'|  ==
        CALL LMTO_OFFSITEXBONDU(ISPA,ISPB,DIS,LMNXA,LMNXB,BONDU,DBONDU)
!!$PRINT*,'IATA,IATB ',IATA,IATB,DIS
!!$PRINT*,'DA    ',DA
!!$PRINT*,'DB    ',DB
!!$PRINT*,'D     ',D
!!$PRINT*,'U22   ',U22
!!$PRINT*,'U3A1B ',U3A1B
!!$PRINT*,'U3A1B ',U3B1A
!!$PRINT*,'BONDU ',BONDU

!!$SVAR=0.D0
!!$DO LMN1A=1,LMNXA
!!$  DO LMN1B=1,LMNXB
!!$    DO LMN2A=1,LMNXA
!!$      DO LMN2B=1,LMNXB
!!$        SVAR1=U22(LMN1A,LMN2A,LMN1B,LMN2B)
!!$        SVAR=MAX(SVAR,ABS(SVAR1-U22(LMN1A,LMN2A,LMN2B,LMN1B)))
!!$        SVAR=MAX(SVAR,ABS(SVAR1-U22(LMN2A,LMN1A,LMN1B,LMN2B)))
!!$        SVAR=MAX(SVAR,ABS(SVAR1-U22(LMN2A,LMN1A,LMN2B,LMN1B)))
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
!!$IF(SVAR.GT.1.D-8) THEN
!!$  CALL ERROR$MSG('SYMMETRY DEVIATION U22 ')
!!$  CALL ERROR$R8VAL('DEVIATION ',SVAR)
!!$  CALL ERROR$STOP('LMTO_OFFSITEXEVAL')
!!$END IF
!!$SVAR=0.D0
!!$DO LMN1A=1,LMNXA
!!$  DO LMN2A=1,LMNXA
!!$    DO LMN3A=1,LMNXA
!!$      DO LMN1B=1,LMNXB
!!$        SVAR1=U3A1B(LMN1A,LMN2A,LMN3A,LMN1B)
!!$        SVAR=MAX(SVAR,ABS(SVAR1-U3A1B(LMN2A,LMN1A,LMN3A,LMN1B)))
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
!!$IF(SVAR.GT.1.D-8) THEN
!!$  CALL ERROR$MSG('SYMMETRY DEVIATION U31 ')
!!$  CALL ERROR$R8VAL('DEVIATION ',SVAR)
!!$  CALL ERROR$STOP('LMTO_OFFSITEXEVAL')
!!$END IF
!!$SVAR=0.D0
!!$DO LMN1A=1,LMNXA
!!$  DO LMN1B=1,LMNXB
!!$    DO LMN2A=1,LMNXA
!!$      DO LMN2B=1,LMNXB
!!$        SVAR1=BONDU(LMN1A,LMN1B,LMN2A,LMN2B)
!!$        SVAR=MAX(SVAR,ABS(SVAR1-BONDU(LMN2A,LMN2B,LMN1A,LMN1B)))
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
!!$IF(SVAR.GT.1.D-8) THEN
!!$  CALL ERROR$MSG('SYMMETRY DEVIATION BONDU ')
!!$  CALL ERROR$R8VAL('DEVIATION ',SVAR)
!!$  CALL ERROR$STOP('LMTO_OFFSITEXEVAL')
!!$END IF
!
!       ========================================================================
!       == SWITCH SELECTED TERMS OFF                                          ==
!       ========================================================================
        IF(.NOT.(HYBRIDSETTING(ISPA)%TNDDO.OR.HYBRIDSETTING(ISPB)%TNDDO)) THEN
          U22=0.D0
          DU22=0.D0
        END IF
        IF(.NOT.(HYBRIDSETTING(ISPA)%T31.OR.HYBRIDSETTING(ISPB)%T31)) THEN
!         == FOR EACH BOND ALL 31 TERMS OR NONE ARE CONSIDERED TO ENSURE =======
!         == THAT THE HAMILTONIAN IS HERMITEAN =================================
          U3A1B=0.D0
          DU3A1B=0.D0
          U3B1A=0.D0
          DU3B1A=0.D0
        END IF
        IF(.NOT.(HYBRIDSETTING(ISPA)%TBONDX.OR.HYBRIDSETTING(ISPB)%TBONDX)) THEN
!         == FOR EACH BOND ALL 31 TERMS OR NONE ARE CONSIDERED TO ENSURE =======
!         == THAT THE HAMILTONIAN IS HERMITEAN =================================
          BONDU=0.D0
          DBONDU=0.D0
        END IF
!
!       ========================================================================
!       == ADD UP EXCHANGE ENERGY                                             ==
!       ========================================================================
        ALLOCATE(H(LMNXA,LMNXB,NDIMD))
        ALLOCATE(HA(LMNXA,LMNXA,NDIMD))
        ALLOCATE(HB(LMNXB,LMNXB,NDIMD))
        H=0.D0
        HA=0.D0
        HB=0.D0
!!$ALLOCATE(H1A(LMNXA,LMNXA,NDIMD))
!!$ALLOCATE(H1B(LMNXB,LMNXB,NDIMD))
!!$H1A=0.D0
!!$H1B=0.D0
!PRINT*,'EX  1',EX,IATA,IATB,DENMAT(NN)%IT
!
!       ========================================================================
!       == NDDO TERM:                                                         ==
!       == U22(1,2,3,4)=INT DX1 INT DX2: A1(X1)*A2(X1)][B3(X2)*B4(X2)]/|R1-R2|==
!       == LMN1A,LMN2A TIED TO COORDINATE X, LMN1B,LMN2B TO X'                ==
!       ========================================================================
        DO LMN1B=1,LMNXB
          DO LMN2B=1,LMNXB
            DO LMN1A=1,LMNXA
              DO LMN2A=1,LMNXA
                SVAR=-0.25D0*U22(LMN1A,LMN2A,LMN2B,LMN1B)
                EX=EX+SVAR*SUM(D(LMN1A,LMN1B,:)*D(LMN2A,LMN2B,:))
                H(LMN1A,LMN1B,:)=H(LMN1A,LMN1B,:)+SVAR*D(LMN2A,LMN2B,:)
                H(LMN2A,LMN2B,:)=H(LMN2A,LMN2B,:)+SVAR*D(LMN1A,LMN1B,:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!PRINT*,'EX  2',EX
!
!       ========================================================================
!       == 3-1 TERMS:                                                         ==
!       == U3A1B(1,2,3,4)=INT DX1 INT DX2:                                    ==
!       ==                          * [A1(X1)*A2(X1)][A3(X2)*B4(X2)]/|X1-X2|  ==
!       == LMN1A,LMN2A TIED TO COORDINATE X, LMN3A,LMN1B TO X'                ==
!       ========================================================================
        DO LMN1A=1,LMNXA
          DO LMN2A=1,LMNXA
            DO LMN3A=1,LMNXA 
              DO LMN1B=1,LMNXB
                SVAR=-0.25D0*U3A1B(LMN1A,LMN2A,LMN3A,LMN1B)
                EX=EX+SVAR*SUM(D(LMN2A,LMN1B,:)*DA(LMN1A,LMN3A,:))
                HA(LMN1A,LMN3A,:)=HA(LMN1A,LMN3A,:)+SVAR*D(LMN2A,LMN1B,:)
                H(LMN2A,LMN1B,:) =H(LMN2A,LMN1B,:) +SVAR*DA(LMN1A,LMN3A,:)
                EX=EX+SVAR*SUM(D(LMN2A,LMN1B,:)*DA(LMN3A,LMN1A,:))
                HA(LMN3A,LMN1A,:)=HA(LMN3A,LMN1A,:)+SVAR*D(LMN2A,LMN1B,:)
                H(LMN2A,LMN1B,:) =H(LMN2A,LMN1B,:) +SVAR*DA(LMN3A,LMN1A,:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO               
!PRINT*,'EX  3',EX
        DO LMN1B=1,LMNXB
          DO LMN2B=1,LMNXB
            DO LMN3B=1,LMNXB 
              DO LMN1A=1,LMNXA
                SVAR=-0.25D0*U3B1A(LMN1B,LMN2B,LMN3B,LMN1A)
                EX=EX+SVAR*SUM(D(LMN1A,LMN2B,:)*DB(LMN1B,LMN3B,:))
                HB(LMN1B,LMN3B,:)=HB(LMN1B,LMN3B,:)+SVAR*D(LMN1A,LMN2B,:)
                H(LMN1A,LMN2B,:) =H(LMN1A,LMN2B,:) +SVAR*DB(LMN1B,LMN3B,:)
                EX=EX+SVAR*SUM(D(LMN1A,LMN2B,:)*DB(LMN3B,LMN1B,:))
                HB(LMN3B,LMN1B,:)=HB(LMN3B,LMN1B,:)+SVAR*D(LMN1A,LMN2B,:)
                H(LMN1A,LMN2B,:) =H(LMN1A,LMN2B,:) +SVAR*DB(LMN3B,LMN1B,:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!PRINT*,'EX  4',EX
!
!       ========================================================================
!       == 2ND ORDER DIFFERENTIAL OVERLAP TERMS:                              ==
!       == BONDU(1,2,3,4)=INT DX1 INT DX2:                                    ==
!       ==                          * [A1(X1)*B2(X1)][A3(X2)*B4(X2)]/|X1-X2|  ==
!       == LMN1A,LMN1B TIED TO COORDINATE X1, LMN2A,LMN2B TO X2               ==
!       ========================================================================
        DO LMN1A=1,LMNXA
          DO LMN1B=1,LMNXB
            DO LMN2A=1,LMNXA 
              DO LMN2B=1,LMNXB
                SVAR=-0.25D0*BONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                EX=EX+SVAR*SUM(DA(LMN1A,LMN2A,:)*DB(LMN1B,LMN2B,:))
                HA(LMN1A,LMN2A,:)=HA(LMN1A,LMN2A,:)+SVAR*DB(LMN1B,LMN2B,:)
                HB(LMN1B,LMN2B,:)=HB(LMN1B,LMN2B,:)+SVAR*DA(LMN1A,LMN2A,:)
!!$SVAR=-0.25D0*BONDU(LMN1A,LMN1B,LMN2A,LMN2B)
!!$EX1=EX1+SVAR*SUM(DA(LMN1A,LMN2A,:)*DB(LMN1B,LMN2B,:))
!!$H1A(LMN1A,LMN2A,:)=H1A(LMN1A,LMN2A,:)+SVAR*DB(LMN1B,LMN2B,:)
!!$H1B(LMN1B,LMN2B,:)=H1B(LMN1B,LMN2B,:)+SVAR*DA(LMN1A,LMN2A,:)
                SVAR=-0.25D0*BONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                EX=EX+SVAR*SUM(D(LMN1A,LMN2B,:)*D(LMN2A,LMN1B,:))
                H(LMN1A,LMN2B,:)=H(LMN1A,LMN2B,:)+SVAR*D(LMN2A,LMN1B,:)
                H(LMN2A,LMN1B,:)=H(LMN2A,LMN1B,:)+SVAR*D(LMN1A,LMN2B,:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!PRINT*,'EX  5',EX
        DEALLOCATE(U22)
        DEALLOCATE(DU22)
        DEALLOCATE(U3A1B)
        DEALLOCATE(DU3A1B)
        DEALLOCATE(U3B1A)
        DEALLOCATE(DU3B1A)
        DEALLOCATE(BONDU)
        DEALLOCATE(DBONDU)
!
!       ========================================================================
!       == ROTATE HAMILTONIAN BACK                                            ==
!       ========================================================================
        DO I=1,NDIMD
          HA(:,:,I)=MATMUL(UROTA,MATMUL(HA(:,:,I),TRANSPOSE(UROTA)))
          H(:,:,I) =MATMUL(UROTA,MATMUL(H(:,:,I) ,TRANSPOSE(UROTB)))
          HB(:,:,I)=MATMUL(UROTB,MATMUL(HB(:,:,I),TRANSPOSE(UROTB)))
!!$H1A(:,:,I)=MATMUL(UROTA,MATMUL(H1A(:,:,I),TRANSPOSE(UROTA)))
!!$H1B(:,:,I)=MATMUL(UROTB,MATMUL(H1B(:,:,I),TRANSPOSE(UROTB)))
        ENDDO
        DEALLOCATE(UROTA)
        DEALLOCATE(UROTB)
!
!       ========================================================================
!       == SHRINK DOWN HAMILTONIAN                                            ==
!       ========================================================================
!       __D,DA IS RE-USED AS WORK ARRAY TO COLLECT HAMILTONIANS______________
        D=0.D0
        DA=0.D0
        DB=0.D0
        CALL LMTO_SHRINKDOWNHTNL(IATA,IATB,NDIMD,LMNXA,LMNXB,H &
     &                                                      ,NA,NB,D(:NA,:NB,:))
        HAMIL(NN)%MAT(:,:,:)=HAMIL(NN)%MAT(:,:,:)+D(:NA,:NB,:)
        CALL LMTO_SHRINKDOWNHTNL(IATA,IATA,NDIMD,LMNXA,LMNXA,HA &
     &                                                     ,NA,NA,DA(:NA,:NA,:))
        HAMIL(NNA)%MAT(:,:,:)=HAMIL(NNA)%MAT(:,:,:)+DA(:NA,:NA,:)
        CALL LMTO_SHRINKDOWNHTNL(IATB,IATB,NDIMD,LMNXB,LMNXB,HB &
     &                                                     ,NB,NB,DB(:NB,:NB,:))
        HAMIL(NNB)%MAT(:,:,:)=HAMIL(NNB)%MAT(:,:,:)+DB(:NB,:NB,:)

!!$DA=0.D0
!!$CALL LMTO_SHRINKDOWNHTNL(IATA,IATA,NDIMD,LMNXA,LMNXA,HA,NA,NA,DA(:NA,:NA,:))
!!$H1(IATA)%MAT(:,:,:)=H1(IATA)%MAT+DA(:NA,:NA,:)
!!$DB=0.D0
!!$CALL LMTO_SHRINKDOWNHTNL(IATB,IATB,NDIMD,LMNXB,LMNXB,HB,NB,NB,DB(:NB,:NB,:))
!!$H1(IATB)%MAT(:,:,:)=H1(IATB)%MAT(:,:,:)+DB(:NB,:NB,:)
!!$DEALLOCATE(H1A)
!!$DEALLOCATE(H1B)

        DEALLOCATE(D)
        DEALLOCATE(H)
        DEALLOCATE(DA)
        DEALLOCATE(HA)
        DEALLOCATE(DB)
        DEALLOCATE(HB)
!!$WRITE(*,*)XDELTA*REAL(IX,8),EX,HAMIL(NN)%MAT(IND1,IND2,IND3)
!!$IF(IX.EQ.3) STOP 'FORCED'
!!$GOTO 1000
      ENDDO
!
!     ==========================================================================
!     == CLEAN UP TO AVOID MEMORY LEAK                                        ==
!     ==========================================================================
!!$PRINT*,'EX1 ',EX1
!!$DO IAT=1,NAT
!!$  LMNXA=H1(IAT)%N
!!$  DO LMN=1,LMNXA
!!$    WRITE(6,FMT='("IAT=",I3," H1=",100F10.5)')IAT,H1(IAT)%MAT(LMN,:,1)
!!$  ENDDO
!!$  DEALLOCATE(H1(IAT)%MAT)
!!$ENDDO
!!$DEALLOCATE(H1)
!!$STOP 'FORCED'

!PRINT*,'EX ',EX
!!$DO IAT=1,NAT
!!$  NN=NNAT(IAT)
!!$  LMNXA=DENMAT(NN)%N1
!!$  DO LMN=1,LMNXA
!!$    WRITE(6,FMT='("IAT=",I3," D=",100F10.5)')IAT,DENMAT(NN)%MAT(LMN,:,1)
!!$  ENDDO
!!$  DO LMN=1,LMNXA
!!$    WRITE(6,FMT='("IAT=",I3," H=",100F10.5)')IAT,HAMIL(NN)%MAT(LMN,:,1)
!!$  ENDDO
!!$ENDDO
      DEALLOCATE(NNAT)

!!$CALL LMTO$REPORTDENMAT(6)
!!$CALL LMTO$REPORTHAMIL(6)
!!$STOP 'FORCED'
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFXINT()
!     **************************************************************************
!     ** COMPUTES OFFSITE MATRIX ELEMENTS 
!     ** WATCH PARALLELIZATION!!!
!     **************************************************************************
      USE LMTO_MODULE,ONLY : POTPAR,OFFSITEX,NSP
      IMPLICIT NONE
      INTEGER(4)        :: GID1,GID2
      INTEGER(4)        :: NR1,NR2
      INTEGER(4)        :: LMNX1,LMNX2
      INTEGER(4)        :: LMN1,LMN2,LMN3,LMN4
      INTEGER(4)        :: NDIS  !#(GRID POINTS)
      INTEGER(4)        :: NF    !#(FIT FUNCTIONS)
      INTEGER(4)        :: ISP1,ISP2,I
      REAL(8)           :: SVAR
      REAL(8)           :: DIS 
      REAL(8)           :: DISMIN,DISMAX !BOUNDS FOR DISTANCE GRID
      REAL(8)  ,PARAMETER :: TOLERANCE=1.D-3
!     **************************************************************************
                                  CALL TRACE$PUSH('LMTO_OFFXINT')
PRINT*,'STARTING INITIALIZATION OF LMTO_OFFSITE'
      ALLOCATE(OFFSITEX(NSP,NSP))
!
!     ==========================================================================
!     == DEFINE ARRAY OF DISTANCES AND DECAY CONSTANTS                        ==
!     == DISTANCE DEPENDENCE WILL BE F(DIS)=SUM_I C_I * EXP(-LAMBDA_I*DIS)    ==
!     ==========================================================================
      DO ISP1=1,NSP
        DO ISP2=1,NSP
          NDIS=5
          OFFSITEX(ISP1,ISP2)%NDIS=NDIS
          ALLOCATE(OFFSITEX(ISP1,ISP2)%DIS(NDIS))
          SVAR=POTPAR(ISP1)%RAD+POTPAR(ISP2)%RAD
          DISMIN=0.5D0*SVAR
!!$PRINT*,'FUDGE WARNING!!!! DISMIN SET TO 1.4 FOR H2 TEST'
!!$DISMIN=1.4D0
          DISMAX=2.D0*SVAR
          DO I=1,NDIS
            DIS=DISMIN+(DISMAX-DISMIN)*REAL(I-1,KIND=8)/REAL(NDIS-1,KIND=8)
            OFFSITEX(ISP1,ISP2)%DIS(I)=DIS
          ENDDO
!
!         == DEFINE DECAY CONSTANTS FOR INTERPOLATING FUNCTIONS ================
          NF=3
          IF(NDIS.LT.NF) THEN
            CALL ERROR$MSG('NUMBER OF INTERPOLATING FUNCTIONS MUST')
            CALL ERROR$MSG('BE EQUAL OR GREATER THAN NUMBER OF GRID POINTS')
            CALL ERROR$STOP('LMTO_OFFSITEXINT')
          END IF
          IF(NF.NE.3) THEN
            CALL ERROR$MSG('NUMBER OF INTERPOLATING FUNCTIONS HARDWIRED TO 3')
            CALL ERROR$STOP('LMTO_OFFSITEXINT')
          END IF
          OFFSITEX(ISP1,ISP2)%NF=NF
          ALLOCATE(OFFSITEX(ISP1,ISP2)%LAMBDA(NF))
          OFFSITEX(ISP1,ISP2)%LAMBDA(1)=1.D0/0.5D0
          OFFSITEX(ISP1,ISP2)%LAMBDA(2)=1.D0/1.D0
          OFFSITEX(ISP1,ISP2)%LAMBDA(3)=1.D0/2.D0
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE OVERLAP MATRIX ELEMENTS                                    ==
!     ==========================================================================
      DO ISP1=1,NSP
        GID1=POTPAR(ISP1)%TAILED%GID
        CALL RADIAL$GETI4(GID1,'NR',NR1)
        DO ISP2=1,NSP
          GID2=POTPAR(ISP2)%TAILED%GID
          CALL RADIAL$GETI4(GID2,'NR',NR2)
PRINT*,'DOING OVERLAP....',ISP1,ISP2
          CALL LMTO_OFFSITEOVERLAPSETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE EXCHANGE INTEGRALS WITH TWO ORBITALS ON EITHER SIDE        ==
!     ==========================================================================
      DO ISP1=1,NSP
        GID1=POTPAR(ISP1)%TAILED%GID
        CALL RADIAL$GETI4(GID1,'NR',NR1)
        DO ISP2=1,NSP
          GID2=POTPAR(ISP2)%TAILED%GID
          CALL RADIAL$GETI4(GID2,'NR',NR2)
PRINT*,'DOING X22 ....',ISP1,ISP2
          CALL LMTO_OFFSITEX22SETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE EXCHANGE INTEGRALS WITH THREE ORBITALS ON THE FIRST SITE   ==
!     ==                               AND ONE ORBITAL ON THE SECOND SITE     ==
!     ==========================================================================
      DO ISP1=1,NSP
        GID1=POTPAR(ISP1)%TAILED%GID
        CALL RADIAL$GETI4(GID1,'NR',NR1)
        DO ISP2=1,NSP
          GID2=POTPAR(ISP2)%TAILED%GID
          CALL RADIAL$GETI4(GID2,'NR',NR2)
PRINT*,'DOING X31 ....',ISP1,ISP2
          CALL LMTO_OFFSITEX31SETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE EXCHANGE INTEGRALS DABAB  ==
!     ==========================================================================
PRINT*,'DOING XABAB....'
      CALL LMTO_TAILEDGAUSSOFFSITEU()   ! ROUTINE PARALLELIZES OVER 'MONOMER'
!
!     ==========================================================================
!     == CONVERT INTEGRALS INTO COEFFICIENTS OF INTERPOLATING FUNCTION        ==
!     ==========================================================================
PRINT*,'CONVERTING....'
      CALL LMTO_OFFSITEXCONVERT()
PRINT*,'INITIALIZATION OF OFFSITE DONE...'
                                  CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEOVERLAPSETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
!     **************************************************************************
!     ** CALCULATES THE OVERLAP OF TAILED ORBITALS ON A RADIAL INTERPOLATION  **
!     ** GRID                                                                 **
!     ** ROUTINE IS PARALLELIZED OVER 'MONOMER'                               **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR,OFFSITEX
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      INTEGER(4),INTENT(IN) :: NR1
      INTEGER(4),INTENT(IN) :: NR2
      REAL(8)   ,INTENT(IN) :: TOLERANCE
      REAL(8)   ,PARAMETER  :: TOLMIN=1.D-8
      INTEGER(4)            :: GID1,GID2
      INTEGER(4)            :: LRX1,LRX2
      INTEGER(4)            :: LNX1,LNX2
      INTEGER(4)            :: IND
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: L1,L2
      INTEGER(4)            :: MABS
      REAL(8)               :: INTEGRAL,DINTEGRAL
      INTEGER(4)            :: LMRX
      REAL(8)               :: PHI1(NR1),PHI2(NR2)
      REAL(8)               :: RGRID1(NR1),RGRID2(NR2)
      INTEGER(4)            :: IDIS
      INTEGER(4)            :: NDIS
      REAL(8)               :: DIS
      REAL(8) ,ALLOCATABLE  :: TOLFAC1(:)
      REAL(8) ,ALLOCATABLE  :: TOLFAC2(:)
      REAL(8)               :: TOL
      INTEGER(4)            :: NTASKS,THISTASK,COUNT
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     ==========================================================================
!     == PREPARATION                                                          ==
!     ==========================================================================
      GID1=POTPAR(ISP1)%TAILED%GID
      GID2=POTPAR(ISP2)%TAILED%GID
      LNX1=POTPAR(ISP1)%TAILED%LNX
      LNX2=POTPAR(ISP2)%TAILED%LNX
      CALL SETUP$ISELECT(ISP1)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX1=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      CALL SETUP$ISELECT(ISP2)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX2=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      NDIS=OFFSITEX(ISP1,ISP2)%NDIS
!
!     ==========================================================================
!     == OBTAIN NORM TO SET TOLERANCE                                         ==
!     ==========================================================================
      ALLOCATE(TOLFAC1(LNX1))
      ALLOCATE(TOLFAC2(LNX2))
      CALL RADIAL$R(GID1,NR1,RGRID1)
      DO LN1=1,LNX1
        CALL RADIAL$INTEGRAL(GID1,NR1 &
     &              ,(RGRID1(:)*POTPAR(ISP1)%TAILED%AEF(:,LN1))**2,TOLFAC1(LN1))
      ENDDO
      CALL RADIAL$R(GID2,NR2,RGRID2)
      DO LN2=1,LNX2
        CALL RADIAL$INTEGRAL(GID2,NR2 &
    &               ,(RGRID2(:)*POTPAR(ISP2)%TAILED%AEF(:,LN2))**2,TOLFAC2(LN2))
      ENDDO
      TOLFAC1=SQRT(TOLFAC1)
      TOLFAC2=SQRT(TOLFAC2)
!
!     ==========================================================================
!     == CALCULATE NUMBER OF INTEGRALS                                        ==
!     ==========================================================================
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        DO LN2=1,LNX2
          L2=POTPAR(ISP2)%TAILED%LOX(LN2)
          DO MABS=0,MIN(L1,L2)
            IND=IND+1
          ENDDO
        ENDDO
      ENDDO
      ALLOCATE(OFFSITEX(ISP1,ISP2)%OVERLAP(NDIS,IND))
      OFFSITEX(ISP1,ISP2)%OVERLAP(:,:)=0.D0
!
!     ==========================================================================
!     == DETERMINE INTEGRALS                                                 ==
!     ==========================================================================
      COUNT=0
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        PHI1(:)=POTPAR(ISP1)%TAILED%AEF(:,LN1) 
        DO LN2=1,LNX2
          L2=POTPAR(ISP2)%TAILED%LOX(LN2)
          PHI2(:)=POTPAR(ISP2)%TAILED%AEF(:,LN2)
          TOL=TOLERANCE*TOLFAC1(LN1)*TOLFAC2(LN2)
          TOL=MAX(TOLMIN,TOL)
          DO MABS=0,MIN(L1,L2)
            IND=IND+1
            DO IDIS=1,NDIS
              COUNT=COUNT+1
              IF(MOD(COUNT-1,NTASKS).NE.THISTASK-1) CYCLE
              DIS=OFFSITEX(ISP1,ISP2)%DIS(IDIS)
!PRINT*,'IND ',IND,LN1,LN2,MABS,DIS,TOL
              CALL LMTO_TWOCENTER(L1,MABS,GID1,NR1,PHI1 &
       &                         ,L2,MABS,GID2,NR2,PHI2 &
       &                         ,DIS,TOL,INTEGRAL)
!PRINT*,'   ',DIS,INTEGRAL
              OFFSITEX(ISP1,ISP2)%OVERLAP(IDIS,IND)=INTEGRAL
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL MPE$COMBINE('MONOMER','+',OFFSITEX(ISP1,ISP2)%OVERLAP)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEX22SETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR,OFFSITEX
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      INTEGER(4),INTENT(IN) :: NR1
      INTEGER(4),INTENT(IN) :: NR2
      REAL(8)   ,INTENT(IN) :: TOLERANCE
      REAL(8)   ,PARAMETER  :: RX=0.1D0
      REAL(8)   ,PARAMETER  :: TOLMIN=1.D-8
      INTEGER(4)            :: GID1,GID2
      INTEGER(4)            :: LRX1,LRX2
      INTEGER(4)            :: LNX1,LNX2
      INTEGER(4)            :: IND
      INTEGER(4)            :: LN1,LN2,LN3,LN4
      INTEGER(4)            :: L1,L2,L3,L4
      INTEGER(4)            :: LR1,LR2
      INTEGER(4)            :: MABS
      REAL(8)               :: INTEGRAL,DINTEGRAL
      INTEGER(4)            :: LMRX
      REAL(8)               :: RHO12(NR1),RHO34(NR2),POT12(NR1)
      REAL(8)               :: RGRID1(NR1),RGRID2(NR2)
      INTEGER(4)            :: IDIS
      INTEGER(4)            :: NDIS
      REAL(8)               :: DIS
      REAL(8) ,ALLOCATABLE  :: TOLFAC1(:)
      REAL(8) ,ALLOCATABLE  :: TOLFAC2(:)
      REAL(8)               :: TOL
      INTEGER(4)            :: IR
      REAL(8)               :: AUX(NR2),SVAR,SVAR1,SVAR2,A,B
      REAL(8)               :: PI,SQ4PI,SQ4PITHIRD
      INTEGER(4)            :: THISTASK,NTASKS,COUNT
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      SQ4PI=SQRT(4.D0*PI)
      SQ4PITHIRD=SQRT(4.D0*PI/3.D0)
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     ==========================================================================
!     == PREPARATION                                                          ==
!     ==========================================================================
      GID1=POTPAR(ISP1)%TAILED%GID
      GID2=POTPAR(ISP2)%TAILED%GID
      LNX1=POTPAR(ISP1)%TAILED%LNX
      LNX2=POTPAR(ISP2)%TAILED%LNX
      CALL SETUP$ISELECT(ISP1)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX1=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      CALL SETUP$ISELECT(ISP2)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX2=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      NDIS=OFFSITEX(ISP1,ISP2)%NDIS
!
!     ==========================================================================
!     == OBTAIN NORM TO SET TOLERANCE                                         ==
!     ==========================================================================
      ALLOCATE(TOLFAC1(LNX1))
      ALLOCATE(TOLFAC2(LNX2))
      CALL RADIAL$R(GID1,NR1,RGRID1)
      DO LN1=1,LNX1
        CALL RADIAL$INTEGRAL(GID1,NR1 &
     &              ,(RGRID1(:)*POTPAR(ISP1)%TAILED%AEF(:,LN1))**2,TOLFAC1(LN1))
      ENDDO
      CALL RADIAL$R(GID2,NR2,RGRID2)
      DO LN2=1,LNX2
        CALL RADIAL$INTEGRAL(GID2,NR2 &
    &               ,(RGRID2(:)*POTPAR(ISP2)%TAILED%AEF(:,LN2))**2,TOLFAC2(LN2))
      ENDDO
      TOLFAC1=SQRT(TOLFAC1)
      TOLFAC2=SQRT(TOLFAC2)
!
!     ==========================================================================
!     == CALCULATE NUMBER OF INTEGRALS                                        ==
!     ==========================================================================
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%TAILED%LOX(LN2)
          DO LN3=1,LNX2
            L3=POTPAR(ISP2)%TAILED%LOX(LN3)
            DO LN4=LN3,LNX2
              L4=POTPAR(ISP2)%TAILED%LOX(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                DO LR2=ABS(L3-L4),MIN(L3+L4,LRX2),2
                  DO MABS=0,MIN(LR1,LR2)
                    IND=IND+1
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ALLOCATE(OFFSITEX(ISP1,ISP2)%X22(NDIS,IND))
!
!     ==========================================================================
!     == DETERMINE INTEGRALS                                                 ==
!     ==========================================================================
      OFFSITEX(ISP1,ISP2)%X22=0.D0
      COUNT=0
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%TAILED%LOX(LN2)
          RHO12(:)=POTPAR(ISP1)%TAILED%AEF(:,LN1) &
       &          *POTPAR(ISP1)%TAILED%AEF(:,LN2)
          DO LN3=1,LNX2
            L3=POTPAR(ISP2)%TAILED%LOX(LN3)
            DO LN4=LN3,LNX2
              L4=POTPAR(ISP2)%TAILED%LOX(LN4)
              RHO34(:)=POTPAR(ISP2)%TAILED%AEF(:,LN3) &
       &              *POTPAR(ISP2)%TAILED%AEF(:,LN4)
              TOL=TOLERANCE*TOLFAC1(LN1)*TOLFAC1(LN2) &
       &                   *TOLFAC2(LN3)*TOLFAC2(LN4)
              TOL=MAX(TOLMIN,TOL)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                CALL RADIAL$POISSON(GID1,NR1,LR1,RHO12,POT12)
!!$!               ==SUBTRACT OUT LONG RANGE PART INCLUDING MONOPOLE AND DIPOLE TERMS
!!$                IF(LR1.EQ.0) THEN
!!$                  SVAR=POTPAR(ISP1)%TAILED%QLN(1,LN1,LN2)*SQ4PI
!!$                  POT12(:)=POT12(:)-SVAR/RGRID1(:)
!!$                ELSE IF(LR1.EQ.1) THEN
!!$                  SVAR=POTPAR(ISP1)%TAILED%QLN(2,LN1,LN2)*SQ4PITHIRD
!!$                  POT12(:)=POT12(:)-SVAR/RGRID1(:)**2
!!$                END IF
!
                DO LR2=ABS(L3-L4),MIN(L3+L4,LRX2),2
!!$                  CALL RADIAL$POISSON(GID2,NR2,LR2,RHO34,POT34)
!!$!                 ==SUBTRACT OUT LONG RANGE PART INCLUDING MONOPOLE AND DIPOLE TERMS
!!$                  IF(LR2.EQ.0) THEN
!!$                    SVAR=POTPAR(ISP2)%TAILED%QLN(1,LN3,LN4)*SQ4PI
!!$                    POT34(:)=POT34(:)-SVAR/RGRID2(:)
!!$                  ELSE IF(LR2.EQ.1) THEN
!!$                    SVAR=POTPAR(ISP1)%TAILED%QLN(2,LN3,LN4)*SQ4PITHIRD
!!$                    POT34(:)=POT34(:)-SVAR/RGRID2(:)**2
!!$                  END IF
!
                  DO MABS=0,MIN(LR1,LR2)
                    IND=IND+1
                    DO IDIS=1,NDIS
                      COUNT=COUNT+1
!                     == OFFSITEX WILL BE ADDED TOGETHER IN THE CALLING ROUTINE
                      IF(MOD(COUNT-1,NTASKS).NE.THISTASK-1) CYCLE
                      DIS=OFFSITEX(ISP1,ISP2)%DIS(IDIS)
                      CALL LMTO_TWOCENTER(LR1,MABS,GID1,NR1,POT12 &
                                         ,LR2,MABS,GID2,NR2,RHO34 &
       &                                 ,DIS,TOL,INTEGRAL)
!                     == SUBTRACT OUT LONG RANGE PART TO ALLOW INTERPOLATION ===
!                     == WILL BE ADDED AGAIN SUBTRACT OUT LONG RANGE PART ======
!                     == TO ALLOW INTERPOLATION =======
                      IF(LR1.LE.1.AND.LR2.LE.1) THEN
                        SVAR=POTPAR(ISP1)%TAILED%QLN(LR1+1,LN1,LN2) &
     &                      *POTPAR(ISP2)%TAILED%QLN(LR2+1,LN3,LN4)
                        SVAR=SVAR/DIS**(LR1+LR2+1)
                        INTEGRAL =INTEGRAL -SVAR
                        SVAR=-REAL(LR1+LR2+1)*SVAR/DIS
                        DINTEGRAL=DINTEGRAL-SVAR
                      END IF
                      OFFSITEX(ISP1,ISP2)%X22(IDIS,IND)=INTEGRAL
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL MPE$COMBINE('MONOMER','+',OFFSITEX(ISP1,ISP2)%X22)
!!$PRINT*,'X22  FOR ',ISP1,ISP2
!!$DO IDIS=1,NDIS
!!$PRINT*,OFFSITEX(ISP1,ISP2)%DIS(IDIS),OFFSITEX(ISP1,ISP2)%X22(IDIS,1)/(4.D0*PI)
!!$ENDDO
!!$STOP 'FORCED'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEX31SETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
!     **************************************************************************
!     **  CALC. U-TENSOR FOR THREE ORBITALS ON ONE ATOM AND ON ON THE OTHER   **
!     **  ROUTINE IS PARALLELIZED OVER 'MONOMER'                              **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR,OFFSITEX
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      INTEGER(4),INTENT(IN) :: NR1
      INTEGER(4),INTENT(IN) :: NR2
      REAL(8)   ,INTENT(IN) :: TOLERANCE
      REAL(8)   ,PARAMETER  :: TOLMIN=1.D-8
      INTEGER(4)            :: GID1,GID2
      INTEGER(4)            :: LRX1,LRX2
      INTEGER(4)            :: LNX1,LNX2
      INTEGER(4)            :: IND
      INTEGER(4)            :: LN1,LN2,LN3,LN4
      INTEGER(4)            :: L1,L2,L3,L4
      INTEGER(4)            :: LR1,LR2
      INTEGER(4)            :: MABS
      REAL(8)               :: INTEGRAL,DINTEGRAL
      INTEGER(4)            :: LMRX
      REAL(8)               :: RHO12(NR1),POT12(NR1),A123(NR1),PHI4(NR2)
      REAL(8)               :: RGRID1(NR1),RGRID2(NR2)
      INTEGER(4)            :: IDIS
      INTEGER(4)            :: NDIS
      REAL(8)               :: DIS
      REAL(8) ,ALLOCATABLE  :: TOLFAC1(:)
      REAL(8) ,ALLOCATABLE  :: TOLFAC2(:)
      REAL(8)               :: TOL
      INTEGER(4)            :: NTASKS,THISTASK,COUNT
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     ==========================================================================
!     == PREPARATION                                                          ==
!     ==========================================================================
      GID1=POTPAR(ISP1)%TAILED%GID
      GID2=POTPAR(ISP2)%TAILED%GID
      LNX1=POTPAR(ISP1)%TAILED%LNX
      LNX2=POTPAR(ISP2)%TAILED%LNX
      CALL SETUP$ISELECT(ISP1)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX1=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      CALL SETUP$ISELECT(ISP2)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX2=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      NDIS=OFFSITEX(ISP1,ISP2)%NDIS
!
!     ==========================================================================
!     == OBTAIN NORM TO SET TOLERANCE                                         ==
!     ==========================================================================
      ALLOCATE(TOLFAC1(LNX1))
      ALLOCATE(TOLFAC2(LNX2))
      CALL RADIAL$R(GID1,NR1,RGRID1)
      DO LN1=1,LNX1
        CALL RADIAL$INTEGRAL(GID1,NR1 &
     &              ,(RGRID1(:)*POTPAR(ISP1)%TAILED%AEF(:,LN1))**2,TOLFAC1(LN1))
      ENDDO
      CALL RADIAL$R(GID2,NR2,RGRID2)
      DO LN2=1,LNX2
        CALL RADIAL$INTEGRAL(GID2,NR2 &
    &               ,(RGRID2(:)*POTPAR(ISP2)%TAILED%AEF(:,LN2))**2,TOLFAC2(LN2))
      ENDDO
      TOLFAC1=SQRT(TOLFAC1)
      TOLFAC2=SQRT(TOLFAC2)
!
!     ==========================================================================
!     == CALCULATE NUMBER OF INTEGRALS                                        ==
!     ==========================================================================
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%TAILED%LOX(LN2)
          DO LN3=1,LNX1
            L3=POTPAR(ISP1)%TAILED%LOX(LN3)
            DO LN4=1,LNX2
              L4=POTPAR(ISP2)%TAILED%LOX(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
                    IND=IND+1
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ALLOCATE(OFFSITEX(ISP1,ISP2)%X31(NDIS,IND))
      OFFSITEX(ISP1,ISP2)%X31(:,:)=0.D0
!
!     ==========================================================================
!     == DETERMINE INTEGRALS                                                 ==
!     ==========================================================================
      COUNT=0
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%TAILED%LOX(LN2)
          RHO12(:)=POTPAR(ISP1)%TAILED%AEF(:,LN1) &
       &          *POTPAR(ISP1)%TAILED%AEF(:,LN2)
          DO LN3=1,LNX1
            L3=POTPAR(ISP1)%TAILED%LOX(LN3)
            DO LN4=1,LNX2
              L4=POTPAR(ISP2)%TAILED%LOX(LN4)
              PHI4=POTPAR(ISP2)%TAILED%AEF(:,LN4)
              TOL=TOLERANCE*TOLFAC1(LN1)*TOLFAC1(LN2) &
       &                   *TOLFAC1(LN3)*TOLFAC2(LN4) 
              TOL=MAX(TOLMIN,TOL)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                CALL RADIAL$POISSON(GID1,NR1,LR1,RHO12,POT12)
                A123(:)=POT12(:)*POTPAR(ISP1)%TAILED%AEF(:,LN3)
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
                    IND=IND+1
                    DO IDIS=1,NDIS
                      DIS=OFFSITEX(ISP1,ISP2)%DIS(IDIS)
                      COUNT=COUNT+1
                      IF(MOD(COUNT-1,NTASKS).NE.THISTASK-1) CYCLE
                      CALL LMTO_TWOCENTER(LR2,MABS,GID1,NR1,A123 &
                                         ,L4,MABS,GID2,NR2,PHI4 &
       &                                 ,DIS,TOL,INTEGRAL)
                      OFFSITEX(ISP1,ISP2)%X31(IDIS,IND)=INTEGRAL
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL MPE$COMBINE('MONOMER','+',OFFSITEX(ISP1,ISP2)%X31)
!
!!$PRINT*,'X31  FOR ',ISP1,ISP2
!!$DO IDIS=1,NDIS
!!$PRINT*,OFFSITEX(ISP1,ISP2)%DIS(IDIS),OFFSITEX(ISP1,ISP2)%X31(IDIS,1) !/(4.D0*PI)
!!$ENDDO
!!$STOP 'FORCED'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEXCONVERT()
!     **************************************************************************
!     ** CONVERT THE DISTANCE-DEPENDENT MATRIX ELEMENTS OF THE U-TENSOR       **
!     ** INTO EXPANSION COEFFICIENTS FOR THE INTERPOLATING FUNCTION           **
!     ** F_J(X)=SUM_I=1^NDIS X22(I,J)*EXP(-LAMBDA(I)*X)                       **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : OFFSITEX,POTPAR,NSP
      IMPLICIT NONE
      INTEGER(4)            :: ISP1,ISP2
      INTEGER(4)            :: NDIS   !#(GRID POINTS)
      INTEGER(4)            :: NF     !#(INTERPOLATING FUNCTIONS)
      INTEGER(4)            :: I
      INTEGER(4)            :: NIND   !#(INTEGRALS)
      REAL(8)   ,ALLOCATABLE:: AMAT(:,:)
      REAL(8)   ,ALLOCATABLE:: AINV(:,:)
      REAL(8)   ,ALLOCATABLE:: AMATIN(:,:)
      REAL(8)   ,ALLOCATABLE:: G(:,:)  !INTERPOLATING FUNCTIONS
REAL(8) :: Y(10),R
INTEGER(4) :: J
      LOGICAL(4)           :: TFIT
!     **************************************************************************
      DO ISP1=1,NSP
        DO ISP2=1,NSP
          NDIS=OFFSITEX(ISP1,ISP2)%NDIS
          NF=OFFSITEX(ISP1,ISP2)%NF
          ALLOCATE(G(NDIS,NF))  ! INTERPOLATING FUNCTIONS
          DO I=1,NF
            G(:,I)=EXP(-OFFSITEX(ISP1,ISP2)%LAMBDA(I) &
     &                 *OFFSITEX(ISP1,ISP2)%DIS(:))
          ENDDO
          TFIT=OFFSITEX(ISP1,ISP2)%NF.NE.OFFSITEX(ISP1,ISP2)%NDIS
          ALLOCATE(AMAT(NF,NF))
          IF(TFIT) THEN
            DO I=1,NF
              DO J=I,NF
                AMAT(I,J)=SUM(G(:,I)*G(:,J))
                AMAT(J,I)=AMAT(I,J)
              ENDDO
            ENDDO
          ELSE
            DO I=1,NF
              AMAT(:,I)=G(:,I)
            ENDDO
          END IF
!
          ALLOCATE(AINV(NF,NF))
          CALL LIB$INVERTR8(NF,AMAT,AINV)
          DEALLOCATE(AMAT)
!
          ALLOCATE(AMATIN(NF,NDIS))
          IF(TFIT) THEN
            AMATIN(:,:)=MATMUL(AINV,TRANSPOSE(G))
          ELSE
            AMATIN=AINV
          END IF
          DEALLOCATE(AINV)
          DEALLOCATE(G)
!
!         ======================================================================
!         == OFF SITE OVERLAP MATRIX 
!         ======================================================================
          NIND=SIZE(OFFSITEX(ISP1,ISP2)%OVERLAP(1,:))
!!$PRINT*,'OFFSITE%DIS    ',OFFSITEX(ISP1,ISP2)%DIS(:)
!!$DO I=1,1
!!$  PRINT*,'OFFSITE%OVERLAP ',OFFSITEX(ISP1,ISP2)%OVERLAP(:,I)
!!$ENDDO
          DO I=1,NIND
            OFFSITEX(ISP1,ISP2)%OVERLAP(:,I) &
     &                         =MATMUL(AMATIN,OFFSITEX(ISP1,ISP2)%OVERLAP(:,I))
          ENDDO
          OFFSITEX(ISP1,ISP2)%OVERLAP(NF+1:,:)=0.D0
!
!         ======================================================================
!         == COULOMB TERMS X22
!         ======================================================================
          NIND=SIZE(OFFSITEX(ISP1,ISP2)%X22(1,:))
!!$DO I=1,1
!!$  PRINT*,'OFFSITE%X22 ',OFFSITEX(ISP1,ISP2)%X22(:,I) &
!!$&   +POTPAR(ISP1)%TAILED%QLN(1,1,1)*POTPAR(ISP2)%TAILED%QLN(1,1,1) &
!!$&   /OFFSITEX(ISP1,ISP2)%DIS(:)
!!$ENDDO
!!$OPEN(UNIT=11,FILE='X22.DAT')
!!$DO I=1,OFFSITEX(ISP1,ISP2)%NDIS
!!$ R=OFFSITEX(ISP1,ISP2)%DIS(I)
!!$ DO J=1,10
!!$   Y(J)=OFFSITEX(ISP1,ISP2)%X22(I,J)
!!$ ENDDO
!!$ WRITE(11,*)R,Y
!!$ENDDO
          DO I=1,NIND
            OFFSITEX(ISP1,ISP2)%X22(:,I)=MATMUL(AMATIN,OFFSITEX(ISP1,ISP2)%X22(:,I))
          ENDDO
          OFFSITEX(ISP1,ISP2)%X22(NF+1:,:)=0.D0
!!$DO I=200,1,-1
!!$  R=REAL(I-1)/REAL(200-1)*15.D0
!!$  IF(R.LT.OFFSITEX(ISP1,ISP2)%DIS(1)) CYCLE
!!$  DO J=1,10
!!$    Y(J)=SUM(OFFSITEX(ISP1,ISP2)%X22(:NF,J)*EXP(-OFFSITEX(ISP1,ISP2)%LAMBDA(:)*R))
!!$  ENDDO
!!$ WRITE(11,*)R,Y
!!$ENDDO
!!$CLOSE(11)
!!$PRINT*,'OFFSITE%LAMBDA ',OFFSITEX(ISP1,ISP2)%LAMBDA(:)
!!$DO I=1,10
!!$  PRINT*,'OFFSITE%X22 ',OFFSITEX(ISP1,ISP2)%X22(:,I)
!!$ENDDO
!
!         ======================================================================
!         == 31 TERMS  X31
!         ======================================================================
          NIND=SIZE(OFFSITEX(ISP1,ISP2)%X31(1,:))
!!$PRINT*,'OFFSITE%DIS    ',OFFSITEX(ISP1,ISP2)%DIS(:)
!!$DO I=1,1
!!$  PRINT*,'OFFSITE%X31 ',OFFSITEX(ISP1,ISP2)%X31(:,I)
!!$ENDDO
!!$OPEN(UNIT=11,FILE='X31.DAT')
!!$DO I=1,NDIS
!!$ R=OFFSITEX(ISP1,ISP2)%DIS(I)
!!$ DO J=1,10
!!$   Y(J)=OFFSITEX(ISP1,ISP2)%X31(I,J)
!!$ ENDDO
!!$ WRITE(11,*)R,Y
!!$ENDDO
          DO I=1,NIND
            OFFSITEX(ISP1,ISP2)%X31(:,I)=MATMUL(AMATIN,OFFSITEX(ISP1,ISP2)%X31(:,I))
          ENDDO
          OFFSITEX(ISP1,ISP2)%X31(NF+1:,:)=0.D0
!!$DO I=200,1,-1
!!$  R=REAL(I-1)/REAL(200-1)*15.D0
!!$  IF(R.LT.OFFSITEX(ISP1,ISP2)%DIS(1)) CYCLE
!!$  DO J=1,10
!!$    Y(J)=SUM(OFFSITEX(ISP1,ISP2)%X31(:NF,J)*EXP(-OFFSITEX(ISP1,ISP2)%LAMBDA(:)*R))
!!$  ENDDO
!!$ WRITE(11,*)R,Y
!!$ENDDO
!!$CLOSE(11)
!!$PRINT*,'OFFSITE%LAMBDA ',OFFSITEX(ISP1,ISP2)%LAMBDA(:)
!!$DO I=1,10
!!$  PRINT*,'OFFSITE%X31 ',OFFSITEX(ISP1,ISP2)%X31(:NF,I)
!!$ENDDO
!!$STOP 'FORCED'
!
!         ======================================================================
!         == BOND OVERLAP INTERACTIONS BONDU                                  ==
!         ======================================================================
          NIND=SIZE(OFFSITEX(ISP1,ISP2)%BONDU(1,:))
!!$DO I=1,1
!!$  PRINT*,'OFFSITE%BONDU ',OFFSITEX(ISP1,ISP2)%BONDU(:,I)
!!$ENDDO
!!$OPEN(UNIT=11,FILE='BONDU.DAT')
!!$DO I=1,NDIS
!!$ R=OFFSITEX(ISP1,ISP2)%DIS(I)
!!$ DO J=1,MIN(10,NIND)
!!$   Y(J)=OFFSITEX(ISP1,ISP2)%BONDU(I,J)
!!$ ENDDO
!!$ WRITE(11,*)R,Y
!!$ENDDO

          DO I=1,NIND
            OFFSITEX(ISP1,ISP2)%BONDU(:,I)=MATMUL(AMATIN &
     &                           ,OFFSITEX(ISP1,ISP2)%BONDU(:,I))
          ENDDO
          OFFSITEX(ISP1,ISP2)%BONDU(NF+1:,:)=0.D0
!!$DO I=200,1,-1
!!$  R=REAL(I-1)/REAL(200-1)*15.D0
!!$  IF(R.LT.OFFSITEX(ISP1,ISP2)%DIS(1)) CYCLE
!!$  DO J=1,MIN(10,NIND)
!!$    Y(J)=SUM(OFFSITEX(ISP1,ISP2)%BONDU(:NF,J)*EXP(-OFFSITEX(ISP1,ISP2)%LAMBDA(:)*R))
!!$  ENDDO
!!$ WRITE(11,*)R,Y
!!$ENDDO
!!$CLOSE(11)
!!$PRINT*,'OFFSITE%LAMBDA  ',OFFSITEX(ISP1,ISP2)%LAMBDA(:)
!!$DO I=1,10
!!$  PRINT*,'OFFSITE%BONDU ',OFFSITEX(ISP1,ISP2)%BONDU(:NF,I)
!!$ENDDO
!!$STOP 'FORCED'

!
!         ======================================================================
!         ==                                                                  ==
!         ======================================================================
          DEALLOCATE(AMATIN)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PRBONDU()
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR,OFFSITEX
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: ISP1=1
      INTEGER(4),PARAMETER :: ISP2=1
      INTEGER(4)           :: IND
      INTEGER(4)           :: LMN1,LMN2,LMN3,LMN4
      INTEGER(4)           :: LMNXA,LMNXB
      INTEGER(4)           :: IAB,ICD
      REAL(8)   ,ALLOCATABLE :: UABCD(:,:,:,:)
!     **************************************************************************
      LMNXA=POTPAR(ISP1)%TAILED%LMNX
      LMNXB=POTPAR(ISP2)%TAILED%LMNX
      ALLOCATE(UABCD(LMNXA,LMNXB,LMNXA,LMNXB))
      UABCD=0.D0
      IND=0
      DO LMN2=1,LMNXB
        DO LMN1=1,LMNXA
          IAB=LMN1+LMNXA*(LMN2-1)
          DO LMN4=1,LMNXB
            DO LMN3=1,LMNXA
              ICD=LMN3+LMNXA*(LMN4-1)
              IF(ICD.GT.IAB) EXIT
              IND=IND+1
              UABCD(LMN1,LMN2,LMN3,LMN4)=OFFSITEX(ISP1,ISP2)%BONDU(1,IND)
              UABCD(LMN3,LMN4,LMN1,LMN2)=OFFSITEX(ISP1,ISP2)%BONDU(1,IND)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO LMN2=1,LMNXB
        DO LMN1=1,LMNXA
          WRITE(*,FMT='(82("="),T20," BONDU FOR LMN1=",I3," AND LMN2=",I3)')LMN1,LMN2
          DO LMN3=1,LMNXA
            WRITE(*,FMT='(I3,50F10.5)')LMN3,UABCD(LMN1,LMN2,LMN3,:)
          ENDDO
        ENDDO
      ENDDO
STOP 'FORCED'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEOVERLAP(ISP1,ISP2,DIS,LMNX1,LMNX2,O,DO)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS                                               **
!     ** TWO ORBITALS ARE ON THE FIRST AND TWO ORBITALS ON THE SECOND ATOM    **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      REAL(8)   ,INTENT(IN) :: DIS
      INTEGER(4),INTENT(IN) :: LMNX1
      INTEGER(4),INTENT(IN) :: LMNX2
      REAL(8)   ,INTENT(OUT):: O(LMNX1,LMNX2)
      REAL(8)   ,INTENT(OUT):: DO(LMNX1,LMNX2)
      INTEGER(4)            :: LRX1,LRX2
      INTEGER(4)            :: LNX1,LNX2
      INTEGER(4)            :: IND
      INTEGER(4)            :: LN1,LN2,LN3
      INTEGER(4)            :: L1,L2,L3
      INTEGER(4)            :: MABS,PM,M
      INTEGER(4),ALLOCATABLE:: LMN0A(:),LMN0B(:)
      INTEGER(4)            :: LMN1,LMN2
      REAL(8)               :: INTEGRAL,DINTEGRAL
      INTEGER(4)            :: LMRX
!     **************************************************************************
!
!     ==========================================================================
!     == PREPARATION                                                          ==
!     ==========================================================================
      LNX1=POTPAR(ISP1)%TAILED%LNX
      LNX2=POTPAR(ISP2)%TAILED%LNX
      CALL SETUP$ISELECT(ISP1)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX1=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      CALL SETUP$ISELECT(ISP2)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX2=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      ALLOCATE(LMN0A(LNX1))
      ALLOCATE(LMN0B(LNX2))
      LMN0A(1)=0
      DO LN1=1,LNX1-1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        LMN0A(LN1+1)=LMN0A(LN1)+2*L1+1
      ENDDO                
      LMN0B(1)=0
      DO LN3=1,LNX2-1
        L3=POTPAR(ISP2)%TAILED%LOX(LN3)
        LMN0B(LN3+1)=LMN0B(LN3)+2*L3+1
      ENDDO                
!
!     ==========================================================================
!     == OBTAIN U-TENSOR                                                      ==
!     ==========================================================================
      O(:,:)=0.D0
      DO(:,:)=0.D0
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        DO LN2=1,LNX2
          L2=POTPAR(ISP2)%TAILED%LOX(LN2)
          DO MABS=0,MIN(L1,L2)
            IND=IND+1
!
!           == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE... ================
            CALL LMTO_OFFSITEXVALUE('OV',ISP1,ISP2,IND,DIS,INTEGRAL,DINTEGRAL)
!
            DO PM=-1,1,2
              IF(MABS.EQ.0.AND.PM.EQ.1) CYCLE
              M=MABS*PM
              LMN1=LMN0A(LN1)+L1+1+M
              LMN2=LMN0B(LN2)+L2+1+M
              O(LMN1,LMN2) =O(LMN1,LMN2) +INTEGRAL
              DO(LMN1,LMN2)=DO(LMN1,LMN2)+DINTEGRAL
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEX22U(ISP1,ISP2,DIS,LMNX1,LMNX2,U,DU)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS                                               **
!     ** TWO ORBITALS ARE ON THE FIRST AND TWO ORBITALS ON THE SECOND ATOM    **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      REAL(8)   ,INTENT(IN) :: DIS
      INTEGER(4),INTENT(IN) :: LMNX1
      INTEGER(4),INTENT(IN) :: LMNX2
      REAL(8)   ,INTENT(OUT):: U(LMNX1,LMNX1,LMNX2,LMNX2)
      REAL(8)   ,INTENT(OUT):: DU(LMNX1,LMNX1,LMNX2,LMNX2)
      INTEGER(4)            :: LRX1,LRX2
      INTEGER(4)            :: LNX1,LNX2
      INTEGER(4)            :: IND
      INTEGER(4)            :: LN1,LN2,LN3,LN4
      INTEGER(4)            :: L1,L2,L3,L4
      INTEGER(4)            :: LR1,LR2
      INTEGER(4)            :: LMR1,LMR2
      INTEGER(4)            :: MABS,PM,M
      INTEGER(4),ALLOCATABLE:: LMN0A(:),LMN0B(:)
      INTEGER(4)            :: LM1,LM2,LM3,LM4
      INTEGER(4)            :: LMN1,LMN2,LMN3,LMN4
      INTEGER(4)            :: M1,M2,M3,M4
      REAL(8)               :: INTEGRAL,DINTEGRAL
      REAL(8)               :: CG1,CG2 ! GAUNT COEFFICIENTS
      INTEGER(4)            :: LMRX
      REAL(8)               :: PI,SQ4PI,SQ4PITHIRD
      REAL(8)               :: SVAR
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      SQ4PI=SQRT(4.D0*PI)
      SQ4PITHIRD=SQRT(4.D0*PI/3.D0)
!
!     ==========================================================================
!     == PREPARATION                                                          ==
!     ==========================================================================
      LNX1=POTPAR(ISP1)%TAILED%LNX
      LNX2=POTPAR(ISP2)%TAILED%LNX
      CALL SETUP$ISELECT(ISP1)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX1=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      CALL SETUP$ISELECT(ISP2)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX2=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      ALLOCATE(LMN0A(LNX1))
      ALLOCATE(LMN0B(LNX2))
      LMN0A(1)=0
      DO LN1=1,LNX1-1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        LMN0A(LN1+1)=LMN0A(LN1)+2*L1+1
      ENDDO                
      LMN0B(1)=0
      DO LN3=1,LNX2-1
        L3=POTPAR(ISP2)%TAILED%LOX(LN3)
        LMN0B(LN3+1)=LMN0B(LN3)+2*L3+1
      ENDDO                
!
!     ==========================================================================
!     == OBTAIN U-TENSOR                                                      ==
!     ==========================================================================
      U(:,:,:,:)=0.D0
      DU(:,:,:,:)=0.D0
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%TAILED%LOX(LN2)
          DO LN3=1,LNX2
            L3=POTPAR(ISP2)%TAILED%LOX(LN3)
            DO LN4=LN3,LNX2
              L4=POTPAR(ISP2)%TAILED%LOX(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                DO LR2=ABS(L3-L4),MIN(L3+L4,LRX2),2
                  DO MABS=0,MIN(LR1,LR2)
                    IND=IND+1
!
!                   == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE...
                    CALL LMTO_OFFSITEXVALUE('22',ISP1,ISP2,IND,DIS &
     &                                          ,INTEGRAL,DINTEGRAL)
                    IF(LR1.LE.1.AND.LR2.LE.1) THEN
                      SVAR=POTPAR(ISP1)%TAILED%QLN(LR1+1,LN1,LN2) &
     &                    *POTPAR(ISP2)%TAILED%QLN(LR2+1,LN3,LN4)
                      SVAR=SVAR/DIS**(LR1+LR2+1)
                      INTEGRAL =INTEGRAL +SVAR
                      SVAR=-REAL(LR1+LR2+1)*SVAR/DIS
                      DINTEGRAL=DINTEGRAL+SVAR
                    END IF
!PRINT*,'INTEGRAL ',INTEGRAL/(4.D0*PI)
!
                    DO PM=-1,1,2
                      IF(MABS.EQ.0.AND.PM.EQ.1) CYCLE
                      M=MABS*PM
                      LMR1=LR1**2+LR1+1+M
                      LMR2=LR2**2+LR2+1+M

                      LM1=L1**2
                      LMN1=LMN0A(LN1)
                      DO M1=1,2*L1+1
                        LM1=LM1+1
                        LMN1=LMN1+1
                        LM2=L2**2
                        LMN2=LMN0A(LN2)
                        DO M2=1,2*L2+1
                          LM2=LM2+1
                          LMN2=LMN2+1
                          CALL SPHERICAL$GAUNT(LM1,LM2,LMR1,CG1)
                          IF(CG1.EQ.0.D0) CYCLE
!
                          LM3=L3**2
                          LMN3=LMN0B(LN3)
                          DO M3=1,2*L3+1
                            LM3=LM3+1
                            LMN3=LMN3+1
                            LM4=L4**2
                            LMN4=LMN0B(LN4)
                            DO M4=1,2*L4+1
                              LM4=LM4+1
                              LMN4=LMN4+1
!
                              CALL SPHERICAL$GAUNT(LM3,LM4,LMR2,CG2)
                              U(LMN1,LMN2,LMN3,LMN4) =U(LMN1,LMN2,LMN3,LMN4) &
       &                                             +CG1*CG2*INTEGRAL
                              DU(LMN1,LMN2,LMN3,LMN4)=DU(LMN1,LMN2,LMN3,LMN4) &
       &                                             +CG1*CG2*DINTEGRAL
                            ENDDO
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!PRINT*,'U(1,1,1,1)',ISP1,ISP2,U(1,1,1,1)
!
!     ==========================================================================
!     == COMPLETE MISSING MATRIX ELEMENTS USING SYMMETRY =======================
!     ==========================================================================
      DO LMN1=1,LMNX1
        DO LMN2=LMN1+1,LMNX1
          U(LMN2,LMN1,:,:)=U(LMN1,LMN2,:,:)
        ENDDO
      ENDDO
      DO LMN3=1,LMNX2
        DO LMN4=LMN3+1,LMNX2
          U(:,:,LMN4,LMN3)=U(:,:,LMN3,LMN4)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEX31U(ISP1,ISP2,DIS,LMNX1,LMNX2,U,DU)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS 
!     ** TWO ORBITALS ARE ON THE FIRST AND TWO ORBITALS ON THE SECOND ATOM    **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      REAL(8)   ,INTENT(IN) :: DIS
      INTEGER(4),INTENT(IN) :: LMNX1
      INTEGER(4),INTENT(IN) :: LMNX2
      REAL(8)   ,INTENT(OUT):: U(LMNX1,LMNX1,LMNX1,LMNX2)
      REAL(8)   ,INTENT(OUT):: DU(LMNX1,LMNX1,LMNX1,LMNX2)
      INTEGER(4)            :: LRX1,LRX2
      INTEGER(4)            :: LNX1,LNX2
      INTEGER(4)            :: IND
      INTEGER(4)            :: LN1,LN2,LN3,LN4
      INTEGER(4)            :: L1,L2,L3,L4
      INTEGER(4)            :: LR1,LR2
      INTEGER(4)            :: LMR1,LMR2
      INTEGER(4)            :: MABS,PM,M
      INTEGER(4),ALLOCATABLE:: LMN0A(:),LMN0B(:)
      INTEGER(4)            :: LM1,LM2,LM3,LM4
      INTEGER(4)            :: LMN1,LMN2,LMN3,LMN4
      INTEGER(4)            :: M1,M2,M3,M4,MR1
      REAL(8)               :: INTEGRAL,DINTEGRAL
      REAL(8)               :: SVAR
      INTEGER(4)            :: LMRX
      REAL(8)               :: CG1,CG2 ! GAUNT COEFFICIENTS
!     **************************************************************************
!
!     ==========================================================================
!     == PREPARATION                                                          ==
!     ==========================================================================
      LNX1=POTPAR(ISP1)%TAILED%LNX
      LNX2=POTPAR(ISP2)%TAILED%LNX
      CALL SETUP$ISELECT(ISP1)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX1=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      CALL SETUP$ISELECT(ISP2)
      CALL SETUP$GETI4('LMRX',LMRX)
      CALL SETUP$UNSELECT()
      LRX2=INT(SQRT(REAL(LMRX-1,KIND=8)+1.D-9))
      ALLOCATE(LMN0A(LNX1))
      ALLOCATE(LMN0B(LNX2))
      LMN0A(1)=0
      DO LN1=1,LNX1-1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        LMN0A(LN1+1)=LMN0A(LN1)+2*L1+1
      ENDDO                
      LMN0B(1)=0
      DO LN3=1,LNX2-1
        L3=POTPAR(ISP2)%TAILED%LOX(LN3)
        LMN0B(LN3+1)=LMN0B(LN3)+2*L3+1
      ENDDO                
!
!     ==========================================================================
!     == OBTAIN U-TENSOR                                                      ==
!     ==========================================================================
      U(:,:,:,:)=0.D0
      DU(:,:,:,:)=0.D0
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%TAILED%LOX(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%TAILED%LOX(LN2)
          DO LN3=1,LNX1
            L3=POTPAR(ISP1)%TAILED%LOX(LN3)
            DO LN4=1,LNX2
              L4=POTPAR(ISP2)%TAILED%LOX(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
                    IND=IND+1
!
!                   == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE...
                    CALL LMTO_OFFSITEXVALUE('31',ISP1,ISP2,IND,ABS(DIS) &
     &                                          ,INTEGRAL,DINTEGRAL)
                    IF(DIS.LT.0.D0) THEN
                      SVAR=(-1.D0)**(L1+L2+L3+L4)
                      INTEGRAL=SVAR*INTEGRAL
                      DINTEGRAL=-SVAR*DINTEGRAL
                    END IF
!
                    DO PM=-1,1,2
                      IF(MABS.EQ.0.AND.PM.EQ.1) CYCLE
                      M=MABS*PM

                      LM1=L1**2
                      LMN1=LMN0A(LN1)
                      DO M1=1,2*L1+1
                        LM1=LM1+1
                        LMN1=LMN1+1
                        LM2=L2**2
                        LMN2=LMN0A(LN2)
                        DO M2=1,2*L2+1
                          LM2=LM2+1
                          LMN2=LMN2+1
                          LMR1=LR1**2
                          DO MR1=1,2*LR1+1
                            LMR1=LMR1+1
                            CALL SPHERICAL$GAUNT(LM1,LM2,LMR1,CG1)
                            IF(CG1.EQ.0.D0) CYCLE
!
                            LM3=L3**2
                            LMN3=LMN0A(LN3)
                            DO M3=1,2*L3+1
                              LM3=LM3+1
                              LMN3=LMN3+1
                              LMR2=LR2**2+LR2+1+M
                              CALL SPHERICAL$GAUNT(LMR1,LM3,LMR2,CG2)
                              IF(CG2.EQ.0.D0) CYCLE
!
                              LM4=L4**2+L4+1+M
                              LMN4=LMN0B(LN4)+L4+1+M
!
                              U(LMN1,LMN2,LMN3,LMN4) =U(LMN1,LMN2,LMN3,LMN4) &
       &                                             +CG1*CG2*INTEGRAL
                              DU(LMN1,LMN2,LMN3,LMN4)=DU(LMN1,LMN2,LMN3,LMN4) &
       &                                             +CG1*CG2*DINTEGRAL
                            ENDDO
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == COMPLETE MISSING MATRIX ELEMENTS USING SYMMETRY =======================
!     ==========================================================================
      DO LMN1=1,LMNX1
        DO LMN2=LMN1+1,LMNX1
          U(LMN2,LMN1,:,:)=U(LMN1,LMN2,:,:)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEXBONDU(ISP1,ISP2,DIS,LMNX1,LMNX2,U,DU)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS 
!     ** TWO ORBITALS ARE ON THE FIRST AND TWO ORBITALS ON THE SECOND ATOM    **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      REAL(8)   ,INTENT(IN) :: DIS
      INTEGER(4),INTENT(IN) :: LMNX1
      INTEGER(4),INTENT(IN) :: LMNX2
      REAL(8)   ,INTENT(OUT):: U(LMNX1,LMNX2,LMNX1,LMNX2)
      REAL(8)   ,INTENT(OUT):: DU(LMNX1,LMNX2,LMNX1,LMNX2)
      INTEGER(4)            :: IND
      INTEGER(4)            :: LMN1,LMN2,LMN3,LMN4
      INTEGER(4)            :: I12,I34
      REAL(8)               :: INTEGRAL,DINTEGRAL
      INTEGER(4)            :: LMRX
!     **************************************************************************
!
!     ==========================================================================
!     == OBTAIN U-TENSOR                                                      ==
!     ==========================================================================
      U(:,:,:,:)=0.D0
      DU(:,:,:,:)=0.D0
      IND=0
      DO LMN2=1,LMNX2
        DO LMN1=1,LMNX1
          I12=LMN1+LMNX1*(LMN2-1)
          I34=0
          DO LMN4=1,LMNX2
            DO LMN3=1,LMNX1
              I34=LMN3+LMNX1*(LMN4-1)
              IF(I34.GT.I12) EXIT
              IND=IND+1
!
!             == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE...
              CALL LMTO_OFFSITEXVALUE('BONDU',ISP1,ISP2,IND,ABS(DIS) &
     &                                          ,INTEGRAL,DINTEGRAL)
              U(LMN1,LMN2,LMN3,LMN4) =INTEGRAL
              U(LMN3,LMN4,LMN1,LMN2) =INTEGRAL
              DU(LMN1,LMN2,LMN3,LMN4)=DINTEGRAL
              DU(LMN3,LMN4,LMN1,LMN2)=DINTEGRAL
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEXVALUE(ID,ISP1,ISP2,IND,DIS,INTEGRAL,DINTEGRAL)
!     **************************************************************************
!     ** INTERPOLATE VALUES TO ACTUAL DISTANCE
!     **************************************************************************
      USE LMTO_MODULE, ONLY  : OFFSITEX
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: ISP1
      INTEGER(4)  ,INTENT(IN) :: ISP2
      INTEGER(4)  ,INTENT(IN) :: IND
      REAL(8)     ,INTENT(IN) :: DIS
      REAL(8)     ,INTENT(OUT):: INTEGRAL
      REAL(8)     ,INTENT(OUT):: DINTEGRAL
      INTEGER(4)  ,SAVE       :: NF
      INTEGER(4)  ,PARAMETER  :: LF=20
      INTEGER(4)  ,SAVE       :: ISP1SAVE=0,ISP2SAVE=0
      REAL(8)     ,SAVE       :: DISSAVE=0.D0
      REAL(8)     ,SAVE       :: G(LF)
!     **************************************************************************
      INTEGRAL=0.D0
      IF(ISP1.NE.ISP1SAVE.OR.ISP2.NE.ISP2SAVE.OR.DIS.NE.DISSAVE) THEN
        ISP1SAVE=ISP1
        ISP2SAVE=ISP2
        NF=OFFSITEX(ISP1,ISP2)%NF
        IF(NF.GT.LF) THEN
          CALL ERROR$MSG('INTERNAL ARRAY SIZE EXCEEDED')
          CALL ERROR$STOP('LMTO_OFFSITEXVALUE')
        END IF
        DISSAVE=DIS      
        G(:NF)=EXP(-OFFSITEX(ISP1,ISP2)%LAMBDA(:)*DIS)
      END IF
      IF(ID.EQ.'22') THEN
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%X22(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%X22(:NF,IND)*G(:NF))
      ELSE IF (ID.EQ.'31') THEN
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%X31(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%X31(:NF,IND)*G(:NF))
      ELSE IF (ID.EQ.'BONDU') THEN
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%BONDU(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%BONDU(:NF,IND)*G(:NF))
      ELSE IF (ID.EQ.'OV') THEN
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%OVERLAP(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%OVERLAP(:NF,IND)*G(:NF))
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('LMTO_OFFSITEXVALUE')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TWOCENTER_TEST()
!     **************************************************************************
!     ** TEST SUBROUTINE "LMTO_TWOCENTER"
!     ** TESTS FROM Z. ROMANOWSKI, INT.J.QUANT. CHEM. 108, 249 (2008)         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)           :: GID
      INTEGER(4),PARAMETER :: NR=250
      REAL(8)              :: R(NR)
      REAL(8)              :: F(NR),G(NR)
      REAL(8)              :: AUX(NR)
      REAL(8)              :: VAL
      INTEGER(4)           :: N1,N2
      REAL(8)              :: SIGMA1,SIGMA2
      INTEGER(4)           :: L1,L2,M1,M2
      REAL(8)              :: DIS
      REAL(8)              :: OVERLAP,SVAR
      REAL(8)              :: TOL
      INTEGER(4)           :: I
!     **************************************************************************
      CALL RADIAL$NEW('SHLOG',GID)
      CALL RADIAL$SETI4(GID,'NR',NR)
      CALL RADIAL$SETR8(GID,'DEX',0.05D0)
      CALL RADIAL$SETR8(GID,'R1',1.D-4)
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
      N1=1; SIGMA1=3.D0 ;L1=0; M1=0
      N2=1; SIGMA2=1.D0 ;L2=0; M2=0
      CALL RADIAL$STO(GID,NR,N1,SIGMA1,F)
      CALL RADIAL$STO(GID,NR,N2,SIGMA2,G)
!     == CHECK NORM
      AUX(:)=R(:)**2*F(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
PRINT*,'NORM OF F=',VAL
      AUX(:)=R(:)**2*G(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
PRINT*,'NORM OF G=',VAL
      TOL=1.D-6
      DO I=1,20
        DIS=0.2D0*REAL(I,KIND=8)
        CALL LMTO_TWOCENTER(L1,M1,GID,NR,F,L2,M2,GID,NR,G,DIS,TOL,OVERLAP)
        SVAR=3.D0*SQRT(3.D0)/(16.D0*DIS)*EXP(-3.D0*DIS)*(3.D0+2.D0*DIS-EXP(2.D0*DIS)*(3.D0-6.D0*DIS))
        PRINT*,'DIS=',DIS,'OVERLAP ',OVERLAP,SVAR,OVERLAP-SVAR
      ENDDO
!     ==========================================================================
      N1=3; SIGMA1=4.D0 ;L1=2; M1=0
      N2=3; SIGMA2=4.D0 ;L2=2; M2=0
      CALL RADIAL$STO(GID,NR,N1,SIGMA1,F)
      CALL RADIAL$STO(GID,NR,N2,SIGMA2,G)
!     == CHECK NORM
      AUX(:)=R(:)**2*F(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
PRINT*,'NORM OF F=',VAL
      AUX(:)=R(:)**2*G(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
PRINT*,'NORM OF G=',VAL
      TOL=1.D-6
      DO I=1,20
        DIS=0.2D0*REAL(I,KIND=8)
        CALL LMTO_TWOCENTER(L1,M1,GID,NR,F,L2,M2,GID,NR,G,DIS,TOL,OVERLAP)
        SVAR=1.D0/315.D0*EXP(-4.D0*DIS)*(315.D0+4.D0*DIS*(315.D0+4.D0*DIS &
    &       *(75.D0+8.D0*DIS*(-15.D0+4.D0*DIS*(-9.D0-2.D0*DIS+8.D0*DIS**2)))))
        PRINT*,'DIS=',DIS,'OVERLAP ',OVERLAP,SVAR,OVERLAP-SVAR
      ENDDO
!
STOP 'FORCED'
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE LMTO_TWOCENTER_MODULE
LOGICAL(4)         :: TPR=.FALSE.
REAL(8)            :: DIS
INTEGER(4)         :: L1
INTEGER(4)         :: M1
INTEGER(4)         :: GID1
INTEGER(4)         :: NR1
REAL(8)   ,POINTER :: F1(:)
INTEGER(4)         :: L2
INTEGER(4)         :: M2
INTEGER(4)         :: GID2
INTEGER(4)         :: NR2
REAL(8)   ,POINTER :: F2(:)
REAL(8)   ,POINTER :: PLMWORK1(:)  ! WORK ARRAYS
REAL(8)   ,POINTER :: PLMWORK2(:)
REAL(8)            :: SLEN ! PARAMETER DEFINING INTEGRATION AREA
ENDMODULE LMTO_TWOCENTER_MODULE
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE LMTO_TWOCENTER(L1_,M1_,GID1_,NR1_,F1_ &
      &                         ,L2_,M2_,GID2_,NR2_,F2_,DIS_,TOLERANCE,OVERLAP)
!      *************************************************************************
!      *************************************************************************
       USE LMTO_TWOCENTER_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: L1_
       INTEGER(4),INTENT(IN) :: M1_
       INTEGER(4),INTENT(IN) :: GID1_
       INTEGER(4),INTENT(IN) :: NR1_
       REAL(8)   ,INTENT(IN) :: F1_(NR1_)
       INTEGER(4),INTENT(IN) :: L2_
       INTEGER(4),INTENT(IN) :: M2_
       INTEGER(4),INTENT(IN) :: GID2_
       INTEGER(4),INTENT(IN) :: NR2_
       REAL(8)   ,INTENT(IN) :: F2_(NR2_)
       REAL(8)   ,INTENT(IN) :: DIS_
       REAL(8)   ,INTENT(IN) :: TOLERANCE
       REAL(8)   ,INTENT(OUT):: OVERLAP
       REAL(8)               :: SVAR1,SVAR2
!      *************************************************************************
       IF(ABS(M1_).GT.L1_.OR.ABS(M2_).GT.L2_) THEN
         CALL ERROR$MSG('MAGNETIC QUANTUM NUMBER OUT OF RANGE (|M| > L)')
         CALL ERROR$STOP('LMTO$TWOCENTER')
       END IF
!
!      =========================================================================
!      == RESULT IS DIAGONAL IN M                                             ==
!      =========================================================================
       OVERLAP=0.D0
       IF(M1_.NE.M2_) RETURN
!
!      =========================================================================
!      == COPY DATA INTO MODULE M                                             ==
!      =========================================================================
       ALLOCATE(F1(NR1_))
       ALLOCATE(F2(NR2_))
       ALLOCATE(PLMWORK1((L1_+1)**2))  ! WORK ARRAY
       ALLOCATE(PLMWORK2((L2_+1)**2))  ! WORK ARRAY
!
!      =========================================================================
!      == COPY DATA INTO MODULE                                               ==
!      =========================================================================
       L1=L1_
       M1=M1_
       GID1=GID1_
       NR1=NR1_
       F1(:)=F1_(:)
!
       L2=L2_
       M2=M2_
       GID2=GID2_
       NR2=NR2_
       F2(:)=F2_(:)
!
       DIS=DIS_ 
       CALL RADIAL$GETR8(GID1,'RMAX',SVAR1)
       CALL RADIAL$GETR8(GID2,'RMAX',SVAR2)
       SLEN=0.5D0*(SVAR1+SVAR2-DIS)
!
!      =========================================================================
!      == PERFORM INTEGRATIONLE M                                             ==
!      =========================================================================
       CALL ADAPT$EVALUATE(TOLERANCE,OVERLAP)
!
       DEALLOCATE(PLMWORK1)
       DEALLOCATE(PLMWORK2)
       DEALLOCATE(F1)
       DEALLOCATE(F2)
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE LMTO_TWOCENTER_MYFUNC(P,VALUE)
!      *************************************************************************
!      *************************************************************************
       USE LMTO_TWOCENTER_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: P(2)  ! POINT [-1,1]X[-1,1]
       REAL(8),INTENT(OUT):: VALUE ! RESULTING OVERLAP
       REAL(8)            :: PI,FPI,SQ2
       REAL(8)            :: COSTHETA1,COSTHETA2
       REAL(8)            :: PLM1,PLM2
       REAL(8)            :: R0(2),T1(2),T2(2),R(2)
       REAL(8)            :: R1,R2
       REAL(8)            :: DRDP
       REAL(8)            :: VAL1,VAL2
!      *************************************************************************
       PI=4.D0*ATAN(1.D0)
       FPI=4.D0*PI
       SQ2=SQRT(2.D0)
!
       R0(:)=(/DIS,0.D0/)
       T1(:)=(/-DIS,DIS/)
       T2(:)=(/SLEN,SLEN/)
       R(:)=R0(:)+T1(:)*0.5D0*(1.D0+P(1))+T2(:)*0.5D0*(1.D0+P(2))
       DRDP=0.25D0*ABS(T1(1)*T2(2)-T1(2)*T2(1))
       R1=R(1)
       R2=R(2)
!
       CALL RADIAL$VALUE(GID1,NR1,F1,R1,VAL1)
       CALL RADIAL$VALUE(GID2,NR2,F2,R2,VAL2)
!
       COSTHETA1=(R1**2-R2**2+DIS**2)/(2.D0*DIS*R1)
       COSTHETA2=(R1**2-R2**2-DIS**2)/(2.D0*DIS*R2)
!      == ROUTINES ARE IN PSHERICAL OBJECT
       CALL PLGNDR((L1+1)**2,L1,COSTHETA1,PLMWORK1)
       CALL PLGNDR((L2+1)**2,L2,COSTHETA2,PLMWORK2)
       PLM1=PLMWORK1(L1*(L1+1)+M1+1)
       PLM2=PLMWORK2(L2*(L2+1)+M2+1)
       PLM1=PLM1*SQRT(REAL(2*L1+1,KIND=8)/FPI)
       PLM2=PLM2*SQRT(REAL(2*L2+1,KIND=8)/FPI)
       IF(M1.NE.0)PLM1=PLM1*SQ2
       IF(M1.NE.0)PLM2=PLM1*SQ2
!
       VALUE=2.D0*PI/DIS*DRDP * R1*VAL1*R2*VAL2*PLM1*PLM2
!IF(TPR)WRITE(*,FMT='("==",10F10.5)')DIS,VALUE,R1,VAL1,R2,VAL2,PLM1*PLM2,DRDP
       RETURN
       END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE ADAPT_MODULE
 TYPE SEGMENT_TYPE
   REAL(8) :: ERR
   REAL(8) :: VALUE
   INTEGER :: DIVIDEAXIS
   REAL(8) :: CENTER(2)
   REAL(8) :: WIDTH(2)
END TYPE SEGMENT_TYPE
LOGICAL(4)           :: TINI=.FALSE.
INTEGER(4),PARAMETER :: NP=17
REAL(8)              :: XY(2,NP)
REAL(8)              :: VALW(NP)
REAL(8)              :: ERRW(NP)
REAL(8)              :: DIRW(2,NP)
END MODULE ADAPT_MODULE
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE ADAPTINI()
!      *************************************************************************
!      ** INITIALIZE ADAPT. CALCULATE GRID POINTS AND WEIGHTS                 **
!      *************************************************************************
       USE ADAPT_MODULE
       IMPLICIT NONE
       REAL(8)           :: WP(NP)
       REAL(8)           :: SVAR
       REAL(8),PARAMETER :: LAMBDA2=SQRT(9.D0/70.D0)
       REAL(8),PARAMETER :: LAMBDA3=SQRT(9.D0/10.D0)
       REAL(8),PARAMETER :: LAMBDA4=SQRT(9.D0/10.D0)
       REAL(8),PARAMETER :: LAMBDA5=SQRT(9.D0/19.D0)
       REAL(8),PARAMETER :: W1=-15264.D0/19683.D0
       REAL(8),PARAMETER :: W2=3920.D0/6561.D0    !=2940/19683
       REAL(8),PARAMETER :: W3=4080.D0/19683.D0
       REAL(8),PARAMETER :: W4=800.D0/19683.D0
       REAL(8),PARAMETER :: W5=6859.D0/19683.D0
       REAL(8),PARAMETER :: W1P=-3884.D0/729.D0
       REAL(8),PARAMETER :: W2P=980.D0/486.D0
       REAL(8),PARAMETER :: W3P=260.D0/1458.D0
       REAL(8),PARAMETER :: W4P=100.D0/729.D0
!      *************************************************************************
       IF(TINI) RETURN
       TINI=.TRUE.
       XY(:, 1)=(/    0.D0,    0.D0/); VALW( 1)=W1 ; WP( 1)=W1P  
       XY(:, 2)=(/+LAMBDA2,    0.D0/); VALW( 2)=W2 ; WP( 2)=W2P 
       XY(:, 3)=(/-LAMBDA2,    0.D0/); VALW( 3)=W2 ; WP( 3)=W2P
       XY(:, 4)=(/0.D0    ,+LAMBDA2/); VALW( 4)=W2 ; WP( 4)=W2P
       XY(:, 5)=(/0.D0    ,-LAMBDA2/); VALW( 5)=W2 ; WP( 5)=W2P
       XY(:, 6)=(/+LAMBDA3,    0.D0/); VALW( 6)=W3 ; WP( 6)=W3P
       XY(:, 7)=(/-LAMBDA3,    0.D0/); VALW( 7)=W3 ; WP( 7)=W3P
       XY(:, 8)=(/0.D0    ,+LAMBDA3/); VALW( 8)=W3 ; WP( 8)=W3P
       XY(:, 9)=(/0.D0    ,-LAMBDA3/); VALW( 9)=W3 ; WP( 9)=W3P
       XY(:,10)=(/+LAMBDA4,+LAMBDA4/); VALW(10)=W4 ; WP(10)=W4P
       XY(:,11)=(/+LAMBDA4,-LAMBDA4/); VALW(11)=W4 ; WP(11)=W4P
       XY(:,12)=(/-LAMBDA4,+LAMBDA4/); VALW(12)=W4 ; WP(12)=W4P
       XY(:,13)=(/-LAMBDA4,-LAMBDA4/); VALW(13)=W4 ; WP(13)=W4P
       XY(:,14)=(/+LAMBDA5,+LAMBDA5/); VALW(14)=W5 ; WP(14)=0.D0
       XY(:,15)=(/+LAMBDA5,-LAMBDA5/); VALW(15)=W5 ; WP(15)=0.D0
       XY(:,16)=(/-LAMBDA5,+LAMBDA5/); VALW(16)=W5 ; WP(16)=0.D0
       XY(:,17)=(/-LAMBDA5,-LAMBDA5/); VALW(17)=W5 ; WP(17)=0.D0
!
       ERRW=VALW-WP   ! WEIGHT FOR ERROR CALCULATION
!
       SVAR=(LAMBDA2/LAMBDA3)**2
       DIRW(:,1)=(/1,1/)*(-2.D0)*(1.D0-SVAR)
       DIRW(:,2)=(/1.D0,0.D0/)
       DIRW(:,3)=(/1.D0,0.D0/)
       DIRW(:,4)=(/0.D0,1.D0/)
       DIRW(:,5)=(/0.D0,1.D0/)
       DIRW(:,6)=(/-SVAR,0.D0/)
       DIRW(:,7)=(/-SVAR,0.D0/)
       DIRW(:,8)=(/0.D0,-SVAR/)
       DIRW(:,9)=(/0.D0,-SVAR/)
       DIRW(:,10:)=0.D0
       RETURN
       END       
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE ADAPT_BASICRULE(SEGMENT)
!      *************************************************************************
!      **  PERFORM GAUSSIAN QUADRATURE FOR ONE SEGMENT                        **
!      *************************************************************************
       USE ADAPT_MODULE
       IMPLICIT NONE
       TYPE(SEGMENT_TYPE), INTENT(INOUT) :: SEGMENT
       REAL(8)                           :: SUMVAL
       REAL(8)                           :: SUMERR
       REAL(8)                           :: SUMDIR(2)
       REAL(8)                           :: P(2)
       INTEGER(4)                        :: IP
       REAL(8)                           :: AREA
       REAL(8)                           :: VAL
       LOGICAL(4),SAVE                   :: TCHOICE=.TRUE.
!      *************************************************************************
       SUMVAL=0.D0
       SUMERR=0.D0
       SUMDIR(:)=0.D0
       DO IP=1,NP
         P(:)=SEGMENT%CENTER(:)+SEGMENT%WIDTH(:)*XY(:,IP)
         CALL ADAPT_INTEGRAND(P,VAL) 
         SUMVAL   =SUMVAL   +VAL*VALW(IP)
         SUMERR   =SUMERR   +VAL*ERRW(IP)
         SUMDIR(:)=SUMDIR(:)+VAL*DIRW(:,IP)
       ENDDO
!      == FACTOR [-1,1]*[-1,2]=4 IS ALREADY ABORBED IN THE WEIGHTS =============
       AREA=SEGMENT%WIDTH(1)*SEGMENT%WIDTH(2)
       SEGMENT%ERR=ABS(SUMERR)*AREA
       SEGMENT%VALUE=SUMVAL*AREA
!      == SELECT DIRECTION WITH LARGER VARIATION. (NEXT DIVISION) ==============
       IF(ABS(SUMDIR(1)).GT.ABS(SUMDIR(2))) THEN
         SEGMENT%DIVIDEAXIS=1
       ELSE IF(ABS(SUMDIR(1)).LT.ABS(SUMDIR(2))) THEN
         SEGMENT%DIVIDEAXIS=2
       ELSE  ! IF THERE IS NO CLEAR DECISION, ALTERNATE
         IF(TCHOICE) THEN
            SEGMENT%DIVIDEAXIS=1
         ELSE
            SEGMENT%DIVIDEAXIS=2
         END IF
       END IF
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE ADAPT_INTEGRAND(P,RES)
!      *************************************************************************
!      **
!      *************************************************************************
       IMPLICIT NONE
       REAL(8)  ,INTENT(IN)  :: P(2)
       REAL(8)  ,INTENT(OUT) :: RES
       REAL(8)               :: PI
       REAL(8)               :: R0(2)
       REAL(8)               :: T1(2)
       REAL(8)               :: T2(2)
       REAL(8)               :: R(2)
       REAL(8)               :: F
       REAL(8)               :: DRDP
       INTEGER(4),PARAMETER  :: TYPE=0
!      *************************************************************************
       IF(TYPE.EQ.0) THEN
         CALL LMTO_TWOCENTER_MYFUNC(P,RES)
         RETURN
!
       ELSE IF(TYPE.EQ.1) THEN
!        == AS SIMPLE CONSTANT =================================================
         R0   =(/0.D0,0.D0/)
         T1(:)=(/1.D0,0.D0/)
         T2(:)=(/0.D0,1.D0/)
         DRDP=0.25D0*ABS(T1(1)*T2(2)-T1(2)*T2(1))
         R(:)=R0(:)+T1(:)*0.5D0*(1.D0+P(1))+T2(:)*0.5D0*(1.D0+P(2))
         F=1.D0
         RES=F*DRDP
!
       ELSE IF (TYPE.EQ.2) THEN
!        == TEST 6 FROM VAN DOOREN AND DE RIDDER ===============================
!        == TARGET RESULT IS -4 ================================================
         PI=4.D0*ATAN(1.D0)
         R0   =(/0.D0,0.D0/)
         T1(:)=(/3.D0*PI,0.D0/)
         T2(:)=(/0.D0,3.D0*PI/)
         DRDP=0.25D0*ABS(T1(1)*T2(2)-T1(2)*T2(1))
         R(:)=R0(:)+T1(:)*0.5D0*(1.D0+P(1))+T2(:)*0.5D0*(1.D0+P(2))
         F=COS(R(1)+R(2))
         RES=F*DRDP
!
       ELSE IF (TYPE.EQ.3) THEN
!        == TEST 8 FROM VAN DOOREN AND DE RIDDER ===============================
!        == TARGET RESULT IS 1.047591113142868 =================================
         R0   =(/0.D0,0.D0/)
         T1(:)=(/0.D0,1.D0/)
         T2(:)=(/1.D0,0.D0/)
         R(:)=R0(:)+T1(:)*0.5D0*(1.D0+P(1))+T2(:)*0.5D0*(1.D0+P(2))
         DRDP=0.25D0*ABS(T1(1)*T2(2)-T1(2)*T2(1))
         F=605.D0*R(2)/( (1.D0+120.D0*(1.D0-R(2))) &
      &                 *((1.D0+120.D0*(1.D0-R(2)))**2+25.D0*R(1)**2*R(2)**2))
         RES=F*DRDP
!
       ELSE IF (TYPE.EQ.4) THEN
!        == TEST 9 FROM VAN DOOREN AND DE RIDDER ===============================
!        == TARGET RESULT IS 499.1249442241215 =================================
         R0   =(/0.D0,0.D0/)
         T1(:)=(/0.D0,1.D0/)
         T2(:)=(/1.D0,0.D0/)
         R(:)=R0(:)+T1(:)*0.5D0*(1.D0+P(1))+T2(:)*0.5D0*(1.D0+P(2))
         DRDP=0.25D0*ABS(T1(1)*T2(2)-T1(2)*T2(1))
         F=1.D0/((R(1)**2+1.D-4)*((R(2)+0.25D0)**2+1.D-4))
         RES=F*DRDP
!
       ELSE IF (TYPE.EQ.5) THEN
!        == TEST 10 FROM VAN DOOREN AND DE RIDDER ==============================
!        == TARGET RESULT IS 1.436563656918090   ===============================
         R0   =(/0.D0,0.D0/)
         T1(:)=(/0.D0,1.D0/)
         T2(:)=(/1.D0,0.D0/)
         R(:)=R0(:)+T1(:)*0.5D0*(1.D0+P(1))+T2(:)*0.5D0*(1.D0+P(2))
         DRDP=0.25D0*ABS(T1(1)*T2(2)-T1(2)*T2(1))
         F=EXP(ABS(R(1)+R(2)-1.D0))
         RES=F*DRDP
       ELSE
         STOP 'ERROR: UNKNOWN TYPE'
       END IF
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE ADAPT$EVALUATE(TOLERANCE,VALUE)
!      *************************************************************************
!      ** ADAPTIVE, TWO-DIMENSIONAL INTEGRATION                               **
!      ** OF THE FUNCTION ADAPT$MYFUNCTION OVER THE AREA [-1,1]X[-1,1].       **
!      ** INTEGRAL IS EVALUATED TO AN ESTIMATED ACCURACY SMALLER THAN         **
!      ** THE SPECIFIED TOLERANCE.                                            **
!      **                                                                     **
!      ** SEE: A.C. GENZ AND A.A. MALIK, J. COMPUT. APPL. MATH. 6,295 (1980)  **
!      *******************************P.E. BLOECHL, GOSLAR 2011*****************
       USE ADAPT_MODULE, ONLY : SEGMENT_TYPE
USE LMTO_TWOCENTER_MODULE, ONLY : TPR
       IMPLICIT NONE
       REAL(8),INTENT(IN)     :: TOLERANCE ! REQUIRED ACCURACY
       REAL(8),INTENT(OUT)    :: VALUE     ! INTEGRAL VALUE
       INTEGER(4),PARAMETER   :: LSTACKX=10000
       TYPE(SEGMENT_TYPE)     :: STACK(LSTACKX)
       TYPE(SEGMENT_TYPE)     :: SEGMENT1
       TYPE(SEGMENT_TYPE)     :: SEGMENT2
       INTEGER(4)             :: NSEGMENTS
       REAL(8)                :: ERROR ! CURRENT ERROR ESTIMATE
       REAL(8)                :: ERRMAX,ERRMIN
       INTEGER(4)             :: IAXIS
       INTEGER(4)             :: I
!      *************************************************************************
       CALL ADAPTINI()
!
!      =========================================================================
!      == INITIALIZE WITH FIRST SEGMENT                                       ==
!      =========================================================================
       NSEGMENTS=1
       STACK(1)%CENTER(:)=0.D0
       STACK(1)%WIDTH(:)=1.D0
       CALL ADAPT_BASICRULE(STACK(1))
       VALUE=STACK(1)%VALUE
       ERROR=STACK(1)%ERR
IF(TPR)PRINT*,'1ST ',VALUE,ERROR
!
!      =========================================================================
!      == BIG LOOP                                                            ==
!      =========================================================================
       DO 
!
!        =======================================================================
!        == DIVIDE SEGMENT WITH LARGEST ERROR (TOP OF THE STACK)              ==
!        =======================================================================
         IAXIS=STACK(1)%DIVIDEAXIS
!
         SEGMENT1%CENTER=STACK(1)%CENTER
         SEGMENT1%WIDTH=STACK(1)%WIDTH
         SEGMENT1%WIDTH(IAXIS)=SEGMENT1%WIDTH(IAXIS)*0.5D0
         SEGMENT1%CENTER(IAXIS)=SEGMENT1%CENTER(IAXIS)-SEGMENT1%WIDTH(IAXIS)
         SEGMENT1%DIVIDEAXIS=0
         SEGMENT1%VALUE=0.D0
         SEGMENT1%ERR=0.D0
!
         SEGMENT2%CENTER=STACK(1)%CENTER
         SEGMENT2%WIDTH=STACK(1)%WIDTH
         SEGMENT2%WIDTH(IAXIS)=SEGMENT2%WIDTH(IAXIS)*0.5D0
         SEGMENT2%CENTER(IAXIS)=SEGMENT2%CENTER(IAXIS)+SEGMENT2%WIDTH(IAXIS)
         SEGMENT2%DIVIDEAXIS=0
         SEGMENT2%VALUE=0.D0
         SEGMENT2%ERR=0.D0
!
!        =======================================================================
!        == INTEGRATE OVER THE TWO SEGMENTS                                   ==
!        =======================================================================
         CALL ADAPT_BASICRULE(SEGMENT1)
         CALL ADAPT_BASICRULE(SEGMENT2)
!
!        =======================================================================
!        == UPDATE TOTALS                                                     ==
!        =======================================================================
         VALUE=VALUE - STACK(1)%VALUE + SEGMENT1%VALUE + SEGMENT2%VALUE
         ERROR=ERROR - STACK(1)%ERR   + SEGMENT1%ERR   + SEGMENT2%ERR
IF(TPR.AND.MODULO(NSEGMENTS,1000).EQ.0)PRINT*,'NEXT ',NSEGMENTS,VALUE,ERROR,STACK(1)%CENTER,STACK(1)%WIDTH
         IF(ERROR.LT.TOLERANCE) RETURN
!
!        =======================================================================
!        == ORDER THE TWO SEGMENTS SO THAT SEGMENT1 HAS THE LARGER ERROR      ==
!        =======================================================================
         IF(SEGMENT2%ERR.GT.SEGMENT1%ERR) THEN
           STACK(1)=SEGMENT1
           SEGMENT1=SEGMENT2
           SEGMENT2=STACK(1)
         END IF
!
!        =======================================================================
!        == INSERT SEGMENT1 SEARCHING FROM TOP OF THE STACK                   ==
!        =======================================================================
         ERRMAX=SEGMENT1%ERR
         I=1
         DO WHILE(ERRMAX.LT.STACK(I+1)%ERR) 
           STACK(I)=STACK(I+1)
           I=I+1
           IF(I.GE.NSEGMENTS) EXIT
         ENDDO
         STACK(I)=SEGMENT1
!
!        =======================================================================
!        == INSERT SEGMENT2 SEARCHING FROM BOTTOM OF THE STACK                ==
!        =======================================================================
         ERRMIN=SEGMENT2%ERR
         NSEGMENTS=NSEGMENTS+1
         IF(NSEGMENTS.GT.LSTACKX) THEN
           CALL ERROR$MSG('NUMBER OF SEGMENTS EXCEEDS MAX')
           CALL ERROR$I4VAL('NSEGMENTS',NSEGMENTS)
           CALL ERROR$STOP('ADAPT$EVALUATE')
         END IF
         I=NSEGMENTS
         DO WHILE(ERRMIN.GT.STACK(I-1)%ERR) 
           STACK(I)=STACK(I-1)
           I=I-1
           IF(I.LE.1) EXIT
         ENDDO
         STACK(I)=SEGMENT2
       ENDDO
       STOP 'LOOP NOT CONVERGED'
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PLOTTAILED()
!     **************************************************************************
!     **  PLOTS THE LOCAL ORBITALS REPRESENTED BY TAILED ORBITALS,            **
!     **  THAT IS USING THE ONSITE STRUCTURE CONSTANTS AND EXTRAPOLATING      **
!     **  TAILS.                                                              **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY: ISPECIES,LNX,LOX,POTPAR
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: LMNXT
      INTEGER(4)             :: LMNXS
      INTEGER(4)             :: LX
      INTEGER(4)             :: LMX
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)   ,ALLOCATABLE :: C(:)
      REAL(8)   ,ALLOCATABLE :: CT(:)
      REAL(8)   ,ALLOCATABLE :: F(:,:)
      REAL(8)                :: SVAR
      INTEGER(4)             :: ISP,IAT,LM
      INTEGER(4)             :: LN,L,LMN,IM
      INTEGER(4)             :: LNT,LT,LMNT,IMT
      INTEGER(4)             :: LNXT
      CHARACTER(5)           :: CHIAT,CHORB
      CHARACTER(128)         :: STRING
!     **************************************************************************

!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        GID=POTPAR(ISP)%TAILED%GID
        CALL RADIAL$GETI4(GID,'NR',NR)
        LX=MAXVAL(LOX(:LNX(ISP),ISP))
        LMX=(LX+1)**2
        LMNX=SUM(2*LOX(:LNX(ISP),ISP)+1)
        LNXT=POTPAR(ISP)%TAILED%LNX
        LMNXT=POTPAR(ISP)%TAILED%LMNX
        ALLOCATE(C(LMNX))
        ALLOCATE(CT(LMNXT))
        ALLOCATE(F(NR,LMX))
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          DO IM=1,2*L+1
            LMN=LMN+1
            C(:)=0.D0
            C(LMN)=1.D0
            CALL LMTO_BLOWUPPSINL(IAT,LMNX,C,LMNXT,CT)
WRITE(*,FMT='("ORBITAL ",I3," WITH L=",I3," AN M=",I3)')LMN,L,IM-L-1
WRITE(*,FMT='(25F10.3)')CT
            F(:,:)=0.D0
            LMNT=0
            DO LNT=1,LNXT
              LT=POTPAR(ISP)%TAILED%LOX(LNT)
              LM=LT**2
              DO IMT=1,2*LT+1
                LMNT=LMNT+1
                LM=LM+1
                IF(ABS(CT(LMNT)).LT.1.D-8) CYCLE
                F(:,LM)=F(:,LM)+POTPAR(ISP)%TAILED%AEF(:,LNT)*CT(LMNT)
              ENDDO
            ENDDO
            WRITE(CHIAT,FMT='(I5)')IAT
            WRITE(CHORB,FMT='(I5)')LMN
            STRING='CHI_FORATOM'//TRIM(ADJUSTL(CHIAT))//'_'//TRIM(ADJUSTL(CHORB))//'.DAT'
            CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LMX,F)
          ENDDO  ! END LOOP OVER IM OF ORBITALS
        ENDDO    ! END LOOP OVER LN OF ORBITALS
        DEALLOCATE(C)
        DEALLOCATE(CT)
        DEALLOCATE(F)
      ENDDO  !END LOOP OVER ATOMS
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PLOTTAILED_OLD()
!     **************************************************************************
!     **  PLOTS THE LOCAL ORBITALS REPRESENTED BY TAILED ORBITALS,            **
!     **  THAT IS USING THE ONSITE STRUCTURE CONSTANTS AND EXTRAPOLATING      **
!     **  TAILS.                                                              **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY: K2,ISPECIES,LNX,LOX,POTPAR,SBAR,SBARLI1
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: LMNXT
      INTEGER(4)             :: LMNXS
      INTEGER(4)             :: LMX
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      INTEGER(4)             :: LX
      INTEGER(4)             :: NNS
      REAL(8)   ,ALLOCATABLE :: SBARLOC(:,:)
      REAL(8)   ,ALLOCATABLE :: F(:,:)
      INTEGER(4),ALLOCATABLE :: LMARR(:)
      REAL(8)                :: SVAR
      INTEGER(4)             :: IAT,ISP,LN,L,NN,LMN,IM,I1,IORB,LM,I0
      INTEGER(4)             :: LNDOT,LMNDOT
      CHARACTER(5)           :: CHIAT,CHORB
      CHARACTER(128)         :: STRING
!     **************************************************************************

!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        GID=POTPAR(ISP)%TAILED%GID
        CALL RADIAL$GETI4(GID,'NR',NR)
        LX=MAXVAL(LOX(:LNX(ISP),ISP))
        LMX=(LX+1)**2
!
!       ========================================================================
!       ==  DETERMINE SIZE OF STRUCTURE CONSTANT ARRAY                        ==
!       ========================================================================
        LMNXS=0 ! #(SCATTERING STATES)
        LMNX=0  ! #(VALENCE STATES)
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          LMNX=LMNX+2*L+1
          IF(POTPAR(ISP)%LNSCATT(LN).EQ.LN) LMNXS=LMNXS+2*L+1
        ENDDO
        LMNXT=LMNX+LMNXS ! #(VALENCE + SCATTERING STATES)
        ALLOCATE(SBARLOC(LMNX,LMNXS))           !C
        ALLOCATE(LMARR(LMNXT))
!    
!       ========================================================================
!       ==  COLLECT LOCAL STRUCTURE CONSTANTS                                 ==
!       ========================================================================
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          DO IM=1,2*L+1
            LMN=LMN+1
            LMARR(LMN)=L**2+IM           
          ENDDO
        ENDDO
        DO LN=1,LNX(ISP)
          IF(POTPAR(ISP)%LNSCATT(LN).NE.LN) CYCLE
          L=LOX(LN,ISP)
          DO IM=1,2*L+1
            LMN=LMN+1
            LMARR(LMN)=L**2+IM           
          ENDDO
        ENDDO
!    
!       ========================================================================
!       ==  COLLECT LOCAL STRUCTURE CONSTANTS                                 ==
!       ========================================================================
        NNS=SIZE(SBAR)
        DO NN=1,NNS
          IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
          IF(SBAR(NN)%IAT2.NE.IAT) CYCLE
          IF(SUM(SBAR(NN)%IT(:)**2).NE.0) CYCLE
          IF(SBAR(NN)%N1.NE.LMNXS.OR.SBAR(NN)%N2.NE.LMNXS) THEN
            CALL ERROR$MSG('INCONSISTENT ARRAY SIZES N1,N2')
            CALL ERROR$I4VAL('N1',LMNXS)
            CALL ERROR$I4VAL('SBAR%N1',SBAR(NN)%N1)
            CALL ERROR$I4VAL('SBAR%N2',SBAR(NN)%N2)
            CALL ERROR$STOP('LMTO_BLOWUPDENMATNL')
          END IF
          LMN=0
          DO LN=1,LNX(ISP)
            L=LOX(LN,ISP)
            I1=SBARLI1(L+1,ISP)-1
            DO IM=1,2*L+1 
              SBARLOC(LMN+IM,:)=SBAR(NN)%MAT(I1+IM,:)  !C
            ENDDO
            LMN=LMN+2*L+1
          ENDDO
          EXIT
        ENDDO
!
!       ========================================================================
!       ==  CALCULATE ORBITALS                                                ==
!       ========================================================================
        ALLOCATE(F(NR,LMX))
PRINT*,'LMARR ',LMARR
        DO IORB=1,LMNX
WRITE(*,FMT='("SBARLOC",10F10.5)')SBARLOC(IORB,:)  !C
          F(:,:)=0.D0
          LMN=0
          DO LN=1,LNX(ISP)
            L=LOX(LN,ISP)
            DO IM=1,2*L+1
              LMN=LMN+1
              LM=LMARR(LMN)
!             == ADD HEAD FUNCTION
              IF(LMN.EQ.IORB) THEN
                F(:,LMARR(LMN))=F(:,LMARR(LMN))+POTPAR(ISP)%TAILED%AEF(:,LN)
              END IF
!             == ADD TAIL FUNCTION FUNCTION
              IF(POTPAR(ISP)%LNSCATT(LN).NE.LN) CYCLE
              LNDOT=POTPAR(ISP)%TAILED%LNDOT(LN)
              LMNDOT=POTPAR(ISP)%TAILED%LMNDOT(LMN)
              F(:,LMARR(LMN))=F(:,LMARR(LMN)) &
      &               -POTPAR(ISP)%TAILED%AEF(:,LNDOT)*SBARLOC(IORB,LMNDOT-LMNX) !C
            ENDDO
          ENDDO
          WRITE(CHIAT,FMT='(I5)')IAT
          WRITE(CHORB,FMT='(I5)')IORB
          STRING='CHI_FORATOM'//TRIM(ADJUSTL(CHIAT))//'_'//TRIM(ADJUSTL(CHORB))//'.DAT'
          CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LMX,F)
          F(:,:)=0.D0
          LMN=0
          DO LN=1,LNX(ISP)
            L=LOX(LN,ISP)
            DO IM=1,2*L+1
              LMN=LMN+1
              LM=LMARR(LMN)
!             == ADD HEAD FUNCTION
              IF(LMN.EQ.IORB) THEN
                F(:,LMARR(LMN))=F(:,LMARR(LMN))+POTPAR(ISP)%TAILED%AEF(:,LN)
              END IF
            ENDDO
          ENDDO
          WRITE(CHIAT,FMT='(I5)')IAT
          WRITE(CHORB,FMT='(I5)')IORB
          STRING='XCHI_FORATOM'//TRIM(ADJUSTL(CHIAT))//'_'//TRIM(ADJUSTL(CHORB))//'.DAT'
          CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LMX,F)
        ENDDO
        DEALLOCATE(F)
!
!       ========================================================================
!       ==  REPORT SOME OTHER DATA                                            ==
!       ========================================================================
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I0=SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            LMN=LMN+1
            SVAR=POTPAR(ISP)%KTOPHIDOT(LN) &
       &      -POTPAR(ISP)%JBARTOPHIDOT(LN)*SBARLOC(LMN,I0+IM)
            SVAR=SVAR/POTPAR(ISP)%KTOPHI(LN) 
            WRITE(*,FMT='("CPHIDOT:",I5,2F10.5)')LMN,K2,SVAR
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  CLEAN UP AFTER ITERATION                                          ==
!       ========================================================================
        DEALLOCATE(LMARR)
        DEALLOCATE(SBARLOC)               
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SIMPLEENERGYTEST()
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,DENMAT,HAMIL,LNX,LOX
      IMPLICIT NONE
      INTEGER(4)            :: NNU
      INTEGER(4)            :: NND
      INTEGER(4)            :: NNH
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IND,INU,INH
      INTEGER(4)            :: NORB
      INTEGER(4)            :: NDIMD
      INTEGER(4)            :: IT(3)
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: LNX1,LMRX,LRX
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      LOGICAL(4),ALLOCATABLE:: TORB(:)
      REAL(8)   ,ALLOCATABLE:: CHI(:,:)
      REAL(8)   ,ALLOCATABLE:: AECORE(:)
      REAL(8)   ,ALLOCATABLE:: CHIPHI(:,:)
      REAL(8)   ,ALLOCATABLE:: ULITTLE(:,:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: U(:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: D(:,:,:)
      REAL(8)   ,ALLOCATABLE:: H(:,:,:)
      REAL(8)               :: EH,EX,EXTOT,EHTOT,EHDC,EXDC
      INTEGER(4)            :: NN,IAT,I,J,K,L,IS,IAT1,IAT2,ISP
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,PARAMETER  :: HFWEIGHT=0.25D0
CHARACTER(128) :: STRING
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_SIMPLEENERGYTEST')
PRINT*,'========================= ENERGYTEST ==============================='
      NAT=SIZE(ISPECIES)
      NND=SIZE(DENMAT)
      IF(.NOT.ALLOCATED(HAMIL)) THEN
        NNH=NND
        ALLOCATE(HAMIL(NND))
        DO NN=1,NND
          HAMIL(NN)%IAT1=DENMAT(NN)%IAT1
          HAMIL(NN)%IAT2=DENMAT(NN)%IAT2
          HAMIL(NN)%IT=DENMAT(NN)%IT
          N1=DENMAT(NN)%N1
          N2=DENMAT(NN)%N2
          N3=DENMAT(NN)%N3
          HAMIL(NN)%N1=N1
          HAMIL(NN)%N2=N2
          HAMIL(NN)%N3=N3
          ALLOCATE(HAMIL(NN)%MAT(N1,N2,N3))
          HAMIL(NN)%MAT=0.D0
        ENDDO
      END IF
      NNH=SIZE(HAMIL)

      EXTOT=0.D0
      EHTOT=0.D0
      DO IAT=1,NAT
!
!       == FIND LOCAL DENSITY MATRIX ===========================================
        IND=-1
        DO NN=1,NND
          IF(DENMAT(NN)%IAT1.NE.IAT) CYCLE
          IF(DENMAT(NN)%IAT2.NE.IAT) CYCLE
          IF(SUM(DENMAT(NN)%IT**2).NE.0) CYCLE
          IND=NN
          EXIT
        ENDDO
        IF(IND.LE.0) THEN
          CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED: IND<0')
          CALL ERROR$STOP('LMTO_ENERGYTEST')
        END IF
!
!       == FIND LOCAL HAMILTONIAN   ============================================
        INH=-1
        DO NN=1,NNH
          IF(HAMIL(NN)%IAT1.NE.IAT) CYCLE
          IF(HAMIL(NN)%IAT2.NE.IAT) CYCLE
          IF(SUM(HAMIL(NN)%IT**2).NE.0) CYCLE
          INH=NN
          EXIT
        ENDDO
        IF(INH.LE.0) THEN
          CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED: INH<0')
          CALL ERROR$STOP('LMTO_ENERGYTEST')
        END IF
!
        NORB=DENMAT(IND)%N1
        NDIMD=DENMAT(IND)%N3
!
!       ========================================================================
!       == CALCULATE U-TENSOR                                                 ==
!       ========================================================================
        ALLOCATE(U(NORB,NORB,NORB,NORB))
        ISP=ISPECIES(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(AECORE(NR))
        CALL SETUP$GETR8A('AECORE',NR,AECORE)
        CALL SETUP$GETI4('LMRX',LMRX)
        CALL SETUP$UNSELECT()
        LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
        LNX1=LNX(ISP)
        ALLOCATE(TORB(LNX1))
        TORB=.TRUE.
        ALLOCATE(CHI(NR,LNX1))
        ALLOCATE(CHIPHI(LNX1,LNX1))
        ALLOCATE(ULITTLE(LRX+1,LNX1,LNX1,LNX1,LNX1))
!        CALL LMTO$DOLOCORB(IAT,ISP,GID,NR,LNX1,LNX1,TORB,CHIPHI,CHI)
        CALL LMTO$DOLOCORB_2(IAT,ISP,GID,NR,LNX1,LNX1,TORB,CHIPHI,CHI)
        CALL LMTO_ULITTLE(GID,NR,LRX,LNX1,LOX(:LNX1,ISP),CHI,ULITTLE)
        CALL LMTO_UTENSOR(LRX,NORB,LNX1,LOX(:LNX1,ISP),ULITTLE,U)
        DEALLOCATE(ULITTLE)
        DEALLOCATE(CHIPHI)
        DEALLOCATE(TORB)
!
WRITE(STRING,FMT='(I5)')IAT
STRING='CHI_FORATOM'//TRIM(ADJUSTL(STRING))//'.DAT'
CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LNX1,CHI)
!
!       ========================================================================
!       == WORK OUT TOTAL ENERGY AND HAMILTONIAN                              ==
!       ========================================================================
        ALLOCATE(D(NORB,NORB,NDIMD))
        ALLOCATE(H(NORB,NORB,NDIMD))
        D=DENMAT(IND)%MAT
        EH=0.D0
        EX=0.D0
        H(:,:,:)=0.D0
        DO I=1,NORB
          DO J=1,NORB
            DO K=1,NORB
              DO L=1,NORB
!               ================================================================
!               == HARTREE TERM (NOT CONSIDERED)                              ==
!               == EH=EH+0.5D0*U(I,J,K,L)*D(K,I,1)*D(L,J,1)                   ==
!               == H(K,I,1)=H(K,I,1)+0.5D0*U(I,J,K,L)*D(L,J,1) !DE/DRHO(K,I,1)==
!               == H(L,J,1)=H(L,J,1)+0.5D0*U(I,J,K,L)*D(K,I,1) !DE/DRHO(L,J,1)==
!               ================================================================
!               ================================================================
!               == EXCHANGE ENERGY =============================================
!               ================================================================
                DO IS=1,NDIMD
                  EX=EX-0.25D0*U(I,J,K,L)*D(K,J,IS)*D(L,I,IS)
                  H(K,J,IS)=H(K,J,IS)-0.25D0*U(I,J,K,L)*D(L,I,IS) 
                  H(L,I,IS)=H(L,I,IS)-0.25D0*U(I,J,K,L)*D(K,J,IS) 
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        HAMIL(INH)%MAT=H
        EXTOT=EXTOT+EX
PRINT*,'ONSITE EXACT EXCHANGE ENERGY FOR ATOM=',IAT,EX
!
!       ========================================================================
!       == ADD CORE VALENCE EXCHANGE                                          ==
!       ========================================================================
        CALL LMTO_CVX(ISPECIES(IAT),NORB,EX,D(:,:,1),H(:,:,1))
        EXTOT=EXTOT+EX
        HAMIL(INH)%MAT(:,:,1)=HAMIL(INH)%MAT(:,:,1)+H(:,:,1)
PRINT*,'CORE VALENCE EXCHANGE ENERGY FOR ATOM=',IAT,EX
!
!       ========================================================================
!       == DOUBLE COUNTING CORRECTION (EXCHANGE ONLY)                         ==
!       ========================================================================
! THIS IS THE TIME CONSUMIN PART OF ENERGYTEST
CALL TIMING$CLOCKON('ENERGYTEST:DC')      
!!$        CALL LMTO_DOUBLECOUNTING(IAT,NDIMD,NORB,D,EX,H)
        CALL LMTO_SIMPLEDC_OLD(GID,NR,NORB,LNX1,LOX(:LNX1,ISP),CHI,LRX,AECORE &
     &                        ,D,EX,H)
        EXTOT=EXTOT-EX
        HAMIL(INH)%MAT=HAMIL(INH)%MAT-H
PRINT*,'DOUBLE COUNTING ',IAT,EX
CALL TIMING$CLOCKOFF('ENERGYTEST:DC')      
        DEALLOCATE(U)
        DEALLOCATE(H)
        DEALLOCATE(D)
        DEALLOCATE(CHI)
        DEALLOCATE(AECORE)
      ENDDO
      EH=EHTOT
      EX=EXTOT
!
!     ==========================================================================
!     == RESCALE WITH HFWEIGHT                                                ==
!     ==========================================================================
      DO NN=1,NNH
        HAMIL(NN)%MAT=HAMIL(NN)%MAT*HFWEIGHT
      ENDDO
      EXTOT=EXTOT*HFWEIGHT
      EHTOT=EHTOT*HFWEIGHT
!
!     ==========================================================================
!     == COMMUNICATE ENERGY TO ENERGYLIST                                     ==
!     ==========================================================================
      CALL ENERGYLIST$SET('LMTO INTERFACE',EXTOT)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EXTOT)
!
!     ==========================================================================
!     == PRINT FOR TEST                                                       ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,FMT='(82("="),T10," DENSITY MATRIX IN A NTBO BASIS ")')
        WRITE(*,FMT='(82("="),T10," ONLY ONSITE TERMS WRITTEN ")')
        DO NN=1,NND
          IAT1=DENMAT(NN)%IAT1
          IAT2=DENMAT(NN)%IAT2
          IT=DENMAT(NN)%IT
          IF(IAT1.NE.IAT2.OR.IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
     &                                                           IAT1,IAT2,IT
          N1=DENMAT(NN)%N1
          N2=DENMAT(NN)%N2
          DO I=1,1 !DENMAT(NN)%N3
            DO J=1,N2 
              WRITE(*,FMT='(I3,30F10.3)')I,DENMAT(NN)%MAT(J,:,I)
            ENDDO
            WRITE(*,FMT='(82("-"))')
          ENDDO
        ENDDO
        WRITE(*,FMT='("XC ENERGY ",F10.5)')EXTOT
        WRITE(*,FMT='(82("="),T10," HAMILTON MATRIX IN A NTBO BASIS ")')
        WRITE(*,FMT='(82("="),T10," ONLY ONSITE TERMS WRITTEN ")')
        DO NN=1,NND
          IAT1=HAMIL(NN)%IAT1
          IAT2=HAMIL(NN)%IAT2
          IT=HAMIL(NN)%IT
          IF(IAT1.NE.IAT2.OR.IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
     &                                                            IAT1,IAT2,IT
          N1=HAMIL(NN)%N1
          N2=HAMIL(NN)%N2
          DO I=1,1 !HAMIL(NN)%N3
            DO J=1,N2 
              WRITE(*,FMT='(I3,30F10.3)')I,HAMIL(NN)%MAT(J,:,I)
            ENDDO
            WRITE(*,FMT='(82("-"))')
          ENDDO
        ENDDO
      END IF
!
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_ENERGYTEST()
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,DENMAT,UTENSOR,HAMIL
      IMPLICIT NONE
      INTEGER(4)            :: NNU
      INTEGER(4)            :: NND
      INTEGER(4)            :: NNH
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IND,INU,INH
      INTEGER(4)            :: NORB
      INTEGER(4)            :: NDIMD
      INTEGER(4)            :: IT(3)
      INTEGER(4)            :: N1,N2,N3
      REAL(8)   ,ALLOCATABLE:: U(:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: D(:,:,:)
      REAL(8)   ,ALLOCATABLE:: H(:,:,:)
      REAL(8)               :: EH,EX,EXTOT,EHTOT,EHDC,EXDC
      INTEGER(4)            :: NN,IAT,I,J,K,L,IS,IAT1,IAT2
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,PARAMETER  :: HFWEIGHT=0.25D0
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_ENERGYTEST')
      NAT=SIZE(ISPECIES)
      NNU=SIZE(UTENSOR)
      NND=SIZE(DENMAT)
      IF(.NOT.ALLOCATED(HAMIL)) THEN
        NNH=NND
        ALLOCATE(HAMIL(NND))
        DO NN=1,NND
          HAMIL(NN)%IAT1=DENMAT(NN)%IAT1
          HAMIL(NN)%IAT2=DENMAT(NN)%IAT2
          HAMIL(NN)%IT=DENMAT(NN)%IT
          N1=DENMAT(NN)%N1
          N2=DENMAT(NN)%N2
          N3=DENMAT(NN)%N3
          HAMIL(NN)%N1=N1
          HAMIL(NN)%N2=N2
          HAMIL(NN)%N3=N3
          ALLOCATE(HAMIL(NN)%MAT(N1,N2,N3))
          HAMIL(NN)%MAT=0.D0
        ENDDO
      END IF
      NNH=SIZE(HAMIL)

      EXTOT=0.D0
      EHTOT=0.D0
      DO IAT=1,NAT
!
!       == FIND LOCAL DENSITY MATRIX ===========================================
        IND=-1
        DO NN=1,NND
          IF(DENMAT(NN)%IAT1.NE.IAT) CYCLE
          IF(DENMAT(NN)%IAT2.NE.IAT) CYCLE
          IF(SUM(DENMAT(NN)%IT**2).NE.0) CYCLE
          IND=NN
          EXIT
        ENDDO
        IF(IND.LE.0) THEN
          CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED: IND<0')
          CALL ERROR$STOP('LMTO_ENERGYTEST')
        END IF
!
!       == FIND LOCAL HAMILTONIAN   ===========================================
        INH=-1
        DO NN=1,NNH
          IF(HAMIL(NN)%IAT1.NE.IAT) CYCLE
          IF(HAMIL(NN)%IAT2.NE.IAT) CYCLE
          IF(SUM(HAMIL(NN)%IT**2).NE.0) CYCLE
          INH=NN
          EXIT
        ENDDO
        IF(INH.LE.0) THEN
          CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED: INH<0')
          CALL ERROR$STOP('LMTO_ENERGYTEST')
        END IF
!
!       == FIND LOCAL U-TENSOR =================================================
        INU=-1
        DO NN=1,NNU
          IF(UTENSOR(NN)%IAT1.NE.IAT) CYCLE
          IF(UTENSOR(NN)%IAT2.NE.IAT) CYCLE
          IF(UTENSOR(NN)%IAT3.NE.IAT) CYCLE
          IF(UTENSOR(NN)%IAT4.NE.IAT) CYCLE
          IT=UTENSOR(NN)%IT2
          IF(IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
          IT=UTENSOR(NN)%IT3
          IF(IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
          IT=UTENSOR(NN)%IT4
          IF(IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
          INU=NN
          EXIT
        ENDDO
        IF(INU.LE.0) THEN
          CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED: INU<0')
          CALL ERROR$STOP('LMTO_ENERGYTEST')
        END IF
!
!       ========================================================================
!       == WORK OUT TOTAL ENERGY AND HAMILTONIAN                              ==
!       ========================================================================
        NORB=DENMAT(IND)%N1
        NDIMD=DENMAT(IND)%N3
        ALLOCATE(U(NORB,NORB,NORB,NORB))
        ALLOCATE(D(NORB,NORB,NDIMD))
        ALLOCATE(H(NORB,NORB,NDIMD))
        U=UTENSOR(INU)%U
        D=DENMAT(IND)%MAT
        EH=0.D0
        EX=0.D0
        H(:,:,:)=0.D0
        DO I=1,NORB
          DO J=1,NORB
            DO K=1,NORB
              DO L=1,NORB
!               ================================================================
!               == HARTREE TERM (NOT CONSIDERED)                              ==
!               == EH=EH+0.5D0*U(I,J,K,L)*D(K,I,1)*D(L,J,1)                   ==
!               == H(K,I,1)=H(K,I,1)+0.5D0*U(I,J,K,L)*D(L,J,1) !DE/DRHO(K,I,1)==
!               == H(L,J,1)=H(L,J,1)+0.5D0*U(I,J,K,L)*D(K,I,1) !DE/DRHO(L,J,1)==
!               ================================================================
!               ================================================================
!               == EXCHANGE ENERGY =============================================
!               ================================================================
                DO IS=1,NDIMD
                  EX=EX-0.25D0*U(I,J,K,L)*D(K,J,IS)*D(L,I,IS)
                  H(K,J,IS)=H(K,J,IS)-0.25D0*U(I,J,K,L)*D(L,I,IS) 
                  H(L,I,IS)=H(L,I,IS)-0.25D0*U(I,J,K,L)*D(K,J,IS) 
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        HAMIL(INH)%MAT=H
        EXTOT=EXTOT+EX
!
!       ========================================================================
!       == ADD CORE VALENCE EXCHANGE                                          ==
!       ========================================================================
        CALL LMTO_CVX(ISPECIES(IAT),NORB,EX,D(:,:,1),H(:,:,1))
        EXTOT=EXTOT+EX
        HAMIL(INH)%MAT(:,:,1)=HAMIL(INH)%MAT(:,:,1)+H(:,:,1)
!
!       ========================================================================
!       == DOUBLE COUNTING CORRECTION (EXCHANGE ONLY)                         ==
!       ========================================================================
! THIS IS THE TIME CONSUMIN PART OF ENERGYTEST
CALL TIMING$CLOCKON('ENERGYTEST:DC')      
        CALL LMTO_DOUBLECOUNTING(IAT,NDIMD,NORB,D,EX,H)
        EXTOT=EXTOT-EX
        HAMIL(INH)%MAT=HAMIL(INH)%MAT-H
CALL TIMING$CLOCKOFF('ENERGYTEST:DC')      
        DEALLOCATE(U)
        DEALLOCATE(H)
        DEALLOCATE(D)
      ENDDO
      EH=EHTOT
      EX=EXTOT
!
!     ==========================================================================
!     == RESCALE WITH HFWEIGHT                                                ==
!     ==========================================================================
      DO NN=1,NNH
        HAMIL(NN)%MAT=HAMIL(NN)%MAT*HFWEIGHT
      ENDDO
      EXTOT=EXTOT*HFWEIGHT
      EHTOT=EHTOT*HFWEIGHT
!
!     ==========================================================================
!     == COMMUNICATE ENERGY TO ENERGYLIST                                     ==
!     ==========================================================================
      CALL ENERGYLIST$SET('LMTO INTERFACE',EXTOT)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EXTOT)
!
!     ==========================================================================
!     == PRINT FOR TEST                                                       ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,FMT='(82("="),T10," DENSITY MATRIX IN A NTBO BASIS ")')
        WRITE(*,FMT='(82("="),T10," ONLY ONSITE TERMS WRITTEN ")')
        DO NN=1,NND
          IAT1=DENMAT(NN)%IAT1
          IAT2=DENMAT(NN)%IAT2
          IT=DENMAT(NN)%IT
          IF(IAT1.NE.IAT2.OR.IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
     &                                                           IAT1,IAT2,IT
          N1=DENMAT(NN)%N1
          N2=DENMAT(NN)%N2
          DO I=1,1 !DENMAT(NN)%N3
            DO J=1,N2 
              WRITE(*,FMT='(I3,30F10.3)')I,DENMAT(NN)%MAT(J,:,I)
            ENDDO
            WRITE(*,FMT='(82("-"))')
          ENDDO
        ENDDO
        WRITE(*,FMT='("XC ENERGY ",F10.5)')EXTOT
        WRITE(*,FMT='(82("="),T10," HAMILTON MATRIX IN A NTBO BASIS ")')
        WRITE(*,FMT='(82("="),T10," ONLY ONSITE TERMS WRITTEN ")')
        DO NN=1,NND
          IAT1=HAMIL(NN)%IAT1
          IAT2=HAMIL(NN)%IAT2
          IT=HAMIL(NN)%IT
          IF(IAT1.NE.IAT2.OR.IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
     &                                                            IAT1,IAT2,IT
          N1=HAMIL(NN)%N1
          N2=HAMIL(NN)%N2
          DO I=1,1 !HAMIL(NN)%N3
            DO J=1,N2 
              WRITE(*,FMT='(I3,30F10.3)')I,HAMIL(NN)%MAT(J,:,I)
            ENDDO
            WRITE(*,FMT='(82("-"))')
          ENDDO
        ENDDO
      END IF
!
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DOUBLECOUNTING(IAT,NDIMD,NORB,D,EX,H)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : ISPECIES,LMORB
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: NORB
      INTEGER(4),INTENT(IN) :: NDIMD
      REAL(8)   ,INTENT(IN) :: D(NORB,NORB,NDIMD)
      REAL(8)   ,INTENT(OUT):: EX
      REAL(8)   ,INTENT(OUT):: H(NORB,NORB,NDIMD)
      INTEGER(4)            :: NND
      INTEGER(4)            :: NNH
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IND,INH
      INTEGER(4)            :: IT(3)
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: LMX
      INTEGER(4)            :: LMRX
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)   ,ALLOCATABLE:: AECORE(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:)
      REAL(8)   ,ALLOCATABLE:: AUX1(:)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)  !(NR,LMR,IDIMD) DENSITY
      REAL(8)   ,ALLOCATABLE:: RHOWC(:,:,:)  !(NR,LMR,IDIMD) DENSITY W CORE
      REAL(8)   ,ALLOCATABLE:: POT(:,:,:)  !(NR,LMR,IDIMD) POTENTIAL
      REAL(8)               :: EH,ETOTC
      REAL(8)               :: SVAR
      REAL(8)               :: CG ! GAUNT COEFFICIENT
      INTEGER(4)            :: I,J,ISP,L
      INTEGER(4)            :: IORB1,IORB2,LM1,LM2,LM3,IDIMD
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      LOGICAL(4),PARAMETER  :: TCV=.TRUE. ! CORE-VALENCE CONTRIBUTION
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_DOUBLECOUNTING')
      CALL DFT$SETL4('XCONLY',.TRUE.)
!
!     ==========================================================================
!     == COLLECT FURTHER INFORMATION                                          ==
!     ==========================================================================
      GID=LMORB(IAT)%GID
      NR=LMORB(IAT)%NR
      LMX=LMORB(IAT)%LMX
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
!
      ISP=ISPECIES(IAT)
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('LMRX',LMRX)
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
!
!     ==========================================================================
!     == CALCULATE DENSITY                                                    ==
!     ==========================================================================
!     == CALCULATE LOCAL ELECTRON DENSITY ======================================
      ALLOCATE(RHO(NR,LMRX,NDIMD))
      ALLOCATE(POT(NR,LMRX,NDIMD))
      RHO(:,:,:)=0.D0
      DO IORB1=1,NORB
        DO IORB2=1,NORB
          DO LM1=1,LMX
            DO LM2=1,LMX
              AUX(:)=LMORB(IAT)%F(:,LM1,IORB1)*LMORB(IAT)%F(:,LM2,IORB2)
              DO IDIMD=1,NDIMD
                SVAR=D(IORB1,IORB2,IDIMD)
                DO LM3=1,LMRX
                  CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
                  RHO(:,LM3,IDIMD)=RHO(:,LM3,IDIMD)+AUX(:)*CG*SVAR
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == ADD CORE DENSITY IF NECESSARY                                        ==
!     ==========================================================================
      ALLOCATE(RHOWC(NR,LMRX,NDIMD))  !WITH CORE
      RHOWC(:,:,:)=RHO(:,:,:)
      IF(TCV) THEN
        CALL SETUP$GETR8A('AECORE',NR,AUX)
        CALL AUGMENTATION_XC(GID,NR,1,1,AUX,ETOTC,POT)
        RHOWC(:,1,1)=RHO(:,1,1)+AUX(:)
      ELSE 
        ETOTC=0.D0
      END IF
!
!     ========================================================================
!     == CALCULATE TOTAL ENERGY                                             ==
!     ========================================================================
!
!== CALCULATE HARTREE ENERGY FOR TEST ==========================================
AUX1=0.D0
DO LM1=1,LMRX
  L=INT(SQRT(REAL(LM1-1,KIND=8))+1.D-5)
  CALL RADIAL$POISSON(GID,NR,L,RHO(:,LM1,1),AUX)
  AUX1(:)=AUX1(:)+0.5D0*AUX(:)*RHO(:,LM1,1)
ENDDO
CALL RADIAL$INTEGRAL(GID,NR,AUX1*R**2,EH)
!
!     == EXCHANGE CORRELATION ENERGY AND POTENTIAL =============================
      POT(:,:,:)=0.D0
      CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOWC,EX,POT)
      EX=EX-ETOTC  ! SUBTRACT CORE EXCHANGE
!
!     ==========================================================================
!     == HAMILTON MATRIX ELEMENTS                                             ==
!     ==========================================================================
      H(:,:,:)=0.D0
      DO IORB1=1,NORB
        DO IORB2=1,NORB
          DO IDIMD=1,NDIMD
            AUX1(:)=0.D0
            DO LM1=1,LMX
              DO LM2=1,LMX
                AUX(:)=LMORB(IAT)%F(:,LM1,IORB1)*LMORB(IAT)%F(:,LM2,IORB2)
                DO LM3=1,LMRX
                  CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
                  AUX1(:)=AUX1(:)+POT(:,LM3,IDIMD)*AUX(:)*CG
                ENDDO
              ENDDO
            ENDDO
            AUX1(:)=AUX1(:)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX1,H(IORB1,IORB2,IDIMD))
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(RHO)
      DEALLOCATE(RHOWC)
      DEALLOCATE(POT)
      DEALLOCATE(R)
      DEALLOCATE(AUX)
      DEALLOCATE(AUX1)
      CALL DFT$SETL4('XCONLY',.FALSE.)
      CALL SETUP$UNSELECT()
!
!     ==========================================================================
!     == PRINT FOR TEST                                                       ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,FMT='(82("="),T10," NON-LOCAL HAMILTON MATRIX IN A NTBO BASIS ")')
        PRINT*,'LMTO_DOUBLECOUNTING TOTAL ENERGY=',EX+EH,' HARTREE=',EH &
     &                         ,' EXC=',EX+EH,' EX(CORE)=',ETOTC
        DO IORB1=1,NORB 
          WRITE(*,FMT='(I3,30F10.3)')IORB1,H(IORB1,:,1)
        ENDDO
        WRITE(*,FMT='(82("-"))')
      END IF

                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TESTDENMAT_1CDENMAT(LMNXX_,NDIMD_,NAT,DENMAT_)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : DENMAT,SBAR,ISPECIES,POTPAR,LNX,LOX
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMNXX_
      INTEGER(4),INTENT(IN) :: NDIMD_
      INTEGER(4),INTENT(IN) :: NAT
      COMPLEX(8),INTENT(IN) :: DENMAT_(LMNXX_,LMNXX_,NDIMD_,NAT)
      INTEGER(4)            :: LMNXX
      INTEGER(4)            :: NND
      INTEGER(4)            :: NNS
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: IND,INS
      INTEGER(4)            :: NN,I,IDIM,IAT,ISP,LN,LMN,IM,I1,I2
      REAL(8)   ,ALLOCATABLE:: KTOPHI(:)
      REAL(8)   ,ALLOCATABLE:: KTOPHIDOT(:)
      REAL(8)   ,ALLOCATABLE:: JBARTOPHIDOT(:)
      REAL(8)   ,ALLOCATABLE:: MAT11(:,:,:)
      REAL(8)   ,ALLOCATABLE:: MAT12(:,:,:)
      REAL(8)   ,ALLOCATABLE:: MAT22(:,:,:)
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
!     **************************************************************************
      IF(.NOT.TPR) RETURN
                                  CALL TRACE$PUSH('LMTO_TESTDENMAT_1CDENMAT')
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      IF(NAT.NE.SIZE(ISPECIES)) THEN
        CALL ERROR$MSG('INCONSISTENT DATA: NAT DIFFERS FROM NAT_')
        CALL ERROR$STOP('LMTO_TESTDENMAT_1CDENMAT')
      END IF
!
!     ==========================================================================
!     == PREPARE WRONSKIANS                                                   ==
!     ==========================================================================
!
!     ==========================================================================
!     == CALCULATE ONE-CENTER DENSITY MATRIX                                  ==
!     ==========================================================================
      NND=SIZE(DENMAT)
      NNS=SIZE(SBAR)
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
!       == EXPAND POTENTIAL PARAMETER ARRAYS ===================================
        LMNXX=SUM(2*LOX(:,ISP)+1)
        ALLOCATE(KTOPHI(LMNXX))
        ALLOCATE(KTOPHIDOT(LMNXX))
        ALLOCATE(JBARTOPHIDOT(LMNXX))
        ALLOCATE(MAT11(LMNXX,LMNXX,4))
        ALLOCATE(MAT12(LMNXX,LMNXX,4))
        ALLOCATE(MAT22(LMNXX,LMNXX,4))
        LMN=0
        DO LN=1,LNX(ISP)
          DO IM=1,2*LOX(LN,ISP)+1
            LMN=LMN+1
            KTOPHI(LMN)      =POTPAR(ISP)%KTOPHI(LN)
            KTOPHIDOT(LMN)   =POTPAR(ISP)%KTOPHIDOT(LN)
            JBARTOPHIDOT(LMN)=POTPAR(ISP)%JBARTOPHIDOT(LN)
          ENDDO
        ENDDO
!       == WRITE ORIGINAL DENSITY MATRIX =======================================
        DO IDIM=1,NDIMD_
          WRITE(*,FMT='(82("="),T30," IAT=",I3,"  AND IDIM=",I1,"  ")') &
    &                                         IAT,IDIM
          DO I=1,LMNXX_
            WRITE(*,FMT='(100F10.5)')REAL(DENMAT_(I,:,IDIM,IAT))
          ENDDO
        ENDDO
!
!       == CALCULATE FROM NTBO DENSITY MATRIX ==================================
        DO NN=1,NNS
          IF(DENMAT(NN)%IAT1.EQ.IAT) CYCLE
          IF(DENMAT(NN)%IAT2.EQ.IAT) CYCLE
          IF(SUM(DENMAT(NN)%IT**2).NE.0) CYCLE
          N1=DENMAT(NN)%N1  
          N2=DENMAT(NN)%N2  
          N3=DENMAT(NN)%N3  
          MAT11(:,:,:)=DENMAT(NN)%MAT(:,:,:)
          MAT12(:,:,:)=MAT11(:,:,:)
          MAT22(:,:,:)=MAT11(:,:,:)
          DO IDIM=1,N3
            DO I2=1,N2
              MAT11(:,I2,IDIM)=KTOPHI(:)*MAT11(:,I2,IDIM)
              MAT12(:,I2,IDIM)=KTOPHI(:)*MAT12(:,I2,IDIM)
              MAT12(:,22,IDIM)=KTOPHIDOT(:)*MAT22(:,I2,IDIM)
            ENDDO
            DO I1=1,N1
              MAT11(I1,:,IDIM)=MAT11(I1,:,IDIM)*KTOPHI(:)
              MAT12(I1,:,IDIM)=MAT12(I1,:,IDIM)*KTOPHIDOT(:)
              MAT22(I1,:,IDIM)=MAT22(I1,:,IDIM)*KTOPHIDOT(:)
            ENDDO
          ENDDO
        ENDDO
!
        DO IND=1,NND
          IF(DENMAT(IND)%IAT1.EQ.IAT) CYCLE
          DO INS=1,NNS
            IF(SBAR(INS)%IAT2.EQ.IAT) CYCLE
            IF(SUM(DENMAT(IND)%IT+SBAR(INS)%IT)**2.NE.0) CYCLE
          ENDDO
        ENDDO


!       == WRITE RESULT ========================================================
        DO IDIM=1,N3
          WRITE(*,FMT='(82("-"),T30," IAT=",I3,"  AND IDIM=",I1,"  ")') &
    &                                       IAT,IDIM
          DO I=1,N1
            WRITE(*,FMT='(100F10.5)')MAT11(I,:,IDIM)
          ENDDO
        ENDDO
!
        DEALLOCATE(KTOPHI)
        DEALLOCATE(KTOPHIDOT)
        DEALLOCATE(JBARTOPHIDOT)
        DEALLOCATE(MAT11)
        DEALLOCATE(MAT12)
        DEALLOCATE(MAT22)
      ENDDO
STOP 'FORCED STOP IN LMTO_TESTDENMAT_1CENTER'
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TESTDENMAT()
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : OVERLAP,DENMAT
      IMPLICIT NONE
      INTEGER(4)            :: NND
      INTEGER(4)            :: NNO
      INTEGER(4)            :: IAT1,IAT2,IT(3)
      INTEGER(4)            :: N1,N2
      INTEGER(4)            :: NN,MM,I,J
      REAL(8)               :: SVAR1,SVAR2
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
!     **************************************************************************
      IF(.NOT.TPR) RETURN
                                           CALL TRACE$PUSH('LMTO_TESTDENMAT')
      NND=SIZE(DENMAT)
      NNO=SIZE(OVERLAP)
      SVAR1=0.D0
      DO NN=1,NND
        IAT1=DENMAT(NN)%IAT1
        IAT2=DENMAT(NN)%IAT2
        IT=DENMAT(NN)%IT
        WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
     &                                                             IAT1,IAT2,IT
        DO MM=1,NNO
          IF(OVERLAP(MM)%IAT1.NE.IAT1) CYCLE
          IF(OVERLAP(MM)%IAT2.NE.IAT2) CYCLE
          IF(OVERLAP(MM)%IT(1).NE.IT(1)) CYCLE
          IF(OVERLAP(MM)%IT(2).NE.IT(2)) CYCLE
          IF(OVERLAP(MM)%IT(3).NE.IT(3)) CYCLE
          N1=DENMAT(NN)%N1
          N2=DENMAT(NN)%N2
PRINT*,'N1,N2 ',N1,N2
          SVAR2=0.D0
          DO J=1,N1 
            WRITE(*,FMT='("D ",200F10.3)')DENMAT(NN)%MAT(J,:,1)
            WRITE(*,FMT='("O ",200F10.3)')OVERLAP(MM)%MAT(:,J)
            DO I=1,N2
              SVAR2=SVAR2+DENMAT(NN)%MAT(I,J,1)*OVERLAP(MM)%MAT(J,I)
            ENDDO
          ENDDO
          SVAR1=SVAR1+SVAR2
          PRINT*,'#PARTICLES FROM THIS PAIR ',SVAR2,MM,NN,SVAR1
        ENDDO
      ENDDO
      PRINT*,'TOTAL #PARTICLES ',SVAR1
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$DOLOCORB_2(IAT,ISP,GID,NR,LNXCHI,LNXPHI,TORB,CHIPHI,CHI)
!     **************************************************************************
!     ** NEW VERSION!!!!!
!     **  CONSTRUCTS ONSITE-MAPPING FROM PARTIAL WAVES TO LOCAL ORBITALS      **
!     **                                                                      **
!     **  TORB SELECTS LOCAL ORBITALS FROM PARTIAL WAVES                      **
!     **                                                                      **
!     **  TAILS ARE DEFINED BY SCATTERING WAVE FUNCTION                       **
!     **      FOR THE VALENCE STATE (I.E.ISCATT=0)                            **
!     **                                                                      **
!     **  OFF-SITE TERMS ARE REPLACED BY EXPONENTIAL TAIL MATCHED DIFFERENTIABLY
!     **                                                                      **
!     **  ON SITE STRUCTURE CONSTANTS ARE SPHERICALLY AVERAGED                **
!     **                                                                      **
!     **   |PSI>= SUM_{I,J} |CHI_I>*CHIPHI(I,J)*<PTILDE_J|\PSITILDE>          **
!     **       WHERE J INCLUDES ONLY ONSITE TERMS                             **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2011 ************
      USE LMTO_MODULE, ONLY : SBAR,POTPAR,SBARLI1,K2
      USE STRINGS_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT     ! ATOM INDEX
      INTEGER(4),INTENT(IN)  :: GID     ! GRID ID
      INTEGER(4),INTENT(IN)  :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN)  :: LNXPHI  ! #(PARTIAL WAVES W/O M-MULTIPLICITY)
      INTEGER(4),INTENT(IN)  :: LNXCHI  ! #(PARTIAL WAVES W/O M-MULTIPLICITY)
      LOGICAL(4),INTENT(IN)  :: TORB(LNXPHI)          ! SELECTS LOCAL ORBITALS
      REAL(8)   ,INTENT(OUT) :: CHIPHI(LNXCHI,LNXPHI) !<PI_I|=\SUM_J CHIPHI(I,J)<P_J|
      REAL(8)   ,INTENT(OUT) :: CHI(NR,LNXCHI)
      INTEGER(4)             :: ISP         ! ATOM TYPE INDEX
      INTEGER(4)             :: LNX         ! #(PARTIAL WAVES)
      INTEGER(4)             :: NNS         ! 
      INTEGER(4)             :: LOX(LNXPHI) ! ANGULAR MOMENTA
      REAL(8)                :: PRO(NR,LNXPHI)    ! PROJECTOR FUNCTIONS
      REAL(8)   ,ALLOCATABLE :: AECHI(:,:)  ! ALL-ELECTRON LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: PSCHI(:,:)  ! PSEUDO LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: NLCHI(:,:)  ! NODELESS LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: SBARAV(:)
      REAL(8)   ,ALLOCATABLE :: AMAT(:,:)   !(LNXCHI1,LNXPHI) 
      REAL(8)   ,ALLOCATABLE :: BMAT(:,:)   !(LNXCHI1,LNXCHI1) 
      REAL(8)   ,ALLOCATABLE :: XMAT(:,:)   !(LNXPHI,LNXCHI) 
      REAL(8)                :: RAD              ! COVALENT RADIUS
      REAL(8)                :: AEZ               ! ATOMIC NUMBER
      REAL(8)                :: AUX(NR)
      REAL(8)                :: R(NR)
      REAL(8)                :: SVAR,SVAR1,SVAR2,VAL,DER
      REAL(8)                :: KVAL,KDER,JVAL,JDER
      REAL(8)                :: QBAR
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: LN,LN1,L,I,J,IIB,LM,LNCHI,IR,IM
      INTEGER(4)             :: LX
      INTEGER(4)             :: IORB
      INTEGER(4)             :: IRAD  ! GRID INDEX JUST BEYOND RAD
      CHARACTER(64)          :: STRING
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO$DOLOCORB_2')
!
!     ==========================================================================
!     == CHECK CONSISTENCY OF INPUT                                           ==
!     ==========================================================================
      LNCHI=0
      DO LN=1,LNXPHI
        IF(TORB(LN))LNCHI=LNCHI+1
      ENDDO
      IF(LNCHI.NE.LNXCHI) THEN
        CALL ERROR$MSG('LOCAL-ORBITAL SELECTION TORB INCONSISTENT WITH LNXCHI')
        CALL ERROR$L4VAL('TORB',TORB)
        CALL ERROR$I4VAL('LNXCHI',LNXCHI)
        CALL ERROR$STOP('LMTO$DOLOCORB_2')
      END IF
!
!     ==========================================================================
!     == RADIAL GRID                                                          ==
!     ==========================================================================
      CALL RADIAL$R(GID,NR,R)
      RAD=POTPAR(ISP)%RAD
      DO IR=1,NR
        IRAD=IR
        IF(R(IR).GT.RAD) EXIT
      ENDDO
!
!     ==========================================================================
!     == COLLECT DATA                                                         ==
!     ==========================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('LNX',LNX)
      IF(LNXPHI.NE.LNX) THEN
        CALL ERROR$STOP('INCONSISTENT #(PARTIAL WAVES)')
        CALL ERROR$STOP('LMTO$DOLOCORB')
      END IF
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      CALL SETUP$GETR8A('PRO',NR*LNX,PRO)
      LX=MAXVAL(LOX(:))
      CALL SETUP$UNSELECT()
!
!     ==========================================================================
!     == FIND ONSITE STRUCTURE CONSTANTS                                      ==
!     ==========================================================================
      ALLOCATE(SBARAV(LX+1))      
      SBARAV(:)=0.D0
      NNS=SIZE(SBAR)
      TCHK=.FALSE.
      DO IIB=1,NNS
        IF(SBAR(IIB)%IAT1.NE.IAT) CYCLE
        IF(SBAR(IIB)%IAT2.NE.IAT) CYCLE
        IF(SUM(SBAR(IIB)%IT**2).NE.0) CYCLE
        DO L=0,LX
          IORB=SBARLI1(L+1,ISP)
          IF(IORB.LE.0) CYCLE
          DO IM=1,2*L+1 
            SBARAV(L+1)=SBARAV(L+1)+SBAR(IIB)%MAT(IORB-1+IM,IORB-1+IM)
          ENDDO
          SBARAV(L+1)=SBARAV(L+1)/REAL(2*L+1,KIND=8)
        ENDDO
!       ------------------------------------------------------------------------
        TCHK=.TRUE.
        EXIT
      ENDDO
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('ONSITE TERM OF STRUCTURE CONSTANTS NOT FOUND')
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$STOP('LMTO$DOLOCORB')
      END IF
!
!     ==========================================================================
!     == CONSTRUCT LOCAL ORBITALS                                             ==
!     ==========================================================================
      ALLOCATE(AECHI(NR,LNXPHI))
      ALLOCATE(PSCHI(NR,LNXPHI))
      ALLOCATE(NLCHI(NR,LNXPHI))
      DO LN=1,LNXPHI
        L=LOX(LN)
        LN1=POTPAR(ISP)%TAILED%LNDOT(LN)
        AECHI(:,LN)=POTPAR(ISP)%TAILED%AEF(:,LN) &
    &              -POTPAR(ISP)%TAILED%AEF(:,LN1)*SBARAV(L+1)
        PSCHI(:,LN)=POTPAR(ISP)%TAILED%PSF(:,LN) &
    &              -POTPAR(ISP)%TAILED%PSF(:,LN1)*SBARAV(L+1)
        NLCHI(:,LN)=POTPAR(ISP)%TAILED%NLF(:,LN) &
    &              -POTPAR(ISP)%TAILED%NLF(:,LN1)*SBARAV(L+1)
      ENDDO
!!$WRITE(STRING,FMT='(I5)')IAT
!!$STRING='LOCORB_FORATOM'//TRIM(ADJUSTL(STRING))//'.DAT'
!!$CALL LMTO_WRITEPHI(TRIM(STRING),GID,NR,LNXPHI,AECHI)
!
!     ==ORTHONORMALIZE LOCAL ORBITALS ==========================================
!     == ORTHONORMALIZATION IS NOT REQUIRED AND SERVES ONLY ESTAETICAL PURPOSES
!!$      DO LN=1,LNX
!!$        L=LOX(LN)
!!$        DO LN1=1,LN-1
!!$          IF(LOX(LN1).NE.L) CYCLE
!!$          AUX(:)=R(:)**2*AECHI(:,LN)*AECHI(:,LN1)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!!$          AECHI(:,LN)=AECHI(:,LN)-AECHI(:,LN1)*VAL
!!$          PSCHI(:,LN)=PSCHI(:,LN)-PSCHI(:,LN1)*VAL
!!$        ENDDO
!!$        AUX(:)=R(:)**2*AECHI(:,LN)**2
!!$        CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!!$        VAL=1.D0/SQRT(VAL)
!!$        AECHI(:,LN)=AECHI(:,LN)*VAL
!!$        PSCHI(:,LN)=PSCHI(:,LN)*VAL
!!$      ENDDO
!
!     ==========================================================================
!     == TRANSFORMATION OF PROJECTORS                                         ==
!     ==========================================================================
      ALLOCATE(AMAT(LNX,LNX))
      ALLOCATE(BMAT(LNX,LNX))
      ALLOCATE(XMAT(LNX,LNX))
      AMAT(:,:)=0.D0
      DO LN=1,LNX
        DO LN1=1,LNX
          IF(LOX(LN).NE.LOX(LN1)) CYCLE
          AUX(:)=R(:)**2*PRO(:,LN)*PSCHI(:,LN1)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,AMAT(LN,LN1))
        ENDDO
      ENDDO
      AMAT=TRANSPOSE(AMAT)
      BMAT(:,:)=0.D0
      DO LN=1,LNX
        BMAT(LN,LN)=1.D0
      ENDDO
      CALL LIB$MATRIXSOLVER8(LNX,LNX,LNX,AMAT,XMAT,BMAT)
      AMAT=TRANSPOSE(XMAT)
      DEALLOCATE(XMAT)
      DEALLOCATE(BMAT)
!
!     ==========================================================================
!     == DELETE ORBITALS NOT IN THE SET                                       ==
!     ==========================================================================
      LNCHI=0
      DO LN=1,LNX
        IF(.NOT.TORB(LN)) CYCLE
        LNCHI=LNCHI+1
        CHI(:,LNCHI)=AECHI(:,LN)
        CHIPHI(LNCHI,:)=AMAT(LN,:)   ! MATCHING COEFFICIENTS
      ENDDO
      DEALLOCATE(AMAT)
!
!     ==========================================================================
!     == PLOT LOCAL ORBITALS                                                  ==
!     ==========================================================================
!!$      CALL SETUP$ISELECT(ISP)
!!$      CALL SETUP$GETR8('AEZ',AEZ)
!!$      WRITE(STRING,FMT='(F3.0)')AEZ
!!$      STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
!!$      CALL SETUP_WRITEPHI(-'CHI'//TRIM(STRING),GID,NR,LNCHI,CHI)
!!$      CALL SETUP$ISELECT(0)
!
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$DOLOCORB(IAT,ISP,GID,NR,LNXCHI,LNXPHI,TORB,CHIPHI,CHI)
!     **************************************************************************
!     **  CONSTRUCTS ONSITE-MAPPING FROM PARTIAL WAVES TO LOCAL ORBITALS      **
!     **                                                                      **
!     **  TORB SELECTS LOCAL ORBITALS FROM PARTIAL WAVES                      **
!     **                                                                      **
!     **  TAILS ARE DEFINED BY SCATTERING WAVE FUNCTION                       **
!     **      FOR THE VALENCE STATE (I.E.ISCATT=0)                            **
!     **                                                                      **
!     **  OFF-SITE TERMS ARE REPLACED BY EXPONENTIAL TAIL MATCHED DIFFERENTIABLY
!     **                                                                      **
!     **  ON SITE STRUCTURE CONSTANTS ARE SPHERICALLY AVERAGED                **
!     **                                                                      **
!     **   |PSI>= SUM_{I,J} |CHI_I>*CHIPHI(I,J)*<PTILDE_J|\PSITILDE>          **
!     **       WHERE J INCLUDES ONLY ONSITE TERMS                             **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2011 ************
      USE LMTO_MODULE, ONLY : SBAR,POTPAR,SBARLI1,K2
      USE STRINGS_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT     ! ATOM INDEX
      INTEGER(4),INTENT(IN)  :: GID     ! GRID ID
      INTEGER(4),INTENT(IN)  :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN)  :: LNXPHI  ! #(PARTIAL WAVES W/O M-MULTIPLICITY)
      INTEGER(4),INTENT(IN)  :: LNXCHI  ! #(PARTIAL WAVES W/O M-MULTIPLICITY)
      LOGICAL(4),INTENT(IN)  :: TORB(LNXPHI)          ! SELECTS LOCAL ORBITALS
      REAL(8)   ,INTENT(OUT) :: CHIPHI(LNXCHI,LNXPHI) !<PI_I|=\SUM_J CHIPHI(I,J)<P_J|
      REAL(8)   ,INTENT(OUT) :: CHI(NR,LNXCHI)
      INTEGER(4)             :: ISP         ! ATOM TYPE INDEX
      INTEGER(4)             :: LNX         ! #(PARTIAL WAVES)
      INTEGER(4)             :: NNS         ! 
      INTEGER(4)             :: LOX(LNXPHI) ! ANGULAR MOMENTA
      INTEGER(4)             :: ISCATT(LNXPHI)    ! COUNTER RELATIVE TO HOMO
      REAL(8)                :: AEPHI(NR,LNXPHI)  ! AE PARTIAL WAVES
      REAL(8)                :: PSPHI(NR,LNXPHI)  ! PSEUDO PARTIAL WAVES
      REAL(8)                :: AEPHIDOT(NR,LNXPHI) ! AE PARTIAL WAVES
      REAL(8)                :: PSPHIDOT(NR,LNXPHI) ! PARTIAL WAVES
      REAL(8)                :: NLPHI(NR,LNXPHI) ! PARTIAL WAVES
      REAL(8)                :: NLPHIDOT(NR,LNXPHI) ! PARTIAL WAVES
      REAL(8)                :: PRO(NR,LNXPHI)    ! PROJECTOR FUNCTIONS
      REAL(8)   ,ALLOCATABLE :: AECHI(:,:)  ! ALL-ELECTRON LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: PSCHI(:,:)  ! PSEUDO LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: NLCHI(:,:)  ! NODELESS LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: SBARAV(:)
      REAL(8)   ,ALLOCATABLE :: AMAT(:,:)   !(LNXCHI1,LNXPHI) 
      REAL(8)   ,ALLOCATABLE :: BMAT(:,:)   !(LNXCHI1,LNXCHI1) 
      REAL(8)   ,ALLOCATABLE :: XMAT(:,:)   !(LNXPHI,LNXCHI) 
      REAL(8)                :: RAD              ! COVALENT RADIUS
      REAL(8)                :: AEZ               ! ATOMIC NUMBER
      REAL(8)                :: AUX(NR)
      REAL(8)                :: R(NR)
      REAL(8)                :: SVAR,SVAR1,SVAR2,VAL,DER
      REAL(8)                :: KVAL,KDER,JVAL,JDER
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: LN,LN1,L,IIB,LNCHI,IR,IM
      INTEGER(4)             :: LX
      INTEGER(4)             :: IORB
      INTEGER(4)             :: IRAD  ! GRID INDEX JUST BEYOND RAD
      CHARACTER(64)          :: STRING
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO$DOLOCORB')
!
!     ==========================================================================
!     == CHECK CONSISTENCY OF INPUT                                           ==
!     ==========================================================================
      LNCHI=0
      DO LN=1,LNXPHI
        IF(TORB(LN))LNCHI=LNCHI+1
      ENDDO
      IF(LNCHI.NE.LNXCHI) THEN
        CALL ERROR$MSG('LOCAL-ORBITAL SELECTION TORB INCONSISTENT WITH LNXCHI')
        CALL ERROR$L4VAL('TORB',TORB)
        CALL ERROR$I4VAL('LNXCHI',LNXCHI)
        CALL ERROR$STOP('LMTO$DOLOCORB')
      END IF
!
!     ==========================================================================
!     == RADIAL GRID                                                          ==
!     ==========================================================================
      CALL RADIAL$R(GID,NR,R)
      RAD=POTPAR(ISP)%RAD
      DO IR=1,NR
        IRAD=IR
        IF(R(IR).GT.RAD) EXIT
      ENDDO
!
!     ==========================================================================
!     == COLLECT DATA                                                         ==
!     ==========================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('LNX',LNX)
      IF(LNXPHI.NE.LNX) THEN
        CALL ERROR$STOP('INCONSISTENT #(PARTIAL WAVES)')
        CALL ERROR$STOP('LMTO$DOLOCORB')
      END IF
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUP$GETR8A('PSPHI',NR*LNX,PSPHI)
      CALL SETUP$GETR8A('QPHI',NR*LNX,NLPHI)
      CALL SETUP$GETR8A('AEPHIDOT',NR*LNX,AEPHIDOT)
      CALL SETUP$GETR8A('PSPHIDOT',NR*LNX,PSPHIDOT)
      CALL SETUP$GETR8A('QPHIDOT',NR*LNX,NLPHIDOT)
      CALL SETUP$GETR8A('PRO',NR*LNX,PRO)
      CALL SETUP$ISELECT(0)
!     == SELECT THE CORRECT PHIDOT FUNCTION (ONLY ONE PER L) ===================
      DO LN=1,LNX
        LN1=POTPAR(ISP)%LNSCATT(LN)
        IF(LN1.EQ.LN) CYCLE
        AEPHIDOT(:,LN)=AEPHIDOT(:,LN1)
        PSPHIDOT(:,LN)=PSPHIDOT(:,LN1)
        NLPHIDOT(:,LN)=NLPHIDOT(:,LN1)
      ENDDO
!
      LX=MAXVAL(LOX(:))
!
!     ==========================================================================
!     == FIND ONSITE STRUCTURE CONSTANTS                                      ==
!     ==========================================================================
      ALLOCATE(SBARAV(LX+1))      
      SBARAV(:)=0.D0
      NNS=SIZE(SBAR)
      TCHK=.FALSE.
      DO IIB=1,NNS
        IF(SBAR(IIB)%IAT1.NE.IAT) CYCLE
        IF(SBAR(IIB)%IAT2.NE.IAT) CYCLE
        IF(SUM(SBAR(IIB)%IT**2).NE.0) CYCLE
        DO L=0,LX
          IORB=SBARLI1(L+1,ISP)
          IF(IORB.LE.0) CYCLE
          DO IM=1,2*L+1 
            SBARAV(L+1)=SBARAV(L+1)+SBAR(IIB)%MAT(IORB-1+IM,IORB-1+IM)
          ENDDO
          SBARAV(L+1)=SBARAV(L+1)/REAL(2*L+1,KIND=8)
        ENDDO
!       ------------------------------------------------------------------------
        TCHK=.TRUE.
        EXIT
      ENDDO
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('ONSITE TERM OF STRUCTURE CONSTANTS NOT FOUND')
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$STOP('LMTO$DOLOCORB')
      END IF
!
!     ==========================================================================
!     == COUNT ONSITE ORBITALS BEFORE EXCLUSION                               ==
!     ==========================================================================
      ALLOCATE(AECHI(NR,LNXPHI))
      ALLOCATE(PSCHI(NR,LNXPHI))
      ALLOCATE(NLCHI(NR,LNXPHI))
!
!     ==========================================================================
!     == CONSTRUCT LOCAL ORBITALS                                             ==
!     ==========================================================================
      LNCHI=0
      DO LN=1,LNXPHI
        L=LOX(LN)
        SVAR1=POTPAR(ISP)%KTOPHI(LN)
        SVAR2=POTPAR(ISP)%KTOPHIDOT(LN) &
    &        -POTPAR(ISP)%JBARTOPHIDOT(LN)*SBARAV(L+1)
        AECHI(:,LN)=AEPHI(:,LN)*SVAR1+AEPHIDOT(:,LN)*SVAR2
        PSCHI(:,LN)=PSPHI(:,LN)*SVAR1+PSPHIDOT(:,LN)*SVAR2
        NLCHI(:,LN)=NLPHI(:,LN)*SVAR1+NLPHIDOT(:,LN)*SVAR2
!
!       == ATTACH EXPONENTIAL TAIL AT THE MATCHING RADIUS ======================
        CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
        CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
        VAL=KVAL-(JVAL-KVAL*POTPAR(ISP)%QBAR(LN))*SBARAV(L+1)
        DER=KDER-(JDER-KDER*POTPAR(ISP)%QBAR(LN))*SBARAV(L+1)
        SVAR=DER/VAL
        IF(SVAR.GT.0.D0) THEN
          CALL SETUP_WRITEPHI(-'FAILEDNLORBITAL.DAT',GID,NR,1,NLCHI(:,LN))
          CALL ERROR$MSG('ORBITAL DOES NOT DECAY')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$I4VAL('ISP',ISP)
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('LN',LN)
          CALL ERROR$R8VAL('MATCHING RADIUS',RAD)
          CALL ERROR$R8VAL('LOGARITMIC DERIVATIVE ',SVAR)
          CALL ERROR$STOP('LMTO$DOLOCORB')
        END IF
        AUX(:)=0.D0
        AUX(IRAD:)=-NLCHI(IRAD:,LN)+VAL*EXP(SVAR*(R(IRAD:)-RAD))
        AECHI(IRAD:,LN)=AECHI(IRAD:,LN)+AUX(IRAD:)
        PSCHI(IRAD:,LN)=PSCHI(IRAD:,LN)+AUX(IRAD:)
        NLCHI(IRAD:,LN)=NLCHI(IRAD:,LN)+AUX(IRAD:)
      ENDDO
!
!     ==ORTHONORMALIZE LOCAL ORBITALS ==========================================
!     == ORTHONORMALIZATION IS NOT REQUIRED AND SERVES ONLY ESTAETICAL PURPOSES
!!$      DO LN=1,LNX
!!$        L=LOX(LN)
!!$        DO LN1=1,LN-1
!!$          IF(LOX(LN1).NE.L) CYCLE
!!$          AUX(:)=R(:)**2*AECHI(:,LN)*AECHI(:,LN1)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!!$          AECHI(:,LN)=AECHI(:,LN)-AECHI(:,LN1)*VAL
!!$          PSCHI(:,LN)=PSCHI(:,LN)-PSCHI(:,LN1)*VAL
!!$        ENDDO
!!$        AUX(:)=R(:)**2*AECHI(:,LN)**2
!!$        CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!!$        VAL=1.D0/SQRT(VAL)
!!$        AECHI(:,LN)=AECHI(:,LN)*VAL
!!$        PSCHI(:,LN)=PSCHI(:,LN)*VAL
!!$      ENDDO
!
!     ==========================================================================
!     == TRANSFORMATION OF PROJECTORS                                         ==
!     ==========================================================================
      ALLOCATE(AMAT(LNX,LNX))
      ALLOCATE(BMAT(LNX,LNX))
      ALLOCATE(XMAT(LNX,LNX))
      AMAT(:,:)=0.D0
      DO LN=1,LNX
        DO LN1=1,LNX
          IF(LOX(LN).NE.LOX(LN1)) CYCLE
          AUX(:)=R(:)**2*PRO(:,LN)*PSCHI(:,LN1)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,AMAT(LN,LN1))
        ENDDO
      ENDDO
      AMAT=TRANSPOSE(AMAT)
      BMAT(:,:)=0.D0
      DO LN=1,LNX
        BMAT(LN,LN)=1.D0
      ENDDO
      CALL LIB$MATRIXSOLVER8(LNX,LNX,LNX,AMAT,XMAT,BMAT)
      AMAT=TRANSPOSE(XMAT)
      DEALLOCATE(XMAT)
      DEALLOCATE(BMAT)
!
!     ==========================================================================
!     == DELETE ORBITALS NOT IN THE SET                                       ==
!     ==========================================================================
      LNCHI=0
      DO LN=1,LNX
        IF(.NOT.TORB(LN)) CYCLE
        LNCHI=LNCHI+1
        CHI(:,LNCHI)=AECHI(:,LN)
        CHIPHI(LNCHI,:)=AMAT(LN,:)   ! MATCHING COEFFICIENTS
      ENDDO
!
!     ==========================================================================
!     == PLOT LOCAL ORBITALS                                                  ==
!     ==========================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETR8('AEZ',AEZ)
      WRITE(STRING,FMT='(F3.0)')AEZ
      STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
      CALL SETUP_WRITEPHI(-'CHI'//TRIM(STRING),GID,NR,LNCHI,CHI)
      CALL SETUP$ISELECT(0)
!
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$REPORTOVERLAP(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NNS
      INTEGER(4)            :: NN,LM1
      REAL(8)               :: Q,Q1,S,S1
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTOVERLAP')
      NNS=SIZE(OVERLAP)
      WRITE(NFIL,FMT='(82("="))')
      WRITE(NFIL,FMT='(82("="),T10,"   OVERLAP MATRIX ELEMENTS   ")')
      WRITE(NFIL,FMT='(82("="))')
      Q=0.D0
      DO NN=1,NNS
        WRITE(NFIL,FMT='(82("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 OVERLAP(NN)%IAT1,OVERLAP(NN)%IAT2,OVERLAP(NN)%IT
        Q1=0.D0
        S1=0.D0
        DO LM1=1,OVERLAP(NN)%N1
          WRITE(NFIL,FMT='(20E10.2)')OVERLAP(NN)%MAT(LM1,:)
          Q1=Q1+SUM(OVERLAP(NN)%MAT(LM1,:)*DENMAT(NN)%MAT(LM1,:,1))
          S1=S1+SUM(OVERLAP(NN)%MAT(LM1,:)*DENMAT(NN)%MAT(LM1,:,4))
        ENDDO
        WRITE(NFIL,*)'NUMBER OF ELECTRONS ATTRIBUTED TO THIS PAIR ',Q1," SPIN=",S1
        Q=Q+Q1
        S=S+S1
        WRITE(NFIL,FMT='(82("="))')
      ENDDO
      WRITE(NFIL,*)'NUMBER OF ELECTRONS',Q," SPIN=",S
                             CALL TRACE$POP()
      RETURN
      END SUBROUTINE LMTO$REPORTOVERLAP
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$REPORTSBAR(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NNS
      INTEGER(4)            :: NN,LM1
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTSBAR')
      NNS=SIZE(SBAR)
      WRITE(NFIL,FMT='(82("="))')
      WRITE(NFIL,FMT='(82("="),T10,"   SCREENED STRUCTURE CONSTANTS   ")')
      WRITE(NFIL,FMT='(82("="))')
      DO NN=1,NNS
        WRITE(NFIL,FMT='(82("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 SBAR(NN)%IAT1,SBAR(NN)%IAT2,SBAR(NN)%IT
        DO LM1=1,SBAR(NN)%N1
          WRITE(NFIL,FMT='(20F10.5)')SBAR(NN)%MAT(LM1,:)
        ENDDO
        WRITE(NFIL,FMT='(82("="))')
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END SUBROUTINE LMTO$REPORTSBAR
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$REPORTORTHODENMAT(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : DENMAT,PERIODICMAT2_TYPE,ISPECIES,LOX,LNX
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NND
      INTEGER(4)            :: NAT
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: LMNX
      REAL(8)   ,ALLOCATABLE:: T(:,:),UNT(:,:)
      INTEGER(4)            :: NN,LM1,IAT,I3,ISP
      TYPE(PERIODICMAT2_TYPE),ALLOCATABLE :: DENMAT1(:)
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTDENMAT')
      NND=SIZE(DENMAT)
      NAT=SIZE(ISPECIES)
!
!     ==========================================================================
!     == CREATE A NEW TEMPORARY ARRAY FOR THE DENSITY MATRIX                  ==
!     ==========================================================================
      ALLOCATE(DENMAT1(NND))
      DO NN=1,NND
        N1=DENMAT(NN)%N1
        N2=DENMAT(NN)%N2
        N3=DENMAT(NN)%N3
        ALLOCATE(DENMAT1(NN)%MAT(N1,N2,N3))
        DENMAT1(NN)%MAT=DENMAT(NN)%MAT
      ENDDO
!
!     ==========================================================================
!     == TRANSFORM TO ONSITE-ORTHOGONALIZED BASOS                             ==
!     ==========================================================================
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LMNX=SUM(2*LOX(:LNX(ISP),ISP)+1)
        ALLOCATE(T(LMNX,LMNX))
        ALLOCATE(UNT(LMNX,LMNX))
        CALL LMTO_ONSORTHO(IAT,LMNX,T,UNT)
!!$DO LM1=1,LMNX
!!$  WRITE(*,FMT='("T  =",I5,30F10.4)')LM1,T(LM1,:)
!!$ENDDO
!!$DO LM1=1,LMNX
!!$  WRITE(*,FMT='("UNT=",I5,30F10.4)')LM1,UNT(LM1,:)
!!$ENDDO
        DO NN=1,NND
          N3=DENMAT(NN)%N3
          IF(DENMAT(NN)%IAT1.EQ.IAT) THEN
            DO I3=1,N3
              DENMAT1(NN)%MAT(:,:,I3)=MATMUL(UNT,DENMAT1(NN)%MAT(:,:,I3))
            ENDDO
          END IF
          IF(DENMAT(NN)%IAT2.EQ.IAT) THEN
            DO I3=1,N3
              DENMAT1(NN)%MAT(:,:,I3)=MATMUL(DENMAT1(NN)%MAT(:,:,I3),TRANSPOSE(UNT))
            ENDDO
          END IF
        ENDDO
        DEALLOCATE(T)
        DEALLOCATE(UNT)
      ENDDO 
!
!     ==========================================================================
!     == WRITE DENSITY MATRIX                                                 ==
!     ==========================================================================
      WRITE(NFIL,FMT='(82("="))')
      WRITE(NFIL,FMT='(82("="),T10," ONSITE ORTHOGONALIZED DENSITY MATRIX ")')
      WRITE(NFIL,FMT='(82("="))')
      DO NN=1,NND
        WRITE(NFIL,FMT='(82("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 DENMAT(NN)%IAT1,DENMAT(NN)%IAT2,DENMAT(NN)%IT
        DO LM1=1,DENMAT(NN)%N1
          WRITE(NFIL,FMT='(20F10.5)')DENMAT1(NN)%MAT(LM1,:,1)
        ENDDO
        WRITE(NFIL,FMT='(82("-"),T10," SPIN SZ CONTRIBTION ")')

        DO LM1=1,DENMAT(NN)%N1
          WRITE(NFIL,FMT='(20F10.5)')DENMAT1(NN)%MAT(LM1,:,4)
        ENDDO

        WRITE(NFIL,FMT='(82("="))')
      ENDDO
!
!     ==========================================================================
!     == CLEAN UP DENMAT1                                                     ==
!     ==========================================================================
      DO NN=1,NND
        DEALLOCATE(DENMAT1(NN)%MAT)
      ENDDO
      DEALLOCATE(DENMAT1)

                             CALL TRACE$POP()
      RETURN
      END SUBROUTINE LMTO$REPORTORTHODENMAT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$REPORTDENMAT(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : DENMAT
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NND
      INTEGER(4)            :: NN,LM1
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTDENMAT')
      NND=SIZE(DENMAT)
      WRITE(NFIL,FMT='(82("="))')
      WRITE(NFIL,FMT='(82("="),T10,"   DENSITY MATRIX   ")')
      WRITE(NFIL,FMT='(82("="))')
      DO NN=1,NND
        WRITE(NFIL,FMT='(82("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 DENMAT(NN)%IAT1,DENMAT(NN)%IAT2,DENMAT(NN)%IT
        DO LM1=1,DENMAT(NN)%N1
          WRITE(NFIL,FMT='(20F10.5)')DENMAT(NN)%MAT(LM1,:,1)
        ENDDO
        WRITE(NFIL,FMT='(82("-"),T10," SPIN SZ CONTRIBTION ")')

        DO LM1=1,DENMAT(NN)%N1
          WRITE(NFIL,FMT='(20F10.5)')DENMAT(NN)%MAT(LM1,:,4)
        ENDDO

        WRITE(NFIL,FMT='(82("="))')
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END SUBROUTINE LMTO$REPORTDENMAT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$REPORTHAMIL(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : HAMIL
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NND
      INTEGER(4)            :: NN,LM1
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTHAMIL')
      NND=SIZE(HAMIL)
      WRITE(NFIL,FMT='(82("="))')
      WRITE(NFIL,FMT='(82("="),T10,"   HAMILTONIAN   ")')
      WRITE(NFIL,FMT='(82("="))')
      DO NN=1,NND
        WRITE(NFIL,FMT='(82("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 HAMIL(NN)%IAT1,HAMIL(NN)%IAT2,HAMIL(NN)%IT
        DO LM1=1,HAMIL(NN)%N1
          WRITE(NFIL,FMT='(20E12.2)')HAMIL(NN)%MAT(LM1,:,1)
        ENDDO
        WRITE(NFIL,FMT='(82("="))')
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END SUBROUTINE LMTO$REPORTHAMIL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SYMMETRIZEHAMIL()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : HAMIL
      IMPLICIT NONE
      INTEGER(4)            :: NND
      INTEGER(4)            :: NN1,NN2,IAT1,IAT2,IT(3),I
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$SYMMETRIZEHAMIL')
      NND=SIZE(HAMIL)
      DO NN1=1,NND
        IAT1=HAMIL(NN1)%IAT1
        IAT2=HAMIL(NN1)%IAT2
        IT(:)=HAMIL(NN1)%IT(:)
        DO NN2=NN1,NND
          IF(HAMIL(NN2)%IAT1.NE.IAT2) CYCLE
          IF(HAMIL(NN2)%IAT2.NE.IAT1) CYCLE
          IF(SUM((HAMIL(NN2)%IT+IT)**2).GT.0) CYCLE
          DO I=1,HAMIL(NN1)%N3
            HAMIL(NN1)%MAT(:,:,I)=0.5D0 &
     &                 *(HAMIL(NN1)%MAT(:,:,I)+TRANSPOSE(HAMIL(NN2)%MAT(:,:,I)))
            HAMIL(NN2)%MAT(:,:,I)=TRANSPOSE(HAMIL(NN1)%MAT(:,:,I))
          ENDDO
        ENDDO
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END SUBROUTINE LMTO$SYMMETRIZEHAMIL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$REPORTPOTBAR(NFIL)
!     **************************************************************************
!     ** WRITE INFORMATION RELATED TO POTENTIAL PARAMETERS                    **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR,NSP,LNX,LOX
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      CHARACTER(64)         :: TITLE
      INTEGER(4)            :: ISP,LN
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTPOTPAR')
      DO ISP=1,NSP
        DO LN=1,LNX(ISP)
          WRITE(TITLE,FMT='(I3," LN=",I3," L=",I2)')ISP,LN,LOX(LN,ISP)
          TITLE='POTENTIAL PARAMETERS FOR ATOM TYPE '//TRIM(ADJUSTL(TITLE))
          CALL REPORT$TITLE(NFIL,TITLE)
          CALL REPORT$R8VAL(NFIL,'RAD',POTPAR(ISP)%RAD,'ABOHR')
          CALL REPORT$I4VAL(NFIL,'LNSCATT',POTPAR(ISP)%LNSCATT(LN),'')
          CALL REPORT$R8VAL(NFIL,'QBAR',POTPAR(ISP)%QBAR(LN),'')
          CALL REPORT$R8VAL(NFIL,'KTOPHI',POTPAR(ISP)%KTOPHI(LN),'')
          CALL REPORT$R8VAL(NFIL,'KTOPHIDOT',POTPAR(ISP)%KTOPHIDOT(LN),'')
          CALL REPORT$R8VAL(NFIL,'JBARTOPHIDOT',POTPAR(ISP)%JBARTOPHIDOT(LN),'')
          CALL REPORT$R8VAL(NFIL,'<PRO|PSPHIDOT>',POTPAR(ISP)%PHIDOTPROJ(LN),'')
        ENDDO
        WRITE(NFIL,FMT='(82("="))')
      ENDDO
                                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_WRITEPHI(FILE,GID,NR,NPHI,PHI)
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
!===============================================================================
!===============================================================================
!===============================================================================
!===================    PLOT ORBITAL ROUTINES       ============================
!===============================================================================
!===============================================================================
!===============================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GRIDPLOT_TAILED(IAT0)
!     **************************************************************************
!     **  WRITES THE ENVELOPE FUNCTIONS CENTERED AT ATOM IAT0 TO FILE         **
!     **                                                                      **
!     **  SET N1,N2,N3 EQUAL TO THE NUMBER OF GRIDPOINTS IN EACH DIRECTION    **
!     **                                                                      **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2009 ************
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      USE LMTO_MODULE, ONLY : ISPECIES,SBAR,POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT0
      INTEGER(4),PARAMETER  :: N1=40,N2=40,N3=40 !GRID (1D?)
      INTEGER(4),PARAMETER  :: NRAD=200
      INTEGER(4),PARAMETER  :: NDIRX=100
      REAL(8)   ,PARAMETER  :: RANGE=8.D0
      REAL(8)               :: DIR(3,NDIRX)
      REAL(8)               :: ORIGIN(3)
      REAL(8)               :: TVEC(3,3)    !BOX
      REAL(8)               :: RBAS(3,3)
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IORB
      INTEGER(4)            :: NORB         !#YLM
      REAL(8)   ,ALLOCATABLE:: R0(:,:)      !(3,NAT)
      REAL(8)   ,ALLOCATABLE:: RAD(:)      !(NAT)
      REAL(8)               :: DR(3)
      INTEGER(4)            :: NNS
      INTEGER(4)            :: NN,IAT,IAT2,ISP,I1,I2,I3,LN,I,J
      INTEGER(4)            :: L
      INTEGER(4)            :: LNX
      INTEGER(4)            :: LMNXS,LMNX
      REAL(8)               :: AEZ
      REAL(8)   ,ALLOCATABLE:: RCLUSTER(:,:) ! POSITIONS OF NEIGBORS
      REAL(8)   ,ALLOCATABLE:: ZCLUSTER(:)   ! ATOMIC NUMBER OF NEIGBORS
      REAL(8)   ,ALLOCATABLE:: ORB(:,:)
      INTEGER(4),ALLOCATABLE:: LOX(:)
      INTEGER(4),ALLOCATABLE:: ISCATT(:)
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: NATCLUSTER !#(ATOMS ON THE CLUSTER)
      CHARACTER(64)         :: FILE
      CHARACTER(15)         :: STRING         !CONTAINS IATO
      CHARACTER(16)         :: FORMATTYPE
      REAL(8)               :: X1D(NRAD)
      INTEGER(4)            :: NDIR,IDIR
      INTEGER(4)            :: NP,IP
      REAL(8)  ,ALLOCATABLE :: P(:,:)
      LOGICAL(4),PARAMETER  :: TGAUSS=.TRUE.
      CHARACTER(2),PARAMETER :: ID='3D'
!     **************************************************************************
                                         CALL TRACE$PUSH('LMTO_GRIDPLOT_TAILED')
!
!     ==========================================================================
!     == COLLECT INFORMATION                                                  ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      NAT=SIZE(ISPECIES)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      ALLOCATE(RAD(NAT))
      NORB=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        RAD(IAT)=POTPAR(ISP)%RAD
!
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(ISCATT(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
        DO LN=1,LNX
          IF(ISCATT(LN).NE.0) CYCLE
          L=LOX(LN)
          IF(IAT.EQ.IAT0)NORB=MAX(NORB,(L+1)**2)
        ENDDO
        DEALLOCATE(LOX)
        DEALLOCATE(ISCATT)
        CALL SETUP$ISELECT(0)
      ENDDO
!
!     ==========================================================================
!     == DEFINE ATOMS ON THE CLUSTER OF NEIGHBORS                             ==
!     ==========================================================================
      NATCLUSTER=0
      NNS=SIZE(SBAR)  !SBAR STRUCCONS. (GLOBAL)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.EQ.IAT0) NATCLUSTER=NATCLUSTER+1
      ENDDO
      ALLOCATE(ZCLUSTER(NATCLUSTER))
      ALLOCATE(RCLUSTER(3,NATCLUSTER))
      IAT=0
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT0) CYCLE
        IAT=IAT+1
        IAT2=SBAR(NN)%IAT2
        RCLUSTER(:,IAT)=R0(:,IAT2) &
     &                 +RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
     &                 +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
     &                 +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
        ISP=ISPECIES(IAT2)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL SETUP$ISELECT(0)
        ZCLUSTER(IAT)=AEZ
      ENDDO
!
!     ==========================================================================
!     == DEFINE GRID POINTS                                                   ==
!     ==========================================================================
      IF(ID.EQ.'3D') THEN
        NP=N1*N2*N3
        ALLOCATE(P(3,NP))
        CALL LMTO_GRIDORB_CUBEGRID(R0(:,IAT0),RANGE,N1,N2,N3,ORIGIN,TVEC,P)
      ELSE IF(ID.EQ.'1D') THEN
        NDIR=0
        DO IAT=1,NATCLUSTER
          NDIR=NDIR+1
          IF(NDIR.GT.NDIRX) THEN
            CALL ERROR$MSG('#(NEIGHBORS EXCEEDS MAXIMUM')
            CALL ERROR$STOP('LMTO_GRIDPLOT_TAILED')
          END IF
          DIR(:,NDIR)=RCLUSTER(:,IAT)-R0(:,IAT0)
          IF(SQRT(SUM(DIR(:,NDIR)**2)).LT.1.D-3) THEN
            NDIR=NDIR-1
            CYCLE   
          END IF
        ENDDO
        NP=NRAD*NDIR
        ALLOCATE(P(3,NP))
        CALL LMTO_GRIDORB_STARGRID(R0(:,IAT0),RANGE,NDIR,DIR,NRAD,X1D,P)
      ELSE
        CALL ERROR$MSG('INVALID ID')
        CALL ERROR$MSG('MUST BE "1D", "2D", OR "3D"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO_GRIDPLOT_TAILED')
      END IF
!
!     ==========================================================================
!     == DETERMINE ENVELOPE FUNCTION AT THE GRID POINTS                       ==
!     ==========================================================================
      ALLOCATE(ORB(NP,NORB))
      CALL LMTO_LMTO_GRIDPLOT_TAILEDINNER(NAT,R0,IAT0,NP,P,NORB,ORB)
!
!     ==========================================================================
!     == WRITE ORBITALS TO FILE                                               ==
!     ==========================================================================
      IF(ID.EQ.'3D') THEN
        DO IORB=1,NORB
          FILE='NTB3D'
          WRITE(STRING,*)IAT0 
          FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))
          WRITE(STRING,*)IORB 
          FILE=TRIM(ADJUSTL(FILE))//'_'//TRIM(ADJUSTL(STRING))//'.CUB'
          CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
          CALL FILEHANDLER$UNIT('HOOK',NFIL)
!
          CALL LMTO_WRITECUBEFILE(NFIL,NATCLUSTER,ZCLUSTER,RCLUSTER &
     &                           ,ORIGIN,TVEC,N1,N2,N3,ORB(:,IORB))
          CALL FILEHANDLER$CLOSE('HOOK')
          CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
        ENDDO
!
!     ==========================================================================
      ELSE IF(ID.EQ.'1D') THEN
        DO IORB=1,NORB
          FILE='NTB1D'
          WRITE(STRING,*)IAT0 
          FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))
          WRITE(STRING,*)IORB 
          FILE=TRIM(ADJUSTL(FILE))//'_'//TRIM(ADJUSTL(STRING))//'.DAT'
          CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
          CALL FILEHANDLER$UNIT('HOOK',NFIL)
          DO I=1,NRAD
            WRITE(NFIL,*)X1D(I),(ORB(NRAD*(IDIR-1)+I,IORB),IDIR=1,NDIR)
          ENDDO
          CALL FILEHANDLER$CLOSE('HOOK')
          CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
        ENDDO
!
!     ==========================================================================
      ELSE IF(ID.EQ.'2D') THEN
        DO IORB=1,NORB
          FILE='NTB2D'
          WRITE(STRING,*)IAT0 
          FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))
          WRITE(STRING,*)IORB 
          FILE=TRIM(ADJUSTL(FILE))//'_'//TRIM(ADJUSTL(STRING))//'.DAT'
          CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
          CALL FILEHANDLER$UNIT('HOOK',NFIL)
          DO I=1,NP
            WRITE(NFIL,*)P(1:2,I),ORB(I,IORB)
          ENDDO
          CALL FILEHANDLER$CLOSE('HOOK')
          CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
        ENDDO
      ELSE
        CALL ERROR$MSG('INVALID ID')
        CALL ERROR$MSG('MUST BE "1D", "2D", OR "3D"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO_GRIDPLOT_TAILED')
      END IF
!
!     ==========================================================================
!     == CLOSE DOWN                                                           ==
!     ==========================================================================
      DEALLOCATE(P)
      DEALLOCATE(ORB)
                                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_LMTO_GRIDPLOT_TAILEDINNER(NAT,R0,IAT1,NP,P,NORB,ORB)
!     **************************************************************************
!     ** MAPS THE SCREENED ENVELOPE FUNCTION AT ATOM IAT1 ONTO THE GRID P     **
!     ** THE SCREENED ENVELOPE FUNCTION IS DESCRIBED BY THE ONSITE            **
!     ** CONTRIBUTION WITH MATCHING EXPONENTIAL TAILS                         **
!     **                                                                      **
!     ** ATTENTION: SCREENING IS DETERMINED BY ISCATT                         **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2010 ************
      USE LMTO_MODULE, ONLY : SBAR,POTPAR,LNX,LOX,ISPECIES,SBARLI1,K2
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT       ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT) ! ATOMIC POSITIONS
      INTEGER(4),INTENT(IN) :: IAT1      ! INDEX OF CENTRAL ATOM
      INTEGER(4),INTENT(IN) :: NP        ! #(GRID POINTS)
      REAL(8)   ,INTENT(IN) :: P(3,NP)   ! GRIDPOINT POSITIONS
      INTEGER(4),INTENT(IN) :: NORB      !#(ORBITALS ON CENTRAL ATOM)
      REAL(8)   ,INTENT(OUT):: ORB(NP,NORB)  ! SCREENED ENVELOPE FUNCTIONS
      INTEGER(4)            :: NNS
      REAL(8)               :: SBARMAT(NORB,NORB)
      INTEGER(4)            :: ISP  ! SPECIES ID OF ATOM IAT1
      INTEGER(4)            :: LMNXS,LMNX
      INTEGER(4)            :: LX
      INTEGER(4)            :: LMX
      REAL(8)   ,ALLOCATABLE:: YLM(:)
      REAL(8)   ,ALLOCATABLE:: SBARLOC(:,:)
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)               :: RX
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      REAL(8)               :: DR(3),DRLEN
      REAL(8)               :: VAL
      INTEGER(4)            :: NN,L,IM,IP,LN,LM,LMN,I0
!     **************************************************************************
      ISP=ISPECIES(IAT1)
      GID=POTPAR(ISP)%TAILED%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      RX=R(NR)
      DEALLOCATE(R)
      LX=MAXVAL(LOX(:LNX(ISP),ISP))
      LMX=(LX+1)**2
!
!     ==========================================================================
!     == SIZE OF SBARLOC                                                      **
!     ==========================================================================
      LMNXS=0
      LMNX=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        LMNX=LMNX+2*L+1
        IF(POTPAR(ISP)%LNSCATT(LN).EQ.LN) LMNXS=LMNXS+2*L+1
      ENDDO
      ALLOCATE(SBARLOC(LMNX,LMNXS))  !C
!
!     ==========================================================================
!     == COLLECT SCREENED STRUCTURE CONSTANTS FOR ATOM IAT1 ====================
!     ==========================================================================
      NNS=SIZE(SBAR)  !SBAR STRUCCONS. (GLOBAL)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT1) CYCLE
        IF(SBAR(NN)%IAT2.NE.IAT1) CYCLE
        IF(MAXVAL(ABS(SBAR(NN)%IT(:))).NE.0) CYCLE
        IF(SBAR(NN)%N1.NE.LMNXS.OR.SBAR(NN)%N2.NE.LMNXS) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF ORBITALS ')
          CALL ERROR$STOP('LMTO_LMTO_GRIDPLOT_TAILEDINNER')
        END IF
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I0=SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1 
            SBARLOC(LMN+IM,:)=SBAR(NN)%MAT(I0+IM,:)  !C
          ENDDO
          LMN=LMN+2*L+1
        ENDDO
        EXIT
      ENDDO
!
!     ==========================================================================
!     == CALCULATE ORBITALS                                                   ==
!     ==========================================================================
      ORB(:,:)=0.D0
      ALLOCATE(YLM(LMX))
      DO IP=1,NP
!       == DR IS THE DISTANCE FROM THE SECOND ATOM =============================
        DR(:)=P(:,IP)-R0(:,IAT1)
        DRLEN=SQRT(SUM(DR**2))
        IF(DRLEN.GT.RX) CYCLE
        CALL SPHERICAL$YLM(LMX,DR,YLM)
        LMN=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          CALL RADIAL$VALUE(GID,NR,POTPAR(ISP)%TAILED%AEF(:,LN),DRLEN,VAL)
          DO IM=1,2*L+1
            LMN=LMN+1
            LM=L**2+IM
            ORB(IP,LMN)=ORB(IP,LMN)+VAL*YLM(LM)
          ENDDO
        ENDDO
        DO LN=LNX(ISP)+1,POTPAR(ISP)%TAILED%LNX
          CALL RADIAL$VALUE(GID,NR,POTPAR(ISP)%TAILED%AEF(:,LN),DRLEN,VAL)
          L=POTPAR(ISP)%TAILED%LOX(LN)
          I0=SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            LM=L**2+IM
            ORB(IP,:)=ORB(IP,:)+VAL*YLM(LM)*SBARLOC(:,I0+IM)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GRIDPLOT_UNTAILED(IAT0)
!     **************************************************************************
!     **  WRITES THE ENVELOPE FUNCTIONS CENTERED AT ATOM IAT0 TO FILE         **
!     **                                                                      **
!     **  SET N1,N2,N3 EQUAL TO THE NUMBER OF GRIDPOINTS IN EACH DIRECTION    **
!     **                                                                      **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2009 ************
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      USE LMTO_MODULE, ONLY : ISPECIES,SBAR,POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT0
      INTEGER(4),PARAMETER  :: N1=40,N2=40,N3=40 !GRID (1D?)
      INTEGER(4),PARAMETER  :: NRAD=200
      INTEGER(4),PARAMETER  :: NDIRX=100
      LOGICAL(4) ,PARAMETER :: T2D=.TRUE.
      LOGICAL(4) ,PARAMETER :: T3D=.FALSE.
      REAL(8)   ,PARAMETER  :: RANGE=8.D0
      REAL(8)               :: DIR(3,NDIRX)
      REAL(8)               :: ORIGIN(3)
      REAL(8)               :: TVEC(3,3)    !BOX
      REAL(8)               :: RBAS(3,3)
      INTEGER(4)            :: NAT
      INTEGER(4)            :: LM1
      INTEGER(4)            :: LMX         !#YLM
      INTEGER(4)            :: NORB
      REAL(8)   ,ALLOCATABLE:: R0(:,:)      !(3,NAT)
      REAL(8)   ,ALLOCATABLE:: RAD(:)      !(NAT)
      INTEGER(4)            :: NNS
      INTEGER(4)            :: L
      INTEGER(4)            :: LNX
      INTEGER(4)            :: NDIR
      REAL(8)               :: AEZ
      REAL(8)   ,ALLOCATABLE:: RCLUSTER(:,:) ! POSITIONS OF NEIGBORS
      REAL(8)   ,ALLOCATABLE:: ZCLUSTER(:)   ! ATOMIC NUMBER OF NEIGBORS
      REAL(8)   ,ALLOCATABLE:: ORB(:,:)
      REAL(8)   ,ALLOCATABLE:: ORB1(:,:) 
      REAL(8)   ,ALLOCATABLE:: ENV(:,:) 
      REAL(8)   ,ALLOCATABLE:: ENV1(:,:) 
      REAL(8)   ,ALLOCATABLE:: ORBG(:,:)  !GAUSS 
      INTEGER(4),ALLOCATABLE:: LOX(:)
      INTEGER(4),ALLOCATABLE:: ISCATT(:)
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: NATCLUSTER !#(ATOMS ON THE CLUSTER)
      CHARACTER(64)         :: FILE
      CHARACTER(15)         :: STRING         !CONTAINS IATO
      REAL(8)               :: X1D(NRAD)
      INTEGER(4)            :: NP
      REAL(8)  ,ALLOCATABLE :: P(:,:)
      LOGICAL(4),PARAMETER  :: TGAUSS=.TRUE.
      CHARACTER(2),PARAMETER :: ID='3D'
      INTEGER(4)            :: NN,IAT,IAT2,ISP,LN,I,IORB,IDIR
!     **************************************************************************
                                      CALL TRACE$PUSH('LMTO_GRIDPLOT_UNTAILED')
      IF(ID.NE.'1D'.AND.ID.NE.'2D'.AND.ID.NE.'3D') THEN
        CALL ERROR$MSG('INVALID ID')
        CALL ERROR$MSG('MUST BE "1D", "2D", OR "3D"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO_GRIDPLOT_TAILED')
      END IF
!
!     ==========================================================================
!     == COLLECT INFORMATION                                                  ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      NAT=SIZE(ISPECIES)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      ALLOCATE(RAD(NAT))
      LMX=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        RAD(IAT)=POTPAR(ISP)%RAD
!
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(ISCATT(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
        NORB=0
        DO LN=1,LNX
         
          NORB=NORB+2*L+1
          IF(ISCATT(LN).NE.0) CYCLE
          L=LOX(LN)
          IF(IAT.EQ.IAT0)LMX=MAX(LMX,(L+1)**2)
        ENDDO
        DEALLOCATE(LOX)
        DEALLOCATE(ISCATT)
        CALL SETUP$ISELECT(0)
      ENDDO
!
!     ==========================================================================
!     == DEFINE ATOMS ON THE CLUSTER OF NEIGHBORS                             ==
!     ==========================================================================
      NATCLUSTER=0
      NNS=SIZE(SBAR)  !SBAR STRUCCONS. (GLOBAL)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.EQ.IAT0) NATCLUSTER=NATCLUSTER+1
      ENDDO
      ALLOCATE(ZCLUSTER(NATCLUSTER))
      ALLOCATE(RCLUSTER(3,NATCLUSTER))
      IAT=0
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT0) CYCLE
        IAT=IAT+1
        IAT2=SBAR(NN)%IAT2
        RCLUSTER(:,IAT)=R0(:,IAT2) &
     &                 +RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
     &                 +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
     &                 +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
        ISP=ISPECIES(IAT2)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL SETUP$ISELECT(0)
        ZCLUSTER(IAT)=AEZ
      ENDDO
!
!     ==========================================================================
!     == DEFINE GRID POINTS                                                   ==
!     ==========================================================================
      IF(ID.EQ.'3D') THEN
        NP=N1*N2*N3
        ALLOCATE(P(3,NP))
        CALL LMTO_GRIDORB_CUBEGRID(R0(:,IAT0),RANGE,N1,N2,N3,ORIGIN,TVEC,P)
      ELSE IF(ID.EQ.'1D') THEN
        NDIR=0
        DO IAT=1,NATCLUSTER
          NDIR=NDIR+1
          IF(NDIR.GT.NDIRX) THEN
            CALL ERROR$MSG('#(NEIGHBORS EXCEEDS MAXIMUM')
            CALL ERROR$STOP('LMTO_GRIDPLOT_TAILED')
          END IF
          DIR(:,NDIR)=RCLUSTER(:,IAT)-R0(:,IAT0)
          IF(SQRT(SUM(DIR(:,NDIR)**2)).LT.1.D-3) THEN
            NDIR=NDIR-1
            CYCLE   
          END IF
        ENDDO
        NP=NRAD*NDIR
        ALLOCATE(P(3,NP))
        CALL LMTO_GRIDORB_STARGRID(R0(:,IAT0),RANGE,NDIR,DIR,NRAD,X1D,P)
      END IF
!
!     ==========================================================================
!     == DETERMINE ENVELOPE FUNCTION AT THE GRID POINTS                       ==
!     ==========================================================================
      NORB=LMX
      ALLOCATE(ORB(NP,NORB))
      ALLOCATE(ENV(NP,NORB))
      ALLOCATE(ORB1(NP,NORB))
      ALLOCATE(ENV1(NP,NORB))
      ALLOCATE(ORBG(NP,NORB))
      CALL LMTO_GRIDENVELOPE(RBAS,NAT,R0,IAT0,LMX,NP,P,ENV,ENV1)
      CALL LMTO_GRIDAUGMENT(RBAS,NAT,R0,IAT0,LMX,NP,P,ORB1,ENV1)
      ORB=ENV+ORB1-ENV1
      ORBG=0.D0
      IF(TGAUSS)CALL LMTO_GRIDGAUSS(RBAS,NAT,R0,IAT0,LMX,NP,P,ORBG)
!
!     ==========================================================================
!     == WRITE ORBITALS TO FILE                                               ==
!     ==========================================================================
      IF(ID.EQ.'3D') THEN
        IF(TGAUSS)ORB=ORBG
        DO IORB=1,NORB
          FILE='NTB3D'
          WRITE(STRING,*)IAT0 
          FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))
          WRITE(STRING,*)IORB 
          FILE=TRIM(ADJUSTL(FILE))//'_'//TRIM(ADJUSTL(STRING))//'.CUB'
          CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
          CALL FILEHANDLER$UNIT('HOOK',NFIL)
!
LM1=0
CALL ERROR$MSG('CODING ERROR: LM1 NOT DEFINED')
CALL ERROR$STOP('LMTO_GRIDPLOT_UNTAILED')
          CALL LMTO_WRITECUBEFILE(NFIL,NATCLUSTER,ZCLUSTER,RCLUSTER &
     &                         ,ORIGIN,TVEC,N1,N2,N3,ORB(:,LM1))
          CALL FILEHANDLER$CLOSE('HOOK')
          CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
        ENDDO
      ELSE IF(ID.EQ.'1D') THEN
        DO LM1=1,LMX
          FILE='NTB1D'
          WRITE(STRING,*)IAT0 
          FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))
          WRITE(STRING,*)IORB 
          FILE=TRIM(ADJUSTL(FILE))//'_'//TRIM(ADJUSTL(STRING))//'.DAT'
          CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
          CALL FILEHANDLER$UNIT('HOOK',NFIL)
          DO I=1,NRAD
            WRITE(NFIL,*)X1D(I),(ORB(NRAD*(IDIR-1)+I,LM1) &
      &                         ,ORB1(NRAD*(IDIR-1)+I,LM1) &
      &                         ,ENV(NRAD*(IDIR-1)+I,LM1) &
      &                         ,ENV1(NRAD*(IDIR-1)+I,LM1) &
      &                         ,ORBG(NRAD*(IDIR-1)+I,LM1),IDIR=1,NDIR)
          ENDDO
          CALL FILEHANDLER$CLOSE('HOOK')
          CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
        ENDDO
      ELSE IF(ID.EQ.'2D') THEN
        DO LM1=1,LMX
          FILE='NTB2D'
          WRITE(STRING,*)IAT0 
          FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))
          WRITE(STRING,*)IORB 
          FILE=TRIM(ADJUSTL(FILE))//'_'//TRIM(ADJUSTL(STRING))//'.DAT'
          CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
          CALL FILEHANDLER$UNIT('HOOK',NFIL)
          DO I=1,NP
            WRITE(NFIL,*)P(1:2,I),ORB(I,LM1)
          ENDDO
          CALL FILEHANDLER$CLOSE('HOOK')
          CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
        ENDDO
      END IF
!
!     ==========================================================================
!     == CLOSE DOWN                                                           ==
!     ==========================================================================
      DEALLOCATE(P)
      DEALLOCATE(ORB)
      DEALLOCATE(ORB1)
      DEALLOCATE(ORBG)
      DEALLOCATE(ENV)
      DEALLOCATE(ENV1)
                                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PLOTLOCORB(IAT0)
!     **************************************************************************
!     **  WRITES THE ENVELOPE FUNCTIONS CENTERED AT ATOM IAT0 TO FILE         **
!     **                                                                      **
!     **  SET N1,N2,N3 EQUAL TO THE NUMBER OF GRIDPOINTS IN EACH DIRECTION    **
!     **                                                                      **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2009 ************
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      USE LMTO_MODULE, ONLY : K2,SBAR,POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT0
      INTEGER(4),PARAMETER  :: N1=40,N2=40,N3=40 !GRID (1D?)
      INTEGER(4),PARAMETER  :: N1D=1000
      LOGICAL(4) ,PARAMETER :: T2D=.TRUE.
      LOGICAL(4) ,PARAMETER :: T3D=.FALSE.
      REAL(8)               :: ORIGIN(3)
      REAL(8)               :: TVEC(3,3)    !BOX
      REAL(8)               :: TLITTLE(3,3)
      REAL(8)               :: RBAS(3,3)
      INTEGER(4)            :: NAT
      INTEGER(4)            :: LM1
      INTEGER(4)            :: LM1X         !#YLM
      REAL(8)   ,ALLOCATABLE:: R0(:,:)      !(3,NAT)
      INTEGER(4),ALLOCATABLE:: ISPECIES1(:)  !(NAT)
      REAL(8)   ,ALLOCATABLE:: RAD(:)      !(NAT)
      REAL(8)               :: XI(3)
      REAL(8)               :: DR(3)
      INTEGER(4)            :: NNS
      INTEGER(4)            :: NN,IAT,IAT2,ISP,I1,I2,I3,LN,I,J
      INTEGER(4)            :: L
      INTEGER(4)            :: LNX
      REAL(8)               :: AEZ
      INTEGER(4),PARAMETER  :: LMXX=36
      REAL(8)   ,ALLOCATABLE:: RCLUSTER(:,:) ! POSITIONS OF NEIGBORS
      REAL(8)   ,ALLOCATABLE:: ZCLUSTER(:)   ! ATOMIC NUMBER OF NEIGBORS
      REAL(8)   ,ALLOCATABLE:: ORB(:,:)
      REAL(8)   ,ALLOCATABLE:: ORB1(:,:) 
      REAL(8)   ,ALLOCATABLE:: ENV(:,:) 
      REAL(8)   ,ALLOCATABLE:: ENV1(:,:) 
      REAL(8)   ,ALLOCATABLE:: ORBG(:,:)  !GAUSS 
      REAL(8)   ,ALLOCATABLE:: QBARVEC(:,:)
      INTEGER(4),ALLOCATABLE:: LOX(:)
      INTEGER(4),ALLOCATABLE:: ISCATT(:)
      INTEGER(4)            :: NFIL
      INTEGER(4)            :: NATCLUSTER !#(ATOMS ON THE CLUSTER)
      CHARACTER(64)         :: FILE
      CHARACTER(15)         :: STRING         !CONTAINS IATO
      REAL(8)               :: X1D(N1D)
      INTEGER(4)            :: NP,IP
      REAL(8)  ,ALLOCATABLE :: P(:,:)
      LOGICAL(4),PARAMETER  :: TGAUSS=.TRUE.
!     **************************************************************************
                                              CALL TRACE$PUSH('LMTO_PLOTLOCORB')
      NNS=SIZE(SBAR)  !SBAR STRUCCONS. (GLOBAL)
!
!     ==========================================================================
!     == COLLECT INFORMATION                                                  ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      ALLOCATE(ISPECIES1(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES1)
      ALLOCATE(RAD(NAT))
      ALLOCATE(QBARVEC(LMXX,NAT))
      QBARVEC(:,:)=0.D0
      LM1X=0
      DO IAT=1,NAT
        ISP=ISPECIES1(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        RAD(IAT)=POTPAR(ISP)%RAD
!
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(ISCATT(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
        DO LN=1,LNX
          IF(ISCATT(LN).NE.0) CYCLE
          L=LOX(LN)
          QBARVEC(L**2+1:(L+1)**2,IAT)=POTPAR(ISP)%QBAR(LN)
          IF(IAT.EQ.IAT0)LM1X=MAX(LM1X,(L+1)**2)
        ENDDO
        DEALLOCATE(LOX)
        DEALLOCATE(ISCATT)
        CALL SETUP$ISELECT(0)
      ENDDO
!
!     ==========================================================================
!     == DEFINE ATOMS ON THE CLUSTER OF NEIGHBORS                             ==
!     ==========================================================================
      NATCLUSTER=0
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.EQ.IAT0) NATCLUSTER=NATCLUSTER+1
      ENDDO
      ALLOCATE(ZCLUSTER(NATCLUSTER))
      ALLOCATE(RCLUSTER(3,NATCLUSTER))
      IAT=0
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT0) CYCLE
        IAT=IAT+1
        IAT2=SBAR(NN)%IAT2
        RCLUSTER(:,IAT)=R0(:,IAT2) &
     &                 +RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
     &                 +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
     &                 +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
        ISP=ISPECIES1(IAT2)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL SETUP$ISELECT(0)
        ZCLUSTER(IAT)=AEZ
      ENDDO
!
!     ==========================================================================
!     == DEFINE CENTERED GRID                                                 ==
!     ==========================================================================
      TVEC(:,:)=0.D0
      TVEC(1,1)=2*7.17D0 !2*CAMNO3
      TVEC(2,2)=2*7.17D0
      TVEC(3,3)=2*7.17D0
      IF(N1.EQ.1)TVEC(:,1)=0.D0
      IF(N2.EQ.1)TVEC(:,2)=0.D0
      IF(N3.EQ.1)TVEC(:,3)=0.D0
      TLITTLE=TVEC
      IF(N1.GT.1)TLITTLE(:,1)=TLITTLE(:,1)/REAL(N1-1,KIND=8)
      IF(N2.GT.1)TLITTLE(:,2)=TLITTLE(:,2)/REAL(N2-1,KIND=8)
      IF(N3.GT.1)TLITTLE(:,3)=TLITTLE(:,3)/REAL(N3-1,KIND=8)
      ORIGIN(:)=-0.5D0*(TVEC(:,1)+TVEC(:,2)+TVEC(:,3)) !ALSO 2D
      ORIGIN(:)=ORIGIN(:)+R0(:,IAT0)
!
!     ==========================================================================
!     == DEFINE GRID POINTS                                                   ==
!     ==========================================================================
      IF(T3D) THEN
      NP=N1*N2*N3
      ALLOCATE(P(3,NP))
      IP=0
      DO I3=1,N3
        XI(3)=REAL(I3-1,KIND=8)
        DO I2=1,N2
          XI(2)=REAL(I2-1,KIND=8)
          DO I1=1,N1
            XI(1)=REAL(I1-1,KIND=8)
            IP=IP+1
            P(:,IP)=ORIGIN(:)+MATMUL(TLITTLE,XI)
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE ENVELOPE FUNCTION AT THE GRID POINTS                       ==
!     ==========================================================================
      ALLOCATE(ORB(NP,LM1X))
      ALLOCATE(ENV(NP,LM1X))
      ALLOCATE(ORB1(NP,LM1X))
      ALLOCATE(ENV1(NP,LM1X))
      ALLOCATE(ORBG(NP,LM1X))
      CALL LMTO_GRIDENVELOPE(RBAS,NAT,R0,IAT0,LM1X,NP,P,ENV,ENV1)
      CALL LMTO_GRIDAUGMENT(RBAS,NAT,R0,IAT0,LM1X,NP,P,ORB1,ENV1)
      ORB=ENV+ORB1-ENV1
      ORBG=0.D0
      IF(TGAUSS)CALL LMTO_GRIDGAUSS(RBAS,NAT,R0,IAT0,LM1X,NP,P,ORBG)
!
!     ==========================================================================
!     == WRITE CUBE FILES                                                     ==
!     ==========================================================================
      IF(TGAUSS)ORB=ORBG
      DO LM1=1,LM1X
        FILE='NTB'
        WRITE(STRING,*)IAT0 
        FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))
        CALL SPHERICAL$YLMNAME(LM1,STRING)
        FILE=TRIM(ADJUSTL(FILE))//TRIM(ADJUSTL(STRING))//'.CUB'
        CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
        CALL FILEHANDLER$UNIT('HOOK',NFIL)
!
        CALL LMTO_WRITECUBEFILE(NFIL,NATCLUSTER,ZCLUSTER,RCLUSTER &
     &                         ,ORIGIN,TVEC,N1,N2,N3,ORB(:,LM1))
        CALL FILEHANDLER$CLOSE('HOOK')
        CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
      ENDDO
      DEALLOCATE(P)
      DEALLOCATE(ORB)
      DEALLOCATE(ORB1)
      DEALLOCATE(ORBG)
      DEALLOCATE(ENV)
      DEALLOCATE(ENV1)
      END IF
!
!     ==========================================================================
!     == DEFINE 1-DIMENSIONAL GRIDS                                           ==
!     == GRIDS POINT TOWARDS NEAREST NEIGHBORS. (OR IN X-DIRECTION)           ==
!     ==========================================================================
      NP=N1D*MAX((NATCLUSTER-1),1)
      ALLOCATE(P(3,NP))
      DO I=1,N1D
        X1D(I)=(-1.D0+2.D0*REAL(I-1)/REAL(N1D-1))*50.D0
      ENDDO
      IP=0
      DO IAT=1,NATCLUSTER
        DR(:)=RCLUSTER(:,IAT)-R0(:,IAT0)
        IF(SQRT(SUM(DR**2)).LT.1.D-3) CYCLE   !AVID ONSITE
        DR(:)=DR(:)/SQRT(SUM(DR**2))
        DO I=1,N1D
          IP=IP+1
          P(:,IP)=R0(:,IAT0)+DR(:)*X1D(I)
        ENDDO
      ENDDO
      IF(NATCLUSTER.EQ.1) THEN
        P(:,:)=0.D0
        P(1,:)=X1D
      END IF
!
!     ==========================================================================
!     == DETERMINE ENVELOPE FUNCTION AT THE GRID POINTS                       ==
!     ==========================================================================
      ALLOCATE(ENV(NP,LM1X))
      ALLOCATE(ORB1(NP,LM1X))
      ALLOCATE(ORB(NP,LM1X))
      ALLOCATE(ENV1(NP,LM1X))
      ALLOCATE(ORBG(NP,LM1X))
      CALL LMTO_GRIDENVELOPE(RBAS,NAT,R0,IAT0,LM1X,NP,P,ENV,ENV1)
      CALL LMTO_GRIDAUGMENT(RBAS,NAT,R0,IAT0,LM1X,NP,P,ORB1,ENV1)
      ORB=ENV+ORB1-ENV1
      ORBG=0.D0
      IF(TGAUSS)CALL LMTO_GRIDGAUSS(RBAS,NAT,R0,IAT0,LM1X,NP,P,ORBG)
!      CALL LMTO_GRIDTAILS(NAT,R0,IAT0,NP,P,LM1X,ORBG)
!
!     ==========================================================================
!     == WRITE ORBITALS TO FILE                                               ==
!     ==========================================================================
      DO LM1=1,LM1X
        FILE='NTB'
        WRITE(STRING,*)IAT0 
        FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))
        CALL SPHERICAL$YLMNAME(LM1,STRING)
        FILE=TRIM(ADJUSTL(FILE))//TRIM(ADJUSTL(STRING))//'.DAT'
        CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
        CALL FILEHANDLER$UNIT('HOOK',NFIL)
!
        IF(NATCLUSTER.EQ.1) THEN
          DO I=1,N1D
            WRITE(NFIL,*)X1D(I),ORB(I,LM1),ORB1(I,LM1),ENV(I,LM1) &
      &                      ,ENV1(I,LM1),ORBG(I,LM1)
          ENDDO
        ELSE
          DO I=1,N1D
            WRITE(NFIL,*)X1D(I),(ORB(N1D*(IAT-1)+I,LM1) &
      &                         ,ORB1(N1D*(IAT-1)+I,LM1) &
      &                         ,ENV(N1D*(IAT-1)+I,LM1) &
      &                         ,ENV1(N1D*(IAT-1)+I,LM1) &
      &                         ,ORBG(N1D*(IAT-1)+I,LM1),IAT=1,MIN(NATCLUSTER-1,5))
          ENDDO
        END IF
        CALL FILEHANDLER$CLOSE('HOOK')
        CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
      ENDDO
      DEALLOCATE(P)
      DEALLOCATE(ORB)
      DEALLOCATE(ORBG)
      DEALLOCATE(ORB1)
      DEALLOCATE(ENV)
      DEALLOCATE(ENV1)

!ATTENTIONATTENTION START FUDGING
IF(T2D) THEN
!
!     ==========================================================================
!     == DEFINE 2-DIMENSIONAL GRID                                           ==
!     ==========================================================================
      NP=N1*N2
      ALLOCATE(P(3,NP))
      TVEC(:,1)=(/10.D0,0.D0,0.D0/)
      TVEC(:,2)=(/0.D0,10.D0,0.D0/)
      ORIGIN(:)=(/1.D0,0.D0,0.D0/)-0.5D0*(TVEC(:,1)+TVEC(:,2))
      TLITTLE(:,1)=TVEC(:,1)/REAL(N1-1,KIND=8)
      TLITTLE(:,2)=TVEC(:,2)/REAL(N2-1,KIND=8)
      IP=0
      DO I=1,N1
        DO J=1,N2
          IP=IP+1
          P(:,IP)=ORIGIN+TLITTLE(:,1)*REAL(I-1,KIND=8) &
     &                  +TLITTLE(:,2)*REAL(J-1,KIND=8) 
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE ENVELOPE FUNCTION AT THE GRID POINTS                       ==
!     ==========================================================================
      ALLOCATE(ENV(NP,LM1X))
      ALLOCATE(ORB1(NP,LM1X))
      ALLOCATE(ORB(NP,LM1X))
      ALLOCATE(ENV1(NP,LM1X))
      ALLOCATE(ORBG(NP,LM1X))
      CALL LMTO_GRIDENVELOPE(RBAS,NAT,R0,IAT0,LM1X,NP,P,ENV,ENV1)
      CALL LMTO_GRIDAUGMENT(RBAS,NAT,R0,IAT0,LM1X,NP,P,ORB1,ENV1)
      ORB=ENV+ORB1-ENV1
      ORBG=0.D0
      IF(TGAUSS)CALL LMTO_GRIDGAUSS(RBAS,NAT,R0,IAT0,LM1X,NP,P,ORBG)
!
!     ==========================================================================
!     == WRITE ORBITALS TO FILE                                               ==
!     ==========================================================================
      IF(TGAUSS)ORB=ORBG
      DO LM1=1,LM1X
        FILE='NTB2D'
        WRITE(STRING,*)IAT0 
        FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))
        CALL SPHERICAL$YLMNAME(LM1,STRING)
        FILE=TRIM(ADJUSTL(FILE))//TRIM(ADJUSTL(STRING))//'.DAT'
        CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
        CALL FILEHANDLER$UNIT('HOOK',NFIL)
!
        DO I=1,NP
          WRITE(NFIL,*)P(1:2,I),ORB(I,LM1)
        ENDDO
        CALL FILEHANDLER$CLOSE('HOOK')
        CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
      ENDDO
      DEALLOCATE(P)
      DEALLOCATE(ORB)
      DEALLOCATE(ORB1)
      DEALLOCATE(ORBG)
      DEALLOCATE(ENV)
      DEALLOCATE(ENV1)
      CALL ERROR$MSG('PROGRAM STOPS AFTER WRITING DATA FOR 2-D PLOT')
      CALL ERROR$MSG('OF LOCAL ORBITALS')
      CALL ERROR$MSG('AVOID THIS OPTION BY SETTING T2D=.FALSE.')
      CALL ERROR$STOP('LMTO_PLOTLOCORB')
END IF

!CALL SPECIALFUNCTION$TEST()
CALL TEST_LMTO$STRUCTURECONSTANTS()
                                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GRIDORB_CUBEGRID(CENTER,RAD,N1,N2,N3,ORIGIN,TVEC,P)
!     **************************************************************************
!     **  DEFINES A GRID CENTERED AT "CENTER" ENCLOSING A SPHERE WITH RADIUS  **
!     **  RAD WITH N1,N2,N3 GRID POINTS IN EACH DIRECTION.                    **
!     **  THE RETURNED DRAWING BOX IS DEFINED BY ORIGIN AND SIDE VECTORS TVEC **
!     **  THE GRID POINTS ARE RETURNED IN P                                   **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2009 ************
      IMPLICIT NONE
      REAL(8)    ,INTENT(IN) :: CENTER(3)
      REAL(8)    ,INTENT(IN) :: RAD 
      INTEGER(4) ,INTENT(IN) :: N1,N2,N3
      REAL(8)    ,INTENT(OUT):: ORIGIN(3)
      REAL(8)    ,INTENT(OUT):: TVEC(3,3)
      REAL(8)    ,INTENT(OUT):: P(3,N1*N2*N3)
      REAL(8)                :: XI(3)
      REAL(8)                :: TLITTLE(3,3)
      INTEGER(4)             :: I1,I2,I3,IP
!     **************************************************************************
!
!     ==========================================================================
!     == DEFINE DRAWING BOX                                                   ==
!     ==========================================================================
      TVEC(:,:)=0.D0
      TVEC(1,1)=2*RAD
      TVEC(2,2)=2*RAD
      TVEC(3,3)=2*RAD
      IF(N1.EQ.1)TVEC(:,1)=0.D0
      IF(N2.EQ.1)TVEC(:,2)=0.D0
      IF(N3.EQ.1)TVEC(:,3)=0.D0
      ORIGIN(:)=CENTER(:)-0.5D0*(TVEC(:,1)+TVEC(:,2)+TVEC(:,3)) !ALSO 2D
!
!     ==========================================================================
!     == DEFINE GRID-STEPS                                                    ==
!     ==========================================================================
      TLITTLE=TVEC
      IF(N1.GT.1)TLITTLE(:,1)=TLITTLE(:,1)/REAL(N1-1,KIND=8)
      IF(N2.GT.1)TLITTLE(:,2)=TLITTLE(:,2)/REAL(N2-1,KIND=8)
      IF(N3.GT.1)TLITTLE(:,3)=TLITTLE(:,3)/REAL(N3-1,KIND=8)
!
!     ==========================================================================
!     == DEFINE GRID POINTS                                                   ==
!     ==========================================================================
      IP=0
      DO I3=1,N3
        XI(3)=REAL(I3-1,KIND=8)
        DO I2=1,N2
          XI(2)=REAL(I2-1,KIND=8)
          DO I1=1,N1
            XI(1)=REAL(I1-1,KIND=8)
            IP=IP+1
            P(:,IP)=ORIGIN(:)+MATMUL(TLITTLE,XI)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GRIDORB_STARGRID(CENTER,RAD,NDIR,DIR,NR,X1D,P)
!     **************************************************************************
!     **  THE GRID POINTS ARE RETURNED IN P                                   **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2009 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: CENTER(3)
      REAL(8)   ,INTENT(IN) :: RAD
      REAL(8)   ,INTENT(IN) :: DIR(3,NDIR)
      INTEGER(4),INTENT(IN) :: NDIR
      REAL(8)   ,INTENT(OUT):: X1D(NR)
      REAL(8)   ,INTENT(OUT):: P(3,NR*NDIR)
      REAL(8)               :: DR(3)
      INTEGER(4)            :: I,IDIR,IP
!     **************************************************************************
      DO I=1,NR
        X1D(I)=RAD*(-1.D0+2.D0*REAL(I-1,KIND=8)/REAL(NR-1,KIND=8))
      ENDDO
      IP=0
      DO IDIR=1,NDIR
        DR(:)=DIR(:,IDIR)
        IF(SQRT(SUM(DR**2)).GT.1.D-5)  THEN
          DR(:)=DR(:)/SQRT(SUM(DR**2))
        END IF
        DO I=1,NR
          IP=IP+1
          P(:,IP)=CENTER(:)+DR(:)*X1D(I)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PLOTNTBO(TYPE,IATORB,LMNORB)
!     **************************************************************************
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2012 ************
      USE LMTO_MODULE, ONLY : ISPECIES,NSP,LNX,LOX,SBAR,K2,POTPAR,SBARLI1
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: TYPE    !
      INTEGER(4)  ,INTENT(IN) :: IATORB
      INTEGER(4)  ,INTENT(IN) :: LMNORB
      LOGICAL(4)              :: TTAILED
      CHARACTER(8)            :: GRIDTYPE
      REAL(8)                 :: RBAS(3,3)
      INTEGER(4)              :: NAT
      REAL(8)   ,ALLOCATABLE  :: R0(:,:)
      INTEGER(4)              :: NNS
      INTEGER(4)              :: NATCLUSTER
      REAL(8)   ,ALLOCATABLE  :: ZCLUSTER(:)
      REAL(8)   ,ALLOCATABLE  :: RCLUSTER(:,:)
      INTEGER(4)              :: NN,ISP,IAT,IAT2
      REAL(8)                 :: CENTER(3)
      REAL(8)   ,PARAMETER    :: RADIUS=10.D0
      REAL(8)                 :: ORIGIN(3),TVEC(3,3) !VIEWBOX
      INTEGER(4),PARAMETER    :: N1=40,N2=40,N3=40
      INTEGER(4),PARAMETER    :: NR=200,NDIR=13
      REAL(8)                 :: DIR(3,NDIR)
      REAL(8)   ,ALLOCATABLE  :: X1D(:)
      INTEGER(4)              :: NP
      REAL(8)   ,ALLOCATABLE  :: P(:,:)
      REAL(8)   ,ALLOCATABLE  :: ORB(:)
      CHARACTER(64)           :: STRING
      CHARACTER(64)           :: FILE
      CHARACTER(8)            :: EXT
      INTEGER(4)              :: NFIL
      INTEGER(4)              :: I,J
!     **************************************************************************
                                      CALL TRACE$PUSH('LMTO_PLOTNTBO')
      IF(INDEX(TYPE,'TAILED').NE.0) THEN
        TTAILED=.TRUE.
      ELSE IF(INDEX(TYPE,'FULL').NE.0) THEN
        TTAILED=.FALSE.
      ELSE 
        CALL ERROR$MSG('ORBITAL TYPE NOT RECOGNIZED')
        CALL ERROR$MSG('MUST BE "FULL" OR "TAILED"')
        CALL ERROR$STOP('LMTO_PLOTNTBO')
      END IF

      IF(INDEX(TYPE,'CUBE').NE.0) THEN
        GRIDTYPE='CUBE'
      ELSE IF(INDEX(TYPE,'STAR').NE.0) THEN
        GRIDTYPE='STAR'
      ELSE 
        CALL ERROR$MSG('GRID-TYPE NOT RECOGNIZED')
        CALL ERROR$MSG('MUST BE "CUBE" OR "STAR"')
        CALL ERROR$STOP('LMTO_PLOTNTBO')
      END IF
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      NAT=SIZE(ISPECIES)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     == DEFINE ATOMS ON THE CLUSTER OF NEIGHBORS                             ==
!     ==========================================================================
      NNS=SIZE(SBAR)
      NATCLUSTER=0
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.EQ.IATORB) NATCLUSTER=NATCLUSTER+1
      ENDDO
      ALLOCATE(ZCLUSTER(NATCLUSTER))
      ALLOCATE(RCLUSTER(3,NATCLUSTER))
      IAT=0
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IATORB) CYCLE
        IAT=IAT+1
        IAT2=SBAR(NN)%IAT2
        RCLUSTER(:,IAT)=R0(:,IAT2) &
     &                 +RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
     &                 +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
     &                 +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
        ISP=ISPECIES(IAT2)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',ZCLUSTER(IAT))
        CALL SETUP$ISELECT(0)
      ENDDO
!
!     ==========================================================================
!     == DEFINE GRID                                                          ==
!     ==========================================================================
      CENTER=R0(:,IATORB)
      IF(GRIDTYPE.EQ.'CUBE') THEN
        NP=N1*N2*N3
        ALLOCATE(P(3,N1*N2*N3))
        CALL LMTO_GRIDORB_CUBEGRID(CENTER,RADIUS,N1,N2,N3,ORIGIN,TVEC,P)
      ELSE IF(GRIDTYPE.EQ.'STAR') THEN
        NP=NR*NDIR
!       == (100) DIRECTIONS ====================================================
        DIR(:,1)=(/1.D0,0.D0,0.D0/)
        DIR(:,2)=(/0.D0,1.D0,0.D0/)
        DIR(:,3)=(/0.D0,0.D0,1.D0/)
!       == (110) DIRECTIONS ====================================================
        DIR(:,4)=(/0.D0,1.D0,1.D0/)
        DIR(:,5)=(/1.D0,0.D0,1.D0/)
        DIR(:,6)=(/1.D0,1.D0,0.D0/)
        DIR(:,7)=(/0.D0,-1.D0,1.D0/)
        DIR(:,8)=(/1.D0,0.D0,-1.D0/)
        DIR(:,9)=(/1.D0,-1.D0,0.D0/)
!       == (111) DIRECTIONS ====================================================
        DIR(:,10)=(/1.D0,1.D0,1.D0/)
        DIR(:,11)=(/-1.D0,1.D0,1.D0/)
        DIR(:,12)=(/1.D0,-1.D0,1.D0/)
        DIR(:,13)=(/1.D0,1.D0,-1.D0/)
        ALLOCATE(P(3,NR*NDIR))
        ALLOCATE(X1D(NR))
        CALL LMTO_GRIDORB_STARGRID(CENTER,RADIUS,NDIR,DIR,NR,X1D,P)
      ELSE
        CALL ERROR$MSG('GRIDTYPE NOT RECOGNIZED. MAY BE "CUBE" OR "STAR"')
        CALL ERROR$STOP('LMTO_PLOTNTBO')
      END IF
!
!     ==========================================================================
!     == CALCULATE ORBITALE                                                   ==
!     ==========================================================================
      ALLOCATE(ORB(NP))
      IF(TTAILED) THEN
        CALL LMTO_TAILED_NTBOOFR(TYPE,IATORB,LMNORB,NP,P,ORB)
      ELSE
        CALL LMTO_NTBOOFR(TYPE,IATORB,LMNORB,NP,P,ORB)
      END IF
!
!     ==========================================================================
!     == WRITE DATA TO FILE                                                   ==
!     ==========================================================================
      FILE='NTB'
      IF(TTAILED) FILE=TRIM(FILE)//'_TAILED'
      WRITE(STRING,*)IATORB
      FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))
!      CALL SPHERICAL$YLMNAME(LM1,STRING)
      WRITE(STRING,*)LMNORB
      IF(GRIDTYPE.EQ.'CUBE') THEN
        EXT=-'.CUB'
      ELSE
        EXT=-'.DAT'
      END IF
      FILE=TRIM(ADJUSTL(FILE))//'_'//TRIM(ADJUSTL(STRING))//TRIM(EXT)
      CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
      CALL FILEHANDLER$UNIT('HOOK',NFIL)
      IF(GRIDTYPE.EQ.'CUBE') THEN
        CALL LMTO_WRITECUBEFILE(NFIL,NATCLUSTER,ZCLUSTER,RCLUSTER &
     &                       ,ORIGIN,TVEC,N1,N2,N3,ORB(:))
      ELSE IF(GRIDTYPE.EQ.'STAR') THEN
        REWIND(NFIL)
        DO I=1,NR
          WRITE(NFIL,*)X1D(I),(ORB(I+NR*(J-1)),J=1,NDIR)
        ENDDO
      ELSE 
        CALL ERROR$MSG('GRIDTYPE NOT RECOGNIZED. MAY BE "CUBE" OR "STAR"')
        CALL ERROR$STOP('LMTO_PLOTNTBO')
      END IF
      CALL FILEHANDLER$CLOSE('HOOK')
      CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
                                      CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TESTORBITALS()
!     **************************************************************************
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2012 ************
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: IATORB=1
      INTEGER(4),PARAMETER :: LMNORB=1
      REAL(8)   ,PARAMETER :: RAD=2.D0
      INTEGER(4),PARAMETER :: NP=38
      REAL(8)              :: DIR(3,NP)
      CHARACTER(8)         :: TYPE='='
      REAL(8)              :: P(3,NP)
      REAL(8)              :: ORB(NP)
      REAL(8)              :: ORBTAILED(NP)
      INTEGER(4)           :: J,LM
      INTEGER(8),PARAMETER :: LX=4
      INTEGER(8),PARAMETER :: LMX=(LX+1)**2
      REAL(8)              :: YLM(LMX)
      REAL(8)              :: FLM(LMX)
      REAL(8)              :: FLMTAILED(LMX)
      REAL(8)              :: CLM(LMX,LMX)
      REAL(8)              :: CLMINV(LMX,LMX)
      REAL(8)              :: SVAR
!     **************************************************************************
!
!     ==========================================================================
!     == DEFINE STARBURST MESH                                                ==
!     ==========================================================================
!     == (100) DIRECTIONS ======================================================
      DIR(:,1)=(/1.D0,0.D0,0.D0/)
      DIR(:,2)=(/0.D0,1.D0,0.D0/)
      DIR(:,3)=(/0.D0,0.D0,1.D0/)
!     == (110) DIRECTIONS ======================================================
      DIR(:,4)=(/0.D0,1.D0,1.D0/)
      DIR(:,5)=(/1.D0,0.D0,1.D0/)
      DIR(:,6)=(/1.D0,1.D0,0.D0/)
      DIR(:,7)=(/0.D0,-1.D0,1.D0/)
      DIR(:,8)=(/1.D0,0.D0,-1.D0/)
      DIR(:,9)=(/1.D0,-1.D0,0.D0/)
!     == (111) DIRECTIONS ======================================================
      DIR(:,10)=(/1.D0,1.D0,1.D0/)
      DIR(:,11)=(/-1.D0,1.D0,1.D0/)
      DIR(:,12)=(/1.D0,-1.D0,1.D0/)
      DIR(:,13)=(/1.D0,1.D0,-1.D0/)
!     == (120) DIRECTIONS ======================================================
      DIR(:,14)=(/2.D0,1.D0,0.D0/)
      DIR(:,15)=(/0.D0,2.D0,1.D0/)
      DIR(:,16)=(/1.D0,0.D0,2.D0/)
      DIR(:,17)=(/0.D0,1.D0,2.D0/)
      DIR(:,18)=(/1.D0,2.D0,0.D0/)
      DIR(:,19)=(/2.D0,0.D0,1.D0/)
!     ==========================================================================
      DIR(:,NP/2+1:)=-DIR(:,:NP/2)
      DO J=1,NP
        DIR(:,J)=DIR(:,J)/SQRT(SUM(DIR(:,J)**2))
      ENDDO
!
!     ==========================================================================
!     == CALCULATE ORBITAL
!     ==========================================================================
      P(:,:)=DIR(:,:)*RAD
      CALL LMTO_TAILED_NTBOOFR(TYPE,IATORB,LMNORB,NP,P,ORB)
      CALL LMTO_NTBOOFR(TYPE,IATORB,LMNORB,NP,P,ORBTAILED)
      DO J=1,NP
        WRITE(*,FMT='(I5," R=",3F5.2," ORB(FULL)=",F20.10," ORB(TAIL)=",F20.10)') &
     &              J,DIR(:,J),ORB(J),ORBTAILED(J)
      ENDDO
!
!     ==========================================================================
!     == CALCULATE ORBITAL
!     ==========================================================================
      FLM(:)=0.D0
      FLMTAILED(:)=0.D0
      CLM(:,:)=0.D0
      DO J=1,NP
        CALL SPHERICAL$YLM(LMX,P(:,J),YLM)
        WRITE(*,FMT='(" YLM(LMX)=",F15.5)')YLM(LMX)
        FLM(:)=FLM(:)+YLM(:)*ORB(J)
        FLMTAILED(:)=FLMTAILED(:)+YLM(:)*ORBTAILED(J)
        DO LM=1,LMX
          CLM(:,LM)=CLM(:,LM)+YLM(:)*YLM(LM)
        ENDDO
      ENDDO
      SVAR=1.D0/REAL(NP,KIND=8)
      FLM(:)=FLM(:)*SVAR
      FLMTAILED(:)=FLMTAILED(:)*SVAR
      CLM(:,:)=CLM(:,:)*SVAR
      DO LM=1,LMX
        WRITE(*,FMT='("LM=",I3," CLM=",30F8.3)')LM,CLM(:,LM)
      ENDDO
      CALL LIB$INVERTR8(LMX,CLM,CLMINV)
      FLM(:)=MATMUL(CLMINV,FLM)
      FLMTAILED(:)=MATMUL(CLMINV,FLMTAILED)
      DO LM=1,LMX
        WRITE(*,FMT='("LM",I3," FLM(F)=",F15.5," FLM(T)=",F15.5)')LM,FLM(LM),FLMTAILED(LM)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBOOFR(TYPE,IATORB,LMNORB,NP,P,ORB)
!     **************************************************************************
!     **  MAPS THE SCREENED ORBITAL LMNORB AT ATOM IATORB                     **
!     **  ONTO THE GRID P                                                     **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2012 ************
      USE LMTO_MODULE, ONLY : ISPECIES,NSP,LNX,LOX,SBAR,K2,POTPAR,SBARLI1
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: TYPE    !
      INTEGER(4)  ,INTENT(IN) :: IATORB  ! INDEX OF CENTRAL ATOM
      INTEGER(4)  ,INTENT(IN) :: LMNORB  ! ORBITAL INDEX
      REAL(8)     ,INTENT(OUT):: ORB(NP) ! SCREENED ENVELOPE FUNCTIONS
      INTEGER(4)  ,INTENT(IN) :: NP      ! #(GRID POINTS)
      REAL(8)     ,INTENT(IN) :: P(3,NP) ! GRIDPOINT POSITIONS
      REAL(8)               :: RBAS(3,3) ! LATTICE VECTORS
      INTEGER(4)            :: NAT       ! #(ATOMS)
      REAL(8)   ,ALLOCATABLE:: R0(:,:)   ! ATOMIC POSITIONS
      INTEGER(4)            :: LMORB,LNORB,LORB,MORB
      INTEGER(4)            :: IAT2,L2X,LM2X
      INTEGER(4)            :: NNS 
      INTEGER(4)            :: NN,IP,ISP,ISP2,L,LM,IM,LN,LMN  !LOOP INDICES
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: LMXX
      REAL(8)               :: R2(3) ! NEIGHBOR ATOM POSITION
      REAL(8)               :: DR(3) ! DISTANCE VECTOR OF GRID POINT TO ATOM
      REAL(8)               :: DIS   ! DISTANCE OF GRID POINT TO ATOM
      REAL(8)               :: VAL,VALDOT
      LOGICAL(4)            :: TONSITE   
      LOGICAL(4)            :: TSPHERE   ! INSIDE AN AUGMENTATION SPHERE?
      REAL(8)  ,ALLOCATABLE :: QBARVEC(:,:)
      REAL(8)  ,ALLOCATABLE :: CVEC(:)   ! COEFFICIENTS OF BARE HANKEL FUNCTNS.
      REAL(8)  ,ALLOCATABLE :: K0(:)     ! BARE SOLID HANKEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: JBAR(:)   ! SCREENED SOLID BESSEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: J0(:)     ! BARE SOLID BESSEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: YLM(:)    ! SPHERICAL HARMONICS
      REAL(8)  ,ALLOCATABLE :: R(:)      ! RADIAL GRID
      REAL(8)  ,ALLOCATABLE :: AEPHI(:,:)   ! ALL-ELECTRON PARTIAL WAVE
      REAL(8)  ,ALLOCATABLE :: NLPHI(:,:)   ! NODELESS PARTIAL WAVE
      REAL(8)  ,ALLOCATABLE :: NLPHIDOT(:,:)! NODELESS SCATTERING PARTIAL WAVE
      REAL(8)  ,ALLOCATABLE :: K0ARR(:)     ! BARE AUGMENTED HANKEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: DK0ARR(:)    ! 
      REAL(8)  ,ALLOCATABLE :: JBARARR(:,:) ! SCREENED AUGMENTED BESSEL FUNCTION
!     **************************************************************************
PRINT*,'MARKE 1A'
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      NAT=SIZE(ISPECIES)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     == DETERMINE LM-, LN-, L- AND M-INDEX OF THE ORBITAL                    ==
!     ==========================================================================
PRINT*,'MARKE 1B'
      LMORB=0
      ISP=ISPECIES(IATORB)
      LMN=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        DO IM=1,2*L+1
          LMN=LMN+1
          IF(LMN.EQ.LMNORB) THEN
            LORB=LOX(LN,ISP)
            LNORB=LN
            LMORB=SBARLI1(LORB+1,ISP)-1+IM
            EXIT
          END IF
        ENDDO
      ENDDO
      IF(LMORB.EQ.0) THEN
        CALL ERROR$MSG('INDEX ERROR')
        CALL ERROR$MSG('NTBOOFR')
      END IF
PRINT*,'MARKE 1C'
!
!     ==========================================================================
!     == DETERMINE QBAR OF THE SCATTERING CHANNEL                             ==
!     ==========================================================================
      LMXX=MAXVAL(LOX(:,:)+1)**2
      ALLOCATE(QBARVEC(LMXX,NSP))
      QBARVEC(:,:)=0.D0
      DO ISP=1,NSP
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          IF(LN.EQ.POTPAR(ISP)%LNSCATT(LN)) THEN
            LM=SBARLI1(L+1,ISP)-1
            DO IM=1,2*L+1
              LM=LM+1
              QBARVEC(LM,ISP)=POTPAR(ISP)%QBAR(LN)        
            ENDDO
          END IF
        ENDDO       
      ENDDO    
PRINT*,'MARKE 1D'
!
!     ==========================================================================
!     == 
!     ==========================================================================
      NNS=SIZE(SBAR)  !SBAR STRUCCONS. (GLOBAL)
      ALLOCATE(CVEC(LMXX))
      ALLOCATE(K0(LMXX))
      ALLOCATE(J0(LMXX))
      ALLOCATE(JBAR(LMXX))
PRINT*,'MARKE 1E'
!
!     ==========================================================================
!     == LOOP OVER ALL NEIGHBORS AND GRID POINTS ===============================
!     ==========================================================================
      ORB(:)=0.D0
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IATORB) CYCLE
        IAT2=SBAR(NN)%IAT2
        ISP2=ISPECIES(IAT2)
        R2(:)=R0(:,IAT2)+RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
     &                  +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
     &                  +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
!
!       == CVEC=1+QBAR*SBAR ====================================================
        TONSITE=(IAT2.EQ.IATORB).AND.(MAXVAL(ABS(SBAR(NN)%IT(:))).EQ.0)
        LM2X=SBAR(NN)%N2
        CVEC(:LM2X)=QBARVEC(:LM2X,ISP2)*SBAR(NN)%MAT(LMORB,:)    !C
        IF(TONSITE)CVEC(LMORB)=CVEC(LMORB)+1.D0
!
!       == LOOP OVER REAL SPACE GRID ===========================================
        DO IP=1,NP
!         == DR IS THE DISTANCE FROM THE SECOND ATOM ===========================
          DR(:)=P(:,IP)-R2(:)
          TSPHERE=(DOT_PRODUCT(DR,DR).LT.POTPAR(ISP2)%RAD**2)
!
!         == DETERMINE BARE HANKEL AND SCREENED BESSEL FUNCTION AT IAT2 ========
          CALL  LMTO$SOLIDHANKEL(DR,POTPAR(ISP2)%RAD,K2,LM2X,K0(1:LM2X))
          IF(TSPHERE) THEN
            CALL  LMTO$SOLIDBESSEL(DR,K2,LM2X,J0(1:LM2X))
            JBAR(:LM2X)=J0(:LM2X)-K0(:LM2X)*QBARVEC(:LM2X,ISP2)
          END IF
!
!         == ENVELOPE FUNCTION =================================================
!         == |KBAR>=|K0>*(1+QBAR*SBART) ========================================
          ORB(IP)=ORB(IP)+DOT_PRODUCT(K0(1:LM2X),CVEC(1:LM2X))
          IF(TSPHERE) THEN
!           == SUBTRACT MULTI-CENTER EXPANSION IN THE SPHERE ===================
!           == -|JBAR>SBAR=-(|J0>-|K0>QBAR)*SBAR ===============================
            ORB(IP)=ORB(IP)+DOT_PRODUCT(JBAR(1:LM2X),SBAR(NN)%MAT(LMORB,:)) 
            IF(TONSITE) ORB(IP)=ORB(IP)-K0(LMORB)
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(CVEC)
      DEALLOCATE(K0)
      DEALLOCATE(J0)
      DEALLOCATE(JBAR)
PRINT*,'MARKE 1F'
!
!     ==========================================================================
!     ==  NOW DO HE AUGMENTATION                                              ==
!     ==========================================================================
      ALLOCATE(YLM(LMXX))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        L2X=MAXVAL( LOX(:LNX(ISP),ISP) )
        ALLOCATE(R(NR))
        ALLOCATE(AEPHI(NR,LNX(ISP)))
        ALLOCATE(NLPHI(NR,LNX(ISP)))
        ALLOCATE(NLPHIDOT(NR,LNX(ISP)))
        ALLOCATE(K0ARR(NR))
        ALLOCATE(DK0ARR(NR))
        ALLOCATE(JBARARR(NR,L2X+1))
        CALL RADIAL$R(GID,NR,R)
        CALL SETUP$GETR8A('AEPHI',NR*LNX(ISP),AEPHI)
        CALL SETUP$GETR8A('NLPHI',NR*LNX(ISP),NLPHI)
        CALL SETUP$GETR8A('NLPHIDOT',NR*LNX(ISP),NLPHIDOT)
        CALL SETUP$ISELECT(0)
        K0ARR(:)=NLPHI(:,LNORB)*POTPAR(ISP)%KTOPHI(LNORB) &
       &        +NLPHIDOT(:,LNORB)*POTPAR(ISP)%KTOPHIDOT(LNORB)
        DK0ARR(:)=(AEPHI(:,LNORB)-NLPHI(:,LNORB))*POTPAR(ISP)%KTOPHI(LNORB) 
        LM2X=(L2X+1)**2
        DO L=0,L2X
          DO LN=1,LNX(ISP)
            IF(LOX(LN,ISP).EQ.L) THEN
              JBARARR(:,L+1)=NLPHIDOT(:,LN)*POTPAR(ISP)%JBARTOPHIDOT(LN)
              EXIT
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(AEPHI)
        DEALLOCATE(NLPHI)
        DEALLOCATE(NLPHIDOT)
        DO NN=1,NNS
          IF(SBAR(NN)%IAT1.NE.IATORB) CYCLE
          IAT2=SBAR(NN)%IAT2
          ISP2=ISPECIES(IAT2)
          IF(ISP2.NE.ISP) CYCLE
          R2(:)=R0(:,IAT2)+RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
       &                  +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
       &                  +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
          TONSITE=(IAT2.EQ.IATORB).AND.(MAXVAL(ABS(SBAR(NN)%IT(:))).EQ.0)
          IF(LM2X.NE.SBAR(NN)%N2) THEN
            CALL ERROR$MSG('LM2X INCONSISTENT WITH ARRAY DIMENSION OF SBAR')
            CALL ERROR$STOP('LMTO_NTBOOFR')
          END IF
!  
!         == LOOP OVER REAL SPACE GRID =========================================
          DO IP=1,NP
!           == DR IS THE DISTANCE FROM THE SECOND ATOM =========================
            DR(:)=P(:,IP)-R2(:)
            TSPHERE=(DOT_PRODUCT(DR,DR).LT.POTPAR(ISP2)%RAD**2)
            DIS=SQRT(SUM(DR**2))
            CALL SPHERICAL$YLM(LM2X,DR,YLM)
!           == INCLUDE DIFFERENCE BETWEEN AEPHI AND NLPHI ======================
            IF(TONSITE) THEN
              CALL RADIAL$VALUE(GID,NR,DK0ARR,DIS,VAL)
              ORB(IP)=ORB(IP)+VAL*YLM(LMORB)
            END IF
            IF(.NOT.TSPHERE) CYCLE
            IF(TONSITE) THEN
              CALL RADIAL$VALUE(GID,NR,K0ARR,DIS,VAL)
              ORB(IP)=ORB(IP)+VAL*YLM(LMORB)
            END IF
            LM=0
            DO L=0,L2X
              CALL RADIAL$VALUE(GID,NR,JBARARR(:,L+1),DIS,VALDOT)
              DO IM=1,2*L+1
                LM=LM+1
                ORB(IP)=ORB(IP)-VALDOT*YLM(LM)*SBAR(NN)%MAT(LMORB,LM) 
              ENDDO
            ENDDO
          ENDDO  ! END LOOP IP
        ENDDO  ! END LOOP NN
        DEALLOCATE(R)
        DEALLOCATE(K0ARR)
        DEALLOCATE(DK0ARR)
        DEALLOCATE(JBARARR)
      ENDDO   ! END LOOP NSP
      DEALLOCATE(YLM)
PRINT*,'MARKE 1G'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILED_NTBOOFR(TYPE,IATORB,LMNORB,NP,P,ORB)
!     **************************************************************************
!     **  MAPS THE SCREENED ORBITAL LMNORB AT ATOM IATORB                     **
!     **  ONTO THE GRID P                                                     **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2012 ************
      USE LMTO_MODULE, ONLY : ISPECIES,NSP,LNX,LOX,SBAR,K2,POTPAR,SBARLI1
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: TYPE    !
      INTEGER(4)  ,INTENT(IN) :: IATORB  ! INDEX OF CENTRAL ATOM
      INTEGER(4)  ,INTENT(IN) :: LMNORB  ! ORBITAL INDEX
      REAL(8)     ,INTENT(OUT):: ORB(NP) ! SCREENED ENVELOPE FUNCTIONS
      INTEGER(4)  ,INTENT(IN) :: NP      ! #(GRID POINTS)
      REAL(8)     ,INTENT(IN) :: P(3,NP) ! GRIDPOINT POSITIONS
      INTEGER(4)            :: NAT       ! #(ATOMS)
      REAL(8)   ,ALLOCATABLE:: R0(:,:)   ! ATOMIC POSITIONS
      INTEGER(4)            :: LMORB,LNORB,LORB,MORB
      INTEGER(4)            :: IAT2,L2X,LM2X
      INTEGER(4)            :: LNDOT
      INTEGER(4)            :: NNS 
      INTEGER(4)            :: NN,IP,ISP,ISP2,L,LM,IM,LN,LMN  !LOOP INDICES
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: LMXX
      REAL(8)               :: R2(3) ! NEIGHBOR ATOM POSITION
      REAL(8)               :: DR(3) ! DISTANCE VECTOR OF GRID POINT TO ATOM
      REAL(8)               :: DIS   ! DISTANCE OF GRID POINT TO ATOM
      REAL(8)               :: VAL,VALDOT
      LOGICAL(4)            :: TONSITE   
      LOGICAL(4)            :: TSPHERE   ! INSIDE AN AUGMENTATION SPHERE?
      REAL(8)  ,ALLOCATABLE :: YLM(:)    ! SPHERICAL HARMONICS
      REAL(8)  ,ALLOCATABLE :: R(:)      ! RADIAL GRID
      REAL(8)  ,ALLOCATABLE :: K0ARR(:)     ! BARE AUGMENTED HANKEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: JBARARR(:,:) ! SCREENED AUGMENTED BESSEL FUNCTION
!     **************************************************************************
PRINT*,'MARKE 1A'
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
PRINT*,'MARKE 1B'
!
!     ==========================================================================
!     == DETERMINE LM-, LN-, L- AND M-INDEX OF THE ORBITAL                    ==
!     ==========================================================================
      LMORB=0
      ISP=ISPECIES(IATORB)
      LMN=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        DO IM=1,2*L+1
          LMN=LMN+1
          IF(LMN.EQ.LMNORB) THEN
            LORB=LOX(LN,ISP)
            LNORB=LN
            LMORB=SBARLI1(LORB+1,ISP)
            EXIT
          END IF
        ENDDO
      ENDDO
      IF(LMORB.EQ.0) THEN
        CALL ERROR$MSG('INDEX ERROR')
        CALL ERROR$MSG('NTBOOFR')
      END IF
PRINT*,'MARKE 1C'
!
!     ==========================================================================
!     == 
!     ==========================================================================
      NNS=SIZE(SBAR)  !SBAR STRUCCONS. (GLOBAL)
PRINT*,'MARKE 1D'
!
!     ==========================================================================
!     ==  NOW DO HE AUGMENTATION                                              ==
!     ==========================================================================
      LMXX=MAXVAL(LOX(:,:)+1)**2
      ALLOCATE(YLM(LMXX))
      ISP=ISPECIES(IATORB)
      L2X=MAXVAL( LOX(:LNX(ISP),ISP) )
      GID=POTPAR(ISP)%TAILED%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      ALLOCATE(K0ARR(NR))
      ALLOCATE(JBARARR(NR,L2X+1))
      CALL RADIAL$R(GID,NR,R)
      K0ARR(:)=POTPAR(ISP)%TAILED%AEF(:,LNORB) 
      LM2X=(L2X+1)**2
      DO L=0,L2X
        DO LN=1,LNX(ISP)
          IF(LOX(LN,ISP).EQ.L) THEN
            LNDOT=POTPAR(ISP)%TAILED%LNDOT(LN)
            JBARARR(:,L+1)=POTPAR(ISP)%TAILED%NLF(:,LNDOT)
            EXIT
          END IF
        ENDDO
      ENDDO
!
      ORB(:)=0.D0
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IATORB) CYCLE
        IF(SBAR(NN)%IAT2.NE.IATORB) CYCLE
        IF(MAXVAL(ABS(SBAR(NN)%IT(:))).NE.0) CYCLE
        IAT2=SBAR(NN)%IAT2
        ISP2=ISPECIES(IAT2)
        R2(:)=R0(:,IAT2)
        IF(LM2X.NE.SBAR(NN)%N2) THEN
          CALL ERROR$MSG('LM2X INCONSISTENT WITH ARRAY DIMENSION OF SBAR')
          CALL ERROR$STOP('LMTO_NTBOOFR')
        END IF
!  
!       == LOOP OVER REAL SPACE GRID ===========================================
        DO IP=1,NP
!         == DR IS THE DISTANCE FROM THE SECOND ATOM ===========================
          DR(:)=P(:,IP)-R2(:)
          DIS=SQRT(SUM(DR**2))
          CALL SPHERICAL$YLM(LM2X,DR,YLM)
          CALL RADIAL$VALUE(GID,NR,K0ARR,DIS,VAL)
          ORB(IP)=ORB(IP)+VAL*YLM(LMORB)
          LM=0
          DO L=0,L2X
            CALL RADIAL$VALUE(GID,NR,JBARARR(:,L+1),DIS,VALDOT)
            DO IM=1,2*L+1
              LM=LM+1
              ORB(IP)=ORB(IP)-VALDOT*YLM(LM)*SBAR(NN)%MAT(LMORB,LM) 
            ENDDO
          ENDDO
        ENDDO  ! END LOOP IP
      ENDDO  ! END LOOP NN
      DEALLOCATE(R)
      DEALLOCATE(K0ARR)
      DEALLOCATE(JBARARR)
      DEALLOCATE(YLM)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GRIDAUGMENT(RBAS,NAT,R0,IAT1,LM1X,NP,P,ORB,ORBI)
!     **************************************************************************
!     ** MAPS THE AUGMENTATION FOR THE SCREENED ENVELOPE FUNCTION             **
!     **  AT ATOM IAT1 ONTO THE GRID P                                        **
!     **                                                                      **
!     ** ATTENTION: SCREENING IS DETERMINED BY ISCATT                         **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2010 ************
      USE LMTO_MODULE, ONLY : SBAR,K2,POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NP        ! #(GRID POINTS)
      REAL(8)   ,INTENT(IN) :: P(3,NP)   ! GRIDPOINT POSITIONS
      REAL(8)   ,INTENT(IN) :: RBAS(3,3) ! LATTICE VECTORS
      INTEGER(4),INTENT(IN) :: NAT       ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT) ! ATOMIC POSITIONS
      INTEGER(4),INTENT(IN) :: IAT1      ! INDEX OF CENTRAL ATOM
      INTEGER(4),INTENT(IN) :: LM1X      !#(ORBITALS ON CENTRAL ATOM)
      REAL(8)   ,INTENT(OUT):: ORB(NP,LM1X)  ! SCREENED ENVELOPE FUNCTIONS
      REAL(8)   ,INTENT(OUT):: ORBI(NP,LM1X) ! MULTICENTER EXPANSION OF ORB
      INTEGER(4)            :: NNS
      INTEGER(4)            :: LNX
      INTEGER(4)            :: LX
      INTEGER(4)            :: ISPECIES1(NAT)
      INTEGER(4)            :: GID,NR
      REAL(8)   ,ALLOCATABLE:: R(:)
      INTEGER(4),ALLOCATABLE:: LOX(:)
      INTEGER(4),ALLOCATABLE:: ISCATT(:)
      REAL(8)   ,ALLOCATABLE:: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE:: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE:: K0AUGARR(:,:)
      REAL(8)   ,ALLOCATABLE:: JBARAUGARR(:,:)
      REAL(8)   ,ALLOCATABLE:: K0(:)
      REAL(8)   ,ALLOCATABLE:: J0(:)
      REAL(8)   ,ALLOCATABLE:: JBAR(:)
      REAL(8)   ,ALLOCATABLE:: JBARAUG(:)
      REAL(8)   ,ALLOCATABLE:: K0AUG(:)
      REAL(8)   ,ALLOCATABLE:: QBARVEC(:)
      REAL(8)   ,ALLOCATABLE:: YLM(:)
      REAL(8)               :: K0VAL,K0DER,J0VAL,J0DER,JBARVAL
      REAL(8)               :: PHIVAL,PHIDER,PHIDOTVAL,PHIDOTDER
      REAL(8)               :: WK0PHI,WK0PHIDOT,WPHIPHIDOT
      REAL(8)               :: QBAR
      REAL(8)               :: SVAR
      REAL(8)               :: RAD
      REAL(8)               :: R2(3)
      REAL(8)               :: DR(3),DIS
      INTEGER(4)            :: NSP
      LOGICAL(4)            :: TONSITE
      INTEGER(4)            :: ISP,LN,L,NN,IAT2,LM1,LM2,LM2X,IP!LOOP INDICES ETC
      INTEGER(4)            :: IR,NFIL
REAL(8) :: X(10)
!     **************************************************************************
      NSP=SIZE(POTPAR)
      NNS=SIZE(SBAR)  !SBAR STRUCCONS. (GLOBAL)
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES1)
      ORB(:,:)=0.D0
      ORBI(:,:)=0.D0
      DO ISP=1,NSP
        RAD=POTPAR(ISP)%RAD
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(ISCATT(LNX))
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        ALLOCATE(AEPHI(NR,LNX))
        ALLOCATE(NLPHIDOT(NR,LNX))
        CALL SETUP$GETR8A('NLPHIDOT',NR*LNX,NLPHIDOT)
        CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
        CALL SETUP$ISELECT(0)
!
!       ========================================================================
!       ==  MATCH PARTIAL WAVES TO K AND JBAR                                 ==
!       ========================================================================
        LX=MAXVAL(LOX)
        ALLOCATE(K0AUGARR(NR,LX+1))
        ALLOCATE(JBARAUGARR(NR,LX+1))
        K0AUGARR(:,:)=0.D0
        JBARAUGARR(:,:)=0.D0
        ALLOCATE(QBARVEC((LX+1)**2))
        QBARVEC=0.D0
        DO LN=1,LNX
          IF(ISCATT(LN).GT.0) CYCLE 
          L=LOX(LN)
          QBAR=POTPAR(ISP)%QBAR(LN)
          CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,J0VAL,J0DER)
          CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,K0VAL,K0DER)
          CALL RADIAL$VALUE(GID,NR,NLPHIDOT(:,LN),RAD,PHIDOTVAL)
          CALL RADIAL$DERIVATIVE(GID,NR,NLPHIDOT(:,LN),RAD,PHIDOTDER)
          CALL RADIAL$VALUE(GID,NR,AEPHI(:,LN),RAD,PHIVAL)
          CALL RADIAL$DERIVATIVE(GID,NR,AEPHI(:,LN),RAD,PHIDER)
          JBARVAL=J0VAL-K0VAL*QBAR
!         == NLPHIDOT MATCHES TO JBAR ==========================================
          JBARAUGARR(:,L+1)=NLPHIDOT(:,LN)/PHIDOTVAL*JBARVAL
          WK0PHIDOT=K0VAL*PHIDOTDER-K0DER*PHIDOTVAL
          WK0PHI=K0VAL*PHIDER-K0DER*PHIVAL
          WPHIPHIDOT=PHIVAL*PHIDOTDER-PHIDER*PHIDOTVAL
!         == AEPHI MATCHES TO K0 ===============================================
          K0AUGARR(:,L+1)=(AEPHI(:,LN)*WK0PHIDOT-NLPHIDOT(:,LN)*WK0PHI) &
      &                  /WPHIPHIDOT
          QBARVEC(L**2+1:(L+1)**2)=QBAR
!PRINT*,'WKJ',WK0PHIDOT,WK0PHI,WPHIPHIDOT
        ENDDO
!PRINT*,'ISCATT ',ISCATT(1:LNX)
!PRINT*,'LOX    ',LOX(1:LNX)
!!$#PRINT*,'QBARVEC ',QBARVEC
!!$NFIL=12
!!$OPEN(UNIT=NFIL,FILE='XX.DAT')
!!$DO IR=10,NR
!!$  CALL LMTO$SOLIDBESSELRAD(0,R(IR),K2,X(1),J0DER)
!!$  CALL LMTO$SOLIDHANKELRAD(0,R(IR),K2,X(2),K0DER)
!!$  CALL LMTO$SOLIDBESSELRAD(1,R(IR),K2,X(3),J0DER)
!!$  CALL LMTO$SOLIDHANKELRAD(1,R(IR),K2,X(4),K0DER)
!!$  X(1)=X(1)-X(2)*QBARVEC(1)
!!$  X(3)=X(3)-X(4)*QBARVEC(2)
!!$  WRITE(NFIL,*)R(IR),K0AUGARR(IR,:),JBARAUGARR(IR,:),X(1:4)
!!$ENDDO
!!$CLOSE(NFIL)
!!$STOP
!
!
!       ========================================================================
!       ==                                                                    ==
!       ========================================================================
        ALLOCATE(K0AUG((LX+1)**2))
        ALLOCATE(JBARAUG((LX+1)**2))
        ALLOCATE(K0((LX+1)**2))
        ALLOCATE(J0((LX+1)**2))
        ALLOCATE(JBAR((LX+1)**2))
        ALLOCATE(YLM((LX+1)**2))
!
        DO NN=1,NNS
          IF(SBAR(NN)%IAT1.NE.IAT1) CYCLE
!IF(SBAR(NN)%IAT2.EQ.IAT1) CYCLE
          IAT2=SBAR(NN)%IAT2
          IF(ISPECIES1(IAT2).NE.ISP) CYCLE
          LM2X=SBAR(NN)%N2
          R2(:)=R0(:,IAT2)+RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
     &                    +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
     &                    +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
          TONSITE=(IAT2.EQ.IAT1).AND.(MAXVAL(ABS(SBAR(NN)%IT(:))).EQ.0)
          DO IP=1,NP
            DR(:)=P(:,IP)-R2(:)
            DIS=SQRT(SUM(DR**2))
            IF(DIS.GT.RAD) CYCLE
!
!           == DETERMINE HANKEL AND SCREENED BESSEL FUNCTION AT IAT2 ===========
            CALL  LMTO$SOLIDHANKEL(DR,RAD,K2,LM2X,K0(1:LM2X))
            CALL  LMTO$SOLIDBESSEL(DR,K2,LM2X,J0(1:LM2X))
            JBAR(:LM2X)=J0(:LM2X)-K0(:LM2X)*QBARVEC(:LM2X)
            CALL SPHERICAL$YLM(LM2X,DR,YLM)
            DO L=0,LX
              LM1=L**2+1
              LM2=(L+1)**2
              CALL RADIAL$VALUE(GID,NR,K0AUGARR(:,L+1),DIS,SVAR)
              K0AUG(LM1:LM2)=SVAR*YLM(LM1:LM2)
              CALL RADIAL$VALUE(GID,NR,JBARAUGARR(:,L+1),DIS,SVAR)
              JBARAUG(LM1:LM2)=SVAR*YLM(LM1:LM2)
            ENDDO
!
            DO LM1=1,LM1X
              ORBI(IP,LM1)=ORBI(IP,LM1) &
      &                           -DOT_PRODUCT(JBAR(1:LM2X),SBAR(NN)%MAT(LM1,:)) !C
              ORB(IP,LM1)=ORB(IP,LM1) &
      &                        -DOT_PRODUCT(JBARAUG(1:LM2X),SBAR(NN)%MAT(LM1,:))  !C
              IF(TONSITE) THEN
                ORBI(IP,LM1)=ORBI(IP,LM1)+K0(LM1)
                ORB(IP,LM1) =ORB(IP,LM1) +K0AUG(LM1)
              END IF
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(K0AUG)
        DEALLOCATE(JBARAUG)
        DEALLOCATE(K0)
        DEALLOCATE(J0)
        DEALLOCATE(JBAR)
        DEALLOCATE(LOX)
        DEALLOCATE(ISCATT)
        DEALLOCATE(NLPHIDOT)
        DEALLOCATE(AEPHI)
        DEALLOCATE(R)
        DEALLOCATE(K0AUGARR)
        DEALLOCATE(JBARAUGARR)
        DEALLOCATE(QBARVEC)
        DEALLOCATE(YLM)
      ENDDO    
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GRIDENVELOPE(RBAS,NAT,R0,IAT1,LM1X,NP,P,ORB,ORBI)
!     **************************************************************************
!     ** MAPS THE SCREENED ENVELOPE FUNCTION AT ATOM IAT1 ONTO THE GRID P     **
!     **                                                                      **
!     ** ATTENTION: SCREENING IS DETERMINED BY ISCATT                         **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2010 ************
      USE LMTO_MODULE, ONLY : SBAR,K2,POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NP        ! #(GRID POINTS)
      REAL(8)   ,INTENT(IN) :: P(3,NP)   ! GRIDPOINT POSITIONS
      REAL(8)   ,INTENT(IN) :: RBAS(3,3) ! LATTICE VECTORS
      INTEGER(4),INTENT(IN) :: NAT       ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT) ! ATOMIC POSITIONS
      INTEGER(4),INTENT(IN) :: IAT1      ! INDEX OF CENTRAL ATOM
      INTEGER(4),INTENT(IN) :: LM1X      !#(ORBITALS ON CENTRAL ATOM)
      REAL(8)   ,INTENT(OUT):: ORB(NP,LM1X)  ! SCREENED ENVELOPE FUNCTIONS
      REAL(8)   ,INTENT(OUT):: ORBI(NP,LM1X) ! MULTICENTER EXPANSION OF ORB
      INTEGER(4)            :: NNS
      INTEGER(4)            :: NN,LM1,IP,IAT,ISP,LN,L  !LOOP INDICES
      INTEGER(4)            :: LNX
      INTEGER(4)            :: LMXX
      INTEGER(4)            :: IAT2,LM2X
      REAL(8)               :: RAD(NAT)  ! ATOMIC RADII( RCV)
      REAL(8)               :: R2(3) ! NEIGHBOR ATOM POSITION
      REAL(8)               :: DR(3) ! DISTANCE OF GRID POINT TO ATOM
      LOGICAL(4)            :: TONSITE   
      LOGICAL(4)            :: TSPHERE   ! INSIDE AN AUGMENTATION SPHERE?
      REAL(8)  ,ALLOCATABLE :: CVEC(:,:) ! COEFFICIENTS OF BARE HANKEL FUNCTNS.
      REAL(8)  ,ALLOCATABLE :: K0(:)     ! BARE SOLID HANKEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: JBAR(:)   ! SCREENED SOLID BESSEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: J0(:)     ! BARE SOLID BESSEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: QBARVEC(:,:)
      INTEGER(4),ALLOCATABLE::LOX(:),ISCATT(:)
      INTEGER(4)            :: ISPECIES1(NAT)
!     **************************************************************************
      NNS=SIZE(SBAR)  !SBAR STRUCCONS. (GLOBAL)
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES1)
!     == DETERMINE LMXX #(ANGULAR MOMENTA ON NEIGHBORING SITE) =================
      LMXX=0
      DO NN=1,NNS
         IF(SBAR(NN)%IAT1.NE.IAT1) CYCLE
        LMXX=MAX(LMXX,SBAR(NN)%N2)
      ENDDO
      ALLOCATE(CVEC(LMXX,LM1X))
      ALLOCATE(K0(LMXX))
      ALLOCATE(JBAR(LMXX))
      ALLOCATE(J0(LMXX))
!
!     == DETERMINE QBAR ========================================================
      ALLOCATE(QBARVEC(LMXX,NAT))
      QBARVEC(:,:)=0.D0
      DO IAT=1,NAT
        ISP=ISPECIES1(IAT)
        CALL SETUP$ISELECT(ISP)
        RAD(IAT)=POTPAR(ISP)%RAD
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(ISCATT(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
        DO LN=1,LNX
          IF(ISCATT(LN).NE.0) CYCLE
          L=LOX(LN)
          QBARVEC(L**2+1:(L+1)**2,IAT)=POTPAR(ISP)%QBAR(LN)
        ENDDO
        DEALLOCATE(LOX)
        DEALLOCATE(ISCATT)
        CALL SETUP$ISELECT(0)
      ENDDO
!
!     == LOOP OVER ALL NEIGHBORS AND GRID POINTS ===============================
      ORB(:,:)=0.D0
      ORBI(:,:)=0.D0 !SPHERE CONTRIBUTION FROM MULTICENTER EXPANSION
      DO NN=1,NNS
         IF(SBAR(NN)%IAT1.NE.IAT1) CYCLE
         IF(SBAR(NN)%N1.NE.LM1X) THEN
           CALL ERROR$MSG('INTERNAL ERROR')
           CALL ERROR$STOP('LMTO$PLOTLOCORB')
         END IF
         IAT2=SBAR(NN)%IAT2
         LM2X=SBAR(NN)%N2
         R2(:)=R0(:,IAT2)+RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
     &                   +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
     &                   +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
!        == CVEC=1+QBAR*SBAR ===================================================
         TONSITE=(IAT2.EQ.IAT1).AND.(MAXVAL(ABS(SBAR(NN)%IT(:))).EQ.0)
         DO LM1=1,LM1X
           CVEC(:LM2X,LM1)=QBARVEC(:LM2X,IAT2)*SBAR(NN)%MAT(LM1,:)    !C
           IF(TONSITE)CVEC(LM1,LM1)=CVEC(LM1,LM1)+1.D0
         ENDDO
!        == LOOP OVER REAL SPACE GRID ==========================================
         DO IP=1,NP
!          == DR IS THE DISTANCE FROM THE SECOND ATOM ==========================
           DR(:)=P(:,IP)-R2(:)
           TSPHERE=(DOT_PRODUCT(DR,DR).LT.RAD(IAT2)**2)
!
!          == DETERMINE BARE HANKEL AND SCREENED BESSEL FUNCTION AT IAT2 =======
           CALL  LMTO$SOLIDHANKEL(DR,RAD(IAT2),K2,LM2X,K0(1:LM2X))
           IF(TSPHERE) THEN
             CALL  LMTO$SOLIDBESSEL(DR,K2,LM2X,J0(1:LM2X))
             JBAR(:LM2X)=J0(:LM2X)-K0(:LM2X)*QBARVEC(:LM2X,IAT2)
           END IF
!          ==  
           DO LM1=1,LM1X
!            == |KBAR>=|K0>*(1+QBAR*SBART) =====================================
             ORB(IP,LM1)=ORB(IP,LM1)+DOT_PRODUCT(K0(1:LM2X),CVEC(1:LM2X,LM1))
!            == -|JBAR>SBAR=-(|J0>-|K0>QBAR)*SBAR ==============================
             IF(TSPHERE) THEN
               ORBI(IP,LM1)=ORBI(IP,LM1) &
      &                        -DOT_PRODUCT(JBAR(1:LM2X),SBAR(NN)%MAT(LM1,:)) !C
               IF(TONSITE) ORBI(IP,LM1)=ORBI(IP,LM1)+K0(LM1)
             END IF
           ENDDO
         ENDDO
      ENDDO
      DEALLOCATE(CVEC)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GRIDGAUSS(RBAS,NAT,R0,IAT0,NORB,NP,P,ORB)
!     **************************************************************************
!     ** MAPS THE SCREENED ENVELOPE FUNCTION AT ATOM IAT1 ONTO THE GRID P     **
!     ** THE FUNCTION IS EXPRESSED BY GAUSS FUNCTIONS                         **
!     **                                                                      **
!     ** ATTENTION: SCREENING IS DETERMINED BY ISCATT                         **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2010 ************
      USE LMTO_MODULE, ONLY : GAUSSORB
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NP        ! #(GRID POINTS)
      REAL(8)   ,INTENT(IN) :: P(3,NP)   ! GRIDPOINT POSITIONS
      REAL(8)   ,INTENT(IN) :: RBAS(3,3) ! LATTICE VECTORS
      INTEGER(4),INTENT(IN) :: NAT       ! #(ATOMS)
      REAL(8)   ,INTENT(IN) :: R0(3,NAT) ! ATOMIC POSITIONS
      INTEGER(4),INTENT(IN) :: IAT0      ! INDEX OF CENTRAL ATOM
      INTEGER(4),INTENT(IN) :: NORB      !#(ORBITALS ON CENTRAL ATOM)
      REAL(8)   ,INTENT(OUT):: ORB(NP,NORB)  ! SCREENED ENVELOPE FUNCTIONS
      REAL(8)               :: DR(3) ! DISTANCE OF GRID POINT TO ATOM
      INTEGER(4)            :: NIJK
      INTEGER(4)            :: NE
      REAL(8)               :: SVAR
      INTEGER(4)            :: IP,IORB,IE
!     **************************************************************************
                                         CALL TRACE$PUSH('LMTO_GRIDGAUSS')
!
!     == LOOP OVER REAL SPACE GRID ==========================================
      NE=GAUSSORB(IAT0)%NE
      NIJK=GAUSSORB(IAT0)%NIJK
      ORB(:,:)=0.D0
      DO IP=1,NP
        
!       == DR IS THE DISTANCE FROM THE CENTRAL ATOM ============================
        DR(:)=P(:,IP)-R0(:,IAT0)
        DO IORB=1,NORB
          ORB(IP,IORB)=0.D0
          DO IE=1,NE
            CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,GAUSSORB(IAT0)%E(IE) &
     &                   ,GAUSSORB(IAT0)%C(:,IE,IORB),DR,SVAR)
            ORB(IP,IORB)=ORB(IP,IORB)+SVAR
          ENDDO
        ENDDO
      ENDDO
                                        CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_WRITECUBEFILE(NFIL,NAT,Z,R,ORIGIN,BOX,N1,N2,N3,DATA)
!     **************************************************************************
!     **  WRITE A GAUSSIAN CUBE FILE (EXTENSION .CUB) WITH VOLUMINETRIC DATA  **
!     **                                                                      **
!     ** REMARK:                                                              **
!     ** UNITS WRITTEN ARE ABOHR, CONSISTENT WITH AVOGADROS IMPLEMENTATION.   **
!     ** THE SPECS REQUIRE N1,N2,N3 TO BE MULTIPLIED BY -1 IF ABOHR ARE USED  **
!     ** AND ANGSTROM IS THE UNIT IF THEY ARE POSITIVE.                       **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2010 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: NAT       ! NUMBER OF ATOMS
      REAL(8)   ,INTENT(IN) :: Z(NAT)    !ATOMIC NUMBER
      REAL(8)   ,INTENT(IN) :: R(3,NAT)
      REAL(8)   ,INTENT(IN) :: ORIGIN(3)
      REAL(8)   ,INTENT(IN) :: BOX(3,3)
      INTEGER(4),INTENT(IN) :: N1,N2,N3
      REAL(8)   ,INTENT(IN) :: DATA(N1,N2,N3)
      REAL(8)               :: ANGSTROM
      REAL(8)               :: SCALE
      INTEGER(4)            :: IAT,I,J,K
!     **************************************************************************
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      SCALE=1.D0
!      SCALE=1.D0/ANGSTROM
      WRITE(NFIL,FMT='("CP-PAW CUBE FILE")')
      WRITE(NFIL,FMT='("NOCHN KOMMENTAR")')
      WRITE(NFIL,FMT='(I5,3F12.6)')NAT,ORIGIN*SCALE
      WRITE(NFIL,FMT='(I5,3F12.6)')-N1,BOX(:,1)/REAL(N1,KIND=8)*SCALE
      WRITE(NFIL,FMT='(I5,3F12.6)')-N2,BOX(:,2)/REAL(N2,KIND=8)*SCALE
      WRITE(NFIL,FMT='(I5,3F12.6)')-N3,BOX(:,3)/REAL(N3,KIND=8)*SCALE
      DO IAT=1,NAT
        WRITE(NFIL,FMT='(I5,4F12.6)')NINT(Z(IAT)),0.D0,R(:,IAT)*SCALE
      ENDDO  
      WRITE(NFIL,FMT='(6(E12.5," "))')(((DATA(I,J,K),K=1,N3),J=1,N2),I=1,N1)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PLOTRADIAL()
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY: GAUSSORBAUG
      IMPLICIT NONE
      REAL(8)    ,PARAMETER :: RAD=5.D0
      INTEGER(4)            :: NIJK
      INTEGER(4)            :: NE
      INTEGER(4)            :: NORB
      INTEGER(4)            :: NAT
      REAL(8)   ,ALLOCATABLE:: E(:)
      REAL(8)   ,ALLOCATABLE:: C(:,:)
      INTEGER(4)            :: IAT,IORB
      CHARACTER(128)        :: FILE
      CHARACTER(128)        :: STRING
!     **************************************************************************
      CALL ATOMLIST$NATOM(NAT)
      DO IAT=1,NAT
        NORB=GAUSSORBAUG(IAT)%NORB
        NIJK=GAUSSORBAUG(IAT)%NIJK
        NE=GAUSSORBAUG(IAT)%NE
        ALLOCATE(E(NE))
        ALLOCATE(C(NIJK,NE))
        E=GAUSSORBAUG(IAT)%E
        DO IORB=1,NORB
          WRITE(STRING,*)IAT
          FILE='PLOTRADIAL_'//TRIM(ADJUSTL(STRING))
          WRITE(STRING,*)IORB
          FILE=TRIM(ADJUSTL(FILE))//'_'//TRIM(ADJUSTL(STRING))//'.DAT'
          C(:,:)=GAUSSORBAUG(IAT)%C(:,:,IORB)
!4*PI*R^2*F^2 IS THE NORM OF THE ORBITAL AS IN THE OVERLAP MATRIX.
          CALL GAUSSIAN$PLOTRADIAL(FILE,RAD,NIJK,NE,E,C)
        ENDDO
        DEALLOCATE(E)
        DEALLOCATE(C)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$PLOTWAVE(NFIL,IDIM0,IB0,IKPT0,ISPIN0,NR1,NR2,NR3)
!     **************************************************************************
!     ** WRAPPER CALLED FROM PAW_GRAPHICS TO PLOT A WAVE FUNCTION IN NTBO     **
!     ** REPRESENTATION                                                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)  :: NFIL
      INTEGER(4)   ,INTENT(IN)  :: IDIM0
      INTEGER(4)   ,INTENT(IN)  :: IB0
      INTEGER(4)   ,INTENT(IN)  :: IKPT0
      INTEGER(4)   ,INTENT(IN)  :: ISPIN0
      INTEGER(4)   ,INTENT(IN)  :: NR1,NR2,NR3
!     **************************************************************************
      CALL LMTO$PLOTWAVE_TAILED(NFIL,IDIM0,IB0,IKPT0,ISPIN0,NR1,NR2,NR3)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$PLOTWAVE_TAILED(NFIL,IDIM0,IB0,IKPT0,ISPIN0,NR1,NR2,NR3)
!     **************************************************************************
!     **************************************************************************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : TON,GAUSSORB,ISPECIES
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)  :: NFIL
      INTEGER(4)   ,INTENT(IN)  :: IDIM0
      INTEGER(4)   ,INTENT(IN)  :: IB0
      INTEGER(4)   ,INTENT(IN)  :: IKPT0
      INTEGER(4)   ,INTENT(IN)  :: ISPIN0
      INTEGER(4)   ,INTENT(IN)  :: NR1,NR2,NR3
      REAL(8)                   :: RBAS(3,3)
      INTEGER(4)                :: NAT
      REAL(8)       ,ALLOCATABLE:: R0(:,:)
      CHARACTER(32),ALLOCATABLE :: NAME(:)
      REAL(8)      ,ALLOCATABLE :: ZAT(:)
      REAL(8)      ,ALLOCATABLE :: Q(:)
      INTEGER(4)                :: NP
      REAL(8)      ,ALLOCATABLE :: P(:,:)
      COMPLEX(8)   ,ALLOCATABLE :: COEFF(:,:)
      COMPLEX(8)   ,ALLOCATABLE :: COEFF1(:)
      COMPLEX(8)   ,ALLOCATABLE :: CVEC(:)
      CHARACTER(64)             :: TITLE='LMTO_WAVEFILE'
      COMPLEX(8)  ,ALLOCATABLE  :: WAVE(:)
      COMPLEX(8)  ,ALLOCATABLE  :: WAVE1(:,:,:)
      REAL(8)                   :: PI
      COMPLEX(8)                :: CI2PI
      COMPLEX(8)                :: EIKT
      COMPLEX(8)                :: CSVAR
      REAL(8)                   :: T(3)
      REAL(8)                   :: XK(3)
      REAL(8)     ,ALLOCATABLE  :: XKARR(:,:)
      INTEGER(4)                :: LMNX,LM1X
      INTEGER(4)                :: NKPT
      INTEGER(4)                :: IP,IR1,IR2,IR3,IAT,IT1,IT2,IT3,IPRO,IORB,IE
      INTEGER(4)                :: IND,I,J,K,ISP
      INTEGER(4)                :: IB1,IB2,IBH,IB
      REAL(8)                   :: X,Y,Z,R2
      INTEGER(4)                :: NIJK,NE,NORB
      INTEGER(4)                :: NPRO,NB,NBH
      COMPLEX(8)                :: CSVAR1,CSVAR2
      COMPLEX(8)  ,PARAMETER    :: CI=(0.D0,1.D0)
      LOGICAL(4)                :: TINV
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO$PLOTWAVE')
      IF(.NOT.TON) RETURN
      PI=4.D0*ATAN(1.D0)
      CI2PI=(0.D0,1.D0)*2.D0*PI
!
!     ==========================================================================
!     ==  COLLECT DATA                                                        ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      ALLOCATE(Q(NAT))
      ALLOCATE(NAME(NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      CALL ATOMLIST$GETCHA('NAME',0,NAT,NAME)
      CALL ATOMLIST$GETR8A('Q',0,NAT,Q)
      CALL DYNOCC$GETI4('NKPT',NKPT)
      ALLOCATE(XKARR(3,NKPT))
      CALL DYNOCC$GETR8A('XK',3*NKPT,XKARR)
      XK(:)=XKARR(:,IKPT0)
      DEALLOCATE(XKARR)
!
!     ==========================================================================
!     ==  PREPARE GRID                                                        ==
!     ==========================================================================
      NP=NR1*NR2*NR3
      ALLOCATE(P(3,NP))
      IP=0
      DO IR3=1,NR3
        DO IR2=1,NR2
          DO IR1=1,NR1
            IP=IP+1
            P(:,IP)=RBAS(:,1)*REAL(IR1-1,KIND=8)/REAL(NR1,KIND=8) &
    &              +RBAS(:,2)*REAL(IR2-1,KIND=8)/REAL(NR2,KIND=8) &
    &              +RBAS(:,3)*REAL(IR3-1,KIND=8)/REAL(NR3,KIND=8) 
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
      CALL WAVES_SELECTWV(IKPT0,ISPIN0)
      CALL PLANEWAVE$SELECT(GSET%ID)
      CALL PLANEWAVE$GETL4('TINV',TINV)
      NB=THIS%NB
      NBH=THIS%NBH
      NPRO=MAP%NPRO
!
!     ==========================================================================
!     == COEFFICIENTS FOR THE NTBO BASIS                                      ==
!     ==========================================================================
!!$WRITE(*,FMT='(100("="),T10," EIGENVECTORS ")')
!!$DO I=1,THIS%NB
!!$  WRITE(*,FMT='(100("(",2F10.5,") "))')THIS%EIGVEC(:,I)
!!$ENDDO
!!$WRITE(*,FMT='(100("="),T10," TBC ")')
!!$DO I=1,THIS%NBH
!!$  WRITE(*,FMT='(100("(",2F10.5,") "))')THIS%TBC(1,I,:)
!!$ENDDO
      ALLOCATE(CVEC(NPRO))
      CVEC(:)=(0.D0,0.D0)
      IF(TINV) THEN
        DO IBH=1,NBH
          IB1=2*IBH-1
          IB2=2*IBH
          CSVAR1=0.5D0*(THIS%EIGVEC(IB1,IB0)-CI*THIS%EIGVEC(IB2,IB0))
          CSVAR2=0.5D0*(THIS%EIGVEC(IB1,IB0)+CI*THIS%EIGVEC(IB2,IB0))
          CVEC(:)=CVEC(:)+THIS%TBC(IDIM0,IBH,:)*CSVAR1 &
      &            +CONJG(THIS%TBC(IDIM0,IBH,:))*CSVAR2
        ENDDO
      ELSE
        DO IB1=1,NB
          CSVAR1=THIS%EIGVEC(IB1,IB0)
          CVEC(:)=CVEC(:)+THIS%TBC(IDIM0,IB1,:)*CSVAR1
        ENDDO             
      END IF
!!$WRITE(*,FMT='(100("="),T10," WAVE FUNCTION COEFFICIENTS ")')
!!$WRITE(*,FMT='(100("(",2F10.5,") "))')CVEC
!
!     ==========================================================================
!     == MAP WACE FUNCTION ONTO GRID                                          ==
!     ==========================================================================
      ALLOCATE(WAVE(NP))
      WAVE=(0.D0,0.D0)
      DO IT1=-1,1
        DO IT2=-1,1
          DO IT3=-1,1
            T(:)=RBAS(:,1)*REAL(IT1,KIND=8) &
     &          +RBAS(:,2)*REAL(IT2,KIND=8) &
     &          +RBAS(:,3)*REAL(IT3,KIND=8)
            EIKT=EXP(CI2PI*(XK(1)*REAL(IT1)+XK(2)*REAL(IT2)+XK(3)*REAL(IT3)))
            IPRO=0
            DO IAT=1,NAT
CALL ERROR$MSG('CODING INCOMPLETE')
CALL ERROR$MSG('REWRITE FROM GAUSSIAN TO TAILED REPRESENTATION')
CALL ERROR$STOP('LMTO$PLOTWAVE_TAILED')
!!$              NIJK=GAUSSORB(IAT)%NIJK
!!$              NE=GAUSSORB(IAT)%NE
!!$              NORB=GAUSSORB(IAT)%NORB
!!$              ALLOCATE(COEFF(NIJK,NE))
!!$              ALLOCATE(COEFF1(NIJK))
!!$              COEFF(:,:)=(0.D0,0.D0)
!!$              DO IORB=1,NORB
!!$                IPRO=IPRO+1
!!$                COEFF(:,:)=COEFF(:,:)+GAUSSORB(IAT)%C(:,:,IORB)*CVEC(IPRO)
!!$              ENDDO
!!$              COEFF(:,:)=COEFF(:,:)*EIKT
!!$              DO IP=1,NP
!!$                X=P(1,IP)-R0(1,IAT)-T(1)
!!$                Y=P(2,IP)-R0(2,IAT)-T(2)
!!$                Z=P(3,IP)-R0(3,IAT)-T(3)
!!$                R2=X**2+Y**2+Z**2
!!$                IF(R2.GT.10.D0**2) CYCLE
!!$
!!$                COEFF1(:)=(0.D0,0.D0)
!!$                DO IE=1,NE
!!$                  COEFF1(:)=COEFF1(:)+COEFF(:,IE)*EXP(-GAUSSORB(IAT)%E(IE)*R2)
!!$                ENDDO
!!$                CSVAR=(0.D0,0.D0)
!!$                DO IND=1,NIJK
!!$                  CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I,J,K)
!!$                  CSVAR=CSVAR+(X**I)*(Y**J)*(Z**K)*COEFF1(IND)
!!$                ENDDO
!!$                WAVE(IP)=WAVE(IP)+CSVAR
!!$              ENDDO
!!$              DEALLOCATE(COEFF)
!!$              DEALLOCATE(COEFF1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(CVEC)
!
!     ==========================================================================
!     == WRITE DATA TO FILE                                                   ==
!     ==========================================================================
      ALLOCATE(ZAT(NAT))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',ZAT(IAT))
        CALL SETUP$ISELECT(0)
      ENDDO
      CALL WRITEWAVEPLOTC(NFIL,TITLE,RBAS,NAT,R0,ZAT,Q,NAME,XK,NR1,NR2,NR3,WAVE)
PRINT*,'IKPT ',IKPT0,' ISPIN=',ISPIN0,' IB=',IB,' IB0 ',IB0,' XK ',XK
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TESTPLOTRADIALGAUSSALONE(ID,L,NPOW,NE,E,C)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE        
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: L
      INTEGER(4)  ,INTENT(IN) :: NPOW
      INTEGER(4)  ,INTENT(IN) :: NE
      REAL(8)     ,INTENT(IN) :: E(NE)
      REAL(8)     ,INTENT(IN) :: C(NPOW,NE)
      INTEGER(4)              :: NFIL
      INTEGER(4)  ,SAVE       :: GID=0
      REAL(8)     ,ALLOCATABLE:: R(:)
      REAL(8)     ,ALLOCATABLE:: G(:)
      REAL(8)                 :: R1,DEX
      INTEGER(4)              :: NR
      INTEGER(4)              :: IR,IE,I
!     **************************************************************************
      IF(GID.EQ.0) THEN
        CALL RADIAL$NEW('SHLOG',GID)
        CALL RADIAL$GRIDPARAMETERS(0.1D0,0.2D0,50.D0,R1,DEX,NR)
        CALL RADIAL$SETI4(GID,'NR',NR)
        CALL RADIAL$SETR8(GID,'DEX',DEX)
        CALL RADIAL$SETR8(GID,'R1',R1)
      END IF
      CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,TRIM(ID)//'.DAT')
      CALL FILEHANDLER$UNIT('HOOK',NFIL)
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      ALLOCATE(G(NR))
      CALL RADIAL$R(GID,NR,R)
      G(:)=0.D0
      DO IE=1,NE
        DO I=0,NPOW-1
          G(:)=G(:)+R(:)**(L+2*I)*EXP(-E(IE)*R(:)**2)*C(I+1,IE)
        ENDDO  
      ENDDO
      DO IR=1,NR
        WRITE(NFIL,FMT='(10F20.5)')R(IR),G(IR)
      ENDDO
      CALL FILEHANDLER$CLOSE('HOOK')
      CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOK')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TESTPLOTRADIALGAUSS(ID,L,GID,NR,F,NPOW,NE,E,C)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE        
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: L
      REAL(8)     ,INTENT(IN) :: F(NR)
      INTEGER(4)  ,INTENT(IN) :: NPOW
      INTEGER(4)  ,INTENT(IN) :: NE
      REAL(8)     ,INTENT(IN) :: E(NE)
      REAL(8)     ,INTENT(IN) :: C(NPOW,NE)
      INTEGER(4)              :: NFIL
      REAL(8)                 :: R(NR)
      REAL(8)                 :: G(NR)
      INTEGER(4)              :: IR,IE,I
!     **************************************************************************
      CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,TRIM(ID)//'.DAT')
      CALL FILEHANDLER$UNIT('HOOK',NFIL)
      CALL RADIAL$R(GID,NR,R)
      G(:)=0.D0
      DO IE=1,NE
        DO I=0,NPOW-1
          G(:)=G(:)+R(:)**(L+2*I)*EXP(-E(IE)*R(:)**2)*C(I+1,IE)
        ENDDO  
      ENDDO
      DO IR=1,NR
        WRITE(NFIL,FMT='(10F20.5)')R(IR),F(IR),G(IR)
      ENDDO
      CALL FILEHANDLER$CLOSE('HOOK')
      CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOK')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TESTTAILEDP(NX)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY: POTPAR,NSP
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      INTEGER(4) :: ISP,LN1,LN2,LN3,LN4
      INTEGER(4) :: L1,L2,L3,L4
      INTEGER(4) :: LR1,LR2
      REAL(8)    :: VAL1,VAL2,SVAR0,SVAR1,SVAR2
      INTEGER(4) :: NE,NIJK
      INTEGER(4) :: IE1,IE2,I,J
      INTEGER(4) :: IP1,IP2,IS,IT
      INTEGER(4) :: LNX
      INTEGER(4) :: COUNT
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: E(:)
!     **************************************************************************
      DO ISP=1,NSP
        LNX=POTPAR(ISP)%TAILED%LNX
        ALLOCATE(LOX(LNX))
        LOX=POTPAR(ISP)%TAILED%LOX
        NE=POTPAR(ISP)%TAILED%SINGLE%NE
        NIJK=POTPAR(ISP)%TAILED%SINGLE%NIJK
        ALLOCATE(E(NE))
        E=POTPAR(ISP)%TAILED%SINGLE%E
!
      COUNT=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO LN2=LN1,LNX
          L2=LOX(LN2)
          DO LN3=1,LNX
            L3=LOX(LN3)
            DO LN4=LN3,LNX
              L4=LOX(LN4)  
              DO LR1=ABS(L1-L2),L1+L2,2
!               ================================================================
                COUNT=COUNT+1
                DO LR2=ABS(L3-L4),L3+L4,2
                  IF(LR2.EQ.LR1) THEN
                    CALL LMTO_TAILEDINDEX_P(NX,LNX,LOX,LN1,LN2,LR1,IP1)
                    CALL LMTO_TAILEDINDEX_P(NX,LNX,LOX,LN3,LN4,LR2,IP2)
                    CALL LMTO_TESTPLOTRADIALGAUSSALONE('RHO',LR1,NIJK,NE,E &
                                        ,POTPAR(ISP)%TAILED%PRODRHO%C(:,:,IP1))
                    CALL LMTO_TESTPLOTRADIALGAUSSALONE('POT',LR2,NIJK,NE,E &
                                        ,POTPAR(ISP)%TAILED%PRODPOT%C(:,:,IP2))
                    VAL1=0.D0
                    DO IE1=1,NE
                      DO IE2=1,NE
                        DO I=1,NIJK
                          DO J=1,NIJK
!                           == R^2 * R^[LR1+2(I-1)]* R^[LR2+2(J-1)]
                            CALL EXPINTEGRAL(2*(LR2+I+J-1),E(IE1)+E(IE2),SVAR0)
                            SVAR1=POTPAR(ISP)%TAILED%PRODRHO%C(I,IE1,IP1)
                            SVAR2=POTPAR(ISP)%TAILED%PRODPOT%C(J,IE2,IP2)
                            VAL1=VAL1+SVAR1*SVAR2*SVAR0
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  END IF
                ENDDO
                DO LR2=ABS(LR1-L3),LR1+L3,2
                  IF(LR2.EQ.L4) THEN                                   
                    CALL LMTO_TAILEDINDEX_T(NX,LNX,LOX,LN1,LN2,LR1,LN3,LR2,IT)
                    CALL LMTO_TAILEDINDEX_S(NX,LNX,LOX,LN4,IS)
                    CALL LMTO_TESTPLOTRADIALGAUSSALONE('TRIPLE',LR2,NIJK,NE,E &
                                        ,POTPAR(ISP)%TAILED%TRIPLE%C(:,:,IT))
                    CALL LMTO_TESTPLOTRADIALGAUSSALONE('SINGLE',L4,NIJK,NE,E &
                                        ,POTPAR(ISP)%TAILED%SINGLE%C(:,:,IS))
                    VAL2=0.D0
                    DO IE1=1,NE
                      DO IE2=1,NE
                        DO I=1,NIJK
                          DO J=1,NIJK
                            CALL EXPINTEGRAL(2*(LR2+I+J-1),E(IE1)+E(IE2),SVAR0)
                            SVAR1=POTPAR(ISP)%TAILED%TRIPLE%C(I,IE1,IT)
                            SVAR2=POTPAR(ISP)%TAILED%SINGLE%C(J,IE2,IS)
                            VAL2=VAL2+SVAR1*SVAR2*SVAR0
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  END IF
                ENDDO
!
!!$DO I=1,NIJK
!!$ WRITE(*,FMT='("PRODRHO",I5,4F20.10)')I,POTPAR(ISP)%TAILED%PRODRHO%C(I,:,IP1)
!!$ENDDO
!!$DO I=1,NIJK
!!$ WRITE(*,FMT='("PRODPOT",I5,4F20.10)')I,POTPAR(ISP)%TAILED%PRODPOT%C(I,:,IP2)
!!$ENDDO
!!$DO I=1,NIJK
!!$ WRITE(*,FMT='("SINGKE",I5,4F20.10)')I,POTPAR(ISP)%TAILED%SINGLE%C(I,:,IS)
!!$ENDDO
!!$DO I=1,NIJK
!!$ WRITE(*,FMT='("TRIPLE",I5,4F20.10)')I,POTPAR(ISP)%TAILED%TRIPLE%C(I,:,IT)
!!$ENDDO
                PRINT*,'IP1,IP2,IS,IT  = ',IP1,IP2,IS,IT
                PRINT*,'LN1-4   = ',LN1,LN2,LN3,LN4
                PRINT*,'L1-4    = ',L1,L2,L3,L4
                PRINT*,'LR1     = ',LR1
                PRINT*,'++++  ',COUNT,VAL1,VAL2,VAL1-VAL2
!IF(COUNT.EQ.2) STOP
!                 ==============================================================
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
        DEALLOCATE(E)
        DEALLOCATE(LOX)
      ENDDO
STOP 'FORCED IN LMTO_TESTTAILEDP'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXPINTEGRAL(N,E,RES)
!     **************************************************************************
!     **************************************************************************
      implicit none
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: E
      LOGICAL(4)            :: TEVEN
      REAL(8)               :: PI
      REAL(8)   ,INTENT(OUT):: RES
      INTEGER(4)            :: K,I
!     **************************************************************************
      K=INT(N/2)
      TEVEN=(2*K.EQ.N)
      PI=4.D0*ATAN(1.D0)
      IF(TEVEN) THEN
        RES=0.5D0*SQRT(PI/E)
        DO I=1,2*K-1,2
          RES=RES*REAL(I,KIND=8)/(2.D0*E)
        ENDDO
      ELSE
        RES=1.D0/(2.D0*E)
        DO I=1,K
          RES=RES*REAL(I,KIND=8)/E
        ENDDO
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILEDINDEX_T(NX,LNX,LOX,LN1,LN2,LR1,LN3,LR2,IT)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LN1,LN2,LR1,LN3,LR2
      INTEGER(4),INTENT(OUT) :: IT
      INTEGER(4)             :: L1,L2,L3,LN1A,LN2A,LN3A,LR1A,LR2A
      INTEGER(4)             :: IS,IP
      INTEGER(4)             :: NPOW2
!     **************************************************************************
      IT=0   !#(TRIPLES)
      IP=0   !#(PRODUCTS)
      IS=0 !#(SINGLES)
      DO LN1A=1,LNX
        L1=LOX(LN1A)
        NPOW2=INT(0.5D0*REAL(NX-L1)+1.000001D0)  !R^(L+2N), N=0,NPOW-1
        IF(NPOW2.LT.1) CYCLE
        IS=IS+1
        DO LN2A=LN1A,LNX
          L2=LOX(LN2A)
          DO LR1A=ABS(L1-L2),L1+L2,2  ! TRIANGLE RULE
            NPOW2=INT(0.5D0*REAL(NX-LR1A)+1.000001D0)  !R^(L+2N), N=0,NPOW2-1
            IF(NPOW2.LT.1) CYCLE
            IP=IP+1
            DO LN3A=1,LNX
              L3=LOX(LN3A)
              DO LR2A=ABS(LR1A-L3),LR1A+L3,2 ! TRIANGLE RULE
                NPOW2=INT(0.5D0*REAL(NX-LR2A)+1.000001D0) !R^(L+2N), N=0,NPOW-1
                IF(NPOW2.LT.1) CYCLE
                IT=IT+1
                IF(LN1A.LT.MIN(LN1,LN2)) CYCLE
                IF(LN2A.LT.MAX(LN1,LN2)) CYCLE
                IF(LR1A.LT.LR1) CYCLE
                IF(LN3A.LT.LN3) CYCLE
                IF(LR2A.LT.LR2) CYCLE
                IF(LN1A.NE.MIN(LN1,LN2).OR.LN2A.NE.MAX(LN1,LN2) &
     &             .OR.LR1A.NE.LR1.OR.LN3A.NE.LN3.OR.LR2A.NE.LR2) THEN
                  CALL ERROR$STOP('LMTO_TAILEDINDEX_T')
                END IF
                RETURN
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END    
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILEDINDEX_P(NX,LNX,LOX,LN1,LN2,LR1,IP)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LN1,LN2,LR1
      INTEGER(4),INTENT(OUT) :: IP
      INTEGER(4)             :: L1,L2,LN1A,LN2A,LR1A
      INTEGER(4)             :: NPOW2
!     **************************************************************************
      IP=0   !#(PRODUCTS)
      DO LN1A=1,LNX
        L1=LOX(LN1A)
        NPOW2=INT(0.5D0*REAL(NX-L1)+1.000001D0)  !R^(L+2N), N=0,NPOW-1
        IF(NPOW2.LT.1) CYCLE
        DO LN2A=LN1A,LNX
          L2=LOX(LN2A)
          DO LR1A=ABS(L1-L2),L1+L2,2  ! TRIANGLE RULE
            NPOW2=INT(0.5D0*REAL(NX-LR1A)+1.000001D0)  !R^(L+2N), N=0,NPOW2-1
            IF(NPOW2.LT.1) CYCLE
            IP=IP+1
            IF(LN1A.LT.MIN(LN1,LN2)) CYCLE
            IF(LN2A.LT.MAX(LN1,LN2)) CYCLE
            IF(LR1A.LT.LR1) CYCLE
            IF(LN1A.NE.MIN(LN1,LN2).OR.LN2A.NE.MAX(LN1,LN2) &
     &         .OR.LR1A.NE.LR1) THEN
              CALL ERROR$I4VAL('LN1',LN1)
              CALL ERROR$I4VAL('LN2',LN2)
              CALL ERROR$I4VAL('LN1A',LN1A)
              CALL ERROR$I4VAL('LN2A',LN2A)
              CALL ERROR$I4VAL('LR1',LR1)
              CALL ERROR$I4VAL('LR1A',LR1A)
              CALL ERROR$STOP('LMTO_TAILEDINDEX_P')
            END IF
            RETURN
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END    
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILEDINDEX_S(NX,LNX,LOX,LN1,IS)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LN1
      INTEGER(4),INTENT(OUT) :: IS
      INTEGER(4)             :: LN1A,L1,NPOW2
!     **************************************************************************
      IS=0
      DO LN1A=1,LNX
        L1=LOX(LN1A)
        NPOW2=INT(0.5D0*REAL(NX-L1)+1.000001D0)  !R^(L+2N), N=0,NPOW-1
        IF(NPOW2.LT.1) CYCLE
        IS=IS+1
        IF(LN1A.LT.LN1) CYCLE
        IF(LN1A.NE.LN1) THEN
          CALL ERROR$I4VAL('LN1',LN1)
          CALL ERROR$I4VAL('LN1A',LN1A)
          CALL ERROR$STOP('LMTO_TAILEDINDEX_S')
        END IF
        RETURN
      ENDDO
      RETURN
      END    
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE LMTO$OVERLAP(NFIL)
!!$!     **************************************************************************      
!!$!     **                                                                      **
!!$!     **                                                                      **
!!$!     **************************************************************************      
!!$      USE LMTO_MODULE
!!$      IMPLICIT NONE
!!$      TYPE OVKJ_TYPE
!!$        REAL(8),POINTER :: KK
!!$        REAL(8),POINTER :: KJBAR
!!$        REAL(8),POINTER :: JBARJBAR
!!$      END TYPE OVKJ_TYPE
!!$      INTEGER(4),INTENT(IN) :: NFIL
!!$      INTEGER(4)            :: NNS
!!$      REAL(8)               :: A(2,2)
!!$      TYPE(OVKJ_TYPE)       :: OVKJ(NSP)
!!$      INTEGER(4)            :: ISVAR
!!$      INTEGER(4)            :: NN,NN1,NN2,LM1,ISP,I
!!$      INTEGER(4)            :: NNS
!!$!     **************************************************************************      
!!$      DO ISP=1,NSP
!!$        ISVAR=(LX(ISP)+1)**2
!!$        ALLOCATE(OVKJ(ISP)%KK(ISVAR))
!!$        ALLOCATE(OVKJ(ISP)%KJBAR(ISVAR))
!!$        ALLOCATE(OVKJ(ISP)%JBARJBAR(ISVAR))
!!$        LM1=0
!!$        DO L=0,LX(ISP)
!!$          A(1,1)=POTPAR(L+1,ISP)%OVUU
!!$          A(1,2)=POTPAR(L+1,ISP)%OVUQ
!!$          A(2,1)=POTPAR(L+1,ISP)%OVUQ
!!$          A(2,2)=POTPAR(L+1,ISP)%OVQQ
!!$          A(:,:)=MATMUL(A,POTPAR(L+1,ISP)%KJBARTOUQ)
!!$          A(:,:)=MATMUL(TRANSPOSE(POTPAR(L+1,ISP)%KJBARTOUQ),A)
!!$          DO IM=1,2*L+1
!!$            LM1=LM1+1
!!$             OVKJ(ISP)%KK(LM1)=A(1,1)
!!$             OVKJ(ISP)%KJBAR(LM1)=A(1,2)
!!$             OVKJ(ISP)%JBARJBAR(LM1)=A(2,2)
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$      NNSS=SIZE(SBAR)
!!$!
!!$!     ===========================================================================
!!$!     == ONSITE TERMS                                                          ==
!!$!     ===========================================================================
!!$      DO NN=1,NNS
!!$        IAT1=OVERLAP(NN)%IAT1
!!$        IAT2=OVERLAP(NN)%IAT2
!!$        IF(IAT1.NE.IAT2) CYCLE
!!$        IT(:)=OVERLAP(NN)%IT(:)
!!$        IF(IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
!!$        ISP1=ISPECIES(IAT1)        
!!$        ISP2=ISPECIES(IAT2)        
!!$        OVERLAP(NN)%MAT(:,:)=0.D0
!!$        DO LM=1,OVERLAP(NN)%N1
!!$          OVERLAP(NN)%MAT(LM,LM)=OVKJ(ISP)%KK(LM)
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ===========================================================================
!!$!     == PAIR TERMS                                                            ==
!!$!     ===========================================================================
!!$      DO NN=1,NNS
!!$        IAT1=OVERLAP(NN)%IAT1
!!$        IAT2=OVERLAP(NN)%IAT2
!!$        IT(:)=OVERLAP(NN)%IT(:)
!!$        ISP1=ISPECIES(IAT1)
!!$        ISP1=ISPECIES(IAT2)
!!$        DO NN1=1,NNSS
!!$          IF(IAT1.NE.SBAR(NN1)%IAT1) CYCLE
!!$          IF(IAT2.NE.SBAR(NN1)%IAT2) CYCLE
!!$          IF(MAXVAL(ABS(IT-SBAR(NN1)%IT)).NE.0) CYCLE
!!$          DO I=1,N1
!!$            DO J=1,N1
!!$              OVERLAP(NN)%MAT(I,J)=OVERLAP(NN)%MAT(I,J) &
!!$     &                       +OVKJ(ISP1)%KJBAR(I)*SBAR(NN1)%MAT(I,J)
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$        DO NN1=1,NNSS
!!$          IF(IAT2.NE.SBAR(NN1)%IAT1) CYCLE
!!$          IF(IAT1.NE.SBAR(NN1)%IAT2) CYCLE
!!$          IF(MAXVAL(ABS(IT+SBAR(NN1)%IT)).NE.0) CYCLE
!!$          DO I=1,N1
!!$            DO J=1,N1
!!$              OVERLAP(NN)%MAT(I,J)=OVERLAP(NN)%MAT(I,J) &
!!$     &                            +SBAR(NN1)%MAT(J,I)*OVKJ(ISP2)%KJBAR(J) 
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ===========================================================================
!!$!     == THREE-CENTER TERMS                                                    ==
!!$!     ===========================================================================
!!$      DO NN=1,NNS
!!$        IAT1=OVERLAP(NN)%IAT1
!!$        IAT2=OVERLAP(NN)%IAT2
!!$        IT(:)=OVERLAP(NN)%IT(:)
!!$        ISP1=ISPECIES(IAT1)
!!$        ISP1=ISPECIES(IAT2)
!!$        DO NN1=1,NNSS
!!$          IF(IAT1.NE.SBAR(NN1)%IAT1) CYCLE
!!$          IT2(:)=SBAR(NN1)%IT
!!$          DO NN2=1,NNSS
!!$            IF(IAT2.NE.SBAR(NN1)%IAT2) CYCLE
!!$            IT3(:)=SBAR(NN2)%IT
!!$            IF(MAXVAL(ABS(IT2-IT3-IT)).NE.0) CYCLE
!!$            DO I=1,N1
!!$               DO J=1,N2
!!$                 DO K=1,N3
!!$                   OVERLAP(NN)%MAT(I,J)=OVERLAP(NN)%MAT(I,J) &
!!$     &               +SBAR(NN1)%MAT(I,K)*OVKJ(ISP3)%JBARJBAR(K)*SBAR(NN3)%MAT(J,K)
!!$                 ENDDO
!!$              ENDDO
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$      RETURN
!!$      END
