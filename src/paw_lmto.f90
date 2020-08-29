! TAILED_TYPE REFERS TO THE ORIGINAL TAILED TYPE BEFORE
!COMMIT   C4E8FD5B386D229804A4B4E068F50BD4337A8EB0
!FROM APRIL 4, 2014 AT 13:25:01 GMT+2
! 1) TYPE POTPAR1_TYPE AND TAILED1_TYPE SHALL REPLACE THE 
!    TYPES POTPAR_TYPE AND TAILED.
!
! THE SUBROUTINES
!      LMTO_MAKEPOTPAR1()
!      LMTO_MAKETAILEDPARTIALWAVES_WITHPOTPAR1()
!      LMTO_MAKETAILEDMATRIXELEMENTS_WITHPOTPAR1()
!
!*******************************************************************************
!*******************************************************************************
!**
!**  LMTO$MAKESTRUCTURECONSTANTS()
!**       ->LMTO_INITIALIZE()
!**            ->LMTO_MAKEPOTPAR1()
!**            ->LMTO_MAKETAILEDPARTIALWAVES_WITHPOTPAR1()
!**            ->LMTO_MAKETAILEDMATRIXELEMENTS_WITHPOTPAR1()
!**            ->LMTO_TAILEDGAUSSFIT()
!**            -> LMTO_OFFXINT()  (TIME CONSUMING, INTERNALLY PARALLEL)
!**       -> LMTO$CLUSTERSTRUCTURECONSTANTS (PARALLEL)
!**       -> LMTO_PCHI
!**       -> LMTO_TAILEDEXPANSION
!**  LMTO$PROJTONTBO_NEW (TRANSFORM PROJECTIONS TO AND FROM LOCAL ORBITALS)
!**  LMTO$DEDF(NB,NKPTL,NSPIN,DEIG)
!**  LMTO$ETOT(NB,NKPTL,NSPIN,DEIG) 
!** 
!**
!** PARALLELIZATION MODEL:
!** 
!*******************************************************************************
!*******************************************************************************
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
  REAL(8)     :: TAILEDLAMBDA1  ! LARGER DECAY CONSTANT FOR THE NTBO TAILS
  REAL(8)     :: TAILEDLAMBDA2  ! SMALLER DECAY CONSTANT FOR THE NTBO TAILS
  REAL(8)     :: RAUG       ! AUGMENTATION RADIUS
  REAL(8)     :: RTAIL      ! RADIUS FOR MATCHING TAILS
  REAL(8)     :: RTAILCUT   ! RADIUS FOR CUTTING TAILS OFF
  INTEGER(4),POINTER :: NORBOFL(:) ! #(LOCAL ORBITALS PER ANGULAR MOMENTUM)
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
TYPE ORBITALSPHHARM_TYPE
  INTEGER(4)         :: GID
  INTEGER(4)         :: NR
  INTEGER(4)         :: LMX
  INTEGER(4)         :: NORB
  REAL(8)   ,POINTER :: F(:,:,:) !(NR,LM,IORB)
END TYPE ORBITALSPHHARM_TYPE
!
TYPE TAILEDEXPANSION_TYPE
  INTEGER(4)         :: LMNXT    !SIZE OF LEFT INDEX, #(TAILED ORBITALS)
  INTEGER(4)         :: LMNXH    !SIZE OF RIGHT INDEX, #(TIGHT-BINDING ORBITALS)
  REAL(8)   ,POINTER :: MAT(:,:) !(LMNXT,LMNXH) EXPANSION COEFFICIENTS
END TYPE TAILEDEXPANSION_TYPE

TYPE TAILED1_TYPE
  INTEGER(4)         :: GID
  INTEGER(4)         :: LNX
  INTEGER(4)         :: LMNX
  INTEGER(4),POINTER :: LOX(:)           ! (LNX)
  REAL(8)   ,POINTER :: AEF(:,:)         ! (NR,LNX)
  REAL(8)   ,POINTER :: PSF(:,:)         ! (NR,LNX)
  REAL(8)   ,POINTER :: NLF(:,:)         ! (NR,LNX)
  REAL(8)   ,POINTER :: U(:,:,:,:)       ! (LMNX,LMNX,LMNX,LMNX)
  REAL(8)   ,POINTER :: OVERLAP(:,:) !(LMNX,LMNX) OVERLAP MATRIX ELEMENTS
  REAL(8)   ,POINTER :: QLN(:,:,:) !(2,LNX,LNX) MONO- AND DIPOLE MATRIX ELEMENTS
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSNLF
END TYPE TAILED1_TYPE

!
!== HOLDS THE POTENTIAL PARAMETER FOR ONE ATOM TYPE ============================
TYPE POTPAR1_TYPE
  ! THE NUMBER OF ACTIVE HEAD FUNCTIONS IS NPHI. 
  ! THE NUMBER OF TAIL FUNCTIONS IS LX.
  ! EACH HEAD FUNCTION IS RELATED TO A PARTIAL WAVE IDENTIFIED BY LN.
  REAL(8)            :: RAD             ! MATCHING RADIUS
  INTEGER(4)         :: NHEAD           ! #(HEAD FUNCTIONS)
  INTEGER(4)         :: NTAIL           ! #(TAIL FUNCTIONS)
  INTEGER(4),POINTER :: ITAIL(:)        !(NHEAD) POINTER TO TAIL FUNCTION
  INTEGER(4),POINTER :: LOFH(:)         !(NHEAD) MAIN ANGULAR MOMENTUM
  INTEGER(4),POINTER :: LNOFH(:)        !(NHEAD) PARTIAL WAVE ID
  REAL(8)   ,POINTER :: KTOPHI(:)       !(NHEAD) |K> = |PHI>    * KTOPHI
  REAL(8)   ,POINTER :: KTOPHIDOT(:)    !(NHEAD)     + |PHIDOT> * KTOPHIDOT
  REAL(8)   ,POINTER :: PHIDOTPROJ(:)   !(LNX)   <P(LN)|PHIDOT(ITAIL(IHEAD))>
  INTEGER(4),POINTER :: LOFT(:)         !(NTAIL)  MAIN ANGULAR MOMENTUM
  INTEGER(4),POINTER :: LNOFT(:)        !(NTAIL)  PARTIAL WAVE ID FOR PHIDOT
  REAL(8)   ,POINTER :: QBAR(:)         !(NTAIL)  |JBAR>=|J>-|K>QBAR
  REAL(8)   ,POINTER :: JBARTOPHIDOT(:) !(NTAIL)|JBAR> = |PHIBARDOT> JBARTOPHIDOT
  REAL(8)   ,POINTER :: PROK(:,:)       !(LNX,NHEAD) <P|K_AUG>
  REAL(8)   ,POINTER :: PROJBAR(:,:)    !(LNX,NTAIL) <P|JBAR_AUG> 
  REAL(8)   ,POINTER :: PHIOV(:,:)      !(LNX,LNX) <AEPHI|THETA_OMEGA|AEPHI>
  TYPE(TAILED1_TYPE) :: TAILED
END TYPE POTPAR1_TYPE
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
LOGICAL(4)            :: TOFFSITE=.FALSE.  !INCLUDE OFFSITE EXCHANGE
LOGICAL(4)            :: TDROP=.FALSE. ! WRITE THE WAVE FUNCTIONS TO FILE
LOGICAL(4)            :: TPICK=.FALSE. ! REAL HAMILTON CORRECTION FROM FILE
REAL(8)               :: K2=-0.25D0    ! 0.5*K2 IS THE KINETIC ENERGY
!REAL(8)               :: K2=0.D0    ! 0.5*K2 IS THE KINETIC ENERGY
REAL(8)               :: RCSCALE=2.D0  !RADIUS SCALE FACTOR FOR NEIGHBORLIST
! RCSCALE=5. IS GOOD FOR THE HUBBARD MODEL WITH LATTICE CONSTANT=3\AA
!REAL(8)               :: RCSCALE=5.D0  !RADIUS SCALE FACTOR FOR NEIGHBORLIST
!         RCSCALE TIMES THE SUM OF COVALENT RADII DEFINES CUTOFF FOR NEIGBORLIST
REAL(8)               :: HFWEIGHT=0.25D0
!
!===============================================================================
!== VARIABLE SECTION                                                          ==
!===============================================================================
LOGICAL(4)              :: TINI=.FALSE.
LOGICAL(4)              :: THTBC=.FALSE. ! HTBC CALCULATED
LOGICAL(4)              :: TOTBC=.FALSE. ! OTBC CALCULATED
CHARACTER(32)           :: MODUS='NONE'
INTEGER(4)              :: NSP=-1
INTEGER(4)              :: ISPSELECTOR=-1 ! USED ONLY FOR HYBRIDSETTING
TYPE(HYBRIDSETTING_TYPE),ALLOCATABLE :: HYBRIDSETTING(:)
INTEGER(4),ALLOCATABLE  :: LNX(:)      !(NSP)
INTEGER(4),ALLOCATABLE  :: LOX(:,:)    !(LNXX,NSP)
INTEGER(4),ALLOCATABLE  :: ISPECIES(:) !(NAT)
REAL(8)   ,ALLOCATABLE  :: ORBRAD(:,:) !(LXX+1,NAT) NODE-POSITION OF THE ORBITAL
TYPE(POTPAR1_TYPE)    ,ALLOCATABLE :: POTPAR1(:) !POTENTIAL PARAMETERS (NEW)
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
LOGICAL(4)                          :: TCTE=.TRUE. 
TYPE(TAILEDEXPANSION_TYPE), POINTER :: CTE(:)=>NULL() ! TAILED EXPANSION
!===============================================================================
!=====  STRUCTURE DEPENDENT DATA  ==============================================
!===============================================================================
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: SBAR_NEW(:)!(NNS) SCREENED STRUCT. CONST.
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: PCHI(:)    !(NNS) <PTILDE|CHITILDE>
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DENMAT_NEW(:)  !(NND) DENSITY MATRIX
! DERIVATIVE OF THE ENERGY WITH RESPECT TO THE DENSITY MATRIX
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: HAMIL_NEW(:)   !(NND) DERIVATIVE OF ENERGY
! DERIVATIGVE OF THE ENERGY WITH RESPECT TO THE INVERSE OVERLAP MATRIX
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DEDOI(:)   !(NND) DERIVATIVE OF ENERGY
! OVERLAP IS ONLY USED FOR TESTDENMAT
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: OVERLAP(:) !(NNS) OVERLAP MATRIX ONLY MAIN
!== INVOVERLAP CALCULATED AS SUM_N <PI|PSI_N><PSI_N|PI>
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: INVOVERLAP(:) !(NNS) OVERLAP MATRIX ONLY MAIN
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
      ELSE IF(ID.EQ.'THTBC') THEN
        THTBC=VAL
!
!     ==========================================================================
!     == IF ACTIVE THE HYBRID CONTRIBTIONS ON THIS ATOM ARE CONSIDERED        ==
!     ==========================================================================
      ELSE IF(ID.EQ.'ACTIVE') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%ACTIVE=VAL
!
      ELSE IF(ID.EQ.'COREVALENCE') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TCV=VAL
!
      ELSE IF(ID.EQ.'FOCKSETUP') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TFOCKSETUP=VAL
!
      ELSE IF(ID.EQ.'NDDO') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TNDDO=VAL
!
      ELSE IF(ID.EQ.'31') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%T31=VAL
!
      ELSE IF(ID.EQ.'BONDX') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
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
      USE LMTO_MODULE, ONLY : TON &
     &                       ,THTBC &
     &                       ,TOFFSITE &
     &                       ,HYBRIDSETTING &
     &                       ,ISPSELECTOR
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      LOGICAL(4)  ,INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'ON') THEN
        VAL=TON
!
      ELSE IF(ID.EQ.'THTBC') THEN
        VAL=THTBC
!
      ELSE IF(ID.EQ.'TOFFSITE') THEN
        VAL=TOFFSITE
!
      ELSE IF(ID.EQ.'FOCKSETUP') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
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
      SUBROUTINE LMTO$GETI4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES & !(NAT)
     &                       ,POTPAR1   !(NSP)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(OUT) :: VAL
      INTEGER(4)               :: IAT,ISP
!     **************************************************************************
      IF(ID.EQ.'NLOCORB') THEN
        VAL=0
        DO IAT=1,SIZE(ISPECIES)
          ISP=ISPECIES(IAT)
          VAL=VAL+SUM(2*POTPAR1(ISP)%LOFH+1)
        ENDDO
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SETI4A(ID,LEN,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES & !(NAT)
     &                       ,POTPAR1 &  !(NSP)
     &                       ,ISPSELECTOR &
     &                       ,HYBRIDSETTING
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      INTEGER(4)  ,INTENT(IN)  :: VAL(LEN)
      INTEGER(4)               :: IAT,ISP
!     **************************************************************************
      IF(ID.EQ.'NORBOFL') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETI4A')
        END IF
        ALLOCATE(HYBRIDSETTING(ISPSELECTOR)%NORBOFL(LEN))
        HYBRIDSETTING(ISPSELECTOR)%NORBOFL=VAL
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO$SETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SETI4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : HYBRIDSETTING &
     &                       ,ISPSELECTOR &
     &                       ,NSP
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
      INTEGER(4)              :: I
!     **************************************************************************
      IF(ID.EQ.'ISP') THEN
!       ========================================================================
!       == THE ATOM TYPE SELECTOR ISP IS USED ONLY FOR SETTINGS               ==
!       ========================================================================
        ISPSELECTOR=VAL
        CALL SETUP$GETI4('NSP',NSP)
        IF(NSP.LE.0) THEN
          CALL ERROR$MSG('NUMBER OF ATOM TYPES UNKNOWN BY SETUP OBJECT')
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
          HYBRIDSETTING(:)%ACTIVE    =.FALSE.
          HYBRIDSETTING(:)%TCV       =.TRUE.
          HYBRIDSETTING(:)%TFOCKSETUP=.FALSE.
          HYBRIDSETTING(:)%LHFWEIGHT =-1.D0
          HYBRIDSETTING(:)%TNDDO     =.FALSE.
          HYBRIDSETTING(:)%T31       =.FALSE.
          HYBRIDSETTING(:)%TBONDX    =.FALSE.
          HYBRIDSETTING(:)%TAILEDLAMBDA1=4.D0
          HYBRIDSETTING(:)%TAILEDLAMBDA2=2.D0
          HYBRIDSETTING(:)%RAUG      =-1.D0
          HYBRIDSETTING(:)%RTAIL     =-1.D0
          DO I=1,NSP
            NULLIFY(HYBRIDSETTING(I)%NORBOFL)
          ENDDO
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
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
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
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%LHFWEIGHT=VAL
!
      ELSE IF(ID.EQ.'TAILLAMBDA2') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TAILEDLAMBDA2=VAL
!
      ELSE IF(ID.EQ.'RAUG') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%RAUG=VAL
!
      ELSE IF(ID.EQ.'RTAIL') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%RTAIL=VAL
!
      ELSE IF(ID.EQ.'RTAILCUT') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%RTAILCUT=VAL
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
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
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
      SUBROUTINE LMTO$SETCH(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : MODUS
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'MODUS') THEN
        MODUS=VAL
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO$SETCH')
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
     &                       ,TOFFSITE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: THISTASK,NTASKS
      CHARACTER(32)         :: ID
      INTEGER(4)            :: ISP,I
      LOGICAL(4)            :: TCHK
      INTEGER(4)            :: LNX1
      INTEGER(4),ALLOCATABLE :: LOX1(:)
      REAL(8)                :: AEZ,RCOV
      REAL(8)                :: LHFWEIGHT
!     **************************************************************************
      IF(.NOT.TON) RETURN
      CALL LMTO_INITIALIZE()
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      CALL REPORT$TITLE(NFIL,'LMTO OBJECT:  GENERIC VARIABLES')
      CALL REPORT$R8VAL(NFIL,'OFFSITE EXCHANGE ADMIXTURE ',HFWEIGHT,' ')
      CALL REPORT$R8VAL(NFIL,'RANGESCALE ',RCSCALE,'*(RCOV(A)+RCOV(B))')
      CALL REPORT$R8VAL(NFIL,'K2 ',K2,'A.U.')
      IF(TOFFSITE) THEN
        CALL REPORT$CHVAL(NFIL,'OFFSITE MATRIX ELEMENTS ARE','INCLUDED')
      ELSE
        CALL REPORT$CHVAL(NFIL,'OFFSITE MATRIX ELEMENTS ARE','EXLUDED')
      END IF
!
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETCH('ID',ID)
        CALL SETUP$GETI4('LNX',LNX1)
        ALLOCATE(LOX1(LNX1))
        CALL SETUP$GETI4A('LOX',LNX1,LOX1)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
        WRITE(NFIL,*)
        CALL REPORT$TITLE(NFIL,'LMTO OBJECT: '//TRIM(ID))
        LHFWEIGHT=HYBRIDSETTING(ISP)%LHFWEIGHT
        IF(LHFWEIGHT.LT.0.D0) LHFWEIGHT=HFWEIGHT
        CALL REPORT$R8VAL(NFIL,'LOCAL EXCHANGE ADMIXTURE ',LHFWEIGHT,' ')
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
        CALL REPORT$R8VAL(NFIL,'AUGMENTATION RADIUS ' &
     &                        ,HYBRIDSETTING(ISP)%RAUG/RCOV,'*R(COV) ')
        CALL REPORT$R8VAL(NFIL,'TAIL MATCHING RADIUS ' &
     &                        ,HYBRIDSETTING(ISP)%RTAIL/RCOV,'*R(COV)')
        CALL REPORT$R8VAL(NFIL,'TRUNCATION RADIUS FOR EXPONENTIAL TAILS' &
     &                        ,HYBRIDSETTING(ISP)%RTAILCUT/RCOV,'*R(COV)')
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
      USE LMTO_MODULE, ONLY : K2,RCSCALE &
     &                       ,POTPAR1 &
     &                       ,ISPECIES,NSP &
     &                       ,SBAR_NEW
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
      REAL(8)   ,ALLOCATABLE :: QBAR(:)    !(N) SCREENING PARAMETER
      REAL(8)   ,ALLOCATABLE :: SBAR1(:,:) !
      REAL(8)                :: SVAR
      INTEGER(4)             :: IAT,IAT1,IAT2,ISP,ISP1,ISP2,LM1X,LM2X 
      INTEGER(4)             :: L,NN,NN1,NN2,NN0,I,I1,I2,LX,L1,L2
      INTEGER(4)             :: I11,I12,I21,I22
      INTEGER(4)             :: IT,IH,LMN1,LMN2
      INTEGER(4)             :: NHEAD,NTAIL
      LOGICAL(4)             :: TCHK
      REAL(8)   ,ALLOCATABLE :: RC(:)
      INTEGER(4)             :: LMX
      INTEGER(4)             :: NTASKS,THISTASK,FROMTASK
LOGICAL(4),SAVE:: TFIRST=.TRUE.
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
                              CALL TRACE$PUSH('LMTO$MAKESTRUCTURECONSTANTS')
                              CALL TIMING$CLOCKON('LMTO STRUCTURECONSTANTS')
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
      ALLOCATE(SBAR_NEW(NNS))
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
          LX=-1
          IF(POTPAR1(ISP)%NHEAD.NE.0)LX=MAX(LX,MAXVAL(POTPAR1(ISP)%LOFH))
          IF(POTPAR1(ISP)%NTAIL.NE.0)LX=MAX(LX,MAXVAL(POTPAR1(ISP)%LOFT))
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
        DO NN=NN1,NN2 !LOOP OVER ALL ATOMS ON THIS CLUSTER
          IAT2=NNLIST(2,NN)
          ISP=ISPECIES(IAT2)
          DO IT=1,POTPAR1(ISP)%NTAIL
            L=POTPAR1(ISP)%LOFT(IT)
            QBAR(I+L**2+1:I+(L+1)**2)=POTPAR1(ISP)%QBAR(IT)
          ENDDO
          LX=LX1(NN-NN1+1)
          I=I+(LX+1)**2
        ENDDO
!
!       ========================================================================
!       == DETERMINE STRUCTURE CONSTANTS                                      ==
!       == HERE, THE STRUCTURE CONSTANTS USE A COMPLETE SET OF ANGULAR MOMENTA==
!       == UP TO A MAXIMUM ANGULAR MOMENTUM. NOT ALL WILL BE USED LATER ON    ==
!       ========================================================================
!NORB=(LX1(1)+1)**2
!N=SUM((LX(1)+1)**2)
        ALLOCATE(SBAR1(NORB,N))
        IF(MOD(IAT1-1,NTASKS).EQ.THISTASK-1) THEN
          CALL TIMING$CLOCKON('STRUCTURE CONSTANTS')
          CALL LMTO$CLUSTERSTRUCTURECONSTANTS(K2,NAT1,RPOS,LX1,QBAR &
       &                                                          ,NORB,N,SBAR1)
          CALL TIMING$CLOCKOFF('STRUCTURE CONSTANTS')
        ELSE
          SBAR1=0.D0
        END IF
!
!       ========================================================================
!       == MAP ONTO SBAR                                                      ==
!       ========================================================================
        I=0
        DO NN=NN1,NN2
          IAT2=NNLIST(2,NN)
          SBAR_NEW(NN)%IAT1=NNLIST(1,NN)
          SBAR_NEW(NN)%IAT2=NNLIST(2,NN)
          SBAR_NEW(NN)%IT(:)=NNLIST(3:5,NN)
          ISP1=ISPECIES(IAT1)
          ISP2=ISPECIES(IAT2)
          NTAIL=POTPAR1(ISP2)%NTAIL
          NHEAD=POTPAR1(ISP1)%NHEAD
          LM1X=SUM(2*POTPAR1(ISP1)%LOFH+1)
          LM2X=SUM(2*POTPAR1(ISP2)%LOFT+1)
          SBAR_NEW(NN)%N1=LM1X
          SBAR_NEW(NN)%N2=LM2X
          ALLOCATE(SBAR_NEW(NN)%MAT(LM1X,LM2X))
          LMN2=0
          DO IT=1,POTPAR1(ISP2)%NTAIL
            L2=POTPAR1(ISP2)%LOFT(IT)
            I21=I+L2**2+1
            I22=I+(L2+1)**2
            LMN1=0
            DO IH=1,POTPAR1(ISP1)%NHEAD
              L1=POTPAR1(ISP1)%LOFH(IH)
              I11=L1**2+1
              I12=(L1+1)**2
              SBAR_NEW(NN)%MAT(LMN1+1:LMN1+2*L1+1,LMN2+1:LMN2+2*L2+1) &
      &                        =SBAR1(I11:I12,I21:I22)
              LMN1=LMN1+2*L1+1
            ENDDO
            LMN2=LMN2+2*L2+1
          ENDDO
          LX=LX1(NN-NN1+1)
          I=I+(LX+1)**2
        ENDDO
        DEALLOCATE(LX1)
        DEALLOCATE(RPOS)
        DEALLOCATE(QBAR)
        DEALLOCATE(SBAR1)
      ENDDO
!
!     == DISTRIBUTE STRUCTURE CONSTANTS
      DO NN=1,NNS
        IAT=SBAR_NEW(NN)%IAT1
        FROMTASK=1+MOD(IAT-1,NTASKS)
        CALL MPE$BROADCAST('MONOMER',FROMTASK,SBAR_NEW(NN)%MAT)
      ENDDO
!!$CALL LMTO$REPORTPERIODICMAT(6,'STRUCTURE CONSTANTS',NNS,SBAR_NEW)
!!$STOP 'FORCED'
!
!     ==========================================================================
!     == REPORT NEIGHBORLIST FOR STRUCTURE CONSTANTS                          ==
!     ==========================================================================
!!$      WRITE(*,FMT='(82("="),T5,"  NEIGHBORLIST FOR SBAR  ")') 
!!$      DO NN=1,NNS
!!$        WRITE(*,FMT='("IAT1=",I5," IAT2=",I5," IT=",3I3)') &
!!$     &          SBAR_NEW(NN)%IAT1,SBAR_NEW(NN)%IAT2,SBAR_NEW(NN)%IT
!!$      ENDDO
!
!     ==========================================================================
!     == CALCULATE <PRO|CHI> FROM STRUCTURE CONSTANTS                         ==
!     ==========================================================================
      CALL LMTO_PCHI()
!
!     ==========================================================================
!     == CALCULATE CTE, THE EXPANSION COEFFICIENTS IN TAILED COMPONENTS       ==
!     ==========================================================================
!IF(TFIRST) THEN
      CALL LMTO_TAILEDEXPANSION()
!  TFIRST=.FALSE.
!END IF
!
!     ==========================================================================
!     == PLOT WAVE FUNCTIONS FOR TESTING (USUALLY INACTIVE)                   ==
!     ==========================================================================
      CALL LMTO_PLOTORBS()
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
      USE LMTO_MODULE, ONLY : SBAR_NEW &
     &                       ,PCHI &
     &                       ,DENMAT=>DENMAT_NEW &
     &                       ,HAMIL=>HAMIL_NEW 
      IMPLICIT NONE
      INTEGER(4)     :: NN
      INTEGER(4)     :: I
!     **************************************************************************
!
!     == CLEAN STRUCTURE CONSTANTS =============================================
      IF(ALLOCATED(SBAR_NEW)) THEN
        NN=SIZE(SBAR_NEW)
        DO I=1,NN
          DEALLOCATE(SBAR_NEW(I)%MAT)
        ENDDO
        DEALLOCATE(SBAR_NEW)
      END IF
!!$!
!!$!     == CLEAN STRUCTURE CONSTANTS =============================================
!!$      IF(ALLOCATED(SBAR)) THEN
!!$        NN=SIZE(SBAR)
!!$        DO I=1,NN
!!$          DEALLOCATE(SBAR(I)%MAT)
!!$        ENDDO
!!$        DEALLOCATE(SBAR)
!!$      END IF
!
!     == CLEAN PROJECTIONS OF LOCAL ORBITALS ===================================
      IF(ALLOCATED(PCHI)) THEN
        NN=SIZE(PCHI)
        DO I=1,NN
          DEALLOCATE(PCHI(I)%MAT)
        ENDDO
        DEALLOCATE(PCHI)
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
      USE LMTO_MODULE, ONLY : TINI &
     &                       ,TON &
     &                       ,TOFFSITE 
      IMPLICIT NONE
      INTEGER(4) :: NAT,ISP
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
                          CALL TRACE$PUSH('LMTO_INITIALIZE')
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
      CALL LMTO_CONSOLIDATETAILEDPARMS()
      CALL LMTO_MAKEPOTPAR1()
!
!     ==========================================================================
!     == ATTACH EXPONENTIAL TAILS TO AUGMENTED HANKEL AND BESSEL FUNCTIONS    ==
!     ==========================================================================
      CALL LMTO_MAKETAILEDPARTIALWAVES_WITHPOTPAR1()
      CALL LMTO_MAKETAILEDMATRIXELEMENTS_WITHPOTPAR1()
      IF(.NOT.TON) RETURN
!
!     ==========================================================================
!     ==  CONSTRUCT OFFSITE INTEGRALS OF TAILED ORBITALS                      ==
!     ==========================================================================
      IF(TOFFSITE) THEN 
                            CALL TIMING$CLOCKON('OFFSITE U-TENSOR')
!       == GAUSSIAN FIT OF TAILED ORBITAL FOR ABAB-TYPE U-TENSOR ===============
                            CALL TIMING$CLOCKON('OFFS.U. GAUSSFIT')
        CALL LMTO_TAILEDGAUSSFIT()
                            CALL TIMING$CLOCKOFF('OFFS.U. GAUSSFIT')
!
!       == OFF-SITE MATRIX ELEMENTS OF OVERLAP MATRIX AND U-TENSOR =============
!       == OVERLAP MATRIX ELEMENTS ARE ALWAYS COMPUTED, U-TENSOR ELEMENTS ======
!       == ONLY WHEN REQUESTED =================================================
                            CALL TIMING$CLOCKON('OFFS.U. INTEGRALS')
        CALL LMTO_OFFXINT()
                            CALL TIMING$CLOCKOFF('OFFS.U. INTEGRALS')
                            CALL TIMING$CLOCKOFF('OFFSITE U-TENSOR')
      END IF
                          CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NULLIFYOFFSITEX()
!     **  THE CODE WILL CHECK IF THE COMPONENTS OF OFFSITEX HAVE BEEN         **
!     **  ASSOCIATED. NULLIFYING IS REQUIRED TO OBTAINED A DEFINED STATUS     **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : NSP,OFFSITEX
      IMPLICIT NONE
      INTEGER(4) :: ISP1,ISP2
!     **************************************************************************
      DO ISP1=1,NSP
        DO ISP2=1,NSP
          NULLIFY(OFFSITEX(ISP1,ISP2)%OVERLAP)
          NULLIFY(OFFSITEX(ISP1,ISP2)%X22)
          NULLIFY(OFFSITEX(ISP1,ISP2)%X31)
          NULLIFY(OFFSITEX(ISP1,ISP2)%BONDU)
        ENDDO
      ENDDO
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
     &                      ,SBARLI1
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      INTEGER(4)             :: LX1
      INTEGER(4)             :: IPOS,IAT,ISP,L,LN
!     **************************************************************************
      IF(ALLOCATED(SBARLI1)) RETURN
                              CALL TRACE$PUSH('LMTO$SBARINDICES')
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
      SUBROUTINE LMTO_MAKEPOTPAR1()
!     **************************************************************************
!     **  SIMILAR TO LMTO_MAKEPOTPAR BUT WITH ANOTHER DATA STRUCTURE          **
!     **  ORGANIZED ACCORDING TO HEADS AND TAILS INSTEAD OF PARTIAL WAVES.    **
!     **  IT IS INTENDED THAT THIS ROUTINE REPLACES LMTO_MAKEPOTPAR.          **
!     **                                                                      **
!     **  POTPAR(ISP)%RAD              : ASA RADIUS                           **
!     **  POTPAR(ISP)%PHIDOTPROJ(LN)   : <PS-PRO|PS-PHIBARDOT>                **
!     **  POTPAR(ISP)%QBAR(LN)         : Q-BAR                                **
!     **  POTPAR(ISP)%KTOPHI(LN)       : K -> |PHI>KTOPHI                     **
!     **  POTPAR(ISP)%KTOPHIDOT(LN)    :    + |PHIBARDOT>KTOPHIDOT            **
!     **  POTPAR(ISP)%JBARTOPHIDOT(LN) : JBAR -> |PHIBARDOT>JBARTOPHIDOT      **
!     **                                                                      **
!     **  POTPAR(ISP)%LNOFH(IHEAD)     : POINTS TO THE PARTIAL WAVE AEPHI(LN) **
!     **                                 OF THE HEAD FUNCTION                 **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : K2,POTPAR1,NSP,LNX,LOX &
     &                       ,HYBRIDSETTING
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: LNX1
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)   ,ALLOCATABLE :: NLPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: PRO(:,:)
      LOGICAL(4),ALLOCATABLE :: TORB(:)
      INTEGER(4),ALLOCATABLE :: ISCATT(:)
      LOGICAL(4),PARAMETER   :: TPRINT=.TRUE.
      REAL(8)                :: RAD
      REAL(8)                :: PHIVAL,PHIDER
      REAL(8)                :: PHIDOTVAL,PHIDOTDER
      REAL(8)                :: KVAL,KDER
      REAL(8)                :: JVAL,JDER
      REAL(8)                :: WJPHI,WJPHIDOT,WKPHI,WKPHIDOT,WJBARPHI,WKJ
      REAL(8)                :: WPHIPHIDOT
      REAL(8)                :: QBAR
      REAL(8)                :: SVAR
      INTEGER(4)             :: NHEAD !#(HEAD FUNCTIONS)
      INTEGER(4)             :: NTAIL !#(TAIL FUNCTIONS)
      INTEGER(4)             :: LX    !X(ANGULAR MOMENTUM)
      INTEGER(4)             :: ISP,LN,L,LNOFH,LNOFT,IHEAD,ITAIL
      INTEGER(4)             :: ISVAR
      LOGICAL(4)             :: TCHK
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO_MAKEPOTPAR1')
      IF(.NOT.ALLOCATED(HYBRIDSETTING)) THEN
        CALL ERROR$MSG('HYBRIDSETTING NOT ALLOCATED')
        CALL ERROR$STOP('LMTO_MAKEPOTPAR1')
      END IF

      ALLOCATE(POTPAR1(NSP))
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
        ALLOCATE(AUX(NR))
        ALLOCATE(ISCATT(LNX1))
        CALL SETUP$GETI4A('ISCATT',LNX1,ISCATT)

!       == MATCHING RADIUS =====================================================
        IF(HYBRIDSETTING(ISP)%RAUG.GE.0.D0) THEN
          RAD=HYBRIDSETTING(ISP)%RAUG
        ELSE
          CALL SETUP$GETR8('AEZ',SVAR)
          CALL PERIODICTABLE$GET(SVAR,'R(COV)',RAD)
        END IF
        POTPAR1(ISP)%RAD=RAD
!
!       == SELECTION OF LOCAL ORBITALS CONSIDERED IN THE U-TENSOR ==============
        ALLOCATE(TORB(LNX1))
        TORB=.FALSE.
        IF(ASSOCIATED(HYBRIDSETTING(ISP)%NORBOFL)) THEN
          DO L=0,SIZE(HYBRIDSETTING(ISP)%NORBOFL)-1  
            ISVAR=HYBRIDSETTING(ISP)%NORBOFL(L+1)
            DO LN=1,LNX1
              IF(LOX(LN,ISP).NE.L) CYCLE
              IF(ISVAR.EQ.0) EXIT  ! NO MORE ORBITALS REQUESTED
              TORB(LN)=.TRUE.
              ISVAR=ISVAR-1
            ENDDO  ! END OF LN-LOOP
            IF(ISVAR.GT.0) THEN
              CALL ERROR$MSG('NUMBER OF NTBOS EXCEEDS NUMBER OF PARTIAL WAVES.')
              CALL ERROR$MSG('INCREASE !STRUCTURE!SPECIES:NPRO OR')
              CALL ERROR$MSG('DECREASE !STRUCTURE!SPECIES!NTBO:NOFL')
              CALL ERROR$I4VAL('ISP',ISP)
              CALL ERROR$I4VAL('L',L)
              CALL ERROR$I4VAL('NOFL',HYBRIDSETTING(ISP)%NORBOFL(L+1))
              CALL ERROR$STOP('LMTO_MAKEPOTPAR1')
            END IF
          ENDDO    ! END OF L-LOOP
        END IF
!
!       == PARTIAL WAVES AND PROJECTORS ========================================
        ALLOCATE(NLPHI(NR,LNX1))
        ALLOCATE(NLPHIDOT(NR,LNX1))
        ALLOCATE(PSPHI(NR,LNX1))
        ALLOCATE(PSPHIDOT(NR,LNX1))
        ALLOCATE(PRO(NR,LNX1))
        CALL SETUP$GETR8A('QPHI',NR*LNX1,NLPHI)
        CALL SETUP$GETR8A('QPHIDOT',NR*LNX1,NLPHIDOT)
        CALL SETUP$GETR8A('PSPHI',NR*LNX1,PSPHI)
        CALL SETUP$GETR8A('PSPHIDOT',NR*LNX1,PSPHIDOT)
        CALL SETUP$GETR8A('PRO',NR*LNX1,PRO)
!
!       ========================================================================
!       ==  COUNT NUMBER OF ACTIVE ORBITALS
!       ========================================================================
        LX=-1
        NHEAD=0
        DO LN=1,LNX1
          IF(.NOT.TORB(LN))CYCLE
          NHEAD=NHEAD+1
          LX=MAX(LX,LOX(LN,ISP))
        ENDDO
        POTPAR1(ISP)%NHEAD=NHEAD
        ALLOCATE(POTPAR1(ISP)%LNOFH(NHEAD))
        ALLOCATE(POTPAR1(ISP)%LOFH(NHEAD))
        ALLOCATE(POTPAR1(ISP)%ITAIL(NHEAD))
        ALLOCATE(POTPAR1(ISP)%KTOPHI(NHEAD))
        ALLOCATE(POTPAR1(ISP)%KTOPHIDOT(NHEAD))
        ALLOCATE(POTPAR1(ISP)%PROK(LNX1,NHEAD))
        ALLOCATE(POTPAR1(ISP)%PHIDOTPROJ(NHEAD))
        POTPAR1(ISP)%PROK(:,:)=0.D0
        IHEAD=0
        DO LN=1,LNX1
          IF(.NOT.TORB(LN))CYCLE
          IHEAD=IHEAD+1
          POTPAR1(ISP)%LNOFH(IHEAD)=LN
          POTPAR1(ISP)%LOFH(IHEAD)=LOX(LN,ISP)
        ENDDO
!       == COUNT NUMBER OF TAIL FUNCTIONS ======================================
        NTAIL=0
        LX=MAXVAL(LOX(:LNX(ISP),ISP))
        DO L=0,LX
          DO IHEAD=1,NHEAD
            IF(POTPAR1(ISP)%LOFH(IHEAD).NE.L) CYCLE
            NTAIL=NTAIL+1
            EXIT
          ENDDO
        ENDDO
        POTPAR1(ISP)%NTAIL=NTAIL
        ALLOCATE(POTPAR1(ISP)%LOFT(NTAIL))
        ALLOCATE(POTPAR1(ISP)%LNOFT(NTAIL))
        ALLOCATE(POTPAR1(ISP)%QBAR(NTAIL))
        ALLOCATE(POTPAR1(ISP)%JBARTOPHIDOT(NTAIL))
        ALLOCATE(POTPAR1(ISP)%PROJBAR(LNX1,NTAIL))
        POTPAR1(ISP)%PROJBAR(:,:)=0.D0
!       == CONNECT EACH TAIL FUNCTION TO A SPECIFIC HEAD FUNCTION FROM WHICH ===
!       == THE SCATTERING FUNCTION IS TAKEN TO DEFINE THE SCREENING CHARGE QBAR.
!       == THE SELECTED HEAD FUNCTION IS THE LAST ONE FOR THE GIVEN L-VALUE.====
!       == IT IS IDENTIFIED BY LNOFT
        ITAIL=0
        DO L=0,LX
          ITAIL=ITAIL+1
          TCHK=.FALSE.
          DO IHEAD=1,NHEAD
            IF(POTPAR1(ISP)%LOFH(IHEAD).NE.L) CYCLE
            TCHK=.TRUE.
            POTPAR1(ISP)%LOFT(ITAIL)=L
            POTPAR1(ISP)%LNOFT(ITAIL)=POTPAR1(ISP)%LNOFH(IHEAD)  !PLACEHOLDER
          ENDDO
          IF(.NOT.TCHK) ITAIL=ITAIL-1  ! NO TAIL FORT THIS L/ DO NOT COUNT UP 
        ENDDO
!       == LINK CORRESPONDING PHIDOT TO EACH PHI ===============================
        DO IHEAD=1,NHEAD
          L=POTPAR1(ISP)%LOFH(IHEAD)
          DO ITAIL=1,NTAIL
            IF(POTPAR1(ISP)%LOFT(ITAIL).NE.L) CYCLE
            POTPAR1(ISP)%ITAIL(IHEAD)=ITAIL
            EXIT
          ENDDO
        ENDDO
!
!       ========================================================================
!       == LNDOT POINTS TO THE PARTIAL WAVE WITH THE PHIDOT FUNCTION          ==
!       ========================================================================
!       == ISCATT=-1 FOR HIGHEST SEMI-CORE STATE WITH THIS L
!       == ISCATT=0 FOR HIGHEST VALENCE STATE WITH THIS L
!       == ISCATT=1 FOR FIRST SCATTERING STATE WITH THIS L
        DO ITAIL=1,NTAIL
          L=POTPAR1(ISP)%LOFT(ITAIL)
          LNOFT=POTPAR1(ISP)%LNOFT(ITAIL)
          DO IHEAD=1,NHEAD
            IF(POTPAR1(ISP)%ITAIL(IHEAD).NE.ITAIL) CYCLE
            LN=POTPAR1(ISP)%LNOFH(IHEAD)
            IF(ISCATT(LN).GT.ISCATT(LNOFT)) LNOFT=LN
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  DETERMINE POTENTIAL PARAMETERS                                    ==
!       ========================================================================
        DO ITAIL=1,NTAIL
          L=POTPAR1(ISP)%LOFT(ITAIL)
          LNOFT=POTPAR1(ISP)%LNOFT(ITAIL)
!         ======================================================================
!         == VALUE AND DERIVATIVE OF PARTIAL WAVES AND ENVELOPE FUNCTIONS     ==
!         == PHIDOT IS PHIBARDOT                                              ==
!         ======================================================================
          CALL RADIAL$VALUE(GID,NR,NLPHIDOT(:,LNOFT),RAD,PHIDOTVAL)
          CALL RADIAL$DERIVATIVE(GID,NR,NLPHIDOT(:,LNOFT),RAD,PHIDOTDER)
          CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
          CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
          WJPHIDOT=JVAL*PHIDOTDER-JDER*PHIDOTVAL
          WKPHIDOT=KVAL*PHIDOTDER-KDER*PHIDOTVAL
          QBAR=WJPHIDOT/WKPHIDOT
!         == |JBAR> = |J> - |K>QBAR ============================================
          POTPAR1(ISP)%QBAR(ITAIL)    = QBAR
!         == JBAR -> |PHIBARDOT> JBARTOPHIDOT ==================================
          POTPAR1(ISP)%JBARTOPHIDOT(ITAIL)=(JVAL-KVAL*QBAR)/PHIDOTVAL
!
!         ==  <PRO|JBAR_AUG> ===================================================
          AUX=R**2*PSPHIDOT(:,LNOFT)*POTPAR1(ISP)%JBARTOPHIDOT(ITAIL)
          DO LN=1,LNX1
            IF(LOX(LN,ISP).NE.L) CYCLE
            CALL RADIAL$INTEGRAL(GID,NR,PRO(:,LN)*AUX,SVAR)
            POTPAR1(ISP)%PROJBAR(LN,ITAIL)=SVAR
          ENDDO
!
          DO IHEAD=1,POTPAR1(ISP)%NHEAD
            IF(POTPAR1(ISP)%LOFH(IHEAD).NE.L) CYCLE
            LNOFH=POTPAR1(ISP)%LNOFH(IHEAD)
!           ====================================================================
!           == VALUE AND DERIVATIVE OF PARTIAL WAVES                          ==
!           == PHIDOT IS PHIBARDOT                                            ==
!           ====================================================================
            CALL RADIAL$VALUE(GID,NR,NLPHI(:,LNOFH),RAD,PHIVAL)
            CALL RADIAL$DERIVATIVE(GID,NR,NLPHI(:,LNOFH),RAD,PHIDER)
            WJPHI=JVAL*PHIDER-JDER*PHIVAL
            WKPHI=KVAL*PHIDER-KDER*PHIVAL
            WPHIPHIDOT=PHIVAL*PHIDOTDER-PHIDER*PHIDOTVAL
            WJBARPHI=WJPHI-WKPHI*QBAR
!           == K    -> |PHI>KTOPHI+|PHIBARDOT> KTOPHIDOT =======================
            POTPAR1(ISP)%KTOPHIDOT(IHEAD)=-WKPHI/WPHIPHIDOT
            POTPAR1(ISP)%KTOPHI(IHEAD)   = WKPHIDOT/WPHIPHIDOT
!
!           ==  <PRO|PSPHIDOT> =================================================
            AUX=R**2*PRO(:,LNOFH)*PSPHIDOT(:,LNOFT)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            POTPAR1(ISP)%PHIDOTPROJ(IHEAD)=SVAR
!
!           ==  <PRO|K_AUG> ===================================================
            AUX=R**2*(PSPHI(:,LNOFH)   *POTPAR1(ISP)%KTOPHI(IHEAD) &
      &              +PSPHIDOT(:,LNOFT)*POTPAR1(ISP)%KTOPHIDOT(IHEAD) )
            DO LN=1,LNX1
              IF(LOX(LN,ISP).NE.L) CYCLE
              CALL RADIAL$INTEGRAL(GID,NR,PRO(:,LN)*AUX,SVAR)
              POTPAR1(ISP)%PROK(LN,IHEAD)=SVAR
            ENDDO
!
!!$PRINT*,'LN=',LN,'================================='
!!$PRINT*,'PHIVAL ',PHIVAL
!!$PRINT*,'PHIDER ',PHIDER
!!$PRINT*,'PHIDOTVAL ',PHIDOTVAL
!!$PRINT*,'PHIDOTDER ',PHIDOTDER
!!$PRINT*,'JVAL ',JVAL
!!$PRINT*,'JDER ',JDER
!!$PRINT*,'KVAL ',KVAL
!!$PRINT*,'KDER ',KDER
!!$PRINT*,'WJPHI     ',WJPHI     
!!$PRINT*,'WJPHIDOT  ',WJPHIDOT  
!!$PRINT*,'WKPHI     ',WKPHI     
!!$PRINT*,'WKPHIDOT  ',WKPHIDOT  
!!$PRINT*,'WPHIPHIDOT',WPHIPHIDOT
!!$PRINT*,'QBAR      ',QBAR      
!!$PRINT*,'WJBARPHI  ',WJBARPHI  
!!$PRINT*,'KTOPHI       ',POTPAR1(ISP)%KTOPHI(IHEAD)    
!!$PRINT*,'KTOPHIDOT    ',POTPAR1(ISP)%KTOPHIDOT(IHEAD)
!!$PRINT*,'JBARTOPHIDOT ',POTPAR1(ISP)%JBARTOPHIDOT(IHEAD)
!
          ENDDO !END OF LOOP OVER HEAD FUNCTIONS
        ENDDO !END OF LOOP OVER ANGULAR MOMENTA
!
        DEALLOCATE(AUX)
        DEALLOCATE(ISCATT)
        DEALLOCATE(TORB)
        DEALLOCATE(NLPHI)
        DEALLOCATE(NLPHIDOT)
        DEALLOCATE(PSPHI)
        DEALLOCATE(PSPHIDOT)
        DEALLOCATE(PRO)
        DEALLOCATE(R)
        CALL SETUP$UNSELECT()
      ENDDO
!
!     ==========================================================================
!     == OVERLAP MATRIX OF PARTIAL WAVES WITHIN SPHERE                        ==
!     == WILL BE USED AS WEIGHT TO ADJUST THE LOCAL ORBITAL EXPANSION         ==
!     ==========================================================================
      CALL LMTO_OVERLAPPHI()
!
!     ==========================================================================
!     == PRINTOUT                                                             ==
!     ==========================================================================
      IF(TPRINT) THEN
        WRITE(*,FMT='(80("="),T10," LMTO_MAKEPOTPAR1  ")')
        DO ISP=1,NSP
          WRITE(*,FMT='(80("-"),T10," ATOM TYPE ",I2," ")')ISP
          NTAIL=POTPAR1(ISP)%NTAIL
          DO ITAIL=1,NTAIL
            L=POTPAR1(ISP)%LOFT(ITAIL)
            WRITE(*,FMT='("L=",I2," QBAR=",F10.4," JBARTOPHIDOT=",F10.4)') &
     &                   L,POTPAR1(ISP)%QBAR(ITAIL) &
     &                    ,POTPAR1(ISP)%JBARTOPHIDOT(ITAIL)
          ENDDO
!
          NHEAD=POTPAR1(ISP)%NHEAD
          DO IHEAD=1,NHEAD
            L=POTPAR1(ISP)%LOFH(IHEAD)
            WRITE(*,FMT='("L=",I2," K->PHI=",F10.4," K->PHIDOT=",F10.4' &
     &                  //'," <PHIDOT|P>=",F10.4)') &
     &              L,POTPAR1(ISP)%KTOPHI(IHEAD),POTPAR1(ISP)%KTOPHIDOT(IHEAD) &
     &               ,POTPAR1(ISP)%PHIDOTPROJ(IHEAD)
          ENDDO
        ENDDO ! END OF LOOP OVER ATOM TYPES
        WRITE(*,FMT='(82("="),T10," FINISHED  ")')
      END IF
                             CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OVERLAPPHI()
!     **************************************************************************
!     **  OVERLAP BETWEEN ALL-ELECTRON PARTIAL WAVES WITHIN SPHERE            **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR1,NSP,LNX,LOX
      IMPLICIT NONE
      INTEGER(4)             :: GID        ! GRID ID
      INTEGER(4)             :: NR         ! #(RADIAL GRID POINTS)
      INTEGER(4)             :: LNX1       ! #(PARTIAL WAVES)
      REAL(8)   ,ALLOCATABLE :: R(:)       ! RADIAL GRID
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:) ! ALL-ELECTRON PARTIAL WAVES
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)                :: VAL
      INTEGER(4)             :: ISP,LN1,LN2 
      LOGICAL(4),PARAMETER   :: TOLD=.FALSE.
      LOGICAL(4),ALLOCATABLE :: TORB(:)
!     **************************************************************************
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
        ALLOCATE(AUX(NR))
!
!       == PARTIAL WAVES AND PROJECTORS ========================================
        ALLOCATE(AEPHI(NR,LNX1))
        CALL SETUP$GETR8A('AEPHI',NR*LNX1,AEPHI)
        ALLOCATE(POTPAR1(ISP)%PHIOV(LNX1,LNX1))
        DO LN1=1,LNX1
          DO LN2=LN1,LNX1
            IF(LOX(LN1,ISP).NE.LOX(LN2,ISP)) CYCLE
            CALL RADIAL$INTEGRATE(GID,NR,R**2*AEPHI(:,LN1)*AEPHI(:,LN2),AUX)
            CALL RADIAL$VALUE(GID,NR,AUX,POTPAR1(ISP)%RAD,VAL)
            POTPAR1(ISP)%PHIOV(LN1,LN2)=VAL
            POTPAR1(ISP)%PHIOV(LN2,LN1)=VAL
          ENDDO
        ENDDO
        DEALLOCATE(R)
        DEALLOCATE(AUX)
        DEALLOCATE(AEPHI)
        CALL SETUP$UNSELECT()
      ENDDO
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
      SUBROUTINE LMTO_CONSOLIDATETAILEDPARMS()
!     **************************************************************************
!     ** CONSOLIDATE THE PARAMETERS OF HYBRIDSETTING                          **
!     ** I.E. INITIALIZE UNSET VALUES TO ACCOUNT FOR ELEMENT SPECIFIC DEFAULTS**
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2015******************
      USE LMTO_MODULE, ONLY : NSP &
     &                       ,HYBRIDSETTING 
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)   :: ISP
      REAL(8)      :: AEZ
      REAL(8)      :: RCOV
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_CONSOLIDATETAILEDPARMS')
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
        CALL SETUP$UNSELECT()
!
!       ========================================================================
!       ==  SET DEFAULTS FOR PARAMETERS FOR THE TAILED REPRESENTATION         ==
!       ========================================================================

!       == MATCHING RADIUS =====================================================
        IF(HYBRIDSETTING(ISP)%RAUG.LE.0.D0) THEN
!         __THE RADIUS MUST CIRCUMSCRIBE SEMI-CORE STATES, BECAUSE THE DECAY____
!         __OF THE HANKEL FUNCTIONS IS USUALLY ADAPTED TO THE VALENCE STATES____
!         __AND THEIR DECAY IS USUALLY MUCH TO SLOW FOR SEMI-CORE STATES________
          HYBRIDSETTING(ISP)%RAUG=1.1D0*RCOV
        END IF
!
!       == RADIUS FOR MATCHING EXPONENTIAL TAILS TO HANKEL FKT.=================
        IF(HYBRIDSETTING(ISP)%RTAIL.LE.0.D0) THEN
          HYBRIDSETTING(ISP)%RTAIL=1.2*RCOV
        END IF
!
!       == RAPIDLY DECAYING EXPONENTIAL ========================================
        IF(HYBRIDSETTING(ISP)%TAILEDLAMBDA1.LE.0.D0) THEN
          HYBRIDSETTING(ISP)%TAILEDLAMBDA1=4.D0
        END IF
!
!       == SLOWLY DECAYING EXPONENTIAL =========================================
        IF(HYBRIDSETTING(ISP)%TAILEDLAMBDA2.LE.0.D0) THEN
          HYBRIDSETTING(ISP)%TAILEDLAMBDA2=2.D0
        END IF
      ENDDO
                              CALL TRACE$POP()
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_MAKETAILEDPARTIALWAVES_WITHPOTPAR1()
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
      USE LMTO_MODULE, ONLY : K2 &
     &                       ,POTPAR1 &
     &                       ,NSP &
     &                       ,LNX &
     &                       ,LOX &
     &                       ,HYBRIDSETTING &
     &                       ,TCTE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      INTEGER(4),PARAMETER   :: LXHIGHER=4
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
      INTEGER(4)             :: NHIGHER
      INTEGER(4)             :: LRX,LMRX
      INTEGER(4)             :: LNXT
      INTEGER(4)             :: LMNXT
      INTEGER(4)             :: NHEAD ! #(HEAD FUNCTIONS)
      INTEGER(4)             :: NTAIL ! #(TAIL FUNCTIONS)
      INTEGER(4)             :: ISP,LN,LN1,LN2,LNT,LMN,LMN1,LMN2,IM,IR,LNDOT
      INTEGER(4)             :: IH,IT
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)                :: RTAIL     ! TAIL-MATCHING RADIUS
      REAL(8)                :: RTAILCUT  ! TAIL-TRUNCATION RADIUS 
      LOGICAL(4)             :: TCHKHIGHER(LXHIGHER+1)
      CHARACTER(64)           :: STRING
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_MAKETAILEDPARTIALWAVES')
      DO ISP=1,NSP
        LAMBDA1=HYBRIDSETTING(ISP)%TAILEDLAMBDA1
        LAMBDA2=HYBRIDSETTING(ISP)%TAILEDLAMBDA2
        RTAIL=HYBRIDSETTING(ISP)%RTAIL
        RTAILCUT=HYBRIDSETTING(ISP)%RTAILCUT
!
        CALL SETUP$ISELECT(ISP)
!       == RADIAL GRID =========================================================
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        RAD=POTPAR1(ISP)%RAD    !AUGMENTATION RADIUS
        DO IR=1,NR
          IRAD=IR
          IF(R(IR).GT.RAD) EXIT
        ENDDO
        POTPAR1(ISP)%TAILED%GID=GID
        NHEAD=POTPAR1(ISP)%NHEAD
        NTAIL=POTPAR1(ISP)%NTAIL
!       __DETERMINE NUMBER OF HIGHER COMPONENTS PURE BESSEL
        NHIGHER=LXHIGHER+1
        DO L=0,LXHIGHER
          TCHKHIGHER(L+1)=.TRUE.
          DO IT=1,NTAIL
            IF(POTPAR1(ISP)%LOFT(IT).EQ.L) THEN
              TCHKHIGHER(L+1)=.FALSE.
              NHIGHER=NHIGHER-1
              EXIT
            END IF
          ENDDO
        ENDDO
        IF(.NOT.TCTE) THEN
          TCHKHIGHER=.FALSE.
          NHIGHER=0
        END IF
!
!       == COUNT PARTIAL WAVES HEAD+TAIL =======================================
        LNXT=NHEAD+NTAIL+NHIGHER
        POTPAR1(ISP)%TAILED%LNX=LNXT
!
!       == FIX ANGULAR MOMENTA FOR HEAD+TAIL REPRESENTATION ====================
        ALLOCATE(POTPAR1(ISP)%TAILED%LOX(LNXT))
        POTPAR1(ISP)%TAILED%LOX(1:NHEAD) =POTPAR1(ISP)%LOFH
        POTPAR1(ISP)%TAILED%LOX(NHEAD+1:NHEAD+NTAIL)=POTPAR1(ISP)%LOFT
        IH=NHEAD+NTAIL
        DO L=0,LXHIGHER
          IF(.NOT.TCHKHIGHER(L+1)) CYCLE
          IH=IH+1
          POTPAR1(ISP)%TAILED%LOX(IH)=L
        ENDDO
        LMNXT=SUM(2*POTPAR1(ISP)%TAILED%LOX+1)
        POTPAR1(ISP)%TAILED%LMNX=LMNXT
!
!       ========================================================================
!       == AUGMENTED HANKEL AND BESSEL FUNCTIONS WITH TAILS ATTACHED          ==
!       ========================================================================
        ALLOCATE(POTPAR1(ISP)%TAILED%AEF(NR,LNXT))
        ALLOCATE(POTPAR1(ISP)%TAILED%PSF(NR,LNXT))
        ALLOCATE(POTPAR1(ISP)%TAILED%NLF(NR,LNXT))
!
!       ========================================================================
!       == INITIALIZE NLF FOR R<RTAIL AS THE ENVELOPE FUNCTION                ==
!       == THE HEAD FUNCTIONS ARE PURE HANKEL FUNCTIONS                       ==
!       == THE TAIL FUNCTIONS ARE SCREENED BESSEL FUNCTIONS                   ==
!       ========================================================================
        POTPAR1(ISP)%TAILED%NLF(:,:)=0.D0
        DO IH=1,NHEAD
          LNT=IH
          L=POTPAR1(ISP)%LOFH(IH)
          DO IR=IRAD,NR   ! WILL BE AUGMENTED FROM 1 TO IRAD-1
            IF(R(IR).GT.RTAIL) EXIT
            CALL LMTO$SOLIDHANKELRAD(L,R(IR),K2,KVAL,KDER)
            POTPAR1(ISP)%TAILED%NLF(IR,LNT)=KVAL
          ENDDO
        ENDDO
        DO IT=1,NTAIL
          LNT=NHEAD+IT
          L=POTPAR1(ISP)%LOFT(IT)
          QBAR=POTPAR1(ISP)%QBAR(IT)
          DO IR=IRAD,NR ! WILL BE AUGMENTED FROM 1 TO IRAD-1
            IF(R(IR).GT.RTAIL) EXIT
            CALL LMTO$SOLIDHANKELRAD(L,R(IR),K2,KVAL,KDER)
            CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,JVAL,JDER)
            POTPAR1(ISP)%TAILED%NLF(IR,LNT)=JVAL-KVAL*QBAR
          ENDDO
        ENDDO
        DO IH=1,NHIGHER
          LNT=NHEAD+NTAIL+IH
          L=POTPAR1(ISP)%TAILED%LOX(LNT)
          DO IR=1,NR ! WILL NOT BE AUGMENTED FROM 1 TO IRAD-1
            IF(R(IR).GT.RTAIL) EXIT
            CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,JVAL,JDER)
            POTPAR1(ISP)%TAILED%NLF(IR,LNT)=JVAL
          ENDDO
        ENDDO
!
!       ========================================================================
!       == TAIL PART: CONSTRUCT EXPONENTIAL TAILS OF THE TAILED ORBITALS      ==
!       == THE HEAD FUNCTIONS ARE PURE HANKEL FUNCTIONS                       ==
!       == THE TAIL FUNCTIONS ARE SCREENED BESSEL FUNCTIONS                   ==
!       ========================================================================
!       == HEAD FUNCTIONS ======================================================
        DO IH=1,NHEAD
          LNT=IH
          L=POTPAR1(ISP)%LOFH(IH)
          CALL LMTO$SOLIDHANKELRAD(L,RTAIL,K2,KVAL,KDER)
!         -- DETERMINE VALUE AND LOGARITHMIC DERIVATIVE OF PHI AND PHIBARDOT----
          A1=KVAL*(KDER/KVAL+LAMBDA2)/(LAMBDA2-LAMBDA1)
          A2=KVAL*(KDER/KVAL+LAMBDA1)/(LAMBDA1-LAMBDA2)
          DO IR=1,NR
            IF(R(IR).LE.RTAIL) CYCLE
            SVAR1=EXP(-LAMBDA1*(R(IR)-RTAIL))
            SVAR2=EXP(-LAMBDA2*(R(IR)-RTAIL))
            POTPAR1(ISP)%TAILED%NLF(IR,LNT)=A1*SVAR1+A2*SVAR2
          ENDDO
        ENDDO
!
!       == TAIL FUNCTIONS ======================================================
        DO IT=1,NTAIL
          LNT=NHEAD+IT
          L=POTPAR1(ISP)%LOFT(IT)
          CALL LMTO$SOLIDHANKELRAD(L,RTAIL,K2,KVAL,KDER)
          CALL LMTO$SOLIDBESSELRAD(L,RTAIL,K2,JVAL,JDER)
!         -- TRANSFORM UNSCREENED BESSEL FUNCTION TO SCREENED BESSEL FUNCTION
          QBAR=POTPAR1(ISP)%QBAR(IT)
          JVAL=JVAL-KVAL*QBAR
          JDER=JDER-KDER*QBAR
!         -- DETERMINE VALUE AND LOGARITHMIC DERIVATIVE OF PHI AND PHIBARDOT----
          B1=JVAL*(JDER/JVAL+LAMBDA2)/(LAMBDA2-LAMBDA1)
          B2=JVAL*(JDER/JVAL+LAMBDA1)/(LAMBDA1-LAMBDA2)
          DO IR=1,NR
            IF(R(IR).LE.RTAIL) CYCLE
            SVAR1=EXP(-LAMBDA1*(R(IR)-RTAIL))
            SVAR2=EXP(-LAMBDA2*(R(IR)-RTAIL))
            POTPAR1(ISP)%TAILED%NLF(IR,LNT)=B1*SVAR1+B2*SVAR2
          ENDDO
        ENDDO  
!
!       == HIGHER PARTIAL WAVES (PURE BESSEL FUNCTIONS) ========================
        DO IH=1,NHIGHER
          LNT=NHEAD+NTAIL+IH
          L=POTPAR1(ISP)%TAILED%LOX(LNT)
          CALL LMTO$SOLIDBESSELRAD(L,RTAIL,K2,JVAL,JDER)
!         -- DETERMINE VALUE AND LOGARITHMIC DERIVATIVE OF PHI AND PHIBARDOT----
          IF(JVAL.EQ.0.D0) THEN
            CALL ERROR$MSG('DIVIDE BY ZERO AHEAD. FIX CODE!')
            CALL ERROR$I4VAL('IH',IH)
            CALL ERROR$I4VAL('NHIGHER',NHIGHER)
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$R8VAL('RTAIL',RTAIL)
            CALL ERROR$R8VAL('K2',K2)
            CALL ERROR$R8VAL('JVAL',JVAL)
            CALL ERROR$STOP('LMTO_MAKETAILEDPARTIALWAVES_WITHPOTPAR1')
          END IF
          B1=JVAL*(JDER/JVAL+LAMBDA2)/(LAMBDA2-LAMBDA1)
          B2=JVAL*(JDER/JVAL+LAMBDA1)/(LAMBDA1-LAMBDA2)
          DO IR=1,NR
            IF(R(IR).LE.RTAIL) CYCLE
            SVAR1=EXP(-LAMBDA1*(R(IR)-RTAIL))
            SVAR2=EXP(-LAMBDA2*(R(IR)-RTAIL))
            POTPAR1(ISP)%TAILED%NLF(IR,LNT)=B1*SVAR1+B2*SVAR2
          ENDDO
        ENDDO  
!
!       ========================================================================
!       == AUGMENTATION ========================================================
!       ========================================================================
!       == NODELESS SPHERE PART. RESULTS IN CONTINUOUS ORBITAL =================
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
        DO IH=1,NHEAD
          IT=POTPAR1(ISP)%ITAIL(IH)
          LN=POTPAR1(ISP)%LNOFH(IH)
          L=POTPAR1(ISP)%LOFH(IH)
          LNDOT=POTPAR1(ISP)%LNOFT(IT)  ! PARTIAL WAVE INDEX FOR PHIDOT
          A1=POTPAR1(ISP)%KTOPHI(IH)
          A2=POTPAR1(ISP)%KTOPHIDOT(IH)
!         == NODESLESS WAVE FUNCTIONS INSERTED ONLY UP TO MATCHING RADIUS
          POTPAR1(ISP)%TAILED%NLF(:IRAD-1,IH)=   NLPHI(:IRAD-1,LN)   *A1 &
       &                                     +NLPHIDOT(:IRAD-1,LNDOT)*A2
!         ==  ADD DIFFERENCE TO FULL WAVE FUNCTIONS ============================
          POTPAR1(ISP)%TAILED%AEF(:,IH)=POTPAR1(ISP)%TAILED%NLF(:,IH) &
       &                             +(   AEPHI(:,LN)   -   NLPHI(:,LN)   )*A1 &
       &                             +(AEPHIDOT(:,LNDOT)-NLPHIDOT(:,LNDOT))*A2
          POTPAR1(ISP)%TAILED%PSF(:,IH)=POTPAR1(ISP)%TAILED%NLF(:,IH) &
       &                             +(   PSPHI(:,LN)   -   NLPHI(:,LN)   )*A1 &
       &                             +(PSPHIDOT(:,LNDOT)-NLPHIDOT(:,LNDOT))*A2
        ENDDO
        DO IT=1,NTAIL
          LNT=NHEAD+IT
          LNDOT=POTPAR1(ISP)%LNOFT(IT)
          A2=POTPAR1(ISP)%JBARTOPHIDOT(IT)
!         == NODESLESS WAVE FUNCTIONS INSERTED ONLY UP TO MATCHING RADIUS
          POTPAR1(ISP)%TAILED%NLF(:IRAD-1,LNT)=NLPHIDOT(:IRAD-1,LNDOT)*A2
!         ==  ADD DIFFERENCE TO FULL WAVE FUNCTIONS ============================
          POTPAR1(ISP)%TAILED%AEF(:,LNT)=POTPAR1(ISP)%TAILED%NLF(:,LNT) &
      &                              +(AEPHIDOT(:,LNDOT)-NLPHIDOT(:,LNDOT))*A2
          POTPAR1(ISP)%TAILED%PSF(:,LNT)=POTPAR1(ISP)%TAILED%NLF(:,LNT) &
      &                              +(PSPHIDOT(:,LNDOT)-NLPHIDOT(:,LNDOT))*A2
        ENDDO
!       == NO AUGMENTATION FOR HIGHER PARTIAL WAVES (VIDE SUPRA) ===============
        DO IH=1,NHIGHER
          LNT=NHEAD+NTAIL+IH
          POTPAR1(ISP)%TAILED%AEF(:,LNT)=POTPAR1(ISP)%TAILED%NLF(:,LNT)
          POTPAR1(ISP)%TAILED%PSF(:,LNT)=POTPAR1(ISP)%TAILED%NLF(:,LNT)
        ENDDO
!
!       ========================================================================
!       ==  CUT OFF TAILS FOR STABILITY                                       ==
!       ========================================================================
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
        DO IR=1,NR
          IF(R(IR).GT.RTAILCUT) THEN  
            POTPAR1(ISP)%TAILED%NLF(IR:,:)=0.D0
            POTPAR1(ISP)%TAILED%AEF(IR:,:)=0.D0
            POTPAR1(ISP)%TAILED%PSF(IR:,:)=0.D0
            EXIT
          END IF
        ENDDO

        DEALLOCATE(AEPHI)
        DEALLOCATE(AEPHIDOT)
        DEALLOCATE(PSPHI)
        DEALLOCATE(PSPHIDOT)
        DEALLOCATE(NLPHI)
        DEALLOCATE(NLPHIDOT)
        DEALLOCATE(R)
!
        CALL SETUP$UNSELECT()
!
!       ========================================================================
!       === WRITE TAILED FUNCTIONS TO FILE                                    ==
!       ========================================================================
        IF(TPR) THEN
          WRITE(STRING,*)ISP
          STRING=ADJUSTL(STRING) 
          CALL LMTO_WRITEPHI(TRIM(STRING)//'_AEF.DAT',GID,NR &
       &                    ,LNXT,POTPAR1(ISP)%TAILED%AEF)
          CALL LMTO_WRITEPHI(TRIM(STRING)//'_PSF.DAT',GID,NR &
       &                    ,LNXT,POTPAR1(ISP)%TAILED%PSF)
          CALL LMTO_WRITEPHI(TRIM(STRING)//'_NLF.DAT',GID,NR &
       &                    ,LNXT,POTPAR1(ISP)%TAILED%NLF)
        END IF
      ENDDO ! END OF LOOP OVER ATOM TYPES
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_MAKETAILEDMATRIXELEMENTS_WITHPOTPAR1()
!     **************************************************************************
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : POTPAR1,NSP
!      USE LMTO_MODULE, ONLY : K2,POTPAR1,NSP,LNX,LOX,HYBRIDSETTING
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      INTEGER(4)             :: GID   ! GRID ID
      INTEGER(4)             :: NR    ! #(RADIAL GRID POINTS)
      INTEGER(4)             :: LMRX  ! #(ANGULAR MOMENTUM FOR DENSITY)
      INTEGER(4)             :: LRX   ! X(ANGULAR MOMENTUM FOR DENSITY)
      INTEGER(4)             :: LMNXT
      INTEGER(4)             :: LNXT
      INTEGER(4),ALLOCATABLE :: LOXT(:)
      REAL(8)   ,ALLOCATABLE :: ULITTLE(:,:,:,:,:)
      INTEGER(4)             :: ISP,LMNT,LNT
!     **************************************************************************
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        CALL SETUP$GETI4('LMRX',LMRX)
        LMNXT=POTPAR1(ISP)%TAILED%LMNX
        LNXT=POTPAR1(ISP)%TAILED%LNX
        ALLOCATE(LOXT(LNXT))
        LOXT(:)=POTPAR1(ISP)%TAILED%LOX(:)
!
!       ========================================================================
!       == ONSITE U-TENSOR OF TAILED PARTIAL WAVES                            ==
!       ========================================================================
        LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
        ALLOCATE(ULITTLE(LRX+1,LNXT,LNXT,LNXT,LNXT))
        CALL LMTO_ULITTLE(GID,NR,LRX,LNXT,LOXT,POTPAR1(ISP)%TAILED%AEF,ULITTLE)

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
        ALLOCATE(POTPAR1(ISP)%TAILED%U(LMNXT,LMNXT,LMNXT,LMNXT))
        CALL LMTO_UTENSOR(LRX,LMNXT,LNXT,LOXT,ULITTLE,POTPAR1(ISP)%TAILED%U)
        DEALLOCATE(ULITTLE)
!
!       ========================================================================
!       == ONSITE OVERLAP MATRIX                                              ==
!       ========================================================================
        ALLOCATE(POTPAR1(ISP)%TAILED%OVERLAP(LMNXT,LMNXT))
        CALL LMTO_ONECENTEROVERLAP(GID,NR,LNXT,LOXT,POTPAR1(ISP)%TAILED%AEF &
     &                            ,LMNXT,POTPAR1(ISP)%TAILED%OVERLAP)
!
!       ========================================================================
!       == MONO- AND DIPOLE MATRIX ELEMENTS                                   ==
!       == USED FOR LONG-DISTANCE MATRIX ELEMENTS                             ==
!       ========================================================================
        ALLOCATE(POTPAR1(ISP)%TAILED%QLN(2,LNXT,LNXT))
        CALL LMTO_ONECENTERQLN(GID,NR,LNXT,LOXT,POTPAR1(ISP)%TAILED%AEF &
     &                            ,POTPAR1(ISP)%TAILED%QLN)
        DEALLOCATE(LOXT)
        CALL SETUP$UNSELECT()
      ENDDO ! END OF LOOP OVER ATOM TYPES
!
!     ==========================================================================
!     == PRINTOUT                                                             ==
!     ==========================================================================
      IF(TPRINT) THEN
        WRITE(*,FMT='(82("="),T10," LMTO_MAKETAILEDMATRIXELEMENTS  ")')
        DO ISP=1,NSP
          WRITE(*,FMT='(82("-"),T10," ATOM TYPE=",I3,"  ")')ISP
          DO LMNT=1,POTPAR1(ISP)%TAILED%LMNX
            WRITE(*,FMT='("LMN=",I3," O=",100F10.5)') &
     &                     LMNT,POTPAR1(ISP)%TAILED%OVERLAP(:,LMNT)
          ENDDO
          DO LNT=1,POTPAR1(ISP)%TAILED%LNX
            WRITE(*,FMT='("LN=",I3," MONOPOLE=",100F10.5)') &
     &                     LNT,POTPAR1(ISP)%TAILED%QLN(1,:,LNT)
          ENDDO
          DO LNT=1,POTPAR1(ISP)%TAILED%LNX
            WRITE(*,FMT='("LN=",I3," DIPOLE=",100F10.5)') &
     &                     LNT,POTPAR1(ISP)%TAILED%QLN(2,:,LNT)
          ENDDO
        ENDDO
!!$        CALL ERROR$MSG('STOPPING AFTER PRINTOUT')
!!$        CALL ERROR$STOP('LMTO_MAKETAILEDMATRIXELEMENTS_WITHPOTPAR1')
      END IF
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
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)            :: LN1,LN2,LN3,LN4
      INTEGER(4)            :: L
      INTEGER(4)            :: LMIN,LMAX,ISVAR1,ISVAR2
      REAL(8)               :: RHO(NR)
      REAL(8)               :: POT(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR
      REAL(8)               :: R(NR)
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
      DO L=0,LRX
        ULITTLE(L+1,:,:,:,:)=ULITTLE(L+1,:,:,:,:)*REAL(2*L+1,KIND=8)/(4.D0*PI)
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
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: FOURPI=4.D0*PI
      INTEGER(4)            :: LN1,LN2,LN3,LN4
      INTEGER(4)            :: IORB1,IORB2,IORB3,IORB4
      INTEGER(4)            :: L1,L2,L3,L4
      INTEGER(4)            :: LM1,LM2,LM3,LM4
      INTEGER(4)            :: M1,M2,M3,M4
      INTEGER(4)            :: L,M,LM,LX
      REAL(8)               :: CG1,CG2
      REAL(8)               :: SVAR
      REAL(8)               :: FOURPIBY2LPLUS1
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_UTENSOR')
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
!     **************************************************************************
!     **  ONSITE OVERLAP MATRIX IN THE TAILED REPRESENTATION                  **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LMNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: CHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: OVERLAP(LMNX,LMNX)
      LOGICAL(4),PARAMETER  :: TPRINT=.FALSE.
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: LMN10,LMN20,LMN
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
!
!     ==========================================================================
!     == DIAGNOSTIC PRINTOUT                                                  ==
!     ==========================================================================
      IF(TPRINT) THEN
        WRITE(*,FMT='(80("="),T20," TAILED OVERLAP MATRIX ")')
        DO LMN=1,LMNX
          WRITE(*,FMT='(10F10.5)')OVERLAP(LMN,:)
        ENDDO
      END IF      
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
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: SQ4PI=SQRT(4.D0*PI)
      REAL(8)   ,PARAMETER  :: SQ4PITHIRD=SQRT(4.D0*PI/3.D0)
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: L1,L2,IM
      REAL(8)               :: AUX(NR),SVAR
      REAL(8)               :: R(NR)
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_ONECENTERMULTIPOLE')
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
      USE LMTO_MODULE, ONLY : NSP &
     &                       ,POTPAR=>POTPAR1
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
      USE LMTO_MODULE, ONLY : POTPAR=>POTPAR1,NSP
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
        IF(POTPAR(ISP)%TAILED%LNX.GT.0) THEN
          L=MAXVAL(POTPAR(ISP)%TAILED%LOX(:))
        ELSE
          L=-1
        END IF
        NPOW=POTPAR(ISP)%TAILED%GAUSSNLF%NIJK  !#(DOUBLE POWERS)
        NPOWX=MAX(NPOW,NPOWX)
        LX=MAX(LX,L)
        NX=MAX(NX,2*(NPOW-1))   !HIGHEST POWER
      ENDDO
      IF(LX.EQ.-1) THEN 
        CALL ERROR$MSG('LX MUST NOT BE -1')
        CALL ERROR$MSG('PROBABLY POTPAR%TAILED%LNX=0, WHICH IS NOT SUPPORTED')
        CALL ERROR$MSG('POSSIBLY NO LOCAL ORBITALS ARE SPECIFIED')
        CALL ERROR$MSG('IN !STRUCTURE!NTBO:NOFL')
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('LMTO_TAILEDGAUSSORBTOYLM')
      END IF
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
      USE LMTO_MODULE, ONLY : NSP &
     &                       ,POTPAR=>POTPAR1 &
     &                       ,OFFSITEX &
     &                       ,HYBRIDSETTING
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
        IF(.NOT.HYBRIDSETTING(ISPA)%TBONDX) CYCLE
        NIJKA=POTPAR(ISPA)%TAILED%GAUSSNLF%NIJK
        LMNXA=POTPAR(ISPA)%TAILED%GAUSSNLF%NORB
        NEA  =POTPAR(ISPA)%TAILED%GAUSSNLF%NE
        ALLOCATE(EA(NEA))
        EA(:)=POTPAR(ISPA)%TAILED%GAUSSNLF%E(:)
        ALLOCATE(ORBA(NIJKA,NEA,LMNXA))
        ORBA(:,:,:)=POTPAR(ISPA)%TAILED%GAUSSNLF%C(:,:,:)
        DO ISPB=1,NSP
          IF(.NOT.HYBRIDSETTING(ISPB)%TBONDX) CYCLE
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
          IF(.NOT.(HYBRIDSETTING(ISPA)%TBONDX.AND.HYBRIDSETTING(ISPB)%TBONDX)) &
    &       CYCLE
          CALL MPE$COMBINE('MONOMER','+',OFFSITEX(ISPA,ISPB)%BONDU)
        ENDDO
      ENDDO      
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PCHI()
!     **************************************************************************
!     ** CALCULATE PROJECTIONS <PTILDE|CHITILDE> ONTO LOCAL ORBITALS          **
!     ** WILL BE USED TO EXTRACT THE LOCAL-ORBITAL PROJECTIONS OF THE         **
!     ** OF THE WAVE FUNCTIONS
!     **************************************************************************
      USE LMTO_MODULE, ONLY : NSP &  ! #(ATOM TYPES)
     &                       ,LNX &  !(NSP) 
     &                       ,LOX &  !(LNXX,NSP)
     &                       ,ISPECIES & !(NAT) ATOM TYPE OF ATOM IAT
     &                       ,POTPAR1  & !(NSP)
     &                       ,SBAR_NEW & !(NNS) STRUCTURE CONSTANTS
     &                       ,PCHI       !(NNS) <PTILDE|CHITILDE>
      IMPLICIT NONE
      TYPE PKJ_TYPE
        REAL(8),ALLOCATABLE :: PROK(:,:)
        REAL(8),ALLOCATABLE :: PROJBAR(:,:)
      END TYPE PKJ_TYPE
      TYPE(PKJ_TYPE)          :: PKJ(NSP)
      INTEGER(4),ALLOCATABLE  :: LMNP0(:)
      INTEGER(4)              :: LMNPX
      INTEGER(4)              :: LMNHX
      INTEGER(4)              :: LMNTX
      INTEGER(4)              :: ISP,L,M
      INTEGER(4)              :: LMNP,LNP,LMNH,IH,LH,LMNT,IT,LT
      INTEGER(4)              :: NN,NNS
      INTEGER(4)              :: IAT1,IAT2,ISP1,ISP2
!     **************************************************************************
!
!     ==========================================================================
!     == <PTILDE|CHI>                                                         ==
!     ==========================================================================
      DO ISP=1,NSP
        ALLOCATE(LMNP0(LNX(ISP)))
        LMNP=0
        DO LNP=1,LNX(ISP)
          LMNP0(LNP)=LMNP
          LMNP=LMNP+2*LOX(LNP,ISP)+1
        ENDDO
        LMNPX=SUM(2*LOX(:LNX(ISP),ISP)+1)  ! #(PROJECTOR FUNCTIONS INCLUDING M)
        LMNHX=SUM(2*POTPAR1(ISP)%LOFH+1)   ! #(HEAD FUNCTIONS INCLUDING M)
        LMNTX=SUM(2*POTPAR1(ISP)%LOFT+1)   ! #(TAIL FUNCTIONS INCLUDING M)

        ALLOCATE(PKJ(ISP)%PROK(LMNPX,LMNHX))
        PKJ(ISP)%PROK=0.D0
        DO LNP=1,LNX(ISP)
          L=LOX(LNP,ISP)
          LMNH=0
          DO IH=1,POTPAR1(ISP)%NHEAD
            LH=POTPAR1(ISP)%LOFH(IH)
            IF(LH.EQ.L) THEN
              LMNP=LMNP0(LNP)
              DO M=1,2*LH+1
                LMNP=LMNP+1
                LMNH=LMNH+1
                PKJ(ISP)%PROK(LMNP,LMNH)=POTPAR1(ISP)%PROK(LNP,IH)
              ENDDO
            ELSE
              LMNH=LMNH+2*LH+1
            END IF
          ENDDO
        ENDDO

        ALLOCATE(PKJ(ISP)%PROJBAR(LMNPX,LMNTX))
        PKJ(ISP)%PROJBAR=0.D0
        DO LNP=1,LNX(ISP)
          L=LOX(LNP,ISP)
          LMNT=0
          DO IT=1,POTPAR1(ISP)%NTAIL
            LT=POTPAR1(ISP)%LOFT(IT)
            IF(LT.EQ.L) THEN
              LMNP=LMNP0(LNP)
              DO M=1,2*LT+1
                LMNP=LMNP+1
                LMNT=LMNT+1
                PKJ(ISP)%PROJBAR(LMNP,LMNT)=POTPAR1(ISP)%PROJBAR(LNP,IT)
              ENDDO
            ELSE
              LMNT=LMNT+2*LT+1
            END IF
          ENDDO
        ENDDO
        DEALLOCATE(LMNP0)
      ENDDO
!
!     ==========================================================================
!     == PROJECTIONS ONTO LOCAL ORBITALS 
!     ==========================================================================
      NNS=SIZE(SBAR_NEW)
      ALLOCATE(PCHI(NNS))
      DO NN=1,NNS
!       == SVAR ENTERS TRANSPOSED. THE FIRST ATOM OF SBAR IS THE SECOND OF PCHI.
        IAT1=SBAR_NEW(NN)%IAT2
        IAT2=SBAR_NEW(NN)%IAT1
        ISP1=ISPECIES(IAT1)
        ISP2=ISPECIES(IAT2)
        LMNPX=SUM(2*LOX(:LNX(ISP1),ISP1)+1)
        LMNHX=SUM(2*POTPAR1(ISP2)%LOFH+1)
        ALLOCATE(PCHI(NN)%MAT(LMNPX,LMNHX))
        PCHI(NN)%IAT1=IAT1
        PCHI(NN)%IAT2=IAT2
        PCHI(NN)%IT=-SBAR_NEW(NN)%IT
!       == <P_I|CHI_J>=<P_I|PHIK_J>-SUM_K <P_I|PHIJBAR_K> SBAR_JK ==============
        PCHI(NN)%MAT=-MATMUL(PKJ(ISP1)%PROJBAR,TRANSPOSE(SBAR_NEW(NN)%MAT))
        IF(IAT1.EQ.IAT2.AND.SUM(ABS(SBAR_NEW(NN)%IT)).EQ.0) THEN
          PCHI(NN)%MAT=PCHI(NN)%MAT+PKJ(ISP1)%PROK
        END IF
      ENDDO
!
!     ==========================================================================
!     == CLEAN UP
!     ==========================================================================
      DO ISP=1,NSP
        DEALLOCATE(PKJ(ISP)%PROK)
        DEALLOCATE(PKJ(ISP)%PROJBAR)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$PROJTONTBO_NEW(ID,XK,NDIM,NBH,NPRO,PROJ,NORB,PIPSI)
!     **************************************************************************
!     ** TRANSFORMS THE PROJECTIONS ONTO COEFFICIENTS FOR                     **
!     ** NATURAL TIGHT-BINDING ORBITALS                                       **
!     **                                                                      **
!     ** REMARK:                                                              **
!     **   NOTE THAT NPRO REFERS TO THE NUMBER OF PARTIAL WAVES (R,L,M,N),    **
!     **   WHILE NRL REFERS TO THE NUMBER OF DIFFERENT (R,L,M) SETS           **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : NSP      & ! #(ATOM TYPES)
     &                       ,ISPECIES & !(NAT)
     &                       ,POTPAR1  & !(NSP)
     &                       ,LNX      & !(NSP)
     &                       ,LOX      & !(LNXX,NSP) (MAIN ANGULAR MOMENTA)
     &                       ,PCHI
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID !MAY BE 'FWRD' OR 'BACK'
      REAL(8)     ,INTENT(IN) :: XK(3) ! # K-POINT IN RELATIVE COORDINATES
      INTEGER(4)  ,INTENT(IN) :: NDIM
      INTEGER(4)  ,INTENT(IN) :: NBH
      INTEGER(4)  ,INTENT(IN) :: NPRO  ! #(PARTIAL WAVES)
      INTEGER(4)  ,INTENT(IN) :: NORB  ! #(LOCAL ORBITALS)
      COMPLEX(8),INTENT(INOUT):: PROJ(NDIM,NBH,NPRO)
      COMPLEX(8),INTENT(INOUT):: PIPSI(NDIM,NBH,NORB)
      REAL(8)   ,PARAMETER    :: PI=4.D0*ATAN(1.D0)
      LOGICAL(4),PARAMETER    :: TPR=.FALSE.
      COMPLEX(8),ALLOCATABLE  :: PCHIOFK(:,:)
      COMPLEX(8),ALLOCATABLE  :: OPCHIOFK(:,:)
      COMPLEX(8),ALLOCATABLE  :: MAT(:,:)
      COMPLEX(8),ALLOCATABLE  :: MATIN(:,:)
      INTEGER(4)              :: NAT                !#(ATOMS)
      INTEGER(4),ALLOCATABLE  :: NPRO1(:),NPRO2(:)  !(NAT)
      INTEGER(4),ALLOCATABLE  :: NORB1(:),NORB2(:)  !(NAT)
      INTEGER(4)              :: NNS
      REAL(8)                 :: KR
      COMPLEX(8)              :: EIKR
      REAL(8)                 :: SVAR  !SUPPORT VARIABLE
      INTEGER(4)              :: IPRO,IORB
      INTEGER(4)              :: IAT1,IAT2,I1,I2,J1,J2,I,J
      INTEGER(4)              :: ISP,NN,IAT,LN1,LN2,LMN1,LMN2,M,IBH,IDIM,L
!     **************************************************************************
      NAT=SIZE(ISPECIES)
!
!     ==========================================================================
!     == DETERMINE START AND END OF THE INDICES FOR A GIVEN ATOM
!     ==========================================================================
CALL TIMING$CLOCKON('LMTO:TO-1')
      ALLOCATE(NPRO1(NAT))
      ALLOCATE(NPRO2(NAT))
      ALLOCATE(NORB1(NAT))
      ALLOCATE(NORB2(NAT))
      IPRO=0
      IORB=0
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        NPRO1(IAT)=IPRO+1
        NORB1(IAT)=IORB+1
        IPRO=IPRO+SUM(2*LOX(:LNX(ISP),ISP)+1)
        IORB=IORB+SUM(2*POTPAR1(ISP)%LOFH+1)
        NPRO2(IAT)=IPRO
        NORB2(IAT)=IORB
      ENDDO
CALL TIMING$CLOCKOFF('LMTO:TO-1')
!
!     ==========================================================================
!     ==  BUILD <PRO|CHI> IN K-SPACE
!     ==========================================================================
CALL TIMING$CLOCKON('LMTO:TO-2')
      ALLOCATE(PCHIOFK(NPRO,NORB))
      PCHIOFK(:,:)=(0.D0,0.D0)
      NNS=SIZE(PCHI)
      DO NN=1,NNS 
        IAT1=PCHI(NN)%IAT1
        IAT2=PCHI(NN)%IAT2
        KR=2.D0*PI*SUM(REAL(PCHI(NN)%IT(:))*XK(:))
        EIKR=CMPLX(COS(KR),SIN(KR),KIND=8) 
        I1=NPRO1(IAT1)
        I2=NPRO2(IAT1)
        J1=NORB1(IAT2)
        J2=NORB2(IAT2)
        PCHIOFK(I1:I2,J1:J2)=PCHIOFK(I1:I2,J1:J2)+PCHI(NN)%MAT*EIKR
      ENDDO
CALL TIMING$CLOCKOFF('LMTO:TO-2')
!
!     ==========================================================================
!     == CALCULAT OPCHIOFK
!     ==========================================================================
CALL TIMING$CLOCKON('LMTO:TO-3')
      ALLOCATE(OPCHIOFK(NPRO,NORB))
      OPCHIOFK=(0.D0,0.D0)
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        DO LN1=1,LNX(ISP)
          L=LOX(LN1,ISP)
          DO LN2=1,LNX(ISP)
            IF(LOX(LN2,ISP).NE.L) CYCLE
            LMN1=NPRO1(IAT)-1+SUM(2*LOX(:LN1-1,ISP)+1)
            LMN2=NPRO1(IAT)-1+SUM(2*LOX(:LN2-1,ISP)+1)
            SVAR=POTPAR1(ISP)%PHIOV(LN1,LN2)
            DO M=1,2*L+1
              LMN1=LMN1+1
              LMN2=LMN2+1
              OPCHIOFK(LMN1,:)=OPCHIOFK(LMN1,:)+SVAR*PCHIOFK(LMN2,:)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
CALL TIMING$CLOCKOFF('LMTO:TO-3')
!
!     ==========================================================================
!     == TRANSFORMATION FROM PARTIAL WAVE PROJECTIONS TO ORBITAL PROJECTIONS  ==
!     ==========================================================================
CALL TIMING$CLOCKON('LMTO:TO-4')
      ALLOCATE(MAT(NORB,NORB))
      MAT=MATMUL(CONJG(TRANSPOSE(PCHIOFK)),OPCHIOFK)
      ALLOCATE(MATIN(NORB,NORB))
      CALL LIB$INVERTC8(NORB,MAT,MATIN)
      DEALLOCATE(MAT)
      ALLOCATE(MAT(NORB,NPRO))
      MAT=MATMUL(MATIN,CONJG(TRANSPOSE(OPCHIOFK)))
      DEALLOCATE(MATIN)
      DEALLOCATE(PCHIOFK)
      DEALLOCATE(OPCHIOFK)
CALL TIMING$CLOCKOFF('LMTO:TO-4')
!
!     ==========================================================================
!     == TRANSFORM PARTIAL WAVE PROJECTIONS TO ORBITAL PROJECTIONS            ==
!     ==========================================================================
CALL TIMING$CLOCKON('LMTO:TO-5')
      IF(ID.EQ.'FWRD') THEN
        PIPSI(:,:,:)=(0.D0,0.D0)
        DO J=1,NPRO
          DO I=1,NORB
             PIPSI(:,:,I)=PIPSI(:,:,I)+MAT(I,J)*PROJ(:,:,J)
          ENDDO
        ENDDO
!
      ELSE IF(ID.EQ.'BACK') THEN
        PROJ(:,:,:)=(0.D0,0.D0)
        DO I=1,NORB
          DO J=1,NPRO
            PROJ(:,:,J)=PROJ(:,:,J)+PIPSI(:,:,I)*CONJG(MAT(I,J))
          ENDDO
        ENDDO
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED. (ALLOWED IS "FWRD" OR "BACK")')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO$PROJTONTBO_NEW')
      END IF
      DEALLOCATE(MAT)
      DEALLOCATE(NPRO1)
      DEALLOCATE(NPRO2)
      DEALLOCATE(NORB1)
      DEALLOCATE(NORB2)
CALL TIMING$CLOCKOFF('LMTO:TO-5')
!
!     ==========================================================================
!     == PRINT FOR TESTING                                                    ==
!     ==========================================================================
      IF(TPR) THEN
        IF(ID.EQ.'FWRD') THEN
          WRITE(*,FMT='(80("="),T10," PROJECTIONS (NEW) ")')
        ELSE
          WRITE(*,FMT='(80("="),T10," DE/D-PROJECTIONS (NEW)")')
        END IF
        DO IBH=1,NBH
          DO IDIM=1,NDIM
            WRITE(*,FMT='(80("-"),T10," IBH=",I5," IDIM=",I2,"  ")')IBH,IDIM
            WRITE(*,FMT='("RE:",10F20.5)')REAL(PROJ(IDIM,IBH,:))
            WRITE(*,FMT='("IM:",10F20.5)')AIMAG(PROJ(IDIM,IBH,:))
          ENDDO
        ENDDO
        IF(ID.EQ.'FWRD') THEN
          WRITE(*,FMT='(80("="),T10," TIGHT-BINDIG COEFFICIENTS (NEW) ")')
        ELSE
          WRITE(*,FMT='(80("="),T10," DE/D-TIGHT-BINDIG COEFFICIENTS (NEW) ")')
        END IF
        DO IBH=1,NBH
          DO IDIM=1,NDIM
            WRITE(*,FMT='(80("-"),T10," IBH=",I5," IDIM=",I2,"  ")')IBH,IDIM
            WRITE(*,FMT='("RE:",10F20.5)')REAL(PIPSI(IDIM,IBH,:))
            WRITE(*,FMT='("IM:",10F20.5)')AIMAG(PIPSI(IDIM,IBH,:))
          ENDDO
        ENDDO
!!$        CALL ERROR$MSG('FORCED STOP AFTER WRITING DIAGNOSTIC RESULTS')
!!$        CALL ERROR$STOP('LMTO$PROJTONTBO_NEW')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$DEDF(NB,NKPTL,NSPIN,DEIG)
!     **************************************************************************
!     ** RETURNS THE CORRECTION FOR DEDF, WHICH IS NOT TAKEN CARE OF BY       **
!     ** <PI|HPSI>.                                                           **
!     **                                                                      **
!     ** REMARK: ACTS ONLY ON THE LOCAL SET OF K-POINTS                       **
!     ** REMARK: SET TO ZERO IF THERE IS NO CORRECTION!                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NKPTL
      INTEGER(4),INTENT(IN) :: NSPIN
      REAL(8)   ,INTENT(OUT):: DEIG(NB,NKPTL,NSPIN)
!     **************************************************************************
      DEIG=0.D0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBODENMAT_NEW()
!     **************************************************************************
!     **  CONSTRUCT DENSITY MATRIX IN A NTBO BASIS                            **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : DENMAT_NEW &
     &                       ,SBAR_NEW &
     &                       ,POTPAR1  &
     &                       ,LOX &
     &                       ,LNX &
     &                       ,ISPECIES
      USE MPE_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
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
COMPLEX(8)  :: PHASE
      INTEGER(4)             :: NTASKS,THISTASK,ICOUNT
!     **************************************************************************
                                          CALL TRACE$PUSH('LMTO_NTBODENMAT_NEW')
      IF(.NOT.ASSOCIATED(THIS%TBC_NEW)) THEN
        CALL ERROR$MSG('THIS%TBC_NEW NOT ASSOCIATED (1)')
        CALL ERROR$STOP('LMTO_NTBODENMAT_NEW')
      END IF
      IF(NDIM.EQ.1) THEN
        NDIMD=NSPIN
      ELSE IF(NDIM.EQ.2) THEN
        NDIMD=4
      END IF
!
!     ==========================================================================
!     == ALLOCATE DENSITY MATRIX
!     ==========================================================================
      NNS=SIZE(SBAR_NEW)
      NND=NNS
      ALLOCATE(DENMAT_NEW(NND))
      DO NN=1,NNS
        IAT1=SBAR_NEW(NN)%IAT1
        IAT2=SBAR_NEW(NN)%IAT2
        IT(:)=SBAR_NEW(NN)%IT(:)
        ISP=ISPECIES(IAT1)
        N1=SUM(2*POTPAR1(ISP)%LOFH+1)
        ISP=ISPECIES(IAT2)
        N2=SUM(2*POTPAR1(ISP)%LOFH+1)
        DENMAT_NEW(NN)%IAT1=IAT1
        DENMAT_NEW(NN)%IAT2=IAT2
        DENMAT_NEW(NN)%IT=IT
        DENMAT_NEW(NN)%N1=N1
        DENMAT_NEW(NN)%N2=N2
        DENMAT_NEW(NN)%N3=NDIMD  !(TOTAL,X,Y,Z)
        ALLOCATE(DENMAT_NEW(NN)%MAT(N1,N2,NDIMD))
        DENMAT_NEW(NN)%MAT(:,:,:)=0.D0
      ENDDO
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
        NPROAT(IAT)=SUM(2*POTPAR1(ISP)%LOFH+1)
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
!
!     ==========================================================================
!     ==  ADD UP DENSITY MATRIX                                               ==
!     ==========================================================================
      NND=SIZE(DENMAT_NEW)
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
          DO NN=1,NND
            ICOUNT=ICOUNT+1
            IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) CYCLE
            IAT1=DENMAT_NEW(NN)%IAT1
            IAT2=DENMAT_NEW(NN)%IAT2
            IT  =DENMAT_NEW(NN)%IT
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
                    C1(:)=THIS%TBC_NEW(:,IBH,I0+I)
                    C2(:)=THIS%TBC_NEW(:,IBH,J0+J)*EIKR   ! EXP(-I*K*T)
                    DO JDIM=1,NDIM
                      CSVAR22(:,JDIM)=CSVAR22(:,JDIM) &
     &                               +0.5D0*((F1+F2)*C1(:)*CONJG(C2(JDIM)) &
     &                                      +(F1-F2)*C1(:)*C2(JDIM))
                    ENDDO
                  ENDDO
                  CSVAR22=REAL(CSVAR22) ! IMAG(CSVAR) CONTAINS CRAP 
                                        !  DUE TO SUPER WAVE FUNCTIONS
                ELSE
                  CSVAR22=(0.D0,0.D0)
                  DO IB=1,NB
                    F1=OCC(IB,IKPT,ISPIN)
                    C1(:)=THIS%TBC_NEW(:,IB,I0+I)
                    C2(:)=THIS%TBC_NEW(:,IB,J0+J)*EIKR ! EXP(-I*K*T)
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
                    DENMAT_NEW(NN)%MAT(I,J,1)=DENMAT_NEW(NN)%MAT(I,J,1) &
     &                                       +REAL(CSVAR22(1,1))
                  ELSE ! NONCOLLINEAR
                    DENMAT_NEW(NN)%MAT(I,J,1)=DENMAT_NEW(NN)%MAT(I,J,1) &
     &                                  +REAL(CSVAR22(1,1)+CSVAR22(2,2))
                    DENMAT_NEW(NN)%MAT(I,J,2)=DENMAT_NEW(NN)%MAT(I,J,2) &
     &                                  +REAL(CSVAR22(1,2)+CSVAR22(2,1))
                    DENMAT_NEW(NN)%MAT(I,J,3)=DENMAT_NEW(NN)%MAT(I,J,3) &
     &                                  -AIMAG(CSVAR22(1,2)-CSVAR22(2,1))
                    DENMAT_NEW(NN)%MAT(I,J,4)=DENMAT_NEW(NN)%MAT(I,J,4) &
     &                                  +REAL(CSVAR22(1,1)-CSVAR22(2,2))
                  END IF
                ELSE IF(NSPIN.EQ.2) THEN !SPIN POLARIZED
                  IF(ISPIN.EQ.1) THEN
                    DENMAT_NEW(NN)%MAT(I,J,1)=DENMAT_NEW(NN)%MAT(I,J,1) &
     &                                       +REAL(CSVAR22(1,1))
                    DENMAT_NEW(NN)%MAT(I,J,2)=DENMAT_NEW(NN)%MAT(I,J,2) &
     &                                       +REAL(CSVAR22(1,1))
                  ELSE
                    DENMAT_NEW(NN)%MAT(I,J,1)=DENMAT_NEW(NN)%MAT(I,J,1) &
     &                                       +REAL(CSVAR22(1,1))
                    DENMAT_NEW(NN)%MAT(I,J,2)=DENMAT_NEW(NN)%MAT(I,J,2) &
     &                                       -REAL(CSVAR22(1,1))
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
        CALL MPE$COMBINE('MONOMER','+',DENMAT_NEW(NN)%MAT)
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,FMT='(82("="),T10," NON-LOCAL DENSITY MATRIX IN A NTBO BASIS ")')
        WRITE(*,FMT='(82("="),T10," FROM LMTO_NTBODENMAT_NEW  ")')
        DO NN=1,NND
          IAT1=DENMAT_NEW(NN)%IAT1
          IAT2=DENMAT_NEW(NN)%IAT2
          IT=DENMAT_NEW(NN)%IT
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
     &                                                         IAT1,IAT2,IT
          N1=DENMAT_NEW(NN)%N1
          N2=DENMAT_NEW(NN)%N2
          DO I=1,1 !DENMAT_NEW(NN)%N3
            DO J=1,N1 
              WRITE(*,FMT='(I3,300F10.3)')I,DENMAT_NEW(NN)%MAT(J,:,I)
            ENDDO
            WRITE(*,FMT='(82("-"))')
          ENDDO
        ENDDO
        CALL ERROR$MSG('FORCED STOP AFTER PRINTING DENSITY MATRIX')
        CALL ERROR$STOP('LMTO_NTBODENMAT_NEW')
      END IF
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBOINVOVERLAP()
!     **************************************************************************
!     **  CONSTRUCT THE INVERSE OVERLAP FROM THE NTBO PROJECTIONS AS          **
!     **    SUM_N <PI_A|\PSI_N><\PSI_N|\PI_B>                                 **
!     **  THE OVERLAP CAN LATER BE OBTAINED BY INVERSION IN K-SPACE           **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2014******************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : INVOVERLAP &
     &                       ,SBAR_NEW &
     &                       ,POTPAR1  &
     &                       ,ISPECIES
      USE MPE_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
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
      REAL(8)   ,ALLOCATABLE :: WKPTL(:)
      INTEGER(4)             :: IAT,NN,II,ISP,IPRO,IKPT,ISPIN,I,J,IBH,IB
      REAL(8)                :: SVAR
      REAL(8)                :: F1,F2
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IAT1,IAT2,IT(3),I0,J0,IDIM,JDIM
      COMPLEX(8)             :: EIKR,C1(NDIM),C2(NDIM),CSVAR22(NDIM,NDIM)
COMPLEX(8)  :: PHASE
      INTEGER(4)             :: NTASKS,THISTASK,ICOUNT
!     **************************************************************************
                                          CALL TRACE$PUSH('LMTO_NTBOINVOVERLAP')
      IF(.NOT.ASSOCIATED(THIS%TBC_NEW)) THEN
        CALL ERROR$MSG('THIS%TBC_NEW NOT ASSOCIATED (1)')
        CALL ERROR$STOP('LMTO_NTBOINVOVERLAP')
      END IF
      IF(NDIM.EQ.1) THEN
        NDIMD=NSPIN
      ELSE IF(NDIM.EQ.2) THEN
        NDIMD=4
      END IF
!
!     ==========================================================================
!     == ALLOCATE DENSITY MATRIX
!     ==========================================================================
      NNS=SIZE(SBAR_NEW)
      NND=NNS
      IF(ALLOCATED(INVOVERLAP)) THEN
        CALL ERROR$MSG('INVOVERLAP ALREADY ALLOCATED')
        CALL ERROR$STOP('LMTO_NTBOINVOVERLAP')
      END IF
      ALLOCATE(INVOVERLAP(NND))
      DO NN=1,NNS
        IAT1=SBAR_NEW(NN)%IAT1
        IAT2=SBAR_NEW(NN)%IAT2
        IT(:)=SBAR_NEW(NN)%IT(:)
        ISP=ISPECIES(IAT1)
        N1=SUM(2*POTPAR1(ISP)%LOFH+1)
        ISP=ISPECIES(IAT2)
        N2=SUM(2*POTPAR1(ISP)%LOFH+1)
        INVOVERLAP(NN)%IAT1=IAT1
        INVOVERLAP(NN)%IAT2=IAT2
        INVOVERLAP(NN)%IT=IT
        INVOVERLAP(NN)%N1=N1
        INVOVERLAP(NN)%N2=N2
        ALLOCATE(INVOVERLAP(NN)%MAT(N1,N2))
        INVOVERLAP(NN)%MAT(:,:)=0.D0
      ENDDO
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      CALL DYNOCC$GETI4('NB',NBX)
      ALLOCATE(WKPTL(NKPTL))
      CALL WAVES_DYNOCCGETR8A('WKPT',NKPTL,WKPTL)
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
        NPROAT(IAT)=SUM(2*POTPAR1(ISP)%LOFH+1)
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
!
!     ==========================================================================
!     ==  ADD UP DENSITY MATRIX                                               ==
!     ==========================================================================
      NND=SIZE(INVOVERLAP)
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
          DO NN=1,NND
            ICOUNT=ICOUNT+1
            IF(MOD(ICOUNT-1,NTASKS).NE.THISTASK-1) CYCLE
            IAT1=INVOVERLAP(NN)%IAT1
            IAT2=INVOVERLAP(NN)%IAT2
            IT  =INVOVERLAP(NN)%IT
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
                    C1(:)=THIS%TBC_NEW(:,IBH,I0+I)
                    C2(:)=THIS%TBC_NEW(:,IBH,J0+J)*EIKR   ! EXP(-I*K*T)
                    DO JDIM=1,NDIM
                      CSVAR22(:,JDIM)=CSVAR22(:,JDIM)+C1(:)*CONJG(C2(JDIM)) 
                    ENDDO
                  ENDDO
                  CSVAR22=WKPTL(IKPT)*CSVAR22
                  CSVAR22=REAL(CSVAR22) ! IMAG(CSVAR) CONTAINS CRAP 
                                        !  DUE TO SUPER WAVE FUNCTIONS
                ELSE
                  CSVAR22=(0.D0,0.D0)
                  DO IB=1,NB
                    C1(:)=THIS%TBC_NEW(:,IB,I0+I)
                    C2(:)=THIS%TBC_NEW(:,IB,J0+J)*EIKR ! EXP(-I*K*T)
                    DO JDIM=1,NDIM
                      CSVAR22(:,JDIM)=CSVAR22(:,JDIM)+C1(:)*CONJG(C2(JDIM))
                    ENDDO
                  ENDDO
                  CSVAR22=WKPTL(IKPT)*CSVAR22
                END IF
!
!           == DISTRIBUTE ONTO DENSITY MATRIX ENTRIES ==========================
!           == D(IDIMD)=SUM_{IDIM,JDIM} D(IDIM,JDIM)*PAULI_{IDIMD}(JDIM,IDIM) ==
!           == IDIMD IN {TOTAL,X,Y,Z}; IDIM IN {UP,DOWN} =======================
!           == TRANSFORMATION MUST BE CONSISTENT WITH WAVES_DENMAT =============
                IF(NSPIN.EQ.1) THEN
                  IF(NDIM.EQ.1) THEN !NON-SPIN-POLARIZED
                    INVOVERLAP(NN)%MAT(I,J)=INVOVERLAP(NN)%MAT(I,J) &
     &                                       +REAL(CSVAR22(1,1))
                  ELSE ! NONCOLLINEAR
                    INVOVERLAP(NN)%MAT(I,J)=INVOVERLAP(NN)%MAT(I,J) &
     &                                      +REAL(CSVAR22(1,1)+CSVAR22(2,2))
                  END IF
                ELSE IF(NSPIN.EQ.2) THEN !SPIN POLARIZED
                  IF(ISPIN.EQ.1) THEN
                    INVOVERLAP(NN)%MAT(I,J)=INVOVERLAP(NN)%MAT(I,J) &
     &                                       +REAL(CSVAR22(1,1))
                  ELSE
                    INVOVERLAP(NN)%MAT(I,J)=INVOVERLAP(NN)%MAT(I,J) &
     &                                       +REAL(CSVAR22(1,1))
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
        CALL MPE$COMBINE('MONOMER','+',INVOVERLAP(NN)%MAT)
      ENDDO
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,FMT='(82("="),T10," INVERSE OVERLAP MATRIX IN A NTBO BASIS ")')
        WRITE(*,FMT='(82("="),T10," FROM LMTO_NTBOINVOVERLAP  ")')
        DO NN=1,NND
          IAT1=INVOVERLAP(NN)%IAT1
          IAT2=INVOVERLAP(NN)%IAT2
          IT=INVOVERLAP(NN)%IT
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
       &                                                       IAT1,IAT2,IT
          N1=INVOVERLAP(NN)%N1
          N2=INVOVERLAP(NN)%N2
          DO J=1,N1 
            WRITE(*,FMT='(I3,300F10.3)')I,INVOVERLAP(NN)%MAT(J,:)
          ENDDO
        ENDDO
        CALL ERROR$MSG('FORCED STOP AFTER PRINTING INVERSE OVERLAP MATRIX')
        CALL ERROR$STOP('LMTO_NTBOINVOVERLAP')
      END IF
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBODENMATDER_NEW()
!     **************************************************************************
!     **  CONSTRUCT HTBC FROM THE DERIVATIVE OF THE ENERGY WITH RESPECT TO    **
!     **  THE DENSITY MATRIX                                                  **
!     **                                                                      **
!     **************************************************************************
      USE WAVES_MODULE, ONLY: WVSET_TYPE,NKPTL,NSPIN,NDIM &
     &                       ,THIS &
     &                       ,MAP &
     &                       ,WAVES_SELECTWV &
     &                       ,GSET
      USE LMTO_MODULE, ONLY : HAMIL=>HAMIL_NEW &
     &                       ,POTPAR1 &
     &                       ,LOX,LNX,ISPECIES &
     &                       ,THTBC !LOGICAL VARIABLE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
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
!     **************************************************************************
                                 CALL TRACE$PUSH('LMTO_NTBODENMATDER_NEW')
      IF(NDIM.EQ.1) THEN
        NDIMD=NSPIN
      ELSE IF(NDIM.EQ.2) THEN
        NDIMD=4
      END IF
      THTBC=.TRUE.   ! HTBC WILL BE CALCULATED
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
        NPROAT(IAT)=SUM(2*POTPAR1(ISP)%LOFH+1)
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
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          NBH=THIS%NBH
          NB=THIS%NB
          IF(.NOT.ASSOCIATED(THIS%HTBC_NEW)) &
     &              ALLOCATE(THIS%HTBC_NEW(NDIM,NBH,NPRO))
          THIS%HTBC_NEW(:,:,:)=(0.D0,0.D0)
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
                    CSVAR22(1,1)=CMPLX(HAMIL(NN)%MAT(I,J,1),0.D0,KIND=8)
                  ELSE ! NONCOLLINEAR
                    CSVAR22(1,1)=CMPLX(HAMIL(NN)%MAT(I,J,1) &
     &                                +HAMIL(NN)%MAT(I,J,4),0.D0,KIND=8)
                    CSVAR22(1,2)=CMPLX(HAMIL(NN)%MAT(I,J,2) &
     &                               ,-HAMIL(NN)%MAT(I,J,3),KIND=8)
                    CSVAR22(2,1)=CMPLX(HAMIL(NN)%MAT(I,J,2) &
     &                               ,+HAMIL(NN)%MAT(I,J,3),KIND=8)
                    CSVAR22(2,2)=CMPLX(HAMIL(NN)%MAT(I,J,1) &
     &                                -HAMIL(NN)%MAT(I,J,4),0.D0,KIND=8)
                  END IF
                ELSE IF(NSPIN.EQ.2) THEN !SPIN POLARIZED
                  IF(ISPIN.EQ.1) THEN
                    CSVAR22(1,1)=CMPLX(HAMIL(NN)%MAT(I,J,1) &
     &                                +HAMIL(NN)%MAT(I,J,2),0.D0,KIND=8)
                  ELSE
                    CSVAR22(1,1)=CMPLX(HAMIL(NN)%MAT(I,J,1) &
     &                                -HAMIL(NN)%MAT(I,J,2),0.D0,KIND=8)
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
                      THIS%HTBC_NEW(IDIM,IB,I0+I)=THIS%HTBC_NEW(IDIM,IB,I0+I) &
      &                           +CSVAR22(IDIM,JDIM)*THIS%TBC_NEW(JDIM,IB,J0+J)
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
!     == PRINT FOR TESTING                                                    ==
!     ==========================================================================
      IF(TPR) THEN
        DO IKPT=1,NKPTL
          DO ISPIN=1,NSPIN
            CALL WAVES_SELECTWV(IKPT,ISPIN)
            NBH=THIS%NBH
            WRITE(*,FMT='(80("="),T10," HTBC(OLD) FOR XK=",3F10.5,"  ")') &
     &                                                                XK(:,IKPT)
            DO IBH=1,NBH
              DO IDIM=1,NDIM
                WRITE(*,FMT='(80("-"),T10," IBH=",I5," IDIM=",I2,"  ")')IBH,IDIM
                WRITE(*,FMT='("RE:",10F20.5)')REAL(THIS%HTBC_NEW(IDIM,IBH,:))
                WRITE(*,FMT='("IM:",10F20.5)')AIMAG(THIS%HTBC_NEW(IDIM,IBH,:))
              ENDDO
            ENDDO
            WRITE(*,FMT='(80("="),T10," HTBC(NEW) FOR XK=",3F10.5,"  ")') &
     &                                                                XK(:,IKPT)
            DO IBH=1,NBH
              DO IDIM=1,NDIM
                WRITE(*,FMT='(80("-"),T10," IBH=",I5," IDIM=",I2,"  ")')IBH,IDIM
                WRITE(*,FMT='("RE:",10F20.5)')REAL(THIS%HTBC_NEW(IDIM,IBH,:))
                WRITE(*,FMT='("IM:",10F20.5)')AIMAG(THIS%HTBC_NEW(IDIM,IBH,:))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CALL ERROR$MSG('FORCED STOP AFTER WRITING DIAGNOSTIC RESULTS')
        CALL ERROR$STOP('LMTO_NTBODENMATDER_NEW')
      END IF
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SETHTBCTOZERO()
!     **************************************************************************
!     **************************************************************************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,THIS,WAVES_SELECTWV
      USE LMTO_MODULE, ONLY: THTBC
      IMPLICIT NONE
      INTEGER(4) :: IKPT,ISPIN
!     **************************************************************************
      THTBC=.FALSE. ! HTBC WILL BE RESET SET TO ZERO
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          IF(ASSOCIATED(THIS%HTBC_NEW))THIS%HTBC_NEW(:,:,:)=(0.D0,0.D0)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_CLEANDENMAT_NEW()
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : DENMAT=>DENMAT_NEW &
     &                       ,HAMIL=>HAMIL_NEW &
     &                       ,INVOVERLAP &
     &                       ,DEDOI
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
!
!     ==========================================================================
!     ==  CLEAN UP INVOVERLAP
!     ==========================================================================
      IF(ALLOCATED(INVOVERLAP)) THEN
        NND=SIZE(INVOVERLAP)
        DO NN=1,NND
          DEALLOCATE(INVOVERLAP(NN)%MAT)
        ENDDO
        DEALLOCATE(INVOVERLAP)
      END IF
!
!     ==========================================================================
!     ==  CLEAN UP DEDOI
!     ==========================================================================
      IF(ALLOCATED(DEDOI)) THEN
        NND=SIZE(DEDOI)
        DO NN=1,NND
          DEALLOCATE(DEDOI(NN)%MAT)
        ENDDO
        DEALLOCATE(DEDOI)
      END IF
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE LMTOAUGMENTATION_MODULE
!*******************************************************************************
!** THE LMTOAUGMENTATION IS A SUBOBJECT OF LMTO, WHICH IS USED TO CONSTRUCT   **
!** PART OF THE DOUBLE COUNTING. IN THE DOUBLE COUNTING THE TOTAL ONE-CENTER  **
!** EXPANSION OF THE DENSITY ENTERS, WHICH MAY DIFFER FROM THE DENSITY OF THE **
!** LOCAL ORBITALS.                                                           **
!**                                                                           **
!** TACTIVE: TACTIVE IS SET TO TRUE WHEN THE DENSITY OF STATES IS SUPPLIED    **
!**          AND IT IS SET TO FALSE BY LMTOAUGMENTATION$CLEAN. DURING ONE     **
!**          SUCH LMTOAGMENTATION$ADD... IS ADDING                            **
!**                                                                           **
!*******************************************************************************
IMPLICIT NONE
LOGICAL(4),SAVE        :: TINI=.FALSE.
LOGICAL(4),SAVE        :: TACTIVE=.FALSE.
INTEGER(4)             :: NAT=0
INTEGER(4)             :: IAT=0           ! ATOM FOCUS
INTEGER(4)             :: NDIMD=0
INTEGER(4)             :: LMNXX=0
COMPLEX(8),ALLOCATABLE :: DENMAT(:,:,:,:) !(LMNXX,LMNXX,NDIMD,NAT)
COMPLEX(8),ALLOCATABLE :: DATH(:,:,:,:)   !(LMNXX,LMNXX,NDIMD,NAT)
END MODULE LMTOAUGMENTATION_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTOAUGMENTATION_INITIALIZE()
!     **************************************************************************
!     **************************************************************************
      USE LMTOAUGMENTATION_MODULE, ONLY : TINI,NAT,NDIMD,LMNXX
      IMPLICIT NONE
!     **************************************************************************
      IF(TINI) RETURN
      CALL ATOMLIST$NATOM(NAT)
      IF(NDIMD.EQ.0) THEN
        CALL ERROR$MSG('NDIMD HAS NOT BEEN SET')
        CALL ERROR$STOP('LMTOAUGMENTATION_INITIALIZE')
      END IF
      IF(LMNXX.EQ.0) THEN
        CALL ERROR$MSG('LMNXX HAS NOT BEEN SET')
        CALL ERROR$STOP('LMTOAUGMENTATION_INITIALIZE')
      END IF
      TINI=.TRUE.
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTOAUGMENTATION$CLEAN()
!     **************************************************************************
!     ** FREE MEMORY. NEEDS TO BE ACTIVATED BY SETTING A DENSITY MATRIX
!     **************************************************************************
      USE LMTOAUGMENTATION_MODULE, ONLY : TACTIVE,IAT,DENMAT,DATH
      IMPLICIT NONE
!     **************************************************************************
      IF(ALLOCATED(DENMAT))DEALLOCATE(DENMAT)
      IF(ALLOCATED(DATH))  DEALLOCATE(DATH)
      IAT=0                      ! REMOVE ATOM FOCUS
      TACTIVE=.FALSE.

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTOAUGMENTATION$GETI4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTOAUGMENTATION_MODULE, ONLY : TINI,NAT,LMNXX,NDIMD,IAT
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ID.EQ.'NAT') THEN
        VAL=NAT   !USE ONLY TO INSPECT THIS OBJECT
      ELSE IF(ID.EQ.'IAT') THEN
        VAL=IAT   !USE ONLY TO INSPECT THIS OBJECT
      ELSE IF(ID.EQ.'LMNXX') THEN
        VAL=LMNXX   !USE ONLY TO INSPECT THIS OBJECT
      ELSE IF(ID.EQ.'NDIMD') THEN
        VAL=NDIMD !USE ONLY TO INSPECT THIS OBJECT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTOAUGMENTATION$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTOAUGMENTATION$SETI4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTOAUGMENTATION_MODULE, ONLY : TINI,NAT,LMNXX,NDIMD,IAT
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'NAT') THEN
        IF(TINI.AND.NAT.NE.VAL) THEN
          CALL ERROR$MSG('NAT MUST NOT BE CHANGED AFTER INITIALIZATION')
          CALL ERROR$STOP('LMTOAUGMENTATION$SETI4')
        END IF
        NAT=VAL
      ELSE IF(ID.EQ.'IAT') THEN
        IAT=VAL
      ELSE IF(ID.EQ.'LMNXX') THEN
        IF(TINI.AND.LMNXX.NE.VAL) THEN
          CALL ERROR$MSG('LMNXX MUST NOT BE CHANGED AFTER INITIALIZATION')
          CALL ERROR$STOP('LMTOAUGMENTATION$SETI4')
        END IF
        LMNXX=VAL
      ELSE IF(ID.EQ.'NDIMD') THEN
        IF(TINI.AND.NDIMD.NE.VAL) THEN
          CALL ERROR$MSG('NDIMD MUST NOT BE CHANGED AFTER INITIALIZATION')
          CALL ERROR$STOP('LMTOAUGMENTATION$SETI4')
        END IF
        NDIMD=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTOAUGMENTATION$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTOAUGMENTATION$SETC8A(ID,LEN,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTOAUGMENTATION_MODULE, ONLY : TINI,TACTIVE &
     &                                   ,NAT,LMNXX,NDIMD,DENMAT,DATH
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      COMPLEX(8)  ,INTENT(IN) :: VAL(LEN)
      INTEGER(4)              :: I
!     **************************************************************************
      IF(.NOT.TINI) CALL LMTOAUGMENTATION_INITIALIZE()
      IF(ID.EQ.'DENMAT') THEN
        IF(LMNXX*LMNXX*NDIMD*NAT.NE.LEN) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('LMNXX',LMNXX)
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTOAUGMENTATION$SETC8A')
        END IF
        IF(TACTIVE) THEN
          CALL ERROR$MSG('DENSITY MATRIX MUST NOT BE CHANGED IN ACTIVE STATE')
          CALL ERROR$MSG('CALL LMTOAUGMENTATION$CLEAN TO RESET OBJECT')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTOAUGMENTATION$SETC8A')
        END IF 
        ALLOCATE(DENMAT(LMNXX,LMNXX,NDIMD,NAT))
        ALLOCATE(DATH(LMNXX,LMNXX,NDIMD,NAT))
        DATH=(0.D0,0.D0)
        DENMAT=RESHAPE(VAL,(/LMNXX,LMNXX,NDIMD,NAT/))
        TACTIVE=.TRUE.
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTOAUGMENTATION$SETC8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTOAUGMENTATION$ADDC8A(ID,LEN,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTOAUGMENTATION_MODULE, ONLY : TINI,TACTIVE,IAT &
     &                                   ,NAT,LMNXX,NDIMD,DENMAT,DATH
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      COMPLEX(8)  ,INTENT(IN) :: VAL(LEN)
      INTEGER(4)              :: ISP,LNX,LMNX
      INTEGER(4)  ,ALLOCATABLE:: LOX(:)
      INTEGER(4)              :: I
!     **************************************************************************
      IF(.NOT.TINI) CALL LMTOAUGMENTATION_INITIALIZE()
      IF(ID.EQ.'DH') THEN
        IF(.NOT.TACTIVE) THEN
          CALL ERROR$MSG('LMTOAUGMENTATION OBJECT IS NOT ACTIVE')
          CALL ERROR$STOP('LMTOAUGMENTATION$ADDC8A')
        END IF 
        IF(IAT.EQ.0) THEN
          CALL ERROR$MSG('NO ATOM FOCUS')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTOAUGMENTATION$ADDC8A')
        END IF
        CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$UNSELECT()
        LMNX=SUM(2*LOX+1)
        IF(LMNX*LMNX*NDIMD.NE.LEN) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$I4VAL('LMNX',LMNX)
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTOAUGMENTATION$ADDC8A')
        END IF
        DATH(:LMNX,:LMNX,:,IAT)=DATH(:LMNX,:LMNX,:,IAT) &
     &                                         +RESHAPE(VAL,(/LMNX,LMNX,NDIMD/))
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTOAUGMENTATION$ADDC8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTOAUGMENTATION$GETC8A(ID,LEN,VAL)
!     **************************************************************************
!     **************************************************************************
      USE LMTOAUGMENTATION_MODULE, ONLY : TINI,TACTIVE,IAT &
     &                                    ,NAT,LMNXX,NDIMD,DATH,DENMAT
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      COMPLEX(8)  ,INTENT(OUT):: VAL(LEN)
      INTEGER(4)              :: ISP,LNX,LMNX
      INTEGER(4)  ,ALLOCATABLE:: LOX(:)
      INTEGER(4)              :: I
!     **************************************************************************
      IF(.NOT.TINI) CALL LMTOAUGMENTATION_INITIALIZE()
      IF(ID.EQ.'DH') THEN
        IF(.NOT.TACTIVE) THEN
          CALL ERROR$MSG('LMTOAUGMENTATION OBJECT IS NOT ACTIVE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTOAUGMENTATION$GETC8A')
        END IF 
        IF(LMNXX*LMNXX*NDIMD*NAT.NE.LEN) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('LMNXX',LMNXX)
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTOAUGMENTATION$GETC8A')
        END IF
        VAL=RESHAPE(DATH,(/LMNXX*LMNXX*NDIMD*NAT/))
!
      ELSE IF(ID.EQ.'DENMAT') THEN
        IF(.NOT.TACTIVE) THEN
          CALL ERROR$MSG('LMTOAUGMENTATION OBJECT IS NOT ACTIVE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTOAUGMENTATION$GETC8A')
        END IF 
        IF(IAT.EQ.0) THEN  ! NO ATOM FOCUS: DENSITY MATRIX FOR ALL ATOMS IS PROVIDED
          IF(LMNXX*LMNXX*NDIMD*NAT.NE.LEN) THEN
            CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
            CALL ERROR$I4VAL('LMNXX',LMNXX)
            CALL ERROR$I4VAL('NDIMD',NDIMD)
            CALL ERROR$I4VAL('NAT',NAT)
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('LMTOAUGMENTATION$GETC8A')
          END IF
          VAL=RESHAPE(DENMAT,(/LMNXX*LMNXX*NDIMD*NAT/))
        ELSE 
          CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
          CALL SETUP$ISELECT(ISP)
          CALL SETUP$GETI4('LNX',LNX)
          ALLOCATE(LOX(LNX))
          CALL SETUP$GETI4A('LOX',LNX,LOX)
          CALL SETUP$UNSELECT()
          LMNX=SUM(2*LOX+1)
          IF(LMNX*LMNX*NDIMD.NE.LEN) THEN
            CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
            CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
            CALL ERROR$I4VAL('LMNX',LMNX)
            CALL ERROR$I4VAL('NDIMD',NDIMD)
            CALL ERROR$I4VAL('IAT',IAT)
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$STOP('LMTOAUGMENTATION$GETC8A')
          END IF
          VAL=RESHAPE(DENMAT(:LMNX,:LMNX,:,IAT),(/LMNX*LMNX*NDIMD/))
        END IF
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTOAUGMENTATION$GETC8A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTOAUGMENTATION$GETRHO(GID,NR,LMRX,NDIMD_,RHO)
!     **************************************************************************
!     ** VALENCE DENSITY FROM PARTIAL WAVES
!     **************************************************************************
      USE LMTOAUGMENTATION_MODULE, ONLY : TINI,TACTIVE,IAT &
     &                                   ,NAT,LMNXX,NDIMD,DENMAT
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: GID
      INTEGER(4),INTENT(IN)  :: NR
      INTEGER(4),INTENT(IN)  :: LMRX
      INTEGER(4),INTENT(IN)  :: NDIMD_
      REAL(8)   ,INTENT(OUT) :: RHO(NR,LMRX,NDIMD)
      INTEGER(4)             :: ISP
      INTEGER(4)             :: GID1
      INTEGER(4)             :: NR1
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)    !(NR,LNX)
      COMPLEX(8),ALLOCATABLE :: DENMAT1(:,:)  !(LMNX,LMNX)
      INTEGER(4)             :: IDIMD
!     **************************************************************************
      IF(.NOT.TINI) CALL LMTOAUGMENTATION_INITIALIZE()
      IF(IAT.EQ.0) THEN
        CALL ERROR$MSG('ATOM INDEX NOT SET')
        CALL ERROR$STOP('LMTOAUGMENTATION$GETRHO')
      END IF
      IF(.NOT.TACTIVE) THEN
        CALL ERROR$MSG('LMTOAUGMENTATION OBJECT IS NOT ACTIVE')
        CALL ERROR$STOP('LMTOAUGMENTATION$GETRHO')
      END IF 
      IF(NDIMD_.NE.NDIMD) THEN
        CALL ERROR$MSG('INCONSISTENT VALUES OF NDIMD')
        CALL ERROR$I4VAL('NDIMD ON INPUT',NDIMD_)
        CALL ERROR$I4VAL('NDIMD FROM MODULE',NDIMD)
        CALL ERROR$STOP('LMTOAUGMENTATION$GETRHO')
      END IF

      CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID1)
      IF(GID.NE.GID1) THEN
        CALL ERROR$MSG('INCONSISTENT VALUE OF GID')
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$I4VAL('ISP',ISP)
        CALL ERROR$I4VAL('GID-ON-INPUT',GID)
        CALL ERROR$I4VAL('SET-GID',GID1)
        CALL ERROR$STOP('LMTOAUGMENTATION$GETRHO')
      END IF
      CALL RADIAL$GETI4(GID,'NR',NR1)
      IF(NR.NE.NR1) THEN
        CALL ERROR$MSG('INCONSISTENT VALUE OF NR')
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$I4VAL('ISP',ISP)
        CALL ERROR$I4VAL('NR-ON-INPUT',NR)
        CALL ERROR$I4VAL('SET-NR',NR1)
        CALL ERROR$STOP('LMTOAUGMENTATION$GETRHO')
      END IF
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      ALLOCATE(AEPHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUP$UNSELECT()
      LMNX=SUM(2*LOX(:LNX)+1)
      ALLOCATE(DENMAT1(LMNX,LMNX))
      DO IDIMD=1,NDIMD
        DENMAT1(:,:)=DENMAT(:LMNX,:LMNX,IDIMD,IAT)
        CALL AUGMENTATION_RHO(NR,LNX,LOX,AEPHI,LMNX,DENMAT1,LMRX,RHO(:,:,IDIMD))
      ENDDO
      DEALLOCATE(AEPHI)
      DEALLOCATE(LOX)
      DEALLOCATE(DENMAT1)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTOAUGMENTATION$MPECOMBINEPOT()
!     **************************************************************************
!     ** PARALLELIZATION:
!     ** SUMS DATH OVER ALL PROCESSORS IN THE MONOMER GROUP                   **
!     ** ASSUMES THAT DATH OF A GIVEN ATOM IS ONLY CALCULATED ON ONE TASK     **
!     **************************************************************************
      USE LMTOAUGMENTATION_MODULE, ONLY : TINI,TACTIVE,DATH
      USE MPE_MODULE
      IMPLICIT NONE
!     **************************************************************************
      IF(.NOT.TINI) CALL LMTOAUGMENTATION_INITIALIZE()
      IF(.NOT.TACTIVE) THEN
        CALL ERROR$MSG('LMTOAUGMENTATION OBJECT IS NOT ACTIVE')
        CALL ERROR$STOP('LMTOAUGMENTATION$MPECOMBINEPOT')
      END IF 
      CALL MPE$COMBINE('MONOMER','+',DATH)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTOAUGMENTATION$ADDPOT(GID,NR,LMRX,NDIMD_,POT)
!     **************************************************************************
!     ** CALCULATES MATRIX ELEMENT OF THE POTENTIAL WITH PARTIAL WAVES        **
!     ** AND ADDS THEM TO THE MODULE-ARRAY DATH                               **
!     **************************************************************************
      USE LMTOAUGMENTATION_MODULE, ONLY : TINI,TACTIVE,NDIMD,DATH,IAT
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: GID
      INTEGER(4),INTENT(IN)  :: NR
      INTEGER(4),INTENT(IN)  :: LMRX
      INTEGER(4),INTENT(IN)  :: NDIMD_
      REAL(8)   ,INTENT(IN)  :: POT(NR,LMRX,NDIMD)
      INTEGER(4)             :: ISP
      INTEGER(4)             :: GID1
      INTEGER(4)             :: NR1
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)    !(NR,LNX)
      REAL(8)                :: AECORE(NR)
      COMPLEX(8),ALLOCATABLE :: DATH1(:,:,:)  !(LMNX,LMNX)
!     **************************************************************************
      IF(.NOT.TINI) CALL LMTOAUGMENTATION_INITIALIZE()
      IF(IAT.EQ.0) THEN
        CALL ERROR$MSG('ATOM INDEX NOT SET')
        CALL ERROR$STOP('LMTOAUGMENTATION$ADDPOT')
      END IF
      IF(.NOT.TACTIVE) THEN
        CALL ERROR$MSG('LMTOAUGMENTATION OBJECT IS NOT ACTIVE')
        CALL ERROR$STOP('LMTOAUGMENTATION$ADDPOT')
      END IF 
      IF(NDIMD_.NE.NDIMD) THEN
        CALL ERROR$MSG('INCONSISTENT VALUES OF NDIMD')
        CALL ERROR$I4VAL('NDIMD ON INPUT',NDIMD_)
        CALL ERROR$I4VAL('NDIMD FROM MODULE',NDIMD)
        CALL ERROR$STOP('LMTOAUGMENTATION$ADDPOT')
      END IF
      CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID1)
      IF(GID.NE.GID1) THEN
        CALL ERROR$MSG('INCONSISTENT VALUE OF GID')
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$I4VAL('ISP',ISP)
        CALL ERROR$I4VAL('GID-ON-INPUT',GID)
        CALL ERROR$I4VAL('SET-GID',GID1)
        CALL ERROR$STOP('LMTOAUGMENTATION$ADDPOT')
      END IF
      CALL RADIAL$GETI4(GID,'NR',NR1)
      IF(NR.NE.NR1) THEN
        CALL ERROR$MSG('INCONSISTENT VALUE OF NR')
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$I4VAL('ISP',ISP)
        CALL ERROR$I4VAL('NR-ON-INPUT',NR)
        CALL ERROR$I4VAL('SET-NR',NR1)
        CALL ERROR$STOP('LMTOAUGMENTATION$ADDPOT')
      END IF
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      ALLOCATE(AEPHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUP$GETR8A('AECORE',NR,AECORE)
      CALL SETUP$UNSELECT()
      LMNX=SUM(2*LOX(:LNX)+1)
      ALLOCATE(DATH1(LMNX,LMNX,NDIMD))
      CALL LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,AEPHI,DATH1)
      DATH(:LMNX,:LMNX,:,IAT)=DATH(:LMNX,:LMNX,:,IAT)+DATH1(:,:,:)
      DEALLOCATE(AEPHI)
      DEALLOCATE(LOX)
      DEALLOCATE(DATH1)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &                          ,AEPOT,AEPHI,DATH)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE EXPECTATION VALUE OF                                 **
!     **  THE ONE-CENTER POTENTIALS WITH THE LOCAL ORBITALS                   **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
      COMPLEX(8),INTENT(OUT):: DATH(LMNX,LMNX,NDIMD)
      INTEGER(4)            :: LMN1,LMN2
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: LM1,LM2,LM3
      INTEGER(4)            :: L1,L2
      INTEGER(4)            :: IM1,IM2
      INTEGER(4)            :: ISPIN
      REAL(8)               :: AEDMU(NR,NDIMD)
      REAL(8)               :: DWORK1(NR)
      REAL(8)               :: CG
      REAL(8)               :: SVAR
      REAL(8)               :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      DATH(:,:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LMN2=0
          LM1=L1**2+IM1
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
!     
!             ==================================================================
!             ==  SUM ALL POTENTIALS THAT ACT ON THE GIVEN PAIR               ==
!             ==  OF PARTIAL WAVES                                            ==
!             ==================================================================
              AEDMU(:,:)=0.D0
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.NE.0.D0) THEN
                  DO ISPIN=1,NDIMD
                    AEDMU(:,ISPIN)=AEDMU(:,ISPIN)+CG*AEPOT(:,LM3,ISPIN)
                  ENDDO
                END IF
              ENDDO
!     
!             ==================================================================
!             ==  PERFORM NOW THE INTEGRATION                                 ==
!             ==================================================================
              DO ISPIN=1,NDIMD
                DWORK1(:)=AEDMU(:,ISPIN)*AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,DWORK1,SVAR)
                DATH(LMN1,LMN2,ISPIN)=DATH(LMN1,LMN2,ISPIN) &
    &                                +CMPLX(SVAR,0.D0,KIND=8)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$ETOT(LMNXX_,NDIMD_,NAT_,DENMAT_,DH_)
!     **************************************************************************
!     **  DENMAT_ ON INPUT IS CALCULATED DIRECTLY FROM THE PROJECTIONS AND    **
!     **  IS USED IN THE AUGMENTATION                                         **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : TON &
     &                       ,MODUS
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMNXX_
      INTEGER(4),INTENT(IN) :: NDIMD_
      INTEGER(4),INTENT(IN) :: NAT_
      COMPLEX(8),INTENT(IN) :: DENMAT_(LMNXX_,LMNXX_,NDIMD_,NAT_)
      COMPLEX(8),INTENT(OUT):: DH_(LMNXX_,LMNXX_,NDIMD_,NAT_)
!     **************************************************************************
      DH_=(0.D0,0.D0)
      IF(.NOT.TON) RETURN
                                    CALL TRACE$PUSH('LMTO$ETOT')
      WRITE(*,FMT='(82("="),T30," LMTO$ENERGY START. MODUS=",A," ")')TRIM(MODUS)
!
!     ==========================================================================
!     ==  PASS DENSITY MATRIX ON INTO LMTOAUGMENTATION OBJECT FOR USE IN THE  ==
!     ==  DOUBLE COUNTING TERM                                                ==
!     ==========================================================================
      CALL LMTOAUGMENTATION$SETI4('NAT',NAT_)
      CALL LMTOAUGMENTATION$SETI4('NDIMD',NDIMD_)
      CALL LMTOAUGMENTATION$SETI4('LMNXX',LMNXX_)
!     == SETC8A('DENMAT'...) ALLOCATES INTERNAL ARRAYS AN AND SETS TACTIVE=TRUE
      CALL LMTOAUGMENTATION$SETC8A('DENMAT',LMNXX_*LMNXX_*NDIMD_*NAT_,DENMAT_)
!
!     ==========================================================================
!     ==  SELECT CHOICES                                                      ==
!     ==========================================================================
      CALL LMTO$SETHTBCTOZERO()
      IF(MODUS.EQ.'DMFT') THEN
        CALL DMFT$GREEN()
      ELSE IF(MODUS.EQ.'HYBRID') THEN
        CALL LMTO_HYBRID()
      ELSE IF(MODUS.EQ.'ROBERT') THEN
        CALL LMTO_ROBERT()
      ELSE
        CALL ERROR$MSG('MODUS NOT RECOGNIZED')
        CALL ERROR$MSG('ALLOWED VALUES ARE "DMFT", "OLDDMFT", "HYBRID"')
        CALL ERROR$CHVAL('MODUS',MODUS)
        CALL ERROR$STOP('LMTO$ETOT')
      END IF
!
!     ==========================================================================
!     ==  COLLECT ONE-CENTER HAMILTONIAN FROM LMTOAUGMENTATION OBJECT         ==
!     ==========================================================================
      CALL LMTOAUGMENTATION$GETC8A('DH',LMNXX_*LMNXX_*NDIMD_*NAT_,DH_)
!
!     == DEALLOCATE INTERNAL ARRAYS AND SET OBJECT INACTIVE=====================
      CALL LMTOAUGMENTATION$CLEAN()
!
      WRITE(*,FMT='(82("="),T30," LMTO$ENERGY DONE ")')
                                    CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_ROBERT()
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : DENMAT=>DENMAT_NEW &
     &                       ,HAMIL=>HAMIL_NEW &
     &                       ,DEDOI &
     &                       ,INVOVERLAP &
     &                       ,POTPAR=>POTPAR1  &
     &                       ,ISPECIES &
     &                       ,HFWEIGHT &
     &                       ,HYBRIDSETTING
      IMPLICIT NONE
      TYPE NLIST_TYPE
        INTEGER(4) :: IAT
        INTEGER(4) :: IT(3)
        INTEGER(4) :: I1
      END TYPE NLIST_TYPE
      REAL(8)                :: ETOT,ETOT1
      INTEGER(4)             :: ISP ! ATOM TYPE
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: NND
      INTEGER(4)             :: NDIMD
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NATCL    ! #(ATOMS IN THE CLUSTER)
      INTEGER(4)             :: NORBCL     !#(ORBITALS W/O SPIN ON THE CLUSTER)
      REAL(8)   ,ALLOCATABLE :: UNS(:,:,:,:) ! NON-SPIN U-TENSOR
      REAL(8)   ,ALLOCATABLE :: U(:,:,:,:)   ! SPINOR U-TENSOR
      REAL(8)   ,ALLOCATABLE :: DEDU(:,:,:,:)! DERIVATIVE W.R.T.SPINOR U-TENSOR
      COMPLEX(8),ALLOCATABLE :: D(:,:)       ! SPINOR CLUSTER DENSITY MATRIX
      COMPLEX(8),ALLOCATABLE :: OINV(:,:)    ! SPINOR CLUSTER DENSITY MATRIX
      COMPLEX(8),ALLOCATABLE :: DEDOINV(:,:)    
      COMPLEX(8),ALLOCATABLE :: H(:,:)       ! SPINOR CLUSTER HAMILTONIAN
      REAL(8)                :: LHFWEIGHT
      TYPE(NLIST_TYPE),ALLOCATABLE :: NLIST(:)
      INTEGER(4)             :: I1UP,F1UP,I1DN,F1DN !INITIAL AND FINAL INDEX
      INTEGER(4)             :: I2UP,F2UP,I2DN,F2DN
      INTEGER(4)             :: N1,N2,N3
      INTEGER(4)             :: IND,IND1,IND2
      INTEGER(4)             :: IAT,IATA,IATB,IAT1,IAT2,NN,IT(3)
!     **************************************************************************
                               CALL TRACE$PUSH('LMTO_ROBERT')
      NAT=SIZE(ISPECIES)
!
!     ==========================================================================
!     == CALCULATE DENSITY MATRIX                                             ==
!     ==========================================================================
      CALL LMTO_NTBODENMAT_NEW()
      CALL LMTO_NTBOINVOVERLAP()!CALCULATE INVOVERLAP=SUM_N<PI|PSI_N><PSI_N|PI>
      NND=SIZE(DENMAT)
      NDIMD=DENMAT(1)%N3
!CALL LMTO$REPORTDENMAT(6)
CALL LMTO$REPORTPERIODICMAT(6,'INVERSE OVERLAP',NND,INVOVERLAP)
!STOP
!
!     ==========================================================================
!     == ALLOCATE HAMILTONIAN                                                 ==
!     ==========================================================================
      IF(ALLOCATED(HAMIL)) THEN
        CALL ERROR$MSG('HAMIL IS ALLOCATED')
        CALL ERROR$STOP('LMTO_ROBERT')
      END IF
      IF(ALLOCATED(DEDOI)) THEN
        CALL ERROR$MSG('DEDOI IS ALLOCATED')
        CALL ERROR$STOP('LMTO_ROBERT')
      END IF
      ALLOCATE(HAMIL(NND))
      ALLOCATE(DEDOI(NND))
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
!
        DEDOI(NN)%IAT1=DENMAT(NN)%IAT1
        DEDOI(NN)%IAT2=DENMAT(NN)%IAT2
        DEDOI(NN)%IT=DENMAT(NN)%IT
        N1=DENMAT(NN)%N1
        N2=DENMAT(NN)%N2
        N3=DENMAT(NN)%N3
        DEDOI(NN)%N1=N1
        DEDOI(NN)%N2=N2
        DEDOI(NN)%N3=N3
        ALLOCATE(DEDOI(NN)%MAT(N1,N2,N3))
        DEDOI(NN)%MAT=0.D0
      ENDDO

      ETOT=0.D0
      DO IAT=1,NAT      
        ISP=ISPECIES(IAT)
        IF(.NOT.(POTPAR(ISP)%NHEAD.GT.0)) CYCLE  ! NO CORRELATED ORBITALS
!
!       ========================================================================
!       == ADJUSTMENT IF LOCAL HFWEIGHT IS DIFFERENT FROM GLOBAL HFWEIGHT ======
!       ========================================================================
        IF(HFWEIGHT.GT.0.D0) THEN
          LHFWEIGHT=HYBRIDSETTING(ISP)%LHFWEIGHT
        ELSE
          LHFWEIGHT=HFWEIGHT
        END IF
!
!       ========================================================================
!       ==  DETERMINE U-TENSOR                                                ==
!       ========================================================================
        LMNX=SUM(2*POTPAR(ISP)%LOFH+1)
        ALLOCATE(UNS(LMNX,LMNX,LMNX,LMNX))
        CALL DMFT_ULOCAL(IAT,LMNX,UNS)
        ALLOCATE(U(2*LMNX,2*LMNX,2*LMNX,2*LMNX))
        U=0.D0
        U(  :LMNX,  :LMNX,  :LMNX,  :LMNX)=UNS
        U(LMNX+1:,LMNX+1:,LMNX+1:,LMNX+1:)=UNS
        U(LMNX+1:,  :LMNX,LMNX+1:,  :LMNX)=UNS
        U(  :LMNX,LMNX+1:,  :LMNX,LMNX+1:)=UNS
        DEALLOCATE(UNS)
!
!       ========================================================================
!       == DEFINE CLUSTER AND SET UP MAPPING                                  ==
!       ========================================================================
        NATCL=0
        DO NN=1,NND
          IF(DENMAT(NN)%IAT1.NE.IAT) CYCLE
          NATCL=NATCL+1
        ENDDO
!
        ALLOCATE(NLIST(NATCL))
        IND=0
        NORBCL=0
        DO NN=1,NND
          IF(DENMAT(NN)%IAT1.NE.IAT) CYCLE
          IND=IND+1
          NLIST(IND)%IAT=DENMAT(NN)%IAT2
          NLIST(IND)%IT=DENMAT(NN)%IT
          NLIST(IND)%I1=NORBCL+1       ! POSITION ON DENSITY MATRIX ARRAY
          NORBCL=NORBCL+DENMAT(NN)%N2
        ENDDO
!
!       ========================================================================
!       ==  EXTRACT DENSITY MATRIX ON THE CLUSTER                             ==
!       ========================================================================
!        ALLOCATE(D(2*LMNX,2*NORBCL))  ! THIS IS FOR STAR-LIKE CONSTRAINTS
        ALLOCATE(D(2*NORBCL,2*NORBCL))
        ALLOCATE(OINV(2*NORBCL,2*NORBCL))
        D(:,:)=(0.D0,0.D0)
        OINV(:,:)=(0.D0,0.D0)
!
!       == ADD THE TERMS CONNECTING THE IMPURITY WITH BATH SITES ==============
!       == THIS WILL BE OVERWRITTEN BY THE LOOP OVER BATH SITES. I LEAVE IT  ==
!       == IN FOR TESTING PURPOSES AND FOR LATER USE, IN CASE THE STAR-LIKE  ==
!       == CHOICE OF DENSITY MATRIX CONSTRAINTS TURNS OUT TO BE BETTER. =======
!       == THE STAR-LIKE CHOICE IS IMPLEMENTED MUCH MORE EFFICIENTLY ==========
!!$        IND1=0
!!$        DO NN=1,NND
!!$          IF(DENMAT(NN)%IAT1.NE.IAT) CYCLE
!!$!         = ASSUMES THAT CLUSTER ELEMENTS ARE TOGETHER
!!$          N1=DENMAT(NN)%N1
!!$          N2=DENMAT(NN)%N2
!!$          I1UP=1
!!$          I1DN=LMNX+1
!!$          I2UP=IND1+1
!!$          I2DN=NORBCL+IND1+1
!!$          CALL LMTO_ROBERT_MAP('FWRD',N1,N2,NDIMD,DENMAT(NN)%MAT &
!!$      &                       ,I1UP,I1DN,I2UP,I2DN,2*NORBCL,2*NORBCL,D)
!!$          CALL LMTO_ROBERT_MAP('FWRD',N1,N2,1,INVOVERLAP(NN)%MAT &
!!$      &                       ,I1UP,I1DN,I2UP,I2DN,2*NORBCL,2*NORBCL,OINV)
!!$          IND1=IND1+N2
!!$        ENDDO
!
!       == CONNECT BATH SITES WITH BATH SITES ==================================
        DO IATA=1,NATCL
          DO IATB=1,NATCL
            IAT1=NLIST(IATA)%IAT
            IAT2=NLIST(IATB)%IAT
            IT=NLIST(IATB)%IT-NLIST(IATA)%IT
            DO NN=1,NND
              IF(DENMAT(NN)%IAT1.NE.IAT1) CYCLE
              IF(DENMAT(NN)%IAT2.NE.IAT2) CYCLE
              IF(SUM(ABS(DENMAT(NN)%IT-IT)).NE.0) CYCLE
              N1=DENMAT(NN)%N1
              N2=DENMAT(NN)%N2
              I1UP=NLIST(IATA)%I1
              I2UP=NLIST(IATB)%I1
              I1DN=I1UP+NORBCL
              I2DN=I2UP+NORBCL
              CALL LMTO_ROBERT_MAP('FWRD',N1,N2,NDIMD,DENMAT(NN)%MAT &
      &                       ,I1UP,I1DN,I2UP,I2DN,2*NORBCL,2*NORBCL,D)
              CALL LMTO_ROBERT_MAP('FWRD',N1,N2,1,INVOVERLAP(NN)%MAT &
      &                       ,I1UP,I1DN,I2UP,I2DN,2*NORBCL,2*NORBCL,OINV)
              EXIT
            ENDDO ! NN
          ENDDO   ! IATB
        ENDDO     ! IATA
!
!       ========================================================================
!       ==  DETERMINE DENSITY MATRIX FUNCTIONAL                               ==
!       ========================================================================
!        ALLOCATE(H(2*LMNX,2*NORBCL))
        ALLOCATE(H(2*NORBCL,2*NORBCL))
        ALLOCATE(DEDOINV(2*NORBCL,2*NORBCL))
        ALLOCATE(DEDU(2*LMNX,2*LMNX,2*LMNX,2*LMNX))
        CALL LMTO_CLUSTERRDMFT(2*LMNX,2*NORBCL,U,D,OINV,ETOT1,H,DEDOINV,DEDU)
        ETOT=ETOT+ETOT1
        DEALLOCATE(D)
        DEALLOCATE(U)
!
!       ========================================================================
!       ==  TRANSFORM HAMILTONIAN
!       ========================================================================
!
!       == ADD THE TERMS CONNECTING THE IMPURITY WITH BATH SITES ==============
!       == THIS WILL BE OVERWRITTEN BY THE LOOP OVER BATH SITES. I LEAVE IT  ==
!       == IN FOR TESTING PURPOSES AND FOR LATER USE, IN CASE THE STAR-LIKE  ==
!       == CHOICE OF DENSITY MATRIX CONSTRAINTS TURNS OUT TO BE BETTER. =======
!       == THE STAR-LIKE CHOICE IS IMPLEMENTED MUCH MORE EFFICIENTLY ==========
!!$        IND1=0
!!$        DO NN=1,NND
!!$          IF(HAMIL(NN)%IAT1.NE.IAT) CYCLE
!!$!         = ASSUMES THAT CLUSTER ELEMENTS ARE TOGETHER
!!$          N1=HAMIL(NN)%N1
!!$          N2=HAMIL(NN)%N2
!!$          I1UP=1
!!$          I1DN=LMNX+1
!!$          I2UP=IND1+1
!!$          I2DN=NORBCL+IND1+1
!!$          CALL LMTO_ROBERT_MAP('BACK',N1,N2,NDIMD,HAMIL(NN)%MAT &
!!$      &                       ,I1UP,I1DN,I2UP,I2DN,2*NORBCL,2*NORBCL,H)
!!$          IND1=IND1+N2
!!$        ENDDO
!
!       == CONNECT BATH SITES WITH BATH SITES ==================================
        DO IATA=1,NATCL
          DO IATB=1,NATCL
            IAT1=NLIST(IATA)%IAT
            IAT2=NLIST(IATB)%IAT
            IT=NLIST(IATB)%IT-NLIST(IATA)%IT
            DO NN=1,NND
              IF(HAMIL(NN)%IAT1.NE.IAT1) CYCLE
              IF(HAMIL(NN)%IAT2.NE.IAT2) CYCLE
              IF(SUM(ABS(HAMIL(NN)%IT-IT)).NE.0) CYCLE
              N1=HAMIL(NN)%N1
              N2=HAMIL(NN)%N2
              I1UP=NLIST(IATA)%I1
              I2UP=NLIST(IATB)%I1
              I1DN=I1UP+NORBCL
              I2DN=I2UP+NORBCL
              CALL LMTO_ROBERT_MAP('BACK',N1,N2,NDIMD,HAMIL(NN)%MAT &
      &                       ,I1UP,I1DN,I2UP,I2DN,2*NORBCL,2*NORBCL,H)
              CALL LMTO_ROBERT_MAP('BACK',N1,N2,NDIMD,DEDOI(NN)%MAT &
      &                       ,I1UP,I1DN,I2UP,I2DN,2*NORBCL,2*NORBCL,DEDOINV)
              EXIT
            ENDDO ! NN
          ENDDO   ! IATB
        ENDDO     ! IATA
!
!       ========================================================================
!       ==  TRANSFORM DERIVATIVE OF U-TENSOR                                  ==
!       ========================================================================
!       == ATTENTION: THIS TERM IS NOT IMPLEMENTED!!!!!        
!
!       ========================================================================
!       ==  CLOSE DOWN                                                        ==
!       ========================================================================
        DEALLOCATE(NLIST)
        DEALLOCATE(H)
        DEALLOCATE(OINV)
        DEALLOCATE(DEDOINV)
        DEALLOCATE(DEDU)
      ENDDO
!
!     ==========================================================================
!     ==  MAP HAMILTONIAN CONTRIBUTION ON THIS%HTBC OF WAVES OBJECT           ==
!     ==========================================================================
      CALL LMTO_NTBODENMATDER_NEW()
!      CALL LMTO_NTBOINVOVERLAPDER_()
      CALL LMTO_CLEANDENMAT_NEW()
!
!     ==========================================================================
!     ==  REPORT ENERGY TO ENERGY LIST                                        ==
!     ==========================================================================
PRINT*,'ENERGY FROM LMTO INTERFACE ',ETOT
      CALL ENERGYLIST$SET('LMTO INTERFACE',ETOT)
      CALL ENERGYLIST$ADD('LOCAL CORRELATION',ETOT)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',ETOT)
                               CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_ROBERT_MAP(ID,N1,N2,NDIMD,MAT &
     &                          ,I1UP,I1DN,I2UP,I2DN,NB1,NB2,CLMAT)
!     **************************************************************************
!     **  ID='FWRD': MAPS A BLOCK OF THE DENSITY MATRIX INTO THE MATRIX       **
!     **             OF THE ENTIRE CLUSTER.                                   **
!     **  ID='BACK': MAPS A BLOCK OF THE HAMILTONIAN MATRIX INTO THE BLOCK    **
!     **             ON THE NEIGHBOR LIST                                     **
!     **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4),INTENT(IN) :: N1  ! SIZE OF 1ST DIMENSION OF THE BLOCK
      INTEGER(4),INTENT(IN) :: N2  ! SIZE OF 2ND DIMENSION OF THE BLOECK
      INTEGER(4),INTENT(IN) :: NDIMD !1,2,4 FOR NON-SPIN, COLLINEAR AND NONCOLL/
      REAL(8)   ,INTENT(INOUT) :: MAT(N1,N2,NDIMD) ! DENMAT OR H ON NGHBORLIST
      INTEGER(4),INTENT(IN) :: I1UP  !INSERTION POINT FOR 1ST INDEX UP-SPIN
      INTEGER(4),INTENT(IN) :: I1DN  !INSERTION POINT FOR 1ST INDEX DOWN-SPIN
      INTEGER(4),INTENT(IN) :: I2UP  !INSERTION POINT FOR 2ND INDEX UP-SPIN
      INTEGER(4),INTENT(IN) :: I2DN  !INSERTION POINT FOR 2ND INDEX DOWN-SPIN
      INTEGER(4),INTENT(IN) :: NB1
      INTEGER(4),INTENT(IN) :: NB2
      COMPLEX(8),INTENT(INOUT):: CLMAT(NB1,NB2) ! DENMAT OR H ON THE CLUSTER
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      INTEGER(4)            :: F1UP,F2UP,F1DN,F2DN,I,J
!     **************************************************************************
      F1UP=I1UP-1+N1
      F1DN=I1DN-1+N1
      F2UP=I2UP-1+N2
      F2DN=I2DN-1+N2
      IF(F1UP.GT.NB1.OR.F1DN.GT.NB1.OR.F2UP.GT.NB2.OR.F2DN.GT.NB2) THEN
        CALL ERROR$MSG('INCONSISTEN ARRAY SIZES')
        CALL ERROR$I4VAL('N1',N1)
        CALL ERROR$I4VAL('N2',N2)
        CALL ERROR$I4VAL('NB1',NB1)
        CALL ERROR$I4VAL('NB2',NB2)
        CALL ERROR$I4VAL('I1UP',I1UP)
        CALL ERROR$I4VAL('I2UP',I2UP)
        CALL ERROR$I4VAL('I1DN',I1DN)
        CALL ERROR$I4VAL('I2DN',I2DN)
        CALL ERROR$I4VAL('F1UP',F1UP)
        CALL ERROR$I4VAL('F2UP',F2UP)
        CALL ERROR$I4VAL('F1DN',F1DN)
        CALL ERROR$I4VAL('F2DN',F2DN)
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO_ROBERT_MAP')
      END IF
!
!     ==========================================================================
!     == FORWARD TRANSFORM OF THE DENSITY MATRIX                              ==
!     ==========================================================================
      IF(ID.EQ.'FWRD') THEN
        IF(NDIMD.EQ.1) THEN
          CLMAT(I1UP:F1UP,I2UP:F2UP)=0.5D0*MAT(:,:,1)
          CLMAT(I1DN:F1DN,I2DN:F2DN)=0.5D0*MAT(:,:,1)
        ELSE IF(NDIMD.EQ.2) THEN
          CLMAT(I1UP:F1UP,I2UP:F2UP)=0.5D0*( MAT(:,:,1)+MAT(:,:,2) )
          CLMAT(I1DN:F1DN,I2DN:F2DN)=0.5D0*( MAT(:,:,1)-MAT(:,:,2) )
        ELSE IF(NDIMD.EQ.4) THEN
          CLMAT(I1UP:F1UP,I2UP:F2UP)=0.5D0*( MAT(:,:,1)+   MAT(:,:,4) )
          CLMAT(I1UP:F1UP,I2UP:F2UP)=0.5D0*( MAT(:,:,1)+   MAT(:,:,4) )
          CLMAT(I1UP:F1UP,I2DN:F2DN)=0.5D0*( MAT(:,:,2)-CI*MAT(:,:,3) )
          CLMAT(I1DN:F1DN,I2UP:F2UP)=0.5D0*( MAT(:,:,2)+CI*MAT(:,:,3) )
        ELSE
          CALL ERROR$MSG('INVALID VALUE OF NDIMD')
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO_ROBERT_MAP')
        END IF
!
!     ==========================================================================
!     == BACK TRANSFORM OF THE HAMULTONIAN                                    ==
!     ==========================================================================
      ELSE IF(ID.EQ.'BACK') THEN
        IF(NDIMD.EQ.1) THEN
          DO J=1,N2
            DO I=1,N1
              MAT(I,J,1)=MAT(I,J,1)+0.5D0*REAL(CLMAT(I1UP-1+I,I2UP-1+J) &
       &                                      +CLMAT(I1DN-1+I,I2DN-1+J))
            ENDDO
          ENDDO
        ELSE IF(NDIMD.EQ.2) THEN
          DO J=1,N2
            DO I=1,N1
              MAT(I,J,1)=MAT(I,J,1)+0.5D0*REAL( CLMAT(I1UP-1+I,I2UP-1+J) &
       &                                       +CLMAT(I1DN-1+I,I2DN-1+J))
              MAT(I,J,2)=MAT(I,J,2)+0.5D0*REAL( CLMAT(I1UP-1+I,I2UP-1+J) &
       &                                       -CLMAT(I1DN-1+I,I2DN-1+J))
            ENDDO
          ENDDO
        ELSE IF(NDIMD.EQ.4) THEN
          DO J=1,N2
            DO I=1,N1
              MAT(I,J,1)=MAT(I,J,1)+0.5D0*REAL( CLMAT(I1UP-1+I,I2UP-1+J) &
       &                                       +CLMAT(I1DN-1+I,I2DN-1+J))
              MAT(I,J,4)=MAT(I,J,4)+0.5D0*REAL( CLMAT(I1UP-1+I,I2UP-1+J) &
       &                                       -CLMAT(I1DN-1+I,I2DN-1+J))
              MAT(I,J,2)=MAT(I,J,2)+0.5D0*REAL( CLMAT(I1UP-1+I,I2DN-1+J) &
       &                                       +CLMAT(I1DN-1+I,I2UP-1+J))
              MAT(I,J,3)=MAT(I,J,3)+0.5D0*AIMAG(CLMAT(I1UP-1+I,I2DN-1+J) &
       &                                       -CLMAT(I1DN-1+I,I2UP-1+J))
            ENDDO
          ENDDO
        ELSE
          CALL ERROR$MSG('INVALID VALUE OF NDIMD')
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LMTO_ROBERT_MAP')
        END IF 
      ELSE
        CALL ERROR$MSG('INVALID VALUE OF ID')
        CALL ERROR$MSG('MUST BE "FWRD" OR "BACK" (UPPERCASE)')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO_ROBERT_MAP')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_CLUSTERRDMFT(N1,N2,U,D,OINV,E,H,DEDOINV,DEDU)
!     **************************************************************************
!     **  PROVIDES THE NON-HARTREE FOCK CONTRIBUTION OF THE 1-PARTICLE REDUCED**
!     **  DENSITY MATRIX FUNCTIONAL FOR AN IMPURITY IN A CLUSTER              **
!     **                                                                      **
!     **  THE DENSITY-MATRIX ELEMENTS FROM THE IMPURITY SITE TO ALL ORBITALS  **
!     **  ON THE CLUSTER ARE USED. THE ORBITALS ARE ORDERED SUCH THAT ALL THE **
!     **  SPIN-UP ORBITALS ARE ARRANGED FIRST AND FOLLOWED BY THE SPIN-DOWN   **
!     **  ORBITALS                                                            **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: N1               ! #(IMPURITY ORBITALS)
      INTEGER(4),INTENT(IN)  :: N2               ! #(CLUSTER ORBITALS)
      REAL(8)   ,INTENT(IN)  :: U(N1,N1,N1,N1)   ! U-TENSOR
      COMPLEX(8),INTENT(IN)  :: D(N2,N2)         ! DENSITY MATRIX
      COMPLEX(8),INTENT(IN)  :: OINV(N2,N2)      ! INVERSE OVERLAP
      REAL(8)   ,INTENT(OUT) :: E                ! ENERGY
      COMPLEX(8),INTENT(OUT) :: H(N2,N2)         ! DEDD
      COMPLEX(8),INTENT(OUT) :: DEDOINV(N2,N2)   ! DEDOINV
      REAL(8)   ,INTENT(OUT) :: DEDU(N1,N1,N1,N1)! DEDU
      INTEGER(4) :: I,J
      COMPLEX(8) :: CSVAR
      REAL(8)    :: EW(N2)
!     **************************************************************************
      E=0.D0
      H=(0.D0,0.D0)
      DEDOINV=(0.D0,0.D0)
      DEDU=0.D0

PRINT*,'CHECK IF HERMITEAN...'
DO I=1,N2
  DO J=I,N2
    CSVAR=OINV(I,J)-CONJG(OINV(J,I)) 
    IF(ABS(CSVAR).GT.1.D-12) PRINT*,'I,J',I,J,OINV(I,J),OINV(J,I) 
  ENDDO
ENDDO    
DO I=1,N2
  DO J=I,N2
    CSVAR=D(I,J)-CONJG(D(J,I)) 
    IF(ABS(CSVAR).GT.1.D-12) PRINT*,'I,J',I,J,D(I,J),D(J,I) 
  ENDDO
ENDDO    
PRINT*,'...CHECK COMPLETED.'
DO I=1,N2/2
PRINT*,'OINV ',I,OINV(:N2/2,I)
ENDDO

      !EIGENVALUES OF OINV (POSITIVE DEFINITE?)
      CALL LIB$DIAGC8(N2,OINV,EW,H)
      DO J=1,N2
        PRINT*,"EW(OINV)",J,EW(J)
      ENDDO
PRINT*,'H1 ',H(:,1)
PRINT*,'H2 ',H(:,2)
STOP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_HYBRID()
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : TOFFSITE
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TTEST1=.FALSE.
      LOGICAL(4),PARAMETER  :: TTEST2=.FALSE.
      INTEGER(4)            :: SWITCH
!     **************************************************************************
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
      CALL LMTO_NTBODENMAT_NEW()
      CALL TIMING$CLOCKOFF('NTBODENMAT')
IF(TTEST1) THEN
  PRINT*,'MARKE BEFORE LMTO_LOCNATORB'
  CALL LMTO_LOCNATORB()
  PRINT*,'MARKE AFTER LMTO_LOCNATORB'
!  STOP
END IF



!!$      CALL LMTO_TESTDENMAT_1CDENMAT(LMNXX_,NDIMD_,NAT_,DENMAT_)
        PRINT*,'TOFFSITE',TOFFSITE
        IF(TOFFSITE) CALL LMTO_TESTDENMAT()
!STOP 'FORCED'

!!$      CALL LMTO$REPORTSBAR(6)
!!$      CALL LMTO$REPORTOVERLAP(6)
!STOP 'FORCED'
!
!     ==========================================================================
!     ==  SOME INFO                                                           ==
!     ==========================================================================
!!$      CALL LMTO$REPORTDENMAT(6)
!!$STOP 'FORCED'
!
!     ==========================================================================
!     ==  CALCULATE ENERGY                                                    ==
!     ==========================================================================
      IF(TTEST2) CALL LMTO_TESTHYBRIDENERGY() !STOPS AFTER TEST
!
      CALL TIMING$CLOCKON('NTBOETOT')
      CALL LMTO_HYBRIDENERGY()
      CALL TIMING$CLOCKOFF('NTBOETOT')
!
!     ==========================================================================
!     ==  CONVERT HAMIL INTO HTBC
!     ==========================================================================
      CALL TIMING$CLOCKON('NTBODENMATDER')
      CALL LMTO_NTBODENMATDER_NEW()
      CALL TIMING$CLOCKOFF('NTBODENMATDER')
!
!!$IF(TOFFSITE)CALL LMTO$REPORTOVERLAP(6)
!!$CALL LMTO$REPORTDENMAT(6)
!!$CALL LMTO$REPORTHAMIL(6)
!STOP 'FORCED'
!
!     ==========================================================================
!     ==  CLEAN DENMAT AND HAMIL                                              ==
!     ==========================================================================
      CALL LMTO_CLEANDENMAT_NEW()
      RETURN
      END

!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE  LMTO_TESTHYBRIDENERGY()
!     **************************************************************************
!     ** TEST ROUTINE FOR LMTO_HYBRIDENERGY                                   **
!     ** LMTO HYBRIDENERGY CALCULATES THE HARTREE-FOCK CORRECTION OF THE      **
!     ** HYBRID FUNCTIONALS IN THE LOCAL BASIS.                               **
!     ** LMTO_HYBRIDENERGY SUMS OVER ALL ATOMS.                               **
!     **                                                                      **
!     ** THE MATRIX ELEMENTS TO BE TESTED NEED TO BE SELECTED BY HAND         **
!     **************************************************************************
      USE LMTO_MODULE             ,ONLY: DENMAT=>DENMAT_NEW &
     &                                  ,HAMIL=>HAMIL_NEW
      USE WAVES_MODULE            ,ONLY: NSPIN,NDIM
!     __ DENMAT,DATH ARE COMPLEX(8),DIMENSION (LMNXX,LMNXX,NDIMD,NAT)
      USE LMTOAUGMENTATION_MODULE, ONLY: DENMAT_PARTWAVE=>DENMAT &
     &                                  ,DATH_PARTWAVE=>DATH
      IMPLICIT NONE
      INTEGER(4),PARAMETER  :: NIX=3
      REAL(8)   ,PARAMETER  :: XDELTA=1.D-2
      REAL(8)   ,PARAMETER  :: YDELTA=1.D-2
      INTEGER(4)            :: LMNXX,NDIMD,NAT
      INTEGER(4)            :: LMN1,LMN2,IDIMD,IAT
      INTEGER(4)            :: NN,IND1,IND2,IND3
      REAL(8)               :: XSVAR,YSVAR
      COMPLEX(8)            :: YCSVAR,YCDELTA
      REAL(8)               :: XVAL(-NIX:NIX)
      REAL(8)               :: XDER(-NIX:NIX)
      REAL(8)               :: YDER(-NIX:NIX)
      INTEGER(4)            :: IX
!     **************************************************************************
      CALL LMTOAUGMENTATION$GETI4('NAT',NAT)
      CALL LMTOAUGMENTATION$GETI4('LMNXX',LMNXX)
      CALL LMTOAUGMENTATION$GETI4('NDIMD',NDIMD)
!
!     ==========================================================================
!     == SELECT MATRIX ELEMENTS TO BE SCANNED                                 ==
!     ==========================================================================
!     == LOCAL-ORBITAL DENSITY MATRIX
      NN=1
      IND1=2
      IND2=2
      IND3=1
!     == PARTIAL WAVE DENSITY MATRIX
      LMN1=1
      LMN2=2
      IDIMD=IND3
      IAT=1
!
!     ==========================================================================
!     == REPORT                                                               ==
!     ==========================================================================
      WRITE(*,FMT='(82("="))')
      WRITE(*,FMT='(82("="),T20," TEST LMTO_HYBRIDENERGY  ")')
      WRITE(*,FMT='(82("="))')
      WRITE(*,FMT='("LOCAL-ORBITAL DENSITY MATRIX")')
      WRITE(*,FMT='(49("."),":",ES10.2,T1,"DISPLACEMENT")')XDELTA
      WRITE(*,FMT='(49("."),":",I3,T1,"ATOM INDEX")')IAT
      WRITE(*,FMT='(49("."),":",2I3,T1,"LMN1,LMN2")')IND1,IND2
      WRITE(*,FMT='(49("."),":",I3,T1,"IDIMD")')IND3
      WRITE(*,FMT='("PARTIAL-WAVE DENSITY MATRIX")')
      WRITE(*,FMT='(49("."),":",ES10.2,T1,"DISPLACEMENT")')YDELTA
      WRITE(*,FMT='(49("."),":",I3,T1,"ATOM INDEX")')IAT
      WRITE(*,FMT='(49("."),":",2I3,T1,"LMN1,LMN2")')LMN1,LMN2
      WRITE(*,FMT='(49("."),":",I3,T1,"IDIMD")')IDIMD
!
!     ==========================================================================
!     == CHECK ARRAY SIZE                                                     ==
!     ==========================================================================
      IF(    LMN1 .GT.LMNXX &
     &   .OR.LMN2 .GT.LMNXX &
     &   .OR.IDIMD.GT.NDIMD &
     &   .OR.IAT  .GT.NAT) THEN
        CALL ERROR$MSG('ERROR STOP DURING TESTING')
        CALL ERROR$MSG('LMN1,LMN2,IDIMD, OR IAT OUT OF RANGE')
        CALL ERROR$I4VAL('LMN1',LMN1)
        CALL ERROR$I4VAL('MAX(LMN1)',LMNXX)
        CALL ERROR$I4VAL('LMN2',LMN2)
        CALL ERROR$I4VAL('MAX(LMN2)',LMNXX)
        CALL ERROR$I4VAL('IDIMD',IDIMD)
        CALL ERROR$I4VAL('MAX(IDIMD)',NDIMD)
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$I4VAL('MAX(IAT)',NAT)
        CALL ERROR$STOP('LMTO_TESTHYBRIDENERGY')
      END IF
      IF(NN.GT.SIZE(DENMAT)) THEN
        CALL ERROR$MSG('ERROR STOP DURING TESTING')
        CALL ERROR$MSG('NN OUT OF RANGE')
        CALL ERROR$STOP('LMTO_TESTHYBRIDENERGY')
      END IF
      IF(    IND1.GT.SIZE(DENMAT(NN)%MAT(:,1,1)) &
     &   .OR.IND2.GT.SIZE(DENMAT(NN)%MAT(1,:,1)) &
     &   .OR.IND3.GT.SIZE(DENMAT(NN)%MAT(1,1,:))) THEN
        CALL ERROR$MSG('ERROR STOP DURING TESTING (TTEST2)')
        CALL ERROR$MSG('IND1,IND2 OR IND3 OUT OF RANGE')
        CALL ERROR$I4VAL('IND1',IND1)
        CALL ERROR$I4VAL('MAX(IND1)',SIZE(DENMAT(NN)%MAT(:,1,1)))
        CALL ERROR$I4VAL('IND2',IND2)
        CALL ERROR$I4VAL('MAX(IND2)',SIZE(DENMAT(NN)%MAT(1,:,1)))
        CALL ERROR$I4VAL('IND3',IND3)
        CALL ERROR$I4VAL('MAX(IND3)',SIZE(DENMAT(NN)%MAT(1,1,:)))
        CALL ERROR$STOP('LMTO_TESTHYBRIDENERGY')
      END IF
!
!     ==========================================================================
!     == CALCULATE ENERGIES AND DERIBATIVES                                   ==
!     ==========================================================================
      XSVAR=DENMAT(NN)%MAT(IND1,IND2,IND3)
      YCSVAR=DENMAT_PARTWAVE(LMN1,LMN2,IDIMD,IAT)
      YSVAR=REAL(YCSVAR,KIND=8)
      YCDELTA=CMPLX(YDELTA,0.D0,KIND=8)
      DO IX=-NIX,NIX,1
        DENMAT(NN)%MAT(IND1,IND2,IND3)      = XSVAR+ XDELTA*REAL(IX,KIND=8)
        DENMAT_PARTWAVE(LMN1,LMN2,IDIMD,IAT)=YCSVAR+YCDELTA*REAL(IX,KIND=8)
        DATH_PARTWAVE(:,:,:,:)=(0.D0,0.D0)
!
        CALL LMTO_HYBRIDENERGY()
!
        CALL ENERGYLIST$GET('LMTO INTERFACE',XVAL(IX))
        XDER(IX)=HAMIL(NN)%MAT(IND1,IND2,IND3)
        YDER(IX)=REAL(DATH_PARTWAVE(LMN1,LMN2,IDIMD,IAT),KIND=8)
        DEALLOCATE(HAMIL) ! IS ALLOCATED IN LMTO_HYBRIDENERGY
      ENDDO
!
!     ==========================================================================
!     == REPORT RESULT                                                        ==
!     ==========================================================================
      DO IX=1,NIX
        XSVAR=(XVAL(IX)-XVAL(-IX))/(2.D0*REAL(IX,KIND=8))
        WRITE(*,FMT='("-->VDER ",3F30.20)')REAL(IX,8),XSVAR &
 &          ,0.5D0*(XDER(IX)+XDER(-IX))*XDELTA+0.5D0*(YDER(IX)+YDER(-IX))*YDELTA
      ENDDO
!
!     ==========================================================================
!     == STOPPING IS REQUIRED BECAUSE INTERNAL DATA HAVE BEEN CHANGED         ==
!     ==========================================================================
      CALL ERROR$MSG('FORCED STOP DURING TESTING (TTEST2)')
      CALL ERROR$STOP('LMTO_TESTHYBRIDENERGY')
      STOP
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
      SUBROUTINE LMTO_LOCNATORB()
!     **************************************************************************
!     ** CALCULATE LOCAL NATURAL ORBITALS AND THEIR OCCUPATIONS               **
!     ** AND PRINT THEM. THE OVERLAP MATRIX IS CALCULATED FROM TAILED ORBITALS**
!     **************************************************************************
      USE LMTO_MODULE, ONLY: ISPECIES &
     &                      ,DENMAT_NEW &
     &                      ,SBAR_NEW &
     &                      ,POTPAR1 &
     &                      ,HYBRIDSETTING &
     &                      ,HFWEIGHT &
     &                      ,TCTE
      IMPLICIT NONE
      INTEGER(4)              :: NAT
      INTEGER(4)              :: NNX
      INTEGER(4)              :: NND
      INTEGER(4)              :: NNS
      INTEGER(4)              :: LMNH
      INTEGER(4)              :: LMNT
      INTEGER(4)              :: LMNXT
      INTEGER(4)              :: IAT1,IAT2,IT(3)
      LOGICAL(4)              :: TONSITE
      INTEGER(4)              :: IAT,ISP,NN,I,J
      REAL(8)                 :: LHFWEIGHT
      REAL(8)                 :: EV
      REAL(8)   ,ALLOCATABLE  :: OVERLAP(:,:)
      REAL(8)   ,ALLOCATABLE  :: OINV(:,:)
      REAL(8)   ,ALLOCATABLE  :: F(:)
      REAL(8)   ,ALLOCATABLE  :: F1(:)
      REAL(8)   ,ALLOCATABLE  :: ORB(:,:)
      REAL(8)   ,ALLOCATABLE  :: UNS(:,:,:,:)
!     **************************************************************************
                                   CALL TRACE$PUSH('LMTO_LOCNATORB')
      CALL CONSTANTS$GET('EV',EV)
      NAT=SIZE(ISPECIES)
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LMNH=SUM(2*POTPAR1(ISP)%LOFH+1) ! #(LOCAL ORBITALS FOR THIS ATOM)
        IF(LMNH.EQ.0) CYCLE
!
!       == IDENTIFY ONSITE DENSITY MATRIX INDEX ================================
        NNX=SIZE(DENMAT_NEW)
        DO NN=1,NNX
          IAT1=DENMAT_NEW(NN)%IAT1
          IAT2=DENMAT_NEW(NN)%IAT2
          IT=DENMAT_NEW(NN)%IT
          TONSITE=(IAT1.EQ.IAT2).AND.(SUM(ABS(IT)).EQ.0)
          IF(.NOT.TONSITE) CYCLE
          IF(IAT1.NE.IAT) CYCLE
          NND=NN
          EXIT
        ENDDO
!
!       == IDENTIFY ONSITE STRUCTURE CONSTANT INDEX ============================
        NNX=SIZE(SBAR_NEW)
        DO NN=1,NNX
          IAT1=SBAR_NEW(NN)%IAT1
          IAT2=SBAR_NEW(NN)%IAT2
          IT=SBAR_NEW(NN)%IT
          TONSITE=(IAT1.EQ.IAT2).AND.SUM(ABS(IT)).EQ.0
          IF(.NOT.TONSITE) CYCLE
          IF(IAT1.NE.IAT) CYCLE
          NNS=NN
          EXIT
        ENDDO
!
!SBAR_NEW(NNS)%MAT=0.D0
        ALLOCATE(OVERLAP(LMNH,LMNH))
        IF(TCTE) THEN
          LMNXT=POTPAR1(ISP)%TAILED%LMNX
          CALL LMTO_EXPANDLOCALWITHCTE('BACK',IAT,1,LMNH,LMNXT &
       &                              ,OVERLAP,POTPAR1(ISP)%TAILED%OVERLAP)
        ELSE
          OVERLAP=POTPAR1(ISP)%TAILED%OVERLAP(:LMNH,:LMNH)                 &
     &         -MATMUL(POTPAR1(ISP)%TAILED%OVERLAP(:LMNH,LMNH+1:)          &
     &                ,TRANSPOSE(SBAR_NEW(NNS)%MAT))                       &
     &         -MATMUL(SBAR_NEW(NNS)%MAT                                   &
     &                ,POTPAR1(ISP)%TAILED%OVERLAP(LMNH+1:,:LMNH))         &
     &         +MATMUL(SBAR_NEW(NNS)%MAT                                   &
     &                ,MATMUL(POTPAR1(ISP)%TAILED%OVERLAP(LMNH+1:,LMNH+1:) &
     &                       ,TRANSPOSE(SBAR_NEW(NNS)%MAT)))
       END IF
!!$PRINT*,'POTPAR1/OVERLAP(HH) ',POTPAR1(ISP)%TAILED%OVERLAP(:LMNH,:LMNH)
!!$PRINT*,'POTPAR1/OVERLAP(HT) ',POTPAR1(ISP)%TAILED%OVERLAP(:LMNH,LMNH+1:)
!!$PRINT*,'POTPAR1/OVERLAP(TH) ',POTPAR1(ISP)%TAILED%OVERLAP(LMNH+1:,:LMNH)
!!$PRINT*,'POTPAR1/OVERLAP(TT) ',POTPAR1(ISP)%TAILED%OVERLAP(LMNH+1:,LMNH+1:)
        WRITE(*,FMT='(80("="),T10," ONSITE OVERLAP FOR ATOM ",I3,"  ")')IAT
        DO I=1,LMNH
          WRITE(*,FMT='(10F10.5)')OVERLAP(I,:)
        ENDDO
        WRITE(*,FMT='(80("="),T10,"  DENSITY MATRIX FOR ATOM ",I3,"  ")')IAT
        DO I=1,LMNH
          WRITE(*,FMT='(10F10.5)')DENMAT_NEW(NND)%MAT(I,:,1)
        ENDDO
!
!       ========================================================================
!       == DETERMINE ONSITE NATURAL ORBITALS ===================================
!       ========================================================================
        ALLOCATE(OINV(LMNH,LMNH))
        ALLOCATE(F(LMNH))
        ALLOCATE(ORB(LMNH,LMNH))
        CALL LIB$INVERTR8(LMNH,OVERLAP,OINV)
        CALL LIB$GENERALEIGENVALUER8(LMNH,DENMAT_NEW(NND)%MAT(:,:,1),OINV,F,ORB)
        ORB=MATMUL(OINV,ORB)
!
!       == REORDER SO THAT OCCUPATIONS DECREASE ================================
        ALLOCATE(F1(LMNH))
        OINV=ORB
        F1=F
        DO I=1,LMNH
          F(LMNH+1-I)=F1(I)
          ORB(:,LMNH+1-I)=OINV(:,I)
        ENDDO          
        DEALLOCATE(OINV)
        DEALLOCATE(F1)
!
        WRITE(*,FMT='(80("="),T10,"  OCCUPATIONS FOR ATOM ",I3,"  ")')IAT
        WRITE(*,FMT='(10F10.5)')F
        WRITE(*,FMT='(80("="),T10,"  NATURAL ORBITALS FOR ATOM ",I3,"  ")')IAT
        DO I=1,LMNH
          WRITE(*,FMT='(10F10.5)')ORB(:,I)
        ENDDO
        DEALLOCATE(F)
!
        DEALLOCATE(OVERLAP)
!
!       ========================================================================
!       == DETERMINE ONSITE U-TENSOR OF LOCAL NATURAL ORBITALS                ==
!       ========================================================================
        LHFWEIGHT=HYBRIDSETTING(ISP)%LHFWEIGHT
!       __LHFWEIGHT=-1 IF NOT SET. IN THAT CASE TAKE THE GLOBAL DEFAULT_________
        IF(LHFWEIGHT.LT.0.D0) LHFWEIGHT=HFWEIGHT
!        
        ALLOCATE(UNS(LMNH,LMNH,LMNH,LMNH))
        CALL DMFT_ULOCAL(IAT,LMNH,UNS)
        CALL LMTO_MAPUTONATORB(LMNH,ORB,UNS)
!
        WRITE(*,FMT='(80("="),T10," BARE U-PARAMETER IN EV FOR ATOM "' &
       &          //',I3,"  ")')IAT
        DO I=1,LMNH
          WRITE(*,FMT='(10F12.1)')(UNS(J,I,J,I)/EV,J=1,LMNH)
        ENDDO
        WRITE(*,FMT='(80("="),T10,"  SCREENED U-PARAMETER IN EV FOR ATOM "' &
       &          //',I3,"  ")')IAT
        DO I=1,LMNH
          WRITE(*,FMT='(10F12.1)')(UNS(J,I,J,I)*LHFWEIGHT/EV,J=1,LMNH)
        ENDDO
        WRITE(*,FMT='(80("="),T10,"  BARE J-PARAMETER IN EV FOR ATOM "' &
       &          //',I3,"  ")')IAT
        DO I=1,LMNH
          WRITE(*,FMT='(10F12.1)')(UNS(J,I,I,J)/EV,J=1,LMNH)
        ENDDO
        WRITE(*,FMT='(80("="),T10,"  SCREENED J-PARAMETER IN EV FOR ATOM "' &
       &          //',I3,"  ")')IAT
        DO I=1,LMNH
          WRITE(*,FMT='(10F12.1)')(UNS(J,I,I,J)*LHFWEIGHT/EV,J=1,LMNH)
        ENDDO
!
        DEALLOCATE(UNS)
        DEALLOCATE(ORB)
      ENDDO
!!$CALL ERROR$MSG('FORCED STOP')
!!$CALL ERROR$STOP('LMTO_LOCNATORB')
                                   CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_MAPUTONATORB(LMNX,ORB,U)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: LMNX
      REAL(8)   ,INTENT(IN)    :: ORB(LMNX,LMNX)
      REAL(8)   ,INTENT(INOUT) :: U(LMNX,LMNX,LMNX,LMNX)
      REAL(8)                  :: U1(LMNX,LMNX,LMNX,LMNX)
      INTEGER(4)               :: I,J
!     **************************************************************************
      U1(:,:,:,:)=0.D0
      DO I=1,LMNX
        DO J=1,LMNX
          U1(:,:,:,I)=U1(:,:,:,I)+U(:,:,:,J)*ORB(J,I)
        ENDDO
      ENDDO
!
      U(:,:,:,:)=0.D0
      DO I=1,LMNX
        DO J=1,LMNX
          U(:,:,I,:)=U(:,:,I,:)+U1(:,:,J,:)*ORB(J,I)
        ENDDO
      ENDDO
!
      U1(:,:,:,:)=0.D0
      DO I=1,LMNX
        DO J=1,LMNX
          U1(:,I,:,:)=U1(:,I,:,:)+U(:,J,:,:)*ORB(J,I)
        ENDDO
      ENDDO
!
      U(:,:,:,:)=0.D0
      DO I=1,LMNX
        DO J=1,LMNX
          U(I,:,:,:)=U(I,:,:,:)+U1(J,:,:,:)*ORB(J,I)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_HYBRIDENERGY()
!     **************************************************************************
!     **  WORK OUT THE ENERGY USING THE LOCAL APPROXIMATION                   **
!     **  TAILED PARTIAL WAVES                                                **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES &
     &                       ,DENMAT=>DENMAT_NEW &
     &                       ,SBAR=>SBAR_NEW &
     &                       ,HAMIL=>HAMIL_NEW,LNX,LOX &
     &                       ,POTPAR1 &
     &                       ,TOFFSITE,HYBRIDSETTING,HFWEIGHT &
     &                       ,TCTE
      USE MPE_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      LOGICAL(4),PARAMETER  :: TPR2=.FALSE.
      INTEGER(4)            :: NND
      INTEGER(4)            :: NNS
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IND,INH,INS
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
      REAL(8)               :: EH,EX,EXTOT,EHTOT,EHDC,EXDC,Q,EAT
      INTEGER(4)            :: NN,IAT,I,J,K,L,IS,IAT1,IAT2,ISP
      INTEGER(4)            :: LMN,LN,IM
      REAL(8)               :: QSPIN(4)
      REAL(8)               :: SVAR
      REAL(8)               :: LHFWEIGHT
      INTEGER(4)            :: IDFTTYPE
INTEGER(4)            :: IDIMD
CHARACTER(128) :: STRING
REAL(8)   ,ALLOCATABLE:: T(:,:),UNT(:,:),MYMAT(:,:)
      REAL(8)               :: HFSCALE
      LOGICAL(4)            :: TACTIVE ! DOES THIS ATOM CONTRIBUTE?
      LOGICAL(4),PARAMETER  :: TPARALLEL=.FALSE.
      INTEGER(4)            :: THISTASK,NTASKS
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_HYBRIDENERGY')
IF(TPR2)PRINT*,'============ LMTO_HYBRIDENERGY ============================='
      NAT=SIZE(ISPECIES)
      IF(ALLOCATED(HAMIL)) THEN
        CALL ERROR$MSG('HAMIL IS ALLOCATED')
        CALL ERROR$STOP('LMTO_HYBRIDENERGY')
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
!
!     ==========================================================================
!     == PARALLELIZATION MODEL                                                ==
!     == THE OUTCOME OF THIS ROUTINE IS                                       ==
!     == -- HAMIL (LMTO_MODULE; HAMIL_NEW)                                    ==
!     == -- EXTOT (PASSED TO ENERGLIST)                                       ==
!     == -- DATH (VIA MODULE LMTOAUGMENTATIUON)                               ==
!     ==                                                                      ==
!     == ONLY ONE OF THE TASKS IN THE MONOMER GROUP ADDS TO HAMIL. AT THE END ==
!     == OF THIS ROUTINE THE RESILT IS SUMMED OVER ALL TASKS.                 ==
!     ==                                                                      ==
!     == THE PARAMETER TPARALLEL DECIDES WHETHER THE PARALLELIZATION IS DONE  ==
!     == OR WETHER EVERYTHING IS CALCULATED ON EACH TASK                      ==
!     ==                                                                      ==
!     ==========================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!
!     ==========================================================================
!     == LOOP OVER ALL ATOMS                                                  ==
!     ==========================================================================
      EXTOT=0.D0
      EHTOT=0.D0
      DO IAT=1,NAT
        IF(TPARALLEL.AND.MOD(IAT-1,NTASKS).NE.THISTASK-1) CYCLE
!
        ISP=ISPECIES(IAT)
        TACTIVE=(POTPAR1(ISP)%NHEAD.GT.0)
        IF(.NOT.TACTIVE) CYCLE
!
!       == COLLECT INFORMATION FROM SETUPS OBJECT ==============================
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(AECORE(NR))
        CALL SETUP$GETR8A('AECORE',NR,AECORE)
        CALL SETUP$GETI4('LMRX',LMRX)
        CALL SETUP$UNSELECT()
        LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
!
!       == ADJUSTMENT IF LOCAL HFWEIGHT IS DIFFERENT FROM GLOBAL HFWEIGHT ======
        LHFWEIGHT=HYBRIDSETTING(ISP)%LHFWEIGHT
        IF(LHFWEIGHT.LT.0.D0) LHFWEIGHT=HFWEIGHT
        HFSCALE=1.D0
        IF(HFWEIGHT.GT.0.D0)HFSCALE=LHFWEIGHT/HFWEIGHT
IF(TPR2)PRINT*,'LHFW=',LHFWEIGHT,' GHFW=',HFWEIGHT
!
!       == FIND ELEMENTS FOR DENSITY MATRIX, HAMILTONIAN AND STRUCTURE CONSTANTS
        NND=SIZE(DENMAT)
        NNS=SIZE(SBAR)
        CALL LMTO_INDEXLOCAL2(IAT,NND,DENMAT,IND)
        CALL LMTO_INDEXLOCAL2(IAT,NND,HAMIL,INH)
        CALL LMTO_INDEXLOCAL1(IAT,NNS,SBAR,INS)
        LMNX=DENMAT(IND)%N1
        NDIMD=DENMAT(IND)%N3
        LNXT=POTPAR1(ISP)%TAILED%LNX
        LMNXT=POTPAR1(ISP)%TAILED%LMNX
        ALLOCATE(LOXT(LNXT))
        LOXT=POTPAR1(ISP)%TAILED%LOX
!
!       ========================================================================
!       == CALCULATE U-TENSOR                                                 ==
!       ========================================================================
        ALLOCATE(U(LMNXT,LMNXT,LMNXT,LMNXT))
        U=POTPAR1(ISP)%TAILED%U
!
!       ========================================================================
!       == 
!       ========================================================================
        ALLOCATE(D(LMNX,LMNX,NDIMD))
        ALLOCATE(H(LMNX,LMNX,NDIMD))
        D=DENMAT(IND)%MAT
        H(:,:,:)=0.D0
!
        ALLOCATE(DT(LMNXT,LMNXT,NDIMD))
        ALLOCATE(HT(LMNXT,LMNXT,NDIMD))
!       == BLOW UP DENSITY MATRIX INTO A BASIS OF HEAD AND TAIL STATES
        IF(TCTE) THEN
          CALL LMTO_EXPANDLOCALWITHCTE('FWRD',IAT,NDIMD,LMNX,LMNXT,D,DT)
        ELSE
          CALL LMTO_EXPANDLOCAL('FWRD',NDIMD,LMNX,LMNXT,SBAR(INS)%MAT,D,DT)
        END IF

IF(TPR2) THEN
  IF(TACTIVE) THEN
    DO I=1,NDIMD
      WRITE(*,FMT='(82("="),T30," DENSITY MATRIX D FOR IAT=",I5,"  IDIM=",I5," ")')IAT,I
      DO LMN=1,LMNX
        WRITE(*,FMT='(200F10.5)')D(LMN,:,I) 
      ENDDO
    ENDDO
    DO I=1,NDIMD
      WRITE(*,FMT='(82("="),T30," DENSITY MATRIX DT FOR IAT=",I5,"  IDIM=",I5," ")')IAT,I
      DO LMN=1,LMNXT
        WRITE(*,FMT='(200F10.5)')DT(LMN,:,I)
      ENDDO
    ENDDO
  END IF
END IF
!
!       ========================================================================
!       == WORK OUT TOTAL ENERGY AND HAMILTONIAN                              ==
!       ========================================================================
        EH=0.D0
        EX=0.D0
        EAT=0.D0
        QSPIN=0.D0
        HT(:,:,:)=0.D0
        DO I=1,LMNXT
          DO J=1,LMNXT
            QSPIN(:NDIMD)=QSPIN(:NDIMD) &
       &                 +POTPAR1(ISP)%TAILED%OVERLAP(I,J)*DT(J,I,:)
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
!               == AN ADDITIONAL FACTOR COMES FROM THE REPRESENTATION ==========
!               == INTO TOTAL AND SPIN
                SVAR=-0.25D0*U(I,J,K,L)
                EX=EX+SVAR*SUM(DT(K,J,:)*DT(L,I,:))
                HT(K,J,:)=HT(K,J,:)+SVAR*DT(L,I,:) 
                HT(L,I,:)=HT(L,I,:)+SVAR*DT(K,J,:) 
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(TCTE) THEN
          CALL LMTO_EXPANDLOCALWITHCTE('BACK',IAT,NDIMD,LMNX,LMNXT,H,HT)
        ELSE
          CALL LMTO_EXPANDLOCAL('BACK',NDIMD,LMNX,LMNXT,SBAR(INS)%MAT,H,HT)
        END IF
        HAMIL(INH)%MAT=HAMIL(INH)%MAT+H*HFSCALE
        EXTOT=EXTOT+EX*HFSCALE
IF(TPR2) THEN
  PRINT*,'TOTAL VALENCE CHARGE ON ATOM=..............',IAT,QSPIN(1)
  PRINT*,'TOTAL SPIN[HBAR/2] ON ATOM=................',IAT,QSPIN(2:NDIMD)
  PRINT*,'EXACT VALENCE EXCHANGE ENERGY FOR ATOM=....',IAT,EX
END IF
        EAT=EAT+EX
        EX=0.D0
!       ========================================================================
!       == ADD CORE-VALENCE EXCHANGE                                          ==
!       ========================================================================
        IF(HYBRIDSETTING(ISP)%TCV) THEN
IF(.TRUE.) THEN
          CALL LMTO_CVX_ACTONPHI(IAT,HFWEIGHT*HFSCALE,EX)
          H=0.D0
          IF(HFSCALE.NE.0.D0) THEN
            EX=EX/(HFSCALE*HFWEIGHT)
          ELSE
            EX=0.D0
          END IF
          EXTOT=EXTOT+EX*HFSCALE
          EAT=EAT+EX
!!$
!!$          CALL LMTO_EXPANDLOCALWITHCTE('FWRD',IAT,NDIMD,LMNX,LMNXT,D,DT)
!!$          CALL LMTO_CVX_NEW(ISP,LMNXT,EX,DT(:,:,1),HT(:,:,1))
!!$          HT(:,:,2:)=0.D0
!!$          CALL LMTO_EXPANDLOCALWITHCTE('BACK',IAT,NDIMD,LMNX,LMNXT,H,HT)
!!$PRINT*,'--TEST--',IAT,ISP,EX,SUM(H(:,:,1)*D(:,:,1)),SUM(H(:,:,1)*TRANSPOSE(D(:,:,1)))

ELSE
!         CALL TEST_LMTO_CVX_NEW(IAT,ISP,NDIMD,LMNX,LMNXT,D)
          IF(TCTE) THEN
            CALL LMTO_EXPANDLOCALWITHCTE('FWRD',IAT,NDIMD,LMNX,LMNXT,D,DT)
          ELSE
            CALL ERROR$STOP('LMTO_HYBRIDENERGY')
            CALL LMTO_EXPANDLOCAL('FWRD',1,LMNX,LMNXT,SBAR(INS)%MAT,D,DT)
          END IF
!  
          CALL LMTO_CVX_NEW(ISP,LMNXT,EX,DT(:,:,1),HT(:,:,1))
          HT(:,:,2:)=0.D0
!!$EX=0.D0
!!$HT=0.D0
PRINT*,'TCTE ',TCTE,SUM(DT*HT)-EX,EX,IAT,ISP,LMNX,LMNXT,NDIMD &
&       ,SUM(ABS(HT(:,:,1)-TRANSPOSE(HT(:,:,1)))) &
&       ,SUM(ABS(DT(:,:,1)-TRANSPOSE(DT(:,:,1))))
!
          H=0.D0
          IF(TCTE) THEN
            CALL LMTO_EXPANDLOCALWITHCTE('BACK',IAT,NDIMD,LMNX,LMNXT,H,HT)
          ELSE
            CALL LMTO_EXPANDLOCAL('BACK',1,LMNX,LMNXT,SBAR(INS)%MAT,H,HT)
          END IF
          HAMIL(INH)%MAT(:,:,1)=HAMIL(INH)%MAT(:,:,1)+H(:,:,1)*HFSCALE
          EXTOT=EXTOT+EX*HFSCALE
          EAT=EAT+EX
END IF
IF(TPR2)PRINT*,'CORE-VALENCE EXCHANGE ENERGY FOR ATOM=.....',IAT,EX
          EX=0.D0
        END IF
!
!       ========================================================================
!       == DOUBLE COUNTING CORRECTION (EXCHANGE ONLY)                         ==
!       == THIS IS THE TIME CONSUMING PART                                    ==
!       ========================================================================
        CALL TIMING$CLOCKON('ENERGYTEST:DC')      
!!$IF(TACTIVE) THEN 
!!$PRINT*,'IAT=',IAT,NDIMD
!!$WRITE(*,FMT='(82("="),T10,"  H BEFORE DC ")')
!!$DO IDIMD=1,NDIMD
!!$  DO LMN=1,LMNX
!!$    WRITE(*,FMT='("IDIMD=",I1,":",100F10.5)')IDIMD,H(LMN,:,IDIMD)
!!$  ENDDO
!!$ENDDO
!!$END IF
!
IF(.TRUE.) THEN
        CALL DFT$SETL4('XCONLY',.TRUE.)
        CALL LMTOAUGMENTATION$SETI4('IAT',IAT)
        CALL LMTO_SIMPLEDC_NEW_NEW(GID,NR,NDIMD,LMNXT,LNXT,LOXT &
     &                    ,POTPAR1(ISP)%TAILED%AEF &
     &                    ,LRX,AECORE,DT,HFSCALE*HFWEIGHT,EX,HT)

!!$CALL TESTSIMPLEDC_NEW_NEW(GID,NR,NDIMD,LMNXT,LNXT,LOXT &
!!$     &                    ,POTPAR1(ISP)%TAILED%AEF &
!!$     &                    ,LRX,AECORE,DT,HFSCALE*HFWEIGHT,EX,HT)
!!$
        CALL LMTOAUGMENTATION$SETI4('IAT',0) !UNSET ATOM INDEX
        CALL DFT$SETL4('XCONLY',.FALSE.)
ELSE
        ALLOCATE(DTALL(LMNXT,LMNXT,NDIMD))
        ALLOCATE(HTALL(LMNXT,LMNXT,NDIMD))
        IF(TCTE) THEN
          CALL LMTO_EXPANDLOCALWITHCTE('FWRD',IAT,NDIMD,LMNX,LMNXT,D,DTALL)
        ELSE
          CALL LMTO_EXPANDLOCAL('FWRD',NDIMD,LMNX,LMNXT,SBAR(INS)%MAT,D,DTALL)
        END IF
        CALL LMTO_SIMPLEDC_NEW(GID,NR,NDIMD,LMNXT,LNXT,LOXT &
     &                    ,POTPAR1(ISP)%TAILED%AEF &
     &                    ,LRX,AECORE,DT,DTALL,EX,HT,HTALL)
        IF(TCTE) THEN
          CALL LMTO_EXPANDLOCALWITHCTE('BACK',IAT,NDIMD,LMNX,LMNXT,H,HTALL)
        ELSE
          CALL LMTO_EXPANDLOCAL('BACK',NDIMD,LMNX,LMNXT,SBAR(INS)%MAT,H,HTALL)
        END IF
        HAMIL(INH)%MAT=HAMIL(INH)%MAT-H*HFSCALE
        DEALLOCATE(DTALL)
        DEALLOCATE(HTALL)
END IF
        IF(TCTE) THEN
          CALL LMTO_EXPANDLOCALWITHCTE('BACK',IAT,NDIMD,LMNX,LMNXT,H,HT)
        ELSE
          CALL LMTO_EXPANDLOCAL('BACK',NDIMD,LMNX,LMNXT,SBAR(INS)%MAT,H,HT)
        END IF
        HAMIL(INH)%MAT=HAMIL(INH)%MAT-H*HFSCALE
        EXTOT=EXTOT-EX*HFSCALE  !HFWEIGT IS MULTIPLIED ON LATER
IF(TPR2) THEN
  IF(.NOT.TPARALLEL) THEN
    PRINT*,'DOUBLE COUNTING CORRECTION ENERGY FOR ATOM=',IAT,-EX
    PRINT*,'EXACT EXCHANGE ENERGY FOR ATOM........... =',IAT,EAT
    PRINT*,'EXCHANGE-XORRECTION FOR ATOM............. =',IAT,EAT-EX
  END IF
END IF
!
!!$IF(TACTIVE) THEN 
!!$PRINT*,'IAT=',IAT,NDIMD
!!$WRITE(*,FMT='(82("="),T10,"  H FROM SIMPLEDC ")')
!!$DO IDIMD=1,NDIMD
!!$  DO LMN=1,LMNX
!!$    WRITE(*,FMT='("IDIMD=",I1,":",100F10.5)')IDIMD,H(LMN,:,IDIMD)
!!$  ENDDO
!!$ENDDO
!!$END IF
!
        CALL TIMING$CLOCKOFF('ENERGYTEST:DC')      
        DEALLOCATE(HT)
        DEALLOCATE(DT)
        DEALLOCATE(U)
        DEALLOCATE(H)
        DEALLOCATE(D)
        DEALLOCATE(AECORE)
        DEALLOCATE(LOXT)
      ENDDO !END OF LOOP OVER ATOMS
!
!     ==========================================================================
!     == OFFSITE EXCHANGE CONTRIBUTION                                        ==
!     ==========================================================================
      IF(TOFFSITE) THEN
        CALL LMTO_OFFSITEXEVAL_NEW(TPARALLEL,EX) !USES NUMERICAL MATRIX ELEMENTS
        PRINT*,'+-+-+-+  OFFSITE EX=',EX
        EXTOT=EXTOT+EX
      END IF
!
!     ==========================================================================
!     == WRAP UP PARALLELIZATION                                              ==
!     ==========================================================================
      IF(TPARALLEL) THEN
        CALL MPE$COMBINE('MONOMER','+',EXTOT)
!       __ COMBINE DATH THE HAMILTONIAN IN TERMS OF PARTIAL WAVES_______________
        CALL LMTOAUGMENTATION$MPECOMBINEPOT()
        DO NN=1,NND
          CALL MPE$COMBINE('MONOMER','+',HAMIL(NN)%MAT)
        ENDDO
      END IF
!
!     ==========================================================================
!     == MAKE HAMILTONIAN HERMITIAN                                           ==
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
IF(TPR2)PRINT*,'ENERGY FROM LMTO INTERFACE ',EXTOT
      CALL ENERGYLIST$SET('LMTO INTERFACE',EXTOT)
      CALL ENERGYLIST$ADD('LOCAL CORRELATION',EXTOT)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',EXTOT)
!
!     ==========================================================================
!     == PRINT FOR TEST                                                       ==
!     ==========================================================================
      IF(TPR.AND.THISTASK.EQ.1) THEN
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
!!$PRINT*,'STOPPING AFTER PRINTING'
!!$STOP 'FORCED STOP'
      END IF
!
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_INDEXLOCAL2(IAT,NND,DENMAT,IND)
!     **************************************************************************
!     ** IDENTIFY THE ONSIDE ELEMENT FOR ATOM IAT ON THE DENSITY MATRIX       **
!     **************************************************************************
      USE LMTO_MODULE , ONLY : PERIODICMAT2_TYPE
      IMPLICIT NONE
      INTEGER(4)             ,INTENT(IN) :: IAT
      INTEGER(4)             ,INTENT(IN) :: NND
      TYPE(PERIODICMAT2_TYPE),INTENT(IN) :: DENMAT(NND)
      INTEGER(4)             ,INTENT(OUT):: IND
      INTEGER(4)                         :: NN
!     **************************************************************************
      IND=-1
      DO NN=1,NND
        IF(DENMAT(NN)%IAT1.NE.IAT) CYCLE
        IF(DENMAT(NN)%IAT2.NE.IAT) CYCLE
        IF(SUM(ABS(DENMAT(NN)%IT)).NE.0) CYCLE
        IND=NN
        EXIT
      ENDDO
      IF(IND.LE.0) THEN
        CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED: IND<0')
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$I4VAL('NND',NND)
        CALL ERROR$STOP('LMTO_INDEXLOCAL2')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_INDEXLOCAL1(IAT,NNS,SBAR,INS)
!     **************************************************************************
!     ** IDENTIFY THE ONSIDE ELEMENT FOR ATOM IAT ON THE DENSITY MATRIX       **
!     **************************************************************************
      USE LMTO_MODULE , ONLY : PERIODICMAT_TYPE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: IAT
      INTEGER(4)            ,INTENT(IN) :: NNS
      TYPE(PERIODICMAT_TYPE),INTENT(IN) :: SBAR(NNS)
      INTEGER(4)            ,INTENT(OUT):: INS
      INTEGER(4)                        :: NN
!     **************************************************************************
      INS=-1
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
        IF(SBAR(NN)%IAT2.NE.IAT) CYCLE
        IF(SUM(ABS(SBAR(NN)%IT)).NE.0) CYCLE
        INS=NN
        EXIT
      ENDDO
      IF(INS.LE.0) THEN
        CALL ERROR$MSG('INTERNAL CONSISTENCY CHECK FAILED: INS<0')
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$I4VAL('NNS',NNS)
        CALL ERROR$STOP('LMTO_INDEXLOCAL1')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILEDEXPANSION()
!     **************************************************************************
!     **  DETERMINES THE EXPANSION COEFFICIENT OF THE NTBOS IN ONE-CENTER    **
!     **  EXPANSION TERMS OF THE TAILED REPRESENTATION                        **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR=>POTPAR1 &
     &                       ,SBAR=>SBAR_NEW &
     &                       ,ISPECIES &
     &                       ,K2 &
     &                       ,CTE       !EXPANSION IN TAILED ORBITALS
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPRINT=.FALSE.
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NNS
      INTEGER(4)             :: LMNXH
      INTEGER(4)             :: LMNXT
      INTEGER(4)             :: IAT1,IAT2
      INTEGER(4)             :: ISP1,ISP2
      INTEGER(4)             :: IT(3)
      LOGICAL(4)             :: TONSITE
      INTEGER(4)             :: L1X,L2X
      REAL(8)                :: RBAS(3,3)
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      REAL(8)   ,ALLOCATABLE :: S(:,:)
      REAL(8)                :: CENTER(3)
      INTEGER(4)             :: NN,IAT,ISP
      INTEGER(4)             :: LMNH,LMNT,LMNA,LMNOFA,LMNOFT
      INTEGER(4)             :: NH,NT,NA
      INTEGER(4)             :: LH,LT,LA
      INTEGER(4)             :: MH,MT,MA
      INTEGER(4)             :: ISVAR
      INTEGER(4)             :: LN,L,M,LMN
      REAL(8)                :: SVAR
!     **************************************************************************
!
!     ==========================================================================
!     == COLLECT ATOMIC STRCTURE                                              ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     == ALLOCATE ARRAY FOR EXPANSION COEFFIENTS                              ==
!     ==========================================================================
      IF(.NOT.ASSOCIATED(CTE)) THEN
        ALLOCATE(CTE(NAT))
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          LMNXH=SUM(2*POTPAR(ISP)%LOFH(:)+1)
          LMNXT=SUM(2*POTPAR(ISP)%TAILED%LOX(:)+1)
          CTE(IAT)%LMNXT=LMNXT    
          CTE(IAT)%LMNXH=LMNXH
          ALLOCATE(CTE(IAT)%MAT(LMNXT,LMNXH))
        ENDDO
      END IF
!
!     ==========================================================================
!     == 
!     ==========================================================================
      DO IAT=1,NAT
        CTE(IAT)%MAT=0.D0
      ENDDO
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IAT1=SBAR(NN)%IAT1
        IAT2=SBAR(NN)%IAT2
        ISP1=ISPECIES(IAT1)
        ISP2=ISPECIES(IAT2)
        IT=SBAR(NN)%IT
        LMNXH=SUM(2*POTPAR(ISP1)%LOFH+1)
        TONSITE=(IAT1.EQ.IAT2).AND.(SUM(IT**2).EQ.0)
        IF(TONSITE) THEN
          DO LMNH=1,LMNXH
            CTE(IAT1)%MAT(LMNH,LMNH)=CTE(IAT1)%MAT(LMNH,LMNH)+1.D0
          ENDDO
          ISVAR=SBAR(NN)%N2
          CTE(IAT1)%MAT(LMNXH+1:LMNXH+ISVAR,:) &
   &               =CTE(IAT1)%MAT(LMNXH+1:LMNXH+ISVAR,:)-TRANSPOSE(SBAR(NN)%MAT)
        ELSE
!         == VARIABLES H REFER TO THE LOCAL ORBITALS AT THE CENTRAL SITE
!         == VARIABLES T REFER TO THE TAILS ON THE NEIGHBORING SITE
!         == VARIABLES A REFER TO THE TAILED COMPONENTS ON THE CENTRAL SITE
!         == |KINFTY_A>=-SUM_B |J_B> S_BA
!         == |KBAR_A>=|K_A> - ONSITESUM_B |JBAR_B> SBARDAGGER_BA
!         == FOR THE HIGHER ANGULAR MOMENTA USE
!         ==  SBARDAGGER_BA = [OFFSITESUM_C S_BC * QBAR_C * SBARDAGGER_CA]
          CENTER=R0(:,IAT1)-R0(:,IAT2) &
   &            -RBAS(:,1)*IT(1)-RBAS(:,2)*IT(2)-RBAS(:,3)*IT(3)
          L1X=MAXVAL(POTPAR(ISP2)%LOFT,DIM=1)
          L2X=MAXVAL(POTPAR(ISP1)%TAILED%LOX,DIM=1)
!         == CATCH ZERO-SIZED ARRAYS
          IF(SIZE(POTPAR(ISP2)%LOFT).EQ.0)      L1X=0
          IF(SIZE(POTPAR(ISP2)%TAILED%LOX).EQ.0)L2X=0
          ALLOCATE(S((L1X+1)**2,(L2X+1)**2))
!         == CENTER IS THE SITE OF THE BESSEL-FUNCTION EXPANSION RELATIVE TO ===
!         == THE POSITION OF THE HANKELFUNCTION ================================
          CALL LMTO$STRUCTURECONSTANTS(CENTER,K2,L1X,L2X,S)
          LMNH=0
          DO NH=1,POTPAR(ISP1)%NHEAD
            LH=POTPAR(ISP1)%LOFH(NH)
            DO MH=1,2*LH+1
              LMNH=LMNH+1
              LMNT=0
              DO NT=1,POTPAR(ISP2)%NTAIL
                LT=POTPAR(ISP2)%LOFT(NT)
                LMNOFT=LT**2
                DO MT=1,2*LT+1
                  LMNT=LMNT+1
                  LMNOFT=LMNOFT+1
                  LMNA=SUM(2*POTPAR(ISP1)%LOFH+1)+SUM(2*POTPAR(ISP1)%LOFT+1)
                  ISVAR=POTPAR(ISP1)%NHEAD+POTPAR(ISP1)%NTAIL+1
                  DO NA=ISVAR,POTPAR(ISP1)%TAILED%LNX
                    LA=POTPAR(ISP1)%TAILED%LOX(NA)
                    LMNOFA=LA**2
                    DO MA=1,2*LA+1 
                      LMNA=LMNA+1
                      LMNOFA=LMNOFA+1
                      CTE(IAT1)%MAT(LMNA,LMNH)=CTE(IAT1)%MAT(LMNA,LMNH) &
     &                                        -SBAR(NN)%MAT(LMNH,LMNT) &
     &                                         *POTPAR(ISP2)%QBAR(NT) &
     &                                         *S(LMNOFT,LMNOFA) 
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(S)
        END IF
      ENDDO
!
!     ==========================================================================
!     == 
!     ==========================================================================
      IF(TPRINT) THEN
        DO IAT=1,NAT
          LMNXT=CTE(IAT)%LMNXT
          LMNXH=CTE(IAT)%LMNXH
          WRITE(*,FMT='(82("-"),T5," EXPANSION OF NTBO IN TAILED COMPONENTS "' &
     &        //'" FOR ATOM ",I3," ")')IAT
!!$          DO LMNH=1,LMNXH
!!$            WRITE(*,FMT='("LMN=",I3," CTE=",100F10.5)')LMNH,CTE(IAT)%MAT(:,LMNH)
!!$          ENDDO
          LMN=0
          WRITE(*,FMT='(4A4,10A20)')'LN','L','M','LMN','CTE'
          ISP=ISPECIES(IAT)
          DO LN=1,POTPAR(ISP)%TAILED%LNX
            L=POTPAR(ISP)%TAILED%LOX(LN)
            DO M=1,2*L+1
              LMN=LMN+1
              SVAR=MAXVAL(ABS(CTE(IAT)%MAT(LMN,:)))
              IF(SVAR.LT.1.D-10) CYCLE
              IF(SVAR.LT.1.D+4) THEN
                WRITE(*,FMT='(4I4,100F12.6)')LN,L,M,LMN,CTE(IAT)%MAT(LMN,:)
              ELSE
                WRITE(*,FMT='(4I4,100E12.3)')LN,L,M,LMN,CTE(IAT)%MAT(LMN,:)
              END IF
            ENDDO
          ENDDO
        ENDDO
!!$        CALL ERROR$MSG('PLANNED STOP AFTER PRINTING')
!!$        CALL ERROR$STOP('LMTO_TAILEDEXPANSION')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_EXPANDLOCALWITHCTE(ID,IAT,NDIMD,LMNX,LMNXT,X,XT)
!     **************************************************************************
!     **  BLOW UP DENSITY MATRIX UP INTO A REPRESENTATION OF                  **
!     **  PURE ANGULAR MOMENTUM COMPONENTS IN THE ONE-CENTER TAILED EXPANSION **
!     **  AND SHRINK HAMILTOINIAN DOWN TO THE NTBO REPRESENTATION             **
!     **                                                                      **
!     ** ID='FWRD' X=INTENT(IN),XT=INTENT(OUT)                                **
!     ** ID='BACK' X=INTENT(OUT),XT=INTENT(IN)                                **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : CTE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID     ! CAN BE 'FWRD' AND 'BACK'
      INTEGER(4)  ,INTENT(IN)    :: IAT    ! ATOM INDEX
      INTEGER(4)  ,INTENT(IN)    :: NDIMD
      INTEGER(4)  ,INTENT(IN)    :: LMNX   ! DIMENSION NTBO-REPRESENTATION
      INTEGER(4)  ,INTENT(IN)    :: LMNXT  ! DIMENSION ONE-CENTER TERMS
      REAL(8)     ,INTENT(INOUT) :: X(LMNX,LMNX,NDIMD)
      REAL(8)     ,INTENT(INOUT) :: XT(LMNXT,LMNXT,NDIMD)
      INTEGER(4)                 :: IDIMD
!     **************************************************************************
      IF(ID.EQ.'FWRD') THEN
        DO IDIMD=1,NDIMD
          XT(:,:,IDIMD)=MATMUL(CTE(IAT)%MAT &
     &                        ,MATMUL(X(:,:,IDIMD),TRANSPOSE(CTE(IAT)%MAT)))
        ENDDO
!
      ELSE IF(ID.EQ.'BACK') THEN
        DO IDIMD=1,NDIMD
          X(:,:,IDIMD)=MATMUL(TRANSPOSE(CTE(IAT)%MAT) &
     &                       ,MATMUL(XT(:,:,IDIMD),CTE(IAT)%MAT))
        ENDDO
      ELSE
        CALL ERROR$STOP('LMTO_EXPANDLOCALWITHCTE')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_EXPANDNONLOCALWITHCTE(ID,IAT1,IAT2,NDIMD &
     &                                     ,LMNX1,LMNX2,LMNXT1,LMNXT2,X,XT)
!     **************************************************************************
!     **  BLOW UP DENSITY MATRIX UP INTO A REPRESENTATION OF                  **
!     **  PURE ANGULAR MOMENTUM COMPONENTS IN THE ONE-CENTER TAILED EXPANSION **
!     **  AND SHRINK HAMILTOINIAN DOWN TO THE NTBO REPRESENTATION             **
!     ** 
!     **************************************************************************
      USE LMTO_MODULE, ONLY : CTE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID     ! CAN BE 'FWRD' AND 'BACK'
      INTEGER(4)  ,INTENT(IN)    :: IAT1   ! ATOM CONNECTED TO 1ST INDEX
      INTEGER(4)  ,INTENT(IN)    :: IAT2   ! ATOM CONNECTED TO 2ND INDEX
      INTEGER(4)  ,INTENT(IN)    :: NDIMD
      INTEGER(4)  ,INTENT(IN)    :: LMNX1  ! 1ST DIMENSION NTBO-REPRESENTATION
      INTEGER(4)  ,INTENT(IN)    :: LMNX2  ! 2ND DIMENSION NTO-REPRESENTATION
      INTEGER(4)  ,INTENT(IN)    :: LMNXT1 ! 1ST DIMENSION ONE-CENTER TERMS
      INTEGER(4)  ,INTENT(IN)    :: LMNXT2 ! 2ND DIMENSION ONE-CENTER TERMS
      REAL(8)     ,INTENT(INOUT) :: X(LMNX1,LMNX2,NDIMD)
      REAL(8)     ,INTENT(INOUT) :: XT(LMNXT1,LMNXT2,NDIMD)
      INTEGER(4)                 :: IDIMD
!     **************************************************************************
      IF(ID.EQ.'FWRD') THEN
        DO IDIMD=1,NDIMD
          XT(:,:,IDIMD)=MATMUL(CTE(IAT1)%MAT &
     &                        ,MATMUL(X(:,:,IDIMD),TRANSPOSE(CTE(IAT2)%MAT)))
        ENDDO
!
      ELSE IF(ID.EQ.'BACK') THEN
        DO IDIMD=1,NDIMD
          X(:,:,IDIMD)=MATMUL(TRANSPOSE(CTE(IAT1)%MAT) &
     &                       ,MATMUL(XT(:,:,IDIMD),CTE(IAT2)%MAT))
        ENDDO
      ELSE
        CALL ERROR$STOP('LMTO_EXPANDNONLOCALWITHCTE')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_EXPANDLOCAL(ID,NDIMD,LMNX,LMNXT,SBAR,X,XT)
!     **************************************************************************
!     ** BLOW DENSITY MATRIX UP INTO A REPRESENTATION OF HEAD AND TAIL STATES **
!     ** AND SHRINK HAMILTOINIAN DOWN TO HEAD-ONLY REPRESENTATION             **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID     ! CAN BE 'FWRD' AND 'BACK'
      INTEGER(4)  ,INTENT(IN)    :: NDIMD
      INTEGER(4)  ,INTENT(IN)    :: LMNX
      INTEGER(4)  ,INTENT(IN)    :: LMNXT
      REAL(8)     ,INTENT(IN)    :: SBAR(LMNX,LMNXT-LMNX) !STRUCTURE CONSTANTS
      REAL(8)     ,INTENT(INOUT) :: X(LMNX,LMNX,NDIMD)
      REAL(8)     ,INTENT(INOUT) :: XT(LMNXT,LMNXT,NDIMD)
      INTEGER(4)                 :: IDIMD
!     **************************************************************************
      IF(ID.EQ.'FWRD') THEN
        DO IDIMD=1,NDIMD
          XT(:LMNX,:LMNX,IDIMD)    =X(:,:,IDIMD)
          XT(LMNX+1:,:LMNX,IDIMD)  =-MATMUL(TRANSPOSE(SBAR),X(:,:,IDIMD))
          XT(:LMNX,LMNX+1:,IDIMD)  =-MATMUL(X(:,:,IDIMD),SBAR)
          XT(LMNX+1:,LMNX+1:,IDIMD)=MATMUL(TRANSPOSE(SBAR) &
     &                                    ,MATMUL(X(:,:,IDIMD),SBAR))
        ENDDO
!
      ELSE IF(ID.EQ.'BACK') THEN
        DO IDIMD=1,NDIMD
          X(:,:,IDIMD)=XT(:LMNX,:LMNX,IDIMD)                           &
      &               -MATMUL(SBAR,XT(LMNX+1:,:LMNX,IDIMD))            &
      &               -MATMUL(XT(:LMNX,LMNX+1:,IDIMD),TRANSPOSE(SBAR)) &
      &               +MATMUL(SBAR,MATMUL(XT(LMNX+1:,LMNX+1:,IDIMD)    &
      &                                  ,TRANSPOSE(SBAR)))
        ENDDO
      ELSE
        CALL ERROR$STOP('LMTO_EXPANDLOCAL')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_EXPANDNONLOCAL(ID,NDIMD,LMNX1,LMNX2,LMNXT1,LMNXT2 &
     &                              ,SBAR1,SBAR2,X,XT)
!     **************************************************************************
!     ** BLOW DENSITY MATRIX UP INTO A REPRESENTATION OF HEAD AND TAIL STATES **
!     ** AND SHRINK HAMILTOINIAN DOWN TO HEAD-ONLY REPRESENTATION             **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID     ! CAN BE 'FWRD' AND 'BACK'
      INTEGER(4)  ,INTENT(IN)    :: NDIMD
      INTEGER(4)  ,INTENT(IN)    :: LMNX1
      INTEGER(4)  ,INTENT(IN)    :: LMNX2
      INTEGER(4)  ,INTENT(IN)    :: LMNXT1
      INTEGER(4)  ,INTENT(IN)    :: LMNXT2
      REAL(8)     ,INTENT(IN)    :: SBAR1(LMNX1,LMNXT1-LMNX1) !STRUCTURE CNSTNTS
      REAL(8)     ,INTENT(IN)    :: SBAR2(LMNX2,LMNXT2-LMNX2) !STRUCTURE CNSTNTS
      REAL(8)     ,INTENT(INOUT) :: X(LMNX1,LMNX2,NDIMD)
      REAL(8)     ,INTENT(INOUT) :: XT(LMNXT1,LMNXT2,NDIMD)
      INTEGER(4)                 :: IDIMD
!     **************************************************************************
      IF(ID.EQ.'FWRD') THEN
        DO IDIMD=1,NDIMD
          XT(:LMNX1,:LMNX2,IDIMD)    =X(:,:,IDIMD)
          XT(LMNX1+1:,:LMNX2,IDIMD)  =-MATMUL(TRANSPOSE(SBAR1),X(:,:,IDIMD))
          XT(:LMNX1,LMNX2+1:,IDIMD)  =-MATMUL(X(:,:,IDIMD),SBAR2)
          XT(LMNX1+1:,LMNX2+1:,IDIMD)=MATMUL(TRANSPOSE(SBAR1) &
     &                                    ,MATMUL(X(:,:,IDIMD),SBAR2))
        ENDDO
!
      ELSE IF(ID.EQ.'BACK') THEN
        DO IDIMD=1,NDIMD
          X(:,:,IDIMD)=XT(:LMNX1,:LMNX2,IDIMD)                           &
      &               -MATMUL(SBAR1,XT(LMNX1+1:,:LMNX2,IDIMD))            &
      &               -MATMUL(XT(:LMNX1,LMNX2+1:,IDIMD),TRANSPOSE(SBAR2)) &
      &               +MATMUL(SBAR1,MATMUL(XT(LMNX1+1:,LMNX2+1:,IDIMD)    &
      &                                  ,TRANSPOSE(SBAR2)))
        ENDDO
      ELSE
        CALL ERROR$STOP('LMTO_EXPANDNONLOCAL')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TEST_LMTO_CVX_NEW(IAT,ISP,NDIMD,LMNX,LMNXT,D_IN)
!     **************************************************************************
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: ISP
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: LMNX
      INTEGER(4),INTENT(IN) :: LMNXT
      REAL(8)   ,INTENT(IN) :: D_IN(LMNX,LMNX,NDIMD)
      REAL(8)               :: D(LMNX,LMNX,NDIMD)
      REAL(8)               :: DT(LMNXT,LMNXT,NDIMD)
      REAL(8)               :: HT(LMNXT,LMNXT,NDIMD)
      REAL(8)               :: H(LMNX,LMNX,NDIMD)
      REAL(8)               :: DELTA(LMNX,LMNX,NDIMD)
      INTEGER(4),PARAMETER  :: NDIS=4
      REAL(8)               :: EARR(-NDIS:NDIS)
      REAL(8)               :: DARR(-NDIS:NDIS)
      REAL(8)               :: EX
      REAL(8)               :: RAN
      REAL(8)               :: Y1,Y2
      INTEGER(4)            :: LMN1,LMN2,IDIS
!     **************************************************************************
      DELTA=0.D0
      DO LMN1=1,LMNX
        DO LMN2=LMN1,LMNX
          CALL RANDOM_NUMBER(RAN)
          DELTA(LMN1,LMN2,1)=-0.5D0+RAN
        ENDDO
      ENDDO
      DELTA(:,:,1)=0.5D0*(DELTA(:,:,1)+TRANSPOSE(DELTA(:,:,1)))
      DELTA=DELTA*1.D-3
!
      DO IDIS=-NDIS,NDIS
        D=D_IN+DELTA*REAL(IDIS,KIND=8)
        CALL LMTO_EXPANDLOCALWITHCTE('FWRD',IAT,NDIMD,LMNX,LMNXT,D,DT)
        CALL LMTO_CVX_NEW(ISP,LMNXT,EX,DT(:,:,1),HT(:,:,1))
        HT(:,:,2:)=0.D0
        CALL LMTO_EXPANDLOCALWITHCTE('BACK',IAT,NDIMD,LMNX,LMNXT,H,HT)
        EARR(IDIS)=EX
        DARR(IDIS)=SUM(H*DELTA)
      ENDDO

      DO IDIS=1,NDIS
        Y1=0.5D0*(EARR(IDIS)-EARR(-IDIS))/REAL(IDIS,KIND=8)
        Y2=0.5D0*(DARR(IDIS)+DARR(-IDIS))
        PRINT*,'IDIS,Y1,Y2 ',IDIS,Y1,Y2
      ENDDO
      STOP
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_CVX_NEW(ISP,LMNXT,EX,DT,HT)
!     **************************************************************************
!     **  CORE VALENCE EXCHANGE ENERGY                                        **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY: POTPAR1
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: ISP             ! ATOM-TYPE INDEX
      INTEGER(4),INTENT(IN)  :: LMNXT           ! #(LOCAL ORBITALS)
      REAL(8)   ,INTENT(IN)  :: DT(LMNXT,LMNXT)  ! TOTAL DENSITY MATRIX
      REAL(8)   ,INTENT(OUT) :: EX              ! CORE-VALENCE ENERGY
      REAL(8)   ,INTENT(OUT) :: HT(LMNXT,LMNXT)  ! HAMILTONIAN CONTRIBUTION
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: CVXMAT(:,:)
      REAL(8)   ,ALLOCATABLE :: CVX(:,:)
      INTEGER(4)             :: NHEAD
      INTEGER(4)             :: NTAIL
      INTEGER(4)             :: LN1,LN2,L1,L2,LMN1,LMN2,IM
      INTEGER(4)             :: IH1,IH2,IT1,IT2
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
!     == MAP CVXMATIX ELEMENTS ONTO HEAD AND TAIL FUNCTIONS                   ==
!     ==========================================================================
      NHEAD=POTPAR1(ISP)%NHEAD
      NTAIL=POTPAR1(ISP)%NTAIL
      ALLOCATE(CVX(NHEAD+NTAIL,NHEAD+NTAIL))
      CVX=0.D0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO LN2=1,LNX
          L2=LOX(LN2)
          IF(L1.NE.L2) CYCLE
          DO IH1=1,NHEAD
            DO IH2=1,NHEAD
              CVX(IH1,IH2)=CVX(IH1,IH2) &
     &                    +POTPAR1(ISP)%PROK(LN1,IH1)*CVXMAT(LN1,LN2) &
     &                                               *POTPAR1(ISP)%PROK(LN2,IH2)
            ENDDO
            DO IT2=1,NTAIL
              CVX(IH1,NHEAD+IT2)=CVX(IH1,NHEAD+IT2) &
     &                          +POTPAR1(ISP)%PROK(LN1,IH1)*CVXMAT(LN1,LN2) &
     &                                            *POTPAR1(ISP)%PROJBAR(LN2,IT2)
            ENDDO
          ENDDO
          DO IT1=1,NTAIL
            DO IH2=1,NHEAD
              CVX(NHEAD+IT1,IH2)=CVX(NHEAD+IT1,IH2) &
     &                          +POTPAR1(ISP)%PROJBAR(LN1,IT1)*CVXMAT(LN1,LN2) &
     &                                               *POTPAR1(ISP)%PROK(LN2,IH2)
            ENDDO
            DO IT2=1,NTAIL
              CVX(NHEAD+IT1,NHEAD+IT2)=CVX(NHEAD+IT1,NHEAD+IT2) &
     &                          +POTPAR1(ISP)%PROJBAR(LN1,IT1)*CVXMAT(LN1,LN2) &
     &                                            *POTPAR1(ISP)%PROJBAR(LN2,IT2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(CVXMAT)
      DEALLOCATE(LOX)
!
!     ==========================================================================
!     == CALCULATE CORE-VALENCE EXCHANGE ENERGY                               ==
!     ==========================================================================
      LNX=NHEAD+NTAIL
      ALLOCATE(LOX(LNX))
      LOX(:NHEAD)=POTPAR1(ISP)%LOFH
      LOX(NHEAD+1:)=POTPAR1(ISP)%LOFT
!
      EX=0.D0
      HT(:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LMN2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
          IF(L1.EQ.L2) THEN
            DO IM=1,2*L1+1
              EX=EX+CVX(LN2,LN1)*DT(LMN1+IM,LMN2+IM)
              HT(LMN1+IM,LMN2+IM)=HT(LMN1+IM,LMN2+IM)+CVX(LN2,LN1)
            ENDDO
          END IF
          LMN2=LMN2+2*L2+1
        ENDDO
        LMN1=LMN1+2*L1+1
      ENDDO      
!     == HT=0.5D0*(HT+TRANSPOSE(HT))
!
!     ==========================================================================
!     == CLOSE DOWN                                                           ==
!     ==========================================================================
      DEALLOCATE(CVX)
      DEALLOCATE(LOX)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_CVX_ACTONPHI(IAT,LHFWEIGHT,ETOT)
!     **************************************************************************
!     **  CORE VALENCE EXCHANGE ENERGY ACTING ON PARTIAL WAVES                **
!     *****************************PETER BLOECHL, GOSLAR 20116******************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT          ! ATOM INDEX
      REAL(8)   ,INTENT(IN)  :: LHFWEIGHT
      REAL(8)   ,INTENT(OUT) :: ETOT              ! CORE-VALENCE ENERGY
      INTEGER(4)             :: ISP
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: NDIMD
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: CVXMAT(:,:)
      COMPLEX(8),ALLOCATABLE :: DENMAT(:,:,:) !(LMNX,LMNX,NDIMD)
      COMPLEX(8),ALLOCATABLE :: DH(:,:,:) !(LMNX,LMNX,NDIMD)
      INTEGER(4)             :: LN1,LN2,L1,L2,LMN1,LMN2,IM
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_CVX_ACTONPHI')
      CALL LMTOAUGMENTATION$SETI4('IAT',IAT)
!
!     ==========================================================================
!     == COLLECT DATA                                                         ==
!     ==========================================================================
      CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      ALLOCATE(CVXMAT(LNX,LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      CALL SETUP$GETR8A('CVX',LNX*LNX,CVXMAT)
      CALL SETUP$UNSELECT()
      LMNX=SUM(2*LOX+1)
      CALL LMTOAUGMENTATION$GETI4('NDIMD',NDIMD)
      ALLOCATE(DENMAT(LMNX,LMNX,NDIMD))
      CALL LMTOAUGMENTATION$GETC8A('DENMAT',LMNX*LMNX*NDIMD,DENMAT)
!
!     ==========================================================================
!     == SCALE ENERGY TERM                                                    ==
!     ==========================================================================
      CVXMAT=CVXMAT*LHFWEIGHT
!
!     ==========================================================================
!     == SUM UP TOTAL ENERGY                                                  ==
!     ==========================================================================
      ALLOCATE(DH(LMNX,LMNX,NDIMD))
      DH=(0.D0,0.D0)
      ETOT=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LMN2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
          IF(L1.EQ.L2) THEN
            DO IM=1,2*L1+1
              ETOT=ETOT+CVXMAT(LN1,LN2)*REAL(DENMAT(LMN1+IM,LMN2+IM,1),KIND=8)
              DH(LMN1+IM,LMN2+IM,1)=CMPLX(CVXMAT(LN1,LN2),KIND=8)
            ENDDO
          END IF
          LMN2=LMN2+2*L2+1
        ENDDO
        LMN1=LMN1+2*L1+1
      ENDDO        
      CALL LMTOAUGMENTATION$ADDC8A('DH',LMNX*LMNX*NDIMD,DH)
!
!     ==========================================================================
!     == CLOSE DOWN                                                           ==
!     ==========================================================================
      DEALLOCATE(DENMAT)
      DEALLOCATE(DH)
      DEALLOCATE(CVXMAT)
      DEALLOCATE(LOX)
      CALL LMTOAUGMENTATION$SETI4('IAT',0)
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SIMPLEDC_NEW(GID,NR,NDIMD,LMNX,LNX,LOX,CHI,LRX,AECORE &
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
!     **  THE DENSITY MATRIX IS IN A (T,X,Y,Z) REPRESENTATION                 **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: LRX
      INTEGER(4)  ,INTENT(IN) :: NDIMD
      INTEGER(4)  ,INTENT(IN) :: LMNX       ! #(LOCAL ORBITALS)
      INTEGER(4)  ,INTENT(IN) :: LNX        ! #(RADIAL FUNCTIONS)
      INTEGER(4)  ,INTENT(IN) :: LOX(LNX)   !MAIN ANGULAR MOMENTUM OF LOCAL ORB.
      REAL(8)     ,INTENT(IN) :: CHI(NR,LNX)
      REAL(8)     ,INTENT(IN) :: AECORE(NR)
      REAL(8)     ,INTENT(IN) :: DENMAT(LMNX,LMNX,NDIMD) ! DENSITY MATRIX
      REAL(8)     ,INTENT(IN) :: DENMATB(LMNX,LMNX,NDIMD) ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT       ! DOUBLE COUNTINNG ENERGY
      REAL(8)     ,INTENT(OUT):: HAM(LMNX,LMNX,NDIMD)  ! DETOT/D(RHO^*)        
      REAL(8)     ,INTENT(OUT):: HAMB(LMNX,LMNX,NDIMD)  ! DETOT/D(RHO^*)        
      INTEGER(4)  ,PARAMETER  :: METHOD=3
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)     ,PARAMETER  :: FOURPI=4.D0*PI
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
      INTEGER(4)              :: LMRX,L
      INTEGER(4)              :: IDIM,LM,LMN
      REAL(8)                 :: ETOTC,ETOTV
!     **************************************************************************
      LMRX=(LRX+1)**2
      ETOT=0.D0
!
!     ==========================================================================
!     ==  TRANSFORM DENSITY MATRIX FROM UP/DOWN TO TOTAL/SPIN                 ==
!     ==========================================================================
      DENMAT1=CMPLX(DENMAT,KIND=8)
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
!     ==  CONSTRUCT FULL DENSITY   =====================================
      DENMAT1=CMPLX(DENMATB,KIND=8)
      ALLOCATE(RHOWC(NR,LMRX,NDIMD))  !WITH CORE
      DO IDIM=1,NDIMD
        CALL AUGMENTATION_RHO(NR,LNX,LOX,CHI &
     &                     ,LMNX,DENMAT1(:,:,IDIM),LMRX,RHOWC(:,:,IDIM))
      ENDDO
      RHOWC(:,1,1)=RHOWC(:,1,1)+AECORE(:)
!
!     ==  POTENTIAL FOR THE FULL DENSITY =====================================
      ALLOCATE(POT(NR,LMRX,NDIMD))
      CALL LMTO_RADXC(GID,NR,LMRX,NDIMD,RHOWC,FXC,POT)
      AUX(:)=(RHO(:,1,1)/(RHOWC(:,1,1)+1.D-6))**2
      DO IDIM=1,NDIMD
        DO LM=1,LMRX
          POT(:,LM,IDIM)=AUX(:)*POT(:,LM,IDIM)
        ENDDO
      ENDDO
      POT(:,1,1)=POT(:,1,1)-2.D0*FXC(:)*AUX(:)/(RHOWC(:,1,1)*Y0+1.D-6)/Y0
!     == POTENTIAL FOR THE CORRELATED DENSITY =================================
      ALLOCATE(POT2(NR,LMRX,NDIMD))
      POT2(:,:,:)=0.D0
      POT2(:,1,1)=2.D0*FXC*AUX(:)/(RHO(:,1,1)*Y0+1.D-6)/Y0
      CALL RADIAL$R(GID,NR,R)
      AUX(:)=FOURPI*R(:)**2*FXC*AUX(:)
      CALL RADIAL$INTEGRAL(GID,NR,AUX,ETOT)
!      PRINT*,'----EXC  ',ETOT
      CALL LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT2,CHI,HAM1)
      HAM=REAL(HAM1)
      CALL LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
      HAMB=REAL(HAM1)
      DEALLOCATE(POT)
      DEALLOCATE(POT2)
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TEST_RADXC_WITHCUT(GID,NR,LMRX,NDIMD,RHO)
!     **************************************************************************
!     ** TEST ROUTINE FOR LMTO_RADXC                                          **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: LMRX
      INTEGER(4)  ,INTENT(IN) :: NDIMD
      REAL(8)     ,INTENT(IN) :: RHO(NR,LMRX,NDIMD)
      INTEGER(4)  ,PARAMETER  :: NDIS=4
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER  :: FOURPI=4.D0*PI
      REAL(8)                 :: R(NR)
      REAL(8)                 :: AUX(NR)
      REAL(8)                 :: RHO1(NR,LMRX,NDIMD)
      REAL(8)                 :: POT1(NR,LMRX,NDIMD)
      REAL(8)                 :: DRHO(NR,LMRX,NDIMD)
      REAL(8)                 :: SVAR
      REAL(8)                 :: ETOT(-NDIS:NDIS)
      REAL(8)                 :: DETOT(-NDIS:NDIS)
      REAL(8)                 :: CUT1(NR)
      REAL(8)                 :: CUT(NR)
      REAL(8)                 :: DCUT(NR)
      REAL(8)                 :: VCUT(NR)
      INTEGER(4)              :: I,LM,IDIMD,L
!     **************************************************************************
      WRITE(*,FMT='(80("=")/80("="),T10,"  LMTO_RADXC TEST  "/80("="))')
      CALL RADIAL$R(GID,NR,R)
      CUT(:)=EXP(-R**2)
!
!     ==========================================================================
!     ==  DEFINE DISPLACEMENTS                                                ==
!     ==========================================================================
      LM=2
      L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
      DRHO(:,:,:)=0.D0
      DRHO(:,LM,1)=1.D-2*R**L*EXP(-R**2)
      DCUT(:)=0.D-6*EXP(-2.D0*R**2)
!
!     ==========================================================================
!     ==  DATA COLLECTION                                                     ==
!     ==========================================================================
      DO I=-NDIS,NDIS
        RHO1=RHO+DRHO*REAL(I,KIND=8)
        CUT1=CUT+DCUT*REAL(I,KIND=8)
        CALL LMTO_RADXC_WITHCUT(GID,NR,LMRX,NDIMD,RHO1,CUT1,ETOT(I),POT1,VCUT)
        DETOT(I)=0.D0
        DO IDIMD=1,NDIMD
          DO LM=1,LMRX
            AUX=R**2*DRHO(:,LM,IDIMD)*POT1(:,LM,IDIMD)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            DETOT(I)=DETOT(I)+SVAR
          ENDDO
        ENDDO
        AUX=R**2*VCUT*DCUT
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        DETOT(I)=DETOT(I)+SVAR
      ENDDO
!
!     ==========================================================================
!     ==  ANALYZE                                                             ==
!     ==========================================================================
      ETOT=ETOT-ETOT(0)

      WRITE(*,FMT='("NDIMD   ",I5)')NDIMD
      WRITE(*,FMT='("LMRX    ",I5)')LMRX
!
      WRITE(*,FMT='(I5,E15.5,E15.5)')0,0.D0,DETOT(0)
      DO I=1,NDIS
        WRITE(*,FMT='(I5,E15.5,E15.5)')I &
     &             ,(ETOT(I)-ETOT(-I))/(2.D0*REAL(I,KIND=8)) &
     &             ,(DETOT(I)+DETOT(-I))/2.D0
      ENDDO
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_RADXC_WITHCUT(GID,NR,LMRX,NDIMD,RHOIN,CUT,EXC,VXC,VCUT)
!     **************************************************************************
!     **  CALCULATES THE EXCHANGE AND CORRELATION ENERGY                      **
!     **  FOR A DENSITY GIVEN ON A RADIAL LOGARITHMIC GRID                    **
!     **  TIMES REAL SPHERICAL HARMONICS                                      **
!     **                                                                      **
!     **  THIS ROUTINE IS ALMOST IDENTICAL TO THE ROUTINE AUGMENTATION_XC     **
!     **  OF THE AUGMENTATION OBJECT PAW_AUGMENTATION.F90. IT DIFFERS IN THAT **
!     **  THE ENERGY DENSITY IS MULTIPLIED WITH A CUTOFF FUNCTION CUT(R)      **
!     **  SO THAT                                                             **
!     **       EXC=4\PI\INT_0^\INFTY DR R^2 CUT(R)*FXC(R)                     **
!     **  IS THE EXCHANGE-CORRELATION ENERGY.                                 **
!     **                                                                      **
!     **  THE TOTAL ENERGY IS AN EXPANSION ABOUT THE SPHERICAL CONTRIBUTION   **
!     **  OF THE DENSITY UP TO QUADRATIC ORDER IN THE NON-SPHERICAL CONTRIBS. **
!     **                                                                      **
!     **       FXC = FXC(XVAL(L=0)*Y0)                                        **
!     **           + 0.5 * D2[FXC]/DXVAL(L=0)*Y0]**2 * XVAL(L>0)**2           **
!     **                                                                      **
!     **  WHERE XVAL=(/RHOT,RHOS,GRHOT**2,GRHOS**2,GRHOT*GRHOS/)              **
!     **  IS A SPHERICAL HARMONICS EXPANSION ON THE RADIAL GRID.              **
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
!     **          ANGULAR-MOMENTUM EXPANSIONS                                 **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1996 ************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TNS=.TRUE. ! NON-SPHERICAL CONTRIBUTIONS ON
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: NDIMD      ! CAN BE 1,2,4
      REAL(8)   ,INTENT(IN) :: RHOIN(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: CUT(NR)             ! CUTOFF FOR DENSITY
      REAL(8)   ,INTENT(OUT):: EXC                 ! INT:CUT*FXC[RHOIN]
      REAL(8)   ,INTENT(OUT):: VXC(NR,LMRX,NDIMD)  ! DEXC/DRHOIN
      REAL(8)   ,INTENT(OUT):: VCUT(NR)            ! DEXC/DCUT(:)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: FOURPI=4.D0*PI
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)   ,PARAMETER  :: CG0LL=Y0
      LOGICAL(4)            :: TGRA   ! SWITCH FOR GRADIENT CORRECTION
      INTEGER(4)            :: NSPIN
      REAL(8)               :: EXC1
      REAL(8)               :: R(NR)
      REAL(8)               :: FXC(NR)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: GRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VGRHO(:,:,:)
      REAL(8)               :: VAL5(5),VXC5(5),V2XC5(5,5),V3XC5(5,5,5)
      REAL(8)               :: XVAL(NR,5,LMRX)
      REAL(8)               :: XDER(NR,5,LMRX)
      INTEGER(4)            :: IR,L,II,ISPIN,ISPIN1,ISPIN2,I,J
      INTEGER(4)            :: LM
      INTEGER(4)            :: IMAX
      REAL(8)               :: FAC
      REAL(8)               :: WORK(NR)
      REAL(8)               :: WORK1(NR)
      REAL(8)               :: WORK2(NR)
      REAL(8)               :: SVEC(5)
      REAL(8),PARAMETER     :: XX=1.D0  !BUGFIX 160509 XX=0.5->1.0
!     **************************************************************************
                                                  CALL TRACE$PUSH('LMTO_RADXC')
      EXC=0.D0
      VXC(:,:,:)=0.D0
      VCUT(:)=0.D0
!
!     ==========================================================================
!     ==   CALCULATE SOME CONSTANTS NEEDED LATER                              ==
!     ==========================================================================
      CALL DFT$GETL4('GC',TGRA)
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
!     == IMAX ALLOWS TO RESTRICT SOME LOOPS (1:5) TO (1:IMAX) ==================
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
              FAC=REAL(L*(L+1),KIND=8)
              XVAL(:,II,1)=XVAL(:,II,1) &
        &                 +CG0LL*(GRHO(:,LM,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                        +FAC*RHO(:,LM,ISPIN1)*RHO(:,LM,ISPIN2)/R(:)**2)
            ENDDO
            DO LM=2,LMRX 
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              XVAL(:,II,LM)=XVAL(:,II,LM) &
        &                  +XX*CG0LL*(GRHO(:,1,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                            +GRHO(:,LM,ISPIN1)*GRHO(:,1,ISPIN2))
            ENDDO
          ENDDO
        ENDDO
!        XVAL(1,:,1)=XVAL(2,:,1) !AVOID DIVIDEBYZERO
      END IF
!
!     ==========================================================================
!     ==  CALCULATE EXCHANGE ENERGY FOR THE SPHERICAL DENSITY                 ==
!     ==========================================================================
      CALL TRACE$PASS('BEFORE DFT')
      FXC(:)=0.D0
      XDER(:,:,:)=0.D0
      VCUT(:)=0.D0
      DO IR=1,NR
!       ==  CYCLE IF THE TOTAL DENSITY VANISHES ================================
        IF(XVAL(IR,1,1).LE.0.D0) CYCLE
!       == NOW CALL DFT ROUTINE ================================================
        VAL5(:)=XVAL(IR,:,1)*Y0
        CALL DFT3(VAL5,EXC1,VXC5,V2XC5,V3XC5)

!       == NOW CALCULATE ENERGY DENSITY AND DERIVATIVES ========================
        FXC(IR)     =FOURPI*EXC1
        XDER(IR,:,1)=FOURPI*VXC5(:)*Y0
        DO LM=2,LMRX
          DO I=1,IMAX        ! IMAX=<5 
            DO J=1,IMAX
              FXC(IR)=FXC(IR)+0.5D0*XVAL(IR,I,LM)*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,:,1)=XDER(IR,:,1) &
       &                  +0.5D0*Y0*XVAL(IR,I,LM)*V3XC5(:,I,J)*XVAL(IR,J,LM)
              XDER(IR,I,LM)=XDER(IR,I,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,J,LM)=XDER(IR,J,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,I,LM)
            ENDDO
          ENDDO
        ENDDO
        VCUT(IR)    =FXC(IR)
        FXC(IR)     =FXC(IR)*CUT(IR)
        XDER(IR,:,:)=XDER(IR,:,:)*CUT(IR)
      ENDDO
      CALL RADIAL$INTEGRAL(GID,NR,FXC(:)*R(:)**2,EXC)
      CALL TRACE$PASS('AFTER DFT')
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
!           ==  RESOLVE SPHERICAL PART =========================================
            VGRHO(:,1,ISPIN1)=VGRHO(:,1,ISPIN1) &
     &                         +Y0*XDER(:,II,1)*GRHO(:,1,ISPIN2)
            VGRHO(:,1,ISPIN2)=VGRHO(:,1,ISPIN2) &
     &                         +Y0*XDER(:,II,1)*GRHO(:,1,ISPIN1)
!           == NONSPHERICAL CONTRIBUTIONS ======================================
            DO LM=2,LMRX
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=REAL(L*(L+1),KIND=8)
              VGRHO(:,1,ISPIN1)=VGRHO(:,1,ISPIN1) &
     &                         +XX*Y0*XDER(:,II,LM)*GRHO(:,LM,ISPIN2)
              VGRHO(:,1,ISPIN2)=VGRHO(:,1,ISPIN2) &
     &                         +XX*Y0*XDER(:,II,LM)*GRHO(:,LM,ISPIN1)
              VRHO(:,LM,ISPIN1)=VRHO(:,LM,ISPIN1) &
     &                         +FAC*Y0*XDER(:,II,1)*RHO(:,LM,ISPIN2)/R**2
              VRHO(:,LM,ISPIN2)=VRHO(:,LM,ISPIN2) &
     &                         +FAC*Y0*XDER(:,II,1)*RHO(:,LM,ISPIN1)/R**2
              VGRHO(:,LM,ISPIN1)=VGRHO(:,LM,ISPIN1) &
     &                         +Y0*XDER(:,II,1)*GRHO(:,LM,ISPIN2) &
     &                         +XX*Y0*XDER(:,II,LM)*GRHO(:,1,ISPIN2)
              VGRHO(:,LM,ISPIN2)=VGRHO(:,LM,ISPIN2) &
     &                         +Y0*XDER(:,II,1)*GRHO(:,LM,ISPIN1) &
     &                         +XX*Y0*XDER(:,II,LM)*GRHO(:,1,ISPIN1)
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
      SUBROUTINE TESTSIMPLEDC_NEW_NEW(GID,NR,NDIMD,LMNX,LNX,LOX,CHI,LRX &
     &                                ,AECORE,DENMAT,HFSCALE,ETOT,HAM)
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: LRX
      INTEGER(4)  ,INTENT(IN) :: NDIMD
      INTEGER(4)  ,INTENT(IN) :: LMNX       ! #(LOCAL ORBITALS)
      INTEGER(4)  ,INTENT(IN) :: LNX        ! #(RADIAL FUNCTIONS)
      INTEGER(4)  ,INTENT(IN) :: LOX(LNX)   !MAIN ANGULAR MOMENTUM OF LOCAL ORB.
      REAL(8)     ,INTENT(IN) :: CHI(NR,LNX)  ! LOCAL ORBITALS
      REAL(8)     ,INTENT(IN) :: AECORE(NR)
      REAL(8)     ,INTENT(IN) :: DENMAT(LMNX,LMNX,NDIMD) ! DENSITY MATRIX
      REAL(8)     ,INTENT(IN) :: HFSCALE
      REAL(8)     ,INTENT(OUT):: ETOT       ! DOUBLE COUNTINNG ENERGY
      REAL(8)     ,INTENT(OUT):: HAM(LMNX,LMNX,NDIMD)  ! DETOT/D(RHO^*)        
      INTEGER(4)  ,PARAMETER  :: NDIS=4
      INTEGER(4)              :: NAT
      INTEGER(4)              :: IAT
      INTEGER(4)              :: I,IDIMD
      INTEGER(4)              :: ISP
      INTEGER(4)              :: LNX_S
      INTEGER(4)              :: LMNX_S
      INTEGER(4)  ,ALLOCATABLE:: LOX_S(:)
      INTEGER(4)              :: LMNXX
      COMPLEX(8)  ,ALLOCATABLE:: DENMAT_BIG(:,:,:,:)      
      COMPLEX(8)  ,ALLOCATABLE:: DH_BIG(:,:,:,:)      
      COMPLEX(8)  ,ALLOCATABLE:: DENMAT_TOT_SAVE(:,:,:) !LMNX_S,LMNX_S,NDIMD
      COMPLEX(8)  ,ALLOCATABLE:: DIS_TOT(:,:,:)         !LMNX_S,LMNX_S,NDIMD
      COMPLEX(8)  ,ALLOCATABLE:: DENMAT_TOT(:,:,:)      !LMNX_S,LMNX_S,NDIMD
      COMPLEX(8)  ,ALLOCATABLE:: DH_TOT(:,:,:)          !LMNX_S,LMNX_S,NDIMD
      REAL(8)                 :: DENMAT_LOC_SAVE(LMNX,LMNX,NDIMD)     
      REAL(8)                 :: DIS_LOC(LMNX,LMNX,NDIMD)      
      REAL(8)                 :: DENMAT_LOC(LMNX,LMNX,NDIMD)     
      REAL(8)                 :: DH_LOC(LMNX,LMNX,NDIMD)      
      REAL(8)                 :: ETOTARR(-NDIS:NDIS)       
      REAL(8)                 :: DETOTARR(-NDIS:NDIS)       
      REAL(8)                 :: RAN,RAN1,RAN2
      REAL(8)                 :: HFSCALE1=1.D0
      REAL(8)                 :: R(NR)
!     **************************************************************************
      ETOT=0.D0
      HAM=0.D0
      CALL LMTOAUGMENTATION$GETI4('IAT',IAT)
IF(IAT.NE.1) RETURN
      CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('LNX',LNX_S)
      ALLOCATE(LOX_S(LNX_S))
      CALL SETUP$GETI4A('LOX',LNX_S,LOX_S)
      CALL RADIAL$R(GID,NR,R)
      CALL SETUP$UNSELECT()
      LMNX_S=SUM(2*LOX_S+1)
      ALLOCATE(DENMAT_TOT_SAVE(LMNX_S,LMNX_S,NDIMD))
      ALLOCATE(DIS_TOT(LMNX_S,LMNX_S,NDIMD))
      ALLOCATE(DENMAT_TOT(LMNX_S,LMNX_S,NDIMD))
      ALLOCATE(DH_TOT(LMNX_S,LMNX_S,NDIMD))

      CALL LMTOAUGMENTATION$GETI4('NAT',NAT)
      CALL LMTOAUGMENTATION$GETI4('LMNXX',LMNXX)
      ALLOCATE(DENMAT_BIG(LMNXX,LMNXX,NDIMD,NAT))
      ALLOCATE(DH_BIG(LMNXX,LMNXX,NDIMD,NAT))
      CALL LMTOAUGMENTATION$GETC8A('DENMAT',LMNXX*LMNXX*NDIMD*NAT,DENMAT_BIG)
      DENMAT_TOT_SAVE=DENMAT_BIG(:LMNX_S,:LMNX_S,:,IAT)
      DENMAT_LOC_SAVE=DENMAT
!
!     ==========================================================================
!     ==  SET UP DISPLACEMENTS
!     ==========================================================================
      DIS_LOC=0.D0
      DIS_TOT=(0.D0,0.D0)
!      DIS_LOC(1,2,1)=1.D-3
!      DIS_TOT(2,2,1)=1.D-2*(1.D0,0.D0)
      DIS_TOT(1,2,1)=1.D-1*(1.D0,0.D0)
      DO IDIMD=1,NDIMD
        DIS_LOC(:,:,IDIMD)=TRANSPOSE(DIS_LOC(:,:,IDIMD))
        DIS_TOT(:,:,IDIMD)=TRANSPOSE(CONJG(DIS_TOT(:,:,IDIMD)))
      ENDDO
!
!     ==========================================================================
!     ==  TEST CALCULATIONS
!     ==========================================================================
      DO I=-NDIS,NDIS
        DENMAT_LOC=DENMAT_LOC_SAVE+DIS_LOC*REAL(I,KIND=8)
        DENMAT_TOT=DENMAT_TOT_SAVE+DIS_TOT*REAL(I,KIND=8)
        DENMAT_BIG(:,:,:,IAT)=(0.D0,0.D0)
        DENMAT_BIG(:LMNX_S,:LMNX_S,:,IAT)=DENMAT_TOT
        CALL LMTOAUGMENTATION$SETC8A('DENMAT',LMNXX*LMNXX*NDIMD*NAT,DENMAT_BIG)
        CALL LMTO_SIMPLEDC_NEW_NEW(GID,NR,NDIMD,LMNX,LNX,LOX,CHI,LRX &
     &                          ,AECORE,DENMAT_LOC,HFSCALE1,ETOT,DH_LOC)
        CALL LMTOAUGMENTATION$GETC8A('DH',LMNXX*LMNXX*NDIMD*NAT,DH_BIG)
        DH_TOT=DH_BIG(:LMNX_S,:LMNX_S,:,IAT)
        ETOTARR(I)=ETOT
!       ==MINUS SIGN BECAUSE DC IS SUBTRACTED
        DETOTARR(I)=SUM(DIS_LOC*DH_LOC)-REAL(SUM(DIS_TOT*DH_TOT),KIND=8)
      ENDDO
!
!     ==========================================================================
!     ==  ANALYZE
!     ==========================================================================
      ETOTARR=ETOTARR-ETOTARR(0)

      WRITE(*,FMT='("ATOM ID ",I5)')IAT
      WRITE(*,FMT='("LMNX    ",I5)')LMNX
      WRITE(*,FMT='("LNX     ",I5)')LNX
      WRITE(*,FMT='("LOX     ",10I5)')LOX(:)
      WRITE(*,FMT='("NDIMD   ",I5)')NDIMD
      WRITE(*,FMT='("LMNXX   ",I5)')LMNXX
      WRITE(*,FMT='("LMNX_S  ",I5)')LMNX_S
      WRITE(*,FMT='("LNX_S   ",I5)')LNX_S
      WRITE(*,FMT='("LOX_S   ",10I5)')LOX_S(:)
      WRITE(*,FMT='("HFSCALE1",F10.5)')HFSCALE1
!
      WRITE(*,FMT='(I5,E15.5,E15.5)')0,0.D0,DETOTARR(0)
      DO I=1,NDIS
        WRITE(*,FMT='(I5,E15.5,E15.5)')I &
     &             ,(ETOTARR(I)-ETOTARR(-I))/(2.D0*REAL(I,KIND=8)) &
     &             ,(DETOTARR(I)+DETOTARR(-I))/2.D0
      ENDDO
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SIMPLEDC_NEW_NEW(GID,NR,NDIMD,LMNX,LNX,LOX,CHI,LRX &
     &                                ,AECORE,DENMAT,HFSCALE,ETOT,HAM)
!     **************************************************************************
!     **  DOUBLE COUNTING CORRECTION FOR THE HYBRID FUNCTIONAL                **
!     **                                                                      **
!     **  EXPRESSES THE TOTAL DENSITY IN TERMS OF A PARTIAL-WAVE EXPANSION.   **
!     **                                                                      **
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
!     **  THE DENSITY MATRIX IS IN A (T,X,Y,Z) REPRESENTATION                 **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: LRX
      INTEGER(4)  ,INTENT(IN) :: NDIMD
      INTEGER(4)  ,INTENT(IN) :: LMNX       ! #(LOCAL ORBITALS)
      INTEGER(4)  ,INTENT(IN) :: LNX        ! #(RADIAL FUNCTIONS)
      INTEGER(4)  ,INTENT(IN) :: LOX(LNX)   !MAIN ANGULAR MOMENTUM OF LOCAL ORB.
      REAL(8)     ,INTENT(IN) :: CHI(NR,LNX)  ! LOCAL ORBITALS
      REAL(8)     ,INTENT(IN) :: AECORE(NR)
      REAL(8)     ,INTENT(IN) :: DENMAT(LMNX,LMNX,NDIMD) ! DENSITY MATRIX
      REAL(8)     ,INTENT(IN) :: HFSCALE
      REAL(8)     ,INTENT(OUT):: ETOT       ! DOUBLE COUNTINNG ENERGY
      REAL(8)     ,INTENT(OUT):: HAM(LMNX,LMNX,NDIMD)  ! DETOT/D(RHO^*)        
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)     ,PARAMETER  :: DELTA=1.D-8
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER  :: FOURPI=4.D0*PI
      REAL(8)     ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      COMPLEX(8)              :: DENMAT1(LMNX,LMNX,NDIMD)
      COMPLEX(8)              :: HAM1(LMNX,LMNX,NDIMD)
      REAL(8)                 :: R(NR)
      REAL(8)                 :: CUT(NR)
      REAL(8)                 :: VCUT(NR)
      REAL(8)     ,ALLOCATABLE:: RHO_LOC(:,:,:)
      REAL(8)     ,ALLOCATABLE:: POT_LOC(:,:,:)
      REAL(8)     ,ALLOCATABLE:: RHO_ALL(:,:,:)
      REAL(8)     ,ALLOCATABLE:: POT_ALL(:,:,:)
      REAL(8)                 :: AUX(NR),AUX2(NR),SVAR
      INTEGER(4)              :: LMRX,L
      INTEGER(4)              :: IDIM,LM,LMN,IR
!     **************************************************************************
      ETOT=0.D0
      HAM=0.D0
!
      LMRX=(LRX+1)**2
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     ==  DENSITY OF LOCAL ORBITALS (INCLUDING CORE)                          ==
!     ==========================================================================
!     == TAKING THE REAL PART IS NOT APPROXIMATION: ONLY THE REAL PART OF THE 
!     == DENSITY MATRIX IN (TXYZ) REPRESENTATION CONTRIBUTES TO THE DENSITY. 
      DENMAT1=CMPLX(DENMAT,KIND=8)
      ALLOCATE(RHO_LOC(NR,LMRX,NDIMD))
      DO IDIM=1,NDIMD
        CALL AUGMENTATION_RHO(NR,LNX,LOX,CHI &
     &                       ,LMNX,DENMAT1(:,:,IDIM),LMRX,RHO_LOC(:,:,IDIM))
      ENDDO
      RHO_LOC(:,1,1)=RHO_LOC(:,1,1)+AECORE(:)
!
!     ==========================================================================
!     == ONE-CENTER EXPANSION OF THE TOTAL DENSITY                            ==
!     ==========================================================================
      ALLOCATE(RHO_ALL(NR,LMRX,NDIMD))
      CALL LMTOAUGMENTATION$GETRHO(GID,NR,LMRX,NDIMD,RHO_ALL)
      RHO_ALL(:,1,1)=RHO_ALL(:,1,1)+AECORE(:)
!
!     ==========================================================================
!     ==  CUTOFF FUNCTION FOR EXCHANGE-CORRELATION INTEGRAL                   ==
!     ==  CUT IS CLOSE TO UNITY IN THE CENTER AND IS ZERO BEYOND THE ATOM     ==
!     ==========================================================================
      CUT(:)=(RHO_LOC(:,1,1)/(RHO_ALL(:,1,1)+DELTA))**2 
!
!     ==========================================================================
!     ==  CALCULATE ENERGY AND POTENTIAL                                      ==
!     ==========================================================================
      ALLOCATE(POT_ALL(NR,LMRX,NDIMD))  ! POTENTIAL FOR PARTIAL-WAVE DENSITY 
      ALLOCATE(POT_LOC(NR,LMRX,NDIMD))  ! POTENTIAL FOR LOCAL ORBITALS
!
!     == EXCHANGE CORRELATION OF THE TOTAL DENSITY INCLUDING CORE ==============
!     == CALL TEST_RADXC_WITHCUT(GID,NR,LMRX,NDIMD,RHO_ALL)
      CALL LMTO_RADXC_WITHCUT(GID,NR,LMRX,NDIMD,RHO_ALL,CUT,ETOT,POT_ALL,VCUT)
!
!     == SUBTRACT FROZEN-CORE ==================================================
      CALL LMTO_RADXC_WITHCUT(GID,NR,1,1,AECORE,CUT,SVAR,POT_LOC,AUX)
      ETOT=ETOT-SVAR
      VCUT(:)=VCUT(:)-AUX(:)
      POT_LOC=0.D0
!      
!     == POTENTIAL FOR PARTIAL-WAVE DENSITY N_T ================================
      POT_ALL(:,1,1)=POT_ALL(:,1,1)-VCUT(:)*2.D0*CUT(:)/(RHO_ALL(:,1,1)+DELTA) 
!
!     == POTENTIAL FOR THE CORRELATED DENSITY ==================================
      POT_LOC(:,:,:)=0.D0
      POT_LOC(:,1,1)=VCUT(:)*2.D0*RHO_LOC(:,1,1)/(RHO_ALL(:,1,1)+DELTA)**2 
!
!     ==========================================================================
!     ==  EXTRACT HAMILTON CONTRIBUTIONS                                      ==
!     ==========================================================================
      HAM1=(0.D0,0.D0)
      CALL LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT_LOC,CHI,HAM1)
      HAM=REAL(HAM1)
!
!     == POTENTIAL ENTERS WITH NEGATIVE SIGN, BECAUSE THE LOCAL EXCHANGE =======
!     == CONTRIBUTION IS TO BE SUBTRACTED FROM THE ONE-CENTER HAMILTONIAN. =====
!     == SIMILARLY THE SCALE FACTOR HFWEIGHT IS APPLIED AT THIS POINT. =========
      CALL LMTOAUGMENTATION$ADDPOT(GID,NR,LMRX,NDIMD,-POT_ALL*HFSCALE)
!!$PRINT*,'HFSCALE',HFSCALE
!!$PRINT*,'LMRX',LMRX
!!$CALL LMTO_WRITEPHI('POT_LOC.DAT',GID,NR,LMRX,POT_LOC)
!!$CALL LMTO_WRITEPHI('POT_ALL.DAT',GID,NR,LMRX,POT_ALL)
!!$CALL LMTO_WRITEPHI('RHO_LOC.DAT',GID,NR,LMRX,RHO_LOC)
!!$CALL LMTO_WRITEPHI('RHO_ALL.DAT',GID,NR,LMRX,RHO_ALL)
!!$CALL LMTO_WRITEPHI('CUT.DAT',GID,NR,1,CUT)
!!$CALL RADIAL$INTEGRAL(GID,NR,FOURPI*R**2*(RHO_ALL(:,1,1)-AECORE(:))*Y0,SVAR)
!!$PRINT*,'RADXC VALENCE CHARGE (PARTIAL WAVE EXPANSION)=',SVAR
!!$CALL RADIAL$INTEGRAL(GID,NR,FOURPI*R**2*(RHO_LOC(:,1,1)-AECORE(:))*Y0,SVAR)
!!$PRINT*,'RADXC VALENCE CHARGE (LOCAL ORBITALS)=        ',SVAR
!!$CALL RADIAL$INTEGRAL(GID,NR,FOURPI*R**2*AECORE*Y0,SVAR)
!!$PRINT*,'RADXC CORE CHARGE                    =        ',SVAR
!!$STOP 'FORCED'
!
      DEALLOCATE(POT_ALL)
      DEALLOCATE(POT_LOC)
!
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
!     **  THE DENSITY MATRIX IS IN A (T,X,Y,Z) REPRESENTATION                 **
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
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER  :: FOURPI=4.D0*PI
      REAL(8)     ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)                 :: EDENSITY(NR)
      REAL(8)                 :: AUX(NR),SVAR
      REAL(8)                 :: FXC(NR)
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
      ETOT=0.D0
!
!     ==========================================================================
!     ==  TRANSFORM DENSITY MATRIX FROM UP/DOWN TO TOTAL/SPIN                 ==
!     ==========================================================================
      DENMAT1=CMPLX(DENMAT,KIND=8)
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
        CALL LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
        DEALLOCATE(POT)
        HAM=REAL(HAM1)

      ELSE IF(METHOD.EQ.2) THEN 
        DENMAT1=CMPLX(DENMATB,KIND=8)
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
        CALL LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
        HAM=REAL(HAM1)
        POT=POT2-POT
        CALL LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
        HAMB=REAL(HAM1)
        DEALLOCATE(POT)
        DEALLOCATE(POT2)
!
      ELSE IF(METHOD.EQ.3) THEN
!       ==  CONSTRUCT FULL DENSITY   =====================================
        DENMAT1=CMPLX(DENMATB,KIND=8)
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
        CALL LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT2,CHI,HAM1)
        HAM=REAL(HAM1)
        CALL LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
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
        CALL LMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
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
!     **  CALCULATES THE EXCHANGE AND CORRELATION ENERGY                      **
!     **  FOR A DENSITY GIVEN ON A RADIAL LOGARITHMIC GRID                    **
!     **  TIMES REAL SPHERICAL HARMONICS                                      **
!     **                                                                      **
!     **  THIS ROUTINE IS ALMOST IDENTICAL TO THE ROUTINE AUGMENTATION_XC     **
!     **  OF THE AUGMENTATION OBJECT PAW_AUGMENTATION.F90. IT DIFFERS IN THAT **
!     **  THE ENERGY DENSITY IS PROVIDED ON A RADIAL GRID, SO THAT            **
!     **       EXC=4\PI\INT_0^\INFTY DR R^2 FXC(R)                            **
!     **  IS THE EXCHANGE-CORRELATION ENERGY.                                 **
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
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: FOURPI=4.D0*PI
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      REAL(8)   ,PARAMETER  :: CG0LL=Y0
      LOGICAL(4)            :: TGRA   ! SWITCH FOR GRADIENT CORRECTION
      INTEGER(4)            :: NSPIN
      REAL(8)               :: EXC1
      REAL(8)               :: R(NR)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: GRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VRHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: VGRHO(:,:,:)
      REAL(8)               :: VAL5(5),VXC5(5),V2XC5(5,5),V3XC5(5,5,5)
      REAL(8)               :: XVAL(NR,5,LMRX)
      REAL(8)               :: XDER(NR,5,LMRX)
      INTEGER(4)            :: IR,L,II,ISPIN,ISPIN1,ISPIN2,I,J
      INTEGER(4)            :: LM
      INTEGER(4)            :: IMAX
      REAL(8)               :: FAC
      REAL(8)               :: WORK(NR)
      REAL(8)               :: WORK1(NR)
      REAL(8)               :: WORK2(NR)
!     **************************************************************************
                                                  CALL TRACE$PUSH('LMTO_RADXC')
      FXC(:)=0.D0
      VXC(:,:,:)=0.D0
!
!     ==========================================================================
!     ==   CALCULATE SOME CONSTANTS NEEDED LATER                              ==
!     ==========================================================================
      CALL DFT$GETL4('GC',TGRA)
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
PRINT*,'TGRA=',TGRA,' NDIMD=',NDIMD,' NSPIN=',NSPIN,' IMAX=',IMAX

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
              FAC=REAL(L*(L+1),KIND=8)
              XVAL(:,II,1)=XVAL(:,II,1) &
        &                 +CG0LL*(GRHO(:,LM,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                        +FAC*RHO(:,LM,ISPIN1)*RHO(:,LM,ISPIN2)/R(:)**2)
            ENDDO
            DO LM=2,LMRX 
              IF(.NOT.TNS) EXIT ! USED TO RESTORE PREVIOUS STATE
              XVAL(:,II,LM)=XVAL(:,II,LM) &
!        &         +0.5D0*CG0LL*(GRHO(:,1,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &         +CG0LL*(GRHO(:,1,ISPIN1)*GRHO(:,LM,ISPIN2) &
        &                +GRHO(:,LM,ISPIN1)*GRHO(:,1,ISPIN2))
            ENDDO
          ENDDO
        ENDDO
        XVAL(1,:,1)=XVAL(2,:,1) !AVOID DIVIDEBYZERO
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
!       == NOW CALCULATE ENERGY DENSITY AND DERIVATIVES ========================
        FXC(IR)=EXC1
        WORK1(IR)=FOURPI*EXC1
        XDER(IR,:,1)  =FOURPI*VXC5(:)*Y0
        DO LM=2,LMRX
          DO I=1,IMAX        ! IMAX=<5 
            DO J=1,IMAX
              WORK1(IR)=WORK1(IR) &
       &               +0.5D0*XVAL(IR,I,LM)*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,:,1)=XDER(IR,:,1) &
       &                  +0.5D0*Y0*XVAL(IR,I,LM)*V3XC5(:,I,J)*XVAL(IR,J,LM)
              XDER(IR,I,LM)=XDER(IR,I,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,J,LM)
              XDER(IR,J,LM)=XDER(IR,J,LM)+0.5D0*V2XC5(I,J)*XVAL(IR,I,LM)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL TRACE$PASS('AFTER DFT')
!
!     CALL RADIAL$INTEGRAL(GID,NR,WORK1(:)*R(:)**2,EXC)
      FXC=WORK1/FOURPI
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
IF(1.EQ.1) THEN
!           ==  RESOLVE SPHERICAL PART =========================================
            DO LM=1,LMRX
              VGRHO(:,1,ISPIN1)=VGRHO(:,1,ISPIN1) &
     &                         +Y0*XDER(:,II,LM)*GRHO(:,LM,ISPIN2)
              VGRHO(:,1,ISPIN2)=VGRHO(:,1,ISPIN2) &
     &                         +Y0*XDER(:,II,LM)*GRHO(:,LM,ISPIN1)
            ENDDO
!           == NONSPHERICAL CONTRIBUTIONS ======================================
            DO LM=2,LMRX
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=REAL(L*(L+1),KIND=8)
              VRHO(:,LM,ISPIN1)=VRHO(:,LM,ISPIN1) &
     &                         +FAC*Y0*XDER(:,II,1)*RHO(:,LM,ISPIN2)
              VRHO(:,LM,ISPIN2)=VRHO(:,LM,ISPIN2) &
     &                         +FAC*Y0*XDER(:,II,1)*RHO(:,LM,ISPIN1)
              VGRHO(:,1,ISPIN1)=VGRHO(:,1,ISPIN1) &
     &                         +Y0*XDER(:,II,1)*GRHO(:,LM,ISPIN2) &
     &                         +Y0*XDER(:,II,LM)*GRHO(:,1,ISPIN2)
              VGRHO(:,1,ISPIN2)=VGRHO(:,1,ISPIN2) &
     &                         +Y0*XDER(:,II,1)*GRHO(:,LM,ISPIN1) &
     &                         +Y0*XDER(:,II,LM)*GRHO(:,1,ISPIN1)
            ENDDO
ELSE           
!           == FIRST RESOLVE XVAL(:,II,1) ======================================
            DO LM=1,LMRX
              IF(LM.NE.1.AND.(.NOT.TNS)) EXIT ! USED TO RESTORE PREVIOUS STATE
              L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
              FAC=REAL(L*(L+1),KIND=8)
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
END IF
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
!!$           CALL RADIAL$DERIVE(GID,NR,VGRHO(:,LM,ISPIN),WORK2)   !NOT SO 
!!$           WORK1(2:)=2.D0/R(2:)*VGRHO(2:,LM,ISPIN)+WORK2(2:)         !GOOD
!!$           WORK1(1)=WORK1(2)
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
      SUBROUTINE LMTO_OFFSITEXEVAL_NEW(TPARALLEL,EX)
!     **************************************************************************
!     **  OFFSITE EXCHANGE ENERGY WITH TAILED PARTIAL WAVES                   **
!     **                                                                      **
!     **  OUTPUT IS                                                           **
!     **   -- EX                                                              **
!     **   -- HAMIL                                                           **
!     **                                                                      **
!     **  PARALLELIZIZATION MODEL:                                            **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES &
     &                       ,POTPAR=>POTPAR1 &
     &                       ,DENMAT=>DENMAT_NEW &
     &                       ,HAMIL=>HAMIL_NEW &
     &                       ,SBAR=>SBAR_NEW &
     &                       ,NSP &
     &                       ,HYBRIDSETTING &
     &                       ,TCTE
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)  :: TPARALLEL
      REAL(8)   ,INTENT(OUT) :: EX
      REAL(8)                :: RBAS(3,3)
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      INTEGER(4)             :: ISP
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4)             :: NND
      INTEGER(4)             :: NNS
      INTEGER(4)             :: NN,NNA,NNB
      INTEGER(4)             :: ISPA,ISPB
      INTEGER(4)             :: LMN1A,LMN2A,LMN1B,LMN2B,LMN3A,LMN3B
      INTEGER(4)             :: LMNXA,LMNXB
      INTEGER(4)             :: LMNXTA,LMNXTB
      INTEGER(4)             :: LNXA,LNXB
      INTEGER(4),ALLOCATABLE :: NNAT(:)    !(NAT) POINTER TO ONSITE DENMAT
      INTEGER(4),ALLOCATABLE :: INS(:)     !(NAT) POINTER TO ONSITE SBAR
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
      INTEGER(4)             :: LMRX
      INTEGER(4)             :: LRX
      INTEGER(4)             :: LX
      LOGICAL(4)             :: TNDDO,T31,TBONDX
      INTEGER(4)            :: THISTASK,NTASKS
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
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
PRINT*,'============ OFFSITEXEVAL ============================='
      EX=0.D0
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
      ALLOCATE(INS(NAT))
      INS=0
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IATA=SBAR(NN)%IAT1
        IATB=SBAR(NN)%IAT2
        IF(IATA.EQ.IATB.AND.SUM(ABS(SBAR(NN)%IT)).EQ.0) INS(IATA)=NN
      ENDDO
      ALLOCATE(NNAT(NAT))
      NNAT=0
      NND=SIZE(DENMAT)
      DO NN=1,NND
        IATA=DENMAT(NN)%IAT1
        IATB=DENMAT(NN)%IAT2
        IF(IATA.EQ.IATB.AND.SUM(DENMAT(NN)%IT**2).EQ.0) NNAT(IATA)=NN
      ENDDO
!
      DO IAT=1,NAT
        IF(INS(IAT).LE.0) THEN
          CALL ERROR$MSG('INDEX ARRAY FOR ONSITE STRUCTURE CONSTANTS')
          CALL ERROR$MSG('NO ONSITE TERMS FOUND FOR ATOM')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$STOP('LMTO_OFFSITEX')
        END IF
        IF(NNAT(IAT).LE.0) THEN
          CALL ERROR$MSG('INDEX ARRAY FOR ONSITE DENSITY MATRIX')
          CALL ERROR$MSG('NO ONSITE TERMS FOUND FOR ATOM')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$STOP('LMTO_OFFSITEX')
        END IF
      ENDDO
      IF(SIZE(HAMIL).NE.NND) THEN
        CALL ERROR$MSG('SIZE OF DENSITY MATRIX AND HAMILTONIAN INONSISTENT')
        CALL ERROR$STOP('LMTO_OFFSITEX')
      END IF
!
!     ==========================================================================
!     == LOOP OVER PAIRS                                                      ==
!     ==========================================================================
      EX=0.D0
      DO NN=1,NND
        IF(TPARALLEL.AND.MOD(NN-1,NTASKS).NE.THISTASK-1) CYCLE
!
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
        TNDDO=HYBRIDSETTING(ISPA)%TNDDO.AND.HYBRIDSETTING(ISPB)%TNDDO
        T31=HYBRIDSETTING(ISPA)%T31.AND.HYBRIDSETTING(ISPB)%T31
        TBONDX=HYBRIDSETTING(ISPA)%TBONDX.AND.HYBRIDSETTING(ISPB)%TBONDX
!       == DO NOTHING UNLESS NEEDED
        IF(.NOT.(TNDDO.OR.T31.OR.TBONDX)) CYCLE
CALL TIMING$CLOCKON('LOOP:OFFX')
!
!
!       ========================================================================
!       == BLOW UP DENSITY MATRIX                                             ==
!       ========================================================================
CALL TIMING$CLOCKON('OFFX:BLOWUP')
        LMNXTA=POTPAR(ISPA)%TAILED%LMNX
        LMNXTB=POTPAR(ISPB)%TAILED%LMNX
        NDIMD=DENMAT(NN)%N3
        LMNXA=DENMAT(NN)%N1
        LMNXB=DENMAT(NN)%N2
        ALLOCATE(D(LMNXTA,LMNXTB,NDIMD))
        ALLOCATE(DA(LMNXTA,LMNXTA,NDIMD))
        ALLOCATE(DB(LMNXTB,LMNXTB,NDIMD))
        IF(TCTE) THEN
          CALL LMTO_EXPANDNONLOCALWITHCTE('FWRD',IATA,IATB,NDIMD,LMNXA,LMNXB &
     &                          ,LMNXTA,LMNXTB,DENMAT(NN)%MAT,D)
          CALL LMTO_EXPANDLOCALWITHCTE('FWRD',IATA,NDIMD,LMNXA,LMNXTA &
     &                          ,DENMAT(NNA)%MAT,DA)
          CALL LMTO_EXPANDLOCALWITHCTE('FWRD',IATB,NDIMD,LMNXB,LMNXTB &
     &                          ,DENMAT(NNB)%MAT,DB)
        ELSE
          CALL LMTO_EXPANDNONLOCAL('FWRD',NDIMD,LMNXA,LMNXB,LMNXTA,LMNXTB &
     &                          ,SBAR(INS(IATA))%MAT,SBAR(INS(IATB))%MAT &
     &                          ,DENMAT(NN)%MAT,D)
          CALL LMTO_EXPANDLOCAL('FWRD',NDIMD,LMNXA,LMNXTA,SBAR(INS(IATA))%MAT &
     &                          ,DENMAT(NNA)%MAT,DA)
          CALL LMTO_EXPANDLOCAL('FWRD',NDIMD,LMNXB,LMNXTB,SBAR(INS(IATB))%MAT &
     &                          ,DENMAT(NNB)%MAT,DB)
        END IF
CALL TIMING$CLOCKOFF('OFFX:BLOWUP')
 !
!       ========================================================================
!       == ROTATE DENSITY MATRIX SO THAT DISTANCE VECTOR POINTS IN Z-DIRECTION=
!       ========================================================================
CALL TIMING$CLOCKON('OFFX:ROTATE1')
!       == CONSTRUCT ROTATION MATRIX ===========================================
!       == DISTANCE VECTOR WILL BE NEW Z-DIRECTION =============================
        DR(:)=RB(:)-RA(:)
        DIS=SQRT(SUM(DR**2))
        ROT(:,3)=DR(:)/DIS  ! ONSITE TERMS ARE ALREADY EXCLUDED
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
        ALLOCATE(UROTA(LMNXTA,LMNXTA))
        ALLOCATE(UROTB(LMNXTB,LMNXTB))
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
CALL TIMING$CLOCKOFF('OFFX:ROTATE1')
!
!       ========================================================================
!       == ADD UP EXCHANGE ENERGY                                             ==
!       ========================================================================
        ALLOCATE(H (LMNXTA,LMNXTB,NDIMD))
        ALLOCATE(HA(LMNXTA,LMNXTA,NDIMD))
        ALLOCATE(HB(LMNXTB,LMNXTB,NDIMD))
        H=0.D0
        HA=0.D0
        HB=0.D0
!
!       ========================================================================
!       == NDDO TERM:                                                         ==
!       == U22(1,2,3,4)=INT DX1 INT DX2: A1(X1)*A2(X1)][B3(X2)*B4(X2)]/|R1-R2|==
!       == LMN1A,LMN2A TIED TO COORDINATE X1, LMN1B,LMN2B TO X2               ==
!       ========================================================================
        IF(TNDDO) THEN
CALL TIMING$CLOCKON('OFFX:NDDO')
          ALLOCATE(U22   (LMNXTA,LMNXTA,LMNXTB,LMNXTB))
          ALLOCATE(DU22  (LMNXTA,LMNXTA,LMNXTB,LMNXTB))
          DIS=SQRT(SUM((RB-RA)**2))
          CALL LMTO_OFFSITEX22U(ISPA,ISPB, DIS,LMNXTA,LMNXTB,U22,DU22)
          DO LMN1B=1,LMNXTB
            DO LMN2B=1,LMNXTB
              DO LMN2A=1,LMNXTA
                DO LMN1A=1,LMNXTA
                  SVAR=-0.25D0*U22(LMN1A,LMN2A,LMN2B,LMN1B)
                  EX=EX+SVAR*SUM(D(LMN1A,LMN1B,:)*D(LMN2A,LMN2B,:))
                  H(LMN1A,LMN1B,:)=H(LMN1A,LMN1B,:)+SVAR*D(LMN2A,LMN2B,:)
                  H(LMN2A,LMN2B,:)=H(LMN2A,LMN2B,:)+SVAR*D(LMN1A,LMN1B,:)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(U22)
          DEALLOCATE(DU22)
CALL TIMING$CLOCKOFF('OFFX:NDDO')
        END IF
!
!       ========================================================================
!       == 3-1 TERMS:                                                         ==
!       == U3A1B(1,2,3,4)=INT DX1 INT DX2:                                    ==
!       ==                          * [A1(X1)*A2(X1)][A3(X2)*B4(X2)]/|X1-X2|  ==
!       == LMN1A,LMN2A TIED TO COORDINATE X1, LMN3A,LMN1B TO X2               ==
!       ========================================================================
        IF(T31) THEN
CALL TIMING$CLOCKON('OFFX:31')
          ALLOCATE(U3A1B (LMNXTA,LMNXTA,LMNXTA,LMNXTB))
          ALLOCATE(DU3A1B(LMNXTA,LMNXTA,LMNXTA,LMNXTB))
          DIS=SQRT(SUM((RB-RA)**2))
          CALL LMTO_OFFSITEX31U(ISPA,ISPB, DIS,LMNXTA,LMNXTB,U3A1B,DU3A1B)
          DO LMN1B=1,LMNXTB
            DO LMN3A=1,LMNXTA 
              DO LMN2A=1,LMNXTA
                DO LMN1A=1,LMNXTA
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
          DEALLOCATE(U3A1B)
          DEALLOCATE(DU3A1B)
          ALLOCATE(U3B1A (LMNXTB,LMNXTB,LMNXTB,LMNXTA))
          ALLOCATE(DU3B1A(LMNXTB,LMNXTB,LMNXTB,LMNXTA))
          CALL LMTO_OFFSITEX31U(ISPB,ISPA,-DIS,LMNXTB,LMNXTA,U3B1A,DU3B1A)
          DO LMN1A=1,LMNXTA
            DO LMN3B=1,LMNXTB 
              DO LMN2B=1,LMNXTB
                DO LMN1B=1,LMNXTB
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
          DEALLOCATE(U3B1A)
          DEALLOCATE(DU3B1A)
CALL TIMING$CLOCKOFF('OFFX:31')
        END IF
!
!       ========================================================================
!       == 2ND ORDER DIFFERENTIAL OVERLAP TERMS:                              ==
!       == BONDU(1,2,3,4)=INT DX1 INT DX2:                                    ==
!       ==                          * [A1(X1)*B2(X1)][A3(X2)*B4(X2)]/|X1-X2|  ==
!       == LMN1A,LMN1B TIED TO COORDINATE X1, LMN2A,LMN2B TO X2               ==
!       ========================================================================
        IF(TBONDX) THEN
CALL TIMING$CLOCKON('OFFX:BONDX')
          ALLOCATE(BONDU (LMNXTA,LMNXTB,LMNXTA,LMNXTB))
          ALLOCATE(DBONDU(LMNXTA,LMNXTB,LMNXTA,LMNXTB))
          DIS=SQRT(SUM((RB-RA)**2))
!         == BONDU(1,2,3,4)=INT DX INT DX': A1(X)*B2(X)][A3(X')*B4(X')]/|R-R'|==
          CALL LMTO_OFFSITEXBONDU(ISPA,ISPB,DIS,LMNXTA,LMNXTB,BONDU,DBONDU)
          DO LMN2B=1,LMNXTB
            DO LMN2A=1,LMNXTA 
              DO LMN1B=1,LMNXTB
                DO LMN1A=1,LMNXTA
                  SVAR=-0.25D0*BONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                  EX=EX+SVAR*SUM(DA(LMN1A,LMN2A,:)*DB(LMN1B,LMN2B,:))
                  HA(LMN1A,LMN2A,:)=HA(LMN1A,LMN2A,:)+SVAR*DB(LMN1B,LMN2B,:)
                  HB(LMN1B,LMN2B,:)=HB(LMN1B,LMN2B,:)+SVAR*DA(LMN1A,LMN2A,:)
                  SVAR=-0.25D0*BONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                  EX=EX+SVAR*SUM(D(LMN1A,LMN2B,:)*D(LMN2A,LMN1B,:))
                  H(LMN1A,LMN2B,:)=H(LMN1A,LMN2B,:)+SVAR*D(LMN2A,LMN1B,:)
                  H(LMN2A,LMN1B,:)=H(LMN2A,LMN1B,:)+SVAR*D(LMN1A,LMN2B,:)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(BONDU)
          DEALLOCATE(DBONDU)
CALL TIMING$CLOCKOFF('OFFX:BONDX')
        END IF
!
!       ========================================================================
!       == ROTATE HAMILTONIAN BACK                                            ==
!       ========================================================================
        DO I=1,NDIMD
          HA(:,:,I)=MATMUL(UROTA,MATMUL(HA(:,:,I),TRANSPOSE(UROTA)))
          H(:,:,I) =MATMUL(UROTA,MATMUL(H(:,:,I) ,TRANSPOSE(UROTB)))
          HB(:,:,I)=MATMUL(UROTB,MATMUL(HB(:,:,I),TRANSPOSE(UROTB)))
        ENDDO
        DEALLOCATE(UROTA)
        DEALLOCATE(UROTB)
!
!       ========================================================================
!       == SHRINK DOWN HAMILTONIAN                                            ==
!       ========================================================================
!       __D,DA IS RE-USED AS WORK ARRAY TO COLLECT HAMILTONIANS______________
CALL TIMING$CLOCKON('OFFX:SHRINKDOWN')
        D=0.D0
        DA=0.D0
        DB=0.D0
        IF(TCTE) THEN
          CALL LMTO_EXPANDNONLOCALWITHCTE('BACK',IATA,IATB,NDIMD &
     &                          ,LMNXA,LMNXB,LMNXTA,LMNXTB &
     &                          ,D(:LMNXA,:LMNXB,:),H)
          CALL LMTO_EXPANDLOCALWITHCTE('BACK',IATA,NDIMD,LMNXA,LMNXTA &
     &                          ,DA(:LMNXA,:LMNXA,:),HA)
          CALL LMTO_EXPANDLOCALWITHCTE('BACK',IATB,NDIMD,LMNXB,LMNXTB &
     &                          ,DB(:LMNXB,:LMNXB,:),HB)
        ELSE
          CALL LMTO_EXPANDNONLOCAL('BACK',NDIMD,LMNXA,LMNXB,LMNXTA,LMNXTB &
     &                          ,SBAR(INS(IATA))%MAT,SBAR(INS(IATB))%MAT &
     &                          ,D(:LMNXA,:LMNXB,:),H)
          CALL LMTO_EXPANDLOCAL('BACK',NDIMD,LMNXA,LMNXTA,SBAR(INS(IATA))%MAT &
     &                          ,DA(:LMNXA,:LMNXA,:),HA)
          CALL LMTO_EXPANDLOCAL('BACK',NDIMD,LMNXB,LMNXTB,SBAR(INS(IATB))%MAT &
     &                          ,DB(:LMNXB,:LMNXB,:),HB)
        END IF
        HAMIL(NN)%MAT(:,:,:) =HAMIL(NN)%MAT(:,:,:) +D(:LMNXA,:LMNXB,:)
        HAMIL(NNA)%MAT(:,:,:)=HAMIL(NNA)%MAT(:,:,:)+DA(:LMNXA,:LMNXA,:)
        HAMIL(NNB)%MAT(:,:,:)=HAMIL(NNB)%MAT(:,:,:)+DB(:LMNXB,:LMNXB,:)
        DEALLOCATE(D)
        DEALLOCATE(H)
        DEALLOCATE(DA)
        DEALLOCATE(HA)
        DEALLOCATE(DB)
        DEALLOCATE(HB)
CALL TIMING$CLOCKOFF('OFFX:SHRINKDOWN')
!!$WRITE(*,*)XDELTA*REAL(IX,8),EX,HAMIL(NN)%MAT(IND1,IND2,IND3)
!!$IF(IX.EQ.3) STOP 'FORCED'
!!$GOTO 1000
CALL TIMING$CLOCKOFF('LOOP:OFFX')
      ENDDO
!
!     ==========================================================================
!     == CLEAN UP TO AVOID MEMORY LEAK                                        ==
!     ==========================================================================
      DEALLOCATE(NNAT)
      DEALLOCATE(INS)
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFXINT()
!     **************************************************************************
!     ** COMPUTES OFFSITE MATRIX ELEMENTS OFFSITEX
!     ** WATCH PARALLELIZATION!!!
!     **************************************************************************
      USE LMTO_MODULE,ONLY : POTPAR=>POTPAR1 &
     &                      ,OFFSITEX &
     &                      ,HYBRIDSETTING &
     &                      ,NSP
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
      DO ISP1=1,NSP
        DO ISP2=1,NSP
          NULLIFY(OFFSITEX(ISP1,ISP2)%OVERLAP)
          NULLIFY(OFFSITEX(ISP1,ISP2)%X22)
          NULLIFY(OFFSITEX(ISP1,ISP2)%X31)
          NULLIFY(OFFSITEX(ISP1,ISP2)%BONDU)
        ENDDO
      ENDDO
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
        IF(.NOT.HYBRIDSETTING(ISP1)%TNDDO) CYCLE
        GID1=POTPAR(ISP1)%TAILED%GID
        CALL RADIAL$GETI4(GID1,'NR',NR1)
        DO ISP2=1,NSP
          IF(.NOT.HYBRIDSETTING(ISP2)%TNDDO) CYCLE
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
        IF(.NOT.HYBRIDSETTING(ISP1)%T31) CYCLE
        GID1=POTPAR(ISP1)%TAILED%GID
        CALL RADIAL$GETI4(GID1,'NR',NR1)
        DO ISP2=1,NSP
          IF(.NOT.HYBRIDSETTING(ISP2)%T31) CYCLE
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
CALL TRACE$PASS('DOING XABAB')
      CALL LMTO_TAILEDGAUSSOFFSITEU()   ! ROUTINE PARALLELIZES OVER 'MONOMER'
CALL TRACE$PASS('XABAB DONE')
!
!     ==========================================================================
!     == CONVERT INTEGRALS INTO COEFFICIENTS OF INTERPOLATING FUNCTION        ==
!     ==========================================================================
PRINT*,'CONVERTING....'
CALL TRACE$PASS('CONVERTING')
      CALL LMTO_OFFSITEXCONVERT()
CALL TRACE$PASS('CONVERSION DONE')
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
      USE LMTO_MODULE, ONLY : POTPAR=>POTPAR1 &
     &                       ,OFFSITEX
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      INTEGER(4),INTENT(IN) :: NR1
      INTEGER(4),INTENT(IN) :: NR2
      REAL(8)   ,INTENT(IN) :: TOLERANCE
      REAL(8)   ,PARAMETER  :: TOLMIN=1.D-8
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
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
REAL(8)::SVAR
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
              CALL LMTO_TWOCENTER(L1,MABS,GID1,NR1,PHI1 &
       &                         ,L2,MABS,GID2,NR2,PHI2 &
       &                         ,DIS,TOL,INTEGRAL)
              OFFSITEX(ISP1,ISP2)%OVERLAP(IDIS,IND)=INTEGRAL
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL MPE$COMBINE('MONOMER','+',OFFSITEX(ISP1,ISP2)%OVERLAP)
!
!     ==========================================================================
!     == PRINT FOR TESTING                                                    ==
!     ==========================================================================
      IF(TPR) THEN
        IND=0
        DO LN1=1,LNX1
          L1=POTPAR(ISP1)%TAILED%LOX(LN1)
          DO LN2=1,LNX2
            L2=POTPAR(ISP2)%TAILED%LOX(LN2)
            DO MABS=0,MIN(L1,L2)
              IND=IND+1
              DO IDIS=1,NDIS
                DIS=OFFSITEX(ISP1,ISP2)%DIS(IDIS)
                WRITE(*,*)ISP1,ISP2,LN1,LN2,MABS,IDIS, &
     &                      OFFSITEX(ISP1,ISP2)%OVERLAP(IDIS,IND)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEX22SETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR=>POTPAR1,OFFSITEX
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      INTEGER(4),INTENT(IN) :: NR1
      INTEGER(4),INTENT(IN) :: NR2
      REAL(8)   ,INTENT(IN) :: TOLERANCE
      REAL(8)   ,PARAMETER  :: RX=0.1D0
      REAL(8)   ,PARAMETER  :: TOLMIN=1.D-8
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: SQ4PI=SQRT(4.D0*PI)
      REAL(8)   ,PARAMETER  :: SQ4PITHIRD=SQRT(4.D0*PI/3.D0)
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
      INTEGER(4)            :: THISTASK,NTASKS,COUNT
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
      USE LMTO_MODULE, ONLY : POTPAR=>POTPAR1,OFFSITEX
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
      USE LMTO_MODULE, ONLY : OFFSITEX,POTPAR=>POTPAR1,NSP
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
      LOGICAL(4)           :: TOVERLAP,TNDDO,T31,TBONDU
!     **************************************************************************
      DO ISP1=1,NSP
        DO ISP2=1,NSP
!         == CHECK WHETHER THERE IS SOMETHING TO CONVERT =======================
          TOVERLAP=ASSOCIATED(OFFSITEX(ISP1,ISP2)%OVERLAP)
          TNDDO   =ASSOCIATED(OFFSITEX(ISP1,ISP2)%X22)
          T31     =ASSOCIATED(OFFSITEX(ISP1,ISP2)%X31)
          TBONDU  =ASSOCIATED(OFFSITEX(ISP1,ISP2)%BONDU)
          IF(.NOT.(TOVERLAP.OR.TNDDO.OR.T31.OR.TBONDU)) CYCLE
!
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
          IF(TOVERLAP) THEN
            NIND=SIZE(OFFSITEX(ISP1,ISP2)%OVERLAP(1,:))
            DO I=1,NIND
              OFFSITEX(ISP1,ISP2)%OVERLAP(:NF,I) &
     &                         =MATMUL(AMATIN,OFFSITEX(ISP1,ISP2)%OVERLAP(:,I))
            ENDDO
            OFFSITEX(ISP1,ISP2)%OVERLAP(NF+1:,:)=0.D0
          END IF
!
!         ======================================================================
!         == COULOMB TERMS X22
!         ======================================================================
          IF(TNDDO) THEN
            NIND=SIZE(OFFSITEX(ISP1,ISP2)%X22(1,:))
            DO I=1,NIND
              OFFSITEX(ISP1,ISP2)%X22(:NF,I)=MATMUL(AMATIN &
      &                                          ,OFFSITEX(ISP1,ISP2)%X22(:,I))
            ENDDO
            OFFSITEX(ISP1,ISP2)%X22(NF+1:,:)=0.D0
          END IF
!
!         ======================================================================
!         == 31 TERMS  X31
!         ======================================================================
          IF(T31) THEN
            NIND=SIZE(OFFSITEX(ISP1,ISP2)%X31(1,:))
            DO I=1,NIND
              OFFSITEX(ISP1,ISP2)%X31(:NF,I)=MATMUL(AMATIN &
       &                                         ,OFFSITEX(ISP1,ISP2)%X31(:,I))
            ENDDO
            OFFSITEX(ISP1,ISP2)%X31(NF+1:,:)=0.D0
          END IF
!
!         ======================================================================
!         == BOND OVERLAP INTERACTIONS BONDU                                  ==
!         ======================================================================
          IF(TBONDU) THEN
            NIND=SIZE(OFFSITEX(ISP1,ISP2)%BONDU(1,:))
            DO I=1,NIND
              OFFSITEX(ISP1,ISP2)%BONDU(:NF,I)=MATMUL(AMATIN &
     &                           ,OFFSITEX(ISP1,ISP2)%BONDU(:,I))
            ENDDO
            OFFSITEX(ISP1,ISP2)%BONDU(NF+1:,:)=0.D0
          END IF
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
      USE LMTO_MODULE, ONLY : POTPAR=>POTPAR1,OFFSITEX
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: ISP1=1
      INTEGER(4),PARAMETER :: ISP2=1
      INTEGER(4)           :: IND
      INTEGER(4)           :: LMN1,LMN2,LMN3,LMN4
      INTEGER(4)           :: LMNXA,LMNXB
      INTEGER(4)           :: IAB,ICD
      REAL(8)   ,ALLOCATABLE :: UABCD(:,:,:,:)
!     **************************************************************************
      IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%BONDU)) THEN
        CALL ERROR$MSG('OFFSITEX(ISP1,ISP2)%BONDU NOT INITIALIZED')
        CALL ERROR$I4VAL('ISP1',ISP1)
        CALL ERROR$I4VAL('ISP2',ISP2)
        CALL ERROR$STOP('LMTO_PRBONDU')
      END IF
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
          WRITE(*,FMT='(82("="),T20," BONDU FOR LMN1=",I3," AND LMN2=",I3)') &
    &                         LMN1,LMN2
          DO LMN3=1,LMNXA
            WRITE(*,FMT='(I3,50F10.5)')LMN3,UABCD(LMN1,LMN2,LMN3,:)
          ENDDO
        ENDDO
      ENDDO
      CALL ERROR$MSG('INTENDED STOP AFTER PRINTING OFFSITEX#BONDU')
      CALL ERROR$STOP('LMTO_PRBONDU')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEOVERLAP(ISP1,ISP2,DIS,LMNX1,LMNX2,O,DO)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS                                               **
!     ** TWO ORBITALS ARE ON THE FIRST AND TWO ORBITALS ON THE SECOND ATOM    **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR=>POTPAR1
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
                                  CALL TRACE$PUSH('LMTO_OFFSITEOVERLAP')
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
!     == OBTAIN OVERLAP MATRIX ELEMENTS                                       ==
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
                                  CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEX22U(ISP1,ISP2,DIS,LMNX1,LMNX2,U,DU)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS                                               **
!     ** TWO ORBITALS ARE ON THE FIRST AND TWO ORBITALS ON THE SECOND ATOM    **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR=>POTPAR1
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      REAL(8)   ,INTENT(IN) :: DIS
      INTEGER(4),INTENT(IN) :: LMNX1
      INTEGER(4),INTENT(IN) :: LMNX2
      REAL(8)   ,INTENT(OUT):: U(LMNX1,LMNX1,LMNX2,LMNX2)
      REAL(8)   ,INTENT(OUT):: DU(LMNX1,LMNX1,LMNX2,LMNX2)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: SQ4PI=SQRT(4.D0*PI)
      REAL(8)   ,PARAMETER  :: SQ4PITHIRD=SQRT(4.D0*PI/3.D0)
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
      REAL(8)               :: SVAR
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
      USE LMTO_MODULE, ONLY : POTPAR=>POTPAR1
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
      REAL(8)               :: FAC,DFAC
      INTEGER(4),PARAMETER  :: VERSION=2
!         VERSION=0  ORIGINAL VERSION
!         VERSION=2  OWN SUBROUTINE VERSION
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
!!$WRITE(*,FMT='(80("="),T20,"  OFFX31  ")')
!!$WRITE(*,FMT='("ISP1/2   ",T20,2I5)')ISP1,ISP2
!!$WRITE(*,FMT='("LNX1/2   ",T20,2I5)')LNX1,LNX2
!!$WRITE(*,FMT='("LOX(ISP1)",T20,20I5)')POTPAR(ISP1)%TAILED%LOX(:)
!!$WRITE(*,FMT='("LOX(ISP2)",T20,20I5)')POTPAR(ISP2)%TAILED%LOX(:)
!!$WRITE(*,FMT='("LMN0A    ",T20,20I5)')LMN0A
!!$WRITE(*,FMT='("LMN0B    ",T20,20I5)')LMN0B
IF(VERSION.EQ.0) THEN
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
                      LM4=L4**2+L4+1+M
                      LMN4=LMN0B(LN4)+L4+1+M
                      LMR2=LR2**2+LR2+1+M

                      LM2=L2**2
                      LMN2=LMN0A(LN2)
                      DO M2=1,2*L2+1
                        LM2=LM2+1
                        LMN2=LMN2+1

                        LM1=L1**2
                        LMN1=LMN0A(LN1)
                        DO M1=1,2*L1+1
                          LM1=LM1+1
                          LMN1=LMN1+1

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
                              CALL SPHERICAL$GAUNT(LMR1,LM3,LMR2,CG2)
                              IF(CG2.EQ.0.D0) CYCLE
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
ELSE IF(VERSION.EQ.1) THEN
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
                      LM4=L4**2+L4+1+M
                      LMN4=LMN0B(LN4)+L4+1+M
                      LMR2=LR2**2+LR2+1+M

                      LM3=L3**2
                      LMN3=LMN0A(LN3)
                      DO M3=1,2*L3+1
                        LM3=LM3+1
                        LMN3=LMN3+1
!
                        LMR1=LR1**2
                        DO MR1=1,2*LR1+1
                          LMR1=LMR1+1
                         CALL SPHERICAL$GAUNT(LMR1,LM3,LMR2,CG2)
                          IF(CG2.EQ.0.D0) CYCLE
                          FAC=CG2*INTEGRAL
                          DFAC=CG2*DINTEGRAL
!
                          LM2=L2**2
                          LMN2=LMN0A(LN2)
                          DO M2=1,2*L2+1
                            LM2=LM2+1
                            LMN2=LMN2+1
!
                            LM1=L1**2
                            LMN1=LMN0A(LN1)
                            DO M1=1,2*L1+1
                              LM1=LM1+1
                              LMN1=LMN1+1
                              CALL SPHERICAL$GAUNT(LM1,LM2,LMR1,CG1)
                              IF(CG1.EQ.0.D0) CYCLE
!
                              U(LMN1,LMN2,LMN3,LMN4) =U(LMN1,LMN2,LMN3,LMN4) &
       &                                             +CG1*FAC
                              DU(LMN1,LMN2,LMN3,LMN4)=DU(LMN1,LMN2,LMN3,LMN4) &
       &                                             +CG1*DFAC
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
ELSE IF(VERSION.EQ.2) THEN
      CALL LMTO_OFFSITEX31U_TEST1(ISP1,LNX1,POTPAR(ISP1)%TAILED%LOX,LMNX1 &
     &                           ,ISP2,LNX2,POTPAR(ISP2)%TAILED%LOX,LMNX2 &
     &                           ,LRX1,DIS,U,DU)
ELSE
  CALL ERROR$MSG('VALUE OF VERSION NOT RECOGNIZED')
  CALL ERROR$STOP('VALUE OF VERSION NOT RECOGNIZED')
END IF
!!$  END IF
!!$END IF

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
      SUBROUTINE LMTO_OFFSITEX31U_TEST1(ISPA,LNXA,LOXA,LMNXA &
     &                                 ,ISPB,LNXB,LOXB,LMNXB &
     &                                 ,LRX,DIS,U,DU)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISPA
      INTEGER(4),INTENT(IN) :: ISPB
      INTEGER(4),INTENT(IN) :: LNXA
      INTEGER(4),INTENT(IN) :: LNXB
      INTEGER(4),INTENT(IN) :: LOXA(LNXA)
      INTEGER(4),INTENT(IN) :: LOXB(LNXB)
      INTEGER(4),INTENT(IN) :: LMNXA
      INTEGER(4),INTENT(IN) :: LMNXB
      INTEGER(4),INTENT(IN) :: LRX !X(ANGULAR MOMENTUM FOR DENSITY)
      REAL(8)   ,INTENT(IN) :: DIS
      REAL(8)   ,INTENT(OUT):: U(LMNXA,LMNXA,LMNXA,LMNXB)
      REAL(8)   ,INTENT(OUT):: DU(LMNXA,LMNXA,LMNXA,LMNXB)
      TYPE X31LIST_TYPE
        INTEGER(4) :: M1
        INTEGER(4) :: M2
        INTEGER(4) :: M3
        INTEGER(4) :: M4
        REAL(8)    :: VAL
      END TYPE X31LIST_TYPE
      TYPE(X31LIST_TYPE),ALLOCATABLE :: LIST(:) !(IIX)
      REAL(8)   ,ALLOCATABLE :: INTEGRAL(:)  !(INDX)
      REAL(8)   ,ALLOCATABLE :: DINTEGRAL(:) !(INDX)
      REAL(8)   ,ALLOCATABLE :: CGMAT(:,:,:) !GAUND COEFFICIENTS
      INTEGER(4),ALLOCATABLE :: IIA1(:,:,:,:,:,:,:) !(LXA+1,LXA+1,LXA+1,...
      INTEGER(4),ALLOCATABLE :: IIA2(:,:,:,:,:,:,:) !(LXA+1,LXA+1,LXA+1,...
      INTEGER(4),ALLOCATABLE :: IIBNDS(:,:) !(2,IIX)
      INTEGER(4)             :: LXA,LXB
      INTEGER(4)             :: LMXA,LMXB
      INTEGER(4)             :: LMN0A(LNXA),LMN0B(LNXB)
      INTEGER(4)             :: MABSX
      INTEGER(4)             :: INDX
      INTEGER(4)             :: LR1X,LR2X
      INTEGER(4)             :: IND
      INTEGER(4)             :: LN1,LN2,LN3,LN4
      INTEGER(4)             :: L1,L2,L3,L4
      INTEGER(4)             :: M1,M2,M3,M4
      INTEGER(4)             :: LM1,LM2,LM3,LM4
      INTEGER(4)             :: LMN1,LMN2,LMN3,LMN4
      INTEGER(4)             :: LR1,LR2,MABS
      INTEGER(4)             :: LMR1,LMR2
      INTEGER(4)             :: M,PM,MR1
      INTEGER(4)             :: LMR1I,LMR1F
      INTEGER(4)             :: IIX
      INTEGER(4)             :: II1,II2,II
      REAL(8)                :: SVAR
      REAL(8)                :: CG1,CG2
      LOGICAL(4)             :: TCHK
!     **************************************************************************
!
!     ==========================================================================
!     == CALCULATE ARRAY SIZES                                                ==
!     ==========================================================================
      LXA=MAXVAL(LOXA)
      LXB=MAXVAL(LOXB)
      LMXA=LXA**2
      LMXB=LXB**2
      INDX=0
      LR1X=-1
      LR2X=-1
      MABSX=-1
      DO LN1=1,LNXA
        L1=LOXA(LN1)
        DO LN2=LN1,LNXA
          L2=LOXA(LN2)
          DO LN3=1,LNXA
            L3=LOXA(LN3)
            DO LN4=1,LNXB
              L4=LOXB(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX),2
                LR1X=MAX(LR1X,LR1)
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  LR2X=MAX(LR2X,LR2)
                  DO MABS=0,MIN(LR2,L4)
                    MABSX=MAX(MABSX,MIN(LR2,L4))
                    INDX=INDX+1
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     THIS ESTIMATE IS NEARLY A FACTOR 10 TOO LARGE!!! WHY?
      IIX=0
      DO L1=0,LXA
        DO L2=0,LXA
          DO L3=0,LXA
            DO L4=0,LXB
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX),2
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
                    IIX=IIX+(2*L1+1)*(2*L2+1)*(2*L3+1)*2
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     == SET UP INDEX ARRAY ====================================================
      LMN1=0
      DO LN1=1,LNXA
        L1=LOXA(LN1)
        LMN0A(LN1)=LMN1
        LMN1=LMN1+2*L1+1
      ENDDO
      LMN1=0
      DO LN1=1,LNXB
        L1=LOXB(LN1)
        LMN0B(LN1)=LMN1
        LMN1=LMN1+2*L1+1
      ENDDO
!!$WRITE(*,FMT='("ISPA,ISPB",T20,10I5)')ISPA,ISPB
!!$WRITE(*,FMT='("LXA,LXB",T20,10I5)')LXA,LXB
!!$WRITE(*,FMT='("LMXA,LMXB",T20,10I5)')LMXA,LMXB
!!$WRITE(*,FMT='("LNXA,LNXB",T20,10I5)')LNXA,LNXB
!!$WRITE(*,FMT='("LOXA",T20,20I5)')LOXA
!!$WRITE(*,FMT='("LOXB",T20,20I5)')LOXB
!!$WRITE(*,FMT='("LMN0A",T20,20I5)')LMN0A
!!$WRITE(*,FMT='("LMN0B",T20,20I5)')LMN0B
!!$WRITE(*,FMT='("INDX,IIX,",T20,10I10)')INDX,IIX
!!$WRITE(*,FMT='("LR1X,LR2X,MABSX",T20,10I5)')LR1X,LR2X,MABSX

!
!     ==========================================================================
!     == EXTRACT FOUR SLATER INTEGRALS                                        ==
!     ==========================================================================
      ALLOCATE(INTEGRAL(INDX))
      ALLOCATE(DINTEGRAL(INDX))
      IND=0
      DO LN1=1,LNXA
        L1=LOXA(LN1)
        DO LN2=LN1,LNXA
          L2=LOXA(LN2)
          DO LN3=1,LNXA
            L3=LOXA(LN3)
            DO LN4=1,LNXB
              L4=LOXB(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX),2
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
                    IND=IND+1
!                   == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE...
                    CALL LMTO_OFFSITEXVALUE('31',ISPA,ISPB,IND,ABS(DIS) &
     &                                          ,INTEGRAL(IND),DINTEGRAL(IND))
                    IF(DIS.LT.0.D0) THEN
                      SVAR=REAL((-1)**(L1+L2+L3+L4))
                      INTEGRAL(IND)=SVAR*INTEGRAL(IND)
                      DINTEGRAL(IND)=-SVAR*DINTEGRAL(IND)
                    END IF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == EXTRACT PRODUCTS OF GAUNT COEFFICIENTS                               ==
!     ==========================================================================
      LM1=MAX(MAX(LMXA,(LR1X+1)**2),(LR2X+1)**2)
      ALLOCATE(CGMAT(LM1,LM1,LM1))
      CGMAT(:,:,:)=0.D0
      DO LM3=1,LMXA
        DO LM2=1,MAX((LR2X+1)**2,LMXA)
          DO LM1=1,(LR1X+1)**2
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG1)
            CGMAT(LM1,LM2,LM3)=CG1
            CGMAT(LM2,LM1,LM3)=CG1
            CGMAT(LM1,LM3,LM2)=CG1
            CGMAT(LM3,LM1,LM2)=CG1
            CGMAT(LM2,LM3,LM1)=CG1
            CGMAT(LM3,LM2,LM1)=CG1
          ENDDO
        ENDDO
      ENDDO

      ALLOCATE(LIST(IIX))
      ALLOCATE(IIA1(LXA+1,LXA+1,LXA+1,LXB+1,LR1X+1,LR2X+1,MABSX+1))
      ALLOCATE(IIA2(LXA+1,LXA+1,LXA+1,LXB+1,LR1X+1,LR2X+1,MABSX+1))
CALL TIMING$CLOCKON('X31-A')
      II=0
      DO L1=0,LXA
        DO L2=0,LXA
          DO L3=0,LXA
            DO L4=0,LXB
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX),2
                LMR1I=(LR1**2)+1
                LMR1F=(LR1+1)**2
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
!
                    IIA1(L1+1,L2+1,L3+1,L4+1,LR1+1,LR2+1,MABS+1)=II+1

                    DO PM=-1,1,2
                      IF(MABS.EQ.0.AND.PM.EQ.1) CYCLE
                      M=MABS*PM
                      M4=L4+1+M
                      LM4=L4**2+M4
                      LMR2=LR2**2+LR2+1+M
                      LM3=L3**2
                      DO M3=1,2*L3+1
                        LM3=LM3+1
                        LM2=L2**2
                        DO M2=1,2*L2+1
                          LM2=LM2+1
                          LM1=L1**2
                          DO M1=1,2*L1+1
                            LM1=LM1+1
                            SVAR=SUM(CGMAT(LMR1I:LMR1F,LM1,LM2) &
                   &                *CGMAT(LMR1I:LMR1F,LM3,LMR2))
                            IF(SVAR.EQ.0.D0) CYCLE
                            II=II+1
                            IF(II.GT.IIX) THEN
                              CALL ERROR$STOP('INDEX II EXCEEDS ARRAY SIZE')
                              CALL ERROR$STOP('LMTO_OFFSITEX31U_TEST1')
                            END IF
                            LIST(II)%M1=M1
                            LIST(II)%M2=M2
                            LIST(II)%M3=M3
                            LIST(II)%M4=M4
                            LIST(II)%VAL=SVAR
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                    IIA2(L1+1,L2+1,L3+1,L4+1,LR1+1,LR2+1,MABS+1)=II
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
CALL TIMING$CLOCKOFF('X31-A')
DEALLOCATE(CGMAT)
!
      ALLOCATE(IIBNDS(2,INDX))
      IND=0
      DO LN1=1,LNXA
        L1=LOXA(LN1)
        DO LN2=LN1,LNXA
          L2=LOXA(LN2)
          DO LN3=1,LNXA
            L3=LOXA(LN3)
            DO LN4=1,LNXB
              L4=LOXB(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX),2
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
                    IND=IND+1
                    IIBNDS(1,IND)=IIA1(L1+1,L2+1,L3+1,L4+1,LR1+1,LR2+1,MABS+1)
                    IIBNDS(2,IND)=IIA2(L1+1,L2+1,L3+1,L4+1,LR1+1,LR2+1,MABS+1)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(IIA1)
      DEALLOCATE(IIA2)
!
!     ==========================================================================
!     == CALCULATE FOUR CENTER INTEGRALS
!     ==========================================================================
CALL TIMING$CLOCKON('X31-B')
      U(:,:,:,:)=0.D0
      DU(:,:,:,:)=0.D0
      IND=0
      DO LN1=1,LNXA
        L1=LOXA(LN1)
        DO LN2=LN1,LNXA
          L2=LOXA(LN2)
          DO LN3=1,LNXA
            L3=LOXA(LN3)
            DO LN4=1,LNXB
              L4=LOXB(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LR1X),2
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
                    IND=IND+1
                    II1=IIBNDS(1,IND)
                    II2=IIBNDS(2,IND)
                    DO II=II1,II2
                      LMN1=LMN0A(LN1)+LIST(II)%M1
                      LMN2=LMN0A(LN2)+LIST(II)%M2
                      LMN3=LMN0A(LN3)+LIST(II)%M3
                      LMN4=LMN0B(LN4)+LIST(II)%M4
                      U(LMN1,LMN2,LMN3,LMN4)=U(LMN1,LMN2,LMN3,LMN4) &
     &                                      +LIST(II)%VAL*INTEGRAL(IND)
                      DU(LMN1,LMN2,LMN3,LMN4)=DU(LMN1,LMN2,LMN3,LMN4) &
     &                                      +LIST(II)%VAL*DINTEGRAL(IND)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
CALL TIMING$CLOCKOFF('X31-B')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OFFSITEXBONDU(ISP1,ISP2,DIS,LMNX1,LMNX2,U,DU)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS 
!     ** TWO ORBITALS ARE ON THE FIRST AND TWO ORBITALS ON THE SECOND ATOM    **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR=>POTPAR1
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
        IF(.NOT.ALLOCATED(OFFSITEX)) THEN
          CALL ERROR$MSG('ARRAY OFFSITEX IN LMTO_MODULE HAS NOT BEEN ALLOCATED')
          CALL ERROR$MSG('INTERNAL CODING ERROR')
          CALL ERROR$STOP('LMTO_OFFSITEXVALUE')
        END IF
        NF=OFFSITEX(ISP1,ISP2)%NF
        IF(NF.GT.LF) THEN
          CALL ERROR$MSG('INTERNAL ARRAY SIZE EXCEEDED')
          CALL ERROR$STOP('LMTO_OFFSITEXVALUE')
        END IF
        DISSAVE=DIS      
        G(:NF)=EXP(-OFFSITEX(ISP1,ISP2)%LAMBDA(:)*DIS)
      END IF
      IF(ID.EQ.'22') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%X22)) THEN
          CALL ERROR$MSG('NDDO TERM REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('LMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%X22(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%X22(:NF,IND)*G(:NF))
      ELSE IF (ID.EQ.'31') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%X31)) THEN
          CALL ERROR$MSG('3-1 EXCHANGE TERM REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('LMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%X31(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%X31(:NF,IND)*G(:NF))
      ELSE IF (ID.EQ.'BONDU') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%BONDU)) THEN
          CALL ERROR$MSG('U-TEM OF BOND OVERLAPS REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('LMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%BONDU(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%BONDU(:NF,IND)*G(:NF))
      ELSE IF (ID.EQ.'OV') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%OVERLAP)) THEN
          CALL ERROR$MSG('OFFSITE OVERLAP REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('LMTO_OFFSITEXVALUE')
        END IF
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
STOP 'FORCED IN LMTO_TWOCENTER_TEST'
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
!      ** DETERMINES THE TWO-CENTER INTEGRAL OF TWO FUNCTIONS SPECIFIED       **
!      ** ON A RADIAL GRID AND REAL SPHERICAL HARMONICS. THE TWO FUNCTIONS    **
!      ** ARE DISPLACED IN Z-DIRECTION                                        **
!      **                                                                     **
!      ** SEE Z. ROMANOWSKI, INT.J.QUANT.CHEM.108, 249 (2008),                **
!      **     Z. ROMANOWSKI, INT.J.QUANT.CHEM.108, 487 (2008) AND             **
!      **     Z. ROMANOWSKI AND A. F. JALBOT, J. MATH. CHEM. 46, 97 (2009)    **
!      *************************************************************************
       USE LMTO_TWOCENTER_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: L1_       ! MAIN ANGULAR MOMENTUM OF F1
       INTEGER(4),INTENT(IN) :: M1_       ! M-INDEX OF REAL SPH. HARMONICS OF F1
       INTEGER(4),INTENT(IN) :: GID1_     ! GRID ID OF F1 (SEE PAW_RADIAL.F90)
       INTEGER(4),INTENT(IN) :: NR1_      ! #(RADIAL GRID POINTS FOR F1)
       REAL(8)   ,INTENT(IN) :: F1_(NR1_) ! RADIAL PART OF F1
       INTEGER(4),INTENT(IN) :: L2_       ! MAIN ANGULAR MOMENTUM OF F2
       INTEGER(4),INTENT(IN) :: M2_       ! M-INDEX OF REAL SPH. HARMONICS OF F2
       INTEGER(4),INTENT(IN) :: GID2_     ! GRID ID OF F2 (SEE PAW_RADIAL.F90)
       INTEGER(4),INTENT(IN) :: NR2_      ! #(RADIAL GRID POINTS FOR F2)
       REAL(8)   ,INTENT(IN) :: F2_(NR2_) ! RADIAL PART OF F2
       REAL(8)   ,INTENT(IN) :: DIS_      ! DISTANCE IN Z-DIRECTION
       REAL(8)   ,INTENT(IN) :: TOLERANCE ! ALLOWED DEVIATION FROM EXACT RESULT
       REAL(8)   ,INTENT(OUT):: OVERLAP   ! RESULTING TWO-CENTER INTEGRAL
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
       REAL(8),PARAMETER  :: PI=4.D0*ATAN(1.D0)
       REAL(8),PARAMETER  :: FPI=4.D0*PI
       REAL(8),PARAMETER  :: SQ2=SQRT(2.D0)
       REAL(8)            :: COSTHETA1,COSTHETA2
       REAL(8)            :: PLM1,PLM2
       REAL(8)            :: R0(2),T1(2),T2(2),R(2)
       REAL(8)            :: R1,R2
       REAL(8)            :: DRDP
       REAL(8)            :: VAL1,VAL2
!      *************************************************************************
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
!      == SPHERICAL_PLGNDR RETURNS THE ASSOCIATED LEGENDRE POLYNOMIALS 
!      == MULTIPLIED WITH  SQRT( (L-M)! / (L+M)! ) FOR LM(L,M)=L(L+1)+1-M 
       CALL SPHERICAL_PLGNDR((L1+1)**2,L1,COSTHETA1,PLMWORK1)
       CALL SPHERICAL_PLGNDR((L2+1)**2,L2,COSTHETA2,PLMWORK2)
       PLM1=PLMWORK1(L1*(L1+1)+M1+1)
       PLM2=PLMWORK2(L2*(L2+1)+M2+1)
       PLM1=PLM1*SQRT(REAL(2*L1+1,KIND=8)/FPI)
       PLM2=PLM2*SQRT(REAL(2*L2+1,KIND=8)/FPI)
       IF(M1.NE.0)PLM1=PLM1*SQ2
       IF(M2.NE.0)PLM2=PLM2*SQ2  ! FIXED ERROR FROM PLM2=PLM1*SQ2 150128
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
       REAL(8)  ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
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
       TYPE(SEGMENT_TYPE),ALLOCATABLE :: STACK(:) !(LSTACKX)
       TYPE(SEGMENT_TYPE)     :: SEGMENT1
       TYPE(SEGMENT_TYPE)     :: SEGMENT2
       INTEGER(4)             :: NSEGMENTS
       REAL(8)                :: ERROR ! CURRENT ERROR ESTIMATE
       REAL(8)                :: ERRMAX,ERRMIN
       INTEGER(4)             :: IAXIS
       INTEGER(4)             :: I
!      *************************************************************************
       ALLOCATE(STACK(LSTACKX))
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
         IF(TPR.AND.MODULO(NSEGMENTS,1000).EQ.0) &
       &     PRINT*,'NEXT ',NSEGMENTS,VALUE,ERROR,STACK(1)%CENTER,STACK(1)%WIDTH
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
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: LMX
      INTEGER(4)            :: LMRX
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:)
      REAL(8)   ,ALLOCATABLE:: AUX1(:)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)  !(NR,LMR,IDIMD) DENSITY
      REAL(8)   ,ALLOCATABLE:: RHOWC(:,:,:)  !(NR,LMR,IDIMD) DENSITY W CORE
      REAL(8)   ,ALLOCATABLE:: POT(:,:,:)  !(NR,LMR,IDIMD) POTENTIAL
      REAL(8)               :: EH,ETOTC
      REAL(8)               :: SVAR
      REAL(8)               :: CG ! GAUNT COEFFICIENT
      INTEGER(4)            :: ISP,L
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
!     ==========================================================================
!     == CALCULATE TOTAL ENERGY                                               ==
!     ==========================================================================
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
      USE LMTO_MODULE, ONLY : ISPECIES,LNX,LOX &
     &                       ,DENMAT=>DENMAT_NEW &
     &                       ,SBAR=>SBAR_NEW     &
     &                       ,POTPAR=>POTPAR1
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
CALL ERROR$MSG('ROUTINE NOT YET ADJUSTED TO NEW DATA STRUCTURE')
CALL ERROR$MSG('DENMAT AND POTPAR NOW CONTAIN ONLY LOCAL ORBITALS')
CALL ERROR$STOP('LMTO_TESTDENMAT_1CDENMAT')
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
!     **  CALCULATE TRACE[DENMAT*OBERLAP] TO CHECK THE SUM-RULE               **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : TOFFSITE &
     &                       ,OVERLAP &
     &                       ,DENMAT=>DENMAT_NEW
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TON=.TRUE.
      INTEGER(4)            :: NND
      INTEGER(4)            :: NNO
      INTEGER(4)            :: IAT1,IAT2,IT(3)
      INTEGER(4)            :: N1,N2
      INTEGER(4)            :: NN,MM,I,J
      REAL(8)               :: SVAR1,SVAR2
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                           CALL TRACE$PUSH('LMTO_TESTDENMAT')
      IF(.NOT.TOFFSITE) THEN
        CALL ERROR$MSG('OVERLAP CAN ONLY BE CALCULATED WITH TOFFSITE=TRUE')
        CALL ERROR$STOP('LMTO_TESTDENMAT')
      END IF
      CALL LMTO_OVERLAPEVAL()
      WRITE(*,FMT='(80("="),T20,"  TR[DENMAT*OVERLAP] ")')
      WRITE(*,FMT='(80("="),T20,"  CHARGE CONTRIBUTION DQ FROM EACH PAIR  ")')
      NND=SIZE(DENMAT)
      NNO=SIZE(OVERLAP)
      SVAR1=0.D0
      DO NN=1,NND
        IAT1=DENMAT(NN)%IAT1
        IAT2=DENMAT(NN)%IAT2
        IT=DENMAT(NN)%IT
        DO MM=1,NNO
          IF(OVERLAP(MM)%IAT1.NE.IAT2) CYCLE
          IF(OVERLAP(MM)%IAT2.NE.IAT1) CYCLE
          IF(OVERLAP(MM)%IT(1).NE.-IT(1)) CYCLE
          IF(OVERLAP(MM)%IT(2).NE.-IT(2)) CYCLE
          IF(OVERLAP(MM)%IT(3).NE.-IT(3)) CYCLE
          N1=DENMAT(NN)%N1
          N2=DENMAT(NN)%N2
          SVAR2=0.D0
          DO J=1,N2
!!$            WRITE(*,FMT='("D ",200F10.3)')DENMAT(NN)%MAT(J,:,1)
!!$            WRITE(*,FMT='("O ",200F10.3)')OVERLAP(MM)%MAT(:,J)
            DO I=1,N1
              SVAR2=SVAR2+DENMAT(NN)%MAT(I,J,1)*OVERLAP(MM)%MAT(J,I)
            ENDDO
          ENDDO
          SVAR1=SVAR1+SVAR2
!         __SVAR2=NUMBER OF ELECTRONS FROM THIS PAIR____________________________
!         __SVAR1=SUM OF CHARGES OVER ALL PAIRS PROCESSED SO FAR________________
          WRITE(*,FMT='("IAT1=",I4," IAT2=",I4," IT=",3I3' &
     &                //'," DQ=",F10.5," SUMQ=",F10.5)')IAT1,IAT2,IT,SVAR2,SVAR1
        ENDDO
      ENDDO
      PRINT*,'TOTAL #PARTICLES ',SVAR1
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OVERLAPEVAL()
!     **************************************************************************
!     **  CALCULATEDS OFF-SITE OVERLAP MATRIX "OVERLAP"                       **
!     **  OVERLAP USES THE SAME DATA STRUCTURE AS DENMAT AND INHERITS         **
!     **  THE SAME NEIGHBORLIST.                                              **
!     **                                                                      **
!     **  CAUTION! ROUTINE IS USED ALSO BY PAW_DMFT OBJECT!                   **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES &
     &                       ,POTPAR=>POTPAR1 &
     &                       ,SBAR=>SBAR_NEW &
     &                       ,NSP &
     &                       ,OFFSITEX &
     &                       ,OVERLAP &
     &                       ,TCTE
      IMPLICIT NONE
      REAL(8)                :: RBAS(3,3)
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4)             :: NNS
      INTEGER(4)             :: NN
      INTEGER(4)             :: ISPA,ISPB
      INTEGER(4)             :: LMNXA,LMNXB
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
      INTEGER(4),ALLOCATABLE :: INS(:)
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_OVERLAPEVAL')
!
!     ==========================================================================
!     == CHECKS
!     ==========================================================================
      IF(.NOT.ALLOCATED(OFFSITEX)) THEN
!       -- THE ARRAY OFFSITEX IS USED IN THE ROUTINE LMTO_OFFSITEOVERLAP BELOW--
        CALL ERROR$MSG('ARRAY OFFSITEX NOT ALLOCATED')
        CALL ERROR$MSG('LMTO_OFFXINT MUST BE CALLED TO SET IT UP')
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
!
!     == DEALLOCATE FIRST ======================================================
      IF(ALLOCATED(OVERLAP)) THEN
        NNS=SIZE(OVERLAP)
        DO NN=1,NNS
          DEALLOCATE(OVERLAP(NN)%MAT)
        ENDDO
        DEALLOCATE(OVERLAP)
      END IF

!     == ALLOCATE OVERLAP MATRIX NOW ===========================================
      NNS=SIZE(SBAR)
      ALLOCATE(OVERLAP(NNS))
      DO NN=1,NNS
        IATA=SBAR(NN)%IAT1
        IATB=SBAR(NN)%IAT2
        ISPA=ISPECIES(IATA)
        ISPB=ISPECIES(IATB)
        NA=SUM(2*POTPAR(ISPA)%LOFH+1)
        NB=SUM(2*POTPAR(ISPB)%LOFH+1)
        OVERLAP(NN)%IAT1=IATA
        OVERLAP(NN)%IAT2=IATB
        OVERLAP(NN)%IT=SBAR(NN)%IT(:)
        OVERLAP(NN)%N1=NA
        OVERLAP(NN)%N2=NB
        ALLOCATE(OVERLAP(NN)%MAT(NA,NB))
        OVERLAP(NN)%MAT(:,:)=0.D0
      ENDDO
!
!     ==========================================================================
!     == COLLECT POINTERS INS(IAT) TO ONSITE STRUCTURE CONSTANTS              ==
!     ==========================================================================
      ALLOCATE(INS(NAT))
      INS=0
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IATA=SBAR(NN)%IAT1
        IATB=SBAR(NN)%IAT2
        IF(IATA.EQ.IATB.AND.SUM(ABS(SBAR(NN)%IT)).EQ.0) INS(IATA)=NN
      ENDDO
!
      DO IAT=1,NAT
        IF(INS(IAT).LE.0) THEN
          CALL ERROR$MSG('INDEX ARRAY FOR ONSITE STRUCTURE CONSTANTS')
          CALL ERROR$MSG('NO ONSITE TERMS FOUND FOR ATOM')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$STOP('LMTO_OVERLAPEVAL')
        END IF
      ENDDO
!
!     ==========================================================================
!     == LOOP OVER PAIRS                                                      ==
!     ==========================================================================
      DO NN=1,NNS
        IATA=OVERLAP(NN)%IAT1
        IATB=OVERLAP(NN)%IAT2
        ISPA=ISPECIES(IATA)
        ISPB=ISPECIES(IATB)
        RA(:)=R0(:,IATA)
        RB(:)=R0(:,IATB)+MATMUL(RBAS,REAL(OVERLAP(NN)%IT(:),KIND=8))
        LMNXA=POTPAR(ISPA)%TAILED%LMNX
        LMNXB=POTPAR(ISPB)%TAILED%LMNX
        NA=OVERLAP(NN)%N1
        NB=OVERLAP(NN)%N2
!
!       == ONSITE ELEMENT=======================================================
        IF(IATA.EQ.IATB.AND.SUM(OVERLAP(NN)%IT**2).EQ.0)  THEN
          ALLOCATE(OV(LMNXA,LMNXB))
          OV=POTPAR(ISPA)%TAILED%OVERLAP
          IF(TCTE) THEN
            CALL LMTO_EXPANDLOCALWITHCTE('BACK',IATA,1,NA,LMNXA &
     &                          ,OVERLAP(NN)%MAT,OV)
          ELSE
            CALL LMTO_EXPANDLOCAL('BACK',1,NA,LMNXA,SBAR(INS(IATA))%MAT &
     &                          ,OVERLAP(NN)%MAT,OV)
          END IF
          DEALLOCATE(OV)
          CYCLE
        END IF
!
!       == CONSIDER ONLY OFFSITE TERMS =========================================
!
!       ========================================================================
!       == CALCULATE OVERLAP MATRIX WITH Z ATONG THE BOND AXIS                ==
!       ========================================================================
        ALLOCATE(OV(LMNXA,LMNXB))
        ALLOCATE(DOV(LMNXA,LMNXB))
        DIS=SQRT(SUM((RB-RA)**2))
        CALL LMTO_OFFSITEOVERLAP(ISPA,ISPB,DIS,LMNXA,LMNXB,OV,DOV)
!!$WRITE(*,FMT='(82("="),T20,"  OVERLAP IN Z DIRECTION ",2I2)')ISPA,ISPB
!!$DO I=1,LMNXA
!!$  WRITE(*,FMT='(200F10.5)')OV(I,:)
!!$ENDDO
!
!       ========================================================================
!       == ROTATE DENSITY MATRIX SO THAT DISTANCE VECTOR POINTS IN Z-DIRECTION=
!       ========================================================================
!       == CONSTRUCT ROTATION MATRIX ===========================================
!       == DISTANCE VECTOR WILL BE NEW Z-DIRECTION =============================
        DR(:)=RB(:)-RA(:)
        DIS=SQRT(SUM(DR**2))
        ROT(:,3)=DR(:)/DIS
!       == FIRST VECTOR IS VECTOR PRODUCT OF THE MOST ORTHOGONAL UNIT VECTOR ===
        I=MINLOC(ABS(DR),1)
        ROT(:,2)=0.D0
        ROT(I,2)=1.D0
        ROT(1,1)=ROT(2,2)*ROT(3,3)-ROT(3,2)*ROT(2,3)
        ROT(2,1)=ROT(3,2)*ROT(1,3)-ROT(1,2)*ROT(3,3)
        ROT(3,1)=ROT(1,2)*ROT(2,3)-ROT(2,2)*ROT(1,3)
        ROT(:,1)=ROT(:,1)/SQRT(SUM(ROT(:,1)**2))
!       == SECOND VECTOR AS VECTOR PRODUCT BETWEEN FIRST AND THIRD VECTOR ======
        ROT(1,2)=ROT(2,3)*ROT(3,1)-ROT(3,3)*ROT(2,1)
        ROT(2,2)=ROT(3,3)*ROT(1,1)-ROT(1,3)*ROT(3,1)
        ROT(3,2)=ROT(1,3)*ROT(2,1)-ROT(2,3)*ROT(1,1)
!       == REMOVE INVERSION, IF PRESENT ========================================
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
        IF(TCTE) THEN
          CALL LMTO_EXPANDNONLOCALWITHCTE('BACK',IATA,IATB,1,NA,NB,LMNXA,LMNXB &
     &                          ,OVERLAP(NN)%MAT,OV)
        ELSE
          CALL LMTO_EXPANDNONLOCAL('BACK',1,NA,NB,LMNXA,LMNXB &
     &                            ,SBAR(INS(IATA))%MAT,SBAR(INS(IATB))%MAT &
     &                            ,OVERLAP(NN)%MAT,OV)
        END IF
        DEALLOCATE(OV)
        DEALLOCATE(DOV)
      ENDDO
!
!     ==========================================================================
!     == CLEAN UP TO AVOID MEMORY LEAK                                        ==
!     ==========================================================================
!!$CALL LMTO$REPORTOVERLAP(6)
!!$STOP 'FORCED'
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$REPORTOVERLAP(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : OVERLAP,DENMAT=>DENMAT_NEW
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
      USE LMTO_MODULE, ONLY : SBAR=>SBAR_NEW
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
      SUBROUTINE LMTO$REPORTPERIODICMAT(NFIL,NAME,NNS,SBAR)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE,ONLY: PERIODICMAT_TYPE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: NFIL
      CHARACTER(*)          ,INTENT(IN) :: NAME
      INTEGER(4)            ,INTENT(IN) :: NNS
      TYPE(PERIODICMAT_TYPE),INTENT(IN) :: SBAR(NNS)
      INTEGER(4)                        :: NN,LM1
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTPERIODICMAT')
      WRITE(NFIL,FMT='(80("="))')
      WRITE(NFIL,FMT='(80("="),T10,"  ",A,"   ")')TRIM(NAME)
      WRITE(NFIL,FMT='(80("="))')
      DO NN=1,NNS
        IF(SBAR(NN)%N1*SBAR(NN)%N2.EQ.0) CYCLE
        WRITE(NFIL,FMT='(80("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 SBAR(NN)%IAT1,SBAR(NN)%IAT2,SBAR(NN)%IT
        DO LM1=1,SBAR(NN)%N1
          WRITE(NFIL,FMT='(20F10.5)')SBAR(NN)%MAT(LM1,:)
        ENDDO
        WRITE(NFIL,FMT='(80("="))')
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END SUBROUTINE LMTO$REPORTPERIODICMAT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$REPORTDENMAT(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : DENMAT=>DENMAT_NEW
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NND
      INTEGER(4)            :: NN,LM1,NDIMD
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTDENMAT')
      NND=SIZE(DENMAT)
      WRITE(NFIL,FMT='(82("="))')
      WRITE(NFIL,FMT='(82("="),T10,"   DENSITY MATRIX   ")')
      WRITE(NFIL,FMT='(82("="))')
      DO NN=1,NND
        IF(DENMAT(NN)%N1*DENMAT(NN)%N2.EQ.0) CYCLE
        WRITE(NFIL,FMT='(82("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 DENMAT(NN)%IAT1,DENMAT(NN)%IAT2,DENMAT(NN)%IT
        DO LM1=1,DENMAT(NN)%N1
          WRITE(NFIL,FMT='(20F10.5)')DENMAT(NN)%MAT(LM1,:,1)
        ENDDO
        NDIMD=DENMAT(NN)%N3
        IF(NDIMD.GT.1) THEN
          WRITE(NFIL,FMT='(82("-"),T10," SPIN SZ CONTRIBTION ")')
          DO LM1=1,DENMAT(NN)%N1
            WRITE(NFIL,FMT='(20F10.5)')DENMAT(NN)%MAT(LM1,:,NDIMD)
          ENDDO
        END IF
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
      USE LMTO_MODULE, ONLY : HAMIL=>HAMIL_NEW
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
      USE LMTO_MODULE, ONLY : HAMIL=>HAMIL_NEW
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
        WRITE(100,FMT='(F15.10,2X,100(F20.10,2X))')R(IR),PHI(IR,:)
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
      SUBROUTINE LMTO_PLOTORBS()
!     **************************************************************************
!     ** PLOT TIGHT-BINDING ORBITALS                                          **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES &
     &                        ,POTPAR=>POTPAR1
      IMPLICIT NONE
!      LOGICAL(4),PARAMETER :: TON=.TRUE.
      LOGICAL(4),PARAMETER :: TON=.FALSE.
      LOGICAL(4),SAVE      :: TINI=.FALSE. !EXECUTE ONLY ONCE
      INTEGER(4) :: IAT,ISP,IORB,IH,L,IM
!     **************************************************************************
      IF(.NOT.TON) RETURN
      IF(TINI) RETURN  !EXECUTE ONLY ONCE
      TINI=.TRUE.
!
!     ==========================================================================
!     == CALCULATE ALL ORBITALS INTO CUBE FILES                               ==
!     ==========================================================================
      PRINT*,'BEFORE LMTO_PLOTTAILED'
      CALL LMTO_PLOTTAILED()
!
!     ==========================================================================
!     == CALCULATE ALL ORBITALS INTO CUBE FILES                               ==
!     ==========================================================================
      PRINT*,'BEFORE LMTO_GRIDPLOT'
      DO IAT=1,SIZE(ISPECIES)
        ISP=ISPECIES(IAT)
        IORB=0
        DO IH=1,POTPAR(ISP)%NHEAD
          L=POTPAR(ISP)%LOFH(IH)
          DO IM=1,2*L+1
            IORB=IORB+1
            CALL LMTO_PLOTNTBO(IAT,IORB,'1C','1D')
            CALL LMTO_PLOTNTBO(IAT,IORB,'MC','1D')
          ENDDO
        ENDDO
      ENDDO
      PRINT*,'PLOTORBS FINISHED'
      CALL ERROR$MSG('PLANNED EXIT AFTER PLOTTING FOR ANALYSIS')
      CALL ERROR$STOP('LMTO_PLOTORBS')
      RETURN
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
      USE LMTO_MODULE, ONLY: ISPECIES &
     &                      ,POTPAR=>POTPAR1
      IMPLICIT NONE
      INTEGER(4)             :: NAT  !#(ATOMS)
      INTEGER(4)             :: LX   !X(ANGULAR MOMENTUM)
      INTEGER(4)             :: LMX  !#(ANGULAR MOMENTA FOR THE ORBITAL)
      INTEGER(4)             :: GID  !GRID ID FOR RADIAL GRID
      INTEGER(4)             :: NR   !#(RADIAL GRID POINTS)
      REAL(8)   ,ALLOCATABLE :: ORB(:,:)  !(NR,LMX)
      INTEGER(4)             :: ISP ! ATOM-TYPE INDEX
      INTEGER(4)             :: NHEAD ! #(RADIAL HEAD FUNCTIONS)
      INTEGER(4)             :: IAT
      INTEGER(4)             :: IORB  ! ORBITAL INDEX
      INTEGER(4)             :: IH,LH,MH
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
        NHEAD=POTPAR(ISP)%NHEAD
        LX=MAXVAL(POTPAR(ISP)%TAILED%LOX)
        LMX=(LX+1)**2
        ALLOCATE(ORB(NR,LMX))
        IORB=0
        DO IH=1,NHEAD
          LH=POTPAR(ISP)%LOFH(IH)
          DO MH=1,2*LH+1
            IORB=IORB+1
            CALL LMTO_TAILEDORBLM(IAT,IORB,NR,LMX,ORB)
            WRITE(CHIAT,FMT='(I5)')IAT
            WRITE(CHORB,FMT='(I5)')IORB
            STRING='CHI_FORATOM'//TRIM(ADJUSTL(CHIAT)) &
     &                          //'_'//TRIM(ADJUSTL(CHORB))//'.DAT'
            CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LMX,ORB)
          ENDDO !END OF LOOP OVER MH (MAGNETIC QUANTUM NUMBER OF HEAD FUNCTION)
        ENDDO !END OF LOOP OVER HEAD FUNCTIONS (WITHOUT M-MULTIPLICITY)
        DEALLOCATE(ORB)
      ENDDO  !END OF LOOP OVER ATOMS
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PLOTNTBO(IAT0,IORB,TYPE,ID)
!     **************************************************************************
!     **  WRITES THE ORBITAL IORB ON ATOM IAT0 TO FILE.                       **
!     **                                                                      **
!     **  SELECT TAILED REPRESENTATION USING TYPE='1C' (ONE-CENTER)           **
!     **  OR THE MULTICENTER EXPANSION USING BARE HANKEL FUNCTION USING 'MC'  **
!     **                                                                      **
!     **  SELECT A CUBE FILE FOR SURFACE PLOTS WITH ID='3D'                   **
!     **                                                                      **
!     **  SET N1,N2,N3 EQUAL TO THE NUMBER OF GRIDPOINTS IN EACH DIRECTION    **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2009-2014********
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      USE LMTO_MODULE, ONLY : ISPECIES &
     &                       ,SBAR=>SBAR_NEW &
     &                       ,POTPAR=>POTPAR1
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: IAT0
      INTEGER(4)  ,INTENT(IN) :: IORB
      CHARACTER(*),INTENT(IN) :: TYPE  ! '1C' OR 'MC'
      CHARACTER(*),INTENT(IN) :: ID    ! '1D', '2D', OR '3D'
      INTEGER(4),PARAMETER    :: N1=40,N2=40,N3=40 !GRID (1D?)
      INTEGER(4),PARAMETER    :: NRAD=200
      INTEGER(4),PARAMETER    :: NDIRX=100
      REAL(8)   ,PARAMETER    :: RANGE=8.D0
      CHARACTER(32),PARAMETER :: STARTYPE='CUBEAXES'
!      CHARACTER(32),PARAMETER :: STARTYPE='NEIGHBORS'
      REAL(8)               :: DIR(3,NDIRX)
      REAL(8)               :: ORIGIN(3)
      REAL(8)               :: TVEC(3,3)    !BOX
      REAL(8)               :: RBAS(3,3)
      INTEGER(4)            :: NAT
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
      REAL(8)   ,ALLOCATABLE:: ORB(:)
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
!     **************************************************************************
                                         CALL TRACE$PUSH('LMTO_PLOTNTBO')
!
!     ==========================================================================
!     == COLLECT INFORMATION                                                  ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      NAT=SIZE(ISPECIES)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     == DEFINE ATOMS ON THE CLUSTER OF NEIGHBORS                             ==
!     == NEEDED FOR BALL STICK MODEL IN CUBE FILE                             ==
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
        IF(STARTYPE.EQ.'NEIGHBORS') THEN
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
        ELSE IF(STARTYPE.EQ.'CUBEAXES') THEN
          NDIR=13
!         == (100) DIRECTIONS ==================================================
          DIR(:,1)=(/1.D0,0.D0,0.D0/)
          DIR(:,2)=(/0.D0,1.D0,0.D0/)
          DIR(:,3)=(/0.D0,0.D0,1.D0/)
!         == (110) DIRECTIONS ==================================================
          DIR(:,4)=(/0.D0,1.D0,1.D0/)
          DIR(:,5)=(/1.D0,0.D0,1.D0/)
          DIR(:,6)=(/1.D0,1.D0,0.D0/)
          DIR(:,7)=(/0.D0,-1.D0,1.D0/)
          DIR(:,8)=(/1.D0,0.D0,-1.D0/)
          DIR(:,9)=(/1.D0,-1.D0,0.D0/)
!         == (111) DIRECTIONS ==================================================
          DIR(:,10)=(/1.D0,1.D0,1.D0/)
          DIR(:,11)=(/-1.D0,1.D0,1.D0/)
          DIR(:,12)=(/1.D0,-1.D0,1.D0/)
          DIR(:,13)=(/1.D0,1.D0,-1.D0/)
        ELSE
          CALL ERROR$MSG('STARTYPE NOT RECOGNIZED')
          CALL ERROR$MSG('MAY BE EITHER "NEIGHBORS" OR "CUBEAXES"')
          CALL ERROR$CHVAL('STARTYPE',STARTYPE)
          CALL ERROR$STOP('LMTO_PLOTNTBO')
        END IF
        NP=NRAD*NDIR
        ALLOCATE(P(3,NP))
        CALL LMTO_GRIDORB_STARGRID(R0(:,IAT0),RANGE,NDIR,DIR,NRAD,X1D,P)
      ELSE
        CALL ERROR$MSG('INVALID ID')
        CALL ERROR$MSG('MUST BE "1D", "2D", OR "3D"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO_PLOTNTBO')
      END IF
!
!     ==========================================================================
!     == DETERMINE ORBITAL ON THE SPECIFIED GRID POINTS                       ==
!     ==========================================================================
      ISP=ISPECIES(IAT0)
      ALLOCATE(ORB(NP))
      IF(TYPE.EQ.'1C') THEN
        CALL LMTO_TAILED_NTBOOFR(IAT0,IORB,NP,P,ORB)
      ELSE IF(TYPE.EQ.'MC') THEN
        CALL LMTO_NTBOOFR(IAT0,IORB,NP,P,ORB)
      ELSE
        CALL ERROR$MSG('TYPE NOT RECOGNIZED')
        CALL ERROR$CHVAL('TYPE',TYPE)
        CALL ERROR$STOP('LMTO_PLOTNTBO')
      END IF
!
!     ==========================================================================
!     == CONSTRUCT FILE NAME                                                  ==
!     ==========================================================================
      FILE='NTBO'
!
!     == ATOM AND ORBITAL ======================================================
      WRITE(STRING,*)IAT0 
      FILE=TRIM(ADJUSTL(FILE))//'_AT'//TRIM(ADJUSTL(STRING))
      WRITE(STRING,*)IORB 
      FILE=TRIM(ADJUSTL(FILE))//'_ORB'//TRIM(ADJUSTL(STRING))
!
!     ==  TAILED OR MULTICENTER EXPANSION ======================================
      IF(TYPE.EQ.'MC') THEN
        FILE=TRIM(ADJUSTL(FILE))//'_'//'MC'  !MULTI CENTER
      ELSE IF(TYPE.EQ.'1C') THEN
        FILE=TRIM(ADJUSTL(FILE))//'_'//'1C'  !1-CENTER (TAILED REPRESENTATION)
      ELSE
        CALL ERROR$MSG('TYPE NOT RECOGNIZED')
        CALL ERROR$CHVAL('TYPE',TYPE)
        CALL ERROR$STOP('LMTO_PLOTNTBO')
      END IF
!
!     == DIMENSION OF PLOT =====================================================
      IF(ID.EQ.'3D') THEN
        FILE=TRIM(ADJUSTL(FILE))//'_'//'3D.CUB'
      ELSE IF(ID.EQ.'2D') THEN
        FILE=TRIM(ADJUSTL(FILE))//'_'//'2D.DAT'
      ELSE IF(ID.EQ.'1D') THEN
        FILE=TRIM(ADJUSTL(FILE))//'_'//'1D.DAT'
      ELSE 
        CALL ERROR$MSG('INVALID ID')
        CALL ERROR$MSG('MUST BE "1D", "2D", OR "3D"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO_PLOTNTBO')
      END IF
!
!     ==========================================================================
!     == WRITE ORBITALS TO FILE                                               ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
      CALL FILEHANDLER$UNIT('HOOK',NFIL)
      IF(ID.EQ.'3D') THEN
        CALL LMTO_WRITECUBEFILE(NFIL,NATCLUSTER,ZCLUSTER,RCLUSTER &
     &                         ,ORIGIN,TVEC,N1,N2,N3,ORB)
      ELSE IF(ID.EQ.'1D') THEN
        DO I=1,NRAD
          WRITE(NFIL,*)X1D(I),(ORB(NRAD*(IDIR-1)+I),IDIR=1,NDIR)
        ENDDO
      ELSE IF(ID.EQ.'2D') THEN
        DO I=1,NP
          WRITE(NFIL,*)P(1:2,I),ORB(I)
        ENDDO
      ELSE
        CALL ERROR$MSG('INVALID ID')
        CALL ERROR$MSG('MUST BE "1D", "2D", OR "3D"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO_PLOTNTBO')
      END IF
      CALL FILEHANDLER$CLOSE('HOOK')
      CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOK')
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
      SUBROUTINE LMTO_TAILEDORBLM(IAT,IORB,NR,LMX,ORB)
!     **************************************************************************
!     ** CALCULATE SPHERICAL HARMONICS CONTRIBUTION OF THE THE TAILED ORBITAL **
!     ** NUMBER IORB ON ATOM IAT                                              **
!     **                                                                      **
!     ** CHI(R)=SUM_LM=1^LMX ORB(R,LM)*YLM(R)                                 **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY: ISPECIES &
     &                      ,POTPAR=>POTPAR1 &
     &                      ,SBAR=>SBAR_NEW &
     &                      ,CTE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: IORB
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMX
      REAL(8)   ,INTENT(OUT):: ORB(NR,LMX)
      INTEGER(4)            :: ISP      
      INTEGER(4)            :: LMN,LN,L,LM,M,LX
!     **************************************************************************
      ISP=ISPECIES(IAT)
      IF(SIZE(ORB(:,1)).NE.SIZE(POTPAR(ISP)%TAILED%AEF(:,1))) THEN
        CALL ERROR$MSG('INCONSISTENT ARGUMENT NR')
        CALL ERROR$STOP('LMTO_TAILEDORBLM')
      END IF
      LX=MAXVAL(POTPAR(ISP)%TAILED%LOX)
      IF(LMX.LT.(LX+1)**2) THEN
        CALL ERROR$MSG('DIMENSION LMX SMALLER THAN REQUIRED')
        CALL ERROR$I4VAL('LMX ON INPUT',LMX)
        CALL ERROR$I4VAL('LX ON POTPAR',LX)
        CALL ERROR$STOP('LMTO_TAILEDORBLM')
      END IF

      ORB=0.D0
      LMN=0    
      DO LN=1,POTPAR(ISP)%TAILED%LNX
        L=POTPAR(ISP)%TAILED%LOX(LN)
        DO M=1,2*L+1
          LMN=LMN+1
          LM=L**2+M
          ORB(:,LM)=ORB(:,LM) &
     &             +POTPAR(ISP)%TAILED%AEF(:,LN)*CTE(IAT)%MAT(LMN,IORB)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILEDORBLM_OLD(IAT,IORB,NR,LMX,ORB)
!     **************************************************************************
!     ** CALCULATE SPHERICAL HARMONICS CONTRIBUTION OF THE THE TAILED ORBITAL **
!     ** NUMBER IORB ON ATOM IAT                                              **
!     **                                                                      **
!     ** CHI(R)=SUM_LM=1^LMX ORB(R,LM)*YLM(R)                                 **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY: ISPECIES &
     &                      ,POTPAR=>POTPAR1 &
     &                      ,SBAR=>SBAR_NEW
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: IORB
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMX
      REAL(8)   ,INTENT(OUT):: ORB(NR,LMX)
      INTEGER(4)            :: ISP      
      INTEGER(4)            :: LX
      INTEGER(4)            :: NNS
      INTEGER(4)            :: INS
      INTEGER(4)            :: NN
      INTEGER(4)            :: NHEAD,NTAIL
      INTEGER(4)            :: I,L,M,LM,LMN
!     **************************************************************************
      ISP=ISPECIES(IAT)
      IF(SIZE(ORB(:,1)).NE.SIZE(POTPAR(ISP)%TAILED%AEF(:,1))) THEN
        CALL ERROR$MSG('INCONSISTENT ARGUMENT NR')
        CALL ERROR$STOP('LMTO_TAILEDORBLM')
      END IF
      LX=MAX(MAXVAL(POTPAR(ISP)%LOFH),MAXVAL(POTPAR(ISP)%LOFT))
      IF(LMX.LT.(LX+1)**2) THEN
        CALL ERROR$MSG('DIMENSION LMX SMALLER THAN REQUIRED')
        CALL ERROR$I4VAL('LMX ON INPUT',LMX)
        CALL ERROR$I4VAL('LX ON POTPAR',LX)
        CALL ERROR$STOP('LMTO_TAILEDORBLM')
      END IF
!
!     ==========================================================================
!     == DETERMINE INDEX OF ONSITE STRUCTURE CONSTANTS                        ==
!     ==========================================================================
      NNS=SIZE(SBAR)
      INS=0
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
        IF(SBAR(NN)%IAT2.NE.IAT) CYCLE
        IF(SUM(SBAR(NN)%IT(:)**2).NE.0) CYCLE
        INS=NN
        EXIT
      ENDDO
!
!     ==========================================================================
!     == COLLECT HEAD CONTRIBUTION                                            ==
!     ==========================================================================
      ORB(:,:)=0.D0
      NHEAD=POTPAR(ISP)%NHEAD
      LMN=0
      DO I=1,NHEAD
        L=POTPAR(ISP)%LOFH(I)
        DO M=1,2*L+1
          LMN=LMN+1
          IF(LMN.NE.IORB) CYCLE
          LM=L**2+M
          ORB(:,LM)=POTPAR(ISP)%TAILED%AEF(:,I)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  ADD TAIL CONTRIBUTIONS                                              ==
!     ==========================================================================
      NTAIL=POTPAR(ISP)%NTAIL
      LMN=0
      DO I=1,NTAIL
        L=POTPAR(ISP)%LOFT(I)
        DO M=1,2*L+1
          LMN=LMN+1
          LM=L**2+M
          ORB(:,LM)=ORB(:,LM)-POTPAR(ISP)%TAILED%AEF(:,NHEAD+I) &
     &                       *SBAR(INS)%MAT(IORB,LMN)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TAILED_NTBOOFR(IAT,IORB,NP,P,CHI)
!     **************************************************************************
!     **  MAPS THE SCREENED ORBITAL LMNORB AT ATOM IATORB                     **
!     **  ONTO THE GRID P                                                     **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2012 ************
      USE LMTO_MODULE, ONLY : ISPECIES &
     &                       ,POTPAR=>POTPAR1
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: IAT     ! INDEX OF CENTRAL ATOM
      INTEGER(4)  ,INTENT(IN) :: IORB    ! ORBITAL INDEX
      INTEGER(4)  ,INTENT(IN) :: NP      ! #(GRID POINTS)
      REAL(8)     ,INTENT(OUT):: CHI(NP) ! SCREENED ENVELOPE FUNCTIONS
      REAL(8)     ,INTENT(IN) :: P(3,NP) ! GRIDPOINT POSITIONS
      INTEGER(4)              :: NAT     ! #(ATOMS)
      REAL(8)   ,ALLOCATABLE  :: R0(:,:) ! ATOMIC POSITIONS
      INTEGER(4)              :: GID     ! GRID ID
      INTEGER(4)              :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4)              :: LX,LMX
      REAL(8)                 :: DR(3) ! DISTANCE VECTOR OF GRID POINT TO ATOM
      REAL(8)                 :: DIS   ! DISTANCE OF GRID POINT TO ATOM
      REAL(8)  ,ALLOCATABLE   :: YLM(:)    ! SPHERICAL HARMONICS
      REAL(8)  ,ALLOCATABLE   :: ORB(:,:)  ! RADIAL PARTS OF THE ORBITAL
      REAL(8)                 :: VAL
      INTEGER(4)              :: IP,LM
      INTEGER(4)              :: ISP
!     **************************************************************************
      ISP=ISPECIES(IAT)
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     == CALCULATE SPHERICAL HARMONICS REPRESENTATION OF THE ORBITAL          ==
!     ==========================================================================
      LX=MAXVAL(POTPAR(ISP)%TAILED%LOX)
      LMX=(LX+1)**2
      GID=POTPAR(ISP)%TAILED%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(ORB(NR,LMX))
      CALL LMTO_TAILEDORBLM(IAT,IORB,NR,LMX,ORB)
!
!     ==========================================================================
!     == LOOP OVER REAL SPACE GRID                                            ==
!     ==========================================================================
      ALLOCATE(YLM(LMX))
      CHI(:)=0.D0
      DO IP=1,NP
        DR(:)=P(:,IP)-R0(:,IAT)
        CALL SPHERICAL$YLM(LMX,DR,YLM)
        DIS=SQRT(SUM(DR**2))
        DO LM=1,LMX
          CALL RADIAL$VALUE(GID,NR,ORB(:,LM),DIS,VAL)
          CHI(IP)=CHI(IP)+VAL*YLM(LM)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBOOFR(IATORB,LMNORB,NP,P,ORB)
!     **************************************************************************
!     **  MAPS THE SCREENED ORBITAL LMNORB AT ATOM IATORB                     **
!     **  ONTO THE GRID P                                                     **
!     **                                                                      **
!     **  THE ONSITE TERMS ARE AUGMENTED WITH THE ALL-ELECTRON PARTIAL WAVES. **
!     **  THE OFF-SITE TERMS ARE OAUGMENTED WITH NODELESS PARTIAL WAVES       **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2012 ************
      USE LMTO_MODULE, ONLY : ISPECIES,NSP,K2,LNX &
     &                       ,SBAR=>SBAR_NEW &
     &                       ,POTPAR=>POTPAR1
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: IATORB  ! INDEX OF CENTRAL ATOM
      INTEGER(4)  ,INTENT(IN) :: LMNORB  ! ORBITAL INDEX
      INTEGER(4)  ,INTENT(IN) :: NP      ! #(GRID POINTS)
      REAL(8)     ,INTENT(OUT):: ORB(NP) ! SCREENED ENVELOPE FUNCTIONS
      REAL(8)     ,INTENT(IN) :: P(3,NP) ! GRIDPOINT POSITIONS
      REAL(8)               :: RBAS(3,3) ! LATTICE VECTORS
      INTEGER(4)            :: NAT       ! #(ATOMS)
      REAL(8)   ,ALLOCATABLE:: R0(:,:)   ! ATOMIC POSITIONS
      INTEGER(4)            :: IHORB     ! HEAD INDEX OF SELECTED ORBITAL
      INTEGER(4)            :: LMORB     ! MAJOR ANGULAR MOMENTUM OF ORBITAL
      INTEGER(4)            :: NNS 
      INTEGER(4)            :: GID       ! GRID ID
      INTEGER(4)            :: NR        ! #(RADIAL GRID POINTS)
      REAL(8)               :: R2(3) ! NEIGHBOR ATOM POSITION
      REAL(8)               :: DR(3) ! DISTANCE VECTOR OF GRID POINT TO ATOM
      REAL(8)               :: DIS   ! DISTANCE OF GRID POINT TO ATOM
      REAL(8)               :: VAL,VALDOT
      LOGICAL(4)            :: TONSITE   
      LOGICAL(4)            :: TSPHERE      ! INSIDE AN AUGMENTATION SPHERE?
      REAL(8)  ,ALLOCATABLE :: QBARVEC(:,:) !(LMNX,ISP)
      INTEGER(4),ALLOCATABLE:: LMOFLMN(:,:) !(LMNX,NSP) ANGULAR MOMENTUM INDEX
      REAL(8)  ,ALLOCATABLE :: SBARVEC(:)   !(LMX)
      REAL(8)  ,ALLOCATABLE :: K0(:)     ! BARE SOLID HANKEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: J0(:)     ! BARE SOLID BESSEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: YLM(:)    ! SPHERICAL HARMONICS
      REAL(8)  ,ALLOCATABLE :: R(:)      ! RADIAL GRID
      REAL(8)  ,ALLOCATABLE :: AEPHI(:,:)   ! ALL-ELECTRON PARTIAL WAVE
      REAL(8)  ,ALLOCATABLE :: AEPHIDOT(:,:)! ALL-ELECTRON SCATTERING WAVE
      REAL(8)  ,ALLOCATABLE :: NLPHI(:,:)   ! NODELESS PARTIAL WAVE
      REAL(8)  ,ALLOCATABLE :: NLPHIDOT(:,:)! NODELESS SCATTERING PARTIAL WAVE
      REAL(8)  ,ALLOCATABLE :: K0ARR(:)     ! BARE AUGMENTED HANKEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: DK0ARR(:)    ! 
      REAL(8)  ,ALLOCATABLE :: JBARARR(:,:) ! SCREENED AUGMENTED BESSEL FUNCTION
      REAL(8)  ,ALLOCATABLE :: DJBARARR(:,:)!
      INTEGER(4)            :: NHEAD,NTAIL ! #(HEAD/TAIL FUNCTIONS)
      INTEGER(4)            :: LMNX  ! X#(TAIL FUNCTIONS) OF ALL ATOM TYPES
      INTEGER(4)            :: LMX   ! X#(ANGULAR MOMENTA) OF ALL ATOM TYPES
      INTEGER(4)            :: NN,IP,ISP,ISP2,L,LM,IM,LN,LMN,IH,IT !LOOP INDICES
      INTEGER(4)            :: IAT2,L2X,LM2X,LMN2X,LM2,LMN2,LNDOT
!     **************************************************************************
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
!     == DETERMINE LM-, LN-, L- INDEX OF THE ORBITAL                    ==
!     ==========================================================================
      ISP=ISPECIES(IATORB)
      LMORB=0  
      LMN=0
      NHEAD=POTPAR(ISP)%NHEAD
      DO IH=1,NHEAD
        L=POTPAR(ISP)%LOFH(IH)
        DO IM=1,2*L+1
          LMN=LMN+1
          IF(LMN.EQ.LMNORB) THEN
            LMORB=L**2+IM   ! NEEDED FOR HANKEL FUNCTIONS
            IHORB=IH
          END IF
        ENDDO
      ENDDO
      IF(LMORB.EQ.0) THEN  
        CALL ERROR$MSG('ORBITAL NOT IDENTIFIED')
        CALL ERROR$I4VAL('IATORB',IATORB)
        CALL ERROR$I4VAL('LMNORB',LMNORB)
        CALL ERROR$STOP('LMTO_NTBOOFR')
      END IF
!
!     ==========================================================================
!     == DETERMINE QBAR OF THE SCATTERING CHANNEL                             ==
!     == THE ARRAY MATCHES THE SECOND (TAIL) INDEX OF SBAR                    ==
!     ==========================================================================
      LMX=0
      LMNX=0
      DO ISP=1,NSP
        LMNX=MAX(LMNX,SUM(2*POTPAR(ISP)%LOFT+1))
        LMX=MAX(LMX,(MAXVAL(POTPAR(ISP)%LOFT)+1)**2)
      ENDDO
      ALLOCATE(QBARVEC(LMNX,NSP))
      ALLOCATE(LMOFLMN(LMNX,NSP))
      QBARVEC(:,:)=0.D0
      LMOFLMN(:,:)=0
      DO ISP=1,NSP
        LMN=0
        DO IT=1,POTPAR(ISP)%NTAIL
          L=POTPAR(ISP)%LOFT(IT)
          DO IM=1,2*L+1
            LMN=LMN+1
            QBARVEC(LMN,ISP)=POTPAR(ISP)%QBAR(IT)
            LMOFLMN(LMN,ISP)=L**2+IM
          ENDDO
        ENDDO
      ENDDO
      ALLOCATE(SBARVEC(LMNX))
!
!     ==========================================================================
!     == CALCULATE INTERSTITIAL ORBITAL                                       ==
!     ==   |KBAR> - (|K0^SPHERE>-|JBAR^SPHERE>*SBART)                         ==
!     ==    = (|K0>-|K0^SPHERE>)*(1+QBAR*SBART) - |J0^SPHERE>SBART            ==
!     ==   WHERE ^SPHERE INDICATES THAT THE FUNCTION IS CUT OF OUTSIDE THE    ==
!     ==   AUGMENTATION SPHERE CONNECTED TO THE INDEX                         ==
!     ==========================================================================
      ALLOCATE(K0(LMX))
      ALLOCATE(J0(LMX))
      ORB(:)=0.D0
      NNS=SIZE(SBAR)  !SBAR STRUCCONS. (GLOBAL)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IATORB) CYCLE
        IAT2=SBAR(NN)%IAT2
        ISP2=ISPECIES(IAT2)
        R2(:)=R0(:,IAT2)+RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
     &                  +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
     &                  +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
        TONSITE=(IAT2.EQ.IATORB).AND.(MAXVAL(ABS(SBAR(NN)%IT(:))).EQ.0)
!       == SELECT STRUCTURE CONSTANT ARRAY =====================================
        LMN2X=SBAR(NN)%N2
        LM2X=(MAXVAL(POTPAR(ISP2)%LOFT)+1)**2
        SBARVEC(:LMN2X)=SBAR(NN)%MAT(LMNORB,:)
!
!       == LOOP OVER REAL SPACE GRID ===========================================
        DO IP=1,NP
!         == DR IS THE DISTANCE FROM THE SECOND ATOM ===========================
          DR(:)=P(:,IP)-R2(:)
          TSPHERE=(DOT_PRODUCT(DR,DR).LT.POTPAR(ISP2)%RAD**2)
          IF(TSPHERE) THEN
            CALL  LMTO$SOLIDBESSEL(DR,K2,LM2X,J0(1:LM2X))
!           == SUBTRACT THE BESSEL CONTRIBUTION OF THE SPHERE TERM
!           == THE HANKEL CONTRIBUTION HAS NOT BEEN ADDED IN THE FIRST PLACE
            DO LMN2=1,LMN2X
              LM2=LMOFLMN(LMN2,ISP2)
              ORB(IP)=ORB(IP)+J0(LM2)*SBARVEC(LMN2)
            ENDDO
          ELSE
!           == EXCLUDE SPHERE TERM OF HANKEL CONTRIBUTION ======================
            CALL LMTO$SOLIDHANKEL(DR,POTPAR(ISP2)%RAD,K2,LM2X,K0(1:LM2X))
            IF(TONSITE) ORB(IP)=ORB(IP)+K0(LMORB)
            DO LMN2=1,LMN2X
              LM2=LMOFLMN(LMN2,ISP2)
              ORB(IP)=ORB(IP)+K0(LM2)*QBARVEC(LMN2,ISP2)*SBARVEC(LMN2)
            ENDDO
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(K0)
      DEALLOCATE(J0)
!
!     ==========================================================================
!     ==  NOW DO HE AUGMENTATION                                              ==
!     ==========================================================================
      ALLOCATE(YLM(LMX))
      DO ISP2=1,NSP
        CALL SETUP$ISELECT(ISP2)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
!       == COLLECT PARTIAL WAVES ===============================================
        ALLOCATE(AEPHI(NR,LNX(ISP2)))
        ALLOCATE(AEPHIDOT(NR,LNX(ISP2)))
        ALLOCATE(NLPHI(NR,LNX(ISP2)))
        ALLOCATE(NLPHIDOT(NR,LNX(ISP2)))
        CALL SETUP$GETR8A('AEPHI',NR*LNX(ISP2),AEPHI)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX(ISP2),AEPHIDOT)
        CALL SETUP$GETR8A('NLPHI',NR*LNX(ISP2),NLPHI)
        CALL SETUP$GETR8A('NLPHIDOT',NR*LNX(ISP2),NLPHIDOT)
        CALL SETUP$ISELECT(0)
!
!       == CONSTRUCT HEAD AND TAIL FUNCTIONS FROM PARTIAL WAVES ===============
        NTAIL=POTPAR(ISP2)%NTAIL
        ALLOCATE(K0ARR(NR))
        ALLOCATE(DK0ARR(NR))
        ALLOCATE(JBARARR(NR,NTAIL))
        ALLOCATE(DJBARARR(NR,NTAIL))
        TONSITE=(ISP2.EQ.ISPECIES(IATORB))
        DO IT=1,NTAIL
          L=POTPAR(ISP2)%LOFT(IT)
          LNDOT=POTPAR(ISP2)%LNOFT(IT)
          JBARARR(:,IT)=NLPHIDOT(:,LNDOT)*POTPAR(ISP2)%JBARTOPHIDOT(IT)
          DJBARARR(:,IT)=(AEPHIDOT(:,LNDOT)-NLPHIDOT(:,LNDOT)) &
      &                 *POTPAR(ISP2)%JBARTOPHIDOT(IT)
          IF(TONSITE) THEN
            IF(POTPAR(ISP2)%LOFH(IHORB).EQ.L) THEN
              LN=POTPAR(ISP2)%LNOFH(IHORB)
              K0ARR(:)=NLPHI(:,LN)   *POTPAR(ISP2)%KTOPHI(IHORB) &
                      +NLPHIDOT(:,LNDOT)*POTPAR(ISP2)%KTOPHIDOT(IHORB) 
              DK0ARR(:)=(AEPHI(:,LN)-NLPHI(:,LN))*POTPAR(ISP2)%KTOPHI(IHORB) &
     &            +(AEPHIDOT(:,LN)-NLPHIDOT(:,LN))*POTPAR(ISP2)%KTOPHIDOT(IHORB)
            END IF
          END IF
        ENDDO
        DEALLOCATE(AEPHI)
        DEALLOCATE(AEPHIDOT)
        DEALLOCATE(NLPHI)
        DEALLOCATE(NLPHIDOT)
!
        DO NN=1,NNS
          IF(SBAR(NN)%IAT1.NE.IATORB) CYCLE
          IAT2=SBAR(NN)%IAT2
          IF(ISPECIES(IAT2).NE.ISP2) CYCLE
          R2(:)=R0(:,IAT2)+RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
       &                  +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
       &                  +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
          TONSITE=(IAT2.EQ.IATORB).AND.(MAXVAL(ABS(SBAR(NN)%IT(:))).EQ.0)
!         == SELECT STRUCTURE CONSTANT ARRAY ===================================
          LMN2X=SBAR(NN)%N2
          LM2X=(MAXVAL(POTPAR(ISP2)%LOFT)+1)**2
          SBARVEC(:LMN2X)=SBAR(NN)%MAT(LMNORB,:)  
!  
!         == LOOP OVER REAL SPACE GRID =========================================
          DO IP=1,NP
!           == DR IS THE DISTANCE FROM THE SECOND ATOM =========================
            DR(:)=P(:,IP)-R2(:)
            DIS=SQRT(SUM(DR**2))
            CALL SPHERICAL$YLM(LM2X,DR,YLM)
!
!           == INCLUDE DIFFERENCE BETWEEN AEPHI AND NLPHI ======================
!           == (THEY MAY DIFFER OUTSIDE THE SPHERE DUE TO CORE TAILS) ==========
!           == AUGMENTATION WITH ALL-ELECTRON PARTIAL WAVES ONLY ONSITE ========
            IF(TONSITE) THEN
              CALL RADIAL$VALUE(GID,NR,DK0ARR,DIS,VAL)
              ORB(IP)=ORB(IP)+VAL*YLM(LMORB)
!
              LMN2=0
              DO IT=1,NTAIL
                CALL RADIAL$VALUE(GID,NR,DJBARARR(:,IT),DIS,VALDOT)
                L=POTPAR(ISP2)%LOFT(IT)
                LM2=L**2
                DO IM=1,2*L+1
                  LM2=LM2+1
                  LMN2=LMN2+1
                  ORB(IP)=ORB(IP)-VALDOT*YLM(LM2)*SBARVEC(LMN2) 
                ENDDO
              ENDDO
            END IF
!
            TSPHERE=(DIS.LT.POTPAR(ISP2)%RAD)
            IF(.NOT.TSPHERE) CYCLE
!           == NOW THE AUGMENTATION ============================================
            IF(TONSITE) THEN
              CALL RADIAL$VALUE(GID,NR,K0ARR,DIS,VAL)
              ORB(IP)=ORB(IP)+VAL*YLM(LMORB)
            END IF
            LMN2=0
            DO IT=1,NTAIL
              CALL RADIAL$VALUE(GID,NR,JBARARR(:,IT),DIS,VALDOT)
              L=POTPAR(ISP2)%LOFT(IT)
              LM2=L**2
              DO IM=1,2*L+1
                LM2=LM2+1
                LMN2=LMN2+1
                ORB(IP)=ORB(IP)-VALDOT*YLM(LM2)*SBARVEC(LMN2) 
              ENDDO
            ENDDO
          ENDDO  ! END LOOP IP
        ENDDO  ! END LOOP NN
        DEALLOCATE(R)
        DEALLOCATE(K0ARR)
        DEALLOCATE(DK0ARR)
        DEALLOCATE(JBARARR)
        DEALLOCATE(DJBARARR)
      ENDDO   ! END LOOP NSP
      DEALLOCATE(YLM)
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
      INTEGER(4),INTENT(IN) :: NDIR
      REAL(8)   ,INTENT(IN) :: DIR(3,NDIR)
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
      REAL(8)              :: P(3,NP)
      REAL(8)              :: ORB(NP)
      REAL(8)              :: ORBTAILED(NP)
      INTEGER(4)           :: J,LM
      INTEGER(4),PARAMETER :: LX=4
      INTEGER(4),PARAMETER :: LMX=(LX+1)**2
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
      CALL LMTO_TAILED_NTBOOFR(IATORB,LMNORB,NP,P,ORB)
      CALL LMTO_NTBOOFR(IATORB,LMNORB,NP,P,ORBTAILED)
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
      REAL(8)      ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      COMPLEX(8)   ,PARAMETER   :: CI=(0.D0,1.D0)
      COMPLEX(8)   ,PARAMETER   :: CI2PI=CI*2.D0*PI
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
      LOGICAL(4)                :: TINV
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO$PLOTWAVE')
      IF(.NOT.TON) RETURN
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
      ALLOCATE(CVEC(NPRO))  ! THE LENGTH SHOULD BE NCHI INSTEAD OF NPRO
      CVEC(:)=(0.D0,0.D0)
      IF(TINV) THEN
        DO IBH=1,NBH
          IB1=2*IBH-1
          IB2=2*IBH
          CSVAR1=0.5D0*(THIS%EIGVEC(IB1,IB0)-CI*THIS%EIGVEC(IB2,IB0))
          CSVAR2=0.5D0*(THIS%EIGVEC(IB1,IB0)+CI*THIS%EIGVEC(IB2,IB0))
          CVEC(:)=CVEC(:)+THIS%TBC_NEW(IDIM0,IBH,:)*CSVAR1 &
      &            +CONJG(THIS%TBC_NEW(IDIM0,IBH,:))*CSVAR2
        ENDDO
      ELSE
        DO IB1=1,NB
          CSVAR1=THIS%EIGVEC(IB1,IB0)
          CVEC(:)=CVEC(:)+THIS%TBC_NEW(IDIM0,IB1,:)*CSVAR1
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
      SUBROUTINE LMTO_EXPANDSBAR()
!     **************************************************************************
!     **  THE SCREENED STRUCTURE CONSTANTS ARE INITIALLY GIVEN IN THE SPACE   **
!     **  WITH DIMENSIONS RELATED TO THE ANGULAR MOMENTA.                     **
!     **  HERE THEY ARE CONVERTED SO THAT THE RIGHT INDEX IS THE NUMBER OF    **
!     **  HEAD FUNCTIONS AND THE LEFT IS THE NUMBER OF TAIL FUNCTIONS         **
!     **                                                                      **
!     **  P_{R'}|K_{RL}> = - \SUM_{L'} |JBAR_{R'L'}> S_{RL,R'L'}              **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR1  &
     &                       ,SBAR_NEW & ! SCREENED STRUCTURE CONSTANTS
     &                       ,SBARLI1  &
     &                       ,ISPECIES & ! ISP=ISPECIES(IAT)
     &                       ,NSP      & ! NUMBER OF ATOM TYPES
     &                       ,LNX,LOX  
      IMPLICIT NONE
      LOGICAL(4),PARAMETER :: TPR=.FALSE. 
      INTEGER(4),PARAMETER :: NFIL=6
      INTEGER(4) :: NNS  ! #(NEIGBOR PAIRS)
      INTEGER(4) :: IAT1,IAT2
      INTEGER(4) :: ISP1,ISP2
      INTEGER(4) :: LMX(NSP)
      INTEGER(4) :: LX
      INTEGER(4) :: NHEAD,NTAIL
      INTEGER(4) :: N1,N2
      INTEGER(4) :: NN,ISP,L
      INTEGER(4) :: I1,L1,M1,LM1A,LM1B,LM1X
      INTEGER(4) :: I2,L2,M2,LM2A,LM2B,LM2X
      REAL(8)   ,ALLOCATABLE :: MAT(:,:)
!     **************************************************************************
      DO ISP=1,NSP
        LX=MAXVAL(LOX(:LNX(ISP),ISP))
        LMX(ISP)=0
        DO L=0,LX
          LMX(ISP)=LMX(ISP)+2*L+1
        ENDDO
      ENDDO
      NNS=SIZE(SBAR_NEW)
      DO NN=1,NNS
        IAT1=SBAR_NEW(NN)%IAT1
        IAT2=SBAR_NEW(NN)%IAT2
        N1=SBAR_NEW(NN)%N1
        N2=SBAR_NEW(NN)%N2
        ALLOCATE(MAT(N1,N2))
        MAT=SBAR_NEW(NN)%MAT
        DEALLOCATE(SBAR_NEW(NN)%MAT)
!
        ISP1=ISPECIES(IAT1)
        ISP2=ISPECIES(IAT2)
        NTAIL=POTPAR1(ISP2)%NTAIL
        NHEAD=POTPAR1(ISP1)%NHEAD
        LM1X=SUM(2*POTPAR1(ISP1)%LOFH+1)
        LM2X=SUM(2*POTPAR1(ISP2)%LOFT+1)
        SBAR_NEW(NN)%N1=LM1X
        SBAR_NEW(NN)%N2=LM2X
        ALLOCATE(SBAR_NEW(NN)%MAT(LM1X,LM2X))
        LM1A=0
        DO I1=1,NHEAD
          L1=POTPAR1(ISP1)%LOFH(I1)
          LM1B=SBARLI1(L1+1,ISP1)-1
          DO M1=1,2*L1+1
            LM1A=LM1A+1
            LM1B=LM1B+1
            LM2A=0
            DO I2=1,NTAIL
              L2=POTPAR1(ISP2)%LOFT(I2)
              LM2B=SBARLI1(L2+1,ISP2)-1
              DO M2=1,2*L2+1
                LM2A=LM2A+1
                LM2B=LM2B+1
                SBAR_NEW(NN)%MAT(LM1A,LM2A)=MAT(LM1B,LM2B)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(MAT)
      ENDDO
!
!     ==========================================================================
!     ==  PRINT RESULT                                                        ==
!     ==========================================================================
      IF(TPR) THEN
        CALL LMTO$REPORTSBAR(NFIL)
        WRITE(NFIL,FMT='(82("="))')
        WRITE(NFIL,FMT='(82("="),T10,"  SCREENED STRUCTURE CONSTANTS(NEW)  ")')
        WRITE(NFIL,FMT='(82("="))')
        DO NN=1,NNS
          WRITE(NFIL &
     &         ,FMT='(82("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 SBAR_NEW(NN)%IAT1,SBAR_NEW(NN)%IAT2,SBAR_NEW(NN)%IT
          DO LM1A=1,SBAR_NEW(NN)%N1
            WRITE(NFIL,FMT='(20F10.5)')SBAR_NEW(NN)%MAT(LM1A,:)
          ENDDO
          WRITE(NFIL,FMT='(82("="))')
        ENDDO 
        CALL ERROR$MSG('PLANNED STOP AFTER WRITING SBAR_NEW')
        CALL ERROR$STOP('LMTO_EXPANDSBAR')
      END IF
      RETURN
      END
