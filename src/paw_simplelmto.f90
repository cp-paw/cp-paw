!*******************************************************************************
!*******************************************************************************
!**  INTERFACE FOR ELECTRONIC CORRELATIONS AND HYBRID FUNCTIONALS             **
!**  EXPRESSED IN TERMS OF LOCAL ORBITALS                                     **
!**                                                                           **
!**  SCHEMATIC FLOW DIAGRAM: AS SEEN FROM PAW_WAVES OBJECT                    **
!**     SIMPLELMTO$MAKE$POTPAR1                                               **
!**                                                                           **
!**     WAVES$OFFSITEDENMAT               OSDENMAT   <- THIS%PROJ;            **
!**                                       (RE)ALLOCATE(OSDENMAT)              **
!**                                       DEALLOCATE(OSHAMIL)                 **
!**     SIMPLELMTO$MAKE$ETOT              ETOT,OSHAMIL <- OSDENMAT            **
!**         WAVES$GETI4('NND')                <- SIZE(OSDENMAT)               **
!**         WAVES$GETRSPACEMATA('DEMAT')     <- OSDENMAT                     **
!**         WAVES$SETRSPACEMATA('HAMIL')      -> OSHAMIL (ALLOCATE(OSHAMIL)   **
!**     WAVES$OFFSITEHAMIL                THIS%HPROJ <- OSHAMIL               **
!**     WAVES$FORCE|WAVES_FORCE_ADDHTBC   FORCE      <- THIS%HPROJ            **
!**     WAVES$HPSI                        HPSI       <- THIS%HPROJ            **
!**                                                  (DEALLOCATE THIS%HPROJ)  **
!**                                                                           **
!**     WAVES$OFFSITEHAMIL AND WAVES$HPROJ ARE ANALOGOUS                      **
!**                                                                           **
!*******************************************************************************
!****************************************P.BLOECHL, GOSLAR 2019 ****************
MODULE SIMPLELMTO_MODULE
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
  REAL(8)     :: RAUG       ! AUGMENTATION RADIUS
  INTEGER(4),POINTER :: NORBOFL(:) ! #(LOCAL ORBITALS PER ANGULAR MOMENTUM)
END TYPE HYBRIDSETTING_TYPE
!
!== HOLDS THE POTENTIAL PARAMETER FOR ONE ATOM TYPE ============================
TYPE POTPAR_TYPE
  ! THE NUMBER OF ACTIVE HEAD FUNCTIONS IS NPHI. 
  ! THE NUMBER OF TAIL FUNCTIONS IS LX.
  ! EACH HEAD FUNCTION IS RELATED TO A PARTIAL WAVE IDENTIFIED BY LN.
  REAL(8)            :: RAUG            ! MATCHING RADIUS
  INTEGER(4)         :: LNXH            ! #(HEAD FUNCTIONS)
  INTEGER(4)         :: LMNXH           ! #(HEAD FUNCTIONS INCLUDING M-VALUES)
  INTEGER(4)         :: NTAIL           ! #(TAIL FUNCTIONS)
  INTEGER(4),POINTER :: ITAIL(:)        !(LNXH) POINTER TO TAIL FUNCTION
  INTEGER(4),POINTER :: LOXH(:)         !(LNXH) MAIN ANGULAR MOMENTUM
  INTEGER(4),POINTER :: LNOFH(:)        !(LNXH) PARTIAL WAVE ID
  REAL(8)   ,POINTER :: KTOPHI(:)       !(LNXH) |K> = |PHI>    * KTOPHI
  REAL(8)   ,POINTER :: KTOPHIDOT(:)    !(LNXH)     + |PHIDOT> * KTOPHIDOT
  REAL(8)   ,POINTER :: PHIDOTPROJ(:)   !(LNX)   <P(LN)|PHIDOT(ITAIL(IHEAD))>
  INTEGER(4),POINTER :: LOFT(:)         !(NTAIL)  MAIN ANGULAR MOMENTUM
  INTEGER(4),POINTER :: LNOFT(:)        !(NTAIL)  PARTIAL WAVE ID FOR PHIDOT
  REAL(8)   ,POINTER :: QBAR(:)         !(NTAIL)  |JBAR>=|J>-|K>QBAR
  REAL(8)   ,POINTER :: JBARTOPHIDOT(:) !(NTAIL)|JBAR> = |PHIBARDOT> JBARTOPHIDOT
  REAL(8)   ,POINTER :: PROK(:,:)       !(LNX,LNXH) <P|K_AUG>
  REAL(8)   ,POINTER :: PROJBAR(:,:)    !(LNX,NTAIL) <P|JBAR_AUG> 
  REAL(8)   ,POINTER :: PHIOV(:,:)      !(LNX,LNX) <AEPHI|THETA_OMEGA|AEPHI>
  INTEGER(4)         :: GID
  REAL(8)   ,POINTER :: AECHI(:,:)      !(NR,LNXH)
  REAL(8)   ,POINTER :: PSCHI(:,:)      !(NR,LNXH)
  REAL(8)   ,POINTER :: NLCHI(:,:)      !(NR,LNXH)
  REAL(8),ALLOCATABLE :: PIPHI(:,:)      !(LMNXH,LMNXPHI) <PI|PHI>
  REAL(8),ALLOCATABLE :: OVERLAP(:,:)    !(LMNXH,LMNXH)    <CHI|CHI>
  REAL(8),ALLOCATABLE :: ONSITEU(:,:,:,:)!(LMNXH,LMNXH,LMNXH,LMNXH)
  REAL(8),ALLOCATABLE :: QLN(:,:,:) !(2,LNXH,LNXH) MONO-&DIPOLE MATRIX ELEMENTS
END TYPE POTPAR_TYPE
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
REAL(8)               :: SCREENL=-1.D0  !SCREENING LENGTH
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
TYPE(POTPAR_TYPE)  ,ALLOCATABLE :: POTPAR(:) !POTENTIAL PARAMETERS (NEW)
INTEGER(4)         ,ALLOCATABLE :: SBARLI1(:,:)
TYPE(UTENSOR_TYPE) ,ALLOCATABLE :: UTENSOR(:)
TYPE(OFFSITEX_TYPE),ALLOCATABLE :: OFFSITEX(:,:)
LOGICAL(4)         ,PARAMETER   :: TSPHERICAL=.FALSE.
END MODULE SIMPLELMTO_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO$SETL4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : TON &
     &                             ,TOFFSITE &
     &                             ,ISPSELECTOR &
     &                             ,HYBRIDSETTING
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'ON') THEN
        TON=VAL
!
      ELSE IF(ID.EQ.'OFFSITE') THEN
        TOFFSITE=VAL
!
!     ==========================================================================
!     == IF ACTIVE, THE HYBRID CONTRIBTIONS ON THIS ATOM ARE CONSIDERED       ==
!     ==========================================================================
      ELSE IF(ID.EQ.'ACTIVE') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%ACTIVE=VAL
!
      ELSE IF(ID.EQ.'COREVALENCE') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TCV=VAL
!
      ELSE IF(ID.EQ.'FOCKSETUP') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TFOCKSETUP=VAL
!
      ELSE IF(ID.EQ.'NDDO') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TNDDO=VAL
!
      ELSE IF(ID.EQ.'31') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%T31=VAL
!
      ELSE IF(ID.EQ.'BONDX') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TBONDX=VAL
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SIMPLELMTO$SETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO$GETL4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : TON &
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
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$GETL4')
        END IF
        VAL=HYBRIDSETTING(ISPSELECTOR)%TFOCKSETUP
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SIMPLELMTO$GETL4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO$GETI4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : ISPECIES & !(NAT)
     &                       ,POTPAR   !(NSP)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(OUT) :: VAL
      INTEGER(4)               :: IAT,ISP
!     **************************************************************************
      IF(ID.EQ.'NLOCORB') THEN
        VAL=0
        DO IAT=1,SIZE(ISPECIES)
          ISP=ISPECIES(IAT)
          VAL=VAL+SUM(2*POTPAR(ISP)%LOXH+1)
        ENDDO
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SIMPLELMTO$GETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO$SETI4A(ID,LEN,VAL)
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : ISPECIES & !(NAT)
     &                       ,POTPAR &  !(NSP)
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
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETI4A')
        END IF
        ALLOCATE(HYBRIDSETTING(ISPSELECTOR)%NORBOFL(LEN))
        HYBRIDSETTING(ISPSELECTOR)%NORBOFL=VAL
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SIMPLELMTO$SETI4A')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO$SETI4(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : HYBRIDSETTING &
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
          CALL ERROR$STOP('SIMPLELMTO$SETI4')
        END IF
        IF(VAL.GT.NSP) THEN
          CALL ERROR$MSG('ATOM TYPE SPECIFIER OUT OF RANGE')
          CALL ERROR$I4VAL('ISP',ISPSELECTOR)
          CALL ERROR$I4VAL('NSP',NSP)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETI4')
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
          HYBRIDSETTING(:)%RAUG      =-1.D0
          DO I=1,NSP
            NULLIFY(HYBRIDSETTING(I)%NORBOFL)
          ENDDO
        END IF
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SIMPLELMTO$SETI4')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO$SETR8(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : ISPSELECTOR &
     &                             ,HYBRIDSETTING &
     &                             ,HFWEIGHT &
     &                             ,K2 &
     &                             ,SCREENL &
     &                             ,RCSCALE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
      INTEGER(4)              :: I
!     **************************************************************************
!
      IF(ID.EQ.'LHFWEIGHT') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%LHFWEIGHT=VAL
!
      ELSE IF(ID.EQ.'RAUG') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$CHVAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%RAUG=VAL
!
      ELSE IF(ID.EQ.'SCALERCUT') THEN
        RCSCALE=VAL
!
      ELSE IF(ID.EQ.'HFWEIGHT') THEN
        HFWEIGHT=VAL
!
      ELSE IF(ID.EQ.'SCREENL') THEN
        SCREENL=VAL
!
      ELSE IF(ID.EQ.'K2') THEN
        K2=VAL

      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SIMPLELMTO$SETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO$GETR8(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : HFWEIGHT &
     &                             ,NSP &
     &                             ,HYBRIDSETTING &
     &                             ,ISPSELECTOR
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
          CALL ERROR$STOP('SIMPLELMTO$GETR8')
        END IF
        VAL=HYBRIDSETTING(ISPSELECTOR)%LHFWEIGHT
        IF(VAL.LT.0.D0) VAL=HFWEIGHT

      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SIMPLELMTO$GETR8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO$SETCH(ID,VAL)
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : MODUS
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
        CALL ERROR$STOP('SIMPLELMTO$SETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO$INITIALIZE()
!     **************************************************************************
!     **  PREPARES POTENTIAL PARAMETERS AND SIMILAR BASIC DATA.               **
!     **  IT IS CALLED BY SIMPLELMTO$MAKESTRUCTURECONSTANTS                   **
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : TINI &
     &                             ,TON &
     &                             ,TOFFSITE 
      IMPLICIT NONE
      INTEGER(4) :: NAT,ISP
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
                          CALL TRACE$PUSH('SIMPLELMTO$INITIALIZE')
!
!     ==========================================================================
!     == COLLECT NSP,LNX,LOX,ISPECIES                                         ==
!     ==========================================================================
      CALL SIMPLELMTO_COLLECTMAPARRAYS()
!
      CALL SIMPLELMTO_CONSOLIDATEHYBRIDSETTING()
!
!     ==========================================================================
!     == DETERMINE POTENTIAL PARAMETERS                                       ==
!     == RAD,LNSCATT,PHIDOTPROJ,QBAR,KTOPHI,KTOPHIDOT,JBARTOPHIDOT            ==
!     ==========================================================================
      CALL SIMPLELMTO_MAKEPOTPAR()
!
!     ==========================================================================
!     == ATTACH EXPONENTIAL TAILS TO AUGMENTED HANKEL AND BESSEL FUNCTIONS    ==
!     ==========================================================================
!!$      CALL SIMPLELMTO_MAKETAILEDPARTIALWAVES_WITHPOTPAR()
!!$      CALL SIMPLELMTO_MAKETAILEDMATRIXELEMENTS_WITHPOTPAR()
!!$      IF(.NOT.TON) RETURN
!!$!
!!$!     ==========================================================================
!!$!     ==  CONSTRUCT OFFSITE INTEGRALS OF TAILED ORBITALS                      ==
!!$!     ==========================================================================
!!$      IF(TOFFSITE) THEN 
!!$                            CALL TIMING$CLOCKON('OFFSITE U-TENSOR')
!!$!       == GAUSSIAN FIT OF TAILED ORBITAL FOR ABAB-TYPE U-TENSOR ===============
!!$                            CALL TIMING$CLOCKON('OFFS.U. GAUSSFIT')
!!$        CALL SIMPLELMTO_TAILEDGAUSSFIT()
!!$                            CALL TIMING$CLOCKOFF('OFFS.U. GAUSSFIT')
!!$!
!!$!       == OFF-SITE MATRIX ELEMENTS OF OVERLAP MATRIX AND U-TENSOR =============
!!$!       == OVERLAP MATRIX ELEMENTS ARE ALWAYS COMPUTED, U-TENSOR ELEMENTS ======
!!$!       == ONLY WHEN REQUESTED =================================================
!!$                            CALL TIMING$CLOCKON('OFFS.U. INTEGRALS')
!!$        CALL SIMPLELMTO_OFFXINT()
!!$                            CALL TIMING$CLOCKOFF('OFFS.U. INTEGRALS')
!!$                            CALL TIMING$CLOCKOFF('OFFSITE U-TENSOR')
!!$      END IF
                          CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_COLLECTMAPARRAYS()
!     **************************************************************************
!     **  STORES A LOCAL COPY OF                                              **
!     **     NSP                                                              **
!     **     LNX(NSP)                                                         **
!     **     LOX(LNXX,NSP)                                                    **
!     **     ISPECIES(NAT)                                                    **
!     **  IN THE SIMPLELMTO_MODULE                                            **
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY: NSP &
     &                      ,LNX &
     &                      ,LOX &
     &                      ,ISPECIES &
     &                      ,ISPSELECTOR
      IMPLICIT NONE
      INTEGER(4)    :: ISP,NAT
!     **************************************************************************
                              CALL TRACE$PUSH('SIMPLELMTO_COLLECTMAPARRAYS')
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
        CALL SIMPLELMTO$SETI4('ISP',1)
        CALL SIMPLELMTO$SETI4('ISP',0)
      END IF
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_CONSOLIDATEHYBRIDSETTING()
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
                          CALL TRACE$PUSH('SIMPLELMTO_CONSOLIDATEHYBRIDSETTING')
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
          CALL ERROR$MSG('RAUG NOT SET. STOPPING')
          CALL ERROR$I4VAL('ISP',ISP)
          CALL ERROR$R8VAL('HYBRIDSETTING(ISP)%RAUG',HYBRIDSETTING(ISP)%RAUG)
          CALL ERROR$STOP('SIMPLELMTO_CONSOLIDATEHYBRIDSETTING')
          HYBRIDSETTING(ISP)%RAUG=1.1D0*RCOV
        END IF
!
        IF(HYBRIDSETTING(ISP)%TBONDX) THEN
          CALL ERROR$MSG('BOND EXCHANGE (BONDX) IS NOT IMPLEMENTED')
          CALL ERROR$I4VAL('ISP',ISP)
         CALL ERROR$L4VAL('HYBRIDSETTING(ISP)%TBONDX',HYBRIDSETTING(ISP)%TBONDX)
          CALL ERROR$STOP('SIMPLELMTO_CONSOLIDATEHYBRIDSETTING')
        END IF
      ENDDO
                              CALL TRACE$POP()
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO$REPORT(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : TON,NSP,HYBRIDSETTING,RCSCALE,HFWEIGHT,K2 &
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
      CALL SIMPLELMTO$INITIALIZE()
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      CALL REPORT$TITLE(NFIL,'SIMPLELMTO OBJECT:  GENERIC VARIABLES')
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
        CALL REPORT$TITLE(NFIL,'SIMPLELMTO OBJECT: '//TRIM(ID))
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
        DEALLOCATE(LOX1)
        CALL SETUP$UNSELECT()
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_MAKEPOTPAR()
!     **************************************************************************
!     **  SIMILAR TO SIMPLELMTO_MAKEPOTPAR BUT WITH ANOTHER DATA STRUCTURE    **
!     **  ORGANIZED ACCORDING TO HEADS AND TAILS INSTEAD OF PARTIAL WAVES.    **
!     **  IT IS INTENDED THAT THIS ROUTINE REPLACES SIMPLELMTO_MAKEPOTPAR.    **
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
      USE SIMPLELMTO_MODULE, ONLY : K2 &
     &                             ,POTPAR &
     &                             ,NSP &
     &                             ,LNX &
     &                             ,LOX &
     &                             ,HYBRIDSETTING
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
      REAL(8)                :: RAUG
      REAL(8)                :: PHIVAL,PHIDER
      REAL(8)                :: PHIDOTVAL,PHIDOTDER
      REAL(8)                :: KVAL,KDER
      REAL(8)                :: JVAL,JDER
      REAL(8)                :: WJPHI,WJPHIDOT,WKPHI,WKPHIDOT,WJBARPHI,WKJ
      REAL(8)                :: WPHIPHIDOT
      REAL(8)                :: QBAR
      REAL(8)                :: SVAR
      INTEGER(4)             :: LNXH !#(HEAD FUNCTIONS)
      INTEGER(4)             :: NTAIL !#(TAIL FUNCTIONS)
      INTEGER(4)             :: LX    !X(ANGULAR MOMENTUM)
      INTEGER(4)             :: ISP,LN,L,LNOFH,LNOFT,IHEAD,ITAIL
      INTEGER(4)             :: ISVAR
      LOGICAL(4)             :: TCHK
!     **************************************************************************
                             CALL TRACE$PUSH('SIMPLELMTO_MAKEPOTPAR')
      IF(.NOT.ALLOCATED(HYBRIDSETTING)) THEN
        CALL ERROR$MSG('HYBRIDSETTING NOT ALLOCATED')
        CALL ERROR$STOP('SIMPLELMTO_MAKEPOTPAR')
      END IF

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
        ALLOCATE(AUX(NR))
        POTPAR(ISP)%GID=GID    !RADIAL GRID IS NEEDED FOR THE LOCAL ORBITALS
!
!       ==
        ALLOCATE(ISCATT(LNX1))
        CALL SETUP$GETI4A('ISCATT',LNX1,ISCATT)

!       == MATCHING RADIUS =====================================================
        IF(HYBRIDSETTING(ISP)%RAUG.GE.0.D0) THEN
          RAUG=HYBRIDSETTING(ISP)%RAUG
        ELSE
          CALL SETUP$GETR8('AEZ',SVAR)
          CALL PERIODICTABLE$GET(SVAR,'R(COV)',RAUG)
        END IF
        POTPAR(ISP)%RAUG=RAUG
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
              CALL ERROR$STOP('SIMPLELMTO_MAKEPOTPAR')
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
        LNXH=0
        DO LN=1,LNX1
          IF(.NOT.TORB(LN))CYCLE
          LNXH=LNXH+1
          LX=MAX(LX,LOX(LN,ISP))
        ENDDO
        POTPAR(ISP)%LNXH=LNXH 
       ALLOCATE(POTPAR(ISP)%LNOFH(LNXH))
        ALLOCATE(POTPAR(ISP)%LOXH(LNXH))
        ALLOCATE(POTPAR(ISP)%ITAIL(LNXH))
        ALLOCATE(POTPAR(ISP)%KTOPHI(LNXH))
        ALLOCATE(POTPAR(ISP)%KTOPHIDOT(LNXH))
        ALLOCATE(POTPAR(ISP)%PROK(LNX1,LNXH))
        ALLOCATE(POTPAR(ISP)%PHIDOTPROJ(LNXH))
        POTPAR(ISP)%PROK(:,:)=0.D0
        IHEAD=0
        DO LN=1,LNX1
          IF(.NOT.TORB(LN))CYCLE
          IHEAD=IHEAD+1
          POTPAR(ISP)%LNOFH(IHEAD)=LN
          POTPAR(ISP)%LOXH(IHEAD)=LOX(LN,ISP)
        ENDDO
        POTPAR(ISP)%LMNXH=SUM(2*POTPAR(ISP)%LOXH(:)+1)
!
!       ========================================================================
!       == COUNT NUMBER OF TAIL FUNCTIONS ======================================
!       == ONE TAIL FUNCTION FOR EACH ANGULAR MOMENTUM WITH AN ORBITAL        ==
!       ========================================================================
        NTAIL=0
        LX=MAXVAL(LOX(:LNX(ISP),ISP))
        DO L=0,LX
          DO IHEAD=1,LNXH
            IF(POTPAR(ISP)%LOXH(IHEAD).NE.L) CYCLE
            NTAIL=NTAIL+1
            EXIT
          ENDDO
        ENDDO
        POTPAR(ISP)%NTAIL=NTAIL
        ALLOCATE(POTPAR(ISP)%LOFT(NTAIL))
        ALLOCATE(POTPAR(ISP)%LNOFT(NTAIL))
        ALLOCATE(POTPAR(ISP)%QBAR(NTAIL))
        ALLOCATE(POTPAR(ISP)%JBARTOPHIDOT(NTAIL))
        ALLOCATE(POTPAR(ISP)%PROJBAR(LNX1,NTAIL))
        POTPAR(ISP)%PROJBAR(:,:)=0.D0
!       == CONNECT EACH TAIL FUNCTION TO A SPECIFIC HEAD FUNCTION FROM WHICH ===
!       == THE SCATTERING FUNCTION IS TAKEN TO DEFINE THE SCREENING CHARGE QBAR.
!       == THE SELECTED HEAD FUNCTION IS THE LAST ONE FOR THE GIVEN L-VALUE.====
!       == IT IS IDENTIFIED BY LNOFT
        ITAIL=0
        DO L=0,LX
          ITAIL=ITAIL+1
          TCHK=.FALSE.
          DO IHEAD=1,LNXH
            IF(POTPAR(ISP)%LOXH(IHEAD).NE.L) CYCLE
            TCHK=.TRUE.
            POTPAR(ISP)%LOFT(ITAIL)=L
            POTPAR(ISP)%LNOFT(ITAIL)=POTPAR(ISP)%LNOFH(IHEAD)  !PLACEHOLDER
          ENDDO
          IF(.NOT.TCHK) ITAIL=ITAIL-1  ! NO TAIL FORT THIS L/ DO NOT COUNT UP 
        ENDDO
!       == LINK CORRESPONDING PHIDOT TO EACH PHI ===============================
        DO IHEAD=1,LNXH
          L=POTPAR(ISP)%LOXH(IHEAD)
          DO ITAIL=1,NTAIL
            IF(POTPAR(ISP)%LOFT(ITAIL).NE.L) CYCLE
            POTPAR(ISP)%ITAIL(IHEAD)=ITAIL
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
          L=POTPAR(ISP)%LOFT(ITAIL)
          LNOFT=POTPAR(ISP)%LNOFT(ITAIL)
          DO IHEAD=1,LNXH
            IF(POTPAR(ISP)%ITAIL(IHEAD).NE.ITAIL) CYCLE
            LN=POTPAR(ISP)%LNOFH(IHEAD)
            IF(ISCATT(LN).GT.ISCATT(LNOFT)) LNOFT=LN
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  DETERMINE POTENTIAL PARAMETERS                                    ==
!       ========================================================================
        DO ITAIL=1,NTAIL
          L=POTPAR(ISP)%LOFT(ITAIL)
          LNOFT=POTPAR(ISP)%LNOFT(ITAIL)
!         ======================================================================
!         == VALUE AND DERIVATIVE OF PARTIAL WAVES AND ENVELOPE FUNCTIONS     ==
!         == PHIDOT IS PHIBARDOT                                              ==
!         ======================================================================
          CALL RADIAL$VALUE(GID,NR,NLPHIDOT(:,LNOFT),RAUG,PHIDOTVAL)
          CALL RADIAL$DERIVATIVE(GID,NR,NLPHIDOT(:,LNOFT),RAUG,PHIDOTDER)
          CALL LMTO$SOLIDBESSELRAD(L,RAUG,K2,JVAL,JDER)
          CALL LMTO$SOLIDHANKELRAD(L,RAUG,K2,KVAL,KDER)
          WJPHIDOT=JVAL*PHIDOTDER-JDER*PHIDOTVAL
          WKPHIDOT=KVAL*PHIDOTDER-KDER*PHIDOTVAL
          QBAR=WJPHIDOT/WKPHIDOT
!         == |JBAR> = |J> - |K>QBAR ============================================
          POTPAR(ISP)%QBAR(ITAIL)    = QBAR
!         == JBAR -> |PHIBARDOT> JBARTOPHIDOT ==================================
          POTPAR(ISP)%JBARTOPHIDOT(ITAIL)=(JVAL-KVAL*QBAR)/PHIDOTVAL
!
!         ==  <PRO|JBAR_AUG> ===================================================
          AUX=R**2*PSPHIDOT(:,LNOFT)*POTPAR(ISP)%JBARTOPHIDOT(ITAIL)
          DO LN=1,LNX1
            IF(LOX(LN,ISP).NE.L) CYCLE
            CALL RADIAL$INTEGRAL(GID,NR,PRO(:,LN)*AUX,SVAR)
            POTPAR(ISP)%PROJBAR(LN,ITAIL)=SVAR
          ENDDO
!
          DO IHEAD=1,POTPAR(ISP)%LNXH
            IF(POTPAR(ISP)%LOXH(IHEAD).NE.L) CYCLE
            LNOFH=POTPAR(ISP)%LNOFH(IHEAD)
!           ====================================================================
!           == VALUE AND DERIVATIVE OF PARTIAL WAVES                          ==
!           == PHIDOT IS PHIBARDOT                                            ==
!           ====================================================================
            CALL RADIAL$VALUE(GID,NR,NLPHI(:,LNOFH),RAUG,PHIVAL)
            CALL RADIAL$DERIVATIVE(GID,NR,NLPHI(:,LNOFH),RAUG,PHIDER)
            WJPHI=JVAL*PHIDER-JDER*PHIVAL
            WKPHI=KVAL*PHIDER-KDER*PHIVAL
            WPHIPHIDOT=PHIVAL*PHIDOTDER-PHIDER*PHIDOTVAL
            WJBARPHI=WJPHI-WKPHI*QBAR
!           == K    -> |PHI>KTOPHI+|PHIBARDOT> KTOPHIDOT =======================
            POTPAR(ISP)%KTOPHIDOT(IHEAD)=-WKPHI/WPHIPHIDOT
            POTPAR(ISP)%KTOPHI(IHEAD)   = WKPHIDOT/WPHIPHIDOT
!
!           ==  <PRO|PSPHIDOT> =================================================
            AUX=R**2*PRO(:,LNOFH)*PSPHIDOT(:,LNOFT)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            POTPAR(ISP)%PHIDOTPROJ(IHEAD)=SVAR
!
!           ==  <PRO|K_AUG> ===================================================
            AUX=R**2*(PSPHI(:,LNOFH)   *POTPAR(ISP)%KTOPHI(IHEAD) &
      &              +PSPHIDOT(:,LNOFT)*POTPAR(ISP)%KTOPHIDOT(IHEAD) )
            DO LN=1,LNX1
              IF(LOX(LN,ISP).NE.L) CYCLE
              CALL RADIAL$INTEGRAL(GID,NR,PRO(:,LN)*AUX,SVAR)
              POTPAR(ISP)%PROK(LN,IHEAD)=SVAR
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
!!$PRINT*,'KTOPHI       ',POTPAR(ISP)%KTOPHI(IHEAD)    
!!$PRINT*,'KTOPHIDOT    ',POTPAR(ISP)%KTOPHIDOT(IHEAD)
!!$PRINT*,'JBARTOPHIDOT ',POTPAR(ISP)%JBARTOPHIDOT(IHEAD)
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
!     == PRINTOUT                                                             ==
!     ==========================================================================
      IF(TPRINT) THEN
        WRITE(*,FMT='(80("="),T10," SIMPLELMTO_MAKEPOTPAR  ")')
        DO ISP=1,NSP
          WRITE(*,FMT='(80("-"),T10," ATOM TYPE ",I2," ")')ISP
          NTAIL=POTPAR(ISP)%NTAIL
          DO ITAIL=1,NTAIL
            L=POTPAR(ISP)%LOFT(ITAIL)
            WRITE(*,FMT='("L=",I2," QBAR=",F10.4," JBARTOPHIDOT=",F10.4)') &
     &                   L,POTPAR(ISP)%QBAR(ITAIL) &
     &                    ,POTPAR(ISP)%JBARTOPHIDOT(ITAIL)
          ENDDO
!
          LNXH=POTPAR(ISP)%LNXH
          DO IHEAD=1,LNXH
            L=POTPAR(ISP)%LOXH(IHEAD)
            WRITE(*,FMT='("L=",I2," K->PHI=",F10.4," K->PHIDOT=",F10.4' &
     &                  //'," <PHIDOT|P>=",F10.4)') &
     &              L,POTPAR(ISP)%KTOPHI(IHEAD),POTPAR(ISP)%KTOPHIDOT(IHEAD) &
     &               ,POTPAR(ISP)%PHIDOTPROJ(IHEAD)
          ENDDO
        ENDDO ! END OF LOOP OVER ATOM TYPES
        WRITE(*,FMT='(82("="),T10," FINISHED  ")')
      END IF
!
!     ==========================================================================
!     == CONSTRUCT LOCAL ORBITALS                                             ==
!     ==========================================================================
      CALL SIMPLELMTO_MAKECHI()
!
!     ==========================================================================
!     == <PITILDE|PHITILDE>     
!     ==========================================================================
      CALL SIMPLELMTO_OVERLAPPHI()
      CALL SIMPLELMTO_PIPHI()
!
!     ==========================================================================
!     == CONSTRUCT ONSITE MATRIX ELEMENTS WITH LOCAL ORBITALS                 ==
!     ==========================================================================
      CALL SIMPLELMTO_ONSITEMATRIXELEMENTS()
!
!     ==========================================================================
!     == CONSTRUCT OFF-SITE MATRIX ELEMENTS                                   ==
!     ==========================================================================
      CALL SIMPLELMTO_OFFXINT()

                             CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_MAKECHI()
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
      USE SIMPLELMTO_MODULE, ONLY : K2 &
     &                       ,POTPAR &
     &                       ,NSP &
     &                       ,LNX &
     &                       ,LOX &
     &                       ,HYBRIDSETTING 
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.TRUE.
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)                :: RAUG
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
      INTEGER(4)             :: LNXH ! #(HEAD FUNCTIONS)
      INTEGER(4)             :: NTAIL ! #(TAIL FUNCTIONS)
      INTEGER(4)             :: ISP,LN,LN1,LN2,LNT,LMN,LMN1,LMN2,IM,IR,LNDOT
      INTEGER(4)             :: IH,IT
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      CHARACTER(64)           :: STRING
!     **************************************************************************
                                           CALL TRACE$PUSH('SIMPLELMTO_MAKECHI')
      DO ISP=1,NSP
!       == RADIAL GRID =========================================================
        GID=POTPAR(ISP)%GID
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        RAUG=POTPAR(ISP)%RAUG    !AUGMENTATION RADIUS
        DO IR=1,NR
          IRAD=IR
          IF(R(IR).GT.RAUG) EXIT
        ENDDO
        LNXH=POTPAR(ISP)%LNXH
!
!       ========================================================================
!       == AUGMENTED HANKEL AND BESSEL FUNCTIONS WITH TAILS ATTACHED          ==
!       ========================================================================
        ALLOCATE(POTPAR(ISP)%AECHI(NR,LNXH))
        ALLOCATE(POTPAR(ISP)%PSCHI(NR,LNXH))
        ALLOCATE(POTPAR(ISP)%NLCHI(NR,LNXH))
!
!       ========================================================================
!       == INITIALIZE NLF FOR R<RTAIL AS THE ENVELOPE FUNCTION                ==
!       == THE HEAD FUNCTIONS ARE PURE HANKEL FUNCTIONS                       ==
!       == THE TAIL FUNCTIONS ARE SCREENED BESSEL FUNCTIONS                   ==
!       ========================================================================
        POTPAR(ISP)%NLCHI(:,:)=0.D0
        DO IH=1,LNXH
          LNT=IH
          L=POTPAR(ISP)%LOXH(IH)
          DO IR=IRAD,NR   ! WILL BE AUGMENTED FROM 1 TO IRAD-1
            CALL LMTO$SOLIDHANKELRAD(L,R(IR),K2,KVAL,KDER)
            POTPAR(ISP)%NLCHI(IR,LNT)=KVAL
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
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8A('AEPHI',NR*LNX(ISP),AEPHI)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX(ISP),AEPHIDOT)
        CALL SETUP$GETR8A('QPHI',NR*LNX(ISP),NLPHI)
        CALL SETUP$GETR8A('QPHIDOT',NR*LNX(ISP),NLPHIDOT)
        CALL SETUP$GETR8A('PSPHI',NR*LNX(ISP),PSPHI)
        CALL SETUP$GETR8A('PSPHIDOT',NR*LNX(ISP),PSPHIDOT)
        CALL SETUP$UNSELECT()
        DO IH=1,LNXH
          IT=POTPAR(ISP)%ITAIL(IH)
          LN=POTPAR(ISP)%LNOFH(IH)
          L=POTPAR(ISP)%LOXH(IH)
          LNDOT=POTPAR(ISP)%LNOFT(IT)  ! PARTIAL WAVE INDEX FOR PHIDOT
          A1=POTPAR(ISP)%KTOPHI(IH)
          A2=POTPAR(ISP)%KTOPHIDOT(IH)
!         == NODESLESS WAVE FUNCTIONS INSERTED ONLY UP TO MATCHING RADIUS
          POTPAR(ISP)%NLCHI(:IRAD-1,IH)=   NLPHI(:IRAD-1,LN)   *A1 &
       &                                +NLPHIDOT(:IRAD-1,LNDOT)*A2
!         ==  ADD DIFFERENCE TO FULL WAVE FUNCTIONS ============================
          POTPAR(ISP)%AECHI(:,IH)=POTPAR(ISP)%NLCHI(:,IH) &
       &                          +(   AEPHI(:,LN)   -   NLPHI(:,LN)   )*A1 &
       &                          +(AEPHIDOT(:,LNDOT)-NLPHIDOT(:,LNDOT))*A2
          POTPAR(ISP)%PSCHI(:,IH)=POTPAR(ISP)%NLCHI(:,IH) &
       &                          +(   PSPHI(:,LN)   -   NLPHI(:,LN)   )*A1 &
       &                          +(PSPHIDOT(:,LNDOT)-NLPHIDOT(:,LNDOT))*A2
        ENDDO
        DEALLOCATE(AEPHI)
        DEALLOCATE(AEPHIDOT)
        DEALLOCATE(PSPHI)
        DEALLOCATE(PSPHIDOT)
        DEALLOCATE(NLPHI)
        DEALLOCATE(NLPHIDOT)
        DEALLOCATE(R)
!
!
!       ========================================================================
!       === WRITE TAILED FUNCTIONS TO FILE                                    ==
!       ========================================================================
        IF(TPR) THEN
          WRITE(STRING,*)ISP
          STRING=ADJUSTL(STRING) 
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_AECHI.DAT',GID,NR &
       &                    ,LNXH,POTPAR(ISP)%AECHI)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_PSCHI.DAT',GID,NR &
       &                    ,LNXH,POTPAR(ISP)%PSCHI)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_NLCHI.DAT',GID,NR &
       &                    ,LNXH,POTPAR(ISP)%NLCHI)
        END IF
      ENDDO ! END OF LOOP OVER ATOM TYPES
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_PIPHI()
!     **************************************************************************
!     ** |PSI>APPROX |CHI>[<CHI|P><PHI|THETA|PHI><P|CHI>]**(-1)               **
!     **                  *<CHI|P><PHI|THETA|PHI><P|PSI>                      **
!     **            =|CHI><PI|PHI><P|PSI>                                     **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2019******************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR &
     &                       ,NSP &
     &                       ,LNX &
     &                       ,LOX &
     &                       ,HYBRIDSETTING 
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4) :: ISP,L1,L2,LN1,LN2,LMN1,LMN2,M
      INTEGER(4) :: LMN1X,LMN2X
      INTEGER(4) :: LN1X,LN2X
      LOGICAL(4),PARAMETER   :: TPR=.TRUE.
      REAL(8)   ,ALLOCATABLE :: AMAT(:,:),AMATINV(:,:)
      REAL(8)   ,ALLOCATABLE :: PHIOV(:,:)
      REAL(8)   ,ALLOCATABLE :: PROK(:,:)
      REAL(8)   ,ALLOCATABLE :: PIPHI(:,:)
!     **************************************************************************
                                           CALL TRACE$PUSH('SIMPLELMTO_PIPHI')
      DO ISP=1,NSP
        LN1X=POTPAR(ISP)%LNXH   !LOCAL ORBITALS
        LN2X=LNX(ISP)             !PARTIAL WAVES
!
!       ========================================================================
!       == <PI|PHI>                                                           ==
!       ========================================================================
        ALLOCATE(PHIOV(LN2X,LN2X))
        ALLOCATE(PROK(LN2X,LN1X))
        ALLOCATE(AMAT(LN1X,LN1X))
        ALLOCATE(AMATINV(LN1X,LN1X))
        ALLOCATE(PIPHI(LN1X,LN2X))      !<PI|CHI>
        PHIOV=POTPAR(ISP)%PHIOV     !<PHI|THETA_OMEGA|PHI>
WRITE(*,*)ISP,LN1X,LN2X,'PHIOV=',POTPAR(ISP)%PHIOV          
        PROK=POTPAR(ISP)%PROK
WRITE(*,*)ISP,'PROK=',POTPAR(ISP)%PROK          
        AMAT=MATMUL(TRANSPOSE(PROK),MATMUL(PHIOV,PROK))
WRITE(*,*)ISP,'AMAT',AMAT
        CALL LIB$INVERTR8(LN1X,AMAT,AMATINV)
WRITE(*,*)ISP,'AMATINV',AMATINV
        PIPHI=MATMUL(AMATINV,MATMUL(TRANSPOSE(PROK),PHIOV))
        DEALLOCATE(PHIOV)
        DEALLOCATE(PROK)
        DEALLOCATE(AMAT)
        DEALLOCATE(AMATINV)
!
!       ========================================================================
!       == EXPAND <PI|PHI>                                                    ==
!       ========================================================================
        LMN1X=SUM(2*POTPAR(ISP)%LOXH(:)+1)
        LMN2X=SUM(2*LOX(:LN2X,ISP)+1)
        ALLOCATE(POTPAR(ISP)%PIPHI(LMN1X,LMN2X))
        POTPAR(ISP)%PIPHI=0.D0

        LMN2=0
        DO LN2=1,LN2X
          L2=LOX(LN2,ISP)
!
          LMN1=0
          DO LN1=1,LN1X
            L1=POTPAR(ISP)%LOXH(LN1)
            IF(L1.EQ.L2) THEN
              DO M=1,2*L1+1
                POTPAR(ISP)%PIPHI(LMN1+M,LMN2+M)=PIPHI(LN1,LN2)
              ENDDO   
            END IF 
            LMN1=LMN1+2*L1+1
          ENDDO
          LMN2=LMN2+2*L2+1
        ENDDO
        DEALLOCATE(PIPHI)

      ENDDO ! END OF LOOP OVER ATOM TYPES
!
!     ==========================================================================
!     == DIAGNOSTIC PRINTOUT                                                  ==
!     ==========================================================================
      IF(TPR) THEN
        DO ISP=1,NSP
          WRITE(*,*)'PIPHI',POTPAR(ISP)%PIPHI          
        ENDDO
!!$        CALL ERROR$MSG('REGULAR STOP AFTER DIAGNOSTIC PRINTOUT')
!!$        CALL ERROR$STOP('SIMPLELMTO_PIPHI')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OVERLAPPHI()
!     **************************************************************************
!     **  OVERLAP BETWEEN ALL-ELECTRON PARTIAL WAVES WITHIN SPHERE            **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR &
     &                             ,NSP &
     &                             ,LNX &
     &                             ,LOX
      IMPLICIT NONE
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
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
        ALLOCATE(POTPAR(ISP)%PHIOV(LNX1,LNX1))
        POTPAR(ISP)%PHIOV=0.D0
        DO LN1=1,LNX1
          DO LN2=LN1,LNX1
            IF(LOX(LN1,ISP).NE.LOX(LN2,ISP)) CYCLE
            CALL RADIAL$INTEGRATE(GID,NR,R**2*AEPHI(:,LN1)*AEPHI(:,LN2),AUX)
            CALL RADIAL$VALUE(GID,NR,AUX,POTPAR(ISP)%RAUG,VAL)
            POTPAR(ISP)%PHIOV(LN1,LN2)=VAL
            POTPAR(ISP)%PHIOV(LN2,LN1)=VAL
          ENDDO
        ENDDO
        DEALLOCATE(R)
        DEALLOCATE(AUX)
        DEALLOCATE(AEPHI)
        CALL SETUP$UNSELECT()
!
!       ========================================================================
!       ==  DIAGNOSTIC PRINTOUT                                               ==
!       ========================================================================
        IF(TPR) THEN
          WRITE(*,FMT='(82("="),T10," LOCAL ORBITALS OVERLAP FOR ISP=",I3)')ISP
          DO LN1=1,LNX1
            WRITE(*,FMT='(20F10.5)')POTPAR(ISP)%PHIOV(LN1,:)
          ENDDO
          CALL ERROR$MSG('REGULAR STOP AFTER DIAGNOSTIC PRINTOUT') 
          CALL ERROR$STOP('SIMPLELMTO_OVERLAPPHI')
        END IF
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_ONSITEMATRIXELEMENTS()
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR &
     &                             ,NSP &
     &                             ,SCREENL
      IMPLICIT NONE
      INTEGER(4)             :: ISP
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4)             :: LMNX
      INTEGER(4)             :: LRX
      REAL(8)   ,ALLOCATABLE :: ULITTLE(:,:,:,:,:) ! SLATER INTEGRALS
!     **************************************************************************
      DO ISP=1,NSP
        GID=POTPAR(ISP)%GID
        CALL RADIAL$GETI4(GID,'NR',NR)
        LNX=POTPAR(ISP)%LNXH
        ALLOCATE(LOX(LNX))
        LOX=POTPAR(ISP)%LOXH
        LMNX=SUM(2*LOX+1)
        LRX=2*MAXVAL(LOX)
!
        ALLOCATE(ULITTLE(LRX+1,LNX,LNX,LNX,LNX))
        CALL SIMPLELMTO_ULITTLE(GID,NR,LRX,LNX,LOX,POTPAR(ISP)%AECHI &
     &                         ,SCREENL,ULITTLE)
        IF(.NOT.ALLOCATED(POTPAR(ISP)%ONSITEU)) THEN
          ALLOCATE(POTPAR(ISP)%ONSITEU(LMNX,LMNX,LMNX,LMNX))
        END IF
        CALL SIMPLELMTO_UTENSOR(LRX,LMNX,LNX,LOX,ULITTLE,POTPAR(ISP)%ONSITEU)
        DEALLOCATE(ULITTLE)

        IF(.NOT.ALLOCATED(POTPAR(ISP)%OVERLAP)) THEN
          ALLOCATE(POTPAR(ISP)%OVERLAP(LMNX,LMNX))
        END IF
        CALL SIMPLELMTO_ONECENTEROVERLAP(GID,NR,LNX,LOX,POTPAR(ISP)%AECHI &
     &                                  ,LMNX,POTPAR(ISP)%OVERLAP)
!
!       ========================================================================
        IF(.NOT.ALLOCATED(POTPAR(ISP)%QLN)) THEN
          ALLOCATE(POTPAR(ISP)%QLN(2,LNX,LNX))
        END IF
        CALL SIMPLELMTO_ONECENTERQLN(GID,NR,LNX,LOX,POTPAR(ISP)%AECHI &
     &                            ,POTPAR(ISP)%QLN)
        DEALLOCATE(LOX)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_ULITTLE(GID,NR,LRX,LNX,LOX,CHI,SCREENL,ULITTLE)
!     **************************************************************************
!     ** SLATER INTEGRALS.                                                    **
!     **                                                                      **
!     ** A NEGATIVE VALUE FOR THE SCREENING LENTH SCREENL IS INTERPRETED AS   **
!     ** NO SCREENING                                                         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LRX
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: CHI(NR,LNX)
      REAL(8)   ,INTENT(IN) :: SCREENL ! SCREENING LENGTH FOR THE INTERACTION
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
                            CALL TRACE$PUSH('SIMPLELMTO_ULITTLE')
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
            IF(SCREENL.LT.0.D0) THEN
              CALL RADIAL$POISSON(GID,NR,L,RHO,POT)
            ELSE
              CALL RADIAL$YUKAWA(GID,NR,L,1.D0/SCREENL,RHO,POT) 
           END IF
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
      SUBROUTINE SIMPLELMTO_UTENSOR(LRX,NORB,LNX,LOX,ULITTLE,U)
!     **************************************************************************
!     ** EXPANDS SLATER INTEGRALS FROM SIMPLELMTO_ULITTLE TO THE FULL U-TENSOR      **
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
                            CALL TRACE$PUSH('SIMPLELMTO_UTENSOR')
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
                          CALL SPHERICAL$GAUNT(LM2,LM4,LM,CG1)
                          CALL SPHERICAL$GAUNT(LM3,LM1,LM,CG2)
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
      SUBROUTINE SIMPLELMTO_ONECENTEROVERLAP(GID,NR,LNX,LOX,CHI,LMNX,OVERLAP)
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
                            CALL TRACE$PUSH('SIMPLELMTO_ONECENTEROVERLAP')
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
      SUBROUTINE SIMPLELMTO_ONECENTERQLN(GID,NR,LNX,LOX,CHI,QLN)
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
                            CALL TRACE$PUSH('SIMPLELMTO_ONECENTERMULTIPOLE')
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
      SUBROUTINE SIMPLELMTO$ETOT()
!     **************************************************************************
!     **  DENMAT_ ON INPUT IS CALCULATED DIRECTLY FROM THE PROJECTIONS AND    **
!     **  IS USED IN THE AUGMENTATION                                         **
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : TON &
     &                             ,MODUS
      IMPLICIT NONE
      INTEGER(4)            :: IAT
!     **************************************************************************
      CALL SIMPLELMTO$INITIALIZE()
      IF(.NOT.TON) RETURN
                                    CALL TRACE$PUSH('SIMPLELMTO$ETOT')
 
!!$DO  IAT=1,NAT_
!!$  WRITE(*,FMT='("DONSITE 1: ",20F10.5)')DENMAT_(:,:,:,IAT)
!!$ENDDO
!
!     ==========================================================================
!     ==  SELECT CHOICES                                                      ==
!     ==========================================================================
      WRITE(*,FMT='(82("="),T30," SIMPLELMTO$ENERGY START. MODUS=",A," ")') &
     &        TRIM(MODUS)
      IF(MODUS.EQ.'DMFT') THEN
        CALL DMFT$GREEN()
      ELSE IF(MODUS.EQ.'HYBRID') THEN
        CALL SIMPLELMTO_HYBRID()
      ELSE
        CALL ERROR$MSG('MODUS NOT RECOGNIZED')
        CALL ERROR$MSG('ALLOWED VALUES ARE "DMFT" AND "HYBRID"')
        CALL ERROR$CHVAL('MODUS',MODUS)
        CALL ERROR$STOP('SIMPLELMTO$ETOT')
      END IF
!
      WRITE(*,FMT='(82("="),T30," LMTO$ENERGY DONE ")')
                                    CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_HYBRID()
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : TOFFSITE
      USE RSPACEOP_MODULE, ONLY : RSPACEMAT_TYPE &
     &                           ,RSPACEOP$DELETE &
     &                           ,RSPACEOP$COPY &
     &                           ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      LOGICAL(4)          ,PARAMETER   :: TPR=.FALSE.
      LOGICAL(4)          ,PARAMETER   :: TTEST=.FALSE.
      LOGICAL(4)          ,PARAMETER   :: TTESTA=.FALSE.
      INTEGER(4)                       :: NAT
      INTEGER(4)                       :: NND
      REAL(8)                          :: RBAS(3,3)
      REAL(8)                          :: STRESS(3,3)
      REAL(8)             ,ALLOCATABLE :: R0(:,:) !(3,NAT)
      REAL(8)             ,ALLOCATABLE :: FORCE(:,:) !(3,NAT)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: DENMAT(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: HAMIL(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: DONSITE(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: HONSITE(:)
      REAL(8)             ,ALLOCATABLE :: FORCET(:,:) !(3,NAT)
      REAL(8)                          :: STRESST(3,3)
      REAL(8)                          :: ETOT
      INTEGER(4)                       :: I,IAT
      INTEGER(4)                       :: NFILO
      LOGICAL(4)                       :: TFIRST=.TRUE.
!     **************************************************************************
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      ALLOCATE(FORCE(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!
!     ==========================================================================
!     ==  COLLECT DENSITY MATRIX                                              ==
!     ==========================================================================
      CALL WAVES$GETI4('NND',NND)
PRINT*,'NND ',NND
      ALLOCATE(DENMAT(NND))
      CALL WAVES$GETRSPACEMATA('DENMAT',NND,DENMAT)
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'DENSITY MATRIX IN PHIS',NND,DENMAT)
!
!     ==========================================================================
!     ==  EXTRACT ON-SITE DENSITY MATRIX IN PARTIAL WAVE EXPANSION            ==
!     ==========================================================================
      ALLOCATE(DONSITE(NAT))
      CALL SIMPLELMTO_ONSITEDENMATEXTRACTADD('EXTRACT',NND,DENMAT,NAT,DONSITE)
!!$DO  IAT=1,NAT
!!$  WRITE(*,FMT='("DONSITE 2: ",20F10.5)')DONSITE(IAT)%MAT
!!$ENDDO
!!$STOP 'FORCED'
!
!     ==========================================================================
!     ==  CONVERT DENSITY MATRIX INTO THE LOCAL-ORBITAL BASIS                 ==
!     ==========================================================================
      CALL SIMPLELMTO_DENMATSHRINKEXPAND('SHRINK',NND,DENMAT) 
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'DENSITY MATRIX IN CHIS',NND,DENMAT)
!
!     ==========================================================================
!     ==  CALCULATE ENERGY                                                    ==
!     ==========================================================================
      ETOT=0.D0
      ALLOCATE(HAMIL(NND))
      DO I=1,NND
        CALL RSPACEOP$COPY(DENMAT(I),HAMIL(I))
        HAMIL(I)%MAT=0.D0
      ENDDO
      ALLOCATE(HONSITE(NAT))
      DO IAT=1,NAT
        CALL RSPACEOP$COPY(DONSITE(IAT),HONSITE(IAT))
        HONSITE(IAT)%MAT=0.D0
      ENDDO
!
IF(.NOT.TTEST.OR.TFIRST) THEN
TFIRST=.FALSE.
      CALL TIMING$CLOCKON('SIMPLELMTO_HYBRIDENERGY')
      CALL SIMPLELMTO_HYBRIDENERGY(NAT,NND,RBAS,R0,DENMAT,DONSITE &
     &                                  ,ETOT,STRESS,FORCE,HAMIL,HONSITE)
      CALL TIMING$CLOCKOFF('SIMPLELMTO_HYBRIDENERGY')
END IF
!
!     ==========================================================================
!     == 
!     ==========================================================================
IF(TTEST.AND.TTESTA) THEN
  CALL FILEHANDLER$UNIT('PROT',NFILO)
  WRITE(NFILO,*)'WARNING FROM SIMPLELMTO_HYBRID! TEST ENVIRONMENT IS ON'
  DO I=1,NND !NEIGBORLIST IS DIRECTIONAL
    CALL RSPACEOP$DELETE(DENMAT(I))
  ENDDO
  CALL WAVES$GETRSPACEMATA('DENMAT',NND,DENMAT)
  CALL SIMPLELMTO_DENMATSHRINKEXPAND('SHRINK',NND,DENMAT) 
  CALL SIMPLELMTO_FAKEETOT(NND,DENMAT,ETOT,HAMIL)
  FORCE=0.D0
  STRESS=0.D0
  DO IAT=1,NAT
    HONSITE(IAT)%MAT=0.D0
  ENDDO
END IF
!
!     ==========================================================================
!     ==  CONVERT HAMIL INTO THE PARTIAL-WAVE EXPANSION                       ==
!     ==========================================================================
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'HAMILTON IN CHI-S',NND,HAMIL)
      CALL SIMPLELMTO_DENMATSHRINKEXPAND('EXPAND',NND,HAMIL)
!
!     ==========================================================================
!     ==  ADD ONSITE HAMILTONIAN TO OFF-SITE HAMILTONIAN                      ==
!     ==========================================================================
      CALL SIMPLELMTO_ONSITEDENMATEXTRACTADD('ADDBACK',NND,HAMIL,NAT,HONSITE)
!!$DO  IAT=1,NAT
!!$  WRITE(*,FMT='("HONSITE 2: ",20F10.5)')HONSITE(IAT)%MAT
!!$ENDDO
IF(TTEST.AND.(.NOT.TTESTA)) THEN
  CALL FILEHANDLER$UNIT('PROT',NFILO)
  WRITE(NFILO,*)'WARNING FROM SIMPLELMTO_HYBRID! TEST ENVIRONMENT IS ON'
  DO I=1,NND !NEIGBORLIST IS DIRECTIONAL
    CALL RSPACEOP$DELETE(DENMAT(I))
  ENDDO
  CALL WAVES$GETRSPACEMATA('DENMAT',NND,DENMAT)
  CALL SIMPLELMTO_FAKEETOT(NND,DENMAT,ETOT,HAMIL)
  FORCE=0.D0
  STRESS=0.D0
  DO IAT=1,NAT
    HONSITE(IAT)%MAT=0.D0
  ENDDO
END IF
!
!     ==========================================================================
!     == COMMUNICATE ENERGY AND DERIVATIVES                                   ==
!     ==========================================================================
!     == ENERGY -> ENERGYLIST ==================================================
      IF(TPR)PRINT*,'ENERGY FROM SIMPLELMTO INTERFACE ',ETOT
      CALL ENERGYLIST$SET('LMTO INTERFACE',ETOT)
      CALL ENERGYLIST$ADD('LOCAL CORRELATION',ETOT)
      CALL ENERGYLIST$ADD('TOTAL ENERGY',ETOT)
!
!     == FORCE -> ATOMLIST =====================================================
      ALLOCATE(FORCET(3,NAT))
      CALL ATOMLIST$GETR8A('FORCE',0,3*NAT,FORCET)
      FORCET=FORCET+FORCE
      CALL ATOMLIST$SETR8A('FORCE',0,3*NAT,FORCET)
      DEALLOCATE(FORCET)
!
!     == STRESS -> CELL ========================================================
      CALL CELL$GETR8A('STRESS_I',9,STRESST)
      STRESST=STRESST-STRESS  ! IN THE PRESENT ROUTINE STRESS=+DE/DEPSILON!
      CALL CELL$SETR8A('STRESS_I',9,STRESST)
!
!     ==  HAMILTONIAN -> WAVES OBJECT ==========================================
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'HAMILTON IN PHIS',NND,HAMIL)
      CALL WAVES$SETRSPACEMATA('HAMIL',NND,HAMIL)
!
!     ==========================================================================
!     ==  CLEAN DENMAT AND HAMIL                                              ==
!     ==========================================================================
      DO I=1,NND
        CALL RSPACEOP$DELETE(DENMAT(I))
        CALL RSPACEOP$DELETE(HAMIL(I))
      ENDDO
      DEALLOCATE(DENMAT)
      DEALLOCATE(HAMIL)
      DO IAT=1,NAT
        CALL RSPACEOP$DELETE(DONSITE(IAT))
        CALL RSPACEOP$DELETE(HONSITE(IAT))
      ENDDO
      DEALLOCATE(DONSITE)
      DEALLOCATE(HONSITE)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_FAKEETOT(NND,DENMAT,ETOT,HAMIL)
!     **************************************************************************
!     **                                                                      **
!     ********************************************P. BLOECHL, GOSLAR 2019*******
      USE RSPACEOP_MODULE, ONLY : RSPACEMAT_TYPE &
     &                           ,RSPACEOP$DELETE &
     &                           ,RSPACEOP$COPY &
     &                           ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN)   :: NND
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: DENMAT(NND)
      REAL(8)             ,INTENT(INOUT):: ETOT
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: HAMIL(NND)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE,SAVE:: HAMILSAVE(:)
      REAL(8)             ,SAVE         :: ETOTSAVE
      INTEGER(4)                        :: IND,J
      INTEGER(4)                        :: IND1,IND2
      REAL(8)                           :: SVAR
      LOGICAL(4)                        :: TCHK
      LOGICAL(4)          ,SAVE         :: TINI=.FALSE.
!     **************************************************************************
!
!     ==========================================================================
!     == STORE ENERGY AND HAMILTONIAN IN THE FIRST CALL
!     ==========================================================================
      IF(.NOT.TINI) THEN
        TINI=.TRUE.
        ALLOCATE(HAMILSAVE(NND))
        DO IND=1,NND
          CALL RSPACEOP$COPY(HAMIL(IND),HAMILSAVE(IND))
        ENDDO
        ETOTSAVE=ETOT
        DO IND=1,NND
          DO J=1,HAMILSAVE(IND)%N3
            ETOTSAVE=ETOTSAVE-SUM(DENMAT(IND)%MAT(:,:,J)*HAMIL(IND)%MAT(:,:,J))
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == OVERWRITE ENERGY AND HAMILTONIAN
!     ==========================================================================
      ETOT=ETOTSAVE
      DO IND=1,NND
!PRINT*,'==',HAMIL(IND)%IAT1,HAMIL(IND)%IAT2,HAMIL(IND)%IT
        CALL RSPACEOP$COPY(HAMILSAVE(IND),HAMIL(IND))
        DO J=1,HAMILSAVE(IND)%N3
          ETOT=ETOT+SUM(DENMAT(IND)%MAT(:,:,J)*HAMIL(IND)%MAT(:,:,J))
        ENDDO
!HAMIL(IND)%MAT=0.D0
      ENDDO
!      ETOT=2.D0*ETOT
!ETOT=0.D0

!
!     ==========================================================================
!     == CHECK WHETHER DENSITY MATRIX IS HERMITEAN
!     ==========================================================================
      DO IND1=1,NND
        TCHK=.FALSE.
        DO IND2=1,NND
          IF(DENMAT(IND2)%IAT2.NE.DENMAT(IND1)%IAT1) CYCLE
          IF(DENMAT(IND2)%IAT1.NE.DENMAT(IND1)%IAT2) CYCLE
          IF(SUM((DENMAT(IND2)%IT+DENMAT(IND1)%IT)**2).NE.0) CYCLE
          IF(TCHK) THEN
            CALL ERROR$MSG('ERROR 4')
            CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
          END IF
          TCHK=.TRUE.
          SVAR=0.D0
          DO J=1,DENMAT(IND1)%N3
            SVAR=SVAR+SUM((DENMAT(IND2)%MAT(:,:,J) &
     &          -TRANSPOSE(DENMAT(IND1)%MAT(:,:,J)))**2)
          ENDDO
          IF(SVAR.GT.1.D-10) THEN
            CALL RSPACEOP$WRITEMAT(6,'DENMAT',NND,DENMAT)
            CALL ERROR$MSG('ERROR 5')
            CALL ERROR$I4VAL('IND1',IND1)
            CALL ERROR$I4VAL('IND2',IND2)
            CALL ERROR$I4VAL('IND1-IAT1',DENMAT(IND1)%IAT1)
            CALL ERROR$I4VAL('IND1-IAT2',DENMAT(IND1)%IAT2)
            CALL ERROR$I4VAL('IND1-IT',DENMAT(IND1)%IT)
            CALL ERROR$I4VAL('IND2-IAT1',DENMAT(IND2)%IAT1)
            CALL ERROR$I4VAL('IND2-IAT2',DENMAT(IND2)%IAT2)
            CALL ERROR$I4VAL('IND2-IT',DENMAT(IND2)%IT)
            CALL ERROR$R8VAL('DEV',SVAR)
            CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
          END IF
        ENDDO
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('ERROR 6')
          CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
        END IF
      ENDDO
!
!     ==========================================================================
!     == MAKE HAMILTONIAN  HERMITEAN
!     ==========================================================================
      DO IND1=1,NND
        TCHK=.FALSE.
        DO IND2=1,NND
          IF(HAMIL(IND2)%IAT2.NE.HAMIL(IND1)%IAT1) CYCLE
          IF(HAMIL(IND2)%IAT1.NE.HAMIL(IND1)%IAT2) CYCLE
          IF(SUM((HAMIL(IND2)%IT+HAMIL(IND1)%IT)**2).NE.0) CYCLE
          IF(TCHK) THEN
            CALL ERROR$MSG('ERROR 1')
            CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
          END IF
          TCHK=.TRUE.
          SVAR=0.D0
          DO J=1,HAMIL(IND1)%N3
            HAMIL(IND2)%MAT(:,:,J)=0.5D0*(HAMIL(IND2)%MAT(:,:,J) &
     &                         +TRANSPOSE(HAMIL(IND1)%MAT(:,:,J)))
            HAMIL(IND1)%MAT(:,:,J)=TRANSPOSE(HAMIL(IND2)%MAT(:,:,J))
          ENDDO
          IF(SVAR.GT.1.D-10) THEN
            CALL RSPACEOP$WRITEMAT(6,'HAMIL',NND,HAMIL)
            CALL ERROR$MSG('ERROR 2')
            CALL ERROR$I4VAL('IND1',IND1)
            CALL ERROR$I4VAL('IND2',IND2)
            CALL ERROR$I4VAL('IND1-IAT1',HAMIL(IND1)%IAT1)
            CALL ERROR$I4VAL('IND1-IAT2',HAMIL(IND1)%IAT2)
            CALL ERROR$I4VAL('IND1-IT',HAMIL(IND1)%IT)
            CALL ERROR$I4VAL('IND2-IAT1',HAMIL(IND2)%IAT1)
            CALL ERROR$I4VAL('IND2-IAT2',HAMIL(IND2)%IAT2)
            CALL ERROR$I4VAL('IND2-IT',HAMIL(IND2)%IT)
            CALL ERROR$R8VAL('DEV',SVAR)
            CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
          END IF
        ENDDO
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('ERROR 3')
          CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
        END IF
      ENDDO
!
!     ==========================================================================
!     == CHECK WHETHER HAMILTONIAN IS HERMITEAN
!     ==========================================================================
      DO IND1=1,NND
        TCHK=.FALSE.
        DO IND2=1,NND
          IF(HAMIL(IND2)%IAT2.NE.HAMIL(IND1)%IAT1) CYCLE
          IF(HAMIL(IND2)%IAT1.NE.HAMIL(IND1)%IAT2) CYCLE
          IF(SUM((HAMIL(IND2)%IT+HAMIL(IND1)%IT)**2).NE.0) CYCLE
          IF(TCHK) THEN
            CALL ERROR$MSG('ERROR 1')
            CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
          END IF
          TCHK=.TRUE.
          SVAR=0.D0
          DO J=1,HAMIL(IND1)%N3
            SVAR=SVAR+SUM((HAMIL(IND2)%MAT(:,:,J) &
     &          -TRANSPOSE(HAMIL(IND1)%MAT(:,:,J)))**2)
          ENDDO
          IF(SVAR.GT.1.D-10) THEN
            CALL RSPACEOP$WRITEMAT(6,'HAMIL',NND,HAMIL)
            CALL ERROR$MSG('ERROR 2')
            CALL ERROR$I4VAL('IND1',IND1)
            CALL ERROR$I4VAL('IND2',IND2)
            CALL ERROR$I4VAL('IND1-IAT1',HAMIL(IND1)%IAT1)
            CALL ERROR$I4VAL('IND1-IAT2',HAMIL(IND1)%IAT2)
            CALL ERROR$I4VAL('IND1-IT',HAMIL(IND1)%IT)
            CALL ERROR$I4VAL('IND2-IAT1',HAMIL(IND2)%IAT1)
            CALL ERROR$I4VAL('IND2-IAT2',HAMIL(IND2)%IAT2)
            CALL ERROR$I4VAL('IND2-IT',HAMIL(IND2)%IT)
            CALL ERROR$R8VAL('DEV',SVAR)
            CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
          END IF
        ENDDO
        IF(.NOT.TCHK) THEN
          CALL ERROR$MSG('ERROR 3')
          CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
        END IF
      ENDDO
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_ONSITEDENMATEXTRACTADD(ID,NND,DENMAT,NAT,DONSITE)
!     **************************************************************************
!     ** COLLECTS THE ONSITE TERMS OFG THE DENSITY MATRIX                     **
!     **                                                                      **
!     ********************************************P. BLOECHL, GOSLAR 2019*******
      USE SIMPLELMTO_MODULE, ONLY : ISPECIES &
     &                             ,POTPAR &
     &                             ,LOX &
     &                             ,LNX
      USE RSPACEOP_MODULE, ONLY : RSPACEMAT_TYPE &
     &                           ,RSPACEOP$DELETE &
     &                           ,RSPACEOP$COPY &
     &                           ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      CHARACTER(*)        ,INTENT(IN)   :: ID
      INTEGER(4)          ,INTENT(IN)   :: NAT
      INTEGER(4)          ,INTENT(IN)   :: NND
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: DENMAT(NND)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: DONSITE(NAT)
      INTEGER(4)                        :: I,IAT
!     ********************************************P. BLOECHL, GOSLAR 2019*******
      IF(ID.EQ.'EXTRACT') THEN
        DO I=1,NND
          IF(DENMAT(I)%IAT1.NE.DENMAT(I)%IAT2) CYCLE
          IF(SUM(DENMAT(I)%IT**2).NE.0) CYCLE
          IAT=DENMAT(I)%IAT1
          CALL RSPACEOP$COPY(DENMAT(I),DONSITE(IAT))
        ENDDO
      ELSE IF(ID.EQ.'ADDBACK') THEN
        DO I=1,NND
          IF(DENMAT(I)%IAT1.NE.DENMAT(I)%IAT2) CYCLE
          IF(SUM(DENMAT(I)%IT**2).NE.0) CYCLE
          IAT=DENMAT(I)%IAT1
          IF(DENMAT(I)%N1.NE.DONSITE(IAT)%N1.OR. &
     &       DENMAT(I)%N2.NE.DONSITE(IAT)%N2.OR. &
     &       DENMAT(I)%N3.NE.DONSITE(IAT)%N3) THEN
            CALL ERROR$MSG('INCONSISTENT ARRAY SIZES')
            CALL ERROR$CHVAL('ID',ID)
            CALL ERROR$I4VAL('ONSITE N1',DONSITE(IAT)%N1)
            CALL ERROR$I4VAL('ONSITE N2',DONSITE(IAT)%N2)
            CALL ERROR$I4VAL('ONSITE N3',DONSITE(IAT)%N3)
            CALL ERROR$I4VAL('OFFSITE N1',DENMAT(I)%N1)
            CALL ERROR$I4VAL('OFFSITE N2',DENMAT(I)%N2)
            CALL ERROR$I4VAL('OFFSITE N3',DENMAT(I)%N3)
            CALL ERROR$STOP('SIMPLELMTO_ONSITEDENMATEXTRACTADD')
          END IF
          DENMAT(I)%MAT=DENMAT(I)%MAT+DONSITE(IAT)%MAT
        ENDDO
      ELSE
        CALL ERROR$MSG('ILLEGAL VALUE OF ID')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SIMPLELMTO_ONSITEDENMATEXTRACTADD')
      END IF
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_DENMATSHRINKEXPAND(ID,NND,DENMAT)
!     **************************************************************************
!     ** TRANSFORMS THE DENSITY MATRIX FROM A PARTIAL-WAVE REPRESENTATION     **
!     ** TO A LOCAL ORBITAL REPRESENTATION (SHRINK) AND IT                    **
!     ** TRANSFORMES THE HAMILTONIAN FROM A LOCAL-ORBITAL REPRESENTATION      **
!     ** TO A PARTIAL-WAVE REPRESENTATION (EXPAND).                           **
!     **                                                                      **
!     **            |PSI> APPROX |CHI><PI|PHI><PTILDE|PSITILDE>               **
!     **                                                                      **
!     ** LNX,LOX REFER TO THE PARTIAL-WAVE REPRESENTATION                     **
!     **                                                                      **
!     ********************************************P. BLOECHL, GOSLAR 2019*******
      USE SIMPLELMTO_MODULE, ONLY : ISPECIES &
     &                             ,POTPAR &
     &                             ,LOX &
     &                             ,LNX
      USE RSPACEOP_MODULE, ONLY : RSPACEMAT_TYPE &
     &                           ,RSPACEOP$DELETE &
     &                           ,RSPACEOP$COPY &
     &                           ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      CHARACTER(*)        ,INTENT(IN)   :: ID
      INTEGER(4)          ,INTENT(IN)   :: NND
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: DENMAT(NND)
      REAL(8)             ,ALLOCATABLE  :: MAT(:,:,:)
      INTEGER(4)                        :: IAT1,IAT2
      INTEGER(4)                        :: ISP1,ISP2
      INTEGER(4)                        :: NCHI1,NCHI2
      INTEGER(4)                        :: NPHI1,NPHI2
      INTEGER(4)                        :: N3
      INTEGER(4)                        :: I,J
!     **************************************************************************
      IF(ID.EQ.'SHRINK') THEN
        DO I=1,NND
          IAT1=DENMAT(I)%IAT1
          IAT2=DENMAT(I)%IAT2
          ISP1=ISPECIES(IAT1)
          ISP2=ISPECIES(IAT2)
          NCHI1=SUM(2*POTPAR(ISP1)%LOXH+1)
          NPHI1=SUM(2*LOX(:LNX(ISP1),ISP1)+1)
          NCHI2=SUM(2*POTPAR(ISP2)%LOXH+1)
          NPHI2=SUM(2*LOX(:LNX(ISP2),ISP2)+1)
          N3=DENMAT(I)%N3
          ALLOCATE(MAT(NCHI1,NCHI2,N3))
          DO J=1,N3
            MAT(:,:,J)=MATMUL(POTPAR(ISP1)%PIPHI &
     &             ,MATMUL(DENMAT(I)%MAT(:,:,J),TRANSPOSE(POTPAR(ISP2)%PIPHI)))
          ENDDO
          DEALLOCATE(DENMAT(I)%MAT)
          ALLOCATE(DENMAT(I)%MAT(NCHI1,NCHI2,N3))
          DENMAT(I)%MAT=MAT
          DEALLOCATE(MAT)
          DENMAT(I)%N1=NCHI1
          DENMAT(I)%N2=NCHI2
        ENDDO
      ELSE IF(ID.EQ.'EXPAND') THEN
        DO I=1,NND
          IAT1=DENMAT(I)%IAT1
          IAT2=DENMAT(I)%IAT2
          ISP1=ISPECIES(IAT1)
          ISP2=ISPECIES(IAT2)
          NCHI1=SUM(2*POTPAR(ISP1)%LOXH+1)
          NPHI1=SUM(2*LOX(:LNX(ISP1),ISP1)+1)
          NCHI2=SUM(2*POTPAR(ISP2)%LOXH+1)
          NPHI2=SUM(2*LOX(:LNX(ISP2),ISP2)+1)
          N3=DENMAT(I)%N3
          ALLOCATE(MAT(NCHI1,NCHI2,N3))
          MAT=DENMAT(I)%MAT
          DEALLOCATE(DENMAT(I)%MAT)
          ALLOCATE(DENMAT(I)%MAT(NPHI1,NPHI2,N3))
          DO J=1,N3
            DENMAT(I)%MAT(:,:,J)=MATMUL(TRANSPOSE(POTPAR(ISP1)%PIPHI) &
      &                                ,MATMUL(MAT(:,:,J),POTPAR(ISP2)%PIPHI))
          ENDDO
          DEALLOCATE(MAT)
          DENMAT(I)%N1=NPHI1
          DENMAT(I)%N2=NPHI2
        ENDDO
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$MSG('ID MUST BE EITHER "SHRINK" OR "EXPAND"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('SIMPLELMTO_DENMATSHRINKEXPAND')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_HYBRIDENERGY(NAT,NND,RBAS,R0,DENMAT,DONSITE &
     &                                  ,ETOT,STRESS,FORCE,HAMIL,HONSITE)
!     **************************************************************************
!     **  WORK OUT THE ENERGY USING THE LOCAL APPROXIMATION                   **
!     **  TAILED PARTIAL WAVES                                                **
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY: ISPECIES &
     &                            ,LNXPHI=>LNX &
     &                            ,LOXPHI=>LOX &
     &                            ,POTPAR &
     &                            ,TOFFSITE &
     &                            ,HYBRIDSETTING &
     &                            ,HFWEIGHT
      USE MPE_MODULE
      USE RSPACEOP_MODULE, ONLY : RSPACEMAT_TYPE &
     &                           ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN)   :: NND
      INTEGER(4)          ,INTENT(IN)   :: NAT
      REAL(8)             ,INTENT(IN)   :: RBAS(3,3)
      REAL(8)             ,INTENT(IN)   :: R0(3,NAT)
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: DENMAT(NND)
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: DONSITE(NAT)
      REAL(8)             ,INTENT(OUT)  :: ETOT
      REAL(8)             ,INTENT(OUT)  :: STRESS(3,3)
      REAL(8)             ,INTENT(OUT)  :: FORCE(3,NAT)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: HAMIL(NND)   !INTENT(OUT)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: HONSITE(NAT) !INTENT(OUT)
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      INTEGER(4)            :: IND
      INTEGER(4)            :: LNX
      INTEGER(4),ALLOCATABLE:: LOX(:)
      INTEGER(4)            :: LMNX
      INTEGER(4)            :: LMNXPHI
      INTEGER(4)            :: NDIMD
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: LMRX,LRX
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      REAL(8)   ,ALLOCATABLE:: AECORE(:)
      REAL(8)   ,ALLOCATABLE:: U(:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: D(:,:,:)
      REAL(8)   ,ALLOCATABLE:: DON(:,:,:)
      REAL(8)   ,ALLOCATABLE:: H(:,:,:)
      REAL(8)   ,ALLOCATABLE:: HON(:,:,:)
      REAL(8)               :: EH,EX,Q,EAT
      INTEGER(4)            :: NN,IAT,I,J,K,L,IS,IAT1,IAT2,ISP
      INTEGER(4)            :: LMN,LN,IM
      REAL(8)               :: QSPIN(4)
      REAL(8)               :: SVAR
      REAL(8)               :: LHFWEIGHT
      INTEGER(4)            :: IDFTTYPE
!REAL(8)   ,ALLOCATABLE:: T(:,:)
      REAL(8)               :: HFSCALE
      LOGICAL(4)            :: TACTIVE ! DOES THIS ATOM CONTRIBUTE?
      LOGICAL(4),PARAMETER  :: TPARALLEL=.FALSE.
      INTEGER(4)            :: THISTASK,NTASKS
!     **************************************************************************
                            CALL TRACE$PUSH('SIMPLELMTO_HYBRIDENERGY')
      IF(TPR)WRITE(*,FMT='(82("="),T10,"  ",A,"  ")')'SIMPLELMTO_HYBRIDENERGY'
!
!     ==========================================================================
!     == PARALLELIZATION MODEL                                                ==
!     == THE OUTCOME OF THIS ROUTINE IS                                       ==
!     == -- HAMIL (SIMPLELMTO_MODULE; HAMIL_NEW)                              ==
!     == -- ETOT (PASSED TO ENERGLIST)                                        ==
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
      ETOT=0.D0
      DO I=1,NND
        HAMIL(I)%MAT=0.D0
      ENDDO
      DO IAT=1,NAT
        HONSITE(IAT)%MAT=0.D0
      ENDDO
      FORCE=0.D0
      STRESS=0.D0
      DO IAT=1,NAT
        IF(TPARALLEL.AND.MOD(IAT-1,NTASKS).NE.THISTASK-1) CYCLE
!
        ISP=ISPECIES(IAT)
        TACTIVE=(POTPAR(ISP)%LNXH.GT.0)
        IF(.NOT.TACTIVE) CYCLE
!
!       ========================================================================
!       == COLLECT INFORMATION FROM SETUPS OBJECT ==============================
!       ========================================================================
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(AECORE(NR))
        CALL SETUP$GETR8A('AECORE',NR,AECORE)
        CALL SETUP$GETI4('LMRX',LMRX)
        CALL SETUP$UNSELECT()
        LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
!
!       ========================================================================
!       == ADJUSTMENT IF LOCAL HFWEIGHT IS DIFFERENT FROM GLOBAL HFWEIGHT ======
!       ========================================================================
        LHFWEIGHT=HYBRIDSETTING(ISP)%LHFWEIGHT
        IF(LHFWEIGHT.LT.0.D0) LHFWEIGHT=HFWEIGHT
        HFSCALE=1.D0
        IF(HFWEIGHT.GT.0.D0)HFSCALE=LHFWEIGHT/HFWEIGHT
        IF(TPR)PRINT*,'LHFW=',LHFWEIGHT,' GHFW=',HFWEIGHT
!
!       ========================================================================
!       == FIND INDEX TO THE ONSITE ELEMENTS OF THE DENSITY MATRIX
!       ========================================================================
        CALL SIMPLELMTO_INDEXLOCAL2(IAT,NND,DENMAT,IND)
        LMNX=DENMAT(IND)%N1
        NDIMD=DENMAT(IND)%N3
        LNX=SIZE(POTPAR(ISP)%LOXH)
        ALLOCATE(LOX(LNX))
        LOX=POTPAR(ISP)%LOXH
!
!       ========================================================================
!       == CALCULATE U-TENSOR                                                 ==
!       ========================================================================
        U=POTPAR(ISP)%ONSITEU
!PRINT*,'U ',U
!PRINT*,'O ',POTPAR(ISP)%OVERLAP
!
!       ========================================================================
!       == 
!       ========================================================================
        ALLOCATE(D(LMNX,LMNX,NDIMD))
        ALLOCATE(H(LMNX,LMNX,NDIMD))
        LMNXPHI=DONSITE(IAT)%N1
        ALLOCATE(DON(LMNXPHI,LMNXPHI,NDIMD))
        ALLOCATE(HON(LMNXPHI,LMNXPHI,NDIMD))
        DON=DONSITE(IAT)%MAT
        D=DENMAT(IND)%MAT
!
!       ========================================================================
!       == WORK OUT TOTAL ENERGY AND HAMILTONIAN                              ==
!       ========================================================================
        EX=0.D0
        EAT=0.D0
        QSPIN=0.D0
        H(:,:,:)=0.D0
        HON(:,:,:)=0.D0
        DO I=1,LMNX
          DO J=1,LMNX
            QSPIN(:NDIMD)=QSPIN(:NDIMD) &
       &                 +POTPAR(ISP)%OVERLAP(I,J)*D(J,I,:)
            DO K=1,LMNX
              DO L=1,LMNX
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
                EX=EX+SVAR*SUM(D(K,J,:)*D(L,I,:))
                H(K,J,:)=H(K,J,:)+SVAR*D(L,I,:) 
                H(L,I,:)=H(L,I,:)+SVAR*D(K,J,:) 
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        HAMIL(IND)%MAT=HAMIL(IND)%MAT+H*HFSCALE
        ETOT=ETOT+EX*HFSCALE
        IF(TPR) THEN
          PRINT*,'TOTAL VALENCE CHARGE ON ATOM=..............',IAT,QSPIN(1)
          PRINT*,'TOTAL SPIN[HBAR/2] ON ATOM=................',IAT,QSPIN(2:NDIMD)
          PRINT*,'EXACT VALENCE EXCHANGE ENERGY FOR ATOM=....',IAT,EX
        END IF
        EAT=EAT+EX
        EX=0.D0
!
!       ========================================================================
!       == ADD CORE-VALENCE EXCHANGE                                          ==
!       ========================================================================
        IF(HYBRIDSETTING(ISP)%TCV) THEN
          CALL SIMPLELMTO_CVX_ACTONPHI(IAT,LMNXPHI,NDIMD,DON,EX,HON)
          HONSITE(IAT)%MAT=HONSITE(IAT)%MAT+HON*HFSCALE
          ETOT=ETOT+EX*HFSCALE
          EAT=EAT+EX
          IF(TPR)PRINT*,'CORE-VALENCE EXCHANGE ENERGY FOR ATOM=.....',IAT,EX
          EX=0.D0
          HON=0.D0
        END IF
!
!       ========================================================================
!       == DOUBLE COUNTING CORRECTION (EXCHANGE ONLY)                         ==
!       == THIS IS THE TIME CONSUMING PART                                    ==
!       ========================================================================
        CALL TIMING$CLOCKON('ENERGYTEST:DC')      
!
EX=0.D0
H=0.D0
HON=0.D0
IF(.TRUE.) THEN
        CALL DFT$SETL4('XCONLY',.TRUE.)
        CALL SIMPLELMTO_DC_NEW(ISP,NDIMD,LMNX,D,LMNXPHI,DON,EX,H,HON)
        CALL DFT$SETL4('XCONLY',.FALSE.)
!HONSITE(IAT)%MAT=0.D0
        HONSITE(IAT)%MAT=HONSITE(IAT)%MAT-HON*HFSCALE
        HAMIL(IND)%MAT  =HAMIL(IND)%MAT-H*HFSCALE
        ETOT            =ETOT-EX*HFSCALE  !HFWEIGHT IS MULTIPLIED ON LATER
END IF
!
        IF(TPR) THEN
          IF(.NOT.TPARALLEL) THEN
            PRINT*,'DOUBLE COUNTING CORRECFTION ENERGY FOR ATOM=',IAT,-EX
            PRINT*,'EXACT EXCHANGE ENERGY FOR ATOM........... =',IAT,EAT
            PRINT*,'EXCHANGE-XORRECTION FOR ATOM............. =',IAT,EAT-EX
          END IF
        END IF
!
        CALL TIMING$CLOCKOFF('ENERGYTEST:DC')      
        DEALLOCATE(U)
        DEALLOCATE(H)
        DEALLOCATE(HON)
        DEALLOCATE(DON)
        DEALLOCATE(D)
        DEALLOCATE(AECORE)
        DEALLOCATE(LOX)
      ENDDO !END OF LOOP OVER ATOMS
!
!     ==========================================================================
!     == OFFSITE EXCHANGE CONTRIBUTION                                        ==
!     ==========================================================================
      IF(TOFFSITE) THEN
!       ==  USES NUMERICAL MATRIX ELEMENTS
        CALL SIMPLELMTO_OFFSITEXEVAL(TPARALLEL,NAT,NND,RBAS,R0,DENMAT &
     &                                             ,EX,STRESS,FORCE,HAMIL) 
        PRINT*,'+-+-+-+  OFFSITE EX=',EX
        ETOT=ETOT+EX
      END IF
!
!     ==========================================================================
!     == WRAP UP PARALLELIZATION                                              ==
!     ==========================================================================
      IF(TPARALLEL) THEN
        CALL MPE$COMBINE('MONOMER','+',ETOT)
        CALL MPE$COMBINE('MONOMER','+',STRESS)
        CALL MPE$COMBINE('MONOMER','+',FORCE)
!       __ COMBINE DATH THE HAMILTONIAN IN TERMS OF PARTIAL WAVES_______________
        DO NN=1,NND
          CALL MPE$COMBINE('MONOMER','+',HAMIL(NN)%MAT)
        ENDDO
        DO IAT=1,NAT
          CALL MPE$COMBINE('MONOMER','+',HONSITE(IAT)%MAT)
        ENDDO
      END IF
!
!     ==========================================================================
!     == MAKE HAMILTONIAN HERMITIAN                                           ==
!     ==========================================================================
!      CALL SIMPLELMTO$SYMMETRIZEHAMIL()
!
!     ==========================================================================
!     == RESCALE WITH HFWEIGHT                                                ==
!     ==========================================================================
      ETOT  =ETOT  *HFWEIGHT
      STRESS=STRESS*HFWEIGHT
      FORCE =FORCE *HFWEIGHT
      DO NN=1,NND
        HAMIL(NN)%MAT=HAMIL(NN)%MAT*HFWEIGHT
      ENDDO
      DO IAT=1,NAT
        HONSITE(IAT)%MAT=HONSITE(IAT)%MAT*HFWEIGHT
      ENDDO
!
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_INDEXLOCAL2(IAT,NND,DENMAT,IND)
!     **************************************************************************
!     ** IDENTIFY THE ONSIDE ELEMENT FOR ATOM IAT ON THE DENSITY MATRIX       **
!     **************************************************************************
      USE RSPACEOP_MODULE, ONLY: RSPACEMAT_TYPE

      IMPLICIT NONE
      INTEGER(4)             ,INTENT(IN) :: IAT
      INTEGER(4)             ,INTENT(IN) :: NND
      TYPE(RSPACEMAT_TYPE)   ,INTENT(IN) :: DENMAT(NND)
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
      SUBROUTINE SIMPLELMTO_CVX_ACTONPHI(IAT,LMNX,NDIMD,DENMAT,ETOT,DH)
!     **************************************************************************
!     **  CORE VALENCE EXCHANGE ENERGY ACTING ON PARTIAL WAVES                **
!     *****************************PETER BLOECHL, GOSLAR 2011-2019**************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT          ! ATOM INDEX
      INTEGER(4),INTENT(IN)  :: LMNX
      INTEGER(4),INTENT(IN)  :: NDIMD
      REAL(8)   ,INTENT(OUT) :: ETOT              ! CORE-VALENCE ENERGY
      REAL(8)   ,INTENT(IN)  :: DENMAT(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT) :: DH(LMNX,LMNX,NDIMD)
      INTEGER(4)             :: ISP
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      REAL(8)   ,ALLOCATABLE :: CVXMAT(:,:)
      INTEGER(4)             :: LN1,LN2,L1,L2,LMN1,LMN2,IM
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_CVX_ACTONPHI')
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
      IF(LMNX.NE.SUM(2*LOX+1)) THEN
        CALL ERROR$MSG('INCONSISTENT INPUT VALUE LMNX')
        CALL ERROR$STOP('SIMPLELMTO_CVX_ACTONPHI')
      END IF
!
!     ==========================================================================
!     == SUM UP TOTAL ENERGY                                                  ==
!     ==========================================================================
      DH=0.D0
      ETOT=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        LMN2=0
        DO LN2=1,LNX
          L2=LOX(LN2)
          IF(L1.EQ.L2) THEN
            DO IM=1,2*L1+1
              ETOT=ETOT+CVXMAT(LN1,LN2)*DENMAT(LMN1+IM,LMN2+IM,1)
              DH(LMN1+IM,LMN2+IM,1)=CVXMAT(LN1,LN2)
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
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_DC_NEW(ISP,NDIMD,LMNX_CHI,D_CHI,LMNX_PHI,D_PHI &
     &                            ,ETOT,H_CHI,H_PHI)
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
      USE SIMPLELMTO_MODULE, ONLY : POTPAR 
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: ISP
      INTEGER(4)  ,INTENT(IN) :: NDIMD
      INTEGER(4)  ,INTENT(IN) :: LMNX_CHI     ! #(LOCAL ORBITALS)
      INTEGER(4)  ,INTENT(IN) :: LMNX_PHI     ! #(PARTIAL WAVES)
      REAL(8)     ,INTENT(IN) :: D_CHI(LMNX_CHI,LMNX_CHI,NDIMD) !DENSITY MATRIX
      REAL(8)     ,INTENT(IN) :: D_PHI(LMNX_PHI,LMNX_PHI,NDIMD) !DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT       ! DOUBLE COUNTINNG ENERGY
      REAL(8)     ,INTENT(OUT):: H_CHI(LMNX_CHI,LMNX_CHI,NDIMD) !HAMILTONIAN
      REAL(8)     ,INTENT(OUT):: H_PHI(LMNX_PHI,LMNX_PHI,NDIMD) !HAMILTONIAN
      LOGICAL(4)  ,PARAMETER  :: TPR=.TRUE.
      REAL(8)     ,PARAMETER  :: DELTA=1.D-8
      REAL(8)     ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)     ,PARAMETER  :: FOURPI=4.D0*PI
      REAL(8)     ,PARAMETER  :: Y0=1.D0/SQRT(4.D0*PI)
      INTEGER(4)              :: GID
      INTEGER(4)              :: NR
      REAL(8)     ,ALLOCATABLE:: AECORE(:)
      INTEGER(4)              :: LNX_CHI    !#(RADIAL FUNCTIONS)
      INTEGER(4)  ,ALLOCATABLE:: LOX_CHI(:) !(LNX_CHI) MAIN ANGULAR MOMENTUM 
      REAL(8)     ,ALLOCATABLE:: CHI(:,:)   !(NR,LNX_CHI) LOCAL ORBITALS
      INTEGER(4)              :: LNX_PHI    !#(RADIAL FUNCTIONS)
      INTEGER(4)  ,ALLOCATABLE:: LOX_PHI(:) !(LNX_CHI) MAIN ANGULAR MOMENTUM 
      REAL(8)     ,ALLOCATABLE:: PHI(:,:)   !(NR,LNX_CHI) LOCAL ORBITALS
      COMPLEX(8)  ,ALLOCATABLE:: DENMAT1(:,:,:)
      COMPLEX(8)  ,ALLOCATABLE:: HAM1(:,:,:)
      REAL(8)     ,ALLOCATABLE:: R(:)       !(NR)
      REAL(8)     ,ALLOCATABLE:: CUT(:)     !(NR)
      REAL(8)     ,ALLOCATABLE:: VCUT(:)    !(NR)
      REAL(8)     ,ALLOCATABLE:: RHO_CHI(:,:,:)
      REAL(8)     ,ALLOCATABLE:: POT_CHI(:,:,:)
      REAL(8)     ,ALLOCATABLE:: RHO_PHI(:,:,:)
      REAL(8)     ,ALLOCATABLE:: POT_PHI(:,:,:)
      REAL(8)     ,ALLOCATABLE:: AUX(:)
      REAL(8)                 :: SVAR
      REAL(8)                 :: AEZ
      REAL(8)                 :: RCOV
      INTEGER(4)              :: LMRX
!      INTEGER(4)              :: LMRX,L
!      INTEGER(4)              :: IDIM,LM,LMN,IR
      INTEGER(4)              :: IDIM,IR
      REAL(8),ALLOCATABLE,SAVE :: CUTSAVE(:)
!     **************************************************************************
      ETOT=0.D0
      H_CHI=0.D0
      H_PHI=0.D0
!
      GID=POTPAR(ISP)%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
!
      LNX_CHI=POTPAR(ISP)%LNXH
      ALLOCATE(LOX_CHI(LNX_CHI))
      LOX_CHI=POTPAR(ISP)%LOXH
      ALLOCATE(CHI(NR,LNX_CHI))
      CHI=POTPAR(ISP)%AECHI
     
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETR8('AEZ',AEZ)
      CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV)
      CALL SETUP$GETI4('LMRX',LMRX)
      ALLOCATE(AECORE(NR))
      CALL SETUP$GETR8A('AECORE',NR,AECORE)
!
      CALL SETUP$GETI4('LNX',LNX_PHI)
      ALLOCATE(LOX_PHI(LNX_PHI))
      CALL SETUP$GETI4A('LOX',LNX_PHI,LOX_PHI)
      ALLOCATE(PHI(NR,LNX_PHI))
      CALL SETUP$GETR8A('AEPHI',NR*LNX_PHI,PHI)
      CALL SETUP$UNSELECT()
!
!     ==========================================================================
!     ==  DENSITY OF LOCAL ORBITALS (INCLUDING CORE)                          ==
!     ==========================================================================
      ALLOCATE(RHO_CHI(NR,LMRX,NDIMD))
      ALLOCATE(RHO_PHI(NR,LMRX,NDIMD))
      ALLOCATE(CUT(NR))
!
!     == TAKING THE REAL PART IS NOT APPROXIMATION: ONLY THE REAL PART OF THE 
!     == DENSITY MATRIX IN (TXYZ) REPRESENTATION CONTRIBUTES TO THE DENSITY. 
      ALLOCATE(DENMAT1(LMNX_CHI,LMNX_CHI,NDIMD))
      DENMAT1=CMPLX(D_CHI,KIND=8)
      DO IDIM=1,NDIMD
        CALL AUGMENTATION_RHO(NR,LNX_CHI,LOX_CHI,CHI &
     &                   ,LMNX_CHI,DENMAT1(:,:,IDIM),LMRX,RHO_CHI(:,:,IDIM))
      ENDDO
      DEALLOCATE(DENMAT1)
      RHO_CHI(:,1,1)=RHO_CHI(:,1,1)+AECORE(:)
!
!     ==========================================================================
!     == ONE-CENTER EXPANSION OF THE TOTAL DENSITY                            ==
!     ==========================================================================
      ALLOCATE(DENMAT1(LMNX_PHI,LMNX_PHI,NDIMD))
      DENMAT1=CMPLX(D_PHI,KIND=8)
      DO IDIM=1,NDIMD
        CALL AUGMENTATION_RHO(NR,LNX_PHI,LOX_PHI,PHI &
     &                   ,LMNX_PHI,DENMAT1(:,:,IDIM),LMRX,RHO_PHI(:,:,IDIM))
      ENDDO
      DEALLOCATE(DENMAT1)
      RHO_PHI(:,1,1)=RHO_PHI(:,1,1)+AECORE(:)
!
!     ==========================================================================
!     ==  CUTOFF FUNCTION FOR EXCHANGE-CORRELATION INTEGRAL                   ==
!     ==  CUT IS CLOSE TO UNITY IN THE CENTER AND IS ZERO BEYOND THE ATOM     ==
!     ==========================================================================
!!$      CUT(:)=(RHO_CHI(:,1,1)/(RHO_PHI(:,1,1)+DELTA))**2 
!!$!     == FUDGE FACTOR: CUT OFF THE CUTOFF FUNCTION TO AVOID AN INCREASE 
!!$!     == AT LARGE DISTANCES. 6 ABOHR IS ABOUT A BOND DISTANCE
!!$      DO IR=1,NR
!!$        IF(R(IR).GT.6.D0) THEN
!!$          CUT(IR:)=0.D0
!!$          EXIT
!!$        END IF
!!$      ENDDO
!!$IF(ALLOCATED(CUTSAVE)) THEN
!!$  CUT=CUTSAVE
!!$ELSE
!!$  ALLOCATE(CUTSAVE(NR))
!!$  CUTSAVE=CUT
!!$END IF
!CUT(:)=EXP(-R**2)
      CUT=1.D0
      DO IR=1,NR
        IF(R(IR).GT.RCOV) THEN
          CUT(IR:)=0.D0
          EXIT
        END IF
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE ENERGY AND POTENTIAL                                      ==
!     ==========================================================================
      ALLOCATE(POT_PHI(NR,LMRX,NDIMD))  ! POTENTIAL FOR PARTIAL-WAVE DENSITY 
      ALLOCATE(POT_CHI(NR,LMRX,NDIMD))  ! POTENTIAL FOR LOCAL ORBITALS
      ALLOCATE(VCUT(NR))
      ALLOCATE(AUX(NR))
!
!     == EXCHANGE CORRELATION OF THE TOTAL DENSITY INCLUDING CORE ==============
!     == CALL TEST1_RADXC_WITHCUT(GID,NR,LMRX,NDIMD,RHO_ALL)
      CALL SIMPLELMTO_RADXC_WITHCUT(GID,NR,LMRX,NDIMD,RHO_PHI,CUT &
     &                             ,ETOT,POT_PHI,VCUT)
!
!     == SUBTRACT FROZEN-CORE ==================================================
      CALL SIMPLELMTO_RADXC_WITHCUT(GID,NR,1,1,AECORE,CUT,SVAR,POT_CHI,AUX)
      ETOT=ETOT-SVAR
      VCUT(:)=VCUT(:)-AUX(:)
      POT_CHI=0.D0
      DO IR=1,NR
        IF(R(IR).GT.6.D0) THEN
          VCUT(IR:)=0.D0
          EXIT
        END IF
      ENDDO
VCUT=0.D0
!VCUT=VCUT*Y0
!      
!     == POTENTIAL FOR PARTIAL-WAVE DENSITY N_T ================================
      POT_PHI(:,1,1)=POT_PHI(:,1,1)-VCUT(:)*2.D0*CUT(:)/(RHO_PHI(:,1,1)+DELTA)
!
!     == POTENTIAL FOR THE CORRELATED DENSITY ==================================
      POT_CHI(:,:,:)=0.D0
      POT_CHI(:,1,1)=VCUT(:)*2.D0*CUT(:)/(RHO_CHI(:,1,1)+DELTA)
!
!     ==========================================================================
!     ==  EXTRACT HAMILTON CONTRIBUTIONS                                      ==
!     ==========================================================================
      ALLOCATE(HAM1(LMNX_CHI,LMNX_CHI,NDIMD))
      HAM1=(0.D0,0.D0)
      CALL SIMPLELMTO_EXPECT(GID,NR,NDIMD,LNX_CHI,LOX_CHI,LMNX_CHI &
     &                     ,LMRX,POT_CHI,CHI,HAM1)
      H_CHI=REAL(HAM1)
      DEALLOCATE(HAM1)
!
      ALLOCATE(HAM1(LMNX_PHI,LMNX_PHI,NDIMD))
      HAM1=(0.D0,0.D0)

      CALL SIMPLELMTO_EXPECT(GID,NR,NDIMD,LNX_PHI,LOX_PHI,LMNX_PHI &
     &                     ,LMRX,POT_PHI,PHI,HAM1)
      H_PHI=REAL(HAM1)
      DEALLOCATE(HAM1)
!!$PRINT*,'H_PHI',H_PHI
!!$PRINT*,'H_CHI',H_CHI
!
!     ==========================================================================
!     ==  PRINT DIAGNOSTIC INFORMATION                                        ==
!     ==========================================================================
      IF(TPR) THEN
        PRINT*,'LMRX',LMRX
        CALL LMTO_WRITEPHI('POT_CHI.DAT',GID,NR,LMRX,POT_CHI)
        CALL LMTO_WRITEPHI('POT_PHI.DAT',GID,NR,LMRX,POT_PHI)
        CALL LMTO_WRITEPHI('RHO_CHI.DAT',GID,NR,LMRX,RHO_CHI)
        CALL LMTO_WRITEPHI('RHO_PHI.DAT',GID,NR,LMRX,RHO_PHI)
        CALL LMTO_WRITEPHI('CUT.DAT',GID,NR,1,CUT)
        CALL LMTO_WRITEPHI('VCUT.DAT',GID,NR,1,VCUT)
        CALL RADIAL$INTEGRAL(GID,NR &
    &                      ,FOURPI*R**2*(RHO_PHI(:,1,1)-AECORE(:))*Y0,SVAR)
        PRINT*,'RADXC VALENCE CHARGE (PARTIAL WAVE EXPANSION)=',SVAR
        CALL RADIAL$INTEGRAL(GID,NR &
    &                       ,FOURPI*R**2*(RHO_CHI(:,1,1)-AECORE(:))*Y0,SVAR)
        PRINT*,'RADXC VALENCE CHARGE (LOCAL ORBITALS)=        ',SVAR
        CALL RADIAL$INTEGRAL(GID,NR,FOURPI*R**2*AECORE*Y0,SVAR)
        PRINT*,'RADXC CORE CHARGE                    =        ',SVAR
!!$        CALL ERROR$MSG('REGULAR STOP AFTER PRINTING DIAGNOSTIC INFORMATION')
!!$        CALL ERROR$STOP('SIMPLELMTO_DC')
      END IF
!
      DEALLOCATE(RHO_PHI)
      DEALLOCATE(RHO_CHI)
      DEALLOCATE(POT_PHI)
      DEALLOCATE(POT_CHI)
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TEST1_RADXC_WITHCUT(GID,NR,LMRX,NDIMD,RHO)
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
      LM=1
      L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
      DRHO(:,:,:)=0.D0
      DRHO(:,LM,1)=1.D-2*R**L*EXP(-R**2)
      DCUT(:)=1.D-1*EXP(-2.D0*R**2)
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
      CALL ERROR$MSG('REGULAR STOP IN TESTING ROUTINE')
      CALL ERROR$STOP('TEST1_RADXC_WITHCUT')
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_RADXC_WITHCUT(GID,NR,LMRX,NDIMD,RHOIN,CUT &
     &                                   ,EXC,VXC,VCUT)
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
        VCUT(IR)    =FXC(IR)/FOURPI
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
      SUBROUTINE SIMPLELMTO_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX &
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
                CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
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
      SUBROUTINE SIMPLELMTO_WRITEPHI(FILE,GID,NR,NPHI,PHI)
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
        WRITE(100,FMT='(F15.10,2X,20(F25.15,2X))')R(IR) &
     &               ,MIN(MAX(PHI(IR,:),-1.D8),1.D8)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!*******************************************************************************
!*******************************************************************************
!****                         OFF-SITE MATRIX ELEMENTS                     *****
!**** DATE ARE PREPARED IN SIMPLELMTO_OFFXINT()                            *****
!**** ENERGIES ARE EVALUEATED IN SIMPLELMTO_OFFSITEXEVAL                   *****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFXINT()
!     **************************************************************************
!     ** COMPUTES OFFSITE MATRIX ELEMENTS OFFSITEX
!     ** WATCH PARALLELIZATION!!!
!     **************************************************************************
      USE SIMPLELMTO_MODULE,ONLY : POTPAR &
     &                      ,OFFSITEX &
     &                      ,HYBRIDSETTING &
     &                      ,NSP
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      REAL(8)   ,PARAMETER :: TOLERANCE=1.D-3
      INTEGER(4),PARAMETER :: NDIS=5  !#(DISTANCE GRID POINTS)
      REAL(8)   ,PARAMETER :: DSMN=0.5D0,DSMX=5.D0 !PARMS FOR DISTANCE GRID
      INTEGER(4),PARAMETER :: NF=3    !#(FIT FUNCTIONS)
!OLD  REAL(8)   ,PARAMETER :: DCAYMN=0.8D0,DCAYMX=2.D0 !PARMS FOR FIT FUNCTIONS
      REAL(8)   ,PARAMETER :: DCAYMN=0.8D0,DCAYMX=2.D0 !PARMS FOR FIT FUNCTIONS
!      REAL(8)   ,PARAMETER :: DCAYMN=0.5D0,DCAYMX=3.D0 !PARMS FOR FIT FUNCTIONS
      INTEGER(4)           :: GID1,GID2
      INTEGER(4)           :: NR1,NR2
      INTEGER(4)           :: LMNX1,LMNX2
      INTEGER(4)           :: LMN1,LMN2,LMN3,LMN4
      INTEGER(4)           :: ISP,ISP1,ISP2,I
      REAL(8)              :: SVAR,SVAR1,AEZ
      REAL(8)              :: RCOV(NSP)
      REAL(8)              :: DIS 
!     **************************************************************************
                                  CALL TRACE$PUSH('SIMPLELMTO_OFFXINT')
PRINT*,'STARTING INITIALIZATION OF SIMPLELMTO_OFFSITE'

      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV(ISP))
        CALL SETUP$UNSELECT()
      ENDDO

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
          OFFSITEX(ISP1,ISP2)%NDIS=NDIS
          ALLOCATE(OFFSITEX(ISP1,ISP2)%DIS(NDIS))
! RAUG IS A POOR CHOICE FOR ESTIMATING THE GRID BOUNDARIES
!          SVAR=POTPAR(ISP1)%RAUG+POTPAR(ISP2)%RAUG
! CHANGED TO RCOV 22.JAN.2019
          SVAR=RCOV(ISP1)+RCOV(ISP2)
          DO I=1,NDIS
            DIS=SVAR*(DSMN+(DSMX-DSMN)*REAL(I-1,KIND=8)/REAL(NDIS-1,KIND=8))
            OFFSITEX(ISP1,ISP2)%DIS(I)=DIS
          ENDDO
!
!         == DEFINE DECAY CONSTANTS FOR INTERPOLATING FUNCTIONS ================
          IF(NDIS.LT.NF) THEN
            CALL ERROR$MSG('NUMBER OF INTERPOLATING FUNCTIONS MUST')
            CALL ERROR$MSG('BE EQUAL OR GREATER THAN NUMBER OF GRID POINTS')
            CALL ERROR$STOP('SIMPLELMTO_OFFSITEXINT')
          END IF
          OFFSITEX(ISP1,ISP2)%NF=NF
          ALLOCATE(OFFSITEX(ISP1,ISP2)%LAMBDA(NF))
          SVAR=(RCOV(ISP1)+RCOV(ISP2))
          DO I=1,NF
            SVAR1=1.D0/(SVAR*(DCAYMN+REAL(I-1,KIND=8)*(DCAYMX-DCAYMN)))
            OFFSITEX(ISP1,ISP2)%LAMBDA(I)=SVAR1
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE OVERLAP MATRIX ELEMENTS                                    ==
!     ==========================================================================
      DO ISP1=1,NSP
        GID1=POTPAR(ISP1)%GID
        CALL RADIAL$GETI4(GID1,'NR',NR1)
        DO ISP2=1,NSP
          GID2=POTPAR(ISP2)%GID
          CALL RADIAL$GETI4(GID2,'NR',NR2)
PRINT*,'DOING OVERLAP....',ISP1,ISP2
          CALL SIMPLELMTO_OFFSITEOVERLAPSETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE EXCHANGE INTEGRALS WITH TWO ORBITALS ON EITHER SIDE        ==
!     ==========================================================================
      DO ISP1=1,NSP
        IF(.NOT.HYBRIDSETTING(ISP1)%TNDDO) CYCLE
        GID1=POTPAR(ISP1)%GID
        CALL RADIAL$GETI4(GID1,'NR',NR1)
        DO ISP2=1,NSP
          IF(.NOT.HYBRIDSETTING(ISP2)%TNDDO) CYCLE
          GID2=POTPAR(ISP2)%GID
          CALL RADIAL$GETI4(GID2,'NR',NR2)
PRINT*,'DOING X22 ....',ISP1,ISP2
          CALL SIMPLELMTO_OFFSITEX22SETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE EXCHANGE INTEGRALS WITH THREE ORBITALS ON THE FIRST SITE   ==
!     ==                               AND ONE ORBITAL ON THE SECOND SITE     ==
!     ==========================================================================
      DO ISP1=1,NSP
        IF(.NOT.HYBRIDSETTING(ISP1)%T31) CYCLE
        GID1=POTPAR(ISP1)%GID
        CALL RADIAL$GETI4(GID1,'NR',NR1)
        DO ISP2=1,NSP
          IF(.NOT.HYBRIDSETTING(ISP2)%T31) CYCLE
          GID2=POTPAR(ISP2)%GID
          CALL RADIAL$GETI4(GID2,'NR',NR2)
PRINT*,'DOING X31 ....',ISP1,ISP2
          CALL SIMPLELMTO_OFFSITEX31SETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE EXCHANGE INTEGRALS DABAB  ==
!     ==========================================================================
!!$PRINT*,'DOING XABAB....'
!!$CALL TRACE$PASS('DOING XABAB')
!!$      CALL SIMPLELMTO_TAILEDGAUSSOFFSITEU()   ! ROUTINE PARALLELIZES OVER 'MONOMER'
!!$CALL TRACE$PASS('XABAB DONE')
!
!     ==========================================================================
!     == CONVERT INTEGRALS INTO COEFFICIENTS OF INTERPOLATING FUNCTION        ==
!     ==========================================================================
PRINT*,'CONVERTING....'
CALL TRACE$PASS('CONVERTING')
      CALL SIMPLELMTO_OFFSITEXCONVERT()
CALL TRACE$PASS('CONVERSION DONE')
PRINT*,'INITIALIZATION OF OFFSITE DONE...'
                                  CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEOVERLAPSETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
!     **************************************************************************
!     ** CALCULATES THE OVERLAP OF TAILED ORBITALS ON A RADIAL INTERPOLATION  **
!     ** GRID                                                                 **
!     ** ROUTINE IS PARALLELIZED OVER 'MONOMER'                               **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR &
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
      GID1=POTPAR(ISP1)%GID
      GID2=POTPAR(ISP2)%GID
      LNX1=POTPAR(ISP1)%LNXH
      LNX2=POTPAR(ISP2)%LNXH
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
     &              ,(RGRID1(:)*POTPAR(ISP1)%AECHI(:,LN1))**2,TOLFAC1(LN1))
      ENDDO
      CALL RADIAL$R(GID2,NR2,RGRID2)
      DO LN2=1,LNX2
        CALL RADIAL$INTEGRAL(GID2,NR2 &
    &               ,(RGRID2(:)*POTPAR(ISP2)%AECHI(:,LN2))**2,TOLFAC2(LN2))
      ENDDO
      TOLFAC1=SQRT(TOLFAC1)
      TOLFAC2=SQRT(TOLFAC2)
!
!     ==========================================================================
!     == CALCULATE NUMBER OF INTEGRALS                                        ==
!     ==========================================================================
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=1,LNX2
          L2=POTPAR(ISP2)%LOXH(LN2)
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
        L1=POTPAR(ISP1)%LOXH(LN1)
        PHI1(:)=POTPAR(ISP1)%AECHI(:,LN1) 
        DO LN2=1,LNX2
          L2=POTPAR(ISP2)%LOXH(LN2)
          PHI2(:)=POTPAR(ISP2)%AECHI(:,LN2)
          TOL=TOLERANCE*TOLFAC1(LN1)*TOLFAC2(LN2)
          TOL=MAX(TOLMIN,TOL)
          DO MABS=0,MIN(L1,L2)
            IND=IND+1
            DO IDIS=1,NDIS
              COUNT=COUNT+1
              IF(MOD(COUNT-1,NTASKS).NE.THISTASK-1) CYCLE
              DIS=OFFSITEX(ISP1,ISP2)%DIS(IDIS)
              CALL SIMPLELMTO_TWOCENTER(L1,MABS,GID1,NR1,PHI1 &
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
          L1=POTPAR(ISP1)%LOXH(LN1)
          DO LN2=1,LNX2
            L2=POTPAR(ISP2)%LOXH(LN2)
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
      SUBROUTINE SIMPLELMTO_OFFSITEX22SETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR &
     &                             ,OFFSITEX &
     &                             ,SCREENL
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      INTEGER(4),INTENT(IN) :: NR1
      INTEGER(4),INTENT(IN) :: NR2
      REAL(8)   ,INTENT(IN) :: TOLERANCE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
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
      GID1=POTPAR(ISP1)%GID
      GID2=POTPAR(ISP2)%GID
      LNX1=POTPAR(ISP1)%LNXH
      LNX2=POTPAR(ISP2)%LNXH
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
     &                   ,(RGRID1(:)*POTPAR(ISP1)%AECHI(:,LN1))**2,TOLFAC1(LN1))
      ENDDO
      CALL RADIAL$R(GID2,NR2,RGRID2)
      DO LN2=1,LNX2
        CALL RADIAL$INTEGRAL(GID2,NR2 &
    &                    ,(RGRID2(:)*POTPAR(ISP2)%AECHI(:,LN2))**2,TOLFAC2(LN2))
      ENDDO
      TOLFAC1=SQRT(TOLFAC1)
      TOLFAC2=SQRT(TOLFAC2)
!
!     ==========================================================================
!     == CALCULATE NUMBER OF INTEGRALS                                        ==
!     ==========================================================================
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%LOXH(LN2)
          DO LN3=1,LNX2
            L3=POTPAR(ISP2)%LOXH(LN3)
            DO LN4=LN3,LNX2
              L4=POTPAR(ISP2)%LOXH(LN4)
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
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%LOXH(LN2)
          RHO12(:)=POTPAR(ISP1)%AECHI(:,LN1) &
       &          *POTPAR(ISP1)%AECHI(:,LN2)
          DO LN3=1,LNX2
            L3=POTPAR(ISP2)%LOXH(LN3)
            DO LN4=LN3,LNX2
              L4=POTPAR(ISP2)%LOXH(LN4)
              RHO34(:)=POTPAR(ISP2)%AECHI(:,LN3) &
       &              *POTPAR(ISP2)%AECHI(:,LN4)
              TOL=TOLERANCE*TOLFAC1(LN1)*TOLFAC1(LN2) &
       &                   *TOLFAC2(LN3)*TOLFAC2(LN4)
              TOL=MAX(TOLMIN,TOL)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                IF(SCREENL.LT.0.D0) THEN
                  CALL RADIAL$POISSON(GID1,NR1,LR1,RHO12,POT12)
                ELSE
                  CALL RADIAL$YUKAWA(GID1,NR1,LR1,1.D0/SCREENL,RHO12,POT12)
                END IF
!!$!               ==SUBTRACT OUT LONG RANGE PART INCLUDING MONOPOLE AND DIPOLE TERMS
!!$                IF(LR1.EQ.0) THEN
!!$                  SVAR=POTPAR(ISP1)%QLN(1,LN1,LN2)*SQ4PI
!!$                  POT12(:)=POT12(:)-SVAR/RGRID1(:)
!!$                ELSE IF(LR1.EQ.1) THEN
!!$                  SVAR=POTPAR(ISP1)%QLN(2,LN1,LN2)*SQ4PITHIRD
!!$                  POT12(:)=POT12(:)-SVAR/RGRID1(:)**2
!!$                END IF
!
                DO LR2=ABS(L3-L4),MIN(L3+L4,LRX2),2
!!$                  CALL RADIAL$POISSON(GID2,NR2,LR2,RHO34,POT34)
!!$!                 ==SUBTRACT OUT LONG RANGE PART INCLUDING MONOPOLE AND DIPOLE TERMS
!!$                  IF(LR2.EQ.0) THEN
!!$                    SVAR=POTPAR(ISP2)%QLN(1,LN3,LN4)*SQ4PI
!!$                    POT34(:)=POT34(:)-SVAR/RGRID2(:)
!!$                  ELSE IF(LR2.EQ.1) THEN
!!$                    SVAR=POTPAR(ISP1)%QLN(2,LN3,LN4)*SQ4PITHIRD
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
                      CALL SIMPLELMTO_TWOCENTER(LR1,MABS,GID1,NR1,POT12 &
                                         ,LR2,MABS,GID2,NR2,RHO34 &
       &                                 ,DIS,TOL,INTEGRAL)
!                     == SUBTRACT OUT LONG RANGE PART TO ALLOW INTERPOLATION ===
!                     == WILL BE ADDED AGAIN ===================================
                      IF(LR1.LE.1.AND.LR2.LE.1) THEN
                        SVAR=POTPAR(ISP1)%QLN(LR1+1,LN1,LN2) &
     &                      *POTPAR(ISP2)%QLN(LR2+1,LN3,LN4)
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
!
!     ==========================================================================
!     == DIAGNOSTIC PRINTOUT                                                  ==
!     ==========================================================================
      IF(TPR) THEN
        PRINT*,'X22  FOR ',ISP1,ISP2
        DO IDIS=1,NDIS
          WRITE(*,FMT='(10F10.5)')OFFSITEX(ISP1,ISP2)%DIS(IDIS) &
     &                           ,OFFSITEX(ISP1,ISP2)%X22(IDIS,1)/(4.D0*PI)
        ENDDO
        CALL ERROR$MSG('REGULAR STOP AFTER PRINTING DIAGNOSTIC INFORMATION')
        CALL ERROR$STOP('SIMPLELMTO_OFFSITEX22SETUP')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEX31SETUP(ISP1,ISP2,NR1,NR2,TOLERANCE)
!     **************************************************************************
!     **  CALC. U-TENSOR FOR THREE ORBITALS ON ONE ATOM AND ON ON THE OTHER   **
!     **  ROUTINE IS PARALLELIZED OVER 'MONOMER'                              **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR &
     &                             ,OFFSITEX &
     &                             ,SCREENL
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
      GID1=POTPAR(ISP1)%GID
      GID2=POTPAR(ISP2)%GID
      LNX1=POTPAR(ISP1)%LNXH
      LNX2=POTPAR(ISP2)%LNXH
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
     &              ,(RGRID1(:)*POTPAR(ISP1)%AECHI(:,LN1))**2,TOLFAC1(LN1))
      ENDDO
      CALL RADIAL$R(GID2,NR2,RGRID2)
      DO LN2=1,LNX2
        CALL RADIAL$INTEGRAL(GID2,NR2 &
    &               ,(RGRID2(:)*POTPAR(ISP2)%AECHI(:,LN2))**2,TOLFAC2(LN2))
      ENDDO
      TOLFAC1=SQRT(TOLFAC1)
      TOLFAC2=SQRT(TOLFAC2)
!
!     ==========================================================================
!     == CALCULATE NUMBER OF INTEGRALS                                        ==
!     ==========================================================================
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%LOXH(LN2)
          DO LN3=1,LNX1
            L3=POTPAR(ISP1)%LOXH(LN3)
            DO LN4=1,LNX2
              L4=POTPAR(ISP2)%LOXH(LN4)
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
!     == DETERMINE INTEGRALS                                                  ==
!     ==========================================================================
      COUNT=0
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%LOXH(LN2)
          RHO12(:)=POTPAR(ISP1)%AECHI(:,LN1) &
       &          *POTPAR(ISP1)%AECHI(:,LN2)
          DO LN3=1,LNX1
            L3=POTPAR(ISP1)%LOXH(LN3)
            DO LN4=1,LNX2
              L4=POTPAR(ISP2)%LOXH(LN4)
              PHI4=POTPAR(ISP2)%AECHI(:,LN4)
              TOL=TOLERANCE*TOLFAC1(LN1)*TOLFAC1(LN2) &
       &                   *TOLFAC1(LN3)*TOLFAC2(LN4) 
              TOL=MAX(TOLMIN,TOL)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                IF(SCREENL.LT.0.D0) THEN
                  CALL RADIAL$POISSON(GID1,NR1,LR1,RHO12,POT12)
                ELSE
                  CALL RADIAL$YUKAWA(GID1,NR1,LR1,1.D0/SCREENL,RHO12,POT12)
                END IF
                A123(:)=POT12(:)*POTPAR(ISP1)%AECHI(:,LN3)
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
                    IND=IND+1
                    DO IDIS=1,NDIS
                      DIS=OFFSITEX(ISP1,ISP2)%DIS(IDIS)
                      COUNT=COUNT+1
                      IF(MOD(COUNT-1,NTASKS).NE.THISTASK-1) CYCLE
                      CALL SIMPLELMTO_TWOCENTER(LR2,MABS,GID1,NR1,A123 &
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
!     ==========================================================================
!     == DIAGNOSTIC PRINTOUT                                                  ==
!     ==========================================================================
      IF(TPR) THEN
        PRINT*,'X31  FOR ',ISP1,ISP2
        DO IDIS=1,NDIS
          PRINT*,OFFSITEX(ISP1,ISP2)%DIS(IDIS) &
     &          ,OFFSITEX(ISP1,ISP2)%X31(IDIS,:) !/(4.D0*PI)
        ENDDO
        STOP 'FORCED'
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEXCONVERT()
!     **************************************************************************
!     ** CONVERT THE DISTANCE-DEPENDENT MATRIX ELEMENTS OF THE U-TENSOR       **
!     ** INTO EXPANSION COEFFICIENTS FOR THE INTERPOLATING FUNCTION           **
!     ** F_J(X)=SUM_I=1^NDIS X22(I,J)*EXP(-LAMBDA(I)*X)                       **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : OFFSITEX &
     &                             ,POTPAR &
     &                             ,NSP
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
                                   CALL TRACE$PUSH('SIMPLELMTO_OFFSITEXCONVERT')
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
!!$PRINT*,'X31-A:',OFFSITEX(ISP1,ISP2)%X31(:NDIS,:)
!!$            CALL SIMPLELMTO_TESTOFFSITEXCONVERT_A(NDIS &
!!$        &                                     ,OFFSITEX(ISP1,ISP2)%DIS &
!!$        &                                     ,OFFSITEX(ISP1,ISP2)%X31(:,1))
            DO I=1,NIND
              OFFSITEX(ISP1,ISP2)%X31(:NF,I)=MATMUL(AMATIN &
       &                                         ,OFFSITEX(ISP1,ISP2)%X31(:,I))
            ENDDO
            OFFSITEX(ISP1,ISP2)%X31(NF+1:,:)=0.D0
!!$PRINT*,'X31-B:',OFFSITEX(ISP1,ISP2)%X31(:NF,:)
!!$            CALL SIMPLELMTO_TESTOFFSITEXCONVERT_B(NF &
!!$       &                   ,OFFSITEX(ISP1,ISP2)%LAMBDA(:NF) &
!!$       &                   ,OFFSITEX(ISP1,ISP2)%X31(:NF,1))
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
                                                                CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_TESTOFFSITEXCONVERT_A(NDIS,DIS,X)
!     **************************************************************************
!     ** THIS AND THE PARTNER ROUTINE ..._B IS USED TO INSPECT THE FIT OF THE **
!     ** DISTANCE-DEPENDENT U-TENSOR MATRIX ELEMENTS BY A SET OF EXPONENTIALS.**
!     ** THE FIRST ROUTINE ..._A WRITES THE DISTANCE-DEPENDENT DATA TO A FILE **
!     ** TEST_A.DAT, WHILE THE SECOND WRITES THE INTERPOLATING FUNCTION TO    **
!     ** FILE TEST_B.DAT. THE SECOND ROUTINE STOPS EXECUTION.
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NDIS
      REAL(8)   ,INTENT(IN) :: DIS(NDIS)
      REAL(8)   ,INTENT(IN) :: X(NDIS)
      INTEGER(4)            :: I
      INTEGER(4)            :: NFIL
!     **************************************************************************
      CALL FILEHANDLER$SETFILE('TEST',.FALSE.,-'TEST_A.DAT')
      CALL FILEHANDLER$SETSPECIFICATION('TEST','FORM','FORMATTED')
      CALL FILEHANDLER$SETSPECIFICATION('TEST','STATUS','REPLACE')
      CALL FILEHANDLER$UNIT('TEST',NFIL)
PRINT*,'X(D)  ',X
      DO I=1,NDIS
        WRITE(NFIL,FMT='(F10.5,10F15.5)')DIS(I),X(I)
      ENDDO
      CALL FILEHANDLER$CLOSE('TEST')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_TESTOFFSITEXCONVERT_B(NF,LAMBDA,X)
!     **************************************************************************
!     ** THIS AND THE PARTNER ROUTINE ..._A IS USED TO INSPECT THE FIT OF THE **
!     ** DISTANCE-DEPENDENT U-TENSOR MATRIX ELEMENTS BY A SET OF EXPONENTIALS.**
!     ** THE FIRST ROUTINE ..._A WRITES THE DISTANCE-DEPENDENT DATA TO A FILE **
!     ** TEST_A.DAT, WHILE THE SECOND WRITES THE INTERPOLATING FUNCTION TO    **
!     ** FILE TEST_B.DAT. THE SECOND ROUTINE STOPS EXECUTION.
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NF
      REAL(8)   ,INTENT(IN) :: LAMBDA(NF)
      REAL(8)   ,INTENT(IN) :: X(NF)
      REAL(8)   ,PARAMETER  :: DMIN=0.D0
      REAL(8)   ,PARAMETER  :: DMAX=5.D0
      INTEGER(4),PARAMETER  :: ND=100
      INTEGER(4)            :: I
      INTEGER(4)            :: NFIL
      REAL(8)               :: D
      REAL(8)               :: FINT
!     **************************************************************************
      CALL FILEHANDLER$SETFILE('TEST',.FALSE.,-'TEST_B.DAT')
      CALL FILEHANDLER$SETSPECIFICATION('TEST','FORM','FORMATTED')
      CALL FILEHANDLER$SETSPECIFICATION('TEST','STATUS','REPLACE')
      CALL FILEHANDLER$UNIT('TEST',NFIL)
PRINT*,'LAMBDA     ',LAMBDA
PRINT*,'X(LAMBDA)  ',X
      DO I=1,ND
        D=DMIN+(DMAX-DMIN)*REAL(I-1,KIND=8)/REAL(ND-1,KIND=8)
        FINT=SUM(X(:)*EXP(-LAMBDA(:)*D))
        WRITE(NFIL,FMT='(F10.5,10F15.5)')D,FINT
      ENDDO
      CALL FILEHANDLER$CLOSE('TEST')
      STOP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_PRBONDU()
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR &
     &                             ,OFFSITEX
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
        CALL ERROR$STOP('SIMPLELMTO_PRBONDU')
      END IF
      LMNXA=POTPAR(ISP1)%LMNXH
      LMNXB=POTPAR(ISP2)%LMNXH
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
      CALL ERROR$STOP('SIMPLELMTO_PRBONDU')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEOVERLAP(ISP1,ISP2,DIS,LMNX1,LMNX2,O,DO)
!     **************************************************************************
!     ** OVERLAP MATRIX ELEMENTS FOR ORBITALS ON TWO ATOMS                    **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR
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
                                  CALL TRACE$PUSH('SIMPLELMTO_OFFSITEOVERLAP')
!
!     ==========================================================================
!     == PREPARATION                                                          ==
!     ==========================================================================
      LNX1=POTPAR(ISP1)%LNXH
      LNX2=POTPAR(ISP2)%LNXH
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
        L1=POTPAR(ISP1)%LOXH(LN1)
        LMN0A(LN1+1)=LMN0A(LN1)+2*L1+1
      ENDDO                
      LMN0B(1)=0
      DO LN3=1,LNX2-1
        L3=POTPAR(ISP2)%LOXH(LN3)
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
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=1,LNX2
          L2=POTPAR(ISP2)%LOXH(LN2)
          DO MABS=0,MIN(L1,L2)
            IND=IND+1
!
!           == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE... ================
            CALL SIMPLELMTO_OFFSITEXVALUE('OV',ISP1,ISP2,IND,DIS,INTEGRAL,DINTEGRAL)
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
      SUBROUTINE SIMPLELMTO_OFFSITEX22U(ISP1,ISP2,DIS,LMNX1,LMNX2,U,DU)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS                                               **
!     ** TWO ORBITALS ARE ON THE FIRST AND TWO ORBITALS ON THE SECOND ATOM    **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR
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
      LNX1=POTPAR(ISP1)%LNXH
      LNX2=POTPAR(ISP2)%LNXH
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
        L1=POTPAR(ISP1)%LOXH(LN1)
        LMN0A(LN1+1)=LMN0A(LN1)+2*L1+1
      ENDDO                
      LMN0B(1)=0
      DO LN3=1,LNX2-1
        L3=POTPAR(ISP2)%LOXH(LN3)
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
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%LOXH(LN2)
          DO LN3=1,LNX2
            L3=POTPAR(ISP2)%LOXH(LN3)
            DO LN4=LN3,LNX2
              L4=POTPAR(ISP2)%LOXH(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                DO LR2=ABS(L3-L4),MIN(L3+L4,LRX2),2
                  DO MABS=0,MIN(LR1,LR2)
                    IND=IND+1
!
!                   == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE...
                    CALL SIMPLELMTO_OFFSITEXVALUE('22',ISP1,ISP2,IND,DIS &
     &                                          ,INTEGRAL,DINTEGRAL)
                    IF(LR1.LE.1.AND.LR2.LE.1) THEN
                      SVAR=POTPAR(ISP1)%QLN(LR1+1,LN1,LN2) &
     &                    *POTPAR(ISP2)%QLN(LR2+1,LN3,LN4)
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
      SUBROUTINE SIMPLELMTO_OFFSITEX31U(ISP1,ISP2,DIS,LMNX1,LMNX2,U,DU)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS 
!     ** TWO ORBITALS ARE ON THE FIRST AND TWO ORBITALS ON THE SECOND ATOM    **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR
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
      INTEGER(4),PARAMETER  :: VERSION=0
!         VERSION=0  ORIGINAL VERSION
!         VERSION=2  OWN SUBROUTINE VERSION
!     **************************************************************************
!
!     ==========================================================================
!     == PREPARATION                                                          ==
!     ==========================================================================
      LNX1=POTPAR(ISP1)%LNXH
      LNX2=POTPAR(ISP2)%LNXH
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
        L1=POTPAR(ISP1)%LOXH(LN1)
        LMN0A(LN1+1)=LMN0A(LN1)+2*L1+1
      ENDDO                
      LMN0B(1)=0
      DO LN3=1,LNX2-1
        L3=POTPAR(ISP2)%LOXH(LN3)
        LMN0B(LN3+1)=LMN0B(LN3)+2*L3+1
      ENDDO                
!
!     ==========================================================================
!     == OBTAIN U-TENSOR                                                      ==
!     ==========================================================================
!!$WRITE(*,FMT='(80("="),T20,"  OFFX31  ")')
!!$WRITE(*,FMT='("ISP1/2   ",T20,2I5)')ISP1,ISP2
!!$WRITE(*,FMT='("LNX1/2   ",T20,2I5)')LNX1,LNX2
!!$WRITE(*,FMT='("LOX(ISP1)",T20,20I5)')POTPAR(ISP1)%LOXH(:)
!!$WRITE(*,FMT='("LOX(ISP2)",T20,20I5)')POTPAR(ISP2)%LOXH(:)
!!$WRITE(*,FMT='("LMN0A    ",T20,20I5)')LMN0A
!!$WRITE(*,FMT='("LMN0B    ",T20,20I5)')LMN0B
IF(VERSION.EQ.0) THEN
      U(:,:,:,:)=0.D0
      DU(:,:,:,:)=0.D0
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%LOXH(LN2)
          DO LN3=1,LNX1
            L3=POTPAR(ISP1)%LOXH(LN3)
            DO LN4=1,LNX2
              L4=POTPAR(ISP2)%LOXH(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
                    IND=IND+1
!
!                   == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE...
                    CALL SIMPLELMTO_OFFSITEXVALUE('31',ISP1,ISP2,IND,ABS(DIS) &
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
!PRINT*,'++TESTU ',DIS,LMN1,LMN2,LMN3,LMN4,U(LMN1,LMN2,LMN3,LMN4),CG1,CG2,INTEGRAL
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
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%LOXH(LN2)
          DO LN3=1,LNX1
            L3=POTPAR(ISP1)%LOXH(LN3)
            DO LN4=1,LNX2
              L4=POTPAR(ISP2)%LOXH(LN4)
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
                DO LR2=ABS(L3-LR1),L3+LR1,2
                  DO MABS=0,MIN(LR2,L4)
                    IND=IND+1
!
!                   == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE...
                    CALL SIMPLELMTO_OFFSITEXVALUE('31',ISP1,ISP2,IND,ABS(DIS) &
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
      CALL SIMPLELMTO_OFFSITEX31U_TEST1(ISP1,LNX1,POTPAR(ISP1)%LOXH,LMNX1 &
     &                           ,ISP2,LNX2,POTPAR(ISP2)%LOXH,LMNX2 &
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
          DU(LMN2,LMN1,:,:)=DU(LMN1,LMN2,:,:)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEX31U_TEST1(ISPA,LNXA,LOXA,LMNXA &
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
      LMXA=(LXA+1)**2
      LMXB=(LXB+1)**2
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
                    CALL SIMPLELMTO_OFFSITEXVALUE('31',ISPA,ISPB,IND,ABS(DIS) &
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
PRINT*,'PM,M3,M2,M1,LMR1I,LMR1F',PM,M3,M2,M1,LMR1I,LMR1F,SVAR
PRINT*,'LM1,LM2,LM3,LMR2,LMR1I,LMR1F,SVAR',LM1,LM2,LM3,LMR2,LMR1I,LMR1F,SVAR,CGMAT(1,1,1)
                            IF(SVAR.EQ.0.D0) CYCLE
                            II=II+1
                            IF(II.GT.IIX) THEN
                              CALL ERROR$STOP('INDEX II EXCEEDS ARRAY SIZE')
                              CALL ERROR$STOP('SIMPLELMTO_OFFSITEX31U_TEST1')
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
PRINT*,'START31LOOP ',LNXA,LNXB,LOXA,LOXB
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
PRINT*,'L1,L2,L3,L4,LR1,LR2,MABS',L1,L2,L3,L4,LR1,LR2,MABS
PRINT*,'II1,II2',II1,II2
                    DO II=II1,II2
                      LMN1=LMN0A(LN1)+LIST(II)%M1
                      LMN2=LMN0A(LN2)+LIST(II)%M2
                      LMN3=LMN0A(LN3)+LIST(II)%M3
                      LMN4=LMN0B(LN4)+LIST(II)%M4
                      U(LMN1,LMN2,LMN3,LMN4)=U(LMN1,LMN2,LMN3,LMN4) &
     &                                      +LIST(II)%VAL*INTEGRAL(IND)
                      DU(LMN1,LMN2,LMN3,LMN4)=DU(LMN1,LMN2,LMN3,LMN4) &
     &                                      +LIST(II)%VAL*DINTEGRAL(IND)
PRINT*,'++',LIST(II)%VAL,INTEGRAL(IND),U(LMN1,LMN2,LMN3,LMN4),DU(LMN1,LMN2,LMN3,LMN4)
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
      SUBROUTINE SIMPLELMTO_OFFSITEXBONDU(ISP1,ISP2,DIS,LMNX1,LMNX2,U,DU)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS 
!     ** TWO ORBITALS ARE ON THE FIRST AND TWO ORBITALS ON THE SECOND ATOM    **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR
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
              CALL SIMPLELMTO_OFFSITEXVALUE('BONDU',ISP1,ISP2,IND,ABS(DIS) &
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
      SUBROUTINE SIMPLELMTO_OFFSITEXVALUE(ID,ISP1,ISP2,IND,DIS &
     &                                   ,INTEGRAL,DINTEGRAL)
!     **************************************************************************
!     ** INTERPOLATE VALUES TO ACTUAL DISTANCE
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY  : OFFSITEX
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
          CALL ERROR$MSG('ARRAY OFFSITEX IN SIMPLELMTO_MODULE HAS NOT BEEN ALLOCATED')
          CALL ERROR$MSG('INTERNAL CODING ERROR')
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        NF=OFFSITEX(ISP1,ISP2)%NF
        IF(NF.GT.LF) THEN
          CALL ERROR$MSG('INTERNAL ARRAY SIZE EXCEEDED')
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        DISSAVE=DIS      
        G(:NF)=EXP(-OFFSITEX(ISP1,ISP2)%LAMBDA(:)*DIS)
      END IF
!
      IF(ID.EQ.'22') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%X22)) THEN
          CALL ERROR$MSG('NDDO TERM REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%X22(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%X22(:NF,IND)*G(:NF))
!
      ELSE IF (ID.EQ.'31') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%X31)) THEN
          CALL ERROR$MSG('3-1 EXCHANGE TERM REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%X31(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%X31(:NF,IND)*G(:NF))
PRINT*,'X31 EVAL A: ',INTEGRAL,DINTEGRAL
PRINT*,'X31 EVAL B: ',OFFSITEX(ISP1,ISP2)%X31
PRINT*,'X31 EVAL C: ',OFFSITEX(ISP1,ISP2)%LAMBDA
PRINT*,'X31 EVAL D: ',G(:NF)
!
      ELSE IF (ID.EQ.'BONDU') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%BONDU)) THEN
          CALL ERROR$MSG('U-TERM OF BONDOVERLAPS REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%BONDU(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%BONDU(:NF,IND)*G(:NF))
!
      ELSE IF (ID.EQ.'OV') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%OVERLAP)) THEN
          CALL ERROR$MSG('OFFSITE OVERLAP REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%OVERLAP(:NF,IND)*G(:NF))
        DINTEGRAL=-SUM(OFFSITEX(ISP1,ISP2)%LAMBDA(:) &
     &                *OFFSITEX(ISP1,ISP2)%OVERLAP(:NF,IND)*G(:NF))
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEXEVAL(TPARALLEL,NAT,NND,RBAS,R0,DENMAT &
     &                                  ,EX,STRESS,FORCE,HAMIL)
!     **************************************************************************
!     **  OFFSITE EXCHANGE ENERGY WITH PARTIAL WAVES                          **
!     **                                                                      **
!     **  OUTPUT IS: EX,HAMIL,FORCE,STRESS                                    **
!     **                                                                      **
!     **  THE HAMILTONIAN IS NOT RESET TO ZERO BUT THE TERMS FROM THIS        **
!     **  ROUTINE ARE ADDED ADDED.                                            **
!     **                                                                      **
!     **  PARALLELIZATION:                                                    **
!     **  TASKS ARE DISTRIBUTED ONTO NODES. RESULTS WILL BE SUMMED OVER ALL   **
!     **  TASKS IN THE CALLING ROUTINE.                                       **
!     **************************************************************************
      USE SIMPLELMTO_MODULE,ONLY : ISPECIES &
     &                            ,POTPAR &
     &                            ,NSP &
     &                            ,HYBRIDSETTING 
      USE RSPACEOP_MODULE  ,ONLY : RSPACEMAT_TYPE &
     &                            ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      LOGICAL(4)          ,INTENT(IN)   :: TPARALLEL
      INTEGER(4)          ,INTENT(IN)   :: NAT
      INTEGER(4)          ,INTENT(IN)   :: NND
      REAL(8)             ,INTENT(IN)   :: RBAS(3,3)
      REAL(8)             ,INTENT(IN)   :: R0(3,NAT)
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: DENMAT(NND)
      REAL(8)             ,INTENT(OUT)  :: EX
      REAL(8)             ,INTENT(OUT)  :: STRESS(3,3)
      REAL(8)             ,INTENT(OUT)  :: FORCE(3,NAT)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: HAMIL(NND)
      LOGICAL(4)          ,PARAMETER    :: TPR=.FALSE.
      LOGICAL(4)          ,PARAMETER    :: TOV=.FALSE.
      INTEGER(4)             :: INDM(NND)
      REAL(8)                :: RBASINV(3,3) ! 1/RBAS
      REAL(8)                :: DEDRBAS(3,3)
      INTEGER(4)             :: ISP
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4)             :: NN,NNA,NNB,NN2
      INTEGER(4)             :: ISPA,ISPB
      INTEGER(4)             :: LMN1A,LMN2A,LMN1B,LMN2B,LMN3A,LMN3B
      INTEGER(4)             :: LMNXA,LMNXB
      INTEGER(4)             :: LNXA,LNXB
      INTEGER(4),ALLOCATABLE :: NNAT(:)    !(NAT) POINTER TO ONSITE DENMAT
      INTEGER(4)             :: IAT,IATA,IATB
      REAL(8)                :: RA(3),RB(3)
      INTEGER(4)             :: NA,NB
      INTEGER(4)             :: NDIMD
      INTEGER(4)             :: I,L,LN,LM,LMN,LMX
      REAL(8)                :: SVAR,DSVAR,SVAR2
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
      REAL(8)   ,ALLOCATABLE :: DAB(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DBA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DAA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DBB(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HAB(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HBA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HAA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HBB(:,:,:)
      REAL(8)   ,ALLOCATABLE :: OV(:,:)
      REAL(8)   ,ALLOCATABLE :: DOV(:,:)
      REAL(8)                :: DEDD
      INTEGER(4)             :: LMRX
      INTEGER(4)             :: LRX
      INTEGER(4)             :: LX
      LOGICAL(4)             :: TNDDO,T31,TBONDX
      INTEGER(4)            :: THISTASK,NTASKS
INTEGER(4) :: J,K,I1,J1,K1
REAL(8) :: SVAR1
INTEGER(4) :: IX,IND1,IND2,IND3
REAL(8)    :: XDELTA,XSVAR
!     **************************************************************************
                            CALL TRACE$PUSH('SIMPLELMTO_OFFSITEXEVAL')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
PRINT*,'============ OFFSITEXEVAL ============================='
!
!     ==========================================================================
!     == COLLECT POINTERS NNAT TO ONSITE DENSITY MATRIX ELEMENTS              ==
!     ==========================================================================
      ALLOCATE(NNAT(NAT))
      NNAT=0    
      DO NN=1,NND
        IATA=DENMAT(NN)%IAT1
        IATB=DENMAT(NN)%IAT2
        IF(IATA.EQ.IATB.AND.SUM(DENMAT(NN)%IT**2).EQ.0) NNAT(IATA)=NN
      ENDDO
!
      DO IAT=1,NAT
        IF(NNAT(IAT).LE.0) THEN
          CALL ERROR$MSG('INDEX ARRAY FOR ONSITE DENSITY MATRIX')
          CALL ERROR$MSG('NO ONSITE TERMS FOUND FOR ATOM')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEX')
        END IF
      ENDDO
!
!     ==========================================================================
!     == FIND THE NEIGHBORLIST PAIRS WITH REVERSED INDICES                    ==
!     ==========================================================================
      INDM(:)=0
      DO NN=1,NND
        IATA=DENMAT(NN)%IAT1
        IATB=DENMAT(NN)%IAT2
        IF(INDM(NN).NE.0) CYCLE
        DO NN2=NN,NND
          IF(DENMAT(NN2)%IAT2.NE.IATA) CYCLE
          IF(DENMAT(NN2)%IAT1.NE.IATB) CYCLE
          IF(SUM((DENMAT(NN)%IT+DENMAT(NN2)%IT)**2).NE.0) CYCLE
          INDM(NN)=NN2          
          INDM(NN2)=NN          
          EXIT
        ENDDO
        IF(INDM(NN).EQ.0) THEN
          CALL ERROR$MSG('ERROR LOCATING INVERSE PAIR')
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXEVAL')
        END IF
      ENDDO
!
!     ==========================================================================
!     == LOOP OVER PAIRS                                                      ==
!     ==========================================================================
      EX=0.D0
      FORCE(:,:)=0.D0
      DEDRBAS(:,:)=0.D0
      DO NN=1,NND
        IF(TPARALLEL.AND.MOD(NN-1,NTASKS).NE.THISTASK-1) CYCLE
!
        IATA=DENMAT(NN)%IAT1
        IATB=DENMAT(NN)%IAT2
!       == CONSIDER ONLY OFFSITE TERMS =========================================
        IF(IATA.EQ.IATB.AND.SUM(DENMAT(NN)%IT**2).EQ.0) CYCLE
!!$CALL SIMPLELMTO$REPORTDENMAT(6)
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
        NNA=NNAT(IATA)
        NNB=NNAT(IATB)
        RA(:)=R0(:,IATA)
        RB(:)=R0(:,IATB)+MATMUL(RBAS,REAL(DENMAT(NN)%IT(:),KIND=8))
        DR(:)=RB(:)-RA(:)
        DIS=SQRT(SUM(DR**2))
!
        TNDDO=HYBRIDSETTING(ISPA)%TNDDO.AND.HYBRIDSETTING(ISPB)%TNDDO
        T31=HYBRIDSETTING(ISPA)%T31.AND.HYBRIDSETTING(ISPB)%T31
        TBONDX=HYBRIDSETTING(ISPA)%TBONDX.AND.HYBRIDSETTING(ISPB)%TBONDX
!       == DO NOTHING UNLESS NEEDED
        IF(.NOT.(TNDDO.OR.T31.OR.TBONDX)) CYCLE
CALL TIMING$CLOCKON('LOOP:OFFX')
!
!       ========================================================================
!       == BLOW UP DENSITY MATRIX                                             ==
!       ========================================================================
CALL TIMING$CLOCKON('OFFX:BLOWUP')
        LMNXA=POTPAR(ISPA)%LMNXH
        LMNXB=POTPAR(ISPB)%LMNXH
        NDIMD=DENMAT(NN)%N3
        IF(LMNXA.NE.DENMAT(NN)%N1.OR.LMNXB.NE.DENMAT(NN)%N2) THEN
          CALL ERROR$MSG('LMNX NOT EQUAL LMNXT') 
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXEVAL')
        END IF
        ALLOCATE(DAB(LMNXA,LMNXB,NDIMD))
        ALLOCATE(DBA(LMNXB,LMNXA,NDIMD))
        ALLOCATE(DAA(LMNXA,LMNXA,NDIMD))
        ALLOCATE(DBB(LMNXB,LMNXB,NDIMD))
        DAB(:,:,:) =DENMAT(NN )%MAT(:,:,:)
        DBA(:,:,:) =DENMAT(INDM(NN))%MAT(:,:,:)
        DAA(:,:,:)=DENMAT(NNA)%MAT(:,:,:)
        DBB(:,:,:)=DENMAT(NNB)%MAT(:,:,:)
CALL TIMING$CLOCKOFF('OFFX:BLOWUP')
 !
!       ========================================================================
!       == ROTATE DENSITY MATRIX SO THAT DISTANCE VECTOR POINTS IN Z-DIRECTION=
!       ========================================================================
CALL TIMING$CLOCKON('OFFX:ROTATE1')
!       == CONSTRUCT ROTATION MATRIX ===========================================
!       == DISTANCE VECTOR WILL BE NEW Z-DIRECTION =============================
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
        LMX=MAX(MAXVAL(POTPAR(ISPA)%LOXH),MAXVAL(POTPAR(ISPB)%LOXH))
        LMX=(LMX+1)**2
        ALLOCATE(YLMROT(LMX,LMX))
        CALL SPHERICAL$ROTATEYLM(LMX,ROT,YLMROT)
!
!       == CONSTRUCT ROTATION MATRIX OF TAILED ORBITALS ========================
        ALLOCATE(UROTA(LMNXA,LMNXA))
        ALLOCATE(UROTB(LMNXB,LMNXB))
        LNXA=POTPAR(ISPA)%LNXH
        UROTA(:,:)=0.D0
        LMN=1
        DO LN=1,LNXA
          L=POTPAR(ISPA)%LOXH(LN)
          LM=L**2+1
          UROTA(LMN:LMN+2*L,LMN:LMN+2*L)=YLMROT(LM:LM+2*L,LM:LM+2*L)
          LMN=LMN+2*L+1
        ENDDO
        LNXB=POTPAR(ISPB)%LNXH
        UROTB(:,:)=0.D0
        LMN=1
        DO LN=1,LNXB
          L=POTPAR(ISPB)%LOXH(LN)
          LM=L**2+1
          UROTB(LMN:LMN+2*L,LMN:LMN+2*L)=YLMROT(LM:LM+2*L,LM:LM+2*L)
          LMN=LMN+2*L+1
        ENDDO
        DEALLOCATE(YLMROT)
!
!       == ROTATE DENSITY MATRIX ===============================================
        !SUGGESTION: SPEED UP BY EXPLOITING THAT UROT IS SPARSE
        DO I=1,NDIMD
          DAA(:,:,I)=MATMUL(TRANSPOSE(UROTA),MATMUL(DAA(:,:,I),UROTA))
          DAB(:,:,I)=MATMUL(TRANSPOSE(UROTA),MATMUL(DAB(:,:,I),UROTB))
          DBA(:,:,I)=MATMUL(TRANSPOSE(UROTB),MATMUL(DBA(:,:,I),UROTA))
          DBB(:,:,I)=MATMUL(TRANSPOSE(UROTB),MATMUL(DBB(:,:,I),UROTB))
        ENDDO
CALL TIMING$CLOCKOFF('OFFX:ROTATE1')
!
!       ========================================================================
!       == ADD UP EXCHANGE ENERGY                                             ==
!       ========================================================================
        ALLOCATE(HAA(LMNXA,LMNXA,NDIMD))
        ALLOCATE(HAB(LMNXA,LMNXB,NDIMD))
        ALLOCATE(HBA(LMNXB,LMNXA,NDIMD))
        ALLOCATE(HBB(LMNXB,LMNXB,NDIMD))
        HAA=0.D0
        HAB=0.D0
        HBA=0.D0
        HBB=0.D0
        DEDD=0.D0  ! ENERGY DERIVATIVE WITH RESPECT TO DISTANCE
!
!       ========================================================================
!       == OVERLAP MATRIX ELEMENTS ARE NOT USED. THEY ARE EVALUATED FOR       ==
!       == TESTING PURPOSES. THE PARAMETER TOV SHALL NORMALLY BE FALSE.       ==
!       ========================================================================
        IF(TOV) THEN
          ALLOCATE(OV (LMNXA,LMNXB))
          ALLOCATE(DOV(LMNXA,LMNXB))
          CALL SIMPLELMTO_OFFSITEOVERLAP(ISPA,ISPB,DIS,LMNXA,LMNXB,OV,DOV)
PRINT*,'OV-REPORT OFFSITE ',IATA,IATB,OV,DOV
PRINT*,'OV-REPORT ONSITE A',POTPAR(ISPA)%OVERLAP
PRINT*,'OV-REPORT ONSITE B',POTPAR(ISPB)%OVERLAP
          DEALLOCATE(OV)
          DEALLOCATE(DOV)
        END IF
!
!       ========================================================================
!       == NDDO TERM:                                                         ==
!       == U22(1,2,3,4)=INT DX1 INT DX2: A1(X1)*A2(X1)][B3(X2)*B4(X2)]/|R1-R2|==
!       == LMN1A,LMN2A TIED TO COORDINATE X1, LMN1B,LMN2B TO X2               ==
!       ========================================================================
        IF(TNDDO) THEN
CALL TIMING$CLOCKON('OFFX:NDDO')
          ALLOCATE(U22   (LMNXA,LMNXA,LMNXB,LMNXB))
          ALLOCATE(DU22  (LMNXA,LMNXA,LMNXB,LMNXB))
          CALL SIMPLELMTO_OFFSITEX22U(ISPA,ISPB, DIS,LMNXA,LMNXB,U22,DU22)
          DO LMN1B=1,LMNXB
            DO LMN2B=1,LMNXB
              DO LMN2A=1,LMNXA
                DO LMN1A=1,LMNXA
!                 ==  W(LMN1A,LMN1B,LMN2A,LMN2B)
                  SVAR =-0.25D0* U22(LMN1A,LMN2A,LMN2B,LMN1B)
                  DSVAR=-0.25D0*DU22(LMN1A,LMN2A,LMN2B,LMN1B)
                  EX  =EX  + SVAR*SUM(DAB(LMN1A,LMN1B,:)*DBA(LMN2B,LMN2A,:))
                  DEDD=DEDD+DSVAR*SUM(DAB(LMN1A,LMN1B,:)*DBA(LMN2B,LMN2A,:))
                  HAB(LMN1A,LMN1B,:)=HAB(LMN1A,LMN1B,:)+SVAR*DBA(LMN2B,LMN2A,:)
                  HBA(LMN2B,LMN2A,:)=HBA(LMN2B,LMN2A,:)+SVAR*DAB(LMN1A,LMN1B,:)

!!$                  SVAR2=SUM(D(LMN1A,LMN1B,:)*D(LMN2A,LMN2B,:))
!!$PRINT*,'--NN--',DIS,SVAR,DSVAR,SVAR2,EX,DEDD,NN,IATA,IATB
!!$                  SVAR2=SUM(D(LMN1A,LMN1B,:)*D(LMN2A,LMN2B,:))
!!$                  EX  =EX  + SVAR*SUM(D(LMN1A,LMN1B,:)*D(LMN2A,LMN2B,:))
!!$                  DEDD=DEDD+DSVAR*SUM(D(LMN1A,LMN1B,:)*D(LMN2A,LMN2B,:))
!!$                  H(LMN1A,LMN1B,:)=H(LMN1A,LMN1B,:)+SVAR*D(LMN2A,LMN2B,:)
!!$                  H(LMN2A,LMN2B,:)=H(LMN2A,LMN2B,:)+SVAR*D(LMN1A,LMN1B,:)
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
          ALLOCATE(U3A1B (LMNXA,LMNXA,LMNXA,LMNXB))
          ALLOCATE(DU3A1B(LMNXA,LMNXA,LMNXA,LMNXB))
          CALL SIMPLELMTO_OFFSITEX31U(ISPA,ISPB, DIS,LMNXA,LMNXB,U3A1B,DU3A1B)
PRINT*,'X31REPORT A:',IATA,IATB,DIS,LMNXA,LMNXB
PRINT*,'X31REPORT A U=',U3A1B,DU3A1B
PRINT*,'X31REPORT A DAB =',DAB
PRINT*,'X31REPORT A DAA=',DAA
PRINT*,'X31REPORT A DBB=',DBB
         DO LMN1B=1,LMNXB
            DO LMN3A=1,LMNXA 
              DO LMN2A=1,LMNXA
                DO LMN1A=1,LMNXA
                  SVAR =-0.25D0* U3A1B(LMN1A,LMN2A,LMN3A,LMN1B)
                  DSVAR=-0.25D0*DU3A1B(LMN1A,LMN2A,LMN3A,LMN1B)
                  EX  =EX   +SVAR*SUM(DAB(LMN2A,LMN1B,:)*DAA(LMN1A,LMN3A,:))
PRINT*,'(1)',SVAR*SUM(DAB(LMN2A,LMN1B,:)*DAA(LMN1A,LMN3A,:))
                  DEDD=DEDD+DSVAR*SUM(DAB(LMN2A,LMN1B,:)*DAA(LMN1A,LMN3A,:))
                  HAA(LMN1A,LMN3A,:)=HAA(LMN1A,LMN3A,:)+SVAR*DAB(LMN2A,LMN1B,:)
                  HAB(LMN2A,LMN1B,:)=HAB(LMN2A,LMN1B,:)+SVAR*DAA(LMN1A,LMN3A,:)
                  EX  =EX   +SVAR*SUM(DAB(LMN2A,LMN1B,:)*DAA(LMN3A,LMN1A,:))
PRINT*,'(2)',SVAR*SUM(DAB(LMN2A,LMN1B,:)*DAA(LMN3A,LMN1A,:))
                  DEDD=DEDD+DSVAR*SUM(DAB(LMN2A,LMN1B,:)*DAA(LMN3A,LMN1A,:))
                  HAA(LMN3A,LMN1A,:)=HAA(LMN3A,LMN1A,:)+SVAR*DAB(LMN2A,LMN1B,:)
                  HAB(LMN2A,LMN1B,:)=HAB(LMN2A,LMN1B,:)+SVAR*DAA(LMN3A,LMN1A,:)
                ENDDO
              ENDDO
            ENDDO
          ENDDO               
          DEALLOCATE(U3A1B)
          DEALLOCATE(DU3A1B)
          ALLOCATE(U3B1A (LMNXB,LMNXB,LMNXB,LMNXA))
          ALLOCATE(DU3B1A(LMNXB,LMNXB,LMNXB,LMNXA))
          CALL SIMPLELMTO_OFFSITEX31U(ISPB,ISPA,-DIS,LMNXB,LMNXA,U3B1A,DU3B1A)
PRINT*,'X31REPORT B:',-DIS,LMNXB,LMNXA,U3B1A,DU3B1A
          DO LMN1A=1,LMNXA
            DO LMN3B=1,LMNXB 
              DO LMN2B=1,LMNXB
                DO LMN1B=1,LMNXB
                  SVAR =-0.25D0* U3B1A(LMN1B,LMN2B,LMN3B,LMN1A)
                  DSVAR=-0.25D0*DU3B1A(LMN1B,LMN2B,LMN3B,LMN1A)
                  EX=EX+SVAR*SUM(DAB(LMN1A,LMN2B,:)*DBB(LMN1B,LMN3B,:))
PRINT*,'(3)',SVAR*SUM(DAB(LMN1A,LMN2B,:)*DBB(LMN1B,LMN3B,:))
                  DEDD=DEDD+DSVAR*SUM(DAB(LMN1A,LMN2B,:)*DBB(LMN1B,LMN3B,:))
                  HBB(LMN1B,LMN3B,:)=HBB(LMN1B,LMN3B,:)+SVAR*DAB(LMN1A,LMN2B,:)
                  HAB(LMN1A,LMN2B,:)=HAB(LMN1A,LMN2B,:)+SVAR*DBB(LMN1B,LMN3B,:)
                  EX  =EX  + SVAR*SUM(DAB(LMN1A,LMN2B,:)*DBB(LMN3B,LMN1B,:))
PRINT*,'(4)',SVAR*SUM(DAB(LMN1A,LMN2B,:)*DBB(LMN3B,LMN1B,:))
                  DEDD=DEDD+DSVAR*SUM(DAB(LMN1A,LMN2B,:)*DBB(LMN3B,LMN1B,:))
                  HBB(LMN3B,LMN1B,:)=HBB(LMN3B,LMN1B,:)+SVAR*DAB(LMN1A,LMN2B,:)
                  HAB(LMN1A,LMN2B,:)=HAB(LMN1A,LMN2B,:)+SVAR*DBB(LMN3B,LMN1B,:)
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
          ALLOCATE(BONDU (LMNXA,LMNXB,LMNXA,LMNXB))
          ALLOCATE(DBONDU(LMNXA,LMNXB,LMNXA,LMNXB))
!         == BONDU(1,2,3,4)=INT DX INT DX': A1(X)*B2(X)][A3(X')*B4(X')]/|R-R'|==
          CALL SIMPLELMTO_OFFSITEXBONDU(ISPA,ISPB,DIS,LMNXA,LMNXB,BONDU,DBONDU)
          DO LMN2B=1,LMNXB
            DO LMN2A=1,LMNXA 
              DO LMN1B=1,LMNXB
                DO LMN1A=1,LMNXA
                  SVAR =-0.25D0* BONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                  DSVAR=-0.25D0*DBONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                  EX  =EX  + SVAR*SUM(DAA(LMN1A,LMN2A,:)*DBB(LMN1B,LMN2B,:))
                  DEDD=DEDD+DSVAR*SUM(DAA(LMN1A,LMN2A,:)*DBB(LMN1B,LMN2B,:)) 
                  HAA(LMN1A,LMN2A,:)=HAA(LMN1A,LMN2A,:)+SVAR*DBB(LMN1B,LMN2B,:)
                  HBB(LMN1B,LMN2B,:)=HBB(LMN1B,LMN2B,:)+SVAR*DAA(LMN1A,LMN2A,:)
                  SVAR =-0.25D0* BONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                  DSVAR=-0.25D0*DBONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                  EX  =EX  + SVAR*SUM(DAB(LMN1A,LMN2B,:)*DAB(LMN2A,LMN1B,:))
                  DEDD=DEDD+DSVAR*SUM(DAB(LMN1A,LMN2B,:)*DAB(LMN2A,LMN1B,:))
                  HAB(LMN1A,LMN2B,:)=HAB(LMN1A,LMN2B,:)+SVAR*DAB(LMN2A,LMN1B,:)
                  HAB(LMN2A,LMN1B,:)=HAB(LMN2A,LMN1B,:)+SVAR*DAB(LMN1A,LMN2B,:)
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
!       == ACCUMULATE FORCES                                                  ==
!       ========================================================================
        FORCE(:,IATB)=FORCE(:,IATB)-DEDD*DR(:)/DIS
        FORCE(:,IATA)=FORCE(:,IATA)+DEDD*DR(:)/DIS
        DO J=1,3
          DEDRBAS(:,J)=DEDRBAS(:,J)+REAL(DENMAT(NN)%IT(:))*DR(J)/DIS
        ENDDO
!
!       ========================================================================
!       == ROTATE HAMILTONIAN BACK                                            ==
!       ========================================================================
        DO I=1,NDIMD
          HAA(:,:,I)=MATMUL(UROTA,MATMUL(HAA(:,:,I),TRANSPOSE(UROTA)))
          HAB(:,:,I)=MATMUL(UROTA,MATMUL(HAB(:,:,I),TRANSPOSE(UROTB)))
          HBA(:,:,I)=MATMUL(UROTB,MATMUL(HBA(:,:,I),TRANSPOSE(UROTA)))
          HBB(:,:,I)=MATMUL(UROTB,MATMUL(HBB(:,:,I),TRANSPOSE(UROTB)))
        ENDDO
        DEALLOCATE(UROTA)
        DEALLOCATE(UROTB)
!
!       ========================================================================
!       == MAP HAMILTONIAN BACK                                               ==
!       ========================================================================
        HAMIL(NN )%MAT(:,:,:)=HAMIL(NN )%MAT(:,:,:)+HAB(:,:,:)
        HAMIL(INDM(NN))%MAT(:,:,:)=HAMIL(INDM(NN))%MAT(:,:,:)+HBA(:,:,:)
        HAMIL(NNA)%MAT(:,:,:)=HAMIL(NNA)%MAT(:,:,:)+HAA(:,:,:)
        HAMIL(NNB)%MAT(:,:,:)=HAMIL(NNB)%MAT(:,:,:)+HBB(:,:,:)
        DEALLOCATE(DAB)
        DEALLOCATE(DBA)
        DEALLOCATE(HAB)
        DEALLOCATE(HBA)
        DEALLOCATE(DAA)
        DEALLOCATE(HAA)
        DEALLOCATE(DBB)
        DEALLOCATE(HBB)
!!$WRITE(*,*)XDELTA*REAL(IX,8),EX,HAMIL(NN)%MAT(IND1,IND2,IND3)
!!$IF(IX.EQ.3) STOP 'FORCED'
!!$GOTO 1000
CALL TIMING$CLOCKOFF('LOOP:OFFX')
      ENDDO
!
!     ==========================================================================
!     == CONVERT DERIVATIVE INTO STRESSES                                     ==
!     ==========================================================================
      CALL LIB$INVERTR8(3,RBAS,RBASINV)
      DEDRBAS=MATMUL(DEDRBAS,RBASINV)
      STRESS=0.5D0*(DEDRBAS+TRANSPOSE(DEDRBAS))
!
!     ==========================================================================
!     == CLEAN UP TO AVOID MEMORY LEAK                                        ==
!     ==========================================================================
      DEALLOCATE(NNAT)
!
!     ==========================================================================
!     == DIAGNOSTIC PRINTOUT                                                  ==
!     ==========================================================================
      IF(TPR) THEN
        CALL SIMPLELMTO_TESTOFFSITEX()
        CALL ERROR$MSG('REGULAR STOP AFTER PRINTING DIAGNOSTIC INFORMATION')
        CALL ERROR$STOP('SIMPLELMTO_OFFSITEXEVAL')
      END IF
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_TESTOFFSITEX()
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE,ONLY : NSP &
     &                            ,POTPAR 
      IMPLICIT NONE
      INTEGER(4),PARAMETER   :: ND=100
      REAL(8)   ,PARAMETER   :: DMIN=0.2D0
      REAL(8)   ,PARAMETER   :: DMAX=5.D0
      REAL(8)                :: DIS
      REAL(8)   ,ALLOCATABLE :: U22(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: DU22(:,:,:,:)
      INTEGER(4)             :: ISPA,ISPB,ID
      INTEGER(4)             :: LMNXA,LMNXB
      INTEGER(4)             :: NFIL=8
!     **************************************************************************
      OPEN(NFIL,FILE='TESTOFFSITEX.DAT')
      DO ISPA=1,NSP
        DO ISPB=ISPA,NSP
          LMNXA=POTPAR(ISPA)%LMNXH
          LMNXB=POTPAR(ISPB)%LMNXH
          ALLOCATE(U22   (LMNXA,LMNXA,LMNXB,LMNXB))
          ALLOCATE(DU22  (LMNXA,LMNXA,LMNXB,LMNXB))
          U22=0.D0
          DU22=0.D0
          WRITE(NFIL,*)DMIN,U22,DU22
          DO ID=1,ND
            DIS=DMIN+REAL(ID-1,KIND=8)/REAL(ND-1,KIND=8)*(DMAX-DMIN)
            CALL SIMPLELMTO_OFFSITEX22U(ISPA,ISPB,DIS,LMNXA,LMNXB,U22,DU22)
            WRITE(NFIL,*)DIS,U22,DU22
          ENDDO
          U22=0.D0
          DU22=0.D0
          WRITE(NFIL,*)DMAX,U22,DU22
          DEALLOCATE(U22)
          DEALLOCATE(DU22)
        ENDDO
      ENDDO
      CALL ERROR$STOP('SIMPLELMTO_TESTOFFSITEX')
      STOP
      END

!
!*******************************************************************************
!*******************************************************************************
!****                         TWO-CENTER INTEGRALS                         *****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_TWOCENTER_TEST()
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
        CALL SIMPLELMTO_TWOCENTER(L1,M1,GID,NR,F,L2,M2,GID,NR,G,DIS,TOL,OVERLAP)
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
        CALL SIMPLELMTO_TWOCENTER(L1,M1,GID,NR,F,L2,M2,GID,NR,G,DIS,TOL,OVERLAP)
        SVAR=1.D0/315.D0*EXP(-4.D0*DIS)*(315.D0+4.D0*DIS*(315.D0+4.D0*DIS &
    &       *(75.D0+8.D0*DIS*(-15.D0+4.D0*DIS*(-9.D0-2.D0*DIS+8.D0*DIS**2)))))
        PRINT*,'DIS=',DIS,'OVERLAP ',OVERLAP,SVAR,OVERLAP-SVAR
      ENDDO
!
      CALL ERROR$MSG('REGULAR STOP AFTER PRINTING DIAGNOSTIC INFORMATION')
      CALL ERROR$STOP('TWO-CENTER INTEGRALS')
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE SIMPLELMTO_TWOCENTER_MODULE
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
ENDMODULE SIMPLELMTO_TWOCENTER_MODULE
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SIMPLELMTO_TWOCENTER(L1_,M1_,GID1_,NR1_,F1_ &
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
       USE SIMPLELMTO_TWOCENTER_MODULE, ONLY: &
      &             F1,F2 &
      &             ,PLMWORK1,PLMWORK2 &
      &             ,DIS &
      &             ,GID1,GID2 &
      &             ,L1,L2 &
      &             ,M1,M2 &
      &             ,NR1,NR2 &
      &             ,SLEN
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
         CALL ERROR$STOP('SIMPLELMTO$TWOCENTER')
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
       CALL NEWADAPT$EVALUATE(TOLERANCE,OVERLAP)
!
       DEALLOCATE(PLMWORK1)
       DEALLOCATE(PLMWORK2)
       DEALLOCATE(F1)
       DEALLOCATE(F2)
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SIMPLELMTO_TWOCENTER_MYFUNC(P,VALUE)
!      *************************************************************************
!      *************************************************************************
       USE SIMPLELMTO_TWOCENTER_MODULE, ONLY: DIS,SLEN &
      &                               ,PLMWORK1,PLMWORK2 &
      &                               ,GID1,GID2 &
      &                               ,NR1,NR2  &
      &                               ,L1,L2  &
      &                               ,M1,M2  &
      &                               ,F1,F2 
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
!      == ROUTINES ARE IN PSHERICAL OBJECT
       CALL PLGNDR((L1+1)**2,L1,COSTHETA1,PLMWORK1)
       CALL PLGNDR((L2+1)**2,L2,COSTHETA2,PLMWORK2)
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
MODULE NEWADAPT_MODULE
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
END MODULE NEWADAPT_MODULE
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE NEWADAPTINI()
!      *************************************************************************
!      ** INITIALIZE ADAPT. CALCULATE GRID POINTS AND WEIGHTS                 **
!      *************************************************************************
       USE NEWADAPT_MODULE, ONLY: TINI &
      &                    ,NP &
      &                    ,XY &
      &                    ,VALW &
      &                    ,ERRW &
      &                    ,DIRW
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
       SUBROUTINE NEWADAPT_BASICRULE(SEGMENT)
!      *************************************************************************
!      **  PERFORM GAUSSIAN QUADRATURE FOR ONE SEGMENT                        **
!      *************************************************************************
       USE NEWADAPT_MODULE, ONLY : SEGMENT_TYPE &
      &                           ,NP &
      &                           ,XY &
      &                           ,DIRW &
      &                           ,VALW &
      &                           ,ERRW 
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
         CALL NEWADAPT_INTEGRAND(P,VAL) 
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
       SUBROUTINE NEWADAPT_INTEGRAND(P,RES)
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
         CALL SIMPLELMTO_TWOCENTER_MYFUNC(P,RES)
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
         CALL ERROR$MSG('UNKNOWN TYPE')
         CALL ERROR$I4VAL('TYPE',TYPE)
         CALL ERROR$STOP('NEWADAPT_INTEGRAND')
       END IF
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE NEWADAPT$EVALUATE(TOLERANCE,VALUE)
!      *************************************************************************
!      ** ADAPTIVE, TWO-DIMENSIONAL INTEGRATION                               **
!      ** OF THE FUNCTION ADAPT$MYFUNCTION OVER THE AREA [-1,1]X[-1,1].       **
!      ** INTEGRAL IS EVALUATED TO AN ESTIMATED ACCURACY SMALLER THAN         **
!      ** THE SPECIFIED TOLERANCE.                                            **
!      **                                                                     **
!      ** SEE: A.C. GENZ AND A.A. MALIK, J. COMPUT. APPL. MATH. 6,295 (1980)  **
!      *******************************P.E. BLOECHL, GOSLAR 2011*****************
       USE NEWADAPT_MODULE, ONLY : SEGMENT_TYPE
USE LMTO_TWOCENTER_MODULE, ONLY : TPR
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN)  :: TOLERANCE ! REQUIRED ACCURACY
       REAL(8)   ,INTENT(OUT) :: VALUE     ! INTEGRAL VALUE
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
       CALL NEWADAPTINI()
!
!      =========================================================================
!      == INITIALIZE WITH FIRST SEGMENT                                       ==
!      =========================================================================
       NSEGMENTS=1
       STACK(1)%CENTER(:)=0.D0
       STACK(1)%WIDTH(:)=1.D0
       CALL NEWADAPT_BASICRULE(STACK(1))
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
         CALL NEWADAPT_BASICRULE(SEGMENT1)
         CALL NEWADAPT_BASICRULE(SEGMENT2)
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
       CALL ERROR$MSG('LOOP NOT CONVERGED')
       CALL ERROR$STOP('NEWADAPT$EVALUATE')
       END
