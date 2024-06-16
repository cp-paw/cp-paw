!*******************************************************************************
!*******************************************************************************
!**  INTERFACE FOR ELECTRONIC CORRELATIONS AND HYBRID FUNCTIONALS             **
!**  EXPRESSED IN TERMS OF LOCAL ORBITALS                                     **
!**                                                                           **
!**  SCHEMATIC FLOW DIAGRAM: AS SEEN FROM PAW_WAVES OBJECT                    **
!**     SIMPLELMTO$MAKE$POTPAR                                                **
!**                                                                           **
!**     WAVES$OFFSITEDENMAT               OSDENMAT   <- THIS%PROJ;            **
!**                                       (RE)ALLOCATE(OSDENMAT)              **
!**                                       DEALLOCATE(OSHAMIL)                 **
!**     SIMPLELMTO$MAKE$ETOT              ETOT,OSHAMIL <- OSDENMAT            **
!**         WAVES$GETI4('NND')                <- SIZE(OSDENMAT)               **
!**         WAVES$GETRSPACEMATA('DEMAT')     <- OSDENMAT                      **
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
  REAL(8)     :: LHFWEIGHT  ! LOCAL EXCHANGE TERMS ARE TREATED WITH LHFWEIGHT
  REAL(8)     :: RAUG       ! AUGMENTATION RADIUS
  INTEGER(4),POINTER :: NORBOFL(:) ! #(LOCAL ORBITALS PER ANGULAR MOMENTUM)
END TYPE HYBRIDSETTING_TYPE
!
!== HOLDS THE POTENTIAL PARAMETER FOR ONE ATOM TYPE ============================
TYPE POTPAR1_TYPE
  !__FROM SIMPLELMTO_MAKECHI1___________________________________________________
  REAL(8)             :: RAUG            ! MATCHING RADIUS
  INTEGER(4)          :: GID             ! GRID ID
  INTEGER(4)          :: GIDE            ! GRID ID FOR EXTENDED GRID (->OFFXINT)
  INTEGER(4)          :: LNXH            ! #(HEAD FUNCTIONS)
  INTEGER(4),POINTER  :: LOXH(:)         !(LNXH) MAIN ANGULAR MOMENTUM ="LOX"
  INTEGER(4),POINTER  :: LNOFH(:)        !(LNXH) PARTIAL WAVE ID
  INTEGER(4)          :: LNXT            ! #(TAIL FUNCTIONS)
  INTEGER(4),POINTER  :: LOXT(:)         !(LNXT) MAIN ANGULAR MOMENTUM ="LOX"
  REAL(8)   ,POINTER  :: AEPHIH(:,:)     !(NR,LNXH)
  REAL(8)   ,POINTER  :: AEPHIT(:,:)     !(NR,LNXH)
  REAL(8)   ,POINTER  :: AECHI(:,:)      !(NR,LNXH)
  REAL(8)   ,POINTER  :: PROPHIH(:,:)    !(LNX,LNXH)
  REAL(8)   ,POINTER  :: PROPHIT(:,:)    !(LNX,LNXT)
  REAL(8),ALLOCATABLE :: CMAT(:,:)       !(LNXH,LNXH) 
  !__FROM SIMPLELMTO_ONSITEMATRIXELEMENTS1______________________________________
  REAL(8),ALLOCATABLE :: ONSITEU(:,:,:,:) !(LMNXH,LMNXH,LMNXH,LMNXH)
  REAL(8),ALLOCATABLE :: QLN(:,:,:)       !(2,LNX,LNX)
  REAL(8),ALLOCATABLE :: PHIOV(:,:)       !(LNX,LNX)   <PHI|THETA_OMEGA|PHI>
  REAL(8),ALLOCATABLE :: OVERLAP(:,:)     !(LNXH,LNXH) <CHI|CHI>
END TYPE POTPAR1_TYPE
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
TYPE(POTPAR1_TYPE)  ,ALLOCATABLE :: POTPAR1(:) !POTENTIAL PARAMETERS (NEW)
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
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%ACTIVE=VAL
!
      ELSE IF(ID.EQ.'COREVALENCE') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TCV=VAL
!
      ELSE IF(ID.EQ.'FOCKSETUP') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TFOCKSETUP=VAL
!
      ELSE IF(ID.EQ.'NDDO') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%TNDDO=VAL
!
      ELSE IF(ID.EQ.'31') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETL4')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%T31=VAL
!
      ELSE IF(ID.EQ.'BONDX') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
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
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
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
     &                       ,POTPAR=>POTPAR1   !(NSP)
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
      USE SIMPLELMTO_MODULE, ONLY : ISPSELECTOR &
     &                             ,HYBRIDSETTING
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      INTEGER(4)  ,INTENT(IN)  :: LEN
      INTEGER(4)  ,INTENT(IN)  :: VAL(LEN)
!     **************************************************************************
      IF(ID.EQ.'NORBOFL') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
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
          HYBRIDSETTING(:)%LHFWEIGHT =0.D0
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
     &                             ,K2 &
     &                             ,SCREENL 
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     **************************************************************************
!
      IF(ID.EQ.'LHFWEIGHT') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%LHFWEIGHT=VAL
!
      ELSE IF(ID.EQ.'RAUG') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$SETR8')
        END IF
        HYBRIDSETTING(ISPSELECTOR)%RAUG=VAL
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
      USE SIMPLELMTO_MODULE, ONLY : HYBRIDSETTING &
     &                             ,ISPSELECTOR
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ID.EQ.'LHFWEIGHT') THEN
        IF(ISPSELECTOR.LE.0) THEN
          CALL ERROR$MSG('ATOM TYPE NOT SELECTED')
          CALL ERROR$MSG('SET VARIABLE "ISP" FIRST')
          CALL ERROR$I4VAL('ISPSELECTOR',ISPSELECTOR)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SIMPLELMTO$GETR8')
        END IF
        VAL=HYBRIDSETTING(ISPSELECTOR)%LHFWEIGHT

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
!     == CONSTRUCT LOCAL ORBITALS, HEAD AND TAIL FUNCTIONS                    ==
!     ==========================================================================
      CALL SIMPLELMTO_MAKECHI1()
!
!     ==========================================================================
!     == CONSTRUCT ONSITE MATRIX ELEMENTS WITH LOCAL ORBITALS                 ==
!     ==========================================================================
      CALL SIMPLELMTO_ONSITEMATRIXELEMENTS1()
!
!     ==========================================================================
!     == CONSTRUCT OFF-SITE MATRIX ELEMENTS                                   ==
!     ==========================================================================
      CALL SIMPLELMTO_OFFXINT()
!!$!
!!$!     =======================================================================
!!$!     ==  CONSTRUCT OFFSITE INTEGRALS OF TAILED ORBITALS                     
!!$!     =======================================================================
!!$      IF(TOFFSITE) THEN 
!!$                            CALL TIMING$CLOCKON('OFFSITE U-TENSOR')
!!$!       == GAUSSIAN FIT OF TAILED ORBITAL FOR ABAB-TYPE U-TENSOR ============
!!$                            CALL TIMING$CLOCKON('OFFS.U. GAUSSFIT')
!!$        CALL SIMPLELMTO_TAILEDGAUSSFIT()
!!$                            CALL TIMING$CLOCKOFF('OFFS.U. GAUSSFIT')
!!$!
!!$!       == OFF-SITE MATRIX ELEMENTS OF OVERLAP MATRIX AND U-TENSOR ==========
!!$!       == OVERLAP MATRIX ELEMENTS ARE ALWAYS COMPUTED, U-TENSOR ELEMENTS ===
!!$!       == ONLY WHEN REQUESTED ==============================================
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
      USE SIMPLELMTO_MODULE, ONLY : TON &
     &                             ,NSP &
     &                             ,HYBRIDSETTING &
     &                             ,K2 &
     &                             ,TOFFSITE &
     &                             ,SCREENL
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4)             :: THISTASK,NTASKS
      CHARACTER(32)          :: ID
      INTEGER(4)             :: ISP,I
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: LNX1
      INTEGER(4),ALLOCATABLE :: LOX1(:)
      REAL(8)                :: AEZ,RCOV
      REAL(8)                :: LHFWEIGHT
      REAL(8)                :: SCALERCUT
!     **************************************************************************
      IF(.NOT.TON) RETURN
      CALL SIMPLELMTO$INITIALIZE()
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
      CALL REPORT$TITLE(NFIL,'SIMPLELMTO OBJECT:  GENERIC VARIABLES')
      CALL WAVES$GETR8('SCALERCUT',SCALERCUT)
      CALL REPORT$R8VAL(NFIL,'RANGESCALE ',SCALERCUT,'*(RCOV(A)+RCOV(B))')
      CALL REPORT$R8VAL(NFIL,'K2 ',K2,'A.U.')
      CALL REPORT$R8VAL(NFIL,'SCREENING LENGTH',SCREENL,'A.U.')
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
      SUBROUTINE SIMPLELMTO_MAKECHI1()
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
     &                             ,POTPAR=>POTPAR1 &
     &                             ,NSP &
     &                             ,LNX &
     &                             ,LOX &
     &                             ,HYBRIDSETTING
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      REAL(8)                :: SVAR
      INTEGER(4)             :: L
!     **************************************************************************
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
      LOGICAL(4),PARAMETER   :: TH1S=.FALSE. ! OVERWRITE AECHI WITH H ORBITAL
      REAL(8)                :: RAUG       ! MATCHING RADIUS
      INTEGER(4)             :: LNXPHI     ! #(PARTIALWAVES)
      INTEGER(4)             :: LNXH       ! #(HEAD FUNCTIONS)
      INTEGER(4)             :: LNXT       ! #(TAIL FUNCTIONS)
      INTEGER(4),ALLOCATABLE :: LNOFH(:,:) !(2,LNXH)
      INTEGER(4),ALLOCATABLE :: LNOFT(:,:) !(2,LNXT)
      INTEGER(4),ALLOCATABLE :: LOXH(:)    !(LNXH) MAIN ANGULAR MOMENTA OF HEAD
      INTEGER(4),ALLOCATABLE :: LOXT(:)    !(LNXT) MAIN ANGULAR MOMENTA OF TAIL
      INTEGER(4)             :: GID        ! GRID ID
      INTEGER(4)             :: NR         ! #(RADIAL GRID POINTS)
      REAL(8)   ,ALLOCATABLE :: R(:)       !(NR) RADIAL GRID
      REAL(8)   ,ALLOCATABLE :: AUX(:)     !(NR) AUXILIARY ARRAY FOR RADIAL GRID
      REAL(8)   ,ALLOCATABLE :: NLPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIH(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIH(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHIH(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIT(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIT(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHIT(:,:)
      REAL(8)   ,ALLOCATABLE :: NLCHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AECHI(:,:)
      REAL(8)   ,ALLOCATABLE :: PSCHI(:,:)
      REAL(8)   ,ALLOCATABLE :: CMAT(:,:)
      REAL(8)   ,ALLOCATABLE :: PRO(:,:)
      REAL(8)   ,ALLOCATABLE :: PROPHIH(:,:)    !(LNXPHI,LNXH)
      REAL(8)   ,ALLOCATABLE :: PROPHIT(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHI1(:)
      REAL(8)   ,ALLOCATABLE :: AEPHI2(:)
      REAL(8)   ,ALLOCATABLE :: NLPHI1(:)
      REAL(8)   ,ALLOCATABLE :: NLPHI2(:)
      REAL(8)   ,ALLOCATABLE :: PSPHI1(:)
      REAL(8)   ,ALLOCATABLE :: PSPHI2(:)
      REAL(8)   ,ALLOCATABLE :: JTAIL(:,:) !(NR,LNXT)
      REAL(8)   ,ALLOCATABLE :: KHEAD(:,:) !(NR,LNXH)
      REAL(8)   ,ALLOCATABLE :: PHIOV(:,:) !<PHI|THETA_OMEGA|PHI>
      REAL(8)                :: DET
      REAL(8)                :: PHI1VAL,PHI1DER,PHI2VAL,PHI2DER
      REAL(8)                :: KVAL,KDER
      REAL(8)                :: JVAL,JDER
      REAL(8)                :: A1,A2
      REAL(8)                :: AEZ,RASA
      INTEGER(4)             :: IRAD   ! FIRST GRIDPOINT BEYOND RAUG
      CHARACTER(64)          :: STRING
      INTEGER(4)             :: ISP,LN,LNH,LNH2,LNT,LN1,LN2,IR
!     **************************************************************************
                             CALL TRACE$PUSH('SIMPLELMTO_MAKECHI')
      IF(.NOT.ALLOCATED(HYBRIDSETTING)) THEN
        CALL ERROR$MSG('HYBRIDSETTING NOT ALLOCATED')
        CALL ERROR$STOP('SIMPLELMTO_MAKECHI')
      END IF

      ALLOCATE(POTPAR(NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        LNXPHI=LNX(ISP)
!
!       == VARIABLE "STRING"  IS USED AS PREFIX FOR WAVE-FUNCTION FILES ========
        WRITE(STRING,*)ISP
        STRING='NEW_'//ADJUSTL(STRING) 
!
!       ========================================================================
!       ==  SELECT LOCAL ORBITALS AND CONNECT HEAD AND TAIL FUNCTIONS TO      ==
!       ==  PARTIAL WAVES                                                     ==
!       ========================================================================
        IF(.NOT.ASSOCIATED(HYBRIDSETTING(ISP)%NORBOFL)) THEN
          CALL ERROR$MSG('NO LOCAL ORBITALS SELECTED FOR THIS ATOM TYPE')
          CALL ERROR$I4VAL('ISP',ISP)
          CALL ERROR$STOP('SIMPLELMTO_MAKEPOTPAR1')
        END IF
!       == ONE TAIL FUNCTION FOR EACH MAIN ANGULAR MOMENTUM ====================
!       == ONE HEAD FUNCTION FOR EACH PARTIAL WAVE =============================
        LNXT=SIZE(HYBRIDSETTING(ISP)%NORBOFL)
        LNXH=SUM(HYBRIDSETTING(ISP)%NORBOFL)
!       == LNOFH POINTS TO THE TWO PARTIALWAVES, WHICH AUGMENT THE HEAD ========
!       == FOR LNOFH.LE.0, THE SCATTERING WAVE WITH L=-LNOFH IS USED ===========
        ALLOCATE(LNOFH(2,LNXH))
        ALLOCATE(LNOFT(2,LNXT))
        CALL SIMPLELMTO_LNHEADTAIL(HYBRIDSETTING(ISP)%NORBOFL &
    &                             ,LNX(ISP),LOX(:LNX(ISP),ISP) &
    &                             ,LNXH,LNOFH,LNXT,LNOFT)
!
!       == COLLECT ANGULAR MOMENTUM OF HEAD AND TAIL FUNCTIONS =================
        ALLOCATE(LOXH(LNXH))
        ALLOCATE(LOXT(LNXT))
        DO LNH=1,LNXH
          IF(LNOFH(1,LNH).GT.0) THEN
            LOXH(LNH)=LOX(LNOFH(1,LNH),ISP)
          ELSE
            LOXH(LNH)=-LNOFH(1,LNH)
          END IF
        ENDDO
        DO LNT=1,LNXT
          IF(LNOFT(1,LNT).GT.0) THEN
            LOXT(LNT)=LOX(LNOFT(1,LNT),ISP)
          ELSE
            LOXT(LNT)=-LNOFT(1,LNT)
          END IF
        ENDDO
!
!       ========================================================================
!       == MATCHING RADIUS                                                    ==
!       ========================================================================
        IF(HYBRIDSETTING(ISP)%RAUG.GE.0.D0) THEN
          RAUG=HYBRIDSETTING(ISP)%RAUG
        ELSE
          CALL SETUP$GETR8('AEZ',SVAR)
          CALL PERIODICTABLE$GET(SVAR,'R(COV)',RAUG)
        END IF
!
!       ========================================================================
!       ==  COLLECT RADIAL GRID
!       ========================================================================
!       == RADIAL GRID =========================================================
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        ALLOCATE(AUX(NR))
!
!       ========================================================================
!       ==  COLLECT PARTIAL WAVES                                             ==
!       ========================================================================
        ALLOCATE(AEPHI(NR,LNXPHI))
        ALLOCATE(AEPHIDOT(NR,LNXPHI))
        ALLOCATE(NLPHI(NR,LNXPHI))
        ALLOCATE(NLPHIDOT(NR,LNXPHI))
        ALLOCATE(PSPHI(NR,LNXPHI))
        ALLOCATE(PSPHIDOT(NR,LNXPHI))
!!$        CALL SETUP$GETR8A('NLPHI',NR*LNXPHI,NLPHI)
!!$        CALL SETUP$GETR8A('NLPHIDOT',NR*LNXPHI,NLPHIDOT)
!ATTENTION: QPHI AND NLPHI SEEM TO BE THE SAME, NAMELY NLPHI. 
!           THE SECOND PARTIAL WAVE SHOULD BE DIFFERENT
!
        CALL SETUP$GETR8A('QPHI',NR*LNXPHI,NLPHI)
        CALL SETUP$GETR8A('QPHIDOT',NR*LNXPHI,NLPHIDOT)
        CALL SETUP$GETR8A('AEPHI',NR*LNXPHI,AEPHI)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNXPHI,AEPHIDOT)
        CALL SETUP$GETR8A('PSPHI',NR*LNXPHI,PSPHI)
        CALL SETUP$GETR8A('PSPHIDOT',NR*LNXPHI,PSPHIDOT)
IF(TPR) THEN
  CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_AEPHI.DAT',GID,NR,LNXPHI,AEPHI)
  CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_PSPHI.DAT',GID,NR,LNXPHI,PSPHI)
  CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_NLPHI.DAT',GID,NR,LNXPHI,NLPHI)
  CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_AEPHIDOT.DAT',GID,NR,LNXPHI,AEPHIDOT)
  CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_PSPHIDOT.DAT',GID,NR,LNXPHI,PSPHIDOT)
  CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_NLPHIDOT.DAT',GID,NR,LNXPHI,NLPHIDOT)
END IF
!
!       ========================================================================
!       ==  CONSTRUCT SOLID HANKEL AND BESSEL FUNCTIONS                       ==
!       ========================================================================
        IF(TPR) THEN
          ALLOCATE(KHEAD(NR,LNXH))
          ALLOCATE(JTAIL(NR,LNXT))
          DO LNH=1,LNXH
            L=LOXH(LNH)
            DO IR=1,NR
              IF(R(IR).GT.1.D-1) THEN
                CALL LMTO$SOLIDHANKELRAD(L,R(IR),K2,KHEAD(IR,LNH),KDER)
              ELSE
                KHEAD(IR,LNH)=0.D0 ! AVOID SINGULARITY AT THE ORIGIN
              END IF
            ENDDO
          ENDDO
          DO LNT=1,LNXT
            L=LOXT(LNT)
            DO IR=1,NR
              CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,JTAIL(IR,LNT),JDER)
            ENDDO
          ENDDO
        END IF
!
!       ========================================================================
!       ==  CONSTRUCT HEAD AND TAIL FUNCTIONS                                 ==
!       ========================================================================
!       == CONSTRUCT PHIH SUCH THAT IT MATCHES SOLID HANKEL FUNCTION ===========
        ALLOCATE(AEPHIH(NR,LNXH))
        ALLOCATE(PSPHIH(NR,LNXH))
        ALLOCATE(NLPHIH(NR,LNXH))
!
        ALLOCATE(NLPHI1(NR))
        ALLOCATE(AEPHI1(NR))
        ALLOCATE(PSPHI1(NR))
        ALLOCATE(NLPHI2(NR))
        ALLOCATE(AEPHI2(NR))
        ALLOCATE(PSPHI2(NR))
        DO LNH=1,LNXH
          L=LOXH(LNH)
          LN1=LNOFH(1,LNH)
          LN2=LNOFH(2,LNH)
          NLPHI1=NLPHI(:,LN1)
          AEPHI1=AEPHI(:,LN1)
          PSPHI1=PSPHI(:,LN1)
          IF(LN2.GT.0) THEN
            NLPHI2=NLPHI(:,LN2)
            AEPHI2=AEPHI(:,LN2)
            PSPHI2=PSPHI(:,LN2)
          ELSE ! USE SCATTERING WAVE FUNCTION
            DO LNT=1,LNXT
              IF(LOXT(LNT).NE.L) CYCLE
              NLPHI2=NLPHIDOT(:,LNOFT(1,LNT))
              AEPHI2=AEPHIDOT(:,LNOFT(1,LNT))
              PSPHI2=PSPHIDOT(:,LNOFT(1,LNT))
              EXIT
            ENDDO
          END IF
!THIS IS TO BE CONSISTENT WITH PRIOR IMPLEMENTATION
!!$DO LN=1,LNXH
!!$  IF(LOXH(LN).NE.L)CYCLE
!!$  LN2=LNOFH(1,LN)
!!$  NLPHI2=NLPHIDOT(:,LN2)
!!$  AEPHI2=AEPHIDOT(:,LN2)
!!$  PSPHI2=PSPHIDOT(:,LN2)
!!$ENDDO

          CALL LMTO$SOLIDHANKELRAD(L,RAUG,K2,KVAL,KDER)
          CALL RADIAL$VALUE(GID,NR,NLPHI1,RAUG,PHI1VAL)
          CALL RADIAL$DERIVATIVE(GID,NR,NLPHI1,RAUG,PHI1DER)
          CALL RADIAL$VALUE(GID,NR,NLPHI2,RAUG,PHI2VAL)
          CALL RADIAL$DERIVATIVE(GID,NR,NLPHI2,RAUG,PHI2DER)
          DET=PHI1VAL*PHI2DER-PHI1DER*PHI2VAL
          A1=(KVAL*PHI2DER-KDER*PHI2VAL)/DET
          A2=-(KVAL*PHI1DER-KDER*PHI1VAL)/DET
          NLPHIH(:,LNH)=NLPHI1*A1+NLPHI2*A2
          AEPHIH(:,LNH)=AEPHI1*A1+AEPHI2*A2
          PSPHIH(:,LNH)=PSPHI1*A1+PSPHI2*A2
        ENDDO
!
!       == CONSTRUCT PHIT SUCH THAT IT MATCHES SOLID BESSEL FUNCTION ===========
        ALLOCATE(AEPHIT(NR,LNXT))
        ALLOCATE(PSPHIT(NR,LNXT))
        ALLOCATE(NLPHIT(NR,LNXT))
        DO LNT=1,LNXT
          L=LOXT(LNT)
          LN1=LNOFT(1,LNT)
          LN2=LNOFT(2,LNT)
          NLPHI1=NLPHI(:,LN1)
          AEPHI1=AEPHI(:,LN1)
          PSPHI1=PSPHI(:,LN1)
          IF(LN2.GT.0) THEN
            NLPHI2=NLPHI(:,LN2)
            AEPHI2=AEPHI(:,LN2)
            PSPHI2=PSPHI(:,LN2)
          ELSE ! USE SCATTERING WAVE FUNCTION
            NLPHI2=NLPHIDOT(:,LNOFT(1,LNT))
            AEPHI2=AEPHIDOT(:,LNOFT(1,LNT))
            PSPHI2=PSPHIDOT(:,LNOFT(1,LNT))
          END IF
!THIS IS TO BE CONSISTENT WITH PRIOR IMPLEMENTATION
!!$DO LN=1,LNXH
!!$  IF(LOXH(LN).NE.L)CYCLE
!!$  LN2=LNOFH(1,LN)
!!$  NLPHI2=NLPHIDOT(:,LN2)
!!$  AEPHI2=AEPHIDOT(:,LN2)
!!$  PSPHI2=PSPHIDOT(:,LN2)
!!$ENDDO
          CALL LMTO$SOLIDBESSELRAD(L,RAUG,K2,JVAL,JDER)
          CALL RADIAL$VALUE(GID,NR,NLPHI1,RAUG,PHI1VAL)
          CALL RADIAL$DERIVATIVE(GID,NR,NLPHI1,RAUG,PHI1DER)
          CALL RADIAL$VALUE(GID,NR,NLPHI2,RAUG,PHI2VAL)
          CALL RADIAL$DERIVATIVE(GID,NR,NLPHI2,RAUG,PHI2DER)
          DET=PHI1VAL*PHI2DER-PHI1DER*PHI2VAL
          A1=(JVAL*PHI2DER-JDER*PHI2VAL)/DET
          A2=-(JVAL*PHI1DER-JDER*PHI1VAL)/DET
          NLPHIT(:,LNT)=NLPHI1*A1+NLPHI2*A2
          AEPHIT(:,LNT)=AEPHI1*A1+AEPHI2*A2
          PSPHIT(:,LNT)=PSPHI1*A1+PSPHI2*A2
        ENDDO
        DEALLOCATE(NLPHI1)
        DEALLOCATE(AEPHI1)
        DEALLOCATE(PSPHI1)
        DEALLOCATE(NLPHI2)
        DEALLOCATE(AEPHI2)
        DEALLOCATE(PSPHI2)
!
!       ========================================================================
!       ==  CONSTRUCT HEAD-AUGMENTED SOLID HANKEL AND BESSEL FUNCTIONS        ==
!       ========================================================================
        IRAD=1
        DO IR=1,NR
          IF(R(IR).LE.RAUG) CYCLE
          IRAD=IR
          EXIT 
        ENDDO
        DO LNH=1,LNXH
          L=LOXH(LNH)
          DO IR=IRAD,NR   ! WILL BE AUGMENTED FROM 1 TO IRAD-1
            CALL LMTO$SOLIDHANKELRAD(L,R(IR),K2,KVAL,KDER)
            AEPHIH(IR,LNH)=AEPHIH(IR,LNH)-NLPHIH(IR,LNH)+KVAL
            PSPHIH(IR,LNH)=PSPHIH(IR,LNH)-NLPHIH(IR,LNH)+KVAL
            NLPHIH(IR,LNH)=KVAL ! DO THIS LAST!!!
          ENDDO
        ENDDO
!
        DO LNT=1,LNXT
          L=LOXT(LNT)
          DO IR=IRAD,NR   ! WILL BE AUGMENTED FROM 1 TO IRAD-1
            CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,JVAL,JDER)
            AEPHIT(IR,LNT)=AEPHIT(IR,LNT)-NLPHIT(IR,LNT)+JVAL
            PSPHIT(IR,LNT)=PSPHIT(IR,LNT)-NLPHIT(IR,LNT)+JVAL
            NLPHIT(IR,LNT)=JVAL ! DO THIS LAST!!!
          ENDDO
        ENDDO
!
!       == MAINTAINING BOTH, AECHI AND AEPHIH, IS UPERFLUOUS AND SHALL BE 
!       == CLEANED UP
        ALLOCATE(AECHI(NR,LNXH))
        ALLOCATE(PSCHI(NR,LNXH))
        ALLOCATE(NLCHI(NR,LNXH))
        AECHI=AEPHIH
        PSCHI=PSPHIH
        NLCHI=NLPHIH
!
!       ========================================================================
!       == PERFORM ON-SITE ORTHONORMALIZATION USING GRAM-SCHMIDT
!       ==  |CHI>=|K-AUG>*CMAT=|PHIH>-|PHIT>SDAGGER*CMAT
!       ========================================================================
        ALLOCATE(CMAT(LNXH,LNXH))
        CMAT(:,:)=0.D0
        DO LNH=1,LNXH
          CMAT(LNH,LNH)=1.D0
        ENDDO
        DO LNH=1,LNXH
!         ==  NORMALIZE ========================================================
          AUX=R(:)**2*AECHI(:,LNH)*AECHI(:,LNH)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          SVAR=1.D0/SQRT(SVAR)
          AECHI(:,LNH) =AECHI(:,LNH) *SVAR
          NLCHI(:,LNH) =NLCHI(:,LNH) *SVAR
          PSCHI(:,LNH) =PSCHI(:,LNH) *SVAR
          AEPHIH(:,LNH)=AEPHIH(:,LNH)*SVAR
          NLPHIH(:,LNH)=NLPHIH(:,LNH)*SVAR
          PSPHIH(:,LNH)=PSPHIH(:,LNH)*SVAR
          CMAT(:,LNH)  =CMAT(:,LNH)  *SVAR
!
!         == ORTHOGONALIZE HIGHER STATES =======================================
          DO LNH2=LNH+1,LNXH
            IF(LOXH(LNH2).NE.LOXH(LNH)) CYCLE          
            AUX=R(:)**2*AECHI(:,LNH)*AECHI(:,LNH2)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            AECHI(:,LNH2) =AECHI(:,LNH2) -AECHI(:,LNH) *SVAR
            NLCHI(:,LNH2) =NLCHI(:,LNH2) -NLCHI(:,LNH) *SVAR
            PSCHI(:,LNH2) =PSCHI(:,LNH2) -PSCHI(:,LNH) *SVAR
            AEPHIH(:,LNH2)=AEPHIH(:,LNH2)-AEPHIH(:,LNH)*SVAR
            NLPHIH(:,LNH2)=NLPHIH(:,LNH2)-NLPHIH(:,LNH)*SVAR
            PSPHIH(:,LNH2)=PSPHIH(:,LNH2)-PSPHIH(:,LNH)*SVAR
            CMAT(:,LNH2)  =CMAT(:,LNH2)  -CMAT(:,LNH)  *SVAR
          ENDDO
        ENDDO
!
!       ========================================================================
!       == PROJECTIONS
!       ========================================================================
        ALLOCATE(PRO(NR,LNXPHI))
        ALLOCATE(PROPHIH(LNXPHI,LNXH))
        ALLOCATE(PROPHIT(LNXPHI,LNXT))
        CALL SETUP$GETR8A('PRO',NR*LNXPHI,PRO)
        PROPHIH(:,:)=0.D0
        DO LNH=1,LNXH
          L=LOXH(LNH)
          DO LN=1,LNXPHI
            IF(LOX(LN,ISP).NE.L) CYCLE
!ATTENTION: REPLACE PSPHIH BY PSCHI
            AUX=R(:)**2*PRO(:,LN)*PSPHIH(:,LNH)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            PROPHIH(LN,LNH)=SVAR
          ENDDO
        ENDDO
        PROPHIT(:,:)=0.D0
        DO LNT=1,LNXT
          L=LOXT(LNT)
          DO LN=1,LNXPHI
            IF(LOX(LN,ISP).NE.L) CYCLE
            AUX=R(:)**2*PRO(:,LN)*PSPHIT(:,LNT)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            PROPHIT(LN,LNT)=SVAR
          ENDDO
        ENDDO
!
!       ========================================================================
!       == PARTIAL-WAVEOVERLAP <PHI|THETA_OMEGA|PHI>                          ==
!       ========================================================================
!       == <PHI|THETA|PHI>======================================================
        ALLOCATE(PHIOV(LNXPHI,LNXPHI))
        CALL SETUP$GETR8('AEZ',AEZ)  !SETUP OBJECT IS ALREADY SELECTED
        CALL PERIODICTABLE$GET(AEZ,'R(ASA)',RASA)
PRINT*,'ASA RADIUS FOR SPECIES ',ISP,' IS ', RASA
        PHIOV=0.D0
        DO LN1=1,LNXPHI
          L=LOX(LN1,ISP)
          DO LN2=LN1,LNXPHI
            IF(L.NE.LOX(LN2,ISP)) CYCLE
            CALL RADIAL$INTEGRATE(GID,NR,R**2*AEPHI(:,LN1)*AEPHI(:,LN2),AUX)
            CALL RADIAL$VALUE(GID,NR,AUX,RASA,SVAR)
            PHIOV(LN1,LN2)=SVAR
            PHIOV(LN2,LN1)=SVAR
          ENDDO
        ENDDO
!
!===============================================================================
!== OVERWRITE LOCAL ORBITAL FOR TESTING PURPOSES
!== USE A HYDROGEN SETUP WITH NPRO=1 .
!===============================================================================
IF(TH1S) THEN
  PRINT*,'WARNING: FUDGE FOR TESTINGIN SIMPLELMTO_MAKECHI1!!!! '
  AECHI(:,1)=EXP(-R)*2.D0   ! 2=1/SQRT(PI)*1/Y0
END IF
!
!       ========================================================================
!       == MAP TO POTPAR                                                     ==
!       ========================================================================
        POTPAR(ISP)%RAUG=RAUG
        POTPAR(ISP)%GID=GID    !RADIAL GRID IS NEEDED FOR THE LOCAL ORBITALS
        POTPAR(ISP)%LNXH=LNXH   
        ALLOCATE(POTPAR(ISP)%LOXH(LNXH))
        POTPAR(ISP)%LOXH=LOXH   
        POTPAR(ISP)%LNXT=LNXT   
        ALLOCATE(POTPAR(ISP)%LOXT(LNXT))
        POTPAR(ISP)%LOXT=LOXT
        ALLOCATE(POTPAR(ISP)%AEPHIH(NR,LNXH))
        POTPAR(ISP)%AEPHIH=AEPHIH
        ALLOCATE(POTPAR(ISP)%AEPHIT(NR,LNXT))
        POTPAR(ISP)%AEPHIT=AEPHIT
        ALLOCATE(POTPAR(ISP)%AECHI(NR,LNXH))
        POTPAR(ISP)%AECHI=AECHI
        ALLOCATE(POTPAR(ISP)%CMAT(LNXH,LNXH))
        POTPAR(ISP)%CMAT=CMAT
        ALLOCATE(POTPAR(ISP)%PROPHIT(LNXPHI,LNXT))
        POTPAR(ISP)%PROPHIT=PROPHIT
        ALLOCATE(POTPAR(ISP)%PROPHIH(LNXPHI,LNXH))
        POTPAR(ISP)%PROPHIH=PROPHIH
        ALLOCATE(POTPAR(ISP)%PHIOV(LNXPHI,LNXPHI))
        POTPAR(ISP)%PHIOV=PHIOV
 !
!       ========================================================================
!       === WRITE LOCAL ORBITALS TO FILE                                      ==
!       ========================================================================
        IF(TPR) THEN
          WRITE(STRING,*)ISP
          STRING='NEW_'//ADJUSTL(STRING) 
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_AECHI.DAT',GID,NR &
       &                    ,LNXH,AECHI)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_PSCHI.DAT',GID,NR &
       &                    ,LNXH,PSCHI)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_NLCHI.DAT',GID,NR &
       &                    ,LNXH,NLCHI)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_AEPHIH.DAT',GID,NR &
       &                    ,LNXH,AEPHIH)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_PSPHIH.DAT',GID,NR &
       &                    ,LNXH,PSPHIH)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_NLPHIH.DAT',GID,NR &
       &                    ,LNXH,NLPHIH)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_AEPHIT.DAT',GID,NR &
       &                    ,LNXT,AEPHIT)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_PSPHIT.DAT',GID,NR &
       &                    ,LNXT,PSPHIT)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_NLPHIT.DAT',GID,NR &
       &                    ,LNXT,NLPHIT)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_KHEAD.DAT',GID,NR &
       &                    ,LNXH,KHEAD)
          CALL SIMPLELMTO_WRITEPHI(TRIM(STRING)//'_JTAIL.DAT',GID,NR &
       &                    ,LNXT,JTAIL)
!
          WRITE(*,FMT='("ISP=",I3," LOXPHI=",20I5)')ISP,LOX(:LNXPHI,ISP)
          WRITE(*,FMT='("ISP=",I3," LOXH  =",20I5)')ISP,LOXH
          WRITE(*,FMT='("ISP=",I3," LOXT  =",20I5)')ISP,LOXT
!
          WRITE(*,FMT='(80("="),T20," <PRO|PHI-HEAD> FOR ISP=",I3," ")')ISP
          DO LN1=1,LNXPHI
            WRITE(*,FMT='(10F10.3)')PROPHIH(LN1,:)
          ENDDO
!
          WRITE(*,FMT='(80("="),T20," <PRO|PHI-TAIL> FOR ISP=",I3," ")')ISP
          DO LN1=1,LNXPHI
            WRITE(*,FMT='(10F10.3)')PROPHIT(LN1,:)
          ENDDO
!
          WRITE(*,FMT='(80("="),T20," CMAT FOR ISP=",I3," ")')ISP
          DO LN1=1,LNXH
            WRITE(*,FMT='(10F10.3)')CMAT(LN1,:)
          ENDDO
!
          WRITE(*,FMT='(80("="),T20," <PHI|THETA-OMEGA|PHI> FOR ISP=",I3," ")')ISP
          DO LN1=1,LNXPHI
            WRITE(*,FMT='(10F10.3)')PHIOV(LN1,:)
          ENDDO
        END IF
!
!       ========================================================================
!       == CLEAN UP                                                           ==
!       ========================================================================
        CALL SETUP$UNSELECT()
        DEALLOCATE(PRO)
        DEALLOCATE(AUX)
        DEALLOCATE(R)
        DEALLOCATE(NLPHI)
        DEALLOCATE(NLPHIDOT)
        DEALLOCATE(PSPHI)
        DEALLOCATE(PSPHIDOT)
        DEALLOCATE(AEPHI)
        DEALLOCATE(AEPHIDOT)
        DEALLOCATE(LNOFH)
        DEALLOCATE(LNOFT)
        DEALLOCATE(LOXH)
        DEALLOCATE(LOXT)
        DEALLOCATE(AEPHIH)
        DEALLOCATE(NLPHIH)
        DEALLOCATE(PSPHIH)
        DEALLOCATE(AEPHIT)
        DEALLOCATE(NLPHIT)
        DEALLOCATE(PSPHIT)
        DEALLOCATE(AECHI)
        DEALLOCATE(NLCHI)
        DEALLOCATE(PSCHI)
        DEALLOCATE(PROPHIH)
        DEALLOCATE(PROPHIT)
        DEALLOCATE(PHIOV)
        DEALLOCATE(CMAT)
        IF(TPR) THEN
          DEALLOCATE(KHEAD)
          DEALLOCATE(JTAIL)
        END IF
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_LNHEADTAIL(NORBOFL,LNX,LOX,LNXH,LNOFH,LNXT,LNOFT)
!     **************************************************************************
!     ** SELECT THE PARTIAL WAVES FOR THE CONSTRUCTION OF HEAD AND TAIL       **
!     ** FUNCTIONS.                                                           **
!     **                                                                      **
!     **    |PHIHEAD_J>=|PHI(LNOFH(1,J))> * A + |PHI(LNOFH(2,J)> * B          **
!     **    |PHITAIL_J>=|PHI(LNOFT(1,J))> * C + |PHI(LNOFT(2,J)> * D          **
!     **                                                                      **
!     ** LNXT IS THE NUMBER OF ANGULAR MOMENTA AND THERE IS ONE TAIL FUNCTION **
!     **      PER ANGULAR MOMENTUM (LNXT=SIZE(NORBOFL)                        **
!     ** LNXH IS THE NUMBER OF HEAD FUNCTIONS WHICH IS EQUAL TO THE NUMBER    **
!     **      OF PARTIAL WAVES (LNXH=SUM(NORBOFL))                            **
!     ** FOR EACH ANGULAR MOMENTUM, THERE IS A SINGLE SCATTERING PARTIAL WAVE **
!     ** FOR LNOFH.LE.0, LNOFH POINTS TO THAT SCATTERING WAVE WITH L=-LNOFH.  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LNXT     ! #(TAIL FUNCTIONS)
      INTEGER(4),INTENT(IN) :: NORBOFL(LNXT) !#(PARTIAL WAVES PER ANGULAR MOM.
      INTEGER(4),INTENT(IN) :: LNX      ! #(PARTIAL WAVES)
      INTEGER(4),INTENT(IN) :: LOX(LNX) ! MAIN ANG.MOM. FOR EACH PARTIAL WAVE
      INTEGER(4),INTENT(IN) :: LNXH     ! #(HEAD FUNCTIONS)
      INTEGER(4),INTENT(OUT):: LNOFH(2,LNXH)
      INTEGER(4),INTENT(OUT):: LNOFT(2,LNXT)
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      INTEGER(4)            :: LX
      INTEGER(4)            :: LN,LNH,L,N
!     **************************************************************************
      IF(LNXH.NE.SUM(NORBOFL)) THEN
        CALL ERROR$MSG('LNXH INCONSISTENT WITH NORBOFL')
        CALL ERROR$STOP('SIMPLELMTO_LNHEADTAIL')
      END IF
      LX=LNXT-1
!
!     ==========================================================================
!     == DETERMINE FIRST PARTIAL WAVE FOR HEAD FUNCTION                       ==
!     ==========================================================================
      LNH=0
      DO L=0,LX
        N=0
        DO LN=1,LNX
          IF(N.EQ.NORBOFL(L+1)) EXIT
          IF(LOX(LN).NE.L) CYCLE
          N=N+1
          LNH=LNH+1
          LNOFH(1,LNH)=LN
        ENDDO
        IF(N.LT.NORBOFL(L+1)) THEN
          CALL ERROR$MSG('LOCAL ORBITAL CONSTRUCTION FAILED')
          CALL ERROR$MSG('NUMBER OF LOCAL ORBITALS PER L EXCEEDS')
          CALL ERROR$MSG('NUMBER OF PARTIAL WAVES')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('TARGET N',NORBOFL(L+1))
          CALL ERROR$I4VAL('ACTUAL N',N)
          CALL ERROR$STOP('SIMPLELMTO_LNHEADTAIL')
        END IF
      ENDDO  
!
!     ==========================================================================
!     == DETERMINE SECOND PARTIAL WAVE FOR HEAD FUNCTION                      ==
!     ==========================================================================
      DO LNH=1,LNXH
        L=LOX(LNOFH(1,LNH))
        LNOFH(2,LNH)=-L   ! SCATTERING WAVE, UNLESS 2. PARTIAL WAVE EXISTS
        DO LN=LNOFH(1,LNH)+1,LNX      
          IF(LOX(LN).NE.L) CYCLE
          LNOFH(2,LNH)=LN
          EXIT
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE PARTIAL WAVES FOR TAIL FUNCTION                            ==
!     ==========================================================================
      LNH=0
      DO L=0,LX
        IF(NORBOFL(L+1).GT.0) THEN
          LNH=SUM(NORBOFL(:L+1))
          LNOFT(:,L+1)=LNOFH(:,LNH)
        ELSE
          LNOFT(2,L+1)=-L
          N=0
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            N=N+1
            LNOFT(N,L+1)=LN   ! IF THERE IS ONLY ONE PARTIAL WAVE, 
                              ! MAKE BOTH IDENTICAL
            IF(N.EQ.2) EXIT
          ENDDO
          IF(N.EQ.0) THEN
            CALL ERROR$MSG('LOCAL-ORBITAL CONSTRUCTION FAILED')
            CALL ERROR$MSG('IT ASSUMES AT LEAST ONE PARTIAL WAVE PER L')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$STOP('SIMPLELMTO_LNHEADTAIL')
          END IF
        END IF
      ENDDO
!
!     ==========================================================================
!     == REPORT                                                               ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(*,FMT='("LOX     =",10I5)')LOX
        WRITE(*,FMT='("NORBOFL =",10I5)')NORBOFL
        WRITE(*,FMT='("LNXH    =",10I5)')LNXH
        WRITE(*,FMT='("LNXT    =",10I5)')LNXT
        WRITE(*,FMT='("LNOFH(1)=",10I5)')LNOFH(1,:)
        WRITE(*,FMT='("LNOFH(2)=",10I5)')LNOFH(2,:)
        WRITE(*,FMT='("LNOFT(1)=",10I5)')LNOFT(1,:)
        WRITE(*,FMT='("LNOFT(2)=",10I5)')LNOFT(2,:)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_ONSITEMATRIXELEMENTS1()
!     **************************************************************************
!     ** CALCULATE ONSITE U-TENSOR AND MULTIPOLE MATRIX ELEMENTS              **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1 &
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
!
REAL(8),ALLOCATABLE :: R(:)
REAL(8),PARAMETER   :: PI=4.D0*ATAN(1.D0)
INTEGER(4)          :: I
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
!       ========================================================================
!       ==  U-TENSOR                                                          ==
!       ========================================================================
        ALLOCATE(ULITTLE(LRX+1,LNX,LNX,LNX,LNX))
!!$ALLOCATE(R(NR))
!!$CALL RADIAL$R(GID,NR,R)
!!$DO I=1,1
!!$  SCREENL=REAL(I)/0.529177D0
!!$  SCREENL=-1.
!!$ !POTPAR(ISP)%AECHI(:,1)=EXP(-R(:))/SQRT(PI)*SQRT(4.D0*PI)
        CALL SIMPLELMTO_ULITTLE(GID,NR,LRX,LNX,LOX,POTPAR(ISP)%AECHI &
     &                         ,SCREENL,ULITTLE)
        IF(.NOT.ALLOCATED(POTPAR(ISP)%ONSITEU)) THEN
          ALLOCATE(POTPAR(ISP)%ONSITEU(LMNX,LMNX,LMNX,LMNX))
        END IF
        CALL SIMPLELMTO_UTENSOR(LRX,LMNX,LNX,LOX,ULITTLE,POTPAR(ISP)%ONSITEU)
!!$  PRINT*,'POTPAR(ISP)%ONSITEU',SCREENL*0.529177D0,' AA U=',POTPAR(ISP)%ONSITEU
!!$ENDDO
!!$STOP 'FORCED'
        DEALLOCATE(ULITTLE)
!
!       ========================================================================
!       ==  ONSITE OVERLAP MATRIX                                             ==
!       == NEEDED ONLY FOR OUTPUT CAN BE REMOVED LATER                        ==
!       ========================================================================
        IF(.NOT.ALLOCATED(POTPAR(ISP)%OVERLAP)) THEN
          ALLOCATE(POTPAR(ISP)%OVERLAP(LMNX,LMNX))
        END IF
        CALL SIMPLELMTO_ONECENTEROVERLAP(GID,NR,LNX,LOX,POTPAR(ISP)%AECHI &
     &                                  ,LMNX,POTPAR(ISP)%OVERLAP)
!
!       ========================================================================
!       ==  DIPOLE MATRIX ELEMENTS                                            ==
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
            IF(SCREENL.LT.0.D0) THEN ! SCREENL<0 IMPLIES SCREENL=INFINITE
              CALL RADIAL$POISSON(GID,NR,L,RHO,POT)
            ELSE
!             __AVOID UNDERFLOW IN YUKAWA FOR TOO SMALL SCREENL_________________
              IF(SCREENL.LT.3.779D-2) THEN
                CALL ERROR$MSG('SET SCREENING LENGTH LARGER THAN 0.02 AA')
                CALL ERROR$MSG('IN !CONTROL!DFT!NTBO:SCREENL[AA]')
                CALL ERROR$MSG('(NUMERICAL PROBLEM IN LIBRARY ROUTINE)')
                CALL ERROR$STOP('SIMPLELMTO_ULITTLE')
              END IF
              SVAR=1.D0/MAX(1.D-8,SCREENL) !AVOID DIVIDE BY ZERO
              CALL RADIAL$YUKAWA(GID,NR,L,SVAR,RHO,POT) 
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
!     **  FOR A DENSITY CHI_LN(R)Y_L(R)CHI_LN'(|R|)*Y_L'(R)                   **
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
      INTEGER(4)            :: L1,L2
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
          IF(ABS(L1-L2).EQ.0) THEN
            AUX(:)=R(:)**2*CHI(:,LN1)*CHI(:,LN2)
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            SVAR=SQ4PI*SVAR
            QLN(1,LN1,LN2)=SVAR
            QLN(1,LN2,LN1)=SVAR
          ELSE IF(ABS(L1-L2).EQ.1) THEN
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
      IF(.NOT.TON) RETURN
                                    CALL TRACE$PUSH('SIMPLELMTO$ETOT')
      CALL SIMPLELMTO$INITIALIZE()
 
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
      CALL TRACE$PASS('SIMPLELMTO$ETOT DONE')
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
      USE SIMPLELMTO_MODULE, ONLY : TON
      USE RSPACEOP_MODULE, ONLY : RSPACEMAT_TYPE &
     &                           ,RSPACEOP$DELETE &
     &                           ,RSPACEOP$COPY &
     &                           ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      LOGICAL(4)          ,PARAMETER   :: TPR=.FALSE.
      CHARACTER(5)        ,PARAMETER   :: TYPE='FINAL'
!      CHARACTER(5)        ,PARAMETER   :: TYPE='PRIOR'
      LOGICAL(4)          ,PARAMETER   :: TTEST=.FALSE.
      LOGICAL(4)          ,PARAMETER   :: TTESTA=.FALSE.
      INTEGER(4)                       :: NAT
      INTEGER(4)                       :: NND
      REAL(8)                          :: RBAS(3,3),RBASINV(3,3)
      REAL(8)                          :: STRESS(3,3)
      REAL(8)             ,ALLOCATABLE :: R0(:,:) !(3,NAT)
      REAL(8)             ,ALLOCATABLE :: FORCE(:,:) !(3,NAT)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: DENMATPHI(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: DENMATCHI(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: HAMILCHI(:)
!     == OVERLAPCHI IS USED ONLY FOR TESTING PURPOSES
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: OVERLAPCHI(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: HAMILPHI(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: DONSITE(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: HONSITE(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: SBARE(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: PIPHI(:)
      REAL(8)             ,ALLOCATABLE :: FORCET(:,:) !(3,NAT)
      REAL(8)                          :: STRESST(3,3)
      REAL(8)                          :: ETOT
      INTEGER(4)                       :: I,IAT,J
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                     CALL TRACE$PUSH('SIMPLELMTO_HYBRID')
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      ALLOCATE(FORCE(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
PRINT*,'POSITIONS',R0
      FORCE(:,:)=0.D0
      STRESS(:,:)=0.D0
!
!     ==========================================================================
!     ==  COLLECT DENSITY MATRIX                                              ==
!     ==========================================================================
      CALL WAVES$GETI4('NND',NND)
      ALLOCATE(DENMATPHI(NND))
      CALL WAVES$GETRSPACEMATA('DENMAT',NND,DENMATPHI)
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'DENSITY MATRIX IN PHIS',NND,DENMATPHI)
!
!     ==========================================================================
!     ==  EXTRACT ON-SITE DENSITY MATRIX IN PARTIAL WAVE EXPANSION            ==
!     ==========================================================================
      ALLOCATE(DONSITE(NAT))
      CALL SIMPLELMTO_ONSITEDENMATEXTRACTADD('EXTRACT',NND,DENMATPHI &
     &                                                ,NAT,DONSITE)
!
!     ==========================================================================
!     ==  CALCULATE BARE STRUCTURECONSTANTS                                   ==
!     ==========================================================================
      ALLOCATE(SBARE(NND))
      DO I=1,NND
        SBARE(I)%IAT1=DENMATPHI(I)%IAT1
        SBARE(I)%IAT2=DENMATPHI(I)%IAT2
        SBARE(I)%IT=DENMATPHI(I)%IT
      ENDDO
      CALL SIMPLELMTO_STRUCTURECONSTANTS(NAT,R0,RBAS,NND,SBARE)
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'BARE STRUCTURE CONSTANTS',NND,SBARE)
!
!     ==========================================================================
!     == CONSTRUCT <PI|PHI> WHICH CONVERTS PARTIAL-WAVE PROJECTIONS INTO      ==
!     == LOCAL-ORBITAL PROJECTIONS
!     ==========================================================================
      ALLOCATE(PIPHI(NND))
      CALL SIMPLELMTO_MAKEPIPHI(NND,SBARE,PIPHI)
      DO I=1,NND
        CALL RSPACEOP$DELETE(SBARE(I))
      ENDDO
      DEALLOCATE(SBARE)
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'<PI|PHI> ',NND,PIPHI)
!
!     ==========================================================================
!     ==  CONVERT DENSITY MATRIX INTO THE LOCAL-ORBITAL BASIS                 ==
!     ==  DENMAT-PHI -> DENMAT-CHI
!     ==  DENMATPHI IS KEPTFORT FORCE CALCULATION
!     ==========================================================================
      ALLOCATE(DENMATCHI(NND))
      CALL SIMPLELMTO_DENMATPHITOCHI(NND,PIPHI,DENMATPHI,DENMATCHI)
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'DENMATCHI ',NND,DENMATCHI)
!
!     ==========================================================================
!     ==  CALCULATE ENERGY                                                    ==
!     ==========================================================================
      ETOT=0.D0
      ALLOCATE(HAMILCHI(NND))
      ALLOCATE(OVERLAPCHI(NND))
      DO I=1,NND
        CALL RSPACEOP$COPY(DENMATCHI(I),HAMILCHI(I))
        HAMILCHI(I)%MAT=0.D0
        OVERLAPCHI(I)%IAT1=DENMATCHI(I)%IAT1
        OVERLAPCHI(I)%IAT2=DENMATCHI(I)%IAT2
        OVERLAPCHI(I)%IT=DENMATCHI(I)%IT
        OVERLAPCHI(I)%N1=DENMATCHI(I)%N1
        OVERLAPCHI(I)%N2=DENMATCHI(I)%N2
        OVERLAPCHI(I)%N3=1
        ALLOCATE(OVERLAPCHI(I)%MAT(OVERLAPCHI(I)%N1,OVERLAPCHI(I)%N2,1))
        OVERLAPCHI(I)%MAT=0.D0
      ENDDO
      ALLOCATE(HONSITE(NAT))
      DO IAT=1,NAT
        CALL RSPACEOP$COPY(DONSITE(IAT),HONSITE(IAT))
        HONSITE(IAT)%MAT=0.D0
      ENDDO
      CALL TIMING$CLOCKON('SIMPLELMTO_HYBRIDENERGY')
      CALL SIMPLELMTO_HYBRIDENERGY(NAT,NND,RBAS,R0,DENMATCHI,DONSITE &
     &                            ,ETOT,STRESS,FORCE,HAMILCHI,HONSITE &
     &                            ,OVERLAPCHI)
      CALL TIMING$CLOCKOFF('SIMPLELMTO_HYBRIDENERGY')
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'OVERLAPCHI ',NND,OVERLAPCHI)
!
                                CALL TRACE$PASS('AFTER SIMPLELMTO_HYBRIDENERGY')
!
!     ==========================================================================
!     ==  FORCE CONTRIBUTION FROM CONVERSION BETWEEN CHI AND PHI              ==
!     ==========================================================================
      ALLOCATE(FORCET(3,NAT))
      CALL SIMPLELMTO_FORCEPHITOCHI(NAT,NND,DENMATPHI,PIPHI,HAMILCHI &
     &                             ,FORCET,RBAS,STRESST)
      FORCE=FORCE+FORCET
      STRESS=STRESS+STRESST
PRINT*,'FORCET',FORCET
      DEALLOCATE(FORCET)
!
!     ==========================================================================
!     ==  CONVERT HAMIL INTO THE PARTIAL-WAVE EXPANSION                       ==
!     ==========================================================================
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'HAMILTON IN CHI-S',NND,HAMILCHI)
      ALLOCATE(HAMILPHI(NND))
      CALL SIMPLELMTO_HAMILCHITOPHI(NND,PIPHI,HAMILCHI,HAMILPHI)
!
!     ==========================================================================
!     ==  ADD ONSITE HAMILTONIAN TO OFF-SITE HAMILTONIAN                      ==
!     ==========================================================================
      CALL SIMPLELMTO_ONSITEDENMATEXTRACTADD('ADDBACK',NND,HAMILPHI,NAT,HONSITE)
!
!     ==========================================================================
!     ==  ADD STRESS FROM FORCES                                              ==
!     ==========================================================================
!!$      CALL LIB$INVERTR8(3,RBAS,RBASINV)
!!$      DO I=1,3
!!$        DO J=1,3
!!$          STRESST(I,J)=SUM(R0(I,:)*FORCE(J,:))
!!$        ENDDO
!!$      ENDDO
!!$      STRESS=STRESS-MATMUL(RBASINV,STRESST)
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
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'HAMILTON IN PHIS',NND,HAMILPHI)
      CALL WAVES$SETRSPACEMATA('HAMIL',NND,HAMILPHI)
!
!     ==========================================================================
!     ==  CLEAN UP                                                            ==
!     ==========================================================================
      DO I=1,NND
        CALL RSPACEOP$DELETE(PIPHI(I))
        CALL RSPACEOP$DELETE(DENMATPHI(I))
        CALL RSPACEOP$DELETE(DENMATCHI(I))
        CALL RSPACEOP$DELETE(HAMILCHI(I))
        CALL RSPACEOP$DELETE(OVERLAPCHI(I))
        CALL RSPACEOP$DELETE(HAMILPHI(I))
      ENDDO
      DEALLOCATE(DENMATPHI)
      DEALLOCATE(DENMATCHI)
      DEALLOCATE(HAMILCHI)
      DEALLOCATE(OVERLAPCHI)
      DEALLOCATE(HAMILPHI)
      DEALLOCATE(PIPHI)
      DO IAT=1,NAT
        CALL RSPACEOP$DELETE(DONSITE(IAT))
        CALL RSPACEOP$DELETE(HONSITE(IAT))
      ENDDO
      DEALLOCATE(DONSITE)
      DEALLOCATE(HONSITE)
      CALL TRACE$PASS('SIMPLELMTO_HYBRID  DONE')
                                     CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_TESTDER()
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : TON
      USE RSPACEOP_MODULE, ONLY : RSPACEMAT_TYPE &
     &                           ,RSPACEOP$DELETE &
     &                           ,RSPACEOP$COPY &
     &                           ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      LOGICAL(4)          ,PARAMETER   :: TPR=.FALSE.
      INTEGER(4)                       :: NAT
      INTEGER(4)                       :: NND
      REAL(8)                          :: RBAS(3,3)
      REAL(8)                          :: STRESS(3,3)
      REAL(8)             ,ALLOCATABLE :: R0(:,:) !(3,NAT)
      REAL(8)             ,ALLOCATABLE :: DIS(:,:) !(3,NAT)
      REAL(8)             ,ALLOCATABLE :: RSAVE(:,:) !(3,NAT)
      REAL(8)             ,ALLOCATABLE :: FORCE(:,:) !(3,NAT)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: DENMATPHI(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: DENMATCHI(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: HAMILCHI(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: HAMILPHI(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: DONSITE(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: HONSITE(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: SBARE(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: SBARESAVE(:,:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: PIPHI(:)
      TYPE(RSPACEMAT_TYPE),ALLOCATABLE :: PIPHISAVE(:,:)
      REAL(8)                          :: RAN
      INTEGER(4)                       :: I,IAT,K,L,II,IDIS
      INTEGER(4)                       :: IAT1,IAT2
      INTEGER(4)                       :: N1,N2
!     **************************************************************************
      IF(.NOT.TON) RETURN
                                     CALL TRACE$PUSH('SIMPLELMTO_HYBRID')
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      ALLOCATE(RSAVE(3,NAT))
      ALLOCATE(DIS(3,NAT))
      ALLOCATE(FORCE(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      RSAVE=R0
      DO IAT=1,NAT
        DO I=1,3
          CALL RANDOM_NUMBER(RAN)
          DIS(I,IAT)=(2.D0*RAN-1.D0)*1.D-3
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  COLLECT DENSITY MATRIX                                              ==
!     ==========================================================================
      CALL WAVES$GETI4('NND',NND)
      ALLOCATE(DENMATPHI(NND))
      CALL WAVES$GETRSPACEMATA('DENMAT',NND,DENMATPHI)
!      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'DENSITY MATRIX IN PHIS',NND,DENMATPHI)
!
!     ==========================================================================
!     ==  EXTRACT ON-SITE DENSITY MATRIX IN PARTIAL WAVE EXPANSION            ==
!     ==========================================================================
      ALLOCATE(DONSITE(NAT))
      CALL SIMPLELMTO_ONSITEDENMATEXTRACTADD('EXTRACT',NND,DENMATPHI &
     &                                                ,NAT,DONSITE)
!
!     =========================================================================
!     =========================================================================
!     =========================================================================
      ALLOCATE(SBARESAVE(NND,-5:5))
      ALLOCATE(PIPHISAVE(NND,-5:5))
      ALLOCATE(SBARE(NND))
      ALLOCATE(PIPHI(NND))
      DO IDIS=-5,5
        R0=RSAVE+DIS*REAL(IDIS,KIND=8)
!
!       ========================================================================
!       ==  CALCULATE BARE STRUCTURECONSTANTS                                 ==
!       ========================================================================
        DO I=1,NND
          SBARE(I)%IAT1=DENMATPHI(I)%IAT1
          SBARE(I)%IAT2=DENMATPHI(I)%IAT2
          SBARE(I)%IT  =DENMATPHI(I)%IT
        ENDDO
        CALL SIMPLELMTO_STRUCTURECONSTANTS(NAT,R0,RBAS,NND,SBARE)
!
        DO I=1,NND
          IAT1=SBARE(I)%IAT1
          IAT2=SBARE(I)%IAT2
          N1=SBARE(I)%N1
          N2=SBARE(I)%N2
          SBARESAVE(I,IDIS)%IAT1=IAT1
          SBARESAVE(I,IDIS)%IAT2=IAT2
          SBARESAVE(I,IDIS)%N1=N1
          SBARESAVE(I,IDIS)%N2=N2
          SBARESAVE(I,IDIS)%N3=1
          ALLOCATE(SBARESAVE(I,IDIS)%MAT(N1,N2,1))
          IF(IDIS.NE.0) THEN
            SBARESAVE(I,IDIS)%MAT(:,:,1)=SBARE(I)%MAT(:,:,1)
          ELSE
            SBARESAVE(I,IDIS)%MAT=0.D0
            DO II=1,3
              SBARESAVE(I,IDIS)%MAT(:,:,1) &
    &                =SBARESAVE(I,IDIS)%MAT(:,:,1)&
    &                +SBARE(I)%MAT(:,:,1+II)*(DIS(II,IAT2)-DIS(II,IAT1))
            ENDDO
!            WRITE(*,*)I,SBARESAVE(I,IDIS)%MAT(:,:,1)
          END IF
        ENDDO
!
!       ========================================================================
!       == CONSTRUCT <PI|PHI> WHICH CONVERTS PARTIAL-WAVE PROJECTIONS INTO    ==
!       == LOCAL-ORBITAL PROJECTIONS
!       ========================================================================
        CALL SIMPLELMTO_MAKEPIPHI(NND,SBARE,PIPHI)
!
        DO I=1,NND
          IAT1=PIPHI(I)%IAT1
          IAT2=PIPHI(I)%IAT2
          N1=PIPHI(I)%N1
          N2=PIPHI(I)%N2
          PIPHISAVE(I,IDIS)%IAT1=IAT1
          PIPHISAVE(I,IDIS)%IAT2=IAT2
          PIPHISAVE(I,IDIS)%N1=N1
          PIPHISAVE(I,IDIS)%N2=N2
          PIPHISAVE(I,IDIS)%N3=1
          ALLOCATE(PIPHISAVE(I,IDIS)%MAT(N1,N2,1))
          IF(IDIS.NE.0) THEN
            PIPHISAVE(I,IDIS)%MAT(:,:,1)=PIPHI(I)%MAT(:,:,1)
          ELSE
            PIPHISAVE(I,IDIS)%MAT=0.D0
            DO II=1,3
              PIPHISAVE(I,IDIS)%MAT(:,:,1) &
    &                =PIPHISAVE(I,IDIS)%MAT(:,:,1)&
    &                +PIPHI(I)%MAT(:,:,1+II)*(DIS(II,IAT2)-DIS(II,IAT1))
            ENDDO
!            WRITE(*,*)I,PIPHISAVE(I,IDIS)%MAT(:,:,1)
          END IF
        ENDDO
!
!       =======================================================================
!       ==
!       =======================================================================
        DO I=1,NND
          DEALLOCATE(SBARE(I)%MAT)
          DEALLOCATE(PIPHI(I)%MAT)
        ENDDO
      ENDDO
      DEALLOCATE(SBARE)
!
!     ==========================================================================
!     ==  ANALYZE
!     ==========================================================================
PRINT*,'ANALYZE SBARE GRADIENT'
      DO I=1,NND
         DO IDIS=1,5
           SBARESAVE(I,IDIS)%MAT(:,:,1) &
     &          =(SBARESAVE(I,IDIS)%MAT(:,:,1)-SBARESAVE(I,-IDIS)%MAT(:,:,1)) &
     &           /REAL(2*IDIS,KIND=8) &
     &           /SBARESAVE(I,0)%MAT(:,:,1)
         ENDDO
         DO K=1,SBARESAVE(I,0)%N1
           DO L=1,SBARESAVE(I,0)%N2
             IF(SBARESAVE(I,0)%MAT(K,L,1).EQ.0.D0) CYCLE
             WRITE(*,FMT='(5F10.5)')(SBARESAVE(I,IDIS)%MAT(K,L,1),IDIS=1,5)
           ENDDO
         ENDDO
      ENDDO      
PRINT*,'ANALYZE PIPHI GRADIENT'
      DO I=1,NND
         DO IDIS=1,5
           PIPHISAVE(I,IDIS)%MAT(:,:,1) &
     &          =(PIPHISAVE(I,IDIS)%MAT(:,:,1)-PIPHISAVE(I,-IDIS)%MAT(:,:,1)) &
     &           /REAL(2*IDIS,KIND=8) &
     &           /PIPHISAVE(I,0)%MAT(:,:,1)
         ENDDO
         DO K=1,PIPHISAVE(I,0)%N1
           DO L=1,PIPHISAVE(I,0)%N2
             IF(PIPHISAVE(I,0)%MAT(K,L,1).EQ.0.D0) CYCLE
             WRITE(*,FMT='(5F10.5)')(PIPHISAVE(I,IDIS)%MAT(K,L,1),IDIS=1,5)
           ENDDO
         ENDDO
      ENDDO      
STOP 'FORCED --'
!
!     ==========================================================================
!     ==  CLEAN UP                                                            ==
!     ==========================================================================
        DEALLOCATE(SBARE)
      DO I=1,NND
        CALL RSPACEOP$DELETE(PIPHI(I))
        CALL RSPACEOP$DELETE(DENMATPHI(I))
        CALL RSPACEOP$DELETE(DENMATCHI(I))
        CALL RSPACEOP$DELETE(HAMILCHI(I))
        CALL RSPACEOP$DELETE(HAMILPHI(I))
      ENDDO
      DEALLOCATE(DENMATPHI)
      DEALLOCATE(DENMATCHI)
      DEALLOCATE(HAMILCHI)
      DEALLOCATE(HAMILPHI)
      DEALLOCATE(PIPHI)
      DO IAT=1,NAT
        CALL RSPACEOP$DELETE(DONSITE(IAT))
        CALL RSPACEOP$DELETE(HONSITE(IAT))
      ENDDO
      DEALLOCATE(DONSITE)
      DEALLOCATE(HONSITE)
                                     CALL TRACE$POP()
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
            CALL ERROR$I4VAL('IND1-IT(1)',DENMAT(IND1)%IT(1))
            CALL ERROR$I4VAL('IND1-IT(2)',DENMAT(IND1)%IT(2))
            CALL ERROR$I4VAL('IND1-IT(3)',DENMAT(IND1)%IT(3))
            CALL ERROR$I4VAL('IND2-IAT1',DENMAT(IND2)%IAT1)
            CALL ERROR$I4VAL('IND2-IAT2',DENMAT(IND2)%IAT2)
            CALL ERROR$I4VAL('IND2-IT(1)',DENMAT(IND2)%IT(1))
            CALL ERROR$I4VAL('IND2-IT(2)',DENMAT(IND2)%IT(2))
            CALL ERROR$I4VAL('IND2-IT(3)',DENMAT(IND2)%IT(3))
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
!!$      DO IND1=1,NND
!!$        TCHK=.FALSE.
!!$        DO IND2=1,NND
!!$          IF(HAMIL(IND2)%IAT2.NE.HAMIL(IND1)%IAT1) CYCLE
!!$          IF(HAMIL(IND2)%IAT1.NE.HAMIL(IND1)%IAT2) CYCLE
!!$          IF(SUM((HAMIL(IND2)%IT+HAMIL(IND1)%IT)**2).NE.0) CYCLE
!!$          IF(TCHK) THEN
!!$            CALL ERROR$MSG('ERROR 1')
!!$            CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
!!$          END IF
!!$          TCHK=.TRUE.
!!$          SVAR=0.D0
!!$          DO J=1,HAMIL(IND1)%N3
!!$            HAMIL(IND2)%MAT(:,:,J)=0.5D0*(HAMIL(IND2)%MAT(:,:,J) &
!!$     &                         +TRANSPOSE(HAMIL(IND1)%MAT(:,:,J)))
!!$            HAMIL(IND1)%MAT(:,:,J)=TRANSPOSE(HAMIL(IND2)%MAT(:,:,J))
!!$          ENDDO
!!$          IF(SVAR.GT.1.D-10) THEN
!!$            CALL RSPACEOP$WRITEMAT(6,'HAMIL',NND,HAMIL)
!!$            CALL ERROR$MSG('ERROR 2')
!!$            CALL ERROR$I4VAL('IND1',IND1)
!!$            CALL ERROR$I4VAL('IND2',IND2)
!!$            CALL ERROR$I4VAL('IND1-IAT1',HAMIL(IND1)%IAT1)
!!$            CALL ERROR$I4VAL('IND1-IAT2',HAMIL(IND1)%IAT2)
!!$            CALL ERROR$I4VAL('IND1-IT',HAMIL(IND1)%IT)
!!$            CALL ERROR$I4VAL('IND2-IAT1',HAMIL(IND2)%IAT1)
!!$            CALL ERROR$I4VAL('IND2-IAT2',HAMIL(IND2)%IAT2)
!!$            CALL ERROR$I4VAL('IND2-IT',HAMIL(IND2)%IT)
!!$            CALL ERROR$R8VAL('DEV',SVAR)
!!$            CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
!!$          END IF
!!$        ENDDO
!!$        IF(.NOT.TCHK) THEN
!!$          CALL ERROR$MSG('ERROR 3')
!!$          CALL ERROR$STOP('SIMPLELMTO_FAKEETOT')
!!$        END IF
!!$      ENDDO
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
            CALL ERROR$I4VAL('IND1-IT(1)',HAMIL(IND1)%IT(1))
            CALL ERROR$I4VAL('IND1-IT(2)',HAMIL(IND1)%IT(2))
            CALL ERROR$I4VAL('IND1-IT(3)',HAMIL(IND1)%IT(3))
            CALL ERROR$I4VAL('IND2-IAT1',HAMIL(IND2)%IAT1)
            CALL ERROR$I4VAL('IND2-IAT2',HAMIL(IND2)%IAT2)
            CALL ERROR$I4VAL('IND2-IT(1)',HAMIL(IND2)%IT(1))
            CALL ERROR$I4VAL('IND2-IT(2)',HAMIL(IND2)%IT(2))
            CALL ERROR$I4VAL('IND2-IT(3)',HAMIL(IND2)%IT(3))
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
      SUBROUTINE SIMPLELMTO_INDBACKON(NND,DENMAT,NAT,INDON,INDBACK)
!     **************************************************************************
!     **                                                                      **
!     ********************************************P. BLOECHL, GOSLAR 2019*******
      USE RSPACEOP_MODULE  , ONLY : RSPACEMAT_TYPE 
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN)   :: NAT         ! #(ATOMS)
      INTEGER(4)          ,INTENT(IN)   :: NND         ! #(NEIGHBORLIST ENTRIES)
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: DENMAT(NND)
      INTEGER(4)          ,INTENT(OUT)  :: INDON(NAT)  ! POINTER TO ONSITE TERM
      INTEGER(4)          ,INTENT(OUT)  :: INDBACK(NND)! POINTER TO BACK HOP
      INTEGER(4)                        :: IAT1,IAT2
      INTEGER(4)                        :: I,J
!     **************************************************************************
                                    CALL TRACE$PUSH('SIMPLELMTO_INDBACKON')
!
!     ==========================================================================
!     == INDBACK: POINTER TO BACK HOP
!     == INDON POINTS FROM EACH ATOM TO THE CORRESPONDIG ONSITE TERM ===========
!     ==========================================================================
      INDBACK=0
      DO I=1,NND
        IF(INDBACK(I).NE.0) CYCLE
        IAT1=DENMAT(I)%IAT1
        IAT2=DENMAT(I)%IAT2
        DO J=I,NND
          IF(DENMAT(J)%IAT2.NE.IAT1) CYCLE
          IF(DENMAT(J)%IAT1.NE.IAT2) CYCLE
          IF(SUM((DENMAT(I)%IT+DENMAT(J)%IT)**2).NE.0) CYCLE
          INDBACK(I)=J
          INDBACK(J)=I
          IF(I.EQ.J) INDON(IAT1)=I
          EXIT
        ENDDO
      ENDDO
                                    CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_DENMATPHITOCHI(NND,PIPHI,DENMATPHI,DENMATCHI)
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
     &                             ,POTPAR=>POTPAR1 &
     &                             ,LOX &
     &                             ,LNX &
     &                             ,NSP
      USE RSPACEOP_MODULE  , ONLY : RSPACEMAT_TYPE &
     &                             ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN)   :: NND
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: DENMATPHI(NND)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: DENMATCHI(NND)
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: PIPHI(NND)
      LOGICAL(4)          ,PARAMETER    :: TPR=.FALSE.
      INTEGER(4)          ,ALLOCATABLE  :: INDBACK(:)  !(NND)
      INTEGER(4)          ,ALLOCATABLE  :: INDON(:)    !(NAT)
      INTEGER(4)                        :: IAT1,IAT2
      INTEGER(4)                        :: ISP,ISP1,ISP2
      INTEGER(4)                        :: N1,N2,N3
      INTEGER(4)                        :: NAT
      INTEGER(4)                        :: I,J
!     **************************************************************************
                                    CALL TRACE$PUSH('SIMPLELMTO_DENMATPHITOCHI')
!
!     ==========================================================================
!     == POINTER TO BACK HOP
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(INDON(NAT))
      ALLOCATE(INDBACK(NND))
      CALL SIMPLELMTO_INDBACKON(NND,DENMATPHI,NAT,INDON,INDBACK)
!
!     ==========================================================================
!     ==  FORWARD TRANSFORMATION FROM RHOPHI TO RHOCHI                        ==
!     ==========================================================================
!       ========================================================================
!       ==  COPY DENMATPHI INTO DENMATCHI
!       ========================================================================
        DO I=1,NND
          IAT1=DENMATPHI(I)%IAT1
          IAT2=DENMATPHI(I)%IAT2
          ISP1=ISPECIES(IAT1)
          ISP2=ISPECIES(IAT2)
          N1  =DENMATPHI(I)%N1
          N2  =DENMATPHI(I)%N2
          N3  =DENMATPHI(I)%N3
          DENMATCHI(I)%IAT1=IAT1
          DENMATCHI(I)%IAT2=IAT2
          DENMATCHI(I)%IT  =DENMATPHI(I)%IT
          N1=SUM(2*POTPAR(ISP1)%LOXH+1)
          N2=SUM(2*POTPAR(ISP2)%LOXH+1)
          N3=DENMATPHI(I)%N3
          DENMATCHI(I)%N1=N1
          DENMATCHI(I)%N2=N2
          DENMATCHI(I)%N3=N3
          ALLOCATE(DENMATCHI(I)%MAT(N1,N2,N3))
          DENMATCHI(I)%MAT=0.D0
        ENDDO
!
!       ========================================================================
!       ==  COMPUTE DENMAT IN PARTIAL WAVE REPRESENTATION                     ==
!       ==  CAUTION: XTRANSPOSE(I)%MAT(I,J)=TRANSPOSE(INDBACK(I)%MAT(J,I)     ==
!       ==  ONLY FIRST ORDER IN THE STRUCTURE CONSTANTS                       ==
!       ==  AND OFFSITE TERMS ONLY ON THE SAME BOND IN A PRODUCT.             ==
!       ========================================================================
        DO I=1,NND
          IAT1=DENMATCHI(I)%IAT1
          IAT2=DENMATCHI(I)%IAT2
          N3=DENMATCHI(I)%N3
          IF(INDBACK(I).EQ.I) THEN  !ONSITE
            DO J=1,N3
              DENMATCHI(I)%MAT(:,:,J)=DENMATCHI(I)%MAT(:,:,J) &
      &                           +MATMUL(PIPHI(I)%MAT(:,:,1) &
      &                           ,MATMUL(DENMATPHI(I)%MAT(:,:,J) &
      &                                  ,TRANSPOSE(PIPHI(I)%MAT(:,:,1)))) 
            ENDDO
          ELSE
            DO J=1,N3
              DENMATCHI(I)%MAT(:,:,J)=DENMATCHI(I)%MAT(:,:,J) &
      &                       + MATMUL(PIPHI(I)%MAT(:,:,1) &
      &                        ,MATMUL(DENMATPHI(INDON(IAT2))%MAT(:,:,J) &
      &                            ,TRANSPOSE(PIPHI(INDON(IAT2))%MAT(:,:,1)))) &
      &                       + MATMUL(PIPHI(INDON(IAT1))%MAT(:,:,1) &
      &                        ,MATMUL(DENMATPHI(I)%MAT(:,:,J) &
      &                            ,TRANSPOSE(PIPHI(INDON(IAT2))%MAT(:,:,1)))) &
      &                       + MATMUL(PIPHI(INDON(IAT1))%MAT(:,:,1) &
      &                        ,MATMUL(DENMATPHI(INDON(IAT1))%MAT(:,:,J) &
      &                             ,TRANSPOSE(PIPHI(INDBACK(I))%MAT(:,:,1))))
!
              DENMATCHI(INDON(IAT1))%MAT(:,:,J) &
      &                                     =DENMATCHI(INDON(IAT1))%MAT(:,:,J) &
      &                       + MATMUL(PIPHI(INDON(IAT1))%MAT(:,:,1) &
      &                        ,MATMUL(DENMATPHI(I)%MAT(:,:,J) &
      &                               ,TRANSPOSE(PIPHI(I)%MAT(:,:,1)))) &
      &                       + MATMUL(PIPHI(I)%MAT(:,:,1) &
      &                        ,MATMUL(DENMATPHI(INDBACK(I))%MAT(:,:,J) &
      &                            ,TRANSPOSE(PIPHI(INDON(IAT1))%MAT(:,:,1)))) 
            ENDDO
          END IF
        ENDDO
!
       IF(TPR)CALL RSPACEOP$WRITEMAT(6,'DENMAT IN CHI-S',NND,DENMATCHI)
                                    CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_HAMILCHITOPHI(NND,PIPHI,HAMILCHI,HAMILPHI)
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
     &                             ,LOX &
     &                             ,LNX
      USE RSPACEOP_MODULE  , ONLY : RSPACEMAT_TYPE &
     &                             ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN)   :: NND
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: PIPHI(NND)
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: HAMILCHI(NND)
      TYPE(RSPACEMAT_TYPE),INTENT(OUT)  :: HAMILPHI(NND)
      INTEGER(4)          ,ALLOCATABLE  :: INDBACK(:)  !(NND)
      INTEGER(4)          ,ALLOCATABLE  :: INDON(:)    !(NAT)
      INTEGER(4)                        :: IAT1,IAT2
      INTEGER(4)                        :: ISP1,ISP2
      INTEGER(4)                        :: N1,N2,N3
      INTEGER(4)                        :: NAT
      INTEGER(4)                        :: I,J
!     **************************************************************************
                                    CALL TRACE$PUSH('SIMPLELMTO_DENMATPHITOCHI')
!
!     ==========================================================================
!     == POINTER TO BACK HOP
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(INDON(NAT))
      ALLOCATE(INDBACK(NND))
      CALL SIMPLELMTO_INDBACKON(NND,HAMILCHI,NAT,INDON,INDBACK)
!
!     ==========================================================================
!     ==  FORWARD TRANSFORMATION FROM RHOPHI TO RHOCHI                        ==
!     ==========================================================================
      DO I=1,NND
        IAT1=HAMILCHI(I)%IAT1
        IAT2=HAMILCHI(I)%IAT2
        ISP1=ISPECIES(IAT1)
        ISP2=ISPECIES(IAT2)
!
        HAMILPHI(I)%IAT1=IAT1
        HAMILPHI(I)%IAT2=IAT2
        HAMILPHI(I)%IT  =HAMILCHI(I)%IT
        N1=SUM(2*LOX(:LNX(ISP1),ISP1)+1)
        N2=SUM(2*LOX(:LNX(ISP2),ISP2)+1)
        N3=HAMILCHI(I)%N3
        HAMILPHI(I)%N1=N1
        HAMILPHI(I)%N2=N2
        HAMILPHI(I)%N3=N3
        ALLOCATE(HAMILPHI(I)%MAT(N1,N2,N3))
        HAMILPHI(I)%MAT=0.D0
      ENDDO
!
!     ========================================================================
!     ==  COMPUTE HAMILTONIAN (DENMAT) IN PARTIAL WAVE  REPRESENTATION      ==
!     ==  CAUTION: XTRANSPOSE(I)%MAT(I,J)=TRANSPOSE(INDBACK(I)%MAT(J,I)     ==
!     ==  ONLY FIRST ORDER IN THE STRUCTURE CONSTANTS                       ==
!     ==  AND OFFSITE TERMS ONLY ON THE SAME BOND IN A PRODUCT.             ==
!     ========================================================================
      DO I=1,NND
        IAT1=HAMILPHI(I)%IAT1
        IAT2=HAMILPHI(I)%IAT2
        N3  =HAMILPHI(I)%N3
        IF(INDBACK(I).EQ.I) THEN  !ONSITE
          DO J=1,N3
            HAMILPHI(I)%MAT(:,:,J)=HAMILPHI(I)%MAT(:,:,J) &
!                           __<PHI(1)|PI(1)>H(1,1)<PI(1)|PHI(1)>________________
     &                          +MATMUL(TRANSPOSE(PIPHI(I)%MAT(:,:,1)) &
     &                          ,MATMUL(HAMILCHI(I)%MAT(:,:,J) &
     &                                 ,PIPHI(I)%MAT(:,:,1)))
          ENDDO
        ELSE
!         == <PHI|PI>_I=TRANSPOSE(<PI|PHI>)_I=TRANSPOSE(<PI|PHI>)_BACK(I) ======
!         == <PHI|PI>_BACK(I)=TRANSPOSE(<PI|PHI>)_I ==== BACK(I) IAT2->IAT1 ====

          DO J=1,N3
            HAMILPHI(I)%MAT(:,:,J)= HAMILPHI(I)%MAT(:,:,J) &
!                           __<PHI(1)|PI(2)>H(2,2)<PI(2)|PHI(2)>________________
     &                    + MATMUL(TRANSPOSE(PIPHI(INDBACK(I))%MAT(:,:,1)) &
     &                     ,MATMUL(HAMILCHI(INDON(IAT2))%MAT(:,:,J) &
     &                            ,PIPHI(INDON(IAT2))%MAT(:,:,1))) &
!                           __<PHI(1)|PI(1)>H(1,1)<PI(1)|PHI(2)>________________
     &                    + MATMUL(TRANSPOSE(PIPHI(INDON(IAT1))%MAT(:,:,1)) &
     &                     ,MATMUL(HAMILCHI(INDON(IAT1))%MAT(:,:,J) &
     &                            ,PIPHI(I)%MAT(:,:,1))) &
!                           __<PHI(1)|PI(1)>H(1,2)<PI(2)|PHI(2)>________________
     &                    + MATMUL(TRANSPOSE(PIPHI(INDON(IAT1))%MAT(:,:,1)) &
     &                     ,MATMUL(HAMILCHI(I)%MAT(:,:,J) &
     &                            ,PIPHI(INDON(IAT2))%MAT(:,:,1))) 
!
            HAMILPHI(INDON(IAT1))%MAT(:,:,J) &
     &                      =HAMILPHI(INDON(IAT1))%MAT(:,:,J) &
!                           __<PHI(1)|PI(1)>H(1,2)<PI(2)|PHI(1)>________________
     &                      + MATMUL(TRANSPOSE(PIPHI(INDON(IAT1))%MAT(:,:,1)) &
     &                       ,MATMUL(HAMILCHI(I)%MAT(:,:,J) &
     &                             ,PIPHI(INDBACK(I))%MAT(:,:,1))) &
!                           __<PHI(1)|PI(2)>H(2,1)<PI(1)|PHI(1)>________________
     &                      + MATMUL(TRANSPOSE(PIPHI(INDBACK(I))%MAT(:,:,1)) &
     &                       ,MATMUL(HAMILCHI(INDBACK(I))%MAT(:,:,J) &
     &                              ,PIPHI(INDON(IAT1))%MAT(:,:,1)))
          ENDDO
        END IF
      ENDDO
                                    CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_FORCEPHITOCHI(NAT,NND,DENMATPHI,PIPHI,HAMILCHI &
     &                                   ,FORCE,RBAS,STRESS)
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
!     ** NOTE THAT PIPHI%MAT(:,:,1) IS <PI|PHI> WHILE                         **
!     **           PIPHI%MAT(:,:,2:4) IS D<PI|PHI>/DR                         **
!     **                                                                      **
!     ********************************************P. BLOECHL, GOSLAR 2020*******
      USE RSPACEOP_MODULE  , ONLY : RSPACEMAT_TYPE &
     &                             ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN)   :: NAT
      INTEGER(4)          ,INTENT(IN)   :: NND
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: DENMATPHI(NND)
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: HAMILCHI(NND)
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: PIPHI(NND)
      REAL(8)             ,INTENT(OUT)  :: FORCE(3,NAT)
      REAL(8)             ,INTENT(OUT)  :: STRESS(3,3)
      REAL(8)             ,INTENT(IN)   :: RBAS(3,3)
      INTEGER(4)          ,ALLOCATABLE  :: INDBACK(:)  !(NND)
      INTEGER(4)          ,ALLOCATABLE  :: INDON(:)    !(NAT)
      REAL(8)             ,ALLOCATABLE  :: MAT(:,:)
      REAL(8)                           :: DEDR(3)
      INTEGER(4)                        :: IAT1,IAT2
      INTEGER(4)                        :: I,J,II
!     **************************************************************************
                                    CALL TRACE$PUSH('SIMPLELMTO_FORCEPHITOCHI')
!
!     ==========================================================================
!     == POINTER TO BACK HOP
!     ==========================================================================
      ALLOCATE(INDON(NAT))
      ALLOCATE(INDBACK(NND))
      CALL SIMPLELMTO_INDBACKON(NND,DENMATPHI,NAT,INDON,INDBACK)
!
!     ==========================================================================
!     ==  CAUTION: XTRANSPOSE(I)%MAT(I,J)=TRANSPOSE(INDBACK(I)%MAT(J,I)       ==
!     ==  ONLY FIRST ORDER IN THE STRUCTURE CONSTANTS                         ==
!     ==  AND OFFSITE TERMS ONLY ON THE SAME BOND IN A PRODUCT.               ==
!     == ASSUMES THAT THE DERIVATIVE OF ONSITE <PI|PHI> WITH RESPECT TO       ==
!     == ATOMIC POSITION VANISHES.                                            ==
!     ==========================================================================
      FORCE=0.D0
      STRESS=0.D0
      DO I=1,NND
        IAT1=PIPHI(I)%IAT1
        IAT2=PIPHI(I)%IAT2
        IF(INDBACK(I).NE.I) THEN  !EXCLUDE ONSITE
          ALLOCATE(MAT(PIPHI(I)%N1,PIPHI(I)%N2))
          MAT(:,:)=0.D0
          DO J=1,DENMATPHI(I)%N3
!           == TERM WITH OFFSITE PHIPHI ========================================
            MAT=MAT+MATMUL(HAMILCHI(INDON(IAT1))%MAT(:,:,J) &
     &                    ,MATMUL(PIPHI(INDON(IAT1))%MAT(:,:,1) &
     &                           ,DENMATPHI(I)%MAT(:,:,J))) &
     &             +MATMUL(HAMILCHI(I)%MAT(:,:,J) &
     &                    ,MATMUL(PIPHI(INDON(IAT2))%MAT(:,:,1) &
     &                           ,DENMATPHI(INDON(IAT2))%MAT(:,:,J))) 
          ENDDO
          DO II=1,3
             DEDR(II)=2.D0*SUM(MAT(:,:)*PIPHI(I)%MAT(:,:,1+II))
          ENDDO   
          FORCE(:,IAT2)=FORCE(:,IAT2)-DEDR(:)
          FORCE(:,IAT1)=FORCE(:,IAT1)+DEDR(:)
!         == THIS IS JUST A SKETCH AS REMINDER TO NOT FORGET THE STRESSES ===
          DO II=1,3
            STRESS(:,II)=STRESS(:,II) &
     &               +MATMUL(RBAS,REAL(PIPHI(I)%IT,KIND=8))*DEDR(II)
          ENDDO
          DEALLOCATE(MAT)
        END IF
      ENDDO
!
                                    CALL TRACE$POP()
      RETURN
      END
MODULE SIMPLELMTO_MYMAT_MODULE
      TYPE MYMAT_TYPE
        REAL(8),ALLOCATABLE :: PROPHIH(:,:)
        REAL(8),ALLOCATABLE :: PROPHIT(:,:)
        REAL(8),ALLOCATABLE :: CMAT(:,:)
        REAL(8),ALLOCATABLE :: PHIOV(:,:)
        REAL(8),ALLOCATABLE :: PHIHHOV(:,:)
        REAL(8),ALLOCATABLE :: PHIHHOVINV(:,:)
        REAL(8),ALLOCATABLE :: PHIHTOV(:,:)
      END TYPE MYMAT_TYPE
END MODULE SIMPLELMTO_MYMAT_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_MAKEPIPHI(NND,SBARE,PIPHI)
!     **************************************************************************
!     **  CALCULATE THE MATRIX <PI|PHI>WHICH CONVERTS PARTIAL-WAVE PROJECTIONS**
!     **  <PTILDE|PSITILDE> INTO LOCAL-ORBITAL PROJECTIONS <PI|PSI>.          **
!     **                                                                      **
!     **            |PSI> APPROX |CHI><PI|PHI><PTILDE|PSITILDE>               **
!     **                                                                      **
!     **  1: OFF-SITE STRUCTURE CONSTANTS ENTER ONLY TO FIRST ORDER.          **
!     **     AS A RESULT PHICHI AND SBARE ARE ON THE SAME NEIGHBORLIST        **
!     **                                                                      **
!     ** LNX,LOX REFER TO THE PARTIAL-WAVE REPRESENTATION                     **
!     **                                                                      **
!     ********************************************P. BLOECHL, GOSLAR 2019*******
      USE SIMPLELMTO_MODULE, ONLY : ISPECIES &
     &                             ,POTPAR=>POTPAR1 &
     &                             ,LOX &
     &                             ,LNX &
     &                             ,NSP
      USE RSPACEOP_MODULE  , ONLY : RSPACEMAT_TYPE &
     &                             ,RSPACEOP$WRITEMAT
      USE SIMPLELMTO_MYMAT_MODULE, ONLY : MYMAT_TYPE
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN)   :: NND
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: SBARE(NND)
      TYPE(RSPACEMAT_TYPE),INTENT(OUT)  :: PIPHI(NND)
      LOGICAL(4)          ,PARAMETER    :: TPR=.FALSE.
      TYPE(MYMAT_TYPE)    ,ALLOCATABLE  :: MYMAT(:)
      INTEGER(4)          ,ALLOCATABLE  :: INDBACK(:)  !(NND)
      INTEGER(4)          ,ALLOCATABLE  :: INDON(:)    !(NAT)
      INTEGER(4)                        :: IAT1,IAT2
      INTEGER(4)                        :: ISP,ISP1,ISP2
      INTEGER(4)                        :: N1,N2,N3
      INTEGER(4)                        :: NAT
      INTEGER(4)                        :: I,II
      REAL(8)                           :: S
!     **************************************************************************
                                    CALL TRACE$PUSH('SIMPLELMTO_MAKEPIPHI')
      NAT=SIZE(ISPECIES)
      ALLOCATE(MYMAT(NSP))
      CALL SIMPLELMTO_EXPANDTOMYMAT(NSP,MYMAT)
!
!     ==========================================================================
!     == POINTER TO BACK HOP
!     ==========================================================================
      ALLOCATE(INDON(NAT))
      ALLOCATE(INDBACK(NND))
      CALL SIMPLELMTO_INDBACKON(NND,SBARE,NAT,INDON,INDBACK)
!
!     ==========================================================================
!     ==  <PI|PHI>            |PSI>=|CHI><PI|PHI><P|PSI>                      ==
!     ==========================================================================
      DO I=1,NND
        IAT1=SBARE(I)%IAT1
        IAT2=SBARE(I)%IAT2
        ISP1=ISPECIES(IAT1)
        ISP2=ISPECIES(IAT2)
        N1=SUM(2*POTPAR(ISP1)%LOXH+1)
        N2=SUM(2*LOX(:LNX(ISP2),ISP2)+1)
        N3=4    !VALUE AND GRADIENT
        PIPHI(I)%IAT1=IAT1
        PIPHI(I)%IAT2=IAT2
        PIPHI(I)%IT  =SBARE(I)%IT
        PIPHI(I)%N1=N1
        PIPHI(I)%N2=N2
        PIPHI(I)%N3=N3
        ALLOCATE(PIPHI(I)%MAT(N1,N2,N3))
        PIPHI(I)%MAT=0.D0
!
        IF(INDBACK(I).EQ.I) THEN
          PIPHI(I)%MAT(:,:,1)=TRANSPOSE(MYMAT(ISP1)%PROPHIH)
        ELSE
          DO II=1,4
!           == S=-1 FOR GRADIENTS TAKING CARE OF INDBACK, ELSE S=+1 ============
            S=REAL(SIGN(1,3-2*II),KIND=8) 
            PIPHI(I)%MAT(:,:,II)=-MATMUL(TRANSPOSE(MYMAT(ISP1)%CMAT) &
     &                        ,MATMUL(SBARE(I)%MAT(:,:,II) &
     &                        ,TRANSPOSE(MYMAT(ISP2)%PROPHIT &
     &                                  -MATMUL(MYMAT(ISP2)%PROPHIH &
     &                                         ,MATMUL(MYMAT(ISP2)%PHIHHOVINV &
     &                                              ,MYMAT(ISP2)%PHIHTOV))))) &
     &                       +S*MATMUL(MYMAT(ISP1)%PHIHTOV &
     &                        ,MATMUL(TRANSPOSE(SBARE(INDBACK(I))%MAT(:,:,II)) &
     &                        ,MATMUL(MYMAT(ISP2)%CMAT &
     &                        ,MATMUL(MYMAT(ISP2)%PHIHHOVINV &
     &                               ,TRANSPOSE(MYMAT(ISP2)%PROPHIH)))))
          ENDDO
        END IF
        DO II=1,4
          PIPHI(I)%MAT(:,:,II)=MATMUL(MYMAT(ISP1)%PHIHHOVINV &
                             ,MATMUL(PIPHI(I)%MAT(:,:,II),MYMAT(ISP2)%PHIOV))
        ENDDO
      ENDDO
      IF(TPR)CALL RSPACEOP$WRITEMAT(6,'<PI|PHI>',NND,PIPHI)
!
!     ==========================================================================
!     ==  CLEAN UP                                                            ==
!     ==========================================================================
      DO ISP=1,NSP
        DEALLOCATE(MYMAT(ISP)%PROPHIH)
        DEALLOCATE(MYMAT(ISP)%PROPHIT)
        DEALLOCATE(MYMAT(ISP)%PHIOV)
        DEALLOCATE(MYMAT(ISP)%CMAT)
        DEALLOCATE(MYMAT(ISP)%PHIHHOV)
        DEALLOCATE(MYMAT(ISP)%PHIHTOV)
      ENDDO
      DEALLOCATE(MYMAT)
                                    CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_EXPANDTOMYMAT(NSP,MYMAT)
!     **************************************************************************
!     ** LNX,LOX REFER TO THE PARTIAL-WAVE REPRESENTATION                     **
!     **                                                                      **
!     ********************************************P. BLOECHL, GOSLAR 2019*******
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1 &
     &                             ,LOX &
     &                             ,LNX 
      USE SIMPLELMTO_MYMAT_MODULE, ONLY : MYMAT_TYPE
      IMPLICIT NONE
      INTEGER(4)      ,INTENT(IN) :: NSP
      TYPE(MYMAT_TYPE),INTENT(OUT):: MYMAT(NSP)
      LOGICAL(4)      ,PARAMETER  :: TPR=.FALSE.
      INTEGER(4)                  :: ISP
      INTEGER(4)                  :: LMNX1,LMNX2
      INTEGER(4)                  :: LMN1I,LMN2I
      INTEGER(4)                  :: LN1,LN2
      INTEGER(4)                  :: L1,L2
      INTEGER(4)                  :: IM
      INTEGER(4)                  :: I
!     **************************************************************************
                                    CALL TRACE$PUSH('SIMPLELMTO_EXPANDTOMYMAT')
!
!     ==========================================================================
!     == EXPAND <P|PHI-H>                                                     ==
!     ==========================================================================
      DO ISP=1,NSP
        LMNX1=SUM(2*LOX(:LNX(ISP),ISP)+1)
        LMNX2=SUM(2*POTPAR(ISP)%LOXH+1)
        ALLOCATE(MYMAT(ISP)%PROPHIH(LMNX1,LMNX2))
        MYMAT(ISP)%PROPHIH=0.D0

        LMN1I=0
        DO LN1=1,LNX(ISP)
          L1=LOX(LN1,ISP)
!
          LMN2I=0 
          DO LN2=1,POTPAR(ISP)%LNXH
            L2=POTPAR(ISP)%LOXH(LN2)
            IF(L2.EQ.L1) THEN
              DO IM=1,2*L2+1
                MYMAT(ISP)%PROPHIH(LMN1I+IM,LMN2I+IM) &
    &                                              =POTPAR(ISP)%PROPHIH(LN1,LN2)
              ENDDO
            END IF
            LMN2I=LMN2I+2*L2+1
          ENDDO
          LMN1I=LMN1I+2*L1+1     
        ENDDO
!
        IF(TPR) THEN
          WRITE(*,FMT='(80("="),T20," <PRO|PHI-HEAD> FOR ISP=",I3," ")')ISP
          DO I=1,LMNX1
            WRITE(*,FMT='(10F10.3)')MYMAT(ISP)%PROPHIH(I,:)
          ENDDO
        END IF
!
      ENDDO
!
!     ==========================================================================
!     == EXPAND <P|PHI-T>                                                     ==
!     ==========================================================================
      DO ISP=1,NSP
        LMNX1=SUM(2*LOX(:LNX(ISP),ISP)+1)
        LMNX2=SUM(2*POTPAR(ISP)%LOXT+1)
        ALLOCATE(MYMAT(ISP)%PROPHIT(LMNX1,LMNX2))
        MYMAT(ISP)%PROPHIT=0.D0

        LMN1I=0
        DO LN1=1,LNX(ISP)
          L1=LOX(LN1,ISP)
!
          LMN2I=0 
          DO LN2=1,POTPAR(ISP)%LNXT
            L2=POTPAR(ISP)%LOXT(LN2)
            IF(L2.EQ.L1) THEN
              DO IM=1,2*L2+1
                MYMAT(ISP)%PROPHIT(LMN1I+IM,LMN2I+IM) &
    &                                              =POTPAR(ISP)%PROPHIT(LN1,LN2)
              ENDDO
            END IF
            LMN2I=LMN2I+2*L2+1
          ENDDO
          LMN1I=LMN1I+2*L1+1     
        ENDDO
!
        IF(TPR) THEN
          WRITE(*,FMT='(80("="),T20," <PRO|PHI-TAIL> FOR ISP=",I3," ")')ISP
          DO I=1,LMNX1
            WRITE(*,FMT='(10F10.3)')MYMAT(ISP)%PROPHIT(I,:)
          ENDDO
        END IF
!
      ENDDO
!
!     ==========================================================================
!     == EXPAND PARTIAL WAVE OVERLAP                                          ==
!     ==========================================================================
      DO ISP=1,NSP
        LMNX1=SUM(2*LOX(:LNX(ISP),ISP)+1)
        LMNX2=LMNX1
        ALLOCATE(MYMAT(ISP)%PHIOV(LMNX1,LMNX2))
        MYMAT(ISP)%PHIOV=0.D0

        LMN1I=0
        DO LN1=1,LNX(ISP)
          L1=LOX(LN1,ISP)
!
          LMN2I=0 
          DO LN2=1,LNX(ISP)
            L2=LOX(LN2,ISP)
            IF(L2.EQ.L1) THEN
              DO IM=1,2*L2+1
                MYMAT(ISP)%PHIOV(LMN1I+IM,LMN2I+IM)=POTPAR(ISP)%PHIOV(LN1,LN2)
              ENDDO
            END IF
            LMN2I=LMN2I+2*L2+1
          ENDDO
          LMN1I=LMN1I+2*L1+1     
        ENDDO
!
        IF(TPR) THEN
          WRITE(*,FMT='(80("="),T20," <PHI|THETA|PHI> FOR ISP=",I3," ")')ISP
          DO I=1,LMNX1
            WRITE(*,FMT='(10F10.3)')MYMAT(ISP)%PHIOV(I,:)
          ENDDO
        END IF
!
      ENDDO
!
!     ==========================================================================
!     == CMAT (RESPONSIBLE FOR NORMALIZATION OF CHI)                          ==
!     ==========================================================================
      DO ISP=1,NSP
        LMNX1=SUM(2*POTPAR(ISP)%LOXH+1)
        LMNX2=LMNX1
        ALLOCATE(MYMAT(ISP)%CMAT(LMNX1,LMNX2))
        MYMAT(ISP)%CMAT=0.D0
        LMN1I=0
        DO LN1=1,POTPAR(ISP)%LNXH
          L1=POTPAR(ISP)%LOXH(LN1)
!
          LMN2I=0 
          DO LN2=1,POTPAR(ISP)%LNXH
            L2=POTPAR(ISP)%LOXH(LN2)
            IF(L2.EQ.L1) THEN
              DO IM=1,2*L2+1
                 MYMAT(ISP)%CMAT(LMN1I+IM,LMN2I+IM)=POTPAR(ISP)%CMAT(LN1,LN2)
              ENDDO
            END IF
            LMN2I=LMN2I+2*L2+1
          ENDDO
          LMN1I=LMN1I+2*L1+1     
        ENDDO
!
        IF(TPR) THEN
          WRITE(*,FMT='(80("="),T20," CMAT FOR ISP=",I3," ")')ISP
          DO I=1,LMNX1
            WRITE(*,FMT='(10F10.3)')MYMAT(ISP)%CMAT(I,:)
          ENDDO
        END IF
!
      ENDDO
!
!     ==========================================================================
!     == <PHI-H|THETA|PHI-H> AND <PHI-H|THETA|PHI_J> 
!     ==========================================================================
      DO ISP=1,NSP
        LMNX1=SUM(2*POTPAR(ISP)%LOXH+1)
        ALLOCATE(MYMAT(ISP)%PHIHHOV(LMNX1,LMNX1))
        MYMAT(ISP)%PHIHHOV=MATMUL(TRANSPOSE(MYMAT(ISP)%PROPHIH) &
     &                           ,MATMUL(MYMAT(ISP)%PHIOV,MYMAT(ISP)%PROPHIH))
!
        ALLOCATE(MYMAT(ISP)%PHIHHOVINV(LMNX1,LMNX1))
        PRINT*,'MARKE 1',ISP,LMNX1 &
       &               ,' LNXH=',POTPAR(ISP)%LNXH &
       &               ,' LOXH=',POTPAR(ISP)%LOXH
        CALL LIB$INVERTR8(LMNX1,MYMAT(ISP)%PHIHHOV,MYMAT(ISP)%PHIHHOVINV)
PRINT*,'MARKE 2'
!
        LMNX1=SUM(2*POTPAR(ISP)%LOXH+1)
        LMNX2=SUM(2*POTPAR(ISP)%LOXT+1)
        ALLOCATE(MYMAT(ISP)%PHIHTOV(LMNX1,LMNX2))
        MYMAT(ISP)%PHIHTOV=MATMUL(TRANSPOSE(MYMAT(ISP)%PROPHIH) &
     &                           ,MATMUL(MYMAT(ISP)%PHIOV,MYMAT(ISP)%PROPHIT))
!
        IF(TPR) THEN
          WRITE(*,FMT='(80("="),T20," <PHI-H|THETA|PHI-H> FOR ISP=",I3," ")')ISP
          DO I=1,LMNX1
            WRITE(*,FMT='(10F10.3)')MYMAT(ISP)%PHIHHOV(I,:)
          ENDDO
          WRITE(*,FMT='(80("="),T20," <PHI-H|THETA|PHI-T> FOR ISP=",I3," ")')ISP
          DO I=1,LMNX1
            WRITE(*,FMT='(10F10.3)')MYMAT(ISP)%PHIHTOV(I,:)
          ENDDO
        END IF
!
      ENDDO
                                    CALL TRACE$POP()
      RETURN
      END
!     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_STRUCTURECONSTANTS(NAT,R0,RBAS,NNS,SBARE)
!     **************************************************************************
!     **  CONSTRUCTS THE STRUCTURE CONSTANTS THAT MEDIATE AN EXPANSION        **
!     **  OF A SOLID HANKEL FUNCTION H_{L,M}(R-R1) CENTERED AT R1             **
!     **  INTO SOLID BESSEL FUNCTIONS  J_{L,M}(R-R2) CENTERED AT R2           **
!     **                                                                      **
!     **    H_{L,M}(R-R1) = - SUM_{L',M'} S_{L,M,L',M'} * J_{L',M'}(R-R2)     **
!     **                                                                      **
!     **  STRUCTURE CONSTANTS AND FORCES ARE CALCULATED                       **
!     **                                                                      **
!     **  SBARE HAS IAT1,IAT2,IT ALREADY SET, SBARE%MAT IS NOT ALLOCATED      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : ISPECIES &
     &                             ,POTPAR=>POTPAR1 &
     &                             ,K2
      USE RSPACEOP_MODULE, ONLY : RSPACEMAT_TYPE &
     &                           ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN)   :: NAT        ! #(ATOMS)
      REAL(8)             ,INTENT(IN)   :: R0(3,NAT)  ! ATOMIC POSITIONS
      REAL(8)             ,INTENT(IN)   :: RBAS(3,3)  ! LATTICE VECTORS
      INTEGER(4)          ,INTENT(IN)   :: NNS        ! #(NEIGHBORLIST ENTRIES)
!     == ON INPUT, SBARE CONTAINS IAT1,IAT2,IT =================================
!     == SBARE HAS N3=4 COMPONENTS, NAMELY VALUE AND THREE DERIVATIVES WITH DR==
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: SBARE(NNS) ! STRUCTURE CONSTANTS
      INTEGER(4)                        :: IAT1,IAT2
      INTEGER(4)                        :: ISP1,ISP2
      INTEGER(4)                        :: N1,N2,N3
      REAL(8)                           :: R21(3) 
      INTEGER(4)                        :: L1X,L2X
      REAL(8)             ,ALLOCATABLE  :: SMAT(:,:,:) !VALUE AND GRADIENT
      INTEGER(4)                        :: I
      INTEGER(4)                        :: LMN1A,LMN1B,LMN2A,LMN2B
      INTEGER(4)                        :: LM1A,LM1B,LM2A,LM2B
      INTEGER(4)                        :: L1,L2
      INTEGER(4)                        :: LN1,LN2
!     **************************************************************************
!
!     ==========================================================================
!     ==  CALCULATE BARE STRUCTURE CONSTANTS                                  ==
!     ==========================================================================
      DO I=1,NNS
        IAT1=SBARE(I)%IAT1
        IAT2=SBARE(I)%IAT2
        ISP1=ISPECIES(IAT1)
        ISP2=ISPECIES(IAT2)
        N1=SUM(2*POTPAR(ISP1)%LOXH+1)
        N2=SUM(2*POTPAR(ISP2)%LOXT+1)
        N3=4   ! VALUE + THREE GRADIENT COMPONENTS
        SBARE(I)%N1=N1
        SBARE(I)%N2=N2
        SBARE(I)%N3=N3
        ALLOCATE(SBARE(I)%MAT(N1,N2,N3)) ! VALUE AND GRADIENT
        SBARE(I)%MAT=0.D0
!       == SKIP ONSITE TERMS. THEY ARE ZERO ====================================
        IF(SBARE(I)%IAT1.EQ.SBARE(I)%IAT2.AND.SUM(SBARE(I)%IT**2).EQ.0) CYCLE
!
!       ========================================================================
!       == CALCULATE BARE STRUCTURE CONSTANTS                                 ==
!       ========================================================================
        R21=R0(:,IAT2)+MATMUL(RBAS,REAL(SBARE(I)%IT,KIND=8))-R0(:,IAT1)
        L1X=MAXVAL(POTPAR(ISP1)%LOXH)
        L2X=MAXVAL(POTPAR(ISP2)%LOXT)
        ALLOCATE(SMAT((L1X+1)**2,(L2X+1)**2,4)) ! VALUE AND GRADIENT
        CALL LMTO$STRUCTURECONSTANTSWGRAD(R21,K2,L1X,L2X,SMAT(:,:,1) &
       &                                                ,SMAT(:,:,2:))
!
!       ========================================================================
!       == MAP STRUCTURE CONSTANTS ON SBARE                                   ==
!       ==   K_I(R) = - S_{I,J}J_J(R)                                         ==
!       ==      <K| = - S * <J|                                               ==
!       == |KINFTY> = |KOMEGA> - |JOMEGA> TRANSPOSE(S)                        ==
!       ========================================================================
        LMN1A=1
        DO LN1=1,POTPAR(ISP1)%LNXH
          L1=POTPAR(ISP1)%LOXH(LN1)
          LMN1B=LMN1A+2*L1
          LM1A=L1**2+1
          LM1B=(L1+1)**2
          LMN2A=1
          DO LN2=1,POTPAR(ISP2)%LNXT
            L2=POTPAR(ISP2)%LOXT(LN2)
            LMN2B=LMN2A+2*L2
            LM2A=L2**2+1
            LM2B=(L2+1)**2
            SBARE(I)%MAT(LMN1A:LMN1B,LMN2A:LMN2B,:)=SMAT(LM1A:LM1B,LM2A:LM2B,:)
            LMN2A=LMN2B+1
          ENDDO
          LMN1A=LMN1B+1
        ENDDO
        DEALLOCATE(SMAT)
      ENDDO
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_HYBRIDENERGY(NAT,NND,RBAS,R0,DENMAT,DONSITE &
     &                                  ,ETOT,STRESS,FORCE,HAMIL,HONSITE &
     &                                  ,OVERLAPCHI)
!     **************************************************************************
!     **  WORK OUT THE ENERGY USING THE LOCAL APPROXIMATION                   **
!     **  TAILED PARTIAL WAVES                                                **
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY: ISPECIES &
     &                            ,LNXPHI=>LNX &
     &                            ,LOXPHI=>LOX &
     &                            ,POTPAR=>POTPAR1 &
     &                            ,TOFFSITE &
     &                            ,HYBRIDSETTING
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
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: HAMIL(NND)      !INTENT(OUT)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: OVERLAPCHI(NND) !INTENT(OUT)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: HONSITE(NAT)    !INTENT(OUT)
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
      LOGICAL(4)            :: TACTIVE ! DOES THIS ATOM CONTRIBUTE?
      LOGICAL(4),PARAMETER  :: TPARALLEL=.FALSE.
      INTEGER(4)            :: THISTASK,NTASKS
LOGICAL(4)            :: TBACK
REAL(8)               :: RBAS1(3,3),EX1,EX2,DEL,RBASINV(3,3),STRESS1(3,3)
REAL(8)               :: STRESS0(3,3)
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
!     == OF THIS ROUTINE, THE RESULT IS SUMMED OVER ALL TASKS.                ==
!     ==                                                                      ==
!     == THE PARAMETER TPARALLEL DECIDES WHETHER THE PARALLELIZATION IS DONE  ==
!     == OR WETHER EVERYTHING IS CALCULATED ON EACH TASK                      ==
!     ==                                                                      ==
!     ==========================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      ETOT=0.D0
      DO I=1,NND
        HAMIL(I)%MAT=0.D0
      ENDDO
      DO IAT=1,NAT
        HONSITE(IAT)%MAT=0.D0
      ENDDO
      FORCE=0.D0
      STRESS=0.D0
!
!     ==========================================================================
!     == ONSITE EXCHANGE CORRECTION                                           ==
!     ==========================================================================
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        TACTIVE=(POTPAR(ISP)%LNXH.GT.0)
        IF(.NOT.TACTIVE) CYCLE
!
        IF(TPARALLEL.AND.MOD(IAT-1,NTASKS).NE.THISTASK-1) CYCLE
PRINT*,'++THISTASK ',TPARALLEL,THISTASK,NTASKS,IAT,MOD(IAT-1,NTASKS)
!
!       ========================================================================
!       ==                                                                    ==
!       ========================================================================
        LHFWEIGHT=HYBRIDSETTING(ISP)%LHFWEIGHT
!
!       ========================================================================
!       == FIND INDEX TO THE ONSITE ELEMENTS OF THE DENSITY MATRIX
!       ========================================================================
        CALL SIMPLELMTO_INDEXLOCAL2(IAT,NND,DENMAT,IND)
        LMNX=DENMAT(IND)%N1
        NDIMD=DENMAT(IND)%N3
        LMNXPHI=DONSITE(IAT)%N1
!
!       ========================================================================
!       == WORK OUT TOTAL ENERGY AND HAMILTONIAN                              ==
!       ========================================================================
        CALL REPORT$TITLE(6,'ONSITE PBE0R CORRECTIONS')
        CALL REPORT$I4VAL(6,'ATOM',IAT,' ')
        ALLOCATE(HON(LMNXPHI,LMNXPHI,NDIMD))
        ALLOCATE(H(LMNX,LMNX,NDIMD))
        CALL SIMPLELMTO_ONSITEX(ISP,HYBRIDSETTING(ISP)%TCV,NDIMD,LMNX,LMNXPHI &
      &                 ,POTPAR(ISP)%ONSITEU,DENMAT(IND)%MAT,DONSITE(IAT)%MAT &
      &                 ,EX,H,HON)
        HONSITE(IAT)%MAT=HONSITE(IAT)%MAT+HON*LHFWEIGHT
        HAMIL(IND)%MAT  =HAMIL(IND)%MAT  +H  *LHFWEIGHT
        ETOT            =ETOT            +EX *LHFWEIGHT  
                                         !HFWEIGHT IS MULTIPLIED ON LATER
        DEALLOCATE(HON)
        DEALLOCATE(H)
!
!       == INCLUDE ONSITE TERM OF THE OVERLAP MATRIX ===========================
        OVERLAPCHI(IND)%MAT(:,:,1)=POTPAR(ISP)%OVERLAP(:,:)
!
!       == CALCULATE CHARGE AND SPIN FOR DIAGONISTIC PURPOSES ==================
        QSPIN=0.D0
        DO I=1,LMNX
          DO J=1,LMNX
            QSPIN(:NDIMD)=QSPIN(:NDIMD) &
       &                 +POTPAR(ISP)%OVERLAP(I,J)*DENMAT(IND)%MAT(J,I,:)
          ENDDO
        ENDDO
        CALL REPORT$R8VAL(6,'CHARGE',QSPIN(1),'(Q_E)')
      ENDDO !END OF LOOP OVER ATOMS
!
!     ==========================================================================
!     == OFFSITE EXCHANGE CONTRIBUTION                                        ==
!     ==========================================================================
      IF(TOFFSITE) THEN
!       == USES NUMERICAL MATRIX ELEMENTS ======================================
!       == UNLIKE THE TOTAL ENERGY, FORCES AND STRESSES ARE NOT ADDED, BUT    ==
!       == DIRECTLY WRITTEN IN FINAL ARRAYS ====================================
        CALL SIMPLELMTO_OFFSITEXEVAL(TPARALLEL,NAT,NND,RBAS,R0,DENMAT &
     &                              ,EX,STRESS,FORCE,HAMIL,OVERLAPCHI) 
        PRINT*,'+-+-+-+  OFFSITE EX=',EX
        ETOT=ETOT+EX

!!$!TESTBEGINTESTBEGINTESTBEGINTESTBEGINTESTBEGINTESTBEGINTESTBEGINTESTBEGINTEST
!!$!****************************************************************************
!!$!** THIS TESTS THE GRADIENTS OF SIMPLELMTO_OFFSITEXEVAL BY NUMERICAL 
!!$!** DIFFERENTIATION. CAUTION! IT TESTS DEDRBAS AND NOT STRESS. CHANGE 
!!$!** SIMPLELMTO_OFFSITEXEVAL SO THAT IT OUTPUTS DEDRBAS INSTEAD OF
!!$!** STRESS=RBAS*DEDRBAS
!!$!****************************************************************************
!!$ CALL SIMPLELMTO_OFFSITEXEVAL(TPARALLEL,NAT,NND,RBAS,R0,DENMAT &
!!$&                              ,EX2,STRESS,FORCE,HAMIL,OVERLAPCHI) 
!!$DO I=1,3
!!$DO J=1,3
!!$ DEL=1.D-3
!!$ RBAS1=RBAS
!!$ RBAS1(I,J)=RBAS(I,J)+DEL
!!$
!!$ CALL SIMPLELMTO_OFFSITEXEVAL(TPARALLEL,NAT,NND,RBAS1,R0,DENMAT &
!!$&                              ,EX2,STRESS0,FORCE,HAMIL,OVERLAPCHI) 
!!$ RBAS1=RBAS
!!$ RBAS1(I,J)=RBAS(I,J)-DEL
!!$ CALL SIMPLELMTO_OFFSITEXEVAL(TPARALLEL,NAT,NND,RBAS1,R0,DENMAT &
!!$&                              ,EX1,STRESS0,FORCE,HAMIL,OVERLAPCHI) 
!!$ STRESS1(I,J)=(EX2-EX1)/(2.D0*DEL)
!!$ WRITE(*,FMT='("STRESSTEST ",2I2,3F10.5)') &
!!$         I,J,STRESS(I,J),STRESS1(I,J),RBAS(I,J)
!!$ENDDO
!!$ENDDO
!!$WRITE(*,FMT='(A10,3(3F10.5," "))')'STRESS',TRANSPOSE(STRESS)
!!$WRITE(*,FMT='(A10,3(3F10.5," "))')'STRESS1',STRESS1
!!$CALL ERROR$STOP('FORCED STOP TEST OF SIMPLELMTO_OFFSITEXEVAL')
!!$!TESTENDTESTENDTESTENDTESTENDTESTENDTESTENDTESTENDTESTENDTESTENDTESTENDTESTEN
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
!     == MAKE HAMILTONIAN HERMITIAN (ONLY FOR TESTING; NOT REQUIRED)          ==
!     ==========================================================================
!     CALL SIMPLELMTO_SYMMETRIZEHAMIL(NAT,NND,HAMIL,HONSITE,OVERLAPCHI)
!
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_SYMMETRIZEHAMIL(NAT,NND,HAMIL,HONSITE,OVERL)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE RSPACEOP_MODULE, ONLY : RSPACEMAT_TYPE &
     &                           ,RSPACEOP$WRITEMAT
      IMPLICIT NONE
      INTEGER(4)          ,INTENT(IN)   :: NND   ! #(NEIGHBORLIST ENTRIES)
      INTEGER(4)          ,INTENT(IN)   :: NAT   ! #(ATOMS)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: HAMIL(NND)      !
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: OVERL(NND)    !
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: HONSITE(NAT)    !
      INTEGER(4)                        :: NN1,NN2,IAT,IAT1,IAT2,IT(3),I
!     **************************************************************************
                             CALL TRACE$PUSH('SIMPLELMTO_SYMMETRIZEHAMIL')
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
            OVERL(NN1)%MAT(:,:,I)=0.5D0 &
     &                 *(OVERL(NN1)%MAT(:,:,I)+TRANSPOSE(OVERL(NN2)%MAT(:,:,I)))
            OVERL(NN2)%MAT(:,:,I)=TRANSPOSE(OVERL(NN1)%MAT(:,:,I))
          ENDDO
        ENDDO
      ENDDO
      DO IAT=1,NAT
        DO I=1,HONSITE(IAT)%N3
          HONSITE(IAT)%MAT(:,:,I)=0.5D0 &
     &             *(HONSITE(IAT)%MAT(:,:,I)+TRANSPOSE(HONSITE(IAT)%MAT(:,:,I)))
        ENDDO
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END SUBROUTINE SIMPLELMTO_SYMMETRIZEHAMIL

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_ONSITEX(ISP,TCV,NDIMD,LMNXCHI,LMNXPHI &
      &                            ,UCHI,DCHI,DPHI,EX,HCHI,HPHI)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP
      LOGICAL(4),INTENT(IN) :: TCV
      INTEGER(4),INTENT(IN) :: NDIMD    ! #(SPIN COMPONENTS (1,2, OR 4)
      INTEGER(4),INTENT(IN) :: LMNXCHI  ! #(LOCAL ORBITALS)
      INTEGER(4),INTENT(IN) :: LMNXPHI  ! #(PARTIAL WAVES)
      REAL(8)   ,INTENT(IN) :: UCHI(LMNXCHI,LMNXCHI,LMNXCHI,LMNXCHI)
      REAL(8)   ,INTENT(IN) :: DCHI(LMNXCHI,LMNXCHI,NDIMD)
      REAL(8)   ,INTENT(IN) :: DPHI(LMNXPHI,LMNXPHI,NDIMD)
      REAL(8)   ,INTENT(OUT):: EX
      REAL(8)   ,INTENT(OUT):: HCHI(LMNXCHI,LMNXCHI,NDIMD)
      REAL(8)   ,INTENT(OUT):: HPHI(LMNXPHI,LMNXPHI,NDIMD)
      LOGICAL(4),PARAMETER  :: TPR=.TRUE.
      REAL(8)   ,ALLOCATABLE:: HPHI1(:,:,:) !(LMNXPHI,LMNXPHI,NDIMD) 
      REAL(8)   ,ALLOCATABLE:: HCHI1(:,:,:) !(LMNXCHI,LMNXCHI,NDIMD) 
      INTEGER(4)            :: I,J,K,L
      REAL(8)               :: EX1
      REAL(8)               :: SVAR
!     **************************************************************************
      EX=0.D0
      HCHI(:,:,:)=0.D0
      HPHI(:,:,:)=0.D0
      ALLOCATE(HCHI1(LMNXCHI,LMNXCHI,NDIMD))
      ALLOCATE(HPHI1(LMNXPHI,LMNXPHI,NDIMD))
!
!     ==========================================================================
!     == EXACT EXCHANGE ENERGY                                                ==
!     ==========================================================================
      DO I=1,LMNXCHI
        DO J=1,LMNXCHI
          DO K=1,LMNXCHI
            DO L=1,LMNXCHI
!             ==================================================================
!             == HARTREE TERM (NOT CONSIDERED)                                ==
!             == EH=EH+0.5D0*U(I,J,K,L)*D(K,I,1)*D(L,J,1)                     ==
!             == H(K,I,1)=H(K,I,1)+0.5D0*U(I,J,K,L)*D(L,J,1) !DE/DRHO(K,I,1)  ==
!             == H(L,J,1)=H(L,J,1)+0.5D0*U(I,J,K,L)*D(K,I,1) !DE/DRHO(L,J,1)  ==
!             ==================================================================
!             == AN ADDITIONAL FACTOR COMES FROM THE REPRESENTATION ============
!             == INTO TOTAL AND SPIN ===========================================
              SVAR=-0.25D0*UCHI(I,J,K,L)
              EX=EX+SVAR*SUM(DCHI(K,J,:)*DCHI(L,I,:))
              HCHI(K,J,:)=HCHI(K,J,:)+SVAR*DCHI(L,I,:) 
              HCHI(L,I,:)=HCHI(L,I,:)+SVAR*DCHI(K,J,:) 
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF(TPR)CALL REPORT$R8VAL(6,'EXACT LOCAL EXCHANGE',EX,'H')
!
!     ========================================================================
!     == ADD CORE-VALENCE EXCHANGE                                          ==
!     ========================================================================
      IF(TCV) THEN
        CALL SIMPLELMTO_CVX_ACTONPHI(ISP,LMNXPHI,NDIMD,DPHI,EX1,HPHI1)
        EX=EX+EX1
        HPHI=HPHI+HPHI1
        IF(TPR)CALL REPORT$R8VAL(6,'CORE-VALENCE EXCHANGE',EX1,'H')
      END IF
!
!     ========================================================================
!     == DOUBLE COUNTING CORRECTION (EXCHANGE ONLY)                         ==
!     == THIS IS THE TIME CONSUMING PART                                    ==
!     ========================================================================
      CALL TIMING$CLOCKON('ENERGYTEST:DC')      
      CALL DFT$SETL4('XCONLY',.TRUE.)
      CALL SIMPLELMTO_DC(ISP,NDIMD,LMNXCHI,DCHI,LMNXPHI,DPHI &
     &                      ,EX1,HCHI1,HPHI1)
      CALL DFT$SETL4('XCONLY',.FALSE.)
      HPHI=HPHI-HPHI1
      HCHI=HCHI-HCHI1
      EX  =EX  -EX1
      IF(TPR)CALL REPORT$R8VAL(6,'DOUBLE COUNTING CORRECTION',-EX1,'H')
      IF(TPR)CALL REPORT$R8VAL(6,'ONSITE EXCHANGE',EX,'H')
      CALL TIMING$CLOCKOFF('ENERGYTEST:DC')      
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
      SUBROUTINE SIMPLELMTO_CVX_ACTONPHI(ISP,LMNX,NDIMD,DENMAT,ETOT,DH)
!     **************************************************************************
!     **  CORE VALENCE EXCHANGE ENERGY ACTING ON PARTIAL WAVES                **
!     *****************************PETER BLOECHL, GOSLAR 2011-2019**************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: ISP          ! ATOM TYPE INDEX
      INTEGER(4),INTENT(IN)  :: LMNX
      INTEGER(4),INTENT(IN)  :: NDIMD
      REAL(8)   ,INTENT(OUT) :: ETOT              ! CORE-VALENCE ENERGY
      REAL(8)   ,INTENT(IN)  :: DENMAT(LMNX,LMNX,NDIMD)
      REAL(8)   ,INTENT(OUT) :: DH(LMNX,LMNX,NDIMD)
      LOGICAL(4),PARAMETER   :: TPR=.FALSE.
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
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      ALLOCATE(CVXMAT(LNX,LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      CALL SETUP$GETR8A('CVX',LNX*LNX,CVXMAT)
      CALL SETUP$UNSELECT()
      IF(LMNX.NE.SUM(2*LOX+1)) THEN
        CALL ERROR$MSG('INCONSISTENT INPUT VALUE LMNX')
        CALL ERROR$I4VAL('ISP',ISP)
        CALL ERROR$I4VAL('LNX',LNX)
        CALL ERROR$I4VAL('LMNX',LMNX)
        CALL ERROR$I4VAL('SUM(2*LOX+1)',SUM(2*LOX+1))
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
!     == DIAGNOSTIC PRINTOUT                                                  ==
!     ==========================================================================
      IF(TPR) THEN
        PRINT*,' LOX=',LOX
        PRINT*,' E_CVX=',ETOT
        WRITE(*,FMT='(82("="),T10,"  DENMAT  ")')
        DO LMN1=1,LMNX
          WRITE(*,FMT='(I3,20F10.5)')LMN1,DENMAT(LMN1,:,1)
        ENDDO
        WRITE(*,FMT='(82("="),T10,"  HAMILTONIAN (W/O SCALING)  ")')
        DO LMN1=1,LMNX
          WRITE(*,FMT='(I3,20F10.5)')LMN1,DH(LMN1,:,1)
        ENDDO
        CALL ERROR$MSG('REGULAR STOP AFTER DIAGNOSTIC OUTPUT')
        CALL ERROR$STOP('SIMPLELMTO_CVX_ACTONPHI')
      END IF
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
      SUBROUTINE SIMPLELMTO_DC(ISP,NDIMD,LMNX_CHI,D_CHI,LMNX_PHI,D_PHI &
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
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1 
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
      LOGICAL(4)  ,PARAMETER  :: TPR=.FALSE.
      REAL(8)     ,PARAMETER  :: DELTA=1.D-2
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
      INTEGER(4)              :: IDIM,IR,LMN
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
      ALLOCATE(CUT(NR))
      CUT(:)=(RHO_CHI(:,1,1)/(RHO_PHI(:,1,1)+DELTA))**2 
!     == PROBLEM 1: A ZERO IN RHO_PHI AT LARGE DISTANCES PRODUCES A DIVERGENCE
!     ==   NEAR THE OUTER BOUNDARY OF THE GRID. THIS CAN BE REMEDIED BY 
!     ==   THE CHOICE OF DELTA, WHICH IS MUCH LARGER THAN RHO_CHI.
!     == PROBLEM 2: PARTIAL WAVES WITH NATURAL BOUNDARY CONDITIONS DO NOT 
!     ==   PRODUCE AN EFFECTIVE CUTOFF, SO THAT THE TAIL REGION DOMINATES.
!     == SOLUTION: CHOOSE DELTA IN THE RANGE OF THE DENSITY IN BETWEEN ATOMS.
!     ==                                                                      ==
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
!
!!$      CUT=1.D0
!!$      DO IR=1,NR
!!$        IF(R(IR).GT.RCOV) THEN
!!$          CUT(IR:)=0.D0
!!$          EXIT
!!$        END IF
!!$      ENDDO
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
!!$      DO IR=1,NR
!!$        IF(R(IR).GT.6.D0) THEN
!!$          VCUT(IR:)=0.D0
!!$          EXIT
!!$        END IF
!!$      ENDDO
!!$VCUT=0.D0
!VCUT=VCUT*Y0
!      
!     == POTENTIAL FOR PARTIAL-WAVE DENSITY N_T ================================
      POT_PHI(:,1,1)=POT_PHI(:,1,1)-VCUT(:)*2.D0*CUT(:)/(RHO_PHI(:,1,1)+DELTA)
!
!     == POTENTIAL FOR THE CORRELATED DENSITY ==================================
      POT_CHI(:,:,:)=0.D0
      POT_CHI(:,1,1)=VCUT(:)*2.D0*RHO_CHI(:,1,1)/(RHO_PHI(:,1,1)+DELTA)**2
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
        CALL LMTO_WRITEPHI('POT_CHI.DAT',GID,NR,LMRX,POT_CHI)
        CALL LMTO_WRITEPHI('POT_PHI.DAT',GID,NR,LMRX,POT_PHI)
        CALL LMTO_WRITEPHI('RHO_CHI.DAT',GID,NR,LMRX,RHO_CHI)
        CALL LMTO_WRITEPHI('RHO_PHI.DAT',GID,NR,LMRX,RHO_PHI)
        CALL LMTO_WRITEPHI('CUT.DAT',GID,NR,1,CUT)
        CALL LMTO_WRITEPHI('VCUT.DAT',GID,NR,1,VCUT)
!
        CALL REPORT$TITLE(6,'DIAGNOSTIC INFO FROM SIMPLELMTO_DC')
        CALL REPORT$I4VAL(6,'LMRX',LMRX,' ')
        CALL RADIAL$INTEGRAL(GID,NR &
    &                      ,FOURPI*R**2*(RHO_PHI(:,1,1)-AECORE(:))*Y0,SVAR)
        CALL REPORT$R8VAL(6,'VALENCE CHARGE (PARTIAL WAVE EXPANSION)',SVAR,' ')
        CALL RADIAL$INTEGRAL(GID,NR &
    &                       ,FOURPI*R**2*(RHO_CHI(:,1,1)-AECORE(:))*Y0,SVAR)
        CALL REPORT$R8VAL(6,'VALENCE CHARGE (LOCAL ORBITALS)',SVAR,' ')
        CALL RADIAL$INTEGRAL(GID,NR,FOURPI*R**2*AECORE*Y0,SVAR)
        CALL REPORT$R8VAL(6,'CORE CHARGE',SVAR,' ')
!
        WRITE(*,FMT='(82("="),T10,"  DENMAT-CHI  ")')
        DO LMN=1,LMNX_CHI
          WRITE(*,FMT='(I3,20F10.5)')LMN,D_CHI(LMN,:,1)
        ENDDO
        WRITE(*,FMT='(82("="),T10,"  DENMAT-PHI  ")')
        DO LMN=1,LMNX_PHI
          WRITE(*,FMT='(I3,20F10.5)')LMN,D_PHI(LMN,:,1)
        ENDDO
        WRITE(*,FMT='(82("="),T10,"  HAMILTONIAN-CHI (W/O SCALING)  ")')
        DO LMN=1,LMNX_CHI
          WRITE(*,FMT='(I3,20F10.5)')LMN,H_CHI(LMN,:,1)
        ENDDO
        WRITE(*,FMT='(82("="),T10,"  HAMILTONIAN-PHI (W/O SCALING)  ")')
        DO LMN=1,LMNX_PHI
          WRITE(*,FMT='(I3,20F10.5)')LMN,H_PHI(LMN,:,1)
        ENDDO
        CALL ERROR$MSG('REGULAR STOP AFTER PRINTING DIAGNOSTIC INFORMATION')
        CALL ERROR$STOP('SIMPLELMTO_DC')
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
        VCUT(IR)    =FXC(IR)  !/FOURPI FACTOR REMOVED BECAUSE NOT PRESENT IN
                                     ! PRIOR VERSION
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
!     ** COMPUTES OFFSITE MATRIX ELEMENTS OFFSITEX                            **
!     ** WATCH PARALLELIZATION!!!                                             **
!     **                                                                      **
!     ** CAUTION: THE CHOICE OF FIT FUNCTIONS FOR THE INTERPOLATION IS NOT    **
!     **    OPTIMIZED. LOGARITHMIC SPACING FOR DECAY AND RADIAL GRID MAY BE   **
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE,ONLY : POTPAR=>POTPAR1 &
     &                      ,OFFSITEX &
     &                      ,HYBRIDSETTING &
     &                      ,NSP
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      REAL(8)   ,PARAMETER :: TOLERANCE=1.D-5
      INTEGER(4),PARAMETER :: NDIS=25  !#(DISTANCE GRID POINTS)
      REAL(8)   ,PARAMETER :: DSMN=0.01D0,DSMX=8.D0 !PARMS FOR DISTANCE GRID
!                             !DISTANCE GRID WILL BE SCALED BY RCOV1+RCOV2 
!!$      INTEGER(4),PARAMETER :: NDIS=10  !#(DISTANCE GRID POINTS)
!!$      REAL(8)   ,PARAMETER :: DSMN=0.01D0,DSMX=10.D0 !PARMS FOR DISTANCE GRID
      INTEGER(4),PARAMETER :: NF=10  !#(FIT FUNCTIONS 1<NF<NDIS!)
      REAL(8)   ,PARAMETER :: DCAYMN=0.01D0,DCAYMX=0.3D0 !PARMS FOR FIT FUNCTIONS
      REAL(8)    ,PARAMETER:: RMAXEX=20.D0 !TARGET EXTENT OF EXTENDED RAD. GRIDS
      CHARACTER(8)         :: GRIDTYPE
      INTEGER(4)           :: GID1,GID2  ! GRID ID
      INTEGER(4)           :: NR1,NR2    ! #(RADIAL GRID POINTS)
      INTEGER(4)           :: LMNX1,LMNX2
      INTEGER(4)           :: LMN1,LMN2,LMN3,LMN4
      INTEGER(4)           :: ISP,ISP1,ISP2,I
      REAL(8)              :: SVAR,SVAR1,AEZ
      REAL(8)              :: RCOV(NSP)  ! COVALENT RADII
      REAL(8)              :: R1,DEX     ! GRID PARAMETERS
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
!     == DEFINE EXTENDED RADIAL GRIDS                                         ==
!     ==========================================================================
!     == THE RADIAL GRIDS SUPPLIED FOR THE SETUPS MAY BE TOO SHORT FOR THE    ==
!     == EVALUATION OF U-TENSOR MATRIX ELEMENTS. AN EXTENDED RADIAL GRID      ==
!     == IS THUS DEFINED FOR THE U-TENSOR. THE FIRST POINTS ON THE EXTENDED   ==
!     == GRID ARE IDENTICAL TO THOSE OF THE ORIGINAL GRID.                    ==
      DO ISP=1,NSP
        GID1=POTPAR(ISP)%GID
        CALL RADIAL$GETCH(GID1,'TYPE',GRIDTYPE)
        CALL RADIAL$NEW(GRIDTYPE,GID2)
        POTPAR(ISP)%GIDE=GID2
        CALL RADIAL$GETR8(GID1,'R1',R1)
        CALL RADIAL$SETR8(GID2,'R1',R1)
        CALL RADIAL$GETR8(GID1,'DEX',DEX)
        CALL RADIAL$SETR8(GID2,'DEX',DEX)
        CALL RADIAL$GETI4(GID1,'NR',NR1)
        IF(GRIDTYPE.EQ.'SHLOG') THEN
          NR2=MAX(NR1,NINT( LOG(1.D0+RMAXEX/R1)/DEX ))
        ELSE IF(GRIDTYPE.EQ.'LOG') THEN
          NR2=MAX(NR1,NINT( LOG(RMAXEX/R1)/DEX ))
        ELSE
          CALL ERROR$MSG('GRIDTYPE NOT RECOGNIZED: MUST BE "LOG" OR "SHLOG"')
          CALL ERROR$CHVAL('GRIDTYPE',GRIDTYPE)
          CALL ERROR$STOP('OFFXINT')
        END IF
        CALL RADIAL$SETI4(GID2,'NR',NR2)
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
        DO ISP2=1,NSP
          IF(.NOT.HYBRIDSETTING(ISP2)%TNDDO) CYCLE
PRINT*,'DOING X22 ....',ISP1,ISP2
          CALL SIMPLELMTO_OFFSITEX22SETUP(ISP1,ISP2,TOLERANCE)
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
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1 &
     &                       ,OFFSITEX
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      INTEGER(4),INTENT(IN) :: NR1
      INTEGER(4),INTENT(IN) :: NR2
      REAL(8)   ,INTENT(IN) :: TOLERANCE
      REAL(8)   ,PARAMETER  :: TOLMIN=1.D-8
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
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
      INTEGER(4)            :: I1,I2,IR
      INTEGER(4)            :: NTASKS,THISTASK,COUNT
      REAL(8)               ::SVAR
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(TPR) OPEN(UNIT=10,FILE='XOV.DAT')
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
!     == DIAGNOSTIC PRINTOUT                                                  ==
!     ==========================================================================
      IF(TPR) THEN
        IF(THISTASK.NE.1) THEN
          CALL ERROR$MSG('REGULAR STOP AFTER PRINTING DIAGNOSTIC INFORMATION')
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEOVERLAPSETUP')
        END IF 
        WRITE(10,*)'#XOV  FOR ',ISP1,ISP2
        WRITE(10,*)'#IND',IND,SHAPE(OFFSITEX(ISP1,ISP2)%X22)
        SVAR=0.D0
        DO I1=1,IND,10
          I2=MIN(IND,I1+9)
          WRITE(10,*)'#====',I1,I2,'====='
          WRITE(10,FMT='(11F10.5)')SVAR+OFFSITEX(ISP1,ISP2)%DIS(1) &
     &                           ,(/(0.D0,IR=1,10)/)
          DO IDIS=1,NDIS
            WRITE(10,FMT='(11F10.5)')SVAR+OFFSITEX(ISP1,ISP2)%DIS(IDIS) &
     &                              ,OFFSITEX(ISP1,ISP2)%OVERLAP(IDIS,I1:I2) &
     &                              ,(/(0.D0,IR=I2+1,I1+9)/)
          ENDDO
          WRITE(10,FMT='(11F10.5)')SVAR+OFFSITEX(ISP1,ISP2)%DIS(NDIS) &
     &                           ,(/(0.D0,IR=1,10)/)
          SVAR=SVAR+25.D0
        ENDDO
        CLOSE(10)
        CALL ERROR$MSG('REGULAR STOP AFTER PRINTING DIAGNOSTIC INFORMATION')
        CALL ERROR$STOP('SIMPLELMTO_OFFSITEOVERLAPSETUP')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_EXPLDISF(SCREENL,DIS,INV,DINV)
!     **************************************************************************
!     ** CALCULATES THE INTERPOLATING FUNCTIONS REPLACING THE LONG-RANGE-PART **
!     ** USED BY SIMPLELMTO_OFFSITEX22U AND SIMPLELMTO_OFFSITEX22SETUP        **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: SCREENL ! SCREENING LENGTH
      REAL(8)   ,INTENT(IN) :: DIS     ! DISTANCE OF THE MULTIPOLE CENTERS
      REAL(8)   ,INTENT(OUT):: INV(3)  ! REPLACEMENT FOR 1/D, 1/D**2, 1/D**3
      REAL(8)   ,INTENT(OUT):: DINV(3) ! DERIVATIVE OF INV(:)
      LOGICAL(4),PARAMETER  :: TOFF=.FALSE.
      REAL(8)   ,PARAMETER  :: LAMBDA=2.D0 
      REAL(8)               :: EXPLDIS
      REAL(8)               :: F1,DF1
!     **************************************************************************
      IF(TOFF) THEN
        INV=0.D0
        DINV=0.D0
        RETURN
      END IF     
      EXPLDIS=EXP(-LAMBDA*DIS)
      INV(1)=(1.D0-(1.D0-0.5D0*(LAMBDA*DIS)**2)*EXPLDIS)/DIS
      INV(2)=(1.D0-(1.D0+(LAMBDA*DIS)-(LAMBDA*DIS)**3/3.D0)*EXPLDIS)/DIS**2
      INV(3)=(1.D0-(1.D0+(LAMBDA*DIS)+(LAMBDA*DIS)**2/2.D0 &
     &                               -(LAMBDA*DIS)**4/8.D0)*EXPLDIS)/DIS**3
      DINV(1)=-INV(1)/DIS &
     &        -LAMBDA*(-LAMBDA*DIS)*EXPLDIS/DIS &
     &        +LAMBDA*(1.D0-0.5D0*(LAMBDA*DIS)**2)*EXPLDIS/DIS
      DINV(2)=-2.D0*INV(2)/DIS &
     &        -LAMBDA*(1.D0-(LAMBDA*DIS)**2)*EXPLDIS/DIS**2 &
     &        +LAMBDA*(1.D0+(LAMBDA*DIS)-(LAMBDA*DIS)**3/3.D0)*EXPLDIS/DIS**2
      DINV(3)=-3.D0*INV(3)/DIS &
     &        -LAMBDA*(1.D0+(LAMBDA*DIS)-(LAMBDA*DIS)**3/2.D0)*EXPLDIS/DIS**3 &
     &        +LAMBDA*(1.D0+(LAMBDA*DIS)+(LAMBDA*DIS)**2/2.D0 &
     &                                  -(LAMBDA*DIS)**4/8.D0)*EXPLDIS/DIS**3
      IF(SCREENL.GE.0.D0) THEN
        EXPLDIS=EXP(-DIS/SCREENL)
        F1=(1.D0+DIS/SCREENL)*EXPLDIS
        DF1=(EXPLDIS-F1)/SCREENL
!       __ATTENTION: DINV MUST BE CALCULATED BEFORE INV IS CHANGED!!____________
        DINV(1)=DINV(1)*F1+INV(1)*DF1
        DINV(2)=DINV(2)*F1+INV(2)*DF1
        DINV(3)=DINV(3)*F1+INV(3)*DF1
        INV(1)=INV(1)*F1
        INV(2)=INV(2)*F1
        INV(3)=INV(3)*F1
      ENDIF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEX22SETUP(ISP1,ISP2,TOLERANCE)
!     **************************************************************************
!     ** NDDO(NEGLECT OF DIFFERENTIAL OVERLAP) CONTRIBUTION TO THE EXCHANGE   **
!     ** THE COULOMB MATRIX ELEMENTS CONSIDERS ONE PAIR OF ORBITALS           **
!     ** ON ONE SITE R, WHICH FORM A DENSITY THAT INTERACTS WITH THE DENSITY  **
!     ** OF A SECOND PAIR OF ORBITALS ON A SECOND CENTER RPRIME               **
!     **                                                                      **
!     ** REMARK: THE LONG-RANGE TERM FROM MONOPOLE AND DIPOLE DENSITIES IS    **
!     **  SUBTRACTED OUT AND WILL BE ADDED IN AFTER INTERPOLATION             **
!     **                                                                      **
!     ** REMARK: USES AN EXTENDED RADIAL GRIDS.                               **
!     **                                                                      **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1 &
     &                             ,OFFSITEX &
     &                             ,SCREENL
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      REAL(8)   ,INTENT(IN) :: TOLERANCE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      LOGICAL(4),PARAMETER  :: TTEST=.FALSE. 
      CHARACTER(2),PARAMETER:: TEST_TYPE='TC' !'TWOCENTER'
!      CHARACTER(2),PARAMETER:: TEST_TYPE='CO' !'CORRECTION'
!      CHARACTER(2),PARAMETER:: TEST_TYPE='LR' !'LONGRANGE'
!      CHARACTER(2),PARAMETER:: TEST_TYPE='MD' !'MONODIPOLE'
      REAL(8)   ,PARAMETER  :: RX=0.1D0
      REAL(8)   ,PARAMETER  :: TOLMIN=1.D-8
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: SQ4PI=SQRT(4.D0*PI)
      REAL(8)   ,PARAMETER  :: SQ4PITHIRD=SQRT(4.D0*PI/3.D0)
      REAL(8)   ,PARAMETER  :: Y0=1.D0/SQ4PI
      INTEGER(4)            :: GID1,GID2   !GRID ID OF NORMAL GRIDS
      INTEGER(4)            :: GID1E,GID2E !GRID ID OF EXTENDED GRIDS
      INTEGER(4)            :: NR1,NR2     ! 
      INTEGER(4)            :: NR1E,NR2E
      REAL(8)   ,ALLOCATABLE:: RGRID1(:)  !(NR1E)
      REAL(8)   ,ALLOCATABLE:: RGRID2(:)  !(NR2E)
      INTEGER(4)            :: LRX1,LRX2
      INTEGER(4)            :: LNX1,LNX2
      INTEGER(4)            :: IND
      INTEGER(4)            :: LN1,LN2,LN3,LN4
      INTEGER(4)            :: L1,L2,L3,L4
      INTEGER(4)            :: LR1,LR2
      INTEGER(4)            :: LM1,LM2,LM3
      INTEGER(4)            :: MABS
      REAL(8)               :: INTEGRAL,DINTEGRAL
      INTEGER(4)            :: LMRX
      REAL(8)   ,ALLOCATABLE:: RHO12(:)   !(NR1E)
      REAL(8)   ,ALLOCATABLE:: RHO34(:)   !(NR2E)
      REAL(8)   ,ALLOCATABLE:: POT12(:)   !(NR1E)
      REAL(8)   ,ALLOCATABLE:: POT34(:)   !(NR2E)
      INTEGER(4)            :: IDIS
      INTEGER(4)            :: NDIS
      REAL(8)               :: DIS
      REAL(8) ,ALLOCATABLE  :: TOLFAC1(:)
      REAL(8) ,ALLOCATABLE  :: TOLFAC2(:)
      REAL(8) ,ALLOCATABLE  :: YLMDIS(:)
      REAL(8)               :: TOL
      INTEGER(4)            :: IR,I1,I2
      REAL(8)               :: SVAR,SVAR1,SVAR2,A,B,YLM2,FAC
      REAL(8)               :: INV(3),DINV(3)
      INTEGER(4)            :: THISTASK,NTASKS,COUNT
LOGICAL(4),PARAMETER  :: TPTCHM=.TRUE. ! POINT-CHARGE MODEL
      REAL(8) :: Q1,Q2,D1(3),D2(3),E,VQ1,VQ2,VD1(3),VD2(3),LAMBDA
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(TPR) OPEN(UNIT=10,FILE='X22.DAT')
!
!     ==========================================================================
!     == PREPARATION: RADIAL GRIDS                                            ==
!     ==========================================================================
      GID1=POTPAR(ISP1)%GID
      GID2=POTPAR(ISP2)%GID
      CALL RADIAL$GETI4(GID1,'NR',NR1)
      CALL RADIAL$GETI4(GID2,'NR',NR2)
      GID1E=POTPAR(ISP1)%GIDE
      GID2E=POTPAR(ISP2)%GIDE
      CALL RADIAL$GETI4(GID1E,'NR',NR1E)
      CALL RADIAL$GETI4(GID2E,'NR',NR2E)
      ALLOCATE(RGRID1(NR1E))
      ALLOCATE(RGRID2(NR2E))
      CALL RADIAL$R(GID1E,NR1E,RGRID1)
      CALL RADIAL$R(GID2E,NR2E,RGRID2)
!
!     ==========================================================================
!     == PREPARATION: ANGULAR MOMENTA                                         ==
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

      ALLOCATE(YLMDIS((LRX2+2)**2))
      CALL SPHERICAL$YLM((LRX2+2)**2,(/0.D0,0.D0,1.D0/),YLMDIS)
!
!     ==========================================================================
!     == OBTAIN NORM TO SET TOLERANCE                                         ==
!     ==========================================================================
      ALLOCATE(TOLFAC1(LNX1))
      ALLOCATE(TOLFAC2(LNX2))
      DO LN1=1,LNX1
        CALL RADIAL$INTEGRAL(GID1,NR1 &
     &                ,(RGRID1(:NR1)*POTPAR(ISP1)%AECHI(:,LN1))**2,TOLFAC1(LN1))
      ENDDO
      DO LN2=1,LNX2
        CALL RADIAL$INTEGRAL(GID2,NR2 &
    &                 ,(RGRID2(:NR2)*POTPAR(ISP2)%AECHI(:,LN2))**2,TOLFAC2(LN2))
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
                    IF(TPR)WRITE(10,*)'# IND',IND,MABS,LR1,LR2,LN1,LN2,LN3,LN4
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      NDIS=OFFSITEX(ISP1,ISP2)%NDIS
      ALLOCATE(OFFSITEX(ISP1,ISP2)%X22(NDIS,IND))
!
!     ==========================================================================
!     == DETERMINE INTEGRALS                                                 ==
!     ==========================================================================
      ALLOCATE(RHO12(NR1E))
      ALLOCATE(POT12(NR1E))
      ALLOCATE(RHO34(NR2E))
      ALLOCATE(POT34(NR2E))
      OFFSITEX(ISP1,ISP2)%X22=0.D0
      COUNT=0
      IND=0
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%LOXH(LN2)
          RHO12(:NR1)=POTPAR(ISP1)%AECHI(:,LN1) &
       &             *POTPAR(ISP1)%AECHI(:,LN2)
          RHO12(NR1+1:)=0.D0
          DO LN3=1,LNX2
            L3=POTPAR(ISP2)%LOXH(LN3)
            DO LN4=LN3,LNX2
              L4=POTPAR(ISP2)%LOXH(LN4)
              RHO34(:NR2)=POTPAR(ISP2)%AECHI(:,LN3) &
       &                 *POTPAR(ISP2)%AECHI(:,LN4)
              RHO34(NR2+1:)=0.D0
              TOL=TOLERANCE*TOLFAC1(LN1)*TOLFAC1(LN2) &
       &                   *TOLFAC2(LN3)*TOLFAC2(LN4)
              TOL=MAX(TOLMIN,TOL)
!
              DO LR1=ABS(L1-L2),MIN(L1+L2,LRX1),2
!               ================================================================
!               == DETERMINE POTENTIAL OF FIRST SITE
!               ================================================================
                IF(SCREENL.LT.0.D0) THEN ! SCREENL<0 IMPLIES SCREENL=INFINITE
                  CALL RADIAL$POISSON(GID1E,NR1E,LR1,RHO12,POT12)
                ELSE
!                 __AVOID UNDERFLOW IN YUKAWA FOR TOO SMALL SCREENL_____________
                  IF(SCREENL.LT.3.779D-2) THEN
                    CALL ERROR$MSG('SET SCREENING LENGTH LARGER THAN 0.02 AA')
                    CALL ERROR$MSG('IN !CONTROL!DFT!NTBO:SCREENL[AA]')
                    CALL ERROR$MSG('(NUMERICAL PROBLEM IN LIBRARY ROUTINE)')
                    CALL ERROR$STOP('SIMPLELMTO_OFFSITEX22SETUP')
                  END IF
                  SVAR=1.D0/MAX(1.D-8,SCREENL) !AVOID DIVIDE BY ZERO
                  CALL RADIAL$YUKAWA(GID1E,NR1E,LR1,SVAR,RHO12,POT12)
                END IF

IF(TTEST) THEN
  POT12=0.D0
  IF(LR1.GT.1) RHO12=0.D0
END IF

!               ==SUBTRACT OUT LONG RANGE PART INCLUDING MONO- AND DIPOLE TERMS
                IF(LR1.EQ.0) THEN
                  SVAR=POTPAR(ISP1)%QLN(1,LN1,LN2)*SQ4PI
                  IF(SCREENL.LE.0.D0) THEN
                    POT12(2:)=POT12(2:)-SVAR/RGRID1(2:)
                    POT12(1)=POT12(2)
                  ELSE
                    POT12(2:)=POT12(2:)-SVAR/RGRID1(2:)*EXP(-RGRID1(2:)/SCREENL)
                    POT12(1)=POT12(2)
                  END IF
                ELSE IF(LR1.EQ.1) THEN
                  SVAR=POTPAR(ISP1)%QLN(2,LN1,LN2)*SQ4PITHIRD
                  IF(SCREENL.LE.0.D0) THEN ! WITHOUT SCREENING
                    POT12(2:)=POT12(2:)-SVAR/RGRID1(2:)**2
                    POT12(1)=POT12(2)
                 ELSE                      ! WITH SCREENING
                    POT12(2:)=POT12(2:)-SVAR/RGRID1(2:)**2 &
      &                      *(1.D0+RGRID1(2:)/SCREENL)*EXP(-RGRID1(2:)/SCREENL)
                    POT12(1)=POT12(2)
                  END IF
                END IF
!
IF(TTEST)POT12=-POT12

!
                DO LR2=ABS(L3-L4),MIN(L3+L4,LRX2),2
!                 ==============================================================
!                 == DETERMINE POTENTIAL OF SECOND SITE
!                 ==============================================================
                  IF(SCREENL.LT.0.D0) THEN ! SCREENL<0 IMPLIES SCREENL=INFINITE
                    CALL RADIAL$POISSON(GID2E,NR2E,LR2,RHO34,POT34)
                  ELSE
!                   __AVOID UNDERFLOW IN YUKAWA FOR TOO SMALL SCREENL___________
                    IF(SCREENL.LT.3.779D-2) THEN
                      CALL ERROR$MSG('SET SCREENING LENGTH LARGER THAN 0.02 AA')
                      CALL ERROR$MSG('IN !CONTROL!DFT!NTBO:SCREENL[AA]')
                      CALL ERROR$MSG('(NUMERICAL PROBLEM IN LIBRARY ROUTINE)')
                      CALL ERROR$STOP('SIMPLELMTO_OFFSITEX22SETUP')
                    END IF
                    SVAR=1.D0/MAX(1.D-8,SCREENL) !AVOID DIVIDE BY ZERO
                    CALL RADIAL$YUKAWA(GID2E,NR2E,LR2,SVAR,RHO34,POT34)
                  END IF
!                 ==SUBTRACT LONG RANGE PART INCLUDING MONO- AND DIPOLE TERMS
IF(TTEST) THEN
  POT34=0.D0
  IF(LR2.GT.1) RHO34=0.D0
END IF
                  IF(LR2.EQ.0) THEN
                    SVAR=POTPAR(ISP2)%QLN(1,LN3,LN4)*SQ4PI
                    IF(SCREENL.LE.0.D0) THEN
                      POT34(2:)=POT34(2:)-SVAR/RGRID2(2:)
                      POT34(1)=POT34(2)
                    ELSE
                      POT34(2:)=POT34(2:) &
       &                       -SVAR/RGRID2(2:)*EXP(-RGRID2(2:)/SCREENL)
                      POT34(1)=POT34(2)
                    END IF
                  ELSE IF(LR2.EQ.1) THEN
                    SVAR=POTPAR(ISP2)%QLN(2,LN3,LN4)*SQ4PITHIRD
                    IF(SCREENL.LE.0.D0) THEN ! WITHOUT SCREENING
                      POT34(2:)=POT34(2:)-SVAR/RGRID2(2:)**2
                      POT34(1)=POT34(2)
                    ELSE                      ! WITH SCREENING
                      POT34(2:)=POT34(2:)-SVAR/RGRID2(2:)**2 &
      &                      *(1.D0+RGRID2(2:)/SCREENL)*EXP(-RGRID2(2:)/SCREENL)
                      POT34(1)=POT34(2)
                    END IF
                  END IF
!
IF(TTEST) POT34=-POT34
!
!                 ==============================================================
!                 == SUM OVER M-VALUES OF                                     ==
!                 == AXIAL ANGULAR MOMENTA OF THE DENSITY ALONG BOND AXIS     ==
!                 ==============================================================
                  DO MABS=0,MIN(LR1,LR2)
! THIS IS A SUM OVER THE ABSOLUTE VALUES OF M. DO I NEED A FACTOR TWO FOR 
! MABS>0 TO ACCOUNT FOR M=-MABS?
                    IND=IND+1
                    DO IDIS=1,NDIS
                      COUNT=COUNT+1
!                     == OFFSITEX WILL BE ADDED TOGETHER IN THE CALLING ROUTINE
                      IF(MOD(COUNT-1,NTASKS).NE.THISTASK-1) CYCLE
                      DIS=OFFSITEX(ISP1,ISP2)%DIS(IDIS)
                      CALL SIMPLELMTO_TWOCENTER(LR1,MABS,GID1E,NR1E,POT12 &
       &                                       ,LR2,MABS,GID2E,NR2E,RHO34 &
       &                                       ,DIS,TOL,INTEGRAL)
!
IF(TTEST.AND.TEST_TYPE.EQ.'TC') GOTO 1000
IF(TTEST.AND.TEST_TYPE.EQ.'CO')    INTEGRAL=0.D0

                      IF(MABS.EQ.0) THEN
                         IF(LR1.EQ.0) THEN
                           CALL RADIAL$VALUE(GID2E,NR2E,POT34,DIS,SVAR)
                           SVAR=SVAR*POTPAR(ISP1)%QLN(LR1+1,LN1,LN2)
                           YLM2=(-1.D0)**LR2 * YLMDIS(LR2**2+LR2+1)
                           INTEGRAL=INTEGRAL+SVAR*YLM2
                         ELSE IF(LR1.EQ.1) THEN
                           CALL RADIAL$DERIVATIVE(GID2E,NR2E,POT34,DIS,SVAR)
                           SVAR=SVAR*POTPAR(ISP1)%QLN(LR1+1,LN1,LN2)
                           YLM2=(-1.D0)**(LR2+1) * YLMDIS(LR2**2+LR2+1)
                           INTEGRAL=INTEGRAL+SVAR*YLM2
                         END IF
                      END IF
                      IF(LR1.EQ.1) THEN
                        CALL RADIAL$VALUE(GID2E,NR2E,POT34,DIS,SVAR)
!!$IF(MABS.EQ.0.AND.LR1.EQ.1) SVAR=-SVAR
!!$! FEHLER FUER MABS=1,LR1=1,LR2=1
!!$IF(MABS.EQ.1.AND.LR1.EQ.1.AND.LR2.EQ.1) SVAR=2.D0*SVAR
                        SVAR=SVAR*POTPAR(ISP1)%QLN(LR1+1,LN1,LN2)
                        SVAR=SVAR*SQ4PITHIRD/DIS
                        LM1=3-MABS
                        LM2=LR2**2+LR2+1-MABS
                        LM3=(LR2+1)**2+(LR2+1)+1
                        CALL SPHERICAL$GAUNT(LM1,LM2,LM3,SVAR1)
                        YLM2=(-1.D0)**(LR2+1) * YLMDIS(LM3)
                        SVAR1=SVAR1*REAL(-LR2,KIND=8)*YLM2
                        IF(LR2.GE.1) THEN
                          LM3=(LR2-1)**2+(LR2-1)+1
                          CALL SPHERICAL$GAUNT(LM1,LM2,LM3,SVAR2)
                          YLM2=(-1.D0)**(LR2-1) * YLMDIS(LM3)
                          SVAR2=SVAR2*REAL(LR2+1,KIND=8)*YLM2
                        ELSE
                          SVAR2=0.D0 
                        END IF
                        INTEGRAL=INTEGRAL+SVAR*(SVAR1+SVAR2)
                      END IF
!
!                     ==========================================================
!                     == ADD IN U_LONGRANGE PART TO COMPLETE MATRIX ELEMENT   ==
!                     == THEN, SUBTRACT OUT LONG RANGE CORRECTION TO ALLOW    ==
!                     == INTERPOLATION. THE CORRECTION                        ==
!                     == WILL BE ADDED AGAIN IN SIMPLELMTO_OFFSITEX22U        ==
!                     ==========================================================
IF(TTEST.AND.TEST_TYPE.EQ.'CO')    GOTO 1000
!IF(TTEST.AND.TEST_TYPE.EQ.'LR') THEN
!INTEGRAL=0.D0
!
                      IF(LR1.LE.1.AND.LR2.LE.1) THEN
                        FAC=POTPAR(ISP1)%QLN(LR1+1,LN1,LN2) &
     &                      *POTPAR(ISP2)%QLN(LR2+1,LN3,LN4)
                        SVAR=FAC
                        CALL SIMPLELMTO_EXPLDISF(SCREENL,DIS,INV,DINV)
                        IF(LR1+LR2.EQ.0) THEN ! MONOPOLE-MONOPOLE
                          IF(SCREENL.LE.0.D0) THEN
                            INTEGRAL =INTEGRAL +SVAR/DIS
                          ELSE
                            INTEGRAL =INTEGRAL +SVAR/DIS*EXP(-DIS/SCREENL)
                          END IF
                          INTEGRAL =INTEGRAL -FAC*INV(1)
                        ELSE IF(LR1+LR2.EQ.1) THEN  !MONOPOLE-DIPOLE
                          SVAR=FAC
                          IF(LR1.EQ.1) SVAR=-SVAR
                          IF(SCREENL.LE.0.D0) THEN
                            SVAR=-SVAR/DIS**2
                            INTEGRAL =INTEGRAL +SVAR
                          ELSE
                            SVAR=-SVAR/DIS**2 &
        &                             *EXP(-DIS/SCREENL)*(1.D0+DIS/SCREENL)
                            INTEGRAL =INTEGRAL +SVAR
                          END IF
                          SVAR=FAC
                          IF(LR1.EQ.1) SVAR=-SVAR
                            INTEGRAL =INTEGRAL -SVAR*INV(2)
                        ELSE IF(LR1+LR2.EQ.2) THEN  !DIPOLE-DIPOLE 
                          SVAR=FAC/DIS**3
                          IF(SCREENL.LE.0.D0) THEN
                            INTEGRAL =INTEGRAL +SVAR
                          ELSE
                            INTEGRAL =INTEGRAL +SVAR*(1.D0+DIS/SCREENL) &
        &                                            *EXP(-DIS/SCREENL)
                          END IF

                          INTEGRAL =INTEGRAL -FAC*INV(3)
                          IF(MABS.EQ.0) THEN
                             IF(SCREENL.LE.0.D0) THEN
                               INTEGRAL =INTEGRAL -SVAR*3.D0
                             ELSE
                               INTEGRAL=INTEGRAL -SVAR &
         &                           *(3.D0+3.D0*DIS/SCREENL+(DIS/SCREENL)**2) &
         &                           *EXP(-DIS/SCREENL)
                             END IF
                            INTEGRAL =INTEGRAL -3.D0*FAC*INV(3)
                          END IF
                        END IF
                      END IF
!!$   GOTO 1000
!!$ END IF
!
IF(TTEST.AND.TEST_TYPE.EQ.'MD') THEN
  Q1=0.D0
  D1=0.D0
  Q2=0.D0
  D2=0.D0
  IF(LR1.EQ.0) THEN
    Q1=POTPAR(ISP1)%QLN(LR1+1,LN1,LN2) 
  ELSE IF(LR1.EQ.1) THEN
    IF(MABS.EQ.0) THEN
      D1=(/0.D0,0.D0,1.D0/)*POTPAR(ISP1)%QLN(LR1+1,LN1,LN2)
    ELSE IF(MABS.EQ.1) THEN
      D1=(/1.D0,0.D0,0.D0/)*POTPAR(ISP1)%QLN(LR1+1,LN1,LN2)
    END IF
  END IF
  IF(LR2.EQ.0) THEN
    Q2=POTPAR(ISP2)%QLN(LR2+1,LN3,LN4) 
  ELSE IF(LR2.EQ.1) THEN
    IF(MABS.EQ.0) THEN
      D2=(/0.D0,0.D0,1.D0/)*POTPAR(ISP2)%QLN(LR2+1,LN3,LN4) 
    ELSE IF(MABS.EQ.1) THEN
      D2=(/1.D0,0.D0,0.D0/)*POTPAR(ISP2)%QLN(LR2+1,LN3,LN4) 
    END IF
  END IF
  LAMBDA=0.D0
  IF(SCREENL.GT.0.D0)LAMBDA=1.D0/SCREENL
  CALL MONOANDDIPOLE(Q1,D1,Q2,D2,LAMBDA,DIS,INTEGRAL,VQ1,VD1,VQ2,VD2)
END IF

1000 CONTINUE
IF(TTEST.AND.ABS(INTEGRAL).LT.1.D-5) INTEGRAL=0.D0


                      OFFSITEX(ISP1,ISP2)%X22(IDIS,IND)=INTEGRAL
                    ENDDO   ! LOOP IDIS
                  ENDDO   ! LOOP MABS
                ENDDO   !LOOP LR2
              ENDDO   ! LOOP LR1
            ENDDO   ! LOOP LN4
          ENDDO   ! LOOP LN3
        ENDDO   ! LOOP LN2
      ENDDO   ! LOOP LN1
      CALL MPE$COMBINE('MONOMER','+',OFFSITEX(ISP1,ISP2)%X22)
!
!     ==========================================================================
!     == DIAGNOSTIC PRINTOUT                                                  ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(10,*)'#X22  FOR ',ISP1,ISP2
        WRITE(10,*)'#IND',IND,SHAPE(OFFSITEX(ISP1,ISP2)%X22)
        SVAR=0.D0
        DO I1=1,IND,10
          I2=MIN(IND,I1+9)
          WRITE(10,*)'#====',I1,I2,'====='
          WRITE(10,FMT='(11F10.5)')SVAR+OFFSITEX(ISP1,ISP2)%DIS(1) &
     &                           ,(/(0.D0,IR=1,10)/)
          DO IDIS=1,NDIS
            WRITE(10,FMT='(11F10.5)')SVAR+OFFSITEX(ISP1,ISP2)%DIS(IDIS) &
     &                          ,OFFSITEX(ISP1,ISP2)%X22(IDIS,I1:I2)/(4.D0*PI) &
     &                             ,(/(0.D0,IR=I2+1,I1+9)/)
          ENDDO
          WRITE(10,FMT='(11F10.5)')SVAR+OFFSITEX(ISP1,ISP2)%DIS(NDIS) &
     &                           ,(/(0.D0,IR=1,10)/)
          SVAR=SVAR+25.D0
        ENDDO
        CLOSE(10)
        CALL ERROR$MSG('REGULAR STOP AFTER PRINTING DIAGNOSTIC INFORMATION')
        CALL ERROR$STOP('SIMPLELMTO_OFFSITEX22SETUP')
      END IF
!
!     ==========================================================================
!     == CLEANUP                                                              ==
!     ==========================================================================
      DEALLOCATE(YLMDIS)
      DEALLOCATE(TOLFAC1)
      DEALLOCATE(TOLFAC2)
      DEALLOCATE(RGRID1)
      DEALLOCATE(RGRID2)
      DEALLOCATE(RHO12)
      DEALLOCATE(RHO34)
      DEALLOCATE(POT12)
      DEALLOCATE(POT34)
      RETURN
      END

!
!    ....1.........2.........3.........4.........5.........6.........7.........8
     SUBROUTINE MONOANDDIPOLE(Q1,D1,Q2,D2,LAMBDA,DIS12,E,VQ1,VD1,VQ2,VD2)
!    ***************************************************************************
!    ** MONOPOLE Q1 AND DIPOLE D1 ARE AT THE ORIGIN R1=(0,0,0)                **
!    ** MONOPOLE Q2 AND DIPOLE D2 ARE AT R2=(0,0,Z)                           **
!    ***************************************************************************
     IMPLICIT NONE
     REAL(8),INTENT(IN)  :: Q1      ! 1ST MONOPOLE AT R1
     REAL(8),INTENT(IN)  :: D1(3)   ! 1ST DIPOLE   AT R1
     REAL(8),INTENT(IN)  :: Q2      ! 2ND MONOPOLE AT R2=R1+(0,0,DIS12)
     REAL(8),INTENT(IN)  :: D2(3)   ! 2ND DIPOLE   AT R2=R1+(0,0,DIS12)
     REAL(8),INTENT(IN)  :: LAMBDA  ! SCREENING PARAMETER
     REAL(8),INTENT(IN)  :: DIS12   ! DISTANCE |R2-R1|
     REAL(8),INTENT(OUT) :: E       ! ENERGY
     REAL(8),INTENT(OUT) :: VQ1     ! POTENTIAL FROM Q2,D2 AT R1
     REAL(8),INTENT(OUT) :: VD1(3)  ! GRADIENT OF THE POTENTIAL FROM Q2,D2 AT R1
     REAL(8),INTENT(OUT) :: VQ2     ! POTENTIAL FROM Q1,D1 AT R2
     REAL(8),INTENT(OUT) :: VD2(3)  ! GRADIENT OF THE POTENTIAL FROM Q1,D1 AT R2
!    ***************************************************************************
     E=Q1*Q2/DIS12 &
    & +(1.D0+LAMBDA*DIS12)*(D1(3)*Q2-Q1*D2(3))/DIS12**2 &
    & +(1.D0+LAMBDA*DIS12)*(SUM(D1*D2)-3.D0*D1(3)*D2(3))/DIS12**3 &
    & -(LAMBDA*DIS12)**2*D1(3)*D2(3)/DIS12**3
     E=E*EXP(-LAMBDA*DIS12)
!    == VQ1=DE/DQ1 =============================================================
     VQ1=( Q2/DIS12-(1.D0+LAMBDA*DIS12)*D2(3)/DIS12**2 )*EXP(-LAMBDA*DIS12)
!
!    == VD1(I)=DE/DD1(I) = GRAD V2(R1) =========================================
     VD1(:)=(1.D0+LAMBDA*DIS12)*D2(:)/DIS12**3
     VD1(3)=VD1(3)+(1.D0+LAMBDA*DIS12)*Q2/DIS12**2 &
    &             -(1.D0+LAMBDA*DIS12)*3.D0*D2(3)/DIS12**3 &
    &             -(LAMBDA*DIS12)**2*D2(3)/DIS12**3
     VD1=VD1*EXP(-LAMBDA*DIS12)
!
!    == VQ2=DE/DQ2 =============================================================
     VQ2=( Q1/DIS12+(1.D0+LAMBDA*DIS12)*D1(3)/DIS12**2 )*EXP(-LAMBDA*DIS12)
!
!    == VD2(I)=DE/DD2(I) =GRAD V1(R2)===========================================
     VD2(:)=(1.D0+LAMBDA*DIS12)*D1(:)/DIS12**3 
     VD2(3)=VD2(3)-(1.D0+LAMBDA*DIS12)*Q1/DIS12**2 &
    &             -(1.D0+LAMBDA*DIS12)*3.D0*D1(3)/DIS12**3 &
    &             -(LAMBDA*DIS12)**2*D1(3)/DIS12**3
     VD2=VD2*EXP(-LAMBDA*DIS12)
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
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1 &
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
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
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
      REAL(8)               :: SVAR
      REAL(8)               :: DIS
      REAL(8) ,ALLOCATABLE  :: TOLFAC1(:)
      REAL(8) ,ALLOCATABLE  :: TOLFAC2(:)
      REAL(8)               :: TOL
      INTEGER(4)            :: NTASKS,THISTASK,COUNT
      INTEGER(4)            :: I1,I2,IR
!     **************************************************************************
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(TPR) OPEN(UNIT=10,FILE='X31.DAT')
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
                IF(SCREENL.LT.0.D0) THEN ! SCREENL<0 IMPLIES SCREENL=INFINITE
                  CALL RADIAL$POISSON(GID1,NR1,LR1,RHO12,POT12)
                ELSE
!                 __AVOID UNDERFLOW IN YUKAWA FOR TOO SMALL SCREENL_____________
                  IF(SCREENL.LT.3.779D-2) THEN
                    CALL ERROR$MSG('SET SCREENING LENGTH LARGER THAN 0.02 AA')
                    CALL ERROR$MSG('IN !CONTROL!DFT!NTBO:SCREENL[AA]')
                    CALL ERROR$MSG('(NUMERICAL PROBLEM IN LIBRARY ROUTINE)')
                    CALL ERROR$STOP('SIMPLELMTO_OFFSITEX31SETUP')
                  END IF
                  SVAR=1.D0/MAX(1.D-8,SCREENL) !AVOID DIVIDE BY ZERO
                  CALL RADIAL$YUKAWA(GID1,NR1,LR1,SVAR,RHO12,POT12)
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
!!$!
!!$!     ==========================================================================
!!$!     == DIAGNOSTIC PRINTOUT                                                  ==
!!$!     ==========================================================================
!!$      IF(TPR) THEN
!!$        PRINT*,'X31  FOR ',ISP1,ISP2
!!$        DO IDIS=1,NDIS
!!$          PRINT*,OFFSITEX(ISP1,ISP2)%DIS(IDIS) &
!!$     &          ,OFFSITEX(ISP1,ISP2)%X31(IDIS,:) !/(4.D0*PI)
!!$        ENDDO
!!$        STOP 'FORCED'
!!$      END IF
!
!     ==========================================================================
!     == DIAGNOSTIC PRINTOUT                                                  ==
!     ==========================================================================
      IF(TPR) THEN
        WRITE(10,*)'#X31  FOR ',ISP1,ISP2
        WRITE(10,*)'#IND',IND,SHAPE(OFFSITEX(ISP1,ISP2)%X31)
        SVAR=0.D0
        DO I1=1,IND,10
          I2=MIN(IND,I1+9)
          WRITE(10,*)'#====',I1,I2,'====='
          WRITE(10,FMT='(11F10.5)')SVAR+OFFSITEX(ISP1,ISP2)%DIS(1) &
     &                           ,(/(0.D0,IR=1,10)/)
          DO IDIS=1,NDIS
            WRITE(10,FMT='(11F10.5)')SVAR+OFFSITEX(ISP1,ISP2)%DIS(IDIS) &
     &                          ,OFFSITEX(ISP1,ISP2)%X31(IDIS,I1:I2)/(4.D0*PI) &
     &                             ,(/(0.D0,IR=I2+1,I1+9)/)
          ENDDO
          WRITE(10,FMT='(11F10.5)')SVAR+OFFSITEX(ISP1,ISP2)%DIS(NDIS) &
     &                           ,(/(0.D0,IR=1,10)/)
          SVAR=SVAR+25.D0
        ENDDO
        CLOSE(10)
        CALL ERROR$MSG('REGULAR STOP AFTER PRINTING DIAGNOSTIC INFORMATION')
        CALL ERROR$STOP('SIMPLELMTO_OFFSITEX31SETUP')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEXCONVERT()
!     **************************************************************************
!     ** CONVERT THE DISTANCE-DEPENDENT MATRIX ELEMENTS OF THE U-TENSOR       **
!     ** INTO EXPANSION COEFFICIENTS FOR THE INTERPOLATING FUNCTION           **
!     ** F_J(X)=SUM_I=1^NDIS X22(I,J)*EXP(-LAMBDA(I)*X)*(1+LAMBDA(I)*X)       **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : OFFSITEX &
     &                             ,POTPAR=>POTPAR1 &
     &                             ,NSP
      IMPLICIT NONE
      INTEGER(4)            :: ISP1,ISP2
      INTEGER(4)            :: NDIS   !#(GRID POINTS)
      INTEGER(4)            :: NF     !#(INTERPOLATING FUNCTIONS)
      INTEGER(4)            :: I
      INTEGER(4)            :: NIND   !#(INTEGRALS)
      LOGICAL(4),PARAMETER  :: TESTNDDOFIT=.FALSE.
      LOGICAL(4),PARAMETER  :: TESTX31FIT=.FALSE.
      REAL(8)   ,ALLOCATABLE:: AMAT(:,:)
      REAL(8)   ,ALLOCATABLE:: AINV(:,:)
      REAL(8)   ,ALLOCATABLE:: AMATIN(:,:)
      REAL(8)   ,ALLOCATABLE:: G(:,:)  !INTERPOLATING FUNCTIONS
      REAL(8)               :: LAMBDA1,LAMBDA2,SVAR1,SVAR2
      REAL(8)               :: DISX
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
     &                 *OFFSITEX(ISP1,ISP2)%DIS(:)) &
     &            *(1.D0+OFFSITEX(ISP1,ISP2)%LAMBDA(I) &
     &                  *OFFSITEX(ISP1,ISP2)%DIS(:))
          ENDDO
          TFIT=OFFSITEX(ISP1,ISP2)%NF.NE.OFFSITEX(ISP1,ISP2)%NDIS
          ALLOCATE(AMAT(NF,NF))
          DISX=MAXVAL(OFFSITEX(ISP1,ISP2)%DIS)
          IF(TFIT) THEN
            DO I=1,NF
              LAMBDA1=OFFSITEX(ISP1,ISP2)%LAMBDA(I) 
              DO J=I,NF
                LAMBDA2=OFFSITEX(ISP1,ISP2)%LAMBDA(J) 
                SVAR1=LAMBDA1*LAMBDA2/(LAMBDA1+LAMBDA2)**2
                SVAR2=(LAMBDA1+LAMBDA2)*DISX
                AMAT(I,J)=SUM(G(:,I)*G(:,J)) &
     &              +(2.D0+SVAR2+SVAR1*(SVAR2**2+2.D0*SVAR2+2.D0)) &
     &               *EXP(-SVAR2)/SVAR2*REAL(NDIS,KIND=8)
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
            IF(TESTNDDOFIT) THEN !----------------------------------------------
              PRINT*,'NDDO-A:',OFFSITEX(ISP1,ISP2)%X22(:NDIS,:)
              CALL SIMPLELMTO_TESTOFFSITEXCONVERT_A(NDIS &
        &                                     ,OFFSITEX(ISP1,ISP2)%DIS &
        &                                     ,OFFSITEX(ISP1,ISP2)%X22(:,1))
            END IF !------------------------------------------------------------

            NIND=SIZE(OFFSITEX(ISP1,ISP2)%X22(1,:))
            DO I=1,NIND
              OFFSITEX(ISP1,ISP2)%X22(:NF,I)=MATMUL(AMATIN &
      &                                          ,OFFSITEX(ISP1,ISP2)%X22(:,I))
            ENDDO
            OFFSITEX(ISP1,ISP2)%X22(NF+1:,:)=0.D0
!
            IF(TESTNDDOFIT) THEN !----------------------------------------------
              PRINT*,'NDDO-B:',OFFSITEX(ISP1,ISP2)%X22(:NF,:)
              CALL SIMPLELMTO_TESTOFFSITEXCONVERT_B(NF &
       &                   ,OFFSITEX(ISP1,ISP2)%LAMBDA(:NF) &
       &                   ,OFFSITEX(ISP1,ISP2)%X22(:NF,1))
            END IF !------------------------------------------------------------
          END IF
!
!         ======================================================================
!         == 31 TERMS  X31
!         ======================================================================
          IF(T31) THEN
            IF(TESTX31FIT) THEN !----------------------------------------------
              PRINT*,'X31-A:',OFFSITEX(ISP1,ISP2)%X22(:NDIS,:)
              CALL SIMPLELMTO_TESTOFFSITEXCONVERT_A(NDIS &
        &                                     ,OFFSITEX(ISP1,ISP2)%DIS &
        &                                     ,OFFSITEX(ISP1,ISP2)%X31(:,1))
            END IF !------------------------------------------------------------
!
            NIND=SIZE(OFFSITEX(ISP1,ISP2)%X31(1,:))
            DO I=1,NIND
              OFFSITEX(ISP1,ISP2)%X31(:NF,I)=MATMUL(AMATIN &
       &                                         ,OFFSITEX(ISP1,ISP2)%X31(:,I))
            ENDDO
            OFFSITEX(ISP1,ISP2)%X31(NF+1:,:)=0.D0
!
            IF(TESTX31FIT) THEN !----------------------------------------------
              PRINT*,'X31-B:',OFFSITEX(ISP1,ISP2)%X22(:NF,:)
              CALL SIMPLELMTO_TESTOFFSITEXCONVERT_B(NF &
       &                   ,OFFSITEX(ISP1,ISP2)%LAMBDA(:NF) &
       &                   ,OFFSITEX(ISP1,ISP2)%X31(:NF,1))
            END IF !------------------------------------------------------------
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
!     ** TEST_A.DAT, WHILE THE SECOND WRIMORE TEST_TES THE INTERPOLATING FUNCTION TO    **
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
      REAL(8)   ,PARAMETER  :: DMAX=50.D0
      INTEGER(4),PARAMETER  :: ND=1000
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
        FINT=SUM(X(:)*EXP(-LAMBDA(:)*D)*(1.D0+LAMBDA(:)*D))
        WRITE(NFIL,FMT='(F10.5,10F15.5)')D,FINT
      ENDDO
      CALL FILEHANDLER$CLOSE('TEST')
      WRITE(NFIL,FMT='(80("+"),T10," STOPPING AFTER TEST ")')
      STOP 'SIMPLELMTO_TESTOFFSITEXCONVERT_B'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_PRBONDU()
!     **************************************************************************
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1 &
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
      LMNXA=SUM(2*POTPAR(ISP1)%LOXH+1)
      LMNXB=SUM(2*POTPAR(ISP2)%LOXH+1)
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
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1
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
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1 &
     &                             ,SCREENL
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
!     -- LAMBDA MUST BE CONSISTENT WITH SIMPLELMTO_OFFSITEX22SETUP
      REAL(8)   ,PARAMETER  :: LAMBDA=-0.5D0 
      REAL(8)               :: EXPLMDADIS
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
      REAL(8)               :: FAC
      REAL(8)               :: INV(3),DINV(3)
!     **************************************************************************
!
!     ==========================================================================
!     == FUNCTIONS CAPTURING THE LONG-RANGE BEHAVIOR OF THE COULOMB INTERACTION=
!     ==========================================================================
      CALL SIMPLELMTO_EXPLDISF(SCREENL,DIS,INV,DINV)
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
!!$CALL SIMPLELMTO_OFFSITEXVALUE('22',ISP1,ISP2,IND,DIS+1.D-3,INTEGRAL,DINTEGRAL)
!!$PRINT*,'DIS+DELTA ',INTEGRAL,DINTEGRAL
!!$CALL SIMPLELMTO_OFFSITEXVALUE('22',ISP1,ISP2,IND,DIS-1.D-3,INTEGRAL,DINTEGRAL)
!!$PRINT*,'DIS-DELTA ',INTEGRAL,DINTEGRAL
!!$STOP

!
!                   == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE...
                    CALL SIMPLELMTO_OFFSITEXVALUE('22',ISP1,ISP2,IND,DIS &
     &                                          ,INTEGRAL,DINTEGRAL)
!PRINT*,'22 FROMGRID ',DIS,INTEGRAL,DINTEGRAL
!                   ============================================================
!                   == ADD IN LONG RANGE PART FROM MONO- AND DIPOLES THAT HAS ==
!                   == BEEN SUBTRACTED IN SIMPLELMTO_OFFSITEX22SETUP          ==
!                   ============================================================
                    IF(LR1.LE.1.AND.LR2.LE.1) THEN
                      FAC=POTPAR(ISP1)%QLN(LR1+1,LN1,LN2) &
     &                   *POTPAR(ISP2)%QLN(LR2+1,LN3,LN4)
                      IF(LR1+LR2.EQ.0) THEN ! MONOPOLE-MONOPOLE
                        INTEGRAL =INTEGRAL +FAC*INV(1)
                        DINTEGRAL=DINTEGRAL+FAC*DINV(1)
!PRINT*,'22 MONOPOLE ',DIS,INTEGRAL,DINTEGRAL
                      ELSE IF(LR1+LR2.EQ.1) THEN !MONOPOLE-DIPOLE
                        SVAR=FAC
                        IF(LR1.EQ.0) SVAR=-SVAR
                        INTEGRAL =INTEGRAL +SVAR*INV(2)
                        DINTEGRAL=DINTEGRAL+SVAR*DINV(2)
                      ELSE IF(LR1+LR2.EQ.2) THEN !DIPOLE-DIPOLE
                        INTEGRAL =INTEGRAL +FAC*INV(3)
                        DINTEGRAL=DINTEGRAL+FAC*DINV(3)
                        IF(MABS.EQ.0) THEN
                          INTEGRAL =INTEGRAL -3.D0*FAC*INV(3)
                          DINTEGRAL=DINTEGRAL-3.D0*FAC*DINV(3)
                        END IF
                      END IF
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
          DU(LMN2,LMN1,:,:)=DU(LMN1,LMN2,:,:)
        ENDDO
      ENDDO
      DO LMN3=1,LMNX2
        DO LMN4=LMN3+1,LMNX2
          U(:,:,LMN4,LMN3)=U(:,:,LMN3,LMN4)
          DU(:,:,LMN4,LMN3)=DU(:,:,LMN3,LMN4)
        ENDDO
      ENDDO
!PRINT*,'22 U        ',DIS,U,DU
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEX31U(ISP1,ISP2,DIS,LMNX1,LMNX2,U,DU)
!     **************************************************************************
!     ** U-TENSOR FOR TWO ATOMS 
!     ** THREE ORBITALS ARE ON THE FIRST AND ONE ORBITALS ON THE SECOND ATOM  **
!     **                                                                      **
!     ** U(I,J,K,L)= INT CHI_I(R1)*CHI_J(R1)*CHI_K(R2)*CHIL(R2)/|R1-R2|       **
!     **    WITH I,J,K ON ATOM A  AND L ON B                                  **
!     **************************************************************************
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1
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
CHARACTER(64) :: STRING
INTEGER(4)    :: IDIS
INTEGER(8)    :: INTBIG
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

!!$INTBIG=LN1+10*(LN2+10*(LN3+10*(LN4+10*(LR1+10*(LR2+10*(MABS+10))))))
!!$WRITE(STRING,*)INTBIG
!!$STRING='Z'//TRIM(ADJUSTL(STRING))//'.DAT'
!!$PRINT*,'STRING ',IND,LN1,LN2,LN3,LN4,LR1,LR2,MABS,INTBIG,TRIM(STRING),HUGE(INTBIG)
!!$OPEN(1237,FILE=TRIM(STRING))
!!$DO IDIS=1,1000
!!$  SVAR=REAL(IDIS)*1.D-2
!!$  CALL SIMPLELMTO_OFFSITEXVALUE('31',ISP1,ISP2,IND,SVAR,INTEGRAL,DINTEGRAL)
!!$  WRITE(1237,FMT='(3F20.10)')SVAR,INTEGRAL,DINTEGRAL
!!$ENDDO
!!$CLOSE(1237)
!
!
!                   == EVALUATE MATRIX ELEMENT FOR THE DISTANCE HERE...
                    CALL SIMPLELMTO_OFFSITEXVALUE('31',ISP1,ISP2,IND,ABS(DIS) &
     &                                          ,INTEGRAL,DINTEGRAL)

!CALL SIMPLELMTO_OFFSITEXVALUE_TEST('31',ISP1,ISP2,IND,ABS(DIS))
                    IF(DIS.LT.0.D0) THEN
                      SVAR=(-1.D0)**(L1+L2+L3+L4)
                      INTEGRAL=SVAR*INTEGRAL
                      DINTEGRAL=-SVAR*DINTEGRAL
CALL ERROR$MSG('DIS<0 SHOULD NOT HAPPEN')
CALL ERROR$STOP('SIMPLELMTO_OFFSITEX31U')
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
!VERSION 1 IS PROBABLY FASTER!
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
      DO LN1=1,LNX1
        L1=POTPAR(ISP1)%LOXH(LN1)
        DO LN2=LN1,LNX1
          L2=POTPAR(ISP1)%LOXH(LN2)
          LMN1=LMN0A(LN1)
          DO M1=1,2*L1+1
            LMN1=LMN1+1
            LMN2=LMN0A(LN2)
            DO M2=1,2*L2+1
              LMN2=LMN2+1
              U(LMN2,LMN1,:,:) =U(LMN1,LMN2,:,:)
              DU(LMN2,LMN1,:,:)=DU(LMN1,LMN2,:,:)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!!$      DO LMN1=1,LMNX1
!!$        DO LMN2=LMN1+1,LMNX1
!!$          U(LMN2,LMN1,:,:) =U(LMN1,LMN2,:,:)
!!$          DU(LMN2,LMN1,:,:)=DU(LMN1,LMN2,:,:)
!!$        ENDDO
!!$      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEXVALUE_TEST(ID,ISP1,ISP2,IND,DIS)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4),INTENT(IN) :: ISP1
      INTEGER(4),INTENT(IN) :: ISP2
      INTEGER(4),INTENT(IN) :: IND
      REAL(8)   ,INTENT(IN) :: DIS
      REAL(8)   ,PARAMETER  :: DELTA=1.D-3
      REAL(8)               :: X1,X2
      REAL(8)               :: Y1,Y2
      REAL(8)               :: DY1,DY2
      REAL(8)               :: DYA,DYB
      REAL(8)               :: S
!     **************************************************************************
      X1=DIS+DELTA
      X2=DIS-DELTA
      CALL SIMPLELMTO_OFFSITEXVALUE(ID,ISP1,ISP2,IND,X1,Y1,DY1)
      CALL SIMPLELMTO_OFFSITEXVALUE(ID,ISP1,ISP2,IND,X2,Y2,DY2)
      DYA=(Y2-Y1)/(X2-X1)
      DYB=(DY1+DY2)/2.D0
      S=(DYA-DY1)/(DY2-DY1)
      IF((S.LT.-1.D-1.OR.S.GT.1.D0+1.D-1).AND.ABS((DYA-DYB)/DYB).GT.1.D-5) THEN
        CALL ERROR$MSG('DERIVATIVE TEST FAILED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$I4VAL('ISP1',ISP1)
        CALL ERROR$I4VAL('ISP2',ISP2)
        CALL ERROR$I4VAL('IND',IND)
        CALL ERROR$R8VAL('DIS ',DIS)
        CALL ERROR$R8VAL('S ',S)
        CALL ERROR$R8VAL('DY/DX(X1) ',DY1)
        CALL ERROR$R8VAL('DY/DX(NUM)',DYA)
        CALL ERROR$R8VAL('DY/DX(X2) ',DY2)
        CALL ERROR$R8VAL('(DYA-DYB)/DYB',DYA/DYB-1.D0)
        CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE_TEST31')
      END IF
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
      USE SIMPLELMTO_MODULE, ONLY : POTPAR=>POTPAR1
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
      REAL(8)     ,SAVE       :: DG(LF)
!     **************************************************************************
      INTEGRAL=0.D0
      IF(ISP1.NE.ISP1SAVE.OR.ISP2.NE.ISP2SAVE.OR.DIS.NE.DISSAVE) THEN
        ISP1SAVE=ISP1
        ISP2SAVE=ISP2
        IF(.NOT.ALLOCATED(OFFSITEX)) THEN
          CALL ERROR$MSG('ARRAY OFFSITEX IN SIMPLELMTO_MODULE')
          CALL ERROR$MSG('HAS NOT BEEN ALLOCATED')
          CALL ERROR$MSG('INTERNAL CODING ERROR')
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        NF=OFFSITEX(ISP1,ISP2)%NF
        IF(NF.GT.LF) THEN
          CALL ERROR$MSG('INTERNAL ARRAY SIZE EXCEEDED')
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        DISSAVE=DIS      
!
!       __INTERPOLATION FUNCTIONS: (1+LAMBDA*R)*EXP(-LAMBDA*R)__________________
!       __EXPONENTIAL TAIL AND VANISHING SLOPE AT THE ORIGIN____________________
        G(:NF)=EXP(-OFFSITEX(ISP1,ISP2)%LAMBDA(:)*DIS) &
     &        *(1.D0+OFFSITEX(ISP1,ISP2)%LAMBDA(:)*DIS) 
        DG(:NF)=-OFFSITEX(ISP1,ISP2)%LAMBDA(:)**2*DIS &
     &          /(1.D0+OFFSITEX(ISP1,ISP2)%LAMBDA(:)*DIS)*G(:NF)
      END IF
!
      IF(ID.EQ.'22') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%X22)) THEN
          CALL ERROR$MSG('NDDO TERM REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL =SUM(OFFSITEX(ISP1,ISP2)%X22(:NF,IND)*G(:NF))
        DINTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%X22(:NF,IND)*DG(:NF))
!
      ELSE IF (ID.EQ.'31') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%X31)) THEN
          CALL ERROR$MSG('3-1 EXCHANGE TERM REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%X31(:NF,IND)*G(:NF))
        DINTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%X31(:NF,IND)*DG(:NF))
!
      ELSE IF (ID.EQ.'BONDU') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%BONDU)) THEN
          CALL ERROR$MSG('U-TERM OF BONDOVERLAPS REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%BONDU(:NF,IND)*G(:NF))
        DINTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%BONDU(:NF,IND)*DG(:NF))
!
      ELSE IF (ID.EQ.'OV') THEN
        IF(.NOT.ASSOCIATED(OFFSITEX(ISP1,ISP2)%OVERLAP)) THEN
          CALL ERROR$MSG('OFFSITE OVERLAP REQUESTED BUT NOT INITIALIZED')
          CALL ERROR$I4VAL('ISP1',ISP1)
          CALL ERROR$I4VAL('ISP2',ISP2)
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
        END IF
        INTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%OVERLAP(:NF,IND)*G(:NF))
        DINTEGRAL=SUM(OFFSITEX(ISP1,ISP2)%OVERLAP(:NF,IND)*DG(:NF))
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('SIMPLELMTO_OFFSITEXVALUE')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SIMPLELMTO_OFFSITEXEVAL(TPARALLEL,NAT,NND,RBAS,R0,DENMAT &
     &                                  ,EX,STRESS,FORCE,HAMIL,OVERLAP)
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
     &                            ,POTPAR=>POTPAR1 &
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
      TYPE(RSPACEMAT_TYPE),INTENT(IN)   :: DENMAT(NND) !INOUT FOR TESTING
      REAL(8)             ,INTENT(OUT)  :: EX
      REAL(8)             ,INTENT(OUT)  :: STRESS(3,3)
      REAL(8)             ,INTENT(OUT)  :: FORCE(3,NAT)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: HAMIL(NND)
      TYPE(RSPACEMAT_TYPE),INTENT(INOUT):: OVERLAP(NND)
      LOGICAL(4)          ,PARAMETER    :: TPR=.FALSE.
      LOGICAL(4)          ,PARAMETER    :: TOV=.TRUE.
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
      REAL(8)                :: ROT(3,3)     ! ROTATION MATRIX
      REAL(8)                :: DROT(3,3,3)  ! D-ROT(I,J)/D-DR(K)
      REAL(8)   ,ALLOCATABLE :: YLMROT(:,:)
      REAL(8)   ,ALLOCATABLE :: DYLMROT(:,:,:)
      REAL(8)   ,ALLOCATABLE :: UROTA(:,:)
      REAL(8)   ,ALLOCATABLE :: UROTB(:,:)
      REAL(8)   ,ALLOCATABLE :: DUROTA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DUROTB(:,:,:)
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
      REAL(8)                :: SQLHFWA,SQLHFWB
      REAL(8)                :: SCALE
      REAL(8)                :: DEDD
      REAL(8)                :: DEDR(3)
      REAL(8)                :: QOFFSITE
      INTEGER(4)             :: LMRX
      INTEGER(4)             :: LRX
      INTEGER(4)             :: LX
      LOGICAL(4)             :: TNDDO,T31,TBONDX
      INTEGER(4)            :: THISTASK,NTASKS
INTEGER(4) :: J,K,I1,J1,K1
REAL(8) :: SVAR1
INTEGER(4) ,PARAMETER :: NX=3
INTEGER(4) :: IX,IND1,IND2,IND3,NNTEST
REAL(8)    :: XDELTA,XSVAR,TESTX(-NX:NX),TESTY(-NX:NX),TESTDYDX(-NX:NX)
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
      QOFFSITE=0.D0
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
!!$NNTEST=-1
!!$IF(NN.EQ.NNTEST) THEN
!!$  IND1=2
!!$  IND2=1
!!$  IND3=1
!!$  OPEN(1753,FILE='TESTSET')
!!$  READ(1753,*)IND1,IND2,IND3
!!$  WRITE(*,*)'TESTSET',IND1,IND2,IND3,' NNTEST=',NN,'NDIMD=',DENMAT(NN)%N3
!!$  CLOSE(1753)
!!$  XSVAR=DENMAT(NN)%MAT(IND1,IND2,IND3)
!!$  XDELTA=1.D-2
!!$  IX=-NX-1
!!$  1000 CONTINUE
!!$  IX=IX+1
!!$  TESTX(IX)=XDELTA*REAL(IX,KIND=8)
!!$  DENMAT(NN)%MAT(IND1,IND2,IND3)=XSVAR+TESTX(IX)
!!$  HAMIL(NN)%MAT(IND1,IND2,IND3)=0.D0
!!$  EX=0.D0
!!$END IF
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
        SQLHFWA=HYBRIDSETTING(ISPA)%LHFWEIGHT**0.25D0
        SQLHFWB=HYBRIDSETTING(ISPB)%LHFWEIGHT**0.25D0
!       == DO NOTHING UNLESS NEEDED
        IF(.NOT.(TNDDO.OR.T31.OR.TBONDX)) CYCLE
        IF(SQLHFWA*SQLHFWB.EQ.0.D0) CYCLE

CALL TIMING$CLOCKON('LOOP:OFFX')
!
!       ========================================================================
!       == BLOW UP DENSITY MATRIX                                             ==
!       ========================================================================
CALL TIMING$CLOCKON('OFFX:BLOWUP')
        LMNXA=SUM(2*POTPAR(ISPA)%LOXH+1)
        LMNXB=SUM(2*POTPAR(ISPB)%LOXH+1)
        NDIMD=DENMAT(NN)%N3
        IF(LMNXA.NE.DENMAT(NN)%N1.OR.LMNXB.NE.DENMAT(NN)%N2) THEN
          CALL ERROR$MSG('LMNX NOT EQUAL LMNXT') 
          CALL ERROR$STOP('SIMPLELMTO_OFFSITEXEVAL')
        END IF
        ALLOCATE(DAB(LMNXA,LMNXB,NDIMD))
        ALLOCATE(DBA(LMNXB,LMNXA,NDIMD))
        ALLOCATE(DAA(LMNXA,LMNXA,NDIMD))
        ALLOCATE(DBB(LMNXB,LMNXB,NDIMD))
        DAB(:,:,:)=DENMAT(NN )%MAT(:,:,:)
        DBA(:,:,:)=DENMAT(INDM(NN))%MAT(:,:,:)
        DAA(:,:,:)=DENMAT(NNA)%MAT(:,:,:)
        DBB(:,:,:)=DENMAT(NNB)%MAT(:,:,:)
CALL TIMING$CLOCKOFF('OFFX:BLOWUP')
 !
!       ========================================================================
!       == ROTATE DENSITY MATRIX SO THAT DISTANCE VECTOR POINTS IN Z-DIRECTION=
!       ========================================================================
CALL TIMING$CLOCKON('OFFX:ROTATE1')
!
!       == ROTATION MATRIX =====================================================
        CALL SPHERICAL$AXISUROT(DR,ROT,DROT)
!
!       == CONSTRUCT TRANSFORMATION MATRIX FOR REAL SPHERICAL HARMONICS ========
        LMX=MAX(MAXVAL(POTPAR(ISPA)%LOXH),MAXVAL(POTPAR(ISPB)%LOXH))
        LMX=(LMX+1)**2
        ALLOCATE(YLMROT(LMX,LMX))
        ALLOCATE(DYLMROT(LMX,LMX,3))
        CALL SPHERICAL$ROTATEYLMWDER(LMX,3,ROT,DROT,YLMROT,DYLMROT)
!
!       == CONSTRUCT ROTATION MATRIX OF TAILED ORBITALS ========================
        ALLOCATE(UROTA(LMNXA,LMNXA))
        ALLOCATE(UROTB(LMNXB,LMNXB))
        ALLOCATE(DUROTA(LMNXA,LMNXA,3))
        ALLOCATE(DUROTB(LMNXB,LMNXB,3))
        LNXA=POTPAR(ISPA)%LNXH
        UROTA(:,:)=0.D0
        DUROTA(:,:,:)=0.D0
        LMN=1
        DO LN=1,LNXA
          L=POTPAR(ISPA)%LOXH(LN)
          LM=L**2+1
          UROTA(LMN:LMN+2*L,LMN:LMN+2*L)   = YLMROT(LM:LM+2*L,LM:LM+2*L)
          DUROTA(LMN:LMN+2*L,LMN:LMN+2*L,:)=DYLMROT(LM:LM+2*L,LM:LM+2*L,:)
          LMN=LMN+2*L+1
        ENDDO
        LNXB=POTPAR(ISPB)%LNXH
        UROTB(:,:)=0.D0
        DUROTB(:,:,:)=0.D0
        LMN=1
        DO LN=1,LNXB
          L=POTPAR(ISPB)%LOXH(LN)
          LM=L**2+1
          UROTB(LMN:LMN+2*L,LMN:LMN+2*L)   = YLMROT(LM:LM+2*L,LM:LM+2*L)
          DUROTB(LMN:LMN+2*L,LMN:LMN+2*L,:)=DYLMROT(LM:LM+2*L,LM:LM+2*L,:)
          LMN=LMN+2*L+1
        ENDDO
        DEALLOCATE(YLMROT)
        DEALLOCATE(DYLMROT)
!
!       == ROTATE DENSITY MATRIX ===============================================
        !SUGGESTION: SPEED UP BY EXPLOITING THAT UROT IS SPARSE
        DO I=1,NDIMD
          DAA(:,:,I)=MATMUL(UROTA,MATMUL(DAA(:,:,I),TRANSPOSE(UROTA)))
          DAB(:,:,I)=MATMUL(UROTA,MATMUL(DAB(:,:,I),TRANSPOSE(UROTB)))
          DBA(:,:,I)=MATMUL(UROTB,MATMUL(DBA(:,:,I),TRANSPOSE(UROTA)))
          DBB(:,:,I)=MATMUL(UROTB,MATMUL(DBB(:,:,I),TRANSPOSE(UROTB)))
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
          DO LMN1A=1,LMNXA
            DO LMN1B=1,LMNXB
              QOFFSITE=QOFFSITE+DBA(LMN1B,LMN1A,1)*OV(LMN1A,LMN1B)
            ENDDO
          ENDDO
!         __OV WILL BE DEALLOCATED BELOW________________________________________
          DEALLOCATE(DOV)
        END IF
!
!       ========================================================================
!       == NDDO TERM:                                                         ==
!       == U22(1,2,3,4)=INT DX1 INT DX2: A1(X1)*A2(X1)][B3(X2)*B4(X2)]/|R1-R2|==
!       == LMN1A,LMN2A TIED TO COORDINATE X1, LMN1B,LMN2B TO X2               ==
!       ==                                                                    ==
!       == CAUTION: H(I,J)=DE/DRHO(I,J), WHICH IS HAMILTONIAN^DAGGE           ==
!       ========================================================================
        IF(TNDDO) THEN
          CALL TIMING$CLOCKON('OFFX:NDDO')
          ALLOCATE(U22   (LMNXA,LMNXA,LMNXB,LMNXB))
          ALLOCATE(DU22  (LMNXA,LMNXA,LMNXB,LMNXB))

!!$CALL SIMPLELMTO_OFFSITEX22U(ISPA,ISPB,DIS-1.D-3,LMNXA,LMNXB,U22,DU22)
!!$PRINT*,'U22-DELTA ',U22,DU22
!!$CALL SIMPLELMTO_OFFSITEX22U(ISPA,ISPB,DIS+1.D-3,LMNXA,LMNXB,U22,DU22)
!!$PRINT*,'U22+DELTA ',U22,DU22
!!$STOP 'FORCED'

          CALL SIMPLELMTO_OFFSITEX22U(ISPA,ISPB, DIS,LMNXA,LMNXB,U22,DU22)
          SCALE=(SQLHFWA*SQLHFWB)**2
          DO LMN1B=1,LMNXB
            DO LMN2B=1,LMNXB
              DO LMN2A=1,LMNXA
                DO LMN1A=1,LMNXA
!                 ==  W(LMN1A,LMN1B,LMN2A,LMN2B) ==============================
                  SVAR =-0.25D0*SCALE* U22(LMN1A,LMN2A,LMN2B,LMN1B)
                  DSVAR=-0.25D0*SCALE*DU22(LMN1A,LMN2A,LMN2B,LMN1B)
                  EX  =EX  + SVAR*SUM(DAB(LMN1A,LMN1B,:)*DBA(LMN2B,LMN2A,:))
                  DEDD=DEDD+DSVAR*SUM(DAB(LMN1A,LMN1B,:)*DBA(LMN2B,LMN2A,:))
                  HAB(LMN1A,LMN1B,:)=HAB(LMN1A,LMN1B,:)+SVAR*DBA(LMN2B,LMN2A,:)
                  HBA(LMN2B,LMN2A,:)=HBA(LMN2B,LMN2A,:)+SVAR*DAB(LMN1A,LMN1B,:)
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
!       ==                                                                    ==
!       == CAUTION: H(I,J)=DE/DRHO(I,J), WHICH IS HAMILTONIAN^DAGGE           ==
!       ========================================================================
        IF(T31) THEN
          CALL TIMING$CLOCKON('OFFX:31')
          ALLOCATE(U3A1B (LMNXA,LMNXA,LMNXA,LMNXB))
          ALLOCATE(DU3A1B(LMNXA,LMNXA,LMNXA,LMNXB))

!!$CALL SIMPLELMTO_OFFSITEX31U(ISPA,ISPB,DIS-1.D-3,LMNXA,LMNXB,U3A1B,DU3A1B)
!!$PRINT*,U3A1B,DU3A1B
!!$CALL SIMPLELMTO_OFFSITEX31U(ISPA,ISPB,DIS+1.D-3,LMNXA,LMNXB,U3A1B,DU3A1B)
!!$PRINT*,U3A1B,DU3A1B
!!$STOP 'FORCED'

          CALL SIMPLELMTO_OFFSITEX31U(ISPA,ISPB, DIS,LMNXA,LMNXB,U3A1B,DU3A1B)
          SCALE=SQLHFWA**3*SQLHFWB
          DO LMN1B=1,LMNXB
            DO LMN3A=1,LMNXA 
              DO LMN2A=1,LMNXA
                DO LMN1A=1,LMNXA
!THE HARTREE TERM CONNECTAS DAA(LMN2A,LMN1A)* DBA(LMN1B,LMN3A)
!THE EXCHANGE COMNECTS     -DAA(LMN2A,LMN3A)* DBA(LMN1B,LMN1A)
!PRINT*,LMN1A,LMN2A,LMN3A,LMN1B,U3A1B(LMN1A,LMN2A,LMN3A,LMN1B),DU3A1B(LMN1A,LMN2A,LMN3A,LMN1B)

                  SVAR =-0.25D0*SCALE* U3A1B(LMN1A,LMN2A,LMN3A,LMN1B)
                  DSVAR=-0.25D0*SCALE*DU3A1B(LMN1A,LMN2A,LMN3A,LMN1B)
!
                  EX  =EX   +SVAR*SUM((DAB(LMN1A,LMN1B,:)+DBA(LMN1B,LMN1A,:)) &
         &                           *(DAA(LMN2A,LMN3A,:)+DAA(LMN3A,LMN2A,:)))
                  DEDD=DEDD+DSVAR*SUM((DAB(LMN1A,LMN1B,:)+DBA(LMN1B,LMN1A,:)) &
         &                           *(DAA(LMN2A,LMN3A,:)+DAA(LMN3A,LMN2A,:)))
                  HAB(LMN1A,LMN1B,:)=HAB(LMN1A,LMN1B,:) &
         &                      +SVAR*(DAA(LMN2A,LMN3A,:)+DAA(LMN3A,LMN2A,:))
                  HBA(LMN1B,LMN1A,:)=HBA(LMN1B,LMN1A,:) &
         &                      +SVAR*(DAA(LMN2A,LMN3A,:)+DAA(LMN3A,LMN2A,:))
                  HAA(LMN2A,LMN3A,:)=HAA(LMN2A,LMN3A,:) &
         &                      +SVAR*(DAB(LMN1A,LMN1B,:)+DBA(LMN1B,LMN1A,:))
                  HAA(LMN3A,LMN2A,:)=HAA(LMN3A,LMN2A,:) &
         &                      +SVAR*(DAB(LMN1A,LMN1B,:)+DBA(LMN1B,LMN1A,:))
                ENDDO
              ENDDO
            ENDDO
          ENDDO               
          DEALLOCATE(U3A1B)
          DEALLOCATE(DU3A1B)
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
          SCALE=(SQLHFWA*SQLHFWB)**2
          DO LMN2B=1,LMNXB
            DO LMN2A=1,LMNXA 
              DO LMN1B=1,LMNXB
                DO LMN1A=1,LMNXA
                  SVAR =-0.25D0*SCALE* BONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                  DSVAR=-0.25D0*SCALE*DBONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                  EX  =EX  + SVAR*SUM(DAA(LMN1A,LMN2A,:)*DBB(LMN1B,LMN2B,:))
                  DEDD=DEDD+DSVAR*SUM(DAA(LMN1A,LMN2A,:)*DBB(LMN1B,LMN2B,:)) 
                  HAA(LMN1A,LMN2A,:)=HAA(LMN1A,LMN2A,:)+SVAR*DBB(LMN1B,LMN2B,:)
                  HBB(LMN1B,LMN2B,:)=HBB(LMN1B,LMN2B,:)+SVAR*DAA(LMN1A,LMN2A,:)
                  SVAR =-0.25D0*SCALE* BONDU(LMN1A,LMN1B,LMN2A,LMN2B)
                  DSVAR=-0.25D0*SCALE*DBONDU(LMN1A,LMN1B,LMN2A,LMN2B)
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
!       == ACCUMULATE FORCES AND STRESSES                                     ==
!       ========================================================================
!       == DAA AND DBB WILL NO MORE BE USED AS DENSITY MATRIX ==================
!       == AND WILL BE USED AS CONTAINER =======================================
        DO I=1,NDIMD
          DAA(:,:,I)=MATMUL(TRANSPOSE(DAA(:,:,I)),HAA(:,:,I)) &
       &            +MATMUL(TRANSPOSE(DBA(:,:,I)),HBA(:,:,I)) &
       &            +MATMUL(DAA(:,:,I),TRANSPOSE(HAA(:,:,I))) &
       &            +MATMUL(DAB(:,:,I),TRANSPOSE(HAB(:,:,I))) 

          DBB(:,:,I)=MATMUL(TRANSPOSE(DBB(:,:,I)),HBB(:,:,I)) &
       &            +MATMUL(TRANSPOSE(DAB(:,:,I)),HAB(:,:,I)) &
       &            +MATMUL(DBB(:,:,I),TRANSPOSE(HBB(:,:,I))) &
       &            +MATMUL(DBA(:,:,I),TRANSPOSE(HBA(:,:,I)))
        ENDDO
!       == SUM UP TOTAL AND SPIN CONTRIBUTIONS =================================
        DO I=2,NDIMD
          DAA(:,:,1)=DAA(:,:,1)+DAA(:,:,I)
          DBB(:,:,1)=DBB(:,:,1)+DBB(:,:,I)
        ENDDO
!       == DEDR IS THE ENERGY DERIVATIVE WITH RESPECT TO THE AXIAL FORCE========
        DEDR(:)=DEDD*DR(:)/DIS
!       == ADD TRANSVERSAL FORCE TO DEDR =======================================
        DO J=1,3
          DUROTA(:,:,J)=MATMUL(DUROTA(:,:,J),TRANSPOSE(UROTA))
          DUROTB(:,:,J)=MATMUL(DUROTB(:,:,J),TRANSPOSE(UROTB))
!!$          DEDR(J)=DEDR(J)+SUM(DAA(:,:,1)*TRANSPOSE(DUROTA(:,:,J))) &
!!$       &                 +SUM(DBB(:,:,1)*TRANSPOSE(DUROTB(:,:,J))) 
          DEDR(J)=DEDR(J)+1.D0*( SUM(DAA(:,:,1)*TRANSPOSE(DUROTA(:,:,J))) &
       &                        +SUM(DBB(:,:,1)*TRANSPOSE(DUROTB(:,:,J))) )
        ENDDO
!       == EXTRACT FORCES AND STRESSES FROM DEDR ===============================
        FORCE(:,IATB)=FORCE(:,IATB)-DEDR(:)
        FORCE(:,IATA)=FORCE(:,IATA)+DEDR(:)
        DO J=1,3
          DEDRBAS(:,J)=DEDRBAS(:,J)+REAL(DENMAT(NN)%IT(:))*DEDR(J)
        ENDDO
!WRITE(*,FMT='("DEDR    ",I5,3F20.10)')NN,DEDR
! 
!
!       ========================================================================
!       == ROTATE HAMILTONIAN BACK                                            ==
!       ========================================================================
        DO I=1,NDIMD
          HAA(:,:,I)=MATMUL(TRANSPOSE(UROTA),MATMUL(HAA(:,:,I),UROTA))
          HAB(:,:,I)=MATMUL(TRANSPOSE(UROTA),MATMUL(HAB(:,:,I),UROTB))
          HBA(:,:,I)=MATMUL(TRANSPOSE(UROTB),MATMUL(HBA(:,:,I),UROTA))
          HBB(:,:,I)=MATMUL(TRANSPOSE(UROTB),MATMUL(HBB(:,:,I),UROTB))
        ENDDO
        IF(TOV) OV(:,:)=MATMUL(TRANSPOSE(UROTA),MATMUL(OV(:,:),UROTB))
        DEALLOCATE(UROTA)
        DEALLOCATE(UROTB)
        DEALLOCATE(DUROTA)
        DEALLOCATE(DUROTB)
!
!       ========================================================================
!       == MAP HAMILTONIAN BACK                                               ==
!       ========================================================================
        HAMIL(NN )%MAT(:,:,:)     =HAMIL(NN )%MAT(:,:,:)+HAB(:,:,:)
        HAMIL(INDM(NN))%MAT(:,:,:)=HAMIL(INDM(NN))%MAT(:,:,:)+HBA(:,:,:)
        HAMIL(NNA)%MAT(:,:,:)     =HAMIL(NNA)%MAT(:,:,:)+HAA(:,:,:)
        HAMIL(NNB)%MAT(:,:,:)     =HAMIL(NNB)%MAT(:,:,:)+HBB(:,:,:)
        DEALLOCATE(DAB)
        DEALLOCATE(DBA)
        DEALLOCATE(HAB)
        DEALLOCATE(HBA)
        DEALLOCATE(DAA)
        DEALLOCATE(HAA)
        DEALLOCATE(DBB)
        DEALLOCATE(HBB)
        IF(TOV) THEN
          OVERLAP(NN)%MAT(:,:,1)=OVERLAP(NN )%MAT(:,:,1)+OV(:,:)
          DEALLOCATE(OV)
        END IF
CALL TIMING$CLOCKOFF('LOOP:OFFX')
!!$IF(NN.EQ.NNTEST) THEN
!!$   TESTY(IX)=EX
!!$   TESTDYDX(IX)=HAMIL(NN)%MAT(IND1,IND2,IND3)
!!$   WRITE(*,*)XDELTA*REAL(IX,8),EX,HAMIL(NN)%MAT(IND1,IND2,IND3)
!!$   IF(IX.EQ.NX) THEN
!!$     WRITE(*,*)0,TESTDYDX(0),TESTDYDX(0)
!!$     DO I=1,NX
!!$       WRITE(*,*)I,(TESTY(I)-TESTY(-I))/(TESTX(I)-TESTX(-I)) &
!!$&                ,0.5D0*(TESTDYDX(I)+TESTDYDX(-I))
!!$     ENDDO
!!$     STOP 'FORCED'
!!$   END IF
!!$   GOTO 1000
!!$END IF
      ENDDO
!
!     ==========================================================================
!     == CONVERT DERIVATIVE INTO STRESSES                                     ==
!     ==========================================================================
!          M_T DDOT(T)=-DE/DT
      STRESS=MATMUL(RBAS,DEDRBAS)

!     STRESS=DEDRBAS !THIS IS FOR TESTING THE NUMERICAL DERIVATIVES
!
!     ==========================================================================
!     == CLEAN UP TO AVOID MEMORY LEAK                                        ==
!     ==========================================================================
      DEALLOCATE(NNAT)
!
!     ==========================================================================
!     == DIAGNOSTIC PRINTOUT                                                  ==
!     ==========================================================================
      WRITE(*,*)'QOFFSITE=',QOFFSITE
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
     &                            ,POTPAR=>POTPAR1 
      IMPLICIT NONE
      INTEGER(4),PARAMETER   :: ND=1000
      REAL(8)   ,PARAMETER   :: DMIN=0.2D0
      REAL(8)   ,PARAMETER   :: DMAX=5.D0
      REAL(8)                :: DIS
      REAL(8)   ,ALLOCATABLE :: U22(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: DU22(:,:,:,:)
      INTEGER(4)             :: ISPA,ISPB,ID
      INTEGER(4)             :: LMNXA,LMNXB
      INTEGER(4)             :: NFIL=8
      INTEGER(4)             :: LMN1=1,LMN2=3,LMN3=3,LMN4=1
!     **************************************************************************
      OPEN(NFIL,FILE='TESTOFFSITEX.DAT')
      REWIND (NFIL)
      DO ISPA=1,NSP
        DO ISPB=ISPA,NSP
          LMNXA=SUM(2*POTPAR(ISPA)%LOXH+1)
          LMNXB=SUM(2*POTPAR(ISPB)%LOXH+1)
          ALLOCATE(U22   (LMNXA,LMNXA,LMNXB,LMNXB))
          ALLOCATE(DU22  (LMNXA,LMNXA,LMNXB,LMNXB))
          U22=0.D0
          DU22=0.D0
  WRITE(NFIL,*)DMIN,U22(LMN1,LMN2,LMN3,LMN4),DU22(LMN1,LMN2,LMN3,LMN4)
!          WRITE(NFIL,*)DMIN,U22,DU22
          DO ID=1,ND
            DIS=DMIN+REAL(ID-1,KIND=8)/REAL(ND-1,KIND=8)*(DMAX-DMIN)
            CALL SIMPLELMTO_OFFSITEX22U(ISPA,ISPB,DIS,LMNXA,LMNXB,U22,DU22)
!            WRITE(NFIL,*)DIS,U22,DU22
     WRITE(NFIL,*)DIS,U22(LMN1,LMN2,LMN3,LMN4),DU22(LMN1,LMN2,LMN3,LMN4)
          ENDDO
          U22=0.D0
          DU22=0.D0
          WRITE(NFIL,*)DIS,U22(LMN1,LMN2,LMN3,LMN4),DU22(LMN1,LMN2,LMN3,LMN4)
!          WRITE(NFIL,*)DMAX,U22,DU22
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
!       == EQ. 19 OF ROMANOWSKI08_IJQC108_249 ==================================
        SVAR=3.D0*SQRT(3.D0)/(16.D0*DIS)*EXP(-3.D0*DIS) &
     &      *(3.D0+2.D0*DIS-EXP(2.D0*DIS)*(3.D0-6.D0*DIS))
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
!       == EQ. 21 OF ROMANOWSKI08_IJQC108_249 ==================================
        SVAR=1.D0/315.D0*EXP(-4.D0*DIS)*(315.D0+4.D0*DIS*(315.D0+4.D0*DIS &
    &       *(75.D0+8.D0*DIS*(-15.D0+4.D0*DIS*(-9.D0-2.D0*DIS+8.D0*DIS**2)))))
        PRINT*,'DIS=',DIS,'OVERLAP ',OVERLAP,SVAR,OVERLAP-SVAR
      ENDDO
!     ==========================================================================
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
!      ** USES THE  TWO-DIMENSIONAL ADAPTIVE INTEGRATION PROPOSED BY          **
!      ** GENZ AND MALIK, WHICH IS IMPLEMENTED IN THE NEWADAPT MODULE         **
!      ** NEWADAPT USES THE FUNCTION CALL SIMPLELMTO_TWOCENTER_MYFUNC(P,VALUE)**
!      ** WHICH IS SUPPLIED HERE.                                             **
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
       CALL RADIAL$VALUE(GID1,NR1,F1,ABS(R1),VAL1)
       CALL RADIAL$VALUE(GID2,NR2,F2,ABS(R2),VAL2)
!
       COSTHETA1=(R1**2-R2**2+DIS**2)/(2.D0*DIS*R1)
       COSTHETA2=-(R2**2-R1**2+DIS**2)/(2.D0*DIS*R2)
!      == SPHERICAL_PLGNDR RETURNS THE ASSOCIATED LEGENDRE POLYNOMIALS 
!      == MULTIPLIED WITH  SQRT( (L-M)! / (L+M)! ) FOR LM(L,M)=L(L+1)+1-M 
       CALL SPHERICAL_PLGNDR((L1+1)**2,L1,COSTHETA1,PLMWORK1)
       CALL SPHERICAL_PLGNDR((L2+1)**2,L2,COSTHETA2,PLMWORK2)
       PLM1=PLMWORK1(L1*(L1+1)-M1+1)
       PLM2=PLMWORK2(L2*(L2+1)-M2+1)
       PLM1=PLM1*SQRT(REAL(2*L1+1,KIND=8)/FPI)
       PLM2=PLM2*SQRT(REAL(2*L2+1,KIND=8)/FPI)
!!$       IF(M1.NE.0)PLM1=PLM1*SQ2
!!$       IF(M2.NE.0)PLM2=PLM2*SQ2  ! FIXED ERROR FROM PLM2=PLM1*SQ2 150128
!
       VALUE=2.D0*PI/DIS*DRDP * R1*VAL1*R2*VAL2*PLM1*PLM2
!IF(TPR)WRITE(*,FMT='("==",10F10.5)')DIS,VALUE,R1,VAL1,R2,VAL2,PLM1*PLM2,DRDP
       RETURN
       END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE NEWADAPT_MODULE
!*******************************************************************************
!**  ADAPTIVE INTEGRATION FOLLOWING                                           **
!**  A.C. GENZ AND A.A. MALIK, J. COMPUT. AND APPL. MATH 6, P.295 (1980)      **
!**                                                                           **
!**  ROUTINE NEWADAPT_INTEGRAND(P,RES) IS THE INTERFACE TO THE EXTERNAL       **
!**  FUNCTION EVALUATION. (THE CALL TO THE  EXTERNAL ROUTINE IS HARDWIRED)    **
!**                                                                           **
!**  THE INTEGRATION IS CALLED AS CALL NEWADAPT$EVALUATE(TOLERANCE,VALUE)     **
!**                                                                           **
!*******************************************************************************
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
!      ** SEE GENZ AND MALIK GENZ80_JCOMPUTAPPLMATH6_295                      **
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
!                           == USE N=2 IN GENZ+MALIK
       REAL(8),PARAMETER :: W1=-15264.D0/19683.D0
       REAL(8),PARAMETER :: W2=3920.D0/6561.D0    !=2940/19683
       REAL(8),PARAMETER :: W3=4080.D0/19683.D0
       REAL(8),PARAMETER :: W4=800.D0/19683.D0
       REAL(8),PARAMETER :: W5=6859.D0/19683.D0
!                           == W1P=W1PRIME
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
!      == FACTOR [-1,1]*[-1,2]=4 IS ALREADY ABSORBED IN THE WEIGHTS ============
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
       INTEGER(4),PARAMETER   :: LSTACKX=50000
       TYPE(SEGMENT_TYPE),ALLOCATABLE  :: STACK(:) !(LSTACKX)
       TYPE(SEGMENT_TYPE)     :: SEGMENT1
       TYPE(SEGMENT_TYPE)     :: SEGMENT2
       TYPE(SEGMENT_TYPE)     :: SEGMENT0
       INTEGER(4)             :: NSEGMENTS
       REAL(8)                :: ERROR ! CURRENT ERROR ESTIMATE
       REAL(8)                :: ERRMAX,ERRMIN
       INTEGER(4)             :: IAXIS
       INTEGER(4)             :: I
!      *************************************************************************
       CALL NEWADAPTINI()
       ALLOCATE(STACK(LSTACKX))
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
IF(STACK(1)%DIVIDEAXIS.EQ.0) THEN
  CALL ERROR$MSG('STACK(1)%DIVIDEAXIS=0')        
  CALL ERROR$I4VAL('NSEGMENTS',NSEGMENTS)
  CALL ERROR$STOP('NEWADAPT$EVALUATE')
END IF
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
           SEGMENT0=SEGMENT1
           SEGMENT1=SEGMENT2
           SEGMENT2=SEGMENT0
         END IF
!
!        =======================================================================
!        == EXPAND STACK AND PLACE NEW SEGMENTS ON DEFAULT POSITIONS          ==
!        =======================================================================
         NSEGMENTS=NSEGMENTS+1
         IF(NSEGMENTS.GT.LSTACKX) THEN
           CALL ERROR$MSG('NUMBER OF SEGMENTS EXCEEDS MAX')
           CALL ERROR$I4VAL('NSEGMENTS',NSEGMENTS)
           CALL ERROR$STOP('ADAPT$EVALUATE')
         END IF
         STACK(1)=SEGMENT1
         STACK(NSEGMENTS)=SEGMENT2
!
!        =======================================================================
!        == ENSURE THAT STACK(I)%ERR IS DESCENDING WITH I                     ==
!        == INSERT SEGMENT1 SEARCHING FROM TOP OF THE STACK                   ==
!        =======================================================================
         ERRMAX=SEGMENT1%ERR
         DO I=2,NSEGMENTS
!          __ENSURE THAT STACK(I-1).GE.STACK(I)_________________________________
           IF(STACK(I)%ERR.GT.ERRMAX) THEN  
             STACK(I-1)=STACK(I)
           ELSE
             STACK(I-1)=SEGMENT1    
             EXIT
           END IF
         ENDDO
!
!        =======================================================================
!        == INSERT SEGMENT2 SEARCHING FROM BOTTOM OF THE STACK                ==
!        =======================================================================
         ERRMIN=SEGMENT2%ERR
         DO I=NSEGMENTS-1,1,-1
!          __ENSURE THAT STACK(I).GE.STACK(I+1)_________________________________
           IF(STACK(I)%ERR.LT.ERRMIN) THEN
             STACK(I+1)=STACK(I)
           ELSE
             STACK(I+1)=SEGMENT1
             EXIT
           END IF
         ENDDO
       ENDDO
       CALL ERROR$MSG('LOOP NOT CONVERGED')
       CALL ERROR$STOP('NEWADAPT$EVALUATE')
       END
