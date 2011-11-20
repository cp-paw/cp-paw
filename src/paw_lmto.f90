MODULE LMTO_MODULE
LOGICAL(4)            :: TON=.false.       
LOGICAL(4)            :: TOFFSITE=.TRUE.  !include offsite exchange
LOGICAL(4)            :: Tdrop=.false. ! write the wave functions to file
REAL(8)   ,PARAMETER  :: K2=-0.25D0    ! 0.5*K2 IS THE KINETIC ENERGY
!REAL(8)   ,PARAMETER :: K2=-0.01D0    ! 0.5*K2 IS THE KINETIC ENERGY
!REAL(8)   ,PARAMETER :: RCSCALE=2.D0  !RADIUS SCALE FACTOR FOR NEIGHBORLIST
REAL(8)   ,PARAMETER  :: RCSCALE=1.2D0  !RADIUS SCALE FACTOR FOR NEIGHBORLIST
real(8)   ,parameter  :: tailed_lambda1=2.d0 ! 1st decay constant for tails
real(8)   ,parameter  :: tailed_lambda2=1.d0 ! 2nd decay constant for tails
!         RCSCALE TIMES THE SUM OF COVALENT RADII DEFINES CUTOFF FOR NEIGBORLIST
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
TYPE ORBITALGAUSSCOEFF_TYPE
  INTEGER(4)         :: NIJK     
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
  INTEGER(4),POINTER :: LOX(:)     ! (LNX)
  INTEGER(4),POINTER :: LNDOT(:)   ! (LNX)
  INTEGER(4),POINTER :: LMNDOT(:)  ! (LMNX)
  REAL(8)   ,POINTER :: AEF(:,:)   ! (NR,LNX)
  REAL(8)   ,POINTER :: PSF(:,:)   ! (NR,LNX)
  REAL(8)   ,POINTER :: NLF(:,:)   ! (NR,LNX)
  REAL(8)   ,POINTER :: U(:,:,:,:) ! (LMNX,LMNX,LMNX,LMNX)
  REAL(8)   ,POINTER :: OVERLAP(:,:) ! (LMNX,LMNX)
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
!
LOGICAL(4)              :: TINI=.FALSE.
LOGICAL(4)              :: TINISTRUC=.FALSE.
INTEGER(4)              :: NSP
INTEGER(4)              :: NRL   ! DIMENSION OF STRUCTURE CONSTANTS IN K-SPACE
INTEGER(4),ALLOCATABLE  :: LNX(:)               !(NSP)
INTEGER(4),ALLOCATABLE  :: LOX(:,:)             !(LNXX,NSP)
INTEGER(4),ALLOCATABLE  :: ISPECIES(:)          !(NAT)
REAL(8)   ,ALLOCATABLE  :: ORBRAD(:,:) !(LXX+1,NAT) NODE-POSITION OF THE ORBITAL
TYPE(POTPAR_TYPE)     ,ALLOCATABLE :: POTPAR(:) !POTENTIAL PARAMETERS
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: SBAR(:)   !(NNS) SCREENED STRUCTURE CONST.
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DENMAT(:) !(NND) DENSITY MATRIX
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: HAMIL(:)  !(NND) DERIVATIVE OF ENERGY
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: OVERLAP(:)!(NNS) OVERLAP MATRIX ONLY MAIN CHANNEL
!!$TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DENMAT_T(:) !(NNS) DENSITY MATRIX
!!$TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: HAMIL_T(:)  !(NNS) DERIVATIVE OF ENERGY
!!$INTEGER(4)                         :: NNUX      !DIMENSION OF UMAT
INTEGER(4)                         :: NNU       !#(ELEMENTS OF UMAT)
TYPE(UMAT_TYPE)       ,ALLOCATABLE :: UMAT(:)   !(NNUX/NNU) U-MATRIX ELEMENTS
INTEGER(4)            ,ALLOCATABLE :: SBARATOMI1(:)
INTEGER(4)            ,ALLOCATABLE :: SBARATOMI2(:)
INTEGER(4)            ,ALLOCATABLE :: SBARLI1(:,:)
!== GAUSSIAN PART OF NTBOS FROM SUPERPOSITION OF HANKEL FUNCTIONS ==============
TYPE(ORBITALGAUSSCOEFF_TYPE),ALLOCATABLE :: GAUSSORB(:) !(NAT)
!== GAUSSIAN PART OF NTBOS FROM TAILED HANKEL AND BESSEL FUNCTIONS =============
TYPE(ORBITALGAUSSCOEFF_TYPE),ALLOCATABLE :: GAUSSORB_T(:) !(NAT)
!== AUGMENTED NTBOS IN TERMS OF GAUSSIANS ======================================
TYPE(ORBITALGAUSSCOEFF_TYPE),ALLOCATABLE :: GAUSSORBAUG(:) !(NAT)
!
TYPE(ORBITALSPHHARM_TYPE),ALLOCATABLE :: LMORB(:)
!
TYPE(UTENSOR_TYPE)    ,ALLOCATABLE :: UTENSOR(:)
!
LOGICAL(4),PARAMETER :: TSPHERICAL=.FALSE.
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
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LMTO$GETL4')
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
      USE LMTO_MODULE, ONLY : TINI
      IMPLICIT NONE
      INTEGER(4) :: NAT
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
!     ==========================================================================
!     == ATTACH EXPONENTIAL TAILS TO AUGMENTED HANKEL AND BESSEL FUNCTIONS    ==
!     ==========================================================================
      CALL LMTO_MAKETAILEDPARTIALWAVES()
!!$!
!!$!     ==========================================================================
!!$!     == DETERMINE GAUSS EXPANSION OF UNSCREENED HANKEL FUNCTIONS KPRIME      ==
!!$!     ==========================================================================
!!$      CALL LMTO_GAUSSFITKPRIME()
!!$      CALL LMTO_GAUSSFITKAUGMENT()
!!$!
!!$!     ==========================================================================
!!$!     == DETERMINE KPRIME AND JBAR WITH EXPONENTIAL TAILS                     ==
!!$!     == USED TO BUILD UP APPROXIMATE NTBO'S IN A ONE-CENTER EXPANSION        ==
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
      USE LMTO_MODULE, ONLY :NSP,LNX,LOX,ISPECIES
      IMPLICIT NONE
      INTEGER(4)    :: ISP,NAT
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_COLLECTMAPARRAYS')
      CALL SETUP$NSPECIES(NSP)
      ALLOCATE(LNX(NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX(ISP))
        CALL SETUP$ISELECT(0)
      ENDDO
!
      ALLOCATE(LOX(MAXVAL(LNX),NSP))
      LOX(:,:)=-1
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4A('LOX',LNX(ISP),LOX(1:LNX(ISP),ISP))
        CALL SETUP$ISELECT(0)
      ENDDO

      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
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
!     **  THE ARRAYS ISPECIES1 AND LX1 HAVE STRANGE NAMES BECAUSE THE    **
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
      USE LMTO_MODULE, ONLY : K2,POTPAR,NSP,LNX,LOX &
     &                       ,tailed_lambda1,tailed_lambda2
      IMPLICIT NONE
      REAL(8)                :: LAMBDA1
      REAL(8)                :: LAMBDA2
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)                :: RAD
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
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_MAKETAILEDPARTIALWAVES')
!     == DECAY CONSTANTS FROM LMTO_MODULE 
      LAMBDA1=TAILED_LAMBDA1
      LAMBDA2=TAILED_LAMBDA2
!
      DO ISP=1,NSP
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
        ALLOCATE(POTPAR(ISP)%TAILED%AEF(NR,LMNXT))
        ALLOCATE(POTPAR(ISP)%TAILED%PSF(NR,LMNXT))
        ALLOCATE(POTPAR(ISP)%TAILED%NLF(NR,LMNXT))
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
        DEALLOCATE(R)
!
WRITE(STRING,FMT='(I5)')ISP
STRING='AETAILS_FORATOMTYPE'//TRIM(ADJUSTL(STRING))//'.DAT'
CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LNXT,POTPAR(ISP)%TAILED%AEF)
WRITE(STRING,FMT='(I5)')ISP
STRING='NLTAILS_FORATOMTYPE'//TRIM(ADJUSTL(STRING))//'.DAT'
CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LNXT,POTPAR(ISP)%TAILED%NLF)
!
!       ========================================================================
!       == ONSITE U-TENSOR OF TAILED PARTIAL WAVES                            ==
!       ========================================================================
        CALL SETUP$GETI4('LMRX',LMRX)
        LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
        ALLOCATE(POTPAR(ISP)%TAILED%U(LMNXT,LMNXT,LMNXT,LMNXT))
        ALLOCATE(ULITTLE(LRX+1,LNXT,LNXT,LNXT,LNXT))
        CALL LMTO_ULITTLE(GID,NR,LRX,LNXT,LOXT,POTPAR(ISP)%TAILED%AEF,ULITTLE)
        CALL LMTO_UTENSOR(LRX,LMNXT,LNXT,LOXT,ULITTLE,POTPAR(ISP)%TAILED%U)
        DEALLOCATE(ULITTLE)
!
!       ========================================================================
!       == ONSITE OVERLAP MATRIX                                              ==
!       ========================================================================
        ALLOCATE(POTPAR(ISP)%TAILED%OVERLAP(LMNXT,LMNXT))
        CALL LMTO_ONECENTEROVERLAP(GID,NR,LNXT,LOXT,POTPAR(ISP)%TAILED%AEF &
     &                            ,LMNXT,POTPAR(ISP)%TAILED%OVERLAP)
        DEALLOCATE(LNDOT)
        DEALLOCATE(LMNDOT)
        DEALLOCATE(LOXT)
      ENDDO
!
!     ==========================================================================
!     ==  CONSTRUCT GAUSSIAN FITS OF PRODUCTS OF TAILD PARTIAL WAVES AND      ==
!     ==  THEIRT POTENTIALS                                                   ==
!     ==========================================================================
      CALL LMTO_TAILEDPRODUCTS()
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
    &                                              *ULITTLE(L+1,LN2,LN4,LN3,LN1)
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
          IF(NPOW2.LT.1) cycle
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
              IF(NPOW2.LT.1) cycle
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
                  IF(NPOW2.LT.1) cycle
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
!CALL LMTO_TESTTAILEDP(nx)
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TESTTAILEDP(nx)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY: POTPAR,NSP
      IMPLICIT NONE
      integer(4),intent(in) :: nx
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
                    CALL LMTO_TAILEDINDEX_P(nx,LNX,LOX,LN1,LN2,LR1,IP1)
                    CALL LMTO_TAILEDINDEX_P(nx,LNX,LOX,LN3,LN4,LR2,IP2)
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
                    CALL LMTO_TAILEDINDEX_T(nx,LNX,LOX,LN1,LN2,LR1,LN3,LR2,IT)
                    CALL LMTO_TAILEDINDEX_S(nx,LNX,LOX,LN4,IS)
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
STOP 'forced in LMTO_TESTTAILEDP'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE EXPINTEGRAL(N,E,RES)
!     **************************************************************************
!     **************************************************************************
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: E
      LOGICAL(4)            :: TEVEN
      REAL(8)               :: PI
      REAL(8)   ,INTENT(OUT):: RES
      INTEGER(4)            :: K
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
      INTEGER(4)             :: npow2
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
                NPOW2=INT(0.5D0*REAL(NX-LR2a)+1.000001D0) !R^(L+2N), N=0,NPOW-1
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
      INTEGER(4)             :: npow2
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
      integer(4)             :: ln1a,l1,npow2
!     **************************************************************************
      IS=0
      DO LN1A=1,LNX
        L1=LOX(LN1A)
        NPOW2=INT(0.5D0*REAL(NX-L1)+1.000001D0)  !R^(L+2N), N=0,NPOW-1
        IF(NPOW2.LT.1) CYCLE
        IS=IS+1
        IF(LN1a.LT.LN1) CYCLE
        IF(LN1A.NE.ln1) then
          CALL ERROR$I4VAL('LN1',LN1)
          CALL ERROR$I4VAL('LN1A',LN1A)
          CALL ERROR$STOP('LMTO_TAILEDINDEX_S')
        END IF
        RETURN
      ENDDO
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
      INTEGER(4)             :: npow2
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
!     **************************************************************************
      CALL SETUP$GETL4('INTERNALSETUPS',TCHK)
      IF(.NOT.TCHK) RETURN
                              CALL TRACE$PUSH('LMTO$MAKESTRUCTURECONSTANTS')
                              CALL TIMING$CLOCKON('LMTO STRUCTURECONSTANTS')
      TINISTRUC=.TRUE.
!
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
        CALL SETUP$ISELECT(0)
      ENDDO
!
!     ==========================================================================
!     == NEIGHBORLIST   NNLIST(:,NN)=(IAT1,IAT2,IT1,IT2,IT3)                  ==
!     ==                     it1,it2,it3 are the lattice translations of iat2 ==
!     ==========================================================================
      ALLOCATE(RC(NAT))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',SVAR)
        CALL PERIODICTABLE$GET(SVAR,'R(COV)',RC(IAT))
      ENDDO
      RC(:)=RC(:)*RCSCALE
      NNX=NNXPERATOM*NAT
      ALLOCATE(NNLIST(5,NNX))
      CALL LMTO$NEIGHBORLIST(RBAS,NAT,R0,RC,NNX,NNS,NNLIST)
      DEALLOCATE(RC)
!
!     ==========================================================================
!     == STRUCTURE CONSTANTS                                                  ==
!     ==========================================================================
      IF(ALLOCATED(SBAR)) THEN
        DO I=1,SIZE(SBAR)
          DEALLOCATE(SBAR(I)%MAT)
        ENDDO
        DEALLOCATE(SBAR)
      END IF
      ALLOCATE(SBAR(NNS))
      DO IAT1=1,NAT
!       == MEMBERS NN1:NN2 ON THE NEIGHBOLIST BuILD THE CLUSTER AROUND ATOM 1 ==
!       == MEMBER NN0 IS THE ONSITE MEMBER                                    ==
        NN1=1
        NN0=0
        DO NN=1,NNS
          IF(NNLIST(1,NN).GT.IAT1)EXIT
          NN2=NN   ! nn2 is the last member with iat1 as first atom 
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
          LX1(NN-NN1+1)=LX ! x(angular momentum for this atom)
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
        ALLOCATE(SBAR1(N,NORB))
!NORB=(LX1(1)+1)**2
!N=(LX(1)+2)**2?
PRINT*,'DOING LMTO$CLUSTERSTRUCTURECONSTANTS.....'
        CALL LMTO$CLUSTERSTRUCTURECONSTANTS(K2,NAT1,RPOS,LX1,QBAR,N,NORB,SBAR1)
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
          ALLOCATE(SBAR(NN)%MAT(LMX2,LMX1))
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
              SBAR(NN)%MAT(J21:J22,J11:J12)=SBAR1(I21:I22,I11:I12)
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
                              CALL TIMING$CLOCKOFF('LMTO STRUCTURECONSTANTS')
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
      USE LMTO_MODULE,ONLY : TON,NRL,ISPECIES,LNX,LOX,SBARATOMI1,SBARLI1
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
      g=conjg(g)
      h=conjg(h)
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
        G(:,I)=C(:)*SBAR(:,I)*F(I)
        G(I,I)=G(I,I)+E(I)
      ENDDO
      CALL LIB$INVERTC8(NRL,G,H)
      G(:,:)=H(:,:)
      DO I=1,NRL
        H(:,I)=H(:,I)*C(I)
      ENDDO
      H=MATMUL(H,SBAR)
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
        KR=-2.D0*PI*SUM(REAL(SBAR(NN)%IT(:))*XK(:))
        I1OFAT1=SBARATOMI1(IAT1)
        I2OFAT1=SBARATOMI2(IAT1)
        I1OFAT2=SBARATOMI1(IAT2)
        I2OFAT2=SBARATOMI2(IAT2)
        EIKR=CMPLX(COS(KR),SIN(KR)) 
        SOFK(I1OFAT2:I2OFAT2,I1OFAT1:I2OFAT1) &
     &                  =SOFK(I1OFAT2:I2OFAT2,I1OFAT1:I2OFAT1)+SBAR(NN)%MAT*EIKR
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_AOFK(NNA,ALIST,XK,N,AOFK)
!     **************************************************************************
!     **  TRANSFORMS THE SCREENED STRUCTURE CONSTANTS INTO K-SPACE            **
!     **************************************************************************
      USE LMTO_MODULE, ONLY: PERIODICMAT_TYPE,POTPAR_TYPE,SBARATOMI1,SBARATOMI2
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)             :: NNA    ! #(ELEMNTS IN NEIGHBORLIST)
      TYPE(PERIODICMAT_TYPE),INTENT(IN) :: ALIST(NNA) !(NNA)
      INTEGER(4)            ,INTENT(IN) :: N
      REAL(8)               ,INTENT(IN) :: XK(3) !K-POINT IN FRACTIONAL COORD.
      COMPLEX(8)            ,INTENT(OUT):: AOFK(N,N)
      REAL(8)                :: KR
      COMPLEX(8)             :: EIKR
      INTEGER(4)             :: IAT1,IAT2
      INTEGER(4)             :: I1OFAT1,I2OFAT1        
      INTEGER(4)             :: I1OFAT2,I2OFAT2        
      INTEGER(4)             :: NN
      REAL(8)                :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      IF(N.NE.MAXVAL(SBARATOMI2)) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$I4VAL('MAXVAL(SBARATOMI2)',MAXVAL(SBARATOMI2))
        CALL ERROR$MSG('LMTO_SOFK')
      END IF
      AOFK(:,:)=(0.D0,0.D0)
      DO NN=1,NNA
        IAT1=ALIST(NN)%IAT1
        IAT2=ALIST(NN)%IAT2
        KR=2.D0*PI*SUM(REAL(ALIST(NN)%IT(:))*XK(:))
        I1OFAT1=SBARATOMI1(IAT1)
        I2OFAT1=SBARATOMI2(IAT1)
        I1OFAT2=SBARATOMI1(IAT2)
        I2OFAT2=SBARATOMI2(IAT2)
        EIKR=CMPLX(COS(KR),SIN(KR)) 
        AOFK(I1OFAT2:I2OFAT2,I1OFAT1:I2OFAT1) &
     &           =AOFK(I1OFAT2:I2OFAT2,I1OFAT1:I2OFAT1)+ALIST(NN)%MAT*EIKR
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_MATRTOK(NNA,ALIST,XK,N,AOFK)
!     **************************************************************************
!     **  TRANSFORMS A REAL SPACE MATRIX TO K-SPACE                           ==
!     **************************************************************************
      USE LMTO_MODULE, ONLY: PERIODICMAT_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)             :: NNA    ! #(ELEMNTS IN NEIGHBORLIST)
      TYPE(PERIODICMAT_TYPE),INTENT(IN) :: ALIST(NNA) !(NNA)
      INTEGER(4)            ,INTENT(IN) :: N
      REAL(8)               ,INTENT(IN) :: XK(3) !K-POINT IN FRACTIONAL COORD.
      COMPLEX(8)            ,INTENT(OUT):: AOFK(N,N)
      REAL(8)                           :: KR
      COMPLEX(8)                        :: EIKR
      INTEGER(4)                        :: NAT
      INTEGER(4)                        :: IAT1,IAT2
      INTEGER(4)            ,ALLOCATABLE:: LENG(:),IP1(:),IP2(:)
      INTEGER(4)                        :: NN,IAT,ISVAR
      REAL(8)                           :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      NAT=MAX(MAXVAL(ALIST(:)%IAT1),MAXVAL(ALIST(:)%IAT2))
      ALLOCATE(LENG(NAT))
      DO NN=1,NNA
        LENG(ALIST(NN)%IAT1)=ALIST(NN)%N1
        LENG(ALIST(NN)%IAT2)=ALIST(NN)%N2
      ENDDO       
      IF(SUM(LENG).NE.N) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$I4VAL('SUM(LENG)',SUM(LENG))
        CALL ERROR$MSG('LMTO_MATRTOK')
      END IF
      ALLOCATE(IP1(NAT))
      ALLOCATE(IP2(NAT))
      ISVAR=1
      DO IAT=1,NAT
        IP1(IAT)=ISVAR
        ISVAR=ISVAR+LENG(IAT)
        IP2(IAT)=ISVAR-1
      ENDDO
!
      AOFK(:,:)=(0.D0,0.D0)
      DO NN=1,NNA
        IAT1=ALIST(NN)%IAT1
        IAT2=ALIST(NN)%IAT2
        KR=2.D0*PI*SUM(REAL(ALIST(NN)%IT(:))*XK(:))
        EIKR=CMPLX(COS(KR),SIN(KR)) 
        AOFK(IP1(IAT1):IP2(IAT1),IP1(IAT2):IP2(IAT2)) &
     &                         =AOFK(IP1(IAT1):IP2(IAT1),IP1(IAT2):IP2(IAT2)) &
     &                         +ALIST(NN)%MAT*EIKR
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_prproj(id)
!     **************************************************************************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,THIS,WAVES_SELECTWV
      USE LMTO_MODULE, ONLY : LOX,LNX,ISPECIES
      implicit none
      character(*),intent(in) :: id
      integer(4)             :: nat
      integer(4),allocatable :: ipro1(:)
      integer(4),allocatable :: nproat(:)
      integer(4)             :: iat,isp,ikpt,ispin,ib,ipro,nb,nbh,i0
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
          PRINT*,"==============IKPT=",IKPT,' iat=',iat,' id=',id
          do iat=1,nat
            I0=IPRO1(IAT)-1
            DO IB=1,NBH
              WRITE(*,FMT='(I3,40("(",2F10.5,")"))')IB,THIS%proj(:,IB,I0+1:I0+NPROAT(IAT))
            ENDDO
          enddo
        enddo
      enddo
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBODENMAT()
!     **************************************************************************
!     **  CONSTRUCT DENSITY MATRIX IN A NTBO BASIS                            **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : SBAR,DENMAT,LOX,LNX,ISPECIES
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
      LOGICAL(4),PARAMETER   :: TPR=.false.
complex(8)  :: phase
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
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          NBH=THIS%NBH
          NB=THIS%NB
          DO NN=1,NND
            IAT1=DENMAT(NN)%IAT1
            IAT2=DENMAT(NN)%IAT2
            IT=DENMAT(NN)%IT
            SVAR=2.D0*PI*SUM(XK(:,IKPT)*REAL(IT,KIND=8))
            EIKR=EXP(CI*SVAR)  !<P_{R+T}|PSI>=<P_R|PSI>*EIKR
            I0=IPRO1(IAT1)-1
            J0=IPRO1(IAT2)-1

!IF(IKPT.NE.2.AND.IKPT.NE.3.AND.IKPT.NE.4.AND.IKPT.NE.7.AND.IKPT.NE.10.AND.IKPT.NE.19) CYCLE
!!$IF(IAT1.EQ.1.AND.IAT2.EQ.1.AND.SUM(IT**2).EQ.0) THEN
!!$  PRINT*,"==============IKPT=",IKPT," XK=",XK(:,IKPT),I0,J0
!!$  DO IB=1,NB
!!$    IF(OCC(IB,IKPT,ISPIN).LT.1.D-5) CYCLE
!!$    PHASE=THIS%TBC(1,IB,I0+1)
!!$    PHASE=CONJG(PHASE)/SQRT((PHASE*CONJG(PHASE)))
!!$    WRITE(*,FMT='(I3,40("(",2F10.5,")"))')IB,THIS%tbc(:,IB,I0+1:I0+NPROAT(IAT1))
!!$  ENDDO
!!$END IF
            DO I=1,NPROAT(IAT1)
              DO J=1,NPROAT(IAT2)
                IF(TINV) THEN
                  CSVAR22=(0.D0,0.D0)
                  DO IBH=1,NBH
                    F1=OCC(2*IBH-1,IKPT,ISPIN)
                    F2=OCC(2*IBH,IKPT,ISPIN)
                    C1(:)=THIS%TBC(:,IBH,I0+I)
                    C2(:)=THIS%TBC(:,IBH,J0+J)*EIKR
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
                    C2(:)=THIS%TBC(:,IB,J0+J)*EIKR
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
      USE LMTO_MODULE, ONLY : SBAR,HAMIL,LOX,LNX,ISPECIES
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
                CSVAR22(:,:)=CSVAR22(:,:)*CONJG(EIKR)
!                CSVAR22(:,:)=CSVAR22(:,:)*EIKR
!
                DO IB=1,NBH
                  DO IDIM=1,NDIM
                    DO JDIM=1,NDIM
                      THIS%HTBC(JDIM,IB,J0+J)=THIS%HTBC(JDIM,IB,J0+J) &
      &                              +CSVAR22(IDIM,JDIM)*THIS%TBC(IDIM,IB,I0+I)
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
      integer(4)  ::nn,nnd
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
      SUBROUTINE LMTO_onsortho(iat,lmnx,t,unt)
!     **************************************************************************
!     **  CONSTRUCT transformation onto onsite orthonormal orbitals           **
!     **                                                                      **
!     **   |chi'_i>=sum_j |chi_j>T_ji;  <chi'_i|chi'_j>=delta_ij              **
!     **   t*unt=1                                                            **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : ISPECIES,potpar
      IMPLICIT NONE
      INTEGER(4),intent(in) :: iat
      INTEGER(4),intent(in) :: lmnx
      real(8)   ,intent(out):: t(lmnx,lmnx)
      real(8)   ,intent(out):: unt(lmnx,lmnx)
      integer(4)            :: isp
      integer(4)            :: lmnxt
      real(8)               :: overlap(lmnx,lmnx)
      real(8)               :: overlap1(lmnx,lmnx)
      real(8)               :: svar
      integer(4)            :: i,j
!     **************************************************************************
      isp=ispecies(iat)
!
!     ==========================================================================
!     == construct overlap matrix                                             ==
!     ==========================================================================
      LMNXT=POTPAR(ISP)%TAILED%LMNX
      call LMTO_SHRINKDOWNHTNL(IAT,IAT,1 &
   &               ,LMNXT,LMNXT,POTPAR(ISP)%TAILED%OVERLAP,LMNX,LMNX,overlap)
      DO I=1,lmnx
        WRITE(*,FMT='(I3," O=",100("(",2F8.4,")"))')I,overlap(i,:)
      ENDDO
!
!     ==========================================================================
!     == orthonormalize                                                       ==
!     ==========================================================================
      T(:,:)=0.D0
      DO I=1,LMNX
        T(I,I)=1.D0
      ENDDO
      DO I=1,LMNX
!       == ORTHOGONALIZE
        DO J=1,I-1
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
!     ==========================================================================
!     == invert                                                               ==
!     ==========================================================================
      call lib$invertr8(lmnx,t,unt)
      return
      end
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
      USE LMTO_MODULE, ONLY : tdrop,SBAR,DENMAT,LOX,LNX,ISPECIES,potpar
      IMPLICIT NONE
      type T_type
        real(8),pointer :: mat(:,:)
        real(8),pointer :: inv(:,:)
      end type T_type
      INTEGER(4)             :: NPRO
      LOGICAL(4)             :: TINV
      integer(4)             :: nfil
      INTEGER(4)             :: NB,NBH
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      INTEGER(4),ALLOCATABLE :: IPRO1(:)
      INTEGER(4),ALLOCATABLE :: NPROAT(:)
      complex(8),allocatable :: psi(:,:,:)
      complex(8),allocatable :: psi1(:,:,:)
      real(8)   ,allocatable :: eigval(:)
      real(8)   ,allocatable :: overlap(:,:)
      real(8)   ,allocatable :: overlap1(:,:)
      real(8)   ,allocatable :: trans(:,:)
      real(8)   ,allocatable :: wkpt(:)
      type(t_type),allocatable :: t(:)
      complex(8),allocatable :: eigvec(:,:)
      character(16)          :: id
      INTEGER(4)             :: NAT
      INTEGER(4)             :: nkpt
      INTEGER(4)             :: lmnxt,lmnx
      INTEGER(4)             :: i,j,i1,i2,idim,ipro,iat,isp,ikpt,ispin,ib,ib1,ibh
      integer(4)             :: iwork16(16)
      real(8)                :: ev
      real(8)                :: svar
!     **************************************************************************
                                           CALL TRACE$PUSH('LMTO_NTBODENMAT')
      IF(.NOT.TDROP) RETURN
      CALL CONSTANTS('EV',EV)
!
!     ==========================================================================
!     ==  ATTACH FILE                                                         ==
!     ==========================================================================
      NFIL=12
      OPEN(NFIL,FILE='NTBOWV.DATA')
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      ALLOCATE(XK(3,NKPTL))
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      nkpt=nkptl
      IF(NKPT.NE.NKPTL) THEN
        CALL ERROR$MSG('ROUTINE NOT SUITED FOR PARALLEL CALCULATIONS')
        CALL ERROR$STOP('LMTO_DROP')
      end if
      allocate(wkpt(nkpt))
      CALL DYNOCC$GETR8A('WKPT',NKPT,WKPT)
!     == multiply with spin multiplicity =======================================      
      if(nspin.eq.2) wkpt=2.d0*wkpt
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
        if(lnx(isp).gt.16) then
          CALL ERROR$MSG('INTERNAL LIMIT FOR VARIABLE LNX EXCEEDED')
          CALL ERROR$STOP('LMTO_DROP')
        end if
        IWORK16(:LNX(ISP))=LOX(:LNX(ISP),ISP)
        ID='PSIINFO'
        WRITE(NFIL,*)ID,IAT,IPRO1(IAT),NPROAT(IAT),LNX(ISP),IWORK16
!
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
!
!     ==========================================================================
!     ==  overlap matrix                                                      ==
!     ==========================================================================
      allocate(t(nat))
      do iat=1,nat
        isp=ispecies(iat)
        lmnx=nproat(iat)
        ALLOCATE(t(iat)%mat(LMNX,LMNX))
        ALLOCATE(t(iat)%inv(LMNX,LMNX))
        call LMTO_onsortho(iat,lmnx,t(iat)%mat,t(iat)%inv)
!       ====================================================================
!!$        LMNXT=POTPAR(ISP)%TAILED%LMNX
!!$        ALLOCATE(overlap(LMNX,LMNX))
!!$        call LMTO_SHRINKDOWNHTNL(IAT,IAT,1 &
!!$     &               ,LMNXT,LMNXT,POTPAR(ISP)%TAILED%OVERLAP,LMNX,LMNX,overlap)
!!$        ALLOCATE(OVERLAP1(LMNX,LMNX))
!!$        ALLOCATE(trans(LMNX,LMNX))
!!$        TRANS(:,:)=0.D0
!!$        DO I=1,LMNX
!!$          TRANS(I,I)=1.D0
!!$        ENDDO
!!$        DO I=1,LMNX
!!$!         == ORTHOGONALIZE
!!$          DO J=1,I-1
!!$            SVAR=OVERLAP(J,I)
!!$            TRANS(:,I)=TRANS(:,I)-TRANS(:,J)*SVAR
!!$            OVERLAP1=OVERLAP
!!$            OVERLAP1(:,I)=OVERLAP1(:,I)-OVERLAP(:,J)*SVAR
!!$            OVERLAP1(I,:)=OVERLAP1(I,:)-OVERLAP(J,:)*SVAR
!!$            OVERLAP(:,:)=OVERLAP1(:,:)
!!$          ENDDO
!!$!         == NORMALIZE
!!$          SVAR=1.D0/SQRT(OVERLAP(I,I))
!!$          TRANS(:,I)=TRANS(:,I)*SVAR
!!$          OVERLAP(:,I)=OVERLAP(:,I)*SVAR
!!$          OVERLAP(I,:)=OVERLAP(I,:)*SVAR
!!$        ENDDO
!!$        t(iat)%mat(:,:)=trans(:,:)
!!$        call lib$invertr8(lmnx,t(iat)%mat,t(iat)%inv)
!!$!
!!$        deallocate(overlap)
!!$        deallocate(overlap1)
!!$        deallocate(trans)
      enddo
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
!         ==checks =============================================================
          IF(.NOT.ASSOCIATED(THIS%TBC)) THEN
            CALL ERROR$MSG('THIS%TBC NOT ASSOCIATED')
            CALL ERROR$STOP('LMTO_DROP')
          END IF
          IF(.NOT.ASSOCIATED(THIS%rlam0)) THEN
            CALL ERROR$MSG('THIS%RLAM0 NOT ASSOCIATED')
            CALL ERROR$STOP('LMTO_DROP')
          END IF
!
!         ==  WRITE WAVE FUNCTIONS =============================================
          ALLOCATE(PSI(NDIM,NPRO,NB))
          ALLOCATE(PSI1(NDIM,NPRO,NB))
          ALLOCATE(eigval(NB))
          ALLOCATE(eigvec(nb,NB))
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
!!$          DO IB=1,NB
!!$            WRITE(*,FMT='(I3," H=",100("(",2F8.4,")"))')IB,THIS%RLAM0(IB,:)/EV
!!$          ENDDO
          CALL LIB$DIAGC8(NB,THIS%RLAM0,EIGVAL,EIGVEC)
!!$          DO IB=1,NB
!!$            WRITE(*,FMT='(I3,"E=",F6.3," U=",100("(",2F8.4,")"))')IB,EIGVAL(IB)/EV,EIGVEC(:,IB)
!!$          ENDDO
!
!         == transform wave function to eigenstates ===========================
          psi1=psi
          psi=(0.d0,0.d0)
          do ib=1,nb
            do ib1=1,nb
              psi(:,:,ib)=psi(:,:,ib)+psi1(:,:,ib1)*eigvec(ib1,ib)
            enddo
          enddo
!
!         == transform to onsite-orthogonal orbitals ==========================
          psi1=psi
          do iat=1,nat
            i1=ipro1(iat)
            i2=i1-1+nproat(iat)
            do ib=1,nb
              do idim=1,ndim
                psi(idim,i1:i2,ib)=matmul(t(iat)%inv,psi1(idim,i1:i2,ib))
              enddo
            enddo
          enddo
!
!         == write wave functions to file =====================================
          ID='KINFO'
          WRITE(NFIL,*)ID,NDIM,NB,NPRO,wkpt(ikpt)
          ID='PSI'
          DO IB=1,NB
            WRITE(NFIL,*)ID,IB,IKPT,ISPIN,EIGVAL(IB),PSI(:,:,IB)
print*,'norm ',ikpt,ispin,ib,dot_product(psi(1,:,ib),psi(1,:,ib))
          ENDDO
          DEALLOCATE(eigval)
          DEALLOCATE(eigvec)
          DEALLOCATE(PSI)
          DEALLOCATE(PSI1)
        ENDDO
      ENDDO
      DO IAT=1,NAT
        DEALLOCATE(T(IAT)%MAT)
        DEALLOCATE(T(IAT)%inv)
      ENDDO
      DEALLOCATE(T)
      CLOSE(NFIL)
!
PRINT*,'TEST READING '
  CALL LMTO_DROPREAD()
      CALL ERROR$MSG('FORCED STOP AFTER WRITING FILE')
      CALL ERROR$STOP('LMTO$DROP')
                                           CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      subroutine lmto_dropread()
!     **************************************************************************
!     ** reads the wave functions in local orbitals, which have been          **
!     ** onsite-orthonormalized.                                              **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      character(16)          :: id
      integer(4)             :: nat       ! #(atoms)
      integer(4)             :: ndim      ! #(spinor-components)
      integer(4)             :: nspin     ! #(spin directions)
      integer(4)             :: nkpt      ! #(k-points)
      integer(4)             :: npro      ! total #(orbitals)
      integer(4)             :: nb        ! #(states)
      integer(4),allocatable :: ipro1(:)  ! first orbital index for this atom
      integer(4),allocatable :: nproat(:) ! #(orbitals) for this atom
      integer(4),allocatable :: lnx(:)    ! #(different orbital shells)
      integer(4),allocatable :: lox(:,:)  ! angular momentum of this orbital shell
      real(8)   ,allocatable :: e(:)      ! energy eigenvalues
      complex(8),allocatable :: c(:,:,:)  ! orbital coefficients
      integer(4)             :: iat,ib,ikpt,ispin
      integer(4)             :: iat1,ib1,ikpt1,ispin1
      integer(4)             :: nfil  ! fortran file unit
!     **************************************************************************
      NFIL=12
      OPEN(NFIL,FILE='NTBOWV.DATA')
      REWIND(nfil)
!
      READ(NFIL,*)ID,NAT,NDIM,NSPIN,NKPT
      IF(ID.NE.'INFO') THEN
        STOP 'INCORRECT ID_1'
      END IF
      ALLOCATE(IPRO1(nAT))
      ALLOCATE(NPROAT(nAT))
      ALLOCATE(LNX(nAT))
      ALLOCATE(LOX(16,nAT))
      DO IAT=1,NAT
        READ(NFIL,*)ID,Iat1,IPRO1(IAT),NPROAT(IAT),LNX(IAT),LOX(:,IAT)
        IF(ID.NE.'PSIINFO') THEN
          STOP 'INCORRECT ID_2'
        END IF
        IF(Iat1.NE.IAT) THEN
          STOP 'ATOM INDEX OUT OF ORDER'
        END IF
      ENDDO
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          READ(NFIL,*)ID,NDIM,NB,NPRO
          IF(ID.NE.'KINFO') THEN
            STOP 'INCORRECT ID_3'
          END IF
          allocate(e(nb))
          allocate(c(ndim,npro,nb))
          DO IB=1,NB
            read(NFIL,*)ID,IB1,IKPT1,ISPIN1,E(IB),c(:,:,IB)
            IF(ID.NE.'PSI') THEN
              print*,'id   =',trim(id)
              print*,'ib   =',ib1,ib
              print*,'ikpt =',ikpt1,ikpt
              print*,'ispin=',ispin1,ispin
              STOP 'INCORRECT ID_4'
            END IF
          ENDDO
!
!         ------> grab wave functions and energies from here <-----
! 
          deallocate(e)
          deallocate(c)
        ENDDO
      ENDDO
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
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMNXX_
      INTEGER(4),INTENT(IN) :: NDIMD_
      INTEGER(4),INTENT(IN) :: NAT_
      COMPLEX(8),INTENT(IN) :: DENMAT_(LMNXX_,LMNXX_,NDIMD_,NAT_)
      LOGICAL(4),SAVE       :: TFIRSTENERGY=.TRUE.
      INTEGER(4)            :: SWITCH
!     **************************************************************************
      IF(.NOT.TON) RETURN
      WRITE(*,FMT='(82("="),T30," LMTO$ENERGY START ")')
!
!     ==========================================================================
!     ==  
!     ==========================================================================
      call lmto_drop()
!
!     ==========================================================================
!     ==  
!     ==========================================================================
      CALL TIMING$CLOCKON('NTBODENMAT')
      CALL LMTO_NTBODENMAT()
      CALL TIMING$CLOCKOFF('NTBODENMAT')
      CALL LMTO_TESTDENMAT_1CDENMAT(LMNXX_,NDIMD_,NAT_,DENMAT_)
      CALL LMTO_TESTDENMAT()
!
!     ==========================================================================
!     ==  some info                                                           ==
!     ==========================================================================
!!$      CALL LMTO$REPORTPOTBAR(6)
!!$      CALL LMTO$REPORTSBAR(6)
!!$      CALL LMTO$REPORTDENMAT(6)
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
!
!     ==========================================================================
!     ==  CONVERT HAMIL INTO HTBC
!     ==========================================================================
      CALL TIMING$CLOCKON('NTBODENMATDER')
      CALL LMTO_NTBODENMATDER()
      CALL TIMING$CLOCKOFF('NTBODENMATDER')
!
!     ==========================================================================
!     ==  CLEAN DENMAT AND HAMIL                                              ==
!     ==========================================================================
      CALL LMTO_CLEANDENMAT()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SIMPLEENERGYTEST2()
!     **************************************************************************
!     **  WORK OUT THE ENERGY USING THE LOCAL APPROXIMATION                   **
!     **  TAILED PARTIAL WAVES                                                **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,DENMAT,HAMIL,LNX,LOX,POTPAR,toffsite
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      LOGICAL(4),PARAMETER  :: TPlot=.false.
      REAL(8)   ,PARAMETER  :: HFWEIGHT=0.25D0
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
      REAL(8)   ,ALLOCATABLE:: AECORE(:)
      REAL(8)   ,ALLOCATABLE:: ULITTLE(:,:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: U(:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: D(:,:,:)
      REAL(8)   ,ALLOCATABLE:: H(:,:,:)
      REAL(8)               :: EH,EX,EXTOT,EHTOT,EHDC,EXDC,Q
      INTEGER(4)            :: NN,IAT,I,J,K,L,IS,IAT1,IAT2,ISP
      real(8)               :: qspin(4)
CHARACTER(128) :: STRING
REAL(8)   ,ALLOCATABLE:: t(:,:),unt(:,:),mymat(:,:)
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
!       ========================================================================
!       == CALCULATE U-TENSOR                                                 ==
!       ========================================================================
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(AECORE(NR))
        CALL SETUP$GETR8A('AECORE',NR,AECORE)
        CALL SETUP$GETI4('LMRX',LMRX)
        LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
        ALLOCATE(D(LMNX,LMNX,NDIMD))
        ALLOCATE(H(LMNX,LMNX,NDIMD))
        D=DENMAT(IND)%MAT
        H(:,:,:)=0.D0
!!$DO I=1,LMNX
!!$  WRITE(*,FMT='("D=",30F10.5)')D(I,:,1)
!!$ENDDO
!!$allocate(t(lmnx,lmnx))
!!$allocate(unt(lmnx,lmnx))
!!$allocate(mymat(lmnx,lmnx))
!!$call LMTO_onsortho(iat,lmnx,t,unt)
!!$do i=1,ndimd
!!$MYMAT=MATMUL(unT,MATMUL(D(:,:,i),TRANSPOSE(UNT)))
!!$print*,'============iat=',iat,'  idim=',i,'==================='
!!$DO j=1,LMNX
!!$  WRITE(*,FMT='("D(MY)=",30F10.5)')MYMAT(j,:)
!!$ENDDO
!!$enddo
!!$deallocate(t)
!!$deallocate(unt)
!!$deallocate(mymat)
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
        HAMIL(INH)%MAT=0.d0
        EH=0.D0
        EX=0.D0
        Qspin=0.D0
        HT(:,:,:)=0.D0
        DO I=1,LMNXT
          DO J=1,LMNXT
            qspin(:ndimd)=qspin(:ndimd)+POTPAR(ISP)%TAILED%OVERLAP(I,J)*DT(J,I,:)
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
                DO IS=1,NDIMD
                  EX=EX-0.25D0*U(I,J,K,L)*DT(K,J,IS)*DT(L,I,IS)
                  HT(K,J,IS)=HT(K,J,IS)-0.25D0*U(I,J,K,L)*DT(L,I,IS) 
                  HT(L,I,IS)=HT(L,I,IS)-0.25D0*U(I,J,K,L)*DT(K,J,IS) 
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,NDIMD,LMNXT,LMNXT,HT,LMNX,LMNX,H)
        HAMIL(INH)%MAT=HAMIL(INH)%MAT+H
        EXTOT=EXTOT+EX
PRINT*,'TOTAL CHARGE ON ATOM=                 ',IAT,QSPIN(1)
PRINT*,'TOTAL SPIN[HBAR/2] ON ATOM=           ',IAT,QSPIN(2:NDIMD)
PRINT*,'EXACT EXCHANGE ENERGY FOR ATOM=       ',IAT,EX
!
!       ========================================================================
!       == ADD CORE VALENCE EXCHANGE                                          ==
!       ========================================================================
        CALL LMTO_CVX(ISP,LMNX,EX,D(:,:,1),H(:,:,1))
        EXTOT=EXTOT+EX
        HAMIL(INH)%MAT(:,:,1)=HAMIL(INH)%MAT(:,:,1)+H(:,:,1)
PRINT*,'CORE VALENCE EXCHANGE ENERGY FOR ATOM=',IAT,EX
!
!       ========================================================================
!       == DOUBLE COUNTING CORRECTION (EXCHANGE ONLY)                         ==
!       ========================================================================
! THIS IS THE TIME CONSUMIN PART OF ENERGYTEST
CALL TIMING$CLOCKON('ENERGYTEST:DC')      
        CALL LMTO_SIMPLEDC(GID,NR,LMNXT,LNXT,LOXT,POTPAR(ISP)%TAILED%AEF &
     &                    ,LRX,AECORE,DT,EX,HT)
!        CALL LMTO_BLOWDOWNHT(IAT,NDIMD,LMNXT,HT,LMNX,H)
!       == THE FOLLOWING IS MORE GENERAL AND SHALL REPLACE THE PREVIOUS ONE
        CALL LMTO_SHRINKDOWNHTNL(IAT,IAT,NDIMD,LMNXT,LMNXT,HT,LMNX,LMNX,H)
        EXTOT=EXTOT-EX
        HAMIL(INH)%MAT=HAMIL(INH)%MAT-H
PRINT*,'DOUBLE COUNTING CORRECTION ENERGY FOR ATOM=',IAT,EX
CALL TIMING$CLOCKOFF('ENERGYTEST:DC')      
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
        CALL LMTO_OFFSITEX(EX)
        EXTOT=EXTOT+EX
      END IF
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
!     == plot wave functions                                                  ==
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
      REAL(8)   ,ALLOCATABLE :: ZA(:)
      REAL(8)   ,ALLOCATABLE :: ZB(:)
      REAL(8)   ,ALLOCATABLE :: D(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: DB(:,:,:)
      REAL(8)   ,ALLOCATABLE :: H(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: HB(:,:,:)
      REAL(8)   ,ALLOCATABLE :: NNAT(:)    !(NAT) POINTER TO ONSITE TERMS
INTEGER(4) :: I,J,K,L,I1,J1,K1
REAL(8) :: SVAR1,SVAR2
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
        ALLOCATE(ZB(NIJKB*NEB))
        H=0.D0
        HA=0.D0
        HB=0.D0
        DO LMN1B=1,LMNXB
          DO LMN2B=1,LMNXB
            ZA(:)=MATMUL(SAB,SPECIAL(ISPB)%APOT(:,LMN1B,LMN2B))
            DO LMN1A=1,LMNXA
              DO LMN2A=1,LMNXA
                SVAR=DOT_PRODUCT(SPECIAL(ISPA)%ARHO(:,LMN1A,LMN2A),ZA)
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
PRINT*,'EX ',EX
!!$CALL LMTO$REPORTDENMAT(6)
!!$CALL LMTO$REPORTHAMIL(6)
!!$STOP
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
          ALLOCATE(SBARLOC1(N1,N2))
        ELSE
          ALLOCATE(SBARLOC2(N1,N2))
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
                SBARLOC1(:,LMN+IM)=SBAR(NN)%MAT(:,I1+IM)
              ELSE
                SBARLOC2(:,LMN+IM)=SBAR(NN)%MAT(:,I1+IM)
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
          IF(.NOT.POTPAR(ISP)%TORB(LN)) THEN
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
        DT(LMNX1+1:,:LMNX2,I)=-MATMUL(SBARLOC1,D1(:,:,I))
        DT(:LMNX1,LMNX2+1:,I)=-MATMUL(D1(:,:,I),TRANSPOSE(SBARLOC2))
        DT(LMNX1+1:,LMNX2+1:,I)=-MATMUL(SBARLOC1,DT(:LMNX1,LMNX2+1:,I))
      ENDDO
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
      INTEGER(4)             :: I,IAT,ISP,NN,LMN1,LMN2,LMNDOT1,LMNDOT2,LMN3
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
          ALLOCATE(SBARLOC1(N1,N2))
        ELSE
          ALLOCATE(SBARLOC2(N1,N2))
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
                SBARLOC1(:,LMN+IM)=SBAR(NN)%MAT(:,I1+IM)
              ELSE
                SBARLOC2(:,LMN+IM)=SBAR(NN)%MAT(:,I1+IM)
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
     &          -MATMUL(HT(:LMNX1,LMNX2+1:,I),SBARLOC2) &
     &          -MATMUL(TRANSPOSE(SBARLOC1),HT(LMNX1+1:,:LMNX2,I)) &
     &          +MATMUL(TRANSPOSE(SBARLOC1) &
     &                 ,MATMUL(HT(LMNX1+1:,LMNX2+1:,I),SBARLOC2))
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
          IF(.NOT.POTPAR(ISP)%TORB(LN)) THEN
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
      SUBROUTINE LMTO_PLOTTAILED()
!     **************************************************************************
!     **  plots the local orbitals represented by tailed orbitals,            **
!     **  that is using the onsite structure constants and extrapolating      **
!     **  tails.                                                              **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY: k2,ISPECIES,LNX,LOX,POTPAR,SBAR,SBARLI1
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
      real(8)                :: svar
      INTEGER(4)             :: IAT,ISP,LN,L,NN,LMN,IM,I1,IORB,LM,i0
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
        ALLOCATE(SBARLOC(LMNXS,LMNX))
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
              SBARLOC(:,LMN+IM)=SBAR(NN)%MAT(:,I1+IM)
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
WRITE(*,FMT='("SBARLOC",10F10.5)')SBARLOC(:,IORB)
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
      &               -POTPAR(ISP)%TAILED%AEF(:,LNDOT)*SBARLOC(LMNDOT-LMNX,IORB)
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
       &      -POTPAR(ISP)%JBARTOPHIDOT(LN)*SBARLOC(I0+IM,LMN)
            SVAR=SVAR/POTPAR(ISP)%KTOPHI(LN) 
            WRITE(*,FMT='("CPHIDOT:",I5,2F10.5)')LMN,K2,SVAR
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  clean up after iteration                                          ==
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
PRINT*,'EXACT EXCHANGE ENERGY FOR ATOM=',IAT,EX
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
!        CALL LMTO_DOUBLECOUNTING(IAT,NDIMD,NORB,D,EX,H)
        CALL LMTO_SIMPLEDC(GID,NR,NORB,LNX1,LOX(:LNX1,ISP),CHI,LRX,AECORE &
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
      SUBROUTINE LMTO_SIMPLEDC(GID,NR,LMNX,LNX,LOX,CHI,LRX,AECORE &
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
!     ==========================================================================
!     ==  CALCULATE HAMILTONIAN IN TOTAL/SPIN REPRESENTATION                  ==
!     ==========================================================================
      CALL LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
      DEALLOCATE(POT)
!
!     ==========================================================================
!     ==  TRANSFORM HAMILTONIAN FROM TOTAL/SPIN TO UP/DOWN                    ==
!     ==========================================================================
      HAM=REAL(HAM1)
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
      INTEGER(4)            :: IAT1,IAT2,IT(3)
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: IND,INS
      INTEGER(4)            :: NN,MM,I,J,IDIM,IAT,ISP,LN,LMN,IM,I1,I2,I3
      REAL(8)               :: SVAR1,SVAR2
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
      SUBROUTINE LMTO_TESTNTBO()
!     **************************************************************************
!     **  DETERMINES THE PROJECTION FROM THE OCCUPATION OF THE NATURAL        **
!     **  TIGHT-BINDING ORBITALS. IF THEY ARE IDENTICAL TO THE ORIGINAL       **
!     **  PAW PROJECTIONS <PTILDE|PSITILDE>, THE NTB COEFFICIENTS ARE CORRECT **
!     **************************************************************************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV
      USE LMTO_MODULE, ONLY : PERIODICMAT_TYPE,SBAR,ISPECIES,POTPAR &
     &                       ,SBARATOMI1,SBARATOMI2,SBARLI1,LNX,LOX 
      IMPLICIT NONE
      REAL(8)       :: XK(3,NKPTL)
      INTEGER(4)    :: NPRO
      INTEGER(4)    :: NBH
      INTEGER(4)    :: IKPT,ISPIN,IBH,IPRO,IDIM
      COMPLEX(8),ALLOCATABLE :: TBC1(:,:)
      COMPLEX(8),ALLOCATABLE :: VEC1(:,:)
      COMPLEX(8),ALLOCATABLE :: VEC2(:,:)
      COMPLEX(8),ALLOCATABLE :: SBAROFK(:,:)
      INTEGER(4)             :: NSMALL,NNS,NAT
      INTEGER(4)             :: IAT,ISP,LN,L,IRL,IM,LN1
      INTEGER(4)             :: NL   !#(PROJECTORS IN THE SAME LM-CHANNEL)
      REAL(8)                :: SVAR
      REAL(8)                :: DEVMAX
!     **************************************************************************
      DEVMAX=0.D0
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      NPRO=MAP%NPRO
!
!     ==========================================================================
!     ==  CHECK PROJECTIONS                                                   ==
!     ==========================================================================
      IF(.NOT.ASSOCIATED(THIS%TBC)) THEN
        CALL ERROR$MSG('THIS%TBC IS NOT ASSOCIATED')
        CALL ERROR$STOP('LMTO_TESTNTBO')
      END IF
      NAT=SIZE(ISPECIES)
      NSMALL=MAXVAL(SBARATOMI2)
      NNS=SIZE(SBAR)
      ALLOCATE(SBAROFK(NSMALL,NSMALL))
      ALLOCATE(TBC1(NDIM,NPRO))
      ALLOCATE(VEC1(NDIM,NSMALL))
      ALLOCATE(VEC2(NDIM,NSMALL))
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          CALL LMTO_AOFK(NNS,SBAR,XK(:,IKPT),NSMALL,SBAROFK)
          NBH=THIS%NBH
          DO IBH=1,NBH
!!$PRINT*,'MARKE 1',IKPT,ISPIN,IBH
!!$PRINT*,'MARKE 1A',SIZE(THIS%TBC)
!!$PRINT*,'MARKE 1B',SIZE(TBC1)
!!$PRINT*,'MARKE 1C',SIZE(VEC1)
!!$PRINT*,'MARKE 1D',SIZE(VEC2)
            TBC1(:,:)=THIS%TBC(:,IBH,:)            
            VEC1(:,:)=(0.D0,0.D0)
            VEC2(:,:)=(0.D0,0.D0)
            IPRO=0
            DO IAT=1,NAT
              ISP=ISPECIES(IAT)
              DO LN=1,LNX(ISP)
                L=LOX(LN,ISP)
                IRL=SBARATOMI1(IAT)-1+SBARLI1(L+1,ISP)-1
                DO IM=1,2*L+1
                  IPRO=IPRO+1
                  IRL=IRL+1
                  VEC1(:,IRL)=VEC1(:,IRL)+TBC1(:,IPRO)
                  VEC2(:,IRL)=VEC2(:,IRL)+POTPAR(ISP)%KTOPHIDOT(LN)*TBC1(:,IPRO)
                ENDDO
              ENDDO
            ENDDO
!
!           == MULTIPLY WITH SBAR ==============================================
            DO IDIM=1,NDIM
              VEC1(IDIM,:)=MATMUL(SBAROFK,VEC1(IDIM,:))
            ENDDO
!
            IPRO=0
            DO IAT=1,NAT
              ISP=ISPECIES(IAT)
              DO LN=1,LNX(ISP)
                L=LOX(LN,ISP)
                NL=0 
                DO LN1=1,LNX(ISP)
                  IF(LOX(LN1,ISP).EQ.L)NL=NL+1
                ENDDO
                IRL=SBARATOMI1(IAT)-1+SBARLI1(L+1,ISP)-1
                DO IM=1,2*L+1
                  IPRO=IPRO+1
                  IRL=IRL+1
!                 == CONTRIBUTION FROM KBAR
                  TBC1(:,IPRO)=TBC1(:,IPRO)*POTPAR(ISP)%KTOPHI(LN)
                  TBC1(:,IPRO)=TBC1(:,IPRO)+VEC2(:,IRL)*POTPAR(ISP)%PHIDOTPROJ(LN)
                  SVAR=POTPAR(ISP)%PHIDOTPROJ(LN)*POTPAR(ISP)%JBARTOPHIDOT(LN)
                  SVAR=SVAR*REAL(NL)
                  TBC1(:,IPRO)=TBC1(:,IPRO)-SVAR*VEC1(:,IRL)
                ENDDO
              ENDDO
            ENDDO
!!$WRITE(*,FMT='("PI",I5,180F10.5)')2*IBH-1,REAL(THIS%PROJ(1,IBH,:))
!!$WRITE(*,FMT='("PF",I5,180F10.5)')2*IBH-1,REAL(TBC1(1,:))
!!$WRITE(*,*)
!!$WRITE(*,FMT='("PI",I5,180F10.5)')2*IBH,AIMAG(THIS%PROJ(1,IBH,:))
!!$WRITE(*,FMT='("PF",I5,180F10.5)')2*IBH,AIMAG(TBC1(1,:))
!!$WRITE(*,*)
!           == CALCULATE DEVIATION (SHOULD BE ZERO) =============================
            TBC1(:,:)=TBC1(:,:)-THIS%PROJ(:,IBH,:)
            DEVMAX=MAX(DEVMAX,MAXVAL(REAL(TBC1(:,:))),MAXVAL(AIMAG(TBC1(:,:))))
          ENDDO
        ENDDO
      ENDDO
      WRITE(*,FMT='("TEST OF NTBO COEFFICIENTS. MAX. DEV.=",E20.5)')DEVMAX
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
!!$      CALL SETUP$GETR8('AEZ',AEZ)
!!$      WRITE(STRING,FMT='(F3.0)')AEZ
!!$      STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
!!$      CALL SETUP_WRITEPHI(-'CHI'//TRIM(STRING),GID,NR,LNCHI,CHI)
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
      REAL(8)                :: QBAR
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: LN,LN1,L,I,J,IIB,LM,LNCHI,IR,IM
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
      CALL SETUP$GETR8('AEZ',AEZ)
      WRITE(STRING,FMT='(F3.0)')AEZ
      STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
      CALL SETUP_WRITEPHI(-'CHI'//TRIM(STRING),GID,NR,LNCHI,CHI)
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
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTOVERLAP')
      NNS=SIZE(OVERLAP)
      WRITE(NFIL,FMT='(82("="))')
      WRITE(NFIL,FMT='(82("="),T10,"   OVERLAP MATRIX ELEMENTS   ")')
      WRITE(NFIL,FMT='(82("="))')
      DO NN=1,NNS
        WRITE(NFIL,FMT='(82("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 OVERLAP(NN)%IAT1,OVERLAP(NN)%IAT2,OVERLAP(NN)%IT
        DO LM1=1,OVERLAP(NN)%N1
          WRITE(NFIL,FMT='(20E10.2)')OVERLAP(NN)%MAT(:,LM1)
        ENDDO
        WRITE(NFIL,FMT='(82("="))')
      ENDDO
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
          WRITE(NFIL,FMT='(20F10.5)')SBAR(NN)%MAT(:,LM1)
        ENDDO
        WRITE(NFIL,FMT='(82("="))')
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END SUBROUTINE LMTO$REPORTSBAR
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
                                              CALL TRACE$PUSH('LMTO_PLOTLOCORB')
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
      call LMTO_LMTO_GRIDPLOT_TAILEDINNER(NAT,R0,IAT0,NP,P,NORB,ORB)
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
      REAL(8)   ,ALLOCATABLE:: r(:)
      REAL(8)               :: RX
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      REAL(8)               :: DR(3),DRLEN
      REAL(8)               :: VAL
      INTEGER(4)            :: NN,L,IM,IP,LN,LM,LMN,I0
!     **************************************************************************
      ISP=ISPECIES(IAT1)
      CALL SETUP$ISELECT(ISP)
      GID=POTPAR(ISP)%TAILED%GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      allocate(r(nr))
      CALL RADIAL$r(GID,NR,R)
      rx=r(nr)
      deallocate(r)
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
      ALLOCATE(SBARLOC(LMNXS,LMNX))
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
            SBARLOC(:,LMN+IM)=SBAR(NN)%MAT(:,I0+IM)
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
            ORB(IP,:)=ORB(IP,:)+VAL*YLM(LM)*SBARLOC(I0+IM,:)
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
      INTEGER(4)            :: LX
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
      CHARACTER(16)         :: FORMATTYPE
      REAL(8)               :: X1D(NRAD)
      REAL(8)               :: VAL
      INTEGER(4)            :: NP,IP
      REAL(8)  ,ALLOCATABLE :: P(:,:)
      LOGICAL(4),PARAMETER  :: TGAUSS=.TRUE.
      CHARACTER(2),PARAMETER :: ID='3D'
      INTEGER(4)            :: NN,IAT,IAT2,ISP,I1,I2,I3,LN,I,J,IORB,IDIR
!     **************************************************************************
                                              CALL TRACE$PUSH('LMTO_PLOTLOCORB')
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
      INTEGER(4)            :: LM2X
      INTEGER(4),PARAMETER  :: LMXX=36
      REAL(8)               :: K0(LMXX)
      REAL(8)               :: J0(LMXX)
      REAL(8)               :: JBAR(LMXX)
      REAL(8)   ,ALLOCATABLE:: CVEC(:,:)  !SCREENPARM
      REAL(8)               :: CVECSUM(LMXX)
      REAL(8)               :: R2(3) !  POSITIONS OF NEIGBORS
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
      CHARACTER(16)         :: FORMATTYPE
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
!CALL TEST_LMTO$STRUCTURECONSTANTS()
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
      &                           -DOT_PRODUCT(JBAR(1:LM2X),SBAR(NN)%MAT(:,LM1))
              ORB(IP,LM1)=ORB(IP,LM1) &
      &                        -DOT_PRODUCT(JBARAUG(1:LM2X),SBAR(NN)%MAT(:,LM1))
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
           CVEC(:LM2X,LM1)=QBARVEC(:LM2X,IAT2)*SBAR(NN)%MAT(:,LM1)
           IF(TONSITE)CVEC(LM1,LM1)=CVEC(LM1,LM1)+1.D0
         ENDDO
!        == LOOP OVER REAL SPACE GRID ==========================================
         DO IP=1,NP
!          == DR IS THE DISTANCE FROM THE SECOND ATOM ==========================
           DR(:)=P(:,IP)-R2(:)
           TSPHERE=(DOT_PRODUCT(DR,DR).LT.RAD(IAT2))
!
!          == DETERMINE BARE HANKEL AND SCREENED BESSEL FUNCTION AT IAT2 =======
           CALL  LMTO$SOLIDHANKEL(DR,RAD(IAT2),K2,LM2X,K0(1:LM2X))
           IF(TSPHERE) THEN
             CALL  LMTO$SOLIDBESSEL(DR,K2,LM2X,J0(1:LM2X))
             JBAR(:LM2X)=J0(:LM2X)-K0(:LM2X)*QBARVEC(:LM2X,IAT2)
           END IF
!          ==  
           DO LM1=1,LM1X
!            == |KBAR>=|K0>*(1+QBAR*SBAR) ======================================
             ORB(IP,LM1)=ORB(IP,LM1)+DOT_PRODUCT(K0(1:LM2X),CVEC(1:LM2X,LM1))
!            == -|JBAR>SBAR=-(|J0>-|K0>QBAR)*SBAR ==============================
             IF(TSPHERE) THEN
               ORBI(IP,LM1)=ORBI(IP,LM1) &
      &                           -DOT_PRODUCT(JBAR(1:LM2X),SBAR(NN)%MAT(:,LM1))
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
!     ** UNITS WRITTEN ARE ABOHR, CONSISTENT WITH AVOGADRO'S IMPLEMENTATION.  **
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
!     ** wrapper called from paw_graphics to plot a wave function in ntbo     **
!     ** representation                                                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN)  :: NFIL
      INTEGER(4)   ,INTENT(IN)  :: IDIM0
      INTEGER(4)   ,INTENT(IN)  :: IB0
      INTEGER(4)   ,INTENT(IN)  :: IKPT0
      INTEGER(4)   ,INTENT(IN)  :: ISPIN0
      INTEGER(4)   ,INTENT(IN)  :: NR1,NR2,NR3
!     **************************************************************************
      call LMTO$PLOTWAVE_tailed(NFIL,IDIM0,IB0,IKPT0,ISPIN0,NR1,NR2,NR3)
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$PLOTWAVE_tailed(NFIL,IDIM0,IB0,IKPT0,ISPIN0,NR1,NR2,NR3)
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
!     == coefficients for the ntbo basis                                      ==
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
!     == map wace function onto grid                                          ==
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
