MODULE LMTO_MODULE
LOGICAL(4)            :: TON=.false.
REAL(8)   ,PARAMETER  :: K2=-0.01D0    ! 0.5*K2 IS THE KINETIC ENERGY
!REAL(8)   ,PARAMETER  :: RCSCALE=2.D0  !RADIUS SCALE FACTOR FOR NEIGHBORLIST
REAL(8)   ,PARAMETER  :: RCSCALE=1.2D0  !RADIUS SCALE FACTOR FOR NEIGHBORLIST
!         rcscale times the sum of covalent radii defines cutoff for neigborlist
!== PARAMETER SET TO DEFINE THE GAUSSIANS FOR THE UNSCREENED HANKEL FUNCTIONS ==
!== IS USED IN LMTO_GAUSSFITKPRIME  ============================================
INTEGER(4),PARAMETER  :: GAUSSFITKPRIME_NPOW=2   ! -1 INDICATES LX
INTEGER(4),PARAMETER  :: GAUSSFITKPRIME_NE=12
REAL(8)   ,PARAMETER  :: GAUSSFITKPRIME_R1=0.6667D0
REAL(8)   ,PARAMETER  :: GAUSSFITKPRIME_SCALER=1.25D0
!== PARAMETER SET TO DEFINE THE GAUSSIANS FOR THE SCREENED HANKEL FUNCTIONS   ==
INTEGER(4),PARAMETER  :: GAUSSFITKBARPRIME_NPOW=4
INTEGER(4),PARAMETER  :: GAUSSFITKBARPRIME_NE=4
REAL(8)   ,PARAMETER  :: GAUSSFITKBARPRIME_R1=1.D0
REAL(8)   ,PARAMETER  :: GAUSSFITKBARPRIME_SCALER=1.5D0
!== PARAMETER SET TO DEFINE THE GAUSSIANS FOR THE AUGMENTATION
INTEGER(4),PARAMETER  :: GAUSSFITKAUGMENT_NPOW=4
INTEGER(4),PARAMETER  :: GAUSSFITKAUGMENT_NE=6
REAL(8)   ,PARAMETER  :: GAUSSFITKAUGMENT_R1=3.D-2
REAL(8)   ,PARAMETER  :: GAUSSFITKAUGMENT_SCALER=2.D0
!
TYPE TAILED_type
  INTEGER(4)         :: GID
  INTEGER(4)         :: lnx
  INTEGER(4)         :: lmnx
  INTEGER(4),pointer :: lox(:)     ! (lnx)
  INTEGER(4),pointer :: lndot(:)   ! (lnx)
  INTEGER(4),pointer :: lmndot(:)  ! (lmnx)
  REAL(8)   ,POINTER :: aef(:,:)   ! (NR,Lnx)
  REAL(8)   ,POINTER :: psf(:,:)   ! (NR,Lnx)
  REAL(8)   ,POINTER :: nlf(:,:)   ! (NR,Lnx)
  real(8)   ,pointer :: u(:,:,:,:) ! (lmnx,lmnx,lmnx,lmnx)
END TYPE tailed_type

TYPE ORBITALSPHHARM_TYPE
  INTEGER(4)         :: GID
  INTEGER(4)         :: NR
  INTEGER(4)         :: LMX
  INTEGER(4)         :: NORB
  REAL(8)   ,POINTER :: F(:,:,:) !(NR,LM,IORB)
END TYPE ORBITALSPHHARM_TYPE
!
TYPE ORBITALGAUSSCOEFF_TYPE
  INTEGER(4)         :: NIJK     
  INTEGER(4)         :: NE
  INTEGER(4)         :: NORB
  REAL(8)   ,POINTER :: E(:)     !(NE)
  REAL(8)   ,POINTER :: C(:,:,:) !(NIJK,NE,NORB)
END TYPE ORBITALGAUSSCOEFF_TYPE
!== POTPARRED CONSIDERS ONLY ONE ANGULAR MOMENTUM CHANNEL PER LM ===============
!== CONSISTENT WITH THE SCREENED STRUCTURE CONSTANTS ===========================
TYPE POTPARRED_TYPE
  REAL(8)   ,POINTER :: DOVERLAPKK(:,:)
  REAL(8)   ,POINTER :: DOVERLAPKJ(:,:)
  REAL(8)   ,POINTER :: DOVERLAPJJ(:,:)
END TYPE POTPARRED_TYPE
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
  REAL(8)   ,POINTER :: DOVERLAPKK(:,:) !(LNX,LNX)
  REAL(8)   ,POINTER :: DOVERLAPKJ(:,:) !(LNX,LNX)
  REAL(8)   ,POINTER :: DOVERLAPJJ(:,:) !(LNX,LNX)
  TYPE(tailed_TYPE)            :: tailed
  TYPE(POTPARRED_TYPE)         :: SMALL
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSKPRIME
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSTAILEDK
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSTAILEDJBAR
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSKAUGMENT
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSJAUGMENT
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
INTEGER(4),ALLOCATABLE  :: LNX(:)              !(NSP)
INTEGER(4),ALLOCATABLE  :: LOX(:,:)            !(LnxX,NSP)
INTEGER(4),ALLOCATABLE  :: ISPECIES(:)         !(NAT)
REAL(8)   ,ALLOCATABLE  :: ORBRAD(:,:) !(LXX+1,NAT) NODE-POSITION OF THE ORBITAL
TYPE(POTPAR_TYPE)     ,ALLOCATABLE :: POTPAR(:)!POTENTIAL PARAMETERS
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: SBAR(:)  !(NNS) SCREENED STRUCTURE CONST.
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DENMAT(:)  !(NNS) DENSITY MATRIX
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: HAMIL(:)  !(NNS) DERIVATIVE OF ENERGY
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: OVERLAP(:) !(NNS) OVERLAP MATRIX ONLY MAIN CHANNEL
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DENMAT_T(:) !(NNS) DENSITY MATRIX
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: HAMIL_T(:)  !(NNS) DERIVATIVE OF ENERGY
INTEGER(4)                         :: NNUX    !DIMENSION OF UMAT
INTEGER(4)                         :: NNU     !#(ELEMENTS OF UMAT)
TYPE(UMAT_TYPE)       ,ALLOCATABLE :: UMAT(:) !(NNUX/NNU) U-MATRIX ELEMENTS
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
logical(4),parameter :: tspherical=.false.

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
      USE LMTO_MODULE
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
      CALL LMTO_SBARINDICES2()   
!
!     ==========================================================================
!     == DETERMINE POTENTIAL PARAMETERS                                       ==
!     == RAD,LNSCATT,PHIDOTPROJ,QBAR,KTOPHI,KTOPHIDOT,JBARTOPHIDOT            ==
!     ==========================================================================
      CALL LMTO_MAKEPOTPAR()
!
!     ==========================================================================
!     == attach exponential tails to augmented Hankel and Bessel functions    ==
!     ==========================================================================
      CALL LMTO_MAKETAILEDPARTIALWAVES()
!
!     ==========================================================================
!     == DETERMINE GAUSS EXPANSION OF UNSCREENED HANKEL FUNCTIONS KPRIME      ==
!     ==========================================================================
      CALL LMTO_GAUSSFITKPRIME()
      CALL LMTO_GAUSSFITKAUGMENT()
!
!     ==========================================================================
!     == DETERMINE KPRIME AND JBAR WITH EXPONENTIAL TAILS                     ==
!     == USED TO BUILD UP APPROXIMATE NTBO'S IN A ONE-CENTER EXPANSION        ==
!     ==========================================================================
      CALL LMTO_GAUSSFITKJTAILS()
      CALL LMTO_ONSITEOVERLAP()
!
!     ==========================================================================
!     == SET UP INDEX ARRAY FOR STRUCTURE CONSTANTS                           ==
!     ==========================================================================
      CALL LMTO_SBARINDICES()
      NAT=SIZE(ISPECIES)
      NRL=SBARATOMI2(NAT)  !#(TIGHT-BINDING ORBITALS)
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
      USE LMTO_MODULE
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
!       ==  collect data                                                      ==
!       ========================================================================
!       == radial grid =========================================================
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
!       == matching radius =====================================================
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(NINT(AEZ),'R(ASA)',RAD) 
        POTPAR(ISP)%RAD=RAD
!       == partial waves and projectors ========================================
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
!       ==  select one phibardot function per l                               ==
!       ========================================================================
!       == GET INFO ON SCATTERING CHANNELS =====================================
!       == ISCATT=0 FOR VALENCE, ISCATT>0 FOR SCATTERING, ISCATT<0 FORSEMICORE =
        ALLOCATE(ISCATT(LNX1))
        CALL SETUP$GETI4A('ISCATT',LNX1,ISCATT)
!       == DETERIMINE SELECTION MAP ============================================
        ALLOCATE(POTPAR(ISP)%LNSCATT(LNX1))
        POTPAR(ISP)%LNSCATT(:)=-1
        DO LN=1,LNX1
          IF(ISCATT(LN).NE.0) CYCLE
          DO LN1=1,LNX1
            IF(LOX(LN1,ISP).NE.LOX(LN,ISP)) CYCLE
            POTPAR(ISP)%LNSCATT(LN1)=LN
          ENDDO
        ENDDO
        DO LN=1,LNX1
          IF(POTPAR(ISP)%LNSCATT(LN).EQ.-1) THEN
            CALL ERROR$MSG('SELECTION OF ENERGY DERIVATIVE PARTIAL WAVE FAILED')
            CALL ERROR$STOP('LMTO_MAKEPOTPAR')
          END IF
        ENDDO
        DEALLOCATE(ISCATT)
!       == MAP =================================================================
        DO LN=1,LNX1
          LN1=POTPAR(ISP)%LNSCATT(ln)
          if(ln1.eq.ln) cycle
          AEPHIDOT(:,LN)=AEPHIDOT(:,LN1)
          PSPHIDOT(:,LN)=PSPHIDOT(:,LN1)
          NLPHIDOT(:,LN)=NLPHIDOT(:,LN1)
        ENDDO
!
!       ========================================================================
!       ==  determine potential parameters                                    ==
!       ========================================================================
        ALLOCATE(POTPAR(ISP)%QBAR(LNX1))
        ALLOCATE(POTPAR(ISP)%PHIDOTPROJ(LNX1))
        ALLOCATE(POTPAR(ISP)%KTOPHI(LNX1))
        ALLOCATE(POTPAR(ISP)%KTOPHIDOT(LNX1))
        ALLOCATE(POTPAR(ISP)%JBARTOPHIDOT(LNX1))
        DO LN=1,LNX(ISP)
          l=LOX(LN,ISP)
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
      USE LMTO_MODULE, ONLY : K2,POTPAR,NSP,LNX,LOX
      IMPLICIT NONE
      REAL(8)   ,PARAMETER   :: LAMBDA1=2.D0
      REAL(8)   ,PARAMETER   :: LAMBDA2=1.D0
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)                :: RAD
      integer(4)             :: irad ! first point beyond rad
      REAL(8)                :: QBAR
      REAL(8)                :: JVAL,JDER,KVAL,KDER
      REAL(8)                :: SVAR1,SVAR2,A1,A2,B1,B2
      INTEGER(4)             :: L
      INTEGER(4)             :: lrx,lmrx
      INTEGER(4)             :: lnxt
      INTEGER(4)             :: lmnxt
      INTEGER(4),allocatable :: loxt(:)
      INTEGER(4),allocatable :: lndot(:)
      INTEGER(4),allocatable :: lmndot(:)
      INTEGER(4)             :: ISP,LN,ln1,ln2,lnt,lmn,lmn1,lmn2,im,IR
      real(8)   ,allocatable :: aephi(:,:)
      real(8)   ,allocatable :: aephidot(:,:)
      real(8)   ,allocatable :: psphi(:,:)
      real(8)   ,allocatable :: psphidot(:,:)
      real(8)   ,allocatable :: nlphi(:,:)
      real(8)   ,allocatable :: nlphidot(:,:)
      real(8)   ,allocatable :: ulittle(:,:,:,:,:)
character(128) :: string
!     **************************************************************************
                              CALL TRACE$PUSH('LMTO_MAKETAILEDPARTIALWAVES')
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
!       == RADIAL GRID =========================================================
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        RAD=POTPAR(ISP)%RAD
        do ir=1,nr
          irad=ir
          if(r(ir).gt.rad) exit
        enddo
print*,'irad1= ',irad,nr,rad,r(nr)
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
        ALLOCATE(loxt(lnxt))
        ALLOCATE(lndot(lnxt))
        ALLOCATE(lmndot(lmnxt))
        LNT=LNX(ISP)
        DO LN=1,LNX(ISP)
          LOXt(:LN)=LOX(:LN,ISP)  ! for the hankel functions
          IF(POTPAR(ISP)%LNSCATT(LN).NE.LN) CYCLE
          LNT=LNT+1
          LOXT(LNT)=LOX(LN,ISP)   ! for the bessel functions
          LNDOT(LN)=LNT
          LNDOT(LNT)=LN
        ENDDO
        DO LN=1,LNX(ISP)          ! COMPLETE LNDOT ARRAY
          LNDOT(LN)=LNDOT(POTPAR(ISP)%LNSCATT(LN))
        ENDDO 
        ALLOCATE(POTPAR(ISP)%TAILED%LOX(LNXT))
        ALLOCATE(POTPAR(ISP)%TAILED%LNDOT(LNXT))
        ALLOCATE(POTPAR(ISP)%TAILED%LMNDOT(LMNXT))
        POTPAR(ISP)%TAILED%LOX(:)=LOXT
        POTPAR(ISP)%TAILED%Lndot=Lndot
        POTPAR(ISP)%TAILED%Lmndot=Lmndot
!
!       ========================================================================
!       == MAPPING "LMNDOT" FROM K TO JBAR FUNCTIONS AND VICE VERSA           ==
!       ==   JBAR(LMN)=F(LMNDOT(LMN)) ;                                       ==
!       ========================================================================
        lmn1=0
        LMN2=SUM(2*LOX(:LNx(ISP),ISP)+1)
        LN2=LNx(ISP)
        DO LN1=1,LNX(ISP)
          IF(POTPAR(ISP)%LNSCATT(LN1).EQ.LN1) THEN
            LN2=LN2+1
            DO IM=1,2*LOX(LN1,ISP)+1
              LMNDOT(LMN2+IM)=LMNDOT(LMN1+IM)
            ENDDO
            LMN=0
            DO LN=1,LNX(ISP)
              IF(POTPAR(ISP)%LNSCATT(LN).EQ.LN1) THEN
                DO IM=1,2*LOX(LN,ISP)+1
                  LMNDOT(LMN+IM)=LMNDOT(LMN2+IM)
                ENDDO
              END IF
              LMN=LMN+2*LOX(LN,ISP)+1
            ENDDO
            LMN2=LMN2+2*Lox(ln1,isp)+1
          END IF
          LMN1=LMN1+2*Lox(ln1,isp)+1
        ENDDO
!
!       ========================================================================
!       == augmented hankel and bessel functions with tails attached          ==
!       ========================================================================
        ALLOCATE(POTPAR(ISP)%TAILED%aeF(NR,LMNXT))
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
        CALL SETUP$GETR8A('NLPHI',NR*LNX(ISP),NLPHI)
        CALL SETUP$GETR8A('NLPHIDOT',NR*LNX(ISP),NLPHIDOT)
        CALL SETUP$GETR8A('PSPHI',NR*LNX(ISP),PSPHI)
        CALL SETUP$GETR8A('PSPHIDOT',NR*LNX(ISP),PSPHIDOT)
!
!       == tail part ===========================================================
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
          DO IR=irad,NR
            SVAR1=EXP(-LAMBDA1*(R(IR)-RAD))
            SVAR2=EXP(-LAMBDA2*(R(IR)-RAD))
            POTPAR(ISP)%TAILED%nlF(IR,LN)       =A1*SVAR1+A2*SVAR2
            POTPAR(ISP)%TAILED%nlF(IR,LNDOT(LN))=B1*SVAR1+B2*SVAR2
          ENDDO
        ENDDO
!
!       == sphere part =========================================================
        DO LN=1,LNX(ISP)
          l=lox(ln,isp)
          CALL LMTO$SOLIDHANKELRAD(L,Rad,K2,KVAL,KDER)
          CALL LMTO$SOLIDBESSELRAD(L,Rad,K2,JVAL,JDER)
          LN1=POTPAR(ISP)%LNSCATT(LN)
          A1=POTPAR(ISP)%KTOPHI(LN)
          A2=POTPAR(ISP)%KTOPHIDOT(LN)
          POTPAR(ISP)%TAILED%nlF(:irad-1,LN)=nlPHI(:irad-1,LN)*A1 &
       &                                    +nlPHIDOT(:irad-1,LN1)*A2
!         == the complex addition of differences in the following is ===========
!         == necessary, if the ae, ps and nl partial waves differ at the =======
!         == matching radius due to the admixed tails of core states ===========
          POTPAR(ISP)%TAILED%aeF(:,LN)=POTPAR(ISP)%TAILED%nlF(:,LN) &
       &      +(aePHI(:,LN)-nlphi(:,ln))*A1+(aephidot(:,ln1)-nlPHIDOT(:,LN1))*A2
          POTPAR(ISP)%TAILED%PSF(:,LN)=POTPAR(ISP)%TAILED%NLF(:,LN) &
          &   +(PSPHI(:,LN)-NLPHI(:,LN))*A1+(PSPHIDOT(:,LN1)-NLPHIDOT(:,LN1))*A2
!
          IF(LN.EQ.LN1) THEN
            A2=POTPAR(ISP)%JBARTOPHIDOT(LN)
            POTPAR(ISP)%TAILED%NLF(:IRAD-1,LNDOT(LN))=NLPHIDOT(:IRAD-1,LN1)*A2
!           == the complex addition of differences in the following is =========
!           == necessary, if the ae, ps and nl partial waves differ at the =====
!           == matching radius due to the admixed tails of core states =========
            POTPAR(ISP)%TAILED%AEF(:,LNDOT(LN)) &
      &                           =POTPAR(ISP)%TAILED%NLF(:,LNDOT(LN)) &
      &                           +(AEPHIDOT(:,LN1)-NLPHIDOT(:,LN1))*A2
            POTPAR(ISP)%TAILED%psF(:,LNDOT(LN)) &
      &                           =POTPAR(ISP)%TAILED%NLF(:,LNDOT(LN)) &
      &                           +(psPHIDOT(:,LN1)-NLPHIDOT(:,LN1))*A2
!
          END IF
        ENDDO
        deallocate(aephi)
        deallocate(aephidot)
        deallocate(psphi)
        deallocate(psphidot)
        deallocate(nlphi)
        deallocate(nlphidot)
        DEALLOCATE(R)
!
WRITE(STRING,FMT='(I5)')ISP
STRING='TAILS_FORATOMTYPE'//TRIM(ADJUSTL(STRING))//'.DAT'
CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LNXT,POTPAR(ISP)%TAILED%aeF)
!
!       ========================================================================
!       == onsite U-TENSOR OF TAILED PARTIAL WAVES                            ==
!       ========================================================================
        CALL SETUP$GETI4('LMRX',LMRX)
        LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
        ALLOCATE(POTPAR(ISP)%TAILED%U(LMNXT,LMNXT,LMNXT,LMNXT))
        ALLOCATE(ULITTLE(LRX+1,LNXT,LNXT,LNXT,LNXT))
        CALL LDAPLUSU_ULITTLE(GID,NR,LRX,LNXT,LOXT,POTPAR(ISP)%TAILED%aeF,ULITTLE)
        CALL LDAPLUSU_UTENSOR(LRX,LMNXT,LNXT,LOXt,ULITTLE,POTPAR(ISP)%TAILED%U)
        DEALLOCATE(ULITTLE)
!
        DEALLOCATE(lndot)
        DEALLOCATE(lmndot)
        DEALLOCATE(loxt)
      ENDDO
                              call trace$pop()
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE LMTO_tailedproducts()
!!$!     **************************************************************************
!!$!     **                                                                      **
!!$!     ******************************PETER BLOECHL, GOSLAR 2011******************
!!$      USE LMTO_MODULE, ONLY : POTPAR,NSP
!!$      IMPLICIT NONE
!!$!     **************************************************************************
!!$      do isp=1,nsp
!!$        lnx=potpar(isp)%tailed%lnx
!!$        allocate(lox(lnx))
!!$        lox=potpar(isp)%tailed%lox
!!$        gid=potpar(isp)%tailed%gid
!!$        call radial$geti4(gid,'nr',nr)
!!$        allocate(aux(nr))
!!$        allocate(w(nr)) ! fitting weighting function
!!$        w(:)=1.d0
!!$!
!!$        
!!$        potpar%tailed%products%ne=
!!$        do ie=1,ne
!!$          potpar%tailed%products%e(ie)=1.d0/(r1*rfac**(ie-1))
!!$        enddo
!!$
!!$        ip=0
!!$        do ln1=1,lnx
!!$          l1=lox(ln1)
!!$          do ln2=ln1,lnx
!!$            l2=lox(ln2)
!!$            aux=potpar(isp)%tailed%aef(:,ln1)*potpar(isp)%tailed%aef(:,ln2)
!!$            do l3=min(l1,l2),l1+l2,2  ! triangle rule
!!$              ip=ip+1
!!$              CALL GAUSSIAN_FITGAUSS(GID,NR,W,L3,aux,NE,NPOW2+1,E(:NE) &
!!$       &                                                    ,CPOWK(:NPOW2+1,:))
!!$
!!$            enddo
!!$          enddo
!!$        enddo
!!$      enddo
!!$      return
!!$      end
!!$
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_KJBARTAILS(GID,NR,L,RAD,K2,QBAR,LAMBDA,KF,JBARF)
!     **************************************************************************
!     **  CONSTRUCTS EXPONENTIAL TAILS FOR HANKEL FUNCTION AND                **
!     **  SCREENED BESSELFUNCTIONS                                            **
!     **  FIT (A+BR^2)*R^L*EXP(-LAMBDA*R) AT THE RADIUS RAD                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: RAD
      REAL(8)   ,INTENT(IN) :: K2
      REAL(8)   ,INTENT(IN) :: QBAR
      REAL(8)   ,INTENT(IN) :: LAMBDA
      REAL(8)   ,INTENT(OUT):: KF(NR)
      REAL(8)   ,INTENT(OUT):: JBARF(NR)
      REAL(8)               :: KVAL,KDER,JVAL,JDER
      REAL(8)               :: SVAR,SVAR1,SVAR2
      REAL(8)               :: AJ,BJ,AK,BK
      REAL(8)               :: R(NR)
      REAL(8)               :: X
      INTEGER(4)            :: IR
!     **************************************************************************
      CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
      CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
      JVAL=JVAL-QBAR*KVAL
      JDER=JDER-QBAR*KDER
      SVAR=0.5D0*((REAL(L)-LAMBDA*RAD)*JVAL-RAD*JDER)
      AJ=SVAR+JVAL
      BJ=-SVAR
      SVAR=0.5D0*((REAL(L)-LAMBDA*RAD)*KVAL-RAD*KDER)
      AK=SVAR+KVAL
      BK=-SVAR
!
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        X=R(IR)/RAD
        SVAR1=X**L*EXP(-LAMBDA*(R(IR)-RAD))
        SVAR2=SVAR1*X**2
        KF(IR)=AK*SVAR1+BK*SVAR2
        JBARF(IR)=AJ*SVAR1+BJ*SVAR2
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GAUSSFITKAUGMENT()
!     **************************************************************************
!     ** CONSTRUCTS THE AUGMENTATION TO THE UNSCREENED HANKEL FUNCTIONS       **
!     ** AND DETERMINES ITS COEFFICIENTS IN A GAUSSIAN EXPANSION.             **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      USE LMTO_MODULE, ONLY : K2,POTPAR,NSP,LNX,LOX,SBARLI1 &
     &                       ,GAUSSFITKAUGMENT_NPOW,GAUSSFITKAUGMENT_NE &
     &                       ,GAUSSFITKAUGMENT_R1,GAUSSFITKAUGMENT_SCALER
      LOGICAL(4)             :: TTEST=.false.
      INTEGER(4),PARAMETER   :: NEX=30
      REAL(8)   ,PARAMETER   :: LAMBDA=1.D0
      REAL(8)                :: E(10)
      INTEGER(4)             :: NPOW
      INTEGER(4)             :: NPOW2
      INTEGER(4)             :: NE
      INTEGER(4)             :: NIJK
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      INTEGER(4)             :: L
      INTEGER(4)             :: LX
      INTEGER(4)             :: LM
      INTEGER(4)             :: LNX1
      INTEGER(4)             :: NORBX
      INTEGER(4)             :: NORBS
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: W(:)
      REAL(8)   ,ALLOCATABLE :: KPRIME(:)
      REAL(8)   ,ALLOCATABLE :: JBAR(:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: CPOWK(:,:)    !(NPOW2,NE)
      REAL(8)   ,ALLOCATABLE :: CPOWJ(:,:)    !(NPOW2,NE)
      REAL(8)   ,ALLOCATABLE :: CK(:,:)       !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: CJ(:,:)       !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: WORKK(:,:)    !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: WORKJ(:,:)    !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: YLMPOL(:,:)
      REAL(8)                :: RAD
      REAL(8)                :: QBAR
      REAL(8)                :: SVAR,SVAR1,SVAR2,JVAL,KVAL,JDER,KDER
      INTEGER(4)             :: ISP,IE,I,J,K,I1,J1,K1,LN,IM,IR,N,IND,IND1,IND2
      INTEGER(4)             :: IND1X,IORB,IORBS
      INTEGER(4)             :: NIJKXPOW2
      INTEGER(4)             :: LN1
      CHARACTER(128)         :: FILE
      CHARACTER(8)           :: STRING1,STRING2
      REAL(8)                :: R0(3)
      INTEGER(4)             :: NFIL
      LOGICAL(4)             :: TSCATT
!     **************************************************************************
                                  CALL TRACE$PUSH('LMTO_GAUSSFITKAUGMENT')
!
!     ==========================================================================
!     == POLYNOMIAL COEFFICIENTS OF SPHERICAL HARMONICS TIMES R**L            ==
!     ==========================================================================
      LX=MAXVAL(LOX)
      ALLOCATE(YLMPOL((LX+1)*(LX+2)*(LX+3)/6,(LX+1)**2))
      CALL GAUSSIAN_YLMPOL(LX,YLMPOL)
!
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        LNX1=LNX(ISP)
        ALLOCATE(NLPHIDOT(NR,LNX1))
        ALLOCATE(AEPHI(NR,LNX1))
        ALLOCATE(AEPHIDOT(NR,LNX1))
        CALL SETUP$GETR8A('QPHIDOT',NR*LNX1,NLPHIDOT)
        CALL SETUP$GETR8A('AEPHI',NR*LNX1,AEPHI)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX1,AEPHIDOT)
!
        NPOW=GAUSSFITKAUGMENT_NPOW !POLYNOMIAL HAS X^NPOW AS HIGHEST POWER
        NE=GAUSSFITKAUGMENT_NE
        DO IE=1,NE
          SVAR=GAUSSFITKAUGMENT_R1*GAUSSFITKAUGMENT_SCALER**(IE-1)
          E(IE)=1.D0/SVAR**2
        ENDDO
        NIJK=((NPOW+1)*(NPOW+2)*(NPOW+3))/6
        LX=MAXVAL(LOX(:LNX(ISP),ISP))
        NORBX=SUM(2*LOX(:LNX1,ISP)+1)
        NORBS=0
        DO L=0,LX
          IF(SBARLI1(L+1,ISP).GE.1)NORBS=MAX(NORBS,SBARLI1(L+1,ISP)+2*L)
        ENDDO
        POTPAR(ISP)%GAUSSKAUGMENT%NE     =NE
        POTPAR(ISP)%GAUSSJAUGMENT%NE     =NE
        POTPAR(ISP)%GAUSSKAUGMENT%NIJK   =NIJK
        POTPAR(ISP)%GAUSSJAUGMENT%NIJK   =NIJK
        POTPAR(ISP)%GAUSSKAUGMENT%NORB   =NORBX
        POTPAR(ISP)%GAUSSJAUGMENT%NORB   =NORBS
        ALLOCATE(POTPAR(ISP)%GAUSSKAUGMENT%E(NE))
        ALLOCATE(POTPAR(ISP)%GAUSSJAUGMENT%E(NE))
        POTPAR(ISP)%GAUSSKAUGMENT%E(:)   =E(:NE)
        POTPAR(ISP)%GAUSSJAUGMENT%E(:)   =E(:NE)
        ALLOCATE(POTPAR(ISP)%GAUSSKAUGMENT%C(NIJK,NE,NORBX))
        ALLOCATE(POTPAR(ISP)%GAUSSJAUGMENT%C(NIJK,NE,NORBS))
!
!       == DETERMINE WEIGHT FOR FIT ============================================
        ALLOCATE(W(NR))
        W(:)=1.D0
!
        RAD=POTPAR(ISP)%RAD
        ALLOCATE(KPRIME(NR))
        ALLOCATE(JBAR(NR))
        ALLOCATE(CPOWK(NPOW,NE))
        ALLOCATE(CPOWJ(NPOW,NE))
        ALLOCATE(WORKK(NIJK,NE))
        ALLOCATE(WORKJ(NIJK,NE))
        ALLOCATE(CK(NIJK,NE))
        ALLOCATE(CJ(NIJK,NE))
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          QBAR=POTPAR(ISP)%QBAR(LN)
          LN1=POTPAR(ISP)%LNSCATT(LN)
          TSCATT=(LN.EQ.LN1)
!
!         ======================================================================
!         == AUGMENTATION FUNCTION                                            ==
!         ======================================================================
          KPRIME=0.D0
          DO IR=1,NR
            IF(R(IR).GT.RAD) EXIT
            CALL LMTO$SOLIDHANKELRAD(L,R(IR),K2,KVAL,KDER)
            CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,JVAL,JDER)
            SVAR=(JVAL-NLPHIDOT(IR,LN1)*POTPAR(ISP)%JBARTOPHIDOT(LN))/QBAR
            KPRIME(IR)=AEPHI(IR,LN)    *POTPAR(ISP)%KTOPHI(LN) &
       &              +AEPHIDOT(IR,LN1)*POTPAR(ISP)%KTOPHIDOT(LN) &
       &              -SVAR
          ENDDO
!
          NPOW2=INT(0.5D0*REAL(NPOW-L))
          CPOWK(:,:)=0.D0
          CALL GAUSSIAN_FITGAUSS(GID,NR,W,L,KPRIME,NE,NPOW2+1,E(:NE) &
       &                                                    ,CPOWK(:NPOW2+1,:))
          IF(TSCATT) THEN
            JBAR(:)=0.D0
            DO IR=1,NR
              IF(R(IR).GT.RAD) EXIT
              JBAR(IR)=(AEPHIDOT(IR,LN1)-NLPHIDOT(IR,LN1)) &
       &                                          *POTPAR(ISP)%JBARTOPHIDOT(LN1)
            ENDDO
            CPOWJ(:,:)=0.D0
            CALL GAUSSIAN_FITGAUSS(GID,NR,W,L,JBAR,NE,NPOW2+1,E(:NE) &
       &                                                    ,CPOWJ(:NPOW2+1,:))
          END IF
!
!         ======================================================================
!         ==  PLOT FIT TO INSPECT QUALITY                                     ==
!         ======================================================================
          IF(TTEST) THEN
            WRITE(STRING1,*)ISP
            WRITE(STRING2,*)LN
            FILE='GAUSSFITKAUGMENT_A_'//TRIM(ADJUSTL(STRING1))
            FILE=TRIM(FILE)//'_'//TRIM(ADJUSTL(STRING2))//'.DAT'
            CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,FILE)
            CALL FILEHANDLER$UNIT('HOOK',NFIL)
            DO IR=1,NR
              KVAL=0.D0
              JVAL=0.D0
              DO IE=1,NE
                DO I=0,NPOW2
                  KVAL=KVAL+R(IR)**(L+2*I)*EXP(-E(IE)*R(IR)**2)*CPOWK(I+1,IE)
                  JVAL=JVAL+R(IR)**(L+2*I)*EXP(-E(IE)*R(IR)**2)*CPOWJ(I+1,IE)
                ENDDO  
              ENDDO
              IF(TSCATT) THEN
                WRITE(NFIL,FMT='(10F20.5)')R(IR),KPRIME(IR),KVAL,JVAL
              ELSE
                WRITE(NFIL,FMT='(10F20.5)')R(IR),KPRIME(IR),KVAL
              END IF
            ENDDO
            CALL FILEHANDLER$CLOSE('HOOK')
            CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOK')
          END IF
!
!         ======================================================================
!         ==  EXPAND                                                          ==
!         ======================================================================
          WORKK(:,:)=0.D0
          WORKJ(:,:)=0.D0
          NIJKXPOW2=0
          DO N=0,NPOW2
!           == MULTIPLY WITH (X^2+Y^2+Z^2)^N ===================================
            DO I=0,N
              CALL BINOMIALCOEFFICIENT(N,I,SVAR1)
              DO J=0,N-I
                CALL BINOMIALCOEFFICIENT(N-I,J,SVAR2)
                K=N-I-J
                CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND1,2*I,2*J,2*K)
                IF(IND1.GT.NIJK) CYCLE
                WORKK(IND1,:)=WORKK(IND1,:)  +CPOWK(N+1,:)*SVAR1*SVAR2
                IF(TSCATT)WORKJ(IND1,:)=WORKJ(IND1,:)+CPOWJ(N+1,:)*SVAR1*SVAR2
                NIJKXPOW2=MAX(IND1,NIJKXPOW2)
              ENDDO
            ENDDO
          ENDDO
!
!         ======================================================================
!         == MULTIPLY WITH SPHERICAL HARMONICS                                ==
!         ======================================================================
          DO IM=1,2*L+1
            CK(:,:)=0.D0
            CJ(:,:)=0.D0
            LM=L**2+IM
            IND1X=SIZE(YLMPOL(:,LM))
            DO IND1=1,IND1X
              IF(YLMPOL(IND1,LM).EQ.0.D0) CYCLE
              CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND1,I,J,K)
              DO IND=1,NIJKXPOW2
                CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I1,J1,K1)
                I1=I1+I
                J1=J1+J
                K1=K1+K
                CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND2,I1,J1,K1)
                IF(IND2.GT.NIJK) CYCLE
                CK(IND2,:)   =CK(IND2,:)   +WORKK(IND,:)*YLMPOL(IND1,LM)
                IF(TSCATT)CJ(IND2,:)=CJ(IND2,:)+WORKJ(IND,:)*YLMPOL(IND1,LM)
              ENDDO
            ENDDO
            IORB=SUM(2*LOX(:LN-1,ISP)+1)+IM
            IORBS=SBARLI1(L+1,ISP)-1+IM
            POTPAR(ISP)%GAUSSKAUGMENT%C(:,:,IORB)=CK(:,:)
            IF(TSCATT)POTPAR(ISP)%GAUSSJAUGMENT%C(:,:,IORBS)=CJ(:,:)
          ENDDO ! END OF LOOP OVER ORBITALS (IM)
        ENDDO  ! END OF LOOP OVER L
!
!       ========================================================================
!       ==  PLOT FIT TO INSPECT QUALITY                                       ==
!       ========================================================================
        IF(TTEST) THEN
          DO IORB=1,NORBX
            WRITE(STRING1,*)ISP
            WRITE(STRING2,*)IORB
            FILE='GAUSSFITKAUGMENT_B_'//TRIM(ADJUSTL(STRING1))
            FILE=TRIM(FILE)//'_'//TRIM(ADJUSTL(STRING2))//'.DAT'
            CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,FILE)
            CALL FILEHANDLER$UNIT('HOOK',NFIL)
            CK(:,:)=POTPAR(ISP)%GAUSSKAUGMENT%C(:,:,IORB)
            DO IR=1,NR
              KVAL=0.D0
              DO IE=1,NE
                R0=(/5.D0,0.D0,0.D0/)*R(IR)
                CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,E(IE),CK(:,IE),R0,SVAR)
                KVAL=KVAL+SVAR
              ENDDO
              WRITE(NFIL,FMT='(10F20.5)')5.D0*R(IR),KVAL,JVAL
            ENDDO
            CALL FILEHANDLER$CLOSE('HOOK')
            CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOK')
          ENDDO
          DO IORB=1,NORBS
            WRITE(STRING1,*)ISP
            WRITE(STRING2,*)IORB
            FILE='GAUSSFITJAUGMENT_B_'//TRIM(ADJUSTL(STRING1))
            FILE=TRIM(FILE)//'_'//TRIM(ADJUSTL(STRING2))//'.DAT'
            CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,FILE)
            CALL FILEHANDLER$UNIT('HOOK',NFIL)
            CJ(:,:)=POTPAR(ISP)%GAUSSJAUGMENT%C(:,:,IORB)
            DO IR=1,NR
              JVAL=0.D0
              DO IE=1,NE
                R0=(/5.D0,0.D0,0.D0/)*R(IR)
                CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,E(IE),CJ(:,IE),R0,SVAR)
                JVAL=JVAL+SVAR
              ENDDO
              WRITE(NFIL,FMT='(10F20.5)')5.D0*R(IR),JVAL
            ENDDO
            CALL FILEHANDLER$CLOSE('HOOK')
            CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOK')
          ENDDO
        END IF
!
!       ========================================================================
!       ==  FINISH OFF                                                        ==
!       ========================================================================
        DEALLOCATE(R)
        DEALLOCATE(W)
        DEALLOCATE(CPOWK)
        DEALLOCATE(CPOWJ)
        DEALLOCATE(WORKK)
        DEALLOCATE(WORKJ)
        DEALLOCATE(CK)
        DEALLOCATE(CJ)
        DEALLOCATE(KPRIME)
        DEALLOCATE(JBAR)
        DEALLOCATE(AEPHI)
        DEALLOCATE(AEPHIDOT)
        DEALLOCATE(NLPHIDOT)
      ENDDO   ! END LOOP OVER SPECIES
                                         CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GAUSSFITKJTAILS()
!     **************************************************************************
!     ** CONSTRUCTS THE UNSCREENED HANKEL FUNCTIONS, WITH THE SIMGULARITY     **
!     ** REPLACED BY A NODELESS SCATTERING FUNCTION AND DETERMINES ITS        **
!     ** COEFFICIENTS IN A GAUSSIAN EXPANSION.                                **
!     **                                                                      **
!     ** SIMILARLY THE UNCREENED HANKEL AND BESSEL FUNCTIONS WITH             **
!     ** AN EXPONENTIAL TAIL R^L*(A+BR^2)*EXP(-LAMBDA*R) ARE EXPANDED         **
!     ** INTO GAUSSIANS.                                                      **
!     **                                                                      **
!     ** CAUTION!!! THE UNSCREENED HANKEL FUNCTION WITH LONG-RANGE TAIL       **
!     **   DOES NOT WORK BECAUSE IT REQUIRES GAUSSIANS WITH A LONGER RANGE    **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      USE LMTO_MODULE, ONLY : K2,POTPAR,NSP,LNX,LOX,SBARLI1 &
     &                       ,GAUSSFITKBARPRIME_NPOW,GAUSSFITKBARPRIME_NE &
     &                       ,GAUSSFITKBARPRIME_R1,GAUSSFITKBARPRIME_SCALER
      INTEGER(4),SAVE        :: GID=-1
      LOGICAL(4)             :: TTEST=.false.
      INTEGER(4),PARAMETER   :: NEX=30
      REAL(8)   ,PARAMETER   :: LAMBDA=1.D0
      REAL(8)                :: E(10)
      INTEGER(4)             :: NPOW
      INTEGER(4)             :: NPOW2
      INTEGER(4)             :: NE
      INTEGER(4)             :: NIJK
      INTEGER(4)             :: GID1
      INTEGER(4)             :: NR,NR1
      INTEGER(4)             :: L
      INTEGER(4)             :: LX
      INTEGER(4)             :: LM
      INTEGER(4)             :: NORBX
      REAL(8)                :: R1,DEX
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: W(:)
      REAL(8)   ,ALLOCATABLE :: KPRIME(:)
      REAL(8)   ,ALLOCATABLE :: JBARPRIME(:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: CPOWK(:,:)    !(NPOW2,NE)
      REAL(8)   ,ALLOCATABLE :: CPOWJBAR(:,:) !(NPOW2,NE)
      REAL(8)   ,ALLOCATABLE :: CK(:,:)       !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: CJBAR(:,:)    !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: WORKK(:,:)    !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: WORKJ(:,:)    !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: YLMPOL(:,:)
      REAL(8)                :: RAD
      REAL(8)                :: QBAR
      REAL(8)                :: SVAR,SVAR1,SVAR2,JVAL,KVAL,JDER
      INTEGER(4)             :: ISP,IE,I,J,K,I1,J1,K1,LN,IM,IR,N,IND,IND1,IND2
      INTEGER(4)             :: IND1X,IORB
      INTEGER(4)             :: NIJKXPOW2
      CHARACTER(128)         :: FILE
      CHARACTER(8)           :: STRING1,STRING2
      REAL(8)                :: R0(3)
      INTEGER(4)             :: NFIL
!     **************************************************************************
                                  CALL TRACE$PUSH('LMTO_GAUSSFITKJTAILS')
!
!     ==========================================================================
!     == POLYNOMIAL COEFFICIENTS OF SPHERICAL HARMONICS TIMES R**L            ==
!     ==========================================================================
      LX=MAXVAL(LOX)
      ALLOCATE(YLMPOL((LX+1)*(LX+2)*(LX+3)/6,(LX+1)**2))
      CALL GAUSSIAN_YLMPOL(LX,YLMPOL)
!
!     ==========================================================================
!     == DEFINE A RADIAL GRID ON WHICH THE FIT IS PERFORMED                   ==
!     ==========================================================================
      IF(GID.EQ.-1) THEN
        CALL RADIAL$NEW('SHLOG',GID)
        CALL RADIAL$GRIDPARAMETERS(0.01D0,0.5D0,25.D0,R1,DEX,NR)
        CALL RADIAL$SETR8(GID,'R1',R1)
        CALL RADIAL$SETR8(GID,'DEX',DEX)
        CALL RADIAL$SETI4(GID,'NR',NR)
      ELSE
        CALL RADIAL$GETI4(GID,'NR',NR)
      END IF
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      ALLOCATE(W(NR))
!
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
      DO ISP=1,NSP
        NPOW=GAUSSFITKBARPRIME_NPOW !POLYNOMIAL HAS X^NPOW AS HIGHEST POWER
        NE=GAUSSFITKBARPRIME_NE
        DO IE=1,NE
          SVAR=GAUSSFITKBARPRIME_R1*GAUSSFITKBARPRIME_SCALER**(IE-1)
          E(IE)=1.D0/SVAR**2
        ENDDO
        NIJK=((NPOW+1)*(NPOW+2)*(NPOW+3))/6
        LX=MAXVAL(LOX(:LNX(ISP),ISP))
        NORBX=0
        DO L=0,LX
          NORBX=MAX(NORBX,SBARLI1(L+1,ISP)+2*L)
        ENDDO
        POTPAR(ISP)%GAUSSTAILEDK%NE     =NE
        POTPAR(ISP)%GAUSSTAILEDJBAR%NE  =NE
        POTPAR(ISP)%GAUSSTAILEDK%NIJK   =NIJK
        POTPAR(ISP)%GAUSSTAILEDJBAR%NIJK=NIJK
        POTPAR(ISP)%GAUSSTAILEDK%NORB   =NORBX
        POTPAR(ISP)%GAUSSTAILEDJBAR%NORB=NORBX
        ALLOCATE(POTPAR(ISP)%GAUSSTAILEDK%E(NE))
        ALLOCATE(POTPAR(ISP)%GAUSSTAILEDJBAR%E(NE))
        POTPAR(ISP)%GAUSSTAILEDK%E(:)=E(:NE)
        POTPAR(ISP)%GAUSSTAILEDJBAR%E(:)=E(:NE)
        ALLOCATE(POTPAR(ISP)%GAUSSTAILEDK%C(NIJK,NE,NORBX))
        ALLOCATE(POTPAR(ISP)%GAUSSTAILEDJBAR%C(NIJK,NE,NORBX))

!
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID1)
        CALL SETUP$GETI4('NR',NR1)
        ALLOCATE(NLPHIDOT(NR1,LNX(ISP)))
        CALL SETUP$GETR8A('QPHIDOT',NR1*LNX(ISP),NLPHIDOT)
!
!       == DETERMINE WEIGHT FOR FIT ============================================
        W(:)=1.D0
!
        RAD=POTPAR(ISP)%RAD
        ALLOCATE(KPRIME(NR))
        ALLOCATE(JBARPRIME(NR))
        ALLOCATE(CPOWK(NPOW,NE))
        ALLOCATE(CPOWJBAR(NPOW,NE))
        ALLOCATE(WORKJ(NIJK,NE))
        ALLOCATE(WORKK(NIJK,NE))
        ALLOCATE(CK(NIJK,NE))
        ALLOCATE(CJBAR(NIJK,NE))
        DO LN=1,LNX(ISP)
          IF(LN.NE.POTPAR(ISP)%LNSCATT(LN)) CYCLE
          L=LOX(LN,ISP)
          QBAR=POTPAR(ISP)%QBAR(LN)
!
!         ======================================================================
!         == KPRIME AND JPRIME WITH TAILS                                     ==
!         ======================================================================
          CALL LMTO_KJBARTAILS(GID,NR,L,RAD,K2,QBAR,LAMBDA,KPRIME,JBARPRIME)
          DO IR=1,NR
            IF(R(IR).LE.RAD) THEN
              CALL RADIAL$VALUE(GID1,NR1,NLPHIDOT(:,LN),R(IR),SVAR)
              JBARPRIME(IR)=SVAR*POTPAR(ISP)%JBARTOPHIDOT(LN)
              CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,JVAL,JDER)
              KPRIME(IR)=(JVAL-JBARPRIME(IR))/QBAR
            END IF
          ENDDO
!
          NPOW2=INT(0.5D0*REAL(NPOW-L))
          CPOWK(:,:)=0.D0
          CPOWJBAR(:,:)=0.D0
          CALL GAUSSIAN_FITGAUSS(GID,NR,W,L,KPRIME,NE,NPOW2+1,E(:NE) &
       &                                                    ,CPOWK(:NPOW2+1,:))
          CALL GAUSSIAN_FITGAUSS(GID,NR,W,L,JBARPRIME,NE,NPOW2+1,E(:NE) &
       &                                                 ,CPOWJBAR(:NPOW2+1,:))
!
!         ======================================================================
!         ==  PLOT FIT TO INSPECT QUALITY                                     ==
!         ======================================================================
          IF(TTEST) THEN
            WRITE(STRING1,*)ISP
            WRITE(STRING2,*)L
            FILE='GAUSSFITKJTAILS_A_'//TRIM(ADJUSTL(STRING1))
            FILE=TRIM(FILE)//'_'//TRIM(ADJUSTL(STRING2))//'.DAT'
            CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,FILE)
            CALL FILEHANDLER$UNIT('HOOK',NFIL)
            DO IR=1,NR
              KVAL=0.D0
              JVAL=0.D0
              DO IE=1,NE
                DO I=0,NPOW2
                  KVAL=KVAL+R(IR)**(L+2*I)*EXP(-E(IE)*R(IR)**2)*CPOWK(I+1,IE)
                  JVAL=JVAL+R(IR)**(L+2*I)*EXP(-E(IE)*R(IR)**2)*CPOWJBAR(I+1,IE)
                ENDDO  
              ENDDO
              WRITE(NFIL,FMT='(10F20.5)')R(IR),KPRIME(IR),JBARPRIME(IR) &
        &                                 ,KVAL,JVAL
            ENDDO
            CALL FILEHANDLER$CLOSE('HOOK')
            CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOK')
          END IF
!
!         ======================================================================
!         ==  EXPAND                                                          ==
!         ======================================================================
          WORKK(:,:)=0.D0
          WORKJ(:,:)=0.D0
          NIJKXPOW2=0
          DO N=0,NPOW2
!           == MULTIPLY WITH (X^2+Y^2+Z^2)^N ===================================
            DO I=0,N
              CALL BINOMIALCOEFFICIENT(N,I,SVAR1)
              DO J=0,N-I
                CALL BINOMIALCOEFFICIENT(N-I,J,SVAR2)
                K=N-I-J
                CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND1,2*I,2*J,2*K)
                IF(IND1.GT.NIJK) CYCLE
                WORKK(IND1,:)=WORKK(IND1,:)  +CPOWK(N+1,:)   *SVAR1*SVAR2
                WORKJ(IND1,:)=WORKJ(IND1,:)  +CPOWJBAR(N+1,:)*SVAR1*SVAR2
                NIJKXPOW2=MAX(IND1,NIJKXPOW2)
              ENDDO
            ENDDO
          ENDDO
!
!         ======================================================================
!         == MULTIPLY WITH SPHERICAL HARMONICS                                ==
!         ======================================================================
          DO IM=1,2*L+1
            CK(:,:)=0.D0
            CJBAR(:,:)=0.D0
            LM=L**2+IM
            IND1X=SIZE(YLMPOL(:,LM))
            DO IND1=1,IND1X
              IF(YLMPOL(IND1,LM).EQ.0.D0) CYCLE
              CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND1,I,J,K)
              DO IND=1,NIJKXPOW2
                CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I1,J1,K1)
                I1=I1+I
                J1=J1+J
                K1=K1+K
                CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND2,I1,J1,K1)
                IF(IND2.GT.NIJK) CYCLE
                CK(IND2,:)   =CK(IND2,:)   +WORKK(IND,:)*YLMPOL(IND1,LM)
                CJBAR(IND2,:)=CJBAR(IND2,:)+WORKJ(IND,:)*YLMPOL(IND1,LM)
              ENDDO
            ENDDO
            IORB=SBARLI1(L+1,ISP)-1+IM
            POTPAR(ISP)%GAUSSTAILEDK%C(:,:,IORB)=CK(:,:)
            POTPAR(ISP)%GAUSSTAILEDJBAR%C(:,:,IORB)=CJBAR(:,:)
          ENDDO ! END OF LOOP OVER ORBITALS (IM)
        ENDDO  ! END OF LOOP OVER L
!
!       ========================================================================
!       ==  PLOT FIT TO INSPECT QUALITY                                       ==
!       ========================================================================
        IF(TTEST) THEN
          DO IORB=1,NORBX
            WRITE(STRING1,*)ISP
            WRITE(STRING2,*)IORB
            FILE='GAUSSFITKJTAILS_B_'//TRIM(ADJUSTL(STRING1))
            FILE=TRIM(FILE)//'_'//TRIM(ADJUSTL(STRING2))//'.DAT'
            CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,FILE)
            CALL FILEHANDLER$UNIT('HOOK',NFIL)
            CK(:,:)=POTPAR(ISP)%GAUSSTAILEDK%C(:,:,IORB)
            CJBAR(:,:)=POTPAR(ISP)%GAUSSTAILEDJBAR%C(:,:,IORB)
            DO IR=1,NR
              KVAL=0.D0
              JVAL=0.D0
              DO IE=1,NE
                R0=(/5.D0,0.D0,0.D0/)*R(IR)
                CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,E(IE),CK(:,IE),R0,SVAR)
                KVAL=KVAL+SVAR
                CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,E(IE),CJBAR(:,IE),R0,SVAR)
                JVAL=JVAL+SVAR
              ENDDO
              WRITE(NFIL,FMT='(10F20.5)')5.D0*R(IR),KVAL,JVAL
            ENDDO
            CLOSE(10001)
            CALL FILEHANDLER$CLOSE('HOOK')
            CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOK')
          ENDDO
        END IF
!
!       ========================================================================
!       ==  FINISH OFF                                                        ==
!       ========================================================================
        DEALLOCATE(CPOWK)
        DEALLOCATE(CPOWJBAR)
        DEALLOCATE(WORKJ)
        DEALLOCATE(WORKK)
        DEALLOCATE(CK)
        DEALLOCATE(CJBAR)
        DEALLOCATE(KPRIME)
        DEALLOCATE(JBARPRIME)
        DEALLOCATE(NLPHIDOT)
      ENDDO   ! END LOOP OVER SPECIES
      DEALLOCATE(R)
      DEALLOCATE(W)
                                         CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GAUSSFITKPRIME()
!     **************************************************************************
!     ** CONSTRUCTS THE UNSCREENED HANKEL FUNCTIONS, WITH THE SINGULARITY     **
!     ** REPLACED BY A NODELESS SCATTERING FUNCTION AND DETERMINS ITS         **
!     ** COEFFICIENTS IN A GAUSSIAN EXPANSION.                                **
!     **                                                                      **
!     ** CAUTION!!! THE UNSCREENED HANKEL FUNCTION WITH LONG-RANGE TAIL       **
!     **   DOES NOT WORK BECAUSE IT REQUIRES GAUSSIANS WITH A LONGER RANGE    **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      USE LMTO_MODULE, ONLY : K2,POTPAR,NSP,LNX,LOX,SBARLI1 &
     &                       ,GAUSSFITKPRIME_NE,GAUSSFITKPRIME_R1 &
     &                       ,GAUSSFITKPRIME_SCALER,GAUSSFITKPRIME_NPOW
      INTEGER(4),SAVE        :: GID=-1
      INTEGER(4),PARAMETER   :: NEX=30
      REAL(8)   ,PARAMETER   :: LAMBDA=1.D0
      LOGICAL(4),PARAMETER   :: TTEST=.false.
      REAL(8)                :: E(NEX)
      INTEGER(4)             :: NPOW
      INTEGER(4)             :: NPOW2
      INTEGER(4)             :: NE
      INTEGER(4)             :: NIJK
      INTEGER(4)             :: GID1
      INTEGER(4)             :: NR,NR1
      INTEGER(4)             :: L
      INTEGER(4)             :: LX
      INTEGER(4)             :: LM
      INTEGER(4)             :: NORBX
      REAL(8)                :: R1,DEX
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: W(:)
      REAL(8)   ,ALLOCATABLE :: KPRIMEL(:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: CPOWKL(:,:)   !(NPOW2,NE)
      REAL(8)   ,ALLOCATABLE :: CKL(:,:)      !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: WORKKL(:,:)   !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: YLMPOL(:,:)
      REAL(8)                :: RAD
      REAL(8)                :: QBAR
      REAL(8)                :: SVAR,SVAR1,SVAR2,JVAL,KVAL,JDER
      INTEGER(4)             :: ISP,IE,I,J,K,I1,J1,K1,LN,IM,IR,N,IND,IND1,IND2
      INTEGER(4)             :: IND1X,IORB
      INTEGER(4)             :: NIJKXPOW2
      CHARACTER(128)         :: FILE
      CHARACTER(8)           :: STRING1,STRING2
      REAL(8)                :: R0(3),RI
      INTEGER(4)             :: NFIL
!     **************************************************************************
                                         CALL TRACE$PUSH('LMTO_GAUSSFITKPRIME')
!
!     ==========================================================================
!     == POLYNOMIAL COEFFICIENTS OF SPHERICAL HARMONICS TIMES R**L            ==
!     ==========================================================================
      LX=MAXVAL(LOX)
      ALLOCATE(YLMPOL((LX+1)*(LX+2)*(LX+3)/6,(LX+1)**2))
      CALL GAUSSIAN_YLMPOL(LX,YLMPOL)
!
!     ==========================================================================
!     == DEFINE A RADIAL GRID ON WHICH THE FIT IS PERFORMED                   ==
!     ==========================================================================
      IF(GID.EQ.-1) THEN
        CALL RADIAL$NEW('SHLOG',GID)
!       == THE GRID MUST REACH FAR OUT TO AVOID SHARP INCREASE =================
!       == OUTSIDE THE REGION CONSIDERED FOR THE FIT ===========================
        CALL RADIAL$GRIDPARAMETERS(0.01D0,5.D0,100.D0,R1,DEX,NR)
        CALL RADIAL$SETR8(GID,'R1',R1)
        CALL RADIAL$SETR8(GID,'DEX',DEX)
        CALL RADIAL$SETI4(GID,'NR',NR)
      END IF
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      ALLOCATE(W(NR))
!
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
      DO ISP=1,NSP
        LX=MAXVAL(LOX(:LNX(ISP),ISP))
        NPOW=MAX(LX,GAUSSFITKPRIME_NPOW)
!       == THE PARAMETERS NE=7 AND E=1/R^2 WITH R0.25*2^(IE-1) IS GOOD =========
!       == UP TO A RADIUS OF ABOUT 10 A.U. =====================================
        NE=GAUSSFITKPRIME_NE
        DO IE=1,NE
          SVAR=GAUSSFITKPRIME_R1*GAUSSFITKPRIME_SCALER**(IE-1)
          E(IE)=1.D0/SVAR**2   
        ENDDO
        NIJK=((NPOW+1)*(NPOW+2)*(NPOW+3))/6
        NORBX=0
        DO L=0,LX
          NORBX=MAX(NORBX,SBARLI1(L+1,ISP)+2*L)
        ENDDO
        POTPAR(ISP)%GAUSSKPRIME%NE=NE
        POTPAR(ISP)%GAUSSKPRIME%NIJK=NIJK
        POTPAR(ISP)%GAUSSKPRIME%NORB=NORBX
        ALLOCATE(POTPAR(ISP)%GAUSSKPRIME%E(NE))
        POTPAR(ISP)%GAUSSKPRIME%E(:)=E(:NE)
        ALLOCATE(POTPAR(ISP)%GAUSSKPRIME%C(NIJK,NE,NORBX))
!
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID1)
        CALL SETUP$GETI4('NR',NR1)
        ALLOCATE(NLPHIDOT(NR1,LNX(ISP)))
        CALL SETUP$GETR8A('QPHIDOT',NR1*LNX(ISP),NLPHIDOT)
        W(:)=1.D0
!
        RAD=POTPAR(ISP)%RAD
        ALLOCATE(KPRIMEL(NR))
        ALLOCATE(CPOWKL(NPOW+1,NE))
        ALLOCATE(WORKKL(NIJK,NE))
        ALLOCATE(CKL(NIJK,NE))
        DO LN=1,LNX(ISP)
          IF(LN.NE.POTPAR(ISP)%LNSCATT(LN)) CYCLE
          L=LOX(LN,ISP)
          QBAR=POTPAR(ISP)%QBAR(LN)
!
!         ======================================================================
!         == PUT KPRIME ON RADIAL GRID AND FIT TO GAUSSIANS                   ==
!         ======================================================================
          DO IR=1,NR
            IF(R(IR).LE.RAD) THEN
              CALL RADIAL$VALUE(GID1,NR1,NLPHIDOT(:,LN),R(IR),SVAR)
              SVAR=SVAR*POTPAR(ISP)%JBARTOPHIDOT(LN)
              CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,JVAL,JDER)
              KPRIMEL(IR)=(JVAL-SVAR)/QBAR
            ELSE
              CALL LMTO$SOLIDHANKELRAD(L,R(IR),K2,KPRIMEL(IR),SVAR)
            END IF
          ENDDO
!
          NPOW2=INT(0.5D0*REAL(NPOW-L))
          CPOWKL(:,:)=0.D0
          CALL GAUSSIAN_FITGAUSS(GID,NR,W,L,KPRIMEL,NE,NPOW2+1,E(:NE) &
       &                                                   ,CPOWKL(:NPOW2+1,:))
!
!         ======================================================================
!         ==  PLOT FIT TO INSPECT QUALITY                                     ==
!         ======================================================================
          IF(TTEST) THEN
            WRITE(STRING1,*)ISP
            WRITE(STRING2,*)L
            FILE='GAUSSFITKPRIME_A_'//TRIM(ADJUSTL(STRING1))
            FILE=TRIM(FILE)//'_'//TRIM(ADJUSTL(STRING2))//'.DAT'
            CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,FILE)
            CALL FILEHANDLER$UNIT('HOOK',NFIL)
            DO IR=1,NR
              RI=R(IR)*2.D0  ! GRID IS EXTENDED ON PURPOSE
              SVAR=0.D0
              DO IE=1,NE
                DO I=0,NPOW2
                  SVAR=SVAR+RI**(L+2*I)*EXP(-E(IE)*RI**2)*CPOWKL(I+1,IE)
                ENDDO  
              ENDDO
              IF(RI.LT.RAD) THEN
                CALL RADIAL$VALUE(GID,NR,KPRIMEL,RI,SVAR1)
              ELSE
                CALL LMTO$SOLIDHANKELRAD(L,RI,K2,SVAR1,SVAR2)
              END IF
              WRITE(NFIL,FMT='(10F20.5)')RI,SVAR1,SVAR
            ENDDO
            CALL FILEHANDLER$CLOSE('HOOK')
            CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOK')
          END IF
!
!         ======================================================================
!         ==  EXPAND                                                          ==
!         ======================================================================
          WORKKL(:,:)=0.D0
          NIJKXPOW2=0
          DO N=0,NPOW2
!           == MULTIPLY WITH (X^2+Y^2+Z^2)^N ===================================
            DO I=0,N
              CALL BINOMIALCOEFFICIENT(N,I,SVAR1)
              DO J=0,N-I
                CALL BINOMIALCOEFFICIENT(N-I,J,SVAR2)
                K=N-I-J
                CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND1,2*I,2*J,2*K)
                IF(IND1.GT.NIJK) CYCLE
                WORKKL(IND1,:)=WORKKL(IND1,:)+CPOWKL(N+1,:) *SVAR1*SVAR2
                NIJKXPOW2=MAX(IND1,NIJKXPOW2)
              ENDDO
            ENDDO
          ENDDO
!
!         ======================================================================
!         == MULTIPLY WITH SPHERICAL HARMONICS                                ==
!         ======================================================================
          DO IM=1,2*L+1
            CKL(:,:)=0.D0
            LM=L**2+IM
            IND1X=SIZE(YLMPOL(:,LM))
            DO IND1=1,IND1X
              IF(YLMPOL(IND1,LM).EQ.0.D0) CYCLE
              CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND1,I,J,K)
              DO IND=1,NIJKXPOW2
                CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I1,J1,K1)
                I1=I1+I
                J1=J1+J
                K1=K1+K
                CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND2,I1,J1,K1)
                IF(IND2.GT.NIJK) CYCLE
                CKL(IND2,:)  =CKL(IND2,:)  +WORKKL(IND,:)*YLMPOL(IND1,LM)
              ENDDO
            ENDDO
            IORB=SBARLI1(L+1,ISP)-1+IM
            POTPAR(ISP)%GAUSSKPRIME%C(:,:,IORB)=CKL(:,:)
          ENDDO ! END OF LOOP OVER ORBITALS (IM)
        ENDDO  ! END OF LOOP OVER L
!
!       ========================================================================
!       ==  PLOT FIT TO INSPECT QUALITY                                       ==
!       ========================================================================
        IF(TTEST) THEN
          DO LM=1,NORBX
            WRITE(STRING1,*)ISP
            WRITE(STRING2,*)LM
            FILE='GAUSSFITKPRIME_B_'//TRIM(ADJUSTL(STRING1))
            FILE=TRIM(FILE)//'_'//TRIM(ADJUSTL(STRING2))//'.DAT'
            CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,FILE)
            CALL FILEHANDLER$UNIT('HOOK',NFIL)
            CKL=POTPAR(ISP)%GAUSSKPRIME%C(:,:,LM)
            DO IR=1,NR
              RI=2.D0*R(IR)
              KVAL=0.D0
              DO IE=1,NE
                R0=(/1.D0,0.D0,0.D0/)*RI
                CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,E(IE),CKL(:,IE),R0,SVAR)
                KVAL=KVAL+SVAR
              ENDDO
              WRITE(NFIL,FMT='(10F20.5)')RI,KVAL
            ENDDO
            CALL FILEHANDLER$CLOSE('HOOK')
            CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOK')
          ENDDO
        END IF
!
!       ========================================================================
!       ==  FINISH OFF                                                        ==
!       ========================================================================
        DEALLOCATE(CPOWKL)
        DEALLOCATE(WORKKL)
        DEALLOCATE(CKL)
        DEALLOCATE(KPRIMEL)
        DEALLOCATE(NLPHIDOT)
      ENDDO   ! END LOOP OVER SPECIES
      DEALLOCATE(R)
      DEALLOCATE(W)
                                         CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBOFROMTAILEDKJ()
!     **************************************************************************
!     **  ONSITE EXPANSION OF TIGHT-BINDING ORBITALS FROM                     **
!     **  TAILED HANKEL AND BESSEL FUNCTIONS                                  **
!     **                                                                      **
!     **  RESULT IS PLACED IN GAUSSORB_T                                      **
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      USE LMTO_MODULE, ONLY : POTPAR,ISPECIES,SBAR,GAUSSORB_T
      LOGICAL           :: TNEW
      INTEGER(4)        :: NAT
      INTEGER(4)        :: NNS
      INTEGER(4)        :: IAT
      INTEGER(4)        :: ISP
      INTEGER(4)        :: NORB
      INTEGER(4)        :: NE
      INTEGER(4)        :: NIJK
      INTEGER(4)        :: NN,IORB1,IORB2
!     **************************************************************************
      TNEW=.NOT.ALLOCATED(GAUSSORB_T)
      IF(TNEW) THEN
        NAT=SIZE(ISPECIES)
        ALLOCATE(GAUSSORB_T(NAT))
      END IF
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.SBAR(NN)%IAT2) CYCLE
        IF(SBAR(NN)%IT(1).NE.0) CYCLE
        IF(SBAR(NN)%IT(2).NE.0) CYCLE
        IF(SBAR(NN)%IT(3).NE.0) CYCLE
        IAT=SBAR(NN)%IAT1
        ISP=ISPECIES(IAT)
        NORB=SBAR(NN)%N1
        IF(TNEW) THEN
          NE  =POTPAR(ISP)%GAUSSTAILEDK%NE
          NIJK=POTPAR(ISP)%GAUSSTAILEDK%NIJK
          GAUSSORB_T(IAT)%NIJK=NIJK
          GAUSSORB_T(IAT)%NE=NE
          GAUSSORB_T(IAT)%NORB=NORB
          ALLOCATE(GAUSSORB_T(IAT)%E(NE))
          GAUSSORB_T(IAT)%E(:)=POTPAR(ISP)%GAUSSTAILEDK%E(:)
          ALLOCATE(GAUSSORB_T(IAT)%C(NIJK,NE,NORB))
        END IF
        GAUSSORB_T(IAT)%C(:,:,:)=POTPAR(ISP)%GAUSSTAILEDK%C(:,:,:)
        DO IORB1=1,NORB
          DO IORB2=1,NORB
            GAUSSORB_T(IAT)%C(:,:,IORB1)=GAUSSORB_T(IAT)%C(:,:,IORB1) &
     &                               -POTPAR(ISP)%GAUSSTAILEDJBAR%C(:,:,IORB2) &
     &                               *SBAR(NN)%MAT(IORB2,IORB1)
          ENDDO
        ENDDO
      ENDDO        
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBOFROMKPRIME()
!     **************************************************************************
!     ** CONSTRUCT TIGHT BINDING ORBITALS AS SUPERPOSITION OF HANKEL FUNCTIONS**
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,GAUSSORB,NSP,LNX,LOX,POTPAR,SBAR
      IMPLICIT NONE
      LOGICAL                :: TNEW
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NE,NEA,NEB
      INTEGER(4)             :: NIJK,NIJKA,NIJKB
      INTEGER(4)             :: NORBX,NORB,NORB1,NORB2
      INTEGER(4)             :: ISP,ISP2
      INTEGER(4)             :: IAT,IAT2
      INTEGER(4)             :: IT(3)
      INTEGER(4)             :: LN,IORB,NN,I,J,IM
      INTEGER(4)             :: NNS
      INTEGER(4)             :: NPOW
      LOGICAL(4)             :: TONSITE
      REAL(8)   ,ALLOCATABLE :: QBAR(:,:)
      REAL(8)   ,ALLOCATABLE :: EA(:),EB(:)
      REAL(8)   ,ALLOCATABLE :: GCB(:,:,:)
      REAL(8)   ,ALLOCATABLE :: MAT(:,:)
      REAL(8)   ,ALLOCATABLE :: MAT1(:,:)
      REAL(8)   ,ALLOCATABLE :: MAT2(:,:)
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      REAL(8)                :: RBAS(3,3)
      REAL(8)                :: RA(3),RB(3)
      REAL(8)   ,ALLOCATABLE :: SCALE1(:)
INTEGER(4) :: IE,IND,K,IX
REAL(8)    :: SVAR,ARR(20),X
REAL(8),ALLOCATABLE ::COO(:,:,:)
INTEGER(4),PARAMETER  :: NX=1000
REAL(8)   ,PARAMETER  :: XMAX=10.D0
REAL(8)               :: RGRID(3,NX)
REAL(8)               :: FGRID(NX,20,2)
!     **************************************************************************
                                         CALL TRACE$PUSH('LMTO_NTBOFROMKPRIME')
DO I=1,NX
  RGRID(:,I)=(/1.D0,1D0,1.D0/)/SQRT(3.D0)*XMAX*2.D0*(REAL(I-1)/REAL(NX-1)-0.5D0)
ENDDO
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      CALL CELL$GETR8A('T0',9,RBAS)
      NAT=SIZE(ISPECIES)
!
!     ==========================================================================
!     == ALLOCATE GAUSSORB, WHICH HOLD THE LOCAL ORBITALS IN GAUSS REPRESENTAT.
!     ==========================================================================
      TNEW=.NOT.ALLOCATED(GAUSSORB)
      IF(TNEW) THEN
        ALLOCATE(GAUSSORB(NAT))
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          NE=4
          NPOW=4 !POLYNOMIAL HAS X^NPOW AS HIGHEST POWER
          NIJK=((NPOW+1)*(NPOW+2)*(NPOW+3))/6
          GAUSSORB(IAT)%NIJK=NIJK
          GAUSSORB(IAT)%NE=NE
          ALLOCATE(GAUSSORB(IAT)%E(NE))
          DO I=1,NE
            GAUSSORB(IAT)%E(I)=1.5D0**(-2*(I-1))
          ENDDO
          NORB=0
          DO LN=1,LNX(ISP)
            IF(POTPAR(ISP)%LNSCATT(LN).NE.LN) CYCLE
            NORB=NORB+2*LOX(LN,ISP)+1
          ENDDO  
          GAUSSORB(IAT)%NORB=NORB
          ALLOCATE(GAUSSORB(IAT)%C(NIJK,NE,NORB))
        ENDDO
      END IF
!
      NORBX=MAXVAL(GAUSSORB(:)%NORB)
      ALLOCATE(QBAR(NORBX,NSP))
      QBAR(:,:)=0.D0
      DO ISP=1,NSP
        IORB=0
        DO LN=1,LNX(ISP)
          IF(POTPAR(ISP)%LNSCATT(LN).NE.LN) CYCLE
          DO IM=1,2*LOX(LN,ISP)+1
            IORB=IORB+1
            QBAR(IORB,ISP)=POTPAR(ISP)%QBAR(LN)
          ENDDO
        ENDDO  
      ENDDO
!
      NEB=MAXVAL(POTPAR(:)%GAUSSKPRIME%NE)
      ALLOCATE(EB(NEB))
!
!     ==========================================================================
!     == ALLOCATE GAUSSORB, WHICH HOLD THE LOCAL ORBITALS IN GAUSS REPRESENTAT.
!     ==========================================================================
      NNS=SIZE(SBAR)
      DO IAT=1,NAT
FGRID(:,:,:)=0.D0
        ISP=ISPECIES(IAT)
        NIJKA=GAUSSORB(IAT)%NIJK
        NEA=GAUSSORB(IAT)%NE
        ALLOCATE(EA(NEA))
        EA(:NEA)=GAUSSORB(IAT)%E
        NORB1=GAUSSORB(IAT)%NORB
        ALLOCATE(GCB(NIJKA,NEA,NORBx))
        ALLOCATE(MAT(NIJKA*NEA,NIJKA*NEA))
        ALLOCATE(MAT1(NIJKA*NEA,NORB1))
        ALLOCATE(MAT2(NIJKA*NEA,NORB1))
        RA(:)=R0(:,IAT)
!
!       == DETERMINE GAUSSORB= <G|KPRIME>(1+QBAR*SBAR) =========================
        GAUSSORB(IAT)%C(:,:,:)=0.D0
        DO NN=1,NNS
          IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
          IAT2=SBAR(NN)%IAT2  
          IT(:)=SBAR(NN)%IT(:)
          TONSITE=IAT2.EQ.IAT.AND.IT(1).EQ.0.AND.IT(2).EQ.0.AND.IT(3).EQ.0
          ISP2=ISPECIES(IAT2)
          NORB2=GAUSSORB(IAT2)%NORB
          NIJKB=POTPAR(ISP2)%GAUSSKPRIME%NIJK
          NEB=POTPAR(ISP2)%GAUSSKPRIME%NE
          EB(:NEB)=POTPAR(ISP2)%GAUSSKPRIME%E
          RB(:)=R0(:,IAT2)+MATMUL(RBAS,REAL(IT))
!!$DO IX=1,NX
!!$  ARR(:)=0.D0
!!$  DO IORB=1,NORB2
!!$    DO IE=1,NEB
!!$      CALL GAUSSIAN_3DORB('CARTESIAN',NIJKB,EB(IE),POTPAR(ISP2)%GAUSSKPRIME%C(:,IE,IORB),RGRID(:,IX)+RA(:)-RB(:),SVAR)
!!$      ARR(IORB)=ARR(IORB)+SVAR
!!$    ENDDO
!!$  ENDDO
!!$  DO I=1,NORB1
!!$    IF(TONSITE)FGRID(IX,I,1)=FGRID(IX,I,1)+ARR(I)
!!$    DO J=1,NORB2
!!$     FGRID(IX,I,1)=FGRID(IX,I,1)+ARR(J)*QBAR(J,ISP2)*SBAR(NN)%MAT(J,I)
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
          CALL GAUSSIAN$GAUSSPSI(NIJKA,NEA,EA,RA,NIJKB,NEB,EB(:NEB),RB &
    &                           ,NORB2,POTPAR(ISP2)%GAUSSKPRIME%C &
    &                           ,GCB(:,:,:NORB2))
          DO I=1,NORB1
            IF(TONSITE) GAUSSORB(IAT)%C(:,:,I)=GAUSSORB(IAT)%C(:,:,I)+GCB(:,:,I)
            DO J=1,NORB2
              GAUSSORB(IAT)%C(:,:,I)=GAUSSORB(IAT)%C(:,:,I) &
    &                               +GCB(:,:,J)*QBAR(J,ISP2)*SBAR(NN)%MAT(J,I)
            ENDDO
          ENDDO
        ENDDO
!
!       ========================================================================
!       == project onto central site                                          ==
!       ========================================================================
        CALL LMTO_GAUSSPROJECT(NIJKA,NEA,EA,MAT)
        DO I=1,NORB1
          IND=0
          DO J=1,NEA
            DO K=1,NIJKA
              IND=IND+1
              MAT1(IND,I)=GAUSSORB(IAT)%C(K,J,I)
            ENDDO
          ENDDO
        ENDDO
        CALL LIB$MATRIXSOLVER8(NIJKA*NEA,NIJKA*NEA,NORB1,MAT,MAT2,MAT1)
        DO I=1,NORB1
          IND=0
          DO J=1,NEA
            DO K=1,NIJKA
              IND=IND+1
              GAUSSORB(IAT)%C(K,J,I)=MAT2(IND,I)
            ENDDO
          ENDDO
        ENDDO
!
!!$OPEN(UNIT=10001,FILE='XY.DAT')
!!$DO IX=1,NX
!!$  ARR(:)=0.D0
!!$  DO IORB=1,NORB1
!!$    DO IE=1,NEA
!!$      CALL GAUSSIAN_3DORB('CARTESIAN',NIJKA,EA(IE),GAUSSORB(IAT)%C(:,IE,IORB),RGRID(:,IX),SVAR)
!!$      ARR(IORB)=ARR(IORB)+SVAR
!!$    ENDDO
!!$  ENDDO
!!$  WRITE(10001,*)RGRID(1,IX),ARR(:NORB1),FGRID(IX,:NORB1,1)
!!$ENDDO
!!$CLOSE(10001)
!!$STOP 'XC'
        DEALLOCATE(EA)
        DEALLOCATE(GCB)
        DEALLOCATE(MAT)
        DEALLOCATE(MAT1)
        DEALLOCATE(MAT2)
      ENDDO
                                         CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GAUSSPROJECT(NIJK,NE,E,U)
!     **************************************************************************
!     **  EVALUATE THE OVERLAP MATRIX U_I,J=<G_I|G_J> FOR A SET OF GAUSSIANS  **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **  1) REPLACE THE EIGENVALUE PROBLEM BY A SINGULAR VALUE DECOMPOSITION **
!     **  2) EXPLOIT THAT THE MATRIX MAT IS BLOCK DIAGONAL WITH 8 BLOCKS      **
!     **     (I1,I2=EVEN,J1,J2=EVEN,K1,K2=EVEN)                               **
!     **     (I1,I2=EVEN,J1,J2=EVEN,K1,K2=ODD)                                **
!     **     (I1,I2=EVEN,J1,J2=ODD ,K1,K2=EVEN)                               **
!     **     (I1,I2=EVEN,J1,J2=ODD ,K1,K2=ODD)                                **
!     **     (I1,I2=ODD ,J1,J2=EVEN,K1,K2=EVEN)                               **
!     **     (I1,I2=ODD ,J1,J2=EVEN,K1,K2=ODD)                                **
!     **     (I1,I2=ODD ,J1,J2=ODD ,K1,K2=EVEN)                               **
!     **     (I1,I2=ODD ,J1,J2=ODD ,K1,K2=ODD)                                **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2010******************
      INTEGER(4),INTENT(IN) :: NIJK
      INTEGER(4),INTENT(IN) :: NE
      REAL(8)   ,INTENT(IN) :: E(NE)
      REAL(8)   ,INTENT(OUT):: U(NIJK,NE,NIJK*NE)
      INTEGER(4),PARAMETER  :: IMAXX=20
      LOGICAL(4)            :: TTEST=.TRUE.
      REAL(8)               :: FACTOR(0:IMAXX)
      REAL(8)               :: EFAC(NE,NE,0:3*IMAXX)
      REAL(8)               :: PI
      REAL(8)               :: SVAR,FAC
      INTEGER(4)            :: IMAX
      INTEGER(4)            :: IPOW
      REAL(8)   ,ALLOCATABLE:: MAT(:,:)
      REAL(8)   ,ALLOCATABLE:: EIG(:)
      INTEGER(4)            :: I,J,IE,JE,IPOS1,IPOS2
      INTEGER(4)            :: IND1,I1,J1,K1,IND2,I2,J2,K2
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
!
!     == INT DX X^N EXP(-A*X^2) = FACTOR(N)/SQRT(A)^(N+1) ======================
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJK,I1,J1,K1)
      IMAX=I1
      IF(IMAX.GT.IMAXX) THEN
        CALL ERROR$STOP('LMTO_GAUSSPROJECT')
      END IF
!
!     ==========================================================================
!     == FACTOR(I)=SQRT(2*PI)*DOUBLE FACTORIAL(I) ==============================
!     ==========================================================================
      FACTOR(:)=0.D0
      SVAR=SQRT(2.D0*PI)
      FACTOR(0)=SVAR
      DO I=1,IMAX
        SVAR=SVAR*REAL(2*I-1,KIND=8)
        FACTOR(2*I)=SVAR
      ENDDO
!
!     ==========================================================================
!     == EFAC(I,J,K)=1/SQRT[2(EI+EJ)]^(3+K)
!     ==========================================================================
      DO IE=1,NE
        DO JE=IE,NE
          FAC=1.D0/SQRT(2.D0*(E(IE)+E(JE)))
          SVAR=FAC**3
          DO I=0,3*IMAX
            EFAC(IE,JE,I)=SVAR
            SVAR=SVAR*FAC
          ENDDO
          EFAC(JE,IE,:)=EFAC(IE,JE,:)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == DETERMINE OVERLAP MATRIX OF GAUSSIANS                                ==
!     ==========================================================================
      ALLOCATE(MAT(NIJK*NE,NIJK*NE))
      MAT(:,:)=0.D0
      DO IND1=1,NIJK
        CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND1,I1,J1,K1)
        DO IND2=1,NIJK
          CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND2,I2,J2,K2)
          SVAR=FACTOR(I1+I2)*FACTOR(J1+J2)*FACTOR(K1+K2)
          IF(SVAR.EQ.0.D0) CYCLE
          IPOW=(I1+I2+J1+J2+K1+K2)
          DO IE=1,NE
            IPOS1=NIJK*(IE-1)+IND1
            DO JE=1,NE
              IPOS2=NIJK*(JE-1)+IND2
              MAT(IPOS1,IPOS2)=SVAR*EFAC(IE,JE,IPOW)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!THIS COPIES THE OVERLAP MATRIX BACK WITHOUT INVERSION 
! INTENDED FOR TESTING ONLY!!
      DO IE=1,NE
        DO IND1=1,NIJK
          IPOS1=NIJK*(IE-1)+IND1
          U(IND1,IE,:)=MAT(IPOS1,:)
        ENDDO
      ENDDO

!!$!
!!$!     ==========================================================================
!!$!     == DETERMINE INVERSE OVERLAP AS SUM_I |I><I|                            ==
!!$!     ==========================================================================
!!$      ALLOCATE(EIG(NIJK*NE))
!!$      CALL LIB$DIAGR8(NIJK*NE,MAT,EIG,U)
!!$      DO I=1,NIJK*NE
!!$        U(:,:,I)=U(:,:,I)/SQRT(ABS(EIG(I)))
!!$      ENDDO
!!$      DEALLOCATE(EIG)
!!$      DEALLOCATE(MAT)
!!$!
!!$!     ==========================================================================
!!$!     == TEST INVERSION                                                       ==
!!$!     ==========================================================================
!!$      IF(TTEST) THEN
!!$        DO IPOS1=1,NIJK*NE
!!$          DO IPOS2=1,NIJK*NE
!!$            SVAR=0.D0
!!$            IND1=0
!!$            DO I=1,NIJK
!!$              DO IE=1,NE
!!$                IND1=IND1+1      
!!$                IND2=0
!!$                DO J=1,NIJK
!!$                  DO JE=1,NE
!!$                    IND2=IND2+1
!!$                    SVAR=SVAR+U(I1,IE,IPOS1)*MAT(IND1,IND2)*U(J,JE,IPOS2)
!!$                  ENDDO
!!$                ENDDO
!!$              ENDDO
!!$            ENDDO
!!$            PRINT*,'INVERSION TEST ',IPOS1,IPOS2,SVAR
!!$          ENDDO
!!$        ENDDO
!!$      END IF
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
      SUBROUTINE LMTO$MAKESTRUCTURECONSTANTS()
!     **************************************************************************
!     **  PRODUCES THE SCREENED STRUCTURE CONSTANTS sbar                      **
!     **                                                                      **
!     **  is needed also by ldaplusu via lmto$dolocorb                        **
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
!
!     == STRUCTURE CONSTANTS ARE DETERMINED ONLY ONCE. THIS NEEDS TO BE CHANGED!
!!$      IF(TINISTRUC) THEN
!!$                              CALL TRACE$POP()
!!$        RETURN
!!$      END IF
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
!!$DO I=1,NNS
!!$WRITE(*,FMT='("NNLIST",5I10)')NNLIST(:,I)
!!$ENDDO
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
!       == MEMBERS NN1:NN2 ON THE NEIGHBOLIST BILD THE CLUSTER AROUND ATOM 1  ==
!       == MEMBER NN0 IS THE ONSITE MEMBER                                    ==
        NN1=1
        NN0=0
        DO NN=1,NNS
          IF(NNLIST(1,NN).GT.IAT1)EXIT
          NN2=NN
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
          N=N+(LX+1)**2
          LX1(NN-NN1+1)=LX
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
      SUBROUTINE LMTO_PERIODICMAT_CP(NNS,MAT1,MAT2)
!     **************************************************************************
!     ** COPIES ALL INFORMATION FROM MAT1 INTO MAT2                           **
!     ** POINTER VARIABLES ALLOCATED                                          **
!     **   (USE LMTO_PERIODICMAT_CLEAN IF STILL ALLOCATED)                    **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : PERIODICMAT_TYPE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN)  :: NNS
      TYPE(PERIODICMAT_TYPE),INTENT(IN)  :: MAT1(NNS) 
      TYPE(PERIODICMAT_TYPE),INTENT(OUT) :: MAT2(NNS) 
      INTEGER(4)                         :: NN,N1,N2
!     **************************************************************************
      DO NN=1,NNS
        MAT2(NN)%IAT1=MAT1(NN)%IAT1
        MAT2(NN)%IAT2=MAT1(NN)%IAT2
        MAT2(NN)%IT=MAT1(NN)%IT
        MAT2(NN)%N1=MAT1(NN)%N1
        MAT2(NN)%N2=MAT1(NN)%N2
        N1=MAT2(NN)%N1
        N2=MAT2(NN)%N2
        ALLOCATE(MAT2(NN)%MAT(N1,N2))
        MAT2(NN)%MAT=MAT1(NN)%MAT
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PERIODICMAT_CLEAN(NNS,MAT)
!     **************************************************************************
!     ** ERASES ALL INFORMATION ON MAT AND DEALLOCATED POINTERS               **
!     ** MAT ITSELF IS NOT DEALLOCATED                                        **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : PERIODICMAT_TYPE
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN)    :: NNS
      TYPE(PERIODICMAT_TYPE),INTENT(INOUT) :: MAT(NNS) 
      INTEGER(4)                           :: NN
!     **************************************************************************
      DO NN=1,NNS
        MAT(NN)%IAT1=0
        MAT(NN)%IAT2=0
        MAT(NN)%IT=0
        MAT(NN)%N1=0
        MAT(NN)%N2=0
        DEALLOCATE(MAT(NN)%MAT)
        NULLIFY(MAT(NN)%MAT)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_BONDORIENTATION(IAT1,IAT2,IT,IORIENT)
!     **************************************************************************
!     ** ASSOCIATES AN ORIENTATION TO A BOND VECTOR                           **
!     ** RETURNS -1,0,1                                                       **
!     **                                                                      **
!     ** CAUTION! IT IS NOT PROTECTED AGAINST |IT(I)| > ITX                   **
!     **          IT IS NOT PROTECTED FOR NAT>2000                            **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT1
      INTEGER(4),INTENT(IN)  :: IAT2
      INTEGER(4),INTENT(IN)  :: IT(3)
      INTEGER(4),INTENT(OUT) :: IORIENT  ! CAN HAVE VALUES -1,0,+1
!     **************************************************************************
      IORIENT=IAT2-IAT1
      IF(IORIENT.EQ.0) THEN
        IORIENT=IT(1)
        IF(IORIENT.EQ.0) THEN
          IORIENT=IT(2)
          IF(IORIENT.EQ.0) THEN
            IORIENT=IT(3)
          END IF
        END IF
      END IF 
      IF(IORIENT.NE.0) IORIENT=IORIENT/ABS(IORIENT)
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_EXPANDQBAR(LMXX,NSP_,QBAR)
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : NSP,LNX,LOX,POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMXX
      INTEGER(4),INTENT(IN) :: NSP_
      REAL(8)   ,INTENT(OUT):: QBAR(LMXX,NSP_)
      INTEGER(4)            :: LX
      INTEGER(4)            :: ISP
      INTEGER(4)            :: L
      INTEGER(4)            :: LN
      INTEGER(4)            :: I1,I2
!     **************************************************************************
      IF(NSP_.NE.NSP) THEN
        CALL ERROR$STOP('LMTO_EXPANDQBAR')
      END IF
      LX=MAXVAL(LOX)
      IF((LX+1)**2.GT.LMXX) THEN
        CALL ERROR$STOP('LMTO_EXPANDQBAR')
      END IF
!
!     ==========================================================================
!     == DETERMINE SCREENING PARAMETER QBAR                                   ==
!     ==========================================================================
      QBAR(:,:)=0.D0
      DO ISP=1,NSP
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I1=L**2+1
          I2=(L+1)**2
          QBAR(I1:I2,ISP)=POTPAR(ISP)%QBAR(POTPAR(ISP)%LNSCATT(LN))
        ENDDO
      ENDDO
      RETURN
      END       
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_KTOKBAR(NNS,SBAR,MAP)
!     **************************************************************************
!     ** DETERMINES THE MATRIX MAP THAT PROVIDES THE SCREENED ENVELOPE        **
!     ** FUNCTION KBAR AS SUPERPOSITION OF UNCREENED SUPERPOSITIONS K         **
!     **    |KBAR_I>=\SUM_J |K_J>MAP_JI=\SUM_J |K_J>(1_JI+QBAR_J*SBAR_JI)     **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : PERIODICMAT_TYPE,LOX,NSP,ISPECIES
      IMPLICIT NONE
      INTEGER(4)            ,INTENT(IN) :: NNS
      TYPE(PERIODICMAT_TYPE),INTENT(IN) :: SBAR(NNS) !(NNS)
      TYPE(PERIODICMAT_TYPE),INTENT(INOUT) :: MAP(NNS) !(NNS)
      REAL(8)              ,ALLOCATABLE :: QBAR(:,:) ! (LMXX,NSP)
      INTEGER(4)                        :: NN,I,IAT2,ISP,N1,N2
      INTEGER(4)                        :: LX
!     **************************************************************************
!
!     ==========================================================================
!     ==  COPY SBAR INTO MAP                                                  ==
!     ==========================================================================
      CALL LMTO_PERIODICMAT_CP(NNS,SBAR,MAP)
!
!     ==========================================================================
!     ==  MULTIPLY WITH QBAR ON THE LEFT                                      ==
!     ==========================================================================
      LX=MAXVAL(LOX)
      ALLOCATE(QBAR((LX+1)**2,NSP))
      CALL LMTO_EXPANDQBAR((LX+1)**2,NSP,QBAR)
      DO NN=1,NNS
        IAT2=MAP(NN)%IAT2
        ISP=ISPECIES(IAT2)
        N1=MAP(NN)%N1
        N2=MAP(NN)%N2
        DO I=1,N1
          MAP(NN)%MAT(:,I)=QBAR(:N2,ISP)*MAP(NN)%MAT(:,I)
        ENDDO
      ENDDO
      DEALLOCATE(QBAR)
!
!     ==========================================================================
!     == ADD IDENTITY                                                         ==
!     ==========================================================================
      DO NN=1,NNS
        IF(MAP(NN)%IAT1.NE.MAP(NN)%IAT2) CYCLE
        IF(MAP(NN)%IT(1).NE.0) CYCLE
        IF(MAP(NN)%IT(2).NE.0) CYCLE
        IF(MAP(NN)%IT(3).NE.0) CYCLE
        N1=MAP(NN)%N1
        DO I=1,N1
          MAP(NN)%MAT(I,I)=MAP(NN)%MAT(I,I)+1.D0
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$OVERLAPFULL()
!     **************************************************************************
!     ** DETERMINES OVERLAP MATRIX OF NATURAL TIGHT-BINDING ORBITALS          **
!     ** USING THE COMPLETE MULTICENTER EXPANSION                             **
!     **                                                                      **
!     ** OVERLAP ONLY CALCULATED FOR PAIRS THAT ARE CONNECTED ALSO BY SBAR    **
!     **                                                                      **
!     ** PROGRAMMED VERY INEFFICIENTLY! OPTIMIZE BEFORE PRODUCTION            **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : PERIODICMAT_TYPE,SBAR,OVERLAP,ISPECIES,POTPAR &
     &                       ,GAUSSORBAUG
      IMPLICIT NONE
      TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: MAP(:) !(NNS)
      INTEGER(4)          :: NNS
      INTEGER(4)          :: NN,NN1,NN2
      INTEGER(4)          :: N1,N2,N2A,N2B
      INTEGER(4)          :: IAT1,IAT2,IAT2A,IAT2B,NEA
      INTEGER(4)          :: ISP1,ISP2,ISP2A,ISP2B,NEB
      INTEGER(4)          :: NIJKA,NIJKB
      REAL(8)             :: RA(3),RB(3)
      REAL(8)             :: RBAS(3,3)
      REAL(8),ALLOCATABLE :: R0(:,:)
      REAL(8)             :: DR(3)
      REAL(8)             :: IT(3)
      INTEGER(4)          :: NAT
      INTEGER(4)          :: N1X
      REAL(8)   ,ALLOCATABLE :: BAREOV(:,:)
!     **************************************************************************
                          CALL TRACE$PUSH('LMTO$OVERLAP_FULL')
!
!     ==========================================================================
!     ==  COLLECT ATOMIC STRUCTURE                                            ==
!     ==========================================================================
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      CALL CELL$GETR8A('T0',9,RBAS)
!
!     ==========================================================================
!     ==  CREATE EMPTY ARRAY OVERLAP (INDICES FROM SBAR, MAT IS ZEROED)       ==
!     ==========================================================================
!     == CLEAR ARRAY, BECAUSE NEIGHBORLIST MAY HAVE CHANGED ====================
      IF(ALLOCATED(OVERLAP)) THEN
        DO NN=1,SIZE(OVERLAP)
          DEALLOCATE(OVERLAP(NN)%MAT)
        ENDDO
        DEALLOCATE(OVERLAP)
      END IF
!
!     == REALLOCATE OVERLAP LIST ===============================================
      NNS=SIZE(SBAR)
      ALLOCATE(OVERLAP(NNS))
      DO NN=1,NNS
        IAT1=SBAR(NN)%IAT1
        IAT2=SBAR(NN)%IAT2
        IT=SBAR(NN)%IT
        N1=GAUSSORBAUG(IAT1)%NORB
        N2=GAUSSORBAUG(IAT2)%NORB
        OVERLAP(NN)%IAT1=IAT1
        OVERLAP(NN)%IAT2=IAT2
        OVERLAP(NN)%IT=IT
        OVERLAP(NN)%N1=N1
        OVERLAP(NN)%N2=N2
        ALLOCATE(OVERLAP(NN)%MAT(N1,N2))
      ENDDO
!
!     ==========================================================================
!     ==  ACCUMULATE NTBO OVERLAP MATRIX:GAUSSIAN CONTRIBUTION                ==
!     ==========================================================================
      DO NN=1,NNS
        IAT1=OVERLAP(NN)%IAT1
        IAT2=OVERLAP(NN)%IAT2
        IT=OVERLAP(NN)%IT
        RA(:)=R0(:,IAT1)
        RB(:)=R0(:,IAT2)+MATMUL(RBAS,REAL(IT))
        NIJKA=GAUSSORBAUG(IAT1)%NIJK
        NIJKB=GAUSSORBAUG(IAT2)%NIJK
        NEA=GAUSSORBAUG(IAT1)%NE
        NEB=GAUSSORBAUG(IAT2)%NE
        N1=GAUSSORBAUG(IAT1)%NORB
        N2=GAUSSORBAUG(IAT2)%NORB
        CALL GAUSSIAN$NEWOVERLAP( &
     &                      NIJKA,N1,NEA,GAUSSORBAUG(IAT1)%E,RA,GAUSSORBAUG(IAT1)%C &
     &                     ,NIJKB,N2,NEB,GAUSSORBAUG(IAT2)%E,RB,GAUSSORBAUG(IAT2)%C &
     &                     ,OVERLAP(NN)%MAT)
      ENDDO
!
!     ==========================================================================
!     ==  ADD TAIL-CORRECTION OF THE GAUSSIAN REPRESENTED ORBITALS
!     ==========================================================================
!      CALL LMTO_DELTAOVERLAPFULL()
CALL TRACE$POP()
RETURN
!!$!
!!$!     ==========================================================================
!!$!     == KBAR=K*MAP;     MAP=(1+QBAR*SBAR)                                    ==
!!$!     ==========================================================================
!!$      ALLOCATE(MAP(NNS))
!!$      CALL LMTO_KTOKBAR(NNS,SBAR,MAP)
!!$!
!!$!     ==========================================================================
!!$!     == LIMIT SPHERE OVERLAP TERMS TO ONE CHANNEL PER ANGULAR MOMENTUM       ==
!!$!     == RESULTING OVERLAP MATRIX IS SMALLER AND THE FULL MATRIX IS READILY   ==
!!$!     == OBTAINED BY ADDING ONSITE TERMS ONLY.                                ==
!!$!     == RESULT LIES IN POTPAR(ISP)%SMALL%DOVERLAPKK, ETC.                    ==
!!$!     ==========================================================================
!!$      CALL LMTO_EXPANDSPHEREDO('LM')
!!$
!!$!
!!$!     ==========================================================================
!!$!     ==  ACCUMULATE OVERLAP MATRIX: SPHERE CORRECTION                        ==
!!$!     ==========================================================================
!!$      DO NN=1,NNS
!!$!       ==  <KBAR(IAT1)|KBAR(IAT2)> ============================================
!!$        N1=OVERLAP(NN)%N1
!!$        N2=OVERLAP(NN)%N2
!!$        IAT1=OVERLAP(NN)%IAT1
!!$        IAT2=OVERLAP(NN)%IAT2
!!$        ISP1=ISPECIES(IAT1)
!!$        ISP2=ISPECIES(IAT2)
!!$!
!!$!       == ONSITE TERM OF HEADS: <KOMEGA(AT1)|KOMEGA(AT2)> =====================
!!$        IF(IAT1.EQ.IAT2) THEN
!!$          IT=OVERLAP(NN)%IT
!!$          IF(IT(1).EQ.0.AND.IT(2).EQ.0.AND.IT(3).EQ.0) THEN
!!$            OVERLAP(NN)%MAT=OVERLAP(NN)%MAT+POTPAR(ISP1)%SMALL%DOVERLAPKK
!!$          END IF
!!$        END IF
!!$!
!!$        DO NN1=1,NNS
!!$          IF(SBAR(NN1)%IAT1.NE.IAT1) CYCLE
!!$!         ==  <KBAR(IAT1)|=SUM(AT2A):MAP(AT1,AT2A)<KBARPRIME(AT2A)| ============
!!$!         ==  <KBAR(IAT1)|=<KBAROMEGA(AT1)|                          ===========
!!$!         ==              +SUM(AT2A):SBAR(AT1,AT2A)<JBAROMEGA(AT2A)| ===========
!!$          IAT2A=SBAR(NN1)%IAT2
!!$          ISP2A=ISPECIES(IAT2A)
!!$          DO NN2=1,NNS
!!$            IF(SBAR(NN2)%IAT1.NE.IAT2) CYCLE
!!$!           ==  <KBAR(IAT2)|=SUM(AT2B):MAP(AT2,AT2B)<KBARPRIME(AT2B)| ==========
!!$!           ==  <KBAR(IAT2)|=<KBAROMEGA(AT2)|                          =========
!!$!           ==              +SUM(AT2B):SBAR(AT2,AT2B)<JBAROMEGA(AT2B)| =========
!!$            IAT2B=SBAR(NN2)%IAT2
!!$            ISP2B=ISPECIES(IAT2B)
!!$!
!!$!           == TAIL-TAIL OVERLAP: SBAR<JBAR|JBAR>SBAR ==========================
!!$            IF(IAT2A.EQ.IAT2B) THEN
!!$              IT=-SBAR(NN1)%IT+OVERLAP(NN)%IT+SBAR(NN2)%IT
!!$              IF(IT(1).EQ.0.AND.IT(2).EQ.0.AND.IT(3).EQ.0) THEN
!!$                OVERLAP(NN)%MAT=OVERLAP(NN)%MAT &
!!$    &                   +MATMUL(TRANSPOSE(SBAR(NN1)%MAT) &
!!$                          ,MATMUL(POTPAR(ISP2A)%SMALL%DOVERLAPJJ,SBAR(NN2)%MAT))
!!$              END IF
!!$            END IF
!!$          ENDDO ! END OF LOOP OVER NN2
!!$!
!!$!         == TAIL-HEAD OVERLAP:  SBAR<JBAR|K> ==================================
!!$          IF(IAT2A.EQ.IAT2) THEN
!!$            IT=-SBAR(NN1)%IT+OVERLAP(NN)%IT
!!$            IF(IT(1).EQ.0.AND.IT(2).EQ.0.AND.IT(3).EQ.0) THEN
!!$              OVERLAP(NN)%MAT=OVERLAP(NN)%MAT &
!!$     &           +TRANSPOSE(MATMUL(POTPAR(ISP2)%SMALL%DOVERLAPKJ,SBAR(NN1)%MAT))
!!$            END IF
!!$          END IF
!!$        ENDDO  ! END OF LOOP OVER NN1
!!$        DO NN2=1,NNS
!!$          IF(SBAR(NN2)%IAT1.NE.IAT2) CYCLE
!!$          IAT2B=SBAR(NN2)%IAT2
!!$!
!!$!         == HEAD-TAIL, OVERLAP:   <K|JBAR>SBAR ================================
!!$          IF(IAT1.EQ.IAT2B) THEN
!!$            IT=OVERLAP(NN)%IT+SBAR(NN2)%IT
!!$            IF(IT(1).EQ.0.AND.IT(2).EQ.0.AND.IT(3).EQ.0) THEN
!!$              OVERLAP(NN)%MAT=OVERLAP(NN)%MAT &
!!$    &                 +MATMUL(POTPAR(ISP1)%SMALL%DOVERLAPKJ,SBAR(NN2)%MAT)
!!$            END IF
!!$          END IF
!!$        ENDDO  ! END OF LOOP OVER NN2
!!$      ENDDO  ! END OF LOOP OVER NN
!!$     
                          CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DELTAOVERLAPFULL()
!     **************************************************************************
!     ** CALCULATES THE CORRECTION FOR THE DIFFERENCE BETWEEN AEPHI AND NLPHI **
!     ** IN THE TAILS OF THE GAUSS EXPANDED WAVE FUNCTIONS                    **
!     **                                                                      **
!     ** THIS CODE IS NOT WELL CHECKED AND MAY HAVE ERRORS. THE RESULTS       **
!     ** INDICATE THAT THE CORRECTION IS NEGLEGIBLE.  THE IMPORTANCE OF THIS  **
!     ** TERM CAN ALSO BE CHECKED FROM AN FCC HYDROGEN CRYSTAL, BECAUSE HERE  **
!     ** AEPHI=NLPHI AND THE CORRECTION IS ZERO.                              **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : NSP,LNX,LOX,SBARLI1 &
     &                       ,PERIODICMAT_TYPE,SBAR,OVERLAP,ISPECIES,POTPAR &
     &                       ,GAUSSORB
      IMPLICIT NONE
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      INTEGER(4)             :: LNXX    !=X(LNX)
      INTEGER(4)             :: LXX    !=X(L)
      INTEGER(4)             :: LNX1    !=LNX(ISP)
      INTEGER(4)             :: NNS
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDELTA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOTDELTA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOTDELTA(:,:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      INTEGER(4),ALLOCATABLE :: LNS(:,:)
      INTEGER(4),ALLOCATABLE :: LNSOFL(:,:)
      INTEGER(4),ALLOCATABLE :: ISCATT(:)
      LOGICAL(4)             :: TONSITEA,TONSITEB
      REAL(8)                :: IT(3),ITA(3),ITB(3)
      REAL(8)                :: SVAR,SVAR1
      INTEGER(4)             :: ISP,ISP1,ISP2,LN1A,LN2A,NN,NNA,NNB,LX
      INTEGER(4)             :: IAT,IAT1,IAT2,IAT1A,IAT2A,IAT1B,IAT2B
      INTEGER(4)             :: LN,L,IM,LM
      INTEGER(4)             :: LMN1,LN1,L1,IM1,LM1
      INTEGER(4)             :: LMN2,LN2,L2,IM2,LM2
!     **************************************************************************
      LXX=-1
      DO ISP=1,NSP
        LXX=MAX(LXX,MAXVAL(LOX(:LNX(ISP),ISP)))
      ENDDO
      LNXX=MAXVAL(LNX)
      ALLOCATE(LNSOFL(LXX+1,NSP))
      ALLOCATE(LNS(LNXX,NSP))
      LNSOFL(:,:)=0
      ALLOCATE(ISCATT(LNXX))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4A('ISCATT',LNXX,ISCATT)
!       == CONSTRUCT MAPPING ONTO VALENCE CHANNELS FOR SCATTERING PART
        DO LN1=1,LNX(ISP)
          LNS(LN1,ISP)=0
          L=LOX(LN1,ISP)
          DO LN=1,LNX(ISP)
            IF(LOX(LN,ISP).NE.LOX(LN1,ISP)) CYCLE
            IF(ISCATT(LN).GT.0) CYCLE
            LNS(LN1,ISP)=LN
            LNSOFL(L+1,ISP)=LN
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  DETERMINE MATRIX ELEMENTS OF PARTIAL WAVES                          ==
!     ==========================================================================
      ALLOCATE(AEPHIDELTA(LNXX,LNXX,NSP))
      ALLOCATE(AEPHIDOTDELTA(LNXX,LNXX,NSP))
      ALLOCATE(NLPHIDOTDELTA(LNXX,LNXX,NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        CALL SETUP$GETI4A('ISCATT',LNXX,ISCATT)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        ALLOCATE(AUX(NR))
        LNX1=LNX(ISP)
        ALLOCATE(AEPHI(NR,LNX1))
        ALLOCATE(AEPHIDOT(NR,LNX1))
        ALLOCATE(NLPHIDOT(NR,LNX1))
        CALL SETUP$GETR8A('AEPHI',NR*LNX1,AEPHI)
        CALL SETUP$GETR8A('QPHIDOT',NR*LNX1,NLPHIDOT)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX1,AEPHIDOT)
!        LX=MAXVAL(LOX(:,ISP))
!        LMX=SBARLI1(LX+1,ISP)+2*LX
!       == CONSTRUCT MATRIX ELEMENTS
CALL SETUP_WRITEPHI('XX1.DAT',GID,NR,1,R**2*AEPHI(:,1)*(AEPHIDOT(:,1)-NLPHIDOT(:,1)))
CALL SETUP_WRITEPHI('XX2.DAT',GID,NR,1,R**2*AEPHI(:,2)*(AEPHIDOT(:,1)-NLPHIDOT(:,1)))
CALL SETUP_WRITEPHI('XX3.DAT',GID,NR,1,R**2*AEPHI(:,3)*(AEPHIDOT(:,3)-NLPHIDOT(:,3)))
        AEPHIDELTA=0.D0
        AEPHIDOTDELTA=0.D0
        NLPHIDOTDELTA=0.D0
        DO LN1=1,LNX1
          LN1A=LNS(LN1,ISP)
          DO LN2=1,LNX1
            LN2A=LNS(LN2,ISP)
            IF(LOX(LN1,ISP).NE.LOX(LN2,ISP)) CYCLE
!
            AUX(:)=R(:)**2*AEPHI(:,LN1)*(AEPHIDOT(:,LN2A)-NLPHIDOT(:,LN2A))
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            SVAR=SVAR*POTPAR(ISP)%JBARTOPHIDOT(LN2A)
            AEPHIDELTA(LN1,LN2,ISP)=POTPAR(ISP)%KTOPHI(LN1)*SVAR
!
            AUX(:)=R(:)**2*AEPHIDOT(:,LN1A)*(AEPHIDOT(:,LN2A)-NLPHIDOT(:,LN2A))
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            SVAR=SVAR*POTPAR(ISP)%JBARTOPHIDOT(LN2A)
            AEPHIDELTA(LN1,LN2,ISP)=AEPHIDELTA(LN1,LN2,ISP) &
     &                             +POTPAR(ISP)%KTOPHIDOT(LN1A)*SVAR
            AEPHIDOTDELTA(LN1,LN2,ISP)=POTPAR(ISP)%JBARTOPHIDOT(LN1A)*SVAR
!
            AUX(:)=R(:)**2*NLPHIDOT(:,LN1A)*(AEPHIDOT(:,LN2A)-NLPHIDOT(:,LN2A))
            CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
            SVAR=POTPAR(ISP)%JBARTOPHIDOT(LN1A)*SVAR*POTPAR(ISP)%JBARTOPHIDOT(LN2A)
            NLPHIDOTDELTA(LN1,LN2,ISP)=SVAR
          ENDDO
WRITE(*,FMT='("AEPHIDELTA   : ",30F10.5)') AEPHIDELTA(LN1,:,ISP)
WRITE(*,FMT='("AEPHIDOTDELTA: ",30F10.5)') AEPHIDOTDELTA(LN1,:,ISP)
WRITE(*,FMT='("NLPHIDOTDELTA: ",30F10.5)') NLPHIDOTDELTA(LN1,:,ISP)
        ENDDO
        DEALLOCATE(R)
        DEALLOCATE(AUX)
        DEALLOCATE(AEPHI)
        DEALLOCATE(AEPHIDOT)
        DEALLOCATE(NLPHIDOT)
      ENDDO
!
!     ==========================================================================
!     ==  CONVERT TO K AND JBAR
!     ==========================================================================
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IAT1=OVERLAP(NN)%IAT1
        IAT2=OVERLAP(NN)%IAT2
        ISP1=ISPECIES(IAT1)
        ISP2=ISPECIES(IAT2)
        DO NNA=1,NNS
          IAT1A=SBAR(NNA)%IAT1
          IF(IAT1A.NE.IAT1) CYCLE
          IAT2A=SBAR(NNA)%IAT2
          ITA(:)=SBAR(NNA)%IT
          TONSITEA=(IAT1A.EQ.IAT2A).AND.(ITA(1).EQ.0) &
     &                             .AND.(ITA(2).EQ.0).AND.(ITA(3).EQ.0)
          DO NNB=1,NNS
            IAT1B=SBAR(NNB)%IAT1
            IF(IAT1B.NE.IAT2) CYCLE
            IAT2B=SBAR(NNB)%IAT2
            ITB(:)=SBAR(NNB)%IT
            TONSITEB=(IAT1B.EQ.IAT2B).AND.(ITB(1).EQ.0) &
     &                               .AND.(ITB(2).EQ.0).AND.(ITB(3).EQ.0)
            IT=-ITA+OVERLAP(NN)%IT+ITB
            IF(IAT2A.NE.IAT2B.OR.IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
            IAT=IAT2B
            ISP=ISPECIES(IAT)
!
            IF(TONSITEA.AND.TONSITEB) CYCLE
            IF(TONSITEA) THEN
!              OVERLAP(NN)%MAT(:,:)=MAT-<K-SBAR*JBAR||DELTA_JBAR>SBAR
               LMN1=0
               DO LN1=1,LNX(ISP1)
                 L1=LOX(LN1,ISP1)
                 DO IM1=1,2*L1+1
                   LMN1=LMN1+1
                   LM1=SBARLI1(L1+1,ISP1)+IM1-1
                   LMN2=0
                   DO LN2=1,LNX(ISP2)
                     L2=LOX(LN2,ISP2)
                     DO IM2=1,2*L2+1
                       LMN2=LMN2+1
                       LM2=SBARLI1(L2+1,ISP)+IM2-1
                       LX=MAXVAL(LOX(:LNX(ISP),ISP))
                       SVAR=0.D0
                       DO L=0,LX
                         LN=LNSOFL(L+1,ISP)
                         LM=SBARLI1(L+1,ISP)-1+IM1
                         SVAR=SVAR+AEPHIDELTA(LN1,LN,ISP)*SBAR(NNB)%MAT(LM,LM2)
                       ENDDO
                       OVERLAP(NN)%MAT(LMN1,LMN2)=OVERLAP(NN)%MAT(LMN1,LMN2)+SVAR
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO  
            ELSE IF(TONSITEB) THEN
!              OVERLAP(NN)%MAT(:,:)=MAT-SBAR*<DELTA_JBAR|K-JBAR*SBAR>
               LMN1=0
               DO LN1=1,LNX(ISP1)
                 L1=LOX(LN1,ISP1)
                 DO IM1=1,2*L1+1
                   LMN1=LMN1+1
                   LM1=SBARLI1(L1+1,ISP1)+IM1-1
                   LMN2=0
                   DO LN2=1,LNX(ISP2)
                     L2=LOX(LN2,ISP2)
                     DO IM2=1,2*L2+1
                       LMN2=LMN2+1
                       LM2=SBARLI1(L2+1,ISP)+IM2-1
                       LX=MAXVAL(LOX(:LNX(ISP),ISP))
                       SVAR=0.D0
                       DO L=0,LX
                         LN=LNSOFL(L+1,ISP)
                         LM=SBARLI1(L+1,ISP)-1+IM2
                         SVAR=SVAR+SBAR(NNB)%MAT(LM,LM1)*AEPHIDELTA(LN2,LN,ISP)
                       ENDDO
                       OVERLAP(NN)%MAT(LMN1,LMN2)=OVERLAP(NN)%MAT(LMN1,LMN2)+SVAR
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO  
            ELSE   ! TAIL-TAIL CONTRIBUTION
!              OVERLAP(NN)%MAT(:,:)=MAT-SBAR(<JBAR|DELTA_JBAR>+<DELTA_JBAR|JBAR>
!                                           -<DELTA_JBAR|DELTA_JBAR>)SBAR
               LMN1=0
               DO LN1=1,LNX(ISP1)
                 L1=LOX(LN1,ISP1)
                 DO IM1=1,2*L1+1
                   LMN1=LMN1+1
                   LM1=SBARLI1(L1+1,ISP1)+IM1-1
                   LMN2=0
                   DO LN2=1,LNX(ISP2)
                     L2=LOX(LN2,ISP2)
                     DO IM2=1,2*L2+1
                       LMN2=LMN2+1
                       LM2=SBARLI1(L2+1,ISP)+IM2-1
                       LX=MAXVAL(LOX(:LNX(ISP),ISP))
                       SVAR=0.D0
                       DO L=0,LX
                         LN=LNSOFL(L+1,ISP)
                         LM=SBARLI1(L+1,ISP)-1
                         SVAR1=0.D0
                         DO IM=1,2*L+1
                           LM=LM+1
                           SVAR1=SVAR1+SBAR(NNA)%MAT(LM,LM1)*SBAR(NNB)%MAT(LM,LM2)
                         ENDDO
                         SVAR=SVAR+SVAR1*(AEPHIDOTDELTA(LN,LN,ISP)+NLPHIDOTDELTA(LN,LN,ISP))
                       ENDDO
                       OVERLAP(NN)%MAT(LMN1,LMN2)=OVERLAP(NN)%MAT(LMN1,LMN2)+SVAR
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO  
            END IF
          ENDDO
        ENDDO
      ENDDO
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$TESTENERGY(lmnxx_,ndimd_,nat_,denmat_)
      USE LMTO_MODULE, ONLY : ton,GAUSSORB,GAUSSORB_T,GAUSSORBAUG
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),intent(in) :: lmnxx_
      INTEGER(4),intent(in) :: ndimd_
      INTEGER(4),intent(in) :: nat_
      complex(8),intent(in) :: denmat_(lmnxx_,lmnxx_,ndimd_,nat_)
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IAT
      LOGICAL(4),SAVE       :: TFIRSTENERGY=.TRUE.
!     **************************************************************************
      IF(.NOT.TON) RETURN
!!$      IF(.NOT.TFIRSTENERGY) THEN
!!$        CALL ERROR$MSG('TEST STOP BEFORE SECOND ITERATION')
!!$        CALL ERROR$STOP('LMTO$TESTENERGY')
!!$      END IF
!!$      IF(TFIRSTENERGY) TFIRSTENERGY=.FALSE.
!RETURN
      WRITE(*,FMT='(82("="),T30," TESTENERGY START ")')
!      CALL LMTO$REPORTPOTBAR(6)
!      CALL LMTO$REPORTSBAR(6)
!
!     ==========================================================================
!     == CONSTRUCT TIGHT-BINDING ORBITALS IN GAUSSIAN REPRESENTATION          ==
!     ==========================================================================
PRINT*,'DOING LMTO_NTBFROMTAILEDKJ.....'
CALL TIMING$CLOCKON('NTBOFROMTAILEDKJ')
      CALL LMTO_NTBOFROMTAILEDKJ()
CALL TIMING$CLOCKOFF('NTBOFROMTAILEDKJ')
PRINT*,'.....LMTO_NTBFROMTAILEDKJ DONE'
!
PRINT*,'DOING LMTO_NTBOFROMKPRIME.....'
CALL TIMING$CLOCKON('NTBOFROMKPRIME')
      CALL LMTO_NTBOFROMKPRIME() !THIS IS THE CORRECT SUPERPOSITION
CALL TIMING$CLOCKOFF('NTBOFROMKPRIME')
PRINT*,'..... LMTO_NTBOFROMKPRIME DONE'
!!$ !== The following does not work
!!$nat=size(gaussorb_t)
!!$call LMTO_CPGAUSSORB(NAT,gaussorb_t,gaussorb)


!!$IF(.NOT.TFIRSTENERGY) THEN
!!$  CALL ERROR$MSG('TEST STOP BEFORE SECOND ITERATION')
!!$  CALL ERROR$STOP('LMTO$TESTENERGY')
!!$END IF
!
!     ==========================================================================
!     ==  ADD AUGMENTATION                                                    ==
!     ==========================================================================
PRINT*,'DOING LMTO_NTBOAUGMENT.....'
CALL TIMING$CLOCKON('NTBOAUGMENT')
      CALL LMTO_NTBOAUGMENT()
CALL TIMING$CLOCKOFF('NTBOAUGMENT')
!     == FROM HERE ON USE GAUSSORBAUG...
PRINT*,'..... LMTO_NTBOAUGMENT DONE'
!
!     ==========================================================================
!     ==  TEST COEFFICIENTS OF TIGHT-BINDING ORBITALS                         ==
!     ==========================================================================
PRINT*,' BEFORE TESTNTBO.....'
!      CALL LMTO_TESTNTBO()
PRINT*,'....... TESTNTBO DONE'
!
!     ==========================================================================
!     ==  CALCULATE OVERLAP MATRIX                                            ==
!     ==========================================================================
!!$PRINT*,'DOING LMTO$OVERLAPFULL.....'
!!$CALL TIMING$CLOCKON('OVERLAPPFULL')
!!$      CALL LMTO$OVERLAPFULL()
!!$CALL TIMING$CLOCKOFf('OVERLAPPFULL')
!!$PRINT*,'..... LMTO$OVERLAPFULL DONE'
!      CALL LMTO$REPORTOVERLAP(6)
!      CALL LMTO_PLOTRADIAL()
!
!     ==========================================================================
!     ==  
!     ==========================================================================
!
!     ==========================================================================
!     == TESTS THE OVERLAP MATRIX OF NATURAL TIGHT-BINDING ORBITALS           ==
!     == BY ESTIMATING THE OVERLAP BETWEEN KOHN SHAM WAVE FUNCTIONS USING THE ==
!     == COEFFICIENTS IN NTBS AND THEIR OVERLAP MATRIX                        ==
!     ==========================================================================
PRINT*,' BEFORE TESTOVERLAP.....'
!      CALL LMTO_TESTOVERLAP()
PRINT*,'....... TESTOVERLAP DONE'
!
!     ==========================================================================
!     ==  MAP ORBITALS TO A RADIAL GRID TIME SPHERICAL HARMONICS              ==
!     ==  THIS WILL BE EXPLOITED IN OMNSIDTEU                                 ==
!     ==========================================================================
PRINT*,' BEFORE LMTO_GAUSSORBTOLM.....'
CALL TIMING$CLOCKON('GAUSSORBTOLM')
      CALL LMTO_GAUSSORBTOLM()
CALL TIMING$CLOCKOFF('GAUSSORBTOLM')
PRINT*,'.......LMTO_GAUSSORBTOLM DONE'
!
!     ==========================================================================
!     ==  ONSITE UTENSOR                                                      ==
!     ==========================================================================
PRINT*,'BEFORE LMTO_UTENSORLAYOUT......'
CALL TIMING$CLOCKON('UTENSORLAYOUT')
      CALL LMTO_UTENSORLAYOUT()
CALL TIMING$CLOCKOFF('UTENSORLAYOUT')
PRINT*,'........LMTO_UTENSORLAYOUT DONE'
!
PRINT*,'BEFORE LMTO_ONSITEU............'
CALL TIMING$CLOCKON('ONSITEU')
      CALL LMTO_ONSITEU()
CALL TIMING$CLOCKOFF('ONSITEU')
PRINT*,' .............LMTO_ONSITEU DONE'
!
!     ==========================================================================
!     ==  
!     ==========================================================================
CALL TIMING$CLOCKON('NTBODENMAT')
      CALL LMTO_NTBODENMAT()
CALL TIMING$CLOCKOFF('NTBODENMAT')
!      CALL LMTO_TESTDENMAT()
!      CALL LMTO_TESTDENMAT_1cdenmat(lmnxx_,ndimd_,nat_,denmat_)
!
!     ==========================================================================
!     ==  CALCULATE ENERGY                                                    ==
!     ==========================================================================
PRINT*,'BEFORE LMTO_ENERGYTEST......'
CALL TIMING$CLOCKON('ENERGYTEST')
!      call LMTO_SIMPLEENERGYTEST()
      call LMTO_SIMPLEENERGYTEST2()
!      CALL LMTO_ENERGYTEST()
CALL TIMING$CLOCKOFF('ENERGYTEST')
PRINT*,'.......LMTO_ENERGYTEST DONE'
!
!     ==========================================================================
!     ==  CONVERT HAMIL INTO HTBC
!     ==========================================================================
CALL TIMING$CLOCKON('NTBODENMATDER')
      CALL LMTO_NTBODENMATDER()
CALL TIMING$CLOCKOFF('NTBODENMATDER')
!!$PRINT*,'WARNING!!! LEAVING LMTO$TESTENERGY'
!!$STOP 'FORCED BEFORE LEAVING LMTO$TESTENERGY'
      IF(TFIRSTENERGY) TFIRSTENERGY=.FALSE.
RETURN
!
!     ==========================================================================
!     ==  
!     ==========================================================================
!PRINT*,' BEFORE LMTO_FOURCENTERGAUSS.....'
!      CALL LMTO_FOURCENTERGAUSS()
!PRINT*,'....... LMTO_FOURCENTERGAUSS DONE'
!
!     ==========================================================================
!     ==  
!     ==========================================================================
PRINT*,'DOING LMTO_PLOTLOCORB.....'
      CALL ATOMLIST$NATOM(NAT)
      DO IAT=1,NAT
         CALL LMTO_PLOTLOCORB(IAT)
      ENDDO
PRINT*,'..... LMTO_PLOTLOCORB DONE'
CALL ERROR$MSG('FORCED STOP')
CALL ERROR$STOP('LMTO$TESTENERGY')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$PLOTWAVE(NFIL,IDIM0,IB0,IKPT0,ISPIN0,NR1,NR2,NR3)
!     **************************************************************************
!     **************************************************************************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : ton,GAUSSORB,ISPECIES
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
      if(.not.ton) return
      PI=4.D0*ATAN(1.D0)
      CI2PI=(0.D0,1.D0)*2.D0*PI
!
!     ==========================================================================
!     ==  COLLECT DATA                                                        ==
!     ==========================================================================
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
WRITE(*,FMT='(100("="),T10," EIGENVECTORS ")')
DO I=1,THIS%NB
  WRITE(*,FMT='(100("(",2F10.5,") "))')THIS%EIGVEC(:,I)
ENDDO
WRITE(*,FMT='(100("="),T10," TBC ")')
DO I=1,THIS%NBH
  WRITE(*,FMT='(100("(",2F10.5,") "))')THIS%TBC(1,I,:)
ENDDO
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
              NIJK=GAUSSORB(IAT)%NIJK
              NE=GAUSSORB(IAT)%NE
              NORB=GAUSSORB(IAT)%NORB
              ALLOCATE(COEFF(NIJK,NE))
              ALLOCATE(COEFF1(NIJK))
              COEFF(:,:)=(0.D0,0.D0)
              DO IORB=1,NORB
                IPRO=IPRO+1
                COEFF(:,:)=COEFF(:,:)+GAUSSORB(IAT)%C(:,:,IORB)*CVEC(IPRO)
              ENDDO
              COEFF(:,:)=COEFF(:,:)*EIKT
              DO IP=1,NP
                X=P(1,IP)-R0(1,IAT)-T(1)
                Y=P(2,IP)-R0(2,IAT)-T(2)
                Z=P(3,IP)-R0(3,IAT)-T(3)
                R2=X**2+Y**2+Z**2
                IF(R2.GT.10.D0**2) CYCLE

                COEFF1(:)=(0.D0,0.D0)
                DO IE=1,NE
                  COEFF1(:)=COEFF1(:)+COEFF(:,IE)*EXP(-GAUSSORB(IAT)%E(IE)*R2)
                ENDDO
                CSVAR=(0.D0,0.D0)
                DO IND=1,NIJK
                  CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I,J,K)
                  CSVAR=CSVAR+(X**I)*(Y**J)*(Z**K)*COEFF1(IND)
                ENDDO
                WAVE(IP)=WAVE(IP)+CSVAR
              ENDDO
              DEALLOCATE(COEFF)
              DEALLOCATE(COEFF1)
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
      SUBROUTINE LMTO_SIMPLEENERGYTEST2()
!     **************************************************************************
!     **  work out the energy using the local approximation                   **
!     **  tailed partial waves                                                **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,DENMAT,HAMIL,lnx,lox,potpar
      IMPLICIT NONE
      INTEGER(4)            :: NNU
      INTEGER(4)            :: NND
      INTEGER(4)            :: NNH
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IND,INU,INH
      INTEGER(4)            :: lmnx
      INTEGER(4)            :: NDIMD
      INTEGER(4)            :: IT(3)
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: lnx1,lmrx,lrx
      INTEGER(4)            :: lnxt,lmnxt
      INTEGER(4),allocatable:: loxt(:)
      INTEGER(4)            :: gid
      INTEGER(4)            :: nr
      logical(4),allocatable:: torb(:)
      REAL(8)   ,ALLOCATABLE:: dt(:,:,:)
      REAL(8)   ,ALLOCATABLE:: ht(:,:,:)
      REAL(8)   ,ALLOCATABLE:: aecore(:)
      REAL(8)   ,ALLOCATABLE:: ulittle(:,:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: U(:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: D(:,:,:)
      REAL(8)   ,ALLOCATABLE:: H(:,:,:)
      REAL(8)               :: EH,EX,EXTOT,EHTOT,EHDC,EXDC
      INTEGER(4)            :: NN,IAT,I,J,K,L,IS,IAT1,IAT2,isp
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,PARAMETER  :: HFWEIGHT=0.25D0
character(128) :: string
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_ENERGYTEST2')
print*,'============ energytest2 ============================='
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
        ISP=ISPECIES(IAT)
        LMNX=DENMAT(IND)%N1
        NDIMD=DENMAT(IND)%N3
!
!       ========================================================================
!       == CALCULATE U-TENSOR                                                 ==
!       ========================================================================
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$GETi4(GID,'NR',NR)
        ALLOCATE(AECORE(NR))
        CALL SETUP$GETR8A('AECORE',nr,AECORE)
        CALL SETUP$GETI4('LMRX',LMRX)
        LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
        ALLOCATE(D(LMNX,LMNX,NDIMD))
        ALLOCATE(H(LMNX,LMNX,NDIMD))
        D=DENMAT(IND)%MAT
        h(:,:,:)=0.d0
!
        lnxt=potpar(isp)%tailed%lnx
        lmnxt=potpar(isp)%tailed%lmnx
        allocate(loxt(lnxt))
        loxt=potpar(isp)%tailed%lox
        allocate(dt(lmnxt,lmnxt,ndimd))
        allocate(ht(lmnxt,lmnxt,ndimd))
        call LMTO_blowupdenmat(iat,ndimd,lmnx,d,lmnxt,dt)
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
        EH=0.D0
        EX=0.D0
        Ht(:,:,:)=0.D0
        DO I=1,LMNXt
          DO J=1,LMNXt
            DO K=1,LMNXt
              DO L=1,LMNXt
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
                  EX=EX-0.25D0*U(I,J,K,L)*Dt(K,J,IS)*Dt(L,I,IS)
                  Ht(K,J,IS)=Ht(K,J,IS)-0.25D0*U(I,J,K,L)*Dt(L,I,IS) 
                  Ht(L,I,IS)=Ht(L,I,IS)-0.25D0*U(I,J,K,L)*Dt(K,J,IS) 
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        call LMTO_blowdownht(iat,ndimd,lmnxt,ht,lmnx,h)
        HAMIL(INH)%MAT=H
        EXTOT=EXTOT+EX
print*,'exact exchange energy for atom=',iat,ex
!
!       ========================================================================
!       == ADD CORE VALENCE EXCHANGE                                          ==
!       ========================================================================
        CALL LMTO_CVX(ISP,LMNX,EX,D(:,:,1),H(:,:,1))
        EXTOT=EXTOT+EX
        HAMIL(INH)%MAT(:,:,1)=HAMIL(INH)%MAT(:,:,1)+H(:,:,1)
print*,'core valence exchange energy for atom=',iat,ex
!
!       ========================================================================
!       == DOUBLE COUNTING CORRECTION (EXCHANGE ONLY)                         ==
!       ========================================================================
! this is the time consumin part of energytest
CALL TIMING$CLOCKON('ENERGYTEST:DC')      
        call lmto_simpledc(GID,NR,lmnxt,LNXt,LOXt,potpar(isp)%tailed%aef &
     &                    ,LRX,AECORE,Dt,Ex,Ht)
        call LMTO_blowdownht(iat,ndimd,lmnxt,ht,lmnx,h)
        EXTOT=EXTOT-EX
print*,'double counting correction energy for atom=',iat,ex
        HAMIL(INH)%MAT=HAMIL(INH)%MAT-H
CALL TIMING$CLOCKOFF('ENERGYTEST:DC')      
        deallocate(ht)
        deallocate(dt)
        DEALLOCATE(U)
        DEALLOCATE(H)
        DEALLOCATE(D)
        deallocate(aecore)
        deallocate(loxt)
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
        WRITE(*,FMT='(82("="),T10," ONLY ONSITE TERMS written ")')
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
      SUBROUTINE LMTO_TAILEDU(IAT,LMNX,U)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR,LNX,LOX,SBAR,ISPECIES
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(out):: U(LMNX,LMNX,LMNX,LMNX)
      INTEGER(4)            :: NNS
      INTEGER(4)            :: ISP
      INTEGER(4)            :: NN
      INTEGER(4)            :: LMNXt
      REAL(8)               :: SBARLOC(lmnx,lmnx)
      REAL(8)   ,ALLOCATABLE:: UT1(:,:,:,:)
      INTEGER(4)            :: IORB1,IORB2,IORB2T
!     **************************************************************************
      ISP=ISPECIES(IAT)
!
!     ======================================================================
!     ==  COLLECT LOCAL STRUCTURE CONSTANTS                               ==
!     ======================================================================
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
        IF(SBAR(NN)%IAT2.NE.IAT) CYCLE
        IF(SUM(SBAR(NN)%IT(:)**2).NE.0) CYCLE
        if(LMNX.ne.sbar(nn)%n1) then
          call error$msg('inconsistent array sizes')
          call error$i4val('iat',iat)
          call error$i4val('lmnx',lmnx)
          call error$i4val('sbar%n1',sbar(nn)%n1)
          call error$i4val('sbar%n2',sbar(nn)%n2)
          call error$stop('LMTO_TAILEDU')
        end if
        SBARLOC(:,:)=SBAR(NN)%MAT(:,:)
        EXIT
      ENDDO
!
!     ==========================================================================
!     == U-TENSOR OF ORBITALS                                                 ==
!     ==========================================================================
      LMNXT=POTPAR(ISP)%TAILED%LMNX
      ALLOCATE(UT1(LMNXT,LMNXT,LMNXT,LMNXT))
      UT1=POTPAR(ISP)%TAILED%U
!
!     ======================================================================
!     ==  TRANSFORM U TENSOR TO NON-SPHERICAL ORBITALS                    ==
!     ======================================================================
      DO IORB1=1,LMNX
        DO IORB2=1,LMNX
          IORB2T=POTPAR(ISP)%TAILED%LMNDOT(IORB2)
          UT1(IORB1,:,:,:)=UT1(IORB1,:,:,:) &
     &                    -UT1(IORB2T,:,:,:)*SBARLOC(IORB2,IORB1)
        ENDDO
      ENDDO
      DO IORB1=1,LMNX
        DO IORB2=1,LMNX
          IORB2T=POTPAR(ISP)%TAILED%LMNDOT(IORB2)
          UT1(:,IORB1,:,:)=UT1(:,IORB1,:,:) &
     &                    -UT1(:,IORB2T,:,:)*SBARLOC(IORB2,IORB1)
        ENDDO
      ENDDO
      DO IORB1=1,LMNX
        DO IORB2=1,LMNX
          IORB2T=POTPAR(ISP)%TAILED%LMNDOT(IORB2)
          UT1(:,:,IORB1,:)=UT1(:,:,IORB1,:) &
     &                    -UT1(:,:,IORB2T,:)*SBARLOC(IORB2,IORB1)
        ENDDO
      ENDDO
      DO IORB1=1,LMNX
        DO IORB2=1,LMNX
          IORB2T=POTPAR(ISP)%TAILED%LMNDOT(IORB2)
          UT1(:,:,:,IORB1)=UT1(:,:,:,IORB1) &
     &                    -UT1(:,:,:,IORB2T)*SBARLOC(IORB2,IORB1)
        ENDDO 
      ENDDO
!
!     ======================================================================
!     ==  EXTRACT RELEVANT PART OF U-TENSOR                               ==
!     ======================================================================
      U(:,:,:,:)=UT1(:LMNX,:LMNX,:LMNX,:LMNX)
!
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_blowupdenmat(iat,ndimd,lmnx,d,lmnxt,dt)
!     **************************************************************************
!     **  core valence exchange energy                                        **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      use lmto_module, only: ispecies,lnx,lox,potpar,sbar,sbarli1,tspherical
      IMPLICIT NONE
      integer(4),intent(in)  :: iat    
      integer(4),intent(in)  :: ndimd
      integer(4),intent(in)  :: lmnx
      integer(4),intent(in)  :: lmnxt
      real(8)   ,intent(in)  :: d(lmnx,lmnx,ndimd)
      real(8)   ,intent(out) :: dt(lmnxt,lmnxt,ndimd)
      integer(4)             :: nns
      real(8)   ,allocatable :: sbarloc(:,:)
      integer(4)             :: i,isp,nn,lmn1,lmn2,lmndot1,lmndot2
      integer(4)             :: L,lmn,ln,i1,im,n1,n2
      real(8)                :: svar
!     **************************************************************************
!
!     ==========================================================================
!     ==  DETERMINE SIZE OF STRUCTURE CONSTANT ARRAY                          ==
!     ==========================================================================
      ISP=ISPECIES(IAT)
      N1=0
      N2=0
      DO LN=1,LNX(ISP)
        L=LOX(LN,ISP)
        N2=N2+2*L+1
        IF(potpar(isp)%LNSCATT(LN).EQ.LN) N1=N1+2*L+1
      ENDDO
      IF(N2.NE.LMNX) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE LMNX')
        CALL ERROR$I4VAL('LMNX',LMNX)
        CALL ERROR$I4VAL('N2',N2)
        CALL ERROR$STOP('LMTO_BLOWUPDT')
      END IF
      ALLOCATE(SBARLOC(N1,N2))
!
!     ==========================================================================
!     ==  COLLECT LOCAL STRUCTURE CONSTANTS                                   ==
!     ==========================================================================
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.IAT) CYCLE
        IF(SBAR(NN)%IAT2.NE.IAT) CYCLE
        IF(SUM(SBAR(NN)%IT(:)**2).NE.0) CYCLE
        IF(SBAR(NN)%N1.NE.n1.OR.SBAR(NN)%N2.NE.N1) THEN
          CALL ERROR$MSG('INCONSISTENT ARRAY SIZES N1,N2')
          CALL ERROR$I4VAL('N1',N1)
          CALL ERROR$I4VAL('SBAR%N1',SBAR(NN)%N1)
          CALL ERROR$I4VAL('SBAR%N2',SBAR(NN)%N2)
          CALL ERROR$STOP('LMTO_BLOWupdt')
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
!     ==========================================================================
!     ==  spherical average of structure constants for test   
!     ==========================================================================
      if(tspherical) Then
        lmn=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I1=SBARLI1(L+1,ISP)-1
          svar=0.d0
          DO IM=1,2*L+1 
            svar=svar+SBARLOC(i1+im,LMN+IM)
          ENDDO
          svar=svar/real(2*l+1,kind=8)
          DO IM=1,2*L+1 
            SBARLOC(:,LMN+IM)=0.d0
            SBARLOC(i1+im,LMN+IM)=svar
          ENDDO
          LMN=LMN+2*L+1
        ENDDO
      end if
!
!     ==========================================================================
!     ==  
!     ==========================================================================
      do i=1,ndimd
        dt(:,:,i)=0.d0
        dt(:lmnx,:lmnx,i)=d(:,:,i)
        dt(lmnx+1:,:lmnx,i)=-matmul(sbarloc,d(:,:,i))
        dt(:lmnx,lmnx+1:,i)=-matmul(d(:,:,i),transpose(sbarloc))
        dt(lmnx+1:,lmnx+1:,i)=-matmul(sbarloc,dt(:lmnx,lmnx+1:,i))
      enddo
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_blowdownht(iat,ndimd,lmnxt,ht,lmnx,h)
!     **************************************************************************
!     **  core valence exchange energy                                        **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      use lmto_module, only: ispecies,lnx,lox,potpar,sbar,sbarli1,tspherical
      IMPLICIT NONE
      integer(4),intent(in)  :: iat    
      integer(4),intent(in)  :: ndimd
      integer(4),intent(in)  :: lmnx
      integer(4),intent(in)  :: lmnxt
      real(8)   ,intent(in)  :: ht(lmnxt,lmnxt,ndimd)
      real(8)   ,intent(out) :: h(lmnx,lmnx,ndimd)
      integer(4)             :: nns
      real(8)   ,allocatable :: sbarloc(:,:)
      real(8)                :: svar
      integer(4)             :: i,isp,nn,lmn1,lmn2,lmndot1,lmndot2,lmn3
      integer(4)             :: L,lmn,ln,i1,im,n1,n2
!     **************************************************************************
!
!     ==========================================================================
!     ==  determine size of structure constant array                          ==
!     ==========================================================================
      isp=ispecies(iat)
      n1=0
      n2=0
      do ln=1,lnx(isp)
        l=lox(ln,isp)
        n2=n2+2*l+1
        if(potpar(isp)%lnscatt(ln).eq.ln) n1=n1+2*l+1
      enddo
      if(n2.ne.lmnx) then
        call error$msg('inconsistent array size lmnx')
        call error$i4val('lmnx',lmnx)
        call error$i4val('n2',n2)
        call error$stop('LMTO_blowdownht')
      end if
      allocate(sbarloc(n1,n2))
!
!     ==========================================================================
!     ==  COLLECT LOCAL STRUCTURE CONSTANTS                                   ==
!     ==========================================================================
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
          CALL ERROR$STOP('LMTO_BLOWDOWNHT')
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
!     ==========================================================================
!     ==  spherical average of structure constants for test   
!     ==========================================================================
      if(tspherical) Then
        lmn=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          I1=SBARLI1(L+1,ISP)-1
          svar=0.d0
          DO IM=1,2*L+1 
            svar=svar+SBARLOC(i1+im,LMN+IM)
          ENDDO
          svar=svar/real(2*l+1,kind=8)
          DO IM=1,2*L+1 
            SBARLOC(:,LMN+IM)=0.d0
            SBARLOC(i1+im,LMN+IM)=svar
          ENDDO
          LMN=LMN+2*L+1
        ENDDO
      end if
!
!     ==========================================================================
!     ==  
!     ==========================================================================
      do i=1,ndimd
        h(:,:,i)=ht(:lmnx,:lmnx,i) &
     &          -matmul(ht(:lmnx,lmnx+1:,i),sbarloc) &
     &          -matmul(transpose(sbarloc),ht(lmnx+1:,:lmnx,i)) &
     &          +matmul(transpose(sbarloc),matmul(ht(lmnx+1:,lmnx+1:,i),sbarloc))
      enddo
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SIMPLEENERGYTEST()
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ISPECIES,DENMAT,HAMIL,lnx,lox
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
      INTEGER(4)            :: lnx1,lmrx,lrx
      INTEGER(4)            :: gid
      INTEGER(4)            :: nr
      logical(4),allocatable:: torb(:)
      REAL(8)   ,ALLOCATABLE:: chi(:,:)
      REAL(8)   ,ALLOCATABLE:: aecore(:)
      REAL(8)   ,ALLOCATABLE:: chiphi(:,:)
      REAL(8)   ,ALLOCATABLE:: ulittle(:,:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: U(:,:,:,:)
      REAL(8)   ,ALLOCATABLE:: D(:,:,:)
      REAL(8)   ,ALLOCATABLE:: H(:,:,:)
      REAL(8)               :: EH,EX,EXTOT,EHTOT,EHDC,EXDC
      INTEGER(4)            :: NN,IAT,I,J,K,L,IS,IAT1,IAT2,isp
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      REAL(8)   ,PARAMETER  :: HFWEIGHT=0.25D0
character(128) :: string
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_ENERGYTEST')
print*,'========================= energytest ==============================='
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
        CALL RADIAL$GETi4(GID,'NR',NR)
        ALLOCATE(AECORE(NR))
        CALL SETUP$GETR8A('AECORE',nr,AECORE)
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
        CALL LDAPLUSU_ULITTLE(GID,NR,LRX,LNX1,LOX(:LNX1,ISP),CHI,ULITTLE)
        CALL LDAPLUSU_UTENSOR(LRX,NORB,LNX1,LOX(:LNX1,ISP),ULITTLE,U)
        DEALLOCATE(ULITTLE)
        DEALLOCATE(CHIPHI)
        DEALLOCATE(torb)
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
print*,'exact exchange energy for atom=',iat,ex
!
!       ========================================================================
!       == ADD CORE VALENCE EXCHANGE                                          ==
!       ========================================================================
        CALL LMTO_CVX(ISPECIES(IAT),NORB,EX,D(:,:,1),H(:,:,1))
        EXTOT=EXTOT+EX
        HAMIL(INH)%MAT(:,:,1)=HAMIL(INH)%MAT(:,:,1)+H(:,:,1)
print*,'core valence exchange energy for atom=',iat,ex
!
!       ========================================================================
!       == DOUBLE COUNTING CORRECTION (EXCHANGE ONLY)                         ==
!       ========================================================================
! this is the time consumin part of energytest
CALL TIMING$CLOCKON('ENERGYTEST:DC')      
!        CALL LMTO_DOUBLECOUNTING(IAT,NDIMD,NORB,D,EX,H)
        call lmto_simpledc(GID,NR,norb,LNX1,LOX(:lnx1,isp),CHI,LRX,AECORE &
     &                        ,D,Ex,H)
        EXTOT=EXTOT-EX
        HAMIL(INH)%MAT=HAMIL(INH)%MAT-H
print*,'double counting ',iat,ex
CALL TIMING$CLOCKOFF('ENERGYTEST:DC')      
        DEALLOCATE(U)
        DEALLOCATE(H)
        DEALLOCATE(D)
        deallocate(chi)
        deallocate(aecore)
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
        WRITE(*,FMT='(82("="),T10," ONLY ONSITE TERMS written ")')
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
      SUBROUTINE lmto_simpledc(GID,NR,LMNX,LNX,LOX,CHI,LRX,AECORE &
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
      real(8)     ,INTENT(IN) :: DENMAT(LMNX,LMNX,4) ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT       ! DOUBLE COUNTINNG ENERGY
      real(8)     ,INTENT(OUT):: HAM(LMNX,LMNX,4)  ! DETOT/D(RHO^*)        
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
      denmat1=cmplx(denmat)
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
      HAM=real(ham1)
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
! this is the time consumin part of energytest
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
        WRITE(*,FMT='(82("="),T10," ONLY ONSITE TERMS written ")')
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
!     **  core valence exchange energy                                        **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      use lmto_module, only: potpar
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: ISP ! ATOM-TYPE INDEX
      INTEGER(4),INTENT(IN)  :: norb ! #(local orbitals)
      REAL(8)   ,INTENT(IN)  :: D(NORB,NORB)  ! TOTAL DENSITY MATRIX
      REAL(8)   ,INTENT(OUT) :: EX            ! CORE-VALENCE ENERGY
      REAL(8)   ,INTENT(OUT) :: H(NORB,NORB)  ! hamiltonian contribution
      integer(4)             :: lnx
      integer(4),allocatable :: lox(:)
      real(8)   ,allocatable :: cvxmat(:,:)
      real(8)                :: c1,c2,svar
      integer(4)             :: ln1,ln2,l1,l2,lmn1,lmn2,im
!     **************************************************************************
!
!     ==========================================================================
!     == collect data                                                         ==
!     ==========================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      ALLOCATE(CVXMAT(LNX,LNX))
      CALL SETUP$GETi4A('LOX',LNX,LOX)
      CALL SETUP$GETR8A('CVX',LNX*LNX,CVXMAT)
!
!     ==========================================================================
!     == calculate core-valence exchange energy                               ==
!     ==========================================================================
      EX=0.D0
      H(:,:)=0.D0
      lmn1=0
      do ln1=1,lnx
        l1=lox(ln1)
        c1=potpar(isp)%ktophi(ln1)
        lmn2=0
        do ln2=1,lnx
          l2=lox(ln2)
          if(l2.eq.l1) then
            c2=potpar(isp)%ktophi(ln2)
            svar=c1*cvxmat(ln1,ln2)*c2
            do im=1,2*l1+1
              ex=ex+d(lmn2+im,lmn1+im)*svar
              h(lmn2+im,lmn1+im)=h(lmn2+im,lmn1+im)+svar
            enddo
          end if
          lmn2=lmn2+2*l2+1
        enddo
        lmn1=lmn1+2*l1+1
      enddo
!
!     ==========================================================================
!     == close down                                                           ==
!     ==========================================================================
      DEALLOCATE(CVXMAT)
      DEALLOCATE(LOX)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DOUBLECOUNTING(iat,ndimd,norb,d,Ex,h)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : ISPECIES,LMORB
      IMPLICIT NONE
      integer(4),intent(in) :: iat
      integer(4),intent(in) :: norb
      integer(4),intent(in) :: ndimd
      real(8)   ,intent(in) :: d(norb,norb,ndimd)
      REAL(8)   ,INTENT(OUT):: Ex
      real(8)   ,intent(out):: h(norb,norb,ndimd)
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
      INTEGER(4)            :: I,J,ISP,l
      INTEGER(4)            :: IORB1,IORB2,LM1,LM2,LM3,IDIMD
      LOGICAL(4),PARAMETER  :: TPR=.false.
      LOGICAL(4),PARAMETER  :: TCV=.true. ! CORE-VALENCE CONTRIBUTION
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
!     == calculate density                                                    ==
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
        CALL SETUP$GETR8A('AECORE',NR,aux)
        CALL AUGMENTATION_XC(GID,NR,1,1,aux,ETOTC,POT)
        RHOWC(:,1,1)=RHO(:,1,1)+aux(:)
      else 
        etotc=0.d0
      END IF
!
!     ========================================================================
!     == calculate total energy                                             ==
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
      EX=EX-ETOTC  ! subtract core exchange
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
     &                         ,' EXC=',EX+EH,' ex(core)=',etotc
        DO iorb1=1,norb 
          WRITE(*,FMT='(I3,30F10.3)')Iorb1,H(iorb1,:,1)
        ENDDO
        WRITE(*,FMT='(82("-"))')
      END IF

                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_DOUBLECOUNTING_old(EHTOT,EXTOT)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : ISPECIES,DENMAT,HAMIL,LMORB
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT):: EHTOT
      REAL(8)   ,INTENT(OUT):: EXTOT
      INTEGER(4)            :: NND
      INTEGER(4)            :: NNH
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IND,INH
      INTEGER(4)            :: NORB
      INTEGER(4)            :: IT(3)
      INTEGER(4)            :: N1,N2,N3
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: LMX
      INTEGER(4)            :: LMRX
      INTEGER(4)            :: NDIMD
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)   ,ALLOCATABLE:: AECORE(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:)
      REAL(8)   ,ALLOCATABLE:: AUX1(:)
      REAL(8)   ,ALLOCATABLE:: D(:,:,:)
      REAL(8)   ,ALLOCATABLE:: H(:,:,:)
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)  !(NR,LMR,IDIMD) DENSITY
      REAL(8)   ,ALLOCATABLE:: RHOWC(:,:,:)  !(NR,LMR,IDIMD) DENSITY W CORE
      REAL(8)   ,ALLOCATABLE:: POT(:,:,:)  !(NR,LMR,IDIMD) POTENTIAL
      REAL(8)               :: EX,EH,ETOTC,ETOTV
      REAL(8)               :: SVAR
      REAL(8)               :: CG ! GAUNT COEFFICIENT
      INTEGER(4)            :: NN,IAT,I,J,K,L,IS,IAT1,IAT2,ISP
      INTEGER(4)            :: IORB1,IORB2,LM1,LM2,LM3,IDIMD
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      LOGICAL(4),PARAMETER  :: TCV=.true. ! CORE-VALENCE CONTRIBUTION
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO_DOUBLECOUNTING')
      NAT=SIZE(ISPECIES)
      NND=SIZE(DENMAT)
      NNH=NND
      IF(.NOT.ALLOCATED(HAMIL)) THEN
        CALL ERROR$MSG('ARRAY HAMIL NOT ASSOCIATED')
        CALL ERROR$STOP('LMTO_DOUBLECOUNTING')
      END IF

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
          CALL ERROR$STOP('LMTO_DOUBLECOUNTING')
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
          CALL ERROR$STOP('LMTO_DOUBLECOUNTING')
        END IF
!
!       ========================================================================
!       == COLLECT FURTHER INFORMATION                                        ==
!       ========================================================================
        GID=LMORB(IAT)%GID
        NR=LMORB(IAT)%NR
        LMX=LMORB(IAT)%LMX
        NORB=LMORB(IAT)%NORB
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
!
        ISP=ISPECIES(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LMRX',LMRX)
        ALLOCATE(AECORE(NR))
        IF(TCV) THEN
          CALL SETUP$GETR8A('AECORE',NR,AECORE)
        ELSE
          AECORE(:)=0.D0
        END IF
        ALLOCATE(AUX(NR))
        ALLOCATE(AUX1(NR))
!
        NDIMD=DENMAT(IND)%N3
        IF(NDIMD.NE.4) THEN
          CALL ERROR$MSG('NDIMD DIFFERS FROM 4')
          CALL ERROR$STOP('LMTO_DOUBLECOUNTING')
        END IF
!
!       ========================================================================
!       == WORK OUT TOTAL ENERGY AND HAMILTONIAN                              ==
!       ========================================================================
!       == CALCULATE LOCAL ELECTRON DENSITY ====================================
        ALLOCATE(D(NORB,NORB,NDIMD))
        D=DENMAT(IND)%MAT
        ALLOCATE(RHO(NR,LMRX,NDIMD))
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
        DEALLOCATE(D)
!
!       == ADD CORE DENSITY IF NECESSARY =======================================
        ALLOCATE(RHOWC(NR,LMRX,NDIMD))  !WITH CORE
        RHOWC=RHO
        RHOWC(:,1,1)=RHO(:,1,1)+AECORE(:)
!
!== CALCULATE HARTREE ENERGY FOR TEST ==========================================
AUX1=0.D0
DO LM1=1,LMRX
  L=INT(SQRT(REAL(LM1-1,KIND=8))+1.D-5)
  CALL RADIAL$POISSON(GID,NR,L,RHO(:,LM1,1),AUX)
  AUX1(:)=AUX1(:)+0.5D0*AUX(:)*RHO(:,LM1,1)
ENDDO
CALL RADIAL$INTEGRAL(GID,NR,AUX1*R**2,EH)
PRINT*,'dc EH ',EH
!
!       == EXCHANGE CORRELATION ENERGY AND POTENTIAL ===========================
        ALLOCATE(POT(NR,LMRX,NDIMD))
        POT(:,:,:)=0.D0
        CALL DFT$SETL4('XCONLY',.TRUE.)
        IF(TCV) THEN
          CALL AUGMENTATION_XC(GID,NR,1,1,AECORE,ETOTC,POT)
        ELSE
          ETOTC=0.D0
        END IF
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOWC,EX,POT)
        CALL DFT$SETL4('XCONLY',.FALSE.)
        EX=EX-ETOTC
PRINT*,'dc EX ',EX,'ETOTV ',ETOTV,'ETOTC ',ETOTC,' EXTOT ',EX+ETOTC
        DEALLOCATE(RHO)
        DEALLOCATE(RHOWC)
        DEALLOCATE(AECORE)
!
!       == WORK OUT HAMILTON MATRIX ELEMENTS ===================================
        ALLOCATE(H(NORB,NORB,4))
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
        HAMIL(INH)%MAT=HAMIL(INH)%MAT-H
        DEALLOCATE(H)
        DEALLOCATE(POT)

        EXTOT=EXTOT+EX
        EHTOT=EHTOT+EH
        DEALLOCATE(R)
        DEALLOCATE(AUX)
        DEALLOCATE(AUX1)
      ENDDO
!
!     ==========================================================================
!     == PRINT FOR TEST                                                       ==
!     ==========================================================================
      IF(TPR) THEN
        PRINT*,'LMTO_DOUBLECOUNTING TOTAL ENERGY=',EXTOT+EHTOT,' HARTREE=',EHTOT,' EXC=',EXTOT
        WRITE(*,FMT='(82("="),T10," NON-LOCAL HAMILTON MATRIX IN A NTBO BASIS ")')
        DO NN=1,NND
          IAT1=HAMIL(NN)%IAT1
          IAT2=HAMIL(NN)%IAT2
          IT=HAMIL(NN)%IT
          WRITE(*,FMT='(82("="),T10," IAT1 ",I4," IAT2=",I4," IT=",3I3)') &
     &                                                               IAT1,IAT2,IT
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

                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TESTDENMAT_1cdenmat(lmnxx_,ndimd_,nat,denmat_)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      USE LMTO_MODULE, ONLY : DENMAT,sbar,ispecies,potpar,lnx,lox
      IMPLICIT NONE
      integer(4),intent(in) :: lmnxx_
      integer(4),intent(in) :: ndimd_
      integer(4),intent(in) :: nat
      complex(8),intent(in) :: denmat_(lmnxx_,lmnxx_,ndimd_,nat)
      INTEGER(4)            :: lmnxx
      INTEGER(4)            :: NND
      INTEGER(4)            :: NNs
      INTEGER(4)            :: IAT1,IAT2,IT(3)
      INTEGER(4)            :: N1,N2,n3
      INTEGER(4)            :: ind,ins
      INTEGER(4)            :: NN,MM,I,J,idim,iat,isp,ln,lmn,im,i1,i2,i3
      REAL(8)               :: SVAR1,SVAR2
      real(8)   ,allocatable:: ktophi(:)
      real(8)   ,allocatable:: ktophidot(:)
      real(8)   ,allocatable:: jbartophidot(:)
      real(8)   ,allocatable:: mat11(:,:,:)
      real(8)   ,allocatable:: mat12(:,:,:)
      real(8)   ,allocatable:: mat22(:,:,:)
      logical(4),parameter  :: tpr=.false.
!     **************************************************************************
      if(.not.tpr) return
                                  CALL TRACE$PUSH('LMTO_TESTDENMAT_1cdenmat')
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      if(nat.ne.SIZE(ISPECIES)) then
        call error$msg('inconsistent data: nat differs from nat_')
        call error$stop('LMTO_TESTDENMAT_1cdenmat')
      end if
!
!     ==========================================================================
!     == prepare wronskians                                                   ==
!     ==========================================================================
!
!     ==========================================================================
!     == calculate one-center density matrix                                  ==
!     ==========================================================================
      NND=SIZE(DENMAT)
      NNs=SIZE(sbar)
      do iat=1,nat
        isp=ispecies(iat)
!       == expand potential parameter arrays ===================================
        lmnxx=sum(2*lox(:,isp)+1)
        allocate(ktophi(lmnxx))
        allocate(ktophidot(lmnxx))
        allocate(jbartophidot(lmnxx))
        allocate(mat11(lmnxx,lmnxx,4))
        allocate(mat12(lmnxx,lmnxx,4))
        allocate(mat22(lmnxx,lmnxx,4))
        lmn=0
        do ln=1,lnx(isp)
          do im=1,2*lox(ln,isp)+1
            lmn=lmn+1
            ktophi(lmn)      =potpar(isp)%ktophi(ln)
            ktophidot(lmn)   =potpar(isp)%ktophidot(ln)
            jbartophidot(lmn)=potpar(isp)%jbartophidot(ln)
          enddo
        enddo
!       == write original density matrix =======================================
        do idim=1,ndimd_
          write(*,fmt='(82("="),t30," iat=",i3,"  and idim=",i1,"  ")') &
    &                                         iat,idim
          do i=1,lmnxx_
            write(*,fmt='(100f10.5)')real(denmat_(i,:,idim,iat))
          enddo
        enddo
!
!       == calculate from ntbo density matrix ==================================
        do nn=1,nns
          if(denmat(nn)%iat1.eq.iat) cycle
          if(denmat(nn)%iat2.eq.iat) cycle
          if(sum(denmat(nn)%it**2).ne.0) cycle
          n1=denmat(nn)%n1  
          n2=denmat(nn)%n2  
          n3=denmat(nn)%n3  
          mat11(:,:,:)=denmat(nn)%mat(:,:,:)
          mat12(:,:,:)=mat11(:,:,:)
          mat22(:,:,:)=mat11(:,:,:)
          do idim=1,n3
            do i2=1,n2
              mat11(:,i2,idim)=ktophi(:)*mat11(:,i2,idim)
              mat12(:,i2,idim)=ktophi(:)*mat12(:,i2,idim)
              mat12(:,22,idim)=ktophidot(:)*mat22(:,i2,idim)
            enddo
            do i1=1,n1
              mat11(i1,:,idim)=mat11(i1,:,idim)*ktophi(:)
              mat12(i1,:,idim)=mat12(i1,:,idim)*ktophidot(:)
              mat22(i1,:,idim)=mat22(i1,:,idim)*ktophidot(:)
            enddo
          enddo
        enddo
!
        do ind=1,nnd
          if(denmat(ind)%iat1.eq.iat) cycle
          do ins=1,nns
            if(sbar(ins)%iat2.eq.iat) cycle
            if(sum(denmat(ind)%it+sbar(ins)%it)**2.ne.0) cycle
          enddo
        enddo


!       == write result ========================================================
        do idim=1,n3
          write(*,fmt='(82("-"),t30," iat=",i3,"  and idim=",i1,"  ")') &
    &                                       iat,idim
          do i=1,n1
            write(*,fmt='(100f10.5)')mat11(i,:,idim)
          enddo
        enddo
!
        deallocate(ktophi)
        deallocate(ktophidot)
        deallocate(jbartophidot)
        deallocate(mat11)
        deallocate(mat12)
        deallocate(mat22)
      enddo
stop 'forced stop in lmto_testdenmat_1center'
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
      logical(4),parameter  :: tpr=.false.
!     **************************************************************************
      if(.not.tpr) return
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
      INTEGER(4)             :: NPRO
      INTEGER(4)             :: NB,NBH,NBX
      INTEGER(4)             :: NDIMD
      REAL(8)   ,ALLOCATABLE :: XK(:,:)
      INTEGER(4),ALLOCATABLE :: IPRO1(:)
      INTEGER(4),ALLOCATABLE :: NPROAT(:)
      REAL(8)   ,ALLOCATABLE :: OCC(:,:,:)
      INTEGER(4)             :: IAT,NN,ISP,IPRO,IKPT,ISPIN,I,J,IBH,IB
      REAL(8)                :: SVAR
      REAL(8)                :: F1,F2
      LOGICAL(4)             :: TINV
      INTEGER(4)             :: IAT1,IAT2,IT(3),I0,J0,IDIM,JDIM
      COMPLEX(8)             :: EIKR,C1(ndim),C2(ndim),CSVAR,csvar22(ndim,ndim)
      COMPLEX(8),PARAMETER   :: CI=(0.D0,1.D0)
      REAL(8)                :: PI
      logical(4),parameter   :: tpr=.false.
!     **************************************************************************
                                           CALL TRACE$PUSH('LMTO_NTBODENMAT')
      PI=4.D0*ATAN(1.D0)
      if(.not.associated(this%tbc)) then
        call error$msg('this%tbc not associated')
        call error$stop('LMTO_NTBODENMAT')
      end if
!
!     ==========================================================================
!     ==========================================================================
      NNS=SIZE(SBAR)
      IF(.NOT.ALLOCATED(DENMAT))ALLOCATE(DENMAT(NNS))
      DO NN=1,NNS
        IAT1=SBAR(NN)%IAT1
        IAT2=SBAR(NN)%IAT2
        ISP=ISPECIES(IAT1)
        N1=SUM(2*LOX(:LNX(ISP),ISP)+1)
        ISP=ISPECIES(IAT2)
        N2=SUM(2*LOX(:LNX(ISP),ISP)+1)
        DENMAT(NN)%IAT1=IAT1
        DENMAT(NN)%IAT2=IAT2
        DENMAT(NN)%IT=SBAR(NN)%IT
        DENMAT(NN)%N1=N1
        DENMAT(NN)%N2=N2
        DENMAT(NN)%N3=4  !(total,x,y,z)
        ALLOCATE(DENMAT(NN)%MAT(N1,N2,4))
        DENMAT(NN)%MAT(:,:,:)=0.D0
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
        NPROAT(IAT)=SUM(2*LOX(:LNX(ISP),ISP)+1)
        IPRO=IPRO+NPROAT(IAT)
      ENDDO
!
!     ==========================================================================
!     ==  ADD UP DENSITY MATRIX                                               ==
!     ==========================================================================
      NPRO=MAP%NPRO
      DO IKPT=1,NKPTL
        DO ISPIN=1,NSPIN
          CALL WAVES_SELECTWV(IKPT,ISPIN)
          NBH=THIS%NBH
          NB=THIS%NB
          CALL PLANEWAVE$SELECT(GSET%ID)
          CALL PLANEWAVE$GETL4('TINV',TINV)
          DO NN=1,NNS
            IAT1=DENMAT(NN)%IAT1
            IAT2=DENMAT(NN)%IAT2
            IT=DENMAT(NN)%IT
            SVAR=2.D0*PI*SUM(XK(:,IKPT)*REAL(IT,kind=8))
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
                    C2(:)=THIS%TBC(:,IBH,J0+J)*EIKR
                    do jdim=1,ndim
                      CSVAR22(:,jdim)=CSVAR22(:,jdim) &
        &                            +0.5D0*((F1+F2)*C1(:)*CONJG(C2(jdim)) &
        &                                   +(F1-F2)*C1(:)*C2(jdim))
                    enddo
                  ENDDO
                  CSVAR22=REAL(CSVAR22) ! IMAG(CSVAR) CONTAINS CRAP DUE TO SUPER WAVE FUNCTIONS
                ELSE
                  CSVAR=(0.D0,0.D0)
                  DO IB=1,NB
                    F1=OCC(IB,IKPT,ISPIN)
                    C1(:)=THIS%TBC(:,IB,I0+I)
                    C2(:)=THIS%TBC(:,IB,J0+J)*EIKR
                    do jdim=1,ndim
                      CSVAR22(:,jdim)=CSVAR22(:,jdim)+F1*C1(:)*CONJG(C2(jdim))
                    enddo
                  ENDDO
                END IF
!
!           == DISTRIBUTE ONTO DENSITY MATRIX ENTRIES ==========================
!           == D(IDIMD)=SUM_{IDIM,JDIM} D(IDIM,JDIM)*PAULI_{IDIMD}(JDIM,IDIM) ==
!           == idimd in {total,x,y,z}; idim in {up,down} =======================
!           == TRANSFORMATION MUST BE CONSISTENT WITH WAVES_DENMAT =============
                IF(NSPIN.EQ.1) THEN
                  IF(NDIM.EQ.1) THEN !NON-SPIN-POLARIZED
                    DENMAT(NN)%MAT(I,J,1)=DENMAT(NN)%MAT(I,J,1) &
          &                              +REAL(CSVAR22(1,1))
                  ELSE ! NONCOLLINEAR
                    DENMAT(NN)%MAT(I,J,1)=DENMAT(NN)%MAT(I,J,1) &
          &                              +REAL(CSVAR22(1,1)+csvar22(2,2))
                    DENMAT(NN)%MAT(I,J,2)=DENMAT(NN)%MAT(I,J,2) &
          &                              +REAL(CSVAR22(1,2)+csvar22(2,1))
                    DENMAT(NN)%MAT(I,J,3)=DENMAT(NN)%MAT(I,J,3) &
          &                              -aimag(CSVAR22(1,2)-csvar22(2,1))
                    DENMAT(NN)%MAT(I,J,4)=DENMAT(NN)%MAT(I,J,4) &
          &                              +REAL(CSVAR22(1,1)-csvar22(2,2))
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
      if(tpr) then
        WRITE(*,FMT='(82("="),T10," NON-LOCAL DENSITY MATRIX IN A NTBO BASIS ")')
        DO NN=1,NNS
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
!!$        CALL ERROR$MSG('FORCED STOP AFTER PRINTING DENSITY MATRIX')
!!$        CALL ERROR$STOP('LMTO_NTBODENMAT')
      end if
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
      USE WAVES_MODULE, ONLY: wvset_type,NKPTL,NSPIN,NDIM,THIS,MAP &
     &                       ,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : SBAR,HAMIL,LOX,LNX,ISPECIES
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      INTEGER(4)             :: N1,N2
      INTEGER(4)             :: NNS
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
      COMPLEX(8)             :: EIKR,C1,C2,CSVAR,csvar22(ndim,ndim)
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
!
!     ==========================================================================
!     ==  ADD UP DENSITY MATRIX                                               ==
!     ==========================================================================
      NNS=SIZE(HAMIL)
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
          DO NN=1,NNS
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
                csvar22(:,:)=csvar22(:,:)*CONJG(EIKR)
!
                DO IB=1,NBH
                  DO IDIM=1,NDIM
                    DO jDIM=1,NDIM
                      THIS%HTBC(jdim,IB,J0+J)=THIS%HTBC(jdim,IB,J0+J) &
      &                              +CSVAR22(IDIM,jdim)*THIS%TBC(IDIM,IB,I0+I)
                    enddo
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
      SUBROUTINE LMTO_UTENSORLAYOUT()
!     **************************************************************************
!     **  SELECTS THE ENTRIES FOR THE  THE U-TENSOR                           **
!     **                                                                      **
!     **  PRELIMINARY VERSION: ONLY ONSITE MATRIX ELEMENTS                    **
!     **  PRELIMINARY VERSION: NO AUGMENTATION CONTRIBUTIONS INCLUDED         **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : UTENSOR,ISPECIES,LNX,LOX
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      INTEGER(4)             :: N1
      INTEGER(4)             :: IAT,NN,ISP
!     **************************************************************************
      IF(ALLOCATED(UTENSOR)) RETURN
                                   CALL TRACE$PUSH('LMTO_UTENSORLAYOUT')
!
!     ==========================================================================
!     == ALLOCATE UTENSOR AND CHOSE MATRIX ELEMENTS TO BE EVALUATED           ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(UTENSOR(NAT))
      NN=0
      DO IAT=1,NAT
        NN=NN+1
        ISP=ISPECIES(IAT)
        N1=SUM(2*LOX(:LNX(ISP),ISP)+1)
        UTENSOR(NN)%IAT1=IAT
        UTENSOR(NN)%IAT2=IAT
        UTENSOR(NN)%IAT3=IAT
        UTENSOR(NN)%IAT4=IAT
        UTENSOR(NN)%IT2=(/0,0,0/)
        UTENSOR(NN)%IT3=(/0,0,0/)
        UTENSOR(NN)%IT4=(/0,0,0/)
        UTENSOR(NN)%N1=N1
        UTENSOR(NN)%N2=N1
        UTENSOR(NN)%N3=N1
        UTENSOR(NN)%N4=N1
        ALLOCATE(UTENSOR(NN)%U(N1,N1,N1,N1))
      ENDDO
                                   CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_ONSITEU()
!     **************************************************************************
!     **                                                                      **
!     ** W=0.5 SUM_{IJKL} W_IJKL CDAGGER_I CDAGGER_J C_L C_K                  **
!     ** WIJKL=INT DR INT DR' CHI_I^*(R)CHI^*_J(R')CHI_K(R)CHI_L(R')/|R-R'|   **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY: LMORB,ISPECIES,UTENSOR
      INTEGER(4)            :: NAT
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: LMX
      INTEGER(4)            :: NORB
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:)
      INTEGER(4)            :: NPAIR
      REAL(8)   ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)   ,ALLOCATABLE:: POT(:,:,:)
      INTEGER(4)            :: IORB1,IORB2,IORB3,IORB4
      INTEGER(4)            :: LM,LM1,LM2,LM3,IPAIR,IPAIR1,IPAIR2,NN
      INTEGER(4)            :: L
      INTEGER(4)            :: ISP
      INTEGER(4)            :: IAT
      INTEGER(4)            :: LMRX
      INTEGER(4)            :: IT(3)
      INTEGER(4)            :: NNU  ! SIZE OF UTENSOR ARRAY
      REAL(8)               :: CG   ! GAUNT COEFFICIENT
      REAL(8)   ,ALLOCATABLE:: U(:,:)
      CHARACTER(64)         :: FILE
      LOGICAL(4),PARAMETER  :: TPR=.false.
!     **************************************************************************
      NNU=SIZE(UTENSOR)
      DO NN=1,NNU
        IAT=UTENSOR(NN)%IAT1
        IF(UTENSOR(NN)%IAT2.NE.IAT) CYCLE
        IF(UTENSOR(NN)%IAT3.NE.IAT) CYCLE
        IF(UTENSOR(NN)%IAT4.NE.IAT) CYCLE
        IT=UTENSOR(NN)%IT2
        IF(IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
        IT=UTENSOR(NN)%IT3
        IF(IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
        IT=UTENSOR(NN)%IT4
        IF(IT(1).NE.0.OR.IT(2).NE.0.OR.IT(3).NE.0) CYCLE
!
        ISP=ISPECIES(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LMRX',LMRX)
        GID=LMORB(IAT)%GID
        NR=LMORB(IAT)%NR
        LMX=LMORB(IAT)%LMX
        NORB=LMORB(IAT)%NORB
        ALLOCATE(R(NR))
        ALLOCATE(AUX(NR))
        CALL RADIAL$R(GID,NR,R)

!!$DO IORB1=1,NORB
!!$  WRITE(FILE,*)IORB1
!!$  FILE='ORB'//TRIM(ADJUSTL(FILE))//'.DAT'
!!$  CALL SETUP_WRITEPHI(FILE,GID,NR,LMRX,LMORB(IAT)%F(:,:LMRX,IORB1))
!!$ENDDO
!
!       ========================================================================
!       ==  CONSTRUCT DENSITIES                                               ==
!       ========================================================================
        NPAIR=NORB*(NORB+1)/2
        ALLOCATE(RHO(NR,LMRX,NPAIR))
        RHO(:,:,:)=0.D0
        IPAIR=0
        DO IORB1=1,NORB
          DO IORB2=IORB1,NORB
            IPAIR=IPAIR+1
            DO LM1=1,LMX
              DO LM2=1,LMX
                AUX(:)=LMORB(IAT)%F(:,LM1,IORB1)*LMORB(IAT)%F(:,LM2,IORB2)
                DO LM3=1,LMRX
                  CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
                  RHO(:,LM3,IPAIR)=RHO(:,LM3,IPAIR)+AUX(:)*CG
                ENDDO
              ENDDO
            ENDDO
            IF(IORB1.EQ.IORB2) THEN
              WRITE(FILE,*)IORB1
              FILE='RHO_'//TRIM(ADJUSTL(FILE))//'.DAT'
              CALL SETUP_WRITEPHI(FILE,GID,NR,LMRX,RHO(:,:,IPAIR))
            END IF
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  CONSTRUCT POTENTIALS                                              ==
!       ========================================================================
        ALLOCATE(POT(NR,LMRX,NPAIR))
        DO IPAIR=1,NPAIR
          DO LM=1,LMRX
            L=INT(SQRT(LM-1+1.D-10))
            CALL RADIAL$POISSON(GID,NR,L,RHO(:,LM,IPAIR),POT(:,LM,IPAIR))
          ENDDO
        ENDDO
!
!       ========================================================================
!       ==  ONSITE U-TENSOR                                                   ==
!       ========================================================================
IF(TPR) THEN
  WRITE(*,FMT='(82("="),T20," ONSITE U-TENSOR  FOR ATOM ",I5,"  ")')iat
END IF
        ALLOCATE(U(NPAIR,NPAIR))
        U(:,:)=0.D0
        DO IPAIR1=1,NPAIR
          DO IPAIR2=IPAIR1,NPAIR
            AUX(:)=0.D0
            DO LM=1,LMRX
              AUX(:)=AUX(:)+RHO(:,LM,IPAIR1)*POT(:,LM,IPAIR2)
            ENDDO
            AUX(:)=AUX(:)*R(:)**2
            CALL RADIAL$INTEGRAL(GID,NR,AUX,U(IPAIR1,IPAIR2))
            U(IPAIR2,IPAIR1)=U(IPAIR1,IPAIR2)
          ENDDO
        ENDDO
!
!       == WRITE ===============================================================
        IPAIR1=0
        DO IORB1=1,NORB
          DO IORB2=IORB1,NORB
            IPAIR1=IPAIR1+1
            IPAIR2=0
            DO IORB3=1,NORB
              DO IORB4=IORB3,NORB
                IPAIR2=IPAIR2+1
                UTENSOR(NN)%U(IORB1,IORB3,IORB2,IORB4)=U(IPAIR1,IPAIR2) 
                UTENSOR(NN)%U(IORB2,IORB3,IORB1,IORB4)=U(IPAIR1,IPAIR2) 
                UTENSOR(NN)%U(IORB1,IORB4,IORB2,IORB3)=U(IPAIR1,IPAIR2) 
                UTENSOR(NN)%U(IORB2,IORB4,IORB1,IORB3)=U(IPAIR1,IPAIR2) 
!
IF(TPR) THEN
  IF(ABS(U(IPAIR1,IPAIR2)).LT.1.D-5) CYCLE
  WRITE(*,*)IORB1,IORB2,IORB3,IORB4,U(IPAIR1,IPAIR2)
END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(AUX)
        DEALLOCATE(R)
        DEALLOCATE(RHO)
        DEALLOCATE(POT)
        DEALLOCATE(U)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GAUSSORBTOLM()
!     **************************************************************************
!     ** MAP ORBITALS FROM GAUSS EXPANSION ONTO RADIAL GRIDS TIMES            **
!     ** SPHERICAL HARMONICS                                                  **
!     **************************************************************************
      USE LMTO_MODULE, ONLY: GAUSSORBAUG,LMORB,ISPECIES
      INTEGER(4)             :: NIJK,LX,INDX,LM,IND,IAT,I,J,K,L,IM
      INTEGER(4)             :: NE,NORB,IE,N,IORB,ISP
      REAL(8)   ,ALLOCATABLE :: YLMFROMPOL(:,:,:)
      INTEGER(4)             :: GID
      REAL(8)                :: R1
      REAL(8)                :: DEX
      INTEGER(4)             :: LMX1
      INTEGER(4)             :: NAT
      INTEGER(4)             :: NR
      INTEGER(4)             :: LMX
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: FEXPR2N(:)
      REAL(8)   ,ALLOCATABLE :: VEC(:)
      REAL(8)   ,ALLOCATABLE :: YLM(:)
      LOGICAL(4),PARAMETER   :: TWRITE=.false.
character(128) :: string,string1
!     **************************************************************************
!
!     ==========================================================================
!     == DETERMINE MATRIX FOR THE TRANSFORMATION TO SPHERICAL HARMONICS       ==
!     ==  IND->(I,J,K) ;                                                      ==
!     ==  (X/R)^I * (Y/R)^J * (Z/R)^K                                         ==
!     ==         =SUM YLM(R) * |R|^(L+2N) * YLMFROMPOL(N+1,LM,IND)            ==
!     ==========================================================================
      NIJK=MAXVAL(GAUSSORBAUG(:)%NIJK)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJK,I,J,K)
      LX=MAX(I,J,K)
      INDX=(LX+1)*(LX+2)*(LX+3)/6
      ALLOCATE(YLMFROMPOL(LX/2+1,(LX+1)**2,INDX))
      CALL GAUSSIAN_YLMFROMPOL(LX,YLMFROMPOL)
!
!     ==========================================================================
!     ==========================================================================
      IF(TWRITE) THEN
        DO IND=1,INDX
          PRINT*,'================== IND=',IND,' =============================='
          DO LM=1,(LX+1)**2
            IF(MAXVAL(ABS(YLMFROMPOL(:,LM,IND))).EQ.0.D0) CYCLE
            WRITE(*,FMT='("LM=",I3," C(:,LM)=",20F10.5)')LM,YLMFROMPOL(:,LM,IND)
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     == DEFINE RADIAL GRID                                                   ==
!     ==========================================================================
!!$      CALL RADIAL$NEW('SHLOG',GID)
!!$      CALL RADIAL$GRIDPARAMETERS(0.01D0,0.05D0,10.D0,R1,DEX,NR)
!!$      CALL RADIAL$SETR8(GID,'R1',R1)
!!$      CALL RADIAL$SETR8(GID,'DEX',DEX)
!!$      CALL RADIAL$SETI4(GID,'NR',NR)
!!$      ALLOCATE(R(NR))
!!$      ALLOCATE(FEXPR2N(NR))
!!$      CALL RADIAL$R(GID,NR,R)
!      
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
      NAT=SIZE(GAUSSORBAUG)
      IF(.NOT.ALLOCATED(LMORB))ALLOCATE(LMORB(NAT))
      ALLOCATE(VEC((LX+1)**2))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL RADIAL$GETI4(GID,'NR',NR)
        ALLOCATE(R(NR))
        ALLOCATE(FEXPR2N(NR))
        CALL RADIAL$R(GID,NR,R)

        NIJK=GAUSSORBAUG(IAT)%NIJK
        CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJK,I,J,K)
        LX=MAX(I,J,K)
        LMX=(LX+1)**2
        NE=GAUSSORBAUG(IAT)%NE
        NORB=GAUSSORBAUG(IAT)%NORB
        LMORB(IAT)%GID=GID
        LMORB(IAT)%NR=NR
        LMORB(IAT)%LMX=LMX
        LMORB(IAT)%NORB=NORB
        ALLOCATE(LMORB(IAT)%F(NR,LMX,NORB))
        LMORB(IAT)%F(:,:,:)=0.D0 
        DO IE=1,NE
          FEXPR2N(:)=EXP(-GAUSSORBAUG(IAT)%E(IE)*R(:)**2)
          DO N=0,LX/2
            LMX1=(LX-2*N+1)**2
!            LMX1=(LX+1)**2
            DO IORB=1,NORB
              VEC(:LMX1)=MATMUL(YLMFROMPOL(N+1,:LMX1,:) &
      &                                             ,GAUSSORBAUG(IAT)%C(:,IE,IORB))
!!$IF(IORB.EQ.1) THEN
!!$  PRINT*,'IAT  ',IAT,' IE=',IE,' N=',N,' LMX1=',LMX1,' LX=',LX
!!$  PRINT*,'C  ',GAUSSORBAUG(IAT)%C(:,IE,IORB)
!!$  PRINT*,'VEC ',VEC(:LMX1)
!!$END IF
              DO LM=1,LMX1
                LMORB(IAT)%F(:,LM,IORB)=LMORB(IAT)%F(:,LM,IORB) &
      &                                                      +FEXPR2N(:)*VEC(LM)
              ENDDO
            ENDDO
            FEXPR2N(:)=FEXPR2N(:)*R(:)**2
          ENDDO
        ENDDO
!
!       == MULTIPLICATION WITH R^L =============================================
        LM=1
        DO L=1,LX  ! STARTS WITH L=1 ON PURPOSE..
          FEXPR2N(:)=R(:)**L
          DO IM=1,2*L+1
            LM=LM+1
            DO IORB=1,NORB
              LMORB(IAT)%F(:,LM,IORB)=LMORB(IAT)%F(:,LM,IORB)*FEXPR2N(:)
            ENDDO
          ENDDO
        ENDDO        
DO IORB=1,NORB
  WRITE(STRING,FMT='(I5)')IAT
  WRITE(STRING1,FMT='(I5)')IORB
  STRING='NSCHI_FORATOM'//TRIM(ADJUSTL(STRING))//'ANDORB'
  STRING=TRIM(ADJUSTL(STRING))//TRIM(ADJUSTL(STRING1))//'.DAT'
  LM=MIN(16,LX**2)
  CALL SETUP_WRITEPHI(TRIM(STRING),GID,NR,LM,LMORB(IAT)%F(:,:LM,IORB))
ENDDO

!CALL SETUP_WRITEPHI('TEST.DAT',GID,NR,16,LMORB(IAT)%F(:,:16,1))
!!$PRINT*,'LMX ',LMX
!!$CALL SETUP_WRITEPHI('TEST.DAT',GID,NR,MIN(LMX,16),LMORB(IAT)%F(:,:MIN(LMX,16),1))
!!$ALLOCATE(YLM((LX+1)**2))
!!$CALL SPHERICAL$YLM((LX+1)**2,(/1.D0,0.D0,0.D0/),YLM)
!!$FEXPR2N(:)=0.D0
!!$DO LM=1,(LX+1)**2
!!$  FEXPR2N(:)=FEXPR2N(:)+LMORB(IAT)%F(:,LM,1)*YLM(LM)
!!$ENDDO    
!!$CALL SETUP_WRITEPHI('TEST2.DAT',GID,NR,1,FEXPR2N)
!!$ FEXPR2N(:)=0.D0
!!$ DO IND=1,NIJK
!!$   DO IE=1,NE
!!$     CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I,J,K)
!!$     IF(J.NE.0.OR.K.NE.0) CYCLE
!!$     FEXPR2N(:)=FEXPR2N(:)+EXP(-GAUSSORBAUG(IAT)%E(IE)*R(:)**2)*R(:)**I &
!!$&                                    *GAUSSORBAUG(IAT)%C(IND,IE,1)
!!$   ENDDO
!!$ ENDDO    
!!$CALL SETUP_WRITEPHI('TEST1.DAT',GID,NR,1,FEXPR2N)
!!$STOP 'FORCED IN LMTO_GAUSSORBAUGTOYLM'
        DEALLOCATE(R)
        DEALLOCATE(FEXPR2N)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_CPGAUSSORB(NAT,FROM,TO)
!     **************************************************************************
!     ** COPY A STRUCTURE OF TYPE ORBITALGAUSSCOEFF_TYPE                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ORBITALGAUSSCOEFF_TYPE
      IMPLICIT NONE
      INTEGER(4)                  ,INTENT(IN)    :: NAT
      TYPE(ORBITALGAUSSCOEFF_TYPE),INTENT(IN)    :: FROM(NAT)
      TYPE(ORBITALGAUSSCOEFF_TYPE),INTENT(INOUT) :: TO(NAT)
      INTEGER(4)                                 :: IAT,NIJK,NE,NORB
!     **************************************************************************
      DO IAT=1,NAT
        NIJK=FROM(IAT)%NIJK
        NE  =FROM(IAT)%NE
        NORB=FROM(IAT)%NORB
        DEALLOCATE(TO(IAT)%E)
        DEALLOCATE(TO(IAT)%C)
        ALLOCATE(TO(IAT)%E(NE))
        ALLOCATE(TO(IAT)%C(NIJK,NE,NORB))
        TO(IAT)%NIJK=NIJK
        TO(IAT)%NE  =NE
        TO(IAT)%NORB=NORB
        TO(IAT)%E=FROM(IAT)%E
        TO(IAT)%C=FROM(IAT)%C
      ENDDO
      RETURN
      END



!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE LMTO$NTBOS(NPRO,NPSI,C,NP,P,PSI)
!!$!     **************************************************************************
!!$      USE LMTO_MODULE
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: NPSI
!!$      INTEGER(4),INTENT(IN) :: NPRO
!!$      REAL(8)   ,INTENT(IN) :: C(NPRO,NPSI)
!!$      INTEGER(4),INTENT(IN) :: NP
!!$      REAL(8)   ,INTENT(IN) :: P(3,NP)   ! GRID POINTS
!!$      REAL(8)   ,INTENT(OUT):: PSI(NP,NPSI)
!!$      REAL(8)   ,ALLOCATABLE:: R0(:,:)
!!$      REAL(8)               :: RBAS(3,3)
!!$      INTEGER(4)            :: NAT
!!$      INTEGER(4)            :: NIJK,NE,NORB
!!$      REAL(8)               :: EI
!!$      REAL(8)               :: RI(3)
!!$      REAL(8)               :: SVAR
!!$!     **************************************************************************
!!$!
!!$!     ==========================================================================
!!$!     ==  COLLECT DATA                                                        ==
!!$!     ==========================================================================
!!$      CALL ATOMLIST$NATOM(NAT)
!!$      ALLOCATE(R0(3,NAT))
!!$      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
!!$      CALL CELL$GETR8A('T0',9,RBAS)
!!$!
!!$!     ==========================================================================
!!$!     ==  ACCUMULATE ENVELOPE FUNCTIONS OF NTBOS                              ==
!!$!     ==========================================================================
!!$      DO IAT=1,NAT
!!$        NIJK=GAUSSORB(IAT)%NIJK
!!$        NE=GAUSSORB(IAT)%NE
!!$        NORB=GAUSSORB(IAT)%NORB
!!$        RI(:)=R0(:,IAT)
!!$        DO IORB=1,NORB
!!$          DO IE=1,NE
!!$            EI=GAUSSORB(IAT)%E(IE)
!!$            CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,EI,GAUSSORB(IAT)%C(:,IE,IORB),RI,SVAR)
!!$            PSI(IP,:)=PSI(IP,:)+SVAR*SMALLVEC(IORB,:,IAT)
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_TESTOVERLAP()
!     **************************************************************************
!     ** TESTS THE OVERLAP MATRIX OF NATURAL TIGHT-BINDING ORBITALS           **
!     ** BY ESTIMATING THE OVERLAP BETWEEN KOHN SHAM WAVE FUNCTIONS USING THE **
!     ** COEFFICIENTS IN NTBS AND THEIR OVERLAP MATRIX                        **
!     **************************************************************************
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV,GSET
      USE LMTO_MODULE, ONLY : PERIODICMAT_TYPE,SBAR,ISPECIES,POTPAR &
     &                       ,SBARATOMI1,SBARATOMI2,SBARLI1,LNX,LOX &
     &                       ,OVERLAP
      IMPLICIT NONE
      REAL(8)       :: XK(3,NKPTL)
      INTEGER(4)    :: NPRO
      INTEGER(4)    :: NBH,NB
      INTEGER(4)    :: NNS
      INTEGER(4)    :: IKPT,ISPIN,IBH,JBH,IPRO,JPRO,IDIM,IB
      COMPLEX(8)    :: CSVAR1,CSVAR2
      COMPLEX(8),ALLOCATABLE :: OVER(:,:)
      REAL(8)   ,ALLOCATABLE :: OVERMAT(:,:)
      LOGICAL(4)             :: TINV
      LOGICAL(4),parameter   :: Tpr=.false.
!     **************************************************************************
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      NPRO=MAP%NPRO
      NNS=SIZE(OVERLAP)
!
!     ==========================================================================
!     ==  DETERMINE OVERLAP MATRIX  IN K-SPACE                                ==
!     ==========================================================================
!!!! PARALLELIZE LOOP WITH RESPECT TO STATES.
      ALLOCATE(OVER(NPRO,NPRO))
      DO IKPT=1,NKPTL
        CALL WAVES_SELECTWV(IKPT,1)
        CALL PLANEWAVE$SELECT(GSET%ID)
        NBH=THIS%NBH
        NB=THIS%NB
        CALL PLANEWAVE$GETL4('TINV',TINV)
        CALL LMTO_MATRTOK(NNS,OVERLAP,XK(:,IKPT),NPRO,OVER)
if(tpr) then
  WRITE(*,FMT='(82("="),T20," XK=",3F10.5,"  ")')XK(:,IKPT)
  PRINT*,'===================== OVERLAP (REAL PART) ===================='
  DO JPRO=1,NPRO
    WRITE(*,FMT='("R(O)",I3,180F10.5)')JPRO,REAL(OVER(:,JPRO))
  ENDDO
  PRINT*,'===================== OVERLAP (IMAGINARY PART) ==============='
  DO JPRO=1,NPRO
    WRITE(*,FMT='("I(O)",180F10.5)')AIMAG(OVER(:,JPRO))
  ENDDO
end if
        DO ISPIN=1,NSPIN
!!$PRINT*,'===================== ISPIN ',ISPIN,'========================='
          CALL WAVES_SELECTWV(IKPT,ISPIN)

PRINT*,'=============WAVE FUNCTION COEFFICIENTS FOR NTBOS  ================='
if(tpr) then
  IF(TINV) THEN
    PRINT*,'SIZE ',SIZE(THIS%TBC),' SHAPE ',SHAPE(THIS%TBC) 
    DO IBH=1,NBH 
      WRITE(*,FMT='("C",I5,180F10.5)')2*IBH-1,REAL(THIS%TBC(1,IBH,:))
      WRITE(*,FMT='("C",I5,180F10.5)')2*IBH,AIMAG(THIS%TBC(1,IBH,:))
    ENDDO
  ELSE 
    DO IB=1,NB
      WRITE(*,FMT='("RE[C]",I5,180F10.5)')IB,REAL(THIS%TBC(1,IB,:))
      WRITE(*,FMT='("IM[C]",I5,180F10.5)')IB,AIMAG(THIS%TBC(1,IB,:))
    ENDDO 
  END IF
end if
          ALLOCATE(OVERMAT(NB,NB))
          OVERMAT(:,:)=0.D0
          DO IBH=1,NBH
            DO JBH=1,NBH
              IF(TINV) THEN
                CSVAR1=(0.D0,0.D0)
                CSVAR2=(0.D0,0.D0)
                DO IPRO=1,NPRO
                  DO JPRO=1,NPRO
                    DO IDIM=1,NDIM
                      CSVAR1=CSVAR1+CONJG(THIS%TBC(IDIM,IBH,IPRO)) &
      &                                 *OVER(IPRO,JPRO)*THIS%TBC(IDIM,JBH,JPRO)
                      CSVAR2=CSVAR2+THIS%TBC(IDIM,IBH,IPRO) &
      &                                 *OVER(IPRO,JPRO)*THIS%TBC(IDIM,JBH,JPRO)
                    ENDDO
                  ENDDO
                ENDDO
!PRINT*,'CSVAR ',IBH,JBH,CSVAR1,CSVAR2
                OVERMAT(2*IBH-1,2*JBH-1)= 0.5D0* REAL(CSVAR1+CSVAR2)
                OVERMAT(2*IBH-1,2*JBH)  = 0.5D0*AIMAG(CSVAR1+CSVAR2)
                OVERMAT(2*IBH  ,2*JBH-1)=-0.5D0*AIMAG(CSVAR1-CSVAR2)
                OVERMAT(2*IBH  ,2*JBH)  = 0.5D0* REAL(CSVAR1-CSVAR2)
              ELSE
                CSVAR1=(0.D0,0.D0)
                DO IPRO=1,NPRO
                  DO JPRO=1,NPRO
                    DO IDIM=1,NDIM
                      CSVAR1=CSVAR1+CONJG(THIS%TBC(IDIM,IBH,IPRO)) &
      &                                 *OVER(IPRO,JPRO)*THIS%TBC(IDIM,JBH,JPRO)
                    ENDDO
                  ENDDO
                ENDDO
                OVERMAT(IBH,JBH)=REAL(CSVAR1)
              END IF
            ENDDO
          ENDDO
if(tpr) then
  DO IB=1,NB
    WRITE(*,FMT='("UNITY ",I5,200F10.5)')IB,OVERMAT(:,IB)
  ENDDO
  WRITE(*,*)
end if
!STOP 'FORCED IN OVERLAPFULL'
          DEALLOCATE(OVERMAT)
        ENDDO
      ENDDO
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
!!$print*,'marke 1',ikpt,ispin,ibh
!!$print*,'marke 1a',size(this%tbc)
!!$print*,'marke 1b',size(tbc1)
!!$print*,'marke 1c',size(vec1)
!!$print*,'marke 1d',size(vec2)
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
      SUBROUTINE LMTO_EXPANDOVERLAPOFK(NSMALL,NBIG,OVS,SBAR,OVB)
!     **************************************************************************
!     ** THE (SMALL) OBERLAP MATRIX OVS IS THE ONE BETWEEN THE MAIN ANGULAR   **
!     ** MOMENTA, THAT IS THERE IS ONE ENTRY PER SITE AND ANGULAR MOMENTUM.   **
!     ** THIS ROUTINE BLOWS THE OVERLAP MATRIX UP TO THE FULL SPACE OF        **
!     ** NATURAL TIGHT-BINDING ORBITALS                                       **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : LNX,LOX,POTPAR,SBARLI1,SBARATOMI1,ISPECIES &
     &                      ,POTPAR_TYPE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSMALL
      INTEGER(4),INTENT(IN) :: NBIG
      COMPLEX(8),INTENT(IN) :: OVS(NSMALL,NSMALL)
      COMPLEX(8),INTENT(IN) :: SBAR(NSMALL,NSMALL)
      COMPLEX(8),INTENT(OUT):: OVB(NBIG,NBIG)
      INTEGER(4)            :: IPRO1,IAT1,LN1,L1,LM1,ISP1,LN1SCATT,LM1A,LM1B
      INTEGER(4)            :: IPRO2,IAT2,LN2,L2,LM2,ISP2,LN2SCATT,LM2A,LM2B
      INTEGER(4)            :: NAT
      REAL(8)               :: SVAR
!     **************************************************************************
      NAT=SIZE(ISPECIES)
      OVB(:,:)=(0.D0,0.D0)
      IPRO1=0
      DO IAT1=1,NAT
        ISP1=ISPECIES(IAT1)          
        DO LN1=1,LNX(ISP1)
          LN1SCATT=POTPAR(ISP1)%LNSCATT(LN1)
          L1=LOX(LN1,ISP1)
          LM1A=SBARATOMI1(IAT1)-1+SBARLI1(L1+1,ISP1)
          LM1B=LM1A+2*L1
          DO LM1=LM1A,LM1B
            IPRO1=IPRO1+1           
!
!           == LOOP OVER SECOND INDEX ==========================================
            IPRO2=0
            DO IAT2=1,NAT
              ISP2=ISPECIES(IAT2)          
              DO LN2=1,LNX(ISP2)
                LN2SCATT=POTPAR(ISP2)%LNSCATT(LN2)
                L2=LOX(LN2,ISP2)
                LM2A=SBARATOMI1(IAT2)-1+SBARLI1(L2+1,ISP2)
                LM2B=LM2A+2*L2
                DO LM2=LM2A,LM2B
                  IPRO2=IPRO2+1           
!                 == COPY FROM OVERLAP =========================================
                  OVB(IPRO1,IPRO2)=OVB(IPRO1,IPRO2)+OVS(LM1,LM2)
!                 == FIX <K|JBAR>SBAR ==========================================
                 SVAR=POTPAR(ISP1)%DOVERLAPKJ(LN1,LN2SCATT)  &
      &              -POTPAR(ISP1)%DOVERLAPKJ(LN1SCATT,LN2SCATT)
                  OVB(IPRO1,IPRO2)=OVB(IPRO1,IPRO2)+SVAR*SBAR(LM1,LM2)
!                 == FIX SBAR<JBAR|K> ==========================================
                 SVAR=POTPAR(ISP1)%DOVERLAPKJ(LN2SCATT,LN1) &
      &               -POTPAR(ISP1)%DOVERLAPKJ(LN2SCATT,LN1SCATT)
                  OVB(IPRO1,IPRO2)=OVB(IPRO1,IPRO2)+CONJG(SBAR(LM2,LM1))*SVAR
!                 == FIX <K|K> =================================================
                  IF(IAT1.EQ.IAT2) THEN
                    IF(L1.EQ.L2) THEN
                      IF(LM2-LM2A.EQ.LM1-LM1A) THEN
                        SVAR=POTPAR(ISP1)%DOVERLAPKK(LN1,LN2) &
      &                      -POTPAR(ISP1)%DOVERLAPKK(LN1SCATT,LN2SCATT)
                         OVB(IPRO1,IPRO2)=OVB(IPRO1,IPRO2)+SVAR
                      END IF
                    END IF
                  END IF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_EXPANDSPHEREDO(ID)
!     **************************************************************************
!     **  CALCULATES  THE CORRECTION OF THE OVERLAP RELATIVE TO THE SMOOTH    **
!     **  PART EXPANDED INTO GAUSSIANS. DONE ONLY FOR THE LEADING CHANNEL     **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : NSP,LNX,LOX,POTPAR,SBARLI1
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4) :: LX
      INTEGER(4) :: LMX
      INTEGER(4) :: L
      INTEGER(4) :: I1,I2
      INTEGER(4) :: ISP,LN,LM
!     **************************************************************************
      IF(ID.EQ.'CLEAN') THEN
        DO ISP=1,NSP
          DEALLOCATE(POTPAR(ISP)%SMALL%DOVERLAPKK)
          DEALLOCATE(POTPAR(ISP)%SMALL%DOVERLAPKJ)
          DEALLOCATE(POTPAR(ISP)%SMALL%DOVERLAPJJ)
        ENDDO
      ELSE IF(ID.EQ.'LM') THEN
      ELSE
        CALL ERROR$MSG('IDENTIFIER NOT RECOGNIZED')
        CALL ERROR$STOP('LMTO_EXPANDSPHEREDO')
      END IF
!
!     == USE THE MAPPING USED IN LMTO$MAKESTRUCTURECONSTANTS
      DO ISP=1,NSP
        LX=MAXVAL(LOX(:,ISP))
        LMX=SBARLI1(LX+1,ISP)+2*LX
        ALLOCATE(POTPAR(ISP)%SMALL%DOVERLAPKK(LMX,LMX))
        ALLOCATE(POTPAR(ISP)%SMALL%DOVERLAPKJ(LMX,LMX))
        ALLOCATE(POTPAR(ISP)%SMALL%DOVERLAPJJ(LMX,LMX))
        POTPAR(ISP)%SMALL%DOVERLAPKK(:,:)=0.D0
        POTPAR(ISP)%SMALL%DOVERLAPKJ(:,:)=0.D0
        POTPAR(ISP)%SMALL%DOVERLAPJJ(:,:)=0.D0
        DO LN=1,LNX(ISP)
!         == LNSCATT POINT TO THE VALENCE CHANNEL FROM WHICH THE CORRESPONDING
!         == SCATTERING WAVE FUNCTION IS TAKEN
          IF(POTPAR(ISP)%LNSCATT(LN).NE.LN) CYCLE
          L=LOX(LN,ISP)
          I1=SBARLI1(L+1,ISP)
          I2=I1+2*L
          DO LM=I1,I2
            POTPAR(ISP)%SMALL%DOVERLAPKK(LM,LM)=POTPAR(ISP)%DOVERLAPKK(LN,LN)
            POTPAR(ISP)%SMALL%DOVERLAPKJ(LM,LM)=POTPAR(ISP)%DOVERLAPKJ(LN,LN)
            POTPAR(ISP)%SMALL%DOVERLAPJJ(LM,LM)=POTPAR(ISP)%DOVERLAPJJ(LN,LN)
          ENDDO
        ENDDO
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_ONSITEOVERLAP()
!     **************************************************************************
!     **  CONSTRUCT THE DIFFERENCE OF THE OVERLAP OF PARTIAL WAVES MATCHING  **
!     **  THE HANKEL AND SCREENED BESSEL FUNCTIONS.                            **
!     **                                                                      **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : POTPAR,NSP,LOX,LNX,K2
      IMPLICIT NONE
      INTEGER(4)             :: LNX1
      INTEGER(4)             :: LMNX,LX
      INTEGER(4)             :: NR        ! #(RADIAL GRID POINTS)
      INTEGER(4)             :: GID       ! GRID-IDENTIFIER
      REAL(8)                :: SVAR      ! AUXILIART VARIABLE
      REAL(8)   ,ALLOCATABLE :: AUX(:)    ! AUXILIARY ARRAY
      REAL(8)   ,ALLOCATABLE :: R(:)      ! RADIAL GRID
      REAL(8)   ,ALLOCATABLE :: KIN(:,:)
      REAL(8)   ,ALLOCATABLE :: JIN(:,:)
      REAL(8)   ,ALLOCATABLE :: KOUT(:,:)
      REAL(8)   ,ALLOCATABLE :: JOUT(:,:)
      REAL(8)   ,ALLOCATABLE :: JBESSEL(:,:)
      REAL(8)   ,ALLOCATABLE :: JBAR(:,:)
      REAL(8)   ,ALLOCATABLE :: KHANKEL(:,:)
      REAL(8)                :: RAD
      REAL(8)                :: KVAL,KDER,VAL,DER
      INTEGER(4)             :: ISP,LN,LN1,LN2,LMN2,L,L1,L2,IR
!     **************************************************************************
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        LNX1=LNX(ISP)
        LMNX=SUM(2*LOX(:LNX1,ISP)+1)
        LX=MAXVAL(LOX(:LNX1,ISP))
        RAD=POTPAR(ISP)%RAD
!       ========================================================================
!       == CONSTRUCT FUNCTIONS MATCHING ONTO BESSEL AND HANKEL FUNCTIONS.     ==
!       == KIN AND JIN ARE SUPERPOSITIONS OF BESSEL FUNCTIONS AND NODELESS    ==
!       == SCATTERING FUNCTIONS MATCHING TO HANKEL AND BESSEL FUNCTIONS.      ==
!       == KOUT AND JOUT ARE MATCHING ALL-ELECTRON PARTIAL WAVES.             ==
!       ========================================================================
        ALLOCATE(KIN(NR,LNX1))
        ALLOCATE(JIN(NR,LNX1))
        ALLOCATE(KOUT(NR,LNX1))
        ALLOCATE(JOUT(NR,LNX1))
        CALL SETUP$GETR8A('QPHI',NR*LNX1,KOUT)    ! NODELESS ATOMIC
        CALL SETUP$GETR8A('AEPHI',NR*LNX1,KIN)    ! NODAL ATOMIC
        CALL SETUP$GETR8A('QPHIDOT',NR*LNX1,JOUT) ! NODELESS SCATTERING
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX1,JIN) ! NODAL SCATTERING
!
!       ========================================================================
!       == ONE SCATTERING FUNCTION PER L IS USED FOR ALL PARTIAL WAVES        ==
!       ========================================================================
        DO LN=1,LNX(ISP)
          JOUT(:,LN)=JOUT(:,POTPAR(ISP)%LNSCATT(LN))
          JIN(:,LN)=JIN(:,POTPAR(ISP)%LNSCATT(LN))
        ENDDO
!
!       ========================================================================
!       == CONSTRUCT BESSEL FUNCTIONS ==========================================
!       ========================================================================
        ALLOCATE(JBESSEL(NR,LX+1))
        ALLOCATE(KHANKEL(NR,LX+1))
        KHANKEL=0.D0
        DO L=0,LX
          DO IR=1,NR
            CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,JBESSEL(IR,L+1),SVAR)
             IF(R(IR).GT.1.D0) CALL LMTO$SOLIDHANKELRAD(L,R(IR),K2,KHANKEL(IR,L+1),SVAR)
          ENDDO 
        ENDDO

        ALLOCATE(JBAR(NR,LNX1))
        DO LN=1,LNX1
          L=LOX(LN,ISP)
          KIN(:,LN)=KIN(:,LN)*POTPAR(ISP)%KTOPHI(LN) &
       &           +JIN(:,LN)*POTPAR(ISP)%KTOPHIDOT(LN)
!!$CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
!!$CALL RADIAL$VALUE(GID,NR,KIN(:,LN),RAD,VAL)
!!$CALL RADIAL$DERIVATIVE(GID,NR,KIN(:,LN),RAD,DER)
!!$PRINT*,'TEST IN  ',LN,KVAL,VAL,KDER,DER
          KOUT(:,LN)=JBESSEL(:,L+1)-JOUT(:,LN)*POTPAR(ISP)%JBARTOPHIDOT(LN)
          KOUT(:,LN)=KOUT(:,LN)/POTPAR(ISP)%QBAR(LN)
!!$CALL RADIAL$VALUE(GID,NR,KOUT(:,LN),RAD,VAL)
!!$CALL RADIAL$DERIVATIVE(GID,NR,KOUT(:,LN),RAD,DER)
!!$PRINT*,'TEST OUT ',LN,KVAL,VAL,KDER,DER
          JIN(:,LN) =JIN(:,LN)*POTPAR(ISP)%JBARTOPHIDOT(LN)
          JOUT(:,LN)=JOUT(:,LN)*POTPAR(ISP)%JBARTOPHIDOT(LN)
          JBAR(:,LN)=JBESSEL(:,L+1)-KHANKEL(:,L+1)*POTPAR(ISP)%QBAR(LN)
        ENDDO
!!$CALL SETUP_WRITEPHI('JBESSEL.DAT',GID,NR,LX+1,JBESSEL)
!!$CALL SETUP_WRITEPHI('KHANKEL.DAT',GID,NR,LX+1,KHANKEL)
!!$CALL SETUP_WRITEPHI('JBAR.DAT',GID,NR,LNX1,JBAR)
        DEALLOCATE(JBAR)
!
!       == FIX TAILS TO AVOID ERRORS IN THE INTEGRATION =====================
        DO IR=1,NR
          IF(R(IR).GT.RAD) THEN
            DO LN=1,LNX1
             L=LOX(LN,ISP)
              KIN(IR:,LN) =KHANKEL(IR:,L+1)
              KOUT(IR:,LN)=KHANKEL(IR:,L+1)
              JIN(IR:,LN) =JBESSEL(IR:,L+1)
              JOUT(IR:,LN)=JBESSEL(IR:,L+1)
            ENDDO
            EXIT
          END IF
        ENDDO
        DEALLOCATE(JBESSEL)
        DEALLOCATE(KHANKEL)
!!$CALL SETUP_WRITEPHI('KIN.DAT',GID,NR,LNX1,KIN)
!!$CALL SETUP_WRITEPHI('KOUT.DAT',GID,NR,LNX1,KOUT)
!!$CALL SETUP_WRITEPHI('JIN.DAT',GID,NR,LNX1,JIN)
!!$CALL SETUP_WRITEPHI('JOUT.DAT',GID,NR,LNX1,JOUT)
!!$PRINT*,'RAD ',RAD
!!$STOP 'FORCED '
        ALLOCATE(POTPAR(ISP)%DOVERLAPKK(LNX1,LNX1))
        ALLOCATE(POTPAR(ISP)%DOVERLAPKJ(LNX1,LNX1))
        ALLOCATE(POTPAR(ISP)%DOVERLAPJJ(LNX1,LNX1))
        POTPAR(ISP)%DOVERLAPKK=0.D0
        POTPAR(ISP)%DOVERLAPKJ=0.D0
        POTPAR(ISP)%DOVERLAPJJ=0.D0
        ALLOCATE(AUX(NR))
        DO LN1=1,LNX1
          L1=LOX(LN1,ISP)
          DO LN2=1,LNX1
            L2=LOX(LN2,ISP)
            IF(L2.EQ.L1) THEN
!             == KKOVERLAP =====================================================
              AUX(:)=(KIN(:,LN1)*KIN(:,LN2)-KOUT(:,LN1)*KOUT(:,LN2))*R(:)**2
              CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
              POTPAR(ISP)%DOVERLAPKK(LN1,LN2)=SVAR
!             == KJOVERLAP =====================================================
              AUX(:)=(KIN(:,LN1)*JIN(:,LN2)-KOUT(:,LN1)*JOUT(:,LN2))*R(:)**2
              CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
              POTPAR(ISP)%DOVERLAPKJ(LN1,LN2)=SVAR
!             == JJOVERLAP =====================================================
              AUX(:)=(JIN(:,LN1)*JIN(:,LN2)-JOUT(:,LN1)*JOUT(:,LN2))*R(:)**2
              CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
              POTPAR(ISP)%DOVERLAPJJ(LN1,LN2)=SVAR
            END IF
          ENDDO
        ENDDO
PRINT*,'============================================='
DO LN1=1,LNX1
  PRINT*,'<KK>', POTPAR(ISP)%DOVERLAPKK(LN1,:)
ENDDO
PRINT*,'============================================='
DO LN1=1,LNX1
  PRINT*,'<KJ>', POTPAR(ISP)%DOVERLAPKJ(LN1,:)
ENDDO
PRINT*,'============================================='
DO LN1=1,LNX1
  PRINT*,'<JJ>', POTPAR(ISP)%DOVERLAPJJ(LN1,:)
ENDDO
        DEALLOCATE(KIN)
        DEALLOCATE(JIN)
        DEALLOCATE(KOUT)
        DEALLOCATE(JOUT)
        DEALLOCATE(R)
        DEALLOCATE(AUX)
      ENDDO  !END LOOP OVER ATOM TYPES
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OVERLAPOFK(XK,NPRO,OVER)
!     **************************************************************************
!     **  CONSTRUCT OVERLAP MATRIX IN K-SPACE
!     **************************************************************************
      USE LMTO_MODULE, ONLY: POTPAR_TYPE,OVERLAP,SBAR,SBARATOMI2
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: XK(3)
      INTEGER(4),INTENT(IN) :: NPRO
      COMPLEX(8),INTENT(OUT):: OVER(NPRO,NPRO)
      INTEGER(4)            :: NNO
      INTEGER(4)            :: NNS
      INTEGER(4)            :: NSMALL
      COMPLEX(8),ALLOCATABLE:: OVERLAPOFK(:,:)
      COMPLEX(8),ALLOCATABLE:: SBAROFK(:,:)
!     **************************************************************************
      NSMALL=MAXVAL(SBARATOMI2)
!
!     ==========================================================================
!     ==  TRANSFORM OVERLAP MATRIX TO K-SPACE                                 ==
!     ==========================================================================
      NNO=SIZE(OVERLAP)
      ALLOCATE(OVERLAPOFK(NSMALL,NSMALL))
      CALL LMTO_AOFK(NNO,OVERLAP,XK,NSMALL,OVERLAPOFK)
!
!     ==========================================================================
!     ==  TRANSFORM SCREENED STRUCTURE CONSTANTS TO K-SPACE                   ==
!     ==========================================================================
      NNS=SIZE(SBAR)
      ALLOCATE(SBAROFK(NSMALL,NSMALL))
      CALL LMTO_AOFK(NNS,SBAR,XK,NSMALL,SBAROFK)
!
!     ==========================================================================
!     ==  SET UP COMPLETE OVERLAP MATRIX                                      ==
!     ==========================================================================
      CALL LMTO_EXPANDOVERLAPOFK(NSMALL,NPRO,OVERLAPOFK,SBAROFK,OVER)
      DEALLOCATE(OVERLAPOFK)
      DEALLOCATE(SBAROFK)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_NTBOAUGMENT()
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : UTENSOR,NSP,ISPECIES,LNX,LOX,GAUSSORB &
     &                       ,GAUSSORBAUG,SBARLI1,POTPAR,SBAR
      IMPLICIT NONE
      INTEGER(4) :: NAT
      INTEGER(4) :: IAT
      INTEGER(4) :: NIJKA,NIJKG,NIJK
      INTEGER(4) :: NEA,NEG,NE
      INTEGER(4) :: NORB
      INTEGER(4) :: ISP,LMN,LN,LN1,LN2,L,L1,L2,LM,LM1,LM2,IORB1,IORB2,NN,NNS
!     **************************************************************************
      NAT=SIZE(ISPECIES)
      IF(.NOT.ALLOCATED(GAUSSORBAUG))ALLOCATE(GAUSSORBAUG(NAT))
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        NIJKG=GAUSSORB(IAT)%NIJK
        NIJKA=POTPAR(ISP)%GAUSSKAUGMENT%NIJK
        NEG=GAUSSORB(IAT)%NE
        NEA=POTPAR(ISP)%GAUSSKAUGMENT%NE
!
        NIJK=MAX(NIJKG,NIJKA)
        NE=NEA+NEG
        NORB=POTPAR(ISP)%GAUSSKAUGMENT%NORB
        GAUSSORBAUG(IAT)%NIJK=NIJK
        GAUSSORBAUG(IAT)%NE=NE
        GAUSSORBAUG(IAT)%NORB=NORB
        ALLOCATE(GAUSSORBAUG(IAT)%E(NE))
        GAUSSORBAUG(IAT)%E(:NEA)=POTPAR(ISP)%GAUSSKAUGMENT%E
        GAUSSORBAUG(IAT)%E(NEA+1:)=GAUSSORB(IAT)%E
        ALLOCATE(GAUSSORBAUG(IAT)%C(NIJK,NE,NORB))
!       == MAP 
        GAUSSORBAUG(IAT)%C(:,:,:)=0.D0
        LMN=0
        DO LN=1,LNX(ISP)
          LN1=POTPAR(ISP)%LNSCATT(LN)
          L=LOX(LN,ISP)
          DO LM=SBARLI1(L+1,ISP),SBARLI1(L+1,ISP)+2*L
            LMN=LMN+1
            GAUSSORBAUG(IAT)%C(:NIJKG,NEA+1:,LMN)=GAUSSORB(IAT)%C(:,:,LM)
          ENDDO
        ENDDO
!       == AUGMENT K ===========================================================
        GAUSSORBAUG(IAT)%C(:NIJKA,:NEA,:)=POTPAR(ISP)%GAUSSKAUGMENT%C
      ENDDO
!
!     ==========================================================================
!     == NOW AUGMENT THE ONSITE SCREENED BESSEL FUNCTIONS                     ==
!     ==========================================================================
      NNS=SIZE(SBAR)
      DO NN=1,NNS
        IF(SBAR(NN)%IAT1.NE.SBAR(NN)%IAT2) CYCLE
        IF(SBAR(NN)%IT(1).NE.0) CYCLE
        IF(SBAR(NN)%IT(2).NE.0) CYCLE
        IF(SBAR(NN)%IT(3).NE.0) CYCLE
        IAT=SBAR(NN)%IAT1
        ISP=ISPECIES(IAT)
        IORB1=0
        DO LN1=1,LNX(ISP)
          L1=LOX(LN1,ISP)
          DO LM1=SBARLI1(L1+1,ISP),SBARLI1(L1+1,ISP)+2*L1
            IORB1=IORB1+1
            DO LN2=1,LNX(ISP)
              IF(POTPAR(ISP)%LNSCATT(LN2).NE.LN2) CYCLE
              L2=LOX(LN2,ISP)
              DO LM2=SBARLI1(L2+1,ISP),SBARLI1(L2+1,ISP)+2*L2
                GAUSSORBAUG(IAT)%C(:NIJKA,:NEA,IORB1) &
     &               =GAUSSORBAUG(IAT)%C(:NIJKA,:NEA,IORB1) &
     &               -POTPAR(ISP)%GAUSSJAUGMENT%C(:,:,LM2)*SBAR(NN)%MAT(LM2,LM1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
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
                            CALL TRACE$PUSH('LDAPLUSU_ULITTLE')
      CALL RADIAL$R(GID,NR,R)
      ULITTLE=0.D0
      DO LN1=1,LNX
        DO LN2=LN1,LNX
          RHO(:)=CHI(:,LN1)*CHI(:,LN2)
!         == USE SELECTION RULE (NOT TO SAVE TIME HERE, BUT LATER FOR THE U-TENSOR)
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
!IF(LOX(LN1).NE.LOX(LN2).OR.LOX(LN2).NE.LOX(LN3).OR.LOX(LN3).NE.LOX(LN4)) SVAR=0.D0
!IF(LOX(LN1)*LOX(LN2)*LOX(LN3)*LOX(LN4).EQ.0) SVAR=0.D0
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
      SUBROUTINE LMTO_FOURCENTERGAUSS()
!     **************************************************************************
!     **  EVALUATES THE U-TENSOR                                              **
!     **                                                                      **
!     **  PRELIMINARY VERSION: ONLY ONSITE MATRIX ELEMENTS                    **
!     **  PRELIMINARY VERSION: NO AUGMENTATION CONTRIBUTIONS INCLUDED         **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : UTENSOR,ISPECIES,LNX,LOX,GAUSSORB,SBARLI1
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      REAL(8)                :: RBAS(3,3)
      INTEGER(4)             :: N1,N2,N3,N4
      INTEGER(4)             :: IAT1,IAT2,IAT3,IAT4
      INTEGER(4)             :: NIJK1,NIJK2,NIJK3,NIJK4
      INTEGER(4)             :: NORB1,NORB2,NORB3,NORB4
      INTEGER(4)             :: NE1,NE2,NE3,NE4
      REAL(8)                :: R1(3),R2(3),R3(3),R4(3)
      REAL(8)   ,ALLOCATABLE :: E1(:),E2(:),E3(:),E4(:)
      INTEGER(4)             :: IAT,NN,ISP,IND,LN,LM,L,IM,I
      REAL(8)   ,ALLOCATABLE :: USMALL(:,:,:,:),UBIG(:,:,:,:)
!     **************************************************************************
!
!     ==========================================================================
!     == ALLOCATE UTENSOR AND CHOSE MATRIX ELEMENTS TO BE EVALUATED           ==
!     ==========================================================================
      NAT=SIZE(ISPECIES)
      ALLOCATE(UTENSOR(NAT))
      NN=0
      DO IAT=1,NAT
        NN=NN+1
        ISP=ISPECIES(IAT)
        N1=SUM(2*LOX(:LNX(ISP),ISP)+1)
        UTENSOR(NN)%IAT1=IAT
        UTENSOR(NN)%IAT2=IAT
        UTENSOR(NN)%IAT3=IAT
        UTENSOR(NN)%IAT4=IAT
        UTENSOR(NN)%IT2=(/0,0,0/)
        UTENSOR(NN)%IT3=(/0,0,0/)
        UTENSOR(NN)%IT4=(/0,0,0/)
        UTENSOR(NN)%N1=N1
        UTENSOR(NN)%N2=N1
        UTENSOR(NN)%N3=N1
        UTENSOR(NN)%N4=N1
        ALLOCATE(UTENSOR(NN)%U(N1,N1,N1,N1))
      ENDDO
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      CALL CELL$GETR8A('T0',9,RBAS)
!
!     ==========================================================================
!     == ALLOCATE HELPER ARRAYS                                               ==
!     ==========================================================================
      NE1=MAXVAL(GAUSSORB(:)%NE)
      ALLOCATE(E1(NE1))
      ALLOCATE(E2(NE1))
      ALLOCATE(E3(NE1))
      ALLOCATE(E4(NE1))
!
!     ==========================================================================
!     == DETERMINE THE GAUSSIAN CONTRIBUTION TO THE OVERLAP                   ==
!     ==========================================================================
      DO NN=1,SIZE(UTENSOR)
        IAT1=UTENSOR(NN)%IAT1
        IAT2=UTENSOR(NN)%IAT2
        IAT3=UTENSOR(NN)%IAT3
        IAT4=UTENSOR(NN)%IAT4
        NIJK1=GAUSSORB(IAT1)%NIJK
        NIJK2=GAUSSORB(IAT2)%NIJK
        NIJK3=GAUSSORB(IAT3)%NIJK
        NIJK4=GAUSSORB(IAT4)%NIJK
        NORB1=GAUSSORB(IAT1)%NORB
        NORB2=GAUSSORB(IAT2)%NORB
        NORB3=GAUSSORB(IAT3)%NORB
        NORB4=GAUSSORB(IAT4)%NORB
        NE1=GAUSSORB(IAT1)%NE
        NE2=GAUSSORB(IAT2)%NE
        NE3=GAUSSORB(IAT3)%NE
        NE4=GAUSSORB(IAT4)%NE
        E1(:NE1)=GAUSSORB(IAT1)%E
        E2(:NE2)=GAUSSORB(IAT2)%E
        E3(:NE3)=GAUSSORB(IAT3)%E
        E4(:NE4)=GAUSSORB(IAT4)%E
        R1=R0(:,IAT1)
        R2=R0(:,IAT2)+MATMUL(RBAS,REAL(UTENSOR(NN)%IT2))
        R3=R0(:,IAT3)+MATMUL(RBAS,REAL(UTENSOR(NN)%IT3))
        R4=R0(:,IAT4)+MATMUL(RBAS,REAL(UTENSOR(NN)%IT4))
        CALL GAUSSIAN$FOURCENTER(NIJK1,NE1,NORB1,E1(:NE1),R1,GAUSSORB(IAT1)%C &
                                ,NIJK2,NE2,NORB2,E2(:NE2),R2,GAUSSORB(IAT2)%C &
                                ,NIJK3,NE3,NORB3,E3(:NE3),R3,GAUSSORB(IAT3)%C &
                                ,NIJK4,NE4,NORB4,E4(:NE4),R4,GAUSSORB(IAT4)%C &
     &                          ,UTENSOR(NN)%U)
PRINT*,'============ UTENSOR ======================'
        WRITE(*,*)NN,(UTENSOR(NN)%U(I,I,I,I),I=1,UTENSOR(NN)%N1)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_OLDFOURCENTERGAUSS()
!     **************************************************************************
!     **  EVALUATES THE U-TENSOR                                              **
!     **                                                                      **
!     **  PRELIMINARY VERSION: ONLY ONSITE MATRIX ELEMENTS                    **
!     **  PRELIMINARY VERSION: NO AUGMENTATION CONTRIBUTIONS INCLUDED         **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : UTENSOR,ISPECIES,LNX,LOX,GAUSSORB,SBARLI1
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      REAL(8)                :: RBAS(3,3)
      INTEGER(4)             :: N1,N2,N3,N4
      INTEGER(4)             :: IAT1,IAT2,IAT3,IAT4
      INTEGER(4)             :: NIJK1,NIJK2,NIJK3,NIJK4
      INTEGER(4)             :: NORB1,NORB2,NORB3,NORB4
      INTEGER(4)             :: NE1,NE2,NE3,NE4
      REAL(8)                :: R1(3),R2(3),R3(3),R4(3)
      REAL(8)   ,ALLOCATABLE :: E1(:),E2(:),E3(:),E4(:)
      INTEGER(4)             :: IAT,NN,ISP,IND,LN,LM,L,IM
      REAL(8)   ,ALLOCATABLE :: USMALL(:,:,:,:),UBIG(:,:,:,:)
!     **************************************************************************
!
!     ==========================================================================
!     == ALLOCATE UTENSOR AND CHOSE MATRIX ELEMENTS TO BE EVALUATED           ==
!     ==========================================================================
      ALLOCATE(UTENSOR(NAT))
      NN=0
      DO IAT=1,NAT
        NN=NN+1
        ISP=ISPECIES(IAT)
        N1=SUM(2*LOX(:LNX(ISP),ISP)+1)
        UTENSOR(NN)%IAT1=IAT
        UTENSOR(NN)%IAT2=IAT
        UTENSOR(NN)%IAT3=IAT
        UTENSOR(NN)%IAT4=IAT
        UTENSOR(NN)%IT2=(/0,0,0/)
        UTENSOR(NN)%IT3=(/0,0,0/)
        UTENSOR(NN)%IT4=(/0,0,0/)
        UTENSOR(NN)%N1=N1
        UTENSOR(NN)%N2=N1
        UTENSOR(NN)%N3=N1
        UTENSOR(NN)%N4=N1
        ALLOCATE(UTENSOR(NN)%U(N1,N1,N1,N1))
      ENDDO
!
!     ==========================================================================
!     == COLLECT ATOMIC STRUCTURE                                             ==
!     ==========================================================================
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      CALL CELL$GETR8A('T0',9,RBAS)
!
!     ==========================================================================
!     == ALLOCATE HELPER ARRAYS                                               ==
!     ==========================================================================
      NE1=MAXVAL(GAUSSORB(:)%NE)
      ALLOCATE(E1(NE1))
      ALLOCATE(E2(NE1))
      ALLOCATE(E3(NE1))
      ALLOCATE(E4(NE1))
      N1=MAXVAL(UTENSOR(:)%N1)
      N2=MAXVAL(UTENSOR(:)%N2)
      N3=MAXVAL(UTENSOR(:)%N3)
      N4=MAXVAL(UTENSOR(:)%N4)
      ALLOCATE(USMALL(N1,N2,N3,N4))
      ALLOCATE(UBIG(N1,N2,N3,N4))
!
!     ==========================================================================
!     == DETERMINE THE GAUSSIAN CONTRIBUTION TO THE OVERLAP                   ==
!     ==========================================================================
      DO NN=1,SIZE(UTENSOR)
        IAT1=UTENSOR(NN)%IAT1
        IAT2=UTENSOR(NN)%IAT2
        IAT3=UTENSOR(NN)%IAT3
        IAT4=UTENSOR(NN)%IAT4
        NIJK1=GAUSSORB(IAT1)%NIJK
        NIJK2=GAUSSORB(IAT2)%NIJK
        NIJK3=GAUSSORB(IAT3)%NIJK
        NIJK4=GAUSSORB(IAT4)%NIJK
        NORB1=GAUSSORB(IAT1)%NORB
        NORB2=GAUSSORB(IAT2)%NORB
        NORB3=GAUSSORB(IAT3)%NORB
        NORB4=GAUSSORB(IAT4)%NORB
        NE1=GAUSSORB(IAT1)%NE
        NE2=GAUSSORB(IAT2)%NE
        NE3=GAUSSORB(IAT3)%NE
        NE4=GAUSSORB(IAT4)%NE
        E1(:NE1)=GAUSSORB(IAT1)%E
        E2(:NE2)=GAUSSORB(IAT2)%E
        E3(:NE3)=GAUSSORB(IAT3)%E
        E4(:NE4)=GAUSSORB(IAT4)%E
        R1=R0(:,IAT1)
        R2=R0(:,IAT2)+MATMUL(RBAS,REAL(UTENSOR(NN)%IT2))
        R3=R0(:,IAT3)+MATMUL(RBAS,REAL(UTENSOR(NN)%IT3))
        R4=R0(:,IAT4)+MATMUL(RBAS,REAL(UTENSOR(NN)%IT4))
        CALL GAUSSIAN$FOURCENTER(NIJK1,NE1,NORB1,E1(:NE1),R1,GAUSSORB(IAT1)%C &
                                ,NIJK2,NE2,NORB2,E2(:NE2),R2,GAUSSORB(IAT2)%C &
                                ,NIJK3,NE3,NORB3,E3(:NE3),R3,GAUSSORB(IAT3)%C &
                                ,NIJK4,NE4,NORB4,E4(:NE4),R4,GAUSSORB(IAT4)%C &
     &                          ,USMALL(:NORB1,:NORB2,:NORB3,:NORB4))
!
!       ========================================================================
!       == EXPAND U-TENSOR TO FULL SIZE                                       ==
!       ========================================================================
        N1=UTENSOR(NN)%N1
        N2=UTENSOR(NN)%N2
        N3=UTENSOR(NN)%N3
        N4=UTENSOR(NN)%N4
!
        ISP=ISPECIES(IAT1)
        IND=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          LM=SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            LM=LM+1
            IND=IND+1
            UBIG(IND,:NORB2,:NORB3,:NORB4)=USMALL(LM,:NORB2,:NORB3,:NORB4)
          ENDDO
        ENDDO
        ISP=ISPECIES(IAT2)
        IND=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          LM=SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            LM=LM+1
            IND=IND+1
             USMALL(:N1,IND,:NORB3,:NORB4)=UBIG(:N1,LM,:NORB3,:NORB4)
          ENDDO
        ENDDO
        ISP=ISPECIES(IAT3)
        IND=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          LM=SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            LM=LM+1
            IND=IND+1
             UBIG(:N1,:N2,IND,:NORB4)=USMALL(:N1,:N2,LM,:NORB4)
          ENDDO
        ENDDO
        ISP=ISPECIES(IAT4)
        IND=0
        DO LN=1,LNX(ISP)
          L=LOX(LN,ISP)
          LM=SBARLI1(L+1,ISP)-1
          DO IM=1,2*L+1
            LM=LM+1
            IND=IND+1
            USMALL(:N1,:N2,:N3,IND)=UBIG(:N1,:N2,:N3,LM)
          ENDDO
        ENDDO
        UTENSOR(NN)%U(:,:,:,:)=USMALL(:N1,:N2,:N3,:N4)
      ENDDO

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_UNNETWORK()
!     **************************************************************************
!     ** SET UP NETWORK OF U-TENSOR ELEMENTS OF SCREENED ORBITALS             **
!     ** A U-TENSER CORRESPONDS TO TWO PAIRS OF NEIGHBORS                     **
!     ** ONE ELEMENT IS IDENTIFIED BY TWO BOND INDICES NN1,NN2 REFERRING      **
!     ** TO AN ATOM PAIR OF SBAR AND A TRANSLATION FOR THE SECOND PAIR        **
!     ** ONLY THOSE U-TENSORS WITH ONE COMMON AT LEAST ONE COMMON ATOM        **
!     ** ARE INCLUDED.                                                        **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : SBAR,ISPECIES,NNUX,NNU,UMAT
      IMPLICIT NONE
      INTEGER(4)             :: NNS
      INTEGER(4)             :: NAT
      INTEGER(4),ALLOCATABLE :: NNONSITE(:) !(NAT)
      INTEGER(4)             :: NN1,NN2
      INTEGER(4)             :: IORIENT
      INTEGER(4)             :: IAT1,IAT2,IAT1B,IAT2B,IAT
      INTEGER(4)             :: NN
      LOGICAL(4)             :: TONSITE1,TONSITE2
!     **************************************************************************
      NNS=SIZE(SBAR)
      NAT=SIZE(ISPECIES)
      NNUX=NAT*100
      ALLOCATE(UMAT(NNUX))
!
!     ==========================================================================
!     == SET UP NETWORK OF U-TENSOR ELEMENTS OF SCREENED ORBITALS             ==
!     == A U-TENSER CORRESPONDS TO TWO PAIRS OF NEIGHBORS                     ==
!     == ONE ELEMENT IS IDENITIFED BY TWO BOND INDICES AND A TRANSLATION      ==
!     == ONLY THOSE U-TENSORS WITH ONE COMMON AT LEAST ONE COMMON ATOM        ==
!     == ARE INCLUDED.                                                        ==
!     ==========================================================================
      NNS=SIZE(SBAR)
      ALLOCATE(NNONSITE(NAT))
!     == IDENTIFY ONSITE TERMS FOR EACH ATOM ===================================
      NNONSITE(:)=0
      DO NN=1,NNS
        CALL LMTO_BONDORIENTATION(SBAR(NN)%IAT1,SBAR(NN)%IAT2,SBAR(NN)%IT &
     &                           ,IORIENT)
        IF(IORIENT.EQ.0) NNONSITE(SBAR(NN)%IAT1)=NN
      ENDDO
!     == CHECK COMPLETENESS OF NNONSITE ARRAY  =================================
      DO IAT=1,NAT
        IF(NNONSITE(IAT).EQ.0) THEN
          CALL ERROR$MSG('MISSING ONSITE TERM. INTERNAL ERROR')
          CALL ERROR$STOP('LMTO$FOURCENTERGAUSS')
        END IF
      ENDDO
!     == SET UP LIST FOR U-TENSOR ==============================================
      NN=0
      DO NN1=1,NNS
        IAT1=SBAR(NN1)%IAT1
        IAT2=SBAR(NN1)%IAT2
        CALL LMTO_BONDORIENTATION(IAT1,IAT2,SBAR(NN1)%IT,IORIENT)
        IF(IORIENT.LT.0) CYCLE   ! REMOVE IDENTICAL BONDS
        TONSITE1=(IORIENT.EQ.0) 
        DO NN2=NN1,NNS           ! AVOID IDENTICAL PAIRS OF BONDS
          IAT1B=SBAR(NN2)%IAT1
          IAT1B=SBAR(NN2)%IAT2
          CALL LMTO_BONDORIENTATION(IAT1B,IAT2B,SBAR(NN2)%IT,IORIENT)
          IF(IORIENT.LT.0) CYCLE ! REMOVE IDENTICAL BONDS
          TONSITE2=(IORIENT.EQ.0)
!
          IF(IAT1.EQ.IAT1B) THEN
!           == IAT1(NN1)=IAT1(NN2) =============================================
            NN=NN+1
            IF(NN.GT.NNUX) THEN
              CALL ERROR$MSG('#(MATRIX ELEMENTS OUT OF RANGE')
              CALL ERROR$STOP('LMTO_UNNETWORK')
            END IF
            UMAT(NN)%NN1=NN1
            UMAT(NN)%NN2=NN2
            UMAT(NN)%IT=(/0,0,0/)
          END IF
          IF(.NOT.TONSITE1.AND.IAT2.EQ.IAT1B) THEN
!           == IAT2(NN1)=IAT1(NN2) =============================================
            NN=NN+1
            IF(NN.GT.NNUX) THEN
              CALL ERROR$MSG('#(MATRIX ELEMENTS OUT OF RANGE')
              CALL ERROR$STOP('LMTO_UNNETWORK')
            END IF
            UMAT(NN)%NN1=NN1
            UMAT(NN)%NN2=NN2
            UMAT(NN)%IT=SBAR(NN1)%IT
          END IF
          IF(.NOT.TONSITE2.AND.IAT2.EQ.IAT1B) THEN
!           == IAT1(NN1)=IAT2(NN2) =============================================
            NN=NN+1
            IF(NN.GT.NNUX) THEN
              CALL ERROR$MSG('#(MATRIX ELEMENTS OUT OF RANGE')
              CALL ERROR$STOP('LMTO_UNNETWORK')
            END IF
            UMAT(NN)%NN1=NN1
            UMAT(NN)%NN2=NN2
            UMAT(NN)%IT=-SBAR(NN2)%IT
          END IF
          IF(.NOT.(TONSITE2.OR.TONSITE2).AND.IAT2.EQ.IAT1B) THEN
!           == IAT1(NN1)=IAT2(NN2) =============================================
            NN=NN+1
            IF(NN.GT.NNUX) THEN
              CALL ERROR$MSG('#(MATRIX ELEMENTS OUT OF RANGE')
              CALL ERROR$STOP('LMTO_UNNETWORK')
            END IF
            UMAT(NN)%NN1=NN1
            UMAT(NN)%NN2=NN2
            UMAT(NN)%IT=SBAR(NN1)%IT-SBAR(NN2)%IT
          END IF
        ENDDO
!       == INCLUDE TWO ONSITE TERMS AT NEIGHBORING SITES =======================
        IF(.NOT.TONSITE1) THEN
          NN=NN+1
            IF(NN.GT.NNUX) THEN
              CALL ERROR$MSG('#(MATRIX ELEMENTS OUT OF RANGE')
              CALL ERROR$STOP('LMTO_UNNETWORK')
            END IF
          UMAT(NN)%NN1=NNONSITE(IAT1)
          UMAT(NN)%NN2=NNONSITE(IAT2)
          UMAT(NN)%IT=SBAR(NN1)%IT
        END IF
      ENDDO
      NNU=NN
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
      USE LMTO_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: NAT
      INTEGER(4)             :: LX1
      INTEGER(4)             :: IPOS,IAT,ISP,L,LN
!     **************************************************************************
      IF(ALLOCATED(SBARATOMI1)) RETURN
                              CALL TRACE$PUSH('LMTO$SBARINDICES')
!
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
                              CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_SBARINDICES2()
!     **************************************************************************
!     **  CREATES AN INDEX ARRAY TO BE USED WITH THE SCREENED STRUCTURE       **
!     **  CONSTANTS.                                                          **
!     **     SBARLI1(L+1,ISP) POINTS TO THE FIRST OF 2*L+1 ENTRIES FOR ATOM   **
!     **                     TYPE ISP AND MAIN ANGULAR MOMENTUM L             **
!     **                     (RELATIVE TO THE FIRST ELEMENT FOR THIS ATOM)    **
!     **                                                                      **
!     **  THE ARRAY LX1 HAVE STRANGE NAMES BECAUSE THE ORIGINAL               **
!     **  NAMES ARE ALREADY USED BY LMTO_MODULE. WE USE SEPARATE ARRAYS,      **
!     **  BECAUSE WE WANT TO REMOVE THESE ARRAYS FROM THE MODULE              **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: LX1
      INTEGER(4)             :: IPOS,ISP,L,LN
!     **************************************************************************
      IF(ALLOCATED(SBARLI1)) RETURN
                              CALL TRACE$PUSH('LMTO$SBARINDICES2')
!
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
!!$WRITE(*,FMT='("A=",20F10.5)')A
!!$WRITE(*,FMT='("B=",20F10.5)')B
!!$WRITE(*,FMT='("C=",20F10.5)')C
!!$WRITE(*,FMT='("D=",20F10.5)')D
!!$WRITE(*,FMT='("E=",20F10.5)')E
!!$WRITE(*,FMT='("F=",20F10.5)')F
      CALL LMTO_PREPARE2(XK,NRL,C,E,F,G,H)
!!$DO I=1,NRL
!!$  WRITE(*,FMT='("G=",20F10.5)')G(I,:)
!!$ENDDO
!!$DO I=1,NRL
!!$  WRITE(*,FMT='("H=",20F10.5)')H(I,:)
!!$ENDDO
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
!!$DO I=1,NBH
!!$  WRITE(*,FMT='("PROJCONTR1=",50F10.5)')REAL(PROJCONTR1(1,I,:))
!!$  WRITE(*,FMT='("PROJCONTR1=",50F10.5)')AIMAG(PROJCONTR1(1,I,:))
!!$ENDDO
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
      USE LMTO_MODULE
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
!!$DO I=1,NRL
!!$  WRITE(*,FMT='("SBAR=",20F10.5)')SBAR(I,:)
!!$ENDDO
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
      USE LMTO_MODULE
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
      SUBROUTINE LMTO$ORBRAD
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : K2,SBAR,ORBRAD
      USE PERIODICTABLE_MODULE
      REAL(8)   ,ALLOCATABLE      :: SBARDIAG(:)
      REAL(8)                     :: RAD1
      REAL(8)                     :: VALPHI,DERPHI
      REAL(8)                     :: VALK,DERK
      REAL(8)                     :: VALJ,DERJ
      REAL(8)                     :: AEZ
      INTEGER(4)                  :: LX1
      INTEGER(4)                  :: LNX
      INTEGER(4),ALLOCATABLE      :: LX(:)
      INTEGER(4),ALLOCATABLE      :: LOX(:)
      INTEGER(4)                  :: NNS
      INTEGER(4)                  :: NAT
      INTEGER(4)                  :: LMX
      INTEGER(4)                  :: IAT,ISP,INB,LM,L,IM
      INTEGER(4)                  :: LXX !X(ANGULAR MOMENTUM)
      INTEGER(4)                  :: NSP !#(ATOM TYPES)
      REAL(8)    ,ALLOCATABLE     :: QBAR(:,:)
      INTEGER(4) ,ALLOCATABLE     :: ISPECIES(:)
!     **************************************************************************
!
!     ==========================================================================
!     == DETERMINE SCREENING PARAMETERS QBAR                                  ==
!     ==========================================================================
      CALL SETUP$GETI4('NSP',NSP)
      ALLOCATE(LX(NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        LX(ISP)=MAXVAL(LOX)
        DEALLOCATE(LOX)
        CALL SETUP$ISELECT(0)
      ENDDO
      LXX=MAXVAL(LX)
      ALLOCATE(QBAR(LXX+1,NSP))
      CALL LMTO_MAKEQBAR(NSP,LXX,K2,QBAR)

!     ==========================================================================
!     == COLLECT DATA AND ALLOCATE ORBRAD                                     ==
!     ==========================================================================
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
!
      IF(.NOT.ALLOCATED(ORBRAD))ALLOCATE(ORBRAD(LXX+1,NAT))
!
!     ==========================================================================
!     == DETERMINE SCREENING PARAMETERS QBAR                                  ==
!     ==========================================================================
      NNS=SIZE(SBAR)
      DO IAT=1,NAT
        ISP=ISPECIES(IAT)
        LX1=LX(ISP)
        LMX=(LX1+1)**2
        ALLOCATE(SBARDIAG(LX1+1))
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(NINT(AEZ),'R(ASA)',RAD1)
!
!       ========================================================================
!       == FIND ONSITE STRUCTURE CONSTANTS                                    ==
!       ========================================================================
        DO INB=1,NNS
          IF(SBAR(INB)%IAT1.NE.IAT) CYCLE
          IF(SBAR(INB)%IAT2.NE.IAT) CYCLE
          IF(SBAR(INB)%IT(1).NE.0) CYCLE
          IF(SBAR(INB)%IT(2).NE.0) CYCLE
          IF(SBAR(INB)%IT(3).NE.0) CYCLE
          IF(SBAR(INB)%N1.NE.LMX) THEN
            CALL ERROR$MSG('DIMENSIONS INCONSISTENT')
            CALL ERROR$STOP('LMTO$ORBRAD')
          END IF
          LM=0
          DO L=0,LX1
            SBARDIAG(L+1)=0.D0
            DO IM=1,2*L+1
              LM=LM+1
              SBARDIAG(L+1)=SBARDIAG(L+1)+SBAR(INB)%MAT(LM,LM)
            ENDDO
            SBARDIAG(L+1)=SBARDIAG(L+1)/REAL(2*L+1,KIND=8)
          ENDDO
          EXIT
        ENDDO
!
!       ========================================================================
!       == DETERMINE RADIUS BY LINEAR EXTRAPOLATION                           ==
!       ========================================================================
        LM=0
        DO L=0,LX1
          CALL LMTO$SOLIDBESSELRAD(L,RAD1,K2,VALJ,DERJ)
          CALL LMTO$SOLIDHANKELRAD(L,RAD1,K2,VALK,DERK)
!         == SCREEN BESSEL FUNCTIONS |JBAR>=|J>-|K>QBAR ========================
          VALJ=VALJ-VALK*QBAR(L+1,ISP)
          DERJ=DERJ-DERK*QBAR(L+1,ISP)
!         ==  |PHI>   <-   |K>-|JBAR>*SBAR =====================================
          VALPHI=VALK-VALJ*SBARDIAG(L+1)
          DERPHI=DERK-DERJ*SBARDIAG(L+1)
!         == DETERMINE APPROXIMATE POSITION OF THE NODE ========================
          ORBRAD(L+1,IAT)=RAD1-VALPHI/DERPHI
          WRITE(6,FMT='("IAT=",I3," Z=",I3," L=",I2," R[ASA]=",F10.5," R[ORB]=",F10.5)') &
     &                IAT,NINT(AEZ),L,RAD1,RAD1-VALPHI/DERPHI
        ENDDO
        DEALLOCATE(SBARDIAG)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_MAKEQBAR(NSP,LXX,K2,QBAR)
!     **************************************************************************
!     ** PARAMETER NEEDED TO  SCREEN THE STRUCTURE CONSTANTS                  **
!     **   |K>-|J>QBAR HAS THE SAME LOGARITHMIC DERIVATIVE AS |PHIDOT>.       **
!     ** VAL AND DER ARE VALUE AND DERIVATIVE OF PHIDOT.                      **
!     **                                                                      **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSP        !#(ATOM-TYPES)
      INTEGER(4),INTENT(IN) :: LXX        ! X(ANGULAR MOMENTUM)
      REAL(8)   ,INTENT(IN) :: K2         ! 
      REAL(8)   ,INTENT(OUT):: QBAR(LXX+1,NSP)  !SCREENING PARAMETER
      INTEGER(4)            :: NSP1       !#(ATOM-TYPES)
      INTEGER(4)            :: LNX        !#(PARTIAL WAVES)
      INTEGER(4)            :: GID        !GRID ID
      INTEGER(4)            :: NR         ! #(RADIAL GRID POINTS)
      INTEGER(4),ALLOCATABLE:: LOX(:)     !(LNX)
      INTEGER(4),ALLOCATABLE:: ISCATT(:)  !(LNX)
      REAL(8)   ,ALLOCATABLE:: NLPHIDOT(:,:)  !(NR,LNX) NODELESS PARTIALWAVE
      REAL(8)               :: AEZ        !ATOMIC NUMBER
      REAL(8)               :: RAD        !ATOMIC RADIUS
      REAL(8)               :: VAL,DER    !VALUE AND DERIAVTIVE OF U
      REAL(8)               :: JVAL,JDER  !VALUE AND DERIAVTIVE OF J
      REAL(8)               :: KVAL,KDER  !VALUE AND DERIAVTIVE OF J
      INTEGER(4)            :: ISP,LN,L
!     **************************************************************************
      CALL SETUP$NSPECIES(NSP1)
      IF(NSP1.NE.NSP) THEN
        CALL ERROR$MSG('INCONSISTENT VALUES')
        CALL ERROR$MSG('NSP ON INPUT DIFFERS FROM THAT OF THE SETUP OBJECT')
        CALL ERROR$STOP('LMTO_MAKEQBAR')
      END IF
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      QBAR(:,:)=0.D0
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(ISCATT(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(NLPHIDOT(NR,LNX))
        CALL SETUP$GETR8A('NLPHIDOT',NR*LNX,NLPHIDOT)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(NINT(AEZ),'R(ASA)',RAD)
        DO L=0,LXX
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            IF(ISCATT(LN).GT.0) CYCLE
            CALL RADIAL$VALUE(GID,NR,NLPHIDOT(:,LN),RAD,VAL)
            CALL RADIAL$DERIVATIVE(GID,NR,NLPHIDOT(:,LN),RAD,DER)
            CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
            CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
            QBAR(L+1,ISP)=(JVAL*DER-VAL*JDER)/(KVAL*DER-VAL*KDER)
          ENDDO
        ENDDO            
        DEALLOCATE(NLPHIDOT)
        DEALLOCATE(LOX)
        DEALLOCATE(ISCATT)
      ENDDO
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE LMTO$OVERLAP(NFIL)
!!$!     **************************************************************************      
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$DOLOCORB_2(IAT,ISP,GID,NR,LNXCHI,LNXPHI,TORB,CHIPHI,CHI)
!     **************************************************************************
!     ** new version!!!!!
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
      USE LMTO_MODULE, ONLY : SBAR,POTPAR,sbarli1,k2
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
      real(8)   ,allocatable :: sbarav(:)
      REAL(8)   ,ALLOCATABLE :: AMAT(:,:)   !(LNXCHI1,LNXPHI) 
      REAL(8)   ,ALLOCATABLE :: BMAT(:,:)   !(LNXCHI1,LNXCHI1) 
      REAL(8)   ,ALLOCATABLE :: XMAT(:,:)   !(LNXPHI,LNXCHI) 
      REAL(8)                :: RAD              ! COVALENT RADIUS
      REAL(8)                :: AEZ               ! ATOMIC NUMBER
      REAL(8)                :: AUX(NR)
      REAL(8)                :: R(NR)
      REAL(8)                :: SVAR,SVAR1,svar2,VAL,DER
      REAL(8)                :: kval,kder,jval,jder
      REAL(8)                :: qbar
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: LMX
      INTEGER(4)             :: LN,ln1,L,I,J,IIB,LM,LNCHI,IR,im
      INTEGER(4)             :: Lx
      INTEGER(4)             :: iorb
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
        IF(TORB(ln))LNCHI=LNCHI+1
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
      CALL SETUP$GETR8A('PRO',NR*LNX,PRO)
      lx=MAXVAL(LOX(:))
!
!     ==========================================================================
!     == FIND ONSITE STRUCTURE CONSTANTS                                      ==
!     ==========================================================================
      allocate(sbarav(lx+1))      
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
print*,'sbarav in dolocorb_2: ',iat,isp,sbarav
!     ==========================================================================
!     == CONSTRUCT LOCAL ORBITALS                                             ==
!     ==========================================================================
      allocate(aechi(nr,lnchi))
      allocate(pschi(nr,lnchi))
      allocate(nlchi(nr,lnchi))
      LNCHI=lnxphi
      DO LN=1,LNXphi
        L=LOX(LN)
        ln1=potpar(isp)%tailed%lndot(ln)
        aechi(:,ln)=potpar(isp)%tailed%aef(:,ln) &
    &              -potpar(isp)%tailed%aef(:,ln1)*sbarav(l+1)
        pschi(:,ln)=potpar(isp)%tailed%psf(:,ln) &
    &              -potpar(isp)%tailed%psf(:,ln1)*sbarav(l+1)
        nlchi(:,ln)=potpar(isp)%tailed%nlf(:,ln) &
    &              -potpar(isp)%tailed%nlf(:,ln1)*sbarav(l+1)
      ENDDO
WRITE(STRING,FMT='(I5)')IAT
STRING='LOCORB_FORATOM'//TRIM(ADJUSTL(STRING))//'.DAT'
CALL LMTO_WRITEPHI(TRIM(STRING),GID,NR,LNxpHI,AECHI)
!
!     ==ORTHONORMALIZE LOCAL ORBITALS ==========================================
!     == orthonormalization is not required and serves only estaetical purposes
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
        CHIPHI(LNCHI,:)=AMAT(LN,:)   ! matching coefficients
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
      USE LMTO_MODULE, ONLY : SBAR,POTPAR,sbarli1,k2
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
      real(8)   ,allocatable :: sbarav(:)
      REAL(8)   ,ALLOCATABLE :: AMAT(:,:)   !(LNXCHI1,LNXPHI) 
      REAL(8)   ,ALLOCATABLE :: BMAT(:,:)   !(LNXCHI1,LNXCHI1) 
      REAL(8)   ,ALLOCATABLE :: XMAT(:,:)   !(LNXPHI,LNXCHI) 
      REAL(8)                :: RAD              ! COVALENT RADIUS
      REAL(8)                :: AEZ               ! ATOMIC NUMBER
      REAL(8)                :: AUX(NR)
      REAL(8)                :: R(NR)
      REAL(8)                :: SVAR,SVAR1,svar2,VAL,DER
      REAL(8)                :: kval,kder,jval,jder
      REAL(8)                :: qbar
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: LMX
      INTEGER(4)             :: LN,ln1,L,I,J,IIB,LM,LNCHI,IR,im
      INTEGER(4)             :: Lx
      INTEGER(4)             :: iorb
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
        IF(TORB(ln))LNCHI=LNCHI+1
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
!     == select the correct phidot function (only one per l) ===================
      DO LN=1,LNX
        LN1=POTPAR(ISP)%LNSCATT(LN)
        IF(LN1.eq.LN) CYCLE
        AEPHIDOT(:,LN)=AEPHIDOT(:,LN1)
        PSPHIDOT(:,LN)=PSPHIDOT(:,LN1)
        NLPHIDOT(:,LN)=NLPHIDOT(:,LN1)
      ENDDO
!
      lx=MAXVAL(LOX(:))
!
!     ==========================================================================
!     == FIND ONSITE STRUCTURE CONSTANTS                                      ==
!     ==========================================================================
      allocate(sbarav(lx+1))      
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
      ALLOCATE(AECHI(NR,LNXphi))
      ALLOCATE(PSCHI(NR,LNXphi))
      ALLOCATE(NLCHI(NR,LNXphi))
!
!     ==========================================================================
!     == CONSTRUCT LOCAL ORBITALS                                             ==
!     ==========================================================================
      LNCHI=0
      DO LN=1,LNXphi
        L=LOX(LN)
        SVAR1=POTPAR(ISP)%KTOPHI(LN)
        SVAR2=POTPAR(ISP)%KTOPHIDOT(LN) &
    &        -POTPAR(ISP)%JBARTOPHIDOT(LN)*SBARAV(L+1)
        AECHI(:,LN)=AEPHI(:,LN)*SVAR1+AEPHIDOT(:,LN)*SVAR2
        PSCHI(:,LN)=PSPHI(:,LN)*SVAR1+PSPHIDOT(:,LN)*SVAR2
        NLCHI(:,LN)=NLPHI(:,LN)*SVAR1+NLPHIDOT(:,LN)*SVAR2
!
!       == ATTACH EXPONENTIAL TAIL AT THE matching radius ======================
        CALL LMTO$SOLIDHANKELRAD(L,Rad,K2,Kval,kder)
        CALL LMTO$SOLIDBESSELRAD(L,Rad,K2,jval,jder)
        val=kval-(jval-kval*potpar(isp)%qbar(ln))*sbarav(l+1)
        der=kder-(jder-kder*potpar(isp)%qbar(ln))*sbarav(l+1)
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
!     == orthonormalization is not required and serves only estaetical purposes
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
        CHIPHI(LNCHI,:)=AMAT(LN,:)   ! matching coefficients
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
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE LMTO_CHIFROMPHI(GID,NR,LXX,RASA,RAD &
!!$     &                           ,LNXPHI,LOXPHI,ISCATT,PHI &
!!$     &                           ,LNXCHI,LOXCHI,AMAT,BMAT)
!!$!     **************************************************************************
!!$!     **  DEFINES PROJECTOR FUNCTIONS FOR LOCAL ORBITALS IN TERMS OF          **
!!$!     **  THE CONVENTIONAL PARTIAL WAVE PROJECTOR FUNCTIONS                   **
!!$!     **                                                                      **
!!$!     **     <P_CHI_I|= SUM_J BMAT(I,J) <P_PHI_J|                             **
!!$!     **     |CHI_I>  = SUM_J |PHI_J> AMAT(J,I)                               **
!!$!     **                                                                      **
!!$!     **  THE LOCAL ORBITALS ARE SELECTED BY NORB WHICH SPECIFIES THE NUMBER  **
!!$!     **  OF LOCAL ORBITALS PER ANGULAR MOMENTUM CHANNEL.                     **
!!$!     **                                                                      **
!!$!     **************************************************************************
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: GID     ! GRID ID
!!$      INTEGER(4),INTENT(IN) :: NR      ! #(RADIAL GRID POINTS)
!!$      INTEGER(4),INTENT(IN) :: LXX     ! X(ANGULAR MOMENTUM)
!!$      REAL(8)   ,INTENT(IN) :: RASA    ! ATOMIC RADIUS             
!!$      REAL(8)   ,INTENT(IN) :: RAD(LXX+1)   ! EXTENT OF HEAD FUNCTIONS
!!$      INTEGER(4),INTENT(IN) :: LNXPHI  ! #(PARTIAL WAVES)
!!$      INTEGER(4),INTENT(IN) :: LOXPHI(LNXPHI)! ANGULAR MOMENTUM PER PARTIAL WAVE
!!$      INTEGER(4),INTENT(IN) :: ISCATT(LNXPHI)! RELATIVE TO HIGHEST VALENCE STATE
!!$      REAL(8)   ,INTENT(IN) :: PHI(NR,LNXPHI)! PARTIAL WAVES
!!$      INTEGER(4),INTENT(IN) :: LNXCHI        ! #(LOCAL ORBITALS)
!!$      INTEGER(4),INTENT(OUT):: LOXCHI(LNXCHI)  ! ANGULAR MOMENTUM PER LOCAL ORB.
!!$      REAL(8)   ,INTENT(OUT):: BMAT(LNXCHI,LNXPHI) ! MAPPING OF PROJECTORS
!!$      REAL(8)   ,INTENT(OUT):: AMAT(LNXPHI,LNXCHI)     ! MAPPING OF WAVE FUNCTIONS
!!$      REAL(8)   ,ALLOCATABLE:: CHI(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: CHI1(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: A(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: B(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: A1(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: MAT(:,:)
!!$      REAL(8)   ,ALLOCATABLE:: MATINV(:,:)
!!$      REAL(8)               :: R(NR)       
!!$      REAL(8)   ,ALLOCATABLE:: G(:)        
!!$      REAL(8)   ,PARAMETER  :: RCG=1.D-2
!!$      REAL(8)   ,ALLOCATABLE:: AUX(:),AUX1(:)
!!$      REAL(8)               :: SVAR1,SVAR2
!!$      INTEGER(4)            :: IR
!!$      INTEGER(4)            :: NX,N,LX,L,LN,LNCHI,NOFL,NOFL2,ISVAR
!!$      INTEGER(4)            :: N1,N2,LN1,LN2,L1,L2
!!$      CHARACTER(64)         :: FILE
!!$      INTEGER(4)            :: NORB(LXX+1)
!!$      REAL(8)               :: VAL1,DER1,VAL2,DER2
!!$      REAL(8)               :: DR
!!$      REAL(8)               :: RCUT(LNXPHI)
!!$!     **************************************************************************
!!$                            CALL TRACE$PUSH('LMTO_CHIFROMPHI')
!!$      CALL RADIAL$R(GID,NR,R)
!!$      LX=MAXVAL(LOXPHI)
!!$      NORB(:)=0
!!$      DO LN=1,LNXCHI
!!$CALL ERROR$MSG('LOXCHI USED BUT NOT DEFINED')
!!$CALL ERROR$STOP('LMTO_CHIFROMPHI')
!!$        L=LOXCHI(LN)
!!$        NORB(L+1)=NORB(L+1)+1
!!$      ENDDO
!!$!
!!$      IF(LX.GT.LXX) THEN
!!$        CALL ERROR$MSG('LX>LXX')
!!$        CALL ERROR$STOP('LMTO_CHIFROMPHI')
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  START WITH PARTIAL WAVES AS LOCAL ORBITALS                          ==
!!$!     == |CHI_I>=SUM_I |PHI_J>A(J,I)                                          ==
!!$!     ==========================================================================
!!$      ALLOCATE(CHI(NR,LNXPHI))     
!!$      CHI(:,:)=0.D0
!!$      ALLOCATE(A(LNXPHI,LNXPHI))        
!!$      A(:,:)=0.D0
!!$      DO LN=1,LNXPHI
!!$        L=LOXPHI(LN)
!!$        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
!!$        A(LN,LN)=1.D0
!!$        CHI(:,LN)=PHI(:,LN)
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     == MAKE HEAD FUNCTION ANTIBONDING WITH NODE AT RCUT ======================
!!$!     ==========================================================================
!!$      ALLOCATE(AUX(NR))
!!$      ALLOCATE(AUX1(NR))
!!$      DO L=0,LX
!!$        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
!!$        NOFL=0
!!$        DO LN=1,LNXPHI
!!$          IF(LOXPHI(LN).NE.L) CYCLE
!!$          NOFL=NOFL+1
!!$          IF(ISCATT(LN).GT.0) CYCLE ! CONSIDER ONLY HEAD FUNCTIONS
!!$!         == SEARCH FOR NEXT PARTIAL WAVE WITH THIS L =========================
!!$          LN1=0
!!$          DO LN2=1,LNXPHI
!!$            IF(LOXPHI(LN2).NE.L) CYCLE  
!!$!ALTERNATIVE 1
!!$            LN1=LN2
!!$            IF(ISCATT(LN2).EQ.1) EXIT
!!$!ALTERNATIVE 2
!!$!            IF(ISCATT(LN2).EQ.1) THEN
!!$!              LN1=LN2
!!$!              EXIT
!!$!            END IF
!!$!END ALTERNATIVES
!!$          ENDDO
!!$          IF(LN1.EQ.0) THEN
!!$PRINT*,'ISCATT ',ISCATT
!!$PRINT*,'LOXPHI ',LOXPHI
!!$            CALL ERROR$MSG('CAN NOT LOCALIZE')
!!$            CALL ERROR$MSG('NO TAIL FUNCTION AVAILABLE')
!!$            CALL ERROR$I4VAL('L',L)
!!$            CALL ERROR$I4VAL('LN',LN)
!!$            CALL ERROR$STOP('LMTO_CHIFROMPHI')
!!$          END IF
!!$!
!!$!         == NOW CONTINUE; LN1 POINTS TO THE NEXT PHI WITH THIS L ==============
!!$!         == IMPOSE NODE CONDITION==============================================
!!$          CALL RADIAL$VALUE(GID,NR,CHI(:,LN),RASA,VAL1)
!!$          CALL RADIAL$DERIVATIVE(GID,NR,CHI(:,LN),RASA,DER1)
!!$          CALL RADIAL$VALUE(GID,NR,CHI(:,LN1),RASA,VAL2)
!!$          CALL RADIAL$DERIVATIVE(GID,NR,CHI(:,LN1),RASA,DER2)
!!$          DR=RAD(L+1)-RASA
!!$          SVAR1=1.D0
!!$          SVAR2=-(VAL1+DR*DER1)/(VAL2+DR*DER2)
!!$          IF(VAL1.EQ.0.D0.AND.VAL2.EQ.0.D0) THEN
!!$            CALL ERROR$MSG('PARTIAL WAVES ARE TRUNCATED INSIDE OF RCUT')
!!$            CALL ERROR$MSG('THIS IS A FLAW OF THE IMPLEMENTATION')
!!$            CALL ERROR$MSG('CHOOSE SMALLER RCUT')
!!$            CALL ERROR$STOP('LDAPLUSU_CHIFROMPHI')
!!$          END IF
!!$!
!!$          CHI(:,LN)=CHI(:,LN)*SVAR1+CHI(:,LN1)*SVAR2
!!$          A(:,LN)=A(:,LN)*SVAR1+A(:,LN1)*SVAR2
!!$!
!!$!         == DETERMINE ACTUAL NODE OF THE ORBITAL ==============================
!!$          DO IR=1,NR
!!$            IF(R(IR).LT.RASA) CYCLE
!!$            IF(CHI(IR,LN)*CHI(IR-1,LN).GT.0.D0) CYCLE
!!$            RCUT(LN)=R(IR-1)-CHI(IR-1,LN)/(CHI(IR,LN)-CHI(IR-1,LN))*(R(IR)-R(IR-1))
!!$            EXIT
!!$          ENDDO          
!!$          CALL RADIAL$VALUE(GID,NR,CHI(:,LN),RCUT(LN),VAL1)
!!$          CALL RADIAL$DERIVATIVE(GID,NR,CHI(:,LN),RCUT(LN),DER1)
!!$          RCUT(LN)=RCUT(LN)-VAL1/DER1
!!$!
!!$!         == ORTHOGONALIZE TO THE LOWER HEAD FUNCTIONS =========================
!!$          DO LN1=1,LN-1
!!$            IF(LOXPHI(LN1).NE.L) CYCLE  
!!$            AUX(:)=CHI(:,LN)*CHI(:,LN1)*R(:)**2
!!$            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$            CALL RADIAL$VALUE(GID,NR,AUX1,MIN(RCUT(LN),RCUT(LN1)),SVAR1)
!!$            CHI(:,LN)=CHI(:,LN)-CHI(:,LN1)*SVAR1
!!$            A(:,LN)=A(:,LN)-A(:,LN1)*SVAR1
!!$          ENDDO
!!$!
!!$!         == NORMALIZE HEAD FUNCTION ===========================================
!!$          AUX(:)=CHI(:,LN)**2*R(:)**2
!!$          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
!!$          CALL RADIAL$VALUE(GID,NR,AUX1,RCUT(LN),SVAR1)
!!$          SVAR1=1.D0/SQRT(SVAR1)
!!$          CHI(:,LN)=CHI(:,LN)*SVAR1
!!$          A(:,LN)=A(:,LN)*SVAR1
!!$        ENDDO
!!$      END DO          
!!$      DEALLOCATE(AUX)
!!$      DEALLOCATE(AUX1)
!!$!
!!$!     ==========================================================================
!!$!     == NOW INVERT THE MATRIX A: B:=A^{-1}                                   ==
!!$!     == |CHI_I>=SUM_J |PHI_J>A(J,I)                                          ==
!!$!     == 1=|PHI><P|=|PHI>A B <P|=|CHI> ( B<P| )                               ==
!!$!     == <P_CHI_I|=SUM_J B(I,J)<P_PHI_J|                                      ==
!!$!     ==========================================================================
!!$      ALLOCATE(B(LNXPHI,LNXPHI))
!!$      B(:,:)=0.D0
!!$      LX=MAXVAL(LOXPHI)
!!$      DO L=0,LX
!!$        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
!!$!
!!$        NX=0
!!$        DO LN=1,LNXPHI
!!$          IF(LOXPHI(LN).EQ.L) NX=NX+1
!!$        ENDDO
!!$        IF(NX.EQ.0) CYCLE   
!!$
!!$        ALLOCATE(MAT(NX,NX))
!!$        ALLOCATE(MATINV(NX,NX))
!!$!        
!!$!       == MAP A INTO MATRIX MAT ===============================================
!!$        N1=0
!!$        DO LN1=1,LNXPHI
!!$          L1=LOXPHI(LN1)
!!$          IF(L1.NE.L) CYCLE
!!$          N1=N1+1
!!$          N2=0
!!$          DO LN2=1,LNXPHI
!!$            L2=LOXPHI(LN2)
!!$            IF(L2.NE.L) CYCLE
!!$            N2=N2+1
!!$            MAT(N2,N1)=A(LN2,LN1)
!!$          ENDDO
!!$        ENDDO
!!$!
!!$!       == INVERT MATRIX =======================================================
!!$        CALL LIB$INVERTR8(NX,MAT,MATINV)
!!$!
!!$!       == MAP MATINV INTO B ===================================================
!!$        N1=0
!!$        DO LN1=1,LNXPHI
!!$          L1=LOXPHI(LN1)
!!$          IF(L1.NE.L) CYCLE
!!$          N1=N1+1
!!$          N2=0
!!$          DO LN2=1,LNXPHI
!!$            L2=LOXPHI(LN2)
!!$            IF(L2.NE.L) CYCLE
!!$            N2=N2+1
!!$            B(LN1,LN2)=MATINV(N1,N2) 
!!$          ENDDO
!!$        ENDDO
!!$        DEALLOCATE(MAT)
!!$        DEALLOCATE(MATINV)
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     === REMOVE TAIL FUNCTIONS                                               ==
!!$!     ==========================================================================
!!$      LNCHI=0
!!$      DO L=0,LX
!!$        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
!!$        NOFL=0
!!$        DO LN=1,LNXPHI
!!$          IF(L.NE.LOXPHI(LN)) CYCLE
!!$          NOFL=NOFL+1
!!$          IF(NOFL.GT.NORB(L+1)) EXIT
!!$          LNCHI=LNCHI+1
!!$          IF(LNCHI.GT.LNXCHI) THEN
!!$            CALL ERROR$MSG('LNCHI OUT OF RANGE')
!!$            CALL ERROR$I4VAL('LNCHI',LNCHI)
!!$            CALL ERROR$I4VAL('LNXCHI',LNXCHI)
!!$            CALL ERROR$STOP('SETUP_CHIFROMPHI')
!!$          END IF
!!$          LOXCHI(LNCHI)=L
!!$          BMAT(LNCHI,:)=B(LN,:)
!!$          AMAT(:,LNCHI)=A(:,LN)
!!$        END DO
!!$      ENDDO
!!$      IF(LNCHI.NE.LNXCHI) THEN
!!$        CALL ERROR$MSG('LNCHI AND LNXCHI ARE INCONSISTENT')
!!$        CALL ERROR$I4VAL('LNCHI',LNCHI)
!!$        CALL ERROR$I4VAL('LNXCHI',LNXCHI)
!!$        CALL ERROR$STOP('SETUP_CHIFROMPHI')
!!$      END IF
!!$!
!!$!     ==========================================================================
!!$!     ==  CLEAN UP                                                            ==
!!$!     ==========================================================================
!!$      DEALLOCATE(CHI)
!!$      DEALLOCATE(A)
!!$      DEALLOCATE(B)
!!$                            CALL TRACE$POP()
!!$      RETURN
!!$      END
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
      SUBROUTINE LMTO_GRIDTAILS(NAT,R0,IAT1,NP,P,NORB,ORB)
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
      INTEGER(4)            :: LX
      INTEGER(4)            :: LMX
      INTEGER(4)            :: LOFN(NORB)
      INTEGER(4)            :: LMOFN(NORB)
      REAL(8)   ,ALLOCATABLE:: YLM(:)
      REAL(8)   ,ALLOCATABLE:: KVAL(:)
      REAL(8)   ,ALLOCATABLE:: JVAL(:)
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)               :: RX
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: LNX1
      INTEGER(4)            :: NN,L,I1,I2,IM,IP,LN,IR
      REAL(8)               :: DR(3),DRLEN
      REAL(8)               :: RAD
      REAL(8)               :: QBAR
      REAL(8)               :: JVAL1,JDER1
      REAL(8)   ,PARAMETER  :: LAMBDA=1.D0
      REAL(8)   ,ALLOCATABLE:: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE:: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE:: KPRIME(:,:)
      REAL(8)   ,ALLOCATABLE:: JBARPRIME(:,:)
!     **************************************************************************
!
!     ==========================================================================
!     == COLLECT SCREENED STRUCTURE CONSTANTS FOR ATOM IAT1 ====================
!     ==========================================================================
      NNS=SIZE(SBAR)  !SBAR STRUCCONS. (GLOBAL)
      DO NN=1,NNS
         IF(SBAR(NN)%IAT1.NE.IAT1) CYCLE
         IF(SBAR(NN)%IAT2.NE.IAT1) CYCLE
         IF(MAXVAL(ABS(SBAR(NN)%IT(:))).NE.0) CYCLE
         IF(NORB.NE.SBAR(NN)%N1) THEN
           CALL ERROR$MSG('INCONSISTENT NUMBER OF ORBITALS ')
           CALL ERROR$STOP('LMTO_GRIDTAILS')
         END IF
         SBARMAT(:,:)=SBAR(NN)%MAT
         EXIT
       ENDDO
!
!     ==========================================================================
!     == CONSTRUCT KBAR AND JBAR WITH TAILS                                   ==
!     ==========================================================================
      ISP=ISPECIES(IAT1)
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      LNX1=LNX(ISP)
      ALLOCATE(NLPHIDOT(NR,LNX1))
      ALLOCATE(AEPHI(NR,LNX1))
      CALL SETUP$GETR8A('QPHIDOT',NR*LNX1,NLPHIDOT)
      CALL SETUP$GETR8A('AEPHI',NR*LNX1,AEPHI)
      LX=MAXVAL(LOX(:LNX(ISP),ISP))
      ALLOCATE(KPRIME(NR,LX+1))
      ALLOCATE(JBARPRIME(NR,LX+1))
      DO LN=1,LNX(ISP)
        IF(POTPAR(ISP)%LNSCATT(LN).NE.LN) CYCLE
        L=LOX(LN,ISP)
        RAD=POTPAR(ISP)%RAD
        QBAR=POTPAR(ISP)%QBAR(LN)
        CALL LMTO_KJBARTAILS(GID,NR,L,RAD,K2,QBAR,LAMBDA &
     &                                          ,KPRIME(:,L+1),JBARPRIME(:,L+1))
        DO IR=1,NR
          IF(R(IR).GT.RAD) EXIT
          JBARPRIME(IR,L+1)=NLPHIDOT(IR,LN)*POTPAR(ISP)%JBARTOPHIDOT(LN)
          CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,JVAL1,JDER1)
          KPRIME(IR,L+1)=(JVAL1-JBARPRIME(IR,L+1))/QBAR
          KPRIME(IR,L+1)=AEPHI(IR,LN)   *POTPAR(ISP)%KTOPHI(LN) &
     &                  +NLPHIDOT(IR,LN)*POTPAR(ISP)%KTOPHIDOT(LN)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == WORK OUT MAPPING                                                     ==
!     ==========================================================================
      RX=R(NR)
      LX=MAXVAL(LOX(:,ISP))
      LMX=(LX+1)**2
      DO L=0,LX
        I1=SBARLI1(L+1,ISP)
        IF(I1.LE.0) CYCLE
        LOFN(I1:I1+2*L)=L
        DO IM=1,2*L+1
          LMOFN(I1-1+IM)=L**2+IM
        ENDDO
      ENDDO
!
!     ==========================================================================
!     == WORK OUT MAPPING                                                     ==
!     ==========================================================================
      ORB(:,:)=0.D0
      ALLOCATE(KVAL(LX+1))
      ALLOCATE(JVAL(LX+1))
      ALLOCATE(YLM(LMX))
      DO IP=1,NP
!       == DR IS THE DISTANCE FROM THE SECOND ATOM ==========================
        DR(:)=P(:,IP)-R0(:,IAT1)
        DRLEN=SQRT(SUM(DR**2))
        IF(DRLEN.GT.RX) CYCLE
        CALL SPHERICAL$YLM(LMX,DR,YLM)
        DO L=0,LX           
          CALL RADIAL$VALUE(GID,NR,KPRIME(:,L+1),DRLEN,KVAL(L+1))
          CALL RADIAL$VALUE(GID,NR,JBARPRIME(:,L+1),DRLEN,JVAL(L+1))
        ENDDO
        DO I1=1,NORB
          ORB(IP,I1)=ORB(IP,I1)+KVAL(LOFN(I1)+1)*YLM(LMOFN(I1))
          DO I2=1,NORB
            ORB(IP,I1)=ORB(IP,I1)-JVAL(LOFN(I2)+1)*YLM(LMOFN(I2))*SBARMAT(I2,I1)
          ENDDO
        ENDDO
      ENDDO      
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
      WRITE(NFIL,FMT='(6(E12.6," "))')(((DATA(I,J,K),K=1,N3),J=1,N2),I=1,N1)
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
