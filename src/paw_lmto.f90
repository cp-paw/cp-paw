MODULE LMTO_MODULE
REAL(8)   ,PARAMETER  :: K2=-0.5D0      ! 0.5*K2 IS THE KINETIC ENERGY
REAL(8)   ,PARAMETER  :: RCSCALE=1.2D0  !2.D0  ! RADIUS SCALE FACTOR FOR NEIGHBORLIST
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
  REAL(8),POINTER    :: QBAR(:)
  REAL(8)   ,POINTER :: PHIDOTPROJ(:)   ! <P|PHIDOT>
  REAL(8)   ,POINTER :: KTOPHI(:)       ! K -> |PHI>KTOPHI+|PHIBARDOT> KTOPHIDOT 
  REAL(8)   ,POINTER :: KTOPHIDOT(:)    ! K -> |PHI>KTOPHI+|PHIBARDOT> KTOPHIDOT 
  REAL(8)   ,POINTER :: JBARTOPHIDOT(:) ! JBAR ->  |PHIBARDOT> JBARTOPHIDOT ====
  INTEGER(4),POINTER :: LNSCATT(:)      ! LN OF CORRESPONDING SCATTERING CHANNEL
  REAL(8)   ,POINTER :: DOVERLAPKK(:,:) !(LNX,LNX)
  REAL(8)   ,POINTER :: DOVERLAPKJ(:,:) !(LNX,LNX)
  REAL(8)   ,POINTER :: DOVERLAPJJ(:,:) !(LNX,LNX)
  TYPE(POTPARRED_TYPE)         :: SMALL
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSKPRIME
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSTAILEDK
  TYPE(ORBITALGAUSSCOEFF_TYPE) :: GAUSSTAILEDJBAR
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
INTEGER(4),ALLOCATABLE  :: LOX(:,:)            !(LXX,NSP)
INTEGER(4),ALLOCATABLE  :: ISPECIES(:)         !(NAT)
REAL(8)   ,ALLOCATABLE  :: ORBRAD(:,:) !(LXX+1,NAT) NODE-POSITION OF THE ORBITAL
TYPE(POTPAR_TYPE)     ,ALLOCATABLE :: POTPAR(:)!POTENTIAL PARAMETERS
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: SBAR(:)  !(NNS) SCREENED STRUCTURE CONST.
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: OVERLAP(:) !(NNS) OVERLAP MATRIX ONLY MAIN CHANNEL
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
TYPE(UTENSOR_TYPE)    ,ALLOCATABLE :: UTENSOR(:)
END MODULE LMTO_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_INITIALIZE()
!     **************************************************************************
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
      CALL LMTO_SBARINDICES2()   ! SBARLI1
!
!     ==========================================================================
!     == DETERMINE POTENTIAL PARAMETERS                                       ==
!     ==========================================================================
      CALL LMTO_MAKEPOTPAR()
!
!     ==========================================================================
!     == DETERMINE GAUSS EXPANDION OF UNSCREENED HANKEL FUNCTIONS KPRIME      ==
!     ==========================================================================
      CALL LMTO_GAUSSFITKPRIME()
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
!     **                                                                      **
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
          CALL REPORT$R8VAL(NFIL,'QBAR',POTPAR(ISP)%QBAR(LN),'')
          CALL REPORT$R8VAL(NFIL,'KTOPHI',POTPAR(ISP)%KTOPHI(LN),'')
          CALL REPORT$R8VAL(NFIL,'KTOPHIDOT',POTPAR(ISP)%KTOPHIDOT(LN),'')
          CALL REPORT$R8VAL(NFIL,'JBARTOPHIDOT',POTPAR(ISP)%JBARTOPHIDOT(LN),'')
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
      REAL(8)   ,ALLOCATABLE :: EOFLN(:)
      REAL(8)   ,ALLOCATABLE :: ESCATT(:)
      REAL(8)   ,ALLOCATABLE :: NLPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: PSPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: PRO(:,:)
      REAL(8)   ,ALLOCATABLE :: KPRIME(:)
      REAL(8)   ,ALLOCATABLE :: W(:)
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
      REAL(8)                :: LAMBDA
      INTEGER(4)             :: ISP,LN,L,LN1,LN2,LX
REAL(8) ::Y(20)
INTEGER(4) :: I,J,IR,L0
      CHARACTER(128)         :: STRING
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO_MAKEPOTPAR')
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ALLOCATE(POTPAR(NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        LNX1=LNX(ISP)
        LX=MAXVAL(LOX(:LNX1,ISP))
        ALLOCATE(ISCATT(LNX1))
        ALLOCATE(EOFLN(LNX1))
        ALLOCATE(ESCATT(LNX1))
        CALL SETUP$GETI4A('ISCATT',LNX1,ISCATT)
        CALL SETUP$GETR8A('EOFLN',LNX1,EOFLN)
        CALL SETUP$GETR8A('ESCATT',LNX1,ESCATT)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        ALLOCATE(NLPHI(NR,LNX1))
        ALLOCATE(AEPHI(NR,LNX1))
        ALLOCATE(NLPHIDOT(NR,LNX1))
        ALLOCATE(AEPHIDOT(NR,LNX1))
        ALLOCATE(PSPHI(NR,LNX1))
        ALLOCATE(PSPHIDOT(NR,LNX1))
        ALLOCATE(PRO(NR,LNX1))
        CALL SETUP$GETR8A('QPHI',NR*LNX1,NLPHI)
        CALL SETUP$GETR8A('AEPHI',NR*LNX1,AEPHI)
        CALL SETUP$GETR8A('QPHIDOT',NR*LNX1,NLPHIDOT)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX1,AEPHIDOT)
        CALL SETUP$GETR8A('PSPHI',NR*LNX1,PSPHI)
        CALL SETUP$GETR8A('PSPHIDOT',NR*LNX1,PSPHIDOT)
        CALL SETUP$GETR8A('PRO',NR*LNX1,PRO)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(NINT(AEZ),'R(ASA)',RAD) 
        POTPAR(ISP)%RAD=RAD
!
        ALLOCATE(POTPAR(ISP)%QBAR(LNX1))
        ALLOCATE(POTPAR(ISP)%PHIDOTPROJ(LNX1))
        ALLOCATE(POTPAR(ISP)%LNSCATT(LNX1))
        ALLOCATE(POTPAR(ISP)%KTOPHI(LNX1))
        ALLOCATE(POTPAR(ISP)%KTOPHIDOT(LNX1))
        ALLOCATE(POTPAR(ISP)%JBARTOPHIDOT(LNX1))
        ALLOCATE(KPRIME(NR))
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        ALLOCATE(NLPHI(NR,LNX1))
        ALLOCATE(AEPHI(NR,LNX1))
        ALLOCATE(NLPHIDOT(NR,LNX1))
        ALLOCATE(AEPHIDOT(NR,LNX1))
        ALLOCATE(PSPHI(NR,LNX1))
        ALLOCATE(PSPHIDOT(NR,LNX1))
        ALLOCATE(PRO(NR,LNX1))
        CALL SETUP$GETR8A('QPHI',NR*LNX1,NLPHI)
        CALL SETUP$GETR8A('AEPHI',NR*LNX1,AEPHI)
        CALL SETUP$GETR8A('QPHIDOT',NR*LNX1,NLPHIDOT)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX1,AEPHIDOT)
        CALL SETUP$GETR8A('PSPHI',NR*LNX1,PSPHI)
        CALL SETUP$GETR8A('PSPHIDOT',NR*LNX1,PSPHIDOT)
        CALL SETUP$GETR8A('PRO',NR*LNX1,PRO)

        ALLOCATE(W(NR))
!
        DO L=0,MAXVAL(LOX(:,ISP))
!         == SELECT PHIBARDOT FROM VALENCE CHANNEL =============================
          LN1=0
          DO LN=1,LNX(ISP)
            IF(LOX(LN,ISP).NE.L) CYCLE
            IF(ISCATT(LN).GT.0) CYCLE
            LN1=LN
          ENDDO
          IF(LN1.EQ.0) THEN
            CALL ERROR$MSG('SELECTION OF ENERGY DERIVATIVE PARTIAL WAVE FAILED')
            CALL ERROR$STOP('LMTO_MAKEPOTPAR')
          END IF
          DO LN=1,LNX(ISP)
            IF(LOX(LN,ISP).NE.L) CYCLE
!LN1=LN ! OLD VERSION RE-ESTABLISHED
!           ====================================================================
!           == VALUE AND DERIVATIVE OF PARTIAL WAVES AND ENVELOPE FUNCTIONS   ==
!           == PHIDOT IS PHIBARDOT                                            ==
!           == THERE IS ONLY A SINGLE PHIDOT PER ANGULAR MOMENTUM             ==
!           ====================================================================
            CALL RADIAL$VALUE(GID,NR,NLPHIDOT(:,LN1),RAD,PHIDOTVAL)
            CALL RADIAL$DERIVATIVE(GID,NR,NLPHIDOT(:,LN1),RAD,PHIDOTDER)
            CALL RADIAL$VALUE(GID,NR,NLPHI(:,LN),RAD,PHIVAL)
            CALL RADIAL$DERIVATIVE(GID,NR,NLPHI(:,LN),RAD,PHIDER)
            CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
            CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
!
!           ====================================================================
!           == CALCULATE POTENTIAL PARAMETERS                                 ==
!           ====================================================================
            WJPHI=JVAL*PHIDER-JDER*PHIVAL
            WJPHIDOT=JVAL*PHIDOTDER-JDER*PHIDOTVAL
            WKPHI=KVAL*PHIDER-KDER*PHIVAL
            WKPHIDOT=KVAL*PHIDOTDER-KDER*PHIDOTVAL
            WPHIPHIDOT=PHIVAL*PHIDOTDER-PHIDER*PHIDOTVAL
            QBAR=WJPHIDOT/WKPHIDOT
            WJBARPHI=WJPHI-WKPHI*QBAR
!
!           ====================================================================
!           == K    -> |PHI>KTOPHI+|PHIBARDOT> KTOPHIDOT =======================
!           == JBAR ->             |PHIBARDOT> JBARTOPHIDOT ====================
!           ====================================================================
            POTPAR(ISP)%LNSCATT(LN)=LN1
            POTPAR(ISP)%QBAR(LN)=QBAR
            POTPAR(ISP)%KTOPHI(LN)=WKPHIDOT/WPHIPHIDOT
            POTPAR(ISP)%KTOPHIDOT(LN)=-WKPHI/WPHIPHIDOT
            POTPAR(ISP)%JBARTOPHIDOT(LN)=-WJBARPHI/WPHIPHIDOT
!
!           ==  <PRO|PSPHIDOT> =================================================
            CALL RADIAL$INTEGRAL(GID,NR,R**2*PRO(:,LN)*PSPHIDOT(:,LN1),SVAR)
            POTPAR(ISP)%PHIDOTPROJ(LN)=SVAR
!           == CROSSCHECK BIOTHOGONALITY =======================================
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
        ENDDO        
!
        DEALLOCATE(EOFLN)
        DEALLOCATE(ESCATT)
        DEALLOCATE(NLPHI)
        DEALLOCATE(AEPHI)
        DEALLOCATE(NLPHIDOT)
        DEALLOCATE(AEPHIDOT)
        DEALLOCATE(PSPHIDOT)
        DEALLOCATE(PSPHI)
        DEALLOCATE(PRO)
        DEALLOCATE(ISCATT)
        DEALLOCATE(R)
        DEALLOCATE(W)
        DEALLOCATE(KPRIME)
      ENDDO
                             CALL TRACE$POP()
      RETURN
      END
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
      LOGICAL(4)             :: TTEST=.TRUE.
      INTEGER(4),PARAMETER   :: NEX=30
      REAL(8)   ,PARAMETER   :: LAMBDA=1.D0
      REAL(8)                :: E(10)
      INTEGER(4)             :: NPOW
      INTEGER(4)             :: NPOW2
      INTEGER(4)             :: NE
      INTEGER(4)             :: LMX
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
      REAL(8)   ,ALLOCATABLE :: KPRIMEL(:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: CPOWK(:,:)    !(NPOW2,NE)
      REAL(8)   ,ALLOCATABLE :: CPOWJBAR(:,:) !(NPOW2,NE)
      REAL(8)   ,ALLOCATABLE :: CPOWKL(:,:)   !(NPOW2,NE)
      REAL(8)   ,ALLOCATABLE :: CK(:,:)       !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: CJBAR(:,:)    !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: CKL(:,:)      !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: WORKK(:,:)    !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: WORKJ(:,:)    !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: WORKKL(:,:)   !(NIJK,NE)
      REAL(8)   ,ALLOCATABLE :: YLMPOL(:,:)
      REAL(8)                :: RAD
      REAL(8)                :: QBAR
      REAL(8)                :: SVAR,SVAR1,SVAR2,JVAL,KVAL,JDER,KDER
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
!     ** CONSTRUCTS THE UNSCREENED HANKEL FUNCTIONS, WITH THE SIMGULARITY     **
!     ** REPLACED BY A NODELESS SCATTERING FUNCTION AND DETERMINES ITS        **
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
      LOGICAL(4),PARAMETER   :: TTEST=.TRUE.
      REAL(8)                :: E(NEX)
      INTEGER(4)             :: NPOW
      INTEGER(4)             :: NPOW2
      INTEGER(4)             :: NE
      INTEGER(4)             :: LMX
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
      REAL(8)                :: SVAR,SVAR1,SVAR2,JVAL,KVAL,JDER,KDER
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
      REAL(8)   ,ALLOCATABLE :: R0(:,:)
      REAL(8)                :: RBAS(3,3)
      REAL(8)                :: RA(3),RB(3)
      REAL(8)   ,ALLOCATABLE :: SCALE1(:)
INTEGER(4) :: IE,IND,K,IX
REAL(8)    :: SVAR,ARR(20),X
REAL(8)    :: EULER
REAL(8),ALLOCATABLE ::COO(:,:,:)
INTEGER(4),PARAMETER  :: NX=100
REAL(8)   ,PARAMETER  :: XMAX=10.D0
REAL(8)               :: RGRID(3,NX)
REAL(8)               :: FGRID(NX,20,2)
!     **************************************************************************
DO I=1,NX
  RGRID(:,I)=(/1.D0,1D0,1.D0/)/SQRT(3.D0)*XMAX*2.D0*(REAL(I-1)/REAL(NX-1)-0.5D0)
ENDDO
      EULER=EXP(1.D0)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      CALL CELL$GETR8A('T0',9,RBAS)
!
!     ==========================================================================
!     == ALLOCATE GAUSSORB, WHICH HOLD THE LOCAL ORBITALS IN GAUSS REPRESENTAT.
!     ==========================================================================
      TNEW=.NOT.ALLOCATED(GAUSSORB)
      IF(TNEW) THEN
        ALLOCATE(GAUSSORB(NAT))
        NAT=SIZE(ISPECIES)
        DO IAT=1,NAT
          ISP=ISPECIES(IAT)
          NE=4
          NPOW=4 !POLYNOMIAL HAS X^NPOW AS HIGHEST POWER
!!$NE=1
!!$NPOW=2
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
        ALLOCATE(GCB(NIJKA,NEA,NORBX))
        ALLOCATE(MAT(NIJKA*NEA,NIJKA*NEA))
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
        CALL LMTO_GAUSSPROJECT(NIJKA,NEA,EA,MAT)
        GCB(:,:,:NORB1)=GAUSSORB(IAT)%C(:,:,:)
        CALL LIB$MATRIXSOLVER8(NIJKA*NEA,NIJKA*NEA,NORB1 &
    &                                      ,MAT,GAUSSORB(IAT)%C,GCB(:,:,:NORB1))
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
!!$!STOP 'XC'
        DEALLOCATE(EA)
        DEALLOCATE(GCB)
        DEALLOCATE(MAT)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_GAUSSPROJECT(NIJK,NE,E,U)
!     **************************************************************************
!     **                                                                      **
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
      FACTOR(:)=0.D0
      SVAR=SQRT(2.D0*PI)
      FACTOR(0)=SVAR
      DO I=1,IMAX
        SVAR=SVAR*REAL(2*I-1,KIND=8)
        FACTOR(2*I)=SVAR
      ENDDO
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
!     **                                                                      **
!     **  PRODUCES A NEIGHBORLIST AND THE SCREENED STRUCTURE CONSTANTS        **
!     **     SBAR                                                             **
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
     &                       ,GAUSSORB
      IMPLICIT NONE
      TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: MAP(:) !(NNS)
      INTEGER(4)          :: NNS
      INTEGER(4)          :: NN,NN1,NN2
      INTEGER(4)          :: N1,N2,N2A,N2B
      INTEGER(4)          :: IAT1,IAT2,IAT2A,IAT2B,NEA
      INTEGER(4)          :: ISP1,ISP2,ISP2A,ISP2B,NEB
      INTEGER(4)          :: NIJKA,NIJKB
      REAL(8)             :: RA(3),RB(3)
      REAL(8),ALLOCATABLE :: EA(:),EB(:)
      REAL(8)             :: RBAS(3,3)
      REAL(8),ALLOCATABLE :: R0(:,:)
      REAL(8)             :: DR(3)
      REAL(8)             :: IT(3)
      INTEGER(4)          :: NAT
      INTEGER(4)          :: N1X
      REAL(8)   ,ALLOCATABLE :: BAREOV(:,:)
!     **************************************************************************
                          CALL TRACE$PUSH('LMTO$OVERLAP_FULL')
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      CALL CELL$GETR8A('T0',9,RBAS)

      NEA=MAXVAL(POTPAR(:)%GAUSSKPRIME%NE)
      ALLOCATE(EA(NEA))
      ALLOCATE(EB(NEA))
!
!     ==========================================================================
!     == KBAR=K*MAP;     MAP=(1+QBAR*SBAR)                                    ==
!     ==========================================================================
      NNS=SIZE(SBAR)
      ALLOCATE(MAP(NNS))
      CALL LMTO_KTOKBAR(NNS,SBAR,MAP)
!
!     ==========================================================================
!     == LIMIT SPHERE OVERLAP TERMS TO ONE CHANNEL PER ANGULAR MOMENTUM       ==
!     == RESULTING OVERLAP MATRIX IS SMALLER AND THE FULL MATRIX IS READILY   ==
!     == OBTAINED BY ADDING ONSITE TERMS ONLY.                                ==
!     == RESULT LIES IN POTPAR(ISP)%SMALL%DOVERLAPKK, ETC.                    ==
!     ==========================================================================
      CALL LMTO_EXPANDSPHEREDO('LM')
!
!     ==========================================================================
!     ==  CREATE EMPTY ARRAY OVERLAP (INDICES FROM SBAR, MAT IS ZEROED)       ==
!     ==  OF TWO VECTORS
!     ==========================================================================
      NNS=SIZE(SBAR)
      IF(.NOT.ALLOCATED(OVERLAP))ALLOCATE(OVERLAP(NNS))
      CALL LMTO_PERIODICMAT_CP(NNS,SBAR,OVERLAP)  !COPY FRAME
      N1X=0
      DO NN=1,NNS
        OVERLAP(NN)%MAT(:,:)=0.D0
        N1X=MAX(N1X,OVERLAP(NN)%N1,OVERLAP(NN)%N2)
      ENDDO
!
!     ==========================================================================
!     ==  ACCUMULATE NTBO OVERLAP MATRIX:GAUSSIAN CONTRIBUTION                ==
!     ==========================================================================
      ALLOCATE(BAREOV(N1X,N1X))
      DO NN=1,NNS
        N1=OVERLAP(NN)%N1
        N2=OVERLAP(NN)%N2
        IAT1=OVERLAP(NN)%IAT1
        IAT2=OVERLAP(NN)%IAT2
        IT=OVERLAP(NN)%IT
        RA(:)=R0(:,IAT1)
        RB(:)=R0(:,IAT2)+MATMUL(RBAS,REAL(IT))
        NIJKA=GAUSSORB(IAT1)%NIJK
        NIJKB=GAUSSORB(IAT2)%NIJK
        NEA=GAUSSORB(IAT1)%NE
        NEB=GAUSSORB(IAT2)%NE
        CALL GAUSSIAN$NEWOVERLAP( &
     &                      NIJKA,N1,NEA,GAUSSORB(IAT1)%E,RA,GAUSSORB(IAT1)%C &
     &                     ,NIJKB,N2,NEB,GAUSSORB(IAT2)%E,RB,GAUSSORB(IAT2)%C &
     &                     ,OVERLAP(NN)%MAT)
      ENDDO
!
!     ==========================================================================
!     ==  ACCUMULATE OVERLAP MATRIX: SPHERE CORRECTION                        ==
!     ==========================================================================
      DO NN=1,NNS
!       ==  <KBAR(IAT1)|KBAR(IAT2)> ============================================
        N1=OVERLAP(NN)%N1
        N2=OVERLAP(NN)%N2
        IAT1=OVERLAP(NN)%IAT1
        IAT2=OVERLAP(NN)%IAT2
        ISP1=ISPECIES(IAT1)
        ISP2=ISPECIES(IAT2)
!
!       == ONSITE TERM OF HEADS: <KOMEGA(AT1)|KOMEGA(AT2)> =====================
        IF(IAT1.EQ.IAT2) THEN
          IT=OVERLAP(NN)%IT
          IF(IT(1).EQ.0.AND.IT(2).EQ.0.AND.IT(3).EQ.0) THEN
            OVERLAP(NN)%MAT=OVERLAP(NN)%MAT+POTPAR(ISP1)%SMALL%DOVERLAPKK
          END IF
        END IF
!
        DO NN1=1,NNS
          IF(SBAR(NN1)%IAT1.NE.IAT1) CYCLE
!         ==  <KBAR(IAT1)|=SUM(AT2A):MAP(AT1,AT2A)<KBARPRIME(AT2A)| ============
!         ==  <KBAR(IAT1)|=<KBAROMEGA(AT1)|                          ===========
!         ==              +SUM(AT2A):SBAR(AT1,AT2A)<JBAROMEGA(AT2A)| ===========
          IAT2A=SBAR(NN1)%IAT2
          ISP2A=ISPECIES(IAT2A)
          DO NN2=1,NNS
            IF(SBAR(NN2)%IAT1.NE.IAT2) CYCLE
!           ==  <KBAR(IAT2)|=SUM(AT2B):MAP(AT2,AT2B)<KBARPRIME(AT2B)| ==========
!           ==  <KBAR(IAT2)|=<KBAROMEGA(AT2)|                          =========
!           ==              +SUM(AT2B):SBAR(AT2,AT2B)<JBAROMEGA(AT2B)| =========
            IAT2B=SBAR(NN2)%IAT2
            ISP2B=ISPECIES(IAT2B)
!
!           == TAIL-TAIL OVERLAP: SBAR<JBAR|JBAR>SBAR ==========================
            IF(IAT2A.EQ.IAT2B) THEN
              IT=-SBAR(NN1)%IT+OVERLAP(NN)%IT+SBAR(NN2)%IT
              IF(IT(1).EQ.0.AND.IT(2).EQ.0.AND.IT(3).EQ.0) THEN
                OVERLAP(NN)%MAT=OVERLAP(NN)%MAT &
    &                   +MATMUL(TRANSPOSE(SBAR(NN1)%MAT) &
                          ,MATMUL(POTPAR(ISP2A)%SMALL%DOVERLAPJJ,SBAR(NN2)%MAT))
              END IF
            END IF
          ENDDO ! END OF LOOP OVER NN2
!
!         == TAIL-HEAD OVERLAP:  SBAR<JBAR|K> ==================================
          IF(IAT2A.EQ.IAT2) THEN
            IT=-SBAR(NN1)%IT+OVERLAP(NN)%IT
            IF(IT(1).EQ.0.AND.IT(2).EQ.0.AND.IT(3).EQ.0) THEN
              OVERLAP(NN)%MAT=OVERLAP(NN)%MAT &
     &           +TRANSPOSE(MATMUL(POTPAR(ISP2)%SMALL%DOVERLAPKJ,SBAR(NN1)%MAT))
            END IF
          END IF
        ENDDO  ! END OF LOOP OVER NN1
        DO NN2=1,NNS
          IF(SBAR(NN2)%IAT1.NE.IAT2) CYCLE
          IAT2B=SBAR(NN2)%IAT2
!
!         == HEAD-TAIL, OVERLAP:   <K|JBAR>SBAR ================================
          IF(IAT1.EQ.IAT2B) THEN
            IT=OVERLAP(NN)%IT+SBAR(NN2)%IT
            IF(IT(1).EQ.0.AND.IT(2).EQ.0.AND.IT(3).EQ.0) THEN
              OVERLAP(NN)%MAT=OVERLAP(NN)%MAT &
    &                 +MATMUL(POTPAR(ISP1)%SMALL%DOVERLAPKJ,SBAR(NN2)%MAT)
            END IF
          END IF
        ENDDO  ! END OF LOOP OVER NN2
      ENDDO  ! END OF LOOP OVER NN
      DEALLOCATE(BAREOV)
     
                          CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$TESTENERGY()
      USE LMTO_MODULE, ONLY : GAUSSORB,GAUSSORB_T
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) :: NAT
      INTEGER(4) :: IAT
!     **************************************************************************
!RETURN
      WRITE(*,FMT='(82("="),T30," TESTENERGY START ")')
      CALL LMTO$REPORTPOTBAR(6)
      CALL LMTO$REPORTSBAR(6)
!
!     ==========================================================================
!     == CONSTRUCT TIGHT-BINDING ORBITALS IN GAUSSIAN REPRESENTATION          ==
!     ==========================================================================
PRINT*,'DOING LMTO_NTBFROMTAILEDKJ.....'
      CALL LMTO_NTBOFROMTAILEDKJ()
PRINT*,'.....LMTO_NTBFROMTAILEDKJ DONE'
!
PRINT*,'DOING LMTO_NTBFROMKPRIME.....'
      CALL LMTO_NTBOFROMKPRIME() 
PRINT*,'..... LMTO_NTBFROMKPRIME DONE'
!GAUSSORB=GAUSSORB_T
!
!     ==========================================================================
!     ==  TEST COEFFICIENTS OF TIGHT-BINDING ORBITALS                         ==
!     ==========================================================================
PRINT*,' BEFORE TESTNTBO.....'
      CALL LMTO_TESTNTBO()
PRINT*,'....... TESTNTBO DONE'
!
!     ==========================================================================
!     ==  CALCULATE OVERLAP MATRIX                                            ==
!     ==========================================================================
PRINT*,'DOING LMTO$OVERLAPFULL.....'
      CALL LMTO$OVERLAPFULL()
PRINT*,'..... LMTO$OVERLAPFULL DONE'
      CALL LMTO$REPORTOVERLAP(6)
!
!     ==========================================================================
!     == TESTS THE OVERLAP MATRIX OF NATURAL TIGHT-BINDING ORBITALS           ==
!     == BY ESTIMATING THE OVERLAP BETWEEN KOHN SHAM WAVE FUNCTIONS USING THE ==
!     == COEFFICIENTS IN NTBS AND THEIR OVERLAP MATRIX                        ==
!     ==========================================================================
PRINT*,' BEFORE TESTOVERLAP.....'
      CALL LMTO_TESTOVERLAP()
PRINT*,'....... TESTOVERLAP DONE'
!
!     ==========================================================================
!     ==  
!     ==========================================================================
PRINT*,' BEFORE LMTO_FOURCENTER.....'
      call lmto_fourcenter()
PRINT*,'....... LMTO_FOURCENTER DONE'
!
!     ==========================================================================
!     ==  
!     ==========================================================================
PRINT*,'DOING LMTO_PLOTLOCORB.....'
      CALL ATOMLIST$NATOM(NAT)
      DO IAT=1,NAT
!        CALL LMTO_PLOTLOCORB(IAT)
      ENDDO
PRINT*,'..... LMTO_PLOTLOCORB DONE'
CALL ERROR$MSG('FORCED STOP')
CALL ERROR$STOP('LMTO$TESTENERGY')
      RETURN
      END
!!$!
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
      USE WAVES_MODULE, ONLY: NKPTL,NSPIN,NDIM,THIS,MAP,WAVES_SELECTWV
      USE LMTO_MODULE, ONLY : PERIODICMAT_TYPE,SBAR,ISPECIES,POTPAR &
     &                       ,SBARATOMI1,SBARATOMI2,SBARLI1,LNX,LOX 
      IMPLICIT NONE
      REAL(8)       :: XK(3,NKPTL)
      INTEGER(4)    :: NPRO
      INTEGER(4)    :: NBH,NB
      INTEGER(4)    :: IKPT,ISPIN,IBH,JBH,IPRO,JPRO,IDIM,IB
      COMPLEX(8)    :: CSVAR1,CSVAR2
      COMPLEX(8),ALLOCATABLE :: OVER(:,:)
      REAL(8)   ,ALLOCATABLE :: OVERMAT(:,:)
      LOGICAL(4)             :: TINV
!     **************************************************************************
!
!     ==========================================================================
!     ==  GET K-POINTS IN RELATIVE COORDINATES                                ==
!     ==========================================================================
      CALL WAVES_DYNOCCGETR8A('XK',3*NKPTL,XK)
      NPRO=MAP%NPRO
!
!     ==========================================================================
!     ==  DETERMINE OVERLAP MATRIX  IN K-SPACE                                ==
!     ==========================================================================
!!!! PARALLELIZE LOOP WITH RESPECT TO STATES.
      ALLOCATE(OVER(NPRO,NPRO))
      DO IKPT=1,NKPTL
        CALL WAVES_SELECTWV(IKPT,1)
        NBH=THIS%NBH
        NB=THIS%NB
        CALL PLANEWAVE$GETL4('TINV',TINV)
PRINT*,'TINV ',TINV
PRINT*,'MARKE 3A',IKPT,NBH,XK(:,IKPT),NPRO
        CALL LMTO_OVERLAPOFK(XK(:,IKPT),NPRO,OVER)
PRINT*,'===================== OVERLAP (REAL PART) ===================='
DO JPRO=1,NPRO
  WRITE(*,FMT='("R(O)",80F10.5)')REAL(OVER(:,JPRO))
ENDDO
PRINT*,'===================== OVERLAP (IMAGINARY PART) ==============='
DO JPRO=1,NPRO
  WRITE(*,FMT='("I(O)",80F10.5)')AIMAG(OVER(:,JPRO))
ENDDO
        DO ISPIN=1,NSPIN
PRINT*,'===================== ISPIN ',ISPIN,'========================='
          CALL WAVES_SELECTWV(IKPT,ISPIN)

PRINT*,'=============WAVE FUNCTION COEFFICIENTS FOR NTBOS  ================='
DO IBH=1,NBH
  WRITE(*,FMT='("C",I5,80F10.5)')2*IBH-1,REAL(THIS%TBC(1,IBH,:))
  WRITE(*,FMT='("C",I5,80F10.5)')2*IBH,AIMAG(THIS%TBC(1,IBH,:))
ENDDO
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
      &                                  *OVER(IPRO,JPRO)*THIS%TBC(IDIM,JBH,JPRO)
                      CSVAR2=CSVAR2+THIS%TBC(IDIM,IBH,IPRO) &
      &                                  *OVER(IPRO,JPRO)*THIS%TBC(IDIM,JBH,JPRO)
                    ENDDO
                  ENDDO
                ENDDO
                OVERMAT(2*IBH-1,2*JBH-1)=0.5D0*REAL(CSVAR1+CSVAR2)
                OVERMAT(2*IBH-1,2*JBH)  =0.5D0*AIMAG(CSVAR1+CSVAR2)
                OVERMAT(2*IBH,2*JBH-1) =-0.5D0*AIMAG(CSVAR1-CSVAR2)
                OVERMAT(2*IBH,2*JBH)    =0.5D0*REAL(CSVAR1-CSVAR2)
              ELSE
                CSVAR1=(0.D0,0.D0)
                DO IPRO=1,NPRO
                  DO JPRO=1,NPRO
                    DO IDIM=1,NDIM
                      CSVAR1=CSVAR1+CONJG(THIS%TBC(IDIM,IBH,IPRO)) &
      &                                  *OVER(IPRO,JPRO)*THIS%TBC(IDIM,JBH,JPRO)
                    ENDDO
                  ENDDO
                ENDDO
                OVERMAT(IBH,JBH)=REAL(CSVAR1)
              END IF
            ENDDO
          ENDDO
DO IB=1,NB
  WRITE(*,FMT='("UNITY ",I5,20F10.5)')IB,OVERMAT(:,IB)
ENDDO
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
!     **  TIGHT-BINDING ORBITALS. IFG THEY ARE IDENTICAL TO THE ORIGINAL      **
!     ** PAW PROJECTIONS <PTILDE|PSITILDE>, THE NTB COEFFICIENTS ARE CORRECT  **
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
WRITE(*,FMT='("PI",I5,80F10.5)')2*IBH-1,REAL(THIS%PROJ(1,IBH,:))
WRITE(*,FMT='("PF",I5,80F10.5)')2*IBH-1,REAL(TBC1(1,:))
WRITE(*,*)
WRITE(*,FMT='("PI",I5,80F10.5)')2*IBH,AIMAG(THIS%PROJ(1,IBH,:))
WRITE(*,FMT='("PF",I5,80F10.5)')2*IBH,AIMAG(TBC1(1,:))
WRITE(*,*)
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
      subroutine lmto_fourcenter()
!     **************************************************************************
!     **************************************************************************
      USE LMTO_MODULE, ONLY : UTENSOR,nsp,ISPECIES,LNX,LOX,GAUSSORB,SBARLI1
      IMPLICIT NONE
      integer(4)             :: gid
      integer(4)             :: nr
      integer(4)             :: lmrx
      real(8)   ,allocatable :: r(:)
      real(8)   ,allocatable :: nlphi(:,:)
      real(8)   ,allocatable :: nlphidot(:,:)
      real(8)   ,allocatable :: aephi(:,:)
      real(8)   ,allocatable :: aephidot(:,:)
      real(8)   ,allocatable :: ulittle(:,:,:,:,:)
      integer(4)             :: isp,lrx,lnx1
!     **************************************************************************
      CALL LMTO_FOURCENTERGAUSS()
!
      do isp=1,nsp
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
!
        CALL SETUP$GETI4('LMRX',LMRX)
        lrx=INT(SQRT(REAL(LMRX)+1.D-8))-1
!
        lnx1=lnx(isp)
        ALLOCATE(NLPHI(NR,LNX1))
        ALLOCATE(AEPHI(NR,LNX1))
        ALLOCATE(NLPHIDOT(NR,LNX1))
        ALLOCATE(AEPHIDOT(NR,LNX1))
        CALL SETUP$GETR8A('QPHI',NR*LNX1,NLPHI)
        CALL SETUP$GETR8A('AEPHI',NR*LNX1,AEPHI)
        CALL SETUP$GETR8A('QPHIDOT',NR*LNX1,NLPHIDOT)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX1,AEPHIDOT)
        ALLOCATE(ULITTLE(LRX+1,LNX1,LNX1,LNX1,LNX1))
!        CALL LDAPLUSU_ULITTLE(GID,NR,LRX,LNX,LOX(:lnx(isp),isp),CHI,ULITTLE)

!        CALL LDAPLUSU_UTENSOR(LRX,NCHI,LNX,LOX,ULITTLE,U)
        deALLOCATE(ULITTLE)
        deALLOCATE(nlphi)
        deALLOCATE(nlphidot)
        deALLOCATE(aephi)
        deALLOCATE(aephidot)
        deALLOCATE(r)
      enddo
      return
      end
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE Lmto_ULITTLE(GID,NR,LRX,LNX,LOX,CHI,ULITTLE)
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
                               CALL TRACE$PUSH('LMTO$PROJTONTBO')
      CALL LMTO_PREPARE1(NPRO,NRL,A,B,C,D,E,F)
WRITE(*,FMT='("A=",20F10.5)')A
WRITE(*,FMT='("B=",20F10.5)')B
WRITE(*,FMT='("C=",20F10.5)')C
WRITE(*,FMT='("D=",20F10.5)')D
WRITE(*,FMT='("E=",20F10.5)')E
WRITE(*,FMT='("F=",20F10.5)')F
      CALL LMTO_PREPARE2(XK,NRL,C,E,F,G,H)
DO I=1,NRL
  WRITE(*,FMT='("G=",20F10.5)')G(I,:)
ENDDO
DO I=1,NRL
  WRITE(*,FMT='("H=",20F10.5)')H(I,:)
ENDDO
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
                               CALL TRACE$PUSH('LMTO$PROJTONTBO')
      CALL LMTO_PREPARE1(NPRO,NRL,A,B,C,D,E,F)
      CALL LMTO_PREPARE2(XK,NRL,C,E,F,G,H)
!
!     ==========================================================================
!     == D * DE/DQ
!     ==========================================================================
      ALLOCATE(PROJCONTR1(NDIM,NBH,NRL))      
      PROJCONTR1(:,:,:)=(0.D0,0.D0)
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
!     **************************************************************************
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
        KR=SUM(REAL(SBAR(NN)%IT(:))*XK(:))
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
!     **************************************************************************
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
        KR=SUM(REAL(ALIST(NN)%IT(:))*XK(:))
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
      SUBROUTINE LMTO$DOLOCORB(IAT,ISP,GID,NR,LNXCHI,LNXPHI,TORB,CHIPHI,CHI)
!     **************************************************************************
!     **  CONSTRUCTS ONSITE-MAPPING FROM PARTIAL WAVES TO LOCAL ORBITALS      **
!     **                                                                      **
!     **  1) SEMI-CORE-LIKE PARTIAL WAVES ARE TREATED AS LOCAL ORBITALS       **
!     **     EVEN THOUGH THEY MAY NOT BE INCLUDED IN THE SET                  **
!     **  2) VALENCE-LIKE LOCAL ORBITALS ARE CONSTRUCTED FROM THE VALENCE-LIKE**
!     **     PARTIAL WAVE AND ITS |Q_{N+1}> PARTNER. THAT ARE SUPER-IMPOSED   **
!     **     AS IN A SCREENED LMTO.                                           **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : SBAR,POTPAR
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
      INTEGER(4),ALLOCATABLE :: LOXCHI(:)    ! ANGULAR MOMENTA
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
      REAL(8)   ,ALLOCATABLE :: SBARONSITE(:,:)   !(LMX,LMX) 
      REAL(8)   ,ALLOCATABLE :: AMAT(:,:)   !(LNXCHI1,LNXPHI) 
      REAL(8)   ,ALLOCATABLE :: BMAT(:,:)   !(LNXCHI1,LNXCHI1) 
      REAL(8)   ,ALLOCATABLE :: XMAT(:,:)   !(LNXPHI,LNXCHI) 
      REAL(8)                :: RCOV              ! COVALENT RADIUS
      REAL(8)                :: AEZ               ! ATOMIC NUMBER
      REAL(8)                :: AUX(NR)
      REAL(8)                :: R(NR)
      REAL(8)                :: SVAR,SVAR1,VAL,DER,VALD,DERD
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: LMX,LNXCHI1
      INTEGER(4)             :: LN,L,I,J,IIB,N1,LM,LNCHI1,LNCHI,IR
      INTEGER(4)             :: LNPHI
      INTEGER(4)             :: IRCOV  ! GRID INDEX JUST BEYOND RCOV
      CHARACTER(64)          :: STRING
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO$DOLOCORB')
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == COLLECT DATA                                                         ==
!     ==========================================================================
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETR8('AEZ',AEZ)
      CALL PERIODICTABLE$GET(NINT(AEZ),'R(COV)',RCOV)
      DO IR=1,NR
        IRCOV=IR
        IF(R(IR).GT.RCOV) EXIT
      ENDDO
      CALL SETUP$GETI4('LNX',LNX)
      IF(LNXPHI.NE.LNX) THEN
        CALL ERROR$STOP('INCONSISTENT #(PARTIAL WAVES)')
        CALL ERROR$STOP('LMTO$DOLOCORB')
      END IF
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUP$GETR8A('PSPHI',NR*LNX,PSPHI)
      CALL SETUP$GETR8A('NLPHI',NR*LNX,NLPHI)
      CALL SETUP$GETR8A('AEPHIDOT',NR*LNX,AEPHIDOT)
      CALL SETUP$GETR8A('PSPHIDOT',NR*LNX,PSPHIDOT)
      CALL SETUP$GETR8A('NLPHIDOT',NR*LNX,NLPHIDOT)
      CALL SETUP$GETR8A('PRO',NR*LNX,PRO)
!
!     ==========================================================================
!     == FIND ONSITE STRUCTURE CONSTANTS                                      ==
!     ==========================================================================
      NNS=SIZE(SBAR)
      TCHK=.FALSE.
      DO IIB=1,NNS
        IF(SBAR(IIB)%IAT1.NE.IAT) CYCLE
        IF(SBAR(IIB)%IAT2.NE.IAT) CYCLE
        N1=SBAR(IIB)%N1
        ALLOCATE(SBARONSITE(N1,N1)) 
        SBARONSITE(:,:)=SBAR(IIB)%MAT
        TCHK=.TRUE.
        EXIT
      ENDDO
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('ONSITE TERM OF STRUCTURE CONSTANTS NOT FOUND')
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$STOP('LMTO$DOLOCORB')
      END IF
!
!     == SPHERICAL AVERAGE =====================================================
      LMX=(MAXVAL(LOX)+1)**2
      DO I=1,LMX
        DO J=1,LMX
          IF(I.EQ.J) CYCLE
          SBARONSITE(I,J)=0.D0
        ENDDO
      ENDDO
      DO L=0,MAXVAL(LOX)
        SVAR=0.D0
        DO LM=L**2+1,(L+1)**2
          SVAR=SVAR+SBARONSITE(LM,LM)
        ENDDO
        SVAR=SVAR/REAL(2*L+1,KIND=8)
        DO LM=L**2+1,(L+1)**2
          SBARONSITE(LM,LM)=SVAR
        ENDDO
PRINT*,'SBARONSITE ',L,SBARONSITE(L**2+1,L**2+1)
      ENDDO
!
!     ==========================================================================
!     == COUNT ONSITE ORBITALS BEFORE EXCLUSION                               ==
!     ==========================================================================
      LNXCHI1=0
      DO LN=1,LNX
        IF(ISCATT(LN).LE.0) LNXCHI1=LNXCHI1+1
      ENDDO
      ALLOCATE(LOXCHI(LNXCHI1))
      ALLOCATE(AECHI(NR,LNXCHI1))
      ALLOCATE(PSCHI(NR,LNXCHI1))
      ALLOCATE(NLCHI(NR,LNXCHI1))
!
!     ==========================================================================
!     == CONSTRUCT LOCAL ORBITALS                                             ==
!     ==========================================================================
      LNCHI=0
      DO LN=1,LNX
        IF(ISCATT(LN).GT.0) CYCLE
        LNCHI=LNCHI+1
        L=LOX(LN)
        LOXCHI(LNCHI)=L
        AECHI(:,LNCHI)=AEPHI(:,LN)          
        PSCHI(:,LNCHI)=PSPHI(:,LN)          
        NLCHI(:,LNCHI)=NLPHI(:,LN)          
        LM=L**2+1
!        SVAR=POTPAR(ISP)%CBAR(LN)-POTPAR(ISP)%ENU(LN) &
!     &        +SBARONSITE(LM,LM)*POTPAR(ISP)%SQDELTABAR(LN)**2
        SVAR=POTPAR(ISP)%KTOPHIDOT(LN) &
     &      +POTPAR(ISP)%JBARTOPHIDOT(LN)*SBARONSITE(LM,LM)
        SVAR=SVAR/POTPAR(ISP)%KTOPHI(LN)
PRINT*,'LN,LNCHI,L, SVAR ',LN,LNCHI,L,SVAR,AEPHIDOT(IRCOV,LN),AECHI(IRCOV,LNCHI)
!IF(AEPHIDOT(IRCOV,LN)*SVAR/AECHI(IRCOV,LNCHI).GT.0.D0) CYCLE
!        IF(SVAR.LT.0.D0) CYCLE
        AECHI(:,LNCHI)=AECHI(:,LNCHI)+AEPHIDOT(:,LN)*SVAR
        PSCHI(:,LNCHI)=PSCHI(:,LNCHI)+PSPHIDOT(:,LN)*SVAR
        NLCHI(:,LNCHI)=NLCHI(:,LNCHI)+NLPHIDOT(:,LN)*SVAR
      ENDDO
PRINT*,'LNX,LOX ',LNX,LOX
PRINT*,'LNXCHI1,LOXCHI ',LNXCHI1,LOXCHI
!
!     == ATTACH EXPONENTIAL TAIL AT THE COVALENT RADIUS ========================
      DO LNCHI=1,LNXCHI1     
        L=LOXCHI(LNCHI)
        CALL RADIAL$VALUE(GID,NR,NLCHI(:,LNCHI),RCOV,VAL)
        CALL RADIAL$DERIVATIVE(GID,NR,NLCHI(:,LNCHI),RCOV,DER)
        CALL RADIAL$VALUE(GID,NR,AECHI(:,LNCHI),RCOV,VALD)
        CALL RADIAL$DERIVATIVE(GID,NR,AECHI(:,LNCHI),RCOV,DERD)
        VALD=VALD-VAL
        DERD=DERD-DER
!
        SVAR=DER/VAL
        IF(ABS(VALD/VAL).GT.1.D-9) THEN ! AVOID DIVIDE-BY-ZERO
          SVAR1=DERD/VALD
        ELSE
          SVAR1=-1.D0
        END IF        
! 
        IF(SVAR.GT.0.D0.OR.SVAR1.GT.0.D0) THEN
          CALL SETUP_WRITEPHI(-'FAILEDAEORBITAL.DAT',GID,NR,LNXCHI1,AECHI)
          CALL SETUP_WRITEPHI(-'FAILEDPSORBITAL.DAT',GID,NR,LNXCHI1,PSCHI)
          CALL SETUP_WRITEPHI(-'FAILEDNLORBITAL.DAT',GID,NR,LNXCHI1,NLCHI)
          CALL SETUP_WRITEPHI(-'FAILEDDIFFORBITAL.DAT',GID,NR,LNXCHI1,AECHI-NLCHI)
          CALL ERROR$MSG('ORBITAL DOES NOT DECAY')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$I4VAL('ISP',ISP)
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('LNCHI',LNCHI)
          CALL ERROR$I4VAL('LOXCHI',LOXCHI(LNCHI))
          CALL ERROR$R8VAL('COVALENT RADIUS',RCOV)
          CALL ERROR$R8VAL('LOGARITMIC DERIVATIVE ',SVAR)
          CALL ERROR$R8VAL('LOGARITMIC DERIVATIVE OF DIFFERENCE',SVAR1)
          CALL ERROR$STOP('LMTO$DOLOCORB')
        END IF
        DO IR=IRCOV,NR
          NLCHI(IR,LNCHI)=VAL*EXP(SVAR*(R(IR)-RCOV))
          AECHI(IR,LNCHI)=VAL*EXP(SVAR*(R(IR)-RCOV)) &
                         +VALD*EXP(SVAR1*(R(IR)-RCOV))
          PSCHI(IR,LNCHI)=VAL*EXP(SVAR*(R(IR)-RCOV)) &
                         +VALD*EXP(SVAR1*(R(IR)-RCOV))
        ENDDO
      ENDDO
!
!     ==ORTHONORMALIZE LOCAL ORBITALS ==========================================
      DO LNCHI=1,LNXCHI1
        L=LOXCHI(LNCHI)
        DO LN=1,LNCHI-1
          IF(LOXCHI(LN).NE.L) CYCLE
          AUX(:)=R(:)**2*AECHI(:,LNCHI)*AECHI(:,LN)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
          AECHI(:,LNCHI)=AECHI(:,LNCHI)-AECHI(:,LN)*VAL
          PSCHI(:,LNCHI)=PSCHI(:,LNCHI)-PSCHI(:,LN)*VAL
        ENDDO
        AUX(:)=R(:)**2*AECHI(:,LNCHI)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
        VAL=1.D0/SQRT(VAL)
        AECHI(:,LNCHI)=AECHI(:,LNCHI)*VAL
        PSCHI(:,LNCHI)=PSCHI(:,LNCHI)*VAL
      ENDDO
!
!     ==========================================================================
!     == TRANSFORMATION OF PROJECTORS                                         ==
!     ==========================================================================
      ALLOCATE(AMAT(LNXCHI1,LNXPHI))
      ALLOCATE(BMAT(LNXCHI1,LNXCHI1))
      ALLOCATE(XMAT(LNXPHI,LNXCHI1))
      AMAT(:,:)=0.D0
      DO LNCHI=1,LNXCHI1
        DO LNPHI=1,LNXPHI
          IF(LOX(LNPHI).NE.LOXCHI(LNCHI)) CYCLE
          AUX(:)=R(:)**2*PRO(:,LNPHI)*PSCHI(:,LNCHI)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,AMAT(LNCHI,LNPHI))
        ENDDO
      ENDDO
      BMAT(:,:)=0.D0
      DO LNCHI=1,LNXCHI1
        BMAT(LNCHI,LNCHI)=1.D0
      ENDDO
      CALL LIB$MATRIXSOLVER8(LNXCHI1,LNXPHI,LNXCHI1,AMAT,XMAT,BMAT)
      AMAT=TRANSPOSE(XMAT)
      DEALLOCATE(XMAT)
      DEALLOCATE(BMAT)
!
!     ==========================================================================
!     == DELETE ORBITALS NOT IN THE SET                                       ==
!     ==========================================================================
      LNCHI=0
      LNCHI1=0
      DO LN=1,LNX
        IF(ISCATT(LN).GT.0) CYCLE
        LNCHI1=LNCHI1+1
        IF(.NOT.TORB(LN)) CYCLE
        LNCHI=LNCHI+1
        CHI(:,LNCHI)=AECHI(:,LNCHI1)
        CHIPHI(LNCHI,:)=AMAT(LNCHI1,:)
PRINT*,'CHIPHI ',CHIPHI(LNCHI,:)
      ENDDO
!
!     ==========================================================================
!     == PLOT LOCAL ORBITALS                                                  ==
!     ==========================================================================
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
      LOGICAL(4) ,PARAMETER :: T2D=.FALSE.
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
      &                         ,ORBG(N1D*(IAT-1)+I,LM1),IAT=1,NATCLUSTER-1)
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
 
