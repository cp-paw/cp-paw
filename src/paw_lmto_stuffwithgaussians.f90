MODULE LMTO_stuff_MODULE
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
!
TYPE ORBITALSPHHARM_TYPE
  INTEGER(4)         :: GID
  INTEGER(4)         :: NR
  INTEGER(4)         :: LMX
  INTEGER(4)         :: NORB
  REAL(8)   ,POINTER :: F(:,:,:) !(NR,LM,IORB)
END TYPE ORBITALSPHHARM_TYPE
!
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
TYPE POTPARg1_TYPE
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
  TYPE(TAILED_TYPE)  :: TAILED
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
LOGICAL(4)            :: TOFFSITE=.FALSE.  !INCLUDE OFFSITE EXCHANGE
LOGICAL(4)            :: TDROP=.FALSE. ! WRITE THE WAVE FUNCTIONS TO FILE
LOGICAL(4)            :: TPICK=.FALSE. ! REAL HAMILTON CORRECTION FROM FILE

REAL(8)               :: K2=-0.25D0    ! 0.5*K2 IS THE KINETIC ENERGY
REAL(8)               :: RCSCALE=1.2D0  !RADIUS SCALE FACTOR FOR NEIGHBORLIST
! RCSCALE=5. IS GOOD FOR THE HUBBARD MODEL WITH LATTICE CONSTANT=3\AA
!REAL(8)               :: RCSCALE=5.D0  !RADIUS SCALE FACTOR FOR NEIGHBORLIST
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
CHARACTER(32)           :: MODUS='NONE'
INTEGER(4)              :: NSP=-1
INTEGER(4)              :: ISPSELECTOR=-1 ! USED ONLY FOR HYBRIDSETTING
TYPE(HYBRIDSETTING_TYPE),ALLOCATABLE :: HYBRIDSETTING(:)
INTEGER(4),ALLOCATABLE  :: LNX(:)      !(NSP)
INTEGER(4),ALLOCATABLE  :: LOX(:,:)    !(LNXX,NSP)
INTEGER(4),ALLOCATABLE  :: ISPECIES(:) !(NAT)
REAL(8)   ,ALLOCATABLE  :: ORBRAD(:,:) !(LXX+1,NAT) NODE-POSITION OF THE ORBITAL
TYPE(POTPAR_TYPE)     ,ALLOCATABLE :: POTPAR(:) !POTENTIAL PARAMETERS
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
!===============================================================================
!=====  STRUCTURE DEPENDENT DATA  ==============================================
!===============================================================================
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: SBAR(:)    !(NNS) SCREEND STRUCTURE CONST.
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: SBAR_NEW(:)!(NNS) SCREENED STRUCT. CONST.
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: PCHI(:)    !(NNS) <PTILDE|CHITILDE>
!TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DENMAT(:)  !(NND) DENSITY MATRIX
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DENMAT_NEW(:)  !(NND) DENSITY MATRIX
TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: HAMIL_NEW(:)   !(NND) DERIVATIVE OF ENERGY
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: OVERLAP(:) !(NNS) OVERLAP MATRIX ONLY MAIN
!INTEGER(4)                         :: NNU       !#(ELEMENTS OF UMAT)
!TYPE(UMAT_TYPE)       ,ALLOCATABLE :: UMAT(:)   !(NNUX/NNU) U-MATRIX ELEMENTS
!!$TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: DENMAT_T(:) !(NNS) DENSITY MATRIX
!!$TYPE(PERIODICMAT2_TYPE),ALLOCATABLE:: HAMIL_T(:)  !(NNS) DERIVATIVE OF ENERGY
!!$INTEGER(4)                         :: NNUX        !DIMENSION OF UMAT
END MODULE LMTO_stuff_MODULE
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
      IMPLICIT NONE
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
