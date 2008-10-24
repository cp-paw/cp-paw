MODULE LMTO_MODULE
REAL(8)   ,PARAMETER  :: K2=0.D0
REAL(8)   ,PARAMETER  :: RC=5.D0  ! CUTOF RADIUS FOR NEIGHBORLIST
TYPE POTPAR_TYPE
   REAL(8) :: QBAR
   REAL(8) :: WKU    ! W[K,U]
   REAL(8) :: WKQ    ! W[K,Q]
   REAL(8) :: WJU    ! W[J,U]
   REAL(8) :: WJQ    ! W[J,Q]
   REAL(8) :: WUQ    ! W[U,Q]
   REAL(8) :: WKJ    ! W[K,J]
   REAL(8) :: KJTOUQ(2,2)    ! (K,J)=(U,Q)*KJTOUQ
   REAL(8) :: KJBARTOUQ(2,2) ! (K,JBAR)=(U,Q)*KJTOUQ
   REAL(8) :: JBARTOQ  ! JBAR=Q*JBARTOQ
   REAL(8) :: VALK     ! K(RAD)
   REAL(8) :: DERK     ! DK/DR(RAD)
   REAL(8) :: VALJ     ! J(RAD)
   REAL(8) :: DERJ     ! DJ/DR(RAD)
   REAL(8) :: VALJBAR  ! JBAR(RAD)
   REAL(8) :: DERJBAR  ! DJBAR/DR(RAD)
   REAL(8) :: VALKBAR  ! KBAR(RAD)
   REAL(8) :: DERKBAR  ! DKBAR/DR(RAD)
   REAL(8) :: OVUU     ! <U|U>
   REAL(8) :: OVUQ     ! <U|Q>
   REAL(8) :: OVQQ     ! <Q|Q>
END TYPE POTPAR_TYPE
TYPE PERIODICMAT_TYPE
INTEGER(4) :: IAT1
INTEGER(4) :: IAT2
INTEGER(4) :: IT(3)
INTEGER(4) :: N1
INTEGER(4) :: N2
REAL(8),POINTER :: MAT(:,:)
END TYPE PERIODICMAT_TYPE
LOGICAL(4)            :: TINIPOTPAR=.FALSE.
LOGICAL(4)            :: TINISTRUC=.FALSE.
INTEGER(4)            :: NSP
INTEGER(4)            :: LXX
INTEGER(4),ALLOCATABLE:: LX(:)               !(NSP)
REAL(8)   ,ALLOCATABLE:: RAD(:)              !(NSP)
integer(4),ALLOCATABLE:: ISPECIES(:)         !(NAT)
TYPE(POTPAR_TYPE),ALLOCATABLE:: POTPAR(:,:)  !(LXX+1,NSP)
TYPE(PERIODICMAT_TYPE),ALLOCATABLE:: SBAR(:) !(NNB)
END MODULE LMTO_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$REPORTSBAR(NFIL)
!     **************************************************************************      
!     **                                                                      **
!     **************************************************************************      
      USE LMTO_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NNB
      INTEGER(4)            :: NN,LM1
!     **************************************************************************      
      NNB=SIZE(SBAR)
      DO NN=1,NNB
        WRITE(NFIL,FMT='(82("="),T10," IAT1=",I5," IAT2=",I5," IT=",3I3," ")') &
     &                 SBAR(NN)%IAT1,SBAR(NN)%IAT2,SBAR(NN)%IT
        DO LM1=1,SBAR(NN)%N1
          WRITE(NFIL,FMT='(20F10.5)')SBAR(NN)%MAT(LM1,:)
        ENDDO
        WRITE(NFIL,FMT='(82("="))')
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
!!$      end TYPE OVKJ_TYPE
!!$      INTEGER(4),INTENT(IN) :: NFIL
!!$      INTEGER(4)            :: NNB
!!$      REAL(8)               :: A(2,2)
!!$      TYPE(OVKJ_TYPE)       :: OVKJ(NSP)
!!$      INTEGER(4)            :: ISVAR
!!$      INTEGER(4)            :: NN,nn1,nn2,LM1,isp,i
!!$      INTEGER(4)            :: NNb
!!$!     **************************************************************************      
!!$      DO ISP=1,NSP
!!$        ISVAR=(LX(ISP)+1)**2
!!$        ALLOCATE(OvKJ(ISP)%KK(ISVAR))
!!$        ALLOCATE(OvKJ(ISP)%KJBAR(ISVAR))
!!$        ALLOCATE(OvKJ(ISP)%JBARJBAR(ISVAR))
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
!!$      NNBs=SIZE(SBAR)
!!$!
!!$!     ===========================================================================
!!$!     == ONSITE TERMS                                                          ==
!!$!     ===========================================================================
!!$      DO NN=1,NNB
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
!!$      DO NN=1,NNB
!!$        IAT1=OVERLAP(NN)%IAT1
!!$        IAT2=OVERLAP(NN)%IAT2
!!$        IT(:)=OVERLAP(NN)%IT(:)
!!$        ISP1=ISPECIES(IAT1)
!!$        ISP1=ISPECIES(IAT2)
!!$        DO NN1=1,NNBS
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
!!$        DO NN1=1,NNBS
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
!!$      DO NN=1,NNB
!!$        IAT1=OVERLAP(NN)%IAT1
!!$        IAT2=OVERLAP(NN)%IAT2
!!$        IT(:)=OVERLAP(NN)%IT(:)
!!$        ISP1=ISPECIES(IAT1)
!!$        ISP1=ISPECIES(IAT2)
!!$        DO NN1=1,NNBS
!!$          IF(IAT1.NE.SBAR(NN1)%IAT1) CYCLE
!!$          IT2(:)=SBAR(NN1)%IT
!!$          DO NN2=1,NNBS
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
      SUBROUTINE LMTO$MAKEPOTENTIALPARAMETERS
!     **************************************************************************      
!     **                                                                      **
!     **************************************************************************      
      USE LMTO_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4),PARAMETER   :: NBX=19
      INTEGER(4)             :: NB
      INTEGER(4)             :: LOFI(NBX)
      INTEGER(4)             :: SOFI(NBX)
      INTEGER(4)             :: NNOFI(NBX)
      REAL(8)                :: FOFI(NBX)
      REAL(8)                :: EOFI(NBX)
      REAL(8)   ,ALLOCATABLE :: POT(:)
      REAL(8)   ,ALLOCATABLE :: PHI(:,:)
      REAL(8)   ,ALLOCATABLE :: UOFI(:,:)
      REAL(8)   ,ALLOCATABLE :: TUOFI(:,:)
      REAL(8)   ,ALLOCATABLE :: UN(:,:)
      REAL(8)   ,ALLOCATABLE :: QNPLUS1(:,:)
      REAL(8)   ,ALLOCATABLE :: PHINU(:,:)
      REAL(8)   ,ALLOCATABLE :: PHINUDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)   ,ALLOCATABLE :: G(:)
      REAL(8)   ,ALLOCATABLE :: DREL(:)
      REAL(8)                :: AEZ
      REAL(8)                :: ETOT
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      INTEGER(4)             :: ISP,IB,L,IR
      INTEGER(4)             :: IC,IV
      INTEGER(4)             :: NC
      INTEGER(4)             :: IROUT
      REAL(8)                :: RBOX
      REAL(8)                :: VALU,DERU,VALQ,DERQ
      REAL(8)                :: VALK,DERK,VALJ,DERJ
      REAL(8)                :: WKU,WJU,WKQ,WJQ,WKJ,WUQ
      REAL(8)                :: CUK,CQK,CUJ,CQJ,CUJBAR,CQJBAR
      REAL(8)                :: E
      REAL(8)                :: QBAR
      REAL(8)                :: VAL
      REAL(8)   ,ALLOCATABLE :: AUX(:),AUX1(:)
      CHARACTER(64),PARAMETER :: KEY='START,REL,NONSO'
!     **************************************************************************      
      IF(TINIPOTPAR) RETURN
      TINIPOTPAR=.TRUE.
!
      CALL SETUP$NSPECIES(NSP)
      ALLOCATE(RAD(NSP))
      ALLOCATE(LX(NSP))
      LXX=-1
      DO ISP=1,NSP
        CALL SETUP$LNX(ISP,LNX)
        ALLOCATE(LOX(LNX))
        CALL SETUP$LOFLN(ISP,LNX,LOX)
        LX(ISP)=MAXVAL(LOX)
        DEALLOCATE(LOX)
        CALL SETUP$AEZ(ISP,AEZ)
        CALL PERIODICTABLE$GET(NINT(AEZ),'R(ASA)',RAD(ISP))
      ENDDO
      LXX=MAXVAL(LX)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ALLOCATE(POTPAR(LXX+1,NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(POT(NR))
        ALLOCATE(PHI(NR,NBX))
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        RBOX=R(NR-3)
        CALL SETUP$AEZ(ISP,AEZ)
        CALL ATOMLIB$AESCF(GID,NR,KEY,RBOX,AEZ,NBX,NB,LOFI,SOFI,FOFI,NNOFI &
    &                     ,ETOT,POT,EOFI,PHI)
!
!       == PRINT ==============================================================
        WRITE(*,FMT='(82("="),T20,"ATOM TYPE ",I3," Z="F5.1)')ISP,AEZ
        DO IB=1,NB
          WRITE(*,FMT='(3I10,2F15.5)')LOFI(IB),SOFI(IB),NNOFI(IB),FOFI(IB),EOFI(IB)        
        ENDDO
!
!       ========================================================================
!       ==  NODELESS WAVE FUNCTIONS                                           ==
!       ========================================================================
        ALLOCATE(UOFI(NR,NBX))
        ALLOCATE(TUOFI(NR,NBX))
        CALL ATOMLIB$NODELESS(GID,NR,RBOX,POT,NB,LOFI,EOFI,UOFI,TUOFI)
!
!       == PRINT ==============================================================
        WRITE(*,FMT='(82("="),T20,"ATOM TYPE ",I3," Z="F5.1)')ISP,AEZ
        DO IB=1,NB
          WRITE(*,FMT='(3I10,2F15.5)')LOFI(IB),SOFI(IB),NNOFI(IB),FOFI(IB),EOFI(IB)        
        ENDDO
!
!       ========================================================================
!       ==  UN,Q_(N+1)                                                        ==
!       ========================================================================
        ALLOCATE(UN(NR,LX(ISP)+1))
        ALLOCATE(QNPLUS1(NR,LX(ISP)+1))
        ALLOCATE(PHINU(NR,LX(ISP)+1))
        ALLOCATE(PHINUDOT(NR,LX(ISP)+1))
        ALLOCATE(G(NR))
        ALLOCATE(DREL(NR))
        ALLOCATE(AUX(NR))
        ALLOCATE(AUX1(NR))
        CALL SETUP$GETI4('NC',NC)
        DO L=0,LX(ISP)
!         == FIND VALENCE STATE FOR THIS L =====================================
          IV=0
          DO IB=NC+1,NB
            IF(LOFI(IB).EQ.L) THEN
              IV=IB 
              EXIT
            END IF
          ENDDO
!         == FIND CORE STATE FOR THIS L ========================================
          IC=0
          DO IB=NC,1,-1
            IF(LOFI(IB).EQ.L) THEN
              IC=IB 
              EXIT
            END IF
          ENDDO
!
!         == MAP VALENCE NODELESS FUNCTION INTO ARRAY ==========================
          IF(IV.NE.0) THEN
            E=EOFI(IV)
            UN(:,L+1)=UOFI(:,IV)
          ELSE
            E=0.D0
            G(:)=0.D0
            IF(IC.NE.0) G(:)=UOFI(:,IC)
            DREL(:)=0.D0
            CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,0,G,L,E,1,UN(:,L+1))
          END IF
!
!         == DETERMINE SCATTERING NODELESS WAVE FUNCTION =======================
          G(:)=UN(:,L+1)
          DREL(:)=0.D0
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT,DREL,0,G,L,E,1,QNPLUS1(:,L+1))
!
!         == CORE-ORTHOGONALIZE NODELESS SCATTERING WAVE FUNCTIONS==============
          DO IR=1,NR
           IF(R(IR).LE.5.D0) CYCLE
           IROUT=IR
           EXIT
          ENDDO
          PHINU(:,L+1)=UN(:,L+1)
          PHINUDOT(:,L+1)=QNPLUS1(:,L+1)
          DO IB=NC,1,-1
            AUX(:)=R(:)**2*PHINU(:,L+1)*UOFI(:,IB)
            AUX(IROUT:)=0.D0    ! AVOID INFTY*0
            CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
            PHINU(:,L+1)=PHINU(:,L+1)-UOFI(:,IB)*VAL
            AUX(:)=R(:)**2*PHINUDOT(:,L+1)*UOFI(:,IB)
            AUX(IROUT:)=0.D0    ! AVOID INFTY*0
            CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
            PHINUDOT(:,L+1)=PHINUDOT(:,L+1)-UOFI(:,IB)*VAL
          ENDDO
!
!         == NORMALIZE PHINU ===================================================
          AUX(:)=R(:)**2*PHINU(:,L+1)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX,RAD(ISP),VAL)
          VAL=1.D0/SQRT(VAL)
          PHINU(:,L+1)=PHINU(:,L+1)*VAL
          PHINUDOT(:,L+1)=PHINUDOT(:,L+1)*VAL
          UN(:,L+1)=UN(:,L+1)*VAL
          QNPLUS1(:,L+1)=QNPLUS1(:,L+1)*VAL
        ENDDO

        DEALLOCATE(G)
        DEALLOCATE(DREL)
!
!       ========================================================================
!       ==  DETERMINE WRONSKI MATRICES AND QBAR                               ==
!       ========================================================================
        DO L=0,LX(ISP)
          CALL RADIAL$VALUE(GID,NR,UN(:,L+1),RAD(ISP),VALU)
          CALL RADIAL$DERIVATIVE(GID,NR,UN(:,L+1),RAD(ISP),DERU)
          CALL RADIAL$VALUE(GID,NR,QNPLUS1(:,L+1),RAD(ISP),VALQ)
          CALL RADIAL$DERIVATIVE(GID,NR,QNPLUS1(:,L+1),RAD(ISP),DERQ)
          CALL LMTO$SOLIDBESSELRAD(L,RAD(ISP),K2,VALJ,DERJ)
          CALL LMTO$SOLIDHANKELRAD(L,RAD(ISP),K2,VALK,DERK)
          WKU=VALK*DERU-DERK*VALU
          WKQ=VALK*DERQ-DERK*VALQ
          WJU=VALJ*DERU-DERJ*VALU
          WJQ=VALJ*DERQ-DERJ*VALQ
          WKJ=VALK*DERJ-DERK*VALJ
!PRINT*,'WKJ ',WKJ,1.D0/RAD(ISP)**2,WKJ-1.D0/RAD(ISP)**2
          WUQ=VALU*DERQ-DERU*VALQ
!         == K -> U*CUK + Q*CQK ================================================
!         == J -> U*CUJ + Q*CQJ ================================================
!         == JBAR -> Q*CQJBAR ==================================================
          CUK=WKQ/WUQ
          CQK=WKU/WUQ
          CUJ=WJQ/WUQ
          CQJ=WJU/WUQ
          CALL LMTO$Q(L,RAD(ISP),VALQ,DERQ,K2,QBAR)
          CUJBAR=CUJ-QBAR*CUK
          CQJBAR=CQJ-QBAR*CQK
          
          POTPAR(L+1,ISP)%QBAR=QBAR
          POTPAR(L+1,ISP)%WKU=WKU
          POTPAR(L+1,ISP)%WKQ=WKQ
          POTPAR(L+1,ISP)%WJU=WJU
          POTPAR(L+1,ISP)%WJQ=WJQ
          POTPAR(L+1,ISP)%WKJ=WKJ
          POTPAR(L+1,ISP)%WUQ=WUQ
          POTPAR(L+1,ISP)%KJTOUQ(1,1)=CUK
          POTPAR(L+1,ISP)%KJTOUQ(1,2)=CQK
          POTPAR(L+1,ISP)%KJTOUQ(2,1)=CUJ
          POTPAR(L+1,ISP)%KJTOUQ(2,2)=CQJ
          POTPAR(L+1,ISP)%KJBARTOUQ(1,1)=CUK
          POTPAR(L+1,ISP)%KJBARTOUQ(1,2)=CQK
          POTPAR(L+1,ISP)%KJBARTOUQ(2,1)=0.D0
          POTPAR(L+1,ISP)%KJBARTOUQ(2,2)=CQJBAR
          POTPAR(L+1,ISP)%JBARTOQ=CQJBAR
          POTPAR(L+1,ISP)%VALK=VALK
          POTPAR(L+1,ISP)%DERK=DERK
          POTPAR(L+1,ISP)%VALJ=VALJ
          POTPAR(L+1,ISP)%DERJ=DERJ
          POTPAR(L+1,ISP)%VALJBAR=VALJ-QBAR*VALK
          POTPAR(L+1,ISP)%DERJBAR=DERJ-QBAR*DERK
!
          AUX=R(:)**2*PHINU(:,L+1)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RAD(ISP),POTPAR(L+1,ISP)%OVUU)
          AUX=R(:)**2*PHINU(:,L+1)*PHINUDOT(:,L+1)
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RAD(ISP),POTPAR(L+1,ISP)%OVUQ)
          AUX=R(:)**2*PHINUDOT(:,L+1)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RAD(ISP),POTPAR(L+1,ISP)%OVQQ)
        ENDDO
!
!       ========================================================================
!       ==  CLEAN UP                                                          ==
!       ========================================================================
        DEALLOCATE(POT)
        DEALLOCATE(PHI)
        DEALLOCATE(UOFI)
        DEALLOCATE(TUOFI)
        DEALLOCATE(PHINU)
        DEALLOCATE(PHINUDOT)
        DEALLOCATE(UN)
        DEALLOCATE(QNPLUS1)
        DEALLOCATE(R)
        DEALLOCATE(AUX)
        DEALLOCATE(AUX1)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$MAKESTRUCTURECONSTANTS
!     **************************************************************************      
!     **                                                                      **
!     **************************************************************************      
      USE LMTO_MODULE
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: NAT       !#(ATOMS)
      REAL(8)                :: RBAS(3,3) !LATTICE VECTORS
      REAL(8)   ,ALLOCATABLE :: R0(:,:)   !(3,NAT) ATOMIC POSITIONS
      INTEGER(4),PARAMETER   :: NNXPERATOM=100
      INTEGER(4)             :: NNX
      INTEGER(4)             :: NNB
      INTEGER(4),ALLOCATABLE :: NNLIST(:,:) !(5,NNX)
      INTEGER(4)             :: NAT1
      INTEGER(4)             :: NORB
      INTEGER(4)             :: N
      INTEGER(4),ALLOCATABLE :: LX1(:)    !(NNB) MAX(ANGULAR MOMENTUM)
      REAL(8)   ,ALLOCATABLE :: RPOS(:,:) !(3,NNB(IAT)) ATOMIC POSITIONS
      REAL(8)   ,ALLOCATABLE :: QBAR(:)  !(N) SCREENING PARAMETER
      REAL(8)   ,ALLOCATABLE :: SBAR1(:,:) !
      REAL(8)   ,ALLOCATABLE :: C(:,:) !
      REAL(8)                :: VAL,DER
      INTEGER(4)             :: IAT,IAT1,IAT2,ISP,ISP1,ISP2,LMX1,LMX2 
      INTEGER(4)             :: L,NN,NN1,NN2,NN0,I,IM
!     **************************************************************************      
      IF(TINISTRUC) RETURN
      TINISTRUC=.TRUE.
      CALL LMTO$MAKEPOTENTIALPARAMETERS
!
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
!
!     ==========================================================================
!     == NEIGHBORLIST                                                         ==
!     ==========================================================================
      NNX=NNXPERATOM*NAT
      ALLOCATE(NNLIST(5,NNX))
      CALL LMTO$NEIGHBORLIST(RBAS,NAT,R0,RC,NNX,NNB,NNLIST)
!!$DO I=1,NNB
!!$WRITE(*,FMT='("NNLIST",5I10)')NNLIST(:,I)
!!$ENDDO
!
!     ==========================================================================
!     == STRUCTURE CONSTANTS                                                  ==
!     ==========================================================================
      ALLOCATE(SBAR(NNB))
      DO IAT1=1,NAT
!       == MEMBERS NN1:NN2 ON THE NEIGHBOLIST BILD THE CLUSTER AROUND ATOM 1  ==
!       == MEMBER NN0 IS THE ONSITE MEMBER                                    ==
        NN1=1
        NN0=0
        DO NN=1,NNB
          IF(NNLIST(1,NN).GT.IAT1)EXIT
          NN2=NN
          IF(NNLIST(1,NN).LT.IAT1)NN1=NN+1
          IF(NNLIST(1,NN).EQ.IAT1) THEN
            IF(NNLIST(2,NN).EQ.IAT1) THEN
              IF(NNLIST(3,NN).EQ.0.AND.NNLIST(4,NN).EQ.0 &
     &                            .AND.NNLIST(5,NN).EQ.0) THEN
                NN0=NN
              END IF
            END IF
          END IF
        ENDDO
!
        NAT1=NN2-NN1+1  ! #(ATOMS ON THE CLUSTER )
        ALLOCATE(LX1(NAT1))
        ALLOCATE(RPOS(3,NAT1))
        DO NN=NN1,NN2
          IAT2=NNLIST(2,NN)
          ISP=ISPECIES(IAT2)
          LX1(NN-NN1+1)=LX(ISP)
          RPOS(:,NN-NN1+1)=R0(:,IAT2)+RBAS(:,1)*REAL(NNLIST(3,NN),KIND=8) &
     &                               +RBAS(:,2)*REAL(NNLIST(4,NN),KIND=8) &
     &                               +RBAS(:,3)*REAL(NNLIST(5,NN),KIND=8)
        ENDDO
        NORB=(LX1(NN0-NN1+1)+1)**2
        N=SUM((LX1(:)+1)**2)
        ALLOCATE(QBAR(N))
        I=0
        DO NN=NN1,NN2
          IAT2=NNLIST(2,NN)
          ISP=ISPECIES(IAT2)
!WRITE(*,FMT='("IAT1=",I5," DR=",3F10.4)')IAT1,RPOS(:,NN-NN1+1)-RPOS(:,1)
          DO L=0,LX(ISP)
            DO IM=1,2*L+1
              I=I+1
              QBAR(I)=POTPAR(L+1,ISP)%QBAR
            ENDDO
          ENDDO
        ENDDO
        ALLOCATE(SBAR1(N,NORB))
        ALLOCATE(C(N,NORB))
        CALL LMTO$CLUSTERSTRUCTURECONSTANTS(K2,NAT1,RPOS,LX1,QBAR,N,NORB,SBAR1)
        CALL LMTO$KBARMULTICENTER(N,NORB,QBAR,SBAR1,C)
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
!WRITE(*,FMT='("IAT1=",I5," DR=",3F10.4)')IAT1,RPOS(:,NN-NN1+1)-RPOS(:,1)
          ISP1=ISPECIES(IAT1)
          ISP2=ISPECIES(IAT2)
          LMX1=(LX(ISP1)+1)**2      
          LMX2=(LX(ISP2)+1)**2      
          SBAR(NN)%N1=LMX1
          SBAR(NN)%N2=LMX2
          ALLOCATE(SBAR(NN)%MAT(LMX1,LMX2))
          SBAR(NN)%MAT(:,:)=TRANSPOSE(SBAR1(I+1:I+LMX2,:))
          I=I+LMX2
        ENDDO
!
!       ========================================================================
!       == DETERMINE "KBAR", I.E. K_I-JBAR_I*SBAR_{II}                        ==
!       ========================================================================
        ISP=ISPECIES(IAT1)
        I=0
        DO L=0,LX(ISP)
          VAL=0.D0
          DER=0.D0
          DO IM=1,2*L+1
            I=I+1
            VAL=VAL+POTPAR(L+1,ISP)%VALK-POTPAR(L+1,ISP)%VALJBAR*SBAR1(I,I)
            DER=DER+POTPAR(L+1,ISP)%DERK-POTPAR(L+1,ISP)%DERJBAR*SBAR1(I,I)
          ENDDO
          POTPAR(L+1,ISP)%VALKBAR=VAL/REAL(2*L+1,KIND=8)
          POTPAR(L+1,ISP)%DERKBAR=DER/REAL(2*L+1,KIND=8)
        ENDDO
        DEALLOCATE(LX1)
        DEALLOCATE(RPOS)
        DEALLOCATE(QBAR)
        DEALLOCATE(SBAR1)
        DEALLOCATE(C)
      ENDDO
      CALL LMTO$REPORTSBAR(6)
      RETURN
      END      
