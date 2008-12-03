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
   REAL(8) :: RAD      ! RADIUS WHERE VALUE AND DERIVATIVE ARE TAKEN
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
INTEGER(4),ALLOCATABLE:: ISPECIES(:)         !(NAT)
TYPE(POTPAR_TYPE),ALLOCATABLE:: POTPAR(:,:)  !(LXX+1,NSP)
TYPE(PERIODICMAT_TYPE),ALLOCATABLE:: SBAR(:) !(NNB)
REAL(8)   ,ALLOCATABLE :: ORBRAD(:,:) !(LXX+1,NAT) NODE-POSITION OF THE ORBITAL
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
!!$      END TYPE OVKJ_TYPE
!!$      INTEGER(4),INTENT(IN) :: NFIL
!!$      INTEGER(4)            :: NNB
!!$      REAL(8)               :: A(2,2)
!!$      TYPE(OVKJ_TYPE)       :: OVKJ(NSP)
!!$      INTEGER(4)            :: ISVAR
!!$      INTEGER(4)            :: NN,NN1,NN2,LM1,ISP,I
!!$      INTEGER(4)            :: NNB
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
!!$      NNBS=SIZE(SBAR)
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
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: UPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: UN(:,:)
      REAL(8)   ,ALLOCATABLE :: QNPLUS1(:,:)
      REAL(8)   ,ALLOCATABLE :: PHINU(:,:)
      REAL(8)   ,ALLOCATABLE :: PHINUDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      INTEGER(4),ALLOCATABLE :: ISCATT(:)
      LOGICAL(4),ALLOCATABLE :: TPHI(:)
      LOGICAL(4),ALLOCATABLE :: TPHIDOT(:)
      REAL(8)                :: AEZ
      REAL(8)                :: ETOT
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      INTEGER(4)             :: ISP,IB,L,IR,LN
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
      logical(4)             :: tchk
!     **************************************************************************
      IF(TINIPOTPAR) RETURN
      TINIPOTPAR=.TRUE.
                                 call trace$push('LMTO$MAKEPOTENTIALPARAMETERS')
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
!       == POTENTIAL PARAMETERS CAN ONLY BE DETERMINED WITH INTERNAL SETUPS ====
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(AUX(NR))
        ALLOCATE(AUX1(NR))
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
!
!       ========================================================================
!       ==  UN,Q_(N+1)                                                        ==
!       ========================================================================
        ALLOCATE(UN(NR,LX(ISP)+1))
        ALLOCATE(QNPLUS1(NR,LX(ISP)+1))
        ALLOCATE(PHINU(NR,LX(ISP)+1))
        ALLOCATE(PHINUDOT(NR,LX(ISP)+1))
!
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(ISCATT(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
        ALLOCATE(AEPHI(NR,LNX))
        ALLOCATE(UPHI(NR,LNX))
        CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
        CALL SETUP$GETR8A('NDLSPHI',NR*LNX,UPHI)
        ALLOCATE(TPHI(LX(ISP)+1))
        ALLOCATE(TPHIDOT(LX(ISP)+1))
        TPHI(:)=.FALSE.
        TPHIDOT(:)=.FALSE.
        DO LN=1,LNX
          L=LOX(LN) 
          IF(ISCATT(LN).EQ.0) THEN
            TPHI(L+1)=.TRUE.
            UN(:,L+1)=UPHI(:,LN)
            PHINU(:,L+1)=AEPHI(:,LN)
          ELSE IF(ISCATT(LN).EQ.1) THEN
            TPHIDOT(L+1)=.TRUE.
            QNPLUS1(:,L+1)=UPHI(:,LN)
            PHINUDOT(:,L+1)=AEPHI(:,LN)
          END IF
        ENDDO
!
!       == NORMALIZE PHINU =====================================================
        DO L=0,LX(ISP)
          IF(TPHI(L+1)) THEN
            IF(.NOT.TPHIDOT(L+1)) THEN
              DO IR=1,NR
                CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,QNPLUS1(IR,L+1),DERJ)
              ENDDO
              PHINUDOT(:,L+1)=QNPLUS1(:,L+1)
            END IF
          ELSE
            IF(TPHIDOT(L+1)) THEN
              CALL ERROR$msg('tphidot=T')
print*,'lox ',lox
print*,'iscatt ',iscatt
              CALL ERROR$i4val('l',l)
              CALL ERROR$STOP('LMTO$MAKEPOTENTIALPARAMETERS')
            END IF
            DO IR=1,NR
              CALL LMTO$SOLIDHANKELRAD(L,R(IR),K2,UN(IR,L+1),DERK)
              CALL LMTO$SOLIDBESSELRAD(L,R(IR),K2,QNPLUS1(IR,L+1),DERJ)
            ENDDO
            PHINU(:,L+1)=UN(:,L+1)
            PHINUDOT(:,L+1)=QNPLUS1(:,L+1)
          END IF
        ENDDO
!
!       == NORMALIZE PHINU =====================================================
        DO L=0,LX(ISP)
          IF(.NOT.TPHI(L+1)) CYCLE
          AUX(:)=R(:)**2*PHINU(:,L+1)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX,RAD(ISP),VAL)
          VAL=1.D0/SQRT(VAL)
          PHINU(:,L+1)=PHINU(:,L+1)*VAL
          PHINUDOT(:,L+1)=PHINUDOT(:,L+1)*VAL
          UN(:,L+1)=UN(:,L+1)*VAL
          QNPLUS1(:,L+1)=QNPLUS1(:,L+1)*VAL
        ENDDO
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
          
          POTPAR(L+1,ISP)%RAD=RAD(ISP)
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
        DEALLOCATE(TPHI)
        DEALLOCATE(TPHIDOT)
        DEALLOCATE(LOX)
        DEALLOCATE(ISCATT)
        DEALLOCATE(AEPHI)
        DEALLOCATE(UPHI)
        DEALLOCATE(PHINU)
        DEALLOCATE(PHINUDOT)
        DEALLOCATE(UN)
        DEALLOCATE(QNPLUS1)
        DEALLOCATE(R)
        DEALLOCATE(AUX)
        DEALLOCATE(AUX1)
      ENDDO
                                                                call trace$pop()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$ORBRAD
!     **************************************************************************      
!     **                                                                      **
!     **************************************************************************      
      USE LMTO_MODULE, ONLY : K2,LX,SBAR,ORBRAD
      USE PERIODICTABLE_MODULE
      REAL(8)   ,ALLOCATABLE      :: SBARDIAG(:)
      REAL(8)                     :: RAD1
      REAL(8)                     :: VALPHI,DERPHI
      REAL(8)                     :: VALK,DERK
      REAL(8)                     :: VALJ,DERJ
      REAL(8)                     :: AEZ
      INTEGER(4)                  :: LX1
      INTEGER(4)                  :: NNB
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
      LXX=MAXVAL(LX)
      ALLOCATE(QBAR(LXX+1,NSP))
      CALL LMTO_MAKEQBAR(NSP,LXX,K2,QBAR)
!
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
      NNB=SIZE(SBAR)
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
        DO INB=1,NNB
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
!     ** VAL AND DER ARE VALUE AND DERIVATIVE OF PHIDOT.                       **
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
      REAL(8)   ,ALLOCATABLE:: UPHI(:,:)  !(NR,LNX) NODELESS PARTIALWAVE
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
        ALLOCATE(UPHI(NR,LNX))
        CALL SETUP$GETR8A('NDLSPHI',NR*LNX,UPHI)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(NINT(AEZ),'R(ASA)',RAD)
        DO L=0,LXX
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            IF(ISCATT(LN).NE.1) CYCLE
            CALL RADIAL$VALUE(GID,NR,UPHI(:,LN),RAD,VAL)
            CALL RADIAL$DERIVATIVE(GID,NR,UPHI(:,LN),RAD,DER)
            CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
            CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
            QBAR(L+1,ISP)=(JVAL*DER-VAL*JDER)/(KVAL*DER-VAL*KDER)
          ENDDO
        ENDDO            
        DEALLOCATE(UPHI)
        DEALLOCATE(LOX)
        DEALLOCATE(ISCATT)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$MAKESTRUCTURECONSTANTS()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : K2,RC,LX,SBAR,POTPAR,TINISTRUC
      IMPLICIT NONE
      INTEGER(4)             :: NAT       !#(ATOMS)
      REAL(8)                :: RBAS(3,3) !LATTICE VECTORS
      REAL(8)   ,ALLOCATABLE :: R0(:,:)   !(3,NAT) ATOMIC POSITIONS
      INTEGER(4),PARAMETER   :: NNXPERATOM=100
      INTEGER(4)             :: NNX
      INTEGER(4),ALLOCATABLE :: NNLIST(:,:) !(5,NNX)
      INTEGER(4)             :: NAT1
      INTEGER(4)             :: NNB
      INTEGER(4)             :: NORB
      INTEGER(4)             :: N
      INTEGER(4),ALLOCATABLE :: LX1(:)     !(NNB) MAX(ANGULAR MOMENTUM)
      REAL(8)   ,ALLOCATABLE :: RPOS(:,:)  !(3,NNB(IAT)) ATOMIC POSITIONS
      REAL(8)   ,ALLOCATABLE :: QBAR1(:,:) !(LXX+1,NSP) SCREENING PARAMETER
      REAL(8)   ,ALLOCATABLE :: QBAR(:)    !(N) SCREENING PARAMETER
      REAL(8)   ,ALLOCATABLE :: SBAR1(:,:) !
      REAL(8)   ,ALLOCATABLE :: C(:,:) !
      REAL(8)                :: VAL,DER
      INTEGER(4)             :: IAT,IAT1,IAT2,ISP,ISP1,ISP2,LMX1,LMX2 
      INTEGER(4)             :: L,NN,NN1,NN2,NN0,I,IM
      INTEGER(4)             :: LXX        !X(ANGULAR MOMENTUM FOR ALL ATOMS)
      INTEGER(4)             :: NSP       !#(ATOM TYPES)
      INTEGER(4),ALLOCATABLE :: ISPECIES(:)
      logical(4)             :: tchk
!     **************************************************************************
      CALL SETUP$GETL4('INTERNALSETUPS',TCHK)
      IF(.NOT.TCHK) return
!
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
!     == DETERMINE SCREENING PARAMETER QBAR                                   ==
!     ==========================================================================
      CALL SETUP$NSPECIES(NSP)
      LXX=MAXVAL(LX)
      ALLOCATE(QBAR1(LXX+1,NSP))
      CALL LMTO_MAKEQBAR(NSP,LXX,K2,QBAR1)
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
!
!       ========================================================================
!       == EXPAND SCREENING PARAMETER QBAR                                    ==
!       ========================================================================
        ALLOCATE(QBAR(N))
        I=0
        DO NN=NN1,NN2
          IAT2=NNLIST(2,NN)
          ISP=ISPECIES(IAT2)
!WRITE(*,FMT='("IAT1=",I5," DR=",3F10.4)')IAT1,RPOS(:,NN-NN1+1)-RPOS(:,1)
          DO L=0,LX(ISP)
            DO IM=1,2*L+1
              I=I+1
              QBAR(I)=QBAR1(L+1,ISP)
!             QBAR(I)=POTPAR(L+1,ISP)%QBAR
            ENDDO
          ENDDO
        ENDDO
!
!       ========================================================================
!       == DETERMINE ASTRUCTURE CONSTANTS                                     ==
!       ========================================================================
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
        DEALLOCATE(LX1)
        DEALLOCATE(RPOS)
        DEALLOCATE(QBAR)
        DEALLOCATE(SBAR1)
        DEALLOCATE(C)
      ENDDO
      DEALLOCATE(QBAR1)
!
      CALL LMTO$REPORTSBAR(6)
!
      CALL LMTO$ORBRAD
      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$DOLOCORB(IAT,LNXCHI,LNXPHI,TORB,PHICHI,CHIPHI)
!     **************************************************************************
!     **  constructs onsite-mapping from partial waves to local orbitals      **
!     **                                                                      **
!     **      |chi_i> =  sum_j |phi_j> * Phichi_ji                            **
!     **                                                                      **
!     **  it uses orbrad, which is the radius obtained by linear extrapolation**
!     **  of the lmTO to zero. orbrad encodes the logarithmic derivative      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ORBRAD
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT     ! atom index
      INTEGER(4),INTENT(IN)  :: LNXCHI  ! #(local orbitals w/o m-multiplicity)
      INTEGER(4),INTENT(IN)  :: LNXPHI  ! #(partial waves w/o m-multiplicity)
      LOGICAL(4),INTENT(IN)  :: TORB(LNXPHI)          ! selects local orbitals
      REAL(8)   ,INTENT(OUT) :: PHICHI(LNXPHI,LNXCHI) !|CHI_I>=\SUM_J|PHI_J>*PHICHI_{JI}
      REAL(8)   ,INTENT(OUT) :: CHIPHI(LNXCHI,LNXPHI)
      INTEGER(4)             :: NAT         ! #(ATOMS)
      INTEGER(4),ALLOCATABLE :: ISPECIES(:) !(nat) atom type for each atom
      INTEGER(4)             :: ISP         ! ATOM TYPE INDEX
      REAL(8)                :: AEZ         ! ATOMIC NUMBER
      REAL(8)                :: RASA        ! ATOMIC RADIUS
      INTEGER(4)             :: LNX         ! #(PARTIAL WAVES)
      INTEGER(4),ALLOCATABLE :: LOX(:)      !(LNX) ANGULAR MOMENTA
      INTEGER(4),ALLOCATABLE :: ISCATT(:)   !(LNX) COUNTER RELATIVE TO HOMO
      INTEGER(4),ALLOCATABLE :: LOXCHI(:)   !(LNXCHI) ANGULAR MOMENTA
      INTEGER(4)             :: GID         ! GRID ID
      INTEGER(4)             :: NR          ! #(RADIAL GRID POINTS)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)  !(NR,LNX) AE PARTIAL WAVES
      REAL(8)   ,ALLOCATABLE :: UPHI(:,:)   !(NR,LNX) NODELESS PARTIAL WAVES
      INTEGER(4)             :: LN,L,I,LXX,ISVAR
      REAL(8)   ,ALLOCATABLE :: RAD(:)      ! orbital radius by linear extrapol.
      CHARACTER(3)           :: ID ! CAN BE OLD OR NEW
!     **************************************************************************
                            CALL TRACE$PUSH('lmto$dolocorb')
!
!     ==========================================================================
!     == COLLECT DATA                                                         ==
!     ==========================================================================
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(ISPECIES(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES)
      ISP=ISPECIES(IAT)
      DEALLOCATE(ISPECIES)
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETR8('AEZ',AEZ)
      CALL PERIODICTABLE$GET(AEZ,'R(ASA)',RASA)
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      ALLOCATE(ISCATT(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
      CALL SETUP$GETI4('GID',GID)
      CALL RADIAL$GETI4(GID,'NR',NR)
      ALLOCATE(AEPHI(NR,LNX))
      ALLOCATE(UPHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
      CALL SETUP$GETR8A('NDLSPHI',NR*LNX,UPHI)
!
!     ==========================================================================
!     == CONSISTENCY CHECKS                                                   ==
!     ==========================================================================
      IF(LNX.NE.LNXPHI) THEN
        CALL ERROR$MSG('INCONSISTENT VALUES OF LNX')
        CALL ERROR$STOP('lmto$DOLOCORB')
      END IF
!
!     == CONSISTENCY CHECK BETWEEN LNXCHI AND TORB =============================
      ISVAR=0
      DO LN=1,LNXPHI
        IF(TORB(LN))ISVAR=ISVAR+1
      ENDDO
      IF(ISVAR.NE.LNXCHI) THEN
        CALL ERROR$MSG('INCONSISTENT VALUES OF LNX')
        CALL ERROR$STOP('lmto$DOLOCORB')
      END IF
!
!     ==========================================================================
!     == COLLECT INFORMATION FOR LOCAL ORBITALS                               ==
!     ==========================================================================
      ALLOCATE(LOXCHI(LNXCHI))
      ISVAR=0
      DO LN=1,LNX
        IF(.NOT.TORB(LN)) CYCLE
        ISVAR=ISVAR+1
        LOXCHI(ISVAR)=LOX(LN)
      ENDDO
      LXX=MAXVAL(LOX)
      ALLOCATE(RAD(LXX+1))
PRINT*,'ISP ',IAT,AEZ
PRINT*,'ORBRAD ',ORBRAD(:,IAT)
      RAD(:)=ORBRAD(:LXX+1,IAT)
CALL LMTO_WRITEPHI('TESTAEPHI',GID,NR,LNX,AEPHI)
CALL LMTO_WRITEPHI('TESTUPHI',GID,NR,LNX,UPHI)
      CALL LMTO_CHIFROMPHI(GID,NR,LXX,RASA,RAD &
     &           ,LNX,LOX,ISCATT,AEPHI &
     &           ,LNXCHI,LOXCHI,PHICHI,CHIPHI)
CALL LMTO_WRITEPHI('TESTAECHI',GID,NR,LNXCHI,MATMUL(AEPHI,PHICHI))
CALL LMTO_WRITEPHI('TESTUCHI',GID,NR,LNXCHI,MATMUL(UPHI,PHICHI))
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_CHIFROMPHI(GID,NR,LXX,RASA,RAD &
     &                           ,LNXPHI,LOXPHI,ISCATT,PHI &
     &                           ,LNXCHI,LOXCHI,AMAT,BMAT)
!     **************************************************************************
!     **  DEFINES PROJECTOR FUNCTIONS FOR LOCAL ORBITALS IN TERMS OF          **
!     **  THE CONVENTIONAL PARTIAL WAVE PROJECTOR FUNCTIONS                   **
!     **                                                                      **
!     **     <P_CHI_I|= SUM_J BMAT(I,J) <P_PHI_J|                             **
!     **     |CHI_I>  = SUM_J |PHI_J> AMAT(J,I)                               **
!     **                                                                      **
!     **  THE LOCAL ORBITALS ARE SELECTED BY NORB WHICH SPECIFIES THE NUMBER  **
!     **  OF LOCAL ORBITALS PER ANGULAR MOMENTUM CHANNEL.                     **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID     ! GRID ID
      INTEGER(4),INTENT(IN) :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: LXX     ! X(ANGULAR MOMENTUM)
      REAL(8)   ,INTENT(IN) :: RASA    ! ATOMIC RADIUS             
      REAL(8)   ,INTENT(IN) :: RAD(LXX+1)   ! EXTENT OF HEAD FUNCTIONS
      INTEGER(4),INTENT(IN) :: LNXPHI  ! #(PARTIAL WAVES)
      INTEGER(4),INTENT(IN) :: LOXPHI(LNXPHI)! ANGULAR MOMENTUM PER PARTIAL WAVE
      INTEGER(4),INTENT(IN) :: ISCATT(LNXPHI)! RELATIVE TO HIGHEST VALENCE STATE
      REAL(8)   ,INTENT(IN) :: PHI(NR,LNXPHI)! PARTIAL WAVES
      INTEGER(4),INTENT(IN) :: LNXCHI        ! #(LOCAL ORBITALS)
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
      REAL(8)               :: R(NR)       
      REAL(8)   ,ALLOCATABLE:: G(:)        
      REAL(8)   ,PARAMETER  :: RCG=1.D-2
      REAL(8)   ,ALLOCATABLE:: AUX(:),AUX1(:)
      REAL(8)               :: SVAR1,SVAR2
      INTEGER(4)            :: IR
      INTEGER(4)            :: NX,N,LX,L,LN,LNCHI,NOFL,NOFL2,ISVAR
      INTEGER(4)            :: N1,N2,LN1,LN2,L1,L2
      CHARACTER(64)         :: FILE
      INTEGER(4)            :: NORB(LXX+1)
      REAL(8)               :: VAL1,DER1,VAL2,DER2
      REAL(8)               :: DR
      REAL(8)               :: RCUT(LNXPHI)
!     **************************************************************************
                            CALL TRACE$PUSH('lmto_CHIFROMPHI')
      CALL RADIAL$R(GID,NR,R)
      LX=MAXVAL(LOXPHI)
      NORB(:)=0
      DO LN=1,LNXCHI
        L=LOXCHI(LN)
        NORB(L+1)=NORB(L+1)+1
      ENDDO
!
      IF(LX.GT.LXX) THEN
        CALL ERROR$MSG('LX>LXX')
        CALL ERROR$STOP('lmto_CHIFROMPHI')
      END IF
!
!     ==========================================================================
!     ==  START WITH PARTIAL WAVES AS LOCAL ORBITALS                          ==
!     == |CHI_I>=SUM_I |PHI_J>A(J,I)
!     ==========================================================================
      ALLOCATE(CHI(NR,LNXPHI))     
      CHI(:,:)=0.D0
      ALLOCATE(A(LNXPHI,LNXPHI))        
      A(:,:)=0.D0
      DO LN=1,LNXPHI
        L=LOXPHI(LN)
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
        DO LN=1,LNXPHI
          IF(LOXPHI(LN).NE.L) CYCLE
          NOFL=NOFL+1
          IF(ISCATT(LN).GT.0) CYCLE ! CONSIDER ONLY HEAD FUNCTIONS
!         == SEARCH FOR NEXT PARTIAL WAVE WITH THIS L =========================
          LN1=0
          DO LN2=1,LNXPHI
            IF(LOXPHI(LN2).NE.L) CYCLE  
!alternative 1
            ln1=ln2
            IF(ISCATT(LN2).EQ.1) exit
!alternative 2
!            IF(ISCATT(LN2).EQ.1) THEN
!              LN1=LN2
!              EXIT
!            END IF
!end alternatives
          ENDDO
          IF(LN1.EQ.0) THEN
print*,'iscatt ',iscatt
print*,'loxphi ',loxphi
            CALL ERROR$MSG('CAN NOT LOCALIZE')
            CALL ERROR$MSG('NO TAIL FUNCTION AVAILABLE')
            CALL ERROR$i4val('l',l)
            CALL ERROR$i4val('ln',ln)
            CALL ERROR$STOP('lmto_CHIFROMPHI')
          END IF
!
!         == NOW CONTINUE; LN1 POINTS TO THE NEXT PHI WITH THIS L ==============
!         == IMPOSE NODE CONDITION==============================================
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN),RASA,VAL1)
          CALL RADIAL$DERIVATIVE(GID,NR,CHI(:,LN),RASA,DER1)
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN1),RASA,VAL2)
          CALL RADIAL$DERIVATIVE(GID,NR,CHI(:,LN1),RASA,DER2)
          DR=RAD(L+1)-RASA
          SVAR1=1.D0
          SVAR2=-(VAL1+DR*DER1)/(VAL2+DR*DER2)
          IF(VAL1.EQ.0.D0.AND.VAL2.EQ.0.D0) THEN
            CALL ERROR$MSG('PARTIAL WAVES ARE TRUNCATED INSIDE OF RCUT')
            CALL ERROR$MSG('THIS IS A FLAW OF THE IMPLEMENTATION')
            CALL ERROR$MSG('CHOOSE SMALLER RCUT')
            CALL ERROR$STOP('LDAPLUSU_CHIFROMPHI')
          END IF
!
          CHI(:,LN)=CHI(:,LN)*SVAR1+CHI(:,LN1)*SVAR2
          A(:,LN)=A(:,LN)*SVAR1+A(:,LN1)*SVAR2
!
!         == DETERMINE ACTUAL NODE OF THE ORBITAL ==============================
          DO IR=1,NR
            IF(R(IR).LT.RASA) CYCLE
            IF(CHI(IR,LN)*CHI(IR-1,LN).GT.0.D0) CYCLE
            RCUT(LN)=R(IR-1)-CHI(IR-1,LN)/(CHI(IR,LN)-CHI(IR-1,LN))*(R(IR)-R(IR-1))
            EXIT
          ENDDO          
          CALL RADIAL$VALUE(GID,NR,CHI(:,LN),RCUT(LN),VAL1)
          CALL RADIAL$DERIVATIVE(GID,NR,CHI(:,LN),RCUT(LN),DER1)
          RCUT(LN)=RCUT(LN)-VAL1/DER1
!
!         == ORTHOGONALIZE TO THE LOWER HEAD FUNCTIONS =========================
          DO LN1=1,LN-1
            IF(LOXPHI(LN1).NE.L) CYCLE  
            AUX(:)=CHI(:,LN)*CHI(:,LN1)*R(:)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,MIN(RCUT(LN),RCUT(LN1)),SVAR1)
            CHI(:,LN)=CHI(:,LN)-CHI(:,LN1)*SVAR1
            A(:,LN)=A(:,LN)-A(:,LN1)*SVAR1
          ENDDO
!
!         == NORMALIZE HEAD FUNCTION =============================================
          AUX(:)=CHI(:,LN)**2*R(:)**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1,RCUT(LN),SVAR1)
          SVAR1=1.D0/SQRT(SVAR1)
          CHI(:,LN)=CHI(:,LN)*SVAR1
          A(:,LN)=A(:,LN)*SVAR1
        ENDDO
      END DO          
      DEALLOCATE(AUX)
      DEALLOCATE(AUX1)
!
!     ==========================================================================
!     == NOW INVERT THE MATRIX A: B:=A^{-1}                                   ==
!     == |CHI_I>=SUM_J |PHI_J>A(J,I)                                          ==
!     == 1=|PHI><P|=|PHI>A B <P|=|CHI> ( B<P| )                               ==
!     == <P_CHI_I|=SUM_J B(I,J)<P_PHI_J|                                      ==
!     ==========================================================================
      ALLOCATE(B(LNXPHI,LNXPHI))
      B(:,:)=0.D0
      LX=MAXVAL(LOXPHI)
      DO L=0,LX
        IF(NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
!
        NX=0
        DO LN=1,LNXPHI
          IF(LOXPHI(LN).EQ.L) NX=NX+1
        ENDDO
        IF(NX.EQ.0) CYCLE   

        ALLOCATE(MAT(NX,NX))
        ALLOCATE(MATINV(NX,NX))
!        
!       == MAP A INTO MATRIX MAT ===============================================
        N1=0
        DO LN1=1,LNXPHI
          L1=LOXPHI(LN1)
          IF(L1.NE.L) CYCLE
          N1=N1+1
          N2=0
          DO LN2=1,LNXPHI
            L2=LOXPHI(LN2)
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
        DO LN1=1,LNXPHI
          L1=LOXPHI(LN1)
          IF(L1.NE.L) CYCLE
          N1=N1+1
          N2=0
          DO LN2=1,LNXPHI
            L2=LOXPHI(LN2)
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
        DO LN=1,LNXPHI
          IF(L.NE.LOXPHI(LN)) CYCLE
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
                            CALL TRACE$POP()
      RETURN
      END
