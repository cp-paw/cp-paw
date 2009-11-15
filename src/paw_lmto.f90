MODULE LMTO_MODULE
REAL(8)   ,PARAMETER  :: K2=-1.D0 ! 0.5*k2 is the kinetic energy
REAL(8)   ,PARAMETER  :: RC=5.D0  ! CUTOF RADIUS FOR NEIGHBORLIST
TYPE POTPAR_TYPE
  REAL(8)          :: RAD
  REAL(8),POINTER  :: QBAR(:)
  REAL(8),POINTER  :: CBAR(:)
  REAL(8),POINTER  :: SQDELTABAR(:)
  REAL(8),POINTER  :: ENU(:)
  REAL(8),POINTER  :: ESCATT(:)
  REAL(8),POINTER  :: A(:)
  REAL(8),POINTER  :: OBAR(:)     !<PHIBARDOT|PHI>/<PHI|PHI>
  REAL(8),POINTER  :: PBAR(:)     !<PHIBARDOT|PHIBARDOT>/<PHI|PHI>
  REAL(8),POINTER  :: PHIPHI(:)   !<PHI|PHI>
END TYPE POTPAR_TYPE
TYPE PERIODICMAT_TYPE
  INTEGER(4)      :: IAT1
  INTEGER(4)      :: IAT2
  INTEGER(4)      :: IT(3)
  INTEGER(4)      :: N1
  INTEGER(4)      :: N2
  REAL(8),POINTER :: MAT(:,:)
END TYPE PERIODICMAT_TYPE
LOGICAL(4)              :: TINISTRUC=.FALSE.
INTEGER(4)              :: NSP
INTEGER(4)              :: LXX
INTEGER(4),ALLOCATABLE  :: LX(:)               !(NSP)
REAL(8)   ,ALLOCATABLE  :: RAD(:)              !(NSP)
INTEGER(4),ALLOCATABLE  :: ISPECIES(:)         !(NAT)
REAL(8)   ,ALLOCATABLE  :: ORBRAD(:,:) !(LXX+1,NAT) NODE-POSITION OF THE ORBITAL
TYPE(POTPAR_TYPE)     ,ALLOCATABLE :: POTPAR(:)
TYPE(PERIODICMAT_TYPE),ALLOCATABLE :: SBAR(:) !(NNB)
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
                             CALL TRACE$PUSH('LMTO$REPORTSBAR')
      NNB=SIZE(SBAR)
      DO NN=1,NNB
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
      USE LMTO_MODULE, ONLY : POTPAR
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4)            :: NSP
      INTEGER(4)            :: LNX
      CHARACTER(64)         :: TITLE
      INTEGER(4)            :: ISP,LN
      INTEGER(4),ALLOCATABLE:: LOX(:)
!     **************************************************************************
                             CALL TRACE$PUSH('LMTO$REPORTPOTPAR')
      NSP=SIZE(POTPAR)
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        LNX=SIZE(POTPAR(ISP)%QBAR)
        ALLOCATE(LOX(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$ISELECT(0)
        DO LN=1,LNX
          WRITE(TITLE,FMT='(I3," LN=",I3," L=",I2)')ISP,LN,LOX(LN)
          TITLE='POTENTIAL PARAMETERS FOR ATOM TYPE '//TRIM(ADJUSTL(TITLE))
          CALL REPORT$TITLE(NFIL,TITLE)
          CALL REPORT$R8VAL(NFIL,'RAD',POTPAR(ISP)%RAD,'ABOHR')
          CALL REPORT$R8VAL(NFIL,'QBAR',POTPAR(ISP)%QBAR(LN),'')
          CALL REPORT$R8VAL(NFIL,'ENU',POTPAR(ISP)%ENU(LN),'H')
          CALL REPORT$R8VAL(NFIL,'ESCATT',POTPAR(ISP)%ESCATT(LN),'H')
          CALL REPORT$R8VAL(NFIL,'CBAR',POTPAR(ISP)%CBAR(LN),'H')
          CALL REPORT$R8VAL(NFIL,'DELTA',POTPAR(ISP)%SQDELTABAR(LN)**2,'H')
          CALL REPORT$R8VAL(NFIL,'OBAR',POTPAR(ISP)%OBAR(LN),'1/H')
          CALL REPORT$R8VAL(NFIL,'PBAR',POTPAR(ISP)%PBAR(LN),'1/H**2')
          CALL REPORT$R8VAL(NFIL,'PHIPHI',POTPAR(ISP)%PHIPHI(LN),'')
          CALL REPORT$R8VAL(NFIL,'A',POTPAR(ISP)%A(LN),'')
        ENDDO
        WRITE(NFIL,FMT='(82("="))')
        DEALLOCATE(LOX)
      ENDDO
                                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$MAKEPOTPAR()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : K2,POTPAR
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4)             :: NSP
      INTEGER(4),ALLOCATABLE :: LX(:)
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4),ALLOCATABLE :: ISCATT(:)
      INTEGER(4)             :: LNX
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      REAL(8)   ,ALLOCATABLE :: R(:)
      INTEGER(4)             :: LNHOMO
      REAL(8)   ,ALLOCATABLE :: EOFLN(:)
      REAL(8)   ,ALLOCATABLE :: ESCATT(:)
      REAL(8)   ,ALLOCATABLE :: NLPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHI(:,:)
      REAL(8)   ,ALLOCATABLE :: NLPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: AEPHIDOT(:,:)
      REAL(8)   ,ALLOCATABLE :: AUX(:),AUX1(:)
      REAL(8)                :: AEZ
      REAL(8)                :: RAD
      REAL(8)                :: PHIVAL,PHIDER
      REAL(8)                :: PHIDOTVAL,PHIDOTDER
      REAL(8)                :: KVAL,KDER
      REAL(8)                :: JVAL,JDER
      REAL(8)                :: WJPHI,WJPHIDOT,WKPHI,WKPHIDOT,WJBARPHI
      REAL(8)                :: QBAR
      REAL(8)                :: OBAR
      REAL(8)                :: PBAR
      REAL(8)                :: CBAR
      REAL(8)                :: SQDELTABAR
      REAL(8)                :: A
      REAL(8)                :: PHIINT
      REAL(8)                :: ENU
      REAL(8)                :: EGAMMA
      INTEGER(4)             :: ISP,LN,L,ln1,ir
real(8):: x1,x2,x3
!     **************************************************************************
      CALL SETUP$NSPECIES(NSP)
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
      DEALLOCATE(LX)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ALLOCATE(POTPAR(NSP))
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(ISCATT(LNX))
        ALLOCATE(LOX(LNX))
        ALLOCATE(EOFLN(LNX))
        ALLOCATE(ESCATT(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
        CALL SETUP$GETR8A('EOFLN',LNX,EOFLN)
        CALL SETUP$GETR8A('ESCATT',LNX,ESCATT)
        CALL SETUP$GETI4('GID',GID)
        CALL SETUP$GETI4('NR',NR)
        ALLOCATE(R(NR))
        CALL RADIAL$R(GID,NR,R)
        ALLOCATE(NLPHI(NR,LNX))
        ALLOCATE(AEPHI(NR,LNX))
        ALLOCATE(NLPHIDOT(NR,LNX))
        ALLOCATE(AEPHIDOT(NR,LNX))
        CALL SETUP$GETR8A('NLPHI',NR*LNX,NLPHI)
        CALL SETUP$GETR8A('AEPHI',NR*LNX,AEPHI)
        CALL SETUP$GETR8A('NLPHIDOT',NR*LNX,NLPHIDOT)
        CALL SETUP$GETR8A('AEPHIDOT',NR*LNX,AEPHIDOT)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(NINT(AEZ),'R(ASA)',RAD) 
        POTPAR(ISP)%RAD=RAD
!
        ALLOCATE(POTPAR(ISP)%QBAR(LNX))
        ALLOCATE(POTPAR(ISP)%CBAR(LNX))
        ALLOCATE(POTPAR(ISP)%SQDELTABAR(LNX))
        ALLOCATE(POTPAR(ISP)%A(LNX))
        ALLOCATE(POTPAR(ISP)%OBAR(LNX))
        ALLOCATE(POTPAR(ISP)%PBAR(LNX))
        ALLOCATE(POTPAR(ISP)%PHIPHI(LNX))
        ALLOCATE(POTPAR(ISP)%ENU(LNX))
        ALLOCATE(POTPAR(ISP)%ESCATT(LNX))
!
        ALLOCATE(AUX(NR))
        ALLOCATE(AUX1(NR))
        DO L=0,MAXVAL(LOX)
!         == select phibardot from valence channel =============================
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
            if(iscatt(ln).gt.0) cycle
            ln1=ln
          enddo
          DO LN=1,LNX
            IF(LOX(LN).NE.L) CYCLE
!ln1=ln ! old version re-established
!           ====================================================================
!           == value and derivative of partial waves and envelope functions   ==
!           ====================================================================
            CALL RADIAL$VALUE(GID,NR,NLPHIDOT(:,LN1),RAD,PHIDOTVAL)
            CALL RADIAL$DERIVATIVE(GID,NR,NLPHIDOT(:,LN1),RAD,PHIDOTDER)
            CALL RADIAL$VALUE(GID,NR,NLPHI(:,LN),RAD,PHIVAL)
            CALL RADIAL$DERIVATIVE(GID,NR,NLPHI(:,LN),RAD,PHIDER)
            CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
            CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
PRINT*,'ISP,L,LN',ISP,L,LN,'=========================================='
PRINT*,'PHIDOT ', PHIDOTVAL,PHIDOTDER,rad-phidotval/phidotder
PRINT*,'PHI    ', PHIVAL,PHIDER
PRINT*,'K      ', KVAL,KDER
PRINT*,'J      ', JVAL,JDER
!
!           ====================================================================
!           == CALCULATE POTENTIAL PARAMETERS                                 ==
!           ====================================================================
            WJPHI=JVAL*PHIDER-JDER*PHIVAL
            WJPHIDOT=JVAL*PHIDOTDER-JDER*PHIDOTVAL
            WKPHI=KVAL*PHIDER-KDER*PHIVAL
            WKPHIDOT=KVAL*PHIDOTDER-KDER*PHIDOTVAL
!           ==  screening charge ===============================================
            QBAR=WJPHIDOT/WKPHIDOT
print*,'qbar ',qbar
            WJBARPHI=WJPHI-WKPHI*QBAR
            ENU=EOFLN(LN)
            EGAMMA=ESCATT(LN1)
!           == band center =====================================================
            CBAR=ENU-WKPHI/WKPHIDOT
!           ==  phiint=<phi|phi> ===============================================
            AUX(:)=R(:)**2*AEPHI(:,LN)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RAD,PHIINT)
!           ==  obar=<phi|phibardot>/<phi|phi>  ================================
            AUX(:)=R(:)**2*AEPHI(:,LN)*AEPHIDOT(:,LN1)
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RAD,OBAR)
            OBAR=OBAR/PHIINT 
!           ==  pbar=<phibardot|phibardot>/<phi|phi> ===========================
            AUX(:)=R(:)**2*AEPHIDOT(:,LN1)**2
            CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
            CALL RADIAL$VALUE(GID,NR,AUX1,RAD,PBAR)
            PBAR=PBAR/PHIINT 
!           == approximate normalization factor for phi ========================
            A=1.D0/SQRT((1.D0-(ENU-EGAMMA)*OBAR)*PHIINT)
!           == band width ======================================================
            SQDELTABAR=-RAD**2*A*WJBARPHI/SQRT(2.D0)
!
!           ====================================================================
!           == NOW MAP ONTO POTPAR                                            ==
!           ====================================================================
            POTPAR(ISP)%QBAR(LN)=QBAR
            POTPAR(ISP)%CBAR(LN)=CBAR
            POTPAR(ISP)%SQDELTABAR(LN)=SQDELTABAR
            POTPAR(ISP)%A(LN)=A
            POTPAR(ISP)%PHIPHI(LN)=PHIINT
            POTPAR(ISP)%OBAR(LN)=OBAR
            POTPAR(ISP)%PBAR(LN)=PBAR
            POTPAR(ISP)%ENU(LN)=ENU
            POTPAR(ISP)%ESCATT(LN)=EGAMMA
          ENDDO
        ENDDO       
        DEALLOCATE(EOFLN)
        DEALLOCATE(ESCATT)
        DEALLOCATE(NLPHI)
        DEALLOCATE(AEPHI)
        DEALLOCATE(NLPHIDOT)
        DEALLOCATE(AEPHIDOT)
        DEALLOCATE(LOX)
        DEALLOCATE(ISCATT)
        DEALLOCATE(R)
        DEALLOCATE(AUX)
        DEALLOCATE(AUX1)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$MAKESTRUCTURECONSTANTS()
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : K2,RC,SBAR,TINISTRUC,POTPAR
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
      REAL(8)   ,ALLOCATABLE :: QBAR1SP(:) !
      INTEGER(4),ALLOCATABLE :: ISCATT(:)
      REAL(8)   ,ALLOCATABLE :: C(:,:) !
      REAL(8)                :: VAL,DER
      INTEGER(4)             :: IAT,IAT1,IAT2,ISP,ISP1,ISP2,LMX1,LMX2 
      INTEGER(4)             :: L,NN,NN1,NN2,NN0,I,IM,LN
      INTEGER(4)             :: LXX        !X(ANGULAR MOMENTUM FOR ALL ATOMS)
      INTEGER(4)             :: LNX
      INTEGER(4),ALLOCATABLE :: LX(:)      !X(ANGULAR MOMENTUM FOR EACH ATOM)
      INTEGER(4),ALLOCATABLE :: LOX(:)     !ANGULAR MOMENTA
      INTEGER(4)             :: NSP       !#(ATOM TYPES)
      INTEGER(4),ALLOCATABLE :: ISPECIES(:)
      LOGICAL(4)             :: TCHK
!     **************************************************************************
      CALL SETUP$GETL4('INTERNALSETUPS',TCHK)
      IF(.NOT.TCHK) RETURN
                              CALL TRACE$PUSH('LMTO$MAKESTRUCTURECONSTANTS')
!
      IF(TINISTRUC) RETURN
      TINISTRUC=.TRUE.
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
      CALL LMTO$MAKEPOTPAR()
      CALL SETUP$NSPECIES(NSP)
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
      ALLOCATE(QBAR1(LXX+1,NSP))
      QBAR1(:,:)=0.D0
      DO ISP=1,NSP
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETI4('LNX',LNX)
        ALLOCATE(LOX(LNX))
        ALLOCATE(ISCATT(LNX))
        CALL SETUP$GETI4A('LOX',LNX,LOX)
        CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
        DO LN=1,LNX
          IF(ISCATT(LN).NE.0) CYCLE
          QBAR1(LOX(LN)+1,ISP)=POTPAR(ISP)%QBAR(LN)
        ENDDO
        DEALLOCATE(ISCATT)
        DEALLOCATE(LOX)
        CALL SETUP$ISELECT(0)
      ENDDO
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
            ENDDO
          ENDDO
        ENDDO
!
!       ========================================================================
!       == DETERMINE STRUCTURE CONSTANTS                                      ==
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
          ALLOCATE(SBAR(NN)%MAT(LMX2,LMX1))
          SBAR(NN)%MAT(:,:)=SBAR1(I+1:I+LMX2,:)
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
      CALL LMTO$REPORTPOTBAR(6)
      CALL LMTO$REPORTSBAR(6)
!     CALL LMTO_plotlocorb(1)
!
                              CALL TRACE$POP()
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
            IF(ISCATT(LN).NE.0) CYCLE
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
      INTEGER(4)             :: NAT         ! #(ATOMS)
      INTEGER(4),ALLOCATABLE :: ISPECIES(:) !(NAT) ATOM TYPE FOR EACH ATOM
      INTEGER(4)             :: ISP         ! ATOM TYPE INDEX
      INTEGER(4)             :: LNX         ! #(PARTIAL WAVES)
      INTEGER(4)             :: NNB         ! 
      INTEGER(4)             :: LOX(LNXPHI) ! ANGULAR MOMENTA
      INTEGER(4),ALLOCATABLE :: LOXCHI(:)    ! ANGULAR MOMENTA
      INTEGER(4)             :: ISCATT(LNXPHI)    ! COUNTER RELATIVE TO HOMO
      REAL(8)                :: AEPHI(NR,LNXPHI)  ! AE PARTIAL WAVES
      REAL(8)                :: PSPHI(NR,LNXPHI)  ! PSEUDO PARTIAL WAVES
      REAL(8)                :: AEPHIDOT(NR,LNXPHI) ! AE PARTIAL WAVES
      REAL(8)                :: PSPHIDOT(NR,LNXPHI) ! PARTIAL WAVES
      REAL(8)                :: nlphi(NR,LNXPHI) ! PARTIAL WAVES
      REAL(8)                :: nlPHIDOT(NR,LNXPHI) ! PARTIAL WAVES
      REAL(8)                :: PRO(NR,LNXPHI)    ! PROJECTOR FUNCTIONS
      REAL(8)   ,ALLOCATABLE :: AECHI(:,:)  ! ALL-ELECTRON LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: PSCHI(:,:)  ! PSEUDO LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: nlCHI(:,:)  ! nodeless LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: SBARONSITE(:,:)   !(LMX,LMX) 
      REAL(8)   ,ALLOCATABLE :: AMAT(:,:)   !(LNXCHI1,LNXPHI) 
      REAL(8)   ,ALLOCATABLE :: BMAT(:,:)   !(LNXCHI1,LNXCHI1) 
      REAL(8)   ,ALLOCATABLE :: XMAT(:,:)   !(LNXPHI,LNXCHI) 
      REAL(8)                :: RCOV              ! COVALENT RADIUS
      REAL(8)                :: AEZ               ! ATOMIC NUMBER
      REAL(8)                :: AUX(NR)
      REAL(8)                :: R(NR)
      REAL(8)                :: SVAR,SVAR1,VAL,der,vald,derd
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: LMX,LNXCHI1
      INTEGER(4)             :: LN,L,I,J,LXX,ISVAR,IIB,N1,LM,LNCHI1,LNCHI,IR
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
      NNB=SIZE(SBAR)
      TCHK=.FALSE.
      DO IIB=1,NNB
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
      ALLOCATE(nlCHI(NR,LNXCHI1))
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
        nlCHI(:,LNCHI)=nlPHI(:,LN)          
        LM=L**2+1
        SVAR=POTPAR(ISP)%CBAR(LN)-POTPAR(ISP)%ENU(LN) &
     &        +SBARONSITE(LM,LM)*POTPAR(ISP)%SQDELTABAR(LN)**2
PRINT*,'LN,LNCHI,L, SVAR ',LN,LNCHI,L,SVAR,AEPHIDOT(IRCOV,LN),AECHI(IRCOV,LNCHI)
PRINT*,'C,ENU   ',POTPAR(ISP)%CBAR(LN),POTPAR(ISP)%ENU(LN)
!IF(AEPHIDOT(IRCOV,LN)*SVAR/AECHI(IRCOV,LNCHI).GT.0.D0) CYCLE
!        IF(SVAR.LT.0.D0) CYCLE
        AECHI(:,LNCHI)=AECHI(:,LNCHI)+AEPHIDOT(:,LN)*SVAR
        PSCHI(:,LNCHI)=PSCHI(:,LNCHI)+PSPHIDOT(:,LN)*SVAR
        nlCHI(:,LNCHI)=nlCHI(:,LNCHI)+nlPHIDOT(:,LN)*SVAR
      ENDDO
print*,'lnx,lox ',lnx,lox
print*,'lnxchi1,loxchi ',lnxchi1,loxchi
!
!     == attach exponential tail at the covalent radius ========================
      DO LNCHI=1,LNXCHI1     
        l=loxchi(lnchi)
        CALL RADIAL$VALUE(GID,NR,NLCHI(:,LNCHI),RCOV,VAL)
        CALL RADIAL$DERIVATIVE(GID,NR,NLCHI(:,LNCHI),RCOV,DER)
        CALL RADIAL$VALUE(GID,NR,AECHI(:,LNCHI),RCOV,VALD)
        CALL RADIAL$DERIVATIVE(GID,NR,AECHI(:,LNCHI),RCOV,DERD)
        VALD=VALD-VAL
        DERD=DERD-DER
!
        SVAR=DER/VAL
        if(abs(vald/val).gt.1.d-9) then ! avoid divide-by-zero
          SVAR1=DERD/VALD
        else
          svar1=-1.d0
        end if        
! 
        IF(SVAR.GT.0.D0.or.svar1.gt.0.d0) THEN
          CALL SETUP_WRITEPHI(-'FAILEDAEORBITAL.DAT',GID,NR,LNXCHI1,AECHI)
          CALL SETUP_WRITEPHI(-'FAILEDPSORBITAL.DAT',GID,NR,LNXCHI1,PSCHI)
          CALL SETUP_WRITEPHI(-'FAILEDNLORBITAL.DAT',GID,NR,LNXCHI1,NLCHI)
          CALL SETUP_WRITEPHI(-'FAILEDdiffORBITAL.DAT',GID,NR,LNXCHI1,aechi-NLCHI)
          CALL ERROR$MSG('ORBITAL DOES NOT DECAY')
          CALL ERROR$I4VAL('IAT',IAT)
          CALL ERROR$I4VAL('ISP',ISP)
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('LNCHI',LNCHI)
          CALL ERROR$I4VAL('LOXCHI',LOXCHI(lnchi))
          CALL ERROR$R8VAL('COVALENT RADIUS',RCOV)
          CALL ERROR$R8VAL('LOGARITMIC DERIVATIVE ',SVAR)
          CALL ERROR$R8VAL('LOGARITMIC DERIVATIVE OF DIFFERENCE',SVAR1)
          CALL ERROR$STOP('LMTO$DOLOCORB')
        END IF
        DO IR=IRCOV,NR
          nlCHI(IR,LNCHI)=VAL*EXP(SVAR*(R(IR)-RCOV))
          AECHI(IR,LNCHI)=VAL*EXP(SVAR*(R(IR)-RCOV)) &
                         +vald*EXP(SVAR1*(R(IR)-RCOV))
          psCHI(IR,LNCHI)=VAL*EXP(SVAR*(R(IR)-RCOV)) &
                         +vald*EXP(SVAR1*(R(IR)-RCOV))
        ENDDO
      enddo
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
      SUBROUTINE LMTO$DOLOCORB_OLD2(IAT,ISP,GID,NR,LNXCHI,LNXPHI,TORB,CHIPHI,CHI)
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
      INTEGER(4)             :: NAT         ! #(ATOMS)
      INTEGER(4),ALLOCATABLE :: ISPECIES(:) !(NAT) ATOM TYPE FOR EACH ATOM
      INTEGER(4)             :: ISP         ! ATOM TYPE INDEX
      INTEGER(4)             :: LNX         ! #(PARTIAL WAVES)
      INTEGER(4)             :: NNB         ! 
      INTEGER(4)             :: LOX(LNXPHI) ! ANGULAR MOMENTA
      INTEGER(4),ALLOCATABLE :: LOXCHI(:)    ! ANGULAR MOMENTA
      INTEGER(4)             :: ISCATT(LNXPHI)    ! COUNTER RELATIVE TO HOMO
      REAL(8)                :: AEPHI(NR,LNXPHI)  ! AE PARTIAL WAVES
      REAL(8)                :: PSPHI(NR,LNXPHI)  ! PSEUDO PARTIAL WAVES
      REAL(8)                :: AEPHIDOT(NR,LNXPHI) ! AE PARTIAL WAVES
      REAL(8)                :: PSPHIDOT(NR,LNXPHI) ! PARTIAL WAVES
      REAL(8)                :: PRO(NR,LNXPHI)    ! PROJECTOR FUNCTIONS
      REAL(8)   ,ALLOCATABLE :: AECHI(:,:)  ! ALL-ELECTRON LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: PSCHI(:,:)  ! PSEUDO LOCAL ORBITALS
      REAL(8)   ,ALLOCATABLE :: SBARONSITE(:,:)   !(LMX,LMX) 
      REAL(8)   ,ALLOCATABLE :: AMAT(:,:)   !(LNXCHI1,LNXPHI) 
      REAL(8)   ,ALLOCATABLE :: BMAT(:,:)   !(LNXCHI1,LNXCHI1) 
      REAL(8)   ,ALLOCATABLE :: XMAT(:,:)   !(LNXPHI,LNXCHI) 
      REAL(8)                :: RCOV              ! COVALENT RADIUS
      REAL(8)                :: AEZ               ! ATOMIC NUMBER
      REAL(8)                :: AUX(NR)
      REAL(8)                :: R(NR)
      REAL(8)                :: SVAR,SVAR1,VAL
      LOGICAL(4)             :: TCHK
      INTEGER(4)             :: LMX,LNXCHI1
      INTEGER(4)             :: LN,L,I,J,LXX,ISVAR,IIB,N1,LM,LNCHI1,LNCHI,IR
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
      CALL SETUP$GETR8A('AEPHIDOT',NR*LNX,AEPHIDOT)
      CALL SETUP$GETR8A('PSPHIDOT',NR*LNX,PSPHIDOT)
      CALL SETUP$GETR8A('PRO',NR*LNX,PRO)
!
!     ==========================================================================
!     == FIND ONSITE STRUCTURE CONSTANTS                                      ==
!     ==========================================================================
      NNB=SIZE(SBAR)
      TCHK=.FALSE.
      DO IIB=1,NNB
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
        LM=L**2+1
        SVAR=POTPAR(ISP)%CBAR(LN)-POTPAR(ISP)%ENU(LN) &
     &        +SBARONSITE(LM,LM)*POTPAR(ISP)%SQDELTABAR(LN)**2
PRINT*,'LN,LNCHI,L, SVAR ',LN,LNCHI,L,SVAR,AEPHIDOT(IRCOV,LN),AECHI(IRCOV,LNCHI)
PRINT*,'C,ENU   ',POTPAR(ISP)%CBAR(LN),POTPAR(ISP)%ENU(LN)
!IF(AEPHIDOT(IRCOV,LN)*SVAR/AECHI(IRCOV,LNCHI).GT.0.D0) CYCLE
        IF(SVAR.LT.0.D0) CYCLE
        AECHI(:,LNCHI)=AECHI(:,LNCHI)+AEPHIDOT(:,LN)*SVAR
        PSCHI(:,LNCHI)=PSCHI(:,LNCHI)+PSPHIDOT(:,LN)*SVAR
      ENDDO
!
!     == CUT OFF THE TAIL OF THE LOCAL ORBITALS ================================
      DO LNCHI=1,LNXCHI1     
        SVAR=AECHI(IRCOV,LNCHI)/AECHI(IRCOV-1,LNCHI)
!       SVAR>1 : CHI INCREASES IN ABSOLUTE VALUE AT THE COVALENT RADIUS
!       0<SVAR<1: CHI DECREASS IN ABSOLUTE VALUE
!       SVAR<0 : ZERO BETWEEN R(IRCOV-1) AND R(IRCOV)
        IF(SVAR.GT.1.D0) THEN
          DO IR=IRCOV-1,2,-1
            SVAR1=AECHI(IR,LNCHI)/AECHI(IR-1,LNCHI)
            IF(SVAR1.LE.0.D0) THEN
PRINT*,'ORBITAL CUTOFF/RCOV',LNCHI,R(IR)/RCOV,' R=',R(IR)
               AECHI(IR:,LNCHI)=0.D0
               PSCHI(IR:,LNCHI)=0.D0
               EXIT
            END IF
            IF(SVAR1.LT.1.D0) THEN
              CALL LMTO_WRITEPHI('FAILEDLOCALORBITAL.DAT',GID,NR,LNXCHI1,AECHI)
              CALL LMTO_WRITEPHI('FAILEDPHI.DAT',GID,NR,LNXCHI1,AEPHI)
              CALL LMTO_WRITEPHI('FAILEDPHIDOT.DAT',GID,NR,LNXCHI1,AEPHIDOT)
              CALL ERROR$MSG('LOCAL ORBITAL CONSTRUCTION FAILED')
              CALL ERROR$MSG('ABSOLUTE VALUE OF LOCAL ORBITAL ...')
              CALL ERROR$MSG('... HAS A MINIMUM INSIDE RCOV')
              CALL ERROR$MSG('LOCAL ORBITALS ARE PRINTED INTO FILE')
              CALL ERROR$CHVAL('FILE','FAILEDLOCALORBITAL.DAT')
              CALL ERROR$I4VAL('IAT',IAT)
              CALL ERROR$I4VAL('ISP',ISP)
              CALL ERROR$I4VAL('LNCHI ',LNCHI)
              CALL ERROR$I4VAL('L ',LOXCHI(LNCHI))
              CALL ERROR$R8VAL('R(IRCOV)',R(IRCOV))
              CALL ERROR$R8VAL('R',R(IR))
              CALL ERROR$STOP('LMTO$DOLOCORB')
            END IF
          ENDDO
        ELSE  
          DO IR=IRCOV,NR
            SVAR1=AECHI(IR,LNCHI)/AECHI(IR-1,LNCHI)
            IF(SVAR1.LE.0.D0) THEN
PRINT*,'ORBITAL CUTOFF/RCOV',LNCHI,R(IR)/RCOV,' R=',R(IR)
              AECHI(IR:,LNCHI)=0.D0
              PSCHI(IR:,LNCHI)=0.D0
              EXIT
            END IF
!           == CHECK IF A CHI HAS A NODE AT ALL... ===============================
            IF(SVAR1.GT.1.D0) THEN
              CALL LMTO_WRITEPHI('FAILEDLOCALORBITAL.DAT',GID,NR,LNXCHI1,AECHI)
              CALL LMTO_WRITEPHI('FAILEDPHI.DAT',GID,NR,LNXCHI1,AEPHI)
              CALL LMTO_WRITEPHI('FAILEDPHIDOT.DAT',GID,NR,LNXCHI1,AEPHIDOT)
              CALL ERROR$MSG('LOCAL ORBITAL CONSTRUCTION FAILED')
              CALL ERROR$MSG('ORBITAL DOES NOT HAVE A NODE BEFORE INCREASING AGAIN')
              CALL ERROR$CHVAL('LOCAL ORBITALS ARE PRINTED TO FILE','FAILEDLOCALORBITAL.DAT')
              CALL ERROR$I4VAL('IAT',IAT)
              CALL ERROR$I4VAL('ISP',ISP)
              CALL ERROR$I4VAL('LNCHI ',LNCHI)
              CALL ERROR$I4VAL('L ',LOXCHI(LNCHI))
              CALL ERROR$R8VAL('R(IR)',R(IR))
              CALL ERROR$R8VAL('R(IRCOV)',R(IRCOV))
              CALL ERROR$STOP('LMTO$DOLOCORB')
            END IF
          ENDDO
         END IF
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
      SUBROUTINE LMTO$DOLOCORB_OLD(IAT,LNXCHI,LNXPHI,TORB,PHICHI,CHIPHI)
!     **************************************************************************
!     **  CONSTRUCTS ONSITE-MAPPING FROM PARTIAL WAVES TO LOCAL ORBITALS      **
!     **                                                                      **
!     **      |CHI_I> =  SUM_J |PHI_J> * PHICHI_JI                            **
!     **                                                                      **
!     **  IT USES ORBRAD, WHICH IS THE RADIUS OBTAINED BY LINEAR EXTRAPOLATION**
!     **  OF THE LMTO TO ZERO. ORBRAD ENCODES THE LOGARITHMIC DERIVATIVE      **
!     **                                                                      **
!     **************************************************************************
      USE LMTO_MODULE, ONLY : ORBRAD
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT     ! ATOM INDEX
      INTEGER(4),INTENT(IN)  :: LNXCHI  ! #(LOCAL ORBITALS W/O M-MULTIPLICITY)
      INTEGER(4),INTENT(IN)  :: LNXPHI  ! #(PARTIAL WAVES W/O M-MULTIPLICITY)
      LOGICAL(4),INTENT(IN)  :: TORB(LNXPHI)          ! SELECTS LOCAL ORBITALS
      REAL(8)   ,INTENT(OUT) :: PHICHI(LNXPHI,LNXCHI) !|CHI_I>=\SUM_J|PHI_J>*PHICHI_{JI}
      REAL(8)   ,INTENT(OUT) :: CHIPHI(LNXCHI,LNXPHI)
      INTEGER(4)             :: NAT         ! #(ATOMS)
      INTEGER(4),ALLOCATABLE :: ISPECIES(:) !(NAT) ATOM TYPE FOR EACH ATOM
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
      REAL(8)   ,ALLOCATABLE :: RAD(:)      ! ORBITAL RADIUS BY LINEAR EXTRAPOL.
      CHARACTER(3)           :: ID ! CAN BE OLD OR NEW
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO$DOLOCORB')
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
      CALL SETUP$GETR8A('NLPHI',NR*LNX,UPHI)
!
!     ==========================================================================
!     == CONSISTENCY CHECKS                                                   ==
!     ==========================================================================
      IF(LNX.NE.LNXPHI) THEN
        CALL ERROR$MSG('INCONSISTENT VALUES OF LNX')
        CALL ERROR$STOP('LMTO$DOLOCORB')
      END IF
!
!     == CONSISTENCY CHECK BETWEEN LNXCHI AND TORB =============================
      ISVAR=0
      DO LN=1,LNXPHI
        IF(TORB(LN))ISVAR=ISVAR+1
      ENDDO
      IF(ISVAR.NE.LNXCHI) THEN
        CALL ERROR$MSG('INCONSISTENT VALUES OF LNX')
        CALL ERROR$STOP('LMTO$DOLOCORB')
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
                            CALL TRACE$PUSH('LMTO_CHIFROMPHI')
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
        CALL ERROR$STOP('LMTO_CHIFROMPHI')
      END IF
!
!     ==========================================================================
!     ==  START WITH PARTIAL WAVES AS LOCAL ORBITALS                          ==
!     == |CHI_I>=SUM_I |PHI_J>A(J,I)                                          ==
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
!ALTERNATIVE 1
            LN1=LN2
            IF(ISCATT(LN2).EQ.1) EXIT
!ALTERNATIVE 2
!            IF(ISCATT(LN2).EQ.1) THEN
!              LN1=LN2
!              EXIT
!            END IF
!END ALTERNATIVES
          ENDDO
          IF(LN1.EQ.0) THEN
PRINT*,'ISCATT ',ISCATT
PRINT*,'LOXPHI ',LOXPHI
            CALL ERROR$MSG('CAN NOT LOCALIZE')
            CALL ERROR$MSG('NO TAIL FUNCTION AVAILABLE')
            CALL ERROR$I4VAL('L',L)
            CALL ERROR$I4VAL('LN',LN)
            CALL ERROR$STOP('LMTO_CHIFROMPHI')
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
!         == NORMALIZE HEAD FUNCTION ===========================================
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
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO_PLOTLOCORB(IAT0)
!     **************************************************************************
!     **  WRITES THE ENVELOPE FUNCTIONS CENTERED AT ATOM IAT0 TO FILE         **
!     **                                                                      **
!     **  SET N1,N2,N3 EQUAL TO THE NUMBER OF GRIDPOINTS IN EACH DIRECTION    **
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE PERIODICTABLE_MODULE
      USE STRINGS_MODULE
      USE LMTO_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT0
      INTEGER(4),PARAMETER  :: N1=50,N2=50,N3=50 !GRID (1D?)
      REAL(8)               :: ORIGIN(3)
      REAL(8)               :: TVEC(3,3)    !BOX
      REAL(8)               :: TLITTLE(3,3)
      REAL(8)               :: RBAS(3,3)
      INTEGER(4)            :: NAT
      INTEGER(4)            :: LM1
      INTEGER(4)            :: LM1X         !#Ylm
      REAL(8)   ,ALLOCATABLE:: R0(:,:)      !(3,NAT)
      INTEGER(4),ALLOCATABLE:: ISPECIES1(:)  !(NAT)
      REAL(8)   ,ALLOCATABLE:: RCOV(:)      !(NAT)
      REAL(8)               :: XI(3)
      REAL(8)               :: DR(3)
      INTEGER(4)            :: NNB
      INTEGER(4)            :: NN,IAT,IAT2,ISP,I1,I2,I3,LN
      INTEGER(4)            :: L
      INTEGER(4)            :: LNX
      REAL(8)               :: AEZ
      INTEGER(4)            :: LM2X
      INTEGER(4),PARAMETER  :: LMXX=36
      REAL(8)               :: K0(LMXX)
      REAL(8)               :: J0(LMXX)
      REAL(8)               :: JBAR(LMXX)
      REAL(8)   ,ALLOCATABLE:: CVEC(:,:)  !screenParm
      REAL(8)               :: CVECSUM(LMXX)
      REAL(8)               :: R2(3)
      REAL(8)   ,ALLOCATABLE:: ORB(:,:,:,:),resO(:,:) !locOrb.
      REAL(8)   ,ALLOCATABLE:: ORBI(:,:,:,:) 
      REAL(8)   ,ALLOCATABLE:: QBARVEC(:,:)
      INTEGER(4),ALLOCATABLE:: LOX(:)
      INTEGER(4),ALLOCATABLE:: ISCATT(:)
      LOGICAL               :: TONSITE      !add head function onsite
      LOGICAL               :: TSPHERE
      INTEGER(4)            :: NFIL
      CHARACTER(64)         :: FILE
      CHARACTER(15)         :: STRING         !contains IATO
!     **************************************************************************
      TVEC(:,:)=0.D0
      TVEC(1,1)=2*7.17D0 !2*CaMnO3
      TVEC(2,2)=2*7.17D0
      TVEC(3,3)=2*7.17D0
!
!     ==========================================================================
!     == DEFINE CENTERED GRID                                                 ==
!     ==========================================================================
      IF(N1.EQ.1)TVEC(:,1)=0.D0
      IF(N2.EQ.1)TVEC(:,2)=0.D0
      IF(N3.EQ.1)TVEC(:,3)=0.D0
      TLITTLE=TVEC
      IF(N1.GT.1)TLITTLE(:,1)=TLITTLE(:,1)/REAL(N1-1,KIND=8)
      IF(N2.GT.1)TLITTLE(:,2)=TLITTLE(:,2)/REAL(N2-1,KIND=8)
      IF(N3.GT.1)TLITTLE(:,3)=TLITTLE(:,3)/REAL(N3-1,KIND=8)
      ORIGIN(:)=-0.5D0*(TVEC(:,1)+TVEC(:,2)+TVEC(:,3)) !also 2D
!
      CALL CELL$GETR8A('T0',9,RBAS)
      CALL ATOMLIST$NATOM(NAT)
      ALLOCATE(R0(3,NAT))
      CALL ATOMLIST$GETR8A('R(0)',0,3*NAT,R0)
      ALLOCATE(ISPECIES1(NAT))
      CALL ATOMLIST$GETI4A('ISPECIES',0,NAT,ISPECIES1)
      ALLOCATE(RCOV(NAT))
      ALLOCATE(QBARVEC(LMXX,NAT))
      QBARVEC(:,:)=0.D0
      LM1X=0
      DO IAT=1,NAT
        ISP=ISPECIES1(IAT)
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETR8('AEZ',AEZ)
        CALL PERIODICTABLE$GET(AEZ,'R(COV)',RCOV(IAT))
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
      ORIGIN(:)=ORIGIN(:)+R0(:,IAT0)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      ALLOCATE(CVEC(LMXX,LM1X))
      ALLOCATE(ORB(N1,N2,N3,LM1X))
      ALLOCATE(resO(N1*N2*N3,LM1X))
      ALLOCATE(ORBI(N1,N2,N3,LM1X))
      ORB(:,:,:,:)=0.D0
      ORBI(:,:,:,:)=0.D0
      NNB=SIZE(SBAR)  !SBAR strucCons. (global)
      DO NN=1,NNB
         IF(SBAR(NN)%IAT1.NE.IAT0) CYCLE
         IF(SBAR(NN)%N1.NE.LM1X) THEN
           CALL ERROR$MSG('INTERNAL ERROR')
           CALL ERROR$STOP('LMTO$PLOTLOCORB')
         END IF
         IAT2=SBAR(NN)%IAT2
         LM2X=SBAR(NN)%N2
         IF(LM2X.GT.LMXX) THEN
           CALL ERROR$MSG('LM2X GREATER THAN LMXX')
           CALL ERROR$MSG('INCREASE PARAMETER LMXX')
           CALL ERROR$STOP('LMTO_PLOTLOCORB')
         END IF
         R2(:)=R0(:,IAT2)+RBAS(:,1)*REAL(SBAR(NN)%IT(1),KIND=8) &
     &                   +RBAS(:,2)*REAL(SBAR(NN)%IT(2),KIND=8) &
     &                   +RBAS(:,3)*REAL(SBAR(NN)%IT(3),KIND=8) 
         TONSITE=(IAT2.EQ.IAT0).AND.(MAXVAL(ABS(SBAR(NN)%IT(:))).EQ.0)
         DO LM1=1,LM1X
           CVEC(:LM2X,LM1)=QBARVEC(:LM2X,IAT2)*SBAR(NN)%MAT(:,LM1)
           IF(TONSITE)CVEC(LM1,LM1)=CVEC(LM1,LM1)+1.D0
         ENDDO
         DO I1=1,N1
           XI(1)=REAL(I1-1,KIND=8)
           DO I2=1,N2
             XI(2)=REAL(I2-1,KIND=8)
             DO I3=1,N3
               XI(3)=REAL(I3-1,KIND=8)
!              == DR IS THE DISTANCE FROM SECOND ATOM ==========================
               DR(:)=ORIGIN(:)+MATMUL(TLITTLE,XI)-R2(:)
               TSPHERE=(DOT_PRODUCT(DR,DR).LT.RCOV(IAT2))
!
               CALL  LMTO$SOLIDHANKEL(DR,RCOV(IAT2),K2,LM2X,K0(1:LM2X))
               IF(TSPHERE) THEN
                 CALL  LMTO$SOLIDBESSEL(DR,K2,LM2X,J0(1:LM2X))
                 JBAR(:LM2X)=J0(:LM2X)-K0(:LM2X)*QBARVEC(:LM2X,IAT2)
               END IF
               DO LM1=1,LM1X
!                == |KBAR>=|K0>*(1+QBAR*SBAR) ==================================
                 ORB(I1,I2,I3,LM1)=ORB(I1,I2,I3,LM1) &
      &                        +DOT_PRODUCT(K0(1:LM2X),CVEC(1:LM2X,LM1))
!                == -|JBAR>SBAR=-(|J0>-|K0>QBAR)*SBAR ==========================
                 IF(TSPHERE) THEN
                   ORBI(I1,I2,I3,LM1)=ORBI(I1,I2,I3,LM1) &
      &                        -DOT_PRODUCT(JBAR(1:LM2X),SBAR(NN)%MAT(:,LM1))
                   IF(TONSITE) ORBI(I1,I2,I3,LM1)=ORBI(I1,I2,I3,LM1)+K0(LM1)
                 END IF
               ENDDO
             ENDDO
           ENDDO
         ENDDO
      ENDDO
      DEALLOCATE(CVEC)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      FILE='TEST_LOCORB'
      WRITE(STRING,*)IAT0 
      FILE=TRIM(ADJUSTL(FILE))//'_IAT'//TRIM(ADJUSTL(STRING))//'.DAT'
      CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
      CALL FILEHANDLER$UNIT('HOOK',NFIL)
      If(1.eq.1) then
!       == writes file for the data explorer ===================================
        WRITE(NFIL,FMT='("class gridpositions counts ",3(1X,I3))')N1,N2,N3
        WRITE(NFIL,FMT='("origin ",3F10.7)')ORIGIN
        WRITE(NFIL,FMT='("#LM1X= ",I2)')LM1X                !dx/2Dcntr 
        WRITE(NFIL,FMT='("#delta(A) ",3F10.7)')TVEC(1,:)/N1 !dx/2Dcntr
        WRITE(NFIL,FMT='("#delta(A) ",3F10.7)')TVEC(2,:)/N2 !dx/2Dcntr
        WRITE(NFIL,FMT='("#delta(A) ",3F10.7)')TVEC(3,:)/N3 !dx/2Dcntr
        resO=RESHAPE(ORB,(/I1*I2*I3,LM1X/))
        PRINT*,"start write data"
        DO LM1=1,LM1X  !dx/2Dcntr
          WRITE(NFIL,FMT='("#LM1= ",I2)')LM1 !dx/2Dcntr
          DO I1=1,N1*N2*N3,10 
             WRITE(NFIL,FMT='(10F10.5)')resO(I1,LM1),resO(I1+1,LM1) &
     &            ,resO(I1+2,LM1),resO(I1+3,LM1), resO(I1+4,LM1) &
     &            ,resO(I1+5,LM1),resO(I1+6,LM1),resO(I1+7,LM1) &
     &            ,resO(I1+8,LM1),resO(I1+9,LM1)
          ENDDO 
        ENDDO
        PRINT*,"done write data"
      else
        WRITE(NFIL,*)ORIGIN(:)
        WRITE(NFIL,*)TVEC(:,:)
        WRITE(NFIL,*)N1,N2,N3,LM1X
        NN=0
        DO I1=1,N1
          DO I2=1,N2
            DO I3=1,N3
              NN=NN+1
!              WRITE(NFIL,*)NN,ORB(I1,I2,I3),BAREORB(I1,I2,I3),ORBI(I1,I2,I3)
              WRITE(NFIL,FMT='(3I5,100F10.5)')I1,I2,I3,ORB(I1,I2,I3,:) &
     &                                             ,ORBI(I1,I2,I3,:)
            ENDDO
          ENDDO
        ENDDO
      end if
      CALL FILEHANDLER$CLOSE('HOOK')
      CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOKERROR')

!CALL SPECIALFUNCTION$TEST()
!CALL TEST_LMTO$STRUCTURECONSTANTS()
      STOP 'FORCED IN LMTO_PLOTLOCORB'
      RETURN
      END
