!
!.......................................................................
!     ******************************************************************
!     ** THERE ARE DIFFERENT VERSIONS OF THE PDOS FILE:               **
!     **   FLAG=LEGACY                                                ** 
!     **        - NO FLAG PRESENT IN PDOS FILE, BT SET IN PDOS$READ   **
!     **        - PRODUCED BY PAW DIECTLY                             **
!     **        - HAS NO OCCUPATIONS                                  **
!     **        - NOT COMPATIBLE WITH TETRAHEDRON METHOD              **    
!     **   FLAG=011004                                                **
!     **        - PRODUCED BY PAW DIECTLY                             **
!     **        - HAS OCCUPATIONS                                     **
!     **        - NOT COMPATIBLE WITH TETRAHEDRON METHOD              **    
!     **   FLAG=181213                                                **
!     **        - PRODUCED BY KPOINT-DIAGONALISATION WITH PAW_BANDS OR**
!     **             OR PAW DIRECTLY                                  **
!     **        - HAS OCCUPATIONS                                     **
!     **        - COMPATIBLE WITH TETRAHEDRON METHOD                  **
!     **        - HAS SYMMETRY INFORMATIONS (SPACEGROUP)              **
!     ******************************************************************

      MODULE PDOS_MODULE
        TYPE STATE_TYPE
          INTEGER(4)                             :: NB
          REAL(8)   ,POINTER                     :: EIG(:)
          REAL(8)   ,POINTER                     :: OCC(:)
          COMPLEX(8),POINTER                     :: VEC(:,:,:)
        END TYPE STATE_TYPE 
        CHARACTER(6)                             :: FLAG
        INTEGER(4)                               :: NAT
        INTEGER(4)                               :: NSP
        INTEGER(4)                               :: NKPT
        INTEGER(4)                               :: NSPIN
        INTEGER(4)                               :: NDIM
        INTEGER(4)                               :: NPRO
        INTEGER(4)                               :: NKDIV(3)
        INTEGER(4)                               :: ISHIFT(3)
        REAL(8)                                  :: RNTOT
        REAL(8)                                  :: NEL
        LOGICAL(4)                               :: TINV
        INTEGER(4)                               :: LNXX
        REAL(8)                                  :: RBAS(3,3)
        REAL(8)   ,ALLOCATABLE                   :: R(:,:)
        INTEGER(4),ALLOCATABLE                   :: LNX(:)
        INTEGER(4),ALLOCATABLE                   :: LOX(:,:)
        INTEGER(4),ALLOCATABLE                   :: ISPECIES(:)
        CHARACTER(16),ALLOCATABLE                :: ATOMID(:)
        INTEGER(4),ALLOCATABLE                   :: IZ(:)
        REAL(8)   ,ALLOCATABLE                   :: RAD(:)
        REAL(8)   ,ALLOCATABLE                   :: PHIOFR(:,:)
        REAL(8)   ,ALLOCATABLE                   :: DPHIDR(:,:)
        REAL(8)   ,ALLOCATABLE                   :: OV(:,:,:)
        REAL(8)   ,ALLOCATABLE                   :: XK(:,:)
        REAL(8)   ,ALLOCATABLE                   :: WKPT(:)
        TYPE(STATE_TYPE),ALLOCATABLE,TARGET      :: STATEARR(:,:)
        TYPE(STATE_TYPE),POINTER                 :: STATE
        INTEGER(4)                               :: SPACEGROUP
        LOGICAL(4)                               :: TSHIFT
      END MODULE PDOS_MODULE
!
!     ..................................................................
      SUBROUTINE PDOS$GETI4(ID,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'NAT') THEN
        VAL=NAT
      ELSE IF(ID.EQ.'NSP') THEN
        VAL=NSP
      ELSE IF(ID.EQ.'NKPT') THEN
        VAL=NKPT
      ELSE IF(ID.EQ.'NSPIN') THEN
        VAL=NSPIN
      ELSE IF(ID.EQ.'NDIM') THEN
        VAL=NDIM
      ELSE IF(ID.EQ.'LNXX') THEN
        VAL=LNXX
      ELSE IF(ID.EQ.'NPRO') THEN
        VAL=NPRO
      ELSE IF(ID.EQ.'SPACEGROUP') THEN
        VAL=SPACEGROUP
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$GETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$SETI4(ID,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'NAT') THEN
        NAT=VAL
      ELSE IF(ID.EQ.'NSP') THEN
        NSP=VAL
      ELSE IF(ID.EQ.'NKPT') THEN
        NKPT=VAL
        IF(ALLOCATED(XK))DEALLOCATE(XK)
      ELSE IF(ID.EQ.'NSPIN') THEN
        NSPIN=VAL
      ELSE IF(ID.EQ.'NDIM') THEN
        NDIM=VAL
      ELSE IF(ID.EQ.'LNXX') THEN
        LNXX=VAL
      ELSE IF(ID.EQ.'NPRO') THEN
        NPRO=VAL
      ELSE IF(ID.EQ.'SPACEGROUP') THEN
        SPACEGROUP=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$SETI4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$SETI4A(ID,LEN,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: VAL(LEN)
      INTEGER(4)              :: ISPIN,IKPT,I
!     ******************************************************************
      IF(ID.EQ.'IZ') THEN
        IF(NSP.EQ.0)NSP=LEN
        IF(LEN.NE.NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NSP',NSP)
          CALL ERROR$STOP('PDOS$SETI4')
        END IF
        IF(.NOT.ALLOCATED(IZ))ALLOCATE(IZ(NSP))
        IZ=VAL
      ELSE IF(ID.EQ.'LNX') THEN
        IF(NSP.EQ.0) NSP=LEN
        IF(LEN.NE.NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NSP',NSP)
          CALL ERROR$STOP('PDOS$SETI4')
        END IF
        IF(.NOT.ALLOCATED(LNX))ALLOCATE(LNX(NSP))
        LNX=VAL
      ELSE IF(ID.EQ.'LOX') THEN
        IF(NSP*LNXX.EQ.0) THEN
          IF(NSP.NE.0) LNXX=LEN/NSP
          IF(LNXX.NE.0) NSP=LEN/LNXX
        END IF
        IF(LEN.NE.NSP*LNXX) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NSP',NSP)
          CALL ERROR$I4VAL('LNXX',LNXX)
          CALL ERROR$STOP('PDOS$SETI4')
        END IF
        IF(.NOT.ALLOCATED(LOX))ALLOCATE(LOX(LNXX,NSP))
        LOX=RESHAPE(VAL,(/LNXX,NSP/))
      ELSE IF(ID.EQ.'ISPECIES') THEN
        IF(NAT.EQ.0) NAT=LEN
        IF(LEN.NE.NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NAT',NAT)
          CALL ERROR$STOP('PDOS$SETI4')
        END IF
        IF(.NOT.ALLOCATED(ISPECIES))ALLOCATE(ISPECIES(NAT))
        ISPECIES=VAL
      ELSE IF(ID.EQ.'NB') THEN
        IF(NKPT*NSPIN.EQ.0) THEN
          IF(NKPT.NE.0) NSPIN=LEN/NKPT
          IF(NSPIN.NE.0) NKPT=LEN/NSPIN
        END IF
        IF(LEN.NE.NKPT*NSPIN) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETI4')
        END IF
        IF(.NOT.ALLOCATED(STATEARR))ALLOCATE(STATEARR(NKPT,NSPIN))
        I=0
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            I=I+1
            STATEARR(IKPT,ISPIN)%NB=VAL(I)
          ENDDO
        ENDDO
      ELSE IF(ID.EQ.'NKDIV') THEN
        IF(NKPT.EQ.0) NKPT=VAL(1)*VAL(2)*VAL(3)
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('PDOS$SETI4')
        END IF
        NKDIV(1:3)=VAL(1:3)
      ELSE IF(ID.EQ.'ISHIFT') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('PDOS$SETI4')
        END IF
        ISHIFT(1:3)=VAL(1:3)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$SETI4A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$GETI4A(ID,LEN,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(OUT):: VAL(LEN)
      INTEGER(4)              :: ISPIN,IKPT,I
!     ******************************************************************
      IF(ID.EQ.'LNX') THEN
        IF(LEN.NE.NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$I4VAL('NSP',NSP)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        IF(.NOT.ALLOCATED(LNX)) THEN
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        VAL=LNX
      ELSE IF(ID.EQ.'LOX') THEN
        IF(LEN.NE.NSP*LNXX) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        IF(.NOT.ALLOCATED(LOX)) THEN 
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        VAL=RESHAPE(LOX,(/LNXX*NSP/))
      ELSE IF(ID.EQ.'ISPECIES') THEN
        IF(LEN.NE.NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        IF(.NOT.ALLOCATED(ISPECIES)) THEN
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        VAL=ISPECIES
      ELSE IF(ID.EQ.'NB') THEN
        IF(NKPT*NSPIN.EQ.0) THEN
          CALL ERROR$MSG('NKPT AND NSPIN MUST BE DEFINED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('NKPT',NKPT)
          CALL ERROR$I4VAL('NSPIN',NSPIN)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        IF(LEN.NE.NKPT*NSPIN) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        IF(.NOT.ALLOCATED(STATEARR)) THEN
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        I=0
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            I=I+1
            VAL(I)=STATEARR(IKPT,ISPIN)%NB
          ENDDO
        ENDDO
      ELSE IF(ID.EQ.'NKDIV') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        VAL=NKDIV
      ELSE IF(ID.EQ.'ISHIFT') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETI4')
        END IF
        VAL=ISHIFT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$GETI4A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$GETL4(ID,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'TINV') THEN
        VAL=TINV
      ELSE IF(ID.EQ.'TSHIFT') THEN
        VAL=TSHIFT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$GETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$SETL4(ID,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'TINV') THEN
        TINV=VAL
      ELSE IF(ID.EQ.'TSHIFT') THEN
        TSHIFT=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$SETL4')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$GETR8(ID,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'RNTOT') THEN
        VAL=RNTOT
      ELSE IF(ID.EQ.'NEL') THEN
        VAL=NEL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$GETR84')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$GETCH(ID,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(6),INTENT(OUT):: VAL
!     ******************************************************************
      IF(ID.EQ.'FLAG') THEN
        VAL=FLAG
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$GETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$SETR8(ID,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'RNTOT') THEN
        RNTOT=VAL
      ELSE IF(ID.EQ.'NEL') THEN
        NEL=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$SETR8')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$GETR8A(ID,LEN,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(OUT):: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'RBAS') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        VAL=RESHAPE(RBAS,(/9/))
      ELSE IF(ID.EQ.'R') THEN
        IF(LEN.NE.3*NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        IF(.NOT.ALLOCATED(R)) THEN
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        VAL=RESHAPE(R,(/3*NAT/))
      ELSE IF(ID.EQ.'RAD') THEN
        IF(LEN.NE.NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        IF(.NOT.ALLOCATED(RAD)) THEN
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        VAL=RESHAPE(RAD,(/NSP/))
      ELSE IF(ID.EQ.'PHI') THEN
        IF(LEN.NE.LNXX*NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        IF(.NOT.ALLOCATED(PHIOFR)) THEN
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        VAL=RESHAPE(PHIOFR,(/LNXX*NSP/))
      ELSE IF(ID.EQ.'DPHIDR') THEN
        IF(LEN.NE.LNXX*NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        IF(.NOT.ALLOCATED(DPHIDR)) THEN
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        VAL=RESHAPE(DPHIDR,(/LNXX*NSP/))
      ELSE IF(ID.EQ.'OVERLAP') THEN
        IF(LEN.NE.LNXX*LNXX*NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        IF(.NOT.ALLOCATED(OV)) THEN
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        VAL=RESHAPE(OV,(/LNXX*LNXX*NSP/))
!
      ELSE IF(ID.EQ.'XK') THEN
        IF(LEN.NE.3*NKPT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        IF(.NOT.ALLOCATED(XK)) THEN
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        VAL=RESHAPE(XK,(/3*NKPT/))
!
      ELSE IF(ID.EQ.'WKPT') THEN
        IF(LEN.NE.NKPT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        IF(.NOT.ALLOCATED(WKPT)) THEN
          CALL ERROR$MSG('REQUESTED QUANTITY HAS NOT BEEN SET')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETR8A')
        END IF
        VAL=WKPT(:)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$GETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$SETR8A(ID,LEN,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'RBAS') THEN
        IF(LEN.NE.9) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
!        RBAS=RESHAPE(VAL,(/3,3/))
!CLEMENS : IFC COMPILER BUGFIX
        RBAS(:,1)=VAL(1:3)
        RBAS(:,2)=VAL(4:6)
        RBAS(:,3)=VAL(7:9)
      ELSE IF(ID.EQ.'R') THEN
        IF(NAT.EQ.0) NAT=LEN/3
        IF(LEN.NE.3*NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(R)) ALLOCATE(R(3,NAT))
        R=RESHAPE(VAL,(/3,NAT/))
!
      ELSE IF(ID.EQ.'RAD') THEN
        IF(NSP.EQ.0) NSP=LEN
        IF(LEN.NE.NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(RAD)) ALLOCATE(RAD(NSP))
        RAD=VAL
!
      ELSE IF(ID.EQ.'PHI') THEN
        IF(NSP*LNXX.EQ.0) THEN
          IF(NSP.EQ.0)  NSP=LEN/LNXX
          IF(LNXX.EQ.0) LNXX=LEN/NSP
        END IF  
        IF(LEN.NE.LNXX*NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(PHIOFR)) ALLOCATE(PHIOFR(LNXX,NSP))
        PHIOFR=RESHAPE(VAL,(/LNXX,NSP/))
!
      ELSE IF(ID.EQ.'DPHIDR') THEN
        IF(NSP*LNXX.EQ.0) THEN
          IF(NSP.EQ.0)  NSP=LEN/LNXX
          IF(LNXX.EQ.0) LNXX=LEN/NSP
        END IF  
        IF(LEN.NE.LNXX*NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(DPHIDR))ALLOCATE(DPHIDR(LNXX,NSP))
        DPHIDR=RESHAPE(VAL,(/LNXX,NSP/))
!
      ELSE IF(ID.EQ.'OVERLAP') THEN
        IF(NSP*LNXX.EQ.0) THEN
          IF(NSP.EQ.0)  NSP=LEN/LNXX**2
          IF(LNXX.EQ.0) LNXX=NINT(SQRT(REAL(LEN/NSP)))
        END IF  
        IF(LEN.NE.LNXX*LNXX*NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(OV))ALLOCATE(OV(LNXX,LNXX,NSP))
        OV=RESHAPE(VAL,(/LNXX,LNXX,NSP/))
!
      ELSE IF(ID.EQ.'XK') THEN
        IF(NKPT.EQ.0) NKPT=LEN/3
        IF(LEN.NE.3*NKPT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(XK)) ALLOCATE(XK(3,NKPT))
        XK=RESHAPE(VAL,(/3,NKPT/))
!
!     == K-POINT WEIGHT ==================================================
      ELSE IF(ID.EQ.'WKPT') THEN
        IF(NKPT.EQ.0) NKPT=LEN
        IF(LEN.NE.NKPT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(WKPT)) ALLOCATE(WKPT(NKPT))
        WKPT=VAL(:)
!
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$SETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$SETCHA(ID,LEN,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      CHARACTER(16),INTENT(IN):: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'ATOMID') THEN
        IF(LEN.NE.NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETCHA')
        END IF
        IF(.NOT.ALLOCATED(ATOMID)) ALLOCATE(ATOMID(NAT))
        ATOMID(:)=VAL(:)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$SETCHA')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$GETCHA(ID,LEN,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN) :: ID
      INTEGER(4)   ,INTENT(IN) :: LEN
      CHARACTER(16),INTENT(OUT):: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'ATOMID') THEN
        IF(LEN.NE.NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETCHA')
        END IF
        IF(.NOT.ALLOCATED(ATOMID)) THEN
          CALL ERROR$MSG('ATOMID NOT ALLOCATED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$GETCHA')
        END IF
        VAL(:)=ATOMID(:)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$GETCHA')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$SETCH(ID,VAL)
      USE PDOS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(6),INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'FLAG') THEN
        FLAG=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$SETCH')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$READ(NFIL)
!     ******************************************************************
!     ** READS PDOS FILE FROM NFIL INTO THE PDOS MODULE               **
!     ******************************************************************
      USE PDOS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4)             :: ISP,IKPT,ISPIN,IB
      INTEGER(4)             :: LNX1,NB
      INTEGER(4)             :: IOS
      CHARACTER(82)          :: IOSTATMSG
      LOGICAL(4)             :: TCHK
      REAL(8)                :: OCCSUM
      INTEGER(4)             :: ILOGICAL
!     ******************************************************************
                             CALL TRACE$PUSH('PDOS$READ')
!
!     ==================================================================
!     == GENERAL QUANTITIES                                           ==
!     ==================================================================
      TCHK=.FALSE.
      READ(NFIL,ERR=100)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,FLAG
      TCHK=.TRUE.
 100  CONTINUE
      IF(.NOT.TCHK) THEN
        PRINT*,'WARNING: NO OCCUPATIONS PRESENT IN PDOS FILE'
        PRINT*,'            OCCUPATIONS WILL BE SET TO 0'
        FLAG='LEGACY'
        REWIND(NFIL)
        READ(NFIL)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX
      END IF
      
      ALLOCATE(LNX(NSP))
      ALLOCATE(LOX(LNXX,NSP))
      ALLOCATE(ISPECIES(NAT))
      READ(NFIL)LNX(:),LOX(:,:),ISPECIES(:)
      
      IF(FLAG.EQ.'181213')THEN
!         == GFORTRAN LOGICAL REPRESENTATION DEFINED WITH TRUE=1, FALSE=0     ==
!         https://gcc.gnu.org/onlinedocs/gfortran/compiler-characteristics/
!         internal-representation-of-logical-variables.html
!         == IFORT LOGICAL REPRESENTATION DEFINED WITH VALUE OF LAST BIT      ==
!         https://www.intel.com/content/www/us/en/docs/fortran-compiler/
!         developer-guide-reference/2024-2/logical-data-representations.html
!         == BOTH SHARE MEANING OF LAST BIT 1=TRUE, 0=FALSE                   ==
!         == ENSURES BACKWARDS COMPATIBILITY WITH OLD PDOS FILES              ==
        READ(NFIL)NKDIV(:),ISHIFT(:),RNTOT,NEL,ILOGICAL
        TINV=BTEST(ILOGICAL,0)
        READ(NFIL)SPACEGROUP,ILOGICAL
        TSHIFT=BTEST(ILOGICAL,0)
      ENDIF
!
!     ==================================================================
!     == ATOMIC STRUCTURE                                             ==
!     ==================================================================
      ALLOCATE(R(3,NAT))
      ALLOCATE(ATOMID(NAT))
      READ(NFIL)RBAS(:,:),R(:,:),ATOMID(:)
!
!     ==================================================================
!     == ELEMENT SPECIFIC QUANTITIES                                  ==
!     ==================================================================
      ALLOCATE(IZ(NSP))
      ALLOCATE(RAD(NSP)); RAD=0.D0
      ALLOCATE(PHIOFR(LNXX,NSP)); PHIOFR=0.D0
      ALLOCATE(DPHIDR(LNXX,NSP)); DPHIDR=0.D0
      ALLOCATE(OV(LNXX,LNXX,NSP)); OV=0.D0
      DO ISP=1,NSP
        LNX1=LNX(ISP)
        READ(NFIL)IZ(ISP),RAD(ISP),PHIOFR(1:LNX1,ISP) &
     &            ,DPHIDR(1:LNX1,ISP),OV(1:LNX1,1:LNX1,ISP)
      ENDDO
!
!     ==================================================================
!     ==  NOW READ PROJECTIONS                                       ==
!     ==================================================================
      OCCSUM=0.0D0
      ALLOCATE(XK(3,NKPT))
      ALLOCATE(WKPT(NKPT))
      ALLOCATE(STATEARR(NKPT,NSPIN))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          READ(NFIL,ERR=9998,END=9998)XK(:,IKPT),NB,WKPT(IKPT)
          STATE%NB=NB
          ALLOCATE(STATE%EIG(NB))
          ALLOCATE(STATE%VEC(NDIM,NPRO,NB))
          ALLOCATE(STATE%OCC(NB))
          DO IB=1,NB
            IF(FLAG.EQ.'LEGACY') THEN
              STATE%OCC(:)=0.D0
              READ(NFIL,ERR=9999,IOSTAT=IOS)STATE%EIG(IB),STATE%VEC(:,:,IB)
            ELSE
              READ(NFIL,ERR=9999,IOSTAT=IOS)STATE%EIG(IB) &
    &                          ,STATE%OCC(IB),STATE%VEC(:,:,IB)
              OCCSUM=OCCSUM+STATE%OCC(IB)
            END IF
          ENDDO
        ENDDO
      ENDDO
PRINT*,"OCCSUM",OCCSUM
                             CALL TRACE$POP
      RETURN
 9999 CONTINUE
      CALL FILEHANDLER$IOSTATMESSAGE(IOS,IOSTATMSG)
      CALL ERROR$MSG('ERROR READING PDOS FILE')
      CALL ERROR$I4VAL('IOS',IOS)
      CALL ERROR$CHVAL('IOSTATMSG',IOSTATMSG)
      CALL ERROR$I4VAL('IB',IB)
      CALL ERROR$I4VAL('IKPT',IKPT)
      CALL ERROR$I4VAL('ISPIN',ISPIN)
      CALL ERROR$I4VAL('NPRO',NPRO)
      CALL ERROR$STOP('PDOS$READ')
      STOP
 9998 CONTINUE
      CALL ERROR$MSG('ERROR READING PDOS FILE')
      CALL ERROR$MSG('OLD VERSION: VARIABLE WKPT IS NOT PRESENT')
      CALL ERROR$MSG('PRODUCE NEW PDOS FILE')
      CALL ERROR$STOP('PDOS$READ')
      STOP
      END SUBROUTINE PDOS$READ
!
!     ..................................................................
      SUBROUTINE PDOS$WRITE(NFIL,FLAG_)
!     ******************************************************************
!     ** WRITES FIRST PART OF PDOS FILE TO NFIL. THE DATA TO BE       **
!     ** WRITTEN, E.G. NAT,NSP,... HAS TO BE TRANSFERED TO THE        **
!     ** PDOS-MODULE BEFORE CALLING THIS FUNCTION. THE KPOINTS,       **
!     ** ENERGIES, PROJECTIONS AND OCCUPATIONS ARE NOT WRITTEN WITH   **
!     ** THIS FUNCTION, BUT WITH PDOS$WRITEK                          **
!     **                                                              **
!     ******************************************************************
      USE PDOS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)        :: NFIL
      CHARACTER(6),INTENT(IN)      :: FLAG_
      INTEGER(4)                   :: ISP,IKPT
      INTEGER(4)                   :: LNX1
      INTEGER(4)                   :: ILOGICAL
!     ******************************************************************
                             CALL TRACE$PUSH('PDOS$WRITE')
      FLAG=FLAG_
!
!     ==================================================================
!     == GENERAL QUANTITIES                                           ==
!     ==================================================================
      WRITE(NFIL)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,FLAG_
      WRITE(NFIL)LNX(:),LOX(:,:),ISPECIES(:)

      IF(FLAG_.EQ.'181213')THEN
!       == PDOS FILE WITH GNU STANDARD OF LOGICAL BIT REPRESENTATION ==
        ILOGICAL=1
        IF(.NOT.TINV)ILOGICAL=0
        WRITE(NFIL)NKDIV(:),ISHIFT(:),RNTOT,NEL,ILOGICAL
        ILOGICAL=1
        IF(.NOT.TSHIFT)ILOGICAL=0
        WRITE(NFIL)SPACEGROUP,ILOGICAL
      ENDIF
!
!     ==================================================================
!     == ATOMIC STRUCTURE                                             ==
!     ==================================================================
      WRITE(NFIL)RBAS,R,ATOMID
!
!     ==================================================================
!     == ELEMENT SPECIFIC QUANTITIES                                  ==
!     ==================================================================
      DO ISP=1,NSP
        LNX1=LNX(ISP)
        WRITE(NFIL)IZ(ISP),RAD(ISP),PHIOFR(1:LNX1,ISP) &
     &       ,DPHIDR(1:LNX1,ISP),OV(1:LNX1,1:LNX1,ISP)
      ENDDO
                             CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$WRITEK(NFIL,XK,NB,NDIM,NPRO,WKPT,EIG,OCC,VECTOR)
!     ******************************************************************
!     **  WRITE EIG,OCC,PROJ FOR ONE KPOINT TO NFIL                   **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)        :: NFIL
      REAL(8)   ,INTENT(IN)        :: XK(3)
      INTEGER(4),INTENT(IN)        :: NB
      INTEGER(4),INTENT(IN)        :: NDIM
      INTEGER(4),INTENT(IN)        :: NPRO
      REAL(8)   ,INTENT(IN)        :: WKPT
      REAL(8)   ,INTENT(IN)        :: EIG(NB)
      REAL(8)   ,INTENT(IN)        :: OCC(NB)
      COMPLEX(8),INTENT(IN)        :: VECTOR(NDIM,NPRO,NB)
      INTEGER(4)                   :: IB1
!     ******************************************************************
                             CALL TRACE$PUSH('PDOS$WRITEK')
      WRITE(NFIL)XK,NB,WKPT
      DO IB1=1,NB
        WRITE(NFIL)EIG(IB1),OCC(IB1),VECTOR(:,:,IB1)
      ENDDO
                             CALL TRACE$POP
      RETURN
      END


