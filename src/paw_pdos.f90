!
!.......................................................................
MODULE PDOS_MODULE
TYPE STATE_TYPE
 INTEGER(4)         :: NB
 REAL(8)   ,POINTER :: EIG(:)
 REAL(8)   ,POINTER :: occ(:)
 COMPLEX(8),POINTER :: VEC(:,:,:)
END TYPE STATE_TYPE 
INTEGER(4)             :: NAT
INTEGER(4)             :: NSP
INTEGER(4)             :: NKPT
INTEGER(4)             :: NSPIN
INTEGER(4)             :: NDIM
INTEGER(4)             :: NPRO
INTEGER(4)             :: LNXX
REAL(8)                :: RBAS(3,3)
REAL(8)   ,ALLOCATABLE :: R(:,:)
INTEGER(4),ALLOCATABLE :: LNX(:)
INTEGER(4),ALLOCATABLE :: LOX(:,:)
INTEGER(4),ALLOCATABLE :: ISPECIES(:)
CHARACTER(16),ALLOCATABLE :: ATOMID(:)
INTEGER(4),ALLOCATABLE :: IZ(:)
REAL(8)   ,ALLOCATABLE :: RAD(:)
REAL(8)   ,ALLOCATABLE :: PHIOFR(:,:)
REAL(8)   ,ALLOCATABLE :: DPHIDR(:,:)
REAL(8)   ,ALLOCATABLE :: OV(:,:,:)
REAL(8)   ,ALLOCATABLE :: XK(:,:)
TYPE(STATE_TYPE),ALLOCATABLE,TARGET :: STATEARR(:,:)
TYPE(STATE_TYPE),POINTER :: STATE
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
      ELSE IF(ID.EQ.'NSPIN') THEN
        NSPIN=VAL
      ELSE IF(ID.EQ.'NDIM') THEN
        NDIM=VAL
      ELSE IF(ID.EQ.'LNXX') THEN
        LNXX=VAL
      ELSE IF(ID.EQ.'NPRO') THEN
        NPRO=VAL
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
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$GETI4A')
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
        VAL=RESHAPE(RAD,(/NAT/))
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
!clemens : ifc compiler bugfix
        rbas(:,1)=val(1:3)
        rbas(:,2)=val(4:6)
        rbas(:,3)=val(7:9)
      ELSE IF(ID.EQ.'R') THEN
        IF(NAT.EQ.0) NAT=LEN/3
        IF(LEN.NE.3*NAT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(R)) ALLOCATE(R(3,NAT))
        R=RESHAPE(VAL,(/3,NAT/))
      ELSE IF(ID.EQ.'RAD') THEN
        IF(NSP.EQ.0) NSP=LEN
        IF(LEN.NE.NSP) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(RAD)) ALLOCATE(RAD(NSP))
        RAD=VAL
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
      ELSE IF(ID.EQ.'XK') THEN
        IF(NKPT.EQ.0) NKPT=LEN/3
        IF(LEN.NE.3*NKPT) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('PDOS$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(XK)) ALLOCATE(XK(3,NKPT))
        XK=RESHAPE(VAL,(/3,NKPT/))
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('PDOS$SETR8A')
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE PDOS$READ(NFIL)
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE PDOS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4)             :: ISP,IKPT,ISPIN,IB
      INTEGER(4)             :: LNX1,NB
      INTEGER(4)             :: NTASKS,THISTASK
      INTEGER(4)             :: IOS
      CHARACTER(82)          :: IOSTATMSG
      character(32)          :: flag
      logical(4)             :: tchk
!     ******************************************************************
                             CALL TRACE$PUSH('PDOS%READ')
      CALL MPE$QUERY(NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
!
!     ==================================================================
!     == GENERAL QUANTITIES                                           ==
!     ==================================================================
      tchk=.false.
      READ(NFIL,ERR=100)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX,flag
      tchk=.true.
 100  CONTINUE
      IF(.NOT.TCHK) THEN
        PRINT*,'WARNING: NO OCCUPATIONS PRESENT IN PDOS FILE'
        PRINT*,'            OCCUPATIONS WILL BE SET TO 0'
        FLAG='OLD VERSION'
        REWIND(NFIL)
        READ(NFIL)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX
      END IF
      ALLOCATE(LNX(NSP))
      ALLOCATE(LOX(LNXX,NSP))
      ALLOCATE(ISPECIES(NAT))
      READ(NFIL)LNX(:),LOX(:,:),ISPECIES(:)
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
!     ==  NOW read PROJECTIONS                                       ==
!     ==================================================================
      ALLOCATE(XK(3,NKPT))
      ALLOCATE(STATEARR(NKPT,NSPIN))
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          READ(NFIL)XK(:,IKPT),NB
          STATE=>STATEARR(IKPT,ISPIN)
          STATE%NB=NB
          ALLOCATE(STATE%EIG(NB))
          ALLOCATE(STATE%VEC(NDIM,NPRO,NB))
          ALLOCATE(state%OCC(NB))
          DO IB=1,NB
            IF(FLAG.EQ.'011004') THEN
              READ(NFIL,ERR=9999,IOSTAT=IOS)STATE%EIG(IB) &
    &                          ,state%OCC(IB),STATE%VEC(:,:,IB)
            ELSE
              state%occ(:)=0.d0
              READ(NFIL,ERR=9999,IOSTAT=IOS)STATE%EIG(IB),STATE%VEC(:,:,IB)
            END IF
          ENDDO
        ENDDO
      ENDDO
                             CALL TRACE$POP
      RETURN
 9999 CONTINUE
      CALL FILEHANDLER$IOSTATMESSAGE(IOS,IOSTATMSG)
      CALL ERROR$MSG('ERROR reading pdos FILE')
      CALL ERROR$I4VAL('IOS',IOS)
      CALL ERROR$CHVAL('IOSTATMSG',IOSTATMSG)
      CALL ERROR$I4VAL('IB',IB)
      CALL ERROR$I4VAL('IKPT',IKPT)
      CALL ERROR$I4VAL('ISPIN',ISPIN)
      CALL ERROR$I4VAL('NPRO',NPRO)
      CALL ERROR$STOP('PDOS$READ')
      END SUBROUTINE PDOS$READ
!
!     ..................................................................
      SUBROUTINE PDOS$WRITE(NFIL)
!     ******************************************************************
!     **  WAVES$GET                                                   **
!     **  GET                                                         **
!     ******************************************************************
      USE PDOS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      INTEGER(4)             :: ISP,IKPT,ISPIN,IB
      INTEGER(4)             :: LNX1,NB
      INTEGER(4)             :: NTASKS,THISTASK
!     ******************************************************************
                             CALL TRACE$PUSH('PDOS$WRITE')
      CALL MPE$QUERY(NTASKS,THISTASK)
      IF(THISTASK.NE.1) RETURN
!
!     ==================================================================
!     == GENERAL QUANTITIES                                           ==
!     ==================================================================
      WRITE(NFIL)NAT,NSP,NKPT,NSPIN,NDIM,NPRO,LNXX
      WRITE(NFIL)LNX(:),LOX(:,:),ISPECIES(:)
!
!     ==================================================================
!     == ATOMIC STRUCTURE                                             ==
!     ==================================================================
      WRITE(NFIL)RBAS(:,:),R(:,:),ATOMID(:)
!
!     ==================================================================
!     == ELEMENT SPECIFIC QUANTITIES                                  ==
!     ==================================================================
      DO ISP=1,NSP
        LNX1=LNX(ISP)
        WRITE(NFIL)IZ(ISP),RAD(ISP),PHIOFR(1:LNX1,ISP) &
     &            ,DPHIDR(1:LNX1,ISP),OV(1:LNX1,1:LNX1,ISP)
      ENDDO
!
!     ==================================================================
!     ==  NOW WRITE PROJECTIONS                                       ==
!     ==================================================================
      DO IKPT=1,NKPT
        DO ISPIN=1,NSPIN
          STATE=>STATEARR(IKPT,ISPIN)
          NB=STATE%NB
          WRITE(NFIL)XK(:,IKPT),NB
          DO IB=1,NB
            WRITE(NFIL)STATE%EIG(NB),STATE%VEC(:,:,IB)
          ENDDO
        ENDDO
      ENDDO
                             CALL TRACE$POP
      RETURN
      END

