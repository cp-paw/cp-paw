MODULE GW_MODULE
INTEGER(4)             :: NSPIN
INTEGER(4)             :: NRL
REAL(8)   ,ALLOCATABLE :: PSPOT(:,:)    !(NR,NSPIN)
REAL(8)   ,ALLOCATABLE :: H1C(:,:,:,:)    !(LMNXX,LMNXX,NSPIN,NAT)
REAL(8)   ,ALLOCATABLE :: O1C(:,:,:,:)    !(LMNXX,LMNXX,NSPIN,NAT)
REAL(8)                :: NAT
REAL(8)                :: IATP
REAL(8)                :: LMNXX
LOGICAL(4)             :: TON=.FALSE.
END MODULE GW_MODULE
!     ..................................................................
      SUBROUTINE GW$SETI4(ID,,VAL)
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     ******************************************************************
      IF(ID.EQ.'NRL') THEN
        NRL=VAL
      ELSE IF(ID.EQ.'NSPIN') THEN
        NSPIN=VAL
      ELSE IF(ID.EQ.'NAT') THEN
        NSPIN=VAL
      ELSE IF(ID.EQ.'ATOMSELECTOR') THEN
        IATP=VAL
      ELSE IF(ID.EQ.'LMNXX') THEN
        LMNXX=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('GW$SETRI4')
      END IF
      RETURN
      END
!     ..................................................................
      SUBROUTINE GW$SETR8A(ID,LEN,VAL)
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      REAL(8)     ,INTENT(IN) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'PSPOT') THEN
        IF(.NOT.TON) RETURN
        IF(LEN.NE.NRL*NSPIN) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('LEN',LEN)
          CALL ERROR$CHVAL('LEN EXPECTED',NRL*NSPIN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GW$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(PSPOT))ALLOCATE(PSPOT(NRL,NSPIN))
        PSPOT=RESHAPE(VAL,(/NRL,NSPIN/))
      ELSE IF(ID.EQ.'1C-HAMILTON') THEN
        IF(.NOT.TON) RETURN
        IF(LEN.NE.LMNXX*NSPIN*NAT.LE.0) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('LEN',LEN)
          CALL ERROR$CHVAL('LMNXX',LMNXX)
          CALL ERROR$CHVAL('NSPIN',NSPIN)
          CALL ERROR$CHVAL('NAT',NAT)
          CALL ERROR$CHVAL('LEN EXPECTED',LMNXX*LMNXX*NSPIN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GW$SETR8A')
        END IF
        IF(IATP.EQ.0) THEN
          CALL ERROR$MSG('NO ATOM SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GW$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(H1C)) THEN
          ALLOCATE(H1C(LMNXX,LMNXX,NSPIN,NAT)
          H1C=0.D0
        END IF
        H1C(:,:,:,IATP)=RESHAPE(VAL,(/LMNXX,LMNXX,NSPIN/))
      ELSE IF(ID.EQ.'1C-OVERLAP') THEN
        IF(.NOT.TON) RETURN
        IF(LEN.NE.LMNXX*NSPIN*NAT.LE.0) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('LEN',LEN)
          CALL ERROR$CHVAL('LMNXX',LMNXX)
          CALL ERROR$CHVAL('NSPIN',NSPIN)
          CALL ERROR$CHVAL('NAT',NAT)
          CALL ERROR$CHVAL('LEN EXPECTED',LMNXX*LMNXX*NSPIN)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GW$SETR8A')
        END IF
        IF(IATP.EQ.0) THEN
          CALL ERROR$MSG('NO ATOM SELECTED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('GW$SETR8A')
        END IF
        IF(.NOT.ALLOCATED(O1C)) THEN
          ALLOCATE(O1C(LMNXX,LMNXX,NSPIN,NAT)
          O1C=0.D0
        END IF
        O1C(:,:,:,IATP)=RESHAPE(VAL,(/LMNXX,LMNXX,NSPIN/))
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('GW$SETR8A')
      END IF
      RETURN
      END
!     ..................................................................
      SUBROUTINE GW$INPUTCNTL(LL_CNTL_)
      USE LLIST_MODULE
      TYPE(LL_TYPE),INTENT(IN) :: LL_CNTL_
      TYPE(LL_TYPE)            :: LL_CNTL
      LOGICAL(4)               :: TCHK
!     ******************************************************************
      LL_CNTL=LL_CNTL_
      CALL LINKEDLIST$EXISTD(LL_CNTL,'GW',TON)
      IF(.NOT.TON) RETURN
      CALL LINKEDLIST$SELECT(LL_CNTL,'GW')
      RETURN
      END
