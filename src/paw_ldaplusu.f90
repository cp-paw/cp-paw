!TODO :
! DATH IS STILL REAL AND SHOULD PROBABLY BE COMPLEX LIKE DENMAT
! DENMAT AND DATH ARE BOTH REAL, IF THE CALCULATION IS COLLINEAR
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE LDAPLUSU_MODULE
TYPE THISISP_TYPE
LOGICAL(4)             :: TINI=.FALSE.
LOGICAL(4)             :: TON=.FALSE.
CHARACTER(8)           :: ATOMTYPE=''
INTEGER(4)             :: GID          !GRID ID FOR RADIAL GRID
INTEGER(4)             :: NR           !#(RADIAL GRID POINTS)
INTEGER(4)             :: LNXCHI       !#(RADIAL FUNCTIONS FOR LOCAL ORBITALS)
INTEGER(4),POINTER     :: LOXCHI(:)    !MAIN ANGULAR MOMENTUM OF LOCAL ORBITAL
INTEGER(4),POINTER     :: NORB(:)      !X(# LOCAL FUNCTIONS PER ANGULAR MOMENTUM)
INTEGER(4)             :: LRX          !X(ANGULAR MOMENTUM IN THE DENSITY)
REAL(8)                :: RCUT=0.D0    !RADIUS OF LOCAL ORBITAL
CHARACTER(16)          :: FUNCTIONALID !CAN BE 'LDA+U','LDA+U(OLD)' OR 'HYBRID'
!== SETTINGS SPECIFICALLY FOR HYBRID FUNCTIONAL ================================
REAL(8)                :: HFWEIGHT=0.D0!CONTRIBUTION OF EXACT EXCHANGE
!== SETTINGS SPECIFICALLY FOR LDA+U
INTEGER(4)             :: MAINLN(2)=(/0,0/) !SHELL TO WHICH UPAR AND JPAR REFER TO
LOGICAL(4)             :: USEDIEL=.FALSE.
LOGICAL(4)             :: USEUPAR=.FALSE.
LOGICAL(4)             :: USEJPAR=.FALSE.
LOGICAL(4)             :: USEFRATIO42=.FALSE.
LOGICAL(4)             :: USEFRATIO62=.FALSE.
REAL(8)                :: DIEL=1            !DIELECTRIC CONSTANT
REAL(8)                :: UPAR=0.D0         !U-PARAMETER
REAL(8)                :: JPAR=0.D0         !J-PARAMETER
REAL(8)                :: FRATIO42=0.D0     ! RATIO OF SLATER INTEGRALS F4 AND F2
REAL(8)                :: FRATIO62=0.D0     ! RATIO OF SLATER INTEGRALS F6 AND F2
LOGICAL(4)             :: TCV=.FALSE.
END TYPE THISISP_TYPE
!
TYPE THIS_TYPE
TYPE(THISISP_TYPE),POINTER :: SP
INTEGER(4)                 :: NCHI              !#(LOCAL ORBITALS)
REAL(8)   ,POINTER         :: CHI(:,:)          !LOCAL (HEAD) ORBITALS
REAL(8)   ,POINTER         :: ULITTLE(:,:,:,:,:)!SLATER INTEGRALS 
REAL(8)   ,POINTER         :: DOWNFOLD(:,:)     !MAPS PARTIAL WAVES TO LOCAL ORBITALS
REAL(8)   ,POINTER         :: CVX(:,:)          !CORE VALENCE EXCHANGE 
END TYPE THIS_TYPE
!
TYPE(THIS_TYPE)   ,ALLOCATABLE,TARGET :: THISARRAY(:)
TYPE(THISISP_TYPE),ALLOCATABLE,TARGET :: THISISPARRAY(:)
TYPE(THIS_TYPE)              ,POINTER :: THIS
TYPE(THISISP_TYPE)           ,POINTER :: THISISP
LOGICAL(4)             :: TON=.FALSE.
INTEGER(4)             :: NTHIS=0      ! #(ATOMS)
INTEGER(4)             :: ITHIS=0      ! SELECTED ATOM 
INTEGER(4)             :: NTHISISP=0   ! #(ATOM TYPES)
INTEGER(4)             :: ITHISISP=0   ! SELECTED ATOM TYPE
CHARACTER(8)           :: DCTYPE='FLL' ! SPECIFIES TYPE OF DOUBLE COUNTING CORRECTION
END MODULE LDAPLUSU_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$NEW(NSP)
!     **************************************************************************
!     **                                                                      **
!     **  SETS UP THE ARRAYS FOR OPERATION. THIS ROUTINE MUST BE CALLED ONCE  **
!     **  BEFORE ANY OTHER CALL TO THE LDAPLUSU OBJECT                        **
!     **                                                                      **
!     **  ATTENTION: CHANGES THE SETTING OF SETUP OBJECT                      **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE, ONLY : NTHISISP,THISISPARRAY
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSP    ! #(ATOMS)
      INTEGER(4)            :: ISP
!     **************************************************************************
!
!     ==========================================================================
!     == DO NOTHING IF THISARRAY IS ALREADY ALLOCATED                         ==
!     ==========================================================================
      IF(NTHISISP.NE.0) THEN
        CALL ERROR$MSG('LDAPLUSU$NEW CANNOT BE CALLED TWICE')
        CALL ERROR$STOP('LDAPLUSU$NEW')
      END IF
!
!     ==========================================================================
!     == CREATE THISARRAY                                                     ==
!     ==========================================================================
      NTHISISP=NSP
      ALLOCATE(THISISPARRAY(NTHISISP))
      DO ISP=1,NSP
        ALLOCATE(THISISPARRAY(ISP)%NORB(1))
        THISISPARRAY(ISP)%NORB(:)=0  ! PER DEFAULT CORRELATE NOTHING
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SELECTTYPE(ISP)
!     **************************************************************************
!     **                                                                      **
!     **  SELECT AN ATOM TYPE OR UNSELECT WITH ISP=0                          **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP  ! ATOM-TYPE INDEX
!     **************************************************************************
!
!     ==========================================================================
!     == TEST FOR CONFLICTS                                                   ==
!     ==========================================================================
      IF(NTHISISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT INITIALIZED')
        CALL ERROR$STOP('LDAPLUSU$SELECTTYPE')
      END IF
!
      IF(ISP.NE.0.AND.ITHISISP.NE.0) THEN
        CALL ERROR$MSG('CANNOT SELECT ATOM TYPE WHILE ANOTHER ONE IS SELECTED')
        CALL ERROR$MSG('USE "CALL LDAPLUSU$SELECTTYPE(0)"')
        CALL ERROR$STOP('LDAPLUSU$SELECTTYPE')
      END IF
!
      IF(ITHIS.NE.0) THEN
        CALL ERROR$MSG('CANNOT SELECT ATOM TYPE WHILE AN ATOM IS SELECTED')
        CALL ERROR$MSG('USE "CALL LDAPLUSU$SELECT(0)"')
        CALL ERROR$STOP('LDAPLUSU$SELECTTYPE')
      END IF
!
!     ==========================================================================
!     == MAKE SELECTION                                                       ==
!     ==========================================================================
      ITHISISP=ISP
      IF(ITHISISP.GT.NTHISISP) THEN
        CALL ERROR$MSG('ISP OUT OF RANGE')
        CALL ERROR$STOP('LDAPLUSU$SELECTTYPE')
      END IF
      IF(ITHISISP.EQ.0) THEN
        NULLIFY(THISISP)
      ELSE
        THISISP=>THISISPARRAY(ITHISISP)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SELECT(IAT)
!     **************************************************************************
!     **                                                                      **
!     **  SELECT AN ATOM OR UNSELECT WITH IAT=0                               **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT  ! UNIQUE INDEX OF THE ATOM TYPE
      INTEGER(4)            :: NAT
      INTEGER(4)            :: IAT1
      INTEGER(4)            :: ISP
!     **************************************************************************
!
!     ==========================================================================
!     == INITIALIZATION                                                       ==
!     ==========================================================================
      IF(NTHISISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT INITIALIZED')
        CALL ERROR$STOP('LDAPLUSU$SELECT')
      END IF
!
      IF(NTHIS.EQ.0) THEN
        CALL ATOMLIST$NATOM(NAT)
        NTHIS=NAT
        ALLOCATE(THISARRAY(NTHIS))
        DO IAT1=1,NAT
          CALL ATOMLIST$GETI4('ISPECIES',IAT1,ISP)
          THISARRAY(IAT1)%SP=>THISISPARRAY(ISP)
          THIS=>THISARRAY(IAT1)
          NULLIFY(THIS%CHI)
          NULLIFY(THIS%ULITTLE)
          NULLIFY(THIS%DOWNFOLD)
          NULLIFY(THIS%CVX)
        ENDDO
        ITHIS=0
      END IF
!
!     ==========================================================================
!     == TEST FOR CONFLICTS                                                   ==
!     ==========================================================================
      IF(IAT.NE.0.AND.ITHIS.NE.0) THEN
        CALL ERROR$MSG('CANNOT SELECT ATOM WHILE ANOTHER ONE IS SELECTED')
        CALL ERROR$MSG('USE "CALL LDAPLUSU$SELECT(0)"')
        CALL ERROR$I4VAL('ATOM TO BE SELECTED ',IAT)
        CALL ERROR$I4VAL('ATOM ALREADY SELECTED ',ITHIS)
        CALL ERROR$STOP('LDAPLUSU$SELECT')
      END IF
!
      IF(ITHISISP.NE.0) THEN
        CALL ERROR$MSG('CANNOT SELECT ATOM WHILE AN ATOM TYPE IS SELECTED')
        CALL ERROR$MSG('USE "CALL LDAPLUSU$SELECTTYPE(0)"')
        CALL ERROR$I4VAL('ATOM TO BE SELECTED ',IAT)
        CALL ERROR$I4VAL('ATOM TYPE SELECTED ',ITHISISP)
        CALL ERROR$STOP('LDAPLUSU$SELECT')
      END IF
!
!     ==========================================================================
!     == MAKE SELECTION                                                       ==
!     ==========================================================================
      ITHIS=IAT
      IF(ITHIS.GT.NTHIS) THEN
        CALL ERROR$MSG('IAT OUT OF RANGE')
        CALL ERROR$I4VAL('NTHIS',NTHIS)
        CALL ERROR$I4VAL('IAT',IAT)
        CALL ERROR$STOP('LDAPLUSU$SELECT')
      END IF
      IF(ITHIS.EQ.0) THEN
        NULLIFY(THIS)
      ELSE
        THIS=>THISARRAY(ITHIS)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SETR8(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ITHISISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETR8')
      END IF

      IF(ID.EQ.'UPAR') THEN
        THISISP%UPAR=VAL
        THISISP%USEUPAR=.TRUE.
        IF(THISISP%USEDIEL) THEN
          CALL ERROR$MSG('DIEL HAS ALREADY BEEN SET')
          CALL ERROR$MSG('UPAR AND DIEL CANNOT BE USED SIMULTANEOUSLY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('SELECTION (ATOM TYPE INDEX)',ITHISISP)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
      ELSE IF(ID.EQ.'JPAR') THEN
        THISISP%JPAR=VAL
        THISISP%USEJPAR=.TRUE.
        IF(THISISP%USEDIEL) THEN
          CALL ERROR$MSG('DIEL HAS ALREADY BEEN SET')
          CALL ERROR$MSG('JPAR AND DIEL CANNOT BE USED SIMULTANEOUSLY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('SELECTION (ATOM INDEX)',ITHISISP)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
      ELSE IF(ID.EQ.'F4/F2') THEN
        THISISP%FRATIO42=VAL
        THISISP%USEFRATIO42=.TRUE.
        IF(THISISP%USEDIEL) THEN
          CALL ERROR$MSG('DIEL HAS ALREADY BEEN SET')
          CALL ERROR$MSG('FRATIO42 AND DIEL CANNOT BE USED SIMULTANEOUSLY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('SELECTION (ATOM INDEX)',ITHISISP)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
      ELSE IF(ID.EQ.'F6/F2') THEN
        THISISP%FRATIO62=VAL
        THISISP%USEFRATIO62=.TRUE.
        IF(THISISP%USEDIEL) THEN
          CALL ERROR$MSG('DIEL HAS ALREADY BEEN SET')
          CALL ERROR$MSG('FRATIO62 AND DIEL CANNOT BE USED SIMULTANEOUSLY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('SELECTION (ATOM INDEX)',ITHISISP)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
      ELSE IF(ID.EQ.'DIEL') THEN
        THISISP%DIEL=VAL
        THISISP%USEDIEL=.TRUE.
        IF(THISISP%USEUPAR.OR.THISISP%USEJPAR) THEN
          CALL ERROR$MSG('UPAR OR JPAR HAVE ALREADY BEEN SET')
          CALL ERROR$MSG('DIEL AND UPAR OR JPAR CANNOT BE USED SIMULTANEOUSLY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('SELECTION (ATOM INDEX)',ITHISISP)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
      ELSE IF(ID.EQ.'RCUT') THEN
        THISISP%RCUT=VAL
      ELSE IF(ID.EQ.'HFWEIGHT') THEN
        THISISP%HFWEIGHT=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$I4VAL('SELECTION (ATOM INDEX)',ITHISISP)
        CALL ERROR$STOP('LDAPLUSU$SETR8')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$GETR8(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
!     **************************************************************************
      IF(ITHISISP.EQ.0) THEN
        CALL ERROR$MSG(' ATOM TYPE IN LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$GETR8')
      END IF

      IF(ID.EQ.'RCUT') THEN
        VAL=THIS%SP%RCUT
      ELSE IF(ID.EQ.'HFWEIGHT') THEN
        VAL=THISISP%HFWEIGHT
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$GETR8')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SETL4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'ON') THEN
        TON=VAL
        RETURN
      END IF
!
!     ==========================================================================
!     == SET ATOM SPECIFIC INFORMATION                                        ==
!     ==========================================================================
      IF(ITHISISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETR8')
      END IF
!
      IF(ID.EQ.'ACTIVE') THEN
        THISISP%TON=VAL
      ELSE IF(ID.EQ.'COREVALENCEEXCHANGE') THEN
        THISISP%TCV=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETL4')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$GETL4(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      LOGICAL(4)  ,INTENT(OUT) :: VAL
!     **************************************************************************
!
!     ==========================================================================
!     == SET ATOM SPECIFIC INFORMATION                                        ==
!     ==========================================================================
      IF(ITHISISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$GETR8')
      END IF
!
      IF(ID.EQ.'ACTIVE') THEN
        VAL=TON.AND.THISISP%TON
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$GETL4')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SETI4A(ID,LENG,VAL)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LENG
      INTEGER(4)  ,INTENT(IN) :: VAL(LENG)
!     **************************************************************************
      IF(ITHISISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETI4A')
      END IF
      IF(ID.EQ.'NCORROFL') THEN
        DEALLOCATE(THISISP%NORB)
        ALLOCATE(THISISP%NORB(LENG))
        THISISP%NORB=VAL
!
      ELSE IF(ID.EQ.'MAINLN') THEN
        THISISP%MAINLN=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETI4A')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$GETI4A(ID,LENG,VAL)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LENG
      INTEGER(4)  ,INTENT(OUT) :: VAL(LENG)
      INTEGER(4)               :: ISVAR
!     **************************************************************************
      IF(ITHIS.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETI4A')
      END IF
      IF(ID.EQ.'NCORROFL') THEN
        ISVAR=SIZE(THIS%SP%NORB)
        IF(ISVAR.LT.LENG) THEN
          CALL ERROR$MSG('INSUFFICIENT ARRAY LENGTH')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LENG(EXTERN)',LENG)
          CALL ERROR$I4VAL('LENG(INTERN)',ISVAR)
          CALL ERROR$STOP('LDAPLUSU$GETI4A')
        END IF
        VAL(:)=0
        VAL(1:ISVAR)=THIS%SP%NORB
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$GETI4A')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SETCH(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     **************************************************************************
      IF(ITHISISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETCH')
      END IF
      IF(ID.EQ.'FUNCTIONALID') THEN
        IF(VAL.NE.'LDA+U'.AND.VAL.NE.'LDA+U(OLD)'.AND.VAL.NE.'HYBRID') THEN
          CALL ERROR$MSG('ILLEGAL VALUE FOR FUNCTIONALID')
          CALL ERROR$CHVAL('FUNCTIONALID',VAL)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETCH')
        END IF
        THISISP%FUNCTIONALID=VAL
      ELSE IF(ID.EQ.'DCTYPE') THEN
        DCTYPE=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETCH')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$GETCH(ID,VAL)
!     **************************************************************************
!     **                                                                      **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: ID
      CHARACTER(*),INTENT(OUT) :: VAL
!     **************************************************************************
      IF(ITHISISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$GETCH')
      END IF
      IF(ID.EQ.'FUNCTIONALID') THEN
        VAL=THISISP%FUNCTIONALID
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$GETCH')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$REPORT(NFIL)
!     **************************************************************************
!     **  REPORTS THE SETTINGS OF THE LDAPLUSU OBJECT                         **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      USE CONSTANTS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      TYPE(THISISP_TYPE),POINTER :: THISISP1
      CHARACTER(64)          :: STRING,NAME
      INTEGER(4)             :: L,ISP,IAT,ISP1
      REAL(8)                :: EV
!     **************************************************************************
      IF(.NOT.TON) RETURN
      CALL REPORT$TITLE(NFIL,'HYBRID FUNCTIONAL')
      DO ISP=1,NTHISISP
        THISISP1=>THISISPARRAY(ISP)
        IF(.NOT.THISISP1%TON) CYCLE
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETCH('ID',STRING)
        CALL SETUP$UNSELECT()
        THISISP1%ATOMTYPE=STRING
        CALL REPORT$CHVAL(NFIL,'ATOM TYPE',TRIM(STRING))
        CALL REPORT$CHVAL(NFIL,'  FUNCTIONAL TYPE',THISISP1%FUNCTIONALID)
        CALL REPORT$R8VAL(NFIL,'  EXTENT OF LOCAL ORBITALS',THISISP1%RCUT,'A_0')
        CALL REPORT$I4VAL(NFIL,'  MAX. ANG.MOMENTUM OF THE DENSITY',THISISP1%LRX,' ')
        DO L=0,SIZE(THISISP1%NORB)-1
          IF(THISISP1%NORB(L+1).EQ.0) CYCLE
          IF(L.EQ.0) THEN
            CALL REPORT$I4VAL(NFIL,'  NUMBER OF CORRELATED S-SHELLS',THISISP1%NORB(L+1),' ')
          ELSE IF(L.EQ.1) THEN
            CALL REPORT$I4VAL(NFIL,'  NUMBER OF CORRELATED P-SHELLS',THISISP1%NORB(L+1),' ')
          ELSE IF(L.EQ.2) THEN
            CALL REPORT$I4VAL(NFIL,'  NUMBER OF CORRELATED D-SHELLS',THISISP1%NORB(L+1),' ')
          ELSE IF(L.EQ.3) THEN
            CALL REPORT$I4VAL(NFIL,'  NUMBER OF CORRELATED F-SHELLS',THISISP1%NORB(L+1),' ')
          ELSE 
            WRITE(STRING,FMT='("  NUMBER OF CORRELATED SHELLS WITH L=",I2)')L+1
            CALL REPORT$I4VAL(NFIL,STRING,THISISP1%NORB(L+1),' ')
          END IF
        ENDDO
        IF(THISISP1%FUNCTIONALID.EQ.'HYBRID') THEN
          CALL REPORT$R8VAL(NFIL,'  HARTREE-FOCK CONTRIBUTION' &
     &                           ,THISISP1%HFWEIGHT*100,'PERCENT')
        ELSE IF(THISISP1%FUNCTIONALID.EQ.'LDA+U' &
     &      .OR.THISISP1%FUNCTIONALID.EQ.'LDA+U(OLD)') THEN
          CALL CONSTANTS('EV',EV)
          CALL REPORT$R8VAL(NFIL,'  U-PARAMETER',THISISP1%UPAR/EV,'EV')
          CALL REPORT$R8VAL(NFIL,'  J-PARAMETER',THISISP1%JPAR/EV,'EV')
          CALL REPORT$R8VAL(NFIL,'  F4/F2',THISISP1%FRATIO42,'')
          CALL REPORT$R8VAL(NFIL,'  F6/F2',THISISP1%FRATIO62,'')
          CALL REPORT$CHVAL(NFIL,'  DOUBLE COUNTING TYPE',DCTYPE)
        END IF
        DO IAT=1,NTHIS
          CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP1)
          IF(ISP1.NE.ISP) CYCLE
          CALL ATOMLIST$GETCH('NAME',IAT,NAME)
          CALL REPORT$CHVAL(NFIL,'ATOM',TRIM(NAME))
        ENDDO
      ENDDO
      NULLIFY(THISISP1)
      RETURN
      END SUBROUTINE LDAPLUSU$REPORT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$ETOT(IAT,LMNXX,NDIMD,DENMAT_,ETOT,DATH_)
!     **************************************************************************
!     **  THIS IS THE MAIN DRIVER ROUTINE FOR THE LDA+U CORRECTION            **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT     ! ATOM INDEX
      INTEGER(4),INTENT(IN)  :: LMNXX
      INTEGER(4),INTENT(IN)  :: NDIMD
      COMPLEX(8),INTENT(IN)  :: DENMAT_(LMNXX,LMNXX,NDIMD)  ! DENSITY MATRIX
      COMPLEX(8)  :: DENMAT(LMNXX,LMNXX,NDIMD)  ! DENSITY MATRIX
      REAL(8)   ,INTENT(OUT) :: ETOT
      REAL(8)   ,INTENT(OUT) :: DATH_(LMNXX,LMNXX,NDIMD)
      CHARACTER(8),PARAMETER :: CHITYPE='FROMPHI'
      INTEGER(4)             :: LRX
      COMPLEX(8),ALLOCATABLE :: DATH(:,:,:)
      INTEGER(4)             :: GID
      INTEGER(4)             :: NR
      INTEGER(4)             :: LNX
      INTEGER(4)             :: LMRX
      INTEGER(4),ALLOCATABLE :: LOX(:)
      INTEGER(4)             :: NCHI
      REAL(8)   ,ALLOCATABLE :: PHITOCHI(:,:)
      INTEGER(4)             :: LNXPHI
      INTEGER(4),ALLOCATABLE :: LOXPHI(:)
      INTEGER(4)             :: NPHI
      REAL(8)   ,ALLOCATABLE :: U(:,:,:,:)
      COMPLEX(8),ALLOCATABLE :: RHO(:,:,:,:)
      COMPLEX(8),ALLOCATABLE :: HAM(:,:,:,:)
      COMPLEX(8),ALLOCATABLE :: HAM1(:,:,:,:)
      COMPLEX(8),ALLOCATABLE :: MATSS(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: AECORE(:)
      INTEGER(4)             :: IS1,IS2,I,LN,M
      INTEGER(4)             :: ISP        ! ATOM-TYPE INDEX
      REAL(8)                :: SVAR,ETOT1
INTEGER(4)             :: LN1,LN2,LN3,LN4
      LOGICAL(4),PARAMETER   :: TCI=.FALSE.
      LOGICAL(4)             :: TCHK
CHARACTER(1) :: SWITCH
!     **************************************************************************
      ETOT=0.D0
      DATH_=0.D0
      IF(.NOT.TON) RETURN
PRINT*,'==== ENTERING LDAPLUSU$ETOT ',IAT,TON
      CALL LDAPLUSU$SELECT(IAT)
      TCHK=THIS%SP%TON
PRINT*,'==== CHECKING IF ATOM ',IAT,' IS ACTIVE ',TCHK
      IF(.NOT.TCHK) THEN
        CALL LDAPLUSU$SELECT(0)  ! UNSELECTING IS REQUIRED
        RETURN
      END IF
                            CALL TRACE$PUSH('LDAPLUSU$ETOT')
!      
!     ==========================================================================
!     ==  CONSTRUCT LOCAL ORBITALS                                            ==
!     ==========================================================================
      CALL ATOMLIST$GETI4('ISPECIES',IAT,ISP)
      CALL SETUP$ISELECT(ISP)
      IF(.NOT.THIS%SP%TINI) THEN
        CALL SETUP$GETI4('LMRX',LMRX)
        THIS%SP%LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
        LRX=THIS%SP%LRX
!
        CALL LDAPLUSU_CHIFROMPHI()   ! USES SETTING OF SETUP OBJECT
PRINT*,'TINI ',THIS%SP%TINI
PRINT*,'TON ',THIS%SP%TON
PRINT*,'LNXCHI ',THIS%SP%LNXCHI
PRINT*,'NR     ',THIS%SP%NR
PRINT*,'LOXCHI ',THIS%SP%LOXCHI
PRINT*,'NORB   ',THIS%SP%NORB
PRINT*,'NCHI   ',THIS%NCHI
PRINT*,'DIEL   ',THIS%SP%DIEL
PRINT*,'UPAR   ',THIS%SP%UPAR
PRINT*,'JPAR   ',THIS%SP%JPAR
!       ========================================================================
!       ==  CALCULATE SMALL U-TENSOR                                          ==
!       ========================================================================
        GID=THIS%SP%GID
        NR=THIS%SP%NR
        LNX=THIS%SP%LNXCHI
        ALLOCATE(LOX(LNX))
        LOX=THIS%SP%LOXCHI
        NCHI=THIS%NCHI
        ALLOCATE(THIS%ULITTLE(LRX+1,LNX,LNX,LNX,LNX))
        CALL LDAPLUSU_ULITTLE(GID,NR,LRX,LNX,LOX,THIS%CHI,THIS%ULITTLE)
        IF(THIS%SP%USEDIEL.AND.(THIS%SP%USEUPAR.OR.THIS%SP%USEJPAR)) THEN
          CALL ERROR$MSG('DIEL AND (UPAR.OR.JPAR) ARE INCOMPATIBLE')
          CALL ERROR$L4VAL('USEDIEL',THIS%SP%USEDIEL)
          CALL ERROR$L4VAL('USEUPAR',THIS%SP%USEUPAR)
          CALL ERROR$L4VAL('USEJPAR',THIS%SP%USEJPAR)
          CALL ERROR$STOP('LDAPLUSU$ETOT')
        END IF
        IF(THIS%SP%USEDIEL) THEN
          THIS%ULITTLE=THIS%ULITTLE/THIS%SP%DIEL
        ELSE IF(THIS%SP%USEUPAR.OR.THIS%SP%USEJPAR) THEN
          CALL LDAPLUSU_MODULITTLEWITHPARMS(LNX,LOX,LRX &
     &         ,THIS%SP%USEUPAR,THIS%SP%UPAR,THIS%SP%USEJPAR,THIS%SP%JPAR &
     &         ,THIS%SP%USEFRATIO42,THIS%SP%FRATIO42,THIS%SP%USEFRATIO62,THIS%SP%FRATIO62 &
     &         ,THIS%SP%MAINLN,THIS%ULITTLE)
        END IF
      ELSE
        GID=THIS%SP%GID
        NR=THIS%SP%NR
        LRX=THIS%SP%LRX
        LNX=THIS%SP%LNXCHI
        ALLOCATE(LOX(LNX))
        LOX=THIS%SP%LOXCHI
        NCHI=THIS%NCHI
      END IF
      CALL SETUP$GETI4('LNX',LNXPHI)
      ALLOCATE(LOXPHI(LNXPHI))
      CALL SETUP$GETI4A('LOX',LNXPHI,LOXPHI)
      NPHI=SUM(2*LOXPHI+1)
      CALL SETUP$UNSELECT()
!
!     ==========================================================================
!     ==  DOWNFOLD                                                            ==
!     ==========================================================================
      ALLOCATE(PHITOCHI(NCHI,NPHI))
      ALLOCATE(RHO(NCHI,NCHI,2,2))
      ALLOCATE(HAM(NCHI,NCHI,2,2))
      ALLOCATE(MATSS(NPHI,NPHI,2,2))
      ALLOCATE(DATH(NPHI,NPHI,NDIMD))
      ALLOCATE(U(NCHI,NCHI,NCHI,NCHI))
SWITCH='0'
1221 CONTINUE
DENMAT=DENMAT_
IF(SWITCH.EQ.'+') THEN
  DENMAT(1,2,2)=DENMAT(1,2,2)+1.D-3
  DENMAT(2,1,2)=DENMAT(2,1,2)+1.D-3
ELSE IF(SWITCH.EQ.'-') THEN
  DENMAT(1,2,2)=DENMAT(1,2,2)-1.D-3
  DENMAT(2,1,2)=DENMAT(2,1,2)-1.D-3
END IF
!
!     == TRANSFORM FROM TOTAL/SPIN TO UP/DOWN REPRESENTATION ===================

      CALL LDAPLUSU_SPINDENMAT('FORWARD',NDIMD,NPHI,DENMAT(1:NPHI,1:NPHI,:),MATSS)
!
PRINT*,'THIS%SP%FUNCTIONALID',THIS%SP%FUNCTIONALID
      IF(THIS%SP%FUNCTIONALID.EQ.'LDA+U(OLD)') THEN
        CALL SETUP$ISELECT(ISP)
        CALL LDAPLUSU_DENMATFLAPW('FORWARD',NPHI,MATSS,NCHI,RHO)
        CALL SETUP$UNSELECT()
      ELSE
!       == TRANSFORM FROM PARTIAL WAVES PHI TO LOCAL ORBITALS CHI ==============
        CALL LDAPLUSU_MAPTOCHI(LNX,LOX,NCHI,LNXPHI,LOXPHI,NPHI,PHITOCHI)
        DO IS1=1,2
          DO IS2=1,2
            RHO(:,:,IS1,IS2)=MATMUL(PHITOCHI &
     &                      ,MATMUL(MATSS(:,:,IS1,IS2),TRANSPOSE(PHITOCHI)))
          ENDDO
        ENDDO
      END IF
!
! PRINTOUT FOR TESTING==========================================================
DO IS1=1,NDIMD
  PRINT*,'===================== DENMAT FOR SPIN',IS1,' ======================'
  I=0
  DO LN=1,LNXPHI
    DO M=1,2*LOXPHI(LN)+1
      I=I+1
      WRITE(*,FMT='(I3,100F8.3)')LOXPHI(LN),REAL(DENMAT(I,:,IS1))
    ENDDO
  ENDDO
ENDDO
PRINT*,'===================== PHITOCHI ======================'
DO I=1,NCHI
  WRITE(*,FMT='(I3,100F8.3)')I,PHITOCHI(I,:)
ENDDO
DO IS1=1,2
  DO IS2=1,2
    IF(SUM(ABS(RHO(:,:,IS1,IS2))).LT.1.D-3) CYCLE
    PRINT*,'===================== RHO FOR SPIN',IS1,IS2,' ====================='
    I=0
    DO LN=1,LNX
      DO M=1,2*LOX(LN)+1
        I=I+1
        WRITE(*,FMT='("RE",I3,100F8.3)')LOX(LN),REAL(RHO(I,:,IS1,IS2))
        WRITE(*,FMT='("IM",I3,100F8.3)')LOX(LN),AIMAG(RHO(I,:,IS1,IS2))
      ENDDO
    ENDDO
  ENDDO
ENDDO
!
DO IS1=1,2
  SVAR=0.D0
  DO LN=1,NCHI
    SVAR=SVAR+REAL(RHO(LN,LN,IS1,IS1))
  ENDDO
  PRINT*,'CHARGE= ',SVAR,' FOR SPIN ',IS1
ENDDO
!
!!$DO LN1=1,LNX
!!$  DO LN2=1,LNX
!!$    DO LN3=1,LNX
!!$      DO LN4=1,LNX
!!$        WRITE(*,*)'ULITTLE',LN1,LN2,LN3,LN4,THIS%ULITTLE(:,LN1,LN2,LN3,LN4)
!!$      ENDDO
!!$    ENDDO
!!$  ENDDO
!!$ENDDO
!
!     ==========================================================================
!     ==  CALCULATE U-TENSOR                                                  ==
!     ==========================================================================
      CALL LDAPLUSU_UTENSOR(LRX,NCHI,LNX,LOX,THIS%ULITTLE,U)
!
!     ==========================================================================
!     ==  HARTREE FOCK INTERACTION ENERGY                                     ==
!     ==========================================================================
      IF(TCI) THEN
        CALL LDAPLUSU_CI(NCHI,U,RHO,ETOT,HAM)
      ELSE
        CALL LDAPLUSU_INTERACTION(NCHI,U,RHO,ETOT,HAM)
      END IF
!
!     ==========================================================================
!     ==  SUBTRACT HARTREE ENERGY                                             ==
!     ==========================================================================
      ALLOCATE(HAM1(NCHI,NCHI,2,2))
      CALL LDAPLUSU_EHARTREE(GID,NR,NCHI,LNX,LOX,THIS%CHI,LRX,RHO,ETOT1,HAM1)
      ETOT=ETOT-ETOT1
      HAM=HAM-HAM1
      DEALLOCATE(HAM1)
!
!     ==========================================================================
!     ==  CORE VALENCE EXCHANGE INTERACTION                                   ==
!     ==========================================================================
      IF(THIS%SP%TCV) THEN
        ALLOCATE(HAM1(NCHI,NCHI,2,2))
        CALL LDAPLUSU_CVX(NCHI,LNX,LOX,RHO,THIS%CVX,ETOT1,HAM1)
        ETOT=ETOT+ETOT1
        HAM=HAM+HAM1
        DEALLOCATE(HAM1)
      END IF
!
!     ==========================================================================
!     ==  DOUBLE COUNTING CORRECTION                                          ==
!     ==========================================================================
      IF(THIS%SP%FUNCTIONALID.EQ.'LDA+U' &
     &              .OR.THIS%SP%FUNCTIONALID.EQ.'LDA+U(OLD)') THEN
        ALLOCATE(HAM1(NCHI,NCHI,2,2))
        CALL LDAPLUSU_EHARTREE(GID,NR,NCHI,LNX,LOX,THIS%CHI,LRX,RHO,ETOT1,HAM1)
        ETOT=ETOT+ETOT1
        HAM=HAM+HAM1
!
        CALL LDAPLUSU_DCLDAPLUSU(DCTYPE,LNX,LOX,NCHI,U,RHO,ETOT1,HAM1)
        ETOT=ETOT-ETOT1
        HAM=HAM-HAM1
        DEALLOCATE(HAM1)
!
      ELSE IF(THIS%SP%FUNCTIONALID.EQ.'HYBRID') THEN
        CALL SETUP$ISELECT(ISP)
        ALLOCATE(AECORE(NR))
        IF(THIS%SP%TCV) THEN
          CALL SETUP$GETR8A('AECORE',NR,AECORE)
        ELSE
          AECORE(:)=0.D0
        END IF
        CALL SETUP$UNSELECT()
        ALLOCATE(HAM1(NCHI,NCHI,2,2))
        CALL LDAPLUSU_EDFT(GID,NR,NCHI,LNX,LOX,THIS%CHI,LRX,AECORE,RHO &
    &                     ,ETOT1,HAM1)

        ETOT=ETOT-ETOT1
        HAM=HAM-HAM1
        DEALLOCATE(HAM1)
        DEALLOCATE(AECORE)
!       == SCALE CORRECTION WITH 0.25 ACCORDING TO PBE0
        ETOT=ETOT*THIS%SP%HFWEIGHT
        HAM=HAM*THIS%SP%HFWEIGHT
      ELSE
        CALL ERROR$MSG('FUNCTIONALID NOT RECOGNIZED')
        CALL ERROR$CHVAL('FUNCTIONALID',THIS%SP%FUNCTIONALID)
        CALL ERROR$STOP('LDAPLUSU$ETOT')
      END IF
PRINT*,'ETOT AFTER SCALING ',ETOT
DO I=1,NCHI
  WRITE(*,FMT='(I3,100F12.5)')I,REAL(HAM(I,:,1,1))
ENDDO
!
!     ==========================================================================
!     ==  UPFOLD                                                              ==
!     ==========================================================================
      IF(THIS%SP%FUNCTIONALID.EQ.'LDA+U(OLD)') THEN
        CALL LDAPLUSU_DENMATFLAPW('BACK',NPHI,MATSS,NCHI,HAM)
      ELSE
!       == TRANSFORM FROM CHI TO PHI ===========================================
        DO IS1=1,2
          DO IS2=1,2
            MATSS(:,:,IS1,IS2)=MATMUL(TRANSPOSE(PHITOCHI),MATMUL(HAM(:,:,IS1,IS2),PHITOCHI))
          ENDDO
        ENDDO
      END IF
!
!     == TRANSFORM FROM (SPIN,SPIN) TO (TOTAL,SPIN) ============================
      CALL LDAPLUSU_SPINDENMAT('BACK',NDIMD,LMNXX,DATH,MATSS)
!
!     == MAKE REAL (THIS IS A FUDGE TO BE FIXED IN AUGMENTATION!)
      DATH_(:,:,:)=(0.D0,0.D0)
      DATH_(:NPHI,:NPHI,:)=REAL(DATH)

PRINT*,'============ TOTAL LDAPLUSU ====IAT=',IAT,'============================'
PRINT*,'IAT=',IAT,' LDA+U ETOT ',ETOT,' SWITCH=',SWITCH
DO IS1=1,NDIMD
PRINT*,'===================== DENMAT FOR SPIN',IS1,' ======================'
I=0
DO LN=1,LNXPHI
DO M=1,2*LOXPHI(LN)+1
I=I+1
WRITE(*,FMT='(I3,100F12.5)')LOXPHI(LN),REAL(DENMAT(I,:,IS1))
ENDDO
ENDDO
ENDDO
IF(SWITCH.EQ.'0') THEN
DO IS1=1,NDIMD
PRINT*,'===================== DATH FOR SPIN',IS1,' ========================'
I=0
DO LN=1,LNXPHI
DO M=1,2*LOXPHI(LN)+1
I=I+1
WRITE(*,FMT='(I3,100F12.5)')LOXPHI(LN),REAL(DATH_(I,:,IS1))
ENDDO
ENDDO
ENDDO
END IF

IF(SWITCH.EQ.'0') THEN
 SWITCH='+'
ELSE IF(SWITCH.EQ.'+') THEN
 SWITCH='-'
ELSE
 STOP 'FORCED'
END IF
!GOTO 1221
!
!     ==========================================================================
!     ==  UNSELECT LDAPLUSU                                                   ==
!     ==========================================================================
      DEALLOCATE(U)
      DEALLOCATE(LOXPHI)
      DEALLOCATE(LOX)
      DEALLOCATE(MATSS)
      DEALLOCATE(RHO)
      DEALLOCATE(HAM)
      CALL LDAPLUSU$SELECT(0)
!CALL ERROR$MSG('DENMAT IS CHANGED INTERNALLY FOR TESTING PURPOSES!!!!')
!CALL ERROR$STOP('FORCED IN LDAPLUS')
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_CHIFROMPHI()
!     **************************************************************************
!     **                                                                      **
!     **  DEFINES TRANSFORMATION FROM PARTIAL WAVES TO LOCAL ORBITALS.        **
!     **  THE RESULT IS STORED WITHOUT MAGNETIC QUANTUM NUMBERS.              **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: LNX
      INTEGER(4)            :: ISP
      INTEGER(4),ALLOCATABLE:: LOX(:)
      INTEGER(4),ALLOCATABLE:: ISCATT(:)
      REAL(8)   ,ALLOCATABLE:: PHI(:,:)
      INTEGER(4)            :: LNXCHI 
      INTEGER(4),ALLOCATABLE:: LOXCHI(:)
      REAL(8)   ,ALLOCATABLE:: MAT(:,:)
      REAL(8)   ,ALLOCATABLE:: R(:)        
      REAL(8)   ,ALLOCATABLE:: AUX(:)        
      REAL(8)               :: SVAR1,SVAR2
      INTEGER(4)            :: IR
      INTEGER(4)            :: IAT
      INTEGER(4)            :: N,LX,L,LN,LNCHI,LN0,NOFL,ISVAR
      INTEGER(4)            :: N1,N2,LN1,LN2,L1,L2
      REAL(8)   ,ALLOCATABLE:: AMAT(:,:),BMAT(:,:)
      LOGICAL(4),ALLOCATABLE:: TORB(:)
!     **************************************************************************
                            CALL TRACE$PUSH('LDAPLUSU_CHIFROMPHI')
!     == SETUP IS STILL SELECTED BY PARENT ROUTINE =============================
      CALL SETUP$GETI4('GID',GID)
      THIS%SP%GID=GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      THIS%SP%NR=NR
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      CALL SETUP$GETI4('LNX',LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$GETI4A('LOX',LNX,LOX)
      ALLOCATE(PHI(NR,LNX))
      CALL SETUP$GETR8A('AEPHI',NR*LNX,PHI)

!     == COLLECT THE SELECTOR OF LOCAL ORBITALS ================================
      ALLOCATE(TORB(LNX))
      IF(SIZE(THIS%SP%NORB).LT.MAXVAL(LOX)+1) THEN
        CALL ERROR$MSG('CONFLICT OF ARRAY DIMENSIONS: NORB < LX')
        CALL ERROR$I4VAL('SIZE(NORB)',SIZE(THIS%SP%NORB))
        CALL ERROR$I4VAL('MAX(LOX)+1',MAXVAL(LOX)+1)
        CALL ERROR$CHVAL('ATOMTYPE',THIS%SP%ATOMTYPE)
        CALL ERROR$STOP('LDAPLUSU_CHIFROMPHI')
      END IF
      DO L=0,MAXVAL(LOX)
        ISVAR=THIS%SP%NORB(L+1)
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE
          TORB(LN)=ISVAR.GT.0
          ISVAR=ISVAR-1
        ENDDO
      ENDDO
!
!     == CHECK CONSISTENCY WITH ISCATT =========================================
      ALLOCATE(ISCATT(LNX))
      CALL SETUP$GETI4A('ISCATT',LNX,ISCATT)
      DO LN=1,LNX
        IF(TORB(LN).AND.ISCATT(LN).GT.0) THEN
          CALL ERROR$MSG('SCATTERING STATES MUST NOT AMONG CORRELATED ORBITALS')
          CALL ERROR$STOP('LDAPLUSU_CHIFROMPHI')
        END IF
      ENDDO
!
      LNXCHI=0
      DO LN=1,LNX
        IF(TORB(LN))LNXCHI=LNXCHI+1
      ENDDO
      THIS%SP%LNXCHI=LNXCHI
      ALLOCATE(LOXCHI(LNXCHI))
      LNCHI=0
      DO LN=1,LNX
        IF(.NOT.TORB(LN))CYCLE
        LNCHI=LNCHI+1
        LOXCHI(LNCHI)=LOX(LN)
      ENDDO
      ALLOCATE(THIS%SP%LOXCHI(LNXCHI))
      THIS%SP%LOXCHI=LOXCHI
      THIS%NCHI=SUM(2*LOXCHI(1:LNXCHI)+1)
      DEALLOCATE(LOXCHI)
!      ALLOCATE(AMAT(LNX,LNXCHI))
!      ALLOCATE(BMAT(LNXCHI,LNX))
!      IAT=ITHIS
!      CALL LMTO$DOLOCORB_OLD(IAT,LNXCHI,LNX,TORB,AMAT,BMAT)
!
!     == STORE DOWNFOLD MATRIX =================================================
!      IF(ASSOCIATED(THIS%DOWNFOLD))DEALLOCATE(THIS%DOWNFOLD)
!      ALLOCATE(THIS%DOWNFOLD(LNXCHI,LNX))
!      DO LN=1,LNXCHI
!        THIS%DOWNFOLD(LN,:)=BMAT(LN,:)
!      ENDDO
!      DEALLOCATE(BMAT)
!
!     == CONSTRUCT LOCAL ORBITALS ==============================================
!      IF(ASSOCIATED(THIS%CHI))DEALLOCATE(THIS%CHI)
!      ALLOCATE(THIS%CHI(NR,LNXCHI))
!      THIS%CHI=M≤ATMUL(PHI,AMAT) 
!      DEALLOCATE(AMAT)

      IAT=ITHIS
      IF(ASSOCIATED(THIS%DOWNFOLD))DEALLOCATE(THIS%DOWNFOLD)
      ALLOCATE(THIS%DOWNFOLD(LNXCHI,LNX))
      ALLOCATE(THIS%CHI(NR,LNXCHI))
      CALL SETUP$GETI4('ISP',ISP)
PRINT*,'IAT  ',IAT,' LNXCHI=',LNXCHI,' LNX=',LNX,' ISP=',ISP
PRINT*,'TORB ',TORB
!      CALL LMTO$DOLOCORB(IAT,ISP,GID,NR,LNXCHI,LNX,TORB,THIS%DOWNFOLD,THIS%CHI)
      CALL LMTO$DOLOCORB_2(IAT,ISP,GID,NR,LNXCHI,LNX,TORB,THIS%DOWNFOLD,THIS%CHI)
!
!     == DETERMINE CORE-VALENCE EXCHANGE =======================================
      IF(THIS%SP%TCV) THEN
        IF(ASSOCIATED(THIS%CVX))DEALLOCATE(THIS%CVX)
        ALLOCATE(THIS%CVX(LNXCHI,LNXCHI))
!        ALLOCATE(MAT(LNX,LNX))
!        CALL SETUP$GETR8A('CVX',LNX**2,MAT)
!        THIS%CVX(:,:)=MATMUL(TRANSPOSE(AMAT),MATMUL(MAT,AMAT))
!        DEALLOCATE(MAT)
        CALL LDAPLUSU_CVXSETUP(ISP,GID,NR,LNXCHI,THIS%SP%LOXCHI,THIS%CHI,THIS%CVX)
      ELSE
        NULLIFY(THIS%CVX)
      END IF
      DEALLOCATE(PHI)
!DO LN=1,LNXCHI
!WRITE(*,FMT='("A  ",I5,10F10.5)')LN,THIS%DOWNFOLD(LN,:)
!ENDDO
!!$
!!$IF(IAT.EQ.1)CALL ATOMLIB_WRITEPHI('CHI_1.DAT',GID,NR,LNXCHI,THIS%CHI)
!!$IF(IAT.EQ.2)CALL ATOMLIB_WRITEPHI('CHI_2.DAT',GID,NR,LNXCHI,THIS%CHI)
!!$ALLOCATE(AUX(NR))
!!$DO LN=1,LNXCHI
!!$  AUX(:)=R(:)**2*THIS%CHI(:,LN)**2
!!$  CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$  PRINT*,"NORM ISP=",ISP," LN= ",LN," NORM=",SVAR1
!!$ENDDO
!!$!
      DEALLOCATE(R)
                            CALL TRACE$POP()
      RETURN      
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_CVXSETUP(ISP,GID,NR,LNX,LOX,AECHI,MAT)
!     **************************************************************************
!     **  CORE-VALENCE EXCHANGE MATRIX ELEMENTS                               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      REAL(8)   ,INTENT(IN) :: AECHI(NR,LNX)
      REAL(8)   ,INTENT(OUT):: MAT(LNX,LNX)
      INTEGER(4)            :: NB,NC
      INTEGER(4)            :: LMRX
      INTEGER(4)            :: LRHOX
      INTEGER(4),ALLOCATABLE:: LOFI(:)
      REAL(8)   ,ALLOCATABLE:: AEPSI(:,:)
!     **************************************************************************
      MAT(:,:)=0.D0
      CALL SETUP$GETI4('NB',NB)
      CALL SETUP$GETI4('NC',NC)
      CALL SETUP$GETI4('LMRX',LMRX)
      LRHOX=INT(SQRT(REAL(LMRX-1)+1.D-8))
      ALLOCATE(LOFI(NB))
      ALLOCATE(AEPSI(NR,NB))
      CALL SETUP$GETI4A('LB',NB,LOFI)
      CALL SETUP$GETR8A('AEPSI',NR*NB,AEPSI)
      CALL LDAPLUSU_CVXMAT(GID,NR,LNX,LOX,AECHI,NC,LOFI(:NC),AEPSI(:,:NC),LRHOX,MAT)
      DEALLOCATE(AEPSI)
      DEALLOCATE(LOFI)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_CVXMAT(GID,NR,LNX,LOX,AEPHI,NC,LOFC,PSIC,LRHOX,MAT)
!     **************************************************************************
!     **  CORE-VALENCE EXCHANGE MATRIX ELEMENTS                               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID           ! GRID ID
      INTEGER(4),INTENT(IN) :: NR            ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: LNX           ! #(PARTIAL WAVES W/O M,SIGMA)
      INTEGER(4),INTENT(IN) :: LOX(LNX)      ! ANGULAR MOMENTA OF PARTIAL WAVES
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX) ! AE PARTIAL WAVES
      INTEGER(4),INTENT(IN) :: NC            ! #(CORE STATES)
      INTEGER(4),INTENT(IN) :: LOFC(NC)      ! ANGULAR MOMENTA OF CORE STATES
      REAL(8)   ,INTENT(IN) :: PSIC(NR,NC)   ! CORE STATES
      INTEGER(4),INTENT(IN) :: LRHOX         ! MAX ANGULAR MOMENTUM OF DENSITY
      REAL(8)   ,INTENT(OUT):: MAT(LNX,LNX)  ! CORE VALENCE X MATRIX ELEMENTS
      INTEGER(4)            :: LX  ! MAX ANGULAR MOMENTUM PARTIAL WAVES
      INTEGER(4)            :: LMX ! MAX #(ANGULAR MOMENTA) OF PARTIAL WAVES
      INTEGER(4)            :: LCX ! HIGHEST ANGULAR MOMENTUM OF CORE STATES
      REAL(8)               :: CG
      REAL(8)               :: RHO1(NR)
      REAL(8)               :: AUX(NR),POT(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VAL
      REAL(8)   ,ALLOCATABLE:: FACTOR(:,:,:)
      INTEGER(4)            :: LMCA
      INTEGER(4)            :: LM1,LC,LRHO,IMC,LMC,LMRHOA,IMRHO,LMRHO
      INTEGER(4)            :: LN1,L1,IC,LN2,L2,LM2
      LOGICAL(4),PARAMETER  :: TPRINT=.TRUE.
      REAL(8)               :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      LX=MAXVAL(LOX)
      LMX=(LX+1)**2
      LCX=MAXVAL(LOFC)
      ALLOCATE(FACTOR(LX+1,LRHOX+1,LCX+1))
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == INCLUDE ANGULAR INTEGRATIONS                                         ==
!     ==========================================================================
!     == SPHERICAL SYMMETRY EXPLOITED: 
      FACTOR=0.D0
      DO L1=0,LX           ! ANGULAR MOMENTUM OF LOCAL ORBITAL
        LM1=L1**2+1
        DO LC=0,LCX        ! ANGULAR MOMENTUM OF CORE STATE
          DO LRHO=0,LRHOX  ! ANGULAR MOMENTUM OF DENSITY
            LMCA=LC**2
            DO IMC=1,2*LC+1
              LMC=LMCA+IMC
              LMRHOA=LRHO**2
              DO IMRHO=1,2*LRHO+1
                LMRHO=LMRHOA+IMRHO
                CALL SPHERICAL$GAUNT(LM1,LMC,LMRHO,CG)
                FACTOR(L1+1,LRHO+1,LC+1)=FACTOR(L1+1,LRHO+1,LC+1)+CG**2
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     = THIS FACTOR IS INCLUDED TO HAVE A PROPER DEFINITION OF SLATER INTEGRALS
      DO LRHO=0,LRHOX
        FACTOR(:,LRHO+1,:)=FACTOR(:,LRHO+1,:)*4.D0*PI/REAL(2*LRHO+1,KIND=8)
      ENDDO
!
!     ==========================================================================
!     == WORK OUT RADIAL INTEGRATIONS                                         ==
!     ==========================================================================
      MAT(:,:)=0.D0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IC=1,NC
          LC=LOFC(IC)
          RHO1(:)=AEPHI(:,LN1)*PSIC(:,IC)
          DO LRHO=0,LRHOX
            CALL RADIAL$POISSON(GID,NR,LRHO,RHO1,POT)
            DO LN2=1,LN1  ! MATRIX IS SYMMETRIC
              L2=LOX(LN2)
              IF(L2.NE.L1) CYCLE
              AUX(:)=R(:)**2*POT(:)*PSIC(:,IC)*AEPHI(:,LN2)
              CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!!$PRINT*,LN1,IC,LC,LRHO,LN2,'VAL ',VAL
!!$CALL ATOMLIB_WRITEPHI('X0.DAT',GID,NR,1,AUX)
              VAL=VAL*REAL(2*LRHO+1,KIND=8)/(4.D0*PI)  !SLATER INTEGRAL
              MAT(LN1,LN2)=MAT(LN1,LN2)-FACTOR(L1+1,LRHO+1,LC+1)*VAL
              MAT(LN2,LN1)=MAT(LN1,LN2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(FACTOR)
!
!     ==========================================================================
!     ==  WRITE FOR TEST PURPOSES                                             ==
!     ==========================================================================
      IF(TPRINT) THEN
        PRINT*,'LOFC ',LOFC
        PRINT*,'NOW THE MATRIX FOR THE CORE VALENCE EXCHANGE INTERACTION'
        WRITE(*,FMT='(4A3,10A18)')'LN1','L1','LN2','L2','MAT(LN1,LN2)'
        DO LN1=1,LNX
          L1=LOX(LN1)
          DO LN2=1,LNX
            L2=LOX(LN2)
            WRITE(*,FMT='(4I3,10F18.10)')LN1,L1,LN2,L2,MAT(LN1,LN2)
          ENDDO
        ENDDO
!        STOP
      END IF
!!$
!!$      DO LN1=1,LNX
!!$        L1=LOX(LN1)
!!$        LM1=L1**2+1
!!$        DO LN2=1,LNX
!!$          L2=LOX(LN2)
!!$          LM2=L2**2+1   
!!$          MAT(LN1,LN2)=0.D0 
!!$          DO IC=1,NC
!!$            LC=LOFC(IC)
!!$            DO IMC=1,2*LC+1
!!$              LMC=LC**2+IMC
!!$              CALL LDAPLUSU_SINGLEU(GID,NR,LRHOX,LM1,AEPHI(:,LN1),LMC,PSIC(:,IC) &
!!$    &                   ,LMC,PSIC(:,IC),LM2,AEPHI(:,LN2),VAL)
!!$              MAT(LN1,LN2)=MAT(LN1,LN2)+VAL
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$
!!$      IF(TPRINT) THEN
!!$        PRINT*,'LOFC ',LOFC
!!$        PRINT*,'NOW THE MATRIX FOR THE CORE VALENCE EXCHANGE INTERACTION'
!!$        WRITE(*,FMT='(4A3,10A18)')'LN1','L1','LN2','L2','MAT(LN1,LN2)'
!!$        DO LN1=1,LNX
!!$          L1=LOX(LN1)
!!$          DO LN2=1,LNX
!!$            L2=LOX(LN2)
!!$            WRITE(*,FMT='(4I3,10F18.10)')LN1,L1,LN2,L2,MAT(LN1,LN2)
!!$          ENDDO
!!$        ENDDO
!!$        STOP
!!$      END IF

      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_SPINDENMAT(ID,NDIMD,LMNXX,MAT1,MAT2)
!     **                                                                      **
!     ** IF="FORWARD": CONVERTS DENSITY MATRIX MAT1 FROM                      **
!     **   (TOTAL,SPIN) REPRESENTATION INTO (SPIN,SPIN) REPRESENTATION MAT2   **
!     ** IF="BACK": CONVERTS HAMILTON MATRIX MAT2 FROM                        **
!     **   (TOTAL,SPIN) REPRESENTATION INTO (SPIN,SPIN) REPRESENTATION MAT2   **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4),INTENT(IN)   :: NDIMD
      INTEGER(4),INTENT(IN)   :: LMNXX
      COMPLEX(8),INTENT(INOUT):: MAT1(LMNXX,LMNXX,NDIMD)
      COMPLEX(8),INTENT(INOUT):: MAT2(LMNXX,LMNXX,2,2)
      COMPLEX(8),PARAMETER    :: CI=(0.D0,1.D0)
!     **************************************************************************
                            CALL TRACE$PUSH('LDAPLUSU_SPINDENMAT')
      IF(ID.EQ.'FORWARD') THEN
        MAT2(:,:,:,:)=(0.D0,0.D0)
        IF(NDIMD.EQ.1) THEN
          MAT2(:,:,1,1)=0.5D0*MAT1(:,:,1)
          MAT2(:,:,2,2)=0.5D0*MAT1(:,:,1)
        ELSE IF(NDIMD.EQ.2) THEN
          MAT2(:,:,1,1)=0.5D0*(MAT1(:,:,1)+MAT1(:,:,2))
          MAT2(:,:,2,2)=0.5D0*(MAT1(:,:,1)-MAT1(:,:,2))
        ELSE IF(NDIMD.EQ.4) THEN
          MAT2(:,:,1,1)=0.5D0*(MAT1(:,:,1)+MAT1(:,:,4))
          MAT2(:,:,2,2)=0.5D0*(MAT1(:,:,1)-MAT1(:,:,4))
          MAT2(:,:,1,2)=0.5D0*(MAT1(:,:,2)-CI*MAT1(:,:,3))
          MAT2(:,:,2,1)=0.5D0*(MAT1(:,:,2)+CI*MAT1(:,:,3))
        ELSE
          CALL ERROR$MSG('NDIMD OUT OF RANGE')
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$STOP('LDAPLUSU_SPINDENMAT')
        END IF
      ELSE IF(ID.EQ.'BACK') THEN
        MAT1(:,:,:)=(0.D0,0.D0)
        IF(NDIMD.EQ.1) THEN
          MAT1(:,:,1)=0.5D0*(MAT2(:,:,1,1)+MAT2(:,:,2,2))
        ELSE IF(NDIMD.EQ.2) THEN
          MAT1(:,:,1)=0.5D0*(MAT2(:,:,1,1)+MAT2(:,:,2,2))
          MAT1(:,:,2)=0.5D0*(MAT2(:,:,1,1)-MAT2(:,:,2,2))
        ELSE IF(NDIMD.EQ.4) THEN
          MAT1(:,:,1)=0.5D0*(MAT2(:,:,1,1)+MAT2(:,:,2,2))
          MAT1(:,:,2)=0.5D0*(MAT2(:,:,1,2)+MAT2(:,:,2,1))
          MAT1(:,:,3)=-0.5D0*CI*(MAT2(:,:,1,2)-MAT2(:,:,2,1))
          MAT1(:,:,4)=0.5D0*(MAT2(:,:,1,1)-MAT2(:,:,2,2))
        ELSE
          CALL ERROR$MSG('NDIMD OUT OF RANGE')
          CALL ERROR$I4VAL('NDIMD',NDIMD)
          CALL ERROR$STOP('LDAPLUSU_SPINDENMAT')
        END IF
      ELSE
        CALL ERROR$MSG('ID MUST BE EITHER "FORWARD" OR "BACK"')
        CALL ERROR$STOP('LDAPLUSU_SPINDENMAT')
      END IF
                            CALL TRACE$POP()
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_MAPTOCHI(LNXCHI,LOXCHI,LMNXCHI &
     &                            ,LNXPHI,LOXPHI,LMNXPHI,DOWNFOLD)
!     **                                                                      **
!     **  EXPANDS TRANSFORMATION FROM PARTIAL WAVES TO LOCAL ORBITALS         **
!     **  TO FULL SIZE                                                        **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: LNXCHI
      INTEGER(4),INTENT(IN)   :: LOXCHI(LNXCHI)
      INTEGER(4),INTENT(IN)   :: LMNXCHI
      INTEGER(4),INTENT(IN)   :: LNXPHI
      INTEGER(4),INTENT(IN)   :: LOXPHI(LNXPHI)
      INTEGER(4),INTENT(IN)   :: LMNXPHI
      REAL(8)   ,INTENT(OUT)  :: DOWNFOLD(LMNXCHI,LMNXPHI)
      INTEGER(4)              :: LMN1,LN1,L1,LMN2,LN2,L2,M
!     **************************************************************************
                            CALL TRACE$PUSH('LDAPLUSU_MAPTOCHI')
      DOWNFOLD=0.D0
      LMN1=0
      DO LN1=1,LNXCHI
        L1=LOXCHI(LN1)
        LMN2=0
        DO LN2=1,LNXPHI
          L2=LOXPHI(LN2)
          IF(L1.EQ.L2) THEN
            DO M=1,2*L1+1
              DOWNFOLD(LMN1+M,LMN2+M)=THIS%DOWNFOLD(LN1,LN2)
            ENDDO
          END IF
          LMN2=LMN2+2*L2+1
        ENDDO
        LMN1=LMN1+2*L1+1
      ENDDO
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_ULITTLE(GID,NR,LRX,LNX,LOX,CHI,ULITTLE)
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
      SUBROUTINE LDAPLUSU_MODULITTLEWITHPARMS(LNX,LOX,LRX,USEUPAR,UPAR &
     &   ,USEJPAR,JPAR,USEFRATIO42,FRATIO42,USEFRATIO62,FRATIO62,MAINLN,ULITTLE)
!     **                                                                      **
!     ** CALCULATES THE INTERACTION ENERGY                                    **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LRX
      LOGICAL(4),INTENT(IN) :: USEUPAR
      LOGICAL(4),INTENT(IN) :: USEJPAR
      LOGICAL(4),INTENT(IN) :: USEFRATIO42
      LOGICAL(4),INTENT(IN) :: USEFRATIO62
      INTEGER(4),INTENT(IN) :: MAINLN(2)
      REAL(8)   ,INTENT(INOUT):: UPAR
      REAL(8)   ,INTENT(INOUT):: JPAR
      REAL(8)   ,INTENT(INOUT):: FRATIO42
      REAL(8)   ,INTENT(INOUT):: FRATIO62
      REAL(8)   ,INTENT(INOUT):: ULITTLE(LRX+1,LNX,LNX,LNX,LNX)
      REAL(8)   ,PARAMETER  :: FIVEEIGTH=0.625D0
      REAL(8)               :: PI,FOURPI
      INTEGER(4)            :: L,LN,LNPROBE,N
      REAL(8)               :: RAWJPAR,RAWUPAR
      REAL(8)               :: SVAR
      REAL(8)               :: XL(LRX+1)
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
!
!     ==========================================================================
!     == DETERMINE SHELL TO WHICH UPAR AND JPAR BELONG                        ==
!     ==========================================================================
      LNPROBE=-1
      N=0
      DO LN=1,LNX
        L=LOX(LNX)
        IF(L.NE.MAINLN(1)) CYCLE
        N=N+1
        IF(N.EQ.MAINLN(2))THEN
          LNPROBE=LN  
          EXIT
        END IF
      ENDDO
      IF(LNX.EQ.1) LNPROBE=1
      IF(LNPROBE.EQ.-1) THEN
        CALL ERROR$MSG('LDAPLUSU_MODULITTLEWITHPARMS')
      END IF
!
!     ==========================================================================
!     == SCALE UP ULITTLE TO SATISFY UPAR                                     ==
!     ==========================================================================
      IF(USEUPAR) THEN
        RAWUPAR=ULITTLE(1,LNPROBE,LNPROBE,LNPROBE,LNPROBE)
        SVAR=UPAR/RAWUPAR
!       == ALL SHELLS OF ORBITALS ARE SCALED SO THAT THE SPECIFIED SHELL HAS THE
!       == SPECIFIED VALUE OF U
        ULITTLE=ULITTLE*SVAR
      ELSE
        UPAR=ULITTLE(1,LNPROBE,LNPROBE,LNPROBE,LNPROBE)
      END IF
!
!     ==========================================================================
!     == CALCULATE J-PARAMETER                                                ==
!     ==========================================================================
      CALL LDAPLUSU_JEXPANSION(MAINLN(1),LRX,XL)
      IF(.NOT.USEJPAR) THEN
!       == CALCULATE J-PARAMETER
        JPAR=0.D0
        DO L=1,LRX   !EXCLUDE L=0
          JPAR=JPAR+ULITTLE(L+1,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*XL(L+1)
        ENDDO
      END IF
!
!     ==========================================================================
!     == FIX UP RATIOS OF SLATER INTEGRALS  F4/F2 AND F6/F2                   ==
!     ==========================================================================
      IF(MAINLN(1).GE.2.AND.LRX.GT.4) THEN
        IF(USEFRATIO42) THEN
          ULITTLE(5,LNPROBE,LNPROBE,LNPROBE,LNPROBE)=FRATIO42 &
    &             *ULITTLE(3,LNPROBE,LNPROBE,LNPROBE,LNPROBE)
        ELSE
          FRATIO42=ULITTLE(5,LNPROBE,LNPROBE,LNPROBE,LNPROBE) &
    &             /ULITTLE(3,LNPROBE,LNPROBE,LNPROBE,LNPROBE)
        END IF
      ELSE
        FRATIO42=0.D0
      END IF
      IF(MAINLN(1).GE.3.AND.LRX.GT.6) THEN
        IF(USEFRATIO62) THEN
          ULITTLE(7,LNPROBE,LNPROBE,LNPROBE,LNPROBE)=FRATIO62 &
    &             *ULITTLE(3,LNPROBE,LNPROBE,LNPROBE,LNPROBE)
        ELSE
          FRATIO62=ULITTLE(7,LNPROBE,LNPROBE,LNPROBE,LNPROBE) &
    &             /ULITTLE(3,LNPROBE,LNPROBE,LNPROBE,LNPROBE)
        END IF
      ELSE
        FRATIO62=0.D0
      END IF
!
!     ==========================================================================
!     == FIX JPAR                                                             ==
!     ==========================================================================
      RAWJPAR=0.D0
      DO L=1,LRX   !EXCLUDE L=0
        RAWJPAR=RAWJPAR+ULITTLE(L+1,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*XL(L+1)
      ENDDO
      IF(RAWJPAR.NE.0.D0) THEN
        SVAR=JPAR/RAWJPAR
        ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)=SVAR &
    &                 *ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)
      END IF
!
!     ==========================================================================
!     == SET UP ULITTLE                                                       ==
!     ==========================================================================
!!$
!!$      ULITTLE=0.D0
!!$      DO L=0,LRX
!!$        FAC=0.D0
!!$        IF(L.EQ.0)  FAC=UPAR
!!$        IF(LRX.GE.2)FAC=JPAR*14.D0/(1.D0+FIVEEIGTH)
!!$        IF(LRX.GE.4)FAC=FIVEEIGTH*JPAR*14.D0/(1.D0+FIVEEIGTH)
!!$        ULITTLE(L+1,:,:,:,:)=FAC*4.D0*PI/REAL(2*L+1,KIND=8)
!!$      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_MODULITTLEWITHPARMS_O(LNX,LOX,LRX,USEUPAR,UPAR &
     &                                       ,USEJPAR,JPAR,MAINLN,ULITTLE)
!     **                                                                      **
!     ** BUGGY VERSION SHALL BE REMOVED !                                     **
!     ** CALCULATES THE INTERACTION ENERGY                                    **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LRX
      LOGICAL(4),INTENT(IN) :: USEUPAR
      LOGICAL(4),INTENT(IN) :: USEJPAR
      INTEGER(4),INTENT(IN) :: MAINLN(2)
      REAL(8)   ,INTENT(IN) :: UPAR
      REAL(8)   ,INTENT(IN) :: JPAR
      REAL(8)   ,INTENT(INOUT):: ULITTLE(LRX+1,LNX,LNX,LNX,LNX)
      REAL(8)   ,PARAMETER  :: FIVEEIGTH=0.625D0
      REAL(8)               :: PI,FOURPI
      INTEGER(4)            :: L,LN,LNPROBE,N
      REAL(8)               :: RAWJPAR,RAWUPAR
      REAL(8)               :: SVAR
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
!
!     ==========================================================================
!     == DETERMIN SHELL TO WHICH UPAR AND JPAR BELONG                         ==
!     ==========================================================================
      LNPROBE=-1
      N=0
      DO LN=1,LNX
        L=LOX(LNX)
        IF(L.NE.MAINLN(1)) CYCLE
        N=N+1
        IF(N.EQ.MAINLN(2))THEN
          LNPROBE=LN  
          EXIT
        END IF
      ENDDO
      IF(LNX.EQ.1) LNPROBE=1
      IF(LNPROBE.EQ.-1) THEN
        CALL ERROR$MSG('LDAPLUSU_MODULITTLEWITHPARMS')
      END IF
!
!     ==========================================================================
!     == SCALE UP ULITTLE TO SATISFY UPAR                                     ==
!     ==========================================================================
      IF(USEUPAR) THEN
        RAWUPAR=ULITTLE(1,LNPROBE,LNPROBE,LNPROBE,LNPROBE)
        SVAR=UPAR/RAWUPAR
        ULITTLE=ULITTLE*SVAR
      END IF
!
!     ==========================================================================
!     == SCALE UP JPAR OF THE MAIN SHELL AND SET OTHERS TO ZERO               ==
!     ==========================================================================
      IF(USEJPAR) THEN
        IF(LOX(LNPROBE).EQ.2) THEN
          RAWJPAR=0.D0
          IF(LRX+1.GT.3) THEN
            SVAR=1.D0/(14.D0*FOURPI)
            RAWJPAR=RAWJPAR+ULITTLE(3,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
          END IF 
          IF(LRX+1.GT.5) THEN
            SVAR=1.D0/(14.D0*FOURPI)
            RAWJPAR=RAWJPAR+ULITTLE(5,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
          END IF
          IF(RAWJPAR.GT.0.D0) THEN
            SVAR=JPAR/RAWJPAR
            ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE) &
     &                         =ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
          END IF
        ELSE IF(LOX(LNPROBE).EQ.3) THEN
          RAWJPAR=0.D0
          IF(LRX+1.GT.3) THEN
            SVAR=268.D0/(6435.D0*FOURPI)
            RAWJPAR=RAWJPAR+ULITTLE(3,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
          END IF 
          IF(LRX+1.GT.5) THEN
            SVAR=195.D0/(6435.D0*FOURPI)
            RAWJPAR=RAWJPAR+ULITTLE(5,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
          END IF
          IF(LRX+1.GT.7) THEN
            SVAR=250.D0/(6435.D0*FOURPI)
            RAWJPAR=RAWJPAR+ULITTLE(7,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
          END IF
          IF(RAWJPAR.GT.0.D0) THEN
            SVAR=JPAR/RAWJPAR
            ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE) &
     &                       =ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
          END IF
        ELSE
          IF(JPAR.EQ.0) THEN
            ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)=0.D0
          ELSE
            CALL ERROR$MSG('JPAR.NE.0 CAN ONLY BE SET OF D AND F SHELLS')
            CALL ERROR$MSG('LDAPLUSU_MODULITTLEWITHPARMS')
          END IF
        END IF
      END IF
      
!
!     ==========================================================================
!     == SET UP ULITTLE                                                       ==
!     ==========================================================================
!!$
!!$      ULITTLE=0.D0
!!$      DO L=0,LRX
!!$        FAC=0.D0
!!$        IF(L.EQ.0)  FAC=UPAR
!!$        IF(LRX.GE.2)FAC=JPAR*14.D0/(1.D0+FIVEEIGTH)
!!$        IF(LRX.GE.4)FAC=FIVEEIGTH*JPAR*14.D0/(1.D0+FIVEEIGTH)
!!$        ULITTLE(L+1,:,:,:,:)=FAC*4.D0*PI/REAL(2*L+1,KIND=8)
!!$      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_UTENSOR(LRX,NORB,LNX,LOX,ULITTLE,U)
!     **                                                                      **
!     ** CALCULATES THE U-TENSOR                                              **
!     **                                                                      **
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
                            CALL TRACE$PUSH('LDAPLUSU_UTENSOR')
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
      SUBROUTINE LDAPLUSU_DCLDAPLUSU(ID,LNX,LOX,NCHI,U,RHO,ETOT,HAM)
!     **                                                                      **
!     **  CALCULATES THE CORRELATION ENERGY FROM UTENSOR AND DENSITY MATRIX   **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID     ! ID FOR DOUBLE COUNTING CORRECTION
      INTEGER(4)  ,INTENT(IN) :: LNX    !                     
      INTEGER(4)  ,INTENT(IN) :: LOX(LNX) ! ANGULAR MOMENTUM  
      INTEGER(4)  ,INTENT(IN) :: NCHI   ! BASIS-SET SIZE              
      REAL(8)     ,INTENT(IN) :: U(NCHI,NCHI,NCHI,NCHI) ! U TENSOR
      COMPLEX(8)  ,INTENT(IN) :: RHO(NCHI,NCHI,2,2) ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT    ! DOUBLE COUNTINNG ENERGY
      COMPLEX(8)  ,INTENT(OUT):: HAM(NCHI,NCHI,2,2)  ! DE/D(RHO^*)        
      REAL(8)                 :: UPAR   ! U-PARAMETER
      REAL(8)                 :: JPAR   ! J-PARAMETER
      REAL(8)                 :: E
      REAL(8)                 :: F(2)
      REAL(8)                 :: V(2)
      INTEGER(4)              :: IS,I,J,LN
      INTEGER(4)              :: NCHI1,NCHI2
      INTEGER(4)              :: L
!     **************************************************************************
                            CALL TRACE$PUSH('LDAPLUSU_ETOT')
      ETOT=0.D0
      HAM(:,:,:,:)=(0.D0,0.D0)
      NCHI1=1
      DO LN=1,LNX
        L=LOX(LNX)
        NCHI2=NCHI1-1+(2*L+1)
!
!       ========================================================================
!       ==  CALCULATE U AND J PARAMETERS                                      ==
!       ========================================================================
        UPAR=0.D0
        JPAR=0.D0
        DO I=NCHI1,NCHI2
          DO J=NCHI1,NCHI2
            UPAR=UPAR+U(I,J,I,J)
            JPAR=JPAR+U(I,J,I,J)-U(I,J,J,I)
          ENDDO
        ENDDO
        UPAR=UPAR/REAL((2*L+1)**2,KIND=8)
        JPAR=UPAR-JPAR/REAL(2*L*(2*L+1),KIND=8)
PRINT*,'UPARAMETER[EV]    ',UPAR*27.211D0 ,'UPARAMETER    ',UPAR
PRINT*,'JPARAMETER[EV](1) ',JPAR*27.211D0 ,'JPARAMETER(1) ',JPAR
!
!       ========================================================================
!       ==  DOUBLE COUNTING CORRECTION                                        ==
!       ========================================================================
        DO IS=1,2
          F(IS)=0.D0
          DO I=NCHI1,NCHI2
            F(IS)=F(IS)+REAL(RHO(I,I,IS,IS))
          ENDDO
        ENDDO
        CALL LDAPLUSU_DOUBLECOUNTING(ID,L,UPAR,JPAR,F,E,V)
        ETOT=ETOT+E
        DO IS=1,2
          DO I=NCHI1,NCHI2
            HAM(I,I,IS,IS)=HAM(I,I,IS,IS)+CMPLX(V(IS),0.D0,KIND=8)
          ENDDO
        ENDDO
        NCHI1=NCHI2+1
      ENDDO
                            CALL TRACE$POP()
      RETURN
      END        
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_INTERACTION(NORB,U,RHO,ETOT,HAM)
!     **************************************************************************
!     ** CALCULATES THE HARTREE AND EXCHANGE ENERGY FROM THE U-TENSOR U       **
!     ** AND THE DENSITY MATRIX RHO.                                          **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NORB                   ! BASIS-SET SIZE    
      REAL(8)     ,INTENT(IN) :: U(NORB,NORB,NORB,NORB) ! U TENSOR
      COMPLEX(8)  ,INTENT(IN) :: RHO(NORB,NORB,2,2)     ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT                   ! ENERGY
      COMPLEX(8)  ,INTENT(OUT):: HAM(NORB,NORB,2,2)     ! DE/D(DENMAT)        
      REAL(8)                 :: UIJKL
      COMPLEX(8)              :: RHOT(NORB,NORB)
      COMPLEX(8)              :: HAMT(NORB,NORB)
      REAL(8)                 :: SVAR,EBLOCK
      INTEGER(4)              :: I,J,K,L,IS1,IS2
!     **************************************************************************
                            CALL TRACE$PUSH('LDAPLUSU_INTERACTION')
!     == SPLIT OFF TOTAL DENSITY FOR HARTREE TERM  =============================
      RHOT(:,:)=RHO(:,:,1,1)+RHO(:,:,2,2)
!     == DETERMINE INTERACTION ENERGY ==========================================
      ETOT=0.D0
      HAMT(:,:)=(0.D0,0.D0)
      HAM(:,:,:,:)=(0.D0,0.D0)
      DO L=1,NORB
        DO K=1,NORB
          EBLOCK=0.D0
          DO J=1,NORB
            DO I=1,NORB
              UIJKL=0.5D0*U(I,J,K,L)
!             == HARTREE TERM ==================================================
              SVAR=REAL(RHOT(K,I)*RHOT(L,J))     
              HAMT(K,I)=HAMT(K,I)+UIJKL*RHOT(L,J)
              HAMT(L,J)=HAMT(L,J)+UIJKL*RHOT(K,I)
!             == EXCHANGE TERM =================================================
              DO IS1=1,2
                DO IS2=1,2
                  SVAR=SVAR-REAL(RHO(L,I,IS2,IS1)*RHO(K,J,IS1,IS2))  !
                  HAM(L,I,IS2,IS1)=HAM(L,I,IS2,IS1)-UIJKL*RHO(K,J,IS1,IS2)
                  HAM(K,J,IS1,IS2)=HAM(K,J,IS1,IS2)-UIJKL*RHO(L,I,IS2,IS1)
                ENDDO
              ENDDO
              EBLOCK=EBLOCK+SVAR*UIJKL
            ENDDO
          ENDDO
          ETOT=ETOT+EBLOCK
        ENDDO
      ENDDO
!     ==  ADD CONTRIBUTION FROM HARTREE TERM ===================================
      HAM(:,:,1,1)=HAM(:,:,1,1)+HAMT(:,:)
      HAM(:,:,2,2)=HAM(:,:,2,2)+HAMT(:,:)
                            CALL TRACE$POP()
      RETURN
      END
!!$!
!!$!     ...1.........2.........3.........4.........5.........6.........7.........8
!!$      SUBROUTINE LDAPLUSU_INTERACTION_TRUEHF(NORB,U,RHO,ETOT,HAM)
!!$!     **                                                                      **
! THIS WAS AN ATTEMPT TO DETERMINE THE INTERACTRION ENERGY FROM A 
! UNCORRELATED MULTIDETERMINANT WAVEFUNTION. IT IS UNFINISHED!!!!
!!$!     ** CALCULATES THE HARTREE AND EXCHANGE ENERGY FROM THE U-TENSOR U       **
!!$!     ** AND THE DENSITY MATRIX RHO.                                          **
!!$!     **                                                                      **
!!$      IMPLICIT NONE
!!$      INTEGER(4)  ,INTENT(IN) :: NORB                   ! BASIS-SET SIZE              
!!$      REAL(8)     ,INTENT(IN) :: U(NORB,NORB,NORB,NORB) ! U TENSOR
!!$      COMPLEX(8)  ,INTENT(IN) :: RHO(NORB,NORB,2,2)     ! DENSITY MATRIX
!!$      REAL(8)     ,INTENT(OUT):: ETOT                   ! ENERGY
!!$      COMPLEX(8)  ,INTENT(OUT):: HAM(NORB,NORB,2,2)     ! DE/D(DENMAT)        
!!$      REAL(8)                 :: UIJKL
!!$      COMPLEX(8)              :: RHOT(NORB,NORB)
!!$      COMPLEX(8)              :: HAMT(NORB,NORB)
!!$      REAL(8)                 :: SVAR,EBLOCK
!!$      INTEGER(4)              :: I,J,K,L,IS1,IS2
!!$!     **************************************************************************
!!$                            CALL TRACE$PUSH('LDAPLUSU_INTERACTION')
!!$!
!!$!     ==========================================================================
!!$!     == DIAGONALIZE DENSITY MATRIX                                           ==
!!$!     ==========================================================================
!!$      DENMAT(1:NORB,1:NORB)=RHO(:,:,1,1)
!!$      DENMAT(1:NORB,NORB+1:2*NORB)=RHO(:,:,1,2)
!!$      DENMAT(NORB+1:2*NORB,1:NORB)=RHO(:,:,2,1)
!!$      DENMAT(NORB+1:2*NORB,NORB+1:2*NORB)=RHO(:,:,2,2)
!!$      CALL LIB$DIAGC8(2*NORB,DENMAT,EIG,UT)
!!$      DO I=1,NCHI
!!$        FN(I)=EIG(NCHI+1-I)
!!$      ENDDO
!!$      FN(NCHI+1)=0.D0
!!$      FN(0)=1.D0
!!$!
!!$!     ==========================================================================
!!$!     == TWO PARTICLE DENSITY                                                 ==
!!$!     ==========================================================================
!!$      RHO2(:,:,:,:)=0.D0
!!$      DO N=1,2*NORB
!!$        DO M=1,2*NORB
!!$          WEIGHT=MIN(FN(N),FN(M))-FN(N)*FM(M)
!!$          IF(ABS(WEIGHT).LT.1.D-7) CYCLE
!!$          DO I=1,2*NORB
!!$            ISI=1+MOD(I-1,NORB)
!!$            DO K=1,2*NORB
!!$              ISK=1+MOD(K-1,NORB)
!!$              CSVAR1=WEIGHT*UT(K,N)*UT(I,N)
!!$              DO J=1,2*NORB
!!$                ISJ=1+MOD(J-1,NORB)
!!$                DO L=1,2*NORB
!!$                  ISL=1+MOD(L-1,NORB)
!!$                  RHO2(I,J,K,L)=U(I,J,K,L)+CSVAR*UT(L,M)*UT(J,M)
!!$                ENDDO
!!$              ENDDO
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     ==  CORRECTION TO THE INTERACTION ENERGY                                ==
!!$!     ==========================================================================
!!$      ETOT=0.D0
!!$      HAMMAT(:,:)=(0.D0,0.D0)
!!$      DO I=1,NORB
!!$        DO J=1,NORB
!!$          EBLOCK=0.D0
!!$          DO K=1,NORB
!!$            DO L=1,NORB
!!$              DU=U(I,J,K,L)-U(I,J,L,K)
!!$              DO IS1=0,1
!!$                DO IS2=0,1
!!$                  I1=I+NORB*IS1
!!$                  I2=J+NORB*IS2
!!$                  I3=K+NORB*IS2
!!$                  I4=L+NORB*IS1
!!$                  ETOT=ETOT+DU*RHO2(I1,I2,I3,I4)
!!$                  HAM2(I1,I2,I3,I4)=HAM2(I1,I2,I3,I4)+DU
!!$                ENDDO
!!$              ENDDO
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     ==  TRANSFORM TO DERIVATIVES FO OCCUPATIONS ANE EIGENVECTORS            ==
!!$!     ==========================================================================
!!$      RHO2(:,:,:,:)=0.D0
!!$      DO N=1,2*NORB
!!$        DO M=1,2*NORB
!!$          WEIGHT=MIN(FN(N),FN(M))-FN(N)*FM(M)
!!$          IF(ABS(WEIGHT).LT.1.D-7) CYCLE
!!$          DO I=1,2*NORB
!!$            ISI=1+MOD(I-1,NORB)
!!$            DO K=1,2*NORB
!!$              ISK=1+MOD(K-1,NORB)
!!$              CSVAR1=WEIGHT*UT(K,N)*UT(I,N)
!!$              DO J=1,2*NORB
!!$                ISJ=1+MOD(J-1,NORB)
!!$                DO L=1,2*NORB
!!$                  ISL=1+MOD(L-1,NORB)
!!$                  RHO2(I,J,K,L)=U(I,J,K,L)+CSVAR*UT(L,M)*UT(J,M)
!!$!----                  DEDF(N)
!!$                ENDDO
!!$              ENDDO
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!!$!
!!$!     ==========================================================================
!!$!     ==  TRANSFORM TWO-PARTICLE HAMILTONIAN INTO A ONE-PARTICLE HAMILTONIAN  ==
!!$!     ==========================================================================
!!$                            CALL TRACE$POP()
!!$      RETURN
!!$      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_DOUBLECOUNTING(ID,L,U,J,F,E,V)
!     **                                                                      **
!     ** DOUBLE COUNTING CORRECTION TO THE LDA+U TOTAL ENERGY                 **
!     ** THE ENERGY PRODUCED SHALL BE ADDED TO THE TOTAL ENERGY               **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID   ! SWITCH BETWEEN DIFFERENT FORMULATIONS
      INTEGER(4)  ,INTENT(IN) :: L    ! MAIN ANGULAR MOMENTUM
      REAL(8)     ,INTENT(IN) :: U    ! U-PARAMETER
      REAL(8)     ,INTENT(IN) :: J    ! J-PARAMETER
      REAL(8)     ,INTENT(IN) :: F(2) ! MEAN OCCUPATION/PER SPIN
      REAL(8)     ,INTENT(OUT):: E    ! DOUBLE COUNTINNG ENERGY
      REAL(8)     ,INTENT(OUT):: V(2)  ! DE/DF
      REAL(8)                 :: FTOT,VTOT,SVAR
!     **************************************************************************
                            CALL TRACE$PUSH('LDAPLUSU_DOUBLECOUNTING')
!     ==========================================================================
!     ==  OPTION                                                              ==
!     ==========================================================================
      IF(ID.EQ.'FLL') THEN
        FTOT=F(1)+F(2)
        E=0.5D0*U*FTOT*(FTOT-1.D0) &
     &   -0.5D0*J*F(1)*(F(1)-1.D0) &
     &   -0.5D0*J*F(2)*(F(2)-1.D0)
        VTOT=0.5D0*U*(2.D0*FTOT-1.D0)
        V(1)=VTOT-0.5D0*J*(2.D0*F(1)-1.D0)
        V(2)=VTOT-0.5D0*J*(2.D0*F(2)-1.D0)
!
!     ==========================================================================
!     ==  OPTION APPROXIMATE MEAN FIELD (AMF)                                 ==
!     ==========================================================================
      ELSE IF(ID.EQ.'AMF') THEN
        SVAR=REAL(L)/REAL(2*L+1)*(U-J)
        E=U*F(1)*F(2)+SVAR*(F(1)**2+F(2)**2)
        V(1)=U*F(2)+2.D0*SVAR*F(1)
        V(2)=U*F(1)+2.D0*SVAR*F(2)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('LDAPLUSU_DOUBLECOUNTING')
      END IF
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_EHARTREE(GID,NR,LMNX,LNX,LOX,CHI,LRX,D,E,H)
!     **                                                                      **
!     **  DOUBLE COUNTING CORRECTION FOR THE HYBRID FUNCTIONAL                **
!     **                                                                      **
!     **  DETERMINES THE HARTREE AND EXCHANGE-ONLY ENERGY FROM THE            **
!     **  DFT FUNCTIONAL                                                      **
!     **  FOR THE DENSITY BUILT FROM THE LOCAL ORBITALS AND THE CORE DENSITY  **
!     **  THIS ENERGY NEEDS TO BE SUBTRACTED FROM THE TOTAL ENERGY            **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: LRX
      INTEGER(4)  ,INTENT(IN) :: LMNX       ! #(LOCAL ORBITALS)
      INTEGER(4)  ,INTENT(IN) :: LNX        ! #(RADIAL FUNCTIONS)
      INTEGER(4)  ,INTENT(IN) :: LOX(LNX)   !MAIN ANGULAR MOMENTUM OF LOCAL ORB.
      REAL(8)     ,INTENT(IN) :: CHI(NR,LNX)
      COMPLEX(8)  ,INTENT(IN) :: D(LMNX,LMNX,2,2) ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: E                ! HARTREE ENERGY
      COMPLEX(8)  ,INTENT(OUT):: H(LMNX,LMNX,2,2) ! DETOT/D(RHO^*)        
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)              :: DENMAT1(LMNX,LMNX)
      COMPLEX(8)              :: HAM1(LMNX,LMNX)
      REAL(8)                 :: R(NR)
      REAL(8)     ,ALLOCATABLE:: RHO(:,:)
      REAL(8)     ,ALLOCATABLE:: POT(:,:)
      REAL(8)                 :: EDENSITY(NR)
      INTEGER(4)              :: LMRX,L
      INTEGER(4)              :: IDIM,LM
!     **************************************************************************
      LMRX=(LRX+1)**2
!
!     ==========================================================================
!     ==  TRANSFORM DENSITY MATRIX FROM UP/DOWN TO TOTAL/SPIN                 ==
!     ==========================================================================
      DENMAT1(:,:)=D(:,:,1,1)+D(:,:,2,2)
!
!     ==========================================================================
!     ==  CALCULATE DENSITY                                                   ==
!     ==========================================================================
      ALLOCATE(RHO(NR,LMRX))
      CALL AUGMENTATION_RHO(NR,LNX,LOX,CHI,LMNX,DENMAT1,LMRX,RHO)
!
!     ==========================================================================
!     == HARTREE ENERGY AND POTENTIAL ==========================================
!     == CORE CONTRIBUTION IS NOT INCLUDED BECAUSE IT IS NOT REPRESENTED IN   ==
!     == THE U-TENSOR AND ONLY THE EXCHANGE PART OF THE CORE-VALENCE IS INCLUDED
!     ==========================================================================
      ALLOCATE(POT(NR,LMRX))
      EDENSITY=0.D0
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        CALL RADIAL$POISSON(GID,NR,L,RHO(:,LM),POT(:,LM))
        EDENSITY(:)=EDENSITY(:)+0.5D0*POT(:,LM)*RHO(:,LM)
      ENDDO
      CALL RADIAL$R(GID,NR,R)
      EDENSITY=EDENSITY*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,EDENSITY,E)
!
!     ==========================================================================
!     ==  CALCULATE HAMILTONIAN IN TOTAL/SPIN REPRESENTATION                  ==
!     ==========================================================================
      CALL LDAPLUSU_EXPECT(GID,NR,1,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
      DEALLOCATE(POT)
!
!     ==========================================================================
!     ==  TRANSFORM HAMILTONIAN FROM TOTAL/SPIN TO UP/DOWN                    ==
!     ==========================================================================
      H(:,:,:,:)=(0.D0,0.D0)
      H(:,:,1,1)=HAM1(:,:)
      H(:,:,2,2)=HAM1(:,:)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_EDFT(GID,NR,LMNX,LNX,LOX,CHI,LRX,AECORE,DENMAT,ETOT,HAM)
!     **                                                                      **
!     **  DOUBLE COUNTING CORRECTION FOR THE HYBRID FUNCTIONAL                **
!     **                                                                      **
!     **  DETERMINES THE HARTREE AND EXCHANGE-ONLY ENERGY FROM THE            **
!     **  DFT FUNCTIONAL                                                      **
!     **  FOR THE DENSITY BUILT FROM THE LOCAL ORBITALS AND THE CORE DENSITY  **
!     **  THIS ENERGY NEEDS TO BE SUBTRACTED FROM THE TOTAL ENERGY            **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: LRX
      INTEGER(4)  ,INTENT(IN) :: LMNX       ! #(LOCAL ORBITALS)
      INTEGER(4)  ,INTENT(IN) :: LNX        ! #(RADIAL FUNCTIONS)
      INTEGER(4)  ,INTENT(IN) :: LOX(LNX)   !MAIN ANGULAR MOMENTUM OF LOCAL ORB.
      REAL(8)     ,INTENT(IN) :: CHI(NR,LNX)
      REAL(8)     ,INTENT(IN) :: AECORE(NR)
      COMPLEX(8)  ,INTENT(IN) :: DENMAT(LMNX,LMNX,2,2) ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT       ! DOUBLE COUNTINNG ENERGY
      COMPLEX(8)  ,INTENT(OUT):: HAM(LMNX,LMNX,2,2)  ! DETOT/D(RHO^*)        
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
      DENMAT1(:,:,1)=DENMAT(:,:,1,1)+DENMAT(:,:,2,2)
      DENMAT1(:,:,2)=DENMAT(:,:,1,2)+DENMAT(:,:,2,1)
      DENMAT1(:,:,3)=CI*(DENMAT(:,:,1,2)-DENMAT(:,:,2,1))
      DENMAT1(:,:,4)=DENMAT(:,:,1,1)-DENMAT(:,:,2,2)
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
PRINT*,'TOTAL        EXCHANGE ENERGY (LOCAL) ',ETOT
PRINT*,'VALENCE      EXCHANGE ENERGY (LOCAL) ',ETOTV
PRINT*,'CORE         EXCHANGE ENERGY (LOCAL) ',ETOTC
PRINT*,'CORE-VALENCE EXCHANGE ENERGY (LOCAL) ',ETOT-ETOTV-ETOTC
      ETOT=ETOT-ETOTC
IF(ETOT.LT.-3.145D0) THEN
PRINT*,'FILE RHOWC.DAT WRITTEN'
CALL ATOMLIB_WRITEPHI('RHOWC1.DAT',GID,NR,LMRX,RHOWC(:,:,1))
CALL ATOMLIB_WRITEPHI('RHOWC2.DAT',GID,NR,LMRX,RHOWC(:,:,2))
CALL ATOMLIB_WRITEPHI('RHOWC3.DAT',GID,NR,LMRX,RHOWC(:,:,3))
CALL ATOMLIB_WRITEPHI('RHOWC4.DAT',GID,NR,LMRX,RHOWC(:,:,4))
END IF

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
PRINT*,'EDFT: EXC ',ETOT
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
      HAM(:,:,1,1)=HAM1(:,:,1)+   HAM1(:,:,4)
      HAM(:,:,1,2)=HAM1(:,:,2)+CI*HAM1(:,:,3)
      HAM(:,:,2,1)=HAM1(:,:,2)-CI*HAM1(:,:,3)
      HAM(:,:,2,2)=HAM1(:,:,1)-   HAM1(:,:,4)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &                          ,AEPOT,AEPHI,DATH)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE EXPECTATION VALUE OF                                 **
!     **  THE ONE-CENTER POTENTIALS WITH THE LOCAL ORBITALS                   **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NDIMD
      INTEGER(4),INTENT(IN) :: LMRX
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      INTEGER(4),INTENT(IN) :: LMNX
      REAL(8)   ,INTENT(IN) :: AEPOT(NR,LMRX,NDIMD)
      REAL(8)   ,INTENT(IN) :: AEPHI(NR,LNX)
      COMPLEX(8),INTENT(OUT):: DATH(LMNX,LMNX,NDIMD)
      INTEGER(4)            :: LMN1,LMN2
      INTEGER(4)            :: LN1,LN2
      INTEGER(4)            :: LM1,LM2,LM3
      INTEGER(4)            :: L1,L2
      INTEGER(4)            :: IM1,IM2
      INTEGER(4)            :: ISPIN
      REAL(8)               :: AEDMU(NR,NDIMD)
      REAL(8)               :: DWORK1(NR)
      REAL(8)               :: CG
      REAL(8)               :: SVAR
      REAL(8)               :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      DATH(:,:,:)=0.D0
      LMN1=0
      DO LN1=1,LNX
        L1=LOX(LN1)
        DO IM1=1,2*L1+1
          LMN1=LMN1+1
          LMN2=0
          LM1=L1**2+IM1
          DO LN2=1,LNX
            L2=LOX(LN2)
            DO IM2=1,2*L2+1
              LMN2=LMN2+1
              LM2=L2**2+IM2
!     
!             ==================================================================
!             ==  SUM ALL POTENTIALS THAT ACT ON THE GIVEN PAIR               ==
!             ==  OF PARTIAL WAVES                                            ==
!             ==================================================================
              AEDMU(:,:)=0.D0
              DO LM3=1,LMRX
                CALL CLEBSCH(LM1,LM2,LM3,CG)
                IF(CG.NE.0.D0) THEN
                  DO ISPIN=1,NDIMD
                    AEDMU(:,ISPIN)=AEDMU(:,ISPIN)+CG*AEPOT(:,LM3,ISPIN)
                  ENDDO
                END IF
              ENDDO
!     
!             ==================================================================
!             ==  PERFORM NOW THE INTEGRATION                                 ==
!             ==================================================================
              DO ISPIN=1,NDIMD
                DWORK1(:)=AEDMU(:,ISPIN)*AEPHI(:,LN1)*AEPHI(:,LN2)*R(:)**2
                CALL RADIAL$INTEGRAL(GID,NR,DWORK1,SVAR)
                DATH(LMN1,LMN2,ISPIN)=DATH(LMN1,LMN2,ISPIN)+CMPLX(SVAR,0.D0,8)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_JEXPANSION(LSHELL,LRHOX,XL)
!     **                                                                      **
!     **  DETERMINES THE EXPANSION COEFFICIENTS OF THE J-PARAMETER            **
!     **  IN SLATER INTEGRALS. THE J-PARAMETER REFERS TO THE SHELL WITH       **
!     **  MAIN ANGULAR MOMETUM "LSHELL". THE EXPANSION COEFFICIENTS ARE       **
!     **  DETERMINED UP TO ANGULAR MOMENTUM LRHOX.                            **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LSHELL
      INTEGER(4),INTENT(IN) :: LRHOX
      REAL(8)   ,INTENT(OUT):: XL(LRHOX+1)
      REAL(8)               :: PI
      INTEGER(4)            :: LPRIME,M1,M2,MPRIME,LM1,LM2,LMPRIME
      REAL(8)               :: SVAR,CG
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      XL(:)=0.D0
      DO LPRIME=1,LRHOX
        SVAR=0.D0
        DO M1=1,2*LSHELL+1
          DO M2=1,2*LSHELL+1
            LM1=LSHELL**2+M1
            LM2=LSHELL**2+M2
            DO MPRIME=1,2*LPRIME+1
              LMPRIME=LPRIME**2+MPRIME
              CALL CLEBSCH(LM1,LM2,LMPRIME,CG)
              SVAR=SVAR+CG**2
            ENDDO
          ENDDO
        ENDDO
        XL(LPRIME+1)=4.D0*PI/REAL(2*LPRIME+1,KIND=8)*SVAR
      ENDDO
      XL(:)=XL(:)/REAL(2*LSHELL*(2*LSHELL+1),KIND=8)
!
!     ==========================================================================
!     == REPORT RESULT                                                        ==
!     ==========================================================================
!!$      WRITE(*,FMT='("J-PARAMETER EXPANSION COEFFICIENTS IN SLATER INTEGRALS")')
!!$      WRITE(*,FMT='("ANGULAR MOMENTUM SHELL FOR J: L=",I5)')LSHELL
!!$      DO LPRIME=1,LRHOX
!!$        WRITE(*,FMT='("X(",I2,")=",F10.5))')LPRIME,XL(LPRIME+1)
!!$      ENDDO
!!$STOP
      RETURN
      END

!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_DENMATFLAPW(ID,NPHI,UPFOLDED,NCHI,DOWNFOLDED)
!     **************************************************************************
!     **                                                                      **
!     **  DETERMINES THE ORBITAL OCCUPATIONS FOR A CONVENTIONAL  LDA+U TYPE   **
!     **  CALCULATION ACCORDING TO EQ. 9 OF BENGONE ET AL. PRB62, 16392 (2000)**
!     **  AND ITS BACK TRANSFORM                                              **
!     **                                                                      **
!     **************************************************************************
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID   ! CAN BE 'FORWARD' OR 'BACK'
      INTEGER(4)  ,INTENT(IN)    :: NPHI
      INTEGER(4)  ,INTENT(IN)    :: NCHI 
      COMPLEX(8)  ,INTENT(INOUT) :: UPFOLDED(NPHI,NPHI,2,2)
      COMPLEX(8)  ,INTENT(INOUT) :: DOWNFOLDED(NCHI,NCHI,2,2)
      INTEGER(4)                 :: GID
      INTEGER(4)                 :: NR
      INTEGER(4)                 :: LNXCHI
      INTEGER(4)                 :: LNXPHI
      INTEGER(4)  ,ALLOCATABLE   :: LOXPHI(:)
      INTEGER(4)  ,ALLOCATABLE   :: LOXCHI(:)
      REAL(8)     ,ALLOCATABLE   :: R(:)       ! RADIAL GRID
      REAL(8)     ,ALLOCATABLE   :: PHI(:,:)   ! PARTIAL WAVES
      REAL(8)                    :: RCUT
      REAL(8)     ,ALLOCATABLE   :: OVER(:,:)
      REAL(8)     ,ALLOCATABLE   :: AUX1(:),AUX2(:)
      REAL(8)                    :: VAL
      INTEGER(4)                 :: LN1,LN2      
      INTEGER(4)                 :: L1,L2      
      INTEGER(4)                 :: M1,M2     
      INTEGER(4)                 :: IS1,IS2
      INTEGER(4)                 :: LMN1,LMN2
      INTEGER(4)                 :: LCHI
!     **************************************************************************
      IF(ID.NE.'FORWARD'.AND. ID.NE.'BACK') THEN
        CALL ERROR$MSG('ILLEGAL VALUE FOR ID')
        CALL ERROR$MSG('ALLOWED VALUES ARE "FORWARD" AND "BACK"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU_DENMATFLAPW')
      END IF
      GID=THIS%SP%GID
      NR=THIS%SP%NR
      LNXCHI=THIS%SP%LNXCHI
      ALLOCATE(LOXCHI(LNXCHI))
      LOXCHI=THIS%SP%LOXCHI
      RCUT=THIS%SP%RCUT
!
!     ==========================================================================
!     ==  CHECK IF ONLY ONE CORRELATED ORBITAL PER L                          ==
!     ==========================================================================
      IF(LNXCHI.NE.1) THEN
        CALL ERROR$STOP('LDAPLUSU_DENMATFLAPW')
      END IF
!
!     ==========================================================================
!     ==  DETERMINE PROJECTED OVERLAP MATRIX ELEMENTS                         ==
!     ==========================================================================
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R) 
      CALL SETUP$GETI4('LNX',LNXPHI)
      ALLOCATE(LOXPHI(LNXPHI))
      CALL SETUP$GETI4A('LOX',LNXPHI,LOXPHI)
      ALLOCATE(PHI(NR,LNXPHI))
      CALL SETUP$GETR8A('AEPHI',NR*LNXPHI,PHI)
      ALLOCATE(OVER(LNXPHI,LNXPHI))
      ALLOCATE(AUX1(NR))
      ALLOCATE(AUX2(NR))
      OVER(:,:)=0.D0
      DO LN1=1,LNXPHI
        L1=LOXPHI(LN1)
        DO LN2=LN1,LNXPHI
          L2=LOXPHI(LN2)
          IF(L1.EQ.L2) THEN
            AUX1(:)=R(:)**2*PHI(:,LN1)*PHI(:,LN2)
            CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
            CALL RADIAL$VALUE(GID,NR,AUX2,RCUT,VAL)
            OVER(LN1,LN2)=VAL
            OVER(LN2,LN1)=VAL
          END IF
        ENDDO
      ENDDO
      DEALLOCATE(AUX1)
      DEALLOCATE(AUX2)
      DEALLOCATE(PHI)
      DEALLOCATE(R)
!
!     ==========================================================================
!     ==                                                                      ==
!     ==========================================================================
      IF(ID.EQ.'FORWARD') THEN
        DOWNFOLDED(:,:,:,:)=0.D0
      ELSE
        UPFOLDED(:,:,:,:)=0.D0
      END IF
      LCHI=LOXCHI(1)
      DO IS1=1,2
        DO IS2=1,2
          LMN1=0
          DO LN1=1,LNXPHI                
            L1=LOXPHI(LN1)
            DO M1=1,2*L1+1
              LMN1=LMN1+1
              LMN2=0
              DO LN2=1,LNXPHI                
                L2=LOXPHI(LN2)
                DO M2=1,2*L2+1
                  LMN2=LMN2+1
                  IF(L1.NE.LCHI.OR.L2.NE.LCHI) CYCLE
                  IF(ID.EQ.'FORWARD') THEN
                    DOWNFOLDED(M1,M2,IS1,IS2)=DOWNFOLDED(M1,M2,IS1,IS2) &
    &                                 +UPFOLDED(LMN1,LMN2,IS1,IS2)*OVER(LN1,LN2)
                  ELSE
                    UPFOLDED(LMN1,LMN2,IS1,IS2)=DOWNFOLDED(M1,M2,IS1,IS2) &
    &                                          *OVER(LN1,LN2)
                  END IF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO                    
      DEALLOCATE(OVER)
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_CVX(NCHI,LNX,LOX,RHO,DH,ETOT,HAM)
!     **************************************************************************
!     ** ADD THE CORE VALENCE EXCHANGE POTENTIAL                              **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      INTEGER(4),INTENT(IN) :: LNX
      INTEGER(4),INTENT(IN) :: LOX(LNX)
      COMPLEX(8),INTENT(IN) :: RHO(NCHI,NCHI,2,2)
      REAL(8)   ,INTENT(IN) :: DH(LNX,LNX)
      REAL(8)   ,INTENT(OUT):: ETOT
      COMPLEX(8),INTENT(OUT):: HAM(NCHI,NCHI,2,2)
      INTEGER(4)            :: LN1,LN2,LMN1,LMN2,LMN1A,LMN2A
      INTEGER(4)            :: L,IM
!     **************************************************************************
      ETOT=0.D0
      HAM(:,:,:,:)=(0.D0,0.D0)
      DO LN1=1,LNX
        L=LOX(LN1)
        DO LN2=1,LNX
          IF(LOX(LN2).NE.L) CYCLE
          LMN1A=SUM(2*LOX(1:LN1-1)+1)
          LMN2A=SUM(2*LOX(1:LN2-1)+1)
          DO IM=1,2*L+1
            LMN1=LMN1A+IM
            LMN2=LMN2A+IM
            ETOT=ETOT+(RHO(LMN1,LMN2,1,1)+RHO(LMN1,LMN2,2,2))*DH(LN2,LN1)
            HAM(LMN1,LMN2,1,1)=HAM(LMN1,LMN2,1,1)+DH(LN1,LN2)
            HAM(LMN1,LMN2,2,2)=HAM(LMN1,LMN2,2,2)+DH(LN1,LN2)
         ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_CI(NCHI,U,RHO,ETOT,POT)
!     **************************************************************************
!     **   ACTUAL MAIN DRIVER ROUTINE FOR CI 
!     **************************************************************************
      USE CI_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NCHI
      REAL(8)   ,INTENT(IN) :: U(NCHI,NCHI,NCHI,NCHI)      
      COMPLEX(8),INTENT(IN) :: RHO(NCHI,NCHI,2,2)
      REAL(8)   ,INTENT(OUT):: ETOT
      COMPLEX(8),INTENT(OUT):: POT(NCHI,NCHI,2,2)
      COMPLEX(8)            :: RHO2(2*NCHI,2*NCHI)
      COMPLEX(8)            :: POT2(2*NCHI,2*NCHI)
      COMPLEX(8)            :: CSVAR
      INTEGER(4)            :: I,J,K,L,N
      TYPE(CIHAMIL_TYPE)    :: CIHAM
      TYPE(CISTATE_TYPE)    :: CIPSI
      REAL(8)               :: FN(2*NCHI)
      COMPLEX(8)            :: TRANSF(2*NCHI,2*NCHI)
      COMPLEX(8)            :: RHO2A(2*NCHI,2*NCHI)
      COMPLEX(8),ALLOCATABLE:: U1(:,:,:,:)
      COMPLEX(8),ALLOCATABLE:: USVAR(:,:,:,:)
      COMPLEX(8),ALLOCATABLE:: POTSVAR(:,:)
!     **************************************************************************
      ETOT=0.D0
      POT(:,:,:,:)=(0.D0,0.D0)
!
!     ==========================================================================
!     == MAP DENSITY MATRIX ONTO RANK 2                                       ==
!     ==========================================================================
      RHO2(1:NCHI,1:NCHI)=RHO(:,:,1,1)
      RHO2(1:NCHI,NCHI+1:2*NCHI)=RHO(:,:,1,2)
      RHO2(NCHI+1:2*NCHI,1:NCHI)=RHO(:,:,2,1)
      RHO2(NCHI+1:2*NCHI,NCHI+1:2*NCHI)=RHO(:,:,2,2)
!
!     ==========================================================================
!     == MAP UTENSOR                                                          ==
!     ==========================================================================
      ALLOCATE(U1(2*NCHI,2*NCHI,2*NCHI,2*NCHI))
      U1=(0.D0,0.D0)
      U1(1:NCHI,1:NCHI,1:NCHI,1:NCHI)=U(:,:,:,:)
      U1(1:NCHI,NCHI+1:2*NCHI,1:NCHI,NCHI+1:2*NCHI)=U(:,:,:,:)
      U1(NCHI+1:2*NCHI,1:NCHI,NCHI+1:2*NCHI,1:NCHI)=U(:,:,:,:)
      U1(NCHI+1:2*NCHI,NCHI+1:2*NCHI,NCHI+1:2*NCHI,NCHI+1:2*NCHI)=U(:,:,:,:)
!
!     ==========================================================================
!     == DIAGONALIZE DENSITY MATRIX                                           ==
!     ==THE CI-EXPANSION CONVERGES FASTER IF FORMULATED IN NATURAL ORBITALS   ==
!     ==========================================================================
      CALL LIB$DIAGC8(2*NCHI,RHO2,FN,TRANSF)
!
!     == SET UP TRANSFORMED DENSITY MATRIX =====================================
      RHO2(:,:)=(0.D0,0.D0)
      DO I=1,2*NCHI
        RHO2(I,I)=FN(I)
      ENDDO
!
!     == TRANSFORM UTENSOR =====================================================
      ALLOCATE(USVAR(2*NCHI,2*NCHI,2*NCHI,2*NCHI))
      USVAR(:,:,:,:)=(0.D0,0.D0)
      DO N=1,2*NCHI
        DO I=1,2*NCHI
          USVAR(:,:,:,N)=USVAR(:,:,:,N)+U1(:,:,:,I)*TRANSF(I,N)
        ENDDO
      ENDDO
      U1(:,:,:,:)=(0.D0,0.D0)
      DO N=1,2*NCHI
        DO I=1,2*NCHI
          U1(:,:,N,:)=U1(:,:,N,:)+USVAR(:,:,I,:)*TRANSF(I,N)
        ENDDO
      ENDDO
      USVAR(:,:,:,:)=(0.D0,0.D0)
      DO N=1,2*NCHI
        DO I=1,2*NCHI
          USVAR(:,N,:,:)=USVAR(:,N,:,:)+U1(:,I,:,:)*CONJG(TRANSF(I,N))
        ENDDO
      ENDDO
      U1(:,:,:,:)=(0.D0,0.D0)
      DO N=1,2*NCHI
        DO I=1,2*NCHI
          U1(N,:,:,:)=U1(N,:,:,:)+USVAR(I,:,:,:)*CONJG(TRANSF(I,N))
        ENDDO
      ENDDO
      DEALLOCATE(USVAR)
!
!     ==========================================================================
!     == SET UP HAMILTON OPERATOR                                             ==
!     ==========================================================================
      DO I=1,2*NCHI
        DO J=1,2*NCHI
          DO K=1,2*NCHI
            DO L=1,2*NCHI
              CSVAR=CMPLX(U1(I,J,K,L),KIND=8)
              IF(ABS(CSVAR).LT.1.D-6) CYCLE
              CALL CI$SETU(CIHAM,I,J,K,L,CSVAR)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL CI$CLEANHAMILTONIAN(CIHAM)
      CALL CI$WRITEHAMILTONIAN(CIHAM,6)
!
!     ==========================================================================
!     == SET UP HAMILTON OPERATOR                                             ==
!     ==========================================================================
      CALL CI$DYNWITHFIXEDDENMAT(2*NCHI,RHO2,CIHAM,CIPSI,POT2,ETOT)
      CALL CI$WRITEPSI(CIPSI,6) 
!
!     ==========================================================================
!     == TRANSFORM FROM NATURAL ORBITALS BACK                                 ==
!     ==========================================================================
      ALLOCATE(POTSVAR(2*NCHI,2*NCHI))
      POTSVAR(:,:)=(0.D0,0.D0)
      DO N=1,2*NCHI
        DO I=1,2*NCHI
          POTSVAR(:,I)=POTSVAR(:,I)+POT2(:,N)*CONJG(TRANSF(N,I))
        ENDDO
      ENDDO
      POT2(:,:)=(0.D0,0.D0)
      DO N=1,2*NCHI
        DO I=1,2*NCHI
          POT2(I,:)=POT2(I,:)+POTSVAR(N,:)*TRANSF(N,I)
        ENDDO
      ENDDO
      DEALLOCATE(POTSVAR)
!
!     ==========================================================================
!     == MAP POTENTIAL BACK                                                   ==
!     ==========================================================================
      POT(:,:,1,1)=POT2(1:NCHI,1:NCHI)
      POT(:,:,1,2)=POT2(1:NCHI,NCHI+1:2*NCHI)
      POT(:,:,2,1)=POT2(NCHI+1:2*NCHI,1:NCHI)
      POT(:,:,2,2)=POT2(NCHI+1:2*NCHI,NCHI+1:2*NCHI)
!
!     ==========================================================================
!     == CLEAN UP                                                             ==
!     ==========================================================================
      CALL CI$DELETEPSI(CIPSI)
      CALL CI$DELETEHAMILTONIAN(CIHAM)
CALL ERROR$STOP('FORCED STOP IN LDAPLUSU_CI')
      RETURN 
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_SINGLEU(GID,NR,LRHOX,LM1,PHI1,LM2,PHI2,LM3,PHI3,LM4,PHI4,U)
!     **************************************************************************
!     **  CALCULATES A SINGLE MATRIX ELEMENT OF THE U-TENSOR                  **
!     **  NO OPTIMIZATION. TO BE USED ONLY FOR TESTING PURPOSES               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LRHOX
      INTEGER(4),INTENT(IN) :: LM1
      REAL(8)   ,INTENT(IN) :: PHI1(NR)
      INTEGER(4),INTENT(IN) :: LM2
      REAL(8)   ,INTENT(IN) :: PHI2(NR)
      INTEGER(4),INTENT(IN) :: LM3
      REAL(8)   ,INTENT(IN) :: PHI3(NR)
      INTEGER(4),INTENT(IN) :: LM4
      REAL(8)   ,INTENT(IN) :: PHI4(NR)
      REAL(8)   ,INTENT(OUT):: U
      REAL(8)               :: PI
      REAL(8)               :: R(NR)
      REAL(8)               :: RHO13(NR),RHO24(NR),POT(NR)
      REAL(8)               :: AUX(NR)
      INTEGER(4)            :: LMRHO,LRHO,IM
      REAL(8)               :: CG1,CG2,VAL
!     *************************************************************************
      PI=4.D0*ATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
      RHO24=PHI2(:)*PHI4(:)
      RHO13=PHI1(:)*PHI3(:)
      U=0.D0
      LMRHO=0
      DO LRHO=0,LRHOX
        CALL RADIAL$POISSON(GID,NR,LRHO,RHO24,POT)
        AUX=R(:)**2*POT(:)*RHO13(:)
        CALL RADIAL$INTEGRAL(GID,NR,AUX,VAL)
!       == SLATER INTEGRAL IS (2*LRHO+1)/(4PI)*VAL
        DO IM=1,2*LRHO+1
          LMRHO=LMRHO+1
          CALL SPHERICAL$GAUNT(LM1,LM3,LMRHO,CG1)
          CALL SPHERICAL$GAUNT(LM2,LM4,LMRHO,CG2)
          U=U+CG1*CG2*VAL
        ENDDO
      ENDDO
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
      CALL SETUP$UNSELECT()
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
!!$      CALL SETUP$ISELECT(ISP)
!!$      CALL SETUP$GETR8('AEZ',AEZ)
!!$      WRITE(STRING,FMT='(F3.0)')AEZ
!!$      STRING=-'_FORZ'//TRIM(ADJUSTL(STRING))//-'DAT'
!!$      CALL SETUP_WRITEPHI(-'CHI'//TRIM(STRING),GID,NR,LNCHI,CHI)
!!$      CALL SETUP$ISELECT(0)
!
                            CALL TRACE$POP()
      RETURN
      END
