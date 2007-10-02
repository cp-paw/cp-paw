!TODO :
! DATH IS STILL REAL AND SHOULD PROBABLY BE COMPLEX LIKE DENMAT
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE LDAPLUSU_MODULE
TYPE THISTYPE
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
CHARACTER(16)          :: FUNCTIONALID !CAN BE 'LDA+U' OR 'HYBRID'
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
!==  INTERMEDIATE STORAGE FROM INITUALIZATION ================================== 
INTEGER(4)             :: NCHI              !#(LOCAL ORBITALS)
REAL(8)   ,POINTER     :: CHI(:,:)          !LOCAL (HEAD) ORBITALS
REAL(8)   ,POINTER     :: ULITTLE(:,:,:,:,:)!SLATER INTEGRALS 
REAL(8)   ,POINTER     :: DOWNFOLD(:,:)     !MAPS PARTIAL WAVES TO LOCAL ORBITALS
END TYPE THISTYPE
TYPE(THISTYPE),ALLOCATABLE,TARGET :: THISARRAY(:)
TYPE(THISTYPE),POINTER :: THIS
LOGICAL(4)             :: TON=.FALSE.
INTEGER(4)             :: NSP=0    ! #(ATOM TYPES)
INTEGER(4)             :: ISP=0    ! SELECTED ATOM TYPE
CHARACTER(8)           :: DCTYPE='FLL' ! SPECIFIES TYPE OF DOUBLE COUNTING CORRECTION
END MODULE LDAPLUSU_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$NEW(NSP_)
!     **                                                                      **
!     **  SETS UP THE ARRAYS FOR OPERATION. THIS ROUTINE MUST BE CALLED ONCE  **
!     **  BEFORE ANY OTHER CALL TO THE LDAPLUSU OBJECT                        **
!     **                                                                      **
!     **  ATTENTION: CHANGES THE SETTING OF SETUP OBJECT                      **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSP_    ! #(DIFFERENT ATOM TYPES)
!     **************************************************************************
!
!     ==========================================================================
!     == DO NOTHING IF THISARRAY IS ALREADY ALLOCATED                         ==
!     ==========================================================================
      IF(NSP.NE.0) THEN
        IF(NSP_.NE.NSP) THEN
          CALL ERROR$MSG('LDAPLUSU$NEW CANNOT BE CALLED TWICE')
          CALL ERROR$STOP('LDAPLUSU$NEW')
        END IF
        RETURN
      END IF
!
!     ==========================================================================
!     == CREATE THISARRAY                                                     ==
!     ==========================================================================
      NSP=NSP_
      ISP=0
      ALLOCATE(THISARRAY(NSP))
      DO ISP=1,NSP
        ALLOCATE(THISARRAY(ISP)%NORB(1))
        THISARRAY(ISP)%NORB(:)=0  ! PER DEFAULT CORRELATE NOTHING
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SELECT(ISP_)
!     **                                                                      **
!     **  SELECT AN ATOM TYPE OR UNSELECT WITH ISP=0                          **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ISP_  ! UNIQUE INDEX OF THE ATOM TYPE
!     **************************************************************************
      ISP=ISP_
      IF(ISP.GT.NSP) THEN
        IF(NSP.EQ.0) THEN
          CALL ERROR$STOP('LDAPLUSU NOT INITIALIZED')
          CALL ERROR$STOP('LDAPLUSU$SELECT')
        END IF
        CALL ERROR$STOP('ISP OUT OF RANGE')
        CALL ERROR$STOP('LDAPLUSU$SELECT')
      END IF
      IF(ISP.EQ.0) THEN
        NULLIFY(THIS)
      ELSE
        THIS=>THISARRAY(ISP)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SETR8(ID,VAL)
!     **                                                                      **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETR8')
      END IF

      IF(ID.EQ.'UPAR') THEN
        THIS%UPAR=VAL
        THIS%USEUPAR=.TRUE.
        IF(THIS%USEDIEL) THEN
          CALL ERROR$MSG('DIEL HAS ALREADY BEEN SET')
          CALL ERROR$MSG('UPAR AND DIEL CANNOT BE USED SIMULTANEOUSLY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
      ELSE IF(ID.EQ.'JPAR') THEN
        THIS%JPAR=VAL
        THIS%USEJPAR=.TRUE.
        IF(THIS%USEDIEL) THEN
          CALL ERROR$MSG('DIEL HAS ALREADY BEEN SET')
          CALL ERROR$MSG('JPAR AND DIEL CANNOT BE USED SIMULTANEOUSLY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
      ELSE IF(ID.EQ.'F4/F2') THEN
        THIS%FRATIO42=VAL
        THIS%USEFRATIO42=.TRUE.
        IF(THIS%USEDIEL) THEN
          CALL ERROR$MSG('DIEL HAS ALREADY BEEN SET')
          CALL ERROR$MSG('FRATIO42 AND DIEL CANNOT BE USED SIMULTANEOUSLY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
      ELSE IF(ID.EQ.'F6/F2') THEN
        THIS%FRATIO62=VAL
        THIS%USEFRATIO62=.TRUE.
        IF(THIS%USEDIEL) THEN
          CALL ERROR$MSG('DIEL HAS ALREADY BEEN SET')
          CALL ERROR$MSG('FRATIO62 AND DIEL CANNOT BE USED SIMULTANEOUSLY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
      ELSE IF(ID.EQ.'DIEL') THEN
        THIS%DIEL=VAL
        THIS%USEDIEL=.TRUE.
        IF(THIS%USEUPAR.OR.THIS%USEJPAR) THEN
          CALL ERROR$MSG('UPAR OR JPAR HAVE ALREADY BEEN SET')
          CALL ERROR$MSG('DIEL AND UPAR OR JPAR CANNOT BE USED SIMULTANEOUSLY')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETR8')
        END IF
      ELSE IF(ID.EQ.'RCUT') THEN
        THIS%RCUT=VAL
      ELSE IF(ID.EQ.'HFWEIGHT') THEN
        THIS%HFWEIGHT=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETR8')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SETL4(ID,VAL)
!     **                                                                      **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      LOGICAL(4)  ,INTENT(IN) :: VAL
!     **************************************************************************
      IF(ID.EQ.'ON') THEN
        TON=VAL
        RETURN
      END IF
      IF(ISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETR8')
      END IF
      IF(ID.EQ.'ACTIVE') THEN
        THIS%TON=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETL4')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SETI4A(ID,LENG,VAL)
!     **                                                                      **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LENG
      INTEGER(4)  ,INTENT(IN) :: VAL(LENG)
!     **************************************************************************
      IF(ISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETI4A')
      END IF
      IF(ID.EQ.'NCORROFL') THEN
        DEALLOCATE(THISARRAY(ISP)%NORB)
        ALLOCATE(THISARRAY(ISP)%NORB(LENG))
        THIS%NORB=VAL
!
      ELSE IF(ID.EQ.'MAINLN') THEN
        THIS%MAINLN=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETI4A')
      END IF
      RETURN
      END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$SETCH(ID,VAL)
!     **                                                                      **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN) :: VAL
!     **************************************************************************
      IF(ISP.EQ.0) THEN
        CALL ERROR$MSG('LDAPLUSU NOT SELECTED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LDAPLUSU$SETCH')
      END IF
      IF(ID.EQ.'FUNCTIONALID') THEN
        IF(VAL.NE.'LDA+U'.AND.VAL.NE.'HYBRID') THEN
          CALL ERROR$MSG('ILLEGAL VALUE FOR FUNCTIONALID')
          CALL ERROR$CHVAL('FUNCTIONALID',VAL)
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('LDAPLUSU$SETCH')
        END IF
        THIS%FUNCTIONALID=VAL
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
      SUBROUTINE LDAPLUSU$REPORT(NFIL)
!     **                                                                      **
!     **  REPORTS THE SETTINGS OF THE LDAPLUSU OBJECT                         **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      USE CONSTANTS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      TYPE(THISTYPE),POINTER :: THIS1
      CHARACTER(64)          :: STRING
      INTEGER(4)             :: L
      REAL(8)                :: EV
!     **************************************************************************
      IF(.NOT.TON) RETURN
      CALL REPORT$TITLE(NFIL,'HYBRID FUNCTIONAL')
      DO ISP=1,NSP
        THIS1=>THISARRAY(ISP)
        IF(.NOT.THIS1%TON) CYCLE
        CALL SETUP$ISELECT(ISP)
        CALL SETUP$GETCH('ID',STRING)
        THISARRAY(ISP)%ATOMTYPE=STRING
        CALL REPORT$CHVAL(NFIL,'ATOM TYPE',TRIM(STRING))
        CALL REPORT$CHVAL(NFIL,'  FUNCTIONAL TYPE',THIS1%FUNCTIONALID)
        CALL REPORT$R8VAL(NFIL,'  EXTENT OF LOCAL ORBITALS',THIS1%RCUT,'A_0')
        CALL REPORT$I4VAL(NFIL,'  MAX. ANG.MOMENTUM OF THE DENSITY',THIS1%LRX,' ')
        DO L=0,SIZE(THIS1%NORB)-1
          IF(THIS1%NORB(L+1).EQ.0) CYCLE
          IF(L.EQ.0) THEN
            CALL REPORT$I4VAL(NFIL,'  NUMBER OF CORRELATED S-SHELLS',THIS1%NORB(L+1),' ')
          ELSE IF(L.EQ.1) THEN
            CALL REPORT$I4VAL(NFIL,'  NUMBER OF CORRELATED P-SHELLS',THIS1%NORB(L+1),' ')
          ELSE IF(L.EQ.2) THEN
            CALL REPORT$I4VAL(NFIL,'  NUMBER OF CORRELATED D-SHELLS',THIS1%NORB(L+1),' ')
          ELSE IF(L.EQ.3) THEN
            CALL REPORT$I4VAL(NFIL,'  NUMBER OF CORRELATED F-SHELLS',THIS1%NORB(L+1),' ')
          ELSE 
            WRITE(STRING,FMT='("  NUMBER OF CORRELATED SHELLS WITH L=",I2)')L+1
            CALL REPORT$I4VAL(NFIL,STRING,THIS1%NORB(L+1),' ')
          END IF
        ENDDO
        IF(THIS1%FUNCTIONALID.EQ.'HYBRID') THEN
          CALL REPORT$R8VAL(NFIL,'  HARTREE-FOCK CONTRIBUTION' &
     &                           ,THIS1%HFWEIGHT*100,'PERCENT')
        ELSE IF(THIS1%FUNCTIONALID.EQ.'LDA+U') THEN
          CALL CONSTANTS('EV',EV)
          CALL REPORT$R8VAL(NFIL,'  U-PARAMETER',THIS1%UPAR/EV,'EV')
          CALL REPORT$R8VAL(NFIL,'  J-PARAMETER',THIS1%JPAR/EV,'EV')
          CALL REPORT$R8VAL(NFIL,'  F4/F2',THIS1%FRATIO42,'')
          CALL REPORT$R8VAL(NFIL,'  F6/F2',THIS1%FRATIO62,'')
          CALL REPORT$CHVAL(NFIL,'  DOUBLE COUNTING TYPE',DCTYPE)
        END IF
      ENDDO
      RETURN
      END SUBROUTINE LDAPLUSU$REPORT
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU$ETOT(ISP_,LMNXX,NDIMD,DENMAT,ETOT,DATH_)
!     **                                                                      **
!     **  THIS IS THE MAIN DRIVER ROUTINE FOR THE LDA+U CORRECTION            **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: ISP_
      INTEGER(4),INTENT(IN)  :: LMNXX
      INTEGER(4),INTENT(IN)  :: NDIMD
      COMPLEX(8),INTENT(IN)  :: DENMAT(LMNXX,LMNXX,NDIMD)
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
      INTEGER(4)             :: IS1,IS2,I,LN,M
      REAL(8)                :: SVAR,ETOT1
INTEGER(4)             :: LN1,LN2,LN3,LN4
!     **************************************************************************
      ETOT=0.D0
      DATH_=0.D0
      IF(.NOT.TON) RETURN
      CALL LDAPLUSU$SELECT(ISP_)
      IF(.NOT.THIS%TON) RETURN
                            CALL TRACE$PUSH('LDAPLUSU$ETOT')
!      
!     ==========================================================================
!     ==  CONSTRUCT LOCAL ORBITALS                                            ==
!     ==========================================================================
      IF(.NOT.THIS%TINI) THEN
        CALL SETUP$LMRX(ISP_,LMRX)
        THIS%LRX=INT(SQRT(REAL(LMRX)+1.D-8))-1
        LRX=THIS%LRX
!
        CALL LDAPLUSU_CHIFROMPHI()
PRINT*,'TINI ',THIS%TINI
PRINT*,'TON ',THIS%TON
PRINT*,'LNXCHI ',THIS%LNXCHI
PRINT*,'NR     ',THIS%NR
PRINT*,'LOXCHI ',THIS%LOXCHI
PRINT*,'NORB   ',THIS%NORB
PRINT*,'NCHI   ',THIS%NCHI
PRINT*,'DIEL   ',THIS%DIEL
PRINT*,'UPAR   ',THIS%UPAR
PRINT*,'JPAR   ',THIS%JPAR
!       ========================================================================
!       ==  CALCULATE SMALL U-TENSOR                                          ==
!       ========================================================================
        GID=THIS%GID
        NR=THIS%NR
        LNX=THIS%LNXCHI
        ALLOCATE(LOX(LNX))
        LOX=THIS%LOXCHI
        NCHI=THIS%NCHI
        ALLOCATE(THIS%ULITTLE(LRX+1,LNX,LNX,LNX,LNX))
        CALL LDAPLUSU_ULITTLE(GID,NR,LRX,LNX,LOX,THIS%CHI,THIS%ULITTLE)
        IF(THIS%USEDIEL.AND.(THIS%USEUPAR.OR.THIS%USEJPAR)) THEN
          CALL ERROR$MSG('DIEL AND (UPAR.OR.JPAR) ARE INCOMPATIBLE')
          CALL ERROR$L4VAL('USEDIEL',THIS%USEDIEL)
          CALL ERROR$L4VAL('USEUPAR',THIS%USEUPAR)
          CALL ERROR$L4VAL('USEJPAR',THIS%USEJPAR)
          CALL ERROR$STOP('LDAPLUSU$ETOT')
        END IF
        IF(THIS%USEDIEL) THEN
          THIS%ULITTLE=THIS%ULITTLE/THIS%DIEL
        ELSE IF(THIS%USEUPAR.OR.THIS%USEJPAR) THEN
          CALL LDAPLUSU_MODULITTLEWITHPARMS(LNX,LOX,LRX &
     &         ,THIS%USEUPAR,THIS%UPAR,THIS%USEJPAR,THIS%JPAR &
     &         ,THIS%USEFRATIO42,THIS%FRATIO42,THIS%USEFRATIO62,THIS%FRATIO62 &
     &         ,THIS%MAINLN,THIS%ULITTLE)
        END IF
      ELSE
        GID=THIS%GID
        NR=THIS%NR
        LRX=THIS%LRX
        LNX=THIS%LNXCHI
        ALLOCATE(LOX(LNX))
        LOX=THIS%LOXCHI
        NCHI=THIS%NCHI
      END IF
      CALL SETUP$SELECT(ISP)
      CALL SETUP$LNX(ISP,LNXPHI)
      ALLOCATE(LOXPHI(LNXPHI))
      CALL SETUP$LOFLN(ISP,LNXPHI,LOXPHI)
      NPHI=SUM(2*LOXPHI+1)
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
!     == TRANSFORM FROM TOTAL/SPIN TO UP/DOWN REPRESENTATION ===================
!
      CALL LDAPLUSU_SPINDENMAT('FORWARD',NDIMD,NPHI,DENMAT(1:NPHI,1:NPHI,:),MATSS)
!
!     == TRANSFORM FROM PARTIAL WAVES PHI TO LOCAL ORBITALS CHI ================
      CALL LDAPLUSU_MAPTOCHI(LNX,LOX,NCHI,LNXPHI,LOXPHI,NPHI,PHITOCHI)
      DO IS1=1,2
        DO IS2=1,2
          RHO(:,:,IS1,IS2)=MATMUL(PHITOCHI,MATMUL(MATSS(:,:,IS1,IS2),TRANSPOSE(PHITOCHI)))
        ENDDO
      ENDDO
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
  WRITE(*,FMT='(I3,100F8.3)')NCHI,PHITOCHI(I,:)
ENDDO
DO IS1=1,2
  DO IS2=1,2
    IF(SUM(ABS(RHO(:,:,IS1,IS2))).LT.1.D-3) CYCLE
    PRINT*,'===================== RHO FOR SPIN',IS1,IS2,' ======================'
    I=0
    DO LN=1,LNX
      DO M=1,2*LOX(LN)+1
        I=I+1
        WRITE(*,FMT='(I3,100F8.3)')LOX(LN),REAL(RHO(I,:,IS1,IS2))
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
DO LN1=1,LNX
  DO LN2=1,LNX
    DO LN3=1,LNX
      DO LN4=1,LNX
        WRITE(*,*)'ULITTLE',LN1,LN2,LN3,LN4,THIS%ULITTLE(:,LN1,LN2,LN3,LN4)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!
!     ==========================================================================
!     ==  CALCULATE U-TENSOR                                                  ==
!     ==========================================================================
      CALL LDAPLUSU_UTENSOR(LRX,NCHI,LNX,LOX,THIS%ULITTLE,U)
!
!     ==========================================================================
!     ==  HARTREE FOCK INTERACTION ENERGY                                     ==
!     ==========================================================================
      CALL LDAPLUSU_INTERACTION(NCHI,U,RHO,ETOT,HAM)
PRINT*,'E(U) ',ETOT
!
!     ==========================================================================
!     ==  DOUBLE COUNTING CORRECTION                                          ==
!     ==========================================================================
      IF(THIS%FUNCTIONALID.EQ.'LDA+U') THEN
        ALLOCATE(HAM1(NCHI,NCHI,2,2))
        CALL LDAPLUSU_DCLDAPLUSU(DCTYPE,LNX,LOX,NCHI,U,RHO,ETOT1,HAM1)
        ETOT=ETOT-ETOT1
        HAM=HAM-HAM1
        DEALLOCATE(HAM1)
      ELSE IF(THIS%FUNCTIONALID.EQ.'HYBRID') THEN
        ALLOCATE(HAM1(NCHI,NCHI,2,2))
        CALL LDAPLUSU_EDFT(GID,NR,NCHI,LNX,LOX,THIS%CHI,LRX,RHO,ETOT1,HAM1)
PRINT*,'E(DC) ',ETOT1
        ETOT=ETOT-ETOT1
        HAM=HAM-HAM1
        DEALLOCATE(HAM1)
!       == SCALE CORRECTION WITH 0.25 ACCORDING TO PBE0
        ETOT=ETOT*THIS%HFWEIGHT
        HAM=HAM*THIS%HFWEIGHT
      ELSE
        CALL ERROR$MSG('FUNCTIONALID NOT RECOGNIZED')
        CALL ERROR$CHVAL('FUNCTIONALID',THIS%FUNCTIONALID)
        CALL ERROR$STOP('LDAPLUSU$ETOT')
      END IF
print*,'marke 1'
!
!     ==========================================================================
!     ==  UPFOLD                                                              ==
!     ==========================================================================
!     == TRANSFORM FROM CHI TO PHI =============================================
      DO IS1=1,2
        DO IS2=1,2
          MATSS(:,:,IS1,IS2)=MATMUL(TRANSPOSE(PHITOCHI),MATMUL(HAM(:,:,IS1,IS2),PHITOCHI))
        ENDDO
      ENDDO
!
!     == TRANSFORM FROM (SPIN,SPIN) TO (TOTAL,SPIN) ============================
      CALL LDAPLUSU_SPINDENMAT('BACK',NDIMD,LMNXX,DATH,MATSS)
!
!     == MAKE REAL (THIS IS A FUDGE TO BE FIXED IN AUGMENTATION!)
      DATH_(:,:,:)=(0.D0,0.D0)
      DATH_(:NPHI,:NPHI,:)=REAL(DATH)
!!$DO IS1=1,NDIMD
!!$PRINT*,'===================== DATH FOR SPIN',IS1,IS2,' ======================'
!!$I=0
!!$DO LN=1,LNXPHI
!!$DO M=1,2*LOXPHI(LN)+1
!!$I=I+1
!!$WRITE(*,FMT='(I3,100F8.3)')LOXPHI(LN),REAL(DATH(I,:,IS1))
!!$ENDDO
!!$ENDDO
!!$ENDDO
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
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_CHIFROMPHI()
!     **                                                                      **
!     **  DEFINES TRANSFORMATION FROM PARTIAL WAVES TO LOCAL ORBITALS.        **
!     **  THE RESULT IS STORED WITHOUT MAGNETIC QUANTUM NUMBERS.              **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      REAL(8)               :: RCUT
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: LNX
      INTEGER(4),ALLOCATABLE:: LOX(:)
      REAL(8)   ,ALLOCATABLE:: PHI(:,:)
      INTEGER(4)            :: LNXCHI 
      INTEGER(4),ALLOCATABLE:: LOXCHI(:)
      REAL(8)   ,ALLOCATABLE:: CHI(:,:)
      REAL(8)   ,ALLOCATABLE:: CHI1(:,:)
      REAL(8)   ,ALLOCATABLE:: A(:,:)
      REAL(8)   ,ALLOCATABLE:: A1(:,:)
      REAL(8)   ,ALLOCATABLE:: MAT(:,:)
      REAL(8)   ,ALLOCATABLE:: MATINV(:,:)
      REAL(8)   ,ALLOCATABLE:: R(:)        
      REAL(8)   ,ALLOCATABLE:: G(:)        
      REAL(8)   ,PARAMETER  :: RCG=0.3D0
      REAL(8)   ,ALLOCATABLE:: AUX(:)
      REAL(8)               :: SVAR1,SVAR2
      INTEGER(4)            :: IR
      INTEGER(4)            :: NX,N,LX,L,LN,LNCHI,LN0,NOFL,ISVAR
      INTEGER(4)            :: N1,N2,LN1,LN2,L1,L2
!     **************************************************************************
                            CALL TRACE$PUSH('LDAPLUSU_CHIFROMPHI')
      CALL SETUP$ISELECT(ISP)
      CALL SETUP$GETI4('GID',GID)
      THIS%GID=GID
      CALL RADIAL$GETI4(GID,'NR',NR)
      THIS%NR=NR
      CALL SETUP$LNX(ISP,LNX)
      ALLOCATE(LOX(LNX))
      CALL SETUP$LOFLN(ISP,LNX,LOX)
      ALLOCATE(PHI(NR,LNX))
      CALL SETUP$AEPARTIALWAVES(ISP,NR,LNX,PHI)
!
!     ==========================================================================
!     ==  DIVIDE IN HEAD AN TAIL FUNCTIONS                                    ==
!     ==========================================================================
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      ALLOCATE(AUX(NR))
!     == COUNT NUMBER OF PARTIAL WAVES PER L AS AS ESTIMATE FOR THE NUMBER OF CHI
!     == COUNT ONLY THOSE ANGULAR MOMENTUM CHANNELS THAT HAVE CORRELATED ORBITALS
      LX=MAXVAL(LOX)
      LNXCHI=0
PRINT*,'CHIFROMPHI: ',THIS%NORB(1:LX+1)
      DO L=0,LX
        IF(THIS%NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
        NOFL=0
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE !CONSIDER ONLY PARTIAL WAVES WITH THE CORRECT L
          NOFL=NOFL+1
          LNXCHI=LNXCHI+1
        ENDDO
        IF(THIS%NORB(L+1).GT.NOFL) THEN
          CALL ERROR$MSG('#(CORRELATED ORBITALS) EXCEEDS #(PARTIAL WAVES)')
          CALL ERROR$STOP('LDAPLUSU_CHIFROMPHI')
        END IF
      ENDDO
PRINT*,'CHIFROMPHI: ',L,NOFL,LNXCHI
!
!     == ORDER ACCORDING TO L ==================================================
      ALLOCATE(LOXCHI(LNXCHI))        
      ALLOCATE(CHI(NR,LNXCHI))        
!     == |CHI_I>=SUM_I |PHI_J>A(J,I)
      ALLOCATE(A(LNX,LNXCHI))        
      A(:,:)=0.D0
      LNCHI=0
      DO L=0,LX
        IF(THIS%NORB(L+1).EQ.0) CYCLE !IGNORE CHANNELS WITHOUT CORRELATED ORBITALS
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE !CONSIDER ONLY PARTIAL WAVES WITH THE CORRECT L
          LNCHI=LNCHI+1
          LOXCHI(LNCHI)=LOX(LN)
          CHI(:,LNCHI)=PHI(:,LN)
          A(LN,LNCHI)=1.D0
        ENDDO
      ENDDO
PRINT*,'CHIFROMPHI: LNXPHI',LNX
PRINT*,'CHIFROMPHI: LOXPHI',LOX
PRINT*,'CHIFROMPHI: LNXCHI',LNXCHI
PRINT*,'CHIFROMPHI: LOXCHI',LOXCHI
!
!     ==========================================================================
!     == MAKE HEAD FUNCTION ANTIBONDING WITH NODE AT RCUT ======================
!     ==========================================================================
      RCUT=THIS%RCUT
      L=-1
      DO LN=1,LNXCHI
        IF(LOXCHI(LN).NE.L) NOFL=0    ! RESET ORBITAL COUNTER FOR EACH L
        L=LOXCHI(LN)
        NOFL=NOFL+1
        IF(NOFL.GT.THIS%NORB(L+1)) CYCLE  !ONLY CORRELATED ORBITALS WILL BE LOCALIZED
        IF(LN+1.GT.LNXCHI) CYCLE          !CHECK IF THERE IS A TAIL FUNCTION LEFT
        IF(LOXCHI(LN+1).NE.L) CYCLE       !CHECK IF THERE IS A TAIL FUNCTION LEFT
!       == IMPOSE NODE CONDITION================================================
        CALL RADIAL$VALUE(GID,NR,CHI(:,LN),RCUT,SVAR1)
        CALL RADIAL$VALUE(GID,NR,CHI(:,LN+1),RCUT,SVAR2)
        IF(SVAR1.EQ.0.D0.AND.SVAR2.EQ.0.D0) THEN
          CALL ERROR$MSG('PARTIAL WAVES ARE TRUNCATED INSIDE OF RCUT')
          CALL ERROR$MSG('THIS IS A FLAW OF THE IMPLEMENTATION')
          CALL ERROR$MSG('CHOOSE SMALLER RCUT')
          CALL ERROR$STOP('LDAPLUSU_CHIFROMPHI')
        END IF
        CHI(:,LN)=CHI(:,LN)*SVAR2-CHI(:,LN+1)*SVAR1
        A(:,LN)=A(:,LN)*SVAR2-A(:,LN+1)*SVAR1
PRINT*,'A ',LN,A(:,LN)
!       == CUT AT RCUT           ===============================================
        DO IR=1,NR
          IF(R(IR).GT.RCUT) CHI(IR,LN)=0.D0
        ENDDO
!       == ORTHOGONALIZE TO THE LOWER HEAD FUNCTIONS ===========================
PRINT*,'CHIFROMPHI: LN',LN,LN-NOFL+1,LN-1
        DO LN1=LN-NOFL+1,LN-1
          AUX(:)=CHI(:,LN)*CHI(:,LN1)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          CHI(:,LN)=CHI(:,LN)-CHI(:,LN1)*SVAR1
          A(:,LN)=A(:,LN)-A(:,LN1)*SVAR1
        ENDDO
!       == NORMALIZE HEAD FUNCTION =============================================
        AUX(:)=CHI(:,LN)**2*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
        SVAR1=1.D0/SQRT(SVAR1)
        CHI(:,LN)=CHI(:,LN)*SVAR1
        A(:,LN)=A(:,LN)*SVAR1
      END DO          
!
!     ==========================================================================
!     == DELOCALIZE TAIL FUNCTIONS                                            ==
!     ==========================================================================
      ALLOCATE(CHI1(NR,LNXCHI))
      CHI1(:,:)=CHI(:,:)
      ALLOCATE(A1(LNX,LNXCHI))        
      A1(:,:)=A(:,:)
      ALLOCATE(G(NR))
      G(:)=EXP(-(R(:)/RCG)**2)
      L=-1
      DO LN=1,LNXCHI
        IF(LOXCHI(LN).NE.L) THEN
          L=LOXCHI(LN)
          LN0=LN
          CYCLE
        END IF  
!       == MINIMIZE CONTRIBUTION NEAR THE CENTER 
        DO LN2=LN0,LN-1
          AUX(:)=G(:)*CHI1(:,LN)*CHI1(:,LN2)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          AUX(:)=G(:)*CHI1(:,LN2)**2*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
          SVAR1=SVAR1/SVAR2
          CHI1(:,LN)=CHI1(:,LN)-CHI1(:,LN2)*SVAR1
          A1(:,LN)=A1(:,LN)-A1(:,LN2)*SVAR1
        ENDDO
      ENDDO
!     = NOW MAP ONLY TAIL FUNCTIONS BACK AND LEAVE HEAD FUNCTIONS UNTOUCHED        
      L=-1
      DO LN=1,LNXCHI
        IF(LOXCHI(LN).NE.L) THEN
          L=LOXCHI(LN)
          NOFL=0
        END IF  
        NOFL=NOFL+1
        IF(NOFL.LE.THIS%NORB(L+1)) CYCLE ! LEAVE HEAD FUNCTIONS ALONE
        CHI(:,LN)=CHI1(:,LN)
        A(:,LN)=A1(:,LN)
      ENDDO
!
!     === CONSTRUCT TRANSFORMATION MATRIX FROM A ===============================
      LX=MAXVAL(LOXCHI)
      DO L=1,LX
!
        NX=0
        DO LN=1,LNXCHI
          IF(LOXCHI(LN).EQ.L) NX=NX+1
        ENDDO
        IF(NX.EQ.0) CYCLE   

        ALLOCATE(MAT(NX,NX))
        ALLOCATE(MATINV(NX,NX))
!        
        N1=0
        DO LN1=1,LNXCHI
          L1=LOXCHI(LN1)
          IF(L1.NE.L) CYCLE
          N1=N1+1
          N2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            IF(L2.NE.L) CYCLE
            N2=N2+1
            MAT(N2,N1)=A(LN2,LN1)
          ENDDO
        ENDDO
        CALL LIB$INVERTR8(NX,MAT,MATINV)
        N1=0
        DO LN1=1,LNXCHI
          L1=LOXCHI(LN1)
          IF(L1.NE.L) CYCLE
          N1=N1+1
          N2=0
          DO LN2=1,LNX
            L2=LOX(LN2)
            IF(L2.NE.L) CYCLE
            N2=N2+1
            A(LN2,LN1)=MATINV(N1,N2) ! A IT TRANSPOSED SO THAT THE INDICES MATCH
          ENDDO
        ENDDO
        DEALLOCATE(MAT)
        DEALLOCATE(MATINV)
      ENDDO
!
!     === REMOVE TAIL FUNCTIONS ================================================
      LN1=0
      L=-1
      DO LN=1,LNXCHI
        IF(L.NE.LOXCHI(LN)) THEN
          L=LOXCHI(LN)
          N=0
          NX=THIS%NORB(L+1)
        END IF
        N=N+1
        IF(N.LE.NX) THEN
          LN1=LN1+1
          LOXCHI(LN1)=LOXCHI(LN)
          CHI(:,LN1)=CHI(:,LN)
!         -- A HAS BEEN INVERTED. THEREFORE THE CHI-INDEX IS LEFT 
!         -- AND THE PHI-INDEX IN ON THE RIGHT HAND SIDE
          A(:,LN1)=A(:,LN)  
        END IF
      ENDDO
      LNXCHI=LN1

!PRINT*,'== WRITE CHI.DAT'
!OPEN(UNIT=109,FILE='CHI.DAT',FORM='FORMATTED')
!DO IR=1,NR
!  WRITE(109,*)R(IR),CHI(IR,1:LNXCHI)
!ENDDO
!CLOSE(109)
!STOP 'FORCED AFTER WRITING CHI'
!
!     ==========================================================================
!     ==  CUT OF SUPPORT FUNCTIONS                                            ==
!     ==========================================================================
      THIS%LNXCHI=LNXCHI
      ALLOCATE(THIS%LOXCHI(LNXCHI))
      THIS%LOXCHI=LOXCHI(1:LNXCHI)
      THIS%NCHI=SUM(2*LOXCHI(1:LNXCHI)+1)
      ALLOCATE(THIS%CHI(NR,LNXCHI))
      THIS%CHI=CHI(:,1:LNXCHI)
      ALLOCATE(THIS%DOWNFOLD(LNXCHI,LNX))
      DO LN=1,LNXCHI
         THIS%DOWNFOLD(LN,:)=A(:,LN)
      ENDDO
!!$OPEN(10,FILE='CHI_'//TRIM(THIS%ATOMTYPE)//'.DAT')
!!$REWIND 10
!!$DO IR=1,NR
!!$WRITE(10,*)R(IR),CHI(IR,:)
!!$ENDDO
!!$CLOSE(10)
!!$STOP
!
!     ==========================================================================
!     ==  CLEAN UP                                                            ==
!     ==========================================================================
      DEALLOCATE(CHI)
      DEALLOCATE(R)
                            CALL TRACE$POP()
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
      PI=4.D0*DATAN(1.D0)
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
      PI=4.D0*DATAN(1.D0)
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
            ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)=ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
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
            ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)=ULITTLE(2:,LNPROBE,LNPROBE,LNPROBE,LNPROBE)*SVAR
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
                          SVAR=SVAR+FOURPIBY2LPLUS1*CG1*CG2*ULITTLE(L+1,LN2,LN4,LN3,LN1)
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
        upar=upar/real((2*l+1)**2,kind=8)
        jpar=upar-jpar/real(2*l*(2*l+1),kind=8)
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
            HAM(I,I,IS,IS)=HAM(I,I,IS,IS)+CMPLX(V(IS),0.D0)
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
!     **                                                                      **
!     ** CALCULATES THE HARTREE AND EXCHANGE ENERGY FROM THE U-TENSOR U       **
!     ** AND THE DENSITY MATRIX RHO.                                          **
!     **                                                                      **
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
      SUBROUTINE LDAPLUSU_EDFT(GID,NR,LMNX,LNX,LOX,CHI,LRX,DENMAT,ETOT,HAM)
!     **                                                                      **
!     ** DOUBLE COUNTING CORRECTION FOR THE HYBRID FUNCTIONAL                 **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      INTEGER(4)  ,INTENT(IN) :: NR
      INTEGER(4)  ,INTENT(IN) :: LRX
      INTEGER(4)  ,INTENT(IN) :: LMNX       ! #(LOCAL ORBITALS)
      INTEGER(4)  ,INTENT(IN) :: LNX        ! #(RADIAL FUNCTIONS)
      INTEGER(4)  ,INTENT(IN) :: LOX(LNX)   ! MAIN ANGULAR MOMENTUM OF LOCAL ORB.
      REAL(8)     ,INTENT(IN) :: CHI(NR,LNX)
      COMPLEX(8)  ,INTENT(IN) :: DENMAT(LMNX,LMNX,2,2) ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT       ! DOUBLE COUNTINNG ENERGY
      COMPLEX(8)  ,INTENT(OUT):: HAM(LMNX,LMNX,2,2)  ! DETOT/D(RHO^*)        
      COMPLEX(8)              :: DENMAT1(LMNX,LMNX,4)
      COMPLEX(8)              :: HAM1(LMNX,LMNX,4)
      REAL(8)                 :: R(NR)
      REAL(8)     ,ALLOCATABLE:: RHO(:,:,:)
      REAL(8)     ,ALLOCATABLE:: POT(:,:,:)
      REAL(8)                 :: EDENSITY(NR)
      REAL(8)                 :: AUX(NR),SVAR
      INTEGER(4)              :: LMRX,L
      INTEGER(4)              :: IDIM,LM
      COMPLEX(8)  ,PARAMETER  :: CI=(0.D0,1.D0)
      INTEGER(4)  ,PARAMETER  :: NDIMD=4
integer(4) :: lmrx1,ir
integer(4) :: imethod
 REAL(8)     ,ALLOCATABLE:: RHOtest(:,:,:)
 REAL(8)     ,ALLOCATABLE:: POTtest(:,:,:)
 REAL(8)     ,ALLOCATABLE:: RHOtest2(:,:,:)
 REAL(8)     ,ALLOCATABLE:: POTtest2(:,:,:)
 REAL(8)                 :: etot2
!     **************************************************************************
      LMRX=(LRX+1)**2
      ETOT=0.D0
!
!     ==========================================================================
!     ==  TRANSFORM DENSITY MATRIX FROM UP/DOWN TO TOTAL/SPIN                 ==
!     ==========================================================================
      DENMAT1(:,:,1)=DENMAT(:,:,1,1)+DENMAT(:,:,2,2)
      DENMAT1(:,:,2)=DENMAT(:,:,1,2)+DENMAT(:,:,2,1)
      DENMAT1(:,:,3)=-CI*(DENMAT(:,:,1,2)-DENMAT(:,:,2,1))
      DENMAT1(:,:,4)=DENMAT(:,:,1,1)-DENMAT(:,:,2,2)
!
!     ==========================================================================
!     ==  CALCULATE DENSITY                                                   ==
!     ==========================================================================
      ALLOCATE(RHO(NR,LMRX,ndimd))
      DO IDIM=1,NDIMD
        CALL AUGMENTATION_RHO(NR,LNX,LOX,CHI &
     &                       ,LMNX,DENMAT1(:,:,IDIM),LMRX,RHO(:,:,IDIM))
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE ENERGY AND POTENTIAL                                      ==
!     ==========================================================================
      ALLOCATE(POT(NR,LMRX,NDIMD))
      POT(:,:,:)=0.D0
!     == EXCHANGE ENERGY AND POTENTIAL =========================================
      CALL DFT$SETL4('XCONLY',.TRUE.)
!
!===============================================================================
!==== dangerous code!!!!                                                    ====
!==== this formulation is based on a noncollinear formulation, which        ====
!==== yields different results from a collinear formulation even for        ====
!==== a collinear density                                                   ====
!====                                                                       ====
!==== different versions are implemented                                    ====
!==== 1) default noncollinear method                                        ====
!==== 2) a collinear method (which is compared on the fly with the          ====
!====    noncollinear method                                                ====
!===============================================================================
      CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,ETOT,POT)
imethod=0
!imethod=1
      if(imethod.eq.1) then
!       == collinear method with collinear density
        ALLOCATE(RHOtest(NR,LMRX,2))
        ALLOCATE(pottest(NR,LMRX,2))
        pottest(:,:,1)=0.d0
        rhotest(:,:,1)=rho(:,:,1)
        rhotest(:,:,2)=rho(:,:,4)
        CALL AUGMENTATION_XC(GID,NR,LMRX,2,RHOtest,ETOT,POTtest)
        pot(:,:,:)=0.d0
        pot(:,:,1)=pottest(:,:,1)
        pot(:,:,4)=pottest(:,:,2)
        deallocate(rhotest)
        deallocate(pottest)
!
!      else if(imethod.eq.2) then
!       == noncollinear method with collinear density ==========================
        ALLOCATE(RHOtest2(NR,LMRX,ndimd))
        ALLOCATE(pottest2(NR,LMRX,ndimd))
        rhotest2(:,:,:)=0.d0
        pottest2(:,:,:)=0.d0
        rhotest2(:,:,1)=rho(:,:,1)
        rhotest2(:,:,4)=rho(:,:,4)
        CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHOtest2,ETOT2,POTtest2)
!print*,'ldaplusutest',etot2-etot,maxval(abs(pottest2-pot)),maxloc(abs(pottest2-pot))
!        etot=etot2
!        pot(:,:,:)=pottest2(:,:,:)
        deallocate(rhotest2)
        deallocate(pottest2)
!
      else if(imethod.eq.3) then
!       == comparison ==========================================================

      end if
      CALL DFT$SETL4('XCONLY',.FALSE.)
PRINT*,'EXC ',ETOT
!
!     == HARTREE ENERGY AND POTENTIAL ==========================================
      EDENSITY=0.D0
      DO LM=1,LMRX
        L=INT(SQRT(REAL(LM-1,KIND=8))+1.D-5)
        CALL RADIAL$POISSON(GID,NR,L,RHO(:,LM,1),AUX)
        POT(:,LM,1)=POT(:,LM,1)+AUX(:)
        EDENSITY(:)=EDENSITY(:)+0.5D0*AUX(:)*RHO(:,LM,1)
      ENDDO
      CALL RADIAL$R(GID,NR,R)
      EDENSITY=EDENSITY*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,EDENSITY,SVAR)
PRINT*,'EH ',SVAR
      ETOT=ETOT+SVAR
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
      HAM(:,:,1,1)=HAM1(:,:,1)+HAM1(:,:,4)
      HAM(:,:,1,2)=HAM1(:,:,2)-CI*HAM1(:,:,3)
      HAM(:,:,2,1)=HAM1(:,:,2)+CI*HAM1(:,:,3)
      HAM(:,:,2,2)=HAM1(:,:,1)-HAM1(:,:,4)
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
                DATH(LMN1,LMN2,ISPIN)=DATH(LMN1,LMN2,ISPIN)+CMPLX(SVAR,0.D0)
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

