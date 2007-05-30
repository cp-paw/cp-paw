!TODO :
! DATH IS STILL REAL AND SHOULD PROBABLY BE COMPLEX LIKE DENMAT
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE LDAPLUSU_MODULE
TYPE THISTYPE
LOGICAL(4)             :: TINI=.FALSE.
LOGICAL(4)             :: TON=.FALSE.
character(8)           :: atomtype=''
INTEGER(4)             :: GID          !grid id for radial grid
INTEGER(4)             :: NR           !#(radial grid points)
INTEGER(4)             :: LNXCHI       !#(radial functions for local orbitals)
INTEGER(4),POINTER     :: LOXCHI(:)    !main angular momentum of local orbital
INTEGER(4),pointer     :: NORB(:)      !X(# LOCAL FUNCTIONS PER ANGULAR MOMENTUM)
INTEGER(4)             :: LRX          !x(angular momentum in the density)
REAL(8)                :: RCUT=0.D0    !RADIUS of local orbital
character(16)          :: functionalid !can be 'lda+u' or 'hybrid'
!== settings specifically for hybrid functional ================================
real(8)                :: hfweight=0.d0!contribution of exact exchange
!== settings specifically for lda+u
INTEGER(4)             :: MAINLN(2)=(/0,0/) !SHELL TO WHICH UPAR AND JPAR REFER TO
LOGICAL(4)             :: USEDIEL=.FALSE.
LOGICAL(4)             :: USEUPAR=.FALSE.
LOGICAL(4)             :: USEJPAR=.FALSE.
REAL(8)                :: DIEL=1            !DIELECTRIC CONSTANT
REAL(8)                :: UPAR=0.D0         !U-PARAMETER
REAL(8)                :: JPAR=0.D0         !J-PARAMETER
!==  intermediate storage from initualization ================================== 
INTEGER(4)             :: NCHI              !#(local orbitals)
REAL(8)   ,POINTER     :: CHI(:,:)          !local (HEAD) ORBITALS
REAL(8)   ,POINTER     :: ULITTLE(:,:,:,:,:)!SLATER INTEGRALS (EXCEPT FACTOR)
REAL(8)   ,POINTER     :: DOWNFOLD(:,:)     !MAPS PARTIAL WAVES TO local ORBITALS
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
!     **  attention: changes the setting of setup object                      **
!     **                                                                      **
      USE LDAPLUSU_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSP_    ! #(different atom types)
      character(8)          :: string
!     **************************************************************************
!
!     ==========================================================================
!     == do nothing if thisarray is already allocated                         ==
!     ==========================================================================
      IF(NSP.NE.0) THEN
        IF(NSP_.NE.NSP) THEN
          CALL ERROR$MSG('LDAPLUSU$NEW CANNOT BE CALLED TWICE')
          CALL ERROR$STOP('LDAPLUSU$NEW')
        END IF
        return
      END IF
!
!     ==========================================================================
!     == create thisarray                                                     ==
!     ==========================================================================
      NSP=NSP_
      ISP=0
      ALLOCATE(THISARRAY(NSP))
      DO ISP=1,NSP
        ALLOCATE(THISARRAY(ISP)%NORB(1))
        THISARRAY(ISP)%NORB(:)=0  ! PER DEFAULT CORRELATE nothing
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
      INTEGER(4),INTENT(IN) :: ISP_  ! unique index of the atom type
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
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NFIL
      TYPE(THISTYPE),POINTER :: THIS1
      CHARACTER(64)          :: STRING
      INTEGER(4)             :: L
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
      INTEGER(4)             :: Lmrx
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
      REAL(8)                :: SVAR,etot1
INTEGER(4)             :: ln1,ln2,ln3,ln4
!     **************************************************************************
      etot=0.d0
      dath_=0.d0
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
        this%lrx=INT(SQRT(REAL(LMrX)+1.D-8))-1
        lrx=this%lrx
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
          CALL LDAPLUSU_MODULITTLEWITHPARMS(LNX,LOX,LRX,THIS%USEUPAR,THIS%UPAR &
     &                           ,THIS%USEJPAR,THIS%JPAR,THIS%MAINLN,THIS%ULITTLE)
        END IF
      ELSE
        GID=THIS%GID
        NR=THIS%NR
        lrx=this%lrx
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
!     == TRANSFORM FROM TOTAL/SPIN TO up/down representation ===================
!
      CALL LDAPLUSU_SPINDENMAT('FORWARD',NDIMD,NPHI,DENMAT(1:NPHI,1:NPHI,:),MATSS)
!
!     == TRANSFORM FROM partial waves phi to local orbitals CHI ================
      CALL LDAPLUSU_MAPTOCHI(LNX,LOX,NCHI,LNXPHI,LOXPHI,NPHI,PHITOCHI)
      DO IS1=1,2
        DO IS2=1,2
          RHO(:,:,IS1,IS2)=MATMUL(PHITOCHI,MATMUL(MATSS(:,:,IS1,IS2),TRANSPOSE(PHITOCHI)))
        ENDDO
      ENDDO
!
! printout for testing==========================================================
DO IS1=1,ndimd
  PRINT*,'===================== denmat FOR SPIN',IS1,' ======================'
  I=0
  DO LN=1,LNXphi
    DO M=1,2*LOXphi(LN)+1
      I=I+1
      WRITE(*,FMT='(I3,100F8.3)')LOXphi(LN),REAL(denmat(I,:,IS1))
    ENDDO
  ENDDO
ENDDO
PRINT*,'===================== phitochi ======================'
do i=1,nchi
  WRITE(*,FMT='(I3,100F8.3)')nchi,phitochi(I,:)
enddo
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
do is1=1,2
  svar=0.d0
  do ln=1,nchi
    svar=svar+real(rho(ln,ln,is1,is1))
  enddo
  print*,'charge= ',svar,' for spin ',is1
enddo
!
do ln1=1,lnx
  do ln2=1,lnx
    do ln3=1,lnx
      do ln4=1,lnx
        write(*,*)'ulittle',ln1,ln2,ln3,ln4,this%ulittle(:,ln1,ln2,ln3,ln4)
      enddo
    enddo
  enddo
enddo
!
!     ==========================================================================
!     ==  CALCULATE U-TENSOR                                                  ==
!     ==========================================================================
      CALL LDAPLUSU_UTENSOR(LRX,NCHI,LNX,LOX,THIS%ULITTLE,U)
!
!     ==========================================================================
!     ==  Hartree Fock interaction ENERGY                                     ==
!     ==========================================================================
      CALL LDAPLUSU_INTERACTION(Nchi,U,RHO,ETOT,HAM)
print*,'e(u) ',etot
!
!     ==========================================================================
!     ==  double counting correction                                          ==
!     ==========================================================================
      IF(THIS%FUNCTIONALID.EQ.'LDA+U') THEN
        ALLOCATE(HAM1(NCHI,NCHI,2,2))
        CALL LDAPLUSU_DCLDAPLUSU(DCTYPE,lnx,lox,NCHI,U,RHO,ETOT1,HAM1)
        ETOT=ETOT-ETOT1
        HAM=HAM-HAM1
        DEALLOCATE(HAM1)
      ELSE IF(THIS%FUNCTIONALID.EQ.'HYBRID') THEN
        ALLOCATE(HAM1(NCHI,NCHI,2,2))
        call LDAPLUSU_edft(gid,nr,nchi,lnx,lox,this%chi,lrx,rho,etot1,ham1)
print*,'e(dc) ',etot1
        ETOT=ETOT-ETOT1
        HAM=HAM-HAM1
        DEALLOCATE(HAM1)
!       == scale correction with 0.25 according to PBE0
        etot=etot*this%hfweight
        ham=ham*this%hfweight
      ELSE
        CALL ERROR$MSG('FUNCTIONALID NOT RECOGNIZED')
        CALL ERROR$CHVAL('FUNCTIONALID',THIS%FUNCTIONALID)
        CALL ERROR$STOP('LDAPLUSU$ETOT')
      END IF
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
!     **  defines transformation from partial waves to local orbitals.        **
!     **  the result is stored without magnetic quantum numbers.              **
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
print*,'chifromphi: ',this%norb(1:lx+1)
      DO L=0,LX
        IF(THIS%NORB(L+1).EQ.0) CYCLE !ignore CHANNELS WITHOUT CORRELATED ORBITALS
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
print*,'chifromphi: ',l,nofl,lnxchi
!
!     == ORDER ACCORDING TO L ==================================================
      ALLOCATE(LOXCHI(LNXCHI))        
      ALLOCATE(CHI(NR,LNXCHI))        
!     == |CHI_I>=SUM_I |PHI_J>A(J,I)
      ALLOCATE(A(LNX,LNXCHI))        
      A(:,:)=0.D0
      LNCHI=0
      DO L=0,LX
        IF(THIS%NORB(L+1).EQ.0) CYCLE !ignore CHANNELS WITHOUT CORRELATED ORBITALS
        DO LN=1,LNX
          IF(LOX(LN).NE.L) CYCLE !CONSIDER ONLY PARTIAL WAVES WITH THE CORRECT L
          LNCHI=LNCHI+1
          LOXCHI(LNCHI)=LOX(LN)
          CHI(:,LNCHI)=PHI(:,LN)
          A(LN,LNCHI)=1.D0
        ENDDO
      ENDDO
print*,'chifromphi: lnxphi',lnx
print*,'chifromphi: loxphi',lox
print*,'chifromphi: lnxchi',lnxchi
print*,'chifromphi: loxchi',loxchi
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
print*,'a ',ln,a(:,ln)
!       == CUT AT RCUT           ===============================================
        DO IR=1,NR
          IF(R(IR).GT.RCUT) CHI(IR,LN)=0.D0
        ENDDO
!       == ORTHOGONALIZE TO THE LOWER HEAD FUNCTIONS ===========================
print*,'chifromphi: ln',ln,ln-nofl+1,ln-1
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
            A(LN2,LN1)=MATINV(N1,N2) ! a it transposed so that the indices match
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
      do ln=1,lnxchi
         THIS%DOWNFOLD(ln,:)=A(:,ln)
      enddo
!!$OPEN(10,FILE='CHI_'//TRIM(THIS%ATOMTYPE)//'.DAT')
!!$rewind 10
!!$DO Ir=1,NR
!!$WRITE(10,*)R(Ir),CHI(Ir,:)
!!$eNDDO
!!$CLOSE(10)
!!$stop
!
!     ==========================================================================
!     ==  CLEAN UP                                                            ==
!     ==========================================================================
      DEALLOCATE(chi)
      DEALLOCATE(R)
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_SPINDENMAT(ID,NDIMD,LMNXX,MAT1,MAT2)
!     **                                                                      **
!     ** IF="FORWARD": CONVERTS DENSITY MATRIX MAT1 from                      **
!     **   (TOTAL,SPIN) REPRESENTATION INTO (SPIN,SPIN) REPRESENTATION MAT2   **
!     ** IF="BACK": CONVERTS HAMILTON MATRIX MAT2 from                        **
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
          CALL ERROR$i4val('NDIMD',ndimd)
          CALL ERROR$stop('LDAPLUSU_SPINDENMAT')
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
          CALL ERROR$i4val('NDIMD',ndimd)
          CALL ERROR$stop('LDAPLUSU_SPINDENMAT')
        END IF
      ELSE
        CALL ERROR$MSG('ID MUST BE EITHER "FORWARD" OR "BACK"')
        CALL ERROR$stop('LDAPLUSU_SPINDENMAT')
      END IF
                            CALL TRACE$POP()
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LDAPLUSU_MAPTOCHI(LNXCHI,LOXCHI,LMNXCHI &
     &                            ,LNXPHI,LOXPHI,LMNXPHI,DOWNFOLD)
!     **                                                                      **
!     **  expands transformation from partial waves to local orbitals         **
!     **  to full size                                                        **
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
!     ** Slater integrals. (Attention: different pre-factors)                 **
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
!if(lox(ln1).ne.lox(ln2).or.lox(ln2).ne.lox(ln3).or.lox(ln3).ne.lox(ln4)) svar=0.d0
!if(lox(ln1)*lox(ln2)*lox(ln3)*lox(ln4).eq.0) svar=0.d0
                ULITTLE(L+1,LN1,LN2,LN3,LN4)=SVAR
                ULITTLE(L+1,LN2,LN1,LN3,LN4)=SVAR
                ULITTLE(L+1,LN1,LN2,LN4,LN3)=SVAR
                ULITTLE(L+1,LN2,LN1,LN4,LN3)=SVAR
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
      SUBROUTINE LDAPLUSU_MODULITTLEWITHPARMS(LNX,LOX,LRX,USEUPAR,UPAR &
     &                                       ,USEJPAR,JPAR,MAINLN,ULITTLE)
!     **                                                                      **
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
        RAWUPAR=ULITTLE(1,LNPROBE,LNPROBE,LNPROBE,LNPROBE)/FOURPI
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
!     ** CALCULATES THE U-tensor                                              **
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
!     **************************************************************************
                            CALL TRACE$PUSH('LDAPLUSU_UTENSOR')
      U(:,:,:,:)=0.d0
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
                      if(maxval(abs(ulittle(:,ln2,ln4,ln3,ln1))).eq.0.d0) cycle
                      LX=MIN(LRX,L2+L4,L1+L3)
                      SVAR=0.D0
                      LM=0
                      DO L=0,LX
                        DO M=1,2*L+1
                          LM=LM+1
                          CALL CLEBSCH(LM2,LM4,LM,CG1)
                          CALL CLEBSCH(LM3,LM1,LM,CG2)
                          SVAR=SVAR+CG1*CG2*ULITTLE(L+1,LN2,LN4,LN3,LN1)
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
      SUBROUTINE LDAPLUSU_dcldaplusu(ID,lnx,lox,Nchi,U,RHO,ETOT,HAM)
!     **                                                                      **
!     **  calculates the correlation energy from utensor and density matrix   **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID     ! ID FOR DOUBLE COUNTING CORRECTION
      INTEGER(4)  ,INTENT(IN) :: Lnx    !                     
      INTEGER(4)  ,INTENT(IN) :: Lox(lnx) ! angular momentum  
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
      INTEGER(4)              :: IS,I,j,ln
      INTEGER(4)              :: nchi1,nchi2
      INTEGER(4)              :: l
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
            IF(I.NE.J)JPAR=JPAR+U(I,J,J,I)
          ENDDO
        ENDDO
        UPAR=UPAR/REAL(2*L+1)**2
        JPAR=JPAR/REAL(2*L*(2*L+1))*7.D0/5.D0
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
!     ** CALCULATES THE hartree and exchange energy from the u-tensor u       **
!     ** and the density matrix rho.                                          **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NORB                   ! bASIS-SET SIZE              
      REAL(8)     ,INTENT(IN) :: U(NORB,NORB,NORB,NORB) ! U TENSOR
      COMPLEX(8)  ,INTENT(IN) :: RHO(NORB,NORB,2,2)     ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT                   ! energy
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
      SUBROUTINE LDAPLUSU_edft(gid,nr,lmnx,lnx,lox,chi,lrx,denmat,etot,ham)
!     **                                                                      **
!     ** DOUBLE COUNTING CORRECTION for the hybrid functional                 **
!     **                                                                      **
      IMPLICIT NONE
      integer(4)  ,intent(in) :: gid
      integer(4)  ,intent(in) :: nr
      integer(4)  ,intent(in) :: lrx
      INTEGER(4)  ,INTENT(IN) :: lmnx       ! #(local orbitals)
      INTEGER(4)  ,INTENT(IN) :: lnx        ! #(radial functions)
      INTEGER(4)  ,INTENT(IN) :: lox(lnx)   ! main angular momentum of local orb.
      real(8)     ,intent(in) :: chi(nr,lnx)
      COMPLEX(8)  ,INTENT(IN) :: denmat(lmnx,lmnx,2,2) ! DENSITY MATRIX
      REAL(8)     ,INTENT(OUT):: ETOT       ! DOUBLE COUNTINNG ENERGY
      COMPLEX(8)  ,INTENT(OUT):: HAM(lmnx,lmnx,2,2)  ! DEtot/D(RHO^*)        
      complex(8)              :: denmat1(lmnx,lmnx,4)
      complex(8)              :: ham1(lmnx,lmnx,4)
      real(8)                 :: r(nr)
      real(8)     ,allocatable:: rho(:,:,:)
      real(8)     ,allocatable:: pot(:,:,:)
      real(8)                 :: edensity(nr)
      real(8)                 :: aux(nr),svar
      integer(4)              :: lmrx,l
      integer(4)              :: idim,lm
      complex(8)  ,parameter  :: ci=(0.d0,1.d0)
      integer(4)  ,parameter  :: ndimd=4
!     **************************************************************************
      lmrx=(lrx+1)**2
      etot=0.d0
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
      allocate(rho(nr,lmrx,4))
      DO IDIM=1,ndimd
        CALL AUGMENTATION_RHO(NR,lnx,LOX,CHI &
     &                       ,LMNX,DENMAT1(:,:,IDIM),LMRX,RHO(:,:,IDIM))
      ENDDO
!
!     ==========================================================================
!     ==  CALCULATE ENERGY AND POTENTIAL                                      ==
!     ==========================================================================
      ALLOCATE(POT(NR,LMRX,NDIMD))
      pot(:,:,:)=0.d0
!     == EXCHANGE ENERGY AND POTENTIAL =========================================
      CALL DFT$SETL4('XCONLY',.TRUE.)
      CALL AUGMENTATION_XC(GID,NR,LMRX,NDIMD,RHO,ETOT,POT)
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
      edensity=edensity*r(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,EDENSITY,SVAR)
print*,'eh ',svar
      ETOT=ETOT+SVAR
!
!     ==========================================================================
!     ==  CALCULATE HAMILTONIAN IN TOTAL/SPIN REPRESENTATION                  ==
!     ==========================================================================
      CALL LDAPLUSU_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX,POT,CHI,HAM1)
      deALLOCATE(POT)
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
      SUBROUTINE ldaplusu_EXPECT(GID,NR,NDIMD,LNX,LOX,LMNX,LMRX &
     &                          ,AEPOT,AEPHI,DATH)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE EXPECTATION VALUE OF                                 **
!     **  THE ONE-CENTER POTENTIALS WITH THE local orbitals                   **
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
      complex(8),INTENT(OUT):: DATH(LMNX,LMNX,NDIMD)
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
                DATH(LMN1,LMN2,ISPIN)=DATH(LMN1,LMN2,ISPIN)+cmplx(SVAR,0.d0)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
