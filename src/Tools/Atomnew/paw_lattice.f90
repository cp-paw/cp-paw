!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE ONEATOM_MODULE
TYPE THISTYPE
  INTEGER(4)  :: GID
  INTEGER(4)  :: NR
  REAL(8)     :: AEZ=-1.D0
  REAL(8)     :: EREF
  INTEGER(4)  :: NB
  INTEGER(4)  :: NC
  INTEGER(4)  :: LOFI(19)
  INTEGER(4)  :: NNOFI(19)
  REAL(8)     :: FOFI(19)
  REAL(8)     :: EOFI(19)
  REAL(8)     :: EKINC
  REAL(8),POINTER :: AEPOT(:)
  REAL(8),POINTER :: RHOCORE(:)
  REAL(8),POINTER :: DREL(:)
  REAL(8),POINTER :: vpauli(:,:)   ! REPULSIVE POTENTIAL
  REAL(8),POINTER :: Udot(:,:)     ! HIGHEST NODLESS PHIDOT FUNCTION
  REAL(8),POINTER :: UN(:,:)       ! HIGHEST NODLESS VALENCE STATE PER L
  REAL(8)         :: Q
  REAL(8)         :: RAD
  REAL(8)         :: DEDRAD
  REAL(8)         :: DEDQ
END TYPE THISTYPE
INTEGER(4)    ,PARAMETER   :: NATX=10
TYPE(THISTYPE),TARGET,SAVE :: THISARR(NATX)
TYPE(THISTYPE),POINTER     :: THIS
LOGICAL(4)    ,PARAMETER   :: TREL=.TRUE.
LOGICAL(4)    ,SAVE        :: TSPECIAL=.FALSE.
END MODULE ONEATOM_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      PROGRAM TEST
!     **************************************************************************
!     ** NEW ATOMIC PROGRAM FOR THE CALCULATIONS OF ATOMIC SETUPS FOR PAW     **
!     **************************************************************************
      USE ONEATOM_MODULE
      USE PERIODICTABLE_MODULE
!      USE LINKEDLIST_MODULE
!      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NAT=2
      REAL(8)              :: AEZARR(NAT)
      REAL(8)              :: DT=10.D0
      REAL(8)              :: MQ=1.D+5
      REAL(8)              :: ANNEQ=1.D-2
      REAL(8)              :: MRAD=5.D+4
      REAL(8)              :: ANNERAD=1.D-2
      INTEGER(4)           :: NITER=1000
INTEGER(4)           :: NP
      REAL(8)              :: ETOT,DETOT,EKINQ,EKINRAD
      REAL(8)              :: DELTA
      INTEGER(4)           :: I,J,IAT,ITER,ir
      REAL(8)              :: RADP(NAT),RAD0(NAT),RADM(NAT)
      REAL(8)              :: QP(NAT),Q0(NAT),QM(NAT)
      REAL(8)              :: DEDQ(NAT)
      REAL(8)              :: DEDRAD(NAT)
      REAL(8)              :: SVAR,SVAR1,SVAR2,SVAR3
      REAL(8)              :: VTOT,VOL0
      REAL(8)              :: VMAD(NAT),EMAD,MADPRESSURE,PRESSURE
      REAL(8)              :: ETOTLAST
      REAL(8)              :: PI
real(8),allocatable :: r(:)
integer(4)          :: nr
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
!      AEZARR(:)=(/8.D0/)   ! O
!     AEZARR(:)=(/14.D0,8.D0,8.D0/)   ! SIO2
      AEZARR(:)=(/100.D0,8.D0/)        ! MGO
      PRINT*,'STARTIN...'
      CALL ATTACHFILES()
      CALL TRACE$SETL4('ON',.FALSE.)
      CALL DFT$SETI4('TYPE',2)
      CALL DFT$SETL4('SPIN',.FALSE.)
      DO IAT=1,NAT
        CALL ONEATOM$NEW(IAT,AEZARR(IAT))
        CALL PERIODICTABLE$GET(NINT(AEZARR(IAT)),'R(COV)',RAD0(IAT))
        Q0(IAT)=0.D0
      ENDDO
Q0(1)=-2.D0
Q0(2)=2.D0
VTOT=(4.212D0/0.529177)**3/4.D0
RAD0=(3.D0/(4.D0*PI)*VTOT/REAL(NAT))**(1.D0/3.D0)
      QM(:)=Q0(:)
      RADM(:)=RAD0(:)
!
!     ==========================================================================
!     == PRINT                                                                ==
!     ==========================================================================
      CALL REPORT$TITLE(6,'INPUT DATA')
      CALL REPORT$R8VAL(6,'MASS FOR CHARGE DYNAMICS',MQ,'A.U.')
      CALL REPORT$R8VAL(6,'FRICTION FOR CHARGE DYNAMICS',ANNEQ,'A.U.')
      CALL REPORT$R8VAL(6,'MASS FOR RADIUS DYNAMICS',MRAD,'A.U.')
      CALL REPORT$R8VAL(6,'FRICTION FOR RADIUS DYNAMICS',ANNERAD,'A.U.')
      CALL REPORT$I4VAL(6,'NUMBER OF ATOMS',NAT,'')
      DO IAT=1,NAT
        CALL ONEATOM$REPORTATOM(6,IAT)
      ENDDO
!
!     ==========================================================================
!     == this is to construct vpauli and the nodeless scattering function udot==
!     ==========================================================================
      CALL ONEATOM(1,0.d0,3.d0,DETOT,DEDQ(1),DEDRAD(1))
      nr=thisarr(1)%nr
      allocate(r(nr))
      call radial$r(Thisarr%gid,nr,r)
      OPEN(1264,FILE='DAT')
      DO IR=1,NR
        WRITE(1264,fmt='(30e15.5)')R(IR),THISARR(1)%UN(IR,1:2),THISARR(1)%UDOT(IR,1:2)
      ENDDO
      Close(1264)
      OPEN(1265,FILE='vpauli')
      DO IR=1,NR
        WRITE(1265,fmt='(30e15.5)')R(IR),THISARR(1)%vpauli(IR,:)
      ENDDO
      Close(1265)
      deallocate(r)
stop 'ok'
!
!     ==========================================================================
!     == PRINT                                                                ==
!     ==========================================================================
!!$TSPECIAL=.TRUE.
      VTOT=4.D0*PI/3.D0*SUM(RAD0(:)**3)
      NP=10
      OPEN(UNIT=8,FILE='XX_0.DAT')
      OPEN(UNIT=9,FILE='XX_M2.DAT')
      OPEN(UNIT=14,FILE='XX_P2.DAT')
      OPEN(UNIT=11,FILE='XX_P1.DAT')
      OPEN(UNIT=12,FILE='XX_M1.DAT')
      DO I=1,NP
        Q0(1)=-0.D0
        Q0(2)=-Q0(1)
        RAD0(1)=3.0D0+4.0D0*(-0.5D0+REAL(I-1,KIND=8)/REAL(NP-1,KIND=8))
        if(vtot-4.d0*pi/3*rad0(1)**3.le.0.d0) cycle
        RAD0(2)=(3.D0/(4.D0*PI)*VTOT-RAD0(1)**3)**(1.D0/3.D0)
        ETOT=0.D0
        CALL ONEATOM(1,Q0(1),RAD0(1),DETOT,DEDQ(1),DEDRAD(1))
        ETOT=ETOT+DETOT
        CALL ONEATOM(2,Q0(2),RAD0(2),DETOT,DEDQ(2),DEDRAD(2))
        ETOT=ETOT+DETOT
        CALL ESTRUCTURE(NAT,Q0,VTOT,EMAD,VMAD,MADPRESSURE)
        ETOT=ETOT+EMAD
        DEDRAD=DEDRAD/(4.D0*PI*RAD0**2)
        WRITE(8,FMT='(10F20.15)')RAD0(1)**3/SUM(RAD0**3),ETOT,ETOT+EMAD
!=====
        Q0(1)=-2.D0
        Q0(2)=-Q0(1)
        RAD0(1)=3.0D0+4.0D0*(-0.5D0+REAL(I-1,KIND=8)/REAL(NP-1,KIND=8))
        RAD0(2)=(3.D0/(4.D0*PI)*VTOT-RAD0(1)**3)**(1.D0/3.D0)
RAD0=RAD0*1.3D0
        ETOT=0.D0
        CALL ONEATOM(1,Q0(1),RAD0(1),DETOT,DEDQ(1),DEDRAD(1))
        ETOT=ETOT+DETOT
        CALL ONEATOM(2,Q0(2),RAD0(2),DETOT,DEDQ(2),DEDRAD(2))
        ETOT=ETOT+DETOT
        CALL ESTRUCTURE(NAT,Q0,VTOT,EMAD,VMAD,MADPRESSURE)
        DEDRAD=DEDRAD/(4.D0*PI*RAD0**2)
        WRITE(9,FMT='(10F20.15)')RAD0(1)**3/SUM(RAD0**3),ETOT,ETOT+EMAD
!=====
        Q0(1)=2.D0
        Q0(2)=-Q0(1)
        RAD0(1)=3.0D0+4.0D0*(-0.5D0+REAL(I-1,KIND=8)/REAL(NP-1,KIND=8))
        RAD0(2)=(3.D0/(4.D0*PI)*VTOT-RAD0(1)**3)**(1.D0/3.D0)
RAD0=RAD0*1.3D0
        ETOT=0.D0
        CALL ONEATOM(1,Q0(1),RAD0(1),DETOT,DEDQ(1),DEDRAD(1))
        ETOT=ETOT+DETOT
        CALL ONEATOM(2,Q0(2),RAD0(2),DETOT,DEDQ(2),DEDRAD(2))
        ETOT=ETOT+DETOT
        CALL ESTRUCTURE(NAT,Q0,VTOT,EMAD,VMAD,MADPRESSURE)
        DEDRAD=DEDRAD/(4.D0*PI*RAD0**2)
        WRITE(14,FMT='(10F20.15)')RAD0(1)**3/SUM(RAD0**3),ETOT,ETOT+EMAD
!=====
        Q0(1)=1.D0
        Q0(2)=-Q0(1)
        RAD0(1)=3.0D0+4.0D0*(-0.5D0+REAL(I-1,KIND=8)/REAL(NP-1,KIND=8))
        RAD0(2)=(3.D0/(4.D0*PI)*VTOT-RAD0(1)**3)**(1.D0/3.D0)
RAD0=RAD0*1.3D0
        ETOT=0.D0
        CALL ONEATOM(1,Q0(1),RAD0(1),DETOT,DEDQ(1),DEDRAD(1))
        ETOT=ETOT+DETOT
        CALL ONEATOM(2,Q0(2),RAD0(2),DETOT,DEDQ(2),DEDRAD(2))
        ETOT=ETOT+DETOT
        CALL ESTRUCTURE(NAT,Q0,VTOT,EMAD,VMAD,MADPRESSURE)
        DEDRAD=DEDRAD/(4.D0*PI*RAD0**2)
        WRITE(11,FMT='(10F20.15)')RAD0(1)**3/SUM(RAD0**3),ETOT,ETOT+EMAD
!=====
        Q0(1)=-1.D0
        Q0(2)=-Q0(1)
        RAD0(1)=3.0D0+4.0D0*(-0.5D0+REAL(I-1,KIND=8)/REAL(NP-1,KIND=8))
        RAD0(2)=(3.D0/(4.D0*PI)*VTOT-RAD0(1)**3)**(1.D0/3.D0)
RAD0=RAD0*1.5D0
        ETOT=0.D0
        CALL ONEATOM(1,Q0(1),RAD0(1),DETOT,DEDQ(1),DEDRAD(1))
        ETOT=ETOT+DETOT
        CALL ONEATOM(2,Q0(2),RAD0(2),DETOT,DEDQ(2),DEDRAD(2))
        ETOT=ETOT+DETOT
        CALL ESTRUCTURE(NAT,Q0,VTOT,EMAD,VMAD,MADPRESSURE)
        DEDRAD=DEDRAD/(4.D0*PI*RAD0**2)
        WRITE(12,FMT='(10F20.15)')RAD0(1)**3/SUM(RAD0**3),ETOT,ETOT+EMAD
      ENDDO
      CLOSE(8)
      CLOSE(9)
      CLOSE(14)
      CLOSE(11)
      CLOSE(12)
!==================
      STOP
!
      VTOT=4.D0*PI/3.D0*SUM(RAD0(:)**3)
      ETOTLAST=1.D+5
      DO ITER=1,NITER
!
!       ========================================================================
!       == PROPAGATE
!       ========================================================================
        ETOT=0.D0
        DO IAT=1,NAT
          CALL ONEATOM(IAT,Q0(IAT),RAD0(IAT),DETOT,DEDQ(IAT),DEDRAD(IAT))
          ETOT=ETOT+DETOT
        ENDDO
        PRESSURE=-4.D0*PI*SUM(DEDRAD(:)*RAD0(:)**2)
!       == MADELUNG CONTRIBUTION ===============================================
        VOL0=4.D0*PI/3.D0*SUM(RAD0(:)**3)
        CALL ESTRUCTURE(NAT,Q0,VOL0,EMAD,VMAD,MADPRESSURE)
        ETOT=ETOT+EMAD
        DEDQ(:)=DEDQ(:)+VMAD(:)
        PRESSURE=PRESSURE+MADPRESSURE
!
!       ========================================================================
!       == PROPAGATE                                                          ==
!       ========================================================================
        SVAR1=2.D0/(1.D0+ANNEQ)
        SVAR2=1.D0-SVAR1
        SVAR3=-DT**2/MQ/(1.D0+ANNEQ)
        QP=SVAR1*Q0+SVAR2*QM+SVAR3*DEDQ
        SVAR=SUM(QP(:))/REAL(NAT,KIND=8)  !CHARGE CONSERVATION
        QP=QP-SVAR
        SVAR1=2.D0/(1.D0+ANNERAD)
        SVAR2=1.D0-SVAR1
        SVAR3=-DT**2/MRAD/(1.D0+ANNERAD)
        RADP=SVAR1*RAD0+SVAR2*RADM+SVAR3*DEDRAD
        DO I=1,300
          SVAR1=VTOT-4.D0*PI/3.D0*SUM(RADP(:)**3)
          IF(ABS(SVAR1/VTOT).LT.1.D-10) EXIT
          SVAR1=SVAR1/((4.D0*PI)**2*SUM(RADP(:)**2*RAD0(:)**2))
          RADP(:)=RADP(:)+4.D0*PI*RAD0(:)**2*SVAR1
          IF(I.EQ.300) PRINT*,'CONVERGENCE PROBLEM VOLUME CONSTRAINT'
        ENDDO
!WRITE(*,FMT='("DVOL/VOL ",10E20.5)')SVAR1/VTOT,DEDQ/(4.D0*PI*RAD0**2)
!
!       ========================================================================
!       == REPORT ENERGIES                                                    ==
!       ========================================================================
        EKINQ=0.5D0*MQ*SUM((QP-QM)**2)/(2.D0*DT)**2          
        EKINRAD=0.5D0*MRAD*SUM((RADP-RADM)**2)/(2.D0*DT)**2          
        WRITE(*,FMT='("!>",I5,15F10.3)')ITER,EKINQ+EKINRAD,ETOT,EKINQ+EKINRAD+ETOT,Q0,RAD0,PRESSURE
        IF(MOD(ITER,100).EQ.0) THEN
          DO IAT=1,NAT
            CALL ONEATOM$REPORTATOM(6,IAT)
          ENDDO
        END IF
!
!       ========================================================================
!       == SWITCH                                                             ==
!       ========================================================================
        QM=Q0
        Q0=QP
        RADM=RAD0
        RAD0=RADP
      ENDDO
      DO IAT=1,NAT
        CALL ONEATOM$REPORTATOM(6,IAT)
      ENDDO
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONEATOM$REPORTATOM(NFIL,IAT)
!     **************************************************************************
      USE ONEATOM_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      INTEGER(4),INTENT(IN) :: IAT
      INTEGER(4)            :: IB
      REAL(8)               :: EV
      REAL(8)               :: ANGSTROM
      REAL(8)               :: BAR
      REAL(8)               :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      CALL CONSTANTS('EV',EV)
      CALL CONSTANTS('ANGSTROM',ANGSTROM)
      CALL CONSTANTS('BAR',BAR)
      THIS=>THISARR(IAT)
      WRITE(NFIL,FMT='("ATOMIC NUMBER ",F5.1)')THIS%AEZ
      CALL REPORT$R8VAL(NFIL,'Q',THIS%Q,'-E')
      CALL REPORT$R8VAL(NFIL,'RAD',THIS%RAD/ANGSTROM,'ANGSTROM')
      CALL REPORT$R8VAL(NFIL,'DEDQ ',THIS%DEDQ/EV,'EV')
      CALL REPORT$R8VAL(NFIL,'DEDRAD ',THIS%DEDRAD/EV,'EV/A0')
      CALL REPORT$R8VAL(NFIL,'PRESSURE ',-THIS%DEDRAD/(4.D0*PI*THIS%RAD**2)/BAR,'BAR')
      DO IB=THIS%NC+1,THIS%NB
        WRITE(NFIL,FMT='("L=",I1," E[EV]=",F10.3," F=",F10.5)') &
     &             THIS%LOFI(IB),THIS%EOFI(IB)/EV,THIS%FOFI(IB)
      ENDDO       
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ATTACHFILES()
!     **************************************************************************
!     **                                                                      **
!     **  READS THE NAME OF THE CONTROL FILE FROM THE ARGUMENT LIST AND       **
!     **  DEFINES THE FILES IN THE FILEHANDLER OBJECT                         **
!     **                                                                      **
!     **  THE FILES DEFINED ARE                                               **
!     **    PROT     PROTOCOLL FILE                                           **
!     **    INPUT    CONTROL FILE X.ACNTL                                     **
!     **    STP      SETU OUTPUT FILE                                         **
!     **                                                                      **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      CHARACTER(128)         :: CNTLNAME
      CHARACTER(128)         :: ROOTNAME
      INTEGER(4)             :: I
!     **************************************************************************

!     ==========================================================================
!     ==  COLLECT INPUT ARGUMENT AND DEFINES ROOTNAME                         ==
!     ==========================================================================
      CALL GETARG(1,CNTLNAME)
      I=INDEX(CNTLNAME,'.',BACK=.TRUE.)
      ROOTNAME=CNTLNAME(1:I-1)
      CALL FILEHANDLER$SETROOT(ROOTNAME)
!
!     ==========================================================================
!     ==  ATTACH PROTOCOLL FILE  X.PROT                                       ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('PROT',.FALSE.,TRIM(ROOTNAME)//-'.PROT')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('PROT','FORM','FORMATTED')
!
!     ==========================================================================
!     == DEFINE AND READ INPUT FILE X.ACNTL                                   ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('INPUT',.FALSE.,TRIM(ROOTNAME)//-'.ACNTL')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('INPUT','FORM','FORMATTED')
!
!     ==========================================================================
!     == DEFINE SETUP FILE                                                    ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('SETUP',.FALSE.,TRIM(ROOTNAME)//-'.STP')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','STATUS','UNKNOWN')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','ACTION','WRITE')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP','FORM','FORMATTED')
!
!     ==========================================================================
!     == ATTACH SETUP FILE AS INPUT FILE FOR TESTING                          ==
!     ==========================================================================
      CALL FILEHANDLER$SETFILE('SETUP_IN',.FALSE.,TRIM(ROOTNAME)//-'.STP')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP_IN','STATUS','OLD')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP_IN','POSITION','REWIND')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP_IN','ACTION','READ')
      CALL FILEHANDLER$SETSPECIFICATION('SETUP_IN','FORM','FORMATTED')
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE WRITEPHI(FILE,GID,NR,NPHI,PHI)
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
      SUBROUTINE WRITECOMPARE(FILE,GID,NR,NPHI,PHI1,PHI2)
!     **                                                                      **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE
      INTEGER(4)  ,INTENT(IN) :: GID      ! GRID ID
      INTEGER(4)  ,INTENT(IN) :: NR       ! #(RDIAL GRID POINTS)
      INTEGER(4)  ,INTENT(IN) :: NPHI        
      REAL(8)     ,INTENT(IN) :: PHI1(NR,NPHI)
      REAL(8)     ,INTENT(IN) :: PHI2(NR,NPHI)
      INTEGER(4)              :: IR
      REAL(8)                 :: R(NR)
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R)
      OPEN(100,FILE=FILE)
      DO IR=1,NR
        WRITE(100,FMT='(F15.10,2X,8(F25.15,2X))')R(IR),PHI1(IR,:),PHI2(IR,:)
      ENDDO
      CLOSE(100)
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE AESCF(GID,NR,NB,TREL,LOFI,SO,F,NN,Z,RHOADD,RBOX,DREL,POT,EOFI)
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID
      INTEGER(4) ,INTENT(IN)     :: NR
      INTEGER(4) ,INTENT(IN)     :: NB        ! #(STATES)
      LOGICAL(4) ,INTENT(IN)     :: TREL
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)  !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO(NB)    !SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: F(NB)     !OCCUPATION
      INTEGER(4) ,INTENT(IN)     :: NN(NB)    !#(NODES)
      REAL(8)    ,INTENT(IN)     :: Z         !ATOMIC NUMBER
      REAL(8)    ,INTENT(IN)     :: RHOADD(NR)! FIXED CONTRIBUTION OF DENSITY
      REAL(8)    ,INTENT(IN)     :: RBOX      !ATOM ENCLOSED IN A BOX 
      REAL(8)    ,INTENT(OUT)    :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(INOUT)  :: POT(NR)   !POTENTIAL
      REAL(8)    ,INTENT(OUT)    :: EOFI(NB)  !ONE-PARTICLE ENERGIES
      REAL(8)                    :: R(NR)
      REAL(8)                    :: RHO(NR)
      REAL(8)                    :: EREF
      INTEGER(4)                 :: ITER
      INTEGER(4)                 :: NITER=50
      REAL(8)                    :: XAV,XMAX
      LOGICAL(4)                 :: CONVG
      REAL(8)   ,PARAMETER       :: TOL=1.D-5
      INTEGER(4)                  :: NFILO
      REAL(8)                    :: EH,EXC
      LOGICAL(4),PARAMETER       :: TBROYDEN=.true.
      REAL(8)                    :: POTIN(NR)
      INTEGER(4)                 :: IRBOX
      INTEGER(4)                 :: IR
      CHARACTER(32)              :: ID
real(8)                    :: pots(nr,100)
!     ***********************************************************************
      CALL RADIAL$R(GID,NR,R)
      IRBOX=0
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GT.RBOX) EXIT
      ENDDO
      EOFI=0.D0
      EREF=0.D0
      XAV=0.D0
      XMAX=0.D0
      CONVG=.FALSE.
      IF(TBROYDEN) THEN
        CALL BROYDEN$NEW(NR,4,1.D-1)
      ELSE
        CALL MYMIXPOT('ON',GID,NR,POT,XMAX,XAV)
      END IF
      DO ITER=1,NITER
        IF(TREL) THEN
          CALL SCHROEDINGER$DREL(GID,NR,POT,EREF,DREL)
          ID='FULL'
        ELSE
          DREL=0.D0
          ID='NONREL'
        END IF
pots(:,iter)=pot
        CALL AERHO(ID,GID,NR,NB,LOFI,SO,F,NN,RBOX,DREL,POT,RHO,EOFI)
        RHO(:)=RHO(:)+RHOADD(:)
        EREF=EOFI(NB)
!
!       ================================================================
!       ==  EXIT IF CONVERGED                                         ==
!       ================================================================
        IF(CONVG) THEN
          CALL FILEHANDLER$UNIT('PROT',NFILO)
          CALL REPORT$I4VAL(NFILO,"SCF LOOP CONVERGED AFTER ",ITER," ITERATIONS")
          CALL REPORT$R8VAL(NFILO,'AV. DIFF BETWEEN IN- AND OUTPUT POTENTIAL',XAV,'H')
          CALL REPORT$R8VAL(NFILO,'MAXIMUM  OF (VIN-VOUT)*R^2',XMAX,'H')
          IF(TBROYDEN) THEN
            CALL BROYDEN$CLEAR
          ELSE
            CALL MYMIXPOT('OFF',GID,NR,POT,XMAX,XAV)
          END IF
          RETURN
        END IF

!       ====================================================================
!       == CALCULATE OUTPUT POTENTIAL                                     ==
!       ====================================================================
        IF(TBROYDEN) POTIN=POT
        CALL MYVOFRHO(GID,NR,Z,RHO,POT,EH,EXC)
!
!       ================================================================
!       ==  GENERATE NEXT ITERATION USING D. G. ANDERSON'S METHOD     ==
!       ================================================================
        IF(TBROYDEN) THEN
          XAV=SQRT(DOT_PRODUCT(POT-POTIN,POT-POTIN)/REAL(NR,KIND=8))
          XMAX=MAXVAL(ABS(POT-POTIN)) 
          CALL BROYDEN$STEP(NR,POTIN,(POT-POTIN))
          POT=POTIN
       ELSE
          CALL MYMIXPOT('GO',GID,NR,POT,XMAX,XAV)
        END IF
PRINT*,'XMAX ',XMAX,XAV
        CONVG=(XMAX.LT.TOL)
      ENDDO
      CALL WRITEPHI('POTS.DAT',GID,NR,10,POTS(:,1:10))
      CALL ERROR$MSG('SELFCONSISTENCY LOOP NOT CONVERGED')
      CALL ERROR$STOP('AESCF')
      RETURN
      END
!
!     .....................................................VOUT.........
      SUBROUTINE MYVOFRHO(GID,NR,AEZ,RHO,POT,EH,EXC)
!     ******************************************************************
!     **                                                              **
!     **  ELECTROSTATIC AND EXCHANGE-CORRELATION POTENTIAL            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: AEZ
      REAL(8),   INTENT(IN) :: RHO(NR)
      REAL(8),   INTENT(OUT):: POT(NR)
      REAL(8),   INTENT(OUT):: EH
      REAL(8),   INTENT(OUT):: EXC
      REAL(8)               :: PI
      REAL(8)               :: FOURPI
      REAL(8)               :: Y0
      REAL(8)               :: AUX(NR)
      REAL(8)               :: EDEN(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      INTEGER(4)            :: IR
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CALL RADIAL$R(GID,NR,R)
!
!     ==================================================================
!     ==  TOTAL POTENTIAL                                             ==
!     ==================================================================
      EDEN=0.D0
      CALL RADIAL$POISSON(GID,NR,0,RHO,POT)
      EDEN(:)=EDEN(:)+0.5D0*RHO(:)*POT(:)
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AUX)
      POT(:)=POT(:)+AUX(:)
      EDEN(:)=EDEN(:)+RHO(:)*AUX(:)
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,EDEN,EH)
!
!     ==================================================================
!     ==  EXCHANGE CORRELATION
!     ==================================================================
      CALL RADIAL$DERIVE(GID,NR,RHO(:),GRHO)
      DO IR=1,NR
        RH=RHO(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC1,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        EDEN(IR)=4.D0*PI*EXC1   ! ANGULAR INTEGRATION ALREADY INCLUDED
        POT(IR)=POT(IR)+VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      CALL RADIAL$DERIVE(GID,NR,GRHO(:),AUX)
      POT(:)=POT(:)-AUX
      IF(R(1).GT.1.D+10) THEN
         POT(:)=POT(:)-2.D0/R(:)*GRHO(:)
      ELSE
        POT(2:)=POT(2:)-2.D0/R(2:)*GRHO(2:)
        POT(1)=POT(1)-2.D0/R(2)*GRHO(2)
        POT(1)=POT(2)
      END IF
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,EDEN,EXC)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE MYMIXPOT(SWITCH,GID,NR,POT,XMAX,XAV)
!     ******************************************************************
!     **                                                              **
!     **  MIX POTENTIAL USING D.G.ANDERSEN'S METHOD                   **
!     **                                                              **
!     **  1) INITIALIZE WITH SWITCH='ON'                              **
!     **     STORES THE FIRST INPUT POTENTIAL <-POT                   **
!     **     ALLOCATES  INTERNAL ARRAYS                               **
!     **                                                              **
!     **  2) ITERATE WITH WITH SWITCH='GO'                            **
!     **     RECEIVES THE OUTPUT POTENTIAL    <-POT           E       **
!     **     CALCULATES NEW INPUT POTENTIAL   ->POT                   **
!     **     STORES THE NEW INPUT POTENTIAL                           **
!     **                                                              **
!     **  3) CLEAR MEMORY WITH SWITCH='OFF'                           **
!     **                                                              **
!     **  WARNING! DO NOT USE SIMULTANEOUSLY FOR TOW DIFFERENT SCHEMES**
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)   :: SWITCH !CAN BE 'ON', 'GO' OR 'OFF'
      INTEGER(4)  ,INTENT(IN)   :: GID
      INTEGER(4)  ,INTENT(IN)   :: NR
      REAL(8)     ,INTENT(INOUT):: POT(NR) ! 
      REAL(8)     ,INTENT(OUT)  :: XMAX    ! MAX.|R**2*(VOUT-VIN)|
      REAL(8)     ,INTENT(OUT)  :: XAV     ! <R**2(VOUT-VIN)**2>
      REAL(8)     ,PARAMETER    :: ALPHA=1.D-1 !MIXING PARAMETER
      CHARACTER(8)       ,SAVE :: STATUS='OFF'
      LOGICAL(4)         ,SAVE :: TSTART
      INTEGER(4)         ,SAVE :: NRSAVE
      REAL(8),ALLOCATABLE,SAVE :: OLDPOTIN(:)
      REAL(8),ALLOCATABLE,SAVE :: OLDPOTOUT(:)
      REAL(8),ALLOCATABLE,SAVE :: NEWPOTOUT(:)
      REAL(8),ALLOCATABLE,SAVE :: NEWPOTIN(:)
      REAL(8)                  :: BETA
      REAL(8)                  :: SVAR1,SVAR2
      REAL(8)                  :: R(NR)
      REAL(8)                  :: AUX(NR),AUX1(NR),AUX2(NR)
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
      IF(SWITCH.EQ.'GO') THEN
        IF(NR.NE.NRSAVE) THEN
          CALL ERROR$MSG('INCONSISTENT NUMBER OF GRIDPOINTS')
          CALL ERROR$STOP('MIXPOT')
        END IF
!       ================================================================
!       == COPY POT INTO NEWPOTIN                                     ==
!       ================================================================
        NEWPOTOUT(:)=POT(:)
!       ================================================================
!       == CALCULATE MAX AND VARIANCE OF POTOUT-POTIN                 ==
!       ================================================================
        AUX(:)=NEWPOTOUT(:)-NEWPOTIN(:)
        XMAX=MAXVAL(ABS(AUX(:)))
        XAV=SUM(AUX(:)**2)
        SVAR1=REAL(NR,KIND=8)
        XAV=SQRT(XAV/SVAR1)
!       ================================================================
!       ==  CALCULATE MIXING FACTOR BETA                              ==
!       ================================================================
        IF(TSTART) THEN
          TSTART=.FALSE.
          BETA=0.D0
        ELSE
          AUX1(:)=NEWPOTOUT(:)-NEWPOTIN(:)
          AUX2(:)=OLDPOTOUT(:)-OLDPOTIN(:)
          AUX(:)=AUX1(:)-AUX2(:)
          SVAR1=SUM(AUX1(:)*AUX(:))
          SVAR2=SUM(AUX(:)**2)
          BETA=SVAR1/SVAR2
        END IF
!       ================================================================
!       == MIX POTENTIALS                                             ==
!       ================================================================
        AUX1(:)=(1.D0-BETA)*NEWPOTIN(:) + BETA*OLDPOTIN(:)
        AUX2(:)=(1.D0-BETA)*NEWPOTOUT(:)+ BETA*OLDPOTOUT(:)
        POT(:)=AUX1(:) + ALPHA*(AUX2(:)-AUX1(:))
        OLDPOTIN(:) =NEWPOTIN(:)
        OLDPOTOUT(:)=NEWPOTOUT(:)
        NEWPOTIN(:) =POT(:)
!
!     ==================================================================
!     == INITIALIZE MIXING                                            ==
!     ==================================================================
      ELSE IF(SWITCH.EQ.'ON') THEN
        IF(TRIM(STATUS).NE.'OFF') THEN
          DEALLOCATE(OLDPOTIN)
          DEALLOCATE(OLDPOTOUT)
          DEALLOCATE(NEWPOTIN)
          DEALLOCATE(NEWPOTOUT)
        END IF
        NRSAVE=NR
        ALLOCATE(OLDPOTIN(NRSAVE))
        ALLOCATE(OLDPOTOUT(NRSAVE))
        ALLOCATE(NEWPOTIN(NRSAVE))
        ALLOCATE(NEWPOTOUT(NRSAVE))
        TSTART=.TRUE.
        STATUS='ON'
        OLDPOTOUT(:)=0.D0
        OLDPOTIN(:) =0.D0
        NEWPOTOUT(:)=0.D0
        NEWPOTIN(:) =POT(:)
!
!     ==================================================================
!     == CLEAR ARRAYS                                                 ==
!     ==================================================================
      ELSE IF(SWITCH.EQ.'OFF') THEN
        DEALLOCATE(OLDPOTIN)
        DEALLOCATE(OLDPOTOUT)
        DEALLOCATE(NEWPOTIN)
        DEALLOCATE(NEWPOTOUT)
        TSTART=.FALSE.
        STATUS='OFF'
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE AERHO(ID,GID,NR,NB,LOFI,SO,F,NN,RBOX,DREL,POT,RHO,E)
!     **                                                                      **
!     **  DETERMINES THE DENSITY FOR A GIVEN POTENTIAL AND A SPECIFIED SET    **
!     **  OF WAVE FUNCTIONS.                                                  **
!     **  UPDATES THE ONE-PARTICLE ENERGIES AND RELATIVISTIC CORRECTION "DREL"**
!     **                                                                      **
!     **  THE ORBITALS ARE SPECIFIED BY ANGULAR MOMENTUM "LOFI",              **
!     **                             THE NUMBER OF NODES "NN",                **
!     **                             AND THE SPIN ORBIT PARAMETER "SO"        **
!     **                                                                      **
!     **  VALUES OF ID:                                                       **
!     **    NONREL: NONRELATIVISTIC CALCULATION, DREL=0 IS RETURNED           **
!     **    EFFZORA: THE RELATIVISTIC CORRECTIONS ARE TAKEN FROM INPUT        **
!     **            AND REMAIN UNCHANED                                       **
!     **    FULL: PERFORMS RELATIVISTIC CALCULATION USING DREL FROM INPUT     **
!     **          AND RECALCULATES DREL FROM THE MEAN KINETIC ENERGY DENSITY  **
!     **          DIVIDED BY THE DENSITY.                                     **
!     **                                                                      **
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID
      INTEGER(4) ,INTENT(IN)     :: GID       ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR        ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: NB        ! #(STATES)
      INTEGER(4) ,INTENT(IN)     :: LOFI(NB)  !ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO(NB)    !SWITCH FOR SPIN-ORBIT COUP.
      INTEGER(4) ,INTENT(IN)     :: NN(NB)    !#(NODES)
      REAL(8)    ,INTENT(IN)     :: F(NB)     !OCCUPATION
      REAL(8)    ,INTENT(IN)     :: RBOX      !ATOM ENCLOSED IN A BOX
      REAL(8)    ,INTENT(INOUT)  :: DREL(NR)  !RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: POT(NR)   !POTENTIAL
      REAL(8)    ,INTENT(OUT)    :: RHO(NR)   !DENSITY
      REAL(8)    ,INTENT(INOUT)  :: E(NB)     !ONE-PARTICLE ENERGIES
      REAL(8)                    :: R(NR)     !RADIAL GRID POINTS
      REAL(8)                    :: G(NR)     !INHOMOGENEITY (NOT USED)
      REAL(8)                    :: AUX(NR),aux1(nr)
      REAL(8)                    :: PHI(NR)
      REAL(8)                    :: EKIN(NR)
      REAL(8)                    :: EKIN2(NR)
      REAL(8)                    :: ARRAY(NR,5)
      INTEGER(4)                 :: IB,IR
      REAL(8)                    :: SVAR
      REAL(8)                    :: C0LL,PI,Y0
!     **************************************************************************
      IF(ID.NE.'FULL'.AND.ID.NE.'EFFZORA'.AND.ID.NE.'NONREL') THEN
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$MSG('ALLOWED VALUES ARE "FULL",EFFZORA","NONREL"')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('AESCF')
      ENDIF
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!
!     ==========================================================================
!     == DETERMINE BOUND STATES FOR A GIVEN POTENTIAL AND ADD TO DENSITY      ==
!     ==========================================================================
      RHO(:)=0.D0
      EKIN(:)=0.D0
      EKIN2(:)=0.D0
      DO IB=1,NB
         G(:)=0.D0
         IF(ID.EQ.'FULL') THEN
           CALL SCHROEDINGER$DREL(GID,NR,POT,E(IB),DREL)
         ELSE IF(ID.EQ.'NONREL') THEN
           DREL(:)=0.D0 
         END IF
         CALL BOUNDSTATE(GID,NR,LOFI(IB),SO(IB),RBOX,DREL,G,NN(IB),POT,E(IB),PHI)
         AUX(:)=(R(:)*PHI(:))**2
         CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
         CALL RADIAL$VALUE(GID,NR,AUX1,RBOX,SVAR)
         PHI(:)=PHI(:)/SQRT(SVAR)
         RHO(:)  =RHO(:)  +F(IB)*C0LL         *PHI(:)**2
         EKIN(:) =EKIN(:) +F(IB)*C0LL*E(IB)   *PHI(:)**2
         EKIN2(:)=EKIN2(:)+F(IB)*C0LL*E(IB)**2*PHI(:)**2
      ENDDO
      do ir=1,nr
        if(r(ir).gt.rbox) rho(ir)=0.d0
      enddo
!
!     ==========================================================================
!     == DETERMINE RELATIVISTIC CORRECTION                                    ==
!     ==========================================================================
      IF(ID.EQ.'FULL') THEN
        AUX(:)=MAX(RHO(:),1.D-6)
        AUX=1.D0/AUX
        EKIN(:)=EKIN(:)*AUX(:)
        EKIN2(:)=EKIN2(:)*AUX(:)
        EKIN2(:)=MAX(EKIN2(:)-EKIN(:)**2,0.D0)
        EKIN2(:)=SQRT(EKIN2(:))
        DO IR=1,NR
          IF(RHO(IR).LT.1.D-5) THEN
            EKIN(IR)=0.D0
            EKIN2(IR)=0.D0
          END IF
        ENDDO
        EKIN(:)=EKIN(:)/Y0-POT(:)
        CALL SCHROEDINGER$DREL(GID,NR,-EKIN,0.D0,DREL)
      END IF
      IF(ID.EQ.'FULL') THEN
        ARRAY(:,1)=POT(:)*Y0
        ARRAY(:,2)=(EKIN+POT(:))*Y0
        ARRAY(:,3)=(EKIN+POT(:))*Y0+EKIN2(:)
        ARRAY(:,4)=(EKIN+POT(:))*Y0-EKIN2(:)
        CALL WRITEPHI('EKIN.DAT',GID,NR,4,ARRAY)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE BOUNDSTATE(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,PHI)
!     **************************************************************************
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND         **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                     **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G                **
!     **                                                                      **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE                 **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.         **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.         **
!     **                                                                      **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND           **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT           **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                         **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: RBOX    ! BOX RADIUS
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-12
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: DER,DERO
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIINHOM(NR),GHOM(NR)
      REAL(8)                    :: PHI1(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20 ! MAXIMUM FACTOR IN THE WAVE FUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      REAL(8)                    :: ROUT
      REAL(8)                    :: VAL,VAL1,VAL2,r1,r2
      INTEGER(4)                 :: IROUT,IRCL,IRBOX,irend
!     *********************************************************************
                                 CALL TRACE$PUSH('BOUNDSTATE')
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
!     ==  r(irbox) is the first gridpoint just ioutside the box
      IRBOX=1
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GT.RBOX) EXIT
      ENDDO
!          
      ISTART=1
      X0=E
      DX=1.D-2
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO I=1,NITER
        E=X0
!        E=MIN(X0,EMAX)
!
!       =======================================================================
!       ==  CUT OFF THE POTENTIAL                                            ==
!       =======================================================================
        CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)
!       == boundary condition phi(rout)=0 =====================================
        if(r(irout).lt.rbox) then
          rout=r(irout)-1.d-5    !ensure that rout<r(irout)
        else
          rout=rbox
          irout=irbox
        end if
!       ==  set kinetic energy to zero beyond rout
        POT1(:)=POT(:)
        POT1(IROUT:)=POT(IROUT)
!
!       =======================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       =======================================================================
        IDIR=1
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(IROUT).GT.0.OR.PHI(IROUT).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NN',NN)
          CALL ERROR$STOP('BOUNDSTATE')
        END IF
!
!       =======================================================================
!       == ESTIMATE PHASE SHIFT                                              ==
!       =======================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,ROUT,Z0)
        Z0=Z0-REAL(NN+1)
        IF(ABS(2.D0*DX).LE.TOL) EXIT
!       =====================================================================
!       ==  BISECTION                                                      ==
!       =====================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('BOUNDSTATE')
      END IF
!
!     =======================================================================
!     ==  DETERMINE MATCHING POINT                                         ==
!     =======================================================================
      IRMATCH=IRCL
      IF(R(IRCL).GT.5.D0) THEN
        CALL RADIAL$XOFR(GID,5.D0,SVAR)
        IRMATCH=INT(SVAR)
      END IF
!
!     =======================================================================
!     ==  INTEGRATE INWARD                                                 ==
!     =======================================================================
      IF(IRMATCH.LT.IROUT) THEN
        THOM=MAXVAL(ABS(G(:))).EQ.0.D0
        IDIR=-1
!       ==  homogeneous solution that fulfills the outer boundary condition
        irend=min(irout,nr-3)
        GHOM(:)=0.D0
        GHOM(IRend+1)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
        phihom(:)=phihom(:)/phihom(irmatch)
        GHOM(:)=0.D0
        GHOM(Irend+2)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHI1)
        PHI1(:)=PHI1(:)/PHI1(IRMATCH)
!       == extrapolate =======================
        if(irend.lt.irout) then
          R1=R(IREND)
          R2=R(IREND+1)
          VAL1=PHIHOM(IREND)
          VAL2=PHIHOM(IREND+1)
          PHIHOM(IREND:)=VAL1+(VAL2-VAL1)/(R2-R1)*(R(IREND:)-R1)
          VAL1=PHI1(IREND)
          VAL2=PHI1(IREND+1)
          PHI1(IREND:)=VAL1+(VAL2-VAL1)/(R2-R1)*(r(IREND:)-R1)
        end if
!       == fulfill outer boundary condition =====================================
        CALL RADIAL$VALUE(GID,NR,PHIHOM,Rout,VAL1)
        CALL RADIAL$VALUE(GID,NR,PHI1,Rout,VAL2)
        SVAR=VAL1+VAL2
        VAL1=VAL1/SVAR
        VAL2=VAL2/SVAR
        PHIHOM(:)=VAL2*PHIHOM(:)-VAL1*PHI1(:)
        phihom(:)=phihom(:)/phihom(irmatch)
!       == inhomogeneous solution with correct boundary conditions
        IF(.NOT.THOM) THEN     
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
          CALL RADIAL$VALUE(GID,NR,PHIINHOM,ROUT,VAL1)
          CALL RADIAL$VALUE(GID,NR,PHI1,ROUT,VAL2)
          PHIINHOM(:)=PHIINHOM(:)-VAL1/VAL2*PHI1(:)
         ELSE
          PHIINHOM(:)=0.D0
        END IF
!
!       =======================================================================
!       ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!       =======================================================================
        SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
        PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
        CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
        CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
        SVAR=(DERO-DER)/PHI(IRMATCH)
        PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
      END IF
!
!     =======================================================================
!     ==  SET WAVE FUNCTION TO ZERO BEYOND RBOX                            ==
!     =======================================================================
DO IR=1,NR
  IF(.NOT.(PHI(IR).GT.0.D0.OR.PHI(IR).LE.0.D0)) THEN
    PRINT*,'ERROR'
!    PRINT*,'PHIIN',PHI(:IRMATCH-1)
!    PRINT*,'PHIOUT',PHI(IRMATCH:)
    PRINT*,'SVAR ',SVAR
    PRINT*,'IROUT,IRCL ',IROUT,IRCL,NR
    PRINT*,'R,PHI ',R(IR),PHI(IR)
    OPEN(UNIT=8,FILE='XXX.DAT')
    DO I=1,NR
      WRITE(8,*)I,R(I),PHIHOM(I),POT1(I)*Y0
    ENDDO
    CLOSE(8)      
    CALL ERROR$MSG('PHI CONTAINS NANS')
    CALL ERROR$STOP('BOUNDSTATE')
  END IF
ENDDO
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE BOUNDSTATEX(GID,NR,L,SO,RBOX,DREL,G,NN,POT,E,PHI)
!     **                                                                  **
!     **  FINDS A BOUND STATE OF THE RADIAL SCHROEDINGER EQUATION AND     **
!     **  ITS ENERGY FOR A  SPECIFIED NUMBER OF NODES NN.                 **
!     **  SCHROEDINGER EQUATION MAY INVOLVE AN INHOMOGENEITY G            **
!     **                                                                  **
!     **  FIRST, THE ENERGY IS DETERMINED BY BISECTION ON THE             **
!     **  GENERALIZED PHASE SHIFT AT THE OUTERMOST RADIAL GRID POINT.     **
!     **  THIS WAVE FUNCTION MAY HOWEVER STILL DIVERGE EXPONENTIALLY.     **
!     **                                                                  **
!     **  SECONDLY, THE SCHROEDINGER EQUATION IS SOLVED INWARD, AND       **
!     **  MATCHED WITH VALUE, EITHER AT THE CLASSICAL TURNING POINT       **
!     **  OR A SPECIFIED RADIUS, WHATEVER IS SMALLER.                     **
!     **                                                                  **
!     **                                                                  **
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! ANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
      REAL(8)    ,INTENT(IN)     :: RBOX    ! BOX RADIUS
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      INTEGER(4) ,INTENT(IN)     :: NN      !#(NODES)
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL
      REAL(8)    ,INTENT(INOUT)  :: E       !ENERGY
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      INTEGER(4)                 :: ISTART,IBI
      REAL(8)                    :: X0,DX,XM,ZM,Z0
      REAL(8)    ,PARAMETER      :: TOL=1.D-12
      REAL(8)                    :: PI,Y0
      REAL(8)                    :: DER,DERO
      REAL(8)                    :: R(NR)
      REAL(8)                    :: PHIHOM(NR),PHIHOM2(NR),PHIINHOM(NR),GHOM(NR)
      INTEGER(4) ,PARAMETER      :: NITER=100
      INTEGER(4)                 :: I,IR
      INTEGER(4)                 :: IDIR ! SWITCH FOR OUT/INWARD INTEGRATION 
      INTEGER(4)                 :: IRMATCH 
      REAL(8)                    :: SVAR
      REAL(8)                    :: POT1(NR)
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20 ! MAXIMUM FACTOR IN THE WAVE FUNCTION
      REAL(8)   ,PARAMETER       :: EMAX=100.D0 ! MAXIMUM ENERGY
      LOGICAL(4)                 :: THOM
      REAL(8)                    :: ROUT
      REAL(8)                    :: VAL,VAL1,VAL2
      INTEGER(4)                 :: IROUT,IRCL,IRBOX
!     *********************************************************************
                                 CALL TRACE$PUSH('BOUNDSTATE')
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      IRBOX=1
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GT.RBOX) EXIT
      ENDDO
!          
      ISTART=1
      X0=E
      DX=1.D-2
      CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      DO I=1,NITER
        E=X0
!        E=MIN(X0,EMAX)
!
!       =======================================================================
!       ==  CUT OFF THE POTENTIAL                                            ==
!       =======================================================================
        CALL SCHROEDINGER_SPECIALRADS(GID,NR,L,XMAX,POT,E,IRCL,IROUT)
        ROUT=MIN(R(IROUT),RBOX)
        IROUT=MIN(IROUT,IRBOX)
        POT1(:)=POT(:)
        POT1(IROUT:)=POT(IROUT)
!
!       =======================================================================
!       == INTEGRATE RADIAL SCHROEDINGER EQUATION OUTWARD                     ==
!       =======================================================================
        IDIR=1
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHI)
!       == CHECK FOR OVERFLOW
        IF(.NOT.(PHI(IROUT).GT.0.OR.PHI(IROUT).LE.0)) THEN
          CALL ERROR$MSG('OVERFLOW AFTER OUTWARD INTEGRATION')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$I4VAL('NN',NN)
          CALL ERROR$STOP('BOUNDSTATE')
        END IF
!
!       =======================================================================
!       == ESTIMATE PHASE SHIFT                                              ==
!       =======================================================================
        CALL SCHROEDINGER$PHASESHIFT(GID,NR,PHI,ROUT,Z0)
        Z0=Z0-REAL(NN+1)
        IF(ABS(2.D0*DX).LE.TOL) EXIT
!       =====================================================================
!       ==  BISECTION                                                      ==
!       =====================================================================
        CALL BISEC(ISTART,IBI,X0,Z0,DX,XM,ZM)
      ENDDO
      IF(ABS(DX).GT.TOL) THEN
        CALL ERROR$MSG('BISECTION LOOP NOT CONVERGED')
        CALL ERROR$MSG('BOUND STATE NOT FOUND')
        CALL ERROR$STOP('BOUNDSTATE')
      END IF
!
!     =======================================================================
!     ==  DETERMINE MATCHING POINT                                         ==
!     =======================================================================
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
      IRMATCH=IRCL
      IF(R(IRCL).GT.5.D0) THEN
        CALL RADIAL$XOFR(GID,5.D0,SVAR)
        IRMATCH=INT(SVAR)
      END IF
!
!     =======================================================================
!     ==  INTEGRATE INWARD                                                 ==
!     =======================================================================
      IF(IRMATCH.LT.IROUT) THEN
        IDIR=-1
        GHOM(:)=0.D0
        IR=MIN(IROUT+3,NR)
        GHOM(IR)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM)
        GHOM=0.D0
        IR=MIN(IROUT+10,NR)
        GHOM(IR)=1.D-8
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,GHOM,L,E,IDIR,PHIHOM2)
        CALL RADIAL$VALUE(GID,NR,PHIHOM,RBOX,VAL1)
        CALL RADIAL$VALUE(GID,NR,PHIHOM2,RBOX,VAL2)
        IF(VAL2.NE.0.D0) THEN
          SVAR=-VAL1/VAL2
          PHIHOM(:)=PHIHOM(:)+SVAR*PHIHOM2(:)
        END IF
        IF(.NOT.THOM) THEN     
          CALL SCHROEDINGER$SPHERICAL(GID,NR,POT1,DREL,SO,G,L,E,IDIR,PHIINHOM)
          IF(VAL2.NE.0.D0) THEN
            CALL RADIAL$VALUE(GID,NR,PHIINHOM,RBOX,VAL1)
            SVAR=-VAL1/VAL2
            PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM2(:)
          END IF
        ELSE
          PHIINHOM(:)=0.D0
        END IF
!
!       =======================================================================
!       ==  MATCH SOLUTION INSIDE AND OUTSIDE WITH VALUE                     ==
!       =======================================================================
        SVAR=(PHI(IRMATCH)-PHIINHOM(IRMATCH))/PHIHOM(IRMATCH)
        PHIINHOM(:)=PHIINHOM(:)+SVAR*PHIHOM(:)
        CALL RADIAL$DERIVATIVE(GID,NR,PHI,R(IRMATCH),DER)
        CALL RADIAL$DERIVATIVE(GID,NR,PHIINHOM,R(IRMATCH),DERO)
        SVAR=(DERO-DER)/PHI(IRMATCH)
        PHI(IRMATCH:)=PHIINHOM(IRMATCH:)
      END IF
!
!     =======================================================================
!     ==  SET WAVE FUNCTION TO ZERO BEYOND RBOX                            ==
!     =======================================================================
!      PHI(IROUT:)=0.D0
!
DO IR=1,NR
  IF(.NOT.(PHI(IR).GT.0.D0.OR.PHI(IR).LE.0.D0)) THEN
    PRINT*,'ERROR'
!    PRINT*,'PHIIN',PHI(:IRMATCH-1)
!    PRINT*,'PHIOUT',PHI(IRMATCH:)
    PRINT*,'SVAR ',SVAR
    PRINT*,'IROUT,IRCL ',IROUT,IRCL,NR
    PRINT*,'R,PHI ',R(IR),PHI(IR)
    OPEN(UNIT=8,FILE='XXX.DAT')
    DO I=1,NR
      WRITE(8,*)I,R(I),PHIHOM(I),POT1(I)*Y0
    ENDDO
    CLOSE(8)      
    CALL ERROR$MSG('PHI CONTAINS NANS')
    CALL ERROR$STOP('BOUNDSTATE')
  END IF
ENDDO
                                 CALL TRACE$POP()
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE GRIDPARAMETERS(DMIN,DMAX,RX,R1,DEX,NR)
!     **                                                               **
!     **  DETERMINES THE GRID PARAMETERS FOR THE SHIFTED LOGARITHMIC   **
!     **  GRID FROM A SPECIFIED MINIMUM SPACING DMIN, A MAXIMUM        **
!     **  SPACING DMAX AND A MAXIMUM RADIUS RX                         **
!     **                                                               **
      REAL(8),   INTENT(IN) :: DMIN
      REAL(8),   INTENT(IN) :: DMAX
      REAL(8),   INTENT(IN) :: RX
      REAL(8),   INTENT(OUT):: R1
      REAL(8),   INTENT(OUT):: DEX
      INTEGER(4),INTENT(OUT):: NR
      REAL(8)               :: RN
      REAL(8)               :: Q   ! EXP(DEX)
!     *******************************************************************
      RN=2.D0+LOG(DMAX/DMIN)/LOG((RX-DMIN)/(RX-DMAX))
      Q=(DMAX/DMIN)**(1.D0/(RN-2.D0))
      DEX=LOG(Q)
      R1=DMIN/(Q-1)
!
!      PRINT*,'DMIN ',DMIN,R1*(EXP(DEX)-1.D0)            
!      PRINT*,'DMAX ',DMAX,R1*(EXP(DEX*REAL(RN-1))-EXP(DEX*REAL(RN-2)))
!      PRINT*,'RX   ',RX,R1*(EXP(DEX*REAL(RN-1))-1.D0)
      NR=NINT(RN)
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE AEETOT(GID,NR,NB,Z,NC,EOFI,FOFI,AERHOC,AERHOV,AEPOT)
!     **                                                               **
!     **                                                               **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NB
      INTEGER(4),INTENT(IN) :: NC
      REAL(8),   INTENT(IN) :: Z
      REAL(8),   INTENT(IN) :: AERHOC(NR)
      REAL(8),   INTENT(IN) :: AERHOV(NR)
      REAL(8)   ,INTENT(IN) :: EOFI(NB)
      REAL(8)   ,INTENT(IN) :: FOFI(NB)
      REAL(8)   ,INTENT(IN) :: AEPOT(NR)
      REAL(8)               :: EKINC,EKINV
      REAL(8)               :: EHC,EHV
      REAL(8)               :: EXCC,EXCV
      REAL(8)               :: R(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR
      INTEGER(4)            :: NFILO
!     *******************************************************************
      CALL RADIAL$R(GID,NR,R)
!     == CORE KINETIC ENERGY ============================================
      AUX(:)=AERHOC(:)*AEPOT(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      EKINC=SUM(EOFI(1:NC)*FOFI(1:NC))-SVAR
!     == VALENCE KINETCI ENERGY ==========================================
      AUX=AERHOV(:)*AEPOT(:)*R(:)**2
      CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
      EKINV=SUM(EOFI(NC+1:)*FOFI(NC+1:))-SVAR
!     == CORE POTENTIAL ENERGY  ==========================================
      CALL MYVOFRHO(GID,NR,Z,AERHOC,AUX,EHC,EXCC)
!     == CORE POTENTIAL ENERGY  ==========================================
      CALL MYVOFRHO(GID,NR,Z,AERHOC+AERHOV,AUX,EHV,EXCV)
      EHV=EHV-EHC
      EXCV=EXCV-EXCC
!     == WRITE ===========================================================
      CALL FILEHANDLER$UNIT('PROT',NFILO)
      CALL REPORT$TITLE(NFILO,'AE-ENERGIES')
      CALL REPORT$R8VAL(NFILO,'ETOT(VALENCE ONLY)',EKINV+EHV+EXCV,'H')
      CALL REPORT$R8VAL(NFILO,'ETOT',EKINV+EKINC+EHV+EHC+EXCV+EXCC,'H')
      CALL REPORT$R8VAL(NFILO,'E-BAND',SUM(EOFI(NC+1:)*FOFI(NC+1:)),'H')
      CALL REPORT$R8VAL(NFILO,'E-KIN',EKINV,'H')
      CALL REPORT$R8VAL(NFILO,'E-HARTREE',EHV,'H')
      CALL REPORT$R8VAL(NFILO,'E-XC',EXCV,'H')
      CALL REPORT$R8VAL(NFILO,'E-XC-CORE',EXCC,'H')
      CALL REPORT$R8VAL(NFILO,'E-XC-VALENCE',EXCV,'H')
      CALL REPORT$R8VAL(NFILO,'E-CORE',EKINC+EHC+EXCC,'H')
      RETURN
      END


!!$.. this is now in paw_generalpurpose and therefore in the paw-library
!!$!..................................................................................
!!$MODULE BROYDEN_MODULE
!!$LOGICAL(4)         :: TON=.FALSE.
!!$INTEGER(4)         :: NSTEPX=0
!!$INTEGER(4)         :: NSTEP=0
!!$INTEGER(4)         :: NX=0
!!$REAL(8)            :: ALPHA
!!$REAL(8),ALLOCATABLE :: XPREV(:,:)
!!$REAL(8),ALLOCATABLE :: YPREV(:,:)
!!$END MODULE BROYDEN_MODULE
!!$!      .............................................................................
!!$       SUBROUTINE BROYDEN$NEW(NX_,NSTEPX_,ALPHA_)
!!$       USE BROYDEN_MODULE
!!$       IMPLICIT NONE
!!$       INTEGER(4),INTENT(IN)    :: NX_
!!$       INTEGER(4),INTENT(IN)    :: NSTEPX_
!!$       REAL(8)   ,INTENT(IN)    :: ALPHA_
!!$!      *****************************************************************************
!!$       IF(TON) THEN
!!$         CALL ERROR$MSG('BROYDEN OBJECT ALREADY IN USE')
!!$         CALL ERROR$STOP('BROYDEN$NEW')
!!$       END IF
!!$       TON=.TRUE.
!!$       NSTEP=0
!!$       NX=NX_
!!$       NSTEPX=NSTEPX_
!!$       ALPHA=ALPHA_
!!$       ALLOCATE(XPREV(NX,NSTEPX))
!!$       ALLOCATE(YPREV(NX,NSTEPX))
!!$       RETURN
!!$       END
!!$!      .............................................................................
!!$       SUBROUTINE BROYDEN$CLEAR
!!$       USE BROYDEN_MODULE
!!$       IMPLICIT NONE
!!$!      *****************************************************************************
!!$       IF(.NOT.TON) THEN
!!$         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
!!$         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
!!$         CALL ERROR$STOP('BROYDEN$CLEAR')
!!$       END IF
!!$       TON=.FALSE.
!!$       NSTEPX=0
!!$       NSTEP=0
!!$       NX=0
!!$       ALPHA=0.D0
!!$       DEALLOCATE(XPREV)
!!$       DEALLOCATE(YPREV)
!!$       RETURN
!!$       END
!!$!      .............................................................................
!!$       SUBROUTINE BROYDEN$STEP(NX_,X,Y)
!!$       USE BROYDEN_MODULE
!!$       IMPLICIT NONE
!!$       INTEGER(4),INTENT(IN)    :: NX_
!!$       REAL(8)   ,INTENT(INOUT) :: X(NX_)
!!$       REAL(8)   ,INTENT(IN)    :: Y(NX_)
!!$       REAL(8)   ,ALLOCATABLE   :: DX(:,:)
!!$       REAL(8)   ,ALLOCATABLE   :: DY(:,:)
!!$       REAL(8)   ,ALLOCATABLE   :: B(:,:)
!!$       REAL(8)   ,ALLOCATABLE   :: BINV(:,:)
!!$ REAL(8)   ,ALLOCATABLE   :: W(:,:)
!!$       INTEGER(4)               :: I
!!$!      *****************************************************************************
!!$       IF(.NOT.TON) THEN
!!$         CALL ERROR$MSG('BROYDEN OBJECT NOT ACTIVE')
!!$         CALL ERROR$MSG('CALL BROYDEN$NEW FIRST')
!!$         CALL ERROR$STOP('BROYDEN$STEP')
!!$       END IF
!!$       IF(NX_.NE.NX) THEN
!!$         CALL ERROR$MSG('SIZE INCONSISTENT')
!!$         CALL ERROR$STOP('BROYDEN$STEP')
!!$       END IF
!!$!PRINT*,'NSTEP',NSTEP
!!$!
!!$!      =================================================================
!!$!      == SIMPLE MIXING IN THE FIRST STEP                             ==
!!$!      =================================================================
!!$       IF(NSTEP.EQ.0) THEN
!!$         IF(NSTEPX.GT.0)THEN
!!$           NSTEP=1
!!$           XPREV(:,1)=X(:)     
!!$           YPREV(:,1)=Y(:)     
!!$         END IF
!!$         X=X+ALPHA*Y
!!$         RETURN
!!$       END IF
!!$!
!!$!      =================================================================
!!$!      == DETERMINE INVERSE HESSIAN ALPHA+DX OTIMES DY                ==
!!$!      =================================================================
!!$       ALLOCATE(DX(NX,NSTEP))
!!$       ALLOCATE(DY(NX,NSTEP))
!!$       DO I=1,NSTEP
!!$         DY(:,I)=YPREV(:,I)-Y(:)  
!!$         DX(:,I)=XPREV(:,I)-X(:)+ALPHA*DY(:,I)
!!$       ENDDO
!!$       ALLOCATE(B(NSTEP,NSTEP))
!!$       ALLOCATE(BINV(NSTEP,NSTEP))
!!$       B=MATMUL(TRANSPOSE(DY),DY)   !OVERLAP MATRIX OF DY
!!$!PRINT*,'B',B
!!$       CALL LIB$INVERTR8(NSTEP,B,BINV)           
!!$ALLOCATE(W(NX,NSTEP))
!!$W=MATMUL(DY,BINV)            !NEW DY IS BIORTHONORMAL TO OLD DY
!!$!PRINT*,'W ',MATMUL(TRANSPOSE(W),DY)
!!$DEALLOCATE(W)
!!$       DY=MATMUL(DY,BINV)            !NEW DY IS BIORTHONORMAL TO OLD DY
!!$       DEALLOCATE(B)
!!$       DEALLOCATE(BINV)
!!$!
!!$!      =================================================================
!!$!      == STORE HISTORY                                               ==
!!$!      =================================================================
!!$       IF(NSTEP.LT.NSTEPX)NSTEP=NSTEP+1
!!$       DO I=NSTEP,2,-1
!!$         YPREV(:,I)=YPREV(:,I-1)
!!$         XPREV(:,I)=XPREV(:,I-1)
!!$       ENDDO
!!$       XPREV(:,1)=X(:)     
!!$       YPREV(:,1)=Y(:)     
!!$!
!!$!      =================================================================
!!$!      == PREDICT NEW VECTOR                                          ==
!!$!      =================================================================
!!$       X=X+ALPHA*Y-MATMUL(DX,MATMUL(TRANSPOSE(DY),Y))
!!$       DEALLOCATE(DX)
!!$       DEALLOCATE(DY)
!!$       RETURN
!!$       END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONEATOM$NEW(IAT,AEZ)
!     **************************************************************************
!     ** INITIALIZES AN ATOM                                                  **
!     **************************************************************************
      USE ONEATOM_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: IAT
      REAL(8)   ,INTENT(IN)  :: AEZ          ! ATOMIC NUMBER
      REAL(8)   ,PARAMETER   :: R1=1.056D-4  ! INNERMOST POINT OF RADIAL GRID
      REAL(8)   ,PARAMETER   :: DEX=0.05D0   ! LOG SPACING OF RADIAL GRID
      INTEGER(4),PARAMETER   :: NR=250       ! #(RADIAL GRID POINTS)
      INTEGER(4)             :: GID          ! GRID ID
      INTEGER(4)             :: NS(10),NP(10),ND(10),NF(10)
      INTEGER(4)             :: SOFI(19)
      INTEGER(4)             :: LOFI(19)
      INTEGER(4)             :: NNOFI(19)
      REAL(8)                :: FOFI(19)
      REAL(8)                :: EOFI(19)
      REAL(8)                :: RBOX=3.D0
      REAL(8)                :: EV
      REAL(8)                :: SVAR
      INTEGER(4)             :: NB       ! #(STATES IN OCCUPIED SHELLS)
      INTEGER(4)             :: NC       ! #(CORE STATES)
      INTEGER(4)             :: L,IB,ISVAR
      REAL(8)   ,ALLOCATABLE :: AEPOT(:)
      REAL(8)   ,ALLOCATABLE :: DREL(:)
      REAL(8)   ,ALLOCATABLE :: RHOADD(:)
      REAL(8)   ,ALLOCATABLE :: RHOCORE(:)
      REAL(8)   ,ALLOCATABLE :: AUX(:)
      REAL(8)   ,ALLOCATABLE :: PHI(:)
      REAL(8)   ,ALLOCATABLE :: R(:)
      REAL(8)                :: ETOT,DEDQ,DEDRAD,EKINC
      REAL(8)   ,PARAMETER   :: Q=0.D0
      REAL(8)                :: RAD
      REAL(8)                :: PI,C0LL
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      C0LL=1.D0/SQRT(4.D0*PI)
      CALL CONSTANTS('EV',EV)
      IF(IAT.GT.NATX) THEN
        CALL ERROR$MSG('ATOM INDEX OUT OF RANGE')
        CALL ERROR$STOP('ONEATOM$NEW')
      END IF
      THIS=>THISARR(IAT)
      THIS%AEZ=AEZ
!
!     ==========================================================================
!     ==  DEFINE RADIAL GRID                                                  ==
!     ==========================================================================
      CALL RADIAL$NEW('SHLOG',GID)
      THIS%GID=GID
      CALL RADIAL$SETR8(GID,'R1',R1)
      CALL RADIAL$SETR8(GID,'DEX',DEX)
      CALL RADIAL$SETI4(GID,'NR',NR)
      THIS%NR=NR
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      RAD=R(NR-3)
!
!     ==========================================================================
!     ==  DEFINE BANDS                                                        ==
!     ==========================================================================
      LOFI( 1)=0   ! 1S
      LOFI( 2)=0   ! 2S
      LOFI( 3)=1   ! 2P
      LOFI( 4)=0   ! 3S
      LOFI( 5)=1   ! 3P
      LOFI( 6)=0   ! 4S
      LOFI( 7)=2   ! 3D
      LOFI( 8)=1   ! 4P
      LOFI( 9)=0   ! 5S
      LOFI(10)=2   ! 4D
      LOFI(11)=1   ! 5P
      LOFI(12)=0   ! 6S
      LOFI(13)=3   ! 4F
      LOFI(14)=2   ! 5D
      LOFI(15)=1   ! 6P
      LOFI(16)=0   ! 7S
      LOFI(17)=3   ! 5F
      LOFI(18)=2   ! 6D
      LOFI(19)=1   ! 7P
!
!     == DETERMINE NUMBER OF NODES ===================================
      DO L=0,3
        ISVAR=0
        DO IB=1,19
          IF(LOFI(IB).NE.L) CYCLE
          NNOFI(IB)=ISVAR
          ISVAR=ISVAR+1
        ENDDO
      ENDDO
!
!     == OCCUPATIONS  AND NUMBER OF BANDS===============================
      FOFI(:)=0.D0
      SVAR=AEZ
      DO IB=1,19
        FOFI(IB)=MIN(SVAR,REAL(2*(2*LOFI(IB)+1),KIND=8))
        SVAR=SVAR-FOFI(IB)
      ENDDO
!
!     == DETERMINE THE NUMBER OF CORE STATES =============================
      NB=19
      DO IB=1,19
        IF(LOFI(IB).NE.0) CYCLE
        IF(FOFI(IB).NE.0.D0) NC=IB-1
        IF(FOFI(IB).EQ.0.D0) THEN
          NB=IB-1   ! include at least one empty orbital to facilitate fermi level search
          EXIT
        END IF
      ENDDO
!
!     == MAP ON THIS ==================================================
      THIS%NB=NB
      THIS%NC=NC
      THIS%LOFI=LOFI
      THIS%NNOFI=NNOFI
      THIS%FOFI=FOFI
!
!     ==========================================================================
!     ==  DETERMINE ATOMIC POTENTIAL                                          ==
!     ==========================================================================
print*,'aez',aez
print*,'core'
do ib=1,nc
  print*,ib,lofi(ib),nnofi(ib),fofi(ib)
enddo
print*,'valence'
do ib=nc+1,nb
  print*,ib,lofi(ib),nnofi(ib),fofi(ib)
enddo
      ALLOCATE(AEPOT(NR))
      CALL RADIAL$NUCPOT(GID,NR,AEZ,AEPOT)
      ALLOCATE(DREL(NR))
      DREL(:)=0.D0
      ALLOCATE(RHOADD(NR))
      RHOADD=0.D0
      SOFI(:)=0
      EOFI(:)=0.D0
      CALL AESCF(GID,NR,NB,TREL,LOFI,SOFI,FOFI,NNOFI,AEZ,RHOADD,RBOX,DREL,AEPOT,EOFI)
!
!     == CORE DENSITY
      ALLOCATE(PHI(NR))
      ALLOCATE(AUX(NR))
      ALLOCATE(RHOCORE(NR))
      RHOCORE(:)=0.D0
      DO IB=1,NC
        AUX(:)=0.D0 !INHOMOGENITY
        CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RBOX,DREL,AUX,NNOFI(IB),AEPOT,EOFI(IB),PHI)
        AUX(:)=(R(:)*PHI(:))**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        PHI(:)=PHI(:)/SQRT(SVAR)
        RHOCORE(:)=RHOCORE(:)+FOFI(IB)*C0LL*PHI(:)**2
      ENDDO
      DEALLOCATE(PHI)
      DEALLOCATE(AUX)
!
!     == CORE KINETIC ENERGY ====================================================
      ALLOCATE(AUX(NR))
      AUX(:)=R(:)**2*AEPOT(:)*RHOCORE(:)
      CALL RADIAL$INTEGRAL(GID,NR,AUX,EKINC)
      EKINC=-EKINC
      DO IB=1,NC
        EKINC=EKINC+FOFI(IB)*EOFI(IB)
      ENDDO
      DEALLOCATE(AUX)
!
!     == MAP ONTO THIS
      THIS%EKINC=EKINC
      THIS%EOFI=1.D0
      THIS%EOFI(1:NB)=EOFI(1:NB)
      ALLOCATE(THIS%AEPOT(NR))
      THIS%AEPOT(:)=AEPOT(:)
      DEALLOCATE(AEPOT)
      ALLOCATE(THIS%RHOCORE(NR))
      THIS%RHOCORE(:)=RHOCORE(:)
      DEALLOCATE(RHOCORE)
      ALLOCATE(THIS%DREL(NR))
      THIS%DREL=DREL
      DEALLOCATE(DREL)
!
!     == REPORT ================================================================
      WRITE(*,FMT='("EIGENSTATES AFTER NON-SCF CALCULATION IN SCF POT")')
      DO IB=1,NB
        WRITE(*,FMT='("L=",I2," NN=",I2," SO=",I2," F=",F5.2 &
     &                ," E[H]=",F20.9," E[EV]=",F20.9)') &
     &              LOFI(IB),NNOFI(IB),SOFI(IB),FOFI(IB),EOFI(IB),EOFI(IB)/EV
      ENDDO
!
      THIS%EREF=0.D0
      ALLOCATE(this%vpauli(NR,MAXVAL(LOFI(1:nb))+1))
      ALLOCATE(this%UDOT(NR,MAXVAL(LOFI(1:nb))+1))
      ALLOCATE(this%UN(NR,MAXVAL(LOFI(1:nb))+1))
      CALL ONEATOM(IAT,Q,RAD,ETOT,DEDQ,DEDRAD)
      THIS%EREF=ETOT
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONEATOM(IAT,Q,RAD,ETOT,DEDQ,DEDRAD)
!     **************************************************************************
!     **  THE TEMPERATURE MUST NOT BE TOO SMALL. OTHERWISE DEDQ CANNOT BE     **
!     **  ACCURATELY DETERMINED                                               **
!     **************************************************************************
      USE ONEATOM_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: IAT
      REAL(8)   ,INTENT(IN) :: Q
      REAL(8)   ,INTENT(IN) :: RAD
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,INTENT(OUT):: DEDQ
      REAL(8)   ,INTENT(OUT):: DEDRAD
      INTEGER(4)            :: GID
      INTEGER(4)            :: NR
      INTEGER(4)            :: NB
      INTEGER(4)            :: NC
      REAL(8)               :: AEZ
      REAL(8)               :: QVALENCE
      INTEGER(4)            :: LOFI(19)
      INTEGER(4)            :: SOFI(19)
      INTEGER(4)            :: NNOFI(19)
      REAL(8)               :: EOFI(19)
      REAL(8)               :: FOFI(19)
      REAL(8)   ,ALLOCATABLE:: AEPOT(:)
      REAL(8)   ,ALLOCATABLE:: POTIN(:)
      REAL(8)   ,ALLOCATABLE:: RHOADD(:)
      REAL(8)   ,ALLOCATABLE:: RHO(:)
      REAL(8)   ,ALLOCATABLE:: DREL(:)
      REAL(8)   ,ALLOCATABLE:: PHI(:,:)
      REAL(8)   ,ALLOCATABLE:: PHIDOT(:)
      REAL(8)   ,ALLOCATABLE:: G(:)
      REAL(8)   ,ALLOCATABLE:: AUX(:),AUX1(:)
      REAL(8)   ,ALLOCATABLE:: R(:)
      REAL(8)   ,ALLOCATABLE:: vpauli(:,:)
      REAL(8)   ,ALLOCATABLE:: UDOT(:,:)
      REAL(8)   ,ALLOCATABLE:: UN(:,:)
      REAL(8)               :: SVAR,DER,VAL,DERDOT,VALDOT
      REAL(8)               :: PI,C0LL
      REAL(8)               :: EKIN,EH,EXC
      REAL(8)               :: XMAX
      LOGICAL(4)            :: CONVG
      INTEGER(4)            :: IB,IR,ITER,IRBOX,IL,L,LMAX
      INTEGER(4),PARAMETER  :: NITER=100
      REAL(8)   ,PARAMETER  :: TOL=1.D-5
      REAL(8)   ,PARAMETER  :: KBT=1.D-2
      REAL(8)               :: TS,F
INTEGER(4) :: I,J
REAL(8)    :: DEG,SVAR1
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      C0LL=1.D0/SQRT(4.D0*PI)
      THIS=>THISARR(IAT)
      GID=THIS%GID
      NR=THIS%NR
      AEZ=THIS%AEZ      
      NB=THIS%NB
      NC=THIS%NC
      LOFI(:)=THIS%LOFI(:)
      NNOFI(:)=THIS%NNOFI(:)
      EOFI(:)=THIS%EOFI(:)
      FOFI(:)=THIS%FOFI(:)
      ALLOCATE(R(NR))
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GT.RAD) EXIT
      ENDDO
!
!     ==========================================================================
!     ==  OBTAIN SELF-CONSISTENT POTENTIAL                                    ==
!     ==========================================================================
      SOFI(:)=0.D0
      ALLOCATE(AEPOT(NR));       AEPOT(:)=THIS%AEPOT(:)
      ALLOCATE(POTIN(NR));       POTIN(:)=AEPOT(:)
      ALLOCATE(RHOADD(NR));      RHOADD=0.D0
      ALLOCATE(DREL(NR));        DREL=THIS%DREL
      ALLOCATE(RHO(NR))
      ALLOCATE(PHI(NR,19))
      ALLOCATE(PHIDOT(NR))
      ALLOCATE(G(NR))
      ALLOCATE(AUX(NR))
      ALLOCATE(AUX1(NR))
      XMAX=0.D0
      CONVG=.FALSE.
      CALL BROYDEN$NEW(NR,4,1.D0)
      POTIN=AEPOT
      DO ITER=1,NITER
        DO IB=NC+1,NB
          G(:)=0.D0
          CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RAD,DREL,G,NNOFI(IB),POTIN,EOFI(IB),PHI(:,IB))
          AUX(:)=(R(:)*PHI(:,IB))**2
          CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
          CALL RADIAL$VALUE(GID,NR,AUX1(:),RAD,SVAR)
          PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
        ENDDO
!
!       == DETERMINE OCCUPATIONS ===============================================
        QVALENCE=AEZ+Q-REAL(SUM(2*(2*LOFI(1:NC)+1)),KIND=8)
        CALL ONEATOM_OPTF(NB-NC,QVALENCE,KBT,LOFI(NC+1:NB),EOFI(NC+1:NB),FOFI(NC+1:NB),DEDQ)
!
!       == DETERMINE DENSITY ===================================================
        RHO(:)=THIS%RHOCORE(:)
        DO IB=NC+1,NB
          RHO(:)=RHO(:)+FOFI(IB)*C0LL*PHI(:,IB)**2
        ENDDO
!
!       == DETERMINE OUTPUT POTENTIAL ==========================================
        CALL MYVOFRHOWITHRAD(GID,NR,RAD,AEZ,RHO,AEPOT,EH,EXC)
!
!       == EXIT IF CONVERGED ===================================================
        IF(CONVG) THEN
          CALL BROYDEN$CLEAR
          EXIT
        END IF
!
!       == POTENTIAL MIXING ====================================================
!       == THE POTENTIAL SIS SET EQUAL BEYOND THE BOX RADIUS, BECAUSE OTHERWISE
!       == THIS REGION DETERMINES THE CONVERGENCE ==============================
        POTIN(IRBOX:)=AEPOT(IRBOX:)
        XMAX=MAXVAL(ABS(AEPOT-POTIN)) 
!PRINT*,'XMAX ',XMAX,MAXLOC(ABS(AEPOT-POTIN))
        CALL BROYDEN$STEP(NR,POTIN,(AEPOT-POTIN))
        CONVG=(XMAX.LT.TOL)
      ENDDO
      IF(.NOT.CONVG) THEN
        CALL BROYDEN$CLEAR
PRINT*,'LOOP IN ONEATOM NOT CONVERGED!:',XMAX,MAXLOC(ABS(AEPOT-POTIN)),IAT
        CALL ERROR$MSG('SELFCONSISTENCY LOOP NOT CONVERGED')
        CALL ERROR$STOP('ONEATOM')
      END IF
OPEN(10,FILE='XX.DAT') 
REWIND(10)
DO IR=1,NR
  WRITE(10,FMT='(10E20.5)')R(IR),AEPOT(IR),POTIN(IR),RHO(IR),DREL(IR),PHI(IR,NC+1:NB)
ENDDO
CLOSE(10)
!
!     ==========================================================================
!     == CALCULATE PRESSURE                                                   ==
!     ==========================================================================
      DEDRAD=0.D0
      DO IB=NC+1,NB
        IF(FOFI(IB).EQ.0.D0) CYCLE
        CALL RADIAL$DERIVATIVE(GID,NR,PHI(:,IB),RAD,DER)
        IF(DER.EQ.0.D0) CYCLE  ! AVOID CORE WAVE FUNCTIONS
        G(:)=0.D0
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,SOFI(IB),G,LOFI(IB),EOFI(IB),1,PHI(:,IB))
! NORMALIZATION NOT REQUIRED
 AUX=R (:)**2*PHI(:,IB)**2
 CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
 CALL RADIAL$VALUE(GID,NR,AUX1(:),RAD,SVAR)
 PHI(:,IB)=PHI(:,IB)/SQRT(SVAR)
!
        CALL RADIAL$VALUE(GID,NR,PHI(:,IB),RAD,VAL)
        CALL RADIAL$DERIVATIVE(GID,NR,PHI(:,IB),RAD,DER)
        G(:)=PHI(:,IB)         
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,SOFI(IB),G,LOFI(IB),EOFI(IB),1,PHIDOT)
! == ORTHONORMALIZATION IS NOT NECESSARY =======
  AUX=R (:)**2*PHI(:,IB)*PHIDOT(:)
  CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
  CALL RADIAL$VALUE(GID,NR,AUX1(:),RAD,SVAR)
  PHIDOT(:)=PHIDOT(:)-PHI(:,IB)*SVAR
        CALL RADIAL$VALUE(GID,NR,PHIDOT,RAD,VALDOT)
        CALL RADIAL$DERIVATIVE(GID,NR,PHIDOT,RAD,DERDOT)
        DEDRAD=DEDRAD-FOFI(IB)*DER/VALDOT
      ENDDO
!
!     ==========================================================================
!     == CALCULATE TOTAL ENERGY                                               ==
!     ==========================================================================
      AUX(:)=(RHO(:)-THIS%RHOCORE(:))*AEPOT(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,SVAR)
      EKIN=THIS%EKINC+SUM(EOFI*FOFI)-SVAR
      CALL MYVOFRHOWITHRAD(GID,NR,RAD,AEZ,RHO,AUX,EH,EXC)
      ETOT=EKIN+EH+EXC-THIS%EREF
!     == ENTROPY TERM OF THE ELECTRONS
      TS=0.D0
      DO IB=NC+1,NB
        SVAR=REAL(2*(2*LOFI(IB)+1),KIND=8)
        F=FOFI(IB)/SVAR
        IF(F.LT.1.D-50) CYCLE
        IF(1.D0-F.LT.1.D-50) CYCLE
        TS=TS+KBT*SVAR*(F*LOG(F)+(1.D0-F)*LOG(1.D0-F))
      ENDDO   
      ETOT=ETOT+TS   
!
!     ==========================================================================
!     == CALCULATE REPULSIVE POTENTIAL FOR EACH ANGULAR MOMENTUM              ==
!     ==========================================================================
      LMAX=MAXVAL(LOFI(1:nb))
      ALLOCATE(VPAULI(NR,LMAX+1))
      ALLOCATE(UDOT(NR,LMAX+1))
      ALLOCATE(UN(NR,LMAX+1))
      DO IL=1,LMAX+1
        L=IL-1
        G(:)=0.D0
        DO IB=1,NB
          IF(LOFI(IB).NE.L) CYCLE
          CALL BOUNDSTATE(GID,NR,LOFI(IB),SOFI(IB),RAD,DREL,G,0 &
     &                   ,POTIN,EOFI(IB),UN(:,IL))
          J=IB
          G(:)=UN(:,IL)
        ENDDO
        IB=J
        CALL SCHROEDINGER$SPHERICAL(GID,NR,POTIN,DREL,SOFI(IB),G,LOFI(IB) &
     &                             ,EOFI(IB),1,UDOT(:,IL))
!       ==  avoid divide by zero if udot=0 near the nucleus
        i=1
        do ir=1,nr
          if(udot(ir,il).ne.0.d0) then
            i=ir
            exit
          end if
        enddo
!       == now construct pauli potential
        VPAULI(i:,IL)=-UN(i:,IL)/UDOT(i:,IL)
        VPAULI(1:i,IL)=VPAULI(i,IL)
        vpauli(:,il)=vpauli(:,il)+potin(:)
      ENDDO
      vpauli(irbox:,:)=0.d0
      THIS%VPAULI=VPAULI
      THIS%UDOT=UDOT
      THIS%UN=UN
      DEALLOCATE(VPAULI)
      DEALLOCATE(UN)
      DEALLOCATE(UDOT)
!
!     ==========================================================================
!     == PUT EIGENVALUES AND OCCUPATIONS BACK                                 ==
!     ==========================================================================
      THIS%EOFI(:)=EOFI(:)
      THIS%FOFI(:)=FOFI(:)
      THIS%RAD=RAD
      THIS%Q=Q
      THIS%DEDRAD=DEDRAD
      THIS%DEDQ=DEDQ
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ONEATOM_OPTF(N,Q,KBT,L,E,F,MU)
!     **************************************************************************
!     **  determine occupations for a given charge and temperature            **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N       ! #(angular momentum shells)
      REAL(8)   ,INTENT(IN) :: Q       ! total number of electrons
      REAL(8)   ,INTENT(IN) :: KBT     ! temperature as K_B*T
      INTEGER(4),INTENT(IN) :: L(N)    ! main angular momentum quantum number
      REAL(8)   ,INTENT(IN) :: E(N)    ! energy eigenvalue
      REAL(8)   ,INTENT(OUT):: F(N)    ! number of electrons per shell
      REAL(8)   ,INTENT(OUT):: MU      ! chemical potential
      REAL(8)               :: MUMIN,MUMAX
      REAL(8)               :: QMAX
      REAL(8)               :: DQ
      REAL(8)               :: X
      INTEGER(4),PARAMETER  :: NITER=1000
      REAL(8)   ,PARAMETER  :: QTOL=1.D-10
      INTEGER(4)            :: ITER,I
!     **************************************************************************
      MUMIN=MINVAL(E)-10.D0*KBT
      MUMAX=MAXVAL(E)+10.D0*KBT
      QMAX=2.D0*REAL(SUM(2*L(:)+1),KIND=8)
      IF(Q.LT.0.D0.OR.Q.GT.QMAX) THEN
        CALL ERROR$MSG('CHARGE OUT OF RANGE')
        CALL ERROR$R8VAL('Q',Q)
        CALL ERROR$R8VAL('QMIN',0.D0)
        CALL ERROR$R8VAL('QMAX',QMAX)
        CALL ERROR$STOP('ONEATOM_OPTF')
      END IF
      DO ITER=1,NITER
        MU=0.5D0*(MUMIN+MUMAX)
        DO I=1,N
          X=(E(I)-MU)/KBT
          IF(X.GT.50.D0) THEN
            F(I)=0.D0
          ELSE IF(X.LT.-50.D0) THEN
            F(I)=REAL(2*(2*L(I)+1),KIND=8)
          ELSE
            F(I)=REAL(2*(2*L(I)+1),KIND=8)/(1.D0+EXP(X))
          END IF
        ENDDO
        DQ=SUM(F)-Q
        IF(ABS(DQ).LT.QTOL) EXIT
        IF(DQ.GT.0.D0) THEN
          MUMAX=MU
        ELSE
          MUMIN=MU
        END IF
      ENDDO
      IF(ABS(DQ).GT.QTOL) THEN
        CALL ERROR$MSG('LOOP NOT CONVERGED')
        CALL ERROR$MSG('ONEATOM_OPTF')
      END IF
      RETURN
      END
!
!     .....................................................VOUT.........
      SUBROUTINE MYVOFRHOWITHRAD(GID,NR,RAD,AEZ,RHO,POT,EH,EXC)
!     ******************************************************************
!     **                                                              **
!     **  ELECTROSTATIC AND EXCHANGE-CORRELATION POTENTIAL            **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8),   INTENT(IN) :: AEZ
      REAL(8),   INTENT(IN) :: RAD
      REAL(8),   INTENT(IN) :: RHO(NR)
      REAL(8),   INTENT(OUT):: POT(NR)
      REAL(8),   INTENT(OUT):: EH
      REAL(8),   INTENT(OUT):: EXC
      REAL(8)               :: ALPHA
      REAL(8),   PARAMETER  :: RHOMIN=1.D-2
      REAL(8)               :: PI
      REAL(8)               :: FOURPI
      REAL(8)               :: Y0
      REAL(8)               :: RHO1(NR)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: POTH(NR)
      REAL(8)               :: POTXC(NR)  
      REAL(8)               :: RHOPLUS(NR)
      REAL(8)               :: DRHOPLUSDRHO(NR)
      REAL(8)               :: EDEN(NR)
      REAL(8)               :: GRHO(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: VGXC,VXC,EXC1,RH,GRHO2
      REAL(8)               :: DUMMY1,DUMMY2,DUMMY3
      REAL(8)               :: EVAL,MUVAL
      REAL(8)               :: SVAR
      REAL(8)               :: F,X,DF,YP,YM,NBYF,DFDN
      INTEGER(4)            :: IR,IRBOX
      REAL(8)               :: Q
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      FOURPI=4.D0*PI
      Y0=1.D0/SQRT(FOURPI)
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IRBOX=IR
        IF(R(IR).GT.RAD) EXIT
      ENDDO
!
!     ==================================================================
!     ==  TOTAL CHARGE                                                ==
!     ==================================================================
      AUX(:)=4.D0*PI*RHO(:)*R(:)**2*Y0
      CALL RADIAL$INTEGRATE(GID,NR,AUX,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,Q)
      Q=Q-AEZ
!
!     ==================================================================
!     ==  TOTAL POTENTIAL                                             ==
!     ==================================================================
      RHO1=RHO
      IR=MIN(IRBOX+3,NR)
      RHO1(IR:)=0.D0
      CALL RADIAL$NUCPOT(GID,NR,AEZ,POTH)
      EDEN(:)=0.5D0*RHO1(:)*POTH(:)
      CALL RADIAL$POISSON(GID,NR,0,RHO1,AUX)
      POTH(:)=POTH(:)+AUX(:)
      CALL RADIAL$VALUE(GID,NR,POTH,RAD,SVAR)
      SVAR=Q/RAD/Y0-SVAR
      POTH(1:IRBOX-1)=POTH(1:IRBOX-1)+SVAR
      POTH(IRBOX:NR)=Q/R(IRBOX:NR)/Y0
      EDEN(:)=EDEN(:)+0.5D0*RHO1(:)*POTH(:)
!
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,EDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EH)
!
!     ==================================================================
!     ==  EXCHANGE CORRELATION
!     ==================================================================
      RHO1(:IRBOX)=RHO(:)
      RHO1(IRBOX:)=0.D0     
      ALPHA=LOG(2.D0)/RHOMIN
      DO IR=1,NR
        X=RHO1(IR)*Y0
        IF(R(IR).GT.RAD) X=0.D0
        IF(ALPHA*X.GT.30.D0) THEN
          F=X
          DF=1.D0
        ELSE
          YP=EXP(ALPHA*X)
          YM=1.D0/YP
          F=LOG(YP+YM)/ALPHA
          DF=(YP-YM)/(YP+YM)
        END IF
        RHOPLUS(IR)=F/Y0
        DRHOPLUSDRHO(IR)=DF/Y0
      ENDDO
!
      CALL RADIAL$DERIVE(GID,NR,RHOPLUS(:),GRHO)
      DO IR=1,NR
        RH=RHOPLUS(IR)*Y0
        GRHO2=(Y0*GRHO(IR))**2
        CALL DFT(RH,0.D0,GRHO2,0.D0,0.D0,EXC1,VXC,DUMMY1,VGXC,DUMMY2,DUMMY3)
        EDEN(IR)=4.D0*PI*EXC1   ! ANGULAR INTEGRATION ALREADY INCLUDED
        POTXC(IR)=VXC/Y0
        GRHO(IR)=VGXC*2.D0*GRHO(IR)
      ENDDO
      CALL RADIAL$DERIVE(GID,NR,GRHO(:),AUX)
      POTXC(:)=POTXC(:)-AUX
      IF(R(1).GT.1.D+10) THEN
         POTXC(:)=POTXC(:)-2.D0/R(:)*GRHO(:)
      ELSE
        POTXC(2:)=POTXC(2:)-2.D0/R(2:)*GRHO(2:)
        POTXC(1)=POTXC(1)-2.D0/R(2)*GRHO(2)
        POTXC(1)=POTXC(2)
      END IF
!
      DO IR=1,NR
        EVAL=EDEN(IR)/(4.D0*PI)
        MUVAL=POTXC(IR)*Y0
!       ==  N/F(N) =============================
        NBYF=RHO1(IR)/RHOPLUS(IR)
        DFDN=DRHOPLUSDRHO(IR)*Y0
        SVAR=(1.D0-NBYF*DFDN)/RHOPLUS(IR)
        EDEN(IR)=NBYF*EVAL
        POTXC(IR)=NBYF*MUVAL+SVAR*EVAL
        EDEN(IR)=EDEN(IR)*4.D0*PI
        POTXC(IR)=POTXC(IR)/Y0
      ENDDO
      EDEN(:)=EDEN(:)*R(:)**2
      CALL RADIAL$INTEGRATE(GID,NR,EDEN,AUX1)
      CALL RADIAL$VALUE(GID,NR,AUX1,RAD,EXC)
      POT=POTH+POTXC
      RETURN
      END
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE STRUCTURE_MODULE
LOGICAL(4),SAVE      :: TINI=.FALSE.
INTEGER(4),PARAMETER :: NAT=2
REAL(8)              :: RBAS(3,3)
REAL(8)              :: RPOS(3,NAT)
REAL(8)              :: BAREVOL
REAL(8)              :: MADMAT(NAT,NAT)
END MODULE STRUCTURE_MODULE
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ESTRUCTURE_INITIALIZE()
!     **************************************************************************
!     **************************************************************************
      USE STRUCTURE_MODULE
      IMPLICIT NONE
      REAL(8)         ::GBAS(3,3)
      INTEGER(4)      :: IAT
      REAL(8)         :: Q(NAT)
      REAL(8)         :: ETOT
      REAL(8)         :: POT(NAT)
      REAL(8)         :: FORCE(3,NAT)
!     **************************************************************************
      IF(TINI) RETURN
      TINI=.TRUE.
!     == FCC ===================================================================
      RBAS(:,1)=(/0.0D0,0.5D0,0.5D0/)
      RBAS(:,2)=(/0.5D0,0.0D0,0.5D0/)
      RBAS(:,3)=(/0.5D0,0.5D0,0.0D0/)
      RPOS(:,1)=(/0.0D0,0.0D0,0.0D0/)
      RPOS(:,2)=(/0.5D0,0.5D0,0.5D0/)
!
!     == UNIT CELL VOLUME ======================================================
      CALL GBASS(RBAS,GBAS,BAREVOL)
!     == DETERMINE MADELUNG MATRIX =============================================
      DO IAT=1,NAT
        Q(:)=0.D0
        Q(IAT)=1.D0
        CALL MADELUNG(NAT,RBAS,RPOS,Q,ETOT,POT,FORCE)
        MADMAT(:,IAT)=POT(:)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE ESTRUCTURE(NAT_,Q,VOL,ETOT,POT,PRESSURE)
!     **************************************************************************
!     **************************************************************************
      USE STRUCTURE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT_
      REAL(8)   ,INTENT(IN) :: Q(NAT_)
      REAL(8)   ,INTENT(IN) :: VOL
      REAL(8)   ,INTENT(OUT):: ETOT
      REAL(8)   ,INTENT(OUT):: POT(NAT_)
      REAL(8)   ,INTENT(OUT):: PRESSURE
      REAL(8)               :: FORCE(3,NAT)
      REAL(8)               :: SCALE
!     **************************************************************************
      CALL ESTRUCTURE_INITIALIZE()
      IF(NAT_.NE.NAT) THEN
        CALL ERROR$MSG('NAT INCONSISTENT')
        CALL ERROR$STOP('ESTRUCTURE')
      END IF
      SCALE=(VOL/BAREVOL)**(1.D0/3.D0)
      CALL MADELUNG(NAT,RBAS*SCALE,RPOS*SCALE,Q,ETOT,POT,FORCE)
! POT(:)=MATMUL(MADMAT,Q)/SCALE
! ETOT(:)=0.5D0*DOT_PRODUCT(POT,Q)
      PRESSURE=ETOT/(3.D0*VOL)
      RETURN
      END

