
!*******************************************************************************
!**  TODO:                                                                    **
!**     - DYNOCC$GETR8: WHEN TOTCHA IS REQUESTED THE DIFFERENCE BETWEEN THE   **
!**           NUMBER OF ELECTRON (INTERNALLY ALSO CALLED TOTCHA) AND SUMOFZ   **
!**           IS RETURNED, THIS IS THE TOTAL CHARGE OF THE SYSTEM, NOT THE    **
!**           NUMBER OF ELECTRONS                                             **
!*******************************************************************************
!
!........1.........2.........3.........4.........5.........6.........7.........8
MODULE DYNOCC_MODULE
!*******************************************************************************
!**                                                                           **
!**  NAME: DYNOCC                                                             **
!**                                                                           **
!**  PURPOSE: ORGANIZES THE OCCUPATIONS OF THE ONE-PARTICLE STATES            **
!**                                                                           **
!**  FUNCTIONS:                                                               **
!**    DYNOCC$CREATE                                                          **
!**    DYNOCC$SET                                                             **
!**    DYNOCC$GET                                                             **
!**    DYNOCC$STOP                                                            **
!**    DYNOCC$INIOCC                                                          **
!**    DYNOCC$MODOCC                                                          **
!**    DYNOCC$PROPAGATE                                                       **
!**    DYNOCC$SWITCH                                                          **
!**    DYNOCC$REPORT                                                          **
!**                                                                           **
!**  REMARKS:                                                                 **
!**    THE OCCUPATIONS F ARE DIRECTLY RELATED TO THE DYNAMICAL                **
!**    VARIABLES X BY F=3*X**2-2*X**3                                         **
!**                                                                           **
!**  LOGICAL INPUT HIRARCHY:                                                  **
!**    NB,NKPT,NSPIN                                                          **
!**    SUMOFZ                                                                 **
!**    TDYN                                                                   **
!**    IF(TDYN) THEN                                                          **
!**      TEMP,MX,DELAT,ANNEX,TSTOP,FREEZE                                     **
!**      TFIXTOT                                                              **
!**      TFIXSPIN                                                             **
!**      IF(TFIXTOT) THEN                                                     **
!**        TOTCHA   (ON INPUT AS EXCESS CHARGE)                               **
!**      ELSE                                                                 **
!**        TOTPOT                                                             **
!**      END IF                                                               **
!**      IF(TFIXSPIN) THEN                                                    **
!**        SPINCHA                                                            **
!**      ELSE                                                                 **
!**        SPINPOT                                                            **
!**      END IF                                                               **
!**    ELSE                                                                   **
!**      OCC                                                                  **
!**    END IF                                                                 **
!**    STARTTYPE                                                              **
!**                                                                           **
!**    FOUR DIFFERENT STARTTYPE ARE POSSIBLE:                                 **
!**      'X' READS THE OCCUPATION VARIABLES FROM RESTART FILE                 **
!**      'E' READS THE ENERGIES FROM RESTART FILE                             **
!**      'N' FILLS STATES ACCORDING TO BAND NUMBER AND DYNOCC$MODOCC          **
!**                                                                           **
!**    THE DEFAULT OCCUPATIONS CORRESPOND TO STARTTYPE='N'.                   **
!**    IT FILLS THE BANDS CONSISTENT WITH THE TOTAL CHARGE AND SPIN           **
!**    ASSUMING FLAT BANDS. THEN THEY ARE MODIFIED BY DYNOCC$MODOCC           **
!**    (DYNOCC$MODOCC ALSO ADJUSTS THE TOTAL CHARGE AND SPIN!)                **
!**    FOR STARTTYPE='E' OR 'X' THE RESTART FILE IS READ AND THE              **
!**    OCCUPATIONS ARE OVERWRITTEN. IF STARTTYPE='X' THE OCCUPATIONS          **
!**    ARE READ FROM FILE. FOR STARTTYPE='E' THE ENERGIES ARE READ            **
!**    FROM FILE AND THE OCCUPATIONS ARE CREATED NEW USING THE                **
!**    TEMPERATURE, THE TOTAL CHARGE OR SPIN (OR SPINPOT).                    **
!**                                                                           **
!******************************************* PETER E. BLOECHL, 1996 ************
CHARACTER(1):: STARTTYPE       ! CAN BE 'N','E','X'
INTEGER(4)  :: NB=0            ! #(BANDS)
INTEGER(4)  :: NKPT=0          ! #(K-POINTS)
INTEGER(4)  :: NKDIV(3)        ! #NKDIV
INTEGER(4)  :: ISHIFT(3)       ! #ISHIFT
INTEGER(4)  :: NSPIN=0         ! #(SPINS)
LOGICAL(4)  :: TDYN=.FALSE.    ! DYNAMICAL/STATIC OCCUPATION
LOGICAL(4)  :: RESET=.TRUE.    ! SETS OUTPUT MODE OF DYNOCC$REPORT
LOGICAL(4)  :: TSTOP=.FALSE.   ! SET VELOCITY TO ZERO
LOGICAL(4)  :: TFIXSPIN=.FALSE.! FIXED SPIN/ MAGNETIC FIELD 
LOGICAL(4)  :: TFIXTOT=.FALSE. ! FIXED CHARGE/CHEMICAL POTENTIAL
LOGICAL(4)  :: TRESTARTFILEPRESENT=.FALSE. ! USE INFORMATION FROM RESTART FILE 
LOGICAL(4)  :: TADIABATIC=.FALSE.  ! QUASI-ADIABTIC CALCULATION
REAL(8)     :: RETARD=0.D0     ! TIME SCALE FOR RETARDATION OF ENERGY LEVELS.
                               ! USED ONLY WITH TADIABATIC=T
LOGICAL(4)  :: TPROPAGATED=.FALSE. ! TOTAL AND KINETIC ENERGY AVAILABLE
REAL(8)     :: FMAX=1.D0       ! MAX OCCUPATION OF A SINGLE STATE
REAL(8)     :: SUMOFZ=-1.D0    ! SUM OF NUCLEAR CHARGES
REAL(8)     :: CHARGE=0.D0     ! EXCESS NUMBER OF ELECTRONS
REAL(8)     :: TOTCHA=0.D0     ! NUMBER OF ELECTRONS (=SUMOFZ+CHARGE)
REAL(8)     :: SPINCHA=0.D0    ! TOTAL SPIN [ELECTRON SPINS]
REAL(8)     :: MX=800.D0       ! MASS OF OCCUPATION DYNAMICS
REAL(8)     :: TOTPOT=0.D0     ! FERMI LEVEL
REAL(8)     :: SPINPOT=0.D0    ! MAGNETIC FIELD
REAL(8)     :: TEMP=0.D0       ! TEMPERATURE
REAL(8)     :: ANNEX=0.D0      ! FRICTION 
REAL(8)     :: DELTAT=0.D0     ! TIME STEP
CHARACTER(32) :: BZITYPE='SAMP'  ! TYPE OF BRILLOUIN-ZONE INTEGRATION:
                           ! CAN BE 'SAMP'   FOR MONCKHORST-PACK SAMPLING
                           !     OR 'TETRA+' FOR THE IMPROVED TETRAHEDRON METHOD
REAL(8)   ,ALLOCATABLE :: XM(:,:,:)      ! DYNAMICAL VARIABLES FOR OCCUPATIONS
REAL(8)   ,ALLOCATABLE :: X0(:,:,:)
REAL(8)   ,ALLOCATABLE :: XP(:,:,:)
LOGICAL(4),ALLOCATABLE :: TFROZEN(:,:,:)  ! STATES ARE FROZEN/DYNAMIC
LOGICAL(4)             :: TEPSILON=.FALSE.
REAL(8)   ,ALLOCATABLE :: EPSILON(:,:,:)  ! DE/DF=<PSI|H|PSI>
!==== INTEGRATION WEIGHTS (INCLUDING SPIN-DEGENERACY AND K-POINT WEIGHT)
REAL(8)   ,ALLOCATABLE :: WGHT(:,:,:)    
!==== ONE-PARTICLE ENERGIES USED TO DETERMINE INTEGRATION WEIGHTS
REAL(8)   ,ALLOCATABLE :: EPS0(:,:,:)
REAL(8)   ,ALLOCATABLE :: EPSP(:,:,:)
REAL(8)   ,ALLOCATABLE :: EPSM(:,:,:)
! REMARK: MPSIDOT2 IS NOT PROPERLY TREATED, BECAUSE IT IS CALCULATED IN  
!         TWO STEPS
REAL(8)   ,ALLOCATABLE :: MPSIDOT2(:,:,:) ! M_PSI<PSIDOT||PSIDOT> 
!====== K-POINT RELATED
REAL(8)   ,ALLOCATABLE :: XK(:,:) !(3,NKPT)
REAL(8)   ,ALLOCATABLE :: WKPT(:) !(NKPT)
REAL(8)                :: EKIN
REAL(8)                :: EPOT
REAL(8)                :: DEVENERGY  ! ESTIMATED DEVIATION FROM THE GROUND STATE ENERGY
REAL(8)                :: MAXDEVOCC  ! MAX DEVIATION OF STATE OCCUPATIONS
END MODULE DYNOCC_MODULE
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC_FOFX(X,F,DF)
!      **                                                                     **
!      **  EVALUATES THE OCCUPATIONS FROM THE DYNAMICAL VARIABLES             **
!      **                                                                     **
!      *************************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: X
       REAL(8),INTENT(OUT):: F
       REAL(8),INTENT(OUT):: DF   !DF/DX
       REAL(8),PARAMETER  :: PI=4.D0*ATAN(1.D0)
!      *************************************************************************
       F =0.5D0*(1.D0-COS(X*PI))
       DF=0.5D0*PI*SIN(X*PI)
       RETURN
       END 
!
!      .................................................................
       SUBROUTINE DYNOCC_XOFF(F,X,DX)
!      **                                                             **
!      **  EVALUATES THE OCCUPATIONS FROM THE DYNAMICAL VARIABLES     **
!      **                                                             **
!      **  REMARK: THERE ARE SEVERAL SOLUTIONS. THIS ROUTINE          **
!      **          SELECTS THE ONE IN THE INTERVAL [0,1]              **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: F
       REAL(8)   ,INTENT(INOUT):: X
       REAL(8)   ,INTENT(OUT):: DX    
       REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
       REAL(8)               :: DFDX
!      *****************************************************************
!      *****************************************************************
       X=ACOS(1.D0-2.D0*F)/PI
!      == 0<X<1  =========================================================
       DFDX=0.5D0*PI*SIN(PI*X)
       DX=1.D0/DFDX
       RETURN
       END 
!
!      .................................................................
       SUBROUTINE DYNOCC_SOFX(X,S,DS)
!      **                                                             **
!      **  EVALUATES THE ENTROPY FROM THE DYNAMICAL VARIABLES         **
!      **                                                             **
!      **  REMARK: ASSUMTION IS THAT F(X) APPROACHES  THE VALUES      **
!      **          0 AND ONE QUADRATICALLY             ]              **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN)  :: X
       REAL(8),INTENT(OUT) :: S
       REAL(8),INTENT(OUT) :: DS
       REAL(8)             :: X1
       REAL(8)             :: F
       REAL(8)             :: DF
       REAL(8),PARAMETER   :: FMIN=1.D-8
       LOGICAL(4),SAVE     :: TINI_SOFX=.FALSE.
       REAL(8)   ,SAVE     :: AFAC
       REAL(8)   ,SAVE     :: XMIN
       REAL(8)             :: SVAR
!      *****************************************************************
!      =================================================================
!      == EXTRAPOLATION FORMULA FOR ARGUMENTS CLOSE TO F=0 OR F=1     ==
!      == REASON: THE NUMERICAL IMPLEMENTATION OF THE LOGARITHM       ==
!      ==        LIMITED AND DOES NOT EXTEND SUFFICIENTLY CLOSE TO ZERO=
!      == THE INTERPOLATION ASSUMES THAT FOR SMALL ARGUMENTS          ==
!      == F(X) IS PROPORTIONAL TO X**2                                ==
!      ==         S(X)=AFAC*X**2     FOR X<<1                         ==
!      ==         S(X)=AFAC*(1-X)**2 FOR (1-X)<<1                     ==
!      =================================================================
       IF(.NOT.TINI_SOFX) THEN
!        == EXTRAPOLATE S TO ZERO WITH S(X)=A*X**B
         F=FMIN
         CALL DYNOCC_XOFF(F,XMIN,SVAR)
         S=-(F*LOG(F)+(1.D0-F)*LOG(1.D0-F))
         AFAC=S/XMIN**2
         TINI_SOFX=.TRUE.
       END IF
!
!      =================================================================
!      == FOLD X BACK INTO THE INTERVAL [0,1]                         ==
!      =================================================================
       X1=X
       IF(X1.LT.-0.5D0.OR.X1.GT.1.5D0) THEN
         IF(X1.LT.0) THEN
           X1=X1+2*NINT(0.5D0*X1)+2.D0
         END IF
         X1=X1-2.D0*INT(0.5D0*X1)
         IF(X1.GT.1.D0) X1=2.D0-X1
       END IF
!
!      =================================================================
!      == NOW DETERMINE ENTROPY                                       ==
!      =================================================================
       IF(ABS(X1).LT.XMIN) THEN
         S=AFAC*X1**2
         DS=2.D0*AFAC*X1
       ELSE IF(ABS(X1-1.D0).LT.XMIN) THEN
         S=AFAC*(1.D0-X1)**2
         DS=-2.D0*AFAC*(1.D0-X1)
       ELSE
         CALL DYNOCC_FOFX(X1,F,DF)
         S=-(F*LOG(F)+(1.D0-F)*LOG(1.D0-F))
         DS=-(LOG(F)-LOG(1.D0-F))*DF
       END IF
       RETURN
       END
!!$!
!!$!      .................................................................
!!$       SUBROUTINE DYNOCC_SOFX(X,S,DS)
!!$!      **                                                             **
!!$!      **  EVALUATES THE ENTROPY FROM THE DYNAMICAL VARIABLES         **
!!$!      **                                                             **
!!$!      **  REMARK: ASSUMTION IS THAT F(X) APPROACHES  THE VALUES      **
!!$!      **          0 AND ONE QUADRATICALLY             ]              **
!!$!      **                                                             **
!!$!      *****************************************************************
!!$       IMPLICIT NONE
!!$       REAL(8),INTENT(IN)  :: X
!!$       REAL(8),INTENT(OUT) :: S
!!$       REAL(8),INTENT(OUT) :: DS
!!$       REAL(8)             :: F
!!$       REAL(8)             :: DF
!!$!      *****************************************************************
!!$       CALL DYNOCC_FOFX(X,F,DF)
!!$       IF(F.EQ.0.D0.OR.F.EQ.1.D0) THEN
!!$         S=0.D0
!!$         DS=0.D0
!!$         RETURN
!!$       END IF
!!$       S=-(F*LOG(F)+(1.D0-F)*LOG(1.D0-F))
!!$       DS=-(LOG(F)-LOG(1.D0-F))*DF
!!$       RETURN
!!$       END
!      .................................................................
       SUBROUTINE DYNOCC$TEST
       INTEGER(4),PARAMETER :: NB=5
       INTEGER(4),PARAMETER :: NKPT=1
       INTEGER(4),PARAMETER :: NSPIN=1
       REAL(8)              :: WKPT(NKPT)
       REAL(8)              :: EPSILON(NB,NKPT,NSPIN)
       REAL(8)              :: MPSIDOT2(NB,NKPT,NSPIN)
       REAL(8)              :: WGHT(NB,NKPT,NSPIN)
       REAL(8)   ,PARAMETER :: SUMOFZ=2.D0
       REAL(8)              :: MASS=1.D+3
       REAL(8)              :: DELTAT=10.D0
       REAL(8)              :: FRIC=0.D-3
       REAL(8)              :: TEMP=1.D-3
       INTEGER(4),PARAMETER :: NITER=10000
       INTEGER(4)           :: IB,IKPT,ISPIN,ITER
       INTEGER(4)           :: NFILO
       REAL(8)              :: SVAR
       REAL(8)              :: EKIN
       REAL(8)              :: EPOT
       REAL(8)              :: EBAND
       REAL(8)              :: XK(3,NKPT)
!      ******************************************************************
       CALL FILEHANDLER$UNIT('PROT',NFILO)
       DO IKPT=1,NKPT
         WKPT=1.D0/REAL(NKPT,KIND=8)
         XK(:,IKPT)=0.D0
       ENDDO
       DO ISPIN=1,NSPIN
         DO IKPT=1,NKPT
           DO IB=1,NB
             CALL LIB$RANDOM(SVAR)
             EPSILON(IB,IKPT,ISPIN)=SVAR
           ENDDO
           WRITE(NFILO,FMT='("ENERGIES FOR K=",I5," AND SPIN=",I5)')IKPT,ISPIN
           WRITE(NFILO,FMT='("EIG",10F8.3)')EPSILON(:,IKPT,ISPIN)
         ENDDO
       ENDDO
       MPSIDOT2(:,:,:)=0.D0
       CALL DYNOCC$SETL4('DYN',.TRUE.)
       CALL DYNOCC$SETR8A('XK',3*NKPT,XK)
       CALL DYNOCC$SETR8A('WKPT',NKPT,WKPT)
       CALL DYNOCC$SETR8('SUMOFZ',SUMOFZ)
       CALL DYNOCC$SETR8('MASS',MASS)
       CALL DYNOCC$SETR8('TIMESTEP',DELTAT)
       CALL DYNOCC$SETR8('TEMP',1.D0)
       CALL DYNOCC$SETL4('FIXQ',.TRUE.)
       CALL DYNOCC$SETL4('FIXS',.TRUE.)
! CALL DYNOCC$SETL4('FIXQ',.FALSE.)
! CALL DYNOCC$SETL4('FIXS',.FALSE.)
! CALL DYNOCC$SETR8('EFERMI',0.33D0)
! CALL DYNOCC$SETR8('MAGNETICFIELD',0.0D0)
       CALL DYNOCC$SETCH('STARTTYPE','N')
       CALL DYNOCC$CREATE(NB,NKPT,NSPIN)
!REMARK WHY IS XK NEEDED?
       CALL DYNOCC$INIOCC
       CALL DYNOCC$SETR8('TEMP',TEMP)
       CALL DYNOCC$REPORT(NFILO)
       CALL DYNOCC$SETR8('FRICTION',FRIC)
       WRITE(NFILO,FMT='(A10,4A15)')'ITER','ECONS','EKIN','EPOT','EBAND'
       DO ITER=1,NITER
         CALL DYNOCC$SETR8A('EPSILON',NB*NKPT*NSPIN,EPSILON)
         CALL DYNOCC$SETR8A('M<PSIDOT|PSIDOT>',NB*NKPT*NSPIN,MPSIDOT2)
         CALL DYNOCC$PROPAGATE
         CALL DYNOCC$GETR8('EKIN',EKIN)
         CALL DYNOCC$GETR8('EPOT',EPOT)
         CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,WGHT)
         EBAND=0.D0
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             SVAR=0.D0
             DO IB=1,NB
!              == NOTE THAT WGHT CONTAINS ALREADY FMAX AND WKPT
               SVAR=SVAR+WGHT(IB,IKPT,ISPIN)*EPSILON(IB,IKPT,ISPIN)
             ENDDO
             EBAND=EBAND+SVAR
           ENDDO
         ENDDO
         WRITE(NFILO,FMT='(I10,4E15.5)')ITER,EKIN+EPOT+EBAND,EKIN,EPOT,EBAND
         CALL DYNOCC$SWITCH
       ENDDO
        CALL DYNOCC$SETR8A('EPSILON',NB*NKPT*NSPIN,EPSILON)
       CALL DYNOCC$REPORT(NFILO) 
       CALL ERROR$MSG('REGULAR STOP IN TEST ROUTINE')
       CALL ERROR$STOP('DYNOCC$TEST')
       STOP
       END
!
!      .................................................................
       SUBROUTINE DYNOCC$CREATE(NB_,NKPT_,NSPIN_)
!      *****************************************************************
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NB_
       INTEGER(4),INTENT(IN) :: NKPT_
       INTEGER(4),INTENT(IN) :: NSPIN_
!      *****************************************************************
       IF(.NOT.ALLOCATED(XK)) THEN
         CALL ERROR$MSG('K-POINTS DO NOT EXIST')
         CALL ERROR$STOP('DYNOCC$CREATE')
       END IF
       IF(SUMOFZ.LT.0.D0) THEN
         CALL ERROR$MSG('SUMOFZ NOT DEFINED')
         CALL ERROR$STOP('DYNOCC$CREATE')
       END IF
       IF(.NOT.ALLOCATED(WKPT)) THEN
         CALL ERROR$MSG('K-POINT WEIGHTS DO NOT EXIST')
         CALL ERROR$STOP('DYNOCC$CREATE')
       END IF
       NB=NB_
       IF(NKPT.NE.NKPT_.AND.NKPT.NE.0) THEN
         CALL ERROR$MSG('NKPT MUST NOT BE CHANGED')
         CALL ERROR$STOP('DYNOCC$CREATE')
       END IF 
       NKPT=NKPT_
       NSPIN=NSPIN_
       ALLOCATE(XM(NB,NKPT,NSPIN))
       ALLOCATE(X0(NB,NKPT,NSPIN))
       ALLOCATE(XP(NB,NKPT,NSPIN))
       ALLOCATE(EPSILON(NB,NKPT,NSPIN))
       ALLOCATE(TFROZEN(NB,NKPT,NSPIN))
       XM(:,:,:)=0.5D0
       X0(:,:,:)=0.5D0
       XP(:,:,:)=0.5D0
       TFROZEN(:,:,:)=.FALSE.
       EPSILON(:,:,:)=-1.D0
       RETURN
       END
!      .................................................................
       SUBROUTINE DYNOCC$SETL4(ID_,VAL)
!      *****************************************************************
!      **  SET LOGICAL PARAMETERS                                     **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID_
       LOGICAL(4)  ,INTENT(IN) :: VAL
!      *****************************************************************
       IF(ID_.EQ.'DYN') THEN
         TDYN=VAL
       ELSE IF(ID_.EQ.'FIXS') THEN
         TFIXSPIN=VAL
       ELSE IF(ID_.EQ.'FIXQ') THEN
         TFIXTOT=VAL
       ELSE IF(ID_.EQ.'STOP') THEN
         TSTOP=VAL
       ELSE IF(ID_.EQ.'ADIABATIC') THEN
         TADIABATIC=VAL
         IF(BZITYPE.EQ.'TETRA+'.AND.(.NOT.TADIABATIC)) THEN
           CALL ERROR$MSG('ADIABATIC=TRUE IS NOT COMPATIBLE WITH BZITYPE=TETRA+')
           CALL ERROR$CHVAL('ID_',ID_)
           CALL ERROR$STOP('DYNOCC$SETL4')
         END IF
       ELSE IF(ID_.EQ.'TETRAHEDRON') THEN
         BZITYPE='TETRA+'
         TADIABATIC=.TRUE.
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$SETL4')
       END IF
       RESET=.TRUE.
       RETURN
       END
!      .................................................................
       SUBROUTINE DYNOCC$SETI4(ID_,VAL)
!      *****************************************************************
!      **  SET LOGICAL PARAMETERS                                     **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID_
       INTEGER(4)  ,INTENT(IN) :: VAL
!      *****************************************************************
       IF(ID_.EQ.'NKPT') THEN
         IF(NKPT.NE.0) THEN
           CALL ERROR$MSG('NKPT IS ALREADY SET')
           CALL ERROR$STOP('DYNOCC$SETI4')
         END IF
         NKPT=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$SETI4')
       END IF
       RESET=.TRUE.
       RETURN
       END
!
!     ..................................................................
      SUBROUTINE DYNOCC$SETI4A(ID,LEN,VAL)
      USE DYNOCC_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: VAL(LEN)
!     ******************************************************************
      IF(ID.EQ.'NKDIV') THEN
        IF(NKPT.EQ.0) NKPT=VAL(1)*VAL(2)*VAL(3)
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('DYNOCC$SETI4')
        END IF
        NKDIV(1:3)=VAL(1:3)
      ELSE IF(ID.EQ.'ISHIFT') THEN
        IF(LEN.NE.3) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('DYNOCC$SETI4')
        END IF
        ISHIFT(1:3)=VAL(1:3)
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('DYNOCC$SETI4A')
      END IF
      RETURN
      END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC$SETR8A(ID_,LEN_,VAL)
!      *************************************************************************
!      **  SET REAL(8) ARRAY                                                  **
!      *************************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID_
       INTEGER(4)  ,INTENT(IN) :: LEN_
       REAL(8)     ,INTENT(IN) :: VAL(LEN_)
       REAL(8)                 :: SVAR
       INTEGER(4)              :: I,IB,IKPT,ISPIN,IND
!      *************************************************************************
       IF(ID_.EQ.'XK') THEN
         IF(NKPT.EQ.0) THEN
           NKPT=LEN_/3
         END IF
         IF(LEN_.NE.3*NKPT) THEN
           CALL ERROR$MSG('SIZE INCONSISTENT')
           CALL ERROR$CHVAL('ID',ID_)
           CALL ERROR$STOP('DYNOCC$SETR8A')
         END IF
         IF(.NOT.ALLOCATED(XK))ALLOCATE(XK(3,NKPT))
!== RESHAPE HAS BEEN DISABLED: GFORTRAN PRODUCES UNPREDICTABLE RESULTS
!XK(:,:)=RESHAPE(VAL,(/3,NKPT/))
         IND=0
         DO IKPT=1,NKPT
           DO I=1,3
             IND=IND+1
             XK(I,IKPT)=VAL(IND)
           ENDDO
         ENDDO
!
!      =========================================================================
       ELSE IF(ID_.EQ.'WKPT') THEN
         IF(NKPT.EQ.0) THEN
           NKPT=LEN_
         END IF
         IF(LEN_.NE.NKPT) THEN
           CALL ERROR$MSG('SIZE INCONSISTENT')
           CALL ERROR$CHVAL('ID',ID_)
           CALL ERROR$STOP('DYNOCC$SETR8A')
         END IF
         IF(.NOT.ALLOCATED(WKPT))ALLOCATE(WKPT(NKPT))
         WKPT(:)=VAL(:)
         IF(ABS(SUM(WKPT)-1.D0).GT.1.D-6) THEN
           CALL ERROR$MSG('K-POINT WEIGHTS DO NOT SUM UP TO ONE!')
           CALL ERROR$CHVAL('ID',ID_)
           CALL ERROR$STOP('DYNOCC$SETR8A')
         END IF
!
!      =========================================================================
       ELSE IF(ID_.EQ.'EPSILON') THEN
         IF(LEN_.NE.NB*NKPT*NSPIN) THEN
           CALL ERROR$MSG('DIMENSIONS INCONSISTENT')
           CALL ERROR$CHVAL('ID_',ID_)
           CALL ERROR$I4VAL('LEN_',LEN_)
           CALL ERROR$I4VAL('NB*NKPT*NSPIN',NB*NKPT*NSPIN)
           CALL ERROR$STOP('DYNOCC$SETR8A')
         END IF
!== RESHAPE HAS BEEN DISABLED: GFORTRAN PRODUCES UNPREDICTABLE RESULTS
!PRINT*,'ALLOCATED(EPSILON) ',ALLOCATED(EPS0),ALLOCATED(EPSM),ALLOCATED(EPSP),ALLOCATED(EPSILON)
!EPSILON=RESHAPE(VAL,(/NB,NKPT,NSPIN/))
!PRINT*,'DYNOCC_SETR8A',TRIM(ID_),NB,NB*NKPT*NSPIN,LEN_
!PRINT*,'DYNOCC_SETR8A',MINVAL(VAL),MAXVAL(VAL),MINVAL(EPSILON),MAXVAL(EPSILON)
         IND=0
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               IND=IND+1
               EPSILON(IB,IKPT,ISPIN)=VAL(IND)
             ENDDO
           ENDDO
         ENDDO
         TEPSILON=.TRUE.
         IF(.NOT.ALLOCATED(EPS0)) THEN
           ALLOCATE(EPSM(NB,NKPT,NSPIN))
           ALLOCATE(EPS0(NB,NKPT,NSPIN))
           ALLOCATE(EPSP(NB,NKPT,NSPIN))
           EPSM=EPSILON
           EPS0=EPSILON
           EPSP=EPSILON
         END IF
! 
!      =========================================================================
       ELSE IF(ID_.EQ.'M<PSIDOT|PSIDOT>') THEN
         IF(LEN_.NE.NB*NKPT*NSPIN) THEN
           CALL ERROR$MSG('DIMENSIONS INCONSISTENT')
           CALL ERROR$CHVAL('ID_',ID_)
           CALL ERROR$I4VAL('LEN_',LEN_)
           CALL ERROR$I4VAL('NB*NKPT*NSPIN',NB*NKPT*NSPIN)
           CALL ERROR$STOP('DYNOCC$SETR8A')
         END IF
         IF(.NOT.ALLOCATED(MPSIDOT2)) THEN
           ALLOCATE(MPSIDOT2(NB,NKPT,NSPIN))
           MPSIDOT2(:,:,:)=0.D0
           SVAR=1.0D0
         ELSE
           MPSIDOT2(:,:,:)=0.5D0*MPSIDOT2(:,:,:)
           SVAR=0.5D0
         END IF
!== RESHAPE HAS BEEN DISABLED: GFORTRAN PRODUCES UNPREDICTABLE RESULTS
!MPSIDOT2=MPSIDOT2+SVAR*RESHAPE(VAL,(/NB,NKPT,NSPIN/))
!PRINT*,'DYNOCC_SETR8A',TRIM(ID_),NB,NB*NKPT*NSPIN,LEN_
!PRINT*,'DYNOCC_SETR8A',MINVAL(VAL),MAXVAL(VAL),MINVAL(MPSIDOT2),MAXVAL(MPSIDOT2)
         IND=0
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               IND=IND+1
               MPSIDOT2(IB,IKPT,ISPIN)=MPSIDOT2(IB,IKPT,ISPIN)+SVAR*VAL(IND)
             ENDDO
           ENDDO
         ENDDO
! 
!      =========================================================================
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$SETR8A')
       END IF
       RETURN
       END
!
!      .................................................................
       SUBROUTINE DYNOCC$SETCH(ID_,VAL)
!      *****************************************************************
!      **  SET CHARACTER                                              **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID_
       CHARACTER(*),INTENT(IN) :: VAL
!      *****************************************************************
       IF(ID_.EQ.'STARTTYPE') THEN
         STARTTYPE=VAL
         IF(STARTTYPE.NE.'N'.AND.STARTTYPE.NE.'E'.AND.STARTTYPE.NE.'X') THEN
           CALL ERROR$MSG('ILLEGAL VALUE OF STARTTYPE')
           CALL ERROR$MSG('ALLOWED VALUES ARE "N", "E", AND "X"')
           CALL ERROR$CHVAL('ID_',ID_)
           CALL ERROR$STOP('DYNOCC$SETCH')
         END IF
! 
!      ===================================================================
       ELSE IF(ID_.EQ.'BZITYPE') THEN
         IF(VAL.NE.'SAMP'.AND.VAL.NE.'TETRA+') THEN
           CALL ERROR$MSG('ILLEGAL VALUE OF BZITYPE')
           CALL ERROR$MSG('ALLOWED VALUES ARE "SAMP", "TETRA+"')
           CALL ERROR$CHVAL('ID_',ID_)
           CALL ERROR$STOP('DYNOCC$SETCH')
         END IF
         BZITYPE=VAL
! 
!      ===================================================================
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$SETCH')
       END IF
       RETURN
       END
!
!      .................................................................
       SUBROUTINE DYNOCC$SETR8(ID_,VAL)
!      *****************************************************************
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE 
       CHARACTER(*),INTENT(IN) :: ID_
       REAL(8)     ,INTENT(IN) :: VAL
!      *****************************************************************
       IF(ID_.EQ.'SUMOFZ') THEN
         SUMOFZ=VAL
         TOTCHA=SUMOFZ+CHARGE
         RESET=.TRUE.
       ELSE IF(ID_.EQ.'TOTCHA') THEN
         IF(SUMOFZ.LT.0.D0) THEN
           CALL ERROR$MSG('SUMOFZ MUST BE SET BEFORE TOTCHA')
           CALL ERROR$CHVAL('ID',ID_)
           CALL ERROR$STOP('DYNOCC$SETR8')
         END IF
         CHARGE=VAL
         TOTCHA=SUMOFZ+CHARGE
         RESET=.TRUE.
       ELSE IF(ID_.EQ.'SPIN') THEN
         SPINCHA=VAL
         RESET=.TRUE.
         PRINT*,'DYNOCC$SETR8: SPIN HAS BEEN SET TO ',SPINCHA
       ELSE IF(ID_.EQ.'FMAX') THEN
         FMAX=VAL
         RESET=.TRUE.
       ELSE IF(ID_.EQ.'EFERMI') THEN
         TOTPOT=VAL
         RESET=.TRUE.
       ELSE IF(ID_.EQ.'MAGNETICFIELD') THEN
         SPINPOT=VAL
         RESET=.TRUE.
       ELSE IF(ID_.EQ.'MASS') THEN
         MX=VAL
         RESET=.TRUE.
       ELSE IF(ID_.EQ.'FRICTION') THEN
         ANNEX=VAL
         RESET=.TRUE.
       ELSE IF(ID_.EQ.'TIMESTEP') THEN
         DELTAT=VAL
       ELSE IF(ID_.EQ.'TEMP') THEN
         TEMP=VAL
         RESET=.TRUE.
       ELSE IF(ID_.EQ.'RETARD') THEN
         RETARD=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$SETR8')
       END IF
       RETURN
       END
!      .................................................................
       SUBROUTINE DYNOCC$GETL4(ID_,VAL)
!      *****************************************************************
!      **  SET LOGICAL PARAMETERS                                     **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID_
       LOGICAL(4)  ,INTENT(OUT):: VAL
!      *****************************************************************
       IF(ID_.EQ.'DYN') THEN
         VAL=TDYN
       ELSE IF(ID_.EQ.'PROPAGATED') THEN
         VAL=TPROPAGATED
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$GETL4')
       END IF
       RETURN
       END
!      .................................................................
       SUBROUTINE DYNOCC$GETI4(ID_,VAL)
!      *****************************************************************
!      **  SET INTEGER PARAMETERS                                     **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID_
       INTEGER(4)  ,INTENT(OUT):: VAL
!      *****************************************************************
       IF(ID_.EQ.'NKPT') THEN
         VAL=NKPT
       ELSE IF(ID_.EQ.'NSPIN') THEN
         VAL=NSPIN
       ELSE IF(ID_.EQ.'NB') THEN
         VAL=NB
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$GETI4')
       END IF
       RETURN
       END
!     .................................................................
      SUBROUTINE DYNOCC$GETI4A(ID_,LEN_,VAL)
!     *****************************************************************
!     **  SET INTEGER PARAMETERS                                     **
!     *****************************************************************
      USE DYNOCC_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID_
      INTEGER(4)  ,INTENT(IN) :: LEN_
      INTEGER(4)  ,INTENT(OUT):: VAL(LEN_)
!     *****************************************************************
      IF(ID_.EQ.' ') THEN
         VAL(:)=0
      ELSE IF(ID_.EQ.'NKDIV') THEN
        IF(LEN_.NE.3) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$I4VAL('LEN',LEN_)
          CALL ERROR$STOP('DYNOCC$GETI4')
        END IF
        VAL(1:3)=NKDIV(1:3)
      ELSE IF(ID_.EQ.'ISHIFT') THEN
        IF(LEN_.NE.3) THEN
          CALL ERROR$MSG('INCONSISTENT SIZE')
          CALL ERROR$CHVAL('ID',ID_)
          CALL ERROR$I4VAL('LEN',LEN_)
          CALL ERROR$STOP('DYNOCC$GETI4')
        END IF
        VAL(1:3)=ISHIFT(1:3)
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$GETI4A')
       END IF
       RETURN
       END
!
!      .................................................................
       SUBROUTINE DYNOCC$GETR8(ID_,VAL)
!      *****************************************************************
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID_
       REAL(8)     ,INTENT(OUT):: VAL
!      *****************************************************************
       IF(ID_.EQ.'SPIN') THEN
         VAL=SPINCHA    
       ELSE IF(ID_.EQ.'TOTCHA') THEN
         VAL=TOTCHA-SUMOFZ
       ELSE IF(ID_.EQ.'NEL') THEN
         VAL=TOTCHA
       ELSE IF(ID_.EQ.'FMAX') THEN
         VAL=FMAX
       ELSE IF(ID_.EQ.'EKIN') THEN
         CALL DYNOCC_EKIN(VAL)
       ELSE IF(ID_.EQ.'EPOT') THEN
         CALL DYNOCC_EPOT(VAL)
       ELSE IF(ID_.EQ.'TEMP') THEN
         VAL=TEMP
       ELSE IF(ID_.EQ.'EFERMI') THEN
         VAL=TOTPOT
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$GETR8')
       END IF
       RETURN
       END
!
!      .................................................................
       SUBROUTINE DYNOCC$GETR8A(ID,LEN,VAL)
!      *****************************************************************
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       INTEGER(4)  ,INTENT(IN) :: LEN
       REAL(8)     ,INTENT(OUT):: VAL(LEN)
       INTEGER(4)              :: IND,I,IB,IKPT,ISPIN
!      *****************************************************************
!
!      =================================================================
!      ==  K-POINT RELATIVE COORDINATES                               ==
!      =================================================================
       IF(ID.EQ.'XK') THEN
         IF(LEN.NE.3*NKPT) THEN
           CALL ERROR$MSG('SIZE INCONSISTENT')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$I4VAL('LEN',LEN)
           CALL ERROR$I4VAL('NKPT',NKPT)
           CALL ERROR$STOP('DYNOCC$GETR8A')
         END IF
!  VAL=RESHAPE(XK,(/3*NKPT/))
         IND=0
         DO IKPT=1,NKPT
           DO I=1,3
             IND=IND+1
             VAL(IND)=XK(I,IKPT)
           ENDDO
         ENDDO
!
!      =================================================================
!      ==  K-POINT GEOMETRIC INTEGRATION WEIGHTS                      ==
!      =================================================================
       ELSE IF(ID.EQ.'WKPT') THEN
         IF(LEN.NE.NKPT) THEN
           CALL ERROR$MSG('SIZE INCONSISTENT')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$I4VAL('LEN',LEN)
           CALL ERROR$I4VAL('NKPT',NKPT)
           CALL ERROR$STOP('DYNOCC$GETR8A')
         END IF
         VAL=WKPT
!
!      =================================================================
!      ==  OCCUPATIONS                                                ==
!      =================================================================
       ELSE IF(ID.EQ.'OCC') THEN
         IF(LEN.NE.NB*NKPT*NSPIN) THEN
           CALL ERROR$MSG('DIMENSIONS INCONSISTENT')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$I4VAL('LEN',LEN)
           CALL ERROR$I4VAL('NB*NKPT*NSPIN',NB*NKPT*NSPIN)
           CALL ERROR$STOP('DYNOCC$GETR8A')
         END IF
         CALL DYNOCC_WGHT()  ! RECALCULATES WEIGHTS WHEN NEEDED
! VAL=RESHAPE(WGHT,(/NB*NKPT*NSPIN/))
         IND=0
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               IND=IND+1
               VAL(IND)=WGHT(IB,IKPT,ISPIN)
             ENDDO
           ENDDO
         ENDDO
!
!      =========================================================================
!      ==  ONE-PARTICLE ENERGIES 'EPSILON' (FORCE ON OCCUPATION)              ==
!      =========================================================================
       ELSE IF(ID.EQ.'EPSILON') THEN
         IF(LEN.NE.NB*NKPT*NSPIN) THEN
           CALL ERROR$MSG('DIMENSIONS INCONSISTENT')
           CALL ERROR$CHVAL('ID',ID)
           CALL ERROR$I4VAL('LEN',LEN)
           CALL ERROR$I4VAL('NB*NKPT*NSPIN',NB*NKPT*NSPIN)
           CALL ERROR$STOP('DYNOCC$GETR8A')
         END IF
         IND=0
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               IND=IND+1
               VAL(IND)=EPSILON(IB,IKPT,ISPIN)
             ENDDO
           ENDDO
         ENDDO
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('DYNOCC$GETR8A')
       END IF
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DYNOCC$WRITE(NFIL,NFILO,TCHK)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE DYNOCC_MODULE
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)       ,INTENT(IN) :: NFIL       ! RESTART-FILE UNIT
      INTEGER(4)       ,INTENT(IN) :: NFILO
      LOGICAL(4)       ,INTENT(OUT):: TCHK
      TYPE(SEPARATOR_TYPE),PARAMETER :: MYSEPARATOR &
          =SEPARATOR_TYPE(6,'OCCUPATIONS','NONE','MAR2008',' ')
      TYPE(SEPARATOR_TYPE),PARAMETER :: OLD2SEPARATOR &
          =SEPARATOR_TYPE(4,'OCCUPATIONS','NONE','OCT2003',' ')
      TYPE(SEPARATOR_TYPE),PARAMETER :: OLDSEPARATOR &
          =SEPARATOR_TYPE(3,'OCCUPATIONS','NONE','AUG1996',' ')
      INTEGER(4)                   :: NTASKS,THISTASK
      INTEGER(4)       ,PARAMETER  :: FORMATTYPE=1
!     **************************************************************************
                          CALL TRACE$PUSH('DYNOCC$WRITE')
      TCHK=.FALSE.
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)

      IF(FORMATTYPE.EQ.1) THEN
        IF(THISTASK.EQ.1) THEN
          CALL RESTART$WRITESEPARATOR(MYSEPARATOR,NFIL,NFILO,TCHK)
          WRITE(NFIL)NB,NKPT,NSPIN
          WRITE(NFIL)X0(:,:,:)
          WRITE(NFIL)XM(:,:,:)
          IF(TDYN) THEN
            WRITE(NFIL)EPS0(:,:,:)
            WRITE(NFIL)EPSM(:,:,:)
            WRITE(NFIL)EPSILON(:,:,:)
          ELSE  ! EPSILON MUST BE AVAILABLE TO BE ABLE TO START UP
            WRITE(NFIL)EPSILON(:,:,:)
            WRITE(NFIL)EPSILON(:,:,:)
            WRITE(NFIL)EPSILON(:,:,:)
          END IF
        END IF
      ELSE IF(FORMATTYPE.EQ.2) THEN
        IF(THISTASK.EQ.1) THEN
          CALL RESTART$WRITESEPARATOR(OLD2SEPARATOR,NFIL,NFILO,TCHK)
          WRITE(NFIL)NB,NKPT,NSPIN
          WRITE(NFIL)X0(:,:,:)
          WRITE(NFIL)XM(:,:,:)
          WRITE(NFIL)EPSILON(:,:,:)
        END IF
      ELSE IF(FORMATTYPE.EQ.3) THEN
        IF(THISTASK.EQ.1) THEN
          CALL RESTART$WRITESEPARATOR(OLDSEPARATOR,NFIL,NFILO,TCHK)
          WRITE(NFIL)NB,NKPT,NSPIN
          WRITE(NFIL)X0(:,:,:)
          WRITE(NFIL)XM(:,:,:)
        END IF
      ELSE
        CALL ERROR$MSG('INVALID VALUE OF FORMATTYPE')
        CALL ERROR$STOP('DYNOCC$WRITE')
      END IF
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DYNOCC$READ(NFIL,NFILO,TCHK)
!     **************************************************************************
!     **                                                                      **
!     **************************************************************************
      USE DYNOCC_MODULE, ONLY : STARTTYPE &
     &                         ,BZITYPE &
     &                         ,TRESTARTFILEPRESENT &
     &                         ,FMAX &
     &                         ,TEMP & !TEMPERATURE
     &                         ,TFIXSPIN,SPINCHA,SPINPOT &
     &                         ,TFIXTOT,TOTCHA,TOTPOT &
     &                         ,NB,NKPT,NSPIN &
     &                         ,WKPT &
     &                         ,XM,X0,XP &
     &                         ,EPSILON &
     &                         ,EPSM,EPS0,EPSP
      USE RESTART_INTERFACE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)         ,INTENT(IN) :: NFIL       ! RESTART-FILE UNIT
      INTEGER(4)         ,INTENT(IN) :: NFILO
      LOGICAL(4)         ,INTENT(OUT):: TCHK
      TYPE(SEPARATOR_TYPE),PARAMETER :: MYSEPARATOR &
          =SEPARATOR_TYPE(6,'OCCUPATIONS','NONE','MAR2008',' ')
      TYPE(SEPARATOR_TYPE),PARAMETER :: OLD2SEPARATOR &
          =SEPARATOR_TYPE(4,'OCCUPATIONS','NONE','OCT2003',' ')
      TYPE(SEPARATOR_TYPE),PARAMETER :: OLDSEPARATOR &
          =SEPARATOR_TYPE(3,'OCCUPATIONS','NONE','AUG1996',' ')
      TYPE(SEPARATOR_TYPE)           :: SEPARATOR 
      REAL(8)            ,ALLOCATABLE:: TMP0(:,:,:)
      REAL(8)            ,ALLOCATABLE:: TMPM(:,:,:)
      REAL(8)            ,ALLOCATABLE:: TMPE0(:,:,:)
      REAL(8)            ,ALLOCATABLE:: TMPEM(:,:,:)
      REAL(8)            ,ALLOCATABLE:: TMPE(:,:,:)
      INTEGER(4)                     :: NB1,NKPT1,NSPIN1
      INTEGER(4)                     :: NTASKS,THISTASK
      INTEGER(4)                     :: ISPIN,IKPT,IB
      REAL(8)                        :: SVAR,DSVAR
      LOGICAL(4)                     :: TOLD
      INTEGER(4)                     :: NFIL1
      CHARACTER(8)                   :: FORMATTYPE
!     **************************************************************************
                          CALL TRACE$PUSH('DYNOCC$READ')
      IF(STARTTYPE.EQ.' ') THEN
        CALL ERROR$MSG('STARTTYPE NOT SPECIFIED')
        CALL ERROR$STOP('DYNOCC$READ')
      END IF
      IF(STARTTYPE.EQ.'N') THEN
        CALL TRACE$POP
        RETURN
      END IF
!
!     ==========================================================================
!     == CHECK IF FILE IS POSITIONED AT THE CORRECT POSITION                  ==
!     ==========================================================================
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      SEPARATOR=MYSEPARATOR
      IF(THISTASK.EQ.1)CALL RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL TRACE$POP ;RETURN
      END IF
      TRESTARTFILEPRESENT=.TRUE.
!
!     ==========================================================================
!     == IDENTIFY VERSION OF THE RESTART FILE                                 ==
!     ==========================================================================
      FORMATTYPE=SEPARATOR%VERSION
      IF(FORMATTYPE.EQ.'MAR2008') THEN
      ELSE IF(FORMATTYPE.EQ.'OCT2003') THEN
      ELSE IF(FORMATTYPE.EQ.'AUG1996') THEN
      ELSE 
        CALL ERROR$MSG('COULD NOT IDENTIFY VALUE OF FORMATTYPE')
        CALL ERROR$STOP('DYNOCC$READ')
      END IF
!
!     ==========================================================================
!     == CATCH INCONSISTENCIES                                                ==
!     ==========================================================================
!     == CHECK IF RESTART FILE FORMAT FROM AUGUST 1996 =================
      TOLD=(SEPARATOR%VERSION.EQ.'AUG1996')
      IF ((STARTTYPE.EQ.'E').AND.(FORMATTYPE.EQ.'AUG1996')) THEN
        CALL FILEHANDLER$UNIT('PROT',NFIL1)
        WRITE(NFIL1,*)
        WRITE(NFIL1,'("**********************************************")')
        WRITE(NFIL1,'("WARNING:")')
        WRITE(NFIL1,'("IT IS NOT POSSIBLE TO RESTART FROM A RSTRT ")')
        WRITE(NFIL1,'("FILE OF PREVIOUS RELEASES USING ")')
        WRITE(NFIL1,'("!MERMIN START=T. THE CODE AUTOMATICALLY")')
        WRITE(NFIL1,'("ASSUMES !MERMIN START=F")')
        WRITE(NFIL1,'("**********************************************")')
        WRITE(NFIL1,*)
        CALL FILEHANDLER$CLOSE('PROT')
      END IF
!
!     ==========================================================================
!     == READ DATA                                                            ==
!     ==========================================================================
      IF(.NOT.ALLOCATED(EPS0)) THEN
        ALLOCATE(EPSM(NB,NKPT,NSPIN))
        ALLOCATE(EPS0(NB,NKPT,NSPIN))
        ALLOCATE(EPSP(NB,NKPT,NSPIN))
      END IF
      X0(:,:,:)=0.D0
      XM(:,:,:)=0.D0
      EPS0=0.D0
      EPSM=0.D0
      EPSP=0.D0
      EPSILON=0.D0
      IF(THISTASK.EQ.1) THEN
        READ(NFIL)NB1,NKPT1,NSPIN1
        ALLOCATE(TMP0(NB1,NKPT1,NSPIN1))
        ALLOCATE(TMPM(NB1,NKPT1,NSPIN1))
        ALLOCATE(TMPE(NB1,NKPT1,NSPIN1))
        ALLOCATE(TMPE0(NB1,NKPT1,NSPIN1))
        ALLOCATE(TMPEM(NB1,NKPT1,NSPIN1))
        NB1=MIN(NB1,NB)
        NKPT1=MIN(NKPT1,NKPT)
        NSPIN1=MIN(NSPIN1,NSPIN)
        IF(FORMATTYPE.EQ.'MAR2008') THEN
          READ(NFIL)TMP0(:,:,:)
          READ(NFIL)TMPM(:,:,:)
          READ(NFIL)TMPE0(:,:,:)
          READ(NFIL)TMPEM(:,:,:)
          READ(NFIL)TMPE(:,:,:)
          X0(1:NB1,1:NKPT1,1:NSPIN1)=TMP0(1:NB1,1:NKPT1,1:NSPIN1)
          XM(1:NB1,1:NKPT1,1:NSPIN1)=TMPM(1:NB1,1:NKPT1,1:NSPIN1)
!         == ENSURE THAT ADDED BANDS ARE TREATED AS UNOCCUPIED==================
          EPS0(NB1+1:,:NKPT1,:NSPIN1)=MAXVAL(TMPE0)
          EPSM(NB1+1:,:NKPT1,:NSPIN1)=MAXVAL(TMPEM)
          EPSILON(NB1+1:,:NKPT1,:NSPIN1)=MAXVAL(TMPE)
!
          EPS0(1:NB1,1:NKPT1,1:NSPIN1)=TMPE0(1:NB1,1:NKPT1,1:NSPIN1)
          EPSM(1:NB1,1:NKPT1,1:NSPIN1)=TMPEM(1:NB1,1:NKPT1,1:NSPIN1)
          EPSP(:,:,:)=EPS0(:,:,:)
          EPSILON(1:NB1,1:NKPT1,1:NSPIN1)=TMPE(1:NB1,1:NKPT1,1:NSPIN1)
        ELSE IF(FORMATTYPE.EQ.'OCT2003') THEN
          READ(NFIL)TMP0(:,:,:)
          READ(NFIL)TMPM(:,:,:)
          READ(NFIL)TMPE(:,:,:)
          X0(1:NB1,1:NKPT1,1:NSPIN1)=TMP0(1:NB1,1:NKPT1,1:NSPIN1)
          XM(1:NB1,1:NKPT1,1:NSPIN1)=TMPM(1:NB1,1:NKPT1,1:NSPIN1)
          EPSILON(1:NB1,1:NKPT1,1:NSPIN1)=TMPE(1:NB1,1:NKPT1,1:NSPIN1)
          EPSM(:,:,:)=EPSILON(:,:,:)
          EPS0(:,:,:)=EPSILON(:,:,:)
          EPSP(:,:,:)=EPSILON(:,:,:)
        ELSE IF(FORMATTYPE.EQ.'AUG1996') THEN
          READ(NFIL)TMP0(:,:,:)
          READ(NFIL)TMPM(:,:,:)
          X0(1:NB1,1:NKPT1,1:NSPIN1)=TMP0(1:NB1,1:NKPT1,1:NSPIN1)
          XM(1:NB1,1:NKPT1,1:NSPIN1)=TMPM(1:NB1,1:NKPT1,1:NSPIN1)
          EPSM(:,:,:)=EPSILON(:,:,:)
          EPS0(:,:,:)=EPSILON(:,:,:)
          EPSP(:,:,:)=EPSILON(:,:,:)
        ELSE 
          CALL ERROR$MSG('FORMATTYPE SPECIFIED IN RESTART FILE NOT RECOGNIZED')
          CALL ERROR$CHVAL('FORMATTYPE',FORMATTYPE)
          CALL ERROR$STOP('DYNOCC$READ')
        END IF
        DEALLOCATE(TMP0)
        DEALLOCATE(TMPM)
        DEALLOCATE(TMPEM)
        DEALLOCATE(TMPE0)
        DEALLOCATE(TMPE)
!
!       ========================================================================
!       == AUGMENT MISSING DATA                                               ==
!       ========================================================================
!       - THE FOLLING IF-CONDITION IS A WORKAROUND FOR A BUG IN THE LOOP 
!       - VECTORIZER OF THE INTEL FORTAN COMPILER 16, THAT WOULD GENERATE 
!       - AN INFINITE LOOP HERE.
        IF(NKPT1+1.LE.NKPT)THEN
          DO IKPT=NKPT1+1,NKPT
            X0(:,IKPT,:)=X0(:,1,:)
            XM(:,IKPT,:)=XM(:,1,:)
            EPSILON(:,IKPT,:)=EPSILON(:,1,:)
            EPS0(:,IKPT,:)=EPS0(:,1,:)
            EPSM(:,IKPT,:)=EPSM(:,1,:)
          ENDDO
        ENDIF

        IF(NSPIN1.EQ.1.AND.NSPIN.EQ.2) THEN
          X0(:,:,2)=X0(:,:,1)
          XM(:,:,2)=XM(:,:,1)
          EPSILON(:,:,2)=EPSILON(:,:,1)
          EPS0(:,:,2)=EPS0(:,:,1)
          EPSM(:,:,2)=EPSM(:,:,1)
        END IF
      END IF
!
!     ==========================================================================
!     ==  BROADCAST                                                           ==
!     ==========================================================================
      CALL MPE$BROADCAST('MONOMER',1,X0)
      CALL MPE$BROADCAST('MONOMER',1,XM)
      CALL MPE$BROADCAST('MONOMER',1,EPSILON)
      CALL MPE$BROADCAST('MONOMER',1,EPS0)
      CALL MPE$BROADCAST('MONOMER',1,EPSM)
!
!     ==========================================================================
!     == CONVERT ENERGIES INTO OCCUPATIONS                                    ==
!     ==========================================================================
      IF(STARTTYPE.EQ.'E'.AND.BZITYPE.EQ.'SAMP' &
     &             .AND.(FORMATTYPE.NE.'AUG1996')) THEN
!       == NOT SUITABLE FOR TETRAHEDRON METHOD BECAUSE IT FAILS FOR T=0
        CALL DYNOCC_INIOCCBYENERGY(NB,NKPT,NSPIN,FMAX &
      &          ,TEMP,TFIXTOT,TOTCHA,TOTPOT,TFIXSPIN,SPINCHA,SPINPOT &
      &          ,WKPT,EPSILON,X0)
        XM(:,:,:)=X0(:,:,:)
      END IF
!
!     ==========================================================================
!     ==  RESET THERMODYNAMIC VARIABLES                                       ==
!     ==========================================================================
      IF(STARTTYPE.EQ.'X') THEN
        IF(NSPIN.EQ.1) THEN
          SPINCHA=0.D0
        ELSE
          SPINCHA=0.D0
          DO IKPT=1,NKPT
            DO IB=1,NB
              CALL DYNOCC_FOFX(X0(IB,IKPT,1),SVAR,DSVAR)
              SPINCHA=SPINCHA+FMAX*SVAR*WKPT(IKPT)
              CALL DYNOCC_FOFX(X0(IB,IKPT,2),SVAR,DSVAR)
              SPINCHA=SPINCHA-FMAX*SVAR*WKPT(IKPT)
            ENDDO
          ENDDO        
        END IF
        TOTCHA=0.D0
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB=1,NB
              CALL DYNOCC_FOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
              TOTCHA=TOTCHA+SVAR*FMAX*WKPT(IKPT)
            ENDDO
          ENDDO        
        ENDDO
      END IF
                           CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DYNOCC$REPORT(NFIL)
!     **************************************************************************
!     **                                                                      **
!     **  THE DATA OF THE DYNOCC OBJECT ARE DIVIDED IN STATIC                 **
!     **  DYNAMIC DATA.                                                       **
!     **  THE PARAMETER RESET IS CHANGED TO TRUE WHENEVER STATIC              **
!     **  DATA ARE CHANGED AND TO FALSE AFTER EACH REPORT.                    **
!     **  STATIC DATA ARE THEREFORE REPORTED, WHENEVER ONE OF THEM            **
!     **  HAS BEEN MODIFIED.                                                  **
!     **                                                                      **
!     **  DYNAMIC DATA ARE ALWAYS REPORTED                                    **
!     **                                                                      **
!     **************************************************************************
      USE DYNOCC_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NFIL
      REAL(8)               :: SVAR
      REAL(8)               :: OCC(NB,NKPT,NSPIN)
      INTEGER(4)            :: ISPIN,IKPT
      REAL(8)               :: EV
      REAL(8)               :: KELVIN
!     **************************************************************************
                         CALL TRACE$PUSH('DYNOCC$REPORT')
      CALL CONSTANTS('EV',EV)
      CALL CONSTANTS('KB',KELVIN)
!     == 1EL<->0.5HBAR => HBAR=2*EL
!
!     ==========================================================================
!     == TITLE OF THE REPORT                                                  ==
!     ==========================================================================
      CALL REPORT$TITLE(NFIL,'OCCUPATIONS')
!
!     ==========================================================================
!     == REPORT BASIC SETTINGS, IF THEY HAVE CHANGED SINCE THE LAST CALL
!     ==========================================================================
      IF(RESET) THEN
        CALL DYNOCC_REPORTSETTING(NFIL)
        RESET=.FALSE.
      END IF
      IF(.NOT.TDYN) THEN 
         CALL TRACE$POP()
         RETURN
      END IF
!
!     ==========================================================================
!     ==  PRINT ENSEMBLE INFORMATION                                          ==
!     ==========================================================================
      IF(TFIXTOT) THEN
        CALL REPORT$R8VAL(NFIL,'CHEMICAL POTENTIAL',TOTPOT/EV,'EV')
      ELSE
        CALL REPORT$R8VAL(NFIL,'CHARGE',-(TOTCHA-SUMOFZ),'E')
        CALL REPORT$R8VAL(NFIL,'CHARGE TERM (-MU*N)?' &
    &                    ,-TOTPOT*(TOTCHA-SUMOFZ)/EV,'EV')
      END IF
!
      IF(NSPIN.EQ.2) THEN
        IF(TFIXSPIN) THEN
          CALL REPORT$R8VAL(NFIL,'MAGNETIC FIELD',SPINPOT/(EV),'EV/(HBAR/2)')
        ELSE
          CALL REPORT$R8VAL(NFIL,'SPIN',SPINCHA/2.D0,'HBAR')
          CALL REPORT$R8VAL(NFIL,'MAGNETIZATION TERM (-S*B)?' &
    &                           ,-SPINPOT*SPINCHA/EV,'EV')
        END IF
      END IF

      IF(BZITYPE.EQ.'SAMP') THEN 
        CALL DYNOCC$GETR8('EPOT',SVAR)
        IF(TFIXTOT) THEN
          IF(TFIXSPIN) THEN
            CALL REPORT$R8VAL(NFIL,'ENTROPY TERM -TS',SVAR/EV,'EV')
          ELSE
            CALL REPORT$R8VAL(NFIL,'ENTROPY TERM -TS-B*S',SVAR/EV,'EV')
          END IF
        ELSE
          IF(TFIXSPIN) THEN
            CALL REPORT$R8VAL(NFIL,'ENTROPY TERM -TS-MU*N',SVAR/EV,'EV')
          ELSE
            CALL REPORT$R8VAL(NFIL,'ENTROPY TERM -TS-B*S-MU*N',SVAR/EV,'EV')
          END IF
        END IF
      END IF
!
!     ==========================================================================
!     ==  PRINT INFORMATION ABOUT CONVERGENCE                                 ==
!     ==========================================================================
      IF(TPROPAGATED) THEN
        CALL REPORT$R8VAL(NFIL,'MAX DEVIATION STATE OCC.',MAXDEVOCC,'E')
        CALL REPORT$R8VAL(NFIL,'ESTIMATED ENERGY DEVIATION',DEVENERGY,'H')
      END IF
!
!     ==========================================================================
!     == PRINT ENERGIES AND OCCUPATIONS                                       ==
!     ==========================================================================
      IF(TEPSILON) THEN
        CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
        WRITE(NFIL,*)'OCCUPATIONS AND ENERGY EXPECTATION VALUES [EV]'
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            WRITE(NFIL,FMT='("FOR K-POINT: ",I5," AND SPIN ",I1)')IKPT,ISPIN
            CALL DYNOCC_PREIG('EIG',NFIL,NB,EPSILON(:,IKPT,ISPIN),EV)
            CALL DYNOCC_PREIG('OCC',NFIL,NB,OCC(:,IKPT,ISPIN),WKPT(IKPT))
          ENDDO
        ENDDO
!
!       == CHECK IF THE NUMBER OF BAND IS SUFFICIENT ===========================
        IF(SUM(OCC(NB,:,:)).GT.1.D-6) THEN
          CALL REPORT$STRING(NFIL,'WARNING: HIGHEST BAND IS NOT EMPTY')
          CALL REPORT$STRING(NFIL,'         RESULTS ARE POBABLY INCORRECT')
          CALL REPORT$STRING(NFIL,'         INCREASE NUMBER OF BANDS')
        END IF
      END IF
                         CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC_REPORTSETTING(NFIL)
!      *************************************************************************
!      **  REPORTS THE SETTING OF THE DYNOCC OBJECT                           **
!      **  DYNAMIC DATA.                                                      **
!      **  THE PARAMETER RESET IS CHANGED TO TRUE WHENEVER STATIC             **
!      **  DATA ARE CHANGED AND TO FALSE AFTER EACH REPORT.                   **
!      **  STATIC DATA ARE THEREFORE REPORTED, WHENEVER ONE OF THEM           **
!      **  HAS BEEN MODIFIED.                                                 **
!      **                                                                     **
!      **  DYNAMIC DATA ARE ALWAYS REPORTED                                   **
!      **                                                                     **
!      *************************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NFIL
       REAL(8)               :: EV
       REAL(8)               :: KELVIN
       REAL(8)               :: OCC(NB,NKPT,NSPIN)
       INTEGER(4)            :: ISPIN,IKPT
       CHARACTER(80)         :: STRING
!      *************************************************************************
       CALL CONSTANTS('EV',EV)
       CALL CONSTANTS('KB',KELVIN)
!      == 1EL<->0.5HBAR => HBAR=2 EL
!
!      =========================================================================
!      ==                                                                     ==
!      =========================================================================
       IF(.NOT.TDYN) THEN
         CALL REPORT$STRING(NFIL,'FIXED OCCUPATIONS')
         CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
         WRITE(NFIL,*)'OCCUPATIONS'
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             WRITE(NFIL,FMT='("FOR K-POINT: ",I5," AND SPIN ",I1)')IKPT,ISPIN
             CALL DYNOCC_PREIG('OCC',NFIL,NB,OCC(:,IKPT,ISPIN),WKPT(IKPT))
           ENDDO
         ENDDO
         RETURN
       ELSE
         CALL REPORT$STRING(NFIL,'OCCUPATIONS VARIABLE')
       END IF
       IF(TFIXTOT) THEN
         CALL REPORT$R8VAL(NFIL,'FIXED CHARGE',-(TOTCHA-SUMOFZ),'E')
         CALL REPORT$R8VAL(NFIL,'SUM OF NUCLEAR CHARGES',SUMOFZ,'E')
         CALL REPORT$R8VAL(NFIL,'#(ELECTRONS)',TOTCHA,' ')
       ELSE
         CALL REPORT$R8VAL(NFIL,'FIXED CHEMICAL POTENTIAL',TOTPOT/EV,'EV')
         CALL REPORT$R8VAL(NFIL,'SUM OF NUCLEAR CHARGES',SUMOFZ,'E')
       END IF
!
       IF(NSPIN.EQ.2) THEN
         IF(TFIXSPIN) THEN
           CALL REPORT$R8VAL(NFIL,'FIXED SPIN',SPINCHA*0.5D0,'HBAR')
         ELSE
           CALL REPORT$R8VAL(NFIL,'FIXED MAGNETIC FIELD',SPINPOT/EV &
      &                                                 ,'EV/(HBAR/2)')
         END IF
       END IF
!
       IF(TRESTARTFILEPRESENT) THEN
         IF(STARTTYPE.EQ.'N') THEN
           STRING='INITIAL STATES ARE FILLED BOTTOM UP ASSUMING FLAT BANDS'
           CALL REPORT$STRING(NFIL,TRIM(STRING))
         ELSE IF(STARTTYPE.EQ.'X') THEN
           STRING='INITIAL OCCUPATIONS ARE READ FROM RESTART FILE'
           CALL REPORT$STRING(NFIL,TRIM(STRING))
         ELSE IF(STARTTYPE.EQ.'E') THEN
           STRING='INITIAL OCCUPATIONS CONSTRUCTED FROM ENERGIES'
           STRING=TRIM(STRING)//' READ FROM RESTART FILE'
           CALL REPORT$STRING(NFIL,TRIM(STRING))
         ELSE
           CALL ERROR$MSG('VALUE OF STARTTYPE NOT RECOGNIZED')
           CALL ERROR$CHVAL('STARTTYPE',STARTTYPE)
           CALL ERROR$STOP('DYNOCC_REPORTSETTING')
         END IF
       ELSE
         STRING='NO MERMIN INFORMATION READ FROM RESTART FILE'
         CALL REPORT$STRING(NFIL,TRIM(STRING))
         STRING='INITIAL STATES ARE FILLED BOTTOM UP ASSUMING FLAT BANDS'
         CALL REPORT$STRING(NFIL,TRIM(STRING))
       END IF
!
       IF(BZITYPE.EQ.'SAMP') THEN
         CALL REPORT$STRING(NFIL,'BRILLOUIN-ZONE INTEGRATION BY SAMPLING') 
         CALL REPORT$R8VAL(NFIL,'TEMPERATURE',TEMP/KELVIN,'K')
         IF(TADIABATIC) THEN
           STRING='QUASI-ADIABATIC OCCUPATIONS USING MERMIN FUNCTIONAL'
           CALL REPORT$STRING(NFIL,TRIM(STRING))
           CALL REPORT$R8VAL(NFIL,'RETARDATION',RETARD,'DELTAT')
         ELSE 
           STRING='DYNAMICAL OCCUPATIONS USING MERMIN FUNCTIONAL'
           CALL REPORT$STRING(NFIL,TRIM(STRING))
           CALL REPORT$R8VAL(NFIL,'MASS OF X',MX,'A.U.')
           IF(ANNEX.NE.0.D0) THEN
             CALL REPORT$R8VAL(NFIL,'FRICTION',ANNEX,' ')
           END IF
           IF(TSTOP) THEN
             CALL REPORT$CHVAL(NFIL,'INITIAL VELOCITIES ARE SET TO ','ZERO')
           END IF
         END IF
!
!      == NOE TETRAHEDRON METHOD ===============================================
       ELSE IF(BZITYPE.EQ.'TETRA+') THEN
         STRING='BRILLOUIN-ZONE INTEGRATION BY TETRAHEDRON METHOD'
         CALL REPORT$STRING(NFIL,TRIM(STRING))
         STRING='QUASI-ADIABATIC OCCUPATIONS USING MERMIN FUNCTIONAL'
         CALL REPORT$STRING(NFIL,TRIM(STRING))
         CALL REPORT$R8VAL(NFIL,'RETARDATION',RETARD,'DELTAT')
         CALL REPORT$R8VAL(NFIL,'TEMPERATURE',TEMP/KELVIN,'K')
         IF(TEMP.GT.1.D-4) THEN
           CALL ERROR$MSG('TETRAHEDRON METHOD ONLY WORKS FOR T=0')
           CALL ERROR$STOP('DYNOCC_REPORTSETTING')
         END IF
       END IF
       RETURN
       END 
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC_PREIG(ID,NFIL,NB,EIG,UNIT)
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       INTEGER(4)  ,INTENT(IN) :: NFIL
       INTEGER(4)  ,INTENT(IN) :: NB
       REAL(8)     ,INTENT(IN) :: EIG(NB)
       REAL(8)     ,INTENT(IN) :: UNIT
       INTEGER(4)              :: ITEN,IB,I1,I2
       CHARACTER(32)           :: FORMAT
!      *****************************************************************
       IF(ID.EQ.'EIG') THEN
         FORMAT='("EIG",I3,":",10F8.3)'
       ELSE IF(ID.EQ.'OCC') THEN
         FORMAT='("OCC",I3,":",10F8.3)'
       ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$STOP('DYNOCC_PREIG')
       END IF
       ITEN=0
       DO WHILE (NB.GT.ITEN)
         I1=ITEN+1
         I2=MIN(ITEN+10,NB)
         WRITE(NFIL,FMT=FORMAT)ITEN,(EIG(IB)/UNIT,IB=I1,I2)
         ITEN=ITEN+10
       ENDDO
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC$ORDER(IKPT,ISPIN)
!      *****************************************************************
!      **                                                             **
!      **  REORDERS BANDS IN THE DYNOCC OBJECT ACCORDING TO           **
!      **  A PREDEFINED SETTING OF THE SORT OBJECT                    **
!      **                                                             **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: IKPT
       INTEGER(4),INTENT(IN) :: ISPIN
!      *****************************************************************
       IF(IKPT.LT.1.OR.IKPT.GT.NKPT) THEN
         CALL ERROR$MSG('IKPT OUT OF RANGE')
         CALL ERROR$I4VAL('IKPT',IKPT)
         CALL ERROR$STOP('DYNOCC$ORDER')
       END IF
       IF(ISPIN.LT.1.OR.ISPIN.GT.NSPIN) THEN
         CALL ERROR$MSG('SPIN OUT OF RANGE')
         CALL ERROR$I4VAL('ISPIN',ISPIN)
         CALL ERROR$STOP('DYNOCC$ORDER')
       END IF
!      == REORDER INPUT ================================================
       IF(TEPSILON) THEN
         CALL SORT$ORDERR8(1,NB,EPSILON(:,IKPT,ISPIN))
       END IF
       IF(ALLOCATED(MPSIDOT2)) THEN
         CALL SORT$ORDERR8(1,NB,MPSIDOT2(:,IKPT,ISPIN))
       END IF
!      == OCCUPATIONS ARE ONLY REORDERED IF OCCUPATIONS ARE DYNAMICAL
       IF(TDYN) THEN
         CALL SORT$ORDERR8(1,NB,XM(:,IKPT,ISPIN))
         CALL SORT$ORDERR8(1,NB,X0(:,IKPT,ISPIN))
         CALL SORT$ORDERR8(1,NB,XP(:,IKPT,ISPIN))
       END IF
       RETURN
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC_PRINTEPSILON(TEPS0,TEPSILON1,TWGHT)
!      *************************************************************************
!      ** ROUTINE IS USED AS HELPER FOR TESTING                               **
!      ** PRINTS EPSILON OR WGHT                                              **
!      **                                                                     **
!      **                                                                     **
!      *************************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       LOGICAL(4),INTENT(IN):: TEPS0
       LOGICAL(4),INTENT(IN):: TEPSILON1
       LOGICAL(4),INTENT(IN):: TWGHT
       REAL(8)              :: EV
       INTEGER(4)           :: IKPT,ISPIN
       INTEGER(4)           :: NFIL=6
!      *************************************************************************
       CALL CONSTANTS('EV',EV)
       DO ISPIN=1,NSPIN
         DO IKPT=1,NKPT
           WRITE(NFIL,FMT='("FOR K-POINT: ",I5," AND SPIN ",I1)')IKPT,ISPIN
           IF(TEPSILON1) THEN
             CALL DYNOCC_PREIG('EIG',NFIL,NB,EPSILON(:,IKPT,ISPIN),EV)
           END IF
           IF(TEPS0) THEN
             CALL DYNOCC_PREIG('EIG',NFIL,NB,EPS0(:,IKPT,ISPIN),EV)
           END IF
           IF(TWGHT) THEN
             CALL DYNOCC_PREIG('OCC',NFIL,NB,WGHT(:,IKPT,ISPIN),WKPT(IKPT))
           END IF
         ENDDO
       ENDDO
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC_WGHT()
!      *************************************************************************
!      ** CALCULATES INTEGRATION WEIGHTS FOR THE CURRENT TIME STEP            **
!      ** USES X0 FOR TADIABATIC=FALSE.                                       **
!      ** USES EPS0 FOR TADIABATIC=TRUE                                       **
!      **                                                                     **
!      ** IS ONLY EXECUTED ONCE PER TIME STEP.                                **
!      ** IS ONLY EXECUTED IF WGHT IS NOT ALLOCATED. WGHT IS                  **
!      ** DEALLOCATED IN DYNOCC$SWITCH                                        **
!      **                                                                     **
!      *************************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       INTEGER(4)          :: ISPIN,IKPT,IB
       REAL(8)             :: SVAR,DSVAR
       REAL(8)             :: XBAR(NB,NKPT,NSPIN)
!      *************************************************************************
       IF(ALLOCATED(WGHT)) RETURN
                            CALL TRACE$PUSH('DYNOCC_WGHT')
       ALLOCATE(WGHT(NB,NKPT,NSPIN))

!      =========================================================================
!      == CATCH CASE WITHOUT KNOWN ENERGIES                                   ==
!      == OCCURS IN THE FIRST TIME STEP, WHEN STARTING FROM SCRATCH           ==
!      =========================================================================
       IF(.NOT.ALLOCATED(EPS0)) THEN
!PRINT*,'WEIGHTS WITHOUT ALLOCATED EPS0========================='
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               CALL DYNOCC_FOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
               WGHT(IB,IKPT,ISPIN)=FMAX*WKPT(IKPT)*SVAR
             ENDDO
!PRINT*,'-------------------------WGHT ',IKPT,ISPIN,'-------------------------'
!WRITE(*,FMT='(10F8.3)')WGHT(:,IKPT,ISPIN)
           ENDDO
         ENDDO
                                   CALL TRACE$POP
         RETURN
       END IF
!
!      =========================================================================
!      ==                                                                     ==
!      =========================================================================
       IF(TADIABATIC) THEN
         IF(BZITYPE.EQ.'TETRA+') THEN
!!$PRINT*,'ENERGIES WITH ALLOCATED EPS0=== MIN=',MINVAL(EPS0),' MAX=',MAXVAL(EPS0)
!!$DO ISPIN=1,NSPIN
!!$  DO IKPT=1,NKPT
!!$    PRINT*,'-------------------------EPS0 ',IKPT,ISPIN,'-------------------------'
!!$    WRITE(*,FMT='(10F8.3)')EPS0(:,IKPT,ISPIN)
!!$  ENDDO
!!$ENDDO
           CALL DYNOCC_TETRAINTERFACE(NSPIN,NKPT,NB,TFIXTOT,TFIXSPIN &
      &               ,FMAX,TOTCHA,SPINCHA,TOTPOT,SPINPOT,EPS0,WGHT)
!!$PRINT*,'WEIGHTS WITH ALLOCATED EPS0========================='
!!$DO ISPIN=1,NSPIN
!!$  DO IKPT=1,NKPT
!!$PRINT*,'-------------------------WGHT ',IKPT,ISPIN,'-------------------------'
!!$WRITE(*,FMT='(10F8.3)')WGHT(:,IKPT,ISPIN)/WKPT(IKPT)
!!$  ENDDO
!!$ENDDO
         ELSE IF(BZITYPE.EQ.'SAMP') THEN
           CALL DYNOCC_INIOCCBYENERGY(NB,NKPT,NSPIN,FMAX &
      &          ,TEMP,TFIXTOT,TOTCHA,TOTPOT,TFIXSPIN,SPINCHA,SPINPOT &
      &          ,WKPT,EPS0,XBAR)
         ELSE 
           CALL ERROR$MSG('BZITYPE NOT RECOGNIZED')
           CALL ERROR$CHVAL('BZITYPE',BZITYPE)
           CALL ERROR$STOP('DYNOCC_WGHT')
         END IF
!
       ELSE 
         IF(BZITYPE.EQ.'TETRA+') THEN
           CALL ERROR$MSG('BZITYPE=TETRA+ WORKS ONLY WITH TADIABATIC=TRUE')
           CALL ERROR$CHVAL('BZITYPE',BZITYPE)
           CALL ERROR$STOP('DYNOCC_WGHT')
         ELSE IF(BZITYPE.EQ.'SAMP') THEN
           XBAR(:,:,:)=X0(:,:,:)
         ELSE 
           CALL ERROR$MSG('BZITYPE NOT RECOGNIZED')
           CALL ERROR$CHVAL('BZITYPE',BZITYPE)
           CALL ERROR$STOP('DYNOCC_WGHT')
         END IF
       END IF
!     
!
!      =================================================================
!      ==                                                             ==
!      =================================================================
       IF(BZITYPE.EQ.'SAMP') THEN
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               CALL DYNOCC_FOFX(XBAR(IB,IKPT,ISPIN),SVAR,DSVAR)
               WGHT(IB,IKPT,ISPIN)=FMAX*WKPT(IKPT)*SVAR
             ENDDO
           ENDDO
         ENDDO
       END IF
                                   CALL TRACE$POP
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC_EKIN(EKIN_)
!      *****************************************************************
!      **                                                             **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(OUT) :: EKIN_
       INTEGER(4)          :: ISPIN,IKPT
       REAL(8)             :: SVAR
!      *****************************************************************
       EKIN_=0.D0
       IF(BZITYPE.EQ.'SAMP'.AND.(.NOT.TADIABATIC)) THEN
         IF(.NOT.TPROPAGATED) THEN
           CALL ERROR$MSG('KINETIC ENERGY AVAILABLE ONLY AFTER PROPAGATION')
           CALL ERROR$STOP('DYNOCC_EKIN')
         END IF
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             SVAR=FMAX*WKPT(IKPT)*0.5D0*MX/(2.D0*DELTAT)**2
             EKIN_=EKIN_+SVAR*SUM((XP(:,IKPT,ISPIN)-XM(:,IKPT,ISPIN))**2)
           ENDDO
         ENDDO
       END IF
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC_EPOT(EPOT_)
!      *****************************************************************
!      **                                                             **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       REAL(8),INTENT(OUT) :: EPOT_
       INTEGER(4)          :: ISPIN,IKPT,IB
       REAL(8)             :: SVAR,DSVAR
       REAL(8)             :: Q0,S0,SIGMA
!      *****************************************************************
       EPOT=0.D0
       Q0=0.D0
       S0=0.D0
       IF(BZITYPE.EQ.'SAMP') THEN
         DO ISPIN=1,NSPIN
           SIGMA=DBLE(3-2*ISPIN)   ! SPIN DIRECTION       
           DO IKPT=1,NKPT
             DO IB=1,NB
               CALL DYNOCC_FOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
               Q0=Q0+FMAX*WKPT(IKPT)*SVAR
               S0=S0+FMAX*WKPT(IKPT)*SVAR*SIGMA
               CALL DYNOCC_SOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
               EPOT=EPOT-TEMP*FMAX*WKPT(IKPT)*SVAR
             ENDDO
           ENDDO
         ENDDO
       ELSE IF(BZITYPE.EQ.'TETRA+') THEN
         EPOT=0.D0
         IF(ALLOCATED(WGHT)) THEN
           Q0=SUM(WGHT)
           S0=0.D0
           IF(NSPIN.EQ.2) THEN
             S0=SUM(WGHT(:,:,1)-WGHT(:,:,2))
           END IF
         ELSE    ! BEST GUESS
           Q0=TOTCHA
           S0=0.D0
         END IF
       ELSE
         CALL ERROR$MSG('BZITYPE NOT RECOGNIZED')
         CALL ERROR$STOP('DYNOCC_EPOT')
       END IF
       IF(TFIXTOT) THEN
         EPOT=EPOT-TOTPOT*(Q0-TOTCHA)
       ELSE
         EPOT=EPOT-TOTPOT*(Q0-SUMOFZ)
       END IF
       IF(.NOT.TFIXSPIN) THEN
         EPOT=EPOT-SPINPOT*(S0-SPINCHA)
       ELSE
         EPOT=EPOT-SPINPOT*S0
       END IF
       EPOT_=EPOT
       RETURN 
       END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC$PROPAGATE()
!      *************************************************************************
!      ** PROPAGATE OCCUPATION VARIABLES                                      **
!      **                                                                     **
!      ** TWO DIFFERENT MODES OF OPERATION:                                   **
!      ** 1)TADIABATIC=TRUE: USES EPSM/0/P AS VARIABLES AND PROPAGATES        **
!      **     USING RETARDATION                                               **
!      **                                                                     **
!      ** 2)TADIABATIC=.FALSE: USES XM/0/P AS VARIABLES AND PROPAGATES        **
!      **     DYNAMICALLY USING MERMIN FUNCTIONAL                             **
!      **                                                                     **
!      *************************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       REAL(8)   ,PARAMETER :: DSMALL=1.D-5
       REAL(8)   ,PARAMETER :: QTOL=1.D-10
       INTEGER(4),PARAMETER :: NITER=1000
       REAL(8)              :: FX(NB,NKPT,NSPIN)
       REAL(8)              :: WGHT2(NB,NKPT,NSPIN)
       REAL(8)              :: XBAR(NB,NKPT,NSPIN)
       INTEGER(4)           :: ISPIN,IKPT,IB
       REAL(8)              :: SVAR,OCC,DOCCDX,FORCE,DSVAR
       REAL(8)              :: SVAR1,SVAR2,SVAR3
       REAL(8)              :: SIGMA    ! SPIN DIRECTION, I.E.: (+/-)1.D0
       REAL(8)              :: S0,S1
       REAL(8)              :: Q0,Q1
       LOGICAL(4)           :: TMPSIDOT2
       LOGICAL(4)           :: TCONV
       INTEGER(4)           :: ITER
       REAL(8)              :: EREL
       REAL(8)              :: DEKIN
       REAL(8)              :: DF(NB,NKPT,NSPIN)
!      *************************************************************************
       EKIN=0.D0
       EPOT=0.D0
       TPROPAGATED=.TRUE.
       IF(.NOT.TDYN) THEN 
         XP(:,:,:)=X0(:,:,:)
         RETURN
       END IF
                              CALL TRACE$PUSH('DYNOCC$PROPAGATE')
!
!      =========================================================================
!      ==  CHECK DATA                                                         ==
!      =========================================================================
       IF(NSPIN.EQ.1) THEN
         TFIXSPIN=.FALSE.
         SPINPOT=0.D0
       END IF
       IF(.NOT.TEPSILON) THEN
         CALL ERROR$MSG('ONE-PARTICLE ENERGIES HAVE NOT BEEN SET')
         CALL ERROR$STOP('DYNOCC$PROPAGATE')
       END IF
       TMPSIDOT2=ALLOCATED(MPSIDOT2)
!
!      == FILL IN ENERGIES IF IT DID NOT HAPPEN ================================
       IF(.NOT.ALLOCATED(EPS0)) THEN
         ALLOCATE(EPSM(NB,NKPT,NSPIN))
         ALLOCATE(EPS0(NB,NKPT,NSPIN))
         ALLOCATE(EPSP(NB,NKPT,NSPIN))
         EPS0=EPSILON
         EPSM=EPSILON
         EPSP=EPSILON
         IF(TMPSIDOT2) THEN
           EPSM=EPSM+MPSIDOT2
           EPS0=EPS0+MPSIDOT2
           EPSP=EPSP+MPSIDOT2
         END IF
       END IF
!
!      =========================================================================
!      ==  TRANSFORM TRAJECTORY BACK INTO THE INTERVAL [0,1] BY               ==
!      ==  TRANSLATION BY 2 AND MIRROR AT X=0 AND X=1                         ==
!      =========================================================================
       DO ISPIN=1,NSPIN
         DO IKPT=1,NKPT
           DO IB=1,NB
             SVAR=X0(IB,IKPT,ISPIN)
             IF(SVAR.GE.0.D0.AND.SVAR.LT.1.D0) CYCLE
!            == SHIFT TRAJECTORY INTO THE INTERVAL [0,2]
             IF(SVAR.LT.0.D0.OR.SVAR.GE.2.D0) THEN
               SVAR=REAL(2*(INT(0.5D0*SVAR+100.D0)-100),KIND=8)
               X0(IB,IKPT,ISPIN)=X0(IB,IKPT,ISPIN)-SVAR
               XM(IB,IKPT,ISPIN)=XM(IB,IKPT,ISPIN)-SVAR
             END IF
!            == MIRROR TRAJECTORY AT X=1
             IF(X0(IB,IKPT,ISPIN).GT.1.D0) THEN
               X0(IB,IKPT,ISPIN)=2.D0-X0(IB,IKPT,ISPIN)
               XM(IB,IKPT,ISPIN)=2.D0-XM(IB,IKPT,ISPIN)
             END IF
           ENDDO
         ENDDO
       ENDDO
!
!      =========================================================================
!      ==  ESTIMATE TOTPOT FROM T=0                                           ==
!      =========================================================================
       CALL DYNOCC_WGHT()
       IF(BZITYPE.EQ.'TETRA+') THEN
         CALL DYNOCC_TETRAINTERFACE(NSPIN,NKPT,NB,TFIXTOT,TFIXSPIN &
      &               ,FMAX,TOTCHA,SPINCHA,TOTPOT,SPINPOT,EPSILON,WGHT2)
       ELSE IF(BZITYPE.EQ.'SAMP') THEN
         CALL DYNOCC_INIOCCBYENERGY(NB,NKPT,NSPIN,FMAX &
      &          ,TEMP,TFIXTOT,TOTCHA,TOTPOT,TFIXSPIN,SPINCHA,SPINPOT &
      &          ,WKPT,EPSILON,WGHT2)
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               CALL DYNOCC_FOFX(WGHT2(IB,IKPT,ISPIN),SVAR,DSVAR)
               WGHT2(IB,IKPT,ISPIN)=FMAX*WKPT(IKPT)*SVAR
             ENDDO
           ENDDO
         ENDDO
       ELSE 
         CALL ERROR$MSG('BZITYPE NOT RECOGNIZED')
         CALL ERROR$CHVAL('BZITYPE',BZITYPE)
         CALL ERROR$STOP('DYNOCC$PROPAGATE')
       END IF
       DF(:,:,:)=WGHT2(:,:,:)-WGHT(:,:,:)
       MAXDEVOCC=MAXVAL(ABS(DF))
       DEVENERGY=SUM(DF*EPSILON)
!
!      =========================================================================
!      ==  PROPAGATE ENERGIES FOR ADIABATIC CALCULATION                       ==
!      =========================================================================
       EPSP=EPS0+(EPSILON-EPS0)/RETARD
!
IF(MAXVAL(EPSP).GT.1.D+2.OR.MINVAL(EPSP).LT.-1.D+2) THEN
  CALL ERROR$MSG('EXTREME VALUES ARE INDICATIVE OF ERROR')
  CALL ERROR$R8VAL('RETARD',RETARD)
  CALL ERROR$R8VAL('MIN(MPSIDOT2)',MINVAL(MPSIDOT2))
  CALL ERROR$R8VAL('MAX(MPSIDOT2)',MAXVAL(MPSIDOT2))
  CALL ERROR$R8VAL('MIN(EPSILON)',MINVAL(EPSILON))
  CALL ERROR$R8VAL('MAX(EPSILON)',MAXVAL(EPSILON))
  CALL ERROR$R8VAL('MIN(EPS0)',MINVAL(EPS0))
  CALL ERROR$R8VAL('MAX(EPS0)',MAXVAL(EPS0))
  CALL ERROR$R8VAL('MIN(EPSP)',MINVAL(EPSP))
  CALL ERROR$R8VAL('MAX(EPSP)',MAXVAL(EPSP))
  CALL ERROR$STOP('DYNOCC$PROPAGATE')
END IF
!      IF(TMPSIDOT2)EPSP=EPSP+MPSIDOT2/RETARD
!
!      =========================================================================
!      ==  RETURN FOR QUASI ADIABATIC CALCULATION                             ==
!      =========================================================================
       IF(TADIABATIC) THEN
!
!        =======================================================================
!        ==  DETERMINE OCCUPATION VARIABLES FOR THE NEXT TIME STEP            ==
!        =======================================================================
         IF(BZITYPE.EQ.'SAMP') THEN
           CALL DYNOCC_INIOCCBYENERGY(NB,NKPT,NSPIN,FMAX &
      &          ,TEMP,TFIXTOT,TOTCHA,TOTPOT,TFIXSPIN,SPINCHA,SPINPOT &
      &          ,WKPT,EPSP,XP)
         END IF
                              CALL TRACE$POP()
         RETURN
       END IF
!
!      =========================================================================
!      == HERE THE ADIABATIC OPTION IS FINISHED.THE FOLLOWING IS FOR THE      ==
!      == DYNAMICAL OCCUPATIONS USING THE MERMIN FUNCTIONAL                   ==
!      =========================================================================
!
!      =========================================================================
!      ==  FREEZE/UNFREEZE OCCUPATIONS                                        ==
!      =========================================================================
       DEKIN=0.D0
       DO ISPIN=1,NSPIN
         SIGMA=DBLE(3-2*ISPIN)   ! SPIN DIRECTION       
         DO IKPT=1,NKPT
!!$PRINT*,'TF',IKPT,ISPIN,':',TFROZEN(:,IKPT,ISPIN)
!!$PRINT*,'X0',IKPT,ISPIN,':',X0(:,IKPT,ISPIN)
           DO IB=1,NB
!  THE FACTOR 5 IS A CHOSEN PARAMETER
             EREL=(EPSILON(IB,IKPT,ISPIN)-TOTPOT-SIGMA*SPINPOT)/(5.D0*TEMP)
             IF(TFROZEN(IB,IKPT,ISPIN)) THEN
               IF(ABS(EREL).LT.1.D0) THEN
!                == UNFREEZE ALL STATES IN THE WINDOW AROUND EFERMI     
                 TFROZEN(IB,IKPT,ISPIN)=.FALSE.
!!$PRINT*,'UNFREEZE STATE A',IB,IKPT,ISPIN,X0(IB,IKPT,ISPIN),XM(IB,IKPT,ISPIN)
               ELSE
!                == UNFREEZE OCCUPIED STATES ABOVE EFERMI OR UNOCCUPIED BELOW =====
                 SVAR=(X0(IB,IKPT,ISPIN)-0.5D0)*EREL
!!$IF(SVAR.GT.0.D0)PRINT*,'UNFREEZE STATE B',IB,IKPT,ISPIN,X0(IB,IKPT,ISPIN),XM(IB,IKPT,ISPIN)
                 IF(SVAR.GT.0.D0)TFROZEN(IB,IKPT,ISPIN)=.FALSE.
               END IF
             ELSE
!              == DO NOT FREEZE STATES IN THE WINDOW
               IF(ABS(EREL).LT.1.D0) CYCLE
!              == DO NOT FREEZE STATES ON THE WRONG SIDE OF THE THE WINDOW
               SVAR=(X0(IB,IKPT,ISPIN)-0.5D0)*EREL
               IF(SVAR.GT.0.D0) CYCLE
!              == FREEZE STATES WHEN THEY PASS THROUGH INTEGER OCCUPATION
!              == EMPTY STATES  ===========================================
               SVAR=X0(IB,IKPT,ISPIN)*XM(IB,IKPT,ISPIN)
               IF(SVAR.LT.0.AND.EREL.GT.1.D0) THEN
!                == CHECK SIZE OF X0 TO AVOID STEPS IN POTENTIAL ENERGY =====
                 IF(ABS(X0(IB,IKPT,ISPIN)).LT.1.D-5) THEN 
                   TFROZEN(IB,IKPT,ISPIN)=.TRUE.
                   DEKIN=DEKIN+FMAX*WKPT(IKPT)*0.5D0*MX &
     &                        *((X0(IB,IKPT,ISPIN)-XM(IB,IKPT,ISPIN))/DELTAT)**2
!!$PRINT*,'FREEZE STATE A',IB,IKPT,ISPIN,X0(IB,IKPT,ISPIN),XM(IB,IKPT,ISPIN)
                   X0(IB,IKPT,ISPIN)=0.D0
                   XM(IB,IKPT,ISPIN)=0.D0
                 END IF
               END IF
!              == FILLED STATES ===========================================
               SVAR=(X0(IB,IKPT,ISPIN)-1.D0)*(XM(IB,IKPT,ISPIN)-1.D0)
               IF(SVAR.LT.0.AND.EREL.LE.-1.D0) THEN
!                == CHECK SIZE OF X0 TO AVOID STEPS IN POTENTIAL ENERGY =====
                 IF(ABS(X0(IB,IKPT,ISPIN)-1.D0).LT.1.D-5) THEN 
                    TFROZEN(IB,IKPT,ISPIN)=.TRUE.
                    DEKIN=DEKIN+FMAX*WKPT(IKPT)*0.5D0*MX &
     &                      *((X0(IB,IKPT,ISPIN)-XM(IB,IKPT,ISPIN))/DELTAT)**2
!!$PRINT*,'FREEZE STATE B',IB,IKPT,ISPIN,X0(IB,IKPT,ISPIN)-1.D0,XM(IB,IKPT,ISPIN)-1.D0
                    X0(IB,IKPT,ISPIN)=1.D0
                    XM(IB,IKPT,ISPIN)=1.D0
                  END IF
               END IF
             END IF
           ENDDO
         ENDDO
       ENDDO
!
!      == DETERMINE KINETIC ENERGY OF UNFROZEN OCCUPATIONS =====================
       EKIN=0.D0
       DO ISPIN=1,NSPIN
         DO IKPT=1,NKPT
           SVAR=FMAX*WKPT(IKPT)*0.5D0*MX/DELTAT**2
           DO IB=1,NB
             IF(TFROZEN(IB,IKPT,ISPIN)) CYCLE
             EKIN=EKIN+SVAR*(X0(IB,IKPT,ISPIN)-XM(IB,IKPT,ISPIN))**2
           ENDDO
         ENDDO
       ENDDO
!
!      == FEED KINETIC ENERGY BACK  ==========================================
       IF(DEKIN.GT.0.D0.AND.EKIN.GT.0.D0) THEN
         SVAR=SQRT(1.D0+DEKIN/EKIN)
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               IF(TFROZEN(IB,IKPT,ISPIN)) CYCLE
               XM(IB,IKPT,ISPIN)=X0(IB,IKPT,ISPIN) &
                                +(XM(IB,IKPT,ISPIN)-X0(IB,IKPT,ISPIN))*SVAR
             ENDDO
           ENDDO
         ENDDO
         DEKIN=0.D0
       END IF
!
!      =================================================================
!      ==  OFFSET OCCUPATIONS SLIGHTLY AWAY FROM ZERO                 ==
!      =================================================================
       DO ISPIN=1,NSPIN
         DO IKPT=1,NKPT
           DO IB=1,NB
             IF(TFROZEN(IB,IKPT,ISPIN)) CYCLE
             SVAR=0.5D0*ABS(X0(IB,IKPT,ISPIN)-XM(IB,IKPT,ISPIN))
             IF(SVAR.GT.1.D-6) CYCLE ! NO CORRECTION IF NONZERO VELOCITY 
             IF(ABS(X0(IB,IKPT,ISPIN)).LT.DSMALL**2) THEN
               X0(IB,IKPT,ISPIN)=DSMALL
               XM(IB,IKPT,ISPIN)=DSMALL
             ELSE IF(ABS(X0(IB,IKPT,ISPIN)-1.D0).LT.DSMALL**2) THEN
               X0(IB,IKPT,ISPIN)=1.D0-DSMALL
               XM(IB,IKPT,ISPIN)=1.D0-DSMALL
             END IF
           ENDDO
         ENDDO
       ENDDO
!
!      =================================================================
!      ==  SET VELOCITY OF OCCUPATIONS TO ZERO                        ==
!      =================================================================
       IF(TSTOP) THEN
         XM(:,:,:)=X0(:,:,:)
         TSTOP=.FALSE.
       END IF
!
!      =================================================================
!      ==  CALCULATE TOTAL ENERGY AND FORCE ON X                      ==
!      ==  REMARK:  THE FACTOR FMAX IS NOT USED BECAUSE THE SAME TERM ==
!      ==         APPEARS ALSO IN THE KINETIC ENERGY                  ==
!      =================================================================
       EPOT=0.D0
       FX(:,:,:)=0.D0
       DO ISPIN=1,NSPIN 
         DO IKPT=1,NKPT
           DO IB=1,NB
             IF(TFROZEN(IB,IKPT,ISPIN)) CYCLE
!            == FORCE FROM BANDS ======================================
             CALL DYNOCC_FOFX(X0(IB,IKPT,ISPIN),OCC,DOCCDX)
             FORCE=-EPSILON(IB,IKPT,ISPIN)*DOCCDX
!            == NOTE, THAT THIS IS NOT THE EULER LAGRANGE EQUATION   ==
!            == THEREFORE THE MINUS SIGN IN THE FOLLOWING LINE       ==
             IF(TMPSIDOT2) FORCE=FORCE-MPSIDOT2(IB,IKPT,ISPIN)*DOCCDX
!            == FORCE FROM ENTROPY TERM (-TS) =========================
             CALL DYNOCC_SOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
             EPOT=EPOT-FMAX*WKPT(IKPT)*TEMP*SVAR
             FORCE=FORCE-TEMP*DSVAR
             FX(IB,IKPT,ISPIN)=FORCE
           ENDDO
         ENDDO
       ENDDO
!
!      =================================================================
!      ==  PROPAGATE X WITOUT CONSTRAINTS                             ==
!      ==  CONSTRAINT FORCE FX                                        ==
!      ==  XP=XBAR+FX*(TOTPOT+SIGMA*SPINPOT)                          ==
!      =================================================================
       XBAR(:,:,:)=X0(:,:,:)
       SVAR1=2.D0/(1.D0+ANNEX)
       SVAR2=1.D0-SVAR1
       SVAR3=(SVAR1/2.D0)*DELTAT**2/MX
       DO ISPIN=1,NSPIN
         DO IKPT=1,NKPT
           DO IB=1,NB
             IF(TFROZEN(IB,IKPT,ISPIN)) CYCLE
             XBAR(IB,IKPT,ISPIN)=SVAR1*X0(IB,IKPT,ISPIN) &
       &                        +SVAR2*XM(IB,IKPT,ISPIN) &
       &                        +SVAR3*FX(IB,IKPT,ISPIN)
!            == FORCE OF CONSTRAINT ==================================
             CALL DYNOCC_FOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
             FX(IB,IKPT,ISPIN)=DSVAR*SVAR3
           ENDDO
         ENDDO
       ENDDO
!
!      =================================================================
!      ==  APPLY CONSTRAINTS                                          ==
!      =================================================================
       DO ITER=1,NITER
!        == EVALUATE Q(X)=Q0+Q1*FX*ALPHA;  S(X)=Q0+S1*FX*BETA ==========
         Q0=0.D0
         S0=0.D0
         Q1=0.D0 
         S1=0.D0
         DO ISPIN=1,NSPIN
           SIGMA=DBLE(3-2*ISPIN)
           DO IKPT=1,NKPT
             DO IB=1,NB
               FORCE=FX(IB,IKPT,ISPIN)
               XP(IB,IKPT,ISPIN)=XBAR(IB,IKPT,ISPIN) &
     &                         +(TOTPOT+SIGMA*SPINPOT)*FORCE
               CALL DYNOCC_FOFX(XP(IB,IKPT,ISPIN),SVAR,DSVAR)
               SVAR1=SVAR*FMAX*WKPT(IKPT)
               Q0=Q0+SVAR1 
               S0=S0+SVAR1*SIGMA
               SVAR1=DSVAR*FMAX*WKPT(IKPT)*FORCE
               Q1=Q1+SVAR1 
               S1=S1+SVAR1*SIGMA
              ENDDO
           ENDDO
         ENDDO
         IF(NSPIN.EQ.1)THEN
            S0=0.D0
            S1=0.D0
         END IF
!
!        ===============================================================
!        == CHECK CONVERGENCE ==========================================
!        ===============================================================
         TCONV=.TRUE.
         IF(TFIXTOT) TCONV=TCONV.AND.(ABS(Q0-TOTCHA).LT.QTOL) 
         IF(TFIXSPIN)TCONV=TCONV.AND.(ABS(S0-SPINCHA).LT.QTOL)
         IF(TCONV) EXIT
!
!        ===============================================================
!        == FIND ZERO OF LINARIZED Q(X)-TOTCHA; S(X)-SPINCHA ===========
!        ==  Q0+Q1*DTOTPOT+S1*DSPINPOT=0                              ==
!        ==  S0+S1*DTOTPOT+Q1*DSPINPOT=0                              ==
!        ===============================================================
         IF(TFIXTOT.AND.TFIXSPIN) THEN
           TOTPOT =TOTPOT -(+Q1*(Q0-TOTCHA)-S1*(S0-SPINCHA))/(Q1**2-S1**2)
           SPINPOT=SPINPOT-(-S1*(Q0-TOTCHA)+Q1*(S0-SPINCHA))/(Q1**2-S1**2)
         ELSE IF(TFIXTOT.AND.(.NOT.TFIXSPIN)) THEN
           TOTPOT=TOTPOT+(TOTCHA-Q0)/Q1
           SPINPOT=0.D0
         ELSE IF((.NOT.TFIXTOT).AND.TFIXSPIN) THEN
           SPINPOT=SPINPOT+(SPINCHA-S0)/S1
           TOTPOT=0.D0
         END IF
       ENDDO
       IF(.NOT.TCONV) THEN
         CALL ERROR$MSG('LOOP NOT CONVERGED') 
         CALL ERROR$R8VAL('TOTAL CHARGE',Q0)
         CALL ERROR$R8VAL('DEVIATION TOTAL CHARGE',Q0-TOTCHA)
         CALL ERROR$MSG('IF Q1 TOO SMALL, INSTABILITY BECAUSE BANDGAP/T TOO LARGE')
         CALL ERROR$MSG('IN THIS CASE START WITH !MERMIN!START=T OR INCREASE T[K]')
         CALL ERROR$R8VAL('Q1',Q1)
         CALL ERROR$L4VAL('TFIXTOT',TFIXTOT)
         CALL ERROR$L4VAL('TFIXSPIN',TFIXSPIN)
         IF(NSPIN.NE.1)CALL ERROR$R8VAL('DEVIATION TOTAL SPIN',S0)
         CALL ERROR$STOP('DYNOCC$PROPAGATE')
       END IF
       IF(.NOT.TFIXTOT)TOTCHA=Q0
       IF(.NOT.TFIXSPIN)SPINCHA=S0
!
!      ===============================================================
!      == DETERMINE KINETIC ENERGY                                  ==
!      ===============================================================
       CALL DYNOCC_EKIN(EKIN)
                             CALL TRACE$POP
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DYNOCC$SWITCH()
!     ******************************************************************
!     ******************************************************************
      USE DYNOCC_MODULE
      IMPLICIT NONE
!     ******************************************************************
      TPROPAGATED=.FALSE.
      TEPSILON=.FALSE.
      IF(ALLOCATED(MPSIDOT2)) DEALLOCATE(MPSIDOT2)
      IF(ALLOCATED(WGHT)) DEALLOCATE(WGHT)
!
!     == PROPAGATE OCCUPATION VARIABLES
      IF(TDYN) THEN
        XM(:,:,:)=X0(:,:,:)
        X0(:,:,:)=XP(:,:,:)
        XP(:,:,:)=0.D0
        EPSM(:,:,:)=EPS0(:,:,:)
        EPS0(:,:,:)=EPSP(:,:,:)
        EPSP(:,:,:)=0.D0
      END IF
      EKIN=0.D0
      EPOT=0.D0
      RETURN
      END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC_TETRAINTERFACE(NSPIN,NKPT,NB,TFIXTOT,TFIXSPIN &
      &                ,FMAX,TOTCHA,SPINCHA,TOTPOT,SPINPOT,EPSILON,WGHT)
!      *************************************************************************
!      **  EVALUATES THE INTEGRATION WEIGHTS FROM KNOWN ENERGY EIGENVALUES    **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)   :: NSPIN
       INTEGER(4),INTENT(IN)   :: NKPT
       INTEGER(4),INTENT(IN)   :: NB
       LOGICAL(4),INTENT(IN)   :: TFIXTOT
       LOGICAL(4),INTENT(IN)   :: TFIXSPIN
       REAL(8)   ,INTENT(IN)   :: FMAX
       REAL(8)   ,INTENT(INOUT):: TOTCHA
       REAL(8)   ,INTENT(INOUT):: SPINCHA
       REAL(8)   ,INTENT(INOUT):: TOTPOT
       REAL(8)   ,INTENT(INOUT):: SPINPOT
       REAL(8)   ,INTENT(IN)   :: EPSILON(NB,NKPT,NSPIN)
       REAL(8)   ,INTENT(OUT)  :: WGHT(NB,NKPT,NSPIN)
       REAL(8)                 :: WORK1(NB*NSPIN,NKPT)
       REAL(8)                 :: WORK2(NB*NSPIN,NKPT)
       REAL(8)                 :: SIGMA
       REAL(8)                 :: Q0,V0
       INTEGER(4)              :: ISPIN
!      *************************************************************************
                           CALL TRACE$PUSH('DYNOCC_TETRAINTERFACE')
                           CALL TIMING$CLOCKON('TETRAHEDRON METHOD')
       IF(NSPIN.EQ.1) THEN
         IF(TFIXTOT) THEN
           SPINCHA=0.D0
           SPINPOT=0.D0
           CALL BRILLOUIN$DOS(NB,NKPT,EPSILON,WGHT,TOTCHA/FMAX,TOTPOT)
           WGHT=WGHT*FMAX
         ELSE
           CALL ERROR$MSG('OPTION NOT IMPLEMENTED')
           CALL ERROR$MSG('IF NSPIN=1, TFIXTOT MUST BE TRUE')
           CALL ERROR$STOP('DYNOCC_TETRAINTERFACE')
         END IF
       ELSE IF(NSPIN.EQ.2) THEN
         IF(TFIXTOT.AND.TFIXSPIN) THEN
           TOTPOT=0.D0
           SPINPOT=0.D0
           DO ISPIN=1,NSPIN
             SIGMA=REAL(3-2*ISPIN,KIND=8)
             WORK1(:NB,:)=EPSILON(:,:,ISPIN)
             Q0=0.5D0*(TOTCHA+SIGMA*SPINCHA)
             CALL BRILLOUIN$DOS(NB,NKPT,WORK1(1:NB,:),WORK2(:NB,:),Q0,V0)
             TOTPOT=TOTPOT+0.5D0*V0
             SPINPOT=SPINPOT-0.5D0*SIGMA*V0
             WGHT(:,:,ISPIN)=FMAX*WORK2(:NB,:)
           ENDDO
!
         ELSE IF(TFIXTOT.AND..NOT.TFIXSPIN) THEN
           WORK1(1:NB,:)     =EPSILON(:,:,1)+0.5D0*SPINPOT
           WORK1(NB+1:2*NB,:)=EPSILON(:,:,2)-0.5D0*SPINPOT
           CALL BRILLOUIN$DOS(NB*NSPIN,NKPT,WORK1,WORK2,TOTCHA,TOTPOT)
           WGHT(:,:,1)=WORK2(1:NB,:)
           WGHT(:,:,2)=WORK2(NB+1:2*NB,:)
           SPINCHA=SUM(WGHT(:,:,1)-WGHT(:,:,2))!

         ELSE IF(.NOT.TFIXTOT.AND.TFIXSPIN) THEN
           CALL ERROR$MSG('OPTION NOT IMPLEMENTED')
           CALL ERROR$MSG('OPTION: NSPIN=2, TFIXTOT=F TFIXSPIN=T')
           CALL ERROR$STOP('DYNOCC_TETRAINTERFACE')
!
         ELSE IF(.NOT.TFIXTOT.AND..NOT.TFIXSPIN) THEN
           CALL ERROR$MSG('OPTION NOT IMPLEMENTED')
           CALL ERROR$MSG('OPTION: NSPIN=2, TFIXTOT=F TFIXSPIN=F')
           CALL ERROR$STOP('DYNOCC_TETRAINTERFACE')
         END IF
       ELSE 
         CALL ERROR$MSG('NSPIN MUST HAVE VALUE 1 OR 2')
         CALL ERROR$STOP('DYNOCC_TETRAINTERFACE')
       END IF
                           CALL TIMING$CLOCKOFF('TETRAHEDRON METHOD')
                           CALL TRACE$POP()
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC$INIOCC
!      *****************************************************************
!      ** CHOOSE INITIAL CONDITIONS FOR THE OCCUPATIONS               **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
!      *****************************************************************
                                  CALL TRACE$PUSH('DYNOCC$INIOCC')
       CALL DYNOCC_INIOCCBYNUMBER(NB,NKPT,NSPIN,FMAX,TOTCHA,SPINCHA,WKPT,X0)
       XM(:,:,:)=X0(:,:,:)
                                   CALL TRACE$POP
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC$MODOCC(IB_,IKPT_,ISPIN_,OCC_)
!      *****************************************************************
!      **  MODIFY INDIVIDUAL OCCUPATIONS                              **
!      **                                                             **
!      **  OBTAINS THE NEW OCCUPATION FOR A PARTICULAR STATE          **
!      **  AND ADJUSTS TOTCHA, SPINCHA AND THE OCCUPATION VARIABLE    **
!      ** FOR THAT STATE                                              **
!      **                                                             **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: IB_,IKPT_,ISPIN_
       REAL(8)   ,INTENT(IN) :: OCC_
       REAL(8)               :: SVAR,DSVAR
!      *****************************************************************
       IF(IB_.GT.NB.OR.IKPT_.GT.NKPT.OR.ISPIN_.GT.NSPIN) THEN
         CALL ERROR$MSG('DIMENSIONS EXCEEDED')
         CALL ERROR$I4VAL('IB_',IB_)
         CALL ERROR$I4VAL('IKPT_',IKPT_)
         CALL ERROR$I4VAL('ISPIN_',ISPIN_)
         CALL ERROR$I4VAL('NB',NB)
         CALL ERROR$I4VAL('NKPT',NKPT)
         CALL ERROR$I4VAL('NSPIN',NSPIN)
         CALL ERROR$STOP('DYNOCC$MODOCC')
       END IF
!
!      =================================================================
!      ==  ADJUST TOTCHA AND SPINCHA                                  ==
!      =================================================================
       CALL DYNOCC_FOFX(X0(IB_,IKPT_,ISPIN_),SVAR,DSVAR)
       TOTCHA=TOTCHA+OCC_-SVAR*FMAX
       IF(NSPIN.EQ.2) THEN
         SPINCHA=SPINCHA+(OCC_-SVAR*FMAX)*DBLE(3-2*ISPIN_)
       END IF
!
!      =================================================================
!      ==  RESET OCCUPATIONS                                          ==
!      =================================================================
       CALL DYNOCC_XOFF(OCC_/FMAX,SVAR,DSVAR)
       X0(IB_,IKPT_,ISPIN_)=SVAR
       XM(IB_,IKPT_,ISPIN_)=SVAR
       RETURN 
       END
!PB031019START
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC_INIOCCBYNUMBER(NB,NKPT,NSPIN,FMAX &
      &                                ,TOTCHA,SPINCHA,WKPT,X)
!      *****************************************************************
!      ** CHOOSE INITIAL CONDITIONS FOR THE OCCUPATIONS               **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4) ,INTENT(IN) :: NB           ! #(BANDS),MAX
       INTEGER(4) ,INTENT(IN) :: NKPT         ! #(K-POINTS),MAX
       INTEGER(4) ,INTENT(IN) :: NSPIN        ! #(SPINS),MAX
       REAL(8)    ,INTENT(IN) :: TOTCHA       ! TOTAL CHARGE
       REAL(8)    ,INTENT(IN) :: SPINCHA      ! SPIN CHARGE
       REAL(8)    ,INTENT(IN) :: FMAX         ! MAX(OCCUPATION PER STATE)
       REAL(8)    ,INTENT(IN) :: WKPT(NKPT)   ! K-POINT WEIGHTS
       REAL(8)    ,INTENT(OUT):: X(NB,NKPT,NSPIN) ! OCCUPATION VARIABLE
       REAL(8)                :: TOT,SVAR
       INTEGER(4)             :: ISPIN,ISVAR
       REAL(8)                :: DX
!      *****************************************************************
       IF(TOTCHA.LT.0.D0.OR.ABS(SPINCHA).GT.TOTCHA) THEN
         CALL ERROR$MSG('NEGATIVE AMOUNT OF VALENCE ELECTRONS IS NOT ALLOWED')
         CALL ERROR$MSG('AND SPIN DENSITY MUST NOT EXCEED TOTAL DENSITY')
         CALL ERROR$R8VAL('TOTCHA',TOTCHA)
         CALL ERROR$R8VAL('SPINCHA',SPINCHA)
         CALL ERROR$STOP('DYNOCC_INIOCCBYNUMBER')
       END IF
!
       IF(NSPIN.EQ.1) THEN ! NON-SPIN-POLARIZED AND NON-COLLINEAR CASE
         ISVAR=INT(TOTCHA/FMAX)
         X(1:ISVAR,:,1)=1.D0
         X(ISVAR+1:NB,:,1)=0.D0
         TOT=(TOTCHA-FMAX*REAL(ISVAR,KIND=8))/FMAX
         IF(TOT.NE.0.D0) THEN
           CALL DYNOCC_XOFF(TOT,SVAR,DX)
           IF(ISVAR+1.GT.NB) THEN
             CALL ERROR$MSG('INSUFFICIENT NUMBER OF BANDS')
             CALL ERROR$I4VAL('NB',NB)
             CALL ERROR$R8VAL('TOTCHA',TOTCHA)
             CALL ERROR$STOP('DYNOCC_INIOCCBYNUMBER')
           END IF
           X(ISVAR+1,:,1)=SVAR
         END IF
       ELSE        ! SPIN-POLARIZED CASE
         DO ISPIN=1,NSPIN
           TOT=0.5D0*TOTCHA+0.5D0*SPINCHA*REAL(3-2*ISPIN,8)         
           ISVAR=INT(TOT)
           X(1:ISVAR,:,ISPIN)=1.D0
           X(ISVAR+1:NB,:,ISPIN)=0.D0
           TOT=TOT/FMAX-REAL(ISVAR,KIND=8)
           IF(TOT.NE.0.D0) THEN
             CALL DYNOCC_XOFF(TOT,SVAR,DX)
             IF(ISVAR+1.GT.NB) THEN
               CALL ERROR$MSG('INSUFFICIENT NUMBER OF BANDS')
               CALL ERROR$I4VAL('NB',NB)
               CALL ERROR$I4VAL('NSPIN',NSPIN)
               CALL ERROR$R8VAL('TOTCHA',TOTCHA)
               CALL ERROR$STOP('DYNOCC_INIOCCBYNUMBER')
             END IF
             X(ISVAR+1,:,ISPIN)=SVAR
           END IF
         ENDDO
       END IF
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE DYNOCC_INIOCCBYENERGY(NB,NKPT,NSPIN,FMAX &
      &          ,TEMP,TFIXTOT,TOTCHA,TOTPOT,TFIXSPIN,SPINCHA,SPINPOT &
      &          ,WKPT,EPSILON,X0)
!      *****************************************************************
!      ** CHOOSE INITIAL CONDITIONS FOR THE OCCUPATIONS               **
!      *****************************************************************
       IMPLICIT NONE
       INTEGER(4),INTENT(IN)    :: NB
       INTEGER(4),INTENT(IN)    :: NKPT
       INTEGER(4),INTENT(IN)    :: NSPIN
       REAL(8)   ,INTENT(IN)    :: FMAX    ! MAX(OCCUPATION PER STATE)
       REAL(8)   ,INTENT(IN)    :: TEMP
       LOGICAL(4),INTENT(IN)    :: TFIXTOT
       LOGICAL(4),INTENT(IN)    :: TFIXSPIN
       REAL(8)   ,INTENT(INOUT) :: TOTCHA
       REAL(8)   ,INTENT(INOUT) :: TOTPOT
       REAL(8)   ,INTENT(INOUT) :: SPINCHA
       REAL(8)   ,INTENT(INOUT) :: SPINPOT
       REAL(8)   ,INTENT(IN)    :: EPSILON(NB,NKPT,NSPIN)
       REAL(8)   ,INTENT(IN)    :: WKPT(NKPT)
       REAL(8)   ,INTENT(OUT)   :: X0(NB,NKPT,NSPIN)
       REAL(8)                  :: SVAR
       INTEGER(4)               :: IB,IKPT,ISPIN
       INTEGER(4)               :: ISPINDEG
       REAL(8)                  :: EMERMN
       REAL(8)                  :: SIGMA
       REAL(8)                  :: DSVAR
       REAL(8)                  :: WORK(NB,NKPT)
!      *****************************************************************
                                  CALL TRACE$PUSH('DYNOCC_INIOCCBYENERGY')
       X0(:,:,:)=0.D0
       ISPINDEG=NINT(FMAX)
!
!      =================================================================
!      ==  DETERMINE INITIAL OCCUPATIONS FROM ONE-ELECTRON ENERGIES   == 
!      =================================================================
       IF(TFIXTOT) THEN
         IF(NSPIN.EQ.1) THEN 
           CALL DYNOCC_MERMIN(NB,NKPT,NSPIN &
     &            ,TOTCHA,ISPINDEG,TEMP,WKPT,EPSILON,X0,TOTPOT,EMERMN)
         ELSE IF(NSPIN.EQ.2) THEN ! NOW THE SPIN POLARIZED CASE
           IF(TFIXSPIN) THEN
             DO ISPIN=1,NSPIN
               SIGMA=DBLE(3-2*ISPIN)
               CALL DYNOCC_MERMIN(NB,NKPT,1 &
     &            ,0.5D0*(TOTCHA+SIGMA*SPINCHA),ISPINDEG,TEMP &
     &            ,WKPT,EPSILON(:,:,ISPIN),X0(:,:,ISPIN),TOTPOT,EMERMN)
             ENDDO
           ELSE 
             DO ISPIN=1,NSPIN
               SIGMA=DBLE(3-2*ISPIN)
               WORK(:,:)=EPSILON(:,:,ISPIN)+SPINPOT*SIGMA
               CALL DYNOCC_MERMIN(NB,NKPT,1 &
     &           ,0.5D0*(TOTCHA+SIGMA*SPINCHA),ISPINDEG,TEMP &
     &           ,WKPT,WORK(:,:),X0(:,:,ISPIN),TOTPOT,EMERMN)
             ENDDO
           END IF
         ELSE
           CALL ERROR$MSG('NSPIN MUST BE EITHER 1  OR 2')
           CALL ERROR$STOP('DYNOCC$INIOCC')
         END IF
       ELSE IF(.NOT.TFIXTOT) THEN
         IF(.NOT.TFIXSPIN) THEN
           DO ISPIN=1,NSPIN
             SIGMA=DBLE(3-2*ISPIN)
             DO IKPT=1,NKPT
               DO IB=1,NB
                 SVAR=(EPSILON(IB,IKPT,ISPIN)-(TOTPOT+SIGMA*SPINPOT))/TEMP
                 IF(ABS(SVAR).LT.55.D0) THEN
                   X0(IB,IKPT,ISPIN)=1.D0/(1.D0+EXP(SVAR))
                 ELSE
                   X0(IB,IKPT,ISPIN)=0.5D0*(1.D0+DSIGN(1.D0,SVAR))
                 END IF
               ENDDO
             ENDDO
           ENDDO
           X0(:,:,:)=X0(:,:,:)*2.D0/DBLE(NSPIN)
         ELSE 
           CALL DYNOCC_MERMIN(NB,NKPT,NSPIN &
     &          ,TOTCHA,ISPINDEG,TEMP,WKPT,EPSILON,X0,TOTPOT,EMERMN)
         END IF
       END IF
!
!      =========================================================================
!      ==  CONVERT OCCUPATIONS INTO X-VARIABLES                               ==
!      =========================================================================
       DO ISPIN=1,NSPIN
         DO IKPT=1,NKPT
           DO IB=1,NB
             SVAR=X0(IB,IKPT,ISPIN)/FMAX
             CALL DYNOCC_XOFF(SVAR,X0(IB,IKPT,ISPIN),DSVAR)
           ENDDO
         ENDDO
       ENDDO
                                  CALL TRACE$POP()
       RETURN
       END
!PB031019END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE DYNOCC_MERMIN(NBANDS,NKPT,NSPIN &
     &                 ,TOTCHA,ISPINDEG,TEMP,WKPT,EIG,F,CHMPOT,EMERMN)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE OCCUPATIONS OF THE ELECTRONIC LEVELS         **
!     **  ACCORDING TO THE FERMI DISTRIBUTION;                        **
!     **  AND CALCULATES THE ENERGY -T*S RELATED TO THE ENTROPY OF    **
!     **  THE ELECTRONS.                                              **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1991 ****
      IMPLICIT NONE
      LOGICAL(4) ,PARAMETER  :: TPR=.FALSE.
      INTEGER(4) ,PARAMETER  :: ITERX=1000   ! MAX #(ITERATIONS)
      REAL(8)    ,PARAMETER  :: TOL=1.D-10    ! TOLERANCE IN #(ELECTRONS)
      INTEGER(4) ,INTENT(IN) :: NBANDS       ! #(BANDS)
      INTEGER(4) ,INTENT(IN) :: NKPT         ! #(K-POINTS)    
      INTEGER(4) ,INTENT(IN) :: NSPIN        ! #(SPINS)    
      REAL(8)    ,INTENT(IN) :: TOTCHA       ! TOTAL CHARGE
      INTEGER(4) ,INTENT(IN) :: ISPINDEG     ! SPIN DEGENERACY (1 OR 2)
      REAL(8)    ,INTENT(IN) :: TEMP         ! K_B*TEMPERATURE IN HARTREE
      REAL(8)    ,INTENT(IN) :: WKPT(NKPT)         ! K-POINT WEIGHTS
      REAL(8)    ,INTENT(IN) :: EIG(NBANDS,NKPT,NSPIN) ! EIGENVALUES
      REAL(8)    ,INTENT(OUT):: F(NBANDS,NKPT,NSPIN)   ! OCCUPATIONS
      REAL(8)    ,INTENT(OUT):: CHMPOT               ! CHEMICAL POTENTIAL
      REAL(8)    ,INTENT(OUT):: EMERMN               ! HEAT OF THE ELECTRONS
      INTEGER(4)             :: ISTART
      INTEGER(4)             :: IB,I,ISPIN,IKPT
      REAL(8)                :: X0,DX,Y0,XM,YM
      REAL(8)                :: SVAR
      REAL(8)                :: SUMV
      INTEGER(4)             :: ITER       ! ITERATION COUNT
      REAL(8)                :: DQ         ! DEVIATION IN TOTAL CHARGE
      REAL(8)                :: EV         ! ELECTRON VOLT
      REAL(8)                :: DE
      REAL(8)                :: F1
      INTEGER(4)             :: IBI
      REAL(8)                :: FMAX         ! #(ELECTRONS PER STATE) 1 OR 2
!     ******************************************************************
                           CALL TRACE$PUSH('DYNOCC_MERMIN')
!
      IF(ISPINDEG.NE.1.AND.ISPINDEG.NE.2) THEN
        CALL ERROR$MSG('SPIN-DEGENERACY CAN BE EITHER ONE OR TWO')
        CALL ERROR$STOP('DYNOCC_MERMIN')
      END IF
      FMAX=REAL(ISPINDEG,KIND=8)
!
!     ==================================================================
!     ==  ESTIMATE CHEMICAL POTENTIAL BY AVERAGING THE ONE-PARTICLE   ==
!     ==  ENERGIES OF THE HIGHEST OCCUPIED BAND                       ==
!     ==================================================================
      IB=INT(TOTCHA/FMAX)  
      IF(IB.GE.NBANDS) THEN
        CALL ERROR$MSG('TO FEW BANDS FOR THE NUMBER OF ELECTRONS')
        CALL ERROR$STOP('DYNOCC_MERMIN')
      END IF
      IB=MAX(IB,1)
      I=0
      CHMPOT=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          I=I+1
          CHMPOT=CHMPOT+EIG(IB,IKPT,ISPIN)
        ENDDO
      ENDDO
      CHMPOT=CHMPOT/DBLE(I) 
                           CALL TRACE$PASS('A')
!
!     ==================================================================
!     ==  FIND CHEMICAL POTENTIAL BY BISECTION                        ==
!     ==================================================================
      X0=CHMPOT
      DX=1.D-2
      ISTART=1
      CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
      CHMPOT=X0
      DO ITER=1,ITERX
        SUMV=0.D0
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB=1,NBANDS
              SVAR=(EIG(IB,IKPT,ISPIN)-CHMPOT)/TEMP
              IF(SVAR.GT.+50.D0)SVAR=+50.D0
              IF(SVAR.LT.-50.D0)SVAR=-50.D0
              F(IB,IKPT,ISPIN)=1.D0/(1.D0+EXP(SVAR))*FMAX
              SUMV=SUMV+F(IB,IKPT,ISPIN)*WKPT(IKPT)
            ENDDO
          ENDDO
        ENDDO
        DQ=SUMV-TOTCHA
        IF(ABS(DQ).LT.TOL) GOTO 110
        X0=CHMPOT
        Y0=DQ
        CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
        CHMPOT=X0
      ENDDO
      CALL ERROR$MSG('OCCUPATIONS NOT CONVERGED')
      CALL ERROR$MSG('PROBABLY THE TEMPERATURE IS ZERO')
      CALL ERROR$STOP('DYNOCC_MERMIN')
 110  CONTINUE
                           CALL TRACE$PASS('B')
!
!     ==================================================================
!     ==  CALCULATE HEAT OF THE ELECTRONS                             ==
!     ==================================================================
      EMERMN=0.D0
      DO ISPIN=1,NSPIN
        DO IKPT=1,NKPT
          DO IB=1,NBANDS
            F1=F(IB,IKPT,ISPIN)/FMAX
            IF(F1.NE.0.D0.AND.1.D0-F1.NE.0.D0) THEN
              DE=+TEMP*(F1*LOG(F1)+(1.D0-F1)*LOG(1.D0-F1))
              EMERMN=EMERMN+DE*WKPT(IKPT)*FMAX
            END IF
          ENDDO
        ENDDO
      ENDDO
                           CALL TRACE$PASS('C')
!
!     ==================================================================
!     ==  PRINT FOR CHECK                                             ==
!     ==================================================================
      IF(TPR) THEN
        CALL CONSTANTS('EV',EV)
        WRITE(*,FMT='("#ELECTRONS( IN)=",F10.5' &
     &               //'," CHEMICAL POTENTIAL=",F10.3' &
     &               //'/"# ELECTRONS(OUT)=",F10.5)')TOTCHA,CHMPOT/EV,TOTCHA+DQ
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            WRITE(*,FMT='(5("(",F8.3,";",F5.2,")"))') &
     &         (EIG(IB,IKPT,ISPIN)/EV,F(IB,IKPT,ISPIN),IB=1,NBANDS)
          ENDDO
        ENDDO
      END IF
                         CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
MODULE OCCUPATION_MODULE
USE LINKEDLIST_MODULE 
LOGICAL(4)    :: TINI=.FALSE.
TYPE(LL_TYPE) :: LL_OCC
CONTAINS
   SUBROUTINE OCCUPATION_NEWLIST
   IF(TINI) RETURN
   CALL LINKEDLIST$NEW(LL_OCC)
   TINI=.TRUE.
   RETURN
   END SUBROUTINE OCCUPATION_NEWLIST
END MODULE OCCUPATION_MODULE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE OCCUPATION$SET(STRING_,NBYTE_,VAL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE OCCUPATION_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRING_
      INTEGER(4)  ,INTENT(IN) :: NBYTE_
      REAL(8)     ,INTENT(IN) :: VAL(NBYTE_)
!     ******************************************************************
      CALL TRACE$PUSH('OCCUPATION$SET')
      CALL OCCUPATION_NEWLIST
      CALL LINKEDLIST$SET(LL_OCC,STRING_,0,VAL)
      CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE OCCUPATION$GET(STRING_,NBYTE_,VAL)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE OCCUPATION_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRING_
      INTEGER(4)  ,INTENT(IN) :: NBYTE_
      REAL(8)     ,INTENT(OUT):: VAL(NBYTE_)
      LOGICAL(4)              :: TCHK
!     ******************************************************************
      CALL TRACE$PUSH('OCCUPATION$GET')
      CALL OCCUPATION_NEWLIST
      CALL LINKEDLIST$EXISTD(LL_OCC,STRING_,1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$MSG('ITEM HAS NOT BEEN STORED')
        CALL ERROR$CHVAL('STRING',STRING_)
        CALL ERROR$STOP('OCCUPATION$GET')
      END IF
      CALL LINKEDLIST$GET(LL_OCC,STRING_,1,VAL)
      CALL TRACE$POP
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE OCCUPATION$REPORT(NFIL,STRING_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE DYNOCC_MODULE, ONLY : NKPT,NKDIV,XK
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRING_
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4)              :: IKPT,I,J
      REAL(8)                 :: DWORK(3)
      REAL(8)                 :: CELLVOL
      REAL(8)                 :: RBAS(3,3)
      REAL(8)                 :: GBAS(3,3)
!     ******************************************************************
      IF(STRING_.EQ.'KPOINTS') THEN
        CALL REPORT$TITLE(NFIL,'K-POINTS')
        CALL REPORT$I4VAL(NFIL,'KDIV(1)',NKDIV(1),' ')        
        CALL REPORT$I4VAL(NFIL,'KDIV(2)',NKDIV(2),' ')        
        CALL REPORT$I4VAL(NFIL,'KDIV(3)',NKDIV(3),' ')        
!
        CALL CELL$GETR8A('T(0)',9,RBAS)
        CALL GBASS(RBAS,GBAS,CELLVOL)
        DO IKPT=1,NKPT
          DO I=1,3
            DWORK(I)=0.D0
            DO J=1,3
              DWORK(I)=DWORK(I)+GBAS(I,J)*XK(J,IKPT)
            ENDDO
          ENDDO
          WRITE(NFIL,FMT='(" K",I4' &
     &             //'," = (",F7.4," * G1, ",F7.4," * G2, ",F7.4,"* G3)"' &
     &             //'," = (", F8.5,", ",F8.5,",", F8.5,");")') &
     &             IKPT,(XK(I,IKPT),I=1,3) &
     &             ,(DWORK(I),I=1,3)
        ENDDO
      ELSE 
        CALL ERROR$MSG('STRING_ NOT RECOGNIZED')
        CALL ERROR$CHVAL('STRING_',STRING_)
        CALL ERROR$STOP('OCCUPATION$REPORT')
      END IF
      RETURN
      END
      
