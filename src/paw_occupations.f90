!.......................................................................
MODULE DYNOCC_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: DYNOCC                                                     **
!**                                                                   **
!**  PURPOSE: ORGANIZES THE OCCUPATIONS OF THE ONE-PARTICLE STATES    **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    DYNOCC$CREATE                                                  **
!**    DYNOCC$SET                                                     **
!**    DYNOCC$GET                                                     **
!**    DYNOCC$STOP                                                    **
!**    DYNOCC$INIOCC                                                  **
!**    DYNOCC$MODOCC                                                  **
!**    DYNOCC$PROPAGATE                                               **
!**    DYNOCC$SWITCH                                                  **
!**    DYNOCC$REPORT                                                  **
!**                                                                   **
!**  REMARKS:                                                         **
!**    THE OCCUPATIONS F ARE DIRECTLY RELATED TO THE DYNAMICAL        **
!**    VARIABLES X BY F=3*X**2-2*X**3                                 **
!**                                                                   **
!**  LOGICAL INPUT HIRARCHY:                                          **
!**    NB,NKPT,NSPIN                                                  **
!**    SUMOFZ                                                         **
!**    TDYN                                                           **
!**    IF(TDYN) THEN                                                  **
!**      TEMP,MX,DELAT,ANNEX,TSTOP,FREEZE                             **
!**      TFIXTOT                                                      **
!**      TFIXSPIN                                                     **
!**      IF(TFIXTOT) THEN                                             **
!**        TOTCHA   (ON INPUT AS EXCESS CHARGE)                       **
!**      ELSE                                                         **
!**        TOTPOT                                                     **
!**      END IF                                                       **
!**      IF(TFIXSPIN) THEN                                            **
!**        SPINCHA                                                    **
!**      ELSE                                                         **
!**        SPINPOT                                                    **
!**      END IF                                                       **
!**    ELSE                                                           **
!**      OCC                                                          **
!**    END IF                                                         **
!**    STARTTYPE                                                      **
!**                                                                   **
!**    FOUR DIFFERENT STARTTYPE ARE POSSIBLE:                         **
!**      'X' READS THE OCCUPATION VARIABLES FROM RESTART FILE         **
!**      'E' READS THE ENERGIES FROM RESTART FILE                     **
!**      'N' FILLS STATES ACCORDING TO BAND NUMBER AND DYNOCC$MODOCC  **
!**                                                                   **
!**    THE DEFAULT OCCUPATIONS CORRESPOND TO STARTTYPE='N'.           **
!**    IT FILLS THE BANDS CONSISTENT WITH THE TOTAL CHARGE AND SPIN   **
!**    ASSUMING FLAT BANDS. THEN THEY ARE MODIFIED BY DYNOCC$MODOCC   **
!**    (DYNOCC$MODOCC ALSO ADJUSTS THE TOTAL CHARGE AND SPIN!)        **
!**    FOR STARTTYPE='E' OR 'X' THE RESTART FILE IS READ AND THE      **
!**    OCCUPATIONS ARE OVERWRITTEN. IF STARTTYPE='X' THE OCCUPATIONS  **
!**    ARE READ FROM FILE. FOR STARTTYPE='E' THE ENERGIES ARE READ    **
!**    FROM FILE AND THE OCCUPATIONS ARE CREATED NEW USING THE        **
!**    TEMPERATURE, THE TOTAL CHARGE OR SPIN (OR SPINPOT).            **
!**                                                                   **
!******************************************* PETER E. BLOECHL, 1996 *****
CHARACTER(1):: STARTTYPE       ! CAN BE 'N','E','X'
INTEGER(4)  :: NB=0            ! #(BANDS)
INTEGER(4)  :: NKPT=0          ! #(K-POINTS)
INTEGER(4)  :: NSPIN=0         ! #(SPINS)
LOGICAL(4)  :: TDYN=.FALSE.    ! DYNAMICAL/STATIC OCCUPATION
LOGICAL(4)  :: RESET=.TRUE.    ! SETS OUTPUT MODE OF DYNOCC$REPORT
LOGICAL(4)  :: TSTOP=.FALSE.   ! SET VELOCITY TO ZERO
LOGICAL(4)  :: TFIXSPIN=.FALSE.! FIXED SPIN/ MAGNETIC FIELD 
LOGICAL(4)  :: TFIXTOT=.FALSE. ! FIXED CHARGE/CHEMICAL POTENTIAL
REAL(8)     :: FMAX=1.D0       ! MAX OCCUPATION OF A SINGLE STATE
REAL(8)     :: SUMOFZ=-1.D0    ! SUM OF NUCLEAR CHARGES
REAL(8)     :: CHARGE=0.D0     ! EXCESS NUMBER
REAL(8)     :: TOTCHA=0.D0     ! NUMBER OF ELECTRONS (=SUMOFZ+CHARGE)
REAL(8)     :: SPINCHA=0.D0    ! TOTAL SPIN [ELECTRON SPINS]
REAL(8)     :: MX=800.D0       ! MASS OF OCCUPATION DYNAMICS
REAL(8)     :: TOTPOT=0.D0     ! FERMI LEVEL
REAL(8)     :: SPINPOT=0.D0    ! MAGNETIC FIELD
REAL(8)     :: TEMP=0.D0       ! TEMPERATURE
REAL(8)     :: ANNEX=0.D0      ! FRICTION 
REAL(8)     :: DELTAT=0.D0     ! TIME STEP
CHARACTER(32) :: BZITYPE='MP'  ! TYPE OF BRILLOUIN-ZONE INTEGRATION METHOD
                               ! CAN BE 'MP' FOR MONCKHORST-PACK SAMPLING
                               !     OR 'TETRA+' FOR THE IMPROVED TETRAHEDRON METHOD
REAL(8)   ,ALLOCATABLE :: XM(:,:,:)     ! DYNAMICAL VARIABLES FRO OCCUPATIONS
REAL(8)   ,ALLOCATABLE :: X0(:,:,:)
REAL(8)   ,ALLOCATABLE :: XP(:,:,:)
LOGICAL(4)             :: TEPSILON=.FALSE.
REAL(8)   ,ALLOCATABLE :: EPSILON(:,:,:)  ! DE/DF=<PSI|H|PSI>
!! REMARK MPSIDOT2 IS NOT PROPERLY TREATED, BECAUSE IT IS CALCULATED IN TWO STEPS
REAL(8)   ,ALLOCATABLE :: MPSIDOT2(:,:,:) ! M_PSI<PSIDOT||PSIDOT> 
!====== K-POINT RELATED
REAL(8)   ,ALLOCATABLE :: XK(:,:) !(3,NKPT)
REAL(8)   ,ALLOCATABLE :: WKPT(:) !(NKPT)
CONTAINS
!      .................................................................
       SUBROUTINE CUBPOLYNOMROOT(A0_,A1_,A2_,A3_,X1,X2,X3)
!      **                                                             **
!      **  SEARCHES THE THREE ZEROS OF THE 3RD ORDER POLYNOMIAL       **
!      **          A0_+A1_*X+A2_*X**2+A3_*X**3                        **
!      **  IF ONLY ONE ZERO EXISTS, X1=X2=X3                          **
!      **                                                             **
!      **  SEE "NUMERICAL RECIPES"; CAMBRIDGE UNIVERSITY PRESS        **
!      **     SECTION 5.5 QUADRATIC AND CUBIC EQUATIONS, P 145        **
!      **                                                             **
!      **  ATTENTION NOT ALL ROOOTS AR EVALUATED                      **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: A0_,A1_,A2_,A3_
       REAL(8)   ,INTENT(OUT):: X1,X2,X3
       REAL(8)               :: A0,A1,A2,A3
       REAL(8)               :: Q,R,Q3,Q3MR2,THETA,SVAR
       REAL(8)               :: PI
       LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
!      *****************************************************************
       CALL CUBNEWTONRAPHSON(A0_,A1_,A2_,A3_,X1)
       X2=X1
       X3=X1
       RETURN
!      =================================================================
       PI=4.D0*DATAN(1.D0)
       A3=1.D0
       A2=A2_/A3_
       A1=A1_/A3_
       A0=A0_/A3_
       Q=(A2**2-3.D0*A1)/9.D0
       R=(2.D0*A2**3-9.D0*A2*A1+27.D0*A0)/54.D0
       Q3=Q**3
       Q3MR2=Q3-R**2
       IF(Q3MR2.GE.0) THEN
         THETA=ACOS(R/DSQRT(Q3))
         SVAR=-2.D0*DSQRT(Q)
         X1=SVAR*COS(THETA/3.D0)
         X2=SVAR*COS((THETA+2.D0*PI)/3.D0)
         X3=SVAR*COS((THETA+4.D0*PI)/3.D0)
         IF(X1.GT.X2) THEN
           SVAR=X1
           IF(X1.GT.X3) THEN
             X1=X3 ; X3=SVAR
           ELSE
             X1=X2 ; X2=SVAR
           END IF
         END IF
         IF(X2.GT.X3) THEN
           SVAR=X2 ; X2=X3 ; X3=SVAR
         END IF
       ELSE
         SVAR=(DSQRT(-Q3MR2)+DABS(R))**(1.0/3.0)
         X1=-SIGN(1.D0,R)*(SVAR+Q/SVAR)
         X2=X1
         X3=X1
       END IF
       SVAR=A2/3.D0
       X1=X1-SVAR
       X2=X2-SVAR
       X3=X3-SVAR
       IF(TTEST) THEN
         SVAR=ABS(A0_+X1*(A1_+X1*(A2_+X1*A3_)))
         IF(SVAR.GT.1.D-8) THEN
           CALL ERROR$STOP('CUBPOLYNOMROOT1')
         ENDIF
         SVAR=ABS(A0_+X2*(A1_+X2*(A2_+X2*A3_)))
         IF(SVAR.GT.1.D-8) THEN
           CALL ERROR$STOP('CUBPOLYNOMROOT2')
         ENDIF
         SVAR=ABS(A0_+X3*(A1_+X3*(A2_+X3*A3_)))
         IF(SVAR.GT.1.D-8) THEN
           CALL ERROR$STOP('CUBPOLYNOMROOT3')
         ENDIF
       END IF
       RETURN
       END SUBROUTINE CUBPOLYNOMROOT
!      .................................................................
       SUBROUTINE CUBNEWTONRAPHSON(A0,A1,A2,A3,X1)
!      **                                                             **
!      **  SEARCHES THE THREE ZEROS OF THE 3RD ORDER POLYNOMIAL       **
!      **          A0_+A1_*X+A2_*X**2+A3_*X**3                        **
!      **  IF ONLY ONE ZERO EXISTS, X1=X2=X3                          **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: A0,A1,A2,A3
       REAL(8)   ,INTENT(OUT):: X1
       REAL(8)               :: VAL,DER,SVAR
       INTEGER(4),PARAMETER  :: NITER=1000
       INTEGER(4)            :: ITER
       REAL(8)   ,PARAMETER  :: TOL=1.D-8
!      *****************************************************************
       X1=0.D0
       DO ITER=1,NITER
         VAL=A0+X1*(A1+X1*(A2+X1*A3))
         DER=A1+X1*(2.D0*A2+X1*3.D0*A3)
         SVAR=-VAL/DER
         IF(ABS(VAL).LT.TOL) RETURN
         X1=X1+SVAR
       ENDDO
       CALL ERROR$MSG('LOOP NOT CONVERGED')
       CALL ERROR$STOP('CUBNEWTONRAPHSON(DYNOCC_MODULE)')
       RETURN
       END SUBROUTINE CUBNEWTONRAPHSON
END MODULE DYNOCC_MODULE
!
!      .................................................................
       SUBROUTINE DYNOCC_FOFX(X,F,DF)
!      **                                                             **
!      **  EVALUATES THE OCCUPATIONS FROM THE DYNAMICAL VARIABLES     **
!      **  OCC(X)=3X^2-2X^3  A CUBE POLYNOM                          **
!      **  WITH A MIN AT (0,0) AND A MAX AT (1,1)                    **
!      **  OCCOFX IS MADE PERIODIC SO THAT THE OCCUPATIONS ARE       **
!      **  STRICTLY BETWEEN ZERO AND ONE                             **
!      **                                                             **
!      *****************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: X
       REAL(8),INTENT(OUT):: F
       REAL(8),INTENT(OUT):: DF   !DF/DX
       LOGICAL(4),SAVE    :: TNEW=.TRUE.
       REAL(8)   ,SAVE    :: XMIN,XMAX
       REAL(8)            :: X1
       REAL(8)            :: DI
!      *****************************************************************
       IF(TNEW) THEN
         XMIN=0.5D0-DSQRT(0.75D0)
         XMAX=0.5D0+DSQRT(0.75D0)
         TNEW=.FALSE.
       END IF
       DI=(X-XMIN)/(XMAX-XMIN)
       DI=DI-INT(DI)
       IF(DI.LT.0.D0) DI=DI+1.D0
       X1=XMIN+DI*(XMAX-XMIN)
       F=X1*X1*(3.D0-2.D0*X1)
       DF=6.D0*X1*(1.D0-X1)
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
       REAL(8)             :: F
       REAL(8)             :: DF
!      *****************************************************************
       CALL DYNOCC_FOFX(X,F,DF)
       IF(F.EQ.0.D0.OR.F.EQ.1.D0) THEN
         S=0.D0
         DS=0.D0
         RETURN
       END IF
       S=F*LOG(F)+(1.D0-F)*LOG(1.D0-F)
       DS=(LOG(F)-LOG(1.D0-F))*DF
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
       REAL(8)   ,INTENT(OUT):: X
       REAL(8)   ,INTENT(OUT):: DX    
       INTEGER(4),PARAMETER  :: NITER=100
       REAL(8)   ,SAVE       :: TOL=1.D-12
       REAL(8)   ,PARAMETER  :: DFMIN=1.D-3
       INTEGER(4)            :: ITER
       REAL(8)               :: FI,DF
!      *****************************************************************
       IF(F.LT.0.D0.OR.F.GT.1.D0) THEN
         CALL ERROR$MSG('OCCUPATIONS MUST LIE BETWEEN ZERO AND ONE')
         CALL ERROR$STOP('DYNOCC_XOFF')
       END IF
       X=0.5D0
       DO ITER=1,100
         CALL DYNOCC_FOFX(X,FI,DF)
         IF(ABS(DF).LT.DFMIN) THEN
           DF=DFMIN
           IF(DF.LT.0.D0) DF=-DFMIN
         END IF  
         DX=(F-FI)/DF
         IF(ABS(DX).LT.TOL) EXIT
         X=X+DX
       ENDDO
       DX=1/DF
       RETURN
       END 
!      .................................................................
       SUBROUTINE DYNOCC$TEST
       INTEGER(4),PARAMETER :: NB=5
       INTEGER(4),PARAMETER :: NKPT=2
       INTEGER(4),PARAMETER :: NSPIN=2
       REAL(8)              :: WKPT(NKPT)
       REAL(8)              :: EPSILON(NB,NKPT,NSPIN)
       REAL(8)              :: WGHT(NB,NKPT,NSPIN)
       REAL(8)   ,PARAMETER :: SUMOFZ=5.D0
       REAL(8)              :: FMAX=1.D0
       REAL(8)              :: MASS=1.D+4
       REAL(8)              :: DELTAT=10.
       REAL(8)              :: FRIC=1.D-3
       REAL(8)              :: TEMP=1.D-3
       INTEGER(4),PARAMETER :: NITER=1000
       INTEGER(4)           :: IB,IKPT,ISPIN,ITER
       INTEGER(4)           :: NFILO
       REAL(8)              :: SVAR
       REAL(8)              :: EKIN
       REAL(8)              :: EPOT
       REAL(8)              :: EBAND
       REAL(8)              :: XK(3,NKPT)
!      ******************************************************************
       CALL FILEHANDLER$UNIT('PROT',NFILO)
       FMAX=2.D0/REAL(NSPIN,KIND=8)
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
         ENDDO
       ENDDO
       CALL DYNOCC$SETR8A('XK',3*NKPT,XK)
       CALL DYNOCC$SETR8A('WKPT',NKPT,WKPT)
       CALL DYNOCC$SETR8('SUMOFZ',SUMOFZ)
       CALL DYNOCC$SETR8('MASS',MASS)
       CALL DYNOCC$SETR8('TIMESTEP',DELTAT)
       CALL DYNOCC$SETR8('TEMP',1.D0)
       CALL DYNOCC$SETL4('FIXQ',.TRUE.)
       CALL DYNOCC$SETL4('FIXS',.TRUE.)
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
         CALL DYNOCC$PROPAGATE
         CALL DYNOCC$GETR8('EKIN',EKIN)
         CALL DYNOCC$GETR8('EPOT',EPOT)
         CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,WGHT)
         EBAND=0.D0
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             SVAR=0.D0
             DO IB=1,NB
               SVAR=SVAR+WGHT(IB,IKPT,ISPIN)*EPSILON(IB,IKPT,ISPIN)
             ENDDO
             EBAND=EBAND+SVAR
           ENDDO
         ENDDO
         WRITE(NFILO,FMT='(I10,4E15.5)')ITER,EKIN+EPOT+EBAND,EKIN,EPOT,EBAND
         CALL DYNOCC$SWITCH
       ENDDO
       CALL DYNOCC$REPORT(NFILO)
       RETURN
       END

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
       XM(:,:,:)=0.5D0
       X0(:,:,:)=0.5D0
       XP(:,:,:)=0.5D0
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
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$SETL4')
       END IF
       RESET=.TRUE.
       RETURN
       END
!      .................................................................
       SUBROUTINE DYNOCC$SETI4A(ID_,LEN_,VAL)
!      *****************************************************************
!      **  SET LOGICAL PARAMETERS                                     **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID_
       INTEGER(4)  ,INTENT(IN) :: LEN_
       INTEGER(4)  ,INTENT(IN) :: VAL(LEN_)
!      *****************************************************************
       IF(ID_.EQ.'+-+-+-+-') THEN
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID_',ID_)
         CALL ERROR$STOP('DYNOCC$SETI4A')
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
!      .................................................................
       SUBROUTINE DYNOCC$SETR8A(ID_,LEN_,VAL)
!      *****************************************************************
!      **  SET REAL(8) ARRAY                                          **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID_
       INTEGER(4)  ,INTENT(IN) :: LEN_
       REAL(8)     ,INTENT(IN) :: VAL(LEN_)
       REAL(8)                 :: SVAR
!      *****************************************************************
       IF(ID_.EQ.'XK') THEN
         IF(NKPT.EQ.0) THEN
           NKPT=LEN_/3
         END IF
         IF(LEN_.NE.3*NKPT) THEN
           CALL ERROR$MSG('SIZE INCONSISTENT')
           CALL ERROR$STOP('DYNOCC$SETR8A')
         END IF
         IF(.NOT.ALLOCATED(XK))ALLOCATE(XK(3,NKPT))
         XK(:,:)=RESHAPE(VAL,(/3,NKPT/))
       ELSE IF(ID_.EQ.'WKPT') THEN
         IF(NKPT.EQ.0) THEN
           NKPT=LEN_
         END IF
         IF(LEN_.NE.NKPT) THEN
           CALL ERROR$MSG('SIZE INCONSISTENT')
           CALL ERROR$STOP('DYNOCC$SETR8A')
         END IF
         IF(.NOT.ALLOCATED(WKPT))ALLOCATE(WKPT(NKPT))
         WKPT(:)=VAL(:)
!PB031019START
         IF(ABS(SUM(WKPT)-1.D0).GT.1.D-6) THEN
           CALL ERROR$MSG('OCCUPATIONS DO NOT SUM UP TO ONE!')
           CALL ERROR$STOP('DYNOCC$SETR8A')
         END IF
!PB031019END
!      ===================================================================
       ELSE IF(ID_.EQ.'EPSILON') THEN
         IF(LEN_.NE.NB*NKPT*NSPIN) THEN
           CALL ERROR$MSG('DIMENSIONS INCONSISTENT')
           CALL ERROR$CHVAL('ID_',ID_)
           CALL ERROR$I4VAL('LEN_',LEN_)
           CALL ERROR$I4VAL('NB*NKPT*NSPIN',NB*NKPT*NSPIN)
           CALL ERROR$STOP('DYNOCC$SETR8A')
         END IF
         EPSILON=RESHAPE(VAL,(/NB,NKPT,NSPIN/))
         TEPSILON=.TRUE.
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
         MPSIDOT2=MPSIDOT2+SVAR*RESHAPE(VAL,(/NB,NKPT,NSPIN/))
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
       ELSE IF(ID_.EQ.'BZINTEGRATION') THEN
         IF(VAL.NE.'MP'.AND.VAL.NE.'TETRA+') THEN
           CALL ERROR$MSG('ILLEGAL VALUE OF BZINTEGRATION')
           CALL ERROR$MSG('ALLOWED VALUES ARE "MP", "TETRA+"')
           CALL ERROR$CHVAL('ID_',ID_)
           CALL ERROR$STOP('DYNOCC$SETCH')
         END IF
         BZITYPE=VAL
IF(BZITYPE.EQ.'TETRA+') THEN
 CALL BRILLOUIN$TESTING
 CALL ERROR$MSG('FORCED STOP AFTER TESTING BRILLOUIN')
 CALL ERROR$MSG('OPTION NOT IMPLEMENTED')
 CALL ERROR$STOP('DYNOCC$SETCH')
END IF
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
       INTEGER(4)              :: IB,ISPIN,IKPT,IND
       REAL(8)                 :: SVAR
!      *****************************************************************
       IF(ID_.EQ.'SUMOFZ') THEN
         SUMOFZ=VAL
         TOTCHA=SUMOFZ+CHARGE
         RESET=.TRUE.
       ELSE IF(ID_.EQ.'TOTCHA') THEN
         IF(SUMOFZ.LT.0.D0) THEN
           CALL ERROR$MSG('SUMOFZ MUST BE SET BEFORE TOTCHA')
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
!      .................................................................
       SUBROUTINE DYNOCC$GETI4A(ID_,LEN_,VAL)
!      *****************************************************************
!      **  SET INTEGER PARAMETERS                                     **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID_
       INTEGER(4)  ,INTENT(IN)::  LEN_
       INTEGER(4)  ,INTENT(OUT):: VAL(LEN_)
!      *****************************************************************
       IF(ID_.EQ.'') THEN
         VAL(:)=0
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
       REAL(8)                 :: SPEED,SVAR,DSVAR,SUM
       INTEGER(4)              :: IB,IKPT,ISPIN,IND
!      *****************************************************************
       IF(ID_.EQ.'SPIN') THEN
         VAL=SPINCHA    
       ELSE IF(ID_.EQ.'TOTCHA') THEN
         VAL=TOTCHA-SUMOFZ
       ELSE IF(ID_.EQ.'FMAX') THEN
         VAL=FMAX
       ELSE IF(ID_.EQ.'EKIN') THEN
         SVAR=0.D0
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               SPEED=(XP(IB,IKPT,ISPIN)-XM(IB,IKPT,ISPIN))/(2.D0*DELTAT)
               SVAR=SVAR+SPEED**2*WKPT(IKPT)
             ENDDO
           ENDDO
         ENDDO
         VAL=0.5D0*MX*SVAR*FMAX
         IF(.NOT.TDYN) VAL=0.D0
       ELSE IF(ID_.EQ.'HEAT'.OR.ID_.EQ.'EPOT') THEN
!        == REMARK: THE KEYWORD 'HEAT' SHOULD NOT BE USED ANY MORE
!        == ENERGY OF RESERVOIRS: -T*S-MU*N
         SUM=0.D0
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               CALL DYNOCC_SOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
               SUM=SUM+SVAR*WKPT(IKPT)
             ENDDO
           ENDDO
         ENDDO
         VAL=TEMP*SUM*FMAX-TOTPOT*(TOTCHA-SUMOFZ)
!        == THE CHEMICAL POTENTIAL FOR SPIN DIRECTION ISPIN IS
!        == TOTPOT-DBLE(3-2*ISPIN)*SPINPOT
         DO ISPIN=1,NSPIN
           SUM=0.D0
           DO IKPT=1,NKPT
             DO IB=1,NB
               CALL DYNOCC_FOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
               SUM=SUM+SVAR*FMAX*WKPT(IKPT)
             ENDDO
           ENDDO        
           IF(.NOT.TFIXTOT) THEN
             VAL=VAL-TOTPOT*(SUM-SUMOFZ)
           END IF
           IF(.NOT.TFIXSPIN) THEN
             VAL=VAL-DBLE(3-2*ISPIN)*SPINPOT*(SUM-SUMOFZ)
           END IF
         ENDDO
       ELSE IF(ID_.EQ.'TEMP') THEN
         VAL=TEMP
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
       REAL(8)                 :: SPEED,SVAR,DSVAR,SUM
       INTEGER(4)              :: IB,IKPT,ISPIN,IND,I
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
         VAL=RESHAPE(XK,(/3*NKPT/))
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
         IND=0
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             DO IB=1,NB
               IND=IND+1
               CALL DYNOCC_FOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
               VAL(IND)=SVAR*FMAX*WKPT(IKPT)
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
!     ..................................................................
      SUBROUTINE DYNOCC$WRITE(NFIL,NFILO,TCHK)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE DYNOCC_MODULE
      USE RESTART_INTERFACE
      IMPLICIT NONE
      INTEGER(4)       ,INTENT(IN) :: NFIL       ! RESTART-FILE UNIT
      INTEGER(4)       ,INTENT(IN) :: NFILO
      LOGICAL(4)       ,INTENT(OUT):: TCHK
      TYPE(SEPARATOR_TYPE),PARAMETER :: MYSEPARATOR &
          =SEPARATOR_TYPE(4,'OCCUPATIONS','NONE','OCT2003',' ')
      TYPE(SEPARATOR_TYPE),PARAMETER :: OLDSEPARATOR &
          =SEPARATOR_TYPE(3,'OCCUPATIONS','NONE','AUG1996',' ')
      INTEGER(4)                   :: NTASKS,THISTASK
      LOGICAL(4)                   :: TNEW=.TRUE.
!     ******************************************************************
                          CALL TRACE$PUSH('OCCUPATIONS$WRITE')
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      IF(TNEW) THEN
        IF(THISTASK.EQ.1) THEN
          CALL RESTART$WRITESEPARATOR(MYSEPARATOR,NFIL,NFILO,TCHK)
          WRITE(NFIL)NB,NKPT,NSPIN
          WRITE(NFIL)X0(:,:,:)
          WRITE(NFIL)XM(:,:,:)
          WRITE(NFIL)EPSILON(:,:,:)
        END IF
      ELSE
        IF(THISTASK.EQ.1) THEN
          CALL RESTART$WRITESEPARATOR(OLDSEPARATOR,NFIL,NFILO,TCHK)
          WRITE(NFIL)NB,NKPT,NSPIN
          WRITE(NFIL)X0(:,:,:)
          WRITE(NFIL)XM(:,:,:)
        END IF
      END IF
                           CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE DYNOCC$READ(NFIL,NFILO,TCHK)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE DYNOCC_MODULE
      USE RESTART_INTERFACE
      USE MPE_MODULE
      IMPLICIT NONE
      INTEGER(4)         ,INTENT(IN) :: NFIL       ! RESTART-FILE UNIT
      INTEGER(4)         ,INTENT(IN) :: NFILO
      LOGICAL(4)         ,INTENT(OUT):: TCHK
      TYPE(SEPARATOR_TYPE),PARAMETER :: MYSEPARATOR &
          =SEPARATOR_TYPE(4,'OCCUPATIONS','NONE','OCT2003',' ')
      TYPE(SEPARATOR_TYPE),PARAMETER :: OLDSEPARATOR &
          =SEPARATOR_TYPE(3,'OCCUPATIONS','NONE','AUG1996',' ')
      TYPE(SEPARATOR_TYPE)           :: SEPARATOR 
      REAL(8)            ,ALLOCATABLE:: TMP0(:,:,:)
      REAL(8)            ,ALLOCATABLE:: TMPM(:,:,:)
      REAL(8)            ,ALLOCATABLE:: TMPE(:,:,:)
      INTEGER(4)                     :: NB1,NKPT1,NSPIN1
      INTEGER(4)                     :: NTASKS,THISTASK
      INTEGER(4)                     :: ISPIN,IKPT,IB
      REAL(8)                        :: SVAR,DSVAR
      LOGICAL(4)                     :: TOLD
      INTEGER(4)                     :: NFIL1
!     ******************************************************************
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
      CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
      SEPARATOR=MYSEPARATOR
      IF(THISTASK.EQ.1)CALL RESTART$READSEPARATOR(SEPARATOR,NFIL,NFILO,TCHK)
      CALL MPE$BROADCAST('MONOMER',1,TCHK)
      IF(.NOT.TCHK) THEN
        CALL TRACE$POP ;RETURN
      END IF
!     == CHECK IF RESTART FILE FORMAT FROM AUGUST 1996 =================
      TOLD=(SEPARATOR%VERSION.EQ.'AUG1996')
      IF ((STARTTYPE.EQ.'E').AND.(TOLD)) THEN
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
!     ==================================================================
!     == READ DATA                                                    ==
!     ==================================================================
      X0(:,:,:)=0.D0
      XM(:,:,:)=0.D0
      IF(THISTASK.EQ.1) THEN
        IF((SEPARATOR%VERSION.NE.MYSEPARATOR%VERSION).AND.(.NOT.TOLD)) THEN
          CALL ERROR$MSG('VERSION NOT RECOGNIZED')
          CALL ERROR$CHVAL('VERSION',SEPARATOR%VERSION)
          CALL ERROR$STOP('DYNOCC$READ')
        END IF
        READ(NFIL)NB1,NKPT1,NSPIN1
        ALLOCATE(TMP0(NB1,NKPT1,NSPIN1))
        ALLOCATE(TMPM(NB1,NKPT1,NSPIN1))
        IF(.NOT.TOLD)ALLOCATE(TMPE(NB1,NKPT1,NSPIN1))
        READ(NFIL)TMP0(:,:,:)
        READ(NFIL)TMPM(:,:,:)
        IF(.NOT.TOLD)READ(NFIL)TMPE(:,:,:)
!
!       ================================================================
!       == DISCARD SUPERFLOUS DATA AND MAP INTO ARRAY                 ==
!       ================================================================
        NB1=MIN(NB1,NB)
        NKPT1=MIN(NKPT1,NKPT)
        NSPIN1=MIN(NSPIN1,NSPIN)
        X0(1:NB1,1:NKPT1,1:NSPIN1)=TMP0(1:NB1,1:NKPT1,1:NSPIN1)
        XM(1:NB1,1:NKPT1,1:NSPIN1)=TMPM(1:NB1,1:NKPT1,1:NSPIN1)
        IF(.NOT.TOLD)EPSILON(1:NB1,1:NKPT1,1:NSPIN1)=TMPE(1:NB1,1:NKPT1,1:NSPIN1)
        DEALLOCATE(TMP0)
        DEALLOCATE(TMPM)
        IF(.NOT.TOLD)DEALLOCATE(TMPE)
!
!       ================================================================
!       == AUGMENT MISSING DATA                                       ==
!       ================================================================
        DO IKPT=NKPT1+1,NKPT
          X0(:,IKPT,:)=X0(:,1,:)
          XM(:,IKPT,:)=XM(:,1,:)
          IF(.NOT.TOLD)EPSILON(:,IKPT,:)=EPSILON(:,1,:)
        ENDDO
        IF(NSPIN1.EQ.1.AND.NSPIN.EQ.2) THEN
          X0(:,:,2)=X0(:,:,1)
          XM(:,:,2)=XM(:,:,1)
          IF(.NOT.TOLD)EPSILON(:,:,2)=EPSILON(:,:,1)
        END IF
      END IF
!
!     ==================================================================
!     ==  BROADCAST                                                   ==
!     ==================================================================
      CALL MPE$BROADCAST('MONOMER',1,X0)
      CALL MPE$BROADCAST('MONOMER',1,XM)
      IF(.NOT.TOLD)CALL MPE$BROADCAST('MONOMER',1,EPSILON)
!
!     ==================================================================
!     == CONVERT ENERGIES INTO OCCUPATIONS                            ==
!     ==================================================================
      IF(STARTTYPE.EQ.'E'.AND.(.NOT.TOLD)) THEN
        CALL DYNOCC_INIOCCBYENERGY(NB,NKPT,NSPIN,FMAX &
      &          ,TEMP,TFIXTOT,TOTCHA,TOTPOT,TFIXSPIN,SPINCHA,SPINPOT &
      &          ,WKPT,EPSILON,X0)
        XM(:,:,:)=X0(:,:,:)
      END IF
!
!     ==================================================================
!     ==  RESET THERMODYNAMIC VARIABLES                               ==
!     ==================================================================
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
!      .................................................................
       SUBROUTINE DYNOCC$REPORT(NFIL)
!      *****************************************************************
!      **                                                             **
!      **  THE DATA OF THE DYNOCC OBJECT ARE DIVIDED IN STATIC        **
!      **  DYNAMIC DATA.                                              **
!      **  THE PARAMETER RESET IS CHANGED TO TRUE WHENEVER STATIC     **
!      **  DATA ARE CHANGED AND TO FALSE AFTER EACH REPORT.           **
!      **  STATIC DATA ARE THEREFORE REPORTED, WHENEVER ONE OF THEM   **
!      **  HAS BEEN MODIFIED.                                         **
!      **                                                             **
!      **  DYNAMIC DATA ARE ALWAYS REPORTED                           **
!      **                                                             **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       INTEGER(4),INTENT(IN) :: NFIL
       REAL(8)               :: SVAR
       REAL(8)               :: OCC(NB,NKPT,NSPIN)
       INTEGER(4)            :: ISPIN,IKPT
       REAL(8)               :: EV
       REAL(8)               :: KELVIN
!      *****************************************************************
       CALL CONSTANTS('EV',EV)
       CALL CONSTANTS('KB',KELVIN)
!      == 1EL<->0.5HBAR => HBAR=2 EL
       IF(.NOT.RESET.AND.TDYN) THEN
         CALL REPORT$TITLE(NFIL,'OCCUPATIONS')
         CALL REPORT$R8VAL(NFIL,'TEMPERATURE',TEMP/KELVIN,'K')
         CALL REPORT$R8VAL(NFIL,'CHEMICAL POTENTIAL',TOTPOT/EV,'EV')
         CALL REPORT$R8VAL(NFIL,'CHARGE',-(TOTCHA-SUMOFZ),'E')
!
         IF(NSPIN.EQ.2) THEN
           CALL REPORT$R8VAL(NFIL,'MAGNETIC FIELD',SPINPOT/(EV),'EV/(HBAR/2)')
           CALL REPORT$R8VAL(NFIL,'SPIN',SPINCHA/2.D0,'HBAR')
         END IF
         CALL DYNOCC$GETR8('HEAT',SVAR)
         CALL REPORT$R8VAL(NFIL,'HEAT',SVAR,'H')
         IF(TEPSILON) THEN
           CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
           WRITE(NFIL,*)'OCCUPATIONS AND ENERGY EXPECTATION VALUES [EV]'
           DO ISPIN=1,NSPIN
             DO IKPT=1,NKPT
               WRITE(NFIL,*)'FOR K-POINT:',IKPT,' AND SPIN:',ISPIN
               CALL DYNOCC_PREIG('EIG',NFIL,NB,EPSILON(:,IKPT,ISPIN),EV)
               CALL DYNOCC_PREIG('OCC',NFIL,NB,OCC(:,IKPT,ISPIN),WKPT(IKPT))
             ENDDO
           ENDDO
         END IF
       END IF
!
!      =================================================================
!      ==                                                             ==
!      =================================================================
       IF(RESET.AND.(.NOT.TDYN)) THEN
         CALL REPORT$TITLE(NFIL,'OCCUPATIONS')
         CALL REPORT$CHVAL(NFIL,'OCCUPATIONS ARE','FIXED')
         CALL REPORT$R8VAL(NFIL,'#(ELECTRONS)',TOTCHA,' ')
         CALL REPORT$R8VAL(NFIL,'CHARGE',-(TOTCHA-SUMOFZ),'E')
         IF(NSPIN.EQ.2) THEN
           CALL REPORT$R8VAL(NFIL,'SPIN',SPINCHA/2.D0,'HBAR')
         END IF
         CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
         DO ISPIN=1,NSPIN
           DO IKPT=1,NKPT
             WRITE(NFIL,*)'OCCUPATIONS FOR K-POINT:',IKPT &
      &                                           ,' AND SPIN:',ISPIN
             CALL DYNOCC_PREIG('OCC',NFIL,NB,OCC(:,IKPT,ISPIN),WKPT(IKPT))
           ENDDO
         ENDDO
       END IF
!
!      =================================================================
!      ==                                                             ==
!      =================================================================
       IF(RESET.AND.TDYN) THEN
         CALL REPORT$TITLE(NFIL,'OCCUPATIONS')
         CALL REPORT$STRING(NFIL,'DYNAMICAL OCCUPATIONS USING MERMIN FUNCTIONAL')
         CALL REPORT$R8VAL(NFIL,'TEMPERATURE',TEMP/KELVIN,'K')
         CALL REPORT$R8VAL(NFIL,'MASS OF X',MX,'A.U.')
         IF(ANNEX.NE.0.D0) THEN
           CALL REPORT$R8VAL(NFIL,'FRICTION',ANNEX,' ')
         END IF
         IF(TSTOP) THEN
           CALL REPORT$CHVAL(NFIL,'INITIAL VELOCITIES ARE SET TO ','ZERO')
         END IF
         IF(TFIXTOT) THEN
           CALL REPORT$R8VAL(NFIL,'FIXED CHARGE',-(TOTCHA-SUMOFZ),'E')
           CALL REPORT$R8VAL(NFIL,'#(ELECTRONS)',TOTCHA,' ')
           IF(TEPSILON) THEN
             CALL REPORT$R8VAL(NFIL,'CHEMICAL POTENTIAL',TOTPOT/EV,'EV')
           END IF
         ELSE
           CALL REPORT$R8VAL(NFIL,'SUM OF NUCLEAR CHARGES',SUMOFZ,'E')
           CALL REPORT$R8VAL(NFIL,'FIXED CHEMICAL POTENTIAL',TOTPOT/EV,'EV')
           IF(TEPSILON) THEN
             CALL REPORT$R8VAL(NFIL,'CHARGE',-(TOTCHA-SUMOFZ),'E')
           END IF
         END IF
!
         IF(NSPIN.EQ.2) THEN
           IF(TFIXSPIN) THEN
             CALL REPORT$R8VAL(NFIL,'FIXED SPIN',SPINCHA*0.5D0,'HBAR')
             IF(TEPSILON) THEN
               CALL REPORT$R8VAL(NFIL,'MAGNETIC FIELD',SPINPOT/EV,'EV/(HBAR/2)')
             END IF
           ELSE
             CALL REPORT$R8VAL(NFIL,'FIXED MAGNETIC FIELD',SPINPOT/EV,'EV/(HBAR/2)')
             IF(TEPSILON) THEN
               CALL REPORT$R8VAL(NFIL,'SPIN',SPINCHA/2.D0,'HBAR')
             END IF
           END IF
         END IF
         CALL DYNOCC$GETR8('HEAT',SVAR)
         CALL REPORT$R8VAL(NFIL,'HEAT',SVAR,'H')
         IF(TEPSILON) THEN
           CALL DYNOCC$GETR8A('OCC',NB*NKPT*NSPIN,OCC)
           WRITE(NFIL,*)'OCCUPATIONS AND ENERGY EXPECTATION VALUES [EV]'
           DO ISPIN=1,NSPIN
             DO IKPT=1,NKPT
               WRITE(NFIL,*)'FOR K-POINT:',IKPT,' AND SPIN:',ISPIN
               CALL DYNOCC_PREIG('EIG',NFIL,NB,EPSILON(:,IKPT,ISPIN),EV)
               CALL DYNOCC_PREIG('OCC',NFIL,NB,OCC(:,IKPT,ISPIN),WKPT(IKPT))
             ENDDO
           ENDDO
         END IF
       END IF
       IF((.NOT.RESET).AND.(.NOT.TDYN)) THEN
         IF(TEPSILON) THEN
!          CALL DYNOCC$GETR8A('EPSILON',NB*NKPT*NSPIN,EPSILON)
           WRITE(NFIL,*)'ENERGY EXPECTATION VALUES [EV]'
           DO ISPIN=1,NSPIN
             DO IKPT=1,NKPT
               WRITE(NFIL,*)'FOR K-POINT:',IKPT,' AND SPIN:',ISPIN 
               CALL DYNOCC_PREIG('EIG',NFIL,NB,EPSILON(:,IKPT,ISPIN),EV)
             ENDDO
           ENDDO
         END IF
       END IF
       RESET=.FALSE.
       RETURN
       END
!
!      .................................................................
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
!      .................................................................
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
!      .................................................................
       SUBROUTINE DYNOCC$PROPAGATE()
!      *****************************************************************
!      **                                                             **
!      *****************************************************************
       USE DYNOCC_MODULE
       IMPLICIT NONE
       REAL(8)   ,PARAMETER :: TOL=1.D-10
       REAL(8)   ,PARAMETER :: DSMALL=1.D-5
       REAL(8)              :: FX(NB,NKPT,NSPIN)
       REAL(8)              :: XBAR(NB,NKPT,NSPIN)
       INTEGER(4)           :: ISPIN,IKPT,IB
       REAL(8)              :: SVAR,OCC,DOCCDX,FORCE,DSVAR
       REAL(8)              :: SVAR1,SVAR2,SVAR3
       REAL(8)              :: X1,X2,X3
       REAL(8)              :: XMAX
       REAL(8)              :: XMIN
       REAL(8)              :: SIGMA    ! SPIN DIRECTION, I.E.: (+/-)1.D0
       REAL(8)              :: DTOTPOT,DSPINPOT
       REAL(8)              :: S0,S1,S2,S3
       REAL(8)              :: Q0,Q1,Q2,Q3
       LOGICAL(4)           :: TMPSIDOT2
       INTEGER(4)           :: ISVAR
REAL(8)::EV
       LOGICAL(4)           :: TCONV
       REAL(8)              :: QTOL=1.D-10
       INTEGER(4),PARAMETER :: NITER=1000
       INTEGER(4)           :: ITER
INTEGER(4) :: NTASKS,THISTASK
INTEGER(4) :: COUNT = 0
!      *****************************************************************
!!!!!!!
!!DEBUG CLEMENS
!!!!!!!
!OPEN(10,FILE='DYNOCC_STATUS.DAT',STATUS='UNKNOWN',FORM='FORMATTED')
!WRITE(10,*) STARTTYPE! ,' STARTTYPE'
!WRITE(10,*) NB!,' NB'
!WRITE(10,*) NKPT!,' NKPT'
!WRITE(10,*) NSPIN!,' NSPIN'
!WRITE(10,*) TDYN!,' TDYN'
!WRITE(10,*) RESET!,' RESET'
!WRITE(10,*) TSTOP!,' TSTOP'
!WRITE(10,*) TFIXSPIN!,' TFIXSPIN'
!WRITE(10,*) TFIXTOT!,' TFIXTOT'
!WRITE(10,*) FMAX!,' FMAX' 
!WRITE(10,*) SUMOFZ!,' SUMOFZ'
!WRITE(10,*) CHARGE!,' CHARGE'
!WRITE(10,*) TOTCHA!,' TOTCHA'
!WRITE(10,*) SPINCHA!,' SPINCHA'
!WRITE(10,*) MX!,' MX'
!WRITE(10,*) TOTPOT!,' TOTPOT'
!WRITE(10,*) SPINPOT!,' SPINPOT'
!WRITE(10,*) TEMP!,' TEMP'
!WRITE(10,*) ANNEX!,' ANNEX'
!WRITE(10,*) DELTAT!,' DELTAT'
!WRITE(10,*) TEPSILON!,' TEPSILON'
!DO ISPIN = 1, NSPIN
!  DO IKPT = 1, NKPT
!    WRITE(10,*) XM(:,IKPT,ISPIN)!,' XM FOR IKPT/ISPIN',IKPT,ISPIN
!    WRITE(10,*) X0(:,IKPT,ISPIN)!,' XM FOR IKPT/ISPIN',IKPT,ISPIN
!    WRITE(10,*) XP(:,IKPT,ISPIN)!,' XM FOR IKPT/ISPIN',IKPT,ISPIN
!    WRITE(10,*) EPSILON(:,IKPT,ISPIN)!,' EPSILON FOR IKPT/ISPIN',IKPT,ISPIN
!    WRITE(10,*) MPSIDOT2(:,IKPT,ISPIN)!,' MPSIDOT2 FOR IKPT,ISPIN',IKPT,ISPIN
!    
!  END DO
!END DO
!DO IKPT = 1, NKPT
!  WRITE(10,*) XK(:,IKPT)!,' XK FOR IKPT',IKPT
!  WRITE(10,*) WKPT(IKPT)!,' WKPT FOR IKPT',IKPT
!END DO
!CLOSE(10)
!!!!!!!
!DEBUG
!!!!!!!
       IF(.NOT.TDYN) THEN 
         XP(:,:,:)=X0(:,:,:)
         RETURN
       END IF
                              CALL TRACE$PUSH('DYNOCC$PROPAGATE')
       IF(NSPIN.EQ.1) THEN
         TFIXSPIN=.FALSE.
         SPINPOT=0.D0
       END IF
       IF(.NOT.TEPSILON) THEN
         CALL ERROR$MSG('ONE-PARTICLE ENERGIES HAVE NOT BEEN SET')
         CALL ERROR$STOP('DYNOCC$PROPAGATE')
       END IF
       TMPSIDOT2=ALLOCATED(MPSIDOT2)
!      IF(TMPSIDOT2) PRINT*,'MPSIDOT2',MPSIDOT2
!      IF(TMPSIDOT2) PRINT*,'E+MPSIDOT2',EPSILON+MPSIDOT2
!
!      =================================================================
!      ==  AVOID ESCAPING OCCUPATIONS BY IMPOSING PERIODIC            ==
!      ==  BOUNDARY CONDISTIONS AT XMIN AND XMAX                      ==
!      ==  EVERY PASS THROUGH THE BOUNDARY X IS STOPPED               ==
!      ==  F(X)=3X^2-2X^3=0.5 FOR X=0.5+-SQRT(3/4)                    ==
!      =================================================================
       XMIN=0.5D0-DSQRT(0.75D0)
       XMAX=0.5D0+DSQRT(0.75D0)
       DO ISPIN=1,NSPIN
         DO IKPT=1,NKPT
           DO IB=1,NB
             SVAR=(X0(IB,IKPT,ISPIN)-XMIN)/(XMAX-XMIN)
             IF(SVAR.GE.0.D0.AND.SVAR.LT.1.D0) CYCLE
             ISVAR=INT(SVAR+100.D0)-100
             SVAR=(XMAX-XMIN)*REAL(ISVAR,KIND=8)
             X0(IB,IKPT,ISPIN)=X0(IB,IKPT,ISPIN)-SVAR
             XM(IB,IKPT,ISPIN)=XM(IB,IKPT,ISPIN)-SVAR
           ENDDO
         ENDDO
       ENDDO
!
!      =================================================================
!      ==  OFFSET OCCUPATIONS SLIGTLY AWAY FROM ZERO                  ==
!      =================================================================
       DO ISPIN=1,NSPIN
         DO IKPT=1,NKPT
           DO IB=1,NB
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
!      ==  CALCULATE FORCE ON X                                       ==
!      ==  REMARK:  THE FACTOR FMAX IS NOT USED BECAUSE THE SAME TERM ==
!      ==         APPEARS ALSO IN THE KINETIC ENERGY                  ==
!      =================================================================
!1000 CONTINUE
       DO ISPIN=1,NSPIN 
         SIGMA=DBLE(3-2*ISPIN)   ! SPIN DIRECTION       
         DO IKPT=1,NKPT
           DO IB=1,NB
!            == FORCE FROM BANDS ======================================
             CALL DYNOCC_FOFX(X0(IB,IKPT,ISPIN),OCC,DOCCDX)
             FORCE=-EPSILON(IB,IKPT,ISPIN)*DOCCDX
!            FORCE=FORCE+(TOTPOT+SIGMA*SPINPOT)*DOCCDX
             IF(TMPSIDOT2) FORCE=FORCE+MPSIDOT2(IB,IKPT,ISPIN)*DOCCDX
             FX(IB,IKPT,ISPIN)=FORCE
!            == FORCE FROM ENTROPY TERM ===============================
             CALL DYNOCC_SOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
             FX(IB,IKPT,ISPIN)=FX(IB,IKPT,ISPIN)-TEMP*DSVAR
           ENDDO
         ENDDO
       ENDDO
!CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!IF (THISTASK.EQ.1) THEN
!  OPEN(10,FILE='FOFX.XY',STATUS='UNKNOWN',FORM='FORMATTED')
!  DO IKPT = 1, NKPT
!    DO IB = 1, NB
!      WRITE(10,'(2F25.15)') EPSILON(IB,IKPT,1),FX(IB,IKPT,1)
!    END DO
!  END DO
!  CLOSE(10)
!END IF
!
!      =================================================================
!      ==  PROPAGATE X WITOUT CONSTRAINTS                             ==
!      =================================================================
!       SVAR1=2.D0/(1.D0+ANNEX)
!       SVAR2=1.D0-SVAR1
!       SVAR3=(SVAR1/2.D0)*DELTAT**2/MX
!       XP(:,:,:)=SVAR1*X0(:,:,:)+SVAR2*XM(:,:,:)+SVAR3*FX(:,:,:)
!
!      =================================================================
!      ==  PROPAGATE X WITOUT CONSTRAINTS                             ==
!      ==  CONSTRAINT FORCE FX                                        ==
!      ==  XP=XBAR+FX*(TOTPOT+SIGMA*SPINPOT)                          ==
!      =================================================================
       SVAR1=2.D0/(1.D0+ANNEX)
       SVAR2=1.D0-SVAR1
       SVAR3=(SVAR1/2.D0)*DELTAT**2/MX
       DO ISPIN=1,NSPIN
         SIGMA=DBLE(3-2*ISPIN)   ! SPIN DIRECTION       
         DO IKPT=1,NKPT
           DO IB=1,NB
             XBAR(IB,IKPT,ISPIN)=SVAR1*X0(IB,IKPT,ISPIN) &
       &                        +SVAR2*XM(IB,IKPT,ISPIN) &
       &                        +SVAR3*FX(IB,IKPT,ISPIN)
             CALL DYNOCC_FOFX(X0(IB,IKPT,ISPIN),SVAR,DSVAR)
             FX(IB,IKPT,ISPIN)=DSVAR*SVAR3
!            XBAR(IB,IKPT,ISPIN)=XP(IB,IKPT,ISPIN) &
!      &           -FX(IB,IKPT,ISPIN)*(TOTPOT+SIGMA*SPINPOT)
           ENDDO
         ENDDO
       ENDDO
!IF (THISTASK.EQ.1) THEN
!  OPEN(10,FILE='FX_B4_CONSTRAINTS.XY',STATUS='UNKNOWN',FORM='FORMATTED')
!  DO IKPT = 1, NKPT
!    DO IB = 1, NB
!      WRITE(10,'(2F25.15)') EPSILON(IB,IKPT,1),XBAR(IB,IKPT,1)
!    END DO
!  END DO
!  CLOSE(10)
!END IF
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
!IF (THISTASK.EQ.1) THEN
!  OPEN(10,FILE='FX_AFTER_CONSTRAINTS.XY',STATUS='UNKNOWN',FORM='FORMATTED')
!  DO IKPT = 1, NKPT
!    DO IB = 1, NB
!      WRITE(10,'(2F25.15)') EPSILON(IB,IKPT,1),XP(IB,IKPT,1)
!    END DO
!  END DO
!  CLOSE(10)
!END IF

!IF(Q1.LT.0.D0) THEN
!  PRINT*,'CLEMENS: Q1 SMALLER THAN ZERO:',ITER
!  COUNT = COUNT + 1
!  IF (COUNT.GT.1000) THEN
!    CALL ERROR$MSG('PETERS FIX RESULTED IN ENDLESS LOOP')
!    CALL ERROR$STOP('PETERS FIX')
!  END IF
!  GOTO 1000
!ENDIF
!COUNT = 0
!CALL MPE$QUERY('MONOMER',NTASKS,THISTASK)
!WRITE(*,FMT='("DYNOCC",2I4,6F20.15)')THISTASK,ITER,TOTPOT,Q0-TOTCHA,Q1
!
!        == CHECK CONVERGENCE ==========================================
         TCONV=.TRUE.
         IF(TFIXTOT) TCONV=TCONV.AND.(ABS(Q0-TOTCHA).LT.QTOL) 
         IF(TFIXSPIN)TCONV=TCONV.AND.(ABS(S0-SPINCHA).LT.QTOL)
         IF(TCONV) EXIT
!        ===============================================================
!        == FIND ZERO OF LINARIZED Q(X)-TOTCHA; S(X)-SPINCHA ===========
!        ==  Q0+Q1*DTOTPOT+S1*DSPINPOT=0                              ==
!        ==  S0+S1*DTOTPOT+Q1*DSPINPOT=0                              ==
!        ===============================================================
         IF(TFIXTOT.AND.TFIXSPIN) THEN
!           DTOTPOT=-(Q0+S0)/(Q1+S1)
!           IF(DABS(Q0+S0).LT.TOL)DTOTPOT=0.D0
!           DSPINPOT=-(Q0-S0)/(Q1-S1)
!           IF(DABS(Q0-S0).LT.TOL)DSPINPOT=0.D0
!           SVARXSXB=0.5D0*(DTOTPOT+DSPINPOT)
!           DSPINPOT=0.5D0*(DTOTPOT-DSPINPOT)
!           DTOTPOT=SVAR
           TOTPOT =TOTPOT -(+Q1*(Q0-TOTCHA)-S1*(S0-SPINCHA))/(Q1**2-S1**2)
           SPINPOT=SPINPOT-(-S1*(Q0-TOTCHA)+Q1*(S0-SPINCHA))/(Q1**2-S1**2)
         ELSE IF(TFIXTOT.AND.(.NOT.TFIXSPIN)) THEN
           TOTPOT=TOTPOT-(Q0-TOTCHA)/Q1
!           IF(DABS(Q0).LT.TOL)DTOTPOT=0.D0
!           DSPINPOT=0.D0
         ELSE IF((.NOT.TFIXTOT).AND.TFIXSPIN) THEN
            SPINPOT=SPINPOT-(S0-SPINCHA)/S1
!           IF(DABS(S0).LT.TOL)DSPINPOT=0.D0
!           DTOTPOT=0.D0
         ELSE IF((.NOT.TFIXTOT).AND.(.NOT. TFIXSPIN)) THEN
!           DTOTPOT =0.D0
!           DSPINPOT=0.D0
         END IF
!        == UPDATE CHEMICAL POTENTIALS =================================
!         IF(TFIXTOT) THEN
!           TOTPOT=TOTPOT+DTOTPOT
!         ELSE
!           TOTCHA=Q0
!         END IF
!         IF(TFIXSPIN) THEN
!           SPINPOT=SPINPOT+DSPINPOT
!         ELSE
!           SPINCHA=S0
!         END IF
       ENDDO
       IF(.NOT.TCONV) THEN
!CALL MPE$QUERY('NONE',NTASKS,THISTASK)
!WRITE(*,FMT='("DYNOCC",2I4,6F20.15)')THISTASK,ITER,Q0-TOTCHA,TOTPOT,Q1
         CALL ERROR$MSG('LOOP NOT CONVERGED') 
         CALL ERROR$R8VAL('DEVIATION TOTAL CHARGE',Q0)
         CALL ERROR$R8VAL('DEVIATION TOTAL SPIN',S0)
         CALL ERROR$STOP('DYNOCC$PROPAGATE')
       END IF
       IF(.NOT.TFIXTOT)TOTCHA=Q0
       IF(.NOT.TFIXSPIN)SPINCHA=S0
!
!      =================================================================
!      == ADJUST CHEMICAL POTENTIALS AND CHARGES                      ==
!      == AND TEST CHARGE AND SPIN CONSERVATION                       ==
!      =================================================================
!
!PRINT*,'SPIN   ',SPINCHA,' SPINPOT ',SPINPOT
!PRINT*,'CHARGE ',TOTCHA ,' EFERMI  ',TOTPOT
!
                              CALL TRACE$POP
       RETURN
       END
!     ..................................................................
      SUBROUTINE DYNOCC$SWITCH()
!     ******************************************************************
!     ******************************************************************
      USE DYNOCC_MODULE
      IMPLICIT NONE
!     ******************************************************************
      TEPSILON=.FALSE.
      IF(ALLOCATED(MPSIDOT2)) DEALLOCATE(MPSIDOT2)
      IF(TDYN) THEN
        XM(:,:,:)=X0(:,:,:)
        X0(:,:,:)=XP(:,:,:)
        XP(:,:,:)=0.D0
      END IF
      RETURN
      END
!
!      .................................................................
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
!      .................................................................
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
!      .................................................................
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
!      .................................................................
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
       REAL(8)                  :: SVAR,X1,X2,X3
       INTEGER(4)               :: IB,IKPT,ISPIN
       INTEGER(4)               :: ISPINDEG
       REAL(8)                  :: EMERMN
       REAL(8)                  :: SIGMA
       REAL(8)                  :: DSVAR,SVAR1
       REAL(8)                  :: WORK(NB,NKPT)
       REAL(8)                  :: EV
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
           CALL DYNOCC_MERMIN(NB,NKPT,NSPIN,NB,NKPT,NSPIN &
     &            ,TOTCHA,ISPINDEG,TEMP,WKPT,EPSILON,X0,TOTPOT,EMERMN)
         ELSE IF(NSPIN.EQ.2) THEN ! NOW THE SPIN POLARIZED CASE
           IF(TFIXSPIN) THEN
             DO ISPIN=1,NSPIN
               SIGMA=DBLE(3-2*ISPIN)
               CALL DYNOCC_MERMIN(NB,NKPT,1,NB,NKPT,1 &
     &            ,0.5D0*(TOTCHA+SIGMA*SPINCHA),ISPINDEG,TEMP &
     &            ,WKPT,EPSILON(:,:,ISPIN),X0(:,:,ISPIN),TOTPOT,EMERMN)
             ENDDO
           ELSE 
             DO ISPIN=1,NSPIN
               SIGMA=DBLE(3-2*ISPIN)
               WORK(:,:)=EPSILON(:,:,ISPIN)+SPINPOT*SIGMA
               CALL DYNOCC_MERMIN(NB,NKPT,1,NB,NKPT,1 &
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
                 IF(DABS(SVAR).LT.55.D0) THEN
                   X0(IB,IKPT,ISPIN)=1.D0/(1.D0+DEXP(SVAR))
                 ELSE
                   X0(IB,IKPT,ISPIN)=0.5D0*(1.D0+DSIGN(1.D0,SVAR))
                 END IF
               ENDDO
             ENDDO
           ENDDO
           X0(:,:,:)=X0(:,:,:)*2.D0/DBLE(NSPIN)
         ELSE 
           CALL DYNOCC_MERMIN(NB,NKPT,NSPIN,NB,NKPT,NSPIN &
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
!WRITE(*,FMT='(3I3,2F10.5)')IB,IKPT,ISPIN,EPSILON(IB,IKPT,ISPIN)*27.211D0,X0(IB,IKPT,ISPIN)
             SVAR=X0(IB,IKPT,ISPIN)/FMAX
             CALL DYNOCC_XOFF(SVAR,X0(IB,IKPT,ISPIN),DSVAR)
           ENDDO
         ENDDO
       ENDDO
       RETURN
       END
!PB031019END
!
!     .....................................................MERMIN.......
      SUBROUTINE DYNOCC_MERMIN(NX,NKPTX,NSPINX,NBANDS,NKPT,NSPIN &
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
      LOGICAL(4) ,PARAMETER  :: TPR=.TRUE.
      INTEGER(4) ,PARAMETER  :: ITERX=1000   ! MAX #(ITERATIONS)
      REAL(8)    ,PARAMETER  :: TOL=1.D-6    ! TOLERANCE IN #(ELECTRONS)
      INTEGER(4) ,INTENT(IN) :: NBANDS,NX    ! #(BANDS),MAX
      INTEGER(4) ,INTENT(IN) :: NKPT,NKPTX   ! #(K-POINTS),MAX
      INTEGER(4) ,INTENT(IN) :: NSPIN,NSPINX ! #(SPINS),MAX
      REAL(8)    ,INTENT(IN) :: TOTCHA       ! TOTAL CHARGE
      INTEGER(4) ,INTENT(IN) :: ISPINDEG     ! SPIN DEGENERACY (1 OR 2)
      REAL(8)    ,INTENT(IN) :: TEMP         ! K_B*TEMPERATURE IN HARTREE
      REAL(8)    ,INTENT(IN) :: WKPT(NKPTX)         ! K-POINT WEIGHTS
      REAL(8)    ,INTENT(IN) :: EIG(NX,NKPTX,NSPINX) ! EIGENVALUES
      REAL(8)    ,INTENT(OUT):: F(NX,NKPTX,NSPINX)   ! OCCUPATIONS
      REAL(8)    ,INTENT(OUT):: CHMPOT               ! CHEMICAL POTENTIAL
      REAL(8)    ,INTENT(OUT):: EMERMN               ! HEAT OF THE ELECTRONS
      INTEGER(4)             :: ISTART
      INTEGER(4)             :: IB,I,ISPIN,IKPT
      REAL(8)                :: X0,DX,Y0,XM,YM
      REAL(8)                :: SVAR
      REAL(8)                :: SUM
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
      IF(ISPINDEG.NE.1.AND.ISPINDEG.NE.2.D0) THEN
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
        SUM=0.D0
        DO ISPIN=1,NSPIN
          DO IKPT=1,NKPT
            DO IB=1,NBANDS
              SVAR=(EIG(IB,IKPT,ISPIN)-CHMPOT)/TEMP
              IF(SVAR.GT.+50.D0)SVAR=+50.D0
              IF(SVAR.LT.-50.D0)SVAR=-50.D0
              F(IB,IKPT,ISPIN)=1.D0/(1.D0+DEXP(SVAR))*FMAX
              SUM=SUM+F(IB,IKPT,ISPIN)*WKPT(IKPT)
            ENDDO
          ENDDO
        ENDDO
        DQ=SUM-TOTCHA
        IF(DABS(DQ).LT.TOL) GOTO 110
        X0=CHMPOT
        Y0=DQ
        CALL BISEC(ISTART,IBI,X0,Y0,DX,XM,YM)
        CHMPOT=X0
      ENDDO
      CALL ERROR$MSG('OCCUPATIONS NOT CONVERGED')
      CALL ERROR$STOP('MERMIN')
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
              DE=+TEMP*(F1*DLOG(F1)+(1.D0-F1)*DLOG(1.D0-F1))
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
            WRITE(*,FMT='(5("(",F8.3,";",F4.2,")"))') &
     &         (EIG(IB,IKPT,ISPIN)/EV,F(IB,IKPT,ISPIN),IB=1,NBANDS)
          ENDDO
        ENDDO
      END IF
                         CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
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
!     ..................................................................
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
!     ..................................................................
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
!     ..................................................................
      SUBROUTINE OCCUPATION$REPORT(NFIL,STRING_)
!     ******************************************************************
!     **                                                              **
!     ******************************************************************
      USE DYNOCC_MODULE, ONLY : NSPIN,NKPT,NB,XK
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: STRING_
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4)              :: ISPIN,IKPT,IB,I,J
      REAL(8)                 :: EV
      REAL(8)                 :: DWORK(3)
      REAL(8)                 :: CELLVOL
      REAL(8)                 :: RBAS(3,3)
      REAL(8)                 :: GBAS(3,3)
      REAL(8)    ,ALLOCATABLE :: EIG(:,:,:)  !(NB,NKPT,NSPIN)
      REAL(8)    ,ALLOCATABLE :: OCC(:,:,:)  !(NB,NKPT,NSPIN)
!     ******************************************************************
      IF(STRING_.EQ.'KPOINTS') THEN
        WRITE(NFIL,FMT='(/"K-POINTS"/"========")')
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
     &             //'," = (",F5.2,"*G1,",F5.2,"*G2,",F5.2,"*G3)"' &
     &             //'," = (",F7.5,",",F7.5,",",F7.5,");")') &
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
      
