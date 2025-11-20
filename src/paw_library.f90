!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****  INTERFACES-ROUTINES FOR EXTERNAL LIBRARIES                          *****
!****                                                                      *****
!****  ALL EXTERNAL LIBRARIES ARE ONLY CALLED VIA THESE INTERFACES.        *****
!****                                                                      *****
!****  SPECIFIC LIBRARIES ARE SELECTED BY THE PREPROCESSOR STATEMENTS      *****
!****                                                                      *****
!****  CONTAINS:                                                           *****
!****    1) FORTRAN UTILITY LIBRARIES                                      *****
!****    2) LAPACK/BLAS LIBRARIES                                          *****
!****    3) FOURIER TRANSFORM LIBRARIES                                    *****
!****    4) SOME OTHER                                                     *****
!****                                                                      *****
!****  THE INTERFACES SHALL BE TESTED USING LIB$TEST                       *****
!****                                                                      *****
!****  NOTE: THERE IS A LAPACK INTERFACE TO ESSL                           *****
!****         HTTP://WWW.NETLIB.ORG/LAPACK/LAWNSPDF/LAWN82.PDF             *****
!****         HTTP://WWW.NETLIB.ORG/LAPACK/ESSL/                           *****
!****                                                                      *****
!****  NOTE: THERE ARE WRAPPER ROUTINES FOR LAPACK AND FOR FFTW3 TO ESSL   *****
!****        SEARCH FOR 'MIGRATING FROM OTHER LIBRARIES TO ESSL'           *****
!****        ON HTTPS://WWW.IBM.COM/DOCS/EN/ESSL/                          *****
!****                                                                      *****
!*******************************************************************************
!*******************************************************************************
!****         
!****     FFTW3 IS THE DEFAULT    
!****     CPPVAR_BLAS_ESSL          
!****     CPPVAR_FFT_ESSL           
!****     CPPVAR_FFT_ACML           
!****     CPPVAR_FFT_PACK           
!****                               
!****     CPPVAR_USAGE_EXIST        
!****         
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$GETOPTS(OPTIONSTRING,NAME,OPTARG,OPTIND,EXITCODE)
!     **************************************************************************
!     **  INTERFACE FOR THE ROUTINE LIB_GETOPTS (ORIGINALLY GETOPTS)          **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*) ,INTENT(IN)   :: OPTIONSTRING
      CHARACTER(*) ,INTENT(OUT)  :: NAME
      CHARACTER(*) ,INTENT(OUT)  :: OPTARG
      INTEGER(4)   ,INTENT(INOUT):: OPTIND
      INTEGER(4)   ,INTENT(OUT)  :: EXITCODE
!     **************************************************************************
      CALL LIB_GETOPTS(OPTIONSTRING,NAME,OPTARG,OPTIND,EXITCODE)
      RETURN
      END SUBROUTINE LIB$GETOPTS
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_GETOPTS(OPTIONSTRING,NAME,OPTARG,OPTIND,EXITCODE)
!     **************************************************************************
!     **                                                                      **
!     ** $$$  SUBPROGRAM DOCUMENTATION BLOCK                                  **
!     **                                                                      **
!     **  SUBPROGRAM:  GETOPTS                                                **
!     **               PROCESS COMMAND LINE ARGUMENTS FOR VALID OPTIONS       **
!     **    PRGMMR: IREDELL       ORG: W/NP23        DATE: 2000-08-22         **
!     **                                                                   
!     **  ABSTRACT: THIS SUBPROGRAM PROCESSES COMMAND-LINE ARGUMENTS FOR      **
!     **     VALID OPTIONS. IT IS THE FORTRAN EQUIVALENT OF THE BUILT-IN      **
!     **     SHELL COMMAND GETOPTS.                                           **
!     **     OPTIONS ON THE COMMAND LINE COME BEFORE THE POSITIONAL ARGUMENTS.**
!     **     OPTIONS ARE PRECEDED BY A - (MINUS SIGN) OR A + (PLUS SIGN).     **
!     **     OPTIONS ARE SINGLE CASE-SENSITIVE ALPHANUMERIC CHARACTERS.       **
!     **     OPTIONS EITHER DO OR DO NOT HAVE AN EXPECTED ARGUMENT.           **
!     **     OPTIONS WITHOUT AN ARGUMENT MAY BE IMMEDIATELY SUCCEEDED BY      **
!     **     FURTHER OPTIONS WITHOUT THE ACCOMPANYING - OR + PREFIX.          **
!     **     OPTIONS WITH AN ARGUMENT MAY BE SEPARATED FROM THEIR ARGUMENT    **
!     **     BY ZERO OR MORE BLANKS.  THE ARGUMENT CANNOT END WITH A BLANK.   **
!     **     OPTIONS END WHEN NOT PRECEDED BY A - OR A + OR AFTER -- OR ++.   **
!     **     THIS SUBPROGRAM PROCESSES ONE OPTION PER INVOCATION.             **
!     **     THIS SUBPROGRAM IS NOT THREAD-SAFE.                              **
!     **                                                                      **
!     **  PROGRAM HISTORY LOG:                                                **
!     **    2000-08-22  IREDELL                                               **
!     **                                                                      **
!     **  USAGE:    CALL GETOPTS(OPTIONSTRING,NAME,OPTARG,OPTIND,EXITCODE)    **
!     **                                                                      **
!     **    INPUT ARGUMENT LIST:                                              **
!     **      OPTIONSTRING                                                    **
!     **        CHARACTER STRING CONTAINING A LIST OF ALL VALID OPTIONS;      **
!     **        OPTIONS SUCCEEDED BY A : REQUIRE AN ARGUMENT                  **
!     **                                                                      **
!     **     INPUT AND OUTPUT ARGUMENT LIST:                                  **
!     **       OPTIND                                                         **
!     **         INTEGER INDEX OF THE NEXT ARGUMENT TO BE PROCESSED;          **
!     **         SET TO 0 BEFORE INITIAL CALL OR TO RESTART PROCESSING        **
!     **                                                                      **
!     **     OUTPUT ARGUMENT LIST:                                            **
!     **       NAME                                                           **
!     **         CHARACTER STRING CONTAINING THE NAME OF THE NEXT OPTION      **
!     **         OR ? IF NO OPTION OR AN UNKNOWN OPTION IS FOUND              **
!     **         OR : IF AN OPTION HAD A MISSING REQUIRED ARGUMENT;           **
!     **         A + IS PREPENDED TO THE VALUE IN NAME IF THE OPTION BEGINS   **
!     **                                                             WITH A + **
!     **       OPTARG                                                         **
!     **         CHARACTER STRING CONTAINING THE OPTION ARGUMENT IF REQUIRED  **
!     **         OR NULL IF NOT REQUIRED OR NOT FOUND;                        **
!     **         OPTARG CONTAINS THE OPTION FOUND IF NAME IS ? OR :.          **
!     **       EXITCODE                                                       **
!     **         INTEGER RETURN CODE (0 IF AN OPTION WAS FOUND, 1 IF END OF   **
!     **         OPTIONS)                                                     **
!     **       
!     **   SUBPROGRAMS CALLED:
!     **     IARGC
!     **       RETRIEVE NUMBER OF COMMAND-LINE ARGUMENTS
!     **     GETARG
!     **       RETRIEVE A COMMAND-LINE ARGUMENT
!     **     INDEX 
!     **       RETRIEVE THE STARTING POSITION OF A SUBSTRING WITHIN A STRING
!     **       
!     **   REMARKS:
!     **     HERE IS AN EXAMPLE OF HOW TO USE THIS SUBPROGRAM.
!     **       IMPLICIT NONE
!     **       CHARACTER*8 COPT,CARG,CB,CPOS
!     **       INTEGER IA,IB,IOPT,IRET,NARG,NPOS,IPOS,IARGC
!     **       IA=0
!     **       IB=0
!     **       IOPT=0
!     **       DO
!     **         CALL GETOPTS('AB:',COPT,CARG,IOPT,IRET)
!     **         IF(IRET.NE.0) EXIT
!     **         SELECT CASE(COPT)
!     **         CASE('A','+A')
!     **           IA=1
!     **         CASE('B','+B')
!     **           IB=1
!     **           CB=CARG
!     **         CASE('?',':')
!     **           PRINT *,'INVALID OPTION ',CARG(1:1)
!     **           STOP 1
!     **         END SELECT
!     **       ENDDO
!     **       IF(IA.EQ.1) PRINT *,'OPTION A SELECTED'
!     **       IF(IB.EQ.1) PRINT *,'OPTION B SELECTED; ARGUMENT=',CB
!     **       NARG=IARGC()
!     **       NPOS=NARG-IOPT+1
!     **       DO IPOS=1,NPOS
!     **         CALL GETARG(IPOS+IOPT-1,CPOS)
!     **         PRINT *,'POSITIONAL ARGUMENT ',IPOS,' IS ',CPOS
!     **       ENDDO
!     **       END
!     **  
!     **   ATTRIBUTES:
!     **     LANGUAGE: FORTRAN 90
!     **  
!     **  $$$
!     **                                                                      **
!     **  RETRIEVED FROM HTTP://WWW.NCO.NCEP.NOAA.GOV/PMB/CODES/NWPROD/       **
!     **                      CFS.V2.1.4/SORC/CFS_SIGAVG.FD/GETOPTS.F         **
!     **  ON AUG 21, 2014                                                     **
!     **************************************************************************
      IMPLICIT NONE
!     ==  PASSED DATA
      CHARACTER(*) ,INTENT(IN)   :: OPTIONSTRING
      CHARACTER(*) ,INTENT(OUT)  :: NAME
      CHARACTER(*) ,INTENT(OUT)  :: OPTARG
      INTEGER(4)   ,INTENT(INOUT):: OPTIND
      INTEGER(4)   ,INTENT(OUT)  :: EXITCODE
!     ==  SAVED DATA
      CHARACTER(256),SAVE        :: CARG
      CHARACTER(1)  ,SAVE        :: CONE
      INTEGER       ,SAVE        :: NARG,LARG,LCUR
!     ==  LOCAL DATA
      CHARACTER(1)               :: COPT
      INTEGER                    :: LNAME,LOPT
      INTEGER                    :: LENG
      INTEGER                    :: ST
!     **************************************************************************
!     ==========================================================================
!     == INITIALLY SET SAVED DATA.
!     ==========================================================================
      IF(OPTIND.LE.0) THEN
        OPTIND=0
        NARG=COMMAND_ARGUMENT_COUNT()
        CARG=' '
        CONE=' '
        LARG=0
        LCUR=1
      ENDIF
!     ==========================================================================
!     ==  RETRIEVE NEXT COMMAND-LINE ARGUMENT IF NECESSARY;
!     ==  EXIT IF AT END OF OPTIONS
!     ==========================================================================
      IF(LCUR.GT.LARG) THEN
        OPTIND=OPTIND+1
        IF(OPTIND.GT.NARG) THEN
          NAME='?'
          OPTARG=' '
          EXITCODE=1
          RETURN
        END IF
!       == FORMERLY:  CALL GETARG(OPTIND,CARG)
        CALL GET_COMMAND_ARGUMENT(OPTIND,CARG,LENG,ST)
        IF(ST.NE.0) THEN
          CALL ERROR$MSG('FAILURE COLLECTING COMMAND LINE ARGUMENT')
          CALL ERROR$I4VAL('STATUS',ST)
          CALL ERROR$I4VAL('POSITION',OPTIND)
          CALL ERROR$I4VAL('ACTUAL LENGTH OF ARGUMENT',LENG)
          CALL ERROR$STOP('LIB_GETOPTS')
        END IF
!
        CONE=CARG(1:1)
        LARG=LEN_TRIM(CARG)
        LCUR=2
        IF(LARG.EQ.1.OR.(CONE.NE.'-'.AND.CONE.NE.'+')) THEN
          NAME='?'
          OPTARG=' '
          EXITCODE=1
          RETURN
        ELSE IF(LARG.EQ.2.AND.CARG(2:2).EQ.CONE) THEN
          OPTIND=OPTIND+1
          NAME='?'
          OPTARG=' '
          EXITCODE=1
          RETURN
        ENDIF
      ENDIF
!
!     ==========================================================================
!     ==  FIND NEXT OPTION IN THE LIST; EXIT IF OPTION IS UNKNOWN
!     ==========================================================================
      EXITCODE=0
      COPT=CARG(LCUR:LCUR)
      LCUR=LCUR+1
      LOPT=INDEX(OPTIONSTRING,COPT)
      IF(LOPT.EQ.0) THEN
        NAME='?'
        OPTARG=COPT
        RETURN
      END IF
!
!     ==========================================================================
!     == OPTION FOUND; RETRIEVE ITS ARGUMENT IF REQUESTED
!     ==========================================================================
      IF(CONE.EQ.'-') THEN
        NAME=""
        LNAME=1
      ELSE
        NAME="+"
        LNAME=2
      END IF
      NAME(LNAME:LNAME)=COPT
      OPTARG=' '
      IF(LOPT.LT.LEN(OPTIONSTRING).AND.OPTIONSTRING(LOPT+1:LOPT+1).EQ.':') THEN
      IF(LCUR.GT.LARG) THEN
        OPTIND=OPTIND+1
        IF(OPTIND.GT.NARG) THEN
          NAME=':'
          OPTARG=COPT
          RETURN
        ENDIF
!       == FORMERLY: CALL GETARG(OPTIND,CARG)
        CALL GET_COMMAND_ARGUMENT(OPTIND,CARG,LENG,ST)
        IF(ST.NE.0) THEN
          CALL ERROR$MSG('FAILURE COLLECTING COMMAND LINE ARGUMENT')
          CALL ERROR$I4VAL('STATUS',ST)
          CALL ERROR$I4VAL('POSITION',OPTIND)
          CALL ERROR$I4VAL('ACTUAL LENGTH OF ARGUMENT',LENG)
          CALL ERROR$STOP('LIB_GETOPTS')
        END IF
!
        LARG=LEN_TRIM(CARG)
        LCUR=1
      ENDIF
      OPTARG=CARG(LCUR:LARG)
      LCUR=LARG+1
    ENDIF
    RETURN
    END SUBROUTINE LIB_GETOPTS
! 
!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****              INTERFACES TO THE LAPACK SUBROUTINES                    *****
!****                                                                      *****
!****  THE LINEAR ALGEBRA PACKAGE (LAPACK) ADRESSES THE FOLLOWING PROBLEMS *****
!****  LINEAR EQUATIONS                                                    *****
!****  LINEAR LEAST SQUARES                                                *****
!****  GENERALIZED LINEAR LEAST SQUARES                                    *****
!****  SYMMETRIC EIGENPROBLEMS                                             *****
!****  NONSYMMETRIC EIGENPROBLEMS                                          *****
!****  SINGULAR VALUE DECOMPOSITION (SVD)                                  *****
!****  GENERALIZED SYMMETRIC DEFINITE EIGENPROBLEMS                        *****
!****  GENERALIZED NONSYMMETRIC EIGENPROBLEMS                              *****
!****  GENERALIZED SINGULAR VALUE DECOMPOSITION                            *****
!****                                                                      *****
!****  HERE WE COVER ONLY A SUBSET OF THESE                                *****
!****    MATRIX INVERSION                                                  *****
!****    LINEAR EQUATIONS                                                  *****
!****    SYMMETRIC EIGENVALUE PROBLEMS                                     *****
!****                                                                      *****
!****                                                                      *****
!*******************************************************************************
!*******************************************************************************
      MODULE LAPACKOPTIONS_MODULE
        CHARACTER(6)               :: GENERALEIGENVALUEC8_MODE='ZHEGV'
      END MODULE LAPACKOPTIONS_MODULE
!
!     ..................................................................
      SUBROUTINE LAPACKOPTIONS$SETCH(ID,VAL)
      USE LAPACKOPTIONS_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID
      CHARACTER(*),INTENT(IN):: VAL
!     ******************************************************************
      IF(ID.EQ.'GENERALEIGENVALUEC8_MODE') THEN
        GENERALEIGENVALUEC8_MODE=VAL
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('LAPACKOPTIONS$SETCH')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$INVERTR8(N,A,AINV)
!     **************************************************************************
!     **  INVERTS THE REAL, SQUARE MATRIX A                                   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N         ! DIMENSION OF THE MATRIX
      REAL(8)   ,INTENT(IN) :: A(N,N)    ! MATRIX TO BE INVERTED
      REAL(8)   ,INTENT(OUT):: AINV(N,N) ! INVERTED MATRIX
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)   ,ALLOCATABLE:: RES(:,:)
      REAL(8)               :: DEV
      INTEGER(4)            :: I
!     **************************************************************************
!      == GENERAL MATRIX INVERSE ===============================================
      CALL LIB_LAPACK_DGETRI(N,A,AINV)
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(RES(N,N))
        RES=MATMUL(A,AINV)
        DO I=1,N
          RES(I,I)=RES(I,I)-1.D0
        ENDDO
        DEV=MAXVAL(ABS(RES))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$INVERTR8')
        END IF
        DEALLOCATE(RES)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$INVERTC8(N,A,AINV)
!     **************************************************************************
!     **  INVERTS THE COMPLEX, SQUARE MATRIX A                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N         ! DIMENSION OF THE MATRIX
      COMPLEX(8),INTENT(IN) :: A(N,N)    ! MATRIX TO BE INVERTED
      COMPLEX(8),INTENT(OUT):: AINV(N,N) ! INVERTED MATRIX
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      COMPLEX(8),ALLOCATABLE:: RES(:,:)
      REAL(8)               :: DEV
      INTEGER(4)            :: I
!     **************************************************************************
!      == GENERAL MATRIX INVERSE ===============================================
      CALL LIB_LAPACK_ZGETRI(N,A,AINV)
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(RES(N,N))
        RES=MATMUL(A,AINV)
        DO I=1,N
          RES(I,I)=RES(I,I)-(1.D0,0.D0)
        ENDDO
        DEV=MAXVAL(ABS(RES))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$INVERTR8')
        END IF
        DEALLOCATE(RES)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$PSEUDOINVERTR8(N,A,SMIN,AINV)
!     **************************************************************************
!     **  CONSTRUCTS PSEUDO INVERSE OF A USING SINGULAR VALUE DECOMPOSITION   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N         ! DIMENSION OF THE MATRIX
      REAL(8)   ,INTENT(IN) :: A(N,N)    ! MATRIX TO BE INVERTED
      REAL(8)   ,INTENT(IN) :: SMIN      ! MIN SINGULAR VALUE TO BE CONSIDERED
      REAL(8)   ,INTENT(OUT):: AINV(N,N) ! INVERTED MATRIX
      LOGICAL   ,PARAMETER  :: TTEST=.TRUE.
      COMPLEX(8),ALLOCATABLE:: RES(:,:)
      REAL(8)               :: DEV
      INTEGER(4)            :: I
      REAL(8)               :: U(N,N)
      REAL(8)               :: V(N,N)
      REAL(8)               :: S(N)
!     **************************************************************************
      CALL LIB_LAPACK_DGESVD(N,N,A,U,S,V)
      DO I=1,N
        IF(ABS(S(I)).GT.SMIN) THEN
          S(I)=1.D0/S(I)
        ELSE
          S(I)=0.D0
        END IF
      ENDDO
      DO I=1,N
        U(:,I)=U(:,I)*S(I)
      ENDDO
      AINV(:,:)=TRANSPOSE(MATMUL(U,V))
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(RES(N,N))
        RES=MATMUL(A,AINV)
        DO I=1,N
          RES(I,I)=RES(I,I)-(1.D0,0.D0)
        ENDDO
        DEV=MAXVAL(ABS(RES))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$INVERTR8')
        END IF
        DEALLOCATE(RES)
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$SVDR8(M,N,A,U,S,VT)
!     **************************************************************************
!     **  SINGULAR VALUE DECOMPOSITION OF THE MXN MATRIX A                    **
!     **  A = U * S * VT
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: M         ! DIMENSION OF THE MATRIX
      INTEGER(4),INTENT(IN) :: N         ! DIMENSION OF THE MATRIX
      REAL(8)   ,INTENT(IN) :: A(M,N)    ! MATRIX TO BE INVERTED
      REAL(8)   ,INTENT(OUT):: U(M,M)    ! LEFT HAND ORTHOGONAL VECTORS
      REAL(8)   ,INTENT(OUT):: S(M)      ! SINGULAR VECTORS
      REAL(8)   ,INTENT(OUT):: VT(N,N)   ! TRANSPOSED RIGHT-HAND ORTHGNAL VECTRS
!     **************************************************************************
      CALL LIB_LAPACK_DGESVD(M,N,A,U,S,VT)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$SVDC8(M,N,A,U,S,VT)
!     **************************************************************************
!     **  SINGULAR VALUE DECOMPOSITION OF THE MXN MATRIX A                    **
!     **  A = U * S * VT
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: M         ! DIMENSION OF THE MATRIX
      INTEGER(4),INTENT(IN) :: N         ! DIMENSION OF THE MATRIX
      COMPLEX(8),INTENT(IN) :: A(M,N)    ! MATRIX TO BE INVERTED
      COMPLEX(8),INTENT(OUT):: U(M,M)    ! LEFT HAND ORTHOGONAL VECTORS
      REAL(8)   ,INTENT(OUT):: S(M)      ! SINGULAR VECTORS
      COMPLEX(8),INTENT(OUT):: VT(N,N)   ! ADJUNCT RIGHT-HAND  VECTORS
!     **************************************************************************
      CALL LIB_LAPACK_ZGESVD(M,N,A,U,S,VT)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$MATRIXSOLVER8(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  SOLVES THE LINEAR EQUATION SYSTEM AX=B                      **
!     **  WHERE A IS A(N,M)                                           **
!     **  IF A IS NOT SQUARE, THE EQUATION IS SOLVED IN A LEAST       **
!     **  SQUARE SENSE                                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: SVAR1,SVAR2
!     ******************************************************************
!
!     ==================================================================
!     == SOLVE SPECIAL CASE WITH DIMENSIONS 1 FIRST                   ==
!     ==================================================================
      IF(N.EQ.1.AND.M.EQ.1) THEN
        X(1,:)=B(1,:)/A(1,1)
        RETURN
      END IF
!
!     ==================================================================
!     == CALL LAPACK DRIVER ROUTINES                                  ==
!     ==================================================================
      IF(N.EQ.M) THEN
        CALL LIB_LAPACK_DGESV(N,M,NEQ,A,X,B)
      ELSE IF(N.LT.M) THEN
        CALL LIB_LAPACK_DGELSD(N,M,NEQ,A,X,B)
      ELSE
        CALL ERROR$MSG('SYSTEM OF EQUATIONS IS OVER DETERMINED')
        CALL ERROR$STOP(' LIB$MATRIXSOLVER8')
      END IF
!
!     ==================================================================
!     ==  TEST                                                        ==
!     ==================================================================
      IF(TTEST) THEN
        SVAR1=MAXVAL(ABS(MATMUL(A,X)-B))
        SVAR2=MAXVAL(ABS(A))/MAXVAL(ABS(B))
        IF(SVAR1/SVAR2.GT.1.D-4) THEN
          CALL ERROR$R8VAL('SVAR1',SVAR1)
          CALL ERROR$R8VAL('SVAR2',SVAR2)
          CALL ERROR$R8VAL('SVAR1/SVAR2',SVAR1/SVAR2)
          CALL ERROR$STOP('LIB$MATRIXSOLVER8')
        END IF
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$MATRIXSOLVEC8(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  SOLVES THE LINEAR EQUATION SYSTEM AX=B                      **
!     **  WHERE A IS A(N,M)                                           **
!     **  IF A IS NOT SQUARE, THE EQUATION IS SOLVED IN A LEAST       **
!     **  SQUARE SENSE                                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(OUT):: X(M,NEQ)
      COMPLEX(8),INTENT(IN) :: B(N,NEQ)
!     ******************************************************************
      IF(N.EQ.1.AND.M.EQ.1) THEN
        X(1,:)=B(1,:)/A(1,1)
        RETURN
      END IF
      IF(N.EQ.M) THEN
        CALL LIB_LAPACK_ZGESV(N,M,NEQ,A,X,B)
      ELSE IF(N.LT.M) THEN
        CALL LIB_LAPACK_ZGELSD(N,M,NEQ,A,X,B)
      ELSE
        CALL ERROR$MSG('SYSTEM OF EQUATIONS IS OVER DETERMINED')
        CALL ERROR$STOP(' LIB$MATRIXSOLVEC8')
      END IF
      RETURN
      END SUBROUTINE LIB$MATRIXSOLVEC8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DIAGR8(N,H,E,U)
!     **************************************************************************
!     **                                                                      **
!     **  DIAGONALIZES THE REAL, SQUARE MATRIX H AFTER SYMMETRIZATION         **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                           **
!     **                                                                      **
!     **         U(K,I)*H(K,L)*U(L,J)=DELTA(I,J)*E(I)                         **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE              **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE                  **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES                     **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: DEV
      REAL(8)  ,ALLOCATABLE :: EMAT(:,:)
      INTEGER(4)            :: I
!     **************************************************************************
!
!     ==========================================================================
!     == DIAGONALIZE                                                          ==
!     ==========================================================================
      CALL LIB_LAPACK_DSYEV(N,H,E,U)
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N))
!       == TEST EIGENVALUE EQUATION ============================================
        EMAT(:,:)=0.D0
        DO I=1,N
          EMAT(I,I)=E(I)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(U,EMAT)))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$DIAGR8')
        END IF
!       == TEST ORTHONORMALITY OF EIGENVECTORS =================================
        EMAT=MATMUL(TRANSPOSE(U),U)
        DO I=1,N
          EMAT(I,I)=EMAT(I,I)-1.D0
        ENDDO
        DEV=MAXVAL(ABS(EMAT))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('ORTHONORMALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$DIAGR8')
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DIAGC8(N,H,E,U)
!     **************************************************************************
!     **                                                                      **
!     **  DIAGONALIZES THE HERMITEAN, SQUARE MATRIX H                         **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                           **
!     **                                                                      **
!     **      CONJG(U(K,I))*H(K,L)*U(L,J)=DELTA(I,J)*E(I)                     **
!     **                                                                      **
!     **  REMARKS:                                                            **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE              **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE                  **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES                     **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      COMPLEX(8),INTENT(OUT):: U(N,N)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: DEV
      COMPLEX(8),ALLOCATABLE:: EMAT(:,:)
      INTEGER(4)            :: I
!     **************************************************************************
      IF(N.EQ.1) THEN
        E(1)=REAL(H(1,1),KIND=8)
        U(1,1)=(1.D0,0.D0)
        RETURN
      END IF
!
!     ==========================================================================
!     == DIAGONALIZE                                                          ==
!     ==========================================================================
      CALL LIB_LAPACK_ZHEEV(N,H,E,U)
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N))
!       == TEST EIGENVALUE EQUATION ============================================
        EMAT(:,:)=CMPLX(0.D0,0.D0,KIND=8)
        DO I=1,N
          EMAT(I,I)=CMPLX(E(I),0.D0,KIND=8)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(U,EMAT)))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$DIAGC8')
        END IF
!       == TEST ORTHONORMALITY OF EIGENVECTORS =================================
        EMAT=MATMUL(TRANSPOSE(CONJG(U)),U)
        DO I=1,N
          EMAT(I,I)=EMAT(I,I)-CMPLX(1.D0,0.D0,KIND=8)
        ENDDO
        DEV=MAXVAL(ABS(EMAT))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('ORTHONORMALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$DIAGC8')
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$EIGVALNONHERMITEANC8(N,H,E)
!     **************************************************************************
!     **                                                                      **
!     **  DIAGONALIZES THE NONHERMITEAN, SQUARE MATRIX H                      **
!     **  AND RETURNS EIGENVALUES ONLY
!     **                                                                      **
!     **      U(K,I)^(-1)*H(K,L)*U(L,J)=E(I)                                  **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      COMPLEX(8),INTENT(OUT):: E(N)
      COMPLEX(8)            :: H1(N,N)
      COMPLEX(8)            :: WORK(2*N)
      REAL(8)               :: RWORK(2*N)
      COMPLEX(8)            :: VL(1,1),VR(1,1)
      INTEGER(4)            :: INFO
!     **************************************************************************
      IF(N.EQ.1) THEN
        E(1)=H(1,1)
        RETURN
      END IF
!
!     ==========================================================================
!     == DIAGONALIZE                                                          ==
!     ==========================================================================
      H1=H
      CALL ZGEEV('N','N',N,H1,N,E,VL,1,VR,1,WORK,2*N,RWORK,INFO)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('THE I-TH ARGUMENT HAD AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB$EIGVALNONHERMITEANC8')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('THE QR ALGORITHM FAILED TO COMPUTE ALL THE EIGENVALUES')
        CALL ERROR$MSG('AND NO EIGENVECTORS HAVE BEEN COMPUTED;')
        CALL ERROR$MSG('ELEMENTS AND I+1:N OF W CONTAIN EIGENVALUES,')
        CALL ERROR$MSG('WHICH HAVE CONVERGED.')
        CALL ERROR$I4VAL('I',INFO)
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$STOP('LIB$EIGVALNONHERMITEANC8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$NONHERMITEANDIAGC8(N,H,E,VR)
!     **************************************************************************
!     **  DIAGONALIZES THE NONHERMITEAN, SQUARE MATRIX H                      **
!     **  AND RETURNS EIGENVALUE AND RIGHT EIGENVECTORS                       **
!     **                                                                      **
!     **      U(K,I)^(-1)*H(K,L)*U(L,J)=E(I)                                  **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      COMPLEX(8),INTENT(OUT):: E(N)
      COMPLEX(8),INTENT(OUT):: VR(N,N)
      COMPLEX(8)            :: VL(1,1)
      COMPLEX(8)            :: H1(N,N)
      COMPLEX(8),ALLOCATABLE:: WORK(:)
      COMPLEX(8)            :: WORK1(1)
      REAL(8)               :: RWORK(2*N)
      INTEGER(4)            :: LWORK
      INTEGER(4)            :: INFO
!     **************************************************************************
      IF(N.EQ.1) THEN
        E(1)=H(1,1)
        VR(1,1)=(1.D0,0.D0)
        RETURN
      END IF
!
!     ==========================================================================
!     == DIAGONALIZE                                                          ==
!     ==========================================================================
!     ==========================================================================
!     == MAKE COPY OF MATRIX, BECAUSE IT WILL BE OVERWRITTEN                  ==
!     ==========================================================================
      H1=H
!
!     ==========================================================================
!     == ALLOCATE WORK ARRAY                                                  ==
!     ==========================================================================
      LWORK=-1
      CALL ZGEEV('N','V',N,H1,N,E,VL,N,VR,N,WORK1,LWORK,RWORK,INFO)
      LWORK=INT(WORK1(1))
      ALLOCATE(WORK(LWORK))
!
!     ==========================================================================
!     == SOLVE EIGENVALUE PROBLEM                                             ==
!     ==========================================================================
      CALL ZGEEV('N','V',N,H1,N,E,VL,N,VR,N,WORK,LWORK,RWORK,INFO)
!
!     ==========================================================================
!     == CHECK ERROR MESSAGE                                                  ==
!     ==========================================================================
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('THE I-TH ARGUMENT HAD AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB$EIGVALNONHERMITEANC8')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('THE QR ALGORITHM FAILED TO COMPUTE ALL THE EIGENVALUES')
        CALL ERROR$MSG('AND NO EIGENVECTORS HAVE BEEN COMPUTED;')
        CALL ERROR$MSG('ELEMENTS AND I+1:N OF W CONTAIN EIGENVALUES,')
        CALL ERROR$MSG('WHICH HAVE CONVERGED.')
        CALL ERROR$I4VAL('I',INFO)
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$STOP('LIB$EIGVALNONHERMITEANC8')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$GENERALEIGENVALUER8(N,H,S,E,U)
!     **************************************************************************
!     **                                                                      **
!     ** SOLVES THE GENERALIZED, REAL, SYMMETRIC EIGENVALUE PROBLEM           **
!     **                                                                      **
!     **      H*U = S*U*E                                                     **
!     **                                                                      **
!     ** WITH EIGENVECTORS U THAT ARE ORTHOGONAL IN THE SENSE                 **
!     **                                                                      **
!     **      U^T*S*U=IDENTITY                                                **
!     **                                                                      **
!     ** REMARK: H AND S MUST BE SYMMETRIC                                    **
!     **         S MUST BE POSITIVE DEFINITE                                  **
!     **         EIGENVECTORS ARE ORTHONORMAL IN THE SENSE                    **
!     **             MATMUL(TRANSPOSE(U),MATMUL(S,U))=IDENTITY                **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)    ! HAMILTON MATRIX
      REAL(8)   ,INTENT(IN) :: S(N,N)    ! OVERLAP MATRIX
      REAL(8)   ,INTENT(OUT):: E(N)      ! EIGENVALUES
      REAL(8)   ,INTENT(OUT):: U(N,N)    ! EIGENVECTORS
      REAL(8)               :: B(N,N)    ! COPY OF OVERLAP MATRIX
      LOGICAL               :: TSYM
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE. ! IF TRUE TEST RESULT
      REAL(8)               :: DEV       ! DEVIATION
      INTEGER               :: I 
!     **************************************************************************
!
!     ==========================================================================
!     == TAKE CARE OF TRIVIAL CASES                                           ==
!     ==========================================================================
      IF(N.EQ.1) THEN
        E(1)=H(1,1)/S(1,1)
        U(1,1)=1.D0/SQRT(S(1,1))
        RETURN
      ELSE IF(N.EQ.0) THEN
        RETURN
      END IF
!
!     ==========================================================================
!     == TEST IF INPUT MATRICES ARE SYMMETRIC                                 ==
!     ==========================================================================
      DEV=MAXVAL(ABS(H-TRANSPOSE(H)))
      TSYM=(DEV.LT.1.D-5)
      IF(.NOT.TSYM) THEN
        CALL ERROR$MSG('HAMILTONIAN NOT SYMMETRIC')
        CALL ERROR$R8VAL('DEV',DEV)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
      END IF
      DEV=MAXVAL(ABS(S-TRANSPOSE(S)))
      TSYM=TSYM.AND.(DEV.LT.1.D-5)
      IF(.NOT.TSYM) THEN
        CALL ERROR$MSG('OVERLAP MATRIX NOT SYMMETRIC')
        CALL ERROR$R8VAL('DEV',DEV)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
      END IF
!
!     ==========================================================================
!     == DIAGONALIZE                                                          ==
!     ==========================================================================
      CALL LIB_LAPACK_DSYGV(N,H,S,E,U)
!
!     ==========================================================================
!     == TEST RESULT OF THE ROUTINE                                           ==
!     ==========================================================================
      IF(TTEST) THEN
!       == CHECK ORTHONORMALITY OF EIGENVECTORS ================================
        B=MATMUL(TRANSPOSE(U),MATMUL(S,U))
        DO I=1,N
          B(I,I)=B(I,I)-1.D0
        ENDDO
        DEV=SUM(ABS(B))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('EIGENSTATES NOT ORTHONORMAL')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
        END IF
!       == CHECK EIGENVALUE PROBLEM ============================================
        DO I=1,N
          DEV=SUM(ABS(MATMUL(H-E(I)*S,U(:,I))))
          IF(DEV.GT.1.D-7) THEN
            CALL ERROR$MSG('EIGENVALUE PROBLEM FAILED')
            CALL ERROR$R8VAL('DEV',DEV)
            CALL ERROR$STOP('LIB$GENERALEIGENVALUER8')
          END IF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE LIB$GENERALEIGENVALUER8
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$GENERALEIGENVALUEC8(N,H,S,E,U)
!     **************************************************************************
!     **                                                                      **
!     ** SOLVES THE GENERALIZED, COMPLEX NON-SYMMETRIC EIGENVALUE PROBLEM     **
!     **      [H(:,:)-E(I)*S(:,:)]*U(:,I)=0                                   **
!     **                                                                      **
!     ** REMARK: H AND S MUST BE HERMITEAN                                    **
!     **         S MUST BE POSITIVE DEFINITE                                  **
!     **         EIGENVECTORS ARE ORTHONORMAL IN THE SENSE                    **
!     **             MATMUL(TRANSPOSE(U),MATMUL(S,U))=IDENTITY                **
!     **                                                                      **
!     **************************************************************************
      USE LAPACKOPTIONS_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)               :: N         ! DIMENSION
      COMPLEX(8),INTENT(IN)               :: H(N,N)    ! HAMITON MATRIX
      COMPLEX(8),INTENT(IN)               :: S(N,N)    ! OVERLAP MATRIX
      REAL(8)   ,INTENT(OUT)              :: E(N)      ! EIGENVALUES
      COMPLEX(8),INTENT(OUT)              :: U(N,N)    ! EIGENVECTORS
      COMPLEX(8)                          :: S1(N,N)
      LOGICAL   ,PARAMETER                :: TTEST=.FALSE.
      REAL(8)                             :: DEV
      INTEGER                             :: I
!     **************************************************************************
!
!     ==========================================================================
!     == TEST IF INPUT MATRICES ARE SYMMETRIC                                 ==
!     ==========================================================================
      DEV=SUM(ABS(H-TRANSPOSE(CONJG(H))))
      IF(DEV.GT.1.D-8) THEN
        CALL ERROR$MSG('HAMILTON MATRIX NOT HERMITEAN')
        CALL ERROR$R8VAL('DEV',DEV)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
      END IF
      DEV=SUM(ABS(S-TRANSPOSE(CONJG(S))))
      IF(DEV.GT.1.D-8) THEN
        CALL ERROR$MSG('OVERLAP MATRIX NOT HERMITEAN')
        CALL ERROR$R8VAL('DEV',DEV)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
      END IF
!
!     ========================================================================
!     == CALL LAPACK ROUTINE                                                ==
!     ========================================================================
      IF(GENERALEIGENVALUEC8_MODE.EQ.'ZHEGV')THEN
        CALL LIB_LAPACK_ZHEGV(N,H,S,E,U)
      ELSE IF(GENERALEIGENVALUEC8_MODE.EQ.'ZHEGVD')THEN
        CALL LIB_LAPACK_ZHEGVD(N,H,S,E,U)
      ELSE
        CALL ERROR$MSG('GENERALEIGENVALUEC8_MODE UNKNOWN')
        CALL ERROR$CHVAL('MODE',GENERALEIGENVALUEC8_MODE)
        CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
      ENDIF
!
!     ========================================================================
!     == TEST RESULT OF THE ROUTINE                                         ==
!     ========================================================================
      IF(TTEST) THEN
        S1=MATMUL(TRANSPOSE(CONJG(U)),MATMUL(S,U))
        DEV=0.D0
        DO I=1,N
          S1(I,I)=S1(I,I)-(1.D0,0.D0)
          DEV=MAX(DEV,MAXVAL(ABS(MATMUL(H-E(I)*S,U(:,I)))))
        ENDDO
        IF(DEV.GT.1.D-6) THEN
          CALL ERROR$MSG('GENERAL EIGENVALUE TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
        END IF
        DEV=SUM(ABS(S1))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('EIGENSTATES NOT ORTHONORMAL')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB$GENERALEIGENVALUEC8')
        END IF
      END IF
      RETURN
      END
!***********************************************************************
!***********************************************************************
!****                                                               ****
!****  DRIVER ROUTINES FOR LAPACK                                   ****
!****                                                               ****
!***********************************************************************
!***********************************************************************
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_DGETRI(N,A,AINV)
!     ******************************************************************
!     **                                                              **
!     **  INVERTS THE REAL, SQUARE MATRIX A                           **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: DGEICD                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: A(N,N)
      REAL(8)   ,INTENT(OUT):: AINV(N,N)
      INTEGER(4)            :: NAUX
      REAL(8)               :: AUX(100*N)
      INTEGER(4)            :: IPIV(N)
      INTEGER(4)            :: INFO
!     ******************************************************************
      NAUX=100*N
      AINV(1:N,1:N)=A(1:N,1:N)
!
!     ==================================================================
!     == PERFORM LU FACTORIZATION OF A                                ==
!     ==================================================================
      CALL DGETRF(N,N,AINV,N,IPIV,INFO) !LAPACK
!
!     ==================================================================
!     == INVERT A USING THE LU FACTORIZATION                          ==
!     ==================================================================
      CALL DGETRI(N,AINV,N,IPIV,AUX,NAUX,INFO) !LAPACK
!
!     ==================================================================
!     == CHECK ERROR CODE                                             ==
!     ==================================================================
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_LAPACK_DGETRI')
        ELSE
          CALL ERROR$MSG('U(I,I) IS EXACTLY ZERO')
          CALL ERROR$MSG('MATRIX IS SINGULAR: ITS INVERSE COULDNOT BE COMPUTED')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_DGETRI')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_DGETRI
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_ZGETRI(N,A,AINV)
!     ******************************************************************
!     **                                                              **
!     **  INVERTS THE REAL, SQUARE MATRIX A                           **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: DGEICD                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: A(N,N)
      COMPLEX(8),INTENT(OUT):: AINV(N,N)
      INTEGER(4)            :: NAUX
      COMPLEX(8)            :: AUX(100*N)
      INTEGER(4)            :: IPIV(N)
      INTEGER               :: INFO
!     ******************************************************************
      IF(N.LE.0) RETURN
      NAUX=100*N
      AINV(1:N,1:N)=A(1:N,1:N)
!
!     ==================================================================
!     == PERFORM LU FACTORIZATION OF A                                ==
!     ==================================================================
      CALL ZGETRF(N,N,AINV,N,IPIV,INFO) !LAPACK
!
!     ==================================================================
!     == INVERT A USING THE LU FACTORIZATION                          ==
!     ==================================================================
      CALL ZGETRI(N,AINV,N,IPIV,AUX,NAUX,INFO) !LAPACK
!
!     ==================================================================
!     == CHECK ERROR CODE                                             ==
!     ==================================================================
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_LAPACK_DGETRI')
        ELSE
          CALL ERROR$MSG('U(I,I) IS EXACTLY ZERO')
          CALL ERROR$MSG('MATRIX IS SINGULAR: ITS INVERSE COULDNOT BE COMPUTED')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZGETRI')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_ZGETRI
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_LAPACK_DGESVD(M,N,A,U,S,VT)
!     **************************************************************************
!     ** SINGULAR VALUE DECOMPOSITION OF THE NON-SQUARE MATRIX A              **
!     ** A=U*SIGMA*VT  
!     **  SIGMA IS AN M-TIMES-N DIAGONAL MATRIX WITH DIAGONAL ELEMENTS S(I)   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: M
      INTEGER(4),INTENT(IN)  :: N
      REAL(8)   ,INTENT(IN)  :: A(M,N)
      REAL(8)   ,INTENT(OUT) :: U(M,M)
      REAL(8)   ,INTENT(OUT) :: VT(N,N)
      REAL(8)   ,INTENT(OUT) :: S(M)     ! SINGULAR VALUES IN DESCENDING ORDER
      REAL(8)   ,ALLOCATABLE :: WORK(:)
      REAL(8)                :: ACOPY(M,N)
      INTEGER(4)             :: LWORK
      INTEGER(4)             :: INFO
      INTEGER(4)             :: I
!     **************************************************************************
      LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
      ALLOCATE(WORK(LWORK))
      ACOPY=A  ! ACOPY WILL BE OVERWRITTEN
      S=0.D0
      U=0.D0
      VT=0.D0
      WORK=0.D0
      CALL DGESVD('A','A',M,N,ACOPY,M,S,U,M,VT,N,WORK,LWORK,INFO)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('THE I-TH ARGUMENT JHAS AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$I4VAL('1ST INDEX OF MATRIX',M)
        CALL ERROR$I4VAL('2ND INDEX OF MATRIX',N)
        CALL ERROR$STOP('LIB_LAPACK_DGESVD')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('DBDSQR DID NOT CONVERGE')
        CALL ERROR$I4VAL('NUMBER OF UNCONVERGED SUPERDIAGONALS',INFO)
        CALL ERROR$I4VAL('1ST INDEX OF MATRIX',M)
        CALL ERROR$I4VAL('2ND INDEX OF MATRIX',N)
        CALL ERROR$STOP('LIB_LAPACK_DGESVD')
      END IF
      DEALLOCATE(WORK)
      S(N+1:M)=0.D0
      RETURN
      END SUBROUTINE LIB_LAPACK_DGESVD
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_LAPACK_ZGESVD(M,N,A,U,S,VT)
!     **************************************************************************
!     ** SINGULAR VALUE DECOMPOSITION OF THE NON-SQUARE MATRIX A              **
!     ** A=U*SIGMA*VT  
!     **  SIGMA IS AN M-TIMES-N DIAGONAL MATRIX WITH DIAGONAL ELEMENTS S(I)   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: M
      INTEGER(4),INTENT(IN)  :: N
      COMPLEX(8),INTENT(IN)  :: A(M,N)
      COMPLEX(8),INTENT(OUT) :: U(M,M)
      COMPLEX(8),INTENT(OUT) :: VT(N,N)
      REAL(8)   ,INTENT(OUT) :: S(M)     ! SINGULAR VALUES IN DESCENDING ORDER
      COMPLEX(8),ALLOCATABLE :: WORK(:)
      REAL(8)   ,ALLOCATABLE :: RWORK(:)
      COMPLEX(8)             :: ACOPY(M,N)
      INTEGER(4)             :: LWORK
      INTEGER(4)             :: INFO
      INTEGER(4)             :: I
!     **************************************************************************
      LWORK=MAX(1,2*MIN(M,N)+MAX(M,N))
      ALLOCATE(WORK(LWORK))
      ALLOCATE(RWORK(5*MIN(M,N)))
      ACOPY=A  ! ACOPY WILL BE OVERWRITTEN
      CALL ZGESVD('A','A',M,N,ACOPY,M,S,U,M,VT,N,WORK,LWORK,RWORK,INFO)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('THE I-TH ARGUMENT JHAS AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$I4VAL('1ST INDEX OF MATRIX',M)
        CALL ERROR$I4VAL('2ND INDEX OF MATRIX',N)
        CALL ERROR$STOP('LIB_LAPACK_DGESVD')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('ZBDSQR DID NOT CONVERGE')
        CALL ERROR$I4VAL('NUMBER OF UNCONVERGED SUPERDIAGONALS',INFO)
        CALL ERROR$I4VAL('1ST INDEX OF MATRIX',M)
        CALL ERROR$I4VAL('2ND INDEX OF MATRIX',N)
        CALL ERROR$STOP('LIB_LAPACK_DGESVD')
      END IF
      DEALLOCATE(WORK)
      S(N+1:M)=0.D0
      RETURN
      END SUBROUTINE LIB_LAPACK_ZGESVD
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_DGESV(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  DRIVER ROUTINE FOR LAPACK ROUTINE DGESV                     **
!     **                                                              **
!     **  COMPUTES THE SOLUTION TO A REAL SYSTEM OF LINEAR EQUATIONS  **
!     **        A * X = B,                                            **
!     **  WHERE A IS A N-BY-N MATRIX AND X AND B ARE N-BY-NEQ MATRICES**
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      REAL(8)   ,ALLOCATABLE:: A1(:,:)
      INTEGER               :: INFO
      INTEGER(4)            :: IPIV(N)
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: SVAR1,SVAR2
!     ******************************************************************
      IF(N.NE.M) THEN
        CALL ERROR$MSG('WORKS ONLY FOR SQUARE MATRICES')
        CALL ERROR$STOP('LIB_LAPACK_DGESV')
      END IF
!
!     ==================================================================
!     == SOLVE SPECIAL CASE WITH DIMENSIONS 1 FIRST                   ==
!     ==================================================================
      IF(N.EQ.1.AND.M.EQ.1) THEN
        X(1,:)=B(1,:)/A(1,1)
        RETURN
      END IF
!
!     ==================================================================
!     == NOW CALL LAPACK ROUTINE                                      ==
!     ==================================================================
      ALLOCATE(A1(N,M))
      A1=A
      X=B
      CALL DGESV(N,NEQ,A1,N,IPIV,X,N,INFO )
      DEALLOCATE(A1)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('ERROR EXIT FROM DGESV')
        CALL ERROR$MSG('THE I-TH ARGUMENT HAD AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB_LAPACK_DGESV')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('ERROR EXIT FROM DGESV')
        CALL ERROR$MSG('U(I,I) IS EXACTLY ZERO.')
        CALL ERROR$MSG('THE FACTORIZATION HAS BEEN COMPLETED,')
        CALL ERROR$MSG('BUT THE FACTOR U IS EXACTLY SINGULAR,')
        CALL ERROR$MSG('SO THE SOLUTION COULD NOT BE COMPUTED.')
        CALL ERROR$I4VAL('I',INFO)
        CALL ERROR$STOP('LIB_LAPACK_DGESV')
      END IF
!
!     ==================================================================
!     ==  TEST                                                        ==
!     ==================================================================
      IF(TTEST) THEN
        SVAR1=MAXVAL(ABS(MATMUL(A,X)-B))
        SVAR2=MAXVAL(ABS(A))/MAXVAL(ABS(B))
        IF(SVAR1/SVAR2.GT.1.D-4) THEN
          CALL ERROR$R8VAL('SVAR1',SVAR1)
          CALL ERROR$R8VAL('SVAR2',SVAR2)
          CALL ERROR$R8VAL('SVAR1/SVAR2',SVAR1/SVAR2)
          CALL ERROR$STOP('LIB_LAPACK_DGESV')
        END IF
      END IF

      RETURN
      END SUBROUTINE LIB_LAPACK_DGESV
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_LAPACK_DGELSD(N,M,NEQ,A,X,B)
!     **************************************************************************
!     **  DRIVER ROUTINE FOR LAPACK ROUTINE DGELSD                            **
!     **                                                                      **
!     **  COMPUTES THE MINIMUM-NORM SOLUTION TO A COMPLEX LINEAR LEAST        **
!     **  SQUARES PROBLEM:                                                    **
!     **              MINIMIZE 2-NORM(| B - A*X |)                            **
!     **  USING THE SINGULAR VALUE DECOMPOSITION (SVD) OF A.                  **
!     **  A IS AN M-BY-N MATRIX WHICH MAY BE RANK-DEFICIENT.                  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(OUT):: X(M,NEQ)
      REAL(8)   ,INTENT(IN) :: B(N,NEQ)
      DOUBLE PRECISION            :: A1(N,M)
      DOUBLE PRECISION,ALLOCATABLE:: B1(:,:)
      DOUBLE PRECISION            :: S(N)        ! SINGULAR VALUES
      DOUBLE PRECISION,PARAMETER  :: RCOND=-1.D0 ! USE MACHINE PRECISION
      INTEGER                     :: RANK
      INTEGER                     :: LWORK
      DOUBLE PRECISION,ALLOCATABLE:: WORK(:)
      INTEGER         ,ALLOCATABLE:: IWORK(:)
      INTEGER                     :: LIWORK
      INTEGER                     :: INFO
      INTEGER                     :: NRHS
      INTEGER                     :: LDA,LDB
      INTEGER                     :: N1,M1,SMLSIZ,MINMN,MAXMN,NLVL
      INTEGER                     :: ILAENV
      EXTERNAL                    :: ILAENV
      LOGICAL         ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)         ,PARAMETER  :: TOL=1.D-5
      REAL(8)                     :: SVAR1,SVAR2
      CHARACTER(6)    ,PARAMETER  :: TYPE='DGELS'
!     **************************************************************************
      IF(N.LT.1.OR.M.LT.1) THEN
        CALL ERROR$MSG('DIMENSIONS MUST BE NONZERO AND POSITIVE')
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$I4VAL('M',M)
        CALL ERROR$STOP('LIB_LAPACK_DGELSD')
      END IF
      M1=N    ! LAPACK USES M AND N OPPOSITE 
      N1=M
      NRHS=NEQ
      MINMN=MIN(M1,N1)
      MAXMN=MAX(M1,N1)
!
!     ==========================================================================
!     == INTERFACE TO DGELSD                                                  ==
!     ==========================================================================
      IF(TYPE.EQ.'DGELSD') THEN
        SMLSIZ=ILAENV(9,'DGELSD',' ',0,0,0,0 )
        NLVL = MAX(0,INT(LOG(REAL(MINMN,KIND=8)/REAL(SMLSIZ+1,KIND=8))/LOG(2.D0))+1)
        LIWORK=MAX(1,3*MINMN*NLVL+11*MINMN)
        LWORK=12*MINMN+2*MINMN*SMLSIZ+8*MINMN*NLVL+MINMN*NEQ+(SMLSIZ+1)**2
        ALLOCATE(WORK(LWORK))
        ALLOCATE(IWORK(LIWORK))
        LDB=MAX(1,MAXMN)
        ALLOCATE(B1(LDB,NEQ))
        B1=0.D0
        B1(1:M1,:)=B(:,:)
        A1=A
! THE CALL TO DGELSD FAILED USING MKL
        CALL DGELSD(M1,N1,NRHS,A1,M1,B1,LDB,S,RCOND,RANK,WORK,LWORK,IWORK,INFO )
        DEALLOCATE(IWORK)
        DEALLOCATE(WORK)
        X=B1(1:N1,:)
        DEALLOCATE(B1)
        IF(INFO.NE.0) THEN
          IF(INFO.LT.0) THEN
            CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
            CALL ERROR$I4VAL('I',-INFO)
            CALL ERROR$STOP('LIB_LAPACK_DGELSD')
          ELSE
            CALL ERROR$MSG('ALGORITHM FOR COMPUTING SVD FAILED TO CONVERGE')
            CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
            CALL ERROR$MSG('BIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
            CALL ERROR$I4VAL('I',INFO)
            CALL ERROR$STOP('LIB_LAPACK_DGELSD')
          END IF
        END IF
!
!     ==========================================================================
!     == INTERFACE TO DGELS                                                   ==
!     ==========================================================================
      ELSE IF(TYPE.EQ.'DGELS') THEN
        LDA=MAX(1,M1)
        LDB=MAX(1,MAXMN)
        ALLOCATE(B1(LDB,NEQ))
        B1=0.D0
        B1(1:M1,:)=B(:,:)
        A1=A
        LWORK=MINMN+MAX(1,MAX(NRHS,MAXMN))
        ALLOCATE(WORK(LWORK))
        CALL DGELS('N',M1,N1,NRHS,A1,LDA,B1,LDB,WORK,LWORK,INFO)
        DEALLOCATE(WORK)
        X=B1(1:N1,:)
        DEALLOCATE(B1)
!
!     ==========================================================================
!     == INTERFACE TO DGELSS                                                  ==
!     ==========================================================================
      ELSE IF(TYPE.EQ.'DGELSS') THEN
        LDA=MAX(1,M1)
        LDB=MAX(1,MAX(M1,N1))
        ALLOCATE(B1(LDB,NEQ))
        B1=0.D0
        B1(1:M1,:)=B(:,:)
        A1=A
        LWORK=3*MINMN+MAX(2*MINMN,MAX(NRHS,MAXMN))
        ALLOCATE(WORK(LWORK))
        CALL DGELSS(M1,N1,NRHS,A1,LDA,B1,LDB,S,RCOND,RANK,WORK,LWORK,INFO)
        DEALLOCATE(WORK)
        X=B1(1:N1,:)
        DEALLOCATE(B1)
      ELSE
        CALL ERROR$MSG('TYPE NOT RECOGNIZED') 
        CALL ERROR$STOP('LIB_LAPACK_DGELSD')
      END IF
!
!     ==========================================================================
!     == TEST                                                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        SVAR1=MAXVAL(ABS(MATMUL(A,X)-B))
        SVAR2=TOL*MAXVAL(ABS(B))/MAXVAL(ABS(B))
        IF(SVAR1.GT.SVAR2) THEN
          CALL ERROR$MSG('TEST FAILED')
          CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('BIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
          CALL ERROR$R8VAL('MAX DEV',SVAR1)
          CALL ERROR$R8VAL('ALLOWED DEV',SVAR2)
          CALL ERROR$STOP('LIB_LAPACK_DGELSD')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_DGELSD
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_ZGESV(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  OMPLEX SYSTEM  OF  LINEAR EQUATIONS                         **
!     **              A * X = B,                                      **
!     **  WHERE A IS A(N,M) IS A SQUARE MATRIX                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(OUT):: X(M,NEQ)
      COMPLEX(8),INTENT(IN) :: B(N,NEQ)
      COMPLEX(8)            :: A1(N,M)
      INTEGER               :: INFO
      INTEGER               :: LDWORK
      INTEGER               :: IPIVOT(N)
      REAL(8)   ,ALLOCATABLE:: WORK(:)
      INTEGER               :: N1,M1,NEQ1
!     ******************************************************************
      IF(N.NE.M) THEN
        CALL ERROR$MSG('ONLY SYMMETRIC MATRICES ALLOWED')
        CALL ERROR$STOP('LIB_LAPACK_ZGESV')
      END IF

      LDWORK=3*MIN(M,N)+MAX(2*MIN(M,N),MAX(M,N),NEQ)
!     -- USE 3*M+3*N+NEQ
      N1=N
      M1=M
      NEQ1=NEQ
!
!     ===================================================================
!     == SOLVE EQUATION SYSTEM                                         ==
!     ===================================================================
      A1=A 
      X=B
      ALLOCATE(WORK(LDWORK))
      CALL ZGESV(N,NEQ,A1,N,IPIVOT,X,N,INFO)
      DEALLOCATE(WORK)
!
!     ===================================================================
!     == CHECK ERROR CODE                                              ==
!     ===================================================================
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('ITH ARGUMENT HAS ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LLIB_LAPACK_ZGESV')
        ELSE IF(INFO.GT.0) THEN
          CALL ERROR$MSG('PROBLEM IS SINGULAR. NO SOLUTION CAN BE COMPUTED')
          CALL ERROR$MSG('U(I,I) IS  EXACTLY  ZERO. THE FACTORIZATION HAS BEEN') 
          CALL ERROR$MSG('COMPLETED, BUT THE FACTOR U IS EXACTLY SINGULAR,')
          CALL ERROR$MSG('SO THE SOLUTION COULD NOT  BE COMPUTED.')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZGESV')
        END IF
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_ZGELSD(N,M,NEQ,A,X,B)
!     ******************************************************************
!     **  DRIVER ROUTINE FOR LAPACK ROUTINE ZGELSD                    **
!     **                                                              **
!     **  COMPUTES THE MINIMUM-NORM SOLUTION TO A COMPLEX LINEAR LEAST**
!     **  SQUARES PROBLEM:                                            **
!     **              MINIMIZE 2-NORM(| B - A*X |)                    **
!     **  USING THE SINGULAR VALUE DECOMPOSITION (SVD) OF A.          **
!     **  A IS AN M-BY-N MATRIX WHICH MAY BE RANK-DEFICIENT.          **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: NEQ
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(OUT):: X(M,NEQ)
      COMPLEX(8),INTENT(IN) :: B(N,NEQ)
      COMPLEX(8)            :: A1(N,M)
      COMPLEX(8),ALLOCATABLE:: B1(:,:)
      REAL(8)               :: S(N)        ! SINGULAR VALUES
      REAL(8)   ,PARAMETER  :: RCOND=-1.D0 ! USE MACHINE PRECISION
      INTEGER(4)            :: RANK
      INTEGER(4)            :: LCWORK      ! SIZE OF CWORK
      INTEGER(4)            :: LRWORK      ! SIZE OF RWORK
      INTEGER(4)            :: LIWORK      ! SIZE OF IWORK
      COMPLEX(8),ALLOCATABLE:: CWORK(:)
      REAL(8)   ,ALLOCATABLE:: RWORK(:)
      INTEGER(4),ALLOCATABLE:: IWORK(:)
      INTEGER(4)            :: INFO
      INTEGER               :: LDB
      INTEGER               :: N1,M1,SMLSIZ,MINMN,NLVL
      INTEGER               :: ILAENV
      EXTERNAL              :: ILAENV
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)   ,PARAMETER  :: TOL=1.D-5
      REAL(8)               :: SVAR1,SVAR2
!     ******************************************************************
      M1=N    ! LAPACK USES M AND N OPPOSITE 
      N1=M
      SMLSIZ=ILAENV(9,'ZGELSD',' ',0,0,0,0 )
      MINMN=MIN(M1,N1)
      NLVL = MAX(0,INT(LOG(REAL(MINMN/(SMLSIZ+1),KIND=8))/LOG(2.D0))+1)
      LCWORK=2*MINMN+MINMN*NEQ
      LRWORK=10*MINMN+2*MINMN*SMLSIZ+8*MINMN*NLVL+3*SMLSIZ*NEQ+(SMLSIZ+1)**2
      LIWORK=MAX(1,3*MINMN*NLVL+11*MINMN)
      ALLOCATE(CWORK(LCWORK))
      ALLOCATE(RWORK(LRWORK))
      ALLOCATE(IWORK(LIWORK))
      LDB=MAX(1,MAX(M1,N1))
      ALLOCATE(B1(LDB,NEQ))
      B1=0.D0
      B1(1:M1,:)=B(:,:)
      A1=A
      CALL ZGELSD(M1,N1,NEQ,A1,M1,B1,LDB,S,RCOND,RANK &
     &           ,CWORK,LCWORK,RWORK,IWORK,INFO )
      X=B1(1:N1,:)
      DEALLOCATE(CWORK)
      DEALLOCATE(RWORK)
      DEALLOCATE(IWORK)
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZGELSD')
        ELSE
          CALL ERROR$MSG('ALGORITHM FOR COMPUTING SVD FAILED TO CONVERGE')
          CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('BIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZGELSD')
        END IF
      END IF
!     ==================================================================
!     ==                                                              ==
!     ==================================================================
      IF(TTEST) THEN
        SVAR1=MAXVAL(ABS(MATMUL(A,X)-B))
        SVAR2=TOL*MAXVAL(ABS(B))/MAXVAL(ABS(B))
        IF(SVAR1.GT.SVAR2) THEN
          CALL ERROR$MSG('TEST FAILED')
          CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('BIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
          CALL ERROR$R8VAL('MAX DEV',SVAR1)
          CALL ERROR$R8VAL('ALLOWED DEV',SVAR2)
          CALL ERROR$STOP('LIB_LAPACK_ZGELSD')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_ZGELSD
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_DSYEV(N,H,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE REAL, SQUARE MATRIX H AFTER SYMMETRIZATION **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **         H(I,K)*U(K,J)=U(I,J)*E(J)                            **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      REAL(8)               :: WORK(3*N)
      INTEGER(4)            :: INFO
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: DEV
      REAL(8)  ,ALLOCATABLE :: EMAT(:,:)
      INTEGER(4)            :: I
!     ******************************************************************
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
      U=0.5D0*(H+TRANSPOSE(H))
      CALL DSYEV('V','U',N,U,N,E,WORK,3*N,INFO )
!
!     ==================================================================
!     == CHECK ERROR CODE                                             ==
!     ==================================================================
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_LAPACK_DSYEV')
        ELSE
          CALL ERROR$MSG('THE ALGORITHM FAILED TO CONVERGE')
          CALL ERROR$MSG('OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('TRIDIAGONAL FORM DID NOT CONVERGE TO ZERO')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_DSYEV')
        END IF
      END IF
!
!     ==================================================================
!     == TEST                                                         ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N))
!       == TEST EIGENVALUE EQUATION ====================================
        EMAT(:,:)=0.D0
        DO I=1,N
          EMAT(I,I)=E(I)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(U,EMAT)))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_DSYEV')
        END IF
!       == TEST ORTHONORMALITY OF EIGENVECTORS =========================
        EMAT=MATMUL(TRANSPOSE(U),U)
        DO I=1,N
          EMAT(I,I)=EMAT(I,I)-1.D0
        ENDDO
        DEV=MAXVAL(ABS(EMAT))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('ORTHONORMALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_DSYEV')
        END IF
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_DSYEV
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_ZHEEV(N,H,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE HERMITEAN, SQUARE MATRIX H                 **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **      CONJG(U(K,I))*H(K,L)*U(L,J)=DELTA(I,J)*E(I)             **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL: ZHPEV :  MATRIX DIAGONALIZATION  P727               **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      COMPLEX(8),INTENT(OUT):: U(N,N)
      INTEGER(4),PARAMETER  :: LWMAX=4000
      COMPLEX(8)            :: CWORK(LWMAX)
      COMPLEX(8),ALLOCATABLE:: CWORK1(:)
      REAL(8)               :: RWORK(3*N-2)
      INTEGER(4)            :: I
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      COMPLEX(8),ALLOCATABLE:: RES(:,:)
      INTEGER(4)            :: LWORK
      INTEGER(4)            :: INFO
      EXTERNAL ZHEEV
!     ******************************************************************
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
!!$      IF(LWMAX.LT.2*N) THEN
!!$        CALL ERROR$MSG('HARDWIRED LIMIT LWMAX SMALLER THAN 2*N')
!!$        CALL ERROR$I4VAL('LWMAX',LWMAX)
!!$        CALL ERROR$I4VAL('2*N',2*N)
!!$        CALL ERROR$STOP('LIB_LAPACK_ZHEEV')
!!$      END IF
      U=0.5D0*(H+TRANSPOSE(CONJG(H)))
      LWORK=-1
      CALL ZHEEV('V','L',N,U,N,E,CWORK,LWORK,RWORK,INFO) !LAPACK
      LWORK=INT(CWORK(1))
      IF(LWORK.LT.LWMAX) THEN
!       LWORK=MIN(LWMAX,INT(CWORK(1)))
        CALL ZHEEV('V','L',N,U,N,E,CWORK,LWORK,RWORK,INFO) !LAPACK
      ELSE
        ALLOCATE(CWORK1(LWORK))
        CALL ZHEEV('V','L',N,U,N,E,CWORK1,LWORK,RWORK,INFO) !LAPACK
        DEALLOCATE(CWORK1)
      END IF
!
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT TO ZHEEV HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZHEEV')
        ELSE
          CALL ERROR$MSG('THE ALGORITHM FAILED TO CONVERGE')
          CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
          CALL ERROR$MSG('TRIDIAGONAL FORM DID NOT CONVERGE TO ZERO.')
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_LAPACK_ZHEEV')
        END IF
      END IF
!
!     ==================================================================
!     ====  OPTIONAL TEST                                             ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(RES(N,N))
        RES=0.5D0*(H+TRANSPOSE(CONJG(H)))
        RES=MATMUL(RES,U)
        DO I=1,N
          RES(:,I)=RES(:,I)-U(:,I)*E(I)
        ENDDO
        IF(MAXVAL(ABS(RES)).GT.1.D-10) THEN
          CALL ERROR$MSG('DIAGONALIZATION TEST FAILED')
          CALL ERROR$R8VAL('DEV ',MAXVAL(ABS(RES)))
          CALL ERROR$STOP('LIB_LAPACK_ZHEEV')
        END IF
        RES=MATMUL(CONJG(TRANSPOSE(U)),U)
        DO I=1,N
          RES(I,I)=RES(I,I)-CMPLX(1.D0,0.D0,KIND=8)
        ENDDO
        IF(MAXVAL(ABS(RES)).GT.1.D-10) THEN
          CALL ERROR$MSG('ORTHONORMALITY  TEST FAILED')
          CALL ERROR$R8VAL('DEV ',MAXVAL(ABS(RES)))
          CALL ERROR$STOP('LIB_LAPACK_ZHEEV')
        END IF
        DEALLOCATE(RES)
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_LAPACK_DSYGV(N,H,S,E,U)
!     ******************************************************************
!     **                                                              **
!     **  DIAGONALIZES THE REAL, SQUARE MATRIX H AFTER SYMMETRIZATION **
!     **  AND RETURNS EIGENVALUES, AND EIGENVECTORS                   **
!     **                                                              **
!     **         U(K,I)*H(K,L)*U(L,J)=DELTA(I,J)*E(I)                 **
!     **                                                              **
!     **  DEPENDENCIES:                                               **
!     **    ESSL DSPEV  MATRIX DIAGONALIZATION P727                   **
!     **                                                              **
!     **  REMARKS:                                                    **
!     **   1) THE EIGENVECTORS ARE REAL BECAUSE IN CASE THEY ARE      **
!     **      COMPLEX REAL AND IMAGINARY PART ARE DEGENERATE          **
!     **      CAN THUS CAN ACT AS EIGENVECTORS THEMSELVES             **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: H(N,N)
      REAL(8)   ,INTENT(IN) :: S(N,N)
      REAL(8)   ,INTENT(OUT):: E(N)
      REAL(8)   ,INTENT(OUT):: U(N,N)
      REAL(8)               :: B(N,N)      
      REAL(8)               :: WORK(3*N)
      INTEGER(4)            :: INFO
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)   ,ALLOCATABLE:: EMAT(:,:)
      REAL(8)               :: DEV
      INTEGER(4)            :: I
!     ******************************************************************
!
!     ==================================================================
!     == DIAGONALIZE                                                  ==
!     ==================================================================
      U=H
      B=S
      CALL DSYGV(1,'V','U',N,U,N,B,N,E,WORK,3*N,INFO )
!
!     ==================================================================
!     == CHECK ERROR CODE                                             ==
!     ==================================================================
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP(' LIB_LAPACK_DSYGV')
        ELSE
          IF(INFO.LE.N) THEN
            CALL ERROR$MSG('DSYEV FAILED TO CONVERGE')
            CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE')
            CALL ERROR$MSG('ITRIDIAGONAL FORM DID NOT CONVERGE TO ZERO')
            CALL ERROR$I4VAL('I',INFO)
            CALL ERROR$STOP(' LIB_LAPACK_DSYGV')
          ELSE
            CALL ERROR$MSG('LEADNG MINOR OF ORDER I OF B IS NOT POSTVE DEFNITE')
            CALL ERROR$MSG('THE FACTORIZATION OF B COULD NOT BE COMPLETED')
            CALL ERROR$MSG('NO EIGENVALUES OR EIGENVECTORS WERE COMPUTED')
            CALL ERROR$I4VAL('I',INFO-N)
            CALL ERROR$STOP(' LIB_LAPACK_DSYGV')
          END IF
        END IF
      END IF
!
!     ==================================================================
!     == TEST                                                         ==
!     ==================================================================
      IF(TTEST) THEN
        ALLOCATE(EMAT(N,N))
        EMAT=0.D0
        DO I=1,N
          EMAT(I,I)=E(I)
        ENDDO
        DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(S,MATMUL(U,EMAT))))
        IF(DEV.GT.1.D-5) THEN
          CALL ERROR$MSG('EIGENVALUE EQUATION TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP(' LIB_LAPACK_DSYGV')
        END IF
        DEALLOCATE(EMAT)
      END IF
      RETURN
      END SUBROUTINE LIB_LAPACK_DSYGV
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_LAPACK_ZHEGV(N,H,S,E,VEC)
!     **************************************************************************
!     **                                                                      **
!     ** SOLVES THE GENERALIZED, COMPLEX EIGENVALUE PROBLEM                   **
!     **      [H(:,:)-E(I)*S(:,:)]*VEC(:,I)=0                                 **
!     **                                                                      **
!     ** REMARK: H AND S MUST BE HERMITEANC                                   **
!     **         S MUST BE POSITIVE DEFINITE                                  **
!     **         EIGENVECTORS ARE ORTHONORMAL IN THE SENSE                    **
!     **             MATMUL(TRANSPOSE(VEC),MATMUL(S,VEC))=IDENTITY            **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: H(N,N)    ! HAMITON MATRIX
      COMPLEX(8),INTENT(IN) :: S(N,N)    ! OVERLAP MATRIX
      REAL(8)   ,INTENT(OUT):: E(N)      ! EIGENVALUES
      COMPLEX(8),INTENT(OUT):: VEC(N,N)  ! EIGENVECTORS
      INTEGER               :: LDWORK
      COMPLEX(8),ALLOCATABLE:: WORK(:)
      COMPLEX(8)            :: S1(N,N)
      REAL(8)               :: RWORK(3*N-2)
      INTEGER               :: INFO
      LOGICAL   ,PARAMETER  :: TTEST=.FALSE.
      REAL(8)               :: DEV
      INTEGER               :: I
!     **************************************************************************
!
!     ==========================================================================
!     == TEST IF INPUT MATRICES ARE HERMITEAN                                 ==
!     ==========================================================================
      IF(TTEST) THEN
        DEV=SUM(ABS(H-TRANSPOSE(CONJG(H))))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('HAMILTON MATRIX NOT HERMITEAN')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
        END IF
        DEV=SUM(ABS(S-TRANSPOSE(CONJG(S))))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('OVERLAP MATRIX NOT HERMITEAN')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
        END IF
      END IF
!
!     ==========================================================================
!     == CALL LAPACK ROUTINE                                                 ==
!     ==========================================================================
      ! LAPACK ROUTINE OVERWRITES HAMILTONIAN WITH EIGENVECTORS
      VEC=0.5D0*(H+TRANSPOSE(CONJG(H)))  
      S1=S
      !DO WORKSPACE QUERY
      ALLOCATE(WORK(1))
      LDWORK=-1
      CALL ZHEGV(1,'V','U',N,VEC,N,S1,N,E,WORK,LDWORK,RWORK,INFO)
      LDWORK=INT(WORK(1))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LDWORK)) 
      CALL ZHEGV(1,'V','U',N,VEC,N,S1,N,E,WORK,LDWORK,RWORK,INFO)
!
!     ==========================================================================
!     == ERROR MESSAGES                                                       ==
!     ==========================================================================
      IF(INFO.NE.0) THEN
        CALL ERROR$MSG('FAILURE SOLVING THE COMPLEX GENERALIZED')
        CALL ERROR$MSG('EIGENVALUE PROBLEM A*X=(LAMBDA)*B*X USING ZHEGV:')
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('ITH ARGUMENT OF ZHGEV HAS ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
        ELSE IF(INFO.GT.N) THEN
          CALL ERROR$MSG('THE LEADING MINOR OF ORDER I OF B IS NOT POSITIVE')
          CALL ERROR$MSG('DEFINITE. THE FACTORIZATION OF B COULD NOT BE')
          CALL ERROR$MSG('COMPLETED AND NO EIGENVALUES OR EIGENVECTORS ')
          CALL ERROR$MSG('WERE COMPUTED.')
          CALL ERROR$I4VAL('I',INFO-N)
        ELSE   !(0.LT.INFO.LE.N)
          CALL ERROR$MSG('ZHEEV FAILED  TO  CONVERGE')
          CALL ERROR$MSG('I OFF-DIAGONAL ELEMENTS OF AN INTERMEDIATE') 
          CALL ERROR$MSG('TRIDIAGONAL FORM DID NOT CONVERGE TO ZERO')
          CALL ERROR$I4VAL('I',INFO)
        END IF
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
      END IF
!
!     ========================================================================
!     == TEST RESULT OF THE ROUTINE                                         ==
!     ========================================================================
      IF(TTEST) THEN
        S1=MATMUL(TRANSPOSE(CONJG(VEC)),MATMUL(S,VEC))
        DEV=0.D0
        DO I=1,N
          S1(I,I)=S1(I,I)-(1.D0,0.D0)
          DEV=MAX(DEV,MAXVAL(ABS(MATMUL(H-E(I)*S,VEC(:,I)))))
        ENDDO
        IF(DEV.GT.1.D-6) THEN
          CALL ERROR$MSG('GENERAL EIGENVALUE TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
        END IF
        DEV=SUM(ABS(S1))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('EIGENSTATES NOT ORTHONORMAL')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGV')
        END IF
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_LAPACK_ZHEGVD(N,H,S,E,VEC)
!     *********************************************************************
!     **                                                                 **
!     ** SOLVES THE GENERALIZED, REAL NON-SYMMETRIC EIGENVALUE PROBLEM   **
!     **      [H(:,:)-E(I)*S(:,:)]*VEC(:,I)=0                            **
!     **                                                                 **
!     ** REMARK: H AND S MUST BE HERMITEANC                              **
!     **         S MUST BE POSITIVE DEFINITE                             **
!     **         EIGENVECTORS ARE ORTHONORMAL IN THE SENSE               **
!     **             MATMUL(TRANSPOSE(VEC),MATMUL(S,VEC))=IDENTITY       **
!     **                                                                 **
!     *********************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)        :: N
      COMPLEX(8),INTENT(IN)        :: H(N,N)    ! HAMITON MATRIX
      COMPLEX(8),INTENT(IN)        :: S(N,N)    ! OVERLAP MATRIX
      REAL(8)   ,INTENT(OUT)       :: E(N)      ! EIGENVALUES
      COMPLEX(8),INTENT(OUT)       :: VEC(N,N)  ! EIGENVECTORS
      COMPLEX(8)                   :: S1(N,N)
      LOGICAL   ,PARAMETER         :: TTEST=.FALSE.
      REAL(8)                      :: DEV
      INTEGER                      :: I
      INTEGER(4)                   :: LWORK,LRWORK,LIWORK,INFO
      COMPLEX(8),ALLOCATABLE       :: WORK(:)
      REAL(8)   ,ALLOCATABLE       :: RWORK(:)
      INTEGER(4),ALLOCATABLE       :: IWORK(:)
      CHARACTER(1)                 :: JOBZ
!     *********************************************************************
!
!     ========================================================================
!     == TEST IF INPUT MATRICES ARE SYMMETRIC                               ==
!     ========================================================================
      IF(TTEST) THEN
        DEV=SUM(ABS(H-TRANSPOSE(CONJG(H))))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('HAMILTON MATRIX NOT HERMITEAN')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGVD')
        END IF
        DEV=SUM(ABS(S-TRANSPOSE(CONJG(S))))
        IF(DEV.GT.1.D-8) THEN
          CALL ERROR$MSG('OVERLAP MATRIX NOT HERMITEAN')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGVD')
        END IF
      END IF
!
!     ========================================================================
!     == CALL LAPACK ROUTINE                                                ==
!     ========================================================================
      JOBZ='V' !COMPUTE EIGENVALUES AND EIGENVECTORS

      ! LAPACK ROUTINE OVERWRITES HAMILTONIAN WITH EIGENVECTORS
      !VEC=0.5D0*(H+TRANSPOSE(CONJG(H)))  
      VEC=H
      S1=S
      !DO WORKSPACE QUERY
      LWORK=-1
      LRWORK=-1
      LIWORK=-1
      ALLOCATE(WORK(1))
      ALLOCATE(RWORK(1))
      ALLOCATE(IWORK(1))
      
      CALL ZHEGVD(1,JOBZ,'U',N,VEC,N,S1,N,E,WORK,LWORK,RWORK,&
     &        LRWORK,IWORK,LIWORK,INFO)
      IF(INFO.NE.0)THEN
        CALL ERROR$MSG('ZHEGVD WORKSPACE QUERY FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('LAPACK_ZHEGVD')
      ENDIF
      LWORK=INT(WORK(1))
      LRWORK=INT(RWORK(1))
      LIWORK=INT(IWORK(1))
      DEALLOCATE(WORK)
      DEALLOCATE(RWORK)
      DEALLOCATE(IWORK)
      ALLOCATE(WORK(LWORK)) 
      ALLOCATE(RWORK(LRWORK)) 
      ALLOCATE(IWORK(LIWORK)) 
      CALL ZHEGVD(1,JOBZ,'U',N,VEC,N,S1,N,E,WORK,LWORK,RWORK,&
     &               LRWORK,IWORK,LIWORK,INFO)
      DEALLOCATE(WORK)
      DEALLOCATE(RWORK)
      DEALLOCATE(IWORK)

      IF(INFO.NE.0)THEN
        CALL ERROR$MSG('ZHEGVD FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('LAPACK_ZHEGVD')
      ENDIF
      

      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('ITH ARGUMENT OF ZHGEVD HAS ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB_LAPACK_ZHEGVD')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$MSG('FAILED')
        CALL ERROR$I4VAL('INFO',INFO)
        CALL ERROR$STOP('LIB_LAPACK_ZHEGVD')
      END IF
!
!     ========================================================================
!     == TEST RESULT OF THE ROUTINE                                         ==
!     ========================================================================
      IF(TTEST) THEN
        S1=MATMUL(TRANSPOSE(CONJG(VEC)),MATMUL(S,VEC))
        DEV=0.D0
        DO I=1,N
          S1(I,I)=S1(I,I)-(1.D0,0.D0)
          DEV=MAX(DEV,MAXVAL(ABS(MATMUL(H-E(I)*S,VEC(:,I)))))
        ENDDO
        IF(DEV.GT.1.D-6) THEN
          CALL ERROR$MSG('GENERAL EIGENVALUE TEST FAILED')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGVD')
        END IF
        DEV=SUM(ABS(S1))
        IF(DEV.GT.1.D-7) THEN
          CALL ERROR$MSG('EIGENSTATES NOT ORTHONORMAL')
          CALL ERROR$R8VAL('DEV',DEV)
          CALL ERROR$STOP('LIB_LAPACK_ZHEGVD')
        END IF
      END IF
      RETURN
      END
!!$#ENDIF
! 
!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****              INTERFACES TO THE BLAS SUBROUTINES                      *****
!****                                                                      *****
!****                                                                      *****
!*******************************************************************************
!*******************************************************************************
!
!DGEMUL: MATRIX MULTIPLICATION P441
!     ..................................................................
      SUBROUTINE LIB$MATMULR8(N,M,L,A,B,C)
!     ******************************************************************
!     **  MATRIX MULTPLICATION   A*B=C                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: A(N,M)
      REAL(8)   ,INTENT(IN) :: B(M,L)
      REAL(8)   ,INTENT(OUT):: C(N,L)
      REAL(8)     ,PARAMETER:: ONE=1.D0
      REAL(8)     ,PARAMETER:: ZERO=0.D0
!     ******************************************************************
      CALL DGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
      RETURN
      END
!
!ZGEMUL: MATRIX MULTIPLICATION P441
!     . ................................................................
      SUBROUTINE LIB$MATMULC8(N,M,L,A,B,C)
!     ******************************************************************
!     **  MATRIX MULTPLICATION   A*B=C                                **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: M
      INTEGER(4),INTENT(IN) :: L
      COMPLEX(8),INTENT(IN) :: A(N,M)
      COMPLEX(8),INTENT(IN) :: B(M,L)
      COMPLEX(8),INTENT(OUT):: C(N,L)
      CHARACTER(8),PARAMETER:: LIB='ESSL'
      COMPLEX(8)  ,PARAMETER:: ONE=(1.D0,0.D0)
      COMPLEX(8)  ,PARAMETER:: ZERO=(0.D0,0.D0)
!     ******************************************************************
      C(:,:)=(0.D0,0.D0)
      CALL ZGEMM('N','N',N,L,M,ONE,A,N,B,M,ZERO,C,N)
      RETURN
      END
!
!ZGEMM(TRANSA,TRANSB,L,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC):
!C=BETA*C+ALPHA*A*B P456
!     ..................................................................
      SUBROUTINE LIB$ADDPRODUCTC8(TID,N,M,L,A,B,C)
!     ******************************************************************
!     **  ADD THE PRODUCT OF TWO MATRICES                             **
!     **  C=C+MATMUL(A,B)                                             **
!     **  FOR TID=.TRUE., A AND C MAY USE IDENTICAL STORAGE           **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN)   :: TID
      INTEGER(4),INTENT(IN)   :: N
      INTEGER(4),INTENT(IN)   :: M
      INTEGER(4),INTENT(IN)   :: L
      COMPLEX(8),INTENT(IN)   :: A(N,M)
      COMPLEX(8),INTENT(IN)   :: B(M,L)
      COMPLEX(8),INTENT(INOUT):: C(N,L)
      COMPLEX(8),ALLOCATABLE  :: WORK(:,:)
!     ******************************************************************
      IF(TID) THEN
        ALLOCATE(WORK(N,M))
        WORK=A
!       ==  C=C+MATMUL(WORK,B) 
        CALL ZGEMM('N','N',N,L,M,(1.D0,0.D0),WORK,N,B,M,(1.D0,0.D0),C,N)
        DEALLOCATE(WORK)
      ELSE
!       ==  C=C+MATMUL(A,B) 
        CALL ZGEMM('N','N',N,L,M,(1.D0,0.D0),A,N,B,M,(1.D0,0.D0),C,N)
      END IF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$DYADSUMR8(LEN1,LEN2,N,PSI1,PSI2,OPERATOR)
!     ******************************************************************
!     **   OPERATOR = SUM_I |PSI1(I)><PSI2(I)|                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN1
      INTEGER(4),INTENT(IN) :: LEN2
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: PSI1(LEN1,N)
      REAL(8)   ,INTENT(IN) :: PSI2(LEN2,N)
      REAL(8)   ,INTENT(OUT):: OPERATOR(LEN1,LEN2)
!     ******************************************************************
!     == OPERATOR=MATMUL(PSI1,TRANSPOSE(PSI2))
      CALL DGEMM('N','T',LEN1,LEN2,N,1.D0,PSI1,LEN1,PSI2,LEN2,0.D0 &
     &          ,OPERATOR,LEN1)
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$DYADSUMC8(LEN1,LEN2,N,PSI1,PSI2,OPERATOR)
!     ******************************************************************
!     **   OPERATOR = SUM_I |PSI1(I)><PSI2(I)|                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LEN1
      INTEGER(4),INTENT(IN) :: LEN2
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: PSI1(LEN1,N)
      COMPLEX(8),INTENT(IN) :: PSI2(LEN2,N)
      COMPLEX(8),INTENT(OUT):: OPERATOR(LEN1,LEN2)
!     ******************************************************************
!     == OPERATOR=MATMUL(PSI1,TRANSPOSE(PSI2))
      OPERATOR(:,:)=(0.D0,0.D0)
      CALL ZGEMM('N','C',LEN1,LEN2,N,(1.D0,0.D0) &
     &          ,PSI1(:,:),LEN1,PSI2(:,:),LEN2,(0.D0,0.D0),OPERATOR,LEN1)
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$SCALARPRODUCTR8(TID,LEN,N1,PSI1,N2,PSI2,OVERLAP)
!     ******************************************************************
!     ** PERFORMS THE SCALAR PRODUCTS OF TWO ARRAYS OF COMPLEX        **
!     ** STATE VECTORS                                                **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TID
      INTEGER(4),INTENT(IN) :: LEN
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      REAL(8)   ,INTENT(IN) :: PSI1(LEN,N1)
      REAL(8)   ,INTENT(IN) :: PSI2(LEN,N2)
      REAL(8)   ,INTENT(OUT):: OVERLAP(N1,N2)
      INTEGER(4)            :: I,J
!     ******************************************************************
      IF(TID.AND.N1.NE.N2) THEN
        CALL ERROR$MSG('PSI2 AND PSI1 DIFFER FOR TID=.TRUE.')
        CALL ERROR$STOP('LIB$SCALARPRODUCTR8')
      END IF
      IF(TID) THEN
!       ==  OVERLAP(I,J) = 0.D0*OVERLAP+1.D0*SUM_K:PSI1(K,I)*PSI1(K,J) =
        CALL DSYRK('U','T',N1,LEN,1.D0,PSI1,LEN,0.D0,OVERLAP,N1)
        DO I=1,N1
          DO J=I+1,N2
            OVERLAP(J,I)=OVERLAP(I,J)
          ENDDO
        ENDDO
      ELSE
!
      CALL DGEMM('T','N',N1,N2,LEN,1.D0,PSI1(:,:),LEN,PSI2(:,:),LEN &
     &             ,0.D0,OVERLAP,N1)
      END IF
      RETURN
      END
!
!ZHERK  C=BETA*C+ ALPHA*A^TA P477
!     ..................................................................
      SUBROUTINE LIB$SCALARPRODUCTC8(TID,LEN,N1,PSI1,N2,PSI2,OVERLAP)
!     ******************************************************************
!     ** PERFORMS THE SCALAR PRODUCTS OF TWO ARRAYS OF COMPLEX        **
!     ** STATE VECTORS                                                **
!     ******************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TID
      INTEGER(4),INTENT(IN) :: LEN
      INTEGER(4),INTENT(IN) :: N1
      INTEGER(4),INTENT(IN) :: N2
      COMPLEX(8),INTENT(IN) :: PSI1(LEN,N1)
      COMPLEX(8),INTENT(IN) :: PSI2(LEN,N2)
      COMPLEX(8),INTENT(OUT):: OVERLAP(N1,N2)
      INTEGER(4)            :: I,J
!     ******************************************************************
      IF(TID.AND.N1.NE.N2) THEN
        CALL ERROR$MSG('PSI2 AND PSI1 DIFFER FOR TID=.TRUE.')
        CALL ERROR$STOP('LIB$SCALARPRODUCTC8')
      END IF
      OVERLAP(:,:)=(0.D0,0.D0)
      IF(TID) THEN
!       == ATTENTION: SCALAR FACTORS ARE SUPPOSED TO BE REAL AS THEY ARE
        CALL ZHERK('U','C',N1,LEN,1.D0,PSI1,LEN,0.D0,OVERLAP,N1)
        DO I=1,N1
          DO J=I+1,N2
            OVERLAP(J,I)=CONJG(OVERLAP(I,J))
          ENDDO
        ENDDO
      ELSE

      CALL ZGEMM('C','N',N1,N2,LEN,(1.D0,0.D0),PSI1,LEN,PSI2,LEN &
     &            ,(0.D0,0.D0),OVERLAP,N1)
      END IF
      RETURN
      END
!
!ZAXPY(N,ALPHA,X,INCX,Y,INCY)  Y=Y+ALPHA*X P267
!     ..................................................................
      SUBROUTINE LIB$VECTORADDC8(N,X,FAC,Y)
!     ******************************************************************
!     **  ADD TWO VECTORS                                             **
!     **  X=X+FAC*Y  WHERE X,Y ARE VECTORS AND FAC IS A SCALAR        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      COMPLEX(8),INTENT(INOUT) :: X(N)
      COMPLEX(8),INTENT(IN)    :: FAC
      COMPLEX(8),INTENT(IN)    :: Y(N)
!     ******************************************************************
      CALL ZAXPY(N,FAC,X,1,Y,1)
      RETURN
      END
!
!DAXPY(N,ALPHA,X,INCX,Y,INCY)  Y=Y+ALPHA*X P267
!DVES(N,X,INCX,Y,INCY,X,INCZ)  Z=X-Y       P313
!     ..................................................................
      SUBROUTINE LIB$VECTORADDR8(N,X,FAC,Y)
!     ******************************************************************
!     **  ADD TWO VECTORS                                             **
!     **  X=X+FAC*Y  WHERE X,Y ARE VECTORS AND FAC IS A SCALAR        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N
      REAL(8)   ,INTENT(INOUT) :: X(N)
      REAL(8)   ,INTENT(IN)    :: FAC
      REAL(8)   ,INTENT(IN)    :: Y(N)
!     ******************************************************************
      CALL DAXPY(N,FAC,X,1,Y,1)
      RETURN
      END
!
!ISORT  SORTS ELEMENTS IN A SEQUENCE P904
!     .................................................................
      SUBROUTINE LIB$SORTI4(N,X)
!     ******************************************************************
!     **  SORTS THE ARRAY X IN ASCENDING ORDER                        **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N
      INTEGER(4),INTENT(INOUT):: X(N)
      INTEGER(4)              :: IHELP
      INTEGER(4)              :: I,J
!     ******************************************************************
#IF DEFINED(CPPVAR_BLAS_ESSL)
      CALL ISORT(X,1,N)
#ELSE
      DO I=1,N-1
        DO J=I+1,N
          IF(X(I).GT.X(J)) THEN
            IHELP=X(I)
            X(I) =X(J)
            X(J) =IHELP
          ENDIF
        ENDDO
      ENDDO
#ENDIF
      RETURN
      END

!
!*******************************************************************************
!*******************************************************************************
!****                                                                       ****
!****  EXTERNAL INTERFACES FOR FOURIER TRANSFORM CALLS                      ****
!****  THESE ROUTINES FORK INTO THE LIBRARY SPECIFIC DRIVER ROUTINES        ****
!****                                                                       ****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$FFTADJUSTGRD(NR)
!     **************************************************************************
!     **  THIS ROUTINE RETURNS THE ALLOWED FOURIER TRANSFORM LENGTH           **
!     **  THAT IS EQUAL OR LARGER THAN THE LENGTH SUPPLIED, BUT               **
!     **  BUT OTHERWISE AS SMALL AS POSSIBLE.                                 **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT):: NR
      INTEGER(4),PARAMETER  :: MAXI=300
      LOGICAL(4),SAVE       :: TINIT=.TRUE.
      INTEGER(4),SAVE       :: COUNT
      INTEGER(4),SAVE       :: IFR(MAXI)
      INTEGER(4)            :: H,I,J,K,M
      INTEGER(4)            :: ISVAR
      REAL(8)               :: SVAR
!     **************************************************************************
      IF (TINIT) THEN
        TINIT=.FALSE.
!       == ALLOWED LENGTHS =====================================================
        COUNT=0
        OUTER: DO H=1,25
          DO I=0,2
            DO J=0,1
              DO K=0,1
                DO M=0,1
                  IF(COUNT.GE.MAXI) EXIT OUTER
                  SVAR = 2.D0**H * 3.D0**I * 5.D0**J * 7.D0**K *11.D0**M
                  IF(SVAR.GT.37748736.D0) CYCLE
                  COUNT=COUNT+1
                  IFR(COUNT)=NINT(SVAR)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO OUTER
!
        DO I=1,COUNT
          DO J=I+1,COUNT
            IF(IFR(I).GT.IFR(J)) THEN
              ISVAR=IFR(I)
              IFR(I)=IFR(J)
              IFR(J)=ISVAR
            END IF
          ENDDO
        ENDDO
      ENDIF
!
      DO I=2,COUNT
        IF(IFR(I).GE.NR) THEN
          NR=IFR(I)
          RETURN
        END IF
      ENDDO
      CALL ERROR$MSG('REQUESTED GRIDPOINT OUT OF RANGE')
      CALL ERROR$I4VAL('NR',NR)
      CALL ERROR$STOP('PLANEWAVE$ADJUSTFFTGRD')
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$FFTC8(DIR,LEN,NFFT,X,Y)                  
!     **************************************************************************
!     **  1-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **  PACKAGES THE ESSL ROUTINE DCFT                                      **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                             **
!     **                                                                      **
!     **  USE FFTW AS STANDARD                                                **
!     **  USE FFTESSL IF ESSL IS INSTALLED                                    **
!     **  USE FFTPACK AS BACKUP IF C-ROUTINES CANNOT BE LINKED OR             **
!     **      FFTW IS NOT AVAILABLE                                           **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR !'GTOR' OR 'RTOG'
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
!     **************************************************************************
#IF DEFINED(CPPVAR_FFT_ESSL)
      CALL LIB_FFTESSL(DIR,LEN,NFFT,X,Y)                  
#ELIF DEFINED(CPPVAR_FFT_ACML)
      CALL LIB_ACML_FFT1DC8(DIR,LEN,NFFT,X,Y)
#ELIF DEFINED(CPPVAR_FFT_PACK)
      CALL LIB_FFTPACK(DIR,LEN,NFFT,X,Y)
#ELSE
      CALL LIB_FFTW3(DIR,LEN,NFFT,X,Y)
#ENDIF

      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$3DFFTC8(DIR,N1,N2,N3,X,Y)
!     **************************************************************************
!     **  3-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **    USES THE 3D FFTW ROUTINES                                         **
!     **                                        CLEMENS FOERST, 2001          **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(4)            :: DIR !'GTOR' OR 'RTOG'
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
!     **************************************************************************
#IF DEFINED(CPPVAR_FFT_ESSL)
      CALL LIB_3DFFT_ESSL(DIR,N1,N2,N3,X,Y)
#ELIF DEFINED(CPPVAR_FFT_ACML)
      CALL LIB_ACML_FFT3DC8(DIR,N1,N2,N3,X,Y)
#ELIF DEFINED(CPPVAR_FFT_PACK)
      CALL LIB_3DFFTPACK(DIR,N1,N2,N3,X,Y)
#ELSE
      CALL LIB_3DFFTW3(DIR,N1,N2,N3,X,Y)
#ENDIF
      RETURN
      END

#IF DEFINED(CPPVAR_FFT_ACML)
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB_ACML_FFT1DC8(DIR,LEN,NFFT,X,Y)                  
!     **************************************************************************
!     **  1-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **  PACKAGES THE ACML ROUTINE ZFFT1M                                    **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                             **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      INTEGER                 :: NFFT_1
      INTEGER                 :: LEN_1
      INTEGER                 :: MODE
      INTEGER                 :: INFO
      COMPLEX(8)  ,ALLOCATABLE,SAVE :: COMM(:)
      INTEGER(4)              ,SAVE :: LENPREV=0
      REAL(8)                 ,SAVE :: SCALEFORWARD
      REAL(8)                 ,SAVE :: SCALEBACKWARD
!     **************************************************************************
      IF(LEN.EQ.0) RETURN
      NFFT_1=NFFT
      LEN_1=LEN
      IF(DIR.EQ.'RTOG') THEN
        MODE=-1
      ELSE IF(DIR.EQ.'GTOR') THEN
        MODE=1
      ELSE
        CALL ERROR$MSG('INVALID VALUE OF VARIABLE "DIR"')
        CALL ERROR$MSG('CAN BE "GTOR" OR "RTOG"')
        CALL ERROR$CHVAL('DIR',DIR)
        CALL ERROR$STOP('LIB_ACML_FFT1DC8')
      END IF
! 
!     ==========================================================================
!     == INITIALIZE PLAN                                                      ==
!     ==========================================================================
      IF(LEN.NE.LENPREV) THEN
        IF(LENPREV.NE.0) DEALLOCATE(COMM)
        ALLOCATE(COMM(3*LEN+100))
        LENPREV=LEN
        Y=X
        CALL ZFFT1M(100,NFFT,LEN,Y,COMM,INFO)
        IF(INFO.LT.0) THEN
          CALL ERROR$MSG('I-TH ARGUMENT OF ZFFT1M HAD AN ILLEGAL VALUE')
          CALL ERROR$I4VAL('I',-INFO)
          CALL ERROR$STOP('LIB_ACML_FFT1DC8')
        ELSE IF(INFO.GT.0) THEN
          CALL ERROR$I4VAL('I',INFO)
          CALL ERROR$STOP('LIB_ACML_FFT1DC8')
        END IF
        SCALEFORWARD=1.D0/SQRT(REAL(LEN,KIND=8))
        SCALEBACKWARD=1.D0/SCALEFORWARD
      END IF
! 
!     ==========================================================================
!     == PERFORM FOURIER TRANSFORM                                            ==
!     ==========================================================================
      Y=X
      CALL ZFFT1M(MODE,NFFT,LEN,Y,COMM,INFO)
      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('I-TH ARGUMENT OF ZFFT1M HAD AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB_ACML_FFT1DC8')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$I4VAL('I',INFO)
        CALL ERROR$STOP('LIB_ACML_FFT1DC8')
      END IF
! 
!     ==========================================================================
!     == SCALE RESULT                                                         ==
!     ==========================================================================
      IF(DIR.EQ.'RTOG') THEN
         Y=Y*SCALEFORWARD
      ELSE 
         Y=Y*SCALEBACKWARD
      END IF
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE LIB_ACML_FFT3DC8(DIR,N1,N2,N3,X,Y)
!     **************************************************************************
!     **  3-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **    USES THE 3D FFTW ROUTINES                                         **
!     **                                        CLEMENS FOERST, 2001          **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      INTEGER                 :: N1_1,N2_1,N3_1
      INTEGER                 :: MODE
      INTEGER                 :: INFO
      COMPLEX(8)  ,ALLOCATABLE,SAVE :: COMM(:)
      INTEGER(4)              ,SAVE :: N1PREV=0
      INTEGER(4)              ,SAVE :: N2PREV=0
      INTEGER(4)              ,SAVE :: N3PREV=0
      REAL(8)                 ,SAVE :: SCALEFORWARD
      REAL(8)                 ,SAVE :: SCALEBACKWARD
!     **************************************************************************
      IF(N1*N2*N3.EQ.0) RETURN
      N1_1=N1
      N2_1=N2
      N3_1=N3
      IF(DIR.EQ.'RTOG') THEN
        MODE=-1
      ELSE IF(DIR.EQ.'GTOR') THEN
        MODE=1
      ELSE
        CALL ERROR$MSG('INVALID VALUE OF VARIABLE "DIR"')
        CALL ERROR$MSG('CAN BE "GTOR" OR "RTOG"')
        CALL ERROR$CHVAL('DIR',DIR)
        CALL ERROR$STOP('LIB_ACML_3DFFTC8')
      END IF
! 
!     ==========================================================================
!     == INITIALIZE PLAN                                                      ==
!     ==========================================================================
      IF(N1.NE.N1PREV.OR.N2.NE.N2PREV.OR.N3.NE.N3PREV) THEN
        IF(N1PREV.NE.0) DEALLOCATE(COMM)
        ALLOCATE(COMM(N1*N2*N3+3*(N1+N2+N3)))
        N1PREV=N1
        N2PREV=N2
        N3PREV=N3
        SCALEFORWARD=1.D0/SQRT(REAL(N1*N2*N3,KIND=8))
        SCALEBACKWARD=1.D0/SCALEFORWARD
      END IF
! 
!     ==========================================================================
!     == PERFORM FOURIER TRANSFORM                                            ==
!     ==========================================================================
      Y=X
      CALL ZFFT3D(MODE,N1_1,N2_1,N3_1,Y,COMM,INFO)

      IF(INFO.LT.0) THEN
        CALL ERROR$MSG('I-TH ARGUMENT OF ZFFT1M HAD AN ILLEGAL VALUE')
        CALL ERROR$I4VAL('I',-INFO)
        CALL ERROR$STOP('LIB_ACML_3DFFTC8')
      ELSE IF(INFO.GT.0) THEN
        CALL ERROR$I4VAL('I',INFO)
        CALL ERROR$STOP('LIB_ACML_3DFFTC8')
      END IF
! 
!     ==========================================================================
!     == SCALE RESULT                                                         ==
!     ==========================================================================
      IF(DIR.EQ.'RTOG') THEN
         Y=Y*SCALEFORWARD
      ELSE 
         Y=Y*SCALEBACKWARD
      END IF
      RETURN
      END
#ENDIF
!
!*******************************************************************************
!**  DRIVER ROUTINES FOR THE FOURIER TRANSFORMS FROM THE FFTW3 PACKAGE        **
!**                                                                           **
!**  NOTE THAT THE CALLS TO FFTW3 ARE NOT COMPATIBLE WITH FFTW2,              **
!**  WHICH IS CALLED HERE AS FFTW                                             **
!**                                                                           **
!**               HTTP://WWW.FFTW.ORG/                                        **
!*******************************************************************************
!     ..........................................................................
      SUBROUTINE LIB_FFTW3(DIR,LEN,NFFT,X,Y)                  
!     **************************************************************************
!     **  1-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                             **
!     **                                                                      **
!     **************************************************************************
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      CHARACTER(4),SAVE       :: DIRSAVE=' '
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      REAL(8)     ,SAVE       :: SCALE
      INTEGER     ,SAVE       :: ISIGN
      COMPLEX(8)              :: XDUMMY(LEN)
      COMPLEX(8)              :: YDUMMY(LEN)
      INTEGER(4),SAVE         :: NP=0
      INTEGER(4),PARAMETER    :: NPX=100 ! #(DIFFERENT FFT PLANS)
      TYPE(C_PTR),SAVE        :: PLANS2(NPX),PLAN
      INTEGER(4),SAVE         :: PLANS1(NPX)
      LOGICAL                 :: DEF
      INTEGER(4)              :: I
      INCLUDE 'FFTW3.F03'
!     **************************************************************************
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE) THEN
        IF (DIR.EQ.'GTOR') THEN
          ISIGN=1
        ELSE IF (DIR.EQ.'RTOG') THEN
          ISIGN=-1
        ELSE
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('1D-FFTW')
        END IF
!
!       == FIND PLAN IN THE LIST ===============================================
        DEF=.FALSE.
        DO I=1,NP
          IF((LEN*ISIGN).EQ.PLANS1(I)) THEN
            DEF=.TRUE.
            PLAN=PLANS2(I)
            EXIT
          END IF
        END DO
!
!       == CREATE NEW PLAN IF NOT IN THE LIST ==================================
        IF(.NOT.DEF) THEN
          WRITE(*,*) 'LIB_FFTW: CREATE PLAN FOR: ',TRIM(DIR),ISIGN,LEN,NP
          NP=NP+1
          IF(NP.GE.NPX) NP=NPX ! ALLOW ONLY NPX PLANS
          IF(DIR.EQ.'RTOG') THEN
            PLANS2(NP) = FFTW_PLAN_DFT_1D(LEN,XDUMMY,YDUMMY,FFTW_FORWARD &
     &                ,IOR(FFTW_DESTROY_INPUT,IOR(FFTW_MEASURE,FFTW_UNALIGNED)))
          ELSE IF (DIR.EQ.'GTOR') THEN
            PLANS2(NP) = FFTW_PLAN_DFT_1D(LEN,XDUMMY,YDUMMY,FFTW_BACKWARD &
     &                ,IOR(FFTW_DESTROY_INPUT,IOR(FFTW_MEASURE,FFTW_UNALIGNED)))
          ELSE
            CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
            CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
            CALL ERROR$CHVAL('DIR',TRIM(DIR))
            CALL ERROR$STOP('LIB_FFTW3')
          END IF 

          PLANS1(NP)=ISIGN*LEN
          PLAN=PLANS2(NP)
        END IF
        LENSAVE=LEN
        DIRSAVE=DIR
      END IF
!
!     ==========================================================================
!     ==  NOW PERFORM FFT                                                     ==
!     ==========================================================================
      DO I=1,NFFT
        XDUMMY=X(:,I)
        CALL FFTW_EXECUTE_DFT(PLAN, XDUMMY, Y(:,I))
      ENDDO
!
!     ==========================================================================
!     ==  SCALE RESULT                                                        ==
!     ==========================================================================
      IF (DIR.EQ.'RTOG') THEN
        SCALE=1.D0/REAL(LEN,KIND=8)
        Y(:,:)=Y(:,:)*SCALE
      END IF
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE LIB_3DFFTW3(DIR,N1,N2,N3,X,Y)
!     **************************************************************************
!     **  3-D FFT                                                             **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)                    **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)                    **
!     **                                                                      **
!     **    USES THE 3D FFTW ROUTINES                                         **
!     **                                        CLEMENS FOERST, 2001          **
!     **************************************************************************
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      TYPE(C_PTR) :: PLAN
      REAL(8)     ,SAVE       :: SCALE
      INCLUDE 'FFTW3.F03' ! FILENAME MADE LOWERCASE BY F90PP 
!     **************************************************************************
      PRINT*,"3DFFTW3",N1,N2,N3
      IF(DIR.EQ.'RTOG') THEN
        PLAN = FFTW_PLAN_DFT_3D(N3,N2,N1,X,Y,FFTW_FORWARD,FFTW_ESTIMATE)
      ELSE IF (DIR.EQ.'GTOR') THEN
        PLAN = FFTW_PLAN_DFT_3D(N3,N2,N1,X,Y,FFTW_BACKWARD,FFTW_ESTIMATE)
      ELSE
        CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
        CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
        CALL ERROR$CHVAL('DIR',TRIM(DIR))
        CALL ERROR$STOP('LIB_3DFFTW3')
      END IF  
!
!     ==========================================================================
!     ==  EXECUTE FOURIER TRANSFORM                                           ==
!     ==========================================================================
      CALL FFTW_EXECUTE_DFT(PLAN, X, Y)
      CALL FFTW_DESTROY_PLAN(PLAN)
!
!     ==========================================================================
!     ==  SCALE RESULT                                                        ==
!     ==========================================================================
      IF (DIR.EQ.'RTOG') THEN
        SCALE=1.D0/REAL(N1*N2*N3,KIND=8)
        Y=Y*SCALE
      END IF
      RETURN
      END
!
!*******************************************************************************
!**  DRIVER ROUTINES FOR THE FOURIER TRANSFORMS FROM THE                      **
!**     ENGINEERING AND SCIENTIFIC SUBROUTINE LIBRARY (ESSL)                  **
!**                                                                           **
!**  ESSL IS A COMMERCIAL PACKAGE FROM IBM                                    **
!*******************************************************************************
!DCFT:  1-D FFT P765
#IF DEFINED(CPPVAR_FFT_ESSL)
!     ..................................................................
      SUBROUTINE LIB_FFTESSL(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      INTEGER(4)  ,PARAMETER  :: NAUX1=20000
      INTEGER(4)  ,PARAMETER  :: NAUX2=20000
      REAL(8)     ,SAVE       :: AUX2(NAUX2)
      REAL(8)     ,SAVE       :: AUX1(NAUX1)
      CHARACTER(4),SAVE       :: DIRSAVE=' '
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER(4)  ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)              :: IFFT,I
!     ******************************************************************
!  
!     ==================================================================
!     == INITIALIZATION PHASE                                         ==
!     ==================================================================
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE.OR.NFFT.NE.NFFTSAVE) THEN
        IF(LEN.GT.8192) THEN
          CALL ERROR$MSG('FFT TOO LONG')
          CALL ERROR$I4VAL('LEN',LEN)
          CALL ERROR$STOP('FFT')
        END IF
        IF(DIR.EQ.'GTOR') THEN
          ISIGN=-1
          SCALE=1.D0
        ELSE IF(DIR.EQ.'RTOG') THEN
          ISIGN=1
          SCALE=1.D0/REAL(LEN,KIND=8)
        ELSE 
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('FFT')
        END IF
        CALL DCFT(1,X,1,LEN,Y,1,LEN,LEN,NFFT,ISIGN,SCALE,AUX1,NAUX1,AUX2,NAUX2)
        DIRSAVE=DIR
        LENSAVE=LEN
        NFFTSAVE=NFFT
      END IF
!  
!     ==================================================================
!     == FOURIER TRANSFORM                                            ==
!     ==================================================================
      CALL DCFT(0,X,1,LEN,Y,1,LEN,LEN,NFFT,ISIGN,SCALE,AUX1,NAUX1,AUX2,NAUX2)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_3DFFT_ESSL(DIR,N1,N2,N3,X,Y)
!     ******************************************************************
!     **  3-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **    USES THE 3D FFT ROUTINE OF ESSL DCFT3                     **
!     **                                        PETER BLOECHL, 2001   **
!     ******************************************************************
!     **    NOT TESTED                                                **
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: DIM(3)
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      INTEGER(4)              :: NAUX       
      REAL(8)   ,ALLOCATABLE  :: AUX(:)
      INTEGER(4)              :: ISIGN
      REAL(8)                 :: SCALE
      INTEGER(4)              :: S,PSI,LAMBDA
!     ******************************************************************
      NAUX=60000
      IF(N1.GT.2048) NAUX=NAUX+4.56*N1
      IF(N3.LT.252) THEN
        IF(N2.GE.252) THEN
          S=MIN(64,N1)
          LAMBDA=(2*N2+256)*(S+4.56)
          NAUX=NAUX+LAMBDA
        END IF
      ELSE
        IF(N2.GE.252) THEN
          S=MIN(64,N1*N2)
          PSI=(2*N3+256)*(S+4.56)
          NAUX=NAUX+PSI
        ELSE
          S=MIN(64,N1*N2)
          PSI=(2*N3+256)*(S+4.56)
          S=MIN(64,N1)
          LAMBDA=(2*N2+256)*(S+4.56)
          NAUX=NAUX+MAX(PSI,LAMBDA)
        END IF
      END IF
      IF (DIR.EQ.'GTOR') THEN
        ISIGN=1
        SCALE=1.D0
      ELSE IF (DIR.EQ.'RTOG') THEN
        ISIGN=-1
        SCALE=1.D0/REAL(N1*N2*N3,KIND=8)
      ELSE
        CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
        CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
        CALL ERROR$CHVAL('DIR',TRIM(DIR))
        CALL ERROR$STOP('LIB_3DFFT_ESSL')
      END IF
      ALLOCATE(AUX(NAUX))
      CALL DCFT3(X,N1,N1*N2,Y,N1,N1*N2,N1,N2,N3,ISIGN,SCALE,AUX,NAUX)
      DEALLOCATE(AUX)
      RETURN
      END
#ENDIF
!
!
!*******************************************************************************
!**  DRIVER ROUTINES FOR THE FOURIER TRANSFORMS FROM THE                      **
!**     FFTPACK PACKAGE                                                       **
!**                                                                           **
!**  FFTPACK IS A FORTRAN SUBROUTINE LIBRARY OF FAST FOURIER TRANSFORM        **
!**  DEVELOPED AT THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH.               **
!**                                                                           **
!**   FFTPACK5 IS DISTRIBUTED UNDER THE GNU GENERAL PUBLIC LICENSE            **
!**   FROM HTTP://WWW.CISL.UCAR.EDU/CSS/SOFTWARE/FFTPACK5/INDEX.HTML          **
!**                                                                           **
!*******************************************************************************
#IF DEFINED(CPPVAR_FFT_PACK)
!     ..................................................................
      SUBROUTINE LIB_FFTPACK(DIR,LEN,NFFT,X,Y)                  
!     ******************************************************************
!     **  1-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **  PACKAGES THE ESSL ROUTINE DCFT                              **
!     **  REMARK: X AND Y MAY BE IDENTICAL ARRAYS                     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: DIR
      INTEGER(4)  ,INTENT(IN) :: LEN
      INTEGER(4)  ,INTENT(IN) :: NFFT
      COMPLEX(8)  ,INTENT(IN) :: X(LEN,NFFT)
      COMPLEX(8)  ,INTENT(OUT):: Y(LEN,NFFT)
      INTEGER(4)  ,PARAMETER  :: NAUX2=2000
      REAL(8)     ,SAVE       :: AUX2(NAUX2)
      INTEGER(4)  ,SAVE       :: IAUX(30)
      REAL(8)                 :: SEQUENCE(2*LEN)
      CHARACTER(4),SAVE       :: DIRSAVE=' '
      INTEGER(4)  ,SAVE       :: LENSAVE=0
      INTEGER(4)  ,SAVE       :: NFFTSAVE=0
      INTEGER(4)  ,SAVE       :: ISIGN
      REAL(8)     ,SAVE       :: SCALE
      INTEGER(4)              :: IFFT,I
!     ******************************************************************
!  
!     ==================================================================
!     == INITIALIZATION PHASE                                         ==
!     ==================================================================
      IF(DIR.NE.DIRSAVE.OR.LEN.NE.LENSAVE) THEN
        IF(DIR.EQ.'GTOR') THEN
          ISIGN=-1
          SCALE=1.D0
        ELSE IF(DIR.EQ.'RTOG') THEN
          ISIGN=1
          SCALE=1.D0/REAL(LEN,KIND=8)
        ELSE 
          CALL ERROR$MSG('DIRECTION ID NOT RECOGNIZED')
          CALL ERROR$MSG('DIR MUST BE "GTOR" OR "RTOG"')
          CALL ERROR$CHVAL('DIR',TRIM(DIR))
          CALL ERROR$STOP('LIB_FFTPACK')
        END IF
        DIRSAVE=DIR
        LENSAVE=LEN
        IF(LEN.EQ.1) RETURN
        IF(NAUX2.LT.2*LEN) THEN
          CALL ERROR$MSG('AUXILIARY ARRAY TOO SMALL: INCREASE NAUX2')
          CALL ERROR$I4VAL('NAUX2',NAUX2)
          CALL ERROR$I4VAL('2*LEN',2*LEN)
          CALL ERROR$STOP('LIB_FFTPACK')
        END IF
        CALL CFFTI(LEN,AUX2,IAUX)
      END IF
!  
!     ==================================================================
!     == FOURIER TRANSFORM                                            ==
!     ==================================================================
      IF(LEN.EQ.1) RETURN
      IF(ISIGN.EQ.1) THEN
        DO IFFT=1,NFFT
          DO I=1,LEN
            SEQUENCE(2*I-1)= REAL(X(I,IFFT),KIND=8)
            SEQUENCE(2*I  )=AIMAG(X(I,IFFT))
          ENDDO
          CALL CFFTF(LEN,SEQUENCE,AUX2,IAUX)
          DO I=1,LEN
            Y(I,IFFT)=CMPLX(SEQUENCE(2*I-1),SEQUENCE(2*I),KIND=8)*SCALE
          ENDDO
        ENDDO
      ELSE
        DO IFFT=1,NFFT
          DO I=1,LEN
            SEQUENCE(2*I-1)= REAL(X(I,IFFT),KIND=8)
            SEQUENCE(2*I  )=AIMAG(X(I,IFFT))
          ENDDO
          CALL CFFTB(LEN,SEQUENCE,AUX2,IAUX)
          DO I=1,LEN
            Y(I,IFFT)=CMPLX(SEQUENCE(2*I-1),SEQUENCE(2*I),KIND=8)*SCALE
          ENDDO
        ENDDO
      END IF
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB_3DFFTPACK(DIR,N1,N2,N3,X,Y)
!     ******************************************************************
!     **  3-D FFT                                                     **
!     **    DIR='GTOR' => Y(R)=     SUM_G X(G) EXP( I*G*R)            **
!     **    DIR='RTOG' => Y(G)=1/NR SUM_R X(R) EXP(-I*G*R)            **
!     **                                                              **
!     **    USES THE 3D FFT ROUTINE OF ESSL DCFT3                     **
!     **                                        PETER BLOECHL, 2001   **
!     ******************************************************************
!     **    NOT TESTED                                                **
      IMPLICIT NONE
      CHARACTER(4)            :: DIR
      INTEGER(4)              :: DIM(3)
      INTEGER(4)              :: N1,N2,N3
      COMPLEX(8)              :: X(N1,N2,N3)
      COMPLEX(8)              :: Y(N1,N2,N3)
      COMPLEX(8)              :: WORK1(N1*N2*N3),WORK2(N1*N2*N3)
      INTEGER(4)              :: I,J,K,IND
!     ******************************************************************
      CALL LIB_FFTPACK(DIR,N1,N2*N3,X,WORK2)
      IND=0
      DO K=1,N3
        DO J=1,N2
          DO I=1,N1
            IND=IND+1
            WORK1(J+N2*(K-1+N1*(I-1)))=WORK2(IND)
          ENDDO
        ENDDO
      ENDDO
      CALL LIB_FFTPACK(DIR,N2,N1*N3,WORK1,WORK2)
      IND=0
      DO I=1,N1
        DO K=1,N3
          DO J=1,N2
            IND=IND+1
            WORK1(K+N3*(I-1+N1*(J-1)))=WORK2(IND)
          ENDDO
        ENDDO
      ENDDO
      CALL LIB_FFTPACK(DIR,N3,N1*N2,WORK1,WORK2)
      IND=0
      DO J=1,N2
        DO I=1,N1
          DO K=1,N3
            IND=IND+1
            Y(I,J,K)=WORK2(IND)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
#ENDIF
!
!     .................................................................
      SUBROUTINE LIB$GETUSAGE(ID,VALUE)
!     *****************************************************************
!     **                                                             **
!     **  PROVIDES INFORMATION  ON THE USAGE OF SYSTEM RESOURCES     **
!     **                                                             **
!     **  USES STANDARD C-LIBRARY ROUTINE GETRUSAGE                  **
!     **                                                             **
!     **                                                             **
!     **                                                             **
!     *****************************************************************
      IMPLICIT NONE
      TYPE USG_TYPE 
        SEQUENCE
        INTEGER(4) :: UTIME(2)   ! USER TIME    (SECONDS,MICROSECONDS)
        INTEGER(4) :: STIME(2)   ! SYSTEM TIME  (SECONDS,MICROSECONDS)
        INTEGER(4) :: MAXRSS     ! MAX SIZE [KBYTE] 
        INTEGER(4) :: IXRSS      ! INTEGRAL SHARED MEMORY SIZE [KBYTE*SEC]         
        INTEGER(4) :: IDRSS      ! INTEGRAL UNSHARED DATA [KBYTE*SEC]    
        INTEGER(4) :: ISRSS      ! INTEGRAL UNSHARED STACK [KBYTE*SEC]    
        INTEGER(4) :: MINFLT     ! #(PAGE RECLAIMS)
        INTEGER(4) :: MAJFLT     ! #(PAGE FAULTS)
        INTEGER(4) :: NSWAP      ! #(SWAPPING PROCESS OUT OF MAIN)            
        INTEGER(4) :: INBLOCK    ! #(FILE INPUTS)                             
        INTEGER(4) :: OUBLOCK    ! #(FILE OUTPUTS)                            
        INTEGER(4) :: MSGSND     ! #(IPC MESSAGES SEND)                       
        INTEGER(4) :: MSGRCV     ! #(IPC MESSAGES RECEIVED)                   
        INTEGER(4) :: NSIGNALS   ! #(SIGNALS SENT)                            
        INTEGER(4) :: NVCSW      ! #(VOLUNTARY CONTEXT SWITCHES)
        INTEGER(4) :: NIVCSW     ! #(INVOLUNTARY CONTEXT SWITCHES)
      END TYPE USG_TYPE
      CHARACTER(*),INTENT(IN)  :: ID
      REAL(8)     ,INTENT(OUT) :: VALUE
      TYPE(USG_TYPE)           :: USG
      REAL(8)                  :: KBYTE=2.D0**10
      INTEGER(4)               :: GETRUSAGE
      INTEGER(4)               :: RC
      REAL(8)                  :: CPUTIME
#IF DEFINED(CPPVAR_USAGE_EXIST)
      EXTERNAL GETRUSAGE
#ENDIF
!     ******************************************************************
      USG%UTIME=(/0,0/)
      USG%STIME=(/0,0/)
      USG%MAXRSS=0     
      USG%IXRSS=0      
      USG%IDRSS=0      
      USG%ISRSS=0      
      USG%MINFLT=0     
      USG%MAJFLT=0     
      USG%NSWAP=0      
      USG%INBLOCK=0    
      USG%OUBLOCK=0    
      USG%MSGSND=0     
      USG%MSGRCV=0     
      USG%NSIGNALS=0   
      USG%NVCSW=0      
      USG%NIVCSW=0     
!     ==================================================================
!     ==================================================================
#IF DEFINED(CPPVAR_USAGE_EXIST)
      RC=GETRUSAGE(%VAL(0),USG)    ! C-ROUTINE
      IF(RC.NE.0) THEN
        CALL ERROR$MSG('ERROR CALLING MYGETRUSAGE')
        CALL ERROR$STOP('MEMORY$GET')
      END IF
#ENDIF
!
      CPUTIME=(REAL(USG%UTIME(1)+USG%STIME(1),KIND=8) &
     &        +REAL(USG%UTIME(2)+USG%STIME(2),KIND=8))*1.D-6
      CPUTIME=MAX(CPUTIME,1.D-6)
      IF(ID.EQ.'MAXMEM') THEN
        VALUE=REAL(USG%MAXRSS,KIND=8)*KBYTE
      ELSE IF(ID.EQ.'USRTIME') THEN
        VALUE=(REAL(USG%UTIME(1),KIND=8)+REAL(USG%UTIME(2),KIND=8))*1.D-6
      ELSE IF(ID.EQ.'SYSTIME') THEN
        VALUE=(REAL(USG%STIME(1),KIND=8)+REAL(USG%STIME(2),KIND=8))*1.D-6
      ELSE IF(ID.EQ.'CPUTIME') THEN
        VALUE=CPUTIME
      ELSE IF(ID.EQ.'PAGEFAULTRATE') THEN
        VALUE=REAL(USG%MAJFLT,KIND=8)/CPUTIME
      ELSE IF(ID.EQ.'SWAPRATE') THEN
        VALUE=REAL(USG%NSWAP,KIND=8)/CPUTIME
      ELSE IF(ID.EQ.'SWITCHRATE') THEN
        VALUE=REAL(USG%NVCSW+USG%NIVCSW,KIND=8)/CPUTIME
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$STOP('MEMORY$GET')
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$ERFR8(X,Y)
!     **************************************************************************
!     **  COPYRIGHT(C) 1996 TAKUYA OOURA                                      **
!     **  (EMAIL: OOURA@MMM.T.U-TOKYO.AC.JP).                                 **
!     **  YOU MAY USE, COPY, MODIFY THIS CODE FOR ANY PURPOSE AND             **
!     **  WITHOUT FEE. YOU MAY DISTRIBUTE THIS ORIGINAL PACKAGE.              **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
      REAL(8)            :: W
      REAL(8)            :: T
      INTEGER(4)         :: K,I
      REAL(8)            :: A(0:64)
      REAL(8)            :: B(0:64)
      DATA (A(I), I = 0, 12) / &
     &    0.00000000005958930743D0, -0.00000000113739022964D0, &
     &    0.00000001466005199839D0, -0.00000016350354461960D0, &
     &    0.00000164610044809620D0, -0.00001492559551950604D0, &
     &    0.00012055331122299265D0, -0.00085483269811296660D0, &
     &    0.00522397762482322257D0, -0.02686617064507733420D0, &
     &    0.11283791670954881569D0, -0.37612638903183748117D0, &
     &    1.12837916709551257377D0 / 
      DATA (A(I), I = 13, 25) /                                &
     &    0.00000000002372510631D0, -0.00000000045493253732D0, &
     &    0.00000000590362766598D0, -0.00000006642090827576D0, &
     &    0.00000067595634268133D0, -0.00000621188515924000D0, &
     &    0.00005103883009709690D0, -0.00037015410692956173D0, &
     &    0.00233307631218880978D0, -0.01254988477182192210D0, &
     &    0.05657061146827041994D0, -0.21379664776456006580D0, &
     &    0.84270079294971486929D0 / 
      DATA (A(I), I = 26, 38) /                                &
     &    0.00000000000949905026D0, -0.00000000018310229805D0, &
     &    0.00000000239463074000D0, -0.00000002721444369609D0, &
     &    0.00000028045522331686D0, -0.00000261830022482897D0, &
     &    0.00002195455056768781D0, -0.00016358986921372656D0, &
     &    0.00107052153564110318D0, -0.00608284718113590151D0, &
     &    0.02986978465246258244D0, -0.13055593046562267625D0, &
     &    0.67493323603965504676D0 / 
      DATA (A(I), I = 39, 51) /                                &
     &    0.00000000000382722073D0, -0.00000000007421598602D0, &
     &    0.00000000097930574080D0, -0.00000001126008898854D0, &
     &    0.00000011775134830784D0, -0.00000111992758382650D0, &
     &    0.00000962023443095201D0, -0.00007404402135070773D0, &
     &    0.00050689993654144881D0, -0.00307553051439272889D0, &
     &    0.01668977892553165586D0, -0.08548534594781312114D0, &
     &    0.56909076642393639985D0 / 
      DATA (A(I), I = 52, 64) /                                &
     &    0.00000000000155296588D0, -0.00000000003032205868D0, &
     &    0.00000000040424830707D0, -0.00000000471135111493D0, &
     &    0.00000005011915876293D0, -0.00000048722516178974D0, &
     &    0.00000430683284629395D0, -0.00003445026145385764D0, &
     &    0.00024879276133931664D0, -0.00162940941748079288D0, &
     &    0.00988786373932350462D0, -0.05962426839442303805D0, &
     &    0.49766113250947636708D0 / 
      DATA (B(I), I = 0, 12) /                                 &
     &   -0.00000000029734388465D0,  0.00000000269776334046D0, &
     &   -0.00000000640788827665D0, -0.00000001667820132100D0, & 
     &   -0.00000021854388148686D0,  0.00000266246030457984D0, &
     &    0.00001612722157047886D0, -0.00025616361025506629D0, &
     &    0.00015380842432375365D0,  0.00815533022524927908D0, &
     &   -0.01402283663896319337D0, -0.19746892495383021487D0, & 
     &    0.71511720328842845913D0 / 
      DATA (B(I), I = 13, 25) /                                 &
     &   -0.00000000001951073787D0, -0.00000000032302692214D0,  &
     &    0.00000000522461866919D0,  0.00000000342940918551D0,  &
     &   -0.00000035772874310272D0,  0.00000019999935792654D0,  &
     &    0.00002687044575042908D0, -0.00011843240273775776D0,  &
     &   -0.00080991728956032271D0,  0.00661062970502241174D0,  &
     &    0.00909530922354827295D0, -0.20160072778491013140D0,  &
     &    0.51169696718727644908D0 / 
      DATA (B(I), I = 26, 38) /                                 &
     &    0.00000000003147682272D0, -0.00000000048465972408D0,  &
     &    0.00000000063675740242D0,  0.00000003377623323271D0,  &
     &   -0.00000015451139637086D0, -0.00000203340624738438D0,  &
     &    0.00001947204525295057D0,  0.00002854147231653228D0,  &
     &   -0.00101565063152200272D0,  0.00271187003520095655D0,  &
     &    0.02328095035422810727D0, -0.16725021123116877197D0,  &
     &    0.32490054966649436974D0 / 
      DATA (B(I), I = 39, 51) /                                &
     &    0.00000000002319363370D0, -0.00000000006303206648D0, & 
     &   -0.00000000264888267434D0,  0.00000002050708040581D0, & 
     &    0.00000011371857327578D0, -0.00000211211337219663D0, & 
     &    0.00000368797328322935D0,  0.00009823686253424796D0, & 
     &   -0.00065860243990455368D0, -0.00075285814895230877D0, & 
     &    0.02585434424202960464D0, -0.11637092784486193258D0, & 
     &    0.18267336775296612024D0 / 
      DATA (B(I), I = 52, 64) /                                &
     &   -0.00000000000367789363D0,  0.00000000020876046746D0, & 
     &   -0.00000000193319027226D0, -0.00000000435953392472D0, & 
     &    0.00000018006992266137D0, -0.00000078441223763969D0, & 
     &   -0.00000675407647949153D0,  0.00008428418334440096D0, & 
     &   -0.00017604388937031815D0, -0.00239729611435071610D0, & 
     &    0.02064129023876022970D0, -0.06905562880005864105D0, & 
     &    0.09084526782065478489D0 / 
!     **************************************************************************
      Y=ERF(X)    !FORTRAN 2008 INTRINSIC
      RETURN
!
      W = ABS(X)
      IF (W .LT. 2.2D0) THEN
          T = W * W
          K = INT(T)
          T = T - REAL(K,KIND=8)
          K = K * 13
          Y = (((((((((((( A(K    )*T + A(K+ 1)) * T + &
     &        A(K+ 2))*T + A(K+ 3))*T + A(K+ 4)) * T + &
     &        A(K+ 5))*T + A(K+ 6))*T + A(K+ 7)) * T + &
     &        A(K+ 8))*T + A(K+ 9))*T + A(K+10)) * T + &
     &        A(K+11))*T + A(K+12))*W
      ELSE IF (W .LT. 6.9D0) THEN
          K = INT(W)
          T = W - REAL(K,KIND=8)
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T + &
     &        B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T + &
     &        B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T + &
     &        B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T + &
     &        B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      ELSE
          Y = 1
      END IF
      IF(X.LT.0) Y=-Y
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$ERFCR8(X,Y)
!     **************************************************************************
!     **  COPYRIGHT(C) 1996 TAKUYA OOURA                                      **
!     **  (EMAIL: OOURA@MMM.T.U-TOKYO.AC.JP).                                 **
!     **  YOU MAY USE, COPY, MODIFY THIS CODE FOR ANY PURPOSE AND             **
!     **  WITHOUT FEE. YOU MAY DISTRIBUTE THIS ORIGINAL PACKAGE.              **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: Y
      REAL(8)     ,PARAMETER :: PA  = 3.97886080735226000D+00
      REAL(8)     ,PARAMETER :: P0  = 2.75374741597376782D-01
      REAL(8)     ,PARAMETER :: P1  = 4.90165080585318424D-01 
      REAL(8)     ,PARAMETER :: P2  = 7.74368199119538609D-01 
      REAL(8)     ,PARAMETER :: P3  = 1.07925515155856677D+00 
      REAL(8)     ,PARAMETER :: P4  = 1.31314653831023098D+00 
      REAL(8)     ,PARAMETER :: P5  = 1.37040217682338167D+00 
      REAL(8)     ,PARAMETER :: P6  = 1.18902982909273333D+00 
      REAL(8)     ,PARAMETER :: P7  = 8.05276408752910567D-01 
      REAL(8)     ,PARAMETER :: P8  = 3.57524274449531043D-01 
      REAL(8)     ,PARAMETER :: P9  = 1.66207924969367356D-02 
      REAL(8)     ,PARAMETER :: P10 =-1.19463959964325415D-01 
      REAL(8)     ,PARAMETER :: P11 =-8.38864557023001992D-02
      REAL(8)     ,PARAMETER :: P12 = 2.49367200053503304D-03 
      REAL(8)     ,PARAMETER :: P13 = 3.90976845588484035D-02 
      REAL(8)     ,PARAMETER :: P14 = 1.61315329733252248D-02 
      REAL(8)     ,PARAMETER :: P15 =-1.33823644533460069D-02 
      REAL(8)     ,PARAMETER :: P16 =-1.27223813782122755D-02 
      REAL(8)     ,PARAMETER :: P17 = 3.83335126264887303D-03 
      REAL(8)     ,PARAMETER :: P18 = 7.73672528313526668D-03 
      REAL(8)     ,PARAMETER :: P19 =-8.70779635317295828D-04 
      REAL(8)     ,PARAMETER :: P20 =-3.96385097360513500D-03 
      REAL(8)     ,PARAMETER :: P21 = 1.19314022838340944D-04 
      REAL(8)     ,PARAMETER :: P22 = 1.27109764952614092D-03
      REAL(8)                :: T,U
!     **************************************************************************
      Y=ERFC(X)    !FORTRAN 2008 INTRINSIC
      RETURN
!
      T = PA / (PA + ABS(X))
      U = T - 0.5D0
      Y = (((((((((P22 * U + P21) * U + P20) * U + &
     &    P19) * U + P18) * U + P17) * U + P16) * U + &
     &    P15) * U + P14) * U + P13) * U + P12 
      Y = ((((((((((((Y * U + P11) * U + P10) * U + &
     &    P9) * U + P8) * U + P7) * U + P6) * U + P5) * U + &
     &    P4) * U + P3) * U + P2) * U + P1) * U + P0) * T * &
     &    EXP(-X * X)
      IF (X .LT. 0) Y = 2 - Y
      RETURN
      END
!
!.......................................................................
MODULE RANDOM_MODULE
  INTEGER(8),SAVE        :: SEED=1
  INTEGER(8),PARAMETER   :: RANFAC1=22222221
  INTEGER(8),PARAMETER   :: RANFAC2=2**24
! CHOICE EXPLICIT CPPVARIABLE_XLF RNG
! INTEGER(8),SAVE        :: SEED=1_8
! INTEGER(8),PARAMETER   :: RANFAC1=44485709377909_8
! INTEGER(8),PARAMETER   :: RANFAC2=2_8**48
END MODULE RANDOM_MODULE
!
!     ..................................................................
      SUBROUTINE LIB$RANDOM(X)
!     ****************************************************************** 
!     **  RETURNS A RANDOM NUMBER                                     **
!     ******************************************************************
      USE RANDOM_MODULE
      IMPLICIT NONE
      REAL(8)   ,INTENT(OUT) :: X
!     ******************************************************************
      SEED=MODULO(RANFAC1*SEED,RANFAC2)
      X=REAL(SEED,KIND=8)/REAL(RANFAC2,KIND=8)
!     == CHOICE EXPLICIT CPPVARIABLE_XLF RNG
!     CALL RANDOM_NUMBER(Y)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE LIB$RANDOMSEED
!     ****************************************************************** 
!     **  CHOOSE A SEED  FOR THE NADOM NUMBER GENERATOR               **
!     ******************************************************************
      USE RANDOM_MODULE
      IMPLICIT NONE
!     ******************************************************************
      SEED=1
!     == CHOICE EXPLICIT CPPVARIABLE_XLF RNG
!     SEED=1_8
!     == CHOICE ORIGINAL CPPVARIABLE_XLF RNG
!     CALL RANDOM_SEED(GENERATOR=2)
!     CALL RANDOM_SEED
      RETURN
      END
!
!***********************************************************************
!***********************************************************************
!****                                                               ****
!****  TEST ROUTINES                                                ****
!****                                                               ****
!***********************************************************************
!***********************************************************************
!
!     ................................................................
      SUBROUTINE LIB$TEST()
      CALL LIB_TEST_INVERTR8()
      CALL LIB_TEST_DIAGR8()
      CALL LIB_TEST_DIAGC8()
      CALL LIB_TEST_GENERALEIGENVALUER8()
!     == COMPLEX GENERAL EIGENVALUE PROBLEM NOT IMPLEMENTED FOR ESSL
      CALL LIB_TEST_GENERALEIGENVALUEC8()
      CALL LIB_TEST_MATRIXSOLVER8()
      CALL LIB_TEST_MATRIXSOLVEC8()
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_INVERTR8()
      INTEGER(4),PARAMETER :: N=5
      REAL(8)              :: A(N,N)
      REAL(8)              :: AINV(N,N)
      REAL(8)              :: RES1(N,N)
      REAL(8)              :: RES2(N,N)
      INTEGER(4)           :: I
      LOGICAL   ,PARAMETER :: TPR=.FALSE.
!     ****************************************************************
      WRITE(*,FMT='("TEST LIB$INVERTR8")')
!
!     == MAKE INPUT DATA =============================================
      CALL RANDOM_NUMBER(A)
!
!     == WRITE INPUT DATA ============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("A:",10F10.3)')A(I,:)
        ENDDO
      END IF
!
!     == SOLVE PROBLEM ===============================================
      CALL LIB$INVERTR8(N,A,AINV)
!
!     == WRITE OUT DATA ==============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("AINV:",10F10.3)')AINV(I,:)
        ENDDO
      END IF
!
!     == TEST RESULT  ================================================
      RES1=MATMUL(A,AINV)
      RES2=MATMUL(AINV,A)
      DO I=1,N
        RES1(I,I)=RES1(I,I)-1.D0
        RES2(I,I)=RES2(I,I)-1.D0
      ENDDO
!      PRINT*,'TEST ',RES1
      IF(MAXVAL(ABS(RES1)+ABS(RES2)).LT.1.D-10) THEN
        WRITE(*,FMT='(T5,"OK")')
      ELSE
        WRITE(*,FMT='(T5,"FAILED")')
      END IF        
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_DIAGR8()
      INTEGER(4),PARAMETER :: N=5
      REAL(8)              :: H(N,N)
      REAL(8)              :: E(N)
      REAL(8)              :: U(N,N)
      REAL(8)              :: EMAT(N,N)
      REAL(8)              :: RES(N,N)
      INTEGER(4)           :: I
!     ****************************************************************
      WRITE(*,FMT='("TEST LIB$DIAGR8")')
      CALL RANDOM_NUMBER(H)
      H=H+TRANSPOSE(H)
!!$      DO I=1,N
!!$        WRITE(*,FMT='("H:",10F10.3)')H(I,:)
!!$      ENDDO
      CALL LIB$DIAGR8(N,H,E,U)
      EMAT=0.D0
      DO I=1,N
        EMAT(I,I)=E(I)
      ENDDO
      RES=MATMUL(H,U)-MATMUL(U,EMAT)
!      PRINT*,'TEST ',RES1
      IF(MAXVAL(ABS(RES)).LT.1.D-10) THEN
        WRITE(*,FMT='(T5,"OK")')
      ELSE
        WRITE(*,FMT='(T5,"FAILED")')
      END IF        
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_DIAGC8()
      INTEGER(4),PARAMETER :: N=5
      COMPLEX(8)           :: H(N,N)
      REAL(8)              :: E(N)
      COMPLEX(8)           :: U(N,N)
      COMPLEX(8)           :: EMAT(N,N)
      COMPLEX(8)           :: RES(N,N)
      REAL(8)              :: RE,IM
      INTEGER(4)           :: I,J
!     ****************************************************************
      WRITE(*,FMT='("TEST LIB$DIAGC8")')
      DO I=1,N
        DO J=1,N
          CALL RANDOM_NUMBER(RE)
          CALL RANDOM_NUMBER(IM)
          H(I,J)=CMPLX(RE,IM,KIND=8)
        ENDDO
      ENDDO
      H=H+TRANSPOSE(CONJG(H))
!!$      DO I=1,N
!!$        WRITE(*,FMT='("H:",10F10.3)')H(I,:)
!!$      ENDDO
      CALL LIB$DIAGC8(N,H,E,U)
      EMAT=(0.D0,0.D0)
      DO I=1,N
        EMAT(I,I)=CMPLX(E(I),0.D0,KIND=8)
      ENDDO
      RES=MATMUL(H,U)-MATMUL(U,EMAT)
!      PRINT*,'TEST ',RES
      IF(MAXVAL(ABS(RES)).LT.1.D-7) THEN
        WRITE(*,FMT='(T5,"OK")')
      ELSE
        WRITE(*,FMT='(T5,"FAILED")')
      END IF        
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_MATRIXSOLVER8()
      INTEGER(4),PARAMETER :: N=5
      INTEGER(4),PARAMETER :: M=7
      INTEGER(4),PARAMETER :: NEQ=2
      REAL(8)              :: ASQ(N,N)
      REAL(8)              :: BSQ(N,NEQ)
      REAL(8)              :: XSQ(N,NEQ)
      REAL(8)              :: A(N,M)
      REAL(8)              :: B(N,NEQ)
      REAL(8)              :: X(M,NEQ)
      INTEGER(4)           :: I
      REAL(8)              :: DEV
      REAL(8)  ,PARAMETER  :: TOL=1.D-10
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
!     ****************************************************************
      WRITE(*,FMT='("LIB_TEST_MATRIXSOLVER8")')
!
!     ================================================================
!     ================================================================
!     == TEST WITH SQUARE MATRIX                                    ==
!     ================================================================
!     ================================================================
!
!     == MAKE INPUT DATA =============================================
      CALL RANDOM_NUMBER(ASQ)
      CALL RANDOM_NUMBER(BSQ)
!
!     == WRITE INPUT DATA ============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("A:",10F10.3)')ASQ(I,:)
        ENDDO
        DO I=1,N
          WRITE(*,FMT='("B:",10F10.3)')BSQ(I,:)
        ENDDO
      END IF
!
!     == SOLVE PROBLEM ===============================================
      CALL LIB$MATRIXSOLVER8(N,N,NEQ,ASQ,XSQ,BSQ)
!
!     == WRITE OUTPUT DATA ===========================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("X:",10F10.3)')XSQ(I,:)
        ENDDO
      END IF
!
!     == TEST RESULT =================================================
      DEV=MAXVAL(ABS(MATMUL(ASQ,XSQ)-BSQ))
      IF(DEV.LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",E12.5)')DEV
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",E12.5)')DEV
      END IF        
!
!     ================================================================
!     ================================================================
!     == TEST WITH GENERAL RECTANGULAR MATRIX                       ==
!     ================================================================
!     ================================================================
!
!     == MAKE INPUT DATA =============================================
      CALL RANDOM_NUMBER(A)
      CALL RANDOM_NUMBER(B)
!
!     == WRITE INPUT DATA ============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("A:",10F10.3)')ASQ(I,:)
        ENDDO
        DO I=1,N
          WRITE(*,FMT='("B:",10F10.3)')BSQ(I,:)
        ENDDO
      END IF
!
!     == SOLVE PROBLEM ===============================================
      CALL LIB$MATRIXSOLVER8(N,M,NEQ,A,X,B)
!
!     == WRITE OUTPUT DATA ===========================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("X:",10F10.3)')XSQ(I,:)
        ENDDO
      END IF
!
!     == TEST RESULT =================================================
      DEV=MAXVAL(ABS(MATMUL(A,X)-B))
      IF(DEV.LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",E12.5)')DEV
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",E12.5)')DEV
      END IF        
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_MATRIXSOLVEC8()
      INTEGER(4),PARAMETER :: N=5
      INTEGER(4),PARAMETER :: M=7
      INTEGER(4),PARAMETER :: NEQ=2
      COMPLEX(8)           :: ASQ(N,N)
      COMPLEX(8)           :: XSQ(N,NEQ)
      COMPLEX(8)           :: A(N,M)
      COMPLEX(8)           :: B(N,NEQ)
      COMPLEX(8)           :: X(M,NEQ)
      REAL(8)              :: RASQ(N,N)
      REAL(8)              :: RB(N,NEQ)
      REAL(8)              :: RA(N,M)
      INTEGER(4)           :: I
      REAL(8)              :: DEV
      REAL(8)  ,PARAMETER  :: TOL=1.D-10
      COMPLEX(8),PARAMETER :: CI=(0.D0,1.D0)
      LOGICAL(4),PARAMETER :: TPR=.FALSE.
!     ****************************************************************
      WRITE(*,FMT='("LIB_TEST_MATRIXSOLVEC8")')
!
!     ================================================================
!     ================================================================
!     == TEST WITH SQUARE MATRIX                                    ==
!     ================================================================
!     ================================================================
!
!     == MAKE INPUT DATA =============================================
      CALL RANDOM_NUMBER(RASQ)
      ASQ=RASQ
      CALL RANDOM_NUMBER(RASQ)
      ASQ=ASQ+RASQ*CI
      CALL RANDOM_NUMBER(RB)
      B=RB
      CALL RANDOM_NUMBER(RB)
      B=B+RB*CI

!
!     == WRITE INPUT DATA ============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("A:",10F10.3)')ASQ(I,:)
        ENDDO
        DO I=1,N
          WRITE(*,FMT='("B:",10F10.3)')B(I,:)
        ENDDO
      END IF
!
!     == SOLVE PROBLEM ===============================================
      CALL LIB$MATRIXSOLVEC8(N,N,NEQ,ASQ,XSQ,B)
!
!     == WRITE OUTPUT DATA ===========================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("X:",10F10.3)')XSQ(I,:)
        ENDDO
      END IF
!
!     == TEST RESULT =================================================
      DEV=MAXVAL(ABS(MATMUL(ASQ,XSQ)-B))
      IF(DEV.LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",E12.5)')DEV
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",E12.5)')DEV
      END IF        
!
!     ================================================================
!     ================================================================
!     == TEST WITH GENERAL RECTANGULAR MATRIX                       ==
!     ================================================================
!     ================================================================
!
!     == MAKE INPUT DATA =============================================
      CALL RANDOM_NUMBER(RA)
      A=RA-0.5D0
      CALL RANDOM_NUMBER(RA)
      A=A+(RA-0.5D0)*CI
      CALL RANDOM_NUMBER(RB)
      B=(RB-0.5D0)
      CALL RANDOM_NUMBER(RB)
      B=B+(RB-0.5D0)*CI
!
!     == WRITE INPUT DATA ============================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("A:",10F10.3)')A(I,:)
        ENDDO
        DO I=1,N
          WRITE(*,FMT='("B:",10F10.3)')B(I,:)
        ENDDO
      END IF
!
!     == SOLVE PROBLEM ===============================================
      CALL LIB$MATRIXSOLVEC8(N,M,NEQ,A,X,B)
!
!     == WRITE OUTPUT DATA ===========================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("X:",10F10.3)')X(I,:)
        ENDDO
      END IF
!
!     == TEST RESULT =================================================
      DEV=MAXVAL(ABS(MATMUL(A,X)-B))
!      DEV=MAX(MAXVAL(ABS(REAL(MATMUL(A,X)-B))),MAXVAL(ABS(REAL(CI*(MATMUL(A,X)-B)))))
      IF(DEV.LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",E12.5)')DEV
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",E12.5)')DEV
      END IF  
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_GENERALEIGENVALUER8()
      INTEGER(4),PARAMETER :: N=5
      REAL(8)              :: H(N,N)
      REAL(8)              :: S(N,N)
      REAL(8)              :: E(N)
      REAL(8)              :: U(N,N)
      REAL(8)              :: EMAT(N,N)
      REAL(8)              :: DEV
      INTEGER(4)           :: I
      LOGICAL   ,PARAMETER :: TPR=.FALSE.
      REAL(8)   ,PARAMETER :: TOL=1.D-8
!     ****************************************************************
      WRITE(*,FMT='("LIB_TEST_GENERALEIGENVALUEC8")')
!
!     ================================================================
!     == SET UP INPUT DATA                                          ==
!     ================================================================
      CALL RANDOM_NUMBER(H)
      H=0.5D0*(H+TRANSPOSE(H))
      CALL RANDOM_NUMBER(S)
      S=MATMUL(S,TRANSPOSE(S))
      
!!$H=0.D0
!!$H(1,2)=1.D0
!!$H(2,1)=1.D0
!!$S=0.D0
!!$DO I=1,N
!!$  S(I,I)=1.D0
!!$ENDDO
!
!     ================================================================
!     == WRITE INPUT                                                ==
!     ================================================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("H:",10F10.3)')H(I,:)
        ENDDO
        WRITE(*,*)'----------------- '
        DO I=1,N
          WRITE(*,FMT='("S:",10F10.3)')S(I,:)
        ENDDO
      END IF
!
!     ================================================================
!     == SOLVE PROBLEM                                              ==
!     ================================================================
      CALL LIB$GENERALEIGENVALUER8(N,H,S,E,U)
!
!     ================================================================
!     == WRITE OUTPUT                                               ==
!     ================================================================
      IF(TPR) THEN
        WRITE(*,*)
        WRITE(*,FMT='("E:",10F10.3)')E
        WRITE(*,*)
        DO I=1,N
          WRITE(*,FMT='("U:",10F10.3)')U(I,:)
        ENDDO
      END IF
!
!     ================================================================
!     == TEST RESULT                                                ==
!     ================================================================
      EMAT=0.D0
      DO I=1,N
        EMAT(I,I)=E(I)
      ENDDO
      DEV=MAXVAL(ABS(MATMUL(H,U)-MATMUL(S,MATMUL(U,EMAT))))
      IF(DEV.LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",E12.3)')DEV
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",E12.3)')DEV
      END IF        
      RETURN
      END
!
!     ...............................................................
      SUBROUTINE LIB_TEST_GENERALEIGENVALUEC8()
      INTEGER(4),PARAMETER :: N=5
      COMPLEX(8)           :: H(N,N)
      COMPLEX(8)           :: S(N,N)
      REAL(8)              :: E(N)
      COMPLEX(8)           :: U(N,N)
      COMPLEX(8)           :: EMAT(N,N)
      REAL(8)              :: RANMAT(N,N)
      REAL(8)              :: DEV1,DEV2
      INTEGER(4)           :: I
      LOGICAL   ,PARAMETER :: TPR=.FALSE.
      REAL(8)   ,PARAMETER :: TOL=1.D-6
      COMPLEX(8),PARAMETER :: CI=(0.D0,1.D0)
!     ****************************************************************
      WRITE(*,FMT='("LIB_TEST_GENERALEIGENVALUER8")')
!
!     ================================================================
!     == SET UP INPUT DATA                                          ==
!     ================================================================
      CALL RANDOM_NUMBER(RANMAT)
      H=RANMAT
      CALL RANDOM_NUMBER(RANMAT)
      H=H+RANMAT*CI
      H=0.5D0*(H+TRANSPOSE(CONJG(H)))
      CALL RANDOM_NUMBER(RANMAT)
      S=RANMAT
      CALL RANDOM_NUMBER(RANMAT)
      S=S+RANMAT*CI
      S=MATMUL(S,TRANSPOSE(CONJG(S)))
!
!     ================================================================
!     == WRITE INPUT                                                ==
!     ================================================================
      IF(TPR) THEN
        DO I=1,N
          WRITE(*,FMT='("H:",10("(",F10.3,",",F10.3,")"))')H(I,:)
        ENDDO
        WRITE(*,*)
        DO I=1,N
          WRITE(*,FMT='("S:",10("(",F10.3,",",F10.3,")"))')S(I,:)
        ENDDO
        WRITE(*,*)
      END IF
!
!     ================================================================
!     == SOLVE PROBLEM                                              ==
!     ================================================================
      CALL LIB$GENERALEIGENVALUEC8(N,H,S,E,U)
!
!     ================================================================
!     == WRITE OUTPUT                                               ==
!     ================================================================
      IF(TPR) THEN
        WRITE(*,FMT='("E:",10F10.3)')E
        WRITE(*,*)
        DO I=1,N
          WRITE(*,FMT='("U:",10("(",F10.3,",",F10.3,")"))')U(I,:)
        ENDDO
      END IF
!
!     ================================================================
!     == TEST RESULT                                                ==
!     ================================================================
!     == TEST EIGENVALUE PROBLEM
      EMAT=0.D0
      DO I=1,N
        EMAT(I,I)=CMPLX(E(I),0.D0,KIND=8)
      ENDDO
      DEV1=MAXVAL(ABS(MATMUL(H,U)-MATMUL(S,MATMUL(U,EMAT))))
!     == TEST ORTHONORMALITY
      EMAT=0.D0
      DO I=1,N
        EMAT(I,I)=CMPLX(1.D0,0.D0,KIND=8)
      ENDDO
      DEV2=MAXVAL(ABS(MATMUL(TRANSPOSE(CONJG(U)),MATMUL(S,U))-EMAT))
!     ==
      IF(MAX(DEV1,DEV2).LT.TOL) THEN
        WRITE(*,FMT='(T5,"OK: DEV=",2E10.3)')DEV1,DEV2
      ELSE
        WRITE(*,FMT='(T5,"FAILED: DEV=",2E10.3)')DEV1,DEV2
      END IF        
      RETURN
      END
! 
!*******************************************************************************
!*******************************************************************************
!****                                                                      *****
!****              INTERFACES FOR SPECIAL FUNCTIONS                        *****
!****                                                                      *****
!****     INTERFACE FOR BASIC SPECIAL FUNCTIONS FROM SLATEC                *****
!****     SEE SLATEC.F FOR FURTHER DETAILS                                 *****
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DBESJ(L,X,Y)
!     **************************************************************************
!     **  BESSEL FUNCTION OF FIRST KIND                                       **
!     **  SEE SLATEC.F FOR FURTHER DETAILS                                    **
!     **                                                                      **
!     **  CAUTION: THE FORTRAN INTRINSICS ALLOW ONLY INTEGER ORDER L,         **
!     **           WHILE WE MOSTLY NEED HALF-INTEGER ORDER                    **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)           :: L   ! ORDER 
      REAL(8),INTENT(IN)           :: X   ! ARGUMENT
      REAL(8),INTENT(OUT)          :: Y   ! BESSEL FUNCTIOHN OF FIRST KIND
      INTEGER(4)                   :: NZ  ! NZ.NEQ.0 IN OUTPUT: UNDERFLOW
!     **************************************************************************
!     __fortran 2008 function bessel_jn
!     __bessel function of the first kind of order n
!     __this routine is an elemental function and can be called for arrays
!     __there is also a version bessel_jn(n1,n2,x) which allows a range of 
!     __orders to be specified
!     y=bessel_jn(l,x) !fortran 2008 intrinsic function

      CALL DBESJ(X,L,1,Y,NZ)
      IF(NZ.NE.0)THEN
        Y=0.D0
      END IF
      RETURN
      END SUBROUTINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DBESY(L,X,Y)
!     **************************************************************************
!     **  BESSEL FUNCTION OF SECOND KIND                                      **
!     **  SEE SLATEC.F FOR FURTHER DETAILS                                    **
!     **                                                                      **
!     **  CAUTION: THE FORTRAN INTRINSICS ALLOW ONLY INTEGER ORDER L,         **
!     **           WHILE WE MOSTLY NEED HALF-INTEGER ORDER                    **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)           :: L
      REAL(8),INTENT(IN)           :: X
      REAL(8),INTENT(OUT)          :: Y
!     **************************************************************************
      CALL DBESY(X,L,1,Y)

!     y=bessel_yn(l,x)  !fortran 2008 intrinsic function
      RETURN
      END SUBROUTINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DBESI(L,X,Y)
!     **************************************************************************
!     **  MODIFIED BESSEL FUNCTION OF FIRST KIND                              **
!     **  SEE SLATEC.F FOR FURTHER DETAILS                                    **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)           :: L
      REAL(8),INTENT(IN)           :: X
      REAL(8),INTENT(OUT)          :: Y
      INTEGER(4)                   :: NZ
!     ================================================================
      CALL DBESI(X,L,1,1,Y,NZ)
      IF(NZ.NE.0)THEN
        Y=0.0D0
!        CALL ERROR$MSG('UNDERFLOW OR OVERFLOW')
!        CALL ERROR$STOP('LIB$DBESI')
      END IF
      RETURN
      END SUBROUTINE
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LIB$DBESK(L,X,Y)
!     **************************************************************************
!     **  MODIFIED BESSEL FUNCTION OF THIRD KIND                              **
!     **  SEE SLATEC.F FOR FURTHER DETAILS                                    **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN)           :: L
      REAL(8),INTENT(IN)           :: X
      REAL(8),INTENT(OUT)          :: Y
      INTEGER(4)                   :: NZ
!     ================================================================
      CALL DBESK(X,L,1,1,Y,NZ)
      IF(NZ.NE.0)THEN
        Y=0.0D0
!        CALL ERROR$MSG('UNDERFLOW OR OVERFLOW')
!        CALL ERROR$STOP('LIB$DBESK')
      END IF
      RETURN
      END SUBROUTINE
