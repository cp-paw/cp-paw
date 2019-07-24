!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SPECIALFUNCTION$TEST()
!      *************************************************************************
!      **                                                                     **
!      **                                                                     **
!      **                                                                     **
!      *************************************************************************
       IMPLICIT NONE
       INTEGER(4),PARAMETER :: LX=3
       INTEGER(4),PARAMETER :: NX=200
       INTEGER(4),PARAMETER :: NFIL=101
       REAL(8)   ,PARAMETER :: XMAX=10.D0
       REAL(8)   ,PARAMETER :: NXREAL=NX
       REAL(8)   ,PARAMETER :: DX=XMAX/(NXREAL-1.D0) 
       REAL(8)              :: X
       REAL(8)              :: Y(LX+1)
       REAL(8)              :: DYDX(LX+1)
       INTEGER(4)           :: IX,L
!      *************************************************************************
!
!      =========================================================================
!      ==  BESSEL FUNCTIONS                                                   ==
!      =========================================================================
       OPEN(NFIL,FILE='TEST_BESSEL.DAT')
       DO IX=1,NX
         X=DX*REAL(IX-1,KIND=8)
         DO L=0,LX
           CALL SPFUNCTION$BESSEL(L,X,Y(L+1),DYDX(L+1))       
         ENDDO
         WRITE(NFIL,*)X,Y,DYDX
       ENDDO
       CLOSE(NFIL)
!
!      =========================================================================
!      ==  MODIFIED BESSEL FUNCTIONS                                          ==
!      =========================================================================
       OPEN(NFIL,FILE='TEST_MODBESSEL.DAT')
       DO IX=1,NX
         X=DX*REAL(IX-1,KIND=8)
         DO L=0,LX
           CALL SPFUNCTION$MODBESSEL(L,X,Y(L+1),DYDX(L+1))       
         ENDDO
         WRITE(NFIL,*)X,Y,DYDX
       ENDDO
       CLOSE(NFIL)
!
!      =========================================================================
!      ==  BESSEL FUNCTIONS FOR KAPPA=0                                       ==
!      =========================================================================
       OPEN(NFIL,FILE='TEST_BESSEL0.DAT')
       DO IX=1,NX
         X=DX*REAL(IX-1,KIND=8)
         DO L=0,LX
           CALL SPFUNCTION$BESSEL0(L,X,Y(L+1),DYDX(L+1))       
         ENDDO
         WRITE(NFIL,*)X,Y,DYDX
       ENDDO
       CLOSE(NFIL)
!
!      =========================================================================
!      ==  MODIFIED HANKEL FUNCTIONS                                          ==
!      =========================================================================
       OPEN(NFIL,FILE='TEST_MODHANKEL.DAT')
       DO IX=2,NX
         X=DX*REAL(IX-1,KIND=8)
         DO L=0,LX
           CALL SPFUNCTION$MODHANKEL(L,X,Y(L+1),DYDX(L+1))       
         ENDDO
         WRITE(NFIL,*)X,Y,DYDX
       ENDDO
       CLOSE(NFIL)
!
!      =========================================================================
!      ==  NEUMANN FUNCTIONS                                                  ==
!      =========================================================================
       OPEN(NFIL,FILE='TEST_NEUMANN.DAT')
       DO IX=2,NX
         X=DX*REAL(IX-1,KIND=8)
         DO L=0,LX
           CALL SPFUNCTION$NEUMANN(L,X,Y(L+1),DYDX(L+1))       
         ENDDO
         WRITE(NFIL,*)X,Y,DYDX
       ENDDO
       CLOSE(NFIL)
!
!      =========================================================================
!      ==  MODIFIED NEUMANN FUNCTIONS                                         ==
!      =========================================================================
       OPEN(NFIL,FILE='TEST_MODNEUMANN.DAT')
       DO IX=2,NX
         X=DX*REAL(IX-1,KIND=8)
         DO L=0,LX
           CALL SPFUNCTION$MODNEUMANN(L,X,Y(L+1),DYDX(L+1))       
         ENDDO
         WRITE(NFIL,*)X,Y,DYDX
       ENDDO
       CLOSE(NFIL)
!
!      =========================================================================
!      ==  NEUMANN FUNCTIONS FOR KAPPA=0                                      ==
!      =========================================================================
       OPEN(NFIL,FILE='TEST_NEUMANN0.DAT')
       DO IX=2,NX
         X=DX*REAL(IX-1,KIND=8)
         DO L=0,LX
           CALL SPFUNCTION$NEUMANN0(L,X,Y(L+1),DYDX(L+1))
         ENDDO
         WRITE(NFIL,*)X,Y,DYDX
       ENDDO
       CLOSE(NFIL)
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPECIALFUNCTION$ERF(X,VAL)
!     **************************************************************************
!     **  RETURNS THE ERROR FUNCTION ERF(X)                                   **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: X
      REAL(8),INTENT(OUT):: VAL
      REAL(8),EXTERNAL ::DERF
!     **************************************************************************
      CALL LIB$ERFR8(X,VAL)
      RETURN
      END
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SPECIALFUNCTION$GAMMP(X,A,VAL)
!      *************************************************************************
!      **  RETURNS THE INCOMPLETE GAMMA FUNCTION P(A,X)
!      **  SEE NUMERICAL RECIPES
!      *************************************************************************
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: X
       REAL(8),INTENT(IN) :: A
       REAL(8),INTENT(OUT):: VAL
!      *************************************************************************
       CALL ERROR$MSG('FUNCTION NOT USED. IT IS MARKED FOR DELETION')
       CALL ERROR$STOP('SPECIALFUNCTION$GAMMP')
!
       IF(X.LT.0.D0.OR.A.LE.0.D0) THEN
         STOP 'ERROR STOP IN SPECIALFUNCTION$GAMMP'
       END IF
       IF(X.LT.A+1.D0) THEN
!        == USE SERIES EXPANSION ===================================
         CALL SPECIALFUNCTION_GSER(X,A,VAL)
       ELSE
!        == USE CONTINUED FRACTION REPRESENTATION ==================
         CALL SPECIALFUNCTION_GCF(X,A,VAL)
         VAL=1.D0-VAL
       END IF
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SPECIALFUNCTION_GAMMLN(X,VAL)
!      ** LN(GAMMA(XX))
!      ** SEE NUMERICAL RECIPES
       IMPLICIT NONE
       REAL(8),INTENT(IN) :: X
       REAL(8),INTENT(OUT):: VAL
       REAL(8),PARAMETER  :: COF(6)=(/76.18009173D0,-86.50532033D0 &
      &                             ,24.01409822D0,-1.231739516D0 &
      &                             ,0.120858003D-2,-0.536382D-5/)
       REAL(8),PARAMETER  :: STP=2.50662827465D0
       REAL(8)            :: X1,TMP,SER
       INTEGER(4)         :: J
!      *************************************************************
       CALL ERROR$MSG('FUNCTION NOT USED. IT IS MARKED FOR DELETION')
       CALL ERROR$STOP('SPECIALFUNCTION$GAMMP')
!
       X1=X-1.D0
       TMP=X1+5.5D0
       TMP=(X1+0.5D0)*LOG(TMP)-TMP
       SER=1.D0
       DO J=1,6
         X1=X1+1.D0
         SER=SER+COF(J)/X1
       ENDDO
       VAL=TMP+LOG(STP*SER)
       RETURN
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SPECIALFUNCTION_GSER(X,A,VAL)
!      ** 
!      ** INCOMPLETE GAMMA FUNCTION P(A,X) EVALUATED BY ITS SERIES EXPANSION
!      ** 
!      ** SEE NUMERICAL RECIPES
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: X
       REAL(8)   ,INTENT(IN) :: A
       REAL(8)   ,INTENT(OUT):: VAL
       INTEGER(4),PARAMETER  :: ITMAX=100
       REAL(8)   ,PARAMETER  :: EPS=3.D-7
       REAL(8)               :: AP,SUM,DEL,GLN
       INTEGER(4)            :: N
!      *************************************************************
       CALL ERROR$MSG('FUNCTION NOT USED. IT IS MARKED FOR DELETION')
       CALL ERROR$STOP('SPECIALFUNCTION$GAMMP')
!
       IF(X.LT.0.D0) THEN
         STOP 'INVALID ARGUMENT FOR GSER'
       END IF
       IF(X.EQ.0.D0) THEN
         VAL=0.D0
         RETURN
       END IF
       AP=A
       SUM=1.D0/A
       DEL=SUM
       DO N=1,ITMAX
         AP=AP+1.D0
         DEL=DEL*X/AP
         SUM=SUM+DEL
         IF(ABS(DEL).LT.ABS(SUM)*EPS) THEN
           CALL SPECIALFUNCTION_GAMMLN(A,GLN)
           VAL=SUM*EXP(-X+A*LOG(X)-GLN)
           RETURN
         END IF
       ENDDO
       STOP 'A TOO LARGE, ITMAX TOO SMALL; STOP IN GSER'
       END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE SPECIALFUNCTION_GCF(X,A,VAL)
!      ** 
!      ** INCOMPLETE GAMMA FUNCTION Q(A,X) EVALUATED BY ITS CONTINUED 
!      **  FRACTION REPRESENTATION
!      ** 
!      ** SEE NUMERICAL RECIPES
       IMPLICIT NONE
       REAL(8)   ,INTENT(IN) :: X
       REAL(8)   ,INTENT(IN) :: A
       REAL(8)   ,INTENT(OUT):: VAL
       INTEGER(4),PARAMETER  :: ITMAX=100
       REAL(8)   ,PARAMETER  :: EPS=3.D-7
       REAL(8)               :: GOLD,GLN
       REAL(8)               :: A0,A1,B0,B1
       REAL(8)               :: FAC,AN,ANA,ANF,G
       INTEGER(4)            :: N
!      *************************************************************
       CALL ERROR$MSG('FUNCTION NOT USED. IT IS MARKED FOR DELETION')
       CALL ERROR$STOP('SPECIALFUNCTION$GAMMP')
!
       CALL SPECIALFUNCTION_GAMMLN(A,GLN)
       GOLD=0.D0
       A0=1.D0
       A1=X
       B0=0.D0
       B1=1.D0
       FAC=1.D0
       DO N=1,ITMAX
         AN=REAL(N,KIND=8)
         ANA=AN-A
         A0=(A1+A0*ANA)*FAC
         B0=(B1+B0*ANA)*FAC
         ANF=AN*FAC
         A1=X*A0+ANF*A1
         B1=X*B0+ANF*B1
         IF(A1.NE.0.D0) THEN
           FAC=1.D0/A1
           G=B1*FAC
           IF(ABS((G-GOLD)/G).LT.EPS) THEN
             VAL=EXP(-X+A*LOG(X)-GLN)*G
             RETURN
            END IF
            GOLD=G
          ENDIF
        ENDDO
        STOP 'A TOO LARGE, ITMAX SOO SMALL, STOP IN GCF'
        END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPECIALFUNCTION$BESSELOLD(L,X,Y)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION                            **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **    FORMULA 10.1.2              FOR   X < L                           **
!     **    FORMULA 10.1.8 AND  10.1.9  FOR   X > L                           **
!     **                                                                      **
!     ** TODO: THIS ROUTINE IS DUPLICATED BY SPFUNCTION$BESSEL                **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     ****************************************** P.E. BLOECHL, 1991 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN AGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! BESSL FUNCTION AT X
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: TRIG(4)
      REAL(8)               :: FACUL(0:100)
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: ARG
      REAL(8)               :: XSQ
      REAL(8)               :: FAC
      INTEGER(4)            :: I,K,IL,II,ISVAR
!     ******************************************************************
      IF(X.GT.REAL(L,KIND=8)) THEN
        ARG=X-0.5D0*REAL(L,KIND=8)*PI
        TRIG(1)=SIN(ARG)/X
        TRIG(2)=COS(ARG)/X
        TRIG(3)=-TRIG(1)
        TRIG(4)=-TRIG(2)
        Y=TRIG(1)
        IF(L.EQ.0) RETURN
!       ==  DOUBLE FACULTY FACUL(L)=(2*L)!!
        FACUL(0)=1.D0
        DO I=1,2*L
          FACUL(I)=FACUL(I-1)*REAL(I,KIND=8)
        ENDDO
        XSQ=0.5D0/X
        FAC=1.D0
        DO K=1,L
          II=MOD(K,4)+1
          FAC=FACUL(K+L)/FACUL(K)/FACUL(L-K)*XSQ**K
!         FAC=FAC*XSQ*DBLE(L+K)/DBLE(K*(L-K))
          Y=Y+FAC*TRIG(II)
        ENDDO
!       II=MOD(L,4)+1
!       FAC=FAC*XSQ*DBLE(2*L)/DBLE(L)
!       Y=Y+FAC*TRIG(II)
        RETURN
      END IF
!     ==================================================================
!     ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS                        ==
!     ==================================================================
      ISVAR=1
      DO IL=1,L
        ISVAR=ISVAR*(2*IL+1)
      ENDDO
      IF(L.NE.0.D0) THEN
        FAC=X**L/DBLE(ISVAR)
      ELSE
        FAC=1.D0/DBLE(ISVAR)
      END IF
      Y=FAC
      XSQ=-0.5D0*X*X
      ISVAR=2*L+1
      DO I=1,1000
        ISVAR=ISVAR+2
        FAC=FAC*XSQ/DBLE(I*ISVAR)
        Y=Y+FAC
        IF(ABS(FAC).LT.TOL) GOTO 9999
      ENDDO
      CALL ERROR$MSG('Y NOT CONVERGED')
      CALL ERROR$I4VAL('L',L)
      CALL ERROR$R8VAL('X',X)
      CALL ERROR$STOP('SPECIALFUNCTION$BESSEL')
9999  CONTINUE
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPECIALFUNCTION$BESSEL(L,X,Y)
!     **************************************************************************
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION                            **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 1991 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN AGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! BESSL FUNCTION AT X
      REAL(8)               :: DYDX
!     ******************************************************************
      CALL SPFUNCTION$BESSEL(L,X,Y,DYDX)
      RETURN
      END
!
!===============================================================================
!===============================================================================
!===============================================================================
!=====     THE FOLLOWING ROUTINES NEED TO BE TESTED.                          ==
!===============================================================================
!===============================================================================
!===============================================================================
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$BESSELOLD(L,X,Y,DYDX)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION                            **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **    FORMULA 10.1.2              FOR   X < L                           **
!     **    FORMULA 10.1.8 AND  10.1.9  FOR   X > L                           **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     ****************************************** P.E. BLOECHL, 2009 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN AGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! BESSL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: TRIG(4)
      REAL(8)               :: DTRIG(4)
      REAL(8)               :: FACUL(0:2*L)
      REAL(8)               :: ARG
      REAL(8)               :: XSQ
      REAL(8)               :: FAC,DFAC
      INTEGER(4)            :: I,K,IL,II,ISVAR
      LOGICAL(4)            :: CONVG
!     **************************************************************************
      IF(X.LT.0.D0) THEN
        CALL ERROR$MSG('BESSEL FUNCTION FOR NEGATIVE ARG UNDEFINED')
        CALL ERROR$STOP('SPFUNCTION$BESSEL')
      ELSE IF(X.LE.REAL(L,KIND=8)) THEN
!       ========================================================================
!       ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS  X<L                         ==
!       ========================================================================
        ISVAR=1
        DO IL=1,L
          ISVAR=ISVAR*(2*IL+1)
        ENDDO
        IF(L.EQ.0) THEN
          FAC=1.D0/REAL(ISVAR,KIND=8)
          DFAC=0.D0
        ELSE IF(L.EQ.1) THEN
          FAC=X/REAL(ISVAR,KIND=8)
          DFAC=1.D0/REAL(ISVAR,KIND=8)
        ELSE 
          FAC=X**L/REAL(ISVAR,KIND=8)
          DFAC=REAL(L,KIND=8)*X**(L-1)/REAL(ISVAR,KIND=8)
        END IF
        Y=FAC
        DYDX=DFAC
        XSQ=-0.5D0*X*X
        ISVAR=2*L+1
        DO I=1,1000
          ISVAR=ISVAR+2
!         = DO DERIVATIVE BEFORE VALUE BECAUSE VALUE IS MODIFIED ===============
          DFAC=(DFAC*XSQ-FAC*X)/REAL(I*ISVAR,KIND=8)
          FAC=FAC*XSQ/REAL(I*ISVAR,KIND=8)
          Y=Y+FAC
          DYDX=DYDX+DFAC
          CONVG=ABS(FAC).LT.TOL
          IF(CONVG) EXIT
        ENDDO
        IF(.NOT.CONVG) THEN 
          CALL ERROR$MSG('Y NOT CONVERGED')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$R8VAL('X',X)
          CALL ERROR$STOP('SPFUNCTION$BESSEL')
        END IF
      ELSE 
!       ========================================================================
!       ==  EXPANSION FOR LARGE ARGUMENTS  X>L                                ==
!       ========================================================================
        ARG=X-0.5D0*REAL(L,KIND=8)*PI
        TRIG(1)=SIN(ARG)/X
        TRIG(2)=COS(ARG)/X
        DTRIG(1)=COS(ARG)/X-SIN(ARG)/X**2
        DTRIG(2)=-SIN(ARG)/X-COS(ARG)/X**2
        TRIG(3)=-TRIG(1)
        TRIG(4)=-TRIG(2)
        DTRIG(3)=-DTRIG(1)
        DTRIG(4)=-DTRIG(2)
        Y=TRIG(1)
        DYDX=DTRIG(1)
        IF(L.EQ.0) RETURN
!       ==  DOUBLE FACULTY FACUL(L)=(2*L)!!
        FACUL(0)=1.D0
        DO I=1,2*L
          FACUL(I)=FACUL(I-1)*REAL(I,KIND=8)
        ENDDO
        XSQ=0.5D0/X
        FAC=1.D0
        DO K=1,L
          II=MOD(K,4)+1
          FAC=FACUL(K+L)/FACUL(K)/FACUL(L-K)*XSQ**K
          DFAC=-REAL(K,KIND=8)*FAC/X
!         FAC=FAC*XSQ*DBLE(L+K)/DBLE(K*(L-K))
          Y=Y+FAC*TRIG(II)
          DYDX=DYDX+DFAC*TRIG(II)+FAC*DTRIG(II)
        ENDDO
!       II=MOD(L,4)+1
!       FAC=FAC*XSQ*DBLE(2*L)/DBLE(L)
!       Y=Y+FAC*TRIG(II)
        RETURN
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$BESSELOLD
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$BESSEL(L,X,Y,DYDX)
!     **************************************************************************
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION                            **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     ******************************************R. SCHADE, 2013*****************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN AGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! BESSL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE
      REAL(8)               :: JLPHALB,JLMHALB,JLPDREIHALB
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: XMIN=1.D-30
!     **************************************************************************
      IF(X.LT.0.D0) THEN
        CALL ERROR$MSG('BESSEL FUNCTION FOR NEGATIVE ARG UNDEFINED')
        CALL ERROR$STOP('SPFUNCTION$BESSEL')
      ELSE IF(X.LT.XMIN)THEN
        IF(L.EQ.0)THEN
          Y=1.0D0
        ELSE
          Y=0.0D0
        ENDIF
        IF(L.EQ.1)THEN
          DYDX=1.0D0/3.0D0
        ELSE
          DYDX=0.0D0
        ENDIF
      ELSE
        CALL LIB$DBESJ(REAL(L,KIND=8)+0.5D0,X,JLPHALB)
        Y=0.5D0*SQRT(2.D0*PI/X)*JLPHALB
        IF(L.EQ.0)THEN
          CALL LIB$DBESJ(REAL(L,KIND=8)+1.5D0,X,JLPDREIHALB)
          DYDX=-SQRT(0.5D0*PI/X)*JLPDREIHALB
          !DYDX=(COS(X)-SIN(X)/X)/X
        ELSE
          CALL LIB$DBESJ(REAL(L,KIND=8)-0.5D0,X,JLMHALB)
          CALL LIB$DBESJ(REAL(L,KIND=8)+1.5D0,X,JLPDREIHALB)
          DYDX=-0.25D0*SQRT(2.0D0*PI/X)*(JLPHALB/X-JLMHALB+JLPDREIHALB) 
        ENDIF
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$BESSEL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$NEUMANNOLD(L,X,Y,DYDX)
!     **************************************************************************
!     **                                                                      **
!     **  SPHERICAL BESSEL FUNCTION OF THE SECOND KIND Y_L(X)                 **
!     **    Y = -X^L*(-1/X D/DX)^L [COS(X)/X]  FORMULA 10.1.26                **
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL NEUMANN FUNCTION                           **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **    FORMULA 10.1.3   FOR   X < L                                      **
!     **    FORMULA 10.1.9   FOR   X > L                                      **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! NEUMANN FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE OF NEUMANN FUNCTION AT X
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: TRIG(4)
      REAL(8)               :: DTRIG(4)
      REAL(8)               :: FACUL(0:2*L)
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: ARG
      REAL(8)               :: XSQ,DXSQ
      REAL(8)               :: FAC,DFAC
      INTEGER(4)            :: I,K,IL,II,ISVAR
      REAL(8)               :: M1POWERL
!     **************************************************************************
      IF(X.LT.0.D0) THEN
        CALL ERROR$MSG('NEUMANN FUNCTION NOT DEEFINED FOR NEGATIVE ARG')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$NEUMANN')
      ELSE IF(X.EQ.0.D0) THEN
        CALL ERROR$MSG('NEUMANN FUNCTION WITH ZERO ARGUMENT DIVERGES')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$NEUMANN')
      ELSE IF(X.LE.REAL(L,KIND=8)) THEN
!       ========================================================================
!       ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS   X<L                        ==
!       ========================================================================
        ISVAR=-1
        DO IL=2,L
          ISVAR=ISVAR*(2*IL-1)
        ENDDO
        FAC=REAL(ISVAR,KIND=8)/X**(L+1)
        DFAC=-REAL(L+1,KIND=8)*FAC/X
        Y=FAC
        DYDX=DFAC
        XSQ=-0.5D0*X*X
        DXSQ=-X
        ISVAR=-(2*L+1)
        DO I=1,1000
          ISVAR=ISVAR+2
          DFAC=(DFAC*XSQ-X*FAC)/REAL(I*ISVAR,KIND=8)
          FAC=FAC*XSQ/REAL(I*ISVAR,KIND=8)
          Y=Y+FAC
          DYDX=DYDX+DFAC
          IF(ABS(FAC).LT.TOL) RETURN
        ENDDO
        CALL ERROR$MSG('Y NOT CONVERGED')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$NEUMANN')
      ELSE 
!       ========================================================================
!       ==  EXPANSION FOR LARGE ARGUMENTS   X>L                               ==
!       ========================================================================
        ARG=X+0.5D0*REAL(L,KIND=8)*PI
        M1POWERL=(-1)**L
        TRIG(1)=-M1POWERL*COS(ARG)/X
        DTRIG(1)=M1POWERL*(SIN(ARG)+COS(ARG)/X)/X
        TRIG(2)=M1POWERL*SIN(ARG)/X
        DTRIG(2)=M1POWERL*(COS(ARG)-SIN(ARG)/X)/X
        TRIG(3)=-TRIG(1)
        DTRIG(3)=-DTRIG(1)
        TRIG(4)=-TRIG(2)
        DTRIG(4)=-DTRIG(2)
        Y=TRIG(1)
        DYDX=DTRIG(1)
        IF(L.EQ.0) RETURN
!       ==  DOUBLE FACULTY FACUL(L)=(2*L)!!
        FACUL(0)=1.D0
        DO I=1,2*L
          FACUL(I)=FACUL(I-1)*REAL(I,KIND=8)
        ENDDO
        XSQ=0.5D0/X
        FAC=1.D0
        DO K=1,L
          II=MOD(K,4)+1
          FAC=FACUL(K+L)/FACUL(K)/FACUL(L-K)*XSQ**K
          DFAC=-FAC*REAL(K,KIND=8)/X
          Y=Y+FAC*TRIG(II)
          DYDX=DYDX+FAC*DTRIG(II)+DFAC*TRIG(II)
        ENDDO
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$NEUMANNOLD
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$NEUMANN(L,X,Y,DYDX)
!     **************************************************************************
!     **  SPHERICAL BESSEL FUNCTION OF THE SECOND KIND Y_L(X)                 **
!     **    Y = -X^L*(-1/X D/DX)^L [COS(X)/X]  FORMULA 10.1.26                **
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL NEUMANN FUNCTION                           **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     **  CAUTION: THE NEUMANN FUNCTION IS DEFINED WITH AN OPPOSITE SIGN!!!   **
!     **    THIS IS THE NEUMANN FUNCTION MULTIPLIED WITH -1!!                 **
!     **    NEUMANN FUNCTION AND BESSEL FUNCTION OF THE SECOND KIND ARE       **
!     **    DEFINED WITH OPPOSITE SIGN.  (REMARK PBLOECHL FEB. 15,2014)       **
!     ******************************************R. SCHADE, 2013*****************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! NEUMANN FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE OF NEUMANN FUNCTION AT X
      REAL(8)               :: YLPHALB,YLMHALB,YLPDREIHALB
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8),PARAMETER     :: XMIN=1.D-30
!     **************************************************************************
      IF(X.LT.0.D0) THEN
        CALL ERROR$MSG('NEUMANN FUNCTION NOT DEEFINED FOR NEGATIVE ARG')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$NEUMANN')
      ELSE IF(X.EQ.0.D0) THEN
        CALL ERROR$MSG('NEUMANN FUNCTION WITH ZERO ARGUMENT DIVERGES')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$NEUMANN')
      ELSE
        CALL LIB$DBESY(REAL(L,KIND=8)+0.5D0,X,YLPHALB)
        Y=0.5D0*SQRT(2.D0*PI/X)*YLPHALB
        IF(L.EQ.0)THEN
          DYDX=(SIN(X)+COS(X)/X)/X
        ELSE
          CALL LIB$DBESY(REAL(L,KIND=8)-0.5D0,X,YLMHALB)
          CALL LIB$DBESY(REAL(L,KIND=8)+1.5D0,X,YLPDREIHALB)
          DYDX=-0.25D0*SQRT(2.0D0*PI/X)*(YLPHALB/X-YLMHALB+YLPDREIHALB)  
        ENDIF
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$NEUMANN
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$MODBESSELOLD(L,X,Y,DYDX)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE MODIFIED SPHERICAL BESSEL FUNCTION                   **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **    FORMULA 10.2.5   FOR   X < L+1                                    **
!     **    FORMULA 10.2.9   FOR   X > L+1                                    **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! MODIFIED BESSEL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE
      REAL(8)               :: TRIG(2),DTRIG(2)
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: XSQ  !X-SQUARE
      REAL(8)               :: FAC,DFAC
      INTEGER(4)            :: I,K,IL,ISVAR
!     **************************************************************************
      IF(X.GT.REAL(L+1,KIND=8)) THEN
        TRIG(:)=1.D0
        DTRIG(:)=0.D0
        XSQ=0.5D0/X
        Y=XSQ*(EXP(X)-(-1.D0)**L*EXP(-X))
        DYDX=-Y/X+XSQ*(EXP(X)+(-1.D0)**L*EXP(-X))
        DO K=1,L
           FAC=REAL((L+K)*(L-K+1),KIND=8)/REAL(K,KIND=8)*XSQ
           DFAC=-FAC/X
           DTRIG(1)=-DTRIG(1)*FAC-TRIG(1)*DFAC
           TRIG(1)=-TRIG(1)*FAC
           DTRIG(2)=DTRIG(2)*FAC+TRIG(2)*DFAC
           TRIG(2)=+TRIG(2)*FAC
           DYDX=DYDX-XSQ/X*(TRIG(1)*EXP(X)-(-1)**L*TRIG(2)*EXP(-X)) &
    &                +XSQ*((DTRIG(1)+TRIG(1))*EXP(X) &
               -(-1.D0)**L*(DTRIG(2)-TRIG(2))*EXP(-X))
           Y=Y+XSQ*(TRIG(1)*EXP(X)-(-1.D0)**L*TRIG(2)*EXP(-X))
        ENDDO
!
!     ==========================================================================
!     ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS                                ==
!     ==========================================================================
      ELSE
        ISVAR=1
        DO IL=1,L
          ISVAR=ISVAR*(2*IL+1)
        ENDDO
        IF(L.EQ.0) THEN
          FAC=1.D0/REAL(ISVAR,KIND=8)
          DFAC=0.D0
        ELSE IF(L.EQ.1) THEN
          FAC=X/REAL(ISVAR,KIND=8)
          DFAC=1.D0/REAL(ISVAR,KIND=8)
        ELSE
          FAC=X**L/REAL(ISVAR,KIND=8)
          DFAC=REAL(L,KIND=8)*X**(L-1)/REAL(ISVAR,KIND=8)
        END IF
        Y=FAC
        DYDX=DFAC
        XSQ=0.5D0*X*X
        ISVAR=2*L+1
        DO I=1,1000
          ISVAR=ISVAR+2
          DFAC=(DFAC*XSQ+FAC*X)/REAL(I*ISVAR,KIND=8)
          FAC=FAC*XSQ/REAL(I*ISVAR,KIND=8)
          Y=Y+FAC
          DYDX=DYDX+DFAC
          IF(DABS(FAC).LT.TOL) RETURN
        ENDDO
        CALL ERROR$MSG('Y NOT CONVERGED')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$MODBESSEL')
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$MODBESSELOLD
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$MODBESSEL(L,X,Y,DYDX)
!     **************************************************************************
!     **  CALCULATES THE MODIFIED SPHERICAL BESSEL FUNCTION                   **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     ******************************************R. SCHADE, 2013*****************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! MODIFIED BESSEL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE
      REAL(8)               :: ILPHALB,ILMHALB,ILPDREIHALB
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8),PARAMETER     :: XMIN=1.D-30
!     **************************************************************************
      IF(X.LT.XMIN)THEN
        IF(L.EQ.0)THEN
          Y=1.0D0
          DYDX=0.0D0
        ELSE IF (L.EQ.1) THEN
          Y=0.0D0
          DYDX=1.0D0/3.0D0
        ELSE
          Y=0.0D0
          DYDX=0.0D0
        ENDIF
      ELSE
        CALL LIB$DBESI(REAL(L,KIND=8)+0.5D0,X,ILPHALB)
        Y=0.5D0*SQRT(2.D0*PI/X)*ILPHALB
        IF(L.EQ.0)THEN
          CALL LIB$DBESI(1.5D0,X,ILMHALB)
          DYDX=0.5D0*SQRT(2.0D0*PI/X)*ILMHALB
        ELSE
          CALL LIB$DBESI(REAL(L,KIND=8)-0.5D0,X,ILMHALB)
          CALL LIB$DBESI(REAL(L,KIND=8)+1.5D0,X,ILPDREIHALB)
          DYDX=-0.25D0*SQRT(2.0D0*PI/X)*(ILPHALB/X-ILMHALB-ILPDREIHALB)  
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SPFUNCTION$MODBESSEL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$MODNEUMANNOLD(L,X,Y,DYDX)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE MODIFIED SPHERICAL NEUMANN FUNCTION                  **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **    FORMULA 10.2.6   FOR   X < L                                      **
!     **    FORMULA 10.2.10  FOR   X > L                                      **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! MODIFIED NEUMANN FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE
      REAL(8)               :: TRIG(2),DTRIG(2)
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: XSQ
      REAL(8)               :: FAC,DFAC
      INTEGER(4)            :: I,K,IL,ISVAR
      LOGICAL(4)            :: CONVG
      REAL(8)               :: SVAR,DSVAR
!     **************************************************************************
      IF(X.LE.0.D0) THEN
        CALL ERROR$MSG('NOT DEFINED FOR ZERO OR NEGATIVE ARGUMENTS')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$MODNEUMANN')
      ELSE IF(X.LE.REAL(L,KIND=8)) THEN
!       ========================================================================
!       ==  TAYLOR EXPANSION FOR SMALL ARGUMENTS                              ==
!       ========================================================================
        ISVAR= 1
        DO IL=2,L
          ISVAR=ISVAR*(2*IL-1)
        ENDDO
        FAC=REAL(ISVAR,KIND=8)/X**(L+1)/(-1.D0)**L
        DFAC=-REAL(L+1,KIND=8)*FAC/X
        Y=FAC
        DYDX=DFAC
        XSQ=0.5D0*X*X
        ISVAR=-(2*L+1)
        CONVG=.FALSE.
        DO I=1,1000
          ISVAR=ISVAR+2
          DFAC=(DFAC*XSQ+FAC*X)/REAL(I*ISVAR,KIND=8)
          FAC=FAC*XSQ/REAL(I*ISVAR,KIND=8)
          Y=Y+FAC
          DYDX=DYDX+DFAC
          CONVG=ABS(FAC).LT.TOL
          IF(CONVG) EXIT
        ENDDO
        IF(.NOT.CONVG) THEN
          CALL ERROR$MSG('Y NOT CONVERGED')
          CALL ERROR$I4VAL('L',L)
          CALL ERROR$R8VAL('X',X)
          CALL ERROR$STOP('SPFUNCTION$MODNEUMANN')
        END IF
      ELSE 
!       ========================================================================
!       ==  LARGE ARGUMENTS                                                   ==
!       ========================================================================
        TRIG(:)=1.D0
        DTRIG(:)=0.D0
        XSQ=0.5D0/X
        Y=XSQ*(EXP(X)+(-1.D0)**L*EXP(-X))
        DYDX=-Y/X+XSQ*(EXP(X)-(-1.D0)**L*EXP(-X))
        DO K=1,L
           FAC=REAL((L+K)*(L-K+1),KIND=8)/REAL(K, KIND=8)*XSQ
           DFAC=-FAC/X
           DTRIG(1)=-DTRIG(1)*FAC-TRIG(1)*DFAC
           TRIG(1)=-TRIG(1)*FAC
           DTRIG(2)= DTRIG(2)*FAC+TRIG(2)*DFAC
           TRIG(2)= TRIG(2)*FAC
           SVAR=XSQ*(TRIG(1)*DEXP(X)+(-1.D0)**L*TRIG(2)*DEXP(-X))
           DSVAR=-SVAR/X+XSQ*(           (DTRIG(1)+TRIG(1))*EXP(X) &
       &                     +(-1.D0)**L*(DTRIG(2)-TRIG(2))*EXP(-X))
           Y=Y+SVAR
           DYDX=DYDX+DSVAR
        ENDDO
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$MODNEUMANNOLD
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$MODNEUMANN(L,X,Y,DYDX)
!     **************************************************************************
!     **  CALCULATES THE MODIFIED SPHERICAL NEUMANN FUNCTION                  **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                                  **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     ******************************************R. SCHADE, 2013*****************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! MODIFIED NEUMANN FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8),PARAMETER     :: XMIN=1.D-30
      REAL(8)               :: MODHANKEL_Y,MODHANKEL_DYDX
      REAL(8)               :: MODBESSEL_Y,MODBESSEL_DYDX
!     *****************************************R.SCHADE, 2013*******************
      IF(X.LT.0.D0) THEN
        CALL ERROR$MSG('NEUMANN FUNCTION NOT DEEFINED FOR NEGATIVE ARG')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$NEUMANN')
      ELSE IF(X.EQ.0.D0) THEN
        CALL ERROR$MSG('NEUMANN FUNCTION WITH ZERO ARGUMENT DIVERGES')
        CALL ERROR$I4VAL('L',L)
        CALL ERROR$R8VAL('X',X)
        CALL ERROR$STOP('SPFUNCTION$NEUMANN')
      ELSE
        !STEGUN 10.2.4
        CALL SPFUNCTION$MODHANKEL(L,X,MODHANKEL_Y,MODHANKEL_DYDX)
        CALL SPFUNCTION$MODBESSEL(L,X,MODBESSEL_Y,MODBESSEL_DYDX)
        Y   =-2.0D0/PI*(-1)**(L+1)*MODHANKEL_Y   +MODBESSEL_Y
        DYDX=-2.0D0/PI*(-1)**(L+1)*MODHANKEL_DYDX+MODBESSEL_DYDX
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$MODNEUMANN
!!$!
!!$!     .......................................................................
!!$      SUBROUTINE SPFUNCTION$HANKEL(L,X,Y,DYDX)
!!$!     ***********************************************************************
!!$!     **                                                                   **
!!$!     **  CALCULATES THE SPHERICAL HANKEL FUNCTION.                        **
!!$!     **  CONVENTION: ABRAMOWITZ, EQUATION 10.1.1                          **
!!$!     **                                                                   **
!!$!     ***********************************************************************
!!$      IMPLICIT NONE
!!$      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
!!$      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
!!$      COMPLEX(8),INTENT(OUT):: Y ! HANKEL FUNCTION AT X
!!$      COMPLEX(8),INTENT(OUT):: DYDX ! HANKEL FUNCTION AT X
!!$      REAL(8)               :: YREAL, YIMAG,DYREAL, DYIMAG
!!$!     ***********************************************************************
!!$      CALL SPFUNCTION$BESSEL (L,X,YREAL,DYREAL)
!!$      CALL SPFUNCTION$NEUMANN(L,X,YIMAG,DYIMAG)
!!$      Y    = CMPLX(YREAL,YIMAG)
!!$      DYDX = CMPLX(DYREAL,DYIMAG)
!!$      END SUBROUTINE SPFUNCTION$HANKEL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$MODHANKELOLD(L,X,Y,DYDX)
!     **************************************************************************
!     **  MODIFIED SPHERICAL BESSEL FUNCTION OF THE THIRD KIND K_L(X)         **
!     **  AS DEFINED IN ABRAMOWITZ/STEGUN EQ. 10.2.4                          **
!     **                                                                      **
!     **    Y = 0.5*PI*X^L*(-1/X D/DX)^L [EXP(-X)/X]                          **
!     **                                                                      **
!     **  CALCULATES THE MODIFIED SPHERICAL HANKEL FUNCTION                   **
!     **  AS DEFINED BY ABRAMOWITZ AND STEGUN (EQUATION 10.2.4)               **
!     **  USING EQUATION 10.2.15 FOR X>L                                      **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! MODIFIED HANKEL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! MODIFIED HANKEL FUNCTION AT X
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: TRIG,DTRIG
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: XSQ
      REAL(8)               :: FAC,DFAC
      REAL(8)               :: JVAL,JDER
      REAL(8)               :: NVAL,NDER
      REAL(8)               :: SVAR,DSVAR
      INTEGER(4)            :: K
!     **************************************************************************
      IF(X.LE.0.D0) THEN
        CALL ERROR$MSG('UNDEFINED FOR NEGATIVE OR ZERO ARGUMENTS')
        CALL ERROR$STOP('SPFUNCTION$MODHANKEL')
      ELSE IF(X.LE.REAL(L,KIND=8)) THEN
         CALL SPFUNCTION$MODBESSEL (L,X,JVAL,JDER)
         CALL SPFUNCTION$MODNEUMANN(L,X,NVAL,NDER)
         Y = 0.5D0*PI * (-1.D0)**(L+1) * (JVAL - NVAL)
         DYDX = 0.5D0*PI * (-1.D0)**(L+1) * (JDER - NDER)
      ELSE 
        TRIG=1.D0
        DTRIG=0.D0
        XSQ=0.5D0/X
        Y=XSQ*PI*EXP(-X)
        DYDX=-Y/X-Y
        DO K=1,L
           FAC=REAL((L+K)*(L-K+1),KIND=8)/REAL(K,KIND=8)*XSQ
           DFAC=-FAC/X
           DTRIG=DTRIG*FAC+TRIG*DFAC
           TRIG=TRIG*FAC
           SVAR=XSQ*PI*TRIG*EXP(-X)
           DSVAR=-SVAR/X-SVAR+SVAR*DTRIG/TRIG
           Y=Y+SVAR
           DYDX=DYDX+DSVAR
        ENDDO
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$MODHANKELOLD
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$MODHANKEL(L,X,Y,DYDX)
!     **************************************************************************
!     **  MODIFIED SPHERICAL BESSEL FUNCTION OF THE THIRD KIND K_L(X)         **
!     **  AS DEFINED IN ABRAMOWITZ/STEGUN EQ. 10.2.4                          **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     ******************************************R. SCHADE, 2013*****************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! MODIFIED HANKEL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! MODIFIED HANKEL FUNCTION AT X
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)   ,PARAMETER  :: XMIN=1.D-30
      REAL(8)               :: KLPHALB,KLMHALB,KLPDREIHALB
!     **************************************************************************
      IF(X.LE.0.D0) THEN
        CALL ERROR$MSG('UNDEFINED FOR NEGATIVE OR ZERO ARGUMENTS')
        CALL ERROR$STOP('SPFUNCTION$MODHANKEL')
      ELSE 
        CALL LIB$DBESK(REAL(L,KIND=8)+0.5D0,X,KLPHALB)
        Y=0.5D0*SQRT(2.D0*PI/X)*KLPHALB
        IF(L.EQ.0)THEN
          DYDX=-0.5D0*PI*EXP(-X)*(1.0D0/X+1.0D0/X**2)
        ELSE
          CALL LIB$DBESK(REAL(L,KIND=8)-0.5D0,X,KLMHALB)
          CALL LIB$DBESK(REAL(L,KIND=8)+1.5D0,X,KLPDREIHALB)
          DYDX=-0.25D0*SQRT(2.0D0*PI/X)*(KLPHALB/X+KLMHALB+KLPDREIHALB)  
        ENDIF
      END IF
      RETURN
      END SUBROUTINE SPFUNCTION$MODHANKEL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$BESSEL0OLD(L,X,Y,DYDX)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION FOR K=0.                   **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! BESSEL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE
      INTEGER(4)            :: FAC, I
!     **************************************************************************
      IF(L.EQ.0) THEN
        Y = 1.D0
        DYDX=1.D0
      ELSE
        FAC=1
        DO I=1, 2*L+1, 2
          FAC = FAC*I
        END DO
        Y = X**L / REAL(FAC,KIND=8)
        DYDX=REAL(L,KIND=8)*X**(L-1)/REAL(FAC,KIND=8)
      END IF
      END SUBROUTINE SPFUNCTION$BESSEL0OLD
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$BESSEL0(L,X,Y,DYDX)
!     **************************************************************************
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION FOR K=0.                   **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     ******************************************R. SCHADE, 2013*****************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! BESSEL FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE
      REAL(8)               :: FAC
      INTEGER(4)            :: I
      REAL(8)   ,PARAMETER  :: XMIN=1.D-30
!     **************************************************************************
      IF(L.EQ.0) THEN
        Y = 1.D0
        DYDX=0.D0
      ELSE
        IF(X.LT.XMIN)THEN
          IF(L.EQ.0)THEN
            Y=1.0D0
            DYDX=0.0D0
          ELSE IF(L.EQ.1) THEN
            Y=0.0D0
            DYDX=1.0D0/3.0D0
          ELSE
            Y=0.0D0
            DYDX=0.0D0
          ENDIF
        ELSE
          FAC=1.D0
          DO I=1,2*L+1,2
            FAC = FAC*REAL(I,KIND=8)
          END DO
          Y = X**L/FAC
          DYDX=REAL(L,KIND=8)*X**(L-1)/FAC
        ENDIF
      END IF
      END SUBROUTINE SPFUNCTION$BESSEL0
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$NEUMANN0OLD(L,X,Y,DYDX)
!     **************************************************************************
!     **                                                                      **
!     **  CALCULATES THE SPHERICAL NEUMANN FUNCTION FOR K=0.                  **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! NEUMANN FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE
      INTEGER(4)            :: FAC, I
!     **************************************************************************
      IF(X.LE.0.D0) THEN
        CALL ERROR$MSG('NOT DEFINED FOR ZERO OR NEGATIVE ARGUMENTS')
        CALL ERROR$STOP('SPFUNCTION$NEUMANN0')
      END IF
      FAC=1
      DO I=1,2*L-1,2
         FAC=FAC*I
      END DO
      Y = -REAL(FAC,KIND=8) / X**(L+1)
      DYDX = -REAL(L+1,KIND=8)*Y/X
      RETURN
      END SUBROUTINE SPFUNCTION$NEUMANN0OLD
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE SPFUNCTION$NEUMANN0(L,X,Y,DYDX)
!     **************************************************************************
!     **  CALCULATES THE SPHERICAL NEUMANN FUNCTION FOR KAPPA=0.              **
!     **  ALSO SEE NOTES IN 13_04_29_TEST_SPECIAL_FUNCTIONS.PDF               **
!     **                                                                      **
!     ******************************************R. SCHADE, 2013*****************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)   ,INTENT(OUT):: Y ! NEUMANN FUNCTION AT X
      REAL(8)   ,INTENT(OUT):: DYDX ! DERIVATIVE
      REAL(8)               :: FAC
      INTEGER(4)            :: I
!     **************************************************************************
      IF(X.LE.0.D0) THEN
        CALL ERROR$MSG('NOT DEFINED FOR ZERO OR NEGATIVE ARGUMENTS')
        CALL ERROR$STOP('SPFUNCTION$NEUMANN0')
      END IF
      FAC=1.D0
      DO I=1,2*L-1,2
        FAC=FAC*REAL(I,KIND=8)
      END DO
      Y = -FAC/ X**(L+1)
      DYDX = -REAL(L+1,KIND=8)*Y/X
      RETURN
      END SUBROUTINE SPFUNCTION$NEUMANN0
