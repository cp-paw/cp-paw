!
!*************************************************************************
!*************************************************************************
!**                                                                     **
!**  OBJECT: POLYNOM                                                    **
!**                                                                     **
!**  A POLYNOM IS DEFINED BY THE ORDER NI-1 AN ORIGIN X0 AND            **
!**  AN ARRAY OF N-COEFFICIENTS CI(I)                                   **
!**     P(X)=SUM_I=1^N: CI(I)*X**(I-1)                                  **
!**                                                                     **
!**  THE OBJECT POLYNOM ALLOWS TO OBTAIN AN INTERPOLATING POLYNOM FROM  **
!**  A SET OF POINTS (X,Y) THROUGH WHICH THE POLYNOM PASSES.            **
!**  THEN THE ORIGIN CAN BE CHANGED, AND VALUES, DERIVATIVES AND        **
!**  INTEGRALS CAN BE EVALUATED AT ARBITRARY POINTS.                    **
!**                                                                     **
!*************************************************************************
!*************************************************************************
!
!     ....................................................................
      SUBROUTINE POLYNOM$COEFF(NX,X0,COEFF,X,Y)
!     **                                                                **
!     **  OBTAINES THE COEFFICIENTS OF A POLYNOMIAL                     **
!     **  FROM DISCRETE POINTS                                          **
!     **     Y(R)=SUM_I=1^NX: COEFF(I)*(X-X0)**(I-1)                    **
!     **  SO THAT                                                       **
!     **     Y(X(I))=Y(I)                                               **
!     **  X0 IS CHOSEN TO BE THE FIRST GRID POINT                       **
!     **                                                                **
!     ********************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      REAL(8)   ,INTENT(OUT):: X0
      REAL(8)   ,INTENT(OUT):: COEFF(NX)
      REAL(8)   ,INTENT(IN) :: X(NX)
      REAL(8)   ,INTENT(IN) :: Y(NX)
      REAL(8)               :: YREL(NX)  ! Y-Y(1)
      REAL(8)               :: XREL(NX)  ! X-X(1)
      REAL(8)               :: DC(NX),A,B
      REAL(8)               :: SVAR,SVAR1,SVAR2
      INTEGER(4)            :: I,J,K
      LOGICAL(4)            :: TTEST=.FALSE.
!     ********************************************************************
      X0=X(1)
      COEFF(:)=0.D0
      SVAR=Y(1)
      COEFF(1)=SVAR
      DO I=1,NX
        YREL(I)=Y(I)-SVAR
        XREL(I)=X(I)-X0
      ENDDO
!
      DO I=2,NX
!       ** F(I)*PROD^(I-1):[A(J)+B(J)*(R(I)-R(1)]
        DC(1)=YREL(I)
        DO J=1,I-1
          B=1.D0/(XREL(I)-XREL(J))
          A=-B*XREL(J)
          SVAR1=0.D0
          DO K=1,J
            SVAR2=DC(K)
            DC(K)=SVAR2*A+SVAR1
            SVAR1=SVAR2*B
          ENDDO
          DC(J+1)=SVAR1
        ENDDO
!       == UPDATE COEFFICIENTS =========================================
        DO J=1,I
          COEFF(J)=COEFF(J)+DC(J)
        ENDDO
!       == UPDATE FUNCTION VALUES ======================================
        DO J=I+1,NX
          SVAR=1.D0
          DO K=1,I
            YREL(J)=YREL(J)-DC(K)*SVAR
            SVAR=SVAR*XREL(J)
          ENDDO
        ENDDO
      ENDDO
!     ==================================================================      
!     ==  TEST                                                        ==      
!     ==================================================================      
      IF(TTEST) THEN
        DO I=1,NX
          SVAR=0.D0
          DO J=1,NX
            SVAR=SVAR+COEFF(J)*(X(I)-X(1))**(J-1)
          ENDDO
          WRITE(*,FMT='(I5,3E15.5)')I,SVAR,Y(I),Y(I)-SVAR
        ENDDO
        WRITE(*,FMT='("COEFF",10E15.5)')COEFF
      END IF
      RETURN
      END SUBROUTINE POLYNOM$COEFF
!
!     ....................................................................
      SUBROUTINE POLYNOM$SHIFTORIGIN(NX,X0,COEFF,X0NEU)
!     **                                                                **
!     **  SHIFTS THE EXPANSION POINT X0 OF A POLYNOMIAL TO X0NEU        **
!     **     Y(R)=SUM_I=1^NX: COEFF(I)*(X-X0)**(I-1)                    **
!     **                                                                **
!     ********************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NX
      REAL(8)   ,INTENT(IN)    :: X0NEU
      REAL(8)   ,INTENT(INOUT) :: X0
      REAL(8)   ,INTENT(INOUT) :: COEFF(NX)
      REAL(8)                  :: C0(NX)
      REAL(8)                  :: DX
      INTEGER(4)               :: I,J
!     ********************************************************************
      DX=X0NEU-X0
      DO I=NX,2,-1
        C0=COEFF
        DO J=I,NX
          COEFF(J-1)=COEFF(J-1)+DX*C0(J)
        ENDDO
      ENDDO
      X0=X0NEU
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE POLYNOM$PRODUCT(N1,X1,C1,N2,X2,C2,N3,X3,C3)
!     **                                                                **
!     **  FORMS THE PRODUCT OF TWO POLYNOMS                             **
!     **                                                                **
!     **  THE RECOMMENDED VALUE FOR N3 IS N3=N1+N2-1. FOR SMALLER       **
!     **  VALUES HIGHER ORDER WILL BE TRUNCATED AND THE RESULT WILL     **
!     **  DEPEND ON THE VALUE OF X3.                                    **
!     **                                                                **
!     ********************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: N1
      REAL(8)   ,INTENT(IN)    :: X1
      REAL(8)   ,INTENT(IN)    :: C1(N1)
      INTEGER(4),INTENT(IN)    :: N2
      REAL(8)   ,INTENT(IN)    :: X2
      REAL(8)   ,INTENT(IN)    :: C2(N2)
      INTEGER(4),INTENT(IN)    :: N3
      REAL(8)   ,INTENT(IN)    :: X3
      REAL(8)   ,INTENT(OUT)   :: C3(N3)
      REAL(8)                  :: X1S,C1S(N1)
      REAL(8)                  :: X2S,C2S(N2)
      REAL(8)                  :: SUM
      INTEGER(4)               :: I,J,J1,J2
!     ********************************************************************
      X1S=X1
      C1S(:)=C1(:)
      CALL POLYNOM$SHIFTORIGIN(N1,X1S,C1S,X3)
      X2S=X2
      C2S(:)=C2(:)
      CALL POLYNOM$SHIFTORIGIN(N2,X2S,C2S,X3)
      DO I=1,N3
        J1=MAX(1,I+1-N2)
        J2=MIN(N1,I)
        SUM=0.D0
        DO J=J1,J2
          SUM=SUM+C1S(J)*C2S(I-J+1)
        ENDDO
        C3(I)=SUM
      ENDDO
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE POLYNOM$VALUE(NX,X0,COEFF,X,Y)
!     **                                                                **
!     **  EVALUATES THE VALUE Y OF THE POLYNOMIAL AT POINT X            **
!     **                                                                **
!     ********************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NX
      REAL(8)   ,INTENT(IN)    :: X0
      REAL(8)   ,INTENT(IN)    :: COEFF(NX)
      REAL(8)   ,INTENT(IN)    :: X
      REAL(8)   ,INTENT(OUT)   :: Y
      REAL(8)                  :: FAC
      REAL(8)                  :: DX
      INTEGER(4)               :: I
!     ********************************************************************
      DX=X-X0
      Y=0.D0
      FAC=1.D0
      DO I=1,NX
        Y=Y+COEFF(I)*FAC
        FAC=FAC*DX
      ENDDO
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE POLYNOM$DERIVATIVE(NX,X0,COEFF,X,Y)
!     **                                                                **
!     **  EVALUATES THE DERIVATIVE OF THE POLYNOMIAL AT POINT X         **
!     **                                                                **
!     ********************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NX
      REAL(8)   ,INTENT(IN)    :: X0
      REAL(8)   ,INTENT(IN)    :: COEFF(NX)
      REAL(8)   ,INTENT(IN)    :: X
      REAL(8)   ,INTENT(OUT)   :: Y
      REAL(8)                  :: FAC
      REAL(8)                  :: DX
      INTEGER(4)               :: I
!     ********************************************************************
      DX=X-X0
      Y=0.D0
      FAC=1.D0
      DO I=1,NX-1
        Y=Y+COEFF(I+1)*FAC*REAL(I)
        FAC=FAC*DX
      ENDDO
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE POLYNOM$INTEGRAL(NX,X0,COEFF,X,Y)
!     **                                                                **
!     **  EVALUATES THE INTEGRAL OF THE POLYNOMIAL AT POINT X           **
!     **                                                                **
!     ********************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)    :: NX
      REAL(8)   ,INTENT(IN)    :: X0
      REAL(8)   ,INTENT(IN)    :: COEFF(NX)
      REAL(8)   ,INTENT(IN)    :: X
      REAL(8)   ,INTENT(OUT)   :: Y
      REAL(8)                  :: FAC
      REAL(8)                  :: DX
      INTEGER(4)               :: I
!     ********************************************************************
      DX=X-X0
      Y=0.D0
      FAC=DX
      DO I=1,NX
        Y=Y+COEFF(I)*FAC/REAL(I)
        FAC=FAC*DX
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POLYNOM$ZEROS(NX,COEFF,Z)
!     **************************************************************************
!     ** DETERMINE THE ZEROS OF A POLYNOMIAL OF ORDER N WITH N=1,2,3          **
!     ** WITH REAL COEFFICIENTS                                               **
!     **                                                                      **
!     ** SEE ROUTINE POLYNOM$ZEROSC8 FOR POLYNOMIALS OF HIGHER ORDER WITH     **
!     ** COMPLEX COEFFICIENTS                                                 **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      REAL(8)   ,INTENT(IN) :: COEFF(NX)
      COMPLEX(8),INTENT(OUT):: Z(NX-1)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)   ,PARAMETER  :: BY3=1.D0/3.D0
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)            :: N ! ORDER OF THE POLYNOM
      REAL(8)               :: A,B,C,D,DD,P,Q,SVAR,U,V
      COMPLEX(8)            :: CSVAR
!     **************************************************************************
      N=NX-1
      Z(:)=(0.D0,0.D0)
!
!     ==========================================================================
!     == LINEAR EQUATION                                                      ==
!     ==========================================================================
      IF(N.EQ.1) THEN
!       == Y(X)=A*X+B ==========================================================
        A=COEFF(2)
        B=COEFF(1)
        Z(1)=-B/A
!
!     ==========================================================================
!     == QUADRATIC EQUATION                                                   ==
!     =========================================================================
      ELSE IF(N.EQ.2) THEN
!       == AX^2+BX+C ===========================================================
        A=COEFF(3)
        B=COEFF(2)
        C=COEFF(1)
        DD=B**2-4.D0*A*C
        Z(1)=SQRT(CMPLX(DD,KIND=8))
        Z(2)=-Z(1)
        Z(:)=(Z(:)-B)/(2.D0*A)
!
!     ==========================================================================
!     == CARDANOS EQUATION FOR CUBIC POLYNOMIAL
!     ==========================================================================
      ELSE IF(N.EQ.3) THEN
!       == AX^3+BX^2+CX+D ======================================================
        A=COEFF(4)
        B=COEFF(3)
        C=COEFF(2)
        D=COEFF(1)
        P=C/A-BY3*(B/A)**2
        Q=2.D0*(B*BY3/A)**3-BY3*B*C/A**2+D/A
        DD=(P*BY3)**3+(0.5D0*Q)**2        
        IF(DD.LE.0.D0) THEN
          SVAR=BY3 * ACOS(-0.5D0*Q*(-BY3*P)**(-1.5D0))
          Z(1)=COS(SVAR)
          Z(2)=COS(SVAR+BY3*2.D0*PI)
          Z(3)=COS(SVAR-BY3*2.D0*PI)
          Z=Z*SQRT(-BY3*4.D0*P)
        ELSE IF(DD.GT.0.D0) THEN
          U=( -0.5D0*Q+SQRT(DD) )**BY3
          V=-P*BY3/U
          CSVAR=EXP(2.D0*PI*CI*BY3)
         Z(1)=CMPLX(U+V,KIND=8)
         Z(2)=U*CSVAR+V/CSVAR
         Z(3)=U/CSVAR+V*CSVAR
       END IF
       Z(:)=Z(:)-B*BY3/A
!
!     ==========================================================================
!     == NO SOLUTION FOR POLYNOMIALS WITH ORDER HIGHER THAN CUBIC             ==
!     ==========================================================================
      ELSE
        CALL ERROR$MSG('ORDER OF THE POLYNOM OUT OF RANGE')
        CALL ERROR$MSG('ZEROS CAN ONLY BE CALCULATED UP TO CUBIC POLYNOMS')
        CALL ERROR$STOP('POLYNOM$ZEROS')
      END IF     
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE POLYNOM$ZEROSC8(NX,COEFF,Z)
!     **************************************************************************
!     ** FIND THE ROOTS OF A COMPLEX POLYNOMIAL OF ORDER NX-1                 **
!     ** THE POLYNOMIAL IS SUM_{J=1}^NX A_J * X**(J-1)                        **
!     ** USES METHOD OF VIETA (CODE BY ROBERT SCHADE, MODIFIED PETER BLOECHL) **
!     **                                                                      **
!     ** ROBERT SUGGESTS TO LOOK INTO THE METHOD OF JENKINS-TRAUB FOR A       **
!     ** NUMERICALLY MORE STABLE METHOD, SHOULD THIS BE AN ISSUE              **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)               :: NX        ! #(COEFFICIENTS)
      COMPLEX(8),INTENT(IN)               :: COEFF(NX) ! COEFFICIENT ARRAY
      COMPLEX(8),INTENT(OUT)              :: Z(NX-1)   ! ZEROS
      LOGICAL(4),PARAMETER                :: TTEST=.FALSE.
      INTEGER(4)                          :: I,J
      COMPLEX(8),ALLOCATABLE              :: M(:,:)
      COMPLEX(8)                         :: CSVAR  !KIND=16 BEFORE
!     **************************************************************************
!     ==========================================================================
!     ==  VIETA FORMULAS
!     ==========================================================================
      ALLOCATE(M(NX-1,NX-1))
      M(:,:)=(0.D0,0.D0)
      DO I=1,NX-1
        M(1,I)=-COEFF(I+1)/COEFF(1)
        IF(I+1.LE.NX-1)M(I+1,I)=(1.D0,0.D0)
      ENDDO
!
!     ==========================================================================
!     == DETERMINE ZEROS FROM EIGENVALUE PROBLEM                              ==
!     ==========================================================================
      CALL LIB__EIGVALNONHERMITEANC8(NX-1,M,Z)
      Z=(1.D0,0.D0)/Z
!
!     ==========================================================================
!     == TEST RESULT                                                          ==
!     ==========================================================================
      IF(TTEST) THEN
         WRITE(*,FMT='(" P(Z)=(",2F15.5,") * Z^",I2)')COEFF(1),0
        DO I=2,NX
          WRITE(*,FMT='("     +(",2F15.5,") * Z^",I2)')COEFF(I),I-1
        ENDDO
        DO I=1,NX-1  ! LOOP OVER ZEROS
          CSVAR=CMPLX(0.D0,0.D0,KIND=8)   ! KIND=16 BEFORE
          DO J=1,NX
            CSVAR=CSVAR+CMPLX(COEFF(J),KIND=8)*Z(I)**(J-1)  !KIND=16 BEFORE  
          ENDDO
          WRITE(*,FMT='(I3," ZERO Z=  ",2F15.5," ABS(P(Z)= ",E15.5)') &
     &            I,Z(I),ABS(CSVAR)
        ENDDO
        CALL ERROR$MSG('REGULAR STOP AFTER TEST')
        CALL ERROR$STOP('POLYNOM$ZEROSC8')
      END IF
      RETURN
      END SUBROUTINE POLYNOM$ZEROSC8
!
!     ..................................................................
      SUBROUTINE POLYNOM$BINOM(N,B)
!     ******************************************************************
!     ** EVALUATES THE BINOMIAL COEFFICIENTS B(I)=(N-1;I)             **
!     ** USING PASCALS TRIANGLE CONSTRUCTION                          **
!     **                                                              **
!     ** DOES NOT WORK FOR N=68 OR HIGHER BECAUSE OF THE LIMITED      **
!     ** PRECISION OF INTEGER(8) (ASSUMED TO BE 8-BYTE INTEGER)       **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(OUT):: B(N)
      INTEGER(4)            :: IVEC(N)
      INTEGER(4)            :: ISVAR,JSVAR
      INTEGER(4)            :: I,J
!     ****************************************************************** 
     IF(N.GE.68) THEN
        PRINT*,'BINOMIAL COEFFICIENTS TOO LARGE FOR N>68'
        STOP
      END IF
      IVEC(:)=0
      IVEC(1)=1
      DO I=2,N
        ISVAR=IVEC(1)
        DO J=2,I
          JSVAR=IVEC(J)
          IVEC(J)=ISVAR+JSVAR
          ISVAR=JSVAR
        ENDDO
      ENDDO
      DO I=1,N
        B(I)=REAL(IVEC(I))
      ENDDO
      RETURN
      END


