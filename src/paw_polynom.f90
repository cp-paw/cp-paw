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
!     **  FROM DESCRETE POINTS                                          **
!     **     Y(R)=SUM_I=1^NX: COEFF(I)*(X-X0)**(I-1)                    **
!     **  SO THAT                                                       **
!     **     Y(X(I))=Y(I)                                               **
!     **  X0 IS CHOSEN TO BE THE FIRST GRID POINT                       **
!     **                                                                **
!     ********************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      REAL(8)   ,INTENT(out):: X0
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
      REAL(8)   ,INTENT(IN)    :: c1(n1)
      INTEGER(4),INTENT(IN)    :: N2
      REAL(8)   ,INTENT(IN)    :: X2
      REAL(8)   ,INTENT(IN)    :: c2(n2)
      INTEGER(4),INTENT(IN)    :: N3
      REAL(8)   ,INTENT(in)    :: X3
      REAL(8)   ,INTENT(out)   :: c3(n3)
      REAL(8)                  :: x1s,C1s(N1)
      REAL(8)                  :: x2s,C2s(N2)
      REAL(8)                  :: sum
      INTEGER(4)               :: I,J,j1,j2
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
!     **  EVALUATES THE VALUE y oF THE POLYNOMIAL AT POINT X            **
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
!     ..................................................................
      SUBROUTINE POLYNOM$BINOM(N,B)
!     ******************************************************************
!     ** EVALUATES THE BINOMIAL COEFFICIENTS B(I)=(N-1;I)             **
!     ** USING PASCAL'S TRIANGLE CONSTRUCTION                         **
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


