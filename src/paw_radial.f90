!***********************************************************************
!**                                                                   **
!**   TEST ROUTINE FOR RADIAL OBJECT                                  **
!**                                                                   **
!**   NOTE THAT THE INTEGRATE ROUTINE OF THE OLD OBJECT WAS INACCURATE**
!**   BECAUSE OF ONLY QUADRATIC INERPOLATION TO THE ORIGIN            **
!**                                                                   **
!***********************************************************************
!    
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIAL$TEST()
!     **************************************************************************
!     **                                                                      **
!     **   TEST ROUTINE FOR RADIAL OBJECT                                     **
!     **                                                                      **
!     **   NOTE THAT THE INTEGRATE ROUTINE OF THE OLD OBJECT WAS INACCURATE   **
!     **   BECAUSE OF ONLY QUADRATIC INERPOLATION TO THE ORIGIN               **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)           :: GID,GID1,GID2
      INTEGER(4),PARAMETER :: NR=250
      REAL(8)              :: R1=1.056D-4
      REAL(8)              :: DEX=0.05D0
      REAL(8)              :: RARR(NR)
      REAL(8)              :: F(NR)
      REAL(8)              :: DFDR(NR)
      REAL(8)              :: INTF(NR)
      REAL(8)              :: XVAL=1.D0
      REAL(8)              :: FVAL
      REAL(8)              :: INTFVAL
      REAL(8)              :: DERFVAL
      REAL(8)              :: NUMINTF(NR),NUMDFDR(NR)
      REAL(8)              :: VAL
      INTEGER(4)           :: IR
      INTEGER(4)           :: IGRID,IFUNC
      CHARACTER(10)        :: STRING
!     **************************************************************************
!
!     ==========================================================================
!     == DEFINE TWO GRIDS                                                     ==
!     ==========================================================================
      CALL RADIAL$NEW('LOG',GID1)
      CALL RADIAL$SETI4(GID1,'NR',NR)
      CALL RADIAL$SETR8(GID1,'DEX',DEX)
      CALL RADIAL$SETR8(GID1,'R1',R1)
      CALL RADIAL$NEW('SHLOG',GID2)
      CALL RADIAL$SETI4(GID2,'NR',NR)
      CALL RADIAL$SETR8(GID2,'DEX',DEX)
      CALL RADIAL$SETR8(GID2,'R1',R1)
      DO IGRID=1,2
        IF(IGRID.EQ.1) THEN
          GID=GID2
          STRING='SHLOG'
          WRITE(*,FMT='(80("=")/80("="),T20,"  SHIFTED LOG. GRID  "/80("="))')
        ELSE IF(IGRID.EQ.2) THEN
          GID=GID1
          STRING='LOG'
          WRITE(*,FMT='(80("=")/80("="),T20,"  LOG. GRID  "/80("="))')
        END IF
        CALL RADIAL$R(GID,NR,RARR)
!
!       ==================================================================
!       == TEST INTEGRATE AND DERIVATIVE FOR THE LOGARITHMIC GRID      ===
!       ==================================================================
        DO IFUNC=1,3
!
!         ================================================================
!         ==  SELECT TRIAL FUNCTION                                     ==
!         ================================================================
          IF(IFUNC.EQ.1) THEN
            WRITE(*,FMT='(80("=")/80("="),T20,"  GAUSS FUNCTION  "/80("="))')
            DO IR=1,NR
              CALL RADIAL_TESTF1(RARR(IR),F(IR),DFDR(IR),INTF(IR))
            ENDDO
            CALL RADIAL_TESTF1(XVAL,FVAL,DERFVAL,INTFVAL)
          ELSE IF(IFUNC.EQ.2) THEN
            WRITE(*,FMT='(80("=")/80("="),T20,"  SINUS FUNCTION  "/80("="))')
            DO IR=1,NR
              CALL RADIAL_TESTF2(RARR(IR),F(IR),DFDR(IR),INTF(IR))
            ENDDO
            CALL RADIAL_TESTF2(XVAL,FVAL,DERFVAL,INTFVAL)
          ELSE IF(IFUNC.EQ.3) THEN
            WRITE(*,FMT='(80("=")/80("="),T20,"  COSINUS FUNCTION  "/80("="))')
            DO IR=1,NR
              CALL RADIAL_TESTF3(RARR(IR),F(IR),DFDR(IR),INTF(IR))
            ENDDO
            CALL RADIAL_TESTF3(XVAL,FVAL,DERFVAL,INTFVAL)
          END IF
!
!         ================================================================
!         ==  TEST RADIAL$INTEGRATE                                     ==
!         ================================================================
          CALL RADIAL$INTEGRATE(GID,NR,F,NUMINTF)
          PRINT*,'RADIAL$INTEGRATE: MAX. DEV.=', MAXVAL(ABS(NUMINTF-INTF))
!
!         ================================================================
!         ==  TEST RADIAL$DERIVE                                        ==
!         ================================================================
          CALL RADIAL$DERIVE(GID,NR,F,NUMDFDR)
          PRINT*,'RADIAL$DERIVE: MAX. DEV.=', MAXVAL(ABS(NUMDFDR-DFDR))
!
          OPEN(100,FILE='TEST-GAUSS-'//TRIM(STRING)//'.DAT')
          DO IR=1,NR
            WRITE(100,FMT='(6F25.15)')RARR(IR),F(IR),DFDR(IR),NUMDFDR(IR)-DFDR(IR) &
       &                            ,INTF(IR),NUMINTF(IR)-INTF(IR)
          ENDDO
          CLOSE(100)
!
!         ================================================================
!         ==  TEST RADIAL$INTEGRAL                                      ==
!         ================================================================
          CALL RADIAL$INTEGRAL(GID,NR,F,VAL)
          PRINT*,'RADIAL$INTEGRAL: DEV.=',VAL,VAL-INTF(NR)
!
!         ================================================================
!         ==  TEST RADIAL$VALUE                                         ==
!         ================================================================
          CALL RADIAL$VALUE(GID,NR,F,XVAL,VAL)
          PRINT*,'RADIAL$VALUE: VALUE,DEV ',VAL,VAL-FVAL
!
!         ================================================================
!         ==  TEST RADIAL$DERIVATIVE                                    ==
!         ================================================================
          CALL RADIAL$DERIVATIVE(GID,NR,F,XVAL,VAL)
          PRINT*,'RADIAL$DERIVATIVE: VALUE,DEV ',VAL,VAL-DERFVAL
        ENDDO
      ENDDO
      RETURN
      CONTAINS
!  
!       ................................................................
        SUBROUTINE RADIAL_TESTF1(R,F,DFDR,INTF)
        REAL(8),INTENT(IN)  :: R
        REAL(8),INTENT(OUT) :: F
        REAL(8),INTENT(OUT) :: DFDR
        REAL(8),INTENT(OUT) :: INTF
        REAL(8)             :: PI,SVAR
!       ****************************************************************
        PI=4.D0*ATAN(1.D0)
        F=EXP(-R**2)
        DFDR=-2.D0*R*F
        CALL SPECIALFUNCTION$ERF(R,SVAR)
        INTF=0.5D0*SQRT(PI)*SVAR
        RETURN
        END SUBROUTINE RADIAL_TESTF1
!
!       ................................................................
        SUBROUTINE RADIAL_TESTF2(R,F,DFDR,INTF)
        REAL(8),INTENT(IN)  :: R
        REAL(8),INTENT(OUT) :: F
        REAL(8),INTENT(OUT) :: DFDR
        REAL(8),INTENT(OUT) :: INTF
!       ****************************************************************
        F=SIN(2.D0*R)
        DFDR=2.D0*COS(2.D0*R)
        INTF=-0.5D0*(COS(2.D0*R)-1.D0)
        RETURN
        END SUBROUTINE RADIAL_TESTF2
!
!       ................................................................
        SUBROUTINE RADIAL_TESTF3(R,F,DFDR,INTF)
        REAL(8),INTENT(IN)  :: R
        REAL(8),INTENT(OUT) :: F
        REAL(8),INTENT(OUT) :: DFDR
        REAL(8),INTENT(OUT) :: INTF
!       ****************************************************************
        F=COS(2.D0*R)
        DFDR=-2.D0*SIN(2.D0*R)
        INTF=0.5D0*SIN(2.D0*R)
        RETURN
        END SUBROUTINE RADIAL_TESTF3
      END

!***********************************************************************
!**                                                                   **
!**  NAME: RADIAL                                                     **
!**                                                                   **
!**  PURPOSE: OPERATIONS ON A RADIAL LOGARITHMIC GRID                 **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    RADIAL$VALUE(R1,DEX,NR,FUNC,R0,F0)                             **
!**    RADIAL$DERIVE(R1,DEX,NR,F,DFDR)                                **
!**    RADIAL$DERIVATIVE(R1,DEX,NR,F,R0,DFDR0)                        **
!**    RADIAL$INTEGRATE(R1,DEX,NR,F,FINT)                             **
!**    RADIAL$INTEGRAL(R1,DEX,NR,F,R0,FINT0)                          **
!**    RADIAL$FFT(ID,R1,DEX,NR,F1,G1,F2)                              **
!**    RADIAL$POISSON(R1,DEX,NR,RHO,V)                                **
!**    RADIAL$MOMENT(R1,DEX,NR,L,RHO,QLM)                             **
!**    BESSELTRANSFORM$CLEAR                                          **
!**    BESSELTRANSFORM(L,NP,R1,G1,DEX,F,G,DISC)                       **
!**    BESSOV(R1,DEX,NR,F,L,G,RMAX,RES)                               **
!**                                                                   **
!**  REMARKS:                                                         **
!**    1) DOES NOT CONTAIN PERMANET DATA                              **
!**    2) NOT YET IMPLEMENTED                                         **
!**    3) USES BESSELTRANSFORM_MODULE                                 **
!**                                                                   **
!***********************************************************************
!***********************************************************************
!**  RADIAL TOOLBOX                                                   **
!**                                                                   **
!**  CONTAINS HELPER ROUTINES THAT ARE USED ALSO BY THE               **
!**  CONTAINED OBJECTS SUCH AS LOGRADIAL AND SHLOGRADIAL              **
!**                                                                   **
!**                                                                   **
!**                                                                   **
!**                                                                   **
!***********************************************************************
!***********************************************************************
!***********************************************************************

!     ....................................................................
      SUBROUTINE RADIAL_POLYNOMIALCOEFFICIENTS(NP,RI_,FI_,R0,CN)
!     **                                                              **
!     **  OBTAINS THE COEFFICIENTS CN FOR A POWERSERIES EXPANSION OF  **
!     **  ORDER NP THROUGH NP DATA POINTS. R0 IS THE EXPANSION POINT  **
!     **  THE INTERPOLATED FUNCTION IS GIVEN BY                       **
!     **     F(R)=\SUM_{J=0}^{N-1} C_J (R-R_0)^{J}                    **
!     **                                                              **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NP
      REAL(8)   ,INTENT(IN) :: RI_(NP)
      REAL(8)   ,INTENT(IN) :: FI_(NP)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: CN(NP)
      REAL(8)               :: FI(NP)
      REAL(8)               :: RI(NP)
      REAL(8)               :: DCN(NP)
      REAL(8)               :: DCN1(NP)
      REAL(8)               :: SVAR,VAL,A,B
      INTEGER(4)            :: I,J,K
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
!     ****************************************************************
      FI(:)=FI_(:)
      RI(:)=RI_(:)-R0
      CN(:)=0.D0
      DO I=1,NP
!       == CONSTRUCT POLYNOMIAL OF ORDER I-1
        DCN(:)=0.D0
        DCN(1)=FI(I)
        DO J=1,I-1
          A=1.D0/(RI(I)-RI(J))
          B=-A*RI(J)
          DCN1=0.D0
          DO K=1,J
            DCN1(K+1)=DCN1(K+1)+A*DCN(K)
            DCN1(K)=DCN1(K)+B*DCN(K)
          ENDDO
          DCN(:)=DCN1(:)
        ENDDO
        CN(:)=CN(:)+DCN(:)
!       == SUBTRACT POLYNOMIAL FROM DATA
        DO J=1,NP
          VAL=0.D0
          SVAR=1.D0
          DO K=1,I
            VAL=VAL+DCN(K)*SVAR
            SVAR=SVAR*RI(J)
          ENDDO
          FI(J)=FI(J)-VAL
        ENDDO
      ENDDO
!
!     ================================================================
!     == TEST ACCURACY OF THE POLYNOMIAL                            ==
!     ================================================================
      IF(TPR) THEN
        DO I=1,NP
          VAL=0.D0
          DO J=1,NP
            VAL=VAL+CN(J)*RI(I)**(J-1)
          ENDDO
          PRINT*,'VAL ',VAL,FI_(I),VAL-FI_(I)
        ENDDO
      END IF
      RETURN
      END SUBROUTINE RADIAL_POLYNOMIALCOEFFICIENTS
!
!     ....................................................................
      SUBROUTINE RADIAL_POLYNOMIALVALUE(NP,RI,FI_,R0,F0)
!     **                                                              **
!     **  DETERMINES THE VALUE F0 AT R0 FROM A POLYNOM OF ORDER NP    **
!     **  THAT PASSES THROUGH NP DATA POINTS (RI,FI_)                 **
!     **                                                              **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NP
      REAL(8)   ,INTENT(IN) :: RI(NP)
      REAL(8)   ,INTENT(IN) :: FI_(NP)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      REAL(8)               :: FI(NP)
      REAL(8)               :: SVAR
      INTEGER(4)            :: I,J,IP
!     ****************************************************************
      FI(:)=FI_(:)
      F0=0.D0
      DO I=1,NP
!       == CONSTRUCT A POLYNOMIAL WITH VALUE 1 AT THE ITH POINT
!       == A ZERO ON ALL POINTS LEFT TO THE ITH POINT
        SVAR=1.D0
        DO J=1,I-1
          SVAR=SVAR*(R0-RI(J))/(RI(I)-RI(J))
        ENDDO
!       == SVAR IS THE VALUE OF THE POLYNOMIAL AT R0, 
!       == THAT HAS A ZERO ON THE FIRST I-1 GRID POINTS AND THAT HAS
!       == THE VALUE 1 AT THE I-TH GRID POINT
        F0=F0+SVAR*FI(I)
!       == NOW SUBTRACT THAT POLYNOMIAL FROM ALL GRID POINTS FROM I+1 TO NP
        DO IP=I+1,NP
          SVAR=1.D0
          DO J=1,I-1
            SVAR=SVAR*(RI(IP)-RI(J))/(RI(I)-RI(J))
          ENDDO
          FI(IP)=FI(IP)-SVAR*FI(I)
        ENDDO
        FI(I)=0.D0
      ENDDO
      RETURN
      END SUBROUTINE RADIAL_POLYNOMIALVALUE
!
!     ..................................................................
      SUBROUTINE RADIAL_POLYNOMIALDERIVATIVE(NP,RI,FI_,R0,DF0)
!     **                                                              **
!     **  DETERMINES THE GRADIENT DF0 AT R0 FROM A POLYNOM OF ORDER NP**
!     **  THAT PASSES THROUGH NP DATA POINTS (RI,FI_)                 **
!     **                                                              **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NP
      REAL(8)   ,INTENT(IN) :: RI(NP)
      REAL(8)   ,INTENT(IN) :: FI_(NP)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: DF0
      REAL(8)               :: F0(NP)
      REAL(8)               :: FI(NP)
      REAL(8)               :: SVAR,DSVAR,FAC
      INTEGER(4)            :: I,J,IP
!     ****************************************************************
      FI(:)=FI_(:)
      F0=0.D0
      DF0=0.D0
      DO I=1,NP
        SVAR=1.D0
        DSVAR=0.D0
        DO J=1,I-1
!         == DO NOT CHANGE THE ORDER OF THE NEXT TWO STATEMENTS ======
!         == SVAR=PROD_J=1^I-1 (R-RJ)/(RI-RJ)
!         == SVAR=D/DR SVAR=PROD_J=1^I-1 [1/(RI-RJ)+(R-RJ)/(RI-RJ)*D/DR)
          FAC=1.D0/(RI(I)-RI(J))
          DSVAR=(DSVAR*(R0-RI(J))+SVAR)*FAC
          SVAR=SVAR*(R0-RI(J))*FAC
        ENDDO
        F0=F0+SVAR*FI(I)
        DF0=DF0+DSVAR*FI(I)
        DO IP=I+1,NP
          SVAR=1.D0
          DO J=1,I-1
            SVAR=SVAR*(RI(IP)-RI(J))/(RI(I)-RI(J))
          ENDDO
          FI(IP)=FI(IP)-SVAR*FI(I)
        ENDDO
        FI(I)=0.D0
      ENDDO
      RETURN
      END SUBROUTINE RADIAL_POLYNOMIALDERIVATIVE
!
!     .....................................................INTRAD.......
      SUBROUTINE RADIAL_INTEGRATEEQUISPACED(NX,F,G)
!     **                                                              **
!     **  INTEGRATES THE FUNCTION F(X) ON AN EQUISPACED GRID          **
!     **  X(I)=I WITH UNIT STEP SIZE FROM X=1 TO X(I) TO OBTAIN G(X)  **
!     **                                                              **
!     ** INPUT :                                                      **
!     **   NX    NUMBER OF GRID POINTS                                **
!     **   F            FUNCTION TO BE INTEGRATED                     **
!     ** OUTPUT :                                                     **
!     **   G            INTEGRAL OF F (FROM THE FIRST GRID POINT)     **
!     **  REMARKS :                                                   **
!     **  (SEE: NUMERICAL RECIPES EQ: 4.1.14;                         **
!     **  INTEGRATES POLYNOMIALS UP TO 3. ORDER EXACTLY)              **
!     **                                                              **
!     **  THE FIRST GRID POINTS ARE CALCULATED DIRECTLY FROM A CUBE   **
!     **  POLYNOMIAL INTERPOLATION OF THE FIRST FOUR POINTS           **
!     **                                                              **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX
      REAL(8)   ,INTENT(IN) :: F(NX)
      REAL(8)   ,INTENT(OUT):: G(NX)
      REAL(8)   ,PARAMETER  :: C0=-31.D0/48.D0
      REAL(8)   ,PARAMETER  :: C1=+11.D0/48.D0
      REAL(8)   ,PARAMETER  :: C2= -5.D0/48.D0
      REAL(8)   ,PARAMETER  :: C3= +1.D0/48.D0
      INTEGER(4)            :: I
      REAL(8)               :: X,A,B,C,D
!     ******************************************************************
      IF(NX.LT.10) THEN
         CALL ERROR$MSG('NUMBER OF MESHPOINTS SMALLLER THAN 10')
         CALL ERROR$STOP('RADIAL_INTEGRATEEQUISPACED')
      END IF
!     ==================================================================
!     ==  SUMMATION                                                   ==
!     ==================================================================
      G(1)=C0*F(1)+C1*F(2)+C2*F(3)+C3*F(4)+F(1)
      DO I=2,NX
        G(I)=G(I-1)+F(I)
      ENDDO
!     ==================================================================
!     ==  FIX ENDPOINT                                                ==
!     ==================================================================
      DO I=4,NX
        G(I)=G(I)+C0*F(I)+C1*F(I-1)+C2*F(I-2)+C3*F(I-3)
      ENDDO
!     ==================================================================
!     ==  FIX FIRST THREE GRID POINTS                                 ==
!     ==================================================================
!     ==  F(X)=A+B*X+C*X^2+D*X^3
      A=(+24.D0*F(1)-36.D0*F(2)+24.D0*F(3) -6.D0*F(4))/6.D0
      B=(-26.D0*F(1)+57.D0*F(2)-42.D0*F(3)+11.D0*F(4))/6.D0
      C=(  9.D0*F(1)-24.D0*F(2)+21.D0*F(3) -6.D0*F(4))/6.D0
      D=(      -F(1) +3.D0*F(2) -3.D0*F(3)      +F(4))/6.D0
      DO I=1,3
        X=REAL(I,KIND=8)
        G(I)=A*(X-1.D0)+0.5D0*B*(X**2-1.D0) &
     &      +C*(X**3-1.D0)/3.D0+0.25D0*D*(X**4-1.D0)
      ENDDO
      RETURN
      END
!
!     .......................................................SPHLSD.....
      SUBROUTINE RADIAL_DERIVEEQUISPACED(NX,F,G)
!     **                                                              **
!     **  TAKES THE RADIAL DERIVATIVE OF A FUNCTION F GIVEN ON A      **
!     **  LOGARITHMIC GRID                                            **
!     **                                                              **
!     **  BASED ON A INTERPOLATION BY A FOURTH ORDER POLYNOM          **
!     **                                                              **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: NX
      REAL(8)    ,INTENT(IN) :: F(NX)
      REAL(8)    ,INTENT(OUT):: G(NX)
      REAL(8)    ,PARAMETER  :: C21=- 2.D0  
      REAL(8)    ,PARAMETER  :: C22=-14.D0  
      REAL(8)    ,PARAMETER  :: C23=+24.D0  
      REAL(8)    ,PARAMETER  :: C24=-10.D0  
      REAL(8)    ,PARAMETER  :: C25=+ 2.D0  
      REAL(8)    ,PARAMETER  :: C11=-11.D0  
      REAL(8)    ,PARAMETER  :: C12=+18.D0  
      REAL(8)    ,PARAMETER  :: C13=- 9.D0  
      REAL(8)    ,PARAMETER  :: C14=+ 2.D0  
      REAL(8)    ,PARAMETER  :: CI1=1.D0/12.D0
      REAL(8)    ,PARAMETER  :: CI2=-2.D0/3.D0
      INTEGER(4)             :: I
!     ******************************************************************
!     ==================================================================
!     == FORM DERIVATIVE ON THE EQUI SPACED X-GRID                    ==
!     ==================================================================
!CHECK FOUR END POINTS!!
      G(1)=(C11*F(1)+C12*F(2)+C13*F(3)+C14*F(4))/6.D0
      G(2)=(C21*F(1)+C22*F(2)+C23*F(3)+C24*F(4)+C25*F(5))/12.D0
      DO I=3,NX-2
        G(I)=CI1*F(I-2)+CI2*F(I-1)-CI2*F(I+1)-CI1*F(I+2)
      ENDDO
      G(NX-1)=-(C25*F(NX-4)+C24*F(NX-3) &
     &         +C23*F(NX-2)+C22*F(NX-1)+C21*F(NX))/12.D0
      G(NX)=-(C14*F(NX-3)+C13*F(NX-2)+C12*F(NX-1)+C11*F(NX))/6.D0
      RETURN
      END
!
!     .......................................................SPHLSD.....
      SUBROUTINE RADIAL_DGLEQUISPACED(IDIR,NX,A,B,C,D,F)
!     **                                                              **
!     **  SOLVES THE DGL SECOND ORDER ON AN EQUISPACED GRID           **
!     **                                                              **
!     **  [A(X)\PARTIAL^2_X+B(X)\PARTIAL_X+C(X)]F(X)=D(X)             **
!     **                                                              **
!     **  FOR IDIR=1, F(1) AND F(2) MUST BE SUPPLIED ON INPUT         **
!     **  FOR IDIR=-1, F(NX) AND F(NX-1) MUST BE SUPPLIED ON INPUT    **
!     **                                                              **
!     **  CAUTION! THERE IS NO CATCH AGAINST OVERFLOW                 **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: IDIR  ! DIRECTION
      INTEGER(4) ,INTENT(IN) :: NX
      REAL(8)    ,INTENT(IN) :: A(NX)
      REAL(8)    ,INTENT(IN) :: B(NX)
      REAL(8)    ,INTENT(IN) :: C(NX)
      REAL(8)    ,INTENT(IN) :: D(NX)
      REAL(8)    ,INTENT(INOUT):: F(NX)
      REAL(8)                :: AP(NX),A0(NX),AM(NX)
      INTEGER(4)             :: I
!     ******************************************************************
!     == AP*F(+) + A0*F(0) + AM*F(-) = D
      AP(:)=A(:)+0.5D0*B(:)
      A0(:)=-2.D0*A(:)+C(:)
      AM(:)=A(:)-0.5D0*B(:)
      IF(IDIR.GE.0) THEN
        DO I=2,NX-1
          F(I+1)=( -(A0(I)+AM(I))*F(I) +AM(I)*(F(I)-F(I-1)) +D(I) )/AP(I)
        ENDDO
      ELSE IF(IDIR.LT.0) THEN
        DO I=NX-1,2,-1
          F(I-1)=( -(A0(I)+AP(I))*F(I) +AP(I)*(F(I)-F(I+1)) +D(I) )/AM(I)
        ENDDO
      ELSE
         CALL ERROR$MSG('INVALID VALUE OF IDIR')
         CALL ERROR$STOP('RADIAL_DGLEQUISPACED')
      END IF
      RETURN
      END
!
!     .......................................................SPHLSD.....
      SUBROUTINE RADIAL_DGLEQUISPACEDGEN(NX,NF,I1,I2,A,B,C,D,F)
!     **                                                              **
!     **  SOLVES THE DGL SECOND ORDER ON AN EQUISPACED GRID           **
!     **                                                              **
!     **  [A(X)\PARTIAL^2_X+B(X)\PARTIAL_X+C(X)]F(X)=D(X)             **
!     **                                                              **
!     **  FOR IDIR=1, F(1) AND F(2) MUST BE SUPPLIED ON INPUT         **
!     **  FOR IDIR=-1, F(NX) AND F(NX-1) MUST BE SUPPLIED ON INPUT    **
!     **                                                              **
!     **  CAUTION! THERE IS NO CATCH AGAINST OVERFLOW                 **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: NX
      INTEGER(4) ,INTENT(IN) :: NF
      INTEGER(4) ,INTENT(IN) :: I1
      INTEGER(4) ,INTENT(IN) :: I2
      REAL(8)    ,INTENT(IN) :: A(NX)
      REAL(8)    ,INTENT(IN) :: B(NX)
      REAL(8)    ,INTENT(IN) :: C(NX,NF,NF)
      REAL(8)    ,INTENT(IN) :: D(NX,NF)
      REAL(8)    ,INTENT(INOUT):: F(NX,NF)
      REAL(8)                :: AP(NX),A0(NX),AM(NX)
      INTEGER(4)             :: I
      INTEGER(4)             :: IDIR
!     ******************************************************************
      IF(MIN(I1,I2).LT.1.OR.MAX(I1,I2).GT.NX) THEN
        CALL ERROR$MSG('INTEGRATION BOUNDS OUT OF RANGE')
        CALL ERROR$I4VAL('I1',I1)
        CALL ERROR$I4VAL('I2',I2)
        CALL ERROR$I4VAL('NX',NX)
        CALL ERROR$STOP('RADIAL_DGLEQUISPACEDGENC')
      END IF
!     == IDIR IS THE DIRECTION OF THE INTEGRATION =======================
      IDIR=1
      IF(I2.LT.I1) IDIR=-1
!     == AP*F(+) + A0*F(0) + AM*F(-) = D
      AP(:)=A(:)+0.5D0*B(:)
      A0(:)=-2.D0*A(:)
      AM(:)=A(:)-0.5D0*B(:)
      IF(IDIR.GE.0) THEN
        DO I=I1+1,I2-1
          F(I+1,:)=( -(A0(I)+AM(I))*F(I,:) +AM(I)*(F(I,:)-F(I-1,:)) &
     &               -MATMUL(C(I,:,:),F(I,:)) +D(I,:) )/AP(I)
        ENDDO
        F(:I1-1,:)=0.D0
        F(I2+1:,:)=0.D0
      ELSE IF(IDIR.LT.0) THEN
        DO I=I1-1,I2+1,-1
          F(I-1,:)=( -(A0(I)+AP(I))*F(I,:) +AP(I)*(F(I,:)-F(I+1,:)) &
     &               -MATMUL(C(I,:,:),F(I,:)) +D(I,:) )/AM(I)
        ENDDO
        F(:I2-1,:)=0.D0
        F(I1+1:,:)=0.D0
      ELSE
         CALL ERROR$MSG('INVALID VALUE OF IDIR')
         CALL ERROR$STOP('RADIAL_DGLEQUISPACED')
      END IF
      RETURN
      END
!
!     .......................................................SPHLSD.....
      SUBROUTINE RADIAL_DGLEQUISPACEDGENC(NX,NF,I1,I2,A,B,C,D,F)
!     **                                                              **
!     **  SOLVES THE DGL SECOND ORDER ON AN EQUISPACED GRID           **
!     **                                                              **
!     **  [A(X)\PARTIAL^2_X+B(X)\PARTIAL_X+C(X)]F(X)=D(X)             **
!     **                                                              **
!     **  FOR I2>I1, F(I1) AND F(I1+1) MUST BE SUPPLIED ON INPUT      **
!     **  FOR I2<I1, F(I2-1) AND F(I2) MUST BE SUPPLIED ON INPUT      **
!     **                                                              **
!     **  CAUTION! THERE IS NO CATCH AGAINST OVERFLOW                 **
!     **                                                              **
!     **                                                              **
!     **                                                              **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: I1
      INTEGER(4) ,INTENT(IN) :: I2
      INTEGER(4) ,INTENT(IN) :: NX
      INTEGER(4) ,INTENT(IN) :: NF
      REAL(8)    ,INTENT(IN) :: A(NX)
      REAL(8)    ,INTENT(IN) :: B(NX)
      COMPLEX(8) ,INTENT(IN) :: C(NX,NF,NF)
      COMPLEX(8) ,INTENT(IN) :: D(NX,NF)
      COMPLEX(8) ,INTENT(INOUT):: F(NX,NF)
      REAL(8)                :: AP(NX),A0(NX),AM(NX)
      COMPLEX(8)             :: C1(NF,NF,NX)
      COMPLEX(8)             :: F1(NF,NX)
      COMPLEX(8)             :: D1(NF,NX)
      INTEGER(4)             :: I
      INTEGER(4)             :: IMIN,IMAX
      INTEGER(4)             :: IDIR
!     ******************************************************************
      IF(MIN(I1,I2).LT.1.OR.MAX(I1,I2).GT.NX) THEN
        CALL ERROR$MSG('INTEGRATION BOUNDS OUT OF RANGE')
        CALL ERROR$I4VAL('I1',I1)
        CALL ERROR$I4VAL('I2',I2)
        CALL ERROR$I4VAL('NX',NX)
        CALL ERROR$STOP('RADIAL_DGLEQUISPACEDGENC')
      END IF
      IMIN=MIN(I1,I2)
      IMAX=MAX(I1,I2)
!     == IDIR IS THE DIRECTION OF THE INTEGRATION =======================
      IDIR=1
      IF(I2.LT.I1) IDIR=-1
!     == AP*F(+) + A0*F(0) + AM*F(-) = D
      AP(:)=A(:)+0.5D0*B(:)
      A0(:)=-2.D0*A(:)
      AM(:)=A(:)-0.5D0*B(:)
!     == REARRANGE C TO AVOID CASH PROBLEMS
      DO I=IMIN,IMAX
         C1(:,:,I)=C(I,:,:)
         D1(:,I)=D(I,:)
         F1(:,I)=F(I,:)
      ENDDO
      IF(IDIR.GE.0) THEN
        DO I=I1+1,I2-1
!!$          F(I+1,:)=( -(A0(I)+AM(I))*F(I,:) +AM(I)*(F(I,:)-F(I-1,:)) &
!!$    &               -MATMUL(C(I,:,:),F(I,:)) +D(I,:) )/AP(I)
!         == VERSION 2
          F1(:,I+1)=( -(A0(I)+AM(I))*F1(:,I) +AM(I)*(F1(:,I)-F1(:,I-1)) &
     &               -MATMUL(C1(:,:,I),F1(:,I)) +D1(:,I) )/AP(I)
!         == VERSION 3
        ENDDO
      ELSE IF(IDIR.LT.0) THEN
        DO I=I1-1,I2+1,-1
!           F(I-1,:)=( -(A0(I)+AP(I))*F(I,:) +AP(I)*(F(I,:)-F(I+1,:)) &
!    &               -MATMUL(C(I,:,:),F(I,:)) +D(I,:) )/AM(I)
          F1(:,I-1)=( -(A0(I)+AP(I))*F1(:,I) +AP(I)*(F1(:,I)-F1(:,I+1)) &
     &               -MATMUL(C1(:,:,I),F1(:,I)) +D1(:,I) )/AM(I)
        ENDDO
      ELSE
        CALL ERROR$MSG('INVALID VALUE OF IDIR')
        CALL ERROR$STOP('RADIAL_DGLEQUISPACEDGENC')
      END IF
      F(:IMIN-1,:)=(0.D0,0.D0)
      F(IMAX+1:,:)=(0.D0,0.D0)
      RETURN
      END
!
!     .......................................................SPHLSD.....
      SUBROUTINE RADIAL_VERLETD1EQUISPACED(NX,F,DFDX)
!     **                                                              **
!     **                                                              **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: NX
      REAL(8)    ,INTENT(IN) :: F(NX)
      REAL(8)    ,INTENT(OUT):: DFDX(NX)
      INTEGER(4)             :: I
!     ******************************************************************
      DFDX(1)=-1.5D0*F(1)+2.D0*F(2)-0.5D0*F(3)
      DO I=2,NX-1
        DFDX(I)=0.5D0*(F(I+1)-F(I-1))
      ENDDO
      DFDX(NX)=-1.5D0*F(NX)+2.D0*F(NX-1)-0.5D0*F(NX-2)
      RETURN
      END
!
!     .......................................................SPHLSD.....
      SUBROUTINE RADIAL_VERLETD2EQUISPACED(NX,F,D2FDX2)
!     **                                                              **
!     **                                                              **
!     *****************************************************************
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: NX
      REAL(8)    ,INTENT(IN) :: F(NX)
      REAL(8)    ,INTENT(OUT):: D2FDX2(NX)
      INTEGER(4)             :: I
!     ******************************************************************
      D2FDX2(1)=F(1)-2.D0*F(2)+F(3)
      DO I=2,NX-1
        D2FDX2(I)=F(I+1)-2.D0*F(I)+F(I-1)
      ENDDO
!      D2FDX2(NX)=F(NX)-2.D0*F(NX-1)+F(NX-2)
!     == EXTRAPOLATE SECOND DERIVATIVE FROM PREVIOUS TWO GRID POINTS
      D2FDX2(NX)=2.D0*D2FDX2(NX-1)-D2FDX2(NX-2)
      RETURN
      END
!***********************************************************************
!***********************************************************************
!**                                                                   **
!**  RADIAL INTERFACE                                                 **
!**                                                                   **
!**  RESOLVES THE GRID ID GID AND PASSES THE REQUEST ON TO THE        **
!**  LOWER LEVEL OBJECTS SUCH AS LOGRADIAL AND SHLOGRADIAL            **
!**                                                                   **
!*********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
!.......................................................................
MODULE RADIAL_MODULE
INTEGER(4) ,PARAMETER   :: NGRIDTYPES=2
INTEGER(4)              :: NGIDX=0
INTEGER(4)              :: NGID=0 ! #(GID)
INTEGER(4) ,ALLOCATABLE :: GRIDARRAY(:,:) !GRIDTYPE, GRIDNUMBER
CHARACTER(8) ,PARAMETER :: GRIDTYPE(NGRIDTYPES)=(/'LOG   ','SHLOG '/)
END MODULE RADIAL_MODULE
!
!     .................................................................
      SUBROUTINE RADIAL$NEW(TYPEID,GID)
!     *****************************************************************
!     ** DEFINE A NEW GRID  AND RETURN THE GRID ID "GID"             **
!     *****************************************************************
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      USE RADIAL_MODULE
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)  :: TYPEID
      INTEGER(4)  ,INTENT(OUT) :: GID
      INTEGER(4)  ,PARAMETER   :: NGIDBLOCKSIZE=20
      INTEGER(4)  ,ALLOCATABLE :: NEWARRAY(:,:)
      INTEGER(4)               :: TYPE,GIDS
      INTEGER(4)               :: NEWNGIDX
!     *****************************************************************
!      
!     ==================================================================
!     == REALLOCATE GRIDARRAY, IF IT IS TO SMALL                      ==
!     ==================================================================
      IF(NGID.GE.NGIDX) THEN
!       == SAVE CONTENTS OF GRIDARRAY AND DEALLOCATE IT
        IF(ALLOCATED(GRIDARRAY)) THEN
          ALLOCATE(NEWARRAY(2,NGIDX))
          NEWARRAY(:,:)=GRIDARRAY(:,:)
          DEALLOCATE(GRIDARRAY)
        END IF
!       == ALLOCATE GRIDARRAY WITH LARGER SIZE ========================
        NEWNGIDX=NGIDX+NGIDBLOCKSIZE
        ALLOCATE(GRIDARRAY(2,NEWNGIDX))
        GRIDARRAY(:,:)=0
!       == COPY CONTENTS SAVED FROM GRIDARRAY BACK ====================
        IF(ALLOCATED(NEWARRAY)) THEN
          GRIDARRAY(:,:NGIDX)=NEWARRAY(:,:NGIDX)
          DEALLOCATE(NEWARRAY)
        END IF
        NGIDX=NEWNGIDX
      END IF
!
!     ==================================================================
!     == DETERMINE GRID TYPE AND MAKE NEW GRID                        ==
!     ==================================================================
      IF(TYPEID.EQ.GRIDTYPE(1)) THEN 
        TYPE=1
        CALL LOGRADIAL_NEW(GIDS)
      ELSE IF(TYPEID.EQ.GRIDTYPE(2)) THEN 
        TYPE=2
        CALL SHLOGRADIAL_NEW(GIDS)
      ELSE
        CALL ERROR$MSG('ERROR: GRIDTYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$NEW')
      END IF
!
!     =================================================================
!     ==  NOW SAVE NEW GRID                                          ==
!     =================================================================
      NGID=NGID+1
      GID=NGID
      GRIDARRAY(1,GID)=TYPE
      GRIDARRAY(2,GID)=GIDS
      RETURN
      END SUBROUTINE RADIAL$NEW
!
!     .................................................................
      SUBROUTINE RADIAL_RESOLVE(GID,GIDS,TYPE)
!     *****************************************************************
!     **  RESOLVE GRID-ID "GID" AND DETERMINE  GRID TYPE AND GRID-ID **
!     **  FOR THAT GRID TYPE                                         **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      USE RADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)  :: GID
      INTEGER(4) ,INTENT(OUT) :: GIDS
      INTEGER(4) ,INTENT(OUT) :: TYPE
!     ***************************************************************
      IF(GID.GT.NGID) THEN
        CALL ERROR$MSG('GRID ID OUT OF RANGE')
        CALL ERROR$I4VAL('GID',GID)
        CALL ERROR$I4VAL('NGID',NGID)
        CALL ERROR$STOP('RADIAL$RESOLVE')
      END IF
      TYPE=GRIDARRAY(1,GID)
      GIDS=GRIDARRAY(2,GID)
      RETURN
      END SUBROUTINE RADIAL_RESOLVE
!
!     .................................................................
      SUBROUTINE RADIAL$REPORT(NFIL)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      USE RADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: NFIL
      INTEGER(4)              :: GIDS,TYPE,I,NR
      REAL(8)                 :: R1,DEX
!     ***************************************************************
      DO I=1,NGID
        GIDS=GRIDARRAY(2,I)
        TYPE=GRIDARRAY(1,I)
        CALL RADIAL$GETI4(I,'NR',NR)
        CALL RADIAL$GETR8(I,'DEX',DEX)
        CALL RADIAL$GETR8(I,'R1',R1)
        WRITE(NFIL,FMT=*)I,GRIDTYPE(TYPE),GIDS,R1,DEX,NR
      ENDDO
      RETURN
      END SUBROUTINE RADIAL$REPORT
!
!     .................................................................
      SUBROUTINE RADIAL$SETR8(GID,ID,VAL)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      USE RADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(IN) :: VAL
      INTEGER(4)              :: GIDS,TYPE
!     ***************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN 
        CALL LOGRADIAL$SETR8(GIDS,ID,VAL)
      ELSE IF(TYPE.EQ.2) THEN 
        CALL SHLOGRADIAL$SETR8(GIDS,ID,VAL)
      END IF
      RETURN
      END SUBROUTINE RADIAL$SETR8
!
!     .................................................................
      SUBROUTINE RADIAL$GETR8(GID,ID,VAL)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      USE RADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      CHARACTER(*),INTENT(IN) :: ID
      REAL(8)     ,INTENT(OUT):: VAL
      INTEGER(4)              :: GIDS,TYPE
!     ***************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN 
        CALL LOGRADIAL$GETR8(GIDS,ID,VAL)
      ELSE IF(TYPE.EQ.2) THEN 
        CALL SHLOGRADIAL$GETR8(GIDS,ID,VAL)
      END IF
      RETURN
      END SUBROUTINE RADIAL$GETR8
!
!     .................................................................
      SUBROUTINE RADIAL$SETI4(GID,ID,VAL)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      USE RADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(IN) :: VAL
      INTEGER(4)              :: GIDS,TYPE
!     ***************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN 
        CALL LOGRADIAL$SETI4(GIDS,ID,VAL)
      ELSE IF(TYPE.EQ.2) THEN 
        CALL SHLOGRADIAL$SETI4(GIDS,ID,VAL)
      END IF
      RETURN
      END SUBROUTINE RADIAL$SETI4
!
!     .................................................................
      SUBROUTINE RADIAL$GETI4(GID,ID,VAL)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      USE RADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN) :: GID
      CHARACTER(*),INTENT(IN) :: ID
      INTEGER(4)  ,INTENT(OUT):: VAL
      INTEGER(4)              :: GIDS,TYPE
!     ***************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN 
        CALL LOGRADIAL$GETI4(GIDS,ID,VAL)
      ELSE IF(TYPE.EQ.2) THEN 
        CALL SHLOGRADIAL$GETI4(GIDS,ID,VAL)
      END IF
      RETURN
      END SUBROUTINE RADIAL$GETI4
!
!     .......................................................SPHLSD.....
      SUBROUTINE RADIAL$CHANGEGRID(GID1,NR1,F,GID2,NR2,G)
!     **                                                              **
!     **                                                              **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: GID1
      INTEGER(4) ,INTENT(IN) :: GID2
      INTEGER(4) ,INTENT(IN) :: NR1
      INTEGER(4) ,INTENT(IN) :: NR2
      REAL(8)    ,INTENT(IN) :: F(NR1)
      REAL(8)    ,INTENT(OUT):: G(NR2)
      INTEGER(4)             :: IR
      REAL(8)                :: R2(NR2)
      REAL(8)                :: R1(NR1)
!     ******************************************************************
      CALL RADIAL$R(GID1,NR1,R1)
      CALL RADIAL$R(GID2,NR2,R2)
      DO IR=1,NR2
!       == AVOID EXTRAPOLATION BEYOND END OF GRID
        IF(R2(IR).GT.R1(NR1)) THEN
          G(IR:)=0.D0
          EXIT
        END IF
!       == INTERPOLATE
        CALL RADIAL$VALUE(GID1,NR1,F,R2(IR),G(IR))
      ENDDO       
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE RADIAL$R(GID,NR,R)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(OUT):: R(NR)
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
INTEGER(4) :: IR
REAL(8)    :: XEXP,RI
!     ******************************************************************
IF(GID.EQ.0) THEN
  IF(NR.NE.250) STOP 'ERROR IN RADIAL$R/ GID=0'
  XEXP=EXP(0.05D0)
  RI=1.056D-4/XEXP
  DO IR=1,NR
    RI=RI*XEXP
    R(IR)=RI
  ENDDO
  RETURN
END IF
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$R(GIDS,NR,R)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$R(GIDS,NR,R) 
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$R')
      END IF  
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$VALUE(GID,NR,F,R0,F0)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
USE RADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F(NR)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$VALUE(GIDS,NR,F,R0,F0)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$VALUE(GIDS,NR,F,R0,F0) 
      ELSE
PRINT*,'GRIDARRAY ',GRIDARRAY(:,:NGID)
PRINT*,'GIDS ',GIDS
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$I4VAL('GID',GID)
        CALL ERROR$I4VAL('TYPE',TYPE)
        CALL ERROR$STOP('RADIAL$VALUE')
      END IF  
      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIAL$XOFR(GID,R0,X0)
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2007 ************
      USE RADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: X0
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$XOFR(GIDS,R0,X0)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$XOFR(GIDS,R0,X0) 
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$I4VAL('GID',GID)
        CALL ERROR$I4VAL('TYPE',TYPE)
        CALL ERROR$STOP('RADIAL$XOFR')
      END IF  
      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE RADIAL$INTEGRATE(GID,NR,F,G)
!     **************************************************************************
!     **                                                                      **
!     **  determines the indefinite integral (antiderivative) of F            **
!     **                                                                      **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F(NR)
      REAL(8)   ,INTENT(OUT):: G(NR)
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$INTEGRATE(GIDS,NR,F,G)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$INTEGRATE(GIDS,NR,F,G)
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$INTEGRATE')
      END IF  
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$INTEGRAL(GID,NR,F,VAL)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F(NR)
      REAL(8)   ,INTENT(OUT):: VAL
      REAL(8)               :: G(NR)
!     ******************************************************************
      CALL RADIAL$INTEGRATE(GID,NR,F,G)
      VAL=G(NR)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$DERIVE(GID,NR,F,G)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F(NR)
      REAL(8)   ,INTENT(OUT):: G(NR)
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$DERIVE(GIDS,NR,F,G)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$DERIVE(GIDS,NR,F,G)
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$DERIVE')
      END IF  
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$DERIVATIVE(GID,NR,F,R0,F0)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F(NR)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$DERIVATIVE(GIDS,NR,F,R0,F0)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$DERIVATIVE(GIDS,NR,F,R0,F0) 
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$DERIVATIVE')
      END IF  
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$DGL(GID,IDIR,NR,A,B,C,D,F)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: IDIR
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: A(NR)
      REAL(8)   ,INTENT(IN) :: B(NR)
      REAL(8)   ,INTENT(IN) :: C(NR)
      REAL(8)   ,INTENT(IN) :: D(NR)
      REAL(8)   ,INTENT(INOUT):: F(NR)
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$DGL(GIDS,IDIR,NR,A,B,C,D,F)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$DGL(GIDS,IDIR,NR,A,B,C,D,F)
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$DGL')
      END IF  
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$DGLGEN(GID,NR,NF,I1,I2,A,B,C,D,F)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NF
      INTEGER(4),INTENT(IN) :: I1
      INTEGER(4),INTENT(IN) :: I2
      REAL(8)   ,INTENT(IN) :: A(NR)
      REAL(8)   ,INTENT(IN) :: B(NR)
      REAL(8)   ,INTENT(IN) :: C(NR,NF,NF)
      REAL(8)   ,INTENT(IN) :: D(NR,NF)
      REAL(8)   ,INTENT(INOUT):: F(NR,NF)
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$DGLGEN(GIDS,NR,NF,I1,I2,A,B,C,D,F)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$DGLGEN(GIDS,NR,NF,I1,I2,A,B,C,D,F)
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$DGLGEN')
      END IF  
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$DGLGENC(GID,NR,NF,I1,I2,A,B,C,D,F)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NF
      INTEGER(4),INTENT(IN) :: I1
      INTEGER(4),INTENT(IN) :: I2
      REAL(8)   ,INTENT(IN) :: A(NR)
      REAL(8)   ,INTENT(IN) :: B(NR)
      COMPLEX(8),INTENT(IN) :: C(NR,NF,NF)
      COMPLEX(8),INTENT(IN) :: D(NR,NF)
      COMPLEX(8),INTENT(INOUT):: F(NR,NF)
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$DGLGENC(GIDS,NR,NF,I1,I2,A,B,C,D,F)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$DGLGENC(GIDS,NR,NF,I1,I2,A,B,C,D,F)
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$DGLGENC')
      END IF  
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$VERLETD1(GID,NR,F,DFDR)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F(NR)
      REAL(8)   ,INTENT(OUT):: DFDR(NR)
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$VERLETD1(GIDS,NR,F,DFDR)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$VERLETD1(GIDS,NR,F,DFDR)
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$VERLETD1')
      END IF  
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$VERLETD2(GID,NR,F,D2FDR2)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F(NR)
      REAL(8)   ,INTENT(OUT):: D2FDR2(NR)
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL LOGRADIAL$VERLETD2(GIDS,NR,F,D2FDR2)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$VERLETD2(GIDS,NR,F,D2FDR2)
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$VERLETD2')
      END IF  
      RETURN
      END      
!
!     .....................................................POISON.......
      SUBROUTINE RADIAL$POISSON(GID,NR,L,RHO,V)
!     ******************************************************************
!     **                                                              **
!     ** SOLVES THE RADIAL POISSON EQUALTION FOR A GIVEN              **
!     ** ANGULAR MOMENTUM COMPONENT OF THE CHARGE DENSITY             **
!     **                                                              **
!     ** INPUT :                                                      **
!     **   L            MAIN ANGULAR MOMENTUM QUANTUM NUMBER          **
!     **   RHO          INPUT CHARGE DENSITY                          **
!     ** OUTPUT :                                                      *
!     **   V            ELECTROSTATIC POTENTIAL                       **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID     ! GRID ID
      INTEGER(4),INTENT(IN) :: NR      ! NUMBER OF RADIAL GRID POINTS
      INTEGER(4),INTENT(IN) :: L       ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: RHO(NR) ! CHARGE DENSITY
      REAL(8)   ,INTENT(OUT):: V(NR)   ! ELECTROSTATIC POTENTIAL
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: AUX2(NR)
      REAL(8)               :: R(NR)
      REAL(8)               :: PI
      REAL(8)               :: FAC
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
      AUX1(:)=RHO(:)*R(:)**(L+2) 
      CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
      FAC=4.D0*PI/real(2*L+1,kind=8)
!     ==  R(1)=0 THERE IS A PRODUCT OF ZERO AND INFINITY ===============
      IF(R(1).EQ.0.D0) THEN
        V(1)=0.D0
        V(2:)=FAC*R(2:)**(-L-1)*AUX2(2:)
        AUX1(1)=0.D0  ! ASSUMING RHO \SIM R^L
        AUX1(2:)=RHO(2:)*R(2:)**(-L+1)
      ELSE
        V(:)=FAC*R(:)**(-L-1)*AUX2(:)
        AUX1(:)=RHO(:)*R(:)**(-L+1)
      END IF
      CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
      V(:)=V(:)+FAC*R(:)**L*(AUX2(NR)-AUX2(:))
      RETURN
      END
!
!     ......................................................................
      SUBROUTINE radial$NUCPOT(GID,NR,Z,POT)
!     **                                                                  **
!     **  ELECTROSTATIC POTENTIAL OF A NUCLEUS WITH FINITE RADIUS         **
!     **  THE NUCLEUS IS CONSIDERED AS A HOMOGENEOUSLY CHARGED SPHERE.    **
!     **  THE RADIUS IS RELATED  TO THE TOTAL MASS (NUMBER OF NUCLEONS),  **
!     **  WHICH CAN BE LOOKED UP KNOWING THE ATOMIC NUMBER Z.             **
!     **                                                                  **
!     **  the potential is the radial part only and needs to be           **
!     **  with the spherical harmonic Y0                                  **
!     **                                                                  **
!     **                                                                  **
      USE PERIODICTABLE_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID      ! GRID ID
      INTEGER(4),INTENT(IN) :: NR       ! #(RADIAL GRID POINTS)
      REAL(8)   ,INTENT(IN) :: Z        ! ATOMIC NUMBER
      REAL(8)   ,INTENT(OUT):: POT(NR)  ! POTENTIAL
      REAL(8)               :: RNUC     ! NUCLEAR RADIUS
      REAL(8)               :: R(NR)    ! RADIAL GRID
      REAL(8)               :: PI,Y0
      INTEGER(4)            :: IR
!     ***********************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL PERIODICTABLE$GET(z,'RNUC',RNUC)
      CALL RADIAL$R(GID,NR,R)
      DO IR=1,NR
        IF(R(IR).GT.RNUC) THEN
           POT(IR)=-Z/R(IR)
        ELSE
          POT(IR)=-Z/RNUC*(1.5D0-0.5D0*(R(IR)/RNUC)**2)
        END IF
      ENDDO
      POT(:)=POT(:)/Y0
      RETURN
      END
!
!     ...........................................MLTPOL.................
      SUBROUTINE RADIAL$MOMENT(GID,NR,L,RHO,QLM)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE MULTIPOLE MOMENT OF A (CHARGE) DISTRIBUTION  **
!     **  RHO(|R|)*Y_L(R).                                            **
!     **                                                              **
!     **    QLM = INT{D^3R: RHO(R) * |R|**L *Y_L(R) }                 **
!     **                                                              **
!     ****************************************** P.E. BLOECHL, 1997*****
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID      ! FIRST POINT ON RADIAL GRID
      INTEGER(4),INTENT(IN) :: NR      ! NUMBER OF RADIAL GRID POINTS
      INTEGER(4),INTENT(IN) :: L       ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: RHO(NR) ! CHARGE DENSITY
      REAL(8)   ,INTENT(OUT):: QLM     ! MULTIPOLE MOMENT
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: R(NR)
!     ******************************************************************
      CALL RADIAL$R(GID,NR,R)
      AUX1(:)=R(:)**(L+2)*RHO(:)
      CALL RADIAL$INTEGRAL(GID,NR,AUX1,QLM)
      RETURN
      END
!
!..........................................................................
MODULE SHLOGRADIAL_MODULE
TYPE SHLOGGRID_TYPE
  REAL(8)    :: R1
  REAL(8)    :: DEX
  INTEGER(4) :: NR
  REAL(8)   ,POINTER :: R(:)
END TYPE SHLOGGRID_TYPE
TYPE(SHLOGGRID_TYPE), ALLOCATABLE :: GRIDARRAY(:)
INTEGER(4)                        :: NGIDX=0
INTEGER(4)                        :: NGID=0
REAL(8)                           :: R1
REAL(8)                           :: DEX
INTEGER(4)                        :: NR
REAL(8)                           :: XEXP
REAL(8)   ,PARAMETER              :: RMIN=1.D-8  ! OFFSET OF INNERMOST POINT 
                                                 !TO AVOID DIVIDE BY ZERO
END MODULE SHLOGRADIAL_MODULE
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL_NEW(GID)
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT) :: GID
      INTEGER(4),PARAMETER ::NGIDBLOCKSIZE=10
      TYPE(SHLOGGRID_TYPE), ALLOCATABLE :: NEWARRAY(:)
      INTEGER(4)           :: NEWNGIDX
!     ********************************************************************
!       
!     ==================================================================
!     == REALLOCATE GRIDARRAY, IF IT IS TO SMALL                      ==
!     ==================================================================
      IF(NGID.GE.NGIDX) THEN
!       == SAVE CONTENTS OF GRIDARRAY AND DEALLOCATE IT
        IF(ALLOCATED(GRIDARRAY)) THEN
          ALLOCATE(NEWARRAY(NGIDX))
          NEWARRAY(:)=GRIDARRAY(:)
          DEALLOCATE(GRIDARRAY)
        END IF
!       == ALLOCATE GRIDARRAY WITH LARGER SIZE ========================
        NEWNGIDX=NGIDX+NGIDBLOCKSIZE
        ALLOCATE(GRIDARRAY(NEWNGIDX))
        GRIDARRAY(:)%R1=0.D0
        GRIDARRAY(:)%DEX=0.D0
        GRIDARRAY(:)%NR=0
!       == COPY CONTENTS SAVED FROM GRIDARRAY BACK ====================
        IF(ALLOCATED(NEWARRAY)) THEN
          GRIDARRAY(:NGIDX)=NEWARRAY(:)
          DEALLOCATE(NEWARRAY)
        END IF
        NGIDX=NEWNGIDX
      END IF
!
      NGID=NGID+1
      GID=NGID
      GRIDARRAY(NGID)%R1=0.D0
      GRIDARRAY(NGID)%DEX=0.D0
      GRIDARRAY(NGID)%NR=0
      NULLIFY(GRIDARRAY(NGID)%R)
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL_RESOLVE(GID)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      REAL(8)               :: RI
      INTEGER(4)            :: IR
!     *******************************************************************
      IF(GID.GT.NGID) THEN
        CALL ERROR$MSG('GRID ID OUT OF RANGE')
        CALL ERROR$I4VAL('GID',GID)
        CALL ERROR$I4VAL('NGID',NGID)
        CALL ERROR$STOP('SHLOGRADIAL_RESOLVE')
      END IF
      R1=GRIDARRAY(GID)%R1
      DEX=GRIDARRAY(GID)%DEX
      NR=GRIDARRAY(GID)%NR
      XEXP=EXP(DEX)
      IF(.NOT.ASSOCIATED(GRIDARRAY(GID)%R)) THEN
        ALLOCATE(GRIDARRAY(GID)%R(NR))
        RI=R1/XEXP
        DO IR=1,NR 
          RI=RI*XEXP
          GRIDARRAY(GID)%R(IR)=RI-R1
        ENDDO
        GRIDARRAY(GID)%R(1)=RMIN
      END IF
      RETURN
      END
!
!      .................................................................
       SUBROUTINE SHLOGRADIAL$SETR8(GID,ID,VAL)
!      *****************************************************************
!      *****************************************************************
       USE SHLOGRADIAL_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: GID
       CHARACTER(*),INTENT(IN) :: ID
       REAL(8)     ,INTENT(IN) :: VAL
!      ***************************************************************
       IF(GID.GT.NGID) THEN
         CALL ERROR$MSG('GRID ID OUT OF RANGE')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NGID',NGID)
         CALL ERROR$STOP('SHLOGRADIAL$SETR8')
       END IF
       IF(ID.EQ.'R1') THEN
         GRIDARRAY(GID)%R1=VAL
       ELSE IF(ID.EQ.'DEX') THEN
         GRIDARRAY(GID)%DEX=VAL
       ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SHLOGRADIAL$SETR8')
       END IF
       RETURN
     END SUBROUTINE SHLOGRADIAL$SETR8
!
!      .................................................................
       SUBROUTINE SHLOGRADIAL$GETR8(GID,ID,VAL)
!      *****************************************************************
!      *****************************************************************
       USE SHLOGRADIAL_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: GID
       CHARACTER(*),INTENT(IN) :: ID
       REAL(8)     ,INTENT(OUT):: VAL
!      ***************************************************************
       IF(GID.GT.NGID) THEN
         CALL ERROR$MSG('GRID ID OUT OF RANGE')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NGID',NGID)
         CALL ERROR$STOP('SHLOGRADIAL$GETR8')
       END IF
       IF(ID.EQ.'R1') THEN
         VAL=GRIDARRAY(GID)%R1
       ELSE IF(ID.EQ.'DEX') THEN
         VAL=GRIDARRAY(GID)%DEX
       ELSE
          CALL ERROR$MSG('ID NOT RECOGNIZED')
          CALL ERROR$CHVAL('ID',ID)
          CALL ERROR$STOP('SHLOGRADIAL$GETR8')
       END IF
       RETURN
       END SUBROUTINE SHLOGRADIAL$GETR8
!
!      .................................................................
       SUBROUTINE SHLOGRADIAL$SETI4(GID,ID,VAL)
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
       USE SHLOGRADIAL_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: GID
       CHARACTER(*),INTENT(IN) :: ID
       INTEGER(4)  ,INTENT(IN) :: VAL
!      ***************************************************************
       IF(GID.GT.NGID) THEN
         CALL ERROR$MSG('GRID ID OUT OF RANGE')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NGID',NGID)
         CALL ERROR$STOP('SHLOGRADIAL$SETI4')
       END IF
       IF(ID.EQ.'NR') THEN
         GRIDARRAY(GID)%NR=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('SHLOGRADIAL$SETI4')
       END IF
       RETURN
       END SUBROUTINE SHLOGRADIAL$SETI4
!
!      .................................................................
       SUBROUTINE SHLOGRADIAL$GETI4(GID,ID,VAL)
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
       USE SHLOGRADIAL_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: GID
       CHARACTER(*),INTENT(IN) :: ID
       INTEGER(4)  ,INTENT(OUT):: VAL
!      ***************************************************************
       IF(GID.GT.NGID) THEN
         CALL ERROR$MSG('GRID ID OUT OF RANGE')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NGID',NGID)
         CALL ERROR$STOP('SHLOGRADIAL$GETI4')
       END IF
       IF(ID.EQ.'NR') THEN
         VAL=GRIDARRAY(GID)%NR
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('SHLOGRADIAL$GETI4')
       END IF
       RETURN
       END SUBROUTINE SHLOGRADIAL$GETI4
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$R(GID,NR_,R)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(OUT):: R(NR_)
      REAL(8)               :: RI     
      INTEGER(4)            :: IR
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('SHLOGRADIAL$R')
      END IF
      IF(ASSOCIATED(GRIDARRAY(GID)%R)) THEN
        R(:)=GRIDARRAY(GID)%R(:)
        RETURN
      END IF
      RI=R1/XEXP
      DO IR=1,NR 
        RI=RI*XEXP
        R(IR)=RI-R1
      ENDDO
      R(1)=RMIN
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$XOFR(GID,R0,X0)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: X0
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      X0=1.D0+LOG(R0/R1+1.D0)/DEX
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$VALUE(GID,NR_,F,R0,F0)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      REAL(8)               :: RARR(4),FARR(4),RI,XI
      INTEGER(4)            :: IR,IR1
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('SHLOGRADIAL$VALUE')
      END IF
      XI=1.D0+LOG(R0/R1+1.D0)/DEX
      IR1=INT(XI)-1
      IR1=MAX(1,IR1)
      IR1=MIN(IR1,NR-3)
      FARR(:)=F(IR1:IR1+3)
      RI=R1*EXP(DEX*REAL(IR1-1,KIND=8))
      DO IR=1,4 
        RARR(IR)=RI-R1
        RI=RI*XEXP
      ENDDO
      IF(IR1.EQ.1) RARR(1)=RMIN
      CALL RADIAL_POLYNOMIALVALUE(4,RARR,FARR,R0,F0)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$INTEGRATE(GID,NR_,F,G)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(OUT):: G(NR_)
      REAL(8)               :: RI     
      REAL(8)               :: FBAR(NR_)
      INTEGER(4)            :: IR
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('SHLOGRADIAL$INTEGRATE')
      END IF
      RI=R1/XEXP
      DO IR=1,NR 
        RI=RI*XEXP
        FBAR(IR)=F(IR)*(DEX*RI)
      ENDDO
      CALL RADIAL_INTEGRATEEQUISPACED(NR,FBAR,G)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$DERIVE(GID,NR_,F,G)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(OUT):: G(NR_)
      REAL(8)               :: RI     
      INTEGER(4)            :: IR
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('SHLOGRADIAL$DERIVE')
      END IF
      CALL RADIAL_DERIVEEQUISPACED(NR_,F,G)
      RI=R1/XEXP
      DO IR=1,NR 
        RI=RI*XEXP
        G(IR)=G(IR)/(DEX*RI)
      ENDDO
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$DERIVATIVE(GID,NR_,F,R0,F0)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      REAL(8)               :: RARR(4),FARR(4),RI,XI
      INTEGER(4)            :: IR,IR1
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('SHLOGRADIAL$DERIVATIVE')
      END IF
      XI=1.D0+LOG(R0/R1+1.D0)/DEX
      IR1=INT(XI)-1
      IR1=MAX(1,IR1)
      IR1=MIN(IR1,NR-3)
      FARR(:)=F(IR1:IR1+3)
      RI=R1*EXP(DEX*REAL(IR1-1,KIND=8))
      DO IR=1,4 
        RARR(IR)=RI-R1
        RI=RI*XEXP
      ENDDO
      CALL RADIAL_POLYNOMIALDERIVATIVE(4,RARR,FARR,R0,F0)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$DGL(GID,IDIR,NR_,A,B,C,D,F)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: IDIR
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: A(NR_)
      REAL(8)   ,INTENT(IN) :: B(NR_)
      REAL(8)   ,INTENT(IN) :: C(NR_)
      REAL(8)   ,INTENT(IN) :: D(NR_)
      REAL(8)   ,INTENT(INOUT):: F(NR_)
      REAL(8)               :: R(NR_)
      REAL(8)               :: DRDX(NR_)
      REAL(8)               :: B1(NR_)
      REAL(8)               :: C1(NR_)
      REAL(8)               :: D1(NR_)
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('SHLOGRADIAL$DGL')
      END IF
      CALL SHLOGRADIAL$R(GID,NR,R)
      DRDX(:)=DEX*(R(:)+R1)   
      B1(:)=B(:)*DRDX(:)-DEX*A(:)
      C1(:)=C(:)*DRDX(:)**2
      D1(:)=D(:)*DRDX(:)**2
      CALL RADIAL_DGLEQUISPACED(IDIR,NR,A,B1,C1,D1,F)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$DGLGEN(GID,NR_,NF,I1,I2,A,B,C,D,F)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      INTEGER(4),INTENT(IN) :: NF
      INTEGER(4),INTENT(IN) :: I1
      INTEGER(4),INTENT(IN) :: I2
      REAL(8)   ,INTENT(IN) :: A(NR_)
      REAL(8)   ,INTENT(IN) :: B(NR_)
      REAL(8)   ,INTENT(IN) :: C(NR_,NF,NF)
      REAL(8)   ,INTENT(IN) :: D(NR_,NF)
      REAL(8)   ,INTENT(INOUT):: F(NR_,NF)
      REAL(8)               :: R(NR_)
      REAL(8)               :: DRDX(NR_)
      REAL(8)               :: B1(NR_)
      REAL(8)               :: C1(NR_,NF,NF)
      REAL(8)               :: D1(NR_,NF)
      INTEGER(4)            :: I,J
      INTEGER(4)            :: IMIN,IMAX
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('SHLOGRADIAL$DGLGEN')
      END IF
      CALL SHLOGRADIAL$R(GID,NR,R)
      IMIN=MIN(I1,I2)
      IMAX=MAX(I1,I2)
      DRDX(IMIN:IMAX)=DEX*(R(IMIN:IMAX)+R1)   
      B1(IMIN:IMAX)=B(IMIN:IMAX)*DRDX(IMIN:IMAX)-DEX*A(IMIN:IMAX)
      DRDX(IMIN:IMAX)=DRDX(IMIN:IMAX)**2
      DO I=1,NF
        D1(IMIN:IMAX,I)=D(IMIN:IMAX,I)*DRDX(IMIN:IMAX)
        DO J=1,NF
          C1(IMIN:IMAX,I,J)=C(IMIN:IMAX,I,J)*DRDX(IMIN:IMAX)
        ENDDO
      ENDDO
      CALL RADIAL_DGLEQUISPACEDGEN(NR,NF,I1,I2,A,B1,C1,D1,F)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$DGLGENC(GID,NR_,NF,I1,I2,A,B,C,D,F)
!     **                                                               **
!     **  TRANSFORMS THE EQUATION                                      **
!     **    [A(R)\PARTIAL_R^2+B(R)\PARTIAL_R+C(R)]F(R)=D(R)            **
!     **  INTO                                                         **
!     **    [A(R)\PARTIAL^2_X+B1(R)\PARTIAL_X+C1(X)]F(X)=D1(X)         **
!     **  WITH R(X)=R1\EXP(\ALPHA(X-1))                                **
!     **                                                               **
!     **  AND CALLS RADIAL_DGLEQUISPACEDGENC TO SOLVE THE RESULTING    **
!     **  DGL ON THE EQUI-SPACED GRID                                  **
!     **                                                               **
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      INTEGER(4),INTENT(IN) :: NF
      INTEGER(4),INTENT(IN) :: I1
      INTEGER(4),INTENT(IN) :: I2
      REAL(8)   ,INTENT(IN) :: A(NR_)
      REAL(8)   ,INTENT(IN) :: B(NR_)
      COMPLEX(8),INTENT(IN) :: C(NR_,NF,NF)
      COMPLEX(8),INTENT(IN) :: D(NR_,NF)
      COMPLEX(8),INTENT(INOUT):: F(NR_,NF)
      REAL(8)               :: R(NR_)
      REAL(8)               :: DRDX(NR_)
      REAL(8)               :: B1(NR_)
      COMPLEX(8)            :: C1(NR_,NF,NF)
      COMPLEX(8)            :: D1(NR_,NF)
      INTEGER(4)            :: I,J
      INTEGER(4)            :: IMIN,IMAX
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('SHLOGRADIAL$DGLGENC')
      END IF
      CALL SHLOGRADIAL$R(GID,NR,R)
      IMIN=MIN(I1,I2)
      IMAX=MAX(I1,I2)
      DRDX(IMIN:IMAX)=DEX*(R(IMIN:IMAX)+R1)   
      B1(IMIN:IMAX)=B(IMIN:IMAX)*DRDX(IMIN:IMAX)-DEX*A(IMIN:IMAX)
      DRDX(IMIN:IMAX)=DRDX(IMIN:IMAX)**2
      DO I=1,NF
        D1(IMIN:IMAX,I)=D(IMIN:IMAX,I)*DRDX(IMIN:IMAX)
        DO J=1,NF
          C1(IMIN:IMAX,I,J)=C(IMIN:IMAX,I,J)*DRDX(IMIN:IMAX)
        ENDDO
      ENDDO
      CALL RADIAL_DGLEQUISPACEDGENC(NR,NF,I1,I2,A,B1,C1,D1,F)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$VERLETD1(GID,NR_,F,DFDR)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(OUT):: DFDR(NR_)
      REAL(8)               :: R(NR_)
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('SHLOGRADIAL$VERLETD1')
      END IF
      CALL RADIAL_VERLETD1EQUISPACED(NR,F,DFDR)
      CALL SHLOGRADIAL$R(GID,NR,R)
      DFDR(:)=DFDR(:)/(DEX*(R(:)+R1))
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE SHLOGRADIAL$VERLETD2(GID,NR_,F,D2FDR2)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(OUT):: D2FDR2(NR_)
      REAL(8)               :: R(NR_)
      REAL(8)               :: DRDX(NR_)
      REAL(8)               :: DFDR(NR_)
!     ******************************************************************
      CALL SHLOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('SHLOGRADIAL$VERLETD2')
      END IF
      CALL RADIAL_VERLETD2EQUISPACED(NR,F,D2FDR2)
      CALL RADIAL_VERLETD1EQUISPACED(NR,F,DFDR)
      CALL SHLOGRADIAL$R(GID,NR,R)
      DRDX(:)=DEX*(R(:)+R1)
      D2FDR2(:)=(D2FDR2(:)-DEX*DFDR(:))/DRDX(:)**2
      RETURN
      END      

!
!..........................................................................
MODULE LOGRADIAL_MODULE
TYPE LOGGRID_TYPE
  REAL(8)    :: R1
  REAL(8)    :: DEX
  INTEGER(4) :: NR
END TYPE LOGGRID_TYPE
TYPE(LOGGRID_TYPE), ALLOCATABLE :: GRIDARRAY(:)
INTEGER(4)                        :: NGIDX=0
INTEGER(4)                        :: NGID=0
REAL(8)   :: R1
REAL(8)   :: DEX
INTEGER(4):: NR
REAL(8)   :: XEXP
END MODULE LOGRADIAL_MODULE
!
!     ...................................................................
      SUBROUTINE LOGRADIAL_NEW(GID)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT):: GID
      INTEGER(4),PARAMETER  ::NGIDBLOCKSIZE=10
      TYPE(LOGGRID_TYPE), ALLOCATABLE :: NEWARRAY(:)
      INTEGER(4)           :: NEWNGIDX
!     ********************************************************************
!       
!     ==================================================================
!     == REALLOCATE GRIDARRAY, IF IT IS TO SMALL                      ==
!     ==================================================================
      IF(NGID.GE.NGIDX) THEN
!       == SAVE CONTENTS OF GRIDARRAY AND DEALLOCATE IT
        IF(ALLOCATED(GRIDARRAY)) THEN
          ALLOCATE(NEWARRAY(NGIDX))
          NEWARRAY(:)=GRIDARRAY(:)
          DEALLOCATE(GRIDARRAY)
        END IF
!       == ALLOCATE GRIDARRAY WITH LARGER SIZE ========================
        NEWNGIDX=NGIDX+NGIDBLOCKSIZE
        ALLOCATE(GRIDARRAY(NEWNGIDX))
        GRIDARRAY(:)%R1=0.D0
        GRIDARRAY(:)%DEX=0.D0
        GRIDARRAY(:)%NR=0
!       == COPY CONTENTS SAVED FROM GRIDARRAY BACK ====================
        IF(ALLOCATED(NEWARRAY)) THEN
          GRIDARRAY(:NGIDX)=NEWARRAY(:)
          DEALLOCATE(NEWARRAY)
        END IF
        NGIDX=NEWNGIDX
      END IF
!
      NGID=NGID+1
      GID=NGID
      RETURN
      END
!
!     ...................................................................
      SUBROUTINE LOGRADIAL_RESOLVE(GID)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
!     *******************************************************************
      IF(GID.GT.NGID) THEN
        CALL ERROR$MSG('GRID ID OUT OF RANGE')
        CALL ERROR$I4VAL('GID',GID)
        CALL ERROR$I4VAL('NGID',NGID)
        CALL ERROR$STOP('LOGRADIAL_RESOLVE')
      END IF
      R1=GRIDARRAY(GID)%R1
      DEX=GRIDARRAY(GID)%DEX
      NR=GRIDARRAY(GID)%NR
      XEXP=EXP(DEX)
      RETURN
      END
!
!      .................................................................
       SUBROUTINE LOGRADIAL$GETR8(GID,ID,VAL)
!      *****************************************************************
!      *****************************************************************
       USE LOGRADIAL_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: GID
       CHARACTER(*),INTENT(IN) :: ID
       REAL(8)     ,INTENT(OUT):: VAL
!      ***************************************************************
       IF(GID.GT.NGID) THEN
         CALL ERROR$MSG('GRID ID OUT OF RANGE')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NGID',NGID)
         CALL ERROR$STOP('LOGRADIAL$GETR8')
       END IF
       IF(ID.EQ.'R1') THEN
         VAL=GRIDARRAY(GID)%R1
       ELSE IF(ID.EQ.'DEX') THEN
         VAL=GRIDARRAY(GID)%DEX
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('LOGRADIAL$GETR8')
       END IF
       RETURN
       END SUBROUTINE LOGRADIAL$GETR8
!
!      .................................................................
       SUBROUTINE LOGRADIAL$SETR8(GID,ID,VAL)
!      *****************************************************************
!      *****************************************************************
       USE LOGRADIAL_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: GID
       CHARACTER(*),INTENT(IN) :: ID
       REAL(8)     ,INTENT(IN) :: VAL
!      ***************************************************************
       IF(GID.GT.NGID) THEN
         CALL ERROR$MSG('GRID ID OUT OF RANGE')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NGID',NGID)
         CALL ERROR$STOP('LOGRADIAL$SETR8')
       END IF
       IF(ID.EQ.'R1') THEN
         GRIDARRAY(GID)%R1=VAL
       ELSE IF(ID.EQ.'DEX') THEN
         GRIDARRAY(GID)%DEX=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('LOGRADIAL$SETR8')
       END IF
       RETURN
     END SUBROUTINE LOGRADIAL$SETR8
!
!      .................................................................
       SUBROUTINE LOGRADIAL$SETI4(GID,ID,VAL)
!      *****************************************************************
!      *****************************************************************
       USE LOGRADIAL_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: GID
       CHARACTER(*),INTENT(IN) :: ID
       INTEGER(4)  ,INTENT(IN) :: VAL
!      ***************************************************************
       IF(GID.GT.NGID) THEN
         CALL ERROR$MSG('GRID ID OUT OF RANGE')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NGID',NGID)
         CALL ERROR$STOP('LOGRADIAL$SETI4')
       END IF
       IF(ID.EQ.'NR') THEN
         GRIDARRAY(GID)%NR=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('LOGRADIAL$SETI4')
       END IF
       RETURN
       END SUBROUTINE LOGRADIAL$SETI4
!
!      .................................................................
       SUBROUTINE LOGRADIAL$GETI4(GID,ID,VAL)
!      *****************************************************************
!      *****************************************************************
       USE LOGRADIAL_MODULE
       IMPLICIT NONE
       INTEGER(4)  ,INTENT(IN) :: GID
       CHARACTER(*),INTENT(IN) :: ID
       INTEGER(4)  ,INTENT(OUT):: VAL
!      ***************************************************************
       IF(GID.GT.NGID) THEN
         CALL ERROR$MSG('GRID ID OUT OF RANGE')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NGID',NGID)
         CALL ERROR$STOP('LOGRADIAL$GETI4')
       END IF
       IF(ID.EQ.'NR') THEN
         VAL=GRIDARRAY(GID)%NR
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('LOGRADIAL$GETI4')
       END IF
       RETURN
       END SUBROUTINE LOGRADIAL$GETI4
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$R(GID,NR_,R)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(OUT):: R(NR_)
      REAL(8)               :: RI     
      INTEGER(4)            :: IR
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
        CALL ERROR$MSG('GRID SIZE INCONSISTENT')
        CALL ERROR$I4VAL('NR',NR)
        CALL ERROR$I4VAL('NR_',NR_)
        CALL ERROR$STOP('LOGRADIAL$R')
      END IF
      RI=R1/XEXP
      DO IR=1,NR 
        RI=RI*XEXP
        R(IR)=RI
      ENDDO
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$XOFR(GID,R0,X0)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: X0
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      X0=1.D0+LOG(R0/R1)/DEX
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$VALUE(GID,NR_,F,R0,F0)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      REAL(8)               :: RARR(4),FARR(4),RI,XI
      INTEGER(4)            :: IR,IR1
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
        CALL ERROR$MSG('GRID SIZE INCONSISTENT')
        CALL ERROR$I4VAL('NR',NR)
        CALL ERROR$I4VAL('NR_',NR_)
        CALL ERROR$STOP('LOGRADIAL$VALUE')
      END IF
      IF(R0.GT.R1) THEN
        XI=1.D0+LOG(R0/R1)/DEX
        IR1=INT(XI)-1
        IR1=MAX(1,IR1)
        IR1=MIN(IR1,NR-3)
      ELSE
        IR1=1
      END IF
      FARR(:)=F(IR1:IR1+3)
      RI=R1*EXP(DEX*(IR1-1))
      DO IR=1,4
        RARR(IR)=RI
        RI=RI*XEXP
      ENDDO
      CALL RADIAL_POLYNOMIALVALUE(4,RARR,FARR,R0,F0)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$INTEGRATE(GID,NR_,F,G)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(OUT):: G(NR_)
      REAL(8)               :: RI     
      REAL(8)               :: FBAR(NR_)
      REAL(8)               :: CN(4),RARR(4),SVAR
      INTEGER(4)            :: IR
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
        CALL ERROR$MSG('GRID SIZE INCONSISTENT')
        CALL ERROR$I4VAL('NR',NR)
        CALL ERROR$I4VAL('NR_',NR_)
        CALL ERROR$STOP('LOGRADIAL$INTEGRATE')
      END IF
!
!     ==================================================================
!     ==  INTEGRATE FROM R1 TO RI                                     ==
!     ==================================================================
      RI=R1/XEXP
      DO IR=1,NR 
        RI=RI*XEXP
        FBAR(IR)=F(IR)*(DEX*RI)
      ENDDO
      CALL RADIAL_INTEGRATEEQUISPACED(NR,FBAR,G)
!
!     ==================================================================
!     ==  INTEGRATE FROM 0 TO R1                                      ==
!     ==================================================================
      RI=R1/XEXP
      DO IR=1,4
        RI=RI*XEXP
        RARR(IR)=RI
      ENDDO
      CALL RADIAL_POLYNOMIALCOEFFICIENTS(4,RARR,F(1:4),0.D0,CN)
      SVAR=R1*(CN(1)+R1*(CN(2)*0.5D0+R1*(CN(3)/3.D0+CN(4)*0.25D0*R1)))
      G(:)=G(:)+SVAR
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$DERIVE(GID,NR_,F,G)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(OUT):: G(NR_)
      REAL(8)               :: RI     
      INTEGER(4)            :: IR
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
        CALL ERROR$MSG('SIZE INCONSISTENT')
        CALL ERROR$I4VAL('NR',NR)
        CALL ERROR$I4VAL('NR_',NR_)
        CALL ERROR$STOP('LOGRADIAL$DERIVE')
      END IF
      CALL RADIAL_DERIVEEQUISPACED(NR,F,G)
      RI=R1/XEXP
      DO IR=1,NR 
        RI=RI*XEXP
        G(IR)=G(IR)/(DEX*RI)
      ENDDO
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$DERIVATIVE(GID,NR_,F,R0,F0)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      REAL(8)               :: RARR(4),FARR(4),RI,XI
      INTEGER(4)            :: IR,IR1
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
        PRINT*,' ERROR',NR_,NR
        STOP 'ERROR LOGRADIAL$DERIVATIVE'
      END IF
      IF(R0.GT.R1) THEN
        XI=1.D0+LOG(R0/R1)/DEX
        IR1=INT(XI)-1
        IR1=MAX(1,IR1)
        IR1=MIN(IR1,NR-3)
      ELSE
        IR1=1
      END IF
      FARR(:)=F(IR1:IR1+3)
      RI=R1*EXP(DEX*(IR1-1))
      DO IR=1,4
        RARR(IR)=RI
        RI=RI*XEXP
      ENDDO
      CALL RADIAL_POLYNOMIALDERIVATIVE(4,RARR,FARR,R0,F0)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$DGL(GID,IDIR,NR_,A,B,C,D,F)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: IDIR
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: A(NR_)
      REAL(8)   ,INTENT(IN) :: B(NR_)
      REAL(8)   ,INTENT(IN) :: C(NR_)
      REAL(8)   ,INTENT(IN) :: D(NR_)
      REAL(8)   ,INTENT(INOUT):: F(NR_)
      REAL(8)               :: R(NR_)
      REAL(8)               :: DRDX(NR_)
      REAL(8)               :: B1(NR_)
      REAL(8)               :: C1(NR_)
      REAL(8)               :: D1(NR_)
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
        CALL ERROR$MSG('SIZE INCONSISTENT')
        CALL ERROR$I4VAL('NR',NR)
        CALL ERROR$I4VAL('NR_',NR_)
        CALL ERROR$STOP('LOGRADIAL$DGL')
      END IF
      CALL LOGRADIAL$R(GID,NR,R)
      DRDX(:)=DEX*R(:)  ! DR/DX
      B1(:)=B(:)*DRDX(:)-DEX*A(:)
      C1(:)=C(:)*DRDX(:)**2
      D1(:)=D(:)*DRDX(:)**2
      CALL RADIAL_DGLEQUISPACED(IDIR,NR,A,B1,C1,D1,F)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$DGLGEN(GID,NR_,NF,I1,I2,A,B,C,D,F)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      INTEGER(4),INTENT(IN) :: NF
      INTEGER(4),INTENT(IN) :: I1
      INTEGER(4),INTENT(IN) :: I2
      REAL(8)   ,INTENT(IN) :: A(NR_)
      REAL(8)   ,INTENT(IN) :: B(NR_)
      REAL(8)   ,INTENT(IN) :: C(NR_,NF,NF)
      REAL(8)   ,INTENT(IN) :: D(NR_,NF)
      REAL(8)   ,INTENT(INOUT):: F(NR_,NF)
      REAL(8)               :: R(NR_)
      REAL(8)               :: DRDX(NR_)
      REAL(8)               :: B1(NR_)
      REAL(8)               :: C1(NR_,NF,NF)
      REAL(8)               :: D1(NR_,NF)
      INTEGER(4)            :: I,J
      INTEGER(4)            :: IMIN,IMAX
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('LOGRADIAL$DGLGEN')
      END IF
      CALL LOGRADIAL$R(GID,NR,R)
      IMIN=MIN(I1,I2)
      IMAX=MAX(I1,I2)
      DRDX(IMIN:IMAX)=DEX*R(IMIN:IMAX)
      B1(IMIN:IMAX)=B(IMIN:IMAX)*DRDX(IMIN:IMAX)-DEX*A(IMIN:IMAX)
      DRDX(IMIN:IMAX)=DRDX(IMIN:IMAX)**2
      DO I=1,NF
        D1(IMIN:IMAX,I)=D(IMIN:IMAX,I)*DRDX(IMIN:IMAX)
        DO J=1,NF
          C1(IMIN:IMAX,I,J)=C(IMIN:IMAX,I,J)*DRDX(IMIN:IMAX)
        ENDDO
      ENDDO
      CALL RADIAL_DGLEQUISPACEDGEN(NR,NF,I1,I2,A,B1,C1,D1,F)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$DGLGENC(GID,NR_,NF,I1,I2,A,B,C,D,F)
!     **                                                               **
!     **  TRANSFORMS THE EQUATION                                      **
!     **    [A(R)\PARTIAL_R^2+B(R)\PARTIAL_R+C(R)]F(R)=D(R)            **
!     **  INTO                                                         **
!     **    [A(R)\PARTIAL^2_X+B1(R)\PARTIAL_X+C1(X)]F(X)=D1(X)         **
!     **  WITH R(X)=R1\EXP(\ALPHA(X-1))                                **
!     **                                                               **
!     **  AND CALLS RADIAL_DGLEQUISPACEDGENC TO SOLVE THE RESULTING    **
!     **  DGL ON THE EQUI-SPACED GRID                                  **
!     **                                                               **
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      INTEGER(4),INTENT(IN) :: NF
      INTEGER(4),INTENT(IN) :: I1
      INTEGER(4),INTENT(IN) :: I2
      REAL(8)   ,INTENT(IN) :: A(NR_)
      REAL(8)   ,INTENT(IN) :: B(NR_)
      COMPLEX(8),INTENT(IN) :: C(NR_,NF,NF)
      COMPLEX(8),INTENT(IN) :: D(NR_,NF)
      COMPLEX(8),INTENT(INOUT):: F(NR_,NF)
      REAL(8)               :: R(NR_)
      REAL(8)               :: DRDX(NR_)
      REAL(8)               :: B1(NR_)
      COMPLEX(8)            :: C1(NR_,NF,NF)
      COMPLEX(8)            :: D1(NR_,NF)
      INTEGER(4)            :: I,J
      INTEGER(4)            :: IMIN,IMAX
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
         CALL ERROR$MSG('INCONSISTENT NUMBER OF GRID POINTS')
         CALL ERROR$I4VAL('GID',GID)
         CALL ERROR$I4VAL('NR',NR)
         CALL ERROR$I4VAL('NR_',NR_)
         CALL ERROR$STOP('LOGRADIAL$DGLGENC')
      END IF
      CALL LOGRADIAL$R(GID,NR,R)
      IMIN=MIN(I1,I2)
      IMAX=MAX(I1,I2)
      DRDX(IMIN:IMAX)=DEX*R(IMIN:IMAX)
      B1(IMIN:IMAX)=B(IMIN:IMAX)*DRDX(IMIN:IMAX)-DEX*A(IMIN:IMAX)
      DRDX(IMIN:IMAX)=DRDX(IMIN:IMAX)**2
      DO I=1,NF
        D1(IMIN:IMAX,I)=D(IMIN:IMAX,I)*DRDX(IMIN:IMAX)
        DO J=1,NF
          C1(IMIN:IMAX,I,J)=C(IMIN:IMAX,I,J)*DRDX(IMIN:IMAX)
        ENDDO
      ENDDO
      CALL RADIAL_DGLEQUISPACEDGENC(NR,NF,I1,I2,A,B1,C1,D1,F)
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$VERLETD1(GID,NR_,F,DFDR)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(OUT):: DFDR(NR_)
      REAL(8)               :: R(NR_)
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
        CALL ERROR$MSG('SIZE INCONSISTENT')
        CALL ERROR$I4VAL('NR',NR)
        CALL ERROR$I4VAL('NR_',NR_)
        CALL ERROR$STOP('LOGRADIAL$VERLETD1')
      END IF
      CALL RADIAL_VERLETD1EQUISPACED(NR,F,DFDR)
      CALL LOGRADIAL$R(GID,NR,R)
      DFDR(:)=DFDR(:)/(DEX*R(:))
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE LOGRADIAL$VERLETD2(GID,NR_,F,D2FDR2)
      USE LOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      REAL(8)   ,INTENT(OUT):: D2FDR2(NR_)
      REAL(8)               :: R(NR_)
      REAL(8)               :: DRDX(NR_)
      REAL(8)               :: DFDR(NR_)
!     ******************************************************************
      CALL LOGRADIAL_RESOLVE(GID)
      IF(NR_.NE.NR) THEN
        CALL ERROR$MSG('SIZE INCONSISTENT')
        CALL ERROR$I4VAL('NR',NR)
        CALL ERROR$I4VAL('NR_',NR_)
        CALL ERROR$STOP('LOGRADIAL$VERLETD1')
      END IF
      CALL RADIAL_VERLETD2EQUISPACED(NR,F,D2FDR2)
      CALL RADIAL_VERLETD1EQUISPACED(NR,F,DFDR)
      CALL LOGRADIAL$R(GID,NR,R)
      DRDX(:)=DEX*R(:)
      D2FDR2(:)=(D2FDR2(:)-DEX*DFDR(:))/DRDX(:)**2
      RETURN
      END      

!
!.......................................................................
MODULE BESSELTRANSFORM_MODULE
!***********************************************************************
!**                                                                   **
!**  NAME: BESSELTRANSFORM                                            **
!**                                                                   **
!**  PURPOSE: PERFORMS BESSELTRANSFORM BETWEEN LOGARITHMIC RADIAL     **
!**    GRIDS IN REAL AND RECOPROCAL SPACE:                            **
!**                                                                   **
!**    CALCULATES THE SPHERICAL BESSEL TRANSFORM OF ORDER L           **
!**       F_L(G)=INT_0^INFTY : X**2*F_L(R)*J_L(|G|*|R|)               **
!**    FOR A FUNCTION F GIVEN ON THE LOGARITHMIC GRID                 **
!**            R(I)=R1*EXP((I-1)*DEX)  ; I=1,.,NP                     **
!**    IN THE PRESENT IMPLEMENTATION NP MUST BE A POWER OF 2          **
!**    THE RESULT IS RETURNED AS FOFG ON THE LOGARITHMIC GRID         **
!**            G(I)=G1*EXP((I-1)*H)    ;I=1,NP                        **
!**                                                                   **
!**  FUNCTIONS:                                                       **
!**    BESSELTRANSFORM(L,NP,R1,G1,DEX,F,G,DISC)                       **
!**    BESSELTRANSFORM$CLEAR                                          **
!**                                                                   **
!**  REFERENCE:                                                       **
!**    J.D. TALMAN. COMP. PHYS. COMMUN. 30 (1983) 93                  **
!**                                                                   **
!*********************************************** P.E. BLOECHL, 1996*****
COMPLEX(8),ALLOCATABLE :: TA(:)     ! WORK ARRAY, SAVE
COMPLEX(8),ALLOCATABLE :: TRBE(:,:) ! WORK ARRAY,SAVE
COMPLEX(8),ALLOCATABLE :: WW(:)     ! WORK ARRAY,SAVE
INTEGER(4)             :: NPSAVE=0
REAL(8)                :: R1SAVE=0.D0
REAL(8)                :: G1SAVE=0.D0
REAL(8)                :: DEXSAVE=0.D0
INTEGER(4)             :: NSAVE=0
! BESSELTRANSFORM REQUIRES MAPPING FROM SHLOG TO LOG GRID TYPE
INTEGER(4),PARAMETER   :: NMAPX=20
INTEGER(4)             :: NMAP=0
INTEGER(4)             :: MAP(4,NMAPX)
CONTAINS
!     ..................................................................
      SUBROUTINE BTSMALLG(L,R1,G1,DEX,NP,F,G,N,WW,TA)
!     ******************************************************************
!     **                                                              **
!     **  BESSELTRANSFORM ACCURATE FOR SMALL G FOR L>1                **
!     **  THIS ROUTINE CONTAINED IN MODULE BESSELTRANSFORM_MODULE     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L
      INTEGER(4),INTENT(IN) :: NP
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: G1
      REAL(8)   ,INTENT(IN) :: DEX
      REAL(8)   ,INTENT(IN) :: F(NP)
      REAL(8)   ,INTENT(OUT):: G(NP)
      COMPLEX(8),INTENT(IN) :: WW(NP)
      COMPLEX(8),INTENT(IN) :: TA(NP)
      COMPLEX(8)            :: XA(NP)
      REAL(8)               :: PI
      INTEGER(4)            :: NH,NHP
      REAL(8)               :: DT,T,AA,BB,CM
      REAL(8)               :: RIX,XEXPX
      INTEGER(4)            :: JJ,I
!     ******************************************************************
      NH=NP/2
      NHP=NH+1
      PI=4.D0*ATAN(1.D0)
      RIX=R1**(L+1.5)
      XEXPX=EXP((L+1.5)*DEX) 
      DO I=1,NP 
        XA(I)=F(I)*RIX 
        RIX=RIX*XEXPX 
      ENDDO
      CALL NLOGN(N,XA,WW,NP) 
      DO I=1,NH 
        XA(I)=TA(I)*XA(I) 
      ENDDO
      DT=2.D0*PI/(DEX*DBLE(NP)) 
      DO JJ=1,L 
        AA=DBLE(2*JJ)-0.5D0
        T=0.D0
        DO I=1,NH
          XA(I)=XA(I)/CMPLX(AA,T,8)
          T=T+DT
        ENDDO
      ENDDO
      DO I=NHP,NP
        XA(I)=0.D0
      ENDDO
      CALL NLOGN(N,XA,WW,NP)
      CM=SQRT(2.D0*PI)/DBLE(NP) 
      AA=CM*G1**(L-1.5)
      BB=EXP((L-1.5)*DEX) 
      DO I=1,NP 
        G(I)=REAL(AA*XA(I) )
        AA=AA*BB 
      ENDDO
      RETURN
      END SUBROUTINE BTSMALLG
!     ..................................................................
      SUBROUTINE SIEGMANN(L,N,R1,DEX,NP,F,G,WW,TRBE)
!     **                                                              **
!     **  BESSELTRANSFORM USING SIEGMANS METHOD                       **
!     **  ACCURATE FOR SMALL G FOR L=0 OR L=1                         **
!     **  THIS ROUTINE CONTAINED IN MODULE BESSELTRANSFORM_MODULE     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: L
      INTEGER(4) ,INTENT(IN) :: N
      INTEGER(4) ,INTENT(IN) :: NP
      REAL(8)    ,INTENT(IN) :: R1
      REAL(8)    ,INTENT(IN) :: DEX
      REAL(8)    ,INTENT(IN) :: F(NP)
      REAL(8)    ,INTENT(OUT):: G(NP)
      COMPLEX(8),INTENT(IN) :: WW(NP)
      COMPLEX(8),INTENT(IN) :: TRBE(NP,2)
      COMPLEX(8)            :: XA(NP)
      INTEGER(4)             :: LP
      INTEGER(4)             :: NH,NHP
      INTEGER(4)             :: I
      COMPLEX(8)            :: Y1
      COMPLEX(8)            :: Y2
      COMPLEX(8)            :: W1
      COMPLEX(8)            :: W2
      REAL(8)                :: RI3,XEXP3
      REAL(8)                :: XX,YY
      INTEGER(4)             :: IC
      REAL(8)                :: CL
!     ******************************************************************
      NH=NP/2
      NHP=NH+1
      LP=L+1 
      DO I=NH,NP 
        XA(I)=(0.D0,0.D0)
      ENDDO
      RI3=R1**3
      XEXP3=EXP(3.D0*DEX) 
      DO I=1,NH 
        XX=RI3*F(2*I-1) 
        RI3=RI3*XEXP3 
        YY=RI3*F(2*I) 
        XA(I)=CMPLX(XX,YY,8) 
        RI3=RI3*XEXP3 
      ENDDO
      CALL NLOGN(N,XA,WW,NP) 
      Y1=TRBE(1,LP)*XA(1) 
      Y2=TRBE(1,LP)*CONJG(XA(1)) 
      XA(1)=2.D0*(Y1+Y2+CONJG(Y2-Y1)) 
      XA(NHP)=4.D0*CONJG(XA(NHP)*TRBE(NHP,LP)) 
      DO I=2,NH 
        IC=NP-I+2 
        Y1=XA(I) 
        Y2=CONJG(XA(IC)) 
        W1=Y1+Y2 
        W2=WW(I+NH)*(Y1-Y2) 
        Y1=(W1-W2)*TRBE(I,LP) 
        Y2=(W1+W2)*CONJG(TRBE(IC,LP)) 
        W1=Y1+Y2 
        W2=WW(I+NH)*(Y1-Y2) 
        XA(I)=W1+W2 
        XA(IC)=CONJG(W1-W2) 
      ENDDO
!
!     == FFT ===========================================================
      CALL NLOGN(N,XA,WW,NP) 
      CL=0.25D0*DEX/DBLE(NP)
      DO I=1,NH 
        G(2*I-1)=CL*REAL(XA(I),KIND=8) 
        G(2*I)  =CL*AIMAG(XA(I)) 
      ENDDO
      RETURN 
      END SUBROUTINE SIEGMANN
!     ..................................................................
      SUBROUTINE BTLARGEG(L,R1,G1,DEX,NP,F,G,N,WW,TA)
!     ******************************************************************
!     **                                                              **
!     **  BESSELTRANSFORM ACCURATE FOR LARGE G                        **
!     **  THIS ROUTINE CONTAINED IN MODULE BESSELTRANSFORM_MODULE     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L
      INTEGER(4),INTENT(IN) :: NP
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: G1
      REAL(8)   ,INTENT(IN) :: DEX
      REAL(8)   ,INTENT(IN) :: F(NP)
      REAL(8)   ,INTENT(OUT):: G(NP)
      COMPLEX(8),INTENT(IN) :: WW(NP)
      COMPLEX(8),INTENT(IN) :: TA(NP)
      COMPLEX(8)            :: XA(NP)
      COMPLEX(8)            :: Y1
      REAL(8)                :: PI
      REAL(8)                :: RIX,XEXPX
      INTEGER(4)             :: NH,NHP
      INTEGER(4)             :: IJ,IJK
      INTEGER(4)             :: I,JJ
      REAL(8)                :: AA,BB,T,DT,CM
!     ******************************************************************
      PI=4.D0*ATAN(1.D0)
      NH=NP/2
      NHP=NH+1

      RIX=R1**1.5D0
      XEXPX=EXP(1.5D0*DEX) 
      DO I=1,NP 
        XA(I)=RIX*F(I) 
        RIX=RIX*XEXPX 
      ENDDO
      CALL NLOGN(N,XA,WW,NP) 
      IJ=MOD(L,2) 
      IJK=MOD(L,4) 
      DO I=1,NH 
        Y1=XA(I)*TA(I+IJ*NH)
        IF (IJK.GT.1) Y1=-Y1
        XA(I)=Y1 
      ENDDO
      IF (L.NE.0) THEN
        DT=2.D0*PI/(DEX*DBLE(NP)) 
        DO JJ=1,L 
          AA=DBLE(2*JJ-L)-0.5D0
          BB=JJ-0.5D0
          T=0.D0 
          DO I=1,NH 
            XA(I)=XA(I)*CMPLX(BB,-T,8)/CMPLX(AA,T,8) 
            T=T+DT 
          ENDDO
        ENDDO
      END IF
      DO I=NHP,NP 
        XA(I)=(0.D0,0.D0)
      ENDDO
      CALL NLOGN(N,XA,WW,NP) 
      CM=SQRT(2.D0*PI)/DBLE(NP)
      AA=G1**(-1.5D0)*CM
      BB=EXP(-1.5D0*DEX) 
      DO I=1,NP 
        XA(I)=AA*XA(I) 
        AA=AA*BB 
      ENDDO
      DO I=1,NP
        G(I)=REAL(XA(I),KIND=8)
      ENDDO
      RETURN
      END SUBROUTINE BTLARGEG
!     ..................................................................
      SUBROUTINE MATCH(NP,G,GLARGE,DISC)
!     ******************************************************************
!     **                                                              **
!     **  MATCHES G AND GLARGE AT THE POINT OF MINIMUM DISCREPANCY    **
!     **  THIS ROUTINE CONTAINED IN MODULE BESSELTRANSFORM_MODULE     **
!     **                                                              **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NP
      REAL(8)   ,INTENT(IN)   :: GLARGE(NP)
      REAL(8)   ,INTENT(OUT)  :: DISC
      REAL(8)   ,INTENT(INOUT):: G(NP)
      REAL(8)                 :: D1,D2
      REAL(8)                 :: TE,AA
      INTEGER(4)              :: II,I
!     ******************************************************************
      D1=ABS(G(1)-GLARGE(1)) 
      D2=ABS(G(2)-GLARGE(2)) 
      TE=D1+D2 
      II=2 
      D1=D2 
      DO I=3,NP
        D2=ABS(G(I)-GLARGE(I)) 
        AA=D1+D2 
        IF (AA.LT.TE) THEN
          II=I 
          TE=AA 
        ENDIF
        D1=D2 
      ENDDO
      DO I=II,NP
        G(I)=GLARGE(I) 
      ENDDO
      DISC=TE 
      RETURN
      END SUBROUTINE MATCH
!     ..................................................................
      SUBROUTINE MAKETA(R1,G1,H,NP,TA)
!     ******************************************************************
!     **                                                              **
!     **  THE NEXT INSTRUCTIONS INITIALIZE THE ARRAY TA USED IN       **
!     **  THE TRANSFORM AT LARGE K VALUES AND SMALL K VALUES          **
!     **  FOR L GREATER THAN 1                                        **
!     **  THIS ROUTINE CONTAINED IN MODULE BESSELTRANSFORM_MODULE     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN)  :: R1
      REAL(8)   ,INTENT(IN)  :: G1
      REAL(8)   ,INTENT(IN)  :: H
      INTEGER(4),INTENT(IN)  :: NP
      COMPLEX(8),INTENT(OUT) :: TA(NP)
      REAL(8)                :: PI
      INTEGER(4)             :: NH    ! NP/2
      REAL(8)                :: DT,AA,T,S,XX,CC,PHI,RR
      INTEGER(4)             :: I,NA  ! DO LOOP INDICES
      COMPLEX(8)             :: Y1,Y2
      REAL(8)                :: RHOMIN
      REAL(8)                :: KAPMIN
!     ******************************************************************
      RHOMIN=LOG(R1)
      KAPMIN=LOG(G1)
      PI=4.D0*ATAN(1.D0) 
      NH=NP/2
      DT=2.D0*PI/(H*DBLE(NP))
      Y1=1.D0 
      AA=(RHOMIN+KAPMIN)*DT 
      Y2=CMPLX(COS(AA),SIN(AA),8) 
      DO I=1,NH 
        T=(I-1)*DT 
        S=0.D0 
        RR=SQRT(110.25D0+T*T) 
        PHI=ATAN(T/10.5D0) 
        DO NA=1,10 
          S=S+ATAN(T/(NA-0.5D0))
        ENDDO
        S=S-T*LOG(RR)+T-10.D0*PHI+SIN(PHI)/(12.D0*RR) 
        S=S-SIN(3.D0*PHI)/(360.D0*RR**3)+SIN(5.D0*PHI)/(1260.D0*RR**5) &
     &     -SIN(7.D0*PHI)/(1680.D0*RR**7) 
        XX=EXP(PI*T) 
        XX=ATAN((XX-1.D0)/(XX+1.D0)) 
        CC=S-XX 
        TA(I)=Y1*CMPLX(COS(CC),SIN(CC),8) 
        CC=S+XX 
        TA(I+NH)=Y1*CMPLX(COS(CC),SIN(CC),8) 
        Y1=Y1*Y2 
      ENDDO
      TA(1)=TA(1)/2.D0 
      TA(1+NH)=TA(1+NH)/2.D0 
      TA(NH)=TA(NH)/2.D0 
      TA(NP)=TA(NP)/2.D0 
      RETURN
      END SUBROUTINE MAKETA
!     ..................................................................
      SUBROUTINE MAKEWW(NP,WW)
!     ******************************************************************
!     **                                                              **
!     ** INITIALIZE THE ARRAY WW USED BY THE NLOGN SUBROUTINE.        **
!     ** THE ELEMENTS IN THE SECOND HALF OF THE ARRAY WW ARE USED     **
!     ** IN THE IMPLEMENTATION OF SIEGMANS METHOD.                    **
!     **  THIS ROUTINE CONTAINED IN MODULE BESSELTRANSFORM_MODULE     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NP
      COMPLEX(8),INTENT(OUT) :: WW(NP)
      INTEGER(4)             :: NH
      REAL(8)                :: PI
      REAL(8)                :: AN
      REAL(8)                :: XX
      INTEGER(4)             :: I
!     ******************************************************************
      NH=NP/2
      PI=4.D0*ATAN(1.D0)
      AN=DBLE(NP)
      DO I=1,NH 
        XX=(I-1)*PI/AN 
        WW(I+NH)=CMPLX(-SIN(XX),COS(XX),8) 
        XX=2.D0*XX 
        WW(I)=CMPLX(COS(XX),SIN(XX),8) 
      ENDDO
      RETURN 
      END SUBROUTINE MAKEWW
!     ..................................................................
      SUBROUTINE MAKETRBE(R1,G1,DEX,NP,N,WW,TRBE)
!     ******************************************************************
!     **                                                              **
!     **  INITIALIZE THE ARRAY TRBE USED IN THE IMPLEMENTATION        **
!     **  OF SIEGMANS METHOD.                                         **
!     **  THIS ROUTINE CONTAINED IN MODULE BESSELTRANSFORM_MODULE     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NP
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: G1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: N
      COMPLEX(8),INTENT(IN) :: WW(NP)
      COMPLEX(8),INTENT(OUT):: TRBE(NP,2)
      COMPLEX(8)            :: XA(NP)
      INTEGER(4)            :: NH,NHP  ! NP/2; NP/2+2
      REAL(8)               :: XEXP
      REAL(8)               :: XX
      REAL(8)               :: AA,BB
      REAL(8)               :: CC,CD
      INTEGER(4)            :: I
      INTEGER(4)            :: IC
      COMPLEX(8)            :: Y1,Y2,W1,W2
!     ******************************************************************
      NH=NP/2
      NHP=NH+1
! 
      XA(:)=(0.D0,0.D0)
      XEXP=EXP(DEX) 
      XX=R1*G1
      DO I=1,NP 
        AA=SIN(XX)/XX
        XX=XEXP*XX
        BB=SIN(XX)/XX
        XA(I)=CMPLX(AA,BB,8) 
        IF (XX.GT.1.D+8) EXIT
        XX=XX*XEXP
      ENDDO
!
      CALL NLOGN(N,XA,WW,NP) 
!
      TRBE(1,1)=XA(1) 
      TRBE(NHP,1)=CONJG(XA(NHP)) 
      DO I=2,NH 
        IC=NP-I+2 
        Y1=XA(I) 
        Y2=CONJG(XA(IC)) 
        W1=Y1+Y2 
        W2=WW(NH+I)*(Y1-Y2) 
        Y1=W1-W2 
        Y2=CONJG(W1+W2) 
        TRBE(I,1)=0.5D0*CONJG(Y1) 
        TRBE(IC,1)=0.5D0*CONJG(Y2) 
      ENDDO
      XA(:)=(0.D0,0.D0)
      XX=R1*G1
      DO I=1,NP 
        IF (XX.GE.0.1D0) THEN
          AA=(SIN(XX)/XX-COS(XX))/XX 
          XX=XEXP*XX 
          BB=(SIN(XX)/XX-COS(XX))/XX 
          XA(I)=CMPLX(AA,BB,8) 
          IF (XX.GT.1.D+8) EXIT
        ELSE
          CC=XX*XX/2.D0 
          CD=1.D0-CC/5.D0+CC*CC/70.D0-CC*CC*CC/1890.D0+CC**4/83160.D0 
          AA=XX*CD/3.D0 
          XX=XEXP*XX 
          CC=XX*XX/2.D0 
          CD=1.D0-CC/5.D0+CC*CC/70.D0-CC*CC*CC/1890.D0+CC**4/83160.D0 
          BB=XX*CD/3.D0 
          XA(I)=CMPLX(AA,BB,8) 
        END IF
        XX=XX*XEXP 
      ENDDO
!
      CALL NLOGN(N,XA,WW,NP) 
!
      TRBE(1,2)=XA(1) 
      TRBE(NHP,2)=CONJG(XA(NHP)) 
      DO I=2,NH 
        IC=NP-I+2 
        Y1=XA(I) 
        Y2=CONJG(XA(IC)) 
        W1=Y1+Y2 
        W2=WW(NH+I)*(Y1-Y2) 
        Y1=W1-W2 
        Y2=CONJG(W1+W2) 
        TRBE(I,2) =0.5D0*CONJG(Y1)
        TRBE(IC,2)=0.5D0*CONJG(Y2)
      ENDDO
      RETURN
      END SUBROUTINE MAKETRBE
!     ..................................................................
      SUBROUTINE NLOGN(N,X,WW,NP) 
!     ******************************************************************
!     **                                                              **
!     **  FAST FOURIER TRANSFORM ROUTINE                              **
!     **  THIS ROUTINE CONTAINED IN MODULE BESSELTRANSFORM_MODULE     **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: N
      INTEGER(4),INTENT(IN)   :: NP
      COMPLEX(8),INTENT(INOUT):: X(NP)
      INTEGER(4)              :: MM(15) 
      COMPLEX(8)              :: WW(NP)
      COMPLEX(8)              :: WK,HOLD,Q 
      INTEGER(4)              :: I,J,JH,K,II,IBLOCK,L
      INTEGER(4)              :: LX,LBLOCK,NBLOCK,LBHALF,ISTART
!     ******************************************************************
      DO I=1,N 
        MM(I)=2**(N-I) 
      ENDDO
      LX=2*MM(1) 
      DO L=1,N 
        NBLOCK=2**(L-1) 
        LBLOCK=LX/NBLOCK 
        LBHALF=LBLOCK/2 
        K=0
        DO IBLOCK=1,NBLOCK
          WK=WW(K+1) 
          ISTART=LBLOCK*(IBLOCK-1) 
          DO I=1,LBHALF 
            J=ISTART+I 
            JH=J+LBHALF 
            Q=X(JH)*WK 
            X(JH)=X(J)-Q 
            X(J)=X(J)+Q 
          ENDDO 
          DO I=2,N 
            II=I 
            IF (K.LT.MM(I)) EXIT
            K=K-MM(I) 
          ENDDO
          K=K+MM(II) 
        ENDDO
      ENDDO
      K=0
      DO J=1,LX 
        IF (K.GE.J) THEN
          HOLD= X(J) 
          X(J)=X(K+1) 
          X(K+1)= HOLD 
        ENDIF
        DO I=1,N 
          II=I 
          IF (K.LT.MM(I)) EXIT
          K=K-MM(I) 
        ENDDO
        K=K+MM(II) 
      ENDDO
      RETURN 
      END SUBROUTINE NLOGN
END MODULE BESSELTRANSFORM_MODULE
!
!     ..................................................................
      SUBROUTINE BESSELTRANSFORM$CLEAR
!     ******************************************************************
!     **                                                              **
!     **  DEALLOCATES TEMPORARY STORAGE FROM BESSELTRANSFORM          **
!     **                                                              **
!     ******************************************************************
      USE BESSELTRANSFORM_MODULE
      IF(ALLOCATED(TA))   DEALLOCATE(TA)
      IF(ALLOCATED(TRBE)) DEALLOCATE(TRBE)
      IF(ALLOCATED(WW))   DEALLOCATE(WW)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE RADIAL$BESSELTRANSFORM(L,GID1,NR_,F,GID2,NG_,G)
!     **                                                              **
      USE BESSELTRANSFORM_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L      ! ANGULAR MOMENTUM
      INTEGER(4),INTENT(IN) :: GID1
      INTEGER(4),INTENT(IN) :: NR_
      REAL(8)   ,INTENT(IN) :: F(NR_)
      INTEGER(4),INTENT(IN) :: GID2
      INTEGER(4),INTENT(IN) :: NG_
      REAL(8)   ,INTENT(OUT):: G(NG_)
      REAL(8)               :: DISC
      INTEGER(4)            :: GID1A
      INTEGER(4)            :: GID2A
      INTEGER(4)            :: I
      INTEGER(4)            :: NX   ! DIMENSION OF SUPPORT GRIDS
      REAL(8)               :: DEX  ! LOGARITMIC SPACING OF SUPPORT GRIDS
      REAL(8)               :: R1A  ! 
      REAL(8)               :: G1A  ! 
      REAL(8)               :: G1  ! 
      REAL(8)               :: R(NR_) 
      REAL(8)   ,ALLOCATABLE :: FA(:),GA(:)
!     ******************************************************************
                             CALL TRACE$PUSH('RADIAL$BESSELTRANSFORM')
!
!     ==================================================================
!     == DETERMINE SUPPORT GRIDS                                      ==
!     ==================================================================
      GID1A=-1
      GID2A=-1
      DO I=1,NMAP
        IF(MAP(1,I).EQ.GID1.AND.MAP(2,I).EQ.GID2) THEN
          GID1A=MAP(3,I)
          GID2A=MAP(4,I)
        END IF
      ENDDO
!
!     ==================================================================
!     == DEFINE NEW SUPPORT GRIDS                                     ==
!     ==================================================================
      IF(GID1A.LT.0) THEN
        IF(NMAP.GE.NMAPX) THEN
          CALL ERROR$MSG('OUT OF RANGE')
          CALL ERROR$STOP('RADIAL$BESSELTRANSFORM')
        END IF
        CALL RADIAL$NEW('LOG',GID1A)
        CALL RADIAL$NEW('LOG',GID2A)
        NMAP=NMAP+1
        MAP(1,NMAP)=GID1
        MAP(2,NMAP)=GID2
        MAP(3,NMAP)=GID1A
        MAP(4,NMAP)=GID2A
        NX=INT(LOG(REAL(MAX(NR_,NG_),KIND=8))/LOG(2.D0)+1.D0-1.D-8)
        NX=2**NX
        CALL RADIAL$SETI4(GID1A,'NR',NX)
        CALL RADIAL$SETI4(GID2A,'NR',NX)
!
        CALL RADIAL$GETR8(GID2,'DEX',DEX)
        CALL RADIAL$SETR8(GID1A,'DEX',DEX)
        CALL RADIAL$SETR8(GID2A,'DEX',DEX)
!
        CALL RADIAL$GETR8(GID2,'R1',G1)
        CALL RADIAL$SETR8(GID2A,'R1',G1)
        CALL RADIAL$R(GID1,NR_,R)
        R1A=R(NR_)/EXP(DEX*(NX-1))
        CALL RADIAL$SETR8(GID1A,'R1',R1A)
      END IF
!
!     ==================================================================
!     == TRANSFORM TO SUPPORT GRID, BESSELTRANSFORM AND TRANSFORM FROM SUPPORT GRID==
!     ==================================================================
      CALL RADIAL$GETI4(GID1A,'NR',NX)
      ALLOCATE(FA(NX))
      ALLOCATE(GA(NX))
      CALL RADIAL$GETR8(GID1A,'DEX',DEX)
      CALL RADIAL$GETR8(GID1A,'R1',R1A)
      CALL RADIAL$GETR8(GID2A,'R1',G1A)
      CALL RADIAL$CHANGEGRID(GID1,NR_,F,GID1A,NX,FA)
      CALL BESSELTRANSFORM(L,NX,R1A,G1A,DEX,FA,GA,DISC)
      CALL RADIAL$CHANGEGRID(GID2A,NX,GA,GID2,NG_,G)
      DEALLOCATE(FA)
      DEALLOCATE(GA)
!CALL ERROR$STOP('FORCED STOP IN RADIAL$BESSELTRANSFORM')
                                    CALL TRACE$POP
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE BESSELTRANSFORM(L,NP,R1,G1,DEX,F,G,DISC)
!     **                                                              **
!     **  CALCULATES THE SPHERICAL BESSEL TRANSFORM OF ORDER L        **
!     **    F_L(G)=INT_0^INFTY : X**2*F_L(R)*J_L(|G|*|R|)             **
!     **  FOR A FUNCTION F GIVEN ON THE LOGARITHMIC GRID              **
!     **            R(I)=R1*EXP((I-1)*DEX)  ; I=1,.,NP                **
!     ** IN THE PRESENT IMPLEMENTATION NP MUST BE A POWER OF 2        **
!     ** THE RESULT IS RETURNED AS FOFG ON THE LOGARITHMIC GRID       **
!     **            G(I)=G1*EXP((I-1)*H)    ;I=1,NP                   **
!     **                                                              **
!     ** METHOD:  J.D. TALMAN. COMP. PHYS. COMMUN. 30 (1983) 93       **
!     **                                                              **
!     **   TWO TRANSFORMS ARE CALCULATED, ONE OF WHICH IS ACCURATE    **
!     **   SMALL G AND THE OTHER FOR LARGE G. THE TWO SOLUTIONS       **
!     **   ARE JOINT AT THE POINT OF MINIMUM DISCREPANCY.             **
!     **   THE DISCREPANCY AT THE MATCHING POINT IS RETURNED AS       **
!     **   VARIABLE DISC SERVING AS ERROR ESTIMATE.                   **
!     **                                                              **
!     **   FOR L=0,1 THE SMALL-G TRANSFORM USES SIEGMAN'S METHOD      **
!     **                                                              **
!     **   THE FOURIERTRANSFORM ACCORDING TO THE CONVENTION USED      **
!     **   IN PAW IS OBTAINED BY MULTIPLICATION WITH 4*PI*I**L/VOL    **
!     **                                                              **
!     ******************************************************************
      USE BESSELTRANSFORM_MODULE
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN) :: L       ! ANGULAR MOMENTUM
      REAL(8)    ,INTENT(IN) :: R1      ! FIRST R-GRID POINT
      REAL(8)    ,INTENT(IN) :: G1      ! FIRST G-GRID POINT 
      REAL(8)    ,INTENT(IN) :: DEX     ! DEX
      INTEGER(4) ,INTENT(IN) :: NP      ! #(GRID POINTS ON THE RADIAL GRID)
      REAL(8)    ,INTENT(IN) :: F(NP)   ! ARRAY TO BE TRANSFORMED
      REAL(8)    ,INTENT(OUT):: G(NP)   ! TRANSFORMED ARRAY
      REAL(8)    ,INTENT(OUT):: DISC    ! ESTIMATED ABSOLUTE ACCURACY OF G
      INTEGER(4)             :: N       ! 2**N GRID POINTS ON THE RADIAL GRID
      LOGICAL(4)             :: CHANGED
      REAL(8)                :: GLARGE(NP)
      INTEGER(4)             :: ISVAR
!     ******************************************************************
!
!     ==================================================================
!     ==  CHECK IF GRIDS HAVE CHANGED                                 ==
!     ==================================================================
      CHANGED=.FALSE.
      CHANGED=CHANGED.OR.(NP.NE.NPSAVE)
      CHANGED=CHANGED.OR.(R1.NE.R1SAVE)
      CHANGED=CHANGED.OR.(G1.NE.G1SAVE)
      CHANGED=CHANGED.OR.(DEX.NE.DEXSAVE)
!
!     ==================================================================
!     ==  UPDATE PERMANENT INFORMATION                                ==
!     ==================================================================
      IF(CHANGED) THEN
        NPSAVE=NP
        R1SAVE=R1
        G1SAVE=G1
        DEXSAVE=DEX
!       ==================================================================
!       ==  CALCULATE N WITH 2**N=NP                                    ==
!       ==================================================================
        ISVAR=2
        N=1
        DO WHILE (ISVAR.NE.NP)
          N=N+1
          ISVAR=ISVAR*2
          IF(ISVAR.GT.NP) THEN
            CALL ERROR$MSG('NP MUST BE A POWER OF TWO')
            CALL ERROR$I4VAL('NP',NP)
            CALL ERROR$STOP('BESSELTRANSFORM')
          END IF
        ENDDO
        NSAVE=N
!
!       ==================================================================
!       ==  ALLOCATE ARRAYS                                             ==
!       ==================================================================
        IF(ALLOCATED(TA)) DEALLOCATE(TA)
        ALLOCATE(TA(NP))
        IF(ALLOCATED(TRBE)) DEALLOCATE(TRBE)
        ALLOCATE(TRBE(NP,2))
        IF(ALLOCATED(WW)) DEALLOCATE(WW)
        ALLOCATE(WW(NP))
! 
!       ==================================================================
!       ==  INITIALIZE AUXILIARY DATA FOR TRANSFORM                     ==
!       ==================================================================
!       ==  TA IS USED  AT LARGE G AND AT SMALL G FOR L>1 ==============
        CALL MAKETA(R1,G1,DEX,NP,TA)
! 
!       ==  WW IS USED IN THE FFT  =====================================
!       ==  WW(NH:NP) IS USED FOR SIEGMANS METHOD  =====================
        CALL MAKEWW(NP,WW)
! 
!       == TRBE IS USED FOR SIEGMANS METHOD. ===========================
        CALL MAKETRBE(R1,G1,DEX,NP,N,WW,TRBE)
      ELSE
        N=NSAVE
      END IF
! 
!     ==================================================================
!     == CALCULATE RESULTS ACCURATE AT SMALL K VALUES                 ==
!     ==================================================================
      IF (L.GE.2) THEN
        CALL BTSMALLG(L,R1,G1,DEX,NP,F,G,N,WW,TA)
      ELSE
        CALL SIEGMANN(L,N,R1,DEX,NP,F,G,WW,TRBE)
      END IF
! 
!     =================================================================
!     == CALCULATE XA ACCURATE AT LARGE K VALUES                     ==
!     =================================================================
      CALL BTLARGEG(L,R1,G1,DEX,NP,F,GLARGE,N,WW,TA)
! 
!     =================================================================
!     == MATCH THE TWO RESULTS AT THE POINT OF MINIMUM DISCREPANCY   ==
!     =================================================================
      CALL MATCH(NP,G,GLARGE,DISC)
      RETURN
      END
