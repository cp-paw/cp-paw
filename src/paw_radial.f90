!***********************************************************************
!**                                                                   **
!**   TEST ROUTINE FOR RADIAL OBJECT                                  **
!**                                                                   **
!**   NOTE THAT THE INTEGRATE ROUTINE OF THE OLD OBJECT WAS INACCURATE**
!**   BECAUSE OF ONLY QUADRATIC INERPOLATION TO THE ORIGIN            **
!**                                                                   **
!***********************************************************************
!    
!     ................................................................. 
      SUBROUTINE RADIAL$TEST()
      IMPLICIT NONE
      INTEGER(4)           :: GID,GID1,GID2
      INTEGER(4),PARAMETER :: NR=250
      REAL(8)              :: R1=1.056D-4,DEX=0.05D0
      REAL(8)              :: RARR(NR)
      REAL(8)              :: F(NR),DFDR(NR),INTF(NR)
      REAL(8)              :: NUMINTF(NR),NUMDFDR(NR)
      REAL(8)              :: XEXP
      REAL(8)              :: RI
      REAL(8)              :: X,VAL,VAL1,DER,INT
      REAL(8)              :: PI
      INTEGER(4)           :: IR
      INTEGER(4)           :: IGRID
      CHARACTER(10)        :: STRING
!     ****************************************************************
      PI=4.D0*DATAN(1.D0)
      CALL RADIAL$NEW('LOG',GID1)
      CALL RADIAL$SETI4(GID1,'NR',NR)
      CALL RADIAL$SETR8(GID1,'DEX',DEX)
      CALL RADIAL$SETR8(GID1,'R1',R1)
      CALL RADIAL$NEW('SHLOG',GID2)
      CALL RADIAL$SETI4(GID2,'NR',NR)
      CALL RADIAL$SETR8(GID2,'DEX',DEX)
      CALL RADIAL$SETR8(GID2,'R1',R1)
      DO IGRID=1,3
        WRITE(*,FMT='("==============================================")')
        WRITE(*,FMT='("==============================================")')
        IF(IGRID.EQ.1) THEN
          GID=GID2
          STRING='SHLOG'
          WRITE(*,FMT='("SHIFTED LOGARITHMIC GRID ",I5)')GID
        ELSE IF(IGRID.EQ.2) THEN
          GID=GID1
          STRING='LOG'
          WRITE(*,FMT='("LOGARITHMIC GRID ",I5)')GID
        ELSE IF(IGRID.EQ.3) THEN
          GID=0
          STRING='OLD'
          WRITE(*,FMT='("OLD RADIAL OBJECT ",I5)')GID
        END IF
        WRITE(*,FMT='("==============================================")')
        WRITE(*,FMT='("==============================================")')
        CALL RADIAL$R(GID,NR,RARR)
!
!       ==================================================================
!       == TEST INTEGRATE AND DERIVATIVE FOR THE LOGARITHMIC GRID      ===
!       ==================================================================
!       == GAUSS FUNCTION
        WRITE(*,FMT='("==============================================")')
        WRITE(*,FMT='("TESTS ON GAUSS FUNCTUION ")')
        WRITE(*,FMT='("==============================================")')
        DO IR=1,NR
          CALL F1(RARR(IR),F(IR),DFDR(IR),INTF(IR))
        ENDDO
        CALL RADIAL$INTEGRATE(GID,NR,F,NUMINTF)
        CALL RADIAL$DERIVE(GID,NR,F,NUMDFDR)
        PRINT*,' MAX DEV. DERIVATIVE ', MAXVAL(ABS(NUMDFDR-DFDR))
        PRINT*,' MAX DEV. INTEGRATE  ', MAXVAL(ABS(NUMINTF-INTF))
        OPEN(100,FILE='TEST-GAUSS-'//TRIM(STRING)//'.DAT')
        DO IR=1,NR
          WRITE(100,FMT='(6F25.15)')RARR(IR),F(IR),DFDR(IR),NUMDFDR(IR)-DFDR(IR) &
       &                            ,INTF(IR),NUMINTF(IR)-INTF(IR)
        ENDDO
        CLOSE(100)
        CALL RADIAL$INTEGRAL(GID,NR,F,VAL)
        PRINT*,'TEST INTEGRAL        ',VAL,VAL-INTF(NR)
        X=1.D0
        CALL F1(X,VAL,DER,INT)
        CALL RADIAL$VALUE(GID,NR,F,1.D0,VAL1)
        PRINT*,'TEST VALUE           ',VAL1,VAL1-VAL
        CALL RADIAL$DERIVATIVE(GID,NR,F,1.D0,VAL1)
        PRINT*,'TEST DERIVATIVE      ',VAL1,VAL1-DER
!   
!       ================================================================
!       == SINUS FUNCTION
!       ================================================================
        WRITE(*,FMT='("==============================================")')
        WRITE(*,FMT='("TESTS ON SINUS FUNCTUION ")')
        WRITE(*,FMT='("==============================================")')
        DO IR=1,NR
          CALL F2(RARR(IR),F(IR),DFDR(IR),INTF(IR))
        ENDDO
        CALL RADIAL$INTEGRATE(GID,NR,F,NUMINTF)
        CALL RADIAL$DERIVE(GID,NR,F,NUMDFDR)
        PRINT*,' MAX DEV. DERIVE     ', MAXVAL(ABS(NUMDFDR-DFDR))
        PRINT*,' MAX DEV. INTEGRATE  ', MAXVAL(ABS(NUMINTF-INTF))
        OPEN(100,FILE='TEST-SIN-'//TRIM(STRING)//'.DAT')
        DO IR=1,NR
          WRITE(100,*)RARR(IR),F(IR),DFDR(IR),NUMDFDR(IR)-DFDR(IR) &
       &                             ,INTF(IR),NUMINTF(IR)-INTF(IR)
        ENDDO
        CLOSE(100)
        X=1.D0
        CALL F2(X,VAL,DER,INT)
        CALL RADIAL$VALUE(GID,NR,F,1.D0,VAL1)
        PRINT*,'TEST VALUE           ',VAL1,VAL1-VAL
        CALL RADIAL$DERIVATIVE(GID,NR,F,1.D0,VAL1)
        PRINT*,'TEST DERIVATIVE      ',VAL1,VAL1-DER
!   
!       ================================================================
!       == COSINUS FUNCTION
!       ================================================================
        WRITE(*,FMT='("==============================================")')
        WRITE(*,FMT='("TESTS ON COSINUS FUNCTUION ")')
        WRITE(*,FMT='("==============================================")')
        DO IR=1,NR
          CALL F3(RARR(IR),F(IR),DFDR(IR),INTF(IR))
        ENDDO
        CALL RADIAL$INTEGRATE(GID,NR,F,NUMINTF)
        CALL RADIAL$DERIVE(GID,NR,F,NUMDFDR)
        PRINT*,' MAX DEV. DERIVATIVE ', MAXVAL(ABS(NUMDFDR-DFDR))
        PRINT*,' MAX DEV. INTEGRATE  ', MAXVAL(ABS(NUMINTF-INTF))
        OPEN(100,FILE='TEST-COS-'//TRIM(STRING)//'.DAT')
        DO IR=1,NR
          WRITE(100,*)RARR(IR),F(IR),DFDR(IR),NUMDFDR(IR)-DFDR(IR) &
       &                             ,INTF(IR),NUMINTF(IR)-INTF(IR)
        ENDDO
        CLOSE(100)
        X=1.D0
        CALL F3(X,VAL,DER,INT)
        CALL RADIAL$VALUE(GID,NR,F,1.D0,VAL1)
        PRINT*,'TEST VALUE          ',VAL1,VAL,VAL1-VAL
        CALL RADIAL$DERIVATIVE(GID,NR,F,1.D0,VAL1)
        PRINT*,'TEST DERIVATIVE     ',VAL1,DER,VAL1-DER
      ENDDO
      RETURN
      CONTAINS
!  
!       ................................................................
        SUBROUTINE F1(R,F,DFDR,INTF)
        REAL(8),INTENT(IN)  :: R
        REAL(8),INTENT(OUT) :: F
        REAL(8),INTENT(OUT) :: DFDR
        REAL(8),INTENT(OUT) :: INTF
        REAL(8)             :: PI,SVAR
!       ****************************************************************
        PI=4.D0*DATAN(1.D0)
        F=EXP(-R**2)
        DFDR=-2.D0*R*F
        CALL SPECIALFUNCTION$ERF(R,SVAR)
        INTF=0.5D0*SQRT(PI)*SVAR
        RETURN
      END SUBROUTINE F1
!
!       ................................................................
        SUBROUTINE F2(R,F,DFDR,INTF)
        REAL(8),INTENT(IN)  :: R
        REAL(8),INTENT(OUT) :: F
        REAL(8),INTENT(OUT) :: DFDR
        REAL(8),INTENT(OUT) :: INTF
        REAL(8)             :: PI
!       ****************************************************************
        PI=4.D0*DATAN(1.D0)
        F=SIN(2.D0*R)
        DFDR=2.D0*COS(2.D0*R)
        INTF=-0.5D0*(COS(2.D0*R)-1.D0)
        RETURN
      END SUBROUTINE F2
!
!       ................................................................
        SUBROUTINE F3(R,F,DFDR,INTF)
        REAL(8),INTENT(IN)  :: R
        REAL(8),INTENT(OUT) :: F
        REAL(8),INTENT(OUT) :: DFDR
        REAL(8),INTENT(OUT) :: INTF
        REAL(8)             :: PI
!       ****************************************************************
        PI=4.D0*DATAN(1.D0)
        F=COS(2.D0*R)
        DFDR=-2.D0*SIN(2.D0*R)
        INTF=0.5D0*SIN(2.D0*R)
        RETURN
      END SUBROUTINE F3
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
      INTEGER(4)            :: I,J,IP,K
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
      INTEGER(4)             :: ITEST
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
      SUBROUTINE RADIAL_DGLEQUISPACEDgen(NX,nf,idir,A,B,C,D,F)
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
      INTEGER(4) ,INTENT(IN) :: Nf
      REAL(8)    ,INTENT(IN) :: A(NX)
      REAL(8)    ,INTENT(IN) :: B(NX)
      REAL(8)    ,INTENT(IN) :: C(NX,nf,nf)
      REAL(8)    ,INTENT(IN) :: D(NX,nf)
      REAL(8)    ,INTENT(INOUT):: F(NX,nf)
      REAL(8)                :: AP(NX),A0(NX),AM(NX)
      INTEGER(4)             :: I
      INTEGER(4)             :: ITEST
!     ******************************************************************
!     == AP*F(+) + A0*F(0) + AM*F(-) = D
      AP(:)=A(:)+0.5D0*B(:)
      A0(:)=-2.D0*A(:)
      AM(:)=A(:)-0.5D0*B(:)
!print*,'ap ',ap
!print*,'a0 ',a0
!print*,'am ',am
!print*,'d ',d
!print*,'c ',c
      IF(IDIR.GE.0) THEN
        DO I=2,NX-1
          F(I+1,:)=( -(A0(I)+AM(I))*F(I,:) +AM(I)*(F(I,:)-F(I-1,:))-matmul(c(i,:,:),f(i,:)) +D(I,:) )/AP(I)
        ENDDO
      ELSE IF(IDIR.LT.0) THEN
        DO I=NX-1,2,-1
          F(I-1,:)=( -(A0(I)+AP(I))*F(I,:) +AP(I)*(F(I,:)-F(I+1,:))-matmul(c(i,:,:),f(i,:)) +D(I,:) )/AM(I)
        ENDDO
      ELSE
         CALL ERROR$MSG('INVALID VALUE OF IDIR')
         CALL ERROR$STOP('RADIAL_DGLEQUISPACED')
      END IF
!print*,'f ',f
!stop
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
          GRIDARRAY(:,NGIDX)=NEWARRAY(:,NGIDX)
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
!     ******************************************************************
      CALL RADIAL$R(GID2,NR2,R2)
      DO IR=1,NR2
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
  IF(NR.NE.250) STOP 'ERROR IN RADIAL$VALUE/ GID=0'
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
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$I4VAL('GID',GID)
        CALL ERROR$I4VAL('TYPE',TYPE)
        CALL ERROR$STOP('RADIAL$VALUE')
      END IF  
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$INTEGRATE(GID,NR,F,G)
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
        CALL ERROR$STOP('RADIAL$DERIVATIVE')
      END IF  
      RETURN
      END      
!
!     ...................................................................
      SUBROUTINE RADIAL$DGLgen(GID,NR,nf,idir,A,B,C,D,F)
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: nf
      INTEGER(4),INTENT(IN) :: idir
      REAL(8)   ,INTENT(IN) :: A(NR)
      REAL(8)   ,INTENT(IN) :: B(NR)
      REAL(8)   ,INTENT(IN) :: C(NR,nf,nf)
      REAL(8)   ,INTENT(IN) :: D(NR,nf)
      REAL(8)   ,INTENT(INOUT):: F(NR,nf)
      INTEGER(4)            :: GIDS
      INTEGER(4)            :: TYPE
!     ******************************************************************
      CALL RADIAL_RESOLVE(GID,GIDS,TYPE)
      IF(TYPE.EQ.1) THEN
        CALL ERROR$MSG('not implemented for this GRID TYPE')
        CALL ERROR$STOP('RADIAL$DERIVATIVE')
!        CALL LOGRADIAL$DGLgen(GIDS,NR,nf,idir,A,B,C,D,F)
      ELSE IF(TYPE.EQ.2) THEN
        CALL SHLOGRADIAL$DGLgen(GIDS,NR,nf,idir,A,B,C,D,F)
      ELSE
        CALL ERROR$MSG('GRID TYPE NOT RECOGNIZED')
        CALL ERROR$STOP('RADIAL$DERIVATIVE')
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
!     ..................................................................
      SUBROUTINE RADIAL$SCHRODINGER(GID,NR,POT,DREL,SO,G,L,E,IDIR,PHI)
!     **                                                                  **
!     **  SOLVES THE RELATIVISTIC RADIAL DIRAC EQUATION FOR THE           **
!     **  LARGE COMPONENT. DREL=1/MREL-1/M0 IS A MEASURE FOR THE          **
!     **  RELATIVISTIC EFFECTS, WHERE MREL=M0+(E-V)/2C^2 AND M0 IS THE    **
!     **  REST MASS. V=POT*Y0 IS THE POTENTIAL.                           **
!     **  SPIN ORBIT COUPLING IS MEASURED BY SO WHICH CAN HAVE THE VALUES:**
!     **    SO=0 NO-SPIN ORBIT COUPLING                                   **
!     **    SO=1 PARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM               **
!     **    SO=-1 ANTIPARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM          **
!     **  IDIR DESCRIBES THE DIRECTION OF THE INTEGRATION OF THE DIFF.EQ. **
!     **    IDIR=1 OUTWARD INTEGRATION                                    **
!     **    IDIR=-1 INWARD INTEGRATION                                    **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE        **
!     **  HOMOGENEOUS SOLUTION.                                           **
!     **                                                                  **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM   **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                     **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                               **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                         **
!     **                                                                  **
!     **  IN THE PRESENCE OF AN INHOMOGENEITY, THE SOLUTION STARTS        **
!     **  WITH ZERO VALUE AND DERIVATIVE.                                 **
!     **                                                                  **
!     **  IN THE ABSENCE OF AN INHOMOGENEITY, THE SOLUTION STARTS         **
!     **  WITH R**L FROM THE INSIDE AND FROM THE OUTSIDE WITH VALUE ZERO  **
!     **  AND FINITE SLOPE                                                **
!     **                                                                  **
!     **  THE NONRELATIVISTIC SOLUTION IS OBTAINED BY SETTINH DREL=0      **
!     **                                                                  **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE      **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                   **
!     **                                                                  **
!     **  REMARKS:                                                        **
!     **  - POT IS ONLY THE RADIAL PART OF THE POTENTIAL.                 **
!     **    THE POTENTIAL IS POT*Y0 WHERE Y0 IS A SPHERICAL HARMONIC      **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! MAINANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
                 ! SO=0: NO SO; SO=1: L/S PARALLEL; SO=-1: L,S ANTIPARALLEL
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: G(NR)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: E       !ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: IDIR    ! IDIR=1 INTEGRATE OUTWARD
                                            ! IDIR=-1 INTEGRATE INWARD
      REAL(8)    ,INTENT(OUT)    :: PHI(NR) !WAVE-FUNCTION
      REAL(8)                    :: A(NR) 
      REAL(8)                    :: B(NR) 
      REAL(8)                    :: C(NR) 
      REAL(8)                    :: D(NR) 
      REAL(8)                    :: R(NR) 
      REAL(8)                    :: PI
      REAL(8)                    :: Y0
      REAL(8)                    :: SOFACTOR
      REAL(8)                    :: RDPRIME(NR)
      INTEGER(4)                 :: IR
      REAL(8)                    :: RI
      LOGICAL(4)                 :: THOM
!     ************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
!     ==================================================================
!     == SPIN ORBIT COUPLING ===========================================
!     ==================================================================
      IF (SO.EQ.0) THEN
        SOFACTOR=0.D0
      ELSE IF(SO.EQ.1) THEN
        SOFACTOR=REAL(L,KIND=8)
      ELSE IF(SO.EQ.-1) THEN
        SOFACTOR=REAL(-L-1,KIND=8)
      ELSE
         CALL ERROR$MSG('SO CAN ONLY HAVE VALUES -1,0,1')
         CALL ERROR$STOP('RADIAL$SCHRODINGER')
      END IF
!     ==================================================================
!     == PREPARE ARRAYS FOR INTEGRATION ================================
!     ==================================================================
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      CALL RADIAL$R(GID,NR,R)
      A(:)=1.D0+DREL(:)
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0*(1.D0+DREL(2:))/R(2:)+RDPRIME(2:)
      C(2:)=-(1.D0+DREL(2:))*REAL(L*(L+1),KIND=8)/R(2:)**2 &
     &    +RDPRIME(2:)*SOFACTOR/R(2:) &
     &    -2.D0*(POT(2:)*Y0-E)
      B(1)=B(2)
      C(1)=C(2)
      D(:)=-2.D0*G(:)
      THOM=MAXVAL(ABS(G(:))).EQ.0.D0
      IF(IDIR.GE.0) THEN
        PHI(1:2)=0.D0
        IF(THOM)PHI(1:2)=R(1:2)**L
      ELSE
        PHI(NR-1:NR)=0.D0
        IF(THOM)PHI(NR-1)=1.D-8
      END IF
      CALL RADIAL$DGL(GID,IDIR,NR,A,B,C,D,PHI)
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE RADIAL$nonsphbound(GID,NR,nspin,lmrx,POT,DREL,G,E,PHI)
!     **                                                                  **
!     **  SOLVES THE RELATIVISTIC RADIAL DIRAC EQUATION FOR THE           **
!     **  LARGE COMPONENT. DREL=1/MREL-1/M0 IS A MEASURE FOR THE          **
!     **  RELATIVISTIC EFFECTS, WHERE MREL=M0+(E-V)/2C^2 AND M0 IS THE    **
!     **  REST MASS. V=POT*Y0 IS THE POTENTIAL.                           **
!     **  SPIN ORBIT COUPLING IS MEASURED BY SO WHICH CAN HAVE THE VALUES:**
!     **    SO=0 NO-SPIN ORBIT COUPLING                                   **
!     **    SO=1 PARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM               **
!     **    SO=-1 ANTIPARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM          **
!     **  IDIR DESCRIBES THE DIRECTION OF THE INTEGRATION OF THE DIFF.EQ. **
!     **    IDIR=1 OUTWARD INTEGRATION                                    **
!     **    IDIR=-1 INWARD INTEGRATION                                    **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE        **
!     **  HOMOGENEOUS SOLUTION.                                           **
!     **                                                                  **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM   **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                     **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                               **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                         **
!     **                                                                  **
!     **  IN THE PRESENCE OF AN INHOMOGENEITY, THE SOLUTION STARTS        **
!     **  WITH ZERO VALUE AND DERIVATIVE.                                 **
!     **                                                                  **
!     **  IN THE ABSENCE OF AN INHOMOGENEITY, THE SOLUTION STARTS         **
!     **  WITH R**L FROM THE INSIDE AND FROM THE OUTSIDE WITH VALUE ZERO  **
!     **  AND FINITE SLOPE                                                **
!     **                                                                  **
!     **  THE NONRELATIVISTIC SOLUTION IS OBTAINED BY SETTINH DREL=0      **
!     **                                                                  **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE      **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                   **
!     **                                                                  **
!     **  REMARKS:                                                        **
!     **  - POT IS ONLY THE RADIAL PART OF THE POTENTIAL.                 **
!     **    THE POTENTIAL IS POT*Y0 WHERE Y0 IS A SPHERICAL HARMONIC      **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: nspin
      INTEGER(4) ,INTENT(IN)     :: lmrx
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: G(NR,lmrx,nspin)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: E       !ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR,lmrx,nspin*nspin) !POTENTIAL (RADIAL PART ONLY)
      REAL(8)    ,INTENT(OUT)    :: PHI(NR,lmrx,nspin,lmrx*nspin) ! WAVE-FUNCTION
      integer(4)                 :: lx
      INTEGER(4)                 :: lox(lmrx)
      integer(4)                 :: lm,lm1,lm2,lm3,l,m,im,is,is1,is2
      REAL(8)                    :: A(NR)
      REAL(8)                    :: B(NR)
      REAL(8)                    :: C(NR,lmrx,nspin,lmrx,nspin)
      REAL(8)                    :: D(NR,lmrx,nspin) 
      REAL(8)                    :: R(NR) 
      REAL(8)                    :: aux(NR) 
      REAL(8)                    :: PI
      REAL(8)                    :: Y0
      REAL(8)                    :: cg
      REAL(8)                    :: RDPRIME(NR)
      INTEGER(4)                 :: IR
      REAL(8)                    :: RI
      LOGICAL(4)                 :: THOM
      real(8)   ,parameter       :: xmax=1.d+100
      real(8)                    :: rcl,rtest,swkb(nr),aux1(nr),svar
      integer(4)                 :: ircl,irout
!     ************************************************************************
print*,'in radial$nonsphbound'
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/DSQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      lx=int(sqrt(lmrx-1.d0)+0.1d0)
      lm=0
      do l=0,lx
        do m=1,2*l+1
          lm=lm+1
          lox(lm)=l
        enddo
      enddo
!
!     ==================================================================
!     ==  determine classical truning point                           ==
!     ==================================================================
      RCL=R(NR)
      IF(E.LT.POT(NR,1,1)*Y0) THEN
        DO IR=NR-1,1,-1
          IF(E.GT.POT(IR,1,1)*Y0) THEN
            IRCL=IR
            RCL=R(IR)-(POT(IR,1,1)-E/Y0)/(POT(IR+1,1,1)-POT(IR,1,1))*(R(IR+1)-R(IR))
            EXIT
          END IF
        ENDDO
        RTEST=RCL
      END IF
!
!     == USE WKB SOLUTION FOR THE SCHR.GL. FOR A CONSTANT POTENTIAL AND L=0
!     == TO ESTIMATE FACTOR FROM RTEST TO OUTERMOST POINT
      irout=nr
      IF(RTEST.LT.R(NR)) THEN
        AUX(:IRCL)=0.D0
        AUX(IRCL+1:)=SQRT(REAL(L*(L+1),kind=8)/R(ircl+1:)**2+2.D0*(POT(IRCL+1:,1,1)*Y0-E))
        CALL RADIAL$DERIVe(GID,NR,AUX,AUX1)
        CALL RADIAL$INTEGRATE(GID,NR,AUX,Swkb)
!        Swkb(:)=Swkb(:)-0.5D0*LOG(AUX1(:))-LOG(R(:))
        Swkb(:)=Swkb(:)-LOG(R(:))
!       == DETERMINE IROUT WHERE THE WAVE FUNCTION CAN GROW BY A FACTOR 
!       == OF XMAX FROM THE CLASSICAL TURNING POINT
        SVAR=LOG(XMAX)
        DO IR=1,NR
          IF(Swkb(IR).GT.SVAR) THEN
            IROUT=IR-1
            EXIT
          END IF
        ENDDO
      END IF
!
!     ==================================================================
!     ==  prepare potential-independent arrays                        ==
!     ==================================================================
!     == a*d2f/dr2+b*df/dr+c*f=d
      A(:)=1.D0+DREL(:)
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0*(1.D0+DREL(2:))/R(2:)+RDPRIME(2:)
      B(1)=B(2)
      D(:,:,:)=-2.D0*G(:,:,:)
     
!     ==================================================================
!     ==  Coupling between wave function components via potential     ==
!     ==================================================================
print*,'pot ',pot(:100,1,1)
      c(:,:,:,:,:)=0.d0
      do lm1=1,lmrx
        do lm2=1,lmrx
          aux(:)=0.d0
          do lm3=1,lmrx
            call clebsch(lm1,lm2,lm3,cg)
            aux=cg*pot(:,lm3,1)
            if(lm3.eq.1) aux(irout+1:)=aux(irout)
            c(:,lm1,1,lm2,1)=c(:,lm1,1,lm2,1)+aux
            if(nspin.gt.1) then
              c(:,lm1,2,lm2,2)=c(:,lm1,2,lm2,2)+aux
              if(nspin.eq.2) then
                aux=cg*pot(:,lm3,2)
                c(:,lm1,1,lm2,1)=c(:,lm1,1,lm1,1)+aux
                c(:,lm1,2,lm2,2)=c(:,lm1,2,lm1,2)-aux
              else
                call error$stop('in radial_xxx')
                aux=cg*pot(:,lm3,2)
                c(:,lm1,1,lm2,2)=c(:,lm1,1,lm1,2)+aux
                c(:,lm1,2,lm2,1)=c(:,lm1,2,lm1,1)-aux
                aux=cg*pot(:,lm3,3)
                c(:,lm1,1,lm2,2)=c(:,lm1,1,lm1,2)+aux
                c(:,lm1,2,lm2,1)=c(:,lm1,2,lm1,1)-aux
                aux=cg*pot(:,lm3,4)
                c(:,lm1,1,lm2,1)=c(:,lm1,1,lm1,1)+aux
                c(:,lm1,2,lm2,2)=c(:,lm1,2,lm1,2)-aux
              end if
            end if
          enddo
        enddo
      enddo
      c=-2.d0*C
!
      lm=0
      do l=1,lx
        aux(1)=0.d0
        aux(2:)=-(1.D0+DREL(2:))/R(2:)**2 * real(l*(l+1),kind=8)+2.d0*e
        do im=1,2*l+1
          lm=lm+1
          do is=1,nspin
            c(:,lm,is,lm,is)=c(:,lm,is,lm,is)+aux(:)
          enddo
        enddo
      enddo
!     ==  avoid zerobyzero
      c(1,:,:,:,:)=c(2,:,:,:,:)
!     ==================================================================
!     ==  add spin orbit coupling to E                                ==
!     ==================================================================
!     == spin orbit coupling missing
!!$      do l=0,lx
!!$        xl=real(l,kind=8)
!!$        do m=-l,l
!!$          xm=real(m,kind=8)
!!$          lm=lm+1
!!$          lpfac(lm)=sqrt((xl-xm)(xl+xm+1.d0))  ! (L_+)Y_{l,m}=lmpfac*Y_{l,m+1}
!!$          lmfac(lm)=sqrt((xl+xm)(xl-xm+1.d0))  ! (L_-)Y_{l,m}=lmmfac*Y_{l,m-1}
!!$          lzfac(lm)=xm                         ! (L+z)Y_{l,m}=lmmfac*Y_{l,m}
!!$!                                              ! Y_{L,m} are true spherical harmonics
!!$        enddo
!!$      enddo
!print*,'r(irout)',rcl,r(irout),pot1(nr)
!
!     ==================================================================
!     ==  determine bound states                                      ==
!     ==================================================================
print*,'before radial_xxx',nspin,lmrx
      call radial_xxx(gid,nr,lmrx*nspin,ircl,irout,lox,a,b,c,d,phi)
      return
      end
!     ..................................................................
      subroutine radial_xxx(gid,nr,nf,irmatch,irout,lox,a,b,c1,d,phi)
      implicit none
      integer(4),intent(in) :: gid
      integer(4),intent(in) :: nr
      integer(4),intent(in) :: nf
      integer(4),intent(in) :: irmatch
      integer(4),intent(in) :: irout
      integer(4),intent(in) :: lox(nf)
      real(8)   ,intent(in) :: a(nr)
      real(8)   ,intent(in) :: b(nr)
      real(8)   ,intent(in) :: c1(nr,nf,nf)
      real(8)   ,intent(in) :: d(nr,nf)
      real(8)   ,intent(out):: phi(nr,nf,nf)
      logical               :: thom
      real(8)               :: phil(nr,nf,nf)
      real(8)               :: phir(nr,nf,nf)
      real(8)               :: phil_dot(nr,nf,nf)
      real(8)               :: phir_dot(nr,nf,nf)
      real(8)               :: phir_inhom(nr,nf)
      real(8)               :: phil_inhom(nr,nf)
      integer(4)            :: if,if1,if2,ir
      real(8)               :: mat(nf,nf),matinv(nf,nf),vec(nf)
      integer(4)            :: irc
      real(8)               :: kink_inhom(nf)
      real(8)               :: kink_hom(nf,nf)
      real(8)               :: kink_dot(nf,nf)
      real(8)               :: r(nr)
      real(8)               :: dhom(nr,nf)
      real(8) ::bvec(nf),xvec(nf)
      real(8) ::bmat(nf,nf),xmat(nf,nf)
      real(8) ::aux(nr),svar      
real(8) ::c(nr,nf,nf)
integer(4) :: i,j
      real(8):: au(2,2),xu(2),bu(2)
!     ******************************************************************
c=c1
print*,'in radial_xxx'
      call radial$r(gid,nr,r)
      irc=irmatch
!
!     ==================================================================
!     ==  obtain inhomogeneous solution                               ==
!     ==================================================================
      THOM=MAXVAL(ABS(d(:,:))).EQ.0.D0
      phil_inhom(:,:)=0.d0
      phir_inhom(:,:)=0.d0
      if(.not.thom) then
        CALL RADIAL$DGLgen(GID,NR,nf,1,A,B,C,D,PHIl_inhom)
        CALL RADIAL$DGLgen(GID,NR,nf,-1,A,B,C,D,PHIr_inhom)
      end if  
!
!     ==================================================================
!     ==  obtain homogeneous solution and their phidots               ==
!     ==================================================================
print*,'rcl ',irc,r(irc),irout,r(irout),nf
print*,'lox ',lox
!print*,'c',c(:100,nf,nf)
do i=1,nf
  do j=1,nf
    if(i.eq.j) cycle
!    c(:,i,j)=0.d0
  enddo
enddo
      phil(:,:,:)=0.d0
      phir(:,:,:)=0.d0
      phil_dot(:,:,:)=0.d0
      phir_dot(:,:,:)=0.d0
      do if=1,nf
        PHIl(1:2,if,if)=R(1:2)**lox(if)
        dhom(:,:)=0.d0
        CALL RADIAL$DGLgen(GID,NR,nf,1,A,B,C,Dhom,PHIl(:,:,if))
        svar=maxval(abs(phil(irc,:,if)))
        phil(:,:,if)=phil(:,:,if)/svar
        dhom(irout,if)=1.d-8
        CALL RADIAL$DGLgen(GID,NR,nf,-1,A,B,C,Dhom,PHIr(:,:,if))
        svar=maxval(abs(phir(irc,:,if)))
        phir(:,:,if)=phir(:,:,if)/svar
        CALL RADIAL$DGLgen(GID,NR,nf,1,A,B,C,-2.d0*phil(:,:,if),PHIl_dot(:,:,if))
        CALL RADIAL$DGLgen(GID,NR,nf,-1,A,B,C,-2.d0*phir(:,:,if),PHIr_dot(:,:,if))
      enddo
!
!     ==================================================================
!     ==  make solutions continuous                                   ==
!     ==================================================================
      mat(:,:)=phir(irc,:,:)
!     == matrix is not symmetric. thus solve with singular value decomposition
!
!     == make inhomogeneous solution continuous =======================
      if(.not.thom) then
        bvec(:)=phil_inhom(irc,:)-phir_inhom(irc,:)
        call lib$matrixsolve(nf,nf,1,mat,vec,bvec)
        do if=1,nf
          phir_inhom(:,:)=phir_inhom(:,:)+phir(:,:,if)*vec(if)
        enddo
      end if
!
!     == make phidot solutions continuous =============================
      bmat(:,:)=phil_dot(irc,:,:)-phir_dot(irc,:,:)
      call lib$matrixsolvenew(nf,nf,nf,mat,xmat,bmat)
print*,'test_dot ',maxval(abs(matmul(mat,xmat)-bmat))
!!$print*,'test1 ',maxval(abs(matmul(mat,xmat)-bmat))
!!$do if=1,nf
!!$  bvec(:)=phil_dot(irc,:,if)-phir_dot(irc,:,if)
!!$  call lib$matrixsolvenew(nf,nf,1,mat,xvec,bvec)
!!$  print*,'test2 ',if,maxval(abs(matmul(mat,xvec)-bvec))
!!$enddo
!!$stop
      do if1=1,nf
        do if2=1,nf
          phir_dot(:,:,if1)=phir_dot(:,:,if1)+phir(:,:,if2)*xmat(if2,if1)
        enddo
      enddo
print*,'cont_dot',phir_dot(irc,:,:)-phil_dot(irc,:,:)
!
!     == make homogeneous solutions continuous =======================
      bmat(:,:)=phil(irc,:,:)
      call lib$matrixsolve(nf,nf,nf,mat,xmat,bmat)
print*,'test_hom ',maxval(abs(matmul(mat,xmat)-bmat))
      do ir=1,nr        
        phir(ir,:,:)=matmul(phir(ir,:,:),xmat)
      enddo
print*,'cont',phir(irc,:,:)-phil(irc,:,:)
!
!     ==================================================================
!     ==  normalize homogeneous solutions, and scale phidot accordingly=
!     ==================================================================
      do i=1,nf
        aux(:)=0.d0
        do j=1,nf
          aux(1:irc)=aux(1:irc)+phil(1:irc,j,i)**2
          aux(irc+1:)=aux(irc+1:)+phir(irc+1:,j,i)**2
        enddo
        aux(:)=aux(:)*r(:)**2
        call radial$integral(gid,nr,aux,svar)
        svar=1.d0/sqrt(svar)
        phil(:,:,i)=phil(:,:,i)*svar
        phir(:,:,i)=phir(:,:,i)*svar
        phil_dot(:,:,i)=phil_dot(:,:,i)*svar
        phir_dot(:,:,i)=phir_dot(:,:,i)*svar
      enddo
!
!     ==================================================================
!     ==  determine kinks                                             ==
!     ==================================================================
!
!     == determine kinks (verlet expression is correct)  =============
      if(.not.thom) then
        kink_inhom(:)=(phir_inhom(irc+1,:)-phir_inhom(irc-1,:)) &
     &               -(phil_inhom(irc+1,:)-phil_inhom(irc-1,:))
      end if
      kink_hom(:,:)=(phir(irc+1,:,:)-phir(irc-1,:,:)) &
     &             -(phil(irc+1,:,:)-phil(irc-1,:,:))
      kink_dot(:,:)=(phir_dot(irc+1,:,:)-phir_dot(irc-1,:,:)) &
     &             -(phil_dot(irc+1,:,:)-phil_dot(irc-1,:,:))
!
!     == map left and right solutions into one array =================
      phil_inhom(irc:,:)=phir_inhom(irc:,:)
      phil(irc:,:,:)=phir(irc:,:,:)
      phil_dot(irc:,:,:)=phir_dot(irc:,:,:)
!
print*,'phil'
print*,phil_dot(:,nf,nf)
!
!     == diagonalize (Kink_hom+epsilon*kink_dot) Not hermitean!
print*,'nf ',nf
print*,'kink_hom ',kink_hom
print*,'kink_dot ',kink_dot
      call lib$generaleigenvaluer8(nf,kink_hom,kink_dot,vec,mat)
print*,'vec ',vec
stop
!     == solve inhomogeneous problem c=ch+ci(epsilon_n)
!     == ci(epsilon_n)=(kink_hom+epsilon_n*kink_dot)**(-1)*kink_inhom
!     not done 
      call error$stop('radial_xx')
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE RADIAL$KINPHI(GID,NR,DREL,SO,L,PHI,TPHI)
!     **                                                                  **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: L       ! MAINANGULAR MOMENTUM
      INTEGER(4) ,INTENT(IN)     :: SO      ! SWITCH FOR SPIN-ORBIT COUP.
                 ! SO=0: NO SO; SO=1: L/S PARALLEL; SO=-1: L,S ANTIPARALLEL
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: PHI(NR) !WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: TPHI(NR)!KINETIC ENERGY TIMES PHI
      REAL(8)                    :: A(NR) 
      REAL(8)                    :: B(NR) 
      REAL(8)                    :: C(NR) 
      REAL(8)                    :: R(NR) 
      REAL(8)                    :: SOFACTOR
      REAL(8)                    :: RDPRIME(NR)
      REAL(8)                    :: DPHIDR(NR)
      REAL(8)                    :: D2PHIDR2(NR)
!     ************************************************************************
!     ==================================================================
!     == SPIN ORBIT COUPLING ===========================================
!     ==================================================================
      IF (SO.EQ.0) THEN
        SOFACTOR=0.D0
      ELSE IF(SO.EQ.1) THEN
        SOFACTOR=REAL(L,KIND=8)
      ELSE IF(SO.EQ.-1) THEN
        SOFACTOR=REAL(-L-1,KIND=8)
      ELSE
         CALL ERROR$MSG('SO CAN ONLY HAVE VALUES -1,0,1')
         CALL ERROR$STOP('RADIAL$SCHRODINGER')
      END IF
!     ==================================================================
!     == PREPARE ARRAYS FOR INTEGRATION ================================
!     ==================================================================
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      CALL RADIAL$R(GID,NR,R)
      A(:)=1.D0+DREL(:)
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0*(1.D0+DREL(2:))/R(2:)+RDPRIME(2:)
      C(2:)=-(1.D0+DREL(2:))*REAL(L*(L+1),KIND=8)/R(2:)**2 &
     &    +RDPRIME(2:)*SOFACTOR/R(2:) 
      B(1)=B(2)
      C(1)=C(2)
      CALL RADIAL$VERLETD1(GID,NR,PHI,DPHIDR)
      CALL RADIAL$VERLETD2(GID,NR,PHI,D2PHIDR2)
      TPHI(:)=-0.5D0*(A(:)*D2PHIDR2(:)+B(:)*DPHIDR(:)+C(:)*PHI(:))
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
      PI=4.D0*DATAN(1.D0)
      CALL RADIAL$R(GID,NR,R)
      AUX1(:)=RHO(:)*R(:)**(L+2) 
      CALL RADIAL$INTEGRATE(GID,NR,AUX1,AUX2)
      FAC=4.D0*PI/DBLE(2*L+1)
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
!     ...........................................MLTPOL.................
      SUBROUTINE RADIAL$MOMENT(GID,NR,L,RHO,QLM)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE MULTIPOLE MOMENT OF A (CHARGE) DISTRIBUTION  **
!     **  RHO(|R|)*Y_L(R).                                            **
!     **                                                              **
!     **    QLM = INT{D^3R: RHO(R) * |R|**L *Y_L(R) }                 **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1997) **
!     ******************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID      ! FIRST POINT ON RADIAL GRID
      INTEGER(4),INTENT(IN) :: NR      ! NUMBER OF RADIAL GRID POINTS
      INTEGER(4),INTENT(IN) :: L       ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: RHO(NR) ! CHARGE DENSITY
      REAL(8)   ,INTENT(OUT):: QLM     ! MULTIPOLE MOMENT
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: XEXP,RI
      INTEGER(4)            :: IR
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
      XEXP=DEXP(DEX)
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
      INTEGER(4)            :: IR,IR1
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
      SUBROUTINE SHLOGRADIAL$DGLgen(GID,NR_,nf,idir,A,B,C,D,F)
      USE SHLOGRADIAL_MODULE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR_
      INTEGER(4),INTENT(IN) :: Nf
      INTEGER(4),INTENT(IN) :: idir
      REAL(8)   ,INTENT(IN) :: A(NR_)
      REAL(8)   ,INTENT(IN) :: B(NR_)
      REAL(8)   ,INTENT(IN) :: C(NR_,nf,nf)
      REAL(8)   ,INTENT(IN) :: D(NR_,nf)
      REAL(8)   ,INTENT(INOUT):: F(NR_,nf)
      REAL(8)               :: R(NR_)
      REAL(8)               :: DRDX(NR_)
      REAL(8)               :: B1(NR_)
      REAL(8)               :: C1(NR_,nf,nf)
      REAL(8)               :: D1(NR_,nf)
      INTEGER(4)            :: IR,IR1,i,j
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
      DO I=1,NF
        D1(:,I)=D(:,I)*DRDX(:)**2
        DO J=1,NF
          C1(:,I,J)=C(:,I,J)*DRDX(:)**2
        ENDDO
      ENDDO
      CALL RADIAL_DGLEQUISPACEDGEN(NR,NF,idir,A,B1,C1,D1,F)
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
      XEXP=DEXP(DEX)
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
        STOP 'ERROR LOGRADIAL$VALUE'
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
      INTEGER(4)            :: IR,IR1
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
!*************PETER E. BLOECHL, IBM ZURICH RESEARCH LABORATORY(1996)****
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
      PI=4.D0*DATAN(1.D0)
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
        G(I)=AA*XA(I) 
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
      PI=4.D0*DATAN(1.D0)
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
!call error$stop('forced stop in radial$besseltransform')
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
!
!     .....................................................GETF.........
      SUBROUTINE GETF(NR,F,X,RES)
!     **                                                              **
!     **  CALCULATES THE VALUE RES OF THE FUNCTION F AT X             **
!     **  BY INTERPOLATION BY A THIRD ORDER POLYNOM                   **
!     **  F(X=I)=F(I)                                                 **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F(NR)
      REAL(8)   ,INTENT(IN) :: X
      REAL(8)   ,INTENT(OUT):: RES
      INTEGER(4)            :: INCR
      REAL(8)               :: XX1,XX2,XX3,XX4
      REAL(8)               :: P1,P2,P3,P4,P21,P32,P43,P321,P432,P4321
!     ******************************************************************
CALL ERROR$MSG('ROUTINE MARKED FOR DELETION')
CALL ERROR$STOP('BESSOV')
      INCR=INT(X)
      INCR=MIN0(INCR,NR-2)
      INCR=MAX0(INCR,2)
      XX1=X-DBLE(INCR-1)
      XX2=XX1-1.D0
      XX3=XX2-1.D0
      XX4=XX3-1.D0
      P1=F(INCR-1)
      P2=F(INCR)
      P3=F(INCR+1)
      P4=F(INCR+2)
      P21=-XX2*P1+XX1*P2
      P32=-XX3*P2+XX2*P3
      P43=-XX4*P3+XX3*P4
      P321=0.5D0*( -XX3*P21+XX1*P32 )
      P432=0.5D0*( -XX4*P32+XX2*P43 )
      P4321=1.D0/3.D0 * ( -XX4*P321+XX1*P432 )
      RES=P4321
      RETURN
      END
!
!     .....................................................BESSOV.......
      SUBROUTINE BESSOV(R1,DEX,NR,F,L,G,RMAX,RES)
!     **                                                              **
!     **  CALCULATES THE OVERLAP BETWEEN A SPHERICAL BESSEL FUNCTION  **
!     **  WITH THE FUNCTION F GIVEN ON AN LOGARITHMIC GRID            **
!     **  (3-D OVERLAPP)                                              **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1    ! 1. RADIAL GRID POINT
      REAL(8)   ,INTENT(IN) :: DEX   ! FACTOR BETWEEN SUBSEQUENT GRID POINTS
      INTEGER(4),INTENT(IN) :: NR    ! #(GRID POINTS)
      REAL(8)   ,INTENT(IN) :: F(NR) 
      INTEGER(4),INTENT(IN) :: L     ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: G
      REAL(8)   ,INTENT(IN) :: RMAX
      REAL(8)   ,INTENT(OUT):: RES
      REAL(8)               :: FAC(4)
      INTEGER(4)            :: NSTEP
      INTEGER(4)            :: I,IR
      REAL(8)               :: RSTEP
      REAL(8)               :: X,R
      REAL(8)               :: VAL
      REAL(8),EXTERNAL      :: BESSL
!     ******************************************************************
CALL ERROR$MSG('ROUTINE MARKED FOR DELETION')
CALL ERROR$STOP('BESSOV')
      NSTEP=20.D0*G/6.D0
      NSTEP=MAX0(NSTEP,400)
      RSTEP=RMAX/NSTEP
      FAC(1)=RSTEP*17.D0/48.D0
      FAC(2)=RSTEP*59.D0/48.D0
      FAC(3)=RSTEP*43.D0/48.D0
      FAC(4)=RSTEP*49.D0/48.D0
      X=1.D0+DLOG(RMAX/R1) / DEX
      CALL GETF(NR,F,X,VAL)
      RES=FAC(1)*VAL*BESSL(L,RMAX*G)*RMAX**2
      DO I=2,4
        R=RSTEP*DBLE(I-1)
        X=1.D0+DLOG(R/R1) / DEX
        CALL GETF(NR,F,X,VAL)
        RES=RES+FAC(I)*VAL*BESSL(L,R*G)*R**2
        R=RMAX-R
        X=1.D0+DLOG(R/R1) / DEX
        CALL GETF(NR,F,X,VAL)
        RES=RES+FAC(I)*VAL*BESSL(L,R*G)*R**2
      ENDDO
      DO IR=5,NSTEP-3
        R=RSTEP*DBLE(IR-1)
        X=1.D0+DLOG(R/R1) / DEX
        CALL GETF(NR,F,X,VAL)
        RES=RES+RSTEP*VAL*BESSL(L,R*G)*R**2
      ENDDO
      RETURN
      END
!
!     ..............................................BESSL...............
      FUNCTION BESSL(L,X) RESULT(Y)
!     ******************************************************************
!     **                                                              **
!     **  CALCULATES THE SPHERICAL BESSEL FUNCTION                    **
!     **  ACCORDING TO ABRAMOWITZ AND STEGUN                          **
!     **    FORMULA 10.1.2              FOR   X < 8                   **
!     **    FORMULA 10.1.8 AND  10.1.9  FOR   X > 8                   **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
!     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER(4),INTENT(IN) :: L ! MAIN AGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: X ! ARGUMENT
      REAL(8)               :: Y ! BESSL FUNCTION AT X
      REAL(8)               :: TRIG(4)
      REAL(8)               :: FACUL(0:100)
      REAL(8)   ,PARAMETER  :: TOL=1.D-10
      REAL(8)               :: PI
      REAL(8)               :: ARG
      REAL(8)               :: XSQ
      REAL(8)               :: FAC
      INTEGER(4)            :: K,IL,II,ISVAR
!     ******************************************************************
CALL ERROR$MSG('ROUTINE MARKED FOR DELETION')
CALL ERROR$STOP('BESSOV')
      IF(X.GT.DBLE(L)) THEN
        PI=4.D0*DATAN(1.D0)
        ARG=X-0.5D0*DBLE(L)*PI
        TRIG(1)=DSIN(ARG)/X
        TRIG(2)=DCOS(ARG)/X
        TRIG(3)=-TRIG(1)
        TRIG(4)=-TRIG(2)
        Y=TRIG(1)
        IF(L.EQ.0) RETURN
        FACUL(0)=1.D0
        DO I=1,2*L
          FACUL(I)=FACUL(I-1)*DBLE(I)
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
        IF(DABS(FAC).LT.TOL) GOTO 9999
      ENDDO
      CALL ERROR$MSG('Y NOT CONVERGED')
      CALL ERROR$I4VAL('L',L)
      CALL ERROR$R8VAL('X',X)
      CALL ERROR$STOP('FUNCTION BESSL')
9999  CONTINUE
      RETURN
      END

