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
!     ....................................................................
      SUBROUTINE radial_INTERPOLATE(NP,RI,FI_,R0,F0)
!     **                                                              **
!     **  POLYNOMIAL EXTRAPOLATION OF ORDER NP FROM NP POINTS (RI,FI_)**
!     **  TO THE POINT (R0,F0)                                        **
!     **                                                              **
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
        SVAR=1.D0
        DO J=1,I-1
          SVAR=SVAR*(R0-RI(J))/(RI(I)-RI(J))
        ENDDO
        F0=F0+SVAR*FI(I)
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
      END SUBROUTINE radial_INTERPOLATE
!
!     .......................................................SPHLSD.....
      SUBROUTINE RADIAL$DERIVE(R1,DEX,NR,F,G)
!     **                                                              **
!     **  TAKES THE RADIAL DERIVATIVE OF A FUNCTION F GIVEN ON A      **
!     **  LOGARITHMIC GRID                                            **
!     **                                                              **
      IMPLICIT NONE
      REAL(8)    ,INTENT(IN) :: R1
      REAL(8)    ,INTENT(IN) :: DEX
      INTEGER(4) ,INTENT(IN) :: NR
      REAL(8)    ,INTENT(IN) :: F(NR)
      REAL(8)    ,INTENT(OUT):: G(NR)
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
      REAL(8)                :: XEXP
      REAL(8)                :: RI
      INTEGER(4)             :: IR
!     ******************************************************************
      XEXP=DEXP(DEX)
!
!     ==================================================================
!     == FORM DERIVATIVE ON THE EQUI SPACED X-GRID                    ==
!     ==================================================================
      G(1)=(C11*F(1)+C12*F(2)+C13*F(3)+C14*F(4))/6.D0
      G(2)=(C21*F(1)+C22*F(2)+C23*F(3)+C24*F(4)+C25*F(5))/12.D0
      DO IR=3,NR-2
        G(IR)=CI1*F(IR-2)+CI2*F(IR-1)-CI2*F(IR+1)-CI1*F(IR+2)
      ENDDO
      G(NR-1)=-(C25*F(NR-4)+C24*F(NR-3) &
     &         +C23*F(NR-2)+C22*F(NR-1)+C21*F(NR))/12.D0
      G(NR)=-(C14*F(NR-3)+C13*F(NR-2)+C12*F(NR-1)+C11*F(NR))/6.D0
!
!     ==================================================================
!     ==  NOW MULTIPLY WITH DX/DR                                     ==
!     ==================================================================
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        G(IR)=G(IR)/(DEX*RI)
      ENDDO
      RETURN
      END
!
!     .....................................................INTRAD.......
      SUBROUTINE RADIAL$INTEGRATE(R1,DEX,NR,F1,G)
!     **                                                              **
!     **  INTEGRATES THE RADIAL FUNCTION F(R) TO OBTAINE G(R)         **
!     **  GIVEN ON A LOGHARITHMIC MESH                                **
!     **      G(R1)= INTEGRAL OVER R2 FROM 0 TO R1 OF F(R2)           **
!     **                                                              **
!     ** INPUT :                                                      **
!     **   R1,DEX,NR    PARAMETERS FOR THE LOGHARITHMIC MESH          **
!     **                R(I)=R1*EXP(DEX*(I-1)) ; I=1,NR               **
!     **   F            FUNCTION TO BE INTEGRATED                     **
!     ** OUTPUT :                                                     **
!     **   G            INTEGRAL OF F (FROM 0)                        **
!     **  REMARKS :                                                   **
!     **  (SEE: NUMERICAL RECIPES EQ: 4.1.14;                         **
!     **  INTEGRATES POLYNOMIALS UP TO 3. ORDER EXACTLY)              **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F1(NR)
      REAL(8)   ,INTENT(OUT):: G(NR)
      REAL(8)   ,PARAMETER  :: C0=-31.D0/48.D0
      REAL(8)   ,PARAMETER  :: C1=+11.D0/48.D0
      REAL(8)   ,PARAMETER  :: C2= -5.D0/48.D0
      REAL(8)   ,PARAMETER  :: C3= +1.D0/48.D0
      REAL(8)               :: F(NR)
      INTEGER(4)            :: IR
      REAL(8)               :: R2,R3
      REAL(8)               :: S21,S31,SUM0,END1
      REAL(8)               :: A,B,C
      REAL(8)               :: RI,XEXP
!     ******************************************************************
      IF(NR.LT.10) THEN
        CALL ERROR$MSG('NUMBER OF MESHPOINTS SMALLLER THAN 10')
        CALL ERROR$STOP('RADIAL$INTEGRATE')
      END IF
      F(:)=F1(:)
      XEXP=DEXP(DEX)
!     ==================================================================
!     ==  INTEGRATE FROM ZERO TO THE FIRST GRID-POINT                 ==
!     ==  EXTRAPOLATION BY A POLYNOMIAL OF 2. ORDER (QUADRATIC)       ==
!     ==================================================================
      R2=R1*XEXP
      R3=R2*XEXP
      S21=(F(2)-F(1))/(R2-R1)
      S31=(F(3)-F(1))/(R3-R1)
      C=(S21-S31)/(R2-R3)
      B=S21-C*(R1+R2)
      A=F(1)-B*R1-C*R1**2
      SUM0 = ( A*R1 + .5D0*B*R1**2 + 1.D0/3.D0*C*R1**3 ) / DEX
!     ==================================================================
!     ==  TRANSFORMATION ONTO LINEAR GRID                             ==
!     ==================================================================
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        F(IR)=F(IR)*DEX*RI
      ENDDO
!     ==================================================================
!     ==  SUMMATION                                                   ==
!     ==================================================================
      END1=C0*F(1)+C1*F(2)+C2*F(3)+C3*F(4)
      G(1)=F(1)+SUM0+END1
      DO IR=2,NR
        G(IR)=G(IR-1)+F(IR)
      ENDDO
!     ==================================================================
!     ==  FIX ENDPOINT                                                ==
!     ==================================================================
      DO IR=4,NR
        G(IR)=G(IR)+C0*F(IR)+C1*F(IR-1)+C2*F(IR-2)+C3*F(IR-3)
      ENDDO
!     ==================================================================
!     ==  FIX FIRST SEVEN GRID POINTS                                 ==
!     ==================================================================
      IF(NR.GE.6) THEN
        DO IR=1,3
          G(IR)=G(IR)-F(IR) &
     &             -(C0*F(IR)+C1*F(IR+1)+C2*F(IR+2)+C3*F(IR+3))
        ENDDO
      ELSE
        DO IR=1,3
          G(IR)=G(IR)-0.5D0*F(IR)
        ENDDO
      END IF
      RETURN
      END
!
!     .....................................................INTRAD.......
      SUBROUTINE RADIAL$INTEGRAL(R1,DEX,NR,F1,RES)
!     ******************************************************************
!     **                                                              **
!     **  INTEGRATES THE RADIAL FUNCTION F(R) FROM 0 TO GRID POINT NR **
!     **  WHERE F IS GIVEN ON A LOGHARITHMIC MESH                     **
!     **      RES  = INTEGRAL OVER R2 FROM 0 TO R1 OF F(R)            **
!     **                                                              **
!     ** INPUT :                                                      **
!     **   R1,DEX,NR    PARAMETERS FOR THE LOGHARITHMIC MESH          **
!     **                R(I)=R1*EXP(DEX*(I-1)) ; I=1,NR               **
!     **   F            FUNCTION TO BE INTEGRATED                     **
!     ** OUTPUT :                                                     **
!     **   RES          INTEGRAL OF F FROM 0 TO R(NR)                 **
!     **  REMARKS :                                                   **
!     **  (SEE: NUMERICAL RECIPES EQ: 4.1.14;                         **
!     **  INTEGRATES POLYNOMIALS UP TO 3. ORDER EXACTLY)              **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F1(NR)
      REAL(8)   ,INTENT(OUT):: RES
      REAL(8)   ,PARAMETER  :: C0=-31.D0/48.D0
      REAL(8)   ,PARAMETER  :: C1=+11.D0/48.D0
      REAL(8)   ,PARAMETER  :: C2= -5.D0/48.D0
      REAL(8)   ,PARAMETER  :: C3= +1.D0/48.D0
      REAL(8)               :: F(NR)
      INTEGER(4)            :: IR
      REAL(8)               :: R2,R3
      REAL(8)               :: S21,S31,SUM0
      REAL(8)               :: A,B,C
      REAL(8)               :: RI,XEXP
!     ******************************************************************
      IF(NR.LT.10) THEN
        CALL ERROR$MSG('NUMBER OF MESHPOINTS SMALLLER THAN 10')
        CALL ERROR$STOP('INTGRL')
      END IF
      DO IR=1,NR
        F(IR)=F1(IR)
      ENDDO
      XEXP=DEXP(DEX)
!     ==================================================================
!     ==  INTEGRATE FROM ZERO TO THE FIRST GRID-POINT                 ==
!     ==  EXTRAPOLATION BY A POLYNOMIAL OF 2. ORDER (QUADRATIC)       ==
!     ==================================================================
      R2=R1*XEXP
      R3=R2*XEXP
      S21=(F(2)-F(1))/(R2-R1)
      S31=(F(3)-F(1))/(R3-R1)
      C=(S21-S31)/(R2-R3)     
      B=S21-C*(R1+R2)
      A=F(1)-B*R1-C*R1**2
      SUM0 = A*R1 + .5D0*B*R1**2 + 1.D0/3.D0*C*R1**3 
!     ==================================================================
!     ==  TRANSFORMATION ONTO LINEAR GRID                             ==
!     ==================================================================
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        F(IR)=F(IR)*DEX*RI
      ENDDO
!     ==================================================================
!     ==  SUMMATION                                                   ==
!     ==================================================================
      RES=SUM0+C0*F(1)+C1*F(2)+C2*F(3)+C3*F(4)
      DO IR=1,NR
        RES=RES+F(IR)
      ENDDO
      IF(NR.GE.4) THEN
!       ================================================================
!       ==  FIX ENDPOINT                                              ==
!       ================================================================
        RES=RES+C0*F(NR)+C1*F(NR-1)+C2*F(NR-2)+C3*F(NR-3)
      ELSE
!       ================================================================
!       ==  FIX FIRST SEVEN GRID POINTS                               ==
!       ================================================================
        RES=RES-0.5D0*F(NR)
      END IF
      RETURN
      END
!
!     .....................................................INTRAD.......
      SUBROUTINE RADIAL$INTEGRAL1(R1,DEX,NR,F1,RAD,VAL)
!     ******************************************************************
!     **                                                              **
!     **   SUM could be done faster                                   **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F1(NR)
      REAL(8)   ,INTENT(IN) :: RAD
      REAL(8)   ,INTENT(OUT):: VAL
      REAL(8)   ,PARAMETER  :: C0=-31.D0/48.D0
      REAL(8)   ,PARAMETER  :: C1=+11.D0/48.D0
      REAL(8)   ,PARAMETER  :: C2= -5.D0/48.D0
      REAL(8)   ,PARAMETER  :: C3= +1.D0/48.D0
      REAL(8)               :: F(NR),G(NR)
      INTEGER(4)            :: IR,iint
      REAL(8)               :: R2,R3
      REAL(8)               :: S21,S31,SUM0,END1
      REAL(8)               :: A,B,C,RINT(4),GINT(4)
      REAL(8)               :: RI,XEXP
!     ******************************************************************
      IF(NR.LT.10) THEN
        CALL ERROR$MSG('NUMBER OF MESHPOINTS SMALLLER THAN 10')
        CALL ERROR$STOP('RADIAL$INTEGRATE')
      END IF
      F(:)=F1(:)
      XEXP=DEXP(DEX)
!     ==================================================================
!     ==  INTEGRATE FROM ZERO TO THE FIRST GRID-POINT                 ==
!     ==  EXTRAPOLATION BY A POLYNOMIAL OF 2. ORDER (QUADRATIC)       ==
!     ==================================================================
      R2=R1*XEXP
      R3=R2*XEXP
      S21=(F(2)-F(1))/(R2-R1)
      S31=(F(3)-F(1))/(R3-R1)
      C=(S21-S31)/(R2-R3)
      B=S21-C*(R1+R2)
      A=F(1)-B*R1-C*R1**2
      SUM0 = ( A*R1 + .5D0*B*R1**2 + 1.D0/3.D0*C*R1**3 ) / DEX
!     ==================================================================
!     ==  TRANSFORMATION ONTO LINEAR GRID                             ==
!     ==================================================================
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        F(IR)=F(IR)*DEX*RI
      ENDDO
!     ==================================================================
!     ==  SUMMATION                                                   ==
!     ==================================================================
      END1=C0*F(1)+C1*F(2)+C2*F(3)+C3*F(4)
      G(1)=F(1)+SUM0+END1
      DO IR=2,NR
        G(IR)=G(IR-1)+F(IR)
      ENDDO
!     ==================================================================
!     ==  FIX ENDPOINT                                                ==
!     ==================================================================
      DO IR=4,NR
        G(IR)=G(IR)+C0*F(IR)+C1*F(IR-1)+C2*F(IR-2)+C3*F(IR-3)
      ENDDO
!     ==================================================================
!     ==  FIX FIRST SEVEN GRID POINTS                                 ==
!     ==================================================================
      IF(NR.GE.6) THEN
        DO IR=1,3
          G(IR)=G(IR)-F(IR) &
     &             -(C0*F(IR)+C1*F(IR+1)+C2*F(IR+2)+C3*F(IR+3))
        ENDDO
      ELSE
        DO IR=1,3
          G(IR)=G(IR)-0.5D0*F(IR)
        ENDDO
      END IF

      IR=MIN(INT(LOG(RAD/R1)/DEX),NR-3)
      DO iint=1,4
        RINT(iint)=r1*DEXP((IR-1+iint-1)*DEX)
        GINT(iint)=G(IR-1+iint)
      END DO
      CALL RADIAL_INTERPOLATE(4,RINT,GINT,RAD,VAL)      
      RETURN
      END
!
!     .....................................................INTRAD.......
      SUBROUTINE RADIAL$INTEGRAL2(R1,DEX,NR,F1,r,RES)
!     ******************************************************************
!     **                                                              **
!     **  INTEGRATES THE RADIAL FUNCTION F(R) FROM 0 TO GRID POINT NR **
!     **  WHERE F IS GIVEN ON A LOGHARITHMIC MESH                     **
!     **      RES  = INTEGRAL OVER R2 FROM 0 TO R1 OF F(R)            **
!     **                                                              **
!     ** INPUT :                                                      **
!     **   R1,DEX,NR    PARAMETERS FOR THE LOGHARITHMIC MESH          **
!     **                R(I)=R1*EXP(DEX*(I-1)) ; I=1,NR               **
!     **   F            FUNCTION TO BE INTEGRATED                     **
!     ** OUTPUT :                                                     **
!     **   RES          INTEGRAL OF F FROM 0 TO R(NR)                 **
!     **  REMARKS :                                                   **
!     **  (SEE: NUMERICAL RECIPES EQ: 4.1.14;                         **
!     **  INTEGRATES POLYNOMIALS UP TO 3. ORDER EXACTLY)              **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     **                                                              **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR
      REAL(8)   ,INTENT(IN) :: F1(NR)
      REAL(8)   ,INTENT(IN) :: r
      REAL(8)   ,INTENT(OUT):: RES
      REAL(8)   ,PARAMETER  :: C0=-31.D0/48.D0
      REAL(8)   ,PARAMETER  :: C1=+11.D0/48.D0
      REAL(8)   ,PARAMETER  :: C2= -5.D0/48.D0
      REAL(8)   ,PARAMETER  :: C3= +1.D0/48.D0
      REAL(8)               :: F(NR)
      INTEGER(4)            :: IR,ir1
      INTEGER(4)            :: i1   ! first grid point for interpolation
      REAL(8)               :: R2,R3
      REAL(8)               :: S21,S31,SUM0,sum
      REAL(8)               :: A,B,C
      REAL(8)               :: RI,XEXP
      real(8)               :: rp(4),val(4)
!     ******************************************************************
call error$MSG('ROUTINE MARKED FOR DELETION')
CALL ERROR$STOP('RADIAL$INTEGRAL2')
      IF(NR.LT.10) THEN
        CALL ERROR$MSG('NUMBER OF MESHPOINTS SMALLLER THAN 10')
        CALL ERROR$STOP('INTGRL')
      END IF
      DO IR=1,NR
        F(IR)=F1(IR)
      ENDDO
      XEXP=DEXP(DEX)
      i1=min(int(log(r/r1)/dex),1)
!      
!     ==================================================================
!     ==  INTEGRATE FROM ZERO TO THE FIRST GRID-POINT                 ==
!     ==  EXTRAPOLATION BY A POLYNOMIAL OF 2. ORDER (QUADRATIC)       ==
!     ==================================================================
      R2=R1*XEXP
      R3=R2*XEXP
      S21=(F(2)-F(1))/(R2-R1)
      S31=(F(3)-F(1))/(R3-R1)
      C=(S21-S31)/(R2-R3)     
      B=S21-C*(R1+R2)
      A=F(1)-B*R1-C*R1**2
      SUM0 = A*R1 + .5D0*B*R1**2 + 1.D0/3.D0*C*R1**3 
!     ==================================================================
!     ==  TRANSFORMATION ONTO LINEAR GRID                             ==
!     ==================================================================
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        F(IR)=F(IR)*DEX*RI
      ENDDO
!     ==================================================================
!     ==  SUMMATION                                                   ==
!     ==================================================================
      SUM=SUM0+C0*F(1)+C1*F(2)+C2*F(3)+C3*F(4)
      DO IR=1,I1-1
        SUM=SUM+F(IR)
      ENDDO
      IR1=0
      RP=R1*EXP(DEX*REAL(IR1-1,KIND=8))/XEXP
      val(:)=0.d0
      DO IR=I1,I1+3
        IR1=IR1+1
        RI=RI*XEXP
        RP(IR1)=RI
        SUM=SUM+F(IR)
        IF(IR.GE.4) THEN
!         ================================================================
!         ==  FIX ENDPOINT                                              ==
!         ================================================================
          VAL(IR1)=SUM+C0*F(IR)+C1*F(IR-1)+C2*F(IR-2)+C3*F(IR-3)
        ELSE
!         ================================================================
!         ==  FIX FIRST SEVEN GRID POINTS                               ==
!         ================================================================
          VAL(IR1)=VAL(IR1)-0.5D0*F(IR)
        END IF
      ENDDO
      CALL RADIAL_INTERPOLATE(4,RP,VAL,R,RES)
      RETURN
      END
!
!     ....................................................................
      SUBROUTINE RADIAL$VALUE(R1,DEX,NR,FI_,R0,F0)
!     **                                                                **
!     **  POLYNOMIAL EXTRAPOLATION OF ORDER NP FROM NP POINTS (RI,FI_)  **
!     **  TO THE POINT (R0,F0)                                          **
!     **                                                                **
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      INTEGER(4),PARAMETER  :: NP=4
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR 
      REAL(8)   ,INTENT(IN) :: FI_(NR)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      REAL(8)               :: RI(NP)
      REAL(8)               :: XEXP
      INTEGER(4)            :: IP,I
!     ******************************************************************
      XEXP=DEXP(DEX)
      IF(R0.GT.R1) THEN
        IP=INT(1.D0+DLOG(R0/R1)/DEX-0.5D0*DBLE(NP))+1
        IP=MAX(IP,1)
      ELSE
        IP=1
      END IF
      IP=MIN(IP,NR-NP+1)
      RI(1)=R1*DEXP(DEX*DBLE(IP-1))
      DO I=2,NP
        RI(I)=RI(I-1)*XEXP
      ENDDO
      CALL INTERPOLATE(NP,RI,FI_(IP:IP+NP-1),R0,F0)
      IF(TPR) THEN
        PRINT*,'===== TEST INTERP ======='
        DO I=1,NP
          PRINT*,'INTER',RI(I),FI_(I+IP-1)
        ENDDO
        PRINT*,'RES ',R0,F0
        PRINT*,'===== TEST INTERP FINISHED ======='
      END IF
      RETURN
      CONTAINS
!       ....................................................................
        SUBROUTINE INTERPOLATE(NP,RI,FI_,R0,F0)
!       **                                                              **
!       **  POLYNOMIAL EXTRAPOLATION OF ORDER NP FROM NP POINTS (RI,FI_)**
!       **  TO THE POINT (R0,F0)                                        **
!       **                                                              **
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: NP
        REAL(8)   ,INTENT(IN) :: RI(NP)
        REAL(8)   ,INTENT(IN) :: FI_(NP)
        REAL(8)   ,INTENT(IN) :: R0
        REAL(8)   ,INTENT(OUT):: F0
        REAL(8)               :: FI(NP)
        REAL(8)               :: SVAR
        INTEGER(4)            :: I,J,IP
!       ****************************************************************
        FI(:)=FI_(:)
        F0=0.D0
        DO I=1,NP
          SVAR=1.D0
          DO J=1,I-1
            SVAR=SVAR*(R0-RI(J))/(RI(I)-RI(J))
          ENDDO
          F0=F0+SVAR*FI(I)
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
        END SUBROUTINE INTERPOLATE
      END SUBROUTINE RADIAL$VALUE
!
!     ..................................................................
      SUBROUTINE RADIAL$DERIVATIVE(R1,DEX,NR,FI_,R0,F0)
!     ******************************************************************
!     **  DERIVATIVE F0 OF A RADIAL FUNCTION AT THE POINT R0          **
!     ****************************************************************
      IMPLICIT NONE
      LOGICAL(4),PARAMETER  :: TPR=.FALSE.
      INTEGER(4),PARAMETER  :: NP=4
      REAL(8)   ,INTENT(IN) :: R1
      REAL(8)   ,INTENT(IN) :: DEX
      INTEGER(4),INTENT(IN) :: NR 
      REAL(8)   ,INTENT(IN) :: FI_(NR)
      REAL(8)   ,INTENT(IN) :: R0
      REAL(8)   ,INTENT(OUT):: F0
      REAL(8)               :: RI(NP)
      REAL(8)               :: XEXP
      INTEGER(4)            :: IP,I
!     ******************************************************************
      XEXP=DEXP(DEX)
      IF(R0.GT.R1) THEN
        IP=INT(1.D0+DLOG(R0/R1)/DEX-0.5D0*DBLE(NP))+1
        IP=MAX(IP,1)
      ELSE
        IP=1
      END IF
      IP=MIN(IP,NR-NP+1)
      RI(1)=R1*DEXP(DEX*DBLE(IP-1))
      DO I=2,NP
        RI(I)=RI(I-1)*XEXP
      ENDDO
      CALL INTERPOLATEGRAD(NP,RI,FI_(IP:IP+NP-1),R0,F0)
      IF(TPR) THEN
        PRINT*,'===== TEST INTERP ======='
        DO I=1,NP
          PRINT*,'INTER',RI(I),FI_(I+IP-1)
        ENDDO
        PRINT*,'RES ',R0,F0
        PRINT*,'===== TEST INTERP FINISHED ======='
      END IF
      RETURN
      CONTAINS
!       ................................................................
        SUBROUTINE INTERPOLATEGRAD(NP,RI,FI_,R0,DF0)
!       ****************************************************************
!       **  GRADIENT OF A POLYNOMIAL EXTRAPOLATION OF ORDER NP        **
!       **  FROM NP POINTS (RI,FI_) TO POINT (R0,F0)                  **
!       ****************************************************************
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: NP
        REAL(8)   ,INTENT(IN) :: RI(NP)
        REAL(8)   ,INTENT(IN) :: FI_(NP)
        REAL(8)   ,INTENT(IN) :: R0
        REAL(8)   ,INTENT(OUT):: DF0
        REAL(8)               :: F0(NP)
        REAL(8)               :: FI(NP)
        REAL(8)               :: SVAR,dsvar,fac
        INTEGER(4)            :: I,J,IP
!       ****************************************************************
        FI(:)=FI_(:)
        F0=0.D0
        DF0=0.D0
        DO I=1,NP
          SVAR=1.D0
          DSVAR=0.D0
          DO J=1,I-1
!           == DO NOT CHANGE THE ORDER OF THE NEXT TWO STATEMENTS ======
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
        END SUBROUTINE INTERPOLATEgrad
      END SUBROUTINE RADIAL$DERIVATIVE
!
!     .....................................................POISON.......
      SUBROUTINE RADIAL$POISSON(R1,DEX,NR,L,RHO,V)
!     ******************************************************************
!     **                                                              **
!     ** SOLVES THE RADIAL POISSON EQUALTION FOR A GIVEN              **
!     ** ANGULAR MOMENTUM COMPONENT OF THE CHARGE DENSITY             **
!     **                                                              **
!     ** INPUT :                                                      **
!     **   R1,DEX,NR    PARAMETERS FOR THE LOGHARITHMIC MESH          **
!     **                R(I)=R1*EXP(DEX*(I-1)) ; I=1,NR               **
!     **   L            MAIN ANGULAR MOMENTUM QUANTUM NUMBER          **
!     **   RHO          INPUT CHARGE DENSITY                          **
!     ** OUTPUT :                                                      *
!     **   V            ELECTROSTATIC POTENTIAL                       **
!     ** MEMORY ARRAYS:                                               **
!     **   AUX1(NR),AUX2(NR)                                          **
!     ** SUBROUTINES USED:                                            **
!     **   INTRAD                                                     **
!     **                                                              **
!     **          P.E. BLOECHL, IBM RESEARCH LABORATORY ZURICH (1991) **
!     ******************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R1      ! FIRST POINT ON RADIAL GRID
      REAL(8)   ,INTENT(IN) :: DEX     ! LOGARITHMIC SPACING OF THE RADIAL GRID
      INTEGER(4),INTENT(IN) :: NR      ! NUMBER OF RADIAL GRID POINTS
      INTEGER(4),INTENT(IN) :: L       ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: RHO(NR) ! CHARGE DENSITY
      REAL(8)   ,INTENT(OUT):: V(NR)   ! ELECTROSTATIC POTENTIAL
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: AUX2(NR)
      REAL(8)               :: XEXP
      REAL(8)               :: PI
      REAL(8)               :: RI,FAC
      INTEGER(4)            :: IR
!     ******************************************************************
      XEXP=DEXP(DEX)
      PI=4.D0*DATAN(1.D0)
!
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        AUX1(IR)=RHO(IR)*RI**(L+2)
      ENDDO
!
      CALL RADIAL$INTEGRATE(R1,DEX,NR,AUX1,AUX2)
!
      FAC=4.D0*PI/DBLE(2*L+1)
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        V(IR)=FAC*RI**(-L-1)*AUX2(IR)
        AUX1(IR)=RHO(IR)*RI**(-L+1)
      ENDDO
!
      CALL RADIAL$INTEGRATE(R1,DEX,NR,AUX1,AUX2)
!
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        V(IR)=V(IR)+FAC*RI**L*(AUX2(NR)-AUX2(IR))
      ENDDO
!
      RETURN
      END
!
!     ...........................................MLTPOL.................
      SUBROUTINE RADIAL$MOMENT(R1,DEX,NR,L,RHO,QLM)
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
      REAL(8)   ,INTENT(IN) :: R1      ! FIRST POINT ON RADIAL GRID
      REAL(8)   ,INTENT(IN) :: DEX     ! LOGARITHMIC SPACING OF THE RADIAL GRID
      INTEGER(4),INTENT(IN) :: NR      ! NUMBER OF RADIAL GRID POINTS
      INTEGER(4),INTENT(IN) :: L       ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: RHO(NR) ! CHARGE DENSITY
      REAL(8)   ,INTENT(OUT):: QLM     ! MULTIPOLE MOMENT
      REAL(8)               :: AUX1(NR)
      REAL(8)               :: XEXP,RI
      INTEGER(4)            :: IR
!     ******************************************************************
      XEXP=DEXP(DEX)
      RI=R1/XEXP
      DO IR=1,NR
        RI=RI*XEXP
        AUX1(IR)=RI**(L+2)*RHO(IR)
      ENDDO
      CALL RADIAL$INTEGRAL(R1,DEX,NR,AUX1,QLM)
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
        G(2*I-1)=CL*REAL(XA(I)) 
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
        G(I)=REAL(XA(I))
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
      IMPLICIT none
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
      IMPLICIT none
      integer(4),intent(in) :: nr
      real(8)   ,intent(in) :: F(NR)
      real(8)   ,intent(in) :: x
      real(8)   ,intent(out):: res
      integer(4)            :: incr
      real(8)               :: xx1,xx2,xx3,xx4
      real(8)               :: p1,p2,p3,p4,p21,p32,p43,p321,p432,p4321
!     ******************************************************************
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
      IMPLICIT none
      real(8)   ,intent(in) :: r1    ! 1. radial grid point
      real(8)   ,intent(in) :: dex   ! factor between subsequent grid points
      integer(4),intent(in) :: nr    ! #(grid points)
      real(8)   ,intent(in) :: f(nr) 
      integer(4),intent(in) :: l     ! main angular momentum
      real(8)   ,intent(in) :: g
      real(8)   ,intent(in) :: rmax
      real(8)   ,intent(out):: res
      real(8)               :: fac(4)
      integer(4)            :: nstep
      integer(4)            :: i,ir
      real(8)               :: rstep
      real(8)               :: x,r
      real(8)               :: val
      real(8),external      :: bessl
!     ******************************************************************
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
      FUNCTION BESSL(L,X) RESult(Y)
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

