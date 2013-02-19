!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE TEST_LMTO$STRUCTURECONSTANTS()
!     **************************************************************************
!     ** TESTS STRUCTURE CONSTANTS BY COMPARING THE SOLID HANKEL FUNCTION     **
!     ** AGAINST THE EXPANSION INTO SOLID BESSEL FUNCTIONS                    **
!     **                                                                      **
!     ** IF TEST FAILS, IT INTERRUPTS WITH AN ERROR MSG                       **
!     **************************************************************************
      USE STRINGS_MODULE
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: L1X=2
      INTEGER(4),PARAMETER :: L2X=6
      INTEGER(4),PARAMETER :: NP=11
      REAL(8)              :: CENTER(3)  ! POSITION OF EXPANSION CENTER
      REAL(8)              :: K2         ! KAPPA**2
      REAL(8)              :: S((L2X+1)**2,(L1X+1)**2)
      REAL(8)              :: H((L1X+1)**2)
      REAL(8)              :: J((L2X+1)**2)
      REAL(8)              :: MINUSJS((L1X+1)**2)
      INTEGER(4)           :: IP,IK
      INTEGER(4)           :: NFIL
      REAL(8)              :: R(3)
      REAL(8)              :: DR(3)
      REAL(8)              :: DIS
      CHARACTER(64)        :: FILE
      REAL(8)              :: MAXDEV
      LOGICAL    ,PARAMETER:: TPR=.FALSE.
!     **************************************************************************
      CENTER(:)=(/0.D0,1.D0,5.D0/)
      DR(:)=(/0.D0,1.D0,1.D0/)
      MAXDEV=0.D0
      DO IK=-1,1
        K2=REAL(IK,KIND=8)
        IF(TPR) THEN
          WRITE(FILE,*)IK
          FILE='TEST_STRUCTURECONSTANTS_K2EQUALS'//TRIM(ADJUSTL(FILE))//'.DAT'
          CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,-FILE)
          CALL FILEHANDLER$UNIT('HOOK',NFIL)
        END IF
        DO IP=1,NP
          R=CENTER+DR*(-0.5D0+REAL(IP-1,KIND=8)/REAL(NP-1,KIND=8))
          CALL LMTO$SOLIDHANKEL(R,1.D-3,K2,(L1X+1)**2,H)
          CALL LMTO$STRUCTURECONSTANTS(CENTER,K2,L1X,L2X,S)
          CALL LMTO$SOLIDBESSEL(R-CENTER,K2,(L2X+1)**2,J)
          MINUSJS=-MATMUL(J,S) 
          MAXDEV=MAX(MAXDEV,MAXVAL(ABS(H-MINUSJS)))
!
!         ==PRINT ==============================================================
          IF(TPR) THEN
            DIS=SQRT(SUM((R-CENTER)**2)/SUM(CENTER**2))
            IF(2*(IP-1).LT.NP-1)DIS=-DIS
            WRITE(NFIL,*)DIS,MINUSJS,H
          END IF
        ENDDO
        IF(TPR) THEN
          CALL FILEHANDLER$CLOSE('HOOK')
        END IF
      ENDDO
!
!     == INTERRUPT WHEN TEST FAILS =============================================
      IF(MAXDEV.GT.1.D-6) THEN
        CALL ERROR$MSG('TEST FAILED')
        CALL ERROR$R8VAL('MAXDEV',MAXDEV)
        CALL ERROR$STOP('TEST_LMTO$STRUCTURECONSTANTS')
      END IF
!
!     == POINT 'HOOK' TO DEFAULT, WHICH INDICATES AN ERROR =====================
      IF(TPR) THEN
       CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,-'.FORGOTTOASSIGNFILETOHOOKERROR')
      END IF 
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SOLIDBESSEL(R,K2,LMX,J)
!     **************************************************************************
!     ** CONSTRUCTS THE REGULAR SOLUTIONS OF THE HELMHOLTZ EQUATION           **
!     **                   (NABLA^2 + K2)*PSI(R)=0                            **
!     **                                                                      **
!     ** THE SOLUTION BEHAVES AT THE ORIGIN LIKE                              **
!     **    J_{L,M}(R)=1/(2*L+1)!! * |R|^L  Y_{L,M}(R)  +O(R^L+1)             **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R(3)
      REAL(8)   ,INTENT(IN) :: K2   ! SQUARE OF THE WAVE VECTOR
      INTEGER(4),INTENT(IN) :: LMX
      REAL(8)   ,INTENT(OUT):: J(LMX)
      INTEGER(4)            :: LX
      REAL(8)               :: K 
      REAL(8)               :: X,Y,DYDX
      INTEGER(4)            :: LM,M
      INTEGER(4)            :: L
      REAL(8)               :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      CALL SPHERICAL$YLM(LMX,R,J)
      LX=INT(SQRT(REAL(LMX+1.D-5)))-1
      K=SQRT(ABS(K2))     
      X=SQRT(SUM(R**2))
      LM=0
      DO L=0,LX
        IF(K2.GT.0.D0) THEN
          CALL SPFUNCTION$BESSEL(L,K*X,Y,DYDX)  ! ABRAMOWITZ 10.1.25
          Y=Y/K**L
        ELSE IF(K2.LT.0.D0) THEN
          CALL SPFUNCTION$MODBESSEL(L,K*X,Y,DYDX) !ABRAMOWITZ 10.2.4
          Y=Y/K**L
        ELSE
          CALL SPFUNCTION$BESSEL0(L,X,Y,DYDX)  ! ABRAMOWITZ 10.1.2
        END IF
        DO M=1,2*L+1
          LM=LM+1
          J(LM)=J(LM)*Y
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SOLIDHANKEL(R,RAD,K2,LMX,H)
!     **************************************************************************
!     ** CONSTRUCTS THE IRREGULAR SOLUTIONS OF THE HELMHOLTZ EQUATION         **
!     **                   (NABLA^2 + K2)*PSI(R)=0                            **
!     **                                                                      **
!     ** THE SOLUTION BEHAVES AT THE ORIGIN LIKE                              **
!     **    H_{L,M}(R)= (2*L-1)!! * |R|^{-L-1}  Y_{L,M}(R)                    **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R(3) ! POSITION RELATIVE TO THE ORIGIN
      REAL(8)   ,INTENT(IN) :: RAD  ! INSIDE RAD THE FUNCTION IS REGULARIZED
      REAL(8)   ,INTENT(IN) :: K2   ! SQUARE OF THE WAVE VECTOR
      INTEGER(4),INTENT(IN) :: LMX
      REAL(8)   ,INTENT(OUT):: H(LMX)  ! SOLID HANKEL FUNCTION
      INTEGER(4)            :: LX
      REAL(8)               :: K 
      REAL(8)               :: X,XR,Y,DYDX
      INTEGER(4)            :: LM,M
      INTEGER(4)            :: L
      REAL(8)               :: PI
      LOGICAL(4)            :: TCAP
      REAL(8)               :: A,B
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      CALL SPHERICAL$YLM(LMX,R,H)
      LX=INT(SQRT(REAL(LMX+1.D-5)))-1
      K=SQRT(ABS(K2))     
      XR=SQRT(SUM(R**2))
      TCAP=XR.LT.RAD
      X=MAX(XR,RAD)
      LM=0
      DO L=0,LX
        IF(K2.GT.0.D0) THEN
          CALL SPFUNCTION$NEUMANN(L,K*X,Y,DYDX)  ! ABRAMOWITZ 10.1.26
          Y=-Y*K**(L+1)
          DYDX=-DYDX*K**(L+2)
        ELSE IF(K2.LT.0.D0) THEN
          CALL SPFUNCTION$MODHANKEL(L,K*X,Y,DYDX) !ABRAMOWITZ 10.2.4
          Y=Y*2.D0/PI*K**(L+1)
          DYDX=DYDX*2.D0/PI*K**(L+2)
        ELSE
!         ==  Y(X)= 1/(2L-1)!! * X**(-L-1) 
          CALL SPFUNCTION$NEUMANN0(L,X,Y,DYDX)  ! ABRAMOWITZ 10.2.5
          Y=-Y     !
          DYDX=-DYDX 
        END IF 
!
!       == INSIDE RAD, MATCH A PARABOLA TIMES R**L =============================
        IF(TCAP) THEN
          B=0.5D0*(DYDX*X-REAL(L,KIND=8)*Y)/X**(L+2)
          A=Y/X**L-B*X**2
          IF(L.EQ.0) THEN
            Y=A+B*XR**2
            DYDX=2.D0*B*XR
          ELSE
            Y=A*XR**L+B*XR**(L+2)
            DYDX=REAL(L,KIND=8)*A*XR**(L-1)+REAL(L+2)*B*XR**(L+1)
          END IF
        END IF  
        DO M=1,2*L+1
          LM=LM+1
          H(LM)=H(LM)*Y
        ENDDO
      ENDDO
      RETURN
      END
!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$KBARMULTICENTER(N,NORB,QBAR,SBAR,C)
!     **************************************************************************
!     ** DETERMINES THE COEFFICIENTS FOR THE MULTICENTER EXPANSION OF THE     **
!     ** SCREENED HANKELFUNCTIONS IN TERMES OF UNSCREENED ONES.               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: NORB
      REAL(8)   ,INTENT(IN) :: QBAR(N)
      REAL(8)   ,INTENT(IN) :: SBAR(N,NORB)
      REAL(8)   ,INTENT(OUT):: C(N,NORB)
      INTEGER(4)            :: I
!     **************************************************************************
!
!     ==========================================================================
!     == SET UP COEFFICIENTS FOR THE MULTICENTER EXPANSION OF KBAR            ==
!     ==========================================================================
      DO I=1,NORB
        C(:,I)=QBAR(:)*SBAR(:,I)
        C(I,I)=C(I,I)+1.D0
      ENDDO
      RETURN
      END SUBROUTINE LMTO$KBARMULTICENTER

!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$EXPANDQBAR(NAT,LXX,LX,QBAR,N,QBARVEC)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      INTEGER(4),INTENT(IN) :: LXX
      INTEGER(4),INTENT(IN) :: LX(NAT)
      REAL(8)   ,INTENT(IN) :: QBAR(LXX+1,NAT)
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(OUT):: QBARVEC(N)
      INTEGER(4)            :: I,IAT,L,IM
!     **************************************************************************
      IF(N.NE.SUM((LX+1)**2)) THEN
        CALL ERROR$STOP('LMTO$EXPANDQBAR')
      END IF
      I=0
      DO IAT=1,NAT
        DO L=0,LX(IAT)
          DO IM=1,2*L+1
            I=I+1
            QBARVEC(I)=QBAR(L+1,IAT)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

!
!     ..1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SCREENEDSOLIDHANKEL(K2,NAT,RPOS,RAD,LX,N,NORB,C,R,F)
!     **************************************************************************
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: K2
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RPOS(3,NAT)
      REAL(8)   ,INTENT(IN) :: RAD(NAT)
      INTEGER(4),INTENT(IN) :: LX(NAT)
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: NORB
      REAL(8)   ,INTENT(IN) :: C(N,NORB)
      REAL(8)   ,INTENT(IN) :: R(3)
      REAL(8)   ,INTENT(OUT):: F(NORB)
      INTEGER(4)            :: LMXX,LMX
      REAL(8)   ,ALLOCATABLE:: H(:)
      INTEGER(4)            :: I,IAT
!     **************************************************************************
      LMXX=(MAXVAL(LX)+1)**2
      ALLOCATE(H(LMXX))
      F(:)=0.D0
      I=0
      DO IAT=1,NAT
        LMX=(LX(IAT)+1)**2
        CALL LMTO$SOLIDHANKEL(R(:)-RPOS(:,IAT),RAD(IAT),K2,LMX,H)
        F(:)=F(:)+MATMUL(H(:LMX),C(I+1:I+LMX,:))
        I=I+LMX
      ENDDO
      DEALLOCATE(H)
      RETURN
      END SUBROUTINE LMTO$SCREENEDSOLIDHANKEL
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SOLIDBESSELRAD(L,R,K2,J,JDER)
!     **                                                                      **
!     ** CONSTRUCTS THE RADIAL PART OF THE REGULAR SOLUTIONS                  **
!     **  OF THE HELMHOLTZ EQUATION  (NABLA^2 + K2)*PSI(R)=0                  **
!     **                                                                      **
!     ** THE SOLUTION BEHAVES AT THE ORIGIN LIKE                              **
!     **    J_{L,M}(R)=1/(2*L+1)!! * |R|^L  Y_{L,M}(R)  +O(R^L+1)             **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: R
      REAL(8)   ,INTENT(IN) :: K2   ! SQUARE OF THE WAVE VECTOR
      REAL(8)   ,INTENT(OUT):: J
      REAL(8)   ,INTENT(OUT):: JDER
      REAL(8)               :: K 
      REAL(8)               :: X,Y,DYDX
      REAL(8)               :: PI
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      K=SQRT(ABS(K2))     
      X=R
      IF(K2.GT.0.D0) THEN
        CALL SPFUNCTION$BESSEL(L,K*X,Y,DYDX)  ! ABRAMOWITZ 10.1.25
        Y=Y/K**L
        DYDX=DYDX/K**(L-1)
      ELSE IF(K2.LT.0.D0) THEN
        CALL SPFUNCTION$MODBESSEL(L,K*X,Y,DYDX) !ABRAMOWITZ 10.2.4
        Y=Y/K**L
        DYDX=DYDX/K**(L-1)
      ELSE
        CALL SPFUNCTION$BESSEL0(L,X,Y,DYDX)  ! ABRAMOWITZ 10.1.2
      END IF
      J=Y
      JDER=DYDX
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SOLIDHANKELRAD(L,R,K2,H,HDER)
!     **                                                                      **
!     ** CONSTRUCTS RADIAL PART OF THE THE IRREGULAR SOLUTIONS                **
!     ** OF THE HELMHOLTZ EQUATION (NABLA^2 + K2)*PSI(R)=0                    **
!     **                                                                      **
!     ** THE SOLUTION BEHAVES AT THE ORIGIN LIKE                              **
!     **    H_{L,M}(R)=1/(2*L-1)!! * |R|^{-L-1}  Y_{L,M}(R)                   **
!     **                                                                      **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L    ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: R    ! RADIUS
      REAL(8)   ,INTENT(IN) :: K2   ! SQUARE OF THE WAVE VECTOR
      REAL(8)   ,INTENT(OUT):: H    ! RADIAL PART OF THE HANKEL FUNCTION
      REAL(8)   ,INTENT(OUT):: HDER ! RADIAL DERIVATIVE OF THE HANKEL FUNCTION
      REAL(8)               :: K 
      REAL(8)               :: X,Y,DYDX
      REAL(8)               :: PI
      REAL(8)               :: SVAR
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      K=SQRT(ABS(K2))     
      X=R
      IF(K2.GT.0.D0) THEN
        CALL SPFUNCTION$NEUMANN(L,K*X,Y,DYDX)  ! ABRAMOWITZ 10.1.26
        SVAR=-K**(L+1)
        Y=SVAR*Y
        DYDX=SVAR*DYDX*K
      ELSE IF(K2.LT.0.D0) THEN
        CALL SPFUNCTION$MODHANKEL(L,K*X,Y,DYDX) !ABRAMOWITZ 10.2.4
        SVAR=2.D0/PI*K**(L+1)
        Y=SVAR*Y
        DYDX=SVAR*DYDX*K
      ELSE
!       ==  Y(X)= 1/(2L-1)!! * X**(-L-1) 
        CALL SPFUNCTION$NEUMANN0(L,X,Y,DYDX)  ! ABRAMOWITZ 10.2.5
        Y=-Y     ! 
        DYDX=-DYDX
      END IF   
      H=Y
      HDER=DYDX
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$STRUCTURECONSTANTS(R21,K2,L1X,L2X,S)
!     **************************************************************************
!     **  CONSTRUCTS THE STRUCTURE CONSTANTS THAT MEDIATE AN EXPANSION        **
!     **  OF A SOLID HANKEL FUNCTION H_{L,M}(R-R1) CENTERED AT R1             **
!     **  INTO SOLID BESSEL FUNCTIONS  J_{L,M}(R-R2) CENTERED AT R2           **
!     **                                                                      **
!     **    H_{L,M}(R-R1) = SUM_{L',M'} J_{L',M'}(R-R2) * S_{L',M',L,M}       **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      REAL(8)   ,INTENT(IN) :: R21(3) ! EXPANSION CENTER
      INTEGER(4),INTENT(IN) :: L1X
      INTEGER(4),INTENT(IN) :: L2X
      REAL(8)   ,INTENT(IN) :: K2 ! 2ME/HBAR**2
      REAL(8)   ,INTENT(OUT):: S((L2X+1)**2,(L1X+1)**2)
      REAL(8)   ,PARAMETER  :: RAD=1.D-6
      REAL(8)               :: K
      REAL(8)               :: PI
      REAL(8)               :: SVAR
      INTEGER(4)            :: L3X,LM1X,LM2X,LM3X
      INTEGER(4)            :: LM1,LM2,LM3,L,L1,L2,L3,IM,LM3A,LM3B
      INTEGER(4)            :: LOFLM((L1X+L2X+1)**2)
      REAL(8)               :: H((L1X+L2X+1)**2)
      REAL(8)               :: HANKEL ! HANKEL FUNCTION OF THE DISTANCE
      REAL(8)               :: CG     ! GAUNT COEFFICIENT
      COMPLEX(8)            :: KAPPA
!     **************************************************************************
      PI=4.D0*ATAN(1.D0)
      L3X=L1X+L2X
      LM1X=(L1X+1)**2
      LM2X=(L2X+1)**2
      LM3X=(L3X+1)**2
      LM3=0
      DO L=0,L3X
        DO IM=1,2*L+1
          LM3=LM3+1
          LOFLM(LM3)=L
        ENDDO
      ENDDO
      IF(K2.GT.0.D0) THEN
        KAPPA=CMPLX(0.D0,-SQRT(K2))
      ELSE IF(K2.LT.0.D0) THEN
        KAPPA=CMPLX(SQRT(-K2),0.D0)
      ELSE
        KAPPA=(0.D0,0.D0)
      END IF

!     == CALCULATE HANKEL FUNCTION OF THE DISTANCE =============================
      CALL LMTO$SOLIDHANKEL(R21,RAD,K2,LM3X,H)
!
!     ==========================================================================
      S(:,:)=0.D0
      DO LM1=1,LM1X
        L1=LOFLM(LM1)
        DO LM2=1,LM2X
          L2=LOFLM(LM2)
          LM3A=(ABS(L2-L1))**2+1
          LM3B=(L1+L2+1)**2
!          DO LM3=LM3A,LM3B
          DO LM3=1,LM3X
            L3=LOFLM(LM3)
            CALL SPHERICAL$GAUNT(LM1,LM2,LM3,CG)
            IF(CG.EQ.0.D0) CYCLE
            IF(K2.EQ.0.D0) THEN
              IF(L3.NE.L1+L2) CYCLE  ! AVOID 0**0
              SVAR=1.D0
            ELSE
              SVAR=REAL(KAPPA**(L1+L2-L3))
            END IF
            S(LM2,LM1)=S(LM2,LM1)+CG*H(LM3)*SVAR
          ENDDO
        ENDDO
      ENDDO
!
!     == MULTIPLY WITH -4*PI (-1)**L2 ==========================================
      S=-4.D0*PI*S
      DO LM2=1,LM2X
        S(LM2,:)=S(LM2,:)*(-1.D0)**LOFLM(LM2)
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$Q(L,RAD,VAL,DER,K2,QPAR)
!     **************************************************************************
!     ** PARAMETER NEEDED TO  SCREEN THE STRUCTURE CONSTANTS                  **
!     **   |K>-|J>QBAR HAS THE SAME LOGARITHMIC DERIVATIVE AS |PHIDOT>.       **
!     ** VAL AND DER ARE VALUE AND DERIVATIVE OF PHIDOT.                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: L
      REAL(8)   ,INTENT(IN) :: RAD
      REAL(8)   ,INTENT(IN) :: VAL
      REAL(8)   ,INTENT(IN) :: DER
      REAL(8)   ,INTENT(IN) :: K2
      REAL(8)   ,INTENT(OUT):: QPAR
      REAL(8)               :: JVAL,JDER
      REAL(8)               :: KVAL,KDER
!     **************************************************************************
      CALL LMTO$SOLIDBESSELRAD(L,RAD,K2,JVAL,JDER)
      CALL LMTO$SOLIDHANKELRAD(L,RAD,K2,KVAL,KDER)
      QPAR=(JVAL*DER-VAL*JDER)/(KVAL*DER-VAL*KDER)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$SCREEN(TSTART,N,NORB,QBAR,S,SBAR)
!     **************************************************************************
!     **  DETERMINES SCREENED STRUCTURE CONSTANTS FOR A CLUSTER               **
!     **      |KBAR_I>=SUM_J |K_J> (DELTA_JI+QBAR_J*SBAR_JI)                  **
!     **  START WITH SBAR=0 OR GIVE BETTER ESTIMATE                           **
!     **                                                                      **
!     **                                                                      **
!     **      SBAR=                                                           **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      LOGICAL(4),INTENT(IN) :: TSTART
      INTEGER(4),INTENT(IN) :: N            ! #(ORBITALS ON THE CLUSTER)
      INTEGER(4),INTENT(IN) :: NORB         ! #(ORBITALS ON CENTRAL SITE
      REAL(8)   ,INTENT(IN) :: QBAR(N)      !
      REAL(8)   ,INTENT(IN) :: S(N,N)       ! UNSCREENED-STRUCTURE CONSTANTS
      REAL(8)   ,INTENT(INOUT):: SBAR(N,NORB) !SCREENED STRUCTURE CONSTANTS
      REAL(8)   ,PARAMETER  :: TOL=1.D-5    ! TOLERANCE FOR CONVERENCE
      INTEGER(4),PARAMETER  :: NITER=1000   ! X#(ITERATIONS)
      REAL(8)               :: ALPHA=0.5D0  ! MIXING FACTOR
      REAL(8)               :: DSBAR(N,NORB)
      REAL(8)               :: S0(N,NORB)
      REAL(8)               :: A(N,N)
      REAL(8)               :: SBARBIG(N,N)
      INTEGER(4)            :: I
      REAL(8)               :: DELTA
      INTEGER(4)            :: ITER
      LOGICAL(4)            :: CONVG
!     **************************************************************************
                            CALL TRACE$PUSH('LMTO$SCREEN')
!
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
      IF(TSTART) THEN
        DO I=1,N
          A(:,I)=-QBAR(:)*S(:,I)
          A(I,I)=A(I,I)+1.D0
        ENDDO
        CALL LIB$MATRIXSOLVER8(N,N,N,TRANSPOSE(A),SBARBIG,TRANSPOSE(S))
        SBARBIG=TRANSPOSE(SBARBIG)
        SBAR(:,:)=SBARBIG(:,:NORB)
!
!     ==========================================================================
!     ==========================================================================
!     ==========================================================================
      ELSE
        S0(:,:)=S(:,:NORB)
        DO ITER=1,NITER
          DO I=1,NORB
            DSBAR(:,I)=QBAR(:)*SBAR(:,I)
            DSBAR(I,I)=1.D0+DSBAR(I,I)
          ENDDO
          DSBAR=MATMUL(S,DSBAR)-SBAR
          DELTA=MAXVAL(ABS(DSBAR))
          CONVG=DELTA.LT.TOL
          IF(CONVG) EXIT
          SBAR=SBAR+DSBAR*ALPHA
        ENDDO
        IF(.NOT.CONVG) THEN
          CALL ERROR$MSG('LOOP NOT CONVERGED')
          CALL ERROR$STOP('LMTO$SCREEN')
        END IF
      END IF
                            CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$CLUSTERSTRUCTURECONSTANTS(K2,NAT,RPOS,LX,QBAR,N,NORB,SBAR)
!     **************************************************************************
!     **  CONSTRUCTS THE STRUCTURE CONSTANTS THAT MEDIATE AN EXPANSION        **
!     **  OF A SOLID HANKEL FUNCTION H_{L,M}(R-R1) CENTERED AT R1             **
!     **                                                                      **
!     ** REMARK: THE CENTRAL ATOM IS THE FIRST ATOM IN THE LIST               **
!     **                                                                      **
!     ** REMARK: INITIALIZE SBAR WITH ZERO OR A BETTER ESTIMATE               **
!     **                                                                      **
!     **                                                                      **
!     **  CONSTRUCTS THE STRUCTURE CONSTANTS THAT MEDIATE AN EXPANSION        **
!     **  OF A SOLID HANKEL FUNCTION H_{L,M}(R-R1) CENTERED AT R1             **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NAT
      REAL(8)   ,INTENT(IN) :: RPOS(3,NAT)
      INTEGER(4),INTENT(IN) :: LX(NAT)
      REAL(8)   ,INTENT(IN) :: QBAR(N)
      REAL(8)   ,INTENT(IN) :: K2
      INTEGER(4),INTENT(IN) :: N
      INTEGER(4),INTENT(IN) :: NORB
      REAL(8)   ,INTENT(INOUT):: SBAR(N,NORB)
      INTEGER(4)            :: II,IAT,LN,L,IM
      INTEGER(4)            :: I1,I2
      INTEGER(4)            :: IAT1,IAT2
      REAL(8)               :: R1(3),R2(3)
      INTEGER(4)            :: LMN1,LMN2
      INTEGER(4)            :: L1X,L2X
      REAL(8)               :: S0(N,N)
      REAL(8)  ,ALLOCATABLE :: S1(:,:)
      REAL(8)               :: QBARVEC(N)
!     **************************************************************************
                               CALL TRACE$PUSH('LMTO$CLUSTERSTRUCTURECONSTANTS')
!
!     ==========================================================================
!     == CHECK CONSISTENCY OF ARRAY DIMENSIONS                                ==
!     ==========================================================================
      IF(SUM((LX+1)**2).NE.N) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
        CALL ERROR$STOP('LMTO$CLUSTERSTRUCTURECONSTANTS')
      END IF
      IF((LX(1)+1)**2.NE.NORB) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY DIMENSIONS')
        CALL ERROR$STOP('LMTO$CLUSTERSTRUCTURECONSTANTS')
      END IF
!
!     ==========================================================================
!     == SET UP BARE STRUCTURE CONSTANTS                                      ==
!     ==========================================================================
      S0(:,:)=0.D0
      I1=0
      DO IAT1=1,NAT
        R1(:)=RPOS(:,IAT1)
        L1X=LX(IAT1)
        LMN1=(L1X+1)**2
        LMN2=(MAXVAL(LX(:))+1)**2
        ALLOCATE(S1(LMN2,LMN1))
        I2=0
        DO IAT2=1,NAT
          IF(IAT2.EQ.IAT1) THEN
            I2=I2+LMN1
            CYCLE
          END IF
          R2(:)=RPOS(:,IAT2)
!WRITE(*,FMT='("R1=",3F10.4," DR=",3F10.4)')R1,R2-R1
          L2X=LX(IAT2)
          LMN2=(L2X+1)**2
          CALL LMTO$STRUCTURECONSTANTS(R2-R1,K2,L1X,L2X,S1(:LMN2,:))
          S0(I2+1:I2+LMN2,I1+1:I1+LMN1)=S1(:LMN2,:)
          I2=I2+(L2X+1)**2
        ENDDO
        DEALLOCATE(S1)
        I1=I1+(L1X+1)**2
      ENDDO
!
!     ==========================================================================
!     == SCREEN STRUCTURE CONSTANTS                                           ==
!     ==========================================================================
      SBAR=0.D0
!     == THE FIRST PARAMETER SWITCHES BETWEEN AN ITERATIVE AND A DIRECT SOLUTION
!     == OF THE EQUATION. THE ITERATIVE APPROACH MAY BE MORE EFFICIENT IN AN
!     == CAR-PARRINELLO LIKE APPROACH.
      CALL LMTO$SCREEN(.TRUE.,N,NORB,QBAR,S0,SBAR)
                               CALL TRACE$POP()
      RETURN
      END

!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE LMTO$NEIGHBORLIST(RBAS,NAT,R,RC,NNX,NNB,NNLIST)
!     **************************************************************************
!     **  THIS IS A SIMPLE NEIGHBORLIST ROUTINE                               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)   ,INTENT(IN) :: NAT       ! #(ATOMS)
      REAL(8)      ,INTENT(IN) :: RBAS(3,3) ! LATTICE VECTORS
      REAL(8)      ,INTENT(IN) :: R(3,NAT)  ! ATOM POSITIONS
      REAL(8)      ,INTENT(IN) :: RC        ! CUTOFF RADIUS FOR THE NEIGHORLIST
      INTEGER(4)   ,INTENT(IN) :: NNX       ! X#(NEIGHBORS PER ATOM)
      INTEGER(4)   ,INTENT(OUT):: NNB       ! #(NEIGHBORS)
      INTEGER(4)   ,INTENT(OUT):: NNLIST(5,NNX) ! NEIGHBORLIST (IAT1,IAT2,IT(3))
      REAL(8)                  :: RBASINV(3,3)
      REAL(8)                  :: RFOLD(3,NAT)
      REAL(8)                  :: X(3)      ! RELATIVE COORDINATES
      REAL(8)                  :: XFOLD(3)      ! RELATIVE COORDINATES
      INTEGER(4)               :: ITFOLD(3,NAT)   ! SHIFT
      INTEGER(4)               :: IAT,I,IAT1,IAT2,IT1,IT2,IT3
      REAL(8)                  :: RMAX2     ! SQUARED CUTOFF RADIUS
      REAL(8)                  :: TVEC(3)
      INTEGER(4)               :: ITVEC(3)
      INTEGER(4)               :: MIN1,MAX1,MIN2,MAX2,MIN3,MAX3
      REAL(8)                  :: D(3),D2
      REAL(8)                  :: X0,Y0,Z0
!     **************************************************************************
!
!     ==========================================================================
!     == FOLD ATOM POSITIONS INTO THE FIRST UNIT CELL                         ==
!     ==========================================================================
      CALL LIB$INVERTR8(3,RBAS,RBASINV)
      DO IAT=1,NAT
        X(:)=MATMUL(RBASINV,R(:,IAT))
        DO I=1,3
          XFOLD(I)=MODULO(X(I),1.D0)
        ENDDO
        ITFOLD(:,IAT)=NINT(XFOLD-X)
        RFOLD(:,IAT)=MATMUL(RBAS,XFOLD)
      ENDDO
!
!     ==========================================================================
!     == FOLD ATOM POSITIONS INTO THE FIRST UNIT CELL                         ==
!     ==========================================================================
      RMAX2=RC**2
      NNB=0
      DO IAT1=1,NAT
!       == PLACE ONSITE ELEMENT FOR EACH ATOM FIRST IN THE NEIGHBORLIST
        NNB=NNB+1
        IF(NNB.GT.NNX) THEN
          CALL ERROR$MSG('MAXIMUM NUMBER OF NEIGHBORS EXCEEDED')
          CALL ERROR$I4VAL('NNB',NNB)
          CALL ERROR$I4VAL('NNX',NNX)
          CALL ERROR$STOP('LMTO$NEIGHBORLIST')
        END IF
        NNLIST(1,NNB)=IAT1
        NNLIST(2,NNB)=IAT1
        NNLIST(3:5,NNB)=0
        DO IAT2=1,NAT
          X0=RFOLD(1,IAT1)-RFOLD(1,IAT2)
          Y0=RFOLD(2,IAT1)-RFOLD(2,IAT2)
          Z0=RFOLD(3,IAT1)-RFOLD(3,IAT2)
          CALL BOXSPH(RBAS,X0,Y0,Z0,RC,MIN1,MAX1,MIN2,MAX2,MIN3,MAX3)
!         == LOOP OVER BOXES IN THE NEIGHBORHOOD ===============================
          DO IT1=MIN1,MAX1
            DO IT2=MIN2,MAX2
              DO IT3=MIN3,MAX3
                IF(IAT1.EQ.IAT2.AND.IT1.EQ.0.AND.IT2.EQ.0.AND.IT3.EQ.0) CYCLE
                ITVEC(1)=IT1
                ITVEC(2)=IT2
                ITVEC(3)=IT3
                TVEC(:)=MATMUL(RBAS,REAL(ITVEC,KIND=8))
!               == DISTANCE CRITERION ==========================================
                D(:)=RFOLD(:,IAT2)+TVEC(:)-RFOLD(:,IAT1)
                D2=SUM(D(:)**2)
                IF(D2.GT.RMAX2) CYCLE
                NNB=NNB+1
                IF(NNB.GT.NNX) THEN
                  CALL ERROR$MSG('MAXIMUM NUMBER OF NEIGHBORS EXCEEDED')
                  CALL ERROR$I4VAL('NNB',NNB)
                  CALL ERROR$I4VAL('NNX',NNX)
                  CALL ERROR$STOP('LMTO$NEIGHBORLIST')
                END IF
                NNLIST(1,NNB)=IAT1
                NNLIST(2,NNB)=IAT2
                ITVEC(:)=ITVEC(:)+ITFOLD(:,IAT2)-ITFOLD(:,IAT1)
                NNLIST(3:5,NNB)=ITVEC(:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!!$DO IAT1=1,NAT
!!$ D(:)=R(:,IAT1)-RFOLD(:,IAT1)
!!$ WRITE(*,FMT='("D ",3F10.5,3I5)')D(:),ITFOLD(:,IAT1)
!!$ENDDO
!!$DO I=1,NNB
!!$  IAT1=NNLIST(1,I)
!!$  IAT2=NNLIST(2,I)
!!$  ITVEC(:)=NNLIST(3:5,I)
!!$  D(:)=R(:,IAT2)-R(:,IAT1)+MATMUL(RBAS,REAL(ITVEC,KIND=8))
!!$  WRITE(*,FMT='(I5," DIS ",F10.5," D ",3F10.5)')I,SQRT(SUM(D**2)),D(:)
!!$ENDDO
!!$PRINT*,'RC ',RC
!!$STOP
      RETURN
      END
