      SUBROUTINE TEST_SCAN_BLOCK(NFAIL)
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      INTEGER,INTENT(INOUT) :: NFAIL
      INTEGER,PARAMETER :: NP=2
      REAL(8) :: RHO(NP,2),GRAD(NP,3,2),TAU(NP,2),WEIGHT(NP)
      REAL(8) :: DRHO(NP,2),DGRAD(NP,3,2),DTAU(NP,2),DWEIGHT(NP),ENERGY
      REAL(8) :: GR(3,2),V(9),D(9),E,REF,EPS,EP,EM,EXPECTED,ERROR,ACTUAL,BASE
      REAL(8) :: RP(NP,2),GP(NP,3,2),TP(NP,2),WP(NP)
      REAL(8) :: UNUSED_R(NP,2),UNUSED_G(NP,3,2),UNUSED_T(NP,2),UNUSED_W(NP)
      INTEGER :: I,S,K,MODE
      RHO=RESHAPE([0.3D0,0.2D0,0.13D0,0.4D0],[NP,2])
      GRAD(:,:,1)=RESHAPE([0.14D0,-0.13D0,0.08D0,0.21D0,-0.17D0,0.03D0],[NP,3])
      GRAD(:,:,2)=RESHAPE([-0.07D0,0.18D0,0.05D0,-0.11D0,0.03D0,0.12D0],[NP,3])
      TAU=RESHAPE([0.8D0,0.7D0,0.4D0,1.2D0],[NP,2])
      WEIGHT=[0.7D0,1.3D0]
      CALL SCAN$BLOCK(NP,RHO,GRAD,TAU,WEIGHT,ENERGY,DRHO,DGRAD,DTAU,DWEIGHT)
      REF=0.D0
      DO I=1,NP
        GR(:,1)=GRAD(I,:,1)+GRAD(I,:,2)
        GR(:,2)=GRAD(I,:,1)-GRAD(I,:,2)
        V=[SUM(RHO(I,:)),RHO(I,1)-RHO(I,2),SUM(GR(:,1)**2),SUM(GR(:,2)**2) &
     &     ,SUM(GR(:,1)*GR(:,2)),0.D0,0.D0,SUM(TAU(I,:)),TAU(I,1)-TAU(I,2)]
        CALL SCAN$EVAL1(V,E,D)
        REF=REF+WEIGHT(I)*E
      END DO
      CALL TEST_REPORT('SCAN BLOCK WEIGHTED ENERGY',ENERGY,REF,2.D-14,NFAIL)
      ERROR=0.D0
      DO I=1,NP
        DO MODE=1,4
          DO S=1,MERGE(1,2,MODE.EQ.4)
            DO K=1,MERGE(3,1,MODE.EQ.2)
              RP=RHO
              GP=GRAD
              TP=TAU
              WP=WEIGHT
              SELECT CASE(MODE)
              CASE(1)
                BASE=RHO(I,S)
                EXPECTED=DRHO(I,S)
              CASE(2)
                BASE=GRAD(I,K,S)
                EXPECTED=DGRAD(I,K,S)
              CASE(3)
                BASE=TAU(I,S)
                EXPECTED=DTAU(I,S)
              CASE(4)
                BASE=WEIGHT(I)
                EXPECTED=DWEIGHT(I)
              END SELECT
              EPS=1.D-5*MAX(ABS(BASE),1.D-2)
              SELECT CASE(MODE)
              CASE(1)
                RP(I,S)=BASE+EPS
              CASE(2)
                GP(I,K,S)=BASE+EPS
              CASE(3)
                TP(I,S)=BASE+EPS
              CASE(4)
                WP(I)=BASE+EPS
              END SELECT
              CALL SCAN$BLOCK(NP,RP,GP,TP,WP,EP,UNUSED_R,UNUSED_G,UNUSED_T,UNUSED_W)
              SELECT CASE(MODE)
              CASE(1)
                RP(I,S)=BASE-EPS
              CASE(2)
                GP(I,K,S)=BASE-EPS
              CASE(3)
                TP(I,S)=BASE-EPS
              CASE(4)
                WP(I)=BASE-EPS
              END SELECT
              CALL SCAN$BLOCK(NP,RP,GP,TP,WP,EM,UNUSED_R,UNUSED_G,UNUSED_T,UNUSED_W)
              ACTUAL=(EP-EM)/(2.D0*EPS)
              IF(.NOT.IEEE_IS_FINITE(ACTUAL)) ERROR STOP 'Nonfinite SCAN block finite difference'
              ERROR=MAX(ERROR,ABS(ACTUAL-EXPECTED)/MAX(1.D0,ABS(EXPECTED)))
            END DO
          END DO
        END DO
      END DO
      CALL TEST_REPORT('SCAN BLOCK PRIMITIVE AND WEIGHT DERIVATIVES',ERROR,0.D0,2.D-8,NFAIL)
      END SUBROUTINE TEST_SCAN_BLOCK
