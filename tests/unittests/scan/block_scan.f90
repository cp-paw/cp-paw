      SUBROUTINE TEST_SCAN_BLOCK(NFAIL)
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT NONE
      INTEGER,INTENT(INOUT) :: NFAIL
      INTEGER,PARAMETER :: NP=2
      REAL(8) :: RHO(NP,2),GRAD(NP,3,2),TAU(NP,2),WEIGHT(NP)
      REAL(8) :: DRHO(NP,2),DGRAD(NP,3,2),DTAU(NP,2),DWEIGHT(NP),ENERGY
      REAL(8) :: GR(3,2),V(9),D(9),E,REF,EPS,EP,EM,EXPECTED,ERROR,ACTUAL,BASE
      REAL(8) :: RP(NP,2),GP(NP,3,2),TP(NP,2),WP(NP)
      REAL(8) :: UNUSED_R(NP,2),UNUSED_G(NP,3,2),UNUSED_T(NP,2),UNUSED_W(NP),LIMITS(8)
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
      RHO(1,:)=[0.5D0,0.5D0]
      RHO(2,:)=[0.25D0,0.D0]
      TAU(1,:)=[0.25D0,0.25D0]
      TAU(2,:)=[-0.01D0,0.01D0]
      GRAD=0.D0
      GRAD(2,1,1)=0.2D0
      WEIGHT=[2.D0,3.D0]
      CALL SCAN$FIELD_LIMITS(NP,RHO,GRAD,TAU,WEIGHT,LIMITS)
      ERROR=MAXVAL(ABS(LIMITS-[4.D0,1.D0,1.D0,1.D0,0.75D0,2.75D0,0.09D0,0.03D0]))
      CALL TEST_REPORT('SCAN RAW-FIELD LIMIT DIAGNOSTICS',ERROR,0.D0,1.D-14,NFAIL)
      ! A large affected electron count need not imply a large deficit.
      TAU(1,:)=0.25D0-1.D-14
      TAU(2,:)=0.5D0
      GRAD=0.D0
      GRAD(1,1,:)=1.D0
      CALL SCAN$FIELD_LIMITS(NP,RHO,GRAD,TAU,WEIGHT,LIMITS)
      CALL TEST_REPORT('SCAN ROUNDING-LIMIT AFFECTED ELECTRONS',LIMITS(5),2.D0,1.D-14,NFAIL)
      CALL TEST_REPORT('SCAN ROUNDING-LIMIT INTEGRATED DEFICIT',LIMITS(7),4.D-14,2.D-15,NFAIL)
      CALL TEST_SCAN_LARGE_BLOCK(NFAIL)
      END SUBROUTINE TEST_SCAN_BLOCK

      SUBROUTINE TEST_SCAN_LARGE_BLOCK(NFAIL)
      IMPLICIT NONE
      INTEGER,INTENT(INOUT) :: NFAIL
      INTEGER,PARAMETER :: NP=521
      INTEGER :: I
      REAL(8) :: R(NP,2),G(NP,3,2),T(NP,2),W(NP),DR(NP,2),DG(NP,3,2),DT(NP,2),DW(NP)
      REAL(8) :: E,REF,ERR,EI,RR(1,2),RG(1,3,2),RT(1,2),RW(1)
      DO I=1,NP
        R(I,:)=[0.1D0+MOD(I,17)*0.01D0,0.04D0+MOD(I,13)*0.006D0]
        G(I,:,1)=[0.02D0,-0.01D0,0.04D0]*R(I,1)
        G(I,:,2)=[-0.01D0,0.03D0,0.01D0]*R(I,2)
        T(I,:)=R(I,:)*0.76D0
        W(I)=REAL(MOD(I,9),8)/REAL(NP,8)
      END DO
      CALL SCAN$BLOCK(NP,R,G,T,W,E,DR,DG,DT,DW)
      REF=0.D0
      ERR=0.D0
      DO I=1,NP
        CALL SCAN$BLOCK(1,R(I:I,:),G(I:I,:,:),T(I:I,:),W(I:I),EI,RR,RG,RT,RW)
        REF=REF+EI
        ERR=MAX(ERR,ABS(DW(I)-RW(1)),MAXVAL(ABS(RR(1,:)-DR(I,:))), &
     &          MAXVAL(ABS(RG(1,:,:)-DG(I,:,:))),MAXVAL(ABS(RT(1,:)-DT(I,:))))
      END DO
      CALL TEST_REPORT('SCAN LARGE BLOCK ENERGY',E,REF,2.D-14,NFAIL)
      CALL TEST_REPORT('SCAN LARGE BLOCK ROW RESPONSES',ERR,0.D0,2.D-14,NFAIL)
      END SUBROUTINE TEST_SCAN_LARGE_BLOCK
