PROGRAM SKALA_TAU_DENSITY_TEST
  USE PAW_SKALA_PRIMITIVES_MODULE
  USE PAW_SKALA_GRID_MODULE
  USE MPE_MODULE, ONLY: MPE$COMBINE
  USE SETUP_MODULE, ONLY: SETUP_STATE=>THIS
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  REAL(8),PARAMETER :: PI=4.D0*ATAN(1.D0)
  COMPLEX(8),PARAMETER :: CI=(0.D0,1.D0)
  INTEGER :: NTASKS,THISTASK
  INTEGER,ALLOCATABLE :: COLOR(:)

  CALL MPE$INIT
  CALL MPE$QUERY('~',NTASKS,THISTASK)
  ALLOCATE(COLOR(NTASKS))
  COLOR=1
  CALL MPE$NEW('~','MONOMER',NTASKS,COLOR)
  CALL TIMING$START
  CALL TRACE$SETL4('ON',.FALSE.)
  CALL PLANE_WAVES()
  CALL CORE_ORBITALS()
  CALL FIXED_ORBITALS()
  IF(THISTASK.EQ.1) WRITE(*,'(A)') 'INDEPENDENT SKALA TAU CHECKS PASSED'
  CALL MPE$EXIT
CONTAINS
  SUBROUTINE CHECK(LABEL,ERROR,TOLERANCE)
    CHARACTER(*),INTENT(IN) :: LABEL
    REAL(8),INTENT(IN) :: ERROR,TOLERANCE
    REAL(8) :: TOTAL
    IF(.NOT.IEEE_IS_FINITE(ERROR)) ERROR STOP 'Nonfinite tau error'
    TOTAL=ERROR
    CALL MPE$COMBINE('~','MAX',TOTAL)
    IF(THISTASK.EQ.1) WRITE(*,'(A,1X,ES16.8)') LABEL,TOTAL
    IF(TOTAL.GT.TOLERANCE) ERROR STOP 'Independent tau reference failed'
  END SUBROUTINE

  SUBROUTINE PLANE_WAVES()
    INTEGER,PARAMETER :: NG=24,NB=3
    REAL(8) :: CELL(3,3),AINV(3,3),GBAS(3,3),VOLUME,K(3),Q(3),FRAC(3),X(3)
    REAL(8) :: F(NB),RREF,TREF,ERROR_RHO,ERROR_TAU,TOTAL,EXACT,SCALE
    REAL(8),ALLOCATABLE :: RHO(:,:),TAU(:,:),KIN(:,:)
    COMPLEX(8),ALLOCATABLE :: COEFF(:,:,:)
    COMPLEX(8) :: PSI,GRAD(3),PHASE,C
    INTEGER,ALLOCATABLE :: IG(:,:)
    INTEGER :: MODE,NGL,NRL,NR1L,FIRST,I,J,L,IR,G,B,NH,M(3),COUNT_MODE
    LOGICAL :: PACKED
    CHARACTER(32) :: ID
    CELL=RESHAPE([6.D0,0.D0,0.D0,0.4D0,7.D0,0.D0,-0.2D0,0.3D0,8.D0],[3,3])
    VOLUME=6.D0*7.D0*8.D0
    CALL PAW_SKALA_INVERSE3(CELL,AINV)
    GBAS=2.D0*PI*TRANSPOSE(AINV)
    IF(MOD(NG,NTASKS).NE.0) ERROR STOP 'Use a rank count dividing 24'
    NR1L=NG/NTASKS
    FIRST=(THISTASK-1)*NR1L+1
    NRL=NR1L*NG*NG
    ALLOCATE(RHO(NRL,1),TAU(NRL,1),KIN(NRL,1))
    DO MODE=1,4
      PACKED=MODE.GE.3
      K=0.D0
      IF(MODE.EQ.2) K=[0.17D0,-0.11D0,0.23D0]
      IF(MODE.EQ.4) K=[0.5D0,0.D0,0.D0]
      WRITE(ID,'(A,I0)') 'tau-reference-',MODE
      CALL PLANEWAVE$INITIALIZE(ID,'~',CELL,K,PACKED,2.D0,FIRST,NR1L,NG,NG)
      CALL PLANEWAVE$GETI4('NGL',NGL)
      NH=NB
      IF(PACKED) NH=(NB+1)/2
      ALLOCATE(IG(3,NGL),COEFF(NGL,1,NH))
      CALL PLANEWAVE$GETI4A('IGVEC',3*NGL,IG)
      COEFF=0.D0
      COUNT_MODE=0
      DO G=1,NGL
        DO B=1,NH
          C=FOURIER_COEFFICIENT(IG(:,G),B,MODE)
          COEFF(G,1,B)=C/SQRT(VOLUME)
          IF(ABS(C).GT.0.D0) COUNT_MODE=COUNT_MODE+1
        END DO
      END DO
      CALL MPE$COMBINE('~','+',COUNT_MODE)
      IF(COUNT_MODE.EQ.0) ERROR STOP 'Empty analytic plane-wave fixture'
      F=[0.81D0,0.23D0,0.47D0]
      CALL WAVES_DENSITY(NGL,NRL,1,NB,NH,F,COEFF,RHO,.TRUE.,KIN,.TRUE.,TAU)
      IF(ANY(.NOT.IEEE_IS_FINITE(RHO)).OR.ANY(.NOT.IEEE_IS_FINITE(TAU))) ERROR STOP 'Nonfinite plane-wave fields'
      IF(MINVAL(TAU).LT.-1.D-14) ERROR STOP 'Negative plane-wave positive tau'
      ERROR_RHO=0.D0
      ERROR_TAU=0.D0
      IR=0
      DO L=0,NG-1
        DO J=0,NG-1
          DO I=FIRST-1,FIRST+NR1L-2
            IR=IR+1
            FRAC=REAL([I,J,L],8)/REAL(NG,8)
            X=MATMUL(CELL,FRAC)
            RREF=0.D0
            TREF=0.D0
            DO B=1,NB
              IF(PACKED) THEN
                CALL REAL_ORBITAL(X,GBAS,MODE,B,PSI,GRAD)
              ELSE
                PSI=0.D0
                GRAD=0.D0
                DO G=1,3
                  M=0
                  IF(G.EQ.1) M=[1,0,0]
                  IF(G.EQ.2) M=[0,1,0]
                  IF(G.EQ.3) M=[1,1,0]
                  Q=MATMUL(GBAS,REAL(M,8)+K)
                  PHASE=EXP(CI*DOT_PRODUCT(Q,X))
                  C=FOURIER_COEFFICIENT(M,B,MODE)
                  PSI=PSI+C*PHASE
                  GRAD=GRAD+CI*Q*C*PHASE
                END DO
              END IF
              RREF=RREF+F(B)*ABS(PSI)**2/VOLUME
              TREF=TREF+0.5D0*F(B)*SUM(ABS(GRAD)**2)/VOLUME
            END DO
            ERROR_RHO=MAX(ERROR_RHO,ABS(RHO(IR,1)-RREF))
            ERROR_TAU=MAX(ERROR_TAU,ABS(TAU(IR,1)-TREF))
          END DO
        END DO
      END DO
      CALL CHECK(TRIM(ID)//' pointwise density',ERROR_RHO,2.D-12)
      CALL CHECK(TRIM(ID)//' pointwise positive tau',ERROR_TAU,2.D-12)
      TOTAL=SUM(TAU)*VOLUME/REAL(NG**3,8)
      CALL MPE$COMBINE('~','+',TOTAL)
      EXACT=0.D0
      IF(PACKED) THEN
        DO B=1,NB
          Q=MATMUL(GBAS,[1.D0,0.D0,0.D0])
          IF(MODE.EQ.4) Q=MATMUL(GBAS,[0.5D0,0.D0,0.D0])
          SCALE=0.7D0
          IF(B.EQ.2) SCALE=0.9D0
          IF(B.EQ.3) THEN
            SCALE=0.6D0
            Q=MATMUL(GBAS,[0.D0,1.D0,0.D0])
            IF(MODE.EQ.4) Q=MATMUL(GBAS,[0.5D0,1.D0,0.D0])
          END IF
          EXACT=EXACT+0.25D0*F(B)*SCALE**2*SUM(Q**2)
        END DO
      ELSE
        DO B=1,NB
          DO G=1,3
            M=0
            IF(G.EQ.1) M=[1,0,0]
            IF(G.EQ.2) M=[0,1,0]
            IF(G.EQ.3) M=[1,1,0]
            Q=MATMUL(GBAS,REAL(M,8)+K)
            EXACT=EXACT+0.5D0*F(B)*ABS(FOURIER_COEFFICIENT(M,B,MODE))**2*SUM(Q**2)
          END DO
        END DO
      END IF
      CALL CHECK(TRIM(ID)//' kinetic integral',ABS(TOTAL-EXACT),5.D-12)
      DEALLOCATE(IG,COEFF)
    END DO
  END SUBROUTINE

  COMPLEX(8) FUNCTION FOURIER_COEFFICIENT(G,B,MODE) RESULT(C)
    INTEGER,INTENT(IN) :: G(3),B,MODE
    INTEGER :: M(3),PARTNER(3)
    C=0.D0
    IF(MODE.LE.2) THEN
      IF(ALL(G.EQ.[1,0,0])) C=CMPLX(0.2D0*B,0.11D0*(B-2),8)
      IF(ALL(G.EQ.[0,1,0])) C=CMPLX(-0.13D0*(B+1),0.09D0*B,8)
      IF(ALL(G.EQ.[1,1,0])) C=CMPLX(0.07D0*B,-0.17D0,8)
    ELSE
      M=[1,0,0]
      IF(MODE.EQ.4) M=[0,0,0]
      PARTNER=-M
      IF(MODE.EQ.4) PARTNER=PARTNER-[1,0,0]
      IF(B.EQ.1) THEN
        IF(ALL(G.EQ.M)) C=0.5D0*(0.7D0+0.9D0)
        IF(ALL(G.EQ.PARTNER)) C=0.5D0*(0.7D0-0.9D0)
        IF(MODE.EQ.3.AND.ALL(G.EQ.0)) C=CI*0.11D0
      ELSE
        M=[0,1,0]
        PARTNER=-M
        IF(MODE.EQ.4) PARTNER=PARTNER-[1,0,0]
        IF(ALL(G.EQ.M)) C=0.3D0*EXP(CI*0.37D0)
        IF(ALL(G.EQ.PARTNER)) C=0.3D0*EXP(-CI*0.37D0)
      END IF
    END IF
  END FUNCTION

  SUBROUTINE REAL_ORBITAL(X,GBAS,MODE,B,PSI,GRAD)
    REAL(8),INTENT(IN) :: X(3),GBAS(3,3)
    INTEGER,INTENT(IN) :: MODE,B
    COMPLEX(8),INTENT(OUT) :: PSI,GRAD(3)
    REAL(8) :: Q(3),T
    Q=MATMUL(GBAS,[1.D0,0.D0,0.D0])
    IF(MODE.EQ.4) Q=MATMUL(GBAS,[0.5D0,0.D0,0.D0])
    T=DOT_PRODUCT(Q,X)
    IF(B.EQ.1) THEN
      PSI=0.7D0*COS(T)
      GRAD=-0.7D0*SIN(T)*Q
    ELSE IF(B.EQ.2) THEN
      PSI=0.9D0*SIN(T)
      IF(MODE.EQ.3) PSI=PSI+0.11D0
      GRAD=0.9D0*COS(T)*Q
    ELSE
      Q=MATMUL(GBAS,[0.D0,1.D0,0.D0])
      IF(MODE.EQ.4) Q=MATMUL(GBAS,[0.5D0,1.D0,0.D0])
      T=DOT_PRODUCT(Q,X)+0.37D0
      PSI=0.6D0*COS(T)
      GRAD=-0.6D0*SIN(T)*Q
    END IF
  END SUBROUTINE

  SUBROUTINE CORE_ORBITALS()
    INTEGER,PARAMETER :: NR=1001,NS=3
    REAL(8) :: R(NR),PHI(NR,NS),D(NR,NS),SMALL(NR,NS),DS(NR,NS)
    REAL(8) :: TAU(NR),REF(NR),DENSITY(NR),Z,C1,C2,C3,ERROR,WORK(NR),INTEGRAL(NR)
    INTEGER :: I,MODE
    ALLOCATE(SETUP_STATE)
    CALL RADIAL$NEW('LOG',SETUP_STATE%GID)
    CALL RADIAL$SETI4(SETUP_STATE%GID,'NR',NR)
    CALL RADIAL$SETR8(SETUP_STATE%GID,'R1',1.D-5)
    CALL RADIAL$SETR8(SETUP_STATE%GID,'DEX',LOG(4.D6)/REAL(NR-1,8))
    CALL RADIAL$R(SETUP_STATE%GID,NR,R)
    SETUP_STATE%ATOM%NC=NS
    ALLOCATE(SETUP_STATE%ATOM%LOFI(NS),SETUP_STATE%ATOM%SOFI(NS),SETUP_STATE%ATOM%FOFI(NS))
    ALLOCATE(SETUP_STATE%ATOM%AEPSI(NR,NS),SETUP_STATE%ATOM%AEPSISM(NR,NS))
    SETUP_STATE%ATOM%LOFI=[0,0,1]
    SETUP_STATE%ATOM%SOFI=[0,0,-1]
    Z=1.4D0
    C1=2.D0*Z**1.5D0
    C2=Z**1.5D0/(2.D0*SQRT(2.D0))
    C3=Z**2.5D0/(2.D0*SQRT(6.D0))
    PHI(:,1)=C1*EXP(-Z*R)
    D(:,1)=-Z*PHI(:,1)
    PHI(:,2)=C2*(2.D0-Z*R)*EXP(-0.5D0*Z*R)
    D(:,2)=C2*Z*(0.5D0*Z*R-2.D0)*EXP(-0.5D0*Z*R)
    PHI(:,3)=C3*R*EXP(-0.5D0*Z*R)
    D(:,3)=C3*(1.D0-0.5D0*Z*R)*EXP(-0.5D0*Z*R)
    DO MODE=1,2
      SMALL=0.D0
      DS=0.D0
      SETUP_STATE%ATOM%FOFI=[2.D0,0.D0,0.D0]
      IF(MODE.EQ.2) THEN
        SETUP_STATE%ATOM%FOFI=[2.D0,0.7D0,1.3D0]
        SMALL(:,1)=0.03D0*PHI(:,3)
        DS(:,1)=0.03D0*D(:,3)
        SMALL(:,3)=0.02D0*PHI(:,1)
        DS(:,3)=0.02D0*D(:,1)
      END IF
      SETUP_STATE%ATOM%AEPSI=PHI
      SETUP_STATE%ATOM%AEPSISM=SMALL
      CALL SETUP_TAUCORE(NR,TAU)
      IF(ANY(.NOT.IEEE_IS_FINITE(TAU)).OR.ANY(TAU.LT.0.D0)) ERROR STOP 'Invalid positive core tau'
      REF=0.D0
      DO I=1,NS
        REF=REF+0.5D0*SETUP_STATE%ATOM%FOFI(I)*(D(:,I)**2+DS(:,I)**2)/(4.D0*PI)
      END DO
      REF=REF+SETUP_STATE%ATOM%FOFI(3)*PHI(:,3)**2/(4.D0*PI*R**2)
      REF=REF+SETUP_STATE%ATOM%FOFI(1)*SMALL(:,1)**2/(4.D0*PI*R**2)
      ERROR=MAXVAL(ABS(TAU/SQRT(4.D0*PI)-REF))/MAXVAL(REF)
      CALL CHECK('Analytic core orbital pointwise relative error',ERROR,2.D-6)
      IF(MODE.EQ.1) THEN
        DENSITY=2.D0*Z**3/PI*EXP(-2.D0*Z*R)
        CALL CHECK('Hydrogenic 1s tau = Z^2 density / 2', &
 &                 MAXVAL(ABS(TAU/SQRT(4.D0*PI)-0.5D0*Z**2*DENSITY))/MAXVAL(REF),2.D-6)
        WORK=SQRT(4.D0*PI)*TAU*R**2
        CALL RADIAL$INTEGRATE(SETUP_STATE%GID,NR,WORK,INTEGRAL)
        CALL CHECK('Hydrogenic 1s occupied kinetic integral',ABS(INTEGRAL(NR)-Z**2),1.D-7)
      END IF
    END DO
    DEALLOCATE(SETUP_STATE)
  END SUBROUTINE

  SUBROUTINE FIXED_ORBITALS()
    INTEGER,PARAMETER :: NR=1601,NG=24,NP=10
    TYPE(PAW_SKALA_SOURCE) :: S
    REAL(8) :: CELL(3,3),CENTER(3),X(3),POINT(3),AE(5,2),PS(5,2)
    REAL(8) :: VALUE(5,2),REF(5,2),CORE(5),RV(2),GV(3,2),TV(2),GOT(5,2),ERROR,POINT_ERROR,TAU_ERROR
    REAL(8) :: COORDS(3,NP),RHO(NG,NG,NG,2),GRHO(NG,NG,NG,3,2),TAU(NG,NG,NG,2)
    COMPLEX(8) :: C(9,2,2),DUP(9,9),DDN(9,9)
    INTEGER :: I,J,K,B,D,MODE,P,Q
    S%NR=NR
    S%LNX=3
    S%LMNX=9
    S%NDIMD=2
    S%RCUT=5.5D0
    ALLOCATE(S%R(NR),S%PHI(NR,3,2),S%CORE(NR,2),S%LOX(3),S%DM(9,9,2),S%DH(9,9,2))
    S%LOX=[0,1,2]
    S%R=[(0.005D0*REAL(I-1,8)+1.D-5,I=1,NR)]
    DO Q=1,2
      DO I=1,3
        S%PHI(:,I,Q)=AMPLITUDE(Q)*S%R**(I-1)*(1.D0-ALPHA(I,Q)*S%R**2)
      END DO
    END DO
    S%CORE(:,1)=SQRT(4.D0*PI)*2.D0/PI*EXP(-2.D0*S%R)
    S%CORE(:,2)=0.5D0*S%CORE(:,1)
    CALL PAW_SKALA_SOLID_BASIS_INIT(S%BASIS,2)
    DO D=1,2
      DO B=1,2
        DO I=1,9
          C(I,B,D)=CMPLX(0.13D0*SIN(REAL(I+2*B+D,8)),0.09D0*COS(REAL(2*I+B+D,8)),8)
        END DO
      END DO
    END DO
    COORDS(:,1)=[0.37D0,-0.29D0,0.43D0]
    COORDS(:,2)=[-0.51D0,0.27D0,0.11D0]
    COORDS(:,3)=[0.7D0,0.D0,0.D0]
    COORDS(:,4)=[0.D0,0.7D0,0.D0]
    COORDS(:,5)=[0.D0,0.D0,0.7D0]
    COORDS(:,6)=[-0.83D0,-0.41D0,0.31D0]
    COORDS(:,7)=[1.11D0,0.21D0,-0.37D0]
    COORDS(:,8)=[0.003D0,0.002D0,-0.001D0]
    COORDS(:,9)=0.D0
    COORDS(:,10)=[-0.13D0,0.91D0,0.27D0]
    CELL=0.D0
    DO I=1,3
      CELL(I,I)=12.D0
    END DO
    CENTER=6.D0
    DO MODE=1,2
      S%NDIMD=MODE
      DUP=0.D0
      DDN=0.D0
      DO B=1,2
        DO J=1,9
          DO I=1,9
            DUP(I,J)=DUP(I,J)+OCC(B,1)*CONJG(C(I,B,1))*C(J,B,1)
            DDN(I,J)=DDN(I,J)+OCC(B,2)*CONJG(C(I,B,2))*C(J,B,2)
          END DO
        END DO
      END DO
      IF(MODE.EQ.1) DDN=DUP
      S%DM(:,:,1)=DUP+DDN
      S%DM(:,:,2)=DUP-DDN
      DO K=1,NG
        DO J=1,NG
          DO I=1,NG
            X=12.D0*REAL([I-1,J-1,K-1],8)/NG-CENTER
            CALL ORBITAL_FIELDS(X,2,C,MODE,PS)
            RHO(I,J,K,1)=PS(1,1)+PS(1,2)
            RHO(I,J,K,2)=PS(1,1)-PS(1,2)
            GRHO(I,J,K,:,1)=PS(2:4,1)+PS(2:4,2)
            GRHO(I,J,K,:,2)=PS(2:4,1)-PS(2:4,2)
            TAU(I,J,K,1)=PS(5,1)+PS(5,2)
            TAU(I,J,K,2)=PS(5,1)-PS(5,2)
          END DO
        END DO
      END DO
      ERROR=0.D0
      POINT_ERROR=0.D0
      TAU_ERROR=0.D0
      DO P=1,NP
        X=COORDS(:,P)
        CALL ORBITAL_FIELDS(X,1,C,MODE,AE)
        CALL ORBITAL_FIELDS(X,2,C,MODE,PS)
        CALL CORE_FIELDS(X,CORE)
        REF=AE-PS
        REF(:,1)=REF(:,1)+0.5D0*CORE
        REF(:,2)=REF(:,2)+0.5D0*CORE
        CALL PAW_SKALA_SOURCE_POINT(S,X,VALUE)
        IF(ANY(.NOT.IEEE_IS_FINITE(VALUE))) ERROR STOP 'Nonfinite reconstructed partial-wave fields'
        POINT_ERROR=MAX(POINT_ERROR,MAXVAL(ABS(VALUE-REF)))
        POINT=CENTER+X
        CALL PAW_SKALA_GRID_INTERPOLATE_PRIMITIVES(NG,NG,NG,CELL,2,RHO,GRHO,TAU,POINT,RV,GV,TV)
        DO D=1,2
          GOT(1,D)=0.5D0*(RV(1)+REAL(3-2*D,8)*RV(2))
          GOT(2:4,D)=0.5D0*(GV(:,1)+REAL(3-2*D,8)*GV(:,2))
          GOT(5,D)=0.5D0*(TV(1)+REAL(3-2*D,8)*TV(2))
          GOT(:,D)=GOT(:,D)+VALUE(:,D)
          REF(:,D)=AE(:,D)+0.5D0*CORE
        END DO
        ERROR=MAX(ERROR,MAXVAL(ABS(GOT-REF)))
        IF(ANY(.NOT.IEEE_IS_FINITE(GOT))) ERROR STOP 'Nonfinite interpolated PAW fields'
        IF(ANY(GOT(5,:).LT.-1.D-14)) ERROR STOP 'Negative fixed-orbital positive tau'
        TAU_ERROR=MAX(TAU_ERROR,MAXVAL(ABS(GOT(5,:)-REF(5,:))))
      END DO
      CALL CHECK('Fixed complex orbitals AE-PS-core pointwise error',POINT_ERROR,2.D-9)
      CALL CHECK('Fixed orbitals interpolated PAW field error',ERROR,2.D-9)
      CALL CHECK('Fixed orbitals interpolated PAW tau error',TAU_ERROR,2.D-9)
    END DO
  END SUBROUTINE

  REAL(8) FUNCTION AMPLITUDE(Q) RESULT(A)
    INTEGER,INTENT(IN) :: Q
    A=1.D0
    IF(Q.EQ.2) A=0.8D0
  END FUNCTION

  REAL(8) FUNCTION ALPHA(L,Q) RESULT(A)
    INTEGER,INTENT(IN) :: L,Q
    A=0.025D0*REAL(L+Q,8)
  END FUNCTION

  REAL(8) FUNCTION OCC(B,D) RESULT(F)
    INTEGER,INTENT(IN) :: B,D
    F=0.9D0-0.17D0*B-0.13D0*D
  END FUNCTION

  SUBROUTINE CORE_FIELDS(X,V)
    REAL(8),INTENT(IN) :: X(3)
    REAL(8),INTENT(OUT) :: V(5)
    REAL(8) :: R
    R=SQRT(SUM(X**2))
    V(1)=2.D0/PI*EXP(-2.D0*R)
    V(2:4)=0.D0
    IF(R.GT.0.D0) V(2:4)=-2.D0*V(1)*X/R
    ! The cusp has no direction at r=0. Use the production limiting direction.
    IF(R.EQ.0.D0) V(2)=-2.D0*V(1)
    V(5)=0.5D0*V(1)
  END SUBROUTINE

  SUBROUTINE ORBITAL_FIELDS(X,Q,C,MODE,V)
    REAL(8),INTENT(IN) :: X(3)
    INTEGER,INTENT(IN) :: Q,MODE
    COMPLEX(8),INTENT(IN) :: C(9,2,2)
    REAL(8),INTENT(OUT) :: V(5,2)
    REAL(8) :: U(9),GU(3,9),W(9),GW(3,9),R2,A,N0,N1,N2,N3,F
    COMPLEX(8) :: PSI,GRAD(3)
    INTEGER :: I,L,D,B,SPIN
    ! Cartesian s/p/d polynomials are independent of the production harmonics,
    ! radial interpolation, density-matrix contraction and derivative routines.
    N0=1.D0/SQRT(4.D0*PI)
    N1=SQRT(3.D0/(4.D0*PI))
    N2=SQRT(15.D0/(16.D0*PI))
    N3=SQRT(5.D0/(16.D0*PI))
    U=[N0,N1*X(1),N1*X(3),N1*X(2),N2*(X(1)**2-X(2)**2), &
 &     2.D0*N2*X(1)*X(3),N3*(2.D0*X(3)**2-X(1)**2-X(2)**2), &
 &     2.D0*N2*X(2)*X(3),2.D0*N2*X(1)*X(2)]
    GU=0.D0
    GU(1,2)=N1
    GU(3,3)=N1
    GU(2,4)=N1
    GU(:,5)=2.D0*N2*[X(1),-X(2),0.D0]
    GU(:,6)=2.D0*N2*[X(3),0.D0,X(1)]
    GU(:,7)=2.D0*N3*[-X(1),-X(2),2.D0*X(3)]
    GU(:,8)=2.D0*N2*[0.D0,X(3),X(2)]
    GU(:,9)=2.D0*N2*[X(2),X(1),0.D0]
    R2=SUM(X**2)
    DO I=1,9
      L=INT(SQRT(REAL(I-1,8)))+1
      A=ALPHA(L,Q)
      W(I)=AMPLITUDE(Q)*(1.D0-A*R2)*U(I)
      GW(:,I)=AMPLITUDE(Q)*((1.D0-A*R2)*GU(:,I)-2.D0*A*U(I)*X)
    END DO
    V=0.D0
    DO D=1,2
      SPIN=D
      IF(MODE.EQ.1) SPIN=1
      DO B=1,2
        PSI=SUM(C(:,B,SPIN)*W)
        GRAD=MATMUL(GW,C(:,B,SPIN))
        F=OCC(B,SPIN)
        V(1,D)=V(1,D)+F*ABS(PSI)**2
        V(2:4,D)=V(2:4,D)+2.D0*F*REAL(CONJG(PSI)*GRAD,8)
        V(5,D)=V(5,D)+0.5D0*F*SUM(ABS(GRAD)**2)
      END DO
    END DO
  END SUBROUTINE
END PROGRAM
