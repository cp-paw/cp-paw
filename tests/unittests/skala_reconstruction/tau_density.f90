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
  CALL MPE$NEW('~','K',NTASKS,COLOR)
  CALL TIMING$START
  CALL TRACE$SETL4('ON',.FALSE.)
  CALL PLANE_WAVES()
  CALL SPINOR_PLANE_WAVES()
  CALL CORE_ORBITALS()
  CALL FIXED_ORBITALS()
  CALL SPINOR_SOURCES()
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

  SUBROUTINE SPINOR_PLANE_WAVES()
    INTEGER,PARAMETER :: NG=24,NB=3
    COMPLEX(8) :: SIGMA(2,2,4),UNITARY(2,2),C(2),PSI(2),GRAD(3,2),PHASE
    COMPLEX(8) :: VQ(4),BFOURIER(2,2),HREF(2),MATRIX_ELEMENT
    COMPLEX(8),ALLOCATABLE :: COEFF(:,:,:),HPSI(:,:,:),HBASE(:,:,:)
    REAL(8) :: CELL(3,3),AINV(3,3),GBAS(3,3),K(3),Q(3),X(3),VOLUME,F(NB)
    REAL(8) :: AXIS(3),ROTATION(3,3),CROSS(3,3),ANGLE,RREF(4),TREF(4),INTEGRAL(4),EXACT(4)
    REAL(8) :: ERROR_RHO,ERROR_TAU,ERROR_ROTATION,BOUND
    REAL(8),PARAMETER :: VBASE(4)=[0.21D0,0.03D0,-0.04D0,0.07D0]
    REAL(8),PARAMETER :: VCOS(4)=[0.05D0,-0.02D0,0.03D0,0.01D0]
    REAL(8),PARAMETER :: VSIN(4)=[-0.02D0,0.04D0,0.01D0,-0.03D0]
    REAL(8) :: QOUT(3),EOP,EREF,STRESS(3,3),SREF(3,3),ERROR_H
    REAL(8),ALLOCATABLE :: RHO(:,:),TAU(:,:),KIN(:,:),RBASE(:,:),TBASE(:,:),VTAU(:,:)
    INTEGER,ALLOCATABLE :: IG(:,:)
    INTEGER :: NRL,NR1L,FIRST,NGL,MODE,ROT,I,J,L,IR,G,B,D,A,M(3),OUTG,DIFF(3)
    CHARACTER(32) :: ID
    SIGMA=0.D0
    SIGMA(1,1,1)=1.D0
    SIGMA(2,2,1)=1.D0
    SIGMA(1,2,2)=1.D0
    SIGMA(2,1,2)=1.D0
    SIGMA(1,2,3)=-CI
    SIGMA(2,1,3)=CI
    SIGMA(1,1,4)=1.D0
    SIGMA(2,2,4)=-1.D0
    AXIS=[1.D0,2.D0,-3.D0]/SQRT(14.D0)
    ANGLE=0.71D0
    UNITARY=COS(ANGLE/2.D0)*SIGMA(:,:,1)
    DO A=1,3
      UNITARY=UNITARY-CI*SIN(ANGLE/2.D0)*AXIS(A)*SIGMA(:,:,A+1)
    END DO
    CROSS=RESHAPE([0.D0,AXIS(3),-AXIS(2),-AXIS(3),0.D0,AXIS(1), &
 &                 AXIS(2),-AXIS(1),0.D0],[3,3])
    ROTATION=SIN(ANGLE)*CROSS
    DO J=1,3
      DO I=1,3
        ROTATION(I,J)=ROTATION(I,J)+(1.D0-COS(ANGLE))*AXIS(I)*AXIS(J)
      END DO
      ROTATION(J,J)=ROTATION(J,J)+COS(ANGLE)
    END DO
    CELL=RESHAPE([6.D0,0.D0,0.D0,0.4D0,7.D0,0.D0,-0.2D0,0.3D0,8.D0],[3,3])
    VOLUME=6.D0*7.D0*8.D0
    CALL PAW_SKALA_INVERSE3(CELL,AINV)
    GBAS=2.D0*PI*TRANSPOSE(AINV)
    NR1L=NG/NTASKS
    FIRST=(THISTASK-1)*NR1L+1
    NRL=NR1L*NG*NG
    ALLOCATE(RHO(NRL,4),TAU(NRL,4),KIN(NRL,4),RBASE(NRL,4),TBASE(NRL,4),VTAU(NRL,4))
    F=[0.81D0,0.23D0,0.47D0]
    DO MODE=1,2
      K=0.D0
      IF(MODE.EQ.2) K=[0.17D0,-0.11D0,0.23D0]
      WRITE(ID,'(A,I0)') 'tau-spinor-',MODE
      CALL PLANEWAVE$INITIALIZE(ID,'~',CELL,K,.FALSE.,2.D0,FIRST,NR1L,NG,NG)
      CALL PLANEWAVE$GETI4('NGL',NGL)
      ALLOCATE(IG(3,NGL),COEFF(NGL,2,NB),HPSI(NGL,2,NB),HBASE(NGL,2,NB))
      CALL PLANEWAVE$GETI4A('IGVEC',3*NGL,IG)
      DO ROT=0,1
        DO B=1,NB
          DO G=1,NGL
            C=SPINOR_COEFFICIENT(IG(:,G),B)
            IF(ROT.EQ.1) C=MATMUL(UNITARY,C)
            COEFF(G,:,B)=C/SQRT(VOLUME)
          END DO
        END DO
        CALL WAVES_DENSITY(NGL,NRL,2,NB,NB,F,COEFF,RHO,.FALSE.,KIN,.TRUE.,TAU)
        IF(ANY(.NOT.IEEE_IS_FINITE(RHO)).OR.ANY(.NOT.IEEE_IS_FINITE(TAU))) &
 &        ERROR STOP 'Nonfinite spinor fields'
        BOUND=MAXVAL(SQRT(SUM(TAU(:,2:4)**2,DIM=2))-TAU(:,1))
        CALL CHECK('Spinor positive-tau matrix bound',MAX(0.D0,BOUND),2.D-12)
        ERROR_RHO=0.D0
        ERROR_TAU=0.D0
        ERROR_ROTATION=0.D0
        IR=0
        DO L=0,NG-1
          DO J=0,NG-1
            DO I=FIRST-1,FIRST+NR1L-2
              IR=IR+1
              X=MATMUL(CELL,REAL([I,J,L],8)/REAL(NG,8))
              VTAU(IR,:)=VBASE+VCOS*COS(2.D0*PI*REAL(I,8)/NG)+VSIN*SIN(2.D0*PI*REAL(I,8)/NG)
              IF(ROT.EQ.1) VTAU(IR,2:4)=MATMUL(ROTATION,VTAU(IR,2:4))
              RREF=0.D0
              TREF=0.D0
              DO B=1,NB
                PSI=0.D0
                GRAD=0.D0
                DO G=1,3
                  M=[1,0,0]
                  IF(G.EQ.2) M=[0,1,0]
                  IF(G.EQ.3) M=[1,1,0]
                  Q=MATMUL(GBAS,REAL(M,8)+K)
                  PHASE=EXP(CI*DOT_PRODUCT(Q,X))
                  C=SPINOR_COEFFICIENT(M,B)
                  PSI=PSI+C*PHASE
                  DO D=1,3
                    GRAD(D,:)=GRAD(D,:)+CI*Q(D)*C*PHASE
                  END DO
                END DO
                DO A=1,4
                  RREF(A)=RREF(A)+F(B)*REAL(DOT_PRODUCT(PSI,MATMUL(SIGMA(:,:,A),PSI)),8)/VOLUME
                  DO D=1,3
                    TREF(A)=TREF(A)+0.5D0*F(B)*REAL(DOT_PRODUCT(GRAD(D,:), &
 &                            MATMUL(SIGMA(:,:,A),GRAD(D,:))),8)/VOLUME
                  END DO
                END DO
              END DO
              IF(ROT.EQ.1) THEN
                RREF(2:4)=MATMUL(ROTATION,RREF(2:4))
                TREF(2:4)=MATMUL(ROTATION,TREF(2:4))
                ERROR_ROTATION=MAX(ERROR_ROTATION,ABS(RHO(IR,1)-RBASE(IR,1)), &
 &                ABS(TAU(IR,1)-TBASE(IR,1)), &
 &                MAXVAL(ABS(RHO(IR,2:4)-MATMUL(ROTATION,RBASE(IR,2:4)))), &
 &                MAXVAL(ABS(TAU(IR,2:4)-MATMUL(ROTATION,TBASE(IR,2:4)))))
              END IF
              ERROR_RHO=MAX(ERROR_RHO,MAXVAL(ABS(RHO(IR,:)-RREF)))
              ERROR_TAU=MAX(ERROR_TAU,MAXVAL(ABS(TAU(IR,:)-TREF)))
            END DO
          END DO
        END DO
        CALL CHECK('Analytic Pauli spinor density',ERROR_RHO,2.D-12)
        CALL CHECK('Analytic Pauli spinor positive tau',ERROR_TAU,2.D-12)
        IF(ROT.EQ.0) THEN
          RBASE=RHO
          TBASE=TAU
        ELSE
          CALL CHECK('Global SU(2) spinor field covariance',ERROR_ROTATION,2.D-12)
        END IF
        INTEGRAL=SUM(TAU,DIM=1)*VOLUME/REAL(NG**3,8)
        CALL MPE$COMBINE('~','+',INTEGRAL)
        EXACT=0.D0
        DO B=1,NB
          DO G=1,3
            M=[1,0,0]
            IF(G.EQ.2) M=[0,1,0]
            IF(G.EQ.3) M=[1,1,0]
            Q=MATMUL(GBAS,REAL(M,8)+K)
            C=SPINOR_COEFFICIENT(M,B)
            IF(ROT.EQ.1) C=MATMUL(UNITARY,C)
            DO A=1,4
              EXACT(A)=EXACT(A)+0.5D0*F(B)*SUM(Q**2) &
 &                           *REAL(DOT_PRODUCT(C,MATMUL(SIGMA(:,:,A),C)),8)
            END DO
          END DO
        END DO
        CALL CHECK('Analytic integrated Pauli spinor tau',MAXVAL(ABS(INTEGRAL-EXACT)),5.D-12)
        HPSI=0.D0
        EOP=0.D0
        STRESS=0.D0
        CALL WAVES_SKALA_TAUPSI(NGL,NRL,2,NB,NB,COEFF,VTAU,HPSI,.TRUE.,F,EOP,.TRUE.,STRESS)
        ERROR_H=0.D0
        ERROR_ROTATION=0.D0
        EREF=0.D0
        SREF=0.D0
        ! Independent Fourier convolution: H(G',G)=q'.q B(G'-G)/2.
        DO B=1,NB
          DO OUTG=1,NGL
            HREF=0.D0
            QOUT=MATMUL(GBAS,REAL(IG(:,OUTG),8)+K)
            DO G=1,3
              M=[1,0,0]
              IF(G.EQ.2) M=[0,1,0]
              IF(G.EQ.3) M=[1,1,0]
              DIFF=IG(:,OUTG)-M
              VQ=0.D0
              IF(ALL(DIFF.EQ.0)) VQ=VBASE
              IF(ALL(DIFF.EQ.[1,0,0])) VQ=0.5D0*(VCOS-CI*VSIN)
              IF(ALL(DIFF.EQ.[-1,0,0])) VQ=0.5D0*(VCOS+CI*VSIN)
              IF(ROT.EQ.1) VQ(2:4)=MATMUL(ROTATION,VQ(2:4))
              BFOURIER=0.D0
              DO A=1,4
                BFOURIER=BFOURIER+VQ(A)*SIGMA(:,:,A)
              END DO
              C=SPINOR_COEFFICIENT(M,B)/SQRT(VOLUME)
              IF(ROT.EQ.1) C=MATMUL(UNITARY,C)
              Q=MATMUL(GBAS,REAL(M,8)+K)
              HREF=HREF+0.5D0*DOT_PRODUCT(QOUT,Q)*MATMUL(BFOURIER,C)
              MATRIX_ELEMENT=DOT_PRODUCT(COEFF(OUTG,:,B),MATMUL(BFOURIER,C))*VOLUME
              DO J=1,3
                DO I=1,3
                  SREF(I,J)=SREF(I,J)-0.5D0*F(B)*REAL(MATRIX_ELEMENT,8) &
 &                         *(QOUT(I)*Q(J)+QOUT(J)*Q(I))
                END DO
              END DO
            END DO
            EREF=EREF+F(B)*REAL(DOT_PRODUCT(COEFF(OUTG,:,B),HREF),8)*VOLUME
            ERROR_H=MAX(ERROR_H,MAXVAL(ABS(HPSI(OUTG,:,B)-HREF)))
            IF(ROT.EQ.1) ERROR_ROTATION=MAX(ERROR_ROTATION, &
 &              MAXVAL(ABS(HPSI(OUTG,:,B)-MATMUL(UNITARY,HBASE(OUTG,:,B)))))
          END DO
        END DO
        CALL MPE$COMBINE('~','+',EREF)
        CALL MPE$COMBINE('~','+',EOP)
        CALL MPE$COMBINE('~','+',STRESS)
        CALL MPE$COMBINE('~','+',SREF)
        CALL CHECK('Fourier Pauli tau Hamiltonian',ERROR_H,2.D-12)
        CALL CHECK('Fourier Pauli tau expectation',ABS(EOP-EREF),5.D-12)
        CALL CHECK('Fourier Pauli tau strain derivative',MAXVAL(ABS(STRESS-SREF)),5.D-12)
        IF(ROT.EQ.0) THEN
          HBASE=HPSI
        ELSE
          CALL CHECK('Global SU(2) tau Hamiltonian covariance',ERROR_ROTATION,2.D-12)
        END IF
      END DO
      DEALLOCATE(IG,COEFF,HPSI,HBASE)
    END DO
  END SUBROUTINE

  FUNCTION SPINOR_COEFFICIENT(M,B) RESULT(C)
    INTEGER,INTENT(IN) :: M(3),B
    COMPLEX(8) :: C(2)
    C=0.D0
    IF(.NOT.(ALL(M.EQ.[1,0,0]).OR.ALL(M.EQ.[0,1,0]).OR.ALL(M.EQ.[1,1,0]))) RETURN
    C(1)=FOURIER_COEFFICIENT(M,B,2)
    C(2)=CMPLX(0.11D0*(B+SUM(M)),-0.07D0*(2*B-SUM(M)),8)
  END FUNCTION

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

  SUBROUTINE PARTIAL_WAVES(X,Q,W,GW)
    REAL(8),INTENT(IN) :: X(3)
    INTEGER,INTENT(IN) :: Q
    REAL(8),INTENT(OUT) :: W(9),GW(3,9)
    REAL(8) :: U(9),GU(3,9),R2,A,N0,N1,N2,N3
    INTEGER :: I,L
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
  END SUBROUTINE

  SUBROUTINE ORBITAL_FIELDS(X,Q,C,MODE,V)
    REAL(8),INTENT(IN) :: X(3)
    INTEGER,INTENT(IN) :: Q,MODE
    COMPLEX(8),INTENT(IN) :: C(9,2,2)
    REAL(8),INTENT(OUT) :: V(5,2)
    REAL(8) :: W(9),GW(3,9),F
    COMPLEX(8) :: PSI,GRAD(3)
    INTEGER :: D,B,SPIN
    CALL PARTIAL_WAVES(X,Q,W,GW)
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

  SUBROUTINE SPINOR_SOURCES()
    INTEGER,PARAMETER :: NR=801,NP=7
    TYPE(PAW_SKALA_SOURCE) :: S,COLLINEAR
    TYPE(PAW_SKALA_SOURCE_CACHE) :: CACHE
    COMPLEX(8) :: C(9,2,2),SIGMA(2,2,4),DIR(9,9,4),DH(9,9,4),INITIAL(9,9,4)
    COMPLEX(8) :: U(2,2),ROTATED(9,2,2)
    REAL(8) :: F(2),X(3),XP(3),AE(5,4),PS(5,4),V(5,4),VP(5,4),VM(5,4),REF(5,4)
    REAL(8) :: BAR(5,4),CORE(5),DX(3),MOMENT(3,3),CONTRACTION,H,ERROR,FD,FIRSTM(3),CROSS(3),POP(5,2)
    REAL(8) :: CELL(3,3),AINV(3,3),CENTER(3),POINTS(NP,3),BARS(5,4,NP),EXPECT(5,4,NP)
    REAL(8) :: DXREF(3),MREF(3,3),CREF,AXIS(3),ROT(3,3),ANGLE,MOVED(3,3),INV(3,3),XC(3),XX(3)
    REAL(8),ALLOCATABLE :: BATCH(:,:,:)
    INTEGER :: I,J,B,D,K,L,P,MODE,COVERED,STATUS,N1,N2,N3
    LOGICAL :: REQUIRE_ACC
    CHARACTER(8) :: TEXT
    CALL GET_ENVIRONMENT_VARIABLE('CPPAW_TAU_REQUIRE_SOURCE_ACC',TEXT,STATUS=STATUS)
    REQUIRE_ACC=STATUS.EQ.0.AND.TRIM(TEXT).EQ.'1'
    SIGMA=0.D0
    SIGMA(1,1,1)=1.D0
    SIGMA(2,2,1)=1.D0
    SIGMA(1,2,2)=1.D0
    SIGMA(2,1,2)=1.D0
    SIGMA(1,2,3)=-CI
    SIGMA(2,1,3)=CI
    SIGMA(1,1,4)=1.D0
    SIGMA(2,2,4)=-1.D0
    F=[0.81D0,0.37D0]
    DO B=1,2
      DO D=1,2
        DO I=1,9
          C(I,D,B)=CMPLX(0.13D0*SIN(REAL(I+2*B+D,8)),0.09D0*COS(REAL(2*I+B+D,8)),8)
        END DO
      END DO
    END DO
    S%NR=NR
    S%LNX=3
    S%LMNX=9
    S%NDIMD=4
    S%RCUT=1.8D0
    ALLOCATE(S%R(NR),S%PHI(NR,3,2),S%CORE(NR,2),S%LOX(3),S%DM(9,9,4),S%DH(9,9,4))
    S%R=[(0.005D0*(I-1)+1.D-5,I=1,NR)]
    S%LOX=[0,1,2]
    DO D=1,2
      DO L=1,3
        S%PHI(:,L,D)=AMPLITUDE(D)*S%R**(L-1)*(1.D0-ALPHA(L,D)*S%R**2)
      END DO
    END DO
    S%CORE(:,1)=SQRT(4.D0*PI)*2.D0/PI*EXP(-2.D0*S%R)
    S%CORE(:,2)=0.5D0*S%CORE(:,1)
    CALL PAW_SKALA_SOLID_BASIS_INIT(S%BASIS,2)
    DO D=1,4
      DO J=1,9
        DO I=1,9
          DIR(I,J,D)=CMPLX(0.03D0*SIN(REAL(I+J+D,8)),0.01D0*(I-J)*D,8)
          INITIAL(I,J,D)=CMPLX(0.01D0*(I+2*J+D),0.02D0*(I-J-D),8)
        END DO
      END DO
      DO I=1,5
        BAR(I,D)=COS(REAL(I+2*D,8))
      END DO
    END DO
    AXIS=[1.D0,2.D0,-3.D0]/SQRT(14.D0)
    ANGLE=0.71D0
    U=COS(ANGLE/2.D0)*SIGMA(:,:,1)-CI*SIN(ANGLE/2.D0) &
 &     *(AXIS(1)*SIGMA(:,:,2)+AXIS(2)*SIGMA(:,:,3)+AXIS(3)*SIGMA(:,:,4))
    ROT=0.D0
    DO I=1,3
      ROT(I,I)=COS(ANGLE)
      DO J=1,3
        ROT(I,J)=ROT(I,J)+(1.D0-COS(ANGLE))*AXIS(I)*AXIS(J)
      END DO
    END DO
    ROT(1,2)=ROT(1,2)-SIN(ANGLE)*AXIS(3)
    ROT(2,1)=ROT(2,1)+SIN(ANGLE)*AXIS(3)
    ROT(1,3)=ROT(1,3)+SIN(ANGLE)*AXIS(2)
    ROT(3,1)=ROT(3,1)-SIN(ANGLE)*AXIS(2)
    ROT(2,3)=ROT(2,3)-SIN(ANGLE)*AXIS(1)
    ROT(3,2)=ROT(3,2)+SIN(ANGLE)*AXIS(1)
    DO B=1,2
      DO I=1,9
        ROTATED(I,:,B)=MATMUL(U,C(I,:,B))
      END DO
    END DO
    H=1.D-5
    ERROR=0.D0
    DO P=1,NP
      X=[0.17D0+0.07D0*P,-0.31D0+0.04D0*P,0.29D0-0.03D0*P]
      CALL SPINOR_DM(S,C,F,SIGMA)
      CALL SPINOR_ATOM_FIELDS(X,1,C,F,SIGMA,AE)
      CALL SPINOR_ATOM_FIELDS(X,2,C,F,SIGMA,PS)
      IF(P.EQ.1) FIRSTM=AE(1,2:4)
      IF(P.EQ.2) THEN
        CROSS=[FIRSTM(2)*AE(1,4)-FIRSTM(3)*AE(1,3) &
 &            ,FIRSTM(3)*AE(1,2)-FIRSTM(1)*AE(1,4),FIRSTM(1)*AE(1,3)-FIRSTM(2)*AE(1,2)]
        IF(SUM(CROSS**2).LT.1.D-16) ERROR STOP 'Source reference must be genuinely noncollinear'
      END IF
      CALL CORE_FIELDS(X,CORE)
      REF=AE-PS
      REF(:,1)=REF(:,1)+CORE
      S%DH=INITIAL
      CALL PAW_SKALA_SOURCE_POINT(S,X,V,BAR,DX,CONTRACTION)
      DH=S%DH-INITIAL
      ERROR=MAX(ERROR,MAXVAL(ABS(V-REF)),ABS(CONTRACTION-REAL(SUM(CONJG(S%DM)*DH),8)))
      COLLINEAR=S
      COLLINEAR%NDIMD=2
      COLLINEAR%DM(:,:,2)=S%DM(:,:,4)
      CALL PAW_SKALA_SOURCE_POINT(COLLINEAR,X,POP)
      ERROR=MAX(ERROR,MAXVAL(ABS(POP(:,1)-0.5D0*(V(:,1)+V(:,4)))) &
 &                    ,MAXVAL(ABS(POP(:,2)-0.5D0*(V(:,1)-V(:,4)))))
      DO K=1,3
        XP=X
        XP(K)=XP(K)+H
        CALL PAW_SKALA_SOURCE_POINT(S,XP,VP)
        XP(K)=XP(K)-2.D0*H
        CALL PAW_SKALA_SOURCE_POINT(S,XP,VM)
        FD=SUM(BAR*(VP-VM))/(2.D0*H)
        ERROR=MAX(ERROR,ABS(FD-DX(K)),MAXVAL(ABS((VP(1,:)-VM(1,:))/(2.D0*H)-V(K+1,:))))
      END DO
      S%DM=S%DM+H*DIR
      CALL PAW_SKALA_SOURCE_POINT(S,X,VP)
      S%DM=S%DM-2.D0*H*DIR
      CALL PAW_SKALA_SOURCE_POINT(S,X,VM)
      ERROR=MAX(ERROR,ABS(SUM(BAR*(VP-VM))/(2.D0*H)-REAL(SUM(CONJG(DIR)*DH),8)))
      CALL SPINOR_DM(S,ROTATED,F,SIGMA)
      CALL PAW_SKALA_SOURCE_POINT(S,X,VP)
      VM(:,1)=V(:,1)
      VM(:,2:4)=MATMUL(V(:,2:4),TRANSPOSE(ROT))
      ERROR=MAX(ERROR,MAXVAL(ABS(VP-VM)))
    END DO
    CALL CHECK('Spinor PAW source fields, adjoints and SU(2) covariance',ERROR,2.D-9)
    CALL SPINOR_DM(S,C,F,SIGMA)
    CELL=RESHAPE([2.8D0,0.D0,0.D0,0.4D0,3.1D0,0.D0,0.2D0,0.3D0,3.2D0],[3,3])
    CALL PAW_SKALA_INVERSE3(CELL,AINV)
    CENTER=[2.7D0,0.1D0,0.2D0]
    DO P=1,NP
      POINTS(P,:)=[0.17D0+0.07D0*P,-0.31D0+0.04D0*P,0.29D0-0.03D0*P]
      BARS(:,:,P)=(1.D0+0.1D0*P)*BAR
    END DO
    ERROR=0.D0
    S%DH=INITIAL
    DXREF=0.D0
    MREF=0.D0
    CREF=0.D0
    DO P=1,NP
      CALL PAW_SKALA_SOURCE_IMAGES(S,CELL,AINV,CENTER,POINTS(P,:),EXPECT(:,:,P) &
 &                                 ,BARS(:,:,P),DX,MOMENT,CONTRACTION)
      DXREF=DXREF+DX
      MREF=MREF+MOMENT
      CREF=CREF+CONTRACTION
      REF=0.D0
      DO N3=-2,2
        DO N2=-2,2
          DO N1=-2,2
            X=POINTS(P,:)-CENTER-MATMUL(CELL,REAL([N1,N2,N3],8))
            CALL PAW_SKALA_SOURCE_POINT(S,X,V)
            REF=REF+V
          END DO
        END DO
      END DO
      ERROR=MAX(ERROR,MAXVAL(ABS(EXPECT(:,:,P)-REF)))
    END DO
    DH=S%DH
    DO MODE=1,2
      CALL PAW_SKALA_SOURCE_CACHE_PREPARE(CACHE,S,CELL,AINV,CENTER,POINTS,MODE.EQ.2,16777216_8)
      S%DH=INITIAL
      DO P=1,NP
        IF(MODE.EQ.1) THEN
          CALL PAW_SKALA_SOURCE_IMAGES(S,CELL,AINV,CENTER,POINTS(P,:),V,CACHE=CACHE,IP=P)
        ELSE
          CALL PAW_SKALA_SOURCE_IMAGES(S,CELL,AINV,CENTER,POINTS(P,:),V &
 &                                     ,BARS(:,:,P),CACHE=CACHE,IP=P)
        END IF
        ERROR=MAX(ERROR,MAXVAL(ABS(V-EXPECT(:,:,P))))
      END DO
      IF(MODE.EQ.2) ERROR=MAX(ERROR,MAXVAL(ABS(S%DH-DH)))
      CALL PAW_SKALA_SOURCE_FORWARD_BATCH(CACHE,S,BATCH,COVERED)
      IF(REQUIRE_ACC.AND.COVERED.NE.NP) ERROR STOP 'Spinor forward offload did not run'
      IF(COVERED.GT.0) ERROR=MAX(ERROR,MAXVAL(ABS(BATCH-EXPECT(:,:,1:COVERED))))
    END DO
    S%DH=INITIAL
    CALL PAW_SKALA_SOURCE_BACK_BATCH(CACHE,S,BARS,.TRUE.,.TRUE.,DX,MOMENT,CONTRACTION,COVERED)
    IF(REQUIRE_ACC.AND.COVERED.NE.NP) ERROR STOP 'Spinor adjoint offload did not run'
    IF(COVERED.EQ.NP) THEN
      ERROR=MAX(ERROR,MAXVAL(ABS(S%DH-DH)),MAXVAL(ABS(DX-DXREF)) &
 &                    ,MAXVAL(ABS(MOMENT-MREF)),ABS(CONTRACTION-CREF))
    ELSE IF(COVERED.NE.0) THEN
      ERROR STOP 'Unexpected partial spinor source offload'
    END IF
    CALL CHECK('Periodic/cached/accelerated spinor PAW sources',ERROR,2.D-11)
    ! Source images and target points both follow a homogeneous cell strain.
    ERROR=0.D0
    DO K=1,3
      DO L=1,3
        MOVED=CELL
        MOVED(K,:)=CELL(K,:)+H*CELL(L,:)
        CALL PAW_SKALA_INVERSE3(MOVED,INV)
        XC=CENTER
        XC(K)=CENTER(K)+H*CENTER(L)
        XX=POINTS(1,:)
        XX(K)=XX(K)+H*XX(L)
        CALL PAW_SKALA_SOURCE_IMAGES(S,MOVED,INV,XC,XX,VP)
        MOVED(K,:)=CELL(K,:)-H*CELL(L,:)
        CALL PAW_SKALA_INVERSE3(MOVED,INV)
        XC(K)=CENTER(K)-H*CENTER(L)
        XX=POINTS(1,:)
        XX(K)=XX(K)-H*XX(L)
        CALL PAW_SKALA_SOURCE_IMAGES(S,MOVED,INV,XC,XX,VM)
        CALL PAW_SKALA_SOURCE_IMAGES(S,CELL,AINV,CENTER,POINTS(1,:),V,BARS(:,:,1),DX,MOMENT)
        FD=SUM(BARS(:,:,1)*(VP-VM))/(2.D0*H)
        ERROR=MAX(ERROR,ABS(FD-DX(K)*POINTS(1,L)+MOMENT(K,L)))
      END DO
    END DO
    CALL CHECK('Spinor PAW source strain derivative',ERROR,2.D-8)
    CALL LEGACY_SPINOR_FIELDS(C,F,SIGMA,S%DM,DIR)
  END SUBROUTINE

  SUBROUTINE LEGACY_SPINOR_FIELDS(C,F,SIGMA,DM,DIR)
    INTEGER,PARAMETER :: NR=1601,NP=7
    COMPLEX(8),INTENT(IN) :: C(9,2,2),SIGMA(2,2,4),DM(9,9,4),DIR(9,9,4)
    REAL(8),INTENT(IN) :: F(2)
    REAL(8) :: R(NR),PHI(NR,3,2),AECORE(NR),PSCORE(NR),POINTS(NP,3)
    REAL(8) :: VALUE(NP,5,4,2),PLUS(NP,5,4,2),MINUS(NP,5,4,2),BARS(NP,5,4)
    REAL(8) :: REF(5,4),CORE(5),H,ERROR,ANALYTIC,FD
    REAL(8),ALLOCATABLE :: GRID(:,:,:,:),X(:,:),W(:)
    COMPLEX(8) :: DH(9,9,4)
    INTEGER :: GID,I,L,Q,P,D,K,NPOINT,NRUSE,NANG
    CALL RADIAL$NEW('LOG',GID)
    CALL RADIAL$SETI4(GID,'NR',NR)
    CALL RADIAL$SETR8(GID,'R1',1.D-5)
    CALL RADIAL$SETR8(GID,'DEX',LOG(4.D5)/REAL(NR-1,8))
    CALL RADIAL$R(GID,NR,R)
    DO Q=1,2
      DO L=1,3
        PHI(:,L,Q)=AMPLITUDE(Q)*R**(L-1)*(1.D0-ALPHA(L,Q)*R**2)
      END DO
    END DO
    AECORE=SQRT(4.D0*PI)*2.D0/PI*EXP(-2.D0*R)
    PSCORE=0.25D0*AECORE
    DO P=1,NP
      POINTS(P,:)=[0.17D0+0.07D0*P,-0.31D0+0.04D0*P,0.29D0-0.03D0*P]
      DO D=1,4
        DO K=1,5
          BARS(P,K,D)=COS(REAL(P+K+2*D,8))
        END DO
      END DO
    END DO
    CALL PAW_SKALA_ONECENTER_POINTS(GID,NR,3,[0,1,2],9,4,DM,PHI(:,:,1),PHI(:,:,2) &
 &    ,AECORE,PSCORE,1.8D0,NP,POINTS,VALUE(:,1,:,1),VALUE(:,2:4,:,1),VALUE(:,5,:,1) &
 &    ,VALUE(:,1,:,2),VALUE(:,2:4,:,2),VALUE(:,5,:,2))
    ERROR=0.D0
    DO P=1,NP
      CALL CORE_FIELDS(POINTS(P,:),CORE)
      CORE(5)=0.D0 ! The legacy one-center diagnostic reports valence tau only.
      DO Q=1,2
        CALL SPINOR_ATOM_FIELDS(POINTS(P,:),Q,C,F,SIGMA,REF)
        REF(:,1)=REF(:,1)+MERGE(1.D0,0.25D0,Q.EQ.1)*CORE
        ERROR=MAX(ERROR,MAXVAL(ABS(VALUE(P,:,:,Q)-REF)))
      END DO
    END DO
    CALL CHECK('Legacy one-center Pauli fields against Cartesian spinors',ERROR,2.D-8)
    H=1.D-5
    CALL PAW_SKALA_ONECENTER_ADJOINT_POINTS(GID,NR,3,[0,1,2],9,4,PHI(:,:,1),PHI(:,:,2) &
 &    ,1.8D0,NP,POINTS,BARS(:,1,:),BARS(:,2:4,:),BARS(:,5,:),DH)
    CALL PAW_SKALA_ONECENTER_POINTS(GID,NR,3,[0,1,2],9,4,DM+H*DIR,PHI(:,:,1),PHI(:,:,2) &
 &    ,AECORE,PSCORE,1.8D0,NP,POINTS,PLUS(:,1,:,1),PLUS(:,2:4,:,1),PLUS(:,5,:,1) &
 &    ,PLUS(:,1,:,2),PLUS(:,2:4,:,2),PLUS(:,5,:,2))
    CALL PAW_SKALA_ONECENTER_POINTS(GID,NR,3,[0,1,2],9,4,DM-H*DIR,PHI(:,:,1),PHI(:,:,2) &
 &    ,AECORE,PSCORE,1.8D0,NP,POINTS,MINUS(:,1,:,1),MINUS(:,2:4,:,1),MINUS(:,5,:,1) &
 &    ,MINUS(:,1,:,2),MINUS(:,2:4,:,2),MINUS(:,5,:,2))
    ANALYTIC=REAL(SUM(CONJG(DIR)*DH),8)
    FD=SUM(BARS*(PLUS(:,:,:,1)-PLUS(:,:,:,2)-MINUS(:,:,:,1)+MINUS(:,:,:,2)))/(2.D0*H)
    CALL CHECK('Legacy one-center four-component density-matrix adjoint',ABS(FD-ANALYTIC),2.D-9)
    CALL PAW_SKALA_ONECENTER_GRID_SIZE(GID,NR,1.8D0,11,NRUSE,NANG,NPOINT)
    ALLOCATE(GRID(NPOINT,5,4,2),X(NPOINT,3),W(NPOINT))
    CALL PAW_SKALA_ONECENTER_PRIMITIVES(GID,NR,3,[0,1,2],9,4,DM,PHI(:,:,1),PHI(:,:,2) &
 &    ,AECORE,PSCORE,1.8D0,11,NPOINT,X,W,GRID(:,1,:,1),GRID(:,2:4,:,1),GRID(:,5,:,1) &
 &    ,GRID(:,1,:,2),GRID(:,2:4,:,2),GRID(:,5,:,2))
    ERROR=0.D0
    DO P=1,NPOINT,MAX(1,NPOINT/31)
      IF(SUM(X(P,:)**2).LT.1.D-6) CYCLE
      CALL CORE_FIELDS(X(P,:),CORE)
      CORE(5)=0.D0
      DO Q=1,2
        CALL SPINOR_ATOM_FIELDS(X(P,:),Q,C,F,SIGMA,REF)
        REF(:,1)=REF(:,1)+MERGE(1.D0,0.25D0,Q.EQ.1)*CORE
        ERROR=MAX(ERROR,MAXVAL(ABS(GRID(P,:,:,Q)-REF)))
      END DO
    END DO
    CALL CHECK('Radial-angular one-center Pauli diagnostic fields',ERROR,2.D-8)
  END SUBROUTINE

  SUBROUTINE SPINOR_DM(S,C,F,SIGMA)
    TYPE(PAW_SKALA_SOURCE),INTENT(INOUT) :: S
    COMPLEX(8),INTENT(IN) :: C(9,2,2),SIGMA(2,2,4)
    REAL(8),INTENT(IN) :: F(2)
    INTEGER :: I,J,D,B
    S%DM=0.D0
    DO D=1,4
      DO J=1,9
        DO I=1,9
          DO B=1,2
            S%DM(I,J,D)=S%DM(I,J,D)+F(B)*DOT_PRODUCT(C(I,:,B),MATMUL(SIGMA(:,:,D),C(J,:,B)))
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE

  SUBROUTINE SPINOR_ATOM_FIELDS(X,Q,C,F,SIGMA,V)
    REAL(8),INTENT(IN) :: X(3),F(2)
    INTEGER,INTENT(IN) :: Q
    COMPLEX(8),INTENT(IN) :: C(9,2,2),SIGMA(2,2,4)
    REAL(8),INTENT(OUT) :: V(5,4)
    REAL(8) :: W(9),GW(3,9)
    COMPLEX(8) :: PSI(2),GRAD(2,3)
    INTEGER :: B,D,K
    CALL PARTIAL_WAVES(X,Q,W,GW)
    V=0.D0
    DO B=1,2
      PSI=MATMUL(TRANSPOSE(C(:,:,B)),W)
      GRAD=MATMUL(TRANSPOSE(C(:,:,B)),TRANSPOSE(GW))
      DO D=1,4
        V(1,D)=V(1,D)+F(B)*REAL(DOT_PRODUCT(PSI,MATMUL(SIGMA(:,:,D),PSI)),8)
        DO K=1,3
          V(K+1,D)=V(K+1,D)+2.D0*F(B)*REAL(DOT_PRODUCT(PSI,MATMUL(SIGMA(:,:,D),GRAD(:,K))),8)
          V(5,D)=V(5,D)+0.5D0*F(B)*REAL(DOT_PRODUCT(GRAD(:,K),MATMUL(SIGMA(:,:,D),GRAD(:,K))),8)
        END DO
      END DO
    END DO
  END SUBROUTINE
END PROGRAM
