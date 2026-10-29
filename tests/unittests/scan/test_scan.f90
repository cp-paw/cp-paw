!*******************************************************************************
!****  STANDALONE TEST DRIVER FOR THE NATIVE SCAN IMPLEMENTATION            ****
!****  (SCAN$EVAL1/2/3 IN SRC/PAW_SCAN.F90)                                 ****
!****                                                                       ****
!****  BUILD/RUN: SEE RUN.SH IN THIS DIRECTORY. THE KERNEL IS COMPILED      ****
!****  WITH THIS PROGRAM AND INDEPENDENT LIBXC REFERENCE VALUES.           ****
!****                                                                       ****
!****  TESTS:                                                               ****
!****   (1) UNIFORM-GAS LIMIT: FX=1, EPSC=PW92 (VS LITERATURE VALUES)       ****
!****   (2) SINGLE-ORBITAL LIMIT (ALPHA=0,P=0): FX=H0X=1.174                ****
!****   (3) FULLY POLARIZED ONE-ELECTRON DENSITY: CORRELATION VANISHES      ****
!****   (4) SPIN-INVERSION SYMMETRY                                         ****
!****   (5) FINITE-DIFFERENCE CHECK OF ALL ANALYTIC FIRST DERIVATIVES       ****
!****   (6) EVAL2 DIAGONAL AND HESSIAN SYMMETRY                            ****
!****   (7) INDEPENDENT LIBXC ENERGY AND FIRST-DERIVATIVE VALUES           ****
!****   (8) QUADRATURE BLOCK CHAIN RULE AND WEIGHT DERIVATIVES             ****
!*******************************************************************************
      PROGRAM TEST_SCAN
      IMPLICIT NONE
      INTEGER(4)        :: NFAIL
!     **************************************************************************
      NFAIL=0
      CALL TEST_UNIFORMGAS(NFAIL)
      CALL TEST_SINGLEORBITAL(NFAIL)
      CALL TEST_ONEELECTRON(NFAIL)
      CALL TEST_SPINSYMMETRY(NFAIL)
      CALL TEST_FINITEDIFFERENCE(NFAIL)
      CALL TEST_SCAN_CLAMP_DERIVATIVES(NFAIL)
      CALL TEST_EVAL2(NFAIL)
      CALL TEST_SCAN_REFERENCE(NFAIL)
      CALL TEST_SCAN_BLOCK(NFAIL)
      IF(NFAIL.EQ.0) THEN
        WRITE(*,'(A)')'ALL TESTS PASSED'
      ELSE
        WRITE(*,'(A,I5,A)')'FAILURE: ',NFAIL,' TEST(S) FAILED'
        STOP 1
      END IF
      END PROGRAM TEST_SCAN

      SUBROUTINE TEST_SCAN_CLAMP_DERIVATIVES(NFAIL)
      IMPLICIT NONE
      INTEGER,INTENT(INOUT) :: NFAIL
      REAL(8) :: V(9),VP(9),VM(9),E,D(9),EP,EM,DUMMY(9),FD,ERROR,H
      INTEGER :: CASE_,J
      H=1.D-6
      ERROR=0.D0
      DO CASE_=1,4
        V=[0.6D0,0.1D0,0.09D0,0.01D0,0.03D0,0.D0,0.D0,1.2D0,0.3D0]
        SELECT CASE(CASE_)
        CASE(1)
          V(8)=0.37D0
          V(9)=-0.43D0
        CASE(2)
          V(8)=0.37D0
          V(9)=0.43D0
        CASE(3)
          V(2)=0.7D0
        CASE(4)
          V(2)=-0.7D0
        END SELECT
        CALL SCAN$EVAL1(V,E,D)
        DO J=1,9
          VP=V
          VM=V
          VP(J)=VP(J)+H
          VM(J)=VM(J)-H
          CALL SCAN$EVAL1(VP,EP,DUMMY)
          CALL SCAN$EVAL1(VM,EM,DUMMY)
          FD=(EP-EM)/(2.D0*H)
          ERROR=MAX(ERROR,ABS(FD-D(J))/MAX(1.D0,ABS(FD),ABS(D(J))))
        END DO
      END DO
      CALL TEST_REPORT('CLAMPED PRIMITIVE CHAIN RULE',ERROR,0.D0,2.D-8,NFAIL)
      END SUBROUTINE TEST_SCAN_CLAMP_DERIVATIVES
!
!     ..........................................................................
      SUBROUTINE TEST_REPORT(NAME,VALUE,TARGET,TOL,NFAIL)
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: NAME
      REAL(8)     ,INTENT(IN)    :: VALUE
      REAL(8)     ,INTENT(IN)    :: TARGET
      REAL(8)     ,INTENT(IN)    :: TOL
      INTEGER(4)  ,INTENT(INOUT) :: NFAIL
!     **************************************************************************
      IF(ABS(VALUE-TARGET).LE.TOL) THEN
        WRITE(*,'(A,T50,A,2E20.10)')'PASS: '//NAME,' ',VALUE,TARGET
      ELSE
        WRITE(*,'(A,T50,A,3E20.10)')'FAIL: '//NAME,' ',VALUE,TARGET &
     &                             ,VALUE-TARGET
        NFAIL=NFAIL+1
      END IF
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE TEST_UNIFORMGAS(NFAIL)
!     **************************************************************************
!     ** UNIFORM GAS: GRAD=0 AND TAU=TAUUNIF, I.E. P=0 AND ALPHA=1.           **
!     ** THEN FX=1 (LDA EXCHANGE) AND H1=0, SO EPSC = PW92 LSDA VALUES.       **
!     ** REFERENCE CORRELATION ENERGIES (HARTREE) FROM THE PW92/CEPERLEY-     **
!     ** ALDER PARAMETRIZATION: EC(RS=1)=-0.0598, EC(RS=2)=-0.0448,           **
!     ** EC(RS=5)=-0.0281 (TOLERANCE 5.D-4: PARAMETRIZATION ACCURACY)        **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT) :: NFAIL
      REAL(8),PARAMETER :: PI=4.D0*ATAN(1.D0)
      REAL(8)           :: RSTAB(3)=(/1.D0,2.D0,5.D0/)
      REAL(8)           :: ECTAB(3)=(/-0.0598D0,-0.0448D0,-0.0281D0/)
      REAL(8)           :: VAL(9),EXC,DEXC(9)
      REAL(8)           :: N,RS,TAUUNIF,EXLDA,EC
      INTEGER(4)        :: I
      CHARACTER(64)     :: NAME
!     **************************************************************************
      WRITE(*,'(A)')'== TEST 1: UNIFORM-GAS LIMIT (P=0, ALPHA=1) =='
      DO I=1,3
        RS=RSTAB(I)
        N=3.D0/(4.D0*PI*RS**3)
        TAUUNIF=0.3D0*(3.D0*PI**2)**(2.D0/3.D0)*N**(5.D0/3.D0)
        VAL(:)=0.D0
        VAL(1)=N
        VAL(8)=TAUUNIF
        CALL SCAN$EVAL1(VAL,EXC,DEXC)
!       == LDA EXCHANGE ENERGY PER VOLUME ======================================
        EXLDA=-(3.D0/(4.D0*PI))*(3.D0*PI**2*N)**(1.D0/3.D0)*N
        EC=(EXC-EXLDA)/N   ! CORRELATION ENERGY PER PARTICLE
        WRITE(NAME,'(A,F4.1,A)')'EPSC(RS=',RS,',ZETA=0) VS PW92'
        CALL TEST_REPORT(TRIM(NAME),EC,ECTAB(I),5.D-4,NFAIL)
      ENDDO
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE TEST_SINGLEORBITAL(NFAIL)
!     **************************************************************************
!     ** ALPHA=0 AND P=0: FX=H0X*GX(0)=1.174. REALIZED WITH GRAD=0,TAU=0.     **
!     ** CHECK THE EXCHANGE-ONLY PART BY SUBTRACTING THE CORRELATION          **
!     ** OBTAINED FROM A DIRECT CALL TO SCAN_CORR.                            **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT) :: NFAIL
      REAL(8),PARAMETER :: PI=4.D0*ATAN(1.D0)
      REAL(8)           :: VAL(9),EXC,DEXC(9)
      REAL(8)           :: N,EXLDA,FX
      REAL(8)           :: AC,C1,C2,C3,C4
!     **************************************************************************
      WRITE(*,'(A)')'== TEST 2: SINGLE-ORBITAL LIMIT (P=0, ALPHA=0) =='
      N=0.3D0
      VAL(:)=0.D0
      VAL(1)=N
      VAL(8)=0.D0     ! TAU=0 AND GRAD=0 -> TAUW=0 -> ALPHA=0
      CALL SCAN$EVAL1(VAL,EXC,DEXC)
      CALL SCAN_CORR(0.5D0*N,0.5D0*N,0.D0,0.D0,AC,C1,C2,C3,C4)
      EXLDA=-(3.D0/(4.D0*PI))*(3.D0*PI**2*N)**(1.D0/3.D0)*N
      FX=(EXC-AC)/EXLDA
      CALL TEST_REPORT('FX(P=0,ALPHA=0)=H0X',FX,1.174D0,1.D-10,NFAIL)
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE TEST_ONEELECTRON(NFAIL)
!     **************************************************************************
!     ** FULLY SPIN-POLARIZED SINGLE-ORBITAL DENSITY (ZETA=1, TAU=TAUW):      **
!     ** THE SCAN CORRELATION ENERGY VANISHES (NO SELF-CORRELATION).          **
!     ** USE A POINT OF A HYDROGENIC 1S ORBITAL.                              **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT) :: NFAIL
      REAL(8),PARAMETER :: PI=4.D0*ATAN(1.D0)
      REAL(8)           :: N,GN,GN2,TAUW
      REAL(8)           :: AC,C1,C2,C3,C4
!     **************************************************************************
      WRITE(*,'(A)')'== TEST 3: ONE-ELECTRON DENSITY: EPSC=0 =='
!     == HYDROGEN 1S AT R=1: N=EXP(-2R)/PI, |GRAD N|=2N ========================
      N=EXP(-2.D0)/PI
      GN=2.D0*N
      GN2=GN**2
      TAUW=GN2/(8.D0*N)
!     == FULLY POLARIZED: RUP=N, RDN=0 =========================================
      CALL SCAN_CORR(N,0.D0,GN2,TAUW,AC,C1,C2,C3,C4)
      CALL TEST_REPORT('AC(ZETA=1,ALPHA=0)',AC,0.D0,1.D-12,NFAIL)
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE TEST_SPINSYMMETRY(NFAIL)
!     **************************************************************************
!     ** EXC MUST BE INVARIANT UNDER SPIN INVERSION:                          **
!     **   RHOS->-RHOS, GRHOST->-GRHOST, TAUS->-TAUS                          **
!     ** AND THE DERIVATIVES MUST CHANGE SIGN ACCORDINGLY.                    **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT) :: NFAIL
      REAL(8)           :: VAL(9),VALM(9),EXC,DEXC(9),EXCM,DEXCM(9)
!     **************************************************************************
      WRITE(*,'(A)')'== TEST 4: SPIN-INVERSION SYMMETRY =='
      VAL=(/0.4D0,0.15D0,0.5D0,0.2D0,0.25D0,0.D0,0.D0,1.2D0,0.3D0/)
      VALM=VAL
      VALM(2)=-VAL(2)
      VALM(5)=-VAL(5)
      VALM(9)=-VAL(9)
      CALL SCAN$EVAL1(VAL,EXC,DEXC)
      CALL SCAN$EVAL1(VALM,EXCM,DEXCM)
      CALL TEST_REPORT('EXC INVARIANT',EXCM,EXC,1.D-13,NFAIL)
      CALL TEST_REPORT('DEXC(1) EVEN',DEXCM(1),DEXC(1),1.D-12,NFAIL)
      CALL TEST_REPORT('DEXC(2) ODD',DEXCM(2),-DEXC(2),1.D-12,NFAIL)
      CALL TEST_REPORT('DEXC(3) EVEN',DEXCM(3),DEXC(3),1.D-12,NFAIL)
      CALL TEST_REPORT('DEXC(5) ODD',DEXCM(5),-DEXC(5),1.D-12,NFAIL)
      CALL TEST_REPORT('DEXC(8) EVEN',DEXCM(8),DEXC(8),1.D-12,NFAIL)
      CALL TEST_REPORT('DEXC(9) ODD',DEXCM(9),-DEXC(9),1.D-12,NFAIL)
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE TEST_FINITEDIFFERENCE(NFAIL)
!     **************************************************************************
!     ** VALIDATE THE ANALYTIC FIRST DERIVATIVES AGAINST CENTRAL FINITE       **
!     ** DIFFERENCES OF THE ENERGY AT SEVERAL REPRESENTATIVE POINTS.          **
!     ** THE POINTS ARE CONSTRUCTED FROM ACTUAL SPIN GRADIENTS SO THAT        **
!     ** GRHOT2,GRHOS2,GRHOST ARE CONSISTENT, AND TAU IS CHOSEN WELL ABOVE    **
!     ** THE VON-WEIZSAECKER BOUND (NO CLAMPING ACTIVE).                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT) :: NFAIL
      INTEGER(4),PARAMETER :: NP=4
      REAL(8)           :: RUP(NP)=(/0.30D0,0.10D0,0.25D0,1.5D0 /)
      REAL(8)           :: RDN(NP)=(/0.30D0,0.25D0,0.05D0,1.2D0 /)
      REAL(8)           :: GUP(NP)=(/0.20D0,0.30D0,0.40D0,2.0D0 /)
      REAL(8)           :: GDN(NP)=(/0.20D0,0.10D0,0.15D0,1.1D0 /)
      REAL(8)           :: FUP(NP)=(/2.00D0,1.50D0,3.00D0,1.3D0 /) ! TAU FACTOR
      REAL(8)           :: FDN(NP)=(/2.00D0,2.50D0,2.00D0,1.8D0 /)
      REAL(8)           :: VAL(9),VAL1(9),EXC,DEXC(9),EP,EM,D1(9)
      REAL(8)           :: TAUWUP,TAUWDN,TUP,TDN,FD,H,ERR,SCALE
      INTEGER(4)        :: IP,I
      CHARACTER(64)     :: NAME
      LOGICAL(4)        :: TOK
!     **************************************************************************
      WRITE(*,'(A)')'== TEST 5: FINITE-DIFFERENCE CHECK OF DEXC(1:9) =='
      DO IP=1,NP
!       == CONSTRUCT CONSISTENT STATE VARIABLES ================================
!       == GRAD RHOT=GUP+GDN, GRAD RHOS=GUP-GDN (PARALLEL GRADIENTS) ==========
        VAL(1)=RUP(IP)+RDN(IP)
        VAL(2)=RUP(IP)-RDN(IP)
        VAL(3)=(GUP(IP)+GDN(IP))**2
        VAL(4)=(GUP(IP)-GDN(IP))**2
        VAL(5)=(GUP(IP)+GDN(IP))*(GUP(IP)-GDN(IP))
        VAL(6)=0.D0
        VAL(7)=0.D0
        TAUWUP=GUP(IP)**2/(8.D0*RUP(IP))
        TAUWDN=GDN(IP)**2/(8.D0*RDN(IP))
        TUP=FUP(IP)*TAUWUP
        TDN=FDN(IP)*TAUWDN
        VAL(8)=TUP+TDN
        VAL(9)=TUP-TDN
        CALL SCAN$EVAL1(VAL,EXC,DEXC)
        TOK=.TRUE.
        DO I=1,9
          IF(I.EQ.6.OR.I.EQ.7) CYCLE     ! NO LAPLACIAN DEPENDENCE
          SCALE=MAX(ABS(VAL(I)),1.D-2)
          H=1.D-6*SCALE
          VAL1=VAL
          VAL1(I)=VAL(I)+H
          CALL SCAN$EVAL1(VAL1,EP,D1)
          VAL1(I)=VAL(I)-H
          CALL SCAN$EVAL1(VAL1,EM,D1)
          FD=(EP-EM)/(2.D0*H)
          ERR=ABS(FD-DEXC(I))/MAX(ABS(FD),1.D-8)
          WRITE(NAME,'(A,I1,A,I1,A)')'POINT ',IP,' DEXC(',I,') VS FD'
          IF(ERR.LT.1.D-5) THEN
            WRITE(*,'(A,T50,A,2E20.10)')'PASS: '//TRIM(NAME),' ',DEXC(I),FD
          ELSE
            WRITE(*,'(A,T50,A,3E20.10)')'FAIL: '//TRIM(NAME),' ' &
     &                                 ,DEXC(I),FD,ERR
            NFAIL=NFAIL+1
            TOK=.FALSE.
          END IF
        ENDDO
        IF(.NOT.TOK) WRITE(*,'(A,9E12.4)')'  AT VAL=',VAL
      ENDDO
!     == LAPLACIAN DERIVATIVES MUST VANISH =====================================
      CALL TEST_REPORT('DEXC(6)=0',DEXC(6),0.D0,1.D-15,NFAIL)
      CALL TEST_REPORT('DEXC(7)=0',DEXC(7),0.D0,1.D-15,NFAIL)
      RETURN
      END
!
!     ..........................................................................
      SUBROUTINE TEST_EVAL2(NFAIL)
!     **************************************************************************
!     ** CONSISTENCY OF SCAN$EVAL2: THE (NUMERICAL) HESSIAN MUST BE           **
!     ** SYMMETRIC AND ITS DIAGONAL MUST MATCH AN INDEPENDENT FINITE          **
!     ** DIFFERENCE OF THE ENERGY.                                            **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(INOUT) :: NFAIL
      REAL(8)           :: VAL(9),VAL1(9)
      REAL(8)           :: EXC,DEXC(9),D2EXC(9,9)
      REAL(8)           :: E0,EP,EM,D1(9),H,FD2
!     **************************************************************************
      WRITE(*,'(A)')'== TEST 6: SCAN$EVAL2 CONSISTENCY =='
      VAL=(/0.6D0,0.1D0,0.4D0,0.1D0,0.15D0,0.D0,0.D0,1.5D0,0.2D0/)
      CALL SCAN$EVAL2(VAL,EXC,DEXC,D2EXC)
!     == SECOND DERIVATIVE WITH RESPECT TO RHOT FROM ENERGY DIFFERENCES =======
      H=1.D-4*VAL(1)
      VAL1=VAL
      VAL1(1)=VAL(1)+H
      CALL SCAN$EVAL1(VAL1,EP,D1)
      VAL1(1)=VAL(1)-H
      CALL SCAN$EVAL1(VAL1,EM,D1)
      CALL SCAN$EVAL1(VAL,E0,D1)
      FD2=(EP-2.D0*E0+EM)/H**2
      CALL TEST_REPORT('D2EXC(1,1) VS FD OF ENERGY',D2EXC(1,1),FD2 &
     &                ,1.D-4*MAX(ABS(FD2),1.D0),NFAIL)
      CALL TEST_REPORT('HESSIAN SYMMETRY (1,8)',D2EXC(1,8),D2EXC(8,1) &
     &                ,1.D-12,NFAIL)
      RETURN
      END
