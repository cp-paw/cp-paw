!
!     ..................................................................
!schroedinger$spherical
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
      LOGICAL(4)                 :: THOM
!     ************************************************************************
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
!     ==================================================================
!     == SPIN ORBIT COUPLING ===========================================
!     ==================================================================
      IF (SO.EQ.0) THEN
        SOFACTOR=0.D0
      ELSE IF(SO.EQ.1) THEN
        SOFACTOR=REAL(L,KIND=8)       ! PARALLEL SPIN AND ORBIT
      ELSE IF(SO.EQ.-1) THEN
        SOFACTOR=REAL(-L-1,KIND=8)    ! ANTIPARALLELSPIN AND ORBIT
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
     &    -RDPRIME(2:)*SOFACTOR/R(2:) &
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
!schroedinger$nonspherical_bound_c
      SUBROUTINE RADIAL$NONSPHBOUND(GID,NR,NDIMD,LMX,LMRX,POT,DREL,G,E &
     &                             ,NPHI,EB,PHI,TPHI,SPHI,TSPHI,TOK)
!     **                                                                  **
!     **  SOLVES THE RELATIVISTIC RADIAL DIRAC EQUATION FOR THE           **
!     **  LARGE COMPONENT. DREL=1/MREL-1/M0 IS A MEASURE FOR THE          **
!     **  RELATIVISTIC EFFECTS, WHERE MREL=M0+(E-V)/2C^2 AND M0 IS THE    **
!     **  REST MASS. V=POT*Y0 IS THE POTENTIAL.                           **
!     **  SPIN ORBIT COUPLING IS MEASURED BY SO WHICH CAN HAVE THE VALUES:**
!     **    SO=0 NO-SPIN ORBIT COUPLING                                   **
!     **    SO=1 PARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM               **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE        **
!     **  HOMOGENEOUS SOLUTION.                                           **
!     **                                                                  **
!     **  THE SOLUTIONS ARE CALCULATED TO LINEAR ORDER IN DE              **
!     **  WHERE E+DE IS THE NEW ENERGY OF THE WAVE FUNCTIONS              **
!     **                                                                  **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM   **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                     **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                               **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                         **
!     **                                                                  **
!     **  THE NONRELATIVISTIC SOLUTION IS OBTAINED BY SETTINH DREL=0      **
!     **                                                                  **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE      **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                   **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: NDIMD   ! (1,2,4)#(SPINOR COMPONENTS)
      INTEGER(4) ,INTENT(IN)     :: LMX     ! X#(WAVE FUNCTION ANGULAR MOMENTA
      INTEGER(4) ,INTENT(IN)     :: LMRX    ! X#(POTENTIAL ANGULAR MOMENTA)
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: G(NR,LMX,2)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: E       !ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR,LMRX,NDIMD) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NPHI
      COMPLEX(8) ,INTENT(OUT)    :: PHI(NR,LMX,2,NPHI) ! WAVE-FUNCTION
      COMPLEX(8) ,INTENT(OUT)    :: TPHI(NR,LMX,2,NPHI) ! WAVE-FUNCTION
      COMPLEX(8) ,INTENT(OUT)    :: SPHI(NR,LMX,2,NPHI) ! WAVE-FUNCTION
      COMPLEX(8) ,INTENT(OUT)    :: TSPHI(NR,LMX,2,NPHI) ! KINETIC ENERGY * WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: EB(NPHI)
      LOGICAL(4) ,INTENT(OUT)    :: TOK
      REAL(8)                    :: DE(NPHI)
      INTEGER(4)                 :: LX
      INTEGER(4)                 :: LOX(LMX,2) ! ANGULAR MOMENTA
      INTEGER(4)                 :: LM,LM1,LM2,LM3,L,M,IM,IS,IB,IS1,IS2
      REAL(8)                    :: A(NR)
      REAL(8)                    :: B(NR)
      COMPLEX(8)                 :: C(NR,LMX,2,LMX,2)
      COMPLEX(8)                 :: CPOT(NR,LMX,2,LMX,2)
      COMPLEX(8)                 :: CR(LMX,2,LMX,2)
      COMPLEX(8)                 :: CA(LMX,2,LMX,2)
      COMPLEX(8)                 :: D(NR,LMX,2) 
      REAL(8)                    :: DR(NR,LMX,2) 
      REAL(8)                    :: PHIR(NR,LMX,2,2*LMX)
      REAL(8)                    :: DPHI(NR,LMX,2)
      REAL(8)                    :: R(NR) 
      REAL(8)                    :: AUX(NR),AUX1(NR),AUX2(NR) 
      COMPLEX(8)                 :: CSVAR
      REAL(8)                    :: PI
      REAL(8)                    :: Y0
      REAL(8)                    :: CG
      REAL(8)                    :: RDPRIME(NR)
      INTEGER(4)                 :: IR
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20
      REAL(8)                    :: SWKB(NR),SVAR
      INTEGER(4)                 :: IRCL,IROUT
      COMPLEX(8),PARAMETER       :: CI=(0.D0,1.D0)
      COMPLEX(8)                 :: CLS(LMX,2,LMX,2) ! SPIN ORBIT MATRIX (LS)
      LOGICAL(4)                 :: TCHK
!     ************************************************************************
      TOK=.FALSE.
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      LX=INT(SQRT(REAL(LMX))-1.D0)
      IF((LX+1)**2.NE.LMX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('RADIAL$NONSPHBOUND')
      END IF
      LM=0
      DO L=0,LX
        DO M=1,2*L+1
          LM=LM+1
          LOX(LM,:)=L
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DETERMINE CLASSICAL TURNING POINT                           ==
!     ==================================================================
!     == IROUT MUST BE SMALLER OR EQUAL TO NR-1 BECAUSE INHOMOGENEITY
!     == ON THE LAST GRID POINT DOES NOT CONTRIBUTE TO DIFFERENTIAL 
!     == EQUATION. IROUT IS SMALLER OR EQUAL TO IRCL
      IRCL=NR-3   
!       == R(IRCL) IS THE FIRST GRID POINT INSIDE THE CLASSICAL TURNING POINT
      IF(E.LT.POT(IRCL,1,1)*Y0) THEN
        DO IR=IRCL,1,-1
          IF(E.GT.POT(IR,1,1)*Y0) THEN
            IRCL=IR
            EXIT
          END IF
        ENDDO
      END IF
!
!     ==================================================================
!     ==  DETERMINE OUTERMOST POINT 'IROUT' FOR INWARD INTEGRATION    ==
!     ==================================================================
!     == USE WKB SOLUTION FOR THE SCHR.GL. FOR A CONSTANT POTENTIAL AND L=0
!     == TO ESTIMATE FACTOR FROM RCL TO OUTERMOST POINT
      AUX(:IRCL)=0.D0
      AUX(IRCL+1:)=SQRT(MAX(0.D0,2.D0*(POT(IRCL+1:,1,1)*Y0-E)))
      CALL RADIAL$INTEGRATE(GID,NR,AUX,SWKB)
      SWKB(:)=SWKB(:)-LOG(R(:))
!     == DETERMINE IROUT WHERE THE WAVE FUNCTION CAN GROW BY A FACTOR 
!     == OF XMAX FROM THE CLASSICAL TURNING POINT
      SVAR=LOG(XMAX)
      IROUT=NR-1
      DO IR=IRCL,NR
        IF(SWKB(IR).GT.SVAR) THEN
          IROUT=IR-1
          EXIT
        END IF
      ENDDO
!
!     ==================================================================
!     ==  PREPARE POTENTIAL-INDEPENDENT ARRAYS                        ==
!     ==================================================================
!     == A*D2F/DR2+B*DF/DR+C*F=D
      A(:)=1.D0+DREL(:)
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0*(1.D0+DREL(2:))/R(2:)+RDPRIME(2:)
      B(1)=B(2)
      D(:,:,:)=-2.D0*G(:,:,:)
!
!     ==================================================================
!     ==  COUPLING BETWEEN WAVE FUNCTION COMPONENTS VIA POTENTIAL     ==
!     ==================================================================
      C(:,:,:,:,:)=(0.D0,0.D0)
      DO LM1=1,LMX
        DO LM2=1,LMX
          AUX(:)=0.D0
          DO LM3=1,LMRX
            CALL CLEBSCH(LM1,LM2,LM3,CG)
            IF(CG.EQ.0.D0) CYCLE
            AUX=CG*POT(:,LM3,1)
            IF(LM3.EQ.1) AUX(IROUT+1:)=AUX(IROUT)  ! CONSTANT POTENTIAL BEYOND IROUT
            C(:,LM1,1,LM2,1)=C(:,LM1,1,LM2,1)+AUX
            C(:,LM1,2,LM2,2)=C(:,LM1,2,LM2,2)+AUX
            IF(NDIMD.GT.1) THEN
              IF(NDIMD.EQ.2) THEN
                AUX=CG*POT(:,LM3,2)
                C(:,LM1,1,LM2,1)=C(:,LM1,1,LM2,1)+AUX
                C(:,LM1,2,LM2,2)=C(:,LM1,2,LM2,2)-AUX
              ELSE
                AUX=CG*POT(:,LM3,2)
                C(:,LM1,1,LM2,2)=C(:,LM1,1,LM2,2)+AUX
                C(:,LM1,2,LM2,1)=C(:,LM1,2,LM2,1)-AUX
                AUX=CG*POT(:,LM3,3)
                C(:,LM1,1,LM2,2)=C(:,LM1,1,LM2,2)-CI*AUX
                C(:,LM1,2,LM2,1)=C(:,LM1,2,LM2,1)+CI*AUX
                AUX=CG*POT(:,LM3,4)
                C(:,LM1,1,LM2,1)=C(:,LM1,1,LM2,1)+AUX
                C(:,LM1,2,LM2,2)=C(:,LM1,2,LM2,2)-AUX
              END IF
            END IF
          ENDDO
        ENDDO
      ENDDO
      CPOT=C    ! STORE TO EVALUATE KINETIC ENERGY
      C=-2.D0*C
!
!     ==================================================================
!     ==  KINETIC ENERGY TERM TO C                                    ==
!     ==================================================================
      LM=0
      DO L=0,LX
        AUX(1)=0.D0
        AUX(2:)=-(1.D0+DREL(2:))/R(2:)**2 * REAL(L*(L+1),KIND=8)+2.D0*E
        DO IM=1,2*L+1
          LM=LM+1
          DO IS=1,2
            C(:,LM,IS,LM,IS)=C(:,LM,IS,LM,IS)+AUX(:)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  ADD SPIN ORBIT COUPLING TO E                                ==
!     ==================================================================
      CALL RADIAL_LS(LMX,CLS)
!CLS=0.D0
      AUX(2:)=-RDPRIME(2:)/R(2:)
      AUX(1)=AUX(2)
      DO LM1=1,LMX
        DO IS1=1,2
          DO LM2=1,LMX
            DO IS2=1,2
              IF(ABS(CLS(LM1,IS1,LM2,IS2)).LT.1.D-10) CYCLE
              C(:,LM1,IS1,LM2,IS2)=C(:,LM1,IS1,LM2,IS2)+AUX(:)*CLS(LM1,IS1,LM2,IS2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!     ==  AVOID DIVIDE ZEROBYZERO
      C(1,:,:,:,:)=C(2,:,:,:,:)
!
!     ==================================================================
!     ==  DETERMINE BOUND STATES                                      ==
!     ==================================================================
      CALL RADIAL_XXXC(GID,NR,2*LMX,IRCL,IROUT,LOX,A,B,C,D,NPHI,DE,PHI,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$STOP('RADIAL_XXXC FINISHED WITH ERROR')
        CALL ERROR$STOP('RADIAL$NONSPHBOUND')
      END IF
!
!     ==================================================================
!     ==  DETERMINE SMALL COMPONENT                                   ==
!     ==================================================================
       SPHI=(0.D0,0.D0)
!      CALL RADIAL_SMALLCOMPONENT(GID,NR,LMX,NPHI,DREL,PHI,SPHI)
!
!!$      CALL RADIAL_SP(LMX,CR,CA)
!!$      CALL CONSTANTS$GET('C',SVAR)
!!$      A(:)=(1.D0+DREL(:))/(2.D0*SVAR)
!!$      SPHI(:,:,:,:)=(0.D0,0.D0)
!!$      DO IB=1,NPHI
!!$        DO IS=1,2
!!$          DO LM=1,LMX
!!$            CALL RADIAL$DERIVE(GID,NR,REAL(PHI(:,LM,IS,IB)),AUX1)
!!$            CALL RADIAL$DERIVE(GID,NR,AIMAG(PHI(:,LM,IS,IB)),AUX2)
!!$            DPHI(:,LM,IS)=CMPLX(AUX1,AUX2)
!!$          ENDDO
!!$        ENDDO
!!$!
!!$        DO IS2=1,2
!!$          DO LM2=1,LMX
!!$            DO IS1=1,2
!!$              DO LM1=1,LMX
!!$                CSVAR=CR(LM1,IS1,LM2,IS2)
!!$                IF(CSVAR.NE.0.D0) THEN
!!$                  SPHI(:,LM1,IS1,IB)=SPHI(:,LM1,IS1,IB)+CSVAR*DPHI(:,LM2,IS2)
!!$                END IF
!!$                CSVAR=CA(LM1,IS1,LM2,IS2)
!!$                IF(CSVAR.NE.0.D0) THEN
!!$                  SPHI(:,LM1,IS1,IB)=SPHI(:,LM1,IS1,IB)+CSVAR*PHI(:,LM2,IS2,IB)
!!$                END IF
!!$              ENDDO
!!$            ENDDO
!!$          ENDDO
!!$        ENDDO
!!$!
!!$        DO LM=1,LMX
!!$          DO IS=1,2
!!$            SPHI(:,LM,IS,IB)=A(:)*SPHI(:,LM,IS,IB)                  
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!
!     ==================================================================
!     ==  SHIFT ENERGIES                                              ==
!     ==================================================================
      EB(:)=DE(:)+E
!
!     ==================================================================
!     ==  DETERMINE TPHI AND TSPHI                                              ==
!     ==================================================================
!     -- BETTER DIRECTLY WORK OUT THE KINETIC ENERGY BECAUSE  THE 
!     -- SCHROEDINGER EQUATION IS FULFILLED ONLY TO FIRST ORDER IN DE
      DO IB=1,NPHI
        TPHI(:,:,:,IB)=EB(IB)*PHI(:,:,:,IB)
        TSPHI(:,:,:,IB)=EB(IB)*SPHI(:,:,:,IB)
        DO IS1=1,2
          DO LM1=1,LMX
            DO IS2=1,2
              DO LM2=1,LMX
                TPHI(:,LM1,IS1,IB)=TPHI(:,LM1,IS1,IB) &
     &                         -CPOT(:,LM1,IS1,LM2,IS2)*PHI(:,LM2,IS2,IB)
                TSPHI(:,LM1,IS1,IB)=TSPHI(:,LM1,IS1,IB) &
     &                         -CPOT(:,LM1,IS1,LM2,IS2)*SPHI(:,LM2,IS2,IB)
              ENDDO
            ENDDO  
          ENDDO
        ENDDO
      ENDDO
!
      TOK=.TRUE.
      RETURN
      END SUBROUTINE RADIAL$NONSPHBOUND
!
!     ........................................................................
!schroedinger_smallcomponent_c
      SUBROUTINE RADIAL_SMALLCOMPONENT(GID,NR,LMX,NPHI,DREL,PHI,SPHI)
!     **                                                                    **
!     **                                                                    **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: LMX
      INTEGER(4),INTENT(IN) :: NPHI
      REAL(8)   ,INTENT(IN) :: DREL(NR)
      COMPLEX(8),INTENT(IN) :: PHI(NR,LMX,2,NPHI)
      COMPLEX(8),INTENT(OUT):: SPHI(NR,LMX,2,NPHI)
      COMPLEX(8)            :: CR(LMX,2,LMX,2)
      COMPLEX(8)            :: CA(LMX,2,LMX,2)
      REAL(8)               :: A(NR)
      REAL(8)               :: AUX1(NR),AUX2(NR),SVAR
      REAL(8)               :: DPHI(NR,LMX,2)
      COMPLEX(8)            :: CSVAR
      INTEGER(4)            :: IB,IS,LM,IS1,IS2,LM1,LM2
!     ************************************************************************
      CALL RADIAL_SP(LMX,CR,CA)
      CALL CONSTANTS$GET('C',SVAR)
      A(:)=(1.D0+DREL(:))/(2.D0*SVAR)
      SPHI(:,:,:,:)=(0.D0,0.D0)
      DO IB=1,NPHI
        DO IS=1,2
          DO LM=1,LMX
            CALL RADIAL$DERIVE(GID,NR,REAL(PHI(:,LM,IS,IB)),AUX1)
            CALL RADIAL$DERIVE(GID,NR,AIMAG(PHI(:,LM,IS,IB)),AUX2)
            DPHI(:,LM,IS)=CMPLX(AUX1,AUX2)
          ENDDO
        ENDDO
!
        DO IS2=1,2
          DO LM2=1,LMX
            DO IS1=1,2
              DO LM1=1,LMX
                CSVAR=CR(LM1,IS1,LM2,IS2)
                IF(CSVAR.NE.0.D0) THEN
                  SPHI(:,LM1,IS1,IB)=SPHI(:,LM1,IS1,IB)+CSVAR*DPHI(:,LM2,IS2)
                END IF
                CSVAR=CA(LM1,IS1,LM2,IS2)
                IF(CSVAR.NE.0.D0) THEN
                  SPHI(:,LM1,IS1,IB)=SPHI(:,LM1,IS1,IB)+CSVAR*PHI(:,LM2,IS2,IB)
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
        DO LM=1,LMX
          DO IS=1,2
            SPHI(:,LM,IS,IB)=A(:)*SPHI(:,LM,IS,IB)                  
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     .......................................................................
      SUBROUTINE RADIAL_LSOLD(LMX,C)
!     **                                                                    **
!     **  SPIN ORBIT MATRIX ELEMENTS                                        **
!     **                                                                    **
!     **  MATRIX ELEMENTS OF L*SIGMA IN REAL SPHERICAL HARMONICS            **
!     **  WHERE SIGMA ARE THE PAULI MATRICES AND L ARE THE ANGULAR MOMENTA  **
!     **                                                                    **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX
      COMPLEX(8),INTENT(OUT):: C(LMX,2,LMX,2)
      REAL(8)               :: C1(LMX,2,LMX,2)
      INTEGER(4)            :: LX
      INTEGER(4)            :: L,M,M1,M2,LM,LM1P,LM1M,LM2P,LM2M,IS,LM1,LM2,IS1,IS2
      REAL(8)               :: SVAR
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      COMPLEX(8)            :: FAC1,FAC2,FAC3,FAC4
      LOGICAL               :: TPRINT=.FALSE.
      LOGICAL               :: TTEST=.TRUE.
      REAL(8)               :: SQR2IN
      COMPLEX(8),ALLOCATABLE :: H(:,:),U(:,:)
      REAL(8)   ,ALLOCATABLE :: E(:)
!     ************************************************************************
      LX=INT(SQRT(real(LMX))-1.d0) ! lmx=(lx+1)**2
      if((lx+1)**2.ne.lmx) then
        call error$msg('lmx does not correspond to a full shell')
        call error$msg('or rounding errors produced incorrect results')
        call error$i4val('lmx',lmx)
        call error$i4val('lx',lx)
        call error$stop('RADIAL_LSOLD')
      end if
      SQR2IN=1.D0/SQRT(2.D0)
!
!     ========================================================================
!     ==  MATRIX ELEMENTS <L,M,S|\SIGMA*L|L,M,S> WHERE |L,M,S> ARE          ==
!     ==  SPHERICAL HARMONICS (ANGULAR MOMENTUM EIGENSTATES)                ==
!     ========================================================================
      C1=0.D0
      LM=0
      DO L=0,LX 
        DO M=-L,L
          LM=LM+1
          C1(LM,1,LM,1)= REAL(M)
          C1(LM,2,LM,2)=-REAL(M)
          IF(M.EQ.L) CYCLE
          SVAR=SQRT( REAL((L-M)*(L+M+1)) )
          C1(LM,1,LM+1,2)=SVAR
          C1(LM+1,2,LM,1)=SVAR
        ENDDO
      ENDDO
!
!     ==================================================================
!     == TRANSFORM TO REAL SPHERICAL HARMONICS                        ==
!     ==================================================================
      IF(TPRINT) THEN
        DO L=1,LX
          LM=L**2
          DO IS=1,2
            DO M=1,2*L+1
              WRITE(*,FMT='(20F8.3)')C1(LM+M,IS,LM+1:LM+2*L+1,:)
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==================================================================
!     == TRANSFORM TO REAL SPHERICAL HARMONICS                        ==
!     ==================================================================
      C=(0.D0,0.D0)
      DO L=1,LX    !NO S.O. VOUPLING IN S-CHANNEL
        LM1P=L**2+L+1
        LM1M=LM1P
!       == M1=M2=0
        C(LM1M,:,LM1M,:)=CMPLX(C1(LM1M,:,LM1M,:))
!       == (M1=0 AND M2.NEQ.0) OR (N1.NEQ.0 AND M2=0)
        LM2P=L**2+L+1
        LM2M=LM2P
        DO M2=1,L
          LM2P=LM2P+1
          LM2M=LM2M-1
          FAC1=CMPLX(SQR2IN)
          FAC2=CMPLX(SQR2IN*(-1.D0)**M2)
          C(LM1M,:,LM2P,:)=FAC1*C1(LM1M,:,LM2P,:)+FAC2*C1(LM1M,:,LM2M,:)
          C(LM2P,:,LM1M,:)=FAC1*C1(LM2P,:,LM1M,:)+FAC2*C1(LM2M,:,LM1M,:)
          FAC1=-CI*FAC1
          FAC2=CI*FAC2
          C(LM1M,:,LM2M,:)=+FAC1*C1(LM1M,:,LM2P,:)+FAC2*C1(LM1M,:,LM2M,:)
          C(LM2M,:,LM1M,:)=-FAC1*C1(LM2P,:,LM1M,:)-FAC2*C1(LM2M,:,LM1M,:)
        ENDDO
!       == M1.NEQ.0, M2.NEQ.0 ==============================================
        DO M1=1,L
          LM1P=LM1P+1
          LM1M=LM1M-1
          LM2P=L**2+L+1
          LM2M=LM2P
          DO M2=1,L
            LM2P=LM2P+1
            LM2M=LM2M-1
            FAC1=CMPLX(0.5D0)
            FAC2=CMPLX(0.5D0*(-1.D0)**(M1+M2))
            FAC3=CMPLX(0.5D0*(-1.D0)**M2)
            FAC4=CMPLX(0.5D0*(-1.D0)**M1)
            C(LM1P,:,LM2P,:)=+FAC1*C1(LM1P,:,LM2P,:)+FAC2*C1(LM1M,:,LM2M,:) &
    &                        +FAC3*C1(LM1P,:,LM2M,:)+FAC3*C1(LM1M,:,LM2P,:)
            C(LM1M,:,LM2M,:)=+FAC1*C1(LM1P,:,LM2P,:)+FAC2*C1(LM1M,:,LM2M,:) &
    &                        -FAC3*C1(LM1P,:,LM2M,:)-FAC3*C1(LM1M,:,LM2P,:)
            FAC1=FAC1*CI
            FAC2=FAC2*CI
            FAC3=FAC3*CI
            FAC4=FAC4*CI
            C(LM1P,:,LM2M,:)=-FAC1*C1(LM1P,:,LM2P,:)+FAC2*C1(LM1M,:,LM2M,:) &
    &                        +FAC3*C1(LM1P,:,LM2M,:)-FAC3*C1(LM1M,:,LM2P,:)
            C(LM1M,:,LM2P,:)=+FAC1*C1(LM1P,:,LM2P,:)-FAC2*C1(LM1M,:,LM2M,:) &
    &                        +FAC3*C1(LM1P,:,LM2M,:)-FAC3*C1(LM1M,:,LM2P,:)
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     == TRANSFORM TO REAL SPHERICAL HARMONICS                        ==
!     ==================================================================
      IF(TPRINT) THEN
        DO L=1,LX
          ALLOCATE(H(2*(2*L+1),2*(2*L+1)))
          ALLOCATE(U(2*(2*L+1),2*(2*L+1)))
          ALLOCATE(E(2*(2*L+1)))
          LM=L**2
          DO IS=1,2
            DO M=1,2*L+1
              WRITE(*,FMT='(20("(",F10.3,",",F10.3,")"))')C(LM+M,IS,LM+1:LM+2*L+1,:)
            ENDDO
          ENDDO
          H(:2*L+1,:2*L+1)=C(LM+1:LM+2*L+1,1,LM+1:LM+2*L+1,1)
          H(:2*L+1,2*L+2:)=C(LM+1:LM+2*L+1,1,LM+1:LM+2*L+1,2)
          H(2*L+2:,:2*L+1)=C(LM+1:LM+2*L+1,2,LM+1:LM+2*L+1,1)
          H(2*L+2:,2*L+2:)=C(LM+1:LM+2*L+1,2,LM+1:LM+2*L+1,2)
          CALL LIB$DIAGC8(2*(2*L+1),H,E,U)
          PRINT*,'E',L,E
          DEALLOCATE(H)
          DEALLOCATE(E)
          DEALLOCATE(U)
        ENDDO
      END IF
!
!     ==================================================================
!     == TRANSFORM TO REAL SPHERICAL HARMONICS                        ==
!     ==================================================================
      IF(TTEST) THEN
        DO LM1=1,LMX
          DO IS1=1,2
            DO LM2=1,LMX
              DO IS2=1,2
                IF(C(LM1,IS1,LM2,IS2).NE.CONJG(C(LM2,IS2,LM1,IS1))) THEN
                  CALL ERROR$MSG('SPIN ORBIT IS NOT HERMITEAN')
                  CALL ERROR$STOP('RADIAL$LS')
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
      RETURN
      END SUBROUTINE RADIAL_LSOLD
!
!     .......................................................................
      SUBROUTINE RADIAL_LS(LMX,C)
!     **                                                                    **
!     **  SPIN ORBIT MATRIX ELEMENTS                                        **
!     **                                                                    **
!     **  MATRIX ELEMENTS OF L*SIGMA IN REAL SPHERICAL HARMONICS            **
!     **  WHERE SIGMA ARE THE PAULI MATRICES AND L ARE THE ANGULAR MOMENTA  **
!     **                                                                    **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX
      COMPLEX(8),INTENT(OUT):: C(LMX,2,LMX,2)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      LOGICAL               :: TTEST=.FALSE.
      COMPLEX(8)            :: LXMAT(LMX,LMX),LYMAT(LMX,LMX),LZMAT(LMX,LMX)
      REAL(8)               :: SVAR
      COMPLEX(8)            :: CTEST(LMX,2,LMX,2)
!     ************************************************************************
      CALL SPHERICAL$L(LMX,LXMAT,LYMAT,LZMAT)
!
!     ========================================================================
!     ==  MATRIX ELEMENTS <L,M,S|\SIGMA*L|L,M,S> WHERE |L,M,S> ARE          ==
!     ==  SPHERICAL HARMONICS (ANGULAR MOMENTUM EIGENSTATES)                ==
!     ========================================================================
      C(:,1,:,1)=LZMAT(:,:)
      C(:,1,:,2)=LXMAT(:,:)-CI*LYMAT(:,:)
      C(:,2,:,1)=LXMAT(:,:)+CI*LYMAT(:,:)
      C(:,2,:,2)=-LZMAT(:,:)
!
!     ==================================================================
!     == TRANSFORM TO REAL SPHERICAL HARMONICS                        ==
!     ==================================================================
      IF(TTEST) THEN
        SVAR=MAXVAL(ABS(C(:,1,:,1)-TRANSPOSE(CONJG(C(:,1,:,1)))))
        SVAR=MAX(SVAR,MAXVAL(ABS(C(:,1,:,2)-TRANSPOSE(CONJG(C(:,2,:,1))))))
        SVAR=MAX(SVAR,MAXVAL(ABS(C(:,2,:,2)-TRANSPOSE(CONJG(C(:,2,:,2))))))
        IF(SVAR.GT.1.D-10) THEN
          CALL ERROR$MSG('SPIN ORBIT IS NOT HERMITEAN')
          CALL ERROR$R8VAL('DEVIATION',SVAR)
          CALL ERROR$STOP('RADIAL$LS')
        END IF
        CALL RADIAL_LSOLD(LMX,CTEST)
        SVAR=MAXVAL(ABS(C-CTEST))
        IF(SVAR.GT.1.D-7) THEN
          CALL ERROR$MSG('SPIN ORBIT DOES NOT AGREE WITH PREVIOUS')
          CALL ERROR$R8VAL('DEVIATION',SVAR)
          CALL ERROR$STOP('RADIAL$LS')
        END IF
      END IF
      RETURN
      END SUBROUTINE RADIAL_LS
!
!     .......................................................................
      SUBROUTINE RADIAL_SP(LMX,CR,CA)
!     **                                                                   **
!     **  CALCULATES THE MATRICES REQUIRED TO EVALUATE THE SMALL COMPONENT **
!     **  (POSITRONS) FROM THE SOLUTION OF THE LARGE COMPONENT (ELECTRONS) **
!     **                                                                   **
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LMX
      COMPLEX(8),INTENT(OUT):: CR(LMX,2,LMX,2)
      COMPLEX(8),INTENT(OUT):: CA(LMX,2,LMX,2)
      COMPLEX(8),PARAMETER  :: CI=(0.D0,1.D0)
      REAL(8)               :: X(LMX,LMX),Y(LMX,LMX),Z(LMX,LMX)
      COMPLEX(8)            :: LX(LMX,LMX),LY(LMX,LMX),LZ(LMX,LMX)
      COMPLEX(8)            :: RCROSSL(LMX,LMX,3)
!     ***********************************************************************
      CALL SPHERICAL$ER(LMX,X,Y,Z)
      CALL SPHERICAL$L(LMX,LX,LY,LZ)
!
!     =======================================================================
!     == DETERMINE THE RADIAL PART CR                                      ==
!     =======================================================================
      CR(:,1,:,1)=CMPLX(Z(:,:),0.D0)
      CR(:,1,:,2)=CMPLX(X(:,:),-Y(:,:))
      CR(:,2,:,1)=CMPLX(X(:,:),Y(:,:))
      CR(:,2,:,2)=CMPLX(-Z(:,:),0.D0)
      CR(:,:,:,:)=0.5D0*CR(:,:,:,:)
!
!     =======================================================================
!     == DETERMINE THE RADIAL PART CA                                      ==
!     =======================================================================
      RCROSSL(:,:,1)=MATMUL(Y,LZ)-MATMUL(Z,LY)
      RCROSSL(:,:,2)=MATMUL(Z,LX)-MATMUL(X,LZ)
      RCROSSL(:,:,3)=MATMUL(X,LY)-MATMUL(Y,LX)
      CA(:,1,:,1)=RCROSSL(:,:,3)
      CA(:,1,:,2)=RCROSSL(:,:,1)-CI*RCROSSL(:,:,2)
      CA(:,2,:,1)=RCROSSL(:,:,1)+CI*RCROSSL(:,:,2)
      CA(:,2,:,2)=-RCROSSL(:,:,3)
      CA(:,:,:,:)=-0.5D0*CA(:,:,:,:)
      RETURN
      END
!
!     ..................................................................
      SUBROUTINE RADIAL_XXXC(GID,NR,NF,IRMATCH,IROUT,LOX,A,B,C,D,NPHI,DE,PHI,TOK)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID
      INTEGER(4),INTENT(IN) :: NR
      INTEGER(4),INTENT(IN) :: NF
      INTEGER(4),INTENT(IN) :: IRMATCH
      INTEGER(4),INTENT(IN) :: IROUT
      INTEGER(4),INTENT(IN) :: LOX(NF)
      REAL(8)   ,INTENT(IN) :: A(NR)
      REAL(8)   ,INTENT(IN) :: B(NR)
      COMPLEX(8),INTENT(IN) :: C(NR,NF,NF)
      COMPLEX(8),INTENT(IN) :: D(NR,NF)
      INTEGER(4),INTENT(IN) :: NPHI
      REAL(8)   ,INTENT(OUT):: DE(NPHI)
      COMPLEX(8),INTENT(OUT):: PHI(NR,NF,NPHI)
      LOGICAL(4),INTENT(OUT):: TOK
      COMPLEX(8)            :: ALLPHIL(NR,NF,NF)
      COMPLEX(8)            :: ALLPHIR(NR,NF,NF)
      COMPLEX(8)            :: PHIL(NR,NF,NPHI)
      COMPLEX(8)            :: PHIR(NR,NF,NPHI)
      COMPLEX(8)            :: PHIL_DOT(NR,NF,NPHI)
      COMPLEX(8)            :: PHIR_DOT(NR,NF,NPHI)
      COMPLEX(8)            :: PHIWORK2D(NR,NF,NF)
      INTEGER(4)            :: IF,IF1,IF2,IR
      COMPLEX(8)            :: HA(NF-NPHI,NF-NPHI),HX(NF-NPHI,NPHI),HB(NF-NPHI,NPHI)
      COMPLEX(8)            :: ALLKINK_HOM(NF,NF)
      COMPLEX(8)            :: MAT(NF,NF)
      INTEGER(4)            :: IRC
      COMPLEX(8)            :: KINK_HOM(NPHI,NPHI)
      COMPLEX(8)            :: KINK_DOT(NPHI,NPHI)
      COMPLEX(8)            :: KINKC(NPHI,NPHI)
      COMPLEX(8)            :: HAM(NPHI,NPHI),OV(NPHI,NPHI)
      REAL(8)               :: R(NR)
      COMPLEX(8)            :: DHOM(NR,NF)
      LOGICAL   ,PARAMETER  :: TWRITE=.FALSE.
      COMPLEX(8)            :: BVECS(NF,NF),XVECS(NF,NF)
      COMPLEX(8)            :: CAUX(NR)
      COMPLEX(8)            :: CSVAR
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR,SVAR1,SVAR2
      INTEGER(4)            :: I,J
      INTEGER(4)            :: L0
      COMPLEX(8)            :: AMAT(NPHI-1,NPHI-1),BVEC(NPHI-1)
CHARACTER(32):: FILE
!     ******************************************************************
      TOK=.FALSE.
!PRINT*,'NEW RADIAL_XXXC STARTED',NPHI,NF
      IF(IROUT+1.GT.NR) THEN
        CALL ERROR$MSG('IROUT OUT OF RANGE')
        CALL ERROR$STOP('RADIAL_XXXC')
      END IF
      IF(IRMATCH.GT.IROUT) THEN
        CALL ERROR$MSG('IRMATCH OUT OF RANGE')
        CALL ERROR$STOP('RADIAL_XXXC')
      END IF
      L0=(NPHI/2-1)/2
      IF(NPHI.NE.2*(2*L0+1)) THEN
        CALL ERROR$MSG('NPHI MUST BE (2*L0+1)*2')
        CALL ERROR$STOP('RADIAL_XXXC')
      END IF
      CALL RADIAL$R(GID,NR,R)
      IRC=IRMATCH
!PRINT*,'IROUT ',IROUT,R(IROUT),IRC,R(IRC)
!
!     ==================================================================
!     ==  OBTAIN HOMOGENEOUS SOLUTION                                 ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC: DETERMINE PHI',L0
      ALLPHIL(:,:,:)=(0.D0,0.D0)
      ALLPHIR(:,:,:)=(0.D0,0.D0)
      DO IF=1,NF
        ALLPHIL(1:2,IF,IF)=CMPLX(R(1:2)**LOX(IF),0.D0)
        DHOM(:,:)=CMPLX(0.D0,0.D0)
        CALL RADIAL$DGLGENC(GID,NR,NF,1,IRC+1,A,B,C,DHOM,ALLPHIL(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIL(IRC,:,IF)))
        ALLPHIL(1:IRC+1,:,IF)=ALLPHIL(1:IRC+1,:,IF)/SVAR
!
        DHOM(IROUT,IF)=CMPLX(1.D-8,0.D0)
        CALL RADIAL$DGLGENC(GID,NR,NF,IROUT+1,IRC-1,A,B,C,DHOM,ALLPHIR(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIR(IRC,:,IF)))
        ALLPHIR(IRC-1:IROUT+1,:,IF)=ALLPHIR(IRC-1:IROUT+1,:,IF)/SVAR
      ENDDO
!
!     ==================================================================
!     ==  MAKE PHI_HOM SOLUTIONS CONTINUOUS                           ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC: MATCH PHI'
      MAT(:,:)=ALLPHIR(IRC,:,:)
!     == MATRIX IS NOT SYMMETRIC. THUS SOLVE WITH SINGULAR VALUE DECOMPOSITION
      BVECS(:,1:NF)=ALLPHIL(IRC,:,:)
      CALL LIB$MATRIXSOLVENEWC8(NF,NF,NF,MAT,XVECS,BVECS)
      PHIWORK2D(:,:,:)=ALLPHIR(:,:,:)
      ALLPHIR(:,:,:)=(0.D0,0.D0)
      DO IF1=1,NF
        DO IF2=1,NF
          ALLPHIR(:,:,IF1)=ALLPHIR(:,:,IF1)+PHIWORK2D(:,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTA CHANNELS        ==
!     ==================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHI'
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).EQ.L0) CYCLE
        I=I+1
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          CSVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HA(I,J)=CSVAR
        ENDDO
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).NE.L0) CYCLE
          J=J+1
          CSVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HB(I,J)=-CSVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVENEWC8(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      PHIL(:,:,:)=(0.D0,0.D0)
      PHIR(:,:,:)=(0.D0,0.D0)
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).NE.L0) CYCLE
        I=I+1
        PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF1)
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                         +ALLPHIR(IRC-1:IROUT+1,:,IF1)
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ====================================================================
!     ==   ORTHOGONALIZE PHI                                                ==
!     ====================================================================
!PRINT*,'RADIAL_XXXC: ORTHOGONALIZE PHI'
      DO I=1,NPHI
!       == ORTHOGONALIZE
        DO J=1,I-1
          CAUX(:)=(0.D0,0.D0)
          DO IF=1,NF
            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL(:IRC,IF,J))*PHIL(:IRC,IF,I)
            CAUX(IRC+1:IROUT+1)=CAUX(IRC+1:IROUT+1) &
     &             +CONJG(PHIR(IRC+1:IROUT+1,IF,J))*PHIR(IRC+1:IROUT+1,IF,I)
          ENDDO
          CAUX(:)=CAUX(:)*R(:)**2
          AUX(:)=REAL(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          AUX(:)=AIMAG(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
          CSVAR=CMPLX(SVAR1,SVAR2)
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*CSVAR
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*CSVAR
        ENDDO
!       ==  NORMALIZE
        AUX(:)=0.D0
        DO IF=1,NF
          AUX(:IRC)=AUX(:IRC)+ABS(PHIL(:IRC,IF,I))**2
          AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &                      +ABS(PHIR(IRC+1:IROUT+1,IF,I))**2
        ENDDO
        AUX(:)=AUX(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
        SVAR1=1.D0/SQRT(SVAR1)
        PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)*SVAR1
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)*SVAR1
      ENDDO
!
!     ====================================================================
!     ==   TEST CONTINUITY                                              ==
!     ====================================================================
!!$PRINT*,'RADIAL_XXXC: TEST CONTINUITY'
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          SVAR1=ABS(PHIL(IRC,IF2,IF1)-PHIR(IRC,IF2,IF1))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$PRINT*,'DEVIATION IN VALUE',SVAR2
!!$      IF(SVAR2.GT.1.D-6) THEN
!!$        CALL ERROR$MSG('MATCHING OF PHI FAILED')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF   
!!$!
!!$!     == TEST DIFFERENTIABILITY ===================================
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          IF(LOX(IF2).EQ.L0) CYCLE
!!$          SVAR1=ABS( (PHIR(IRC+1,IF2,IF1)-PHIR(IRC-1,IF2,IF1)) &
!!$     &              -(PHIL(IRC+1,IF2,IF1)-PHIL(IRC-1,IF2,IF1)))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$PRINT*,'DEVIATION IN OTHER DERIVATIVES',SVAR2
!!$      IF(SVAR2.GT.1.D-6) THEN
!!$        CALL ERROR$MSG('MATCHING OF DPHI/DR FAILED')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF   
!!$!     ==  TEST OVERLAP  =============================================
!!$      PRINT*,'OVERLAP'
!!$      DO I=1,NPHI
!!$        DO J=I,NPHI
!!$          CAUX(:)=(0.D0,0.D0)
!!$          DO IF=1,NF
!!$            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL(:IRC,IF,I))*PHIL(:IRC,IF,J)
!!$            CAUX(IRC+1:)=CAUX(IRC+1:)+CONJG(PHIR(IRC+1:,IF,I))*PHIR(IRC+1:,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX(:)=REAL(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          AUX(:)=AIMAG(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$          OV(I,J)=CMPLX(SVAR1,SVAR2)
!!$          OV(J,I)=CONJG(OV(I,J))
!!$        ENDDO
!!$      ENDDO
!!$      DO I=1,NPHI
!!$        SVAR=1.D0/SQRT(ABS(OV(I,I)))
!!$        OV(:,I)=OV(:,I)*SVAR
!!$        OV(I,:)=OV(I,:)*SVAR
!!$      ENDDO
!!$PRINT*,'MARKE A',NPHI
!!$      DO I=1,NPHI
!!$        WRITE(*,FMT='(20("(",F10.3,",",F10.3,")"))')OV(I,:)
!!$      ENDDO
!!$PRINT*,'MARKE B',NPHI
!!$!
!!$!     ==  PLOT WAVE FUNCTIONS ======================================
!!$      DO IF1=1,NPHI
!!$        WRITE(FILE,*)IF1
!!$        FILE='PHI'//ADJUSTL(FILE)
!!$        OPEN(8,FILE=FILE,FORM='FORMATTED')
!!$        REWIND 8
!!$        DO IR=1,IRC+1
!!$          WRITE(8,FMT='(80F12.7)')R(IR),REAL(PHIL(IR,:,IF1)),AIMAG(PHIL(IR,:,IF1))
!!$        ENDDO
!!$        DO IR=IRC-1,MIN(IROUT+1,NR)
!!$          WRITE(8,FMT='(80F12.7)')R(IR),REAL(PHIR(IR,:,IF1)),AIMAG(PHIR(IR,:,IF1))
!!$        ENDDO
!!$        CLOSE(8)
!!$      ENDDO
!
!     ==================================================================
!     ==  DETERMINE PHIDOT                                            ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC: DETERMINE PHIDOT'
      PHIL_DOT(:,:,:)=0.D0
      PHIR_DOT(:,:,:)=0.D0
      DO IF=1,NPHI
        CALL RADIAL$DGLGENC(GID,NR,NF,1,IRC+1,A,B,C,-2.D0*PHIL(:,:,IF) &
     &                     ,PHIL_DOT(:,:,IF))
        CALL RADIAL$DGLGENC(GID,NR,NF,IROUT+1,IRC-1,A,B,C,-2.D0*PHIR(:,:,IF) &
     &                     ,PHIR_DOT(:,:,IF))
      ENDDO
!
!     ==================================================================
!     ==  MAKE PHI_DOT CONTINUOUS                                     ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC: MATCH PHIDOT'
      MAT(:,:)=ALLPHIR(IRC,:,:)
      BVECS(:,:NPHI)=PHIL_DOT(IRC,:,:)-PHIR_DOT(IRC,:,:)
      CALL LIB$MATRIXSOLVENEWC8(NF,NF,NPHI,MAT,XVECS(:,:NPHI),BVECS(:,:NPHI))
!
      DO IF1=1,NPHI
        DO IF2=1,NF
          PHIR_DOT(IRC-1:IROUT+1,:,IF1)=PHIR_DOT(IRC-1:IROUT+1,:,IF1) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTUM CHANNELS        ==
!     ==================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHIDOT'
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).EQ.L0) CYCLE
        I=I+1
        DO J=1,NPHI
          CSVAR=(PHIR_DOT(IRC+1,IF1,J)-PHIR_DOT(IRC-1,IF1,J)) &
     &         -(PHIL_DOT(IRC+1,IF1,J)-PHIL_DOT(IRC-1,IF1,J))
          HB(I,J)=-CSVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVENEWC8(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      DO I=1,NPHI
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I) &
     &                               +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ====================================================================
!     ==   ORTHOGONALIZE PHIDOT TO PHI                                  ==
!     ====================================================================
      DO I=1,NPHI
        DO J=1,NPHI
          CAUX(:)=0.D0
          DO IF=1,NF
            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL(:IRC,IF,J))*PHIL_DOT(:IRC,IF,I)
            CAUX(IRC+1:IROUT+1)=CAUX(IRC+1:IROUT+1) &
     &             +CONJG(PHIR(IRC+1:IROUT+1,IF,J))*PHIR_DOT(IRC+1:IROUT+1,IF,I)
          ENDDO
          CAUX(:)=CAUX(:)*R(:)**2
          AUX(:)=REAL(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          AUX(:)=AIMAG(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
          CSVAR=CMPLX(SVAR1,SVAR2)
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*CSVAR
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*CSVAR
        ENDDO
      ENDDO
!
!     ====================================================================
!     ==   TEST CONTINUITY                                              ==
!     ====================================================================
!!$PRINT*,'RADIAL_XXXC: TEST CONTINUITY OF PHIDOT'
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          SVAR1=ABS(PHIL_DOT(IRC,IF2,IF1)-PHIR_DOT(IRC,IF2,IF1))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$PRINT*,'DEVIATION IN VALUE OF PHIDOT ',SVAR2
!!$      IF(SVAR2.GT.1.D-6) THEN
!!$        CALL ERROR$MSG('MATCHING FAILED')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF   
!!$!
!!$!     == TEST DIFFERENTIABILITY ===================================
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          IF(LOX(IF2).EQ.L0) CYCLE
!!$          SVAR1=ABS( (PHIR_DOT(IRC+1,IF2,IF1)-PHIR_DOT(IRC-1,IF2,IF1)) &
!!$     &              -(PHIL_DOT(IRC+1,IF2,IF1)-PHIL_DOT(IRC-1,IF2,IF1)))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$PRINT*,'DEVIATION IN OTHER DERIVATIVES OF PHIDOT',SVAR2
!!$      IF(SVAR2.GT.1.D-4) THEN
!!$        CALL ERROR$MSG('MATCHING OF DPHIDOT/DR FAILED')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF   
!!$!
!!$!     =============================================================
!!$      DO IF1=1,NPHI
!!$        WRITE(FILE,*)IF1
!!$        FILE='PHIDOT'//ADJUSTL(FILE)
!!$        OPEN(8,FILE=FILE,FORM='FORMATTED')
!!$        REWIND 8
!!$        DO IR=1,IRC+1
!!$          WRITE(8,FMT='(80F12.7)')R(IR),REAL(PHIL_DOT(IR,:,IF1)),AIMAG(PHIL_DOT(IR,:,IF1))
!!$        ENDDO
!!$        DO IR=IRC-1,MIN(IROUT+1,NR)
!!$          WRITE(8,FMT='(80F12.7)')R(IR),REAL(PHIR_DOT(IR,:,IF1)),AIMAG(PHIR_DOT(IR,:,IF1))
!!$        ENDDO
!!$        CLOSE(8)
!!$      ENDDO
!
!     ==================================================================
!     ==  REMOVE KINKS BY MIXING PHIDOT INTO PHI                      ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC:MATCH KINKS'
      I=0
      DO IF2=1,NF
        IF(LOX(IF2).NE.L0) CYCLE
        I=I+1
        KINK_HOM(I,:)=(PHIR(IRC+1,IF2,:)-PHIR(IRC-1,IF2,:)) &
     &               -(PHIL(IRC+1,IF2,:)-PHIL(IRC-1,IF2,:))
        KINK_DOT(I,:)=(PHIR_DOT(IRC+1,IF2,:)-PHIR_DOT(IRC-1,IF2,:)) &
     &               -(PHIL_DOT(IRC+1,IF2,:)-PHIL_DOT(IRC-1,IF2,:))
      ENDDO
      CALL LIB$MATRIXSOLVENEWC8(NPHI,NPHI,NPHI,-KINK_DOT,HAM,KINK_HOM)
!!$      SVAR=MAXVAL(ABS(KINK_HOM+MATMUL(KINK_DOT,HAM)))
!!$PRINT*,'REMAINING KINKS ',SVAR
!     == MAKE PHI DIFFERENTIABLE =========================================
      DO I=1,NPHI
        DO J=1,NPHI
          PHIL(:IRC,:,I)         =PHIL(:IRC,:,I)         +PHIL_DOT(:IRC,:,J)  *HAM(J,I)
          PHIR(IRC+1:IROUT+1,:,I)=PHIR(IRC+1:IROUT+1,:,I)+PHIR_DOT(IRC+1:IROUT+1,:,J)*HAM(J,I)
        ENDDO
      ENDDO
!     == MAKE PHI DIFFERENTIABLE =========================================
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("HAM ",50("(",2F10.5")   "))')HAM(I,:)
!!$ENDDO
!SVAR=0.5D0*MAXVAL(ABS(HAM-TRANSPOSE(CONJG(HAM))))
!PRINT*,'DEVIATION FROM HERMITICITY',SVAR
!
!     ==================================================================
!     ==  DETERMINE OVERLAP MATRIX                                    ==
!     ==================================================================
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("H ",50("(",2F10.5")   "))')HAM(I,:)
!!$ENDDO
!!$
!!$      DO I=1,NPHI
!!$        DO J=I,NPHI
!!$          CAUX(:)=0.D0
!!$          DO IF=1,NF
!!$            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL_DOT(:IRC,IF,I))*PHIL_DOT(:IRC,IF,J)
!!$            CAUX(IRC+1:IROUT+1)=CAUX(IRC+1:IROUT+1) &
!!$     &             +CONJG(PHIR_DOT(IRC+1:IROUT+1,IF,I))*PHIR_DOT(IRC+1:IROUT+1,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX(:)=REAL(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          AUX(:)=AIMAG(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$          OV(I,J)=CMPLX(SVAR1,SVAR2)
!!$          OV(J,I)=CONJG(OV(I,J))
!!$        ENDDO
!!$      ENDDO
!!$      OV=MATMUL(TRANSPOSE(CONJG(HAM)),MATMUL(OV,HAM))
!!$      DO I=1,NPHI
!!$        OV(I,I)=OV(I,I)+(1.D0,0.D0)
!!$      ENDDO
!!$!
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("O ",50("(",2F10.5")   "))')OV(I,:)
!!$ENDDO
!
      DO I=1,NPHI
        DO J=I,NPHI
          CAUX(:)=0.D0
          DO IF=1,NF
            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL(:IRC,IF,I))*PHIL(:IRC,IF,J)
            CAUX(IRC+1:IROUT+1)=CAUX(IRC+1:IROUT+1) &
     &             +CONJG(PHIR(IRC+1:IROUT+1,IF,I))*PHIR(IRC+1:IROUT+1,IF,J)
          ENDDO
          CAUX(:)=CAUX(:)*R(:)**2
          AUX(:)=REAL(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
          AUX(:)=AIMAG(CAUX)
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
          OV(I,J)=CMPLX(SVAR1,SVAR2)
          OV(J,I)=CONJG(OV(I,J))
        ENDDO
      ENDDO
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("O ",50("(",2F10.5")   "))')OV(I,:)
!!$ENDDO
!
!     ==================================================================
!     ==  DETERMINE EIGENSTATES                                       ==
!     ==================================================================
      HAM=0.5D0*(HAM+TRANSPOSE(CONJG(HAM)))
!      CALL LIB$DIAGC8(NPHI,HAM,DE,KINKC)
      CALL LIB$GENERALEIGENVALUEC8(NPHI,HAM,OV,DE,KINKC)
      PHI(:,:,:)=0.D0
      DO I=1,NPHI
        DO J=1,NPHI
          PHI(:IRC,:,I)         =PHI(:IRC,:,I)         +PHIL(:IRC,:,J)*KINKC(J,I)      
          PHI(IRC+1:IROUT+1,:,I)=PHI(IRC+1:IROUT+1,:,I)+PHIR(IRC+1:IROUT+1,:,J)*KINKC(J,I)      
        ENDDO
      ENDDO
!PRINT*,'DE',DE
!!$!
!!$!     ==================================================================
!!$!     ==  NORMALIZE SOLUTIONS                                         ==
!!$!     ==================================================================
!!$      DO I=1,NPHI
!!$        AUX(:)=0.D0
!!$        DO J=1,NF
!!$          AUX(:)=AUX(:)+ABS(PHI(:,J,I))**2
!!$        ENDDO
!!$        AUX(:)=AUX(:)*R(:)**2
!!$        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$        SVAR=1.D0/SQRT(SVAR)
!!$        PHI(:,:,I)=PHI(:,:,I)*SVAR
!!$      ENDDO
!
!     ==================================================================
!     ==  TEST OVERLAP                                                ==
!     ==================================================================
!!$      DO IF1=1,NPHI
!!$        WRITE(FILE,*)IF1
!!$        FILE='PHIFIN'//ADJUSTL(FILE)
!!$        OPEN(8,FILE=FILE,FORM='FORMATTED')
!!$        REWIND 8
!!$        DO IR=1,IROUT+1
!!$          WRITE(8,FMT='(80F20.7)')R(IR),REAL(PHI(IR,:,IF1)),AIMAG(PHI(IR,:,IF1))
!!$        ENDDO
!!$        CLOSE(8)
!!$      ENDDO
!!$!
!     == TEST ORTHONORMALITY OVERLAP MATRIX
!!$      DO I=1,NPHI
!!$        HAM(I,I)=(0.D0,0.D0)
!!$        DO J=I,NPHI
!!$          CAUX(:)=0.D0
!!$          DO IF=1,NF
!!$            CAUX(:)=CAUX(:)+CONJG(PHI(:,IF,I))*PHI(:,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX=REAL(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          AUX=AIMAG(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$          HAM(I,J)=CMPLX(SVAR1,SVAR2)
!!$          HAM(J,I)=CONJG(HAM(I,J))
!!$        ENDDO
!!$        HAM(I,I)=HAM(I,I)-(1.D0,0.D0)
!!$      ENDDO
!!$      DO I=1,NPHI
!!$        WRITE(*,FMT='("O ",50("(",2F10.5")   "))')HAM(I,:)
!!$      ENDDO
         SVAR=MAXVAL(ABS(HAM))
         PRINT*,'DEVIATION FROM ORTHONORMALITY',SVAR,NF
!!$      IF(SVAR.GT.0.5D0) THEN
!!$        CALL ERROR$MSG('WAVE FUNCTIONS NOT ORTHOGONAL')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF
      TOK=.TRUE.
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE RADIAL$NONSPHBOUND_nso_sloc_ov(GID,NR,LMX,LMRX,POT,lxsl,potsl,ovsl,G,Enu &
     &                             ,NPHI,EB,PHI,TPHI,TOK)
!     **                                                                  **
!     **  SOLVES THE nonrelATIVISTIC RADIAL schroedinger EQUATION FOR THE **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE        **
!     **  HOMOGENEOUS SOLUTION.                                           **
!     **                                                                  **
!     **  THE SOLUTIONS ARE CALCULATED TO LINEAR ORDER IN DE              **
!     **  WHERE E+DE IS THE NEW ENERGY OF THE WAVE FUNCTIONS              **
!     **                                                                  **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM   **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                     **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                               **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                         **
!     **                                                                  **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE      **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                   **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID-ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: LMX     ! X#(WAVE FUNCTION ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)     :: LMRX    ! X#(POTENTIAL ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)     :: Lxsl    ! #(semilocal potenials, one per l)
      REAL(8)    ,INTENT(IN)     :: potsl(NR,Lxsl)  ! semi-local potentials
      REAL(8)    ,INTENT(IN)     :: ovsl(NR,lxsl)   ! semi-local overlap contribtion
      REAL(8)    ,INTENT(IN)     :: G(NR,LMX)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: Enu     ! expansion ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR,LMRX) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NPHI              ! #(wave functions)
      real(8)    ,INTENT(OUT)    :: PHI(NR,LMX,nphi)  ! WAVE-FUNCTION
      real(8)    ,INTENT(OUT)    :: TPHI(NR,LMX,NPHI) ! p**2/(2m)*WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: EB(NPHI)          ! one-particke energies
      LOGICAL(4) ,INTENT(OUT)    :: TOK               ! error flag
      INTEGER(4)                 :: LX
      INTEGER(4)                 :: LOX(LMX) ! ANGULAR MOMENTA
      INTEGER(4)                 :: LM,LM1,LM2,LM3,L,M,IM,IB,ir
      REAL(8)                    :: A(NR)
      REAL(8)                    :: B(NR)
      reaL(8)                    :: C(NR,LMX,LMX)
      reaL(8)                    :: Cov(NR,LMX)
      REAL(8)                    :: D(NR,LMX) 
      REAL(8)                    :: DPHI(NR,LMX)
      REAL(8)                    :: R(NR)                        !
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: PI                           !
      REAL(8)                    :: Y0                           !
      REAL(8)                    :: CG
      REAL(8)                    :: RDPRIME(NR)                  !
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20
      REAL(8)                    :: SWKB(NR),SVAR
      INTEGER(4)                 :: IRCL,IROUT
      LOGICAL(4)                 :: TCHK
!     ************************************************************************
print*,'warning! overlap not included yet'
      TOK=.FALSE.
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      LX=INT(SQRT(REAL(LMX))-1.D0)
      IF((LX+1)**2.NE.LMX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('RADIAL$NONSPHBOUND_NONSO')
      END IF
      LM=0
      DO L=0,LX
        DO M=1,2*L+1
          LM=LM+1
          LOX(LM)=L
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DETERMINE CLASSICAL TURNING POINT                           ==
!     ==================================================================
!     == IROUT MUST BE SMALLER OR EQUAL TO NR-1 BECAUSE INHOMOGENEITY
!     == ON THE LAST GRID POINT DOES NOT CONTRIBUTE TO DIFFERENTIAL 
!     == EQUATION. IROUT IS SMALLER OR EQUAL TO IRCL
      IRCL=NR-3   
!       == R(IRCL) IS THE FIRST GRID POINT INSIDE THE CLASSICAL TURNING POINT
      IF(Enu.LT.POT(IRCL,1)*Y0) THEN
        DO IR=IRCL,1,-1
          IF(Enu.GT.POT(IR,1)*Y0) THEN
            IRCL=IR
            EXIT
          END IF
        ENDDO
      END IF
!
!     ==================================================================
!     ==  DETERMINE OUTERMOST POINT 'IROUT' FOR INWARD INTEGRATION    ==
!     ==================================================================
!     == USE WKB SOLUTION FOR THE SCHR.GL. FOR A CONSTANT POTENTIAL AND L=0
!     == TO ESTIMATE FACTOR FROM RCL TO OUTERMOST POINT
      AUX(:IRCL)=0.D0
      AUX(IRCL+1:)=SQRT(MAX(0.D0,2.D0*(POT(IRCL+1:,1)*Y0-Enu)))
      CALL RADIAL$INTEGRATE(GID,NR,AUX,SWKB)
      SWKB(:)=SWKB(:)-LOG(R(:))
!     == DETERMINE IROUT WHERE THE WAVE FUNCTION CAN GROW BY A FACTOR 
!     == OF XMAX FROM THE CLASSICAL TURNING POINT
      SVAR=LOG(XMAX)
      IROUT=NR-1
      DO IR=IRCL,NR
        IF(SWKB(IR).GT.SVAR) THEN
          IROUT=IR-1
          EXIT
        END IF
      ENDDO
!
!     ==================================================================
!     ==  PREPARE POTENTIAL-INDEPENDENT ARRAYS                        ==
!     ==================================================================
!     == A*D2F/DR2+B*DF/DR+C*F=D
      A(:)=1.D0
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0/R(2:)
      B(1)=B(2)
      D(:,:)=-2.D0*G(:,:)
!
!     ==================================================================
!     ==  COUPLING BETWEEN WAVE FUNCTION COMPONENTS VIA POTENTIAL     ==
!     ==================================================================
      C(:,:,:)=0.D0
      DO LM1=1,LMX
        DO LM2=1,LMX
          AUX(:)=0.D0
          DO LM3=1,LMRX
            CALL CLEBSCH(LM1,LM2,LM3,CG)
            IF(CG.EQ.0.D0) CYCLE
            AUX=CG*POT(:,LM3)
            IF(LM3.EQ.1) AUX(IROUT+1:)=AUX(IROUT)  ! CONSTANT POTENTIAL BEYOND IROUT
            C(:,LM1,LM2)=C(:,LM1,LM2)+AUX
          ENDDO
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  add semi-lcoal potential                                    ==
!     ==================================================================
      lm=0
      do l=1,lxsl
        do m=1,2*l+1
          lm=lm+1
          c(:,lm,lm)=c(:,lm,lm)+potsl(:,l)
          cov(:,lm)=+ovsl(:,l)
        enddo
      enddo
!
!     ==================================================================
!     ==  scale c                                                     ==
!     ==================================================================
      C=-2.D0*C
!
!     ==================================================================
!     ==  add KINETIC ENERGY TERM TO C and shift energy zero to enu   ==
!     ==================================================================
      LM=0
      DO L=0,LX
        AUX(1)=0.D0
        AUX(2:)=-1.D0/R(2:)**2 * REAL(L*(L+1),KIND=8)+2.D0*Enu
        DO IM=1,2*L+1
          LM=LM+1
          C(:,LM,LM)=C(:,LM,LM)+AUX(:)
        ENDDO
      ENDDO
!     ==  AVOID DIVIDE-BY-ZERO
      C(1,:,:)=C(2,:,:)
!
!     ==================================================================
!     ==  DETERMINE BOUND STATES                                      ==
!     ==================================================================
      CALL RADIAL_XXXR(GID,NR,LMX,IRCL,IROUT,LOX,A,B,C,D,NPHI,eb,PHI,TCHK)
      EB(:)=enu+eb(:) ! change energies relative to eny to absolute energies
      IF(.NOT.TCHK) THEN
        CALL ERROR$STOP('RADIAL_XXXR FINISHED WITH ERROR')
        CALL ERROR$STOP('RADIAL$NONSPHBOUND')
      END IF
!
!     ==================================================================
!     ==  SHIFT ENERGIES                                              ==
!     ==================================================================
!
      TOK=.TRUE.
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE RADIAL$NONSPHBOUND_nonso(GID,NR,LMX,LMRX,POT,DREL,G,Enu &
     &                             ,NPHI,EB,PHI,TPHI,TOK)
!     **                                                                  **
!     **  SOLVES THE RELATIVISTIC RADIAL DIRAC EQUATION FOR THE           **
!     **  LARGE COMPONENT. DREL=1/MREL-1/M0 IS A MEASURE FOR THE          **
!     **  RELATIVISTIC EFFECTS, WHERE MREL=M0+(E-V)/2C^2 AND M0 IS THE    **
!     **  REST MASS. V=POT*Y0 IS THE POTENTIAL.                           **
!     **  SPIN ORBIT COUPLING IS MEASURED BY SO WHICH CAN HAVE THE VALUES:**
!     **    SO=0 NO-SPIN ORBIT COUPLING                                   **
!     **    SO=1 PARALLEL SPIN AND ORBITAL ANGULAR MOMENTUM               **
!     **  G IS AN INHOMOGENEITY, WHICH MUST BE SET TO ZERO FOR THE        **
!     **  HOMOGENEOUS SOLUTION.                                           **
!     **                                                                  **
!     **  THE SOLUTIONS ARE CALCULATED TO LINEAR ORDER IN DE              **
!     **  WHERE E+DE IS THE NEW ENERGY OF THE WAVE FUNCTIONS              **
!     **                                                                  **
!     **  THE DIFFERENTIAL EQUATION IS SOLVED WITH THE VERLET ALGORITHM   **
!     **    DPHI/DX=(PHI(+)-PHI(-))/2                                     **
!     **    D2PHI/DX2=PHI(+)-2PHI(0)+PHI(-)                               **
!     **  WHERE X IS THE VARIABLE WITH R(X=I)=R_I                         **
!     **                                                                  **
!     **  THE NONRELATIVISTIC SOLUTION IS OBTAINED BY SETTINH DREL=0      **
!     **                                                                  **
!     **  ATTENTION! THE ROUTINE IS NOT GUARDED AGAINST OVERFLOW DUE      **
!     **    TO THE EXPONENTIAL INCREASE OF THE SOLUTION                   **
!     **                                                                  **
!     *********************** COPYRIGHT: PETER BLOECHL, GOSLAR 2006 ********
      IMPLICIT NONE
      INTEGER(4) ,INTENT(IN)     :: GID     ! GRID-ID FOR RADIAL GRID
      INTEGER(4) ,INTENT(IN)     :: NR      ! #(RADIAL GRID POINTS)
      INTEGER(4) ,INTENT(IN)     :: LMX     ! X#(WAVE FUNCTION ANGULAR MOMENTA)
      INTEGER(4) ,INTENT(IN)     :: LMRX    ! X#(POTENTIAL ANGULAR MOMENTA)
      REAL(8)    ,INTENT(IN)     :: DREL(NR)!RELATIVISTIC CORRECTION
!                                           ! DREL= M0/MREL(R)-1
      REAL(8)    ,INTENT(IN)     :: G(NR,LMX)   !INHOMOGENITY
      REAL(8)    ,INTENT(IN)     :: Enu     ! expansion ENERGY
      REAL(8)    ,INTENT(IN)     :: POT(NR,LMRX) !POTENTIAL (RADIAL PART ONLY)
      INTEGER(4) ,INTENT(IN)     :: NPHI              ! #(wave functions)
      real(8)    ,INTENT(OUT)    :: PHI(NR,LMX,nphi)  ! WAVE-FUNCTION
      real(8)    ,INTENT(OUT)    :: TPHI(NR,LMX,NPHI) ! p**2/(2m)*WAVE-FUNCTION
      REAL(8)    ,INTENT(OUT)    :: EB(NPHI)          ! one-particke energies
      LOGICAL(4) ,INTENT(OUT)    :: TOK               ! error flag
      INTEGER(4)                 :: LX
      INTEGER(4)                 :: LOX(LMX) ! ANGULAR MOMENTA
      INTEGER(4)                 :: LM,LM1,LM2,LM3,L,M,IM,IB,ir
      REAL(8)                    :: A(NR)
      REAL(8)                    :: B(NR)
      reaL(8)                    :: C(NR,LMX,LMX)
      REAL(8)                    :: D(NR,LMX) 
      REAL(8)                    :: DPHI(NR,LMX)
      REAL(8)                    :: R(NR)                        !
      REAL(8)                    :: AUX(NR)
      REAL(8)                    :: PI                           !
      REAL(8)                    :: Y0                           !
      REAL(8)                    :: CG
      REAL(8)                    :: RDPRIME(NR)                  !
      REAL(8)   ,PARAMETER       :: XMAX=1.D+20
      REAL(8)                    :: SWKB(NR),SVAR
      INTEGER(4)                 :: IRCL,IROUT
      LOGICAL(4)                 :: TCHK
!     ************************************************************************
      TOK=.FALSE.
      PI=4.D0*ATAN(1.D0)
      Y0=1.D0/SQRT(4.D0*PI)
      CALL RADIAL$R(GID,NR,R)
      CALL RADIAL$DERIVE(GID,NR,DREL,RDPRIME) 
      LX=INT(SQRT(REAL(LMX))-1.D0)
      IF((LX+1)**2.NE.LMX) THEN   !LMX=(LX+1)**2
        CALL ERROR$MSG('LMX DOES NOT CORRESPOND TO A FULL SHELL')
        CALL ERROR$MSG('OR ROUNDING ERRORS PRODUCED INCORRECT RESULTS')
        CALL ERROR$I4VAL('LMX',LMX)
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('RADIAL$NONSPHBOUND_NONSO')
      END IF
      LM=0
      DO L=0,LX
        DO M=1,2*L+1
          LM=LM+1
          LOX(LM)=L
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  DETERMINE CLASSICAL TURNING POINT r(ircl)                   ==
!     ==   and outmost point for inward integration r(irout)          ==
!     ==================================================================
      call schroedinger_specialrads(gid,nr,0,xmax,pot(:,1),enu,ircl,irout)
!
!     ==================================================================
!     ==  PREPARE POTENTIAL-INDEPENDENT ARRAYS                        ==
!     ==================================================================
!     == A*D2F/DR2+B*DF/DR+C*F=D
      A(:)=1.D0+DREL(:)
!     == AVOID DIVIDE BY ZERO IF THE FIRST GRID POINT IS THE ORIGIN.
!     == THE FORCES ON THE FIRST GRID POINT ARE NOT USED,
!     == BECAUSE RADIAL$DGL IS BASED ON THE VERLET ALGORITHM
!     == THAT CANNOT USE THE FORCES ON THE FIRST AND LAST GRID POINT
      B(2:)=2.D0*(1.D0+DREL(2:))/R(2:)+RDPRIME(2:)
      B(1)=B(2)
      D(:,:)=-2.D0*G(:,:)
!
!     ==================================================================
!     ==  COUPLING BETWEEN WAVE FUNCTION COMPONENTS VIA POTENTIAL     ==
!     ==================================================================
      C(:,:,:)=0.D0
      DO LM1=1,LMX
        DO LM2=1,LMX
          AUX(:)=0.D0
          DO LM3=1,LMRX
            CALL CLEBSCH(LM1,LM2,LM3,CG)
            IF(CG.EQ.0.D0) CYCLE
            AUX=CG*POT(:,LM3)
            IF(LM3.EQ.1) AUX(IROUT+1:)=AUX(IROUT)  ! CONSTANT POTENTIAL BEYOND IROUT
            C(:,LM1,LM2)=C(:,LM1,LM2)+AUX
          ENDDO
        ENDDO
      ENDDO
      C=-2.D0*C
!
!     ==================================================================
!     ==  add KINETIC ENERGY TERM TO C and shift energy zero to enu   ==
!     ==================================================================
      LM=0
      DO L=0,LX
        AUX(1)=0.D0
        AUX(2:)=-(1.D0+DREL(2:))/R(2:)**2 * REAL(L*(L+1),KIND=8)+2.D0*Enu
        DO IM=1,2*L+1
          LM=LM+1
          C(:,LM,LM)=C(:,LM,LM)+AUX(:)
        ENDDO
      ENDDO
!     ==  AVOID DIVIDE-BY-ZERO
      C(1,:,:)=C(2,:,:)
!
!     ==================================================================
!     ==  DETERMINE BOUND STATES                                      ==
!     ==================================================================
      CALL RADIAL_XXXR(GID,NR,LMX,IRCL,IROUT,LOX,A,B,C,D,NPHI,eb,PHI,TCHK)
      IF(.NOT.TCHK) THEN
        CALL ERROR$STOP('RADIAL_XXXR FINISHED WITH ERROR')
        CALL ERROR$STOP('RADIAL$NONSPHBOUND')
      END IF
!
!     ==================================================================
!     ==  SHIFT ENERGIES                                              ==
!     ==================================================================
      EB(:)=enu+eb(:) ! change energies relative to eny to absolute energies
!
      TOK=.TRUE.
      RETURN
      END 
!
!     ..................................................................
      SUBROUTINE RADIAL_XXXR(GID,NR,NF,IRMATCH,IROUT,LOX,A,B,C,D,NPHI,DE,PHI,TOK)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID             ! GRID-ID FOR RADIAL GRID
      INTEGER(4),INTENT(IN) :: NR              ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: NF              ! #(angular momenta)
      INTEGER(4),INTENT(IN) :: IRMATCH         ! matching point for inside-outside integration
      INTEGER(4),INTENT(IN) :: IROUT           ! outermost point to be considered
      INTEGER(4),INTENT(IN) :: LOX(NF)         ! angular momenta of wave fyunction components
      REAL(8)   ,INTENT(IN) :: A(NR)
      REAL(8)   ,INTENT(IN) :: B(NR)
      reAL(8)   ,INTENT(IN) :: C(NR,NF,NF)
      REAl(8)   ,INTENT(IN) :: D(NR,NF)
      INTEGER(4),INTENT(IN) :: NPHI            ! #(wave functions)
      REAL(8)   ,INTENT(OUT):: DE(NPHI)        ! one-particle eigenvalues
      real(8)   ,INTENT(OUT):: PHI(NR,NF,NPHI) ! wave functions
      LOGICAL(4),INTENT(OUT):: TOK             ! error flag
      REAL(8)               :: ALLPHIL(NR,NF,NF)
      REAL(8)               :: ALLPHIR(NR,NF,NF)
      REAL(8)               :: PHIL(NR,NF,NPHI)
      REAL(8)               :: PHIR(NR,NF,NPHI)
      REAL(8)               :: PHIL_DOT(NR,NF,NPHI)
      REAL(8)               :: PHIR_DOT(NR,NF,NPHI)
      REAL(8)               :: PHIWORK2D(NR,NF,NF)
      INTEGER(4)            :: IF,IF1,IF2,IR
      REAL(8)               :: HA(NF-NPHI,NF-NPHI),HX(NF-NPHI,NPHI),HB(NF-NPHI,NPHI)
      REAL(8)               :: ALLKINK_HOM(NF,NF)
      REAL(8)               :: MAT(NF,NF)

      REAL(8)               :: DE1(NPHI)
      INTEGER(4)            :: IRC
      REAL(8)               :: KINK_HOM(NPHI,NPHI)
      REAL(8)               :: KINK_DOT(NPHI,NPHI)
      REAL(8)               :: KINKC(NPHI,NPHI)
      REAL(8)               :: HAM(NPHI,NPHI),OV(NPHI,NPHI)
      REAL(8)               :: R(NR)
      REAL(8)               :: DHOM(NR,NF)
      LOGICAL   ,PARAMETER  :: TWRITE=.FALSE.
      REAL(8)               :: BVECS(NF,NF),XVECS(NF,NF)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR
      INTEGER(4)            :: I,J
      INTEGER(4)            :: L0
      REAL(8)               :: AMAT(NPHI-1,NPHI-1),BVEC(NPHI-1)
CHARACTER(32):: FILE
!     ******************************************************************
      TOK=.FALSE.
!PRINT*,'NEW RADIAL_XXXR STARTED',NPHI,NF
      IF(IROUT+1.GT.NR) THEN
        CALL ERROR$MSG('IROUT OUT OF RANGE')
        CALL ERROR$STOP('RADIAL_XXXC')
      END IF
      IF(IRMATCH.GT.IROUT) THEN
        CALL ERROR$MSG('IRMATCH OUT OF RANGE')
        CALL ERROR$STOP('RADIAL_XXXC')
      END IF
      L0=(NPHI-1)/2
      IF(NPHI.NE.2*L0+1) THEN
        CALL ERROR$MSG('NPHI MUST BE 2*L0+1')
        CALL ERROR$i4val('NPHI',nphi)
        CALL ERROR$i4val('l0',l0)
        CALL ERROR$STOP('RADIAL_XXXC')
      END IF
      CALL RADIAL$R(GID,NR,R)
      IRC=IRMATCH
!
!     ==================================================================
!     ==  OBTAIN HOMOGENEOUS SOLUTION                                 ==
!     ==================================================================
!PRINT*,'RADIAL_XXXr: DETERMINE PHI',Lox
      ALLPHIL(:,:,:)=0.D0
      ALLPHIR(:,:,:)=0.D0
      DO IF=1,NF
        ALLPHIL(1:2,IF,IF)=R(1:2)**LOX(IF)
        DHOM(:,:)=0.D0
        CALL RADIAL$DGLGEN(GID,NR,NF,1,IRC+1,A,B,C,DHOM,ALLPHIL(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIL(IRC,:,IF)))
        ALLPHIL(1:IRC+1,:,IF)=ALLPHIL(1:IRC+1,:,IF)/SVAR
        DHOM(IROUT,IF)=1.D-8
        CALL RADIAL$DGLGEN(GID,NR,NF,IROUT+1,IRC-1,A,B,C,DHOM,ALLPHIR(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIR(IRC,:,IF)))
        ALLPHIR(IRC-1:IROUT+1,:,IF)=ALLPHIR(IRC-1:IROUT+1,:,IF)/SVAR
      ENDDO
!
!     ==================================================================
!     ==  MAKE PHI_HOM SOLUTIONS CONTINUOUS                           ==
!     ==================================================================
!PRINT*,'RADIAL_XXXr: MATCH PHI'
      MAT(:,:)=ALLPHIR(IRC,:,:)
!     == MATRIX IS NOT SYMMETRIC. THUS SOLVE WITH SINGULAR VALUE DECOMPOSITION
      BVECS(:,1:NF)=ALLPHIL(IRC,:,:)
      CALL LIB$MATRIXSOLVENEW(NF,NF,NF,MAT,XVECS,BVECS)

      PHIWORK2D(:,:,:)=ALLPHIR(:,:,:)
      ALLPHIR(:,:,:)=0.D0
      DO IF1=1,NF
        DO IF2=1,NF
          ALLPHIR(:,:,IF1)=ALLPHIR(:,:,IF1)+PHIWORK2D(:,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTA CHANNELS        ==
!     ==================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHI'
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).EQ.L0) CYCLE
        I=I+1
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          SVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HA(I,J)=SVAR
        ENDDO
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).NE.L0) CYCLE
          J=J+1
          SVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HB(I,J)=-SVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVENEW(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      PHIL(:,:,:)=0.D0
      PHIR(:,:,:)=0.D0
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).NE.L0) CYCLE
        I=I+1
        PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF1)
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                         +ALLPHIR(IRC-1:IROUT+1,:,IF1)
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ====================================================================
!     ==   ORTHOGONALIZE PHI                                                ==
!     ====================================================================
!PRINT*,'RADIAL_XXXC: ORTHOGONALIZE PHI'
      DO I=1,NPHI
!       == ORTHOGONALIZE
        DO J=1,I-1
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)=AUX(:IRC)+PHIL(:IRC,IF,J)*PHIL(:IRC,IF,I)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &                   +PHIR(IRC+1:IROUT+1,IF,J)*PHIR(IRC+1:IROUT+1,IF,I)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*SVAR
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*SVAR
        ENDDO
!       ==  NORMALIZE
        AUX(:)=0.D0
        DO IF=1,NF
          AUX(:IRC)         =AUX(:IRC)         +PHIL(:IRC,IF,I)**2
          AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1)+PHIR(IRC+1:IROUT+1,IF,I)**2
        ENDDO
        AUX(:)=AUX(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        SVAR=1.D0/SQRT(SVAR)
        PHIL(:IRC+1,:,I)       =PHIL(:IRC+1,:,I)*SVAR
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)*SVAR
      ENDDO
!
!     ==================================================================
!     ==  DETERMINE PHIDOT                                            ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC: DETERMINE PHIDOT'
      PHIL_DOT(:,:,:)=0.D0
      PHIR_DOT(:,:,:)=0.D0
      DO IF=1,NPHI
        CALL RADIAL$DGLGEN(GID,NR,NF,1,IRC+1,A,B,C,-2.D0*PHIL(:,:,IF) &
     &                     ,PHIL_DOT(:,:,IF))
        CALL RADIAL$DGLGEN(GID,NR,NF,IROUT+1,IRC-1,A,B,C,-2.D0*PHIR(:,:,IF) &
     &                     ,PHIR_DOT(:,:,IF))
      ENDDO
!
!     ==================================================================
!     ==  MAKE PHI_DOT CONTINUOUS                                     ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC: MATCH PHIDOT'
      MAT(:,:)=ALLPHIR(IRC,:,:)
      BVECS(:,:NPHI)=PHIL_DOT(IRC,:,:)-PHIR_DOT(IRC,:,:)
      CALL LIB$MATRIXSOLVENEW(NF,NF,NPHI,MAT,XVECS(:,:NPHI),BVECS(:,:NPHI))
!
      DO IF1=1,NPHI
        DO IF2=1,NF
          PHIR_DOT(IRC-1:IROUT+1,:,IF1)=PHIR_DOT(IRC-1:IROUT+1,:,IF1) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTUM CHANNELS        ==
!     ==================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHIDOT'
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).EQ.L0) CYCLE
        I=I+1
        DO J=1,NPHI
          SVAR=(PHIR_DOT(IRC+1,IF1,J)-PHIR_DOT(IRC-1,IF1,J)) &
     &         -(PHIL_DOT(IRC+1,IF1,J)-PHIL_DOT(IRC-1,IF1,J))
          HB(I,J)=-SVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVENEW(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      DO I=1,NPHI
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I) &
     &                               +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ====================================================================
!     ==   ORTHOGONALIZE PHIDOT TO PHI                                  ==
!     ====================================================================
      DO I=1,NPHI
        DO J=1,NPHI
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)=AUX(:IRC)+PHIL(:IRC,IF,J)*PHIL_DOT(:IRC,IF,I)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &             +PHIR(IRC+1:IROUT+1,IF,J)*PHIR_DOT(IRC+1:IROUT+1,IF,I)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*SVAR
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*SVAR
        ENDDO
      ENDDO
!
!     ==================================================================
!     ==  REMOVE KINKS BY MIXING PHIDOT INTO PHI                      ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC:MATCH KINKS'
      I=0
      DO IF2=1,NF
        IF(LOX(IF2).NE.L0) CYCLE
        I=I+1
        KINK_HOM(I,:)=(PHIR(IRC+1,IF2,:)-PHIR(IRC-1,IF2,:)) &
     &               -(PHIL(IRC+1,IF2,:)-PHIL(IRC-1,IF2,:))
        KINK_DOT(I,:)=(PHIR_DOT(IRC+1,IF2,:)-PHIR_DOT(IRC-1,IF2,:)) &
     &               -(PHIL_DOT(IRC+1,IF2,:)-PHIL_DOT(IRC-1,IF2,:))
      ENDDO
      CALL LIB$MATRIXSOLVENEW(NPHI,NPHI,NPHI,-KINK_DOT,HAM,KINK_HOM)
!     == MAKE PHI DIFFERENTIABLE =========================================
      DO I=1,NPHI
        DO J=1,NPHI
          PHIL(:IRC,:,I)         =PHIL(:IRC,:,I)         +PHIL_DOT(:IRC,:,J)  *HAM(J,I)
          PHIR(IRC+1:IROUT+1,:,I)=PHIR(IRC+1:IROUT+1,:,I)+PHIR_DOT(IRC+1:IROUT+1,:,J)*HAM(J,I)
        ENDDO
      ENDDO
!     == MAKE PHI DIFFERENTIABLE =========================================
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("HAM ",50("(",2F10.5")   "))')HAM(I,:)
!!$ENDDO
!SVAR=0.5D0*MAXVAL(ABS(HAM-TRANSPOSE(CONJG(HAM))))
!PRINT*,'DEVIATION FROM HERMITICITY',SVAR
!
!     ==================================================================
!     ==  DETERMINE OVERLAP MATRIX                                    ==
!     ==================================================================
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("H ",50("(",2F10.5")   "))')HAM(I,:)
!!$ENDDO
!!$
!!$      DO I=1,NPHI
!!$        DO J=I,NPHI
!!$          CAUX(:)=0.D0
!!$          DO IF=1,NF
!!$            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL_DOT(:IRC,IF,I))*PHIL_DOT(:IRC,IF,J)
!!$            CAUX(IRC+1:IROUT+1)=CAUX(IRC+1:IROUT+1) &
!!$     &             +CONJG(PHIR_DOT(IRC+1:IROUT+1,IF,I))*PHIR_DOT(IRC+1:IROUT+1,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX(:)=REAL(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          AUX(:)=AIMAG(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$          OV(I,J)=CMPLX(SVAR1,SVAR2)
!!$          OV(J,I)=CONJG(OV(I,J))
!!$        ENDDO
!!$      ENDDO
!!$      OV=MATMUL(TRANSPOSE(CONJG(HAM)),MATMUL(OV,HAM))
!!$      DO I=1,NPHI
!!$        OV(I,I)=OV(I,I)+(1.D0,0.D0)
!!$      ENDDO
!!$!
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("O ",50("(",2F10.5")   "))')OV(I,:)
!!$ENDDO
!
      DO I=1,NPHI
        DO J=I,NPHI
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)         =AUX(:IRC)+PHIL(:IRC,IF,I)*PHIL(:IRC,IF,J)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &                         +PHIR(IRC+1:IROUT+1,IF,I)*PHIR(IRC+1:IROUT+1,IF,J)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          OV(I,J)=SVAR
          OV(J,I)=OV(I,J)
        ENDDO
      ENDDO
!call writephi(gid,nr,'homogeneous',nf,nphi,irc,phil,phir)
!
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("O ",50("(",2F10.5")   "))')OV(I,:)
!!$ENDDO
!
!     ==================================================================
!     ==  DETERMINE EIGENSTATES                                       ==
!     ==================================================================
      HAM=0.5D0*(HAM+TRANSPOSE(HAM))
!      CALL LIB$DIAGC8(NPHI,HAM,DE,KINKC)
      CALL LIB$GENERALEIGENVALUER8(NPHI,HAM,OV,DE,KINKC)
      PHI(:,:,:)=0.D0
      DO I=1,NPHI
        DO J=1,NPHI
          PHI(:IRC,:,I)         =PHI(:IRC,:,I)         +PHIL(:IRC,:,J)*KINKC(J,I)      
          PHI(IRC+1:IROUT+1,:,I)=PHI(IRC+1:IROUT+1,:,I)+PHIR(IRC+1:IROUT+1,:,J)*KINKC(J,I)      
        ENDDO
      ENDDO
!PRINT*,'DE',DE
!!$!
!!$!     ==================================================================
!!$!     ==  NORMALIZE SOLUTIONS                                         ==
!!$!     ==================================================================
!!$      DO I=1,NPHI
!!$        AUX(:)=0.D0
!!$        DO J=1,NF
!!$          AUX(:)=AUX(:)+ABS(PHI(:,J,I))**2
!!$        ENDDO
!!$        AUX(:)=AUX(:)*R(:)**2
!!$        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$        SVAR=1.D0/SQRT(SVAR)
!!$        PHI(:,:,I)=PHI(:,:,I)*SVAR
!!$      ENDDO
!
!     ==================================================================
!     ==  TEST OVERLAP                                                ==
!     ==================================================================
!!$      DO IF1=1,NPHI
!!$        WRITE(FILE,*)IF1
!!$        FILE='PHIFIN'//ADJUSTL(FILE)
!!$        OPEN(8,FILE=FILE,FORM='FORMATTED')
!!$        REWIND 8
!!$        DO IR=1,IROUT+1
!!$          WRITE(8,FMT='(80F20.7)')R(IR),REAL(PHI(IR,:,IF1)),AIMAG(PHI(IR,:,IF1))
!!$        ENDDO
!!$        CLOSE(8)
!!$      ENDDO
!!$!
!     == TEST ORTHONORMALITY OVERLAP MATRIX
!!$      DO I=1,NPHI
!!$        HAM(I,I)=(0.D0,0.D0)
!!$        DO J=I,NPHI
!!$          CAUX(:)=0.D0
!!$          DO IF=1,NF
!!$            CAUX(:)=CAUX(:)+CONJG(PHI(:,IF,I))*PHI(:,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX=REAL(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          AUX=AIMAG(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$          HAM(I,J)=CMPLX(SVAR1,SVAR2)
!!$          HAM(J,I)=CONJG(HAM(I,J))
!!$        ENDDO
!!$        HAM(I,I)=HAM(I,I)-(1.D0,0.D0)
!!$      ENDDO
!!$      DO I=1,NPHI
!!$        WRITE(*,FMT='("O ",50("(",2F10.5")   "))')HAM(I,:)
!!$      ENDDO
!!$      SVAR=MAXVAL(ABS(HAM))
!!$      PRINT*,'DEVIATION FROM ORTHONORMALITY',SVAR,NF
!!$      IF(SVAR.GT.0.5D0) THEN
!!$        CALL ERROR$MSG('WAVE FUNCTIONS NOT ORTHOGONAL')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF
      TOK=.TRUE.
call writephi(gid,nr,'final',nf,nphi,irc,phi,phi)
      RETURN
      END SUBROUTINE RADIAL_XXXR
!
!     ..................................................................
      SUBROUTINE RADIAL_XXXR_Ov(GID,NR,NF,IRMATCH,IROUT,LOX,A,B,C,D,NPHI,DE,PHI,TOK)
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID             ! GRID-ID FOR RADIAL GRID
      INTEGER(4),INTENT(IN) :: NR              ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: NF              ! #(angular momenta)
      INTEGER(4),INTENT(IN) :: IRMATCH         ! matching point for inside-outside integration
      INTEGER(4),INTENT(IN) :: IROUT           ! outermost point to be considered
      INTEGER(4),INTENT(IN) :: LOX(NF)         ! angular momenta of wave fyunction components
      REAL(8)   ,INTENT(IN) :: A(NR)
      REAL(8)   ,INTENT(IN) :: B(NR)
      reAL(8)   ,INTENT(IN) :: C(NR,NF,NF)
      REAl(8)   ,INTENT(IN) :: D(NR,NF)
      INTEGER(4),INTENT(IN) :: NPHI            ! #(wave functions)
      REAL(8)   ,INTENT(OUT):: DE(NPHI)        ! one-particle eigenvalues
      real(8)   ,INTENT(OUT):: PHI(NR,NF,NPHI) ! wave functions
      LOGICAL(4),INTENT(OUT):: TOK             ! error flag
      REAL(8)               :: ALLPHIL(NR,NF,NF)
      REAL(8)               :: ALLPHIR(NR,NF,NF)
      REAL(8)               :: PHIL(NR,NF,NPHI)
      REAL(8)               :: PHIR(NR,NF,NPHI)
      REAL(8)               :: PHIL_DOT(NR,NF,NPHI)
      REAL(8)               :: PHIR_DOT(NR,NF,NPHI)
      REAL(8)               :: PHIWORK2D(NR,NF,NF)
      INTEGER(4)            :: IF,IF1,IF2,IR
      REAL(8)               :: HA(NF-NPHI,NF-NPHI),HX(NF-NPHI,NPHI),HB(NF-NPHI,NPHI)
      REAL(8)               :: ALLKINK_HOM(NF,NF)
      REAL(8)               :: MAT(NF,NF)
      REAL(8)               :: DE1(NPHI)
      INTEGER(4)            :: IRC
      REAL(8)               :: KINK_HOM(NPHI,NPHI)
      REAL(8)               :: KINK_DOT(NPHI,NPHI)
      REAL(8)               :: KINKC(NPHI,NPHI)
      REAL(8)               :: HAM(NPHI,NPHI),OV(NPHI,NPHI)
      REAL(8)               :: R(NR)
      REAL(8)               :: DHOM(NR,NF)
      LOGICAL   ,PARAMETER  :: TWRITE=.FALSE.
      REAL(8)               :: BVECS(NF,NF),XVECS(NF,NF)
      REAL(8)               :: AUX(NR)
      REAL(8)               :: SVAR
      INTEGER(4)            :: I,J
      INTEGER(4)            :: L0
      REAL(8)               :: AMAT(NPHI-1,NPHI-1),BVEC(NPHI-1)
CHARACTER(32):: FILE
!     ******************************************************************
      TOK=.FALSE.
!PRINT*,'NEW RADIAL_XXXR STARTED',NPHI,NF
      IF(IROUT+1.GT.NR) THEN
        CALL ERROR$MSG('IROUT OUT OF RANGE')
        CALL ERROR$STOP('RADIAL_XXXC')
      END IF
      IF(IRMATCH.GT.IROUT) THEN
        CALL ERROR$MSG('IRMATCH OUT OF RANGE')
        CALL ERROR$STOP('RADIAL_XXXC')
      END IF
      L0=(NPHI-1)/2
      IF(NPHI.NE.2*L0+1) THEN
        CALL ERROR$MSG('NPHI MUST BE 2*L0+1')
        CALL ERROR$i4val('NPHI',nphi)
        CALL ERROR$i4val('l0',l0)
        CALL ERROR$STOP('RADIAL_XXXC')
      END IF
      CALL RADIAL$R(GID,NR,R)
      IRC=IRMATCH
!PRINT*,'IROUT ',IROUT,R(IROUT),IRC,R(IRC)
!
!     ==================================================================
!     ==  OBTAIN HOMOGENEOUS SOLUTION                                 ==
!     ==================================================================
!PRINT*,'RADIAL_XXXr: DETERMINE PHI',Lox
      ALLPHIL(:,:,:)=0.D0
      ALLPHIR(:,:,:)=0.D0
      DO IF=1,NF
        ALLPHIL(1:2,IF,IF)=R(1:2)**LOX(IF)
        DHOM(:,:)=0.D0
        CALL RADIAL$DGLGEN(GID,NR,NF,1,IRC+1,A,B,C,DHOM,ALLPHIL(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIL(IRC,:,IF)))
        ALLPHIL(1:IRC+1,:,IF)=ALLPHIL(1:IRC+1,:,IF)/SVAR
        DHOM(IROUT,IF)=1.D-8
        CALL RADIAL$DGLGEN(GID,NR,NF,IROUT+1,IRC-1,A,B,C,DHOM,ALLPHIR(:,:,IF))
        SVAR=MAXVAL(ABS(ALLPHIR(IRC,:,IF)))
        ALLPHIR(IRC-1:IROUT+1,:,IF)=ALLPHIR(IRC-1:IROUT+1,:,IF)/SVAR
      ENDDO
!
!     ==================================================================
!     ==  MAKE PHI_HOM SOLUTIONS CONTINUOUS                           ==
!     ==================================================================
!PRINT*,'RADIAL_XXXr: MATCH PHI'
      MAT(:,:)=ALLPHIR(IRC,:,:)
!     == MATRIX IS NOT SYMMETRIC. THUS SOLVE WITH SINGULAR VALUE DECOMPOSITION
      BVECS(:,1:NF)=ALLPHIL(IRC,:,:)
      CALL LIB$MATRIXSOLVENEW(NF,NF,NF,MAT,XVECS,BVECS)
      PHIWORK2D(:,:,:)=ALLPHIR(:,:,:)
      ALLPHIR(:,:,:)=0.D0
      DO IF1=1,NF
        DO IF2=1,NF
          ALLPHIR(:,:,IF1)=ALLPHIR(:,:,IF1)+PHIWORK2D(:,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTA CHANNELS        ==
!     ==================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHI'
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).EQ.L0) CYCLE
        I=I+1
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          SVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HA(I,J)=SVAR
        ENDDO
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).NE.L0) CYCLE
          J=J+1
          SVAR=(ALLPHIR(IRC+1,IF1,IF2)-ALLPHIR(IRC-1,IF1,IF2)) &
     &         -(ALLPHIL(IRC+1,IF1,IF2)-ALLPHIL(IRC-1,IF1,IF2))
          HB(I,J)=-SVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVENEW(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      PHIL(:,:,:)=0.D0
      PHIR(:,:,:)=0.D0
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).NE.L0) CYCLE
        I=I+1
        PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF1)
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                         +ALLPHIR(IRC-1:IROUT+1,:,IF1)
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ====================================================================
!     ==   ORTHOGONALIZE PHI                                                ==
!     ====================================================================
!PRINT*,'RADIAL_XXXC: ORTHOGONALIZE PHI'
      DO I=1,NPHI
!       == ORTHOGONALIZE
        DO J=1,I-1
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)=AUX(:IRC)+PHIL(:IRC,IF,J)*PHIL(:IRC,IF,I)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &                   +PHIR(IRC+1:IROUT+1,IF,J)*PHIR(IRC+1:IROUT+1,IF,I)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PHIL(:IRC+1,:,I)=PHIL(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*SVAR
          PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*SVAR
        ENDDO
!       ==  NORMALIZE
        AUX(:)=0.D0
        DO IF=1,NF
          AUX(:IRC)         =AUX(:IRC)         +PHIL(:IRC,IF,I)**2
          AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1)+PHIR(IRC+1:IROUT+1,IF,I)**2
        ENDDO
        AUX(:)=AUX(:)*R(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        SVAR=1.D0/SQRT(SVAR)
        PHIL(:IRC+1,:,I)       =PHIL(:IRC+1,:,I)*SVAR
        PHIR(IRC-1:IROUT+1,:,I)=PHIR(IRC-1:IROUT+1,:,I)*SVAR
      ENDDO
!
!!$!     ====================================================================
!!$!     ==   TEST CONTINUITY                                              ==
!!$!     ====================================================================
!!$!PRINT*,'RADIAL_XXXr: TEST CONTINUITY'
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          SVAR1=ABS(PHIL(IRC,IF2,IF1)-PHIR(IRC,IF2,IF1))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$!PRINT*,'DEVIATION IN VALUE',SVAR2
!!$      IF(SVAR2.GT.1.D-6) THEN
!!$        CALL ERROR$MSG('MATCHING OF PHI FAILED')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF   
!!$!
!!$!     == TEST DIFFERENTIABILITY ===================================
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          IF(LOX(IF2).EQ.L0) CYCLE
!!$          SVAR1=ABS( (PHIR(IRC+1,IF2,IF1)-PHIR(IRC-1,IF2,IF1)) &
!!$     &              -(PHIL(IRC+1,IF2,IF1)-PHIL(IRC-1,IF2,IF1)))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$!PRINT*,'DEVIATION IN OTHER DERIVATIVES',SVAR2
!!$      IF(SVAR2.GT.1.D-6) THEN
!!$        CALL ERROR$MSG('MATCHING OF DPHI/DR FAILED')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF   
!!$!     ==  TEST OVERLAP  =============================================
!!$!      PRINT*,'OVERLAP'
!!$      DO I=1,NPHI
!!$        DO J=I,NPHI
!!$          CAUX(:)=0.D0
!!$          DO IF=1,NF
!!$            CAUX(:IRC)=CAUX(:IRC)+PHIL(:IRC,IF,I)*PHIL(:IRC,IF,J)
!!$            CAUX(IRC+1:)=CAUX(IRC+1:)+PHIR(IRC+1:,IF,I)*PHIR(IRC+1:,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX(:)=CAUX
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          OV(I,J)=SVAR1
!!$          OV(J,I)=OV(I,J)
!!$        ENDDO
!!$      ENDDO
!!$      DO I=1,NPHI
!!$        SVAR=1.D0/SQRT(ABS(OV(I,I)))
!!$        OV(:,I)=OV(:,I)*SVAR
!!$        OV(I,:)=OV(I,:)*SVAR
!!$      ENDDO
!!$!PRINT*,'MARKE A',NPHI
!!$      DO I=1,NPHI
!!$        WRITE(*,FMT='(20("(",F10.3,",",F10.3,")"))')OV(I,:)
!!$      ENDDO
!!$!PRINT*,'MARKE B',NPHI
!!$!
!!$!     ==  PLOT WAVE FUNCTIONS ======================================
!!$      DO IF1=1,NPHI
!!$        WRITE(FILE,*)IF1
!!$        FILE='PHI'//ADJUSTL(FILE)
!!$        OPEN(8,FILE=FILE,FORM='FORMATTED')
!!$        REWIND 8
!!$        DO IR=1,IRC+1
!!$          WRITE(8,FMT='(80F12.7)')R(IR),PHIL(IR,:,IF1)
!!$        ENDDO
!!$        DO IR=IRC-1,MIN(IROUT+1,NR)
!!$          WRITE(8,FMT='(80F12.7)')R(IR),PHIR(IR,:,IF1)
!!$        ENDDO
!!$        CLOSE(8)
!!$      ENDDO
!
!     ==================================================================
!     ==  DETERMINE PHIDOT                                            ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC: DETERMINE PHIDOT'
      PHIL_DOT(:,:,:)=0.D0
      PHIR_DOT(:,:,:)=0.D0
      DO IF=1,NPHI
        CALL RADIAL$DGLGEN(GID,NR,NF,1,IRC+1,A,B,C,-2.D0*PHIL(:,:,IF) &
     &                     ,PHIL_DOT(:,:,IF))
        CALL RADIAL$DGLGEN(GID,NR,NF,IROUT+1,IRC-1,A,B,C,-2.D0*PHIR(:,:,IF) &
     &                     ,PHIR_DOT(:,:,IF))
      ENDDO
!
!     ==================================================================
!     ==  MAKE PHI_DOT CONTINUOUS                                     ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC: MATCH PHIDOT'
      MAT(:,:)=ALLPHIR(IRC,:,:)
      BVECS(:,:NPHI)=PHIL_DOT(IRC,:,:)-PHIR_DOT(IRC,:,:)
      CALL LIB$MATRIXSOLVENEW(NF,NF,NPHI,MAT,XVECS(:,:NPHI),BVECS(:,:NPHI))
!
      DO IF1=1,NPHI
        DO IF2=1,NF
          PHIR_DOT(IRC-1:IROUT+1,:,IF1)=PHIR_DOT(IRC-1:IROUT+1,:,IF1) &
     &                            +ALLPHIR(IRC-1:IROUT+1,:,IF2)*XVECS(IF2,IF1)
        ENDDO
      ENDDO
!
!     ==================================================================
!     == REMOVE KINKS IN NON-RELEVANT ANGULAR MOMENTUM CHANNELS        ==
!     ==================================================================
!PRINT*,'REMOVE KINKS OF OTHER ANGULAR MOMENTUM CHANNELS OF PHIDOT'
      I=0
      DO IF1=1,NF
        IF(LOX(IF1).EQ.L0) CYCLE
        I=I+1
        DO J=1,NPHI
          SVAR=(PHIR_DOT(IRC+1,IF1,J)-PHIR_DOT(IRC-1,IF1,J)) &
     &         -(PHIL_DOT(IRC+1,IF1,J)-PHIL_DOT(IRC-1,IF1,J))
          HB(I,J)=-SVAR
        ENDDO
      ENDDO
      IF(NF.GT.NPHI) THEN
        CALL LIB$MATRIXSOLVENEW(NF-NPHI,NF-NPHI,NPHI,HA,HX,HB)
      END IF
      DO I=1,NPHI
        J=0
        DO IF2=1,NF
          IF(LOX(IF2).EQ.L0) CYCLE
          J=J+1
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)+ALLPHIL(:IRC+1,:,IF2)*HX(J,I)
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I) &
     &                               +ALLPHIR(IRC-1:IROUT+1,:,IF2)*HX(J,I)
        ENDDO
      ENDDO  
!
!     ====================================================================
!     ==   ORTHOGONALIZE PHIDOT TO PHI                                  ==
!     ====================================================================
      DO I=1,NPHI
        DO J=1,NPHI
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)=AUX(:IRC)+PHIL(:IRC,IF,J)*PHIL_DOT(:IRC,IF,I)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &             +PHIR(IRC+1:IROUT+1,IF,J)*PHIR_DOT(IRC+1:IROUT+1,IF,I)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          PHIL_DOT(:IRC+1,:,I)=PHIL_DOT(:IRC+1,:,I)-PHIL(:IRC+1,:,J)*SVAR
          PHIR_DOT(IRC-1:IROUT+1,:,I)=PHIR_DOT(IRC-1:IROUT+1,:,I)-PHIR(IRC-1:IROUT+1,:,J)*SVAR
        ENDDO
      ENDDO
!
!     ====================================================================
!     ==   TEST CONTINUITY                                              ==
!     ====================================================================
!!$PRINT*,'RADIAL_XXXC: TEST CONTINUITY OF PHIDOT'
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          SVAR1=ABS(PHIL_DOT(IRC,IF2,IF1)-PHIR_DOT(IRC,IF2,IF1))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$PRINT*,'DEVIATION IN VALUE OF PHIDOT ',SVAR2
!!$      IF(SVAR2.GT.1.D-6) THEN
!!$        CALL ERROR$MSG('MATCHING FAILED')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF   
!!$!
!!$!     == TEST DIFFERENTIABILITY ===================================
!!$      SVAR2=0.D0
!!$      DO IF1=1,NPHI
!!$        DO IF2=1,NF
!!$          IF(LOX(IF2).EQ.L0) CYCLE
!!$          SVAR1=ABS( (PHIR_DOT(IRC+1,IF2,IF1)-PHIR_DOT(IRC-1,IF2,IF1)) &
!!$     &              -(PHIL_DOT(IRC+1,IF2,IF1)-PHIL_DOT(IRC-1,IF2,IF1)))
!!$          SVAR2=MAX(SVAR2,SVAR1)
!!$        ENDDO
!!$      ENDDO
!!$PRINT*,'DEVIATION IN OTHER DERIVATIVES OF PHIDOT',SVAR2
!!$      IF(SVAR2.GT.1.D-4) THEN
!!$        CALL ERROR$MSG('MATCHING OF DPHIDOT/DR FAILED')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF   
!!$!
!!$!     =============================================================
!!$      DO IF1=1,NPHI
!!$        WRITE(FILE,*)IF1
!!$        FILE='PHIDOT'//ADJUSTL(FILE)
!!$        OPEN(8,FILE=FILE,FORM='FORMATTED')
!!$        REWIND 8
!!$        DO IR=1,IRC+1
!!$          WRITE(8,FMT='(80F12.7)')R(IR),REAL(PHIL_DOT(IR,:,IF1)),AIMAG(PHIL_DOT(IR,:,IF1))
!!$        ENDDO
!!$        DO IR=IRC-1,MIN(IROUT+1,NR)
!!$          WRITE(8,FMT='(80F12.7)')R(IR),REAL(PHIR_DOT(IR,:,IF1)),AIMAG(PHIR_DOT(IR,:,IF1))
!!$        ENDDO
!!$        CLOSE(8)
!!$      ENDDO
!
!     ==================================================================
!     ==  REMOVE KINKS BY MIXING PHIDOT INTO PHI                      ==
!     ==================================================================
!PRINT*,'RADIAL_XXXC:MATCH KINKS'
      I=0
      DO IF2=1,NF
        IF(LOX(IF2).NE.L0) CYCLE
        I=I+1
        KINK_HOM(I,:)=(PHIR(IRC+1,IF2,:)-PHIR(IRC-1,IF2,:)) &
     &               -(PHIL(IRC+1,IF2,:)-PHIL(IRC-1,IF2,:))
        KINK_DOT(I,:)=(PHIR_DOT(IRC+1,IF2,:)-PHIR_DOT(IRC-1,IF2,:)) &
     &               -(PHIL_DOT(IRC+1,IF2,:)-PHIL_DOT(IRC-1,IF2,:))
      ENDDO
      CALL LIB$MATRIXSOLVENEW(NPHI,NPHI,NPHI,-KINK_DOT,HAM,KINK_HOM)
!!$      SVAR=MAXVAL(ABS(KINK_HOM+MATMUL(KINK_DOT,HAM)))
!!$PRINT*,'REMAINING KINKS ',SVAR
!     == MAKE PHI DIFFERENTIABLE =========================================
      DO I=1,NPHI
        DO J=1,NPHI
          PHIL(:IRC,:,I)         =PHIL(:IRC,:,I)         +PHIL_DOT(:IRC,:,J)  *HAM(J,I)
          PHIR(IRC+1:IROUT+1,:,I)=PHIR(IRC+1:IROUT+1,:,I)+PHIR_DOT(IRC+1:IROUT+1,:,J)*HAM(J,I)
        ENDDO
      ENDDO
!     == MAKE PHI DIFFERENTIABLE =========================================
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("HAM ",50("(",2F10.5")   "))')HAM(I,:)
!!$ENDDO
!SVAR=0.5D0*MAXVAL(ABS(HAM-TRANSPOSE(CONJG(HAM))))
!PRINT*,'DEVIATION FROM HERMITICITY',SVAR
!
!     ==================================================================
!     ==  DETERMINE OVERLAP MATRIX                                    ==
!     ==================================================================
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("H ",50("(",2F10.5")   "))')HAM(I,:)
!!$ENDDO
!!$
!!$      DO I=1,NPHI
!!$        DO J=I,NPHI
!!$          CAUX(:)=0.D0
!!$          DO IF=1,NF
!!$            CAUX(:IRC)=CAUX(:IRC)+CONJG(PHIL_DOT(:IRC,IF,I))*PHIL_DOT(:IRC,IF,J)
!!$            CAUX(IRC+1:IROUT+1)=CAUX(IRC+1:IROUT+1) &
!!$     &             +CONJG(PHIR_DOT(IRC+1:IROUT+1,IF,I))*PHIR_DOT(IRC+1:IROUT+1,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX(:)=REAL(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          AUX(:)=AIMAG(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$          OV(I,J)=CMPLX(SVAR1,SVAR2)
!!$          OV(J,I)=CONJG(OV(I,J))
!!$        ENDDO
!!$      ENDDO
!!$      OV=MATMUL(TRANSPOSE(CONJG(HAM)),MATMUL(OV,HAM))
!!$      DO I=1,NPHI
!!$        OV(I,I)=OV(I,I)+(1.D0,0.D0)
!!$      ENDDO
!!$!
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("O ",50("(",2F10.5")   "))')OV(I,:)
!!$ENDDO
!
      DO I=1,NPHI
        DO J=I,NPHI
          AUX(:)=0.D0
          DO IF=1,NF
            AUX(:IRC)         =AUX(:IRC)+PHIL(:IRC,IF,I)*PHIL(:IRC,IF,J)
            AUX(IRC+1:IROUT+1)=AUX(IRC+1:IROUT+1) &
     &                         +PHIR(IRC+1:IROUT+1,IF,I)*PHIR(IRC+1:IROUT+1,IF,J)
          ENDDO
          AUX(:)=AUX(:)*R(:)**2
          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
          OV(I,J)=SVAR
          OV(J,I)=OV(I,J)
        ENDDO
      ENDDO
!!$DO I=1,NPHI
!!$  WRITE(*,FMT='("O ",50("(",2F10.5")   "))')OV(I,:)
!!$ENDDO
!
!     ==================================================================
!     ==  DETERMINE EIGENSTATES                                       ==
!     ==================================================================
      HAM=0.5D0*(HAM+TRANSPOSE(HAM))
!      CALL LIB$DIAGC8(NPHI,HAM,DE,KINKC)
      CALL LIB$GENERALEIGENVALUER8(NPHI,HAM,OV,DE,KINKC)
      PHI(:,:,:)=0.D0
      DO I=1,NPHI
        DO J=1,NPHI
          PHI(:IRC,:,I)         =PHI(:IRC,:,I)         +PHIL(:IRC,:,J)*KINKC(J,I)      
          PHI(IRC+1:IROUT+1,:,I)=PHI(IRC+1:IROUT+1,:,I)+PHIR(IRC+1:IROUT+1,:,J)*KINKC(J,I)      
        ENDDO
      ENDDO
!PRINT*,'DE',DE
!!$!
!!$!     ==================================================================
!!$!     ==  NORMALIZE SOLUTIONS                                         ==
!!$!     ==================================================================
!!$      DO I=1,NPHI
!!$        AUX(:)=0.D0
!!$        DO J=1,NF
!!$          AUX(:)=AUX(:)+ABS(PHI(:,J,I))**2
!!$        ENDDO
!!$        AUX(:)=AUX(:)*R(:)**2
!!$        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
!!$        SVAR=1.D0/SQRT(SVAR)
!!$        PHI(:,:,I)=PHI(:,:,I)*SVAR
!!$      ENDDO
!
!     ==================================================================
!     ==  TEST OVERLAP                                                ==
!     ==================================================================
!!$      DO IF1=1,NPHI
!!$        WRITE(FILE,*)IF1
!!$        FILE='PHIFIN'//ADJUSTL(FILE)
!!$        OPEN(8,FILE=FILE,FORM='FORMATTED')
!!$        REWIND 8
!!$        DO IR=1,IROUT+1
!!$          WRITE(8,FMT='(80F20.7)')R(IR),REAL(PHI(IR,:,IF1)),AIMAG(PHI(IR,:,IF1))
!!$        ENDDO
!!$        CLOSE(8)
!!$      ENDDO
!!$!
!     == TEST ORTHONORMALITY OVERLAP MATRIX
!!$      DO I=1,NPHI
!!$        HAM(I,I)=(0.D0,0.D0)
!!$        DO J=I,NPHI
!!$          CAUX(:)=0.D0
!!$          DO IF=1,NF
!!$            CAUX(:)=CAUX(:)+CONJG(PHI(:,IF,I))*PHI(:,IF,J)
!!$          ENDDO
!!$          CAUX(:)=CAUX(:)*R(:)**2
!!$          AUX=REAL(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR1)
!!$          AUX=AIMAG(CAUX)
!!$          CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR2)
!!$          HAM(I,J)=CMPLX(SVAR1,SVAR2)
!!$          HAM(J,I)=CONJG(HAM(I,J))
!!$        ENDDO
!!$        HAM(I,I)=HAM(I,I)-(1.D0,0.D0)
!!$      ENDDO
!!$      DO I=1,NPHI
!!$        WRITE(*,FMT='("O ",50("(",2F10.5")   "))')HAM(I,:)
!!$      ENDDO
!!$      SVAR=MAXVAL(ABS(HAM))
!!$      PRINT*,'DEVIATION FROM ORTHONORMALITY',SVAR,NF
!!$      IF(SVAR.GT.0.5D0) THEN
!!$        CALL ERROR$MSG('WAVE FUNCTIONS NOT ORTHOGONAL')
!!$        CALL ERROR$STOP('RADIAL_XXXC')
!!$      END IF
      TOK=.TRUE.
      RETURN
      END SUBROUTINE RADIAL_XXXR_OV
!
!     ...................................................................
      subroutine schroedinger_specialrads(gid,nr,l,xmax,v00,e,ircl,irout)
!     **                                                               **
!     **  estimate the classical turning point r(ircl)                 **
!
!     **  estimate the outermost grid point for inward integration     **
!     **  from a WKB solution of the Schroedinger equation             **
!     **                                                               **
      implicit none
      integer(4),intent(in) :: gid      ! grid id
      integer(4),intent(in) :: nr       ! #(radial grid points)
      integer(4),intent(in) :: l        ! angular momentum
      real(8)   ,intent(in) :: xmax     ! maximum tolerable ratio of phi 
      real(8)   ,intent(in) :: v00(nr)  ! radial spherical potential
      real(8)   ,intent(in) :: e        ! energy
      integer(4),intent(out):: ircl     ! classical turning point
      integer(4),intent(out):: irout    ! outermost grid point
      real(8)               :: pi       ! pi
      real(8)               :: y0       ! spherical harmonic for lm=0
      real(8)               :: r(nr)    ! radial grid
      real(8)               :: fac,svar,sumval,xmaxlog
      integer(4)            :: ir
      logical               :: tchk
!     *******************************************************************
      pi=4.d0*datan(1.d0)
      y0=1.d0/sqrt(4.d0*pi)
      xmaxlog=log(xmax)
      call radial$r(gid,nr,r)
      fac=0.5d0*real(l*(l+1),kind=8)
      sumval=0.d0
      irout=1
      ircl=nr
      tchk=.false.
      do ir=2,nr-1
        svar=v00(ir)*y0+fac/r(ir)**2-e
        if(svar.lt.0.d0) then
          svar=0.d0
          sumval=0.d0
          ircl=ir+1
          cycle
        end if
        sumval=sumval+sqrt(2.d0*svar)*0.5d0*(r(ir+1)-r(ir-1))
        irout=ir-1
        if(sumval.gt.xmaxlog) then
          tchk=.true.
          exit
        end if
      enddo
!     == fix up end of the grid
      if(.not.tchk) then
        svar=v00(nr)*y0+fac/r(nr)**2-e
        sumval=sumval+0.5d0*(r(nr)-r(nr-1))*svar
        if(sumval.gt.xmaxlog) then
          irout=nr-1
        else
          irout=nr
        end if
      end if
      return
      end



!     ....................................................................
      subroutine writephi(gid,nr,file,nf,nphi,irc,phil,phir)
      implicit none
      integer(4)  ,intent(in) :: gid
      character(*),intent(in) :: file
      integer(4)  ,intent(in) :: nr
      integer(4)  ,intent(in) :: nf
      integer(4)  ,intent(in) :: irc
      integer(4)  ,intent(in) :: nphi
      real(8)     ,intent(in) :: phil(nr,nf,nphi)
      real(8)     ,intent(in) :: phir(nr,nf,nphi)
      integer(4)              :: ir,iphi
      character(32)           :: file1
      real(8)                 :: r(nr)
!     **********************************************************************       
      call radial__r(gid,nr,r)
      do iphi=1,nphi
        write(FILE1,fmt='(i1)')IPHI
        FILE1=TRIM(FILE)//'_'//TRIM(FILE1)//'.DAT'
        open(100,file=file1)
        do ir=1,irc
          write(100,*)r(IR),phil(ir,:,iphi)
        enddo
        do ir=irc+1,nr
          write(100,*)R(IR),phir(ir,:,iphi)
        enddo
        close(100)
      enddo
      return 
      end
