!*******************************************************************************
!*******************************************************************************
!**  PAW_GAUSSIAN OBJECT                                                      **
!**                                                                           **
!**  PERFORMS OPERATIONS OF CARTESIAN AND HERMITE GAUSSIANS                   **
!**                                                                           **
!**  INDEXING: SEE ROUTINE GAUSSIAN_INDEX FOR EXPLANATION                     **
!**  
!**  
!**  
!**  
!**  
!*******************************************************************************
!*******************************************************************************
!     
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_INDEX(N,NIJK,IJK)
!     **************************************************************************
!     **  DEMONSTRATES THE INDEXING OF 3-D GAUSSIANS                          **
!     **                                                                      **
!     **  GAUSSIANS ARE ORDERED ACCORDING TO THEIR HIGHEST TOTAL POWER,       **
!     **  BECAUSE THAT IS THE NATURAL CUTOFF FOR AN EXPANSION                 **
!     **                                                                      **
!     **  CARTESIAN GAUSSIANS ARE DEFINED AS                                  **
!     **  G(IND)=X^IJK(1,IND) * Y^IJK(2,IND) * Z^IJK(3,IND) * EXP(-P*R^2)     **
!     **************************************************************************
      INTEGER(4),INTENT(IN)  :: N           ! HIGHEST POWER FOR GAUSSIANS
      INTEGER(4),INTENT(IN)  :: NIJK        ! #(GAUSSIANS)
      INTEGER(4),INTENT(OUT) :: IJK(3,NIJK) ! INDEX ARRAY
      INTEGER(4)             :: IND         ! GAUSSIAN INDEX
      INTEGER(4)             :: M           ! M=I+J+K
      INTEGER(4)             :: I,J,K       ! POWERS DEFINING THE GAUSSIAN
!     **************************************************************************
      IF(NIJK.NE.(N+1)*(N+2)*(N+3)/6) THEN
        CALL ERROR$MSG('INCONSISTENT DIMENSIONS')
        CALL ERROR$STOP('GAUSSIAN_INDEX')
      END IF
!
!     ==========================================================================
!     == LOOP OVER ALL GAUSSIANS WITH I+J+K=0,...,N                           ==
!     ==========================================================================
      IND=0
      DO M=0,N              ! M=I+J+K IS THE ACTUAL TOTAL POWER 
        DO I=0,M            ! I GROWS SLOWEST
          DO J=0,M-I        ! J GROWS
            K=M-I-J         ! K SHRINKS
            IND=IND+1
            IJK(1,IND)=I    ! POWER FOR THE X-DIRECTION
            IJK(2,IND)=J    ! POWER FOR THE Y-DIRECTION
            IJK(3,IND)=K    ! POWER FOR THE Z-DIRECTION
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_GAUSSINDEX(ID,IND,I,J,K)
!     **************************************************************************
!     ** TRANSFORMS GAUSS INDEX IND INTO POWERS (I,J,K) FOR XYZ AND VICE VERSA**
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN)    :: ID    ! SWITCH: 'IJKFROMIND' OR 'INDFROMIJK'
      INTEGER(4)  ,INTENT(INOUT) :: IND   ! GAUSSIAN INDEX
      INTEGER(4)  ,INTENT(INOUT) :: I,J,K ! POWERS FOR XYZ
      INTEGER(4)                 :: N,IND1,I1
!     **************************************************************************
!
!     ==========================================================================
!     == DETERMINE (I,J,K) FROM GAUSSIAN INDEX 'IND'                          ==
!     ==========================================================================
      IF(ID.EQ.'IJKFROMIND') THEN
        N=0
        DO WHILE ((N+1)*(N+2)*(N+3).LT.6*IND)
          N=N+1
        ENDDO
        IND1=(N+1)*(N+2)*(N+3)/6
        IF(IND1.EQ.IND) THEN
          I=N
          J=0
          K=0
          RETURN
        END IF
        DO I1=0,N
          IF(2*IND1-(N-I1+1)*(N-I1+2)+2.GT.2*IND) THEN
            I=I1-1
            J=(2*IND-2*IND1+(N-I+1)*(N-I+2)-2)/2
            K=N-I-J
            EXIT
          END IF
        ENDDO
        RETURN
!
!     ==========================================================================
!     == DETERMINE GAUSSIAN INDEX 'IND' FROM (I,J,K)                          ==
!     ==========================================================================
      ELSE IF(ID.EQ.'INDFROMIJK') THEN  
        N=I+J+K
        IND=(N+1)*(N+2)*(N+3)/6-(N-I+1)*(N-I+2)/2+J+1
      ELSE
        CALL ERROR$MSG('ID NOT RECOGNIZED')
        CALL ERROR$CHVAL('ID',ID)
        CALL ERROR$STOP('GAUSSIAN$GAUSSINDEX')
      END IF
      RETURN
      END   
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_YLMPOL(LX,YLMPOL)
!     **************************************************************************
!     ** R^L*Y_LM = SUM_IJK X^IY^JZ^K * YLMPOL(IJK,LM)                        **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LX
      REAL(8)   ,INTENT(OUT):: YLMPOL((LX+1)*(LX+2)*(LX+3)/6,(LX+1)**2)
!     **************************************************************************
      IF(LX.LT.0) THEN
        CALL ERROR$MSG('LX MUST NOT BE NEGATIVE')
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('GAUSSIAN_YLMPOL')
      END IF
      CALL SPHERICAL$YLMPOLYNOMIALS(LX,(LX+1)*(LX+2)*(LX+3)/6,(LX+1)**2,YLMPOL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_YLMFROMPOL(LX,YLMFROMPOL)
!     **************************************************************************
!     ** EXPANDS THE ANGULAR PART OF CARTESIAN GAUSSIANS INTO REAL SPHERICAL  **
!     ** HARMONICS                                                            **
!     **  IND->(I,J,K) ;                                                      **
!     **  (X/R)^I * (Y/R)^J * (Z/R)^K                                         **
!     **         =SUM YLM(R) * |R|^(L+2N) * YLMFROMPOL(N+1,LM,IND)            **
!     **  (DESCRIPTION IN CHAPTER "WORKING WITH GAUSSIANS")                   **
!     ***************************************PETER BLOECHL, GOSLAR 2010*********
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: LX
      REAL(8)   ,INTENT(OUT):: YLMFROMPOL(LX/2+1,(LX+1)**2 &
     &                                          ,(LX+1)*(LX+2)*(LX+3)/6)
      REAL(8)               :: YLMPOL((LX+1)*(LX+2)*(LX+3)/6,(LX+1)**2)
      REAL(8)               :: AMAT((LX+1)*(LX+2)*(LX+3)/6 &
     &                             ,(LX+1)*(LX+2)*(LX+3)/6)
      REAL(8)               :: AINV((LX+1)*(LX+2)*(LX+3)/6 &
     &                             ,(LX+1)*(LX+2)*(LX+3)/6)
      INTEGER(4)            :: INDX
      INTEGER(4)            :: LMX,LMX1,LX1,INDX1
      INTEGER(4)            :: I0A,I0B
      INTEGER(4)            :: N,IND,IND1,I,J,K,IP2,JP2,KP2
      INTEGER(4)            :: LM,L,M
      LOGICAL(4)            :: TWRITE=.FALSE.
!     **************************************************************************
      IF(LX.LT.0) THEN
        CALL ERROR$MSG('LX MUST NOT BE NEGATIVE')
        CALL ERROR$I4VAL('LX',LX)
        CALL ERROR$STOP('GAUSSIAN_YLMFROMPOL')
      END IF
      INDX=(LX+1)*(LX+2)*(LX+3)/6
      LMX=(LX+1)**2
      CALL SPHERICAL$YLMPOLYNOMIALS(LX,INDX,LMX,YLMPOL)
      AMAT=0.D0
      AMAT(:,:LMX)=YLMPOL(:,:)        
      I0A=0
      I0B=LMX
      DO N=1,LX/2
        LX1=LX-2*N
        LMX1=(LX1+1)**2
        INDX1=(LX-1)*(LX)*(LX+1)/6  ! LEAVE SPACE FOR R^2
!PRINT*,'N=',N,'I0A,I0B=',I0A+1,I0A+LMX1,I0B+1,I0B+LMX1,'LX1,LMX1,INDX1',LX1,LMX1,INDX1
        DO IND=1,INDX1
          CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I,J,K)
          IP2=I+2
          JP2=J+2
          KP2=K+2
          CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND1,IP2,J,K)
          AMAT(IND1,I0B+1:I0B+LMX1)=AMAT(IND1,I0B+1:I0B+LMX1) &
     &                             +AMAT(IND,I0A+1:I0A+LMX1)
          CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND1,I,JP2,K)
          AMAT(IND1,I0B+1:I0B+LMX1)=AMAT(IND1,I0B+1:I0B+LMX1) &
     &                             +AMAT(IND,I0A+1:I0A+LMX1)
          CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND1,I,J,KP2)
          AMAT(IND1,I0B+1:I0B+LMX1)=AMAT(IND1,I0B+1:I0B+LMX1) &
     &                             +AMAT(IND,I0A+1:I0A+LMX1)
        ENDDO
        I0A=I0B
        I0B=I0B+LMX1
      ENDDO
!
!     ==========================================================================
!     == WRITE POWER SERIES EXPANSION                                         ==
!     ==========================================================================
      IF(TWRITE) THEN
        IND1=0
        DO N=0,LX/2
          LM=0
          DO L=0,LX-2*N
            DO M=-L,L
              LM=LM+1
              IND1=IND1+1
              WRITE(*,*)'=== L=',L,' M=',M,' LM=',LM,' N=',N,' IND1=',IND1,'==='
              DO IND=1,INDX
                IF(AMAT(IND,IND1).EQ.0.D0) CYCLE
                CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I,J,K)
                WRITE(*,*)I,J,K,AMAT(IND,IND1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
!
!     ==========================================================================
!     ==  INVERT MATRIX                                                       ==
!     ==========================================================================
      CALL LIB$INVERTR8(INDX,AMAT,AINV)
!
!     ==========================================================================
!     ==  RESOLVE                                                             ==
!     ==========================================================================
      YLMFROMPOL(:,:,:)=0.D0
      YLMFROMPOL(1,:LMX,:)=AINV(:LMX,:)
      I0A=LMX
      DO N=1,LX/2
        LMX1=(LX-2*N+1)**2
        YLMFROMPOL(N+1,:LMX1,:)=AINV(I0A+1:I0A+LMX1,:)
        I0A=I0A+LMX1
      ENDDO              
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_SHIFTCENTER(NX,NIJK,E,R,C)
!     **************************************************************************
!     **  SHIFTS THE CENTER FOR THE EXPANSION OF A FUNCTION EXPRESSED IN      **
!     **  CARTESIAN COORDINATES TO A NEW CENTER                               **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4)  ,INTENT(IN)  :: NX     ! HIGHEST POWER IN EXPANSION
      INTEGER(4)  ,INTENT(IN)  :: NIJK   ! #(COEFFICIENTS IN C)
      REAL(8)     ,INTENT(IN)  :: E      ! EXPONENT
      REAL(8)     ,INTENT(IN)  :: R(3)   ! POSITION RELATIVE TO THE CENTER
      REAL(8)     ,INTENT(OUT) :: C(NIJK,NIJK) ! TRANSF. MATRIX FOR COEFFICIENTS
      REAL(8)                  :: C1(0:NX,0:NX)
      REAL(8)                  :: C2(0:NX,0:NX)
      REAL(8)                  :: C3(0:NX,0:NX)
      REAL(8)                  :: FAC(3),FAC1(3),XJ
      INTEGER(4)               :: J,K,N,IND1,IND2,I1,I2,J1,J2,K1,K2,M1,M2,NMK
      REAL(8)                  :: GARR(0:NX,3)
      REAL(8)                  :: SVAR
!     **************************************************************************
!     == WITHOUT A SHIFT, RETURN THE IDENTITY ==================================
      IF(SUM(R(:)**2).LT.1.D-10) THEN
        C(:,:)=0.D0
        DO N=1,NIJK
          C(N,N)=1.D0
        ENDDO
        RETURN
      END IF
!
!     ==========================================================================
!     == EXPAND GAUSSIAN                                                      ==
!     ==========================================================================
      FAC(:)=EXP(-E*R(:)**2)
      GARR(0,:)=FAC(:)
      XJ=0.D0
      DO J=1,NX
        XJ=XJ+1.D0
        FAC(:)=FAC(:)*2.D0*E*R(:)/XJ
        GARR(J,:)=FAC(:)
      ENDDO
!
!     ==========================================================================
!     == EXPAND POLYNOMIAL ABOUT OTHER CENTER AND ASSEMBLE RESULT             ==
!     ==========================================================================
      C1(:,:)=0.D0
      C2(:,:)=0.D0
      C3(:,:)=0.D0
      FAC(:)=1.D0
      DO K=0,NX
        DO N=K,NX
          CALL BINOMIALCOEFFICIENT(N,K,SVAR)
          FAC1(:)=FAC(:)*SVAR
          NMK=N-K
          DO J=0,NX-NMK
            C1(J+NMK,N)=C1(J+NMK,N)+GARR(J,1)*FAC1(1)
            C2(J+NMK,N)=C2(J+NMK,N)+GARR(J,2)*FAC1(2)
            C3(J+NMK,N)=C3(J+NMK,N)+GARR(J,3)*FAC1(3) 
          ENDDO
        ENDDO
        FAC(:)=-R(:)*FAC(:)
      ENDDO
!!$      == THIS IS THE SAFE VERSION
!!$      C1(:,:)=0.D0
!!$      C2(:,:)=0.D0
!!$      C3(:,:)=0.D0
!!$      DO N=0,NX
!!$        DO J=0,NX
!!$          K1=MAX(N+J-NX,0)
!!$          FAC(:)=(-R(:))**K1
!!$          DO K=K1,N
!!$            CALL BINOMIALCOEFFICIENT(N,K,SVAR)
!!$            C1(N+J-K,N)=C1(N+J-K,N)+GARR(J,1)*SVAR*FAC(1)
!!$            C2(N+J-K,N)=C2(N+J-K,N)+GARR(J,2)*SVAR*FAC(2)
!!$            C3(N+J-K,N)=C3(N+J-K,N)+GARR(J,3)*SVAR*FAC(3) 
!!$            FAC(:)=-R(:)*FAC(:)
!!$          ENDDO
!!$        ENDDO
!!$      ENDDO
!
!     ==========================================================================
!     == COMPOSE RESULT IN THREE DIMENSIONS                                   ==
!     ==========================================================================
      IND1=0
      DO M1=0,NX
        DO I1=0,M1
          DO J1=0,M1-I1
            K1=M1-I1-J1
            IND1=IND1+1
!
            IND2=0
            DO M2=0,NX
              DO I2=0,M2
                DO J2=0,M2-I2
                  K2=M2-I2-J2
                  IND2=IND2+1
! 
                  C(IND2,IND1)=C1(I2,I1)*C2(J2,J1)*C3(K2,K1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_3DORB(ID,NIJK,E,C,R,G)
!     **************************************************************************
!     **  CALCULATES THE VALUE OF A FUNCTION EXPRESSED AS AN EXPANSION INTO   **
!     **  CARTESIAN OR HERMITE GAUSSIANS. THE FUNCTION IS GIVEN BY THE        **
!     **  COEFFICIENT ARRAY C                                                 **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: ID  !SELECTS "CARTESIAN" OR "HERMITE" GAUSSIANS
      INTEGER(4)  ,INTENT(IN) :: NIJK    ! #(COEFFICIENTS IN C)
      REAL(8)     ,INTENT(IN) :: E       ! EXPONENT
      REAL(8)     ,INTENT(IN) :: R(3)    ! POSITION RELATIVE TO THE CENTER
      REAL(8)     ,INTENT(IN) :: C(NIJK) ! COEFFICIENT ARRAY 
      REAL(8)     ,INTENT(OUT):: G       ! VALUE OF THE FUNCTION
      INTEGER(4)              :: N,M,I,J,K,IND,NIJK1
      REAL(8)                 :: GX,GY,GZ
!     **************************************************************************
      IF(ID.EQ.'CARTESIAN') THEN
        G=0.D0
        DO IND=1,NIJK
          CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I,J,K)
          G=G+C(IND)*R(1)**I*R(2)**J*R(3)**K
        ENDDO
        G=G*EXP(-E*SUM(R(:)**2))
        RETURN
      END IF 
!
      NIJK1=NIJK  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJK1,N,J,K)
      G=0.D0
      IND=0
      DO M=0,N
        DO I=0,M
          DO J=0,M-I
            K=M-I-J
            IND=IND+1
            IF(ID.EQ.'CARTESIAN') THEN
              CALL GAUSSIAN_CARTESIANGAUSSIAN1D(I,E,R(1),GX)
              CALL GAUSSIAN_CARTESIANGAUSSIAN1D(J,E,R(2),GY)
              CALL GAUSSIAN_CARTESIANGAUSSIAN1D(K,E,R(3),GZ)
            ELSE IF(ID.EQ.'HERMITE') THEN
              CALL GAUSSIAN_HERMITEGAUSSIAN1D(I,E,R(1),GX)
              CALL GAUSSIAN_HERMITEGAUSSIAN1D(J,E,R(2),GY)
              CALL GAUSSIAN_HERMITEGAUSSIAN1D(K,E,R(3),GZ)
            ELSE
              CALL ERROR$MSG('ID NOT RECOGNIZED')
              CALL ERROR$CHVAL('ID',ID)
              CALL ERROR$STOP('GAUSSIAN_3DORB')
            END IF
            G=G+GX*GY*GZ*C(IND)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_CARTESIANGAUSSIAN1D(N,E,X,G)
!     **************************************************************************
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: E  !EXPONENT
      REAL(8)   ,INTENT(IN) :: X  !POSITION
      REAL(8)   ,INTENT(OUT):: G
!     **************************************************************************
      G=EXP(-E*X**2)*X**N
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_HERMITEGAUSSIAN1D(N,E,X,G)
!     **************************************************************************
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: E  !EXPONENT
      REAL(8)   ,INTENT(IN) :: X  !POSITION
      REAL(8)   ,INTENT(OUT):: G
      REAL(8)               :: G0,GM,GP
      INTEGER(4)            :: J
!     **************************************************************************
      GM=0.D0
      G0=EXP(-E*X*X)
      DO J=1,N
        GP=2.D0*E*(X*G0-REAL(J-1)*GM)
        GM=G0
        G0=GP
      ENDDO
      G=G0
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_FITGAUSS(GID,NR,W,L,F,NEP,NPOW,EP,C)
!     **************************************************************************
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: GID   ! GRID ID
      INTEGER(4),INTENT(IN) :: NR    ! #(RADIAL GRID POINTS)
      INTEGER(4),INTENT(IN) :: L     ! MAIN ANGULAR MOMENTUM
      REAL(8)   ,INTENT(IN) :: W(NR) ! WEIGHTING FOR FIT
      REAL(8)   ,INTENT(IN) :: F(NR) ! FUNCTION TO BE FITTED
      INTEGER(4),INTENT(IN) :: NEP   ! #(EXPONENTS)
      INTEGER(4),INTENT(IN) :: NPOW  ! #(POWERS) HIGHEST POWER=L+2*(NPOW-1)
      REAL(8)   ,INTENT(IN) :: EP(NEP)
      REAL(8)   ,INTENT(OUT):: C(NPOW,NEP)
      REAL(8)               :: R2(NR),RL(NR)
      REAL(8)               :: B(NPOW,NEP)
      REAL(8)               :: A(NPOW,NEP,NPOW,NEP)
      REAL(8)               :: SVAR,AUX(NR)
      REAL(8)               :: Q(2*NPOW,NEP,NEP)
      REAL(8)               :: SCALE(NPOW,NEP)
      REAL(8)               :: WR2(NR)
      INTEGER(4)            :: I,J,I1,J1,I2,J2
      LOGICAL(4),PARAMETER  :: TTEST=.FALSE.
!     **************************************************************************
      CALL RADIAL$R(GID,NR,R2)
      RL(:)=R2(:)**L   ! R2 IS STILL R AND WILL BE SQUARED ONLY IN THE NEXT LINE
      R2(:)=R2(:)**2
      WR2(:)=W(:) !*R2(:)
!
!     ==========================================================================
!     ==  CONSTRUCT B=<F|W|G_I>                                               ==
!     ==========================================================================
      DO I=1,NEP
        AUX(:)=RL(:)*WR2(:)*EXP(-EP(I)*R2(:))*F(:)
        CALL RADIAL$INTEGRAL(GID,NR,AUX,B(1,I))
        DO J=2,NPOW
          AUX(:)=AUX(:)*R2(:)  
          CALL RADIAL$INTEGRAL(GID,NR,AUX,B(J,I))
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  CONSTRUCT B=<GI|W|GJ>                                               ==
!     ==========================================================================
      DO I1=1,NEP
        DO I2=I1,NEP
          AUX(:)=RL(:)**2*WR2(:)*EXP(-(EP(I1)+EP(I2))*R2(:))
          CALL RADIAL$INTEGRAL(GID,NR,AUX,Q(1,I1,I2))
          DO J=2,2*NPOW
            AUX(:)=AUX(:)*R2(:)  
            CALL RADIAL$INTEGRAL(GID,NR,AUX,Q(J,I1,I2))
          ENDDO
          IF(I1.NE.I2) Q(:,I2,I1)=Q(:,I1,I2)
        ENDDO
      ENDDO 
!
      A(:,:,:,:)=0.D0
      DO I2=1,NEP
        DO J2=1,NPOW
          DO I1=1,NEP
            DO J1=1,NPOW
              A(J1,I1,J2,I2)=Q(J1+J2-1,I1,I2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  SCALE TEST FUNCTIONS TO AVOID NUMERICAL PROBLEMS
!     ==========================================================================
      DO I=1,NEP
        SVAR=EXP(1.D0)/(2.D0*EP(I))
        DO J=1,NPOW
          SCALE(J,I)=1.D0/(SVAR*REAL(L+2*J-2))**(L+2*J-2)
        ENDDO
      ENDDO
      DO I2=1,NEP
        DO J2=1,NPOW
          DO I1=1,NEP
            DO J1=1,NPOW
              A(J1,I1,J2,I2)=SCALE(J1,I1)*A(J1,I1,J2,I2)*SCALE(J2,I2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO I1=1,NEP
        DO J1=1,NPOW
          B(J1,I1)=SCALE(J1,I1)*B(J1,I1)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  SOLVE LEAST SQUARES EQUATION                                        ==
!     ==========================================================================
      CALL LIB$MATRIXSOLVER8(NEP*NPOW,NEP*NPOW,1,A,C,B)
!
!     ==========================================================================
!     ==  SCALE COEFFICIENTS
!     ==========================================================================
      DO I1=1,NEP
        DO J1=1,NPOW
          C(J1,I1)=SCALE(J1,I1)*C(J1,I1)
        ENDDO
      ENDDO
!
!     ==========================================================================
!     ==  TEST
!     ==========================================================================
      IF(TTEST) THEN
        DO I1=1,NPOW
          DO J1=1,NEP
            SVAR=0.D0
            DO I2=1,NPOW
              DO J2=1,NEP
                SVAR=SVAR+A(I1,J1,I2,J2)*C(I2,J2)/SCALE(I2,J2)
              ENDDO
            ENDDO
            SVAR=SVAR-B(I1,J1)
            IF(ABS(SVAR).GT.1.D-6) THEN
              PRINT*,'FITTING ERROR ',SVAR,I1,J1
            END IF
          ENDDO
        ENDDO
        AUX(:)=WR2(:)*F(:)**2
        CALL RADIAL$INTEGRAL(GID,NR,AUX,SVAR)
        DO I1=1,NPOW
          DO J1=1,NEP
            DO I2=1,NPOW
              DO J2=1,NEP
                SVAR=SVAR+C(I1,J1)*A(I1,J1,I2,J2)*C(I2,J2) &
      &                  /SCALE(I1,J1)/SCALE(I2,J2)
              ENDDO
            ENDDO
            SVAR=SVAR-2.D0*B(I1,J1)/SCALE(I1,J1)*C(I1,J1)
          ENDDO
        ENDDO
PRINT*,'SQUARE DEVIATION ',SVAR
      END IF
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_CONTRACT2GAUSS(NIJKA,EA,RA,CA,NIJKB,EB,RB,CB &
     &                           ,NIJKP,EP,RP,CP)
!     **************************************************************************
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NIJKA
      REAL(8)   ,INTENT(IN) :: EA
      REAL(8)   ,INTENT(IN) :: RA(3)
      REAL(8)   ,INTENT(IN) :: CA(NIJKA)
      INTEGER(4),INTENT(IN) :: NIJKB
      REAL(8)   ,INTENT(IN) :: EB
      REAL(8)   ,INTENT(IN) :: RB(3)
      REAL(8)   ,INTENT(IN) :: CB(NIJKB)
      INTEGER(4),INTENT(IN) :: NIJKP
      REAL(8)   ,INTENT(OUT):: EP
      REAL(8)   ,INTENT(OUT):: RP(3)
      REAL(8)   ,INTENT(OUT):: CP(NIJKP)
      INTEGER(4)            :: NA,NB,NP
      INTEGER(4)            :: NABX
      INTEGER(4)            :: J,K
      REAL(8)   ,ALLOCATABLE:: HX(:,:,:),HY(:,:,:),HZ(:,:,:)
      INTEGER(4)            :: INDA,MA,IA,JA,KA
      INTEGER(4)            :: INDB,MB,IB,JB,KB
      INTEGER(4)            :: INDP,MP,IP,JP,KP
      INTEGER(4)            :: NIJKA1,NIJKB1,NIJKP1
      REAL(8)               :: SVAR,SVAR1
!     **************************************************************************
!     ==========================================================================
!     ==  TESTS                                                               ==
!     ==========================================================================
      NIJKA1=NIJKA  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKA1,NA,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$I4VAL('NIJKA',NIJKA)
        CALL ERROR$STOP('GAUSSIAN_CONTRACTHGAUSS')
      END IF
      NIJKB1=NIJKB  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKB1,NB,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$I4VAL('NIJKB',NIJKB)
        CALL ERROR$STOP('GAUSSIAN_CONTRACTHGAUSS')
      END IF
      NIJKP1=NIJKP  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKP1,NP,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$I4VAL('NIJKP',NIJKP)
        CALL ERROR$STOP('GAUSSIAN_CONTRACTHGAUSS')
      END IF
!
      EP=EA+EB
      RP=(RA*EA+RB*EB)/EP
      NABX=MAX(NA,NB)
      ALLOCATE(HX(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HY(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HZ(0:NABX,0:NABX,0:2*NABX))
      CALL GAUSSIAN_HERMITEC(NABX,RA(1),RB(1),EA,EB,HX)
      CALL GAUSSIAN_HERMITEC(NABX,RA(2),RB(2),EA,EB,HY)
      CALL GAUSSIAN_HERMITEC(NABX,RA(3),RB(3),EA,EB,HZ)

      INDP=0
      DO MP=0,NP
        DO IP=0,MP
          DO JP=0,MP-IP
            KP=MP-IP-JP
            INDP=INDP+1
!
            SVAR=0.D0
            INDB=0
            DO MB=0,NB
              DO IB=0,MB
                DO JB=0,MB-IB
                  KB=MB-IB-JB
                  INDB=INDB+1
!
                  INDA=0
                  DO MA=0,NA
                    DO IA=0,MA
                      SVAR1=0.D0
                      DO JA=0,MA-IA
                        KA=MA-IA-JA
                        INDA=INDA+1
                        SVAR1=SVAR1+CA(INDA)*HY(JA,JB,JP)*HZ(KA,KB,KP)
                      ENDDO
                      SVAR=SVAR+CB(INDB)*HX(IA,IB,IP)*SVAR1
                    ENDDO
                  ENDDO

                ENDDO
              ENDDO
            ENDDO
            CP(INDP)=SVAR
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_CONTRACT2GAUSSMAT(NIJKA,NEA,EA,RA,NIJKB,NEB,EB,RB &
     &                                     ,NIJKP,NEP,EP,RP,CPAB)
!     **************************************************************************
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NIJKA
      INTEGER(4),INTENT(IN) :: NEA
      REAL(8)   ,INTENT(IN) :: EA(NEA)
      REAL(8)   ,INTENT(IN) :: RA(3)
      INTEGER(4),INTENT(IN) :: NIJKB
      INTEGER(4),INTENT(IN) :: NEB
      REAL(8)   ,INTENT(IN) :: EB(NEB)
      REAL(8)   ,INTENT(IN) :: RB(3)
      INTEGER(4),INTENT(IN) :: NIJKP
      INTEGER(4),INTENT(IN) :: NEP
      REAL(8)   ,INTENT(OUT):: EP(NEP)
      REAL(8)   ,INTENT(OUT):: RP(3,NEP)
      REAL(8)   ,INTENT(OUT):: CPAB(NIJKA,NIJKB,NIJKP,NEA,NEB)
      INTEGER(4)            :: NA,NB,NP
      INTEGER(4)            :: NABX
      INTEGER(4)            :: J,K
      REAL(8)   ,ALLOCATABLE:: HX(:,:,:),HY(:,:,:),HZ(:,:,:)
      INTEGER(4)            :: INDA,MA,IA,JA,KA,IEA
      INTEGER(4)            :: INDB,MB,IB,JB,KB,IEB
      INTEGER(4)            :: INDP,MP,IP,JP,KP,IEP
      INTEGER(4)            :: NIJKA1,NIJKB1,NIJKP1
      REAL(8)               :: SVAR,SVAR1
!     **************************************************************************
!     ==========================================================================
!     ==  TESTS                                                               ==
!     ==========================================================================
      NIJKA1=NIJKA  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKA1,NA,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$I4VAL('NIJKA',NIJKA)
        CALL ERROR$STOP('GAUSSIAN_CONTRACTHGAUSS')
      END IF
      NIJKB1=NIJKB  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKB1,NB,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$I4VAL('NIJKB',NIJKB)
        CALL ERROR$STOP('GAUSSIAN_CONTRACTHGAUSS')
      END IF
      NIJKP1=NIJKP  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKP1,NP,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$I4VAL('NIJKP',NIJKP)
        CALL ERROR$STOP('GAUSSIAN_CONTRACTHGAUSS')
      END IF
!
      NABX=MAX(NA,NB)
      ALLOCATE(HX(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HY(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HZ(0:NABX,0:NABX,0:2*NABX))
      IEP=0
      DO IEA=1,NEA
        DO IEB=1,NEB
          IEP=IEP+1
          EP(IEP)=EA(IEA)+EB(IEB)
          RP(:,IEP)=(RA(:)*EA(IEA)+RB(:)*EB(IEB))/EP(IEP)
          CALL GAUSSIAN_HERMITEC(NABX,RA(1),RB(1),EA(IEA),EB(IEB),HX)
          CALL GAUSSIAN_HERMITEC(NABX,RA(2),RB(2),EA(IEA),EB(IEB),HY)
          CALL GAUSSIAN_HERMITEC(NABX,RA(3),RB(3),EA(IEA),EB(IEB),HZ)

          INDP=0
          DO MP=0,NP
            DO IP=0,MP
              DO JP=0,MP-IP
                KP=MP-IP-JP
                INDP=INDP+1
!         
                SVAR=0.D0
                INDB=0
                DO MB=0,NB
                  DO IB=0,MB
                    DO JB=0,MB-IB
                      KB=MB-IB-JB
                      INDB=INDB+1
!         
                      INDA=0
                      DO MA=0,NA
                        DO IA=0,MA
                          SVAR1=0.D0
                          DO JA=0,MA-IA
                            KA=MA-IA-JA
                            INDA=INDA+1
                            CPAB(INDA,INDB,INDP,IEA,IEB)=HX(IA,IB,IP) &
         &                                              *HY(JA,JB,JP) &
         &                                              *HZ(KA,KB,KP)
                          ENDDO
                        ENDDO
                      ENDDO
!
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN$OVERLAP(NIJKA,NCA,EA,RA,CA,NIJKB,NCB,EB,RB,CB,SAB)
!     **************************************************************************
!     ** EVALUATES THE OVERLAP OF TWO SETS OF FUNCTIONS REPRESENTED BY        **
!     ** CARTESIAN GAUSSIANS.                                                 **
!     ** SET "A" WITH NCA FUNCTIONS HAS THE EXPONENT EA AND IS CENTERED AT RA **
!     ** SET "B" WITH NCB FUNCTIONS HAS THE EXPONENT EB AND IS CENTERED AT RB **
!     ** THE RESULTING OVERLAP MATRIX IS S_IJ=<A_I|B_J>                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NIJKA
      INTEGER(4),INTENT(IN) :: NCA
      REAL(8)   ,INTENT(IN) :: EA
      REAL(8)   ,INTENT(IN) :: RA(3)
      REAL(8)   ,INTENT(IN) :: CA(NIJKA,NCA)
      INTEGER(4),INTENT(IN) :: NIJKB
      INTEGER(4),INTENT(IN) :: NCB
      REAL(8)   ,INTENT(IN) :: EB
      REAL(8)   ,INTENT(IN) :: RB(3)
      REAL(8)   ,INTENT(IN) :: CB(NIJKB,NCB)
      REAL(8)   ,INTENT(OUT) :: SAB(NCA,NCB)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)            :: NA,NB
      INTEGER(4)            :: NABX
      REAL(8)   ,ALLOCATABLE:: HX(:,:,:),HY(:,:,:),HZ(:,:,:)
      REAL(8)   ,ALLOCATABLE:: PRIMOV(:,:)
      INTEGER(4)            :: INDA,MA,IA,JA,KA
      INTEGER(4)            :: INDB,MB,IB,JB,KB
      INTEGER(4)            :: NIJKA1,NIJKB1
      REAL(8)               :: SVAR
!     **************************************************************************
!
!     ==========================================================================
!     ==  TESTS                                                               ==
!     ==========================================================================
      NIJKA1=NIJKA  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKA1,NA,JA,KA)
      IF(JA.NE.0.OR.KA.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_OVERLAP')
      END IF
      NIJKB1=NIJKB  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKB1,NB,JB,KB)
      IF(JB.NE.0.OR.KB.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_OVERLAP')
      END IF
!
!     ==========================================================================
!     ==  DETERMINE HERMITE COEFFICIENTS                                      ==
!     ==========================================================================
      NABX=MAX(NA,NB)
      ALLOCATE(HX(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HY(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HZ(0:NABX,0:NABX,0:2*NABX))
      CALL GAUSSIAN_HERMITEC(NABX,RA(1),RB(1),EA,EB,HX)
      CALL GAUSSIAN_HERMITEC(NABX,RA(2),RB(2),EA,EB,HY)
      CALL GAUSSIAN_HERMITEC(NABX,RA(3),RB(3),EA,EB,HZ)
!
!     ==========================================================================
!     ==  WORK OUT OVERLAP MATRIX                                             ==
!     ==========================================================================
      ALLOCATE(PRIMOV(NIJKA,NIJKB))
      PRIMOV(:,:)=0.D0
      SAB(:,:)=0.D0
      INDA=0
      DO MA=0,NA
        DO IA=0,MA
          DO JA=0,MA-IA
            KA=MA-IA-JA
            INDA=INDA+1
!
            INDB=0
            DO MB=0,NB
              DO IB=0,MB
                DO JB=0,MB-IB
                  KB=MB-IB-JB
                  INDB=INDB+1
!
                  SVAR=HX(IA,IB,0)*HY(JA,JB,0)*HZ(KA,KB,0)
                  PRIMOV(INDA,INDB)=SVAR
!!$                  DO J=1,NCB
!!$                    DO I=1,NCA
!!$                      SAB(I,J)=SAB(I,J)+SVAR*CA(INDA,I)*CB(INDB,J)
!!$                    ENDDO
!!$                  ENDDO
!
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      PRIMOV(:,:)=PRIMOV(:,:)*(PI/(EA+EB))**1.5D0
!      SAB(:,:)=SAB(:,:)*(PI/(EA+EB))**1.5D0
      SAB=MATMUL(TRANSPOSE(CA(:,:)),MATMUL(PRIMOV,CB))
!
!!$DO I=1,NIJKA
!!$  WRITE(*,FMT='("CA",I8,20E10.2)')I,CA(I,:),CB(I,:)
!!$ENDDO
!!$DO J=1,NCB
!!$  WRITE(*,FMT='("S",20E10.2)')SAB(:,J)
!!$ENDDO
!!$STOP 'FORCED IN GAUSSIAN$OVERLAP'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN$NEWOVERLAP(NIJKA,NCA,NEA,EA,RA,CA &
     &                           ,NIJKB,NCB,NEB,EB,RB,CB,SAB)
!     **************************************************************************
!     ** EVALUATES THE OVERLAP OF TWO SETS OF FUNCTIONS REPRESENTED BY        **
!     ** CARTESIAN GAUSSIANS.                                                 **
!     ** SET "A" WITH NCA FUNCTIONS HAS THE EXPONENT EA AND IS CENTERED AT RA **
!     ** SET "B" WITH NCB FUNCTIONS HAS THE EXPONENT EB AND IS CENTERED AT RB **
!     ** THE RESULTING OVERLAP MATRIX IS S_IJ=<A_I|B_J>                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NIJKA
      INTEGER(4),INTENT(IN) :: NCA
      INTEGER(4),INTENT(IN) :: NEA
      REAL(8)   ,INTENT(IN) :: EA(NEA)
      REAL(8)   ,INTENT(IN) :: RA(3)
      REAL(8)   ,INTENT(IN) :: CA(NIJKA,NEA,NCA)
      INTEGER(4),INTENT(IN) :: NIJKB
      INTEGER(4),INTENT(IN) :: NEB
      INTEGER(4),INTENT(IN) :: NCB
      REAL(8)   ,INTENT(IN) :: EB(NEB)
      REAL(8)   ,INTENT(IN) :: RB(3)
      REAL(8)   ,INTENT(IN) :: CB(NIJKB,NEB,NCB)
      REAL(8)   ,INTENT(OUT):: SAB(NCA,NCB)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)            :: NA,NB
      INTEGER(4)            :: NABX
      REAL(8)   ,ALLOCATABLE:: HX(:,:,:),HY(:,:,:),HZ(:,:,:)
      REAL(8)   ,ALLOCATABLE:: PRIMOV(:,:)
      INTEGER(4)            :: INDA,MA,IA,JA,KA,IEA
      INTEGER(4)            :: INDB,MB,IB,JB,KB,IEB
      INTEGER(4)            :: NIJKA1,NIJKB1
!     **************************************************************************
!
!     ==========================================================================
!     ==  TESTS                                                               ==
!     ==========================================================================
      NIJKA1=NIJKA  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKA1,NA,JA,KA)
      IF(JA.NE.0.OR.KA.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_OVERLAP')
      END IF
      NIJKB1=NIJKB  !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKB1,NB,JB,KB)
      IF(JB.NE.0.OR.KB.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_OVERLAP')
      END IF
!
!     ==========================================================================
!     ==  DETERMINE HERMITE COEFFICIENTS                                      ==
!     ==========================================================================
      NABX=MAX(NA,NB)
      ALLOCATE(HX(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HY(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HZ(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(PRIMOV(NIJKA,NIJKB))
      SAB(:,:)=0.D0
      DO IEA=1,NEA
        DO IEB=1,NEB
          CALL GAUSSIAN_HERMITEC(NABX,RA(1),RB(1),EA(IEA),EB(IEB),HX)
          CALL GAUSSIAN_HERMITEC(NABX,RA(2),RB(2),EA(IEA),EB(IEB),HY)
          CALL GAUSSIAN_HERMITEC(NABX,RA(3),RB(3),EA(IEA),EB(IEB),HZ)
!
!         ======================================================================
!         ==  WORK OUT OVERLAP MATRIX                                         ==
!         ======================================================================
          INDA=0
          DO MA=0,NA
            DO IA=0,MA
              DO JA=0,MA-IA
                KA=MA-IA-JA
                INDA=INDA+1
!
                INDB=0
                DO MB=0,NB
                  DO IB=0,MB
                    DO JB=0,MB-IB
                      KB=MB-IB-JB
                      INDB=INDB+1
                      PRIMOV(INDA,INDB)=HX(IA,IB,0)*HY(JA,JB,0)*HZ(KA,KB,0)
                    ENDDO
                  ENDDO
                ENDDO

              ENDDO
            ENDDO
          ENDDO
          PRIMOV(:,:)=PRIMOV(:,:)*(PI/(EA(IEA)+EB(IEB)))**1.5D0
          SAB=SAB+MATMUL(TRANSPOSE(CA(:,IEA,:)),MATMUL(PRIMOV,CB(:,IEB,:)))
        ENDDO
      ENDDO  
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN$GAUSSOVERLAP(NIJKA,NEA,EA,RA,NIJKB,NEB,EB,RB,SAB)
!     **************************************************************************
!     ** EVALUATES THE OVERLAP MATRIX OF TWO SETS OF CARTESIAN GAUSSIANS      **
!     ** SET "A" WITH NEA EXPONENTS EA IS CENTERED AT RA                      **
!     ** SET "B" WITH NEB EXPONENTS EB IS CENTERED AT RB                      **
!     ** THE RESULTING OVERLAP MATRIX IS S_IJ=<G_I|G_J>                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NIJKA
      INTEGER(4),INTENT(IN) :: NEA
      REAL(8)   ,INTENT(IN) :: EA(NEA)
      REAL(8)   ,INTENT(IN) :: RA(3)
      INTEGER(4),INTENT(IN) :: NIJKB
      INTEGER(4),INTENT(IN) :: NEB
      REAL(8)   ,INTENT(IN) :: EB(NEB)
      REAL(8)   ,INTENT(IN) :: RB(3)
      REAL(8)   ,INTENT(OUT):: SAB(NIJKA,NEA,NIJKB,NEB)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)            :: NA,NB
      INTEGER(4)            :: NABX
      REAL(8)   ,ALLOCATABLE:: HX(:,:,:),HY(:,:,:),HZ(:,:,:)
      INTEGER(4)            :: INDA,MA,IA,JA,KA,IEA
      INTEGER(4)            :: INDB,MB,IB,JB,KB,IEB
      INTEGER(4)            :: NIJKA1,NIJKB1
!     **************************************************************************
!
!     ==========================================================================
!     ==  TESTS                                                               ==
!     ==========================================================================
      NIJKA1=NIJKA !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKA1,NA,JA,KA)
      IF(JA.NE.0.OR.KA.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_GAUSSOVERLAP')
      END IF
      NIJKB1=NIJKB !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKB1,NB,JB,KB)
      IF(JB.NE.0.OR.KB.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_GAUSSOVERLAP')
      END IF
!
!     ==========================================================================
!     ==  DETERMINE HERMITE COEFFICIENTS                                      ==
!     ==========================================================================
      NABX=MAX(NA,NB)
      ALLOCATE(HX(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HY(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HZ(0:NABX,0:NABX,0:2*NABX))
      SAB(:,:,:,:)=0.D0
      DO IEA=1,NEA
        DO IEB=1,NEB
          CALL GAUSSIAN_HERMITEC(NABX,RA(1),RB(1),EA(IEA),EB(IEB),HX)
          CALL GAUSSIAN_HERMITEC(NABX,RA(2),RB(2),EA(IEA),EB(IEB),HY)
          CALL GAUSSIAN_HERMITEC(NABX,RA(3),RB(3),EA(IEA),EB(IEB),HZ)
!
!         ======================================================================
!         ==  WORK OUT OVERLAP MATRIX                                         ==
!         ======================================================================
          INDA=0
          DO MA=0,NA
            DO IA=0,MA
              DO JA=0,MA-IA
                KA=MA-IA-JA
                INDA=INDA+1
!
                INDB=0
                DO MB=0,NB
                  DO IB=0,MB
                    DO JB=0,MB-IB
                      KB=MB-IB-JB
                      INDB=INDB+1
                      SAB(INDA,IEA,INDB,IEB)=HX(IA,IB,0)*HY(JA,JB,0)*HZ(KA,KB,0)
                    ENDDO
                  ENDDO
                ENDDO

              ENDDO
            ENDDO
          ENDDO
          SAB(:,IEA,:,IEB)=SAB(:,IEA,:,IEB)*(PI/(EA(IEA)+EB(IEB)))**1.5D0
        ENDDO
      ENDDO  
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN$GAUSSPSI(NIJKA,NEA,EA,RA,NIJKB,NEB,EB,RB,NCB,CB,GCB)
!     **************************************************************************
!     ** EVALUATES THE MATRIX ELEMENTS <G(I,IE)|PSI_N> OF A FUNCTION |PSI_N>  **
!     ** AT RB WITH GAUSSIANS CENTERED AT RA                                  **
!     ** SET "A" WITH NCA FUNCTIONS HAS THE EXPONENT EA AND IS CENTERED AT RA **
!     ** SET "B" WITH NCB FUNCTIONS HAS THE EXPONENT EB AND IS CENTERED AT RB **
!     ** THE RESULTING OVERLAP MATRIX IS S_IJ=<A_I|B_J>                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NIJKA
      INTEGER(4),INTENT(IN) :: NEA
      REAL(8)   ,INTENT(IN) :: EA(NEA)
      REAL(8)   ,INTENT(IN) :: RA(3)
      INTEGER(4),INTENT(IN) :: NIJKB
      INTEGER(4),INTENT(IN) :: NEB
      REAL(8)   ,INTENT(IN) :: EB(NEB)
      REAL(8)   ,INTENT(IN) :: RB(3)
      INTEGER(4),INTENT(IN) :: NCB
      REAL(8)   ,INTENT(IN) :: CB(NIJKB,NEB,NCB)   ! |PSI_N>=SUM_J:|GB_J>CB(J,N)
      REAL(8)   ,INTENT(OUT):: GCB(NIJKA,NEA,NCB)  ! <GA_I|PSI_N>
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)            :: NA,NB
      INTEGER(4)            :: NABX
      INTEGER(4)            :: J,K
      REAL(8)   ,ALLOCATABLE:: HX(:,:,:),HY(:,:,:),HZ(:,:,:)
      INTEGER(4)            :: INDA,MA,IA,JA,KA,IEA
      INTEGER(4)            :: INDB,MB,IB,JB,KB,IEB
      REAL(8)               :: SVAR,FAC
      INTEGER(4)            :: NIJKA1,NIJKB1
!     **************************************************************************
      GCB(:,:,:)=0.D0
!
!     ==========================================================================
!     ==  TESTS                                                               ==
!     ==========================================================================
      NIJKA1=NIJKA !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKA1,NA,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_OVERLAP')
      END IF
      NIJKB1=NIJKB !AVOID MAPPING INTENT(IN) ONTO INTENT(INOUT)
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKB1,NB,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN_OVERLAP')
      END IF
!
!     ==========================================================================
!     ==  DETERMINE HERMITE COEFFICIENTS                                      ==
!     ==========================================================================
      NABX=MAX(NA,NB)
      ALLOCATE(HX(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HY(0:NABX,0:NABX,0:2*NABX))
      ALLOCATE(HZ(0:NABX,0:NABX,0:2*NABX))
      DO IEA=1,NEA
        DO IEB=1,NEB
          CALL GAUSSIAN_HERMITEC(NABX,RA(1),RB(1),EA(IEA),EB(IEB),HX)
          CALL GAUSSIAN_HERMITEC(NABX,RA(2),RB(2),EA(IEA),EB(IEB),HY)
          CALL GAUSSIAN_HERMITEC(NABX,RA(3),RB(3),EA(IEA),EB(IEB),HZ)
!
!         ======================================================================
!         ==  WORK OUT OVERLAP <G|PSI_N>                                      ==
!         ======================================================================
          FAC=(PI/(EA(IEA)+EB(IEB)))**1.5D0
          INDA=0
          DO MA=0,NA
            DO IA=0,MA
              DO JA=0,MA-IA
                KA=MA-IA-JA
                INDA=INDA+1
!
                INDB=0
                DO MB=0,NB
                  DO IB=0,MB
                    DO JB=0,MB-IB
                      KB=MB-IB-JB
                      INDB=INDB+1
                      SVAR=FAC*HX(IA,IB,0)*HY(JA,JB,0)*HZ(KA,KB,0)
                      GCB(INDA,IEA,:)=GCB(INDA,IEA,:)+SVAR*CB(INDB,IEB,:)
                    ENDDO
                  ENDDO
                ENDDO

              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO  
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_HERMITEC(N,XA,XB,EA,EB,H)
!     **************************************************************************
!     ** EVALUATES A SET OF HERMITE COEFFICIENTS                              **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: XA
      REAL(8)   ,INTENT(IN) :: EA
      REAL(8)   ,INTENT(IN) :: XB
      REAL(8)   ,INTENT(IN) :: EB
      REAL(8)   ,INTENT(OUT):: H(0:N,0:N,0:2*N)
      REAL(8)               :: EP,MU,ONEBY2P,XP,PA,PB,TP1
      INTEGER(4)            :: I,J,IJSUM,T,I1,I2
!     **************************************************************************
      IF(N.LT.0) THEN
        CALL ERROR$MSG('N MUST BE A NON-NEGATIVE INTEGER')
        CALL ERROR$I4VAL('N',N)
        CALL ERROR$STOP('GAUSSIAN_HERMITEC')
      END IF
      H(:,:,:)=0.D0
      EP=EA+EB
      MU=EA*EB/EP
      ONEBY2P=1.D0/(2.D0*EP)
      XP=(EA*XA+EB*XB)/EP
      PA=XP-XA
      PB=XP-XB
      H(0,0,0)=EXP(-MU*(XB-XA)**2)
      IF(N.LT.1) RETURN
!     == INITIALIZE UP TO I+J=1 TO AVOID REFERRING TO INDICES OUT OF RANGE
      H(1,0,0)=PA*H(0,0,0)
      H(1,0,1)=ONEBY2P*H(0,0,0)
      H(0,1,0)=PB*H(0,0,0)
      H(0,1,1)=ONEBY2P*H(0,0,0)
      DO IJSUM=2,2*N
        IF(IJSUM.LE.N) THEN
          TP1=1.D0
          H(0,IJSUM,0)=PB*H(0,IJSUM-1,0)+TP1*H(0,IJSUM-1,1)
          DO T=1,IJSUM-2
            TP1=REAL(T+1)
            H(0,IJSUM,T)=ONEBY2P*H(0,IJSUM-1,T-1) &
     &                 +PB*H(0,IJSUM-1,T)+TP1*H(0,IJSUM-1,T+1)
          ENDDO
          H(0,IJSUM,IJSUM-1)=ONEBY2P*H(0,IJSUM-1,IJSUM-2) &
     &                       +PB*H(0,IJSUM-1,IJSUM-1)
          H(0,IJSUM,IJSUM)=ONEBY2P*H(0,IJSUM-1,IJSUM-1)
        END IF
!
        I1=MAX(1,IJSUM-N)
        I2=MIN(IJSUM,N)
        DO I=I1,I2
          J=IJSUM-I
!         == DO (I,J,0) SEPARATELY TO AVOID REFERRING TO (I-1,J,-1)
          TP1=1.D0
          H(I,J,0)=PA*H(I-1,J,0)+TP1*H(I-1,J,1)
!         == NOW NORMAL EXPRESSIONS
          DO T=1,IJSUM-2
            TP1=REAL(T+1)
            H(I,J,T)=ONEBY2P*H(I-1,J,T-1)+PA*H(I-1,J,T)+TP1*H(I-1,J,T+1)
          ENDDO
          H(I,J,IJSUM-1)=ONEBY2P*H(I-1,J,IJSUM-2)+PA*H(I-1,J,IJSUM-1)
          H(I,J,IJSUM)=ONEBY2P*H(I-1,J,IJSUM-1)
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN$FOURCENTER(NIJKA,NEA,NORBA,EA,RA,CA &
     &                              ,NIJKB,NEB,NORBB,EB,RB,CB &
     &                              ,NIJKC,NEC,NORBC,EC,RC,CC &
     &                              ,NIJKD,NED,NORBD,ED,RD,CD,UABCD)
!     **************************************************************************
!     ** 
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NIJKA
      INTEGER(4),INTENT(IN) :: NIJKB
      INTEGER(4),INTENT(IN) :: NIJKC
      INTEGER(4),INTENT(IN) :: NIJKD
      INTEGER(4),INTENT(IN) :: NEA
      INTEGER(4),INTENT(IN) :: NEB
      INTEGER(4),INTENT(IN) :: NEC
      INTEGER(4),INTENT(IN) :: NED
      INTEGER(4),INTENT(IN) :: NORBA
      INTEGER(4),INTENT(IN) :: NORBB
      INTEGER(4),INTENT(IN) :: NORBC
      INTEGER(4),INTENT(IN) :: NORBD
      REAL(8)   ,INTENT(IN) :: EA(NEA)
      REAL(8)   ,INTENT(IN) :: EB(NEB)
      REAL(8)   ,INTENT(IN) :: EC(NEC)
      REAL(8)   ,INTENT(IN) :: ED(NED)
      REAL(8)   ,INTENT(IN) :: RA(3)
      REAL(8)   ,INTENT(IN) :: RB(3)
      REAL(8)   ,INTENT(IN) :: RC(3)
      REAL(8)   ,INTENT(IN) :: RD(3)
      REAL(8)   ,INTENT(IN) :: CA(NIJKA,NEA,NORBA)
      REAL(8)   ,INTENT(IN) :: CB(NIJKB,NEB,NORBB)
      REAL(8)   ,INTENT(IN) :: CC(NIJKC,NEC,NORBC)
      REAL(8)   ,INTENT(IN) :: CD(NIJKD,NED,NORBD)
      REAL(8)   ,INTENT(OUT):: UABCD(NORBA,NORBB,NORBC,NORBD)
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)             :: NIJKX !TEMP.COPY OF NIJKA ETC.
      INTEGER(4)             :: NIJKP,NIJKQ
      REAL(8)                :: EP,EQ
      REAL(8)                :: RP(3),RQ(3)
      REAL(8)   ,ALLOCATABLE :: CP(:),CQ(:)
      INTEGER(4)             :: NA,NB,NC,ND,NP,NQ,N
      INTEGER(4)             :: J,K
      REAL(8)   ,ALLOCATABLE :: HI(:,:,:)
      REAL(8)                :: UABCD1
      INTEGER(4)             :: INDP,INDQ,MP,MQ,IP,IQ,JP,JQ,KP,KQ
      INTEGER(4)             :: IORBA,IORBB,IORBC,IORBD
      INTEGER(4)             :: IEA,IEB,IEC,IED
      REAL(8)                :: SGN
!     **************************************************************************
      NIJKX=NIJKA  ! COPY TO AVOID INTENT(IN/INOUT) CONFUSION
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKX,NA,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN$FOURCENTER')
      END IF
      NIJKX=NIJKB  ! COPY TO AVOID INTENT(IN/INOUT) CONFUSION
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKX,NB,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN$FOURCENTER')
      END IF
      NIJKX=NIJKC  ! COPY TO AVOID INTENT(IN/INOUT) CONFUSION
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKX,NC,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN$FOURCENTER')
      END IF
      NIJKX=NIJKD  ! COPY TO AVOID INTENT(IN/INOUT) CONFUSION
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKX,ND,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN$FOURCENTER')
      END IF
!
      NP=NA+NB 
      NQ=NC+ND 
      J=0        !AVOID INTENT(IN/INOUT) CONFUSION
      K=0        !AVOID INTENT(IN/INOUT) CONFUSION
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKP,NP,J,K)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKQ,NQ,J,K)
      ALLOCATE(CP(NIJKP))
      ALLOCATE(CQ(NIJKQ))
      N=NP+NQ
      ALLOCATE(HI(0:N,0:N,0:N))
      UABCD(:,:,:,:)=0.D0
      DO IORBA=1,NORBA
        DO IORBB=1,NORBB
          DO IEA=1,NEA
            DO IEB=1,NEB
      CALL GAUSSIAN_CONTRACT2GAUSS(NIJKA,EA(IEA),RA,CA(:,IEA,IORBA) &
     &                            ,NIJKB,EB(IEB),RB,CB(:,IEB,IORBB) &
     &                            ,NIJKP,EP,RP,CP)
PRINT*,'MARKE 1',IORBA,IORBB,IEA,IEB
              DO IORBC=1,NORBC
                DO IORBD=1,NORBD
                  DO IEC=1,NEC
                    DO IED=1,NED
      CALL GAUSSIAN_CONTRACT2GAUSS(NIJKC,EC(IEC),RC,CC(:,IEC,IORBC) &
     &                            ,NIJKD,ED(IED),RD,CD(:,IED,IORBD) &
     &                            ,NIJKQ,EQ,RQ,CQ)
      CALL GAUSSIAN_HERMITEINTEGRAL(N,EP*EQ/(EP+EQ),RP-RQ,HI)
      UABCD1=0.D0
      INDP=0
      DO MP=0,NP      
        DO IP=0,MP
          DO JP=0,MP-IP
            KP=MP-IP-JP
            INDP=INDP+1
            IF(ABS(CP(INDP)).LT.1.D-6) CYCLE
!         
            INDQ=0
            SGN=-1.D0
            DO MQ=0,NQ      
              SGN=-SGN     !SGN=(-1)^(IQ+JQ+KQ)=(-1)^MQ
              DO IQ=0,MQ
                DO JQ=0,MQ-IQ
                  KQ=MQ-IQ-JQ
                  INDQ=INDQ+1
                  UABCD1=UABCD1+SGN*HI(IP+IQ,JP+JQ,KP+KQ) &
     &                            *CP(INDP)*CQ(INDQ)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      UABCD1=2.D0*PI**2.5D0/(EP*EQ*SQRT(EP+EQ))*UABCD1
      UABCD(IORBA,IORBB,IORBC,IORBD)=UABCD(IORBA,IORBB,IORBC,IORBD)+UABCD1
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN$ZDIRECTION_FOURCENTER(NIJKA,NEA,EA,NFA,CA &
     &                                         ,NIJKB,NEB,EB,NFB,CB,DIS,UABCD)
!     **************************************************************************
!     ** CALCULATES FOUR-CENTER MATRIX ELEMENTS FOR CARTESIAN GAUSSIANS       **
!     ** LOCATED AT THE ENDS OF A SINGLE BOND ORIENTED IN Z-DIRECTION         **
!     ** U(1,2,3,4)=INT DX INT DX': A1(X)*B2(X)][A3(X')*B4(X')]/|R-R'|        **
!     ** THE ORDER OF INDICES DEVIATES FROM THE CONVENTION FOR THE UTENSOR!   **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)  :: NIJKA
      INTEGER(4),INTENT(IN)  :: NIJKB
      INTEGER(4),INTENT(IN)  :: NEA
      INTEGER(4),INTENT(IN)  :: NEB
      REAL(8)   ,INTENT(IN)  :: EA(NEA)
      REAL(8)   ,INTENT(IN)  :: EB(NEB)
      INTEGER(4),INTENT(IN)  :: NFA
      INTEGER(4),INTENT(IN)  :: NFB
      REAL(8)   ,INTENT(IN)  :: CA(NIJKA,NEA,NFA)
      REAL(8)   ,INTENT(IN)  :: CB(NIJKB,NEB,NFB)
      REAL(8)   ,INTENT(IN)  :: DIS
      REAL(8)   ,INTENT(OUT) :: UABCD(NFA,NFB,NFA,NFB)
      REAL(8)   ,PARAMETER   :: PI=4.D0*ATAN(1.D0)
      INTEGER(4)             :: NIJKX  !COPY OF NIJKA,ETC
      INTEGER(4)             :: NIJKP
      INTEGER(4)             :: NEP
      REAL(8)                :: RA(3)
      REAL(8)                :: RB(3)
      REAL(8)                :: EP(NEA*NEB)
      REAL(8)                :: RP(3,NEA*NEB)
      INTEGER(4)             :: NA,NB,NP,N
      INTEGER(4)             :: J,K
      REAL(8)   ,ALLOCATABLE :: HI(:,:,:)
      REAL(8)   ,ALLOCATABLE :: CPAB(:,:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: CPAB1(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: CPAB2(:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: CPAB3(:,:,:,:,:)
      REAL(8)   ,ALLOCATABLE :: UPQ(:,:,:,:)
      REAL(8)                :: FACTOR,SVAR
      INTEGER(4)             :: INDA,INDB,INDP,INDQ
      INTEGER(4)             :: MP,MQ,IP,IQ,JP,JQ,KP,KQ
      INTEGER(4)             :: IEA,IEB,IEP,IEQ
      INTEGER(4)             :: IFA,IFB,IFC,IFD
      REAL(8)                :: SGN
!     **************************************************************************
                     CALL TRACE$PUSH('GAUSSIAN$ZDIRECTION_FOURCENTER')
      NIJKX=NIJKA
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKX,NA,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN$FOURCENTER')
      END IF
      NIJKX=NIJKB
      CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',NIJKX,NB,J,K)
      IF(J.NE.0.OR.K.NE.0) THEN
        CALL ERROR$MSG('COEFFICIENT ARRAY DOES NOT COVER COMPLETE SHELLS')
        CALL ERROR$STOP('GAUSSIAN$FOURCENTER')
      END IF
!
!     ==========================================================================
!     ==  CONTRACT PAIRS OF CARTESIAN GAUSSIANS INTO HERMITE GAUSSIANS        ==
!     ==========================================================================
      RA(:)=0.D0
      RB(:)=(/0.D0,0.D0,DIS/)
      NP=NA+NB 
      J=0
      K=0
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKP,NP,J,K)
      NEP=NEA*NEB
      ALLOCATE(CPAB(NIJKA,NIJKB,NIJKP,NEA,NEB))
PRINT*,'IN GAUSSIAN$ZDIRECTION_FOURCENTER: BEFORE CONTRACT2GAUSSMAT'
      CALL GAUSSIAN_CONTRACT2GAUSSMAT(NIJKA,NEA,EA,RA,NIJKB,NEB,EB,RB &
     &                               ,NIJKP,NEP,EP,RP,CPAB)
!
!     ==========================================================================
!     ==  CONTRACT INTO ORBITALS
!     ==========================================================================
PRINT*,'IN GAUSSIAN$ZDIRECTION_FOURCENTER: BEFORE FOLD INTO ORBITALS LOOP 1'
      ALLOCATE(CPAB3(NIJKP,NIJKA,NEA,NEB,NFB))
      CPAB3=0.D0
      DO IFB=1,NFB
        DO IEA=1,NEA
          DO INDA=1,NIJKA
            DO IEB=1,NEB
              DO INDB=1,NIJKB
                CPAB3(:,INDA,IEA,IEB,IFB)=CPAB3(:,INDA,IEA,IEB,IFB) &
    &                                    +CPAB(INDA,INDB,:,IEA,IEB) &
    &                                    *CB(INDB,IEB,IFB)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(CPAB)
PRINT*,'IN GAUSSIAN$ZDIRECTION_FOURCENTER: BEFORE FOLD INTO ORBITALS LOOP 2'
      ALLOCATE(CPAB1(NIJKP,NEP,NFA,NFB))
      CPAB1=0.D0
      DO IFA=1,NFA
        DO IFB=1,NFB
          IEP=0
          DO IEA=1,NEA
            DO IEB=1,NEB
              IEP=IEP+1
              DO INDA=1,NIJKA
                CPAB1(:,IEP,IFA,IFB)=CPAB1(:,IEP,IFA,IFB) &
    &                               +CPAB3(:,INDA,IEA,IEB,IFB) &
    &                               *CA(INDA,IEA,IFA)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(CPAB3)
!
!     ==========================================================================
!     ==  COULOMB INTEGRAL OF HERMITE GAUSSIANS                               ==
!     ==========================================================================
PRINT*,'IN GAUSSIAN$ZDIRECTION_FOURCENTER: BEFORE HERMITEINTEGRALS ',NP
      N=2*NP
      ALLOCATE(HI(0:N,0:N,0:N))
      ALLOCATE(UPQ(NIJKP,NIJKP,NEP,NEP))
      DO IEP=1,NEP
        DO IEQ=1,NEP
          CALL GAUSSIAN_HERMITEINTEGRAL(N,EP(IEP)*EP(IEQ)/(EP(IEP)+EP(IEQ)) &
     &                                   ,RP(:,IEP)-RP(:,IEQ),HI)
          FACTOR=2.D0*PI**2.5D0/(EP(IEP)*EP(IEQ)*SQRT(EP(IEP)+EP(IEQ)))
          INDP=0
          DO MP=0,NP      
            DO IP=0,MP
              DO JP=0,MP-IP
                KP=MP-IP-JP
                INDP=INDP+1
!         
                INDQ=0
                SGN=-1.D0
                DO MQ=0,NP      
                  SGN=-SGN     !SGN=(-1)^(IQ+JQ+KQ)=(-1)^MQ
                  DO IQ=0,MQ
                    DO JQ=0,MQ-IQ
                      KQ=MQ-IQ-JQ
                      INDQ=INDQ+1
                      UPQ(INDP,INDQ,IEP,IEQ)=FACTOR*SGN*HI(IP+IQ,JP+JQ,KP+KQ) 
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(HI)
!
!     ==========================================================================
!     ==  FOURCENTER INTEGRALS IN CARTESIAN GAUSSIANS                         ==
!     ==========================================================================
PRINT*,'IN GAUSSIAN$ZDIRECTION_FOURCENTER: MAPPING LOOP 1',NIJKP,NEP,NFA,NFB
      ALLOCATE(CPAB2(NIJKP,NEP,NFA,NFB))
      CPAB2=0.D0
      DO IFC=1,NFA
        DO IFD=1,NFB
          DO IEQ=1,NEP
            DO INDQ=1,NIJKP
              CPAB2(:,:,IFC,IFD)=CPAB2(:,:,IFC,IFD) &
                                +UPQ(:,INDQ,:,IEQ)*CPAB1(INDQ,IEQ,IFC,IFD)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
PRINT*,'IN GAUSSIAN$ZDIRECTION_FOURCENTER: MAPPING LOOP 2'
      UABCD(:,:,:,:)=0.D0
      DO IFC=1,NFA
        DO IFD=1,NFB
          DO IFA=1,NFA
            DO IFB=1,NFB
              IF(IFA+NFA*(IFB-1).GT.IFC+NFA*(IFD-1)) CYCLE
              SVAR=0.D0
              DO IEP=1,NEP
                DO INDP=1,NIJKP
                  SVAR=SVAR+CPAB1(INDP,IEP,IFA,IFB)*CPAB2(INDP,IEP,IFC,IFD)
                ENDDO
              ENDDO
              UABCD(IFA,IFB,IFC,IFD)=SVAR
              UABCD(IFC,IFD,IFA,IFB)=SVAR
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(CPAB1)
      DEALLOCATE(CPAB2)
PRINT*,'IN GAUSSIAN$ZDIRECTION_FOURCENTER: DONE'
!
                     CALL TRACE$POP()
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN$PLOTGAUSSEXPANSION(NIJK,E,C,RAD,N,F)
!     **************************************************************************
!     ** MAP AN FUNCTION EXPRESSED AS GAUSSIAN EXPANSION ONTO A GRID          **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NIJK     ! #(GAUSS COEFFICIENTS)
      REAL(8)   ,INTENT(IN) :: E        ! EXPONENT OF GAUSS EXPANSION
      REAL(8)   ,INTENT(IN) :: C(NIJK)  ! GAUSS COEFFICIENTS
      REAL(8)   ,INTENT(IN) :: RAD      ! RADIUS OF SPHERE ENCLOSED BY GRID
      INTEGER(4),INTENT(IN) :: N        ! DIMENSION OF GRID IN 3 DIRECTIONS
      REAL(8)   ,INTENT(OUT):: F(N,N,N) ! FUNCTION ON THE GRID
      REAL(8)               :: DELTA    ! STEP WIDTH
      REAL(8)               :: R(3)     ! ACTUAL POSITION ON THE GRID
      INTEGER(4)            :: I,J,K
!     **************************************************************************
      DELTA=2*RAD/REAL(N+1)
      R(3)=-RAD-DELTA
      DO K=1,N
        R(3)=R(3)+DELTA
        R(2)=-RAD-DELTA
        DO J=1,N
          R(2)=R(2)+DELTA
          R(1)=-RAD-DELTA
          DO I=1,N
            R(1)=R(1)+DELTA
            CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,E,C,R,F(I,J,K))
          ENDDO
        ENDDO
       ENDDO
      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN$PLOTRADIAL(FILE,RAD,NIJK,NE,E,C)
!     **************************************************************************
!     ** WRITE A GAUSSIAN EXPANSION ONTO A RADIAL GRID                        **
!     ** USE STARBURST MESH AND WRITE MEAN, MAX, AND MIN                      **
!     **************************************************************************
      IMPLICIT NONE
      CHARACTER(*),INTENT(IN) :: FILE   ! FILENAME OF OUTPUT
      REAL(8)   ,INTENT(IN) :: RAD      ! RADIUS OF SPHERE ENCLOSED BY GRID
      INTEGER(4),INTENT(IN) :: NIJK     ! #(GAUSS COEFFICIENTS)
      INTEGER(4),INTENT(IN) :: NE       ! #(GAUSS COEFFICIENTS)
      REAL(8)   ,INTENT(IN) :: E(NE)    ! EXPONENT OF GAUSS EXPANSION
      REAL(8)   ,INTENT(IN) :: C(NIJK,NE) ! GAUSS COEFFICIENTS
      INTEGER(4),PARAMETER  :: NR=200
      INTEGER(4),PARAMETER  :: NDIR=26
      REAL(8)               :: EXARR(NE)
      REAL(8)               :: VEC(3,NDIR)
      REAL(8)               :: F(NDIR)
      REAL(8)               :: SVAR,RI
      INTEGER(4)            :: IND,I,J,K
      INTEGER(4)            :: IR,IDIR,IE
      INTEGER(4)            :: NFIL
!     **************************************************************************
!
!     ==========================================================================
!     == FCC/SIC/BCC-TYPE STARBURST MESH                                      ==
!     ==========================================================================
!     == FCC DIRECTIONS ========================================================
      SVAR=1.D0/SQRT(2.D0)
      VEC(:, 1)=(/+1.D0,+0.D0,+1.D0/)*SVAR
      VEC(:, 2)=(/-1.D0,+0.D0,+1.D0/)*SVAR
      VEC(:, 3)=(/+0.D0,+1.D0,+1.D0/)*SVAR
      VEC(:, 4)=(/+0.D0,-1.D0,+1.D0/)*SVAR
      VEC(:, 5)=(/+1.D0,+1.D0,+0.D0/)*SVAR
      VEC(:, 6)=(/+1.D0,-1.D0,+0.D0/)*SVAR
      VEC(:, 7)=(/-1.D0,+1.D0,+0.D0/)*SVAR
      VEC(:, 8)=(/-1.D0,-1.D0,+0.D0/)*SVAR
      VEC(:, 9)=(/+1.D0,+0.D0,-1.D0/)*SVAR
      VEC(:,10)=(/-1.D0,+0.D0,-1.D0/)*SVAR
      VEC(:,11)=(/+0.D0,+1.D0,-1.D0/)*SVAR
      VEC(:,12)=(/+0.D0,-1.D0,-1.D0/)*SVAR
!     == SIC DIRECTIONS ========================================================
      VEC(:,13)=(/+1.D0,+0.D0,+0.D0/)
      VEC(:,14)=(/-1.D0,+0.D0,+0.D0/)
      VEC(:,15)=(/+0.D0,+1.D0,+0.D0/)
      VEC(:,16)=(/+0.D0,-1.D0,+0.D0/)
      VEC(:,17)=(/+0.D0,+0.D0,+1.D0/)
      VEC(:,18)=(/+0.D0,+0.D0,-1.D0/)
!     == BCC DIRECTIONS ========================================================
      SVAR=1.D0/SQRT(3.D0)
      VEC(:,19)=(/+1.D0,+1.D0,+1.D0/)*SVAR
      VEC(:,20)=(/-1.D0,+1.D0,+1.D0/)*SVAR
      VEC(:,21)=(/+1.D0,-1.D0,+1.D0/)*SVAR
      VEC(:,22)=(/-1.D0,-1.D0,+1.D0/)*SVAR
      VEC(:,23)=(/+1.D0,+1.D0,-1.D0/)*SVAR
      VEC(:,24)=(/-1.D0,+1.D0,-1.D0/)*SVAR
      VEC(:,25)=(/+1.D0,-1.D0,-1.D0/)*SVAR
      VEC(:,26)=(/-1.D0,-1.D0,-1.D0/)*SVAR
!
      CALL FILEHANDLER$SETFILE('HOOK',.FALSE.,FILE)
      CALL FILEHANDLER$UNIT('HOOK',NFIL)
      DO IR=1,NR
        RI=RAD*REAL(IR-1)/REAL(NR-1)
        EXARR(:)=EXP(-E(:)*RI**2)
        F(:)=0.D0
        DO IND=1,NIJK
          CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I,J,K)
          DO IDIR=1,NDIR
            DO IE=1,NE
              F(IDIR)=F(IDIR)+VEC(1,IDIR)**I*VEC(2,IDIR)**J*VEC(3,IDIR)**K &
     &                       *RI**(I+J+K)*EXARR(IE)*C(IND,IE)
            ENDDO
          ENDDO
        ENDDO
        WRITE(NFIL,*)RI,SUM(F)/REAL(NDIR),MAXVAL(F),MINVAL(F)
      ENDDO
      CALL FILEHANDLER$CLOSE('HOOK')
      CALL FILEHANDLER$SETFILE('HOOK',.TRUE.,'.FORGOTTOASSIGNFILETOHOOK')
      RETURN
      END      
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_HERMITEINTEGRAL(N,E,R,H0)
!     **************************************************************************
!     ** EVALUATES A HERMITE COULOMB INTEGRALS                                **
!     ** ONLY INTEGRALS WITH T+U+V.LE.N ARE CALCULATED. OTHERS REMAIN EMPTY.  **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: E
      REAL(8)   ,INTENT(IN) :: R(3)
      REAL(8)   ,INTENT(OUT):: H0(0:N,0:N,0:N)
      REAL(8)               :: HP(0:N,0:N,0:N)
      REAL(8)               :: F(0:3*N)
      INTEGER(4)            :: T,U,V,I
      INTEGER(4)            :: M
      REAL(8)               :: X
      REAL(8)               :: FAC,SVAR
!     **************************************************************************
!
!     ==========================================================================
!     ==  DETERMINE BOYS FUNCTIONS                                            ==
!     ==========================================================================
      X=E*(R(1)**2+R(2)**2+R(3)**2)
      CALL GAUSSIAN_BOYS(3*N,X,F)
      FAC=-2.D0*E
      SVAR=1.D0
      DO M=0,3*N
        F(M)=SVAR*F(M)    !(-2*E)^M * FM
        SVAR=SVAR*FAC     !(-2*E)^M
      ENDDO
!
!     ==========================================================================
!     ==  USE RECURSION TO CONSTRUCT HERMITE COULOMB INTEGRALS                ==
!     ==========================================================================
      H0(:,:,:)=0.D0
      DO I=0,N  ! ONLY DETERMINE ELEMENTS WITH T+U+V < N+1
        M=N-I
!       == ASSUME THAT ALL ELEMENTS WITH UPPER INDEX M+1 ARE THERE
!       == NON-ZERO ELEMENTS OF H0 HAVE T+U+V<I ================================
!       == ALL INDICES T,U,V MUST LIE IN [0,N]
!
!       == RAISE THIRD INDEX V =================================================
        IF(I.GE.1) THEN  ! AVOID ACCESSING ELEMENTS WITH INDICES <0
          V=1
          DO U=0,MIN(I-V,N)
            DO T=0,MIN(I-V-U,N)
              HP(T,U,V)=R(3)*H0(T,U,V-1)
            ENDDO
          ENDDO
        END IF
        DO V=2,MIN(I,N)  
          DO U=0,MIN(I-V,N)
            DO T=0,MIN(I-V-U,N)
              HP(T,U,V)=REAL(V-1)*H0(T,U,V-2)+R(3)*H0(T,U,V-1)
            ENDDO
          ENDDO
        ENDDO
!
!       == RAISE SECOND INDEX ==================================================
        IF(I.GE.1) THEN  ! AVOID ACCESSING ELEMENTS WITH INDICES <0
          U=1
          DO T=0,MIN(I-U,N)
            HP(T,U,0)=R(2)*H0(T,U-1,0)
          ENDDO
        END IF
        DO U=2,MIN(I,N)
          DO T=0,MIN(I-U,N)
            HP(T,U,0)=REAL(U-1)*H0(T,U-2,0)+R(2)*H0(T,U-1,0)
          ENDDO
        ENDDO
!
!       == RAISE FIRST INDEX ===================================================
        IF(I.GE.1) THEN  ! AVOID ACCESSING ELEMENTS WITH INDICES <0
          T=1
          HP(T,0,0)=R(1)*H0(T-1,0,0)
        END IF
        DO T=2,MIN(I,N)
          HP(T,0,0)=REAL(T-1)*H0(T-2,0,0)+R(1)*H0(T-1,0,0)
        ENDDO
!
!       == FILL IN (000) ELEMENT ===============================================
        HP(0,0,0)=F(M)
!
!       == SWAP ARRAYS =========================================================
        H0=HP
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_HERMITEINTEGRALOLD(N,E,R,HI)
!     **************************************************************************
!     ** EVALUATES A HERMITE COULOMB INTEGRALS                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: E
      REAL(8)   ,INTENT(IN) :: R(3)
      REAL(8)   ,INTENT(OUT):: HI(0:N,0:N,0:N)
      REAL(8)               :: F(0:3*N)
      INTEGER(4)            :: T,U,V
      INTEGER(4)            :: M
      REAL(8)               :: X
      REAL(8)               :: FAC,SVAR
!     **************************************************************************
!
!     ==========================================================================
!     ==  DETERMINE BOYS FUNCTIONS                                            ==
!     ==========================================================================
      X=E*(R(1)**2+R(2)**2+R(3)**2)
      CALL GAUSSIAN_BOYS(3*N,X,F)
      FAC=-2.D0*E
      SVAR=1.D0
      DO M=1,3*N
        F(M)=SVAR*F(M)    !(-2*E)^M * FM
        SVAR=SVAR*FAC     !(-2*E)^M
      ENDDO
!
!     ==========================================================================
!     ==  USE RECURSION TO CONSTRUCT HERMITE COULOMB INTEGRALS                ==
!     ==========================================================================
      DO M=1,3*N
        HI(0,0,0)=F(M)       ! LEVEL M PRODUCES ALL INTEGRALS WITH T+U+V=M
        DO T=MIN(M,N),2,-1   ! 2 .LE. T .LE. MIN(M,N)
          DO U=0,M-T         ! 0 .LE. U .LE. M-T 
            V=M-T-U          ! T+U+V=M
            HI(T,U,V)=REAL(T-1)*HI(T-2,U,V)+R(1)*HI(T-1,U,V)
          ENDDO
        ENDDO
        T=1
        DO U=0,M-T
          V=M-T-U
          HI(T,U,V)=R(1)*HI(T-1,U,V)
        ENDDO
        T=0
        DO U=0,M-T
          V=M-T-U
          HI(T,U,V)=0.D0
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_BOYS(N,X,F)
!     **************************************************************************
!     ** EVALUATES A SET OF BOYS FUNCTIONS FOR M=0,...,N                      **
!     **          F_M(X)=INT_0^1 DT: T^(2M)*EXP(-X*T^2)                       **
!     ** RESULT IS ACCURATE ONLY FOR ABOUT X<D (EMPIRICAL OBSERVATION)        **
!     **                                                                      **
!     ** METHOD USED IS A VARIATION OF:                                       **
!     **   B.A. MAMEDOV, ON THE EVALUATION OF BOYS FUNCTIONS USING DOWNWARD   **
!     **   RECURSION RELATION, J. MATH. CHEM. 36, P301                        **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: N
      REAL(8)   ,INTENT(IN) :: X
      REAL(8)   ,INTENT(OUT):: F(0:N)
      REAL(8)   ,PARAMETER  :: D=15.D0   ! #(ACCURATE DIGITS) 
      REAL(8)   ,PARAMETER  :: PI=4.D0*ATAN(1.D0)
      REAL(8)               :: Y,EXPMX,SVAR,F0
      INTEGER(4)            :: MT,M,MBREAK
!     **************************************************************************
!
!     ==========================================================================
!     == CATCH THE CASE WHEN X=0                                              ==
!     ==========================================================================
      IF(ABS(X).LT.1.D-5) THEN
        DO M=0,N
          F(M)=1.D0/REAL(2*M+1,KIND=8) &
     &        -X/REAL(2*M+2,KIND=8) &
     &        -X**2/REAL(2*(2*M+3),KIND=8) &
     &        -X**3/REAL(6*(2*M+3),KIND=4) 
        ENDDO
        RETURN
      END IF
!
!     ==========================================================================
!     == SET UP DATA FOR THE MAIN PART                                        ==
!     ==========================================================================
      IF(X.LT.35.D0) THEN
        CALL LIB$ERFR8(SQRT(X),SVAR)
        F0=SQRT(PI/(4.D0*X))*SVAR  ! BOYS FUNCTION F0
      ELSE
        F0=SQRT(PI/(4.D0*X))              ! LARGE ARGUMENT LIMIT OF F0
      END IF
      MBREAK=NINT(X)
      EXPMX=EXP(-X)
!
!     ==========================================================================
!     == LARGE DISTANCE LIMIT                                                 ==
!     ==========================================================================
      SVAR=1.D0/(2.D0*X)
      Y=F0
      DO M=0,MIN(MBREAK,N)
        F(M)=Y
        Y=(Y*REAL(2*M+1,KIND=8)-EXPMX)*SVAR
      ENDDO
!
!     ==========================================================================
!     == NOW THE NORMAL RECURSION                                             ==
!     ==========================================================================
      IF(ABS(X-REAL(N)).GT.1.D-5) THEN   !EQ.10 OF MAMEDOV04_JCX36_301
        MT=N+INT(D/ABS(LOG10(REAL(N)/X)))
      ELSE
        MT=N+INT(D/ABS(LOG10(REAL(N))))
      END IF   
      MT=2*INT(0.5D0*REAL(MT+2))  !MT MUST BE THE NEXT HIGHER EVEN NUMBER
      Y=EXPMX/REAL(2*M+1,KIND=8)
      DO M=MT,MBREAK+1,-1
        Y=(2.D0*X*Y+EXPMX)/REAL(2*M+1,KIND=8) !=F(M-1)
        IF(M.LE.N)F(M)=Y
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN$YLMTIMESRN(J,LM,NIJK,C)
!     **************************************************************************
!     **  PRODUCE COEFFICIENT ARRAY FOR A POLYNOMIAL R^(L+2*J)*Y_LM           **
!     **                                                                      **
!     ******************************PETER BLOECHL, GOSLAR 2011******************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: J
      INTEGER(4),INTENT(IN)   :: LM
      INTEGER(4),INTENT(IN)   :: NIJK
      REAL(8)   ,INTENT(OUT)  :: C(NIJK)
      LOGICAL(4),SAVE         :: TINI=.FALSE.
      INTEGER(4),SAVE         :: NX=-1
      INTEGER(4),SAVE         :: LMX=-1
      INTEGER(4),SAVE         :: NIJKX=-1
      INTEGER(4)              :: NF
      REAL(8)   ,POINTER,SAVE :: CMAT(:,:)
      INTEGER(4),POINTER,SAVE :: MAP(:,:)
!     **************************************************************************
!     ==========================================================================
!     == CHECK SIZE OF TABLE
!     ==========================================================================
      IF(2*J.GT.NX.OR.LM.GT.LMX.OR.NIJK.GT.NIJKX) THEN
        IF(TINI) THEN
          DEALLOCATE(CMAT)
          DEALLOCATE(MAP)
          TINI=.FALSE.
        END IF
      END IF
!
!     ==========================================================================
!     == REBUILD TABLE IF NEEDED                                              ==
!     ==========================================================================
      IF(.NOT.TINI) THEN
        NX=MAX(4,2*J) ! THE 4 IS TO AVOID BUILDING UNLIKELY SMALL TABLE SIZES
        DO 
          NIJKX=(NX+1)*(NX+2)*(NX+3)/6
          LMX=(NX+1)**2
          IF(NIJKX.GE.NIJK.AND.LMX.GE.LM) EXIT
          NX=NX+1
        ENDDO
        NF=(NX+2)*((NX+2)**2-1)/6
        ALLOCATE(CMAT(NIJKX,NF))
        ALLOCATE(MAP(NX+1,LMX))
        CALL GAUSSIAN_EXPANDMAT(NX,NIJKX,NF,CMAT,MAP)
        TINI=.TRUE.
      END IF
!
!     ==========================================================================
!     == COLLECT RESULT                                                       ==
!     ==========================================================================
      C(:)=CMAT(:NIJK,MAP(J+1,LM))
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_EXPANDMAT(NX,NIJK,NF,C,MAP)
!     **************************************************************************
!     **  EXPRESS FUNCTIONS R^(L+2*J)*Y_L(R) IN A GAUSS REPRESENTATION        **
!     **  (L,M,N)-> C(:NIJK,MAP(J+1,LM))                                      **
!     **                                                                      **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NX                  ! HIGHEST POWER
      INTEGER(4),INTENT(IN) :: NIJK                ! #(TERMS IN GAUSSIAN EXPANS)
      INTEGER(4),INTENT(IN) :: NF                  ! #(ENTRIES IN TRANSF.)
      REAL(8)   ,INTENT(OUT):: C(NIJK,NF)          ! TRANSFORMATION MATRIX
      INTEGER(4),INTENT(OUT):: MAP(NX+1,(NX+1)**2) ! LOOKUPTABLE
      INTEGER(4)            :: LX
      INTEGER(4)            :: INDX
      INTEGER(4)            :: IND1X
      REAL(8)               :: R2N(NIJK)
      REAL(8)               :: YLMPOL(NIJK,(NX+1)**2)
      INTEGER(4)            :: N,I,J,K,I1,J1,K1,IORB,L,IM,LM,IND,IND1,IND2
      INTEGER(4)            :: TWOI,TWOJ,TWOK
      REAL(8)               :: SVAR1,SVAR2
!     ******************************PETER BLOECHL, GOSLAR 2011******************
!                            CALL TRACE$PUSH('GAUSSIAN_EXPANDMAT')
      IF(NIJK.NE.(NX+1)*(NX+2)*(NX+3)/6) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$STOP('GAUSSIAN_EXPANDMAT')
      END IF
      IF(NF.NE.(NX+2)*((NX+2)**2-1)/6) THEN
        CALL ERROR$MSG('INCONSISTENT ARRAY SIZE')
        CALL ERROR$STOP('GAUSSIAN_EXPANDMAT')
      END IF
      LX=NX
      MAP(:,:)=0
      C(:,:)=0.D0
!
!     ==========================================================================
!     == POLYNOMIAL COEFFICIENTS OF SPHERICAL HARMONICS TIMES R**L            ==
!     ==========================================================================
      CALL GAUSSIAN_YLMPOL(LX,YLMPOL)
!
!     ==========================================================================
!     ==  DETERMINE EXPANSION COEFFICIENTS                                    ==
!     ==========================================================================
      IORB=0
      DO N=0,NX/2
!
!       ========================================================================
!       ==  MULTIPLY WITH (X^2+Y^2+Z^2)^N                                     ==
!       ========================================================================
        INDX=0
        R2N(:)=0.D0
        DO I=0,N
          TWOI=2*I  ! AVOID INTENT(IN/OUT) CONFUSION
          CALL BINOMIALCOEFFICIENT(N,I,SVAR1)
          DO J=0,N-I
            TWOJ=2*J  ! AVOID INTENT(IN/OUT) CONFUSION
            CALL BINOMIALCOEFFICIENT(N-I,J,SVAR2)
            K=N-I-J
            TWOK=2*K  ! AVOID INTENT(IN/OUT) CONFUSION
            CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND1,TWOI,TWOJ,TWOK)
            IF(IND1.GT.NIJK) CYCLE
            R2N(IND1)=R2N(IND1)+SVAR1*SVAR2
            INDX=MAX(INDX,IND1)
          ENDDO
        ENDDO
!
!       ========================================================================
!       == MULTIPLY WITH SPHERICAL HARMONICS                                  ==
!       ========================================================================
        DO L=0,NX-2*N
          DO IM=1,2*L+1
            LM=L**2+IM
            IORB=IORB+1
            MAP(N+1,LM)=IORB   ! LOOKUPTABLE
            IND1X=SIZE(YLMPOL(:,LM))
            DO IND1=1,IND1X
              IF(YLMPOL(IND1,LM).EQ.0.D0) CYCLE
              CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND1,I,J,K)
              DO IND=1,INDX
                CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I1,J1,K1)
                I1=I1+I
                J1=J1+J
                K1=K1+K
                CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND2,I1,J1,K1)
                IF(IND2.GT.NIJK) CYCLE
                C(IND2,IORB)=C(IND2,IORB)+R2N(IND)*YLMPOL(IND1,LM)
              ENDDO
            ENDDO
          ENDDO ! END OF LOOP OVER MAGNETIC QUANTUM NUMBERS
        ENDDO  ! END OF LOOP OVER POWERS N
      ENDDO  ! END OF LOOP OVER MAIN ANGULAR MOMENTUM L
!                            CALL TRACE$POP()
      RETURN
      END
!
!*******************************************************************************
!*******************************************************************************
!** TEST ROUTINES FOR THE GAUSSIAN OBJECT                                     **
!*******************************************************************************
!*******************************************************************************
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_TEST_ALL()
      IMPLICIT NONE
!     **************************************************************************
      CALL GAUSSIAN_TEST_FOURCENTER()
      CALL GAUSSIAN_TEST_CONTRACTGAUSS()
      CALL GAUSSIAN_TEST_INDEX()
      CALL GAUSSIAN_TEST_PLOTGAUSS()
      CALL GAUSSIAN_TEST_PRODUCTRULE()
      CALL GAUSSIAN_TEST_BOYS()
      STOP
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_TEST_HERMITEINTEGRAL()
!     **************************************************************************
!     ** TESTROUTINE FOR GAUSSIAN_HERMITEINTEGRAL                             **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: N=2
      REAL(8)              :: E=1.3D0
      REAL(8)              :: R(3)
      REAL(8)              :: HI(0:N,0:N,0:N)
      REAL(8)              :: HIREF1(0:N,0:N,0:N)
      REAL(8)              :: HIREF2(0:N,0:N,0:N)
      REAL(8)              :: F0(0:3*N),FP(0:3*N),FM(0:3*N)
      REAL(8)              :: DF(0:3*N),D2F(0:3*N)
      REAL(8)              :: DF2(0:3*N),D2F2(0:3*N)
      REAL(8)  ,PARAMETER  :: DELTA=1.D-3
      REAL(8)              :: X
      INTEGER(4)           :: I,J,K
!     **************************************************************************
      CALL RANDOM_NUMBER(R(1))
      CALL RANDOM_NUMBER(R(2))
      CALL RANDOM_NUMBER(R(3))
      X=E*SUM(R(:)**2)
      CALL GAUSSIAN_BOYS(3*N,X,F0)
      CALL GAUSSIAN_BOYS(3*N,X+1.D0*DELTA,FP)
      CALL GAUSSIAN_BOYS(3*N,X-1.D0*DELTA,FM)
      DF =(FP-FM)/(2.D0*DELTA)
      D2F=(FP-2.D0*F0+FM)/DELTA**2
      CALL GAUSSIAN_BOYS(3*N,X+2.D0*DELTA,FP)
      CALL GAUSSIAN_BOYS(3*N,X-2.D0*DELTA,FM)
      DF2 =(FP-FM)/(4.D0*DELTA)
      D2F2=(FP-2.D0*F0+FM)/(2.D0*DELTA)**2
!
      HIREF1(0,0,0)=F0(0)
      HIREF1(1,0,0)=2.D0*E*R(1)*DF(0)      
      HIREF1(0,1,0)=2.D0*E*R(2)*DF(0)      
      HIREF1(0,0,1)=2.D0*E*R(3)*DF(0)      
      HIREF1(2,0,0)=2.D0*E*DF(0)+(2.D0*E*R(1))**2*D2F(0)      
      HIREF1(0,2,0)=2.D0*E*DF(0)+(2.D0*E*R(2))**2*D2F(0)      
      HIREF1(0,0,2)=2.D0*E*DF(0)+(2.D0*E*R(3))**2*D2F(0)      
      HIREF1(1,1,0)=(2.D0*E)**2*R(1)*R(2)*D2F(0)      
      HIREF1(1,0,1)=(2.D0*E)**2*R(1)*R(3)*D2F(0)      
      HIREF1(0,1,1)=(2.D0*E)**2*R(2)*R(3)*D2F(0)      
!
      HIREF2(0,0,0)=F0(0)
      HIREF2(1,0,0)=2.D0*E*R(1)*DF2(0)      
      HIREF2(0,1,0)=2.D0*E*R(2)*DF2(0)      
      HIREF2(0,0,1)=2.D0*E*R(3)*DF2(0)      
      HIREF2(2,0,0)=2.D0*E*DF2(0)+(2.D0*E*R(1))**2*D2F2(0)      
      HIREF2(0,2,0)=2.D0*E*DF2(0)+(2.D0*E*R(2))**2*D2F2(0)      
      HIREF2(0,0,2)=2.D0*E*DF2(0)+(2.D0*E*R(3))**2*D2F2(0)      
      HIREF2(1,1,0)=(2.D0*E)**2*R(1)*R(2)*D2F2(0)      
      HIREF2(1,0,1)=(2.D0*E)**2*R(1)*R(3)*D2F2(0)      
      HIREF2(0,1,1)=(2.D0*E)**2*R(2)*R(3)*D2F2(0)      

      CALL GAUSSIAN_HERMITEINTEGRAL(N,E,R,HI)

      DO I=0,N
        DO J=0,N
          DO K=0,N
            IF(I+J+K.GT.2) CYCLE
            WRITE(*,FMT='(3I3,F10.5,3E10.3)')I,J,K,HI(I,J,K) &
    &                    ,HIREF1(I,J,K)-HI(I,J,K),HIREF2(I,J,K)-HI(I,J,K)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_TEST_FOURCENTER()
!     **************************************************************************
!     ** TESTROUTINE FOR GAUSSIAN_FOURCENTER                                  **
!     ** CURRENTLY IT JUST CALCULATES A MATRIX ELEMENT BETWEEN FUNCTIONS      **
!     ** EXPRESSED AS EXPANSIONS IN CARTESIAN GAUSSIANS                       **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) :: NIJKA,NIJKB,NIJKC,NIJKD,NIJKP,NIJKQ
      REAL(8)    :: EA(1),EB(1),EC(1),ED(1)
      REAL(8)    :: RA(3),RB(3),RC(3),RD(3)
      REAL(8),ALLOCATABLE:: CA(:,:,:),CB(:,:,:),CC(:,:,:),CD(:,:,:)
      REAL(8)    :: UABCD(1,1,1,1)
      INTEGER(4) :: NA,NB,NC,ND,NP,NQ
      INTEGER(4) :: I,J,K
!     **************************************************************************
      NA=2
      NB=2
      NC=2
      ND=2
      NP=9  ! RESULTS BECOME VERY POOR IF NP<NA+NB!!
      NQ=9  ! RESULTS BECOME VERY POOR IF NQ<NC+ND!!
      CALL RANDOM_NUMBER(EA(1))
      CALL RANDOM_NUMBER(EB(1))
      CALL RANDOM_NUMBER(EC(1))
      CALL RANDOM_NUMBER(ED(1))
      DO I=1,3
        CALL RANDOM_NUMBER(RA(I))
        CALL RANDOM_NUMBER(RB(I))
        CALL RANDOM_NUMBER(RC(I))
        CALL RANDOM_NUMBER(RD(I))
      ENDDO
      J=0  ! AVOID INTENT(IN/OUT) CONFUSION
      K=0  ! AVOID INTENT(IN/OUT) CONFUSION
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKA,NA,J,K)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKB,NB,J,K)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKC,NC,J,K)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKD,ND,J,K)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKP,NP,J,K)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKQ,NQ,J,K)
      ALLOCATE(CA(NIJKA,1,1))
      ALLOCATE(CB(NIJKB,1,1))
      ALLOCATE(CC(NIJKC,1,1))
      ALLOCATE(CD(NIJKD,1,1))
      DO I=1,NIJKA
        CALL RANDOM_NUMBER(CA(I,1,1))
      ENDDO
      DO I=1,NIJKB
        CALL RANDOM_NUMBER(CB(I,1,1))
      ENDDO
      DO I=1,NIJKC
        CALL RANDOM_NUMBER(CC(I,1,1))
      ENDDO
      DO I=1,NIJKD
        CALL RANDOM_NUMBER(CD(I,1,1))
      ENDDO
 
      CALL GAUSSIAN$FOURCENTER(NIJKA,1,1,EA,RA,CA,NIJKB,1,1,EB,RB,CB &
     &                        ,NIJKC,1,1,EC,RC,CC,NIJKD,1,1,ED,RD,CD,UABCD)
PRINT*,'UABCD',UABCD
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_TEST_CONTRACTGAUSS()
!     **************************************************************************
!     ** TESTROUTINE
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) :: NA
      INTEGER(4) :: NB
      INTEGER(4) :: NP
      INTEGER(4)           :: I,J,K,IPROBE
      REAL(8)              :: EA,EB,EP
      REAL(8)              :: RA(3),RB(3),RP(3)
      INTEGER(4)           :: NIJKA,NIJKB,NIJKP
      INTEGER(4)           :: NPROBE
      REAL(8)              :: R(3)
      REAL(8),ALLOCATABLE  :: CA(:),CB(:),CP(:)
      REAL(8)              :: GA,GB,GP
!     **************************************************************************
      NA=6
      NB=2
      NP=9  ! RESULTS BECOME VERY POOR IF NP<NA+NB!!
      CALL RANDOM_NUMBER(EA)
      CALL RANDOM_NUMBER(EB)
      DO I=1,3
        CALL RANDOM_NUMBER(RA(I))
        CALL RANDOM_NUMBER(RB(I))
      ENDDO
      J=0  !AVOID INTENT(IN/OUT) CONFUSION
      K=0  !AVOID INTENT(IN/OUT) CONFUSION
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKA,NA,J,K)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKB,NB,J,K)
      CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',NIJKP,NP,J,K)
      ALLOCATE(CA(NIJKA))
      ALLOCATE(CB(NIJKB))
      ALLOCATE(CP(NIJKP))
      DO I=1,NIJKA
        CALL RANDOM_NUMBER(CA(I))
      ENDDO
      DO I=1,NIJKB
        CALL RANDOM_NUMBER(CB(I))
      ENDDO

      CALL GAUSSIAN_CONTRACT2GAUSS(NIJKA,EA,RA,CA,NIJKB,EB,RB,CB &
     &                           ,NIJKP,EP,RP,CP)
      NPROBE=10
      DO IPROBE=1,NPROBE
        DO I=1,3
          CALL RANDOM_NUMBER(R(I))
        ENDDO
        CALL GAUSSIAN_3DORB('CARTESIAN',NIJKA,EA,CA,R-RA,GA)
        CALL GAUSSIAN_3DORB('CARTESIAN',NIJKB,EB,CB,R-RB,GB)
        CALL GAUSSIAN_3DORB('HERMITE',NIJKP,EP,CP,R-RP,GP)
        WRITE(*,FMT='("R=",3F10.5," A*B ",F20.10," P=",F20.10," DIFF ",F20.10)')R,GA*GB,GP,GP-GA*GB
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_TEST_INDEX()
!     **************************************************************************
!     ** TESTROUTINE
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4) :: IND,I,J,K,IND1,I1,J1,K1,N,M
!     **************************************************************************
      N=5
      IND=0
      DO M=0,N
        DO I=0,M
          DO J=0,M-I
            K=M-I-J
            IND=IND+1
            CALL GAUSSIAN_GAUSSINDEX('IJKFROMIND',IND,I1,J1,K1)
            CALL GAUSSIAN_GAUSSINDEX('INDFROMIJK',IND1,I1,J1,K1)
            WRITE(*,FMT='("IND=",I5," IND2=",I5," N=",I3," N=",I3," (I,J,K)=(",3I3,"   (I,J,K)=(",3I3,")")') &
     &                 IND,IND1,M,I1+J1+K1,I,J,K,I1,J1,K1
            IF(IND1.NE.IND.OR.M.NE.I1+J1+K1 &
     &                    .OR.I1.NE.I.OR.J1.NE.J.OR.K1.NE.K) THEN
              STOP 'TEST FAILED!!!'
            END IF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_TEST_PRODUCTRULE()
!     **************************************************************************
!     ** TESTROUTINE
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NX=4
      INTEGER(4),PARAMETER :: NFIL=100
      REAL(8)              :: H(0:NX,0:NX,0:2*NX)
      INTEGER(4)           :: IA,IB
      REAL(8)              :: XA,EA
      REAL(8)              :: XB,EB
      REAL(8)              :: XP,EP
      INTEGER(4)           :: NP
      REAL(8)              :: X,Y1,Y2,SVAR
      INTEGER(4)           :: I,M
!     **************************************************************************
      IA=4
      IB=3
      XA=-0.5D0
      EA=1.D0
      XB=0.5D0
      EB=2.D0
!
      EP=EA+EB
      XP=(EA*XA+EB*XB)/EP
      CALL GAUSSIAN_HERMITEC(NX,XA,XB,EA,EB,H)
      OPEN(UNIT=NFIL,FILE='DAT')
      NP=100
      DO I=1,NP
        X=-3.D0+6.D0*REAL(I-1)/REAL(NP-1)
        CALL GAUSSIAN_CARTESIANGAUSSIAN1D(IA,EA,X-XA,Y1)
        CALL GAUSSIAN_CARTESIANGAUSSIAN1D(IB,EB,X-XB,SVAR)
        Y1=Y1*SVAR
        Y2=0.D0
        DO M=0,IA+IB
          CALL GAUSSIAN_HERMITEGAUSSIAN1D(M,EP,X-XP,SVAR)
          Y2=Y2+H(IA,IB,M)*SVAR
        ENDDO
        WRITE(NFIL,*)X,Y1,Y2 !,Y2-Y1
      ENDDO
      CLOSE(NFIL)
PRINT*,'H',H(IA,IB,0:IA+IB)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_TEST_PLOTGAUSS()
!     **************************************************************************
!     ** TESTROUTINE
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NFIL=10
      INTEGER(4),PARAMETER :: NX=100
      INTEGER(4),PARAMETER :: NF=6
      REAL(8)              :: E=1.D0
      REAL(8)              :: X,Y(NF)
      INTEGER(4)           :: I,J
!     **************************************************************************
!
!     ==========================================================================
!     == PLOT GAUSSIANS                                                       ==
!     ==========================================================================
      OPEN(UNIT=NFIL,FILE='DAT')
      DO I=1,NX
        X=-3+6.D0*REAL(I-1)/REAL(NX-1)
        DO J=1,3
          CALL GAUSSIAN_CARTESIANGAUSSIAN1D(J-1,E,X,Y(J))
        ENDDO
        DO J=1,3
          CALL GAUSSIAN_HERMITEGAUSSIAN1D(J-1,E,X,Y(3+J))
        ENDDO
        WRITE(NFIL,*)X,Y(4:NF)
      ENDDO
      CLOSE(NFIL)
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_TEST_BOYS()
!     **************************************************************************
!     ** TESTROUTINE FOR GAUSSIAN_BOYS                                        **
!     ** COMPARISON WITH DATA FROM B.A. MAMEDOV, J. MATH CHEM 36, P301 (2004) **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: N1=12,N2=6,N3=6
      REAL(8)              :: F1(N1),F2(N2),F3(N3)
      REAL(8)              :: X1(N1),X2(N2),X3(N3)
      INTEGER(4)           :: M1(N1),M2(N2),M3(N3)
      REAL(8)              :: F(0:100)
      REAL(8)              :: DEV
      REAL(8)   ,PARAMETER :: TOL=1.D-8
      INTEGER(4)           :: I
!     **************************************************************************
!     == TABLE 1 ===============================================================
      M1=(/8,15,20,25,31,11,42,75,100,20,45,100/)
      X1=(/16.D0,27.D0,30.D0,13.D0,34.D0,38.D0,32.D0,30.D0,33.D0 &
     &   ,1.4D-3,6.4D-5,2.6D-7/)
      F1=(/4.02308592502660D-07,1.08359515555596D-11,1.37585444267909D-13 &
     &    ,8.45734447905704D-08,2.90561943091301D-16,4.04561442253925D-12 &
     &    ,5.02183610419086D-16,1.01429517438537D-15,3.42689684943483D-17 &
     &    ,2.43577075309547D-02,1.09883228385254D-02,4.97512309732144D-03/)
!
!     == TABLE 2 ===============================================================
      M2=(/8,16,21,12,15,18/)
      X2=(/42.D0,50.D0,56.D0,60.D0,53.D0,58.D0/)
      F2=(/1.11826597752251D-10,2.40509456111904D-16,1.43739730342730D-19 &
     &    ,4.05791663779760D-15,3.14434039868936D-16,1.78336953967902D-18/)
!
!     == TABLE 3 ===============================================================
      M3=(/8,14,20,33,36,100/)
      X3=(/63.D0,68.D0,73.D0,85.D0,100.D0,120.D0/)
      F3=(/3.56261924865627D-12,3.09783511327517D-17,1.71295886102040D-21 &
     &    ,1.74268831008018D-29,3.08919970425521D-33,4.97723065221079D-53/)
!
!     ==========================================================================
      DO I=1,N1
        CALL GAUSSIAN_BOYS(M1(I),X1(I),F)
        DEV=ABS((F(M1(I))-F1(I))/F1(I))
        WRITE(*,FMT='("N=",I4," X=",F10.3," F=",2E25.10," DEV= ",E10.1)') &
     &        M1(I),X1(I),F(M1(I)),F1(I),DEV
        IF(DEV.GT.TOL) THEN
PRINT*,'ERROR================'
!!$          CALL ERROR$MSG('TEST OF GAUSSIAN_BOYS FAILED')
!!$          CALL ERROR$STOP('GAUSSIAN_TEST_BOYS')
        END IF
      ENDDO
!
!     ==========================================================================
      DO I=1,N2
        CALL GAUSSIAN_BOYS(M2(I),X2(I),F)
        DEV=ABS((F(M2(I))-F2(I))/F2(I))
        WRITE(*,FMT='("N=",I4," X=",F10.3," F=",2E25.10," DEV= ",E10.1)') &
     &        M2(I),X2(I),F(M2(I)),F2(I),DEV
        IF(DEV.GT.TOL) THEN
PRINT*,'ERROR================'
!!$          CALL ERROR$MSG('TEST OF GAUSSIAN_BOYS FAILED')
!!$          CALL ERROR$STOP('GAUSSIAN_TEST_BOYS')
        END IF
      ENDDO
!
!     ==========================================================================
      DO I=1,N3
        CALL GAUSSIAN_BOYS(M3(I),X3(I),F)
        DEV=ABS((F(M3(I))-F3(I))/F3(I))
        WRITE(*,FMT='("N=",I4," X=",F10.3," F=",2E25.10," DEV= ",E10.1)') &
     &        M3(I),X3(I),F(M3(I)),F3(I),DEV
        IF(DEV.GT.TOL) THEN
PRINT*,'ERROR================'
!!$          CALL ERROR$MSG('TEST OF GAUSSIAN_BOYS FAILED')
!!$          CALL ERROR$STOP('GAUSSIAN_TEST_BOYS')
        END IF
      ENDDO
STOP 'FORCED'
      RETURN
      END
!
!     ...1.........2.........3.........4.........5.........6.........7.........8
      SUBROUTINE GAUSSIAN_TEST_SHIFTCENTER()
!     **************************************************************************
!     **  TESTROUTINE FOR GAUSSIAN_SHIFTCENTER                                **
!     **************************************************************************
      IMPLICIT NONE
      INTEGER(4),PARAMETER :: NX=10
      INTEGER(4),PARAMETER :: NIJK=(NX+1)*(NX+2)*(NX+3)/6
      INTEGER(4),PARAMETER :: NP=100
      REAL(8)              :: E
      REAL(8)              :: DR(3)
      REAL(8)              :: R(3)
      REAL(8) ,ALLOCATABLE :: T(:,:) !(NIJK,NIJK)
      REAL(8)              :: C0(NIJK)
      REAL(8)              :: C1(NIJK)
      REAL(8)              :: F1,F0
      INTEGER(4)           :: NFIL=6
      INTEGER(4)           :: I
!     **************************************************************************
      ALLOCATE(T(NIJK,NIJK))
      E=1.D0
      DR(:)=(/1.D0,0.D0,0.D0/)
      C0(:)=0.D0
      C0(4)=1.D0
!
      CALL GAUSSIAN_SHIFTCENTER(NX,NIJK,E,DR,T)
      C1=MATMUL(T,C0)
!
      OPEN(NFIL,FILE='DAT',FORM='FORMATTED')
      DO I=1,NP
        R(:)=0.D0
        R(1)=-3.D0+6.D0*REAL(I-1)/REAL(NP-1)
        CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,E,C0,R-DR,F0)
        CALL GAUSSIAN_3DORB('CARTESIAN',NIJK,E,C1,R,F1)
        WRITE(NFIL,*)R(1),F0,F1,F1-F0
      ENDDO
      CLOSE(NFIL)
      STOP
      RETURN
      END
